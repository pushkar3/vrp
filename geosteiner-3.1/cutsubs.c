/***********************************************************************

	File:	cutsubs.c
	Rev:	a-1
	Date:	01/16/2001

	Copyright (c) 1996, 2001 by David M. Warme

************************************************************************

	Processing of violated cutsets.

************************************************************************

	Modification Log:

	a-1:	01/16/2001
		: Split off from cutset.c.  Modified to take cut
		:  vertices rather than set of spanning edges.

************************************************************************/

#include "bb.h"
#include "constrnt.h"
#include "cutset.h"
#include "steiner.h"


/*
 * Global Routines
 */

struct constraint * add_cutset_to_list (bitmap_t *,
					struct constraint *,
					double *,
					bitmap_t *,
					bitmap_t *,
					struct cinfo *);


/*
 * Local Routines
 */

static bool		is_subset (bitmap_t *, bitmap_t *, int);

/*
 * This routine handles all of the details of adding a new cutset to the
 * given list of cutsets.  Any existing constraints that are looser
 * (supersets) are deleted.  The new cutset is ignored if it is looser
 * than any existing cutset.
 */

	struct constraint *
add_cutset_to_list (

bitmap_t *		cutset,		/* IN - new cutset to add */
struct constraint *	cutlist,	/* IN - list to add to */
double *		x,		/* IN - current LP solution */
bitmap_t *		vert_mask,	/* IN - set of valid terminals */
bitmap_t *		edge_mask,	/* IN - set of valid hyperedges */
struct cinfo *		cip		/* IN - compatibility info */
)
{
int			i;
int			j;
int			nedges;
int			nmasks;
int			num_in_cut;
int			count;
int *			vp1;
int *			vp2;
struct constraint *	p;
struct constraint **	hookp;
bitmap_t *		cut_edges;
double			z;

	nedges	= cip -> num_edges;
	nmasks	= cip -> num_edge_masks;

	cut_edges = NEWA (nmasks, bitmap_t);
	memset (cut_edges, 0, nmasks * sizeof (*cut_edges));

	count = 0;
	z = 0.0;
	for (i = 0; i < cip -> num_edges; i++) {
		if (NOT BITON (edge_mask, i)) continue;
		num_in_cut = 0;
		vp1 = cip -> edge [i];
		vp2 = cip -> edge [i + 1];
		while (vp1 < vp2) {
			j = *vp1++;
			if (BITON (cutset, j)) {
				++num_in_cut;
			}
		}
		if (num_in_cut <= 0) {
			/* this hyperedge resides entirely	*/
			/* outside of the cut...  doesn't span!	*/
			continue;
		}
		if (num_in_cut >= cip -> edge_size [i]) {
			/* this hyperedge resides entirely	*/
			/* within the cut...  doesn't span!	*/
			continue;
		}
		SETBIT (cut_edges, i);
		++count;
		z += x [i];
	}

	/* Check for an all-zero cutset.  These occasionally	*/
	/* happen because of numeric issues...			*/

	if (count <= 0) {
		/* Empty cutset!  OOOPS! */
#if 1
		tracef (" %% WARNING!  empty cutset!\n");
#endif
		free ((char *) cut_edges);
		return (cutlist);
	}

	if (z >= 1.0 - FUZZ) {
#if 1
		tracef (" %% WARNING!  bogus cutset!\n");
#endif
		free ((char *) cut_edges);
		return (cutlist);
	}

	/* If this new cutset is a superset of an existing one,	*/
	/* then there is nothing to add, and nothing to delete.	*/
	for (p = cutlist; p NE NULL; p = p -> next) {
		if (is_subset (p -> mask, cut_edges, nmasks)) {
			free (cut_edges);
			return (cutlist);
		}
	}

	/* Delete all current cutsets which have this new one	*/
	/* as a subset.						*/
	hookp = &cutlist;
	while ((p = *hookp) NE NULL) {
		if (p -> type NE CT_CUTSET) {
			hookp = &(p -> next);
		}
		else if (is_subset (cut_edges, p -> mask, nmasks)) {
			*hookp = p -> next;
			free ((char *) (p -> mask));
			free ((char *) p);
		}
		else {
			hookp = &(p -> next);
		}
	}

	p = NEW (struct constraint);
	p -> next	= NULL;
	p -> iteration	= 0;
	p -> type	= CT_CUTSET;
	p -> mask	= cut_edges;
	*hookp = p;

	return (cutlist);
}

/*
 * This routine returns TRUE if-and-only-if the first bit-mask is
 * a subset of the second.
 */

	static
	bool
is_subset (

bitmap_t *	bp1,		/* IN - first set. */
bitmap_t *	bp2,		/* IN - second set. */
int		nmasks		/* IN - number of masks in set. */
)
{
int		i;

	for (i = 0; i < nmasks; i++) {
		if ((*bp1 & *bp2) NE *bp1) return (FALSE);
		++bp1;
		++bp2;
	}

	return (TRUE);
}
