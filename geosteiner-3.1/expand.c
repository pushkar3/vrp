/***********************************************************************

	File:	expand.c
	Rev:	a-1
	Date:	01/16/2001

	Copyright (c) 1996, 2001 by David M. Warme

************************************************************************

	Expand logical constraints into physical constraints.

************************************************************************

	Modification Log:

	a-1:	01/16/2001	warme
		: Split off from constrnt.c.

************************************************************************/

#include "bb.h"
#include "constrnt.h"
#include "steiner.h"


/*
 * Global Routines
 */

struct rcoef *	expand_constraint (struct constraint *,
				   struct rcoef *,
				   bitmap_t *,
				   struct cinfo *);

/*
 * This routine expands a "logical" constraint into its coefficient/OP/RHS
 * form.
 */

	struct rcoef *
expand_constraint (

struct constraint *	lcp,		/* IN - logical constraint */
struct rcoef *		cp,		/* OUT - coefficient row */
bitmap_t *		edge_mask,	/* IN - set of valid hyperedges */
struct cinfo *		cip		/* IN - compatibility info */
)
{
int			i;
int			j;
int			t;
int			ssize;
int			isize;
int			nedges;
int			kmasks;
struct rcoef *		orig_cp;
bitmap_t *		bp1;
bitmap_t		mask;
int *			vp1;
int *			vp2;

	nedges = cip -> num_edges;

	switch (lcp -> type) {
	case CT_CUTSET:
		orig_cp = cp;
		/* We are given a set F of full sets... */
		for (i = 0; i < nedges; i++) {
			if (NOT BITON (lcp -> mask, i)) continue;
			if (NOT BITON (edge_mask, i)) continue;
			cp -> var = i + RC_VAR_BASE;
			cp -> val = 1.0;
			++cp;
		}
		if (cp <= orig_cp) {
			/* Empty cutset! */
			fatal ("expand_constraint: Bug 1.");
		}
		cp -> var = RC_OP_GE;
		cp -> val = 1.0;
		break;

	case CT_SUBTOUR:
		/* We are given a set S of terminals...  Get size */
		/* of subtour - 1... */
		ssize = -1;
		bp1 = &(lcp -> mask [0]);
		kmasks = cip -> num_vert_masks;
		for (i = 0; i < kmasks; i++) {
			mask = *bp1++;
			ssize += NBITSON (mask);
		}
		for (j = 0; j < nedges; j++) {
			if (NOT BITON (edge_mask, j)) continue;
			vp1 = cip -> edge [j];
			vp2 = cip -> edge [j + 1];
			isize = -1;
			while (vp1 < vp2) {
				t = *vp1++;
				if (BITON (lcp -> mask, t)) {
					++isize;
				}
			}
			if (isize <= 0) continue;
			cp -> var = j + RC_VAR_BASE;
			cp -> val = isize;
			++cp;
		}
		cp -> var = RC_OP_LE;
		cp -> val = ssize;
		break;

	default:
		fatal ("expand_constraint: Bug 2.");
		break;
	}

	return (cp);
}
