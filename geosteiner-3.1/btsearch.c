/***********************************************************************

	File:	btsearch.c
	Rev:	a-1
	Date:	11/30/2000

	Copyright (c) 1993, 2001 by David M. Warme

************************************************************************

	Routines to perform the backtrack search for a solution.

************************************************************************

	Modification Log:

	a-1:	11/30/2000	warme
		: Created.

************************************************************************/

#include "search.h"
#include "steiner.h"


/*
 * Global Routines
 */

double			backtrack_search (struct cinfo *, bitmap_t *);
void			initialize_btsearch (void);

#if 0
struct full_set *	sort_full_trees (struct full_set *);
#endif


/*
 * External References
 */

	/* none */


/*
 * Local Constants
 */

#define	MAX_TERMINALS	32


/*
 * Local Types
 */

struct sinfo {
	struct cinfo *		cip;

	/* the best solution so far... */
	dist_t			best_length;
	int			best_count;
	bitmap_t		best_solution;
	bool			found;

	/* Current partial solution... */
	bitmap_t		solution;

	/* arrays accessed in stack fashion during backtrack... */
	bitmap_t *		terms_left;
	bitmap_t *		compat;
	bitmap_t *		incompat;

	/* New ones we have to initialize */
	bitmap_t *		edge_vmasks;
	bitmap_t *		vert_emasks;
	bitmap_t *		cmasks;
	bitmap_t *		incmasks;
};


/*
 * Local Routines
 */

static void			allocate_search_stacks (struct sinfo *);
static bool			find_inaccessible_terminal (struct sinfo *,
							    int);
static pnum_t			find_starting_point (struct cinfo *);
static void			free_search_stacks (struct sinfo *);
static void			init_compat_incompat_masks (struct sinfo *);
static void			init_edge_vmasks (struct sinfo *);
static void			init_vert_emasks (struct sinfo *);
static bool			perform_search (struct sinfo *, bitmap_t *);
static void			search_recurse (struct sinfo *,
						int, int, dist_t);


/*
 * Local Variables
 */

static int8u *		one_bits_in_byte [256 + 1];
static int8u		one_bits_in_byte_array [1024];

/*
 * Perform all initialization for the backtrack search that only needs
 * to be done once -- no matter how many times the search is called.
 *
 * Currently all we need to do is initialize the "one-bits-in-byte" table.
 */

	void
initialize_btsearch (void)

{
int		i;
int		j;
int		byte;
int8u *		p;

	p = &one_bits_in_byte_array [0];
	for (i = 0; i < 256; i++) {
		one_bits_in_byte [i] = p;
		byte = i;
		j = 0;
		while (byte > 0) {
			if ((byte & 0x01) NE 0) {
				*p++ = j;
			}
			++j;
			byte >>= 1;
		}
	}

	/* Terminate last table entry. */
	one_bits_in_byte [i] = p;
}

/*
 * This routine performs a backtrack search for a Steiner Minimal Tree of
 * the given point set.  It composes a minimum-length solution from some
 * combination of the given full-sets.
 */

	double
backtrack_search (

struct cinfo *		cip,		/* IN - compatibility info */
bitmap_t *		smt_mask	/* OUT - SMT as a set of full-sets */
)
{
bool			failed;
struct sinfo		sinfo;

	if ((cip -> num_vert_masks NE 1) OR (cip -> num_edge_masks NE 1)) {
		fatal ("backtrack_search: Bug 1.");
	}

	/* First, zero out the result bit-mask... */
	smt_mask [0] = 0;

	/* Create and initialize a search-context structure containing	*/
	/* worst-case memory allocations...				*/
	sinfo.cip = cip;
	allocate_search_stacks (&sinfo);

	init_edge_vmasks (&sinfo);
	init_vert_emasks (&sinfo);
	init_compat_incompat_masks (&sinfo);

	failed = perform_search (&sinfo, smt_mask);

	if (failed) {
		fatal ("backtrack_search: Bug 2.");
	}

	free ((char *) (sinfo.incmasks));
	free ((char *) (sinfo.cmasks));
	free ((char *) (sinfo.vert_emasks));
	free ((char *) (sinfo.edge_vmasks));

	/* Free up the search stacks. */
	free_search_stacks (&sinfo);

	return (sinfo.best_length);
}

/*
 * Initialize an array containing one mask per edge.  The mask indicates
 * which vertices the edge contains.
 */

	static
	void
init_edge_vmasks (

struct sinfo *	sip		/* IN/OUT - search/compatibility info */
)
{
int			i;
int			j;
int			nedges;
struct cinfo *		cip;
int *			vp1;
int *			vp2;
bitmap_t *		array;
bitmap_t		mask;

	cip = sip -> cip;

	nedges	= cip -> num_edges;

	array = NEWA (nedges, bitmap_t);

	for (i = 0; i < nedges; i++) {
		mask = 0;
		vp1 = cip -> edge [i];
		vp2 = cip -> edge [i + 1];
		while (vp1 < vp2) {
			j = *vp1++;
			mask |= (1 << j);
		}
		array [i] = mask;
	}

	sip -> edge_vmasks = array;
}

/*
 * Initialize an array containing one mask per vertex.  The mask indicates
 * which edges the vertex is incident to.
 */

	static
	void
init_vert_emasks (

struct sinfo *	sip		/* IN/OUT - search/compatibility info */
)
{
int			i;
int			j;
int			nverts;
struct cinfo *		cip;
int *			ep1;
int *			ep2;
bitmap_t *		array;
bitmap_t		mask;

	cip = sip -> cip;

	nverts	= cip -> num_verts;

	array = NEWA (nverts, bitmap_t);

	for (i = 0; i < nverts; i++) {
		mask = 0;
		ep1 = cip -> term_trees [i];
		ep2 = cip -> term_trees [i + 1];
		while (ep1 < ep2) {
			j = *ep1++;
			mask |= (1 << j);
		}
		array [i] = mask;
	}

	sip -> vert_emasks = array;
}

/*
 * Initialize an array containing one mask per edge.  The mask indicates
 * which edges are COMPATIBLE with the given edge.  In this case we
 * assume two edges are compatible if they share exactly 1 vertex.
 */

	static
	void
init_compat_incompat_masks (

struct sinfo *	sip		/* IN/OUT - search/compatibility info */
)
{
int			j;
int			k;
int			e1;
int			e2;
int			nedges;
struct cinfo *		cip;
int *			vp1;
int *			vp2;
int *			ep1;
int *			ep2;
bitmap_t *		carray;
bitmap_t *		iarray;
bitmap_t		edges_seen;
bitmap_t		e1_vmask;
bitmap_t		cmask;
bitmap_t		imask;
bitmap_t		mask;
bitmap_t *		vmasks;

	cip = sip -> cip;

	nedges	= cip -> num_edges;

	carray = NEWA (nedges, bitmap_t);
	iarray = NEWA (nedges, bitmap_t);

	vmasks = sip -> edge_vmasks;

	for (e1 = 0; e1 < nedges; e1++) {
		imask = 0;
		cmask = 0;
		edges_seen = 0;
		e1_vmask = vmasks [e1];
		vp1 = cip -> edge [e1];
		vp2 = cip -> edge [e1 + 1];
		while (vp1 < vp2) {
			j = *vp1++;
			ep1 = cip -> term_trees [j];
			ep2 = cip -> term_trees [j + 1];
			while (ep1 < ep2) {
				e2 = *ep1++;
				if ((edges_seen & (1 << e2)) NE 0) continue;
				edges_seen |= (1 << e2);
				mask = e1_vmask & vmasks [e2];
				k = NBITSON (mask);
				if (k EQ 1) {
					cmask |= (1 << e2);
				}
				else {
					/* k should be > 1 here! */
					imask |= (1 << e2);
				}
			}
		}
		carray [e1] = cmask;
		iarray [e1] = imask;
	}

	sip -> cmasks	= carray;
	sip -> incmasks	= iarray;
}

/*
 * This routine pre-sorts the list of full-sets and re-numbers them so that
 * at each node in the search tree, we select candidate sub-trees in an
 * order that heuristically puts the most likely solutions first.  This has
 * a tendency to reduce the search space by producing earlier length cutoffs
 * on the later sub-trees.
 */

#if 0

	struct full_set *
sort_full_trees (

struct full_set *	fsp		/* IN - full-trees to sort. */
)
{
int			i;
int			n;
int			h;
struct full_set *	tmp;
struct full_set **	trees;
struct full_set **	endp;
struct full_set **	p1;
struct full_set **	p2;
struct full_set **	p3;
struct full_set **	p4;
dist_t			dist_1;
dist_t			dist_2;
int			n_1;
int			n_2;
struct full_set *	root;
struct full_set **	hookp;

	trees = put_trees_in_array (fsp, &n);

	endp = trees + n;

	for (h = 1; h <= n; h = 3*h+1) {
	}

	do {
		h = h / 3;
		p4 = trees + h;
		p1 = p4;
		while (p1 < endp) {
			tmp = *p1;
			dist_1 = tmp -> tree_len;
			n_1	= tmp -> terminals -> n - 1;
			p2 = p1;
			while (TRUE) {
				p3 = (p2 - h);
				dist_2	= (*p3) -> tree_len;
				n_2	= (*p3) -> terminals -> n - 1;
				if ((dist_2 * n_1) <= (dist_1 * n_2)) break;
				*p2 = *p3;
				p2 = p3;
				if (p2 < p4) break;
			}
			*p2 = tmp;
			++p1;
		}
	} while (h > 1);

	i = 0;
	root = NULL;
	hookp = &root;
	p1 = trees;
	while (p1 < endp) {
		tmp = *p1++;
		tmp -> tree_num = i++;
		tmp -> next	= NULL;
		*hookp = tmp;
		hookp = &(tmp -> next);
	}

	return (root);
}

#endif

/*
 * This routine handles the top-most level of the search.
 */

	static
	bool
perform_search (

struct sinfo *		sip,		/* IN - search/compatibility info */
bitmap_t *		smt_mask	/* OUT - solution */
)
{
int			nterms;
bitmap_t		omit;
bitmap_t		mask;
pnum_t			first_point;
struct cinfo *		cip;
int			t;
int *			tp1;
int *			endp;
bool			failed;

	cip = sip -> cip;

	/* Set up initial state of best solution so far... */
	sip -> best_length	= INF_DISTANCE;
	sip -> best_count	= 0;
	sip -> best_solution	= 0;
	sip -> found		= FALSE;

	sip -> terms_left [0] = cip -> initial_vert_mask [0];

	mask = sip -> terms_left [0];
	nterms = NBITSON (mask);

	sip -> compat [0]	= 0;
	sip -> incompat [0]	= 0;
	sip -> solution		= 0;

	/* Do the search once for each full-tree that contains the */
	/* starting point...					   */
	first_point = find_starting_point (cip);

	omit = 0;

	tp1	= cip -> term_trees [first_point];
	endp	= cip -> term_trees [first_point + 1];
	while (tp1 < endp) {
		t = *tp1++;

		/* Skip trees that are not part of this component... */
		if (NOT BITON (cip -> initial_edge_mask, t)) continue;

		sip -> solution |= (1 << t);

		sip -> terms_left [1] =
			sip -> terms_left [0] & ~(sip -> edge_vmasks [t]);

		sip -> compat [1] = sip -> cmasks [t];

		sip -> incompat [1] =   sip -> incmasks [t]
				      | omit
				      | ~ (cip -> initial_edge_mask [0]);

		search_recurse (sip,
				1,
				nterms - cip -> edge_size [t],
				cip -> cost [t]);

		sip -> solution &= ~(1 << t);
		omit |= (1 << t);
	}

	failed = NOT (sip -> found);

	smt_mask [0] = sip -> best_solution;

	return (failed);
}

/*
 * This routine allocates all of the search stacks needed.
 */

	static
	void
allocate_search_stacks (

struct sinfo *		sip		/* IN - search/compatibility info */
)
{
int		nedges;
int		n;
struct cinfo *	cip;

	cip	= sip -> cip;
	nedges	= cip -> num_edges;

	n = (nedges + 1);

	sip -> terms_left	= NEWA (n, bitmap_t);
	sip -> compat		= NEWA (n, bitmap_t);
	sip -> incompat		= NEWA (n, bitmap_t);
}

/*
 * This routine frees up the search stacks.
 */

	static
	void
free_search_stacks (

struct sinfo *		sip		/* IN - search info. */
)
{
	free ((char *) sip -> incompat);
	free ((char *) sip -> compat);
	free ((char *) sip -> terms_left);

	sip -> terms_left	= NULL;
	sip -> compat		= NULL;
	sip -> incompat		= NULL;
}

/*
 * This routine selects a single terminal to use as a starting point in
 * the search.  We choose a point that is a member of the least number
 * of full-sets.
 */

	static
	pnum_t
find_starting_point (

struct cinfo *		cip		/* IN - compatibility info */
)
{
int			i;
int			j;
int			nverts;
int			nedges;
int			cnt;
int			fewest;
pnum_t			starting;
int *			vp1;
int *			vp2;
int32u			counts [MAX_TERMINALS];

	nverts = cip -> num_verts;
	nedges = cip -> num_edges;

	(void) memset (counts, 0, nverts * sizeof (counts [0]));

	for (i = 0; i < nedges; i++) {
		if (NOT BITON (cip -> initial_edge_mask, i)) continue;

		vp1 = cip -> edge [i];
		vp2 = cip -> edge [i + 1];
		while (vp1 < vp2) {
			j = *vp1++;
			++(counts [j]);
		}
	}

	starting = -1;
	fewest	 = INT_MAX;
	for (i = 0; i < nverts; i++) {
		cnt = counts [i];
		if (cnt <= 0) {
			/* Ignore terminals that are not part	*/
			/* of this sub-problem...		*/
		}
		else if (cnt < fewest) {
			starting = i;
			fewest = cnt;
		}
	}

	if (starting < 0) {
		fatal ("find_starting_point: Bug 1.");
	}

	return (starting);
}

/*
 * This routine performs the recursive portion of the backtrack search.
 */

	static
	void
search_recurse (

struct sinfo *	sip,		/* IN/OUT - search/compatibility info */
int		level,		/* IN - recursion level */
int		nleft,		/* IN - number of terminals left to connect */
dist_t		length		/* IN - length of current partial soln */
)
{
int			k;
bitmap_t *		terms_left;
bitmap_t *		compat;
bitmap_t *		incompat;
bitmap_t		omit;
bitmap_t		mask;
int			i1;
bitmap_t		word;
int			bitno;
int8u *			p1;
int8u *			p2;
int			byte;
dist_t			budget;
struct cinfo *		cip;

	if (length >= sip -> best_length) {
		/* Current partial solution is already too long	*/
		/* to ever become the best seen so far...	*/
		return;
	}

	cip = sip -> cip;

	if (nleft <= 0) {
		if (nleft < 0) {
			fatal ("search_recurse: Bug 1.");
		}

		/* All terminals have been connected!  We have	*/
		/* a new BEST solution!  Save it off.		*/

		sip -> best_solution = sip -> solution;
		sip -> best_length	= length;
		sip -> found		= TRUE;
		return;
	}

	if ((nleft >= 10) AND
	    find_inaccessible_terminal (sip, level)) {
		/* There is a terminal still be be hooked up such that	*/
		/* none of its full-trees is compatible with the	*/
		/* current partial solution!  Early cutoff!		*/
		return;
	}

	terms_left	= &(sip -> terms_left [level]);

	compat		= &(sip -> compat [level]);
	incompat	= &(sip -> incompat [level]);

#if 0 /* Disabled so that we don't need to sort the edges. */

	/* Do the length/n ratio cutoff test -- if none of the feasible	*/
	/* pieces has a length/n ratio that beats budget/nleft, then	*/
	/* there is now way to get a better solution!			*/
	if ((nleft >= 4) AND (sip -> best_length < INF_DISTANCE)) {

		/* Compute the "length budget" with which we must	*/
		/* connect all remaining points...			*/
		budget = sip -> best_length - length;

		if (nleft <= 0) goto no_ratio_cutoff;

		/* Find piece with the lowest length/terms ratio that	*/
		/* is not incompatible with current partial solution...	*/
		word = ~(incompat [0]);
		bitno = 0;
		byte = 0;
		if (word NE 0) {
			do {
				byte = (word & 0xFF);
				if (byte NE 0) break;
				bitno += 8;
				word >>= 8;
			} while (word NE 0);
		}
		if (byte NE 0) {
			k = bitno + one_bits_in_byte [byte] [0];
			if (k < cip -> num_edges) {
				fsp = cip -> full_trees [k];
				if ((cip -> cost [k] * nleft) >=
				    (budget * (cip -> edge_size [k] - 1))) {
					/* Even the piece with the best	*/
					/* length/n ratio does not beat	*/
					/* whats required for a better	*/
					/* solution...  Cutoff!		*/
#if 0
					tracef ("%% RATIO CUTOFF!\n");
#endif
					return;
				}
			}
		}
		no_ratio_cutoff:	;
	}
#endif

	/* Determine the set of trees we will try to add to the	*/
	/* partial solution.					*/

	word = compat [0] & ~incompat [0];

	omit = 0;

	/* Loop over each feasible full-tree. */

	for (bitno = 0; word NE 0; bitno += 8, word >>= 8) {
		byte = (word & 0xFF);
		if (byte EQ 0) continue;
		p1 = one_bits_in_byte [byte];
		p2 = one_bits_in_byte [byte + 1];
		while (p1 < p2) {
			k = bitno + *p1++;
			/* Full-set K is feasible.  Test to see	*/
			/* if it is truly compatible...		*/

			mask = sip -> edge_vmasks [k] & ~terms_left [0];
			i1 = NBITSON (mask);
			if (i1 NE 1) continue;

			/* Tree K intersects the current partial */
			/* solution at exactly one point!	 */

			compat [1] = compat [0] | sip -> cmasks [k];

			incompat [1] =	  incompat [0]
					| sip -> incmasks [k]
					| omit;

			terms_left [1] =
				terms_left [0] & ~ sip -> edge_vmasks [k];

			sip -> solution |= (1 << k);

			search_recurse (sip,
					level + 1,
					nleft - cip -> edge_size [k] + 1,
					length + cip -> cost [k]);

			sip -> solution &= ~(1 << k);
			omit |= (1 << k);
		}
	}
}

/*
 * This routine checks the given node level to see if there are any
 * inaccessible terminals.
 */

	static
	bool
find_inaccessible_terminal (

struct sinfo *	sip,		/* IN - search/compatibility info */
int		level		/* IN - recursion level */
)
{
int		t;
int		bitno;
int		byte;
bitmap_t	cur_incompat;
bitmap_t	word;
bitmap_t	mask;
int8u *		p1;
int8u *		p2;
struct cinfo *	cip;

	cip = sip -> cip;

	cur_incompat = sip -> incompat [level];

	word = sip -> terms_left [level];

	for (bitno = 0; word NE 0; bitno += 8, word >>= 8) {
		byte = (word & 0xFF);
		if (byte EQ 0) continue;

		p1 = one_bits_in_byte [byte];
		p2 = one_bits_in_byte [byte + 1];
		while (p1 < p2) {
			t = bitno + *p1++;

			mask = sip -> vert_emasks [t];
			if (mask EQ (mask & cur_incompat)) {
				/* EVERY edge containing T has been	*/
				/* found to be incompatible with >= 1	*/
				/* edge in the current partial tree	*/
				/* solution!				*/
				return (TRUE);
			}
		}
	}

	return (FALSE);
}
