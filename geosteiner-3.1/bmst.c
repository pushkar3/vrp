/***********************************************************************

	File:	bmst.c
	Rev:	b-1
	Date:	02/28/2001

	Copyright (c) 1997, 2001 by David M. Warme & Martin Zachariasen

************************************************************************

	Routines to compute Minimum Spanning Trees using the
	Bottleneck Steiner Distance.

************************************************************************

	Modification Log:

	a-1:	09/04/97	warme
		: Created.
	b-1:	02/28/2001	martinz
		: Added new procedures "bmst" and "bmst_terms" that
		: return the edges of the MST (instead of only the
		: length)

************************************************************************/

#include "bsd.h"
#include "steiner.h"


/*
 * Global Routines
 */

dist_t			bmst_length (struct pset *, struct bsd *);
int			bmst (struct pset *, struct bsd *, struct edge *);
dist_t			bmst_terms_length (int *, int, struct bsd *);
int			bmst_terms (int *, int, struct bsd *, struct edge *);


/*
 * External References
 */

	/* none */

/*
 * This routine computes the length of a Minimum Spanning Tree for
 * the given point set as measured using the Bottleneck Steiner Distance.
 */

	dist_t
bmst_length (

struct pset *		pts,		/* IN - point set */
struct bsd *		bsdp		/* IN - BSD data structure */
)
{
int		i;
int		nedges;
struct edge *	ep;
struct edge *	edges;
dist_t		total;

	edges = NEWA (pts -> n - 1, struct edge);

	nedges = bmst (pts, bsdp, &edges [0]);

	if (nedges NE pts -> n - 1) {
		fatal ("bmst_length: bug 1.");
	}

	total = 0;
	ep = edges;
	for (i = 0; i < nedges; i++, ep++) {
		total += ep -> len;
	}

	free ((char *) edges);

	return (total);
}

/*
 * This routine computes a Minimum Spanning Tree for
 * the given point set as measured using the Bottleneck Steiner Distance.
 */

	int
bmst (

struct pset *		pts,		/* IN - point set */
struct bsd *		bsdp,		/* IN - BSD data structure */
struct edge *		edges		/* OUT - edge list */
)
{
int		i;
int		j;
int		n;
int		nedges;
int		mst_edge_count;
struct point *	p1;
struct point *	p2;
struct edge *	ep;
struct edge *	edge_array;

	n	= pts -> n;
	nedges	= (n * (n - 1)) >> 1;

	edge_array = NEWA (nedges, struct edge);

	ep = edge_array;
	p1 = &(pts -> a [0]);
	for (i = 0; i < n; i++, p1++) {
		p2 = p1 + 1;
		for (j = i + 1; j < n; j++, p2++) {
			ep -> len	= bsd (bsdp, p1 -> pnum, p2 -> pnum);
			ep -> p1	= ((pnum_t) i);
			ep -> p2	= ((pnum_t) j);
			++ep;
		}
	}

	mst_edge_count = mst_edge_list (n, nedges, &edge_array [0], edges);

	free ((char *) edge_array);

	return (mst_edge_count);
}

/*
 * This routine computes the length of a Minimum Spanning Tree for
 * the given list of terminals as measured using the
 * Bottleneck Steiner Distance.
 */

	dist_t
bmst_terms_length (

int *			terms,		/* IN - array of terminal numbers */
int			n,		/* IN - number of terminals */
struct bsd *		bsdp		/* IN - BSD data structure */
)
{
int		i;
int		nedges;
struct edge *	ep;
struct edge *	edges;
dist_t		total;

	edges = NEWA (n - 1, struct edge);

	nedges = bmst_terms (terms, n, bsdp, edges);

	if (nedges NE n - 1) {
		fatal ("bmst_terms_length: bug 1.");
	}

	total = 0;
	ep = edges;
	for (i = 0; i < nedges; i++, ep++) {
		total += ep -> len;
	}

	free ((char *) edges);

	return (total);
}

/*
 * This routine computes a Minimum Spanning Tree for
 * the given list of terminals as measured using the
 * Bottleneck Steiner Distance.
 */

	int
bmst_terms (

int *			terms,		/* IN - array of terminal numbers */
int			n,		/* IN - number of terminals */
struct bsd *		bsdp,		/* IN - BSD data structure */
struct edge *		edges		/* OUT - edge list */
)
{
int		i;
int		j;
int		t;
int		nedges;
int		mst_edge_count;
int *		ip;
int *		jp;
struct edge *	ep;
struct edge *	edge_array;

	nedges	= (n * (n - 1)) >> 1;

	edge_array = NEWA (nedges, struct edge);

	ep = edge_array;
	ip = &terms [0];
	for (i = 0; i < n; i++, ip++) {
		t = *ip;
		jp = ip + 1;
		for (j = i + 1; j < n; j++, jp++) {
			ep -> len	= bsd (bsdp, t, *jp);
			ep -> p1	= ((pnum_t) i);
			ep -> p2	= ((pnum_t) j);
			++ep;
		}
	}

	mst_edge_count = mst_edge_list (n, nedges, &edge_array [0], edges);

	free ((char *) edge_array);

	return (mst_edge_count);
}
