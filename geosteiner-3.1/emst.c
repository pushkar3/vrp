/***********************************************************************

	File:	emst.c
	Rev:	b-1
	Date:	02/28/2001

	Copyright (c) 1993, 2001 by Martin Zachariasen & David M. Warme

************************************************************************

	Routines to compute Euclidean Minimum Spanning Trees.

************************************************************************

	Modification Log:

	a-1:	11/10/98	martinz
		: Created.
	b-1:	02/28/2001	warme
		: Renamed build_edges, and eliminated its 3rd "at_least"
		:  parameter.  Build complete graph when N < 10.
		: Use common routines now in mst.c.
		: Fix duplicate points problem.

************************************************************************/

#include "dsuf.h"
#include "steiner.h"

#define ANSI_DECLARATORS
#define REAL coord_t
#include "triangle.h"

/*
 * Global Routines
 */

int			euclidean_mst (struct pset *, struct edge *);
dist_t			euclidean_mst_length (struct pset *);


/*
 * External References
 */

	/* none */


/*
 * Local Routines
 */

static int		build_euclidean_edges (struct pset *,
					       struct edge **);
static int *		heapsort_x (struct pset *);

/*
 * This routine computes the total length of the Euclidean Minimum
 * Spanning Tree of the given set of points.
 */

	dist_t
euclidean_mst_length (

struct pset *		pts	/* IN - point set to find MST length of. */
)
{
int		i;
int		nedges;
dist_t		total;
struct edge *	ep;
struct edge *	edges;

	edges = NEWA (pts -> n - 1, struct edge);

	nedges = euclidean_mst (pts, &edges [0]);

	if (nedges NE pts -> n - 1) {
		fatal ("euclidean_mst_length: Bug 1.");
	}

	total = 0;
	ep = &edges [0];
	for (i = 0; i < nedges; i++) {
		total += ep -> len;
		++ep;
	}

	free ((char *) edges);

	return (total);
}

/*
 * This routine computes an Euclidean Minimum Spanning Tree for the
 * given point set.  The output is a list of edges.
 *
 * The "Triangle" package by Jonathan Richard Shewchuk is used
 * to provide the Delaunay triangulation for the set of points.
 */

	int
euclidean_mst (

struct pset *		pts,		/* IN - point set. */
struct edge *		edges		/* OUT - edge list. */
)
{
int		nedges;
int		mst_edge_count;
struct edge *	edge_array;

	nedges = build_euclidean_edges (pts, &edge_array);

	mst_edge_count = mst_edge_list (pts -> n,
					nedges,
					&edge_array [0],
					edges);

	free ((char *) edge_array);

	return (mst_edge_count);
}

/*
 * This routine builds an edge-list containing all edges in
 * the Delaunay triangulation of the set of points.
 */

	static
	int
build_euclidean_edges (

struct pset *		pts,		/* IN - set of points */
struct edge **		edges_out	/* OUT - edge list */
)
{
int			i, i1;
int			j, j1;
int             	k;
int             	n;
int	 		nedges;
int			ndup;
struct edge *		edges;
struct point *		p1;
struct point *		p2;
int *			order;
int *			orig_vnum;
bool *			dflags;
struct edge *		zedges;
struct triangulateio	in, out, vorout;

	n = pts -> n;

	if (n < 10) {
		/* Build the complete graph... */
		nedges = n * (n - 1) / 2;
		edges = NEWA (nedges, struct edge);
		*edges_out = edges;
		p1 = &(pts -> a [0]);
		for (i = 0; i < n; i++, p1++) {
			p2 = p1 + 1;
			for (j = i + 1; j < n; j++, p2++) {
				edges -> len	= EDIST (p1, p2);
				edges -> p1	= ((pnum_t) i);
				edges -> p2	= ((pnum_t) j);
				++edges;
			}
		}
		return (nedges);
	}

	/* Triangle does a "good" job when given duplicate points -- it	*/
	/* emits edges to only one of the given coincident points.	*/
	/* This is bad for us, however, since we want a fully connected	*/
	/* MST, and cannot get one if there are vertices for which we	*/
	/* have no incident edges.  Therefore we must detect all	*/
	/* duplicate points and manually generate one zero-length edge	*/
	/* for each (copies 2 through K of a point each have an edge to	*/
	/* the first copy of the point).				*/

	order = heapsort_x (pts);

	dflags = NEWA (n, bool);
	memset (dflags, FALSE, n * sizeof (dflags [0]));

	zedges = NEWA (n, struct edge);
	memset (zedges, 0, n * sizeof (zedges));

	ndup = 0;
	for (i = 0; i < n - 1; ) {
		i1 = order [i];
		p1 = &(pts -> a [i1]);
		for (j = i; ; ) {
			++j;
			if (j >= n) break;

			j1 = order [j];
			p2 = &(pts -> a [j1]);
			if (p1 -> x NE p2 -> x) break;
			if (p1 -> y NE p2 -> y) break;

			/* Point j1 is a duplicate of point i1.  The	*/
			/* sort also guarantees that i1 < j1.  Omit	*/
			/* point j1.					*/

			dflags [j1] = TRUE;

			/* Generate zero-length edge (i1,j1). */

			zedges [ndup].len	= 0.0;
			zedges [ndup].p1	= i1;
			zedges [ndup].p2	= j1;
			++ndup;
		}
		i = j;
	}

	free (order);

	/* Get array to map duplicate-free points back to originals. */
	orig_vnum = NEWA (n, int);

	/* Set up data structure to call triangle. */
	in.pointlist	  = NEWA (2 * n, coord_t);
	j = 0;
	for (i = 0; i < n; i++) {
		if (dflags [i]) continue;
		in.pointlist [2*j    ] = pts -> a [i].x;
		in.pointlist [2*j + 1] = pts -> a [i].y;
		orig_vnum [j] = i;
		++j;
	}
	in.numberofpoints = j;

	free (dflags);

	in.numberofpointattributes	= 0;
	in.pointattributelist		= 0;
	in.pointmarkerlist		= 0;
	in.numberofsegments		= 0;
	in.numberofholes		= 0;
	in.numberofregions		= 0;
	in.regionlist			= 0;
	out.pointlist			= 0;
	out.edgelist			= 0;
	out.trianglelist		= 0;
	out.neighborlist		= 0;

	triangulate ("zeNBQ", &in, &out, &vorout);

	free (in.pointlist);
	free (out.pointlist);
	free (out.trianglelist);
	free (out.neighborlist);

	nedges = out.numberofedges;

	edges = NEWA (nedges + ndup, struct edge);
	*edges_out = edges;

	for (k = 0; k < ndup; k++) {
		*edges++ = zedges [k];
	}

	free (zedges);

	for (k = 0; k < nedges; k++) {
		i  = orig_vnum [out.edgelist [2*k    ]];
		j  = orig_vnum [out.edgelist [2*k + 1]];
		p1 = &(pts -> a [i]);
		p2 = &(pts -> a [j]);

		edges -> len	= EDIST (p1, p2);
		edges -> p1	= ((pnum_t) i);
		edges -> p2	= ((pnum_t) j);
		++edges;
	}

	free (out.edgelist);
	free (orig_vnum);

	return (nedges + ndup);
}

/*
 * Use the heapsort algorithm to sort the given terminals in increasing
 * order by the following keys:
 *
 *	1.	X coordinate
 *	2.	Y coordinate
 *	3.	index (i.e., position within input data)
 *
 * Of course, we do not move the points, but rather permute an array
 * of indexes into the points.
 */

	static
	int *
heapsort_x (

struct pset *		pts		/* IN - the terminals to sort */
)
{
int			i, i1, i2, j, k, n;
struct point *		p1;
struct point *		p2;
int *			index;

	n = pts -> n;

	index = NEWA (n, int);
	for (i = 0; i < n; i++) {
		index [i] = i;
	}

	/* Construct the heap via sift-downs, in O(n) time. */
	for (k = n >> 1; k >= 0; k--) {
		j = k;
		for (;;) {
			i = (j << 1) + 1;
			if (i + 1 < n) {
				/* Increment i (to right subchild of j) */
				/* if the right subchild is greater. */
				i1 = index [i];
				i2 = index [i + 1];
				p1 = &(pts -> a [i1]);
				p2 = &(pts -> a [i2]);
				if ((p2 -> x > p1 -> x) OR
				    ((p2 -> x EQ p1 -> x) AND
				     ((p2 -> y > p1 -> y) OR
				      ((p2 -> y EQ p1 -> y) AND
				       (i2 > i1))))) {
					++i;
				}
			}
			if (i >= n) {
				/* Hit bottom of heap, sift-down is done. */
				break;
			}
			i1 = index [j];
			i2 = index [i];
			p1 = &(pts -> a [i1]);
			p2 = &(pts -> a [i2]);
			if ((p1 -> x > p2 -> x) OR
			    ((p1 -> x EQ p2 -> x) AND
			     ((p1 -> y > p2 -> y) OR
			      ((p1 -> y EQ p2 -> y) AND
			       (i1 > i2))))) {
				/* Greatest child is smaller.  Sift-	*/
				/* down is done. */
				break;
			}
			/* Sift down and continue. */
			index [j] = i2;
			index [i] = i1;
			j = i;
		}
	}

	/* Now do actual sorting.  Exchange first/last and sift down. */
	while (n > 1) {
		/* Largest is at index [0], swap with index [n-1],	*/
		/* thereby putting it into final position.		*/
		--n;
		i = index [0];
		index [0] = index [n];
		index [n] = i;

		/* Now restore the heap by sifting index [0] down. */
		j = 0;
		for (;;) {
			i = (j << 1) + 1;
			if (i + 1 < n) {
				/* Increment i (to right subchild of j) */
				/* if the right subchild is greater. */
				i1 = index [i];
				i2 = index [i + 1];
				p1 = &(pts -> a [i1]);
				p2 = &(pts -> a [i2]);
				if ((p2 -> x > p1 -> x) OR
				    ((p2 -> x EQ p1 -> x) AND
				     ((p2 -> y > p1 -> y) OR
				      ((p2 -> y EQ p1 -> y) AND
				       (i2 > i1))))) {
					++i;
				}
			}
			if (i >= n) {
				/* Hit bottom of heap, sift-down is done. */
				break;
			}
			i1 = index [j];
			i2 = index [i];
			p1 = &(pts -> a [i1]);
			p2 = &(pts -> a [i2]);
			if ((p1 -> x > p2 -> x) OR
			    ((p1 -> x EQ p2 -> x) AND
			     ((p1 -> y > p2 -> y) OR
			      ((p1 -> y EQ p2 -> y) AND
			       (i1 > i2))))) {
				/* Greatest child is smaller.  Sift-	*/
				/* down is done. */
				break;
			}
			/* Sift down and continue. */
			index [j] = i2;
			index [i] = i1;
			j = i;
		}
	}

	return (index);
}
