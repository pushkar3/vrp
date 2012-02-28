/***********************************************************************

	File:	sll.c
	Rev:	b-1
	Date:	01/22/2000

	Copyright (c) 1998, 2001 by Martin Zachariasen

************************************************************************

	Implementation of a variant of the Smith, Lee and Liebman
	heuristic for Euclidean Steiner trees (Networks 11, 1981)

************************************************************************

	Modification Log:

	a-1:	11/10/98	martinz
		: Created.  Derived from C++ geosteiner96 program
	b-1:	01/22/2000	martinz
		: Heuristic now returns a non-infinity result when
		:  all points are co-linear.
		: Split off elementary geometric functions to efuncs.h.
		: Translate instance mean of points is at origin to
		:  maximize precision of eq-points and Steiner points.
		: Use BSD information if it is available.
		: Use dist_t instead of double.

************************************************************************/

#include "bsd.h"
#include "dsuf.h"
#include "steiner.h"
#include "efuncs.h"

#define ANSI_DECLARATORS
#define REAL coord_t
#include "triangle.h"


/*
 * Global Routines
 */

dist_t		smith_lee_liebman (struct pset *, struct bsd *);


/*
 * Local Types
 */

struct sll_edge {
	dist_t		len;
	dist_t		bsd;
	pnum_t		p1;
	pnum_t		p2;
	bool		mst;
};

struct sll_fst {
	dist_t		len;
	dist_t		ratio;
	pnum_t		p [4];
};


/*
 * Local Routines
 */

static void		add_fst (struct pset *		pts,
				 int *			terms,
				 struct sll_fst *	sll_fsts,
				 int *			fst_count,
				 dist_t			mst_l);
static int		comp_edges (const void *, const void *);
static int		comp_fsts  (const void *, const void *);
static dist_t		fst_length (struct pset *, int *);


/*
 * Local Variables
 */

static struct sll_edge *	sll_edges;
static struct sll_fst *		sll_fsts;

/*
 * Compute the length of a shortest full Steiner tree
 * for a set of terminals with up to 4 terminals.
 * Special fast version which only computes length.
 * Returns 0.0 if no Steiner tree exists.
 * Assume terminals are ordered as they appear on Steiner polygon.
 */

	static
	dist_t
fst_length (

struct pset *		pts,	/* IN - terminal list */
int *			terms	/* IN - indices of terminals to consider */
)
{
int			i;
int			term_count;
struct point *		a;
struct point *		b;
struct point *		c;
struct point *		d;
struct point		e, ctr, e_ad, e_cb;
dist_t			l;
dist_t			min_length;

	/* What is the number of terminals? (should be four or less) */
	term_count = 0;
	while ((term_count < 4) AND (terms [term_count] >= 0)) {
		++term_count;
	}

	if (term_count EQ 2) {
		return (EDIST (&(pts -> a [terms [0]]),
			       &(pts -> a [terms [1]])));
	}

	if (term_count EQ 3) {
		a = &(pts -> a [terms [0]]);
		b = &(pts -> a [terms [1]]);
		c = &(pts -> a [terms [2]]);

		eq_point (a, b, &e);
		eq_circle_center (a, b, &e, &ctr);

		if (right_turn (&e, a, c) AND
		    left_turn (&e, b, c) AND
		    (sqr_dist (&ctr, c) > sqr_dist (&ctr, a))) {
			return EDIST (&e, c);
		}
	}

	if (term_count EQ 4) {
		min_length = INF_DISTANCE;
		for (i = 0; i <= 1; i++) {
			/* Using Lemma 5.2 p. 64 in Hwang, Richards, Winter */
			a = &(pts -> a [terms [i]]);
			d = &(pts -> a [terms [i+1]]);
			c = &(pts -> a [terms [i+2]]);
			b = &(pts -> a [terms [(i+3) % 4]]);

			/* Find intersetion point between ac and bd.
			   It is the same in both iterations */

			if ((i EQ 0) AND
			    NOT segment_intersection (a, c, b, d, &ctr)) break;

			eq_point (a, d, &e_ad);
			eq_point (c, b, &e_cb);
			if (NOT wedge120 (d, a, &e_cb) AND
			    NOT wedge120 (&e_cb, d, a) AND
			    NOT wedge120 (&e_ad, b, c) AND
			    NOT wedge120 (b, c, &e_ad) AND
			    NOT wedge120 (a, &ctr, d)) {
				l = EDIST (&e_ad, &e_cb);
				if (l < min_length) {
					min_length = l;
				}
			}
		}
		if (min_length < INF_DISTANCE) return (min_length);
	}

	return (0.0);
}

/*
 * For sorting edges by length (should be changed!!)
 */

	static
	int
comp_edges (

const void *		p1,
const void *		p2
)
{
dist_t			l1;
dist_t			l2;

	l1 = sll_edges [ *(int*) p1].len;
	l2 = sll_edges [ *(int*) p2].len;

	if (l1 < l2) return (-1);
	if (l1 > l2) return (1);
	return (0);
}


/*
 * For sorting FSTs by ratio (should be changed!!)
 */

	static
	int
comp_fsts (

const void *		p1,
const void *		p2
)
{
dist_t			l1;
dist_t			l2;

	l1 = sll_fsts [ *(int*) p1].ratio;
	l2 = sll_fsts [ *(int*) p2].ratio;

	if (l1 < l2) return (-1);
	if (l1 > l2) return (1);
	return (0);
}

/*
 * Add generated FST to list of FSTs
 */

	static
	void
add_fst (

struct pset *	pts,		/* IN - terminal list */
int *		terms,		/* IN - indices of terminals in FST */
struct sll_fst* sll_fsts,	/* IN/OUT - list of FSTs */
int *		fst_count,	/* IN/OUT - current FST count */
dist_t		mst_l		/* Length of MST spanning terminals in FST */
)
{
int		j;
dist_t		smt_l;
dist_t		ratio;

	smt_l = fst_length (pts, terms);

	if (smt_l > 0) {
		ratio = smt_l / mst_l;
		if (ratio < 1.0) {
			for (j = 0; j < 4; j++) {
				sll_fsts [*fst_count].p [j] = terms [j];
			}
			sll_fsts [*fst_count].len   = smt_l;
			sll_fsts [*fst_count].ratio = ratio;
			++(*fst_count);
		}
	}
}

/*
 * Variant of heuristic by Smith, Lee and Liebman (Networks 11, 1981)
 * All 4-terminal candidates with three connected MST-edges
 * are put on the priority queue
 */

	dist_t
smith_lee_liebman (

struct pset *	ptss,	/* IN - terminals list for which
				heuristic tree should be constructed */
struct bsd *	bsdp	/* IN - BSD data structure */
)
{
int			i, j, jj, ei, fi, fst_count, mst_count;
int			ni, nj, nt, np1, np2, p1, p2, p4;
bool			neighbour_edge_mst;
bool			nb_edge;
bool			convex_region;
bool			in_same_block;
int			neighbour_edge;
int			neighbour_edge_idx;
int			outgoing_mst_edge;
int			term_count;
int			root [4];
int			terms [4];
int *			triangleedges;
int *			sortededges;
int *			sortedfsts;
dist_t			mst_l, total_mst_l, outgoing_mst_edge_l;
dist_t			total_smt_l, mx, my;
struct point		minp, maxp;
struct dsuf		mstsets, fstsets;
struct triangulateio	in, out, vorout;
struct pset *		pts;

	/* Special cases */

	if (ptss -> n <= 1) return (0.0);
	if (ptss -> n EQ 2) return (EDIST (&(ptss -> a [0]),
					   &(ptss -> a [1])));

	/* Compute mean of all terminals and translate in order */
	/* to improve the precision of computed eq-points.	*/
	mx = 0.0;
	my = 0.0;
	for (i = 0; i < ptss -> n; i++) {
		mx += ptss -> a[i].x;
		my += ptss -> a[i].y;
	}
	mx = floor(mx / ((dist_t) ptss -> n));
	my = floor(my / ((dist_t) ptss -> n));

	/* Transfer to triangulate data structure and call triangle */

	in.numberofpoints = ptss -> n;
	in.pointlist = NEWA (2*in.numberofpoints, REAL);

	pts = NEW_PSET(ptss -> n);
	pts -> n = ptss -> n;
	for (i = 0; i < pts -> n; i++) {
		pts -> a[i].x	 = ptss -> a[i].x - mx;
		pts -> a[i].y	 = ptss -> a[i].y - my;
		pts -> a[i].pnum = ptss -> a[i].pnum;
		in.pointlist [2*i  ] = pts -> a [i].x;
		in.pointlist [2*i+1] = pts -> a [i].y;
	}

	in.numberofpointattributes	= 0;
	in.pointattributelist		= NULL;
	in.pointmarkerlist		= NULL;
	in.numberofsegments		= 0;
	in.numberofholes		= 0;
	in.numberofregions		= 0;
	in.regionlist			= NULL;
	out.pointlist			= NULL;
	out.trianglelist		= NULL;
	out.neighborlist		= NULL;

	triangulate ("znNBQ", &in, &out, &vorout);

	if (out.numberoftriangles EQ 0) {

		/* There are no triangles (all input points co-linear).
		   Compute SMT as straight line segment between points */
		minp.x = INF_DISTANCE; maxp.x = -INF_DISTANCE;
		minp.y = INF_DISTANCE; maxp.y = -INF_DISTANCE;
		for (i = 0; i < pts -> n; i++) {
			minp.x = MIN (minp.x, pts -> a [i].x);
			maxp.x = MAX (maxp.x, pts -> a [i].x);
			minp.y = MIN (minp.y, pts -> a [i].y);
			maxp.y = MAX (maxp.y, pts -> a [i].y);
		}

		free (in.pointlist);
		free (out.pointlist);
		free (out.trianglelist);
		free (out.neighborlist);
		return (EDIST (&minp, &maxp));
	}

	/* Construct edge information */
	sll_edges	= NEWA (out.numberofedges, struct sll_edge);
	triangleedges	= NEWA (3*out.numberoftriangles, int);
	for (i = 0; i < 3*out.numberoftriangles; i++) {
		triangleedges [i] = -1;
	}

	ei = 0;
	for (i = 0; i < out.numberoftriangles; i++) {
		for (j = 0; j < 3; j++) {
			p1 = out.trianglelist [3*i + j];
			p2 = out.trianglelist [3*i + ((j+1) % 3)];
			if (triangleedges [3*i + j] EQ -1) {
				/* only add once */
				sll_edges [ei].len = EDIST (&(pts -> a [p1]),
							    &(pts -> a [p2]));
				/* Get BSD info if it is there */
				if ((bsdp NE NULL) AND
				    (pts -> a [p1].pnum >= 0) AND
				    (pts -> a [p2].pnum >= 0)) {
					sll_edges [ei].bsd = bsd (bsdp,
								  pts -> a [p1].pnum,
								  pts -> a [p2].pnum);
				}
				else {
					sll_edges [ei].bsd = sll_edges [ei].len;
				}

				sll_edges [ei].p1  = p1;
				sll_edges [ei].p2  = p2;
				sll_edges [ei].mst = FALSE;
				triangleedges [3*i + j] = ei;

				/* Now go through neighbouring triangles */
				/* and add information about the new edge */
				for (ni = 0; ni < 3; ni++) {
					nt = out.neighborlist [3*i + ni];
					if (nt EQ -1) continue;
					for (nj = 0; nj < 3; nj++) {
						np1 = out.trianglelist [3*nt + nj];
						np2 = out.trianglelist [3*nt + ((nj+1) % 3)];
						if ((np1 EQ p2) AND
						    (np2 EQ p1)) {
							 /* found reverse edge */
							triangleedges [3*nt + nj] = ei;
						}
					}
				}
				++ei;
			}
		}
	}

	/* Sort edges */

	sortededges = NEWA (out.numberofedges, int);
	for (i = 0; i < out.numberofedges; i++) {
		sortededges [i] = i;
	}
	qsort (sortededges, out.numberofedges, sizeof(int), comp_edges);

	/* Use Kruskal to find MST */

	dsuf_create (&mstsets, pts -> n);
	for (i = 0; i < pts -> n; i++) {
		dsuf_makeset (&mstsets, i);
	}
	total_mst_l = 0.0;
	for (i = 0; i < out.numberofedges; i++) {
		/* go through edges in sorted order */
		ei = sortededges [i];
		root [1] = dsuf_find (&mstsets, sll_edges [ei].p1);
		root [2] = dsuf_find (&mstsets, sll_edges [ei].p2);

		if (root [1] NE root [2]) {
			dsuf_unite (&mstsets, root [1], root [2]);
			total_mst_l += sll_edges [ei].len;
			sll_edges [ei].mst = TRUE; /* remember this edge */
		}
	}

	/* Generate FSTs with 3 and 4 terminals */

	sll_fsts = NEWA (4*out.numberoftriangles, struct sll_fst);
	fst_count = 0;
	for (i = 0; i < out.numberoftriangles; i++) {
		/* Count the number of MST edges and find their length sum */
		mst_count = 0;
		mst_l = 0.0;
		for (j = 0; j < 3; j++) {
			ei = triangleedges [3*i + j];
			if (sll_edges [ei].mst) {
				++mst_count;
				mst_l += sll_edges [ei].len;
			}
		}

		/* If there are two MST edges try to construct an FST */
		if (mst_count NE 2) continue;

		for (j = 0; j < 3; j++) {
			terms [j] = out.trianglelist [3*i + j];
		}
		terms [3] = -1;

		/* Add 3-terminal FST */
		add_fst (pts, terms, sll_fsts, &fst_count, mst_l);

		/* Now go through neighbouring triangles and add */
		/* valid FSTs */
		for (ni = 0; ni < 3; ni++) {
			 /* get triangle index */
			nt = out.neighborlist [3*i + ni];
			if (nt EQ -1) continue;

			/* Find neighbouring edge and outgoing MST */
			/* edge (note that there can be at most one) */
			neighbour_edge		= -1;
			neighbour_edge_idx	= -1;
			outgoing_mst_edge	= -1;
			neighbour_edge_mst	= FALSE;
			outgoing_mst_edge_l	= 0.0;

			for (nj = 0; nj < 3; nj++) {
				ei = triangleedges [3*nt + nj];

				/* Is this the neighouring edge? */
				nb_edge = FALSE;
				for (j = 0; j < 3; j++) {
					if (triangleedges [3*i + j] EQ ei) {
						nb_edge = TRUE;
					}
				}
				if (nb_edge) {
					if (sll_edges [ei].mst) {
						/* neighbour edge is an */
						/* MST edge */
						neighbour_edge_mst = TRUE;
					}
					neighbour_edge	   = ei;
					neighbour_edge_idx = nj;
				}
				else if (sll_edges [ei].mst) {
					/* outgoing MST edge identified */
					outgoing_mst_edge = ei;
					outgoing_mst_edge_l
						= sll_edges [ei].len;
				}
			}

			if (outgoing_mst_edge < 0) continue;

			if (neighbour_edge_mst AND (i >= nt)) continue;
			/* only consider once... */

			/* Idenfify fourth terminal */

			p4 = out.trianglelist [3*nt + ((neighbour_edge_idx + 2) % 3)];

			/* Make list of points on border of triangles */

			j = -1;
			jj = -1;
			do {
				terms [++jj] = out.trianglelist [3*i + (++j)];
			} while (triangleedges [3*i + j] NE neighbour_edge);
			terms [++jj] = p4;
			while (j < 2) {
				terms [++jj] = out.trianglelist [3*i + (++j)];
			}

			/* Make sure that the two triangles make up a */
			/* convex region */

			convex_region = TRUE;
			for (j = 0; j < 4; j++) {
				if (right_turn (&(pts -> a [terms [j]]),
						&(pts -> a [terms [(j+1) % 4]]),
						&(pts -> a [terms [(j+2) % 4]]))) {
					convex_region = FALSE;
					break;
				}
			}

			if (convex_region) {
				/* Add 4-terminal FST */
				add_fst (pts,
					 terms,
					 sll_fsts,
					 &fst_count,
					 mst_l + outgoing_mst_edge_l);
			}
		}
	}

	/* Sort FSTs */

	sortedfsts = NEWA (fst_count, int);
	for (i = 0; i < fst_count; i++) {
		sortedfsts [i] = i;
	}
	qsort (sortedfsts, fst_count, sizeof (int), comp_fsts);

	/* Build heuristic tree */

	dsuf_create (&fstsets, pts -> n);
	for (i = 0; i < pts -> n; i++) {
		dsuf_makeset (&fstsets, i);
	}
	total_smt_l = 0.0;
	for (i = 0; i < fst_count; i++) {
		/* go through FSTs in sorted order */
		fi = sortedfsts [i];
		term_count = 0;
		while ((term_count < 4) AND
		       (sll_fsts [fi].p [term_count] >= 0)) {
			++term_count;
		}

		/* check that no two terminals are in the same block */
		in_same_block = FALSE;
		for (j = 0; j < term_count; j++) {
			root [j] = dsuf_find (&fstsets, sll_fsts [fi].p [j]);
		}
		for (j = 0; j < term_count-1; j++) {
			for (jj = j+1; jj < term_count; jj++) {
				if (root [j] EQ root [jj]) {
					in_same_block = TRUE;
					break;
				}
			}
		}
		if (in_same_block) continue;

		for (j = 1; j < term_count; j++) {
			dsuf_unite (&fstsets,
				    dsuf_find (&fstsets, sll_fsts [fi].p [0]),
				    dsuf_find (&fstsets, sll_fsts [fi].p [j]));
		}
		total_smt_l += sll_fsts [fi].len;
	}

	for (i = 0; i < out.numberofedges; i++) {
		/* go through edges in sorted order */
		ei = sortededges [i];
		root [1] = dsuf_find (&fstsets, sll_edges [ei].p1);
		root [2] = dsuf_find (&fstsets, sll_edges [ei].p2);
		if (root [1] NE root [2]) {
			dsuf_unite (&fstsets, root [1], root [2]);
			/* Use BSD length here */
			total_smt_l += sll_edges [ei].bsd;
		}
	}

	/* Free all allocated arrays, including those allocated by Triangle */

	dsuf_destroy (&mstsets);
	dsuf_destroy (&fstsets);

	free (pts);
	free (in.pointlist);
	free (out.pointlist);
	free (out.trianglelist);
	free (out.neighborlist);
	free (triangleedges);
	free (sll_edges);
	free (sll_fsts);
	free (sortededges);
	free (sortedfsts);

	return (total_smt_l);
}
