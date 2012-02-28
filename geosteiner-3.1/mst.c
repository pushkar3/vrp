/***********************************************************************

	File:	mst.c
	Rev:	a-1
	Date:	02/28/2001

	Copyright (c) 1993, 2001 by David M. Warme

************************************************************************

	Routines to compute Minimum Spanning Trees.

************************************************************************

	Modification Log:

	a-1:	02/28/2001	warme
		: Created.  Gathered all copies of these routines
		:  into this one file.

************************************************************************/

#include "dsuf.h"
#include "steiner.h"


/*
 * Global Routines
 */

int			mst_edge_list (int,
				       int,
				       struct edge *,
				       struct edge *);
void			sort_edge_list (struct edge *, int);


/*
 * Local Routines
 */

	/* none */

/*
 * This routine computes the MST of a given list of edges.
 */

	int
mst_edge_list (

int			n,		/* IN - number of vertices */
int			nedges,		/* IN - number of edges */
struct edge *		edge_list,	/* IN - list of edges */
struct edge *		edges		/* OUT - MST edge list */
)
{
int		i;
int		mst_edge_count;
int		components;
int		max_vert;
struct edge *	ep;
struct edge *	ep_endp;
int		root1;
int		root2;
struct dsuf	sets;

	sort_edge_list (edge_list, nedges);

	/* Don't assume that the vertex numbers are well-behaved,	*/
	/* except that they must be non-negative.  We do a quick scan	*/
	/* to determine the largest vertex number and then allocate	*/
	/* a union-find data structure large enough to handle it.  Note	*/
	/* that we then use this union-find data structure in a		*/
	/* completely sparse way -- we only ever access set items for	*/
	/* vertices that are named by an edge.				*/

	max_vert = 1;		/* avoid zero-size union-find... */
	ep = edge_list;
	for (i = 0; i < nedges; i++, ep++) {
		if (ep -> p1 > max_vert) {
			max_vert = ep -> p1;
		}
		if (ep -> p2 > max_vert) {
			max_vert = ep -> p2;
		}
	}

	dsuf_create (&sets, max_vert + 1);

	/* Note that it is not a problem to "makeset" a vertex more	*/
	/* than once...							*/
	ep = edge_list;
	for (i = 0; i < nedges; i++, ep++) {
		dsuf_makeset (&sets, ep -> p1);
		dsuf_makeset (&sets, ep -> p2);
	}

	components = n;
	mst_edge_count = 0;
	ep = edge_list;
	ep_endp = (ep + nedges);

	while (components > 1) {
		if (ep >= ep_endp) {
			/* Ran out of edges before MST complete! */
			fatal ("mst_edge_list: Bug 1.");
		}
		root1 = dsuf_find (&sets, ep -> p1);
		root2 = dsuf_find (&sets, ep -> p2);
		if (root1 NE root2) {
			dsuf_unite (&sets, root1, root2);
			*edges = *ep;
			++edges;
			++mst_edge_count;
			--components;
		}
		++ep;
	}

	dsuf_destroy (&sets);

	return (mst_edge_count);
}

/*
 * This routine sorts the given edge list in INCREASING order by edge length.
 */

	void
sort_edge_list (

struct edge *		a,	/* IN/OUT - array of edges to be sorted. */
int			n	/* IN - number of elements in array. */
)
{
int		h;
struct edge	tmp;
dist_t		key;
struct edge *	p1;
struct edge *	p2;
struct edge *	p3;
struct edge *	p4;
struct edge *	endp;

	endp = &a [n];

	for (h = 1; h <= n; h = 3*h+1) {
	}

	do {
		h = h / 3;
		p4 = &a [h];
		p1 = p4;
		while (p1 < endp) {
			tmp = *p1;
			key = tmp.len;
			p2 = p1;
			while (TRUE) {
				p3 = (p2 - h);
				if (p3 -> len <= key) break;
				*p2 = *p3;
				p2 = p3;
				if (p2 < p4) break;
			}
			*p2 = tmp;
			++p1;
		}
	} while (h > 1);
}
