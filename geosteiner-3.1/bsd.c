/***********************************************************************

	File:	bsd.c
	Rev:	a-3
	Date:	02/28/2001

	Copyright (c) 1997, 2001 by David M. Warme

************************************************************************

	Bottleneck Steiner Distance stuff.

************************************************************************

	Modification Log:

	a-1:	09/04/97	warme
		: Split off from coarse.c.
	a-2:	10/04/98	warme
		: Completely re-coded to conserve memory and provide
		:  multiple implementations that provide different
		:  space-time tradeoffs.
	a-3:	02/28/2001	warme
		: Changes for 3.1 release.  Migrate data from bsd
		:  structure to local vars.

************************************************************************/

#include "bsd.h"
#include "steiner.h"


/*
 * Global Routines
 */

dist_t			bsd (struct bsd *, int, int);
struct bsd *		compute_bsd (int, struct edge *, int);
void			shutdown_bsd (struct bsd *);


/*
 * Local Types
 */

	/* A structure for representing an MST in adjacency list format. */

struct mstadj {
	int		edge;		/* Edge number */
	int		vnum;		/* Index of other vertex */
};


/*
 * Local Routines
 */

static struct mstadj **	build_adjacency_list (int, int, int, struct edge *);
static void		walk (int,
			      int,
			      struct mstadj **,
			      bool *,
			      int *);

/*
 * This routine returns the Bottleneck Steiner Distance between
 * two terminals, specified by index.
 */

	dist_t
bsd (

struct bsd *	bsdp,		/* IN - BSD data structure */
int		i,		/* IN - first terminal */
int		j		/* IN - second terminal */
)
{
int		n;

	n = bsdp -> n;
	if ((i < 0) OR (i >= n) OR (j < 0) OR (j >= n)) {
		fatal ("bsd: Bug 1.");
	}

	if (i EQ j) return (0.0);

	switch (bsdp -> method) {
	case 1:	/* Complete lower-triangular matrix. */
		if (i > j) {
			i = ((i * (i - 1)) >> 1) + j;
		}
		else {
			i = ((j * (j - 1)) >> 1) + i;
		}
		j = bsdp -> ematrix [i];
		return (bsdp -> mst_edges [j].len);

	case 2:	/* Cache of BSD rows. */
		/* Not yet implemented! */
		fatal ("bsd: Bug 2.");
		break;

	default:
		fatal ("bsd: Bug 1.");
		break;
	}

	return (0.0);
}

/*
 * This routine initializes the BSD data structures according to the
 * given point set and implementation method.  The first thing is to
 * compute the actual minimum spanning tree.  The rest of the initialization
 * is method-specific.
 */

	struct bsd *
compute_bsd (

int			nedges,		/* IN - number of MST edges */
struct edge *		mst_edges,	/* IN - edges of the MST */
int			method		/* IN - implementation method */
)
{
int			i;
int			j;
int			nterms;
struct bsd *		bsdp;
int16u *		rowp;
int *			tmprow;
struct mstadj **	adj_list;
struct edge *		edges;
bool *			mark;

	nterms = nedges + 1;

	/* Allocate and zero the BSD structure. */
	bsdp = NEW (struct bsd);
	memset (bsdp, 0, sizeof (*bsdp));

	if (method EQ 0) {
		/* Default is to use the full matrix for 4000 points or	*/
		/* less.  Use the cache of rows for larger problems.	*/
#if 0
		method = (nterms <= 4000) ? 1 : 2;
#else
		method = 1;
#endif
	}
	bsdp -> method	= method;

	bsdp -> n	= nterms;
	edges		= NEWA (nedges + 1, struct edge); /* 1 extra! */

	/* Add initial zero-length edge, so that we have an edge number	*/
	/* that represents a BSD of zero.				*/
	edges [0].len	= 0.0;
	edges [0].p1	= 0;
	edges [0].p2	= 0;

	memcpy (&edges [1], mst_edges, nedges * sizeof (mst_edges [0]));

	bsdp -> mst_edges = edges;

	adj_list = build_adjacency_list (nterms, 1, nedges + 1, edges);

	mark = NEWA (nterms, bool);
	for (i = 0; i < nterms; i++) {
		mark [i] = FALSE;
	}

	switch (method) {
	case 1:	/* Complete lower triangular matrix. */

		bsdp -> ematrix = NEWA (nterms * (nterms - 1) >> 1, int16u);
		tmprow = NEWA (nterms, int);

		/* Fill in the matrix, one row at a time... */
		rowp = bsdp -> ematrix;
		for (i = 0; i < nterms; i++) {
			walk (i, 0, adj_list, mark, tmprow);
			for (j = 0; j < i; j++) {
				*rowp++ = tmprow [j];
			}
		}
		free ((char *) tmprow);
		break;

	case 2:
		fatal ("compute_bsd: Bug 2.");
		break;

	default:
		fatal ("compute_bsd: Bug 3.");
		break;
	}

	free ((char *) mark);
	free ((char *) adj_list [0]);
	free ((char *) adj_list);

	return (bsdp);
}

/*
 * This routine converts a list of edges into a full graph structure
 * represented in adjacency list form.  The adjacency list is in two parts:
 * an array of pointers indexed by node number, and an array of "adj"
 * structures indexed by these pointers.  To find all neighbors of node K,
 * look at every "mstadj" structure between adj_list [K] and adj_list [K+1].
 */

	static
	struct mstadj **
build_adjacency_list (

int			n,		/* IN - number of nodes */
int			first_edge,	/* IN - first edge number */
int			nedges,		/* IN - number of edges */
struct edge *		edges		/* IN - array of edges */
)
{
int			i;
struct edge *		ep;
struct mstadj *		ap;
int *			count;
struct mstadj **	tmp_ptr;
struct mstadj **	adj_list;

	adj_list	= NEWA (n + 1, struct mstadj *);
	ap		= NEWA (2 * nedges, struct mstadj);
	tmp_ptr		= NEWA (n, struct mstadj *);
	count		= NEWA (n, int);

	for (i = 0; i < n; i++) {
		count [i] = 0;
	}

	/* Count neighbors for each node. */
	ep = &edges [first_edge];
	for (i = first_edge; i < nedges; ep++, i++) {
		++(count [ep -> p1]);
		++(count [ep -> p2]);
	}

	/* Set up the pointers for each node. */
	for (i = 0; i < n; i++) {
		adj_list [i]	= ap;
		tmp_ptr [i]	= ap;
		ap += count [i];
	}
	adj_list [i] = ap;		/* set overall end of list. */

	/* Deposit neighbors into properly allocated lists. */
	ep = &edges [first_edge];
	for (i = first_edge; i < nedges; ep++, i++) {
		ap = tmp_ptr [ep -> p1]++;
		ap -> edge	= i;
		ap -> vnum	= ep -> p2;

		ap = tmp_ptr [ep -> p2]++;
		ap -> edge	= i;
		ap -> vnum	= ep -> p1;
	}

	free ((char *) count);
	free ((char *) tmp_ptr);

	return (adj_list);
}

/*
 * This routine recursively walks through the given Minimum Spanning Tree
 * starting at the given node, and determines the longest edge along the
 * path from the original given root node, to every other node in the tree.
 *
 * Note: this routine assumes that edges are listed in increasing order,
 * just as Kruskal's MST algorithm emits them.
 */

	static
	void
walk (

int			node,		/* IN - current node in tree */
int			longest,	/* IN - longest edge # in cur path */
struct mstadj **	adj_list,	/* IN - MST adjacency list */
bool *			mark,		/* IN - vertices visited */
int *			rowp		/* OUT - array of longest edge #'s */
)
{
struct mstadj *		ap;
struct mstadj *		endp;

	/* We now know the longest edge along the path from the root */
	/* to this current node!  Remember it in the matrix! */
	rowp [node] = longest;

	mark [node] = TRUE;
	ap	= adj_list [node];
	endp	= adj_list [node + 1];
	while (ap < endp) {
		if (NOT mark [ap -> vnum]) {
			walk (ap -> vnum,
			     (longest >= ap -> edge) ? longest : ap -> edge,
			      adj_list,
			      mark,
			      rowp);
		}
		++ap;
	}
	mark [node] = FALSE;
}

/*
 * This routine destroys and deallocates the BSD information.
 */

	void
shutdown_bsd (

struct bsd *	bsdp		/* IN - BSD info to destroy */
)
{
	if (bsdp EQ NULL) {
		return;
	}

	if (bsdp -> ematrix NE NULL) {
		free ((char *) (bsdp -> ematrix));
		bsdp -> ematrix = NULL;
	}

	if (bsdp -> mst_edges NE NULL) {
		free ((char *) (bsdp -> mst_edges));
		bsdp -> mst_edges = NULL;
	}

	free ((char *) bsdp);
}
