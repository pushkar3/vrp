/***********************************************************************

	File:	sec_comp.h
	Rev:	b-1
	Date:	02/28/2001

	Copyright (c) 1996, 2001 by David M. Warme

************************************************************************

	Data structures for decomposition to make separation of
	SEC's easier.

************************************************************************

	Modification Log:

	a-1:	10/05/96	warme
		: Created.
	b-1:	02/28/2001	warme
		: Changes for 3.1 release.
		: Renamed many members of "struct comp" -- we now
		:  prefer "vertex" over "terminal" and "edge" over
		:  "full set" or FST.
		: Replace "tmasks" and "kmasks" with new "rverts".

************************************************************************/

#ifndef SEC_COMP_H
#define	SEC_COMP_H

#include "steiner.h"


/*
 * The following structure contains all of the information that we
 * need to represent a single component of an SEC separation problem.
 * We use a variety of techniques to split one separation problem into
 * several components, and to simplify those components.
 */

struct comp {
	struct comp *	next;		/* Linked list of components */
	int		num_verts;	/* Number of vertices in component */
	int		num_edges;	/* Number of edges in component */
	int8u		flags;		/* See below */
	/* Info about edges in the component... */
	double *	x;		/* LP assigned weight */
	int **		everts;		/* Vertices in each edge */
	/* Info about vertices in the component... */
	int **		vedges;		/* Edges incident to each vertex */
	double *	tviol;		/* Extra violation for each vertex */
	int **		rverts;		/* List of original vertices that */
					/* are represented by each vertex */
					/* (we simplify by merging multiple */
					/* vertices into one...) */
	bitmap_t *	vert_mask;	/* Mask of valid vertices */
	bitmap_t *	edge_mask;	/* Mask of valid edges */
	struct constraint * cp;		/* List of constraints for this	*/
					/* component */
};

/*
 * Various flags.  These are used to keep track of which hurdles a given
 * component has cleared.  By keeping track of this info we apply only
 * those transformations that have not previously been done to the component.
 * When a component is simplified or further divided, however, we can reset
 * these flags on the sub-component(s) -- permitting further decomposing or
 * simplification without fear of infinite looping or recursion.
 */

#define	CFLG_CONG	0x01	/* Known to contain only congested vertices */
#define CFLG_CC		0x02	/* Known to be a connected component */
#define CFLG_BCC	0x04	/* Known to be a bi-connected component */
#define	CFLG_CHAIN	0x08	/* Long integral chains merged */

#define CFLG_ALL	(CFLG_CONG | CFLG_CC | CFLG_BCC | CFLG_CHAIN)


extern struct constraint * check_component_subtour (bitmap_t *,
						    struct comp *,
						    struct constraint *,
						    double *,
						    bitmap_t *,
						    struct bbinfo *);
extern struct comp *	delete_vertex_from_component (int, struct comp *);
extern struct comp *	find_congested_components (double *	x,
						   bitmap_t *	tmap,
						   bitmap_t *	fset_mask,
						   bool		print_flag,
						   struct cinfo * cip);

extern int		find_least_congested_vertex (bitmap_t *,
						     struct comp *);
extern void		free_congested_component (struct comp *);

#endif
