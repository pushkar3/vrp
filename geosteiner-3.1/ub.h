/***********************************************************************

	File:	ub.h
	Rev:	b-1
	Date:	02/28/2001

	Copyright (c) 1997, 2001 by David M. Warme

************************************************************************

	Declarations pertaining to the heuristic upper bound code.

************************************************************************

	Modification Log:

	a-1:	09/06/97	warme
		: Created.
	b-1:	02/28/2001	warme
		: Add new "struct ubinfo" to better encapsulate
		:  global state info.
		: Add startup and shutdown routines.

************************************************************************/


#ifndef UB_H
#define	UB_H

#include "bb.h"
#include "steiner.h"


/*
 * The following structure contains information needed by the upper
 * bounding heuristics -- information that is computed only once, and
 * then used each time the heuristic is called.
 */

struct ubinfo {
	int	num_rankings;	/* Number of valid rankings of the FSTs */
	int *	rankings [2];	/* Various rankings of the FSTs */
	int *	mst_edges;	/* The MST edges, shortest to longest */
	double	best_z;		/* Best solution seen during heuristic */
};


extern void		compute_heuristic_upper_bound (double *,
						       struct bbinfo *);
extern void		shutdown_heuristic_upper_bound (struct ubinfo *);
extern struct ubinfo *	startup_heuristic_upper_bound (struct cinfo *);

#endif
