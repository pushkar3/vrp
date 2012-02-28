/***********************************************************************

	File:	bsd.h
	Rev:	a-2
	Date:	02/28/2001

	Copyright (c) 1998, 2001 by David M. Warme

************************************************************************

	Declarations for the Bottleneck Steiner Distance
	implementation.

************************************************************************

	Modification Log:

	a-1:	10/04/98	warme
		: Created.
	a-2:	02/28/2001	warme
		: Changes for 3.1 release.  Hide certain interfaces
		:  in the .c file where they belong.

************************************************************************/

#ifndef	BSD_H
#define	BSD_H

#include "steiner.h"


/*
 * A structure to encapsulate all the data we store for rapidly
 * computing the BSD of two terminals.
 */

struct bsd {
	int			method;	/* Implementation method to use */

	/* Stuff common to all implementations. */
	int			n;		/* Number of terminals */
	struct edge *		mst_edges;	/* List of MST edges */

	/* Stuff for the Lower Triangular matrix of edge #'s implementation. */
	int16u *	ematrix;	/* The full matrix */

	/* Stuff for the row-cache implementation. */
};


extern dist_t		bsd (struct bsd *, int, int);
extern struct bsd *	compute_bsd (int, struct edge *, int);
extern void		shutdown_bsd (struct bsd *);

#endif
