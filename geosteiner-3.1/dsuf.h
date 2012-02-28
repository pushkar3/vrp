/***********************************************************************

	File:	dsuf.h
	Rev:	a-1
	Date:	11/03/98

	Copyright (c) 1998, 2001 by David M. Warme

************************************************************************

	Declarations for the "disjoint set union-find" data structure.

************************************************************************

	Modification Log:

	a-1:	11/03/98	warme
		: Created.

************************************************************************/

#ifndef	DSUF_H
#define	DSUF_H


#include "steiner.h"


/*
 * This is the so-called "Disjoint Set Union-Find" data structure.
 * The operations are "makeset", "find", and "union".
 *
 * See chapter 2 of "Data Structures and Network Algorithms", by
 * Robert Endre Tarjan, SIAM, 1983 for complete details.
 */


struct dsuf {
	int *	parent;
	int *	rank;
	int	set_size;
};


/*
 * Global Routines
 */

extern void		dsuf_create (struct dsuf *, int);
extern void		dsuf_destroy (struct dsuf *);
extern int		dsuf_find (struct dsuf *, int);
extern void		dsuf_makeset (struct dsuf *, int);
extern void		dsuf_unite (struct dsuf *, int, int);

#endif
