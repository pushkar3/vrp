/***********************************************************************

	File:	emptyr.h
	Rev:	a-1
	Date:	09/28/98

	Copyright (c) 1998, 2001 by David M. Warme and Martin Zachariasen

************************************************************************

	Routines for efficiently determining whether or not two
	terminals define an empty rectangle.  We precompute this
	information and store it compactly.

************************************************************************

	Modification Log:

	a-1:	09/28/98	warme
		: Created.  Implemented Zachariasen's algorithm
		:  using Warme's infrastructure.

************************************************************************/

#ifndef	EMPTYR_H
#define	EMPTYR_H


#include "steiner.h"


/*
 * Global Routines
 */

extern int		count_empty_rectangles (bitmap_t *, int);
extern bitmap_t *	init_empty_rectangles (struct pset *, int *);
extern bool		is_empty_rectangle (bitmap_t *, int, int);
extern void		shutdown_empty_rectangles (bitmap_t *);

#endif
