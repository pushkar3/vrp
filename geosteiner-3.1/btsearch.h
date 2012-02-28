/***********************************************************************

	File:	search.h
	Rev:	a-1
	Date:	11/30/2000

	Copyright (c) 1993, 2001 by David M. Warme

************************************************************************

	Declarations for the special backtrack search used
	by the FST pruning code.

************************************************************************

	Modification Log:

	a-1:	11/30/2000	warme
		: Created.

************************************************************************/

#ifndef BTSEARCH_H
#define	BTSEARCH_H

#include "steiner.h"


/*
 * Function Prototypes
 */

extern double		backtrack_search (struct cinfo *, bitmap_t *);
extern void		initialize_btsearch (void);


#endif
