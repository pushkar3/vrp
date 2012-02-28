/***********************************************************************

	File:	sec2.h
	Rev:	a-1
	Date:	05/16/97

	Copyright (c) 1997, 2001 by David M. Warme

************************************************************************

	Declarations for the deterministic separation procedure
	for the "generalized SEC's" that uses a reduction to min-cut
	in a bipartite network.

************************************************************************

	Modification Log:

	a-1:	05/16/97	warme
		: Created.

************************************************************************/

#ifndef SEC2_H
#define	SEC2_H

#include "bb.h"
#include "constrnt.h"
#include "sec_comp.h"
#include "steiner.h"


/*
 * Function Prototypes
 */

extern struct constraint *	sec_flow_separator (struct comp **,
						    double *,
						    bitmap_t *,
						    struct bbinfo *,
						    struct constraint *);


#endif
