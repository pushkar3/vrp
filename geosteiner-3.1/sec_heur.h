/***********************************************************************

	File:	sec_heur.h
	Rev:	b-1
	Date:	02/28/2001

	Copyright (c) 1996, 2001 by David M. Warme

************************************************************************

	Data structures for the various heuristic separation
	procedures for Strong SEC's (Subtour Elimination
	Constraints).

************************************************************************

	Modification Log:

	a-1:	05/13/96	warme
		: Created.
	b-1:	02/28/2001	warme
		: Changes for 3.1 release.
		: Made check_unique_subtour be global.
		: Change calling convention of enumerate_all_subtours
		:  and find_small_subtours.

************************************************************************/

#ifndef SEC_HEUR_H
#define	SEC_HEUR_H

#include "bb.h"
#include "constrnt.h"
#include "sec_comp.h"
#include "steiner.h"
/*
 * Size of largest subproblem that we separate via exhaustive enumeration.
 */

#define	SEC_ENUM_LIMIT		10


/*
 * Function Prototypes
 */

	/* Pre-declarations of structure tags... */
	struct bbnode;
	struct constraint;

extern struct constraint *	check_subtour (bitmap_t *,
					       struct constraint *,
					       double *,
					       bitmap_t *,
					       struct cinfo *);
extern struct constraint *	check_unique_subtour (bitmap_t *,
						      int,
						      struct constraint *);
extern struct constraint *	enumerate_all_subtours (struct comp *,
							struct constraint *,
							struct bbinfo *);
extern struct constraint *	find_integer_cycles (double *,
						     bitmap_t *,
						     bitmap_t *,
						     struct constraint *,
						     struct cinfo *);
extern struct constraint *	find_small_subtours (struct comp *,
						     struct constraint *,
						     struct bbinfo *);
extern bool			is_equal (bitmap_t *, bitmap_t *, int);
extern struct constraint *	sec_flow_heuristic (struct comp *,
						    double *,
						    bitmap_t *,
						    struct cinfo *,
						    struct constraint *);


#endif
