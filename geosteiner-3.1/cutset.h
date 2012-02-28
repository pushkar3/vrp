/***********************************************************************

	File:	cutset.h
	Rev:	b-2
	Date:	02/28/2001

	Copyright (c) 1996, 2001 by David M. Warme

************************************************************************

	Data structures for separating cutset constraints.

************************************************************************

	Modification Log:

	b-1:	11/14/96	warme
		: Created.
	b-2:	02/28/2001	warme
		: Changes for 3.1 release.

************************************************************************/

#ifndef CUTSET_H
#define	CUTSET_H

#include "constrnt.h"
#include "flow.h"
#include "steiner.h"


/*
 * The following data structure defines the flow graph that we use
 * to separate cutset constraints that are fractionally violated.
 */

struct cs_info {

	/* Data used by the flow solver... */
	struct flow_prob	prob;	/* The network flow formulation */
	struct flow_soln	soln;	/* The network flow solution */
	struct flow_temp	temp;	/* Temporary data structures */

	/* Data used to set the arc capacities and modify the	*/
	/* flow network during cutset separation. */
	int *		arc_to_fset;	/* arc # -> full set #. */
};


extern struct constraint * add_cutset_to_list (bitmap_t *,
					       struct constraint *,
					       double *,
					       bitmap_t *,
					       bitmap_t *,
					       struct cinfo *);
extern void		build_cutset_separation_formulation (bitmap_t *,
							     bitmap_t *,
							     struct cinfo *,
							     struct cs_info *);
extern struct constraint * find_cutset_constraints (double *,
						    bitmap_t *,
						    bitmap_t *,
						    struct cinfo *);
extern struct constraint * find_fractional_cutsets (double *,
						    struct cs_info *,
						    bitmap_t *,
						    bitmap_t *,
						    struct cinfo *);

#endif
