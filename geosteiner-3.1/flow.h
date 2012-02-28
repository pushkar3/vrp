/***********************************************************************

	File:	flow.h
	Rev:	a-1
	Date:	07/05/96

	Copyright (c) 1996, 2001 by David M. Warme

************************************************************************

	Data structures and declarations for the maximum
	flow network solver.

************************************************************************

	Modification Log:

	a-1:	07/05/96	warme
		: Created.

************************************************************************/

#ifndef FLOW_H
#define	FLOW_H

#include "steiner.h"


/*
 * The value used to represent "infinite" flow.
 */

#define	INFINITE_FLOW		10000.0


/*
 * The following structure contains all of the information needed
 * to describe a single network flow problem instance.
 */

struct flow_prob {
	int		num_nodes;	/* Number of nodes */
	int		num_arcs;	/* Number of arcs */
	int		source;		/* Source node number */
	int		sink;		/* Sink node number */
	int **		out;		/* outgoing arcs for each node */
	int **		in;		/* incoming arcs for each node */
	int *		arc_src;	/* source node for each arc */
	int *		arc_dst;	/* destination node for each arc */
	double *	capacity;	/* capacity for each arc */
};


/*
 * The next structure contains all of the information that describes
 * a network max-flow solution.
 */

struct flow_soln {
	double		z;		/* total network src->sink flow */
	double *	flow;		/* flow for each arc */
	bitmap_t *	cut;		/* nodes on source side of cut */
};


/*
 * Finally, this structure contains temporary data structures used
 * internally by the network flow solver.  Its contents are unimportant
 * to the outside world.
 */

struct flow_temp {
	/* Pre-allocated temporary buffers used during a	*/
	/* single run of the flow solver -- may be reused for	*/
	/* several consecutive runs...				*/
	int *		pred_arc;	/* predecessor arc in path */
	double *	delta;		/* flow increment available at node */
	int *		queue;		/* breadth-first queue of nodes */
};


extern void		compute_max_flow (struct flow_prob *,
					  struct flow_temp *,
					  struct flow_soln *);
extern void		create_flow_solution_data (struct flow_prob *,
						   struct flow_soln *);
extern void		create_flow_temp_data (struct flow_prob *,
					       struct flow_temp *);
extern void		free_flow_solution_data (struct flow_soln *);
extern void		free_flow_temp_data (struct flow_temp *);


#endif
