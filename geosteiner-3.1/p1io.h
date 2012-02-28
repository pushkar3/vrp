/***********************************************************************

	File:	p1io.h
	Rev:	b-1
	Date:	02/28/2001

	Copyright (c) 1993, 2001 by David M. Warme

************************************************************************

	Equates for reading and writing the data that is
	output from phase 1.

************************************************************************

	Modification Log:

	a-1:	10/20/96	warme
		: Created.  Split off from p1io.c.
	b-1:	02/28/2001	warme
		: Made two more routines global for 3.1 release.

************************************************************************/

#ifndef P1IO_H
#define	P1IO_H

#include "steiner.h"

/*
 * The following equates define the versions of the phase 1 data
 * format.
 */

/*
 * Version 0: (no version number)
 *
 * OR-library Steiner problem in graph format, extended to handle the
 * Steiner tree in HYPERGRAPH problem by requiring each edge (hyper or
 * otherwise) to reside with its cost on a single line.  The last number
 * on each line is the edge cost, the others are vertex numbers.  It is
 * permissible for EOF to appear instead of the "number of vertices to
 * be connected together" field.  In this case, all vertices are assumed
 * to be terminals -- resulting in an MST in graph/hypergraph problem
 * instance.
 */

#define	P1IO_VERSION_0		0

/*
 * Version 1:	Obsolete!
 */

#define	P1IO_VERSION_1		1

/*
 * Version 2:	The original multi-metric FST format:
 */

#define	P1IO_VERSION_2		2

/*
 * Version 3:
 *
 *	Added Steiner tree in hypergraph functionality, plus
 *	several other items.
 */

#define	P1IO_VERSION_3		3

#define CURRENT_P1IO_VERSION	P1IO_VERSION_3


/*
 * Function Prototypes
 */

extern void	free_phase_1_data (struct cinfo *);
extern void	print_phase_1_data (struct cinfo *, int);
extern void	read_phase_1_data (struct cinfo *);

extern int **	compute_basic_incompat (struct cinfo *);
extern void	init_term_trees (struct cinfo *);

#endif
