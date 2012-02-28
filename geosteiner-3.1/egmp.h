/***********************************************************************

	File:	egmp.h
	Rev:	b-1
	Date:	11/25/2000

	Copyright (c) 2000, 2001 by David M. Warme

************************************************************************

	Support routines for the EFST generator that use the
	GNU Multi-Precision arithmetic library (GMP -- if we have
	it) to compute certain items with high accuracy.

************************************************************************

	Modification Log:

	b-1:	11/25/2000	warme
		: Created.

************************************************************************/


#ifdef HAVE_GMP

#ifndef EGMP_H
#define EGMP_H

#include "gmp.h"


struct einfo;
struct eqp_t;


typedef struct {			/* a + b * sqrt(3) */
	mpq_t		a;
	mpq_t		b;
} qr3_t;

struct qr3_point {
	qr3_t		x;
	qr3_t		y;
};


/*
 * Global Routines
 */

extern double	compute_EFST_length (struct einfo * eip, struct eqp_t * eqpt);
extern void	qr3_clear (qr3_t * p);
extern void	qr3_init (qr3_t * p);
extern void	update_eqpoint_and_displacement (struct einfo *	eip,
						 struct eqp_t *	eqpk);

#endif

#endif
