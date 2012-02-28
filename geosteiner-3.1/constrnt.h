/***********************************************************************

	File:	constrnt.h
	Rev:	b-1
	Date:	02/28/2001

	Copyright (c) 1996, 2001 by David M. Warme

************************************************************************

	Data structures for describing constraints.

************************************************************************

	Modification Log:

	a-1:	07/17/96	warme
		: Created.
	b-1:	02/28/2001	warme
		: Changes for 3.1 release.
		: Made add_constraint_to_pool and expand_constraint
		:  be global.

************************************************************************/

#ifndef CONSTRNT_H
#define	CONSTRNT_H

#include "steiner.h"


/*
 * The following structures represent "raw" (or physical) constraints.
 * These constraints are expressed using the actual row of coefficients,
 * operator and right-hand-side.
 *
 * Negative values of var represent the constraint's OPERATOR and rhs,
 * as well as marking the end.
 */

#define RC_OP_LE	0
#define RC_OP_EQ	1
#define RC_OP_GE	2
#define RC_VAR_BASE	3

struct rcoef {
	int16u		var;	/* variable (var>=0) or operator (var<0) */
	short		val;	/* coefficient value */
};

struct rcon {
	int		len;	/* length of constraint LHS */
	struct rcoef *	coefs;	/* the actual coefficients of the row */
	int		next;	/* next rcon number in hash bucket chain */
	int		lprow;	/* not in LP if <0, current row if >=0 */
	int		biter;	/* most recent iteration during which this */
				/* constraint was binding */
	short		hval;	/* hash value for this entry */
	short		flags;	/* various flags for entry */
	int		uid;	/* unique ID */
	int		refc;	/* reference count: number of *suspended* */
				/* nodes for which this constraint is */
				/* binding */
};

/* flags */

#define	RCON_FLAG_DISCARD	0x0001	/* Discard at next opportunity. */


/*
 * Structures used to store "logical" constraints, which are expressed
 * in terms of their significance to the problem, not their coefficients.
 */

enum ctype {
	CT_CUTSET,
	CT_SUBTOUR
};

struct constraint {
	struct constraint *	next;
	int			iteration;
	enum ctype		type;
	bitmap_t *		mask;
};


/*
 * The constraint pool.  We maintain all constraints that have ever
 * been generated in the pool, along with the hash table, freelists
 * and scratch buffers.
 */

#define	CPOOL_HASH_SIZE	1009

struct cpool {
	int			uid;	/* Bumped when pool changes */
	struct rcon *		rows;	/* All constraints, in sequence */
	int			nrows;	/* Number of constraints in the pool */
	int			maxrows; /* Allocated size of pool */
	int			num_nz;	/* Number of non-zeros in the pool */
	int *			lprows;	/* maps LP row # to constraint # */
	int			nlprows; /* number of LP rows */
	int			npend;	/* # rows pending addition to LP */
	struct rblk *		blocks;	/* List of rcon allocation blocks */
	struct rcoef *		cbuf;	/* Scratch constraint buffer */
	int			iter;	/* LP iteration count - a timestamp */
					/* used to find constraints that */
					/* haven't been useful in a while */
	int			initrows; /* number of initial rows */
	int			nvars;	/* Number of variables - LP columns */
	int			hwmrow;	/* High water mark for LP rows */
	int			hwmnz;	/* High water mark for LP non-zeros */
	int			hash [CPOOL_HASH_SIZE];
};


/*
 * This structure maintains a single block of free rcoef's, from which
 * we allocate sequentially.  To maintain better cache behavior when
 * scanning all constraints for violations, we always begin a new block
 * whenever the constraint we are adding does not fit in the current
 * block.
 */

struct rblk {
	struct rblk *		next;	/* next rblk in pool */
	struct rcoef *		base;	/* base of array of rcons */
	struct rcoef *		ptr;	/* allocation pointer */
	int			nfree;	/* number of free rcons remaining */
};


/*
 * Function Prototypes
 */

struct bbinfo;

extern bool		add_constraint_to_pool (struct cpool *,
						struct rcoef *,
						bool);
extern int		add_constraints (struct bbinfo *, struct constraint *);
extern void		add_pending_rows_to_LP (struct bbinfo *);
extern LP_t *		build_initial_formulation (struct cpool *,
						   bitmap_t *,
						   bitmap_t *,
						   struct cinfo *,
						   struct lpmem *);
extern void		debug_print_constraint (char *,
						char *,
						struct constraint *,
						double *,
						bitmap_t *,
						struct cinfo *);
extern void		delete_slack_rows_from_LP (struct bbinfo *);
extern void		destroy_initial_formulation (struct bbinfo *);
extern void		destroy_node_basis (struct bbnode *, struct bbinfo *);
extern struct rcoef *	expand_constraint (struct constraint *,
					   struct rcoef *,
					   bitmap_t *,
					   struct cinfo *);
extern void		free_constraint_pool (struct cpool *);
extern void		initialize_constraint_pool (struct cpool *,
						    bitmap_t *,
						    bitmap_t *,
						    struct cinfo *);
extern bool		is_violation (struct rcoef *, double *);
extern void		mark_row_pending_to_LP (struct cpool *, int);
extern void		restore_node_basis (struct bbnode *, struct bbinfo *);
extern void		save_node_basis (struct bbnode *, struct bbinfo *);
extern int		solve_LP_over_constraint_pool (struct bbinfo *);


#endif
