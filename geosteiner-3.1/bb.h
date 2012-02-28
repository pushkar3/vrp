/***********************************************************************

	File:	bb.h
	Rev:	a-2
	Date:	02/28/2001

	Copyright (c) 1996, 2001 by David M. Warme

************************************************************************

	Data structure for branch-and-bound (branch-and-cut)
	procedure.

************************************************************************

	Modification Log:

	a-1:	07/17/96	warme
		: Created.
	a-2:	02/28/2001	warme
		: Changes for 3.1 release, CPLEX 7.0, etc.

************************************************************************/

#ifndef BB_H
#define	BB_H

#include "config.h"
#include "steiner.h"

#ifdef CPLEX
 #if CPLEX < 40
  #define PROTOTYPE_MAX
  #include "cpxdefs.inc"
 #else
  #define CPX_PROTOTYPE_ANSI
  #include "cplex.h"
 #endif
#endif

#ifdef LPSOLVE
#include "lpkit.h"
#endif


/*
 * Constants
 */

			/* floating point "fuzz" for comparisons. */
#define	FUZZ		0.000001

	/* Constants for the heaps... */

#define	INIT_HSIZE	128
#define	NUM_BB_HEAPS	2

#define	BEST_NODE_HEAP	0
#define	WORST_NODE_HEAP	1


/*
 * The type used to represent an LP problem instance depends upon
 * whieh LP-solver we are using...
 */

#ifdef CPLEX
typedef struct cpxlp		LP_t;
#endif

#ifdef LPSOLVE
typedef lprec			LP_t;
#endif


/*
 * LP result status codes that are independent of the particular
 * solver used.
 */

#define	BBLP_OPTIMAL		0
#define	BBLP_CUTOFF		1
#define	BBLP_INFEASIBLE		2

/*
 * The following structure is used to represent a single node of the
 * branch-and-bound tree.
 */

struct bbnode {
	double		z;	/* node's current objective value -- this */
				/* is initially either the parent's or an */
				/* estimate from try_branch. */
	bool		optimal; /* optimal LP relaxation found? */
	int		num;	/* node number -- assigned in order of */
				/* NODE CREATION! */
	int		iter;	/* number of LP iterations for this node */
	int		parent;	/* node number of parent node */
	int		index [NUM_BB_HEAPS]; /* index in corresponding */
				/* bbtree heap */
	int		var;	/* var that was branched to make this node */
	int		dir;	/* direction branch var was restricted */
	int		depth;	/* depth in tree -- 0 EQ root node */
	int		br1cnt;	/* count of var=1 branch decisions */
	double *	x;	/* most recent LP solution */
	int		cpiter;	/* constraint pool time stamp, used to see */
				/* if x value is up-to-date w.r.t pool */
	double *	zlb;	/* lower bounds on Xi=0 and Xi=1 branches */
	bitmap_t *	fixed;	/* variables fixed to some value */
	bitmap_t *	value;	/* value variables are fixed at */
	int		n_uids;	/* Number of UIDs in bc_uids list */
	int *		bc_uids; /* Unique IDs of all constraints that are */
				/* binding for this node.  Only non-null */
				/* when node is suspended. */
	int *		bc_row;	/* position of bc_uid row in LP tableaux */
	int *		rstat;	/* basis info for corresponding bc_uids row */
	int *		cstat;	/* basis info for each column */
	double *	bheur;	/* Branch heuristic values */
	struct bbnode *	next;	/* next unprocessed node in LIFO order */
	struct bbnode * prev;	/* previous unprocessed node in LIFO order */
};

/*
 * The following typedef defines the type of function that is used by
 * a heap to determine if the first node belongs higher up (closer to
 * the root) than the second node.  Different heaps use different
 * comparison functions to achieve different sorting orders.
 */

typedef int	bbheap_func_t (struct bbnode *, struct bbnode *);


/*
 * The following structure is used to represent a single heap of
 * branch-and-bound nodes.
 */

struct bbheap {
	struct bbnode ** array;	/* the heap array */
	int		 hsize;	/* current heap array allocation */
	int		 nheap;	/* number of nodes in the heap array */
				/* comparision function to determine if */
				/* first node belongs above second in */
				/* this heap... */
	bbheap_func_t *	 is_parent_funcp;
};

/*
 * The following structure is used to represent the branch-and-bound
 * tree.  It contains the following ways of accessing the nodes:
 *
 *	- a linked-list for processing the nodes in depth-first order
 *	- a heap to access the current best node (lowest objective)
 *	- a heap to access the WORST node -- used to find those nodes
 *	  that must be cut-off whenever a better feasible integer
 *	  solution is found.
 */

struct bbtree {
	struct bbnode *	first;	/* first node in LIFO (depth first) order */
	struct bbnode *	free;	/* node freelist */
	int		snum;	/* node creation serial number counter */
	int		nmasks;	/* Size of "fixed" and "value" bit masks */
	int		node_policy;	/* Next node policy */
	struct bbheap	heap [NUM_BB_HEAPS]; /* heaps used to access nodes */
				/* in various orders */
};

/*
 * Next node policies...
 */

#define	NN_DEPTH_FIRST		0
#define	NN_BEST_NODE		1

/*
 * The following structure contains all of the global information that
 * we pass around between the various components of the branch-and-cut
 * procedure:
 *
 *	- the phase-1 data (original points, full sets, and
 *	    compatibility info
 *	- the main LP problem instance
 *	- the branch-and-bound tree
 *	- info used by the various separation procedures
 *	- the best feasible integer solution seen so far
 *	- the current branch-and-bound node, including:
 *		- its LP objective value and solution
 *		- its local fixed variables and values
 */

struct bbinfo {
	struct cinfo *	cip;	/* original point set, full sets, and */
				/* compatibility info */
	bitmap_t *	tmap;	/* Set of valid terminals in problem */
	bitmap_t *	fset_mask; /* Set of valid full sets */
	LP_t *		lp;	/* the main LP problem instance */
	struct lpmem *	lpmem;	/* memory allocated for LP problem instance */
	struct cpool *	cpool;	/* the global pool of constraints */
	struct bbtree *	bbtree;	/* the branch-and-bound tree */
	struct cs_info * csip;	/* cutset separation info */
	double		preempt_z; /* objective value above which to preempt */
				   /* the current node */
	double		best_z;	/* best feasible integer objective so far */
	bitmap_t *	smt;	/* best feasible integer solution so far */
	struct bbnode *	node;	/* current branch-and-bound node */
	double		_z;	/*   its objective value (obsolete) */
	double *	_x;	/*   its LP solution vector (obsolete) */
	int		slack_size; /* size of slack vector */
	double *	slack;	/*   its LP slack variables */
	double *	dj;	/*   its LP reduced costs */
	bitmap_t *	fixed;	/*   set of fixed variables */
	bitmap_t *	value;	/*   values of fixed variables */
	struct bbstats * statp;	/* statistics structure */
	cpu_time_t	t0;	/* CPU time at start of branch and cut */
	double		prevlb;	/* previous best lower bound */
	struct ubinfo *	ubip;	/* info used by upper bound heuristic */
};

/*
 * This structure is used to hold statistics for the branch-and-cut.
 */

struct cstats {			/* Constraint statistics... */
	int	num_prows;	/* Number of rows in constraint pool */
	int	num_lprows;	/* Number of rows in LP */
	int	num_pnz;	/* Number of non-zeros in constraint pool */
	int	num_lpnz;	/* Number of non-zeros in LP */
};

struct bbstats {
	int		n;		/* Number of terminals */
	int		m;		/* Number of full sets */
	cpu_time_t	p1time;		/* Phase 1 CPU time */
	cpu_time_t	p2time;		/* Phase 2 CPU time */
	double		z;		/* Integer optimal objective value */
	int		num_nodes;	/* Number of b&b nodes */
	int		num_lps;	/* Number of LP's solved */
	struct cstats	cs_init;	/* Constraint stats initially */
	struct cstats	cs_root;	/* Constraint stats for root node */
	struct cstats	cs_final;	/* Final constraint stats */
	/* Root node statistics */
	double		root_z;		/* Objective value for root node */
	bool		root_opt;	/* Is root_z optimal? */
	int		root_lps;	/* Number of LP's solved at root */
	cpu_time_t	root_time;	/* CPU time to finish root node */
};

/*
 * Here are a bunch of macros that we use to insulate us from the
 * calling convention differences between CPLEX versions 3.0 and 4.0.
 */

#ifdef CPLEX
 #if CPLEX >= 40
  /* Version 4.0 or greater...					*/

  /* For simplicity, make the CPLEX environment be a global.	*/
  extern CPXENVptr	cplex_env;

  /* These changed in 5.0...					*/
  #if CPLEX >= 50
   #define _MYCPX_copybase(lp, cstat, rstat) \
		(CPXcopybase (cplex_env, lp, cstat, rstat))
   #define _MYCPX_primopt(lp)	(CPXprimopt (cplex_env, lp))
   #define _MYCPX_INFBOUND	CPX_INFBOUND
  #else
   #define _MYCPX_copybase(lp, cstat, rstat) \
		(CPXloadbase (cplex_env, lp, cstat, rstat))
   #define _MYCPX_primopt(lp)	(CPXoptimize (cplex_env, lp))
   #define _MYCPX_INFBOUND	INFBOUND
  #endif

  /* Unchanged since 4.0...					*/

  #define _MYCPX_addrows(lp, ccnt, rcnt, nzcnt, rhs, sense, \
			 rmatbeg, rmatind, rmatval, colname, rowname) \
		(CPXaddrows (cplex_env, lp, ccnt, rcnt, nzcnt, rhs, sense, \
			     rmatbeg, rmatind, rmatval, colname, rowname))
  #define _MYCPX_chgbds(lp, cnt, index, lu, bd) \
		(CPXchgbds (cplex_env, lp, cnt, index, lu, bd))
  #define _MYCPX_delsetrows(lp, delstat) \
		(CPXdelsetrows (cplex_env, lp, delstat))
  #define _MYCPX_dualopt(lp)	(CPXdualopt (cplex_env, lp))
  #define _MYCPX_freeprob(lpp)	(CPXfreeprob (cplex_env, lpp))
  #define _MYCPX_getbase(lp, cstat, rstat) \
		(CPXgetbase (cplex_env, lp, cstat, rstat))
  #define _MYCPX_getdj(lp, dj, begin, end) \
		(CPXgetdj (cplex_env, lp, dj, begin, end))
  #define _MYCPX_getnumcols(lp)	(CPXgetnumcols (cplex_env, lp))
  #define _MYCPX_getnumnz(lp)	(CPXgetnumnz (cplex_env, lp))
  #define _MYCPX_getnumrows(lp)	(CPXgetnumrows (cplex_env, lp))
  #define _MYCPX_getnzspace(lp)	(CPXgetnzspace (cplex_env, lp))
  #define _MYCPX_getobj(lp,obj,b,e) (CPXgetobj (cplex_env, lp, obj, b, e))
  #define _MYCPX_getrowspace(lp) (CPXgetrowspace (cplex_env, lp))
  #define _MYCPX_getslack(lp, slack, begin, end) \
		(CPXgetslack (cplex_env, lp, slack, begin, end))
  #define _MYCPX_loadlp(probname, numcols, numrows, objsen, obj, rhs, \
			 sense, matbeg, matcnt, matind, matval, lb, ub, \
			 rngval, colspace, rowspace, nzspace) \
		(CPXloadlp (cplex_env, probname, numcols, numrows, objsen, \
			    obj, rhs, sense, matbeg, matcnt, matind, matval, \
			    lb, ub, rngval, colspace, rowspace, nzspace))
  #define _MYCPX_lpwrite(lp, fname) (CPXlpwrite (cplex_env, lp, fname))
  #define _MYCPX_openCPLEX(stp) (CPXopenCPLEX (stp))
  #define _MYCPX_setadvind(value) \
		(CPXsetintparam (cplex_env, CPX_PARAM_ADVIND, value))
  #define _MYCPX_setlogfile(stream) (CPXsetlogfile (cplex_env, stream))
  #define _MYCPX_setobjulim(limit, small, big) \
		(CPXsetdblparam (cplex_env, CPX_PARAM_OBJULIM, limit))
  #define _MYCPX_setscaind(flag, small, big) \
		(CPXsetintparam (cplex_env, CPX_PARAM_SCAIND, flag))
  #define _MYCPX_solution(lp, stat, z, x, pi, slack, dj) \
		(CPXsolution (cplex_env, lp, stat, z, x, pi, slack, dj))

  #define _MYCPX_MIN	CPX_MIN
  #define _MYCPX_MAX	CPX_MAX
 #else
  /* version 3.0 */

  #define _MYCPX_addrows(lp, ccnt, rcnt, nzcnt, rhs, sense, \
			 rmatbeg, rmatind, rmatval, colname, rowname) \
		(addrows (lp, ccnt, rcnt, nzcnt, rhs, sense, \
			  rmatbeg, rmatind, rmatval, colname, rowname))
  #define _MYCPX_chgbds(lp, cnt, index, lu, bd) \
		(chgbds (lp, cnt, index, lu, bd))
  #define _MYCPX_delsetrows(lp, delstat) (delsetrows (lp, delstat))
  #define _MYCPX_dualopt(lp)	(dualopt (lp))
  #define _MYCPX_freeprob(lpp)	(freeprob (lpp), 0)
  #define _MYCPX_getbase(lp, cstat, rstat) \
		(getbase (lp, cstat, rstat))
  #define _MYCPX_getdj(lp, dj, begin, end) \
		(getdj (lp, dj, begin, end))
  #define _MYCPX_getnumcols(lp)	(getmac (lp))
  #define _MYCPX_getnumnz(lp)	(getmat (lp))
  #define _MYCPX_getnumrows(lp) (getmar (lp))
  #define _MYCPX_getnzspace(lp) (getmatsz (lp))
  #define _MYCPX_getobj(lp,obj,b,e) (getobj (lp, obj, b, e))
  #define _MYCPX_getrowspace(lp) (getmarsz (lp))
  #define _MYCPX_getslack(lp, slack, begin, end) \
		(getslack (lp, slack, begin, end))
  #define _MYCPX_copybase(lp, cstat, rstat) \
		(loadbase (lp, cstat, rstat))
  #define _MYCPX_loadlp(probname, numcols, numrows, objsen, obj, rhs, \
			 sense, matbeg, matcnt, matind, matval, lb, ub, \
			 rngval, colspace, rowspace, nzspace) \
		(loadprob (probname, numcols, numrows, 0, objsen, \
			    obj, rhs, sense, matbeg, matcnt, matind, matval, \
			    lb, ub, rngval, \
			    NULL, NULL, NULL, NULL, NULL, \
			    NULL, NULL, NULL, NULL, NULL, \
			    NULL, NULL, NULL, NULL, NULL, \
			    NULL, NULL, \
			    colspace, rowspace, nzspace, 0, 0, 0, 0, 0))
  #define _MYCPX_lpwrite(lp, fname) (lpwrite (lp, fname))
  #define _MYCPX_primopt(lp)	(optimize (lp))
  #define _MYCPX_setadvind(value) (setadvind (value))
  #define _MYCPX_setlogfile(stream) (setlogfile (stream))
  #define _MYCPX_setobjulim(limit, small, big) \
		(setobjulim (limit, small, big))
  #define _MYCPX_setscaind(flag, small, big) \
		(setscaind (flag, small, big))
  #define _MYCPX_solution(lp, stat, z, x, pi, slack, dj) \
		(solution (lp, stat, z, x, pi, slack, dj))

  #define _MYCPX_MIN	1
  #define _MYCPX_MAX	-1
  #define _MYCPX_INFBOUND	INFBOUND
 #endif
#endif

/*
 * Some macros to do common things to LP's
 */

#ifdef CPLEX
#define GET_LP_NUM_COLS(lp)	(_MYCPX_getnumcols (lp))
#define GET_LP_NUM_ROWS(lp)	(_MYCPX_getnumrows (lp))
#define	GET_LP_NUM_NZ(lp)	(_MYCPX_getnumnz (lp))
#endif

#ifdef LPSOLVE
#define	GET_LP_NUM_COLS(lp)	((lp) -> columns)
#define	GET_LP_NUM_ROWS(lp)	((lp) -> rows)
#define	GET_LP_NUM_NZ(lp)	((lp) -> non_zeros)
#endif


/*
 * A structure to keep track of dynamic memory used by an LP.
 */

#ifdef CPLEX
struct lpmem {
	double *	objx;
	double *	rhsx;
	char *		senx;
	int *		matbeg;
	int *		matcnt;
	int *		matind;
	double *	matval;
	double *	bdl;
	double *	bdu;
	int		obj_scale;	/* objective scale factor */
};
#endif

#ifdef LPSOLVE
struct lpmem {
	/* lp_solve_2.0 dynamically manages the LP tableaux memory... */
	int	dummy;
};
#endif


/*
 * Extern declarations
 */

extern bool		Seed_Pool_With_2SECs;


extern dist_t		branch_and_cut (struct bbinfo *);
extern bool		check_for_better_IFS (double *,
					      struct bbinfo *,
					      double *);
extern struct bbinfo *	create_bbinfo (struct cinfo *);
extern void		new_upper_bound (double, struct bbinfo *);


/*
 * This flag is set by a signal handler.  We use it to manually terminate
 * constraint generation for the current node, thereby forcing a branch.
 */

extern volatile bool		force_branch_flag;


#endif
