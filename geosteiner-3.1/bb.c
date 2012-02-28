/***********************************************************************

	File:	bb.c
	Rev:	b-2
	Date:	02/28/2001

	Copyright (c) 1995, 2001 by David M. Warme

************************************************************************

	The guts of the branch-and-cut.

************************************************************************

	Modification Log:

	a-1:	11/17/95	warme
		: Created.
	b-1:	11/14/96	warme
		: Renamed this program.
		: Split off the cutset stuff into another file.
		: Other cleanups for release.
	b-2:	02/28/2001	warme
		: Numerous changes for 3.1 release.
		: Split the main() routine off into bbmain.c.
		: Split certain utility routines off into bbsubs.c.
		: Major improvements to branch var selection.
		: New create_bbinfo routine does setup for
		:  branch_and_cut.
		: Enable variable fixing code.

************************************************************************/

#include "bb.h"
#include "bbsubs.h"
#include "config.h"
#include "constrnt.h"
#include "cra.h"
#include "cutset.h"
#include "genps.h"
#include "p1io.h"
#include "sec2.h"
#include "sec_comp.h"
#include "sec_heur.h"
#include "steiner.h"
#include "ub.h"

#include <signal.h>


/*
 * Global Routines
 */

dist_t			branch_and_cut (struct bbinfo *);
bool			check_for_better_IFS (double *,
					      struct bbinfo *,
					      double *);
struct bbinfo *		create_bbinfo (struct cinfo *);
void			new_upper_bound (double, struct bbinfo *);


int		Branch_Var_Policy = 1;
int		Check_Branch_Vars_Thoroughly = 1;
bool		Check_Root_Constraints = FALSE;
bool		Choose_Branch_Vars_Carefully = TRUE;
double		Initial_Upper_Bound = DBL_MAX;
bool		Print_Root_LP = FALSE;

volatile bool	force_branch_flag;


/*
 * Local Equates
 */


	/* Lower bound outcomes... */
#define	LB_INFEASIBLE		1
#define	LB_CUTOFF		2
#define	LB_INTEGRAL		3
#define	LB_FRACTIONAL		4
#define LB_PREEMPTED		5

	/* Variable fixing outcomes... */
#define	VFIX_NOTHING_FIXED	0
#define	VFIX_VARIABLES_FIXED	1
#define	VFIX_FIXED_FRACTIONAL	2
#define	VFIX_INFEASIBLE		3

#define	UP_FIRST	TRUE


/*
 * Local Types
 */

struct bvar {			/* Info about best branch var seen */
	int	var;		/* Best branch var, or -1 */
	double	z0;		/* Objective when Xj=0 */
	double	z1;		/* Objective when Xj=1 */
	double	test_2nd_val;	/* Only check 2nd branch if 1st > this. */
};

#ifdef CPLEX

struct basis_save {	/* Structure to save basis state for CPLEX */
	int *		cstat;
	int *		rstat;
};

#endif


/*
 * Local Routines
 */

static int		carefully_choose_branching_variable (struct bbinfo *,
							     double *,
							     double *);
static void		change_var_bounds (LP_t *, int, double, double);
static void		check_root_constraints (struct bbinfo *);
static int		choose_branching_variable (struct bbinfo *,
						   double *,
						   double *);
static bool		compare_branch_vars (struct bbinfo *,
					     int,
					     struct bvar *);
static int		compute_good_lower_bound (struct bbinfo *);
static void		cut_off_existing_nodes (double, struct bbtree *);
static struct constraint * do_separations (struct bbinfo *,
					   cpu_time_t **);
static bool		eval_branch_var (struct bbinfo *,
					 int,
					 int,
					 struct basis_save *,
					 double);
static int		fix_variables (struct bbinfo *,
				       int *, int,
				       int *, int);
static RETSIGTYPE	force_branch_signal_handler (int);
static bool		integer_feasible_solution (double *,
						   bitmap_t *,
						   bitmap_t *,
						   struct cinfo *,
						   int *);
static void		new_lower_bound (double, struct bbinfo *);
static void		print_root_lp (struct bbinfo *);
static int		reduced_cost_var_fixing (struct bbinfo *);
static struct bbnode *	select_next_node (struct bbtree *);
static void		sort_branching_vars (int *, int, double *);
static void		trace_node (struct bbinfo *, char, char *);
static void		update_node_preempt_value (struct bbinfo *);

#ifdef CPLEX
static void		destroy_LP_basis (struct basis_save *);
static double		try_branch (LP_t *,
				    int,
				    int,
				    double *,
				    double,
				    struct basis_save *);
static void		save_LP_basis (LP_t *, struct basis_save *);
#endif


/*
 * Local Variables
 */

static FILE *		rcfile;

/*
 * Set up the initial branch-and-bound problem, including the root
 * node and the initial constraint pool.
 */

	struct bbinfo *
create_bbinfo (

struct cinfo *		cip		/* IN - compatibility info */
)
{
struct bbinfo *		bbip;
int			i;
int			j;
int			k;
int			n;
int			nmasks;
int			nedges;
LP_t *			lp;
int			status;
bitmap_t *		vert_mask;
bitmap_t *		edge_mask;
bitmap_t *		req_edges;
bitmap_t *		fixed;
bitmap_t *		value;
bitmap_t *		smt;
struct cpool *		cpool;
struct bbstats *	statp;
struct bbtree *		bbtree;
struct bbnode *		root;
struct lpmem *		lpmem;
struct rcon *		rcp;

	/* Create the global branch-and-bound info structure... */
	bbip = NEW (struct bbinfo);
	memset (bbip, 0, sizeof (*bbip));

	bbip -> t0 = get_cpu_time ();

	nmasks = cip -> num_edge_masks;
	nedges = cip -> num_edges;

	vert_mask	= cip -> initial_vert_mask;
	edge_mask	= cip -> initial_edge_mask;
	req_edges	= cip -> required_edges;

	/* initialize global pool of constraints. */
	cpool = NEW (struct cpool);
	initialize_constraint_pool (cpool, vert_mask, edge_mask, cip);

	/* Build initial formulation. */
	lpmem = NEW (struct lpmem);
	lp = build_initial_formulation (cpool,
					vert_mask,
					edge_mask,
					cip,
					lpmem);

	/* Initialize the branch-and-bound tree... */
	bbtree = create_bbtree (nmasks);

	/* Create vectors to describe the current problem... */
	fixed	= NEWA (nmasks, bitmap_t);
	value	= NEWA (nmasks, bitmap_t);
	smt	= NEWA (nmasks, bitmap_t);

	for (i = 0; i < nmasks; i++) {
		fixed [i] = 0;
		value [i] = 0;
		smt [i]	  = 0;
	}
	for (i = 0; i < nedges; i++) {
		if (NOT BITON (edge_mask, i)) {
			/* variables that are outside of the problem */
			/* are fixed at zero... */
			SETBIT (fixed, i);
			change_var_bounds (lp, i, 0.0, 0.0);
		}
		else if (BITON (req_edges, i)) {
			/* Front-end has determined that this hyperedge	*/
			/* MUST be present in an optimal solution!	*/
			SETBIT (fixed, i);
			SETBIT (value, i);
			change_var_bounds (lp, i, 1.0, 1.0);
		}
	}

	/* Create the root node... */
	root = NEW (struct bbnode);
	memset (root, 0, sizeof (*root));

	root -> z	= -DBL_MAX;
	root -> optimal	= FALSE;
	root -> num	= (bbtree -> snum)++;
	root -> iter	= 0;
	root -> parent	= -1;
	for (i = 0; i < NUM_BB_HEAPS; i++) {
		root -> index [i] = -1;
	}
	root -> var	= -1;
	root -> dir	= 0;
	root -> depth	= 0;
	root -> br1cnt	= 0;
	root -> x	= NEWA (nedges, double);
	root -> cpiter	= -1;		/* x is not current. */
	root -> zlb	= NEWA (2 * nedges, double);
	root -> fixed	= fixed;
	root -> value	= value;
	root -> n_uids	= 0;
	root -> bc_uids	= NULL;
	root -> bc_row	= NULL;
	root -> rstat	= NULL;
	root -> cstat	= NULL;
	root -> bheur	= NEWA (nedges, double);
	root -> next	= NULL;
	root -> prev	= NULL;

	for (i = 0; i < nedges; i++) {
		root -> bheur [i] = 0.0;
	}

	for (i = 0; i < 2*nedges; i++) {
		root -> zlb [i] = -DBL_MAX;
	}

	/* Create the branch-and-bound statistics structure... */
	statp = NEW (struct bbstats);
	memset (statp, 0, sizeof (*statp));

	statp -> n		= cip -> num_verts;
	statp -> m		= nedges;
	statp -> p1time		= cip -> p1time;
	statp -> num_nodes	= 0;
	statp -> num_lps	= 0;

	statp -> cs_init.num_prows	= cpool -> nrows;
	statp -> cs_init.num_lprows	= GET_LP_NUM_ROWS (lp);
	statp -> cs_init.num_pnz	= cpool -> num_nz;
	statp -> cs_init.num_lpnz	= GET_LP_NUM_NZ (lp);

	/* Fill in the global branch-and-bound info structure... */
	bbip -> cip		= cip;
	bbip -> tmap		= vert_mask;
	bbip -> fset_mask	= edge_mask;
	bbip -> lp		= lp;
	bbip -> lpmem		= lpmem;
	bbip -> cpool		= cpool;
	bbip -> bbtree		= bbtree;
	bbip -> csip		= NULL;
	bbip -> preempt_z	= Initial_Upper_Bound;
	bbip -> best_z		= Initial_Upper_Bound;
	bbip -> smt		= smt;
	bbip -> node		= root;
	bbip -> slack_size	= 0;
	bbip -> slack		= NULL;
	bbip -> dj		= NEWA (nedges, double);
	bbip -> fixed		= NULL;
	bbip -> value		= NULL;
	bbip -> statp		= statp;
	bbip -> prevlb		= -DBL_MAX;
	bbip -> ubip		= NULL;

	/* Make the root node inactive by putting it in the bbtree... */

	n = cpool -> nlprows;
	root -> n_uids	= n;
	root -> bc_uids	= NEWA (n, int);
	root -> bc_row	= NEWA (n, int);

	j = 0;
	rcp = &(cpool -> rows [0]);
	for (i = 0; i < cpool -> nrows; i++, rcp++) {
		k = rcp -> lprow;
		if (k < 0) continue;
		++(rcp -> refc);
		root -> bc_uids [j] = rcp -> uid;
		root -> bc_row [j]  = k;
		++j;
	}
	if (j NE n) {
		fatal ("create_bbinfo: Bug 1.");
	}

	root -> next = NULL;
	root -> prev = NULL;
	bbtree -> first = root;

	bbheap_insert (root, bbtree, BEST_NODE_HEAP);
	bbheap_insert (root, bbtree, WORST_NODE_HEAP);

	return (bbip);
}

/*
 * This routine is the top-level of the branch-and-cut.
 */

	dist_t
branch_and_cut (

struct bbinfo *		bbip		/* IN - branch-and-bound info */
)
{
int			i;
int			j;
int			nmasks;
int			nedges;
int			status;
bitmap_t *		fixed;
bitmap_t *		value;
bitmap_t *		delta;
bitmap_t *		smt;
double *		x;
double			z0;
double			z1;
double			tmpz;
cpu_time_t		t1;
RETSIGTYPE		(*save_sigterm) (int);
struct cinfo *		cip;
LP_t *			lp;
struct cpool *		cpool;
struct bbtree *		bbtree;
struct bbstats *	statp;
struct bbnode *		node;
struct bbnode *		node_to_free;
struct bbnode *		node2;

#ifdef CPLEX
int *			b_index;
char *			b_lu;
double *		b_bd;
double			objlim;
double			save_objlim;
#endif

	cip	= bbip -> cip;
	cpool	= bbip -> cpool;
	lp	= bbip -> lp;
	statp	= bbip -> statp;
	smt	= bbip -> smt;
	bbtree	= bbip -> bbtree;

	nmasks = cip -> num_edge_masks;
	nedges = cip -> num_edges;

#if CPLEX
	/* Save the existing objective limit, and set it to infinity. */
	CPXgetdblparam (cplex_env, CPX_PARAM_OBJULIM, &save_objlim);
	objlim = DBL_MAX;
	CPXsetdblparam (cplex_env, CPX_PARAM_OBJULIM, objlim);
#endif

	if (Check_Root_Constraints) {
		rcfile = fopen ("/tmp/lp.x", "w");
		if (rcfile EQ NULL) {
			fprintf (stderr,
				 "Warning: unable to create /tmp/lp.x\n");
			Check_Root_Constraints = FALSE;
		}
	}

	/* Establish handler for the SIGTERM (kill -15) signal... */
	force_branch_flag = FALSE;
	save_sigterm = signal (SIGTERM, force_branch_signal_handler);

	/* Create vectors to describe the current problem... */
	fixed	= NEWA (nmasks, bitmap_t);
	value	= NEWA (nmasks, bitmap_t);
	delta	= NEWA (nmasks, bitmap_t);

	/* No variables have been fixed yet... */
	for (i = 0; i < nmasks; i++) {
		fixed [i] = 0;
		value [i] = 0;
	}
	bbip -> fixed = fixed;
	bbip -> value = value;

#ifdef CPLEX
	/* Create arrays for changing variable bounds... */
	b_index	= NEWA (2 * nedges, int);
	b_lu	= NEWA (2 * nedges, char);
	b_bd	= NEWA (2 * nedges, double);
#endif

#if 0
	/* Build cutset separation formulation. */
	bbip -> csip = NEW (struct cs_info);
	build_cutset_separation_formulation (vert_mask,
					     edge_mask,
					     cip,
					     bbip -> csip);
#endif

	/* Init the heuristic upper bound. */

	bbip -> ubip = startup_heuristic_upper_bound (bbip -> cip);


	/* At this point, all nodes are inactive. */

	for (;;) {

		/* Select the next node to process. */
		node = select_next_node (bbtree);
		if (node EQ NULL) break;

		/* This is perhaps a new lower bound... */
		new_lower_bound (node -> z, bbip);

		if (node -> z > -DBL_MAX) {
			tracef ("%% Resuming node %d at %24.20f\n",
				node -> num,
				UNSCALE (node -> z, &(cip -> scale)));
		}
		else {
			tracef ("%% Resuming node %d\n", node -> num);
		}

		/* Restore the LP tableaux and basis for this node.	*/
		/* Decrement the reference counts on this node's	*/
		/* binding rows.  Since it is now the ACTIVE node,	*/
		/* there is no reason to continue protecting these rows	*/
		/* from deletion by other nodes.			*/
		restore_node_basis (node, bbip);

		/* Determine new preemption value (i.e. the objective	*/
		/* value of the next-best node).			*/
		update_node_preempt_value (bbip);

		/* Modify LP to represent problem from new node. */
		for (i = 0; i < nmasks; i++) {
			delta [i] =   (fixed [i] ^ node -> fixed [i])
				    | (value [i] ^ node -> value [i]);
		}
#ifdef CPLEX
		j = 0;
		for (i = 0; i < nedges; i++) {
			if (NOT BITON (delta, i)) continue;
			/* Force bounds for variable 'i' to be correct... */
			b_index [j]	= i;	/* variable i, */
			b_lu [j]	= 'L';	/*	lower bound */
			b_index [j+1]	= i;	/* variable i, */
			b_lu [j+1]	= 'U';	/*	upper bound */
			if (NOT BITON (node -> fixed, i)) {
				/* new variable is NOT fixed... */
				b_bd [j]	= 0.0;
				b_bd [j+1]	= 1.0;
			}
			else if (NOT BITON (node -> value, i)) {
				/* new variable is fixed to 0 */
				b_bd [j]	= 0.0;
				b_bd [j+1]	= 0.0;
			}
			else {
				/* new variable is fixed to 1 */
				b_bd [j]	= 1.0;
				b_bd [j+1]	= 1.0;
			}
			j += 2;
		}
		if (j > 0) {
			if (_MYCPX_chgbds (bbip -> lp, j, b_index, b_lu, b_bd) NE 0) {
				fatal ("branch_and_cut: Bug 1.");
			}
#if 0
			++(bbip -> cpool -> uid);
#endif
		}
#endif

#ifdef LPSOLVE
		j = 0;
		for (i = 0; i < nedges; i++) {
			if (NOT BITON (delta, i)) continue;
			++j;
			/* Force bounds on variable 'i' to be correct... */
			if (NOT BITON (node -> fixed, i)) {
				/* variable is NOT fixed... */
				set_bounds (lp, i + 1, 0.0, 1.0);
			}
			else if (NOT BITON (node -> value, i)) {
				/* variable is fixed to 0 */
				set_bounds (lp, i + 1, 0.0, 0.0);
			}
			else {
				/* variable is fixed to 1 -- must set	*/
				/* bounds in this order to avoid	*/
				/* lb > ub condition between calls...	*/
				set_bounds (lp, i + 1, 1.0, 1.0);
			}
		}
		if (j > 0) {
#if 0
			++(bbip -> cpool -> uid);
#endif
		}
#endif

		for (i = 0; i < nmasks; i++) {
			fixed [i] = node -> fixed [i];
			value [i] = node -> value [i];
		}

		/* Set up new node to be processed */
		bbip -> node = node;
		if (node -> iter <= 0) {
			/* Haven't processed this node before	*/
			/* -- tally another node...		*/
			++(statp -> num_nodes);
		}

		/* Process the current node... */
		status = compute_good_lower_bound (bbip);

		if (node -> depth EQ 0) {
			/* Finished the root node... */
			statp -> root_z = node -> z;
			statp -> root_lps = statp -> num_lps;

			/* Slack rows should have already been deleted... */
			statp -> cs_root.num_prows = cpool -> nrows;
			statp -> cs_root.num_lprows =
				GET_LP_NUM_ROWS (bbip -> lp);
			statp -> cs_root.num_pnz = cpool -> num_nz;
			statp -> cs_root.num_lpnz =
				GET_LP_NUM_NZ (bbip -> lp);
			statp -> root_opt = node -> optimal;

			t1 = get_cpu_time ();
			statp -> root_time = t1 - bbip -> t0;

			if ((status EQ LB_FRACTIONAL) AND Print_Root_LP) {
				print_root_lp (bbip);
			}

			if (Check_Root_Constraints) {
				check_root_constraints (bbip);
			}
		}

		node_to_free = NULL;	/* Default is no node to free... */

		switch (status) {
		case LB_INFEASIBLE:
			/* Node is fathomed! */
			trace_node (bbip, ' ', "infeasible");
			node_to_free = node;
			break;

		case LB_CUTOFF:
			/* Node is fathomed! */
			trace_node (bbip, ' ', "cutoff");
			node_to_free = node;
			break;

		case LB_INTEGRAL:
			if (node -> z >= bbip -> best_z) {
				trace_node (bbip, ' ', "cutoff");
				break;
			}
			/* We have a new best feasible integer solution! */
			for (i = 0; i < nmasks; i++) {
				smt [i] = 0;
			}
			x = node -> x;
			for (i = 0; i < nedges; i++) {
				if (x [i] >= 0.5) {
					SETBIT (smt, i);
				}
			}
			new_upper_bound (node -> z, bbip);

			trace_node (bbip, '*', NULL);

			node_to_free = node;
			break;

		case LB_FRACTIONAL:
			j = choose_branching_variable (bbip, &z0, &z1);
			if (j < 0) {
				/* At least one variable was fixed due	*/
				/* to cutoff or infeasibility.  It is	*/
				/* possible that the entire node is now	*/
				/* cutoff...				*/
				if (node -> z >= bbip -> best_z) {
					trace_node (bbip, ' ', "cutoff");
					node_to_free = node;
					break;
				}
				goto suspend;
			}

			/* Create two nodes... */
			if (z0 < z1) {
				add_bbnode (bbip, j, 0, z0);
				add_bbnode (bbip, j, 1, z1);
			}
			else if (z1 < z0) {
				add_bbnode (bbip, j, 1, z1);
				add_bbnode (bbip, j, 0, z0);
			}
			else if (UP_FIRST) {	/* To break ties... */
				add_bbnode (bbip, j, 0, z0);
				add_bbnode (bbip, j, 1, z1);
			}
			else {
				add_bbnode (bbip, j, 1, z1);
				add_bbnode (bbip, j, 0, z0);
			}
			trace_node (bbip, ' ', NULL);

			/* This node is done (became 2 children), free it. */
			node_to_free = node;
			break;

		case LB_PREEMPTED:
suspend:
			tracef ("%% suspending node %d at %24.20f\n",
				node -> num,
				UNSCALE (node -> z, &(cip -> scale)));
			/* This node is no longer the best.  Put it	*/
			/* back into the heap and get another one...	*/
			node2 = bbtree -> first;
			if (node2 NE NULL) {
				node2 -> prev = node;
			}
			node -> next = node2;
			node -> prev = NULL;
			bbtree -> first = node;

			bbheap_insert (node, bbtree, BEST_NODE_HEAP);
			bbheap_insert (node, bbtree, WORST_NODE_HEAP);

			/* Deactivating this node -- remember the basis */
			save_node_basis (node, bbip);

			/* Do NOT free this node! */
			break;
		}

		/* If there is a node to free, do so now... */
		if (node_to_free NE NULL) {
			/* Free up saved basis info and decrement	*/
			/* constraint reference counts before freeing.	*/
			destroy_node_basis (node_to_free, bbip);

			node_to_free -> next = bbtree -> free;
			bbtree -> free = node_to_free;
		}
	}

	statp -> cs_final.num_prows	= cpool -> nrows;
	statp -> cs_final.num_lprows	= GET_LP_NUM_ROWS (bbip -> lp);
	statp -> cs_final.num_pnz	= cpool -> num_nz;
	statp -> cs_final.num_lpnz	= GET_LP_NUM_NZ (bbip -> lp);

	/* New lower bound! */
	new_lower_bound (bbip -> best_z, bbip);

	t1 = get_cpu_time ();
	statp -> p2time = t1 - bbip -> t0;
	statp -> z = bbip -> best_z;

	/* Restore normal handling of SIGTERM... */
	(void) signal (SIGTERM, save_sigterm);

#if CPLEX
	free ((char *) b_bd);
	free ((char *) b_lu);
	free ((char *) b_index);
#endif

	free ((char *) delta);
	free ((char *) value);
	free ((char *) fixed);

#if CPLEX
	CPXsetdblparam (cplex_env, CPX_PARAM_OBJULIM, save_objlim);
#endif

	/* Return length of final tree... */

	return ((dist_t) bbip -> best_z);
}

/*
 * A handler for the SIGTERM signal.  When we receive this signal we
 * stop generating constraints for the current node and branch instead.
 */

	static
	RETSIGTYPE
force_branch_signal_handler (

int		sig		/* IN - signal being handled */
)
{
	/* Stop generating constraints and force a branch instead. */
	force_branch_flag = TRUE;

	/* Prepare to handle the signal again... */
	(void) signal (SIGTERM, force_branch_signal_handler);
}

/*
 * This routine selects the next node to process from the given
 * branch-and-bound tree.  This is where we implement the specific
 * search policy.
 */

	static
	struct bbnode *
select_next_node (

struct bbtree *		tp	/* IN - branch-and-bound tree */
)
{
struct bbnode *		p;

	if (tp -> first EQ NULL) {
		/* No more nodes! */
		return (NULL);
	}

	switch (tp -> node_policy) {
	case NN_DEPTH_FIRST:
		/* Get node created most recently... */
		p = tp -> first;
		break;

	case NN_BEST_NODE:
		/* Get node with lowest objective function value... */
		p = tp -> heap [BEST_NODE_HEAP].array [0];
		break;

	default:
		fatal ("select_next_node: Bug 1.");
	}

	delete_node_from_bbtree (p, tp);

	return (p);
}

/*
 * This routine traces the result for a given node.
 */

	static
	void
trace_node (

struct bbinfo *		bbip,	/* IN - the branch-and-bound info */
char			c1,	/* IN - space or * char */
char *			msg	/* IN - message OR NULL */
)
{
struct bbnode *		p;
struct bbtree *		tp;
struct bbheap *		hp;
double			lowz;
char *			p1;
char *			p2;
char			c2;
char			buf1 [32];
char			buf2 [32];
char			buf3 [32];
char			buf4 [32];
char			buf5 [32];
char			buf6 [32];
char			line [136];

	p	= bbip -> node;
	tp	= bbip -> bbtree;

	if (msg EQ NULL) {
		(void) sprintf (buf1, "%14.4f", p -> z);
	}
	else {
		(void) sprintf (buf1, "%14s", msg);
	}
	if (bbip -> best_z EQ DBL_MAX) {
		buf2 [0] = '\0';
	}
	else {
		(void) sprintf (buf2, "%14.4f", bbip -> best_z);
	}
	hp = &(tp -> heap [BEST_NODE_HEAP]);
	if (hp -> nheap > 0) {
		(void) sprintf (buf3,
				"%14.4f",
				hp -> array [0] -> z);
	}
	else {
		buf3 [0] = '\0';
	}
	if (p -> var < 0) {
		buf4 [0] = '\0';
		c2 = ' ';
	}
	else {
		(void) sprintf (buf4, "x%d", p -> var);
		c2 = (p -> dir EQ 0) ? 'D' : 'U';
	}
	if (p -> parent < 0) {
		buf5 [0] = '\0';
	}
	else {
		(void) sprintf (buf5, "%d", p -> parent);
	}
	if (p -> depth <= 0) {
		buf6 [0] = '\0';
	}
	else {
		(void) sprintf (buf6, "%d", p -> depth);
	}

	(void) sprintf (line,
			"%c%6d%6d%14s%14s%14s%6s %c%6s%6s",
			c1,		/* space or * */
			p -> num,	/* node number */
			hp -> nheap,	/* nodes left */
			buf1,		/* objective/cutoff/infeas */
			buf2,		/* best integer soln */
			buf3,		/* best node */
			buf4,		/* variable name */
			c2,		/* branch direction */
			buf5,		/* parent node number */
			buf6);		/* node depth */

	p1 = line;
	for (p2 = p1; *p2 NE '\0'; p2++) {
	}
	while ((p2 > p1) AND (p2 [-1] EQ ' ')) {
		--p2;
	}
	*p2 = '\0';
	(void) tracef ("  %% %s\n", line);
}

/*
 * This routine plots the LP relaxation solution for the root node.  We
 * do this whenever -r is specified and we get an optimal LP relaxation
 * solution that is fractional.
 */

	static
	void
print_root_lp (

struct bbinfo *		bbip		/* IN - branch-and-bound info */
)
{
struct cinfo *		cip;
struct bbnode *		nodep;
struct bbstats *	statp;
char *			descr;
char			buf1 [32];
char			buf2 [32];
char			title [256];

	cip	= bbip -> cip;
	nodep	= bbip -> node;
	statp	= bbip -> statp;

	convert_cpu_time (statp -> root_time, buf1);
	dist_to_string (buf2, nodep -> z, &(cip -> scale));

	descr = "Steiner Minimal Tree";
	if ((cip -> description NE NULL) AND
	    (cip -> description [0] NE '\0')) {
		descr = cip -> description;
	}
	sprintf (title,
		 "%s:  %lu points,  Root Node, LP %d,  length = %s,  %s seconds",
		 descr,
		 cip -> num_verts,
		 nodep -> iter,
		 buf2,
		 buf1);

	plot_lp_solution (title, bbip -> node -> x, cip, BIG_PLOT);
}

/*
 * This routine chooses the next variable to branch on at this
 * node.
 */

	static
	int
choose_branching_variable (

struct bbinfo *		bbip,		/* IN - branch-and-bound info */
double *		z0,		/* OUT - value to give Xi=0 node */
double *		z1		/* OUT - value to give Xi=1 node */
)
{
int			i;
int			nedges;
int			best_var;
struct cinfo *		cip;
struct bbnode *		nodep;
double *		x;
bitmap_t *		fset_mask;
double			xi;
double			max_infeas;
double			infeas;

	if (Choose_Branch_Vars_Carefully) {
		/* Do it very carefully! */
		return (carefully_choose_branching_variable (bbip, z0, z1));
	}

	/* There are LOTS of things that we MIGHT do here in	*/
	/* the FUTURE!  For right now, just take the variable	*/
	/* that is closest to 1/2 -- take larger variables in	*/
	/* the event of a tie...				*/

	cip		= bbip -> cip;
	nodep		= bbip -> node;
	x		= nodep -> x;
	fset_mask	= bbip -> fset_mask;

	nedges	= cip -> num_edges;

	best_var	= -1;
	max_infeas	= 0.0;

	for (i = nedges - 1; i >= 0; i--) {
		if (NOT BITON (fset_mask, i)) continue;
		xi = x [i];
		if (xi <= FUZZ) continue;
		if (xi + FUZZ >= 1.0) continue;
		infeas = 1.0 - xi;
		if (xi < infeas) {
			infeas = xi;
		}
		if (infeas > max_infeas) {
			best_var = i;
			max_infeas = infeas;
		}
	}
	if (best_var < 0) {
		fatal ("choose_branching_variable: Bug 1.");
	}

	/* Give both nodes the same value... */
	*z0	= nodep -> z;
	*z1	= nodep -> z;

	return (best_var);
}

/*
 * This routine does a very careful job of choosing the next variable
 * to branch on.  For each fractional variable Xi, we solve the LP
 * first with Xi=0 and then Xi=1, yielding objective values Zi0 and Zi1
 * correspondingly.  We then choose variable Xi for which the value
 * min(Zi0, Zi1) is MAXIMIZED.  This provides us with the best possible
 * "one-branch/no-cuts" improvement in the lower bound.
 */

	static
	int
carefully_choose_branching_variable (

struct bbinfo *		bbip,		/* IN - branch-and-bound info */
double *		node_z0,	/* OUT - value for Xi=0 node */
double *		node_z1		/* OUT - value for Xi=1 node */
)
{
int			i;
int			j;
int			n;
int			nedges;
int			nfrac;
int			limit;
int			new_limit;
int			logn;
int			failure_limit;
int			num_failures;
struct cinfo *		cip;
struct bbnode *		nodep;
double *		x;
double *		rank;
bitmap_t *		edge_mask;
int *			fvars;
bool			fixed;
bool			cur_var_is_better;
double			xi;
double			z0;
double			z1;
double			z;
double			delta;
double			test_2nd_val;
double			num;
double			den;
struct bvar		best;
struct basis_save	bsave;

	cip		= bbip -> cip;
	nodep		= bbip -> node;
	x		= nodep -> x;
	edge_mask	= bbip -> fset_mask;

	nedges	= cip -> num_edges;

	fvars = NEWA (nedges, int);

start_all_over:

	nfrac = 0;
	for (i = 0; i < nedges; i++) {
		if (NOT BITON (edge_mask, i)) continue;
		xi = x [i];
		if (xi <= FUZZ) continue;
		if (xi + FUZZ >= 1.0) continue;
		fvars [nfrac++] = i;
	}
#if 1
	tracef ("\n %% Carefully choosing branching variable, nfrac = %d\n",
		nfrac);
#endif

	/* Compute final heuristic ranking of the candidate branch	*/
	/* variables.  We multiply the node's "bheur" values by a	*/
	/* "complexity factor" that is 1 for variables sitting at 0.5,	*/
	/* and increases as the variable's "rational" value becomes	*/
	/* more complicated.  Thus, we prefer vars stuck at 1/2, and	*/
	/* then prefer vars stuck at 1/3 or 2/3, vars stuck at 1/4 or	*/
	/* 3/4, etc.							*/

	rank = NEWA (nedges, double);
	memcpy (rank, nodep -> bheur, nedges * sizeof (double));

	for (j = 0; j < nfrac; j++) {
		i = fvars [j];

		/* Compute closest rational approximation. */

		(void) cra (x [i], &num, &den);

		/* Factor is 1 + the denominator's "distance" from 1/2	*/
		/* + the numerator's "distance" from 1/2.		*/

		rank [i] *= ((den - 1.0) + fabs (num - 0.5 * den));
	}

	/* Sort the fractional variables so that good branch choices	*/
	/* appear first with high probability.				*/

	sort_branching_vars (fvars, nfrac, rank);

#if 0
	tracef (" %% Branch Variable Info:\n");
	for (j = 0; j < nfrac; j++) {
		i = fvars [j];
		tracef (" %% %d:\tx%d\t= %.6f,\tbheur = %g,\trank = %g\n",
			j, i, x [i], nodep -> bheur [i], rank [i]);
	}
#endif

	free ((char *) rank);

	/* Snapshot the current basis so that we can quickly	*/
	/* get back to it each time...				*/
	save_LP_basis (bbip -> lp, &bsave);

	/* Compute the non-improvement limit.  When we have tested this	*/
	/* many consecutive variables without finding a better choice,	*/
	/* we punt.  We use 2 * log(N), where N is the number of	*/
	/* fractional vars, and log(x) is the floor of the base-2 log.	*/

	failure_limit = nfrac;
	if (nfrac >= 20) {
		logn = 0;
		for (n = nfrac; n > 1; n >>= 1) {
			++logn;
		}
		failure_limit = 2 * Check_Branch_Vars_Thoroughly * logn;
	}

	best.var	= -1;
	best.z0		= nodep -> z;
	best.z1		= nodep -> z;

	test_2nd_val	= -DBL_MAX;

	/* Do a quick scan without forcing anything to determine the	*/
	/* best initial choice of branch variable.			*/

	for (j = 0; j < nfrac; j++) {
		i = fvars [j];
		if (i < 0) continue;	/* var was fixed! */

		cur_var_is_better = compare_branch_vars (bbip, i, &best);
		if (cur_var_is_better) {
			test_2nd_val = best.test_2nd_val;
		}
	}

#if 1
	if (best.var >= 0) {
		tracef (" %% Initial guess is x%d,"
			" Z0 = %-24.15g, Z1 = %-24.15g\n\n",
			best.var,
			best.z0,
			best.z1);
	}
#endif

	/* Now do the expensive part -- testing one or both branches	*/
	/* of good candidate variables.					*/

	new_limit = -1;
	limit = nfrac;

	num_failures = 0;

again:

	for (j = 0; j < limit; j++) {

		if (force_branch_flag) {
			if (best.var >= 0) goto get_out;

			/* Ignore this interrupt! */
			force_branch_flag = FALSE;
		}

		i = fvars [j];
		if (i < 0) continue;	/* var was fixed! */

		xi = x [i];

		z0 = nodep -> zlb [2 * i + 0];
		z1 = nodep -> zlb [2 * i + 1];

		if (((z0 > nodep -> z) AND (z0 < z1)) OR
		    ((z0 EQ z1) AND (xi <= 0.5))) {
			/* Check the Xi=0 branch, and then the Xi=1 branch. */
			fixed = eval_branch_var (bbip,
						 i,
						 0,	/* Xi=0, then Xi=1 */
						 &bsave,
						 test_2nd_val);
		}
		else {
			/* Check the Xi=1 branch, and then the Xi=0 branch. */
			fixed = eval_branch_var (bbip,
						 i,
						 1,	/* Xi=1, then Xi=0 */
						 &bsave,
						 test_2nd_val);
		}

		if (fixed) {
#if 1
			/* Special return code that says to try */
			/* re-solving the LP again.		*/
			destroy_LP_basis (&bsave);
			free ((char *) fvars);
			return (-1);
#elif 0
			goto start_all_over;
#else
			fvars [j] = -1;
			new_limit = j;
			limit = nfrac;
			num_failures = 0;
#endif
		}

		/* Test the current var to see if it is better than the	*/
		/* best seen so var.					*/

		cur_var_is_better = compare_branch_vars (bbip, i, &best);

		if (cur_var_is_better) {
			/* Establish a new threshold for testing 2nd branch. */
			test_2nd_val = best.test_2nd_val;

			z0 = nodep -> zlb [2 * i + 0];
			z1 = nodep -> zlb [2 * i + 1];

			z = z0;
			if (z1 < z) {
				z = z1;
			}

#if 1
			tracef (" %%   New best:  x%d, Z = %-24.15g\n",
				best.var, z);
#endif

			if (z >= bbip -> best_z) {
				/* Nice deal!  This variable	*/
				/* forces a cutoff on both	*/
				/* branches!  No need to look	*/
				/* at the rest...		*/
				new_limit = -1;
				break;
			}
			num_failures = 0;
		}
		else {
			++num_failures;
			if (num_failures >= failure_limit) {
#if 1
				tracef (" %% %d consecutive failures: giving up.\n",
					num_failures);
#endif
				goto get_out;
			}
		}
	}

	if (best.var < 0) {
		fatal ("carefully_choose_branching_variable: Bug 1.");
	}

	if (new_limit >= 0) {
		/* We have fixed some variables.  Must retest. */
		limit = new_limit;
		new_limit = -1;
		goto again;
	}

get_out:

#if 1
	tracef (" %% Best branch is x%d, Z0 = %-24.15g, Z1 = %-24.15g\n\n",
		best.var, best.z0, best.z1);
#endif

	destroy_LP_basis (&bsave);

	free ((char *) fvars);

	*node_z0 = best.z0;
	*node_z1 = best.z1;

	return (best.var);
}

/*
 * Sort the given list of fractional variables so that the best branching
 * variables appear early in the list with high probability.
 */

	static
	void
sort_branching_vars (

int *			fvars,	/* IN/OUT - fractional variables to sort */
int			n,	/* IN - number of fractional vars */
double *		rank	/* IN - heuristic rank of each variable */
)
{
int			i, i1, i2, j, k;

	/* Construct the heap via sift-downs, in O(n) time. */
	for (k = n >> 1; k >= 0; k--) {
		j = k;
		for (;;) {
			i = (j << 1) + 1;
			if (i + 1 < n) {
				/* Increment i (to right subchild of j) */
				/* if the right subchild is larger. */
				i1 = fvars [i];
				i2 = fvars [i + 1];
				if (rank [i2] > rank [i1]) {
					++i;
				}
			}
			if (i >= n) {
				/* Hit bottom of heap, sift-down is done. */
				break;
			}
			i1 = fvars [j];
			i2 = fvars [i];
			if (rank [i2] < rank [i1]) {
				/* Largest child is smaller.  Sift-	*/
				/* down is done. */
				break;
			}
			/* Sift down and continue. */
			fvars [j] = i2;
			fvars [i] = i1;
			j = i;
		}
	}

	/* Now do actual sorting.  Exchange first/last and sift down. */
	while (n > 1) {
		/* Largest is at fvars [0], swap with fvars [n-1],	*/
		/* thereby putting it into final position.		*/
		--n;
		i = fvars [0];
		fvars [0] = fvars [n];
		fvars [n] = i;

		/* Now restore the heap by sifting fvars [0] down. */
		j = 0;
		for (;;) {
			i = (j << 1) + 1;
			if (i + 1 < n) {
				/* Increment i (to right subchild of j) */
				/* if the right subchild is larger. */
				i1 = fvars [i];
				i2 = fvars [i + 1];
				if (rank [i2] > rank [i1]) {
					++i;
				}
			}
			if (i >= n) {
				/* Hit bottom of heap, sift-down is done. */
				break;
			}
			i1 = fvars [j];
			i2 = fvars [i];
			if (rank [i2] < rank [i1]) {
				/* Largest child is smaller.  Sift-	*/
				/* down is done. */
				break;
			}
			/* Sift down and continue. */
			fvars [j] = i2;
			fvars [i] = i1;
			j = i;
		}
	}
}

/*
 * This routine does the guts of carefully testing a single fractional
 * branching variable.  The caller indicates which direction should
 * be tested first.
 */

	static
	bool
eval_branch_var (

struct bbinfo *		bbip,		/* IN - branch-and-bound info */
int			var,		/* IN - variable to branch */
int			dir1,		/* IN - first branch direction */
struct basis_save *	basp,		/* IN - basis to restore when done */
double			test_2nd_val	/* IN - test 2nd if 1st is > this */
)
{
int			i;
int			nedges;
int			dir2;
struct cinfo *		cip;
LP_t *			lp;
struct bbnode *		nodep;
bool			found;
bool			fixed;
double *		x;
double			z;

	cip	= bbip -> cip;
	lp	= bbip -> lp;
	nodep	= bbip -> node;

	nedges = cip -> num_edges;

	x = NEWA (nedges, double);

	dir2 = 1 - dir1;

	fixed = FALSE;

	/* Try the first branch direction... */
	z = try_branch (lp, var + 1, dir1, x, DBL_MAX, basp);

#if CPLEX
	z = ldexp (z, bbip -> lpmem -> obj_scale);
#endif

	/* Check for a better integer feasible solution... */
	found = check_for_better_IFS (x, bbip, &z);

	/* Update per-variable lower bounds. */
	i = 2 * var + dir1;
	if (z > nodep -> zlb [i]) {
		nodep -> zlb [i] = z;
	}
	else {
		z = nodep -> zlb [i];
	}

#if 1
	tracef (" %%%s\tx%d = %d,\tZ%d = %-24.15g\n",
		found ? " !!!" : "",
		var, dir1, dir1, z);
#endif

	/* Try finding a good heuristic solution on the branched solution. */
	compute_heuristic_upper_bound (x, bbip);

#if 1
	if (z >= bbip -> best_z + 1.0e-8 * fabs (bbip -> best_z)) {
		/* Cutoff or infeasible.  Var must be fixed to	*/
		/* other direction.  Reoptimize and get new	*/
		/* basis.					*/
		SETBIT (bbip -> fixed, var);
		SETBIT (bbip -> node -> fixed, var);
		if (dir1) {
			CLRBIT (bbip -> value, var);
			CLRBIT (bbip -> node -> value, var);
		}
		else {
			SETBIT (bbip -> value, var);
			SETBIT (bbip -> node -> value, var);
		}
		change_var_bounds (lp,
				   var,
				   (double) dir2,
				   (double) dir2);

		nodep -> cpiter = -1;	/* force re-solve of LP */

		solve_LP_over_constraint_pool (bbip);
		lp = bbip -> lp;
		save_LP_basis (lp, basp);

		/* Try finding a good heuristic solution on the	*/
		/* new fixed solution...			*/
		compute_heuristic_upper_bound (bbip -> node -> x, bbip);

		/* The variable has been fixed! */
		fixed = TRUE;
		goto all_done;
	}
#endif

	if (z <= test_2nd_val)
	{
		/* No need to test the second branch. */
		goto all_done;
	}

	/* Try the second branch direction... */
	z = try_branch (lp, var + 1, dir2, x, DBL_MAX, basp);

#if CPLEX
	z = ldexp (z, bbip -> lpmem -> obj_scale);
#endif

	/* Check for better integer feasible solution... */
	found = check_for_better_IFS (x, bbip, &z);

	/* Update per-variable lower bounds. */
	i = 2 * var + dir2;
	if (z > nodep -> zlb [i]) {
		nodep -> zlb [i] = z;
	}
	else {
		z = nodep -> zlb [i];
	}

#if 1
	tracef (" %%%s\tx%d = %d,\tZ%d = %-24.15g\n",
		found ? " !!!" : "",
		var, dir2, dir2, z);
#endif

	/* Try finding a good heuristic solution on the branched solution. */
	compute_heuristic_upper_bound (x, bbip);

#if 1
	if (z >= bbip -> best_z + 1.0e-8 * fabs (bbip -> best_z)) {
		/* Cutoff or infeasible.  Var must be fixed to	*/
		/* other direction.  Reoptimize and get new	*/
		/* basis.					*/
		SETBIT (bbip -> fixed, var);
		SETBIT (bbip -> node -> fixed, var);
		if (dir2) {
			CLRBIT (bbip -> value, var);
			CLRBIT (bbip -> node -> value, var);
		}
		else {
			SETBIT (bbip -> value, var);
			SETBIT (bbip -> node -> value, var);
		}
		change_var_bounds (lp,
				   var,
				   (double) dir1,
				   (double) dir1);

		nodep -> cpiter = -1;	/* force re-solve of LP */

		solve_LP_over_constraint_pool (bbip);
		lp = bbip -> lp;
		save_LP_basis (lp, basp);

		/* Try finding a good heuristic solution on the	*/
		/* new fixed solution...			*/
		compute_heuristic_upper_bound (bbip -> node -> x, bbip);

		/* The variable has been fixed! */
		fixed = TRUE;
	}
#endif

all_done:

	free ((char *) x);

	return (fixed);
}

/*
 * See if one candidate branch variable is better than another.
 * We implement various policies here.
 */

	static
	bool
compare_branch_vars (

struct bbinfo *		bbip,		/* IN - branch and bound info */
int			i1,		/* IN - first branch var */
struct bvar *		bvp		/* IN/OUT - current best branch var */
)
{
struct bbnode *		nodep;
bool			cur_var_is_better;
double			z;
double			z0;
double			z1;
double			zmin;
double			zmax;
double			best_z0;
double			best_z1;
double			best_zmin;
double			best_zmax;
double			ub;
double			gap;
double			delta;
double			prod;
double			best_prod;
double			test2;

#define TOLERANCE	1.0e-10

	nodep	= bbip -> node;

	ub	= bbip -> best_z;
	z	= nodep -> z;

	if (i1 < 0) {
		return (FALSE);
	}

	z0	= nodep -> zlb [2 * i1 + 0];
	z1	= nodep -> zlb [2 * i1 + 1];

	if (z0 < z1) {
		zmin = z0;
		zmax = z1;
	}
	else {
		zmin = z1;
		zmax = z0;
	}

	if (bvp -> var < 0) {
		/* New branch var is better.  Still have to provide a	*/
		/* test_2nd_val threshold, which is policy specific.	*/

		switch (Branch_Var_Policy) {
		case 0:
			test2 = zmin;
			break;

		case 1:
			test2 = zmin - TOLERANCE * fabs (zmin);
			break;

		case 2:
			gap = fabs (ub - z);
			if (gap <= 0.0) {
				fatal ("compare_branch_vars: Bug 1.");
			}
			prod = fabs ((z0 - z) * (z1 - z));
			test2 = z + prod / gap;
			break;

		default:
			fatal ("compare_branch_vars: Bug 2.");
			break;
		}

		bvp -> var		= i1;
		bvp -> z0		= z0;
		bvp -> z1		= z1;
		bvp -> test_2nd_val	= test2;

		return (TRUE);
	}

	best_z0	= bvp -> z0;
	best_z1	= bvp -> z1;

	if (best_z0 < best_z1) {
		best_zmin = best_z0;
		best_zmax = best_z1;
	}
	else {
		best_zmin = best_z1;
		best_zmax = best_z0;
	}

	switch (Branch_Var_Policy) {
	case 0:
		/* Naive max of mins. */
		cur_var_is_better = FALSE;
		if (zmin > best_zmin) {
			cur_var_is_better = TRUE;
			test2 = zmin;
		}
		break;

	case 1:
		/* Smarter lexicographic max of mins.  If the mins	*/
		/* are "about equal", use the max values to decide.	*/
fuzzy_lexical:
		cur_var_is_better = FALSE;

		delta = TOLERANCE * fabs (best_zmin);
		if ((zmin - best_zmin) > delta) {
			cur_var_is_better = TRUE;
		}
		else if ((fabs (zmin - best_zmin) <= delta) AND
			 (zmax > best_zmax)) {
			cur_var_is_better = TRUE;
		}
		test2 = zmin - TOLERANCE * fabs (zmin);
		break;

	case 2:
		/* Product of improvements.  Uses method 1 to break	*/
		/* close ties.						*/
		prod = fabs ((z0 - z) * (z1 - z));
		best_prod = fabs ((best_z0 - z) * (best_z1 - z));

		/* Compute tolerance factor that is a fraction of the	*/
		/* current gap.						*/
		gap = fabs (ub - z);
		if (gap <= 0.0) {
			fatal ("compare_branch_vars: Bug 3.");
		}
		delta = 1.0E-5 * gap;
		if (fabs (prod - best_prod) <= delta) {
			/* Products are nearly equal.  Use fuzzy	*/
			/* lexicographic max of mins.			*/
			goto fuzzy_lexical;
		}
		cur_var_is_better = FALSE;
		if (prod > best_prod) {
			cur_var_is_better = TRUE;
			test2 = z + prod / gap;
		}
		break;

	default:
		fatal ("compare_branch_vars: Bug 4.");
		break;
	}

	if (cur_var_is_better) {

		bvp -> var		= i1;
		bvp -> z0		= z0;
		bvp -> z1		= z1;
		bvp -> test_2nd_val	= test2;
	}

	return (cur_var_is_better);

#undef TOLERANCE
}

/*
 * This routine checks to see if the result of doig a "test-branch" on
 * a variable JUST HAPPENS to result in an integer feasible solution
 * that is the best so far.  Although this would be total serendipity,
 * we do it anyway.  The check takes virtually no time because we almost
 * always locate a fractional variable within the first few probes.
 * And besides, it would be a shame to NOT notice something this important
 * -- especially if we don't have ANY upper bound yet!
 */

	bool
check_for_better_IFS (

double *		x,		/* IN - solution to test */
struct bbinfo *		bbip,		/* IN - branch and bound info */
double *		true_z		/* OUT - true Z value, if integral */
)
{
int			i;
int			nedges;
int			nmasks;
struct cinfo *		cip;
double			z;
double			real_z;
int			num_frac;

	cip = bbip -> cip;
	nedges	= cip -> num_edges;

	/* The special try_branch code for lp_solve can leave us with	*/
	/* an LP solution that isn't primal feasible.  In fact, it may	*/
	/* not even satisfy the variable bounds!  Detect this case	*/
	/* right up front.						*/

	for (i = 0; i < nedges; i++) {
		if (x [i] < -FUZZ) return (FALSE);
		if (x [i] > 1.0 + FUZZ) return (FALSE);
	}

	if (NOT integer_feasible_solution (x,
					   bbip -> tmap,
					   bbip -> fset_mask,
					   cip,
					   &num_frac)) {
		/* Not integer feasible -- get out. */
		return (FALSE);
	}

	/* We literally stumbled across an Integer Feasible Solution!	*/

	/* Re-calculate the final objective function (in order to	*/
	/* eliminate some of the LP solver's numerical errors)...	*/

	z = 0.0;
	for (i = 0; i < nedges; i++) {
		if (x [i] + FUZZ < 1.0) continue;
		z += cip -> cost [i];
	}

	/* Give caller the correct Z value... */
	*true_z = z;

	/* Damn Intel FPU is keeping 80 bits for Z... */
	/* Damn compilers keep getting smarter too!  ;^> */
	store_double (&real_z, z);

	if (real_z >= bbip -> best_z) {
		/* No better than what we have... */
		return (FALSE);
	}

	/* We have a new best integer feasible solution!  Record it!	*/
	nmasks = cip -> num_edge_masks;
	for (i = 0; i < nmasks; i++) {
		bbip -> smt [i] = 0;
	}
	for (i = 0; i < nedges; i++) {
		if (x [i] >= 0.5) {
			SETBIT (bbip -> smt, i);
		}
	}

	new_upper_bound (real_z, bbip);

	return (TRUE);
}

/*
 * This routine saves the current basis of the given LP.
 */

#ifdef CPLEX

	static
	void
save_LP_basis (

LP_t *			lp,		/* IN - LP to save basis for */
struct basis_save *	basp		/* OUT - saved basis info */
)
{
int		rows;
int		cols;

	rows = _MYCPX_getnumrows (lp);
	cols = _MYCPX_getnumcols (lp);

	basp -> cstat = NEWA (cols, int);
	basp -> rstat = NEWA (rows, int);

	if (_MYCPX_getbase (lp, basp -> cstat, basp -> rstat) NE 0) {
		fatal ("save_LP_basis: Bug 1.");
	}
}

/*
 * Destroy the saved basis info...
 */

	static
	void
destroy_LP_basis (

struct basis_save *	basp		/* IN - basis info to free up */
)
{
#if CPLEX >= 40
	free ((char *) (basp -> rstat));
	free ((char *) (basp -> cstat));
#endif
}

#endif

/*
 * This routine tries the given branch by solving the LP.  It
 * returns the resulting objective value, or INF_DISTANCE if something
 * goes wrong (like infeasible).
 */

#ifdef CPLEX

	static
	double
try_branch (

LP_t *			lp,		/* IN - LP to re-optimize */
int			var,		/* IN - variable to try branching */
int			dir,		/* IN - branch direction, 0 or 1 */
double *		x,		/* OUT - LP solution obtained */
double			ival,		/* IN - value to give if infeasible */
struct basis_save *	basp		/* IN - basis to restore when done */
)
{
int		status;
double		z;
int		b_index [2];
char		b_lu [2];
double		b_bd [2];

	--var;		/* vars are zero-origined in CPLEX... */

	b_index [0] = var;	b_lu [0] = 'L';
	b_index [1] = var;	b_lu [1] = 'U';
	if (dir EQ 0) {
		b_bd [0] = 0.0;
		b_bd [1] = 0.0;
	}
	else {
		b_bd [0] = 1.0;
		b_bd [1] = 1.0;
	}
	if (_MYCPX_chgbds (lp, 2, b_index, b_lu, b_bd) NE 0) {
		fatal ("try_branch: Bug 1.");
	}

	/* Solve the current LP instance... */
	status = _MYCPX_dualopt (lp);
	if (status NE 0) {
		tracef (" %%  WARNING dualopt: status = %d\n", status);
	}

	/* Get current LP solution... */
	if (_MYCPX_solution (lp, &status, &z, x, NULL, NULL, NULL) NE 0) {
		fatal ("try_branch: Bug 2.");
	}

	/* Determine type of LP result... */
	switch (status) {
	case CPX_OPTIMAL:
	case CPX_OPTIMAL_INFEAS:
		break;

	case CPX_INFEASIBLE:
	case CPX_UNBOUNDED:
			/* (CPLEX 3.0 sometimes gives us infeasible!) */
	case CPX_OBJ_LIM:	/* Objective limit exceeded in Phase II. */
		z = ival;
		break;

	default:
		tracef (" %% Status = %d\n", status);
		_MYCPX_lpwrite (lp, "core.lp");
		fatal ("try_branch: Bug 2.");
		break;
	}

	b_bd [0] = 0.0;
	b_bd [1] = 1.0;
	if (_MYCPX_chgbds (lp, 2, b_index, b_lu, b_bd) NE 0) {
		fatal ("try_branch: Bug 3.");
	}

	/* Restore the basis... */
	status = _MYCPX_copybase (lp, basp -> cstat, basp -> rstat);
	if (status NE 0) {
		fprintf (stderr, "try_branch: status = %d\n", status);
		fatal ("try_branch: Bug 4.");
	}

	return (z);
}

#endif

/*
 * This routine computes the lower-bound for the current node, which
 * consists of solving the LP and generating violated constraints
 * until either:
 *
 *	- LP becomes infeasible
 *	- LP objective meets or exceeds cutoff value
 *	- LP solution is integral
 *	- separation finds no more violated constraints
 */

	static
	int
compute_good_lower_bound (

struct bbinfo *		bbip		/* IN - branch-and-bound info */
)
{
int			i;
int			j;
int			n;
int			num_const;
int			status;
bitmap_t *		tmap;
bitmap_t *		fset_mask;
struct cinfo *		cip;
LP_t *			lp;
struct bbnode *		nodep;
double *		x;
double			z;
struct constraint *	cp;
struct constraint *	tmp;
int			iteration;
int			fix_status;
int			num_fractional;
cpu_time_t		Tlp;
cpu_time_t *		Tp;
bool			is_int;
char			tbuf1 [16];
char			tbuf2 [256];
char			time_str [256];
char			title [256];
cpu_time_t		Tn [20];

	cip	  = bbip -> cip;
	tmap	  = bbip -> tmap;
	fset_mask = bbip -> fset_mask;
	lp	  = bbip -> lp;
	nodep	  = bbip -> node;
	x	  = nodep -> x;

	Tp = &Tn [0];

	*Tp++ = get_cpu_time ();

	iteration = 1;
	num_const = 0;

	for (;;) {
		status = solve_LP_over_constraint_pool (bbip);
		z = nodep -> z;

		Tlp = get_cpu_time ();

		++(bbip -> node -> iter);
		++(bbip -> statp -> num_lps);

#if 0
		/* Display LP solution vector in machine-readable form... */
		for (i = 0; i < cip -> num_edges; i++) {
			tracef (" %% %08lx %08lx\n",
				((bitmap_t *) &x [i]) [0],
				((bitmap_t *) &x [i]) [1]);
		}
#endif

		if (Check_Root_Constraints AND
		    (bbip -> node -> depth EQ 0)) {
			fwrite (x,
				1,
				cip -> num_edges * sizeof (*x),
				rcfile);
		}

#if 1
		convert_cpu_time (Tlp - *--Tp, time_str);
		while (Tp > &Tn [0]) {
			--Tp;
			convert_cpu_time (Tp [1] - Tp [0], tbuf1);
			(void) sprintf (tbuf2, "%s/%s", tbuf1, time_str);
			strcpy (time_str, tbuf2);
		}
		(void) sprintf (title,
				"Node %d LP %d Solution, length = %f, %s %d",
				bbip -> node -> num, bbip -> node -> iter,
				z, time_str, num_const);
#if 0
		plot_lp_solution (title, x, cip, BIG_PLOT);
#else
		(void) tracef ("  %% %s\n", title);
#endif
#endif

		switch (status) {
		case BBLP_OPTIMAL:
			if (z >= bbip -> best_z) {
				nodep -> z = bbip -> best_z;
				return (LB_CUTOFF);
			}
			break;

		case BBLP_CUTOFF:
			nodep -> z = bbip -> best_z;
			return (LB_CUTOFF);

		case BBLP_INFEASIBLE:
			nodep -> z = bbip -> best_z;
			return (LB_INFEASIBLE);

		default:
			tracef ("%% solve status = %d\n", status);
			fatal ("compute_good_lower_bound: Bug 3.");
		}

#ifdef CPLEX
		/* Now get rid of any rows that have become	*/
		/* slack.  (We don't lose these constraints:	*/
		/* they're still sitting around in the		*/
		/* constraint pool.)				*/
		delete_slack_rows_from_LP (bbip);
#endif

		/* Solution is feasible, check for integer-feasible... */
		is_int = integer_feasible_solution (x,
						    tmap,
						    fset_mask,
						    cip,
						    &num_fractional);

		tracef (" %% %d fractional variables\n", num_fractional);

		if (is_int) {
			/* All vars are either 0 or 1 and the	*/
			/* solution is connected -- we have a	*/
			/* Steiner tree!			*/

			/* Re-calculate the final objective	*/
			/* function to eliminate numerical	*/
			/* errors in the value of Z...		*/
			z = 0.0;
			for (i = 0; i < cip -> num_edges; i++) {
				if (x [i] + FUZZ < 1.0) continue;
				z += cip -> cost [i];
			}
			nodep -> z = z;
			if (z >= bbip -> best_z) {
				/* probably a repeat performance... */
				nodep -> z = bbip -> best_z;
				return (LB_CUTOFF);
			}
			bbip -> node -> optimal = TRUE;
			return (LB_INTEGRAL);
		}

		/* Check to see if this node's objective value	*/
		/* is now high enough to be preempted...	*/
		if (nodep -> z > bbip -> preempt_z) {
			/* Node's value is no longer the lowest...	*/
			/* Preempt this one in favor of another.	*/
			return (LB_PREEMPTED);
		}

		/* Perhaps we have a new lower bound? */
		new_lower_bound (z, bbip);

		Tp = &Tn [0];
		*Tp++ = get_cpu_time ();

		compute_heuristic_upper_bound (x, bbip);

		/* If we have improved the upper bound, it is possible	*/
		/* that this node can now be cutoff...			*/
		if (nodep -> z >= bbip -> best_z) {
			nodep -> z = bbip -> best_z;
			return (LB_CUTOFF);
		}

		/* Try to fix some variables using reduced costs... */
		fix_status = reduced_cost_var_fixing (bbip);
		if (fix_status EQ VFIX_INFEASIBLE) {
			nodep -> z = bbip -> best_z;
			return (LB_INFEASIBLE);
		}
		if (fix_status EQ VFIX_FIXED_FRACTIONAL) {
			continue;
		}

		if (force_branch_flag) {
			/* User kicked us!  Stop separating and branch! */
			force_branch_flag = FALSE;
			break;
		}

		/* Apply all separation algorithms to solution... */
		cp = do_separations (bbip, &Tp);

		if (cp EQ NULL) {
			/* No more violated constraints found! */
			break;
		}

#ifdef LPSOLVE
		/* Now get rid of any rows that have become	*/
		/* slack.  (We don't lose these constraints:	*/
		/* they're still sitting around in the		*/
		/* constraint pool.)				*/
		delete_slack_rows_from_LP (bbip);
#endif

		/* Add new contraints to the constraint pool. */
		num_const = add_constraints (bbip, cp);

		if (num_const <= 0) {
			/* Separation routines found violations, but	*/
			/* the constraint pool disagrees...		*/
			fatal ("compute_good_lower_bound: Bug 4.");
		}

		while (cp NE NULL) {
			tmp = cp;
			cp = tmp -> next;
			free ((char *) (tmp -> mask));
			free ((char *) tmp);
		}
		++iteration;
	}

#if 1
	/* Print execution times of final iteration... */
	convert_cpu_time (0, time_str);
	--Tp;
	while (Tp > &Tn [0]) {
		--Tp;
		convert_cpu_time (Tp [1] - Tp [0], tbuf1);
		(void) sprintf (tbuf2, "%s/%s", tbuf1, time_str);
		strcpy (time_str, tbuf2);
	}
	(void) tracef ("  %% Final iteration: %s\n", time_str);
#endif

	/* Only get here with fractional solution and	*/
	/* no more violated constraints were found.	*/

	return (LB_FRACTIONAL);
}

/*
 * Routine to print out an updated lower bound value.
 */

	static
	void
new_lower_bound (

double			lb,		/* IN - new lower bound value */
struct bbinfo *		bbip		/* IN - branch-and-bound info */
)
{
int			i;
struct cinfo *		cip;
struct scale_info *	sip;
double			prev;
double			old_gap;
double			new_gap;
cpu_time_t		t1;
char			buf1 [32];

	cip = bbip -> cip;
	sip = &(cip -> scale);

	prev = bbip -> prevlb;
	if (lb <= prev) {
		/* There has been no improvement - get out. */
		return;
	}

	if (prev <= -DBL_MAX) {
		/* Don't make lower bound jump from initial value... */
		prev = lb;
	}

	/* Print out the old and new lower bounds, with timestamp. */
	t1 = get_cpu_time ();
	convert_cpu_time (t1 - bbip -> t0, buf1);

	if ((bbip -> best_z >= DBL_MAX) OR (bbip -> best_z EQ 0.0)) {
		old_gap = 99.9;
		new_gap = 99.9;
	}
	else {
		old_gap = 100.0 * (bbip -> best_z - prev) / bbip -> best_z;
		new_gap = 100.0 * (bbip -> best_z - lb) / bbip -> best_z;
	}

	tracef (" %% @LO %s %24.20f %2.10f\n",
		buf1, UNSCALE (prev, sip), old_gap);
	tracef (" %% @LN %s %24.20f %2.10f\n",
		buf1, UNSCALE (lb, sip), new_gap);

	bbip -> prevlb = lb;
}

/*
 * Routine to print out an updated upper bound value.
 */

	void
new_upper_bound (

double			ub,		/* IN - new upper bound value */
struct bbinfo *		bbip		/* IN - branch-and-bound info */
)
{
struct cinfo *		cip;
struct scale_info *	sip;
double			prev;
int			i;
cpu_time_t		t1;
double			old_gap;
double			new_gap;
char			buf1 [64];
char			buf2 [64];

	cip = bbip -> cip;
	sip = &(cip -> scale);

	prev = bbip -> best_z;
	if (ub >= prev) {
		/* Supposed to be an improvement! */
		fatal ("new_upper_bound: Bug 1.");
	}

	/* We HAVE a new best solution! */
	bbip -> best_z = ub;

#if CPLEX
	{ double toobig, toosmall, ulim;
	  /* Set new cutoff value for future LPs... */
	  ulim = ldexp (ub, -(bbip -> lpmem -> obj_scale));
	  if (_MYCPX_setobjulim (ulim, &toosmall, &toobig) NE 0) {
		fatal ("new_upper_bound: Bug 2.");
	  }
	}
#endif

#if LPSOLVE
	/* Set new cutoff value for future LPs... */
	/* (This may not really work in lp_solve.) */
	bbip -> lp -> obj_bound = ub;
#endif

	cut_off_existing_nodes (ub, bbip -> bbtree);

	/* Might want to do this if all other nodes were cut off. */
	update_node_preempt_value (bbip);

	/* Now print out the trace messages. */
	if (prev >= DBL_MAX) {
		/* Don't make upper bound jump from infinity... */
		prev = ub;
	}

	/* Print out the old and new lower and upper bounds, with timestamp. */
	t1 = get_cpu_time ();
	convert_cpu_time (t1 - bbip -> t0, buf1);

	if (bbip -> prevlb <= -DBL_MAX) {
		old_gap = 99.9;
		new_gap = 99.9;
	}
	else {
		old_gap = 100.0 * (prev - bbip -> prevlb) / prev;
		new_gap = 100.0 * (ub - bbip -> prevlb) / ub;
	}

	tracef (" %% @UO %s %24.20f %2.10f\n",
		buf1, UNSCALE (prev, sip), old_gap);
	tracef (" %% @UN %s %24.20f %2.10f\n",
		buf1, UNSCALE (ub, sip), new_gap);
}

/*
 * This routine deletes any existing node whose objective value is
 * cut off by the given latest feasible integer solution.
 */

	static
	void
cut_off_existing_nodes (

double		best_z,		/* IN - new best objective value */
struct bbtree *	tp		/* IN - branch-and-bound tree */
)
{
int		num_cut;
struct bbheap *	hp;
struct bbnode *	p;

	num_cut = 0;

	/* We process the nodes from WORST to best... */
	hp = &(tp -> heap [WORST_NODE_HEAP]);

	while (hp -> nheap > 0) {
		/* Get node with highest objective value... */
		p = hp -> array [0];
		if (p -> index [WORST_NODE_HEAP] NE 0) {
			fatal ("cut_off_existing_nodes: Bug 1.");
		}
		if (p -> z < best_z) {
			/* All remaining nodes are < best_z... */
			break;
		}

		/* This node has been cut off! */
		delete_node_from_bbtree (p, tp);
		p -> next = tp -> free;
		tp -> free = p;
		++num_cut;
	}

	if (num_cut > 0) {
		tracef (" %%	=== %d nodes cut off ===\n", num_cut);
	}
}

/*
 * This routine updates the node preemption value.  The node to be
 * preempted must be active when this routine is called (i.e., it must
 * be removed from the heap so that the "next best" node is at the top
 * of the heap).
 */

	static
	void
update_node_preempt_value (

struct bbinfo *		bbip		/* IN - branch-and-bound info */
)
{
struct bbtree *		tp;
struct bbheap *		hp;
struct bbnode *		node2;

	tp = bbip -> bbtree;
	hp = &(tp -> heap [BEST_NODE_HEAP]);

	if (hp -> nheap <= 0) {
		/* No other nodes.  Preempt only at cutoff. */
		bbip -> preempt_z = bbip -> best_z;
	}
	else {
		/* Preempt current node when next-best is exceeded. */
		node2 = hp -> array [0];
		bbip -> preempt_z = node2 -> z;
	}
}

/*
 * This routine performs most of the separations -- in the proper order.
 */

	static
	struct constraint *
do_separations (

struct bbinfo *		bbip,		/* IN - branch-and-bound info */
cpu_time_t **		Tpp		/* IN/OUT - CPU time vector */
)
{
double *		x;
bitmap_t *		vert_mask;
bitmap_t *		edge_mask;
struct cinfo *		cip;
struct comp *		comp;
struct comp *		p;
struct comp **		hookp;
struct comp *		p2;
cpu_time_t *		Tp;
struct constraint *	cp;
struct constraint *	cp2;
struct constraint *	tmp;
bool			print_flag;
bool			optimal;
bool			any_skipped;

	x		= bbip -> node -> x;
	vert_mask	= bbip -> tmap;
	edge_mask	= bbip -> fset_mask;
	cip		= bbip -> cip;

	Tp = *Tpp;

	/* Find all zero-weight cutsets... */
	cp = find_cutset_constraints (x, vert_mask, edge_mask, cip);
	*Tp++ = get_cpu_time ();

	/* Find solid integer cycles... */
	cp = find_integer_cycles (x, vert_mask, edge_mask, cp, cip);
	*Tp++ = get_cpu_time ();

	/* Break problem up into congested components... */
	print_flag = TRUE;
	comp = find_congested_components (x,
					  vert_mask,
					  edge_mask,
					  print_flag,
					  cip);

	/* Exhaustively enumerate all components that are sufficiently	*/
	/* small...  Delete them from the list when done.		*/
	hookp = &comp;
	while ((p = *hookp) NE NULL) {
		if (p -> num_verts <= SEC_ENUM_LIMIT) {
			cp2 = enumerate_all_subtours (p, NULL, bbip);
			*hookp = p -> next;
			p -> next = NULL;
			free_congested_component (p);
			while (cp2 NE NULL) {
				tmp = cp2 -> next;
				cp2 -> next = cp;
				cp = cp2;
				cp2 = tmp;
			}
		}
		else {
			hookp = &(p -> next);
		}
	}
	*Tp++ = get_cpu_time ();

	/* Find violated SEC's using a heuristic flow	*/
	/* formulation.					*/
	for (p = comp; p NE NULL; p = p -> next) {
		p -> cp = sec_flow_heuristic (p,
					      x,
					      edge_mask,
					      cip,
					      p -> cp);
	}
	*Tp++ = get_cpu_time ();

	/* Find small-cardinality subtour violations	*/
	/* by partial enumeration...			*/
	for (p = comp; p NE NULL; p = p -> next) {
		p -> cp = find_small_subtours (p, p -> cp, bbip);
	}
	*Tp++ = get_cpu_time ();

	/* Discard each component for which we have found at least one	*/
	/* violation.  Gather all constraints onto the main list...	*/
	hookp = &comp;
	for (;;) {
		p = *hookp;
		if (p EQ NULL) break;
		cp2 = p -> cp;
		if (cp2 NE NULL) {
			/* Gather these constraints onto main list... */
			p -> cp = NULL;
			while (cp2 NE NULL) {
				tmp = cp2 -> next;
				cp2 -> next = cp;
				cp = cp2;
				cp2 = tmp;
			}
			*hookp = p -> next;
			p -> next = NULL;
			free_congested_component (p);

		}
		else {
			/* No constraints yet for this component.  We	*/
			/* may want to try the more expensive method...	*/
			/* Retain this component.			*/
			hookp = &(p -> next);
		}
	}

#if 0
	/* Time to use the new-fangled SEC separator... */
	p2 = comp;
	comp = NULL;
	cp = sec_flow_separator (&p2, x, edge_mask, bbip, cp);
	*Tp++ = get_cpu_time ();
#else
	/* Time to use the new-fangled SEC separator... */
	/* Do it one component at a time, so that we can see if */
	/* there are any components for which no violations were found. */
	while (comp NE NULL) {
		p2 = comp;
		comp = comp -> next;
		p2 -> next = NULL;
		cp2 = sec_flow_separator (&p2, x, edge_mask, bbip, NULL);
		while (cp2 NE NULL) {
			tmp = cp2 -> next;
			cp2 -> next = cp;
			cp = cp2;
			cp2 = tmp;
		}
	}
	*Tp++ = get_cpu_time ();
#endif

	/* If this separation routine does not find any SEC violations,	*/
	/* it means that none exist!					*/
	optimal = TRUE;

	/* Note: congested components are all freed now... */

#if 0
	if (cp EQ NULL) {
		/* Nothing else found -- look for fractional cutsets... */
		cp = find_fractional_cutsets (x,
					      bbip -> csip,
					      vert_mask,
					      edge_mask,
					      cip);
		*Tp++ = get_cpu_time ();
	}
#endif

	if ((cp EQ NULL) AND optimal) {
		/* We KNOW that we have NO violations!  The LP	*/
		/* relaxation is now OPTIMAL!			*/
		bbip -> node -> optimal = TRUE;
	}

	*Tpp = Tp;

	return (cp);
}

/*
 * This routine attempts to use LP reduced costs to fix variables.  Any
 * variable whose reduced cost exceeds the current LP/IP gap can be
 * permanently fixed to its current value (either 0 or 1) for the duration
 * of the current bb-node (and all of its children).
 */

	static
	int
reduced_cost_var_fixing (

struct bbinfo *		bbip	/* IN - branch-and-bound info */
)
{
int			i;
int			nedges;
int			nmasks;
int			status;
int			nfix0;
int			nfix1;
struct cinfo *		cip;
LP_t *			lp;
double *		zlb;
double			gap;
double			threshold;
int *			newfix0;
int *			newfix1;
struct bbnode *		nodep;
double *		x;

	if (bbip -> best_z >= DBL_MAX) {
		/* Can't reasonably attempt this until we have	*/
		/* a valid upper bound...			*/
		return (VFIX_NOTHING_FIXED);
	}

	nodep = bbip -> node;
	gap = bbip -> best_z - nodep -> z;

	/* Only fix if we significantly exceed the gap... */
	gap *= (1.0 + FUZZ);

	threshold = nodep -> z + gap;

	lp	= bbip -> lp;
	cip	= bbip -> cip;
	nedges	= cip -> num_edges;
	nmasks	= cip -> num_edge_masks;

	nodep	= bbip -> node;
	x	= nodep -> x;
	zlb	= nodep -> zlb;

	newfix0 = NEWA (nedges, int);
	newfix1 = NEWA (nedges, int);

	nfix0	= 0;
	nfix1	= 0;

	for (i = 0; i < nedges; i++) {
		if (BITON (bbip -> fixed, i)) continue;

		if (zlb [2 * i] > threshold) {
			newfix1 [nfix1++] = i;
		}
		if (zlb [2 * i + 1] > threshold) {
			newfix0 [nfix0++] = i;
		}
	}

	status = VFIX_NOTHING_FIXED;

	if ((nfix0 > 0) OR (nfix1 > 0)) {
		status = fix_variables (bbip, newfix0, nfix0, newfix1, nfix1);
	}

	free ((char *) newfix1);
	free ((char *) newfix0);

	return (status);
}

/*
 * This routine fixes variables to zero and/or one.  We are given two
 * lists of variables to fixed, those to be fixed to zero, and those
 * to be fixed to one.  This routine then iteratively applies a series
 * of deductive steps that can cause additional variables to be fixed
 * based upon connectivity and compatibility criteria.
 */

	static
	int
fix_variables (

struct bbinfo *		bbip,		/* IN - branch-and-bound info */
int *			fix_to_0,	/* IN - vars to fix to 0 */
int			nfix0,		/* IN - number of vars fixed to 0 */
int *			fix_to_1,	/* IN - vars to fix to 1 */
int			nfix1		/* IN - number of vars fixed to 1 */
)
{
int			i;
int			j;
int			k;
int			t;
int			fs;
int			nedges;
int			nmasks;
int			kmasks;
int			status;
int			last_fset;
int			fix0_count;
int			fix1_count;
struct cinfo *		cip;
struct bbnode *		nodep;
bitmap_t *		fixmask0;
bitmap_t *		fixmask1;
bitmap_t *		terms_checked;
int *			ep1;
int *			ep2;
int *			ep3;
int *			ep4;
int *			vp1;
int *			vp2;
int			vars_fixed;
int			fix_frac;

#undef	PRINT_FIXED_VARIABLES

	cip	= bbip -> cip;
	nodep	= bbip -> node;

	nedges	= cip -> num_edges;
	nmasks	= cip -> num_edge_masks;
	kmasks	= cip -> num_vert_masks;

	fixmask0 	= NEWA (2 * nmasks + kmasks, bitmap_t);
	fixmask1 	= fixmask0 + nmasks;
	terms_checked	= fixmask1 + nmasks;

	/* Initialize masks of vars in fix-to-0 and fix-to-1 lists. */
	/* We use these to prevent adding duplicate entries later on. */
	for (i = 0; i < nmasks; i++) {
		fixmask0 [i] = 0;
		fixmask1 [i] = 0;
	}
	for (i = 0; i < nfix0; i++) {
		SETBIT (fixmask0, fix_to_0 [i]);
	}
	for (i = 0; i < nfix1; i++) {
		SETBIT (fixmask1, fix_to_1 [i]);
	}

	status = VFIX_NOTHING_FIXED;

	fix_frac = 0;

	fix0_count = 0;
	fix1_count = 0;

	/* Iteratively fix variables until no more fixing can be done. */
	do {
		vars_fixed = FALSE;

		/* =============== Handle fixing vars to 0 =============== */

		ep1 = fix_to_0;
		ep2 = ep1;
		ep3 = ep1 + nfix0;
		while (ep1 < ep3) {
			fs = *ep1++;
			CLRBIT (fixmask0, fs);
			if (NOT BITON (bbip -> fset_mask, fs)) continue;
			if (BITON (bbip -> fixed, fs)) {
				if (NOT BITON (bbip -> value, fs)) {
					/* Already fixed to zero!  Ignore. */
					continue;
				}
				/* Already fixed to one! */
				status = VFIX_INFEASIBLE;
				goto alldone;
			}
			/* Fix it to zero now! */
			SETBIT (bbip -> fixed, fs);
			CLRBIT (bbip -> value, fs);
			SETBIT (nodep -> fixed, fs);
			CLRBIT (nodep -> value, fs);
			if ((FUZZ < nodep -> x [fs]) AND
			    (nodep -> x [fs] + FUZZ < 1.0)) {
				++fix_frac;
			}
			change_var_bounds (bbip -> lp, fs, 0.0, 0.0);
			++fix0_count;
#ifdef PRINT_FIXED_VARIABLES
			tracef (" %%	Fixed x%-3d = 0\n", fs);
#endif
			/* Save this edge for later check. */
			*ep2++ = fs;
			/* We have fixed at least 1 variable! */
			status = VFIX_VARIABLES_FIXED;
			vars_fixed = TRUE;
		}
		ep1 = fix_to_0;
		if (ep1 < ep2) {
			/* Check if any of the vars we just set to 0	*/
			/* permit us to deduce vars that must be 1...	*/
			for (i = 0; i < kmasks; i++) {
				terms_checked [i] = 0;
			}
			while (ep1 < ep2) {
				i = *ep1++;
				vp1 = cip -> edge [i];
				vp2 = cip -> edge [i + 1];
				while (vp1 < vp2) {
					t = *vp1++;
					if (BITON (terms_checked, t)) continue;
					SETBIT (terms_checked, t);
					ep3 = cip -> term_trees [t];
					ep4 = cip -> term_trees [t + 1];
					k = 0;
					last_fset = -1;
					while (ep3 < ep4) {
						fs = *ep3++;
						if (NOT BITON (bbip -> fset_mask, fs)) continue;
						if (BITON (bbip -> fixed, fs) AND
						    NOT BITON (bbip -> value, fs)) continue;
						/* Full set fs has term t */
						/* and is NOT fixed to zero. */
						last_fset = fs;
						++k;
						if (k > 1) break;
					}
					if (k <= 0) {
						/* disconnected terminal! */
						status = VFIX_INFEASIBLE;
						goto alldone;
					}
					if ((k EQ 1) AND
					    NOT BITON (fixmask1, last_fset)) {
						/* one full set left, it */
						/* must be taken! */
						SETBIT (fixmask1, last_fset);
						fix_to_1 [nfix1++] = last_fset;
					}
				}
			}
		}
		/* Set the Fix-to-0 list to empty.  Fixmask0 should now	*/
		/* have all bits turned off.				*/
		nfix0 = 0;


		/* =============== Handle fixing vars to 1 =============== */

		ep1 = fix_to_1;
		ep2 = ep1 + nfix1;
		while (ep1 < ep2) {
			fs = *ep1++;
			CLRBIT (fixmask1, fs);
			if (NOT BITON (bbip -> fset_mask, fs)) continue;
			if (BITON (bbip -> fixed, fs)) {
				if (BITON (bbip -> value, fs)) {
					/* Already fixed to one!  Ignore. */
					continue;
				}
				/* Already fixed to zero! */
				status = VFIX_INFEASIBLE;
				goto alldone;
			}
			/* Fix it to one now! */
			SETBIT (bbip -> fixed, fs);
			SETBIT (bbip -> value, fs);
			SETBIT (nodep -> fixed, fs);
			SETBIT (nodep -> value, fs);
			if ((FUZZ < nodep -> x [fs]) AND
			    (nodep -> x [fs] + FUZZ < 1.0)) {
				++fix_frac;
			}
			change_var_bounds (bbip -> lp, fs, 1.0, 1.0);
			++fix1_count;
#ifdef PRINT_FIXED_VARIABLES
			tracef (" %%	Fixed x%-3d = 1\n", fs);
#endif
			/* We have fixed at least 1 variable! */
			status = VFIX_VARIABLES_FIXED;
			vars_fixed = TRUE;

			/* Fix every *other* incompatible FST to zero! */
			ep3 = cip -> inc_edges [fs];
			ep4 = cip -> inc_edges [fs + 1];
			while (ep3 < ep4) {
				j = *ep3++;
				if (j EQ fs) continue;
				if (NOT BITON (fixmask0, j)) {
					SETBIT (fixmask0, j);
					fix_to_0 [nfix0++] = j;
				}
			}
		}
		/* Set the Fix-to-1 list to empty.  Fixmask1 should now	*/
		/* have all bits turned off.				*/
		nfix1 = 0;
	} while (vars_fixed);

alldone:

	if ((fix0_count | fix1_count) NE 0) {
		/* Problem has changed -- force re-solve of LP. */
		nodep -> cpiter = -1;
	}

	switch (status) {
	case VFIX_NOTHING_FIXED:
		break;

	case VFIX_VARIABLES_FIXED:
		if (fix_frac > 0) {
			tracef (" %% Fixed %d vars to 0 and %d vars to 1 (%d were fractional).\n",
				fix0_count, fix1_count, fix_frac);
			status = VFIX_FIXED_FRACTIONAL;
		}
		else {
			tracef (" %% Fixed %d vars to 0 and %d vars to 1.\n",
				fix0_count, fix1_count);
		}
		break;

	case VFIX_INFEASIBLE:
		tracef (" %% Variable fixing detected infeasibility!\n");
		break;

	default:
		fatal ("fix_variables: Bug 1.");
		break;
	}

	free ((char *) fixmask0);

	return (status);
}

/*
 * This routine changes the bounds on the given LP variable...
 */

	static
	void
change_var_bounds (

LP_t *			lp,		/* IN - LP to changes bounds of */
int			var,		/* IN - variable to fix */
double			lower,		/* IN - lower bound */
double			upper		/* IN - upper bound */
)
{
#if CPLEX
int			b_index [2];
char			b_lu [2];
double			b_bd [2];

	b_index [0] = var;	b_lu [0] = 'L';		b_bd [0] = lower;
	b_index [1] = var;	b_lu [1] = 'U';		b_bd [1] = upper;
	if (_MYCPX_chgbds (lp, 2, b_index, b_lu, b_bd) NE 0) {
		fatal ("change_var_bounds: Bug 1.");
	}
#endif

#if LPSOLVE
	set_bounds (lp, var + 1, lower, upper);
#endif
}

/*
 * This routine checks to see if we have an integer feasible solution.
 * First we check for integrality, then we check connectedness.
 */

	static
	bool
integer_feasible_solution (

double *		x,		/* IN - LP solution to check. */
bitmap_t *		tmap,		/* IN - subset of terminals. */
bitmap_t *		fset_mask,	/* IN - subset of full-sets. */
struct cinfo *		cip,		/* IN - compatibility info. */
int *			num_fractional	/* OUT - number of fractional vars. */
)
{
int			i;
int			j;
int			t;
int			fs;
int			fs2;
int			kmasks;
int			nedges;
int			nmasks;
int			num_int;
int			num_frac;
int			starting_fset;
struct pset *		terms;
bitmap_t *		integral_fsets;
int *			ep1;
int *			ep2;
int *			vp1;
int *			vp2;
int *			sp;
int *			stack;
bitmap_t		mask;
bitmap_t *		terms_left;

	nedges = cip -> num_edges;
	nmasks = cip -> num_edge_masks;
	kmasks = cip -> num_vert_masks;

	/* First do a quick check of the integrality of the solution... */
	num_frac	= 0;
	for (i = 0; i < nedges; i++) {
#if 0
		/* Disabling this check because the heuristic upper-	*/
		/* bound routine sometimes produces solutions that	*/
		/* include edges that are NOT in the fset_mask!  For	*/
		/* example consider the case when an MST edge gets	*/
		/* pruned.  The heuristic may use it to glue the final	*/
		/* pieces together.  We want to consider the solution	*/
		/* valid if it forms a tree, even if it uses edges that	*/
		/* we know are suboptimal!				*/
		if (NOT BITON (fset_mask, i)) continue;
#endif
		if (x [i] <= FUZZ) continue;
		if (x [i] + FUZZ >= 1.0) continue;

		/* Variable has a fractional value... */
		++num_frac;
	}
	*num_fractional = num_frac;

	if (num_frac > 0) {
		/* We have fractional variables -- solution is	*/
		/* definitely NOT integer feasible...		*/
		return (FALSE);
	}

	/* All solution variables are either 0 or 1 -- integral!  This	*/
	/* case is much less common.  Loop back over the solution,	*/
	/* making a note of all full sets present in the solution.	*/

	integral_fsets = NEWA (nmasks, bitmap_t);
	for (i = 0; i < nmasks; i++) {
		integral_fsets [i] = 0;
	}

	num_int		= 0;
	starting_fset	= -1;
	j = 0;
	for (i = 0; i < nedges; i++) {
#if 0
		/* Disabling this check because the heuristic upper-	*/
		/* bound routine sometimes produces solutions that	*/
		/* include edges that are NOT in the fset_mask!  For	*/
		/* example consider the case when an MST edge gets	*/
		/* pruned.  The heuristic may use it to glue the final	*/
		/* pieces together.  We want to consider the solution	*/
		/* valid if it forms a tree, even if it uses edges that	*/
		/* we know are suboptimal!				*/
		if (NOT BITON (fset_mask, i)) continue;
#endif
		if (x [i] >= 0.5) {
			SETBIT (integral_fsets, i);
			starting_fset = i;
			++num_int;
			j += (cip -> edge_size [i] - 1);
		}
	}

	if (j NE cip -> num_verts - 1) {
		/* Wrong cardinality of edges -- cannot be a tree. */
		free ((char *) integral_fsets);
		return (FALSE);
	}

	if (starting_fset < 0) {
		/* No full sets in solution -- problem must have one or	*/
		/* fewer terminals!  This is connected by default.	*/
		free ((char *) integral_fsets);
		return (TRUE);
	}

	/* Create temporary mask of terminals we have not yet seen... */
	terms_left = NEWA (kmasks, bitmap_t);
	for (i = 0; i < kmasks; i++) {
		terms_left [i] = tmap [i];
	}

	stack = NEWA (num_int, int);
	sp = stack;

	/* Find connected component containing the starting_fset... */
	CLRBIT (integral_fsets, starting_fset);
	--num_int;
	*sp++ = starting_fset;

	while (sp > stack) {
		fs = *--sp;
		vp1 = cip -> edge [fs];
		vp2 = cip -> edge [fs + 1];
		while (vp1 < vp2) {
			t = *vp1++;
			if (NOT BITON (terms_left, t)) continue;
			CLRBIT (terms_left, t);
			ep1 = cip -> term_trees [t];
			ep2 = cip -> term_trees [t + 1];
			while (ep1 < ep2) {
				fs2 = *ep1++;
				if (NOT BITON (integral_fsets, fs2)) continue;
				CLRBIT (integral_fsets, fs2);
				--num_int;
				*sp++ = fs2;
			}
		}
	}

	/* See if any terminals were not reached... */
	mask = 0;
	for (i = 0; i < kmasks; i++) {
		mask |= terms_left [i];
	}

	free ((char *) stack);
	free ((char *) terms_left);
	free ((char *) integral_fsets);

	if (mask NE 0) {
		/* At least one more connected component -- solution	*/
		/* is not connected, and therefore infeasible!  (We	*/
		/* also know we have at least one integer cycle!)	*/
		return (FALSE);
	}

	/* Solution is a Steiner tree!  (Not necessarily minimal.) */

	return (TRUE);
}

/*
 * This routine goes through each of the constraints in the LP tableaux
 * for the root node.  Each constraint is checked against each LP
 * solution (recorded in the rcfile) to find the earliest iteration in
 * which we could have had root LP optimality -- if our separation
 * algorithm magically generated the *right* constraints.
 */

	static
	void
check_root_constraints (

struct bbinfo *		bbip	/* IN - branch-and-bound info */
)
{
int			i;
int			j;
int			nedges;
int			nbytes;
int			iter;
int *			list;
struct cinfo *		cip;
struct cpool *		pool;
int *			ip1;
int *			ip2;
int *			ip3;
double *		x;

	cip	= bbip -> cip;
	pool	= bbip -> cpool;

	/* Develop the list of all (non-initial) constraints... */
	list = NEWA (pool -> nlprows, int);
	ip3 = list;
	for (i = 0; i < pool -> nlprows; i++) {
		j = pool -> lprows [i];
		if (j < pool -> initrows) continue;
		*ip3++ = j;
	}


	fclose (rcfile);

	rcfile = fopen ("/tmp/lp.x", "r");

	nedges = cip -> num_edges;
	nbytes = nedges * sizeof (double);

	x = NEWA (nedges, double);

	iter = 0;
	while (ip3 > list) {
		i = fread (x, 1, nbytes, rcfile);
		if (i NE nbytes) {
			fatal ("check_root_constraints: Bug 1.");
		}

		/* Delete all remaining constraints that violate x. */
		ip1 = ip2 = list;
		while (ip2 < ip3) {
			j = *ip2++;
			if (NOT is_violation (pool -> rows [j].coefs, x)) {
				/* No violation -- keep constraint around. */
				*ip1++ = j;
			}
		}
		ip3 = ip1;
		tracef (" %% @r iter %d, %d constraints left\n",
			iter, ip3 - list);
		if (ip3 <= list) break;
		++iter;
	}

	tracef (" %% @RC Could have gotten root constraints in %d iterations!\n", iter);

	fclose (rcfile);
	rcfile = NULL;

	free ((char *) x);
	free ((char *) list);
}
