/***********************************************************************

	File:	bbsubs.c
	Rev:	a-1
	Date:	02/28/2001

	Copyright (c) 1995, 2001 by David M. Warme

************************************************************************

	Low-level subtroutines of the branch-and-cut.

************************************************************************

	Modification Log:

	a-1:	02/28/2001	warme
		: Split off from bb.c.  Changes for 3.1 release.

************************************************************************/

#include "bb.h"
#include "bbsubs.h"
#include "config.h"
#include "constrnt.h"
#include "steiner.h"
#include "ub.h"


/*
 * Global Routines
 */

void			add_bbnode (struct bbinfo *, int, int, double);
void			append_node_to_tree (struct bbnode *, struct bbtree *);
void			bbheap_insert (struct bbnode *, struct bbtree *, int);
struct bbtree *		create_bbtree (int);
void			delete_node_from_bbtree (struct bbnode *,
						 struct bbtree *);
void			destroy_bbinfo (struct bbinfo *);
void			shutdown_lp_solver (void);
void			startup_lp_solver (void);


#if defined(CPLEX)
 #if CPLEX >= 40
  CPXENVptr		cplex_env;
 #endif
#endif


/*
 * Local Routines
 */

static void		bbheap_delete (struct bbnode *, struct bbtree *, int);
static void		bbheap_free (struct bbheap *);
static void		bbheap_init (struct bbheap *, bbheap_func_t *);
static void		destroy_bbnode (struct bbnode *);
static int		node_is_better (struct bbnode *, struct bbnode *);
static int		node_is_worse (struct bbnode *, struct bbnode *);

#ifdef CPLEX
static void		startup_cplex (void);
#endif

/*
 * This routine does everything needed to start up whichever
 * LP solver we are using.
 */

	void
startup_lp_solver (void)

{
#ifdef CPLEX
	startup_cplex ();
#endif

	/* Nothing for lp_solve... */
}



/*
 * This routine does everything needed to shut down whichever
 * LP solver we are using.
 */

	void
shutdown_lp_solver (void)

{
#if defined (CPLEX) AND (CPLEX >= 40)
	/* Shut down CPLEX... */
	if (CPXcloseCPLEX (&cplex_env) NE 0) {
		fprintf (stderr, "Warning: Unable to close CPLEX.\n");
	}
#endif

	/* Nothing for lp_solve... */
}

/*
 * This routine performs everything needed to start up the
 * newer versions of CPLEX.
 */

#if defined(CPLEX) AND (CPLEX >= 40)

	static
	void
startup_cplex (void)

{
int		status;
FILE *		fp;
CPXCHANNELptr	_res, _warn, _error, _log;
char		msg [512];

#ifdef HAVE_STDERR_IS_LVALUE
    {	FILE * esave;
	/* Flush the #@!$% CPLEX startup banner! */
	esave = stderr;
	stderr = fopen ("/dev/null", "w");
#endif

	cplex_env = _MYCPX_openCPLEX (&status);

#ifdef HAVE_STDERR_IS_LVALUE
	/* Undo flushing of stderr... */
	fp = stderr;
	stderr = esave;
	fclose (fp);
    }
#endif

	if (cplex_env EQ NULL) {
		if (CPXgeterrorstring (NULL, status, msg) EQ NULL) {
			strcpy (msg, "No CPLEX error message.");
		}
		fprintf (stderr, "%s\n", msg);
		goto shutdown;
	}

	/* Get rid of CPLEX's default message destinations. */
	CPXgetchannels (cplex_env, &_res, &_warn, &_error, &_log);
	CPXdisconnectchannel (cplex_env, _res);
	CPXdisconnectchannel (cplex_env, _warn);
	CPXdisconnectchannel (cplex_env, _error);
	CPXdisconnectchannel (cplex_env, _log);

	CPXsetintparam (cplex_env, CPX_PARAM_SCRIND, 0);

	fp = fopen ("cplex.log", "a");
	if (fp EQ NULL) {
		perror ("cplex.log");
		goto shutdown;
	}

	/* Send all log stuff to the cplex.log file. */
	CPXsetlogfile (cplex_env, fp);

	/* But discard the results stuff... */
	CPXdisconnectchannel (cplex_env, _res);

	/* CPLEX is now ready to roll! */
	return;

	/* Something didn't work -- shut down CPLEX and get out. */
shutdown:
	CPXcloseCPLEX (&cplex_env);
	exit (1);
}

#endif

/*
 * This routine performs everything needed to start up older
 * versions of CPLEX.
 */

#if defined(CPLEX) AND (CPLEX < 40)

	static
	void
startup_cplex (void)

{
	/* Older CPLEX library isn't as noisy! */
	_MYCPX_setlogfile (stderr);
}

#endif

/*
 * This routine adds a new node to the branch-and-bound tree.
 */

	void
add_bbnode (

struct bbinfo *		bbip,		/* IN - branch-and-bound info */
int			var,		/* IN - variable to branch on */
int			dir,		/* IN - branch direction */
double			z		/* IN - value to give to node */
)
{
int			i;
int			j;
int			n;
int			nedges;
int			nmasks;
struct bbnode *		parent;
struct bbtree *		tp;
struct bbnode *		p;
struct bbnode *		p1;
struct bbnode *		p2;
struct bbnode **	tmp1;
struct bbnode **	tmp2;

	if (z >= bbip -> best_z) {
		/* Node is already worse than the cutoff value!	*/
		/* Don't bother adding it to the tree!		*/
		return;
	}

	parent	= bbip -> node;		/* current node becomes parent */
	tp	= bbip -> bbtree;	/* tree to insert in */

	tracef (" %% @NC %4d %4d	x%d = %d	%f\n",
		tp -> snum, parent -> num, var, dir, z);

	nmasks	= tp -> nmasks;
	nedges	= bbip -> cip -> num_edges;

	/* Get a new tree node... */
	p = tp -> free;
	if (p NE NULL) {
		tp -> free = p -> next;
	}
	else {
		p = NEW (struct bbnode);
		p -> x	   = NEWA (nedges, double);
		p -> zlb   = NEWA (2 * nedges, double);
		p -> fixed = NEWA (nmasks, bitmap_t);
		p -> value = NEWA (nmasks, bitmap_t);
		p -> bheur = NEWA (nedges, double);
	}
	p -> z		= z;
	p -> optimal	= FALSE;
	p -> num	= (tp -> snum)++;
	p -> iter	= 0;
	p -> parent	= parent -> num;
	p -> var	= var;
	p -> dir	= dir;
	p -> depth	= 1 + parent -> depth;
	if (dir EQ 0) {
		p -> br1cnt = parent -> br1cnt;
	}
	else {
		p -> br1cnt = parent -> br1cnt + 1;
	}
	p -> cpiter = -1;	/* force re-solve of LP. */

	/* Get most up-to-date fixed variables... */
	for (i = 0; i < nmasks; i++) {
		p -> fixed [i] = bbip -> fixed [i];
		p -> value [i] = bbip -> value [i];
	}
	SETBIT (p -> fixed, var);
	if (dir EQ 0) {
		CLRBIT (p -> value, var);
	}
	else {
		SETBIT (p -> value, var);
	}
	p -> n_uids	= 0;
	p -> bc_uids	= NULL;
	p -> bc_row	= NULL;
	p -> rstat	= NULL;
	p -> cstat	= NULL;

	/* Copy parent's LP solution, etc. */
	memcpy (p -> x, parent -> x, nedges * sizeof (p -> x [0]));
	memcpy (p -> zlb, parent -> zlb, nedges * (2 * sizeof (p -> zlb [0])));
	memcpy (p -> bheur, parent -> bheur, nedges * sizeof (p -> bheur [0]));

	/* Save the current basis (actually the parent's basis)	*/
	/* into this node.					*/
	save_node_basis (p, bbip);

	/* Insert node into depth-first list... */
	p1 = tp -> first;
	if (p1 NE NULL) {
		p1 -> prev = p;
	}
	p -> next	= p1;
	p -> prev	= NULL;
	tp -> first = p;

	/* Insert node into "best-node" heap... */
	bbheap_insert (p, tp, BEST_NODE_HEAP);

	/* Insert node into "worst-node" heap... */
	bbheap_insert (p, tp, WORST_NODE_HEAP);
}

/*
 * Append the given branch-and-bound node to the given tree.  It is added
 * to the END of the node list, and sorted into each of the heaps.
 */

	void
append_node_to_tree (

struct bbnode *		p,	/* IN - node to append to tree */
struct bbtree *		tp	/* IN - tree to append it to */
)
{
struct bbnode *		p1;

	p1 = tp -> first;
	if (p1 EQ NULL) {
		tp -> first = p;
	}
	else {
		while (p1 -> next NE NULL) {
			p1 = p1 -> next;
		}
		p1 -> next = p;
	}
	p -> next = NULL;
	p -> prev = p1;

	bbheap_insert (p, tp, BEST_NODE_HEAP);
	bbheap_insert (p, tp, WORST_NODE_HEAP);
}

/*
 * This routine deletes an arbitrary node from an arbitrary position
 * within the given branch-and-bound tree.
 */

	void
delete_node_from_bbtree (

struct bbnode *		p,	/* IN - node to delete */
struct bbtree *		tp	/* IN - tree to delete it from */
)
{
struct bbnode * p1;
struct bbnode * p2;


	/* Delete it from LIFO list... */
	p1 = p -> prev;
	p2 = p -> next;
	if (p1 NE NULL) {
		p1 -> next = p2;
	}
	else {
		tp -> first = p2;
	}
	if (p2 NE NULL) {
		p2 -> prev = p1;
	}
	p -> next = NULL;
	p -> prev = NULL;

	/* Delete it from the best-node heap... */
	bbheap_delete (p, tp, BEST_NODE_HEAP);

	/* Delete it from the worst-node heap... */
	bbheap_delete (p, tp, WORST_NODE_HEAP);
}

/*
 * This routine creates an initial, empty branch-and-bound tree.
 */

	struct bbtree *
create_bbtree (

int		nmasks		/* IN - number of edge masks */
)
{
struct bbtree *		tp;

	tp = NEW (struct bbtree);

	tp -> first		= NULL;
	tp -> free		= NULL;
	tp -> snum		= 0;
	tp -> nmasks		= nmasks;
	tp -> node_policy	= NN_BEST_NODE;

	/* Initialize best and worst order heaps */
	bbheap_init (&(tp -> heap [BEST_NODE_HEAP]), node_is_better);
	bbheap_init (&(tp -> heap [WORST_NODE_HEAP]), node_is_worse);

	return (tp);
}

/*
 * This routine initializes a branch-and-bound node heap.
 */

	static
	void
bbheap_init (

struct bbheap *		hp,	/* IN - the heap to initialize */
bbheap_func_t *		funcp	/* IN - node comparison function */
)
{
	hp -> array		= NEWA (INIT_HSIZE, struct bbnode *);
	hp -> hsize		= INIT_HSIZE;
	hp -> nheap		= 0;
	hp -> is_parent_funcp	= funcp;
}



	static
	void
bbheap_free (

struct bbheap *		hp	/* IN - heap to free up */
)
{
int		i;

	for (i = 0; i < hp -> nheap; i++) {
		destroy_bbnode (hp -> array [i]);
	}

	free ((char *) (hp -> array));
}

/*
 * This routine adds the given branch-and-bound node to the specified
 * heap of the given branch-and-bound tree.
 */

	void
bbheap_insert (

struct bbnode *		p,	/* IN - node to insert */
struct bbtree *		tp,	/* IN - branch-and-bound tree with heaps */
int			heap_no	/* IN - heap number to insert node into */
)
{
int			i;
int			j;
int			n;
struct bbheap *		hp;
struct bbnode **	tmp;
struct bbnode *		p2;
bbheap_func_t *		is_parent;

	/* Use specified heap... */
	hp = &(tp -> heap [heap_no]);

	/* Verify sufficient heap space to hold new node... */
	n = hp -> nheap;
	if (n >= hp -> hsize) {
		tmp = NEWA (2 * n, struct bbnode *);
		for (i = 0; i < n; i++) {
			tmp [i] = hp -> array [i];
		}
		free ((char *) (hp -> array));
		hp -> array = tmp;
		hp -> hsize = 2 * n;
	}

	is_parent = hp -> is_parent_funcp;

	/* Insert node into heap by placing it at the end of	*/
	/* the array and sifting up...				*/
	for (i = n; i > 0;) {
		j = ((i - 1) >> 1);
		p2 = hp -> array [j];
		if (is_parent (p2, p)) break;
		p2 -> index [heap_no] = i;
		hp -> array [i] = p2;
		i = j;
	}
	p -> index [heap_no] = i;
	hp -> array [i] = p;
	hp -> nheap = n + 1;
}

/*
 * This routine deletes a single node from the specified heap of the
 * given branch-and-bound tree.
 */

	static
	void
bbheap_delete (

struct bbnode *		p,	/* IN - the node to delete from the heap */
struct bbtree *		tp,	/* IN - branch-and-bound tree with heaps */
int			heap_no	/* IN - heap number from which to delete */
)
{
int			i;
int			j;
int			n;
struct bbheap *		hp;
struct bbnode *		p1;
struct bbnode *		p2;
struct bbnode *		p3;
bbheap_func_t *		is_parent;

	/* Use proper heap to delete from... */
	hp = &(tp -> heap [heap_no]);
	n = hp -> nheap;

	if ((n <= 0) OR (n > hp -> hsize)) {
		fatal ("delete_bbheap: Bug 1.");
	}

	/* Deleting from a heap requires three steps:		*/
	/*	1. Move the last node into the deleted node's	*/
	/*	   current position.				*/
	/*	2. Sift this replacement node up.		*/
	/*	3. Sift this replacement node down.		*/

	/* Get position of node being deleted... */
	i = p -> index [heap_no];
	if (i >= n) {
		fatal ("delete_bbheap: Bug 2.");
	}
	if (hp -> array [i] NE p) {
		fatal ("delete_bbheap: Bug 3.");
	}

	/* Deleted node is no longer in the array... */
	p -> index [heap_no] = -1;

	/* Heap now has one fewer element in it... */
	--n;
	hp -> nheap = n;
	if (n <= 0) {
		/* Heap is now empty -- done! */
		return;
	}

	/* Move last heap element into vacated spot. */
	p1 = hp -> array [n];
	if (p1 EQ p) {
		/* Deleting last item from heap -- we are done! */
		return;
	}
	if (p1 -> index [heap_no] NE n) {
		fatal ("delete_bbheap: Bug 4.");
	}

	/* We now assume that node "p1" will be in position "i"... */

	is_parent = hp -> is_parent_funcp;

	/* First sift up... */
	while (i > 0) {
		j = ((i - 1) >> 1);
		p2 = hp -> array [j];
		if (is_parent (p2, p1)) break;
		p2 -> index [heap_no] = i;
		hp -> array [i] = p2;
		i = j;
	}

	/* Now sift down... */
	while (i < n) {
		j = (i << 1) + 1;
		if (j >= n) break;
		p2 = hp -> array [j];
		if ((j + 1) < n) {
			p3 = hp -> array [j + 1];
			if (is_parent (p3, p2)) {
				++j;
				p2 = p3;
			}
		}
		if (is_parent (p1, p2)) break;
		p2 -> index [heap_no] = i;
		hp -> array [i] = p2;
		i = j;
	}
	p1 -> index [heap_no] = i;
	hp -> array [i] = p1;

	hp -> nheap = n;
}

/*
 * This routine returns TRUE if-and-only-if node 1 is "better" than
 * node 2 -- in other words, node 1 must be above node 2 in the
 * "best-node" heap.
 */

	static
	int
node_is_better (

struct bbnode *		p1,	/* IN - node 1 */
struct bbnode *		p2	/* IN - node 2 */
)
{
	if (p1 -> z < p2 -> z) return (TRUE);
	if (p1 -> z > p2 -> z) return (FALSE);

	/* Objective values are equal -- use node creation	*/
	/* order to break the tie.  What we want to achieve is	*/
	/* that the "var=0" and "var=1" children of a node	*/
	/* (which will have equal objective values) get taken	*/
	/* from the "best-node" heap in proper order with	*/
	/* respect to the UP_FIRST switch.  (The node that was	*/
	/* created last should be taken from the heap first.)	*/

	if (p1 -> num >= p2 -> num) return (TRUE);
	return (FALSE);
}


/*
 * This routine returns TRUE if-and-only-if node 1 is "worse" than
 * node 2 -- in other words, node 2 must be above node 2 in the
 * "worst-node" heap.
 */

	static
	int
node_is_worse (

struct bbnode *		p1,	/* IN - node 1 */
struct bbnode *		p2	/* IN - node 2 */
)
{
	return (p1 -> z >= p2 -> z);
}

/*
 * This routine destroys the given branch-and-bound info structure,
 * freeing all of the memory it points to.
 */

	void
destroy_bbinfo (

struct bbinfo *		bbip		/* IN - branch-and-bound info */
)
{
struct bbtree *		bbtree;
struct bbnode *		p1;
struct bbnode *		p2;

	shutdown_heuristic_upper_bound (bbip -> ubip);

	if (bbip -> statp NE NULL) {
		free ((char *) (bbip -> statp));
	}
	bbip -> value = NULL;
	bbip -> fixed = NULL;
	free ((char *) (bbip -> dj));
	if (bbip -> slack NE NULL) {
		free ((char *) (bbip -> slack));
	}
	bbip -> node = NULL;
	free ((char *) (bbip -> smt));

	bbtree = bbip -> bbtree;
	if (bbtree NE NULL) {
		p1 = bbtree -> free;
		for (;;) {
			if (p1 EQ NULL) break;
			p2 = p1 -> next;
			destroy_bbnode (p1);
			p1 = p2;
		}
		bbtree -> first = NULL;
		bbtree -> free = NULL;

		bbheap_free (&(bbtree -> heap [BEST_NODE_HEAP]));
		bbheap_free (&(bbtree -> heap [WORST_NODE_HEAP]));
		free ((char *) bbtree);
	}

	if (bbip -> cpool NE NULL) {
		free_constraint_pool (bbip -> cpool);
	}

	if (bbip -> lp NE NULL) {
		destroy_initial_formulation (bbip);
	}

	if (bbip -> lpmem NE NULL) {
		free ((char *) (bbip -> lpmem));
	}

	/* This items all belong to the cinfo.  Just zap them. */
	bbip -> cip		= NULL;
	bbip -> tmap		= NULL;
	bbip -> fset_mask	= NULL;

	free ((char *) bbip);
}

/*
 * Destroy the given branch-and-bound node, freeing all of the memory
 * it refers to.
 */

	static
	void
destroy_bbnode (

struct bbnode *		p		/* IN - node to destroy */
)
{
	free ((char *) (p -> fixed));
	free ((char *) (p -> value));

	if (p -> x NE NULL) {
		free ((char *) (p -> x));
	}
	if (p -> zlb NE NULL) {
		free ((char *) (p -> zlb));
	}
	if (p -> bc_uids NE NULL) {
		free ((char *) (p -> bc_uids));
	}
	if (p -> bc_row NE NULL) {
		free ((char *) (p -> bc_row));
	}
	if (p -> rstat NE NULL) {
		free ((char *) (p -> rstat));
	}
	if (p -> cstat NE NULL) {
		free ((char *) (p -> cstat));
	}
	free ((char *) (p -> bheur));
	free ((char *) p);
}
