/***********************************************************************

	File:	rfst.c
	Rev:	b-1
	Date:	12/22/2000

	Copyright (c) 1998, 2001 by Martin Zachariasen & David M. Warme

************************************************************************

	The main routine for the rectilinear FST generator.  It
	reads a point set from standard input, generates the FSTs,
	and outputs them to standard output.

************************************************************************

	Modification Log:

	a-1:	08/31/98	warme
		: Created.  Derived from Martin's prototype.
	b-1:	12/22/2000	martinz
		: Added -k option to denote the maximum size of FSTs
		: that should be generated.
		: FSTs for which corner-flipped versions split
		: into two or more FSTs are removed.

************************************************************************/

#include "bsd.h"
#include "emptyr.h"
#include "p1io.h"
#include "rfst.h"
#include "steiner.h"


/*
 * Global Routines
 */

int			main (int, char **);


/*
 * Local Constants
 */

#define EMPTY_DIAMOND_PROPERTY	1
#define UB_SHORTLEG		1
#define NOSPLIT_CORNER_FLIPPED	1
#define KAHNG_ROBINS_HEURISTIC	0
#define DO_STATISTICS		0


/*
 * Local Types
 */

struct ibuf {
	struct ibuf *	next;
	int		count;
	int		buf [1];
};


/*
 * Local Routines
 */

static void		add_zero_length_fsts (struct rinfo *, int, int **);
static void		build_fst_list (struct rinfo *);
static void		build_rfst_graph (struct rinfo *,
					  int,
					  int,
					  int,
					  struct pset *,
					  struct pset *,
					  struct edge *,
					  int);
static void		compute_rfsts_for_unique_terminals (struct rinfo *);
static void		compute_successors (struct rinfo *);
static void		compute_ub0 (struct rinfo *);
static void		compute_ub1 (struct rinfo *);
static void		compute_zt (struct rinfo *);
static void		convert_delta_cpu_time (char *);
static void		decode_params (int, char **);
static dist_t		delta_x_func (struct point *, struct point *);
static dist_t		delta_y_func (struct point *, struct point *);
static bool		diamond_empty (struct rinfo *,
				       struct point *,
				       struct point *,
				       int,
				       int);
static int		generate_duplicate_terminal_groups (struct pset *,
							    int *,
							    int ***);
static cpu_time_t	get_delta_cpu_time (void);
static void		grow_RFST (struct rinfo *	rip,
				   int			size,
				   dist_t		length,
				   int			dir,
				   dist_t		ub_length,
				   dist_t *		ub_shortleg,
				   int			longindex);
static int *		heapsort_x (struct pset *);
static int *		heapsort_y (struct pset *);
static int		lrindex_dir_0 (struct point *, struct point *);
static int		lrindex_dir_1 (struct point *, struct point *);
static int		lrindex_dir_2 (struct point *, struct point *);
static int		lrindex_dir_3 (struct point *, struct point *);
static struct full_set ** put_trees_in_array (struct full_set *, int *);
static struct pset *	remove_duplicates (struct pset * pts,
					   int		 ndg,
					   int **	 list,
					   int **	 fwd_map,
					   int **	 rev_map);
static void		renumber_terminals (struct rinfo *,
					    struct pset *,
					    int *);
static dist_t		test_and_save_fst (struct rinfo *,
					   int,
					   dist_t,
					   int,
					   int);
static void		usage (void);


/*
 * Local Variables
 */

static char *		description;
static char *		me;
static int		MaxFSTSize = 0;
static int		output_version = CURRENT_P1IO_VERSION;
static bool		Print_Detailed_Timings = FALSE;
static cpu_time_t	T0;
static cpu_time_t	Tn;

/*
 * Data tables for implementing the DSTDIR and DSTDIRP macros without
 * conditional branches.  This should be faster on machines for which
 * pipeline flushes are an issue.
 * We could really use the pointer-to-member feature of C++ here!
 */

#define X	offsetof(struct point, x)
#define Y	offsetof(struct point, y)

const int8u	dstdir_offsets [] = {
	X, Y, X, Y, X,
};

#undef X
#undef Y

/* Tables for the array-of-function-pointers implementation. */

typedef dist_t	(* dstdir_funcp) (struct point *, struct point *);

const dstdir_funcp	dstdir_funcs [] = {
	delta_x_func,
	delta_y_func,
	delta_x_func,
	delta_y_func,
	delta_x_func,
};


/* Tables for the LRINDEX implementation using an array of functions. */

const lrindex_funcp	lrindex_funcs [] = {
	lrindex_dir_0,
	lrindex_dir_1,
	lrindex_dir_2,
	lrindex_dir_3,
};

/*
 * The main routine for the "rfst" program.  It reads a point set
 * from standard input, generates all of the rectilinear FSTs,
 * and outputs them to standard output in our special "phase 1 I/O"
 * format.
 */

	int
main (

int		argc,
char **		argv
)
{
int			i;
int			j;
int			k;
int			ndg;
int			ntrees;
int			count;
int			fpsave;
int **			dup_grps;
int *			fwd_map;
int *			rev_map;
int *			ip1;
struct full_set *	fsp;
struct pset *		terms;
struct rinfo		rinfo;
struct cinfo		cinfo;
struct pset *		pts;
struct pset *		pts2;
cpu_time_t		Trenum;
struct scale_info	scale;
char			buf1 [32];

	fpsave = set_floating_point_double_precision ();

	setbuf (stdout, NULL);

	decode_params (argc, argv);

	pts = get_points (stdin, &scale);

	init_output_conversion (pts, RECTILINEAR, &scale);

	T0 = get_cpu_time ();
	Tn = T0;

	rinfo.x_order = heapsort_x (pts);

	if (Print_Detailed_Timings) {
		convert_delta_cpu_time (buf1);
		fprintf (stderr, "Sort X:                 %s\n", buf1);
	}

	/* Find all duplicate terminals in the input. */
	ndg = generate_duplicate_terminal_groups (pts,
						  rinfo.x_order,
						  &dup_grps);

	rinfo.num_term_masks = BMAP_ELTS (pts -> n);

	/* Remove all but the first of each duplicate terminal. */
	/* Compute forward and reverse maps to renumber the terminals. */
	pts2 = remove_duplicates (pts, ndg, dup_grps, &fwd_map, &rev_map);

	/* Renumber the x_order list -- instead of re-sorting pts2. */
	j = 0;
	for (i = 0; i < pts -> n; i++) {
		k = rinfo.x_order [i];
		if ((k < 0) OR (pts -> n < k)) {
			fatal ("main: Bug 1.");
		}
		k = fwd_map [k];
		if (k < 0) continue;
		rinfo.x_order [j++] = k;
	}

	if (Print_Detailed_Timings) {
		convert_delta_cpu_time (buf1);
		fprintf (stderr, "Remove Duplicates:      %s\n", buf1);
	}

	rinfo.y_order = heapsort_y (pts2);

	if (Print_Detailed_Timings) {
		convert_delta_cpu_time (buf1);
		fprintf (stderr, "Sort Y:                 %s\n", buf1);
	}

	/* From now on, we work only with the reduced terminal set, and */
	/* we assume that all terminals are unique. */

	rinfo.pts = pts2;

	compute_rfsts_for_unique_terminals (&rinfo);

	/* Now put the terminal numbers back the way they were, */
	/* renumber the terminals within each RFST, etc. */

	renumber_terminals (&rinfo, pts, rev_map);

	/* Link the FSTs together into one long list, and number them. */
	build_fst_list (&rinfo);

	/* Add one FST for each duplicate terminal that was removed. */
	if (ndg > 0) {
		add_zero_length_fsts (&rinfo, ndg, dup_grps);
	}

	/* Measure renumber time.  This also sets Tn so that Tn-T0 is	*/
	/* the total processing time.					*/
	Trenum = get_delta_cpu_time ();

	if (Print_Detailed_Timings) {
		convert_cpu_time (Trenum, buf1);
		fprintf (stderr, "Renumber Terminals:     %s\n", buf1);
		convert_cpu_time (Tn - T0, buf1);
		fprintf (stderr, "Total:                  %s\n", buf1);
	}

	/* Zero out the compatibility info. */
	memset (&cinfo, 0, sizeof (cinfo));

	cinfo.num_edges			= rinfo.ntrees;
	cinfo.num_verts			= rinfo.pts -> n;
	cinfo.num_vert_masks		= BMAP_ELTS (cinfo.num_verts);
	cinfo.num_edge_masks		= BMAP_ELTS (cinfo.num_edges);
	cinfo.edge			= NEWA (rinfo.ntrees + 1, int *);
	cinfo.edge_size			= NEWA (rinfo.ntrees, int);
	cinfo.cost			= NEWA (rinfo.ntrees, dist_t);
	cinfo.tflag			= NEWA (cinfo.num_verts, bool);
	cinfo.metric			= RECTILINEAR;
	cinfo.scale			= scale;
	cinfo.integrality_delta		= 1.0;
	cinfo.pts			= rinfo.pts;
	cinfo.full_trees		= put_trees_in_array (rinfo.full_sets,
							      &ntrees);
	cinfo.mst_length		= rinfo.mst_length;
	cinfo.description		= gst_strdup (description);
	cinfo.p1time			= Tn - T0;

	for (i = 0; i < cinfo.num_verts; i++) {
		cinfo.tflag [i] = TRUE;
	}

	count = 0;
	for (i = 0; i < rinfo.ntrees; i++) {
		fsp = cinfo.full_trees [i];
		k = fsp -> terminals -> n;
		cinfo.edge_size [i]	= k;
		cinfo.cost [i]		= fsp -> tree_len;
		count += k;
	}
	ip1 = NEWA (count, int);
	for (i = 0; i < rinfo.ntrees; i++) {
		cinfo.edge [i] = ip1;
		terms = cinfo.full_trees [i] -> terminals;
		k = terms -> n;
		for (j = 0; j < k; j++) {
			*ip1++ = terms -> a [j].pnum;
		}
	}
	cinfo.edge [i] = ip1;

	/* Print all of the data we have generated... */
	print_phase_1_data (&cinfo, output_version);

	/* Clean up. */

	free ((char *) (rinfo.y_order));
	free ((char *) fwd_map);
	free ((char *) rev_map);
	free ((char *) (rinfo.x_order));

	free ((char *) pts2);
#if 0
	free ((char *) pts);
#endif

	free_phase_1_data (&cinfo);

	restore_floating_point_precision (fpsave);

	exit (0);
}

/*
 * This routine decodes the various command-line arguments.
 */

	static
	void
decode_params (

int		argc,
char **		argv
)
{
char *		ap;
char		c;
int		v;

	--argc;
	me = *argv++;
	while (argc > 0) {
		ap = *argv++;
		if (*ap NE '-') {
			usage ();
		}
		++ap;
		while ((c = *ap++) NE '\0') {
			switch (c) {
			case 'd':
				if (*ap EQ '\0') {
					if (argc <= 0) {
						usage ();
					}
					ap = *argv++;
					--argc;
				}
				if (strlen (ap) >= 80) {
					fprintf (stderr,
						 "Description must be less"
						 " than 80 characters.\n");
					usage ();
				}
				description = ap;
				/* Change newlines to spaces... */
				for (;;) {
					ap = strchr (ap, '\n');
					if (ap EQ NULL) break;
					*ap++ = ' ';
				}
				ap = "";
				break;

			case 'k':
				/* Max number of terminals in FSTs */
				if (*ap EQ '\0') {
					if (argc <= 0) {
						usage ();
					}
					ap = *argv++;
					--argc;
				}
				MaxFSTSize = atoi (ap);
				ap = "";
				break;

			case 't':
				Print_Detailed_Timings = TRUE;
				break;

			case 'v':
				if (*ap EQ '\0') {
					if (argc <= 0) {
						usage ();
					}
					ap = *argv++;
					--argc;
				}
				if ((*ap < '0') OR (*ap > '9')) {
					usage ();
				}
				v = atoi (ap);
				if ((v < 0) OR (v > CURRENT_P1IO_VERSION)) {
					fprintf (stderr,
						 "%s: Bad version `%s'."
						 "  Valid versions range"
						 " from 0 to %d.\n",
						 me,
						 ap,
						 CURRENT_P1IO_VERSION);
					usage ();
				}
				output_version = v;
				ap = "";
				break;

			default:
				usage ();
				break;
			}
		}
		--argc;
	}
}

/*
 * This routine prints out the proper usage and exits.
 */

static char *	arg_doc [] = {
	"",
	"\t-d txt\tDescription of problem instance.",
	"\t-k K\tOnly generate FSTs spanning up to K terminals.",
	"\t-t\tPrint detailed timings on stderr.",
	"\t-v N\tGenerates version N output data format.",
	"",
	NULL
};

	static
	void
usage (void)

{
char **		pp;
char *		p;

	(void) fprintf (stderr,
			"\nUsage: %s [-t] [-d description] [-k K] [-v N]"
			" <points-file\n",
			me);
	pp = &arg_doc [0];
	while ((p = *pp++) NE NULL) {
		(void) fprintf (stderr, "%s\n", p);
	}
	exit (1);
}

/*
 * Use the heapsort algorithm to sort the given terminals in increasing
 * order by the following keys:
 *
 *	1.	X coordinate
 *	2.	Y coordinate
 *	3.	index (i.e., position within input data)
 *
 * Of course, we do not move the points, but rather permute an array
 * of indexes into the points.
 */

	static
	int *
heapsort_x (

struct pset *		pts		/* IN - the terminals to sort */
)
{
int			i, i1, i2, j, k, n;
struct point *		p1;
struct point *		p2;
int *			index;

	n = pts -> n;

	index = NEWA (n, int);
	for (i = 0; i < n; i++) {
		index [i] = i;
	}

	/* Construct the heap via sift-downs, in O(n) time. */
	for (k = n >> 1; k >= 0; k--) {
		j = k;
		for (;;) {
			i = (j << 1) + 1;
			if (i + 1 < n) {
				/* Increment i (to right subchild of j) */
				/* if the right subchild is greater. */
				i1 = index [i];
				i2 = index [i + 1];
				p1 = &(pts -> a [i1]);
				p2 = &(pts -> a [i2]);
				if ((p2 -> x > p1 -> x) OR
				    ((p2 -> x EQ p1 -> x) AND
				     ((p2 -> y > p1 -> y) OR
				      ((p2 -> y EQ p1 -> y) AND
				       (i2 > i1))))) {
					++i;
				}
			}
			if (i >= n) {
				/* Hit bottom of heap, sift-down is done. */
				break;
			}
			i1 = index [j];
			i2 = index [i];
			p1 = &(pts -> a [i1]);
			p2 = &(pts -> a [i2]);
			if ((p1 -> x > p2 -> x) OR
			    ((p1 -> x EQ p2 -> x) AND
			     ((p1 -> y > p2 -> y) OR
			      ((p1 -> y EQ p2 -> y) AND
			       (i1 > i2))))) {
				/* Greatest child is smaller.  Sift-	*/
				/* down is done. */
				break;
			}
			/* Sift down and continue. */
			index [j] = i2;
			index [i] = i1;
			j = i;
		}
	}

	/* Now do actual sorting.  Exchange first/last and sift down. */
	while (n > 1) {
		/* Largest is at index [0], swap with index [n-1],	*/
		/* thereby putting it into final position.		*/
		--n;
		i = index [0];
		index [0] = index [n];
		index [n] = i;

		/* Now restore the heap by sifting index [0] down. */
		j = 0;
		for (;;) {
			i = (j << 1) + 1;
			if (i + 1 < n) {
				/* Increment i (to right subchild of j) */
				/* if the right subchild is greater. */
				i1 = index [i];
				i2 = index [i + 1];
				p1 = &(pts -> a [i1]);
				p2 = &(pts -> a [i2]);
				if ((p2 -> x > p1 -> x) OR
				    ((p2 -> x EQ p1 -> x) AND
				     ((p2 -> y > p1 -> y) OR
				      ((p2 -> y EQ p1 -> y) AND
				       (i2 > i1))))) {
					++i;
				}
			}
			if (i >= n) {
				/* Hit bottom of heap, sift-down is done. */
				break;
			}
			i1 = index [j];
			i2 = index [i];
			p1 = &(pts -> a [i1]);
			p2 = &(pts -> a [i2]);
			if ((p1 -> x > p2 -> x) OR
			    ((p1 -> x EQ p2 -> x) AND
			     ((p1 -> y > p2 -> y) OR
			      ((p1 -> y EQ p2 -> y) AND
			       (i1 > i2))))) {
				/* Greatest child is smaller.  Sift-	*/
				/* down is done. */
				break;
			}
			/* Sift down and continue. */
			index [j] = i2;
			index [i] = i1;
			j = i;
		}
	}

	return (index);
}

/*
 * Use the heapsort algorithm to sort the given terminals in increasing
 * order by the following keys:
 *
 *	1.	Y coordinate
 *	2.	X coordinate
 *	3.	index (i.e., position within input data)
 *
 * Of course, we do not move the points, but rather permute an array
 * of indexes into the points.
 */

	static
	int *
heapsort_y (

struct pset *		pts		/* IN - the terminals to sort */
)
{
int			i, i1, i2, j, k, n;
struct point *		p1;
struct point *		p2;
int *			index;

	n = pts -> n;

	index = NEWA (n, int);
	for (i = 0; i < n; i++) {
		index [i] = i;
	}

	/* Construct the heap via sift-downs, in O(n) time. */
	for (k = n >> 1; k >= 0; k--) {
		j = k;
		for (;;) {
			i = (j << 1) + 1;
			if (i + 1 < n) {
				/* Increment i (to right subchild of j) */
				/* if the right subchild is greater. */
				i1 = index [i];
				i2 = index [i + 1];
				p1 = &(pts -> a [i1]);
				p2 = &(pts -> a [i2]);
				if ((p2 -> y > p1 -> y) OR
				    ((p2 -> y EQ p1 -> y) AND
				     ((p2 -> x > p1 -> x) OR
				      ((p2 -> x EQ p1 -> x) AND
				       (i2 > i1))))) {
					++i;
				}
			}
			if (i >= n) {
				/* Hit bottom of heap, sift-down is done. */
				break;
			}
			i1 = index [j];
			i2 = index [i];
			p1 = &(pts -> a [i1]);
			p2 = &(pts -> a [i2]);
			if ((p1 -> y > p2 -> y) OR
			    ((p1 -> y EQ p2 -> y) AND
			     ((p1 -> x > p2 -> x) OR
			      ((p1 -> x EQ p2 -> x) AND
			       (i1 > i2))))) {
				/* Greatest child is smaller.  Sift-	*/
				/* down is done. */
				break;
			}
			/* Sift down and continue. */
			index [j] = i2;
			index [i] = i1;
			j = i;
		}
	}

	/* Now do actual sorting.  Exchange first/last and sift down. */
	while (n > 1) {
		/* Largest is at index [0], swap with index [n-1],	*/
		/* thereby putting it into final position.		*/
		--n;
		i = index [0];
		index [0] = index [n];
		index [n] = i;

		/* Now restore the heap by sifting index [0] down. */
		j = 0;
		for (;;) {
			i = (j << 1) + 1;
			if (i + 1 < n) {
				/* Increment i (to right subchild of j) */
				/* if the right subchild is greater. */
				i1 = index [i];
				i2 = index [i + 1];
				p1 = &(pts -> a [i1]);
				p2 = &(pts -> a [i2]);
				if ((p2 -> y > p1 -> y) OR
				    ((p2 -> y EQ p1 -> y) AND
				     ((p2 -> x > p1 -> x) OR
				      ((p2 -> x EQ p1 -> x) AND
				       (i2 > i1))))) {
					++i;
				}
			}
			if (i >= n) {
				/* Hit bottom of heap, sift-down is done. */
				break;
			}
			i1 = index [j];
			i2 = index [i];
			p1 = &(pts -> a [i1]);
			p2 = &(pts -> a [i2]);
			if ((p1 -> y > p2 -> y) OR
			    ((p1 -> y EQ p2 -> y) AND
			     ((p1 -> x > p2 -> x) OR
			      ((p1 -> x EQ p2 -> x) AND
			       (i1 > i2))))) {
				/* Greatest child is smaller.  Sift-	*/
				/* down is done. */
				break;
			}
			/* Sift down and continue. */
			index [j] = i2;
			index [i] = i1;
			j = i;
		}
	}

	return (index);
}

/*
 * This routine finds all pairs of points whose coordinates are exactly
 * identical.  This is a degeneracy that can cause successive algorithms
 * much heartburn, so we take care of them immediately by keeping only
 * the earliest (lowest pnum) point, and marking the others as duplicates.
 * We generate a partition of such terminals into subsets -- if terminals
 * A and B reside in the same subset then they have identical coordinates.
 * We refer to such a subset as a Duplicate Terminal Group.
 *
 * For each terminal group we retain ONLY THE FIRST member -- the terminal
 * map bits are turned OFF for all other members of a terminal group,
 * effectively eliminating those terminals from the problem.
 */

	static
	int
generate_duplicate_terminal_groups (

struct pset *		pts,	/* IN - original terminal set */
int *			xorder, /* IN - terminal numbers sorted by X coord */
int ***			grps	/* OUT - duplicate terminal groups */
)
{
int			i;
int			j;
int			n;
int			n_grps;
struct point *		p0;
struct point *		p1;
struct point *		p2;
int *			ip;
int **			real_ptrs;
int *			real_terms;
int **			ptrs;
int *			terms;

	n = pts -> n;

	n_grps = 0;

	ptrs	= NEWA (n + 1, int *);
	terms	= NEWA (n, int);

	ip = &terms [0];
	for (i = 1; i < n; ) {
		p0 = &(pts -> a [xorder [i - 1]]);
		p1 = &(pts -> a [xorder [i]]);

		if ((p0 -> y NE p1 -> y) OR (p0 -> x NE p1 -> x)) {
			/* Not identical. */
			++i;
			continue;
		}

		/* Terminals xorder[i-1] and xorder[i] are identical. */

		for (j = i + 1; j < n; j++) {
			p2 = &(pts -> a [xorder [j]]);
			if (p0 -> y NE p2 -> y) break;
			if (p0 -> x NE p2 -> x) break;
		}
		/* j now points to the first non-equal terminal */

		/* Make a new duplicate terminal group... */
		ptrs [n_grps++] = ip;
		*ip++ = xorder [i - 1];
		while (i < j) {
			*ip++ = xorder [i++];
		}

		/* Skip whole group of coincident points and continue */
	}
	ptrs [n_grps] = ip;

	if (n_grps <= 0) {
		*grps = NULL;
	}
	else {
		/* Transfer to permanent memory of proper size... */
		real_ptrs = NEWA (n_grps + 1, int *);
		real_terms = NEWA (ip - terms, int);
		(void) memcpy ((char *) real_terms,
			       (char *) terms,
			       (ip - terms) * sizeof (int));
		for (i = 0; i <= n_grps; i++) {
			real_ptrs [i] = &real_terms [ptrs [i] - ptrs [0]];
		}
		*grps = real_ptrs;
	}

	free ((char *) terms);
	free ((char *) ptrs);

	return (n_grps);
}

/*
 * This routine removes all but the first terminal in each duplicate
 * terminal group.  We also prepare arrays for mapping forward from
 * old to new terminal numbers, and backward from new terminal numbers
 * to old.
 */

	static
	struct pset *
remove_duplicates (

struct pset *		pts,		/* IN - original point set */
int			ndg,		/* IN - number of duplicate groups */
int **			list,		/* IN - list of duplicate groups */
int **			fwd_map_ptr,	/* OUT - map to renumber old to new */
int **			rev_map_ptr	/* OUT - map to renumber new to old */
)
{
int		i;
int		j;
int		n;
int		kmasks;
int		numdel;
int		new_n;
int *		ip1;
int *		ip2;
int *		fwd;
int *		rev;
bitmap_t *	deleted;
struct pset *	newpts;

	n = pts -> n;

	kmasks = BMAP_ELTS (n);
	deleted = NEWA (kmasks, bitmap_t);
	for (i = 0; i < kmasks; i++) {
		deleted [i] = 0;
	}

	numdel = 0;
	for (i = 0; i < ndg; i++) {
		ip1 = list [i];
		ip2 = list [i + 1];

		/* Retain the first in this group, exclude all	*/
		/* of the others...				*/
		while (++ip1 < ip2) {
			if (BITON (deleted, *ip1)) {
				fatal ("remove_duplicates: Bug 1.");
			}
			++numdel;
			SETBIT (deleted, *ip1);
		}
	}

	new_n = n - numdel;

	fwd = NEWA (n, int);
	rev = NEWA (new_n, int);
	newpts = NEW_PSET (new_n);
	memset (newpts, 0, PSET_SIZE (new_n));
	newpts -> n = new_n;
	j = 0;
	for (i = 0; i < n; i++) {
		if (BITON (deleted, i)) {
			fwd [i] = -1;
		}
		else {
			newpts -> a [j].x		= pts -> a [i].x;
			newpts -> a [j].y		= pts -> a [i].y;
			newpts -> a [j].pnum		= j;
			rev [j] = i;
			fwd [i] = j;
			++j;
		}
	}

	free ((char *) deleted);

	*fwd_map_ptr = fwd;
	*rev_map_ptr = rev;

	return (newpts);
}

/*
 * Compute the RFSTs for the given set of terminals, which are now
 * guaranteed to be unique.
 */

	static
	void
compute_rfsts_for_unique_terminals (

struct rinfo *		rip	/* IN - global RFST info */
)
{
int			i;
int			n;
int			dir;
int			nedges;
struct pset *		pts;
struct edge *		ep;
struct edge *		mst_edges;
struct rlist *		rp;
dist_t			mst_len;
double			ub_shortleg [2];
char			buf1 [32];

	pts = rip -> pts;
	n = pts -> n;

	compute_successors (rip);

	if (Print_Detailed_Timings) {
		convert_delta_cpu_time (buf1);
		fprintf (stderr, "Compute Successors:     %s\n", buf1);
	}

	rip -> empty_rect = init_empty_rectangles (pts, rip -> succ [0]);

	if (Print_Detailed_Timings) {
		convert_delta_cpu_time (buf1);
		fprintf (stderr, "Empty Rectangles:       %s\n", buf1);
	}

	compute_ub0 (rip);

	if (Print_Detailed_Timings) {
		convert_delta_cpu_time (buf1);
		fprintf (stderr, "Compute UB0:            %s\n", buf1);
	}

	mst_edges = NEWA (n - 1, struct edge);
	nedges = rect_mst (pts, mst_edges, rip -> empty_rect);
	if (nedges NE n - 1) {
		fatal ("compute_rfsts_for_unique_terminals: Bug 1.");
	}

	mst_len = 0.0;
	ep = mst_edges;
	for (i = 0; i < nedges; ep++, i++) {
		mst_len += ep -> len;
	}
	rip -> mst_length = mst_len;

	if (Print_Detailed_Timings) {
		convert_delta_cpu_time (buf1);
		fprintf (stderr, "Compute MST:            %s\n", buf1);
	}

	rip -> bsd = compute_bsd (nedges, mst_edges, 0);

	free ((char *) mst_edges);

	if (Print_Detailed_Timings) {
		convert_delta_cpu_time (buf1);
		fprintf (stderr, "Compute BSD:            %s\n", buf1);
	}

	compute_zt (rip);

	if (Print_Detailed_Timings) {
		convert_delta_cpu_time (buf1);
		fprintf (stderr, "Compute ZT:             %s\n", buf1);
	}

	compute_ub1 (rip);

	if (Print_Detailed_Timings) {
		convert_delta_cpu_time (buf1);
		fprintf (stderr, "Compute UB1:            %s\n", buf1);
	}

	/*
	 * Preprocessing finished.  Start growing FSTs!
	 */

	rip -> terms		= NEWA (n, int);
	rip -> longterms	= NEWA (n + 1, int);
	rip -> maxedges		= NEWA (n, dist_t);
	rip -> shortterm	= NEWA (n, int);
	rip -> lrindex		= NEWA (n, int);
	rip -> term_check	= NEWA (n, bool);
	rip -> hash		= NEWA (n, struct rlist *);

	for (i = 0; i < n; i++) {
		rip -> lrindex [i] = 0;
		rip -> term_check [i] = FALSE;
		rip -> hash [i] = NULL;
	}
	rip -> list.forw	= &(rip -> list);
	rip -> list.back	= &(rip -> list);

	rip -> fsts_checked = 0;

	ub_shortleg [0] = INF_DISTANCE;
	ub_shortleg [1] = INF_DISTANCE;
	if (MaxFSTSize EQ 0) MaxFSTSize = n;

	for (dir = 0; dir < 2; dir++) {
		for (i = 0; i < n; i++) {
			rip -> terms [0]	= i;
			rip -> maxedges [i]	= 0.0;
			/* Long leg candidate list is initially empty,	*/
			/* add candidates on demand.			*/
			rip -> longterms [0]	= i;
			rip -> longterms [1]	= -1;
			grow_RFST (rip,
				   1,		/* size */
				   0.0,		/* length */
				   dir,
				   0.0,		/* ub_length */
				   ub_shortleg,
				   0);		/* longindex */

#if DO_STATISTICS
			for (j = 0; longterms [j] >= 0; j++) {
				++long_leg_count;
			}
#endif
		}
	}

	/* Finally add MST-edges */

	ep = &(rip -> bsd -> mst_edges [1]);
	for (i = 1; i < n; i++) {
		rip -> terms [0] = ep -> p1;
		rip -> terms [1] = ep -> p2;
		test_and_save_fst (rip,
				   2,
				   ep -> len,
				   0,
				   TYPE_1);
		++ep;
	}

#if DO_STATISTICS
	fprintf (stderr, "Total FSTs checked: %d\n", rip -> fsts_checked);
#endif

	if (Print_Detailed_Timings) {
		convert_delta_cpu_time (buf1);
		fprintf (stderr, "Grow FSTs:              %s\n", buf1);
	}

	/* Disconnect FSTs from hash table. */

	for (rp = rip -> list.forw;
	     rp NE &(rip -> list);
	     rp = rp -> forw) {
		rp -> next = NULL;
	}

	/* Clean up. */

	free ((char *) (rip -> hash));		rip -> hash = NULL;
	free ((char *) (rip -> term_check));	rip -> term_check = NULL;
	free ((char *) (rip -> lrindex));	rip -> lrindex = NULL;
	free ((char *) (rip -> shortterm));	rip -> shortterm = NULL;
	free ((char *) (rip -> maxedges));	rip -> maxedges = NULL;
	free ((char *) (rip -> longterms));	rip -> longterms = NULL;
	free ((char *) (rip -> terms));		rip -> terms = NULL;

	free ((char *) (rip -> ub1 [0]));
	free ((char *) (rip -> zt [0] [0]));
	free ((char *) (rip -> zt [0]));

	shutdown_bsd (rip -> bsd);

	free ((char *) (rip -> ub0 [0]));
	free ((char *) (rip -> empty_rect));	rip -> empty_rect = NULL;
	free ((char *) (rip -> succ [0]));

	memset (&(rip -> succ), 0, sizeof (rip -> succ));
	memset (&(rip -> ub0),	0, sizeof (rip -> ub0));
	memset (&(rip -> ub1),	0, sizeof (rip -> ub1));
	memset (&(rip -> zt),	0, sizeof (rip -> zt));
}

/*
 * Sort the terminals by both X and Y coordinates, and then create the
 * successor lists.  These permit us to start from a random terminal and
 * traverse all remaining terminals in order by one of the four compass
 * directions.
 */

	static
	void
compute_successors (

struct rinfo *		rip		/* The global RFST info */
)
{
int			i;
int			n;
struct pset *		pts;
int *			index;
int *			succ;

	pts = rip -> pts;
	n = pts -> n;

	/* Allocate permanent successor lists. */
	succ = NEWA (4 * n, int);
	rip -> succ [0] = succ;		succ += n;
	rip -> succ [1] = succ;		succ += n;
	rip -> succ [2] = succ;		succ += n;
	rip -> succ [3] = succ;
	rip -> succ [4] = rip -> succ [0];
	rip -> succ [5] = rip -> succ [1];
	rip -> succ [6] = rip -> succ [2];

	/* Get terminal numbers sorted by X coordinate. */
	index = rip -> x_order;

	/* Set direction 0 (due east) successors. */
	succ = rip -> succ [0];
	for (i = 1; i < n; i++) {
		succ [index [i - 1]] = index [i];
	}
	succ [index [i - 1]] = -1;

	/* Set direction 2 (due west) successors. */
	succ = rip -> succ [2];
	for (i = 1; i < n; i++) {
		succ [index [i]] = index [i - 1];
	}
	succ [index [0]] = -1;

	/* Get terminal numbers sorted by Y coordinate. */
	index = rip -> y_order;

	/* Set direction 1 (due north) successors. */
	succ = rip -> succ [1];
	for (i = 1; i < n; i++) {
		succ [index [i - 1]] = index [i];
	}
	succ [index [i - 1]] = -1;

	/* Set direction 3 (due south) successors. */
	succ = rip -> succ [3];
	for (i = 1; i < n; i++) {
		succ [index [i]] = index [i - 1];
	}
	succ [index [0]] = -1;
}

/*
 * Compute the UB0 bounds.
 */

	static
	void
compute_ub0 (

struct rinfo *		rip		/* IN/OUT - global RFST info */
)
{
int			i;
int			j;
int			n;
int			dir;
struct pset *		pts;
struct point *		p1;
struct point *		p2;
int *			succ;
dist_t *		array;
dstdiroff_t		dir_off;
dstdiroff_t		dirp_off;
double			bound;
double			d1, d2, d3;

	pts = rip -> pts;
	n = pts -> n;

	array = NEWA (4 * n, dist_t);
	rip -> ub0 [0]	= array;	array += n;
	rip -> ub0 [1]	= array;	array += n;
	rip -> ub0 [2]	= array;	array += n;
	rip -> ub0 [3]	= array;
	rip -> ub0 [4]	= rip -> ub0 [0];
	rip -> ub0 [5]	= rip -> ub0 [1];
	rip -> ub0 [6]	= rip -> ub0 [2];

	for (dir = 0; dir < 4; dir++) {
		array =	 rip -> ub0 [dir];
		succ = rip -> succ [dir];
		dir_off = DSTDIR_GETOFFSET (dir);
		dirp_off = DSTDIR_GETOFFSET (dir + 1);

		p1 = &(pts -> a [0]);
		for (i = 0; i < n; i++, p1++) {
			bound = INF_DISTANCE;
			for (j = succ [i]; j >= 0; j = succ [j]) {
				p2 = &(pts -> a [j]);
				d1 = DSTDIR_OFFSET (p1, p2, dir_off);
				if (d1 > bound) break;
				d2 = DSTDIR_OFFSET (p1, p2, dirp_off);
				if (d1 > d2) {
					d3 = d1 + d2;	/* RDIST (p1,p2) */
					if (d3 < bound) {
						bound = d3;
					}
				}
			}
			array [i] = bound;
		}
	}
}

/*
 * Find the short leg candidates for each terminal i and direction.
 */

	static
	void
compute_zt (

struct rinfo *		rip		/* IN/OUT - global RFST info */
)
{
int			i;
int			j;
int			n;
int			dir;
int			dirp;
int			lr;
int			index;
struct pset *		pts;
struct point *		p1;
struct point *		p2;
int *			succ;
int *			ip;
int *			ip1;
int **			zt;
dstdiroff_t		dir_off;
dstdiroff_t		dir1_off;
lrdiroff_t		lrdir_off;
dist_t *		ub0;
dist_t *		ub0p;
dist_t			d1, d2;
dist_t			limit;
dist_t			b;
struct ibuf *		root;
struct ibuf *		ibp;
struct ibuf **		hookp;
int *			bufp;
int			bleft;
int *			offset;
int *			op;

static int		lr_table [] = {0, 1, 1, 0};
static int		dirp_table [] = {3, 2, 3, 2};


	pts = rip -> pts;
	n = pts -> n;

	/* Compute the short leg candidate lists.  We save them in a	*/
	/* dynamic buffer first.  Later we massage the data structure	*/
	/* into something more appropriate for access.			*/

	root = NULL;
	hookp = &root;
	bufp = NULL;
	bleft = 0;
	index = 0;
	offset = NEWA (4 * (n + 1), int);
	op = offset;

	for (dir = 0; dir < 4; dir++) {
		lr	= lr_table [dir];
		dirp	= dirp_table [dir];

		succ = rip -> succ [dir];
		dir_off = DSTDIR_GETOFFSET (dir);
		dir1_off = DSTDIR_GETOFFSET (dir + 1);
		lrdir_off = LRINDEX_GETOFFSET (dir);
		ub0 = rip -> ub0 [dir];
		ub0p = rip -> ub0 [dirp];

		p1 = &(pts -> a [0]);
		for (i = 0; i < n; i++, p1++) {
			*op++ = index;
			limit = ub0 [i];
			for (j = succ [i]; j >= 0; j = succ [j]) {
				p2 = &(pts -> a [j]);
				d1 = DSTDIR_OFFSET (p1, p2, dir_off);
				if (d1 == 0.0) continue;
				if (d1 > limit) break;

				d2 = DSTDIR_OFFSET (p1, p2, dir1_off);
				if (d2 EQ 0.0) break;
				if (LRINDEX_OFFSET (p1, p2, lrdir_off) NE lr) {
					continue;
				}
				if (d2 > ub0p [j]) continue;

				b = bsd (rip -> bsd, i, j);
				if (d1 > b) continue;
				if (d2 > b) continue;
				if (is_empty_rectangle (rip -> empty_rect, i, j)) {
					/* j is a short leg candidate for i */
					/* Store it away! */
					if (bleft <= 0) {
						ibp = (struct ibuf *)
							new (offsetof (struct ibuf,
								       buf [n]));
						ibp -> next = NULL;
						ibp -> count = 0;
						*hookp = ibp;
						hookp = &(ibp -> next);
						bufp = &(ibp -> buf [0]);
						bleft = n;
					}
					*bufp++ = j;
					++(ibp -> count);
					--bleft;
					++index;
				}
			}
		}
		*op++ = index;
	}

	/* Copy list data into nice contiguous memory. */

	ip = NEWA (index, int);
	ip1 = ip;
	for (ibp = root; ibp NE NULL; ibp = ibp -> next) {
		bufp = &(ibp -> buf [0]);
		for (i = 0; i < ibp -> count; i++) {
			*ip1++ = *bufp++;
		}
	}

	if (ip1 NE ip + index) {
		fatal ("compute_zt: Bug 1.");
	}

	/* Allocate arrays of pointers into the lists. */

	rip -> zt [0]	= NEWA (4 * (n + 1), int *);
	rip -> zt [1]	= rip -> zt [0] + (n + 1);
	rip -> zt [2]	= rip -> zt [1] + (n + 1);
	rip -> zt [3]	= rip -> zt [2] + (n + 1);
	rip -> zt [4]	= rip -> zt [0];
	rip -> zt [5]	= rip -> zt [1];
	rip -> zt [6]	= rip -> zt [2];

	op = offset;
	for (dir = 0; dir < 4; dir++) {
		zt = rip -> zt [dir];
		for (i = 0; i <= n; i++, p1++) {
			zt [i] = &(ip [*op++]);
		}
	}
	if (rip -> zt [0] [0] + index NE rip -> zt [3] [n]) {
		fatal ("compute_zt: Bug 2.");
	}

	while (root NE NULL) {
		ibp = root;
		root = ibp -> next;
		free ((char *) ibp);
	}

	free ((char *) offset);
}

/*
 * Compute the upper bounds UB1.
 */

	static
	void
compute_ub1 (

struct rinfo *		rip		/* IN/OUT - global RFST info */
)
{
int			i;
int			j;
int			n;
int			dir;
int			lr;
int			last;
struct pset *		pts;
struct point *		p1;
struct point *		p2;
struct point *		p3;
int *			succ;
int *			ip1;
int *			ip2;
int **			zt;
dstdiroff_t		dir_off;
dstdiroff_t		dirp_off;
lrdiroff_t		lrdir_off;
dist_t *		array;
dist_t *		ub1;
dist_t			d1, d2, d3, d4;
dist_t			bound;
struct point		s;

static int		lr_table [] = {0, 1, 1, 0};

	pts = rip -> pts;
	n = pts -> n;

	array = NEWA (4 * n, dist_t);

	rip -> ub1 [0]	= array;	array += n;
	rip -> ub1 [1]	= array;	array += n;
	rip -> ub1 [2]	= array;	array += n;
	rip -> ub1 [3]	= array;
	rip -> ub1 [4]	= rip -> ub1 [0];
	rip -> ub1 [5]	= rip -> ub1 [1];
	rip -> ub1 [6]	= rip -> ub1 [2];

	for (dir = 0; dir < 4; dir++) {
		lr	= lr_table [dir];

		ub1 = rip -> ub1 [dir];
		zt = rip -> zt [dir];

		succ = rip -> succ [dir];
		dir_off = DSTDIR_GETOFFSET (dir);
		dirp_off = DSTDIR_GETOFFSET (dir + 1);
		lrdir_off = LRINDEX_GETOFFSET (dir);

		p1 = &(pts -> a [0]);
		for (i = 0; i < n; i++, p1++) {
			ip1 = zt [i];
			ip2 = zt [i + 1];
			if (ip1 >= ip2) {
				/* No short leg candidates.  UB1 is zero! */
				ub1 [i] = 0.0;
				continue;
			}

			/* Get last short leg candidate in this direction. */
			last = ip2 [-1];

			bound = INF_DISTANCE;
			p3 = &(pts -> a [last]);
			/* Get Steiner point for this last terminal. */

			SPOINT (&s, p1, p3, dir);

			d3 = DSTDIR_OFFSET (p1, p3, dirp_off);
			for (j = succ [last]; j >= 0; j = succ [j]) {
				p2 = &(pts -> a [j]);
				d1 = DSTDIR_OFFSET (&s, p2, dir_off);
				if (d1 > bound) break;

				d2 = DSTDIR_OFFSET (&s, p2, dirp_off);
				if ((LRINDEX_OFFSET (p1, p2, lrdir_off) EQ lr)
				    AND
				    (d3 > d2)) {
					bound = d1;
					break;
				}
				if (d1 > d2) {
					d4 = d1 + d2; /* RDIST(&s, p2) */
					if (d4 < bound) {
						bound = d4;
					}
				}
			}
			ub1 [i] = DSTDIR_OFFSET (p1, p3, dir_off) + bound;
		}
	}
}

/*
 * Recursively grow a rectilinear FST.
 */

	static
	void
grow_RFST (

struct rinfo *	rip,		/* IN/OUT - The global RFST info */
int		size,		/* IN - number of terms in partial tree */
dist_t		length,		/* IN - length of partial tree */
int		dir,		/* IN - growth direction from root */
dist_t		ub_length,	/* IN - upper bound on partial tree length */
dist_t *	ub_shortleg,	/* IN - left/right short leg upper bounds */
int		longindex	/* IN - current long leg terminal */
)
{
int			i, j, k, l, j2;
int			lr;
int			r;
int			lastlr;
int			dirp;
int			lp;
bool			pass_bsd;
bool			need_to_restore;
struct pset *		pts;
struct point *		p;
struct point *		q;
struct point *		root;
struct point *		last;
int *			terms;
int *			succ;
int *			ip1;
int *			ip2;
dist_t			dstdir, dstdirp;
dist_t			d1, d2, d3;
dist_t			ub1;
dist_t			b;
dist_t			max_backbone;
dist_t			min_bsd;
dist_t			new_ub_length;
dist_t			new_ub_shortleg [2];
dist_t			lastdstdirp;
dist_t			long_leg_max;

#define LONGTERMS	rip -> longterms

	terms		= rip -> terms;
	r		= terms [0];
	l		= terms [size - 1];
	lastlr		= rip -> lrindex [l];
	succ		= rip -> succ [dir];

	pts		= rip -> pts;

	root = &(pts -> a [r]);
	last = &(pts -> a [l]);

	max_backbone = INF_DISTANCE;

	lastdstdirp = DSTDIRP (root, last, dir);

	/* We set this flag whenever we change one of the "maxedges" values. */
	need_to_restore = FALSE;

	for (;;) {
		i = LONGTERMS [++longindex];
		if (i < -1) {
			/* No more candidates, and we have already	*/
			/* determined that no more can be found.	*/
			break;
		}
		if (i EQ -1) {
			/* Dynamically add next candidate to longterms. */
			for (i = succ [LONGTERMS [longindex - 1]];
			     i >= 0;
			     i = succ [i]) {
				p = &(pts -> a [i]);
				dstdirp = DSTDIRP (root, p, dir);
				if (dstdirp EQ 0.0) {
					lr = 2;
					rip -> lrindex [i] = 2;
					rip -> shortterm [i] = 0;
					LONGTERMS [longindex] = i;
					LONGTERMS [longindex + 1] = -2;
					break;
				}

				lr = LRINDEX (root, p, dir);
				if (lr EQ 0) {
					dirp = (dir - 1) & 0x03;
				}
				else {
					dirp = (dir + 1) & 0x03;
				}

				/* Find short leg candidate (if it exists) */
				j = -1;
				ip1 = rip -> zt [dirp] [i];
				ip2 = rip -> zt [dirp] [i + 1];
				while (ip1 < ip2) {
					k = *--ip2;
					if (LRINDEX (root, &(pts -> a [k]), dir)
					    EQ lr) {
						j = k;
						break;
					}
				}
				rip -> shortterm [i] = j;

				/* Check upper bounds */
				ub1 = 0.0;
				if (j >= 0) {
					ub1 = rip -> ub1 [dirp] [i];
				}
				d1 = rip -> ub0 [dirp] [i];
				if (d1 < ub1) {
					d1 = ub1;
				}
				if (dstdirp > d1) continue;

				rip -> lrindex [i] = lr;
				LONGTERMS [longindex]	  = i;
				LONGTERMS [longindex + 1] = -1;
				break;
			}
			if (i < 0) {
				/* This long leg has no further candidates! */
				/* Mark it so we never try to do find any   */
				/* more candidates! */
				LONGTERMS [longindex] = -2;
				break;
			}
		}
		else {
			/* Next long leg candidate available in longterms! */
			p = &(pts -> a [i]);
			lr = rip -> lrindex [i];
			dstdirp = DSTDIRP (root, p, dir);
		}

		dstdir	 = DSTDIR (last, p, dir);

		/* Check if consecutive terminals share Steiner point. */
		if ((size >= 3) AND (dstdir EQ 0.0)) continue;

		/* Update maximum backbone length using empty diamond */
		/* property. */
		if (dstdirp < dstdir) {
			d1 = dstdir + dstdirp;
			if (d1 < max_backbone) {
				max_backbone = d1;
			}
		}

		/* Update maximum backbone length using empty corner	*/
		/* rectangle property.					*/
		if ((size >= 2) AND (lr EQ lastlr) AND (dstdirp < lastdstdirp)) {
			if (dstdir < max_backbone) {
				max_backbone = dstdir;
			}
		}

		/* Check length of new backbone segment */
		if (dstdir > max_backbone) break;

		if (lr EQ 2) {
			/* Terminal is ON the backbone!	 Save	*/
			/* immediately as type (i) and break.	*/
			if (size >= 2) {
				terms [size] = i;
				test_and_save_fst (rip,
						   size + 1,
						   length + dstdir + dstdirp,
						   dir,
						   TYPE_1);
			}
			break;
		}

		/* Terminal on the wrong side of the long leg? */
		if ((size >= 2) AND (lr EQ lastlr)) continue;

		/* Empty rectangle with last terminal? */
		if (NOT is_empty_rectangle (rip -> empty_rect, l, i)) continue;

		/* Check if new backbone segment is longer than any BSD. */
		pass_bsd = TRUE;
		min_bsd	 = INF_DISTANCE;
		for (j = 0; j < size; j++) {
			k = terms [j];
			d1 = rip -> maxedges [k];
			if (dstdir > d1) {
				d1 = dstdir;
				rip -> maxedges [k] = d1;
				need_to_restore = TRUE;
			}
			b = bsd (rip -> bsd, i, k);
			if (d1 > b) {
				pass_bsd = FALSE;
				break;
			}
			if (b < min_bsd) {
				min_bsd = b;
			}
		}
		if (NOT pass_bsd) continue;
		new_ub_length = ub_length + min_bsd;

		/* dirp = (lr EQ 0) ? dir+3 : dir+1; */
		dirp = (dir + (lr + lr - 1)) & 3;

		/* Try to make a type (ii) FST */

		j = rip -> shortterm [i];
		if (j < 0) goto try_typeI;	/* no short leg candidate! */

		/* Is backbone rectangle empty? */
		if (NOT is_empty_rectangle (rip -> empty_rect, r, j)) goto try_typeI;
		/* Check UB1. */
		if (dstdirp > rip -> ub1 [dirp] [i]) goto try_typeI;

		/* Check BSD for each terminal in current tree */
		q = &(pts -> a [j]);
		for (j2 = 0; j2 < size; j2++) {
			k = terms [j2];
			d1 = rip -> maxedges [k];
			d2 = DSTDIRP (root, q, dir);
			d3 = DSTDIRP (p, q, dir);
			if (d2 > d1) {
				d1 = d2;
			}
			/* d1 = MAX (maxedges [k], DSTDIRP (root, q, dir)) */
			if (d1 > d3) {
				d3 = d1;
			}
			/* d3 = MAX (d1, DSTDIRP (p, q, dir)) */
			if (d3 > bsd (rip -> bsd, i, k)) goto try_typeI;
			d3 = DSTDIR (p, q, dir);
			if (d1 > d3) {
				d3 = d1;
			}
			/* d3 = MAX (d1, DSTDIR (p, q, dir)) */
			if (d3 > bsd (rip -> bsd, j, k)) goto try_typeI;
		}

#if UB_SHORTLEG
		/* Check short leg upper bound */
		if (DSTDIRP (root, q, dir) > ub_shortleg [lr]) goto try_growing;
#endif

		/* Perform FST specific tests and save type (ii) tree */
		terms [size]		= i;
		terms [size + 1]	= j;
		test_and_save_fst (rip,
				   size + 2,
				   length + DSTDIR (last, q, dir)
					  + DSTDIRP (root, p, dir),
				   dir,
				   TYPE_2);

		/* Try to make a type (i) FST */
try_typeI:

		/* Check UB0 */
		if (dstdirp > rip -> ub0 [dirp] [i]) continue;

		/* Check BSD to each terminal in previous tree */
		pass_bsd = TRUE;
		for (j2 = 0; j2 < size; j2++) {
			k = terms [j2];
			d1 = rip -> maxedges [k];
			if (dstdirp > d1) {
				d1 = dstdirp;
			}
			if (d1 > bsd (rip -> bsd, k, i)) {
				pass_bsd = FALSE;
				break;
			}
		}
		if (NOT pass_bsd) continue;

		/* Check if BSD upper bound is shorter than partial tree */
		/* (including connecting Steiner point).		 */
		if (length + dstdir > new_ub_length) continue;

		/* Check if BSD upper bound is shorter than partial tree */
		if (length + dstdir + dstdirp > new_ub_length) goto try_growing;

		/* Do not make 2-terminal FSTs (MST edges). */
		if (size <= 1) goto try_growing;

		/* Is backbone rectangle empty? */
		if (NOT is_empty_rectangle (rip -> empty_rect, r, i)) goto try_growing;

#if UB_SHORTLEG
		/* Check short leg upper bound */
		if (dstdirp > ub_shortleg [lr]) goto try_growing;
#endif

		/* Perform FST specific tests and save type (i) tree */
		terms [size] = i;
		new_ub_length = test_and_save_fst (rip,
						   size + 1,
						   length + dstdir + dstdirp,
						   dir,
						   TYPE_1);

		/* Try to grow the current tree. */

try_growing:
		/* Should we generate larger FSTs? */
		if (size >= MaxFSTSize) continue;

		/* Upper bound (A). */
		d1 = ub_shortleg [lr];
		if (dstdirp < d1) {
			d1 = dstdirp;
		}
		new_ub_shortleg [lr] = d1;

		d1 = ub_shortleg [1-lr];
		if (size >= 2) {
			/* Upper bound (B). */
			d2 = rip -> ub0 [dirp] [i];
			if (min_bsd < d2) {
				d2 = min_bsd;
			}
			d2 -= dstdirp;
			if (d2 < d1) {
				d1 = d2;
			}

			/* Upper bound (C). */
			if (dstdir < d1) {
				d1 = dstdir;
			}

			if (size >= 3) {
				/* Upper bound (D). */
				lp = terms [size - 2];
				d2 = DSTDIRP (root, &(pts -> a [lp]), dir);
				if (dstdirp < d2) {
					d2 = dstdirp;
				}
				d2 = DSTDIR (&(pts -> a [lp]), p, dir) - d2;
				if (d2 < d1) {
					d1 = d2;
				}
			}
		}
		new_ub_shortleg [1-lr] = d1;

		terms [size] = i;
		rip -> maxedges [i] = dstdirp;
		grow_RFST (rip,
			   size + 1,
			   length + dstdir + dstdirp,
			   dir,
			   new_ub_length,
			   new_ub_shortleg,
			   longindex);
	}

	if (need_to_restore) {
		/* Restore caller's maxedges values.  We do this by	*/
		/* recomputing them from scratch in a scan back toward	*/
		/* the root.						*/
		long_leg_max = 0.0;
		for (i = size - 1; i > 0; i--) {
			k = terms [i];
			p = &(pts -> a [k]);
			d1 = DSTDIRP (root, p, dir);
			if (long_leg_max > d1) {
				d1 = long_leg_max;
			}
			rip -> maxedges [k] = d1;
			j = terms [i - 1];
			q = &(pts -> a [j]);
			d1 = DSTDIR (q, p, dir);
			if (d1 > long_leg_max) {
				long_leg_max = d1;
			}
		}
		rip -> maxedges [r] = long_leg_max;
	}
}

/*
 * This routine performs all of the FST specific screening tests.
 * If all are passed, the FST is saved.
 */

	static
	dist_t
test_and_save_fst (

struct rinfo *	rip,		/* IN/OUT - The global RFST info */
int		size,		/* IN - number of terminals in RFST */
dist_t		length,		/* IN - length of RFST */
int		dir,		/* IN - backbone growth direction from root */
int		type		/* IN - type of tree */
)
{
int			i, j, k;
int			z12, z13, z23;
int			nstein;
int			nedges;
int			last;
int *			terms;
struct pset *		pts;
struct point *		p1;
struct point *		p2;
struct point *		p3;
struct point *		p4;
struct rlist *		rp;
struct rlist **		hookp;
struct rlist *		rp1;
struct rlist *		rp2;
struct pset *		new_terms;
struct pset *		new_steiners;
dist_t			b;
dist_t			d1;
struct full_set *	fsp;
struct edge *		edges;
struct point		s1;
struct point		s2;

	++(rip -> fsts_checked);

	terms = rip -> terms;

	/* Is this FST too large? */
	if (size > MaxFSTSize) return (length);

	if (size > 2) {
		b = bmst_terms_length (terms, size, rip -> bsd);
		if (length >= b) return (b);
	}

	/* Simple duplicate tests for size 3 and 4 */
	if (dir EQ 1) {
		if (size EQ 3) return (length);
		if ((size EQ 4) AND (type NE TYPE_1)) return (length);
	}

	pts = rip -> pts;

	if (size > 4) {
		/* No two terminals on the long leg should share a	*/
		/* Steiner point.					*/

		for (i = 1; i < size; i++) {
			d1 = DSTDIR (&(pts -> a [terms [i-1]]),
				     &(pts -> a [terms [i]]),
				     dir);
			if (d1 EQ 0.0) return (length);
		}
	}
	else if (size EQ 4) {
		p1 = &(pts -> a [terms [0]]);
		p2 = &(pts -> a [terms [1]]);
		if (DSTDIR (p1, p2, dir) EQ 0.0) return (length);
		p3 = &(pts -> a [terms [2]]);
		p4 = &(pts -> a [terms [3]]);
		if (DSTDIR (p3, p4, dir) EQ 0.0) return (length);

		if (DSTDIR (p2, p3, dir) EQ 0.0) {
			if (DSTDIRP (p1, p4, dir) NE 0.0) {
				return (length);
			}
			type = CROSS;
		}
	}
	else if (size EQ 3) {
		/* Make sure that 3-terminal FST is not degenerate */
		p1 = &(pts -> a [terms [0]]);
		p2 = &(pts -> a [terms [1]]);
		p3 = &(pts -> a [terms [2]]);
		z12 = (DSTDIR  (p1, p2, dir) EQ 0.0) + (DSTDIRP (p1, p2, dir) EQ 0.0);
		z13 = (DSTDIR  (p1, p3, dir) EQ 0.0) + (DSTDIRP (p1, p3, dir) EQ 0.0);
		z23 = (DSTDIR  (p2, p3, dir) EQ 0.0) + (DSTDIRP (p2, p3, dir) EQ 0.0);
		if (z12 + z13 > 1) return (length);
		if (z12 + z23 > 1) return (length);
		if (z13 + z23 > 1) return (length);
	}

#if EMPTY_DIAMOND_PROPERTY

	/* Check empty diamond property for transformed FST */
	i = 0;
	last = size - 1;
	if (type EQ TYPE_2) {
		last = size - 2;
		if ((size & 1) EQ 0) {
			i = 1;		/* type (ii), even */
		}
	}
	else if ((size & 1) NE 0) {
		i = 1;			/* type (i), odd */
	}
	p1 = &(pts -> a [terms [0]]);
	p2 = &(pts -> a [terms [size-1]]);
	while (i < last) {
		/* Check skew diamond... */
		SPOINT (&s1, p1, &(pts -> a [terms [i]]), dir);
		SPOINT (&s2, p2, &(pts -> a [terms [i+1]]), dir);
		if (NOT diamond_empty (rip, &s1, &s2, terms [i], dir)) {
			return (length);
		}
		++i;
		if (i >= last) break;

		/* Check diamond for flat segment... */
		SPOINT (&s1, p2, &(pts -> a [terms [i]]), dir);
		SPOINT (&s2, p2, &(pts -> a [terms [i+1]]), dir);
		if (NOT diamond_empty (rip, &s1, &s2, terms [i], dir)) {
			return (length);
		}
		++i;
	}
#endif

#if NOSPLIT_CORNER_FLIPPED
	/* Test that corner flipped FST does not split into two or more */
	/* FSTs. Such FSTs are not needed.				*/

	d1 = DSTDIRP (&(pts -> a [terms [size-1]]),
		      &(pts -> a [terms [0]]),
		      dir);


	if (type EQ TYPE_1) {
		i = size - 3;
	}
	else {
		i = size - 4;
	}

	while (i > 0) {
		if (DSTDIRP (&(pts -> a [terms [i]]),
			     &(pts -> a [terms [0]]),
			     dir) <= d1) {
			return (length);
		}
		i -= 2;
	}
#endif

#if KAHNG_ROBINS_HEURISTIC
	if (size > 5) {
		struct pset * new_pts;
		new_pts = NEW_PSET (size);
		new_pts -> n = size;
		p1 = &(new_pts -> a [0]);
		for (i = 0; i < size; i++, p1++) {
			p2 = &(pts -> a [terms [i]]);
			p1 -> x		 = p2 -> x;
			p1 -> y		 = p2 -> y;
			p1 -> pnum	 = ((pnum_t) i);
		}
		d1 = kahng_robins_length (new_pts, length);
		free ((char *) new_pts);
		if (length > d1) return (d1);
	}
#endif

	/* General duplicate test.  We use a hash table, for speed.	*/
	/* For correctness, the hash function must not depend upon the	*/
	/* order of the terminals in the FST.  A simple checksum has	*/
	/* this property and tends to avoid favoring any one bucket.	*/

	/* Compute hash and prepare for rapid set comparison. */
	k = 0;
	for (i = 0; i < size; i++) {
		j = terms [i];
		rip -> term_check [j] = TRUE;
		k += j;
	}
	k %= pts -> n;

	hookp = &(rip -> hash [k]);
	for (;;) {
		rp = *hookp;
		if (rp EQ NULL) break;
		if (rp -> size EQ size) {
			p1 = &(rp -> fst -> terminals -> a [0]);
			for (i = 0; ; i++, p1++) {
				if (i >= size) goto found_rfst;
				if (NOT rip -> term_check [p1 -> pnum]) break;
			}
		}
		hookp = &(rp -> next);
	}

found_rfst:

	for (i = 0; i < size; i++) {
		rip -> term_check [terms [i]] = FALSE;
	}

	if (rp NE NULL) {
		/* An FST for these terminals already exists. */
		fsp = rp -> fst;
		if (fsp -> tree_len <= length) {
			return (fsp -> tree_len);
		}
		/* The new one is shorter!  Delete the old one. */
		*hookp = rp -> next;
		rp2 = rp -> forw;
		rp1 = rp -> back;
		rp2 -> back = rp1;
		rp1 -> forw = rp2;
		free ((char *) (fsp -> terminals));
		free ((char *) (fsp -> steiners));
		free ((char *) (fsp -> edges));
		free ((char *) fsp);
		free ((char *) rp);
	}

	/* Build FST graph in edge list form. */

	new_terms = NEW_PSET (size);
	new_terms -> n = size;
	for (i = 0; i < size; i++) {
		new_terms -> a [i] = pts -> a [terms [i]];
	}

	if (size <= 2) {
		nstein = 0;
		new_steiners = NULL;
		nedges = 1;
	}
	else if (type EQ CROSS) {
		nstein = 1;
		new_steiners = NEW_PSET (1);
		nedges = 4;
	}
	else {
		nstein = size - 2;
		new_steiners = NEW_PSET (nstein);
		nedges = 2 * size - 3;
	}

	edges = NEWA (nedges, struct edge);

	build_rfst_graph (rip,
			  size,
			  dir,
			  type,
			  new_terms,
			  new_steiners,
			  edges,
			  nedges);

	fsp = NEW (struct full_set);

	fsp -> next		= NULL;
	fsp -> tree_num		= 0;
	fsp -> tree_len		= length;
	fsp -> terminals	= new_terms;
	fsp -> steiners		= new_steiners;
	fsp -> nedges		= nedges;
	fsp -> edges		= edges;

	rp = NEW (struct rlist);

	rp2 = &(rip -> list);
	rp1 = rp2 -> back;
	rp -> back	= rp1;
	rp -> forw	= rp2;
	rp -> next	= rip -> hash [k];
	rp -> size	= size;
	rp -> fst	= fsp;

	rp1 -> forw	= rp;
	rp2 -> back	= rp;
	rip -> hash [k] = rp;

	return (length);
}

/*
 * Check that the diamond defined by two points is empty of terminals.
 */

	static
	bool
diamond_empty (

struct rinfo *		rip,		/* IN - global RFST info */
struct point *		p,		/* IN - first point */
struct point *		q,		/* IN - second point */
int			i,		/* IN - terminal to start scan from */
int			dir		/* IN - scan direction */
)
{
int			j;
int *			succ;
struct point *		r;
struct pset *		pts;
dstdiroff_t		dstdir_off;
dstdiroff_t		dstdirp_off;
dist_t			d;
dist_t			dstdir;
dist_t			dstdirp;

	pts = rip -> pts;
	succ = rip -> succ [dir];
	dstdir_off  = DSTDIR_GETOFFSET (dir);
	dstdirp_off = DSTDIR_GETOFFSET (dir + 1);

	d = RDIST (p, q);

	for (j = succ [i]; j >= 0; j = succ [j]) {
		r = &(pts -> a [j]);
		dstdir = DSTDIR_OFFSET (r, q, dstdir_off);
		if (dstdir > d) break;	/* finished in this direction */
		dstdirp = DSTDIR_OFFSET (r, q, dstdirp_off);
		if ((RDIST (r, p) < d) AND (dstdir + dstdirp < d)) {
			return (FALSE); /* found terminal in the diamond! */
		}
	}

	succ = rip -> succ [dir + 2];

	for (j = succ [i]; j >= 0; j = succ [j]) {
		r = &(pts -> a [j]);
		dstdir = DSTDIR_OFFSET (r, p, dstdir_off);
		if (dstdir > d) break;	/* finished in this direction */
		dstdirp = DSTDIR_OFFSET (r, p, dstdirp_off);
		if ((RDIST (r, q) < d) AND (dstdir + dstdirp < d)) {
			return (FALSE); /* found terminal in the diamond! */
		}
	}

	/* The diamond is empty! */
	return (TRUE);
}

/*
 * This routine constructs a graph of the RFST in edge list form.
 * We generate either the normal topology, or the "corner-flipped"
 * topology -- which ever yields a "left-most" and "top-most" topology.
 * This routine therefore also fills in the proper Steiner points.
 */

	static
	void
build_rfst_graph (

struct rinfo *		rip,
int			size,
int			dir,
int			type,
struct pset *		terms,
struct pset *		steins,
struct edge *		edges,
int			nedges
)
{
int			i, j, k;
struct point *		p1;
struct point *		p2;
struct point *		p3;
struct edge *		ep;

	ep = edges;

	p1 = &(terms -> a [0]);
	p2 = &(terms -> a [1]);

	if (size <= 2) {
		ep -> p1	= 0;
		ep -> p2	= 1;
		ep -> len	= RDIST (p1, p2);
		return;
	}

	if (type EQ CROSS) {
		steins -> n = 1;
		p1 = &(terms -> a [0]);
		p2 = &(terms -> a [1]);
		p3 = &(steins -> a [0]);
		SPOINT (p3, p1, p2, dir);
		p3 -> pnum = 0;
		for (i = 0; i < 4; i++) {
			ep -> p1 = i;
			ep -> p2 = 4;
			ep -> len = RDIST (&(terms -> a [i]),
					   &(steins -> a [0]));
			++ep;
		}
		return;
	}

	/* Decide whether to build a normal or corner-flipped topology. */
	/* This really only affects the Steiner point locations, not	*/
	/* the structure of the graph.					*/

	steins -> n = size - 2;

	k = size - 1;
	if (type EQ TYPE_2) {
		--k;
	}
	k = terms -> a [k].pnum;
	if ((dir EQ 1) AND
	    ((size <= 2) OR (rip -> lrindex [k] EQ 1))) {
		/* Build normal (unflipped) topology. */
		p1 = (type EQ TYPE_1) ? &(terms -> a [0])
				      : &(terms -> a [size-1]);
		p2 = &(terms -> a [size-2]);
		p3 = &(steins -> a [size - 3]);
		for (i = size - 3; i >= 0; i--) {
			SPOINT (p3, p1, p2, dir);
			p3 -> pnum = i;
			p1 = &(terms -> a [0]);
			--p2;
			--p3;
		}
	}
	else {
		/* Build corner-flipped topology. */
		if (type EQ TYPE_1) {
			k = ((size & 1) EQ 0) ? size-1 : 0;
		}
		else {
			k = ((size & 1) EQ 0) ? 0 : size-1;
		}
		p1 = &(terms -> a [k]);
		p2 = &(terms -> a [1]);
		p3 = &(steins -> a [0]);
		for (i = 0; i < size-2; i++) {
			SPOINT (p3, p1, p2, dir);
			p3 -> pnum = i;
			p1 = &(terms -> a [size - 1]);
			++p2;
			++p3;
		}
	}

	/* Now that the Steiner points are in the corrct places, build	*/
	/* the edges of the FST graph.	They have the same structure	*/
	/* regardless of whether the FST is type (i), type (ii),	*/
	/* flipped or unflipped.					*/

	j = 0;
	for (i = 0; i < size-2; i++) {
		ep -> p1	= j;
		j = size + i;
		ep -> p2	= j;
		++ep;
		ep -> p1	= j;
		ep -> p2	= i + 1;
		++ep;
	}
	ep -> p1 = j;
	ep -> p2 = i + 1;
	++ep;

	if (edges + nedges NE ep) {
		fatal ("build_rfst_graph: Bug 2.");
	}

	ep = edges;
	for (i = 0; i < nedges; i++, ep++) {
		j = ep -> p1;
		p1 = (j < size) ? &(terms -> a [j])
				: &(steins -> a [j - size]);
		k = ep -> p2;
		p2 = (k < size) ? &(terms -> a [k])
				: &(steins -> a [k - size]);
		ep -> len = RDIST (p1, p2);

		if (ep -> len EQ 0.0) {
			fatal ("build_rfst_graph: Bug 3.");
		}

	}
}

/*
 * Map all terminal numbers back to their original values.  (i.e., with
 * the duplicate terminals reinstated.	We must renumber the terminals
 * of each RFST also.
 */

	static
	void
renumber_terminals (

struct rinfo *		rip,		/* IN/OUT - global RFST info */
struct pset *		to_pts,		/* IN - point set to map to (orig) */
int *			rev_map		/* IN - map from new to old terms */
)
{
int			i;
int			j;
int			from_n;
int			to_n;
int			kmasks;
struct pset *		from_pts;
struct rlist *		rp1;
struct rlist *		rp2;
struct full_set *	fsp;
struct pset *		terms;
struct point *		p1;

	from_pts	= rip -> pts;
	from_n		= from_pts -> n;
	to_n		= to_pts -> n;

	kmasks = rip -> num_term_masks;

	/* Restore original set of terminals. */
	rip -> pts = to_pts;

	/* Renumber terminals in each FST. */
	rp2 = &(rip -> list);
	for (rp1 = rp2 -> forw; rp1 NE rp2; rp1 = rp1 -> forw) {
		fsp = rp1 -> fst;
		terms = fsp -> terminals;
		p1 = &(terms -> a [0]);
		for (i = 0; i < terms -> n; i++, p1++) {
			j = p1 -> pnum;
			if ((j < 0) OR (j >= from_n)) {
				fatal ("renumber_terminals: Bug 1.");
			}
			j = rev_map [j];
			if ((j < 0) OR (j >= to_n)) {
				fatal ("renumber_terminals: Bug 2.");
			}
			p1 -> pnum = j;
		}
	}
}

/*
 * Link all of the FSTs together into one long list and number them
 * each sequentially.  Free the doubly-linked rlists as we go.
 */

	static
	void
build_fst_list (

struct rinfo *		rip		/* IN - global RFST info */
)
{
int			i;
struct rlist *		rp1;
struct rlist *		rp2;
struct rlist *		rp3;
struct full_set *	fsp;
struct full_set **	hookp;

	hookp = &(rip -> full_sets);

	rp2 = &(rip -> list);
	i = 0;
	for (rp1 = rp2 -> forw; rp1 NE rp2; ) {
		fsp = rp1 -> fst;
		fsp -> tree_num = i++;
		*hookp = fsp;
		hookp = &(fsp -> next);
		rp3 = rp1;
		rp1 = rp1 -> forw;
		free ((char *) rp3);
	}
	*hookp = NULL;

	rip -> list.forw = rp2;
	rip -> list.back = rp2;

	/* Make it easy to add zero-length FSTs onto the end. */
	rip -> ntrees	= i;
	rip -> hookp	= hookp;
}

/*
 * This routine adds one RFST (of zero length) connecting the first
 * terminal of each duplicate terminal group to each terminal that
 * was deleted during the major processing.
 */

	static
	void
add_zero_length_fsts (

struct rinfo *		rip,		/* IN - global RFST info */
int			ndg,		/* IN - number of duplicate groups */
int **			list		/* IN - list of duplicate groups */
)
{
int			i;
int			t;
int			u;
int			kmasks;
int *			ip1;
int *			ip2;
struct pset *		pts;
struct pset *		terms;
struct edge *		edges;
struct full_set *	fsp;

	pts	= rip -> pts;
	kmasks	= rip -> num_term_masks;

	for (i = 0; i < ndg; i++) {
		ip1 = list [i];
		ip2 = list [i + 1];
		if (ip1 + 2 > ip2) {
			/* Group with fewer than two members! */
			fatal ("add_zero_length_fsts: Bug 1.");
		}
		t = *ip1++;
		while (ip1 < ip2) {
			u = *ip1++;

			terms = NEW_PSET (2);
			terms -> n = 2;
			terms -> a [0]	= pts -> a [t];
			terms -> a [1]	= pts -> a [u];

			edges = NEW (struct edge);
			edges -> p1	= 0;
			edges -> p2	= 1;
			edges -> len	= 0.0;

			fsp = NEW (struct full_set);

			fsp -> next	 = NULL;
			fsp -> tree_num	 = (rip -> ntrees)++;
			fsp -> tree_len	 = 0.0;
			fsp -> terminals = terms;
			fsp -> steiners	 = NULL;
			fsp -> nedges	 = 1;
			fsp -> edges	 = edges;

			*(rip -> hookp) = fsp;
			rip -> hookp	= &(fsp -> next);
		}
	}
}

/*
 * This routine allocates an array of pointers that lets us access a particular
 * tree number directly.
 */

	static
	struct full_set **
put_trees_in_array (

struct full_set *	fsp,		/* IN - list of full-trees. */
int *			acount		/* OUT - count of array elements. */
)
{
struct full_set **	ap;
struct full_set *	p;
int			count;
int			num;
struct full_set **	array;

	count = 0;
	for (p = fsp; p NE NULL; p = p -> next) {
		++count;
	}

	array = NEWA (count, struct full_set *);

	num = 0;
	ap = array;
	for (p = fsp; p NE NULL; p = p -> next) {
		*ap++ = p;
	}

	*acount = count;

	return (array);
}

/*
 * Routine to compute the X distance between two points.
 */

	static
	dist_t
delta_x_func (

struct point *		p1,		/* IN - point 1 */
struct point *		p2		/* IN - point 2 */
)
{
	return (fabs (p1 -> x - p2 -> x));
}


/*
 * Routine to compute the Y distance between two points.
 */

	static
	dist_t
delta_y_func (

struct point *		p1,		/* IN - point 1 */
struct point *		p2		/* IN - point 2 */
)
{
	return (fabs (p1 -> y - p2 -> y));
}

/*
 * Routines to compute whether a point Q is left or right of the ray
 * from P in a given direction.	 The index returned is 0 if Q is on
 * the left, and 1 otherwise.
 */

	static
	int
lrindex_dir_0 (

struct point *	p,
struct point *	q
)
{
	return (q -> y <= p -> y);
}

	static
	int
lrindex_dir_1 (

struct point *	p,
struct point *	q
)
{
	return (q -> x >= p -> x);
}

	static
	int
lrindex_dir_2 (

struct point *	p,
struct point *	q
)
{
	return (q -> y >= p -> y);
}

	static
	int
lrindex_dir_3 (

struct point *	p,
struct point *	q
)
{
	return (q -> x <= p -> x);
}

/*
 * Compute the CPU time used since the last time we called this routine.
 */

	static
	cpu_time_t
get_delta_cpu_time (void)

{
cpu_time_t		now;
cpu_time_t		delta;

	now = get_cpu_time ();

	delta = now - Tn;

	Tn = now;

	return (delta);
}


/*
 * Compute and format the CPU time used since the last time we called
 * this routine.
 */

	static
	void
convert_delta_cpu_time (

char *		buf		/* OUT - ASCII time string, CPU seconds */
)
{
cpu_time_t	delta;

	delta = get_delta_cpu_time ();
	convert_cpu_time (delta, buf);
}
