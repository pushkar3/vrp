/***********************************************************************

	File:	prunefst.c
	Rev:	a-1
	Date:	11/22/2000

	Copyright (c) 2000, 2001 by Martin Zachariasen & David M. Warme

************************************************************************

	Pruning of Euclidean and rectilinear FSTs using method
	originally proposed by Fossmeier & Kaufmann.
	Use implementation similar to the one suggested by Althaus,
	but with significant improvements in the compatibility tests.

************************************************************************

	Modification Log:

	a-1:	11/22/2000	martinz
		: Created.

************************************************************************/

#include "bb.h"
#include "bsd.h"
#include "btsearch.h"
#include "dsuf.h"
#include "efuncs.h"
#include "emptyr.h"
#include "p1io.h"
#include "steiner.h"

/*
 * Global Routines
 */

int			main (int, char **);

/*
 * External References
 */

	/* none */

/*
 * Local Macros
 */

#ifndef MIN
 #define MIN(a,b)	(((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
 #define MAX(a,b)	(((a) > (b)) ? (a) : (b))
#endif


/*
 * Local Types
 */

struct pg_edge {
	dist_t				len; /* Length of edge */
	int				fst; /* Corresponding FST */
	int				p1;  /* First vertex */
	int				p2;  /* Second vertex */
};

struct clt_info {
	dist_t				dist;
	int				term;
	int				aterm1;
	int				aterm2;
};

struct pinfo {
	struct cinfo *			cip;
	int				num_pg_edges;
	int				num_pg_verts;
	struct pg_edge *		pg_edges;
	int *				steiner_index;
	struct clt_info **		clt_info;
	int *				clt_count;
	bitmap_t *			compat_mask;
	double				eps;
};

struct incompat {
	struct incompat *		next;
	int				fst;
};

struct bc3 {
	struct cinfo *	cip;		/* problem data */
	int		kmasks;		/* size of vert_mask */
	int		nmasks;		/* size of edge_mask */
	int *		dfs;		/* DFS number of each vertex */
	int *		low;		/* lowest DFS num in component */
	int *		parent;		/* parents of vertices in DFS tree */
	int		max_stack;	/* size of stack */
	int *		stack;		/* base-address of edge stack */
	int *		sp;		/* current stack pointer */
	int		counter;	/* DFS number generator */
	bitmap_t *	bcc_vmask;	/* scratch buffer for new BCCs */
	bitmap_t *	bcc_emask;	/* scratch buffer for new BCCs */
	bitmap_t *	edges_seen;	/* edges already pushed */
	int *		bcc_vlist;	/* scratch vertex list buffer */
	int *		degree;		/* temp vertex degree counter */
	int *		made_req;	/* caller's list of required edges */
	int		req_count;	/* cur index into made-req */
};


/*
 * Local Routines
 */

static void		add_incompat (struct incompat **, int, int, int *);
static int		bcc_find_required (struct cinfo *, int *, int);
static void		bcc3 (struct bc3 *, int);
static int		comp_ints (const void *, const void *);
static int		comp_pg_edges (const void *, const void *);
static int		comp_clt (const void *, const void *);
static void		compute_incompatibility_info (struct pinfo *, struct bsd *);
static void		compute_pruning_info (struct cinfo *, struct pinfo *);
static void		convert_delta_cpu_time (char *);
static void		decode_params (int, char**);
static cpu_time_t	get_delta_cpu_time (void);
static bool		passes_upper_bound_tests (struct pinfo *, struct bsd *, int, int,
						  struct pset *, bitmap_t *, int *, bitmap_t *);
static void		process_bcc3 (struct bc3 *, int *, int *);
static void		prune_fsts(struct cinfo *, struct bsd *, double);
static bool		prune_this_fst (struct cinfo *, struct pinfo *, int);
static dist_t		terminal_edge_distance(struct cinfo *, struct point *,
						struct point *, struct point *,
						struct point *,
						dist_t *, dist_t *);
static void		test_close_terminal (struct pinfo *,
					     int, int, struct clt_info **);
static void		usage ();
static void		zap_deleted_fsts(struct cinfo *);


/*
 * Local Variables
 */

static char *		description = NULL;
static int		EpsilonFactor = 32;
static char *		me;
static int		output_version = CURRENT_P1IO_VERSION;
static bool		Print_Detailed_Timings = FALSE;
static cpu_time_t	T0;
static cpu_time_t	Tn;


/*
 * Function Prototypes.
 */

extern int		bmst (struct pset *, struct bsd *, struct edge *);
extern dist_t		smith_lee_liebman(struct pset *);
extern dist_t		greedy_heuristic(struct pset *, struct bsd *);

/*
 * The main routine for the "prunefst" program.	 It takes the output from
 * the Euclidean or rectilinear FST generator, prunes FSTs that cannot
 * appear in an optimal solution and outputs the pruned set of FSTs
 * to the backend (FST concatenation procedure).
 */

	int
main (

int		argc,
char **		argv
)
{
int			nedges;
int			fpsave;
char			buf1 [32];
struct cinfo		cinfo;
struct bsd *		bsd;
struct edge *		mst_edges;
cpu_time_t		Tzap;
bitmap_t *		empty_rect;
double			eps;

	fpsave = set_floating_point_double_precision ();

	setbuf (stdout, NULL);

	decode_params (argc, argv);

	init_tables ();

	/* Read FST data */
	read_phase_1_data (&cinfo);

	if ((cinfo.metric NE RECTILINEAR) AND (cinfo.metric NE EUCLIDEAN)) {
		fprintf (stderr, "Can only prune geometric (Euclidean or rectilinear) FSTs\n");
		exit (1);
	}

	startup_lp_solver ();

	initialize_btsearch ();

	T0 = get_cpu_time ();
	Tn = T0;

	/* Compute minimum spanning tree */

	mst_edges = NEWA (cinfo.pts -> n, struct edge);
	eps = 0.0;
	nedges	  = 0;
	if (cinfo.metric EQ EUCLIDEAN) {
		eps	= ((double) EpsilonFactor) * DBL_EPSILON;
		nedges	= euclidean_mst (cinfo.pts, mst_edges);
	}
	else if (cinfo.metric EQ RECTILINEAR) {
		empty_rect = init_empty_rectangles (cinfo.pts, NULL);
		nedges	   = rect_mst (cinfo.pts, mst_edges, empty_rect);
		free ((char *) empty_rect);
	}
	if (nedges NE cinfo.pts -> n - 1) {
		fatal ("prunefst: bug 1.");
	}

	if (Print_Detailed_Timings) {
		convert_delta_cpu_time (buf1);
		fprintf (stderr, "Compute MST:            %s\n", buf1);
	}
	/* Compute bottleneck Steiner distances */

	bsd = compute_bsd (nedges, mst_edges, 0);
	free ((char *) mst_edges);

	if (Print_Detailed_Timings) {
		convert_delta_cpu_time (buf1);
		fprintf (stderr, "Compute BSD:            %s\n", buf1);
	}

	/* Prune FSTs */
	prune_fsts (&cinfo, bsd, eps);

	if (Print_Detailed_Timings) {
		convert_delta_cpu_time (buf1);
		fprintf (stderr, "Pruning FSTs:           %s\n", buf1);
	}

	/* Remove deleted FSTs permanently */
	zap_deleted_fsts (&cinfo);

	/* Measure zap time.  This also sets Tn so that Tn-T0 is   */
	/* the total processing time.				   */
	Tzap = get_delta_cpu_time ();

	if (Print_Detailed_Timings) {
		convert_cpu_time (Tzap, buf1);
		fprintf (stderr, "Zap deleted FSTs:       %s\n", buf1);
		convert_cpu_time (Tn - T0, buf1);
		fprintf (stderr, "Total:                  %s\n", buf1);
	}
	cinfo.p1time += (Tn - T0);
	if (description NE NULL) {
		free ((char *) cinfo.description);
		cinfo.description = gst_strdup (description);
	}

	/* Print pruned FST data */
	print_phase_1_data (&cinfo, output_version);

	shutdown_bsd (bsd);

	shutdown_lp_solver ();

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

			case 'f':
				/* Epsilon multiplication factor */
				if (*ap EQ '\0') {
					if (argc <= 0) {
						usage ();
					}
					ap = *argv++;
					--argc;
				}
				EpsilonFactor = atoi (ap);
				ap = "";
				break;
			case 't':
				Print_Detailed_Timings = TRUE;
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
	"\t-f F\tEpsilon multiplication factor F for floating point number",
	"\t\tcomparisons (default: 32).",
	"\t-t\tPrint detailed timings on stderr.",
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
			"\nUsage: %s [-t] [-d description] [-f F]",
			me);
	pp = &arg_doc [0];
	while ((p = *pp++) NE NULL) {
		(void) fprintf (stderr, "%s\n", p);
	}
	exit (1);
}

/*
 * Main pruning procedure. We use a method proposed by Fossmeier and Kaufmann
 * based on thesing whether is advantageous to extend a given FST with a
 * terminal no currently spanned. If so, the FST is discarded.
 */
	static
	void
prune_fsts (

struct cinfo *		cip,	/* IN/OUT - compatibility info */
struct bsd *		bsd,	/* IN	  - BSD data structure	*/
double			eps	/* IN	  - epsilon value */
)
{
int			i, j, t, r1, r2;
int			tcomp;
int			fsave;
int			pruned_total;
int			required_total;
int			old_pruned_total;
int			scan;
int			nverts;
int			kmasks;
int			nedges;
int			nmasks;
int			adj_edge;
int			del_comps_count;
int *			ep1;
int *			ep2;
int *			vp;
int *			vp1;
int *			vp2;
int *			comps_edge;
int *			made_req;
int			req_count;
bool			first_2edge;
bool			this_2edge;
struct dsuf		comps;
struct dsuf		del_comps;
struct pinfo		pinfo;
struct pset *		ltlist;
bitmap_t *		ltmask;
int *			lflist;
bitmap_t *		lfmask;

	nverts = cip -> num_verts;
	kmasks = cip -> num_vert_masks;
	nmasks = cip -> num_edge_masks;
	nedges = cip -> num_edges;

	pinfo.cip = cip;
	pinfo.eps = eps;
	pinfo.compat_mask = NEWA (nmasks, bitmap_t);
	memset (pinfo.compat_mask, 0, nmasks * sizeof (bitmap_t));

	/* The following are only used for calling upper bound procedure */
	ltlist	  = NEW_PSET (nverts);
	ltmask	  = NEWA (kmasks, bitmap_t);
	lflist	  = NEWA (nedges, int);
	lfmask	  = NEWA (nmasks, bitmap_t);

	/* Initialize masks */
	for (i = 0; i < kmasks; i++) {
		ltmask [i] = 0;
	}
	for (i = 0; i < nmasks; i++) {
		lfmask [i] = 0;
	}

	/* Perform thorough upper bound test for every not-yet pruned FST */
	pruned_total = 0;
	for (i = 0; i < nedges; i++) {

		if (BITON (cip -> initial_edge_mask, i)) {
			if (NOT passes_upper_bound_tests (&pinfo, bsd, i, -1,
							  ltlist, ltmask, lflist, lfmask)) {
				CLRBIT (cip -> initial_edge_mask, i);
				pruned_total++;
			}
			else {
				SETBIT (pinfo.compat_mask, i);
			}
		}
		else {
			pruned_total++;
		}
	}

	/* Compute basic compatibility */
	/* CHANGE: USE COMP INFO THAT IS ALREADY THERE? */
	compute_incompatibility_info (&pinfo, bsd);

	/* Compute pruning information */
	compute_pruning_info(cip, &pinfo);

	/* Build union-find structure */
	dsuf_create (&comps, cip -> num_verts);
	for (t = 0; t < cip -> num_verts; t++) {
		dsuf_makeset (&comps, t);
	}

	comps_edge = NEWA (cip -> num_verts, int);
	made_req   = NEWA (nedges, int);

	/* Perform actual pruning */
	if (Print_Detailed_Timings) {
		fprintf(stderr, "- scan 0 finished. %6d FSTs pruned\n", pruned_total);
	}
	required_total = 0;
	for (scan = 1; scan < nedges; scan++) {

		old_pruned_total   = pruned_total;
		for (i = 0; i < nedges; i++) {
			if ((BITON (cip -> initial_edge_mask, i)) AND
			   (NOT BITON (cip -> required_edges, i))) {

				/* Set up mask of compatible FSTs */
				ep1 = cip -> inc_edges [i];
				ep2 = cip -> inc_edges [i + 1];
				while (ep1 < ep2) {
					j = *ep1++;
					CLRBIT (pinfo.compat_mask, j);
				}

				/* Test if FST can be pruned */
				if (prune_this_fst(cip, &pinfo, i)) {
					CLRBIT (cip -> initial_edge_mask, i);
					CLRBIT (pinfo.compat_mask, i);
					pruned_total++;
				}

				/* Reset mask */
				ep1 = cip -> inc_edges [i];
				while (ep1 < ep2) {
					j = *ep1++;
					if (BITON (cip -> initial_edge_mask, j))
						SETBIT (pinfo.compat_mask, j);
				}
			}
		}

		/* Test if any connected component (initially one
		   for each terminal) only has one adjacent FST */

	try_again:
		req_count = 0;
		for (t = 0; t < cip -> num_verts; t++) {
			comps_edge [t] = -1;
		}

		/* First find the number of adjacent FSTs... */
		for (t = 0; t < cip -> num_verts; t++) {
			tcomp = dsuf_find (&comps, t); /* t's component */
			if (comps_edge [tcomp] EQ -2) continue;
			ep1 = cip -> term_trees [t];
			ep2 = cip -> term_trees [t + 1];
			while (ep1 < ep2) {
				i = *ep1++;
				if ((BITON (cip -> initial_edge_mask, i)) AND
				    (NOT BITON (cip -> required_edges, i))) {
					adj_edge = comps_edge [tcomp];
					if (adj_edge EQ -2) break;
					if (adj_edge EQ -1) {
						/* First FST */
						comps_edge [tcomp] = i;
					}
					else {
						/* Second (or more) FST */
						first_2edge = (cip -> edge_size [adj_edge] EQ 2);
						this_2edge  = (cip -> edge_size [i] EQ 2);

						/* If all adjacent edges have been 2-edges, */
						/* then pick the shortest (if equal then take edge */
						/* with minimum index) */

						if (first_2edge AND this_2edge) {
							if ((cip -> cost [i] < cip -> cost [adj_edge]) OR
							    ((cip -> cost [i] EQ cip -> cost [adj_edge]) AND
							     (i < adj_edge))) {
								comps_edge [tcomp] = i;
							}
						}
						else {
							comps_edge [tcomp] = -2; break;
						}
					}
				}
			}
		}

		/* ... then check the counts */
		for (t = 0; t < cip -> num_verts; t++) {
			fsave = comps_edge [t];
			if ((fsave >= 0) AND
			    (NOT BITON (cip -> required_edges, fsave))) {

				if (NOT BITON (cip -> initial_edge_mask, fsave)) {
					/* fatal: FST already removed */
					fatal("prune_fsts: bug 1.");
				}

				SETBIT (cip -> required_edges, fsave);
				required_total++;
				made_req [req_count++] = fsave;
			}
		}

#if 1
		i = req_count;
		req_count = bcc_find_required (cip, made_req, req_count);
		required_total += (req_count - i);
#else
		/* Test if leaving out an FST disconnects the remaining FSTs */
		/* CHANGE! This is rather time-consuming. Should be implemented using bicomps */

		for (i = 0; i < nedges; i++) {
			if ((BITON (cip -> initial_edge_mask, i)) AND
			   (NOT BITON (cip -> required_edges, i))) {

				/* create union-find data structure */
				dsuf_create (&del_comps, cip -> num_verts);
				for (t = 0; t < cip -> num_verts; t++) {
					dsuf_makeset (&del_comps, t);
				}

				del_comps_count = cip -> num_verts;

				for (j = 0; j < nedges; j++) {
					if (BITON (cip -> initial_edge_mask, j)) {

						if (j EQ i) continue; /* skip FST being tested */

						/* Unite vertices spanned */
						vp1 = cip -> edge [j];
						vp2 = cip -> edge [j + 1] - 1;
						while (vp1 < vp2) {
							r1 = dsuf_find (&del_comps, *vp1);
							r2 = dsuf_find (&del_comps, *vp2);

							if (r1 NE r2) {
								dsuf_unite (&del_comps, r1, r2);
								del_comps_count--;
							}
							vp1++;
						}
						if (del_comps_count EQ 1) break;
					}
				}

				if (del_comps_count > 1) {

					/* FST must be included in SMT */

					SETBIT (cip -> required_edges, i);
					required_total++;
					made_req [req_count++] = i;
				}
				dsuf_destroy (&del_comps);
			}
		}
#endif

		/* Now update data structures for new required FSTs */
		for (j = 0; j < req_count; j++) {

			fsave = made_req [j];

			/* Unite vertices spanned */
			vp1 = cip -> edge [fsave];
			vp2 = cip -> edge [fsave + 1] - 1;
			while (vp1 < vp2) {
				r1 = dsuf_find (&comps, *vp1);
				r2 = dsuf_find (&comps, *vp2);

				if (r1 EQ r2) { /* fatal: cycle created */
					fatal("prune_fsts: bug 2.");
				}
				dsuf_unite (&comps, r1, r2);
				vp1++;
			}

			/* Prune all incompatible FSTs */
			ep1 = cip -> inc_edges [fsave];
			ep2 = cip -> inc_edges [fsave + 1];
			while (ep1 < ep2) {
				i = *ep1++;
				if (BITON (cip -> initial_edge_mask, i)) {

					if (BITON (cip -> required_edges, i)) {
						/* fatal: FST already required */
						fatal("prune_fsts: bug 3.");
					}
					CLRBIT (cip -> initial_edge_mask, i);
					pruned_total++;
				}
			}
		}

		/* Remove FSTs making cycles among required FSTs */
		for (i = 0; i < nedges; i++) {
			if ((BITON (cip -> initial_edge_mask, i)) AND
			    (NOT BITON (cip -> required_edges, i))) {

				/* Check if a pair of vertices span the same component */
				vp1 = cip -> edge [i];
				vp2 = cip -> edge [i + 1];
				while (vp1 < vp2) {
					vp = vp1 + 1;
					while (vp < vp2) {
						r1 = dsuf_find (&comps, *vp1);
						r2 = dsuf_find (&comps, *vp);

						if (r1 EQ r2) { /* cycle created - remove FST */
							CLRBIT (cip -> initial_edge_mask, i);
							pruned_total++;
							vp = vp1 = vp2;
							break;
						}
						vp++;
					}
					vp1++;
				}
			}
		}

		if (req_count > 0) goto try_again;

		if (old_pruned_total EQ pruned_total) break;

		if (Print_Detailed_Timings) {
			fprintf(stderr, "- scan %d finished. %6d FSTs pruned\n",
				scan, pruned_total - old_pruned_total);
		}
	}

	if (Print_Detailed_Timings) {
		fprintf(stderr, "- pruning finished: before: %d  after: %d  required: %d\n",
			nedges, nedges - pruned_total, required_total);
	}

	/* Free pruning data information... */

	for (i = 0; i < nedges; i++) {
		if (pinfo.clt_info[i] NE NULL) {
			free ((char *) pinfo.clt_info [i]);
		}
	}
	free ((char *) pinfo.clt_count);
	free ((char *) pinfo.clt_info);
	free ((char *) pinfo.steiner_index);
	free ((char *) pinfo.pg_edges);

	dsuf_destroy (&comps);

	free ((char *) made_req);
	free ((char *) comps_edge);

	free ((char *) lfmask);
	free ((char *) lflist);
	free ((char *) ltmask);
	free ((char *) ltlist);

	free ((char *) pinfo.compat_mask);
}

/*
 * Remove FSTs permanently that have been marked as not needed.
 * However, MST edges (2-terminal FSTs) are not deleted even if
 * marked as never needed.
 * It should be noted that no arrays are reallocated; data is packed in place.
 */
	static
	void
zap_deleted_fsts (

struct cinfo *	cip		/* IN/OUT - compatibility info */
)
{
int			i;
int			j;
int			k;
int			new_nedges;
int *			ni;
int *			vp;
int *			vp1;
int *			vp2;
int *			ep;
int *			ep1;
int *			ep2;
struct full_set *	fsp;

	/* First we count the number of FSTs that remain */
	/* and set up map from old to new edge index */

	new_nedges = 0;
	ni  = NEWA (cip -> num_edges, int);
	for (i = 0; i < cip -> num_edges; i++) {
		if ((BITON (cip -> initial_edge_mask, i)) OR (cip -> edge_size [i] EQ 2))
			ni[i] = new_nedges++;
		else
			ni[i] = -1;
	}

	/* Pack edge_size and cost arrays */
	for (i = 0; i < cip -> num_edges; i++) {
		if (ni[i] >= 0) {
			cip -> edge_size [ ni[i] ] = cip -> edge_size [i];
			cip -> cost	 [ ni[i] ] = cip -> cost [i];
		}
	}

	/* Pack edge array */
	vp = cip -> edge [0];
	for (i = 0; i < cip -> num_edges; i++) {
		if (ni[i] >= 0) {
			vp1 = cip -> edge [i];
			vp2 = cip -> edge [i + 1];
			cip -> edge [ ni[i] ] = vp;
			while (vp1 < vp2)
				*(vp++) = *(vp1++);
		}
	}
	cip -> edge [ new_nedges ] = vp;

	/* Pack term_trees array */
	ep = cip -> term_trees [0];
	for (j = 0; j < cip -> num_verts; j++) {
		vp1 = cip -> term_trees [j];
		vp2 = cip -> term_trees [j + 1];
		cip -> term_trees [j] = ep;
		while (vp1 < vp2) {
			i = *(vp1++);
			if (ni[i] >= 0)
				*(ep++) = ni[i];
		}
	}
	cip -> term_trees[ cip -> num_verts ] = ep;

	/* Pack inc_edges array */
	ep = cip -> inc_edges [0];
	for (i = 0; i < cip -> num_edges; i++) {
		if (ni[i] >= 0) {
			ep1 = cip -> inc_edges [i];
			ep2 = cip -> inc_edges [i + 1];
			cip -> inc_edges [ ni[i] ] = ep;
			while (ep1 < ep2) {
				k = *(ep1++);
				if (ni[k] >= 0)
					*(ep++) = ni [k];
			}
		}
	}
	cip -> inc_edges [ new_nedges ] = ep;

	/* Pack full_trees array */
	if (cip -> full_trees NE NULL) {
		for (i = 0; i < cip -> num_edges; i++) {
			if (ni[i] >= 0) {
				cip -> full_trees [ ni[i] ] = cip -> full_trees [i];
				cip -> full_trees [ ni[i] ] -> tree_num = ni[i];
			}
			else {

				/* Free FST data */
				fsp = cip -> full_trees [i];
				free ((char *) (fsp -> terminals));
				free ((char *) (fsp -> steiners));
				free ((char *) (fsp -> edges));
				free ((char *) fsp);
			}
		}
	}

	/* Pack bit maps */
	for (i = 0; i < cip -> num_edges; i++) {
		if (ni[i] >= 0) {
			if (BITON (cip -> initial_edge_mask, i))
				SETBIT (cip -> initial_edge_mask, ni[i]);
			else
				CLRBIT (cip -> initial_edge_mask, ni[i]);

			if (BITON (cip -> required_edges, i))
				SETBIT (cip -> required_edges, ni[i]);
			else
				CLRBIT (cip -> required_edges, ni[i]);
		}
	}

	/* Finally set edge count */
	cip -> num_edges      = new_nedges;
	cip -> num_edge_masks = BMAP_ELTS (new_nedges);

	free((char *) ni);
}

/*
 * Check if a given FST can be pruned
 */
	static
	bool
prune_this_fst (

struct cinfo *	cip,		/* IN - compatibility info */
struct pinfo *	pip,		/* IN - pruning data structure */
int		fst
)
{
int			i;
int			clt_count;
int			curr_clt;
int			root;
int			root1;
int			root2;
bool			prune_fst;
struct dsuf		comps;
struct pg_edge *	pg_edge;
struct clt_info *	clt;

	clt_count = pip -> clt_count [fst];
	if (clt_count == 0) return FALSE;

	/* Create disjoint set */
	dsuf_create (&comps, pip -> num_pg_verts + 1);
	for (i = 0; i < pip -> num_pg_verts + 1; i++)
		dsuf_makeset (&comps, i);

	/* Add pruning graph edges (in sorted order) */

	prune_fst = FALSE;
	curr_clt  = 0;
	clt = &(pip -> clt_info [fst][curr_clt]);
	for (i = 0; i < pip -> num_pg_edges; i++) {
		pg_edge = &(pip -> pg_edges[i]);
		while (pg_edge -> len > clt -> dist) {
			root  = dsuf_find (&comps, clt -> term);
			root1 = dsuf_find (&comps, clt -> aterm1);
			root2 = dsuf_find (&comps, clt -> aterm2);

			if ((root NE root1) AND (root NE root2)) {
				prune_fst = TRUE; /* This FST can be pruned! */
				goto prune_exit;
			}

			curr_clt++;
			if (curr_clt >= clt_count)
				goto prune_exit; /* This FST cannot be pruned... */
			clt = &(pip -> clt_info [fst][curr_clt]);
		}

		/* Add edge if FST is not deleted or incompatible */
		if (BITON (pip -> compat_mask, pg_edge -> fst)) {
			root1 = dsuf_find (&comps, pg_edge -> p1);
			root2 = dsuf_find (&comps, pg_edge -> p2);
			if (root1 NE root2)
				dsuf_unite (&comps, root1, root2);
		}

	}
	prune_exit:

	dsuf_destroy (&comps);
	return prune_fst;
}

/*
 * Add a pair of FSTs as incompatible
 */

	static
	void
add_incompat (

struct incompat **	incompat,	/* IN/OUT - incomp. data structure */
int			fst1,		/* IN	  - first FST */
int			fst2,		/* IN	  - second FST */
int*			counts		/* IN/OUT - incomp. counts */
)
{
struct incompat *	icp;

	icp = incompat [fst1];
	incompat [fst1] = NEW (struct incompat);
	incompat [fst1] -> fst	= fst2;
	incompat [fst1] -> next = icp;
	counts [fst1]++;

	icp = incompat [fst2];
	incompat [fst2] = NEW (struct incompat);
	incompat [fst2] -> fst	= fst1;
	incompat [fst2] -> next = icp;
	counts [fst2]++;
}

/*
 * Computes for each FST, a list of those FSTs which are incompatible,
 * that is, fulfill one of the following conditions:
 * 1. The FSTs have two ore more terminals in common
 * 2. The FSTs have one terminal in common and the BSD of their
 * terminals is shorter than the total length of the FSTs
 * 3. An heuristic tree spanning the terminals is shorter.
 * 4. Another pair of FSTs has shorther or equal length.
 * 5. The MSTHG for the FSTs within the terminals spanned
 *    together with BSD-MST edges is shorter.
 */

	static
	void
compute_incompatibility_info (

struct pinfo *		pip,	/* IN/OUT - pruning info */
struct bsd *		bsd	/* IN	  - BSD data structure	*/
)
{
int			i;
int			j;
int			k;
int			t;
int			fs;
int			common;
int			comterm;
int			nverts;
int			nedges;
int			kmasks;
int			nmasks;
int			total;
bitmap_t *		fmask;
bitmap_t *		lfmask;
bitmap_t *		edge_mask;
bitmap_t *		tmask;
bitmap_t *		ltmask;
struct pset *		ltlist;
struct incompat *	icp;
struct incompat *	icpn;
struct incompat **	incompat;
int *			ep1;
int *			ep2;
int *			vp1;
int *			vp2;
int *			flist;
int *			lflist;
int **			inc_edges;
int *			counts;
struct cinfo *		cip;

	/* Initialize and allocate various variables and arrays */
	cip	  = pip -> cip;
	nverts	  = cip -> num_verts;
	kmasks	  = cip -> num_vert_masks;
	nedges	  = cip -> num_edges;
	nmasks	  = cip -> num_edge_masks;
	edge_mask = cip -> initial_edge_mask;

	inc_edges = NEWA (nedges + 1, int *);
	counts	  = NEWA (nedges, int);
	incompat  = NEWA (nedges, struct incompat *);
	flist	  = NEWA (nedges, int);
	fmask	  = NEWA (nmasks, bitmap_t);
	tmask	  = NEWA (kmasks, bitmap_t);
	ltlist	  = NEW_PSET (nverts);
	ltmask	  = NEWA (kmasks, bitmap_t);
	lflist	  = NEWA (nedges, int);
	lfmask	  = NEWA (nmasks, bitmap_t);

	for (i = 0; i < kmasks; i++) {
		tmask [i] = 0;
		ltmask [i] = 0;
	}
	for (i = 0; i < nmasks; i++) {
		fmask [i] = 0;
		lfmask [i] = 0;
	}
	for (i = 0; i < nedges; i++) {
		incompat [i] = NULL;
		counts [i]   = 0;
	}

	/* Compute the list of (lists of) incomatible FSTs... */
	total = 0;

	if (Print_Detailed_Timings) {
		fprintf(stderr, "- computing incompatible FSTs for each FST\n");
	}
	for (i = 0; i < nedges; i++) {
		if (NOT BITON (edge_mask, i)) continue;

		/* Develop list of all FSTs adjacent to FST i... */
		k = 0; j = 0;
		vp1 = cip -> edge [i];
		vp2 = cip -> edge [i + 1];
		while (vp1 < vp2) {
			t = *vp1++;
			SETBIT (tmask, t);
			ep1 = cip -> term_trees [t];
			ep2 = cip -> term_trees [t + 1];
			while (ep1 < ep2) {
				fs = *ep1++;
				if (BITON (fmask, fs)) continue;
				if (NOT BITON (edge_mask, fs)) continue;
				if (fs <= i) continue; /* test pairs once */
				SETBIT (fmask, fs);
				flist [k++] = fs;
			}
		}

		/* Now loop through all adjacent FSTs */
		ep1 = &flist [0];
		ep2 = &flist [k];
		while (ep1 < ep2) {
			fs = *ep1++;
			CLRBIT (fmask, fs);

			/* Count number of vertices in common. */
			common = 0; comterm = 0;
			vp1 = cip -> edge [fs];
			vp2 = cip -> edge [fs + 1];
			while (vp1 < vp2) {
				t = *vp1++;
				if (BITON (tmask, t)) {
					++common; comterm = t;
				}
			}
			if (common >= 2) {
				/* Too many - retain as incompatible */
				add_incompat(incompat, i, fs, counts); total += 2;
				continue;
			}

			/* One terminal in common. Perform thorough upper tests. */
			if (NOT passes_upper_bound_tests (pip, bsd, i, fs,
							  ltlist, ltmask, lflist, lfmask)) {
				/* Did not pass - retain as incompatible */
				add_incompat(incompat, i, fs, counts); total += 2;
				continue;
			}
		}

		vp1 = cip -> edge [i];
		vp2 = cip -> edge [i + 1];
		while (vp1 < vp2) {
			j = *vp1++;
			CLRBIT (tmask, j);
		}
	}

	/* Now allocate and copy into contiguous memory... */
	ep1 = NEWA (total, int);
	for (i = 0; i < nedges; i++) {
		inc_edges [i] = ep1;
		k = counts [i];
		if (k <= 0) continue;
		icp = incompat [i];
		while (icp NE NULL) {
			*ep1++ = icp -> fst;
			icpn = icp -> next;
			free ((char *) icp);
			icp = icpn;
		}
		qsort(inc_edges [i], k, sizeof(int), comp_ints);
	}
	inc_edges [i] = ep1;

	if (ep1 - inc_edges [0] NE total) {
		fatal ("compute_incompatibility_info: bug 1.");
	}
	if (cip -> inc_edges NE NULL) {
		if (cip -> inc_edges [0] NE NULL) {
			free ((char *) (cip -> inc_edges [0]));
		}
		free ((char *) (cip -> inc_edges));
	}
	cip -> inc_edges = inc_edges;

	/* Free allocated memory */

	free ((char *) lfmask);
	free ((char *) lflist);
	free ((char *) ltmask);
	free ((char *) ltlist);
	free ((char *) tmask);
	free ((char *) fmask);
	free ((char *) flist);
	free ((char *) incompat);
	free ((char *) counts);
}

/*
 * Compute pruning information. For every FST identify a list
 * of "close" terminals and find their distance to the FST.
 */

	static
	void
compute_pruning_info (

struct cinfo *		cip,	/* IN	  - compatibility info */
struct pinfo *		pip	/* IN/OUT - pruning data structure  */
)
{
int			i;
int			j;
int			k;
int			t;
int			total;
int			steiner_index;
int			nedges;
int			nverts;
int			kmasks;
int *			vp1;
int *			vp2;
bitmap_t *		tmask;
struct full_set *	fsp;
struct pset *		terms;
struct clt_info*	cli;
struct clt_info*	clip;
struct point *		p1;
struct point *		p2;
dist_t			l;

	nedges = cip -> num_edges;
	nverts = cip -> num_verts;
	kmasks = cip -> num_vert_masks;

	/* Terminal mask */
	tmask = NEWA (kmasks, bitmap_t);
	for (i = 0; i < kmasks; i++) {
		tmask [i] = 0;
	}

	/* Generate list of all edges in all FSTs */
	total = 0;
	for (i = 0; i < nedges; i++) {
		total += cip -> full_trees [i] -> nedges;
	}

	steiner_index = nverts;
	pip -> num_pg_edges	= total;
	pip -> pg_edges		= NEWA (total, struct pg_edge);
	memset (pip -> pg_edges, 0, total * sizeof (pip -> pg_edges [0]));
	pip -> steiner_index	= NEWA (nedges, int);
	k = 0;
	l = 0.0;
	for (i = 0; i < nedges; i++) {
		fsp = cip -> full_trees [i];
		terms = fsp -> terminals;
		pip -> steiner_index [i] = steiner_index;
		for (j = 0; j < fsp -> nedges; j++) {

			/* Compute length of this edge */
			p1 = (fsp -> edges[j].p1 < terms -> n)
			     ? &(terms -> a[ fsp -> edges[j].p1 ])
			     : &(fsp -> steiners -> a [fsp -> edges[j].p1 - terms -> n]);
			p2 = (fsp -> edges[j].p2 < terms -> n)
			     ? &(terms -> a[ fsp -> edges[j].p2 ])
			     : &(fsp -> steiners -> a [fsp -> edges[j].p2 - terms -> n]);
			if (cip -> metric EQ EUCLIDEAN) {
				l = EDIST(p1, p2) *
				    (1.0 - pip -> eps * ((double) cip -> edge_size [i]));
			}
			if (cip -> metric EQ RECTILINEAR) {
				l = RDIST(p1, p2);
			}

			pip -> pg_edges[k].fst = i;
			pip -> pg_edges[k].len = l;
			pip -> pg_edges[k].p1  =
				(fsp -> edges[j].p1 < terms -> n)
				? terms -> a[ fsp -> edges[j].p1 ].pnum
				: fsp -> edges[j].p1 - terms -> n + steiner_index;
			pip -> pg_edges[k].p2  =
				(fsp -> edges[j].p2 < terms -> n)
				? terms -> a[ fsp -> edges[j].p2 ].pnum
				: fsp -> edges[j].p2 - terms -> n + steiner_index;
			k++;
		}
		if (fsp -> steiners NE NULL) {
			steiner_index += fsp -> steiners -> n;
		}
	}
	pip -> num_pg_verts = steiner_index;

	/* Sort edge list */
	qsort(pip -> pg_edges, pip -> num_pg_edges,
	      sizeof(struct pg_edge), comp_pg_edges);

	/* Pruning graph has been constructed.
	   Now construct close terminal lists */

	pip -> clt_info	 = NEWA (nedges, struct clt_info *);
	pip -> clt_count = NEWA (nedges, int);
	cli = NEWA (nverts, struct clt_info);
	memset (cli, 0, nverts * sizeof (cli [0]));
	if (Print_Detailed_Timings) {
		fprintf(stderr, "- constructing close terminal list for each FST\n");
	}
	for (i = 0; i < nedges; i++) {

		/* Mark terminals in current FST */
		vp1 = cip -> edge [i];
		vp2 = cip -> edge [i + 1];
		while (vp1 < vp2) {
			t = *vp1++;
			SETBIT (tmask, t);
		}

		/* Find all close terminals and add information
		   to close terminal data structure */
		clip = cli;
		for (t = 0; t < nverts; t++) {
			if (NOT BITON (tmask, t)) {
				test_close_terminal(pip, i, t, &clip);
			}
		}

		/* Unmark terminals in current FST */
		vp1 = cip -> edge [i];
		while (vp1 < vp2) {
			t = *vp1++;
			CLRBIT (tmask, t);
		}

		/* Sort close terminals */
		pip -> clt_count [i] = clip - cli;
		pip -> clt_info [i]  = NULL;
		if (pip -> clt_count [i] > 0) {
			qsort(cli, pip -> clt_count [i],
			      sizeof(struct clt_info), comp_clt);
			pip -> clt_info [i] = NEWA (pip -> clt_count [i],
						    struct clt_info);
			for (j = 0; j < pip -> clt_count [i]; j++)
				pip -> clt_info [i][j] = cli[j];
		}
	}

	free ((char *) cli);
	free ((char *) tmask);
}
/*
 * Computes upper bounds for single FST or pair of FSTs.
 * Returns TRUE if all upper bound tests are passed
 * and FALSE otherwise.
 */

	static
	bool
passes_upper_bound_tests (

struct pinfo *		pip,	/* IN - pruning info */
struct bsd *		bsd,	/* IN - BSD data structure */
int			fst1,	/* IN - first FST */
int			fst2,	/* IN - second FST */
struct pset *		ltlist, /* IN - local terminal list (should just be allocated) */
bitmap_t *		ltmask, /* IN - local terminal mask (should just be cleared) */
int *			lflist, /* IN - local FST list (should just be allocated) */
bitmap_t *		lfmask	/* IN - local FST mask (should just be cleared) */
)
{
int			i;
int			j;
int			t;
int			hid;
int			nverts;
int			nedges;
int			kmasks;
int			nmasks;
int			lnedges;
int			lfs;
int			lfs2;
int			lsmtcount;
int			fstmaxind;
int			isterms;
bool			all_spanned;
bool			shorter_pair_found;
struct bbinfo *		bbip;
struct cinfo		lcinfo;
struct edge *		bsdmst;
struct edge *		ep;
int *			hep;
int *			ep1;
int *			ep2;
int *			vp1;
int *			vp2;
int *			lterm;
struct cinfo *		cip;
bitmap_t		smt_mask;
dist_t			l;
dist_t			bsdl;
dist_t			msthgl;
dist_t			minl;
dist_t			pairl;

	cip = pip -> cip;

	/* Construct list of terminals and set terminal mask */
	i = 0;
	vp1 = cip -> edge [fst1];
	vp2 = cip -> edge [fst1 + 1];
	while (vp1 < vp2) {
		t = *vp1++;
		ltlist -> a[i++] = cip -> pts -> a[t];
		SETBIT (ltmask, t);
	}

	if (fst2 >= 0) { /* negative means not defined */
		vp1 = cip -> edge [fst2];
		vp2 = cip -> edge [fst2 + 1];
		while (vp1 < vp2) {
			t = *vp1++;
			if (BITON (ltmask, t)) continue;
			ltlist -> a[i++] = cip -> pts -> a[t];
			SETBIT (ltmask, t);
		}
	}
	ltlist -> n = i;

	/* Reset terminal masks (in case we quit) */
	for (i = 0; i < ltlist -> n; i++) {
		CLRBIT (ltmask, ltlist -> a[i].pnum);
	}

	/* Lower bound on total length of these two FSTs (or this single FST)*/
	if (fst2 >= 0) {
		l = cip -> cost [fst1] + cip -> cost [fst2];
	}
	else {
		l = cip -> cost [fst1];
	}

	l *= (1.0 - pip -> eps * ((double) ltlist -> n));

	/* Compute BSD-MST */
	bsdmst = NEWA (ltlist -> n - 1, struct edge);
	if (bmst (ltlist, bsd, bsdmst) NE ltlist -> n - 1)
		fatal ("passes_upper_bound_test: bug 1");

	/* BSD-MST test... */
	bsdl = 0.0;
	ep   = bsdmst;
	for (i = 0; i < ltlist -> n - 1; i++, ep++) {
		bsdl += ep -> len;
	}
	if ((ltlist -> n > 3) AND (bsdl <= l)) {
		free ((char *) bsdmst);
		return FALSE;
	}

	/* Heuristic upper bound test (rectilinear) ... */
	if ((cip -> metric EQ RECTILINEAR) AND (ltlist -> n <= 15)) {
		if (kahng_robins_length (ltlist, l) < l) {
			free ((char *) bsdmst);
			return FALSE;
		}
	}

	/* Heuristic upper bound test (Euclidean) ... */
	if (cip -> metric EQ EUCLIDEAN) {
		if (greedy_heuristic (ltlist, bsd) < l) {
			free ((char *) bsdmst);
			return FALSE;
		}
	}

	/* If we are testing a pair of MST edges, we stop here */
	if (ltlist -> n EQ 3) {
		free ((char *) bsdmst);
		return TRUE;
	}

	/* Set terminal masks again */
	for (i = 0; i < ltlist -> n; i++) {
		SETBIT (ltmask, ltlist -> a[i].pnum);
	}

	/* Prepare for calling the branch-and-cut MSTHG procedure.	   */
	/* Construct list of FSTs spaning a subset of the given terminals. */
	/* We only need FSTs of cardinality 3 or larger.		   */

	lnedges = 0;
	for (i = 0; i < ltlist -> n; i++) {
		ep1 = cip -> term_trees [ ltlist -> a[i].pnum ];
		ep2 = cip -> term_trees [ ltlist -> a[i].pnum + 1 ];
		while (ep1 < ep2) {
			lfs = *ep1++;
			if (BITON (lfmask, lfs)) continue;
			if (NOT BITON (cip -> initial_edge_mask, lfs)) continue;
			if (cip -> edge_size [lfs] <= 2) continue;
			if ((lfs EQ fst1) AND (fst2 < 0)) continue; /* skip FST being tested */

			/* Does FST span a subset of given terminals? */
			all_spanned = TRUE;
			vp1 = cip -> edge [lfs];
			vp2 = cip -> edge [lfs + 1];
			while (vp1 < vp2) {
				t = *vp1++;
				if (NOT BITON (ltmask, t)) {
					all_spanned = FALSE;
					break;
				}
			}
			if (all_spanned) {
				SETBIT (lfmask, lfs);
				lflist [lnedges++] = lfs;
			}
		}
	}

	/* Reset terminal and FST masks */
	for (i = 0; i < ltlist -> n; i++) {
		CLRBIT (ltmask, ltlist -> a[i].pnum);
	}
	for (i = 0; i < lnedges; i++) {
		CLRBIT (lfmask, lflist [i]);
	}

	/* If no large FSTs then return (test cannot be performed) */
	if (lnedges EQ 0) {
		free ((char *) bsdmst);
		return TRUE;
	}

	/* Test if there is a PAIR of FSTs that has smaller or equal total length.   */
	/* In case the length is equal we only need to keep the "canonical"	     */
	/* par, i.e., for which the maximum index of (large) FSTs is minimized.	     */

	if (fst2 >= 0) {
		fstmaxind = -1;
		if ((cip -> edge_size [fst1] NE 2) AND (fst1 > fstmaxind)) fstmaxind = fst1;
		if ((cip -> edge_size [fst2] NE 2) AND (fst2 > fstmaxind)) fstmaxind = fst2;
	}
	else {
		fstmaxind = fst1;
	}

	nverts = ltlist -> n;
	shorter_pair_found = FALSE;
	for (i = 0; i < lnedges; i++) {
		lfs = lflist [i];

		/* Mark terminals in this FST */
		vp1 = cip -> edge [lfs];
		vp2 = cip -> edge [lfs + 1];
		while (vp1 < vp2) {
			t = *vp1++;
			SETBIT (ltmask, t);
		}

		if (cip -> edge_size [lfs] + 1 EQ nverts) {

			/* Find shortest BSD-MST edge to remaining terminal */
			ep   = bsdmst;
			minl = INF_DISTANCE;
			for (j = 0; j < nverts - 1; j++, ep++) {
				if ((NOT BITON (ltmask, ltlist -> a[ep -> p1].pnum)) OR
				    (NOT BITON (ltmask, ltlist -> a[ep -> p2].pnum))) {
					if (ep -> len < minl) minl = ep -> len;
				}
			}
			pairl = cip -> cost [lfs] + minl;
			if  (pairl <  l)			shorter_pair_found = TRUE;
			if ((pairl <= l) AND (lfs < fstmaxind)) shorter_pair_found = TRUE;
		}
		else {
			/* Try to combine with another FST */
			for (j = i+1; j < lnedges; j++) {
				lfs2  = lflist [j];
				if (cip -> edge_size [lfs] + cip -> edge_size [lfs2] - 1 NE nverts) continue;
				pairl = cip -> cost [lfs] + cip -> cost [lfs2];
				if (pairl > l) continue;

				/* Count intersecting terminals */
				isterms = 0;
				vp1 = cip -> edge [lfs2];
				vp2 = cip -> edge [lfs2 + 1];
				while (vp1 < vp2) {
					t = *vp1++;
					if (BITON (ltmask, t)) isterms++;
				}
				if (isterms NE 1) continue; /* not spanning all terminals */

				if  (pairl <  l)	    shorter_pair_found = TRUE;
				if ((pairl <= l)       AND
				    (lfs  < fstmaxind) AND
				    (lfs2 < fstmaxind))	    shorter_pair_found = TRUE;
			}
		}

		/* Reset terminal mask */
		for (j = 0; j < nverts; j++) {
			CLRBIT (ltmask, ltlist -> a[j].pnum);
		}
		if (shorter_pair_found) break;
	}

	if (shorter_pair_found) {
		free ((char *) bsdmst);
		return FALSE;
	}

	/* Allocate and initialize local cinfo structure */
	nedges =  lnedges + ltlist -> n - 1; /* we add BSD-MST edges */
	kmasks =  BMAP_ELTS (nverts);
	nmasks =  BMAP_ELTS (nedges);

	memset (&lcinfo, 0, sizeof (lcinfo));
	lcinfo.num_verts		= nverts;
	lcinfo.num_edges		= nedges;
	lcinfo.num_vert_masks		= kmasks;
	lcinfo.initial_vert_mask	= NEWA (kmasks, bitmap_t);
	lcinfo.num_edge_masks		= nmasks;
	lcinfo.initial_edge_mask	= NEWA (nmasks, bitmap_t);
	lcinfo.required_edges		= NEWA (nmasks, bitmap_t);
	lcinfo.edge			= NEWA (nedges + 1, int *);
	lcinfo.edge[0]			= NEWA (cip -> edge[cip -> num_edges] -
						cip -> edge[0], int);
	lcinfo.edge_size		= NEWA (nedges, int);
	lcinfo.cost			= NEWA (nedges, dist_t);
	lcinfo.tflag			= NEWA (nverts, bool);
	lcinfo.metric			= PURE_GRAPH;
	lcinfo.scale.scale		= 0;
	lcinfo.scale.scale_mul		= 1.0;
	lcinfo.scale.scale_div		= 1.0;
	lcinfo.integrality_delta	= 0;
	lcinfo.mst_length		= 0;
	lcinfo.description		= gst_strdup ("Local MSTHG");
	lterm				= NEWA (cip -> num_verts, int);

	for (i = 0; i < kmasks; i++) {
		lcinfo.initial_vert_mask [i] = 0;
	}
	for (i = 0; i < nverts; i++) {
		SETBIT (lcinfo.initial_vert_mask, i);
		lcinfo.tflag [i] = TRUE; /* we only meed to solve MSTHG's */
		lterm [ ltlist -> a[i].pnum ] = i;
	}
	for (i = 0; i < nmasks; i++) {
		lcinfo.initial_edge_mask [i] = 0;
		lcinfo.required_edges [i] = 0;
	}
	for (i = 0; i < nedges; i++) {
		SETBIT (lcinfo.initial_edge_mask, i);
	}

	hid   = 0;
	hep   = lcinfo.edge [0];

	/* Add all identified FSTs to hyperedge list */
	for (i = 0; i < lnedges; i++) {
		lfs = lflist [i];

		/* Make it more attractive to seek an MST with many FSTs */
		lcinfo.cost [hid]      = cip -> cost [lfs] -
					 cip -> integrality_delta / (double) nverts;
		lcinfo.edge_size [hid] = cip -> edge_size [lfs];
		lcinfo.edge [hid]      = hep;
		vp1 = cip -> edge [lfs];
		vp2 = cip -> edge [lfs + 1];
		while (vp1 < vp2) {
			*hep++ = lterm[ *vp1++ ];
		}
		hid++;
	}

	/* Add all BSD-MST edges to hyperedge list */
	ep = bsdmst;
	for (i = 0; i < ltlist -> n - 1; i++, ep++) {
		lcinfo.cost [hid]      = ep -> len -
					 cip -> integrality_delta / (double) nverts;
		lcinfo.edge_size [hid] = 2;
		lcinfo.edge [hid]      = hep;
		*(hep++) = ep -> p1;
		*(hep++) = ep -> p2;
		hid++;
	}
	lcinfo.edge [hid] = hep;

	/* Call branch-and-cut procedure */

	tracef_control.disabled = TRUE;
	init_term_trees (&lcinfo);
	lcinfo.inc_edges = compute_basic_incompat (&lcinfo);

	if ((lcinfo.num_verts <= 24) AND (lcinfo.num_edges <= 32)) {

		/* Small problem -- use backtrack search. */

		smt_mask = 0;
		msthgl = backtrack_search (&lcinfo, &smt_mask);
		lsmtcount = 0;
		for (i = 0; i < lcinfo.num_edges; i++) {
			if ((smt_mask & (1 << i)) NE 0) {
				++lsmtcount;
			}
		}
	}
	else {

		bbip = create_bbinfo (&lcinfo);

#if 0
		/* Print problem */
		{
		int * ep1;
		int * ep2;
		printf("\n%d %d\n", lcinfo.num_verts, lcinfo.num_edges);
		for (i = 0; i < lcinfo.num_edges; i++) {
			ep1 = lcinfo.edge[i];
			ep2 = lcinfo.edge[i+1];
			while (ep1 < ep2)
				printf("%d ", *(ep1++) + 1);
			printf("%f\n", lcinfo.cost[i]);
		}
		printf("%d\n", lcinfo.num_verts);
		for (i = 0; i < lcinfo.num_verts; i++)
			printf("%d\n", i + 1);
		}
#endif

		/* Find MST in hypergraph. */

		msthgl = branch_and_cut (bbip);

		/* Count number of FSTs in SMT */
		lsmtcount = 0;
		for (i = 0; i < lcinfo.num_edges; i++) {
			if (BITON (bbip->smt, i)) lsmtcount++;
		}

		/* Free allocated memory */
		destroy_bbinfo (bbip);

	}

	free_phase_1_data (&lcinfo);
	free ((char *) lterm);
	free ((char *) bsdmst);

	/* Return result of final test */
	if (fst2 >= 0) {
		if (msthgl < l - cip -> integrality_delta) {
			return FALSE;
		}
		/* If equal length then there should be at least three FSTs */
		else if ((msthgl < l) AND (lsmtcount >= 3)) {
			return FALSE;
		}
	}
	else {
		if (msthgl <= l) {
			return FALSE;
		}
	}

	return TRUE; /* all tests passed */

}

/*
 * Find distance and attachment terminals for given terminal
 * to a specific FST
 */

	static
	void
test_close_terminal (

struct pinfo *		pip,	/* IN	  - pruning data structure  */
int			fst,	/* IN	  - FST index */
int			term,	/* IN	  - Terminal */
struct clt_info**	clip	/* IN/OUT - store close terminal info here */
)
{
int			j;
int			i1;
int			i2;
dist_t			d;
dist_t			d1;
dist_t			d2;
dist_t			pg_longest;
struct point *		pt;
struct point *		p1;
struct point *		p2;
struct point		clp;
struct full_set *	fsp;
struct pset *		fsp_terms;
struct pset *		fsp_steins;
struct cinfo *		cip;
struct clt_info		clt;

	cip	   = pip -> cip;
	pt	   = &(cip -> pts -> a [term]);
	fsp	   = cip -> full_trees [fst];
	fsp_terms  = cip -> full_trees [fst] -> terminals;
	fsp_steins = cip -> full_trees [fst] -> steiners;
	pg_longest = pip -> pg_edges [pip -> num_pg_edges-1].len;

	/* First a rough test to eliminate the terminal */
	if (EDIST(&(fsp_terms -> a[0]), pt) >
	    (fsp_terms -> n - 1) * pg_longest)
		return; /* this terminal is too far away */

	/* Now find the edge that is closest edge to this terminal */
	memset (&clt, 0, sizeof (clt));
	clt.term = term;
	clt.dist = INF_DISTANCE;
	for (j = 0; j < fsp -> nedges; j++) {

		/* Get coordinates and pruning graph indices for endpoints */
		if (fsp -> edges[j].p1 < fsp_terms -> n) {
			p1 = &(fsp_terms  -> a[ fsp -> edges[j].p1 ]);
			i1 = fsp_terms -> a[ fsp -> edges[j].p1 ].pnum;
		}
		else {
			p1 = &(fsp_steins -> a[ fsp -> edges[j].p1 -
						fsp_terms -> n ]);
			i1 = pip -> steiner_index [fst] +
				fsp -> edges[j].p1 - fsp_terms -> n;
		}

		if (fsp -> edges[j].p2 < fsp_terms -> n) {
			p2 = &(fsp_terms  -> a[ fsp -> edges[j].p2 ]);
			i2 = fsp_terms -> a[ fsp -> edges[j].p2 ].pnum;
		}
		else {
			p2 = &(fsp_steins -> a[ fsp -> edges[j].p2 -
						fsp_terms -> n ]);
			i2 = pip -> steiner_index [fst] +
				fsp -> edges[j].p2 - fsp_terms -> n;
		}

		/* Compute closest distance from terminal to edge */
		d = terminal_edge_distance(cip, pt, p1, p2, &clp, &d1, &d2) *
		    (1.0 + pip -> eps * ((double) cip -> edge_size [fst]));
		if (d < clt.dist) {
			clt.dist = d;
			clt.aterm1 = (d1 <= d) ? i1 : pip -> num_pg_verts;
			clt.aterm2 = (d2 <= d) ? i2 : pip -> num_pg_verts;
		}
	}

	/* Is distance smaller than longest edge in pruning graph? */
	if (clt.dist < pg_longest) {
		*((*clip)++) = clt; /* save this terminal */
	}
}

/*
 * Distance from point to edge
 */

	static
	dist_t
terminal_edge_distance (

struct cinfo *		cip,	/* IN  - compatibility info */
struct point *		pt,	/* IN  - point */
struct point *		p1,	/* IN  - first edge point */
struct point *		p2,	/* IN  - second edge point */
struct point *		clp,	/* OUT - closest point */
dist_t *		d1,	/* OUT - distance from clp to p1 */
dist_t *		d2	/* OUT - distance to clp to p2 */
)
{
dist_t			d;
dist_t			l1;
dist_t			l2;
dist_t			l;
struct point *		a;
struct point *		b;
struct point		e;
struct point		ctr;

	d = 0.0;
	if (cip -> metric EQ RECTILINEAR) {

		/* Compute distance to rectangle given by p1 and p2.
		   Also identify (one) closest point. */

		l1 = pt -> x - p1 -> x;
		l2 = pt -> x - p2 -> x;
		clp -> x = pt -> x;
		if ((l1 < 0.0) AND (l2 < 0.0)) {
			clp -> x = MIN(p1 -> x, p2 -> x);
		}
		if ((l1 > 0.0) AND (l2 > 0.0)) {
			clp -> x = MAX(p1 -> x, p2 -> x);
		}

		l1 = pt -> y - p1 -> y;
		l2 = pt -> y - p2 -> y;
		clp -> y = pt -> y;
		if ((l1 < 0.0) AND (l2 < 0.0)) {
			clp -> y = MIN(p1 -> y, p2 -> y);
		}
		if ((l1 > 0.0) AND (l2 > 0.0)) {
			clp -> y = MAX(p1 -> y, p2 -> y);
		}
		d   = RDIST(pt, clp);
		*d1 = RDIST(p1, clp);
		*d2 = RDIST(p2, clp);

	}

	if (cip -> metric EQ EUCLIDEAN) {

		/* Compute Steiner point for the three points and then
		   the distance from pt to this Steiner point */

		/* Make sure points are oriented nicely */
		if (left_turn(p1, p2, pt)) {
			a = p1; b = p2;
		}
		else {
			a = p2; b = p1;
		}

		eq_point (a, b, &e);
		eq_circle_center (a, b, &e, &ctr);

		if (right_turn (&e, a, pt)) {
			if (left_turn (&e, b, pt)) {
				if (sqr_dist(&ctr, pt) > sqr_dist(&ctr, &e)) {
					project_point(&e, &ctr, pt, clp);
					l = EDIST(&e, pt);
				}
				else {
					*clp = *pt;
					l = EDIST(a, pt) + EDIST(b, pt);
				}
			}
			else {
				*clp = *b;
				l = EDIST(a, b) + EDIST(pt, b);
			}
		}
		else {
			*clp = *a;
			l = EDIST(b, a) + EDIST(pt, a);
		}

		d   = l - EDIST(p1, p2);
		*d1 = EDIST(p1, clp);
		*d2 = EDIST(p2, clp);
	}


	return d;
}

/*
 * Find all FSTs whose removal splits the problem.  Make these
 * all required.
 */

	static
	int
bcc_find_required (

struct cinfo *	cip,		/* IN - compatibility info */
int *		made_req,	/* IN/OUT - list of FSTs made REQUIRED */
int		req_count	/* IN - existing required count */
)
{
int			i;
int			j;
int			nverts;
int			nedges;
int			kmasks;
int			nmasks;
int			max_stack;
struct bc3		bc;

	nverts	= cip -> num_verts;
	nedges	= cip -> num_edges;

	if (nverts < 2) {
		/* No BCC's. */
		return (req_count);
	}

	kmasks	= BMAP_ELTS (nverts);
	nmasks	= BMAP_ELTS (nedges);

	bc.cip		= cip;
	bc.kmasks	= kmasks;
	bc.nmasks	= nmasks;
	bc.dfs		= NEWA (nverts, int);
	bc.low		= NEWA (nverts, int);
	bc.parent	= NEWA (nverts, int);

	for (i = 0; i < nverts; i++) {
		bc.dfs	[i] = 0;
		bc.low [i] = 0;
		bc.parent [i] = -1;
	}

	j = (nedges > nverts) ? nedges : nverts;
	bc.max_stack	= j;
	bc.stack	= NEWA (j, int);
	bc.sp		= bc.stack;
	bc.counter	= 0;

	bc.bcc_vmask	= NEWA (kmasks, bitmap_t);
	bc.bcc_emask	= NEWA (nmasks, bitmap_t);
	bc.edges_seen	= NEWA (nmasks, bitmap_t);
	bc.bcc_vlist	= NEWA (nverts, int);
	bc.degree	= NEWA (nverts, int);
	bc.made_req	= made_req;
	bc.req_count	= req_count;

	for (i = 0; i < bc.kmasks; i++) {
		bc.bcc_vmask [i] = 0;
	}
	for (i = 0; i < bc.nmasks; i++) {
		bc.bcc_emask [i] = 0;
		bc.edges_seen [i] = 0;
	}

	/* Traverse each connected component, identifying its BCC's as	*/
	/* we go.							*/
	for (i = 0; i < nverts; i++) {
		if (NOT BITON (cip -> initial_vert_mask, i)) continue;
		if (bc.dfs [i] > 0) continue;

		/* Traverse one connected component, finding	*/
		/* each of its BCCs as we go.			*/

		bcc3 (&bc, i);
	}

	free ((char *) bc.degree);
	free ((char *) bc.bcc_vlist);
	free ((char *) bc.edges_seen);
	free ((char *) bc.bcc_emask);
	free ((char *) bc.bcc_vmask);
	free ((char *) bc.stack);
	free ((char *) bc.parent);
	free ((char *) bc.low);
	free ((char *) bc.dfs);

	return (bc.req_count);
}

/*
 * This is the recursive part of the bi-connected-components algorithm.	 It
 * is the standard method, with a few tweaks to work on hypergraphs instead.
 * We process each bi-connected component individually.
 */

	static
	void
bcc3 (

struct bc3 *		bcp,		/* IN - global BCC data */
int			v		/* IN - current DFS vertex */
)
{
int			i;
int			e;
int			e2;
int			w;
int *			ep1;
int *			ep2;
int *			vp1;
int *			vp2;
int *			vp3;
int *			vp4;
int *			stack_endp;
int *			sp;
int *			stack;
struct cinfo *		cip;

	cip = bcp -> cip;

	if ((v < 0) OR (v >= cip -> num_verts)) {
		fatal ("bcc3:  Bug 1.");
	}

	++(bcp -> counter);
	bcp -> dfs [v] = bcp -> counter;
	bcp -> low [v] = bcp -> counter;
	ep1 = cip -> term_trees [v];
	ep2 = cip -> term_trees [v + 1];
	while (ep1 < ep2) {
		e = *ep1++;
		if ((e < 0) OR (e >= cip -> num_edges)) {
			fatal ("bcc3: Bug 2.");
		}
		if (NOT BITON (cip -> initial_edge_mask, e)) continue;
		if (NOT BITON (bcp -> edges_seen, e)) {
			/* We haven't seen this edge before.  Push	*/
			/* it onto the stack...				*/
			stack_endp = &(bcp -> stack [bcp -> max_stack]);
			if ((bcp -> sp < bcp -> stack) OR
			    (bcp -> sp >= stack_endp)) {
				fatal ("bcc3: Bug 3.");
			}
			*(bcp -> sp)++ = e;
			SETBIT (bcp -> edges_seen, e);
		}
		/* Scan the vertices and process them... */
		vp1 = cip -> edge [e];
		vp2 = cip -> edge [e + 1];
		while (vp1 < vp2) {
			w = *vp1++;
			if ((w < 0) OR (w >= cip -> num_verts)) {
				fatal ("bcc3: Bug 4.");
			}
			if (bcp -> dfs [w] EQ 0) {
				bcp -> parent [w] = v;
				bcc3 (bcp, w);
				if (bcp -> low [w] >= bcp -> dfs [v]) {
					/* We have a new BCC! */
					stack	= bcp -> stack;
					sp	= bcp -> sp;
					do {
						if (sp <= stack) {
							fatal ("bcc3: Bug 5.");
						}
						e2 = *--sp;
					} while (e2 NE e);

					/* Process the bi-connected comp. */
					process_bcc3 (bcp, sp, bcp -> sp);

					/* Pop BCC edges from stack */
					bcp -> sp = sp;
				}
				if (bcp -> low [w] < bcp -> low [v]) {
					bcp -> low [v] = bcp -> low [w];
				}
			}
			else if ((w NE bcp -> parent [v]) AND
				 (bcp -> dfs [w] < bcp -> low [v])) {
				bcp -> low [v] = bcp -> dfs [w];
			}
		}
	}
}

/*
 * Process a single bi-connected component, specified as a list of edges.
 * We look for vertices having degree 1 within the component.  When this
 * happens, the incident edge is required.
 */

	static
	void
process_bcc3 (

struct bc3 *	bcp,		/* IN - global BCC data */
int *		edge_ptr,	/* IN - list of BCC edges */
int *		endp		/* IN - end of BCC edge list */
)
{
int		i;
int		e;
int		v;
int		n;
int *		ep1;
int *		ep2;
int *		ep3;
int *		vp1;
int *		vp2;
int *		vlp;
int *		elist;
struct cslist * cslp;
struct cinfo *	cip;
double		xe;

	cip = bcp -> cip;

	/* Gather a list of all vertices in this BCC.  Compute	*/
	/* their degrees (with respect to the BCC).		*/
	vlp = bcp -> bcc_vlist;

	ep1 = edge_ptr;
	ep2 = endp;
	while (ep1 < ep2) {
		e = *ep1++;
		SETBIT (bcp -> bcc_emask, e);
		vp1 = cip -> edge [e];
		vp2 = cip -> edge [e + 1];
		while (vp1 < vp2) {
			v = *vp1++;
			if (NOT BITON (cip -> initial_vert_mask, v)) continue;
			if (NOT BITON (bcp -> bcc_vmask, v)) {
				*vlp++ = v;
				SETBIT (bcp -> bcc_vmask, v);
				bcp -> degree [v] = 0;
			}
			++(bcp -> degree [v]);
		}
	}

	/* All of the vertices of this BCC are now known, as are their	*/
	/* degrees (relative to the component).	 Now look for vertices	*/
	/* of degree 1.							*/

	vp1 = bcp -> bcc_vlist;
	vp2 = vlp;
	while (vp1 < vp2) {
		v = *vp1++;
		CLRBIT (bcp -> bcc_vmask, v);		/* clean up as we go */

		n = bcp -> degree [v];
		if (n > 1) continue;

		ep1 = cip -> term_trees [v];
		ep2 = cip -> term_trees [v + 1];
		for (;;) {
			if (ep1 >= ep2) {
				fatal ("process_bcc3: Bug 1.");
			}
			e = *ep1++;
			if (BITON (bcp -> bcc_emask, e)) break;
		}
		if (BITON (cip -> required_edges, e)) continue;
		i = (bcp -> req_count)++;
		bcp -> made_req [i] = e;
		SETBIT (cip -> required_edges, e);
	}

	/* Clean up edge mark flags. */
	ep1 = edge_ptr;
	ep2 = endp;
	while (ep1 < ep2) {
		e = *ep1++;
		CLRBIT (bcp -> bcc_emask, e);
	}
}

/*
 * For sorting integers in place
 */

	static
	int
comp_ints (

const void *		p1,
const void *		p2
)
{
int			l1;
int			l2;

	l1 = *((int *) p1);
	l2 = *((int *) p2);

	if (l1 < l2) return (-1);
	if (l1 > l2) return (1);
	return (0);
}


/*
 * For sorting pruning graph in place
 */

	static
	int
comp_pg_edges (

const void *		p1,
const void *		p2
)
{
dist_t			l1;
dist_t			l2;

	l1 = ((struct pg_edge *) p1) -> len;
	l2 = ((struct pg_edge *) p2) -> len;

	if (l1 < l2) return (-1);
	if (l1 > l2) return (1);
	return (0);
}

/*
 * For sorting close terminals in place
 */

	static
	int
comp_clt (

const void *		p1,
const void *		p2
)
{
dist_t			l1;
dist_t			l2;

	l1 = ((struct clt_info *) p1) -> dist;
	l2 = ((struct clt_info *) p2) -> dist;

	if (l1 < l2) return (-1);
	if (l1 > l2) return (1);
	return (0);
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
