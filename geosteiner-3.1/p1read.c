/***********************************************************************

	File:	p1read.c
	Rev:	b-2
	Date:	02/28/2001

	Copyright (c) 1993, 2001 by David M. Warme

************************************************************************

	Routines for reading the data that is output from phase 1.

************************************************************************

	Modification Log:

	b-1:	01/12/97	warme
		: Split p1io.c into reading and writing parts.
	b-2:	02/28/2001	warme
		: Numerous changes for 3.1 release.
		: New input scaling stuff.
		: Fixed uninitialized incompat array bug.
		: Use common sort_ints routine.
		: Free *all* cinfo data structures, and be more
		:  careful about freeing term_trees.

************************************************************************/

#include "p1io.h"
#include "steiner.h"


/*
 * Global Routines
 */

int **		compute_basic_incompat (struct cinfo *);
void		free_phase_1_data (struct cinfo *);
void		init_term_trees (struct cinfo *);
void		read_phase_1_data (struct cinfo *);


/*
 * External References
 */

	/* none */


/*
 * Local Routines
 */

static void		double_to_hex (double, char *);
static int		get_d (void);
static double		get_dec_double (void);
static double		get_hex_double (void);
static void		get_line (char *, int);
static int		get_version (void);
static bitmap_t		get_x (void);
static double		hex_to_double (char *);
static int		hexdig (char);
static void		init_inc_edges (struct cinfo *, int **, int *);
static void		read_duplicate_terminal_groups (struct cinfo *, int);
static void		read_version_0 (struct cinfo *, int);
static void		read_version_2 (struct cinfo *, int);
static void		remove_duplicates (int, int **, bitmap_t *);
static void		skip (void);
static void		verify_symmetric (int **, int);

/*
 * This routine reads in the data as printed by the print-phase-1-data
 * routine.
 */

	void
read_phase_1_data (

struct cinfo *	cip		/* OUT - compatibility info. */
)
{
int		version;

	/* Zero everything initially... */
	(void) memset (cip, 0, sizeof (*cip));

	version = get_version ();

	switch (version) {
	case P1IO_VERSION_0:
		read_version_0 (cip, version);
		break;

	case P1IO_VERSION_2:
	case P1IO_VERSION_3:
		read_version_2 (cip, version);
		break;

	default:
		fatal ("read_phase_1_data: Bug 1.");
		break;
	}

	/* Set up the output conversion stuff. */
	init_output_conversion (cip -> pts, cip -> metric, &(cip -> scale));
}

/*
 * This routine reads in the extended OR-library format.
 */

	static
	void
read_version_0 (

struct cinfo *	cip,		/* OUT - compatibility info. */
int		version		/* IN - version number. */
)
{
int			i;
int			j;
int			k;
char			c;
int			nverts;
int			nedges;
int			kmasks;
int			nmasks;
int			scale_factor;
int			total_edge_cardinality;
int			nt;
char *			s;
struct numlist *	line;
struct numlist *	p;
struct numlist **	hookp;
struct numlist **	prev_hookp;
struct numlist *	costs;
struct numlist **	cost_hookp;
int *			ip1;
int *			ip2;

	nverts = get_d ();
	nedges = get_d ();

	kmasks = BMAP_ELTS (nverts);
	nmasks = BMAP_ELTS (nedges);

	cip -> num_verts	= nverts;
	cip -> num_edges	= nedges;
	cip -> num_vert_masks	= kmasks;
	cip -> num_edge_masks	= nmasks;

	cip -> edge		= NEWA (nedges + 1, int *);
	cip -> edge_size	= NEWA (nedges, int);
	cip -> cost		= NEWA (nedges, dist_t);
	cip -> tflag		= NEWA (nverts, bool);
	cip -> metric		= PURE_GRAPH;

	cip -> initial_vert_mask = NEWA (kmasks, bitmap_t);
	cip -> initial_edge_mask = NEWA (nmasks, bitmap_t);
	cip -> required_edges	 = NEWA (nmasks, bitmap_t);

	memset (cip -> initial_vert_mask, 0, kmasks * sizeof (bitmap_t));
	memset (cip -> initial_edge_mask, 0, nmasks * sizeof (bitmap_t));
	memset (cip -> required_edges,	  0, nmasks * sizeof (bitmap_t));

	for (i = 0; i < nverts; i++) {
		SETBIT (cip -> initial_vert_mask, i);
	}
	for (i = 0; i < nedges; i++) {
		SETBIT (cip -> initial_edge_mask, i);
	}

	costs = NULL;
	cost_hookp = &costs;

	total_edge_cardinality = 0;

	for (i = 0; i < nedges; i++) {
		line = parse_line_of_numbers (stdin);
		if (line EQ NULL) {
			fprintf (stderr, "Unexpected EOF!\n");
			exit (1);
		}

		/* Count the numbers and detach the last one */
		j = 0;
		prev_hookp = NULL;
		hookp = &line;
		for (;;) {
			p = *hookp;
			if (p EQ NULL) break;
			prev_hookp = hookp;
			hookp = &(p -> next);
			++j;
		}
		if (j < 3) {
			fprintf (stderr,
				 "Hyperedge with less than 2 vertices!\n");
			exit (1);
		}

		/* Transfer last number onto cost list. */
		p = *prev_hookp;
		*prev_hookp = NULL;
		--j;
		cip -> edge_size [i] = j;
		total_edge_cardinality += j;

		*cost_hookp = p;
		cost_hookp = &(p -> next);

		ip1 = NEWA (j, int);
		cip -> edge [i] = ip1;

		/* Verify remaining numbers are integers and convert. */
		j = 0;
		for (p = line; p NE NULL; p = p -> next) {
			if (p -> expon < 0) goto nonint;
			if (p -> expon > 0) {
				do {
					p -> mantissa *= 10.0;
				} while (--(p -> expon) > 0);
			}
			if ((p -> mantissa < 1) OR
			    (nverts < p -> mantissa) OR
			    (p -> mantissa NE floor (p -> mantissa))) {
				fprintf (stderr,
					 "Vertex number %g out of range.\n",
					 p -> mantissa);
				exit (1);
			}
			*ip1++ = p -> mantissa - 1;
		}
		if (cip -> edge [i] + cip -> edge_size [i] NE ip1) {
			fatal ("read_version_0: Bug 1.");
		}
	}

	/* Now convert the hyperedge costs. */
	scale_factor = compute_scaling_factor (costs);

	set_scale_info (&(cip -> scale), scale_factor);

	i = 0;
	while (costs NE NULL) {
		p = costs;
		costs = p -> next;
		p -> next = NULL;
		cip -> cost [i++] = read_numlist (p, &(cip -> scale));
	}
	if (i NE nedges) {
		fatal ("read_version_0: Bug 2.");
	}

	/* Read in the list of terminal vertices... */
	if (scanf (" %d", &nt) NE 1) {
		/* EOF -- assume all vertices are terminals. */
		for (i = 0; i < nverts; i++) {
			cip -> tflag [i] = TRUE;
		}
	}
	else {
		/* All are Steiner vertices unless listed! */
		for (i = 0; i < nverts; i++) {
			cip -> tflag [i] = FALSE;
		}
		for (i = 0; i < nt; i++) {
			j = get_d ();
			if ((j < 1) OR (nverts < j)) {
				fatal ("read_version_0: Bug 3.");
			}
			cip -> tflag [j - 1] = TRUE;
		}
	}

	/* Copy hyperedges into contiguous form. */
	ip1 = NEWA (total_edge_cardinality, int);
	for (i = 0; i < nedges; i++) {
		ip2 = cip -> edge [i];
		cip -> edge [i] = ip1;
		k = cip -> edge_size [i];
		for (j = 0; j < k; j++) {
			*ip1++ = ip2 [j];
		}
		free ((char *) ip2);
	}
	cip -> edge [i] = ip1;

	init_term_trees (cip);
	cip -> inc_edges = compute_basic_incompat (cip);
	return;

	/* Error exit. */
nonint:
	fprintf (stderr, "Expected integer!\n");
	exit (1);
}

/*
 * This routine reads in the new data formats -- versions 2 and 3.
 */

	static
	void
read_version_2 (

struct cinfo *	cip,		/* OUT - compatibility info. */
int		version		/* IN - version number. */
)
{
int			i;
int			j;
int			k;
int			nverts;
int			nedges;
int			kmasks;
int			nmasks;
int			nt;
int			ns;
int			scale_factor;
int			fst_edge_count;
int			num_incompat;
int			ndg;
int			count;
int			total_edge_card;
struct pset *		terms;
struct pset *		steins;
struct point *		p1;
struct full_set *	fsp;
struct edge *		ep;
bitmap_t *		bp1;
bitmap_t *		bp2;
bitmap_t *		bp3;
bitmap_t		mask;
int			bitcount;
int *			ip1;
int *			ip2;
int *			icounts;
int **			incompat;
bool			geometric;
dist_t			len;
char			buf1 [128];

	skip ();	/* Skip rest of version line */

	get_line (buf1, sizeof (buf1));
	i = strlen (buf1);
	if (i > 79) {
		fprintf (stderr, "Description field too long.\n");
		exit (1);
	}
	cip -> description = new (i + 1);
	strcpy (cip -> description, buf1);

	cip -> metric = get_d ();
	if ((cip -> metric NE RECTILINEAR) AND
	    (cip -> metric NE EUCLIDEAN) AND
	    (cip -> metric NE PURE_GRAPH)) {
		fprintf (stderr, "Bad metric: %d\n", cip -> metric);
		exit (1);
	}
	geometric = (cip -> metric NE PURE_GRAPH);

	nverts = (int) get_d ();

	if (geometric) {
		(void) get_dec_double ();
		cip -> mst_length = get_hex_double ();
	}
	else {
		cip -> mst_length = 0;
	}

	ndg = 0;
	if ((version <= P1IO_VERSION_2) AND geometric) {
		/* Version 2 has duplicate terminal groups... */
		ndg = (int) get_d ();
	}

	scale_factor = (int) get_d ();

	set_scale_info (&(cip -> scale), scale_factor);

	cip -> integrality_delta = 0;
	if (version >= P1IO_VERSION_3) {
		/* Version 3 has integrality delta... */
		(void) get_dec_double ();
		cip -> integrality_delta = get_hex_double ();
	}

	skip ();	/* skip rest of line before machine description */
	get_line (buf1, sizeof (buf1));	/* read machine description line */

	cip -> p1time = get_d ();

	nedges = (int) get_d ();

	kmasks = BMAP_ELTS (nverts);
	nmasks = BMAP_ELTS (nedges);

	cip -> initial_vert_mask = NEWA (kmasks, bitmap_t);
	for (i = 0; i < kmasks; i++) {
		cip -> initial_vert_mask [i] = 0;
	}
	for (i = 0; i < nverts; i++) {
		SETBIT (cip -> initial_vert_mask, i);
	}
	cip -> initial_edge_mask	= NEWA (nmasks, bitmap_t);
	cip -> required_edges		= NEWA (nmasks, bitmap_t);
	for (i = 0; i < nmasks; i++) {
		cip -> initial_edge_mask [i]	= 0;
		cip -> required_edges [i]	= 0;
	}

	cip -> num_edges	= nedges;
	cip -> num_verts	= nverts;
	cip -> num_vert_masks	= kmasks;
	cip -> num_edge_masks	= nmasks;
	cip -> edge		= NEWA (nedges + 1, int *);
	cip -> edge_size	= NEWA (nedges, int);
	cip -> cost		= NEWA (nedges, dist_t);
	cip -> tflag		= NEWA (nverts, bool);
	if (geometric) {
		cip -> pts	= NEW_PSET (nverts);
		ZERO_PSET (cip -> pts, nverts);
		cip -> pts -> n = nverts;
		cip -> full_trees = NEWA (nedges, struct full_set *);
	}

	incompat		= NEWA (nedges, int *);
	icounts			= NEWA (nedges, int);

	if (geometric) {
		/* Read in the terminals... */
		for (i = 0; i < nverts; i++) {
			p1 = &(cip -> pts -> a [i]);
			(void) get_dec_double ();
			(void) get_dec_double ();
			p1 -> x		 = get_hex_double ();
			p1 -> y		 = get_hex_double ();
			p1 -> pnum	 = i;
		}
	}

	if (version >= P1IO_VERSION_3) {
		/* Version 3 has terminal/Steiner flag for each vertex. */
		for (i = 0; i < nverts; i++) {
			cip -> tflag [i] = get_d ();
		}
	}
	else {
		/* Version 2: assume all vertices are terminals. */
		for (i = 0; i < nverts; i++) {
			cip -> tflag [i] = TRUE;
		}
	}

	if ((version <= P1IO_VERSION_2) AND geometric) {
		read_duplicate_terminal_groups (cip, ndg);
	}

	/* hyperedges... */
	total_edge_card = 0;
	for (i = 0; i < nedges; i++) {
		nt = get_d ();
		total_edge_card += nt;
		cip -> edge_size [i] = nt;
		ip1 = NEWA (nt, int);
		cip -> edge [i] = ip1;
		if (geometric) {
			fsp = NEW (struct full_set);
			(void) memset (fsp, 0, sizeof (*fsp));
			cip -> full_trees [i] = fsp;
			fsp -> tree_num = i;
			terms = NEW_PSET (nt);
			ZERO_PSET (terms, nt);
			terms -> n = nt;
		}
		for (j = 0; j < nt; j++) {
			k = get_d ();
			if ((k < 1) OR (k > nverts)) {
				fprintf (stderr,
					 "Terminal index out of range.\n");
				exit (1);
			}
			--k;
			*ip1++ = k;
			if (geometric) {
				terms -> a [j] = cip -> pts -> a [k];
			}
		}
		(void) get_dec_double ();
		len = get_hex_double ();
		cip -> cost [i] = len;

		if (geometric) {
			fsp -> tree_len = len;
			ns = get_d ();
			steins = NEW_PSET (ns);
			ZERO_PSET (steins, ns);
			steins -> n = ns;
			for (j = 0; j < ns; j++) {
				p1 = &(steins -> a [j]);
				(void) get_dec_double ();
				(void) get_dec_double ();
				p1 -> x = get_hex_double ();
				p1 -> y = get_hex_double ();
				p1 -> pnum = j;
			}
			fsp -> terminals = terms;
			fsp -> steiners = steins;
			fst_edge_count = get_d ();
			fsp -> nedges = fst_edge_count;
			ep = NEWA (fst_edge_count, struct edge);
			fsp -> edges = ep;
			for (j = 0; j < fst_edge_count; j++, ep++) {
				ep -> len = 0;	/* should be unused... */
				k = get_d ();
				if ((-ns <= k) AND (k < 0)) {
					k = nt - k - 1;
				}
				else if ((0 < k) AND (k <= nt)) {
					--k;
				}
				else {
					fprintf (stderr, "Invalid edge endpoint!\n");
					exit (1);
				}
				ep -> p1 = k;
				k = get_d ();
				if ((-ns <= k) AND (k < 0)) {
					k = nt - k - 1;
				}
				else if ((0 < k) AND (k <= nt)) {
					--k;
				}
				else {
					fprintf (stderr, "Invalid edge endpoint!\n");
					exit (1);
				}
				ep -> p2 = k;
			}
		}
		k = get_d ();	/* full set status... */
		switch (k) {
		case 0:
			break;

		case 1:
			SETBIT (cip -> initial_edge_mask, i);
			break;

		case 2:
			SETBIT (cip -> required_edges, i);
			SETBIT (cip -> initial_edge_mask, i);
			break;

		default:
			fprintf (stderr, "Invalid full set status: %d\n", k);
			exit (1);
		}

		num_incompat = get_d ();
		icounts [i] = num_incompat;
		ip1 = NULL;
		if (num_incompat > 0) {
			ip1 = NEWA (num_incompat, int);
			for (j = 0; j < num_incompat; j++) {
				k = get_d ();
				if ((k <= 0) OR (k > nedges)) {
					fprintf (stderr, "Bad incompatible index.\n");
					exit (1);
				}
				ip1 [j] = k - 1;
			}
		}
		incompat [i] = ip1;

		if (version <= P1IO_VERSION_2) {
			/* Version 2 has strongly compatible full sets... */
			k = get_d ();
			for (j = 0; j < k; j++) {
				(void) get_d ();
			}
		}
	}

	/* Copy hyperedges into contiguous form. */
	ip1 = NEWA (total_edge_card, int);
	for (i = 0; i < nedges; i++) {
		ip2 = cip -> edge [i];
		cip -> edge [i] = ip1;
		k = cip -> edge_size [i];
		for (j = 0; j < k; j++) {
			*ip1++ = ip2 [j];
		}
		free ((char *) ip2);
	}
	cip -> edge [i] = ip1;

	init_term_trees (cip);
	init_inc_edges (cip, incompat, icounts);

	free ((char *) icounts);
	free ((char *) incompat);
}

/*
 * This routine reads in the duplicate terminal groups.  We do not know
 * how big this is going to be.  Therefore we must plan for the worst
 * and reallocate when we are done...
 */

	static
	void
read_duplicate_terminal_groups (

struct cinfo *		cip,		/* IN/OUT - compatibility info */
int			ndg		/* IN - number of duplicate groups */
)
{
int		i;
int		j;
int		k;
int		n;
int		t;
int		nverts;
int		kmasks;
int *		ip;
int *		real_terms;
bitmap_t *	mark;
int *		index;
int *		terms;
int **		dup_grps;

	if (ndg <= 0) {
		return;				/* nothing to read */
	}

	nverts = cip -> num_verts;

	if (2 * ndg > nverts) {
		/* Too many duplicate terminal groups! */
		fatal ("read_duplicate_terminal_groups: Bug 1.");
	}

	kmasks = cip -> num_vert_masks;

	mark = NEWA (kmasks, bitmap_t);
	index = NEWA (ndg + 1, int);
	terms = NEWA (nverts, int);

	for (i = 0; i < kmasks; i++) {
		mark [i] = 0;
	}

	k = 0;
	for (i = 0; i < ndg; i++) {
		index [i] = k;
		n = (int) get_d ();
		for (j = 0; j < n; j++) {
			t = get_d ();
			if ((t < 1) OR (cip -> num_verts < t)) {
				fatal ("read_duplicate_terminal_groups: Bug 2.");
			}
			--t;
			if (BITON (mark, t)) {
				fatal ("read_duplicate_terminal_groups: Bug 3.");
			}
			SETBIT (mark, t);
			terms [k++] = t;
		}
	}
	index [i] = k;

	real_terms = NEWA (k, int);
	(void) memcpy ((char *) real_terms,
		       (char *) terms,
		       k * sizeof (int));

	dup_grps = NEWA (ndg + 1, int *);
	for (i = 0; i <= ndg; i++) {
		dup_grps [i] = real_terms + index [i];
	}

	free ((char *) dup_grps);
	free ((char *) real_terms);
	free ((char *) terms);
	free ((char *) index);
	free ((char *) mark);

	/* Remove duplicate terminals from the problem... */
	remove_duplicates (ndg, dup_grps, cip -> initial_vert_mask);
}

/*
 * This routine creates the "term_trees" array that is indexed by point
 * number and gives a list of all tree-numbers involving that point.
 */

	void
init_term_trees (

struct cinfo *	cip		/* IN - compatibility info */
)
{
int			i;
int			j;
int			n;
int			nverts;
int			nedges;
struct full_set *	fsp;
struct pset *		terms;
struct point *		p1;
int			total;
int *			ip;
int *			counts;
int **			ptrs;
int *			vp1;
int *			vp2;

	nverts = cip -> num_verts;
	nedges = cip -> num_edges;

	counts	= NEWA (nverts, int);

	total = 0;
	for (i = 0; i < nverts; i++) {
		counts [i] = 0;
	}

	for (i = 0; i < nedges; i++) {
		vp1 = cip -> edge [i];
		vp2 = cip -> edge [i + 1];
		while (vp1 < vp2) {
			++counts [*vp1++];
			++total;
		}
	}

	ip			= NEWA (total, int);
	cip -> term_trees	= NEWA (nverts + 1, int *);
	ptrs			= NEWA (nverts, int *);

	for (i = 0; i < nverts; i++) {
		cip -> term_trees [i]	= ip;
		ptrs [i]		= ip;
		ip += counts [i];
	}
	cip -> term_trees [i] = ip;

	for (i = 0; i < nedges; i++) {
		vp1 = cip -> edge [i];
		vp2 = cip -> edge [i + 1];
		while (vp1 < vp2) {
			j = *vp1++;
			*(ptrs [j])++ = i;
		}
	}

	free ((char *) ptrs);
	free ((char *) counts);
}

/*
 * This routine merges the FST incompatibility info we read in with
 * the "basic incompatibility" info we generate internally from the
 * FSTs themselves.  The final result is saved as "cinfo.inc_edges".
 */

	static
	void
init_inc_edges (

struct cinfo *		cip,		/* IN - compatibility info */
int **			incompat,	/* IN - input incompatibility lists */
int *			icounts		/* IN - lengths of incompat lists */
)
{
int			i;
int			j;
int			k;
int			nedges;
int			total;
int *			ip1;
int *			ip2;
int *			ip3;
int *			ip4;
int *			ip5;
int *			flist;
int *			newcounts;
int **			inc_edges;
int **			basic_incompat;
bitmap_t *		edge_mask;

	nedges = cip -> num_edges;
	edge_mask = cip -> initial_edge_mask;

	inc_edges = NEWA (nedges + 1, int *);
	newcounts = NEWA (nedges, int);

	basic_incompat = compute_basic_incompat (cip);

	flist = NEWA (nedges, int);

	total = 0;
	for (i = 0; i < nedges; i++) {
		k = icounts [i];
		ip1 = flist;
		ip2 = incompat [i];
		sort_ints (ip2, k);
		ip3 = ip2 + k;
		ip4 = basic_incompat [i];
		ip5 = basic_incompat [i + 1];
		for (;;) {
			for (;;) {
				if (ip2 >= ip3) break;
				j = *ip2;
				if (BITON (edge_mask, j)) break;
				++ip2;
			}
			if (ip2 >= ip3) {
				while (ip4 < ip5) {
					*ip1++ = *ip4++;
				}
				break;
			}
			if (ip4 >= ip5) {
				while (ip2 < ip3) {
					j = *ip2++;
					if (NOT BITON (edge_mask, j)) continue;
					*ip1++ = j;
				}
				break;
			}
			k = *ip4;
			if (j <= k) {
				if (j EQ k) {
					++ip4;
				}
				*ip1++ = j;
				++ip2;
			}
			else {
				*ip1++ = k;
				++ip4;
			}
		}
		k = ip1 - flist;
		newcounts [i] = k;
		total += k;
		ip2 = NULL;
		if (k > 0) {
			ip2 = NEWA (k, int);
			for (j = 0; j < k; j++) {
				ip2 [j] = flist [j];
			}
		}
		inc_edges [i] = ip2;
	}
	free((char *) flist);

	free (basic_incompat [0]);
	free (basic_incompat);

	ip1 = NEWA (total, int);
	for (i = 0; i < nedges; i++) {
		ip2 = inc_edges [i];
		inc_edges [i] = ip1;

		if (ip2 EQ NULL) continue;

		k = newcounts [i];
		for (j = 0; j < k; j++) {
			*ip1++ = ip2 [j];
		}
		free ((char *) ip2);
	}
	inc_edges [i] = ip1;

	free ((char *) newcounts);

	if (ip1 - inc_edges [0] NE total) {
		fatal ("init_inc_edges: Bug 1.");
	}

	cip -> inc_edges = inc_edges;

	verify_symmetric (cip -> inc_edges, nedges);
}

/*
 * This routine computes for each FST, a list of those FSTs having
 * at least 2 vertices in common.
 */

	int **
compute_basic_incompat (

struct cinfo *		cip	/* IN/OUT - compatibility info. */
)
{
int			i;
int			j;
int			k;
int			t;
int			fs;
int			common;
int			nedges;
int			kmasks;
int			nmasks;
int			total;
bitmap_t *		fsmask;
bitmap_t *		edge_mask;
bitmap_t *		tmask;
int *			ep1;
int *			ep2;
int *			vp1;
int *			vp2;
int *			flist;
int **			incompat;
int *			counts;

	kmasks	  = cip -> num_vert_masks;
	nedges	  = cip -> num_edges;
	nmasks	  = cip -> num_edge_masks;
	edge_mask = cip -> initial_edge_mask;

	incompat = NEWA (nedges + 1, int *);
	counts	 = NEWA (nedges, int);
	flist	 = NEWA (nedges, int);

	for (i = 0; i < nedges; i++) {
		incompat [i] = NULL;
		counts [i]   = 0;
	}

	fsmask = NEWA (nmasks, bitmap_t);
	for (i = 0; i < nmasks; i++) {
		fsmask [i] = 0;
	}

	tmask = NEWA (kmasks, bitmap_t);
	for (i = 0; i < kmasks; i++) {
		tmask [i] = 0;
	}

	/* Compute the list of (lists of) basically incomatible FSTs... */
	total = 0;
	for (i = 0; i < nedges; i++) {
		incompat [i] = NULL;
		if (NOT BITON (edge_mask, i)) continue;
		/* Develop list of all FSTs adjacent to FST i... */
		k = 0;
		vp1 = cip -> edge [i];
		vp2 = cip -> edge [i + 1];
		while (vp1 < vp2) {
			t = *vp1++;
			SETBIT (tmask, t);
			ep1 = cip -> term_trees [t];
			ep2 = cip -> term_trees [t + 1];
			while (ep1 < ep2) {
				fs = *ep1++;
				if (BITON (fsmask, fs)) continue;
				if (NOT BITON (edge_mask, fs)) continue;
				if (fs EQ i) continue;
				SETBIT (fsmask, fs);
				flist [k++] = fs;
			}
		}
		ep1 = &flist [0];
		ep2 = &flist [k];
		k = 0;
		while (ep1 < ep2) {
			fs = *ep1++;
			CLRBIT (fsmask, fs);
			/* Count number of vertices in common. */
			common = 0;
			vp1 = cip -> edge [fs];
			vp2 = cip -> edge [fs + 1];
			while (vp1 < vp2) {
				j = *vp1++;
				if (BITON (tmask, j)) {
					++common;
				}
			}
			if (common >= 2) {
				/* Too many in common!  Retain... */
				flist [k++] = fs;
			}
		}
		counts [i] = k;
		total += k;
		if (k > 0) {
			/* Save off sorted list of incompatible FSTs. */
			sort_ints (flist, k);
			ep1 = NEWA (k, int);
			incompat [i] = ep1;
			ep2 = flist;
			for (j = 0; j < k; j++) {
				*ep1++ = *ep2++;
			}
		}
		vp1 = cip -> edge [i];
		vp2 = cip -> edge [i + 1];
		while (vp1 < vp2) {
			j = *vp1++;
			CLRBIT (tmask, j);
		}
	}
	free ((char *) tmask);
	free ((char *) fsmask);
	free ((char *) flist);

	/* Now allocate and copy into contiguous memory... */
	ep1 = NEWA (total, int);
	for (i = 0; i < nedges; i++) {
		ep2 = incompat [i];
		incompat [i] = ep1;
		if (ep2 EQ NULL) continue;
		k = counts [i];
		for (j = 0; j < k; j++) {
			*ep1++ = ep2 [j];
		}
		free ((char *) ep2);
	}
	incompat [i] = ep1;

	free ((char *) counts);

	if (ep1 - incompat [0] NE total) {
		fatal ("compute_basic_incompat: Bug 1.");
	}

	return (incompat);
}

/*
 * Verify that the given "incompatibility" relation is symmetric.
 */

	static
	void
verify_symmetric (

int **		incompat,		/* IN - incompatibility lists */
int		nedges			/* IN - number of FSTs */
)
{
int		i;
int		j;
int **		begp;
int **		endp;

	begp = NEWA (nedges, int *);
	endp = NEWA (nedges, int *);

	for (i = 0; i < nedges; i++) {
		begp [i] = incompat [i];
		endp [i] = incompat [i + 1];
	}

	for (i = 0; i < nedges; i++) {
		for (;;) {
			if (begp [i] >= endp [i]) break;
			j = begp [i] [0];
			if (j > i) break;
			++(begp [i]);
			if ((begp [j] >= endp [j]) OR
			    (begp [j] [0] NE i)) {
				/* Incompatibility array not symmetric! */
				fatal ("verify_symmetric: Bug 1.");
			}
			++(begp [j]);
		}
	}

	for (i = 0; i < nedges; i++) {
		if (begp [i] NE endp [i]) {
			/* Incompatibility array not symmetric! */
			fatal ("verify_symmetric: Bug 2.");
		}
	}

	free ((char *) endp);
	free ((char *) begp);
}

/*
 * This routine determines the version number of the input data.
 */

	static
	int
get_version (void)

{
int		c;
int		n;
int		version;

	/* Strip any white space... */
	do {
		c = getchar ();
	} while ((c EQ ' ') OR
		 (c EQ '\t') OR
		 (c EQ '\n') OR
		 (c EQ '\f'));

	if (c < 0) {
		fprintf (stderr, "Unexpected EOF!\n");
		exit (1);
	}
	if (('0' <= c) AND (c <= '9')) {
		/* No version number input -- version 0. */
		ungetc (c, stdin);
		return (P1IO_VERSION_0);
	}

	if (c NE 'V') {
		fprintf (stderr, "Bad data version number!\n");
		exit (1);
	}
	version = -1;
	n = scanf ("%d", &version);
	if (n NE 1) {
		fprintf (stderr, "Data version number not found!\n");
		exit (1);
	}
	if (version <= P1IO_VERSION_0) {
		fprintf (stderr, "Corrupt version number!\n");
		exit (1);
	}
	if (version > CURRENT_P1IO_VERSION) {
		fprintf (stderr, "Version number out of range!\n");
		exit (1);
	}

	return (version);
}

/*
 * This routine reads a single decimal number from stdin.
 */

	static
	int
get_d (void)

{
int		i, n;

	n = scanf (" %d", &i);
	if (n NE 1) {
		fprintf (stderr, "Unexpected EOF!\n");
		exit (1);
	}

	return (i);
}


/*
 * This routine reads in a 32-bit bitmap element as a hex number
 * from stdin.
 */

	static
	bitmap_t
get_x (void)

{
bitmap_t	i;
int		n;

	n = scanf (" %lx", &i);
	if (n NE 1) {
		fprintf (stderr, "Unexpected EOF!\n");
		exit (1);
	}

	return (i);
}

/*
 * This routine reads in an entire line as a string.
 */

	static
	void
get_line (

char *		buf,		/* OUT - text of line */
int		buflen		/* IN - length of buffer */
)
{
int		n;
char *		p;

	p = fgets (buf, buflen, stdin);
	if (p EQ NULL) {
		fprintf (stderr, "Unexpected EOF!\n");
		exit (1);
	}
	n = strlen (buf);
	if (n > 0) {
		buf [n - 1] = '\0';	/* Discard newline. */
	}
}

/*
 * This routine reads a decimal floating point number.
 */

	static
	double
get_dec_double (void)

{
double		x;

	if (scanf (" %lg", &x) NE 1) {
		fprintf (stderr, "Expected floating point number.\n");
		exit (1);
	}

	return (x);
}

/*
 * This routine reads a HEXIDECIMAL floating point number.
 */

	static
	double
get_hex_double (void)

{
char		buf1 [128];

	if (scanf (" %s", buf1) NE 1) {
		fprintf (stderr, "Unexpected EOF!\n");
		exit (1);
	}

	return (hex_to_double (buf1));
}

/*
 * This routine skips the remainder of the current line.
 */

	static
	void
skip (void)

{
int		c;

	for (;;) {
		c = getchar ();
		if (c < 0) break;
		if (c EQ '\n') break;
	}
}

/*
 * This routine converts the printable ASCII hex format back into
 * a floating point number.
 */

	static
	double
hex_to_double (

char *		s		/* IN - the hex ASCII string to decode */
)
{
char		c;
char *		p1;
char *		p2;
int		digit;
int		expon;
int		esign;
double		msign;
double		mant;
double		value;

	/* Strip leading white space */
	for (;;) {
		c = *s++;
		if (c EQ ' ') continue;
		if (c EQ '\n') continue;
		if (c EQ '\t') continue;
		if (c EQ '\f') continue;
		if (c NE '\v') break;
	}

	msign = 1.0;
	if (c EQ '-') {
		msign = -1.0;
		c = *s++;
	}

	mant = 0.0;

	if (c NE '.') return (msign * mant);

	/* find end of mantissa */
	p1 = s;
	for (;;) {
		c = *s++;
		digit = hexdig (c);
		if (digit < 0) break;
	}

	/* rescan mantissa in reverse */
	for (p2 = s - 1; p2 > p1; ) {
		digit = hexdig (*--p2);
#if 0
		tracef ("	@	%d\n", digit);
#endif
		mant += ((double) digit);
		mant *= 0.0625;		/* divide by 16... */
	}

	expon = 0;
	if ((c EQ 'x') OR (c EQ 'X')) {
		c = *s++;
		esign = 1;
		if (c EQ '-') {
			esign = -1;
			c = *s++;
		}
		for (;;) {
			digit = hexdig (c);
			if (digit < 0) break;
			expon = (expon << 4) + digit;
			c = *s++;
		}
		expon *= esign;
	}

	value = msign * ldexp (mant, expon);

	return (value);
}

/*
 * Routine to identify and convert hexidecimal ASCII digits.
 */

	static
	int
hexdig (

char		c
)
{
	switch (c) {

	case '0':		return (0);
	case '1':		return (1);
	case '2':		return (2);
	case '3':		return (3);
	case '4':		return (4);
	case '5':		return (5);
	case '6':		return (6);
	case '7':		return (7);
	case '8':		return (8);
	case '9':		return (9);
	case 'A': case 'a':	return (10);
	case 'B': case 'b':	return (11);
	case 'C': case 'c':	return (12);
	case 'D': case 'd':	return (13);
	case 'E': case 'e':	return (14);
	case 'F': case 'f':	return (15);

	}

	return (-1);
}

/*
 * This routine turns off the "tmap" bit for all but the first
 * terminal in each duplicate terminal group.  This effectively
 * removes the duplicated terminals from the problem.
 */

	static
	void
remove_duplicates (

int		ndg,		/* IN - number of duplicate terminal groups */
int **		list,		/* IN - list of duplicate terminal groups */
bitmap_t *	tmap		/* IN/OUT - valid subset of terminals */
)
{
int		i;
int *		ip1;
int *		ip2;

	for (i = 0; i < ndg; i++) {
		ip1 = list [i];
		ip2 = list [i + 1];

		/* Retain the first in this group, exclude all	*/
		/* of the others...				*/
		while (++ip1 < ip2) {
			CLRBIT (tmap, *ip1);
		}
	}
}

/*
 * This routine frees up the phase 1 data that was read in by
 * read_phase_1_data.
 */

	void
free_phase_1_data (

struct cinfo *		cip		/* IN - all phase 1 data */
)
{
int			i;
int			nedges;
struct full_set *	fsp;

	nedges	= cip -> num_edges;

	free ((char *) (cip -> edge [0]));
	free ((char *) (cip -> edge));
	free ((char *) (cip -> edge_size));
	free ((char *) (cip -> cost));
	free ((char *) (cip -> tflag));

	if (cip -> full_trees NE NULL) {
		for (i = 0; i < nedges; i++) {
			fsp = cip -> full_trees [i];
			free ((char *) (fsp -> terminals));
			free ((char *) (fsp -> steiners));
			free ((char *) (fsp -> edges));
			free ((char *) fsp);
		}
		free ((char *) (cip -> full_trees));
	}

	if (cip -> pts NE NULL) {
		free ((char *) (cip -> pts));
	}
	free ((char *) (cip -> initial_vert_mask));
	free ((char *) (cip -> initial_edge_mask));
	free ((char *) (cip -> required_edges));

	if (cip -> term_trees NE NULL) {
		if (cip -> term_trees [0] NE NULL) {
			free ((char *) (cip -> term_trees [0]));
		}
		free ((char *) (cip -> term_trees));
	}

	if (cip -> inc_edges NE NULL) {
		if (cip -> inc_edges [0] NE NULL) {
			free ((char *) (cip -> inc_edges [0]));
		}
		free ((char *) (cip -> inc_edges));
	}

	if (cip -> description NE NULL) {
		free ((char *) (cip -> description));
	}
}
