/***********************************************************************

	File:	dumpfst.c
	Rev:	b-1
	Date:	01/18/2001

	Copyright (c) 1993, 2001 by David M. Warme

************************************************************************

	Main routine for a utility to dump a set of FSTs.

************************************************************************

	Modification Log:

	a-1:	10/22/98	warme
		: Created.
	b-1:	01/21/2001	warme
		: Added -d and -h options to output only statistical
		:  summaries of the given FSTs.
		: Used global sort_ints() function.

************************************************************************/

#include "p1io.h"
#include "steiner.h"


/*
 * Global Routines
 */

int			main (int, char **);


/*
 * Local Routines
 */

static void		decode_params (int, char **);
static void		usage (void);


/*
 * Local Variables
 */

static bool		aflag = FALSE;
static bool		dflag = FALSE;
static bool		hflag = FALSE;
static bool		lflag = FALSE;
static char *		me;
static bool		sflag = FALSE;

/*
 * The main routine for the dumpfst utility.
 */

	int
main (

int		argc,
char **		argv
)
{
int			i;
int			j;
int			n;
int			nedges;
int			nmasks;
int			nverts;
int			nignore;
int			nterms;
int			nsteins;
int			npruned;
int			nrequired;
int			nmaybe;
int			fpsave;
int *			bucket;
int *			fterms;
bitmap_t *		vert_mask;
bitmap_t *		edge_mask;
struct cinfo		cinfo;
char			buf [128];
char *			s;

	fpsave = set_floating_point_double_precision ();

	setbuf (stdout, NULL);

	decode_params (argc, argv);

	init_tables ();

	read_phase_1_data (&cinfo);

	vert_mask	= cinfo.initial_vert_mask;
	edge_mask	= cinfo.initial_edge_mask;

	nedges = cinfo.num_edges;
	nmasks = cinfo.num_edge_masks;
	nverts = cinfo.num_verts;

	if (dflag) {

		/* Display statistics on the set of FSTs -- */
		/* summary info only. */

		switch (cinfo.metric) {
		case RECTILINEAR:	s = "Rectlinear";	break;
		case EUCLIDEAN:		s = "Euclidean";	break;
		case PURE_GRAPH:	s = "Graph";		break;
		default:		s = "???";		break;
		}

		printf ("Metric:		%s\n", s);

		nignore = 0;
		nterms	= 0;
		nsteins = 0;
		for (i = 0; i < nverts; i++) {
			if (NOT BITON (vert_mask, i)) {
				++nignore;
			}
			else if (cinfo.tflag [i]) {
				++nterms;
			}
			else {
				++nsteins;
			}
		}

		printf ("Vertices:		%d\n", nverts);
		printf ("  Unused:		%d\n", nignore);
		printf ("  Terminals:		%d\n", nterms);
		printf ("  Steiners:		%d\n", nsteins);

		npruned	  = 0;
		nrequired = 0;
		nmaybe	  = 0;
		for (i = 0; i < nedges; i++) {
			if (NOT BITON (edge_mask, i)) {
				++npruned;
			}
			else if (BITON (cinfo.required_edges, i)) {
				++nrequired;
			}
			else {
				++nmaybe;
			}
		}

		printf ("Edges:			%d\n", nedges);
		printf ("  Pruned:		%d\n", npruned);
		printf ("  Required:		%d\n", nrequired);
		printf ("  Undecided:		%d\n", nmaybe);
	}
	if (hflag) {
		bucket = NEWA (nverts + 1, int);
		for (i = 0; i <= nverts; i++) {
			bucket [i] = 0;
		}
		for (i = 0; i < nedges; i++) {
			if ((NOT aflag) AND (NOT BITON (edge_mask, i))) continue;
			j = cinfo.edge_size [i];
			++(bucket [j]);
		}
		printf ("Size\tCount\n----\t-----\n");
		for (i = 0; i <= nverts; i++) {
			if (bucket [i] <= 0) continue;
			printf ("%d\t%d\n", i, bucket [i]);
		}
		free ((char *) bucket);
	}

	if ((NOT dflag) AND (NOT hflag)) {

		/* Not a summary mode -- dump the FSTs. */

		for (i = 0; i < nedges; i++) {
			if ((NOT aflag) AND (NOT BITON (edge_mask, i))) continue;
			fterms = cinfo.edge [i];
			n = cinfo.edge [i + 1] - fterms;
			if (sflag) {
				sort_ints (fterms, n);
			}
			for (j = 0; j < n; j++) {
				printf (" %d", fterms [j]);
			}
			if (lflag) {
				dist_to_string (buf, cinfo.cost [i], &cinfo.scale);
				printf (" %s", buf);
			}
			printf ("\n");
		}
	}

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
bool		full_dump;

	full_dump = FALSE;

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
			case 'a':
				aflag = TRUE;
				break;

			case 'd':
				dflag = TRUE;
				break;

			case 'h':
				hflag = TRUE;
				break;

			case 'l':
				lflag = TRUE;
				full_dump = TRUE;
				break;

			case 's':
				sflag = TRUE;
				full_dump = TRUE;
				break;

			default:
				usage ();
				break;
			}
		}
		--argc;
	}

	if ((dflag OR hflag) AND full_dump) {
		/* -d and -h are summary modes that contradict */
		/* the flags pertaining to full dumps of the FSTs. */
		usage ();
	}
	if (dflag AND (NOT hflag) AND aflag) {
		/* Attempting "dumpfst -d -a"... */
		usage ();
	}
}

/*
 * This routine prints out the proper usage and exits.
 */

static char *	arg_doc [] = {
	"",
	"\t-a\tDump all FSTs, even those marked as never used.",
	"\t-d\tDisplay statistics about FSTs.",
	"\t-h\tDisplay histogram of FST sizes.",
	"\t\t(-a includes never used FSTs in histogram.)",
	"\t-l\tInclude length of each FST.",
	"\t-s\tSorts the terminals of each FST.",
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
			"\nUsage:\n"
			"\t%s [-d] [-h [-a]] <FST_file\n"
			"\t%s [-als] <FST_file\n",
			me,
			me);
	pp = &arg_doc [0];
	while ((p = *pp++) NE NULL) {
		(void) fprintf (stderr, "%s\n", p);
	}
	exit (1);
}
