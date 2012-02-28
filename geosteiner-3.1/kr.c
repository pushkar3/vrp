/***********************************************************************

	File:	kr.c
	Rev:	b-1
	Date:	02/28/2001

	Copyright (c) 1993, 2001 by David M. Warme

************************************************************************

	Utility program to compute Kahng-Robins Steiner heuristic.

************************************************************************

	Modification Log:

	a-1:	05/21/96	warme
		: Created.
	b-1:	02/28/2001	warme
		: Changes for 3.1 release.
		: Fix Intel floating point precision problems.
		: Changes for new input scaling stuff.

************************************************************************/

#include "genps.h"
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
 * Local Routines
 */

static void		decode_params (int, char **);
static void		do_kahng_robins (struct pset *, struct scale_info *);
static void		do_minimum_spanning_tree (struct pset *,
						  struct scale_info *);
static void		usage (void);


/*
 * Local Variables
 */

static int32u		cpu_time_limit;
static char *		me;
static bool		Skip_Kahng_Robins		= FALSE;
static bool		Skip_Minimum_Spanning_Tree	= FALSE;

/*
 * Utility program to compute Minimum Spanning Tree and
 * 1-Steiner heuristic without doing a full Steiner tree.
 */

	int
main (

int		argc,
char **		argv
)
{
int			i;
int			fpsave;
cpu_time_t		T1;
struct pset *		pts;
struct scale_info	scale;

	fpsave = set_floating_point_double_precision ();

	setbuf (stdout, NULL);

	decode_params (argc, argv);

	pts = get_points (stdin, &scale);

	init_output_conversion (pts, RECTILINEAR, &scale);

	/* Enable CPU time limitation, if any. */
	start_limiting_cpu_time (cpu_time_limit);

	/* Prepare for plotting all terminals. */
	define_Plot_Terminals (pts, &scale);

	if (NOT Skip_Minimum_Spanning_Tree) {
		do_minimum_spanning_tree (pts, &scale);
	}
	if (NOT Skip_Kahng_Robins) {
		do_kahng_robins (pts, &scale);
	}

	free ((char *) pts);

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
			case 'k':
				Skip_Kahng_Robins = TRUE;
				break;

			case 'l':
				if (*ap EQ '\0') {
					if (argc <= 0) {
						usage ();
					}
					ap = *argv++;
					--argc;
				}
				if (NOT decode_cpu_time_limit (ap, &cpu_time_limit)) {
					usage ();
				}
				ap = "";
				break;

			case 'm':
				Skip_Minimum_Spanning_Tree = TRUE;
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
	"\t-k\tDisables computation of Kahng-Robins Tree.",
	"\t-m\tDisables computation of Minimum Spanning Tree.",
	"",
	"\tExample CPU times are:",
	"\t\t-l 3days2hours30minutes15seconds",
	"\t\t-l1000seconds -l1000 -l 2h30m",
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
			"\nUsage: %s [-km] [-l cpu-time-limit]\n",
			me);
	pp = &arg_doc [0];
	while ((p = *pp++) NE NULL) {
		(void) fprintf (stderr, "%s\n", p);
	}
	exit (1);
}

/*
 * This routine computes and plots the Minimum Spanning Tree for the
 * given point set.
 */

	static
	void
do_minimum_spanning_tree (

struct pset *		pts,		/* IN - the set of terminals */
struct scale_info *	sip		/* IN - problem scaling info */
)
{
int			i;
int			n;
int			nedges;
cpu_time_t		T0;
cpu_time_t		T1;
struct edge *		ep;
struct point *		p1;
struct point *		p2;
dist_t			mst_len;
char			tbuf [20];
char			title [80];
struct edge *		m;
char			buf1 [32];

	n = pts -> n;

	m = NEWA (n - 1, struct edge);

	T0 = get_cpu_time ();
	nedges = rect_mst (pts, &m [0], NULL);

	T1 = get_cpu_time ();

	if (nedges NE n - 1) {
		fatal ("do_minimum_spanning_tree: Bug 1.");
	}

	(void) printf ("\n  %%  Minimum Spanning Tree\n");

	begin_plot (BIG_PLOT);
	(void) printf ("\tPlot_Terminals\n");
	mst_len = 0;
	for (i = 0; i < nedges; i++) {
		ep = &m [i];
		p1 = &(pts -> a [ep -> p1]);
		p2 = &(pts -> a [ep -> p2]);
		draw_segment (p1, p2, sip);
		mst_len += ep -> len;
	}
	convert_cpu_time (T1 - T0, tbuf);
	dist_to_string (buf1, mst_len, sip);
	(void) sprintf (title, "Minimum Spanning Tree:  %lu points,  length = %s,  %s seconds",
			(int32u) n, buf1, tbuf);
	end_plot (title);

	free ((char *) m);
}

/*
 * This routine computes an approximate Steiner Minimal Tree using
 * the Kahng-Robins heuristic, and plots the result.
 */

	static
	void
do_kahng_robins (

struct pset *		pts,		/* IN - the set of terminals */
struct scale_info *	sip		/* IN - problem scaling info */
)
{
int			i;
int			n;
int			nedges;
cpu_time_t		T0;
cpu_time_t		T1;
struct pset *		newpts;
struct edge *		ep;
struct point *		p1;
struct point *		p2;
dist_t			kr_len;
char			tbuf [20];
char			title [80];
struct edge *		m;
char			buf1 [32];

	n = pts -> n;

	/* Make a temporary copy, big enough to hold the Steiner points too. */
	newpts = NEW_PSET (n + n - 2);
	COPY_PSET (newpts, pts);

	m = NEWA (n + n - 3, struct edge);

	T0 = get_cpu_time ();

	nedges = kahng_robins (newpts, 0, &m [0]);

	T1 = get_cpu_time ();

	(void) printf ("\n  %%  Kahng Robins\n");

	begin_plot (BIG_PLOT);
	(void) printf ("\tPlot_Terminals\n");
	kr_len = 0;
	for (i = 0; i < nedges; i++) {
		ep = &m [i];
		p1 = &(newpts -> a [ep -> p1]);
		p2 = &(newpts -> a [ep -> p2]);
		draw_segment (p1, p2, sip);
		kr_len += ep -> len;
	}
	convert_cpu_time (T1 - T0, tbuf);
	dist_to_string (buf1, kr_len, sip);
	(void) sprintf (title, "1-Steiner Heuristic:  %lu points,  length = %s,  %s seconds",
			(int32u) n, buf1, tbuf);
	end_plot (title);

	free ((char *) m);
	free ((char *) newpts);
}
