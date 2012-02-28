/***********************************************************************

	File:	plotfst.c
	Rev:	c-1
	Date:	01/21/2001

	Copyright (c) 1993, 2001 by David M. Warme

************************************************************************

	Main routine for a utility to plot the FSTs in
	various ways.

************************************************************************

	Modification Log:

	b-1:	01/10/97	warme
		: Split off from old_bs.c.
		: Reading in the phase 1 data.
	c-1:	01/21/2001	warme
		: Changed name from "fsplot" to "plotfst".
		: Added -p argument to plot the point set.

************************************************************************/

#include "genps.h"
#include "steiner.h"


/*
 * Global Routines
 */

int			main (int, char **);

bool			Print_Grouped_Full_Sets		= FALSE;
bool			Print_Overlaid_Full_Sets	= FALSE;


/*
 * External References
 */

	/* none */


/*
 * Local Routines
 */

static void		decode_params (int, char **);
static void		usage (void);


/*
 * Local Variables
 */

static char *		me;
static bool		Print_Full_Sets = FALSE;
static bool		Print_Points = FALSE;

/*
 * The main routine for the plotfst utility.
 */

	int
main (

int		argc,
char **		argv
)
{
int			i;
int			nedges;
int			nmasks;
int			fpsave;
bitmap_t *		edge_mask;
bitmap_t *		all_fsets_mask;
bitmap_t *		no_fsets_mask;
int			count;
char			tbuf [20];
char			title [128];
struct cinfo		cinfo;

	fpsave = set_floating_point_double_precision ();

	setbuf (stdout, NULL);

	decode_params (argc, argv);

	init_tables ();

	read_phase_1_data (&cinfo);

	edge_mask	= cinfo.initial_edge_mask;

	convert_cpu_time (cinfo.p1time, tbuf);
	printf (" %% Phase 1: %s seconds\n", tbuf);

	/* Prepare for plotting all terminals. */
	define_Plot_Terminals (cinfo.pts, &cinfo.scale);

	nedges = cinfo.num_edges;
	nmasks = cinfo.num_edge_masks;

	all_fsets_mask = NEWA (nmasks, bitmap_t);
	no_fsets_mask = NEWA (nmasks, bitmap_t);
	for (i = 0; i < nmasks; i++) {
		all_fsets_mask [i] = 0;
		no_fsets_mask [i] = 0;
	}
	for (i = 0; i < nedges; i++) {
		SETBIT (all_fsets_mask, i);
	}

	if (Print_Points) {
		if ((cinfo.description NE NULL) AND
		    (cinfo.description [0] NE '\0')) {
			strcpy (title, cinfo.description);
		}
		else {
			sprintf (title, "%lu points", (int32u) cinfo.num_verts);
		}
		overlay_plot_subset (title, no_fsets_mask, &cinfo, BIG_PLOT);
	}
	if (Print_Full_Sets) {
		plot_full_sets (all_fsets_mask, &cinfo, SMALL_PLOT);
	}
	if (Print_Grouped_Full_Sets) {
		plot_full_sets_grouped (all_fsets_mask, &cinfo, SMALL_PLOT);
	}

	if (Print_Overlaid_Full_Sets) {
		sprintf (title,
			 "All FSTs:  %lu points,  %s seconds",
			 (int32u) cinfo.num_verts, tbuf);
		overlay_plot_subset (title, edge_mask, &cinfo, BIG_PLOT);
	}

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
			case 'f':
				Print_Full_Sets = TRUE;
				break;

			case 'g':
				Print_Grouped_Full_Sets = TRUE;
				break;

			case 'o':
				Print_Overlaid_Full_Sets = TRUE;
				break;

			case 'p':
				Print_Points = TRUE;
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
	"\t-f\tPrints all full-sets in \"fly specks\" fashion.",
	"\t-g\tPrints full sets in \"grouped fly specks\" fashion.",
	"\t-o\tPrints all full-sets in overlaid fashion.",
	"\t-p\tPrints the point set, no FSTs.",
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
			"\nUsage: %s [-fgop]\n",
			me);
	pp = &arg_doc [0];
	while ((p = *pp++) NE NULL) {
		(void) fprintf (stderr, "%s\n", p);
	}
	exit (1);
}
