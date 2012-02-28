/***********************************************************************

	File:	genps.c
	Rev:	b-1
	Date:	02/28/2001

	Copyright (c) 1993, 2001 by David M. Warme

************************************************************************

	Routines for outputting stuff.

************************************************************************

	Modification Log:

	a-1:	04/18/93	warme
		: Created.  Moved most output stuff into this file
		:  from other places.
	b-1:	02/28/2001	warme
		: Changes for 3.1 release.  Pass scaling info
		:  explicitly to these routines.  Added plot_subtour.
		:  Always define X and Y scales to be equal.
		:  Fix format in non-geometric case of plot_lp_solution.

************************************************************************/

#include "genps.h"
#include "steiner.h"


/*
 * Global Routines
 */

void		begin_plot (enum plot_size);
void		define_Plot_Terminals (struct pset *, struct scale_info *);
void		draw_segment (struct point *,
			      struct point *,
			      struct scale_info *);
void		end_plot (char *);
void		overlay_plot_subset (char *,
				     bitmap_t *,
				     struct cinfo *,
				     enum plot_size);
void		plot_full_sets (bitmap_t *, struct cinfo *, enum plot_size);
void		plot_full_sets_grouped (bitmap_t *,
					struct cinfo *,
					enum plot_size);
void		plot_lp_solution (char *,
				  double *,
				  struct cinfo *,
				  enum plot_size);
void		plot_subtour (char *,
			      double *,
			      struct cinfo *,
			      bitmap_t *,
			      enum plot_size);


/*
 * External References
 */

	/* none */


/*
 * Local Equates
 */

#define	SMALL_PLOTS_PER_PAGE	12


/*
 * Local Routines
 */

static void		announce_new_page (void);
static void		define_coordinate_axes (struct pset *,
						struct scale_info *);
static void		draw_efst (struct full_set *, struct cinfo *);
static void		draw_fst (struct full_set *, struct cinfo *);
static void		draw_rfst (struct full_set *, struct cinfo *);
static void		draw_fractional_fst (struct full_set *,
					     double,
					     struct cinfo *);
static void		fst_comment (struct full_set *);
static double		get_limit (double);
static void		page_break (void);


/*
 * Local Variables
 */

static int		current_page = 0;
static enum plot_size	current_plot_size;
static int		small_plots_on_page = 0;

/*
 * This routine emits Postscript that defines the coordinates of
 * all the terminals.  The "N DefineTerminals" procedure causes space
 * to be allocated for N terminals.  We then emit the terminals, one
 * per line with "X Y DT" procedure calls.  The "DT" procedure simply
 * stuffs the given X,Y coordinates into the next slot in the TermX and
 * TermY arrays.
 *
 * Once the terminals are defined, the X,Y coordinates of a terminal (e.g.,
 * terminal 57) can be pushed onto the Postscript stack using "57 T".
 */

	void
define_Plot_Terminals (

struct pset *		pts,		/* IN - terminals to plot */
struct scale_info *	sip		/* IN - problem scaling info */
)
{
int			i;
int			n;
struct point *		p1;
char			buf1 [32];
char			buf2 [32];

	if (pts EQ NULL) return;

	n = pts -> n;

	printf ("\n%%%%BeginSetup\n");

	define_coordinate_axes (pts, sip);

	printf ("\n%d DefineTerminals\n", n);

	p1 = &(pts -> a [0]);
	for (i = 0; i < n; i++, p1++) {
		coord_to_string (buf1, p1 -> x, sip);
		coord_to_string (buf2, p1 -> y, sip);
		(void) printf ("\t%s\t%s\tDT\n", buf1, buf2);
	}
	printf ("\n%%%%EndSetup\n\n");
}

/*
 * This routine determines appropriate ranges for the X and Y coordinate
 * axes, and emits the corresponding Postscript definitions for the
 * MinX, MaxX, MinY and MaxY variables.
 */

	static
	void
define_coordinate_axes (

struct pset *		pts,		/* IN - terminals to plot */
struct scale_info *	sip		/* IN - problem scaling info */
)
{
int			i;
int			n;
struct point *		p1;
double			x, y;
double			minxcoord, maxxcoord;
double			minycoord, maxycoord;
double			xspan, yspan, span;
double			axmin, axmax;
double			aymin, aymax;

	n = pts -> n;

	if (n < 1) {
		printf ("\n0 1 0 1 SetAxes\n");
		return;
	}

	p1 = &(pts -> a [0]);

	minxcoord = maxxcoord = p1 -> x;
	minycoord = maxycoord = p1 -> y;
	++p1;

	for (i = 1; i < n; i++, p1++) {
		x = p1 -> x;
		y = p1 -> y;
		if (x < minxcoord) {
			minxcoord = x;
		}
		else if (x > maxxcoord) {
			maxxcoord = x;
		}
		if (y < minycoord) {
			minycoord = y;
		}
		else if (y > maxycoord) {
			maxycoord = y;
		}
	}

	minxcoord = UNSCALE (minxcoord, sip);
	maxxcoord = UNSCALE (maxxcoord, sip);
	minycoord = UNSCALE (minycoord, sip);
	maxycoord = UNSCALE (maxycoord, sip);

	/* We only generate square plots having equal scales on both	*/
	/* axes.  Determine the "span" of the plot, i.e., the length of	*/
	/* each axis in the plot.					*/

	xspan = maxxcoord - minxcoord;
	yspan = maxycoord - minycoord;

	if (xspan EQ 0.0) {
		if (yspan EQ 0.0) {
			/* Single point. */
			if (maxxcoord NE 0.0) {
				if (fabs (maxxcoord) >= fabs (maxycoord)) {
					span = 2.0 * fabs (maxxcoord);
				}
				else {
					span = 2.0 * fabs (maxycoord);
				}
			}
			else if (maxycoord NE 0.0) {
				span = 2.0 * fabs (maxycoord);
			}
			else {
				/* Single point at the origin. */
				span = 2.0;
			}
		}
		else {
			span = get_limit (yspan);
		}
	}
	else if (yspan EQ 0.0) {
		span = get_limit (xspan);
	}
	else if (xspan >= yspan) {
		span = get_limit (xspan);
	}
	else {
		span = get_limit (yspan);
	}

	/* Determine the minimum x axis value. */

	if (xspan EQ 0.0) {
		goto center_x;
	}
	else if ((0.0 <= minxcoord) AND (maxxcoord <= span)) {
		axmin = 0.0;
	}
	else if ((-span <= minxcoord) AND (maxxcoord <= 0.0)) {
		axmin = -span;
	}
	else if ((-0.5 * span <= minxcoord) AND (maxxcoord <= 0.5 * span)) {
		axmin = -0.5 * span;
	}
	else {
center_x:
		/* Center the x coordinates. */
		axmin = 0.5 * (minxcoord + maxxcoord - span);
	}
	axmax = axmin + span;

	/* Determine the minimum y axis value. */

	if (yspan EQ 0.0) {
		goto center_y;
	}
	else if ((0.0 <= minycoord) AND (maxycoord <= span)) {
		aymin = 0.0;
	}
	else if ((-span <= minycoord) AND (maxycoord <= 0.0)) {
		aymin = -span;
	}
	else if ((-0.5 * span <= minycoord) AND (maxycoord <= 0.5 * span)) {
		aymin = -0.5 * span;
	}
	else {
center_y:
		/* Center the y coordinates */
		aymin = 0.5 * (minycoord + maxycoord - span);
	}
	aymax = aymin + span;

	/* Good enough for now... */
	printf ("\n%g %g %g %g SetAxes\n", axmin, axmax, aymin, aymax);
}

/*
 * This routine returns a reasonable scale limit for the given quantity.
 * These are always 1, 2 or 5 times a power of 10.
 */

	static
	double
get_limit (

double		value		/* IN - value to get scale limit for */
)
{
int		i;
double		limit;

	if (value >= 1.0) {
		limit = 1.0;
		for (i = 0; i <= 20; i++) {
			if (limit >= value) return (limit);
			if (2.0 * limit >= value) return (2.0 * limit);
			if (5.0 * limit >= value) return (5.0 * limit);
			limit *= 10.0;
		}
		return (value);
	}

	limit = 1.0;
	for (i = 0; i <= 20; i++) {
		if (0.5 < value * limit) return (1.0 / limit);
		limit *= 10.0;
		if (2.0 < value * limit) return (5.0 / limit);
		if (1.0 < value * limit) return (2.0 / limit);
	}
	return (value);
}

/*
 * This routine generates some Postscript commands to plot each full-set
 * in the given mask.
 */

	void
plot_full_sets (

bitmap_t *		fset_mask,	/* IN - subset of full-sets to plot */
struct cinfo *		cip,		/* IN - compatibility info */
enum plot_size		plot_size	/* IN - size of plot to produce */
)
{
int			i;
int			n;
struct full_set *	fsp;
int *			vp1;
int *			vp2;
char			title [256];
char			buf1 [32];

	n = cip -> num_edges;

	if ((cip -> pts NE NULL) AND (cip -> full_trees NE NULL)) {
		/* Draw the FSTs with postscript. */
		for (i = 0; i < n; i++) {
			if (NOT BITON (fset_mask, i)) continue;
			fsp = cip -> full_trees [i];
			begin_plot (plot_size);
			(void) printf ("\tPlot_Terminals\n");
			fst_comment (fsp);
			draw_fst (fsp, cip);
			dist_to_string (buf1,
					cip -> cost [i],
					&(cip -> scale));
			(void) sprintf (title, "FST %lu,  Length = %s",
					(int32u) i, buf1);
			end_plot (title);
		}
		page_break ();
	}
	else {
		/* Just print the hyperedges */
		for (i = 0; i < n; i++) {
			if (NOT BITON (fset_mask, i)) continue;
			printf ("Edge %d: ", i);
			vp1 = cip -> edge [i];
			vp2 = cip -> edge [i + 1];
			while (vp1 < vp2) {
				printf ("%d ", *vp1++);
			}
			dist_to_string (buf1,
					cip -> cost [i],
					&(cip -> scale));
			printf ("%s\n", buf1);
		}
	}
}

/*
 * This routine generates some Postscript commands to plot each full-set
 * in the given mask.  Note that we "group" as many mutually-disjoint
 * full sets as possible on each plot so as to minimize the amount of
 * paper we spew...
 */

	void
plot_full_sets_grouped (

bitmap_t *		fset_mask,	/* IN - subset of full-sets to plot */
struct cinfo *		cip,		/* IN - compatibility info */
enum plot_size		plot_size	/* IN - size of plot to produce */
)
{
int			i;
int			j;
int			n;
int			nleft;
int			kmasks;
int			nmasks;
int			nplot;
bitmap_t *		left;
bitmap_t *		tmask;
int *			plist;
struct full_set *	fsp;
struct pset *		pts;
char *			cp1;
char *			cp2;
char			fsnums [110];
char			fsnum [8];
char			title [120];

	n = cip -> num_edges;
	nmasks = cip -> num_edge_masks;
	kmasks = cip -> num_vert_masks;

	left = NEWA (nmasks, bitmap_t);
	tmask = NEWA (kmasks, bitmap_t);
	plist = NEWA (n, int);

	/* Make a local copy of all full sets left to plot. */
	for (i = 0; i < nmasks; i++) {
		left [i] = fset_mask [i];
	}
	nleft = n;

	while (nleft > 0) {
		begin_plot (plot_size);

		for (i = 0; i < kmasks; i++) {
			tmask [i] = 0;
		}

		nplot = 0;
		cp1 = &fsnums [sizeof (fsnums)];
		*--cp1 = '\0';

		for (i = n - 1; i >= 0; i--) {
			if (NOT BITON (left, i)) continue;

			/* Skip full set "i" if not disjoint	*/
			/* with all others in this plot...	*/
			fsp = cip -> full_trees [i];
			pts = fsp -> terminals;
			for (j = 0; j < pts -> n; j++) {
				if (BITON (tmask, pts -> a [j].pnum)) break;
			}
			if (j < pts -> n) continue;

			(void) sprintf (fsnum, "%lu", (int32u) i);
			for (cp2 = fsnum; *cp2 NE '\0'; cp2++) {
			}

			/* Stop if label does not fit! */
			if ((cp2 - fsnum) + 2 > (cp1 - fsnums)) break;

			while (cp2 > fsnum) {
				*--cp1 = *--cp2;
			}
			*--cp1 = ' ';
			*--cp1 = ',';

			plist [nplot++] = i;
			CLRBIT (left, i);
			--nleft;

			fsp = cip -> full_trees [i];
			fst_comment (fsp);
			draw_fst (fsp, cip);

			pts = fsp -> terminals;
			for (j = 0; j < pts -> n; j++) {
				SETBIT (tmask, pts -> a [j].pnum);
			}
		}

		(void) printf ("\tPlot_Terminals\n");
		(void) sprintf (title,
				"FST%s %s.",
				(nplot > 1) ? "s" : "",
				cp1 + 2);
		end_plot (title);
	}

	page_break ();

	free ((char *) plist);
	free ((char *) tmask);
	free ((char *) left);
}

/*
 * This routine generates some Postscript commands to plot a SUBSET of
 * the given full-sets in overlaid fashion.
 */

	void
overlay_plot_subset (

char *			title,		/* IN - title to display with. */
bitmap_t *		fset_mask,	/* IN - subset of full-sets to plot */
struct cinfo *		cip,		/* IN - compatibility info */
enum plot_size		plot_size	/* IN - size of plot to produce */
)
{
int			i;
int			n;
int *			vp1;
int *			vp2;
struct full_set *	fsp;
char			buf1 [32];

	n = cip -> num_edges;

	if ((cip -> pts NE NULL) AND (cip -> full_trees NE NULL)) {
		/* Draw the FSTs with postscript. */
		begin_plot (plot_size);
		(void) printf ("\tPlot_Terminals\n");

		for (i = 0; i < n; i++) {
			if (NOT BITON (fset_mask, i)) continue;
			fsp = cip -> full_trees [i];
			fst_comment (fsp);
			draw_fst (fsp, cip);
		}

		end_plot (title);
	}
	else {
		/* Just print the hyperedges */
		for (i = 0; i < n; i++) {
			if (NOT BITON (fset_mask, i)) continue;
			printf ("Edge %d: ", i);
			vp1 = cip -> edge [i];
			vp2 = cip -> edge [i + 1];
			while (vp1 < vp2) {
				printf ("%d ", *vp1++);
			}
			dist_to_string (buf1,
					cip -> cost [i],
					&(cip -> scale));
			printf ("%s\n", buf1);
		}
	}
}

/*
 * This routine draws a single full Steiner tree.
 */

	static
	void
draw_fst (

struct full_set *	fsp,		/* IN - FST to draw */
struct cinfo *		cip		/* IN - compatibility info */
)
{
	switch (cip -> metric) {
	case RECTILINEAR:
		draw_rfst (fsp, cip);
		break;

	case EUCLIDEAN:
		draw_efst (fsp, cip);
		break;

	default:
		fatal ("draw_fst: Bug 1.");
		break;
	}
}

/*
 * This routine draws a single rectilinear FST.
 *
 * NOTE:  We ASSUME that the given FST graph represents a left-most Hwang
 * topology!  (Beyond this, we don't really care about the ordering of
 * terminals, edges or Steiner points.)  If we therefore draw all
 * "diagonal" (i.e., neither vertical nor horizontal) edges in left-most
 * fashion (vertical segment to left of the horizontal segment), then the
 * FST will be correctly drawn as a left-most Hwang topology.
 */

	static
	void
draw_rfst (

struct full_set *	fsp,		/* IN - full set to draw */
struct cinfo *		cip		/* IN - compatibility info */
)
{
int			i;
int			nt;
struct point *		p1;
struct point *		p2;
struct pset *		terms;
struct pset *		steins;
struct edge *		ep;
char			cbufx [32];
char			cbufy [32];
char			buf1 [256];
char			buf2 [256];

	terms = fsp -> terminals;
	steins = fsp -> steiners;

	nt = terms -> n;

	ep = fsp -> edges;
	for (i = 0; i < fsp -> nedges; i++, ep++) {
		if (ep -> p1 < nt) {
			p1 = &(terms -> a [ep -> p1]);
			(void) sprintf (buf1, "%d T", p1 -> pnum);
		}
		else {
			p1 = &(steins -> a [ep -> p1 - nt]);
			coord_to_string (cbufx, p1 -> x, &(cip -> scale));
			coord_to_string (cbufy, p1 -> y, &(cip -> scale));
			(void) sprintf (buf1, "%s\t%s", cbufx, cbufy);
		}

		if (ep -> p2 < nt) {
			p2 = &(terms -> a [ep -> p2]);
			(void) sprintf (buf2, "%d T", p2 -> pnum);
		}
		else {
			p2 = &(steins -> a [ep -> p2 - nt]);
			coord_to_string (cbufx, p2 -> x, &(cip -> scale));
			coord_to_string (cbufy, p2 -> y, &(cip -> scale));
			(void) sprintf (buf2, "%s\t%s", cbufx, cbufy);
		}

		if ((p1 -> x EQ p2 -> x) OR (p1 -> y EQ p2 -> y)) {
			(void) printf ("\t%s\t%s\tS\n",
				       buf1,
				       buf2);
		}
		else if (p1 -> x < p2 -> x) {
			(void) printf ("\t%s\t%s\tC\n",
				       buf1,
				       buf2);
		}
		else {
			(void) printf ("\t%s\t%s\tC\n",
				       buf2,
				       buf1);
		}
	}
}

/*
 * This routine draws a single Euclidean FST.
 */

	static
	void
draw_efst (

struct full_set *	fsp,		/* IN - FST to draw */
struct cinfo *		cip		/* IN - compatibility info */
)
{
int			i;
int			j;
int			nt;
struct point *		p1;
struct pset *		terms;
struct pset *		steins;
struct edge *		ep;
char			buf1 [32];
char			buf2 [32];

	terms = fsp -> terminals;
	steins = fsp -> steiners;

	nt = terms -> n;

	ep = fsp -> edges;
	for (i = 0; i < fsp -> nedges; i++, ep++) {
		j = ep -> p1;
		if (j < nt) {
			p1 = &(terms -> a [j]);
			(void) printf ("\t%d T", p1 -> pnum);
		}
		else {
			p1 = &(steins -> a [j - nt]);
			coord_to_string (buf1, p1 -> x, &(cip -> scale));
			coord_to_string (buf2, p1 -> y, &(cip -> scale));
			(void) printf ("\t%s\t%s", buf1, buf2);
		}
		j = ep -> p2;
		if (j < nt) {
			p1 = &(terms -> a [j]);
			(void) printf ("\t%d T\tS\n", p1 -> pnum);
		}
		else {
			p1 = &(steins -> a [j - nt]);
			coord_to_string (buf1, p1 -> x, &(cip -> scale));
			coord_to_string (buf2, p1 -> y, &(cip -> scale));
			(void) printf ("\t%s\t%s\tS\n", buf1, buf2);
		}
	}
}

/*
 * This routine plots an LP solution.  This is a set of full sets in which
 * full set i has weight Wi, where 0 <= Wi <= 1.  The full sets with
 * weight of 1 are drawn normally.  The fractional ones are drawn as
 * gray-scale "stars" emanating from the center of mass of the terminals.
 */

	void
plot_lp_solution (

char *		title,		/* IN - title for plot. */
double *	weights,	/* IN - weight of each full set. */
struct cinfo *	cip,		/* IN - compatibility info. */
enum plot_size	plot_size	/* IN - size of plot to produce. */
)
{
int			i;
int			n;
struct full_set *	fsp;
int *			vp1;
int *			vp2;

	n = cip -> num_edges;

	if ((cip -> pts NE NULL) AND (cip -> full_trees NE NULL)) {
		/* Draw the FSTs with postscript. */
		begin_plot (plot_size);

		for (i = 0; i < n; i++) {
			if (weights [i] <= 0.000001) continue;
			fsp = cip -> full_trees [i];
			draw_fractional_fst (fsp, weights [i], cip);
		}

		(void) printf ("\tPlot_Terminals\n");

		end_plot (title);
	}
	else {
		/* Just output the weighted hyperedges. */
		printf ("\n");
		for (i = 0; i < n; i++) {
			if (weights [i] <= 0.000001) continue;
			printf ("x^%d = %10.8g\t", i, weights [i]);
			vp1 = cip -> edge [i];
			vp2 = cip -> edge [i + 1];
			while (vp1 < vp2) {
				printf (" %d", *vp1++);
			}
			printf ("\n");
		}
	}
}

/*
 * This routine plots a particular subtour violation S within an LP solution.
 * Only FSTs having at least 2 vertices in the subtour S are displayed.
 * The full sets with weight of 1 are drawn normally.  The fractional ones
 * are drawn as gray-scale "stars" emanating from the center of mass of the
 * terminals.
 */

	void
plot_subtour (

char *		title,		/* IN - title for plot. */
double *	weights,	/* IN - weight of each full set. */
struct cinfo *	cip,		/* IN - compatibility info. */
bitmap_t *	S,		/* IN - subtour to plot. */
enum plot_size	plot_size	/* IN - size of plot to produce. */
)
{
int			i;
int			j;
int			k;
int			n;
struct full_set *	fsp;
int *			vp1;
int *			vp2;

	n = cip -> num_edges;

	if ((cip -> pts NE NULL) AND (cip -> full_trees NE NULL)) {
		/* Draw the FSTs with postscript. */
		begin_plot (plot_size);

		for (i = 0; i < n; i++) {
			if (weights [i] <= 0.000001) continue;
			k = 0;
			vp1 = cip -> edge [i];
			vp2 = cip -> edge [i + 1];
			while (vp1 < vp2) {
				j = *vp1++;
				if (BITON (S, j)) {
					++k;
				}
			}
			if (k < 2) continue;
			fsp = cip -> full_trees [i];
			draw_fractional_fst (fsp, weights [i], cip);
		}

		(void) printf ("\t0.75 setgray\n");
		(void) printf ("\tPlot_Terminals\n");
		(void) printf ("\t0 setgray\n");

		for (j = 0; j < cip -> num_verts; j++) {
			if (NOT BITON (S, j)) continue;
			printf ("\t%d\tPT\n", j);
		}

		end_plot (title);
	}
	else {
		/* Just output the weighted hyperedges. */
		printf ("\n");
		for (i = 0; i < n; i++) {
			if (weights [i] <= 0.000001) continue;
			printf ("x^%d = %10.8g\t", i, weights [i]);
			vp1 = cip -> edge [i];
			vp2 = cip -> edge [i + 1];
			while (vp1 < vp2) {
				printf (" %d", *vp1++);
			}
			printf ("\n");
		}
	}
}

/*
 * This routine draws a single fractional-weight full set.  We draw these
 * as a "star" using the center-of-mass of the terminals as the hub.  The
 * weight is used to determine the darkness of the lines.
 */

	static
	void
draw_fractional_fst (

struct full_set *	fsp,		/* IN - full set to draw */
double			weight,		/* IN - weight for full set */
struct cinfo *		cip		/* IN - compatibility info */
)
{
int			i;
int			n;
coord_t			cx;
coord_t			cy;
struct pset *		terms;
struct point *		p1;
char			buf1 [32];
char			buf2 [32];

	terms = fsp -> terminals;

	if (weight + 0.000001 >= 1.0) {
		/* Draw integral full sets "normally"... */
		fst_comment (fsp);
		(void) printf ("\t0 setgray\n");
		draw_fst (fsp, cip);
		return;
	}

	fst_comment (fsp);

	n = terms -> n;

	/* Determine the coordinates of the "hub". */
	p1 = &(terms -> a [0]);
	cx = 0.0;
	cy = 0.0;
	for (i = 0; i < n; p1++, i++) {
		cx += p1 -> x;
		cy += p1 -> y;
	}
	cx /= n;
	cy /= n;
	coord_to_string (buf1, cx, &(cip -> scale));
	coord_to_string (buf2, cy, &(cip -> scale));

	(void) printf ("\t%f setgray\n", 1.0 - weight);

	p1 = &(terms -> a [0]);
	for (i = 0; i < n; p1++, i++) {
		(void) printf ("\t%d T\t%s\t%s\tS\n",
			       p1 -> pnum, buf1, buf2);
	}
}

/*
 * This routine emits a Postscript comment describing the FST.
 */

	static
	void
fst_comment (

struct full_set *	fsp		/* IN - full set to describe */
)
{
int			i;
int			n;
struct pset *		terms;
struct point *		p1;

	terms = fsp -> terminals;

	(void) printf (" %% fs%d:", fsp -> tree_num);
	n = terms -> n;
	p1 = &(terms -> a [0]);
	for (i = 0; i < n; p1++, i++) {
		(void) printf (" %lu", (int32u) (p1 -> pnum));
	}
	(void) printf ("\n");
}

/*
 * This routine draws a line segment between two points.
 */

	void
draw_segment (

struct point *		p1,	/* IN - first point */
struct point *		p2,	/* IN - second point */
struct scale_info *	sip	/* IN - problem scaling info */
)
{
char		buf1 [32];
char		buf2 [32];
char		buf3 [32];
char		buf4 [32];

	coord_to_string (buf1, p1 -> x, sip);
	coord_to_string (buf2, p1 -> y, sip);
	coord_to_string (buf3, p2 -> x, sip);
	coord_to_string (buf4, p2 -> y, sip);

	(void) printf ("	%s	%s	%s	%s	S\n",
		       buf1, buf2, buf3, buf4);
}

/*
 * This routine emits appropriate Postscript code to begin a plot of
 * the given size.  If this is the first plot on a page, then we also
 * emit DSC comments for the page number.  This helps out ghostview and
 * other Postscript previewers.
 */

	void
begin_plot (

enum plot_size		size		/* IN - size of plot to begin */
)
{
	current_plot_size = size;
	switch (size) {
	case BIG_PLOT:
		page_break ();
		announce_new_page ();
		(void) printf ("BeginPlot\n");
		break;

	case SMALL_PLOT:
		if (small_plots_on_page EQ 0) {
			announce_new_page ();
		}
		(void) printf ("BeginSmallPlot\n");
		break;

	default:
		fatal ("begin_plot: Bug 1.");
		break;
	}
}

/*
 * This routine emits appropriate Postscript code to end the plot that
 * is currently in progress.  We also track the number of finished
 * plots per page here.
 */

	void
end_plot (

char *			title		/* IN - title for plot */
)
{
	(void) printf ("  (%s)\n", title);
	switch (current_plot_size) {
	case BIG_PLOT:
		(void) printf ("EndPlot\n\n");
		break;

	case SMALL_PLOT:
		(void) printf ("EndSmallPlot2\n\n");
		++small_plots_on_page;
		if (small_plots_on_page >= SMALL_PLOTS_PER_PAGE) {
			small_plots_on_page = 0;
		}
		break;

	default:
		fatal ("end_plot: Bug 1.");
		break;
	}
}

/*
 * This routine puts a %%Page: comment into the output to mark a page
 * boundary.  This is for the benefit of Ghostview and other Postscript
 * previewers.
 */

	static
	void
announce_new_page (void)

{
	++current_page;
	(void) printf ("%%%%Page: %lu %lu\n",
		       (int32u) current_page,
		       (int32u) current_page);
}

/*
 * This routine introduces a page break.  This is needed when doing
 * small plots, and the page has not yet filled up.
 */

	static
	void
page_break (void)

{
	if (small_plots_on_page NE 0) {
		(void) printf ("FlushSmallPlot\n");
		small_plots_on_page = 0;
	}
}
