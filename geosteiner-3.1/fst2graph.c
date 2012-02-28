/***********************************************************************

	File:	fst2graph.c
	Rev:	b-1
	Date:	01/19/2001

	Copyright (c) 1995, 2001 by David M. Warme & Martin Zachariasen

************************************************************************

	The main routine to generate ordinary graphs from FSTs.
	There are two possibilities:
	1. Reduced grid-graph for rectilinear problem,
	2. Graph with all terminals and Steiner points as nodes
	   and all line segments as edges (Euclidean and rectilinear
	   problem).

************************************************************************

	Modification Log:

	a-1:	05/18/95	warme
		: Created.
	b-1:	01/19/2001	martinz
		: Changed int16u to int32u in order to be able to handle
		:  bigger instances (also changed macros).
		: Changed name from redg.c to fst2graph.c
		: Support for OR-Library format and SteinLib format
		:  (made old format obsolete).
		: Added non-grid graph generation.

************************************************************************/

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
 * Local Equates
 */

#define	UP_EDGE		0x80000000	/* edge exits point and goes up */
#define	RIGHT_EDGE	0x40000000	/* edge exits point and goes right */
#define	GMASK		0x3FFFFFFF	/* mask of other bits */


/*
 * Local Types
 */

struct grid {
	int		nt;
	int		ns;
	coord_t *	x_coord;
	coord_t *	y_coord;
	int *		xindex;
	int *		yindex;
	int32u *	gridp;
};


/*
 * Local Routines
 */

static void		compute_edge_graph (bitmap_t *, bitmap_t *,
					    struct cinfo *);
static void		compute_grid_graph (bitmap_t *, bitmap_t *,
					    struct cinfo *);
static void		decode_params (int, char **);
static void		draw_bb_grid (int, int, int, int);
static void		draw_bb_grid_horizontal (int, int, int);
static void		draw_bb_grid_vertical (int, int, int);
static void		draw_full_sets_on_grid (bitmap_t *, struct cinfo *);
static void		edge_down (int32u *, int, int, int, dist_t,
				   struct cinfo *);
static void		edge_left (int32u *, int, int, int, dist_t,
				   struct cinfo *);
static void		edge_right (int32u *, int, int, int, dist_t,
				    struct cinfo *);
static void		edge_up (int32u *, int, int, int, dist_t,
				 struct cinfo *);
static void		identify_steiner_points (void);
static int		map_x (coord_t);
static int		map_y (coord_t);
static void		output_edge (int, int, dist_t, struct cinfo *);
static void		sort_x_inc (struct pset *);
static void		sort_y_inc (struct pset *);
static void		usage (void);


/*
 * Local Variables
 */

static bool		count_edges;
static struct grid	grid;
static char *		me;
static int		number_of_edges;
static bool		Print_Grid_Graph	= TRUE;
static bool		Print_ORLibrary_Format	= TRUE;
static bool		Print_SteinLib_Format	= FALSE;
static bool		Print_Unscaled = FALSE;
static char *		description = NULL;

/*
 * This is the main routine for computing Reduced Grid-Graphs.
 */

	int
main (

int		argc,
char **		argv
)
{
int			fpsave;
bitmap_t *		vert_mask;
bitmap_t *		edge_mask;
struct cinfo		cinfo;

	fpsave = set_floating_point_double_precision ();

	setbuf (stdout, NULL);

	decode_params (argc, argv);

	init_tables ();

	read_phase_1_data (&cinfo);
	if (description NE NULL) {
		free ((char *) cinfo.description);
		cinfo.description = gst_strdup (description);
	}

	if ((cinfo.metric NE RECTILINEAR) AND (cinfo.metric NE EUCLIDEAN)) {
		fprintf (stderr, "This only be done for geometric (Euclidean or rectilinear) FSTs\n");
		exit (1);
	}

	vert_mask	= cinfo.initial_vert_mask;
	edge_mask	= cinfo.initial_edge_mask;

	if ((cinfo.metric EQ RECTILINEAR) AND (NOT Print_Unscaled)) {
		/* Force printing of integer data. */
		cinfo.scale.min_precision = 0;
		if (cinfo.scale.scale > 0) {
			cinfo.scale.scale = 0;
		}
	}

	if ((cinfo.metric EQ EUCLIDEAN) OR (!Print_Grid_Graph)) {
		compute_edge_graph (vert_mask, edge_mask, &cinfo);
	}
	else {
		compute_grid_graph (vert_mask, edge_mask, &cinfo);
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

			case 'e':
				Print_Grid_Graph	= FALSE;
				break;

			case 's':
				Print_ORLibrary_Format	= FALSE;
				Print_SteinLib_Format	= TRUE;
				break;

			case 'u':
				Print_Unscaled = TRUE;
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
	"Reads FST info from stdin. Produces an ordinary graph",
	"on stdout which is either:",
	" - graph of all line segments in all FSTs (default for Euclidean problem)",
	" - reduced grid graph (default for rectilinear problem)",
	"Output data is printed in OR-Library format.",
	"Distances in the rectilinear problem are scaled to integers."
	"",
	"\t-d txt\tDescription of problem instance.",
	"\t-s\tPrints data in SteinLib format.",
	"\t-e\tGenerate edge graph for the rectilinear problem.",
	"\t-u\tOutput unscaled (fractional) data.",
	"",
	NULL
};

	static
	void
usage (void)

{
char **		pp;
char *		p;

	(void) fprintf (stderr, "\nUsage: %s [-seu] [-d description]\n", me);
	pp = &arg_doc [0];
	while ((p = *pp++) NE NULL) {
		(void) fprintf (stderr, "%s\n", p);
	}
	exit (1);
}

/*
 * This routine computes the edge-graph obtained by taking the union
 * of all of all line segments in all FSTs
 */

	static
	void
compute_edge_graph (

bitmap_t *	tmap,		/* IN - valid set of terminals */
bitmap_t *	fset_mask,	/* IN - valid FSTs */
struct cinfo *	cip		/* IN - compatibility info */
)
{
int			i, j, i1, i2;
int			sp_index;
int			sp_total;
int			n2edges;
int			nedges;
struct full_set *	fsp;
struct pset *		terms;
struct point *		p1;
struct point *		p2;
dist_t			d = 0.0;
char			buf1 [128];
char			buf2 [128];

	/* Total number of Steiner points and ordinary edges */
	nedges = cip -> num_edges;
	sp_total = 0;
	n2edges	 = 0;
	for (i = 0; i < nedges; i++) {
		if (NOT BITON (fset_mask, i)) continue;
		sp_total += cip -> full_trees [i] -> steiners -> n;
		n2edges	 += cip -> full_trees [i] -> nedges;
	}

	/* Write header */
	if (Print_ORLibrary_Format) {
		printf ("%d %d\n", cip -> num_verts + sp_total, n2edges);
	}

	if (Print_SteinLib_Format) {
		printf ("33d32945 STP File, STP Format Version 1.00\n"
			"Section Comment\n"
			"Name    \"%s\"\n"
			"Creator \"GeoSteiner\"\n"
			"Remark  \"Reduced graph from FST generator\"\n"
			"End\n\n"
			"Section Graph\n"
			"Nodes %d\n"
			"Edges %d\n",
			cip -> description, cip -> num_verts + sp_total, n2edges);
	}

	/* Write (ordinary) graph edges */
	sp_index = cip -> num_verts;
	for (i = 0; i < nedges; i++) {
		if (NOT BITON (fset_mask, i)) continue;
		fsp = cip -> full_trees [i];
		terms = fsp -> terminals;
		for (j = 0; j < fsp -> nedges; j++) {

			/* Compute length of this edge */
			p1 = (fsp -> edges[j].p1 < terms -> n)
				? &(terms -> a[ fsp -> edges[j].p1 ])
				: &(fsp -> steiners -> a [fsp -> edges[j].p1 - terms -> n]);
			p2 = (fsp -> edges[j].p2 < terms -> n)
				? &(terms -> a[ fsp -> edges[j].p2 ])
				: &(fsp -> steiners -> a [fsp -> edges[j].p2 - terms -> n]);

			i1 = (fsp -> edges[j].p1 < terms -> n)
				? terms -> a[ fsp -> edges[j].p1 ].pnum + 1
				: fsp -> edges[j].p1 - terms -> n + sp_index + 1;
			i2 = (fsp -> edges[j].p2 < terms -> n)
				? terms -> a[ fsp -> edges[j].p2 ].pnum + 1
				: fsp -> edges[j].p2 - terms -> n + sp_index + 1;

			if (Print_SteinLib_Format) {
				printf ("E ");
			}

			if (cip -> metric EQ EUCLIDEAN) {
				d = EDIST(p1, p2);
			}
			else {
				d = RDIST(p1, p2);
			}
			dist_to_string (buf1, d, &(cip -> scale));
			printf ("%d %d %s\n", i1, i2, buf1);
		}
		if (fsp -> steiners NE NULL) {
			sp_index += fsp -> steiners -> n;
		}
	}

	/* Write terminals (and coordinates in STP-format) */
	if (Print_ORLibrary_Format) {
		printf ("%d\n", cip -> num_verts);
		for (i = 0; i < cip -> num_verts; i++) {
			printf ("%d\n", i+1);
		}
	}

	if (Print_SteinLib_Format) {
		printf ("End\n\n"
			"Section Terminals\n"
			"Terminals %d\n",
			cip -> num_verts);
		for (i = 0; i < cip -> num_verts; i++) {
			printf ("T %d\n", i+1);
		}
		printf ("End\n\n"
			"Section Coordinates\n");
		for (i = 0; i < cip -> num_verts; i++) { /* terminals */
			dist_to_string (buf1,
					cip -> pts -> a[i].x,
					&(cip -> scale));
			dist_to_string (buf2,
					cip -> pts -> a[i].y,
					&(cip -> scale));
			printf ("DD %d %s %s\n", i+1, buf1, buf2);
		}
		sp_index = cip -> num_verts + 1;
		for (i = 0; i < nedges; i++) { /* Steiner points */
			fsp = cip -> full_trees [i];
			for (j = 0; j < fsp -> steiners -> n; j++) {
				dist_to_string (buf1,
						fsp -> steiners -> a[j].x,
						&(cip -> scale));
				dist_to_string (buf2,
						fsp -> steiners -> a[j].y,
						&(cip -> scale));
				printf ("DD %d %s %s\n", sp_index++, buf1, buf2);
			}
		}
		printf ("End\n\n"
			"EOF\n");
	}
}

/*
 * This routine computes the grid-graph obtained by taking the union
 * of all FSTs
 */

	static
	void
compute_grid_graph (

bitmap_t *	tmap,		/* IN - valid set of terminals */
bitmap_t *	fset_mask,	/* IN - valid FSTs */
struct cinfo *	cip		/* IN - compatibility info */
)
{
int		i;
int		j;
int		k;
int		nverts;
int		nverts2;
int		v1;
int		spi;
int32u *	gridp;
coord_t		coord;
dist_t		prev_coord;
int		prev_index;
struct pset *	tmp;
char		buf1 [128];
char		buf2 [128];

	nverts = cip -> num_verts;

	grid.nt		= nverts;
	grid.ns		= 0;
	grid.x_coord	= NEWA (nverts, coord_t);
	grid.y_coord	= NEWA (nverts, coord_t);
	grid.xindex	= NEWA (nverts, int);
	grid.yindex	= NEWA (nverts, int);

	/* Compute map giving index of each terminal in sequence when	*/
	/* sorted by increasing X coordinate.				*/
	tmp = NEW_PSET (nverts);
	tmp -> n = nverts;
	for (i = 0; i < nverts; i++) {
		tmp -> a [i] = cip -> pts -> a [i];
	}
	sort_x_inc (tmp);
	prev_coord = INF_DISTANCE;
	prev_index = -1;
	for (i = 0; i < nverts; i++) {
		coord = tmp -> a [i].x;
		grid.x_coord [i] = coord;
		j = tmp -> a [i].pnum;
		if (coord NE prev_coord) {
			prev_index = i;
		}
		grid.xindex [j] = prev_index;
		prev_coord = coord;
	}

	/* Compute map giving index of each terminal in sequence when	*/
	/* sorted by increasing Y coordinate.				*/
	tmp -> n = nverts;
	for (i = 0; i < nverts; i++) {
		tmp -> a [i] = cip -> pts -> a [i];
	}
	sort_y_inc (tmp);
	prev_coord = INF_DISTANCE;
	prev_index = -1;
	for (i = 0; i < nverts; i++) {
		coord = tmp -> a [i].y;
		grid.y_coord [i] = tmp -> a [i].y;
		j = tmp -> a [i].pnum;
		if (coord NE prev_coord) {
			prev_index = i;
		}
		grid.yindex [j] = prev_index;
		prev_coord = coord;
	}

	free ((char *) tmp);

	/* Allocate and zero matrix to hold the grid... */
	nverts2 = (nverts + 1) * nverts;
	grid.gridp = NEWA (nverts2, int32u);
	for (i = 0; i < nverts2; i++) {
		grid.gridp [i] = 0;
	}
	grid.gridp += nverts;

	/* Set vertex number for each terminal in the grid... */
	for (i = 0; i < nverts; i++) {
		if (NOT BITON (tmap, i)) {
			/* Should not have any duplicate terminals! */
			fatal ("compute_grid_graph: Bug 1.");
		}
		j = grid.xindex [i];
		k = grid.yindex [i];
		grid.gridp [k * nverts + j] = i + 1;
	}

	draw_full_sets_on_grid (fset_mask, cip);

	identify_steiner_points ();

	/* Count every edge... */
	count_edges = TRUE;
	number_of_edges = 0;
	gridp = grid.gridp;
	for (i = 0; i < nverts; i++) {
		for (j = 0; j < nverts; j++, gridp++) {
			v1 = (*gridp & GMASK);
			if (v1 EQ 0) continue;
			if ((*gridp & RIGHT_EDGE) NE 0) {
				edge_right (gridp, v1, i, j, 0, cip);
			}
			if ((*gridp & UP_EDGE) NE 0) {
				edge_up (gridp, v1, i, j, 0, cip);
			}
			if ((j > 0) AND ((gridp [-1] & RIGHT_EDGE) NE 0)) {
				edge_left (gridp, v1, i, j, 0, cip);
			}
			if ((i > 0) AND ((gridp [-nverts] & UP_EDGE) NE 0)) {
				edge_down (gridp, v1, i, j, 0, cip);
			}
		}
	}

	/* Output total number of vertices, plus number of terminals. */
	if (Print_ORLibrary_Format) {
		printf ("%d %d\n", nverts + grid.ns, number_of_edges);
	}

	if (Print_SteinLib_Format) {
		printf ("33d32945 STP File, STP Format Version 1.00\n"
			"Section Comment\n"
			"Name    \"%s\"\n"
			"Creator \"GeoSteiner\"\n"
			"Remark  \"Reduced graph from FST generator\"\n"
			"End\n\n"
			"Section Graph\n"
			"Nodes %d\n"
			"Edges %d\n",
			cip -> description, nverts + grid.ns, number_of_edges);
	}

	/* Output every edge... */
	count_edges = FALSE;
	gridp = grid.gridp;
	for (i = 0; i < nverts; i++) {
		for (j = 0; j < nverts; j++, gridp++) {
			v1 = (*gridp & GMASK);
			if (v1 EQ 0) continue;
			if ((*gridp & RIGHT_EDGE) NE 0) {
				edge_right (gridp, v1, i, j, 0, cip);
			}
			if ((*gridp & UP_EDGE) NE 0) {
				edge_up (gridp, v1, i, j, 0, cip);
			}
			if ((j > 0) AND ((gridp [-1] & RIGHT_EDGE) NE 0)) {
				edge_left (gridp, v1, i, j, 0, cip);
			}
			if ((i > 0) AND ((gridp [-nverts] & UP_EDGE) NE 0)) {
				edge_down (gridp, v1, i, j, 0, cip);
			}
		}
	}

	/* Write terminals (and coordinates in STP-format) */
	if (Print_ORLibrary_Format) {
		printf ("%d\n", nverts);
		for (i = 0; i < nverts; i++) {
			printf ("%d\n", i+1);
		}
	}

	if (Print_SteinLib_Format) {
		printf ("End\n\n"
			"Section Terminals\n"
			"Terminals %d\n",
			nverts);
		for (i = 0; i < nverts; i++) {
			printf ("T %d\n", i+1);
		}
		printf ("End\n\n"
			"Section Coordinates\n");

		/* Output coordinates of each terminal... */
		for (i = 0; i < nverts; i++) {
			coord_to_string (buf1,
					 cip -> pts -> a [i].x,
					 &(cip -> scale));
			coord_to_string (buf2,
					 cip -> pts -> a [i].y,
					 &(cip -> scale));
			printf ("DD %d %s %s\n", i+1, buf1, buf2);
		}

		/* Output coordinates of each steiner point... */
		spi = nverts+1;
		gridp = grid.gridp;
		for (i = 0; i < nverts; i++) {
			for (j = 0; j < nverts; j++, gridp++) {
				if ((*gridp & GMASK) > nverts) {
					coord_to_string (buf1,
							 grid.x_coord [j],
							 &(cip -> scale));
					coord_to_string (buf2,
							 grid.y_coord [i],
							 &(cip -> scale));
					printf ("DD %d %s %s\n",
						spi, buf1, buf2);
					spi++;
				}
			}
		}
		printf ("End\n\n"
			"EOF\n");
	}


	free ((char *) (grid.gridp - nverts));
	free ((char *) grid.yindex);
	free ((char *) grid.xindex);
	free ((char *) grid.y_coord);
	free ((char *) grid.x_coord);
}

/*
 * This routine sorts the given point set in order by INCREASING X coord.
 */

	static
	void
sort_x_inc (

struct pset *		pts	/* IN/OUT - point set to sort. */
)
{
int		n;
int		h;
struct point	tmp;
coord_t		key;
struct point *	p1;
struct point *	p2;
struct point *	p3;
struct point *	p4;
struct point *	endp;

	n = pts -> n;

	endp = &(pts -> a [n]);

	for (h = 1; h <= n; h = 3*h+1) {
	}

	do {
		h = h / 3;
		p4 = &(pts -> a [h]);
		p1 = p4;
		while (p1 < endp) {
			tmp = *p1;
			key = tmp.x;
			p2 = p1;
			while (TRUE) {
				p3 = (p2 - h);
				if (p3 -> x <= key) break;
				*p2 = *p3;
				p2 = p3;
				if (p2 < p4) break;
			}
			*p2 = tmp;
			++p1;
		}
	} while (h > 1);
}

/*
 * This routine sorts the given point set in order by INCREASING Y coord.
 */

	void
sort_y_inc (

struct pset *		pts	/* IN/OUT - point set to sort. */
)
{
int		n;
int		h;
struct point	tmp;
coord_t		key;
struct point *	p1;
struct point *	p2;
struct point *	p3;
struct point *	p4;
struct point *	endp;

	n = pts -> n;

	endp = &(pts -> a [n]);

	for (h = 1; h <= n; h = 3*h+1) {
	}

	do {
		h = h / 3;
		p4 = &(pts -> a [h]);
		p1 = p4;
		while (p1 < endp) {
			tmp = *p1;
			key = tmp.y;
			p2 = p1;
			while (TRUE) {
				p3 = (p2 - h);
				if (p3 -> y <= key) break;
				*p2 = *p3;
				p2 = p3;
				if (p2 < p4) break;
			}
			*p2 = tmp;
			++p1;
		}
	} while (h > 1);
}

/*
 * This routine draws all valid FSTs onto the given grid in
 * an overlaid fashion.
 */

	static
	void
draw_full_sets_on_grid (

bitmap_t *	fset_mask,	/* IN - mask of valid FSTs */
struct cinfo *	cip		/* IN - compatibility info */
)
{
int			i;
int			j;
int			u, v, ux, uy, vx, vy;
int			nt;
int			ns;
struct full_set *	fsp;
struct pset *		terms;
struct pset *		steins;
struct edge *		ep;
int *			xi;
int *			yi;

	xi = grid.xindex;
	yi = grid.yindex;
	for (i = 0; i < cip -> num_edges; i++) {
		if (NOT BITON (fset_mask, i)) continue;
		fsp = cip -> full_trees [i];

		ep = fsp -> edges;
		terms = fsp -> terminals;
		steins = fsp -> steiners;
		nt = terms -> n;
		ns = steins -> n;
		for (j = 0; j < fsp -> nedges; j++, ep++) {
			u = ep -> p1;
			if (u < 0) {
				fatal ("draw_full_sets_on_grid: Bug 1.");
			}
			if (u < nt) {
				u = terms -> a [u].pnum;
				ux = xi [u];
				uy = yi [u];
			}
			else {
				u -= nt;
				if (u >= ns) {
					fatal ("draw_full_sets_on_grid: Bug 2.");
				}
				ux = map_x (steins -> a [u].x);
				uy = map_y (steins -> a [u].y);
			}

			v = ep -> p2;
			if (v < 0) {
				fatal ("draw_full_sets_on_grid: Bug 3.");
			}
			if (v < nt) {
				v = terms -> a [v].pnum;
				vx = xi [v];
				vy = yi [v];
			}
			else {
				v -= nt;
				if (v >= ns) {
					fatal ("draw_full_sets_on_grid: Bug 4.");
				}
				vx = map_x (steins -> a [v].x);
				vy = map_y (steins -> a [v].y);
			}

			draw_bb_grid (ux, uy, vx, vy);
		}
	}
}

/*
 * This routine determines the grid index of the given X coordinate.
 */

	static
	int
map_x (

coord_t		x
)
{
int		i;

	for (i = 0; i < grid.nt; i++) {
		if (x EQ grid.x_coord [i]) return (i);
		if (x < grid.x_coord [i]) {
			fatal ("map_x: Bug 1.");
		}
	}

	fatal ("map_x: Bug 2.");
}

/*
 * This routine determines the grid index of the given Y coordinate.
 */

	static
	int
map_y (

coord_t		y
)
{
int		i;

	for (i = 0; i < grid.nt; i++) {
		if (y EQ grid.y_coord [i]) return (i);
		if (y < grid.y_coord [i]) {
			fatal ("map_y: Bug 1.");
		}
	}

	fatal ("map_y: Bug 2.");
}

/*
 * This routine is used to draw a backbone onto the grid.  A backbone
 * consists of a vertical line segment from the first point, and a
 * horizontal line from the second point.  Each segment consists only
 * of the points between the given points and the intersections of the
 * horizontal and vertical lines.
 */

	static
	void
draw_bb_grid (

int		x0,
int		y0,
int		x1,
int		y1
)
{
int		nt;

	nt = grid.nt;
	if ((x0 < 0) OR (x0 >= nt)) {
		fatal ("draw_bb_grid: Bug 1.");
	}
	if ((y0 < 0) OR (y0 >= nt)) {
		fatal ("draw_bb_grid: Bug 2.");
	}
	if ((x1 < 0) OR (x1 >= nt)) {
		fatal ("draw_bb_grid: Bug 3.");
	}
	if ((y1 < 0) OR (y1 >= nt)) {
		fatal ("draw_bb_grid: Bug 4.");
	}

	/* Assume RFSTs are left-most topologies.  Plot corner	*/
	/* segments so that the vertical segment is to the left	*/
	/* of the horizontal segment.				*/

	if (x0 < x1) {
		draw_bb_grid_vertical (x0, y0, y1);
		draw_bb_grid_horizontal (x1, y1, x0);
	}
	else {
		draw_bb_grid_vertical (x1, y1, y0);
		draw_bb_grid_horizontal (x0, y0, x1);
	}
}

/*
 * This routine draws a vertical line onto the grid.
 */

	static
	void
draw_bb_grid_vertical (

int		x0,
int		y0,
int		y1
)
{
int		i;
int		n;
int32u *	gridp;

	if (y0 > y1) {
		i = y0;
		y0 = y1;
		y1 = i;
	}

	n = grid.nt;

	if ((x0 < 0) OR (x0 >= n) OR (y0 < 0) OR (y1 > n)) {
		fatal ("draw_bb_grid_vertical: Bug 1.");
	}

	gridp = grid.gridp;
	gridp += (n * y0 + x0);
	while (y0 < y1) {
		*gridp |= UP_EDGE;
		gridp += n;
		++y0;
	}
}

/*
 * This routine draws a horizontal line onto the grid.
 */

	static
	void
draw_bb_grid_horizontal (

int		x0,
int		y0,
int		x1
)
{
int		i;
int		n;
int32u *	gridp;

	if (x0 > x1) {
		i = x0;
		x0 = x1;
		x1 = i;
	}

	n = grid.nt;

	if ((y0 < 0) OR (y0 >= n) OR (x0 < 0) OR (x1 > n)) {
		fatal ("draw_bb_grid_horizontal: Bug 1.");
	}

	gridp = grid.gridp;
	gridp += (n * y0 + x0);
	while (x0 < x1) {
		*gridp |= RIGHT_EDGE;
		++gridp;
		++x0;
	}
}

/*
 * This routine is used to identify the steiner points on the grid.
 */

	static
	void
identify_steiner_points (void)

{
int		i;
int		j;
int		nt;
int		ns;
int		code;
int32u *	gridp;

	nt = grid.nt;
	ns = 0;
	gridp = grid.gridp;
	for (i = 0; i < nt; i++) {
		for (j = 0; j < nt; j++, gridp++) {
			code = 0;
			if ((*gridp & RIGHT_EDGE) NE 0) {
				code |= 1;
			}
			if ((*gridp & UP_EDGE) NE 0) {
				code |= 0x02;
			}
			if ((j > 0) AND ((gridp [-1] & RIGHT_EDGE) NE 0)) {
				code |= 0x04;
			}
			if ((i > 0) AND ((gridp [-nt] & UP_EDGE) NE 0)) {
				code |= 0x08;
			}
			switch (code) {
			case 0x00:
				if ((*gridp & GMASK) NE 0) {
					/* Terminal with no connections! */
					fatal ("identify_steiner_points: Bug 1.");
				}
				break;

			case 0x01:
			case 0x02:
			case 0x04:
			case 0x08:
				if ((*gridp & GMASK) EQ 0) {
					/* degree 1 vertex must be terminal! */
					fatal ("identify_steiner_points: Bug 2.");
				}
				break;

			case 0x03:
			case 0x05:
			case 0x06:
			case 0x09:
			case 0x0A:
			case 0x0C:
				break;

			default:
				/* vertex has degree 3 or 4... */
				if ((*gridp & GMASK) EQ 0) {
					/* This is a new steiner point! */
					++ns;
					*gridp = (*gridp & ~GMASK) |
						 ((nt + ns) & GMASK);
				}
				break;
			}
		}
	}

	grid.ns = ns;
}

/*
 * This routine chases the grid-graph edge going RIGHT from position
 * I, J.  We traverse at least one grid element, but stop once we
 * discover a terminal or steiner node.
 */

	static
	void
edge_right (

int32u *	gridp,		/* IN - current grid element */
int		v1,		/* IN - first vertex of edge */
int		i,		/* IN - current Y index */
int		j,		/* IN - current X index */
dist_t		total,		/* IN - total length of edge so far */
struct cinfo *	cip		/* IN - compatibility info */
)
{
int		nt;
int		v2;
int		code;

	nt = grid.nt;
	for (;;) {
		/* Step once to the right. */
		if ((j + 1) >= nt) {
			/* Ran off edge of grid! */
			fatal ("edge_right: Bug 1.");
		}
		total += (grid.x_coord [j+1] - grid.x_coord [j]);
		++j;
		++gridp;

		code = 0;
		if ((*gridp & RIGHT_EDGE) NE 0) {
			code |= 1;
		}
		if ((*gridp & UP_EDGE) NE 0) {
			code |= 0x02;
		}
		if ((j > 0) AND ((gridp [-1] & RIGHT_EDGE) NE 0)) {
			code |= 0x04;
		}
		if ((i > 0) AND ((gridp [-nt] & UP_EDGE) NE 0)) {
			code |= 0x08;
		}
		if ((code & 0x04) EQ 0) {
			/* No edge going out the way we came in! */
			fatal ("edge_right: Bug 2.");
		}
		v2 = (*gridp & GMASK);

		/* If we have stepped to a vertex, get out! */
		if (v2 > 0) break;

		/* Determine which direction to go... */
		switch (code) {
		case 0x05:	break;		/* keep on going... */
		case 0x06:
			edge_up (gridp, v1, i, j, total, cip);
			return;

		case 0x0C:
			edge_down (gridp, v1, i, j, total, cip);
			return;

		default:
			/* Bad type of intersection... */
			fatal ("edge_right: Bug 3.");
		}
	}

	if (v1 < v2) {
		/* We have an edge of the grid-graph! */
		output_edge (v1, v2, total, cip);
	}
}

/*
 * This routine chases the grid-graph edge going UP from position
 * I, J.  We traverse at least one grid element, but stop once we
 * discover a terminal or steiner node.
 */

	static
	void
edge_up (

int32u *	gridp,		/* IN - current grid element */
int		v1,		/* IN - first vertex of edge */
int		i,		/* IN - current Y index */
int		j,		/* IN - current X index */
dist_t		total,		/* IN - total length of edge so far */
struct cinfo *	cip		/* IN - compatibility info */
)
{
int		nt;
int		v2;
int		code;

	nt = grid.nt;
	for (;;) {
		/* Step once upward. */
		if ((i + 1) >= nt) {
			/* Ran off edge of grid! */
			fatal ("edge_up: Bug 1.");
		}
		total += (grid.y_coord [i+1] - grid.y_coord [i]);
		++i;
		gridp += nt;

		code = 0;
		if ((*gridp & RIGHT_EDGE) NE 0) {
			code |= 1;
		}
		if ((*gridp & UP_EDGE) NE 0) {
			code |= 0x02;
		}
		if ((j > 0) AND ((gridp [-1] & RIGHT_EDGE) NE 0)) {
			code |= 0x04;
		}
		if ((i > 0) AND ((gridp [-nt] & UP_EDGE) NE 0)) {
			code |= 0x08;
		}
		if ((code & 0x08) EQ 0) {
			/* No edge going out the way we came in! */
			fatal ("edge_up: Bug 2.");
		}
		v2 = (*gridp & GMASK);

		/* If we have stepped to a vertex, get out! */
		if (v2 > 0) break;

		/* Determine which direction to go... */
		switch (code) {
		case 0x0A:	break;		/* keep on going... */
		case 0x09:
			edge_right (gridp, v1, i, j, total, cip);
			return;

		case 0x0C:
			edge_left (gridp, v1, i, j, total, cip);
			return;

		default:
			/* Bad type of intersection... */
			fatal ("edge_up: Bug 3.");
		}
	}

	if (v1 < v2) {
		/* We have an edge of the grid-graph! */
		output_edge (v1, v2, total, cip);
	}
}

/*
 * This routine chases the grid-graph edge going LEFT from position
 * I, J.  We traverse at least one grid element, but stop once we
 * discover a terminal or steiner node.
 */

	static
	void
edge_left (

int32u *	gridp,		/* IN - current grid element */
int		v1,		/* IN - first vertex of edge */
int		i,		/* IN - current Y index */
int		j,		/* IN - current X index */
dist_t		total,		/* IN - total length of edge so far */
struct cinfo *	cip		/* IN - compatibility info */
)
{
int		nt;
int		v2;
int		code;

	nt = grid.nt;
	for (;;) {
		/* Step once to the left. */
		if (j <= 0) {
			/* Ran off edge of grid! */
			fatal ("edge_left: Bug 1.");
		}
		total += (grid.x_coord [j] - grid.x_coord [j-1]);
		--j;
		--gridp;

		code = 0;
		if ((*gridp & RIGHT_EDGE) NE 0) {
			code |= 1;
		}
		if ((*gridp & UP_EDGE) NE 0) {
			code |= 0x02;
		}
		if ((j > 0) AND ((gridp [-1] & RIGHT_EDGE) NE 0)) {
			code |= 0x04;
		}
		if ((i > 0) AND ((gridp [-nt] & UP_EDGE) NE 0)) {
			code |= 0x08;
		}
		if ((code & 0x01) EQ 0) {
			/* No edge going out the way we came in! */
			fatal ("edge_left: Bug 2.");
		}
		v2 = (*gridp & GMASK);

		/* If we have stepped to a vertex, get out! */
		if (v2 > 0) break;

		/* Determine which direction to go... */
		switch (code) {
		case 0x05:	break;		/* keep on going... */
		case 0x03:
			edge_up (gridp, v1, i, j, total, cip);
			return;

		case 0x09:
			edge_down (gridp, v1, i, j, total, cip);
			return;

		default:
			/* Bad type of intersection... */
			fatal ("edge_left: Bug 3.");
		}
	}

	if (v1 < v2) {
		/* We have an edge of the grid-graph! */
		output_edge (v1, v2, total, cip);
	}
}

/*
 * This routine chases the grid-graph edge going DOWN from position
 * I, J.  We traverse at least one grid element, but stop once we
 * discover a terminal or steiner node.
 */

	static
	void
edge_down (

int32u *	gridp,		/* IN - current grid element */
int		v1,		/* IN - first vertex of edge */
int		i,		/* IN - current Y index */
int		j,		/* IN - current X index */
dist_t		total,		/* IN - total length of edge so far */
struct cinfo *	cip		/* IN - compatibility info */
)
{
int		nt;
int		v2;
int		code;

	nt = grid.nt;
	for (;;) {
		/* Step once downward. */
		if (i <= 0) {
			/* Ran off edge of grid! */
			fatal ("edge_down: Bug 1.");
		}
		total += (grid.y_coord [i] - grid.y_coord [i-1]);
		--i;
		gridp -= nt;

		code = 0;
		if ((*gridp & RIGHT_EDGE) NE 0) {
			code |= 1;
		}
		if ((*gridp & UP_EDGE) NE 0) {
			code |= 0x02;
		}
		if ((j > 0) AND ((gridp [-1] & RIGHT_EDGE) NE 0)) {
			code |= 0x04;
		}
		if ((i > 0) AND ((gridp [-nt] & UP_EDGE) NE 0)) {
			code |= 0x08;
		}
		if ((code & 0x02) EQ 0) {
			/* No edge going out the way we came in! */
			fatal ("edge_down: Bug 2.");
		}
		v2 = (*gridp & GMASK);

		/* If we have stepped to a vertex, get out! */
		if (v2 > 0) break;

		/* Determine which direction to go... */
		switch (code) {
		case 0x0A:	break;		/* keep on going... */
		case 0x03:
			edge_right (gridp, v1, i, j, total, cip);
			return;

		case 0x06:
			edge_left (gridp, v1, i, j, total, cip);
			return;

		default:
			/* Bad type of intersection... */
			fatal ("edge_down: Bug 3.");
		}
	}

	if (v1 < v2) {
		/* We have an edge of the grid-graph! */
		output_edge (v1, v2, total, cip);
	}
}

/*
 * This routine either outputs the edges, or counts them, as appropriate.
 */

	static
	void
output_edge (

int		v1,		/* IN - first vertex of edge */
int		v2,		/* IN - second vertex of edge */
dist_t		total,		/* IN - total length of edge */
struct cinfo *	cip		/* IN - compatibility info */
)
{
char		buf1 [128];

	if (count_edges) {
		++number_of_edges;
	}
	else {
		if (Print_SteinLib_Format) {
			printf ("E ");
		}

		dist_to_string (buf1, total, &(cip -> scale));
		printf ("%d %d %s\n", v1, v2, buf1);
	}
}
