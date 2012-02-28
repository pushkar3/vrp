/***********************************************************************

	File:	steiner.h
	Rev:	b-1
	Date:	02/28/2001

	Copyright (c) 1993, 2001 by David M. Warme

************************************************************************

	General declarations for the Steiner Tree program.

************************************************************************

	Modification Log:

	a-1:	02/20/93	warme
		: Created.
	a-2:	08/31/98	warme
		: Split off metric dependent stuff.
	b-1:	02/28/2001	warme
		: Changes for 3.1 release.
		: Modified "struct numlist" to be partially scaled
		:  numeric representation instead of textual.
		: Added struct scale_info and UNSCALE macro.
		: Added tracef routine and struct tracef_control.
		: Pass new scale_info to all routines that need it.
		: Added several new routines.

************************************************************************/

#ifndef	STEINER_H
#define	STEINER_H

#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*
 * Typedef's to insulate us from the ugliness of C types.
 */

typedef unsigned char		int8u;
#if	defined(__STDC__) && ((__STDC__+0)==1)
typedef	signed char		int8s;
#else
typedef char			int8s;
#endif
typedef unsigned short		int16u;
typedef short			int16s;
typedef unsigned long		int32u;
typedef long			int32s;
typedef char			bool;


/*
 * Macros to protect us from some of C's sharpest edges.
 */

#define NOT	!
#define AND	&&
#define OR	||
#define FALSE	0
#define	TRUE	1
#define EQ	==
#define NE	!=

/*
 * Type used to represent CPU time usage.
 */

typedef int32u		cpu_time_t;


/*
 * Primitive types associated with points.
 */

typedef double		coord_t;	/* type of X and Y coordinates. */
typedef double		dist_t;		/* distance. */
typedef	int		pnum_t;		/* type of point number indices. */

#define	INF_DISTANCE	((dist_t) DBL_MAX)	/* "infinite" distance. */


/*
 * Bit Maps
 */

typedef	int32u		bitmap_t;	/* Element of Bit-map vector. */
#define	BPW		32		/* Bits per word of a bit-map. */
#define	BMAP_ELTS(n)	(((n) + (BPW-1)) / BPW)

/*
 * The folling represents a single edge in a weighted tree or graph.
 */

struct edge {
	dist_t		len;		/* Length of edge. */
	pnum_t		p1;		/* First endpoint of edge. */
	pnum_t		p2;		/* Second endpoint of edge. */
};


/*
 * A linked list structure, each node of which holds a single real
 * number in special "scaled" form.  We use these during input so that we
 * can examine an entire set of numbers before we commit to an internal
 * representation for them.
 */

struct numlist {
	double			mantissa;	/* number scaled to an */
						/* integer value */
	int			expon;		/* val = mantissa * 10^expon */
	int			nsig;		/* # signif digits in mant */
	struct numlist *	next;		/* next in linked list */
};

/*
 * A structure that describes how external coordinates, lengths, etc.
 * are scaled internally.  We do this to eliminate some of the bad
 * behavior of floating point representation.
 *
 * The internal values are actually 10**scale times BIGGER than the
 * external values read in.  Since scale can be negative, we have two
 * distinct conversion factors -- one to multiply by, and one to divide
 * by.
 */

struct scale_info {
	int		scale;		/* internal vals mpy'd by 10**scale */
	int		min_precision;	/* min decimal places to output */
	double		scale_mul;	/* 1 or 10**-scale, if scale < 0 */
	double		scale_div;	/* 1 or 10**scale, if scale >= 0 */
};


/*
 * Macro to unscale an internal value back to external form.
 */

#define UNSCALE(val,p)	((val) * (p) -> scale_mul / (p) -> scale_div)

/*
 * The various representations for geometric objects:
 *
 * The following represents a single point.
 */

struct point {
	coord_t		x;
	coord_t		y;
	pnum_t		pnum;
};


/*
 * The following represents an input set of points, with the count.
 * We always allocate these dynamically with an appropriate number
 * of points in the "a" array.
 */

struct pset {
	int		n;
	struct point	a [1];
};


/*
 * The structure used to represent a single FST.
 */

struct full_set {
	struct full_set *	next;
	int			tree_num;
	dist_t			tree_len;
	struct pset *		terminals;
	struct pset *		steiners;
	int			nedges;
	struct edge *		edges;
};

/*
 * The "cinfo" describes the basic problem instance, and also contains
 * some additional information (i.e., compatibility/incompatibility info).
 */

struct cinfo {
	/* The basic problem description. */
	int			num_verts;
	int			num_edges;
	int			num_vert_masks;
	int			num_edge_masks;
	int **			edge;		/* Terminals in each edge */
	int *			edge_size;	/* Cardinality of each edge */
	dist_t *		cost;		/* Length of each edge */
	bool *			tflag;		/* TRUE: vertex is terminal, */
						/* FALSE: vertex is Steiner. */
	int			metric;		/* Rectilinear, Euclid, etc. */

	struct scale_info	scale;		/* problem scaling info */

	/* Minimum length difference between two SMTs of different lengths */
	dist_t			integrality_delta;
	/* Initial problem bit-masks. */
	bitmap_t *		initial_vert_mask;
	bitmap_t *		initial_edge_mask;
	bitmap_t *		required_edges;
	/* pre-initialized tables. */
	int **			term_trees;
	int **			inc_edges;

	char *			description;
	cpu_time_t		p1time;

	/* Optional geometric information. */
	struct pset *		pts;
	struct full_set **	full_trees;
	dist_t			mst_length;
};


/*
 * The following equates define the permissible distance metrics.
 */

#define	RECTILINEAR		1
#define	EUCLIDEAN		2
#define	PURE_GRAPH		3

/*
 * Useful macros.
 */

#define DELTAX(p1,p2)	(fabs((p1) -> x - (p2) -> x))
#define	DELTAY(p1,p2)	(fabs((p1) -> y - (p2) -> y))

#define	RDIST(p1,p2)	(DELTAX ((p1), (p2)) + DELTAY ((p1), (p2)))
#define	EDIST(p1,p2)	(hypot (DELTAX ((p1), (p2)), DELTAY ((p1), (p2))))


#define	NEW(type)	((type *) new (sizeof (type)))
#define	NEWA(n, type)	((type *) new ((size_t) ((n) * sizeof (type))))


#define	NULL_PSET	((struct pset *) 0)
#define	PSET_SIZE(n)	(offsetof (struct pset, a [n]))
#define NEW_PSET(n)	((struct pset *) new (PSET_SIZE (n)))
#define	COPY_PSET(dp, sp)	(void) memcpy ((dp), (sp), PSET_SIZE ((sp) -> n))
#define	ZERO_PSET(dp, n)	(void) memset ((dp), 0, PSET_SIZE (n))

#define	SETBIT(bm, n)	((bm) [(n) / BPW] |= (1ul << ((n) % BPW)))
#define	CLRBIT(bm, n)	((bm) [(n) / BPW] &= ~(1ul << ((n) % BPW)))
#define	BITON(bm, n)	(((bm) [(n) / BPW] & (1ul << ((n) % BPW))) NE 0)
#define	NBITSON(m)	(  nbits [ (m)        & 0xFFlu]		\
			 + nbits [((m) >>  8) & 0xFFlu]		\
			 + nbits [((m) >> 16) & 0xFFlu]		\
			 + nbits [((m) >> 24) & 0xFFlu])

/*
 * Stuff for measuring CPU time usage.
 */

#define	TICKS_PER_SEC	100		/* This is the units WE use! */

/*
 * Global Variables
 */

extern int8u		nbits [];


/*
 * A structure for controlling the output of tracef().
 */

struct tracef_control {
	bool		disabled;
};

extern struct tracef_control	tracef_control;


/*
 * Function Prototypes.
 */

extern int		compute_scaling_factor (struct numlist *);
extern void		convert_cpu_time (cpu_time_t, char *);
extern void		coord_to_string (char *, coord_t, struct scale_info *);
extern bool		decode_cpu_time_limit (char *, int32u *);
extern void		dist_to_string (char *, dist_t, struct scale_info *);
extern int		euclidean_mst (struct pset *, struct edge *);
extern dist_t		euclidean_mst_length (struct pset *);
extern void		fatal (char *);
extern cpu_time_t	get_cpu_time (void);
extern struct pset *	get_points (FILE *, struct scale_info *);
extern char *		gst_strdup (const char *);
extern void		init_output_conversion (struct pset *,
						int,
						struct scale_info *);
extern void		init_tables (void);
extern int		kahng_robins (struct pset *, dist_t, struct edge *);
extern dist_t		kahng_robins_length (struct pset *, dist_t);
extern int		mst_edge_list (int, int, struct edge *, struct edge *);
extern void *		new (size_t);
extern struct numlist *	parse_line_of_numbers (FILE *);
extern void		print_mask (char *, bitmap_t *, int);
extern coord_t		read_numlist (struct numlist *, struct scale_info *);
extern int		rect_mst (struct pset *, struct edge *, bitmap_t *);
extern dist_t		rect_mst_length (struct pset *);
extern void		restore_floating_point_precision (int);
extern int		save_floating_point_precision (void);
extern int		set_floating_point_double_precision (void);
extern void		set_scale_info (struct scale_info *, int);
extern void		sort_edge_list (struct edge *, int);
extern void		sort_ints (int *, int);
extern void		start_limiting_cpu_time (int32u);
extern void		store_double (double *, double);
extern void		tracef (const char *, ...);

#endif
