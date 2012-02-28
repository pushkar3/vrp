/***********************************************************************

	File:	p1write.c
	Rev:	b-2
	Date:	02/28/2001

	Copyright (c) 1993, 2001 by David M. Warme

************************************************************************

	Routines for writing the data that is output from phase 1.

************************************************************************

	Modification Log:

	b-1:	01/12/97	warme
		: Split p1io.c into reading and writing parts.
	b-2:	02/28/2001	warme
		: Changes for 3.1 release.
		: Use new scaling stuff.
		: Terminate last tflags line.
		: Don't include pruned FSTs in incompatible lists.
		: Avoid empty line (a single tab) when incompatible
		:  list is empty.

************************************************************************/

#include "config.h"
#include "p1io.h"
#include "steiner.h"


/*
 * Global Routines
 */

void		print_phase_1_data (struct cinfo *, int);


/*
 * External References
 */

extern char * get_machine_string (void);


/*
 * Local Routines
 */

static void		double_to_hex (double, char *);
static void		print_version_0 (struct cinfo *, int);
static void		print_version_2 (struct cinfo *, int);

/*
 * This routine prints out all of the data that is output from
 * phase 1 of the algorithm.  The output is almost readable by
 * humans, but is designed more to be machine-readable.
 */

	void
print_phase_1_data (

struct cinfo *		cip,		/* IN - compatability info. */
int			version		/* IN - version to generate. */
)
{
	switch (version) {
	case P1IO_VERSION_0:
		print_version_0 (cip, version);
		break;

	case P1IO_VERSION_2:
	case P1IO_VERSION_3:
		/* Versions 2 and 3 are different enough that we use	*/
		/* a completely different routine. */
		print_version_2 (cip, version);
		break;

	default:
		fatal ("print_phase_1_data: Bug 1.");
	}
}

/*
 * This routine prints out the extended OR-library format -- version 0.
 */

	static
	void
print_version_0 (

struct cinfo *		cip,		/* IN - compatibility info. */
int			version		/* IN - version to generate. */
)
{
int			i;
int			j;
int			n;
int			m;
int *			vp1;
int *			vp2;
char			buf1 [64];

	n = cip -> num_verts;
	m = cip -> num_edges;

	printf (" %d %d\n", n, m);

	/* hyperedges... */
	for (i = 0; i < m; i++) {
		vp1 = cip -> edge [i];
		vp2 = cip -> edge [i + 1];
		while (vp1 < vp2) {
			j = *vp1++;
			printf (" %d", j + 1);
		}

		dist_to_string (buf1, cip -> cost [i], &(cip -> scale));
		printf (" %s\n", buf1);
	}

	/* Info about terminals. */
	j = 0;
	for (i = 0; i < n; i++) {
		if (NOT (cip -> tflag [i])) continue;
		++j;
	}
	printf ("\t%d\n", j);
	j = 0;
	for (i = 0; i < n; i++) {
		if (NOT (cip -> tflag [i])) continue;
		printf (" %d", i + 1);
		if (++j >= 10) {
			printf ("\n");
			j = 0;
		}
	}
}

/*
 * This routine prints out the new data format -- versions 2 and 3.
 */

	static
	void
print_version_2 (

struct cinfo *		cip,		/* IN - compatibility info. */
int			version		/* IN - version to generate. */
)
{
int			i;
int			j;
int			k;
int			n;
int			m;
int			col;
int			kmasks;
int			nmasks;
int			count;
dist_t			mst_len;
struct point *		p1;
struct full_set *	fsp;
struct pset *		terms;
struct pset *		steins;
int *			ep1;
int *			ep2;
int *			ep3;
int *			flist;
int *			vp1;
int *			vp2;
bitmap_t *		tmask;
bool			geometric;
char *			desc;
char			buf1 [64];
char			buf2 [64];
char			buf3 [64];
char			buf4 [64];

#define	MAXCOL	78

	geometric = ((cip -> metric NE PURE_GRAPH) AND
		     (cip -> full_trees NE NULL) AND
		     (cip -> pts NE NULL));

	n = cip -> num_verts;
	m = cip -> num_edges;
	kmasks = cip -> num_vert_masks;
	nmasks = cip -> num_edge_masks;

	printf ("V%d\n", version);
	desc = cip -> description;
	printf ("%s\n", (desc NE NULL) ? desc : "");
	if (geometric) {
		printf ("%d\n", cip -> metric);
	}
	else {
		/* Not enough info to put out geometric stuff -- pure graph */
		printf ("%d\n", PURE_GRAPH);
	}
	printf ("%d\n", n);

	if (geometric) {
		/* Compute MST length... */
		mst_len = cip -> mst_length;
		dist_to_string (buf1, mst_len, &(cip -> scale));
		double_to_hex ((double) mst_len, buf2);
		printf ("%s %s\n", buf1, buf2);
	}

	if ((version <= P1IO_VERSION_2) AND geometric) {
		/* No duplicate terminal groups... */
		printf ("0\n");
	}
	printf ("%d\n", cip -> scale.scale);
	if (version >= P1IO_VERSION_3) {
		/* Version 3 has integrality delta... */
		dist_to_string (buf1,
				cip -> integrality_delta,
				&(cip -> scale));
		double_to_hex ((double) cip -> integrality_delta, buf2);
		printf ("%s %s\n", buf1, buf2);
	}
	printf ("%s\n", get_machine_string ());
	printf ("%lu\n", cip -> p1time);	/* CPU time */
	printf ("%d\n", m);			/* Number of hyperedges */

	if (geometric) {
		p1 = &(cip -> pts -> a [0]);
		for (i = 0; i < n; i++, p1++) {
			coord_to_string (buf1, p1 -> x, &(cip -> scale));
			coord_to_string (buf2, p1 -> y, &(cip -> scale));
			double_to_hex ((double) (p1 -> x), buf3);
			double_to_hex ((double) (p1 -> y), buf4);
			printf ("\t%s\t%s\t%s\t%s\n", buf1, buf2, buf3, buf4);
		}
	}

	if (version >= P1IO_VERSION_3) {
		/* Print the terminal/Steiner flag for each vertex. */
		j = 0;
		for (i = 0; i < n; i++) {
			if (j EQ 0) {
				printf ("\t");
			}
			printf (" %d", cip -> tflag [i]);
			if (++j >= 10) {
				printf ("\n");
				j = 0;
			}
		}
		if (j > 0) {
			printf ("\n");
		}
	}

	flist = NEWA (m, int);
	tmask = NEWA (kmasks, bitmap_t);
	for (i = 0; i < kmasks; i++) {
		tmask [i] = 0;
	}

	/* hyperedges... */
	for (i = 0; i < m; i++) {
		printf ("\t%d\n", cip -> edge_size [i]);
		printf ("\t");
		vp1 = cip -> edge [i];
		vp2 = cip -> edge [i + 1];
		col = 8;
		while (vp1 < vp2) {
			j = *vp1++;
			sprintf (buf1, "%d", j + 1);
			k = strlen (buf1);
			if (col + 1 + k >= MAXCOL) {
				printf ("\n\t\t%s", buf1);
				col = 16 + k;
			}
			else if (col <= 8) {
				printf ("%s", buf1);
				col += k;
			}
			else {
				printf (" %s", buf1);
				col += (1 + k);
			}
		}
		printf ("\n");

		dist_to_string (buf1, cip -> cost [i], &(cip -> scale));
		double_to_hex ((double) (cip -> cost [i]), buf2);
		printf ("\t%s\t%s\n", buf1, buf2);

		if (geometric) {
			fsp = cip -> full_trees [i];
			terms = fsp -> terminals;
			steins = fsp -> steiners;

			if (steins EQ NULL) {
				printf ("\t0\n");
			}
			else {
				printf ("\t%d\n", steins -> n);
				for (j = 0; j < steins -> n; j++) {
					p1 = &(steins -> a [j]);
					coord_to_string (buf1,
							 p1 -> x,
							 &(cip -> scale));
					coord_to_string (buf2,
							 p1 -> y,
							 &(cip -> scale));
					double_to_hex ((double) (p1 -> x), buf3);
					double_to_hex ((double) (p1 -> y), buf4);
					printf ("\t%s\t%s\t%s\t%s\n", buf1, buf2, buf3, buf4);
				}
			}

			printf ("\t%d\n", fsp -> nedges);
			for (j = 0; j < fsp -> nedges; j++) {
				k = fsp -> edges [j].p1;
				printf ("\t\t%d",
					(k < terms -> n)
					    ? k + 1
					    : terms -> n - k - 1);
				k = fsp -> edges [j].p2;
				printf ("\t%d\n",
					(k < terms -> n)
					    ? k + 1
					    : terms -> n - k - 1);
			}
		}

		if ((cip -> initial_edge_mask NE NULL) AND
		    (NOT BITON (cip -> initial_edge_mask, i))) {
			printf ("\t0\n");	/* edge never needed */
		}
		else if ((cip -> required_edges NE NULL) AND
			 (BITON (cip -> required_edges, i))) {
			printf ("\t2\n");	/* edge always needed */
		}
		else {
			printf ("\t1\n");	/* edge sometimes needed */
		}

		/* List the hyperedges that are incompatible with this one. */
		if (cip -> inc_edges EQ NULL) {
			/* No incompatibility info to give! */
			printf ("\t0\n");
		}
		else if ((cip -> initial_edge_mask NE NULL) AND
			 (NOT BITON (cip -> initial_edge_mask, i))) {
			/* This hyperedge was pruned. */
			/* Do not list all others as incompatible!!! */
			printf ("\t0\n");
		}
		else {
			/* Gather the list of incompatible edges we	*/
			/* choose to mention.  (i.e., eliminate those	*/
			/* that are "basic" incompatibilities.)		*/
			vp1 = cip -> edge [i];
			vp2 = cip -> edge [i + 1];
			while (vp1 < vp2) {
				j = *vp1++;
				SETBIT (tmask, j);
			}
			ep1 = flist;
			ep2 = cip -> inc_edges [i];
			ep3 = cip -> inc_edges [i + 1];
			while (ep2 < ep3) {
				j = *ep2++;
				if (i EQ j) continue;
				if (NOT BITON (cip -> initial_edge_mask, j)) continue;
				count = 0;
				vp1 = cip -> edge [j];
				vp2 = cip -> edge [j + 1];
				while (vp1 < vp2) {
					k = *vp1++;
					if (BITON (tmask, k)) {
						++count;
					}
				}
				if (count >= 2) {
					/* FST's i and j have >= 2 terms */
					/* in common -- don't list this */
					/* obvious incompatibility... */
					continue;
				}
				*ep1++ = j;
			}

			printf ("\t%d\n", ep1 - flist);
			if (ep1 > flist) {
				printf ("\t");
				col = 8;
				for (ep2 = flist; ep2 < ep1; ep2++) {
					j = *ep2;
					sprintf (buf1, "%d", j + 1);
					k = strlen (buf1);
					if (col + 1 + k >= MAXCOL) {
						printf ("\n\t\t%s", buf1);
						col = 16 + k;
					}
					else if (col <= 8) {
						printf ("%s", buf1);
						col += (1 + k);
					}
					else {
						printf (" %s", buf1);
						col += (1 + k);
					}
				}
				printf ("\n");
			}
			vp1 = cip -> edge [i];
			vp2 = cip -> edge [i + 1];
			while (vp1 < vp2) {
				j = *vp1++;
				CLRBIT (tmask, j);
			}
		}

		if (version <= P1IO_VERSION_2) {
			/* Number of strongly compatible full sets... */
			printf ("\t0\n");
		}
	}

	free ((char *) tmask);
	free ((char *) flist);

#undef MAXCOL
}

/*
 * This routine converts a double into a printable ASCII string
 * that represents the exact numeric value in hexidecimal.  The
 * printable string has the following format:
 *
 *	[-].{hexdigits}[x[-]{hexdigits}]
 *
 * The leading minus sign is optional and the exponent part will be
 * left off if it is zero.
 */

	static
	void
double_to_hex (

double		value,		/* IN - the floating point value to encode */
char *		s		/* OUT - the hex ASCII string */
)
{
double		mant;
int		expon;
int		digit;
char *		p;
char		buf [12];

static char	hex [] = "0123456789ABCDEF";

	if (value < 0) {
		*s++ = '-';
		value = - value;
	}

	*s++ = '.';

	mant = frexp (value, &expon);

	if (mant EQ 0.0) {
		*s++ = '0';
	}
	else {
		while (mant > 0.0) {
			mant *= 16.0;
			digit = (int) mant;
			*s++ = hex [digit];
			mant -= ((double) digit);
		}
		if (expon NE 0) {
			*s++ = 'x';
			if (expon < 0) {
				*s++ = '-';
				expon = - expon;
			}
			p = &buf [0];
			while (expon > 0) {
				*p++ = hex [expon & 0x0F];
				expon >>= 4;
			}
			while (p > &buf [0]) {
				*s++ = *--p;
			}
		}
	}

	*s++ = '\0';
}
