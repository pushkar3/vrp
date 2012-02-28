/***********************************************************************

	File:	emptyr.c
	Rev:	b-1
	Date:	01/31/2000

	Copyright (c) 1998, 2001 by David M. Warme and Martin Zachariasen

************************************************************************

	Routines for efficiently determining whether or not two
	terminals define an empty rectangle.  We precompute this
	information and store it compactly.

************************************************************************

	Modification Log:

	a-1:	09/28/98	warme
		: Created.  Implemented Zachariasen's algorithm
		:  using Warme's infrastructure.
	b-1:	01/31/2000	martinz
		: Successor array is computed within initialization
		: if given as a NULL pointer.

************************************************************************/

#include "emptyr.h"
#include "steiner.h"


/*
 * Global Routines
 */

int		count_empty_rectangles (bitmap_t *, int);
bitmap_t *	init_empty_rectangles (struct pset *, int *);
bool		is_empty_rectangle (bitmap_t *, int, int);
void		shutdown_empty_rectangles (bitmap_t *);


/*
 * Local Macros
 */

#define _POS(i,j)	(((i * (i-1)) >> 1) + j)
#define POS(i,j)	((i >= j) ? _POS(i,j) : _POS(j,i))


/*
 * Local Routines
 */

static void		set_bit (bitmap_t *, int, int);
static int *		heapsort_x (struct pset *);

/*
 * Pre-compute an N by N boolean matrix whose (i,j)-th element is
 * TRUE if-and-only-if the interior of the rectangle defined by
 * terminals i and j is devoid of terminals.
 *
 * Since this matrix is symmetric and the diagonal elements are
 * all TRUE, we store only the lower triangle in a compact bit-vector
 * form.
 *
 * The succ0 array may be given as the NULL pointer in which case it is
 * computed by the the procedure.
 */

	bitmap_t *
init_empty_rectangles (

struct pset *		pts,	/* IN - point set to use */
int *			succ0	/* IN - next "dir 0" successor for each term */
)
{
int		i;
int		j;
int		n;
bitmap_t *	bits;
int		nbits;
int		nwords;
int *		x_order;
int *		succ;
struct point *	p;
struct point *	q;
double		dx, dy;
double		x, y;
double		top_dist, bot_dist;
double		old_top_dist, old_bot_dist;
double		top_x, bot_x;

	n = pts -> n;

	/* Do we need to compute succ array? */
	if (succ0 EQ NULL) {
		x_order = heapsort_x (pts);
		succ	= NEWA (n, int);

		for (i = 1; i < n; i++) {
			succ [x_order [i - 1]] = x_order [i];
		}
		succ [x_order [i - 1]] = -1;
	}
	else {
		succ = succ0;
	}

	/* Create the lower-triangular bit-vector and zero it. */
	nbits = (n * (n - 1)) >> 1;
	nwords = BMAP_ELTS (nbits);

	bits = NEWA (nwords, bitmap_t);
	for (i = 0; i < nwords; i++) {
		bits [i] = 0;
	}

	p = &(pts -> a [0]);
	for (i = 0; i < n; i++, p++) {
		x = p -> x;
		y = p -> y;

		top_dist	= INF_DISTANCE;
		bot_dist	= INF_DISTANCE;
		old_top_dist	= INF_DISTANCE;
		old_bot_dist	= INF_DISTANCE;
		top_x		= x;
		bot_x		= x;
		for (j = succ [i]; j >= 0; j = succ [j]) {
			q = &(pts -> a [j]);
			dx = q -> x - x;
			if (dx EQ 0.0) {
				/* Q is exactly on vertical line through P. */
				set_bit (bits, i, j);
				continue;
			}

			dy = q -> y - y;
			if (dy EQ 0.0) {
				/* Q is exactly on horiz line through Q. */
				set_bit (bits, i, j);
				continue;
			}

			if (dy > 0.0) {
				/* Q is on top (above P). */
				if (dy <= top_dist) {
					set_bit (bits, i, j);
					if (q -> x > top_x) {
						old_top_dist = top_dist;
						top_x = q -> x;
					}
					top_dist = dy;
				}
				else if ((q -> x EQ top_x) AND
					 (dy <= old_top_dist)) {
					set_bit (bits, i, j);
				}
			}
			else {
				/* Q is on bottom (below P). */
				dy = - dy;
				if (dy <= bot_dist) {
					set_bit (bits, i, j);
					if (q -> x > bot_x) {
						old_bot_dist = bot_dist;
						bot_x = q -> x;
					}
					bot_dist = dy;
				}
				else if ((q -> x EQ bot_x) AND
					 (dy <= old_bot_dist)) {
					set_bit (bits, i, j);
				}
			}
		}
	}

	if (succ0 EQ NULL) {
		free ((char *) x_order);
		free ((char *) succ);
	}

	return (bits);
}

/*
 * Set bit (i,j) in the bit matrix.
 */

	static
	void
set_bit (

bitmap_t *		bits,	/* IN - bit vector for matrix */
int			i,	/* IN - row number */
int			j	/* IN - column number */
)
{
int			pos;

	/* The diagonal is not stored! */
	if (i EQ j) return;

	pos = POS (i, j);
	SETBIT (bits, pos);
}

/*
 * Clean up the empty rectangles data structure.
 */

	void
shutdown_empty_rectangles (

bitmap_t *		bits	/* IN - empty rectangle bit matrix */
)
{
	free ((char *) bits);
}

/*
 * Determine if the rectangle defined by terminals i and j is empty.
 */

	bool
is_empty_rectangle (

bitmap_t *		bits,	/* IN - empty rectangle bit matrix */
int			i,	/* IN - first terminal number */
int			j	/* IN - second terminal number */
)
{
int			pos;

	if (i EQ j) {
		/* The rectangle is zero by zero, and therefore has	*/
		/* no interior.	 There cannot be any terminals inside!	*/
		return (TRUE);
	}

	pos = POS (i, j);

	return (BITON (bits, pos));
}

/*
 * This routine counts the total number of 1 bits in the bit matrix.
 * This represents the number of pairs of terminals that define
 * rectangles whose INTERIORS are devoid of terminals.
 */

	int
count_empty_rectangles (

bitmap_t *		bits,	/* IN - empty rectangle bit matrix */
int			n	/* IN - number of terminals */
)
{
int		i;
int		nbits;
int		nwords;
int		count;
bitmap_t	mask1;
bitmap_t	mask2;
bitmap_t	limit;

	nbits = (n * (n - 1)) >> 1;

	nwords = nbits / BPW;
	nbits  = nbits % BPW;

	count = 0;
	for (i = 0; i < nwords; i++) {
		mask1 = bits [i];
		while (mask1 NE 0) {
			mask1 ^= (mask1 & -mask1);
			++count;
		}
	}

	if (nbits > 0) {
		/* Count straggling bits in last word. */
		mask1 = bits [nwords];
		limit = 1 << nbits;	/* first bit to EXCLUDE from count */

		while (mask1 NE 0) {
			mask2 = (mask1 & -mask1);
			if (mask2 >= limit) break;
			mask1 ^= mask2;
			++count;
		}
	}

	return (count);
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
