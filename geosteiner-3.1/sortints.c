/***********************************************************************

	File:	sortints.c
	Rev:	a-1
	Date:	01/05/2001

	Copyright (c) 1993, 2001 by David M. Warme

************************************************************************

	Routine to efficiently sort an array of integers.

************************************************************************

	Modification Log:

	a-1:	01/04/2001	warme
		: Created.

************************************************************************/

#include "steiner.h"


/*
 * Global Routines
 */

void			sort_ints (int *, int);

/*
 * Sort an array of integers into increasing order.
 * We assume that the array can be either very small, very large,
 * or somewhere in between.  We also want to efficiently handle the
 * case where the initial array is already sorted.
 *
 * We make a single pass of the bubble sort, remembering the position
 * (K-1,K) of the last exchange.  If K <= 1 upon completion of this scan,
 * then there we no exchanges made and the array is already sorted.
 * Otherwise, items (K,...,N-1) are now fully sorted, so it suffices to
 * sort the first N = K array elements.  If K <= 16, then finish up with
 * the bubble sort -- otherwise, use heapsort.
 */

	void
sort_ints (

int *		array,		/* IN/OUT - the array to sort */
int		n		/* IN - number of ints in the array */
)
{
int		i;
int		j;
int		k;
int		ai;
int		aj;
int		ai1;

	k = 0;
	for (i = 1; i < n; i++) {
		ai1 = array [i - 1];
		ai = array [i];
		if (ai1 > ai) {
			array [i - 1] = ai;
			array [i] = ai1;
			k = i;
		}
	}

	if (k <= 1) {
		/* Array is now completely sorted. */
		return;
	}

	n = k;

	if (n <= 16) {
		/* Remaining problem is pretty small -- polish it off	*/
		/* using bubble sort...					*/
		do {
			k = 0;
			for (i = 1; i < n; i++) {
				ai1 = array [i - 1];
				ai  = array [i];
				if (ai1 > ai) {
					array [i - 1] = ai;
					array [i] = ai1;
					k = i;
				}
			}
			n = k;
		} while (n > 1);
		return;
	}

	/* Problem is bigger.  Use heapsort, which always uses	*/
	/* O(N * LOG(N)) time -- no matter what.		*/

	/* Construct the heap via sift-downs, in O(n) time.	*/
	
	for (k = n >> 1; k >= 0; k--) {
		j = k;
		for (;;) {
			i = (j << 1) + 1;
			if (i >= n) {
				/* Hit bottom of heap, sift-down is done. */
				break;
			}
			ai = array [i];
			if (i + 1 < n) {
				/* Increment i (to right subchild of j) */
				/* if the right subchild is greater. */
				ai1 = array [i + 1];
				if (ai1 > ai) {
					++i;
					ai = ai1;
				}
			}
			aj = array [j];
			if (ai <= aj) {
				/* Greatest child isn't bigger.  Sift-	*/
				/* down is done.			*/
				break;
			}
			/* Sift down and continue. */
			array [j] = ai;
			array [i] = aj;
			j = i;
		}
	}

	/* Now do actual sorting.  Exchange first/last and sift down. */
	while (n > 1) {
		/* Largest is at array [0], swap with array [n - 1],	*/
		/* thereby putting it into final position.		*/
		--n;
		i = array [0];
		array [0] = array [n];
		array [n] = i;

		/* Now restore the heap by sifting array [0] down. */
		j = 0;
		for (;;) {
			i = (j << 1) + 1;
			if (i >= n) {
				/* Hit bottom of heap, sift-down is done. */
				break;
			}
			ai = array [i];
			if (i + 1 < n) {
				/* Increment i (to right subchild of j) */
				/* if the right subchild is greater. */
				ai1 = array [i + 1];
				if (ai1 > ai) {
					++i;
					ai = ai1;
				}
			}
			aj = array [j];
			if (ai <= aj) {
				/* Greatest child isn't bigger.  Sift-	*/
				/* down is done. */
				break;
			}
			/* Sift down and continue. */
			array [j] = ai;
			array [i] = aj;
			j = i;
		}
	}
}
