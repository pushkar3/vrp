/***********************************************************************

	File:	dsuf.c
	Rev:	a-1
	Date:	11/03/98

	Copyright (c) 1998, 2001 by David M. Warme

************************************************************************

	The "disjoint set union-find" data structure.

************************************************************************

	Modification Log:

	a-1:	11/03/98	warme
		: Created.

************************************************************************/

#include "dsuf.h"
#include "steiner.h"


/*
 * Global Routines
 */

void		dsuf_create (struct dsuf *, int);
int		dsuf_find (struct dsuf *, int);
void		dsuf_makeset (struct dsuf *, int);
void		dsuf_unite (struct dsuf *, int, int);


/*
 * Local Routines
 */

	/* none */

/*
 * This routine creates a collection of N disjoint sets.  They are left
 * uninitialized so that a sparse collection can be accessed quickly.
 */

	void
dsuf_create (

struct dsuf *	dsp,		/* IN/OUT - sets to create */
int		n		/* IN - number of disjoint sets */
)
{
	if (n <= 0) {
		fatal ("dsuf_create: Bug 1.");
	}

	dsp -> set_size		= n;
	dsp -> parent		= NEWA (n, int);
	dsp -> rank		= NEWA (n, int);
}


/*
 * Destroy the given collection of disjoint sets.
 */

	void
dsuf_destroy (

struct dsuf *	dsp		/* IN - sets to destroy */
)
{
	free ((char *) (dsp -> rank));
	free ((char *) (dsp -> parent));

	dsp -> set_size	= 0;
	dsp -> parent	= NULL;
	dsp -> rank	= NULL;
}

/*
 * This routine makes a single disjoint set for item "i".
 */

	void
dsuf_makeset (

struct dsuf *	dsp,		/* IN - collection of sets */
int		i		/* IN - item to make into a disjoint set */
)
{
	if ((i < 0) OR (i >= dsp -> set_size)) {
		/* Item out of bounds. */
		fatal ("dsuf_makeset: Bug 1.");
	}
	dsp -> parent [i]	= i;
	dsp -> rank [i]		= 0;
}

/*
 * This routine "unites" two sets that were previously disjoint.  I and J
 * must be the "canonical" member of each disjoint set (i.e. they must
 * each be the output of a "find" operation), and must be distinct.
 *
 * We perform the "union by rank" heuristic here.
 */

	void
dsuf_unite (

struct dsuf *	dsp,		/* IN - collection of sets */
int		i,		/* IN - first set to unite */
int		j		/* IN - second set to unite */
)
{
int		ri;
int		rj;

	if ((i < 0) OR (i >= dsp -> set_size)) {
		/* Item I is out of range. */
		fatal ("dsuf_unite: Bug 1.");
	}
	if ((j < 0) OR (j >= dsp -> set_size)) {
		/* Item J is out of range. */
		fatal ("dsuf_unite: Bug 2.");
	}
	if (i EQ j) {
		/* Attempt to unite I with I. */
		fatal ("dsuf_unite: Bug 3.");
	}

	ri = dsp -> rank [i];
	rj = dsp -> rank [j];

	if (ri EQ rj) {
		/* Both subtrees have the same maximum depth.  We	*/
		/* arbitrarily choose I to be underneath J.  The rank	*/
		/* of J must then increase.				*/
		dsp -> parent [i] = j;
		dsp -> rank [j]   = rj + 1;
	}
	else if (ri > rj) {
		/* Tree I is (probably) deeper.  Putting J underneath	*/
		/* will not increase I's rank.				*/
		dsp -> parent [j] = i;
	}
	else {
		/* Tree J is (probably) deeper... */
		dsp -> parent [i] = j;
	}
}

/*
 * This routine, given a member I of one of the disjoint sets A, will
 * choose a cannonical member J of set A and return it.  Until set A gets
 * united with some other set, find (I) will always return the same J.
 *
 * This routine performs the "path compression" heuristic.
 */

	int
dsuf_find (

struct dsuf *	dsp,		/* IN - collection of sets */
int		i		/* IN - item to find cannonical item for */
)
{
int		j;
int		k;

	/* Yes, I know this routine is very elegent when coded	*/
	/* recursively...  Here's the iterative version.	*/

	j = dsp -> parent [i];
	if (i EQ j) {
		/* A cannonical element has itself as parent. */
		return (i);
	}

	/* We must search up the tree -- and compress when done... */
	while (TRUE) {
		k = dsp -> parent [j];
		if (j EQ k) break;
		j = k;
	}

	/* Now compress the path (make all items in chain point directly */
	/* at the root K) -- we never have to do this search again!	 */
	while (i NE k) {
		j = dsp -> parent [i];
		dsp -> parent [i] = k;
		i = j;
	}

	return (k);
}
