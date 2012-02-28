/***********************************************************************

	File:	utils.c
	Rev:	b-1
	Date:	02/28/2001

	Copyright (c) 1993, 2001 by David M. Warme

************************************************************************

	Various utility routines.

************************************************************************

	Modification Log:

	a-1:	04/18/93	warme
		: Created.  Collected these routines into this file
		:  from other places.
	b-1:	02/28/2001	warme
		: Changes for 3.1 release.
		: Added tracef and struct tracef_control.
		: Added routines for Intel floating point precision fix.
		: Added gst_strdup and store_double.

************************************************************************/

#include "config.h"
#include "steiner.h"
#include <stdarg.h>

#ifdef HAVE_X86_FLOATING_POINT_PRECISION_FIX
#include <fpu_control.h>
#endif


/*
 * Global Routines
 */

void			fatal (char *);
char *			gst_strdup (const char *);
void			init_tables (void);
void *			new (size_t);
void			print_mask (char *, bitmap_t *, int);
void			restore_floating_point_precision (int);
int			save_floating_point_precision (void);
int			set_floating_point_double_precision (void);
void			store_double (double *, double);
void			tracef (const char *, ...);


/*
 * Global Variables
 */

int8u			nbits [256];


/*
 * External References
 */

	/* none */


/*
 * Local Routines
 */

static void		init_nbits (void);

/*
 * This routine displays a fatal message and then dies!
 */

	void
fatal (

char *		msg		/* IN - message to display. */
)
{
	(void) fprintf (stderr, "%s\n", msg);
	(void) fflush (stderr);
	abort ();
}

/*
 * This routine performs all dynamic memory allocation for the program.
 * We test for out of memory condition here.
 */

	void *
new (

size_t		size		/* IN - size of chunk in bytes. */
)
{
void *		p;

	if (size EQ 0) {
		/* Avoid implementation-defined bahavior of malloc! */
		size = 1;
	}

	p = malloc (size);
	if (p EQ NULL) {
		(void) fprintf (stderr, "Out of memory!\n");
		exit (1);
	}

	return (p);
}

/*
 * This routine prints out a given subset.
 */

	void
print_mask (

char *		msg,		/* IN - header message to display */
bitmap_t *	bp1,		/* IN - subset to display */
int		n		/* IN - number of bits in mask */
)
{
int		i;
int		count;

	(void) tracef ("%s", msg);
	count = 0;
	for (i = 0; i < n; i++) {
		if (BITON (bp1, i)) {
			if (++count > 20) {
				tracef ("\n%s	", msg);
				count = 0;
			}
			tracef (" %lu", (int32u) i);
		}
	}
	tracef ("\n");
}

/*
 * This routine initializes various global tables.
 */

	void
init_tables (void)

{
	init_nbits ();
}


/*
 * This routine initializes a translation table that converts an 8-bit
 * byte into the count of ONE bits in the byte.
 */

	static
	void
init_nbits (void)

{
int		i;
int32u		bits;
int		count;

	for (i = 0; i < 256; i++) {
		count = 0;
		for (bits = ((int32u) i); bits > 0; bits >>= 1) {
			if ((bits & 0x01) NE 0) {
				++count;
			}
		}
		nbits [i] = ((int8u) count);
	}
}

/*
 * This routine does CONTROLLABLE printf-ing.
 */

struct tracef_control	tracef_control = { FALSE };


	void
tracef (

const char *	format,		/* IN - format string */
...				/* IN - additional args to format */
)
{
va_list		ap;

	if (tracef_control.disabled) return;

	va_start (ap, format);

	vprintf (format, ap);

	va_end (ap);
}

/*
 * This routine is our own version of the "strdup" function, which
 * is not available everywhere.  This version also uses "new" instead
 * of "malloc".
 */

	char *
gst_strdup (

const char *	s
)
{
size_t		n;
char *		p;

	if (s EQ NULL) return (NULL);

	n = strlen (s) + 1;

	p = new (n);

	strcpy (p, s);

	return (p);
}

/*
 * Special routine to work around problems such as the Intel floating
 * point implementation -- where registers have more precision than
 * memory.  This routine FORCES a value to be stored into memory, from
 * which we can then re-load a value that has exactly the precision of
 * a double, and no more.
 */

	void
store_double (

double *	dp,		/* OUT - double variable to store into */
double		x		/* IN - double value to store */
)
{
	*dp = x;
}

/*
 * On processors that support such a mode (currently only Intel x86):
 * 1 - return the current floating point precision setting,
 * 2 - force the current floating point precision to "double" (returning
 *     the previous setting),
 * 3 - restore a previously sampled floating point precision setting.
 */

/* Somebody should teach these guys how to write header files... */
#define	_FPU_CW_PC	_FPU_EXTENDED

	int
save_floating_point_precision (void)

{
#ifdef HAVE_X86_FLOATING_POINT_PRECISION_FIX
fpu_control_t		cw;

	_FPU_GETCW (cw);
	return (cw);
#else
	return (0);
#endif
}


	int
set_floating_point_double_precision (void)

{
#ifdef HAVE_X86_FLOATING_POINT_PRECISION_FIX
fpu_control_t		orig_cw, cw;

	_FPU_GETCW (orig_cw);
	cw = (orig_cw & ~_FPU_CW_PC) | _FPU_DOUBLE;
	_FPU_SETCW (cw);
	return (orig_cw);
#else
	return (0);
#endif
}


	void
restore_floating_point_precision (

int		prec		/* IN - precision setting to restore */
)
{
#ifdef HAVE_X86_FLOATING_POINT_PRECISION_FIX
fpu_control_t	cw;
	_FPU_GETCW (cw);
	cw = (cw & ~_FPU_CW_PC) | (prec & _FPU_CW_PC);
	_FPU_SETCW (cw);
#endif
}
