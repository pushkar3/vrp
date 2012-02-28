/***********************************************************************

	File:	cputime.c
	Rev:	a-1
	Date:	04/16/93

	Copyright (c) 1993, 2001 by David M. Warme

************************************************************************

	This file contains routines for computing CPU time
	consumption for the program.

************************************************************************

	Modification Log:

	a-1:	04/16/93	warme
		: Created.

************************************************************************/

#include "config.h"
#include "steiner.h"

#ifdef UNIX_CPU_TIME
#include <sys/types.h>
#include <sys/time.h>
#include <sys/times.h>
#endif

#include <time.h>


/*
 * Global Routines
 */

cpu_time_t	get_cpu_time (void);
void		convert_cpu_time (cpu_time_t, char *);


/*
 * Define API_CLOCKS_PER_SECOND to be the tick rate presented to the
 * program by the chosen API.
 */

#ifdef UNIX_CPU_TIME
 #define API_CLOCKS_PER_SECOND	CLK_TCK		/* The Posix value */
#else
 #define API_CLOCKS_PER_SECOND	CLOCKS_PER_SEC	/* The ANSI C value */
#endif

/*
 * This routine gets the current amount of CPU time that the program has
 * used.  It is returned in a format we define that is independent of the
 * operating system (i.e. hundreth's of a second).
 *
 * Because the times we are measuring can be quite long (several days),
 * we can run into overflow problems if the operating system gives us
 * a CPU timer with too high a resolution (i.e. microseconds).  This
 * code can be compiled to use either the standard interfaces defined
 * by ANSI C, or the native UNIX interface.  Choose the one that has
 * a resolution closest to our own TICKS_PER_SEC units.
 *
 * The basic conversion formula is:
 *
 *			   sys_ticks * our_ticks_per_sec    1
 *	our_ticks = floor (----------------------------- + ---)
 *			          sys_ticks_per_sec	    2
 *
 * which can be re-written
 *
 *			   2 * sys_ticks * our_tps + sys_tps
 *	our_ticks = floor (---------------------------------)
 *			             2 * sys_tps
 *
 * We compute this using one of several methods, depending upon the
 * relationship between sys_ticks_per_sec and our_ticks_per_sec.
 * The idea of it all is to minimize the chance of arithmetic overflow.
 */

	cpu_time_t
get_cpu_time (void)

{
clock_t		total;
int32u		seconds;
int32u		ticks;
cpu_time_t	cpu_time;

static clock_t	clocks_per_sec = 0;
static int32u	Q, R, method;

#ifdef UNIX_CPU_TIME
	{ struct tms	t;

		times (&t);

		total	= t.tms_utime  + t.tms_stime
			+ t.tms_cutime + t.tms_cstime;
	}
#else
	/* Using the ANSI C defined interface... */
	total = clock ();
#endif

	if (clocks_per_sec EQ 0) {
		/* On some systems this is a system call.  Do only once! */
		clocks_per_sec = API_CLOCKS_PER_SECOND;

		/* Now pick the conversion method and parameters. */
		if (TICKS_PER_SEC >= clocks_per_sec) {
			Q = TICKS_PER_SEC / clocks_per_sec;
			R = TICKS_PER_SEC % clocks_per_sec;
			method = (R EQ 0) ? 0 : 1;
		}
		else {
			R = clocks_per_sec % TICKS_PER_SEC;
			if (R EQ 0) {
				Q = clocks_per_sec / TICKS_PER_SEC;
				method = 2;
			}
			else {
				method = 3;
			}
		}
	}

	seconds = total / clocks_per_sec;
	ticks   = total % clocks_per_sec;

	/* Convert the fractions of a second ticks into .01 second units. */
	switch (method) {
	case 0:
		ticks *= Q;
		break;

	case 1:
		ticks = ticks * Q
			+ (2 * ticks * R + clocks_per_sec) /
			  (2 * clocks_per_sec);
		break;

	case 2:
		ticks = (2 * ticks + Q) / (2 * Q);
		break;

	case 3:
		ticks = ((2 * TICKS_PER_SEC) * ticks + clocks_per_sec)
			/ (2 * clocks_per_sec);
		break;

	default:
		fatal ("get_cpu_time: Bug 1.");
		break;
	}

	cpu_time = seconds * TICKS_PER_SEC + ticks;

	return (cpu_time);
}

/*
 * This routine will convert the given CPU time into a printable
 * null-terminated ASCII string.  The answer contains two decimal places.
 */

	void
convert_cpu_time (

cpu_time_t	time,
char *		out_buf
)
{
cpu_time_t	secs;
cpu_time_t	ticks;
char *		p;
char		buf [20];

#define	ZERO	((cpu_time_t) 0)
#define	TEN	((cpu_time_t) 10)

	secs	= time / TICKS_PER_SEC;
	ticks	= time % TICKS_PER_SEC;

	p = &buf [0];
	if (secs <= ZERO) {
		*p++ = '0';
	}
	else {
		while (secs > ZERO) {
			*p++ = (secs % TEN) + '0';
			secs /= TEN;
		}
	}
	while (p > &buf [0]) {
		*out_buf++ = *--p;
	}
	*out_buf++ = '.';
	ticks *= TEN;
	*out_buf++ = (ticks / TICKS_PER_SEC) + '0';
	ticks %= TICKS_PER_SEC;
	ticks *= TEN;
	*out_buf++ = (ticks / TICKS_PER_SEC) + '0';
	*out_buf++ = '\0';

#undef ZERO
#undef TEN
}
