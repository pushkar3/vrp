/***********************************************************************

	File:	cpulimit.c
	Rev:	b-1
	Date:	01/10/97

	Copyright (c) 1993, 2001 by David M. Warme

************************************************************************

	Routines to enforce CPU time limits.

************************************************************************

	Modification Log:

	b-1:	01/10/97	warme
		: Split off from old_bs.c

************************************************************************/

#include "config.h"
#include "steiner.h"
#include <signal.h>


/*
 * Global Routines
 */

bool		decode_cpu_time_limit (char *, int32u *);
void		start_limiting_cpu_time (int32u);


/*
 * External References
 */

	/* none */


/*
 * Local Routines
 */

static RETSIGTYPE	check_limit (int);

static int32u		cpu_time_limit;
static int32u		end_time;

/*
 * This routine decodes an argument string that represents an amount of
 * CPU time.  The time is the sum of a sequence of time groups.  Each
 * time group is a sequence of decimal digits followed by an optional
 * letter defining units of minutes, hours, days, etc.
 * The regular expression accepted is:
 *
 * ([0-9]+([dhms][a-z]*)?)*
 */

	bool
decode_cpu_time_limit (

char *		ap,
int32u *	limit
)
{
int32u		total;
int32u		group;
char		c;

	*limit = 0;
	total = 0;
	c = *ap++;
	while (c NE '\0') {
		if ((c < '0') OR (c > '9')) {
			return (FALSE);
		}
		group = c - '0';
		while (((c = *ap++) NE '\0') AND ('0' <= c) AND (c <= '9')) {
			group = (10 * group) + (c - '0');
		}
		switch (c) {
		case 'd':	group *= (24 * 60 * 60);	break;
		case 'h':	group *= (60 * 60);		break;
		case 'm':	group *= 60;			break;
		case 's':					break;
		case '\0':					break;
		default:	return (FALSE);
		}
		/* strip rest of unit name... */
		while ((c > 0) AND ((c < '0') OR (c > '9'))) {
			c = *ap++;
		}
		total += group;
	}

	*limit = total;

	return (TRUE);

}

/*
 * This routine begins limiting CPU time starting right now.  If the given
 * limit is ever exceeded, the process exits with an error message.
 */

	void
start_limiting_cpu_time (

int32u		limit		/* IN - maximum CPU seconds to permit. */
)
{
int32u		now;
int32u		then;

	if (limit > 0) {
		cpu_time_limit = limit;

		now = get_cpu_time ();
		then = now + limit * TICKS_PER_SEC;

		/* Record the time at which we should shut down... */
		end_time = then;

		/* Schedule the first timer interrupt. */
		check_limit (0);
	}
}

/*
 * This routine handles the "ALARM" signal, which goes off after a set
 * amount of *REAL* time.  We check the CPU time usage -- if it has
 * exceeded the limit, we print a message and exit.  Otherwise, we
 * reschedule another alarm signal.
 */

	static
	RETSIGTYPE
check_limit (

int		sig
)
{
int32u		now;
int32u		delta_seconds;

	now = get_cpu_time ();

	if (now >= end_time) {
		(void) printf ("\n\n%% ***** CPU time limit of %lu"
				" seconds exceeded *****\n",
			       cpu_time_limit);
		exit (1);
	}

	delta_seconds = (end_time - now + TICKS_PER_SEC - 1) / TICKS_PER_SEC;

	signal (SIGALRM, check_limit);
	alarm (delta_seconds);
}
