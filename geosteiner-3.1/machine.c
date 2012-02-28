/***********************************************************************

	File:	machine.c
	Rev:	a-1
	Date:	01/12/99

	Copyright (c) 1999, 2001 by David M. Warme

************************************************************************

	This file contains routines for obtaining a string that
	describes the machine running the program.

************************************************************************

	Modification Log:

	a-1:	01/12/99	warme
		: Created.

************************************************************************/

#include "config.h"
#include "steiner.h"


/* Determine how to implement get_machine_string... */
#if defined(MACHDESC)
	/* This method overrides all other options... */
	#define PRINT_CONSTANT_STRING
#elif defined(UNAME_FUNCTION_WORKS)
	#define USE_UNAME_SYSCALL
#elif defined(HAVE_POPEN) AND defined(HAVE_PCLOSE) AND defined(UNAME_PATH)
	#define USE_UNAME_COMMAND
#else
	#define MACHDESC		"Unknown"
	#define PRINT_CONSTANT_STRING
#endif

#ifdef USE_UNAME_SYSCALL
#include <sys/utsname.h>
#endif


/*
 * Global Routines
 */

char *		get_machine_string (void);


/*
 * Local Variables
 */

static char *	machine_string;

/*
 * This routine obtains a machine description string from the "uname" 
 * system call.
 */

#ifdef USE_UNAME_SYSCALL

	char *
get_machine_string (void)

{
size_t		len;
char *		s;
struct utsname	un;

	if (machine_string NE NULL) return (machine_string);

	if (uname (&un) < 0) {
		machine_string = "Unknown";
		return (machine_string);
	}

	len	= strlen (un.sysname)
		+ strlen (un.nodename)
		+ strlen (un.release)
		+ strlen (un.version)
		+ strlen (un.machine);

	s = NEWA (len + 5, char);

	sprintf (s, "%s %s %s %s %s",
		 un.sysname, un.nodename, un.release, un.version, un.machine);

	machine_string = s;

	return (s);
}

#endif

/*
 * This routine obtains a machine description string by invoking the
 * "uname -a" command.
 */

#ifdef USE_UNAME_COMMAND

	char *
get_machine_string (void)

{
int		c;
char *		p;
char *		endp;
size_t		size;
FILE *		fp;
char		command [1024];
char		buffer [80];

	if (machine_string NE NULL) return (machine_string);

	sprintf (command, "%s -a", UNAME_PATH);
	fp = popen (command, "r");
	if (fp EQ NULL) goto dont_know;

	p = buffer;
	endp = p + sizeof (buffer) - 1;
	for (;;) {
		c = getc (fp);
		if (c < 0) break;
		if (p >= endp) break;
		if (c EQ '\n') break;
		*p++ = c;
	}
	*p++ = '\0';
	size = p - buffer;
	p = NEWA (size, char);
	strcpy (p, buffer);
	machine_string = p;

	/* Drain the stream... */
	while (c >= 0) {
		c = getc (fp);
	}

	if (pclose (fp) NE 0) {
dont_know:
		machine_string = "Unknown";
	}
	return (machine_string);
}

#endif

/*
 * This routine returns a fixed machine description string.
 */

#ifdef PRINT_CONSTANT_STRING

	char *
get_machine_string (void)

{
	return (MACHDESC);
}

#endif
