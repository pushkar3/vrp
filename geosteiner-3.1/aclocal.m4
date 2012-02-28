dnl *******************************************************************
dnl
dnl	File:	aclocal.m4
dnl	Rev:	a-2
dnl	Date:	02/28/2001
dnl
dnl	Copyright (c) 1998, 2001 by David M. Warme
dnl
dnl *******************************************************************
dnl
dnl	Customized versions of some of the GNU autoconfig macros.
dnl
dnl *******************************************************************
dnl
dnl	Modification Log:
dnl
dnl	a-1:	12/20/98	warme
dnl		: Created.
dnl		: Removed -g from the default CFLAGS.
dnl	a-2:	02/28/2001	warme
dnl		: Changes for 3.1 release, including tests for
dnl		:  Floating point precision fixes and checks
dnl		:  for -lpthread.
dnl
dnl *******************************************************************

undefine([AC_PROG_CC])

AC_DEFUN(AC_PROG_CC,
[AC_BEFORE([$0], [AC_PROG_CPP])dnl
AC_CHECK_PROG(CC, gcc, gcc)
if test -z "$CC"; then
  AC_CHECK_PROG(CC, cc, cc, , , /usr/ucb/cc)
  test -z "$CC" && AC_MSG_ERROR([no acceptable cc found in \$PATH])
fi

AC_PROG_CC_WORKS
AC_PROG_CC_GNU

if test $ac_cv_prog_gcc = yes; then
  GCC=yes
dnl Check whether -g works, even if CFLAGS is set, in case the package
dnl plays around with CFLAGS (such as to build both debugging and
dnl normal versions of a library), tasteless as that idea is.
  ac_test_CFLAGS="${CFLAGS+set}"
  ac_save_CFLAGS="$CFLAGS"
  CFLAGS=
  AC_PROG_CC_G
  if test "$ac_test_CFLAGS" = set; then
    CFLAGS="$ac_save_CFLAGS"
# elif test $ac_cv_prog_cc_g = yes; then
#   CFLAGS="-g -O2"
  else
    CFLAGS="-O2"
  fi
else
  GCC=
  test "${CFLAGS+set}" = set || CFLAGS="-O"
fi
])

AC_DEFUN(AC_CHECK_GCC_FLOAT_STORE,
[dnl Check for best solution for the Intel x86 FPU on Linux case...
dnl (no need to verify GCC here -- it either fixes the problem or it doesn't.)
AC_CACHE_CHECK(for Intel x86 floating point precision fix,
		ac_cv_x86_fp_precision_fix,
		AC_LANG_SAVE
		AC_LANG_C
		AC_TRY_RUN(dnl
[#if STDC_HEADERS
#include <stdlib.h>
#endif
AC_CHECK_X86_PRECISION_FIX_CODE],
ac_cv_x86_fp_precision_fix=yes,
ac_cv_x86_fp_precision_fix=no,
ac_cv_x86_fp_precision_fix=no)
AC_LANG_RESTORE)
if test "$ac_cv_x86_fp_precision_fix" = "yes"
then
	AC_DEFINE(HAVE_X86_FLOATING_POINT_PRECISION_FIX)
	# This is a better solution.  No need to use -ffloat-store
	ac_cv_gcc_needs_float_store=no
else
 if test "$GCC" != "yes"
 then
	# Not GCC -- can't use -ffloat-store
	ac_cv_gcc_needs_float_store=no
 else
	AC_CACHE_CHECK(if gcc -ffloat-store needed,
		       ac_cv_gcc_needs_float_store,
	AC_LANG_SAVE
	AC_LANG_C
	AC_TRY_RUN(dnl
[#if STDC_HEADERS
#include <stdlib.h>
#endif
AC_CHECK_FLOAT_STORE_CODE],
ac_cv_gcc_needs_float_store=yes,
ac_cv_gcc_needs_float_store=no,
ac_cv_gcc_needs_float_store=no)
AC_LANG_RESTORE)
 fi
fi
GCC_FLOAT_STORE=
if test "$ac_cv_gcc_needs_float_store" = "yes"
then
	GCC_FLOAT_STORE='-ffloat-store'
fi
AC_SUBST(GCC_FLOAT_STORE)
])

AC_DEFUN(AC_CHECK_X86_PRECISION_FIX_CODE,
[
#include <fpu_control.h>
extern void store_double (
#ifdef __STDC__
double *, double
#endif
);

main ()

{
fpu_control_t	cw, orig_cw;
double		x1, eps1;
double		x2, eps2;
double		x3, eps3;
double		z;

	/* Compute an epsilon with results forced through memory. */

	eps1 = 1.0;
	x1 = 1.0;
	for (;;) {
		store_double (&z, x1 + 1.0);
		if (z == 1.0) break;
		eps1 = x1;
		x1 *= 0.5;
	}

/* These guys should learn how to write headers... */
#define	_FPU_CW_PC	_FPU_EXTENDED

	_FPU_GETCW (orig_cw);

	/* Compute an epsilon using extended precision enabled. */

	cw = (orig_cw & ~_FPU_CW_PC) | _FPU_EXTENDED;
	_FPU_SETCW (cw);

	eps2 = 1.0;
	x2 = 1.0;
	for (;;) {
		z = x2 + 1.0;
		if (z == 1.0) break;
		eps2 = x2;
		x2 *= 0.5;
	}

	/* Compute an epsilon with double precision forced. */

	cw = (orig_cw & ~_FPU_CW_PC) | _FPU_DOUBLE;
	_FPU_SETCW (cw);

	eps3 = 1.0;
	x3 = 1.0;
	for (;;) {
		z = x3 + 1.0;
		if (z == 1.0) break;
		eps3 = x3;
		x3 *= 0.5;
	}

	_FPU_SETCW (orig_cw);

	if ((eps1 != eps2) && (eps1 == eps3)) {
		/* Houston, we have a problem... */
		/* Roger, set PC=double to fix it... */
		exit (0);
	}
	/* Any failure here says we can't use this fix -- it doesn't work. */
	exit (1);
}

	void
store_double (

double *	p,
double		val
)
{
	*p = val;
}
])

AC_DEFUN(AC_CHECK_FLOAT_STORE_CODE,
[
extern void store_double (
#ifdef __STDC__
double *, double
#endif
);

main ()

{
double		x, eps1;
double		y, eps2;
double		z;

	eps1 = 1.0;
	x = 1.0;
	for (;;) {
		z = x + 1.0;
		if (z == 1.0) break;
		eps1 = x;
		x *= 0.5;
	}

	eps2 = 1.0;
	y = 1.0;
	for (;;) {
		store_double (&z, y + 1.0);
		if (z == 1.0) break;
		eps2 = y;
		y *= 0.5;
	}

	if (eps1 == eps2) {
		/* Exiting with failure here says we don't */
		/* need -ffloat-store. */
		exit (1);
	}

	/* Houston, we have a problem...  Roger, use -ffloat-store. */
	exit (0);
}

	void
store_double (

double *	p,
double		val
)
{
	*p = val;
}
])

dnl
dnl AC_CPLEX_CHECK_LIBPTHREAD - Check to see if we need to
dnl use -lpthread when linking CPLEX.
dnl
AC_DEFUN(AC_CPLEX_CHECK_LIBPTHREAD,
[dnl
 AC_CACHE_CHECK(if CPLEX needs -lpthread, ac_cv_cplex_libpthread,
	[ dnl
	cpx_save_LIBS=${LIBS}
	cpx_save_CPPFLAGS=${CPPFLAGS}
	LIBS="${LIB} ${ac_cv_cplex_lib}"
	[cpxhdrdir="`echo $ac_cv_cplex_header | sed -e 's!/[^/]*$!!'`"]
	CPPFLAGS="${CPPFLAGS} -I${cpxhdrdir}"

	AC_TRY_LINK(AC_CPLEX_CHECK_PTHREAD_INCLUDES, dnl
		    AC_CPLEX_CHECK_PTHREAD_CODE, dnl
		    ac_cv_cplex_libpthread=no,	dnl -lpthread not needed
		    # Link did not work.  Try adding -lpthread.
		    LIBS="${LIBS} -lpthread"
		    AC_TRY_LINK(AC_CPLEX_CHECK_PTHREAD_INCLUDES, dnl
				AC_CPLEX_CHECK_PTHREAD_CODE, dnl
				ac_cv_cplex_libpthread=yes,
				ac_cv_cplex_libpthread=no))

	LIBS=${cpx_save_LIBS}
	CPPFLAGS=${cpx_save_CPPFLAGS}
 ])
])

dnl
dnl Macro to generate proper include files for -lpthread checking code.
dnl
AC_DEFUN(AC_CPLEX_CHECK_PTHREAD_INCLUDES, dnl
[#if CPLEX < 40
 #include "cpxdefs.inc"
#else
 #include "cplex.h"
#endif
])

dnl
dnl Code to test if CPLEX needs -lpthread
dnl
AC_DEFUN(AC_CPLEX_CHECK_PTHREAD_CODE,
[
#if CPLEX >= 50
 int		status;
/* CPXENVptr	cpx_env; */
 CPXLPptr	lp;

/*	cplex_env = CPXopenCPLEX (&status); */
	lp = NULL;
	CPXprimopt ((void *) 0, lp);
	CPXdualopt ((void *) 0, lp);
#else
 #if CPLEX >= 40
  int		status;
  CPXENVptr	cpx_env;
  CPXLPptr	lp;

	cplex_env = CPXopenCPLEX (&status);
	lp = NULL;
	CPXoptimize (cpx_env, lp);
	CPXdualopt (cpx_env, lp);
 #else
  CPXLPptr	lp;

	lp = NULL;
	optimize (lp);
	dualopt (lp);
 #endif
#endif
])

dnl AC_DEFINE_EXPAND_VALUE(VARIABLE [, VALUE])
AC_DEFUN(AC_DEFINE_EXPAND_VALUE,
[dnl call AC_DEFINE after expanding given value
ac_expand_text1=""
ac_expand_text2="`eval echo \"$2\"`"
while test "$ac_expand_text1" != "$ac_expand_text2"
do
	ac_expand_text3="`eval echo \"$ac_expand_text2\"`"
	ac_expand_text1="$ac_expand_text2"
	ac_expand_text2="$ac_expand_text3"
done
AC_DEFINE_UNQUOTED([$1], "$ac_expand_text1")
])

dnl AC_FIND_FILE(VARIABLE, FILENAME, TESTFLAG, DIRS)
AC_DEFUN(AC_FIND_FILE,
[dnl loop over the dirs, until found
  $1=''
  for cv_ff_dir in $4
  do
	if test "$3" "$cv_ff_dir/$2"
	then
		$1="$cv_ff_dir/$2"
		break
	fi
  done
])
