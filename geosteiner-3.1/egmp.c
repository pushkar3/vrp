/***********************************************************************

	File:	egmp.c
	Rev:	b-1
	Date:	11/25/2000

	Copyright (c) 2000, 2001 by David M. Warme

************************************************************************

	Support routines for the EFST generator that use the
	GNU Multi-Precision arithmetic library (GMP -- if we have
	it) to compute certain items with high accuracy.

************************************************************************

	Modification Log:

	b-1:	11/25/2000	warme
		: Created.

************************************************************************/


#include "config.h"

#ifdef HAVE_GMP

#include "egmp.h"
#include "efst.h"
#include "gmp.h"
#include "steiner.h"


/*
 * Global Routines
 */

double		compute_EFST_length (struct einfo * eip, struct eqp_t * eqpt);
void		qr3_clear (qr3_t * p);
void		qr3_init (qr3_t * p);
void		update_eqpoint_and_displacement (struct einfo *	eip,
						 struct eqp_t *	eqpk);

extern int	Multiple_Precision;

/*
 * Local Equates
 */

#define	DEBUG_PRINT	0


/*
 * Local Routines
 */

static void		compute_eqpoint (struct qr3_point *	p,
					 struct einfo *		eip,
					 struct eqp_t *		eqpk);
static void		qr3_mul (qr3_t * dst, qr3_t * p1, qr3_t * p2);
static double		qr3_to_double (qr3_t *);
static void		r_to_q (mpq_t q_dst, double r_src);

/*
 * A routine to recompute the coordinates of the given eq-point and its
 * corresponding displacement vector.  Using purely floating point
 * arithmetic, errors tend to accumulate as the eq-point order grows.
 * Once we are pretty sure that we are going to save this eq-point, we
 * call this routine that re-computes these quantities correct to within
 * 1/2 ULP of the floating point arithmetic -- thus preventing the
 * accumulation of such errors.
 */

	void
update_eqpoint_and_displacement (

struct einfo *		eip,		/* IN - EFST generation info */
struct eqp_t *		eqpk		/* IN - eq-point to recompute */
)
{
double			nx;
double			ny;
double			ex;
double			ey;
mpq_t			rtmp;
struct point *		tp;
struct qr3_point *	ep;
qr3_t			dv_tmp;

	mpq_init (rtmp);
	qr3_init (&dv_tmp);

	ep = &(eip -> cur_eqp);

	compute_eqpoint (ep, eip, eqpk);

#if DEBUG_PRINT
	printf ("\nPoint %3d: original = (%24.20f, %24.20f)\n",
		eqpk -> E.pnum,
		eqpk -> E.x,
		eqpk -> E.y);
#endif

	nx = qr3_to_double (&(ep -> x));
	ny = qr3_to_double (&(ep -> y));

#if DEBUG_PRINT
	printf ("\tNew:          (%24.20f, %24.20f)\n", nx, ny);

	ex = fabs (nx - eqpk -> E.x) / nx;
	ey = fabs (ny - eqpk -> E.y) / ny;

	printf ("\tErrors:\t%d\t(%14g, %14g)\n", eqpk -> S, ex, ey);

	printf ("\tSymbolic: X = ");
	mpz_out_str (stdout, 10, mpq_numref (ep -> x.a));
	printf (" / ");
	mpz_out_str (stdout, 10, mpq_denref (ep -> x.a));
	printf (" + ");
	mpz_out_str (stdout, 10, mpq_numref (ep -> x.b));
	printf (" * sqrt (3) / ");
	mpz_out_str (stdout, 10, mpq_denref (ep -> x.b));
	printf ("\n");

	printf ("\tSymbolic: Y = ");
	mpz_out_str (stdout, 10, mpq_numref (ep -> y.a));
	printf (" / ");
	mpz_out_str (stdout, 10, mpq_denref (ep -> y.a));
	printf (" + ");
	mpz_out_str (stdout, 10, mpq_numref (ep -> y.b));
	printf (" * sqrt (3) / ");
	mpz_out_str (stdout, 10, mpq_denref (ep -> y.b));
	printf ("\n");
#endif

	eqpk -> E.x = nx;
	eqpk -> E.y = ny;

	tp = &(eip -> eqp [eqpk -> DV.pnum].E);

	r_to_q (rtmp, tp -> x);
	mpq_sub (dv_tmp.a, ep -> x.a, rtmp);
	mpq_set (dv_tmp.b, ep -> x.b);
	eqpk -> DV.x = qr3_to_double (&dv_tmp);

	r_to_q (rtmp, tp -> y);
	mpq_sub (dv_tmp.a, ep -> y.a, rtmp);
	mpq_set (dv_tmp.b, ep -> y.b);
	eqpk -> DV.y = qr3_to_double (&dv_tmp);

	qr3_clear (&dv_tmp);
	mpq_clear (rtmp);
}

/*
 * Compute the exact coordinates of the given eq-point as a member
 * of the field Q(sqrt(3)) -- i.e., as (A + B*sqrt(3), C + D*sqrt(3))
 * where A, B, C and D are exact rational numbers.
 */

	void
compute_eqpoint (

struct qr3_point *	out,		/* OUT - eq-point coordinates */
struct einfo *		eip,		/* IN - EFST generation info */
struct eqp_t *		eqpk		/* IN - eq-point to calculate */
)
{
int			c;
int			n;
struct qr3_point	P, Q, R;
mpq_t			c2, c3, t1, t2;

	if (eqpk -> L EQ NULL) {
		/* Base case -- a terminal. */

		r_to_q (out -> x.a, eqpk -> E.x);
		r_to_q (out -> y.a, eqpk -> E.y);

		mpq_set_ui (out -> x.b, 0, 1);
		mpq_set_ui (out -> y.b, 0, 1);

		return;
	}

	/* Recurse, let P be the right point, and Q the left eq-point. */

	qr3_init (&P.x);
	qr3_init (&P.y);
	qr3_init (&Q.x);
	qr3_init (&Q.y);
	qr3_init (&R.x);
	qr3_init (&R.y);

	compute_eqpoint (&P, eip, eqpk -> R);
	compute_eqpoint (&Q, eip, eqpk -> L);

	mpq_sub (R.x.a, Q.x.a, P.x.a);
	mpq_sub (R.x.b, Q.x.b, P.x.b);
	mpq_sub (R.y.a, Q.y.a, P.y.a);
	mpq_sub (R.y.b, Q.y.b, P.y.b);

	mpz_init_set_ui (mpq_numref (c2), 2);
	mpz_init_set_ui (mpq_denref (c2), 1);

	mpz_init_set_ui (mpq_numref (c3), 3);
	mpz_init_set_ui (mpq_denref (c3), 1);

	mpq_init (t1);
	mpq_init (t2);

	mpq_mul (t1, c3, R.y.b);
	mpq_sub (t2, R.x.a, t1);
	mpq_div (t1, t2, c2);
	mpq_add (out -> x.a, P.x.a, t1);

	mpq_sub (t1, R.x.b, R.y.a);
	mpq_div (t2, t1, c2);
	mpq_add (out -> x.b, P.x.b, t2);

	mpq_mul (t1, c3, R.x.b);
	mpq_add (t2, R.y.a, t1);
	mpq_div (t1, t2, c2);
	mpq_add (out -> y.a, P.y.a, t1);

	mpq_add (t1, R.y.b, R.x.a);
	mpq_div (t2, t1, c2);
	mpq_add (out -> y.b, P.y.b, t2);

	mpq_clear (t2);
	mpq_clear (t1);
	mpq_clear (c3);
	mpq_clear (c2);

	qr3_clear (&R.y);
	qr3_clear (&R.x);
	qr3_clear (&Q.y);
	qr3_clear (&Q.x);
	qr3_clear (&P.y);
	qr3_clear (&P.x);
}

/*
 * Convert a double to a rational number.
 */

	static
	void
r_to_q (

mpq_t		q_dst,		/* OUT - rational number */
double		r_src		/* IN  - double to convert */
)
{
int		mant_size;
int		expon;
mpf_t		x;
mpz_t		tmp;
mpz_ptr		np;
mpz_ptr		dp;

	/* First, convert to MP-floating point. */
	mpf_init2 (x, 64);
	mpf_set_d (x, r_src);

	/* Access mantissa as an mpz_t... */
	tmp [0]._mp_alloc	= x [0]._mp_prec + 1;
	tmp [0]._mp_size	= x [0]._mp_size;
	tmp [0]._mp_d		= x [0]._mp_d;

	np = mpq_numref (q_dst);
	dp = mpq_denref (q_dst);

	/* Since mantissa is fractional (not an integer), we	*/
	/* must include the mantissa size in the exponent...	*/
	mant_size = x [0]._mp_size;
	if (mant_size < 0) {
		mant_size = -mant_size;
	}
	expon = x [0]._mp_exp - mant_size;
	if (expon < 0) {
		/* Set numerator */
		mpz_set (np, tmp);

		/* Divide by limb-base**(- expon) */
		mpz_set_ui (dp, 1);
		mpz_mul_2exp (dp, dp, 32 * (- expon));
	}
	else {
		/* Set denominator */
		mpz_set_ui (dp, 1);

		/* Set numerator to proper shift factor. */
		mpz_set_ui (np, 1);
		mpz_mul_2exp (np, np, 32 * expon);
		/* multiply mantissa in... */
		mpz_mul (np, np, tmp);
	}

	mpq_canonicalize (q_dst);

	mpf_clear (x);
	/* DO NOT clear tmp! */
}

/*
 * Initialize an Q(sqrt(3)) number.
 */

	void
qr3_init (

qr3_t *		p
)
{
	mpq_init (p -> a);
	mpq_init (p -> b);
}


/*
 * Clear a Q(sqrt(3)) number.
 */

	void
qr3_clear (

qr3_t *		p
)
{
	mpq_clear (p -> a);
	mpq_clear (p -> b);
}

/*
 * Accurately translate an element Z of Q(sqrt(3)) into a double.
 * Let Z = A + B * sqrt(3), where A and B are rational.
 *
 *	Z - A = B * sqrt(3)	==>	Z^2 - 2*A*Z + A^2 - 3*B^2 = 0
 *
 * We solve start with a floating point approximation of Z, convert
 * it to a rational and then apply Newton iterations until we
 * achieve the desired precision.  The termination test is as
 * follows:
 *
 *	|Z - (A + B*sqrt(3))| < eps * (A + B*sqrt(3)),
 *
 *	Z^2 - 2*(A + B*sqrt(3))*Z + A^2 + 2*A*B*sqrt(3) + 3*B^2
 *		< eps^2 * (A^2 + 2*A*B*sqrt(3) + 3*B^2),
 *
 *	Z^2 - 2*A*Z + A^2 + 3*B^2 - eps^2*(A^2 + 3*B^2)
 *			< (2*B*Z - (1 - eps^2) 2*A*B) * sqrt(3),
 *
 *	Z^2 - 2*A*Z + (1 - eps^2)*(A^2 + 3*B^2)
 *			< (Z - (1 - eps^2)*A) * 2 * B * sqrt(3),
 *
 * The obvious thing here is to square both sides and compare.
 * However, one must be careful regarding the signs of the two
 * sides when doing this.
 */

	static
	double
qr3_to_double (

qr3_t *		p		/* IN - A + B*sqrt(3) to convert */
)
{
int		flag;
mpq_t		z;
mpq_t		t1;
mpq_t		t2;
mpq_t		t3;
mpq_t		t4;
mpq_t		t5;
mpq_t		t6;
mpq_t		c_1_minus_eps2;
double		xf;

#define EPS	64		/* relative error of 2^(-EPS) */

	if (mpz_sgn (mpq_numref (p -> b)) EQ 0) {
		return (mpq_get_d (p -> a));
	}

	mpq_init (z);
	mpq_init (t1);
	mpq_init (t2);
	mpq_init (t3);
	mpq_init (t4);
	mpq_init (t5);
	mpq_init (t6);
	mpq_init (c_1_minus_eps2);

	/* c_1_minus_eps2 is now 0 / 1.  Set it to be	*/
	/* 1 - 2**( - 2 * EPS).				*/

	mpz_mul_2exp (mpq_denref (c_1_minus_eps2),
		      mpq_denref (c_1_minus_eps2),
		      2*EPS);
	mpz_sub_ui (mpq_numref (c_1_minus_eps2),
		    mpq_denref (c_1_minus_eps2),
		    1);

	/* Compute t1 = (A^2 - 3*B^2) */

	mpq_mul (t3, p -> a, p -> a);
	mpq_mul (t4, p -> b, p -> b);
	mpz_mul_ui (mpq_numref (t4), mpq_numref (t4), 3);
	mpq_canonicalize (t4);
	mpq_sub (t1, t3, t4);

	/* Compute t2 = (1 - eps^2) * (A^2 + 3*B^2).	*/

	mpq_add (t2, t3, t4);
	mpq_mul (t2, t2, c_1_minus_eps2);

	/* Compute t3 = (1 - eps^2)*A.	*/

	mpq_mul (t3, c_1_minus_eps2, p -> a);

	/* Compute t4 = 12 * B^2 */

	mpz_mul_ui (mpq_numref (t4), mpq_numref (t4), 3);
	mpq_canonicalize (t4);

	/* Constants for the termination test are now ready.	*/
	/* Time to start up the Newton iteration.		*/

	/* Approximate solution using floating point as a start... */

	xf = mpq_get_d (p -> a) + mpq_get_d (p -> b) * sqrt(3.0);

	/* Put into rational form. */
	r_to_q (z, xf);

	for (;;) {
		/* Compute Z = new Z, via Newton */
		mpq_mul (t5, z, z);
		mpq_sub (t5, t5, t1);
		mpq_sub (t6, z, p -> a);
		mpz_mul_ui (mpq_numref (t6), mpq_numref (t6), 2);
		mpq_canonicalize (t6);
		mpq_div (z, t5, t6);

#if DEBUG_PRINT
		printf ("	New z is %24.20f\n", mpq_get_d (z));
#endif

		if (Multiple_Precision <= 1) {
			/* Termination test takes time.  At level 1,	*/
			/* we just do 1 Newton iteration and quit.	*/
			break;
		}

		/* Now test termination... */
		mpq_set (t6, p -> a);
		mpz_mul_ui (mpq_numref (t6), mpq_numref (t6), 2);
		mpq_canonicalize (t6);
		mpq_sub (t5, z, t6);
		mpq_mul (t5, t5, z);
		mpq_add (t5, t5, t2);

		mpq_sub (t6, z, t3);
		mpq_mul (t6, t6, p -> b);

		/* Time to check the signs... */

		flag = 1;
		if (mpz_sgn (mpq_numref (t5)) < 0) {
			if (mpz_sgn (mpq_numref (t6)) >= 0) break;
			flag = -1;
		}
		else if (mpz_sgn (mpq_numref (t6)) < 0) {
			continue;
		}

		mpq_mul (t5, t5, t5);
		mpq_mul (t6, t6, t6);
		mpz_mul_ui (mpq_numref (t6), mpq_numref (t6), 12);
		mpq_canonicalize (t6);

#if DEBUG_PRINT
		printf ("	t5 = %24g, t6 = %24g\n",
			mpq_get_d (t5),
			mpq_get_d (t6));
#endif

		if (mpq_cmp (t5, t6) * flag < 0) break;
	}

	xf = mpq_get_d (z);

	mpq_clear (c_1_minus_eps2);
	mpq_clear (t6);
	mpq_clear (t5);
	mpq_clear (t4);
	mpq_clear (t3);
	mpq_clear (t2);
	mpq_clear (t1);
	mpq_clear (z);

	return (xf);
}

/*
 * Routine to compute the length of a given EFST (i.e., Simpson line)
 * to within 1/2 ULP.  (Modulo good behavior of the sqrt() function...)
 *
 * The exact position of the eq-point end of the Simpson line has already
 * been stored in eip -> cur_eqp.
 */

	double
compute_EFST_length (

struct einfo *		eip,		/* IN - EFST generation info */
struct eqp_t *		eqpt		/* IN - terminal end of Simpson line */
)
{
qr3_t			x;
qr3_t			y;
double			len;

	qr3_init (&x);
	qr3_init (&y);

	r_to_q (x.a, - eqpt -> E.x);
	r_to_q (y.a, - eqpt -> E.y);

	mpq_add (x.a, x.a, eip -> cur_eqp.x.a);
	mpq_add (y.a, y.a, eip -> cur_eqp.y.a);
	mpq_set (x.b, eip -> cur_eqp.x.b);
	mpq_set (y.b, eip -> cur_eqp.y.b);

	qr3_mul (&x, &x, &x);
	qr3_mul (&y, &y, &y);

	mpq_add (x.a, x.a, y.a);
	mpq_add (x.b, x.b, y.b);

	len = sqrt (qr3_to_double (&x));

	qr3_clear (&y);
	qr3_clear (&x);

	return (len);
}

/*
 * Multiply two elements of Q(sqrt(3)).
 */

	static
	void
qr3_mul (

qr3_t *		dst,		/* OUT - destination */
qr3_t *		p1,		/* IN - first operand */
qr3_t *		p2		/* IN - second operand */
)
{
mpq_t		res_a;
mpq_t		tmp1;
mpq_t		tmp2;

	mpq_init (res_a);
	mpq_init (tmp1);
	mpq_init (tmp2);

	mpq_mul (tmp1, p1 -> a, p2 -> a);
	mpq_mul (tmp2, p1 -> b, p2 -> b);
	mpz_mul_ui (mpq_numref (tmp2), mpq_numref (tmp2), 3);
	mpq_canonicalize (tmp2);
	mpq_add (res_a, tmp1, tmp2);

	mpq_mul (tmp1, p1 -> a, p2 -> b);
	mpq_mul (tmp2, p1 -> b, p2 -> a);
	mpq_add (dst -> b, tmp1, tmp2);

	mpq_set (dst -> a, res_a);

	mpq_clear (tmp2);
	mpq_clear (tmp1);
	mpq_clear (res_a);
}

#endif
