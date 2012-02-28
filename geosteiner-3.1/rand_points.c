/***********************************************************************

	File:	rand_points.c
	Rev:	b-1
	Date:	02/28/2001

	Copyright (c) 1993, 2001 by David M. Warme

************************************************************************

	This program generates random point sets.

************************************************************************

	Modification Log:

	a-1:	02/21/93	warme
		: Created.
	a-2:	02/03/96	warme
		: Embedded the code that my own system uses for
		:  "rand/srand" so that this program produces
		:  the same (bad) mapping between seeds and
		:  sequences on all systems.
		: Added a much better random number generator
		:  using 64-bit shift register with XOR feedback.
		:  Perhaps now the point sets will start "looking"
		:  more random...
	b-1:	02/28/2001	warme
		: Process command line arguments in a more standard way.
		: Add usage() routine for when crud is given.

************************************************************************/

#include "steiner.h"
#include <time.h>


/*
 * Global Routines
 */

int			main (int, char **);


/*
 * External References
 */

	/* none */


/*
 * Local Equates
 */

#define	MAX_COORD	10000


/*
 * Local Routines
 */

static int16u		__rand (void);
static void		convert_seed (char *);
static int32u		cv_number (char *);
static void		decode_params (int, char **);
static void		generate_remainder_table (void);
static int32u		rand64 (int32u);
static int32u		random_coordinate (void);
static void		seed_to_string (char *, int32u, char *);
static void		usage (void);


/*
 * Local Variables
 */

static int		fflag;
static int32u		hi_seed;
static int32u		lo_seed;
static char *		me;
static int		npoints;
static int		oflag;
static bool		use_new_generator;

/*
 * This routine generates random point sets.
 */

	int
main (

int		argc,		/* IN - argument count. */
char **		argv		/* IN - argument vector. */
)
{
int			i;
int32u			x;
int32u			y;
char			sbuf [32];

	setbuf (stdout, NULL);

	decode_params (argc, argv);

	if (oflag NE 0) {
		if (use_new_generator) {
			/* Display seed in Hex. */
			seed_to_string (sbuf, 16, "0x");
		}
		else {
			/* Display seed in decimal. */
			seed_to_string (sbuf, 10, "");
		}
		fprintf ((oflag EQ 1) ? stdout : stderr,
			 "%% rand_points: seed = %s\n",
			 sbuf);
	}

	for (i = 0; i < npoints; i++) {
		x = random_coordinate ();
		if (NOT use_new_generator) {
			(void) __rand ();	/* Discard another number. */
		}
		y = random_coordinate ();
		(void) printf ("%5lu %5lu\n", x, y);
	}

	if (fflag NE 0) {
		if (use_new_generator) {
			/* Display seed in Hex. */
			seed_to_string (sbuf, 16, "0x");
		}
		else {
			/* Display seed in decimal. */
			seed_to_string (sbuf, 10, "");
		}
		fprintf ((oflag EQ 1) ? stdout : stderr,
			 "%% rand_points: final state = %s\n",
			 sbuf);
	}

	exit (0);
}

/*
 * This routine decodes the various command-line arguments.
 */

	static
	void
decode_params (

int		argc,
char **		argv
)
{
char *		ap;
char		c;

	--argc;
	me = *argv++;

	fflag = 0;
	oflag = 0;

	hi_seed = 1;
	lo_seed = 0;
	npoints = 10;

	while (argc > 0) {
		ap = *argv++;
		if (*ap NE '-') {
			npoints = cv_number (ap);
			break;
		}
		++ap;
		while ((c = *ap++) NE '\0') {
			switch (c) {
			case 'e':
				oflag = 2;
				break;

			case 'f':
				fflag = 1;
				break;

			case 'o':
				oflag = 1;
				break;

			case 'n':
				use_new_generator = TRUE;
				break;

			case 's':
				if (*ap EQ '\0') {
					if (argc <= 0) {
						usage ();
					}
					ap = *argv++;
					--argc;
				}
				convert_seed (ap);
				ap = "";
				break;

			case 'r':
				hi_seed = time (NULL);
				lo_seed = 0;
				break;

			default:
				usage ();
				break;
			}
		}
		--argc;
	}
}

/*
 * This routine prints out the proper usage and exits.
 */

static char *	arg_doc [] = {
	"",
	"\tGenerate N random points.  Default N is 10.",
	"",
	"\t-e\tWrite initial seed to stderr.",
	"\t-f\tWrite final seed to stderr (or stdout if -o given).",
	"\t-o\tWrite initial seed (and final seed if -f given) to stdout.",
	"\t-n\tUse newer (better) random generator.",
	"\t-s K\tSet random seed to K.",
	"\t-r\tRandomize: use current time as seed.",
	NULL
};

	static
	void
usage (void)

{
char **		pp;
char *		p;

	(void) fprintf (stderr,
			"\nUsage: %s [-efonr] [-s seed] [N]\n",
			me);

	pp = &arg_doc [0];
	while ((p = *pp++) NE NULL) {
		(void) fprintf (stderr, "%s\n", p);
	}
	exit (1);
}

/*
 * This routine converts the given string into a number.
 */

	static
	int32u
cv_number (

char *		num_string	/* IN - number to convert. */
)
{
char *		s;
int		c;
int32u		val;

	s = num_string;
	val = 0;
	while ((c = *s++) NE '\0') {
		if (NOT isdigit (c)) {
			(void) fprintf (stderr,
					"%s: `%s' is not a decimal number.\n",
					me, s);
			exit (1);
		}
		val = (val * 10) + (c - '0');
	}

	return (val);
}

/*
 * This routine converts the given ASCII string into a seed.  Either
 * decimal or Hex is accepted.
 */

	static
	void
convert_seed (

char *		string		/* IN - ASCII string to get seed from. */
)
{
char *		p;
int		i;
int		c;
int32u		carry;
int32u		base;
int16u		digs [4];

	p = string;

	base = 10;
	if ((p [0] EQ '0') AND (tolower (p [1]) EQ 'x')) {
		base = 16;
		p += 2;
	}

	for (i = 0; i < 4; i++) {
		digs [i] = 0;
	}

	for (;;) {
		c = *p++;
		if (c EQ '\0') {
			hi_seed = (((int32u) digs [2]) << 16) + digs [3];
			lo_seed = (((int32u) digs [0]) << 16) + digs [1];
			return;
		}

		/* Convert ASCII hex char to value */
		if (c < '0') c = -1;
		else if (c <= '9') c -= '0';
		else if (c < 'A') c = -1;
		else if (c <= 'F') c -= ('A' - 10);
		else if (c < 'a') c = -1;
		else if (c <= 'f') c -= ('a' - 10);
		else c = -1;

		if ((c < 0) OR (c >= base)) {
			fprintf (stderr,
				 "Warning: invalid seed: %s\n",
				 string);
			hi_seed = 1;
			lo_seed = 0;
			return;
		}

		carry = c;
		for (i = 3; i >= 0; i--) {
			carry += digs [i] * base;
			digs [i] = carry;
			carry >>= 16;
		}
	}
}

/*
 * This routine converts the current seed into a printable string
 * in the given base.
 */

	static
	void
seed_to_string (

char *		buf,		/* OUT - buffer to convert into. */
int32u		base,		/* IN - base to display in. */
char *		prefix		/* IN - prefix string for output. */
)
{
int		i;
int32u		rem;
char *		dstp;
int16u		digs [4];
char		tmp [32];

	/* Put prefix into destination buffer first... */
	while (*prefix NE '\0') {
		*buf++ = *prefix++;
	}

	if ((lo_seed EQ 0) AND (hi_seed EQ 0)) {
		*buf++ = '0';
		*buf++ = '\0';
		return;
	}

	digs [3] = (lo_seed >> 16);
	digs [2] = lo_seed;
	digs [1] = (hi_seed >> 16);
	digs [0] = hi_seed;

	dstp = &tmp [sizeof (tmp) / sizeof (tmp [0])];
	*--dstp = '\0';

	while ((digs [0] NE 0) OR (digs [1] NE 0) OR
	       (digs [2] NE 0) OR (digs [3] NE 0)) {
		rem = 0;
		for (i = 3; i >= 0; i--) {
			rem = (rem << 16) + digs [i];
			digs [i] = rem / base;
			rem %= base;
		}
		*--dstp = "0123456789ABCDEF" [rem];
	}
	(void) strcpy (buf, dstp);
}

/*
 * This routine computes a new random coordinate.  This is done in one
 * of two ways:  using the "bad" PDP-11 Unix method, or using a better
 * 64-bit shift register with XOR feedback.  (For backward compatibility,
 * the old way is the default...)
 */

	static
	int32u
random_coordinate (void)

{
int32u		n;

	if (use_new_generator) {
		return (rand64 (MAX_COORD));
	}

	/* Doing it the "bad" old way... */

	n = (((int32u) __rand ()) >> 2) ^ (((int32u) __rand ()) << 3);

	return (n % MAX_COORD);
}

/*
 * This routine produces random numbers the "bad" way (i.e. the way
 * the original PDP-11 Unix did).
 */

	static
	int16u
__rand (void)

{
long		x;

	x = 1103515245L * ((long) hi_seed) + 12345L;
	hi_seed = x;

	return ((x >> 16) & 0x7FFF);
}

/*
 * A better random-number generator package.  It uses a 64-bit shift
 * register to produce a pseudo-random bit stream with a period of
 * 2^64 - 1 bits.  Mathematically, the Kth shift register state is
 * exactly X^K mod G(X), where G(X) is a polynomial (of degree 64) that
 * is irreducible over the integers modulo 2.
 *
 * In this case, the Generator Polynomial is:
 *
 *	X^64 + X^62 + X^61 + X^56 + X^55 + X^54 + X^52 + X^50
 *	  + X^48 + X^44 + X^43 + X^42 + X^40 + X^39 + X^36 + X^35
 *	  + X^31 + X^28 + X^27 + X^26 + X^24 + X^23 + X^19 + X^17
 *	  + X^16 + X^14 + X^13 + X^10 + X^9 + X^7 + X^6 + X^5 + 1
 *
 * These polynomials are packed into 32-bit computer words as
 * follows:  the Most-Significant-Bit of the LO-reg is X^0, and
 * the Least-Significant bit of the HI-reg is X^63.
 */

#define	GEN_POLY_HI	0x19B8AB86
#define	GEN_POLY_LO	0x8766D1B9


/*
 * For run-time efficiency, the following table of remainders is pre-
 * initialized at compile time.  The values depend only upon the particular
 * choice of GEN_POLY.  Should a new polynomial be required, run the
 * "generate_remainder_table" routine to obtain new table initializations.
 *
 * NOTE: if you choose another generator polynomial, it MUST be
 *	 irreducible!!!
 */

static int32u		remainders [512] = {
	0x00000000,	0x00000000,	0x9ac99dbb,	0x11ea8045,
	0x3b5e9804,	0x10a45787,	0xa19705bf,	0x014ed7c2,
	0x76bd3008,	0x2148af0e,	0xec74adb3,	0x30a22f4b,
	0x4de3a80c,	0x31ecf889,	0xd72a35b7,	0x200678cc,
	0xed7a6010,	0x42915e1c,	0x77b3fdab,	0x537bde59,
	0xd624f814,	0x5235099b,	0x4ced65af,	0x43df89de,
	0x9bc75018,	0x63d9f112,	0x010ecda3,	0x72337157,
	0xa099c81c,	0x737da695,	0x3a5055a7,	0x629726d0,
	0xd4396352,	0xb653eb35,	0x4ef0fee9,	0xa7b96b70,
	0xef67fb56,	0xa6f7bcb2,	0x75ae66ed,	0xb71d3cf7,
	0xa284535a,	0x971b443b,	0x384dcee1,	0x86f1c47e,
	0x99dacb5e,	0x87bf13bc,	0x031356e5,	0x965593f9,
	0x39430342,	0xf4c2b529,	0xa38a9ef9,	0xe528356c,
	0x021d9b46,	0xe466e2ae,	0x98d406fd,	0xf58c62eb,
	0x4ffe334a,	0xd58a1a27,	0xd537aef1,	0xc4609a62,
	0x74a0ab4e,	0xc52e4da0,	0xee6936f5,	0xd4c4cde5,
	0xa6bf65d7,	0x5fd68167,	0x3c76f86c,	0x4e3c0122,
	0x9de1fdd3,	0x4f72d6e0,	0x07286068,	0x5e9856a5,
	0xd00255df,	0x7e9e2e69,	0x4acbc864,	0x6f74ae2c,
	0xeb5ccddb,	0x6e3a79ee,	0x71955060,	0x7fd0f9ab,
	0x4bc505c7,	0x1d47df7b,	0xd10c987c,	0x0cad5f3e,
	0x709b9dc3,	0x0de388fc,	0xea520078,	0x1c0908b9,
	0x3d7835cf,	0x3c0f7075,	0xa7b1a874,	0x2de5f030,
	0x0626adcb,	0x2cab27f2,	0x9cef3070,	0x3d41a7b7,
	0x72860685,	0xe9856a52,	0xe84f9b3e,	0xf86fea17,
	0x49d89e81,	0xf9213dd5,	0xd311033a,	0xe8cbbd90,
	0x043b368d,	0xc8cdc55c,	0x9ef2ab36,	0xd9274519,
	0x3f65ae89,	0xd86992db,	0xa5ac3332,	0xc983129e,
	0x9ffc6695,	0xab14344e,	0x0535fb2e,	0xbafeb40b,
	0xa4a2fe91,	0xbbb063c9,	0x3e6b632a,	0xaa5ae38c,
	0xe941569d,	0x8a5c9b40,	0x7388cb26,	0x9bb61b05,
	0xd21fce99,	0x9af8ccc7,	0x48d65322,	0x8b124c82,
	0x43b368dc,	0x8cdc55c3,	0xd97af567,	0x9d36d586,
	0x78edf0d8,	0x9c780244,	0xe2246d63,	0x8d928201,
	0x350e58d4,	0xad94facd,	0xafc7c56f,	0xbc7e7a88,
	0x0e50c0d0,	0xbd30ad4a,	0x94995d6b,	0xacda2d0f,
	0xaec908cc,	0xce4d0bdf,	0x34009577,	0xdfa78b9a,
	0x959790c8,	0xdee95c58,	0x0f5e0d73,	0xcf03dc1d,
	0xd87438c4,	0xef05a4d1,	0x42bda57f,	0xfeef2494,
	0xe32aa0c0,	0xffa1f356,	0x79e33d7b,	0xee4b7313,
	0x978a0b8e,	0x3a8fbef6,	0x0d439635,	0x2b653eb3,
	0xacd4938a,	0x2a2be971,	0x361d0e31,	0x3bc16934,
	0xe1373b86,	0x1bc711f8,	0x7bfea63d,	0x0a2d91bd,
	0xda69a382,	0x0b63467f,	0x40a03e39,	0x1a89c63a,
	0x7af06b9e,	0x781ee0ea,	0xe039f625,	0x69f460af,
	0x41aef39a,	0x68bab76d,	0xdb676e21,	0x79503728,
	0x0c4d5b96,	0x59564fe4,	0x9684c62d,	0x48bccfa1,
	0x3713c392,	0x49f21863,	0xadda5e29,	0x58189826,
	0xe50c0d0b,	0xd30ad4a4,	0x7fc590b0,	0xc2e054e1,
	0xde52950f,	0xc3ae8323,	0x449b08b4,	0xd2440366,
	0x93b13d03,	0xf2427baa,	0x0978a0b8,	0xe3a8fbef,
	0xa8efa507,	0xe2e62c2d,	0x322638bc,	0xf30cac68,
	0x08766d1b,	0x919b8ab8,	0x92bff0a0,	0x80710afd,
	0x3328f51f,	0x813fdd3f,	0xa9e168a4,	0x90d55d7a,
	0x7ecb5d13,	0xb0d325b6,	0xe402c0a8,	0xa139a5f3,
	0x4595c517,	0xa0777231,	0xdf5c58ac,	0xb19df274,
	0x31356e59,	0x65593f91,	0xabfcf3e2,	0x74b3bfd4,
	0x0a6bf65d,	0x75fd6816,	0x90a26be6,	0x6417e853,
	0x47885e51,	0x4411909f,	0xdd41c3ea,	0x55fb10da,
	0x7cd6c655,	0x54b5c718,	0xe61f5bee,	0x455f475d,
	0xdc4f0e49,	0x27c8618d,	0x468693f2,	0x3622e1c8,
	0xe711964d,	0x376c360a,	0x7dd80bf6,	0x2686b64f,
	0xaaf23e41,	0x0680ce83,	0x303ba3fa,	0x176a4ec6,
	0x91aca645,	0x16249904,	0x0b653bfe,	0x07ce1941,
	0x8766d1b9,	0x19b8ab86,	0x1daf4c02,	0x08522bc3,
	0xbc3849bd,	0x091cfc01,	0x26f1d406,	0x18f67c44,
	0xf1dbe1b1,	0x38f00488,	0x6b127c0a,	0x291a84cd,
	0xca8579b5,	0x2854530f,	0x504ce40e,	0x39bed34a,
	0x6a1cb1a9,	0x5b29f59a,	0xf0d52c12,	0x4ac375df,
	0x514229ad,	0x4b8da21d,	0xcb8bb416,	0x5a672258,
	0x1ca181a1,	0x7a615a94,	0x86681c1a,	0x6b8bdad1,
	0x27ff19a5,	0x6ac50d13,	0xbd36841e,	0x7b2f8d56,
	0x535fb2eb,	0xafeb40b3,	0xc9962f50,	0xbe01c0f6,
	0x68012aef,	0xbf4f1734,	0xf2c8b754,	0xaea59771,
	0x25e282e3,	0x8ea3efbd,	0xbf2b1f58,	0x9f496ff8,
	0x1ebc1ae7,	0x9e07b83a,	0x8475875c,	0x8fed387f,
	0xbe25d2fb,	0xed7a1eaf,	0x24ec4f40,	0xfc909eea,
	0x857b4aff,	0xfdde4928,	0x1fb2d744,	0xec34c96d,
	0xc898e2f3,	0xcc32b1a1,	0x52517f48,	0xddd831e4,
	0xf3c67af7,	0xdc96e626,	0x690fe74c,	0xcd7c6663,
	0x21d9b46e,	0x466e2ae1,	0xbb1029d5,	0x5784aaa4,
	0x1a872c6a,	0x56ca7d66,	0x804eb1d1,	0x4720fd23,
	0x57648466,	0x672685ef,	0xcdad19dd,	0x76cc05aa,
	0x6c3a1c62,	0x7782d268,	0xf6f381d9,	0x6668522d,
	0xcca3d47e,	0x04ff74fd,	0x566a49c5,	0x1515f4b8,
	0xf7fd4c7a,	0x145b237a,	0x6d34d1c1,	0x05b1a33f,
	0xba1ee476,	0x25b7dbf3,	0x20d779cd,	0x345d5bb6,
	0x81407c72,	0x35138c74,	0x1b89e1c9,	0x24f90c31,
	0xf5e0d73c,	0xf03dc1d4,	0x6f294a87,	0xe1d74191,
	0xcebe4f38,	0xe0999653,	0x5477d283,	0xf1731616,
	0x835de734,	0xd1756eda,	0x19947a8f,	0xc09fee9f,
	0xb8037f30,	0xc1d1395d,	0x22cae28b,	0xd03bb918,
	0x189ab72c,	0xb2ac9fc8,	0x82532a97,	0xa3461f8d,
	0x23c42f28,	0xa208c84f,	0xb90db293,	0xb3e2480a,
	0x6e278724,	0x93e430c6,	0xf4ee1a9f,	0x820eb083,
	0x55791f20,	0x83406741,	0xcfb0829b,	0x92aae704,
	0xc4d5b965,	0x9564fe45,	0x5e1c24de,	0x848e7e00,
	0xff8b2161,	0x85c0a9c2,	0x6542bcda,	0x942a2987,
	0xb268896d,	0xb42c514b,	0x28a114d6,	0xa5c6d10e,
	0x89361169,	0xa48806cc,	0x13ff8cd2,	0xb5628689,
	0x29afd975,	0xd7f5a059,	0xb36644ce,	0xc61f201c,
	0x12f14171,	0xc751f7de,	0x8838dcca,	0xd6bb779b,
	0x5f12e97d,	0xf6bd0f57,	0xc5db74c6,	0xe7578f12,
	0x644c7179,	0xe61958d0,	0xfe85ecc2,	0xf7f3d895,
	0x10ecda37,	0x23371570,	0x8a25478c,	0x32dd9535,
	0x2bb24233,	0x339342f7,	0xb17bdf88,	0x2279c2b2,
	0x6651ea3f,	0x027fba7e,	0xfc987784,	0x13953a3b,
	0x5d0f723b,	0x12dbedf9,	0xc7c6ef80,	0x03316dbc,
	0xfd96ba27,	0x61a64b6c,	0x675f279c,	0x704ccb29,
	0xc6c82223,	0x71021ceb,	0x5c01bf98,	0x60e89cae,
	0x8b2b8a2f,	0x40eee462,	0x11e21794,	0x51046427,
	0xb075122b,	0x504ab3e5,	0x2abc8f90,	0x41a033a0,
	0x626adcb2,	0xcab27f22,	0xf8a34109,	0xdb58ff67,
	0x593444b6,	0xda1628a5,	0xc3fdd90d,	0xcbfca8e0,
	0x14d7ecba,	0xebfad02c,	0x8e1e7101,	0xfa105069,
	0x2f8974be,	0xfb5e87ab,	0xb540e905,	0xeab407ee,
	0x8f10bca2,	0x8823213e,	0x15d92119,	0x99c9a17b,
	0xb44e24a6,	0x988776b9,	0x2e87b91d,	0x896df6fc,
	0xf9ad8caa,	0xa96b8e30,	0x63641111,	0xb8810e75,
	0xc2f314ae,	0xb9cfd9b7,	0x583a8915,	0xa82559f2,
	0xb653bfe0,	0x7ce19417,	0x2c9a225b,	0x6d0b1452,
	0x8d0d27e4,	0x6c45c390,	0x17c4ba5f,	0x7daf43d5,
	0xc0ee8fe8,	0x5da93b19,	0x5a271253,	0x4c43bb5c,
	0xfbb017ec,	0x4d0d6c9e,	0x61798a57,	0x5ce7ecdb,
	0x5b29dff0,	0x3e70ca0b,	0xc1e0424b,	0x2f9a4a4e,
	0x607747f4,	0x2ed49d8c,	0xfabeda4f,	0x3f3e1dc9,
	0x2d94eff8,	0x1f386505,	0xb75d7243,	0x0ed2e540,
	0x16ca77fc,	0x0f9c3282,	0x8c03ea47,	0x1e76b2c7,
};

/*
 * This function generates the hex numbers used to initialize the previous
 * table of polynomial remainders.  Simply select a polynomial, run this
 * function, and place its output into the initializer above, and recompile.
 */

#if 0

	static
	void
generate_remainder_table (void)

{
int		i;
int		j;
int32u		lo;
int32u		hi;

	for (i = 0; i < 256; i++) {
		lo = 0;
		hi = i;
		for (j = 0; j < 8; j++) {
			if ((hi & 0x00000001) != 0) {
				hi = (lo << 31) | (hi >> 1);
				lo >>= 1;
				hi ^= GEN_POLY_HI;
				lo ^= GEN_POLY_LO;
			}
			else {
				hi = (lo << 31) | (hi >> 1);
				lo >>= 1;
			}
		}
		printf ("\t0x%08lx,\t0x%08lx,", lo, hi);
		if ((i & 0x01) == 1) {
			printf ("\n");
		}
	}
}

#endif

/*
 * This routine returns a random 32-bit integer.  The value of zero is
 * only slightly less likely to occur than all others.
 */

	static
	int32u
rand64 (

int32u		modulus		/* IN - want random from 0..modulus-1 */
)
{
int32u		hi;
int32u		lo;
int32u *	p;
int		i;
int32u		rem;

#define	SHIFT8 \
	p = &remainders [2 * (hi & 0xFF)]; \
	hi = ((lo << 24) | (hi >> 8)) ^ p [1]; \
	lo = (lo >> 8) ^ p [0]

	hi = hi_seed;
	lo = lo_seed;

	SHIFT8;
	SHIFT8;
	SHIFT8;
	SHIFT8;

	hi_seed = hi;
	lo_seed = lo;

	if (modulus < 0x10000) {
		/* Fast code to mod a 64-bit number by a 16-bit number... */
		rem = lo % modulus;
		rem = (rem << 16) + (hi >> 16);
		rem %= modulus;
		rem = (rem << 16) + (hi & 0xFFFF);
		return (rem % modulus);
	}

	/* Must do things the hard way...  Sigh...  Actually, this */
	/* code only works for moduli up to 0x7FFFFFFF... */
	rem = lo % modulus;
	for (i = 0; i < 32; i++) {
		rem <<= 1;
		if ((hi & 0x80000000) NE 0) {
			++rem;
		}
		rem %= modulus;
		hi <<= 1;
	}

	return (rem);
}
