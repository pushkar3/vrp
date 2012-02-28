/***********************************************************************

	File:	genps.h
	Rev:	b-1
	Date:	02/28/2001

	Copyright (c) 1993, 2001 by David M. Warme

************************************************************************

	Stuff for generating plots in postscript.

************************************************************************

	Modification Log:

	a-1:	11/05/98	warme
		: Created.  Split off from steiner.h.
	b-1:	02/28/2001	warme
		: Pass scale_info explicitly to these routines.
		: Added plot_subtour.

************************************************************************/

#ifndef	GENPS_H
#define	GENPS_H

#include "steiner.h"


/*
 * Enumerated type to identify a plot size.  Big is a single plot per
 * page.  Small is 12 plots per page.
 */

enum plot_size {
	BIG_PLOT,
	SMALL_PLOT
};


/*
 * Function Prototypes.
 */

extern void		begin_plot (enum plot_size);
extern void		define_Plot_Terminals (struct pset *,
					       struct scale_info *);
extern void		draw_segment (struct point *,
				      struct point *,
				      struct scale_info *);
extern void		end_plot (char *);
extern void		overlay_plot_subset (char *,
					     bitmap_t *,
					     struct cinfo *,
					     enum plot_size);
extern void		plot_full_sets (bitmap_t *,
					struct cinfo *,
					enum plot_size);
extern void		plot_full_sets_grouped (bitmap_t *,
						struct cinfo *,
						enum plot_size);
extern void		plot_lp_solution (char *,
					  double *,
					  struct cinfo *,
					  enum plot_size);
extern void		plot_subtour (char *,
				      double *,
				      struct cinfo *,
				      bitmap_t *,
				      enum plot_size);

#endif
