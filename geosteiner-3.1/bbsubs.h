/***********************************************************************

	File:	bbsubs.h
	Rev:	a-1
	Date:	02/28/2001

	Copyright (c) 1998, 2001 by David M. Warme

************************************************************************

	Declarations of low-level branch-and-cut routines.

************************************************************************

	Modification Log:

	a-1:	02/28/2001	warme
		: Split off from bb.h.  Changes for 3.1 release.

************************************************************************/

#ifndef	BBSUBS_H
#define	BBSUBS_H


#include "steiner.h"

struct bbnode;
struct bbtree;
struct bbinfo;

extern void		add_bbnode (struct bbinfo *, int, int, double);
extern void		append_node_to_tree (struct bbnode *, struct bbtree *);
extern void		bbheap_insert (struct bbnode *, struct bbtree *, int);
extern struct bbtree *	create_bbtree (int);
extern void		delete_node_from_bbtree (struct bbnode *,
						 struct bbtree *);
extern void		destroy_bbinfo (struct bbinfo *);
extern void		shutdown_lp_solver (void);
extern void		startup_lp_solver (void);

#endif
