/*  t_dens                                                                  */
/*                                                                          */
/*  kdens, kdcomp, kdensf: an unfinished kernel density command that no     */
/*  part of TDA calls (kdens itself exits with "not implemented"); plotd    */
/*  is the density estimator in use.  The bodies have been removed.          */
/*                                                                          */
/*  TDA. Program for Transition Data Analysis, written by Goetz Rohwer.     */
/*  Copyright (C) 1989,1991-97 Goetz Rohwer. All rights reserved.           */
/*                                                                          */
/*  This file is part of TDA.                                               */
/*                                                                          */
/*  TDA is free software; you can redistribute it and/or modify             */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 2 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  TDA is distributed in the hope that it will be useful,                  */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License       */
/*  along with this program (it should be in a file named COPYING),         */
/*  if not, write to the Free Software Foundation, Inc.,                    */
/*  675 Mass Ave, Cambridge, MA 02139, USA.                                 */
/*                                                                          */

#include "tda.h"
#include "t_gen.h"
#include "t_pgen.h"
#include "t_parm.h"
#include "t_gdat.h"
#include "t_var.h"
#include "t_plot.h"
#include "t_alloc.h"
#include "t_cdf.h"
#include "tda_context.h"
#include "tda_compat.h"

/*  functions in t_dens.c */

int kdens(TDAContext *ctx, int idx);
int kdensf(TDAContext *ctx, int n,float *x,int m,float a,float d,float sig,float *crit);
int kdcomp(const void *, const void *, void *);

/* ------------------------------------------------------------------------ */
/*  kdens(idx)  Density estimation.                                         */
/*  ##          kdens(v=,sig=,x=,sc=,pl=) [=fname]                          */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */
