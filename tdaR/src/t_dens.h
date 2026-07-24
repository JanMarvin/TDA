/* t_dens.h */

#ifndef _TDENS_H
#define _TDENS_H
#include "tda_ctx_fwd.h"

/*  functions in t_dens.c */

int kdens(TDAContext *ctx, int idx);
int kdensf(TDAContext *ctx, int n,float *x,int m,float a,float d,float sig,float *crit);

#endif /* _TDENS_H */

/* end of t_dens.h */


