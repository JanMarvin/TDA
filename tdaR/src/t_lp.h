/* t_lp.h */

#ifndef _TLP_H
#define _TLP_H
#include "tda_ctx_fwd.h"

/*  functions in t_lp.c */

int lpf1(TDAContext *ctx, int m,int n,int p,double *tab,double *x,double *y,double *z,double tdd);
int lpi(TDAContext *ctx, int n,int m,int nest,int *a0,int *b0,int *opts,int *itmp,int *vmax);
int lpia(TDAContext *ctx, int n,int m,int nest,int *an,int *an0,int **av,int **ar,
         int *b0,int *opts,int *itmp,int *vmax);

#endif /* _TLP_H */

/* end of t_lp.h */


