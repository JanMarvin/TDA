/* t_areg.h */

#ifndef _TAREG_H
#define _TAREG_H
#include "tda_ctx_fwd.h"

/*  functions in t_areg.c */

int npreg(TDAContext *ctx);
double mkern1(TDAContext *ctx, int opt,int n,double d,double a,double *x);
int lowess(TDAContext *ctx, double *x,double *y,int n,double f,int nsteps,double delta,
    double *ys,double *rw,double *res);
int rmidmean(TDAContext *ctx, int n,double *x,double *y,int r,double *xm,double *ym,
    float *ym1,float *ym2);
int inpreg(TDAContext *ctx);

#endif /* _TAREG_H */

/* end of t_areg.h */


