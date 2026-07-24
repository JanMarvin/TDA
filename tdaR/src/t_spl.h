/* t_spl.h */

#ifndef _TSPL_H
#define _TSPL_H
#include "tda_ctx_fwd.h"

/*  functions in t_spl.c */

int spltena(TDAContext *ctx, int n,double *x,double *y,double *xp,double *yp,
    int mode,double slp1,double slpn,double sigma);
void spltenad(TDAContext *ctx, int n,double *x,double *y,double *xp,double *yp,
    double sigma,double t,double *tx,double *ty,int first);
int spltenc(TDAContext *ctx, int n,double *x,double *y,double *xp,double *yp,double sigma);
void spltencd(TDAContext *ctx, int n,double *x,double *y,double *xp,double *yp,
    double sigma,double t,double *tx,double *ty,int first);
int splakima(TDAContext *ctx, int mode,int l,double *x,double *y,int m,int n,double *u,double *v);
int splf(TDAContext *ctx, int m,double *x,double *y,double *w,double xa,double xb,
    int k,double s,int *nn,double * t,double *c,int maxit);
double spld(TDAContext *ctx, int n, double *t,int k,double *c,int nu,double arg,int l,double *tmp);
int interv(TDAContext *ctx, int n,double *x,double arg, int init);

#endif /* _TSPL_H */

/* end of t_spl.h */


