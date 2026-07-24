/* t_smo.h */

#ifndef _TSMO_H
#define _TSMO_H
#include "tda_ctx_fwd.h"

/*  functions in t_smo.c */

int sma(TDAContext *ctx);
int smd(TDAContext *ctx);
int s_sma(TDAContext *ctx, int n,double *x,int nw,double *w,int r,int opt);
int s_smd(TDAContext *ctx, char *sm,int opt,int n,double *y);
void nnsmooth(TDAContext *ctx, int n,double *x,double *y,int r);
void nneighbor(TDAContext *ctx, int n,double *x,int r,int k,int *k1,int *k2);
int spl(TDAContext *ctx);
int scplot(TDAContext *ctx);

#endif /* _TSMO_H */

/* end of t_smo.h */


