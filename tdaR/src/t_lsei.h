/* t_lsei.h */

#ifndef _TLSEI_H
#define _TLSEI_H
#include "tda_ctx_fwd.h"

/*  functions in t_lsei.c */

int lsei(TDAContext *ctx, double *w, int me, int ma, int mg, int n, double *x, int cov, double *rnorme, double *rnorml, int *ranka, int *ranke);
int wnnls(TDAContext *ctx, double *w, int nw, int me, int ma, int n, int l, double *x, double *rnorm);
int lhhfti(TDAContext *ctx, int m, int n, int na, int nb, double *a, double *b, double *rnorm, int mode, double *d, double *h, int *ip, double tau);
int syminv(TDAContext *ctx, int n, double *x);

#endif /* _TLSEI_H */

/*  end of t_lsei.h */


