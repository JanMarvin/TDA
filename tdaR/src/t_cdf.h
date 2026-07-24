/* t_cdf.h */

#ifndef _TCDF_H
#define _TCDF_H
#include "tda_ctx_fwd.h"

/*  functions in t_cdf.c */

double dnf(TDAContext *ctx, double x);
double cdnf(TDAContext *ctx, double z);
double cdnf1(TDAContext *ctx, double z);
double cdnif1(TDAContext *ctx, double p);
double cdchif(TDAContext *ctx, double x, int df);
double cdff(TDAContext *ctx, double x, int m, int n);
double cdtf(TDAContext *ctx, double t, int df);
double cdnm(TDAContext *ctx, double x);
int dmv(TDAContext *ctx, int m,int k,double *h,double *r,double *prob,double eps,double *err);
double bivn(TDAContext *ctx, double ah,double ak,double r,int *err);

/* ------------------------------------------------------------------------ */

#endif /* _TCDF_H */

/*  end of t_cdf.h */

