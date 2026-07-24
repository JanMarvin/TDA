/* t_gf.h */

#ifndef _TGF_H
#define _TGF_H
#include "tda_ctx_fwd.h"

/*  functions in t_gf.c */

void clear_npflg(TDAContext *ctx);
void prn_npflg(TDAContext *ctx);
int imax(TDAContext *ctx, int i,int j);
int imin(TDAContext *ctx, int i,int j);
int iabs(TDAContext *ctx, int x);
double dmax(TDAContext *ctx, double x,double y);
double dmin(TDAContext *ctx, double x,double y);
double dsign(TDAContext *ctx, double x);
int isign(TDAContext *ctx, int x);
double rexp(TDAContext *ctx, double x);
double rlog(TDAContext *ctx, double x);
double loggam(TDAContext *ctx, double x,int *err);
double icgam(TDAContext *ctx, double x, double alpha,int *err);
double beta(TDAContext *ctx, double a, double b,int *err);
double incbeta(TDAContext *ctx, double x, double aa, double bb,int *err);
double digam(TDAContext *ctx, double x,int *err);
double trigam(TDAContext *ctx, double x,int *err);
int icdgam(TDAContext *ctx, double x,double p,double *d);
double bincoeff(TDAContext *ctx, int n,int  m);
int l_poisson(TDAContext *ctx, double theta,double k,double *val,double *val1,double *val2);
int l_negbin(TDAContext *ctx, double alpha,double gamma,double k,int opt,double *d);

/* ------------------------------------------------------------------------ */

#endif

