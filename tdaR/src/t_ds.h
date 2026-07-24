/* t_ds.h */

#ifndef _TDS_H
#define _TDS_H
#include "tda_ctx_fwd.h"

/*  functions in t_ds.c */

int mean(TDAContext *ctx, int n,short *vinum,double *mean);
int std(TDAContext *ctx, int n,short *vinum,double *mean,double *x,int opt);
int corr(TDAContext *ctx, int n,short *vinum,double *mean,double *std,double *x);
void rcorr(TDAContext *ctx, int n,short *vinum,double *x);
int cov(TDAContext *ctx, int n,short *vinum,double *m,double *x);
double xmean(TDAContext *ctx, int n,double *x,double *wt);
double xstd(TDAContext *ctx, int n,double *x,double *wt,double *std);
double quantf(TDAContext *ctx, int n,double *x,double q);

#endif /* _TDS_H */

/* end of t_ds.h */


