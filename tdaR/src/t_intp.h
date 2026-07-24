/* t_intp.h */

#ifndef _TINTP_H
#define _TINTP_H
#include "tda_ctx_fwd.h"

/*  functions in t_intp.c */

int triang(TDAContext *ctx);
int intp(TDAContext *ctx);
int intp_grid(TDAContext *ctx, int alg,int n,double *x,double *y,double *z,int nu,int nv,
    double *u,double *v,double *w,int ncp);
void intp_dbox2(TDAContext *ctx, int n,double *x,double *y);
void intp_dbox3(TDAContext *ctx, int n,double *x,double *y,double *z);

#endif /* _TINTP_H */

/* end of t_intp.h */


