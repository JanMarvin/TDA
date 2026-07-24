/* t_sort.h */

#ifndef _TSORT_H
#define _TSORT_H
#include "tda_ctx_fwd.h"

/*  functions in t_sort.c */

int sortd(TDAContext *ctx, int n,double *x,int opt);
int sortd2(TDAContext *ctx, int n,double *x,double *y);
int sortd3(TDAContext *ctx, int n,double *x,double *y,double *z);
int vsort(TDAContext *ctx, int n,short *vidx,int opt,int cflg,int pflg);   
int sorti(TDAContext *ctx, int n,int *x,int opt);
int sorti2(TDAContext *ctx, int n,int *x,int *y);
int sortdpi(TDAContext *ctx, int n,int *x,int *ptr);
int sortdpi2(TDAContext *ctx, int n,int *x,int *y,int *ptr);
int sortdp(TDAContext *ctx, int n,double *x,int *ptr);
int sortdp2(TDAContext *ctx, int n,double *x,double *y,int *ptr);
int sortdp2f(TDAContext *ctx, int n,float *x,float *y,int *ptr);
int sortdp4f(TDAContext *ctx, int n,float *xa,float *ya,float *xb,float *yb,int *ptr);
int sortdp2a(TDAContext *ctx, int n,double *x,short *c,int *ptr,int opt);
int sortdpn(TDAContext *ctx, int m,int n,double *x,int nc,int *col,int *ptr);
int sortdpn1(TDAContext *ctx, int m,int n,double *x,int nc,int *col,int *ptr,short *cen);
int sortdp2c(TDAContext *ctx, int n,double *x,double *y,int *ptr,double xmin,double ymin,
    double xmax,double ymax);
int s2ccomp(TDAContext *ctx, const void *arg1,const void *arg2);


#endif /* _TSORT_H */

/* end of t_sort.h */


