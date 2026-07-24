/* t_lin.h */

#ifndef _TLIN_H
#define _TLIN_H
#include "tda_ctx_fwd.h"

/*  functions in t_lin.c */

int decomp(TDAContext *ctx, int n, double *x);
void dclear(TDAContext *ctx, int n, double *da);
void dcopy(TDAContext *ctx, int n,double *dx,int incx, double *dy,int incy);
double ddot(TDAContext *ctx, int n,double *dx,int incx,double *dy,int incy);
void dscal(TDAContext *ctx, int n,double da,double *dx,int incx);
double fsign(TDAContext *ctx, double a,double b);
void daxpy(TDAContext *ctx, int n,double da,double *dx,int incx,double *dy,int incy);
double twonrm(TDAContext *ctx, int n,double *v); 
int dpodi(TDAContext *ctx, double *a,int ma,int n,double *det,int job);
int dpofa(TDAContext *ctx, double *a,int ma,int n);
void drot(TDAContext *ctx, int n,double *dx,int incx,double *dy,int incy,double dc,double ds);
void drotg(TDAContext *ctx, double *da,double *db,double *dc,double *ds);
void dxpy(TDAContext *ctx, int n,int m,double *x,double *y,double *xpy);
void dzero(TDAContext *ctx, int n,int m,double *a,int cda);
double dnrm2(TDAContext *ctx, int n,double *dx,int incx);
void dqrdc(TDAContext *ctx, double *x,int cdx,int n,int p,double *qraux,int *jpvt,      
    double *work,int job);
void dqrsl(TDAContext *ctx, double *x,int cdx,int n,int k,double *qraux,double *y,double *qy,
    double *qty,double *b,double *rsd,double *xb,int job,int *info);
void dtrco(TDAContext *ctx, double *t,int cdt,int n,double *rcond,double *z,int job);
int idamax(TDAContext *ctx, int n,double *dx,int incx);
void dswap(TDAContext *ctx, int n,double *dx,int incx,double *dy,int incy);
double dasum(TDAContext *ctx, int n,double *dx,int incx);
void dtrsl(TDAContext *ctx, double *t,int cdt,int n,double *b,int job,int *info);                                                

#endif /* _TLIN_H */

/*  end of t_lin.h */

