/* t_lin.h */

#ifndef _TLIN_H
#define _TLIN_H

/*  functions in t_lin.c */

int decomp(int n, double *x);
void dclear(int n, double *da);
void dcopy(int n,double *dx,int incx, double *dy,int incy);
double ddot(int n,double *dx,int incx,double *dy,int incy);
void dscal(int n,double da,double *dx,int incx);
double fsign(double a,double b);
void daxpy(int n,double da,double *dx,int incx,double *dy,int incy);
double twonrm(int n,double *v); 
int dpodi(double *a,int ma,int n,double *det,int job);
int dpofa(double *a,int ma,int n);
void drot(int n,double *dx,int incx,double *dy,int incy,double dc,double ds);
void drotg(double *da,double *db,double *dc,double *ds);
void dxpy(int n,int m,double *x,double *y,double *xpy);
void dzero(int n,int m,double *a,int cda);
double dnrm2(int n,double *dx,int incx);
void dqrdc(double *x,int cdx,int n,int p,double *qraux,int *jpvt,      
    double *work,int job);
void dqrsl(double *x,int cdx,int n,int k,double *qraux,double *y,double *qy,
    double *qty,double *b,double *rsd,double *xb,int job,int *info);
void dtrco(double *t,int cdt,int n,double *rcond,double *z,int job);
int idamax(int n,double *dx,int incx);
void dswap(int n,double *dx,int incx,double *dy,int incy);
double dasum(int n,double *dx,int incx);
void dtrsl(double *t,int cdt,int n,double *b,int job,int *info);                                                

#endif /* _TLIN_H */

/*  end of t_lin.h */

