/* t_sort.h */

#ifndef _TSORT_H
#define _TSORT_H

/*  functions in t_sort.c */

int sortd(int n,double *x,int opt);
int sortd2(int n,double *x,double *y);
int sortd3(int n,double *x,double *y,double *z);
int vsort(int n,short *vidx,int opt,int cflg,int pflg);   
int sorti(int n,int *x,int opt);
int sorti2(int n,int *x,int *y);
int sortdpi(int n,int *x,int *ptr);
int sortdpi2(int n,int *x,int *y,int *ptr);
int sortdp(int n,double *x,int *ptr);
int sortdp2(int n,double *x,double *y,int *ptr);
int sortdp2f(int n,float *x,float *y,int *ptr);
int sortdp4f(int n,float *xa,float *ya,float *xb,float *yb,int *ptr);
int sortdp2a(int n,double *x,short *c,int *ptr,int opt);
int sortdpn(int m,int n,double *x,int nc,int *col,int *ptr);
int sortdpn1(int m,int n,double *x,int nc,int *col,int *ptr,short *cen);
int sortdp2c(int n,double *x,double *y,int *ptr,double xmin,double ymin,
    double xmax,double ymax);
int s2ccomp(const void *arg1,const void *arg2);

extern int *VSORTPtr;   /* pointer to sorted data matrix cases */

#endif /* _TSORT_H */

/* end of t_sort.h */


