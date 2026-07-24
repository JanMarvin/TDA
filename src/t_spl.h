/* t_spl.h */

#ifndef _TSPL_H
#define _TSPL_H

/*  functions in t_spl.c */

int spltena(int n,double *x,double *y,double *xp,double *yp,
    int mode,double slp1,double slpn,double sigma);
void spltenad(int n,double *x,double *y,double *xp,double *yp,
    double sigma,double t,double *tx,double *ty,int first);
int spltenc(int n,double *x,double *y,double *xp,double *yp,double sigma);
void spltencd(int n,double *x,double *y,double *xp,double *yp,
    double sigma,double t,double *tx,double *ty,int first);
int splakima (int mode,int l,double *x,double *y,int m,int n,double *u,double *v);
int splf (int m,double *x,double *y,double *w,double xa,double xb,
    int k,double s,int *nn,double * t,double *c,int maxit);
double spld(int n, double *t,int k,double *c,int nu,double arg,int l,double *tmp);
int interv(int n,double *x,double arg, int init);

#endif /* _TSPL_H */

/* end of t_spl.h */


