/* t_areg.h */

#ifndef _TAREG_H
#define _TAREG_H

/*  functions in t_areg.c */

int npreg(void);
double mkern1(int opt,int n,double d,double a,double *x);
int lowess(double *x,double *y,int n,double f,int nsteps,double delta,
    double *ys,double *rw,double *res);
int rmidmean(int n,double *x,double *y,int r,double *xm,double *ym,
    float *ym1,float *ym2);
int inpreg(void);

#endif /* _TAREG_H */

/* end of t_areg.h */


