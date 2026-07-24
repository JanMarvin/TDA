/* t_cdf.h */

#ifndef _TCDF_H
#define _TCDF_H

/*  functions in t_cdf.c */

double dnf(double x);
double cdnf(double z);
double cdnf1(double z);
double cdnif1(double p);
double cdchif(double x, int df);
double cdff(double x, int m, int n);
double cdtf(double t, int df);
double cdnm(double x);
int dmv(int m,int k,double *h,double *r,double *prob,double eps,double *err);
double bivn(double ah,double ak,double r,int *err);

/* ------------------------------------------------------------------------ */
extern double DMVArg[];         /* arguments for dmv()                      */

#endif /* _TCDF_H */

/*  end of t_cdf.h */

