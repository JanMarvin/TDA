/* t_gf.h */

#ifndef _TGF_H
#define _TGF_H

/*  functions in t_gf.c */

void clear_npflg(void);
void prn_npflg(void);
int imax(int i,int j);
int imin(int i,int j);
int iabs(int x);
double dmax(double x,double y);
double dmin(double x,double y);
double dsign(double x);
int isign(int x);
double rexp(double x);
double rlog(double x);
double loggam(double x,int *err);
double icgam(double x, double alpha,int *err);
double beta(double a, double b,int *err);
double incbeta(double x, double aa, double bb,int *err);
double digam(double x,int *err);
double trigam(double x,int *err);
int icdgam(double x,double p,double *d);
double bincoeff(int n,int  m);
int l_poisson(double theta,double k,double *val,double *val1,double *val2);
int l_negbin(double alpha,double gamma,double k,int opt,double *d);

/* ------------------------------------------------------------------------ */
extern short NPFlgs[];      /* Numerical problems flags                     */

#endif

