/* t_lsei.h */

#ifndef _TLSEI_H
#define _TLSEI_H

/*  functions in t_lsei.c */

int lsei(double *w, int me, int ma, int mg, int n, double *x, int cov, double *rnorme, double *rnorml, int *ranka, int *ranke);
int wnnls(double *w, int nw, int me, int ma, int n, int l, double *x, double *rnorm);
int lhhfti(int m, int n, int na, int nb, double *a, double *b, double *rnorm, int mode, double *d, double *h, int *ip, double tau);
int syminv(int n, double *x);

#endif /* _TLSEI_H */

/*  end of t_lsei.h */


