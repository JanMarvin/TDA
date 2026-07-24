/* t_svd.h */

#ifndef _TSVD_H
#define _TSVD_H

/*  functions in t_svd.c */

int eigen(int n,double *a,double *evr,double *evi,double *vec,int bal);
int svdecomp(int m,int n,double *a,double *q,double *u,double *v,int mode);
int ginv(int m,int n,double *a);
int ginv1(int m,int n,double *a,double *sv);
int evecf(int n, double *d, double *z);
int simitz(int n,int p,int km,double eps,int em,double *x,double *d,int *nit);
int eigen2(int n,double *a,double *evr,double *evi,double *vec,int mode,int itmax);
int eigen1(int n,double *a,double *evr,double *evi,double *vecr,double *veci,int *indic);

#endif /* _TSVD_H */

/* end of t_svd.h */


