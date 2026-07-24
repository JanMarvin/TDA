/* t_ds.h */

#ifndef _TDS_H
#define _TDS_H

/*  functions in t_ds.c */

int mean(int n,short *vinum,double *mean);
int std(int n,short *vinum,double *mean,double *x,int opt);
int corr(int n,short *vinum,double *mean,double *std,double *x);
void rcorr(int n,short *vinum,double *x);
int cov(int n,short *vinum,double *m,double *x);
double xmean(int n,double *x,double *wt);
double xstd(int n,double *x,double *wt,double *std);
double quantf(int n,double *x,double q);

#endif /* _TDS_H */

/* end of t_ds.h */


