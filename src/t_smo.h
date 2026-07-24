/* t_smo.h */

#ifndef _TSMO_H
#define _TSMO_H

/*  functions in t_smo.c */

int sma(void);
int smd(void);
int s_sma(int n,double *x,int nw,double *w,int r,int opt);
int s_smd(char *sm,int opt,int n,double *y);
void nnsmooth(int n,double *x,double *y,int r);
void nneighbor(int n,double *x,int r,int k,int *k1,int *k2);
int spl(void);
int scplot(void);

#endif /* _TSMO_H */

/* end of t_smo.h */


