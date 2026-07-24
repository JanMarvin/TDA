/* t_mes.h */

#ifndef _TMES_H
#define _TMES_H

/*  functions in t_mes.c */

int ghd(void);
int ghd1(void);
int conj(void);
int mreg(void);
int mreg_proj(int n,double *x);
int mreg_proj1(int n,int m,int *t,int *p,int *pp,double *x,double *y);
int mreg_proj2(int n,int m,int *t,double *x,double *y);
int mreg_ties(int n,double *y,int *t);
int nmca(void);

#endif /* _TMES_H */

/* end of t_mes.h */


