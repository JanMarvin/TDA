/* t_intp.h */

#ifndef _TINTP_H
#define _TINTP_H

/*  functions in t_intp.c */

int triang(void);
int intp(void);
int intp_grid(int alg,int n,double *x,double *y,double *z,int nu,int nv,
    double *u,double *v,double *w,int ncp);
void intp_dbox2(int n,double *x,double *y);
void intp_dbox3(int n,double *x,double *y,double *z);

#endif /* _TINTP_H */

/* end of t_intp.h */


