/* t_gmin.h */

#ifndef _TGMIN_H
#define _TGMIN_H

/*  functions in t_gmin.c */

int f_min(int typ);
int prn_fsval(int n,double *x,int gmina,double *lb,double *ub,int opt);
int get_func(char *s,int opt,int *n3,int iflag,char *msg);
void prn_feval(void);
int get_flval(double *f,int n,double *x,int deriv,int opt,double *g,double *h,double *d);
int get_filval(int n,double *x,int deriv,int *ni);
int get_fival(int n,double *x,int deriv,int aflag);
int get_ifval(double *fl,double *fu,int n,double *xl,double *xu,int deriv,
    double *gl,double *gu);
int get_icval(double *fl,double *fu,int n,double *xl,double *xu);          

#endif /* _TGMIN_H */

/*  end of t_gmin.h */

