/* t_eval3.h */

#ifndef _TEVAL3_H
#define _TEVAL3_H

/*  functions in t_eval3.c                                                  */

int eval_specf(int atyp,int vn,int noc,int sl,double **buf);
void check_expr(short cnt,int *ptyp,int *nv,int *nc,int *n2,int *n3,int opt);
int evalf(int iflag);
int fnd_alloc(int opt,int deriv,int fnn,int pn,int iflag);
int deriv0(int opt,int deriv,int j,int iflag);
int deriv1(int deriv,int ix,int j,int iflag);
int deriv2(int op,int deriv,int j,double *val,int iflag);
int deriv3(int op,int deriv,int j,double *val,int iflag);
int derivf1(int op,int deriv,int j,double *val,int iflag);
int derivf2(int op,int deriv,int j,double *val,double *dval);
int derivf3(int op,int deriv,int j,double *val,int iflag);
double get_integrand(double x,int *err);
int ia_op(int typ,int j,double *buf);
void ia_mul(double a,double b,double c,double d,double *ra,double *rb);
void ia_iv(double a,double b,double c,double d,double *ra,double *rb);
void ia_sin(double a,double b,double *ra,double *rb);
void ia_cos(double a,double b,double *ra,double *rb);


#endif /* _TEVAL3_H */

/*  end of t_eval3.h */




