/* t_eval3.h */

#ifndef _TEVAL3_H
#define _TEVAL3_H
#include "tda_ctx_fwd.h"

/*  functions in t_eval3.c                                                  */

int eval_specf(TDAContext *ctx, int atyp,int vn,int noc,int sl,double **buf);
void check_expr(TDAContext *ctx, short cnt,int *ptyp,int *nv,int *nc,int *n2,int *n3,int opt);
int evalf(TDAContext *ctx, int iflag);
int fnd_alloc(TDAContext *ctx, int opt,int deriv,int fnn,int pn,int iflag);
int deriv0(TDAContext *ctx, int opt,int deriv,int j,int iflag);
int deriv1(TDAContext *ctx, int deriv,int ix,int j,int iflag);
int deriv2(TDAContext *ctx, int op,int deriv,int j,double *val,int iflag);
int deriv3(TDAContext *ctx, int op,int deriv,int j,double *val,int iflag);
int derivf1(TDAContext *ctx, int op,int deriv,int j,double *val,int iflag);
int derivf2(TDAContext *ctx, int op,int deriv,int j,double *val,double *dval);
int derivf3(TDAContext *ctx, int op,int deriv,int j,double *val,int iflag);
double get_integrand(TDAContext *ctx, double x,int *err);
int ia_op(TDAContext *ctx, int typ,int j,double *buf);
void ia_mul(TDAContext *ctx, double a,double b,double c,double d,double *ra,double *rb);
void ia_iv(TDAContext *ctx, double a,double b,double c,double d,double *ra,double *rb);
void ia_sin(TDAContext *ctx, double a,double b,double *ra,double *rb);
void ia_cos(TDAContext *ctx, double a,double b,double *ra,double *rb);


#endif /* _TEVAL3_H */

/*  end of t_eval3.h */




