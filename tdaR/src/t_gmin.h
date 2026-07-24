/* t_gmin.h */

#ifndef _TGMIN_H
#define _TGMIN_H
#include "tda_ctx_fwd.h"

/*  functions in t_gmin.c */

int f_min(TDAContext *ctx, int typ);
int prn_fsval(TDAContext *ctx, int n,double *x,int gmina,double *lb,double *ub,int opt);
int get_func(TDAContext *ctx, char *s,int opt,int *n3,int iflag,char *msg);
void prn_feval(TDAContext *ctx);
int get_flval(TDAContext *ctx, double *f,int n,double *x,int deriv,int opt,double *g,double *h,double *d);
int get_filval(TDAContext *ctx, int n,double *x,int deriv,int *ni);
int get_fival(TDAContext *ctx, int n,double *x,int deriv,int aflag);
int get_ifval(TDAContext *ctx, double *fl,double *fu,int n,double *xl,double *xu,int deriv,
    double *gl,double *gu);
int get_icval(TDAContext *ctx, double *fl,double *fu,int n,double *xl,double *xu);          

#endif /* _TGMIN_H */

/*  end of t_gmin.h */

