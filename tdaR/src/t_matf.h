/* t_matf.h */

#ifndef _TMATF_H
#define _TMATF_H
#include "tda_ctx_fwd.h"

/*  functions in t_matf.c */

void mp_info(TDAContext *ctx);                        
int mp_alloc(TDAContext *ctx, int t,int m,int n);                 
void mp_putvar(TDAContext *ctx, int i,int n,double *mat,int nv,short *vidx,int icase); 
int mp_putlog(TDAContext *ctx, double x);
int mp_putpar(TDAContext *ctx, int n,double *x);
int mp_putmpar(TDAContext *ctx, int m,int n,double *x);
int mp_putcov(TDAContext *ctx, int n,double *x);
int get_mexpr(TDAContext *ctx, char *exp,int *row,int *col,int *ivflg);
int eval_mexpr(TDAContext *ctx, int off,int row,int col,double *mat,int ivflg,
    double *mat1,int *rsel,int *csel,int idx,double dval);



#endif /* _TMATF_H */

/*  end of t_matf.h */















