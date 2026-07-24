/* t_mes.h */

#ifndef _TMES_H
#define _TMES_H
#include "tda_ctx_fwd.h"

/*  functions in t_mes.c */

int ghd(TDAContext *ctx);
int ghd1(TDAContext *ctx);
int tda_conj(TDAContext *ctx);
int mreg(TDAContext *ctx);
int mreg_proj(TDAContext *ctx, int n,double *x);
int mreg_proj1(TDAContext *ctx, int n,int m,int *t,int *p,int *pp,double *x,double *y);
int mreg_proj2(TDAContext *ctx, int n,int m,int *t,double *x,double *y);
int mreg_ties(TDAContext *ctx, int n,double *y,int *t);
int nmca(TDAContext *ctx);

#endif /* _TMES_H */

/* end of t_mes.h */
