/* t_com.h */

#ifndef _TCOM_H
#define _TCOM_H
#include "tda_ctx_fwd.h"

/*  functions in t_com.c */

int com(TDAContext *ctx);
int comb_n(TDAContext *ctx, int n,int *com,int first); 
int comb_nm(TDAContext *ctx, int n,int m,int *com,int first);
int comb_nm1(TDAContext *ctx, int n,int m,int *com,int first);
int perm(TDAContext *ctx, int n,int *com,int *p,int *d,int first);
int partit(TDAContext *ctx, int *k,int *p,int last);
int partit1(TDAContext *ctx, int n,int k,int *c,int *d,int last);
int pcyc(TDAContext *ctx);
int indep(TDAContext *ctx);

#endif /* _TCOM_H */

/* end of t_com.h */
int partit2(TDAContext *ctx, int n,int r,int *a,int *p,int *t,int first);
