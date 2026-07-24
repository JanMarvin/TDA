/* t_gdd.h */

#ifndef _TGDD_H
#define _TGDD_H
#include "tda_ctx_fwd.h"

/*  functions in t_gdd.c */

int gdd(TDAContext *ctx);
int gdd_alloc(TDAContext *ctx, int n); 
void gdd_free(TDAContext *ctx, int opt);
void gdd_setgt(TDAContext *ctx, int typ,int gt,int gtt);
int gdd_setnd(TDAContext *ctx);
int gdd_setptr(TDAContext *ctx, int opt);
int gdd_setap(TDAContext *ctx);
void gdd_prot(TDAContext *ctx);
int gdd_fndni(TDAContext *ctx, int k);
int check_gdj(TDAContext *ctx, int j);
double gdd_adj(TDAContext *ctx, int i,int j,int gn);
int g_suc(TDAContext *ctx, int i,int gn);
int g_pre(TDAContext *ctx, int i,int gn);
int g_loop(TDAContext *ctx, int i,int gn);
int g_getnd(TDAContext *ctx, int n,int *nd);
int g_getne(TDAContext *ctx, int n,int *ni,int *nj);
int gdd_check(TDAContext *ctx, int gn,int opt);
int gdd_check2(TDAContext *ctx, int gn,int gn1,int opt);
int gdd_tcheck(TDAContext *ctx, int typ,int dir,int val);
int gdd_ei(TDAContext *ctx, int k);             
int gdd_ej(TDAContext *ctx, int k);             
double gdd_ev(TDAContext *ctx, int k,int gn);      
int gdd_node(TDAContext *ctx, int i);             
void n_info(TDAContext *ctx, int i,int n);
void n_info_e(TDAContext *ctx);
int gcd(TDAContext *ctx);


#endif /* _TGDD_H */

/* end of t_gdd.h */


