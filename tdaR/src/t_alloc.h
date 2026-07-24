/* t_alloc.h */

#ifndef _TALLOC_H
#define _TALLOC_H
#include "tda_ctx_fwd.h"

/*  functions in t_alloc.c */

void a_clean(TDAContext *ctx);
int alloc_actmp(TDAContext *ctx, int n);
int alloc_actmp1(TDAContext *ctx, int n);
int alloc_acx(TDAContext *ctx, int n);
int alloc_acy(TDAContext *ctx, int n);
int alloc_acy1(TDAContext *ctx, int n);
int alloc_acz(TDAContext *ctx, int n);
int alloc_act(TDAContext *ctx, int n);
int alloc_acu(TDAContext *ctx, int n);
int alloc_acv(TDAContext *ctx, int n);
int alloc_acw(TDAContext *ctx, int n);
int alloc_acxf(TDAContext *ctx, int n);
int alloc_acyf(TDAContext *ctx, int n);
int alloc_aczf(TDAContext *ctx, int n);
int alloc_acuf(TDAContext *ctx, int n);
int alloc_acvf(TDAContext *ctx, int n);
int alloc_acm(TDAContext *ctx, int n);
int alloc_acn(TDAContext *ctx, int n);
int alloc_aci(TDAContext *ctx, int n);
int alloc_acj(TDAContext *ctx, int n);
int alloc_acr(TDAContext *ctx, int n);
int alloc_acs(TDAContext *ctx, int n);
int alloc_ack(TDAContext *ctx, int n);
int alloc_acl(TDAContext *ctx, int n);
int alloc_acc(TDAContext *ctx, int n);
int alloc_acd(TDAContext *ctx, int n);
int alloc_ace(TDAContext *ctx, int n);
int alloc_acf(TDAContext *ctx, int n);
int alloc_acns(TDAContext *ctx, int n);
int alloc_acms(TDAContext *ctx, int n);
int alloc_acptr(TDAContext *ctx, int n);
int alloc_aciptr(TDAContext *ctx, int n);
int alloc_acjptr(TDAContext *ctx, int n);


#endif /* _TALLOC_H */

/* end of t_alloc.h */


