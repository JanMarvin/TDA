/* t_fnrtb.h */

#ifndef _TFNRTB_H
#define _TFNRTB_H
#include "tda_ctx_fwd.h"

/*  functions in t_fnrtb.c */

int fn_exp(TDAContext *ctx, int resid);
int fn_exp1(TDAContext *ctx, int resid);
int fn_exp2(TDAContext *ctx, int resid);
int fn_pol(TDAContext *ctx, int resid);
int fn_pol1(TDAContext *ctx, int resid);
double get_mpol1(TDAContext *ctx, double t,int *err); 
double get_mpol2(TDAContext *ctx, double t,int *err); 
double get_mpol3(TDAContext *ctx, double t,int *err); 
int fn_gm(TDAContext *ctx, int resid);
int fn_wei(TDAContext *ctx, int resid);
int fn_sic(TDAContext *ctx, int resid);
int fn_ll(TDAContext *ctx, int resid);
int fn_ll2(TDAContext *ctx, int resid);
int fn_ln(TDAContext *ctx, int resid);
int fn_ig(TDAContext *ctx, int resid);
int fn_gam(TDAContext *ctx, int resid);

#endif /* _TFNRTB_H */

/*  end of t_fnrtb.h */




