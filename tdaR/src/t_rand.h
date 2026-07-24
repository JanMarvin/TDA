/* t_rand.h */

#ifndef _TRAND_H
#define _TRAND_H
#include "tda_ctx_fwd.h"

/*  functions in t_rand.c */

double random1(TDAContext *ctx);
double random2(TDAContext *ctx);
double normal(TDAContext *ctx);
double normal1(TDAContext *ctx);
int rdmn_init(TDAContext *ctx, int n,double *a);
void rdmn_free(TDAContext *ctx);
void rdmn(TDAContext *ctx);
int rdp1(TDAContext *ctx, double lambda);

/* ------------------------------------------------------------------------ */



#endif /* _TRAND_H */

/*  end of t_rand.h */

