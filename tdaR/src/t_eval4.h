/* t_eval4.h */

#ifndef _TEVAL4_H
#define _TEVAL4_H
#include "tda_ctx_fwd.h"

/*  functions in t_eval4.c                                                  */

int alloc_mex(TDAContext *ctx, int n,int opt);
int mparse(TDAContext *ctx);          
int alloc_mex_eval(TDAContext *ctx, int opt,int n,int m);
int mex_eval(TDAContext *ctx, char *s);


#endif /* _TEVAL_H */

/*  end of t_eval.h */




