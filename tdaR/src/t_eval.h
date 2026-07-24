/* t_eval.h */

#ifndef _TEVAL_H
#define _TEVAL_H
#include "tda_ctx_fwd.h"

/*  functions in t_eval.c                                                  */

int alloc_est(TDAContext *ctx, int n,int opt);
void parse(TDAContext *ctx, short cnt,int *typ,double *val,int *idx);
int v_parse(TDAContext *ctx, char *s,int iflag);
int v_search(TDAContext *ctx, char *p,int *len,int *narg); 
int v_search1(TDAContext *ctx, char *p,int *len);
char *skip_const(TDAContext *ctx, char *p);                   
void prn_emsg1(TDAContext *ctx, int n);
void prn_emsg2(TDAContext *ctx, int n);

/* ------------------------------------------------------------------------ */



















#endif /* _TEVAL_H */

/*  end of t_eval.h */




