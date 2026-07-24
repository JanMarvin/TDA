/* t_fnqr.h */

#ifndef _TFNQR_H
#define _TFNQR_H
#include "tda_ctx_fwd.h"

/*  functions in t_fnqr.c */

int fn_logit1(TDAContext *ctx);
int fn_logit2(TDAContext *ctx);
int fn_logit3(TDAContext *ctx);
int fn_logit4(TDAContext *ctx);
int fn_probit1(TDAContext *ctx);
int fn_probit2(TDAContext *ctx);
int fn_probit3(TDAContext *ctx);
int fn_prob3(TDAContext *ctx, int icase,int wave,int cat,double *par,double *prob);
int fn_probit4(TDAContext *ctx);
int fn_rmod1(TDAContext *ctx);


#endif /* _TFNQR_H */

/*  end of t_fnqr.h */


