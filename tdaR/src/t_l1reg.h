/* t_l1reg.h */

#ifndef _TL1REG_H
#define _TL1REG_H
#include "tda_ctx_fwd.h"

/*  functions in t_l1reg.c */

int l1reg(TDAContext *ctx);
int l1regf(TDAContext *ctx, int m, int n, double *x, double *y, double *s, double *res, int *rank);

#endif /* _TL1REG_H */

/*  end of t_l1reg.h */

