/* t_brr.h */

#ifndef _TBRR_H
#define _TBRR_H
#include "tda_ctx_fwd.h"

/*  functions in t_brr.c */

#ifndef LINT
int brr(TDAContext *ctx);
int brr_gen(TDAContext *ctx, int ns, int nu,int *nrep,int opt);
#endif

#endif /* _TBRR_H */

/* end of t_brr.h */

