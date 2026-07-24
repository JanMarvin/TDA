/* t_dstat.h */

#ifndef _TDSTAT_H
#define _TDSTAT_H
#include "tda_ctx_fwd.h"

/*  functions in t_dstat.c */

int dstat(TDAContext *ctx);
int quant(TDAContext *ctx);
int mfreq(TDAContext *ctx, int typ);
int pcov(TDAContext *ctx, int typ);
int atab(TDAContext *ctx);
int dma(TDAContext *ctx);

#endif /* _TDSTAT_H */

/* end of t_dstat.h */


