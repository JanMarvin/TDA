/* t_rate.h */

#ifndef _TRATE_H
#define _TRATE_H
#include "tda_ctx_fwd.h"

/*  functions in t_rate.c */

int rate(TDAContext *ctx);
void prn_resid1(TDAContext *ctx, int icase, int org, int des, double ts, double tf,
    double rate, double surv, double resid, double wt);

/* ------------------------------------------------------------------------ */




#endif /* _TRATE_H */

/* end of t_rate.h */


