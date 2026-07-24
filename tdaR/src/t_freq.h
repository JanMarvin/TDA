/* t_freq.h */

#ifndef _TFREQ_H
#define _TFREQ_H
#include "tda_ctx_fwd.h"

/*  functions in t_freq.c */

int cfreq(TDAContext *ctx, int n, int m, int *x, int mcel, int *value, int *freq, int *bf,
    int wflg, double *wt, double *wfreq);
int cfreq1(TDAContext *ctx, int n,short *idx,int mcel,int *value,int *freq,int *bf,
    int widx,double *wfreq);

#endif /* _TFREQ_H */

/*  end of t_freq.h */

