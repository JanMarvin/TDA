/* t_ass.h */

#ifndef _TASS_H
#define _TASS_H
#include "tda_ctx_fwd.h"

/*  functions in t_ass.c */

int gap(TDAContext *ctx);
int gqap(TDAContext *ctx);
int gqapd(TDAContext *ctx, int n,int niter,double alpha,double beta,int look4,int *seed,
    int *f,int *d,int *a,int *b,int *srtf,int *srtif,int *srtd,int *srtid,
    int *srtc,int *srtic,int *indexd,int *indexf,int *cost,int *fdind,
    int *opta,int *bestv,int *iter,int *ito);
int gloc(TDAContext *ctx);
int qap_w(TDAContext *ctx, int n,double *c,double *f,double *d,int *loc3n,double *wtmp);

#endif /* _TASS_H */

/* end of t_ass.h */


