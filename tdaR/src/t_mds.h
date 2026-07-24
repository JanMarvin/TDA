/* t_mds.h */

#ifndef _TMDS_H
#define _TMDS_H
#include "tda_ctx_fwd.h"

/*  functions in t_mds.c */

int dmet(TDAContext *ctx); 
int dmet1(TDAContext *ctx); 
int gpro(TDAContext *ctx);
int uds(TDAContext *ctx);
int sga(TDAContext *ctx);
int unf(TDAContext *ctx);
int mdsc(TDAContext *ctx);
int mdsm(TDAContext *ctx);
int mdsn(TDAContext *ctx);
int mdsn1(TDAContext *ctx);   /* censored counterpart, dispatched first */
int mds4_fn(TDAContext *ctx);
int mdsx(TDAContext *ctx);
int mds3_fn(TDAContext *ctx);
int mdsr(TDAContext *ctx);
int rfit(TDAContext *ctx);
int rfit1(TDAContext *ctx);

#endif /* _TMDS_H */

/* end of t_mds.h */


