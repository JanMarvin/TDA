/* t_pdat.h */

#ifndef _TPDAT_H
#define _TPDAT_H
#include "tda_ctx_fwd.h"

/*  functions in t_pdat.c */

int pdata(TDAContext *ctx);
void dtda(TDAContext *ctx, char *fname,int noc,int nv,short *vidx,int nq);
void dspss(TDAContext *ctx, char *fname,int nv,short *vidx,int nq);
int pdatr(TDAContext *ctx);
int pdatd(TDAContext *ctx);

#endif /* _TPDAT_H */

/*  end of t_pdat.h */











