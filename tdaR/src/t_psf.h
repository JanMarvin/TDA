/* t_psf.h */

#ifndef _TPSF_H
#define _TPSF_H
#include "tda_ctx_fwd.h"

/*  functions in t_psf.c */

int psfile(TDAContext *ctx);
int psetup(TDAContext *ctx);
void setup_2dsys(TDAContext *ctx);
int ps_close(TDAContext *ctx, int opt,int scx,int kflag);
int pl_axis(TDAContext *ctx, int typ,int opt,int opt1,int lt,double lw);
int psetup3(TDAContext *ctx);
void setup_3dsys(TDAContext *ctx);
void setup_3proj(TDAContext *ctx, double lon,double lat);

#endif /* _TPSF_H */

/*  end of t_psf.h */











