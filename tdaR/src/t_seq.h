/* t_seq.h */

#ifndef _TSEQ_H
#define _TSEQ_H
#include "tda_ctx_fwd.h"

/*  functions in t_seq.c */

int t_seq(TDAContext *ctx);
int check_sdj(TDAContext *ctx, int j);
int seq_sget(TDAContext *ctx, int i,int t,int n);
int seq_getsn(TDAContext *ctx, int sn,int opt);
int seq_scheck(TDAContext *ctx, int i,int j);
void seq_afree(TDAContext *ctx, int opt);                   
void prn_ssel(TDAContext *ctx, int n);
int seq_csel(TDAContext *ctx);
void seq_dtda(TDAContext *ctx, int typ,int sn,char *fname,int noc,int nv,short *vidx);



#endif /* _TSEQ_H */

/* end of t_seq.h */


