/* t_mdat.h */

#ifndef _TMDAT_H
#define _TMDAT_H
#include "tda_ctx_fwd.h"

/*  functions in t_mdat.c */

int clear(TDAContext *ctx);
void clear_a(TDAContext *ctx);
int clearnl(TDAContext *ctx);
int tsel(TDAContext *ctx);       
void tsel_off(TDAContext *ctx, int opt);       
int cwt(TDAContext *ctx);       
int wr_sys(TDAContext *ctx);
int rd_sys(TDAContext *ctx);
int dblock(TDAContext *ctx);          
int dblock_alloc(TDAContext *ctx, int n);   
int repsel(TDAContext *ctx);          
int repsel_off(TDAContext *ctx);   
int repsel_alloc(TDAContext *ctx, int n);   


#endif /* _TMDAT_H */

/*  end of t_mdat.h */











