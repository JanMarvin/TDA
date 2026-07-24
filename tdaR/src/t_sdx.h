/* t_sdx.h */

#ifndef _TSDX_H
#define _TSDX_H
#include "tda_ctx_fwd.h"

/*  functions in t_sdx.c */

int sdgen(TDAContext *ctx); 
int sdshp(TDAContext *ctx); 
int rdbf(TDAContext *ctx); 
int rcsv(TDAContext *ctx); 
int rplz(TDAContext *ctx); 
int gtopo(TDAContext *ctx); 
int sddcwp(TDAContext *ctx); 
int sdgshhs(TDAContext *ctx); 

#endif /* _TSDX_H */

/* end of t_sdx.h */


