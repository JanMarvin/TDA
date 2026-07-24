/* t_gcmd.h */

#ifndef _TGCMD_H
#define _TGCMD_H
#include "tda_ctx_fwd.h"

/*  functions in t_gcmd.c */

int mem_use(TDAContext *ctx);
int prntime(TDAContext *ctx);
int prn_txt(TDAContext *ctx);
int machp(TDAContext *ctx, int opt);
int t_parse(TDAContext *ctx);          
int data(TDAContext *ctx, int opt);


#endif /* _TGCMD_H */

/*  end of t_gcmd.h */











