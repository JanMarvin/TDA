/* t_gg.h */

#ifndef _TGG_H
#define _TGG_H
#include "tda_ctx_fwd.h"

/*  functions in t_gg.c */

int giset(TDAContext *ctx);
int gcni(TDAContext *ctx);
int gcset(TDAContext *ctx);
double g_mineval(TDAContext *ctx, int n1,int *nodes1,int n2,int *nodes2,int opt);
int gcliq(TDAContext *ctx); 
int ggcliq(TDAContext *ctx);
int g_gclptr_i(TDAContext *ctx, int max);
int g_gclptr_c(TDAContext *ctx, int n,int *nodes);


#endif /* _TGG_H */

/* end of t_gg.h */


