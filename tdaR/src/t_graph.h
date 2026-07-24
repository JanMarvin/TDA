/* t_graph.h */

#ifndef _TGRAPH_H
#define _TGRAPH_H
#include "tda_ctx_fwd.h"

/*  functions in t_graph.c */

int gni(TDAContext *ctx);
int gdln(TDAContext *ctx);
int gdp(TDAContext *ctx);
int gda(TDAContext *ctx);
int gdu(TDAContext *ctx);
int gde(TDAContext *ctx);
int gdcset(TDAContext *ctx);
int gsym(TDAContext *ctx);
int gdot(TDAContext *ctx);
int gtcl(TDAContext *ctx);
int gsort(TDAContext *ctx);
int gdcyc(TDAContext *ctx);
int gev(TDAContext *ctx); 
void g_op1(TDAContext *ctx, int n,double *z,double *w);
int gsp(TDAContext *ctx);
int tnet(TDAContext *ctx);

#endif /* _TGRAPH_H */

/* end of t_graph.h */


