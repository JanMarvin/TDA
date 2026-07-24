/* t_cplot.h */

#ifndef _TCPLOT_H
#define _TCPLOT_H
#include "tda_ctx_fwd.h"

/*  functions in t_cplot.c */

int pl_plotc(TDAContext *ctx);
void cont_plot(TDAContext *ctx, int nlev,double *flev,int opt);
int pl_plotcm(TDAContext *ctx);
int pl_plotr(TDAContext *ctx);


#endif /* _TCPLOT_H */

/* end of t_cplot.h */









