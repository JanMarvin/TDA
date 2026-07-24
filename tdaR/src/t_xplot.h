/* t_xplot.h */

#ifndef _TXPLOT_H
#define _TXPLOT_H
#include "tda_ctx_fwd.h"

/*  functions in t_xplot.c */

int xplot(TDAContext *ctx);
int alloc_scx(TDAContext *ctx, int n);
void x_psetup(TDAContext *ctx, double xa,double xb,double ya,double yb);
int x_getn(TDAContext *ctx, double x);
int x_getaxval(TDAContext *ctx, double xmin,double xmax,double *xa,double *xb,int opt);
void x_plaxis(TDAContext *ctx, int opt,int nx,double xa,double xb);
void x_scplot(TDAContext *ctx, int n,int ns,int opt,int clip);
int x_cchk(TDAContext *ctx, double xmin,double xmax,double ymin,double ymax);
int xcheck(TDAContext *ctx);
int xlog(TDAContext *ctx);
int xlog1(TDAContext *ctx);
int xlogp(TDAContext *ctx);
int xconh(TDAContext *ctx);
int xconhp(TDAContext *ctx, int g,int lt);
int xreg(TDAContext *ctx);
int xplotf(TDAContext *ctx);
int xopen(TDAContext *ctx);
int xdelete(TDAContext *ctx);
int xdens(TDAContext *ctx);
int xfunc(TDAContext *ctx);

#endif /* _TXPLOT_H */

/*  end of t_xplot.h */











