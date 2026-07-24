/* t_plot.h */

#ifndef _TPLOT_H
#define _TPLOT_H
#include "tda_ctx_fwd.h"

/*  functions in t_plot.c */

int check_pcmd(TDAContext *ctx, int opt,int dim);              
int check_ps(TDAContext *ctx, int dim);              
void upd_bbox(TDAContext *ctx, int opt,double x,double y); 
void set_clip(TDAContext *ctx); 
void ps_ltyp(TDAContext *ctx, int typ);
void ps_lwidth(TDAContext *ctx, double lw);
void ps_fill(TDAContext *ctx, double g); 
void ps_2dplot(TDAContext *ctx, double x,double y, int opt);
double ps_2dx(TDAContext *ctx, double x);
double ps_2dy(TDAContext *ctx, double y);
void plot_str(TDAContext *ctx, double x,double y,char *s, double siz, int opt,int adj,
    int r,int xopt,int cflag);
void ps_nlab(TDAContext *ctx, int n,double x,double y,double fs,int opt);
void ps_sym(TDAContext *ctx, int typ, double px, double py,double siz);
int ps_grid(TDAContext *ctx, int typ);
int pl_label(TDAContext *ctx, int typ);
int pl_text(TDAContext *ctx);
int pl_frame(TDAContext *ctx);
int pl_rec(TDAContext *ctx);
int pl_plotp(TDAContext *ctx, int typ);
int pl_plotm(TDAContext *ctx);
int pl_plotf(TDAContext *ctx);
int pl_ploth(TDAContext *ctx);
int pl_plotd(TDAContext *ctx);
int pl_ploto(TDAContext *ctx, int typ);
int pl_plotk(TDAContext *ctx); 
int pl_plotch(TDAContext *ctx);
int pl_plots(TDAContext *ctx, int typ);
int pl_plg(TDAContext *ctx);


/*--------------------------------------------------------------------------*/




                         







#endif /* _TPLOT_H */

/* end of t_plot.h */









