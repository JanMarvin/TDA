/* t_gm.h */

#ifndef _TGM_H
#define _TGM_H
#include "tda_ctx_fwd.h"

/*  functions in t_gm.c */

void g_spat(TDAContext *ctx, double ax,double ay,double az,double bx,double by,double bz,
    double *rx,double *ry,double *rz);
int g_intersect(TDAContext *ctx, double xia,double yia,double xib,double yib,
    double xja,double yja,double xjb,double yjb,double *xxa,double *yya,
    double *xxb,double *yyb);
int g_inseg(TDAContext *ctx, double px,double py,double xa,double ya,double xb,double yb);
double g_perp(TDAContext *ctx, double ux,double uy,double vx,double vy);
int g_left(TDAContext *ctx, double xa,double ya,double xb,double yb,double x,double y);
int g_ccw(TDAContext *ctx, double x0,double y0,double x1,double y1,double x2,double y2);
int g_isect(TDAContext *ctx, double x1,double y1,double x2,double y2,double x3,double y3,
    double x4,double y4);
double g_line_len(TDAContext *ctx, int n,double *x,double *y);
double g_pol_area(TDAContext *ctx, int n,double *x,double *y);
double g_pol_xmin(TDAContext *ctx, int n,double *x,double *y);
double g_pol_xmax(TDAContext *ctx, int n,double *x,double *y);
double g_pol_ymin(TDAContext *ctx, int n,double *x,double *y);
double g_pol_ymax(TDAContext *ctx, int n,double *x,double *y);
void g_pol_rect(TDAContext *ctx, int n,double *x,double *y,double *xmin,double *xmax,
    double *ymin,double *ymax);
int g_inpoly(TDAContext *ctx, int n, double *x,double *y,double xa,double ya);
int g_l_inpoly(TDAContext *ctx, int n, double *x,double *y,double xa,double ya,double xb,
    double yb);
int g_p_inpoly(TDAContext *ctx, int na, double *xa,double *ya,int nb,double *xb,double *yb);
int g_pol_checkp(TDAContext *ctx, int n,double *x,double *y,double px,double py);
int g_pol_isect(TDAContext *ctx, int n, double *x,double *y,double xa,double ya,
    double xb,double yb);
int g_chull(TDAContext *ctx, float *x,float *y,int m,int *in,int *ih,int *il);
void g_xchxy(TDAContext *ctx, double *x,double *y);
double g_len(TDAContext *ctx, double xa,double ya,double xb,double yb);
void g_poly_ccw(TDAContext *ctx, int n,double *x,double *y);

#endif /* _TGM_H */

/* end of t_gm.h */


