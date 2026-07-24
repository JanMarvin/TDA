/* t_plot3.h */

#ifndef _TPLOT3_H
#define _TPLOT3_H
#include "tda_ctx_fwd.h"

/*  functions in t_plot3.c */

double arc_to_degree(TDAContext *ctx, double x);
double degree_to_arc(TDAContext *ctx, double x);
void ps_3dprj(TDAContext *ctx, double x1,double x2,double x3,double *x,double *y);
void ps_3dprj1(TDAContext *ctx, double x1,double x2,double x3,double *x,double *y,double *z);
void ps_3dprj_inv(TDAContext *ctx, double x1,double x2,double x3,double *x,double *y,double *z);
void ps_3dnorm(TDAContext *ctx, double *x1,double *x2,double *x3);
double ps_3dist(TDAContext *ctx, double x1,double x2,double x3,double y1,double y2,double y3);
void ps_3dgeo(TDAContext *ctx, double *x1,double *x2,double *x3);
void ps_3dgeo_inv(TDAContext *ctx, double *x1,double *x2,double *x3);
void ps_3dplot(TDAContext *ctx, double x1,double x2,double x3,int opt,int geo,int nc);
void ps_3dplot1(TDAContext *ctx, double x,double y,int opt,int nc);

int plotp3(TDAContext *ctx, int typ);
int plcurv3(TDAContext *ctx);
int plsurf3(TDAContext *ctx);
int plsurf3d(TDAContext *ctx);
int pltext3(TDAContext *ctx);
int plcirc3(TDAContext *ctx);
int plglob3(TDAContext *ctx);

void pl_circ3(TDAContext *ctx, double x,double y,double z,double r,double s,double a,double b,
    int hide,double xp,double yp,double zp,int nc);
int pl_lon3(TDAContext *ctx, double x,double y,double z,double r,double lon,double lat1,double lat2,int nc);
int pl_lat3(TDAContext *ctx, double x,double y,double z,double r,double lat,double lon1,double lon2,int nc);
int pl_meridian(TDAContext *ctx, double lon,double lata,double latb,double r,int n,int nc,int lt,double lw);
int pl_parallel(TDAContext *ctx, double lat,double lona,double lonb,double r,int n,int nc,int lt,double lw);

#endif /* _TPLOT3_H */

/* end of t_plot3.h */









