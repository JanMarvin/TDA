/* t_gm.h */

#ifndef _TGM_H
#define _TGM_H

/*  functions in t_gm.c */

void g_spat(double ax,double ay,double az,double bx,double by,double bz,
    double *rx,double *ry,double *rz);
int g_intersect(double xia,double yia,double xib,double yib,
    double xja,double yja,double xjb,double yjb,double *xxa,double *yya,
    double *xxb,double *yyb);
int g_inseg(double px,double py,double xa,double ya,double xb,double yb);
double g_perp(double ux,double uy,double vx,double vy);
int g_left(double xa,double ya,double xb,double yb,double x,double y);
int g_ccw(double x0,double y0,double x1,double y1,double x2,double y2);
int g_isect(double x1,double y1,double x2,double y2,double x3,double y3,
    double x4,double y4);
double g_line_len(int n,double *x,double *y);
double g_pol_area(int n,double *x,double *y);
double g_pol_xmin(int n,double *x,double *y);
double g_pol_xmax(int n,double *x,double *y);
double g_pol_ymin(int n,double *x,double *y);
double g_pol_ymax(int n,double *x,double *y);
void g_pol_rect(int n,double *x,double *y,double *xmin,double *xmax,
    double *ymin,double *ymax);
int g_inpoly(int n, double *x,double *y,double xa,double ya);
int g_l_inpoly(int n, double *x,double *y,double xa,double ya,double xb,
    double yb);
int g_p_inpoly(int na, double *xa,double *ya,int nb,double *xb,double *yb);
int g_pol_checkp(int n,double *x,double *y,double px,double py);
int g_pol_isect(int n, double *x,double *y,double xa,double ya,
    double xb,double yb);
int g_chull(float *x,float *y,int m,int *in,int *ih,int *il);
void g_xchxy(double *x,double *y);
double g_len(double xa,double ya,double xb,double yb);
void g_poly_ccw(int n,double *x,double *y);

#endif /* _TGM_H */

/* end of t_gm.h */


