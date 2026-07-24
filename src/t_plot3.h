/* t_plot3.h */

#ifndef _TPLOT3_H
#define _TPLOT3_H

/*  functions in t_plot3.c */

double arc_to_degree(double x);
double degree_to_arc(double x);
void ps_3dprj(double x1,double x2,double x3,double *x,double *y);
void ps_3dprj1(double x1,double x2,double x3,double *x,double *y,double *z);
void ps_3dprj_inv(double x1,double x2,double x3,double *x,double *y,double *z);
void ps_3dnorm(double *x1,double *x2,double *x3);
double ps_3dist(double x1,double x2,double x3,double y1,double y2,double y3);
void ps_3dgeo(double *x1,double *x2,double *x3);
void ps_3dgeo_inv(double *x1,double *x2,double *x3);
void ps_3dplot(double x1,double x2,double x3,int opt,int geo,int nc);
void ps_3dplot1(double x,double y,int opt,int nc);

int plotp3(int typ);
int plcurv3(void);
int plsurf3(void);
int plsurf3d(void);
int pltext3(void);
int plcirc3(void);
int plglob3(void);

void pl_circ3(double x,double y,double z,double r,double s,double a,double b,
    int hide,double xp,double yp,double zp,int nc);
int pl_lon3(double x,double y,double z,double r,double lon,double lat1,double lat2,int nc);
int pl_lat3(double x,double y,double z,double r,double lat,double lon1,double lon2,int nc);
int pl_meridian(double lon,double lata,double latb,double r,int n,int nc,int lt,double lw);
int pl_parallel(double lat,double lona,double lonb,double r,int n,int nc,int lt,double lw);

#endif /* _TPLOT3_H */

/* end of t_plot3.h */









