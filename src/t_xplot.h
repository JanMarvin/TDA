/* t_xplot.h */

#ifndef _TXPLOT_H
#define _TXPLOT_H

/*  functions in t_xplot.c */

int xplot(void);
int alloc_scx(int n);
void x_psetup(double xa,double xb,double ya,double yb);
int x_getn(double x);
int x_getaxval(double xmin,double xmax,double *xa,double *xb,int opt);
void x_plaxis(int opt,int nx,double xa,double xb);
void x_scplot(int n,int ns,int opt,int clip);
int x_cchk(double xmin,double xmax,double ymin,double ymax);
int xcheck(void);
int xlog(void);
int xlog1(void);
int xlogp(void);
int xconh(void);
int xconhp(int g,int lt);
int xreg(void);
int xplotf(void);
int xopen(void);
int xdelete(void);
int xdens(void);
int xfunc(void);

#endif /* _TXPLOT_H */

/*  end of t_xplot.h */











