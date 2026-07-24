/* t_clip.h */

#ifndef _TCLIP_H
#define _TCLIP_H
#include "tda_ctx_fwd.h"

/*  functions in t_clip.c */

int sdclip(TDAContext *ctx);
void sdclip_err_msg(TDAContext *ctx, int n,int id);
void sdclip_write_header(TDAContext *ctx, FILE *fd,int id,int typ,int n,int sdid,int id1);
void sdclip_write_point(TDAContext *ctx, FILE *fd,double x,double y);

int sdclip_poly_rec(TDAContext *ctx, int n,double *x,double *y,double *sx,double *sy,int *nf,
    int *d,double xmin,double ymin,double xmax,double ymax);

int clip_simp(TDAContext *ctx, int n,double *x,double *y,double xmin,double ymin,
    double xmax,double ymax);
int line_intersect(double xa1,double ya1,double xa2,double ya2,
    double xb1,double yb1,double xb2,double yb2,
    double *sax,double *say,double *sbx,double *sby);

#endif /* _TCLIP_H */

/* end of t_clip.h */
int clip_line(TDAContext *ctx, double xa,double ya,double xb,double yb,double xmin,double ymin,double xmax,double ymax,double *sax,double *say,double *sbx,double *sby,int *ca,int *cb,int *ra,int *rb,int *dir);
