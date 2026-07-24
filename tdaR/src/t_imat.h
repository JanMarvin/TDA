/* t_imat.h */

#ifndef _TIMAT_H
#define _TIMAT_H
#include "tda_ctx_fwd.h"

/*  functions in t_imat.c */

int m_midf(TDAContext *ctx, char *cmd);
int m_midf1(TDAContext *ctx, char *cmd);
int m_midf2(TDAContext *ctx, char *cmd);
int m_midf3(TDAContext *ctx, char *cmd);
int gmin(TDAContext *ctx);
int range(TDAContext *ctx);
int idf(TDAContext *ctx);
int sddf(TDAContext *ctx);
int iddf(TDAContext *ctx);
int imean(TDAContext *ctx);
int ivar(TDAContext *ctx);
int icov2(TDAContext *ctx, int typ);
int igini(TDAContext *ctx);
int ivar1(TDAContext *ctx);
int alloc_par(TDAContext *ctx, int n,int opt);
int alloc_list(TDAContext *ctx, int n,int m);
int igmin(TDAContext *ctx, int opt,int narg,int nbmax,double *par,double *lb,double *ub,int mxit,
    double tolbw,double tolfd,double tolfe,int gc,int typ);
int igmin_res(TDAContext *ctx, double tolfe,int opt);
void igmin_cpar(TDAContext *ctx, int n,double *x,double *par);

void i_add(TDAContext *ctx, double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_sub(TDAContext *ctx, double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_mul(TDAContext *ctx, double xl,double xh,double yl,double yh,double *rl,double *rh);
int  i_div(TDAContext *ctx, double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_abs(TDAContext *ctx, double xl,double xh,double *rl,double *rh);
void i_max(TDAContext *ctx, double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_min(TDAContext *ctx, double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_square(TDAContext *ctx, double xl,double xh,double *rl,double *rh);
void i_sqrt(TDAContext *ctx, double xl,double xh,double *rl,double *rh);
double ivarf(TDAContext *ctx, int n,double x,double *xl,double *xh);






#endif /* _TIMAT_H */

/*  end of t_imat.h */















