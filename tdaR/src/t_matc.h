/* t_matc.h */

#ifndef _TMATC_H
#define _TMATC_H
#include "tda_ctx_fwd.h"

/*  functions in t_matc.c */

char *n_mexpr(TDAContext *ctx, char *p,int n,int idx,int *rmax,int *cmax,int *ivflg,char *sel);
void mx_free(TDAContext *ctx);
double get_mxval(TDAContext *ctx, int idx,int row,int col);
double get_matval(TDAContext *ctx, int idx,int row,int col);
void m_copy(TDAContext *ctx, double *x,double *y,int n);
char *m_check1(TDAContext *ctx, char *p);                  
char *m_check2(TDAContext *ctx, char *p,int opt);                  
char *m_getmat(TDAContext *ctx, char *p,int row,int col,int *idx,char *cmd,int opt);

int m_cent(TDAContext *ctx, char *cmd,int opt);
int m_dcent(TDAContext *ctx, char *cmd);
int m_cross(TDAContext *ctx, char *cmd);
int m_chol(TDAContext *ctx, char *cmd);
int m_agg(TDAContext *ctx, char *cmd); 
int m_vec(TDAContext *ctx, char *cmd,int opt);
int m_ivec(TDAContext *ctx, char *cmd);
int m_mldes(TDAContext *ctx, char *cmd);
int m_sum(TDAContext *ctx, char *cmd,int opt);
int m_drow(TDAContext *ctx, char *cmd,int opt);
int m_diag(TDAContext *ctx, char *cmd,int opt);
int m_transp(TDAContext *ctx, char *cmd);
int m_sqrt(TDAContext *ctx, char *cmd,int opt);
int m_nrow(TDAContext *ctx, char *cmd,int opt);
int m_num(TDAContext *ctx, char *cmd);
int m_mch(TDAContext *ctx, char *cmd);
int m_mpfit(TDAContext *ctx, char *cmd);
int m_mpinv(TDAContext *ctx, char *cmd);
int m_mperm(TDAContext *ctx, char *cmd,int opt);
int m_mqap(TDAContext *ctx, char *cmd);
int m_mcel(TDAContext *ctx, char *cmd);
int m_mnc(TDAContext *ctx, char *cmd);
int m_mpz(TDAContext *ctx, char *cmd);
int m_mpb(TDAContext *ctx, char *cmd,int opt);
int m_mevs(TDAContext *ctx, char *cmd);
int m_mev(TDAContext *ctx, char *cmd);
int m_ginv(TDAContext *ctx, char *cmd);
int m_svd(TDAContext *ctx, char *cmd,int opt);
int m_wvec(TDAContext *ctx, char *cmd);
int m_wvec1(TDAContext *ctx, char *cmd);
int m_scal1(TDAContext *ctx, char *cmd);
int m_mdefc(TDAContext *ctx, char *cmd);
int m_mdefi(TDAContext *ctx, char *cmd);
int m_invs(TDAContext *ctx, char *cmd);
int m_invd(TDAContext *ctx, char *cmd);
int m_ple(TDAContext *ctx, char *cmd);
int m_nvar(TDAContext *ctx, char *cmd);
int m_print(TDAContext *ctx, char *cmd); 
int m_mkp(TDAContext *ctx, char *cmd);
int m_mls(TDAContext *ctx, char *cmd,int opt);
int m_mlp(TDAContext *ctx, char *cmd,int opt);
int m_mlpi(TDAContext *ctx, char *cmd);
int m_mul(TDAContext *ctx, char *cmd);
int m_cat(TDAContext *ctx, char *cmd,int opt);
int m_srow(TDAContext *ctx, char *cmd,int opt);
int m_sort(TDAContext *ctx, char *cmd,int opt);
int m_expr(TDAContext *ctx, char *cmd);
int m_expr1(TDAContext *ctx, char *cmd);
int m_setv(TDAContext *ctx, char *cmd); 
int m_brr(TDAContext *ctx, char *cmd); 
int m_mtrim(TDAContext *ctx, char *cmd); 
int m_mmp(TDAContext *ctx, char *cmd,int opt);
int m_mpit(TDAContext *ctx, char *cmd,int opt);
int m_mqp(TDAContext *ctx, char *cmd,int opt);
int m_mlsei1(TDAContext *ctx, char *cmd);
int m_kmet(TDAContext *ctx, char *cmd);


#endif /* _TMATC_H */

/*  end of t_matc.h */

