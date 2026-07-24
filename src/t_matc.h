/* t_matc.h */

#ifndef _TMATC_H
#define _TMATC_H

/*  functions in t_matc.c */

char *n_mexpr(char *p,int n,int idx,int *rmax,int *cmax,int *ivflg,char *sel);
void mx_free(void);
double get_mxval(int idx,int row,int col);
double get_matval(int idx,int row,int col);
void m_copy(double *x,double *y,int n);
char *m_check1(char *p);                  
char *m_check2(char *p,int opt);                  
char *m_getmat(char *p,int row,int col,int *idx,char *cmd,int opt);

int m_cent(char *cmd,int opt);
int m_dcent(char *cmd);
int m_cross(char *cmd);
int m_chol(char *cmd);
int m_agg(char *cmd); 
int m_vec(char *cmd,int opt);
int m_ivec(char *cmd);
int m_mldes(char *cmd);
int m_sum(char *cmd,int opt);
int m_drow(char *cmd,int opt);
int m_diag(char *cmd,int opt);
int m_transp(char *cmd);
int m_sqrt(char *cmd,int opt);
int m_nrow(char *cmd,int opt);
int m_num(char *cmd);
int m_mch(char *cmd);
int m_mpfit(char *cmd);
int m_mpinv(char *cmd);
int m_mperm(char *cmd,int opt);
int m_mqap(char *cmd);
int m_mcel(char *cmd);
int m_mnc(char *cmd);
int m_mpz(char *cmd);
int m_mpb(char *cmd,int opt);
int m_mevs(char *cmd);
int m_mev(char *cmd);
int m_ginv(char *cmd);
int m_svd(char *cmd,int opt);
int m_wvec(char *cmd);
int m_wvec1(char *cmd);
int m_scal1(char *cmd);
int m_mdefc(char *cmd);
int m_mdefi(char *cmd);
int m_invs(char *cmd);
int m_invd(char *cmd);
int m_ple(char *cmd);
int m_nvar(char *cmd);
int m_print(char *cmd); 
int m_mkp(char *cmd);
int m_mls(char *cmd,int opt);
int m_mlp(char *cmd,int opt);
int m_mlpi(char *cmd);
int m_mul(char *cmd);
int m_cat(char *cmd,int opt);
int m_srow(char *cmd,int opt);
int m_sort(char *cmd,int opt);
int m_expr(char *cmd);
int m_expr1(char *cmd);
int m_setv(char *cmd); 
int m_brr(char *cmd); 
int m_mtrim(char *cmd); 
int m_mmp(char *cmd,int opt);
int m_mpit(char *cmd,int opt);
int m_mqp(char *cmd,int opt);
int m_mlsei1(char *cmd);
int m_kmet(char *cmd);

extern double *MX[];
extern double *MX1[];
extern int MXRow[];
extern int MXCol[];
extern int MXIV[];
extern char MXName[]; 

#endif /* _TMATC_H */

/*  end of t_matc.h */

