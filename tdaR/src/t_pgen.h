/* t_pgen.h */

#ifndef _TPGEN_H
#define _TPGEN_H
#include "tda_ctx_fwd.h"

/*  functions in t_pgen.c */

void p_err(TDAContext *ctx, int typ,int n);
void ps_err(TDAContext *ctx, int typ,char *txt,int n);
void p_warn(TDAContext *ctx, int typ,int n);
void p_serr(TDAContext *ctx, char *s,int n);
void p_wmsg(TDAContext *ctx, int typ,char *fname,int n);
void p_vterr(TDAContext *ctx, int typ);
int check_cmd(TDAContext *ctx, int opt);  
void p_cmd(TDAContext *ctx, char *pcmd,int n);            
void prn_sve(TDAContext *ctx);
int eval_sve(TDAContext *ctx, int i);
int getd1(TDAContext *ctx, double *x,int ix,int ig,int opt);
int getd2(TDAContext *ctx, double *x,double *y,int ix,int iy,int opt);
int getdxy(TDAContext *ctx, double *x,double *y,int nx,int nif,int opt);
int check_nvar(TDAContext *ctx, int opt);
void prn_nwvar(TDAContext *ctx, int opt);
void prn_data(TDAContext *ctx, int m,int nc,int n,double *x,FILE *fd,char *fmt);
#ifdef TDA_R_PACKAGE
void export_prn(TDAContext *ctx, const char *name, int m,int nc,int n,double *x);
#endif
void prn1_data(TDAContext *ctx, double *y,double *x,int m,int ny,int nx,FILE *fd,char *fmt);
void prn1_coeff(TDAContext *ctx, int n,double *x,double *se,int ni,int df,short *vidx,int seflg);
void makefmt(TDAContext *ctx, int *fmt1,int *fmt2,char *fmts,size_t fmtsz,int mlen,char sepc,int opt);
void makenfmt(TDAContext *ctx, int *fmt,char *fmts,size_t fmtsz,int mlen,char sepc);
void pmfmt(TDAContext *ctx, int n,int m);
void pmtfmt(TDAContext *ctx, int n,int m);
void prn_cwt(TDAContext *ctx);        
void prn_sfmt(TDAContext *ctx, char *s,int n,char *fmt,double x);
void prn_f1mat(TDAContext *ctx, FILE *fd,int prn,char *fmt,int n,double *x);

#endif /* _TPGEN_H */

/* end of t_pgen.h */


