/* t_pgen.h */

#ifndef _TPGEN_H
#define _TPGEN_H

/*  functions in t_pgen.c */

void p_err(int typ,int n);
void ps_err(int typ,char *txt,int n);
void p_warn(int typ,int n);
void p_serr(char *s,int n);
void p_wmsg(int typ,char *fname,int n);
void p_vterr(int typ);
int check_cmd(int opt);  
void p_cmd(char *pcmd,int n);            
void prn_sve(void);
int eval_sve(int i);
int getd1(double *x,int ix,int ig,int opt);
int getd2(double *x,double *y,int ix,int iy,int opt);
int getdxy(double *x,double *y,int nx,int nif,int opt);
int check_nvar(int opt);
void prn_nwvar(int opt);
void prn_data(int m,int nc,int n,double *x,FILE *fd,char *fmt);
void prn1_data(double *y,double *x,int m,int ny,int nx,FILE *fd,char *fmt);
void prn1_coeff(int n,double *x,double *se,int ni,int df,short *vidx,int seflg);
void makefmt(int *fmt1,int *fmt2,char *fmts,int mlen,char sepc,int opt);
void makenfmt(int *fmt,char *fmts,int mlen,char sepc);
void pmfmt(int n,int m);
void pmtfmt(int n,int m);
void prn_cwt(void);        
void prn_sfmt(char *s,int n,char *fmt,double x);
void prn_f1mat(FILE *fd,int prn,char *fmt,int n,double *x);

#endif /* _TPGEN_H */

/* end of t_pgen.h */


