/* t_var.h */

#ifndef _TVAR_H
#define _TVAR_H
#include "tda_ctx_fwd.h"

/*  functions in t_var.c */

char get_vnchar(TDAContext *ctx, char c);
int check_vname(TDAContext *ctx, char *p);
int get_vnlen(TDAContext *ctx, char *p);
int get_vnl(TDAContext *ctx, int n,short *vinum);
int alloc_vmax(TDAContext *ctx);
int get_nidx(TDAContext *ctx);
void clear_vidx(TDAContext *ctx, int i);
int get_vidx(TDAContext *ctx, char *p,char *name);
int get_vidx1(TDAContext *ctx, char *p,char *name);
int get_nvidx(TDAContext *ctx, char *p,short *vinum);
int save_var(TDAContext *ctx, char *vd,int opt);
void var_err(TDAContext *ctx, char *vd);
int alloc_vdat(TDAContext *ctx, int idx,int opt);
int alloc_vsdat(TDAContext *ctx, int n,short *ivar);
int get_slen(TDAContext *ctx, int j,int noc);
void prn_var(TDAContext *ctx, int idx);
void prn_vname(TDAContext *ctx, int i);
void prn_vlabel(TDAContext *ctx, int i);
void prn_hvar(TDAContext *ctx);
void prn_hlabel(TDAContext *ctx);
void clear_dm(TDAContext *ctx);
void clear_var(TDAContext *ctx, int idx);
void clear_avar(TDAContext *ctx, int idx);
void free_var(TDAContext *ctx, int i);
int nlist(TDAContext *ctx);
int nl_check(TDAContext *ctx, char *nl);
void nl_free(TDAContext *ctx, int j);
int alloc_vl(TDAContext *ctx, int n);
char *get_nvi(TDAContext *ctx, char *s,int *n,int nv,short *vidx,int *nb);
char *get_nvia(TDAContext *ctx, char *s,int *n,int opt,int *nb);
int recode(TDAContext *ctx);
int ndvar(TDAContext *ctx);
int get_mxvlen(TDAContext *ctx, int n,short *vidx);
void prn_vlist(TDAContext *ctx, int n,short *vidx);
int svb_alloc(TDAContext *ctx, int n);








                        /* flag if at least one var label present.          */



extern char NLName[][VNLMax + 1];   /* names                                */



#endif /* _TVAR_H */

/* end of t_var.h */


