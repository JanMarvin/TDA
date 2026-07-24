/* t_mat.h */

#ifndef _TMAT_H
#define _TMAT_H
#include "tda_ctx_fwd.h"

/*  functions in t_mat.c */

int t_mat(TDAContext *ctx);
void mdefcpy(TDAContext *ctx, char *d,char *s);
int check_local(TDAContext *ctx, char *name);
void m_cmdmsg(TDAContext *ctx);
int alloc_local(TDAContext *ctx, char *cmd);
void free_local(TDAContext *ctx);
int alloc_mat(TDAContext *ctx);
void mat_free(TDAContext *ctx);
int mat_getidx(TDAContext *ctx, char *mname,int opt);   
char *get_mname(TDAContext *ctx, char *p,char *mname,int opt);
int mat_newmat(TDAContext *ctx, char *mname,int m,int n);
int mat_ncheck(TDAContext *ctx, char *p,int opt);
int mat_alloc(TDAContext *ctx, int idx,int opt);
int mat_newidx(TDAContext *ctx, char *mname,int m,int n);
void mat_err(TDAContext *ctx, int opt);   
int mat_ncopy(TDAContext *ctx, int row,int col,double *x,char *name);
int mat_info(TDAContext *ctx);     



#endif /* _TMAT_H */

/*  end of t_mat.h */
int m_exp(TDAContext *ctx, char *cmd);
