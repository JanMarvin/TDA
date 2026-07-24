/* t_spss.h */

#ifndef _TSPSS_H
#define _TSPSS_H
#include "tda_ctx_fwd.h"

/*  functions in t_spss.c */

int rd_spss(TDAContext *ctx); 
int wr_spss(TDAContext *ctx);
int rd_spss1(TDAContext *ctx); 
int wr_spss1(TDAContext *ctx);
int get_dvarp(TDAContext *ctx, int *np,char *vname);
void make_arcd(TDAContext *ctx, int fn,char *fname,int len,int nrec,int vn);

#endif /* _TSPSS_H */

/* end of t_spss.h */


