/* t_edat.h */

#ifndef _TEDAT_H
#define _TEDAT_H
#include "tda_ctx_fwd.h"

/*  functions in t_edat.c */

int check_edj(TDAContext *ctx, int j);
int edef(TDAContext *ctx);
int def_edat(TDAContext *ctx, int mode);
void edat_off(TDAContext *ctx, int opt);
void prn_edef(TDAContext *ctx);
void edef_info(TDAContext *ctx);
int check_edat(TDAContext *ctx);
void prn_edat(TDAContext *ctx);
int get_edat(TDAContext *ctx, int i,int *sn,int *org,int *des,double *ts,double *tf);
int sort_ed(TDAContext *ctx, int opt); 
int get_spell(TDAContext *ctx, int init,int *icase,int *sn,int *org,int *des,
    double *ts,double *tf,int *spl,int *nspl);
int epdat(TDAContext *ctx);
int epsdat(TDAContext *ctx);

/*--------------------------------------------------------------------------*/







/*--------------------------------------------------------------------------*/

#endif /* _TEDAT_H */

/* end of t_edat.h */


