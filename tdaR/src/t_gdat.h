/* t_gdat.h */

#ifndef _TGDAT_H
#define _TGDAT_H
#include "tda_ctx_fwd.h"

/*  functions in t_gdat.c */

int new_var(TDAContext *ctx);
double get_data(TDAContext *ctx, int j,int i);
void put_data(TDAContext *ctx, double x,int j,int i);
void put_str(TDAContext *ctx, char *buf,int blen,int j,int i,int bflag);
void get_str(TDAContext *ctx, char *buf,int j,int i);
double dscan(TDAContext *ctx, char *p,int len,int *mval);
char *skip_sep(TDAContext *ctx, char *p);
char *skip_dval(TDAContext *ctx, char *p);
void make_vfmt(TDAContext *ctx, int nidx);
int sdnvar(TDAContext *ctx);
int check_nvdef(void);
void sdnvar_close(TDAContext *ctx);
int get_nsdxy(TDAContext *ctx, char *buf,double *x,double *y);
int sdnvar_alloc(TDAContext *ctx, int opt,int n);
int check_sd(TDAContext *ctx, int opt);

/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
                        /* depends on the actual tsel command.              */






#endif /* _TGDAT_H */

/*  end of t_gdat.h */











