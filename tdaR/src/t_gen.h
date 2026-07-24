/* t_gen.h */

#ifndef _TGEN_H
#define _TGEN_H
#include "tda_ctx_fwd.h"
#include "tda_attr.h"

/*  functions in t_gen.c */

void gerr_exit(TDAContext *ctx, int n);
void memrq(TDAContext *ctx, int n, int m);
void prn_mem(TDAContext *ctx);      
void newline(TDAContext *ctx);           
char *check_comment(TDAContext *ctx, char *p);
int check_drec(TDAContext *ctx, char *buf);
void prn_message(TDAContext *ctx, int rec,int opt,int wflag);
int get_sline(TDAContext *ctx, char *p); 
void prnchar(TDAContext *ctx, unsigned char c, int n, int mode);
void fprnchar(TDAContext *ctx, FILE *fd, unsigned char c, int n, int mode);
int mpr(TDAContext *ctx, int m,int n,double *a,int ivflg,double *b,char *fmt,char *fn,char *s,int aflag);
int mprf(TDAContext *ctx, int m,int n,float *a,char *fmt,FILE *fd);
int mprd(TDAContext *ctx, int m,int n,double *a,char *fmt,FILE *fd);
char *skip_b(TDAContext *ctx, char *p);
char *skip_c(TDAContext *ctx, char *p);
char *skip_cb(TDAContext *ctx, char *p);
char *skip_int(TDAContext *ctx, char *p);
char *skip_dbl(TDAContext *ctx, char *p);
char *skip_expr(TDAContext *ctx, char *p);
char *skip_blev(TDAContext *ctx, char *p);
char *skip_nc(TDAContext *ctx, char *p);
char *skip_com(TDAContext *ctx, char *p);
char *skip_xa(TDAContext *ctx, char *p);
int get_fname(TDAContext *ctx, char *p,char *q);
int get_float(char *p,int nmax,float *x);
void printf1(TDAContext *ctx, const char *fmt, ...) TDA_PRINTF(2, 3);
/*  x_show: the xshow command.  The X11 front end it needs is not part
    of this tree, so t_xstub.c provides the only definition.  */
int x_show(TDAContext *ctx);
void printf2(TDAContext *ctx, const char *fmt, ...) TDA_PRINTF(2, 3);
void printfe(TDAContext *ctx, const char *fmt, ...) TDA_PRINTF(2, 3);
void fflushe(TDAContext *ctx);
#if TIME_ON
void prn_time(TDAContext *ctx, FILE *fd);      
void prn1_time(TDAContext *ctx, char *buf); 
#endif
  
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/

#endif










                                                

