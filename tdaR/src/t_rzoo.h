/* t_rzoo.h */

#ifndef _TRZOO_H
#define _TRZOO_H
#include "tda_ctx_fwd.h"

/*  functions in t_rzoo.c */

void arcd_off(TDAContext *ctx);
int arcd(TDAContext *ctx);
void prn_afiles(TDAContext *ctx);
int get_drec(TDAContext *ctx, int fn, char *buf, int nmax);
int get_record(TDAContext *ctx, int fn, char *buf, int nmax);
int alloc_avar(TDAContext *ctx, int idx,int n,int opt);
int check_avar(TDAContext *ctx, int n,int adic);
int get_avar(TDAContext *ctx);
void get_astr(TDAContext *ctx, char *buf,int i);
int arcc(TDAContext *ctx);
int arcv(TDAContext *ctx);
int arcvc(TDAContext *ctx);
int get_afmt(TDAContext *ctx, int *w1,int *w2);

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/*  For each file we allocate a read buffer where the decompression is      */
/*  done, and a line buffer.                                                */

/* ------------------------------------------------------------------------ */


#endif /* _TRZOO_H */

/*  end of t_rzoo.h */












