/* t_con.h */

#ifndef _TCON_H
#define _TCON_H
#include "tda_ctx_fwd.h"

/*  functions in t_con.c */

int p_con(TDAContext *ctx, char *pcmd,int nx,int nif,int mw,int nw,double *w,int *ne,int *ni);
int con_proc(TDAContext *ctx, char *pcmd);
void con_free(TDAContext *ctx);


#endif /* _TCON_H */

/* end of t_con.h */


