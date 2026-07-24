/* t_stata.h */

#ifndef _TSTATA_H
#define _TSTATA_H
#include "tda_ctx_fwd.h"

/*  functions in t_stata.c */

int rd_stata(TDAContext *ctx); 
int wr_stata(TDAContext *ctx);
int st_gets(TDAContext *ctx, char *p);
int st_geti(TDAContext *ctx, char *p);
double st_getsf(TDAContext *ctx, char *p,int l); 
double st_getf(TDAContext *ctx, char *p);
double st_getd(TDAContext *ctx, char *p);


#endif /* _TSTATA_H */

/* end of t_stata.h */


