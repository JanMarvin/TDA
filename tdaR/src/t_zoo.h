/* t_zoo.h */

#ifndef _TZOO_H
#define _TZOO_H
#include "tda_ctx_fwd.h"

/*  functions in t_zoo.c */

int get_zoo(TDAContext *ctx);
int dbf_init(TDAContext *ctx, int opt);
int lzs(TDAContext *ctx, int fn);
int lzd(TDAContext *ctx, int fn);
int lzh_decode(TDAContext *ctx, int fn);

#endif /* _TZOO_H */

/*  end of t_zoo.h */




