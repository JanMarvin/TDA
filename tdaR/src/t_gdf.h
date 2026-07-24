/* t_gdf.h */

#ifndef _TGDF_H
#define _TGDF_H
#include "tda_ctx_fwd.h"

/*  functions in t_gdf.c */

int gdf(TDAContext *ctx);
int gdf_dcheck(TDAContext *ctx);
int alloc_gdfytyp(TDAContext *ctx, int n);
int alloc_gdfgen(TDAContext *ctx, int nu,int ndim);
int gdf_gcheck(TDAContext *ctx);
void gdf_pmtyp(TDAContext *ctx);
int lsreg1(TDAContext *ctx);

#endif /* _TGDF_H */

/* end of t_gdf.h */


