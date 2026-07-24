/* t_dmu.h */

#ifndef _TDMU_H
#define _TDMU_H
#include "tda_ctx_fwd.h"

/*  functions in t_dmu.c */

int dm_ccnt(TDAContext *ctx);
int dm_lcnt(TDAContext *ctx);
int dm_dump(TDAContext *ctx);
int dm_dsplit(TDAContext *ctx);
int esort(TDAContext *ctx);
int eskip(TDAContext *ctx);
int eselect(TDAContext *ctx);
int emerge(TDAContext *ctx);
int ejoin(TDAContext *ctx);

#endif /* _TDMU_H */

/*  end of t_dmu.h */

