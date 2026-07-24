/* t_sd.h */

#ifndef _TSD_H
#define _TSD_H
#include "tda_ctx_fwd.h"

/*  functions in t_sd.c */

int sd_getdata(TDAContext *ctx, int i,int order,int cflag,int opt);
int sdplot(TDAContext *ctx); 
int sdplot31(TDAContext *ctx);
int sdplot32(TDAContext *ctx);
int sdplot33(TDAContext *ctx);
int sdpdata(TDAContext *ctx);
int sdinf(TDAContext *ctx);
int sdencl(TDAContext *ctx);
int sdsel(TDAContext *ctx);

#endif /* _TSD_H */

/* end of t_sd.h */


