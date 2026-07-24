/* t_map.h */

#ifndef _TMAP_H
#define _TMAP_H
#include "tda_ctx_fwd.h"

/*  functions in t_map.c */

int psetupg(TDAContext *ctx); 
int sdpgrat(TDAContext *ctx);
int sdpgeo(TDAContext *ctx); 
int sdpmap(TDAContext *ctx); 
double geod_distance(TDAContext *ctx, double lon1,double lat1,double lon2,double lat2);

#endif /* _TMAP_H */

/* end of t_map.h */


