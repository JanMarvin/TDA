/* t_int.h */

#ifndef _TINT_H
#define _TINT_H
#include "tda_ctx_fwd.h"

/*  functions in t_int.c */

int niset(TDAContext *ctx);
void ni_info(TDAContext *ctx);
int t_int(TDAContext *ctx);
int prn_int(TDAContext *ctx, int m,double rerr);
double ni_gen(TDAContext *ctx, double a,double b,int ftyp,int nhp,int htyp,int *err);
int qsniff(TDAContext *ctx, double a,double b,double rerr,double *val,int *nf,int *nfu,int ftyp);


#endif /* _TINT_H */

/* end of t_int.h */


