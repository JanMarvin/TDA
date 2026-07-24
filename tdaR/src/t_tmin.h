/* t_tmin.h */

#ifndef _TTMIN_H
#define _TTMIN_H
#include "tda_ctx_fwd.h"

/*  functions in t_tmin.c */

int ts_min(TDAContext *ctx, int n,double *x,double *typx,double *xpls,double *g,double *gpls,
    double *s,double *d,double *dn,double *e,double *wk1,double *wk2,
    double *h,int *pivot);

                    /* 1 if tensor model                                    */

                    /* 1 if also analytical gradient                        */
                    /* 2 if also analytical hessian                         */



#endif /* _TTMIN_H */

/*  end of t_tmin.h */











