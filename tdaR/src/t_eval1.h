/* t_eval1.h */

#ifndef _TEVAL1_H
#define _TEVAL1_H
#include "tda_ctx_fwd.h"

/*  functions in t_eval1.c                                                  */

int v_eval1(TDAContext *ctx, int n,short cnt,int *typ,double *val,int *idx,double *res,
    int vflag,int deriv,int ip,int ivp,int iflag);

#endif /* _TEVAL1_H */

/*  end of t_eval1.h */




