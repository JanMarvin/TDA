/* t_ml.h */

#ifndef _TML_H
#define _TML_H
#include "tda_ctx_fwd.h"

/*  functions in t_ml.c */

void set_mldef(TDAContext *ctx);
void set_mlopt(TDAContext *ctx);
int get_dsv(TDAContext *ctx, int n,double *par,int gmina,double *lb,double *ub,int nl);
int ml_init(TDAContext *ctx, int opt,int typ,int gmina);
int syminv1(TDAContext *ctx, int n, double *d, double *h);
int syminv2(TDAContext *ctx, int n, double *x);
void fn(TDAContext *ctx, double *x, int mod, int typ, int scal, int opt,int *err);
int checkov(TDAContext *ctx);
void prn_mlres(TDAContext *ctx, int typ);
void prn_ml1res(TDAContext *ctx, int typ);
void prvec(TDAContext *ctx, char *txt,int n,double *x);
void prmat(TDAContext *ctx, char *txt,int n,int m,double *x);
void prhess(TDAContext *ctx, char *txt,int n);
void prval(TDAContext *ctx, char *txt,double x);
void prot_init(TDAContext *ctx, int typ,int mod,int gmina);

/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/





                            /* hessian requested                            */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/





/* ------------------------------------------------------------------------ */

#endif /* _TML_H */

/*  end of t_ml.h */
