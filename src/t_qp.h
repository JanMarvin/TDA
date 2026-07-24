/* t_qp.h */

#ifndef _TQP_H
#define _TQP_H

/*  functions in t_qp.c */

int qld(int n,int m,int meq,int lql,int maxit,double *a,double *b,double *grad,
    double *g,double *x,double *xl,double *xu,int *iact,int *nact,
    double vsmall,double *diag);

#endif /* _TQP_H */

/*  end of t_lqp.h */


