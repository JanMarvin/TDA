/* t_mata.h */

#ifndef _TMATA_H
#define _TMATA_H
#include "tda_ctx_fwd.h"

/*  functions in t_mata.c */

int m_rperm(TDAContext *ctx, int n,int *icn,int licn,int *ip,int *lenr,int *iperm,
    int *pr,int *arp,int *cv,int *out);
int m_perm(TDAContext *ctx, int n,int *icn,int licn,int *ip,int *lenr,int *arp,   
    int *ib,int *lowl,int *numb,int *prev);
int mpfit(TDAContext *ctx, int n,int m,double *a,double *u,double *v,double *b,int *iter,double *eps);

#endif /* _TMATA_H */

/*  end of t_mata.h */















