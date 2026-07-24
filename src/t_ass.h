/* t_ass.h */

#ifndef _TASS_H
#define _TASS_H

/*  functions in t_ass.c */

int gap(void);
int gqap(void);
int gqapd(int n,int niter,double alpha,double beta,int look4,int *seed,
    int *f,int *d,int *a,int *b,int *srtf,int *srtif,int *srtd,int *srtid,
    int *srtc,int *srtic,int *indexd,int *indexf,int *cost,int *fdind,
    int *opta,int *bestv,int *iter,int *ito);
int gloc(void);
int qap_w(int n,double *c,double *f,double *d,int *loc3n,double *wtmp);

#endif /* _TASS_H */

/* end of t_ass.h */


