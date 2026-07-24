/* t_tmin.h */

#ifndef _TTMIN_H
#define _TTMIN_H

/*  functions in t_tmin.c */

int ts_min(int n,double *x,double *typx,double *xpls,double *g,double *gpls,
    double *s,double *d,double *dn,double *e,double *wk1,double *wk2,
    double *h,int *pivot);

extern int TSMETH;  /* 0 if quadratic (Newton) model                        */
                    /* 1 if tensor model                                    */

extern int TSDERIV; /* 0 if only function values                            */
                    /* 1 if also analytical gradient                        */
                    /* 2 if also analytical hessian                         */

extern int TSFNC;   /* 0 use standard models, 1 user-defined functions      */
extern int TSMOD;   /* model number for TSFNC = 0.                          */

extern double *TSX;        
extern int TSXA;
extern double *TSTYPX;        
extern int TSTYPXA;
extern double *TSXPLS;        
extern int TSXPLSA;
extern double *TSG;           
extern int TSGA;
extern double *TSGPLS;        
extern int TSGPLSA;
extern double *TSS;           
extern int TSSA;
extern double *TSD;           
extern int TSDA;
extern double *TSDN;           
extern int TSDNA;
extern double *TSE;           
extern int TSEA;
extern double *TSWK1;         
extern int TSWK1A;
extern double *TSWK2;         
extern int TSWK2A;
extern double *TSH;           
extern int TSHA;
extern int *TSPIVOT;       
extern int TSPIVOTA;

#endif /* _TTMIN_H */

/*  end of t_tmin.h */











