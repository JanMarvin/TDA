/* t_rate.h */

#ifndef _TRATE_H
#define _TRATE_H

/*  functions in t_rate.c */

int rate(void);
void prn_resid1(int icase, int org, int des, double ts, double tf,
    double rate, double surv, double resid, double wt);

/* ------------------------------------------------------------------------ */
extern int MTerm;       /* max number of model terms                        */
extern int RMTyp;       /* type number of rate model                        */
extern int MixTyp;      /* set for mixture models                           */
extern int DisFlg;      /* set for discrete time models                     */
extern int VDAlloc;     /* 1 if model specification allocated               */
extern int XMFlg[];     /* flags for model terms (a,b,c,d)                  */
extern int NDef;        /* number of vector definitions                     */
extern short *VTerm;    /* model specification: term                        */
extern short *VOrg;     /* model specification: origin state                */
extern short *VDes;     /* model specification: destination state           */
extern short *VSNum;    /* model specification: spell number                */
extern short *VNVar;    /* model specification: number of variables         */
extern short **VNIVar;  /* pointer to variable numbers                      */

extern short *PIdx;     /* array with variables numbers of parameters       */
extern char *PIdxM;     /* array with corresponding model terms             */
extern short *PIdxPtr;  /* pointer to transition-specific parameters        */
extern short *TranTVar; /* set for VTyp 5 variables                         */
extern short *VTNum;    /* max number of covariates in transition           */
extern int VTMax;       /* max( VTNum[i] ).                                 */

extern double *POLWrk1; /* work array                                       */
extern double *POLWrk2; /* work array                                       */
extern double *BCBeta;  /* work array                                       */
extern double *BCTG;    /* work array                                       */
extern double *BCGam;   /* work array                                       */
extern double *DRJKW;   /* work array                                       */
extern double *DAJKW;   /* work array                                       */
extern double *DBJKW;   /* work array                                       */


#endif /* _TRATE_H */

/* end of t_rate.h */


