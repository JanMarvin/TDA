/* t_qrmod.h */

#ifndef _TQRMOD_H
#define _TQRMOD_H

/*  functions in t_qrmod.c */

int qreg(void);
int rmod(void);

/*--------------------------------------------------------------------------*/
extern int QRTyp;           /* type of quantal response model               */
extern int NWave;           /* numbe of waves                               */
extern int IntFlg;          /* zero if without an intercept                 */
extern int NPX;             /* number of independent variables              */
extern int NPX1;            /* NPX + IntFlg                                 */

extern short *PYVar;        /* index of dependent variable in wave ...      */
extern short **PXVar;       /* indeces of independent variables in wave ... */

extern int NPZ;             /* number of Z variables                        */
extern short **PZVar;       /* indices of z variables in wave ...           */

extern int NPU;             /* number of U variables                        */

extern int NYCat;           /* number of categories in dependent variables  */
extern int *YCat;           /* categories of dep. variables                 */
extern int *YCatI;          /* inverse of YCat                              */
extern double *PYFreq;      /* weighted frequencies of categories           */
/* ------------------------------------------------------------------------ */
extern int PDMCFlg;         /* if missing value code defined                */
extern char *PDMCIdx;       /* vector with select flags                     */
extern short *PDMCTi;       /* number of waves of participation             */
/* ------------------------------------------------------------------------ */
extern int QR_INIT;         /* set to 1 for first call of fn...().          */
extern double *QRProb;      /* temp. storage for probabilities              */
extern double *QRTmp;       /* temp. storage                                */

extern int *L4_ISET;        /* temp storage used in fn_logit4               */
extern int *L4_IPTR;           
extern double *L4_XB;          
/* ------------------------------------------------------------------------ */
extern int RMPN;            /* number of patterns in Rasch model            */
extern int RMPV;            /* number of variables                          */
extern short *RMPAT;        /* array with patterns                          */
extern int *RMPF;           /* frequencies of patterns                      */

#endif /* _TQRMOD_H */

/*  end of t_qrmod.h */

