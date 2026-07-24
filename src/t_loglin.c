/****************************************************************************/
/*  t_loglin                                                                */
/*                                                                          */
/*  TDA. Program for Transition Data Analysis, written by Goetz Rohwer.     */
/*  Copyright (C) 1989,1991-97 Goetz Rohwer. All rights reserved.           */
/*                                                                          */
/*  This file is part of TDA.                                               */
/*                                                                          */
/*  TDA is free software; you can redistribute it and/or modify             */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 2 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  TDA is distributed in the hope that it will be useful,                  */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License       */
/*  along with this program (it should be in a file named COPYING),         */
/*  if not, write to the Free Software Foundation, Inc.,                    */
/*  675 Mass Ave, Cambridge, MA 02139, USA.                                 */
/*                                                                          */

#include "tda.h"
#include "t_gen.h"
#include "t_pgen.h"
#include "t_parm.h"
#include "t_ml.h"
#include "t_gdat.h"
#include "t_alloc.h"
#include "t_var.h"
#include "t_gf.h"
#include "t_cdf.h"
#include "t_freq.h"
#include "t_con.h"


/*  functions in t_loglin.c */

int loglin(void);
int ll_mcheck(void);
int ll_dcol(void);
int ll_alloc(int opt);
void ll_free(void);
int ll_tab(void);
int loglin1(void);
void ll_ppar(void);
void ll_error(int err,char *mod);
void prn_mdes(FILE *fd, char *desc, int n);
void prn_lldmat(int i);
void prn_llcov(int i);
void prn_ctab(int mode);
int llmodcheck(char *desc);
int expand1(char *s, char *r, int mode);
int expand2(char *s, char *r, int mode);
char *distribute(int n, char *r, int mode, int *len);
int build_dm(char *descr, int col, int *slen, int mode);
int design_comb(int dn, int *cat, int *ptr);
int ll_readdm(int i, int dcol, int mode);
int logl_fit(void);
void logl_solve(int n, double *lu, double *b, double *x);
void logl_inv(int n, double *tri, double *xinv);
int logl_chol(int n, double *sym, double *tri);
/* void logl_gen(double *c, double *obs, double *expv, double *res); */
int logscrn(void);
int logscr_conf(int m, int *iset, int *jset);
int logscr_comb(int *iset, int n, int m, int last);
int logscr_eval(short *iar, int nc, int nv, int ibeg);
void logscr_res(void);
double logscr_lr(void);
void logscr_fit(int ncon);
int e_expm(void);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

#define MaxDX 200           /* max number of model terms                    */
#define MaxABC 200 

int NDIM = 0;               /* number of dimensions (= PMNV)                */
int NDIM1 = 0;              /* number of dimensions minus 1                 */
int *CTNCat;                /* number of categories in i.th dimension       */
int CTNCatA = 0;            /* if allocated                                 */
int **CTCat;                /* value of categories in the i.th dimension    */
int CTCatA = 0;             /* if allocated                                 */
int *CTCatAI;               /* if CTCat[i] is allocated                     */
int CTCatAIA = 0;

int NCTab = 0;              /* actual number of cells of the table          */
int NCTab1 = 0;             /* number of cells in a complete table          */
int CTNum = 0;              /* number of table elements                     */
int CTIFlg = 0;             /* set if table is incomplete                   */

short *CTVal;               /* value combinations of the table              */
int CTValA = 0;
int *CTFreq;                /* frequencies of table rows                    */
int CTFreqA = 0;
double *CTSCAL;             /* scale values                                 */
int CTSCALA = 0;

int LLTyp[MaxLL];           /* Type of model                                */
int LLTyp1 = 0;             /* Counter for models of type 1                 */
int LLTyp2 = 0;             /* Counter for models of type 2                 */

int *LLChk;                 /* Used to expand model descriptions            */
int LLD1Len = 0;            /* max length of first expansion                */
int LLD2Len = 0;            /* max length of first expansion                */
char *LLD1;                 /* first expansion                              */
int LLD1A = 0;
char *LLD2;                 /* second expansion                             */
int LLD2A = 0;
int ExMax;                  /* max size of second expansion                 */
int *ExIdx;
int ExIdxA = 0;
int *ExVc;
int ExVcA = 0;
char *ExV;
int ExVA = 0;
char *ExVl;
int ExVlA = 0;
char *ExVTyp;
int ExVTypA = 0;
char *ExOpl;
int ExOplA = 0;
char *ExMult;
int ExMultA = 0;
char *ExUsed;
int ExUsedA = 0;

int LLDColM = 0;            /* max number of columns of the design matrix   */
int *LLDCol;                /* qrray with design matrix columns             */
int LLDColA = 0;

char *CTVLab;               /* label of model parameters                    */
int CTVLabA = 0;
int CTVLenM = 0;            /* max length of available space for labels     */
int CTVLen = 0;             /* actual length of string with labels          */

short *LLDMat;              /* design matrix                                */
int LLDMatA = 0;
int LLDMCol = 0;            /* actual number of columns of design matrix    */

double *CTRes;              /* residuals                                    */
int CTResA = 0;
double *CTMat;
int CTMatA = 0;
double *CTB;
int CTBA = 0;
double *CTPar;
int CTParA = 0;
double *CTWork;
int CTWorkA = 0;
double *CTFit;              /* fitted table                                 */
int CTFitA = 0;

int CTNPar = 0;             /* number of model parameters                   */
int CTNDF = 0;              /* degrees of freedom                           */
double CTCHI = 0.0;         /* Pearson X square statistic                   */
double CTLR = 0.0;          /* Likelihood ratio statistic                   */

/*  The following variables are used for marginal and partial association   */

double *PART;   
int PARTA = 0;
double *MARG;
int MARGA = 0;
double *GSQ;
int GSQA = 0;
int *DFS;
int DFSA = 0;
int *DGFR;
int DGFRA = 0;
int *LSISet;
int LSISetA = 0;
int *LSJSet;
int LSJSetA = 0;

float *LSx;
int LSxA = 0;
float *LSy;
int LSyA = 0;
short int *LSip;
int LSipA = 0; 
short *LSim;
int LSimA = 0;        
short *LSconf;
int LSconfA = 0;
short *LSdim;
int LSdimA = 0;
short *LScord;
int LScordA = 0;
short *LSsiz;
int  LSsizA = 0;

/* ------------------------------------------------------------------------ */
/*  loglin    log-linear models                                             */
/*                                                                          */
/*  loglin(                                                                 */
/*      maxcat=...,             max cat per dimension, def. 1000            */
/*      w = ...,                variable containing weights                 */
/*      scale = ...,            variable for scaling                        */
/*      ptab=...,               print table                                 */
/*      nfmt=...,               print format table values, def. 2           */
/*      mod=...,                model descriptions                          */
/*      mxit=...,               max number of iterations, def. 20           */
/*      tolf=...,               convergence tolerance, def. 1e-8            */
/*      tfmt=...,               print format parameter, def. 10.4           */
/*      ppar=...,               output file for parameter estimates         */
/*      pres=...,               residuals                                   */
/*      pcov=...,               covariance matrix                           */
/*      mfmt=...,               print format cov matrix, def. 12.4          */
/*      screen,                 screening for associations                  */
/*                                                                          */
/*  ) = X1,...,Xn;              table dimensions                            */
/*                                                                          */
/*  Number of variables = number of dimensions is max 26 (A,..,Z)           */
/*                                                                          */
/*  Return 0 if successful, -1 if error.                                    */

int loglin(void)
{
    int err;

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("Loglinear models. Current memory: %d bytes.\n\n",MemReq);

    MxItFlg = 0;
    TOLF = 1.e-8;

    if (parm(CmdBuf + 6,4,1))     /* get parameters */
        goto LLFin;

    if (PMNFmtF == 0) {
        PMNFmt = 2;
        makenfmt(&PMNFmt,PMNFmtS,1,SEPC);
    }
    if (MxIter < 1 || MxItFlg == 0)
        MxIter = 20;

    NDIM = PMNV;                /* number of table dimensions */
    if (NDIM > 26) {
        printf1("Error: max number of dimensions is 26.\n");
        goto LLFin;
    }
    if (ll_tab())               /* make table */
        goto LLFin;

    if (PMTabFDef) {            /* print table */ 
        prn_ctab(1);
        printf1("Table written to: %s\n",PMTabFName);
    }
    newline();
     
    if (SCRNFlg) {              /* screening */
        if (logscrn())
            goto LLFin;
    }

    /* check model requests */

    if (PMModN == 0) {
        err = 0;
        goto LLFin;
    }
    if (ll_mcheck())    /* check model descriptions */
        goto LLFin;

    if (ll_dcol())      /* get max number of design matrix columns */
        goto LLFin;

    if (ll_alloc(1))    /* allocate memory for model estimation */
        goto LLFin;

    if (loglin1())      /* estimate models */
        goto LLFin;

    err = 0;

LLFin:
    con_free();
    ll_alloc(0);
    ll_free();
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ll_mcheck()   Check model descriptions                                  */
/*                Return 0 if OK, -1 if error.                              */

int ll_mcheck(void)
{
    FILE *fd;
    register int i;
    int err,n;
    char *p;

    err = -1;

    printf1("Check of model requests.\n");

    for (i = 0; i < PMModN; ++i) {
        p = PMModS[i];
        if (*p == 'd' && *(p + 1) == ':') {
            if ((fd = fopen(p + 2,OPEN_RD))) {   
                LLTyp[i] = -1;
                LLTyp2++;
                fclose(fd);
            }
            else {
                ll_error(-9,p);
                goto LLMCFin;
            }
        }
        else {
            n = llmodcheck(p);
            if (n < 0) {
                ll_error(n,p);
                goto LLMCFin;
            }
            LLTyp[i] = n;
            LLTyp1++;
        }
    }
    LLD1Len = 1;

    for (i = 0; i < PMModN; ++i) {      /* get max length of first expansion */
        if (LLTyp[i] > 0) {
            n = expand1(PMModS[i],p,0);
            if (n < 0) {
                ll_error(n,PMModS[i]);
                goto LLMCFin;
            }
            LLTyp[i] = n;
            LLD1Len = imax(LLD1Len,n);
        }
    }   

    /* allocate memory for first expansion */

    if (!(LLD1 = (char *) calloc (LLD1Len + 2,sizeof(char)))) {
        p_err(-2,1);
        goto LLMCFin;
    }
    LLD1A = LLD1Len + 2;
    memrq(LLD1A,sizeof(char));

    /* allocate memory for second expansion */

    ExMax = LLD1Len;

    if (!(ExIdx  = (int *)  calloc(ExMax + 1,sizeof(int)))) {
        p_err(-2,1);
        goto LLMCFin;
    }
    ExIdxA = ExMax + 1;
    memrq(ExIdxA,sizeof(int));

    if (!(ExVc   = (int *)  calloc(ExMax + 1,sizeof(int)))) { 
        p_err(-2,1);
        goto LLMCFin;
    }
    ExVcA = ExMax + 1;
    memrq(ExVcA,sizeof(int));

    if (!(ExV = (char *) calloc(ExMax + 1,sizeof(char)))) {
        p_err(-2,1);
        goto LLMCFin;
    }
    ExVA = ExMax + 1;
    memrq(ExVA,sizeof(char));

    if (!(ExVl   = (char *) calloc(ExMax + 1,sizeof(char)))) {
        p_err(-2,1);
        goto LLMCFin;
    }
    ExVlA = ExMax + 1;
    memrq(ExVlA,sizeof(char));

    if (!(ExVTyp = (char *) calloc(ExMax + 1,sizeof(char)))) {
        p_err(-2,1);
        goto LLMCFin;
    }
    ExVTypA = ExMax + 1;
    memrq(ExVTypA,sizeof(char));

    if (!(ExOpl  = (char *) calloc(ExMax + 1,sizeof(char)))) {
        p_err(-2,1);
        goto LLMCFin;
    }
    ExOplA = ExMax + 1;
    memrq(ExOplA,sizeof(char));

    if (!(ExMult = (char *) calloc(ExMax + 1,sizeof(char)))) {
        p_err(-2,1);
        goto LLMCFin;
    }
    ExMultA = ExMax + 1;
    memrq(ExMultA,sizeof(char));

    if (!(ExUsed = (char *) calloc(ExMax + 1,sizeof(char)))) {
        p_err(-2,1);
        goto LLMCFin;
    }
    ExUsedA = ExMax + 1;
    memrq(ExUsedA,sizeof(char));

    LLD2Len = 1;                    /* get max length of second expansion */

    for (i = 0; i < PMModN; ++i) {
        if (LLTyp[i] > 0) {
            n = expand1(PMModS[i],LLD1,1);
            n = expand2(LLD1,p,0);
            if (n < 0) {
                ll_error(n,PMModS[i]);
                goto LLMCFin;
            }
            LLTyp[i] = n;
            LLD2Len = imax(LLD2Len,n);
        }
    }   

    /* allocate memory for second expansion */

    if (!(LLD2 = (char *) calloc (LLD2Len + 2,sizeof(char)))) {
        p_err(-2,1);
        goto LLMCFin;
    }
    LLD2A = LLD2Len + 2;
    memrq(LLD2A,sizeof(char));

    err = 0;

LLMCFin:
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ll_dcol     Get number of design matrix columns.                        */
/*              Return 0 if OK, -1 if error.                                */

int ll_dcol(void)
{
    register int i;
    int n,err,slen;

    err = -1;

    /*  allocate memory for number of design matrix columns */

    if (!(LLDCol = (int *)calloc(PMModN,sizeof(int)))) {
        p_err(-2,1);
        goto LLDCFin;
    }
    LLDColA = PMModN;        
    memrq(LLDColA,sizeof(int));

    LLDColM = 1;

    for (i = 0; i < PMModN; ++i) {     
        if (LLTyp[i] > 0) {
            n = expand1(PMModS[i],LLD1,1);
            n = expand2(LLD1,LLD2,1);
            n = build_dm(LLD2,0,&slen,0);
            if (n < 0) {
                ll_error(n,PMModS[i]);
                goto LLDCFin;
            }
            LLDCol[i] = n;
            LLDColM = imax(LLDColM,n);
        }
        else if (LLTyp[i] < 0) {
            n = ll_readdm(i,0,0);
            if (n < 0) {
                ll_error(n,PMModS[i]);
                goto LLDCFin;
            }
            LLDCol[i] = n;
            LLDColM = imax(LLDColM,n);
        }
    }   
    err = 0;

LLDCFin:
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ll_alloc(opt)   If opt != 0 allocate memory for model estimation,       */
/*                  else free previously allocated memory.                  */
/*                  Return 0 if OK, -1 if error.                            */

int ll_alloc(int opt)
{
    int err = 0;

    if (opt == 0)
        goto LLAFin;

    err = -1;

    /*  allocate memory for model parameter labels */

    CTVLenM = (LLDColM + 1) * (NDIM * 20 + 1);

    if (!(CTVLab = (char *) calloc(CTVLenM + 10,sizeof(char)))) {
        p_err(-2,1);
        goto LLAFin;
    }
    CTVLabA = CTVLenM + 10;
    memrq(CTVLabA,sizeof(char));

    /*  allocate memory for design matrix LLDMat */

    if (!(LLDMat = (short *)calloc(NCTab * LLDColM + 1,sizeof(short)))) { 
        p_err(-2,1);
        goto LLAFin;
    }
    LLDMatA = NCTab * LLDColM + 1;
    memrq(LLDMatA,sizeof(short));

    /*  memory allocation for model estimation, space is allocated for
        the largest model. */

    if (!(CTRes  = (double *) calloc(NCTab + 1,sizeof(double)))) {  
        p_err(-2,1);
        goto LLAFin;
    }
    CTResA = NCTab + 1;
    memrq(CTResA,sizeof(double));

    if (!(CTMat  = (double *) calloc(LLDColM * LLDColM + 1,sizeof(double)))) {  
        p_err(-2,1);
        goto LLAFin;
    }
    CTMatA = LLDColM + 1;
    memrq(CTMatA,sizeof(double));

    if (!(CTB = (double *) calloc(LLDColM + 1,sizeof(double)))) { 
        p_err(-2,1);
        goto LLAFin;
    }
    CTBA = LLDColM + 1;
    memrq(CTBA,sizeof(double));

    if (!(CTPar  = (double *) calloc(LLDColM + 1,sizeof(double)))) {  
        p_err(-2,1);
        goto LLAFin;
    }
    CTParA = LLDColM + 1;
    memrq(CTParA,sizeof(double));

    if (!(CTWork = (double *) calloc(LLDColM + 1,sizeof(double)))) {
        p_err(-2,1);
        goto LLAFin;
    }
    CTWorkA = LLDColM + 1;
    memrq(CTWorkA,sizeof(double));

    if (!(CTFit = (double *) calloc(NCTab + 1,sizeof(double)))) {
        p_err(-2,1);
        goto LLAFin;
    }
    CTFitA = NCTab + 1;
    memrq(CTFitA,sizeof(double));

    return(0);

LLAFin:
    if (CTVLabA > 0) {
        free((char *)CTVLab);
        memrq(-CTVLabA,sizeof(char));
        CTVLabA = 0;    
    }
    if (LLDMatA > 0) {
        free((char *)LLDMat);
        memrq(-LLDMatA,sizeof(short));
        LLDMatA = 0;    
    }
    if (CTResA > 0) {
        free((char *)CTRes);
        memrq(-CTResA,sizeof(double));
        CTResA = 0;    
    }
    if (CTMatA > 0) {
        free((char *)CTMat);
        memrq(-CTMatA,sizeof(double));
        CTMatA = 0;    
    }
    if (CTBA > 0) {
        free((char *)CTB);
        memrq(-CTBA,sizeof(double));
        CTBA = 0;    
    }
    if (CTParA > 0) {
        free((char *)CTPar);
        memrq(-CTParA,sizeof(double));
        CTParA = 0;    
    }
    if (CTWorkA > 0) {
        free((char *)CTWork);
        memrq(-CTWorkA,sizeof(double));
        CTWorkA = 0;    
    }
    if (CTFitA > 0) {
        free((char *)CTFit);
        memrq(-CTFitA,sizeof(double));
        CTFitA = 0;    
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ll_free. Free allocated memory.                                         */

void ll_free(void)
{
    register int i;

    if (CTNCatA > 0) {
        free((char *)CTNCat);
        memrq(-CTNCatA,sizeof(int));
        CTNCatA = 0;    
    }
    if (CTCatA > 0) {
        for (i = 0; i < NDIM; ++i) {
            if (CTCatAI[i] > 0) {
                free((char *)CTCat[i]);
                memrq(-CTCatAI[i],sizeof(int));
            }
        }
        free((char *)CTCat);
        memrq(-CTCatA,sizeof(int *));
        CTCatA = 0;    
    }
    if (CTCatAIA > 0) {
        free((char *)CTCatAI);
        memrq(-CTCatAIA,sizeof(int));
        CTCatAIA = 0;    
    }
    if (CTValA > 0) {
        free((char *)CTVal);
        memrq(-CTValA,sizeof(short));
        CTValA = 0;    
    }
    if (CTFreqA > 0) {
        free((char *)CTFreq);
        memrq(-CTFreqA,sizeof(int));
        CTFreqA = 0;    
    }
    if (CTSCALA > 0) {
        free((char *)CTSCAL);
        memrq(-CTSCALA,sizeof(double));
        CTSCALA = 0;    
    }
    if (LLD1A > 0) {
        free((char *)LLD1);
        memrq(-LLD1A,sizeof(char));
        LLD1A = 0;    
    }
    if (LLD2A > 0) {
        free((char *)LLD2);
        memrq(-LLD2A,sizeof(char));
        LLD2A = 0;    
    }
    if (ExIdxA > 0) {
        free((char *)ExIdx);
        memrq(-ExIdxA,sizeof(int));
        ExIdxA = 0;    
    }
    if (ExVcA > 0) {
        free((char *)ExVc);
        memrq(-ExVcA,sizeof(int));
        ExVcA = 0;    
    }
    if (ExVA > 0) {
        free((char *)ExV);
        memrq(-ExVA,sizeof(char));
        ExVA = 0;    
    }
    if (ExVlA > 0) {
        free((char *)ExVl);
        memrq(-ExVlA,sizeof(char));
        ExVlA = 0;    
    }
    if (ExVTypA > 0) {
        free((char *)ExVTyp);
        memrq(-ExVTypA,sizeof(char));
        ExVTypA = 0;    
    }
    if (ExOplA > 0) {
        free((char *)ExOpl);
        memrq(-ExOplA,sizeof(char));
        ExOplA = 0;    
    }
    if (ExMultA > 0) {
        free((char *)ExMult);
        memrq(-ExMultA,sizeof(char));
        ExMultA = 0;    
    }
    if (ExUsedA > 0) {
        free((char *)ExUsed);
        memrq(-ExUsedA,sizeof(char));
        ExUsedA = 0;    
    }
    if (LLDColA > 0) {
        free((char *)LLDCol);
        memrq(-LLDColA,sizeof(int));
        LLDColA = 0;    
    }
}

/* ------------------------------------------------------------------------ */
/*  ll_tab.  Create table.                                                  */
/*           Return 0 if OK, -1 if error.                                   */

int ll_tab(void)
{
    register int i,j,k,l;
    int err,m,n,eerr;
    double tmp;

    err = -1;

    printf1("Definition and structure of contingency table.\n");
    if (PMWVar >= 0)  
        printf1("Frequencies defined by: %s\n",VName[PMWVar]);
    if (PMSCAL >= 0)  
        printf1("Scale variable defined by: %s\n",VName[PMSCAL]);
    newline();

    /*  Calculation of categories of dimensions of the table:
        NCTab = number of cells of the complete table.
        MaxYC = maximum number of categories in a dimension.
        CTNCat[i] = number of categories in i.th dimension,
        CTCat[i][j] (j = 1,2,3,...,CTNCat[i]) is the value of the j.th
        categorie in the i.th dimension. */

    if (!(CTNCat = (int *) calloc(NDIM,sizeof(int)))) {   
        p_err(-2,1);
        goto LLTFin;
    }
    CTNCatA = NDIM;
    memrq(CTNCatA,sizeof(int));

    if (!(CTCatAI = (int *) calloc(NDIM,sizeof(int)))) {   
        p_err(-2,1);
        goto LLTFin;
    }
    CTCatAIA = NDIM;
    memrq(CTCatAIA,sizeof(int));

    if (!(CTCat = (int **) calloc(NDIM,sizeof(int *)))) {   
        p_err(-2,1);
        goto LLTFin;
    }
    CTCatA = NDIM;
    memrq(CTCatA,sizeof(int *));

    NCTab = 1;

    if (alloc_acr(NOC + 1))                 /* data */
        goto LLTFin;
    if (alloc_acs(PMMaxCat + 1))            /* values of categories */
        goto LLTFin;
    if (alloc_aci(PMMaxCat + 1))            /* frequencies */
        goto LLTFin;
    if (alloc_acj(PMMaxCat + 1))            /* bf pointer */
        goto LLTFin;

    eerr = 0;

    printf1("Dimension  Variable ");
    prnchar(' ',VNameLen - 8,0);
    printf1("  Categories\n");

    for (i = 0; i < NDIM; ++i) {

        for (j = 0; j < NOC; ++j) {
            AcR[j + 1] = (int)get_data(PMVIdx[i],j);
        }
        n = cfreq(1,NOC,AcR,PMMaxCat,AcS,AcI,AcJ,0,&tmp,&tmp);

        if (n < 1) {
            printf1("Error: can't sort dimension %d.\n",i + 1);
            goto LLTFin;
        }
        if (n == 1)
            eerr++;

        CTNCat[i] = n;
        if (!(CTCat[i] = (int *) calloc(n + 1,sizeof(int)))) {   
            p_err(-2,1);
            goto LLTFin;
        }
        CTCatAI[i] = n + 1;
        memrq(n + 1,sizeof(int));

        NCTab *= n;

        printf1("%3d  %c     %s ",i + 1,'A' + i,VName[PMVIdx[i]]);
        prnchar(' ',VNameLen - strlen(VName[PMVIdx[i]]),0);
        printf1("  %3d : ",n);
        for (j = 1; j <= n; ++j) {
            k = AcS[AcJ[j]];
            CTCat[i][j] = k;
            printf1(PMNFmtS,k);
        }
        newline();
    }
    newline();
    if (eerr) {
        printf1("Error: only one value in at least one dimension.\n");
        goto LLTFin;
    }

    /*  build contingency table, value combinations in CTVal, frequencies
        in CTFreq, NDIM dimensions, NCTab cells */

    if (!(CTVal = (short *) calloc(NCTab * NDIM + 1,sizeof(short)))) {   
        p_err(-2,1);
        goto LLTFin;
    }
    CTValA = NCTab * NDIM + 1;
    memrq(CTValA,sizeof(short));

    if (!(CTFreq = (int *) calloc(NCTab + 1,sizeof(int)))) {   
        p_err(-2,1);
        goto LLTFin;
    }
    CTFreqA = NCTab + 1;
    memrq(CTFreqA,sizeof(int));

    if (PMSCAL >= 0) {
        if (!(CTSCAL = (double *) calloc(NCTab + 1,sizeof(double)))) {   
            p_err(-2,1);
            goto LLTFin;
        }
        CTSCALA = NCTab + 1;
        memrq(CTSCALA,sizeof(double));
    }

    for (k = 0; k < NDIM; ++k)
        AcS[k] = 1;

    i = 0;
    while (1) {
        for (k = 0; k < NDIM; ++k)  
            CTVal[++i] = (short) CTCat[k][AcS[k]];

        j = 0;
        for (k = NDIM - 1; k >= 0; --k) {
            if ((AcS[k] += 1) <= CTNCat[k]) {
                j = 1;
                break;     
            }
            else
                AcS[k] = 1;
        }
        if (j == 0)
            break;
    }
    n = 1;
    for (j = 0; j < NOC; ++j) {

        if (PMWVar >= 0)  
            n = (int)get_data(PMWVar,j);

        for (i = 0; i < NDIM; ++i)  
            AcS[i] = (int)get_data(PMVIdx[i],j);

        l = m = 1;
        for (i = NDIM - 1; i >= 0; --i) {
            for (k = 1; k <= CTNCat[i]; ++k) {
                if (AcS[i] == CTCat[i][k])
                    break;
            }
            l += (k - 1) * m;
            m *= CTNCat[i];
        }
        CTFreq[l] += n;
        CTNum += n;
        if (PMSCAL >= 0)
            CTSCAL[l] += get_data(PMSCAL,j);
    }
    printf1("Table has %d cells, %d counts.\n",NCTab,CTNum);
    if (CTNum < 1)  
        goto LLTFin;   

    CTIFlg = NCTab1 = 0;
    for (i = 1; i <= NCTab; ++i) {
        if (CTFreq[i] > 0)
            NCTab1++;
    }
    if (NCTab1 != NCTab) {       /* check if table is complete */
        CTIFlg = 1;
        printf1("Table is incomplete:");
        printf1(" %d (%5.2f%%) of %d cells.\n",NCTab1,
                                 100.0 * (double)NCTab1/(double)NCTab,NCTab);
    }
    else 
        printf1("Table is complete.\n",NCTab);

    err = 0;

LLTFin:
    alloc_acr(0);
    alloc_acs(0);
    alloc_aci(0);
    alloc_acj(0);
    return(err);
}  

/* ------------------------------------------------------------------------ */
/*  loglin1     Estimation of log-linear models                             */
/*              Return 0 if OK, -1 if error                                 */

int loglin1(void)
{
    register int i,j,k,l;
    register char *q;
    int r,err,slen;
    double tmp;

    err = -1;

    printf1("Begin of model estimation.\n\n");

    for (i = 0; i < PMModN; ++i) {

        LLDMCol = LLDCol[i];        /* columns of the design matrix */
        slen = 1;
        printf1("Model: ");         /* print model description */
          
        if (LLTyp[i] == 0)          /* if null model */
            printf1("0\n");
             
        else if (LLTyp[i] > 0) {

            prn_mdes(stdout,PMModS[i],7);
            expand1(PMModS[i],LLD1,1);     
            r = expand2(LLD1,LLD2,1);     
            if (r > PMModSA[i]) {
                printf1("Expanded: ");
                prn_mdes(stdout,LLD2,10);
            }

            /*  build the design matrix. slen is the max length of the model
                parameter labels. */

            r = build_dm(LLD2,LLDMCol,&slen,1);
            if (r < 0 || r != LLDMCol) { 
                ll_error(r,PMModS[i]);
                goto LL1Fin;
            }
        }
        else if (LLTyp[i] < 0) { 
            slen = 10;
    
            printf1("%s\n",PMModS[i]);

            r = ll_readdm(i,LLDMCol,1);

            if (r < 0) { 
                ll_error(r,PMModS[i]);
                goto LL1Fin;
            }
        }
        newline();     
            
        CTNPar = LLDMCol;   /* number of model parameters = number of
                               columns of the design matrix. */

        if (LLDMCol && PMF1Def)     /* if requested print design matrix */
            prn_lldmat(i);  

        printf1("Number of model parameters: %d\n",CTNPar + 1); 
        printf1("Maximum number of iterations: %d\n",MxIter);
        printf1("Tolerance for convergence: %lg\n\n",TOLF);
            
        r = logl_fit();     /* call function for estimation */
         
        printf1("Convergence ");
        if (r > MxIter)
            printf1("not ");
        printf1("reached in %d iterations.\n\n",r);

        if (CTLR < 1.e-6)
            CTLR = 0.0;
        if (CTCHI < 1.e-6)
            CTCHI = 0.0;

        printf1("Likelihood Ratio Statistic %14.4f",CTLR);
        if (CTNDF > 0)
            printf1("  Prob: %6.4f",1.0 - cdchif(CTLR,CTNDF));
        printf1("\nPearson's Chi Square       %14.4f",CTCHI);
        if (CTNDF > 0)
            printf1("  Prob: %6.4f",1.0 - cdchif(CTCHI,CTNDF));
        if (CTNDF > 0)
            printf1("\nF-Statistic                %14.4f",CTLR / (double)CTNDF);
        printf1("\nDegrees of Freedom %d\n\n",CTNDF);
             
        /*  calculate constant in CTPar[0] */

        tmp = 0.0;
        for (j = 1; j <= CTNPar; ++j)  
            tmp += (double)LLDMat[j] * CTPar[j];

        CTPar[0] = rlog(CTFit[1]) - tmp;

        /*  print estimated parameters and errors */

        if (slen < 10)
            slen = 10;

        printf1("Idx  Parameter");
        prnchar(' ',slen - 9,0);
        prnchar(' ',PMTFmt1 - 10,0); printf1("     Coeff");
        prnchar(' ',PMTFmt1 -  9,0); printf1("     Error");
        prnchar(' ',PMTFmt1 -  9,0); printf1("    T-Stat  Signif\n");
        prnchar('-',12 + slen + 3 * (PMTFmt1 + 1),0);

        q = CTVLab;
        l = 0;
        r = -1;
        for (j = 0; j <= LLDMCol; ++j) {
            printf1("\n%3d  ",j);
            if (LLTyp[i] > 0) {
                for (k = 0; k < slen; ++k) {
                    if (l < CTVLen && *q) {
                        printf1("%c",*q++);
                        l++;
                    }
                    else
                        printf1(" ");
                }
                l++;
                q++;
            }
            else if (j == 0)
                printf1("Constant  ");
            else
                printf1("Col%4d   ",j);

            r++;
            printf1(PMTFmtS,CTPar[r]);

            if (j && (tmp = CTMat[(r - 1) * CTNPar + r]) > 0.0) {
                tmp = sqrt(tmp);
                printf1(PMTFmtS,tmp);
                tmp = CTPar[r] / tmp;
                printf1(PMTFmtS,tmp);  
                tmp = 2.0 * cdnf(fabs(tmp)) - 1.0;
                printf1(" %6.4f",tmp);
            }
            else {
                prnchar(' ',PMTFmt1 - 3,0); printf1("--- ");
                prnchar(' ',PMTFmt1 - 3,0); printf1("---     ---");
            }
        }
        printf1("\n\n");

        if (PMPPFDef)       /* write parameter estimates */
            ll_ppar();

        if (PMResFDef) {    /* if requested print table of residuals */

            fprintf(PMResFd,"Residuals of Model: ");
            if (LLTyp[i] > 0)
                prn_mdes(PMResFd,PMModS[i],20);
            else
                fprintf(PMResFd,"%s\n",PMModS[i]+2);

            prn_ctab(2);
            printf1("Residuals written to: %s\n",PMResFName);
        }
        if (LLDMCol && PMCovFDef) {    /* if requested print covariance matrix */
            prn_llcov(i);
            printf1("Covariance matrix written to: %s\n",PMCovFName);
        }
        printf1("\n");
    }
    err = 0;

LL1Fin:
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  ll_ppar()              Write parameter to PMPPFd.                       */

void ll_ppar(void)
{
    register int i;

    for (i = 0; i <= LLDMCol; ++i) {
        fprintf(PMPPFd,PMTFmtS,CTPar[i]);
        fprintf(PMPPFd,"\n");
    }
    printf1("Parameter estimates written to: %s\n",PMPPFName);
}

/* ------------------------------------------------------------------------ */
/*  ll_error    Print error message.                                        */

void ll_error(int err,char *mod)
{
    printf1("Error: %s\n",mod);

    switch (err) {
        case  -1:  printf1("Syntax error.\n");
                   break;
        case  -2:  printf1("Can't find at least one dimension.\n");
                   break;
        case  -3:  printf1("Can't read the level of a dimension.\n");
                   break;
        case  -4:  printf1("At least one level is out of range.\n");
                   break;
        case  -5:  printf1("Exceeded max number of table dimensions.\n");
                   break;
        case  -6:  printf1("Exceeded max space for model expansion.\n");
                   break;
        case  -7:  printf1("Error in build_dm() algorithm.\n");
                   break;
        case  -8:  printf1("Can't open the file.\n");
                   break;
        case  -9:  printf1("File is empty.\n");
                   break;
        case -10:  printf1("File contains errors.\n");
                   break;
        case -11:  printf1("One or more lines have less than %d entries.\n",LLDMCol);
                   break;
        case -12:  printf1("Found more than %d lines in the file.\n",NCTab);
                   break;
        default:   break;
    }
}

/* ------------------------------------------------------------------------ */
/*  prn_mdes     Print model description desc to file fd.                   */

void prn_mdes(FILE *fd, char *desc, int n)
{
    register char *p;

    p = desc;
    while (*p)  
        fprintf(fd,"%c",*p++);
    fprintf(fd,"\n");
}

/* ------------------------------------------------------------------------ */
/*  prn_lldmat     Print design matrix of model i to PMF1d.                 */

void prn_lldmat(int i)
{
    register int j,k;

    fprintf(PMF1d,"# Design matrix of model: ");
    if (LLTyp[i] > 0)
        prn_mdes(PMF1d,PMModS[i],24);
    else if (LLTyp[i] < 0)
        fprintf(PMF1d,"%s\n",PMModS[i]+2);
   
    for (j = 1; j <= NCTab; ++j) {
        for (k = 1; k <= LLDMCol; ++k)
            fprintf(PMF1d,PMNFmtS,(int)LLDMat[(j - 1) * LLDMCol + k]);
        fprintf(PMF1d,"\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  prn_llcov      Print covariance matrix of model i to PMCovFd.           */

void prn_llcov(int i)
{
    register int j,k;

    fprintf(PMCovFd,"# Covariance matrix of model: ");
    if (LLTyp[i] > 0)
        prn_mdes(PMCovFd,PMModS[i],28);
    else
        fprintf(PMCovFd,"%s\n",PMModS[i]+2);

    for (j = 1; j <= CTNPar; ++j) {
        for (k = 1; k <= CTNPar; ++k)  
            fprintf(PMCovFd,PMMFmtS,CTMat[(j - 1) * CTNPar + k]);
        fprintf(PMCovFd,"\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  prn_ctab (mode)                                                         */
/*      mode 1 : Print contingency table to PMTabFd                         */
/*      mode 2 : Print residuals to PMResFd                                 */

void prn_ctab(int mode)
{
    FILE *ofd;
    register int j,k;
    double s = 1.0;

    if (mode == 1)
        ofd = PMTabFd;
    else
        ofd = PMResFd;

    if (mode == 2) {
        fprintf(ofd,"\nIndex  ");
        for (k = 0; k < NDIM; ++k) {
            fprnchar(ofd,' ',PMNFmt - 1,0);
            fprintf(ofd,"%c ",'A' + k);
        }
        fprintf(ofd,"  Observed ");
        fprintf(ofd,"      Fitted ");
        fprintf(ofd,"    Residual        Scale\n");
        fprnchar(ofd,'-',56 + NDIM * (PMNFmt + 1),1);
    }
    for (k = 1; k <= NCTab; ++k) {
        fprintf(ofd,"%5d  ",k);
        for (j = 1; j <= NDIM; ++j)  
            fprintf(ofd,PMNFmtS,(int)CTVal[(k - 1) * NDIM + j]);
        fprintf(ofd,"%10d ",CTFreq[k]);

        if (mode == 2)  
            fprintf(ofd,"%12.4f %12.4f ",CTFit[k],CTRes[k]);

        if (PMSCAL >= 0)
            s = CTSCAL[k];

        if (PMSCAL >= 0 || mode == 2)
            fprintf(ofd,"%12.4f",s);

        fprintf(ofd,"\n");
    }
    if (mode == 2)  
        fprintf(ofd,"\n");
}

/* ------------------------------------------------------------------------ */
/*  llmodcheck   Check model description.                                   */
/*      Return:  1  if no error                                             */
/*               0  if empty string (null model)                            */
/*              -1  if syntax error                                         */
/*              -2  if one or more dimensions not found                     */
/*              -3  if a level cannot be read                               */
/*              -4  if a level is out of range                              */

int llmodcheck(char *desc)
{
    register int k,l,fnd;
    register char *p;
    int n,llev,rlev,nt;

    if (!*desc)
        return(0);

    p = desc;
    nt = llev = rlev = 0;
    while (*p) {
        if (*p == '(')
            llev++;
        else if (*p == ')')
            rlev++;
        else if (*p == '+' || *p == '.') {
            if (!nt)                    
                return(-1);
        }
        else if (*p >= 'A' && *p <= 'Z') {
            k = (int)*p - (int)'A';
            if (k >= NDIM)
                return(-2);
            nt++;
            if (*++p == '[') {
                if (sscanf(++p,"%d",&n) != 1)   
                    return(-3);
                if (*p == '-')
                    p++;
                while (*p && *p >= '0' && *p <= '9')
                    p++;
                if (*p != ']') 
                    return(-1);

                fnd = 0;
                for (l = 1; l < CTNCat[k]; ++l) {
                    if (n == CTCat[k][l]) {
                        fnd = 1;
                        break;
                    }
                }
                if (!fnd)
                    return(-4);
            }
            else
                p--;
        }
        else
            return(-1);
        p++;
    }
    if (!nt || llev != rlev)
        return(-1);
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  expand1(s,r,mode)                                                       */
/*      The input string s is assumed to be the description of a log linear */
/*      model. If the string contains expressions like AB, ABC, ..., these  */
/*      expressions are expanded: A+B+A.B, ...                              */
/*                                                                          */
/*      If mode = 0 the function returns the length of the expanded string, */
/*      if mode = 1 the string is actually expanded and returned in r.      */
/*      The function returns                                                */
/*          -5, if more than MaxABC characters are used in an ABC.. string  */

int expand1(char *s, char *r, int mode)
{
    int j,i,k,last,len = 0;
    int iset[MaxABC + 1],dim[MaxABC + 1];
    register char *p,*q;

    p = s;
    if (!*p)
        return(0);

    while (*p) {
        j = 0;
        if (*p >= 'A' && *p <= 'Z') {

            dim[++j] = (int)*p - (int)'A';
            q = p + 1;
            while (*q >= 'A' && *q <= 'Z') {
                if (++j > MaxABC)
                    return(-5);
                dim[j] = (int)*q - (int)'A';
                q++;
            }
            if (j > 1) {
                if (mode)
                    *r++ = '(';
                len++;

                for (i = 1; i <= j; ++i) {
                    last = 1;
                    while (!(last = logscr_comb(iset,j,i,last))) {
                        for (k = 1; k <= i; ++k) {
                            if (k > 1) {
                                if (mode)
                                    *r++ = '.';
                                len++;
                            }
                            if (mode)
                                *r++ = 'A' + dim[iset[k]];
                            len++;
                        }
                        if (mode)
                            *r++ = '+';
                        len++;
                    }
                }
                r--;
                if (mode)
                    *r++ = ')';
                p = q;
            }
            else {
                if (mode)
                    *r++ = *p;
                p++;
                len++;
            }
        }
        else {
            if (mode)
                *r++ = *p;
            p++;
            len++;
        }
    }
    if (mode)
        *r = '\0';
    return(len);
}

/* ------------------------------------------------------------------------ */
/*  expand2(s,r,mode)                                                       */
/*      The input string s is assumed to be an "algebraic" expression       */
/*      that consists of upper case letters (A - Z), operands '+' and '.',  */
/*      and brackets. If mode = 1 the string is expanded in r, else only    */
/*      the length of the expanded string is calculated.                    */
/*                                                                          */
/*      The function returns the length of the expanded string, or          */
/*          -1  if a syntax error is deteced                                */
/*          -3  if a level cannot be read                                   */
/*          -6 if max space is exceeded                                     */
/*                                                                          */
/*      The algorithm is adapted from:  M.J. Levine, Symbolic Expansion     */
/*      of Algebraic Expressions, Communications of the ACM, 13 (1970),     */
/*      No. 3, 191 - 192.                                                   */

int expand2(char *s, char *r, int mode)
{
    register int n,lvl,utyp,ttyp,plus = 0;
    register char t,uvar,tvar;
    int len,uvarc,tvarc;

    if (!*s)
        return (0);

    len = utyp = uvarc = lvl = n = 0;
    uvar = '\0';
   
    if (*s >= 'A' && *s <= 'Z') {
        uvar = *s;
        if (*++s == '[') {
            if (sscanf(++s,"%d",&uvarc) != 1)  
                return(-3);
            if (*s == '-')
                s++;
            while (*s && *s >= '0' && *s <= '9')
                s++;
            if (*s != ']') 
                return(-1);
            utyp = 2;
        }
        else {
            s--;
            utyp = 1;
        }
    }
    while (1) {

        t = *s;
        ttyp = utyp;
        tvar = uvar;
        tvarc = uvarc;

        if (!*++s)
            utyp = -1;

        else {
            if (*s == '+' && ttyp > 0)
                utyp = 0;
            else if (*s >= 'A' && *s <= 'Z') {
                uvar = *s;
                if (*++s == '[') {
                    if (sscanf(++s,"%d",&uvarc) != 1)   
                        return(-3);
                    if (*s == '-')
                        s++;
                    while (*s && *s >= '0' && *s <= '9')
                        s++;
                    if (*s != ']') 
                        return(-1);
                    utyp = 2;
                }
                else {
                    s--;
                    utyp = 1;
                }
            }
            else   
                utyp = 0;
        }
        if (ttyp < 0)
            t = '+';

        if (ttyp < 0 || utyp > 0 || (utyp == 0 && *s == '(')) {
  
            if (t == '.') {
                ExMult[n] = 1;
                ExOpl[n] = (char)lvl;
                if (++n > ExMax)
                    return(-6);
            }
            else if (t == '+') {
                ExMult[n] = 0;
                ExOpl[n] = (char)lvl;
                if (lvl == 0) {
                    if (plus) {
                        if (mode)
                            *r++ = '+';
                        len++;
                    }
                    r = distribute(n,r,mode,&len);        
                    plus++;
                    n = 0;
                }
                else if (++n > ExMax)
                    return(-6);
            }
            else if (t == '(')
                lvl++;
            else  
                return (-1);
        }
        else {
            if (t == ')' && lvl > 0)
                lvl--;
            else if (ttyp > 0) {
                ExV[n] = tvar;
                ExVTyp[n] = ttyp;
                ExVc[n] = tvarc;
                ExVl[n] = (char)lvl;
            }
            else
                return (-1);
        }
        if (ttyp < 0)  
            break;
    }
    if (lvl != 0)
        return (-1);
    if (mode)
        *r = '\0';
    return(len);
}

/* ------------------------------------------------------------------------ */
/*  distribute(), used by expand2().                                        */

char *distribute(int n, char *r, int mode, int *len)
{
    register int i,j,k,l,alter,level,product,term,plus = 0;
    char sstmp[80];

    for (k = 0; k <= n; ++k)
        ExUsed[k] = 0;
NEXT:
    alter = 1;
    j = i = -1;
FACTOR:
    if (ExUsed[++i])
        goto FACTOR;
    ExIdx[++j] = i;

    while (1) {
        if (ExMult[i])
            goto FACTOR;
        level = (int)ExOpl[i];
     
        if (level <= 0)
            break;
        else {
            if (alter) {
                l = level;
                level = (int)ExVl[i] + 1;
                ExUsed[i] = product = term = 1;
                alter = 0;
                for (k = i - 1; k >= 0; k--) {
                    if ((int)ExOpl[k] < level) {
                        level = (int)ExOpl[k];
                        product = ExMult[k];
                        if (product)
                            level++;
                        if (level <= l)
                            term = 0;
                    }
                    if (product)
                        ExUsed[k] = (char)term;
                }
            }
            else {
                while (1) {
                    if (level > (int)ExOpl[++i])
                        break;
                }
            }
        }
    }         

    /*  build summands, the factors are ExV[ExIdx[k]], k = 0,j */
   
    if (plus) {
        *len += 1;
        if (mode)
            *r++ = '+';
    }
    for (k = 0; k <= j; ++k) {
        if (k) {
            *len += 1;
            if (mode)
                *r++ = '.';
        }
        *len += 1;
        if (mode)  
            *r++ = ExV[ExIdx[k]];

        if (ExVTyp[ExIdx[k]] == 2) {
            sprintf(sstmp,"[%d]",ExVc[ExIdx[k]]);
            *len += strlen(sstmp);
            if (mode) {
                strcpy(r,sstmp);
                r += strlen(sstmp);
            }
        }
    }
    if (!alter) {
        plus++;
        goto NEXT;
    }
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  build_dm(descr,col,slen,mode)                                           */
/*      If mode == 0 calculate number of columns of the design matrix,      */
/*      else create the design matrix in LLDMat. It is assumed that memory    */
/*      for the design matrix is already allocated and has space for col    */
/*      columns. Also, if mode != 0, labels of the model parameters are     */
/*      created in the string CTVLab, and slen is set to the maximum        */
/*      length of these labels.                                             */
/*                                                                          */
/*      The function returns the number of columns of the design matrix     */
/*      or an error code:                                                   */
/*          -1  if syntax error                                             */
/*          -2  if a dimension is not defined in the table                  */
/*          -3  if a level of a dimension cannot be read                    */
/*          -4  if a level of a dimension is out of range                   */
/*          -5  if exceeded max number of dimensions                        */
/*          -7  if more than col columns are needed.                        */

int build_dm(char *descr, int col, int *slen, int mode)
{
    register int i,j,k,l,m,m1,m2;
    int len,nsl,fnd,r,n,dn,x,ptr[MaxDX + 1];
    int cat[MaxDX + 1],ctyp[MaxDX + 1],clev[MaxDX + 1];
    register char *p,*q,*ctvptr;
    char sstmp[20];

    if (mode) {
        ctvptr = CTVLab;
        *ctvptr = '\0';
        strcat(ctvptr,"Constant");
        ctvptr += 9;
        *ctvptr = '\0';
        CTVLen = 9;
        *slen = 8;
        if (!col)   
            return(0);
    }
    p = descr;
    j = r = 0;

    while (*p) {

        nsl = dn = 0;
        m = 1;
        while (*p) {
            if (*p < 'A' || *p > 'Z')   
                return(-1);
            k = (int)(*p - 'A');
            if (k < 0)  
                return(-2);
            if (++dn > MaxDX)    
                return(-5);
            ptr[dn] = k;
            if (*++p == '[') {
                if (sscanf(++p,"%d",&n) != 1)
                    return(-3);
                fnd = 0;
                for (l = 1; l < CTNCat[k]; ++l) {
                    if (n == CTCat[k][l]) {
                        fnd = 1;
                        break;
                    }
                }
                if (!fnd)
                    return(-4);

                ctyp[dn] = 1;
                clev[dn] = n;
                nsl++;
          
                if (*p == '-')
                    p++;
                while (*p && *p >= '0' && *p <= '9')
                    p++;
                if (*p++ != ']')
                    return(-1);
            }
            else {
                ctyp[dn] = 0;
                m *= CTNCat[k] - 1;
            }
            if (!*p || *p == '+')
                break;
            else if (*p != '.')  
                return(-1);
            p++;
        }
        if (*p == '+' && !*++p)
            return(-1);

        r += m;

        if (mode) {     /* create m columns of design matrix */

            for (k = 1; k <= dn; ++k)
                cat[k] = 1;

            while (1) {

                /*  check for single levels */

                n = 0;
                for (k = 1; k <= dn; ++k) {
                    if (ctyp[k]) {
                        m = ptr[k];
                        m2 = cat[k];
                        if (clev[k] == CTCat[m][m2])  
                            n++;
                    }
                }
                if (!nsl || nsl == n) {

                    j++;
                    if (j > col)  
                        return(-7);
                    len = 0;
    
                    for (i = 1; i <= NCTab; ++i) {

                        n = 1;
                        for (k = 1; k <= dn; ++k) {
                            if (n) {
                                m = ptr[k];
                                x = (int) CTVal[(i - 1) * NDIM + m + 1];
                                m1 = CTNCat[m];
                                m2 = cat[k];
                                if (x == CTCat[m][m1])
                                    n *= -1;
                                else if (x != CTCat[m][m2])
                                    n = 0;
                            }
                                /* create label of parameter */

                            if (i == 1 && CTVLen < CTVLenM) {
                                m = ptr[k];
                                m2 = cat[k];
                                if (k > 1) {
                                    *ctvptr++ = '.';
                                    CTVLen++;
                                    len++;
                                }
                                if (CTVLen < CTVLenM) {
                                    *ctvptr++ = 'A' + (char)m;
                                    len++;
                                    CTVLen++;
                                    sprintf(sstmp,"[%d]",CTCat[m][m2]);
                                    q = sstmp;
                                    while (*q && CTVLen < CTVLenM) {
                                        *ctvptr++ = *q++;
                                        len++;
                                        CTVLen++;
                                    }
                                }
                            }
                        }
                        LLDMat[(i - 1) * col + j] = (short) n;
                    }
                    if (CTVLen < CTVLenM) {
                        *ctvptr++ = '\0';
                        CTVLen++;
                        if (*slen < len)
                            *slen = len;
                    }
                }
                if (design_comb(dn,cat,ptr))   /* next combination */
                    break;
            }
        }
    }
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  design_comb(dn,cat,ptr)                                                 */
/*      Create combinations of the elements of the array cat[1,...,dn].     */
/*      dn is the number of dimensions. The field ptr[1,...,dn] contains    */
/*      pointers so that CTNCat[ptr[j]] is the category associated with     */
/*      this dimension.                                                     */
/*                                                                          */
/*      Return 0 if a new combination is found, else return 1.              */

int design_comb(int dn, int *cat, int *ptr) 
{
    register int j;

    for (j = dn; j >= 1; --j) {
        if ((cat[j] += 1) < CTNCat[ptr[j]])
            return(0);
        else
            cat[j] = 1;
    }
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  ll_readdm(i,dcol,mode)                                                  */
/*      Read a design matrix file. If mode = 0 calculate number of columns  */
/*      as found in the first (non-comment) line of the file. Else read     */
/*      up to NCTab lines into the design matrix LLDMat.                    */
/*                                                                          */
/*      Return:  the number of entries in the first line, if mode = 0, or   */
/*               the number of lines read into the design matrix, or        */
/*               -8  if the file cannot be openend,                         */
/*               -9  if the file is empty,                                  */
/*              -10  if a syntax error is detected                          */
/*              -11  if mode = 1 and less than dcol entries in a line,      */
/*              -12  if mode = 1 and more than NCTab lines in the file.     */

int ll_readdm(int i, int dcol, int mode)
{
    FILE *fd;
    register int col,row;
    int n;
    register char *p;

    if (alloc_acc(RLMaxDef + 1))    
        return(-2);   

    if (!(fd = fopen(PMModS[i] + 2,OPEN_RD)))   
        return(-8);
    else {
        row = col = 0;
        while (fgets(AcC,RLMaxDef,fd)) {
            p = skip_b(AcC);

            col = 0;
            if ((p = check_comment(AcC)) != NULL) {   /* if no comment */

                while (*p && *p != '\n') {    
                    if (sscanf(p,"%d",&n) != 1) { 
                        fclose(fd);
                        return(-10);
                    }
                    if (!col)
                        row++;
                    col++;
                    if (mode && col <= dcol) {
                        if (row > NCTab) {
                            fclose(fd);
                            return(-12);
                        }
                        LLDMat[(row - 1) * dcol + col] = (short) n; 
                    }
                    p = skip_int(p);
                    while (*p && *p == ' ')
                        p++;
                }
                if (mode) {
                    if (col < dcol) {
                        fclose(fd);
                        return(-11);
                    }
                }
                else if (col)
                    break;
            }
        }
        fclose(fd);
        if (!row)
            return(-9);
        else if (mode)
            col = dcol;
        return(col);
    }
}

/* ------------------------------------------------------------------------ */
/*  logl_fit    Fitting a log-linear model.                                 */
/*      The algorithm is adapted from the program FREQ by:                  */
/*      S.J. Haberman, Analysis of Qualitative Data, Vol 1 & 2, NY 1979.    */
/*                                                                          */
/*      Return: Number of iterations. Successful if less than MaxIt.        */

int logl_fit(void) 
{
    register int i,k,l,ii,ll;
    int it,irank;
    double tmp,tmp1,tmp2,e,f,sum,sum1,w;

    /* initialize the algorithm */

    f = DBLMAX;
    CTNDF = NCTab - 1;
    tmp = sum = sum1 = 0.0;
    for (i = 1; i <= NCTab; ++i) {
        CTFit[i] = (double)CTFreq[i] + 0.5;
        CTRes[i] = rlog(CTFit[i]);

        if (PMSCAL >= 0) {
            if (CTSCAL[i] > EPSI1)      
                CTRes[i] -= rlog((double)CTSCAL[i]);
            else {
                CTFit[i] = 0.0;
                CTNDF--;
            }
        }
        sum += CTFit[i];
        sum1 += CTFit[i] * CTRes[i];
    }
    if (sum > TOLF) {
        w = sum1 / sum; 
        for (i = 1; i <= NCTab; ++i)
            CTRes[i] -= w;
    }
    for (k = 1; k <= CTNPar; ++k) {
        CTPar[k] = sum = 0.0;
        for (i = 1; i <= NCTab; ++i)   
            sum += CTFit[i] * CTRes[i] * (double)LLDMat[(i - 1) * LLDMCol + k];
        CTWork[k] = sum;
    }

    /* perform an iteration, maximal MxIter iterations */

    for (it = 1; it <= MxIter; ++it) {
        /*********************
        if (it > 1 || init1) {
        *********************/
            sum = 0.0;
            for (i = 1; i <= NCTab; ++i)
                sum += CTFit[i];
            w = sum;
            for (k = 1; k <= CTNPar; ++k) {
                sum = 0.0;
                for (i = 1; i <= NCTab; ++i)  
                    sum += CTFit[i] * (double)LLDMat[(i - 1) * LLDMCol + k];
                if (w > TOLF)
                    CTB[k] = sum / w;
                else
                    CTB[k] = 0.0;
            }

            /* obtain weighted sums of cross-products */

            for (k = 1; k <= CTNPar; ++k) {
                for (l = 1; l <= k; ++l) {
                    sum = 0.0;
                    for (i = 1; i <= NCTab; ++i)  
                        sum += ((double)LLDMat[(i - 1) * LLDMCol + k] - CTB[k]) *
                               ((double)LLDMat[(i - 1) * LLDMCol + l] - CTB[l]) *
                                                                      CTFit[i]; 
                    CTMat[(k - 1) * CTNPar + l] = sum;
                }
            }

            /* obtain Cholesky decomposition of s */

            irank = logl_chol(CTNPar,CTMat,CTMat);

            /* see if further steps are needed */

            if (f < TOLF) {
                logl_inv(CTNPar,CTMat,CTMat);
                CTNDF -= irank;
                goto LF_Fin;
            }

            /* obtain new value of b */

            logl_solve(CTNPar,CTMat,CTWork,CTWork);
            for (k = 1; k <= CTNPar; ++k) 
                CTPar[k] += CTWork[k];

            /* update U and CTFit and check for convergence */

            if (it > 1)
                f = 0.0;
        /**
        }
        **/
        for (i = 1; i <= NCTab; ++i) {
            sum = 0.0;
            if (CTNPar > 0) {
                for (k = 1; k <= CTNPar; ++k)
                    sum += CTPar[k] * (double)LLDMat[(i - 1) * LLDMCol + k];
                if (it > 1) {
                    e = fabs(sum - CTRes[i]);   
                    if (e > f)
                        f = e;
                }
            }
            CTRes[i] = sum;
            CTFit[i] = rexp(CTRes[i]);        

            if (PMSCAL >= 0) {
                if (CTSCAL[i] > EPSI1)
                    CTFit[i] *= CTSCAL[i];
                else
                    CTFit[i] = 0.0;         
            }
        }
        sum = w = 0.0;
        for (i = 1; i <= NCTab; ++i) {
            if (PMSCAL < 0 || CTSCAL[i] > EPSI1) {
                w += (double)CTFreq[i];
                sum += CTFit[i];
            }
        }
        if (sum > TOLF)
            w /= sum;
        else
            w = 0.0;
        for (i = 1; i <= NCTab; ++i)
            CTFit[i] *= w;

        if (CTNPar <= 0)
            goto LF_Fin;

        /* prepare differences between fitted and observed lin. combinations */

        for (k = 1; k <= CTNPar; ++k) {
            sum = 0.0;

            for (i = 1; i <= NCTab; ++i) {

                if (CTFit[i] > TOLF) {
                    tmp = (double)CTFreq[i];
                    sum += (tmp - CTFit[i]) * (double)LLDMat[(i - 1) * LLDMCol + k];
                }
            }
            CTWork[k] = sum;
        }
    }
    logl_inv(CTNPar,CTMat,CTMat);
    CTNDF -= irank;
    
    /* calc chi-square and likelihood ratio statistics */

LF_Fin:
    CTLR = CTCHI = 0.0;
    for (i = 1; i <= NCTab; ++i) {

        if (CTFit[i] > TOLF) {
            tmp1 = (double)CTFreq[i];
            tmp2 = tmp1 - CTFit[i];
            CTCHI += tmp2 * tmp2 / CTFit[i];
            if (tmp1 > TOLF)  
                CTLR += tmp1 * rlog(tmp1 / CTFit[i]);       
        }
    }
    CTLR *= 2.0;

    /* calculation of adjusted residuals in array CTRes */

    sum1 = 0.0;
    for (i = 1; i <= NCTab; ++i)
        sum1 += CTFit[i];

    for (i = 1; i <= NCTab; ++i) {
        ii = (i - 1) * LLDMCol;
        sum = 0.0;
        if (CTNPar > 1) {
            for (k = 2; k <= CTNPar; ++k) {
                ll = k - 1;
                for (l = 1; l <= ll; ++l)
                sum += CTMat[(k - 1) * CTNPar + l] * 
                                            ((double)LLDMat[ii + k] - CTB[k]) *
                                            ((double)LLDMat[ii + l] - CTB[l]);
            }
            sum *= 2.0;
        }
        for (k = 1; k <= CTNPar; ++k) {
            tmp  = (double)LLDMat[ii + k] - CTB[k];
            sum += CTMat[(k - 1) * CTNPar + k] * tmp * tmp;
        }
        sum = CTFit[i] * (1.0 - CTFit[i] * sum);
 
        if (sum1 > TOLF)
            sum -= CTFit[i] * CTFit[i] / sum1;
        CTRes[i] = 0.0;

        if (sum > TOLF) {
            tmp = (double)CTFreq[i];
            CTRes[i] = (tmp - CTFit[i]) / sqrt(sum);   
        }
    }
    return(it);     /* return number of iterations */
}

/* ------------------------------------------------------------------------ */
/*  logl_solve (n,lu,b,x)                                                   */
/*      Solve lu * x = b with lu(i,i) != 0                                  */
/*      n is the dimension of the matrices.                                 */

void logl_solve(int n, double *lu, double *b, double *x)
{
    register int i,j,k,i1;
    int in;
    double sum;

    x[1] = b[1];
    if (lu[1] <= 0.0)
        x[1] = 0.0;
    for (i = 2; i <= n; ++i) {
        sum = 0.0;
        i1 = i - 1;
        for (j = 1; j <= i1; ++j) 
            sum += lu[(i - 1) * n + j] * x[j];
        x[i] = b[i] - sum;
        if (lu[1] <= 0.0)
            x[i] = 0.0;
    }
    i = n;
    if (lu[n * n] > TOLF) 
        x[i] /= lu[n * n];
    else
        x[i] = 0.0;
    if (n == 1)
        return;

    for (k = 2; k <= n; ++k) {
        sum = 0.0;
        i1 = i;
        i--;
        in = (i - 1) * n;
        if (lu[in + i] <= 0.0)
            x[i] = 0.0;
        else {
            for (j = i1; j <= n; ++j)
                sum += lu[in + j] * x[j];
            x[i] = (x[i] - sum) / lu[in + i];
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  logl_inv (n,tri,xinv)                                                   */
/*      Calculate inverse of tri with a modified Cholesky decomposition.    */
/*      The inverse is returned in xinv. tri and xinv have dimenension n,n  */

void logl_inv(int n, double *tri, double *xinv)
{
    register int i,i1,k,j;
    int k1,in,jn;
    double sum;

    for (i = 1; i <= n; ++i) {
        in = (i - 1) * n;
        if (tri[in + i] > TOLF)
            xinv[in + i] = 1.0 / tri[in + i];
        else
            xinv[in + i] = 0.0;
    }
    if (n == 1)
        return;

    for (i = 2; i <= n; ++i) {
        in = (i - 1) * n;
        i1 = i - 1;
        for (k = 1; k <= i1; ++k) {
            sum = -tri[in + k];
            if (k < i1) {
                k1 = k + 1;
                for (j = k1; j <= i1; ++j)
                    sum -= tri[in + j] * xinv[(j - 1) * n + k];
            }
            xinv[in + k] = sum;
            xinv[(k - 1) * n + i] = sum * xinv[in + i];
        }
    }
    for (i = 1; i <= n; ++i) {
        i1 = i + 1;
        for (j = 1; j <= i; ++j) {
            jn = (j - 1) * n;
            sum = xinv[jn + i];
            if (i != n) {
                for (k = i1; k <= n; ++k)
                    sum += xinv[(k - 1) * n + i] * xinv[jn + k];
            }
            xinv[(i - 1) * n + j] = sum;
        }
    }
    for (i = 2; i <= n; ++i) {
        i1 = i - 1;
        for (j = 1; j <= i1; ++j)
            xinv[(j - 1) * n + i] = xinv[(i - 1) * n + j];
    }
}

/* ------------------------------------------------------------------------ */
/*  logl_chol(n,sym,tri)                                                    */
/*      Modified Cholesky decomposition for n,n -matrix sym (not negative   */
/*      definite. Output matrix is tri, return the rank irank.              */

int logl_chol(int n, double *sym, double *tri)
{
    register int i,j,k,l;
    int in,jn,irank;
    double x,sum,tol;

    tol = 1.0e-4;   

    irank = 0;
    for (i = 1; i <= n; ++i) {
        in = (i - 1) * n;
        x = tol * sym[in + i];
        tri[i] = sym[in + 1];
        if (i > 1) {
            if (tri[1] > TOLF)
                tri[in + 1] = tri[i] / tri[1];
            for (j = 2; j <= i; ++j) {
                jn = (j - 1) * n;
                if (tri[jn + j] > TOLF) {
                    sum = 0.0;
                    k = j - 1;
                    for (l = 1; l <= k; ++l) 
                        sum += tri[in + l] * tri[(l - 1) * n + j];
                    tri[jn + i] = sym[in + j] - sum;
                    if (j < i)
                        tri[in + j] = tri[jn + i] / tri[jn + j];
                }
            }
        }
        if (tri[in + i] > x)
            irank++;
        else {
            for (j = 1; j <= n; ++j)
                tri[in + j] = tri[(j - 1) * n + i] = 0.0;
        }
    }
    return(irank);
}

/* ------------------------------------------------------------------------ */
/*  logl_gen    Calculate generalized residuals.                            */

/************
void logl_gen(double *c, double *obs, double *expv, double *res)
{
    register int i,k,l;
    int ll;
    double sum,sum1,sum2,st;

    sum = sum1 = 0.0;
    for (i = 1; i <= NCTab; ++i) {
        sum1 += CTFit[i] * c[i];
        sum  += CTFit[i];
    }
    sum2 = sum1 * sum1 / sum;

    if (CTNPar > 0) {
        for (k = 1; k <= CTNPar; ++k) {
            sum = 0.0;
            for (i = 1; i <= NCTab; ++i)  
                sum += ((double)LLDMat[(i-1) * LLDMCol + k] - CTB[k]) *
                                                             CTFit[i] * c[i];
            CTWork[k] = sum;
        }
        sum = 0.0;
        for (k = 2; k <= CTNPar; ++k) {
            ll = k - 1;
            for (l = 1; l <= ll; ++l)
                sum += CTMat[(k - 1) * CTNPar + l] * CTWork[k] * CTWork[l];
        }
        sum += sum;
        for (k = 1; k <= CTNPar; ++k) 
            sum += CTMat[(k - 1) * CTNPar + k] * CTWork[k] * CTWork[k];
    }
    sum1 = 0.0;
    for (i = 1; i <= NCTab; ++i)  
        sum1 += CTFit[i] * c[i] * c[i];
    st = sum1 - sum - sum2;
    *res = sum = sum1 = 0.0;
    for (i = 1; i <= NCTab; ++i) {
        sum += (double)CTFreq[i] * c[i];
        sum1 +=  CTFit[i] * c[i];
    }
    *obs = sum;
    *expv = sum1;     
    if (st > 0.0) {
        st = sqrt(st);
        *res = (*obs - *expv) / st;
    }
}
***/

/* -##--------------------------------------------------------------------- */
/*  logscrn    Test of marginal and partial association.                    */
/*                                                                          */
/*      Ref.: E.D. Lustbader, R.K. Stodola, Partial and Marginal            */
/*      Association in Multidimensional Contingency Tables, AS 160          */
/*      Applied Statistics 30 (1981), 97 - 105; M.B. Brown, Screening       */
/*      Effects in Multidimensional Contingency Tables, Applied             */
/*      Statistics 25 (1976), 37 - 46.                                      */
/*                                                                          */
/*      Return 0 if OK, -1 if error.                                        */

#define MaxDev  0.25    /* used in logfit */

int logscrn(void) 
{
int last;
    register int i,j,m,jn,jn1;
    int err,df,iemax,mm,np,np1,ntab2,l3;
    short itp;
    double g21,g22,g23;

    err = -1;

    printf1("Screening the table.\n");

    /*  calculate iemax, i.e. the maximal number of interaction effects
        in the saturated model. */

    NDIM1 = NDIM - 1;
    iemax = 1;
    j = NDIM / 2;
    for (i = 1; i <= j; ++i)
        iemax = (iemax * (NDIM - i + 1)) / i;

    if (!(CTFit = (double *) calloc(NCTab + 1,sizeof(double)))) {
        p_err(-2,1);
        goto LSFin;
    }
    CTFitA = NCTab + 1;
    memrq(CTFitA,sizeof(double));

    if (!(GSQ  = (double *)calloc(NDIM + 1,sizeof(double)))) {  
        p_err(-2,1);
        goto LSFin;
    }
    GSQA = NDIM + 1;
    memrq(GSQA,sizeof(double));

    if (!(PART = (double *)calloc(iemax * NDIM1 + 1,sizeof(double)))) {
        p_err(-2,1);
        goto LSFin;
    }
    PARTA = iemax * NDIM1 + 1;
    memrq(PARTA,sizeof(double));

    if (!(MARG = (double *)calloc(iemax * NDIM1 + 1,sizeof(double)))) {
        p_err(-2,1);
        goto LSFin;
    }
    MARGA = iemax * NDIM1 + 1;
    memrq(MARGA,sizeof(double));

    if (!(DGFR = (int *)calloc(NDIM + 1,sizeof(int)))) {
        p_err(-2,1);
        goto LSFin;
    }
    DGFRA = NDIM + 1;
    memrq(DGFRA,sizeof(int));

    if (!(DFS  = (int *)calloc(iemax * NDIM1 + 1,sizeof(int)))) {
        p_err(-2,1);
        goto LSFin;
    }
    DFSA = iemax * NDIM1 + 1;
    memrq(DFSA,sizeof(int));

    j = iemax * NDIM1;
    for (i = 1; i <= j; ++i) {
        MARG[i] = PART[i] = -77.77;
        DFS[i] = -9999;
    }
    ntab2 = NCTab / 2 + 1;      

    /*  calculate mm, i.e. maximum of (i * NPi) (i = 1,NDIM) to allocate
        memory. */

    mm = 1;
    m = NDIM1 / 2;
    for (i = 1; i <= m; ++i)
        mm = (mm * (NDIM - i)) / i;
    mm *= NDIM;

    /*  allocate memory for some auxiliary arrays */

    if (!(LSISet = (int   *)calloc(NDIM + 1,sizeof(int)))) { 
        p_err(-2,1);
        goto LSFin;
    }
    LSISetA = NDIM + 1;
    memrq(LSISetA,sizeof(int));

    if (!(LSJSet = (int   *)calloc(NDIM + 1,sizeof(int)))) { 
        p_err(-2,1);
        goto LSFin;
    }
    LSJSetA = NDIM + 1;
    memrq(LSJSetA,sizeof(int));

    if (!(LSdim  = (short *)calloc(NDIM + 1,sizeof(short)))) { 
        p_err(-2,1);
        goto LSFin;
    }
    LSdimA = NDIM + 1;
    memrq(LSdimA,sizeof(short));

    if (!(LSsiz  = (short *)calloc(NDIM + 1,sizeof(short)))) {  
        p_err(-2,1);
        goto LSFin;
    }
    LSsizA = NDIM + 1;
    memrq(LSsizA,sizeof(short));

    if (!(LScord = (short *)calloc(NDIM + 1,sizeof(short)))) {
        p_err(-2,1);
        goto LSFin;
    }
    LScordA = NDIM + 1;
    memrq(LScordA,sizeof(short));

    if (!(LSconf = (short *)calloc(NDIM * iemax + 1,sizeof(short)))) {
        p_err(-2,1);
        goto LSFin;
    }
    LSconfA = NDIM * iemax + 1;
    memrq(LSconfA,sizeof(short));

    if (!(LSip = (short *)calloc(NDIM1 * iemax + 1,sizeof(short)))) {
        p_err(-2,1);
        goto LSFin;
    }
    LSipA = NDIM1 * iemax + 1;
    memrq(LSipA,sizeof(short));

    if (!(LSim = (short *)calloc(NDIM1 * mm + 1,sizeof(short)))) {  
        p_err(-2,1);
        goto LSFin;
    }
    LSimA = NDIM1 * mm + 1;
    memrq(LSimA,sizeof(short));

    if (!(LSx = (float *)calloc(ntab2 + 1,sizeof(float)))) {  
        p_err(-2,1);
        goto LSFin;
    }
    LSxA = ntab2 + 1;
    memrq(LSxA,sizeof(float));

    if (!(LSy = (float *)calloc(ntab2 + 1,sizeof(float)))) {
        p_err(-2,1);
        goto LSFin;
    }
    LSyA = ntab2 + 1;
    memrq(LSyA,sizeof(float));

    /*  change order of dimensions in the contingency table to account
        for the different array handling in Fortran and C. */

    for (i = 1; i <= NDIM; ++i)
        LSdim[i] = CTNCat[NDIM - i];

    DGFR[NDIM] = NCTab - 1;

    logscr_res();
    GSQ[1] = logscr_lr();

    for (m = 1; m <= NDIM1; ++m) {
        np = logscr_conf(m,LSISet,LSJSet);
        logscr_res();
        DGFR[m] = logscr_eval(LSip,np,m,1);
        logscr_fit(np);
        GSQ[m + 1] = logscr_lr();

        for (i = 1; i <= m; ++i) {
            itp = (int)LSip[i];
            np1 = np - 1;
            for (j = 1; j <= np1; ++j)
                LSip[(j - 1) * NDIM1 + i] = LSip[j * NDIM1 + i];
            LSip[(np - 1) * NDIM1 + i] = (int)itp;
        }
        l3 = -m + 1;
        for (j = 1; j <= np; ++j) {
            jn = (j - 1) * NDIM1;
            jn1 = (np - j) * NDIM1;   /* change of order */

            logscr_res();
            df = logscr_eval(LSip,np - 1,m,1);
            logscr_fit(np - 1);
            g21 = logscr_lr();

            DFS[jn1 + m] = DGFR[m] - df;
            PART[jn1 + m] = g21 - GSQ[m + 1];

            if (m == 1) {
                MARG[jn1 + 1] = PART[jn1 + 1];
            }
            else {
                logscr_res();
                df = logscr_eval(LSip,1,m,np);
                logscr_fit(1);
                g22 = logscr_lr();
                l3 += m;
                logscr_res();
                df = logscr_eval(LSim,m,m - 1,l3);
                logscr_fit(m);
                g23 = logscr_lr();
                MARG[jn1 + m] = g23 - g22;
            }
            for (i = 1; i <= m; ++i) {
                itp = (int)LSip[(np - 1) * NDIM1 + i];
                LSip[(np - 1) * NDIM1 + i] = LSip[jn + i];
                LSip[jn + i] = (int)itp;
            }
        }
        DGFR[NDIM] -= DGFR[m];
        GSQ[m] -= GSQ[m + 1];
    }

    /*  print test statistics for all interaction of order ... are zero */

    printf1("Test for zero interactions.\n\n");
    printf1("Order   DF  ");
    prnchar(' ',PMTFmt1 - 9,0);  printf1("Statistic");
    prnchar(' ',PMTFmt1 - 4,0); printf1("Prob \n");
    prnchar('-',11 + 2 * (PMTFmt1 + 1),0);
    for (i = 1; i <= NDIM; ++i) {
        printf1("\n%4d %5d  ",i,DGFR[i]);
        printf1(PMTFmtS,GSQ[i]);
        printf1(PMTFmtS,1.0 - cdchif(GSQ[i],DGFR[i]));
    }

    /*  print test statistics for marginal and partial association */

    printf1("\n\nMarginal and partial association.\n");
    printf1("\nModel");
    if ((mm = 3 * NDIM1) < 5)
        mm = 5;
    prnchar(' ',mm - 5,0);
    printf1("  DF");
    prnchar(' ',PMTFmt1 - 8,0); printf1("Marginal");
    prnchar(' ',PMTFmt1 - 3,0); printf1("Prob");
    prnchar(' ',PMTFmt1 - 6,0); printf1("Partial");
    prnchar(' ',PMTFmt1 - 3,0); printf1("Prob\n");
    prnchar('-',3 * NDIM1 + 4 * (PMTFmt1 + 1) + 3,0);

    for (i = 1; i <= NDIM1; ++i) {
        last = 1;
        m = i;
        while (!(last = logscr_comb(LSISet,NDIM,i,last))) {
            printf1("\n");
            for (j = 1; j <= i; ++j)
                printf1(" %c ",'A' + (char)(LSISet[j] - 1));
            prnchar(' ',mm - 3 * i,0);
            printf1("%4d",DFS[m]); 
            printf1(PMTFmtS,MARG[m]);
            printf1(PMTFmtS,1.0 - cdchif(MARG[m],DFS[m]));
            printf1(PMTFmtS,PART[m]);
            printf1(PMTFmtS,1.0 - cdchif(PART[m],DFS[m]));
            m += NDIM1;
        }
    }
    printf1("\n\n");
    err = 0;

LSFin:
    if (CTFitA > 0) {
        free((char *)CTFit);
        memrq(-CTFitA,sizeof(double));
        CTFitA = 0;    
    }
    if (GSQA > 0) {
        free((char *)GSQ);
        memrq(-GSQA,sizeof(double));
        GSQA = 0;    
    }
    if (PARTA > 0) {
        free((char *)PART);
        memrq(-PARTA,sizeof(double));
        PARTA = 0;    
    }
    if (MARGA > 0) {
        free((char *)MARG);
        memrq(-MARGA,sizeof(double));
        MARGA = 0;    
    }
    if (DGFRA > 0) {
        free((char *)DGFR);
        memrq(-DGFRA,sizeof(int));
        DGFRA = 0;    
    }
    if (DFSA > 0) {
        free((char *)DFS);
        memrq(-DFSA,sizeof(int));
        DFSA = 0;    
    }
    if (LSISetA > 0) {
        free((char *)LSISet);
        memrq(-LSISetA,sizeof(int));
        LSISetA = 0;    
    }
    if (LSJSetA > 0) {
        free((char *)LSJSet);
        memrq(-LSJSetA,sizeof(int));
        LSJSetA = 0;    
    }
    if (LSdimA > 0) {
        free((char *)LSdim);
        memrq(-LSdimA,sizeof(short));
        LSdimA = 0;    
    }
    if (LSsizA > 0) {
        free((char *)LSsiz);
        memrq(-LSsizA,sizeof(short));
        LSsizA = 0;    
    }
    if (LScordA > 0) {
        free((char *)LScord);
        memrq(-LScordA,sizeof(short));
        LScordA = 0;    
    }
    if (LSconfA > 0) {
        free((char *)LSconf);
        memrq(-LSconfA,sizeof(short));
        LSconfA = 0;    
    }
    if (LSipA > 0) {
        free((char *)LSip);
        memrq(-LSipA,sizeof(short));
        LSipA = 0;    
    }
    if (LSimA > 0) {
        free((char *)LSim);
        memrq(-LSimA,sizeof(short));
        LSimA = 0;    
    }
    if (LSxA > 0) {
        free((char *)LSx);
        memrq(-LSxA,sizeof(float));
        LSxA = 0;    
    }
    if (LSyA > 0) {
        free((char *)LSy);
        memrq(-LSyA,sizeof(float));
        LSyA = 0;    
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  logscr_conf    Create actual model configuration.                       */

int logscr_conf(int m, int *iset, int *jset)
{
    register int i;
    int l,js,nm,nm1,np,np1,ilast,jlast;

    ilast = 1;
    np = nm = 0;
    
    while (1) {
        while (1) {
            ilast = logscr_comb(iset,NDIM,m,ilast);
            if (ilast)
                return(np);
            np++;
            np1 = (np - 1) * NDIM1;
            for (i = 1; i <= m; ++i)
                LSip[np1 + i] = (short)iset[i];
            if (m != 1)
                break;
        }
        jlast = 1;
        l = m - 1;

        while (1) {
            jlast = logscr_comb(jset,m,l,jlast);
            if (jlast)
                break;
            nm++;
            nm1 = (nm - 1) * NDIM1;
            for (i = 1; i <= l; ++i) {
                js = jset[i];
                LSim[nm1 + i] = (short)iset[js];
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  logscr_comb                                                             */
/*      Create all combinations of length m of iset[1,...,n].               */
/*      If called with last = 1, the function is initialized. Then, as      */
/*      long as a new combination is found, the function return 0. If no    */
/*      further combination can be found, the function returns 1.           */

int logscr_comb(int *iset, int n, int m, int last) 
{
    register int i,k,l;

    if (last) {
        for (i = 1; i <= m; ++i)
            iset[i] = (int)i;
    }
    else {
        k = m;
        while (k > 0) {
            l = (int)iset[k] + 1;
            if ((l + m - k) <= n) {
                for (i = k; i <= m; ++i)  
                    iset[i] = l++;
                return(last);
            }
            k--;
        }
    }
    if (last)
        last = 0;
    else
        last = 1;
    return(last);
}

/* ------------------------------------------------------------------------ */
/*  logscr_eval                                                             */
/*      Note: config has column dimension NDIM, iar has dimension NDIM1.      */
/*      Return: degrees of freedom.                                         */

int logscr_eval(short *iar, int nc, int nv, int ibeg)
{
    register int i,j,jn;
    int k,kk,l,l1,df;

    df = 0;
    for (j = 1; j <= nc; ++j) {
        jn = (j - 1) * NDIM;
        kk = 1;
        l  = ibeg + j - 1;
        l1 = (l - 1) * NDIM1;

        for (i = 1; i <= nv; ++i) {
            k = iar[l1 + i];
            kk *= (LSdim[k] - 1);
            LSconf[jn + i] = (int)k;
        }
        LSconf[jn + nv + 1] = 0;
        df += kk;
    }
    return(df);
}

/* ------------------------------------------------------------------------ */
/*  logscr_res  Initialize CTFit with mean values.                          */

void logscr_res(void)
{
    register int i;
    double avg;

    avg = (double)NOC / (double)NCTab;
    for (i = 1; i <= NCTab; ++i)  
        CTFit[i] = avg;
}

/* ------------------------------------------------------------------------ */
/*  logscr_lr  Calculate and return likelihood ratio Chi square statistic.  */

double logscr_lr(void)
{
    register int i;
    double tmp,gsq;

    gsq = 0.0;
    for (i = 1; i <= NCTab; ++i) {
        tmp = (double)CTFreq[i];
        if (CTFit[i] > EPSI && tmp > EPSI)  
            gsq += tmp * rlog(tmp / CTFit[i]);
    }
    gsq *= 2.0;
    return(gsq);
}

/* ------------------------------------------------------------------------ */
/*  logscr_fit                                                              */
/*      Fit a log-linear model. The algorithm is adapted from               */
/*      Haberman, Log-Linear Fit for Contingency Tables (AS 51),            */
/*      Applied Statistics 21 (1972), 218 - 225.                            */

void logscr_fit(int ncon)
{
    register int i,j,k,ii,kk;
    int iin,l,n,option,isz;
    double xmax,e;

    for (kk = 1; kk <= MxIter; ++kk) {
        xmax = 0.0;
        for (ii = 1; ii <= ncon; ++ii) {
            iin = (ii - 1) * NDIM;
            option = 1;
            LSsiz[1] = 1;
            n = 0;
            for (k = 1; k <= NDIM1; ++k) {
                if (!(l = (int)LSconf[iin + k]))
                    break;
                LSsiz[k + 1] = LSsiz[k] * (short)LSdim[l];
                n++;
            }
            isz = (int)LSsiz[n + 1];
            for (j = 1; j <= isz; ++j)  
                LSx[j] = LSy[j] = 0.0;
L_Cont1:
            for (k = 1; k <= NDIM; ++k)
                LScord[k] = 0;
            i = 1;
L_Cont2:
            j = 1;
            for (k = 1; k <= n; ++k) {
                l = (int)LSconf[iin + k];
                j += (int)LScord[l] * (int)LSsiz[k];
            }
            if (option) {
                LSx[j] += (float)CTFreq[i];
                LSy[j] += (float)CTFit[i];
            }
            else {
                if (LSy[j] <= 0.0)
                    CTFit[i] = 0.0;
                else
                    CTFit[i] *= (double)(LSx[j] / LSy[j]);
            }
            i++;
            for (k = 1; k <= NDIM; ++k) {
                LScord[k] += 1;
                if ((int)LScord[k] < LSdim[k])
                    goto L_Cont2;
                LScord[k] = 0;
            }
            if (option) {
                option = 0;
                goto L_Cont1;
            }
            for (i = 1; i <= isz; ++i) {
                e = fabs((double)(LSx[i] - LSy[i]));
                if (e > xmax)
                    xmax = e;
            }
        }
        if (xmax < MaxDev)
            return;
    }
}

/* ------------------------------------------------------------------------ */
/*  e_expm    Expand an expression                                          */
/*                                                                          */
/*  expm = expression;                                                      */
/*                                                                          */
/*  Return 0 if successful, -1 if error.                                    */

int e_expm(void)
{
    int err;

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("Expand model description. Current memory: %d bytes.\n\n",MemReq);

    if (parm(CmdBuf + 4,8,1))     /* get parameters */
        goto EEFin;

    printf1("expm=%s\n",PMRHSTR);

    err = 0;
EEFin:
    p_clean();
    return(err);
}

