/****************************************************************************/
/*  t_qrmod                                                                 */
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
#include "t_gdat.h"   
#include "t_freq.h"   
#include "t_ml.h"   
#include "t_gf.h"   
#include "t_var.h"   
#include "t_cdf.h"   
#include "t_alloc.h"   
#include "t_sort.h"   
#include "t_con.h"   
#include "t_min.h"   
#include "t_fnqr.h"   

/*  functions in t_qrmod.c */

int qreg(void);
int prn_qrmod(int mod);
int qr_check(void);
void qr_nparm(void);
int qr_alloc(int typ,int opt);
void qr_sval(void);
void prn_qcoeff(int mod);
void prn_qc1(int k);
void qr_options(int mod);
void qr_resid(int mod);
void qr_resida(int k, int iv);
int qr_svar(int opt);
int qr_depcat(int opt);
int pdmc_init(int opt);
void qr_pdata(void);
void qr_prob1(int mod,int icase,int wave,double *par,double *prob);
void qr_prob2(int mod,int icase,int wave,double *par,double *prob);
void qr_log3(int icase,int wave,double *par,double *prob);
void qr_corr(void);
int rmod(void);

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */

int QRTyp = 0;              /* type of quantal response model               */
int NWave = 0;              /* numbe of waves                               */
int NCV = 0;                /* number of valid cases                        */
int IntFlg = 1;             /* zero if without an intercept                 */

short *PYVar;               /* index of dependent variable in wave ...      */
int PYVarA = 0;             /* if allocated                                 */

int NPX = 0;                /* number of independent variables              */
int NPX1;                   /* NPX + IntFlg                                 */
short **PXVar;              /* indices of independent variables in wave ... */
int PXVarA = 0;             /* if allocated                                 */
int PXNVarA = 0;            /* if allocated                                 */

int NPZ = 0;                /* number of Z variables                        */
short **PZVar;              /* indices of z variables in wave ...           */
int PZVarA = 0;             /* if allocated                                 */
int PZNVarA = 0;            /* if allocated                                 */
int PZNVarAL = 0;           /* if allocated                                 */

int NYCat = 0;              /* number of categories in dependent variables  */
int *YCat;                  /* categories of dep. variables                 */
int *YCatI;                 /* inverse of YCat                              */
int YCatA = 0;              /* if allocated                                 */
int YCatIA = 0;             /* if allocated                                 */
double *PYFreq;             /* weighted frequencies of categories           */
int PYFreqA = 0;            /* if allocated                                 */
/* ------------------------------------------------------------------------ */
int PDMCFlg = 0;            /* if missing value code defined                */
char *PDMCIdx;              /* vector with select flags                     */
short *PDMCTi;              /* number of waves of participation             */
int PDMCIdxA = 0;           /* if allocated                                 */
int PDMCTiA = 0;            /* if allocated                                 */
/* ------------------------------------------------------------------------ */
int QR_INIT = 0;            /* set to 1 for first call of fn...().          */

double *QRProb;             /* temp. storage for probabilities              */
int QRProbA = 0;     
double *QRTmp;              /* temp. storage                                */
int QRTmpA = 0;     

int *L4_ISET;               /* temp storage used in fn_logit4               */
int *L4_IPTR;           
double *L4_XB;          
int L4_ISETA = 0;   
int L4_IPTRA = 0;   
int L4_XBA = 0;     
/* ------------------------------------------------------------------------ */
int RMPN = 0;               /* number of patterns in Rasch model            */
int RMPV = 0;               /* number of variables                          */
short *RMPAT;               /* array with patterns                          */
int RMPATA = 0;             /* allocated                                    */
int *RMPF;                  /* frequencies of patterns                      */
int RMPFA = 0;              /* allocated                                    */

/* ------------------------------------------------------------------------ */
/*  qreg()          Quantal response models.                                */
/*                                                                          */
/*                  qreg(                                                   */
/*                      m=...,      model number, def. 1                    */
/*                                   1 = binary logit                       */
/*                                   2 = binary probit                      */
/*                                   3 = ordinal logit                      */
/*                                   4 = oridnal probit                     */
/*                                   5 = multinomial logit                  */
/*                      nw=...,     number of waves, def. 1                 */
/*                      ni=...,      1 if without intercept, def. 0         */
/*                      maxcat=..., max number of categories, def. 1000     */
/*                      pmin                                                */
/*                      dsv=...,    starting values                         */
/*                      ppar=...,   print parameter                         */
/*                      pcov=...,   cov matrix                              */ 
/*                      tfmt=...,   print format parameter, def. 10.4       */
/*                      mfmt=...,   print format pcov, def. 10.6            */
/*                      con=...,    constraints                             */
/*                  ) = varlist;                                            */
/*                                                                          */
/*                  Return 0 if OK, otherwise -1.                           */

int qreg(void)
{
    int err,nv,n;
           
    err = -1;
    if (check_cmd(2))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    printf1("Quantal response models. Current memory: %d bytes.\n",MemReq);

    set_mldef();                /* set defaults for ML estimation */
    PDMCFlg = 0;                      

    if (parm(CmdBuf + 4,4,1))   /* get parameters */
        goto QRFin;

    if (PMFmtF == 0)
        pmfmt(10,4);

    if (PMMaxCatFlg == 0)
        PMMaxCat = MaxYC;

    NWave = PMNW;

    newline();
    if ((PMM == 5 || PMM == 6) && PMNQ < 1) {
        printf1("Error: multinomial logit and probit models needs nq parameter.\n");
        goto QRFin;
    }
    if ((PMM == 7 || PMM == 8) && NWave < 2) {
        printf1("Error: this model needs at least two waves.\n");
        goto QRFin;
    }
    if (PMM == 6 && PMNQ > DMVMax) {
        printf1("Error: max number of categories is %d.\n",DMVMax);
        goto QRFin;
    }
    if (PMM == 8 && NWave > DMVMax) {
        printf1("Error: max number of waves is %d.\n",DMVMax);
        goto QRFin;
    }
    nv = check_nvar(1);         /* check variables */
    if (nv == 0)                /* PMNV is number of variables */       
        goto QRFin;        

    if (prn_qrmod(PMM))         /* print type of model */
        goto QRFin;
    newline();
       
    prn_nwvar(1);               /* print variables */
          
    if (NPZ > 0) {
        if (PMM != 5 && PMM != 6) {
            printf1("Error: z variables only with m=5 or m=6.\n");
            goto QRFin;
        }
    }
    if (PMM == 7)   /* always without intercept */
        PMNI = 1;

    NPX1 = NPX;
    if (PMNI) {
        printf1("Model without intercept.\n");
        PMNI = 1;
        IntFlg = 0;
    }
    else {
        PMNI = 0;
        IntFlg = 1;
        NPX1++;
        if (NPX + NPZ == 0)  
            printf1("Model without independent variables.\n");
    }
    if (NPX1 + NPZ == 0) {
        printf1("Error: no model parameters.\n");
        goto QRFin;
    }
    if (qr_svar(1))         /* save variables */
        goto QRFin;

    NCV = pdmc_init(1);     /* init missing value flags */
    if (NCV < 2) {
        if (NCV >= 0)
            printf1("Error: insufficient number of cases.\n");
        goto QRFin;
    }
    newline();

    if (qr_depcat(1))       /* distribution of categories */
        goto QRFin;

    if (NPZ > 0 && NYCat != PMNQ) {
        printf1("Error: nq=%d inconsistent with data.\n",PMNQ);
        goto QRFin;
    }
    n = 0;
    if (WIVar >= 0) {
        prn_cwt();          /* case weight information */
        newline();
    }
    if (qr_check())         /* check model requirements */
        goto QRFin;

    qr_nparm();             /* set NParm = number of parameters */

    set_mlopt();            /* adjust ML options */         

    if (ml_init(1,0,0))     /* init ML estimation */
        goto QRFin;

    if (QRTyp == QRPROB3 || QRTyp == QRPROB4) {
        if (PMNHP < 2)
            PMNHP = 2;
        else if (PMNHP > 20)
            PMNHP = 20;
        else if (PMNHP == 11)
            PMNHP = 10;

        if (PMOPT < 1 || PMOPT > 3)
            PMOPT = 1;

        printf1("Control of integration (nhp): %d",PMNHP);
        if (PMNHP >= 12)
            printf1(" [eps=%lg]",PMEPS);
        printf1("\nType of parameterization: %d\n",PMOPT);
    }

    if (get_dsv(NParm,Par,0,ParLB,ParUB,1))   /* try to get starting values */
        goto QRFin;
  
    if (DSVFlg == 0)
        qr_sval();          /* default starting values */

    prot_init(2,QRTyp,0);   /* init protocol file */

    if (PMNConS > 0) {      /* process constraints */
        newline();
        if (con_proc(CmdBuf))
            goto QRFin;
    }
    if (NCV < NParm1) {
        printf1("Error: insufficient number of cases for %d parameters.\n",NParm1);
        goto QRFin;
    }

    if (qr_alloc(QRTyp,1))  /* additional storage */
        goto QRFin;

    QR_INIT = 1;            /* flag first call of fn_logit(), fn_probit() */

    if (ffmin(QRTyp,0,1,1))   /* minimization */
        goto QRFin;
           
    prn_mlres(0);           /* print results of min algorithm */

    if (QRTyp == QRLOG4)  
        printf1("Number of cases used for likelihood: %d\n\n",NFLUsed);

    prn_qcoeff(QRTyp);      /* print estimated parameters etc. */   
           
    if (LConv >= 0)  
        prn_ml1res(0);      /* additional results */
      
    if (PMF1Def)            /* print data to output file */
        qr_pdata();       

    qr_options(QRTyp);      /* additional options */

    if ((QRTyp == QRPROB3 || QRTyp == QRPROB4) && PMPPFDef)  /* write correlation matrix */
        qr_corr();       

    err = 0;

QRFin:
    QR_INIT = 0; 
    con_free();
    ml_init(0,0,0);      
    qr_alloc(0,0);
    qr_depcat(0);
    pdmc_init(0);
    qr_svar(0);
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_qrmod(mod)      Print name of quantal response model and set QRTyp. */
/*                                                                          */
/*              1 : QRLOG1      binary logit                                */
/*              2 : QRPROB1     binary probit                               */
/*              3 : QRLOG2      ordinal logit                               */
/*              4 : QRPROB2     ordinal probit                              */
/*              5 : QRLOG3      multinomial logit                           */
/*              6 : QRPROB3     multinomial probit                          */
/*              7 : QRLOG4      conditional logit                           */
/*              8 : QRPROB4     simultaneous binary probit                  */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int prn_qrmod(int mod)
{
    printf1("Model: ");
    switch (mod) {
        case  1:        printf1("binary logit.\n");
                        QRTyp = QRLOG1;
                        break;
        case  2:        printf1("binary probit.\n");
                        QRTyp = QRPROB1;
                        break;
        case  3:        printf1("ordinal logit.\n");
                        PMNI = 1;
                        QRTyp = QRLOG2;
                        break;
        case  4:        printf1("ordinal probit.\n");
                        PMNI = 1;
                        QRTyp = QRPROB2;
                        break;
        case  5:        printf1("multinomial logit.\n");
                        QRTyp = QRLOG3;
                        break;
        case  6:        printf1("multivariate probit.\n");
                        QRTyp = QRPROB3;
                        break;
        case  7:        printf1("conditional logit.\n");
                        QRTyp = QRLOG4;
                        break;
        case  8:        printf1("simultaneous binary probit.\n");
                        QRTyp = QRPROB4;
                        break;

        default:        printf1("%d -- not defined.\n",mod);
                        return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  qr_check()      check model requirements, model is QRTyp                */
/*                  and set PMDOPT according to derivatives.                */
/*  Return 0 if OK, -1 if error.                                            */

int qr_check(void)
{
    if (QRTyp == QRLOG1   || QRTyp == QRPROB1 || QRTyp == QRLOG4 || QRTyp == QRPROB4) {
        if (NYCat != 2) {
            printf1("Error: model requires two categories in dependent variable.\n");
            return(-1);
        }
    }
    else if (QRTyp == QRLOG2 || QRTyp == QRPROB2) {
        if (NYCat < 3) {
            printf1("Error: model requires three or more categories in dependent variable.\n");
            return(-1);
        }
    }
    if (QRTyp == QRPROB3 || QRTyp == QRPROB4) { /* only function values */
        PMDOPT = 0;
        if (MINA < 7 || MINA > 8)
            MINA = 8;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  qr_nparm()      set number of model parameters for QRTyp                */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

void qr_nparm(void)
{
    NParm = 0;       

    switch (QRTyp) {

        case QRLOG1: 
        case QRLOG4:
        case QRPROB1:   NParm = NPX1;
                        break;
        case QRLOG2:
        case QRPROB2:   NParm = NPX + NYCat - 1;
                        break;

        case QRLOG3:    NParm = (NYCat - 1) * NPX1 + NPZ;
                        break;

        case QRPROB3:   NParm = (NYCat - 1) * NPX1 + NPZ + (NYCat * (NYCat - 1)) / 2;
                        break;

        case QRPROB4:   NParm = NWave * NPX1 + (NWave * (NWave - 1)) / 2;
                        break;

        default:        printf1("Error: model %d not available.\n",QRTyp);
                        return;       
    }
}

/* ------------------------------------------------------------------------ */
/*  qr_alloc(typ,opt)   If opt != 0 allocate additional storage, depending  */
/*                      on typ. Otherwise free previously allocated         */
/*                      memory.                                             */
/*                      Return 0 if OK, otherwise -1.                       */

int qr_alloc(int typ,int opt)
{
    int err,n;

    if (opt == 0) {
        err = 0;
        goto QRALFin;
    }
    err = -1;

    n = NYCat;
    if (n < NWave)
        n = NWave;

    if (!(QRTmp = (double *)calloc(n + 1,sizeof(double))))  
        goto QRALFin;  
    memrq(n + 1,sizeof(double));
    QRTmpA = n + 1;

    if (typ == QRPROB3 || typ == QRPROB4)
        n *= n;                

    if (!(QRProb = (double *)calloc(n + 1,sizeof(double))))  
        goto QRALFin;  
    memrq(n + 1,sizeof(double));
    QRProbA = n + 1;

    if (typ == QRLOG4) {       /* conditional logit */

        n = NWave + 1;

        if (!(L4_ISET = (int *)calloc(n,sizeof(int))))  
            goto QRALFin;  
        memrq(n,sizeof(int));
        L4_ISETA = n;
       
        if (!(L4_IPTR = (int *)calloc(n,sizeof(int))))  
            goto QRALFin;  
        memrq(n,sizeof(int));
        L4_IPTRA = n;
         
        if (!(L4_XB = (double *)calloc(n,sizeof(double))))  
            goto QRALFin;  
        memrq(n,sizeof(double));
        L4_XBA = n;
    }
    return(0);

QRALFin:
    if (err)
        p_err(-2,1);

    if (QRProbA > 0) {
        free((char *)QRProb);
        memrq(-QRProbA,sizeof(double));
        QRProbA = 0;        
    }
    if (QRTmpA > 0) {
        free((char *)QRTmp);
        memrq(-QRTmpA,sizeof(double));
        QRTmpA = 0;        
    }
    if (L4_ISETA > 0) {
        free((char *)L4_ISET);
        memrq(-L4_ISETA,sizeof(int));
        L4_ISETA = 0;        
    }
    if (L4_IPTRA > 0) {
        free((char *)L4_IPTR);
        memrq(-L4_IPTRA,sizeof(int));
        L4_IPTRA = 0;        
    }
    if (L4_XBA > 0) {
        free((char *)L4_XB);
        memrq(-L4_XBA,sizeof(double));
        L4_XBA = 0;        
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  qr_sval()       Calculate  default starting values for model QRTyp      */

void qr_sval(void) 
{
    register int j;
    double tmp,tmp1;

    for (j = 1; j <= NParm; ++j)
        Par[j] = 0.0;

    if (QRTyp == QRLOG2 || QRTyp == QRPROB2) {
        tmp = tmp1 = 0.0;
        for (j = 0; j < NYCat; ++j)  
            tmp1 += PYFreq[j];
        for (j = 1; j < NYCat; ++j) {
            tmp += PYFreq[j - 1];
            if (tmp > 0.0 && tmp1 > tmp)
                Par[j] = rlog((tmp1 - tmp) / tmp);           
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  Function: prn_qcoeff(mod)                                               */
/*      Print estimated coefficients for quantal response model mod.        */

void prn_qcoeff(int mod)
{
    register int j,k,l,li;
    int len,nvlen,wave;
    double tsum = 0.0;

    nvlen = PMNVLEN;
    if (nvlen < 10)
        nvlen = 10;

    len = 23;
    printf1("Idx Cat Term   ");
    printf1("Variable ");
    prnchar(' ',nvlen - 8,0);
    prnchar(' ',PMTFmt1 - 5,0); printf1("Coeff"); 
    prnchar(' ',PMTFmt1 - 4,0); printf1("Error"); 
    prnchar(' ',PMTFmt1 - 6,0); printf1("C/Error  Signif");
       
    k = 1;
    wave = 1;

NXT:
    if (NPX1) {
        for (j = 1; j < NYCat; ++j) {

            printf1("\n");
            prnchar('-',len + 3 * (PMTFmt1 + 1) + nvlen,0);

            for (l = 1; l <= NParm; ++l) {

                if (mod == QRLOG2 || mod == QRPROB2) {
                    printf1("\n%3d   - ",k);
                    if (l < NYCat) {
                        printf1(" I %-3d Alpha %-3d ",l,l); 
                        prnchar(' ',nvlen - 9,0);
                    }
                    else {
                        li = l - NYCat;
                        printf1(" X     %s ",VName[PXVar[li][0]]);
                        prnchar(' ',nvlen - strlen(VName[PXVar[li][0]]),0);
                    }
                }
                else {
                    printf1("\n%3d %3d ",k,YCat[j]);
                    if (IntFlg && l == 1) {
                        if (mod == QRPROB4)
                            printf1(" W%2d",wave);
                        else
                            printf1(" I  ");
                        printf1("   Intercept ");
                        prnchar(' ',nvlen - 9,0);
                    }
                    else {
                        li = l - IntFlg - 1;
                        if (mod == QRPROB4)
                            printf1(" W%2d",wave);
                        else
                            printf1(" X  ");
                        printf1("   %s ",VName[PXVar[li][0]]);
                        prnchar(' ',nvlen - strlen(VName[PXVar[li][0]]),0);
                    }
                }
                prn_qc1(k);
                k++;
                if (l >= NPX1 && mod != QRLOG2 && mod != QRPROB2)    
                    break;
            }
            if (mod == QRLOG2 || mod == QRPROB2)
                break;
        }
    }
   
    if (NPZ && (mod == QRLOG3 || mod == QRPROB3)) {     /* generic variables */ 
        printf1("\n");
        prnchar('-',len + 3 * (PMTFmt1 + 1) + nvlen,0);

        for (j = 0; j < NPZ; ++j) {
            printf1("\n%3d   - ",k);
            printf1(" Z%-4d %s ",j + 1,VName[PZVar[j][0]]);
            prnchar(' ',nvlen - strlen(VName[PZVar[j][0]]),0);

            prn_qc1(k);
            k++;
        }
    }
    if (mod == QRPROB4 && wave++ < NWave)
        goto NXT;

    if (mod == QRPROB3 || mod == QRPROB4) {
        printf1("\n");
        prnchar('-',len + 3 * (PMTFmt1 + 1) + nvlen,0);
        for (j = 2; j <= NYCat; ++j) {
            for (l = 1; l < j; ++l) {
                printf1("\n%3d   - ",k);
                printf1(" S     Sigma%2d,%2d ",j,l);
                prnchar(' ',nvlen - 10,0);
                prn_qc1(k);
                k++;
            }
        }
    }
    newline();
    newline();
}

void prn_qc1(int k)     /* print coefficients */
{
    double tmp,tmp1,tmp2;

    printf1(PMTFmtS,Par[k]);
    if (LConv >= 0 && Diag[k] >= 0.0) {
        if (Diag[k] > EPSI1)
            tmp = sqrt(Diag[k]);
        else
            tmp = 0.0;

        printf1(PMTFmtS,tmp);        
    }
    else {
        tmp = 0.0;
        prnchar(' ',PMTFmt1 - 3,0);
        printf1("--- ");
    }
    if (LConv >= 0 && tmp > 0.0) {
        tmp1 = Par[k] / tmp;
        printf1(PMTFmtS,tmp1);          
        tmp2 = 2.0 * cdnf(fabs(tmp1)) - 1.0;
        printf1(" %6.4lf",tmp2);          
    }
    else {
        prnchar(' ',PMTFmt1 - 3,0);
        printf1("---     ---");         
    }
}

/* ------------------------------------------------------------------------ */
/*  qr_options(mod)     process additional options for model mod.           */

void qr_options(int mod)
{
    int i;

    for (i = 0; i < PMResN; ++i) {
        switch (PMRes[i]) {
            case  1:                        /* stand. residuals */
                if (NWave != 1)
                    break;
                if (mod == QRLOG1 || mod == QRPROB1)
                    qr_resid(mod);
                break;

            default:
                break;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  qr_resid()                                                              */
/*      Print standardized coefficients, only for binary and multinomial    */
/*      logit models with cross-sectional data (NWave = 1).                 */
       
void qr_resid(int mod)
{
    register int j,k,l,i,len,nvlen;

    nvlen = PMNVLEN;
    if (nvlen < 9)
        nvlen = 9;
        
    printf1("\nStandardized coefficients.\n\n");

    len = 15;
    printf1("Idx Cat Term   Variable ");
    prnchar(' ',nvlen - 8,0);

    prnchar(' ',PMTFmt1 - 10,0); printf1("     Coeff"); 
    prnchar(' ',PMTFmt1 -  5,0); printf1("Exp(C)"); 
    prnchar(' ',PMTFmt1 -  8,0); printf1("Exp(C*SD)");
    prnchar(' ',PMTFmt1 -  7,0); printf1("Std.Dev.");
    k = 1;

    if (NPX1) {
        for (j = 1; j < NYCat; ++j) {

            newline();           
            prnchar('-',len + 4 * (PMTFmt1 + 1) + nvlen,0);

            for (l = 1; l <= NParm; ++l) {

                printf1("\n%3d %3d ",k,YCat[j]);
                if (IntFlg && l == 1) {
                    printf1(" I     Intercept ");
                    prnchar(' ',nvlen - 9,0);
                    qr_resida(k,-1);
                }
                else {
                    i = l - IntFlg - 1;
                    printf1(" X     %s ",VName[PXVar[i][0]]);           
                    prnchar(' ',nvlen - strlen(VName[PXVar[i][0]]),0);
                    qr_resida(k,PXVar[i][0]);
                }
                k++;
                if (l >= NPX1)
                    break;
            }
        }
    }
    /***
    if (NPZ && mod == QRLOG3) {        generic variables   
        printf1("\n");
        prnchar('-',len + 3 * (PMTFmt1 + 1) + nvlen,0);
        kk = k;
        for (j = 0; j < NPZ; ++j) {
            for (jj = 0; jj < PZN[j]; ++jj) {
                i = PZVar[j][jj];
                printf1("\n%3d %3d ",kk++,YCat[jj]);
                printf1(" Z%-4d %s ",PZNum[j],VName[i]);
                prnchar(' ',nvlen - strlen(VName[i]),0);
                qr_resida(k,i);
            }
            k++;
        }
    }
    **/
    newline();
}

void qr_resida(int k, int iv)     /* print standardized coefficients */
{
    register int i;
    double wt,sdev,sum,ws,tmp;

    sdev = 0.0;

    if (iv >= 0) {
        wt = 1.0;  
        sum = ws = 0.0;

        for (i = 0; i < NOC; ++i) {
                     
            if (PDMCFlg && PDMCTi[i] == 0)  
                continue;

            if (!PDMCFlg || PDMCIdx[i * NWave]) {

                if (WIVar >= 0)
                    wt = get_data(WIVar,i) * WNorm;
                ws += wt;
                tmp = get_data(iv,i);
                sum += tmp * wt;
                sdev += tmp * tmp * wt;
            }
            if (ws > 0.0) {
                tmp = sum / ws;
                sdev -= tmp * tmp * ws;
                if (sdev > 0.0 && ws > 1.0)
                    sdev = sqrt(sdev / (ws - 1.0));
                else
                    sdev = 0.0;
            }
            else
                sdev = 0.0;
        }
    }
    printf1(PMTFmtS,Par[k]);
    printf1(PMTFmtS,rexp(Par[k]));
    printf1(PMTFmtS,rexp(Par[k] * sdev));
    printf1(PMTFmtS,sdev);
}

   
/* ------------------------------------------------------------------------ */
/*  qr_svar(opt)    If opt = 1, set var structure for qr models.            */
/*                  if opt = 0, free previously allocated memory.           */
/*                                                                          */
/*                  NWave               number of waves (= PMNW)            */
/*                  short PYVar[w]      indices of dependent variables      */
/*                                      w = 0,...,NWave - 1.                */
/*                  short PXVar[i][w]   indices of independent variables    */  
/*                                      i = 0,...,NPX-1, w=0,...,NWave-1    */
/*                                                                          */
/*                  short PZVar[i][w]   indices of z variables              */  
/*                                      i = 0,...,NPZ-1                     */
/*                                      w = 0,...,NWave - 1                 */
/*                                      j = 0,...,PMNQ - 1                  */
/*                                      PZVar[i][w * PMNQ + j]              */
/*                                      later: NYCat = PMNQ                 */
/*                                                                          */
/*                  Return 0 if successful, otherwise -1.                   */

int qr_svar(int opt)
{
    register int i,j,k,kk;
    int err;

    err = -1;
    if (opt == 0) {
        err = 0;
        goto QRSVFree;
    }
    if (!(PYVar = (short *)calloc(NWave,sizeof(short))))  
        goto QRSVFree;
    memrq(NWave,sizeof(short));
    PYVarA = NWave;

    for (i = 0; i < NWave; ++i)     /* save indices of dep variable */
        PYVar[i] = PMVIdx[i];

    if (NPX > 0) {                  /* if there are independent variables */

        if (!(PXVar = (short **)calloc(NPX,sizeof(short *))))  
            goto QRSVFree;
        memrq(NPX,sizeof(short *));
        PXVarA = NPX;

        for (i = 0; i < NPX; ++i) {
            if (!(PXVar[i] = (short *)calloc(NWave,sizeof(short))))  
                goto QRSVFree;
            memrq(NWave,sizeof(short));
            PXNVarA++;           
        }
        k = NWave;
        for (i = 0; i < NPX; ++i) {
            for (j = 0; j < NWave; ++j)
                PXVar[i][j] = PMVIdx[k++];
        }
    }

    if (NPZ > 0) {                  /* if there are z variables */

        if (!(PZVar = (short **)calloc(NPZ,sizeof(short *))))  
            goto QRSVFree;
        memrq(NPZ,sizeof(short *));
        PZVarA = NPZ;
        PZNVarAL = NWave * PMNQ;
        for (i = 0; i < NPZ; ++i) {
            if (!(PZVar[i] = (short *)calloc(PZNVarAL,sizeof(short))))  
                goto QRSVFree;
            memrq(PZNVarAL,sizeof(short));
            PZNVarA++;           
        }
        kk = 0;
        for (i = 0; i < NPZ; ++i) {
            for (k = 0; k < NWave; ++k) {
                for (j = 0; j < PMNQ; ++j) 
                    PZVar[i][k * PMNQ + j] = PMZIdx[kk++];
            }
        }
    }
    return(0);

QRSVFree:
    if (PYVarA) {
        free((char *)PYVar);
        memrq(-NWave,sizeof(short));
        PYVarA = 0;       
    }
    if (PXVarA) {
        for (i = 0; i < PXNVarA; ++i) {
            free((char *)PXVar[i]);
            memrq(-NWave,sizeof(short));
        }
        PXNVarA = 0;       
        free((char *)PXVar);
        memrq(-PXVarA,sizeof(short *));
        PXVarA = 0;       
    }
    if (PZVarA) {
        for (i = 0; i < PZNVarA; ++i) {
            free((char *)PZVar[i]);
            memrq(-PZNVarAL,sizeof(short));
        }
        PZNVarAL = PZNVarA = 0;       
        free((char *)PZVar);
        memrq(-PZVarA,sizeof(short *));
        PZVarA = 0;       
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  Function: qr_depcat(opt)                                                */
/*                                                                          */
/*      If opt = 1 calculate and print categories of dependent variable.    */
/*      Also print weights of the categories.                               */
/*      If opt = 0 free previously allocated memory.                        */
/*                                                                          */
/*      NYCat    =   Number of categories                                   */
/*      YCat[i]  =   Array with values of categories                        */
/*                                                                          */
/*      The categories are mapped to internal cat numbers: 0,...,NYCat-1.   */
/*      Mapping is through: YCatI[j], so that: if j is the value of a       */
/*      category (in the input data), then YCatI[j] is the internal cat     */
/*      number. Therefore: If YCat[i] = j, then YCatI[j] = i, i = 0,NYCat-1 */
/*                                                                          */
/*      Generally, the access is as follows:                                */
/*                                                                          */
/*      There are NWave dependent variables with associated variable        */
/*      numbers: PYVar[l], l = 0,NWave-1. Then, for any l (0,NWave-1):      */
/*      PYVar[l] is the internal variable number associated with            */
/*      the l.th dependent variable. The value for the record icase is      */
/*      found by val = get_data(PYVar[l],icase).                            */
/*                                                                          */
/*      Then, the associated internal cat number is found by                */
/*      cat = YCatI[val].                                                   */
/*                                                                          */
/*      The global array PYFreq is build with the weighted frequencies      */
/*      of the categories of the dependent variable, pooled over waves.     */
/*                                                                          */
/*      Note: If panel missing codes are defined, these observations        */
/*      are excluded from the frequency distribution.                       */
/*                                                                          */
/*      Return 0 if successful, -1 if insuff memory, -2 other error.        */

int qr_depcat(int opt)
{
    register int i,j,k,l,kk;
    int err,m,yca,mycat,yfreqa,yptra,jfya,wyfreqa;
    int *yc,*yfreq,*yptr,*jfy;
    double tmp,wt,*wyfreq[MaxW],d1,d2;

    if (opt == 0) {          /* just free memory */
        err = 0;
        goto QRDEPC1Fin;
    }
    wyfreqa = yca = yfreqa = yptra = jfya = 0;
    err = -1;

    printf1("Categories of dependent variable.\n");
    printf1("Maximum number of categories: %d\n\n",PMMaxCat);

    if (!(yc = (int *)calloc (PMMaxCat + 1,sizeof(int))))  
        goto QRDEPCFin;
    memrq(PMMaxCat + 1,sizeof(int));
    yca = 1;

    if (!(yfreq = (int *)calloc (PMMaxCat + 1,sizeof(int))))  
        goto QRDEPCFin;
    memrq(PMMaxCat + 1,sizeof(int));
    yfreqa = 1;

    if (!(yptr = (int *)calloc (PMMaxCat + 1,sizeof(int))))  
        goto QRDEPCFin;
    memrq(PMMaxCat + 1,sizeof(int));
    yptra = 1;

    if (!(jfy = (int *)calloc (NWave * NOC + 1,sizeof(int))))  
        goto QRDEPCFin;
    memrq(NWave * NOC + 1,sizeof(int));
    jfya = 1;

    k = 0;
    for (j = 0; j < NWave; ++j) {
        l = PYVar[j];
        for (i = 0; i < NOC; ++i) { 

            if (!PDMCFlg || PDMCIdx[i * NWave + j]) {
                m = (int)get_data(l,i);
                if (m < 0) {
                    printf1("Error: dependent variable in wave %d has negative categories.\n",j + 1);
                    err = -2;
                    goto QRDEPCFin;
                }
                jfy[++k] = m;
            }
        }
    }

    NYCat = cfreq(1,k,jfy,PMMaxCat,yc,yfreq,yptr,0,&d1,&d2);
    
    if (NYCat <= 0) {
        printf1("Can't create distribution of categories of dependent variables.\n");
        printf1("Possibly too many categories, current maximum is %d.\n",PMMaxCat);
        err = -2;
        goto QRDEPCFin;
    }
    if (NYCat < 2) {
        printf1("Error: dependent variable has less than two categories.\n");
        err = -2;
        goto QRDEPCFin;
    }

    /*  Now we have NYCat categories in the array yc. Build another array
        of categories, YCat, where the categories are sorted in ascending
        order: YCat[j], j = 0, NYCat - 1. */        
                                                    
    if (!(YCat = (int *)calloc (NYCat,sizeof(int))))   
        goto QRDEPCFin; 
    memrq(NYCat,sizeof(int));
    YCatA = NYCat;

    if (!(PYFreq = (double *)calloc (NYCat,sizeof(double))))  
        goto QRDEPCFin; 
    memrq(NYCat,sizeof(double));
    PYFreqA = NYCat;

    mycat = 0;
    for (i = 1; i <= NYCat; ++i) {
        j = yptr[i];
        YCat[i - 1] = m = yc[j];
        if (mycat < m)
            mycat = m;
    }

    /*  build inverse array YCatI to map the categories to internal
        category numbers: 0,..., NYCat - 1. */

    if (!(YCatI = (int *)calloc (mycat + 1,sizeof(int))))  
        goto QRDEPCFin; 
    memrq(mycat + 1,sizeof(int));
    YCatIA = mycat + 1;

    for (i = 0; i < NYCat; ++i) {
        j = YCat[i];
        YCatI[j] = i;
    }

    /*  Build wyfreq[i], weighted frequencies, separately for each wave i */
    /*  PYFreq contains the same but pooled over all waves. */

    k = 0;
    wt = 1.0;            
    k = 0;
    for (i = 0; i < NWave; ++i) {

        if (!(wyfreq[i] = (double *)calloc (NYCat + 1,sizeof(double))))  
            goto QRDEPCFin;
        memrq(NYCat + 1,sizeof(double));
        wyfreqa++;

        for (j = 0; j < NOC; ++j) { 

            if (!PDMCFlg || PDMCIdx[j * NWave + i]) {

                kk = jfy[++k];
                l = YCatI[kk];

                if (WIVar >= 0)                      /* get weights */
                    wt = get_data(WIVar,j) * WNorm;
       
                wyfreq[i][l] += wt;
                PYFreq[l] += wt;
            }
        }
    }
    printf1("Index      ");
    for (j = 0; j < NYCat; ++j)  
        printf1("  %8d",j);
       
    printf1("   (Weighted)\nCategory   ");

    for (j = 0; j < NYCat; ++j)  
        printf1("  %8d",YCat[j]);
       
    printf1("  Observations\n");
                 
    prnchar('-',25 + NYCat * 10,0);
    for (i = 0; i < NWave; ++i) {
        wt = 0.0;
        for (j = 0; j < NYCat; ++j)  
            wt += wyfreq[i][j];


        printf1("\nWave %-3d N ",i + 1);
        for (j = 0; j < NYCat; ++j)  
            printf1("%10.2lf",wyfreq[i][j]);
             
        printf1(" %13.2lf\n         Pct ",wt);
        for (j = 0; j < NYCat; ++j) {
            if (wt > 0.0)
                tmp = 100.0 * wyfreq[i][j] / wt;
            else
                tmp = 0.0;

            printf1("%8.2lf  ",tmp);
        }
    }
    newline();                     
    err = 0;
    newline();

QRDEPCFin:
    if (yca) {
        free((char *)yc);
        memrq(-PMMaxCat - 1,sizeof(int));
    }
    if (yfreqa) {
        free((char *)yfreq);
        memrq(-PMMaxCat - 1,sizeof(int));
    }
    if (yptra) {
        free((char *)yptr);
        memrq(-PMMaxCat - 1,sizeof(int));
    }
    if (jfya) {
        free((char *)jfy);
        memrq(-NWave * NOC - 1,sizeof(int));
    }
    for (i = 0; i < wyfreqa; ++i) {
        free((char *)wyfreq[i]);
        memrq(-NYCat - 1,sizeof(double));
    }

QRDEPC1Fin:
    if (err || opt == 0) {
        if (YCatA) {
            free((char *)YCat);
            memrq(-YCatA,sizeof(int));
            YCatA = 0;
        }
        if (YCatIA) {
            free((char *)YCatI);
            memrq(-YCatIA,sizeof(int));
            YCatIA = 0;
        }
        if (PYFreqA) {
            free((char *)PYFreq);
            memrq(-PYFreqA,sizeof(double));
            PYFreqA = 0;
        }
    }
    if (err == -1)
        p_err(-2,1);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  pdmc_init   If opt = 1 init data structure for panel miss. values.      */
/*              If opt = 0 free previously allocated memory.                */
/*                                                                          */
/*              PDMCIdx of length NOC * NWave is set to 1 if an observation */
/*                      is available                                        */
/*              PDMCTi  of length NOC is the number of waves an individual  */
/*                      has participated.                                   */
/*                                                                          */
/*              A case is missing if the dependent variable has a negative  */
/*              value.                                                      */
/*              PMPMin = min number of waves for participation (pmin)       */
/*              PDMCTi is only set to a nonzero value if an individual has  */
/*              participated at least PMPMin waves with valid values.       */
/*              Otherwise PDMCTi is zero and also PDMCIdx is zero for       */
/*              all waves.                                                  */
/*                                                                          */
/*              If there is at least one missing value then PMMCFlg = 1,    */  
/*              otherwise PDMCFlg = 0.                                      */
/*                                                                          */
/*              Return  number of valid cases.                              */
/*                     -1  if insufficient memory                           */

int pdmc_init(int opt)
{
    register int i,j;
    int err,n,nn,ni;

    if (opt == 0) {
        err = 0;
        goto PDMCFin;
    }
    printf1("Checking available data (pmin=%d)\n",PMPMin);

    err = -1;
    PDMCFlg = 0;         

    if (!(PDMCIdx = (char *) calloc(NOC * NWave + 1,sizeof(char)))) {
        p_err(-2,1);
        goto PDMCFin;
    }
    PDMCIdxA = NOC * NWave + 1;
    memrq(PDMCIdxA,sizeof(char));

    if (!(PDMCTi = (short *) calloc(NOC + 1,sizeof(short)))) {
        p_err(-2,1);
        goto PDMCFin;
    }
    PDMCTiA = NOC + 1;
    memrq(PDMCTiA,sizeof(short));

    nn = ni = 0;
    for (i = 0; i < NOC; ++i) {
        PDMCTi[i] = 0;
        n = 0;
        for (j = 0; j < NWave; ++j) {
            if (get_data(PMVIdx[j],i) >= 0) {
                PDMCIdx[i * NWave + j] = 1;
                n++;
            }
        }
        if (n >= PMPMin) {
            PDMCTi[i] = n;
            ni++;
            nn += n;
        }
        else
            PDMCFlg = 1;
    }
    printf1("Number of cases with valid data: %d\n",ni);
    if (NWave > 1)
        printf1("Number of valid observations: %d\n",nn);

    return(ni);

PDMCFin:
    if (err || opt == 0) {
        if (PDMCIdxA) {
            free(PDMCIdx);
            memrq(-PDMCIdxA,sizeof(char));
            PDMCIdxA = 0;                 
        } 
        if (PDMCTiA) {
            free((char *)PDMCTi);
            memrq(-PDMCTiA,sizeof(short));
            PDMCTiA = 0;                 
        } 
        PDMCFlg = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  qr_pdata()      Print data to PMF1d.                                    */
/*                  Use PMFmtS.                                             */
/*                                                                          */
/*                  Return 0 if successful, otherwise -1.                   */

void qr_pdata(void)
{
    register int i,j,k,l;
    int r,n,nrec,pflag;
    double wt,prob;

    if (PMF1Def == 0)
        return;

    nrec = 0;
    for (i = 0; i < NOC; ++i) { 
 
        if (PDMCFlg && PDMCTi[i] == 0)  
            continue;

        for (j = 0; j < NWave; ++j) {

            if (!PDMCFlg || PDMCIdx[i * NWave + j]) {

                fprintf(PMF1d,"%6d %2d ",i + 1,j + 1);

                fprintf(PMF1d,PMFmtS,get_data(PYVar[j],i));
                for (k = 0; k < NPX; ++k)
                    fprintf(PMF1d,PMFmtS,get_data(PXVar[k][j],i));
   
                for (k = 0; k < NPZ; ++k) {
                    for (l = 0; l < NYCat; ++l) {
                        n = PZVar[k][j * NYCat + l];
                        fprintf(PMF1d,PMFmtS,get_data(n,i));
                    }
                }
   
                if (WIVar >= 0) {   
                    wt = get_data(WIVar,i) * WNorm;
                    fprintf(PMF1d,PMFmtS,wt);
                }
                pflag = 1;

                switch (QRTyp) {
                    case QRLOG1:
                    case QRPROB1:   qr_prob1(QRTyp,i,j,Par,QRTmp);
                                    break;
                    case QRLOG2:
                    case QRPROB2:   qr_prob2(QRTyp,i,j,Par,QRTmp);
                                    break;

                    case QRLOG3:    qr_log3(i,j,Par,QRTmp);
                                    break;

                    case QRPROB3:   for (k = -1; k < NYCat; ++k) {
                                        r = fn_prob3(i,j,k,Par,&prob);
                                        if (r)
                                            prob = -1.0;
                                        QRTmp[k + 1] = prob;
                                    }
                                    break;
                    default:        pflag = 0;
                                    break;
                }
                if (pflag) {
                    for (k = 0; k <= NYCat; ++k)
                        fprintf(PMF1d,PMFmtS,QRTmp[k]);
                }
                fprintf(PMF1d,"\n");
                nrec++;
            }
        }
    }
    printf1("Data (%d records) written to: %s\n",nrec,PMF1dName);

    if (PMTDAFDef) {
        fprintf(PMTDAFd,"# data written by qreg command.\n");
        fprintf(PMTDAFd,"nvar(\n");
        fprintf(PMTDAFd,"  dfile = %s,\n",PMF1dName);
        fprintf(PMTDAFd,"  noc = %d,\n",nrec);
        k = 0;
        fprintf(PMTDAFd,"  CaseID");
        fprnchar(PMTDAFd,' ',VNameLen - 6,0);
        fprintf(PMTDAFd," [6.0] = c%-2d,\n",++k);

        fprintf(PMTDAFd,"  Wave");
        fprnchar(PMTDAFd,' ',VNameLen - 4,0);
        fprintf(PMTDAFd," [2.0] = c%-2d,\n",++k);

        for (j = -1; j < NPX; ++j) {
            if (j < 0)
                i = PYVar[0];
            else if (j < NPX)
                i = PXVar[j][0];

            fprintf(PMTDAFd,"  %s",VName[i]);
            fprnchar(PMTDAFd,' ',VNameLen - strlen(VName[i]),0);
            fprintf(PMTDAFd," [%d.%d] = c%-2d,\n",PMFmt1,PMFmt2,++k);
        }
        for (j = 0; j < NPZ; ++j) {
            for (l = 0; l < NYCat; ++l) {
                i = PZVar[j][l];
                fprintf(PMTDAFd,"  %s",VName[i]);
                fprnchar(PMTDAFd,' ',VNameLen - strlen(VName[i]),0);
                fprintf(PMTDAFd," [%d.%d] = c%-2d,\n",PMFmt1,PMFmt2,++k);
            }
        }
        if (WIVar >= 0) {
            fprintf(PMTDAFd,"  %s",VName[WIVar]);
            fprnchar(PMTDAFd,' ',VNameLen - strlen(VName[WIVar]),0);
            fprintf(PMTDAFd," [%d.%d] = c%-2d,\n",PMFmt1,PMFmt2,++k);
        }
        if (pflag) {
            for (j = 0; j <= NYCat; ++j) {
                if (j == 0)
                    fprintf(PMTDAFd,"  PROB  ");
                else
                    fprintf(PMTDAFd,"  PROB%-2d",j - 1);
                fprnchar(PMTDAFd,' ',VNameLen - 6,0);
                fprintf(PMTDAFd," [%d.%d] = c%-2d,\n",PMFmt1,PMFmt2,++k);
            }
        }
        fprintf(PMTDAFd,");\n");
        printf1("TDA description written to: %s\n",PMTDAFName);
    }
}

/* ------------------------------------------------------------------------ */
/*  qr_prob1(mod,icase,wave,par,prob)                                       */
/*                                                                          */
/*  mod is QRLOG1 or QRPROB1. Return probability for icase and wave:        */
/*  prob[0] = prob for actual choice                                        */
/*  prob[1] = prob for Y = 0                                                */
/*  prob[2] = prob for Y = 1                                                */

void qr_prob1(int mod,int icase,int wave,double *par,double *prob)
{
    register int j,k;
    double xb,tmp;

    tmp = 0.0;
    if (IntFlg)
        xb = par[1];
    else
        xb = 0.0;
    for (j = IntFlg; j < NParm; ++j) {
        k = j - IntFlg;
        xb += par[j + 1] * get_data(PXVar[k][wave],icase);
    }
    if (mod == QRLOG1) {            /* binary logit */
        tmp = rexp(xb);
        tmp /= (tmp + 1.0);
    }
    else if (mod == QRPROB1)        /* binary probit */
        tmp = cdnf(xb);

    prob[2] = tmp;
    prob[1] = 1.0 - tmp;

    k = (int) PYVar[wave];         /* dependent variable */

    j = YCatI[(int)get_data(k,icase)];   /* internal category */
    if (j == 0)  
        prob[0] = prob[1]; 
    else
        prob[0] = prob[2];
}

/* ------------------------------------------------------------------------ */
/*  qr_prob2(mod,icase,wave,par,prob)                                       */
/*                                                                          */
/*  mod is QRLOG2 or QRPROB2. Return probabilitiesfor icase and wave        */
/*  in prob[].                                                              */

void qr_prob2(int mod,int icase,int wave,double *par,double *prob)
{
    register int i,j,k;
    double xb,tmp;

    tmp = 0.0;

    k = (int) PYVar[wave];                /* dependent variable */
    k = YCatI[(int)get_data(k,icase)];    /* internal category */
    
    j = 0;
    xb = 0.0;
    for (i = NYCat; i <= NParm; ++i)   
        xb += par[i] * get_data(PXVar[j++][wave],icase);
    
    for (k = 1; k < NYCat; ++k) {
        if (mod == QRLOG2) {
            tmp = rexp(par[k] + xb);
            prob[k] = tmp / (1.0 + tmp);
        }
        else
            prob[k] = cdnf(par[k] + xb);
    }
    prob[NYCat] = 0.0;

    k = (int) PYVar[wave];         /* dependent variable */

    j = YCatI[(int)get_data(k,icase)];   /* internal category */
    prob[0] = prob[j + 1];
}

/* ------------------------------------------------------------------------ */
/*  qr_log3(icase,wave,par,prob)                                            */
/*                                                                          */
/*  Return probabilities for icase and wave based on                        */
/*  multinomial logit model. Prob for actual choice in prob[0].             */
/*                                                                          */

void qr_log3(int icase,int wave,double *par,double *prob)
{
    register int i,j,l,u,k1,k2;
    double tmp,d;  

    d = 1.0;        /* d is the inclusive value */
    u = 1;
    for (j = 2; j <= NYCat; ++j) {   

        if (IntFlg)
            tmp = par[u++];
        else
            tmp = 0.0;

        for (i = 0; i < NPX; ++i)  
            tmp += par[u++] * get_data(PXVar[i][wave],icase);
           
        l = (NYCat - 1) * NPX1 + 1;
    
        for (i = 0; i < NPZ; ++i) {
            k1 = PZVar[i][wave * NYCat + j - 1];
            k2 = PZVar[i][wave * NYCat];
            tmp += par[l++] * (get_data(k1,icase) - get_data(k2,icase));
        }
        prob[j] = rexp(tmp);
        d += prob[j];
    }
    prob[1] = 1.0;
    for (i = 2; i <= NYCat; ++i) {   
        prob[i] /= d;
        prob[1] -= prob[i];
    }
    u = (int) PYVar[wave];         /* dependent variable */

    j = YCatI[(int)get_data(u,icase)];   /* internal category */
    prob[0] = prob[j + 1];
}

/* ------------------------------------------------------------------------ */
/*  qr_corr()   print correlation matrix for multinomial probit model.      */

void qr_corr(void)
{
    register int k,l,k1,l1,kk;
    int m,jptr,nn;
    double s,ss[DMVMax + 1];
         
    if (QRTyp == QRPROB3) {
        nn = NYCat;
        jptr = (nn - 1) * NPX1 + NPZ;    /* points to sigma in par */
    }
    else {
        nn = NWave;
        jptr = nn * NPX1; 
    }
    Par[0] = 1.0;
    for (k = 1; k <= nn; ++k) {
        for (l = 1; l <= nn; ++l) {
                    
            if (PMOPT == 3) {
                if (k == 1)  
                    k1 = 0;
                else  
                    k1 = jptr + ((k - 1) * (k - 2) / 2) + 1;
                if (l == 1)  
                    l1 = 0;
                else  
                    l1 = jptr + ((l - 1) * (l - 2) / 2) + 1;

                m = imin(k,l);
                s = 0.0;
                for (kk = 1; kk < m; ++kk)  
                    s += Par[l1++] * Par[k1++];
                if (k > m)
                    s += Par[k1];
                else if (l > m)
                    s += Par[l1];
                else
                    s += 1.0;
            }
            else {
                if (k == l)
                    s = 1.0;
                else {
                    if (l < k)      
                        s = Par[jptr + ((k - 1) * (k - 2)) / 2 + l];
                    else               
                        s = Par[jptr + ((l - 1) * (l - 2)) / 2 + k];

                    if (PMOPT == 2) {
                        s = rexp(s);
                        s = 2.0 * s / (1.0 + s) - 1.0;
                    }
                }
            }
            QRProb[(k - 1) * nn + l] = s;           
        }
    }
    for (k = 1; k <= nn; ++k) {
        s = QRProb[(k - 1) * nn + k];
        if (s <= 0.0) {
            printf1("Warning: estimated correlation matrix not positive definite.\n");
            return;
        }
        ss[k] = sqrt(s);
    }
    for (k = 1; k <= nn; ++k) {
        fprintf(PMPPFd,"# ");
        for (l = 1; l <= nn; ++l) {
            s = QRProb[(k - 1) * nn + l] / (ss[k] * ss[l]);     
            fprintf(PMPPFd,PMTFmtS,s);
        }
        fprintf(PMPPFd,"\n");
    }
    printf1("Estimated correlation matrix written to: %s\n",PMPPFName);
}

/* ------------------------------------------------------------------------ */
/*  rmod()          Rasch model.                                            */
/*  ##                                                                      */
/*                  rmod(                                                   */
/*                      m=...,      model number, def. 1                    */
/*                                   1 = simple Rasch model                 */
/*                      maxcat=..., max number of categories, def. 1000     */
/*                      dsv=...,    starting values                         */
/*                      ppar=...,   print parameter                         */
/*                      pcov=...,   cov matrix                              */ 
/*                      tfmt=...,   print format parameter, def. 10.4       */
/*                      mfmt=...,   print format pcov, def. 10.6            */
/*                      df=...,     write combinations to output file       */
/*                  ) = varlist;                                            */
/*                                                                          */
/*                  Return 0 if OK, otherwise -1.                           */

int rmod(void)
{
    register int i,j,k;
    int err,nv,n;
    double d1,d2;
           
    err = -1;
    if (check_cmd(2))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    printf1("Rasch models. Current memory: %d bytes.\n",MemReq);

    set_mldef();                /* set defaults for ML estimation */

    if (parm(CmdBuf + 4,4,1))   /* get parameters */
        goto RMFin;

    if (PMFmtF == 0)
        pmfmt(10,4);

    if (PMMaxCatFlg == 0)
        PMMaxCat = MaxYC;

    printf1("Number of variables: %d\n",PMNV);

    /* make frequency distribution from data */
    /* max number of categories is: PMMaxCat */

    RMPV = PMNV;

    if (alloc_acr(NOC * RMPV + 1))          /* raw data */
        goto RMFin;
    if (alloc_acs(PMMaxCat * RMPV + 1))     /* data combinations */
        goto RMFin;
    if (alloc_aci(PMMaxCat + 1))            /* frequencies */
        goto RMFin;
    if (alloc_acj(PMMaxCat + 1))            /* bf pointer */
        goto RMFin;

    for (i = 0; i < NOC; ++i) {
        for (j = 0; j < RMPV; ++j) {
            k = (int)get_data(PMVIdx[j],i);
            if (k < 0 || k > 1) {
                printf1("Error: all values must be 0 or 1.\n");
                goto RMFin;
            }
            AcR[i * RMPV + j + 1] = k;
        }
    }
    RMPN = cfreq(RMPV,NOC,AcR,PMMaxCat,AcS,AcI,AcJ,0,&d1,&d2);

    if (RMPN <= 0) {
        printf1("Can't create distribution of categories of variables.\n");
        printf1("Possibly too many categories, current maximum is %d.\n",PMMaxCat);
        goto RMFin;
    }
    printf1("Number of patterns: %d\n\n",RMPN);

    /* save in RMPAT, RMPF */

    if (!(RMPAT = (short *)calloc(RMPN * RMPV + 1,sizeof(short))))  
        goto RMFin;  
    RMPATA = RMPN * RMPV + 1;
    memrq(RMPATA,sizeof(short));

    if (!(RMPF = (int *)calloc(RMPN + 1,sizeof(int))))  
        goto RMFin;  
    RMPFA = RMPN + 1;
    memrq(RMPFA,sizeof(int));

    for (i = 1; i <= RMPN; ++i) {
        k = AcJ[i];
        for (j = 1; j <= RMPV; ++j) {
            RMPAT[(i - 1) * RMPV + j] = AcS[(k - 1) * RMPV + j];

            if (PMF1Def)  
                fprintf(PMF1d,"%d ",RMPAT[(i-1) * RMPV + j]);
        }
        RMPF[i] = AcI[k];

        if (PMF1Def)  
            fprintf(PMF1d,"%6d\n",RMPF[i]);
    }
    NParm = RMPV + RMPN;

    alloc_acr(0);
    alloc_acs(0);
    alloc_aci(0);
    alloc_acj(0);




    PMDOPT = 0;
    if (MINA == 5 || MINA == 6)
        MINA = 4;
    else if (MINA == 7 || MINA > 8)
        MINA = 8;

    set_mlopt();            /* adjust ML options */         
            
    if (ml_init(1,0,0))     /* init ML estimation */
        goto RMFin;

    if (get_dsv(NParm,Par,0,ParLB,ParUB,1))   /* try to get starting values */
        goto RMFin;
  
/*  if (DSVFlg == 0)
        qr_sval();     */   /* default starting values */

    prot_init(2,RMOD1,0);   /* init protocol file */

    /**
    if (RMPN < NParm1) {
        printf1("Error: insufficient number of cases for %d parameters.\n",NParm1);
        goto RMFin;
    }
    **/

    if (ffmin(RMOD1,0,1,1))   /* minimization */
        goto RMFin;
goto RMFin;
    prn_mlres(0);           /* print results of min algorithm */

    if (QRTyp == QRLOG4)  
        printf1("Number of cases used for likelihood: %d\n\n",NFLUsed);

    prn_qcoeff(QRTyp);      /* print estimated parameters etc. */   
           
    if (LConv >= 0)  
        prn_ml1res(0);      /* additional results */
      
    if (PMF1Def)            /* print data to output file */
        qr_pdata();       

    qr_options(QRTyp);      /* additional options */

    if ((QRTyp == QRPROB3 || QRTyp == QRPROB4) && PMPPFDef)  /* write correlation matrix */
        qr_corr();       

    err = 0;

RMFin:
    if (RMPATA > 0) {
        free((char *)RMPAT);
        memrq(-RMPATA,sizeof(short));
        RMPATA = 0;        
    }
    if (RMPFA > 0) {
        free((char *)RMPF);
        memrq(-RMPFA,sizeof(int));
        RMPFA = 0;        
    }
    ml_init(0,0,0);      
    p_clean();
    return(err);
}


