/****************************************************************************/
/*  t_seqm                                                                  */
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

#include "tda.h"
#include "t_gen.h"
#include "t_pgen.h"
#include "t_parm.h"
#include "t_seq.h"
#include "t_gdat.h"
#include "t_var.h"
#include "t_rand.h"
#include "t_gf.h"
#include "t_mat.h"
#include "t_alloc.h"
#include "t_lp.h"

/*  functions in t_seqm.c */

int seqm(void);
int get_costf(int opt);
int seqm_init(int opt);
void seq_pvar(int i,int j);
int seq_gst(int i,int opt);
void seq_pps(int i,int j,int opt,int ic1,int ic2);
void seq_pmd(void);
void seq_plcs(void);
void seq_pseq(void);
void seq_dmat(void);
int seq_pmatch(void);
float seq_mdist(void);
void seqm_pcost(void);
void seqm_dtda(int n);
int seq_scost(void);
int seqpm(void);
int seq_search(void);
int seq_sfind(int k);
void seq_ppt(int k);
int seq_pprn(int i);
int seqm_alloc1(int opt);
void seqpm_dtda(int n);
int subm(void);
double subm_cost(int i,int j);   

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */

int SeqMTyp = 0;    /* type of distance measure                             */
int SeqSN = 0;      /* selected sequence                                    */
int SeqSN1 = 0;     /* first sequence                                       */
int SeqSN2 = 0;     /* second sequence                                      */
int SeqMLen = 0;    /* max length of sequences                              */
int SeqNS = 0;      /* number of states                                     */
int SeqTA = 0;      /* start of time axis                                   */  
int SeqTB = 0;      /* end of time axis                                     */  

int SeqSIL = 0;     /* length of first sequence                             */
int SeqSJL = 0;     /* length of second sequence                            */
int *SeqSI;         /* first sequence                                       */
int *SeqSJ;         /* second sequence                                      */
int *SeqSS;         /* work array                                           */
int SeqSIA = 0;
int SeqSJA = 0;
int SeqSSA = 0;

int ICostTyp = 0;   /* type of indel cost                                   */
int SCostTyp = 0;   /* type of substitution cost                            */

float ICAlpha = 0.0;
float ICBeta = 0.0;

float *ICost;       /* optimal matching: insertion cost                     */
float *SCost;       /* optimal matching: substitution cost                  */
int ICostA = 0;
int SCostA = 0;

int SeqMDim = 0;    /* SeqMlen + 1                                          */
float *SeqMD;       /* minimum cost (D matrix)                              */
float *SeqME;       
float *SeqMF;       
int SeqMDA = 0;
int SeqMEA = 0;
int SeqMFA = 0;
int SeqNMD = 0;     /* number of alignments                                 */

int SeqSPL = 0;     /* max length of test pattern strings                   */
char *SeqDP;        /* DP matrix for pattern search                         */
int SeqDPA = 0;

/*--------------------------------------------------------------------------*/
/*  seqm()      Optimal matching                                            */
/*              Return 0 if OK, -1 if error.                                */

int seqm(void)
{
    register int i;
    int err,n,m;

    err = -1;         
    if (check_cmd(0))
        return(-1);
         
    printf1("Sequence proximity measures. Current memory: %d bytes.\n",MemReq);
         
    if (parm(CmdBuf + 4,1,1))    /* get parameters */
        goto SEQMFin;

    if (PMTFmtF == 0)
        pmtfmt(5,2);

    if (PMM < 1)
        PMM = 1;

    switch (PMM) {

        case 1: printf1("Optimal matching.\n");
                SeqSN = seq_getsn(PMSN,1);   /* get sequence number */
                if (SeqSN < 0)  
                    goto SEQMFin;

                SeqSN1 = SeqSN2 = SeqSN;
                SeqTA = SeqTMin[SeqSN];     /* range of time axis */
                SeqTB = SeqTMax[SeqSN];
                break;

        case 2: printf1("Comparing parallel sequences.\n");
                if (PMSN < 1 || PMSN1 < 1 ||
                         (SeqSN1 = seq_getsn(PMSN,0)) < 0 ||
                                  (SeqSN2 = seq_getsn(PMSN1,0)) < 0) { 

                    printf1("Error: need two valid sequence numbers.\n");
                    goto SEQMFin;
                }
                printf1("Using sequence data structures: %d and %d.\n",PMSN,PMSN1);
                if (SeqSN1 == SeqSN2) {
                    printf1("Error: should be different sequences.\n");
                    goto SEQMFin;
                }
                SeqSN = SeqSN1;                     

                if (seq_scheck(SeqSN1,SeqSN2)) {
                    printf1("Error: must have identical state space.\n");
                    goto SEQMFin;
                }
                SeqTA = imin(SeqTMin[SeqSN1],SeqTMin[SeqSN2]);   /* range of time axis */
                SeqTB = imax(SeqTMax[SeqSN1],SeqTMax[SeqSN2]);    

                break;


        default:  printf1("Error: method %d not available.\n",PMM);
                  goto SEQMFin;
    }
    SeqNS = SeqSTN[SeqSN];          /* number of states */
    SeqMLen = SeqTB - SeqTA + 1;    /* max sequence length */
    SeqMDim = SeqMLen + 1;          /* dimension of D matrix */

    printf1("Number of states: %d. Max sequence length: %d\n",SeqNS,SeqMLen);
    if (PMSM[1])  
        printf1("Option (sm=1): skip internal gaps.\n");
    if (PMSM[2])  
        printf1("Option (sm=2): skip identical states.\n");
    if (PMRRN == 1)
        printf1("Using common sequence length (rr=1).\n");

    if (PMF1Def)
        printf1("Test output will be written to: %s\n",PMF1dName);
    newline();

    if (get_costf(1))       /* set up cost functions */
        goto SEQMFin;

    if (PMF1Def)            /* print to test output file */
        seqm_pcost();

    if (ICostTyp == 1)      /* for linear indel cost functions */
        SeqMTyp = 1;
    else
        SeqMTyp = 0;

    if (PMMax >= 1 && SeqMTyp == 0)  
        printf1("Restricted alignment: max=%d\n",PMMax);
    else
        PMMax = 0;

    if (seqm_init(1))       /* init D, E, F matrices */
        goto SEQMFin;

    /* check cn option */

    if (PMM == 1 && PMNCN > 0) {

        printf1("Checking cn option.\n");
        m = 0;
        for (i = 0; i < PMNCN; ++i) {
            if (PMCN[i] < 1 || PMCN[i] > NOC) {
                printf1("Error: cn parameter inconsistent with number of cases (= %d).\n",NOC);
                goto SEQMFin;
            }
            n = PMCN[i];
            if (seq_gst(n - 1,1) < 1) {
                printf1("Warning: sequence %d has zero length and will be ignored.\n",n);
                PMCN[i] = -1;
            }
            else 
                m++;
        }
        printf1("Will use %d reference sequences.\n",m);
        if (m == 0)
            goto SEQMFin;
    }
    n = seq_pmatch();               /* alignment */

    if (PMTDAFDef && n > 0)         /* write TDA description */
        seqm_dtda(n);

    err = 0;
   
SEQMFin:
    seqm_init(0);
    get_costf(0);
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  get_costf(opt)                                                          */
/*                                                                          */
/*  If opt != 0 set up cost functions. We have the following options:       */
/*                                                                          */
/*  ICostTyp = 0    fixed in ICost[]                                        */
/*             1    icost = alpha,beta (linear indel cost function)         */
/*             2    icost = data based                                      */
/*                                                                          */
/*  SCostTyp = 0    fixed in SCost[]                                        */
/*             1    icost = absolute difference                             */
/*             2    icost = data based                                      */
/*                                                                          */
/*  Also adjust PMS.                                                        */
/*  Also allocate SeqSI[], SeqSJ[], and SeqSS[].                            */
/*                                                                          */
/*  If opt = 0 free previously allocated memory.                            */

int get_costf(int opt)
{
    register int i,j;
    int err,ic,sc,n,im,imr,imc,sm,smr,smc;

    if (opt == 0) {
        err = 0;
        goto COSTFFin;
    }
    err = -2;

    n = SeqMLen;
    if (!(SeqSI = (int *)calloc(n,sizeof(int))))   
        goto COSTFFin;
    memrq(n,sizeof(int));
    SeqSIA = n;

    if (!(SeqSJ = (int *)calloc(n,sizeof(int))))   
        goto COSTFFin;
    memrq(n,sizeof(int));
    SeqSJA = n;

    if (!(SeqSS = (int *)calloc(n,sizeof(int))))   
        goto COSTFFin;
    memrq(n,sizeof(int));
    SeqSSA = n;

    if (!(ICost = (float *)calloc(n,sizeof(float))))   
        goto COSTFFin;
    memrq(n,sizeof(float));
    ICostA = n;

    n = SeqNS * SeqNS;  
    if (!(SCost = (float *)calloc(n,sizeof(float))))   
        goto COSTFFin;
    memrq(n,sizeof(float));
    SCostA = n;
            
    err = -1;
    ic = sc = 0;

    /* indel cost */

    if (PMICOSTMAT >= 0) {
        ICostTyp = 99;
        im = PMICOSTMAT;
        printf1("Indel cost defined by matrix: %s\n",MatName[im]);
        imr = MatRow[im];
        imc = MatCol[im];
        if ((imr != 1 || imc != SeqMLen) && (imr != SeqMLen || imc != 1)) {
            printf1("Error: need an (1,%d) or (%d,1) matrix.\n",SeqMLen,SeqMLen);
            goto COSTFFin;
        }
    }
    else if (PMICOSTA >= 0.0 && PMICOSTB >= 0.0) {
        ICostTyp = 1;
        ICAlpha = PMICOSTA;
        ICBeta = PMICOSTB;
        printf1("Linear indel cost function: alpha=%g, beta=%g.\n",ICAlpha,ICBeta);
    }
    else if (PMICOSTA >= 0.0) {
        ICostTyp = 0;
        printf1("Indel cost: %g.\n",PMICOSTA);
    }
    else {
        ic = 1;
        ICostTyp = 0;
        PMICOSTA = 1.0;
        printf1("Default indel cost: 1.\n");
    }

    /* substitution cost */

    if (PMSCOSTMAT >= 0) {
        SCostTyp = 99;
        sm = PMSCOSTMAT;
        printf1("Substitution cost defined by matrix: %s\n",MatName[sm]);
        smr = MatRow[sm];
        smc = MatCol[sm];
        if (smr < SeqNS || smc < SeqNS) {
            printf1("Error: need at least an (%d,%d) matrix.\n",SeqNS,SeqNS);
            goto COSTFFin;
        }
    }
    else if (PMSCOSTM == 1) {
        SCostTyp = 1;
        printf1("Substitution cost defined by absolute difference.\n");
    }
    else if (PMSCOSTM == 2) {
        SCostTyp = PMSCOSTM;
        printf1("Substitution cost based on data, type %d.\n",SCostTyp);
    }
    else {
        sc = 1;
        SCostTyp = 0;
        printf1("Default substitution cost: 2.\n");
    }
    if (ic == 0 || sc == 0) {     /* icost or scost parameters used */
        if (PMS > 1)
            PMS = 0;
    }

    /* set up indel cost in ICost[] */

    for (i = 0; i < SeqMLen; ++i)
        if (ICostTyp == 0)
            ICost[i] = PMICOSTA;
        else if (ICostTyp == 99) {
            ICost[i] = (float)MatVal[im][i + 1];
    }
    if (ICostTyp == 99)
        ICostTyp = 0;

    /* set up substitution cost in SCost[] */

    if (SCostTyp == 2) {
        seq_scost();
        SCostTyp = 0;
    }
    else {
        for (i = 0; i < SeqNS; ++i) {
            for (j = 0; j < SeqNS; ++j) {
                if (j != i) {
                    if (SCostTyp == 0)  
                        SCost[i * SeqNS + j] = 2.0; 
                    else if (SCostTyp == 99)
                        SCost[i * SeqNS + j] = MatVal[sm][i * smc + j + 1];
                }
            }
        }
        if (SCostTyp == 99)
            SCostTyp = 0;
    }
    return(0);

COSTFFin:
    if (err == -2)  
        p_err(-2,1);

    if (SeqSIA > 0) {
        free((char *)SeqSI);
        memrq(-SeqSIA,sizeof(int));
        SeqSIA = 0;
    }
    if (SeqSJA > 0) {
        free((char *)SeqSJ);
        memrq(-SeqSJA,sizeof(int));
        SeqSJA = 0;
    }
    if (SeqSSA > 0) {
        free((char *)SeqSS);
        memrq(-SeqSSA,sizeof(int));
        SeqSSA = 0;
    }
    if (ICostA > 0) {
        free((char *)ICost);
        memrq(-ICostA,sizeof(float));
        ICostA = 0;         
    }
    if (SCostA > 0) {
        free((char *)SCost);
        memrq(-SCostA,sizeof(float));
        SCostA = 0;         
    }
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqm_init(opt)  If opt != 0 allocate and initialize data structure.     */
/*                  Otherwise free previously allocated memory.             */
/*                                                                          */
/*                  Always allocate and initialize: SeqMD[]                 */
/*                  If SeqMTyp = 1, also allocate and initialize            */
/*                  SeqME[] and SeqMF[].                                    */
/*                                                                          */
/*                  If PMF1Def and PMTST[1] print                           */
/*                  matrices to test output file.                           */
/*                                                                          */
/*                  Return 0 if OK, -1 if error.                            */

int seqm_init(int opt)
{
    register int i;
    int err,n;
    float tmp;

    if (opt == 0) {
        err = 0;
        goto SEQMIFin;
    }
    err = -1;

    n = SeqMDim * SeqMDim;
    if (!(SeqMD = (float *)calloc(n,sizeof(float))))   
        goto SEQMIFin;
    memrq(n,sizeof(float));
    SeqMDA = n;

    if (SeqMTyp == 1) {
        if (!(SeqME = (float *)calloc(n,sizeof(float))))   
            goto SEQMIFin;
        memrq(n,sizeof(float));
        SeqMEA = n;

        if (!(SeqMF = (float *)calloc(n,sizeof(float))))   
            goto SEQMIFin;
        memrq(n,sizeof(float));
        SeqMFA = n;
    }
    SeqMD[0] = 0.0;                 /* init SeqMD */

    if (SeqMTyp == 0) {
        for (i = 0; i < SeqMLen; ++i) {
            SeqMD[i + 1] = SeqMD[i] + ICost[i];
            SeqMD[(i + 1) * SeqMDim] = SeqMD[i * SeqMDim] + ICost[i];
        }
    }
    else if (SeqMTyp == 1) {        /* also SeqME and SeqMF */

        SeqME[0] = SeqMF[0] = 0.0;

        for (i = 0; i < SeqMLen; ++i) {
            tmp = (float)ICAlpha + (float)ICBeta * (float)i;
            SeqMD[i + 1] = SeqMD[(i + 1) * SeqMDim] = tmp;
            SeqMF[i + 1] =                                   
            SeqME[(i + 1) * SeqMDim] = tmp + (float)ICAlpha;
        }
    }
    if (PMF1Def && PMTST[1]) {
        fprintf(PMF1d,"Initial D matrix\n");
        mprf(SeqMDim,SeqMDim,SeqMD - 1,PMFmtS,PMF1d);

        if (SeqMTyp == 1) {        /* also SeqME and SeqMF */

            fprintf(PMF1d,"\nInitial E matrix\n");
            mprf(SeqMDim,SeqMDim,SeqME - 1,PMFmtS,PMF1d);

            fprintf(PMF1d,"\nInitial F matrix\n");
            mprf(SeqMDim,SeqMDim,SeqMF - 1,PMFmtS,PMF1d);
        }
        fprintf(PMF1d,"\n");
    }
    return(0);

SEQMIFin:
    if (err)  
        p_err(-2,1);

    if (SeqMDA > 0) {
        free((char *)SeqMD);
        memrq(-SeqMDA,sizeof(float));
        SeqMDA = 0;
    }
    if (SeqMEA > 0) {
        free((char *)SeqME);
        memrq(-SeqMEA,sizeof(float));
        SeqMEA = 0;
    }
    if (SeqMFA > 0) {
        free((char *)SeqMF);
        memrq(-SeqMFA,sizeof(float));
        SeqMFA = 0;
    }
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seq_pvar(i,j) Print additional variables to output file, for cases i,j  */

void seq_pvar(int i,int j)
{
    register int k,iv;

    for (k = 0; k < PMNV; ++k) {
        iv = (int)PMVIdx[k];
        fprintf(PMFd,VPFmtS[iv],get_data(iv,i));
        if (j >= 0)
            fprintf(PMFd,VPFmtS[iv],get_data(iv,j));
    }
}

/*--------------------------------------------------------------------------*/
/*  seq_gst(i,opt)  Get sequence for case i, return length.                 */
/*                  If opt = 1, use SeqSI, else SeqSJ.                      */
/*                  Return 0 for zero length, -1 for internal gaps.         */
/*                  Note: using internal state numbers.                     */
/*                                                                          */
/*  Now this function recognizes PMSM[].                                    */
/*  PMSM[1] : skip internal gaps                                            */
/*  PMSM[2] : skip identical states                                         */

int seq_gst(int i,int opt)
{
    register int t,l,a,b,n,nn;

    a = -1;
    for (t = SeqTA; t <= SeqTB; ++t) {
        n = seq_sget(i,t,SeqSN);
        if (n >= 0) {
            a = t;
            break;
        }
    }
    if (a < 0)
        return(0);

    b = -1;
    for (t = SeqTB; t >= a; --t) {
        n = seq_sget(i,t,SeqSN);
        if (n >= 0) {
            b = t;
            break;
        }
    }
    nn = l = 0;  
    for (t = a; t <= b; ++t) {

        n = seq_sget(i,t,SeqSN);
        if (n < 0) {
            if (PMSM[1])
                continue;
            return(-1);
        }
        if (PMSM[2] && l && n == nn)
            continue;
        nn = n;
        n = SeqSTNII[SeqSN][n];
        if (opt == 1)
            SeqSI[l++] = n;
        else
            SeqSJ[l++] = n;
    }
    return(l);
}

/*--------------------------------------------------------------------------*/
/*  seq_pps(i,j)    Print sequences to test output file. If i != 0, print   */
/*                  first sequence (SeqSI), if j != 0, print second         */
/*                  sequence (SeqSJ).                                       */

void seq_pps(int i,int j,int opt,int ic1,int ic2)
{
    register int k;

    if (PMF1Def && PMTST[2]) {  
        if (i) {
            fprintf(PMF1d,"Sequence A (case number %d)\n",ic1 + 1);
            for (k = 0; k < SeqSIL; ++k)  
                fprintf(PMF1d," %2d",SeqSTNI[SeqSN1][SeqSI[k]]);
            fprintf(PMF1d,"\n");
        }
        if (j) {
            fprintf(PMF1d,"Sequence B (case number %d)\n",ic2 + 1);
            for (k = 0; k < SeqSJL; ++k)  
                fprintf(PMF1d," %2d",SeqSTNI[SeqSN2][SeqSJ[k]]);
            fprintf(PMF1d,"\n");
        }   
        if (opt)
            fprintf(PMF1d,"\n");    
    }
}
   
/*--------------------------------------------------------------------------*/
/*  seq_pmd()   Print number of alignments to stderr.                       */

void seq_pmd(void)
{
    if (SILENTFlg >= 2)
        return;

    SeqNMD++;
    if ((SeqNMD / 100) * 100 == SeqNMD)
        printfe("Done: %7d     %c",SeqNMD,CR);
}

/*--------------------------------------------------------------------------*/
/*  seq_plcs()      Print longest common subsequence.                       */

void seq_plcs(void)
{
    register int i,j,k,ii;
    int dij;

    k = 0;
    i = SeqSIL - 1;
    j = SeqSJL - 1;
    while (i >= 0 && j >= 0) {
        ii = (i + 1) * SeqMDim;
        dij = (int)SeqMD[ii + j + 1];
        if (dij == (int)SeqMD[i * SeqMDim + j + 1] + (int)ICost[i])
            i--;
        else if (dij == (int)SeqMD[ii + j] + (int)ICost[j])
            j--;
        else {
            SeqSS[k++] = SeqSI[i];
            i--;
            j--;
        }
    }
    fprintf(PMFd,"%4d ",k);    
  
    for (i = k - 1; i >= 0; --i)
        fprintf(PMFd,"%2d ",SeqSTNI[SeqSN][SeqSS[i]]);
    for (i = k; i < SeqMLen; ++i)
        fprintf(PMFd,"-1 ");
}

/*--------------------------------------------------------------------------*/
/*  seq_pseq()      Print sequential distances from SeqMD.                  */

void seq_pseq(void)
{
    register int i,k;

    k = imin(SeqSIL,SeqSJL);
    for (i = 1; i <= k; ++i)   
        fprintf(PMFd,PMTFmtS,(double)SeqMD[i * SeqMDim + i]);
       
    if (k < SeqSIL) {
        for (i = k + 1; i <= SeqSIL; ++i)  
            fprintf(PMFd,PMTFmtS,(double)SeqMD[i * SeqMDim + SeqSJL]);
    }
    else if (k < SeqSJL) {
        for (i = k + 1; i <= SeqSJL; ++i)  
            fprintf(PMFd,PMTFmtS,(double)SeqMD[SeqSIL * SeqMDim + i]);
    }
    k = imax(SeqSIL,SeqSJL);
    while (k++ < SeqMLen)
        fprintf(PMFd,PMTFmtS,-1.0);
}

/*--------------------------------------------------------------------------*/
/*  seq_dmat()          Print D matrix to test output file.                 */
/*                      If SeqMTyp = 1 : also E and F matrices              */

void seq_dmat(void)
{
    register int i,j,im;
    int is,js,n;
    float tmp;

    for (im = 0; im < 3; ++im) {
        if (im == 0)  
            fprintf(PMF1d,"D");
        else if (im == 1)
            fprintf(PMF1d,"E");
        else if (im == 2)
            fprintf(PMF1d,"F");

        fprintf(PMF1d," Matrix  ");

        /* fprnchar(PMF1d,' ',10,0); */ 
        for (j = 0; j <= SeqSJL; ++j)  
            fprintf(PMF1d,"%7d ",j);
        fprintf(PMF1d,"\n         ");
        fprnchar(PMF1d,'-',8 + 8 * SeqSJL,1);
        fprintf(PMF1d,"                B ");
        for (j = 0; j < SeqSJL; ++j) {
            js = SeqSJ[j];
            n = SeqSTNI[SeqSN][js];
            fprintf(PMF1d,"%7d ",n);
        }
        fprintf(PMF1d,"\n         ");
        fprnchar(PMF1d,'-',8 + 8 * SeqSJL,1);
    
        for (i = 0; i <= SeqSIL; ++i) {
            if (i == 0)
                fprintf(PMF1d,"%4d   A |",i );
            else {
                is = SeqSI[i - 1];
                n = SeqSTNI[SeqSN][is];
                fprintf(PMF1d,"%4d %3d |",i,n);
            }
            for (j = 0; j <= SeqSJL; ++j) {

                if (im == 0)
                    tmp = SeqMD[i * SeqMDim + j];
                else if (im == 1)
                    tmp = SeqME[i * SeqMDim + j];
                else if (im == 2)
                    tmp = SeqMF[i * SeqMDim + j];

                fprintf(PMF1d,"%7.2f ",tmp);
            }
            fprintf(PMF1d,"\n");
        }
        fprintf(PMF1d,"\n");
        if (SeqMTyp == 0)
            break;
    }
}   
    
/*--------------------------------------------------------------------------*/
/*  seq_pmatch()    Optimal matching, for options:                          */
/*                  PMM = 1 : optimal matching                              */
/*                  PMM = 2 : compare parallel sequences                    */
/*                                                                          */
/*  Return number of records written to output file.                        */

int seq_pmatch(void)
{
    register int i,j,k;
    int nm,m,nn,id,jd,sil,sjl;
    float tmp,dmax;
    double r,rdh;

    printf1("\nStarting alignment procedure.\n");

    SeqNMD = 0;         /* number of alignments */
    nn = m = 0;
    dmax = 0.0;
    id = jd = 0;

    if (PMM == 1) {

        if (PMNCN <= 0) {           /* pairwise comparison */

            nm = NOC * (NOC - 1) / 2;
            if (PMR > 0) {
                printf1("Random selection of approximately %d sequence pairs.\n",PMR);
                if (PMR >= nm) {
                    printf1("Will be ignored.\n");
                    PMR = 0;
                }
                else {
                    rdh = (double)PMR / (double)nm;
                    nm = PMR;
                }
            }
            if (SILENTFlg < 2)
                printfe("Total alignments: %d\n",nm);

            for (i = 0; i < NOC; ++i) {
                sil = seq_gst(i,1);
                if (sil < 1) {
                    m++;
                    continue;
                }
                for (j = 0; j < i; ++j) {
                    sjl = seq_gst(j,2);
                    if (sjl < 1)
                        continue;

                    if (PMRRN == 1)
                        SeqSIL = SeqSJL = imin(sil,sjl);
                    else {
                        SeqSIL = sil;
                        SeqSJL = sjl;
                    }

                    if (PMR > 0) {
                        r = random1();
                        if (r > rdh)
                            continue;
                    }
                    seq_pps(1,1,1,i,j);         /* write to test output */

                    tmp = seq_mdist();          /* calculate distance */
               
                    seq_pmd();                  /* print to stderr */

                    if (dmax < tmp) {
                        dmax = tmp;
                        id = i;
                        jd = j;
                    }
                    fprintf(PMFd,"%6d %6d %6d %6d ",i + 1,j + 1,SeqSIL,SeqSJL);

                    if (PMS == 1)           /* sequential */
                        seq_pseq();
                    else {   
                        fprintf(PMFd,PMTFmtS,(double)tmp);
                        if (PMS == 2)      /* LCS */
                            seq_plcs();
                    }           
                    seq_pvar(i,j);          /* print additional variables */
                    fprintf(PMFd,"\n");
                    nn++;
                }
            }
        }
        else {              /* compare with PMCN[] sequences */

            for (i = 0; i < NOC; ++i) {
                sil = seq_gst(i,1);

                if (sil < 1) {
                    m++;
                    continue;
                }
                fprintf(PMFd,"%5d %3d ",i + 1,sil);

                for (j = 0; j < PMNCN; ++j) {
                    k = PMCN[j] - 1;
                    if (k >= 0) {
                        sjl = seq_gst(k,2);
                        if (sjl < 1)
                            continue;

                        if (PMRRN == 1)
                            SeqSIL = SeqSJL = imin(sil,sjl);
                        else {
                            SeqSIL = sil;
                            SeqSJL = sjl;
                        }
                        seq_pps(1,1,1,i,k);         /* write to test output */

                        tmp = seq_mdist();    
               
                        seq_pmd();                  /* print to stderr */
                     
                        if (PMS == 1)               /* sequential */
                            seq_pseq();
                        else {
                            fprintf(PMFd,PMTFmtS,(double)tmp);
                            if (PMS == 2)          /* LCS */
                                seq_plcs();
                        }           
                    }
                }
                seq_pvar(i,-1);
                fprintf(PMFd,"\n");
                nn++;
            }
        }           
    }
    else if (PMM == 2) {      /* compare two parallel sequences */

        for (i = 0; i < NOC; ++i) {
            SeqSN = SeqSN1;
            SeqSIL = seq_gst(i,1);
            SeqSN = SeqSN2;
            SeqSJL = seq_gst(i,2);

            if (SeqSIL < 1 || SeqSJL < 1) {
                m++;
                continue;
            }
            if (PMRRN == 1)
                SeqSIL = SeqSJL = imin(SeqSIL,SeqSJL);

            seq_pps(1,1,1,i,i);         /* write to test output */
            tmp = seq_mdist();          /* calculate distance */
            seq_pmd();                  /* print to stderr */

            if (dmax < tmp) {
                dmax = tmp;
                id = i;
            }
            fprintf(PMFd,"%5d %3d %3d ",i + 1,SeqSIL,SeqSJL);

            if (PMS == 1)           /* sequential */
                seq_pseq();
            else {   
                fprintf(PMFd,PMTFmtS,(double)tmp);
                if (PMS == 2)      /* LCS */
                    seq_plcs();
            }           
            seq_pvar(i,-1);            
            fprintf(PMFd,"\n");
            nn++;
        }
    }
    printf1("Number of sequences (cases): %d\n",NOC);
    printf1("Sequences with zero length or internal gaps: %d\n",m);
    printf1("Sequences used for alignment: %d\n",NOC - m);
    printf1("\nNumber of alignments: %d\n",SeqNMD);
    printf1("%d record(s) written to output file: %s\n",nn,PMFdName);
    if (nn > 1) {
        if (PMM == 1 && PMNCN <= 0)  
            printf1("Maximum distance between sequences %d and %d: %g\n",id+1,jd+1,dmax);
        else if (PMM == 2)  
            printf1("Maximum distance in case %d: %g\n",id,dmax);
    }
    return(nn);
}

/*--------------------------------------------------------------------------*/
/*  seq_mdist()             Calculate distance between sequences in         */
/*                          SeqSI and SeqSJ.                                */
/*                          Return: distance.                               */

float seq_mdist(void)
{
    register int i,j,k,l,isp;
    int is,js,j1,j2,j2p,nn;
    float itmp,dtmp,stmp;
    float ilen,jlen,u,tmp,tmp1,tmp2;

    nn = SeqNS * SeqNS;

    if (PMMax >= 1) {               /* restricted alignment */
        ilen = (float)(SeqSIL + 1);
        jlen = (float)(SeqSJL + 1);
        u = (ilen / jlen + jlen / ilen) / sqrt(ilen * ilen + jlen * jlen);
        u = (float)PMMax * u / 1.4142;
    }
    j1 = 0;
    j2 = j2p = SeqSJL;

    for (i = 0; i < SeqSIL; ++i) {

        if (PMMax >= 1) {
            tmp = (double)(i + 1) / ilen;    
            tmp1 = jlen * (tmp - u);
            tmp2 = jlen * (tmp + u);
            if (tmp1 <= 1.0)
                j1 = 0;
            else
                j1 = (int)tmp1;
            if (tmp2 >= jlen)
                j2 = SeqSJL;
            else
                j2 = (int)tmp2;
        }
        is = SeqSI[i];                         
        isp = is * SeqNS;
        k = i * SeqMDim + j1;
        l = k + SeqMDim;

        for (j = j1; j < j2; ++j) {

            if (SeqMTyp == 1) {     /* indel functions */

                stmp = SeqMD[l] + (float)ICAlpha;
                if (stmp > SeqME[l] + (float)ICBeta)
                    stmp = SeqME[l] + (float)ICBeta;
                l++;
                SeqME[l] = stmp;

                stmp = SeqMD[k + 1] + (float)ICAlpha;
                if (stmp > SeqMF[k + 1] + (float)ICBeta)
                    stmp = SeqMF[k + 1] + (float)ICBeta;
                SeqMF[l] = stmp;
            }
            js = SeqSJ[j];

            if (is == js)
                stmp = 0.0;

            else if (SCostTyp == 0)     /* time-independent substitution cost */
                stmp = SCost[isp + js];
 
            else if (SCostTyp == 1) {   /* absolute difference */

                stmp = fabs((double)(SeqSTNI[SeqSN][is] - SeqSTNI[SeqSN][js]));
    
            }
            stmp += SeqMD[k];
            k++; 

            if (SeqMTyp == 0) {
                if (j > j1 || j == 0) {
                    itmp = ICost[j] + SeqMD[l];
                    if (stmp > itmp)
                        stmp = itmp;
                }
                l++;
                if (j < j2p) {
                    dtmp = ICost[i] + SeqMD[k];
                    if (stmp > dtmp)
                        stmp = dtmp;
                }
            }
            else if (SeqMTyp == 1) {
                if (stmp > SeqME[l])
                    stmp = SeqME[l];
                if (stmp > SeqMF[l])
                    stmp = SeqMF[l];
            }
            SeqMD[l] = stmp;
        }
        j2p = j2;
    }
    if (PMF1Def && PMTST[3])
        seq_dmat();
    return(stmp);
}

/*--------------------------------------------------------------------------*/
/*  seqm_pcost()    Print cost definitions to test output file.             */

void seqm_pcost(void)
{
    if (PMF1Def == 0)
        return;

    fprintf(PMF1d,"Optimal matching test output file.\n");
    fprintf(PMF1d,"Number of states: %d\n",SeqNS);
    fprintf(PMF1d,"Max sequence lenght: %d\n",SeqMLen);

    if (ICostTyp == 1)
        fprintf(PMF1d,"\nLinear indel cost function: alpha=%g, beta=%g.\n",ICAlpha,ICBeta);
    else {
        fprintf(PMF1d,"\nIndel cost\n");
        mprf(1,SeqMLen,ICost - 1,PMFmtS,PMF1d);
    }
    if (SCostTyp == 0) {
        fprintf(PMF1d,"\nSubstitution cost\n");
        mprf(SeqNS,SeqNS,SCost - 1,PMFmtS,PMF1d);
    }
    else if (SCostTyp == 1)
        fprintf(PMF1d,"\nSubstitution cost: absolute difference.\n");
    fprintf(PMF1d,"\n");
}

/*--------------------------------------------------------------------------*/
/*  seqm_dtda(n)     Write TDA description file, n = number of cases.       */

void seqm_dtda(int n)
{
    register int i,j,k,l;

    fprintf(PMTDAFd,"nvar(\n");
    fprintf(PMTDAFd,"  dfile = %s,\n",PMFdName);
    fprintf(PMTDAFd,"  noc = %d,\n",n);

    k = 1;

    if (PMM == 1) {

        if (PMNCN <= 0) {

            fprintf(PMTDAFd,"  SEQ1N <5>[6.0] = c%-3d, # sequence A case number\n",k++);
            fprintf(PMTDAFd,"  SEQ2N <5>[6.0] = c%-3d, # sequence B case number\n",k++);
            fprintf(PMTDAFd,"  SEQ1L <2>[6.0] = c%-3d, # sequence A length\n",k++);
            fprintf(PMTDAFd,"  SEQ2L <2>[6.0] = c%-3d, # sequence B length\n",k++);

            if (PMS == 1) {
                for (i = 1; i <= SeqMLen; ++i) 
                    fprintf(PMTDAFd,"  DIST%d <4>[%d.%d] = c%-3d, # distance, time %d\n",i,PMTFmt1,PMTFmt2,k++,i);
            }
            else {
                fprintf(PMTDAFd,"  DIST  <4>[%d.%d] = c%-3d, # distance\n",PMTFmt1,PMTFmt2,k++);
                
                if (PMS == 2) {
                    fprintf(PMTDAFd,"  LCSL  <2>[4.0] = c%-3d, # length of LCS\n",k++);
                    for (i = 1; i <= SeqMLen; ++i)
                        fprintf(PMTDAFd,"  LCS%d <2>[2.0] = c%-3d, # LCS t=%d\n",i,k++,i);
                }
            }                 
        }            
        else {
            fprintf(PMTDAFd,"  SEQ1N <5>[6.0] = c%-3d, # sequence A case number\n",k++);
            fprintf(PMTDAFd,"  SEQ1L <2>[6.0] = c%-3d, # sequence A length\n",k++);

            for (i = 0; i < PMNCN; ++i) {
                j = PMCN[i];      
                if (j > 0) {
                    if (PMS == 1) {
                        for (l = 1; l <= SeqMLen; ++l) 
                            fprintf(PMTDAFd,"  DIST%d_%d <4>[%d.%d] = c%-3d, # distance, t=%d, cn=%d\n",j,l,PMTFmt1,PMTFmt2,k++,l,j);
                    }
                    else {
                        fprintf(PMTDAFd,"  DIST%d <4>[%d.%d] = c%-3d, # distance to sequence cn=%d\n",j,PMTFmt1,PMTFmt2,k++,j);
                        if (PMS == 2) {
                            fprintf(PMTDAFd,"  LCSL%d <2>[4.0] = c%-3d, # length of LCS, cn=%d\n",j,k++,j);
                            for (l = 1; l <= SeqMLen; ++l)
                                fprintf(PMTDAFd,"  LCS%d_%d <2>[2.0] = c%-3d, # LCS t=%d, cn=%d\n",j,l,k++,l,j);
                        }
                    }                  
                }
            }
        }
        for (i = 0; i < PMNV; ++i) {
            j = PMVIdx[i];
            fprintf(PMTDAFd,"  %s_A <%d>[%d.%d] = c%-3d, # sequence A\n",       
                                    VName[j],VSLen[j],VPFmt1[j],VPFmt2[j],k++);
            if (PMNCN <= 0)
                fprintf(PMTDAFd,"  %s_B <%d>[%d.%d] = c%-3d, # sequence B\n",         
                                    VName[j],VSLen[j],VPFmt1[j],VPFmt2[j],k++);
        }
    }
    else if (PMM == 2) {

        fprintf(PMTDAFd,"  SEQN  <5>[6.0] = c%-3d, # case number\n",k++);
        fprintf(PMTDAFd,"  SEQ1L <2>[6.0] = c%-3d, # sequence A length\n",k++);
        fprintf(PMTDAFd,"  SEQ2L <2>[6.0] = c%-3d, # sequence B length\n",k++);

        if (PMS == 1) {
            for (l = 1; l <= SeqMLen; ++l) 
                fprintf(PMTDAFd,"  DIST%d_%d <4>[%d.%d] = c%-3d, # distance, t=%d\n",j,l,PMTFmt1,PMTFmt2,k++,l);
        }
        else {         
            fprintf(PMTDAFd,"  DIST  <2>[6.0] = c%-3d, # distance\n",k++);
            if (PMS == 2) {
                fprintf(PMTDAFd,"  LCSL  <2>[4.0] = c%-3d, # length of LCS\n",k++);
                for (i = 1; i <= SeqMLen; ++i)
                    fprintf(PMTDAFd,"  LCS%d <2>[2.0] = c%-3d, # LCS t=%d\n",i,k++,i);
            }
        }
        for (i = 0; i < PMNV; ++i) {
            j = PMVIdx[i];
            fprintf(PMTDAFd,"  %s <%d>[%d.%d] = c%-3d,\n",       
                                    VName[j],VSLen[j],VPFmt1[j],VPFmt2[j],k++);
        }
    }
    fprintf(PMTDAFd,");\n");
    printf1("TDA description written to: %s\n",PMTDAFName);
}

/*--------------------------------------------------------------------------*/
/*  seq_scost()     Calculate substitution costs based on transition        */  
/*                  frequencies.                                            */
/*                                                                          */
/*  Return 0 if successful, otherwise -1.                                   */

int seq_scost(void)
{
    register int i,j,is,js;
    int n,m;
    float tmp;

    m = SeqNS * SeqNS;
    for (i = 0; i < m; ++i)
        SCost[i] = 0.0;
           
    m = 0;
    for (i = 0; i < NOC; ++i) {
        n = seq_gst(i,1);
        if (n < 1) {
            m++;
            continue;
        }
        is = SeqSI[0];
        if (SCostTyp == 2) {
            for (j = 1; j < n; ++j) {
                js = SeqSI[j];
                SCost[is * SeqNS + js] += 1.0;
                is = js;
            }
        }
    }
    if (SCostTyp == 2) {
        for (is = 0; is < SeqNS; ++is) {
            tmp = 0.0;
            j = is * SeqNS;
            for (js = 0; js < SeqNS; ++js)  
                tmp += SCost[j++];

            j = is * SeqNS;
            for (js = 0; js < SeqNS; ++js) {
                if (tmp > 0.0)
                    SCost[j++] /= tmp;
            }
        }
        SCost[0] = 0.0;
        for (is = 1; is < SeqNS; ++is) {
            j = is * SeqNS;
            for (js = 0; js < is; ++js) {
                tmp = 2.0 - SCost[is * SeqNS + js] - SCost[js * SeqNS + is];
                SCost[j++] = tmp;
                SCost[js * SeqNS + is] = tmp;
            }
            SCost[j] = 0.0;
        }
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  seqpm.     Sequence pattern matching.                                   */
/*                                                                          */  
/*          seqpm(                                                          */
/*              sn = ...,       number of sequence data structure, def. 1   */
/*              ps = ...,...,   definition of max 20 patterns               */
/*              df = ...,       test output file                            */
/*              nfmt=...,       print format for number of matches, def. 4  */
/*              v=...,          add ... variables to output file            */
/*              dtda=...,       TDA description file                        */
/*          ) = fname;          name of output file                         */
/*                                                                          */
/*          Patterns must be given as follows:                              */
/*                                                                          */
/*          ps = [a1,a2,...],[b1,b2,...],...                                */
/*                                                                          */
/*          The characters may be:                                          */
/*          nonnegative integers for valid states                           */
/*          ?    matches any character                                      */
/*          *    matches any sequence of characters                         */
/*          +                                                               */
/*          -                                                               */
/*                                                                          */
/*          Return 0 if OK, -1 if error.                                    */

int seqpm(void)
{
    int err;

    err = -1;         
    if (check_cmd(0))
        return(-1);
         
    printf1("Sequence pattern matching. Current memory: %d bytes.\n",MemReq);
         
    if (parm(CmdBuf + 5,1,1))    /* get parameters */
        goto SEQPMFin;

    SeqSN = seq_getsn(PMSN,1);   /* get sequence number */
    if (SeqSN < 0)  
        goto SEQPMFin;

    SeqNS = SeqSTN[SeqSN];          /* number of states */
    SeqTA = SeqTMin[SeqSN];         /* range of time axis */
    SeqTB = SeqTMax[SeqSN];
    SeqMLen = SeqTB - SeqTA + 1;    /* max sequence length */
    SeqMDim = SeqMLen + 1;          /* dimension of D matrix */

    printf1("Number of states: %d. Max sequence length: %d\n",SeqNS,SeqMLen);

    err = seq_search();
   
SEQPMFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seq_search()    Search for patterns.                                    */
/*  #                                                                      */
/*  Return 0 if successful, otherwise -1.                                   */

int seq_search(void)
{
    register int i,j,k;
    int err,m,nn,nm;

    err = -1;
    printf1("\nSearching for patterns.\n");

    if (PMNPS <= 0) {
        p_err(-19,1);
        goto SSRCHFin;
    }

    SeqSPL = 0;                 /* max length of test pattern strings */
    for (i = 0; i < PMNPS; ++i) {
        printf1("Pattern%2d: ",i + 1);
        seq_pprn(i);
        if (SeqSPL < PMPSN[i])
            SeqSPL = PMPSN[i];
    }
    SeqSPL++;   
    if (seqm_alloc1(1)) {    /* allocate memory */
        p_err(-2,1);
        goto SSRCHFin;
    }
    if (PMF1Def) {      /* print to test output file */
        PMTST[2] = 1;
        fprintf(PMF1d,"Pattern search test output\n\n");
    }
    SeqNMD = 0;
    nn = m = 0;
    for (i = 0; i < NOC; ++i) {
        SeqSIL = seq_gst(i,1);
        if (SeqSIL < 1) {
            m++;
            continue;
        }
        fprintf(PMFd,PMNFmtS,i + 1,SeqSIL);
        fprintf(PMFd,PMNFmtS,SeqSIL);

        /* change to original state numbers */

        for (j = 0; j < SeqSIL; ++j)
            SeqSJ[j] = SeqSTNI[SeqSN][SeqSI[j]];
        SeqSJL = SeqSIL;

        seq_pps(1,0,1,i,0);               /* write to test output */

        for (k = 0; k < PMNPS; ++k) {
            nm = seq_sfind(k);
            fprintf(PMFd,PMNFmtS,nm);
   
            if (PMF1Def)            /* write to test output */
                seq_ppt(k);        
        }
        seq_pmd();                  /* print to stderr */
        seq_pvar(i,-1);             /* add variables */
        fprintf(PMFd,"\n");
        nn++;
    }
    printf1("\nNumber of sequences (cases): %d\n",NOC);
    printf1("Sequences with zero length or internal gaps: %d\n",m);
    printf1("%d records written to output file: %s\n",nn,PMFdName);

    if (PMTDAFDef && nn > 0)         /* write TDA description */
        seqpm_dtda(nn);

    if (PMF1Def)  
        printf1("Test output written to: %s\n",PMF1dName);

    err = 0;

SSRCHFin:
    seqm_alloc1(0);     /* free memory */
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seq_sfind(k)    Find test pattern k.                                    */
/*                                                                          */
/*  Return number of matches.                                               */

int seq_sfind(int k)
{
    register int i,j,l,ll,j1,jj;
    int m,a,b,bl,n,nn,nnf,nf,u,fndflg;
    char c;

    m = PMPSN[k];       /* length of test pattern */

    if (PMF1Def) {      /* write test pattern to test output */

        fprintf(PMF1d,"Pattern:");   
        for (i = 0; i < m; ++i) {
            n = (int)PMPS[k][i];
            if (n >= 0) 
                fprintf(PMF1d," %d",n);
            else {
                c = '\0';
                if (n == -1) 
                    c = '?';
                else if (n == -2) 
                    c = '+';
                else if (n == -3) 
                    c = '-';
                else if (n == -4) 
                    c = '*';
                fprintf(PMF1d," %c",c);
            }
        }
        fprintf(PMF1d,"\n");   
    }
    nf = 0;
    n = SeqSPL * SeqMDim;
    for (i = 1; i < n; ++i)
        SeqDP[i] = 0;
    SeqDP[0] = 1;

    j = 0; 
    for (i = 0; i < m; ++i) {   /* init first column */

        if (SeqDP[j] && (int)PMPS[k][i] <= -3) {
            j += SeqMDim;
            SeqDP[j] = 1;
        }
        else
            break;
    }
    bl = j1 = 0;
    while (j1 < SeqSJL) {

        fndflg = nn = nnf = 0;
        for (j = j1; j < SeqSJL; ++j) {

            l = j + 1;
            b = SeqSJ[j];
            n = 0;

            for (i = 0; i < m; ++i) {

                ll = l + SeqMDim;
                a = (int)PMPS[k][i];
                u = 0;
                if (SeqDP[l] && a <= -3)
                    u = 1; 
                else {
                    jj = 1;
                    if (j == j1)
                        jj += j1;

                    if (SeqDP[l - jj]) {
                        if ((a >= 0 && a != b) || (a == -2 && b != bl))
                            ;
                        else
                            u = 1;
                    }
                    if (u == 0 && SeqDP[ll - jj]) {
                        if (a == -4 || (a <= -2 && b == bl))
                            u = 1;
                    }
                }
                if (u) {
                    SeqDP[ll] = 1;
                    n++;
                }
                l = ll;
       
            }
            if (n == 0)  
                break;
                   
            if (SeqDP[m * SeqMDim + j + 1]) {
                fndflg = 1;
                nnf = nn;
            }
            nn++;
            bl = b;
        }
        if (fndflg) {
            nf++;
            if (PMF1Def) {      /* write match to test output */

                fprintf(PMF1d,"  Match:");
                for (l = j1; l <= j1 + nnf; ++l)
                    fprintf(PMF1d," %d",SeqSJ[l]);
                fprintf(PMF1d,"\n");
            }
        }
        if (nn > 0)
            j1 += nn;
        else
            j1++;
    }
    return(nf);
}

/*--------------------------------------------------------------------------*/
/*  seq_ppt(k)      Print DP matrix to test output file.                    */
/*                  For test pattern k.                                     */

void seq_ppt(int k)
{
    register int i,j;
    int n;
    char c;

    if (PMF1Def== 0)
        return;

    fprintf(PMF1d,"\nD Matrix\n");
    fprnchar(PMF1d,' ',9,0);
    for (j = 0; j < SeqSIL; ++j) {
        n = SeqSTNI[SeqSN][SeqSI[j]];
        fprintf(PMF1d,"%3d ",n);
    }
    fprintf(PMF1d,"\n    ");
    fprnchar(PMF1d,'-',4 + 4 * SeqSIL,1);

    for (i = 0; i < SeqSPL; ++i) {
        if (i == 0)
            fprintf(PMF1d,"    |");
        else {
            n = (int)PMPS[k][i - 1];
            if (n >= 0) 
                fprintf(PMF1d,"%3d |",n);
            else {
                c = '\0';
                if (n == -1) 
                    c = '?';
                else if (n == -2) 
                    c = '+';
                else if (n == -3) 
                    c = '-';
                else if (n == -4) 
                    c = '*';
                fprintf(PMF1d,"  %c |",c);
            }
        }
        for (j = 0; j <= SeqSIL; ++j)  
            fprintf(PMF1d,"%3d ",(int)SeqDP[i * SeqMDim + j]);
        fprintf(PMF1d,"\n");
    }
    fprintf(PMF1d,"\n");
}   

/*--------------------------------------------------------------------------*/
/*  seq_pprn(i)     Print test pattern i. Return 1 if pattern contains      */
/*                  special characters, otherwise 0.                        */

int seq_pprn(int i)
{
    register int j,k,n;
    char c;

    n = 0;
    for (j = 0; j < PMPSN[i]; ++j) {
        if (j)  
            printf1(",");
        k = (int)PMPS[i][j];
        if (k >= 0) 
            printf1("%d",k);
        else {
            c = '\0';
            if (k == -1) 
                c = '?';
            else if (k == -2) 
                c = '+';
            else if (k == -3) 
                c = '-';
            else if (k == -4) 
                c = '*';
            printf1("%c",c);
            n = 1;
        }
    }
    newline();
    return(n);
}

/*--------------------------------------------------------------------------*/
/*  seqm_alloc1(opt)    If opt != 0 allocate, else free memory.             */
/*                      Return 0 if OK, -1 if error.                        */

int seqm_alloc1(int opt)
{
    int n,err = 0;

    if (opt) {
        err = -1;
        n = SeqSPL * SeqMDim;
        if (!(SeqDP = (char *)calloc(n,sizeof(char))))   
            goto SEQMA1Fin;
        memrq(n,sizeof(char));
        SeqDPA = n;
    
        n = SeqMLen;
        if (!(SeqSI = (int *)calloc(n,sizeof(int))))   
            goto SEQMA1Fin;
        memrq(n,sizeof(int));
        SeqSIA = n;

        if (!(SeqSJ = (int *)calloc(n,sizeof(int))))   
            goto SEQMA1Fin;
        memrq(n,sizeof(int));
        SeqSJA = n;

        return(0);
    }

SEQMA1Fin:
    if (SeqDPA) {
        free((char *)SeqDP);
        memrq(-SeqDPA,sizeof(char));
        SeqDPA = 0;
    }
    if (SeqSIA) {
        free((char *)SeqSI);
        memrq(-SeqSIA,sizeof(int));
        SeqSIA = 0;
    }
    if (SeqSJA) {
        free((char *)SeqSJ);
        memrq(-SeqSJA,sizeof(int));
        SeqSJA = 0;
    }
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqpm_dtda(n)   Write TDA description file, n = number of cases.        */

void seqpm_dtda(int n)
{
    register int i,j,k;

    fprintf(PMTDAFd,"nvar(\n");
    fprintf(PMTDAFd,"  dfile = %s,\n",PMFdName);
    fprintf(PMTDAFd,"  noc = %d,\n",n);

    fprintf(PMTDAFd,"  SID  <5>[%d.0] = c1  , # sequence case number\n",PMNFmt);
    fprintf(PMTDAFd,"  SLEN <5>[%d.0] = c2  , # sequence length\n",PMNFmt);
    k = 3;
    for (i = 0; i < PMNPS; ++i) {
        fprintf(PMTDAFd,"  SM%-2d <5>[%d.0] = c%-3d, # number of matches, pattern %d\n",k-2,PMNFmt,k,k-2);
        k++;
    }
    for (i = 0; i < PMNV; ++i) {
        j = PMVIdx[i];
        fprintf(PMTDAFd,"  %s <%d>[%d.%d] = c%-3d,\n",       
                                    VName[j],VSLen[j],VPFmt1[j],VPFmt2[j],k++);
    }
    fprintf(PMTDAFd,");\n");
    printf1("TDA description written to: %s\n",PMTDAFName);
}

/*--##----------------------------------------------------------------------*/
/*  subm()      substitution metric for distributions                       */
/*                                                                          */
/*              subm(                                                       */
/*                  scost=...,     name of substitution cost matrix         */
/*                  eps=...,       epsilon for distribution check,          */
/*                                 def. 1.e-6                               */
/*                  fmt=...,       print format, def. 10.4                  */ 
/*                  df=...,        output file for substitution vector      */
/*              ) = X,Y1,...,Yk;                                            */
/*                                                                          */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int subm(void)
{
    register int i,j,k,l,ix,iy;
    int err,m,n,m1,n1,ni,nj,r,sm,smr,smc;
    double x,xsum,delta;

    err = -1;         
    if (check_cmd(0))
        return(-1);
         
    printf1("Substitution metric. Current memory: %d bytes.\n",MemReq);
         
    if (parm(CmdBuf + 4,4,1))    /* get parameters */
        goto SUBMFin;

    if (PMFmtF == 0)
        pmfmt(10,4);

    if (PMNV < 2) {
        printf1("Error: need at least two variables.\n");
        goto SUBMFin;
    }
    if (PMSCOSTMAT >= 0) {
        sm = PMSCOSTMAT;
        printf1("Substitution cost defined by matrix: %s\n",MatName[sm]);
        smr = MatRow[sm];
        smc = MatCol[sm];
        if (smr < NOC || smc < NOC) {
            printf1("Error: need at least an (%d,%d) matrix.\n",NOC,NOC);
            goto SUBMFin;
        }
        for (i = 0; i < NOC; ++i) {
            for (j = 1; j <= NOC; ++j) {
                if (MatVal[sm][i * smc + j] < 0.0) {
                    printf1("Error: substitution cost matrix contains negative values.\n");
                    goto SUBMFin;
                }
            }
        }
    }

    /* check distributions */

    n = 0;
    for (j = 0; j < PMNV; ++j) {
        ix = PMVIdx[j];
        xsum = 0.0;
        n = 0;
        for (i = 0; i < NOC; ++i) {
            x = get_data(ix,i);
            if (x < 0.0 || x > 1.0) {
                n = 1;
                break;
            }
            xsum += x;
        }
/*****  if (n || fabs(xsum - 1.0) > PMEPS) { *****/
        if (n) {
            printf("Error: no distribution in variable %s\n",VName[ix]);
            goto SUBMFin;
        }
    }   
    if (alloc_aci(NOC + 1))      
        goto SUBMFin;
    if (alloc_acj(NOC + 1))      
        goto SUBMFin;
    if (alloc_acu(NOC + 1))      
        goto SUBMFin;
    if (alloc_acv(NOC + 1))      
        goto SUBMFin;

    ix = PMVIdx[0];
    for (k = 1; k < PMNV; ++k) {
        delta = 0.0;

        ni = nj = 0;
        iy = PMVIdx[k];
        for (i = 0; i < NOC; ++i) {
            x = get_data(ix,i) - get_data(iy,i);
            if (x > PMEPS) {
                ni++;
                AcI[ni] = i+ 1;
                AcU[ni] = x;
            }
            else if (x < -PMEPS) {
                nj++;
                AcJ[nj] = i + 1;
                AcV[nj] = -x;
            }              
        }
        if (ni > 0 && nj > 0) {
            n = ni * nj;
            m = ni + nj;
            n1 = n + 1;
            m1 = m + 1;
            if (alloc_acz(n1 * m1 + 1))      
                goto SUBMFin;
            if (alloc_acx(n1 + 1))      
                goto SUBMFin;
            if (alloc_acy(m1 + 1))      
                goto SUBMFin;
    
            l = 1;
            for (i = 1; i <= ni; ++i) {
                for (j = 1; j <= nj; ++j)  
                    AcZ[l++] = -subm_cost(AcI[i],AcJ[j]);
            }
            l = 0;
            for (i = 1; i <= ni; ++i) {
                for (j = 1; j <= nj; ++j)  
                    AcZ[i * n1 + l + j] = 1.0;
                l += nj;   
            }
      
            for (i = 1; i <= nj; ++i) {
                for (j = 0; j < ni; ++j) {
                    AcZ[(ni + i) * n1 + i + j * nj] = 1.0;
                }
            }
            l = 2 * n1;
            for (i = 1; i <= ni; ++i) { 
                AcZ[l] = AcU[i]; 
                l += n1;
            }
            for (j = 1; j <= nj; ++j) {
                AcZ[l] = AcV[j];
                l += n1;
            }
            /********** 
            for (i = 1; i <= m1; ++i) {
                for (j = 1; j <= n1; ++j)  
                    printf("%lf ",AcZ[(i - 1) * n1 + j]);
                newline();
            }
            *****/
            if (PMF1Def) {
                for (i = 1; i <= NOC; ++i) {
                    for (j = 1; j <= NOC; ++j)  
                        fprintf(PMF1d,PMFmtS,subm_cost(i,j));
                    fprintf(PMF1d,"\n");
                }
                fprintf(PMF1d,"\n");
 
                for (i = 0; i < NOC; ++i) {
                    fprintf(PMF1d,"%4d ",i + 1);
                    fprintf(PMF1d,PMFmtS,get_data(ix,i));
                    fprintf(PMF1d,PMFmtS,get_data(iy,i));
                    fprintf(PMF1d,"\n");
                }                                      
                fprintf(PMF1d,"\n");

                for (i = 1; i <= ni; ++i) {
                    fprintf(PMF1d,"%4d %4d ",i,AcI[i]);
                    fprintf(PMF1d,PMFmtS,AcU[i]);
                    fprintf(PMF1d,"\n");
                }
                fprintf(PMF1d,"\n");
                for (j = 1; j <= nj; ++j) {
                    fprintf(PMF1d,"%4d %4d ",j,AcJ[j]);
                    fprintf(PMF1d,PMFmtS,AcV[j]);
                    fprintf(PMF1d,"\n");
                }
                fprintf(PMF1d,"\n");
                for (i = 1; i <= m1; ++i) {
                    for (j = 1; j <= n1; ++j)  
                        fprintf(PMF1d,PMFmtS,AcZ[(i - 1) * n1 + j]);
                    fprintf(PMF1d,"\n");
                }
                fprintf(PMF1d,"\n");
            }

            r = lpf1(m-1,n,m-1,AcZ,AcX,AcY,&delta,PMEPS);
            if (r) {
                switch (r) {
                    case -2:    printf1("No solution.\n");
                                break;
                    case -3:    printf1("No solution. Dual objective function is unbounded.\n");
                                break;
                    case -4:    printf1("No solution. Primal objective function is unbounded.\n");
                                break;
                    case -5:    printf1("Error: insufficient memory.\n");
                                break;
                    default:    printf1("Error (%d).\n",r);
                                break;
                }   
                err = 0;
                goto SUBMFin;
            }
            printf1("Distance: ");
            printf1(PMFmtS,-delta);
            newline();


            if (PMF1Def) {
                l = 1;
                for (i = 1; i <= ni; ++i) {
                    for (j = 1; j <= nj; ++j) {
                        fprintf(PMF1d,"%4d %4d ",AcI[i],AcJ[j]);
                        fprintf(PMF1d,PMFmtS,AcX[l++]);
                        fprintf(PMF1d,"\n");
                    }
                }
                printf1("\nSubstitution vector written to: %s\n",PMF1dName);
            }
        }
    }
    err = 0;
   
SUBMFin:
    p_clean();
    return(err);
}

double subm_cost(int i,int j)    
{
    int sm,smc;
    double c;

    if (PMSCOSTMAT >= 0) { 
        sm = PMSCOSTMAT;
        smc = MatCol[sm];
        c = MatVal[sm][(i - 1) * smc + j];
    }
    else
        c = fabs((double)(i - j));
    return(c);
}


