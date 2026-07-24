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
#include "tda_context.h"

/*  functions in t_seqm.c */

int seqm(TDAContext *ctx);
int get_costf(TDAContext *ctx, int opt);
int seqm_init(TDAContext *ctx, int opt);
void seq_pvar(TDAContext *ctx, int i,int j);
int seq_gst(TDAContext *ctx, int i,int opt);
void seq_pps(TDAContext *ctx, int i,int j,int opt,int ic1,int ic2);
void seq_pmd(TDAContext *ctx);
void seq_plcs(TDAContext *ctx);
void seq_pseq(TDAContext *ctx);
void seq_dmat(TDAContext *ctx);
int seq_pmatch(TDAContext *ctx);
float seq_mdist(TDAContext *ctx);
void seqm_pcost(TDAContext *ctx);
void seqm_dtda(TDAContext *ctx, int n);
int seq_scost(TDAContext *ctx);
int seqpm(TDAContext *ctx);
int seq_search(TDAContext *ctx);
int seq_sfind(TDAContext *ctx, int k);
void seq_ppt(TDAContext *ctx, int k);
int seq_pprn(TDAContext *ctx, int i);
int seqm_alloc1(TDAContext *ctx, int opt);
void seqpm_dtda(TDAContext *ctx, int n);
int subm(TDAContext *ctx);
double subm_cost(TDAContext *ctx, int i,int j);   

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */








/*--------------------------------------------------------------------------*/
/*  seqm()      Optimal matching                                            */
/*              Return 0 if OK, -1 if error.                                */

int seqm(TDAContext *ctx)
{
    register int i;
    int err,n,m;

    err = -1;         
    if (check_cmd(ctx, 0))
        return(-1);
         
    printf1(ctx, "Sequence proximity measures. Current memory: %d bytes.\n",ctx->MemReq);
         
    if (parm(ctx, ctx->CmdBuf + 4,1,1))    /* get parameters */
        goto SEQMFin;

    if (ctx->PMTFmtF == 0)
        pmtfmt(ctx, 5,2);

    if (ctx->PMM < 1)
        ctx->PMM = 1;

    switch (ctx->PMM) {

        case 1: printf1(ctx, "Optimal matching.\n");
                ctx->SeqSN = seq_getsn(ctx, ctx->PMSN,1);   /* get sequence number */
                if (ctx->SeqSN < 0)  
                    goto SEQMFin;

                ctx->SeqSN1 = ctx->SeqSN2 = ctx->SeqSN;
                ctx->SeqTA = ctx->SeqTMin[ctx->SeqSN];     /* range of time axis */
                ctx->SeqTB = ctx->SeqTMax[ctx->SeqSN];
                break;

        case 2: printf1(ctx, "Comparing parallel sequences.\n");
                if (ctx->PMSN < 1 || ctx->PMSN1 < 1 ||
                         (ctx->SeqSN1 = seq_getsn(ctx, ctx->PMSN,0)) < 0 ||
                                  (ctx->SeqSN2 = seq_getsn(ctx, ctx->PMSN1,0)) < 0) { 

                    printf1(ctx, "Error: need two valid sequence numbers.\n");
                    goto SEQMFin;
                }
                printf1(ctx, "Using sequence data structures: %d and %d.\n",ctx->PMSN,ctx->PMSN1);
                if (ctx->SeqSN1 == ctx->SeqSN2) {
                    printf1(ctx, "Error: should be different sequences.\n");
                    goto SEQMFin;
                }
                ctx->SeqSN = ctx->SeqSN1;                     

                if (seq_scheck(ctx, ctx->SeqSN1,ctx->SeqSN2)) {
                    printf1(ctx, "Error: must have identical state space.\n");
                    goto SEQMFin;
                }
                ctx->SeqTA = imin(ctx, ctx->SeqTMin[ctx->SeqSN1],ctx->SeqTMin[ctx->SeqSN2]);   /* range of time axis */
                ctx->SeqTB = imax(ctx, ctx->SeqTMax[ctx->SeqSN1],ctx->SeqTMax[ctx->SeqSN2]);    

                break;


        default:  printf1(ctx, "Error: method %d not available.\n",ctx->PMM);
                  goto SEQMFin;
    }
    ctx->SeqNS = ctx->SeqSTN[ctx->SeqSN];          /* number of states */
    ctx->SeqMLen = ctx->SeqTB - ctx->SeqTA + 1;    /* max sequence length */
    ctx->SeqMDim = ctx->SeqMLen + 1;          /* dimension of D matrix */

    printf1(ctx, "Number of states: %d. Max sequence length: %d\n",ctx->SeqNS,ctx->SeqMLen);
    if (ctx->PMSM[1])  
        printf1(ctx, "Option (sm=1): skip internal gaps.\n");
    if (ctx->PMSM[2])  
        printf1(ctx, "Option (sm=2): skip identical states.\n");
    if (ctx->PMRRN == 1)
        printf1(ctx, "Using common sequence length (rr=1).\n");

    if (ctx->PMF1Def)
        printf1(ctx, "Test output will be written to: %s\n",ctx->PMF1dName);
    newline(ctx);

    if (get_costf(ctx, 1))       /* set up cost functions */
        goto SEQMFin;

    if (ctx->PMF1Def)            /* print to test output file */
        seqm_pcost(ctx);

    if (ctx->ICostTyp == 1)      /* for linear indel cost functions */
        ctx->SeqMTyp = 1;
    else
        ctx->SeqMTyp = 0;

    if (ctx->PMMax >= 1 && ctx->SeqMTyp == 0)  
        printf1(ctx, "Restricted alignment: max=%d\n",ctx->PMMax);
    else
        ctx->PMMax = 0;

    if (seqm_init(ctx, 1))       /* init D, E, F matrices */
        goto SEQMFin;

    /* check cn option */

    if (ctx->PMM == 1 && ctx->PMNCN > 0) {

        printf1(ctx, "Checking cn option.\n");
        m = 0;
        for (i = 0; i < ctx->PMNCN; ++i) {
            if (ctx->PMCN[i] < 1 || ctx->PMCN[i] > ctx->NOC) {
                printf1(ctx, "Error: cn parameter inconsistent with number of cases (= %d).\n",ctx->NOC);
                goto SEQMFin;
            }
            n = ctx->PMCN[i];
            if (seq_gst(ctx, n - 1,1) < 1) {
                printf1(ctx, "Warning: sequence %d has zero length and will be ignored.\n",n);
                ctx->PMCN[i] = -1;
            }
            else 
                m++;
        }
        printf1(ctx, "Will use %d reference sequences.\n",m);
        if (m == 0)
            goto SEQMFin;
    }
    n = seq_pmatch(ctx);               /* alignment */

    if (ctx->PMTDAFDef && n > 0)         /* write TDA description */
        seqm_dtda(ctx, n);

    err = 0;
   
SEQMFin:
    seqm_init(ctx, 0);
    get_costf(ctx, 0);
    p_clean(ctx);
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

int get_costf(TDAContext *ctx, int opt)
{
    register int i = 0,j = 0;
    int err = 0,ic = 0,sc = 0,n = 0,im = 0,imr = 0,imc = 0,sm = 0,smr = 0,smc = 0;

    if (opt == 0) {
        err = 0;
        goto COSTFFin;
    }
    err = -2;

    n = ctx->SeqMLen;
    if (!(ctx->SeqSI = (int *)calloc((size_t)(n),sizeof(int))))   
        goto COSTFFin;
    memrq(ctx, n,sizeof(int));
    ctx->SeqSIA = n;

    if (!(ctx->SeqSJ = (int *)calloc((size_t)(n),sizeof(int))))   
        goto COSTFFin;
    memrq(ctx, n,sizeof(int));
    ctx->SeqSJA = n;

    if (!(ctx->SeqSS = (int *)calloc((size_t)(n),sizeof(int))))   
        goto COSTFFin;
    memrq(ctx, n,sizeof(int));
    ctx->SeqSSA = n;

    if (!(ctx->ICost = (float *)calloc((size_t)(n),sizeof(float))))   
        goto COSTFFin;
    memrq(ctx, n,sizeof(float));
    ctx->ICostA = n;

    n = ctx->SeqNS * ctx->SeqNS;  
    if (!(ctx->SCost = (float *)calloc((size_t)(n),sizeof(float))))   
        goto COSTFFin;
    memrq(ctx, n,sizeof(float));
    ctx->SCostA = n;
            
    err = -1;
    ic = sc = 0;

    /* indel cost */

    if (ctx->PMICOSTMAT >= 0) {
        ctx->ICostTyp = 99;
        im = ctx->PMICOSTMAT;
        printf1(ctx, "Indel cost defined by matrix: %s\n",ctx->MatName[im]);
        imr = ctx->MatRow[im];
        imc = ctx->MatCol[im];
        if ((imr != 1 || imc != ctx->SeqMLen) && (imr != ctx->SeqMLen || imc != 1)) {
            printf1(ctx, "Error: need an (1,%d) or (%d,1) matrix.\n",ctx->SeqMLen,ctx->SeqMLen);
            goto COSTFFin;
        }
    }
    else if (ctx->PMICOSTA >= 0.0 && ctx->PMICOSTB >= 0.0) {
        ctx->ICostTyp = 1;
        ctx->ICAlpha = (float)(ctx->PMICOSTA);
        ctx->ICBeta = (float)(ctx->PMICOSTB);
        printf1(ctx, "Linear indel cost function: alpha=%g, beta=%g.\n",(double)(ctx->ICAlpha),(double)(ctx->ICBeta));
    }
    else if (ctx->PMICOSTA >= 0.0) {
        ctx->ICostTyp = 0;
        printf1(ctx, "Indel cost: %g.\n",ctx->PMICOSTA);
    }
    else {
        ic = 1;
        ctx->ICostTyp = 0;
        ctx->PMICOSTA = 1.0;
        printf1(ctx, "Default indel cost: 1.\n");
    }

    /* substitution cost */

    if (ctx->PMSCOSTMAT >= 0) {
        ctx->SCostTyp = 99;
        sm = ctx->PMSCOSTMAT;
        printf1(ctx, "Substitution cost defined by matrix: %s\n",ctx->MatName[sm]);
        smr = ctx->MatRow[sm];
        smc = ctx->MatCol[sm];
        if (smr < ctx->SeqNS || smc < ctx->SeqNS) {
            printf1(ctx, "Error: need at least an (%d,%d) matrix.\n",ctx->SeqNS,ctx->SeqNS);
            goto COSTFFin;
        }
    }
    else if (ctx->PMSCOSTM == 1) {
        ctx->SCostTyp = 1;
        printf1(ctx, "Substitution cost defined by absolute difference.\n");
    }
    else if (ctx->PMSCOSTM == 2) {
        ctx->SCostTyp = ctx->PMSCOSTM;
        printf1(ctx, "Substitution cost based on data, type %d.\n",ctx->SCostTyp);
    }
    else {
        sc = 1;
        ctx->SCostTyp = 0;
        printf1(ctx, "Default substitution cost: 2.\n");
    }
    if (ic == 0 || sc == 0) {     /* icost or scost parameters used */
        if (ctx->PMS > 1)
            ctx->PMS = 0;
    }

    /* set up indel cost in ICost[] */

    for (i = 0; i < ctx->SeqMLen; ++i)
        if (ctx->ICostTyp == 0)
            ctx->ICost[i] = (float)(ctx->PMICOSTA);
        else if (ctx->ICostTyp == 99) {
            ctx->ICost[i] = (float)ctx->MatVal[im][i + 1];
    }
    if (ctx->ICostTyp == 99)
        ctx->ICostTyp = 0;

    /* set up substitution cost in SCost[] */

    if (ctx->SCostTyp == 2) {
        seq_scost(ctx);
        ctx->SCostTyp = 0;
    }
    else {
        for (i = 0; i < ctx->SeqNS; ++i) {
            for (j = 0; j < ctx->SeqNS; ++j) {
                if (j != i) {
                    if (ctx->SCostTyp == 0)  
                        ctx->SCost[i * ctx->SeqNS + j] = 2.0; 
                    else if (ctx->SCostTyp == 99)
                        ctx->SCost[i * ctx->SeqNS + j] = (float)(ctx->MatVal[sm][i * smc + j + 1]);
                }
            }
        }
        if (ctx->SCostTyp == 99)
            ctx->SCostTyp = 0;
    }
    return(0);

COSTFFin:
    if (err == -2)  
        p_err(ctx, -2,1);

    if (ctx->SeqSIA > 0) {
        free((char *)ctx->SeqSI);
        memrq(ctx, -ctx->SeqSIA,sizeof(int));
        ctx->SeqSIA = 0;
    }
    if (ctx->SeqSJA > 0) {
        free((char *)ctx->SeqSJ);
        memrq(ctx, -ctx->SeqSJA,sizeof(int));
        ctx->SeqSJA = 0;
    }
    if (ctx->SeqSSA > 0) {
        free((char *)ctx->SeqSS);
        memrq(ctx, -ctx->SeqSSA,sizeof(int));
        ctx->SeqSSA = 0;
    }
    if (ctx->ICostA > 0) {
        free((char *)ctx->ICost);
        memrq(ctx, -ctx->ICostA,sizeof(float));
        ctx->ICostA = 0;         
    }
    if (ctx->SCostA > 0) {
        free((char *)ctx->SCost);
        memrq(ctx, -ctx->SCostA,sizeof(float));
        ctx->SCostA = 0;         
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

int seqm_init(TDAContext *ctx, int opt)
{
    register int i;
    int err,n;
    double tmp;

    if (opt == 0) {
        err = 0;
        goto SEQMIFin;
    }
    err = -1;

    n = ctx->SeqMDim * ctx->SeqMDim;
    if (!(ctx->SeqMD = (float *)calloc((size_t)(n),sizeof(float))))   
        goto SEQMIFin;
    memrq(ctx, n,sizeof(float));
    ctx->SeqMDA = n;

    if (ctx->SeqMTyp == 1) {
        if (!(ctx->SeqME = (float *)calloc((size_t)(n),sizeof(float))))   
            goto SEQMIFin;
        memrq(ctx, n,sizeof(float));
        ctx->SeqMEA = n;

        if (!(ctx->SeqMF = (float *)calloc((size_t)(n),sizeof(float))))   
            goto SEQMIFin;
        memrq(ctx, n,sizeof(float));
        ctx->SeqMFA = n;
    }
    ctx->SeqMD[0] = 0.0;                 /* init SeqMD */

    if (ctx->SeqMTyp == 0) {
        for (i = 0; i < ctx->SeqMLen; ++i) {
            ctx->SeqMD[i + 1] = ctx->SeqMD[i] + ctx->ICost[i];
            ctx->SeqMD[(i + 1) * ctx->SeqMDim] = ctx->SeqMD[i * ctx->SeqMDim] + ctx->ICost[i];
        }
    }
    else if (ctx->SeqMTyp == 1) {        /* also SeqME and SeqMF */

        ctx->SeqME[0] = ctx->SeqMF[0] = 0.0;

        for (i = 0; i < ctx->SeqMLen; ++i) {
            tmp = (double)((float)ctx->ICAlpha + (float)ctx->ICBeta * (float)i);
            ctx->SeqMD[(i + 1) * ctx->SeqMDim] = (float)tmp;
            ctx->SeqMD[i + 1] = (float)tmp;
            ctx->SeqMF[i + 1] =                                   
            ctx->SeqME[(i + 1) * ctx->SeqMDim] = (float)((double)(tmp) + (double)((float)ctx->ICAlpha));
        }
    }
    if (ctx->PMF1Def && ctx->PMTST[1]) {
        fprintf(ctx->PMF1d,"Initial D matrix\n");
        mprf(ctx, ctx->SeqMDim,ctx->SeqMDim,ctx->SeqMD - 1,ctx->PMFmtS,ctx->PMF1d);

        if (ctx->SeqMTyp == 1) {        /* also SeqME and SeqMF */

            fprintf(ctx->PMF1d,"\nInitial E matrix\n");
            mprf(ctx, ctx->SeqMDim,ctx->SeqMDim,ctx->SeqME - 1,ctx->PMFmtS,ctx->PMF1d);

            fprintf(ctx->PMF1d,"\nInitial F matrix\n");
            mprf(ctx, ctx->SeqMDim,ctx->SeqMDim,ctx->SeqMF - 1,ctx->PMFmtS,ctx->PMF1d);
        }
        fprintf(ctx->PMF1d,"\n");
    }
    return(0);

SEQMIFin:
    if (err)  
        p_err(ctx, -2,1);

    if (ctx->SeqMDA > 0) {
        free((char *)ctx->SeqMD);
        memrq(ctx, -ctx->SeqMDA,sizeof(float));
        ctx->SeqMDA = 0;
    }
    if (ctx->SeqMEA > 0) {
        free((char *)ctx->SeqME);
        memrq(ctx, -ctx->SeqMEA,sizeof(float));
        ctx->SeqMEA = 0;
    }
    if (ctx->SeqMFA > 0) {
        free((char *)ctx->SeqMF);
        memrq(ctx, -ctx->SeqMFA,sizeof(float));
        ctx->SeqMFA = 0;
    }
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seq_pvar(i,j) Print additional variables to output file, for cases i,j  */

void seq_pvar(TDAContext *ctx, int i,int j)
{
    register int k,iv;

    for (k = 0; k < ctx->PMNV; ++k) {
        iv = (int)ctx->PMVIdx[k];
        rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[iv],get_data(ctx, iv,i));
        if (j >= 0)
            rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[iv],get_data(ctx, iv,j));
#ifdef TDA_R_PACKAGE
        /* v= columns are written at each variable's own print format,
           so they need a producer like everything else; the caller
           names the row being built in ctx->RExportRowKey */
        if (ctx->RExportRowKey[0]) {
            tda_export_cell(ctx, ctx->RExportRowKey, get_data(ctx, iv,i));
            if (j >= 0)
                tda_export_cell(ctx, ctx->RExportRowKey,
                                get_data(ctx, iv,j));
        }
#endif
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

int seq_gst(TDAContext *ctx, int i,int opt)
{
    register int t,l,a,b,n,nn;

    a = -1;
    for (t = ctx->SeqTA; t <= ctx->SeqTB; ++t) {
        n = seq_sget(ctx, i,t,ctx->SeqSN);
        if (n >= 0) {
            a = t;
            break;
        }
    }
    if (a < 0)
        return(0);

    b = -1;
    for (t = ctx->SeqTB; t >= a; --t) {
        n = seq_sget(ctx, i,t,ctx->SeqSN);
        if (n >= 0) {
            b = t;
            break;
        }
    }
    nn = l = 0;  
    for (t = a; t <= b; ++t) {

        n = seq_sget(ctx, i,t,ctx->SeqSN);
        if (n < 0) {
            if (ctx->PMSM[1])
                continue;
            return(-1);
        }
        if (ctx->PMSM[2] && l && n == nn)
            continue;
        nn = n;
        n = ctx->SeqSTNII[ctx->SeqSN][n];
        if (opt == 1)
            ctx->SeqSI[l++] = n;
        else
            ctx->SeqSJ[l++] = n;
    }
    return(l);
}

/*--------------------------------------------------------------------------*/
/*  seq_pps(i,j)    Print sequences to test output file. If i != 0, print   */
/*                  first sequence (SeqSI), if j != 0, print second         */
/*                  sequence (SeqSJ).                                       */

void seq_pps(TDAContext *ctx, int i,int j,int opt,int ic1,int ic2)
{
    register int k;

    if (ctx->PMF1Def && ctx->PMTST[2]) {  
        if (i) {
            fprintf(ctx->PMF1d,"Sequence A (case number %d)\n",ic1 + 1);
#ifdef TDA_R_PACKAGE
            /* the pair's own identity and the two state sequences, so
               the dp-matrix blocks can be assembled without reading
               the test-output file back */
            tda_export_cell(ctx, "seqm.dp.pair", (double)(ic1 + 1));
#endif
            for (k = 0; k < ctx->SeqSIL; ++k) { 
                fprintf(ctx->PMF1d," %2d",ctx->SeqSTNI[ctx->SeqSN1][ctx->SeqSI[k]]);
#ifdef TDA_R_PACKAGE
                tda_export_cell(ctx, "seqm.dp.seqA",
                                (double)ctx->SeqSTNI[ctx->SeqSN1][ctx->SeqSI[k]]);
#endif
            }
            fprintf(ctx->PMF1d,"\n");
#ifdef TDA_R_PACKAGE
            tda_export_endrow(ctx, "seqm.dp.seqA");
            tda_export_flush(ctx, "seqm.dp.seqA");
#endif
        }
        if (j) {
            fprintf(ctx->PMF1d,"Sequence B (case number %d)\n",ic2 + 1);
#ifdef TDA_R_PACKAGE
            tda_export_cell(ctx, "seqm.dp.pair", (double)(ic2 + 1));
            tda_export_endrow(ctx, "seqm.dp.pair");
#endif
            for (k = 0; k < ctx->SeqSJL; ++k) { 
                fprintf(ctx->PMF1d," %2d",ctx->SeqSTNI[ctx->SeqSN2][ctx->SeqSJ[k]]);
#ifdef TDA_R_PACKAGE
                tda_export_cell(ctx, "seqm.dp.seqB",
                                (double)ctx->SeqSTNI[ctx->SeqSN2][ctx->SeqSJ[k]]);
#endif
            }
            fprintf(ctx->PMF1d,"\n");
#ifdef TDA_R_PACKAGE
            tda_export_endrow(ctx, "seqm.dp.seqB");
            tda_export_flush(ctx, "seqm.dp.seqB");
#endif
        }   
        if (opt)
            fprintf(ctx->PMF1d,"\n");    
    }
}
   
/*--------------------------------------------------------------------------*/
/*  seq_pmd()   Print number of alignments to stderr.                       */

void seq_pmd(TDAContext *ctx)
{
    if (ctx->SILENTFlg >= 2)
        return;

    ctx->SeqNMD++;
    if ((ctx->SeqNMD / 100) * 100 == ctx->SeqNMD)
        printfe(ctx, "Done: %7d     %c",ctx->SeqNMD,CR);
}

/*--------------------------------------------------------------------------*/
/*  seq_plcs()      Print longest common subsequence.                       */

void seq_plcs(TDAContext *ctx)
{
    register int i,j,k,ii;
    int dij;

    k = 0;
    i = ctx->SeqSIL - 1;
    j = ctx->SeqSJL - 1;
    while (i >= 0 && j >= 0) {
        ii = (i + 1) * ctx->SeqMDim;
        dij = (int)ctx->SeqMD[ii + j + 1];
        if (dij == (int)ctx->SeqMD[i * ctx->SeqMDim + j + 1] + (int)ctx->ICost[i])
            i--;
        else if (dij == (int)ctx->SeqMD[ii + j] + (int)ctx->ICost[j])
            j--;
        else {
            ctx->SeqSS[k++] = ctx->SeqSI[i];
            i--;
            j--;
        }
    }
    fprintf(ctx->PMFd,"%4d ",k);    
#ifdef TDA_R_PACKAGE
    if (ctx->RExportRowKey[0])
        tda_export_cell(ctx, ctx->RExportRowKey, (double)k);
#endif
  
    for (i = k - 1; i >= 0; --i) {
        fprintf(ctx->PMFd,"%2d ",ctx->SeqSTNI[ctx->SeqSN][ctx->SeqSS[i]]);
#ifdef TDA_R_PACKAGE
        if (ctx->RExportRowKey[0])
            tda_export_cell(ctx, ctx->RExportRowKey,
                            (double)ctx->SeqSTNI[ctx->SeqSN][ctx->SeqSS[i]]);
#endif
    }
    for (i = k; i < ctx->SeqMLen; ++i) {
        fprintf(ctx->PMFd,"-1 ");
#ifdef TDA_R_PACKAGE
        if (ctx->RExportRowKey[0])
            tda_export_cell(ctx, ctx->RExportRowKey, -1.0);
#endif
    }
}

/*--------------------------------------------------------------------------*/
/*  seq_pseq()      Print sequential distances from SeqMD.                  */

void seq_pseq(TDAContext *ctx)
{
    register int i,k;

    k = imin(ctx, ctx->SeqSIL,ctx->SeqSJL);
    for (i = 1; i <= k; ++i) {  
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMTFmtS,(double)ctx->SeqMD[i * ctx->SeqMDim + i]);
#ifdef TDA_R_PACKAGE
        if (ctx->RExportRowKey[0])
            tda_export_cell(ctx, ctx->RExportRowKey,
                            (double)ctx->SeqMD[i * ctx->SeqMDim + i]);
#endif
    }
       
    if (k < ctx->SeqSIL) {
        for (i = k + 1; i <= ctx->SeqSIL; ++i) { 
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMTFmtS,(double)ctx->SeqMD[i * ctx->SeqMDim + ctx->SeqSJL]);
#ifdef TDA_R_PACKAGE
            if (ctx->RExportRowKey[0])
                tda_export_cell(ctx, ctx->RExportRowKey,
                                (double)ctx->SeqMD[i * ctx->SeqMDim + ctx->SeqSJL]);
#endif
        }
    }
    else if (k < ctx->SeqSJL) {
        for (i = k + 1; i <= ctx->SeqSJL; ++i) { 
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMTFmtS,(double)ctx->SeqMD[ctx->SeqSIL * ctx->SeqMDim + i]);
#ifdef TDA_R_PACKAGE
            if (ctx->RExportRowKey[0])
                tda_export_cell(ctx, ctx->RExportRowKey,
                                (double)ctx->SeqMD[ctx->SeqSIL * ctx->SeqMDim + i]);
#endif
        }
    }
    k = imax(ctx, ctx->SeqSIL,ctx->SeqSJL);
    while (k++ < ctx->SeqMLen) {
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMTFmtS,-1.0);
#ifdef TDA_R_PACKAGE
        if (ctx->RExportRowKey[0])
            tda_export_cell(ctx, ctx->RExportRowKey, -1.0);
#endif
    }
}

/*--------------------------------------------------------------------------*/
/*  seq_dmat()          Print D matrix to test output file.                 */
/*                      If SeqMTyp = 1 : also E and F matrices              */

void seq_dmat(TDAContext *ctx)
{
    register int i,j,im;
    int is,js,n;
    double tmp = 0.0;

    for (im = 0; im < 3; ++im) {
        if (im == 0)  
            fprintf(ctx->PMF1d,"D");
        else if (im == 1)
            fprintf(ctx->PMF1d,"E");
        else if (im == 2)
            fprintf(ctx->PMF1d,"F");

        fprintf(ctx->PMF1d," Matrix  ");

        /* fprnchar(PMF1d,' ',10,0); */ 
        for (j = 0; j <= ctx->SeqSJL; ++j)  
            fprintf(ctx->PMF1d,"%7d ",j);
        fprintf(ctx->PMF1d,"\n         ");
        fprnchar(ctx, ctx->PMF1d,'-',8 + 8 * ctx->SeqSJL,1);
        fprintf(ctx->PMF1d,"                B ");
        for (j = 0; j < ctx->SeqSJL; ++j) {
            js = ctx->SeqSJ[j];
            n = ctx->SeqSTNI[ctx->SeqSN][js];
            fprintf(ctx->PMF1d,"%7d ",n);
        }
        fprintf(ctx->PMF1d,"\n         ");
        fprnchar(ctx, ctx->PMF1d,'-',8 + 8 * ctx->SeqSJL,1);
    
        for (i = 0; i <= ctx->SeqSIL; ++i) {
            if (i == 0)
                fprintf(ctx->PMF1d,"%4d   A |",i );
            else {
                is = ctx->SeqSI[i - 1];
                n = ctx->SeqSTNI[ctx->SeqSN][is];
                fprintf(ctx->PMF1d,"%4d %3d |",i,n);
            }
            for (j = 0; j <= ctx->SeqSJL; ++j) {

                if (im == 0)
                    tmp = (double)(ctx->SeqMD[i * ctx->SeqMDim + j]);
                else if (im == 1)
                    tmp = (double)(ctx->SeqME[i * ctx->SeqMDim + j]);
                else if (im == 2)
                    tmp = (double)(ctx->SeqMF[i * ctx->SeqMDim + j]);

                fprintf(ctx->PMF1d,"%7.2f ",tmp);
#ifdef TDA_R_PACKAGE
                tda_export_cell(ctx, im == 0 ? "seqm.dp.D" :
                                     im == 1 ? "seqm.dp.E" : "seqm.dp.F",
                                (double)tmp);
#endif
            }
            fprintf(ctx->PMF1d,"\n");
#ifdef TDA_R_PACKAGE
            tda_export_endrow(ctx, im == 0 ? "seqm.dp.D" :
                                   im == 1 ? "seqm.dp.E" : "seqm.dp.F");
#endif
        }
#ifdef TDA_R_PACKAGE
        /* one flush per printed matrix, so repeated pairs come back as
           seqm.dp.D, seqm.dp.D.2, ... beside their own blocks */
        tda_export_flush(ctx, im == 0 ? "seqm.dp.D" :
                              im == 1 ? "seqm.dp.E" : "seqm.dp.F");
#endif
        fprintf(ctx->PMF1d,"\n");
        if (ctx->SeqMTyp == 0)
            break;
    }
}   
    
/*--------------------------------------------------------------------------*/
/*  seq_pmatch()    Optimal matching, for options:                          */
/*                  PMM = 1 : optimal matching                              */
/*                  PMM = 2 : compare parallel sequences                    */
/*                                                                          */
/*  Return number of records written to output file.                        */

int seq_pmatch(TDAContext *ctx)
{
    register int i,j,k;
#ifdef TDA_R_PACKAGE
    /* seq_pseq(), seq_plcs() and seq_pvar() all add values to this same
       row; naming it here is what routes their cells to it */
    snprintf(ctx->RExportRowKey, sizeof(ctx->RExportRowKey), "seqm.pairs");
#endif
    int nm,m,nn,id,jd,sil,sjl;
    double tmp,dmax;
    double r,rdh = 0.0;

    printf1(ctx, "\nStarting alignment procedure.\n");

    ctx->SeqNMD = 0;         /* number of alignments */
    nn = m = 0;
    dmax = 0.0;
    id = jd = 0;

    if (ctx->PMM == 1) {

        if (ctx->PMNCN <= 0) {           /* pairwise comparison */

            nm = ctx->NOC * (ctx->NOC - 1) / 2;
            if (ctx->PMR > 0) {
                printf1(ctx, "Random selection of approximately %d sequence pairs.\n",ctx->PMR);
                if (ctx->PMR >= nm) {
                    printf1(ctx, "Will be ignored.\n");
                    ctx->PMR = 0;
                }
                else {
                    rdh = (double)ctx->PMR / (double)nm;
                    nm = ctx->PMR;
                }
            }
            if (ctx->SILENTFlg < 2)
                printfe(ctx, "Total alignments: %d\n",nm);

            for (i = 0; i < ctx->NOC; ++i) {
                sil = seq_gst(ctx, i,1);
                if (sil < 1) {
                    m++;
                    continue;
                }
                for (j = 0; j < i; ++j) {
                    sjl = seq_gst(ctx, j,2);
                    if (sjl < 1)
                        continue;

                    if (ctx->PMRRN == 1)
                        ctx->SeqSIL = ctx->SeqSJL = imin(ctx, sil,sjl);
                    else {
                        ctx->SeqSIL = sil;
                        ctx->SeqSJL = sjl;
                    }

                    if (ctx->PMR > 0) {
                        r = random1(ctx);
                        if (r > rdh)
                            continue;
                    }
                    seq_pps(ctx, 1,1,1,i,j);         /* write to test output */

                    tmp = (double)(seq_mdist(ctx));          /* calculate distance */
               
                    seq_pmd(ctx);                  /* print to stderr */

                    if (dmax < tmp) {
                        dmax = tmp;
                        id = i;
                        jd = j;
                    }
                    fprintf(ctx->PMFd,"%6d %6d %6d %6d ",i + 1,j + 1,ctx->SeqSIL,ctx->SeqSJL);
#ifdef TDA_R_PACKAGE
                    /* every value this row writes, in the order it is
                       written: the four identifiers, then whatever the
                       s= option puts in the distance position, then
                       v='s own columns from seq_pvar() */
                    tda_export_cell(ctx, "seqm.pairs", (double)(i + 1));
                    tda_export_cell(ctx, "seqm.pairs", (double)(j + 1));
                    tda_export_cell(ctx, "seqm.pairs", (double)ctx->SeqSIL);
                    tda_export_cell(ctx, "seqm.pairs", (double)ctx->SeqSJL);
#endif

                    if (ctx->PMS == 1)           /* sequential */
                        seq_pseq(ctx);
                    else {   
                        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMTFmtS,(double)tmp);
#ifdef TDA_R_PACKAGE
                        tda_export_cell(ctx, "seqm.pairs", (double)tmp);
#endif
                        if (ctx->PMS == 2)      /* LCS */
                            seq_plcs(ctx);
                    }           
                    seq_pvar(ctx, i,j);          /* print additional variables */
                    fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
                    tda_export_endrow(ctx, "seqm.pairs");
#endif
                    nn++;
                }
            }
        }
        else {              /* compare with PMCN[] sequences */

            for (i = 0; i < ctx->NOC; ++i) {
                sil = seq_gst(ctx, i,1);

                if (sil < 1) {
                    m++;
                    continue;
                }
                fprintf(ctx->PMFd,"%5d %3d ",i + 1,sil);

                for (j = 0; j < ctx->PMNCN; ++j) {
                    k = ctx->PMCN[j] - 1;
                    if (k >= 0) {
                        sjl = seq_gst(ctx, k,2);
                        if (sjl < 1)
                            continue;

                        if (ctx->PMRRN == 1)
                            ctx->SeqSIL = ctx->SeqSJL = imin(ctx, sil,sjl);
                        else {
                            ctx->SeqSIL = sil;
                            ctx->SeqSJL = sjl;
                        }
                        seq_pps(ctx, 1,1,1,i,k);         /* write to test output */

                        tmp = (double)(seq_mdist(ctx));    
               
                        seq_pmd(ctx);                  /* print to stderr */
                     
                        if (ctx->PMS == 1)               /* sequential */
                            seq_pseq(ctx);
                        else {
                            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMTFmtS,(double)tmp);
                            if (ctx->PMS == 2)          /* LCS */
                                seq_plcs(ctx);
                        }           
                    }
                }
                seq_pvar(ctx, i,-1);
                fprintf(ctx->PMFd,"\n");
                nn++;
            }
        }           
    }
    else if (ctx->PMM == 2) {      /* compare two parallel sequences */

        for (i = 0; i < ctx->NOC; ++i) {
            ctx->SeqSN = ctx->SeqSN1;
            ctx->SeqSIL = seq_gst(ctx, i,1);
            ctx->SeqSN = ctx->SeqSN2;
            ctx->SeqSJL = seq_gst(ctx, i,2);

            if (ctx->SeqSIL < 1 || ctx->SeqSJL < 1) {
                m++;
                continue;
            }
            if (ctx->PMRRN == 1)
                ctx->SeqSIL = ctx->SeqSJL = imin(ctx, ctx->SeqSIL,ctx->SeqSJL);

            seq_pps(ctx, 1,1,1,i,i);         /* write to test output */
            tmp = (double)(seq_mdist(ctx));          /* calculate distance */
            seq_pmd(ctx);                  /* print to stderr */

            if (dmax < tmp) {
                dmax = tmp;
                id = i;
            }
            fprintf(ctx->PMFd,"%5d %3d %3d ",i + 1,ctx->SeqSIL,ctx->SeqSJL);

            if (ctx->PMS == 1)           /* sequential */
                seq_pseq(ctx);
            else {   
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMTFmtS,(double)tmp);
                if (ctx->PMS == 2)      /* LCS */
                    seq_plcs(ctx);
            }           
            seq_pvar(ctx, i,-1);            
            fprintf(ctx->PMFd,"\n");
            nn++;
        }
    }
    printf1(ctx, "Number of sequences (cases): %d\n",ctx->NOC);
    printf1(ctx, "Sequences with zero length or internal gaps: %d\n",m);
    printf1(ctx, "Sequences used for alignment: %d\n",ctx->NOC - m);
#ifdef TDA_R_PACKAGE
    ctx->RExportRowKey[0] = '\0';
    tda_export_flush(ctx, "seqm.pairs");
#endif
    printf1(ctx, "\nNumber of alignments: %d\n",ctx->SeqNMD);
    printf1(ctx, "%d record(s) written to output file: %s\n",nn,ctx->PMFdName);
    if (nn > 1) {
        if (ctx->PMM == 1 && ctx->PMNCN <= 0)  
            printf1(ctx, "Maximum distance between sequences %d and %d: %g\n",id+1,jd+1,dmax);
        else if (ctx->PMM == 2)  
            printf1(ctx, "Maximum distance in case %d: %g\n",id,dmax);
    }
    return(nn);
}

/*--------------------------------------------------------------------------*/
/*  seq_mdist()             Calculate distance between sequences in         */
/*                          SeqSI and SeqSJ.                                */
/*                          Return: distance.                               */

float seq_mdist(TDAContext *ctx)
{
    register int i,j,k,l,isp;
    int is,js,j1,j2,j2p;
    double itmp,dtmp,stmp = 0.0;
    double ilen = 0.0,jlen = 0.0,u = 0.0,tmp,tmp1,tmp2;


    if (ctx->PMMax >= 1) {               /* restricted alignment */
        ilen = (double)((float)(ctx->SeqSIL + 1));
        jlen = (double)((float)(ctx->SeqSJL + 1));
        u = (ilen / jlen + jlen / ilen) / sqrt(ilen * ilen + jlen * jlen);
        u = (double)((float)(ctx->PMMax)) * (double)(u) / 1.4142;
    }
    j1 = 0;
    j2 = j2p = ctx->SeqSJL;

    for (i = 0; i < ctx->SeqSIL; ++i) {

        if (ctx->PMMax >= 1) {
            tmp = (double)(i + 1) / ilen;    
            tmp1 = jlen * (tmp - u);
            tmp2 = jlen * (tmp + u);
            if (tmp1 <= 1.0)
                j1 = 0;
            else
                j1 = (int)tmp1;
            if (tmp2 >= jlen)
                j2 = ctx->SeqSJL;
            else
                j2 = (int)tmp2;
        }
        is = ctx->SeqSI[i];                         
        isp = is * ctx->SeqNS;
        k = i * ctx->SeqMDim + j1;
        l = k + ctx->SeqMDim;

        for (j = j1; j < j2; ++j) {

            if (ctx->SeqMTyp == 1) {     /* indel functions */

                stmp = (double)(ctx->SeqMD[l] + (float)ctx->ICAlpha);
                if ((double)(stmp) > (double)((ctx->SeqME[l])) + (double)((float)ctx->ICBeta))
                    stmp = (double)(ctx->SeqME[l] + (float)ctx->ICBeta);
                l++;
                ctx->SeqME[l] = (float)(stmp);

                stmp = (double)(ctx->SeqMD[k + 1] + (float)ctx->ICAlpha);
                if ((double)(stmp) > (double)((ctx->SeqMF[k + 1])) + (double)((float)ctx->ICBeta))
                    stmp = (double)(ctx->SeqMF[k + 1] + (float)ctx->ICBeta);
                ctx->SeqMF[l] = (float)(stmp);
            }
            js = ctx->SeqSJ[j];

            if (is == js)
                stmp = 0.0;

            else if (ctx->SCostTyp == 0)     /* time-independent substitution cost */
                stmp = (double)(ctx->SCost[isp + js]);
 
            else if (ctx->SCostTyp == 1) {   /* absolute difference */

                stmp = fabs((double)(ctx->SeqSTNI[ctx->SeqSN][is] - ctx->SeqSTNI[ctx->SeqSN][js]));
    
            }
            stmp += (double)(ctx->SeqMD[k]);
            k++; 

            if (ctx->SeqMTyp == 0) {
                if (j > j1 || j == 0) {
                    itmp = (double)(ctx->ICost[j] + ctx->SeqMD[l]);
                    if (stmp > itmp)
                        stmp = itmp;
                }
                l++;
                if (j < j2p) {
                    dtmp = (double)(ctx->ICost[i] + ctx->SeqMD[k]);
                    if (stmp > dtmp)
                        stmp = dtmp;
                }
            }
            else if (ctx->SeqMTyp == 1) {
                if ((double)(stmp) > (double)(ctx->SeqME[l]))
                    stmp = (double)(ctx->SeqME[l]);
                if ((double)(stmp) > (double)(ctx->SeqMF[l]))
                    stmp = (double)(ctx->SeqMF[l]);
            }
            ctx->SeqMD[l] = (float)(stmp);
        }
        j2p = j2;
    }
    if (ctx->PMF1Def && ctx->PMTST[3])
        seq_dmat(ctx);
    return ((float)(stmp));
}

/*--------------------------------------------------------------------------*/
/*  seqm_pcost()    Print cost definitions to test output file.             */

void seqm_pcost(TDAContext *ctx)
{
    if (ctx->PMF1Def == 0)
        return;

    fprintf(ctx->PMF1d,"Optimal matching test output file.\n");
    fprintf(ctx->PMF1d,"Number of states: %d\n",ctx->SeqNS);
    fprintf(ctx->PMF1d,"Max sequence lenght: %d\n",ctx->SeqMLen);

    if (ctx->ICostTyp == 1)
        fprintf(ctx->PMF1d,"\nLinear indel cost function: alpha=%g, beta=%g.\n",(double)(ctx->ICAlpha),(double)(ctx->ICBeta));
    else {
        fprintf(ctx->PMF1d,"\nIndel cost\n");
        mprf(ctx, 1,ctx->SeqMLen,ctx->ICost - 1,ctx->PMFmtS,ctx->PMF1d);
    }
    if (ctx->SCostTyp == 0) {
        fprintf(ctx->PMF1d,"\nSubstitution cost\n");
        mprf(ctx, ctx->SeqNS,ctx->SeqNS,ctx->SCost - 1,ctx->PMFmtS,ctx->PMF1d);
    }
    else if (ctx->SCostTyp == 1)
        fprintf(ctx->PMF1d,"\nSubstitution cost: absolute difference.\n");
    fprintf(ctx->PMF1d,"\n");
}

/*--------------------------------------------------------------------------*/
/*  seqm_dtda(n)     Write TDA description file, n = number of cases.       */

void seqm_dtda(TDAContext *ctx, int n)
{
    register int i,j,k,l;

    fprintf(ctx->PMTDAFd,"nvar(\n");
    fprintf(ctx->PMTDAFd,"  dfile = %s,\n",ctx->PMFdName);
    fprintf(ctx->PMTDAFd,"  noc = %d,\n",n);

    k = 1;

    if (ctx->PMM == 1) {

        if (ctx->PMNCN <= 0) {

            fprintf(ctx->PMTDAFd,"  SEQ1N <5>[6.0] = c%-3d, # sequence A case number\n",k++);
            fprintf(ctx->PMTDAFd,"  SEQ2N <5>[6.0] = c%-3d, # sequence B case number\n",k++);
            fprintf(ctx->PMTDAFd,"  SEQ1L <2>[6.0] = c%-3d, # sequence A length\n",k++);
            fprintf(ctx->PMTDAFd,"  SEQ2L <2>[6.0] = c%-3d, # sequence B length\n",k++);

            if (ctx->PMS == 1) {
                for (i = 1; i <= ctx->SeqMLen; ++i) 
                    fprintf(ctx->PMTDAFd,"  DIST%d <4>[%d.%d] = c%-3d, # distance, time %d\n",i,ctx->PMTFmt1,ctx->PMTFmt2,k++,i);
            }
            else {
                fprintf(ctx->PMTDAFd,"  DIST  <4>[%d.%d] = c%-3d, # distance\n",ctx->PMTFmt1,ctx->PMTFmt2,k++);
                
                if (ctx->PMS == 2) {
                    fprintf(ctx->PMTDAFd,"  LCSL  <2>[4.0] = c%-3d, # length of LCS\n",k++);
                    for (i = 1; i <= ctx->SeqMLen; ++i)
                        fprintf(ctx->PMTDAFd,"  LCS%d <2>[2.0] = c%-3d, # LCS t=%d\n",i,k++,i);
                }
            }                 
        }            
        else {
            fprintf(ctx->PMTDAFd,"  SEQ1N <5>[6.0] = c%-3d, # sequence A case number\n",k++);
            fprintf(ctx->PMTDAFd,"  SEQ1L <2>[6.0] = c%-3d, # sequence A length\n",k++);

            for (i = 0; i < ctx->PMNCN; ++i) {
                j = ctx->PMCN[i];      
                if (j > 0) {
                    if (ctx->PMS == 1) {
                        for (l = 1; l <= ctx->SeqMLen; ++l) 
                            fprintf(ctx->PMTDAFd,"  DIST%d_%d <4>[%d.%d] = c%-3d, # distance, t=%d, cn=%d\n",j,l,ctx->PMTFmt1,ctx->PMTFmt2,k++,l,j);
                    }
                    else {
                        fprintf(ctx->PMTDAFd,"  DIST%d <4>[%d.%d] = c%-3d, # distance to sequence cn=%d\n",j,ctx->PMTFmt1,ctx->PMTFmt2,k++,j);
                        if (ctx->PMS == 2) {
                            fprintf(ctx->PMTDAFd,"  LCSL%d <2>[4.0] = c%-3d, # length of LCS, cn=%d\n",j,k++,j);
                            for (l = 1; l <= ctx->SeqMLen; ++l)
                                fprintf(ctx->PMTDAFd,"  LCS%d_%d <2>[2.0] = c%-3d, # LCS t=%d, cn=%d\n",j,l,k++,l,j);
                        }
                    }                  
                }
            }
        }
        for (i = 0; i < ctx->PMNV; ++i) {
            j = ctx->PMVIdx[i];
            fprintf(ctx->PMTDAFd,"  %s_A <%d>[%d.%d] = c%-3d, # sequence A\n",       
                                    ctx->VName[j],ctx->VSLen[j],ctx->VPFmt1[j],ctx->VPFmt2[j],k++);
            if (ctx->PMNCN <= 0)
                fprintf(ctx->PMTDAFd,"  %s_B <%d>[%d.%d] = c%-3d, # sequence B\n",         
                                    ctx->VName[j],ctx->VSLen[j],ctx->VPFmt1[j],ctx->VPFmt2[j],k++);
        }
    }
    else if (ctx->PMM == 2) {

        fprintf(ctx->PMTDAFd,"  SEQN  <5>[6.0] = c%-3d, # case number\n",k++);
        fprintf(ctx->PMTDAFd,"  SEQ1L <2>[6.0] = c%-3d, # sequence A length\n",k++);
        fprintf(ctx->PMTDAFd,"  SEQ2L <2>[6.0] = c%-3d, # sequence B length\n",k++);

        if (ctx->PMS == 1) {
            for (l = 1; l <= ctx->SeqMLen; ++l) 
                fprintf(ctx->PMTDAFd,"  DIST%d <4>[%d.%d] = c%-3d, # distance, t=%d\n",l,ctx->PMTFmt1,ctx->PMTFmt2,k++,l);
        }
        else {         
            fprintf(ctx->PMTDAFd,"  DIST  <2>[6.0] = c%-3d, # distance\n",k++);
            if (ctx->PMS == 2) {
                fprintf(ctx->PMTDAFd,"  LCSL  <2>[4.0] = c%-3d, # length of LCS\n",k++);
                for (i = 1; i <= ctx->SeqMLen; ++i)
                    fprintf(ctx->PMTDAFd,"  LCS%d <2>[2.0] = c%-3d, # LCS t=%d\n",i,k++,i);
            }
        }
        for (i = 0; i < ctx->PMNV; ++i) {
            j = ctx->PMVIdx[i];
            fprintf(ctx->PMTDAFd,"  %s <%d>[%d.%d] = c%-3d,\n",       
                                    ctx->VName[j],ctx->VSLen[j],ctx->VPFmt1[j],ctx->VPFmt2[j],k++);
        }
    }
    fprintf(ctx->PMTDAFd,");\n");
    printf1(ctx, "TDA description written to: %s\n",ctx->PMTDAFName);
}

/*--------------------------------------------------------------------------*/
/*  seq_scost()     Calculate substitution costs based on transition        */  
/*                  frequencies.                                            */
/*                                                                          */
/*  Return 0 if successful, otherwise -1.                                   */

int seq_scost(TDAContext *ctx)
{
    register int i,j,is,js;
    int n,m;
    double tmp;

    m = ctx->SeqNS * ctx->SeqNS;
    for (i = 0; i < m; ++i)
        ctx->SCost[i] = 0.0;
           
    m = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        n = seq_gst(ctx, i,1);
        if (n < 1) {
            m++;
            continue;
        }
        is = ctx->SeqSI[0];
        if (ctx->SCostTyp == 2) {
            for (j = 1; j < n; ++j) {
                js = ctx->SeqSI[j];
                ctx->SCost[is * ctx->SeqNS + js] =
                (float)((double)ctx->SCost[is * ctx->SeqNS + js] + 1.0);
                is = js;
            }
        }
    }
    if (ctx->SCostTyp == 2) {
        for (is = 0; is < ctx->SeqNS; ++is) {
            tmp = 0.0;
            j = is * ctx->SeqNS;
            for (js = 0; js < ctx->SeqNS; ++js)  
                tmp += (double)(ctx->SCost[j++]);

            j = is * ctx->SeqNS;
            for (js = 0; js < ctx->SeqNS; ++js) {
                /* Rohwer: SCost[j++] /= tmp, inside the if -- j advances
                   only when the row sum is positive.  The compound
                   assignment evaluates the subscript once; splitting it
                   for the explicit casts must not put j++ in the
                   subscript, which would read and modify j with no
                   sequence point between them. */
                if (tmp > 0.0) {
                    ctx->SCost[j] = (float)((double)ctx->SCost[j] /
                                            (double)tmp);
                    j++;
                }
            }
        }
        ctx->SCost[0] = 0.0;
        for (is = 1; is < ctx->SeqNS; ++is) {
            j = is * ctx->SeqNS;
            for (js = 0; js < is; ++js) {
                tmp = 2.0 - (double)((ctx->SCost[is * ctx->SeqNS + js])) - (double)(ctx->SCost[js * ctx->SeqNS + is]);
                ctx->SCost[j++] = (float)(tmp);
                ctx->SCost[js * ctx->SeqNS + is] = (float)(tmp);
            }
            ctx->SCost[j] = 0.0;
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

int seqpm(TDAContext *ctx)
{
    int err;

    err = -1;         
    if (check_cmd(ctx, 0))
        return(-1);
         
    printf1(ctx, "Sequence pattern matching. Current memory: %d bytes.\n",ctx->MemReq);
         
    if (parm(ctx, ctx->CmdBuf + 5,1,1))    /* get parameters */
        goto SEQPMFin;

    ctx->SeqSN = seq_getsn(ctx, ctx->PMSN,1);   /* get sequence number */
    if (ctx->SeqSN < 0)  
        goto SEQPMFin;

    ctx->SeqNS = ctx->SeqSTN[ctx->SeqSN];          /* number of states */
    ctx->SeqTA = ctx->SeqTMin[ctx->SeqSN];         /* range of time axis */
    ctx->SeqTB = ctx->SeqTMax[ctx->SeqSN];
    ctx->SeqMLen = ctx->SeqTB - ctx->SeqTA + 1;    /* max sequence length */
    ctx->SeqMDim = ctx->SeqMLen + 1;          /* dimension of D matrix */

    printf1(ctx, "Number of states: %d. Max sequence length: %d\n",ctx->SeqNS,ctx->SeqMLen);

    err = seq_search(ctx);
   
SEQPMFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seq_search()    Search for patterns.                                    */
/*  #                                                                      */
/*  Return 0 if successful, otherwise -1.                                   */

int seq_search(TDAContext *ctx)
{
    register int i,j,k;
    int err,m,nn,nm;

    err = -1;
    printf1(ctx, "\nSearching for patterns.\n");

    if (ctx->PMNPS <= 0) {
        p_err(ctx, -19,1);
        goto SSRCHFin;
    }

    ctx->SeqSPL = 0;                 /* max length of test pattern strings */
    for (i = 0; i < ctx->PMNPS; ++i) {
        printf1(ctx, "Pattern%2d: ",i + 1);
        seq_pprn(ctx, i);
        if (ctx->SeqSPL < ctx->PMPSN[i])
            ctx->SeqSPL = ctx->PMPSN[i];
    }
    ctx->SeqSPL++;   
    if (seqm_alloc1(ctx, 1)) {    /* allocate memory */
        p_err(ctx, -2,1);
        goto SSRCHFin;
    }
    if (ctx->PMF1Def) {      /* print to test output file */
        ctx->PMTST[2] = 1;
        fprintf(ctx->PMF1d,"Pattern search test output\n\n");
    }
    ctx->SeqNMD = 0;
    nn = m = 0;
#ifdef TDA_R_PACKAGE
    snprintf(ctx->RExportRowKey, sizeof(ctx->RExportRowKey), "seqpm.table");
#endif
    for (i = 0; i < ctx->NOC; ++i) {
        ctx->SeqSIL = seq_gst(ctx, i,1);
        if (ctx->SeqSIL < 1) {
            m++;
            continue;
        }
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->SeqSIL);
#ifdef TDA_R_PACKAGE
        tda_export_cell(ctx, "seqpm.table", (double)(i + 1));
        tda_export_cell(ctx, "seqpm.table", (double)ctx->SeqSIL);
#endif

        /* change to original state numbers */

        for (j = 0; j < ctx->SeqSIL; ++j)
            ctx->SeqSJ[j] = ctx->SeqSTNI[ctx->SeqSN][ctx->SeqSI[j]];
        ctx->SeqSJL = ctx->SeqSIL;

        seq_pps(ctx, 1,0,1,i,0);               /* write to test output */

        for (k = 0; k < ctx->PMNPS; ++k) {
            nm = seq_sfind(ctx, k);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nm);
#ifdef TDA_R_PACKAGE
            tda_export_cell(ctx, "seqpm.table", (double)nm);
#endif
   
            if (ctx->PMF1Def)            /* write to test output */
                seq_ppt(ctx, k);        
        }
        seq_pmd(ctx);                  /* print to stderr */
        seq_pvar(ctx, i,-1);             /* add variables */
        fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
        tda_export_endrow(ctx, "seqpm.table");
#endif
        nn++;
    }
#ifdef TDA_R_PACKAGE
    ctx->RExportRowKey[0] = '\0';
    tda_export_flush(ctx, "seqpm.table");
#endif
    printf1(ctx, "\nNumber of sequences (cases): %d\n",ctx->NOC);
    printf1(ctx, "Sequences with zero length or internal gaps: %d\n",m);
    printf1(ctx, "%d records written to output file: %s\n",nn,ctx->PMFdName);

    if (ctx->PMTDAFDef && nn > 0)         /* write TDA description */
        seqpm_dtda(ctx, nn);

    if (ctx->PMF1Def)  
        printf1(ctx, "Test output written to: %s\n",ctx->PMF1dName);

    err = 0;

SSRCHFin:
    seqm_alloc1(ctx, 0);     /* free memory */
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seq_sfind(k)    Find test pattern k.                                    */
/*                                                                          */
/*  Return number of matches.                                               */

int seq_sfind(TDAContext *ctx, int k)
{
    register int i,j,l,ll,j1,jj;
    int m,a,b,bl,n,nn,nnf,nf,u,fndflg;
    char c;

    m = ctx->PMPSN[k];       /* length of test pattern */

    if (ctx->PMF1Def) {      /* write test pattern to test output */

        fprintf(ctx->PMF1d,"Pattern:");   
        for (i = 0; i < m; ++i) {
            n = (int)ctx->PMPS[k][i];
            if (n >= 0) 
                fprintf(ctx->PMF1d," %d",n);
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
                fprintf(ctx->PMF1d," %c",c);
            }
        }
        fprintf(ctx->PMF1d,"\n");   
    }
    nf = 0;
    n = ctx->SeqSPL * ctx->SeqMDim;
    for (i = 1; i < n; ++i)
        ctx->SeqDP[i] = 0;
    ctx->SeqDP[0] = 1;

    j = 0; 
    for (i = 0; i < m; ++i) {   /* init first column */

        if (ctx->SeqDP[j] && (int)ctx->PMPS[k][i] <= -3) {
            j += ctx->SeqMDim;
            ctx->SeqDP[j] = 1;
        }
        else
            break;
    }
    bl = j1 = 0;
    while (j1 < ctx->SeqSJL) {

        fndflg = nn = nnf = 0;
        for (j = j1; j < ctx->SeqSJL; ++j) {

            l = j + 1;
            b = ctx->SeqSJ[j];
            n = 0;

            for (i = 0; i < m; ++i) {

                ll = l + ctx->SeqMDim;
                a = (int)ctx->PMPS[k][i];
                u = 0;
                if (ctx->SeqDP[l] && a <= -3)
                    u = 1; 
                else {
                    jj = 1;
                    if (j == j1)
                        jj += j1;

                    if (ctx->SeqDP[l - jj]) {
                        if ((a >= 0 && a != b) || (a == -2 && b != bl))
                            ;
                        else
                            u = 1;
                    }
                    if (u == 0 && ctx->SeqDP[ll - jj]) {
                        if (a == -4 || (a <= -2 && b == bl))
                            u = 1;
                    }
                }
                if (u) {
                    ctx->SeqDP[ll] = 1;
                    n++;
                }
                l = ll;
       
            }
            if (n == 0)  
                break;
                   
            if (ctx->SeqDP[m * ctx->SeqMDim + j + 1]) {
                fndflg = 1;
                nnf = nn;
            }
            nn++;
            bl = b;
        }
        if (fndflg) {
            nf++;
            if (ctx->PMF1Def) {      /* write match to test output */

                fprintf(ctx->PMF1d,"  Match:");
                for (l = j1; l <= j1 + nnf; ++l)
                    fprintf(ctx->PMF1d," %d",ctx->SeqSJ[l]);
                fprintf(ctx->PMF1d,"\n");
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

void seq_ppt(TDAContext *ctx, int k)
{
    register int i,j;
    int n;
    char c;

    if (ctx->PMF1Def== 0)
        return;

    fprintf(ctx->PMF1d,"\nD Matrix\n");
    fprnchar(ctx, ctx->PMF1d,' ',9,0);
    for (j = 0; j < ctx->SeqSIL; ++j) {
        n = ctx->SeqSTNI[ctx->SeqSN][ctx->SeqSI[j]];
        fprintf(ctx->PMF1d,"%3d ",n);
    }
    fprintf(ctx->PMF1d,"\n    ");
    fprnchar(ctx, ctx->PMF1d,'-',4 + 4 * ctx->SeqSIL,1);

    for (i = 0; i < ctx->SeqSPL; ++i) {
        if (i == 0)
            fprintf(ctx->PMF1d,"    |");
        else if (i - 1 >= ctx->PMPSN[k]) {
            /*  SeqSPL is the LONGEST pattern of the set, but PMPS[k] is
                only PMPSN[k] long, so this loop read past the end of
                every pattern shorter than the longest (seqpm1.cf).  Those positions have no
                symbol, so they print blank.  */
            fprintf(ctx->PMF1d,"    |");
        }
        else {
            n = (int)ctx->PMPS[k][i - 1];
            if (n >= 0) 
                fprintf(ctx->PMF1d,"%3d |",n);
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
                fprintf(ctx->PMF1d,"  %c |",c);
            }
        }
        for (j = 0; j <= ctx->SeqSIL; ++j)  
            fprintf(ctx->PMF1d,"%3d ",(int)ctx->SeqDP[i * ctx->SeqMDim + j]);
        fprintf(ctx->PMF1d,"\n");
    }
    fprintf(ctx->PMF1d,"\n");
}   

/*--------------------------------------------------------------------------*/
/*  seq_pprn(i)     Print test pattern i. Return 1 if pattern contains      */
/*                  special characters, otherwise 0.                        */

int seq_pprn(TDAContext *ctx, int i)
{
    register int j,k,n;
    char c;

    n = 0;
    for (j = 0; j < ctx->PMPSN[i]; ++j) {
        if (j)  
            printf1(ctx, ",");
        k = (int)ctx->PMPS[i][j];
        if (k >= 0) 
            printf1(ctx, "%d",k);
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
            printf1(ctx, "%c",c);
            n = 1;
        }
    }
    newline(ctx);
    return(n);
}

/*--------------------------------------------------------------------------*/
/*  seqm_alloc1(opt)    If opt != 0 allocate, else free memory.             */
/*                      Return 0 if OK, -1 if error.                        */

int seqm_alloc1(TDAContext *ctx, int opt)
{
    int n,err = 0;

    if (opt) {
        err = -1;
        n = ctx->SeqSPL * ctx->SeqMDim;
        if (!(ctx->SeqDP = (char *)calloc((size_t)(n),sizeof(char))))   
            goto SEQMA1Fin;
        memrq(ctx, n,sizeof(char));
        ctx->SeqDPA = n;
    
        n = ctx->SeqMLen;
        if (!(ctx->SeqSI = (int *)calloc((size_t)(n),sizeof(int))))   
            goto SEQMA1Fin;
        memrq(ctx, n,sizeof(int));
        ctx->SeqSIA = n;

        if (!(ctx->SeqSJ = (int *)calloc((size_t)(n),sizeof(int))))   
            goto SEQMA1Fin;
        memrq(ctx, n,sizeof(int));
        ctx->SeqSJA = n;

        return(0);
    }

SEQMA1Fin:
    if (ctx->SeqDPA) {
        free((char *)ctx->SeqDP);
        memrq(ctx, -ctx->SeqDPA,sizeof(char));
        ctx->SeqDPA = 0;
    }
    if (ctx->SeqSIA) {
        free((char *)ctx->SeqSI);
        memrq(ctx, -ctx->SeqSIA,sizeof(int));
        ctx->SeqSIA = 0;
    }
    if (ctx->SeqSJA) {
        free((char *)ctx->SeqSJ);
        memrq(ctx, -ctx->SeqSJA,sizeof(int));
        ctx->SeqSJA = 0;
    }
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqpm_dtda(n)   Write TDA description file, n = number of cases.        */

void seqpm_dtda(TDAContext *ctx, int n)
{
    register int i,j,k;

    fprintf(ctx->PMTDAFd,"nvar(\n");
    fprintf(ctx->PMTDAFd,"  dfile = %s,\n",ctx->PMFdName);
    fprintf(ctx->PMTDAFd,"  noc = %d,\n",n);

    fprintf(ctx->PMTDAFd,"  SID  <5>[%d.0] = c1  , # sequence case number\n",ctx->PMNFmt);
    fprintf(ctx->PMTDAFd,"  SLEN <5>[%d.0] = c2  , # sequence length\n",ctx->PMNFmt);
    k = 3;
    for (i = 0; i < ctx->PMNPS; ++i) {
        fprintf(ctx->PMTDAFd,"  SM%-2d <5>[%d.0] = c%-3d, # number of matches, pattern %d\n",k-2,ctx->PMNFmt,k,k-2);
        k++;
    }
    for (i = 0; i < ctx->PMNV; ++i) {
        j = ctx->PMVIdx[i];
        fprintf(ctx->PMTDAFd,"  %s <%d>[%d.%d] = c%-3d,\n",       
                                    ctx->VName[j],ctx->VSLen[j],ctx->VPFmt1[j],ctx->VPFmt2[j],k++);
    }
    fprintf(ctx->PMTDAFd,");\n");
    printf1(ctx, "TDA description written to: %s\n",ctx->PMTDAFName);
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

int subm(TDAContext *ctx)
{
    register int i,j,k,l,ix,iy;
    int err,m,n,m1,n1,ni,nj,r,sm,smr,smc;
    double x,delta;

    err = -1;         
    if (check_cmd(ctx, 0))
        return(-1);
         
    printf1(ctx, "Substitution metric. Current memory: %d bytes.\n",ctx->MemReq);
         
    if (parm(ctx, ctx->CmdBuf + 4,4,1))    /* get parameters */
        goto SUBMFin;

    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    if (ctx->PMNV < 2) {
        printf1(ctx, "Error: need at least two variables.\n");
        goto SUBMFin;
    }
    if (ctx->PMSCOSTMAT >= 0) {
        sm = ctx->PMSCOSTMAT;
        printf1(ctx, "Substitution cost defined by matrix: %s\n",ctx->MatName[sm]);
        smr = ctx->MatRow[sm];
        smc = ctx->MatCol[sm];
        if (smr < ctx->NOC || smc < ctx->NOC) {
            printf1(ctx, "Error: need at least an (%d,%d) matrix.\n",ctx->NOC,ctx->NOC);
            goto SUBMFin;
        }
        for (i = 0; i < ctx->NOC; ++i) {
            for (j = 1; j <= ctx->NOC; ++j) {
                if (ctx->MatVal[sm][i * smc + j] < 0.0) {
                    printf1(ctx, "Error: substitution cost matrix contains negative values.\n");
                    goto SUBMFin;
                }
            }
        }
    }

    /* check distributions */

    n = 0;
    for (j = 0; j < ctx->PMNV; ++j) {
        ix = ctx->PMVIdx[j];
        n = 0;
        for (i = 0; i < ctx->NOC; ++i) {
            x = get_data(ctx, ix,i);
            if (x < 0.0 || x > 1.0) {
                n = 1;
                break;
            }
        }
/*****  if (n || fabs(xsum - 1.0) > PMEPS) { *****/
        if (n) {
            tda_out("Error: no distribution in variable %s\n",ctx->VName[ix]);
            goto SUBMFin;
        }
    }   
    if (alloc_aci(ctx, ctx->NOC + 1))      
        goto SUBMFin;
    if (alloc_acj(ctx, ctx->NOC + 1))      
        goto SUBMFin;
    if (alloc_acu(ctx, ctx->NOC + 1))      
        goto SUBMFin;
    if (alloc_acv(ctx, ctx->NOC + 1))      
        goto SUBMFin;

    ix = ctx->PMVIdx[0];
    for (k = 1; k < ctx->PMNV; ++k) {
        delta = 0.0;

        ni = nj = 0;
        iy = ctx->PMVIdx[k];
        for (i = 0; i < ctx->NOC; ++i) {
            x = get_data(ctx, ix,i) - get_data(ctx, iy,i);
            if (x > ctx->PMEPS) {
                ni++;
                ctx->AcI[ni] = i+ 1;
                ctx->AcU[ni] = x;
            }
            else if (x < -ctx->PMEPS) {
                nj++;
                ctx->AcJ[nj] = i + 1;
                ctx->AcV[nj] = -x;
            }              
        }
        if (ni > 0 && nj > 0) {
            n = ni * nj;
            m = ni + nj;
            n1 = n + 1;
            m1 = m + 1;
            if (alloc_acz(ctx, n1 * m1 + 1))      
                goto SUBMFin;
            if (alloc_acx(ctx, n1 + 1))      
                goto SUBMFin;
            if (alloc_acy(ctx, m1 + 1))      
                goto SUBMFin;
    
            l = 1;
            for (i = 1; i <= ni; ++i) {
                for (j = 1; j <= nj; ++j)  
                    ctx->AcZ[l++] = -subm_cost(ctx, ctx->AcI[i],ctx->AcJ[j]);
            }
            l = 0;
            for (i = 1; i <= ni; ++i) {
                for (j = 1; j <= nj; ++j)  
                    ctx->AcZ[i * n1 + l + j] = 1.0;
                l += nj;   
            }
      
            for (i = 1; i <= nj; ++i) {
                for (j = 0; j < ni; ++j) {
                    ctx->AcZ[(ni + i) * n1 + i + j * nj] = 1.0;
                }
            }
            l = 2 * n1;
            for (i = 1; i <= ni; ++i) { 
                ctx->AcZ[l] = ctx->AcU[i]; 
                l += n1;
            }
            for (j = 1; j <= nj; ++j) {
                ctx->AcZ[l] = ctx->AcV[j];
                l += n1;
            }
            /********** 
            for (i = 1; i <= m1; ++i) {
                for (j = 1; j <= n1; ++j)  
                    tda_out("%lf ",AcZ[(i - 1) * n1 + j]);
                newline(ctx);
            }
            *****/
            if (ctx->PMF1Def) {
                for (i = 1; i <= ctx->NOC; ++i) {
                    for (j = 1; j <= ctx->NOC; ++j)  
                        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,subm_cost(ctx, i,j));
                    fprintf(ctx->PMF1d,"\n");
                }
                fprintf(ctx->PMF1d,"\n");
 
                for (i = 0; i < ctx->NOC; ++i) {
                    fprintf(ctx->PMF1d,"%4d ",i + 1);
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,get_data(ctx, ix,i));
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,get_data(ctx, iy,i));
                    fprintf(ctx->PMF1d,"\n");
                }                                      
                fprintf(ctx->PMF1d,"\n");

                for (i = 1; i <= ni; ++i) {
                    fprintf(ctx->PMF1d,"%4d %4d ",i,ctx->AcI[i]);
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcU[i]);
                    fprintf(ctx->PMF1d,"\n");
                }
                fprintf(ctx->PMF1d,"\n");
                for (j = 1; j <= nj; ++j) {
                    fprintf(ctx->PMF1d,"%4d %4d ",j,ctx->AcJ[j]);
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcV[j]);
                    fprintf(ctx->PMF1d,"\n");
                }
                fprintf(ctx->PMF1d,"\n");
                for (i = 1; i <= m1; ++i) {
                    for (j = 1; j <= n1; ++j)  
                        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcZ[(i - 1) * n1 + j]);
                    fprintf(ctx->PMF1d,"\n");
                }
                fprintf(ctx->PMF1d,"\n");
            }

            r = lpf1(ctx, m-1,n,m-1,ctx->AcZ,ctx->AcX,ctx->AcY,&delta,ctx->PMEPS);
            if (r) {
                switch (r) {
                    case -2:    printf1(ctx, "No solution.\n");
                                break;
                    case -3:    printf1(ctx, "No solution. Dual objective function is unbounded.\n");
                                break;
                    case -4:    printf1(ctx, "No solution. Primal objective function is unbounded.\n");
                                break;
                    case -5:    printf1(ctx, "Error: insufficient memory.\n");
                                break;
                    default:    printf1(ctx, "Error (%d).\n",r);
                                break;
                }   
                err = 0;
                goto SUBMFin;
            }
            printf1(ctx, "Distance: ");
            rt_printf1_d(ctx, ctx->PMFmtS,-delta);
            newline(ctx);


            if (ctx->PMF1Def) {
                l = 1;
                for (i = 1; i <= ni; ++i) {
                    for (j = 1; j <= nj; ++j) {
                        fprintf(ctx->PMF1d,"%4d %4d ",ctx->AcI[i],ctx->AcJ[j]);
                        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[l++]);
                        fprintf(ctx->PMF1d,"\n");
                    }
                }
                printf1(ctx, "\nSubstitution vector written to: %s\n",ctx->PMF1dName);
            }
        }
    }
    err = 0;
   
SUBMFin:
    p_clean(ctx);
    return(err);
}

double subm_cost(TDAContext *ctx, int i,int j)    
{
    int sm,smc;
    double c;

    if (ctx->PMSCOSTMAT >= 0) { 
        sm = ctx->PMSCOSTMAT;
        smc = ctx->MatCol[sm];
        c = ctx->MatVal[sm][(i - 1) * smc + j];
    }
    else
        c = fabs((double)(i - j));
    return(c);
}


