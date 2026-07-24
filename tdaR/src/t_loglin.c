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
#include "tda_context.h"


/*  functions in t_loglin.c */

int loglin(TDAContext *ctx);
int ll_mcheck(TDAContext *ctx);
int ll_dcol(TDAContext *ctx);
int ll_alloc(TDAContext *ctx, int opt);
void ll_free(TDAContext *ctx);
int ll_tab(TDAContext *ctx);
int loglin1(TDAContext *ctx);
void ll_ppar(TDAContext *ctx);
void ll_error(TDAContext *ctx, int err,char *mod);
void prn_mdes(TDAContext *ctx, FILE *fd, char *desc, int n);
void prn_lldmat(TDAContext *ctx, int i);
void prn_llcov(TDAContext *ctx, int i);
void prn_ctab(TDAContext *ctx, int mode);
int llmodcheck(TDAContext *ctx, char *desc);
int expand1(TDAContext *ctx, char *s, char *r, int mode);
int expand2(TDAContext *ctx, char *s, char *r, int mode);
char *distribute(TDAContext *ctx, int n, char *r, int mode, int *len);
int build_dm(TDAContext *ctx, char *descr, int col, int *slen, int mode);
int design_comb(TDAContext *ctx, int dn, int *cat, int *ptr);
int ll_readdm(TDAContext *ctx, int i, int dcol, int mode);
int logl_fit(TDAContext *ctx);
void logl_solve(TDAContext *ctx, int n, double *lu, double *b, double *x);
void logl_inv(TDAContext *ctx, int n, double *tri, double *xinv);
int logl_chol(TDAContext *ctx, int n, double *sym, double *tri);
/* void logl_gen(double *c, double *obs, double *expv, double *res); */
int logscrn(TDAContext *ctx);
int logscr_conf(TDAContext *ctx, int m, int *iset, int *jset);
int logscr_comb(TDAContext *ctx, int *iset, int n, int m, int last);
int logscr_eval(TDAContext *ctx, short *iar, int nc, int nv, int ibeg);
void logscr_res(TDAContext *ctx);
double logscr_lr(TDAContext *ctx);
void logscr_fit(TDAContext *ctx, int ncon);
int e_expm(TDAContext *ctx);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

#define MaxDX 200           /* max number of model terms                    */
#define MaxABC 200 











/*  The following variables are used for marginal and partial association   */



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

int loglin(TDAContext *ctx)
{
    int err;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Loglinear models. Current memory: %d bytes.\n\n",ctx->MemReq);

    ctx->MxItFlg = 0;
    ctx->TOLF = 1.e-8;

    if (parm(ctx, ctx->CmdBuf + 6,4,1))     /* get parameters */
        goto LLFin;

    if (ctx->PMNFmtF == 0) {
        ctx->PMNFmt = 2;
        makenfmt(ctx, &ctx->PMNFmt,ctx->PMNFmtS,sizeof(ctx->PMNFmtS),1,ctx->SEPC);
    }
    if (ctx->MxIter < 1 || ctx->MxItFlg == 0)
        ctx->MxIter = 20;

    ctx->NDIM = ctx->PMNV;                /* number of table dimensions */
    if (ctx->NDIM > 26) {
        printf1(ctx, "Error: max number of dimensions is 26.\n");
        goto LLFin;
    }
    if (ll_tab(ctx))               /* make table */
        goto LLFin;

    if (ctx->PMTabFDef) {            /* print table */ 
        prn_ctab(ctx, 1);
        printf1(ctx, "Table written to: %s\n",ctx->PMTabFName);
    }
    newline(ctx);
     
    if (ctx->SCRNFlg) {              /* screening */
        if (logscrn(ctx))
            goto LLFin;
    }

    /* check model requests */

    if (ctx->PMModN == 0) {
        err = 0;
        goto LLFin;
    }
    if (ll_mcheck(ctx))    /* check model descriptions */
        goto LLFin;

    if (ll_dcol(ctx))      /* get max number of design matrix columns */
        goto LLFin;

    if (ll_alloc(ctx, 1))    /* allocate memory for model estimation */
        goto LLFin;

    if (loglin1(ctx))      /* estimate models */
        goto LLFin;

    err = 0;

LLFin:
    con_free(ctx);
    ll_alloc(ctx, 0);
    ll_free(ctx);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ll_mcheck()   Check model descriptions                                  */
/*                Return 0 if OK, -1 if error.                              */

int ll_mcheck(TDAContext *ctx)
{
    FILE *fd;
    register int i;
    int err,n;
    char *p = NULL;

    err = -1;

    printf1(ctx, "Check of model requests.\n");

    for (i = 0; i < ctx->PMModN; ++i) {
        p = ctx->PMModS[i];
        if (*p == 'd' && *(p + 1) == ':') {
            if ((fd = fopen(p + 2,OPEN_RD))) {   
                ctx->LLTyp[i] = -1;
                ctx->LLTyp2++;
                fclose(fd);
            }
            else {
                ll_error(ctx, -9,p);
                goto LLMCFin;
            }
        }
        else {
            n = llmodcheck(ctx, p);
            if (n < 0) {
                ll_error(ctx, n,p);
                goto LLMCFin;
            }
            ctx->LLTyp[i] = n;
            ctx->LLTyp1++;
        }
    }
    ctx->LLD1Len = 1;

    for (i = 0; i < ctx->PMModN; ++i) {      /* get max length of first expansion */
        if (ctx->LLTyp[i] > 0) {
            n = expand1(ctx, ctx->PMModS[i],p,0);
            if (n < 0) {
                ll_error(ctx, n,ctx->PMModS[i]);
                goto LLMCFin;
            }
            ctx->LLTyp[i] = n;
            ctx->LLD1Len = imax(ctx, ctx->LLD1Len,n);
        }
    }   

    /* allocate memory for first expansion */

    if (!(ctx->LLD1 = (char *) calloc ((size_t)(ctx->LLD1Len + 2),sizeof(char)))) {
        p_err(ctx, -2,1);
        goto LLMCFin;
    }
    ctx->LLD1A = ctx->LLD1Len + 2;
    memrq(ctx, ctx->LLD1A,sizeof(char));

    /* allocate memory for second expansion */

    ctx->ExMax = ctx->LLD1Len;

    if (!(ctx->ExIdx  = (int *)  calloc((size_t)(ctx->ExMax + 1),sizeof(int)))) {
        p_err(ctx, -2,1);
        goto LLMCFin;
    }
    ctx->ExIdxA = ctx->ExMax + 1;
    memrq(ctx, ctx->ExIdxA,sizeof(int));

    if (!(ctx->ExVc   = (int *)  calloc((size_t)(ctx->ExMax + 1),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto LLMCFin;
    }
    ctx->ExVcA = ctx->ExMax + 1;
    memrq(ctx, ctx->ExVcA,sizeof(int));

    if (!(ctx->ExV = (char *) calloc((size_t)(ctx->ExMax + 1),sizeof(char)))) {
        p_err(ctx, -2,1);
        goto LLMCFin;
    }
    ctx->ExVA = ctx->ExMax + 1;
    memrq(ctx, ctx->ExVA,sizeof(char));

    if (!(ctx->ExVl   = (char *) calloc((size_t)(ctx->ExMax + 1),sizeof(char)))) {
        p_err(ctx, -2,1);
        goto LLMCFin;
    }
    ctx->ExVlA = ctx->ExMax + 1;
    memrq(ctx, ctx->ExVlA,sizeof(char));

    if (!(ctx->ExVTyp = (char *) calloc((size_t)(ctx->ExMax + 1),sizeof(char)))) {
        p_err(ctx, -2,1);
        goto LLMCFin;
    }
    ctx->ExVTypA = ctx->ExMax + 1;
    memrq(ctx, ctx->ExVTypA,sizeof(char));

    if (!(ctx->ExOpl  = (char *) calloc((size_t)(ctx->ExMax + 1),sizeof(char)))) {
        p_err(ctx, -2,1);
        goto LLMCFin;
    }
    ctx->ExOplA = ctx->ExMax + 1;
    memrq(ctx, ctx->ExOplA,sizeof(char));

    if (!(ctx->ExMult = (char *) calloc((size_t)(ctx->ExMax + 1),sizeof(char)))) {
        p_err(ctx, -2,1);
        goto LLMCFin;
    }
    ctx->ExMultA = ctx->ExMax + 1;
    memrq(ctx, ctx->ExMultA,sizeof(char));

    if (!(ctx->ExUsed = (char *) calloc((size_t)(ctx->ExMax + 1),sizeof(char)))) {
        p_err(ctx, -2,1);
        goto LLMCFin;
    }
    ctx->ExUsedA = ctx->ExMax + 1;
    memrq(ctx, ctx->ExUsedA,sizeof(char));

    ctx->LLD2Len = 1;                    /* get max length of second expansion */

    for (i = 0; i < ctx->PMModN; ++i) {
        if (ctx->LLTyp[i] > 0) {
            n = expand1(ctx, ctx->PMModS[i],ctx->LLD1,1);
            n = expand2(ctx, ctx->LLD1,p,0);
            if (n < 0) {
                ll_error(ctx, n,ctx->PMModS[i]);
                goto LLMCFin;
            }
            ctx->LLTyp[i] = n;
            ctx->LLD2Len = imax(ctx, ctx->LLD2Len,n);
        }
    }   

    /* allocate memory for second expansion */

    if (!(ctx->LLD2 = (char *) calloc ((size_t)(ctx->LLD2Len + 2),sizeof(char)))) {
        p_err(ctx, -2,1);
        goto LLMCFin;
    }
    ctx->LLD2A = ctx->LLD2Len + 2;
    memrq(ctx, ctx->LLD2A,sizeof(char));

    err = 0;

LLMCFin:
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ll_dcol     Get number of design matrix columns.                        */
/*              Return 0 if OK, -1 if error.                                */

int ll_dcol(TDAContext *ctx)
{
    register int i;
    int n,err,slen;

    err = -1;

    /*  allocate memory for number of design matrix columns */

    if (!(ctx->LLDCol = (int *)calloc((size_t)(ctx->PMModN),sizeof(int)))) {
        p_err(ctx, -2,1);
        goto LLDCFin;
    }
    ctx->LLDColA = ctx->PMModN;        
    memrq(ctx, ctx->LLDColA,sizeof(int));

    ctx->LLDColM = 1;

    for (i = 0; i < ctx->PMModN; ++i) {     
        if (ctx->LLTyp[i] > 0) {
            n = expand1(ctx, ctx->PMModS[i],ctx->LLD1,1);
            n = expand2(ctx, ctx->LLD1,ctx->LLD2,1);
            n = build_dm(ctx, ctx->LLD2,0,&slen,0);
            if (n < 0) {
                ll_error(ctx, n,ctx->PMModS[i]);
                goto LLDCFin;
            }
            ctx->LLDCol[i] = n;
            ctx->LLDColM = imax(ctx, ctx->LLDColM,n);
        }
        else if (ctx->LLTyp[i] < 0) {
            n = ll_readdm(ctx, i,0,0);
            if (n < 0) {
                ll_error(ctx, n,ctx->PMModS[i]);
                goto LLDCFin;
            }
            ctx->LLDCol[i] = n;
            ctx->LLDColM = imax(ctx, ctx->LLDColM,n);
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

int ll_alloc(TDAContext *ctx, int opt)
{
    int err = 0;

    if (opt == 0)
        goto LLAFin;

    err = -1;

    /*  allocate memory for model parameter labels */

    ctx->CTVLenM = (ctx->LLDColM + 1) * (ctx->NDIM * 20 + 1);

    if (!(ctx->CTVLab = (char *) calloc((size_t)(ctx->CTVLenM + 10),sizeof(char)))) {
        p_err(ctx, -2,1);
        goto LLAFin;
    }
    ctx->CTVLabA = ctx->CTVLenM + 10;
    memrq(ctx, ctx->CTVLabA,sizeof(char));

    /*  allocate memory for design matrix LLDMat */

    if (!(ctx->LLDMat = (short *)calloc((size_t)(ctx->NCTab) * (size_t)(ctx->LLDColM) + 1,sizeof(short)))) { 
        p_err(ctx, -2,1);
        goto LLAFin;
    }
    ctx->LLDMatA = ctx->NCTab * ctx->LLDColM + 1;
    memrq(ctx, ctx->LLDMatA,sizeof(short));

    /*  memory allocation for model estimation, space is allocated for
        the largest model. */

    if (!(ctx->CTRes  = (double *) calloc((size_t)(ctx->NCTab + 1),sizeof(double)))) {  
        p_err(ctx, -2,1);
        goto LLAFin;
    }
    ctx->CTResA = ctx->NCTab + 1;
    memrq(ctx, ctx->CTResA,sizeof(double));

    if (!(ctx->CTMat  = (double *) calloc((size_t)(ctx->LLDColM) * (size_t)(ctx->LLDColM) + 1,sizeof(double)))) {  
        p_err(ctx, -2,1);
        goto LLAFin;
    }
    ctx->CTMatA = ctx->LLDColM + 1;
    memrq(ctx, ctx->CTMatA,sizeof(double));

    if (!(ctx->CTB = (double *) calloc((size_t)(ctx->LLDColM + 1),sizeof(double)))) { 
        p_err(ctx, -2,1);
        goto LLAFin;
    }
    ctx->CTBA = ctx->LLDColM + 1;
    memrq(ctx, ctx->CTBA,sizeof(double));

    if (!(ctx->CTPar  = (double *) calloc((size_t)(ctx->LLDColM + 1),sizeof(double)))) {  
        p_err(ctx, -2,1);
        goto LLAFin;
    }
    ctx->CTParA = ctx->LLDColM + 1;
    memrq(ctx, ctx->CTParA,sizeof(double));

    if (!(ctx->CTWork = (double *) calloc((size_t)(ctx->LLDColM + 1),sizeof(double)))) {
        p_err(ctx, -2,1);
        goto LLAFin;
    }
    ctx->CTWorkA = ctx->LLDColM + 1;
    memrq(ctx, ctx->CTWorkA,sizeof(double));

    if (!(ctx->CTFit = (double *) calloc((size_t)(ctx->NCTab + 1),sizeof(double)))) {
        p_err(ctx, -2,1);
        goto LLAFin;
    }
    ctx->CTFitA = ctx->NCTab + 1;
    memrq(ctx, ctx->CTFitA,sizeof(double));

    return(0);

LLAFin:
    if (ctx->CTVLabA > 0) {
        free((char *)ctx->CTVLab);
        memrq(ctx, -ctx->CTVLabA,sizeof(char));
        ctx->CTVLabA = 0;    
    }
    if (ctx->LLDMatA > 0) {
        free((char *)ctx->LLDMat);
        memrq(ctx, -ctx->LLDMatA,sizeof(short));
        ctx->LLDMatA = 0;    
    }
    if (ctx->CTResA > 0) {
        free((char *)ctx->CTRes);
        memrq(ctx, -ctx->CTResA,sizeof(double));
        ctx->CTResA = 0;    
    }
    if (ctx->CTMatA > 0) {
        free((char *)ctx->CTMat);
        memrq(ctx, -ctx->CTMatA,sizeof(double));
        ctx->CTMatA = 0;    
    }
    if (ctx->CTBA > 0) {
        free((char *)ctx->CTB);
        memrq(ctx, -ctx->CTBA,sizeof(double));
        ctx->CTBA = 0;    
    }
    if (ctx->CTParA > 0) {
        free((char *)ctx->CTPar);
        memrq(ctx, -ctx->CTParA,sizeof(double));
        ctx->CTParA = 0;    
    }
    if (ctx->CTWorkA > 0) {
        free((char *)ctx->CTWork);
        memrq(ctx, -ctx->CTWorkA,sizeof(double));
        ctx->CTWorkA = 0;    
    }
    if (ctx->CTFitA > 0) {
        free((char *)ctx->CTFit);
        memrq(ctx, -ctx->CTFitA,sizeof(double));
        ctx->CTFitA = 0;    
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ll_free. Free allocated memory.                                         */

void ll_free(TDAContext *ctx)
{
    register int i;

    if (ctx->CTNCatA > 0) {
        free((char *)ctx->CTNCat);
        memrq(ctx, -ctx->CTNCatA,sizeof(int));
        ctx->CTNCatA = 0;    
    }
    if (ctx->CTCatA > 0) {
        for (i = 0; i < ctx->NDIM; ++i) {
            if (ctx->CTCatAI[i] > 0) {
                free((char *)ctx->CTCat[i]);
                memrq(ctx, -ctx->CTCatAI[i],sizeof(int));
            }
        }
        free((char *)ctx->CTCat);
        memrq(ctx, -ctx->CTCatA,sizeof(int *));
        ctx->CTCatA = 0;    
    }
    if (ctx->CTCatAIA > 0) {
        free((char *)ctx->CTCatAI);
        memrq(ctx, -ctx->CTCatAIA,sizeof(int));
        ctx->CTCatAIA = 0;    
    }
    if (ctx->CTValA > 0) {
        free((char *)ctx->CTVal);
        memrq(ctx, -ctx->CTValA,sizeof(short));
        ctx->CTValA = 0;    
    }
    if (ctx->CTFreqA > 0) {
        free((char *)ctx->CTFreq);
        memrq(ctx, -ctx->CTFreqA,sizeof(int));
        ctx->CTFreqA = 0;    
    }
    if (ctx->CTSCALA > 0) {
        free((char *)ctx->CTSCAL);
        memrq(ctx, -ctx->CTSCALA,sizeof(double));
        ctx->CTSCALA = 0;    
    }
    if (ctx->LLD1A > 0) {
        free((char *)ctx->LLD1);
        memrq(ctx, -ctx->LLD1A,sizeof(char));
        ctx->LLD1A = 0;    
    }
    if (ctx->LLD2A > 0) {
        free((char *)ctx->LLD2);
        memrq(ctx, -ctx->LLD2A,sizeof(char));
        ctx->LLD2A = 0;    
    }
    if (ctx->ExIdxA > 0) {
        free((char *)ctx->ExIdx);
        memrq(ctx, -ctx->ExIdxA,sizeof(int));
        ctx->ExIdxA = 0;    
    }
    if (ctx->ExVcA > 0) {
        free((char *)ctx->ExVc);
        memrq(ctx, -ctx->ExVcA,sizeof(int));
        ctx->ExVcA = 0;    
    }
    if (ctx->ExVA > 0) {
        free((char *)ctx->ExV);
        memrq(ctx, -ctx->ExVA,sizeof(char));
        ctx->ExVA = 0;    
    }
    if (ctx->ExVlA > 0) {
        free((char *)ctx->ExVl);
        memrq(ctx, -ctx->ExVlA,sizeof(char));
        ctx->ExVlA = 0;    
    }
    if (ctx->ExVTypA > 0) {
        free((char *)ctx->ExVTyp);
        memrq(ctx, -ctx->ExVTypA,sizeof(char));
        ctx->ExVTypA = 0;    
    }
    if (ctx->ExOplA > 0) {
        free((char *)ctx->ExOpl);
        memrq(ctx, -ctx->ExOplA,sizeof(char));
        ctx->ExOplA = 0;    
    }
    if (ctx->ExMultA > 0) {
        free((char *)ctx->ExMult);
        memrq(ctx, -ctx->ExMultA,sizeof(char));
        ctx->ExMultA = 0;    
    }
    if (ctx->ExUsedA > 0) {
        free((char *)ctx->ExUsed);
        memrq(ctx, -ctx->ExUsedA,sizeof(char));
        ctx->ExUsedA = 0;    
    }
    if (ctx->LLDColA > 0) {
        free((char *)ctx->LLDCol);
        memrq(ctx, -ctx->LLDColA,sizeof(int));
        ctx->LLDColA = 0;    
    }
}

/* ------------------------------------------------------------------------ */
/*  ll_tab.  Create table.                                                  */
/*           Return 0 if OK, -1 if error.                                   */

int ll_tab(TDAContext *ctx)
{
    register int i,j,k,l;
    int err,m,n,eerr;
    double tmp;

    err = -1;

    printf1(ctx, "Definition and structure of contingency table.\n");
    if (ctx->PMWVar >= 0)  
        printf1(ctx, "Frequencies defined by: %s\n",ctx->VName[ctx->PMWVar]);
    if (ctx->PMSCAL >= 0)  
        printf1(ctx, "Scale variable defined by: %s\n",ctx->VName[ctx->PMSCAL]);
    newline(ctx);

    /*  Calculation of categories of dimensions of the table:
        NCTab = number of cells of the complete table.
        MaxYC = maximum number of categories in a dimension.
        CTNCat[i] = number of categories in i.th dimension,
        CTCat[i][j] (j = 1,2,3,...,CTNCat[i]) is the value of the j.th
        categorie in the i.th dimension. */

    if (!(ctx->CTNCat = (int *) calloc((size_t)(ctx->NDIM),sizeof(int)))) {   
        p_err(ctx, -2,1);
        goto LLTFin;
    }
    ctx->CTNCatA = ctx->NDIM;
    memrq(ctx, ctx->CTNCatA,sizeof(int));

    if (!(ctx->CTCatAI = (int *) calloc((size_t)(ctx->NDIM),sizeof(int)))) {   
        p_err(ctx, -2,1);
        goto LLTFin;
    }
    ctx->CTCatAIA = ctx->NDIM;
    memrq(ctx, ctx->CTCatAIA,sizeof(int));

    if (!(ctx->CTCat = (int **) calloc((size_t)(ctx->NDIM),sizeof(int *)))) {   
        p_err(ctx, -2,1);
        goto LLTFin;
    }
    ctx->CTCatA = ctx->NDIM;
    memrq(ctx, ctx->CTCatA,sizeof(int *));

    ctx->NCTab = 1;

    if (alloc_acr(ctx, ctx->NOC + 1))                 /* data */
        goto LLTFin;
    if (alloc_acs(ctx, ctx->PMMaxCat + 1))            /* values of categories */
        goto LLTFin;
    if (alloc_aci(ctx, ctx->PMMaxCat + 1))            /* frequencies */
        goto LLTFin;
    if (alloc_acj(ctx, ctx->PMMaxCat + 1))            /* bf pointer */
        goto LLTFin;

    eerr = 0;

    printf1(ctx, "Dimension  Variable ");
    prnchar(ctx, ' ',ctx->VNameLen - 8,0);
    printf1(ctx, "  Categories\n");

    for (i = 0; i < ctx->NDIM; ++i) {

        for (j = 0; j < ctx->NOC; ++j) {
            ctx->AcR[j + 1] = (int)get_data(ctx, ctx->PMVIdx[i],j);
        }
        n = cfreq(ctx, 1,ctx->NOC,ctx->AcR,ctx->PMMaxCat,ctx->AcS,ctx->AcI,ctx->AcJ,0,&tmp,&tmp);

        if (n < 1) {
            printf1(ctx, "Error: can't sort dimension %d.\n",i + 1);
            goto LLTFin;
        }
        if (n == 1)
            eerr++;

        ctx->CTNCat[i] = n;
        if (!(ctx->CTCat[i] = (int *) calloc((size_t)(n + 1),sizeof(int)))) {   
            p_err(ctx, -2,1);
            goto LLTFin;
        }
        ctx->CTCatAI[i] = n + 1;
        memrq(ctx, n + 1,sizeof(int));

        ctx->NCTab *= n;

        printf1(ctx, "%3d  %c     %s ",i + 1,'A' + i,ctx->VName[ctx->PMVIdx[i]]);
        prnchar(ctx, ' ',ctx->VNameLen - (int)strlen(ctx->VName[ctx->PMVIdx[i]]),0);
        printf1(ctx, "  %3d : ",n);
        for (j = 1; j <= n; ++j) {
            k = ctx->AcS[ctx->AcJ[j]];
            ctx->CTCat[i][j] = k;
            rt_printf1_i(ctx, ctx->PMNFmtS,k);
        }
        newline(ctx);
    }
    newline(ctx);
    if (eerr) {
        printf1(ctx, "Error: only one value in at least one dimension.\n");
        goto LLTFin;
    }

    /*  build contingency table, value combinations in CTVal, frequencies
        in CTFreq, NDIM dimensions, NCTab cells */

    if (!(ctx->CTVal = (short *) calloc((size_t)(ctx->NCTab) * (size_t)(ctx->NDIM) + 1,sizeof(short)))) {   
        p_err(ctx, -2,1);
        goto LLTFin;
    }
    ctx->CTValA = ctx->NCTab * ctx->NDIM + 1;
    memrq(ctx, ctx->CTValA,sizeof(short));

    if (!(ctx->CTFreq = (int *) calloc((size_t)(ctx->NCTab + 1),sizeof(int)))) {   
        p_err(ctx, -2,1);
        goto LLTFin;
    }
    ctx->CTFreqA = ctx->NCTab + 1;
    memrq(ctx, ctx->CTFreqA,sizeof(int));

    if (ctx->PMSCAL >= 0) {
        if (!(ctx->CTSCAL = (double *) calloc((size_t)(ctx->NCTab + 1),sizeof(double)))) {   
            p_err(ctx, -2,1);
            goto LLTFin;
        }
        ctx->CTSCALA = ctx->NCTab + 1;
        memrq(ctx, ctx->CTSCALA,sizeof(double));
    }

    for (k = 0; k < ctx->NDIM; ++k)
        ctx->AcS[k] = 1;

    i = 0;
    while (1) {
        for (k = 0; k < ctx->NDIM; ++k)  
            ctx->CTVal[++i] = (short) ctx->CTCat[k][ctx->AcS[k]];

        j = 0;
        for (k = ctx->NDIM - 1; k >= 0; --k) {
            if ((ctx->AcS[k] += 1) <= ctx->CTNCat[k]) {
                j = 1;
                break;     
            }
            else
                ctx->AcS[k] = 1;
        }
        if (j == 0)
            break;
    }
    n = 1;
    for (j = 0; j < ctx->NOC; ++j) {

        if (ctx->PMWVar >= 0)  
            n = (int)get_data(ctx, ctx->PMWVar,j);

        for (i = 0; i < ctx->NDIM; ++i)  
            ctx->AcS[i] = (int)get_data(ctx, ctx->PMVIdx[i],j);

        l = m = 1;
        for (i = ctx->NDIM - 1; i >= 0; --i) {
            for (k = 1; k <= ctx->CTNCat[i]; ++k) {
                if (ctx->AcS[i] == ctx->CTCat[i][k])
                    break;
            }
            l += (k - 1) * m;
            m *= ctx->CTNCat[i];
        }
        ctx->CTFreq[l] += n;
        ctx->CTNum += n;
        if (ctx->PMSCAL >= 0)
            ctx->CTSCAL[l] += get_data(ctx, ctx->PMSCAL,j);
    }
    printf1(ctx, "Table has %d cells, %d counts.\n",ctx->NCTab,ctx->CTNum);
    if (ctx->CTNum < 1)  
        goto LLTFin;   

    ctx->CTIFlg = ctx->NCTab1 = 0;
    for (i = 1; i <= ctx->NCTab; ++i) {
        if (ctx->CTFreq[i] > 0)
            ctx->NCTab1++;
    }
    if (ctx->NCTab1 != ctx->NCTab) {       /* check if table is complete */
        ctx->CTIFlg = 1;
        printf1(ctx, "Table is incomplete:");
        printf1(ctx, " %d (%5.2f%%) of %d cells.\n",ctx->NCTab1,
                                 100.0 * (double)ctx->NCTab1/(double)ctx->NCTab,ctx->NCTab);
    }
    else 
        printf1(ctx, "Table is complete.\n");

    err = 0;

LLTFin:
    alloc_acr(ctx, 0);
    alloc_acs(ctx, 0);
    alloc_aci(ctx, 0);
    alloc_acj(ctx, 0);
    return(err);
}  

/* ------------------------------------------------------------------------ */
/*  loglin1     Estimation of log-linear models                             */
/*              Return 0 if OK, -1 if error                                 */

int loglin1(TDAContext *ctx)
{
    register int i,j,k,l;
    register char *q;
    int r,err,slen;
    double tmp;
#ifdef TDA_R_PACKAGE
    double crow[5];
    char clab[64];
#endif

    err = -1;

    printf1(ctx, "Begin of model estimation.\n\n");

    for (i = 0; i < ctx->PMModN; ++i) {

        ctx->LLDMCol = ctx->LLDCol[i];        /* columns of the design matrix */
        slen = 1;
        printf1(ctx, "Model: ");         /* print model description */
          
        if (ctx->LLTyp[i] == 0)          /* if null model */
            printf1(ctx, "0\n");
             
        else if (ctx->LLTyp[i] > 0) {

            prn_mdes(ctx, TDA_CONSOLE,ctx->PMModS[i],7);
            expand1(ctx, ctx->PMModS[i],ctx->LLD1,1);     
            r = expand2(ctx, ctx->LLD1,ctx->LLD2,1);     
            if (r > ctx->PMModSA[i]) {
                printf1(ctx, "Expanded: ");
                prn_mdes(ctx, TDA_CONSOLE,ctx->LLD2,10);
            }

            /*  build the design matrix. slen is the max length of the model
                parameter labels. */

            r = build_dm(ctx, ctx->LLD2,ctx->LLDMCol,&slen,1);
            if (r < 0 || r != ctx->LLDMCol) { 
                ll_error(ctx, r,ctx->PMModS[i]);
                goto LL1Fin;
            }
        }
        else if (ctx->LLTyp[i] < 0) { 
            slen = 10;
    
            printf1(ctx, "%s\n",ctx->PMModS[i]);

            r = ll_readdm(ctx, i,ctx->LLDMCol,1);

            if (r < 0) { 
                ll_error(ctx, r,ctx->PMModS[i]);
                goto LL1Fin;
            }
        }
        newline(ctx);     
            
        ctx->CTNPar = ctx->LLDMCol;   /* number of model parameters = number of
                               columns of the design matrix. */

        if (ctx->LLDMCol && ctx->PMF1Def)     /* if requested print design matrix */
            prn_lldmat(ctx, i);  

        printf1(ctx, "Number of model parameters: %d\n",ctx->CTNPar + 1); 
        printf1(ctx, "Maximum number of iterations: %d\n",ctx->MxIter);
        printf1(ctx, "Tolerance for convergence: %lg\n\n",ctx->TOLF);
            
        r = logl_fit(ctx);     /* call function for estimation */
         
        printf1(ctx, "Convergence ");
        if (r > ctx->MxIter)
            printf1(ctx, "not ");
        printf1(ctx, "reached in %d iterations.\n\n",r);

        if (ctx->CTLR < 1.e-6)
            ctx->CTLR = 0.0;
        if (ctx->CTCHI < 1.e-6)
            ctx->CTCHI = 0.0;

        printf1(ctx, "Likelihood Ratio Statistic %14.4f",ctx->CTLR);
        if (ctx->CTNDF > 0)
            printf1(ctx, "  Prob: %6.4f",1.0 - cdchif(ctx, ctx->CTLR,ctx->CTNDF));
        printf1(ctx, "\nPearson's Chi Square       %14.4f",ctx->CTCHI);
        if (ctx->CTNDF > 0)
            printf1(ctx, "  Prob: %6.4f",1.0 - cdchif(ctx, ctx->CTCHI,ctx->CTNDF));
        if (ctx->CTNDF > 0)
            printf1(ctx, "\nF-Statistic                %14.4f",ctx->CTLR / (double)ctx->CTNDF);
        printf1(ctx, "\nDegrees of Freedom %d\n\n",ctx->CTNDF);
#ifdef TDA_R_PACKAGE
        {
            double erow[6];
            erow[0] = ctx->CTLR;
            erow[2] = ctx->CTCHI;
            if (ctx->CTNDF > 0) {
                erow[1] = 1.0 - cdchif(ctx, ctx->CTLR,ctx->CTNDF);
                erow[3] = 1.0 - cdchif(ctx, ctx->CTCHI,ctx->CTNDF);
                erow[4] = ctx->CTLR / (double)ctx->CTNDF;
            }
            else
                erow[4] = (double)NAN;
        erow[3] = erow[4];
        erow[1] = erow[3];
            erow[5] = (double)ctx->CTNDF;
            tda_export_mat(ctx, "loglin.stats", erow, 1, 6);
        }
#endif
             
        /*  calculate constant in CTPar[0] */

        tmp = 0.0;
        for (j = 1; j <= ctx->CTNPar; ++j)  
            tmp += (double)ctx->LLDMat[j] * ctx->CTPar[j];

        ctx->CTPar[0] = rlog(ctx, ctx->CTFit[1]) - tmp;

        /*  print estimated parameters and errors */

        if (slen < 10)
            slen = 10;

        printf1(ctx, "Idx  Parameter");
        prnchar(ctx, ' ',slen - 9,0);
        prnchar(ctx, ' ',ctx->PMTFmt1 - 10,0); printf1(ctx, "     Coeff");
        prnchar(ctx, ' ',ctx->PMTFmt1 -  9,0); printf1(ctx, "     Error");
        prnchar(ctx, ' ',ctx->PMTFmt1 -  9,0); printf1(ctx, "    T-Stat  Signif\n");
        prnchar(ctx, '-',12 + slen + 3 * (ctx->PMTFmt1 + 1),0);

        q = ctx->CTVLab;
        l = 0;
        r = -1;
        for (j = 0; j <= ctx->LLDMCol; ++j) {
            printf1(ctx, "\n%3d  ",j);
            if (ctx->LLTyp[i] > 0) {
#ifdef TDA_R_PACKAGE
                size_t lb = 0;
#endif
                for (k = 0; k < slen; ++k) {
                    if (l < ctx->CTVLen && *q) {
#ifdef TDA_R_PACKAGE
                        if (lb < sizeof(clab) - 1)
                            clab[lb++] = *q;
#endif
                        printf1(ctx, "%c",*q++);
                        l++;
                    }
                    else
                        printf1(ctx, " ");
                }
#ifdef TDA_R_PACKAGE
                /* the parameter label the table prints, character by
                   character; it exists only as text, so the R side had
                   no way to name a coefficient without the console */
                clab[lb] = '\0';
                tda_export_str_row(ctx, "loglin.coeff.names", clab);
#endif
                l++;
                q++;
            }
            else if (j == 0) {
                printf1(ctx, "Constant  ");
#ifdef TDA_R_PACKAGE
                tda_export_str_row(ctx, "loglin.coeff.names", "Constant");
#endif
            }
            else {
                printf1(ctx, "Col%4d   ",j);
#ifdef TDA_R_PACKAGE
                snprintf(clab, sizeof(clab), "Col%d", j);
                tda_export_str_row(ctx, "loglin.coeff.names", clab);
#endif
            }

            r++;
            rt_printf1_d(ctx, ctx->PMTFmtS,ctx->CTPar[r]);
#ifdef TDA_R_PACKAGE
            crow[0] = (double)j;
            crow[1] = ctx->CTPar[r];
            crow[4] = (double)NAN;
        crow[3] = crow[4];
        crow[2] = crow[3];
#endif

            if (j && (tmp = ctx->CTMat[(r - 1) * ctx->CTNPar + r]) > 0.0) {
                tmp = sqrt(tmp);
                rt_printf1_d(ctx, ctx->PMTFmtS,tmp);
#ifdef TDA_R_PACKAGE
                crow[2] = tmp;
#endif
                tmp = ctx->CTPar[r] / tmp;
                rt_printf1_d(ctx, ctx->PMTFmtS,tmp);  
#ifdef TDA_R_PACKAGE
                crow[3] = tmp;
#endif
                tmp = 2.0 * cdnf(ctx, fabs(tmp)) - 1.0;
                printf1(ctx, " %6.4f",tmp);
#ifdef TDA_R_PACKAGE
                crow[4] = tmp;
#endif
            }
            else {
                prnchar(ctx, ' ',ctx->PMTFmt1 - 3,0); printf1(ctx, "--- ");
                prnchar(ctx, ' ',ctx->PMTFmt1 - 3,0); printf1(ctx, "---     ---");
            }
#ifdef TDA_R_PACKAGE
            tda_export_row(ctx, "loglin.coeff", crow, 5);
#endif
        }
        printf1(ctx, "\n\n");
#ifdef TDA_R_PACKAGE
        tda_export_flush(ctx, "loglin.coeff");
        tda_export_str_flush(ctx, "loglin.coeff.names");
#endif

        if (ctx->PMPPFDef)       /* write parameter estimates */
            ll_ppar(ctx);

        if (ctx->PMResFDef) {    /* if requested print table of residuals */

            fprintf(ctx->PMResFd,"Residuals of Model: ");
            if (ctx->LLTyp[i] > 0)
                prn_mdes(ctx, ctx->PMResFd,ctx->PMModS[i],20);
            else
                fprintf(ctx->PMResFd,"%s\n",ctx->PMModS[i]+2);

            prn_ctab(ctx, 2);
            printf1(ctx, "Residuals written to: %s\n",ctx->PMResFName);
        }
#ifdef TDA_R_PACKAGE
        /* one matrix per model, loglin.vcov, loglin.vcov.2, ... in mod= order */
        if (ctx->LLDMCol)
            export_prn(ctx, "loglin.vcov", ctx->CTNPar, ctx->CTNPar, ctx->CTNPar,
                       ctx->CTMat);
#endif
        if (ctx->LLDMCol && ctx->PMCovFDef) {    /* if requested print covariance matrix */
            prn_llcov(ctx, i);
            printf1(ctx, "Covariance matrix written to: %s\n",ctx->PMCovFName);
        }
        printf1(ctx, "\n");
    }
    err = 0;

LL1Fin:
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  ll_ppar()              Write parameter to PMPPFd.                       */

void ll_ppar(TDAContext *ctx)
{
    register int i;

    for (i = 0; i <= ctx->LLDMCol; ++i) {
        rt_fprintf_d(ctx, ctx->PMPPFd,ctx->PMTFmtS,ctx->CTPar[i]);
        fprintf(ctx->PMPPFd,"\n");
    }
    printf1(ctx, "Parameter estimates written to: %s\n",ctx->PMPPFName);
}

/* ------------------------------------------------------------------------ */
/*  ll_error    Print error message.                                        */

void ll_error(TDAContext *ctx, int err,char *mod)
{
    printf1(ctx, "Error: %s\n",mod);

    switch (err) {
        case  -1:  printf1(ctx, "Syntax error.\n");
                   break;
        case  -2:  printf1(ctx, "Can't find at least one dimension.\n");
                   break;
        case  -3:  printf1(ctx, "Can't read the level of a dimension.\n");
                   break;
        case  -4:  printf1(ctx, "At least one level is out of range.\n");
                   break;
        case  -5:  printf1(ctx, "Exceeded max number of table dimensions.\n");
                   break;
        case  -6:  printf1(ctx, "Exceeded max space for model expansion.\n");
                   break;
        case  -7:  printf1(ctx, "Error in build_dm() algorithm.\n");
                   break;
        case  -8:  printf1(ctx, "Can't open the file.\n");
                   break;
        case  -9:  printf1(ctx, "File is empty.\n");
                   break;
        case -10:  printf1(ctx, "File contains errors.\n");
                   break;
        case -11:  printf1(ctx, "One or more lines have less than %d entries.\n",ctx->LLDMCol);
                   break;
        case -12:  printf1(ctx, "Found more than %d lines in the file.\n",ctx->NCTab);
                   break;
        default:   break;
    }
}

/* ------------------------------------------------------------------------ */
/*  prn_mdes     Print model description desc to file fd.                   */

void prn_mdes(TDAContext *ctx, FILE *fd, char *desc, int n)
{
    register char *p;

    (void)ctx;
    (void)n;        /* unused: the signature is shared */

    p = desc;
    while (*p)  
        fprintf(fd,"%c",*p++);
    fprintf(fd,"\n");
}

/* ------------------------------------------------------------------------ */
/*  prn_lldmat     Print design matrix of model i to PMF1d.                 */

void prn_lldmat(TDAContext *ctx, int i)
{
    register int j,k;

    fprintf(ctx->PMF1d,"# Design matrix of model: ");
    if (ctx->LLTyp[i] > 0)
        prn_mdes(ctx, ctx->PMF1d,ctx->PMModS[i],24);
    else if (ctx->LLTyp[i] < 0)
        fprintf(ctx->PMF1d,"%s\n",ctx->PMModS[i]+2);
   
    for (j = 1; j <= ctx->NCTab; ++j) {
        for (k = 1; k <= ctx->LLDMCol; ++k)
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,(int)ctx->LLDMat[(j - 1) * ctx->LLDMCol + k]);
        fprintf(ctx->PMF1d,"\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  prn_llcov      Print covariance matrix of model i to PMCovFd.           */

void prn_llcov(TDAContext *ctx, int i)
{
    register int j,k;

    fprintf(ctx->PMCovFd,"# Covariance matrix of model: ");
    if (ctx->LLTyp[i] > 0)
        prn_mdes(ctx, ctx->PMCovFd,ctx->PMModS[i],28);
    else
        fprintf(ctx->PMCovFd,"%s\n",ctx->PMModS[i]+2);

    for (j = 1; j <= ctx->CTNPar; ++j) {
        for (k = 1; k <= ctx->CTNPar; ++k)  
            rt_fprintf_d(ctx, ctx->PMCovFd,ctx->PMMFmtS,ctx->CTMat[(j - 1) * ctx->CTNPar + k]);
        fprintf(ctx->PMCovFd,"\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  prn_ctab (mode)                                                         */
/*      mode 1 : Print contingency table to PMTabFd                         */
/*      mode 2 : Print residuals to PMResFd                                 */

void prn_ctab(TDAContext *ctx, int mode)
{
    FILE *ofd;
    register int j,k;
    double s = 1.0;
#ifdef TDA_R_PACKAGE
    double *erow = (double *)malloc((size_t)(ctx->NDIM + 6) * sizeof(double));
#endif

    if (mode == 1)
        ofd = ctx->PMTabFd;
    else
        ofd = ctx->PMResFd;

    if (mode == 2) {
        fprintf(ofd,"\nIndex  ");
        for (k = 0; k < ctx->NDIM; ++k) {
            fprnchar(ctx, ofd,' ',ctx->PMNFmt - 1,0);
            fprintf(ofd,"%c ",'A' + k);
        }
        fprintf(ofd,"  Observed ");
        fprintf(ofd,"      Fitted ");
        fprintf(ofd,"    Residual        Scale\n");
        fprnchar(ctx, ofd,'-',56 + ctx->NDIM * (ctx->PMNFmt + 1),1);
    }
    for (k = 1; k <= ctx->NCTab; ++k) {
        fprintf(ofd,"%5d  ",k);
        for (j = 1; j <= ctx->NDIM; ++j)  
            rt_fprintf_i(ctx, ofd,ctx->PMNFmtS,(int)ctx->CTVal[(k - 1) * ctx->NDIM + j]);
        fprintf(ofd,"%10d ",ctx->CTFreq[k]);

        if (mode == 2)  
            fprintf(ofd,"%12.4f %12.4f ",ctx->CTFit[k],ctx->CTRes[k]);

        if (ctx->PMSCAL >= 0)
            s = ctx->CTSCAL[k];

        if (ctx->PMSCAL >= 0 || mode == 2)
            fprintf(ofd,"%12.4f",s);

        fprintf(ofd,"\n");
#ifdef TDA_R_PACKAGE
        if (erow != NULL) {
            int e = 0;
            erow[e++] = (double)k;
            for (j = 1; j <= ctx->NDIM; ++j)
                erow[e++] = ctx->CTVal[(k - 1) * ctx->NDIM + j];
            erow[e++] = (double)ctx->CTFreq[k];
            if (mode == 2) {
                erow[e++] = ctx->CTFit[k];
                erow[e++] = ctx->CTRes[k];
            }
            if (ctx->PMSCAL >= 0 || mode == 2)
                erow[e++] = s;
            tda_export_row(ctx, mode == 2 ? "loglin.residuals"
                                          : "loglin.table", erow, e);
        }
#endif
    }
    if (mode == 2)  
        fprintf(ofd,"\n");
#ifdef TDA_R_PACKAGE
    free(erow);
    tda_export_flush(ctx, mode == 2 ? "loglin.residuals" : "loglin.table");
#endif
}

/* ------------------------------------------------------------------------ */
/*  llmodcheck   Check model description.                                   */
/*      Return:  1  if no error                                             */
/*               0  if empty string (null model)                            */
/*              -1  if syntax error                                         */
/*              -2  if one or more dimensions not found                     */
/*              -3  if a level cannot be read                               */
/*              -4  if a level is out of range                              */

int llmodcheck(TDAContext *ctx, char *desc)
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
            if (k >= ctx->NDIM)
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
                for (l = 1; l < ctx->CTNCat[k]; ++l) {
                    if (n == ctx->CTCat[k][l]) {
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

int expand1(TDAContext *ctx, char *s, char *r, int mode)
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
                    while (!(last = logscr_comb(ctx, iset,j,i,last))) {
                        for (k = 1; k <= i; ++k) {
                            if (k > 1) {
                                if (mode)
                                    *r++ = '.';
                                len++;
                            }
                            if (mode)
                                *r++ = (char)('A' + dim[iset[k]]);
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

int expand2(TDAContext *ctx, char *s, char *r, int mode)
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
                ctx->ExMult[n] = 1;
                ctx->ExOpl[n] = (char)lvl;
                if (++n > ctx->ExMax)
                    return(-6);
            }
            else if (t == '+') {
                ctx->ExMult[n] = 0;
                ctx->ExOpl[n] = (char)lvl;
                if (lvl == 0) {
                    if (plus) {
                        if (mode)
                            *r++ = '+';
                        len++;
                    }
                    r = distribute(ctx, n,r,mode,&len);        
                    plus++;
                    n = 0;
                }
                else if (++n > ctx->ExMax)
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
                ctx->ExV[n] = tvar;
                ctx->ExVTyp[n] = (char)(ttyp);
                ctx->ExVc[n] = tvarc;
                ctx->ExVl[n] = (char)lvl;
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

char *distribute(TDAContext *ctx, int n, char *r, int mode, int *len)
{
    register int i,j,k,l,alter,level,product,term,plus = 0;
    char sstmp[80];

    for (k = 0; k <= n; ++k)
        ctx->ExUsed[k] = 0;
NEXT:
    alter = 1;
    j = i = -1;
FACTOR:
    if (ctx->ExUsed[++i])
        goto FACTOR;
    ctx->ExIdx[++j] = i;

    while (1) {
        if (ctx->ExMult[i])
            goto FACTOR;
        level = (int)ctx->ExOpl[i];
     
        if (level <= 0)
            break;
        else {
            if (alter) {
                l = level;
                level = (int)ctx->ExVl[i] + 1;
                ctx->ExUsed[i] = (char)(product = term = 1);
                alter = 0;
                for (k = i - 1; k >= 0; k--) {
                    if ((int)ctx->ExOpl[k] < level) {
                        level = (int)ctx->ExOpl[k];
                        product = ctx->ExMult[k];
                        if (product)
                            level++;
                        if (level <= l)
                            term = 0;
                    }
                    if (product)
                        ctx->ExUsed[k] = (char)term;
                }
            }
            else {
                while (1) {
                    if (level > (int)ctx->ExOpl[++i])
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
            *r++ = ctx->ExV[ctx->ExIdx[k]];

        if (ctx->ExVTyp[ctx->ExIdx[k]] == 2) {
            snprintf(sstmp,sizeof(sstmp),"[%d]",ctx->ExVc[ctx->ExIdx[k]]);
            *len += (int)strlen(sstmp);
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

int build_dm(TDAContext *ctx, char *descr, int col, int *slen, int mode)
{
    register int i,j,k,l,m,m1,m2;
    int len,nsl,fnd,r,n,dn,x,ptr[MaxDX + 1];
    int cat[MaxDX + 1],ctyp[MaxDX + 1],clev[MaxDX + 1];
    register char *p,*q,*ctvptr = NULL;
    char sstmp[20];

    if (mode) {
        ctvptr = ctx->CTVLab;
        *ctvptr = '\0';
        strcat(ctvptr,"Constant");
        ctvptr += 9;
        *ctvptr = '\0';
        ctx->CTVLen = 9;
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
                for (l = 1; l < ctx->CTNCat[k]; ++l) {
                    if (n == ctx->CTCat[k][l]) {
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
                m *= ctx->CTNCat[k] - 1;
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
                        if (clev[k] == ctx->CTCat[m][m2])  
                            n++;
                    }
                }
                if (!nsl || nsl == n) {

                    j++;
                    if (j > col)  
                        return(-7);
                    len = 0;
    
                    for (i = 1; i <= ctx->NCTab; ++i) {

                        n = 1;
                        for (k = 1; k <= dn; ++k) {
                            if (n) {
                                m = ptr[k];
                                x = (int) ctx->CTVal[(i - 1) * ctx->NDIM + m + 1];
                                m1 = ctx->CTNCat[m];
                                m2 = cat[k];
                                if (x == ctx->CTCat[m][m1])
                                    n *= -1;
                                else if (x != ctx->CTCat[m][m2])
                                    n = 0;
                            }
                                /* create label of parameter */

                            if (i == 1 && ctx->CTVLen < ctx->CTVLenM) {
                                m = ptr[k];
                                m2 = cat[k];
                                if (k > 1) {
                                    *ctvptr++ = '.';
                                    ctx->CTVLen++;
                                    len++;
                                }
                                if (ctx->CTVLen < ctx->CTVLenM) {
                                    *ctvptr++ = 'A' + (char)m;
                                    len++;
                                    ctx->CTVLen++;
                                    snprintf(sstmp,sizeof(sstmp),"[%d]",ctx->CTCat[m][m2]);
                                    q = sstmp;
                                    while (*q && ctx->CTVLen < ctx->CTVLenM) {
                                        *ctvptr++ = *q++;
                                        len++;
                                        ctx->CTVLen++;
                                    }
                                }
                            }
                        }
                        ctx->LLDMat[(i - 1) * col + j] = (short) n;
                    }
                    if (ctx->CTVLen < ctx->CTVLenM) {
                        *ctvptr++ = '\0';
                        ctx->CTVLen++;
                        if (*slen < len)
                            *slen = len;
                    }
                }
                if (design_comb(ctx, dn,cat,ptr))   /* next combination */
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

int design_comb(TDAContext *ctx, int dn, int *cat, int *ptr) 
{
    register int j;

    for (j = dn; j >= 1; --j) {
        if ((cat[j] += 1) < ctx->CTNCat[ptr[j]])
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

int ll_readdm(TDAContext *ctx, int i, int dcol, int mode)
{
    FILE *fd;
    register int col,row;
    int n;
    register char *p;

    if (alloc_acc(ctx, RLMaxDef + 1))    
        return(-2);   

    if (!(fd = fopen(ctx->PMModS[i] + 2,OPEN_RD)))   
        return(-8);
    else {
        row = col = 0;
        while (fgets(ctx->AcC,RLMaxDef,fd)) {
            p = skip_b(ctx, ctx->AcC);

            col = 0;
            if ((p = check_comment(ctx, ctx->AcC)) != NULL) {   /* if no comment */

                while (*p && *p != '\n') {    
                    if (sscanf(p,"%d",&n) != 1) { 
                        fclose(fd);
                        return(-10);
                    }
                    if (!col)
                        row++;
                    col++;
                    if (mode && col <= dcol) {
                        if (row > ctx->NCTab) {
                            fclose(fd);
                            return(-12);
                        }
                        ctx->LLDMat[(row - 1) * dcol + col] = (short) n; 
                    }
                    p = skip_int(ctx, p);
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

int logl_fit(TDAContext *ctx) 
{
    register int i,k,l,ii,ll;
    int it,irank = 0;
    double tmp,tmp1,tmp2,e,f,sum,sum1,w;

    /* initialize the algorithm */

    f = ctx->DBLMAX;
    ctx->CTNDF = ctx->NCTab - 1;
    tmp = sum = sum1 = 0.0;
    for (i = 1; i <= ctx->NCTab; ++i) {
        ctx->CTFit[i] = (double)ctx->CTFreq[i] + 0.5;
        ctx->CTRes[i] = rlog(ctx, ctx->CTFit[i]);

        if (ctx->PMSCAL >= 0) {
            if (ctx->CTSCAL[i] > ctx->EPSI1)      
                ctx->CTRes[i] -= rlog(ctx, (double)ctx->CTSCAL[i]);
            else {
                ctx->CTFit[i] = 0.0;
                ctx->CTNDF--;
            }
        }
        sum += ctx->CTFit[i];
        sum1 += ctx->CTFit[i] * ctx->CTRes[i];
    }
    if (sum > ctx->TOLF) {
        w = sum1 / sum; 
        for (i = 1; i <= ctx->NCTab; ++i)
            ctx->CTRes[i] -= w;
    }
    for (k = 1; k <= ctx->CTNPar; ++k) {
        ctx->CTPar[k] = sum = 0.0;
        for (i = 1; i <= ctx->NCTab; ++i)   
            sum += ctx->CTFit[i] * ctx->CTRes[i] * (double)ctx->LLDMat[(i - 1) * ctx->LLDMCol + k];
        ctx->CTWork[k] = sum;
    }

    /* perform an iteration, maximal MxIter iterations */

    for (it = 1; it <= ctx->MxIter; ++it) {
        /*********************
        if (it > 1 || init1) {
        *********************/
            sum = 0.0;
            for (i = 1; i <= ctx->NCTab; ++i)
                sum += ctx->CTFit[i];
            w = sum;
            for (k = 1; k <= ctx->CTNPar; ++k) {
                sum = 0.0;
                for (i = 1; i <= ctx->NCTab; ++i)  
                    sum += ctx->CTFit[i] * (double)ctx->LLDMat[(i - 1) * ctx->LLDMCol + k];
                if (w > ctx->TOLF)
                    ctx->CTB[k] = sum / w;
                else
                    ctx->CTB[k] = 0.0;
            }

            /* obtain weighted sums of cross-products */

            for (k = 1; k <= ctx->CTNPar; ++k) {
                for (l = 1; l <= k; ++l) {
                    sum = 0.0;
                    for (i = 1; i <= ctx->NCTab; ++i)  
                        sum += ((double)ctx->LLDMat[(i - 1) * ctx->LLDMCol + k] - ctx->CTB[k]) *
                               ((double)ctx->LLDMat[(i - 1) * ctx->LLDMCol + l] - ctx->CTB[l]) *
                                                                      ctx->CTFit[i]; 
                    ctx->CTMat[(k - 1) * ctx->CTNPar + l] = sum;
                }
            }

            /* obtain Cholesky decomposition of s */

            irank = logl_chol(ctx, ctx->CTNPar,ctx->CTMat,ctx->CTMat);

            /* see if further steps are needed */

            if (f < ctx->TOLF) {
                logl_inv(ctx, ctx->CTNPar,ctx->CTMat,ctx->CTMat);
                ctx->CTNDF -= irank;
                goto LF_Fin;
            }

            /* obtain new value of b */

            logl_solve(ctx, ctx->CTNPar,ctx->CTMat,ctx->CTWork,ctx->CTWork);
            for (k = 1; k <= ctx->CTNPar; ++k) 
                ctx->CTPar[k] += ctx->CTWork[k];

            /* update U and CTFit and check for convergence */

            if (it > 1)
                f = 0.0;
        /**
        }
        **/
        for (i = 1; i <= ctx->NCTab; ++i) {
            sum = 0.0;
            if (ctx->CTNPar > 0) {
                for (k = 1; k <= ctx->CTNPar; ++k)
                    sum += ctx->CTPar[k] * (double)ctx->LLDMat[(i - 1) * ctx->LLDMCol + k];
                if (it > 1) {
                    e = fabs(sum - ctx->CTRes[i]);   
                    if (e > f)
                        f = e;
                }
            }
            ctx->CTRes[i] = sum;
            ctx->CTFit[i] = rexp(ctx, ctx->CTRes[i]);        

            if (ctx->PMSCAL >= 0) {
                if (ctx->CTSCAL[i] > ctx->EPSI1)
                    ctx->CTFit[i] *= ctx->CTSCAL[i];
                else
                    ctx->CTFit[i] = 0.0;         
            }
        }
        sum = w = 0.0;
        for (i = 1; i <= ctx->NCTab; ++i) {
            if (ctx->PMSCAL < 0 || ctx->CTSCAL[i] > ctx->EPSI1) {
                w += (double)ctx->CTFreq[i];
                sum += ctx->CTFit[i];
            }
        }
        if (sum > ctx->TOLF)
            w /= sum;
        else
            w = 0.0;
        for (i = 1; i <= ctx->NCTab; ++i)
            ctx->CTFit[i] *= w;

        if (ctx->CTNPar <= 0)
            goto LF_Fin;

        /* prepare differences between fitted and observed lin. combinations */

        for (k = 1; k <= ctx->CTNPar; ++k) {
            sum = 0.0;

            for (i = 1; i <= ctx->NCTab; ++i) {

                if (ctx->CTFit[i] > ctx->TOLF) {
                    tmp = (double)ctx->CTFreq[i];
                    sum += (tmp - ctx->CTFit[i]) * (double)ctx->LLDMat[(i - 1) * ctx->LLDMCol + k];
                }
            }
            ctx->CTWork[k] = sum;
        }
    }
    logl_inv(ctx, ctx->CTNPar,ctx->CTMat,ctx->CTMat);
    ctx->CTNDF -= irank;
    
    /* calc chi-square and likelihood ratio statistics */

LF_Fin:
    ctx->CTLR = ctx->CTCHI = 0.0;
    for (i = 1; i <= ctx->NCTab; ++i) {

        if (ctx->CTFit[i] > ctx->TOLF) {
            tmp1 = (double)ctx->CTFreq[i];
            tmp2 = tmp1 - ctx->CTFit[i];
            ctx->CTCHI += tmp2 * tmp2 / ctx->CTFit[i];
            if (tmp1 > ctx->TOLF)  
                ctx->CTLR += tmp1 * rlog(ctx, tmp1 / ctx->CTFit[i]);       
        }
    }
    ctx->CTLR *= 2.0;

    /* calculation of adjusted residuals in array CTRes */

    sum1 = 0.0;
    for (i = 1; i <= ctx->NCTab; ++i)
        sum1 += ctx->CTFit[i];

    for (i = 1; i <= ctx->NCTab; ++i) {
        ii = (i - 1) * ctx->LLDMCol;
        sum = 0.0;
        if (ctx->CTNPar > 1) {
            for (k = 2; k <= ctx->CTNPar; ++k) {
                ll = k - 1;
                for (l = 1; l <= ll; ++l)
                sum += ctx->CTMat[(k - 1) * ctx->CTNPar + l] * 
                                            ((double)ctx->LLDMat[ii + k] - ctx->CTB[k]) *
                                            ((double)ctx->LLDMat[ii + l] - ctx->CTB[l]);
            }
            sum *= 2.0;
        }
        for (k = 1; k <= ctx->CTNPar; ++k) {
            tmp  = (double)ctx->LLDMat[ii + k] - ctx->CTB[k];
            sum += ctx->CTMat[(k - 1) * ctx->CTNPar + k] * tmp * tmp;
        }
        sum = ctx->CTFit[i] * (1.0 - ctx->CTFit[i] * sum);
 
        if (sum1 > ctx->TOLF)
            sum -= ctx->CTFit[i] * ctx->CTFit[i] / sum1;
        ctx->CTRes[i] = 0.0;

        if (sum > ctx->TOLF) {
            tmp = (double)ctx->CTFreq[i];
            ctx->CTRes[i] = (tmp - ctx->CTFit[i]) / sqrt(sum);   
        }
    }
    return(it);     /* return number of iterations */
}

/* ------------------------------------------------------------------------ */
/*  logl_solve (n,lu,b,x)                                                   */
/*      Solve lu * x = b with lu(i,i) != 0                                  */
/*      n is the dimension of the matrices.                                 */

void logl_solve(TDAContext *ctx, int n, double *lu, double *b, double *x)
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
    if (lu[n * n] > ctx->TOLF) 
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

void logl_inv(TDAContext *ctx, int n, double *tri, double *xinv)
{
    register int i,i1,k,j;
    int k1,in,jn;
    double sum;

    for (i = 1; i <= n; ++i) {
        in = (i - 1) * n;
        if (tri[in + i] > ctx->TOLF)
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

int logl_chol(TDAContext *ctx, int n, double *sym, double *tri)
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
            if (tri[1] > ctx->TOLF)
                tri[in + 1] = tri[i] / tri[1];
            for (j = 2; j <= i; ++j) {
                jn = (j - 1) * n;
                if (tri[jn + j] > ctx->TOLF) {
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
void logl_gen(TDAContext *ctx, double *c, double *obs, double *expv, double *res)
{
    register int i,k,l;
    int ll;
    double sum,sum1,sum2,st;

    sum = sum1 = 0.0;
    for (i = 1; i <= ctx->NCTab; ++i) {
        sum1 += ctx->CTFit[i] * c[i];
        sum  += ctx->CTFit[i];
    }
    sum2 = sum1 * sum1 / sum;

    if (ctx->CTNPar > 0) {
        for (k = 1; k <= ctx->CTNPar; ++k) {
            sum = 0.0;
            for (i = 1; i <= ctx->NCTab; ++i)  
                sum += ((double)ctx->LLDMat[(i-1) * ctx->LLDMCol + k] - ctx->CTB[k]) *
                                                             ctx->CTFit[i] * c[i];
            ctx->CTWork[k] = sum;
        }
        sum = 0.0;
        for (k = 2; k <= ctx->CTNPar; ++k) {
            ll = k - 1;
            for (l = 1; l <= ll; ++l)
                sum += ctx->CTMat[(k - 1) * ctx->CTNPar + l] * ctx->CTWork[k] * ctx->CTWork[l];
        }
        sum += sum;
        for (k = 1; k <= ctx->CTNPar; ++k) 
            sum += ctx->CTMat[(k - 1) * ctx->CTNPar + k] * ctx->CTWork[k] * ctx->CTWork[k];
    }
    sum1 = 0.0;
    for (i = 1; i <= ctx->NCTab; ++i)  
        sum1 += ctx->CTFit[i] * c[i] * c[i];
    st = sum1 - sum - sum2;
    *res = sum = sum1 = 0.0;
    for (i = 1; i <= ctx->NCTab; ++i) {
        sum += (double)ctx->CTFreq[i] * c[i];
        sum1 +=  ctx->CTFit[i] * c[i];
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

int logscrn(TDAContext *ctx) 
{
int last;
    register int i,j,m,jn,jn1;
    int err,df,iemax,mm,np,np1,ntab2,l3;
    short itp;
    double g21,g22,g23;

    err = -1;

    printf1(ctx, "Screening the table.\n");

    /*  calculate iemax, i.e. the maximal number of interaction effects
        in the saturated model. */

    ctx->NDIM1 = ctx->NDIM - 1;
    iemax = 1;
    j = ctx->NDIM / 2;
    for (i = 1; i <= j; ++i)
        iemax = (iemax * (ctx->NDIM - i + 1)) / i;

    if (!(ctx->CTFit = (double *) calloc((size_t)(ctx->NCTab + 1),sizeof(double)))) {
        p_err(ctx, -2,1);
        goto LSFin;
    }
    ctx->CTFitA = ctx->NCTab + 1;
    memrq(ctx, ctx->CTFitA,sizeof(double));

    if (!(ctx->GSQ  = (double *)calloc((size_t)(ctx->NDIM + 1),sizeof(double)))) {  
        p_err(ctx, -2,1);
        goto LSFin;
    }
    ctx->GSQA = ctx->NDIM + 1;
    memrq(ctx, ctx->GSQA,sizeof(double));

    if (!(ctx->PART = (double *)calloc((size_t)(iemax) * (size_t)(ctx->NDIM1) + 1,sizeof(double)))) {
        p_err(ctx, -2,1);
        goto LSFin;
    }
    ctx->PARTA = iemax * ctx->NDIM1 + 1;
    memrq(ctx, ctx->PARTA,sizeof(double));

    if (!(ctx->MARG = (double *)calloc((size_t)(iemax) * (size_t)(ctx->NDIM1) + 1,sizeof(double)))) {
        p_err(ctx, -2,1);
        goto LSFin;
    }
    ctx->MARGA = iemax * ctx->NDIM1 + 1;
    memrq(ctx, ctx->MARGA,sizeof(double));

    if (!(ctx->DGFR = (int *)calloc((size_t)(ctx->NDIM + 1),sizeof(int)))) {
        p_err(ctx, -2,1);
        goto LSFin;
    }
    ctx->DGFRA = ctx->NDIM + 1;
    memrq(ctx, ctx->DGFRA,sizeof(int));

    if (!(ctx->DFS  = (int *)calloc((size_t)(iemax) * (size_t)(ctx->NDIM1) + 1,sizeof(int)))) {
        p_err(ctx, -2,1);
        goto LSFin;
    }
    ctx->DFSA = iemax * ctx->NDIM1 + 1;
    memrq(ctx, ctx->DFSA,sizeof(int));

    j = iemax * ctx->NDIM1;
    for (i = 1; i <= j; ++i) {
        ctx->MARG[i] = ctx->PART[i] = -77.77;
        ctx->DFS[i] = -9999;
    }
    ntab2 = ctx->NCTab / 2 + 1;      

    /*  calculate mm, i.e. maximum of (i * NPi) (i = 1,NDIM) to allocate
        memory. */

    mm = 1;
    m = ctx->NDIM1 / 2;
    for (i = 1; i <= m; ++i)
        mm = (mm * (ctx->NDIM - i)) / i;
    mm *= ctx->NDIM;

    /*  allocate memory for some auxiliary arrays */

    if (!(ctx->LSISet = (int   *)calloc((size_t)(ctx->NDIM + 1),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto LSFin;
    }
    ctx->LSISetA = ctx->NDIM + 1;
    memrq(ctx, ctx->LSISetA,sizeof(int));

    if (!(ctx->LSJSet = (int   *)calloc((size_t)(ctx->NDIM + 1),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto LSFin;
    }
    ctx->LSJSetA = ctx->NDIM + 1;
    memrq(ctx, ctx->LSJSetA,sizeof(int));

    if (!(ctx->LSdim  = (short *)calloc((size_t)(ctx->NDIM + 1),sizeof(short)))) { 
        p_err(ctx, -2,1);
        goto LSFin;
    }
    ctx->LSdimA = ctx->NDIM + 1;
    memrq(ctx, ctx->LSdimA,sizeof(short));

    if (!(ctx->LSsiz  = (short *)calloc((size_t)(ctx->NDIM + 1),sizeof(short)))) {  
        p_err(ctx, -2,1);
        goto LSFin;
    }
    ctx->LSsizA = ctx->NDIM + 1;
    memrq(ctx, ctx->LSsizA,sizeof(short));

    if (!(ctx->LScord = (short *)calloc((size_t)(ctx->NDIM + 1),sizeof(short)))) {
        p_err(ctx, -2,1);
        goto LSFin;
    }
    ctx->LScordA = ctx->NDIM + 1;
    memrq(ctx, ctx->LScordA,sizeof(short));

    if (!(ctx->LSconf = (short *)calloc((size_t)(ctx->NDIM) * (size_t)(iemax) + 1,sizeof(short)))) {
        p_err(ctx, -2,1);
        goto LSFin;
    }
    ctx->LSconfA = ctx->NDIM * iemax + 1;
    memrq(ctx, ctx->LSconfA,sizeof(short));

    if (!(ctx->LSip = (short *)calloc((size_t)(ctx->NDIM1) * (size_t)(iemax) + 1,sizeof(short)))) {
        p_err(ctx, -2,1);
        goto LSFin;
    }
    ctx->LSipA = ctx->NDIM1 * iemax + 1;
    memrq(ctx, ctx->LSipA,sizeof(short));

    if (!(ctx->LSim = (short *)calloc((size_t)(ctx->NDIM1) * (size_t)(mm) + 1,sizeof(short)))) {  
        p_err(ctx, -2,1);
        goto LSFin;
    }
    ctx->LSimA = ctx->NDIM1 * mm + 1;
    memrq(ctx, ctx->LSimA,sizeof(short));

    if (!(ctx->LSx = (float *)calloc((size_t)(ntab2 + 1),sizeof(float)))) {  
        p_err(ctx, -2,1);
        goto LSFin;
    }
    ctx->LSxA = ntab2 + 1;
    memrq(ctx, ctx->LSxA,sizeof(float));

    if (!(ctx->LSy = (float *)calloc((size_t)(ntab2 + 1),sizeof(float)))) {
        p_err(ctx, -2,1);
        goto LSFin;
    }
    ctx->LSyA = ntab2 + 1;
    memrq(ctx, ctx->LSyA,sizeof(float));

    /*  change order of dimensions in the contingency table to account
        for the different array handling in Fortran and C. */

    for (i = 1; i <= ctx->NDIM; ++i)
        ctx->LSdim[i] = (short)(ctx->CTNCat[ctx->NDIM - i]);

    ctx->DGFR[ctx->NDIM] = ctx->NCTab - 1;

    logscr_res(ctx);
    ctx->GSQ[1] = logscr_lr(ctx);

    for (m = 1; m <= ctx->NDIM1; ++m) {
        np = logscr_conf(ctx, m,ctx->LSISet,ctx->LSJSet);
        logscr_res(ctx);
        ctx->DGFR[m] = logscr_eval(ctx, ctx->LSip,np,m,1);
        logscr_fit(ctx, np);
        ctx->GSQ[m + 1] = logscr_lr(ctx);

        for (i = 1; i <= m; ++i) {
            itp = ctx->LSip[i];
            np1 = np - 1;
            for (j = 1; j <= np1; ++j)
                ctx->LSip[(j - 1) * ctx->NDIM1 + i] = ctx->LSip[j * ctx->NDIM1 + i];
            ctx->LSip[(np - 1) * ctx->NDIM1 + i] = (short)itp;
        }
        l3 = -m + 1;
        for (j = 1; j <= np; ++j) {
            jn = (j - 1) * ctx->NDIM1;
            jn1 = (np - j) * ctx->NDIM1;   /* change of order */

            logscr_res(ctx);
            df = logscr_eval(ctx, ctx->LSip,np - 1,m,1);
            logscr_fit(ctx, np - 1);
            g21 = logscr_lr(ctx);

            ctx->DFS[jn1 + m] = ctx->DGFR[m] - df;
            ctx->PART[jn1 + m] = g21 - ctx->GSQ[m + 1];

            if (m == 1) {
                ctx->MARG[jn1 + 1] = ctx->PART[jn1 + 1];
            }
            else {
                logscr_res(ctx);
                df = logscr_eval(ctx, ctx->LSip,1,m,np);
                logscr_fit(ctx, 1);
                g22 = logscr_lr(ctx);
                l3 += m;
                logscr_res(ctx);
                df = logscr_eval(ctx, ctx->LSim,m,m - 1,l3);
                logscr_fit(ctx, m);
                g23 = logscr_lr(ctx);
                ctx->MARG[jn1 + m] = g23 - g22;
            }
            for (i = 1; i <= m; ++i) {
                itp = (short)((int)ctx->LSip[(np - 1) * ctx->NDIM1 + i]);
                ctx->LSip[(np - 1) * ctx->NDIM1 + i] = ctx->LSip[jn + i];
                ctx->LSip[jn + i] = (short)((int)itp);
            }
        }
        ctx->DGFR[ctx->NDIM] -= ctx->DGFR[m];
        ctx->GSQ[m] -= ctx->GSQ[m + 1];
    }

    /*  print test statistics for all interaction of order ... are zero */

    printf1(ctx, "Test for zero interactions.\n\n");
    printf1(ctx, "Order   DF  ");
    prnchar(ctx, ' ',ctx->PMTFmt1 - 9,0);  printf1(ctx, "Statistic");
    prnchar(ctx, ' ',ctx->PMTFmt1 - 4,0); printf1(ctx, "Prob \n");
    prnchar(ctx, '-',11 + 2 * (ctx->PMTFmt1 + 1),0);
    for (i = 1; i <= ctx->NDIM; ++i) {
        printf1(ctx, "\n%4d %5d  ",i,ctx->DGFR[i]);
        rt_printf1_d(ctx, ctx->PMTFmtS,ctx->GSQ[i]);
        rt_printf1_d(ctx, ctx->PMTFmtS,1.0 - cdchif(ctx, ctx->GSQ[i],ctx->DGFR[i]));
    }

    /*  print test statistics for marginal and partial association */

    printf1(ctx, "\n\nMarginal and partial association.\n");
    printf1(ctx, "\nModel");
    if ((mm = 3 * ctx->NDIM1) < 5)
        mm = 5;
    prnchar(ctx, ' ',mm - 5,0);
    printf1(ctx, "  DF");
    prnchar(ctx, ' ',ctx->PMTFmt1 - 8,0); printf1(ctx, "Marginal");
    prnchar(ctx, ' ',ctx->PMTFmt1 - 3,0); printf1(ctx, "Prob");
    prnchar(ctx, ' ',ctx->PMTFmt1 - 6,0); printf1(ctx, "Partial");
    prnchar(ctx, ' ',ctx->PMTFmt1 - 3,0); printf1(ctx, "Prob\n");
    prnchar(ctx, '-',3 * ctx->NDIM1 + 4 * (ctx->PMTFmt1 + 1) + 3,0);

    for (i = 1; i <= ctx->NDIM1; ++i) {
        last = 1;
        m = i;
        while (!(last = logscr_comb(ctx, ctx->LSISet,ctx->NDIM,i,last))) {
            printf1(ctx, "\n");
            for (j = 1; j <= i; ++j)
                printf1(ctx, " %c ",'A' + (char)(ctx->LSISet[j] - 1));
            prnchar(ctx, ' ',mm - 3 * i,0);
            printf1(ctx, "%4d",ctx->DFS[m]); 
            rt_printf1_d(ctx, ctx->PMTFmtS,ctx->MARG[m]);
            rt_printf1_d(ctx, ctx->PMTFmtS,1.0 - cdchif(ctx, ctx->MARG[m],ctx->DFS[m]));
            rt_printf1_d(ctx, ctx->PMTFmtS,ctx->PART[m]);
            rt_printf1_d(ctx, ctx->PMTFmtS,1.0 - cdchif(ctx, ctx->PART[m],ctx->DFS[m]));
            m += ctx->NDIM1;
        }
    }
    printf1(ctx, "\n\n");
    err = 0;

LSFin:
    if (ctx->CTFitA > 0) {
        free((char *)ctx->CTFit);
        memrq(ctx, -ctx->CTFitA,sizeof(double));
        ctx->CTFitA = 0;    
    }
    if (ctx->GSQA > 0) {
        free((char *)ctx->GSQ);
        memrq(ctx, -ctx->GSQA,sizeof(double));
        ctx->GSQA = 0;    
    }
    if (ctx->PARTA > 0) {
        free((char *)ctx->PART);
        memrq(ctx, -ctx->PARTA,sizeof(double));
        ctx->PARTA = 0;    
    }
    if (ctx->MARGA > 0) {
        free((char *)ctx->MARG);
        memrq(ctx, -ctx->MARGA,sizeof(double));
        ctx->MARGA = 0;    
    }
    if (ctx->DGFRA > 0) {
        free((char *)ctx->DGFR);
        memrq(ctx, -ctx->DGFRA,sizeof(int));
        ctx->DGFRA = 0;    
    }
    if (ctx->DFSA > 0) {
        free((char *)ctx->DFS);
        memrq(ctx, -ctx->DFSA,sizeof(int));
        ctx->DFSA = 0;    
    }
    if (ctx->LSISetA > 0) {
        free((char *)ctx->LSISet);
        memrq(ctx, -ctx->LSISetA,sizeof(int));
        ctx->LSISetA = 0;    
    }
    if (ctx->LSJSetA > 0) {
        free((char *)ctx->LSJSet);
        memrq(ctx, -ctx->LSJSetA,sizeof(int));
        ctx->LSJSetA = 0;    
    }
    if (ctx->LSdimA > 0) {
        free((char *)ctx->LSdim);
        memrq(ctx, -ctx->LSdimA,sizeof(short));
        ctx->LSdimA = 0;    
    }
    if (ctx->LSsizA > 0) {
        free((char *)ctx->LSsiz);
        memrq(ctx, -ctx->LSsizA,sizeof(short));
        ctx->LSsizA = 0;    
    }
    if (ctx->LScordA > 0) {
        free((char *)ctx->LScord);
        memrq(ctx, -ctx->LScordA,sizeof(short));
        ctx->LScordA = 0;    
    }
    if (ctx->LSconfA > 0) {
        free((char *)ctx->LSconf);
        memrq(ctx, -ctx->LSconfA,sizeof(short));
        ctx->LSconfA = 0;    
    }
    if (ctx->LSipA > 0) {
        free((char *)ctx->LSip);
        memrq(ctx, -ctx->LSipA,sizeof(short));
        ctx->LSipA = 0;    
    }
    if (ctx->LSimA > 0) {
        free((char *)ctx->LSim);
        memrq(ctx, -ctx->LSimA,sizeof(short));
        ctx->LSimA = 0;    
    }
    if (ctx->LSxA > 0) {
        free((char *)ctx->LSx);
        memrq(ctx, -ctx->LSxA,sizeof(float));
        ctx->LSxA = 0;    
    }
    if (ctx->LSyA > 0) {
        free((char *)ctx->LSy);
        memrq(ctx, -ctx->LSyA,sizeof(float));
        ctx->LSyA = 0;    
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  logscr_conf    Create actual model configuration.                       */

int logscr_conf(TDAContext *ctx, int m, int *iset, int *jset)
{
    register int i;
    int l,js,nm,nm1,np,np1,ilast,jlast;

    ilast = 1;
    np = nm = 0;
    
    while (1) {
        while (1) {
            ilast = logscr_comb(ctx, iset,ctx->NDIM,m,ilast);
            if (ilast)
                return(np);
            np++;
            np1 = (np - 1) * ctx->NDIM1;
            for (i = 1; i <= m; ++i)
                ctx->LSip[np1 + i] = (short)iset[i];
            if (m != 1)
                break;
        }
        jlast = 1;
        l = m - 1;

        while (1) {
            jlast = logscr_comb(ctx, jset,m,l,jlast);
            if (jlast)
                break;
            nm++;
            nm1 = (nm - 1) * ctx->NDIM1;
            for (i = 1; i <= l; ++i) {
                js = jset[i];
                ctx->LSim[nm1 + i] = (short)iset[js];
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

int logscr_comb(TDAContext *ctx, int *iset, int n, int m, int last) 
{
    (void)ctx;        /* unused: the signature is shared */
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

int logscr_eval(TDAContext *ctx, short *iar, int nc, int nv, int ibeg)
{
    register int i,j,jn;
    int k,kk,l,l1,df;

    df = 0;
    for (j = 1; j <= nc; ++j) {
        jn = (j - 1) * ctx->NDIM;
        kk = 1;
        l  = ibeg + j - 1;
        l1 = (l - 1) * ctx->NDIM1;

        for (i = 1; i <= nv; ++i) {
            k = iar[l1 + i];
            kk *= (ctx->LSdim[k] - 1);
            ctx->LSconf[jn + i] = (short)((int)k);
        }
        ctx->LSconf[jn + nv + 1] = 0;
        df += kk;
    }
    return(df);
}

/* ------------------------------------------------------------------------ */
/*  logscr_res  Initialize CTFit with mean values.                          */

void logscr_res(TDAContext *ctx)
{
    register int i;
    double avg;

    avg = (double)ctx->NOC / (double)ctx->NCTab;
    for (i = 1; i <= ctx->NCTab; ++i)  
        ctx->CTFit[i] = avg;
}

/* ------------------------------------------------------------------------ */
/*  logscr_lr  Calculate and return likelihood ratio Chi square statistic.  */

double logscr_lr(TDAContext *ctx)
{
    register int i;
    double tmp,gsq;

    gsq = 0.0;
    for (i = 1; i <= ctx->NCTab; ++i) {
        tmp = (double)ctx->CTFreq[i];
        if (ctx->CTFit[i] > ctx->EPSI && tmp > ctx->EPSI)  
            gsq += tmp * rlog(ctx, tmp / ctx->CTFit[i]);
    }
    gsq *= 2.0;
    return(gsq);
}

/* ------------------------------------------------------------------------ */
/*  logscr_fit                                                              */
/*      Fit a log-linear model. The algorithm is adapted from               */
/*      Haberman, Log-Linear Fit for Contingency Tables (AS 51),            */
/*      Applied Statistics 21 (1972), 218 - 225.                            */

void logscr_fit(TDAContext *ctx, int ncon)
{
    register int i,j,k,ii,kk;
    int iin,l,n,option,isz;
    double xmax,e;

    for (kk = 1; kk <= ctx->MxIter; ++kk) {
        xmax = 0.0;
        for (ii = 1; ii <= ncon; ++ii) {
            iin = (ii - 1) * ctx->NDIM;
            option = 1;
            ctx->LSsiz[1] = 1;
            n = 0;
            for (k = 1; k <= ctx->NDIM1; ++k) {
                if (!(l = (int)ctx->LSconf[iin + k]))
                    break;
                ctx->LSsiz[k + 1] = ctx->LSsiz[k] * (short)ctx->LSdim[l];
                n++;
            }
            isz = (int)ctx->LSsiz[n + 1];
            for (j = 1; j <= isz; ++j)  
                ctx->LSx[j] = ctx->LSy[j] = 0.0;
L_Cont1:
            for (k = 1; k <= ctx->NDIM; ++k)
                ctx->LScord[k] = 0;
            i = 1;
L_Cont2:
            j = 1;
            for (k = 1; k <= n; ++k) {
                l = (int)ctx->LSconf[iin + k];
                j += (int)ctx->LScord[l] * (int)ctx->LSsiz[k];
            }
            if (option) {
                ctx->LSx[j] += (float)ctx->CTFreq[i];
                ctx->LSy[j] += (float)ctx->CTFit[i];
            }
            else {
                if ((double)(ctx->LSy[j]) <= 0.0)
                    ctx->CTFit[i] = 0.0;
                else
                    ctx->CTFit[i] *= (double)(ctx->LSx[j] / ctx->LSy[j]);
            }
            i++;
            for (k = 1; k <= ctx->NDIM; ++k) {
                ctx->LScord[k] += 1;
                if ((int)ctx->LScord[k] < ctx->LSdim[k])
                    goto L_Cont2;
                ctx->LScord[k] = 0;
            }
            if (option) {
                option = 0;
                goto L_Cont1;
            }
            for (i = 1; i <= isz; ++i) {
                e = fabs((double)(ctx->LSx[i] - ctx->LSy[i]));
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

int e_expm(TDAContext *ctx)
{
    int err;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Expand model description. Current memory: %d bytes.\n\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,8,1))     /* get parameters */
        goto EEFin;

    printf1(ctx, "expm=%s\n",ctx->PMRHSTR);

    err = 0;
EEFin:
    p_clean(ctx);
    return(err);
}

