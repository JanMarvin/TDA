/****************************************************************************/
/*  t_nlreg                                                                 */
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
#include "t_eval.h"   
#include "t_eval3.h"
#include "t_alloc.h"   
#include "t_ml.h"   
#include "t_gf.h"   
#include "t_var.h"   
#include "t_lin.h"   
#include "t_cdf.h"   
#include "t_gmin.h"   
#include "t_sort.h"   
#include "tda_context.h"

/*  functions in t_nlreg.c */

int nlreg(TDAContext *ctx);
int  nlreg_cen(TDAContext *ctx, int iter);
void nlreg_coeff(TDAContext *ctx, int np);
void nlreg_pcov(TDAContext *ctx, int np);
void nlreg_res(TDAContext *ctx, int np,int nx);
void nlreg_prot(TDAContext *ctx, int np);
int doddrv(TDAContext *ctx, int n,int m,int np,double *x,char *ifixx,int ifixu, double *y,double *beta,char *ifixb, double *wd,int ldwd,double *w,int ndigit,double taufac, double sstol,double partol,double *work,int lwork, int *iwork,int liwork);
int dodmn(TDAContext *ctx, int n,int np,int m,double *x,char *ifixx,int ifixu,double *y, double *betac,char *ifixb,double *beta,double *betan,double *betas, double *s,double *delta,double *deltan,double *deltas,double *t, double *f,double *fn,double *fs,double *fjacb,double *fjacx, double *w,double *wd,int ldwd,double *ssf,double *ss, double *tt,int ldtt,double *xplusd,double *ddelt,double *sss, double *work,int lwork,int *iwork,int liwork);
void dodlm(TDAContext *ctx, int n,int np,int npp,int m,double *f,double *fjacb,             double *fjacx,double *w,double *wd,int ldwd,double *ss, double *tt,int ldtt,double *ddelt,double *alpha2,double *tau,double epsmac, double *sss,double *wrk1,double *tfjacb,double *omega,double *yt, double *u,double *qraux,double *wrk2,int *jpvt,double *s,double *t, int *nlms,double *rcond,int *irank);
void dodstp(TDAContext *ctx, int n,int np,int npp,int m,double *f,double *fjacb,            double *fjacx,double *w,double *wd,int ldwd,double *ss, double *tt,int ldtt,double *ddelt,double alpha,double epsmac, double *sss,double *tfjacb,double *vdtd,double *omega,double *yt, double *u,double *qraux,double *wrk2,int *jpvt,double *s,double *t, double *phi,int *irank,double *rcond);
int devfun(TDAContext *ctx, int n,int np,int m,double *betac,double *beta,char *ifixb, double *x,double *y,double *delta,double *xplusd,double *w,double *f);
int devjac(TDAContext *ctx, int n,int np,int npp,int m,double *betac,double *beta, char *ifixb,char *ifixx,int ifixu,double *x,double *delta, double *xplusd,double *ss,double *tt,int ldtt, int neta,double *pv,double *stp,double *fjacb, double *fjacx,double *w);
int djacfd(TDAContext *ctx, int n,int np,int m,double *beta,double *x,double *delta, double *xplusd,double *ss,double *tt,int ldtt,int neta, double *pv,double *stp,char *ifixb,double *fjacb, char *ifixx,int ifixu,double *fjacx);
void diniwk(TDAContext *ctx, int n,int m,int np,double *work,int *iwork,double *x,                          double *beta,double sstol,double partol,double taufac,int epsmai, int sstoli,int partli,int taufci,int ssfi,int tti,int ldtti);
int dpack(TDAContext *ctx, int n2,double *v1,double *v2,char *ifix);
void dunpac(TDAContext *ctx, int n2,double *v1,double *v2,char *ifix);
void ddiagi(TDAContext *ctx, int n,int m,double *s,int lds,double *v,int cdv,double *sv,int cdsv);
void ddiags(TDAContext *ctx, int n,int m,double *s,int lds,double *v,int cdv,double *sv,int cdsv);
void ddiagw(TDAContext *ctx, int n,int m,double *w,double *v,int cdv,double *wv,int cdwv);
void dwds(TDAContext *ctx, int n,int m,double *w,double *wd,int ldwd,double *t,double *wdt);
void dsetn(TDAContext *ctx, int n,int m,double *x,int *nrow);
int detaf(TDAContext *ctx, int n,int np,int m,double *xplusd,double *beta,double *eta, int *neta,double epsmac,int nrow,double *partmp,double *pvtemp);
void dchex(TDAContext *ctx, double *r,int cdr,int p,int k,int l,double *z,int cdz,int nz, double *c,double *s,int job);
void dsclb(TDAContext *ctx, int np,double *beta,double *ssf);
void dscld(TDAContext *ctx, int n,int m,double *x,double *tt);
void dwinf(TDAContext *ctx, int n,int m,int np,int *deltai,int *epsi,int *wssi,int *wssdei, int *wssepi,int *rvari,int *partli,int *sstoli,int *taufci,int *epsmai, int *olmavi,int *fjacbi,int *fjacxi,int *xplusi,int *betaci,int *betasi, int *betani,int *deltsi,int *deltni,int *ddelti,int *fsi,int *fni,int *si, int *sssi,int *ssi,int *ssfi,int *ti,int *tti,int *taui,int *alphai, int *vcvi,int *omegai,int *yti,int *ui,int *qrauxi,int *wrk1i,int *sei, int *rcondi,int *etai,int *actrsi,int *pnormi,int *prersi,int *rnorsi, int *lwkmn);
void diwinf(TDAContext *ctx, int m,int np,int *jpvti,int *nnzwi,int *nppi,int *idfi, int *nrowi,int *netai,int *int2i,int *iranki,int *ldtti,int *liwkmn);
void dacces(TDAContext *ctx, int n,int m,int np,double *work,int lwork,int *iwork,int liwork, int access,int *jpvt,int *wrk1,int *tfjacb,int *omega,int *yt,int *u,int *qraux, int *wrk2,int *nnzw,int *npp,double *partol,double *sstol, double *taufac,double *epsmac,int *neta, double *wss,double *wssdel,double *wsseps,double *rvar,int *idf, double *tau,double *alpha,int *int2,double *olmavg,double *rcond, int *irank,double *actrs,double *pnorm,double *prers,double *rnorms);
void didts(TDAContext *ctx, int n,int m,double *w,double *wd,int ldwd,double alpha,double *tt, int ldtt,double *t,double *dtt);
int odr_fun(TDAContext *ctx, int n,int np,int m,double *beta,double *xplusd,double *f,int cdf); 
int odr_jac(TDAContext *ctx, int n,int np,int m,double *beta,double *xplusd,            double *fjacb,double *fjacx);

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */


/* ------------------------------------------------------------------------ */
/*  nlreg() Nonlinear regression                                            */  
/*                                                                          */
/*          nlreg(                                                          */
/*              v=...,      varlist: Y,X1,X2,...                            */
/*              ni=...,     1 if without intercept                          */
/*              yw=...,     censoring indicator                             */
/*              opt=...,    1 = OLS, 2 = ODR, def. 1                        */
/*              mxit=...,   max number of iterations, def. 50               */
/*              mxcyc=...,  max number of cycles, def. 20                   */
/*              tolsg=...,  tolerance for squared convergence, def. 1e-8    */
/*              tolsp=...,  tolerance for parameter convergence, def. 1e-10 */
/*              tolp=...,   tolerance for cycles (censored data)            */
/*              xp=...,     starting values                                 */
/*              dsv=...,    starting values                                 */
/*              tfmt=...,   print format for parameters, def. 10.4          */
/*              ppar=...,   print parameter                                 */
/*              pcov=...,   print covariance matrix                         */
/*              mfmt=...,   print format pcov                               */
/*              pres=...,   print residuals                                 */
/*              fmt=...,    print format pres                               */
/*              prot=...,   protocol file, if prot1=... print also          */
/*                          delta scale factors.                            */
/*              pfmt=...,   print format for protocol file, def. -19.11     */
/*                                                                          */
/*          ) = function;       optional                                    */
/*                                                                          */
/*          Return 0 if OK, otherwise -1.                                   */

int nlreg(TDAContext *ctx)
{
    register int i,j;
    int err,nv,np,nx,lwork,liwork,ifixu,ldwd,ndigit,ncen,conv;
    double tmp,wt,wsum,ww[2],wwd[2],*wptr;
    char *ifixx;

    ww[1] = wwd[1] = -1.0; ldwd = 1; ifixx = NULL;          

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Nonlinear regression. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLSG = 1.e-8;              /* default tolerances */
    ctx->TOLSP = 1.e-10;

    if (parm(ctx, ctx->CmdBuf + 5,9,0))   /* get parameters */
        goto NLREGFin;

    if (ctx->PMNW != 1) {
        printf1(ctx, "Error: nlreg command requires cross-section data (nw=1).\n");
        goto NLREGFin;
    }
    nv = check_nvar(ctx, 0);         /* check variables */
    if (nv < 2) {               /* PMNV is number of variables */       
        printf1(ctx, "Error: need at least two variables.\n");
        goto NLREGFin;        
    }
    nx = nv - 1;                /* number of indep. variables */

    if (ctx->PMMXCYC < 0)
        ctx->PMMXCYC = 20;

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 10,4);

    ndigit = 0;
    if (ctx->PMNDIGFlg)
        ndigit = (int)ctx->PMNDIGIT;

    printf1(ctx, "\nType: ");
    if (ctx->PMOPT == 1) {
        printf1(ctx, "ordinary least squares regression.\n");
        ctx->ODRTyp = 0;
    }
    else {
        printf1(ctx, "orthogonal distance regression.\n");
        ctx->ODRTyp = 1;
    }
    printf1(ctx, "Function: ");
    if (ctx->FNFlg) {
        printf1(ctx, "user-defined.\n");
        np = ctx->FNArgN;
        if (ctx->ODRTyp)
            ctx->PMDOPT = 1;     /* always numerical approximation */
    }
    else {
        printf1(ctx, "linear predictor [");
        np = ctx->PMNV - 1;
        if (ctx->PMNI) {
            printf1(ctx, "without");
            ctx->PMNI = 1;
        }
        else {
            printf1(ctx, "including");
            np++;
        }
        printf1(ctx, " intercept].\n");
    }
    printf1(ctx, "Derivatives: ");
    if (ctx->PMDOPT == 1) {
        ctx->ODRAd = 0;
        printf1(ctx, "numerical approximation.\n\n");
    }
    else {
        ctx->ODRAd = 1;
        printf1(ctx, "analytical.\n\n");
    }

    prn_nwvar(ctx, 0);                       /* print variables */
    if (np < 1) {
        printf1(ctx, "Error: no parameters.\n");
        goto NLREGFin;
    }
    ctx->ODRMXIt = 50;
    if (ctx->MxItFlg)
        ctx->ODRMXIt = ctx->MxIter;

    printf1(ctx, "Number of parameters:            %d\n",np);
    printf1(ctx, "Number of cases:                 %d\n",ctx->NOC);
    printf1(ctx, "Max number of iterations:        %d\n",ctx->ODRMXIt);
    printf1(ctx, "Tolerance squared convergence:   %g\n",ctx->TOLSG);
    printf1(ctx, "Tolerance parameter convergence: %g\n",ctx->TOLSP);
    newline(ctx);  

    if (ctx->FNFlg && ctx->ODRAd) {        /* allocate memory for derivatives */

        if (fnd_alloc(ctx, 1,1,ctx->FNArgN,ctx->FNPN,0))   
            goto NLREGFin;
    }

    /* get data into AcY[] and AcX[], optionally weights in AcW[] */
    /* AcU[] is for model parameters */                       

    if (alloc_acy(ctx, ctx->NOC + 1))
        goto NLREGFin;
    if (alloc_acx(ctx, ctx->NOC * nx + 1))
        goto NLREGFin;
    if (alloc_acu(ctx, np + 1))
        goto NLREGFin;

    if (ctx->WIVar >= 0) {
        if (alloc_acw(ctx, ctx->NOC + 1))
            goto NLREGFin;
    }

    /* if yw=... is used to define a censoring indicator save this
       information in AcNS[], also allocate: AcY1 for previous parameter
       values, AcZ for residuals, AcM for sorting pointer, AcT temporary. */

    if (ctx->PMYWVar >= 0) {
        if (alloc_acns(ctx, ctx->NOC + 1))
            goto NLREGFin;
        if (alloc_acz(ctx, ctx->NOC + 1))
            goto NLREGFin;
        if (alloc_act(ctx, ctx->NOC + 1))
            goto NLREGFin;
        if (alloc_acm(ctx, ctx->NOC + 1))
            goto NLREGFin;
        if (alloc_acy1(ctx, np + 1))
            goto NLREGFin;
    }
    ncen = 0;
    wt = 1.0;
    wsum = 0.0;
    for (i = 0; i < ctx->NOC; ++i) {
        ctx->AcY[i + 1] = get_data(ctx, ctx->PMVIdx[0],i);
        for (j = 1; j <= nx; ++j)  
            ctx->AcX[i * nx + j] = get_data(ctx, ctx->PMVIdx[j],i);
             
        if (ctx->WIVar >= 0) {
            wt = get_data(ctx, ctx->WIVar,i);
            if (wt < 0.0) {
                printf1(ctx, "Error: found negative weight in case %d.\n",i + 1);
                goto NLREGFin; 
            }
            wsum += wt;
            ctx->AcW[i + 1] = wt;
        }
        if (ctx->PMYWVar >= 0) {
            if ((get_data(ctx, ctx->PMYWVar,i)) != 0.0)
                ctx->AcNS[i] = 1;
            else {
                ctx->AcNS[i] = 0;
                ncen++;
            }
        }
    }
    if (ctx->WIVar >= 0) {
        printf1(ctx, "Using weights defined by: %s\n",ctx->VName[ctx->WIVar]);
        if (wsum < ctx->EPSI1) {
            printf1(ctx, "Error: sum of weights is almost zero.\n");
            goto NLREGFin;
        }
        wptr = ctx->AcW;
    }
    else {
        ww[1] = -1.0;       /* no weights */
        wptr = ww;
    }
    if (get_dsv(ctx, np,ctx->AcU,0,ctx->ParLB,ctx->ParUB,1))    /* try to get starting values */
        goto NLREGFin;                      /* ParLB,ParUB not used */

    if (ctx->FNFlg) {            /* show parameters and starting values */

        if (prn_fsval(ctx, np,ctx->AcU + 1,0,ctx->ParLB + 1,ctx->ParUB + 1,0))
            goto NLREGFin;
    }
    nlreg_prot(ctx, np);         /* init protocol file */

    /* we don't use ifixx and ifixb */

    ifixu = 0;
    if (alloc_acc(ctx, np + 1))      /* used for ifixb[] */
        goto NLREGFin;

    for (i = 1; i <= np; ++i)   /* all parameters free */
        ctx->AcC[i] = 1;

    /* allocate working memory */

    lwork = 17 + 7 * ctx->NOC + 10 * ctx->NOC * nx + 2 * ctx->NOC * np + 8 * np + 2;
    liwork = 19 + 2 * np + nx;

    if (alloc_acv(ctx, lwork + 1))
        goto NLREGFin;

    if (alloc_acn(ctx, liwork + 1))
        goto NLREGFin;

    if (ctx->PMYWVar >= 0) {
        printf1(ctx, "Number of censored observations: %d [%g pct]\n",ncen,
                                            100.0 * (double)ncen / (double)ctx->NOC);
        printf1(ctx, "Will perform up to %d cycles to adjust for censored cases.\n",ctx->PMMXCYC);
        printf1(ctx, "Tolerance for y-convergence: %g\n\n",ctx->TOLP);
    }

    /* call driver */

    conv = 0;
    for (ctx->ORDCyc = 0; ctx->ORDCyc <= ctx->PMMXCYC; ++ctx->ORDCyc) {

        ctx->ODRCOVPtr = ctx->ODRSEPtr = ctx->ODRDf = ctx->ODRNFJ = ctx->ODRNF = ctx->ODRNIt  = 0;   
        ctx->ODRXPDPtr = ctx->ODRFNPtr = ctx->ODREPSPtr = 0;
        ctx->ODRWEps = ctx->ODRWDel = ctx->ODRRVar = 0.0; 

        err = doddrv(ctx, ctx->NOC,nx,np,ctx->AcX,ifixx,ifixu,ctx->AcY,ctx->AcU,ctx->AcC,wwd,ldwd,wptr,ndigit,
              0.0,ctx->TOLSG,ctx->TOLSP,ctx->AcV,lwork,ctx->AcN,liwork);

        if (ctx->PMYWVar < 0) {
            conv = 1;
            break;
        }
        printf1(ctx, "Cycle%3d : ",ctx->ORDCyc);
        if (err <= 0 || ctx->ODRFNPtr == 0) {
            printf1(ctx, "cannot find solution.\n");
            break;
        }
/* ## */
        if (ctx->ORDCyc > 0) {       /* check convergence */
            wsum = 0.0;
            for (j = 1; j <= np; ++j) {
                tmp = fabs(ctx->AcU[j] - ctx->AcY1[j]) / dmax(ctx, 1.0,ctx->AcU[j]);
                wsum = dmax(ctx, wsum,tmp);
            }
            printf1(ctx, "%17.10f",wsum);
        }
        else
            printf1(ctx, "     Tolerance   ");
        printf1(ctx, "  Iterations: %d  Function calls: %2d (%d)\n",ctx->ODRNIt,ctx->ODRNF,ctx->ODRNFJ);

        for (j = 1; j <= np; ++j)   /* save parameters */
            ctx->AcY1[j] = ctx->AcU[j];

        if (ctx->ORDCyc > 0 && wsum <= ctx->TOLP) {
            conv = 1;
            break;
        }
        if (nlreg_cen(ctx, ctx->ORDCyc))      /* update censored values */
            break;
    }
     
    printf1(ctx, "\nConvergence ");
    if (err <= 0 || conv == 0)
        printf1(ctx, "not ");
    printf1(ctx, "reached. Number of iterations: %d. Function calls: %d (%d)\n",
                                                        ctx->ODRNIt,ctx->ODRNF,ctx->ODRNFJ);
    if (err <= 0) {
        printf1(ctx, "Problem: ");
        switch (err) {
            case  0:    printf1(ctx, "exceeded max number of iterations.\n");
                        break;
            case -5:    printf1(ctx, "found negative weights.\n");        
                        break;
            case -6:    printf1(ctx, "cannot evaluate function.\n");
                        break;
            case -7:    printf1(ctx, "cannot determine function accuracy.\n");
                        break;
            case -8:    printf1(ctx, "step size search failed.\n");
                        break;
            case -9:    printf1(ctx, "cannot evaluate function or derivatives.\n");
                        break;
            default:    printf1(ctx, "not specified (%d).\n",err);
                        break;
        }
        if (err < 0)
            goto NLREGFin;
    }
    else {
        printf1(ctx, "Type: ");
        if (err == 1)
            printf1(ctx, "squared ");
        else if (err == 2)
            printf1(ctx, "parameter ");
        else                 
            printf1(ctx, "squared and parameter ");
        printf1(ctx, "convergence.\n");
    }
    if (err > 0) {
        printf1(ctx, "\nRank deficiency: %d\n",ctx->ODRRDef);
        if (ctx->ODRRDef == 0) {
            printf1(ctx, "Degrees of freedom: %d\n",ctx->ODRDf);
            printf1(ctx, "Sum of squared epsilon: %g\n",ctx->ODRWEps);
            if (ctx->ODRTyp)
                printf1(ctx, "Sum of squared delta: %g\n",ctx->ODRWDel);
            printf1(ctx, "Variance of residuals: %g\n",ctx->ODRRVar);
        }
    }
         
    nlreg_coeff(ctx, np);

    if (ctx->PMCovFDef && ctx->ODRCOVPtr)     /* write cov matrix */
        nlreg_pcov(ctx, np);

    if (ctx->PMResFDef)                  /* write residuals */                 
        nlreg_res(ctx, np,nx);

    /**** pmat("work",lwork,1,AcV); ***/
             
    err = 0;
    
NLREGFin:
    if (ctx->FNFlg)
        fnd_alloc(ctx, 0,0,0,0,0);
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  nlreg_cen()     update dep. variable for censored cases.                */
/*                  return 0 if ok., -1 if error                            */
/*                                                                          */
/*  Semantics:                                                                 */
/*  this is the Buckley-James (1979) estimator.  Residuals are sorted;      */
/*  each censored case's probability mass is redistributed equally to       */
/*  the residuals to its right (Efron's redistribute-to-the-right =         */
/*  the Kaplan-Meier construction), and the censored response is            */
/*  replaced by fit + the KM-weighted mean of the larger uncensored         */
/*  residuals.  The mxcyc= cycles are the standard Buckley-James            */
/*  iteration.  A line-faithful R replication reproduces TDA's              */
/*  estimates to ~4 decimals; pinned in test-r-comparisons.R.               */

int nlreg_cen(TDAContext *ctx, int iter)
{
    (void)iter;        /* unused: the signature is shared */
    register int i,j,k,l;
    int m;
    double tmp,tmp1;
             
    for (i = 0; i < ctx->NOC; ++i) {
        tmp = ctx->AcV[ctx->ODRFNPtr + i + 1];
        ctx->AcZ[i] = get_data(ctx, ctx->PMVIdx[0],i) - tmp;
        if (ctx->AcNS[i] == 0)
            ctx->AcY[i + 1] = tmp;
    }
    if (sortdp2a(ctx, ctx->NOC,ctx->AcZ,ctx->AcNS,ctx->AcM,0))      /* sort */
        return(-1);

    for (i = 0; i < ctx->NOC; ++i)  
        ctx->AcT[i] = 0.0;

    tmp = 1.0 / (double)ctx->NOC;
    m = ctx->NOC - 1;
    for (i = 0; i < ctx->NOC; ++i) {
        j = ctx->AcM[i];
        ctx->AcT[j] += tmp;
        if (ctx->AcNS[j] == 0 && i < ctx->NOC - 1) {
            tmp1 = ctx->AcT[j] / (double)m;
            for (k = i + 1; k < ctx->NOC; ++k)  
                ctx->AcT[ctx->AcM[k]] += tmp1;
            ctx->AcT[j] = 0.0;
        }
        m--;
    }
    tmp = 0.0;
    for (i = 0; i < ctx->NOC; ++i) {
        k = ctx->AcM[i];
        tmp += ctx->AcT[k];
    }
    for (i = 0; i < ctx->NOC; ++i) {
        j = ctx->AcM[i];
        if (ctx->AcNS[j] == 0) {         /* if censored */
            if (i == ctx->NOC - 1) {
                ctx->AcY[j + 1] += ctx->AcZ[j];
            }
            else {
                tmp = tmp1 = 0.0;
                for (k = i + 1; k < ctx->NOC; ++k) {
                    l = ctx->AcM[k];
                    if (ctx->AcNS[l] != 0 || k == ctx->NOC - 1) {
                        tmp += ctx->AcZ[l] * ctx->AcT[l];
                        tmp1 += ctx->AcT[l];
                    }
                }
                if (tmp1 != 0.0)
                    ctx->AcY[j + 1] += tmp / tmp1;
            }
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  nlreg_coeff()   print estimated coefficients.                           */

void nlreg_coeff(TDAContext *ctx, int np)
{
    register int i,k,l; 
    double tmp,tmp1;

    if (ctx->FNFlg == 0) {
        prn1_coeff(ctx, np,ctx->AcU,ctx->AcV + ctx->ODRSEPtr,ctx->PMNI,ctx->ODRDf,ctx->PMVIdx,0);
    }
    else {
        l = 9;
        if (l < ctx->FNMLen)
            l = ctx->FNMLen;

        printf1(ctx, "\nIdx  Parameter ");
        prnchar(ctx, ' ',l - 9,0);
        prnchar(ctx, ' ',ctx->PMTFmt1 - 4,0); printf1(ctx, "Value");
        prnchar(ctx, ' ',ctx->PMTFmt1 - 4,0); printf1(ctx, "Error");
        prnchar(ctx, ' ',ctx->PMTFmt1 - 6,0); printf1(ctx, "Value/E  Signif\n");
        prnchar(ctx, '-',14 + l + 3 * (ctx->PMTFmt1 + 1),1);
    
        for (i = 1; i <= ctx->FNArgN; ++i) {
            k = ctx->FNArgSP[i - 1];
            printf1(ctx, "%3d  %s  ",i,ctx->FNArgDef[k]);
            prnchar(ctx, ' ',l - (int)strlen(ctx->FNArgDef[k]),0);
    
            rt_printf1_d(ctx, ctx->PMTFmtS,ctx->AcU[i]); 

            if (ctx->ODRSEPtr > 0 && ctx->ODRDf > 0) {
                tmp = ctx->AcV[ctx->ODRSEPtr + i];
                rt_printf1_d(ctx, ctx->PMTFmtS,tmp);
                if (tmp > 0.0) {
                    tmp1 = ctx->AcU[i] / tmp;
                    tmp = 2.0 * cdtf(ctx, fabs(tmp1),ctx->ODRDf) - 1.0;
                    rt_printf1_d(ctx, ctx->PMTFmtS,tmp1);
                    printf1(ctx, "%7.4f",tmp);
                }
                else {
                    prnchar(ctx, ' ',ctx->PMTFmt1 - 3,0);
                    printf1(ctx, "---     ---"); 
                }
            }
            else {
                prnchar(ctx, ' ',ctx->PMTFmt1 - 3,0); printf1(ctx, "--- "); 
                prnchar(ctx, ' ',ctx->PMTFmt1 - 3,0); printf1(ctx, "---     ---"); 
            }
            newline(ctx);
        }
    }
    newline(ctx);           
}

/* ------------------------------------------------------------------------ */
/*  nlreg_pcov()    print covariance matrix.                                */

void nlreg_pcov(TDAContext *ctx, int np)
{
    register int i,j;
    double tmp,*cov;

    cov = ctx->AcV + ctx->ODRCOVPtr;
    for (i = 1; i <= np; ++i) {
        for (j = 1; j <= np; ++j) {
            if (j >= i)
                tmp = cov[(i - 1) * np + j];
            else
                tmp = cov[(j - 1) * np + i];

            rt_fprintf_d(ctx, ctx->PMCovFd,ctx->PMMFmtS,tmp);
        }
        fprintf(ctx->PMCovFd,"\n");
    }
    p_wmsg(ctx, 2,ctx->PMCovFName,1);
}

/* ------------------------------------------------------------------------ */
/*  nlreg_res()     print data and residuals.                               */

void nlreg_res(TDAContext *ctx, int np,int nx)
{
    (void)np;        /* unused: the signature is shared */
    register int i,j,k;
    double tmp;

    for (i = 1; i <= ctx->NOC; ++i) {
        fprintf(ctx->PMResFd,"%6d ",i);
        tmp = get_data(ctx, ctx->PMVIdx[0],i - 1);
        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,tmp);
        tmp = -9.999;
        if (ctx->ODRFNPtr)
            tmp = ctx->AcV[ctx->ODRFNPtr + i];
        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,tmp);
        tmp = -9.999;
        if (ctx->ODREPSPtr)
            tmp = ctx->AcV[ctx->ODREPSPtr + i];
        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,tmp);

        tmp = -9.999;
        for (j = 1; j <= nx; ++j) {
            rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,ctx->AcX[(i-1)*nx+j]);
            if (ctx->ODRTyp) {
                if (ctx->ODRXPDPtr)
                    tmp = ctx->AcV[ctx->ODRXPDPtr + (i-1)*nx + j];
                rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,tmp);
            }
        }
        if (ctx->WIVar >= 0)  
            rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,ctx->AcW[i]);

        fprintf(ctx->PMResFd,"\n");
    }
    p_wmsg(ctx, 3,ctx->PMResFName,1);

    if (ctx->PMTDAFDef) {
        fprintf(ctx->PMTDAFd,"# data written by nlreg command.\n");
        fprintf(ctx->PMTDAFd,"nvar(\n");
        fprintf(ctx->PMTDAFd,"  dfile = %s,\n",ctx->PMResFName);
        fprintf(ctx->PMTDAFd,"  noc = %d,\n",ctx->NOC);
        k = 1;
        fprintf(ctx->PMTDAFd,"  Case");
        fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - 4,0);
        fprintf(ctx->PMTDAFd,"[6.0] = c%-2d,\n",k++);
        j = ctx->PMVIdx[0];
        fprintf(ctx->PMTDAFd,"  %s",ctx->VName[j]);
        fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - (int)strlen(ctx->VName[j]),0);
        fprintf(ctx->PMTDAFd,"[%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,k++);

        fprintf(ctx->PMTDAFd,"  E%s",ctx->VName[j]);
        fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - (int)strlen(ctx->VName[j]) - 1,0);
        fprintf(ctx->PMTDAFd,"[%d.%d] = c%-2d, # estimated\n",ctx->PMFmt1,ctx->PMFmt2,k++);

        fprintf(ctx->PMTDAFd,"  EPS");
        fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - 3,0);
        fprintf(ctx->PMTDAFd,"[%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,k++);

        for (i = 1; i <= nx; ++i) {
            j = ctx->PMVIdx[i];
            fprintf(ctx->PMTDAFd,"  %s",ctx->VName[j]);
            fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - (int)strlen(ctx->VName[j]),0);
            fprintf(ctx->PMTDAFd,"[%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,k++);

            if (ctx->ODRTyp) {
                fprintf(ctx->PMTDAFd,"  E%s",ctx->VName[j]);
                fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - (int)strlen(ctx->VName[j]) - 1,0);
                fprintf(ctx->PMTDAFd,"[%d.%d] = c%-2d, # estimated\n",ctx->PMFmt1,ctx->PMFmt2,k++);
            }
        }
        if (ctx->WIVar >= 0) {
            fprintf(ctx->PMTDAFd,"  %s",ctx->VName[ctx->WIVar]);
            fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - (int)strlen(ctx->VName[ctx->WIVar]),0);
            fprintf(ctx->PMTDAFd,"[%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,k++);
        }
        fprintf(ctx->PMTDAFd,");\n");
    
        printf1(ctx, "TDA description written to: %s\n",ctx->PMTDAFName);
    }
}

/* ------------------------------------------------------------------------ */
/*  nlreg_prot()  protocol file.                                            */

void nlreg_prot(TDAContext *ctx, int np)
{
    register int i;

    if (ctx->PMProtFDef == 0)
        return;

    printf1(ctx, "Protocol will be written to: %s\n",ctx->PMProtFName);

    fprintf(ctx->PMProtFd,"Protocol: nonlinear regression.\n");
    if (ctx->ODRTyp == 0)
        fprintf(ctx->PMProtFd,"Ordinary least squares regression.\n");
    else  
        fprintf(ctx->PMProtFd,"Orthogonal distance regression.\n");

    fprintf(ctx->PMProtFd,"Number of model parameters: %d\n",np);
    fprintf(ctx->PMProtFd,"\nStarting values.\n");
    for (i = 1; i <= np; ++i)
        rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,ctx->AcU[i]);
    fprintf(ctx->PMProtFd,"\n\n");
}

/*  ----------------------------------------------------------------------- */    
/*  int doddrv()                                                            */
/*                                                                          */
/*  Program for orthogonal distance regression, adapted from the FORTRAN    */  
/*  source code in ACM algorithm 676. The original program was written by   */
/*  Paul T. Boggs, Richard H. Byrd, Janet R. Donaldson, Robert B. Schnabel. */
/*                                                                          */
/*  REferences:                                                             */
/*  BOGGS, P. T., R. H. BYRD, J. R. DONALDSON, AND R. B. SCHNABEL (1987),   */
/*  "ODRPACK -- SOFTWARE FOR WEIGHTED ORTHOGONAL DISTANCE REGRESSION,"      */
/*  UNIVERSITY OF COLORADO DEPARTMENT OF COMPUTER SCIENCE                   */
/*  TECHNICAL REPORT NUMBER CU-CS-360-87.                                   */
/*  (TO APPEAR IN ACM TRANS. MATH. SOFTWARE.)                               */
/*                                                                          */
/*  BOGGS, P. T., R. H. BYRD, J. R. DONALDSON, AND R. B. SCHNABEL (1989),   */
/*  "REFERENCE GUIDE FOR ODRPACK SOFTWARE FOR WEIGHTED                      */
/*  ORTHOGONAL DISTANCE REGRESSION,"                                        */
/*  ONLINE DOCUMENTATION AVAILABLE FROM AUTHORS                             */
/*                                                                          */
/*  BOGGS, P. T., R. H. BYRD, AND R. B. SCHNABEL (1987),                    */
/*  "A STABLE AND EFFICIENT ALGORITHM FOR NONLINEAR                         */
/*  ORTHOGONAL DISTANCE REGRESSION,"                                        */
/*  SIAM J. SCI. STAT. COMPUT., 8(6):1052-1078.                             */
/*                                                                          */
/*  doddrv(n,m,np,x,*ifixx,ifixu,*y,*beta,*ifixb,*wd,ldwd,*w,ndigit,        */
/*      taufac,sstol,partol,*work,lwork,*iwork,liwork)                      */
/*                                                                          */
/*  Return   1 squared convergence                                          */
/*           2 parameter convergence                                        */
/*           3 both squared and parameter convergence                       */
/*           0 exceeded max number of iterations                            */
/*          -1 error in input parameters                                    */
/*          -2 if number of parameters zero or more than # of cases         */
/*          -3 if lwork insufficient                                        */
/*          -4 if liwork insufficient                                       */
/*          -5 if negative weights                                          */
/*          -6 cannot evaluate function with starting values                */
/*          -7 cannot calculate number of good digits                       */
/*          -8 error in step size search                                    */
/*          -9 cannot evaluate function or derivatives                      */

int doddrv(TDAContext *ctx, int n,int m,int np,double *x,char *ifixx,int ifixu,double *y, double *beta,char *ifixb,double *wd,int ldwd,double *w,int ndigit, double taufac,double sstol,double partol,double *work,int lwork, int *iwork,int liwork)
{
    register int i;
    double epsmac,eta,tmp;      
    int err,actrsi,alphai,betaci,betani,betasi,ddelti,deltai,deltni,deltsi;
    int nn,npp,epsmai,etai,fi,fjacbi,fjacxi,fni,fsi,idfi,int2i;
    int iranki,jpvti,ldtt,ldtti,liwkmn,lwkmn,neta,netai;
    int nnzwi,nppi,nrow,nrowi,olmavi,omegai;
    int partli,pnormi,prersi,qrauxi,rcondi,rnorsi,rvari,si,ssfi,ssi;
    int sssi,sstoli,taufci,taui,tfjaci,ti,tti,ui,wrk1i,wrk2i,wssi;
    int wssdei,wssepi,xplusi,yti;
 
    err = 0; 
    if (n < 1 || m < 1 || np < 1)
        return(-1);

    /* SET STARTING LOCATIONS WITHIN INTEGER WORKSPACE */

    diwinf(ctx, m,np,&jpvti,&nnzwi,&nppi,&idfi,&nrowi,&netai,         
        &int2i,&iranki,&ldtti,&liwkmn);
 
    /* SET STARTING LOCATIONS WITHIN DOUBLE PRECISION WORK SPACE */
   
    dwinf(ctx, n,m,np,&deltai,&fi,&wssi,&wssdei,&wssepi,&rvari,
         &partli,&sstoli,&taufci,&epsmai,&olmavi,
         &fjacbi,&fjacxi,&xplusi,&betaci,&betasi,&betani,&deltsi,
         &deltni,&ddelti,&fsi,&fni,&si,&sssi,&ssi,&ssfi,&ti,&tti,&taui,
         &alphai,&tfjaci,&omegai,&yti,&ui,&qrauxi,&wrk1i,&wrk2i,&rcondi,
         &etai,&actrsi,&pnormi,&prersi,&rnorsi,&lwkmn);

    /* check input parameter */
      
    if (ldwd != 1 && ldwd != n)
        return(-1);

    npp = 0;                     /* # of parameters to be estimated */
    for (i = 1; i <= np; ++i) {
        if (ifixb[i])      
            npp += 1;     
    }
    if (n <= 0 || m <= 0 || npp <= 0 || npp > n)
        return(-2);
 
    if (lwork < lwkmn)
        return(-3);   
    if (liwork < liwkmn)      
        return(-4);
 
    if (w[1] >= 0.0) {    
         nn = 0;
         for (i = 1; i <= n; ++i) {
            if (w[i] < 0.0)       
                return(-5);
            else if (w[i] > 0.0)       
               nn++;    
        }
        if (nn < npp)      
            return(-1);
    }         
    if (wd[1] >= 0.0) {    
        nn = ldwd * m;
        for (i = 1; i <= nn; ++i) {
            if (wd[i] <= 0.0)
                return(-6);
        }
    }
    for (i = 0; i <= liwork; ++i)       /* initialization */
        iwork[i] = 0;

    for (i = 0; i <= lwork; ++i)
        work[i] = 0.0;    
    
    diniwk(ctx, n,m,np,work,iwork,x,beta,sstol,partol,taufac,
           epsmai,sstoli,partli,taufci,ssfi,tti,ldtti);

    work[taui] = -work[taufci];

    /* SET UP FOR PARAMETER ESTIMATION -
       PULL BETA'S TO BE ESTIMATED AND CORRESPONDING SCALE VALUES
       AND STORE IN WORK(betaci) AND WORK(ssi), RESPECTIVELY */
 
    iwork[nppi] = dpack(ctx, np,work+betaci-1,beta,ifixb);

    if (work[ssfi] > 0.0)        
        iwork[nppi] = dpack(ctx, np,work+ssi-1,work+ssfi-1,ifixb);
    else  
        work[ssi] = work[ssfi];

    /* EVALUATE THE WEIGHTED EPSILONS AT THE STARTING POINT */
   
    err = devfun(ctx, n,np,m,work+betaci-1,beta,ifixb,x,y,
              work+deltai-1,work+xplusi-1,w,work+fi-1);
    if (err)      
        return(-6);  

    /* FIND NUMBER OF NONZERO WEIGHTS */
 
    if (w[1] < 0.0)       
        iwork[nnzwi] = n;
    else {
        iwork[nnzwi] = 0;
        for (i = 1; i <= n; ++i) {
            if (w[i] > 0.0)       
                iwork[nnzwi] += 1;                
        }
    }        
 
    /* COMPUTE NORM OF THE INITIAL ESTIMATES */
  
    ddiags(ctx, iwork[nppi],1,work+ssi-1,iwork[nppi],work+betaci-1,1,work+sssi-1,1);           
    ddiags(ctx, n,m,work+tti-1,iwork[ldtti],work+deltai-1,m,work+sssi+iwork[nppi]-1,m);
    work[pnormi] = dnrm2(ctx, iwork[nppi]+n*m,work+sssi-1,1);
    
    /* COMPUTE SUM OF SQUARES OF THE WEIGHTED EPSILONS AND WEIGHTED DELTAS */
 
    dcopy(ctx, n,work+fi-1,1,work+sssi-1,1);
    work[wssepi] = ddot(ctx, n,work+sssi-1,1,work+sssi-1,1);
    dwds(ctx, n,m,w,wd,ldwd,work+deltai-1,work+sssi+n-1);
    work[wssdei] = ddot(ctx, n*m,work+sssi+n-1,1,work+sssi+n-1,1);
    work[wssi] = work[wssepi] + work[wssdei];



    if (ctx->PMProtFDef) {      /* print some information to protocol file */

        prval(ctx, "Factor for trust region radius",work[taufci]);
        prval(ctx, "Initial weighted sum of epsilon",work[wssepi]);
        prvec(ctx, "Parameter scale factors",np,work+ssi-1);

        if (ctx->PMProtFDef > 1 && ctx->ODRTyp) {
            prval(ctx, "Initial weighted sum of delta",work[wssdei]);
            prmat(ctx, "Delta scale factors",n,m,work+tti-1);
        }
    }
   
    /* SELECT FIRST ROW OF X + DELTA THAT CONTAINS NO ZEROS */

    nrow = -1;
    dsetn(ctx, n,m,work+xplusi-1,&nrow); 
    iwork[nrowi] = nrow;

    /* SET NUMBER OF GOOD DIGITS IN FUNCTION RESULTS */
 
    epsmac = work[epsmai];
    tmp = -rlog(ctx, epsmac) / log(10.0);

    if ((ndigit < 2) || (ndigit > (int)tmp)) {     
        iwork[netai] = -1;
        err = detaf(ctx, n,np,m,work+xplusi-1,beta,&eta,&neta,epsmac,
                                   nrow,work+betani-1,work+fni-1);
        if (err) {    
            if (ctx->PMProtFDef) 
                fprintf(ctx->PMProtFd,"\nCannot evaluate function accuracy.\n");
            return(-7);
        }
        else {
           iwork[netai] = neta;
           work[etai] = eta;
        }      
    }
    else {
        iwork[netai] = ndigit;
        tmp = (double)(-ndigit);
        work[etai] = pow(10.0,tmp);
    }         
    if (ctx->PMProtFDef) {
        fprintf(ctx->PMProtFd,"\nNumber of good digits: %d\n",iwork[netai]);
        fprintf(ctx->PMProtFd,"Accuracy: %g\n",work[etai]);
    }        
    ldtt = iwork[ldtti]; 
    err = dodmn(ctx, n,np,m,x,ifixx,ifixu,y,work+betaci-1,ifixb,beta,work+betani-1,
          work+betasi-1,work+si-1,work+deltai-1,work+deltni-1,work+deltsi-1,
          work+ti-1,work+fi-1,work+fni-1,work+fsi-1,work+fjacbi-1,
          work+fjacxi-1,w,wd,ldwd,work+ssfi-1,work+ssi-1,work+tti-1,
          ldtt,work+xplusi-1,work+ddelti-1,work+sssi-1,work,lwork,iwork,liwork);
    

    if (err >= 0) {
        ctx->ODRFNPtr = fni - 1;
        ctx->ODREPSPtr = fi - 1;
        ctx->ODRXPDPtr = xplusi - 1;
    }
    return(err);
}

/*  ----------------------------------------------------------------------- */    
/*     dacces()                                                             */
/*                                                                          */
/*  dacces(n,m,np,work,lwork,iwork,liwork,access,jpvt,wrk1,tfjacb,          */
/*      omega,yt,u,qraux,wrk2,nnzw,npp,partol,sstol,taufac,epsmac,      */
/*      neta,wss,wssdel,wsseps,rvar,idf,tau,alpha,     */
/*      int2,olmavg,rcond,irank,actrs,pnorm,prers,rnorms)                   */
/*                                                                          */
/*  PURPOSE  ACCESS OR STORE VALUES IN THE WORK ARRAYS                      */
/*                                                                          */

void dacces(TDAContext *ctx, int n,int m,int np,double *work,int lwork,int *iwork,int liwork, int access,int *jpvt,int *wrk1,int *tfjacb,int *omega,int *yt,int *u,int *qraux, int *wrk2,int *nnzw,int *npp,double *partol,double *sstol, double *taufac,double *epsmac,int *neta, double *wss,double *wssdel,double *wsseps,double *rvar,int *idf, double *tau,double *alpha,int *int2,double *olmavg,double *rcond, int *irank,double *actrs,double *pnorm,double *prers,double *rnorms)
{
    (void)liwork; (void)lwork;        /* unused: the signature is shared */
    int actrsi,alphai,betaci,betani,betasi,ddelti,deltai,deltni,deltsi;
    int epsmai,etai,fi,fjacbi,fjacxi,fni,fsi,idfi,int2i;               
    int iranki,jpvti,ldtti,liwkmn,lwkmn;
    int netai,nnzwi,nppi,nrowi;
    int olmavi,omegai,partli,pnormi,prersi,qrauxi,rcondi,rnorsi,rvari;
    int si,ssfi,ssi,sssi,sstoli,taufci,taui,tfjaci,ti,tti,ui,wrk1i;
    int wrk2i,wssi,wssdei,wssepi,xplusi,yti;
 
    /* FIND STARTING LOCATIONS WITHIN INTEGER WORKSPACE */
 
    diwinf(ctx, m,np,&jpvti,&nnzwi,&nppi,&idfi,&nrowi,&netai,
        &int2i,&iranki,&ldtti,&liwkmn);
 
    /* FIND STARTING LOCATIONS WITHIN DOUBLE PRECISION WORK SPACE */
 
    dwinf(ctx, n,m,np,&deltai,&fi,&wssi,&wssdei,&wssepi,&rvari,
         &partli,&sstoli,&taufci,&epsmai,&olmavi,
         &fjacbi,&fjacxi,&xplusi,&betaci,&betasi,&betani,&deltsi,
         &deltni,&ddelti,&fsi,&fni,&si,&sssi,&ssi,&ssfi,&ti,&tti,&taui,
         &alphai,&tfjaci,&omegai,&yti,&ui,&qrauxi,&wrk1i,&wrk2i,&rcondi,
         &etai,&actrsi,&pnormi,&prersi,&rnorsi,&lwkmn); 
 
    if (access) {   /* SET STARTING LOCATIONS FOR WORK VECTORS */
 
        *jpvt   = jpvti;
        *wrk1   = wrk1i;
        *tfjacb = tfjaci;
        *omega  = omegai;
        *yt     = yti;
        *u      = ui;
        *qraux  = qrauxi;
        *wrk2   = wrk2i;
 
        *actrs  = work[actrsi];
        *alpha  = work[alphai];
        *epsmac = work[epsmai];
        *olmavg = work[olmavi];
        *partol = work[partli];
        *pnorm  = work[pnormi];
        *prers  = work[prersi];
        *rcond  = work[rcondi];
        *wss    = work[wssi];
        *wssdel = work[wssdei];
        *wsseps = work[wssepi];
        *rvar   = work[rvari];
        *rnorms = work[rnorsi];
        *sstol  = work[sstoli];
        *tau    = work[taui];
        *taufac = work[taufci];
 
        *neta   = iwork[netai];
        *irank  = iwork[iranki];
        *nnzw   = iwork[nnzwi];
        *npp    = iwork[nppi];
        *idf    = iwork[idfi];
        *int2   = iwork[int2i];
                              
    }
    else {  /* STORE VALUES INTO THE WORK VECTORS */
 
        work[actrsi]  = *actrs;
        work[alphai]  = *alpha;
        work[olmavi]  = *olmavg;
        work[partli]  = *partol;
        work[pnormi]  = *pnorm; 
        work[prersi]  = *prers;
        work[rcondi]  = *rcond;
        work[wssi]    = *wss;
        work[wssdei]  = *wssdel;
        work[wssepi]  = *wsseps;
        work[rvari]   = *rvar;
        work[rnorsi]  = *rnorms;
        work[sstoli]  = *sstol;
        work[taui]    = *tau;
 
        iwork[iranki] = *irank;
        iwork[idfi]   = *idf;
        iwork[int2i]  = *int2;
    }
}

/*  ----------------------------------------------------------------------- */    
/*  dchex()                                                                 */
/*                                                                          */ 
/*  dchex(r,cdr,p,k,l,z,cdz,nz,c,s,job)                                     */
/*                                                                          */
/*  PURPOSE  UPDATES THE CHOLESKY FACTORIZATION  A=TRANS(R)*R  OF A         */
/*           POSITIVE DEFINITE MATRIX A OF ORDER P UNDER DIAGONAL           */
/*           PERMUTATIONS OF THE FORM  TRANS(E)*A*E  WHERE E IS A           */
/*           PERMUTATION MATRIX.                                            */
/*                                                                          */

void dchex(TDAContext *ctx, double *r,int cdr,int p,int k,int l,double *z,int cdz,int nz, double *c,double *s,int job)
{
    double t,t1;
    int i,ii,il,iu,j,jj,km1,kp1,lm1,lmk;
 
    km1 = k - 1;
    kp1 = k + 1;
    lmk = l - k;
    lm1 = l - 1;
 
    if (job == 2)
        goto L130;
 
    /* RIGHT CIRCULAR SHIFT. */
 
    for (i = 1; i <= l; ++i) {
        ii = l - i + 1;
        s[i] = r[(ii-1)*cdr+l]; /* R(ii,l); */ 
    }
    for (jj = k; jj <= lm1; ++jj) {
        j = lm1 - jj + k;
        for (i = 1; i <= j; ++i)  
            r[(i-1)*cdr+j+1] = r[(i-1)*cdr+j];
        r[j*cdr+j+1] = 0.0;
    }
    for (i = 1; i < k; ++i) {
        ii = l - i + 1;
        r[(i-1)*cdr+k] = s[ii];
    }
 
    /* CALCULATE THE ROTATIONS. */
 
    t = s[1];
    for (i = 1; i <= lmk; ++i) {
        t1 = s[i];
        drotg(ctx, &s[i+1],&t,&c[i],&t1);
        s[i] = t1;
        t = s[i+1];
    }
    r[(k-1)*cdr+k] = t;
    for (j = kp1; j <= p; ++j) {
        il = imax(ctx, 1,l-j+1); 
        for (ii = il; ii <= lmk; ++ii) {
            i = l - ii;
            t = c[ii] * r[(i-1)*cdr+j] + s[ii] * r[i*cdr+j];
            r[i*cdr+j] = c[ii] * r[i*cdr+j] - s[ii] * r[(i-1)*cdr+j];
            r[(i-1)*cdr+j] = t;
        }
    }
 
    /* IF REQUIRED, APPLY THE TRANSFORMATIONS TO Z. */
 
    if (nz  <  1)             
        return;    

    for (j = 1; j <= nz; ++j) {
        for (ii = 1; ii <= lmk; ++ii) {
            i = l - ii;

            t = c[ii] * z[(i-1)*cdz+j] + s[ii] * z[i*cdz+j];
            z[i*cdz+j] = c[ii] * z[i*cdz+j] - s[ii] * z[(i-1)*cdz+j];
            z[(i-1)*cdz+j] = t;
        }
    }
    return;
 
L130:   /* LEFT CIRCULAR SHIFT */
 
    /* REORDER THE COLUMNS */
 
    for (i = 1; i <= k; ++i) {
        ii = lmk + i;
        s[ii] = r[(i-1)*cdr+k];
    }
    for (j = k; j <= lm1; ++j) {
        for (i = 1; i <= j; ++i)  
            r[(i-1)*cdr+j] = r[(i-1)*cdr+j+1];
        jj = j - km1;
        s[jj] = r[j*cdr+j+1];
    }
    for (i = 1; i <= k; ++i) {
        ii = lmk + i;
        r[(i-1)*cdr+l] = s[ii];
    }
    for (i = kp1; i <= l; ++i)
        r[(i-1)*cdr+l] = 0.0;
 
    /* REDUCTION LOOP. */

    for (j = k; j <= p; ++j) {
        if (j  !=  k) {           
                                /* APPLY THE ROTATIONS. */
            iu = imin(ctx, j-1,l-1);
            for (i = k; i <= iu; ++i) {
                ii = i - k + 1;
                t = c[ii] * r[(i-1)*cdr+j] + s[ii] * r[i*cdr+j];
                r[i*cdr+j] = c[ii] * r[i*cdr+j] - s[ii] * r[(i-1)*cdr+j];
                r[(i-1)*cdr+j] = t;
            }
        }
        if (j  < l) {            
            jj = j - k + 1;
            t = s[jj];
            drotg(ctx, &r[(j-1)*cdr+j],&t,&c[jj],&s[jj]);
        }
    }
          
    /* APPLY THE ROTATIONS TO Z. */
 
    if (nz  <  1)           
        return;   

    for (j = 1; j <= nz; ++j) {
        for (i = k; i <= lm1; ++i) {
            ii = i - km1;
            t = c[ii] * z[(i-1)*cdz+j] + s[ii] * z[i*cdz+j];
            z[i*cdz+j] = c[ii] * z[i*cdz+j] - s[ii] * z[(i-1)*cdz+j];
            z[(i-1)*cdz+j] = t;
        }
    }
}             

/*  ----------------------------------------------------------------------- */    
/*  ddiagi()                                                                */
/*                                                                          */
/*  ddiagi(n,m,s,lds,v,cdv,sv,ddsv)                                         */
/*                                                                          */
/*  PURPOSE  SCALE THE VECTOR V BY THE INVERSE OF THE DIAGONAL MATRIX S     */
/*           AND RETURN THE RESULT IN VECTOR SV                             */
/*                                                                          */

void ddiagi(TDAContext *ctx, int n,int m,double *s,int lds,double *v,int cdv,double *sv,int cdsv)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j;
 
    if (n == 0 || m == 0)        
        return;
 
    if (s[1] < 0.0) {    
        for (j = 1; j <= m; ++j) {
            for (i = 1; i <= n; ++i)
                sv[(i-1)*cdsv+j] = v[(i-1)*cdv+j] / fabs(s[1]);
        }
    }
    else {
        if (lds == 1) {    
            for (j = 1; j <= m; ++j) {
                for (i = 1; i <= n; ++i) 
                    sv[(i-1)*cdsv+j] = v[(i-1)*cdv+j]/s[j];
            }
        }
        else {
            for (j = 1; j <= m; ++j) {
                for (i = 1; i <= n; ++i) 
                    sv[(i-1)*cdsv+j] = v[(i-1)*cdv+j]/s[(i-1)*m+j];
            }
        }          
    }           
}

/*  ----------------------------------------------------------------------- */    
/*  ddiags()                                                                */
/*                                                                          */
/*  ddiags(n,m,s,lds,v,cdv,sv,cdsv)                                         */
/*                                                                          */
/*  PURPOSE  SCALE THE VECTOR V BY THE DIAGONAL MATRIX S                    */
/*           AND RETURN THE RESULT IN VECTOR SV.                            */
/*                                                                          */

void ddiags(TDAContext *ctx, int n,int m,double *s,int lds,double *v,int cdv,double *sv,int cdsv)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j;

    if (n == 0 || m == 0)          
        return;
 
    if (s[1] < 0.0) {    
        for (j = 1; j <= m; ++j) {
            for (i = 1; i <= n; ++i) 
                sv[(i-1)*cdsv+j] = fabs(s[1])*v[(i-1)*cdv+j];
        }
    }
    else { 
        if (lds == 1) {       
            for (j = 1; j <= m; ++j) {
                for (i = 1; i <= n; ++i) 
                    sv[(i-1)*cdsv+j] = s[j]*v[(i-1)*cdv+j];    
            }
        }
        else {
            for (j = 1; j <= m; ++j) {
                for (i = 1; i <= n; ++i)
                    sv[(i-1)*cdsv+j] = s[(i-1)*m+j]*v[(i-1)*cdv+j];
            }
        }           
    }          
}

/*  ----------------------------------------------------------------------- */    
/*  ddiagw()                                                                */
/*                                                                          */
/*  ddiagw(n,m,w,v,cdv,wv,cdwv)                                             */
/*                                                                          */
/*  PURPOSE  SCALE THE N BY M ARRAY V BY THE DIAGONAL OBSERVATIONAL         */
/*           ERROR WEIGHT MATRIX W AND RETURN THE RESULT IN VECTOR WV.      */
/*           N.B.  IF THE FIRST ELEMENT OF W IS NEGATIVE, THE DEFAULT       */
/*           WEIGHTING OF ONE FOR ALL ELEMENTS WILL BE INVOKED, I.E.,       */
/*           THE RESULTS WILL BE "UNWEIGHTED."                              */
/*                                                                          */

void ddiagw(TDAContext *ctx, int n,int m,double *w,double *v,int cdv,double *wv,int cdwv)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j;
 
    if (n == 0 || m == 0)          
        return;
 
    if (w[1] < 0.0) {    
        for (j = 1; j <= m; ++j) {
            for (i = 1; i <= n; ++i)
                wv[(i-1)*cdwv+j] = v[(i-1)*cdv+j];
        }
    }
    else {
        for (j = 1; j <= m; ++j) {
            for (i = 1; i <= n; ++i)
                wv[(i-1)*cdwv+j] = w[i]*v[(i-1)*cdv+j];  
        }
    }
}

/*  ----------------------------------------------------------------------- */    
/*  detaf()                                                                 */
/*                                                                          */
/*  detaf(n,np,m,xplusd,beta,eta,neta,epsmac,nrow,partmp,pvtemp)            */                                                
/*                                                                          */
/*  COMPUTE NOISE AND NUMBER OF GOOD DIGITS IN FUNCTION RESULTS             */
/*  (THIS ROUTINE IS MODELED AFTER STARPAC SUBROUTINE ETAFUN)               */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int detaf(TDAContext *ctx, int n,int np,int m,double *xplusd,double *beta,double *eta, int *neta,double epsmac,int nrow,double *partmp,double *pvtemp)
{
    register int i,k;
    double a,b,fac,j,rsssm,rsssmj,sqrtmp;     
    double rss[6];
    double p1 = 0.1;
    double p2 = 0.2;
 
    sqrtmp = sqrt(epsmac);
    rsssm = 0.0;
    rsssmj = 0.0;
    for (i = 1; i <= 5; ++i) {
        j = (double)(i - 3);
        for (k = 1; k <= np; ++k)  
            partmp[k] = beta[k] * (1.0 + j * sqrtmp);
                        
        ctx->ODRNF += 1;
        if (odr_fun(ctx, n,np,m,partmp,xplusd,pvtemp,1))
            return(-1);
 
        rss[i] = pvtemp[nrow];
 
        rsssm += rss[i];
        rsssmj += j * rss[i];
    }
    a = p2 * rsssm;
    b = p1 * rsssmj;
    if (rss[3] != 0.0)       
        fac = 1.0 / fabs(rss[3]);
    else   
        fac = 1.0;

    for (i = 1; i <= 5; ++i) {
        j = (double)(i - 3);
        rss[i] = fabs((rss[i] - (a + j * b)) * fac);
    }
    *eta = epsmac;         
    if (*eta < rss[1])
        *eta = rss[1];
    if (*eta < rss[2])
        *eta = rss[2];
    if (*eta < rss[3])
        *eta = rss[3];
    if (*eta < rss[4])
        *eta = rss[4];
    if (*eta < rss[5])
        *eta = rss[5];
           
    *neta = (int)(-rlog(ctx, *eta) / log(10.0));
    return(0);
}          

/*  ----------------------------------------------------------------------- */    
/*  devfun()                                                                */
/*                                                                          */
/*  devfun(n,np,m,betac,beta,ifixb,x,y,delta,xplusd,w,f)                    */                                         
/*                                                                          */
/*  COMPUTE THE WEIGHTED EPSILON'S FOR THE CURRENT POINT                    */
/*                                                                          */
/*  Return:   0 if OK, -1 if error in odr_fun().                            */
/*                                                                          */

int devfun(TDAContext *ctx, int n,int np,int m,double *betac,double *beta,char *ifixb, double *x,double *y,double *delta,double *xplusd,double *w,double *f)
{
    int err;   

    /* INSERT CURRENT UNFIXED BETA ESTIMATES INTO BETA */
 
    dunpac(ctx, np,betac,beta,ifixb);
 
    /* COMPUTE XPLUSD = X + DELTA */
   
    dxpy(ctx, n,m,x,delta,xplusd);
 
    /* EVALUATE THE PREDICTED VALUES OF THE FUNCTION FOR THE CURRENT POINT */
 
    err = odr_fun(ctx, n,np,m,beta,xplusd,f,1);
    ctx->ODRNF++;
    if (err)      
        return(err);
 
    /* COMPUTE WEIGHTED EPSILONS FOR CURRENT POINT AND STORE IN F */
    /* f(i) = f(i) - y(i) */
   
    daxpy(ctx, n,-1.0,y,1,f,1);
    ddiagw(ctx, n,1,w,f,1,f,1);
    return(0);
}          

/*  ----------------------------------------------------------------------- */    
/*  devjac()                                                                */
/*                                                                          */
/*  devjac(n,np,npp,m,betac,beta,ifixb,ifixx,ifixu,x,delta,      */
/*      xplusd,ss,tt,ldtt,neta,pv,stp,fjacb,             */
/*      fjacx,w)                                                */
/*                                                                          */
/*  COMPUTE THE WEIGHTED JACOBIANS WRT BETA AND DELTA                       */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */ 

int devjac(TDAContext *ctx, int n,int np,int npp,int m,double *betac,double *beta, char *ifixb,char *ifixx,int ifixu,double *x,double *delta, double *xplusd,double *ss,double *tt,int ldtt, int neta,double *pv,double *stp,double *fjacb, double *fjacx,double *w)
{
    int i,j,jfx,err;
 
    /* INSERT CURRENT UNFIXED BETA ESTIMATES INTO BETA */
  
    dunpac(ctx, np,betac,beta,ifixb);
 
    /* COMPUTE XPLUSD = X + DELTA */
 
    dxpy(ctx, n,m,x,delta,xplusd);          

    /* COMPUTE THE JACOBIAN WRT THE ESTIMATED BETAS (FJACB) AND
       THE JACOBIAN WRT DELTA (FJACX) */

    if (ctx->ODRAd) {    
        err = odr_jac(ctx, n,np,m,beta,xplusd,fjacb,fjacx);
        ctx->ODRNFJ += 1;       
    }
    else {
        err = djacfd(ctx, n,np,m,beta,x,delta,xplusd,ss,tt,ldtt,neta,pv,
              stp,ifixb,fjacb,ifixx,ifixu,fjacx);
    }      
    if (err)
        return(-1);
 
   /* WEIGHT THE JACOBIAN WRT THE ESTIMATED BETAS */
 
    if (ctx->ODRAd) {     
        jfx = 0;
        for (j = 1; j <= np; ++j) {
            if (ifixb[j]) {     
                jfx += 1;    
                ddiagw(ctx, n,1,w,fjacb+j-1,np,fjacb+jfx-1,np);
            }
        }
    }
    else {
        for (j = 1; j <= npp; ++j)  
            ddiagw(ctx, n,1,w,fjacb+j-1,np,fjacb+j-1,np);
    }          

    /* WEIGHT OR ZERO THE JACOBIAN'S WRT X AS APPROPRIATE */
 
    if (ctx->ODRTyp) {   
        if (ifixu) { /* if ifixx is used    
                        CHECK FOR POSSIBLY FIXED COLUMNS OR ELEMENTS OF X */
 
            /* WEIGHT JACOBIAN WRT X(I,J) FOR I=1,N AND
               THEN ZERO APPROPRIATE ELEMENTS */
 
            for (j = 1; j <= m; ++j) {
                ddiagw(ctx, n,1,w,fjacx+j-1,m,fjacx+j-1,m);
                for (i = 1; i <= n; ++i) {
                    if (ifixx[(i-1)*m+j] == 0)        
                        fjacx[(i-1)*m+j] = 0.0;
                }
            }
        }
        else {      /* WEIGHT JACOBIAN WRT X(I,J) FOR I=1,N AND J=1,M */
    
            for (j = 1; j <= m; ++j)  
                ddiagw(ctx, n,1,w,fjacx+j-1,m,fjacx+j-1,m);
        }           
    }
    else  /* zero all of fjacx for OLS */
        dzero(ctx, n,m,fjacx,m);            
    return(0);
}            
          
/*  ----------------------------------------------------------------------- */    
/*  didts()                                                                 */
/*                                                                          */
/*  didts(n,m,w,wd,ldwd,alpha,tt,ldtt,t,dtt)                                */
/*                                                                          */
/*           SCALE MATRIX TT BY THE INVERSE OF DT, I.E., COMPUTE            */
/*           DTT = T * INV(DT) WHERE DT = (W*D)**2 + ALPHA*TT**2,           */
/*           W AND D ARE DEFINED BY EQ.2 OF THE PROLOGUE OF DODR            */
/*           AND DODRC, AND TT IS THE SCALING MATRIX FOR THE DELTA'S,       */
/*           ALSO DEFINED IN THE PROLOGUE OF DODR AND DODRC.                */
/*                                                                          */

void didts(TDAContext *ctx, int n,int m,double *w,double *wd,int ldwd,double alpha,double *tt, int ldtt,double *t,double *dtt)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j;
    double dt,term1,term2;     
 
    if (n == 0 || m == 0)        
        return;

    if (w[1] >= 0.0) {     
        if (wd[1] > 0.0) {   
            if (ldwd >= n) {    
                if (tt[1] > 0.0) {    
                    if (ldtt >= n) {    
                        for (j = 1; j <= m; ++j) {
                            for (i = 1; i <= n; ++i) {
                                if (w[i] != 0.0 || alpha != 0.0) {     
                                    dtt[(i-1)*m+j] = t[(i-1)*m+j] /
                                           (pow(w[i]*wd[(i-1)*m+j],2.0) +
                                            pow(alpha*tt[(i-1)*m+j],2.0));
                                }
                                else {
                                    dtt[(i-1)*m+j] = 0.0;
                                }     
                            }
                        }
                    }
                    else {
                        for (j = 1; j <= m; ++j) {
                            term2 = pow(alpha*tt[j],2.0);
                            for (i = 1; i <= n; ++i) {
                                if (w[i] != 0.0 || alpha != 0.0) {    
                                    dtt[(i-1)*m+j] = t[(i-1)*m+j]/
                                         (pow(w[i]*wd[(i-1)*m+j],2.0)+term2);
                                }
                                else  
                                    dtt[(i-1)*m+j] = 0.0;
                            }
                        }
                    }     
                }
                else {
                    term2 = pow(alpha*tt[1],2.0);
                        for (j = 1; j <= m; ++j) {
                            for (i = 1; i <= n; ++i) {
                                if (w[i] != 0.0 || alpha != 0.0) {    
                                    dtt[(i-1)*m+j] = t[(i-1)*m+j]/
                                        (pow(w[i]*wd[(i-1)*m+j],2.0)+term2);
                                }
                                else  
                                    dtt[(i-1)*m+j] = 0.0;
                            }
                        }
                }     
            }
            else {
                if (tt[1] > 0.0) {    
                    if (ldtt >= n) {   
                        for (j = 1; j <= m; ++j) {
                            for (i = 1; i <= n; ++i) {
                                if (w[i] != 0.0 || alpha != 0.0) {   
                                    dtt[(i-1)*m+j] = t[(i-1)*m+j]/
                                           (pow(w[i]*wd[j],2.0) +
                                            pow(alpha*tt[(i-1)*m+j],2.0));
                                }
                                else  
                                    dtt[(i-1)*m+j] = 0.0;
                            }
                        }
                    }
                    else {
                        for (j = 1; j <= m; ++j) {
                            term2 = pow(alpha*tt[j],2.0);
                            for (i = 1; i <= n; ++i) {
                                if (w[i] != 0.0 || alpha != 0.0) {    
                                    dtt[(i-1)*m+j] = t[(i-1)*m+j]/
                                          (pow(w[i]*wd[j],2.0)+term2);
                                }
                                else  
                                    dtt[(i-1)*m+j] = 0.0;
                            }
                        }   
                    }     
                }
                else {
                    term2 = pow(alpha*tt[1],2.0);
                    for (j = 1; j <= m; ++j) {
                        for (i = 1; i <= n; ++i) {
                            if (w[i] != 0.0 || alpha != 0.0) {    
                                dtt[(i-1)*m+j] = t[(i-1)*m+j]/
                                           (pow(w[i]*wd[j],2.0) + term2);
                            }
                            else  
                                dtt[(i-1)*m+j] = 0.0;
                        }
                    }   
                }    
            }      
        }
        else {
            if (tt[1] > 0.0) {    
                if (ldtt >= n) {   
                    for (j = 1; j <= m; ++j) {
                        for (i = 1; i <= n; ++i) {
                            if (w[i] != 0.0 || alpha != 0.0) {   
                                dtt[(i-1)*m+j] = t[(i-1)*m+j]/
                                            (pow(w[i]*wd[1],2.0) +
                                             pow(alpha*tt[(i-1)*m+j],2.0));
                            }
                            else  
                                dtt[(i-1)*m+j] = 0.0;
                        }
                    }
                }
                else {
                    for (j = 1; j <= m; ++j) {
                        term2 = pow(alpha*tt[j],2.0);
                        for (i = 1; i <= n; ++i) {
                            if (w[i] != 0.0 || alpha != 0.0) {     
                                dtt[(i-1)*m+j] = t[(i-1)*m+j]/
                                         (pow(w[i]*wd[1],2.0)+term2);
                            }
                            else  
                                dtt[(i-1)*m+j] = 0.0;
                        }
                    }
                }     
            }
            else {
                term2 = pow(alpha*tt[1],2.0);
                for (j = 1; j <= m; ++j) {
                    for (i = 1; i <= n; ++i) {
                        if (w[i] != 0.0 || alpha != 0.0) {     
                            dtt[(i-1)*m+j] = t[(i-1)*m+j]/
                                      (pow(w[i]*wd[1],2.0)+term2);
                        }
                        else  
                            dtt[(i-1)*m+j] = 0.0;
                    }
                }
            }      
        }          
    }
    else {
        if (wd[1] > 0.0) {    
            if (ldwd >= n) {   
                if (tt[1] > 0.0) {    
                    if (ldtt >= n) {     
                        for (j = 1; j <= m; ++j) {
                            for (i = 1; i <= n; ++i) {
                                dtt[(i-1)*m+j] = t[(i-1)*m+j]/
                                          (pow(wd[(i-1)*m+j],2.0) +
                                           pow(alpha*tt[(i-1)*m+j],2.0));
                            }   
                        }
                    }
                    else {
                        for (j = 1; j <= m; ++j) {
                            term2 = pow(alpha*tt[j],2.0);
                            for (i = 1; i <= n; ++i) {
                                dtt[(i-1)*m+j] = t[(i-1)*m+j]/
                                           (pow(wd[(i-1)*m+j],2.0) + term2);
                            }   
                        }
                    }     
                }
                else {
                    term2 = pow(alpha*tt[1],2.0);
                    for (j = 1; j <= m; ++j) {
                        for (i = 1; i <= n; ++i) {
                            dtt[(i-1)*m+j] = t[(i-1)*m+j]/
                                          (pow(wd[(i-1)*m+j],20) + term2);
                        }
                    }
                }      
            }
            else {
                if (tt[1] > 0.0) {     
                    if (ldtt >= n) {   
                        for (j = 1; j <= m; ++j) {
                            term1 = pow(wd[j],2.0);
                            for (i = 1; i <= n; ++i) {
                                dtt[(i-1)*m+j] = t[(i-1)*m+j]/
                                       (term1+alpha*pow(tt[(i-1)*m+j],2.0));
                            }
                        }
                    }
                    else {
                        for (j = 1; j <= m; ++j) {
                            dt = 1.0 / (pow(wd[j],2.0)+alpha*pow(tt[j],2.0));
                            for (i = 1; i <= n; ++i) {
                                dtt[(i-1)*m+j] = t[(i-1)*m+j]*dt;
                            }
                        }
                    }      
                }
                else {
                    term2 = pow(alpha*tt[1],2.0);
                    for (j = 1; j <= m; ++j) {
                        term1 = pow(wd[j],2.0);
                        dt = 1.0/(term1+term2);
                        for (i = 1; i <= n; ++i) {
                            dtt[(i-1)*m+j] = t[(i-1)*m+j]*dt;
                        }
                    }
                }      
            }      
        }
        else {
            if (tt[1] > 0.0) {   
                if (ldtt >= n) {    
                    term1 = pow(wd[1],2.0);
                    for (j = 1; j <= m; ++j) {
                        for (i = 1; i <= n; ++i) {
                            dtt[(i-1)*m+j] = t[(i-1)*m+j]/
                                  (term1 + alpha*pow(tt[(i-1)*m+j],2.0));
                        }
                    }
                }
                else {
                    term1 = pow(wd[1],2.0);
                    for (j = 1; j <= m; ++j) {
                        term2 = alpha*pow(tt[j],2.0);
                        dt = 1.0/(term1+term2);
                        for (i = 1; i <= n; ++i) {
                            dtt[(i-1)*m+j] = t[(i-1)*m+j]*dt;
                        }
                    }
                }     
            }
            else {
                dt = 1.0/(pow(wd[1],2.0) +alpha* pow(tt[1],2.0));
                for (j = 1; j <= m; ++j) {
                    for (i = 1; i <= n; ++i) {
                        dtt[(i-1)*m+j] = t[(i-1)*m+j]*dt;
                    }
                }
            }      
        }      
    }       
}
            
/*  ----------------------------------------------------------------------- */    
/*  diniwk()    initialization.                                             */
/*                                                                          */
/*  diniwk(n,m,np,work,iwork,x,beta,sstol,partol,taufac,epsmai,sstoli,      */
/*      partli,taufci,ssfi,tti,ldtti)                                       */
/*                                                                          */
  
void diniwk(TDAContext *ctx, int n,int m,int np,double *work,int *iwork,double *x,double *beta,                                    double sstol,double partol,double taufac,int epsmai,int sstoli,int partli, int taufci,int ssfi,int tti,int ldtti)
{
    work[epsmai] = ctx->EPSI;    /* machine epsilon */
 
    if (partol < work[epsmai] || partol >= 1.0)      
         work[partli] = pow(work[epsmai],(2.0/3.0));
    else  
        work[partli] = partol;
 
    if (sstol < work[epsmai] || sstol >= 1.0)        
        work[sstoli] = sqrt(work[epsmai]);
    else  
        work[sstoli] = sstol;
 
    if (taufac <= 0.0)      
        work[taufci] = 1.0;
    else
        work[taufci] = taufac;

    /* COMPUTE SCALING FOR BETA'S AND DELTA'S */
        
    dsclb(ctx, np,beta,work+ssfi-1);
    iwork[ldtti] = n;
    dscld(ctx, n,m,x,work+tti-1);          
}

/*  ----------------------------------------------------------------------- */    
/*  djacfd()                                                                */
/*                                                                          */
/*  djacfd(n,np,m,beta,x,delta,xplusd,ss,tt,ldtt,neta,pv,stp,         */
/*      ifixb,fjacb,ifixx,ifixu,fjacx)                                      */     
/*                                                                          */
/*  COMPUTE FINITE DIFFERENCE APPROXIMATIONS TO THE                         */
/*  JACOBIAN WRT THE ESTIMATED BETAS AND WRT THE DELTAS                     */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */ 

int djacfd(TDAContext *ctx, int n,int np,int m,double *beta,double *x,double *delta, double *xplusd,double *ss,double *tt,int ldtt,int neta,double *pv, double *stp,char *ifixb,double *fjacb,char *ifixx,int ifixu,double *fjacx)
{
    int i,j,jfx,doit,err;
    double betaj,sqreps,typj;     
 
    /* SET THE RELATIVE STEP SIZE FOR COMPUTING THE JACOBIANS */
 
    sqreps = pow(10.0,(double)(-neta) / 2.0);
 
    /* COMPUTE THE PREDICTED VALUES OF THE FUNCTION AT THE GIVEN POINT */
 
    err = odr_fun(ctx, n,np,m,beta,xplusd,pv,1);
    ctx->ODRNF += 1;      
    if (err)             
        return(-1);

    /* COMPUTE THE JACOBIAN WRT THE ESTIMATED BETAS */
 
    jfx = 0;
    for (j = 1; j <= np; ++j) {
        if (ifixb[j]) {                 /* if free parameter */
            jfx += 1;     
            betaj = beta[j];
            typj = 1.0/ss[jfx];
            stp[j] = betaj + sqreps*fsign(ctx, 1.0,betaj)*dmax(ctx, fabs(betaj),typj);
            stp[j] = stp[j] - betaj;
            beta[j] = betaj + stp[j];
            err = odr_fun(ctx, n,np,m,beta,xplusd,fjacb + jfx-1,np);
            ctx->ODRNF += 1;      
            if (err)      
                return(-1);

            for (i = 1; i <= n; ++i)  
                fjacb[(i-1)*np + jfx] = (fjacb[(i-1)*np+jfx]-pv[i])/stp[j];

            beta[j] = betaj;
        }          
    }
   
    /* COMPUTE THE JACOBIAN WRT THE X'S */
 
    if (ctx->ODRTyp) {    
        for (j = 1; j <= m; ++j) {
            doit = 0;
            if (ifixu == 0)         /* if ifixx is not used */
               doit = 1;    
            else {
                for (i = 1; i <= n; ++i) {
                    if (ifixx[(i-1)*m+j] != 0) {     
                        doit = 1;
                        break;   
                    }        
                }
            }   
            if (doit == 0)  
                dzero(ctx, n,1,fjacx+j-1,m);            
              
            else {
                for (i = 1; i <= n; ++i) {
                    if (tt[1] > 0.0) {     
                        if (ldtt == 1)       
                            typj = 1.0/tt[j];
                        else 
                            typj = 1.0/tt[(i-1)*m+j];
                    }
                    else
                        typj = fabs(1.0/tt[1]);
                         
                    stp[i] = xplusd[(i-1)*m+j] +
                                sqreps * fsign(ctx, 1.0,xplusd[(i-1)*m+j]) *
                                         dmax(ctx, fabs(xplusd[(i-1)*m+j]),typj);

                    stp[i] = xplusd[(i-1)*m+j] + sqreps*fsign(ctx, 1.0,xplusd[(i-1)*m+j])*
                                             dmax(ctx, fabs(xplusd[(i-1)*m+j]),typj);
                    stp[i] -= xplusd[(i-1)*m+j];
                    xplusd[(i-1)*m+j] += stp[i];
                }
                err = odr_fun(ctx, n,np,m,beta,xplusd,fjacx + j-1,m);
                ctx->ODRNF += 1;       
                if (err)        
                    return(-1);

                for (i = 1; i <= n; ++i) {
                    fjacx[(i-1)*m+j] = (fjacx[(i-1)*m+j] - pv[i]) / stp[i];
                    xplusd[(i-1)*m+j] = x[(i-1)*m+j] + delta[(i-1)*m+j];
                }
            }      
        }
    }        
    return(0);
}

/*  ----------------------------------------------------------------------- */    
/*  dodlm()                                                                 */
/*                                                                          */
/*  dodlm(n,np,npp,m,f,fjacb,fjacx,w,wd,ldwd,ss,tt,ldtt,ddelt,              */
/*      alpha2,tau,epsmac,sss,wrk1,tfjacb,omega,yt,u,qraux,wrk2,jpvt,       */
/*      s,t,nlms,rcond,irank)                                               */
/*                                                                          */
/*  COMPUTE LEVENBERG-MARQUARDT PARAMETER AND STEPS S AND T USING           */
/*  ANALOG OF THE TRUST-REGION LEVENBERG-MARQUARDT ALGORITHM.               */
/*                                                                          */

void dodlm(TDAContext *ctx, int n,int np,int npp,int m,double *f,double *fjacb,               double *fjacx,double *w,double *wd,int ldwd,double *ss, double *tt,int ldtt,double *ddelt,double *alpha2,double *tau,double epsmac, double *sss,double *wrk1,double *tfjacb,double *omega,double *yt, double *u,double *qraux,double *wrk2,int *jpvt,double *s,double *t, int *nlms,double *rcond,int *irank)
{
    int i,j;
    double alpha1,alphan,bot,phi1,phi2,sa,top;    
    double p001 = 0.001;
    double p1 = 0.1;
 
    /* COMPUTE FULL GAUSS-NEWTON STEP (ALPHA=0) */
 
    alpha1 = 0.0;
    dodstp(ctx, n,np,npp,m,f,fjacb,fjacx,w,wd,ldwd,ss,tt,ldtt,ddelt,alpha1,
        epsmac,sss,tfjacb,wrk1,omega,yt,u,qraux,wrk2,
        jpvt,s,t,&phi1,irank,rcond);

    /* INITIALIZE TAU IF NECESSARY */

    if (*tau < 0.0)        
        *tau = fabs(*tau)*phi1;
 
    /* CHECK IF FULL GAUSS-NEWTON STEP IS OPTIMAL */
 
    if ((phi1 - *tau) <= p1 * (*tau)) {    
        *nlms = 1;
        *alpha2 = 0.0;
        return;
    }             
 
    /* FULL GAUSS-NEWTON STEP IS OUTSIDE TRUST REGION -
       FIND LOCALLY CONSTRAINED OPTIMAL STEP */
 
    phi1 -= *tau;

    /* INITIALIZE UPPER AND LOWER BOUNDS FOR ALPHA  */
 
    bot = 0.0;
    if (npp >= 1) {    
        for (i = 1; i <= npp; ++i) {
            sss[i] = ddot(ctx, n,fjacb + i - 1,np,f,1);
        }
        ddiagi(ctx, npp,1,ss,npp,sss,1,sss,1);
    }      
    for (j = 1; j <= m; ++j) {
        for (i = 1; i <= n; ++i)  
            wrk1[(i-1)*m+j] = fjacx[(i-1)*m+j]*f[i] + ddelt[(i-1)*m+j];
    }
    ddiagi(ctx, n,m,tt,ldtt,wrk1,m,sss+npp,m);
    top = dnrm2(ctx, npp+n*m,sss,1) / *tau;
    if (*alpha2 > top || *alpha2 == 0.0)       
        *alpha2 = p001*top;

    /* MAIN LOOP */
 
    for (i = 1; i <= 10; ++i) {

        /*  COMPUTE LOCALLY CONSTRAINED STEPS S AND T AND PHI(ALPHA) FOR
            CURRENT VALUE OF ALPHA */
 
        dodstp(ctx, n,np,npp,m,f,fjacb,fjacx,w,wd,ldwd,ss,tt,ldtt,
            ddelt,*alpha2,epsmac,sss,tfjacb,wrk1,omega,yt,u,qraux,wrk2,
            jpvt,s,t,&phi2,irank,rcond);

        phi2 -= *tau;
 
        /* CHECK WHETHER CURRENT STEP IS OPTIMAL */
 
        if (fabs(phi2) <= p1 * *tau ||(*alpha2 == bot && phi2 < 0.0)) {    
            *nlms = i+1;
            return;
        }      
 
        /* CURRENT STEP IS NOT OPTIMAL */
 
        /* UPDATE BOUNDS FOR ALPHA AND COMPUTE NEW ALPHA */
 
        if (phi1-phi2 == 0.0) {   
            *nlms = 12;
            return;
        }                     
        sa = phi2*(alpha1 - *alpha2) / (phi1-phi2);
        if (phi2 < 0.0)      
            top = dmin(ctx, top,*alpha2);
        else
            bot = dmax(ctx, bot,*alpha2);

        if (phi1*phi2 > 0.0)     
            bot = dmax(ctx, bot,*alpha2-sa);
        else
            top = dmin(ctx, top,*alpha2-sa);

        alphan = *alpha2 - sa*(phi1 + *tau) / *tau;

        if (alphan >= top || alphan <= bot)      
            alphan = dmax(ctx, p001*top,sqrt(top*bot));

        /* GET READY FOR NEXT ITERATION */
 
        alpha1 = *alpha2;
        *alpha2 = alphan;
        phi1 = phi2;
    }
 
    /* SET NLMS TO INDICATE AN OPTIMAL STEP COULD NOT BE FOUND IN 10 TRYS */
  
    *nlms = 12;
}                

/*  ----------------------------------------------------------------------- */    
/*  dodmn()                                                                 */
/*                                                                          */
/*  dodmn(n,np,m,x,ifixx,ifixu,y,betac,ifixb,beta,betan,betas,s,            */
/*      delta,deltan,deltas,t,f,fn,fs,fjacb,fjacx,w,wd,ldwd,                */
/*      ssf,ss,tt,ldtt,xplusd,ddelt,sss,work,lwork,iwork,liwork)            */
/*                                                                          */
/*  Iteratively compute least squares solution.                             */
/*                                                                          */
/*  return  1   squared convergence                                         */
/*          2   parameter convergence                                       */
/*          3   both, 1 and 2.                                              */
/*          0   exceeded max number of iterations.                          */
/*         -8   error in step size search                                   */
/*         -9   error in function evaluation                                */
/*                                                                          */

int dodmn(TDAContext *ctx, int n,int np,int m,double *x,char *ifixx,int ifixu,double *y, double *betac,char *ifixb,double *beta,double *betan,double *betas, double *s,double *delta,double *deltan,double *deltas,double *t, double *f,double *fn,double *fs,double *fjacb,double *fjacx, double *w,double *wd,int ldwd,double *ssf,double *ss, double *tt,int ldtt,double *xplusd,double *ddelt,double *sss, double *work,int lwork,int *iwork,int liwork)
{
    int info,r,i,j,idf,int2,irank,jpvt,junfix,neta,nlms,niter;
    int nnzw,npp,omega,qraux,tfjacb,u,wrk1,wrk2;
    int yt,access,cnvpar,cnvss,didvcv,intdbl,lstep;

    double actred,actrs,alpha,dirder,epsmac,olmavg,partol,pnorm,prered;
    double prers,ratio,rcond,rnorm,rnormn = 0.0,rnorms,rvar,sstol,tau,taufac;
    double temp,temp1,temp2,tsnorm,wss,wssdel,wsseps;     
    double p0001 = 0.0001;
    double p1 = 0.1;
    double p25 = 0.25;
    double p5  = 0.5;
    double p75 = 0.75;
    double w2[3];
 
    info = 0;
    w2[1] = -1.0;

    /* INITIALIZE NECESSARY VARIABLES */
 
    access = 1;
    dacces(ctx, n,m,np,work,lwork,iwork,liwork,access,&jpvt,&wrk1,&tfjacb,&omega,
        &yt,&u,&qraux,&wrk2,&nnzw,&npp,&partol,&sstol,&taufac,&epsmac,&neta,
        &wss,&wssdel,&wsseps,&rvar,&idf,&tau,&alpha,&int2,
        &olmavg,&rcond,&irank,&actrs,&pnorm,&prers,&rnorms);
   
    rnorm = sqrt(wss);
    didvcv = intdbl = 0;
    lstep = 1;
 
    if (rnorm  ==  0.0) {   /* stop if initial estimate is solution */
        info = 1;
        olmavg = 0.0;
        goto DODMNFin;
    }         
 
    /* MAIN LOOP */

    if (ctx->ORDCyc == 0)
        if (ctx->SILENTFlg < 2)
            printfe(ctx, "\nCycle  Iter    Norm of Errors      Parameter Change    FCall\n");
 
DODMNNxt:
    ctx->ODRNIt += 1;       
    rnorms = rnorm;
 
    /* EVALUATE JACOBIAN */
  
    r = devjac(ctx, n,np,npp,m,betac,beta,ifixb,ifixx,ifixu,x,delta,
                 xplusd,ss,tt,ldtt,neta,fn,sss,fjacb,fjacx,w);
    if (r)
        return(-9);

    /* COMPUTE DDELT = (W*D)**2 * DELTA */

    dwds(ctx, n,m,w,wd,ldwd,delta,ddelt);
    dwds(ctx, n,m,w,wd,ldwd,ddelt,ddelt);
 
    /* SUB LOOP FOR INTERNAL DOUBLING OR COMPUTING NEW STEP WHEN OLD FAILED  */
 
    niter = 0;

DODMN1:    /* COMPUTE STEPS S AND T */

    if (++niter > ctx->ODRMXIt)  
        return(-8);

    niter++;
    dodlm(ctx, n,np,npp,m,f,fjacb,fjacx,w,wd,ldwd,ss,tt,ldtt,ddelt,&alpha,
        &tau,epsmac,sss,work+wrk1-1,work+tfjacb-1,work+omega-1,work+yt-1,
        work+u-1,work+qraux-1,work+wrk2-1,iwork+jpvt-1,s,t,&nlms,&rcond,&irank);

    olmavg += (double)nlms;
 
    /* COMPUTE BETAN = BETAC + S, DELTAN = DELTA + T */
 
    dxpy(ctx, npp,1,betac,s,betan);
    dxpy(ctx, n,m,delta,t,deltan);

    /* ##  COMPUTE NORM OF SCALED STEPS S AND T (TSNORM) */
 
    if (npp >= 1)      
         ddiags(ctx, npp,1,ss,npp,s,1,sss,1);
          
    ddiags(ctx, n,m,tt,ldtt,t,m,sss+npp,m);
    tsnorm = dnrm2(ctx, npp+n*m,sss,1);
 
    /* COMPUTE SCALED PREDICTED REDUCTION */
   
    for (i = 1; i <= n; ++i) {
        sss[i] = ddot(ctx, npp,fjacb+(i-1)*np,1,s,1) +
                 ddot(ctx, m,fjacx+(i-1)*m,1,t+(i-1)*m,1);
    }
    dwds(ctx, n,m,w,wd,ldwd,t,sss+n);
    temp1 = dnrm2(ctx, n+n*m,sss,1) / rnorm;
    temp2 = sqrt(alpha) * tsnorm / rnorm;
    prered = pow(temp1,2.0) + pow(temp2,2.0) / p5;
    dirder = -(pow(temp1,2.0) + pow(temp2,2.0));
 
    /* EVALUATE WEIGHTED EPSILONS AT NEW POINT */
   
    r = devfun(ctx, n,np,m,betan,beta,ifixb,x,y,deltan,xplusd,w,fn);
    if (r) {   
                 /* SET NORM TO INDICATE STEP SHOULD BE REJECTED */
 
        rnormn /= (p1 * p75);
    }
    else {  /* COMPUTE NORM OF NEW WEIGHTED EPSILONS AND WEIGHTED DELTAS */
 
        dcopy(ctx, n,fn,1,sss,1);
        dwds(ctx, n,m,w,wd,ldwd,deltan,sss+n);
        rnormn = dnrm2(ctx, n+n*m,sss,1);
    }          
 
    /* COMPUTE SCALED ACTUAL REDUCTION  */
 
    if (p1*rnormn < rnorm)         
        actred = 1.0 - pow(rnormn/rnorm,2.0);
    else
        actred = -1.0;
 
    /* COMPUTE RATIO OF ACTUAL REDUCTION TO PREDICTED REDUCTION */
 
    if (prered  ==  0.0)     
        ratio = 0.0;
    else
        ratio = actred/prered;

    /* CHECK ON LACK OF REDUCTION IN INTERNAL DOUBLING CASE */
 
    if (intdbl && (ratio < p0001 || rnormn > rnorms)) {    
        tau *= p5;
        alpha /= p5;
        dcopy(ctx, npp,betas,1,betan,1);
        dcopy(ctx, n*m,deltas,1,deltan,1);
        dcopy(ctx, n,fs,1,fn,1);
        actred = actrs;
        prered = prers;
        rnormn = rnorms;
        ratio = p5;
    }        
 
    /* UPDATE STEP BOUND */
 
    intdbl = 0;
    if (ratio < p25) {    
        if (actred >= 0.0)     
            temp = p5;
        else
            temp = p5*dirder/(dirder+p5*actred);

        if (p1*rnormn >= rnorm || temp < p1)     
            temp = p1;

        tau = temp * dmin(ctx, tau,tsnorm/p1);
        alpha /= temp;
    }
    else if (alpha == 0.0) {    
        tau = tsnorm/p5;
    }
    else if (ratio >= p75 && nlms <= 11) {   
 
        /* STEP QUALIFIES FOR INTERNAL DOUBLING
           - UPDATE TAU AND alpha
           - SAVE INFORMATION FOR CURRENT POINT */
 
        intdbl = 1;
 
        tau = tsnorm / p5;
        alpha = alpha * p5;
 
        dcopy(ctx, npp,betan,1,betas,1);
        dcopy(ctx, n*m,deltan,1,deltas,1);
        dcopy(ctx, n,fn,1,fs,1);
        actrs = actred;
        prers = prered;
        rnorms = rnormn;
    }        
 
    /* IF INTERNAL DOUBLING, SKIP CONVERGENCE CHECKS */

    if (intdbl && tau > 0.0) {    
        int2 += 1;
        goto DODMN1;
    }          
 
    /* CHECK ACCEPTANCE */

    if (ratio >= p0001) {    
        dcopy(ctx, n,fn,1,f,1);
        dcopy(ctx, npp,betan,1,betac,1);
        dcopy(ctx, n*m,deltan,1,delta,1);
        rnorm = rnormn;
        if (npp >= 1)       
           ddiags(ctx, npp,1,ss,npp,betac,1,sss,1);
        ddiags(ctx, n,m,tt,ldtt,delta,m,sss+npp,m);
        pnorm = dnrm2(ctx, npp+n*m,sss,1);
        lstep = 1;
    }
    else  
        lstep = 0;
 
    /* TEST CONVERGENCE */
 
    info = 0;
    cnvpar = cnvss = 0;
    if (rnorm == 0.0 || (fabs(actred) <= sstol && prered <= sstol && p5*ratio <= 1.0))
        cnvss = 1;

    if (tau <= partol*pnorm)
        cnvpar = 1;

    if (cnvss)  
        info = 1;
    if (cnvpar)     
        info = 2;
    if (cnvss && cnvpar) 
        info = 3;
 
    /* info about current interation */
 
    if (ctx->SILENTFlg < 2) {
        printfe(ctx, "%4d %5d  %20.13e %20.13e %6d (%d)\n",ctx->ORDCyc,ctx->ODRNIt,rnorm,0.0,ctx->ODRNF,ctx->ODRNFJ);
        fflushe(ctx);
    }
    if (ctx->PMProtFDef) {
        fprintf(ctx->PMProtFd,"\nIteration: %d\n",ctx->ODRNIt);
        fprintf(ctx->PMProtFd,"Current number of function calls: %d (%d)\n\n",ctx->ODRNF,ctx->ODRNFJ);
        prvec(ctx, "Parameter",np,betac);
        prval(ctx, "Reciprocal condition number",rcond);
        prval(ctx, "Norm of errors",rnorm);
    }
    if (info == 0) {        /* check if finished */
        if (lstep) {       
            if (ctx->ODRNIt >= ctx->ODRMXIt)      
                goto DODMNFin;
            else
                goto DODMNNxt;  /* begin next iteration */
        }      
        else                    /* step failed */
            goto DODMN1;
    }         
 
DODMNFin:
    if (ctx->PMProtFDef) {
        fprintf(ctx->PMProtFd,"\nEnd of iterations.\n");
        switch (info) {
            case  0:    fprintf(ctx->PMProtFd,"Exceeded max number of iterations.\n");
                        break;
            case  1:    fprintf(ctx->PMProtFd,"Squared convergence.\n");
                        break;
            case  2:    fprintf(ctx->PMProtFd,"Parameter convergence.\n");
                        break;
            case  3:    fprintf(ctx->PMProtFd,"Squared and parameter convergence.\n");
                        break;
            default:    fprintf(ctx->PMProtFd,"No convergence.\n");
                        break;
        }             
    }
 
    /* COMPUTE UNWEIGHTED EPSILONS AND X+DELTA TO RETURN TO USER */
 
    r = devfun(ctx, n,np,m,betac,beta,ifixb,x,y,delta,xplusd,w2,f);
    if (r)
        return(-9);
 
    /* COMPUTE VARIANCE COVARIANCE MATRIX OF ESTIMATED PARAMETERS
       IN UPPER TRIANGULAR PORTION OF work[tfjacb] */

    if (irank == 0 && info) {       /* if no rank deficiency */
 
        /* EVALUATE JACOBIANS AT FINAL SOLUTION */

        r = devjac(ctx, n,np,npp,m,betac,beta,ifixb,ifixx,ifixu,x,delta,
                     xplusd,ssf,tt,ldtt,neta,fn,sss,fjacb,fjacx,w);
        if (r)      
            return(-9);  

        idf = 0;
        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= npp; ++j) {
                if (fjacb[(i-1)*np+j] != 0.0) {    
                    idf += 1;
                    goto L70;
                }
            }
            for (j = 1; j <= m; ++j) {
                if (fjacx[(i-1)*m+j] != 0.0) {     
                    idf += 1;
                    goto L70;
                }
            }
L70:        ;
        }
        if (ctx->ODRTyp) {   
 
            /* PROBLEM IS ODR --
               SET UP OMEGA AND TFJACB
               (VDTD = fjacx * INV(DT) WHERE DT = (W*D)**2) */
 
            didts(ctx, n,m,w,wd,ldwd,0.0,tt,ldtt,fjacx,work+wrk1-1);
   
            for (i = 1; i <= n; ++i) {
                work[omega-1+i] = sqrt(1.0+ddot(ctx, m,work+wrk1+i-2,1,fjacx+(i-1)*m,1));

                for (j = 1; j <= npp; ++j)  
                    work[tfjacb+(i-1)*np+j-1] = fjacb[(i-1)*np+j]/work[omega-1+i];
            }
        }
        else {  /* PROBLEM IS OLS */ 
 
            dcopy(ctx, n*npp,fjacb,1,work+tfjacb-1,1);
        }        
        dqrdc(ctx, work+tfjacb-1,np,n,npp,work+qraux-1,iwork+jpvt-1,work+wrk2-1,0);
        dpodi(ctx, work+tfjacb-1,np,npp,work+wrk2-1,1);
               
        if (idf > npp) {    
            idf -= npp;
            rvar = rnorm*rnorm/(double)idf;
        }
        else {
            idf = 0;
            rvar = rnorm*rnorm;
        }       
        dscal(ctx, n*npp,rvar,work+tfjacb-1,1);
        dcopy(ctx, npp,work+tfjacb-1,np+1,work+wrk2-1,1);   /* diagonal */

        didvcv = 1;
        if (np > npp) {   
            junfix = npp-1;
            for (j = np - 1; j >= 0; --j) {
                if (ifixb[j+1] == 0)      
                    work[wrk2+j] = 0.0;
                else {
                    temp = work[wrk2+junfix];
                    if (temp >= 0.0)
                        temp = sqrt(temp);
                    else {
                        temp = 0.0;
                        didvcv = 0;
                    }
                    work[wrk2+j] = temp;                     
                    junfix -= 1;
                }
            }
        }
        else {
            for (j = 0; j < np; ++j) {
                temp = work[wrk2+j];
                if (temp >= 0.0)
                    temp = sqrt(temp);
                else {
                    temp = 0.0;
                    didvcv = 0;
                }
                work[wrk2+j] = temp;                 
            }
        }      
        if (didvcv) {    
            ctx->ODRCOVPtr = tfjacb - 1;
            ctx->ODRDf = idf;
        }
    }        
    ctx->ODRSEPtr = wrk2 - 1;

    /* STORE VARIOUS SCALARS IN WORK ARRAYS FOR RETURN TO USER */
 
    olmavg /= (double)ctx->ODRNIt;
 
    /* COMPUTE WEIGHTED EPSILONS AND WEIGHTED DELTAS FOR RETURN TO USER */

    ddiagw(ctx, n,1,w,f,1,sss,1);
    wsseps = ddot(ctx, n,sss,1,sss,1);
    dwds(ctx, n,m,w,wd,ldwd,delta,sss+n);
    wssdel = ddot(ctx, n*m,sss+n,1,sss+n,1);
    wss = wsseps + wssdel;
 
    /* COMPUTE ESTIMATED RESPONSE VARIABLE RETURN TO USER, I.E.,
       EST<Y> = OBS<Y> + EST<EPSILON> */
 
    dxpy(ctx, n,1,y,f,fn);
    access = 0;
    dacces(ctx, n,m,np,work,lwork,iwork,liwork,access,&jpvt,&wrk1,&tfjacb,
        &omega,&yt,&u,&qraux,&wrk2,&nnzw,&npp,&partol,&sstol,&taufac,
        &epsmac,&neta,&wss,&wssdel,&wsseps,&rvar,&idf,&tau,&alpha,&int2,
        &olmavg,&rcond,&irank,&actrs,&pnorm,&prers,&rnorms);

    ctx->ODRRDef = irank;
    ctx->ODRWEps = wsseps;
    ctx->ODRWDel = wssdel;
    ctx->ODRRVar = rvar;
    return(info);
}

/*  ----------------------------------------------------------------------- */    
/*  dodstp()                                                                */
/*                                                                          */
/*  dodstp(n,np,npp,m,f,fjacb,fjacx,w,wd,ldwd,ss,tt,ldtt,ddelt,alpha,       */
/*      epsmac,sss,tfjacb,vdtd,omega,yt,u,qraux,wrk2,jpvt,s,t,phi,          */
/*      irank,rcond)                                                        */
/*                                                                          */
/*  COMPUTE LOCALLY CONSTRAINED STEPS S AND T, AND PHI(ALPHA)               */
/*                                                                          */

void dodstp(TDAContext *ctx, int n,int np,int npp,int m,double *f,double *fjacb,              double *fjacx,double *w,double *wd,int ldwd,double *ss, double *tt,int ldtt,double *ddelt,double alpha,double epsmac, double *sss,double *tfjacb,double *vdtd,double *omega,double *yt, double *u,double *qraux,double *wrk2,int *jpvt,double *s,double *t, double *phi,int *irank,double *rcond)
{
    int i,iimax,inf,ipvt,j,kp,elim;
    double co,si,temp;     
    double dum[2];
 
    /* COMPUTE LOOP PARAMETERS WHICH DEPEND ON WEIGHT STRUCTURE */
 
    /* SET UP JPVT IF ALPHA = 0 */ 
 
    if (alpha == 0.0) {    
        kp = npp;
        for (i = 1; i <= npp; ++i)  
            jpvt[i] = i;
    }
    else {
        if (npp >= 1)       
            kp = npp - *irank;
        else 
            kp = npp;
    }        

    /* SET UP OMEGA AND TFJACB
       (VDTD = fjacx * INV(DT) WHERE DT = (W*D)**2 + alpha*TT**2) */
  
    didts(ctx, n,m,w,wd,ldwd,alpha,tt,ldtt,fjacx,vdtd);
   
    for (i = 1; i <= n; ++i)  
        omega[i] = sqrt(1.0 + ddot(ctx, m,vdtd+(i-1)*m,1,fjacx+(i-1)*m,1));
     
    for (j = 1; j <= kp; ++j) {
        for (i = 1; i <= n; ++i)  
            tfjacb[(i-1)*np+j] = fjacb[(i-1)*np+jpvt[j]] / omega[i];
    }

    /* SET UP VDTD AND YT
       (VDTD = ddelt * INV(DT) WHERE DT = (W*D)**2 + alpha*TT**2) */
  
    didts(ctx, n,m,w,wd,ldwd,alpha,tt,ldtt,ddelt,vdtd);
  
    for (i = 1; i <= n; ++i) {
        vdtd[(i-1)*m+1] = ddot(ctx, m,fjacx+(i-1)*m,1,vdtd+(i-1)*m,1);
        yt[i] = -(f[i] - vdtd[(i-1)*m+1]) / omega[i];
    }

    /* COMPUTE S */ 
 
    /* DO QR FACTORIZATION (WITH COLUMN PIVOTING OF TRJACB IF ALPHA = 0) */
 
    if (alpha == 0.0) {    
        ipvt = 1;
        for (i = 1; i <= npp; ++i)  
            jpvt[i] = 0;
    }
    else   
        ipvt = 0;

    dqrdc(ctx, tfjacb,np,n,kp,qraux,jpvt,wrk2,ipvt);
 
    /* GET TR(Q)*YT */
 
    dqrsl(ctx, tfjacb,kp,n,kp,qraux,yt,dum,yt,dum,dum,dum,1000,&inf);
 
    /* ELIMINATE ALPHA PART USING GIVENS ROTATIONS */
 
    if (alpha != 0.0) {    
        dzero(ctx, npp,1,s,1);              
        for (i = 1; i <= kp; ++i) {
            dzero(ctx, kp,1,wrk2,1);          
            if (ss[1] > 0.0) {     
               wrk2[i] = sqrt(alpha) * ss[jpvt[i]];
            }
            else
               wrk2[i] = sqrt(alpha) * fabs(ss[1]);
    
            for (j = i; j <= kp; ++j) {
                drotg(ctx, &tfjacb[(j-1)*np+j],&wrk2[j],&co,&si);
                if (kp - j >= 1) {    
                    drot(ctx, kp-j,tfjacb + (j-1)*np+j,np,wrk2+j,1,co,si);
                }     
                temp = co*yt[j] + si*s[jpvt[i]];
                s[jpvt[i]] = -si*yt[j] + co*s[jpvt[i]];
                yt[j] = temp;
            }
        }
    }       
 
    /* COMPUTE SOLUTION - ELIMINATE VARIABLES IF NECESSARY  */
 
    if (npp >= 1) {    
        if (alpha == 0.0) {    
            kp = npp;
 
            /* ESTIMATE RCOND - U WILL CONTAIN APPROX NULL VECTOR */
 
L100:       dtrco(ctx, tfjacb,np,kp,rcond,u,1);
            if (*rcond <= epsmac) {   
                elim = 1;
                iimax = idamax(ctx, kp,u,1);
 
                /* IMAX IS THE COLUMN TO REMOVE - USE DCHEX AND FIX JPVT */
  
                if (iimax != kp) {                  
                    dchex(ctx, tfjacb,np,kp,iimax,kp,yt,1,1,qraux,wrk2,2);
                    j = jpvt[iimax];
                    for (i = iimax; i < kp; ++i) {
                        jpvt[i] = jpvt[i+1];
                    }
                    jpvt[kp] = j;
                }      
                kp--;     
            }
            else 
                elim = 0;

            if (elim && kp >= 1)     
                goto L100;
            else
                *irank = npp-kp;
        }     
         
        /* BACKSOLVE AND UNSCRAMBLE */
 
        for (i = kp + 1; i <= npp; ++i)  
            yt[i] = 0.0;
             
        if (kp >= 1)      
            dtrsl(ctx, tfjacb,kp,kp,yt,1,&inf);     
               
        for (i = 1; i <= npp; ++i)  
            s[jpvt[i]] = yt[i];
    }     

    /* COMPUTE T */
   
    for (i = 1; i <= n; ++i) {
        temp = f[i] + ddot(ctx, npp,fjacb+(i-1)*np,1,s,1);   
        u[i] = (temp - vdtd[(i-1)*m+1]) / (omega[i] * omega[i]);
    }
    for (j = 1; j <= m; ++j) {
        for (i = 1; i <= n; ++i)  
            t[(i-1)*m+j] = -(fjacx[(i-1)*m+j] * u[i] + ddelt[(i-1)*m+j]);
    }
          
    /* (T = T * INV(DT) WHERE DT = (W*D)**2 + alpha*TT**2)  */

    didts(ctx, n,m,w,wd,ldwd,alpha,tt,ldtt,t,t);
        
    /* COMPUTE PHI(alpha) FROM SCALED S AND T */
            
    if (npp >= 1)      
        ddiags(ctx, npp,1,ss,npp,s,1,sss,1);

    ddiags(ctx, n,m,tt,ldtt,t,m,sss+npp,m);
    *phi = dnrm2(ctx, npp+n*m,sss,1);
}             

/*  ----------------------------------------------------------------------- */    
/*  int dpack()                                                             */
/*                                                                          */
/*  dpack(n2,*n1,*v1,*v2,*ifix)                                             */
/*                                                                          */
/*  SELECT THE UNFIXED ELEMENTS OF V2 AND RETURN THEM IN V1                 */
/*  Return number of parameters.                                            */

int dpack(TDAContext *ctx, int n2,double *v1,double *v2,char *ifix)
{
    (void)ctx;        /* unused: the signature is shared */
    int i,n;
 
    n = 0;
    for (i = 1; i <= n2; ++i) {
        if (ifix[i]) {    
            n += 1;  
            v1[n] = v2[i];
        }       
    }
    return(n);
}

/*  ----------------------------------------------------------------------- */    
/*  dsclb()                                                                 */
/*                                                                          */
/*  dsclb(np,beta,ssf)                                                      */
/*                                                                          */
/*  COMPUTE APPROPRIATE SCALE VALUES FOR BETA'S ACCORDING TO                */
/*  THE ALGORITHM GIVEN IN THE PROLOGUES FOR DODR AND DODRC                 */
/*                                                                          */

void dsclb(TDAContext *ctx, int np,double *beta,double *ssf)
{
    double bmax,bmin,tmp1,tmp2;
    int k,bigdif;
 
    bmax = fabs(beta[1]);
    for (k = 2; k <= np; ++k) {
        bmax = dmax(ctx, bmax,fabs(beta[k]));
    }
    if (bmax == 0.0) {     /* ALL INPUT VALUES OF BETA ARE ZERO */
 
        for (k = 1; k <= np; ++k)  
            ssf[k] = 1.0;
    }
    else {  /* SOME OF THE INPUT VALUES ARE NONZERO */
 
        bmin = bmax;
        for (k = 1; k <= np; ++k) {
            if (beta[k] != 0.0) {    
               bmin = dmin(ctx, bmin,fabs(beta[k]));
            }      
        }
        bigdif = 0;
        tmp1 = rlog(ctx, bmax) / log(10.0);
        tmp2 = rlog(ctx, bmin) / log(10.0);

        if (tmp1 - tmp2 >= 1.0)
            bigdif = 1;

        for (k = 1; k <= np; ++k) {
            if (beta[k] == 0.0)      
                ssf[k] =  10.0/bmin;
            else {
                if (bigdif)      
                    ssf[k] = 1.0/fabs(beta[k]);
                else
                    ssf[k] = 1.0/bmax;
            }     
        }
    }        
}

/*  ----------------------------------------------------------------------- */    
/*  dscld()                                                                 */
/*                                                                          */
/*  dscld(n,m,x,tt)                                                         */
/*                                                                          */
/*  COMPUTE APPROPRIATE SCALE VALUES FOR DELTA'S ACCORDING TO               */
/*  THE ALGORITHM GIVEN IN THE PROLOGUES FOR DODR AND DODRC                 */
/*                                                                          */

void dscld(TDAContext *ctx, int n,int m,double *x,double *tt)
{
    double xmax,xmin,tmp1,tmp2;
    int i,j,bigdif;

    for (j = 1; j <= m; ++j) {
        xmax = fabs(x[j]);
        for (i = 2; i <= n; ++i) {
            xmax = dmax(ctx, xmax,fabs(x[(i-1)*m+j]));
        }
        if (xmax == 0.0) {               /* X = 0 */
            for (i = 1; i <= n; ++i)  
                tt[(i-1)*m+j] = 1.0;
        }
        else {
            xmin = xmax;
            for (i = 1; i <= n; ++i) {
                if (x[(i-1)*m+j] != 0.0) {     
                    xmin = dmin(ctx, xmin,fabs(x[(i-1)*m+j]));
                }      
            }
            bigdif = 0;
            tmp1 = rlog(ctx, xmax) / log(10.0);
            tmp2 = rlog(ctx, xmin) / log(10.0);
            if (tmp1 - tmp2 >= 1.0)
                bigdif = 1;

            for (i = 1; i <= n; ++i) {
                if (x[(i-1)*m+j] != 0.0) {    
                    if (bigdif)     
                        tt[(i-1)*m+j] = 1.0/fabs(x[(i-1)*m+j]);
                    else
                        tt[(i-1)*m+j] = 1.0/xmax;
                }
                else {
                  tt[(i-1)*m+j] = 10.0/xmin;
                }     
            }
        }     
    }
}          

/*  ----------------------------------------------------------------------- */    
/*  dsetn()                                                                 */
/*                                                                          */
/*  dsetn(n,m,x,*nrow)                                                      */
/*                                                                          */
/*  SELECT THE ROW AT WHICH THE DERIVATIVE WILL BE CHECKED                  */
/*                                                                          */

void dsetn(TDAContext *ctx, int n,int m,double *x,int *nrow)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j;
 
    if ((*nrow >= 1) && (*nrow <= n))        
        return;

    /* SELECT FIRST ROW OF INDEPENDENT VARIABLES WHICH CONTAINS NO ZEROS
       IF THERE IS ONE, OTHERWISE FIRST ROW IS USED. */
 
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= m; ++j) {
            if (x[(i-1)*m+j] == 0.0)
                goto L20;
        }
        *nrow = i;
        return;
L20: ;
    }
    *nrow = 1;
}

/*  ----------------------------------------------------------------------- */    
/*  dunpac()                                                                */
/*                                                                          */
/*  dunpac(n2,v1,v2,ifix)                                                   */
/*                                                                          */
/*  COPY THE ELEMENTS OF V1 INTO THE LOCATIONS OF V2 WHICH ARE UNFIXED.     */
/*                                                                          */

void dunpac(TDAContext *ctx, int n2,double *v1,double *v2,char *ifix)
{
    (void)ctx;        /* unused: the signature is shared */
    int i,n1;
 
    n1 = 0;
    for (i = 1; i <= n2; ++i) {
        if (ifix[i]) {    
            n1 += 1;
            v2[i] = v1[n1];
        }
    }
}

/*  ----------------------------------------------------------------------- */    
/*  dwds()                                                                  */
/*                                                                          */
/*  dwds(n,m,w,wd,ldwd,t,wdt)                                               */
/*                                                                          */
/*           SCALE MATRIX T USING w*d, I.E., COMPUTE                        */
/*           wdt = w*d*t                                                    */
/*           WHERE w AND d ARE DEFINED BY EQ.2 OF THE PROLOGUES FOR         */
/*           DODR AND DODRC.                                                */
/*                                                                          */
/*  Note: we assume that column dimension of t and wdt, and also of w if    */
/*  w[1] >= 0, is always m.                                                 */  

void dwds(TDAContext *ctx, int n,int m,double *w,double *wd,int ldwd,double *t,double *wdt)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j;
    double temp;     

    if (n == 0 || m == 0)         
        return;

    if (w[1] >= 0.0) {   
        if (wd[1] > 0.0) {    
            if (ldwd >= n) {    
                for (j = 1; j <= m; ++j) {
                    for (i = 1; i <= n; ++i) {
                        wdt[(i-1)*m+j] = w[i]*wd[(i-1)*m+j]*t[(i-1)*m+j];
                    }
                }
            }
            else {
                for (j = 1; j <= m; ++j) {
                    for (i = 1; i <= n; ++i) {
                        wdt[(i-1)*m+j] = w[i]*wd[j]*t[(i-1)*m+j];
                    }
                }
            }     
        }
        else {
            for (j = 1; j <= m; ++j) {
                for (i = 1; i <= n; ++i) {
                    wdt[(i-1)*m+j] = w[i]*fabs(wd[1])*t[(i-1)*m+j];
                }
            }
        }       
    }
    else {
        if (wd[1] > 0.0) {   
            if (ldwd >= n) {    
                for (j = 1; j <= m; ++j) {
                    for (i = 1; i <= n; ++i) {
                        wdt[(i-1)*m+j] = wd[(i-1)*m+j]*t[(i-1)*m+j];
                    }
                }
            }
            else {
                for (j = 1; j <= m; ++j) {
                    temp = wd[j];
                    for (i = 1; i <= n; ++i) {
                        wdt[(i-1)*m+j] = temp*t[(i-1)*m+j];
                    }
                }
            }      
        }
        else {
            temp = fabs(wd[1]);
            for (j = 1; j <= m; ++j) {
                for (i = 1; i <= n; ++i) {
                    wdt[(i-1)*m+j] = temp*t[(i-1)*m+j];
                }
            }
        }       
    }        
}

/*  ----------------------------------------------------------------------- */    
/*  dwinf()                                                                 */
/*                                                                          */
/*  dwinf(n,m,np,*deltai,*epsi,*wssi,*wssdei,*wssepi,*rvari,                */
/*        *partli,*sstoli,*taufci,*epsmai,*olmavi,                          */
/*        *fjacbi,*fjacxi,*xplusi,*betaci,*betasi,*betani,*deltsi,          */
/*        *deltni,*ddelti,*fsi,*fni,*si,*sssi,*ssi,*ssfi,*ti,*tti,*taui,    */
/*        *alphai,*vcvi,*omegai,*yti,*ui,*qrauxi,*wrk1i,*sei,*rcondi,       */
/*        *etai,*actrsi,*pnormi,*prersi,*rnorsi,*lwkmn)                     */
/*                                                                          */
/*  SET STORAGE LOCATIONS WITHIN DOUBLE PRECISION WORK SPACE                */
/*                                                                          */

void dwinf(TDAContext *ctx, int n,int m,int np,int *deltai,int *epsi,int *wssi,int *wssdei, int *wssepi,int *rvari,int *partli,int *sstoli,int *taufci,int *epsmai, int *olmavi,int *fjacbi,int *fjacxi,int *xplusi,int *betaci,int *betasi, int *betani,int *deltsi,int *deltni,int *ddelti,int *fsi,int *fni,int *si, int *sssi,int *ssi,int *ssfi,int *ti,int *tti,int *taui,int *alphai, int *vcvi,int *omegai,int *yti,int *ui,int *qrauxi,int *wrk1i,int *sei, int *rcondi,int *etai,int *actrsi,int *pnormi,int *prersi,int *rnorsi, int *lwkmn)
{
    (void)ctx;        /* unused: the signature is shared */
    *deltai =          1;
    *epsi   = *deltai + n*m;
    *wssi   = *epsi   + n;
    *wssdei = *wssi  + 1;
    *wssepi = *wssdei + 1;
    *rvari  = *wssepi + 1;
    *partli = *rvari  + 1;
    *sstoli = *partli + 1;
    *taufci = *sstoli + 1;
    *epsmai = *taufci + 1;
    *olmavi = *epsmai + 1;
    *fjacbi = *olmavi + 1;
    *fjacxi = *fjacbi + n*np;
    *xplusi = *fjacxi + n*m;
    *betaci = *xplusi + n*m;
    *betasi = *betaci + np;
    *betani = *betasi + np;
    *deltsi = *betani + np;
    *deltni = *deltsi + n*m;
    *ddelti = *deltni + n*m;
    *fsi    = *ddelti + n*m;
    *fni    = *fsi    + n;
    *si     = *fni    + n;
    *sssi   = *si     + np;
    *ssi    = *sssi  + n*m + n;
    *ssfi   = *ssi    + np;
    *ti     = *ssfi   + np;
    *tti    = *ti     + n*m;
    *taui   = *tti    + n*m;
    *alphai = *taui   + 1;
    *vcvi   = *alphai + 1;
    *omegai = *vcvi   + n*np;
    *yti    = *omegai + n;
    *ui     = *yti    + n;
    *qrauxi = *ui     + n;
    *wrk1i  = *qrauxi + np;
    *sei    = *wrk1i  + n*m;
    *rcondi = *sei    + np;
    *etai   = *rcondi + 1;
    *actrsi = *etai   + 1;
    *pnormi = *actrsi + 1;
    *prersi = *pnormi + 1;
    *rnorsi = *prersi + 1;
    *lwkmn  = *rnorsi;
}

/*  ----------------------------------------------------------------------- */    
/*  diwinf()                                                                */
/*                                                                          */
/*  diwinf(m,np,jpvti,nnzwi,nppi,idfi,nrowi,netai,int2i,iranki,             */
/*      ldtti,liwkmn)                                                       */
/*                                                                          */
/*  SET STORAGE LOCATIONS WITHIN INTEGER WORK SPACE                         */
/*                                                                          */

void diwinf(TDAContext *ctx, int m,int np,int *jpvti,int *nnzwi,int *nppi,int *idfi, int *nrowi,int *netai,int *int2i,int *iranki,int *ldtti,int *liwkmn)
{
    (void)ctx; (void)m;        /* unused: the signature is shared */
    *jpvti  = 1;            
    *nnzwi  = *jpvti  + np;
    *nppi   = *nnzwi  + 1;
    *idfi   = *nppi   + 1;
    *nrowi  = *idfi   + 1;
    *netai  = *nrowi  + 1;
    *int2i  = *netai  + 1;
    *iranki = *int2i  + 1;
    *ldtti  = *iranki + 1;
    *liwkmn = *ldtti;
}           

/*  ----------------------------------------------------------------------- */  
/*  odr_fun()   Function evaluation.                                        */
/*  ##                                                                      */
/*  odr_fun(n,np,m,beta,xplusd,f,cdf)                                       */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int odr_fun(TDAContext *ctx, int n,int np,int m,double *beta,double *xplusd,double *f,int cdf)
{
    register int i,j,k;
    double tmp;

    if (ctx->FNFlg) {        /* user-defined function */

        for (j = 0; j < np; ++j)        /* set values for evaluation */
            ctx->FNArgVal[j] = beta[j + 1];

        for (i = 0; i < n; ++i) {

            /* if ODRTyp = 1 we also need to exchange the data for the
               X variables in TDA's data matrix! We use the AVVAL option
               in get_fival(ctx) and v_eval1(ctx), resp. */

            if (ctx->ODRTyp) {
                for (j = 1; j <= m; ++j) {
                    k = ctx->PMVIdx[j];
                    ctx->AVVAL[k] = xplusd[i * m + j];
                }
            }
            k = get_fival(ctx, i,&tmp,0,ctx->ODRTyp);
            if (k) {
                printf1(ctx, "Error (%d): can't evaluate function for case %d.\n",k,i + 1);
                prn_emsg2(ctx, k);
                return(-1);
            }
            f[i * cdf + 1] = tmp;                            
        }
    }
    else {              /* linear predictor */

        for (i = 0; i < n; ++i) {
            tmp = 0.0;
            j = 0;
            if (ctx->PMNI == 0)
                tmp = beta[++j];
     
            for (k = 1; k <= m; ++k) {
                tmp += xplusd[i * m + k] * beta[++j];
            }
            f[i * cdf + 1] = tmp;                            
        }
    }
    return(0);
}

/*  ----------------------------------------------------------------------- */  
/*  odr_jac()   Evaluation of derivatives.                                  */
/*                                                                          */  
/*  odr_jac(n,np,m,beta,xplusd,fjacb,fjacx)                                 */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int odr_jac(TDAContext *ctx, int n,int np,int m,double *beta,double *xplusd,double *fjacb, double *fjacx)
{
    register int i,j,k;
    double tmp;

    if (ctx->FNFlg) {        /* user-defined function */

        if (ctx->ODRTyp) {
            printf1(ctx, "Error in odr_jac!\n");
            exit(0);
        }
        for (j = 0; j < np; ++j)        /* set values for evaluation */
            ctx->FNArgVal[j] = beta[j + 1];

        for (i = 0; i < n; ++i) {

            k = get_fival(ctx, i,&tmp,1,0);
            if (k) {
                printf1(ctx, "Error (%d): can't evaluate derivatives for case %d.\n",k,i + 1);
                prn_emsg2(ctx, k);
                return(-1);
            }
            for (j = 1; j <= np; ++j)  
                fjacb[i * np + j] = ctx->FNGrad[0][j - 1];     
        }
    }
    else {              /* linear predictor */

        for (i = 0; i < n; ++i) {
            j = 0;
            if (ctx->PMNI == 0) {
                fjacb[i * np + 1] = 1.0;                
                j = 1;
            }
            k = 1;
            while (++j <= np) {
                fjacb[i * np + j] = xplusd[i * m + k];
                k++;
            }
        }          
        if (ctx->ODRTyp) {
            for (i = 0; i < n; ++i) {
                if (ctx->PMNI)
                    k = 1;
                else
                    k = 2;
                for (j = 1; j <= m; ++j) {
                    fjacx[i * m + j] = beta[k];
                    k++;
                }
            }
        }          
    }
    return(0);
}

