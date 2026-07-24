/****************************************************************************/
/*  t_glm                                                                   */
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
#include "t_gmin.h"   
#include "t_lsei.h"   
#include "t_gf.h"   
#include "t_var.h"   
#include "t_con.h"   
#include "t_cdf.h"   
#include "tda_context.h"

/*  functions in t_glm.c */

int glm(TDAContext *ctx);
int prn_glmd(TDAContext *ctx, int d);
void glm_prot(TDAContext *ctx, int nx);
int glm_it(TDAContext *ctx, int nx,int intflg,int ne,int ni);
int get_mue(TDAContext *ctx, double eta,double *mue);
int rtbis(TDAContext *ctx, double x1,double x2,double *x,double eta,double *dx);
int glm_lfd(TDAContext *ctx, double mue,double *d);
int glm_dev(TDAContext *ctx, int nx,int intflg,int df);
void glm_dtda(TDAContext *ctx, int intflg,int nx);

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */


/* ------------------------------------------------------------------------ */
/*  glm()   Generalized linear models.                                      */  
/*                                                                          */
/*          glm(                                                            */
/*              d =...,     1 = normal                                      */
/*                          2 = binomial                                    */
/*                          3 = poisson                                     */
/*                          4 = gamma                                       */
/*                          5 = inverse gaussian                            */
/*              link=...,   1 = identity                                    */
/*                          2 = log                                         */
/*                          3 = logit                                       */
/*                          4 = reciprocal                                  */
/*                          5 = probit                                      */
/*                          6 = comp. log-log                               */
/*                          7 = square root                                 */
/*                          8 = quadratic inverse                           */
/*              v=...,      varlist: Y,X1,X2,...                            */
/*              yw=...,     variable for proportions                        */
/*              ni=...,     1 if without intercept                          */
/*              lsecon=     equality constraints                            */
/*              lsicon=     inequality constraints                          */
/*              mxit=...,   max number of iterations, def. 20               */
/*              tolsp=...,  tolerance for convergence, def. 1e-6            */
/*              xp=...,     starting values                                 */
/*              dsv=...,    starting values                                 */
/*              tfmt=...,   print format for parameters, def. 10.4          */
/*              ppar=...,   print parameter                                 */
/*              pcov=...,   print covariance matrix                         */
/*              mfmt=...,   print format pcov                               */
/*              pres=...,   print residuals                                 */
/*              fmt=...,    print format pres                               */
/*              prot=...,   protocol file                                   */
/*              pfmt=...,   print format for protocol file, def. -19.11     */
/*              ab=...,     domain for link function, def. (0,1)            */
/*              tolf=...,   tolerance for inverse link function, 10-12      */
/*                                                                          */
/*          ) = link_function;                                              */
/*                                                                          */
/*          Return 0 if OK, otherwise -1.                                   */

int glm(TDAContext *ctx)
{
    register int i;
    int err,r,nv,nx,nx1,nw,mw,intflg,df,ncon,ne,ni;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Generalized linear models. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLSP = 1.e-6;
    ctx->TOLF  = 1.e-12;

    if (parm(ctx, ctx->CmdBuf + 3,9,0))   /* get parameters */
        goto GLMFin;

    if (ctx->PMNW != 1) {
        printf1(ctx, "Error: glm command requires cross-section data (nw=1).\n");
        goto GLMFin;
    }
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    newline(ctx);
    if (prn_glmd(ctx, (int)ctx->PMD))
        goto GLMFin;

    nv = check_nvar(ctx, 0);         /* check variables */
    if (nv == 0) {              /* PMNV is number of variables */       
        printf1(ctx, "Error: no variables.\n");
        goto GLMFin;        
    }
    prn_nwvar(ctx, 0);               /* print variables */

    nx = nx1 = nv - 1;
    if (ctx->PMNI) {
        printf1(ctx, "Model without intercept.\n");
        ctx->PMNI = 1;
        intflg = 0;
    }
    else {
        ctx->PMNI = 0;
        intflg = 1;
        nx1++;
        if (nx == 0)  
            printf1(ctx, "Model without independent variables.\n");
    }
    if (nx1 == 0) {
        printf1(ctx, "Error: no model parameters.\n");
        goto GLMFin;
    }
    if (ctx->PMYWVar >= 0) {
        printf1(ctx, "Variable used to define proportions: %s",ctx->VName[ctx->PMYWVar]);
        if (ctx->PMD != 2) {
            printf1(ctx, " (will be ignored)");
            ctx->PMYWVar = -1;
        }
        printf1(ctx, "\n\n");
    }   

    /* allocate memory for derivatives of link function */

    if (fnd_alloc(ctx, 1,1,ctx->FNArgN,ctx->FNPN,0))
        goto GLMFin;

    /* data matrix: AcW, parameter vector: AcX */

    ne = ni = ncon = 0;         /* number of constraints */
    nw = nx1 + 1;
    r = mw = ctx->NOC + ctx->NCONSTR;     
    if (r < nw)
        r = nw;

    if (alloc_acw(ctx, r * nw + 1))
        goto GLMFin;

    if (alloc_acx(ctx, nx1 + 1))
        goto GLMFin;

    if (alloc_acy(ctx, nx1 + 1))
        goto GLMFin;

    if (ctx->MxItFlg == 0)
        ctx->MxIter = 20;

    printf1(ctx, "Estimation with iteratively re-weighted least squares.\n");
    if (ctx->WIVar >= 0)  
        printf1(ctx, "Using weights defined by: %s\n",ctx->VName[ctx->WIVar]);

    printf1(ctx, "Number of model parameters: %d\n",nx1);
         
    if (ctx->NCONSTR > 0) {
        newline(ctx);
        if (p_con(ctx, ctx->CmdBuf,nx,ctx->PMNI,mw,nw,ctx->AcW,&ne,&ni))      /* get constraints */
            goto GLMFin;
        ncon = ne + ni;
        printf1(ctx, "Equality constraints: %d\n",ne);        
        printf1(ctx, "Inequality constraints: %d\n\n",ni);          

        /* need some additional memory */

        if (alloc_acu(ctx, ne * nw + 1))
            goto GLMFin;

        if (alloc_acv(ctx, ne * nw + 1))
            goto GLMFin;

        r = ne * nw;
        for (i = 1; i <= r; ++i)
            ctx->AcU[i] = ctx->AcW[i];

        r = ni * nw;
        df = (ne + ctx->NOC) * nw;
        for (i = 1; i <= r; ++i)
            ctx->AcV[i] = ctx->AcW[df + i];
    }
    if (ctx->NOC + ne < nx1 + 1) {
        printf1(ctx, "Error: need at least %d cases.\n",nx1 + 1);
        goto GLMFin;
    }
    printf1(ctx, "Maximum number of iterations: %d\n",ctx->MxIter);
    printf1(ctx, "Tolerance for scaled parameter change: %g\n",ctx->TOLSP);

    if (get_dsv(ctx, nx1,ctx->AcX,0,ctx->ParLB,ctx->ParUB,1))   /* try to get starting values */
        goto GLMFin;                        /* ParLB,ParUB not used */

    glm_prot(ctx, nx1);      /* init protocol file */
    newline(ctx);

    r = glm_it(ctx, nx,intflg,ne,ni);    /* perform iterations */
    if (r < 0) {
        if (r == -2)  
            printf1(ctx, "Convergence not reached.\nExceeded maximum number of iterations.\n");
        goto GLMFin;
    }
    printf1(ctx, "Convergence reached in %d iterations.\n",ctx->GLMIT);
    printf1(ctx, "Final scaled parameter change: %g\n",ctx->GLMPC);

    df = ctx->NOC - ctx->GLMRK;
    if (df < 0)
        df = 0;

    if (glm_dev(ctx, nx,intflg,df))      /* calculate and print deviance etc */
        goto GLMFin;

    if (r > 0)      /* no covariance matrix calculated */
        df = 0;
    else {
        for (i = 1; i <= nx1; ++i)  
            ctx->AcY[i] = ctx->AcW[(i - 1) * nw + i];
    }  
    prn1_coeff(ctx, nx1,ctx->AcX,ctx->AcY,ctx->PMNI,df,ctx->PMVIdx,1);

    newline(ctx);
#ifdef TDA_R_PACKAGE
    if (r == 0)
        export_prn(ctx, "glm.vcov", nx1,nw,nx1,ctx->AcW);
#endif
    if (r == 0 && ctx->PMCovFDef) {                          /* write cov matrix */
        prn_data(ctx, nx1,nw,nx1,ctx->AcW,ctx->PMCovFd,ctx->PMMFmtS);
        p_wmsg(ctx, 2,ctx->PMCovFName,1);
    }
    err = 0;

GLMFin:
    fnd_alloc(ctx, 0,0,0,0,0);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_glmd(d)     Type type of distribution and set GLMD and GLMLF        */

int prn_glmd(TDAContext *ctx, int d)
{
    int nl;

    printf1(ctx, "Distribution: ");
    switch (d) {
        case 2: printf1(ctx, "binomial.\n");
                nl = 3;
                break;
        case 3: printf1(ctx, "Poisson.\n");
                nl = 2;
                break;
        case 4: printf1(ctx, "gamma.\n");
                nl = 4;
                break;
        case 5: printf1(ctx, "inverse Gaussian.\n");
                nl = 8;
                break;
        default: printf1(ctx, "normal.\n");
                d = 1;
                nl = 1;
                break;
    }
    ctx->GLMD = d;

    printf1(ctx, "Link function: ");

    if (ctx->FNFlg) {        /* user-defined link function */
        ctx->GLMLF = 0;
        printf1(ctx, "user-defined.\n");
        if (ctx->FNArgN != 1) {
            printf1(ctx, "Error: link function must have exactly one argument.\n");
            return(-1);  
        }
        if (ctx->FVFlg) {
            printf1(ctx, "Error: link function must not refer to data matrix variables.\n");
            return(-1);   
        }
        if (ctx->PMXYFlg == 0) {
            ctx->PMX = 0.0;
            ctx->PMY = 1.0;
        }
        if (ctx->PMX >= ctx->PMY) {
            printf1(ctx, "Error in definition of range for mue.\n");
            return(-1);  
        }
        printf1(ctx, "Domain of link function: (%g,%g)\n",ctx->PMX,ctx->PMY);
    }
    else {
        if (ctx->PMLINK >= 1 && ctx->PMLINK <= 8)
            ctx->GLMLF = ctx->PMLINK;
        else
            ctx->GLMLF = nl;

        switch (ctx->GLMLF) {
            case  1:    printf1(ctx, "identity.\n"); break;
            case  2:    printf1(ctx, "log.\n"); break;
            case  3:    printf1(ctx, "logit.\n"); break;
            case  4:    printf1(ctx, "reciprocal.\n"); break;
            case  5:    printf1(ctx, "probit.\n"); break;
            case  6:    printf1(ctx, "complementary log-log.\n"); break;
            case  7:    printf1(ctx, "square root.\n"); break;
            case  8:    printf1(ctx, "quadratic inverse.\n"); break;
        }
    }
    newline(ctx);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  glm_prot()  protocol file.                                              */

void glm_prot(TDAContext *ctx, int nx)
{
    register int i;

    if (ctx->PMProtFDef == 0)
        return;

    printf1(ctx, "Protocol will be written to: %s\n",ctx->PMProtFName);

    fprintf(ctx->PMProtFd,"Generalized linear model.\n");
    fprintf(ctx->PMProtFd,"Number of model parameters: %d\n",nx);
    fprintf(ctx->PMProtFd,"\nStarting values.\n");
    for (i = 1; i <= nx; ++i)
        rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,ctx->AcX[i]);
    fprintf(ctx->PMProtFd,"\n\n");
}

/* ------------------------------------------------------------------------ */
/*  glm_it(nx,intflg,ne,ni) Iterated weighted regression for glm.           */
/*                          nx = number of x variables, intflg = 1 if       */
/*                          with intercept. ne and ni are the numbers of    */
/*                          equality and inequality constraints, resp.      */
/*                                                                          */
/*                          AcX contains current parameter vector.          */
/*                          Function requires AcY and AcW as working        */
/*                          space.                                          */
/*                                                                          */
/*  If PMResFDef write data and estimated Y values into output file, and    */  
/*  optionally an TDA description file.                                     */
/*                                                                          */
/*  Return:  0  if successful.                                              */
/*          -1  if error                                                    */
/*          -2  if max number of iterations reached.                        */
/*           1  if covariance matrix not calculated.                        */

int glm_it(TDAContext *ctx, int nx,int intflg,int ne,int ni)
{
    register int i,j,k,iw;
    int iter,err,nw,nx1,ranka,ranke,cov;
    double y,yn,eta,tmp,mue,d,w,rnorme,rnorml,wt,wsum;

    rnorme = rnorml = 0.0;
    ranka = ranke = ctx->GLMIT = ctx->GLMRK = 0;
    nx1 = nx;
    if (intflg)
        nx1++;
    nw = nx1 + 1;
    cov = 0;
    yn = wt = 1.0;
               
    if (ctx->SILENTFlg < 2)
        printfe(ctx, "\n  Iter    Function Value\n");

    for (iter = 0; iter <= ctx->MxIter; ++iter) {

        for (j = 1; j <= nx1; ++j)      /* save current parameters */
            ctx->AcY[j] = ctx->AcX[j];

        if (ctx->PMProtFDef)
            fprintf(ctx->PMProtFd,"\nIteration %d\n",iter);
 
        /* set up matrix for weighted regression */

        iw = ne;
        wsum = 0.0;

        for (i = 0; i < ctx->NOC; ++i) {

            y = get_data(ctx, ctx->PMVIdx[0],i);

            eta = 0.0;                  /* eta = xi * beta */
            j = 1;
            if (intflg)
                eta += ctx->AcX[j++];
            for (k = 1; k <= nx; ++k) {
                tmp = get_data(ctx, ctx->PMVIdx[k],i);
                eta += tmp * ctx->AcX[j++];
            }
            if (get_mue(ctx, eta,&mue))      /* g(mue) = eta */
                return(-1);
   
            if (ctx->PMYWVar >= 0 && ctx->PMD == 2) {
                yn = get_data(ctx, ctx->PMYWVar,i);
                if (y < 0.0 || yn < 1.0 || y > yn) {
                    printf1(ctx, "Error in case %d: y = %g, count = %g\n",i+1,y,yn);
                    return(-1);
                }
                mue *= yn;
            }

            /* ### evaluate link function with mue, and get first derivative in d */
      
            /* tda_out("iter=%d i=%d eta=%g mue=%g yn=%g\n",iter,i,eta,mue,yn); */
   
            tmp = mue / yn;
            if ((err = glm_lfd(ctx, tmp,&d))) {
                printf1(ctx, "Error: cannot evaluate link function, or derivative, with mue = %g (case %d)\n",mue / yn,i + 1);
                if (err > 0)
                    prn_emsg2(ctx, err);
                return(-1);
            }

            /* tda_out("yn=%g d=%g\n",yn,d); */

            if (ctx->PMYWVar >= 0 && ctx->PMD == 2)
                d /= yn;
     
            /* calculate working weights */

            if (ctx->WIVar >= 0) {
                wt = get_data(ctx, ctx->WIVar,i);
                if (wt < 0.0) {
                    printf1(ctx, "Error: found negative weight in case %d.\n",i + 1);
                    return(-1);
                }
            }
            wsum += wt;

            tmp = d * d;     
            switch (ctx->GLMD) {
                case 2:     tmp *= mue * (1.0 - mue / yn);
                            break;
                case 3:     tmp *= mue;                 
                            break;
                case 4:     tmp *= mue * mue;
                            break;
                case 5:     tmp *= mue * mue * mue;
                            break;
                default:    break;
            }
            if (wt > 0.0) {
                if (tmp <= 0.0) {
                    printf1(ctx, "Error: cannot calculate weight in case %d\n",i + 1);
                    return(-1);
                }
                w = sqrt(wt / tmp);
            }
            else     
                w = 0.0;
                 
            j = 1;
            if (intflg) {
                ctx->AcW[iw * nw + j] = w;
                j++;
            }
            for (k = 1; k <= nx; ++k) {
                tmp = get_data(ctx, ctx->PMVIdx[k],i);
                ctx->AcW[iw * nw + j] = tmp * w;
                j++;
            }
            if (w > 0.0)   
                ctx->AcW[iw * nw + j] = w * (eta  + (y - mue) * d);
            else
                ctx->AcW[iw * nw + j] = 0.0;

            if (cov && ctx->PMResFDef) {

                fprintf(ctx->PMResFd,"%6d ",i + 1);
                tmp = get_data(ctx, ctx->PMVIdx[0],i);
                rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMFmtS,y);
                rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMFmtS,mue);
                rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMFmtS,eta);
                for (k = 1; k <= nx1; ++k) {
                    tmp = ctx->AcW[iw * nw + k];           
                    if (w > 0.0)
                        tmp /= w;
                    rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMFmtS,tmp);           
                }
                rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMFmtS,wt);
                rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMFmtS,w * w);
                fprintf(ctx->PMResFd,"\n");
            }
            iw++;
        }
        if (cov && ctx->PMResFDef) {
            printf1(ctx, "Data and estimated values written to: %s\n",ctx->PMResFName);
            if (ctx->PMTDAFDef)  
                glm_dtda(ctx, intflg,nx);
            newline(ctx);
            if (cov == 2)
                return(1);
        }
        if (wsum < ctx->EPSI1) {
            printf1(ctx, "Error: sum of weights is almost zero.\n");
            return(-1);
        }
        if (ne > 0) {       /* add equality constraints */
            j = ne * nw;
            for (i = 1; i <= j; ++i)
                ctx->AcW[i] = ctx->AcU[i];
        }
        if (ni > 0) {       /* add inequality constraints */
            j = ni * nw;
            k = (ne + ctx->NOC) * nw;
            for (i = 1; i <= j; ++i)
                ctx->AcW[k + i] = ctx->AcV[i];
        }

        /**** 
        for (i = 0; i < NOC + ne + ni; ++i) {
            for (j = 1; j <= nw; ++j) 
                printf1(ctx, "%f ",AcW[i * nw + j]);
            printf1(ctx, "\n");
        }
        newline(ctx);
        **********/
         
        /*  solve least squares problem */

        err = lsei(ctx, ctx->AcW,ne,ctx->NOC,ni,nx1,ctx->AcX,cov,&rnorme,&rnorml,&ranka,&ranke);

        if (err == -1) {
            printf1(ctx, "Equality constraints are contradictory.\n");
            printf1(ctx, "Calculated a least squares solution.\n");
        }
        else if (err) {
            printf1(ctx, "Error: cannot solve least squares problem (%d)\n",err);
            if (err == -2)  
                printf1(ctx, "Inequality constraints are incompatible.\n");
            else if (err == -3)  
                printf1(ctx, "Equality and inequality constraints are contradictory.\n");
            else if (err == -4)  
                p_err(ctx, -2,1);
            return(-1);
        }
        if (ranka != nx1 - ranke) {
            printf1(ctx, "Error: rank of least squares data matrix: %d\n",ranka);
            return(-1);
        }
        ctx->GLMPC = 0.0;
        for (j = 1; j <= nx1; ++j) {
            tmp = fabs(ctx->AcX[j] - ctx->AcY[j]) / dmax(ctx, fabs(ctx->AcX[j]),1.0);                  
            ctx->GLMPC = dmax(ctx, ctx->GLMPC,tmp);                                            
        }
        if (ctx->SILENTFlg < 2)
            printfe(ctx, "  %3d  %20.13e\n",iter,ctx->GLMPC);

        if (ctx->PMProtFDef) {
            fprintf(ctx->PMProtFd,"Rank of least squares data matrix: %d\n",ranka);
            fprintf(ctx->PMProtFd,"Norm of least squares residuals: %g\n",rnorml); 
            if (ne > 0) {
                fprintf(ctx->PMProtFd,"Rank of equality constraints: %d\n",ranke);
                fprintf(ctx->PMProtFd,"Norm of residuals of equality constraints: %g\n",rnorme); 
            }
            prvec(ctx, "\nParameter vector",nx1,ctx->AcX);
            prval(ctx, "Scaled parameter change",ctx->GLMPC);
        }
        if (ctx->GLMPC <= ctx->TOLSP) {
            ctx->GLMIT = iter;
            ctx->GLMRK = ranka;

            /* do not calculate cov matrix if ... */
            if (ni > 0 || ctx->NOC <= ranka || fabs(rnorme + rnorml) < ctx->EPSI1) {
                if (ctx->PMResFDef) {
                    if (cov)
                        return(1);
                    else
                        cov = 2;
                }
                else
                    return(1);
            }
            else if (cov == 1) {    /* convergence reached and covariance
                                       matrix successfully calculated */

                /* if not normal distribution, we need to rescale the
                   covariance matrix */

                if (ctx->NOC <= ctx->GLMRK || rnorml < ctx->EPSI)
                    return(1);

                if (ctx->GLMD == 1)
                    return(0);

                tmp = (double)(ctx->NOC - ctx->GLMRK) / (rnorml * rnorml);
                for (i = 1; i <= nx1; ++i) {
                    for (j = 1; j <= nx1; ++j)
                        ctx->AcW[(i - 1) * nw + j] *= tmp;
                }
                return(0);
            }
            else if (cov == 2)
                return(1);
            else  
                cov = 1;
        }
    }
    return(-2);
}

/* ------------------------------------------------------------------------ */
/*  get_mue(eta,mue)    get mue such that g(mue) = eta.                     */
/*                      return 0 if successful, -1 if error.                */

int get_mue(TDAContext *ctx, double eta,double *mue)
{
    int err;
    double dx;

    dx = 0.0;
    err = 0;
    switch (ctx->GLMLF) {  
        case 0:     err = rtbis(ctx, ctx->PMX,ctx->PMY,mue,eta,&dx);   /* user-defined */
                    break;
        case 1:     *mue = eta;                         /* identity */
                    break;
        case 2:     *mue = rexp(ctx, eta);                   /* log */
                    break;
        case 3:     dx = rexp(ctx, eta);                     /* logit */
                    *mue = dx / (1.0 + dx);
                    break;
        case 4:     if (eta != 0.0)                     /* reciprocal */
                        *mue = 1.0 / eta;
                    else
                        err = -1;
                    break;
        case 5:     *mue = cdnf(ctx, eta);                   /* probit */
                    break;
        case 6:     dx = rexp(ctx, eta);                     /* compl. log-log */
                    *mue = 1.0 - rexp(ctx, -dx);
                    break;
        case 7:     *mue = eta * eta;                   /* square root */
                    break;
        case 8:     if (eta <= 0.0)                     /* quadratic inverse */
                        err = -1;
                    *mue = sqrt(1.0 / eta);
                    break;
    }
    if (err) {
        printf1(ctx, "Error: cannot find mue corresponding to eta = %g [dx=%g]\n",eta,dx);
        if (ctx->GLMLF == 0)
            printf1(ctx, "Check link function and range for mue.\n");
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  rtbis(x1,x2,x,eta)      find x such that g(x) = eta.                    */
/*  ###                     g(x) is the link function                       */
/*                          x1,x2 should be brackets such that              */
/*                          g(x1) <= eta <= g(x2) or                        */
/*                          g(x2) <= eta <= g(x1)                           */
/*                                                                          */
/*  Return 0 if successful, otherwise -1.                                   */

#define RTJMAX 500

int rtbis(TDAContext *ctx, double x1,double x2,double *x,double eta,double *dx)
{
    int j,err;
    double f1,f2,fmid,xmid;
    double tol,eps;

    eps = 10.e-15;
    tol = ctx->TOLF;     
    *dx = 0.0;
    ctx->FNArgVal[0] = x1 + eps;
    err = get_fival(ctx, 0,&f1,0,0);
    /* printf1("x1=%g f1=%g err=%d\n",x1,f1,err); */ 
    if (err)
        return(err);
    f1 -= eta;

    ctx->FNArgVal[0] = x2 - eps;
    err = get_fival(ctx, 0,&f2,0,0);
    /* printf1("x2=%g f2=%g err=%d\n",x2,f2,err); */
    if (err)
        return(err);
    f2 -= eta;

    if (f1 * f2 >= 0.0)  
        return(-1);
 
    if (f1 < 0.0) {
        *x = x1;
        *dx = x2 - x1;
    }
    else {
        *x = x2;
        *dx = x1 - x2;
    }
    for (j = 0; j < RTJMAX; ++j) {
        *dx *= 0.5;
        xmid = *x + *dx;
        ctx->FNArgVal[0] = xmid;
        err = get_fival(ctx, 0,&fmid,0,0);
        if (err)
            return(err);
        fmid -= eta;
        if (fmid <= 0.0)
            *x = xmid;
        if (fabs(*dx) < tol || fmid == 0.0) {
            /* printf1("eta=%20.16e x=%20.16e\n",eta,*x); */ 
            return(0);
        }
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  glm_lfd(mue,d)  Evaluate link function with mue, and get first          */
/*                  derivative, return in d.                                */
/*                  If OK return 0, else -1.                                */

int glm_lfd(TDAContext *ctx, double mue,double *d)
{
    int err;
    double tmp;

    *d = 0.0;
    err = 0;
    switch (ctx->GLMLF) {  
        case 0:     ctx->FNArgVal[0] = mue;                  /* user-defined */
                    err = get_fival(ctx, 0,&tmp,1,0);
                    if (err == 0)
                        *d = ctx->FNGrad[0][0];
                    break;
        case 1:     *d = 1.0;                           /* identity */
                    break;
        case 2:     if (mue == 0.0)                     /* log */
                        err = -1;
                    else
                        *d = 1.0 / mue;
                    break;
        case 3:     tmp = mue * (1.0 - mue);            /* logit */
                    if (tmp == 0.0)          
                        err = -1;
                    else
                        *d = 1.0 / tmp;
                    break;
        case 4:     if (mue == 0.0)                     /* reciprocal */
                        err = -1;         
                    else
                        *d = -1.0 / (mue * mue);
                    break;
        case 5:     if (mue <= 0.0 || mue >= 1.0)
                        err = -1;
                    else {
                        tmp = cdnif1(ctx, mue);
                        *d = 1.0 / dnf(ctx, tmp);            /* probit */
                    }
                    break;
        case 6:     tmp = 1.0 - mue;
                    if (tmp <= 0.0)                     /* compl. log-log */
                        err = -1;            
                    else  
                        *d = -1.0 / (tmp * rlog(ctx, tmp));
                    break;
        case 7:     if (mue <= 0.0)                     /* square root */
                        err = -1;
                    else
                        *d = 1.0 / (2.0 * sqrt(mue));
                    break;
        case 8:     tmp = mue * mue * mue;              /* quadratic inverse */
                    if (tmp == 0.0)
                        err = -1;
                    else
                        *d = -2.0 / tmp;
                    break;
    }
    return(err);
}
  
/* ------------------------------------------------------------------------ */
/*  glm_dev(nx,intflg,df)   This function calculates:                       */
/*                                                                          */
/*  GLMDEV = deviance                                                       */
/*  GLMX   = Pearson statistic                                              */
/*  PHI_ML = ML-based scale parameter                                       */
/*  PHI_DB = deviance-based scale parameter                                 */
/*                                                                          */  
/*  Return:  0  if successful, -1 if error.                                 */

int glm_dev(TDAContext *ctx, int nx,int intflg,int df)
{
    register int i,j,k;
    double eta,y,yn,mue,wt,wsum,tmp,tmp1,dev;

    yn = 1.0;
    wt = 1.0;
    wsum = 0.0;

    ctx->GLMDEV = 0.0;   /* deviance */
    ctx->GLMX = 0.0;     /* Pearson statistic */
    ctx->PHI_ML = 0.0;   /* ML-based scale parameter */
    ctx->PHI_DB = 0.0;   /* deviance-based scale parameter */

    for (i = 0; i < ctx->NOC; ++i) {

        y = get_data(ctx, ctx->PMVIdx[0],i);
        eta = 0.0;
        j = 1;
        if (intflg)
            eta += ctx->AcX[j++];
        for (k = 1; k <= nx; ++k) {
            tmp = get_data(ctx, ctx->PMVIdx[k],i);
            eta += tmp * ctx->AcX[j++];
        }
        if (get_mue(ctx, eta,&mue))  
            return(-1);

        if (ctx->PMYWVar >= 0 && ctx->PMD == 2) {
            yn = get_data(ctx, ctx->PMYWVar,i);
            if (y < 0.0 || yn < 1.0 || y > yn) {
                printf1(ctx, "Error in case %d: y = %g, count = %g\n",i+1,y,yn);
                return(-1);
            }
            mue *= yn;
        }
        if (ctx->WIVar >= 0) {
            wt = get_data(ctx, ctx->WIVar,i);
            if (wt < 0.0) {
                printf1(ctx, "Error: found negative weight in case %d.\n",i + 1);
                return(-1);
            }
        }
        wsum += wt;

        dev = 0.0;
        switch (ctx->GLMD) {
            case 2:     if (y > 0.0)                /* binomial */
                            dev += 2.0 * y * rlog(ctx, y / mue);
                        tmp = yn - y;
                        tmp1 = yn - mue;
                        if (tmp > 0.0 && tmp1 > 0.0)
                            dev += 2.0 * tmp * rlog(ctx, tmp / tmp1);
                        break;
            case 3:     if (y > 0.0)                /* Poisson */
                            dev = 2.0 * y * rlog(ctx, y / mue);
                        dev -= 2.0 * (y - mue);
                        break;

            case 4:     if (y > 0.0 && mue > 0.0)              /* Gamma */
                            dev = -2.0 * rlog(ctx, y / mue);
                        if (mue > 0.)
                            dev += 2.0 * (y - mue) / mue;
                        break;

            case 5:     if (y > 0.0 && mue > 0.0)   /* Inverse Gaussian */
                            dev = (y - mue) * (y - mue) / (y * mue * mue);
                        break;

            default:    tmp = y - mue;              /* normal */
                        tmp1 = tmp * tmp;
                        dev = tmp1;
                        ctx->GLMX += tmp1;
                        ctx->PHI_ML += tmp1;
                        break;
        }
        ctx->GLMDEV += wt * dev;
    }
    if (wsum < ctx->EPSI1) {
        printf1(ctx, "Error: sum of weights is almost zero.\n");
        return(-1);
    }
    switch (ctx->GLMD) {
        case  2:    break;
        case  3:    break;
        case  4:    break;
        case  5:    break;
        default:    ctx->PHI_ML /= wsum;
                    break;
    }
    if (df > 0)
        ctx->PHI_DB = ctx->GLMDEV / (double)df;

    printf1(ctx, "\nRank of data matrix: %d\n",ctx->GLMRK);
    printf1(ctx, "Degrees of freedom: %d\n",df);
    prn_sfmt(ctx, "Deviance",31,ctx->PMTFmtS,ctx->GLMDEV);

    if (ctx->GLMD == 1) {
        prn_sfmt(ctx, "Pearson statistic",31,ctx->PMTFmtS,ctx->GLMX);
        prn_sfmt(ctx, "ML-based scaling factor",31,ctx->PMTFmtS,ctx->PHI_ML);
        prn_sfmt(ctx, "Deviance-based scaling factor",31,ctx->PMTFmtS,ctx->PHI_DB);
    }
    else if (ctx->GLMD == 2)
        prn_sfmt(ctx, "Fixed scaling factor",31,ctx->PMTFmtS,1.0);


    return(0);
}

/* ------------------------------------------------------------------------ */
/*  glm_dtda()      write TDA description file                              */

void glm_dtda(TDAContext *ctx, int intflg,int nx)
{
    register int i,j,k;
   
    if (ctx->PMTDAFDef) {
        fprintf(ctx->PMTDAFd,"# data written by glm command.\n");
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

        fprintf(ctx->PMTDAFd,"  Mue");
        fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - 3,0);
        fprintf(ctx->PMTDAFd,"[%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,k++);

        fprintf(ctx->PMTDAFd,"  Eta");
        fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - 3,0);
        fprintf(ctx->PMTDAFd,"[%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,k++);

        if (intflg) {
            fprintf(ctx->PMTDAFd,"  Int");
            fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - 3,0);
            fprintf(ctx->PMTDAFd,"[%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,k++);
        }
        for (i = 1; i <= nx; ++i) {
            j = ctx->PMVIdx[i];
            fprintf(ctx->PMTDAFd,"  %s",ctx->VName[j]);
            fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - (int)strlen(ctx->VName[j]),0);
            fprintf(ctx->PMTDAFd,"[%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,k++);
        }
        fprintf(ctx->PMTDAFd,"  Weight");
        fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - 6,0);
        fprintf(ctx->PMTDAFd,"[%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,k++);

        fprintf(ctx->PMTDAFd,"  WWeight");
        fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - 7,0);
        fprintf(ctx->PMTDAFd,"[%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,k++);
        fprintf(ctx->PMTDAFd,");\n");
    
        printf1(ctx, "TDA description written to: %s\n",ctx->PMTDAFName);
    }
}


