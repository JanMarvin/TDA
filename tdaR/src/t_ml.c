/****************************************************************************/
/*  t_ml                                                                    */
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
#include "t_gf.h"
#include "t_fnrta.h"
#include "t_fnrtb.h"
#include "t_fnrtc.h"
#include "t_fnqr.h"
#include "t_gdat.h"
#include "t_parm.h"
#include "t_con.h"
#include "t_tmin.h"
#include "t_edat.h"
#include "t_rate.h"
#include "t_gmin.h"
#include "t_lin.h"
#include "t_eval.h"
#include "t_int.h"
#include "t_mat.h"
#include "t_matf.h"
#include "t_mds.h"
#include "t_cl.h"
#include "t_ireg.h"
#include "tda_context.h"

/*  functions in t_ml.c */

void set_mldef(TDAContext *ctx);
void set_mlopt(TDAContext *ctx);
int get_dsv(TDAContext *ctx, int n,double *par,int gmina,double *lb,double *ub,int nl);
int ml_init(TDAContext *ctx, int opt,int typ,int gmina);
int syminv1(TDAContext *ctx, int n, double *d, double *h);
int syminv2(TDAContext *ctx, int n, double *x);
void fn(TDAContext *ctx, double *x, int mod, int typ, int scal, int opt,int *err);
int fn_umod(TDAContext *ctx);
void prjvec(TDAContext *ctx, double *x);
void prjhess(TDAContext *ctx);
int checkov(TDAContext *ctx);
void prn_mlres(TDAContext *ctx, int typ);
void prn_ml1res(TDAContext *ctx, int typ);
void prvec(TDAContext *ctx, char *txt, int n, double *x);
void prmat(TDAContext *ctx, char *txt,int n,int m,double *x);
void prhess(TDAContext *ctx, char *txt, int n);
void prval(TDAContext *ctx, char *txt, double x);
void prot_init(TDAContext *ctx, int typ,int mod,int gmina);
void prot_prni(TDAContext *ctx, char *s,int n,int *x);
void prot_prns(TDAContext *ctx, char *s,int n,short *x);
void prot_prnc(TDAContext *ctx, char *s,int n,char *x);

/*--------------------------------------------------------------------------*/
/*  global variables defined in t_ml.c                                      */



/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/


                            /* derivatives in step size search.             */



                            /* hessian requested                            */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/



                            /* used to handle constraints                   */
                            /* used to handle constraints!                  */



/*--------------------------------------------------------------------------*/
/*  set_mldef()     set default options for ML estimation.                  */
         
void set_mldef(TDAContext *ctx)
{
    ctx->MINA = 5;          /* type of min algorithm */
    ctx->MxIter = -1;        /* max number of iterations */
    ctx->MxItFlg = 0;        /* set if mxit command used */
    ctx->MxIt1 = -1;         /* max number of subiterations */
    ctx->MxItR = 100;        /* max number of random search iterations */
    ctx->Crite = 1;          /* type of convergence criterion */
    ctx->CritFlg = 0;        /* Set if user defined criterion */
    ctx->DScal = -1.0;       /* Default scaling factor */
    ctx->DScalFlg = 0;       /* Set to 1 if dscal command */
    ctx->TOLG = 1.e-6;       /* Tolerance for gradient */
    ctx->TOLF = 1.e-12;      /* Tolerance for function value changes */
    ctx->TOLP = 1.e-4;       /* Tolerance for parameter changes */
    ctx->TOLV = 1.e-10;      /* Tolerance for Simplex algorithm */
    ctx->TOLS = 1.e-4;       /* Tolerance for direct search algorithm */
    ctx->TOLSG = 1.e-5;      /* Tolerance for scaled gradient */
    ctx->TOLSP = 1.e-8;      /* Tolerance for scaled parameter change */

    ctx->AMue = 0.2;         /* Mue for Armijo conditions in fmin */
    ctx->SMin = 1.e-10;      /* Minimal step size */
    ctx->SLen = 1.0;         /* Step length in direct search algorithm */
    ctx->SRed = 0.5;         /* Step reduction in direct search algorithm */
    ctx->CCheck = 5;         /* Check convergence interval in Simplex alg. */
    ctx->STFlg = 0;          /* Set if function evaluation without derivatives */
    ctx->CCTyp = 2;          /* Type of covariance matrix calculation */
    ctx->SRLen = 1.0;        /* Step length for random search algorithm */
    ctx->SRRed = 0.1;        /* Step length reduction for random search */
}

/*--------------------------------------------------------------------------*/
/*  set_mlopt()     adjust ML options after the interpreation in parm().    */  
         
void set_mlopt(TDAContext *ctx)
{
    if (ctx->MINA == 7 || ctx->MINA == 8) {

        if (ctx->PMDOPT < 0 || ctx->PMDOPT > 2)
            ctx->PMDOPT = 2;

        ctx->TSDERIV = ctx->PMDOPT;
        if (ctx->MxIter < 0) {
            if (ctx->TSDERIV == 0)
                ctx->MxIter = 200;
            else if (ctx->TSDERIV == 1)
                ctx->MxIter = 100;
            else
                ctx->MxIter = 20;
        }
        ctx->Crite = 6;
    }
    if (ctx->MINA == 1)
        ctx->Crite = 5;
    else if (ctx->MINA == 2) 
        ctx->Crite = 4;
    else if (ctx->MINA == 3 || ctx->MINA == 4)
        ctx->Crite = 1;

    if (ctx->MxIter < 0) {
        if (ctx->MINA <= 4)
            ctx->MxIter = 100;
        else
            ctx->MxIter = 20;
    }
    if (ctx->MxIt1 < 0)
        ctx->MxIt1 = 50;
}

/* ------------------------------------------------------------------------ */
/*  get_dsv(n,par,gmina,lb,ub,nl)                                           */
/*                                                                          */
/*              Get starting values into par[i], i = 1,...,n.               */
/*              First check xp=... , then input file PMDSVName.             */
/*              If gmina != 0 put lower and upper bounds into lb[],ub[].    */
/*              Default is 0 +/- 1. If nl != 0 additional newline.          */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */
/*              If successful set DSVFlg = 1.                               */
/*                                                                          */

int get_dsv(TDAContext *ctx, int n,double *par,int gmina,double *lb,double *ub,int nl) 
{
    FILE *fd;
    register int i;
    int err,m,fflag;
    char buf[202],*p;
    double tmp;

    ctx->DSVFlg = 0;
    err = -1;
    fflag = m = 0;
    if (ctx->PMNX > 0) {
        m = imin(ctx, n,ctx->PMNX);
        for (i = 1; i <= m; ++i) {
            par[i] = ctx->PMXX[i - 1];
            if (gmina) {
                lb[i] = ctx->PMXA[i - 1];
                ub[i] = ctx->PMXB[i - 1];
            }
        }
        printf1(ctx, "%d starting value(s) from xp parameter.\n",m);
        ctx->DSVFlg = 1;
    }
    else if (ctx->PMDSVFlg) {
                     
        printf1(ctx, "Reading starting values from: %s\n",ctx->PMDSVName);
        if (!(fd = fopen(ctx->PMDSVName,OPEN_RD))) {   
            printf1(ctx, "Error: can't open %s\n",ctx->PMDSVName);
            goto GETDSVFin;
        }
        fflag = 1;

        while (fgets(buf,200,fd)) {
               
            if (check_drec(ctx, buf)) {  /* if a data record */
                p = skip_b(ctx, buf);    /* skip blanks */
                if (sscanf(p,"%lg",&tmp) != 1)  
                    goto GETDSVFin;
                par[m + 1] = tmp;
                if (gmina) {
                    p = skip_dbl(ctx, p);
                    p = skip_b(ctx, p);
                    if (sscanf(p,"%lg",&tmp) == 1) {
                        lb[m + 1] = tmp;
                        p = skip_dbl(ctx, p);
                        p = skip_b(ctx, p);
                    }
                    else
                        lb[m + 1] = par[m + 1] - 1.0;

                    if (sscanf(p,"%lg",&tmp) == 1)  
                        ub[m + 1] = tmp;
                    else
                        ub[m + 1] = par[m + 1] + 1.0;
                }
                m++;
            }
            if (m >= n)
                break;
        }
        if (m > 0)  
            ctx->DSVFlg = 1;
    }
    if (ctx->DSVFlg) {
        if (m < n)  
            printf1(ctx, "Warning: found %d, would need %d starting values.\n",m,n);
        if (nl)
            newline(ctx);
    }
    for (i = m + 1; i <= n; ++i) {
        par[i] = 0.0;       
        if (gmina) {
            lb[i] = -1.0;
            ub[i] =  1.0;
        }
    }
    err = 0;   

GETDSVFin:
    if (fflag)
        fclose(fd);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ml_init(opt,typ,gmina)                                                  */  
/*                                                                          */
/*  Initialization for function minimization.                               */
/*                                                                          */
/*  typ = 0 : ML                                                            */
/*        1 : general func minimization                                     */
/*        2 : nonlinear regression                                          */
/*                                                                          */
/*  If opt != 0 allocate required memory for model estimation with          */
/*  NParm parameters. If opt == 0 free previously allocated memory.         */
/*  If gmina != 0 allocate additional memory for parameter bounds.          */
/*                                                                          */
/*  Return 0 if OK, -1 if error (insufficient memory).                      */

int ml_init(TDAContext *ctx, int opt,int typ,int gmina)
{
    int err = 0;

    ctx->DSVFlg = 0;             /* flag for externally supplied starting values */
    clear_npflg(ctx);          /* clear numerical problem flags */
    ctx->NINTMUsed = -1;

    if (opt == 0) {
        if (ctx->MLAlloc == 0)
            return(0);
        goto MLIFin;
    }
    err = -1;

    if (typ == 0)
        printf1(ctx, "Maximum likelihood estimation.");
    else if (typ == 1)
        printf1(ctx, "Function minimization.");
    else
        printf1(ctx, "Nonlinear regression.");

    printf1(ctx, "\nAlgorithm %d: ",ctx->MINA);

    switch (ctx->MINA) {
        case  1:    printf1(ctx, "Direct Search\n\n");
                    break;
        case  2:    printf1(ctx, "Simplex\n\n");
                    break;
        case  3:    printf1(ctx, "conjugate gradients\n\n");
                    break;
        case  4:    printf1(ctx, "BFGS\n\n");
                    break;
        case  6:    printf1(ctx, "Newton (II)\n\n");
                    break;
        case  7:    printf1(ctx, "CES (quadratic model [dopt=%d])\n\n",ctx->TSDERIV);
                    break;
        case  8:    printf1(ctx, "CES (tensor model [dopt=%d])\n\n",ctx->TSDERIV);
                    break;
        default:    printf1(ctx, "Newton (I)\n\n"); /* MINA 5 is default */
                    break;
    }
    ctx->NParmA = ctx->NParm;
    ctx->HSiz = ctx->NParm * (ctx->NParm - 1) / 2;     /* size of hessian lower triangular */
    ctx->NParm1 = ctx->NParm;                     /* dimensions in reduced space */
    ctx->HSiz1 = ctx->HSiz;

    ctx->HSizA = ctx->HSiz;
    if (ctx->MINA == 3) {                  /* Hess is used as a working area */
        if (ctx->HSizA < 5 * ctx->NParm + 2)     /* in ffmin5(ctx). */
            ctx->HSizA = 5 * ctx->NParm + 2;
    }
    else if (ctx->MINA == 4) {
        if (ctx->HSizA < ctx->NParm * (ctx->NParm + 7) / 2)
            ctx->HSizA = ctx->NParm * (ctx->NParm + 7) / 2;
    }

    /* print info about ML estimation */

    printf1(ctx, "Number of model parameters: %d\n",ctx->NParm);
    if (ctx->CCTyp)
        printf1(ctx, "Type of covariance matrix: %d\n",ctx->CCTyp);
    printf1(ctx, "Maximum number of iterations: %d\n",ctx->MxIter);

    if (ctx->MINA == 1 || ctx->MINA == 2)  
        printf1(ctx, "Starting value of step length: %lg\n",ctx->SLen);
       
    if (ctx->MINA == 1)   
        printf1(ctx, "Step length reduction factor: %lg\n",ctx->SRed);
       
    printf1(ctx, "Convergence criterion: %d\n",ctx->Crite);

    if (ctx->Crite == 1)
        printf1(ctx, "Tolerance for norm of final gradient: %lg\n",ctx->TOLG);
    else if (ctx->Crite == 2)
        printf1(ctx, "Tolerance for change in function value: %lg\n",ctx->TOLF);
    else if (ctx->Crite == 3)
        printf1(ctx, "Tolerance for change in parameters: %lg\n",ctx->TOLP);
    else if (ctx->Crite == 4) 
        printf1(ctx, "Tolerance for variance of function values: %lg\n",ctx->TOLV);
    else if (ctx->Crite == 5)  
        printf1(ctx, "Minimum value of step length: %lg\n",ctx->TOLS);
    else if (ctx->Crite == 6) {
        printf1(ctx, "Tolerance for scaled gradient: %lg\n",ctx->TOLSG); 
        printf1(ctx, "Tolerance for scaled parameter change: %lg\n",ctx->TOLSP);
    }
    if (ctx->MINA == 5 || ctx->MINA == 6)  
        printf1(ctx, "Mue of Armijo condition: %lg\n",ctx->AMue);
       
    if (ctx->MINA >= 3 && ctx->MINA <= 6)  
        printf1(ctx, "Minimum of step size value: %lg\n",ctx->SMin);
       
    if (typ > 0 && ctx->DScal <= 0.0)
        ctx->DScal = 1.0;
   
    printf1(ctx, "Scaling factor: %lg\n\n",ctx->DScal);

    /* set initial conditions */

    ctx->NumF = ctx->NumFG = ctx->NumFH = 0;

    /* memory allocation */

    if (!(ctx->Par = (double *)calloc((size_t)(ctx->NParmA + 1),sizeof(double))))   
        goto MLIFin;         
    memrq(ctx, ctx->NParmA + 1,sizeof(double));
    ctx->ParA = ctx->NParmA + 1;
              
    if (gmina) {
        if (!(ctx->ParLB = (double *)calloc((size_t)(ctx->NParmA + 1),sizeof(double))))   
            goto MLIFin;         
        memrq(ctx, ctx->NParmA + 1,sizeof(double));
        ctx->ParLBA = ctx->NParmA + 1;

        if (!(ctx->ParUB = (double *)calloc((size_t)(ctx->NParmA + 1),sizeof(double))))   
            goto MLIFin;         
        memrq(ctx, ctx->NParmA + 1,sizeof(double));
        ctx->ParUBA = ctx->NParmA + 1;
    }         
    if (!(ctx->WrkD = (double *)calloc((size_t)(ctx->NParm + 1),sizeof(double))))   
        goto MLIFin;         
    memrq(ctx, ctx->NParm + 1,sizeof(double));
    ctx->WrkDA = ctx->NParm + 1;

    if (!(ctx->WrkE = (double *)calloc((size_t)(ctx->NParm + 1),sizeof(double))))   
        goto MLIFin;          
    memrq(ctx, ctx->NParm + 1,sizeof(double));
    ctx->WrkEA = ctx->NParm + 1;

    if (!(ctx->WrkH = (double *)calloc((size_t)(ctx->HSiz + 1),sizeof(double))))   
        goto MLIFin;         
    memrq(ctx, ctx->HSiz + 1,sizeof(double));
    ctx->WrkHA = ctx->HSiz + 1;

    if (!(ctx->Par1 = (double *)calloc((size_t)(ctx->NParm + 1),sizeof(double))))   
        goto MLIFin;         
    memrq(ctx, ctx->NParm + 1,sizeof(double));
    ctx->Par1A = ctx->NParm + 1;
              
    if (!(ctx->Grad = (double *)calloc((size_t)(ctx->NParm + 1),sizeof(double))))   
        goto MLIFin;         
    memrq(ctx, ctx->NParm + 1,sizeof(double));
    ctx->GradA = ctx->NParm + 1;
              
    if (!(ctx->Diag = (double *)calloc((size_t)(ctx->NParmA + 1),sizeof(double))))   
        goto MLIFin;         
    memrq(ctx, ctx->NParmA + 1,sizeof(double));
    ctx->DiagA = ctx->NParmA + 1;
              
    if (!(ctx->Hess = (double *)calloc((size_t)(ctx->HSizA + 1),sizeof(double))))   
        goto MLIFin;         
    memrq(ctx, ctx->HSizA + 1,sizeof(double));
    ctx->HessA = ctx->HSizA + 1;
              
    if (ctx->CCTyp != 2 || ctx->MINA == 7 || ctx->MINA == 8) {
        if (!(ctx->GTmp = (double *)calloc((size_t)(ctx->NParm + 1),sizeof(double))))   
            goto MLIFin;         
        memrq(ctx, ctx->NParm + 1,sizeof(double));
        ctx->GTmpA = ctx->NParm + 1;
    }         
    if (ctx->MINA >= 5 && ctx->MINA <= 8) {
        if (!(ctx->Srch = (double *)calloc((size_t)(ctx->NParm + 1),sizeof(double))))   
            goto MLIFin;         
        memrq(ctx, ctx->NParm + 1,sizeof(double));
        ctx->SrchA = ctx->NParm + 1;
    }         
    if (ctx->MINA == 7 || ctx->MINA == 8) {
        if (!(ctx->TSHess = (double *)calloc((size_t)(ctx->NParm) * (size_t)(ctx->NParm) + 1,sizeof(double))))   
            goto MLIFin;         
        ctx->TSHessA = ctx->NParm * ctx->NParm + 1;
        memrq(ctx, ctx->TSHessA,sizeof(double));
              
        if (!(ctx->TSGrad = (double *)calloc((size_t)(ctx->NParm + 1),sizeof(double))))   
            goto MLIFin;         
        ctx->TSGradA = ctx->NParm + 1;
        memrq(ctx, ctx->TSGradA,sizeof(double));
    }         
    ctx->MLAlloc = 1;
    return(0);

MLIFin:
    if (ctx->SrchA) {
        free((char *)ctx->Srch);
        memrq(ctx, -ctx->SrchA,sizeof(double));
        ctx->SrchA = 0;
    }
    if (ctx->GTmpA) {
        free((char *)ctx->GTmp);
        memrq(ctx, -ctx->GTmpA,sizeof(double));
        ctx->GTmpA = 0;
    }
    if (ctx->HessA) {
        free((char *)ctx->Hess);
        memrq(ctx, -ctx->HessA,sizeof(double));
        ctx->HessA = 0;
    }
    if (ctx->DiagA) {
        free((char *)ctx->Diag);
        memrq(ctx, -ctx->DiagA,sizeof(double));
        ctx->DiagA = 0;
    }
    if (ctx->GradA) {
        free((char *)ctx->Grad);
        memrq(ctx, -ctx->GradA,sizeof(double));
        ctx->GradA = 0;
    }
    if (ctx->Par1A) {
        free((char *)ctx->Par1);
        memrq(ctx, -ctx->Par1A,sizeof(double));
        ctx->Par1A = 0;
    }
    if (ctx->WrkHA) {
        free((char *)ctx->WrkH);
        memrq(ctx, -ctx->WrkHA,sizeof(double));
        ctx->WrkHA = 0;
    }
    if (ctx->WrkEA) {
        free((char *)ctx->WrkE);
        memrq(ctx, -ctx->WrkEA,sizeof(double));
        ctx->WrkEA = 0;
    }
    if (ctx->WrkDA) {
        free((char *)ctx->WrkD);
        memrq(ctx, -ctx->WrkDA,sizeof(double));
        ctx->WrkDA = 0;
    }
    if (ctx->ParA) {
        free((char *)ctx->Par);
        memrq(ctx, -ctx->ParA,sizeof(double));
        ctx->ParA = 0;
    }
    if (ctx->ParLBA) {
        free((char *)ctx->ParLB);
        memrq(ctx, -ctx->ParLBA,sizeof(double));
        ctx->ParLBA = 0;
    }
    if (ctx->ParUBA) {
        free((char *)ctx->ParUB);
        memrq(ctx, -ctx->ParUBA,sizeof(double));
        ctx->ParUBA = 0;
    }
    if (ctx->TSHessA) {
        free((char *)ctx->TSHess);
        memrq(ctx, -ctx->TSHessA,sizeof(double));
        ctx->TSHessA = 0;
    }
    if (ctx->TSGradA) {
        free((char *)ctx->TSGrad);
        memrq(ctx, -ctx->TSGradA,sizeof(double));
        ctx->TSGradA = 0;
    }
    ctx->MLAlloc = 0;
    if (err)  
        p_err(ctx, -2,1);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  syminv1 Inversion of the hessian matrix with dimension                  */
/*          n = NParm or NParm1.                                            */
/*          The diagonal elements are assumed to be in d, the lower         */
/*          tiangular is assumed to be in h. The inversion is performed     */
/*          directly in the lower triangular. The function returns the      */
/*          rank of the matrix. For the determination of the pseudo rank    */
/*          the tolerance EPSI1 is used.                                    */
/*                                                                          */
/*  Return -1 if insufficient memory.                                       */

int syminv1(TDAContext *ctx, int n, double *d, double *h)
{
    register int i,j,k,kk,ii,jj;
    int aflag,r;
    double crit,tmp1,tmp2,tmp3,*wrk;

    if (!(wrk = (double *)calloc((size_t)(n + 1),sizeof(double))))
        return(-1); 
    memrq(ctx, n + 1,sizeof(double));
                
    r = n;         
    for (i = 1; i <= n; ++i)            
        wrk[i] = d[i];

    for (k = 1; k <= n; ++k) {

        kk = (k - 1) * (k - 2) / 2;

        aflag = 0;                          /* flag for aliased parameters  */

        /* check for current sweep-column */

        tmp1 = d[k];
        if (tmp1 <= ctx->EPSI1 * wrk[k]) {
            aflag = 1;
        }
        else {
            for (i = 1; i < k; ++i) {
                tmp2 = d[i];
                if (tmp2 != 0.0) {
                    tmp3 = h[kk + i];
                    crit = 1.0 / (tmp2 + (tmp3 * tmp3) / tmp1);
                    if (crit < ctx->EPSI1 * wrk[i]) {
                        aflag = 1;
                        break;
                    }
                }
            }
        }
        if (!aflag) {

            /* collinearity check passed, so updating triangle */
 
            for (j = 1; j < k; ++j) {
                jj = (j - 1) * (j - 2) / 2;
                tmp2 = h[kk + j];
                tmp3 = tmp2 / tmp1;
                for (i = 1; i < j; ++i)  
                    h[jj + i] += h[kk + i] * tmp3;
                d[j] += tmp2 * tmp3;
            }
            for (i = k + 1; i <= n; ++i) {
                ii = (i - 1) * (i - 2) / 2;
                tmp2 = h[ii + k] / tmp1;
                for (j = k + 1; j <= i; ++j) {
                    jj = (j - 1) * (j - 2) / 2;
                    if (j < i)
                        h[ii + j] -= h[jj + k] * tmp2;
                    else
                        d[j] -= h[jj + k] * tmp2;
                }
            }
            for (i = k + 1; i <= n; ++i) {
                ii = (i - 1) * (i - 2) / 2;
                tmp2 = h[ii + k] / tmp1;
                for (j = 1; j < k; ++j)  
                    h[ii + j] -= h[kk + j] * tmp2;
            }
            for (j = 1; j < k; ++j)
                h[kk + j] /= tmp1;
     
            for (i = k + 1; i <= n; ++i) {
                ii = (i - 1) * (i - 2) / 2;
                h[ii + k] /= -tmp1;
            }
            d[k] = 1.0 / tmp1;
        }
        else {
 
            /*  collinearity check not passed, so updating the k column     */
            /*  and k row with zeroes                                       */
 
            for (i = 1; i < k; ++i)  
                h[kk + i] = 0.0;
            d[k] = 0.0;

            for (i = k + 1; i <= n; ++i) {
                ii = (i - 1) * (i - 2) / 2;
                h[ii + k] = 0.0;
            }
            r--;
        }
    }
    memrq(ctx, -n - 1,sizeof(double));
    free((char *)wrk);
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  syminv2 Inversion of a pos. def. symmetric matrix x with dimension n.   */
/*          It is assumed that the array x only contains the lower          */
/*          triangle, including the diagonal elements.                      */
/*                                                                          */
/*          The inversion is performed directly in the lower triangle.      */
/*          The function returns the rank of the matrix. For the            */
/*          determination of the pseudo rank tolerance EPSI1 is used.       */
/*                                                                          */
/*  Return -1 if insufficient memory.                                       */

int syminv2(TDAContext *ctx, int n, double *x)
{
    register int i,j,k,kk,ii,jj;
    int aflag,r;
    double crit,tmp1,tmp2,tmp3,*wrk;

    if (!(wrk = (double *)calloc((size_t)(n + 1),sizeof(double))))
        return(-1); 
    memrq(ctx, n + 1,sizeof(double));
                
    r = n;         
    for (i = 1; i <= n; ++i)            
        wrk[i] = x[(i * i + i) / 2];

    for (k = 1; k <= n; ++k) {

        kk = ((k - 1) * k) / 2;

        aflag = 0;                          /* flag for aliased parameters  */

        /* check for current sweep-column */

        tmp1 = x[kk + k];

        if (tmp1 <= ctx->EPSI1 * wrk[k]) {
            aflag = 1;
        }
        else {
            for (i = 1; i < k; ++i) {
                ii = ((i - 1) * i) / 2;
                tmp2 = x[ii + i];
                if (tmp2 != 0.0) {
                    tmp3 = x[kk + i];
                    crit = 1.0 / (tmp2 + (tmp3 * tmp3) / tmp1);
                    if (crit < ctx->EPSI1 * wrk[i]) {
                        aflag = 1;
                        break;
                    }
                }
            }
        }
        if (!aflag) {

            /* collinearity check passed, so updating triangle */
 
            for (j = 1; j < k; ++j) {
                jj = ((j - 1) * j) / 2;
                tmp2 = x[kk + j];
                tmp3 = tmp2 / tmp1;
                for (i = 1; i <= j; ++i)  
                    x[jj + i] += x[kk + i] * tmp3;
            }
            for (i = k + 1; i <= n; ++i) {
                ii = ((i - 1) * i) / 2;
                tmp2 = x[ii + k] / tmp1;
                for (j = k + 1; j <= i; ++j) {
                    jj = ((j - 1) * j) / 2;
                    x[ii + j] -= x[jj + k] * tmp2;
                }
            }
            for (i = k + 1; i <= n; ++i) {
                ii = ((i - 1) * i) / 2;
                tmp2 = x[ii + k] / tmp1;
                for (j = 1; j < k; ++j)  
                    x[ii + j] -= x[kk + j] * tmp2;
            }
            for (j = 1; j < k; ++j)
                x[kk + j] /= tmp1;
     
            for (i = k + 1; i <= n; ++i) {
                ii = ((i - 1) * i) / 2;
                x[ii + k] /= -tmp1;
            }
            x[kk + k] = 1.0 / tmp1;
        }
        else {
 
            /*  collinearity check not passed, so updating the k column     */
            /*  and k row with zeroes                                       */
 
            for (i = 1; i <= k; ++i)  
                x[kk + i] = 0.0;

            for (i = k + 1; i <= n; ++i) {
                ii = ((i - 1) * i) / 2;
                x[ii + k] = 0.0;
            }
            r--;
        }
    }
    memrq(ctx, -n - 1,sizeof(double));
    free((char *)wrk);
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  fn  Call of model-specific functions for calculation of loglikelihood,  */
/*      gradient and hessian. Parameters are assumed to be in x[], a        */
/*      vector of length NParm1. If constraints are defined the vector is   */
/*      expanded before calling the model-specific functions.               */
/*                                                                          */
/*      The function value is calculated in FTmp. The number of function    */
/*      calls is counted in the global variable NumF, calls with            */
/*      derivatives are counted in NumFG and NumFH.                         */
/*                                                                          */
/*      If scal != 0, the calculated function value, the gradient, and      */
/*      the hessian, are scaled with DScal before returning.                */
/*                                                                          */
/*      If typ = 0  only function value is calculated                       */
/*      If typ = 1  function value and gradient are calculated              */
/*      If typ = 2  function value, gradient, and hessian are calculated    */
/*                                                                          */
/*      If opt == 1, calculation of score test statistics based on log-     */
/*      likelihood in the unrestricted parameter space.                     */
/*      Results saved in FNTStat and FNTRank.                               */
/*                                                                          */
/*      If first call of fn(), and if GCTestFlg = 1, calculation and        */
/*      print of global chi square (score test statistic with actual        */
/*      starting values). Only done if typ = 2.                             */
/*      First call recognized by global variable FN_First.                  */
/*                                                                          */
/*      The actual number of cases used to evaluate the function is         */
/*      counted in NOCUsed.                                                 */
/*                                                                          */
/*      Return err = -1 if insufficient memory, 1 if error in function.     */

void fn(TDAContext *ctx, double *x, int mod, int typ, int scal, int opt,int *err)
{
    (void)opt;        /* unused: the signature is shared */
    register int i,j;
    double tmp; /** ,tmp1, *d,*h **/  

    *err = 0;
    ctx->NumF++;
    ctx->LFunc = 1;
    ctx->NOCUsed = ctx->LGrad = ctx->LSec = 0;
    if (typ == 1) {
        ctx->LGrad = 1;
        ctx->NumFG++;
    }
    else if (typ == 2) {    /* LSec = 1 implies LGrad = 1 */
        ctx->LGrad = ctx->LSec = 1;
        ctx->NumFG++;
        ctx->NumFH++;
    }
    if (ctx->NCon1) {                             /* expand x */
        for (i = 1; i <= ctx->NParm; ++i) {
            tmp = 0.0;
            for (j = 1; j <= ctx->NParm1; ++j)
                tmp += ctx->ConQ[(ctx->CIP[i] - 1) * ctx->NParm1 + j] * x[j];
            ctx->ParC[i] = ctx->ParS[i] + tmp;
        }
        ctx->TPar = ctx->ParC;
    }
    else
        ctx->TPar = x;

    if (ctx->PMProtFDef >= 3) {
        fprintf(ctx->PMProtFd,"\nFunction evaluation %d (%d,%d). ",ctx->NumF,ctx->NumFG,ctx->NumFH);
        fprintf(ctx->PMProtFd,"Func %d Grad %d Sec %d\n",ctx->LFunc,ctx->LGrad,ctx->LSec);
    }
    if (ctx->LFunc)  
        ctx->FTmp = 0.0;                                              
    if (ctx->LGrad)
        dclear(ctx, ctx->NParm,ctx->Grad);                               
    if (ctx->LSec) {
        dclear(ctx, ctx->NParm,ctx->Diag);                                 
        dclear(ctx, ctx->HSiz,ctx->Hess);                                   
    }
    switch (mod) {

        /* quantal response models */

        case QRLOG1:    *err = fn_logit1(ctx);
                        break;
        case QRLOG2:    *err = fn_logit2(ctx);
                        break;
        case QRLOG3:    *err = fn_logit3(ctx);
                        break;
        case QRLOG4:    *err = fn_logit4(ctx);
                        break;

        case QRPROB1:   *err = fn_probit1(ctx);
                        break;
        case QRPROB2:   *err = fn_probit2(ctx);
                        break;
        case QRPROB3:   *err = fn_probit3(ctx);
                        break;
        case QRPROB4:   *err = fn_probit4(ctx);
                        break;

        case RMOD1:     *err = fn_rmod1(ctx);      
                        break;
                        
        /* transition rate models */

        case MCOX:      *err = fn_cox(ctx);      
                        break;
        case MEXP:      *err = fn_exp(ctx, 0);   
                        break;
        case MEXP1:     *err = fn_exp1(ctx, 0);  
                        break;
        case MEXP2:     *err = fn_exp2(ctx, 0);  
                        break;
        case MPOL:      *err = fn_pol(ctx, 0);  
                        break;
        case MPOL1:     *err = fn_pol1(ctx, 0);
                        break;
        case MGM:       *err = fn_gm(ctx, 0);   
                        break;
        case MWEI:      *err = fn_wei(ctx, 0);
                        break;
        case MSIC:      *err = fn_sic(ctx, 0);   
                        break;
        case MLL:       *err = fn_ll(ctx, 0);
                        break;
        case MLL2:      *err = fn_ll2(ctx, 0);
                        break;
        case MLN:       *err = fn_ln(ctx, 0);   
                        break;
        case MIG:       *err = fn_ig(ctx, 0);   
                        break;
        case MGAM:      *err = fn_gam(ctx, 0);      
                        break;
        case DLR:   
        case CLL:       *err = fn_dis(ctx, mod);  
                        break;

        /* general functions */

        case CLMOD1:    *err = clu_fn(ctx);   /* clu */
                        break;

        /*******************************************
        case MDSMOD1:   *err = mds1_fn();   
                        break;
        case MDSMOD2:   *err = mds2_fn();   
                        break;
        ********************************************/
        case MDSMOD3:   *err = mds3_fn(ctx);   /* MDS */
                        break;
        case MDSMOD4:   *err = mds4_fn(ctx);   /* MDS */
                        break;

        case IVRMOD1:   *err = ivr1_fn(ctx);   /* interval regression */
                        break;

        case UMOD:      *err = fn_umod(ctx);   /* user-defined models */
                        break;

        default:        printf1(ctx, "Error: fn %d not available\n",mod);
                        exit(0);
    }
    if (ctx->PMProtFDef >= 3) {
        prval(ctx, "Function",ctx->FTmp);
        prvec(ctx, "Parameter",ctx->NParm,ctx->TPar);
        if (ctx->LGrad)  
            prvec(ctx, "Gradient",ctx->NParm,ctx->Grad);
        if (ctx->LSec)   
            prhess(ctx, "Hessian",ctx->NParm);
    }   
    if (*err)
        return;
           
   
    /* calculate score test statistic if opt == 1, or if first call
       of fn(ctx) and GCTestFlg = 1. */
    /***
    if (opt == 1 || (GCTestFlg && FN_First)) {

        memrq(ctx, NParm + HSiz + 2,sizeof(double));
        if (!(d = (double *)calloc(NParm + 1,sizeof(double))) ||
            !(h = (double *)calloc(HSiz + 1,sizeof(double))))   
            p_err(ctx, -2,1);
         
        for (i = 1; i <= NParm; ++i)    
            d[i] = -Diag[i];
        for (i = 1; i <= HSiz; ++i)
            h[i] = -Hess[i];

        check memory    

        if ((FNTRank = syminv1(ctx, NParm,d,h)) != NParm) {   
            FNTStat = -1.0;
        }
        else {
            FNTStat = 0.0;
            for (i = 1; i <= NParm; ++i) {
                tmp1 = 0.0;
                for (j = 1; j <= NParm; ++j) {

                    if (i == j)
                        tmp = d[j];
                    else if (i < j)
                        tmp = h[(j - 1) * (j - 2) / 2 + i];
                    else
                        tmp = h[(i - 1) * (i - 2) / 2 + j];

                    tmp1 += Grad[j] * tmp;
                }
                FNTStat += tmp1 * Grad[i];
            }
        }
        memrq(ctx, -NParm - HSiz - 2,sizeof(double));
        free((char *)d);
        free((char *)h);

        if (GCTestFlg && FN_First) {
            GCTest = FNTStat;
            GCRank = FNTRank;
        }
    }
    *******/

    if (ctx->LFunc && ctx->FN_First)
        ctx->FMax1 = ctx->FTmp;

    if (ctx->LFunc && scal)
        ctx->FTmp *= ctx->DScal;                                      

    if (ctx->LGrad) {
        if (scal) {
            for (j = 1; j <= ctx->NParm; ++j)  
                ctx->Grad[j] *= ctx->DScal;                              
        }
        if (ctx->NCon1)   
            prjvec(ctx, ctx->Grad);
    }
    if (ctx->LSec) {
        if (scal) {
            for (j = 1; j <= ctx->NParm; ++j)   
                ctx->Diag[j] *= ctx->DScal;                                 
            for (j = 1; j <= ctx->HSiz; ++j)
                ctx->Hess[j] *= ctx->DScal;                                             
        }
        if (ctx->NCon1)  
            prjhess(ctx);
    }
    ctx->FN_First = 0;
}

/* ------------------------------------------------------------------------ */
/*  fn_umod()   calculate function, gradient, hessian for user-defined      */
/*              models. Use parameter values in TPar[] (i=1,NParm).         */
/*              if LFunc  return function value in FTmp.                    */
/*              if LGrad  return gradient in Grad[]                         */
/*              if LSec   return hessian in Hess[] and Diag[]               */
/*                                                                          */
/*              If CCovFlg = 1 return outer product of gradients instead    */
/*              of Hessian.                                                 */
/*                                                                          */
/*              Return 0 if OK, -1 if insuff memory, 1 if error in function */

int fn_umod(TDAContext *ctx)
{
    int i,err,deriv;

    ctx->NOCUsed = deriv = 0;
    if (ctx->LGrad)
        deriv = 1;
    if (ctx->LSec)  
        deriv = 2;

    err = get_flval(ctx, &ctx->FTmp,ctx->NParm,ctx->TPar + 1,deriv,2,ctx->Grad,ctx->Hess,ctx->Diag);

    if (err) {
        printf1(ctx, "\nError in function evaluation.\n");
        prn_emsg2(ctx, err);
        printf1(ctx, "Current function arguments:");
        for (i = 1; i <= ctx->NParm; ++i)
            printf1(ctx, " %lg",ctx->TPar[i]);
        printf1(ctx, "\n");
        err = 1;
    }
    else {
        err = checkov(ctx);          /* check overflow */
        if (err) {
            printf1(ctx, "\nError in function evaluation.\n");
            printf1(ctx, "Numerical overflow.\n");
        }
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prjvec(x)   Project the NParm vector x to NParm1 = NParm - NCon1 by     */
/*              multiplication with the transposed projection matrix ConQ.  */
/*              Use of WrkD as temporary storage.                           */

void prjvec(TDAContext *ctx, double *x)
{
    register int i,j;
    double tmp;

    for (i = 1; i <= ctx->NParm1; ++i) {
        tmp = 0.0;
        for (j = 1; j <= ctx->NParm; ++j)  
            tmp += ctx->ConQ[(ctx->CIP[j] - 1) * ctx->NParm1 + i] * x[j];
        ctx->WrkD[i] = tmp;
    }
    for (i = 1; i <= ctx->NParm1; ++i)  
        x[i] = ctx->WrkD[i];
}

/* ------------------------------------------------------------------------ */
/*  prjhess     Project the Hessian (Hess and Diag) into the NParm1 =       */
/*              NParm - NCon1 space.                                        */
/*              Use of WrkH as temporary storage.                           */

void prjhess(TDAContext *ctx)
{
    register int i,j,k,l;
    double tmp;

    for (i = 1; i <= ctx->NParm1; ++i) {
        for (j = 1; j <= i; ++j) {
            tmp = 0.0;
            for (k = 1; k <= ctx->NParm; ++k) {
                for (l = 1; l <= ctx->NParm; ++l) {
                    if (l == k)
                        tmp += ctx->Diag[k] * ctx->ConQ[(ctx->CIP[l] - 1) * ctx->NParm1 + i] * 
                                         ctx->ConQ[(ctx->CIP[k] - 1) * ctx->NParm1 + j];
                    else if (l < k)
                        tmp += ctx->Hess[(k - 1) * (k - 2) / 2 + l] *
                                         ctx->ConQ[(ctx->CIP[l] - 1) * ctx->NParm1 + i] * 
                                         ctx->ConQ[(ctx->CIP[k] - 1) * ctx->NParm1 + j];
                    else
                        tmp += ctx->Hess[(l - 1) * (l - 2) / 2 + k] *
                                         ctx->ConQ[(ctx->CIP[l] - 1) * ctx->NParm1 + i] * 
                                         ctx->ConQ[(ctx->CIP[k] - 1) * ctx->NParm1 + j];
                }
            }
            if (i == j)
                ctx->WrkD[i] = tmp;
            else
                ctx->WrkH[(i - 1) * (i - 2) / 2 + j] = tmp;
        }
    }
    for (i = 1; i <= ctx->NParm1; ++i)  
        ctx->Diag[i] = ctx->WrkD[i];
    for (i = 1; i <= ctx->HSiz1; ++i)
        ctx->Hess[i] = ctx->WrkH[i];
}

/* ------------------------------------------------------------------------ */
/*  checkov.  Check overflow and underflow. OFMax and OFMin are used as     */
/*            ceilings and floors.                                          */
/*            Return 0 if OK,  1 if error.                                  */

int checkov(TDAContext *ctx)
{
    register int j,k = 0;

    if (ctx->LFunc) {
        if (ctx->FTmp <= ctx->OFMin) {
            k = 1;
            ctx->FTmp  = ctx->OFMin;
        }
        if (ctx->FTmp >= ctx->OFMax) {
            k = 1;  
            ctx->FTmp  = ctx->OFMax;
        }
    }
    if (ctx->LGrad) {
        for (j = 1; j <= ctx->NParm; ++j) {
            if (ctx->Grad[j] <= ctx->OFMin) {
                k = 1;           
                ctx->Grad[j] =  ctx->OFMin;
            }
            if (ctx->Grad[j] >= ctx->OFMax) {
                k = 1;
                ctx->Grad[j]  = ctx->OFMax;
            }
        }
    }
    if (ctx->LSec) {
        for (j = 1; j <= ctx->NParm; ++j) {
            if (ctx->Diag[j] <= ctx->OFMin) {
                k = 1;  
                ctx->Diag[j] =  ctx->OFMin;
            }
            if (ctx->Diag[j] >= ctx->OFMax) {
                k = 1;   
                ctx->Diag[j]  = ctx->OFMax;
            }
        }
        for (j = 1; j <= ctx->HSiz; ++j) {
            if (ctx->Hess[j] <= ctx->OFMin) {
                k = 1;   
                ctx->Hess[j] =  ctx->OFMin;
            }
            if (ctx->Hess[j] >= ctx->OFMax) {
                k = 1;   
                ctx->Hess[j]  = ctx->OFMax;
            }
        }
    }
    ctx->NPFlgs[2] = (short)(ctx->NPFlgs[2] + (k)); 
    return(k); 
}

/* ------------------------------------------------------------------------ */
/*  prn_mlres(typ)      Print results of function minimization.             */
/*                      typ 0 : ML                                          */
/*                          1 : general function                            */
/*                          2 : nonlinear regression                        */

void prn_mlres(TDAContext *ctx, int typ)
{
    printf1(ctx, "\nConvergence ");
    if (ctx->LConv) 
        printf1(ctx, "not ");
    printf1(ctx, "reached in %d iterations.\n",ctx->Iter);
    printf1(ctx, "Number of function evaluations: %d (%d,%d)\n",ctx->NumF,ctx->NumFG,ctx->NumFH);
    if (ctx->NINTMUsed >= 0)
        ni_info(ctx);

    if (ctx->LConv) {
        printf1(ctx, "\nProblem: ");

        switch (ctx->LConv) {
            case  1:    if (typ == 0)
                            printf1(ctx, "cannot increase log likelihood.\n");
                        else
                            printf1(ctx, "cannot decrease function.\n");
                        break;
            case  2:    printf1(ctx, "reached max number of iterations.\n");
                        break;
            case -1:    printf1(ctx, "final hessian (or outer product) not positive definite.\n");
                        break;
            case -2:    printf1(ctx, "step size search failed.\n");
                        break;
            case -3:    printf1(ctx, "cannot calculate eigenvalues and vectors of the hessian.\n");
                        break;
            case -4:    printf1(ctx, "search vector is no descent direction.\n");
                        break;
            case -5:    printf1(ctx, "exceeded max iterations in step size search.\n");
                        break;
            default:    printf1(ctx, "LConv = %d\n",ctx->LConv);
                        break;
        }
    }
    if (typ == 0)
        printf1(ctx, "\nMaximum of log likelihood: %lg",ctx->FMax);
    else
        printf1(ctx, "\nMinimum of function: %lg",ctx->FMax);
#ifdef TDA_R_PACKAGE
    /* each diagnostic beside its own print, full precision where the
       text has %lg's six digits (CONTRIBUTING.md) */
    tda_export_mat(ctx, typ == 0 ? "ml.logLik" : "ml.fmin",
                   &ctx->FMax, 1, 1);
#endif

    if (ctx->MINA == 2)  
        printf1(ctx, "\nVariance of final function values: %lg",ctx->CValV);
       
    if (ctx->LConv >= 0) {

        if ((ctx->CValG) != 0.0)
            printf1(ctx, "\nNorm of final gradient vector: %lg",ctx->CValG);
#ifdef TDA_R_PACKAGE
        if ((ctx->CValG) != 0.0)
            tda_export_mat(ctx, "ml.gradient", &ctx->CValG, 1, 1);
#endif

        if (ctx->MINA == 5 || ctx->MINA == 6) {
            printf1(ctx, "\nLast absolute change of function value: %lg",ctx->CValF);
            printf1(ctx, "\nLast relative change in parameters: %lg",ctx->CValP);
#ifdef TDA_R_PACKAGE
            tda_export_mat(ctx, "ml.change.f", &ctx->CValF, 1, 1);
            tda_export_mat(ctx, "ml.change.p", &ctx->CValP, 1, 1);
#endif
        }
        else if (ctx->MINA == 7 || ctx->MINA == 8) {
            printf1(ctx, "\nFinal scaled gradient: %lg",ctx->TOLSGF);
            printf1(ctx, "\nFinal scaled parameter change: %lg",ctx->TOLSPF);
#ifdef TDA_R_PACKAGE
            tda_export_mat(ctx, "ml.scaled.g", &ctx->TOLSGF, 1, 1);
            tda_export_mat(ctx, "ml.scaled.p", &ctx->TOLSPF, 1, 1);
#endif
        }
    }
    printf1(ctx, "\n\n");
    prn_npflg(ctx);           /* print numerical problems if any */
}

/* ------------------------------------------------------------------------ */
/*  prn_ml1res(typ)     Print additional results of func minimization.      */
/*                      typ 0 : ML                                          */
/*                          1 : general function                            */
/*                          2 : nonlinear regression                        */
/*                                                                          */
/*                      do only if LConv >= 0.                              */
/*                      if CCTyp > 0 and PMCovFDef print cov. matrix.       */

void prn_ml1res(TDAContext *ctx, int typ)
{
    if (typ == 0)
        printf1(ctx, "Log likelihood");
    else
        printf1(ctx, "Function");
    printf1(ctx, " (starting values): ");
    rt_printf1_d(ctx, ctx->PMTFmtS,ctx->FMax1); 

    if (typ == 0)
        printf1(ctx, "\nLog likelihood");
    else
        printf1(ctx, "\nFunction");
    printf1(ctx, " (final estimates): "); 
    rt_printf1_d(ctx, ctx->PMTFmtS,ctx->FMax); 
    newline(ctx);    
    newline(ctx);
#ifdef TDA_R_PACKAGE
    /* the starting/final pair beside its print (tfmt text is usually
       full precision already; the export removes the dependence on it) */
    {
        double erow[2];
        erow[0] = ctx->FMax1;
        erow[1] = ctx->FMax;
        tda_export_mat(ctx, typ == 0 ? "ml.logLik.pair" : "ml.f.pair",
                       erow, 1, 2);
    }
#endif
    if (ctx->PMPPWFlg || ctx->PMCovWFlg) {
        if (ctx->PMPPWFlg)  
            printf1(ctx, "Parameter estimates written to: %s\n",ctx->PMPPFName);
        if (ctx->PMCovWFlg)  
            printf1(ctx, "Covariance matrix written to: %s\n",ctx->PMCovFName);
        newline(ctx);
    }
    mp_info(ctx);              /* print info about matrices */

    /***
    if (GCTestFlg)         score test with starting values   
        gctest();

    if (LMTestFlg)         LM test   
        lmtest();
      
    if (NWTest)            Wald tests   
        wtest();
    ****/
}

/* ------------------------------------------------------------------------ */
/*  prvec(txt,n,x)                                                          */
/*      Prints a vector with n elements. The string txt is printed first.   */
/*      Used for iteration protocol.                                        */

void prvec(TDAContext *ctx, char *txt, int n, double *x)
{
    register int i; 

    fprintf(ctx->PMProtFd,"%s\n",txt); 
    for (i = 1; i <= n; ++i)  
        rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,x[i]);
    fprintf(ctx->PMProtFd,"\n");
}

/* ------------------------------------------------------------------------ */
/*  prmat(txt,n,m,x)                                                        */
/*      Prints a matrix with n,m elements. The string txt is printed first. */
/*      Used for protocol file.                                             */

void prmat(TDAContext *ctx, char *txt,int n,int m,double *x)
{
    register int i,j;

    fprintf(ctx->PMProtFd,"%s\n",txt); 
    for (i = 0; i < n; ++i) {
        for (j = 1; j <= m; ++j) 
            rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,x[i * m + j]);
        fprintf(ctx->PMProtFd,"\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  prhess(txt,n)       Print Hessian to protocol file.                     */

void prhess(TDAContext *ctx, char *txt, int n)
{
    register int i,j,k;

    fprintf(ctx->PMProtFd,"%s\n",txt); 
    k = 1;
    for (i = 1; i <= n; ++i) {
        for (j = 1; j < i; ++j)
            rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,ctx->Hess[k++]);
        rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,ctx->Diag[i]);
        fprintf(ctx->PMProtFd,"\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  prval(txt,x)                                                            */
/*      Print a double value x. The string txt is printed first.            */
/*      Used for iteration protocol.                                        */

void prval(TDAContext *ctx, char *txt, double x)
{
    fprintf(ctx->PMProtFd,"%s\n",txt); 
    rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,x);  
    fprintf(ctx->PMProtFd,"\n");
}

/* ------------------------------------------------------------------------ */
/*  prot_init(typ,mod,gmina)                                                */
/*                                                                          */
/*  Init protocol file                                                      */
/*                      typ 1 = transition rate models                      */
/*                          2 = quantal response models                     */
/*                          3 = user-def. log-likelihood                    */
/*                          4 = general function                            */
/*                          5 = nonlinear regression                        */
/*                          6 = user-def rate model                         */
/*                                                                          */
/*  If typ 1 add some information about transitions and parameter vector.   */
/*  If gmina != 0 print lower/upper parameter bounds.                       */
/*                                                                          */

void prot_init(TDAContext *ctx, int typ,int mod,int gmina)
{
    register int i,j;

    if (ctx->PMProtFDef == 0)
        return;

    printf1(ctx, "Protocol will be written to: %s\n",ctx->PMProtFName);

    switch (typ) {
        case 1:     fprintf(ctx->PMProtFd,"Transition rate model: %d\n\n",mod);
                    break;
        case 2:     fprintf(ctx->PMProtFd,"Quantal response model: %d\n\n",mod);
                    break;
        case 3:     fprintf(ctx->PMProtFd,"User-defined log-likelihood.\n\n");
                    break;
        case 4:     fprintf(ctx->PMProtFd,"General function minimization.\n\n");
                    break;
        case 5:     fprintf(ctx->PMProtFd,"Nonlinear regression.\n\n");
                    break;
        case 6:     fprintf(ctx->PMProtFd,"User-defined rate model.\n\n");
                    break;
        default:    break;
    }
    fprintf(ctx->PMProtFd,"Number of model parameters: %d\n",ctx->NParm);
    fprintf(ctx->PMProtFd,"Starting values");
    if (gmina)
        fprintf(ctx->PMProtFd," and bounding box");
    fprintf(ctx->PMProtFd,".\n");

    for (i = 1; i <= ctx->NParm; ++i) {
        fprintf(ctx->PMProtFd,"%3d  ",i);
        rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,ctx->Par[i]);
        if (gmina) {
            rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,ctx->ParLB[i]);
            rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,ctx->ParUB[i]);
        }
        fprintf(ctx->PMProtFd,"\n");
    }
    fprintf(ctx->PMProtFd,"\n");

    if (typ == 1) {     /* transition rate models */

        fprintf(ctx->PMProtFd,"Number of transitions: %d\n",ctx->NTran1);

        for (i = 1; i <= ctx->MaxSnn; ++i) {
            fprintf(ctx->PMProtFd,"TranPtr:  %2d  :",i);
            for (j = 0; j <= ctx->MaxOrg; ++j)
                fprintf(ctx->PMProtFd," %2d",ctx->TranPtr[(i - 1) * ctx->MaxOrg1 + j]);
            fprintf(ctx->PMProtFd,"\n");
        }
        prot_prni(ctx, "SnTran1: ",ctx->NTran1,ctx->SnTran1);
        prot_prni(ctx, "OrgTran1:",ctx->NTran1,ctx->OrgTran1);
        prot_prni(ctx, "DesTran1:",ctx->NTran1,ctx->DesTran1);

        prot_prns(ctx, "VTNum:   ",ctx->NTran1,ctx->VTNum);
        prot_prns(ctx, "PIdxPtr: ",ctx->NTran1,ctx->PIdxPtr);
        prot_prns(ctx, "PIdx:    ",ctx->NParm,ctx->PIdx + 1);
        prot_prnc(ctx, "PIdxM:   ",ctx->NParm,ctx->PIdxM + 1);
        fprintf(ctx->PMProtFd,"\n");
    }
}

void prot_prni(TDAContext *ctx, char *s,int n,int *x)
{
    register int i;

    fprintf(ctx->PMProtFd,"%s",s);
    for (i = 0; i < n; ++i)
        fprintf(ctx->PMProtFd," %2d",x[i]);
    fprintf(ctx->PMProtFd,"\n");
}

void prot_prns(TDAContext *ctx, char *s,int n,short *x)
{
    register int i;

    fprintf(ctx->PMProtFd,"%s",s);
    for (i = 0; i < n; ++i)
        fprintf(ctx->PMProtFd," %2d",(int)x[i]);
    fprintf(ctx->PMProtFd,"\n");
}

void prot_prnc(TDAContext *ctx, char *s,int n,char *x)
{
    register int i;

    fprintf(ctx->PMProtFd,"%s",s);
    for (i = 0; i < n; ++i)
        fprintf(ctx->PMProtFd," %2d",(int)x[i]);
    fprintf(ctx->PMProtFd,"\n");
}




