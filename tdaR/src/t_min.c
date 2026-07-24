/****************************************************************************/
/*  t_min                                                                   */
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
#include "t_gf.h"
#include "t_ml.h"
#include "t_svd.h"
#include "t_rand.h"
#include "t_con.h"
#include "t_alloc.h"
#include "t_tmin.h"
#include "t_gdat.h"
#include "t_int.h"
#include "t_mat.h"
#include "t_matf.h"
#include "tda_context.h"

/*  functions in t_min.c */

int ffmin(TDAContext *ctx, int mod,int typ,int cov,int msg);
int ffmin1(TDAContext *ctx, int mod);
int ffmin2(TDAContext *ctx, int mod);
int ffmin3(TDAContext *ctx, int mod);
int ffmin5(TDAContext *ctx, int mod);
void exphess(TDAContext *ctx);
int chinv(TDAContext *ctx, int opt);
int ffmin7(TDAContext *ctx, int mod);

/* ------------------------------------------------------------------------ */
/*  ffmin(mod,typ,cov,msg)      Function minimization                       */
/*                                                                          */
/*          typ 0 : ML                                                      */
/*              1 : gen function                                            */
/*              2 : nonlinear regression                                    */
/*                                                                          */
/*          If cov == 0 do not calculate covariance matrix.                 */
/*          If msg == 0 do not print messages.                              */
/*                                                                          */
/*          a) Call specified minimization algorithm. If constraints are    */
/*             defined, this is done in the restricted parameter space.     */
/*          b) If successful and cov != 0 additional call of fn_function to */
/*             calculate the covariance matrix (this depends on the user-   */
/*             specified type, given by CCTyp). If CCTyp = 0 the cov.       */
/*             matrix is not calculated.                                    */
/*                                                                          */
/*             If typ = 2, the cov matrix is scaled with                    */
/*                                                                          */
/*              2 * f / (NOC - NParm)                                       */
/*                                                                          */
/*          c) If constraints are defined, the final parameter vector and   */
/*             the covariance matrix are expanded to the orginal space.     */
/*          d) Print result to the iteration protocol.                      */
/*          e) If requested print covariance matrix to output file.         */
/*                                                                          */
/*          Result:                                                         */
/*              LConv depending on the success of the function.             */
/*              FMax  log likelihood at maximum                             */
/*              Par   parameter vector in the original parameter space      */
/*              Hess/Diag   covariance matrix in the original space.        */
/*                                                                          */
/*          LConv 1 cannot increase log likelihood                          */
/*                2 exceeded max number of iterations                       */
/*               -1 final hessian (or outer product) not positive definite  */
/*               -2 step size search failed                                 */
/*               -3 cannot calculate eigenvalues and vectors of the hessian */
/*               -4 search vector is no descent direction                   */
/*               -5 exceeded max iterations in step size search             */
/*               -6 insufficient memory.                                    */
/*               -7 error in function evaluation                            */
/*                                                                          */
/*      Note: This function also calculates:                                */
/*      a)  LM test statistics, if LMTestFlg = 1.                           */
/*                                                                          */
/*      Return:     0 if successful.                                        */
/*                 -1 if insufficient memory.                               */
/*                 +1 if error in function evaluation.                      */

int ffmin(TDAContext *ctx, int mod,int typ,int cov,int msg)
{
    register int i,j,k,l,ij;
    int r,err,rr,n;
    double tmp,tmp1,tmp2,tmp3;

    err = -1;
    ctx->Iter = ctx->LConv = ctx->CGradFlg = ctx->CCovFlg = ctx->MCovFlg = ctx->NumF = ctx->NumFG = ctx->NumFH = 0;
    ctx->NINTMUsed = -1;
    clear_npflg(ctx);          /* clear numerical problem flags */

    if (ctx->MxIter < 1) {
        ctx->LConv = 2;
        return(0);
    }
    /* tda_out_flush(); */

    ctx->FN_First = 1;   /* set for first call of fn(ctx) */


    /*  Call minimization algorithm specified by MINA.
        The resulting parameters are given in Par with dimension NParm1.
        So, if constraints are defined, Par is in the reduced space.
        The evaluated function value is returned in FMin. */

    if (ctx->PMProtFDef) {
        fprintf(ctx->PMProtFd,"Function minimization with algorithm %d",ctx->MINA);
        if (ctx->MINA == 7 || ctx->MINA == 8)
            fprintf(ctx->PMProtFd," [dopt=%d]",ctx->PMDOPT);
        fprintf(ctx->PMProtFd,".\n");
    }

    switch (ctx->MINA) {

        case 1:     /* Direct search */
                    r = ffmin1(ctx, mod);
                    break;
        case 2:     /* Simplex method */
                    r = ffmin2(ctx, mod);
                    break;
        case 3:     /* conjugate gradient or BFGS methods */
        case 4:     r = ffmin3(ctx, mod);
                    break;

        case 7:     /* CES (qudratic and tensor model) */
        case 8:     r = ffmin7(ctx, mod);
                    break;

        default:    /* Newton algorithms (MINA 5 and 6) */
                    r = ffmin5(ctx, mod);
                    break;
    }
    if (r) {            /* insufficient memory or error in function */
        if (r < 0) {
            err = -1;
            ctx->LConv = -6;
        }
        else {
            err = r;
            ctx->LConv = -7;
        }
        goto FFMFin;
    }
    if (ctx->SILENTFlg < 2) {
        printfe(ctx, "\n");
        fflushe(ctx);
    }
    if (ctx->PMProtFDef) {
        switch (ctx->LConv) {
          case  0:  fprintf(ctx->PMProtFd,"\nConvergence reached with criterion %d.\n\n",ctx->Crite);
                    break;
          case  1:  fprintf(ctx->PMProtFd,"\nCannot increase log likelihood.\n\n");
                    break;
          case  2:  fprintf(ctx->PMProtFd,"\nExceeded max number of iterations.\n\n");
                    break;
          case -1:  fprintf(ctx->PMProtFd,"\nFinal hessian or outer product not positive definite.\n\n");
                    break;
          case -2:  fprintf(ctx->PMProtFd,"\nStep size search failed.\n\n");
                    break;
          case -3:  fprintf(ctx->PMProtFd,"\nCannot calculate eigenvalues and vectors of the hessian.\n\n");
                    break;
          case -4:  fprintf(ctx->PMProtFd,"\nSearch vector not in descent direction.\n\n");
                    break;
          case -5:  fprintf(ctx->PMProtFd,"\nExceeded max iterations in step size search.\n\n");
                    break;
        }
    }

    ctx->FMax = ctx->FMin / ctx->DScal;

    if (ctx->PMMPLogDef == 1)                /* create matrix for loglikelihood */
        mp_putlog(ctx, ctx->FMax);

    /*  If successful calculate covariance matrix. This depends on the
        type defined by CCTyp. All calculations are done in a parameter space
        which is reduced if constraints are present.
        Note: fn(ctx) is called without the scaling option.

        Note: If cov = 0 or CCTyp = 0, cov matrix is not calculated */

    if (ctx->PMMPGradDef == 1) {
        if (cov == 0 || ctx->CCTyp == 0 || ctx->NOCUsed < 1)
            ctx->PMMPGradDef = -1;
        else if (mp_alloc(ctx, 4,ctx->NOCUsed,ctx->NParm + ctx->PM2NV))
            ctx->PMMPGradDef = -1;
    }

    if (cov == 0 || ctx->CCTyp == 0)
        goto FFMCon;

    if (ctx->LConv >= 0 && ctx->CCTyp > 0) {

        if (ctx->PMProtFDef)
            fprintf(ctx->PMProtFd,"Covariance matrix calculation with ccov=%d\n\n",ctx->CCTyp);

        if (ctx->PMMPGradDef == 1) {     /* calculate gradients */
            ctx->CGradFlg = 1;
            fn(ctx, ctx->Par,mod,1,0,0,&rr);
            ctx->CGradFlg = 0;
            if (rr) {
                err = rr;
                goto FFMFin;
            }
        }

        ctx->MCovFlg = 1;    /* signals calculation of covariance matrix */

        /* calculate outer product */

        if (ctx->CCTyp == 1 || ctx->CCTyp == 3) {

            ctx->CCovFlg = 1;
            fn(ctx, ctx->Par,mod,2,0,0,&rr);
            ctx->CCovFlg = 0;
            if (rr) {
                err = rr;
                goto FFMFin;
            }

            if (ctx->CCTyp == 1) {
                if ((i = syminv1(ctx, ctx->NParm1,ctx->Diag,ctx->Hess)) != ctx->NParm1) {

                    if (i < 0) {    /* insufficient memory */
                        ctx->LConv = -6;
                        err = -1;
                        goto FFMFin;
                    }
                    ctx->LConv = -1;
                }
            }
            else {  /* save outer product in Hess1/Diag1 */

                if (!(ctx->Diag1 = (double *)calloc((size_t)(ctx->NParm1 + 1),sizeof(double)))) {
                    ctx->LConv = -6;
                    err = -1;
                    goto FFMFin;
                }
                ctx->Diag1A = ctx->NParm1 + 1;
                memrq(ctx, ctx->Diag1A,sizeof(double));

                if (!(ctx->Hess1 = (double *)calloc((size_t)(ctx->HSiz1 + 1),sizeof(double)))) {
                    ctx->LConv = -6;
                    err = -1;
                    goto FFMFin;
                }
                ctx->Hess1A = ctx->HSiz1 + 1;
                memrq(ctx, ctx->Hess1A,sizeof(double));

                for (i = 1; i <= ctx->NParm1; ++i)
                    ctx->Diag1[i] = ctx->Diag[i];
                for (i = 1; i <= ctx->HSiz1; ++i)
                    ctx->Hess1[i] = ctx->Hess[i];
            }
        }

        /*  calculate hessian in Hess / Diag */

        if (ctx->CCTyp == 2 || ctx->CCTyp == 3) {

            /* if algorithm 7 or 8 and we have only function values, then
               use the Hessian returned from ffmin7(ctx) */

            if ((ctx->MINA == 7 || ctx->MINA == 8) && ctx->PMDOPT == 0)
                goto FFMCONT;

            if (ctx->LMTestFlg && ctx->NCon1) {
                fn(ctx, ctx->Par,mod,2,0,1,&rr);       /* LM test if LMTestFlg = 1 */
                ctx->LMTest = ctx->FNTStat;
                ctx->LMRank = ctx->FNTRank;
            }
            else
                fn(ctx, ctx->Par,mod,2,0,0,&rr);

            if (rr) {
                err = rr;
                goto FFMFin;
            }
FFMCONT:
            if (typ == 0) {
                for (i = 1; i <= ctx->NParm1; ++i)
                    ctx->Diag[i] *= -1.0;
                for (i = 1; i <= ctx->HSiz1; ++i)
                    ctx->Hess[i] *= -1.0;
            }
            if ((i = syminv1(ctx, ctx->NParm1,ctx->Diag,ctx->Hess)) != ctx->NParm1) {
                if (i < 0) {
                    ctx->LConv = -6; /* insufficient memory */
                    err = -1;
                    goto FFMFin;
                }
                else
                    ctx->LConv = -1;

            }

            if (ctx->CCTyp == 3) {   /* calculate cross product */

                if (ctx->LConv >= 0) {

                    for (i = 1; i <= ctx->NParm1; ++i)
                        ctx->WrkD[i] = ctx->Diag[i];
                    for (i = 1; i <= ctx->HSiz1; ++i)
                        ctx->WrkH[i] = ctx->Hess[i];
                    ij = 1;
                    for (i = 1; i <= ctx->NParm1; ++i) {
                        for (j = 1; j <= i; ++j) {
                            tmp = 0.0;
                            for (k = 1; k <= ctx->NParm1; ++k) {
                                tmp1 = 0.0;
                                for (l = 1; l <= ctx->NParm1; ++l) {
                                    if (l == k)
                                        tmp3 = ctx->Diag1[k];
                                    else if (l < k)
                                        tmp3 = ctx->Hess1[(k - 1) * (k - 2) / 2 + l];
                                    else
                                        tmp3 = ctx->Hess1[(l - 1) * (l - 2) / 2 + k];

                                    if (j == l)
                                        tmp2 = ctx->WrkD[l];
                                    else if (j < l)
                                        tmp2 = ctx->WrkH[(l - 1) * (l - 2) / 2 + j];
                                    else
                                        tmp2 = ctx->WrkH[(j - 1) * (j - 2) / 2 + l];

                                    tmp1 += tmp3 * tmp2;
                                }
                                if (k == i)
                                    tmp2 = ctx->WrkD[i];
                                else if (k < i)
                                    tmp2 = ctx->WrkH[(i - 1) * (i - 2) / 2 + k];
                                else
                                    tmp2 = ctx->WrkH[(k - 1) * (k - 2) / 2 + i];

                                tmp += tmp1 * tmp2;
                            }
                            if (i == j)
                                ctx->Diag[i] = tmp;
                            else
                                ctx->Hess[ij++] = tmp;
                        }
                    }
                }
            }
        }
        ctx->MCovFlg = 0;

        tmp = 0.0;                      /* calculate norm of gradient */
        for (i = 1; i <= ctx->NParm1; ++i)
            tmp += ctx->Grad[i] * ctx->Grad[i];
        ctx->CValG = sqrt(tmp);
    }

    /*  If constraints are defined, expand the parameter vector and the
        covariance matrix. */

FFMCon:

    if (ctx->NCon1) {

        for (i = 1; i <= ctx->NParm; ++i) {
            tmp = 0.0;
            for (j = 1; j <= ctx->NParm1; ++j)
                tmp += ctx->ConQ[(ctx->CIP[i] - 1) * ctx->NParm1 + j] * ctx->Par[j];
            ctx->WrkD[i] = tmp;
        }
        for (i = 1; i <= ctx->NParm; ++i)
            ctx->Par[i] = ctx->ParS[i] + ctx->WrkD[i];

        if (ctx->CCTyp > 0)
            exphess(ctx);
    }

    if (ctx->PMMPParDef == 1)                /* create matrix for parameters */
        mp_putpar(ctx, ctx->NParm,ctx->Par);

    if (ctx->PMProtFDef) {
        prvec(ctx, "Parameter vector",ctx->NParm,ctx->Par);

        if (ctx->LConv >= 0 && ctx->CCTyp > 0)
            prhess(ctx, "Covariance matrix",ctx->NParm);

        fprintf(ctx->PMProtFd,"\nEnd of algorithm.\n\n");
    }
    if (cov == 0 || ctx->CCTyp == 0) {
        err = 0;
        goto FFMFin;
    }

    /* scale cov matrix for nonlinear regression */

    n = ctx->NOC - (ctx->NParm - ctx->NCon1);

    if (typ == 2 && ctx->CCTyp == 2 && ctx->LConv >= 0 && n > 0) {

        tmp = 2.0 * ctx->FMax / (double)n;

        for (j = 1; j <= ctx->NParm; ++j)
            ctx->Diag[j] *= tmp;

        for (j = 1; j <= ctx->HSiz; ++j)
            ctx->Hess[j] *= tmp;
    }
#ifdef TDA_R_PACKAGE
    /* the block below is where the covariance matrix is read out of
       Hess/Diag; under R it is entered whether or not a file or matrix
       was asked for, so that the ml.vcov export at its end is always
       made -- every write inside keeps its own flag */
    if (1) {
#else
    if (ctx->PMPPFDef || ctx->PMCovFDef || ctx->PMMPCovDef) {
#endif

        if (ctx->PMMPCovDef == 1) {              /* allocate matrix for cov matrix */
            if (mp_alloc(ctx, 3,ctx->NParm,ctx->NParm))
                ctx->PMMPCovDef = -1;
        }

        for (i = 1; i <= ctx->NParm; ++i) {

            if (ctx->PMPPFDef)
                rt_fprintf_d(ctx, ctx->PMPPFd,ctx->PMTFmtS,ctx->Par[i]);

            if (ctx->LConv >= 0) {

                if (ctx->PMPPFDef) {
                    tmp = ctx->Diag[i];
                    if (tmp > 0.0)
                        tmp = sqrt(tmp);
                    else
                        tmp = 0.0;
                    rt_fprintf_d(ctx, ctx->PMPPFd,ctx->PMTFmtS,tmp);
                }
                if ((ctx->PMCovFDef || ctx->PMMPCovDef == 1) && ctx->CCTyp > 0) {

                    for (j = 1; j <= ctx->NParm; ++j) {
                        if (j == i)
                            tmp = ctx->Diag[i];
                        else if (j < i)
                            tmp = ctx->Hess[(i - 1) * (i - 2) / 2 + j];
                        else
                            tmp = ctx->Hess[(j - 1) * (j - 2) / 2 + i];

                        if (ctx->PMCovFDef)
                            rt_fprintf_d(ctx, ctx->PMCovFd,ctx->PMMFmtS,tmp);

                        if (ctx->PMMPCovDef == 1)
                            ctx->MatVal[ctx->MPCovIdx][(i - 1) * ctx->NParm + j] = tmp;
                    }
                    if (ctx->PMCovFDef)
                        fprintf(ctx->PMCovFd,"\n");
                }
            }
            if (ctx->PMPPFDef)
                fprintf(ctx->PMPPFd,"\n");
        }
        if (ctx->PMPPFDef)
            ctx->PMPPWFlg = 1;
        if (ctx->LConv >= 0 && ctx->PMCovFDef && ctx->CCTyp > 0)
            ctx->PMCovWFlg = 1;
#ifdef TDA_R_PACKAGE
        /* hand the same covariance matrix the pcov= file gets directly
           to R as well (CONTRIBUTING.md); text and file output
           above are untouched */
        if (ctx->LConv >= 0 && ctx->CCTyp > 0 && ctx->NParm > 0) {
            int ei, ej;
            double etmp;
            double *ev = (double *)malloc((size_t)ctx->NParm *
                                          (size_t)ctx->NParm *
                                          sizeof(double));
            if (ev != NULL) {
                for (ei = 1; ei <= ctx->NParm; ++ei) {
                    for (ej = 1; ej <= ctx->NParm; ++ej) {
                        if (ej == ei)
                            etmp = ctx->Diag[ei];
                        else if (ej < ei)
                            etmp = ctx->Hess[(ei - 1) * (ei - 2) / 2 + ej];
                        else
                            etmp = ctx->Hess[(ej - 1) * (ej - 2) / 2 + ei];
                        ev[(ej - 1) * ctx->NParm + (ei - 1)] = etmp;
                    }
                }
                tda_export_mat(ctx, "ml.vcov", ev, ctx->NParm, ctx->NParm);
                free(ev);
            }
        }
#endif
    }
    err = 0;

FFMFin:
    if (msg) {
        if (err < 0)
            printf1(ctx, "\nError: insufficient memory for function minimization.\n");
        else if (err > 0) {
            printf1(ctx, "\nError in function evaluation.\nCannot continue.\n");
            prn_npflg(ctx);
        }
    }
    if (ctx->Diag1A) {
        free((char *)ctx->Diag1);
        memrq(ctx, -ctx->Diag1A,sizeof(double));
        ctx->Diag1A = 0;
    }
    if (ctx->Hess1A) {
        free((char *)ctx->Hess1);
        memrq(ctx, -ctx->Hess1A,sizeof(double));
        ctx->Hess1A = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ffmin1                                                                  */
/*      Function minimization with a direct search method.                  */
/*      References:                                                         */
/*      A. F. Kaupe, Direct Search, Algorithm 178, Communications of the    */
/*      ACM 6 (1963), 313                                                   */
/*      M. Bell, M.C. Pike, Comm. ACM 9 (1966), 684                         */
/*      R. de Vogelaere, Comm. ACM 11 (1968), 498                           */
/*      F.K. Tomlin, L.B. Smith, Comm. ACM 1969, 637                        */
/*      L.B. Smith, Comm. ACM 1969, 638                                     */
/*                                                                          */
/*      The current implementation corresponds to the proposals of Bell     */
/*      and Pike with modifications proposed by Tomlin & Smith and Smith.   */
/*                                                                          */
/*      Return:     0 if successful.                                        */
/*                 -1 if insufficient memory.                               */
/*                 +1 if error in function evaluation.                      */

int ffmin1(TDAContext *ctx, int mod)
{
    register int k;
    int err;
    double delta,fm,fm1,theta,tmp,*step;

    err = 0;

    if (ctx->SILENTFlg < 2)
        printfe(ctx, "\n  Iter    Function Value       Step Length     Par Change   FCall\n");

    if (ctx->PMProtFDef == 1)
        fprintf(ctx->PMProtFd,"\n  Iter    Function Value       Step Length     Par Change   FCall\n");

    delta = ctx->SLen;
    if (!(step   = (double *) calloc((size_t)(ctx->NParm + 1),sizeof(double))))
        return(-1);
    memrq(ctx, ctx->NParm + 1,sizeof(double));

    for (k = 1; k <= ctx->NParm1; ++k) {
        tmp = fabs(ctx->Par[k]);
        if (tmp < ctx->EPSI1)
            step[k] = delta;
        else
            step[k] = delta * tmp;
    }
    fn(ctx, ctx->Par,mod,0,1,0,&err);
    if (err)
        goto FM1Fin;

    ctx->FMin = ctx->FTmp;

FM1L1:
    fm = ctx->FMin;
    for (k = 1; k <= ctx->NParm1; ++k)
        ctx->Par1[k] = ctx->Par[k];

    for (k = 1; k <= ctx->NParm1; ++k) {
        ctx->Par1[k] += step[k];
        fn(ctx, ctx->Par1,mod,0,1,0,&err);
        if (err)
            goto FM1Fin;

        fm1 = ctx->FTmp;

        if (fm1 < fm)
            fm = fm1;
        else {
            step[k] = -step[k];
            ctx->Par1[k] += 2.0 * step[k];

            fn(ctx, ctx->Par1,mod,0,1,0,&err);
            if (err)
                goto FM1Fin;

            fm1 = ctx->FTmp;

            if (fm1 < fm)
                fm = fm1;
            else
                ctx->Par1[k] -= step[k];
        }
    }
    if (fm < ctx->FMin) {

FM1L2:
        for (k = 1; k <= ctx->NParm1; ++k) {
            if ((ctx->Par1[k] > ctx->Par[k] && step[k] <  0.0) ||
                        (ctx->Par1[k] <= ctx->Par[k] && step[k] >= 0.0))
                step[k] = -step[k];
            theta = ctx->Par[k];
            ctx->Par[k] = ctx->Par1[k];
            ctx->Par1[k] = 2.0 * ctx->Par1[k] - theta;
        }
        ctx->FMin = fm;
        fn(ctx, ctx->Par1,mod,0,1,0,&err);
        if (err)
            goto FM1Fin;

        fm = fm1 = ctx->FTmp;

        for (k = 1; k <= ctx->NParm1; ++k) {
            ctx->Par1[k] += step[k];

            fn(ctx, ctx->Par1,mod,0,1,0,&err);
            if (err)
                goto FM1Fin;

            fm1 = ctx->FTmp;

            if (fm1 < fm)
                fm = fm1;
            else {
                step[k] = -step[k];
                ctx->Par1[k] += 2.0 * step[k];

                fn(ctx, ctx->Par1,mod,0,1,0,&err);
                if (err)
                    goto FM1Fin;

                fm1 = ctx->FTmp;
                if (fm1 < fm)
                    fm = fm1;
                else
                    ctx->Par1[k] -= step[k];
            }
        }
        if (fm >= ctx->FMin)
            goto FM1L1;

        for (k = 1; k <= ctx->NParm1; ++k) {
            if (fabs(ctx->Par1[k] - ctx->Par[k]) > 0.5 * fabs(step[k]))
                goto FM1L2;
        }
    }
    ctx->Iter++;

    if (ctx->SILENTFlg < 2) {
        printfe(ctx, "  %3d  %20.13e %17.10e ",ctx->Iter,ctx->FMin,delta);
        printfe(ctx, "        --   %6d\n",ctx->NumF);
    }
    if (ctx->PMProtFDef) {
        if (ctx->PMProtFDef > 1)
            fprintf(ctx->PMProtFd,"\nITER      Function Value       Step Length     Par Change   FCall\n");
        fprintf(ctx->PMProtFd,"  %3d  %20.13e %17.10e ",ctx->Iter,ctx->FMin,delta);
        fprintf(ctx->PMProtFd,"        --   %6d\n",ctx->NumF);
    }

    if (delta > ctx->TOLS) {
        if (ctx->Iter >= ctx->MxIter) {
            ctx->LConv = 2;
            goto FM1Fin;
        }
        delta *= ctx->SRed;
        for (k = 1; k <= ctx->NParm1; ++k)
            step[k] *= ctx->SRed;
        goto FM1L1;
    }
    else
        ctx->LConv = 0;

FM1Fin:
    if (err)
        err = 1;

    free((char *)step);
    memrq(ctx, -ctx->NParm - 1,sizeof(double));
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ffmin2                                                                  */
/*      Function minimization with a Simplex algorithm.                     */
/*                                                                          */
/*      Ref.: R. O'Neill, Function Minimization using a Simplex Procedure,  */
/*      Algorithm AS 47, Applied Statistics 20 (1971), 338 - 345.           */
/*      Nelder and Mead, Computer Journal 7 (1965), 308-313.                */
/*                                                                          */
/*      Remarks:                                                            */
/*      1. J.M. Chambers, J.E. Ertel, Applied Statistics 23 (1974),250      */
/*      2. R. O'Neill, Applied Statistics 23 (1974), 252                    */
/*      3. P.R. Benyon, Applied Statistics 25 (1976), 97                    */
/*      4. I.D. Hill, Applied Statistics 27 (1978), 380                     */
/*                                                                          */
/*      Return:     0 if successful.                                        */
/*                 -1 if insufficient memory.                               */
/*                 +1 if error in function evaluation.                      */

int ffmin2(TDAContext *ctx, int mod)
{
    register int i,j,l;
    int err,rr,nn,ilo,ihi,ccount;
    int stepa,pstara,p2stara,pbara,ya,pa;
    double *step,*p,*pstar,*p2star,*pbar,*y;
    double dn,dnn,z,sum,summ,ylo,rcoeff,ystar,ecoeff,y2star;
    double ccoeff,del,ynewlo;
    double reduc = 0.001;

    stepa = pstara = p2stara = pbara = ya = pa = 0;

    if (ctx->SILENTFlg < 2)
        printfe(ctx, "\n  Iter    Function Value         Variance      Par Change   FCall\n");

    if (ctx->PMProtFDef == 1)
        fprintf(ctx->PMProtFd,"\n  Iter    Function Value         Variance      Par Change   FCall\n");

    err = -1;
    if (!(step = (double *) calloc((size_t)(ctx->NParm + 1),sizeof(double))))
        goto FM2End;
    memrq(ctx, ctx->NParm + 1,sizeof(double));
    stepa = ctx->NParm + 1;

    if (!(pstar = (double *) calloc((size_t)(ctx->NParm + 1),sizeof(double))))
        goto FM2End;
    memrq(ctx, ctx->NParm + 1,sizeof(double));
    pstara = ctx->NParm + 1;

    if (!(p2star = (double *) calloc((size_t)(ctx->NParm + 1),sizeof(double))))
        goto FM2End;
    memrq(ctx, ctx->NParm + 1,sizeof(double));
    p2stara = ctx->NParm + 1;

    if (!(pbar = (double *) calloc((size_t)(ctx->NParm + 1),sizeof(double))))
        goto FM2End;
    memrq(ctx, ctx->NParm + 1,sizeof(double));
    pbara = ctx->NParm + 1;

    if (!(y = (double *) calloc((size_t)(ctx->NParm + 2),sizeof(double))))
        goto FM2End;
    memrq(ctx, ctx->NParm + 2,sizeof(double));
    ya = ctx->NParm + 2;

    if (!(p = (double *) calloc((size_t)((ctx->NParm + 1)) * (size_t)(ctx->NParm) + 1,sizeof(double))))
        goto FM2End;
    pa = (ctx->NParm + 1) * ctx->NParm + 1;
    memrq(ctx, pa,sizeof(double));

    err = 0;

    for (i = 0; i <= ctx->NParm1; ++i)
        step[i] = ctx->SLen;

    rcoeff = 1.0;
    ecoeff = 2.0;
    ccoeff = 0.5;
    ccount = ctx->CCheck;
    dn = ctx->NParm1;
    nn = ctx->NParm1 + 1;
    dnn = nn;
    del = 1.0;

FM2Res:
    for (i = 1; i <= ctx->NParm1; ++i)
        p[(i - 1) * nn + nn] = ctx->Par[i];

    fn(ctx, ctx->Par,mod,0,1,0,&rr);
    if (rr) {
        err = 1;
        goto FM2End;
    }
    z = ctx->FTmp;

    y[nn] = z;
    sum = z;
    summ = z * z;
    for (j = 1; j <= ctx->NParm1; ++j) {
        ctx->Par[j] += step[j] * del;

        for (i = 1; i <= ctx->NParm1; ++i)
            p[(i - 1) * nn + j] = ctx->Par[i];

        fn(ctx, ctx->Par,mod,0,1,0,&rr);
        if (rr) {
            err = 1;
            goto FM2End;
        }
        z = ctx->FTmp;

        y[j] = z;
        sum += z;
        summ += z * z;
        ctx->Par[j] -= step[j] * del;
    }

    while (1) {
        ylo = y[1];
        ynewlo = ylo;
        ilo = ihi = 1;

        for (i = 2; i <= nn; ++i) {
            if (y[i] < ylo) {
                ylo = y[i];
                ilo = i;
            }
            if (y[i] > ynewlo) {
                ynewlo = y[i];
                ihi = i;
            }
        }
        sum -= ynewlo;
        summ -= ynewlo * ynewlo;

        for (i = 1; i <= ctx->NParm1; ++i) {
            z = 0.0;
            for (j = 1; j <= nn; ++j)
                z += p[(i - 1) * nn + j];

            z -= p[(i - 1) * nn + ihi];
            pbar[i] = z / dn;
        }
        for (i = 1; i <= ctx->NParm1; ++i)
           pstar[i] = (1.0 + rcoeff) * pbar[i] - rcoeff * p[(i - 1) * nn + ihi];

        fn(ctx, pstar,mod,0,1,0,&rr);
        if (rr) {
            err = 1;
            goto FM2End;
        }

        ystar = ctx->FTmp;
        if (ystar >= ylo)
            goto FM2L12;

        for (i = 1; i <= ctx->NParm1; ++i)
            p2star[i] = ecoeff * pstar[i] + (1.0 - ecoeff) * pbar[i];

        fn(ctx, p2star,mod,0,1,0,&rr);
        if (rr) {
            err = 1;
            goto FM2End;
        }
        y2star = ctx->FTmp;
        if (y2star >= ystar) {
            for (i = 1; i <= ctx->NParm1; ++i)
                p[(i - 1) * nn + ihi] = pstar[i];
            y[ihi] = ystar;
            sum += y[ihi];
            summ += y[ihi] * y[ihi];
            goto FM2L901;
        }

FM2L10:
        for (i = 1; i <= ctx->NParm1; ++i)
            p[(i - 1) * nn + ihi] = p2star[i];
        y[ihi] = y2star;
        sum += y[ihi];
        summ += y[ihi] * y[ihi];
        goto FM2L901;

FM2L12:
        l = 0;
        for (i = 1; i <= nn; ++i) {
            if (y[i] > ystar)
                l++;
        }
        if (l > 1) {
            for (i = 1; i <= ctx->NParm1; ++i)
                p[(i - 1) * nn + ihi] = pstar[i];
            y[ihi] = ystar;
            sum += y[ihi];
            summ += y[ihi] * y[ihi];
        }
        else {
            if (l != 0) {
                for (i = 1; i <= ctx->NParm1; ++i)
                    p[(i - 1) * nn + ihi] = pstar[i];
                y[ihi] = ystar;
            }
            for (i = 1; i <= ctx->NParm1; ++i)
                p2star[i] = ccoeff * p[(i - 1) * nn + ihi] + (1.0 - ccoeff)
                                                                * pbar[i];
            fn(ctx, p2star,mod,0,1,0,&rr);
            if (rr) {
                err = 1;
                goto FM2End;
            }
            y2star = ctx->FTmp;
            if (y2star <= y[ihi])
                goto FM2L10;

            sum = summ = 0.0;
            for (j = 1; j <= nn; ++j) {
                for (i = 1; i <= ctx->NParm1; ++i) {
                    p[(i - 1) * nn + j] =
                           (p[(i - 1) * nn + j] + p[(i - 1) * nn + ilo]) * 0.5;
                    ctx->Par1[i] = p[(i - 1) * nn + j];
                }
                fn(ctx, ctx->Par1,mod,0,1,0,&rr);
                if (rr) {
                    err = 1;
                    goto FM2End;
                }
                y[j] = ctx->FTmp;

                sum += y[j];
                summ += y[j] * y[j];
            }
        }
FM2L901:
        if (--ccount == 0) {
            ctx->Iter++;
            ccount = ctx->CCheck;
            ctx->CValV = fabs((summ - (sum * sum) / dnn) / dn);

            if (ctx->SILENTFlg < 2) {
                printfe(ctx, "  %3d  %20.13e %17.10e ",ctx->Iter,ctx->FTmp,ctx->CValV);
                printfe(ctx, "        --   %6d\n",ctx->NumF);
            }
            if (ctx->PMProtFDef) {
                if (ctx->PMProtFDef > 1)
                    fprintf(ctx->PMProtFd,"\nITER      Function Value         Variance      Par Change   FCall\n");
                fprintf(ctx->PMProtFd,"  %3d  %20.13e %17.10e ",ctx->Iter,ctx->FTmp,ctx->CValV);
                fprintf(ctx->PMProtFd,"        --   %6d\n",ctx->NumF);
            }
            if (ctx->Iter >= ctx->MxIter || ctx->CValV <= ctx->TOLV)
                break;
        }
    }
    if (y[ihi] > y[ilo])
        ihi = ilo;

    for (i = 1; i <= ctx->NParm1; ++i)
        ctx->Par1[i] = p[(i - 1) * nn + ihi];

    ctx->FMin = ynewlo = y[ihi];

    if (ctx->Iter >= ctx->MxIter) {
        ctx->LConv = 2;
        goto FM2Fin;
    }
    for (i = 1; i <= ctx->NParm1; ++i) {     /* check final result */
        del = step[i] * reduc;
        ctx->Par1[i] += del;

        fn(ctx, ctx->Par1,mod,0,1,0,&rr);
        if (rr) {
            err = 1;
            goto FM2End;
        }

        z = ctx->FTmp;
        if (ynewlo - z > ctx->TOLF)
            goto FM2Nxt;

        ctx->Par1[i] = ctx->Par1[i] - del - del;

        fn(ctx, ctx->Par1,mod,0,1,0,&rr);
        if (rr) {
            err = 1;
            goto FM2End;
        }

        z = ctx->FTmp;
        if (ynewlo - z > ctx->TOLF)
            goto FM2Nxt;

        ctx->Par1[i] += del;
    }
    ctx->FMin = ynewlo;
    ctx->LConv = 0;

FM2Fin:
    for (i = 1; i <= ctx->NParm1; ++i)
        ctx->Par[i] = ctx->Par1[i];
    goto FM2End;

FM2Nxt:
    for (i = 1; i <= ctx->NParm1; ++i)
        ctx->Par[i] = ctx->Par1[i];
    del = reduc;
    if (ctx->SILENTFlg < 2)
        printfe(ctx, "  Restart\n");
    if (ctx->PMProtFDef)
        fprintf(ctx->PMProtFd,"RESTART.\n");

    goto FM2Res;

FM2End:
    if (pa > 0) {
        free((char *)p);
        memrq(ctx, -pa,sizeof(double));
    }
    if (ya > 0) {
        free((char *)y);
        memrq(ctx, -ya,sizeof(double));
    }
    if (pbara > 0) {
        free((char *)pbar);
        memrq(ctx, -pbara,sizeof(double));
    }
    if (p2stara > 0) {
        free((char *)p2star);
        memrq(ctx, -p2stara,sizeof(double));
    }
    if (pstara > 0) {
        free((char *)pstar);
        memrq(ctx, -pstara,sizeof(double));
    }
    if (stepa > 0) {
        free((char *)step);
        memrq(ctx, -stepa,sizeof(double));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ffmin3  Function minimization with conjugate gradient or BFGS           */
/*          Adapted and translated to C from: Algorithm 500 of Collected    */
/*          Algorithms from ACM., Trans. Math. Software 6 (1980), 618 - 22  */
/*          Written by D.F. Shanno and K.H. Phua.                           */
/*                                                                          */
/*          Note: most comments have been stripped. See the original        */
/*          Fortran source for a full documentation.                        */
/*                                                                          */
/*          Hess[] is used as a working area. It must have dimension of     */
/*          at least 5 * NParm + 2 if MINA = 3 (conjugate gradients)        */
/*          at least NParm * (NParm + 7) / 2 if MINA = 4 (BFGS).            */
/*                                                                          */
/*      Return:     0 if successful.                                        */
/*                 -1 if insufficient memory.                               */
/*                 +1 if error in function evaluation.                      */

int ffmin3(TDAContext *ctx, int mod)
{
    register int i = 0,j = 0,ii = 0,ij = 0;
    int nx = 0,ng = 0,nry = 0,nrd = 0,ncons = 0,ncons1 = 0,ncons2 = 0,nrst = 0,rsw = 0,ifun = 0,ncalls = 0;
    int err = 0,it1 = 0,nxpi = 0,ngpi = 0,nrdpi = 0,nrypi = 0,ngpj = 0;
    double tmp = 0.0,fm = 0.0,fp = 0.0,dal = 0.0,dg1 = 0.0,gsq = 0.0,alpha = 0.0,at = 0.0,ap = 0.0,dp = 0.0,dg = 0.0,rtst = 0.0,step = 0.0,u1 = 0.0,u2 = 0.0,u3 = 0.0,u4 = 0.0;

    if (ctx->SILENTFlg < 2)
        printfe(ctx, "\n  Iter    Function Value     Norm of Gradient  Par Change   FCall\n");

    if (ctx->PMProtFDef == 1)
        fprintf(ctx->PMProtFd,"\n  Iter    Function Value     Norm of Gradient  Par Change   FCall\n");

    alpha = 1.0;
    err = ifun = 0;
    nx = ctx->NParm1;
    ng = nx + ctx->NParm1;

    if (ctx->MINA == 3) {       /* conjugate gradient */
        nry = ng + ctx->NParm1;
        nrd = nry + ctx->NParm1;
        ncons = 5 * ctx->NParm1;
        ncons1 = ncons + 1;
        ncons2 = ncons + 2;
    }
    else
        ncons = 3 * ctx->NParm1;

FFL20:
    fn(ctx, ctx->Par,mod,1,1,0,&err);              /* initial function evaluation */
    if (err)
        goto FM3Fin;

    ctx->FMin = ctx->FTmp;
    ifun++;

    nrst = ctx->NParm1;
    rsw = 1;
    dg1 = 0.0;
    for (i = 1; i <= ctx->NParm1; ++i) {
        ctx->Hess[i] = -ctx->Grad[i];
        dg1 -= ctx->Grad[i] * ctx->Grad[i];
    }
    gsq = -dg1;

    dg = dg1;   /* inserted */

    ctx->CValG = 0.0;                     /* calculate norm of gradient */
    for (i = 1; i <= ctx->NParm1; ++i) {
        tmp = ctx->Grad[i];
        ctx->CValG += tmp * tmp;
    }
    ctx->CValG = sqrt(ctx->CValG);
    if (ctx->CValG <= ctx->TOLG)      /* check if initial point is minimizer */
        goto FM3Fin;

    /* entry point of major iteration loop */

    while (1) {

        if (++ctx->Iter > ctx->MxIter) {
            ctx->Iter--;
            ctx->LConv = 2;      /* exceeded max number of iterations (func calls) */
            goto FM3Fin;
        }
        if (ctx->SILENTFlg < 2) {
            printfe(ctx, "  %3d  %20.13e %17.10e ",ctx->Iter,ctx->FMin,ctx->CValG);
            printfe(ctx, "        --   %6d (%d,%d)\n",ctx->NumF,ctx->NumFG,ctx->NumFH);
        }
        if (ctx->PMProtFDef) {
            if (ctx->PMProtFDef > 1)
                fprintf(ctx->PMProtFd,"\nITER      Function Value     Norm of Gradient  Par Change   FCall\n");
            fprintf(ctx->PMProtFd,"  %3d  %20.13e %17.10e ",ctx->Iter,ctx->FMin,ctx->CValG);
            fprintf(ctx->PMProtFd,"        --   %6d (%d,%d)\n",ctx->NumF,ctx->NumFG,ctx->NumFH);

            if (ctx->PMProtFDef > 1) {
                prvec(ctx, "\nParameter",ctx->NParm1,ctx->Par);
                prvec(ctx, "Gradient",ctx->NParm1,ctx->Grad);
            }
        }
        fm = ctx->FMin;
        ncalls = ifun;

        /* begin linear search */

        alpha *= dg / dg1;

        if (ctx->PMProtFDef > 1)
            fprintf(ctx->PMProtFd,"Begin line search: alpha = %g dg=%g dg1=%g\n",
                                               alpha,dg,dg1);
        if (nrst == 1 || ctx->MINA == 4)
            alpha = 1.0;
        if (rsw)
            alpha = 1.0 / sqrt(gsq);
        ap = 0.0;
        fp = fm;
        dg = dp = dg1;

        step = 0.0;
        for (i = 1; i <= ctx->NParm1; ++i) {
            step += ctx->Hess[i] * ctx->Hess[i];
            nxpi = nx + i;
            ngpi = ng + i;
            ctx->Hess[nxpi] = ctx->Par[i];
            ctx->Hess[ngpi] = ctx->Grad[i];
        }
        step = sqrt(step);

        /*  begin of linear search iterations */

        it1 = 0;
FFL80:

        if (ctx->PMProtFDef > 1)
            fprintf(ctx->PMProtFd,"Line search with alpha=%g step=%g\n",alpha,step);

        if (++it1 > ctx->MxIt1) {
            ctx->LConv = -5;     /* linear search failed */
            goto FM3Fin;
        }
        if (alpha * step <= ctx->SMin) {
            if (!rsw)
                goto FFL20;

            if (ctx->PMProtFDef > 1)
                fprintf(ctx->PMProtFd,"Line search failed (smin=%g alpha=%g step=%g).\n",ctx->SMin,alpha,step);

            ctx->LConv = -2;     /* linear search failed */
            goto FM3Fin;
        }
        for (i = 1; i <= ctx->NParm1; ++i) {
            nxpi = nx + i;
            ctx->Par[i] = ctx->Hess[nxpi] + alpha * ctx->Hess[i];
        }
        fn(ctx, ctx->Par,mod,1,1,0,&err);      /* new function evaluation */
        if (err)
            goto FM3Fin;

        ctx->FMin = ctx->FTmp;
        ifun++;

        dal=0.0;
        for (i = 1; i <= ctx->NParm1; ++i) {
            dal += ctx->Grad[i] * ctx->Hess[i];
        }
        if (ctx->FMin > fm && dal < 0.0) {
            alpha /= 3.0;
            ap = 0.0;
            fp = fm;
            dp = dg;
            goto FFL80;
        }
        if (ctx->FMin > (fm + 0.0001 * alpha * dg) || fabs(dal / dg) > 0.9 ||
           ((ifun - ncalls) <= 1 && fabs(dal / dg) > ctx->EPSI1 && ctx->MINA == 3)) {

            u1 = dp + dal - 3.0 * (fp - ctx->FMin) / (ap - alpha);
            u2 = u1 * u1 - dp * dal;
            if (u2 <= 0.0)
                u2 = 0.0;
            else
                u2 = sqrt(u2);
            at = alpha - (alpha - ap) * (dal + u2 - u1) / (dal - dp + 2.0 * u2);

            if ((dal / dp) <= 0.0) {
                if (at < (1.01 * dmin(ctx, alpha,ap)) || at > (0.99 * dmax(ctx, alpha,ap)))
                    at = (alpha + ap) / 2.0;
            }
            else {
                if ((dal > 0.0 && 0.0 < at && at < (0.99 * dmin(ctx, alpha,ap))) ||
                    (dal <= 0.0 && at > (1.01 * dmax(ctx, alpha,ap)))) {
                    ;
                }
                else {
                    if (dal <= 0.0)
                        at = 2.0 * dmax(ctx, alpha,ap);
                    else
                        at = dmin(ctx, alpha,ap) / 2.0;
                }
            }
            ap = alpha;
            fp = ctx->FMin;
            dp = dal;
            alpha = at;
            goto FFL80;
        }

        /*  line search was successful. test convergence */

        gsq = 0.0;
        for (i = 1; i <= ctx->NParm1; ++i) {
            gsq += ctx->Grad[i] * ctx->Grad[i];
        }
        tmp = 0.0;                      /* calculate norm of gradient */
        for (i = 1; i <= ctx->NParm1; ++i)
            tmp += ctx->Grad[i] * ctx->Grad[i];
        ctx->CValG = sqrt(tmp);
        if (ctx->CValG <= ctx->TOLG) {
            ctx->Iter++;
            goto FM3Prn;
        }
        for (i = 1; i <= ctx->NParm1; ++i) {
            ctx->Hess[i] *= alpha;
        }
        if (ctx->MINA == 3) {       /* conjugate gradients */
            rtst = 0.0;
            for (i = 1; i <= ctx->NParm1; ++i) {
                ngpi = ng + i;
                rtst += ctx->Grad[i] * ctx->Hess[ngpi];
            }
            if (fabs(rtst / gsq) > 0.2)
                nrst = ctx->NParm1;

            if (nrst == ctx->NParm1) {
                ctx->Hess[ncons + 1] = 0.0;
                ctx->Hess[ncons + 2] = 0.0;
                for (i = 1; i <= ctx->NParm1; ++i) {
                    nrdpi = nrd + i;
                    nrypi = nry + i;
                    ngpi = ng + i;
                    ctx->Hess[nrypi]  = ctx->Grad[i] - ctx->Hess[ngpi];
                    ctx->Hess[nrdpi]  = ctx->Hess[i];
                    ctx->Hess[ncons1] += ctx->Hess[nrypi] * ctx->Hess[nrypi];
                    ctx->Hess[ncons2] += ctx->Hess[i] * ctx->Hess[nrypi];
                }
            }
            u1 = u2 = 0.0;
            for (i = 1; i <= ctx->NParm1; ++i) {
                nrdpi = nrd + i;
                nrypi = nry + i;
                u1 = u1 - ctx->Hess[nrdpi] * ctx->Grad[i] / ctx->Hess[ncons1];
                u2 = u2 + ctx->Hess[nrdpi] * ctx->Grad[i] * 2.0 / ctx->Hess[ncons2] -
                                               ctx->Hess[nrypi] * ctx->Grad[i] / ctx->Hess[ncons1];
            }
            u3 = ctx->Hess[ncons2] / ctx->Hess[ncons1];
            for (i = 1; i <= ctx->NParm1; ++i) {
                nxpi = nx + i;
                nrdpi = nrd + i;
                nrypi = nry + i;
                ctx->Hess[nxpi] = -u3 * ctx->Grad[i] - u1 * ctx->Hess[nrypi] - u2 * ctx->Hess[nrdpi];
            }
            if (nrst != ctx->NParm1) {
                u1 = u2 = u3 = u4 = 0.0;
                for (i = 1; i <= ctx->NParm1; ++i) {
                    ngpi = ng + i;
                    nrdpi = nrd + i;
                    nrypi = nry + i;
                    u1 = u1 - (ctx->Grad[i] - ctx->Hess[ngpi]) * ctx->Hess[nrdpi] / ctx->Hess[ncons1];
                    u2 = u2 - (ctx->Grad[i] - ctx->Hess[ngpi]) * ctx->Hess[nrypi] / ctx->Hess[ncons1] +
                                      2.0 * ctx->Hess[nrdpi] * (ctx->Grad[i] - ctx->Hess[ngpi]) / ctx->Hess[ncons2];
                    u3 = u3 + ctx->Hess[i] * (ctx->Grad[i] - ctx->Hess[ngpi]);
                }
                step = 0.0;
                for (i = 1; i <= ctx->NParm1; ++i) {
                    ngpi = ng + i;
                    nrdpi = nrd + i;
                    nrypi = nry + i;
                    step = (ctx->Hess[ncons2] / ctx->Hess[ncons1]) * (ctx->Grad[i] - ctx->Hess[ngpi]) +
                                                  u1 * ctx->Hess[nrypi] + u2 * ctx->Hess[nrdpi];
                    u4 = u4 + step * (ctx->Grad[i] - ctx->Hess[ngpi]);
                    ctx->Hess[ngpi] = step;
                }
                u1 = u2 = 0.0;
                for (i = 1; i <= ctx->NParm1; ++i) {
                    u1 = u1 - ctx->Hess[i] * ctx->Grad[i] / u3;
                    ngpi = ng + i;
                    u2 = u2 + (1.0 + u4 / u3) * ctx->Hess[i] * ctx->Grad[i] / u3 -
                                                         ctx->Hess[ngpi] * ctx->Grad[i] / u3;
                }
                for (i = 1; i <= ctx->NParm1; ++i) {
                    ngpi = ng + i;
                    nxpi = nx + i;
                    ctx->Hess[nxpi] = ctx->Hess[nxpi] - u1 * ctx->Hess[ngpi] - u2 * ctx->Hess[i];
                }
            }
            dg1 = 0.0;
            for (i = 1; i <= ctx->NParm1; ++i) {
                nxpi = nx + i;
                ctx->Hess[i] = ctx->Hess[nxpi];
                dg1 += ctx->Hess[i] * ctx->Grad[i];
            }
            if (dg1 > 0.0) {
                ctx->LConv = -4; /* search vector not a descent direction */
                goto FM3Fin;
            }
            if (nrst == ctx->NParm1)
                nrst = 0;
            nrst++;
            rsw = 0;
        }
        else {  /* BFGS */
            u1 = 0.0;
            for (i = 1; i <= ctx->NParm1; ++i) {
                ngpi = ng + i;
                ctx->Hess[ngpi] = ctx->Grad[i] - ctx->Hess[ngpi];
                u1 += ctx->Hess[i] * ctx->Hess[ngpi];
            }
            if (rsw) {
                u2 = 0.0;
                for (i = 1; i <= ctx->NParm1; ++i) {
                    ngpi = ng + i;
                    u2 += ctx->Hess[ngpi] * ctx->Hess[ngpi];
                }
                ij = 1;
                u3 = u1 / u2;
                for (i = 1; i <= ctx->NParm1; ++i) {
                    for (j = i; j <= ctx->NParm1; ++j) {
                        ncons1 = ncons + ij;
                        ctx->Hess[ncons1] = 0.0;
                        if (i == j)
                            ctx->Hess[ncons1] = u3;
                        ij++;
                    }
                    nxpi = nx + i;
                    ngpi = ng + i;
                    ctx->Hess[nxpi] = u3 * ctx->Hess[ngpi];
                }
                u2 *= u3;
            }
            else {
                u2 = 0.0;
                for (i = 1; i <= ctx->NParm1; ++i) {
                    u3 = 0.0;
                    ij = i;
                    if (i != 1) {
                        ii = i - 1;
                        for (j = 1; j <= ii; ++j) {
                            ngpj = ng + j;
                            ncons1 = ncons + ij;
                            u3 += ctx->Hess[ncons1] * ctx->Hess[ngpj];
                            ij = ij + ctx->NParm1 - j;
                        }
                    }
                    for (j = i; j <= ctx->NParm1; ++j) {
                        ncons1 = ncons + ij;
                        ngpj = ng + j;
                        u3 += ctx->Hess[ncons1] * ctx->Hess[ngpj];
                        ij++;
                    }
                    ngpi = ng + i;
                    u2 += u3 * ctx->Hess[ngpi];
                    nxpi = nx + i;
                    ctx->Hess[nxpi] = u3;
                }
            }
            u4 = 1.0 + u2 / u1;
            for (i = 1; i <= ctx->NParm1; ++i) {
                nxpi = nx + i;
                ngpi = ng + i;
                ctx->Hess[ngpi] = u4 * ctx->Hess[i] - ctx->Hess[nxpi];
            }
            ij = 1;

            for (i = 1; i <= ctx->NParm1; ++i) {
                nxpi = nx + i;
                u3 = ctx->Hess[i] / u1;
                u4 = ctx->Hess[nxpi] / u1;
                for (j = i; j <= ctx->NParm1; ++j) {
                    ncons1 = ncons + ij;
                    ngpj = ng + j;
                    ctx->Hess[ncons1] = ctx->Hess[ncons1] + u3 * ctx->Hess[ngpj] - u4 * ctx->Hess[j];
                    ij++;
                }
            }
            dg1 = 0.0;
            for (i = 1; i <= ctx->NParm1; ++i) {
                u3 = 0.0;
                ij = i;
                if (i != 1) {
                    ii = i - 1;
                    for (j = 1; j <= ii; ++j) {
                        ncons1 = ncons + ij;
                        u3 -= ctx->Hess[ncons1] * ctx->Grad[j];
                        ij = ij + ctx->NParm1 - j;
                    }
                }
                for (j = i; j <= ctx->NParm1; ++j) {
                    ncons1 = ncons + ij;
                    u3 -= ctx->Hess[ncons1] * ctx->Grad[j];
                    ij++;
                }
                dg1 += u3 * ctx->Grad[i];
                ctx->Hess[i] = u3;
            }
            if (dg1 > 0.0) {
                ctx->LConv = -4; /* search vector not a descent direction */
                goto FM3Fin;
            }
            rsw = 0;
        }
    }
FM3Prn:
    if (ctx->SILENTFlg < 2) {
        printfe(ctx, "  %3d  %20.13e %17.10e ",ctx->Iter,ctx->FMin,ctx->CValG);
        printfe(ctx, "        --   %6d (%d,%d)\n",ctx->NumF,ctx->NumFG,ctx->NumFH);
    }
    if (ctx->PMProtFDef) {
        if (ctx->PMProtFDef > 1)
            fprintf(ctx->PMProtFd,"\nITER      Function Value     Norm of Gradient  Par Change   FCall\n");
        fprintf(ctx->PMProtFd,"  %3d  %20.13e %17.10e ",ctx->Iter,ctx->FMin,ctx->CValG);
        fprintf(ctx->PMProtFd,"        --   %6d (%d,%d)\n",ctx->NumF,ctx->NumFG,ctx->NumFH);

        if (ctx->PMProtFDef > 1) {
            prvec(ctx, "\nParameter",ctx->NParm1,ctx->Par);
            prvec(ctx, "Gradient",ctx->NParm1,ctx->Grad);
        }
    }

FM3Fin:
    if (err)
        err = 1;
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ffmin5  Function minimization with a modified Newton method.            */
/*          If MINA == 5, modification is with gradient methods, if         */
/*          MINA == 6 use additional information of max eigenvector.        */
/*                                                                          */
/*          Return -1 if insufficient memory, otherwise 0.                  */

int ffmin5(TDAContext *ctx, int mod)
{
    register int i,j,k;
    int err,npflg,ssflg;
    double step,test,test1,test2,lval,tmp,emin,*ev,*evec;

    err = 0;
    lval = 0.0;

    if (ctx->SILENTFlg < 2)
        printfe(ctx, "\n  Iter    Function Value     Norm of Gradient  Par Change   FCall\n");

    if (ctx->PMProtFDef == 1)
        fprintf(ctx->PMProtFd,"\n  Iter    Function Value     Norm of Gradient  Par Change   FCall\n");

    ctx->LConv = 2;      /* return status of this function */
    npflg = 0;      /* set if hessian not positive definite and we use the  */
                    /* modified algorithm.                                  */

                             /* first call of function evaluation, including   */
    fn(ctx, ctx->Par,mod,2,1,0,&err);  /* derivatives in Grad, Hess, Diag. The function  */
    if (err)
        goto FM5Fin;

    ctx->FMin = ctx->FTmp;        /* value is returned in FTmp.                     */

    while (++ctx->Iter <= ctx->MxIter) {   /* Main loop over iterations */
        /******************
        if (FTmp <= 0.0)
            NPFlgs[3] += 1;
        *********************/
        tmp = 0.0;                      /* calculate norm of gradient */
        for (j = 1; j <= ctx->NParm1; ++j)
            tmp += ctx->Grad[j] * ctx->Grad[j];
        ctx->CValG = sqrt(tmp);

        if (ctx->Iter > 1) {
            ctx->CValF = fabs(ctx->FMin - lval);
            if ((lval) != 0.0)
                ctx->CValF /= lval;

            ctx->CValP = 0.0;
            for (j = 1; j <= ctx->NParm1; ++j) {
                tmp = fabs(ctx->Par[j] - ctx->Par1[j]);
                if ((ctx->Par[j]) != 0.0)
                    tmp /= fabs(ctx->Par[j]);
                if (ctx->CValP < tmp)
                    ctx->CValP = tmp;
            }
            if (ctx->PMProtFDef > 1) {
                prval(ctx, "Change in function value",ctx->CValF);
                prval(ctx, "Change in parameter estimates",ctx->CValP);
            }
        }
        if (ctx->SILENTFlg < 2)
            printfe(ctx, "  %3d  %20.13e %17.10e ",ctx->Iter,ctx->FMin,ctx->CValG);

        if (ctx->PMProtFDef) {
            if (ctx->PMProtFDef > 1)
                fprintf(ctx->PMProtFd,"\nITER      Function Value     Norm of Gradient  Par Change   FCall\n");
            fprintf(ctx->PMProtFd,"  %3d  %20.13e %17.10e ",ctx->Iter,ctx->FMin,ctx->CValG);
        }
        if (ctx->Iter == 1) {
            if (ctx->SILENTFlg < 2)
                printfe(ctx, "         --  %6d (%d,%d)",ctx->NumF,ctx->NumFG,ctx->NumFH);
            if (ctx->PMProtFDef)
                fprintf(ctx->PMProtFd,"         --  %6d (%d,%d)\n",ctx->NumF,ctx->NumFG,ctx->NumFH);
        }
        else {
            if (ctx->SILENTFlg < 2)
                printfe(ctx, "%11.4e  %6d (%d,%d)",ctx->CValP,ctx->NumF,ctx->NumFG,ctx->NumFH);
            if (ctx->PMProtFDef)
                fprintf(ctx->PMProtFd,"%11.4e  %6d (%d,%d)\n",ctx->CValP,ctx->NumF,ctx->NumFG,ctx->NumFH);
        }
        if (ctx->PMProtFDef > 1) {
            prvec(ctx, "\nParameter",ctx->NParm1,ctx->Par);
            prvec(ctx, "Gradient",ctx->NParm1,ctx->Grad);
            prhess(ctx, "Hessian",ctx->NParm1);
        }

        /*  Check for convergence. If Crite > 1 begin with second iteration */

        if (ctx->Crite == 1 || ctx->Iter > 1) {

            switch (ctx->Crite) {
                case 1:     if (fabs(ctx->CValG) <= ctx->TOLG) ctx->LConv = 0;
                            break;
                case 2:     if (fabs(ctx->CValF) <= ctx->TOLF) ctx->LConv = 0;
                            break;
                case 3:     if (fabs(ctx->CValP) <= ctx->TOLP) ctx->LConv = 0;
                            break;
                default:    break;
            }
            if (!ctx->CritFlg && ctx->Iter > 1 && (fabs(ctx->CValG) <= ctx->TOLG || fabs(ctx->CValF) <= ctx->TOLF))
                ctx->LConv = 0;

            if (!ctx->LConv)                 /* convergence was reached.         */
                goto FM5Fin;

            /****************************************************************
            if (Iter > 1 && fabs(CValF) <= TOLF) {  check for minimum decrease
                                                    of the function value.
                LConv = 1;
                goto FM5Fin;
            }
            *****************************************************************/
        }
        lval = ctx->FMin;                        /* save actual function value   */
        for (j = 1; j <= ctx->NParm1; ++j)       /* and parameters.              */
            ctx->Par1[j] = ctx->Par[j];

        if (npflg) {     /* we need a copy of the hessian in this case */

            for (j = 1; j <= ctx->NParm1; ++j) {
                evec[(j - 1) * ctx->NParm1 + j] = ctx->Diag[j];
                for (k = 1; k < j; ++k)
                    evec[(j - 1) * ctx->NParm1 + k] =
                    evec[(k - 1) * ctx->NParm1 + j] = ctx->Hess[(j - 1) * (j - 2) / 2 + k];
            }
        }

        /*  Try to calculate the Newton search direction by a Cholesky      */
        /*  decomposition of the hessian used to solve Hess * Srch = -Grad. */
        /*  This is done by chinv(). If the return is zero, is was success- */
        /*  full and Srch contains the Newton search vector. Otherwise the  */
        /*  hessian is not positive definite.                               */

        if (!chinv(ctx, 1)) {

            /*  The hessian is positive definite. Srch contains now the     */
            /*  Newton search vector. If constraints, Srch is the search    */
            /*  vector in the reduced parameter space.                      */

            if (ctx->SILENTFlg < 2)
                printfe(ctx, "\n");
            if (ctx->PMProtFDef > 1)
                prvec(ctx, "Search vector",ctx->NParm1,ctx->Srch);
            ssflg = 1;
        }
        else {  /*  The hessian is not positive definite. */

            ctx->NPFlgs[4] += 1;
            if (ctx->SILENTFlg < 2)
                printfe(ctx, "  (npd)\n");
            if (ctx->PMProtFDef > 1)
                fprintf(ctx->PMProtFd,"Hessian not positive definite.\n");

            /*  First part of the search vector is the negative gradient.   */

            for (j = 1; j <= ctx->NParm1; ++j)
                ctx->Srch[j] = -ctx->Grad[j];

            if (ctx->MINA == 5) {

                if (ctx->PMProtFDef > 1)
                    prvec(ctx, "Search vector",ctx->NParm1,ctx->Srch);

                ssflg = 1;
            }
            else {

                ssflg = 0;

                /*  The second part is the eigenvector of the hessian       */
                /*  associated with the minimal eigenvalue.                 */

                if (!npflg) {   /* we need some additional storage in this case */

                    if (!(ev = (double *)calloc((size_t)(ctx->NParm1 + 1),sizeof(double))))
                        return(-1);

                    if (!(evec = (double *)calloc((size_t)(ctx->NParm1) * (size_t)(ctx->NParm1) + 1,sizeof(double)))) {
                        free((char *)ev);
                        return(-1);
                    }
                    memrq(ctx, ctx->NParm1 + ctx->NParm1 * ctx->NParm1 + 2,sizeof(double));

                    fn(ctx, ctx->Par,mod,2,1,0,&err);
                    if (err)
                        goto FM5Fin;

                    for (j = 1; j <= ctx->NParm1; ++j) {
                        evec[(j - 1) * ctx->NParm1 + j] = ctx->Diag[j];
                        for (k = 1; k < j; ++k)
                            evec[(j - 1) * ctx->NParm1 + k] =
                            evec[(k - 1) * ctx->NParm1 + j] = ctx->Hess[(j - 1) * (j - 2) / 2 + k];
                    }
                    npflg = 1;
                }

                /*  Calculate eigenvalues and eigenvectors of the hessian */

                if (evecf(ctx, ctx->NParm1,ev,evec)) {
                    ctx->LConv = -3;
                    goto FM5Fin;
                }

                /*  determine the minimum eigenvalue emin, and the associated   */
                /*  eigenvector, stored again in evec.                          */

                k = 1;
                emin = ev[1];
                for (j = 2; j <= ctx->NParm1; ++j) {
                    if (emin > ev[j]) {
                        emin = ev[j];
                        k = j;
                    }
                }
                i = 1;
                for (j = 1; j <= ctx->NParm1; ++j)
                    evec[i++] = evec[(j - 1) * ctx->NParm1 + k];

                if (ctx->PMProtFDef > 1) {
                    prval(ctx, "Minimum eigenvalue",emin);
                    prvec(ctx, "Associated eigenvector",ctx->NParm1,evec);
                    fprintf(ctx->PMProtFd,"Step size search.\n");
                }

                /*  Now select the sign so that evec is a nonascent direction   */
                /*  with negative curvature.                                    */

                test1 = 0.0;
                for (j = 1; j <= ctx->NParm1; ++j)
                    test1 += evec[j] * ctx->Grad[j];

                if (test1 > 0.0) {
                    for (j = 1; j <= ctx->NParm1; ++j)
                        evec[j] = -evec[j];
                }

                /*  The step size search is one with the second-order Armijo    */
                /*  condition. Cf. McCormick 1983, p. 135.                      */

                test1 = test2 = 0.0;
                for (j = 1; j <= ctx->NParm1; ++j) {
                    test1 += ctx->Grad[j] * ctx->Srch[j];
                    test2 += evec[j] * evec[j];
                }
                test2 *= emin;
                test = ctx->AMue * (test1 + 0.5 * test2);
                step = 1.0;

                while (1) {
                    for (j = 1; j <= ctx->NParm1; ++j)
                        ctx->Par[j] = ctx->Par1[j] + step * ctx->Srch[j] + sqrt(step) * evec[j];


                    if (ctx->STFlg)          /* depending on STFlg the calculation   */
                        fn(ctx, ctx->Par,mod,0,1,0,&err);      /* includes derivatives or not.         */
                    else
                        fn(ctx, ctx->Par,mod,2,1,0,&err);
                    if (err)
                        goto FM5Fin;

                    if (ctx->PMProtFDef > 1)
                        prval(ctx, "Function value",ctx->FTmp);

                    /*******************
                    if (FTmp <= 0.0)
                        NPFlgs[3] += 1;
                    *********************/
                    if (ctx->FTmp - ctx->FMin <= test * step)
                        break;

                    /* check for minimum change in function value */
                    if (fabs(ctx->FTmp - ctx->FMin) <= ctx->TOLF) {
                        break;
                        /*********** changed
                        LConv = 1;
                        goto FM5Fin;
                        ********************/
                    }
                    if ((step *= 0.5) < ctx->SMin) {    /* check for minimal step size  */
                        ctx->LConv = -2;
                        goto FM5Fin;
                    }
                }
            }
        }
        if (ssflg) {        /* standard line search */

            if (ctx->PMProtFDef > 1)
                fprintf(ctx->PMProtFd,"Step size search.\n");

            /*  Calculation of a suitable step size. We use the first-      */
            /*  order Armijo condition, cf. McCormick 1983, p. 134.         */

            test = 0.0;
            for (j = 1; j <= ctx->NParm1; ++j)       /* test criterion in the    */
                test += ctx->Grad[j] * ctx->Srch[j];      /* reduced space.           */

            test *= ctx->AMue;
            step = 1.0;

            while (1) {
                for (j = 1; j <= ctx->NParm1; ++j)
                    ctx->Par[j] = ctx->Par1[j] + step * ctx->Srch[j];

                if (ctx->STFlg)                    /* depending on STFlg the calculation   */
                    fn(ctx, ctx->Par,mod,0,1,0,&err);   /* includes derivatives or not.         */
                else
                    fn(ctx, ctx->Par,mod,2,1,0,&err);
                if (err)
                    goto FM5Fin;

                if (ctx->PMProtFDef > 1)
                    prval(ctx, "Function value",ctx->FTmp);

                /******************
                if (FTmp <= 0.0)
                    NPFlgs[3] += 1;
                *********************/
                if (ctx->FTmp - ctx->FMin <= test * step)
                    break;

                if (ctx->PMProtFDef > 1) {
                    prval(ctx, "Change in function value",ctx->FTmp - ctx->FMin);
                    prval(ctx, "Armijo test criterion",test * step);
                }

                /* check for minimum change in function value */
                if (fabs(ctx->FTmp - ctx->FMin) <= ctx->TOLF) {
                    break;
                    /************* changed
                    LConv = 1;
                    goto FM5Fin;
                    *********************/
                }
                if ((step *= 0.5) < ctx->SMin) {    /* check for minimal step size  */
                    ctx->LConv = -2;
                    goto FM5Fin;
                }
            }
        }
        if (ctx->STFlg) {                 /* in this case we need an additional function  */
            fn(ctx, ctx->Par,mod,2,1,0,&err);  /* evaluation including derivatives.            */
            if (err)
                goto FM5Fin;
        }
        ctx->FMin = ctx->FTmp;    /* set new function value */

        if (ctx->PMProtFDef > 1)
            prval(ctx, "Step size",step);
    }
    ctx->LConv = 2;                /* exceeded max number of iterations */
    ctx->Iter--;

FM5Fin:
    if (npflg) {
        memrq(ctx, ctx->NParm1 + ctx->NParm1 * ctx->NParm1 + 2,-(int)sizeof(double));
        free((char *)ev);
        free((char *)evec);
    }
    if (err)
        err = 1;
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  exphess   Expand the Hessian matrix (given by Hess and Diag) to NParm   */
/*            by premultiplication with ConQ and by postmultiplication      */
/*            with the transpose of ConQ.                                   */
/*            Use of WrkD and WrkH as temporary storage.                    */

void exphess(TDAContext *ctx)
{
    register int i,j,k,l;
    double tmp;

    for (i = 1; i <= ctx->NParm; ++i) {
        for (j = 1; j <= i; ++j) {
            tmp = 0.0;
            for (k = 1; k <= ctx->NParm1; ++k) {
                for (l = 1; l <= ctx->NParm1; ++l) {
                    if (l == k)
                        tmp += ctx->Diag[k] * ctx->ConQ[(ctx->CIP[i] - 1) * ctx->NParm1 + l] *
                                         ctx->ConQ[(ctx->CIP[j] - 1) * ctx->NParm1 + k];
                    else if (l < k)
                        tmp += ctx->Hess[(k - 1) * (k - 2) / 2 + l] *
                                         ctx->ConQ[(ctx->CIP[i] - 1) * ctx->NParm1 + l] *
                                         ctx->ConQ[(ctx->CIP[j] - 1) * ctx->NParm1 + k];
                    else
                        tmp += ctx->Hess[(l - 1) * (l - 2) / 2 + k] *
                                         ctx->ConQ[(ctx->CIP[i] - 1) * ctx->NParm1 + l] *
                                         ctx->ConQ[(ctx->CIP[j] - 1) * ctx->NParm1 + k];
                }
            }
            if (i == j)
                ctx->WrkD[i] = tmp;
            else
                ctx->WrkH[(i - 1) * (i - 2) / 2 + j] = tmp;
        }
    }
    for (i = 1; i <= ctx->NParm; ++i)
        ctx->Diag[i] = ctx->WrkD[i];
    for (i = 1; i <= ctx->HSiz; ++i)
        ctx->Hess[i] = ctx->WrkH[i];
}

/* ------------------------------------------------------------------------ */
/*  chinv(opt)                                                              */
/*      If opt = 1, this function solves Hess * Srch = -Grad. First, a      */
/*      Cholesky decomposition of Hess is performed and the result is       */
/*      stored again in Hess (and Diag). If Hess is not positive definite,  */
/*      the function return with an positive index. Else the search vector  */
/*      is calculated by solving Hess * Srch = -Grad, and the function      */
/*      returns zero. If opt = 0, it is only checked if the hessian matrix  */
/*      is positive definite.                                               */
/*                                                                          */
/*      Return 0 if positive definite, else 1.                              */

int chinv(TDAContext *ctx, int opt)
{
    register int i,j,k,l,m;
    double tmp,tmp1;

    for (i = 1; i <= ctx->NParm1; ++i) {
        m = (i - 1) * (i - 2) / 2;
        tmp1 = ctx->Diag[i];
        for (j = 1; j < i; ++j) {
            l = (j - 1) * (j - 2) / 2;
            tmp = ctx->Hess[m + j];
            for (k = 1; k < j; ++k)
                tmp -= ctx->Hess[m + k] * ctx->Hess[++l];
            tmp /= ctx->Diag[j];
            ctx->Hess[m + j] = tmp;
            tmp1 -= tmp * tmp;
        }
        if (tmp1 < ctx->EPSI1)
            return(i);
        ctx->Diag[i] = sqrt(tmp1);
    }
    if (!opt)
        return(0);

    /* now solution of Hess * Srch = -Grad */

    for (i = 1; i <= ctx->NParm1; ++i) {      /* forward substitution */
        k = (i - 1) * (i - 2) / 2;
        tmp = -ctx->Grad[i];
        for (j = 1; j < i; ++j)
            tmp -= ctx->Hess[++k] * ctx->WrkD[j];
        ctx->WrkD[i] = tmp / ctx->Diag[i];
    }
    for (i = ctx->NParm1; i >= 1; --i) {      /* backwards substitution */
        tmp = ctx->WrkD[i];
        for (j = i + 1; j <= ctx->NParm1; ++j) {
            k = (j - 1) * (j - 2) / 2;
            tmp -= ctx->Hess[k + i] * ctx->Srch[j];
        }
        ctx->Srch[i] = tmp / ctx->Diag[i];
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  ffmin7  Function minimization with CES (quadratic or tensor model).     */
/*                                                                          */
/*  Return 0 if OK, -1 if insufficient memory, 1 if error in function       */

int ffmin7(TDAContext *ctx, int mod)
{
    register int i,j,k;
    int err,r,wa,na,*acn;
    double *typx,*acw,*acy,*acu,*acv,*w;

    err = -1;
    na = wa = 0;

    if (!(w = (double *)calloc((size_t)(5) * (size_t)(ctx->NParm1) + 15,sizeof(double))))
        goto FFMIN7Fin;

    wa = 5 * ctx->NParm1 + 15;
    memrq(ctx, wa,sizeof(double));

    typx = w;
    acw  = typx + ctx->NParm1 + 1;
    acy  = acw  + ctx->NParm1 + 1;
    acu  = acy  + ctx->NParm1 + 1;
    acv  = acu  + ctx->NParm1 + 1;

    if (!(acn = (int *)calloc((size_t)(ctx->NParm1 + 1),sizeof(int))))
        goto FFMIN7Fin;

    na = ctx->NParm1 + 1;
    memrq(ctx, na,sizeof(int));

    for (i = 1; i <= ctx->NParm1; ++i)
        typx[i] = 1.0;

    ctx->TSFNC = 0;      /* use standard models */
    ctx->TSMOD = mod;
    if (ctx->MINA == 7)
        ctx->TSMETH = 0;
    else
        ctx->TSMETH = 1;

    r = ts_min(ctx, ctx->NParm1,ctx->Par,typx,ctx->Par1,ctx->GTmp,ctx->TSGrad,ctx->Srch,acw,acy,ctx->WrkE,acu,acv,
               ctx->TSHess,acn);

    if (r < 0) {            /* error in function evaluation */
        err = 1;
        goto FFMIN7Fin;
    }
    switch (r) {
        case 1:
        case 2:     ctx->LConv = 0;
                    break;
        case 3:
        case 5:     ctx->LConv = -2;
                    break;
        case 4:     ctx->LConv = 2;
                    break;
        default:    ctx->LConv = 1;
                    break;
    }

    /* if we only have function values then put gradient into Grad[]
       and Hessian into Hess[] and Diag[] */

    if (ctx->PMDOPT == 0) {
        k = 1;
        for (i = 1; i <= ctx->NParm1; ++i) {
            for (j = 1; j < i; ++j)
                ctx->Hess[k++] = ctx->TSHess[(i - 1) * ctx->NParm1 + j] / ctx->DScal;
            ctx->Diag[i] = ctx->TSHess[(i - 1) * ctx->NParm1 + i] / ctx->DScal;
            ctx->Grad[i] = ctx->TSGrad[i] / ctx->DScal;
        }
    }
    err = 0;

FFMIN7Fin:

    if (wa > 0) {
        free((char *)w);
        memrq(ctx, -wa,sizeof(double));
    }
    if (na > 0) {
        free((char *)acn);
        memrq(ctx, -na,sizeof(int));
    }
    return(err);
}
