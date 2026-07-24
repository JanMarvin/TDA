/****************************************************************************/
/*  t_dem                                                                   */
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
#include "t_gf.h"
#include "t_ml.h"
#include "t_alloc.h"
#include "t_mat.h"
#include "t_matc.h"
#include "t_svd.h"
#include "tda_context.h"

/*  functions in t_dem.c */

int spmod(TDAContext *ctx);
double enorm(TDAContext *ctx, int n,double *x);
double r_alph(TDAContext *ctx, int n,double *x,double *y);
int etest(TDAContext *ctx);
int rap(TDAContext *ctx);
double rap_rate(TDAContext *ctx, int i,int j,int n,int na,int ny);

/* ------------------------------------------------------------------------ */
/*  spmod()     Stationary solution of Leslie-matrix                        */
/*                                                                          */
/*              spmod(                                                      */
/*                  tolf=...,       tolerance, def. 1.e-6                   */
/*                  mxit=...,       max number of iterations, def. 100      */
/*                  df=...,         output file                             */
/*                  prot=...,       protocol file, each iteration           */
/*                  fmt=...,        print format, def. 8.6                  */
/*              ) = F,S,X;          variables: fertility, survivor rates    */  
/*                                  and starting values                     */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int spmod(TDAContext *ctx)
{
    register int i,j;
    int err,iter,i0,i1,i2,conv;
    double alpha,a;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Stationary solution of Leslie-matrix. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLF = 1.0e-6;

    if (parm(ctx, ctx->CmdBuf + 5,4,1))     /* get parameters */
        goto DEMFin;
   
    if (ctx->PMNV != 3) {
        printf1(ctx, "Error: need three variables.\n");
        goto DEMFin;
    }
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 8,6);

    i0 = ctx->PMVIdx[0];
    i1 = ctx->PMVIdx[1];
    i2 = ctx->PMVIdx[2];

    /* allocate memory */

    if (alloc_acx(ctx, ctx->NOC + 1))
        goto DEMFin;
    if (alloc_acy(ctx, ctx->NOC + 1))
        goto DEMFin;

    if (ctx->MxItFlg == 0)
        ctx->MxIter = 100;
    printf1(ctx, "Max number of iterations: %d\n",ctx->MxIter);
    printf1(ctx, "Tolerance for convergence: %lg\n\n",ctx->TOLF);

    for (i = 0; i < ctx->NOC; ++i)
        ctx->AcX[i] = get_data(ctx, i2,i);
    enorm(ctx, ctx->NOC,ctx->AcX);

    iter = 0;
    if (ctx->PMProtFDef) {
        fprintf(ctx->PMProtFd,"%6d ",iter);
        for (i = 0; i < ctx->NOC; ++i)
            rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMFmtS,ctx->AcX[i]);
        fprintf(ctx->PMProtFd,"\n");
    }
    conv = 0;
    for (iter = 1; iter <= ctx->MxIter; ++iter) {

        ctx->AcY[0] = 0.0;
        for (i = 0; i < ctx->NOC; ++i)
            ctx->AcY[0] += ctx->AcX[i] * get_data(ctx, i0,i);

        for (i = 1; i < ctx->NOC; ++i)
            ctx->AcY[i] = ctx->AcX[i - 1] * get_data(ctx, i1,i - 1);

        alpha = r_alph(ctx, ctx->NOC,ctx->AcX,ctx->AcY);
        enorm(ctx, ctx->NOC,ctx->AcY);

        if (ctx->PMProtFDef) {
            fprintf(ctx->PMProtFd,"%6d ",iter);
            for (i = 0; i < ctx->NOC; ++i)
                rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMFmtS,ctx->AcY[i]);
            rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMFmtS,alpha);
            fprintf(ctx->PMProtFd,"\n");
        }

        /*  a holds the PREVIOUS alpha and is assigned at the foot of
            this loop, but the loop starts at iter = 1, so "iter > 0"
            was true on the first pass and compared alpha against a
            value never assigned (Valgrind, spmod.cf).  There is no
            previous alpha to compare on the first iteration.  */
        if (iter > 1) {
            if (fabs(alpha - a) <= ctx->TOLF) {
                conv = 1;
                break;
            }
        }
        for (i = 0; i < ctx->NOC; ++i)  
            ctx->AcX[i] = ctx->AcY[i];
        a = alpha;
    }
    if (conv == 0) {        /* if no convergence */
        printf1(ctx, "No convergence.\n");
        goto DEMFin;
    }
    a = 0.0;
    for (i = 0; i < ctx->NOC; ++i)    
        a = dmax(ctx, a,fabs(ctx->AcX[i] - ctx->AcY[i]));

    printf1(ctx, "Convergence after %d iterations.\n",iter);
    printf1(ctx, "Final growth factor: ");
    rt_printf1_d(ctx, ctx->PMFmtS,alpha);
    printf1(ctx, "\nMax difference in stationary vector: %lg\n",a);

#ifdef TDA_R_PACKAGE
    /* direct exports (CONTRIBUTING.md): the growth factor and
       the normalized stationary vector */
    tda_export_cell(ctx, "spmod.growth", alpha);
    tda_export_endrow(ctx, "spmod.growth");
    /* the converged vector lives in AcY (the iteration's ping-pong
       buffer, already normalized); i2 still holds the START values */
    for (i = 0; i < ctx->NOC; ++i) {
        tda_export_cell(ctx, "spmod.stationary", ctx->AcY[i]);
        tda_export_endrow(ctx, "spmod.stationary");
    }
    tda_export_flush(ctx, "spmod.growth");
    tda_export_flush(ctx, "spmod.stationary");
#endif

    if (ctx->PMF1Def) {
        for (i = 0; i < ctx->NOC; ++i)
            ctx->AcX[i] = get_data(ctx, i2,i);
        enorm(ctx, ctx->NOC,ctx->AcX);

        for (i = 0; i < ctx->NOC; ++i) {
            fprintf(ctx->PMF1d,"%4d ",i + 1);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[i]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[i]);
            fprintf(ctx->PMF1d,"\n");
        }
        printf1(ctx, "Results written to: %s\n",ctx->PMF1dName);
    }   
    if (ctx->PMProtFd)
        printf1(ctx, "Protocol written to: %s\n",ctx->PMProtFName);

    a = 0.0;
    for (i = 0; i < ctx->NOC; ++i) {
        a += ctx->AcY[i];
    }
    tda_out("Sum of Y = %g\n",a);


    for (i = 0; i < ctx->NOC; ++i) {
        a = 0.0;
        if (i == 0) {
            for (j = 0; j < ctx->NOC; ++j) 
                a += ctx->AcY[j] * get_data(ctx, i0,j);
        }
        else {
            a = get_data(ctx, i1,i - 1) * ctx->AcY[i - 1];
        }
        ctx->AcX[i] = a;
    }
    a = 0.0;
    for (i = 0; i < ctx->NOC; ++i) {
        a += ctx->AcX[i];
    }
    tda_out("Sum of next Y = %g\n",a);

    for (i = 0; i < ctx->NOC; ++i) {
        tda_out("%4d %16.8lf %16.8lf %16.8lf\n",i+1,ctx->AcY[i],ctx->AcX[i],ctx->AcY[i] * alpha);
    }
                




    err = 0;

DEMFin:
    p_clean(ctx);  
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  enorm(n,x)  return sum of components of x, scale x to                   */
/*              sum of components = 1.                                      */

double enorm(TDAContext *ctx, int n,double *x)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;
    double xs;

    xs = 0.0;
    for (i = 0; i < n; ++i) {
        xs  += x[i];
    }
    if (xs > 0.0) {
        for (i = 0; i < n; ++i)
            x[i] /= xs;
    }
    return(xs);
}

/* ------------------------------------------------------------------------ */
/*  r_alph(n,x,y)   return alpha.                                           */
/*             sum of components = 1.                                      */

double r_alph(TDAContext *ctx, int n,double *x,double *y)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;
    double xs,ys;

    xs = ys = 0.0;
    for (i = 0; i < n; ++i) {
        xs += x[i];
        ys += y[i];
    }
    if (xs > 0.0)
        return(ys / xs);
    return(0.0);
}

/* ------------------------------------------------------------------------ */
/*  etest()     Test of eigenvalue/vector calculations.                     */
/*  ##                                                                      */
/*              etest(                                                      */
/*                  alg=...,        algorithm, def. 1                       */
/*                                  1 =                                     */
/*                                  2 = eigen1                              */
/*                                  3 = eigen2                              */
/*                  mxit=...,       max iterations, alg. 3, def. 100        */  
/*                  fmt=...,        print format, def. -19,11               */
/*              ) = matrix;                                                 */  
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int etest(TDAContext *ctx)
{
    register int i,j,k;
    int err,r,row,col;
    /*  ivflg is an INPUT to n_mexpr -> get_mexpr -> v_parse, which uses
        it to decide whether interval operators are allowed.  It was
        declared and passed uninitialised (Valgrind, etest.cf), so
        whether intervals were accepted depended on stack contents.
        0 = not allowed, which is what a plain etest wants.  */
    int ivflg = 0;
    double dr,di,tmp,tmp1;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Test of eigenvalue/vector calculation. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,8,1))     /* get parameters */
        goto ETFin;
   
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, -19,11);

    if (n_mexpr(ctx, ctx->PMRHSTR,1,0,&row,&col,&ivflg,NULL) == NULL)
        goto ETFin;
   
    if (row != col) {
        mat_err(ctx, 1);
        goto ETFin;
    }
    if (ctx->PMALG < 1 || ctx->PMALG > 3)
        ctx->PMALG = 1;

    printf1(ctx, "Algorithm: %d",ctx->PMALG);
    if (ctx->MxItFlg == 0)
        ctx->MxIter = 100;
    if (ctx->PMALG == 3)
        printf1(ctx, " [max iterations: %d]",ctx->MxIter);
    newline(ctx);

    if (alloc_acw(ctx, row * col + 1))
        goto ETFin;

    for (i = 1; i <= row * col; ++i)
        ctx->AcW[i] = ctx->MX[0][i];

    printf1(ctx, "\n%s\n",ctx->PMRHSTR);
    for (i = 1; i <= row; ++i) {
        for (j = 1; j <= col; ++j)  
            rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcW[(i - 1) * col + j]);
        newline(ctx);
    }

    if (alloc_acn(ctx, col + 1))
        goto ETFin;
    if (alloc_acx(ctx, col + 1))
        goto ETFin;
    if (alloc_acy(ctx, col + 1))
        goto ETFin;
    if (alloc_acu(ctx, row * col + 1))
        goto ETFin;
    if (alloc_acv(ctx, row * col + 1))
        goto ETFin;

    if (ctx->PMALG == 1)  
        r = eigen(ctx, col,ctx->MX[0],ctx->AcX,ctx->AcY,ctx->AcU,1);
    else if (ctx->PMALG == 2)  
        r = eigen1(ctx, col,ctx->MX[0],ctx->AcX,ctx->AcY,ctx->AcU,ctx->AcV,ctx->AcN);
    else {
        r = eigen2(ctx, col,ctx->MX[0],ctx->AcX,ctx->AcY,ctx->AcU,1,ctx->MxIter);
        if (r > 0 && r < ctx->MxIter) {
            printf1(ctx, "Required %d iterations.\n",r);
            r = 0;
        }
    }       
     
    if (r) {
        printf1(ctx, "No success in calculations (err=%d).\n",r);
        goto ETFin;
    }
     
    if (ctx->PMALG != 2) {
        for (k = 1; k <= col; ++k) {
            if (fabs(ctx->AcY[k]) != 0.0) {
                for (i = 1; i <= row; ++i) {
                    ctx->AcV[(i - 1) * col + k] = ctx->AcU[(i - 1) * col + k + 1];
                    ctx->AcV[(i - 1) * col + k + 1] = -ctx->AcU[(i - 1) * col + k + 1];
                    ctx->AcU[(i - 1) * col + k + 1] = ctx->AcU[(i - 1) * col + k];
                }
                k++;
            }
        }
    }
    for (k = 1; k <= col; ++k) {
   
        printf1(ctx, "\n(%3d) ",k);
        rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcX[k]);
        rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcY[k]);
        printf1(ctx, "  [eigenvalue]\n\n");
#ifdef TDA_R_PACKAGE
        /*  The printed form is rounded by PMFmtS, so an R caller that
            parsed it back would compare a fit against the format rather
            than against the numbers.  Exported at full precision: one
            row per eigenvalue, real part then imaginary.  */
        tda_export_cell(ctx, "etest.eigen", (double)k);
        tda_export_cell(ctx, "etest.eigen", ctx->AcX[k]);
        tda_export_cell(ctx, "etest.eigen", ctx->AcY[k]);
        tda_export_endrow(ctx, "etest.eigen");
#endif

        for (i = 1; i <= row; ++i) {
            printf1(ctx, "      ");
            rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcU[(i - 1) * col + k]);
            rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcV[(i - 1) * col + k]);
            if (i == 1)
                printf1(ctx, "  [eigenvector]");
            newline(ctx);
        }

        /* check */

        dr = di = 0.0;
        for (i = 1; i <= row; ++i) {
            tmp = 0.0;
            for (j = 1; j <= col; ++j)  
                tmp += ctx->AcW[(i - 1) * col + j] * ctx->AcU[(j - 1) * col + k];
            tmp1 = ctx->AcX[k] * ctx->AcU[(i - 1) * col + k] -
                               ctx->AcY[k] * ctx->AcV[(i - 1) * col + k];
            dr = dmax(ctx, dr,fabs(tmp - tmp1));
        }
        for (i = 1; i <= row; ++i) {
            tmp = 0.0;
            for (j = 1; j <= col; ++j)  
                tmp += ctx->AcW[(i - 1) * col + j] * ctx->AcV[(j - 1) * col + k];
            tmp1 = ctx->AcX[k] * ctx->AcV[(i - 1) * col + k] +
                               ctx->AcY[k] * ctx->AcU[(i - 1) * col + k];
            di = dmax(ctx, di,fabs(tmp - tmp1));
        }
        printf1(ctx, "\n      ");
        rt_printf1_d(ctx, ctx->PMFmtS,dr);
        rt_printf1_d(ctx, ctx->PMFmtS,di);
        printf1(ctx, "  [check]\n");
    }
    err = 0;

ETFin:
    mx_free(ctx);
    p_clean(ctx);  
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  rap     Rates in age-period form.                                       */
/*                                                                          */
/*          rap(                                                            */
/*              year=y1,y2,     range of years for table construction       */
/*              age=a1,a2,      range of ages for table construction        */
/*              nc=...,         smoothing, def. nc=0                        */
/*              df=...,         output file for tables                      */
/*              nfmt=...,       integer print format, def. 4                */
/*                                                                          */
/*          ) = TS,TC,TF,D;                                                 */
/*                                                                          */
/*  Values of variables will be truncated to integers.                      */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int rap(TDAContext *ctx)
{
    register int i,j;
    int err,its,itc,itf,id,mints,maxts,d,t,ts,tc,tf,nn,ny,na,a,nr,ne;
    double r;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Rates in age-period form. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 3,4,1))     /* get parameters */
        goto RAPFin;
   
    if (ctx->PMNV != 4) {
        printf1(ctx, "Error: need four variables on right-hand side.\n");
        goto RAPFin;
    }
    if (ctx->PMNC < 0)
        ctx->PMNC = 0;   

    its = ctx->PMVIdx[0];
    itc = ctx->PMVIdx[1];
    itf = ctx->PMVIdx[2];
    id  = ctx->PMVIdx[3];

    for (i = 0; i < ctx->NOC; ++i) {
        ts = (int)(get_data(ctx, its,i));
        tc = (int)(get_data(ctx, itc,i));
        tf = (int)(get_data(ctx, itf,i));
        if (ts > tc || tc > tf) {
            printf1(ctx, "Inconsistent data in case %d\n",i + 1);
            goto RAPFin;
        }
    }
    mints = maxts = (int)get_data(ctx, its,0);

    for (i = 1; i < ctx->NOC; ++i) {
        t = (int)(get_data(ctx, its,i));
        mints = imin(ctx, mints,t);
        maxts = imax(ctx, maxts,t);
    }

    if (mints < 0) {
        printf1(ctx, "Error: minimal starting time is negative: %d\n",mints);
        goto RAPFin;
    }
    nn = maxts - mints + 1;

    if (alloc_acn(ctx, nn + 1))
        goto RAPFin;
    if (alloc_aci(ctx, nn + 1))
        goto RAPFin;
    if (alloc_acr(ctx, nn + 1))
        goto RAPFin;
    if (alloc_acs(ctx, nn + 1))
        goto RAPFin;

    for (i = 0; i < nn; ++i)
        ctx->AcR[i] = ctx->AcS[i] = -1;

    for (i = 0; i < ctx->NOC; ++i) {
        ts = (int)(get_data(ctx, its,i));
        tc = (int)(get_data(ctx, itc,i));
        tf = (int)(get_data(ctx, itf,i));
        d = (int)(get_data(ctx, id ,i));
        j = ts - mints;
        ctx->AcN[j] += 1;
        if (d != 0)
            ctx->AcI[j] += 1;

        if (ctx->AcR[j] < 0)
            ctx->AcR[j] = tc;
        else
            ctx->AcR[j] = imin(ctx, ctx->AcR[j],tc);

        if (ctx->AcS[j] < 0)
            ctx->AcS[j] = tf;
        else
            ctx->AcS[j] = imax(ctx, ctx->AcS[j],tf);
    }

    printf1(ctx, "\n  Idx       TS        N      D=1      D=0   Min TC   Max TF\n");
    prnchar(ctx, '-',59,1);

    j = 1;
    for (i = 0; i < nn; ++i) {
        if (ctx->AcN[i] < 1)
            continue;
        printf1(ctx, "%5d %8d %8d %8d %8d %8d %8d\n",j++,mints + i,ctx->AcN[i],
            ctx->AcI[i],ctx->AcN[i] - ctx->AcI[i],ctx->AcR[i],ctx->AcS[i]);
    }
    newline(ctx);

    if (ctx->PMYEAR1 < 0 || ctx->PMAGE1 < 0) {
        err = 0;
        goto RAPFin;
    }
    printf1(ctx, "Table construction for years %d-%d and ages %d-%d.\n\n",
        ctx->PMYEAR1,ctx->PMYEAR2,ctx->PMAGE1,ctx->PMAGE2);

    ny = ctx->PMYEAR2 - ctx->PMYEAR1 + 1;
    na = ctx->PMAGE2 - ctx->PMAGE1 + 1;

    if (alloc_acn(ctx, ny * na + 1))                      
        goto RAPFin;
    if (alloc_acm(ctx, ny * na + 1))
        goto RAPFin;

    nr = ne = 0;

    for (i = 0; i < ctx->NOC; ++i) {
        ts = (int)(get_data(ctx, its,i));
        tc = (int)(get_data(ctx, itc,i));
        tf = (int)(get_data(ctx, itf,i));
        d = (int)(get_data(ctx, id ,i));

        for (j = tc; j <= tf; ++j) {
            if (j < ctx->PMYEAR1 || j > ctx->PMYEAR2)
                continue;
            a = j - ts;
            if (a < ctx->PMAGE1 || a > ctx->PMAGE2)
                continue;

            t = j - ctx->PMYEAR1;
            a -= ctx->PMAGE1;

            if (t < 0 || a < 0 || t >= ny || a >= na) {
                tda_out("t=%d a=%d\n",t,a);
                exit(0);
            }
            ctx->AcN[a * ny + t] += 1;
            nr++;
            if (d != 0 && j == tf) {
                ctx->AcM[a * ny + t] += 1;
                ne++;
            }
        }
    }

    printf1(ctx, "                Entries       Cells       Empty\n");
    prnchar(ctx, '-',47,1);
    t = d = 0;
    a = ny * na;
    for (i = 0; i < a; ++i) {
        if (ctx->AcN[i] == 0)
            t++;
        if (ctx->AcM[i] == 0)
            d++;
    }
    printf1(ctx, "Risk table   %10d  %10d  %10d\n",nr,a,t);
    printf1(ctx, "Event table  %10d  %10d  %10d\n",ne,a,d);
    newline(ctx);


    if (ctx->PMF1Def) {
        fprintf(ctx->PMF1d,"\nRisk table\n");
        fprnchar(ctx, ctx->PMF1d,' ',ctx->PMNFmt + 1,0);
        for (j = 0; j < ny; ++j)
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->PMYEAR1 + j);
        fprintf(ctx->PMF1d,"\n");
        for (i = na - 1; i >= 0; --i) {
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->PMAGE1 + i);
            for (j = 0; j < ny; ++j)
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->AcN[i * ny + j]);
            fprintf(ctx->PMF1d,"\n");
        }

        fprintf(ctx->PMF1d,"\nEvent table\n");
        fprnchar(ctx, ctx->PMF1d,' ',ctx->PMNFmt + 1,0);
        for (j = 0; j < ny; ++j)
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->PMYEAR1 + j);
        fprintf(ctx->PMF1d,"\n");
        for (i = na - 1; i >= 0; --i) {
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->PMAGE1 + i);
            for (j = 0; j < ny; ++j)
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->AcM[i * ny + j]);
            fprintf(ctx->PMF1d,"\n");
        }
        fprintf(ctx->PMF1d,"\nRates\n");
        fprnchar(ctx, ctx->PMF1d,' ',ctx->PMNFmt + 1,0);
        for (j = 0; j < ny; ++j)
            fprintf(ctx->PMF1d," %6d",ctx->PMYEAR1 + j);
        fprintf(ctx->PMF1d,"\n");
        for (i = na - 1; i >= 0; --i) {
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->PMAGE1 + i);
            for (j = 0; j < ny; ++j) {
                r = rap_rate(ctx, i,j,ctx->PMNC,na,ny);
                fprintf(ctx->PMF1d," %6.2f",r);
            }
            fprintf(ctx->PMF1d,"\n");
        }
        printf1(ctx, "Tables written to: %s\n",ctx->PMF1dName);
    }

    err = 0;

RAPFin:
    p_clean(ctx);  
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  rap_rate(i,j,n)                                                         */

double rap_rate(TDAContext *ctx, int i,int j,int n,int na,int ny)
{
    int ii,iii,jj,jjj,nr,ne;
    double r;

    nr = ne = 0;
    for (ii = i - n; ii <= i + n; ++ii) {
        iii = ii;
        if (iii < 0)
            iii = 0;
        else if (iii >= na)
            iii = na - 1;

        for (jj = j - n; jj <= j + n; ++jj) {
                
            jjj = jj;
            if (jjj < 0)
                jjj = 0;
            else if (jjj >= ny)
                jjj = ny - 1;

            nr += ctx->AcN[iii * ny + jjj];
            ne += ctx->AcM[iii * ny + jjj];
        }
    }
    if (nr == 0)
        r = -1.0;
    else 
        r = 100.0 * (double)ne / (double)nr;                  
    return(r);
}


