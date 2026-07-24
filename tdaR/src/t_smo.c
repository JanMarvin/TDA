/****************************************************************************/
/*  t_smo                                                                   */
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
#include "t_var.h"
#include "t_alloc.h"
#include "t_gdat.h"
#include "t_spl.h"
#include "t_sort.h"
#include "t_ml.h"
#include "t_gf.h"
#include "t_areg.h"
#include "t_plot.h"
#include "tda_context.h"

/*  functions in t_smo.c */

int sma(TDAContext *ctx);
int s_sma(TDAContext *ctx, int n,double *x,int nw,double *w,int r,int opt);
int smd(TDAContext *ctx);
int s_smd(TDAContext *ctx, char *sm,int opt,int n,double *y);
int s_smdd(TDAContext *ctx, int n,double *x,double *y,int ns,int opt,int r,int dir);
double med3(TDAContext *ctx, double *x);
double med4(TDAContext *ctx, double *x);
double med5(TDAContext *ctx, double *x);
void nnsmooth(TDAContext *ctx, int n,double *x,double *y,int r);
void nneighbor(TDAContext *ctx, int n,double *x,int r,int k,int *k1,int *k2);
int spl(TDAContext *ctx);
int scplot(TDAContext *ctx);
void scplot1(TDAContext *ctx, int n,double *x,double *y,int ns);
int scplot2(TDAContext *ctx, int n);
void scplot2a(TDAContext *ctx, int n,double x,double y);

/* ------------------------------------------------------------------------ */
/*  sma()           Moving averages.                                        */
/*                                                                          */
/*                  sma(                                                    */
/*                      gss=...,        weights, required                   */
/*                      opt=...,        rules for ends                      */
/*                      r=...,          repeat factor, def. 1               */
/*                      fmt=...,        print format, def. 10.4             */
/*                      df=...,         optional output file                */
/*                      sel=...,        optional case selection             */
/*                  ) = varlist;                                            */
/*                                                                          */
/*                  Return 0 if OK, otherwise -1.                           */

int sma(TDAContext *ctx)
{
    register int i;
    int err,n,iv,ivv,ix,nrec;
    double w;  

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Moving averages. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 3,4,1))     /* get parameters */
        goto SMAFin;

    if (ctx->PMR < 1)
        ctx->PMR = 1;
    if (ctx->PMOPT != 2)
        ctx->PMOPT = 1;
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    if (ctx->PMNTP < 1) {
        printf1(ctx, "Error: gss parameter must specify at least two weights.\n");
        goto SMAFin;
    }

    printf1(ctx, "Weights:");
    for (i = ctx->PMNTP - 1; i >= 0; --i)
        printf1(ctx, " %g",ctx->PMTP[i]);
    for (i = 1; i < ctx->PMNTP; ++i)
        printf1(ctx, " %g",ctx->PMTP[i]);
    newline(ctx);

    n = ctx->PMNTP - 1;
    w = ctx->PMTP[0];
    for (i = 1; i <= n; ++i) {
        if (ctx->PMTP[i] < 0.0) {
            printf1(ctx, "Error: negative weights are not allowed.\n");
            goto SMAFin;
        }
        w += 2.0 * ctx->PMTP[i];
    }
    printf1(ctx, "Sum of weights: %g\n",w);

    printf1(ctx, "Using: ");
    if (ctx->PMOPT == 1)
        printf1(ctx, "copy-on end-value rule.\n");
    else            
        printf1(ctx, "replicate end-value rule.\n");
    printf1(ctx, "Repeat factor: %d\n",ctx->PMR);

    if (alloc_acx(ctx, ctx->NOC))         /* allocate AcX */
        goto SMAFin;

    if (alloc_acy(ctx, ctx->NOC))         /* allocate AcY */
        goto SMAFin;
     
    nrec = ix = 0;
    for (iv = 0; iv < ctx->PMNV; ++iv) {     /* for all variables */

        ivv = ctx->PMVIdx[iv];

        if (ctx->PMF1Def)  
            fprintf(ctx->PMF1d,"# Var: %s\n",ctx->VName[ivv]);
        else
            printf1(ctx, "\nVariable: %s\n",ctx->VName[ivv]);

        n = getd1(ctx, ctx->AcX,ivv,-1,1);
        if (n == 0)  
            continue;    
             
        for (i = 0; i < n; ++i)
            ctx->AcY[i] = ctx->AcX[i];

        s_sma(ctx, n,ctx->AcY,ctx->PMNTP,ctx->PMTP,ctx->PMR,ctx->PMOPT);

        for (i = 0; i < n; ++i) {
#ifdef TDA_R_PACKAGE
            /* the smoothed series, whichever of the two paths below
               prints it -- both write the same four values */
            {
                double erow[4];
                erow[0] = (double)ix;
                erow[1] = (double)(i + 1);
                erow[2] = ctx->AcX[i];
                erow[3] = ctx->AcY[i];
                tda_export_row(ctx, "smooth.table", erow, 4);
            }
#endif
            if (ctx->PMF1Def) {
                fprintf(ctx->PMF1d,"%3d %6d ",ix,i + 1);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[i]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[i]);
                fprintf(ctx->PMF1d,"\n");
                nrec++;
            }
            else {
                printf1(ctx, "%3d %6d ",ix,i + 1);
                rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcX[i]);
                rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcY[i]);
                printf1(ctx, "\n");
            }
        }
        ix++;
    }
    if (ctx->PMF1Def)  
        printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);

    err = 0;

SMAFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  s_sma(n,x,nw,w,r,opt)   moving averages                                 */
/*                                                                          */
/*  n = number of input data points in x[].                                 */
/*  at return, x[] contains the smoothed values.                            */
/*  nw = number of weights in w[]                                           */
/*  r = repeat factor                                                       */
/*  opt = 1 : copy-on end-value rule                                        */
/*  opt = 2 : replicate end-value rule                                      */
/*                                                                          */
/*  Return 0 if OK, -1 if insufficient memory.                              */

int s_sma(TDAContext *ctx, int n,double *x,int nw,double *w,int r,int opt)
{
    register int i,j,k;
    int ii,nn,na;
    double tmp;

    if (n < 1 || nw < 1)
        return(0);
        
    if (alloc_actmp(ctx, n))         /* allocate AcTmp */
        return(-1); 

    for (i = 0; i < n; ++i)
        ctx->AcTmp[i] = x[i];

    nw--;           
    nn = n - 1;
    na = 0;
    if (opt == 1) {
        nn -= nw;
        na = nw;
    }

    ii = 0;
    while (1) {
        for (i = na; i <= nn; ++i) {
            tmp = w[0] * x[i];     
            for (j = 1; j <= nw; ++j) {
                k = i + j;
                if (k >= n)
                    tmp += w[j] * x[n - 1];
                else
                    tmp += w[j] * x[k];
                k = i - j;
                if (k < 0)
                    tmp += w[j] * x[0];
                else
                    tmp += w[j] * x[k];
            }
            ctx->AcTmp[i] = tmp;           
        }
        for (i = 0; i < n; ++i)
            x[i] = ctx->AcTmp[i];

        if (++ii >= r)
            break;
    }
    alloc_actmp(ctx, 0);     /* free AcTmp */
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  smd(idx)        Running median smoother                                 */
/*                                                                          */
/*                  smd(                                                    */
/*                      sm=[...],       a sequence of: 2,3,4,5,H,R,S,t      */
/*                      opt=...,        rules for ends                      */
/*                      fmt=...,        print format, def. 10.4             */
/*                      df=...,         optional output file                */
/*                      sel=...,        optional case selection             */
/*                  ) = varlist;                                            */
/*                                                                          */
/*  sm=[....]   any sequence of: 2,3,4,5,H,R,S,t                            */
/*                                                                          */
/*  opt = 1 : copy-on end-value rule                                        */
/*  opt = 2 : replicate end-value rule                                      */
/*  opt = 3 : end-value smoothing                                           */  
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int smd(TDAContext *ctx)
{
    register int i;
    int err,n,iv,ivv,ix,nrec;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Running Median Smoother. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 3,4,1))     /* get parameters */
        goto SMDFin;

    if (ctx->PMRHSTRA < 1 || ctx->PMOPT < 1 || ctx->PMOPT > 3) {
        p_err(ctx, -1,1);
        goto SMDFin;
    }
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    printf1(ctx, "Using: sm=[%s]\nUsing: ",ctx->PMRHSTR);
    if (ctx->PMOPT == 1)
        printf1(ctx, "copy-on end-value rule.\n");
    else if (ctx->PMOPT == 2)
        printf1(ctx, "replicate end-value rule.\n");
    else
        printf1(ctx, "end-value smoothing.\n");

    if (alloc_acx(ctx, ctx->NOC))         /* allocate AcX */
        goto SMDFin;

    if (alloc_acy(ctx, ctx->NOC))         /* allocate AcY */
        goto SMDFin;

    if (ctx->PMF1Def)  
        fprintf(ctx->PMF1d,"# Running medians [%s]\n",ctx->PMRHSTR);

    nrec = ix = 0;
    for (iv = 0; iv < ctx->PMNV; ++iv) {     /* for all variables */

        ivv = ctx->PMVIdx[iv];

        if (ctx->PMF1Def)  
            fprintf(ctx->PMF1d,"# Var: %s\n",ctx->VName[ivv]);
        else
            printf1(ctx, "\nVariable: %s\n",ctx->VName[ivv]);

        n = getd1(ctx, ctx->AcX,ivv,-1,1);
        if (n == 0)  
            continue;    

        for (i = 0; i < n; ++i)
            ctx->AcY[i] = ctx->AcX[i];
           
        if (s_smd(ctx, ctx->PMRHSTR,ctx->PMOPT,n,ctx->AcY)) {
            printf1(ctx, "Syntax error: [%s]\n",ctx->PMRHSTR);
            goto SMDFin;
        }
        for (i = 0; i < n; ++i) {
#ifdef TDA_R_PACKAGE
            /* the smoothed series, whichever of the two paths below
               prints it -- both write the same four values */
            {
                double erow[4];
                erow[0] = (double)ix;
                erow[1] = (double)(i + 1);
                erow[2] = ctx->AcX[i];
                erow[3] = ctx->AcY[i];
                tda_export_row(ctx, "smooth.table", erow, 4);
            }
#endif
            if (ctx->PMF1Def) {
                fprintf(ctx->PMF1d,"%3d %6d ",ix,i + 1);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[i]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[i]);
                fprintf(ctx->PMF1d,"\n");
                nrec++;
            }
            else {
                printf1(ctx, "%3d %6d ",ix,i + 1);
                rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcX[i]);
                rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcY[i]);
                printf1(ctx, "\n");
            }
        }
        ix++;
    }
    if (ctx->PMF1Def)  
        printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);
    err = 0;

SMDFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  s_smd(sm=,opt,n,y)          running median smoother                     */
/*                                                                          */
/*  sm=[....]   any sequence of: 2,3,4,5,H,R,S,t                            */
/*                                                                          */
/*  opt = 1 : copy-on end-value rule                                        */
/*  opt = 2 : replicate end-value rule                                      */
/*  opt = 3 : end-value smoothing                                           */  
/*                                                                          */
/*  y[] contains n input values.                                            */
/*  On return, y[] contains the smoothed values.                            */
/*  Return 0 if OK, -1 if syntax error, -2 if insufficient memory.          */

int s_smd(TDAContext *ctx, char *sm,int opt,int n,double *y)
{
    register int i;
    int r,dir,tflag,err;
    register char *p;   
    double tmp;

    if (n <= 1)
        return(0);

    err = -2;
    if (alloc_acu(ctx, n))         /* allocate AcU */
        goto SSMDFin;

    if (alloc_acv(ctx, n))         /* allocate AcV */
        goto SSMDFin;

    err = -1;
    for (i = 0; i < n; ++i)  
        ctx->AcV[i] = y[i];

    tflag = 0;
SMDRep:      
    dir = 0;
    p = sm;         
    while (*p) {
        r = 0;
        if (*p == '2') {
            p++;
            s_smdd(ctx, n,y,ctx->AcU,2,opt,r,dir);
            if (dir)
                dir = 0;
            else
                dir = 1;
        }
        else if (*p == '3') {
            p++;
            if (*p == 'R') {
                r = 1;
                p++;
            }
            s_smdd(ctx, n,y,ctx->AcU,3,opt,r,0);
        }
        else if (*p == '4') {
            p++;
            s_smdd(ctx, n,y,ctx->AcU,4,opt,r,dir);
            if (dir)
                dir = 0;
            else
                dir = 1;
        }
        else if (*p == '5') {
            p++;
            if (*p == 'R') {
                r = 1;
                p++;
            }
            s_smdd(ctx, n,y,ctx->AcU,5,opt,r,0);
        }
        else if (*p == 'H') {
            p++;
            s_smdd(ctx, n,y,ctx->AcU,0,opt,0,0);
        }
        else if (*p == 'S') {
            p++;
            s_smdd(ctx, n,y,ctx->AcU,1,opt,0,0);
        }
        else if (*p == 't' && !*(p + 1)) {
            tflag++;
            break;
        }
        else  
            goto SSMDFin;
    }
    if (tflag == 1) {
        for (i = 0; i < n; ++i) {
            tmp = ctx->AcV[i];
            ctx->AcV[i] = y[i];
            y[i] = tmp - y[i];
        }
        goto SMDRep;
    }
    else if (tflag == 2) {
        for (i = 0; i < n; ++i)  
            y[i] += ctx->AcV[i];
    }
    err = 0;

SSMDFin:
    alloc_acu(ctx, 0);       /* free AcU */
    alloc_acv(ctx, 0);       /* free AcV */
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  s_smdd(n,x,y,ns,opt,r,dir)  running median calculation.                 */
/*                                                                          */
/*  n = number of input values in x[].                                      */
/*  y[] work array, output values also in x[].                              */
/*                                                                          */
/*  ns = 0  hanning                                                         */
/*       1  splitting                                                       */
/*       2  2-medians, depends on dir                                       */
/*       3  3-medians                                                       */
/*       4  4-medians, depends on dir                                       */
/*       5  5-medians                                                       */
/*                                                                          */
/*  if dir = 0 use previous, otherwise use next value.                      */
/*                                                                          */
/*  if r = 1 repeat until no change (only with ns = 3 or 5)                 */
/*                                                                          */
/*  opt = 1 : copy-on end-value rule                                        */
/*  opt = 2 : replicate end-value rule                                      */
/*  opt = 3 : end-value smoothing                                           */  

int s_smdd(TDAContext *ctx, int n,double *x,double *y,int ns,int opt,int r,int dir)
{
    register int i,j;
    int nn;
    double s[3];

    if (n <= 2)  
        return(0);

    nn = n - 1;
    while (1) {
        for (i = 0; i < n; ++i) {
            y[i] = x[i];
            if (ns == 0) {                  /* hanning */
                if (i > 0 && i < nn)
                    y[i] = 0.25 * x[i - 1] + 0.5 * x[i] + 0.25 * x[i + 1];
            }
            else if (ns == 1) {             /* splitting */

                if (i >= 2 && i <= nn - 3 && x[i] == x[i + 1]) {
                    s[2] = x[i];
                    s[1] = x[i - 1];
                    s[0] = 3.0 * x[i - 1] - 2.0 * x[i - 2];
                    y[i] = med3(ctx, s);

                    s[0] = x[i + 1];
                    s[1] = x[i + 2];
                    s[2] = 3.0 * x[i + 2] - 2.0 * x[i + 3];
                    y[i + 1] = med3(ctx, s);
                }
            }
            else if (ns == 2) {
                if (dir == 0) {
                    if (i > 0)
                        y[i] = 0.5 * x[i - 1] + 0.5 * x[i];
                }
                else {
                    if (i < nn)
                        y[i] = 0.5 * x[i] + 0.5 * x[i + 1];
                }
            }
            else if (ns == 3) {
                if (i > 0 && i < nn)   
                    y[i] = med3(ctx, x + i - 1);
            }
            else if (ns == 4) {
                if (i > 0 && i < nn) {
                    if (dir == 0) {
                        if (i == 1)
                            y[i] = med3(ctx, x);
                        else             
                            y[i] = med4(ctx, x + i - 2);
                    }
                    else {
                        if (i == nn - 1)
                            y[i] = med3(ctx, x + nn - 2);
                        else
                            y[i] = med4(ctx, x + i - 1);
                    }
                }
            }
            else if (ns == 5) {
                if (i > 0 && i < nn) {
                    if (i == 1)
                        y[i] = med3(ctx, x);
                    else if (i == nn - 1)
                        y[i] = med3(ctx, x + nn - 2);
                    else   
                        y[i] = med5(ctx, x + i - 2);
                }
            }
        }

        /* apply end-value option */

        if (opt == 2) {         /* nothing to do for 3 and 5 */

            if (ns == 0) {
                y[0] = 0.75 * x[0] + 0.25 * x[1];
                y[nn] = 0.75 * x[nn] + 0.25 * x[nn - 1];
            }
        }
        else if (opt == 3) {    /* E rule, nothing to do for H */

            if (ns == 3 || ns == 5) {
                s[0] = 3.0 * y[1] - 2.0 * y[2];
                s[1] = x[0];
                s[2] = y[1];
                y[0] = med3(ctx, s);
                s[0] = 3.0 * y[nn - 1] - 2.0 * y[nn - 2];
                s[1] = x[nn];
                s[2] = y[nn - 1];
                y[0] = med3(ctx, s);
            }
        }
        j = 0;
        for (i = 0; i < n; ++i) {
            if (fabs(x[i] - y[i]) > ctx->EPSI1)
                j = 1;
            x[i] = y[i];
        }
        if (r == 0 || j == 0)
            break;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  med3(x)     return median of x[0],x[1],x[2]                             */

double med3(TDAContext *ctx, double *x)
{
    (void)ctx;        /* unused: the signature is shared */
    double s[2];

    s[0] = x[0];
    s[1] = x[1];           
    if (s[1] < s[0]) {
        s[0] = x[1];
        s[1] = x[0];
    }
    if (x[2] <= s[0])
        return(s[0]);
    if (x[2] <= s[1])
        return(x[2]);
    return(s[1]);
}

/* ------------------------------------------------------------------------ */
/*  med5(x)     return median of x[0],x[1],x[2],x[3],x[4]                   */

double med5(TDAContext *ctx, double *x)
{
    register int i,k,jmin,jmax;
    double xmin,xmax,s[5];

    jmin = jmax = 0;
    xmin = xmax = x[0];
    for (i = 1; i < 5; ++i) {
        if (xmin > x[i]) {
            xmin = x[i];
            jmin = i;
        }
        if (xmax < x[i]) {
            xmax = x[i];
            jmax = i;
        }
    }
    k = 0;
    for (i = 0; i < 5; ++i) {
        if (i != jmin && i != jmax)
            s[k++] = x[i];
    }
    return(med3(ctx, s));
}

/* ------------------------------------------------------------------------ */
/*  med4(x)     return median of x[0],x[1],x[2],x[3]                        */

double med4(TDAContext *ctx, double *x)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,jmin,jmax;
    double xmin,xmax,tmp;

    jmin = jmax = 0;
    xmin = xmax = x[0];
    for (i = 1; i < 4; ++i) {
        if (xmin > x[i]) {
            xmin = x[i];
            jmin = i;
        }
        if (xmax < x[i]) {
            xmax = x[i];
            jmax = i;
        }
    }
    tmp = 0.0;
    for (i = 0; i < 4; ++i) {
        if (i != jmin && i != jmax)
            tmp += x[i];
    }
    return(tmp / 2.0);
}

/*--------------------------------------------------------------------------*/
/*  nnsmooth(n,x,y,r)                                                       */

void nnsmooth(TDAContext *ctx, int n,double *x,double *y,int r)
{
    register int i,j;
    int k1,k2;
    double tmp;

    for (i = 0; i < n; ++i) {
        y[i] = x[i];
        if (r < n) {
            nneighbor(ctx, n,x,r,i,&k1,&k2);
            tmp = 0.0;
            for (j = 0; j < r; ++j)
                tmp += x[k1 + j];
            if (r > 0)
                tmp /= (double)r;
            y[i] = tmp;
        }
    }
}

/*--------------------------------------------------------------------------*/
/*  nneighbor(n,x,r,k,k1,k2)                                                */
/*                                                                          */
/*  find k1 and k2 such that x[i] (i = k1,...,k2) are the r nearest         */
/*  neighbors of x[k].                                                      */

void nneighbor(TDAContext *ctx, int n,double *x,int r,int k,int *k1,int *k2)
{
    (void)ctx;        /* unused: the signature is shared */
    register int m,j1,j2;

    if (r > n)
        r = n;
    if (r == n) {
        *k1 = 0;
        *k2 = n - 1;
        return;
    }
    j1 = 0;
    j2 = 0;
    m = 1;
    while (1) {
        j1++;
        j2++;
        if (k - j1 < 0) {
            j1--;
            j2 = r - j1 - 1;
            break;
        }
        else if (k + j2 >= n) {
            j2--;
            j1 = r - j2 - 1;
        }
        else if (fabs(x[k] - x[k - j1]) <= fabs(x[k] - x[k + j2]))
            j2--;
        else
            j1--;
        if (++m >= r)
            break;
    }
    *k1 = k - j1;
    *k2 = k + j2;
}

/* ------------------------------------------------------------------------ */
/*  spl()           Smoothing spline function.                              */
/*                                                                          */
/*                  spl(                                                    */
/*                      sig=...,        smoothing factor, def. 0.0          */
/*                      deg=...,        degree, def. 2                      */
/*                      max=...,        max number of knots, def. NOC/2     */
/*                      fmt=...,        print format, def. 10.4             */
/*                      df=...,         optional output file                */
/*                      rx=a(d)b        range of interpolation              */
/*                  ) = X,Y [,W];                                           */
/*                                                                          */
/*                  Return 0 if OK, otherwise -1.                           */

int spl(TDAContext *ctx)
{
    register int i,j,l,k;
    int err,r,ix,iy,iw,nrec,nk;
    double xa,xb,x,y,eps;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Smoothing spline function. Current memory: %d bytes.\n",ctx->MemReq);
    ctx->MxIter = 0;

    if (parm(ctx, ctx->CmdBuf + 3,4,1))     /* get parameters */
        goto SPLFin;

    if (ctx->PMNV != 2 && ctx->PMNV != 3) {
        printf1(ctx, "Error: need two or three variables on right-hand side.\n");
        goto SPLFin;
    }
    if (ctx->PMSIG < 0.0)
        ctx->PMSIG = 0.0;
    if (ctx->PMDEG < 2)
        ctx->PMDEG = 2;
    if (ctx->PMMax < 3 * ctx->PMDEG + 1)
        ctx->PMMax = imax(ctx, ctx->NOC / 2 + 1,3 * ctx->PMDEG + 1);
    if (ctx->MxIter < 1)
        ctx->MxIter = 10;
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    printf1(ctx, "Degree: %d. Smoothing factor: %g\n",ctx->PMDEG,ctx->PMSIG);
    printf1(ctx, "Max number of iterations: %d\n",ctx->MxIter);
    printf1(ctx, "Max number of knots: %d\n",ctx->PMMax);

    if (alloc_acx(ctx, ctx->NOC + 1))         /* allocate AcX */
        goto SPLFin;

    if (alloc_acy(ctx, ctx->NOC + 1))         /* allocate AcY */
        goto SPLFin;
     
    if (alloc_acw(ctx, ctx->NOC + 1))         /* allocate AcW */
        goto SPLFin;
     
    if (alloc_acu(ctx, ctx->PMMax + 1))       /* allocate AcU */
        goto SPLFin;
     
    if (alloc_acv(ctx, ctx->PMMax + 1))       /* allocate AcV */
        goto SPLFin;
     
    if (alloc_actmp(ctx, ctx->PMDEG + 2))     /* allocate AcTmp */
        goto SPLFin;
     
    ix = ctx->PMVIdx[0];
    iy = ctx->PMVIdx[1];
    if (ctx->PMNV > 2)
        iw = ctx->PMVIdx[2];
    else
        iw = -1;

    for (i = 0; i < ctx->NOC; ++i) {
        ctx->AcX[i] = get_data(ctx, ix,i);
        ctx->AcY[i] = get_data(ctx, iy,i);
        if (iw >= 0) {
            ctx->AcW[i] = get_data(ctx, iw,i);
            if (ctx->AcW[i] <= 0.0) {
                printf1(ctx, "Error: found zero or negative weight in case %d.\n",i + 1);
                goto SPLFin;
            }
        }
        else
            ctx->AcW[i] = 1.0;
    }
    if (sortd3(ctx, ctx->NOC,ctx->AcX,ctx->AcY,ctx->AcW))        /* sort in ascending order */
        goto SPLFin;

    for (i = 1; i < ctx->NOC; ++i) {
        if (ctx->AcX[i] <= ctx->AcX[i - 1]) {
            printf1(ctx, "Error: cannot sort into strictly ascending order.\n");
            goto SPLFin;
        }   
    }
    nk = ctx->PMMax;
    xa = ctx->AcX[0];
    xb = ctx->AcX[ctx->NOC - 1];

    r = splf(ctx, ctx->NOC,ctx->AcX - 1,ctx->AcY - 1,ctx->AcW - 1,xa,xb,ctx->PMDEG,ctx->PMSIG,&nk,ctx->AcU,ctx->AcV,ctx->MxIter);

    if (r < 0) {
        printf1(ctx, "Error (%d)",r);
        if (r == -2) 
            printf1(ctx, ": exceeded max number of knots");
        else if (r == -3)
            printf1(ctx, ": exceeded max number of iterations");
        printf1(ctx, ".\n");
        goto SPLFin;
    }
    printf1(ctx, "Successful calculation of ");
    if (r == 0)
        printf1(ctx, "smoothing spline.\n");
    else if (r == 1)
        printf1(ctx, "interpolating spline.\n");
    else
        printf1(ctx, "OLS approximation.\n");
    printf1(ctx, "Actual number of knots: %d\n",nk);

    xa = ctx->AcU[ctx->PMDEG + 1];
    xb = ctx->AcU[nk - ctx->PMDEG];

    printf1(ctx, "Interval for interpolation: [%g,%g]\n",xa,xb);

    if (ctx->PMF1Def) {
        eps = 1.e-6;
        nrec = 0;
        if (ctx->PMRXFlg) {
            printf1(ctx, "Request is interpolation in [%g,%g], increment: %g\n",
                                            (double)(ctx->PMRXA),(double)(ctx->PMRXB),(double)(ctx->PMRXD));

            r = (int)((ctx->PMRXB - ctx->PMRXA) / ctx->PMRXD) + 1; 
            x = (double)(ctx->PMRXA);
            interv(ctx, nk,ctx->AcU,x,1);           /* initialization */
            k = 1;
            for (i = 0; i <= r; ++i) {
                x = (double)(ctx->PMRXA) + (double)(i) * (double)(ctx->PMRXD);
                if (x < xa)
                    continue;
                if (x > xb + eps)
                    break;

                fprintf(ctx->PMF1d,"%6d ",k++);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,x);
#ifdef TDA_R_PACKAGE
                tda_export_cell(ctx, "spl.table", (double)(k - 1));
                tda_export_cell(ctx, "spl.table", x);
#endif

                if (fabs(x - xb) <= eps)
                    l = nk - ctx->PMDEG - 1;
                else
                    l = interv(ctx, nk,ctx->AcU,x,0);

                for (j = 0; j <= ctx->PMDEG; ++j) {
                    y = spld(ctx, nk,ctx->AcU,ctx->PMDEG,ctx->AcV,j,x,l,ctx->AcTmp);
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,y);
#ifdef TDA_R_PACKAGE
                    tda_export_cell(ctx, "spl.table", y);
#endif
                }
                /* braces: unbraced, the export lands outside the if and
                   every row gains a weight cell */
                if (iw >= 0) {
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcW[i]);
#ifdef TDA_R_PACKAGE
                    tda_export_cell(ctx, "spl.table", ctx->AcW[i]);
#endif
                }

                fprintf(ctx->PMF1d,"\n");
#ifdef TDA_R_PACKAGE
                tda_export_endrow(ctx, "spl.table");
#endif
                nrec++;
            }
        }
        else {
            interv(ctx, nk,ctx->AcU,xa,1);           /* initialization */
            for (i = 0; i < ctx->NOC; ++i) {
                x = ctx->AcX[i];
                if (x < xa)
                    continue;
                if (x > xb)
                    break;

                fprintf(ctx->PMF1d,"%6d ",i + 1);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,x);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[i]);
#ifdef TDA_R_PACKAGE
                tda_export_cell(ctx, "spl.table", (double)(i + 1));
                tda_export_cell(ctx, "spl.table", x);
                tda_export_cell(ctx, "spl.table", ctx->AcY[i]);
#endif

                if (fabs(x - xb) <= eps)
                    l = nk - ctx->PMDEG - 1;
                else
                    l = interv(ctx, nk,ctx->AcU,x,0);

                for (j = 0; j <= ctx->PMDEG; ++j) {
                    y = spld(ctx, nk,ctx->AcU,ctx->PMDEG,ctx->AcV,j,x,l,ctx->AcTmp);
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,y);
#ifdef TDA_R_PACKAGE
                    tda_export_cell(ctx, "spl.table", y);
#endif
                }
                /* braces: unbraced, the export lands outside the if and
                   every row gains a weight cell */
                if (iw >= 0) {
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcW[i]);
#ifdef TDA_R_PACKAGE
                    tda_export_cell(ctx, "spl.table", ctx->AcW[i]);
#endif
                }

                fprintf(ctx->PMF1d,"\n");
#ifdef TDA_R_PACKAGE
                tda_export_endrow(ctx, "spl.table");
#endif
                nrec++;
            }

        }
        printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);
    }

    err = 0;

SPLFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  scplot()        Scatterplots                                            */
/*  ##                                                                      */
/*                  scplot(                                                 */
/*                      opt=...,        1   only symbols, default           */
/*                                      2   sunflower plot                  */
/*                                      3   lowess                          */
/*                                      4   midmeans                        */
/*                      lt=...,         line type, def. 1                   */
/*                      lw=...,         line width, def. 0.2                */
/*                      s=...,          symbol type, def. 1                 */
/*                      fs=...,         symbol size, def. 2 mm              */
/*                      nc=...,                                             */
/*                      nn=n1,n2,       for sunflower plot, def. 1,1        */
/*                      fss=...,        symbol size sunflower plot, def. 2  */
/*                      r=...,                                              */
/*                      sig=...,        for lowess, def. 0.5                */
/*                      ns=...,                                             */
/*                      d=...,          for lowess                          */
/*                      sel=...,        optional case selection             */
/*                      df=...,         optional output file                */
/*                      fmt=...,        print format, def. 10.4             */
/*                  ) = X,Y;            variables                           */
/*                                                                          */
/*                  Return 0 if OK, -1 if error.                            */

int scplot(TDAContext *ctx)
{
    register int i;
    int nn,err,ix,iy,first,nrec;
    double x,y;
           
    err = - 1;
    nrec = 0;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Scatterplot. Current memory: %d bytes.\n",ctx->MemReq);

    if (check_ps(ctx, 2))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 6,4,1))       /* get parameters */
        goto PLSCFin;

    if (ctx->PMOPT < 1 || ctx->PMOPT > 4)
        ctx->PMOPT = 1;

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    if (ctx->NOC < 2)
        goto PLSCFin;

    if (ctx->PMNV != 2) {
        p_err(ctx, -1,1);
        goto PLSCFin;
    }
    ix = ctx->PMVIdx[0];
    iy = ctx->PMVIdx[1];

    if (alloc_acx(ctx, ctx->NOC + 1))
        goto PLSCFin;
    if (alloc_acy(ctx, ctx->NOC + 1))
        goto PLSCFin;

    nn = getd2(ctx, ctx->AcX,ctx->AcY,ix,iy,1);
    if (nn < 2) {
        p_err(ctx, -26,1);
        goto PLSCFin;
    }
    if (ctx->PMOPT > 1) {
        if (sortd2(ctx, nn,ctx->AcX,ctx->AcY))     /* sort simultaneously according to AcX */
            goto PLSCFin;
    }
    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    fprintf(ctx->PSFd,"gsave\n");
    if (ctx->PMNC != 1)
        set_clip(ctx);
    ps_lwidth(ctx, ctx->PMLW);
    ps_ltyp(ctx, ctx->PMLT);

    if (ctx->PMOPT != 2) {   /* standard symbols */

        if (ctx->PMS >= 1 && ctx->PMS <= 17 && ctx->PMFS > 0.0) {  /* plot symbols */

            fprintf(ctx->PSFd,"gsave\n");
            fprintf(ctx->PSFd,"/ssiz %4.2f def\n",ctx->PtMM * ctx->PMFS / 2.0);

            for (i = 0; i < nn; ++i) {
                x = ps_2dx(ctx, ctx->AcX[i]);       
                y = ps_2dy(ctx, ctx->AcY[i]);       
                ps_sym(ctx, ctx->PMS,x,y,ctx->PtMM * ctx->PMFS / 2.0);
            }
            fprintf(ctx->PSFd,"grestore\n");
        }
        /***
        ps_lwidth(ctx, PMLW);
        ps_ltyp(ctx, PMLT);
        ***/
    }

    if (ctx->PMOPT == 2) {               /* sunflower plot */
        err = scplot2(ctx, nn);
    }
    else if (ctx->PMOPT == 3) {          /* lowess */

        if (alloc_acu(ctx, nn + 1))
            goto PLSCFin;
        if (alloc_acv(ctx, nn + 1))
            goto PLSCFin;
        if (alloc_actmp(ctx, nn + 1))
            goto PLSCFin;

        if (ctx->PMNS < 0)               /* number of steps */
            ctx->PMNS = 2;
        if (ctx->PMSIG <= 0.0 || ctx->PMSIG >= 1.0)
            ctx->PMSIG = 0.5;
        if (ctx->PMD < 0.0)
            ctx->PMD = 0.0;
        printf1(ctx, "Lowess: sig=%g ns=%d d=%g\n",ctx->PMSIG,ctx->PMNS,ctx->PMD);

        if (lowess(ctx, ctx->AcX - 1,ctx->AcY - 1,nn,ctx->PMSIG,ctx->PMNS,ctx->PMD,ctx->AcU,ctx->AcV,ctx->AcTmp))
            goto PLSCFin;

        first = 0;
        for (i = 0; i < nn; ++i) {
            ps_2dplot(ctx, ctx->AcX[i],ctx->AcU[i + 1],first);    
            first = 1;
        }
        fprintf(ctx->PSFd,"stroke\n");

        if (ctx->PMF1Def) {
            for (i = 0; i < nn; ++i) {
                fprintf(ctx->PMF1d,"%6d ",i + 1);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[i]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[i]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcU[i + 1]);
                fprintf(ctx->PMF1d,"\n");
                nrec++;
            }
        }
    }
    else if (ctx->PMOPT == 4) {                  /* midmeans */

        if (alloc_acu(ctx, nn + 1))
            goto PLSCFin;
        if (alloc_acv(ctx, nn + 1))
            goto PLSCFin;
        if (alloc_acxf(ctx, nn + 1))
            goto PLSCFin;
        if (alloc_acyf(ctx, nn + 1))
            goto PLSCFin;

        if (ctx->PMR < 2)
            ctx->PMR = 2;
        if (ctx->PMR >= nn) {
            printf1(ctx, "Error: r must be less than number of cases.\n");
            goto PLSCFin;
        }
        printf1(ctx, "Midmeans: r = %d\n",ctx->PMR);

        if (rmidmean(ctx, nn,ctx->AcX,ctx->AcY,ctx->PMR,ctx->AcU,ctx->AcV,ctx->AcXF,ctx->AcYF)) {
            p_err(ctx, -2,1);
            goto PLSCFin;
        }
        if (ctx->PMF1Def) {
            for (i = 0; i < nn; ++i) {
                fprintf(ctx->PMF1d,"%6d ",i + 1);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[i]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[i]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcU[i]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcV[i]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,(double)ctx->AcXF[i]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,(double)ctx->AcYF[i]);
                fprintf(ctx->PMF1d,"\n");
                nrec++;
            }
        }

        if (ctx->PMNS <= 0 || ctx->PMNS >= nn)
            ctx->PMNS = 0;
        else {
            if (alloc_actmp(ctx, nn + 1))
                goto PLSCFin;
        }
        scplot1(ctx, nn,ctx->AcU,ctx->AcV,ctx->PMNS);

        for (i = 0; i < nn; ++i)
            ctx->AcV[i] = (double)ctx->AcXF[i];
        scplot1(ctx, nn,ctx->AcU,ctx->AcV,ctx->PMNS);

        for (i = 0; i < nn; ++i)
            ctx->AcV[i] = (double)ctx->AcYF[i];
        scplot1(ctx, nn,ctx->AcU,ctx->AcV,ctx->PMNS);

    }
    fprintf(ctx->PSFd,"grestore\n");

    if (ctx->PMF1Def && nrec > 0)
        printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);

    err = 0;

PLSCFin:
    p_clean(ctx);       
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  scplot1(n,x,y,ns)   plot x,y. if ns > 0 then smooth.                    */

void scplot1(TDAContext *ctx, int n,double *x,double *y,int ns)
{
    register int i,first;

    if (ns)  
        nnsmooth(ctx, n,y,ctx->AcTmp,ctx->PMNS);
       
    first = 0;
    for (i = 0; i < n; ++i) {
        if (ns)
            ps_2dplot(ctx, x[i],ctx->AcTmp[i],first);    
        else
            ps_2dplot(ctx, x[i],y[i],first);    
        first = 1;
    }
    fprintf(ctx->PSFd,"stroke\n");
}

/*--------------------------------------------------------------------------*/
/*  scplot2     Sunflower plot                                              */
/*              n = number of points in AcX,AcY                             */
/*              return 0 if OK, -1 if insufficient memory.                  */

int scplot2(TDAContext *ctx, int n)
{
    register int i,j,ix,iy;
    int err,nx,ny,px,py;
    double ax,ay,x,y;

    err = -1;

    nx = ctx->PMNN1;
    ny = ctx->PMNN2;
    if (nx < 1)
        nx = 1;             
    if (ny < 1)
        ny = 1;           
    ax = ctx->UXLen / (double)nx;
    ay = ctx->UYLen / (double)ny;

    printf1(ctx, "Sunflower plot: nx=%d [%g] ny=%d [%g]\n",nx,ax,ny,ay);

    if (alloc_acn(ctx, nx * ny))
        goto SCP2Fin;

    for (i = 0; i < n; ++i) {
        px = (int)((ctx->AcX[i] - ctx->PA1[0]) / ax);
        py = (int)((ctx->AcY[i] - ctx->PA1[1]) / ay);
        if (px < 0)
            px = 0;
        if (px >= nx)
            px = nx - 1;
        if (py < 0)
            py = 0;
        if (py >= ny)
            py = ny - 1;
        ctx->AcN[py * nx + px] += 1;
    }
    fprintf(ctx->PSFd,"/ssiz %4.2f def\n",ctx->PtMM * ctx->PMFS / 5.0);

    /*****
    xs = 0.35 * ax * PSXLen / UXLen;    
    ys = 0.35 * ay * PSYLen / UYLen;
    ******/

    y = ctx->PA1[1] + ay / 2.0;
    for (iy = 0; iy < ny; ++iy) {
        x = ctx->PA1[0] + ax / 2.0;
        for (ix = 0; ix < nx; ++ix) {
            j = ctx->AcN[iy * nx + ix];
            if (j > 0)                  /* plot */
                scplot2a(ctx, j,x,y);
            x += ax;
        }
        y += ay;
    }
    err = 0;

SCP2Fin:
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  scplot2a(n)     create sunflower with n ticks.                          */
/*  ##                                                                      */

void scplot2a(TDAContext *ctx, int n,double x,double y)
{
    int i;
    double a,b,siz;

    x = ps_2dx(ctx, x);
    y = ps_2dy(ctx, y);
    fprintf(ctx->PSFd,"gsave\n%5.2f %5.2f translate\n",x,y);

    if (n == 1)
        ps_sym(ctx, 5,0.0,0.0,ctx->PtMM * ctx->PMFS / 5.0);      /* always symbol 5 */

    else if (n >= 2) {
        siz = ctx->PMFSS * ctx->PtMM / 2.0;    /* size in points */
        a = 2.0 * Pi / (double)n;
        b = Pi / 2.0;
        for (i = 0; i < n; ++i) {
            x = siz * cos(a * (double)i + b);
            y = siz * sin(a * (double)i + b);
            fprintf(ctx->PSFd,"%6.2f %6.2f m\n0 0 l\nstroke\n",x,y);             
        }
        /** fprintf(PSFd,"stroke\n"); **/         
    }
    fprintf(ctx->PSFd,"grestore\n");
}



