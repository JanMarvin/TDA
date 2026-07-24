/****************************************************************************/
/*  t_ireg                                                                  */
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
#include "t_gdat.h"
#include "t_gf.h"
#include "t_alloc.h"
#include "t_mat.h"
#include "t_matc.h"
#include "t_matf.h"
#include "t_sort.h"
#include "t_eval.h"  
#include "t_eval1.h"  
#include "t_eval3.h"  
#include "t_ml.h"  
#include "t_gmin.h"  
#include "t_gdd.h"  
#include "t_svd.h"  
#include "t_lp.h"  
#include "t_lsei.h"  
#include "t_com.h"  
#include "t_imat.h"  
#include "t_plot.h"  
#include "t_cplot.h"  
#include "t_min.h"  
#include "t_con.h"  
#include "tda_context.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_ireg                                                     */

int ilsreg(TDAContext *ctx);
int imreg(TDAContext *ctx);

int ivls(TDAContext *ctx);
void mcross(TDAContext *ctx, int n,int m,double *x,double *a);

int ivreg(TDAContext *ctx);
void ivreg_m(TDAContext *ctx, int opt,int np,double *ml,double *mh);
void ivreg_y(TDAContext *ctx, int opt,int np,double xml,double xmh);
int ivreg_n(TDAContext *ctx, int opt,int np);
int ivreg_upd(TDAContext *ctx, int opt,int np,int i);
double ivreg_b0(TDAContext *ctx, int np,double *xy);              
int ivreg_bb(TDAContext *ctx, int n,double *xl,double *xh,double *bl,double *bh);

int ivreg_min(TDAContext *ctx, int typ,int opt,int mopt,int narg,int nbmax,double *par,double *lb, double *ub,int mxit,double tolbw,double tolfd,double tolfe,int gc);
int ivreg_h(TDAContext *ctx, int opt,int n,double *xl,double *xh, double *rl,double *rh,int deriv,double *gl,double *gh);

int ivar1f(TDAContext *ctx, int n,double *xl,double *xh,double *ip,int *idx,double *v0,double *v1);
int ivreg_z(TDAContext *ctx, int np,int rp,int ns,int opt);

int ivreg1(TDAContext *ctx);
int ivreg1_h(TDAContext *ctx, int mopt,int opt,int n,double *xl,double *xh, double *rl,double *rh,int deriv,double *gl,double *gh,double fminp);

double i_center(TDAContext *ctx, double ac,double ar,double bc,double br);
double i_radius(TDAContext *ctx, double ac,double ar,double bc,double br);
void i_bvar(TDAContext *ctx, double xl,double xh,double bl,double bh,double *ul,double *uh);
double i_norm(TDAContext *ctx, int opt,double xl,double xh,double yl,double yh);
double ivreg1_nrm(TDAContext *ctx, int opt,double xl,double xh,double yl,double yh, double al,double ah,double bl,double bh);

int ivreg2(TDAContext *ctx);
double ivreg2m(TDAContext *ctx, double *bc,double *br);
double ivreg2f(TDAContext *ctx, double bc,double br);
int ivreg2p(TDAContext *ctx, int nlev,double *flev);
int ivr1_fn(TDAContext *ctx);
int ivreg2d(TDAContext *ctx);
void ivreg2dp(TDAContext *ctx, int n,int na,double bc,double br);


/*  Global parameters and variables                                         */



/* ------------------------------------------------------------------------ */
/*  ilsreg()    LS regression with interval-valued dependent variable.      */
/*                                                                          */
/*                  ilsreg(                                                 */
/*                      x=...,          list of x values                    */
/*                      tfmt=...,       print format, def. 10.4             */  
/*                  ) = YL,YH,X;                                            */
/*                                                                          */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int ilsreg(TDAContext *ctx)
{
    register int i,j;
    int err,il,ih,ix;     
    double x,yl,yh,ylm,yhm,v,w,xm,xxm;
    double tmp,bmin,bmax,amin,amax,ni,nw;
         
    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "LS regression with interval-valued dependent variable. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 6,4,1))     /* get parameters */
        goto ILSFin;

    if (ctx->PMOPT != 2)
        ctx->PMOPT = 1;

    if (ctx->PMNV != 3) {            /* need three variables */
        p_err(ctx, -1,1);
        goto ILSFin;
    }
    il = ctx->PMVIdx[0];
    ih = ctx->PMVIdx[1];
    ix = ctx->PMVIdx[2];

    if (alloc_acx(ctx, ctx->NOC + 1))     /* x variable */ 
        goto ILSFin;
    if (alloc_acy(ctx, ctx->NOC + 1))     /* lower bound */ 
        goto ILSFin;
    if (alloc_acw(ctx, ctx->NOC + 1))     /* width */
        goto ILSFin;

    tmp = v = xxm = xm = ylm = yhm = 0.0;
    for (i = 0; i < ctx->NOC; ++i) {
        x  = get_data(ctx, ix,i);
        yl = get_data(ctx, il,i);
        yh = get_data(ctx, ih,i);
        if (yh < yl) {
            printf1(ctx, "Error in case %d, found: %lg,%lg\n",i+1,yl,yh);
            goto ILSFin;
        }
        ctx->AcX[i] = x;
        ctx->AcY[i] = yl;
        ctx->AcW[i] = yh - yl;
        ylm += yl;
        yhm += yh;
        xm += x;
        xxm += x * x;
        v += x * yl;

/*      printf1("i=%d x=%lg  yl=%lg w=%lg \n",i,x,  yl,AcW[i]); */      
    }
    w = (double)ctx->NOC * xxm - xm * xm;
    v = (double)ctx->NOC * v - xm * ylm;
    xm /= (double)ctx->NOC;
    ylm /= (double)ctx->NOC;
    yhm /= (double)ctx->NOC;

    if (w <= 0.0) {
        printf1(ctx, "Error: indep. variable has zero variance.\n");
        goto ILSFin;
    }
    yl = yh = bmin = bmax = 0.0;
          
    for (i = 0; i < ctx->NOC; ++i) {
        tmp = ctx->AcX[i] - xm;   
        if (tmp > 0.0) {
            bmax += tmp * ctx->AcW[i];
            yl += ctx->AcY[i];
            yh += ctx->AcY[i] + ctx->AcW[i];
        }
        else if (tmp < 0.0) {
            bmin += tmp * ctx->AcW[i];
            yh += ctx->AcY[i];
            yl += ctx->AcY[i] + ctx->AcW[i];
        }
    }
    bmax = ((double)ctx->NOC * bmax + v) / w;
    bmin = ((double)ctx->NOC * bmin + v) / w;
    amax = yh / (double)ctx->NOC - xm * bmax;
    amin = yl / (double)ctx->NOC - xm * bmin;

    printf1(ctx, "\nMean of independent variable: ");
    rt_printf1_d(ctx, ctx->PMTFmtS,xm);
    printf1(ctx, "\nRange of mean of dependent variable: ");
    rt_printf1_d(ctx, ctx->PMTFmtS,ylm);
    rt_printf1_d(ctx, ctx->PMTFmtS,yhm);
    printf1(ctx, "\nBeta minimum: ");
    rt_printf1_d(ctx, ctx->PMTFmtS,bmin);
    printf1(ctx, "  corresponding alpha: ");
    rt_printf1_d(ctx, ctx->PMTFmtS,amin);
    printf1(ctx, "\nBeta maximum: ");
    rt_printf1_d(ctx, ctx->PMTFmtS,bmax);
    printf1(ctx, "  corresponding alpha: ");
    rt_printf1_d(ctx, ctx->PMTFmtS,amax);
    newline(ctx);
#ifdef TDA_R_PACKAGE
    {
        double erow[4];
        erow[0] = bmin; erow[1] = bmax;
        erow[2] = amin; erow[3] = amax;
        tda_export_mat(ctx, "iv.bounds", erow, 1, 4);
    }
#endif

    /* calculate bounds for given x values */

    if (ctx->PMNTP > 0) {

        newline(ctx);
        prnchar(ctx, ' ',ctx->PMTFmt1 - 10,0);
        printf1(ctx, "   x-value ");
        prnchar(ctx, ' ',ctx->PMTFmt1 - 10,0);
        printf1(ctx, " y minimum ");
        prnchar(ctx, ' ',ctx->PMTFmt1 - 10,0);
        printf1(ctx, " y maximum\n");
        prnchar(ctx, '-',3 * ctx->PMTFmt1 + 2,1);

        for (j = 0; j < ctx->PMNTP; ++j) {
            x = ctx->PMTP[j];
            ni = 1.0 / (double)ctx->NOC;
            nw = (double)ctx->NOC * (x - xm) / w;

            yl = yh = bmin = bmax = 0.0;

            for (i = 0; i < ctx->NOC; ++i) {
                tmp = nw * (ctx->AcX[i] - xm) + ni;

                if (tmp > 0.0)  
                    bmax += tmp * ctx->AcW[i];
                else if (tmp < 0.0)  
                    bmin += tmp * ctx->AcW[i];
            }
            yh = ylm + (x - xm) * v / w + bmax;
            yl = ylm + (x - xm) * v / w + bmin;

            rt_printf1_d(ctx, ctx->PMTFmtS,x);
            rt_printf1_d(ctx, ctx->PMTFmtS,yl);
            rt_printf1_d(ctx, ctx->PMTFmtS,yh);
            newline(ctx);

            if (ctx->PMF1Def) {
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMTFmtS,x);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMTFmtS,yl);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMTFmtS,yh);
                fprintf(ctx->PMF1d,"\n");
            }   
        }
        if (ctx->PMF1Def)  
            printf1(ctx, "\n%d records written to: %s\n",ctx->PMNTP,ctx->PMF1dName);

    }
    newline(ctx);
    err = 0;

ILSFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  imreg()     conditional means with interval-valued data.                */
/*                                                                          */  
/*                  imreg(                                                  */
/*                      fmt=...,        print format, def. 10.4             */  
/*                  ) = YL,YH,XL,XU;                                        */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int imreg(TDAContext *ctx)
{
    register int i,j;
    int err,ilx,ihx,ily,ihy,n,nn,nrec;
    double xl,xh,yl,yh,ymin,ymax,a,b;
         
    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Conditional means of interval-valued data. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,4,1))     /* get parameters */
        goto IMREGFin;

    if (ctx->PMNV != 4) {            /* need four variables */
        p_err(ctx, -1,1);
        goto IMREGFin;
    }
    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 10,4);

    ily = ctx->PMVIdx[0];
    ihy = ctx->PMVIdx[1];
    ilx = ctx->PMVIdx[2];
    ihx = ctx->PMVIdx[3];

    if (alloc_acx(ctx, ctx->NOC + 1))     /* x lower bound */ 
        goto IMREGFin;
    if (alloc_acu(ctx, ctx->NOC + 1))     /* x upper bound */ 
        goto IMREGFin;
    if (alloc_acy(ctx, ctx->NOC + 1))     /* y lower bound */
        goto IMREGFin;
    if (alloc_acv(ctx, ctx->NOC + 1))     /* y upper bound */
        goto IMREGFin;
    if (alloc_actmp(ctx, 2 * ctx->NOC + 1))     /* induced partition */
        goto IMREGFin;

    n = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        xl = get_data(ctx, ilx,i);
        xh = get_data(ctx, ihx,i);
        if (xh <= xl + ctx->EPSI1) {
            printf1(ctx, "Error in case %d, found: %lg,%lg\n",i+1,xl,xh);
            goto IMREGFin;
        }
        yl = get_data(ctx, ily,i);
        yh = get_data(ctx, ihy,i);
        if (yh <= yl + ctx->EPSI1) {
            printf1(ctx, "Error in case %d, found: %lg,%lg\n",i+1,yl,yh);
            goto IMREGFin;
        }
        ctx->AcX[i] = xl;
        ctx->AcU[i] = xh;
        ctx->AcY[i] = yl;
        ctx->AcV[i] = yh;
        ctx->AcTmp[n++] = xl;
        ctx->AcTmp[n++] = xh;
    }
    if (sortd(ctx, n,ctx->AcTmp,0))
        goto IMREGFin;
   
    nn = 1;
    for (i = 1; i < n; ++i) {           /* AcTmp contains induced partition */
        if (ctx->AcTmp[i] > ctx->AcTmp[i - 1])
            ctx->AcTmp[nn++] = ctx->AcTmp[i];
    }
    nrec = 0;

    for (i = 0; i < nn; ++i) {
        a = ctx->AcTmp[i];
        b = ctx->AcTmp[i + 1];
             
        n = 0;
        ymin = ymax = 0.0;
        for (j = 0; j < ctx->NOC; ++j) {
            if (ctx->AcX[j] <= a && ctx->AcU[j] >= b) {
                ymin += ctx->AcY[j];
                ymax += ctx->AcV[j];
                n++;
            }
        }
        if (n > 1) {
            ymin /= (double)n;
            ymax /= (double)n;
        }
        if (ctx->PMF1Def) {
            fprintf(ctx->PMF1d,"%6d ",i + 1);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,a);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ymin);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ymax);
            fprintf(ctx->PMF1d,"\n");

            fprintf(ctx->PMF1d,"%6d ",i + 1);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,b);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ymin);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ymax);
            fprintf(ctx->PMF1d,"\n");
            nrec += 2;
        }
        if (i >= nn - 2)
            break;
    }
    if (nrec > 0)
        printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);

    newline(ctx);
    err = 0;

IMREGFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ivls        Linear interval equations.                                  */
/*                                                                          */
/*              ivls(                                                       */
/*                  sc=...,         def. 0.01                               */
/*                  mxit=...,       def. 1000                               */
/*                  tolp=...        def. 1.e-8                              */
/*                  tfmt=...,       print format, def. 10.4                 */  
/*              ) = YL,YH,X1L,X1H,...;                                      */
/*                                                                          */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int ivls(TDAContext *ctx)
{
    register int i,j,k;
    int err,n,nv,m,r,iter,ni,nf,m1,nlp;
    int *z,*p,*xlidx,*xhidx;
    double tol,xl,xh,x,d,v,c,c1;
    double *acmat,*admat,*bcvec,*bdvec,*rmat,*gmat,*x0vec,*gvec,*lpmat;
    double *xlvec,*xhvec,*dvec,*d1vec,*lpx,*lpy;
         
    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    ctx->TOLP = 1.0e-8;

    printf1(ctx, "Linear interval equations. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,4,1))     /* get parameters */
        goto IVLSFin;

    nv = ctx->PMNV / 2;
    if (2 * nv != ctx->PMNV) {   
        printf1(ctx, "Error: need even number of variables on right-hand side.\n");
        goto IVLSFin;
    }
    if (nv < 2) {   
        printf1(ctx, "Error: need at least two interval-valued variables.\n");
        goto IVLSFin;
    }
    n = ctx->NOC;
    m = nv - 1;

    if (n < m) {
        printf1(ctx, "Error: number of cases must not be less than number of variables.\n");
        goto IVLSFin;
    }
    if (ctx->PMSC <= ctx->EPSI)  
        ctx->PMSC = 0.01;

    if (ctx->MxIter < 1 || ctx->MxItFlg == 0)
        ctx->MxIter = 1000;

    printf1(ctx, "\nTolerance for preliminary bounds (sc): %lg\n",ctx->PMSC);
    printf1(ctx, "Max number of iterations: %d\n",ctx->MxIter);
    printf1(ctx, "Tolerance (tolp): %lg\n\n",ctx->TOLP);

    if (alloc_acx(ctx, 2 * n * m + 4))   
        goto IVLSFin;
    acmat = ctx->AcX;
    admat = ctx->AcX + n * m + 2;

    if (alloc_acy(ctx, 3 * n + 6))   
        goto IVLSFin;
    bcvec = ctx->AcY;
    bdvec = ctx->AcY + n + 2;
    dvec  = ctx->AcY + 2 * n + 4;

    if (alloc_acz(ctx, m * n + m * m + 4))   
        goto IVLSFin;
    rmat = ctx->AcZ;
    gmat = ctx->AcZ + m * n + 2;

    if (alloc_acw(ctx, 5 * m + 10))       
        goto IVLSFin;
    x0vec = ctx->AcW;
    gvec  = ctx->AcW + m + 2;
    d1vec = ctx->AcW + 2 * m + 6;
    xlvec = ctx->AcW + 3 * m + 8;
    xhvec = ctx->AcW + 4 * m + 10;

    for (i = 0; i < n; ++i) {
        k = 0;
        for (j = 0; j < nv; ++j) {
            xl = get_data(ctx, ctx->PMVIdx[k++],i);
            xh = get_data(ctx, ctx->PMVIdx[k++],i);
            if (xh < xl) {
                printf1(ctx, "Error: no valid interval in record %d.\n",i + 1);
                goto IVLSFin;
            }
            x = (xl + xh) / 2.0;
            d = (xh - xl) / 2.0;
            if (j == 0) {
                bcvec[i] = x;                 
                bdvec[i] = d;
            }
            else {
                acmat[i * m + j] = x;
                admat[i * m + j] = d;
            }
        }
    }
/***
    for (i = 0; i < n; ++i) {
        tda_out("yc=%lf ",bcvec[i]);
        tda_out("yd=%lf  xc: ",bdvec[i]);
        for (j = 1; j <= m; ++j)  
            tda_out("%lf ",acmat[i * m + j]);
        tda_out(" xd: ");
        for (j = 1; j <= m; ++j)  
            tda_out("%lf ",admat[i * m + j]);
        newline(ctx);
    }
***/

    mcross(ctx, n,m,acmat,gmat);
    r = ginv(ctx, m,m,gmat);
    if (r < 0) {
        printf1(ctx, "Error in matrix inversion. Cannot continue.\n");
        goto IVLSFin;
    }
    else if (r != m) {
        printf1(ctx, "Error: insufficient rank=%d of data matrix.\n",r);
        goto IVLSFin;
    }

    /* calculate R matrix (m,n) */

    for (i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
            x = 0.0;
            for (k = 1; k <= m; ++k)
                x += gmat[i * m + k] * acmat[j * m + k];
            rmat[i * n + j + 1] = x;
        }
    }

    /* calculate G matrix (m,m) */

    for (i = 0; i < m; ++i) {
        for (j = 1; j <= m; ++j) {
            x = 0.0;
            for (k = 0; k < n; ++k) 
                x += fabs(rmat[i * n + k + 1]) * admat[k * m + j];
            gmat[i * m + j] = x;
        }
    }

    /* calculate x0 vector (m) in x0vec */

    for (i = 0; i < m; ++i) {
        x = 0.0;
        for (k = 1; k <= n; ++k)
            x += rmat[i * n + k] * bcvec[k - 1];
        x0vec[i] = x;
    }

    /* calculate g vector (m) gvec */

    for (i = 0; i < n; ++i) {
        d = bdvec[i];
        for (j = 1; j <= m; ++j)
            d += admat[i * m + j] * fabs(x0vec[j - 1]);
        dvec[i] = d;
    }
    for (i = 0; i < m; ++i) {
        x = 0.0;
        for (k = 0; k < n; ++k) 
            x += fabs(rmat[i * n + k + 1]) * dvec[k];
        gvec[i] = x;
    }

    /* iteration for calculating solution of Gd + g < d */
      
    tol = ctx->DBLMAX / 1000.0;

    for (i = 0; i < m; ++i)
        d1vec[i] = 0.0;

    for (iter = 0; iter < ctx->MxIter; ++iter) {

        for (i = 0; i < m; ++i)
            dvec[i] = d1vec[i];

        for (i = 0; i < m; ++i) {
            x = 0.0;
            for (j = 1; j <= m; ++j)  
                x += gmat[i * m + j] * dvec[j - 1];
            d1vec[i] = x + gvec[i] + ctx->PMSC;
        }
        d = 0.0;
        for (i = 0; i < m; ++i)  
            d = dmax(ctx, d,fabs(dvec[i] - d1vec[i]));

        if (d < ctx->PMSC || d > tol)
            break;
    }
    if (d >= ctx->PMSC) {
        printf1(ctx, "Error: no convergence (d=%lg).\n",d);
        goto IVLSFin;
    }

    /* calculate preliminary bounds in xlvec and xhvec */

    for (i = 0; i < m; ++i) {
        x = 0.0;
        for (k = 1; k <= n; ++k)
            x += rmat[i * n + k] * bcvec[k - 1];
        xlvec[i] = x - d1vec[i];
        xhvec[i] = x + d1vec[i];
    }
    printf1(ctx, "Preliminary bounds\n");
    for (i = 0; i < m; ++i) {
        rt_printf1_d(ctx, ctx->PMTFmtS,xlvec[i]);
        rt_printf1_d(ctx, ctx->PMTFmtS,xhvec[i]);
        newline(ctx);
    }
    newline(ctx);

    /* try to find optimal enclosure */

    m1 = m + 1;

    if (alloc_acu(ctx, (2 * n + 1) * m1 + 2))   
        goto IVLSFin;
    lpmat = ctx->AcU;

    if (alloc_acv(ctx, m + 2 * n + 4))   
        goto IVLSFin;
    lpx = ctx->AcV;
    lpy = ctx->AcV + m + 2;

    if (alloc_acn(ctx, 2 * m + 4))   
        goto IVLSFin;
    z = ctx->AcN;
    p = ctx->AcN + m + 2;

    if (alloc_acm(ctx, 2 * m + 4))   
        goto IVLSFin;
    xlidx = ctx->AcM;
    xhidx = ctx->AcM + m + 2;
             
    nlp = 1;
    ni = 0;
    for (i = 0; i < m; ++i) {
        if (xlvec[i] > 0.0)
            z[i] = 1;
        else if (xhvec[i] < 0.0)
            z[i] = -1;
        else {
            p[ni++] = i;
            nlp *= 2;
        }
    }
    nlp *= 2 * m;
    printf1(ctx, "Trying to find optimal enclosure.\n");
    printf1(ctx, "Need to solve %d lp problems.\n\n",nlp);
             
    for (i = 0; i < ni; ++i)
        z[p[i]] = -1;
                    
    while (1) {
        /***/
        tda_out("z: ");
        for (i = 0; i < m; ++i)
            tda_out("%2d ",z[i]);
        newline(ctx);
        /***/

        /* solve lp problems */
     
        for (i = 1; i <= n; ++i) {
            lpmat[i * m1 + m1] = bcvec[i - 1] + bdvec[i - 1];
            lpmat[(i + n) * m1 + m1] = bdvec[i - 1] - bcvec[i - 1];
           
            for (j = 1; j <= m; ++j) {
                x = acmat[(i - 1) * m + j];
                if (z[j - 1] < 0)
                    x = -x;

                lpmat[i * m1 + j] = x - admat[(i - 1) * m + j];
                lpmat[(i + n) * m1 + j] = -(x + admat[(i - 1) * m + j]);
            }         
        }
       
        for (k = 0; k < m; ++k) {

            if (z[k] > 0) {
                c = -1.0;
                c1 = 1.0;
            }
            else {
                c = 1.0;
                c1 = -1.0;
            }
            lpmat[k + 1] = c;

            /***
            tda_out("\nlpmat fuer min\n");
            for (i = 0; i <= 2 * n; ++i) {
                for (j = 1; j <= m1; ++j)
                    tda_out("%lf ",lpmat[i * m1 + j]);
                newline(ctx);
            }
            ***/
         
            r = lpf1(ctx, 2 * n,m,0,lpmat,lpx,lpy,&v,ctx->TOLP);       

            /***/
            tda_out("min  r=%d\n",r);
            tda_out("lpx: ");
            for (i = 1; i <= m; ++i)
                tda_out("%g ",lpx[i]);
            newline(ctx);
            /***/

            if (r == 0) {
                if (xlidx[k] == 0) {
                    xlvec[k] = lpx[k + 1] * c1;
                    xlidx[k] = 1;
                }
                else
                    xlvec[k] = dmin(ctx, xlvec[k],lpx[k + 1] * c1);
            }
            lpmat[k + 1] = c1;

            /***
            tda_out("\nlpmat fuer max\n");
            for (i = 0; i <= 2 * n; ++i) {
                for (j = 1; j <= m1; ++j)
                    tda_out("%f ",lpmat[i * m1 + j]);
                newline(ctx);
            }
            ***/

            r = lpf1(ctx, 2 * n,m,0,lpmat,lpx,lpy,&v,ctx->TOLP);       

            /***/
            tda_out("max r=%d\n",r);
            tda_out("lpx: ");
            for (i = 1; i <= m; ++i)
                tda_out("%g ",lpx[i]);
            newline(ctx);
            /***/

            if (r == 0) {
                if (xhidx[k] == 0) {
                    xhvec[k] = lpx[k + 1] * c1;
                    xhidx[k] = 1;
                }
                else
                    xhvec[k] = dmax(ctx, xhvec[k],lpx[k + 1] * c1);
            }
            lpmat[k + 1] = 0.0;
        }
        if (ni == 0)
            break;

        nf = 0;
        for (i = 0; i < ni; ++i) {
            k = p[i];
            if (z[k] == -1) {
                z[k] = 1;
                break;
            }
            else {   
                if (i == ni - 1) {
                    nf = 1;
                    break;
                }
                z[k] = -1;
            }
        }
        if (nf)
            break;
    }
    nf = 1;
    for (i = 0; i < m; ++i) {
        if (xlidx[i] == 0 || xhidx[i] == 0) {
            nf = 0;
            break;
        }
    }
    if (nf == 0) {
        printf1(ctx, "Cannot find a solution.\n");
    }   
    else {
        nf = 1;
        printf1(ctx, "Final bounds\n");
        for (i = 0; i < m; ++i) {
            rt_printf1_d(ctx, ctx->PMTFmtS,xlvec[i]);
            rt_printf1_d(ctx, ctx->PMTFmtS,xhvec[i]);
            newline(ctx);
            if (xlvec[i] > xhvec[i] + ctx->EPSI)
                nf = 0;
        }
        if (nf == 0)
            printf1(ctx, "\nThe solution set seems to be empty.\n");
    }
    newline(ctx);
    err = 0;

IVLSFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  mcross(n,m,x,a)  Return cross product x'x in a.                         */

void mcross(TDAContext *ctx, int n,int m,double *x,double *a)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j,k;
    double tmp;

    for (i = 1; i <= m; ++i) {
        for (j = 1; j <= m; ++j) {
            tmp = 0.0;
            for (k = 0; k < n; ++k)
                tmp += x[k * m + i] * x[k * m + j];
            a[(i - 1) * m + j] = tmp;
        }
    }
}

/* ----------------------------------------------------------------------- */
/*  ivreg       LS regression with interval data, for a single              */
/*              independent variables (intercept added).                    */
/*                                                                          */
/*              ivreg(                                                      */
/*                  opt=...,        option, def. 1                          */
/*                                  1 = only first two steps                */
/*                                  2 = heuristic optimization              */
/*                                  3 = exact optimization                  */
/*                  ns=...,         min or max, def. 0                      */
/*                                  0 = min                                 */
/*                                  1 = max                                 */
/*                  n=...,          block size, opt 2                       */
/*                  nbox=...,       max number of boxes, def. 100           */
/*                  mxit=...,       max number of iterations, def. 100      */
/*                  tolbw=...,      tolerance box widght, def. 1.e-4        */
/*                  tolf=...,       tolerance for beta iteration, 1.e-6     */
/*                  tolfd=...,      tolerance function range, def. 1.e-10   */
/*                  tolbw=...,      tolerance global min/max, def. 1.e-10   */
/*                  fmt=...,        print format, def. 13.6                 */  
/*                  df=...,         output file containing data             */
/*                  prot=...,       protocol file                           */
/*              ) = YL,YH,XL,XH;                                            */
/*                                                                          */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int ivreg(TDAContext *ctx)
{
    register int i,j,k;
    int err,n,np,r,nx = 0,ny,nx0,ny0,iter,rp;
    double xl,xh;
           
    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    ctx->TOLF  = 1.e-6;       /* def tolerance for beta iteration                 */
    ctx->TOLBW = 1.e-4;       /* def tolerance for box length                     */
    ctx->TOLFD = 1.e-10;      /* def tolerance for function range                 */
    ctx->TOLFE = 1.e-10;      /* def tolerance for global minimum                 */

    printf1(ctx, "Regression with interval data. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,4,1))     /* get parameters */
        goto IVREGFin;

    if (ctx->PMOPT < 1 || ctx->PMOPT > 3)
        ctx->PMOPT = 1;

    if (ctx->PMN < 2)
        ctx->PMN = 20;

    if (ctx->PMNV != 4) {   
        printf1(ctx, "Error: need four variables on right-hand side.\n");
        goto IVREGFin;
    }
    if (ctx->PMNS != 1)
        ctx->PMNS = 0;

    if (ctx->PMFmtF == 0)        /* default print format */
        pmfmt(ctx, 13,6);

    n = ctx->NOC;
    np = 2 * n; 

    if (ctx->PMNS == 0)
        printf1(ctx, "\nMinimization.\n");
    else
        printf1(ctx, "\nMaximization.\n");
    printf1(ctx, "Number of data points: %d\n",n);
    if (n < 2) {
        printf1(ctx, "Error: need at least two data points.\n");
        goto IVREGFin;
    }
    printf1(ctx, "Tolerance for box width: %g\n\n",ctx->TOLBW);

    /* allocate RPar, RParL, RParU (i = 0,...,np-1) */

    if (alloc_par(ctx, np,1))
        goto IVREGFin;

    /* temporarily used for variance calculation */

    if (alloc_acx(ctx, n + 1))              
        goto IVREGFin; 
    if (alloc_acy(ctx, n + 1))  
        goto IVREGFin; 
    if (alloc_actmp(ctx, 2 * n + 1))  
        goto IVREGFin; 
    if (alloc_ack(ctx, n + 1))  
        goto IVREGFin; 

    /*  get data into RParL and RParU */

    k = 0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < 2; ++j) {
            xl = get_data(ctx, ctx->PMVIdx[j * 2],i);
            xh = get_data(ctx, ctx->PMVIdx[j * 2 + 1],i);
            if (xh < xl) {
                printf1(ctx, "Error: no valid interval in record %d.\n",i + 1);
                goto IVREGFin;
            }
            ctx->RParL[k] = xl;
            ctx->RParU[k] = xh;
            ctx->RPar[k] = (xl + xh) / 2.0;
            k++;

            if (j == 1) {
                ctx->AcX[i] = xl;
                ctx->AcY[i] = xh;
            }
        }
    }
/**
    tda_out("data RPar, RParl and RParU\n"); 
    for (i = 0; i < np; ++i) {
        tda_out("%3d %f %f %f\n",i,RPar[i],RParL[i],RParU[i]);
    }
**/
    /* ============================================================ */

    ivreg_m(ctx, 0,np,&ctx->IVRXML,&ctx->IVRXMU);      /* calculate mean of X */

    ivreg_y(ctx, ctx->PMNS,np,ctx->IVRXML,ctx->IVRXMU);     /* fix y values */

    ivreg_m(ctx, 1,np,&ctx->IVRYML,&ctx->IVRYMU);      /* calculate mean of Y */

    ny = ivreg_n(ctx, 1,np);                 /* number of fixed y values */

    printf1(ctx, "First step: fixed %d values of Y.\n",ny);
    printf1(ctx, "Mean Y = [ ");
    rt_printf1_d(ctx, ctx->PMFmtS,ctx->IVRYML);  
    printf1(ctx, ", ");
    rt_printf1_d(ctx, ctx->PMFmtS,ctx->IVRYMU);  
    tda_out("]\nMean X = [ ");
    rt_printf1_d(ctx, ctx->PMFmtS,ctx->IVRXML);  
    printf1(ctx, ", ");
    rt_printf1_d(ctx, ctx->PMFmtS,ctx->IVRXMU);  
    tda_out("]\n");

    /* calculate variance of X */

    r = ivar1f(ctx, n,ctx->AcX,ctx->AcY,ctx->AcTmp,ctx->AcK,&ctx->IVRVARL,&ctx->IVRVARU);
    if (r) {
        printf1(ctx, "\nError: cannot find bounds for variance of X.\n");
        goto IVREGFin;
    }
    printf1(ctx, "Var  X = [ ");
    rt_printf1_d(ctx, ctx->PMFmtS,ctx->IVRVARL);  
    printf1(ctx, ", ");
    rt_printf1_d(ctx, ctx->PMFmtS,ctx->IVRVARU);  
    printf1(ctx, "]\n\n");
#ifdef TDA_R_PACKAGE
    /* the three bracketed ranges printed just above, in the order they
       appear: mean of Y, mean of X, variance of X */
    {
        double erow[6];
        erow[0] = ctx->IVRYML;   erow[1] = ctx->IVRYMU;
        erow[2] = ctx->IVRXML;   erow[3] = ctx->IVRXMU;
        erow[4] = ctx->IVRVARL;  erow[5] = ctx->IVRVARU;
        tda_export_mat(ctx, "iv.ranges", erow, 1, 6);
    }
#endif

/**
    tda_out("data RPar, RParl and RParU\n"); 
    for (i = 0; i < np; ++i) {
        tda_out("%3d %f %f %f\n",i,RPar[i],RParL[i],RParU[i]);
    }
**/
    /* ============================================================ */

    printf1(ctx, "Second step: fixing values of X (and Y).\n");

    ctx->IVRBETA = ivreg_b0(ctx, np,ctx->RPar);                /* preliminary beta */
    ivreg_bb(ctx, np,ctx->RParL,ctx->RParU,&ctx->IVRBL,&ctx->IVRBU);     /* bounds for beta */

    if (ctx->PMNS == 0)                              /* minimization */
        ctx->IVRBU = ctx->IVRBETA;
    else    
        ctx->IVRBL = ctx->IVRBETA;

    printf1(ctx, "Iteration     bounds of current beta          bounds of x-mean  x-fixed  y-fixed\n");
    printf1(ctx, "%7d  %14.6f %12.6f %12.6f %12.6f %8d %8d\n",0,ctx->IVRBL,ctx->IVRBU,ctx->IVRXML,ctx->IVRXMU,0,ny);
           
    iter = nx0 = 0;
    ny0 = ny;
            
    while (1) {                                 

        for (i = 0; i < np; i += 2) {
            if (fabs(ctx->RParU[i + 1] - ctx->RParL[i + 1]) <= ctx->TOLBW)
                continue;

            ivreg_upd(ctx, ctx->PMNS,np,i);
        }
        if (nx0 < n) {
            ivreg_m(ctx, 0,np,&ctx->IVRXML,&ctx->IVRXMU);      /* calculate new mean of X */
            nx = ivreg_n(ctx, 0,np);                 /* number of fixed x values */
        }
        if (ny0 < n) {                          /* also update y values */
            ivreg_y(ctx, ctx->PMNS,np,ctx->IVRXML,ctx->IVRXMU);     /* fix y values */
            ivreg_m(ctx, 1,np,&ctx->IVRYML,&ctx->IVRYMU);      /* calculate mean of Y */
            ny = ivreg_n(ctx, 1,np);                 /* number of fixed y values */
        }
        if (nx == nx0 && ny == ny0)  
            break;
        nx0 = nx;
        ny0 = ny;

        /* try to update variance of X */

        k = 0;
        for (i = 0; i < np; i += 2) {
            ctx->AcX[k] = ctx->RParL[i + 1];
            ctx->AcY[k] = ctx->RParU[i + 1];
            k++;
        }
        ivar1f(ctx, n,ctx->AcX,ctx->AcY,ctx->AcTmp,ctx->AcK,&ctx->IVRVARL,&ctx->IVRVARU);            
/*
        tda_out("new var %g %g\n",IVRVARL,IVRVARU);
*/
        /* update beta */

        ctx->IVRBETA = ivreg_b0(ctx, np,ctx->RPar);                
        ivreg_bb(ctx, np,ctx->RParL,ctx->RParU,&ctx->IVRBL,&ctx->IVRBU);  
        if (ctx->PMNS == 0)              
            ctx->IVRBU = ctx->IVRBETA;
        else    
            ctx->IVRBL = ctx->IVRBETA;

        printf1(ctx, "%7d  %14.6f %12.6f %12.6f %12.6f %8d %8d\n",++iter,ctx->IVRBL,ctx->IVRBU,ctx->IVRXML,ctx->IVRXMU,nx,ny);
#ifdef TDA_R_PACKAGE
        /* one row per iteration; the reader takes the last, which is
           the converged beta range */
        {
            double erow[6];
            erow[0] = (double)iter;
            erow[1] = ctx->IVRBL;   erow[2] = ctx->IVRBU;
            erow[3] = ctx->IVRXML;  erow[4] = ctx->IVRXMU;
            erow[5] = (double)nx;
            tda_export_row(ctx, "iv.iterations", erow, 6);
        }
#endif



    }
           
    /* calculate rp = remaining parameters */

    rp = 0; 
    for (i = 0; i < np; i += 2) {
        if (fabs(ctx->RParU[i] - ctx->RParL[i]) > ctx->TOLBW)  
            rp++;
        if (fabs(ctx->RParU[i + 1] - ctx->RParL[i + 1]) > ctx->TOLBW)  
            rp++;
    }
    printf1(ctx, "\nRemaining coordinates: %d\n",rp);
              
    /* ============================================================ */

    if (ctx->PMOPT == 1 || rp == 0)
        goto IVREGCont; 
                   
    printf1(ctx, "\nThird step: additional optimization (opt=%d).\n",ctx->PMOPT);

    if (ctx->PMOPT == 2) {
        ivreg_z(ctx, np,rp,ctx->PMN,ctx->PMNS);
        goto IVREGFin;

    }
    else if (ctx->PMOPT == 3) {

        /* The shipped defaults were a flat mxit=100, nbox=100, which
           made the exact search give up even on problems it certifies
           in well under a second (7 free coordinates need ~1500
           iterations; 13 need ~500000).  When the user has not set the
           limits, scale them to the number of free coordinates instead:
           enough to certify what is certifiable in reasonable time,
           capped so hopeless problems still terminate.  User-set
           values are respected exactly as before.  Deliberate change,
           made on request -- see README. */

        if (ctx->MxIter < 1 || ctx->MxItFlg == 0) {
            ctx->MxIter = 20000;
            for (i = 7; i < rp && ctx->MxIter < 1000000; ++i)
                ctx->MxIter *= 2;
            if (ctx->MxIter > 1000000)
                ctx->MxIter = 1000000;
        }
        if (ctx->PMNBOX < 1) {       /* max number of boxes */
            ctx->PMNBOX = ctx->MxIter / 10;
            if (ctx->PMNBOX < 1000)
                ctx->PMNBOX = 1000;
            if (ctx->PMNBOX > 100000)
                ctx->PMNBOX = 100000;
        }

        /* reorganize x- and y-values */

        if (alloc_acx(ctx, n + 1))  
            goto IVREGFin; 
        if (alloc_acy(ctx, n + 1))  
            goto IVREGFin; 
        if (alloc_aci(ctx, n + 1))  
            goto IVREGFin; 
        if (alloc_acj(ctx, n + 1))  
            goto IVREGFin; 

        ctx->IVRX  = ctx->AcX;
        ctx->IVRY  = ctx->AcY;
        ctx->IVRXP = ctx->AcI;
        ctx->IVRYP = ctx->AcJ;

        j = k = 0;
        for (i = 0; i < np; i += 2) {

            ctx->IVRY[k] = ctx->RPar[i];
            ctx->IVRX[k] = ctx->RPar[i + 1];

            ctx->IVRYP[k] = -1;
            if (fabs(ctx->RParU[i] - ctx->RParL[i]) > ctx->TOLBW) {
                ctx->IVRYP[k] = j;
                ctx->RParL[j] = ctx->RParL[i];
                ctx->RParU[j] = ctx->RParU[i];
                ctx->RPar[j]  = ctx->RPar[i];
                j++;
            }
            ctx->IVRXP[k] = -1;
            if (fabs(ctx->RParU[i + 1] - ctx->RParL[i + 1]) > ctx->TOLBW) {
                ctx->IVRXP[k] = j;
                ctx->RParL[j] = ctx->RParL[i + 1];
                ctx->RParU[j] = ctx->RParU[i + 1];
                ctx->RPar[j]  = ctx->RPar[i + 1];
                j++;
            }
            k++;
        }
        if (alloc_list(ctx, ctx->PMNBOX,rp)) 
            goto IVREGFin;

        printf1(ctx, "Number of parameters: %d\n",rp);              
        printf1(ctx, "Maximum number of iterations: %d\n",ctx->MxIter);              
        printf1(ctx, "Maximum number of boxes: %d\n",ctx->PMNBOX);
        printf1(ctx, "Tolerance for box width: %g\n",ctx->TOLBW);
        printf1(ctx, "Tolerance for function range: %g\n",ctx->TOLFD);
        printf1(ctx, "Tolerance for global minimum: %g\n\n",ctx->TOLFE);
      
        r = ivreg_min(ctx, 0,ctx->PMNS,0,rp,ctx->PMNBOX,ctx->RPar,ctx->RParL,ctx->RParU,ctx->MxIter,ctx->TOLBW,ctx->TOLFD,ctx->TOLFE,1);
        if (r) {
            if (r == 1)
                printf1(ctx, "Exceeded maximum number of boxes.\n");
            printf1(ctx, "The value below is the best point found, not proven\n");
            printf1(ctx, "optimal; raise nbox= (and mxit=) to certify it.\n");
            printf1(ctx, "Current best value of beta:");
            rt_printf1_d(ctx, ctx->PMFmtS,ctx->GO_FMIN);
            newline(ctx);
            goto IVREGFin;
        }
        igmin_res(ctx, ctx->TOLFE,0);
    }    

IVREGCont:
    if (ctx->PMF1Def) {          /* write to output file */

        if (ctx->PMOPT == 1 || rp == 0) {
            for (i = 0; i < np; i += 2) {
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->RPar[i]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->RPar[i + 1]);
                fprintf(ctx->PMF1d,"\n");
            }
        }
        else if (ctx->PMOPT == 3) {   

            for (i = 0; i < n; ++i) {
                if ((k = ctx->IVRYP[i]) >= 0)
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->RPar[k]);
                else
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->IVRY[i]);

                if ((k = ctx->IVRXP[i]) >= 0)
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->RPar[k]);
                else
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->IVRX[i]);
                fprintf(ctx->PMF1d,"\n");
            }                     
        }
        printf1(ctx, "%d records written to: %s\n",n,ctx->PMF1dName);
    }
    newline(ctx);
    err = 0;

IVREGFin:
    alloc_list(ctx, 0,0);
    alloc_par(ctx, 0,0);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ivreg_m(opt,np,ml,mh)                                                   */
/*                                                                          */
/*  Calculate bounds of mean value. If opt = 0 then X, otherwise Y.         */
/*  Return bounds in ml and mh.                                             */

void ivreg_m(TDAContext *ctx, int opt,int np,double *ml,double *mh)
{
    register int i;

    *ml = *mh = 0.0;
    for (i = 0; i < np; i += 2) {
        if (opt == 1) {
            *ml += ctx->RParL[i];
            *mh += ctx->RParU[i];
        }
        else {
            *ml += ctx->RParL[i + 1];
            *mh += ctx->RParU[i + 1];
        }
    }
    *ml /= (double)(np / 2);
    *mh /= (double)(np / 2);
}

/* ------------------------------------------------------------------------ */
/*  ivreg_y(opt,np,xml,xmh)                                                 */
/*                                                                          */
/*  Fix values of y variable. Use global variables RPar, RParL, RParU.      */
/*  xml and xmh contain bounds of mean of X variable.                       */
/*  If opt = 0 minimization, otherwise maximization.                        */

void ivreg_y(TDAContext *ctx, int opt,int np,double xml,double xmh)
{
    register int i;

    for (i = 0; i < np; i += 2) {
        if (ctx->RParU[i + 1] < xml) {
            if (opt == 0)   
                ctx->RPar[i] = ctx->RParL[i] = ctx->RParU[i];
            else
                ctx->RPar[i] = ctx->RParU[i] = ctx->RParL[i];
        }
        else if (ctx->RParL[i + 1] > xmh) {
            if (opt == 0)   
                ctx->RPar[i] = ctx->RParU[i] = ctx->RParL[i];
            else
                ctx->RPar[i] = ctx->RParL[i] = ctx->RParU[i];
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  ivreg_n(opt,np)     If opt = 0 return number of fixed x values,         */
/*                      otherwise number of fixed y values.                 */

int ivreg_n(TDAContext *ctx, int opt,int np)
{
    register int i,n;

    n = 0;
    for (i = 0; i < np; i += 2) {
        if (opt == 1) {
            if (fabs(ctx->RParU[i] - ctx->RParL[i]) <= ctx->TOLBW)
                n++;          
        }
        else {
            if (fabs(ctx->RParU[i + 1] - ctx->RParL[i + 1]) <= ctx->TOLBW)
                n++;          
        }
    }
    return(n);
}

/* ------------------------------------------------------------------------ */
/*  ivreg_upd(opt,np,i)                                                     */
/*                                                                          */
/*  Update x-parameter i+1. If opt=0 minimization, otherwise maximization.  */
/*  Use current bounds for beta: IVRBL and IVRBU.                           */
/*  Return 1 if successful update, otherwise 0.                             */

int ivreg_upd(TDAContext *ctx, int opt,int np,int i)
{
    (void)np;        /* unused: the signature is shared */
    double xl,xh,yl,yh,tl,th;

    i_sub(ctx, ctx->RParL[i],ctx->RParU[i],ctx->IVRYML,ctx->IVRYMU,&yl,&yh);
    i_sub(ctx, ctx->RParL[i + 1],ctx->RParU[i + 1],ctx->IVRXML,ctx->IVRXMU,&xl,&xh);
    i_mul(ctx, 2.0 * ctx->IVRBL,2.0 * ctx->IVRBU,xl,xh,&tl,&th);
       
    if (yh < tl) {
        if (opt == 0)                      
            ctx->RPar[i + 1] = ctx->RParL[i + 1] = ctx->RParU[i + 1];
        else
            ctx->RPar[i + 1] = ctx->RParU[i + 1] = ctx->RParL[i + 1];
        return(1);
    }
    else if (yl > th) {
        if (opt == 0)                      
            ctx->RPar[i + 1] = ctx->RParU[i + 1] = ctx->RParL[i + 1];
        else
            ctx->RPar[i + 1] = ctx->RParL[i + 1] = ctx->RParU[i + 1];
        return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  ivreg_b0(np,xy)     Return beta based on interval centers.              */

double ivreg_b0(TDAContext *ctx, int np,double *xy)               
{
    register int i;   
    double n,tmp,sx,sy,sxx,sxy;
       
    n = (double)(np / 2);

    sx = sy = sxx = sxy = 0.0;
    for (i = 0; i < np; i += 2) {
        sx += xy[i + 1];
        sy += xy[i];
        sxx += xy[i + 1] * xy[i + 1];
        sxy += xy[i] * xy[i + 1];
    }
    tmp = n * sxx - sx * sx;
    if (tmp == 0.0) {
        printf1(ctx, "ERROR in ivreg_b0. Zero variance.\n");
        exit(0);
    }
    return((n * sxy - sx * sy) / tmp);
}

/* ------------------------------------------------------------------------ */
/*  ivreg_bb(n,xl,xh,bl,bh)                                                 */
/*                                                                          */
/*  Calculate bounds for beta. Return in bl and bh.                         */
/*                                                                          */
/*  Return 0 if OK, -1 if error                                             */

int ivreg_bb(TDAContext *ctx, int n,double *xl,double *xh,double *bl,double *bh)
{
    register int i;   
    double nn,sl,sh,tl,th,ul,uh;
       
    nn = (double)(n / 2);
    ul = uh = 0.0;

    for (i = 0; i < n; i += 2) {
        i_sub(ctx, nn * xl[i],nn * xh[i],nn * ctx->IVRYML,nn * ctx->IVRYMU,&sl,&sh);
        i_mul(ctx, sl,sh,xl[i + 1],xh[i + 1],&tl,&th);
        i_add(ctx, ul,uh,tl,th,&sl,&sh);
        ul = sl;
        uh = sh;
    }
    if (i_div(ctx, ul,uh,nn * nn * ctx->IVRVARL,nn * nn * ctx->IVRVARU,bl,bh)) {
        printf1(ctx, "FATAL ERROR in ivreg_bb.\n");
        exit(0);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  ivreg_min(typ,opt,mopt,narg,nbmax,par,lb,ub,mxit,tolbw,tolfd,tolfe,gc)  */
/*                                                                          */
/*  Branch and bound bisection for ivreg command.                           */
/*                                                                          */
/*              typ         0 ivreg                                         */
/*                          1 ivreg1                                        */
/*              opt         0 if minimization                               */
/*                          1 if maximization                               */
/*              mopt        only for typ 1                                  */
/*                          0 absolute, 1 squared                           */
/*              narg        number of arguments                             */
/*              nbmax       max number of boxes                             */
/*              par         initial parameters, par[i], i = 1,...,narg      */
/*              lb          initial lower bounds                            */
/*              ub          initial upper bounds                            */
/*              mxit        max number of iterations                        */
/*              tolbw       tolerance for box width                         */
/*              tolfd       tolerance for function range                    */
/*              tolfe       tolerance for global minimum                    */
/*              gc          0 don't use gradients, 1 use gradients          */
/*                                                                          */
/*  It is asssumed that a list for nbmax boxes is already allocated:        */
/*                                                                          */
/*      GO_LB       lower bounds                                            */
/*      GO_UB       upper bounds                                            */
/*      GO_LBF      lower function values                                   */
/*      GO_UBF      upper function values                                   */
/*      GO_FLG      -1 if not yet processed                                 */
/*                   0 if dropped                                           */
/*                   1 if temporarily accepted                              */
/*                   2 if finally accepted                                  */
/*                                                                          */
/*  The ivreg_min function sets the following global variables:             */
/*                                                                          */
/*  GO_FN       number of function evaluations                              */
/*  GO_IFN      number of inclusion function evalutaions                    */
/*  GO_IT       number of iterations performed                              */
/*  GO_NB       number of boxes in final list                               */
/*  GO_NBU      number of boxes used                                        */
/*  GO_FMIN     best function value                                         */
/*                                                                          */  
/*  The parameters for the best function value are returned in par[].       */
/*                                                                          */
/*  Return  0 if successful                                                 */
/*         -1 if insuff memory                                              */
/*         -2 error in function evaluation                                  */ 
/*          1 if exceeded max number of boxes                               */

int ivreg_min(TDAContext *ctx, int typ,int opt,int mopt,int narg,int nbmax,double *par,double *lb, double *ub,int mxit,double tolbw,double tolfd,double tolfe,int gc)
{
    (void)tolbw;        /* unused: the signature is shared */
    register int i,j;
    int err,r,nb,is,iter,ja,jb,jj,first;
    double fl,fu,tmp,w,f,fminp;
    double *lptra,*lptrb,*uptra,*uptrb;

    fminp = ctx->DBLMAX;
    ctx->GO_IT = 0;              /* number of iterations performed */
    ctx->GO_FN = 0;              /* number of function evaluations */
    ctx->GO_IFN = 0;             /* number of inclusion function evaluations */
    ctx->GO_NBU = 0;             /* number of boxes used */
    ctx->GO_NB = 0;              /* number of boxes in final list */

    err = -1;

    for (i = 0; i < nbmax; ++i)
        ctx->GO_FLG[i] = 0;

    if (alloc_actmp(ctx, narg + 1))  
        goto IGMINFin; 
    if (gc) {
        if (alloc_acu(ctx, narg + 1))    /* lower bounds of gradient */
            goto IGMINFin; 
        if (alloc_acv(ctx, narg + 1))    /* upper bounds of gradient */
            goto IGMINFin; 
    }
    err = 0;

    lptra = ctx->GO_LB;                      /* put initial box on list */
    uptra = ctx->GO_UB;

    for (i = 0; i < narg; ++i) {
        lptra[i] = lb[i];
        uptra[i] = ub[i];
    }
    ctx->GO_FLG[0] = -1;     /* not yet processed */

    /* calculate inclusion function with initial boxes */
   
    ctx->GO_IFN++;

    if (typ == 0)
        r = ivreg_h(ctx, 1,narg,lptra,uptra,&fl,&fu,0,ctx->AcU,ctx->AcV);                                      
    else
        r = ivreg1_h(ctx, mopt,1,narg,lptra,uptra,&fl,&fu,0,ctx->AcU,ctx->AcV,fminp);                                      

    if (r) {
        printf1(ctx, "Error in evaluating inclusion function.\n");
        err = -2;
        goto IGMINFin;
    }
    ctx->GO_LBF[0] = fl;
    ctx->GO_UBF[0] = fu;
    
    /* calculate function at starting values */

    ctx->GO_FN++;
    if (typ == 0)
        r = ivreg_h(ctx, 0,narg,par,&tmp,&fminp,&tmp,0,ctx->AcU,ctx->AcV);                                      
    else
        r = ivreg1_h(ctx, mopt,0,narg,par,&tmp,&fminp,&tmp,0,ctx->AcU,ctx->AcV,fminp);                                      

    if (r) {
        printf1(ctx, "Error in evaluating function.\n");
        err = -2;
        goto IGMINFin;
    }
    nb = 1;     /* number of boxes */
    iter = 0;                                     

    if (ctx->SILENTFlg < 2)
        printfe(ctx, "\n  Iter    Function Value       NBox  FCall  IFCall\n");
    
    while (++iter <= mxit) {       
#ifdef TDA_R_PACKAGE
        /* an interrupt behaves exactly like exhausting the iteration
           limit: the loop ends, the ordinary bookkeeping runs, and the
           result is reported as uncertified with this line saying why */
        if (tda_check_interrupt()) {
            printf1(ctx, "Interrupted by user.\n");
            break;
        }
#endif

        if (ctx->SILENTFlg < 2)
            printfe(ctx, "%5d  %20.13e  %6d %6d %7d\n",iter,fminp,nb,ctx->GO_FN,ctx->GO_IFN);

        /* find box with lowest lower, or largest upper, function bound,
           and check for boxes that can be dropped */

        ja = -1;
        first = 1;
        for (i = 0; i < nb; ++i) {

            if (ctx->GO_FLG[i]) {    /* check all boxes */
  
                if (opt == 0) {
                    if (ctx->GO_LBF[i] > fminp)  
                        ctx->GO_FLG[i] = 0;
 
                    else if (ctx->GO_FLG[i] < 0) {
                        if (first) {      
                            f = ctx->GO_LBF[i];
                            ja = i;
                            first = 0;    
                        }
                        else if (f > ctx->GO_LBF[i]) {
                            f = ctx->GO_LBF[i];
                            ja = i;
                        }
                    }
                }   
                else {

                    if (ctx->GO_UBF[i] < fminp)  
                        ctx->GO_FLG[i] = 0;
                    else if (ctx->GO_FLG[i] < 0) {
                        if (first) {      
                            f = ctx->GO_UBF[i];
                            ja = i;
                            first = 0;
                        }                 
                        else if (f < ctx->GO_UBF[i]) {
                            f = ctx->GO_UBF[i];
                            ja = i;
                        }   
                    }
                }   
            }
        }
        if (ja < 0)      /* here we have no more unprocessed boxes */
            break;
                    
        /* bisect box ja and enter subboxes into list */

        jb = -1;        /* find free box */

        for (i = 0; i < nbmax; ++i) {
            if (ctx->GO_FLG[i] == 0) {
                jb = i;
                break;
            }
        }
        if (jb < 0) {       /* exceeded max number of boxes */
            err = 1;
            break;          
        }
        if (jb >= nb) {
            nb = jb + 1;
            if (ctx->GO_NBU < nb)
                ctx->GO_NBU = nb;
        }
        lptra = ctx->GO_LB + ja * narg;
        uptra = ctx->GO_UB + ja * narg;

        lptrb = ctx->GO_LB + jb * narg;
        uptrb = ctx->GO_UB + jb * narg;
   
        is = 0;                         /* coordinate with largest width */
        w = uptra[0] - lptra[0];

        for (i = 1; i < narg; ++i) {
            tmp = uptra[i] - lptra[i];
            if (w < tmp) {
                w = tmp;
                is = i;
            }
        }

        if (w <= ctx->TOLBW) {               /* temporarily accept this box */
            ctx->GO_FLG[ja] = 1;
            continue;
        }

        for (i = 0; i < narg; ++i) {    /* copy box ja to jb */
            lptrb[i] = lptra[i];
            uptrb[i] = uptra[i];
        }

        tmp = lptra[is] + w / 2.0;
        uptra[is] = lptrb[is] = tmp;
        ctx->GO_FLG[jb] = -1;

        for (jj = 0; jj < 2; ++jj) {           /* for both boxes */

            if (jj)   
                ja = jb;

            lptra = ctx->GO_LB + ja * narg;
            uptra = ctx->GO_UB + ja * narg;

            /* evaluate inclusion function */

            ctx->GO_IFN++;
            if (typ == 0)
                r = ivreg_h(ctx, 1,narg,lptra,uptra,&fl,&fu,gc,ctx->AcU,ctx->AcV);                                      
            else
                r = ivreg1_h(ctx, mopt,1,narg,lptra,uptra,&fl,&fu,gc,ctx->AcU,ctx->AcV,fminp);                                      
            if (r) {
                printf1(ctx, "Error in evaluating inclusion function.\n");
                err = -2;
                goto IGMINFin;
            }
            if (opt == 0) {
                if (fl > fminp) {
                    ctx->GO_FLG[ja] = 0;             /* free box */
                    if (ja == nb - 1)
                        nb--;    
                    continue;
                }
            }
            else {
                if (fu < fminp) {
                    ctx->GO_FLG[ja] = 0;             /* free box */
                    if (ja == nb - 1)
                        nb--;    
                    continue;
                }
            }
            ctx->GO_LBF[ja] = fl;
            ctx->GO_UBF[ja] = fu;

            if (opt == 0) {
                if (fminp > fu) {
                    fminp = fu;
                    igmin_cpar(ctx, narg,uptra,par);
                }
            }
            else {
                if (fminp < fl) {
                    fminp = fl;
                    igmin_cpar(ctx, narg,lptra,par);
                }
            }
            if (fabs(fu - fl) < tolfd)      /* accept for solution */
                ctx->GO_FLG[ja] = 1;
   
            if (gc) {                       /* gradient check */

                r = 0;
                for (j = 0; j < narg; ++j) {
                    if (opt == 0) {
                        if (ctx->AcU[j] > 0.0) {
                            if (lptra[j] > lb[j]) {
                                r = 1;
                                break;
                            }
                            uptra[j] = lptra[j];
                            r = -1;
                        }
                        else if (ctx->AcV[j] < 0.0) {
                            if (uptra[j] < ub[j]) {
                                r = 1;
                                break;
                            }
                            lptra[j] = uptra[j];
                            r = -1;
                        }
                    }
                    else {
                        if (ctx->AcU[j] > 0.0) {
                            if (uptra[j] < ub[j]) {
                                r = 1;
                                break;
                            }
                            lptra[j] = uptra[j];
                            r = -1;
                        }
                        else if (ctx->AcV[j] < 0.0) {
                            if (lptra[j] > lb[j]) {
                                r = 1;
                                break;
                            }
                            uptra[j] = lptra[j];
                            r = -1;
                        }
                    }
                }
                if (r == 1) {
                    ctx->GO_FLG[ja] = 0;             /* free box */
                    if (ja == nb - 1)
                        nb--;    
                    continue;
                }
                else if (r == -1) {         /* update inclusion function */
                    ctx->GO_IFN++;
                    if (typ == 0)
                        r = ivreg_h(ctx, 1,narg,lptra,uptra,&fl,&fu,0,ctx->AcU,ctx->AcV);                                      
                    else
                        r = ivreg1_h(ctx, mopt,1,narg,lptra,uptra,&fl,&fu,0,ctx->AcU,ctx->AcV,fminp);                                      
                    if (r) {
                        printf1(ctx, "Error in evaluating inclusion function.\n");
                        err = -2;
                        goto IGMINFin;
                    }

                    if (opt == 0) {
                        if (fl > fminp) {
                            ctx->GO_FLG[ja] = 0;             /* free box */
                            if (ja == nb - 1)
                                nb--;    
                            continue;
                        }
                    }
                    else {
                        if (fu < fminp) {
                            ctx->GO_FLG[ja] = 0;             /* free box */
                            if (ja == nb - 1)
                                nb--;    
                            continue;
                        }
                    }
                    ctx->GO_LBF[ja] = fl;
                    ctx->GO_UBF[ja] = fu;
                }
            }

            /* update fminp with function value for midpoint of a box with
               least lower bound */
                   
            for (j = 0; j < narg; ++j)
                ctx->AcTmp[j] = (lptra[j] + uptra[j]) / 2.0;

            /* calculate function at midpoint */

            ctx->GO_FN++;
            if (typ == 0)
                r = ivreg_h(ctx, 0,narg,ctx->AcTmp,&tmp,&f,&tmp,0,ctx->AcU,ctx->AcV);                                      
            else
                r = ivreg1_h(ctx, mopt,0,narg,ctx->AcTmp,&tmp,&f,&tmp,0,ctx->AcU,ctx->AcV,fminp);                                      

            if (r) {
                printf1(ctx, "Error in evaluating function.\n");
                err = -2;
                goto IGMINFin;
            }
            if (opt == 0) {
                if (fminp > f) {
                    fminp = f;
                    igmin_cpar(ctx, narg,ctx->AcTmp,par);
                }
            }   
            else {
                if (fminp < f) {
                    fminp = f;
                    igmin_cpar(ctx, narg,ctx->AcTmp,par);
                }
            }   
        }
    }
    ctx->GO_NB = nb;
    ctx->GO_IT = iter;
    ctx->GO_FMIN = fminp;

    if (ctx->PMProtFDef) {
        fprintf(ctx->PMProtFd,"IGMIN optimization.\n");
        prval(ctx, "Best function value",ctx->GO_FMIN);
        prvec(ctx, "Parameters",narg,par - 1);
        fprintf(ctx->PMProtFd,"\n");
        
        if (opt == 0)
            tmp = ctx->GO_FMIN - tolfe;
        else
            tmp = ctx->GO_FMIN + tolfe;

        ja = 1;
        for (i = 0; i < nb; ++i) {
            if (ctx->GO_FLG[i] == 0)  
                continue;

            if (ctx->GO_FLG[i] >= 1) {

                if (opt == 0) {                 /* minimization */
                    if (ctx->GO_LBF[i] >= tmp)         
                        ctx->GO_FLG[i] = 2;
                }
                else {
                    if (ctx->GO_UBF[i] <= tmp)         
                        ctx->GO_FLG[i] = 2;
                }
            }
            lptra = ctx->GO_LB + i * narg;
            uptra = ctx->GO_UB + i * narg;

            w = uptra[0] - lptra[0];
            for (j = 1; j < narg; ++j)  
                w = dmax(ctx, w,uptra[j] - lptra[j]);

            fprintf(ctx->PMProtFd,"Box   Width                 Lower function bound  Upper function bound  Acceptance\n");
            fprintf(ctx->PMProtFd,"%3d  %21.14e %21.14e %21.14e  %6d\n",
                                 ja++,w,ctx->GO_LBF[i],ctx->GO_UBF[i],ctx->GO_FLG[i]);

            fprintf(ctx->PMProtFd,"Parameter vector\nLB: ");
            for (j = 0; j < narg; ++j)  
                rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,lptra[j]);
            fprintf(ctx->PMProtFd,"\nUB: ");
            for (j = 0; j < narg; ++j)  
                rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,uptra[j]);
            fprintf(ctx->PMProtFd,"\n\n");
        }
    }

IGMINFin:
    alloc_acu(ctx, 0);           
    alloc_acv(ctx, 0);           
    alloc_actmp(ctx, 0);           
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ivreg_h    objective and inclusion function for ivreg's exact step.     */
/*                                                                          */
/*  The quantity is the OLS slope                                           */
/*      b = (n*Sxy - Sx*Sy) / (n*Sxx - Sx*Sx)                               */
/*  where observations the fixing steps left undecided are replaced by      */
/*  the current parameters (IVRXP/IVRYP index them; -1 = observed value).   */
/*  opt=0 evaluates b at a point; opt=1 computes an interval enclosure of   */
/*  b over the parameter box via the i_* interval arithmetic; deriv adds    */
/*  bounds for the gradient in gl[i],gh[i].                                 */
/*                                                                          */
/*  ivreg_h1 further below is a preliminary substitute (its own words)      */
/*  using the mean-centred enclosure with precomputed variance bounds; it   */
/*  has NO call sites anywhere -- dead code, kept as Rohwer left it.        */
/*                                                                          */
/*  Return 0 if OK, -1 if error                                             */

int ivreg_h(TDAContext *ctx, int opt,int n,double *xl,double *xh, double *rl,double *rh,int deriv,double *gl,double *gh)
{
    register int i,k;   
    double nn,tmp,sxl,sxh,syl,syh,sxxl,sxxh,sxyl,sxyh;
    double xxl,xxh,yyl,yyh,sl,sh,tl,th,ul,uh,nl,nh;
       
/**    
tda_out("ivreg_h opt=%d\n",opt);

    for (i = 0; i < n; ++i) {
        printf1(ctx, "i=%d xl=%20.12f ",i,xl[i]);
        if (opt)
            printf1(ctx, "xh=%20.12f",xh[i]);
        newline(ctx);
    }     
**/
    nn = (double)ctx->NOC;     

    if (opt == 0) {
        sxl = syl = sxxl = sxyl = 0.0;
        for (i = 0; i < ctx->NOC; ++i) {
            if (ctx->IVRXP[i] >= 0)
                xxl = xl[ctx->IVRXP[i]];
            else
                xxl = ctx->IVRX[i];

            if (ctx->IVRYP[i] >= 0)
                yyl = xl[ctx->IVRYP[i]];
            else
                yyl = ctx->IVRY[i];

            sxl += xxl;
            syl += yyl;
            sxxl += xxl * xxl;
            sxyl += xxl * yyl;
        }
        tmp = nn * sxxl - sxl * sxl;
        if (tmp == 0.0) {
            printf1(ctx, "\nFATAL ERROR in ivreg_h: zero variance.\n");
            exit(0);
        }
        *rl = (nn * sxyl - sxl * syl) / tmp;
        return(0);
    }
    sxl = syl = sxxl = sxyl = 0.0;
    sxh = syh = sxxh = sxyh = 0.0;
    *rl = *rh = 0.0;

    for (i = 0; i < ctx->NOC; ++i) {

        if (ctx->IVRXP[i] >= 0) {
            xxl = xl[ctx->IVRXP[i]];
            xxh = xh[ctx->IVRXP[i]];
        }
        else
            xxl = xxh = ctx->IVRX[i];

        if (ctx->IVRYP[i] >= 0) {
            yyl = xl[ctx->IVRYP[i]];
            yyh = xh[ctx->IVRYP[i]];
        }
        else
            yyl = yyh = ctx->IVRY[i];

        i_add(ctx, sxl,sxh,xxl,xxh,&tl,&th);
        sxl = tl;
        sxh = th;
        i_add(ctx, syl,syh,yyl,yyh,&tl,&th);
        syl = tl;
        syh = th;
        i_square(ctx, xxl,xxh,&sl,&sh);
        i_add(ctx, sxxl,sxxh,sl,sh,&tl,&th);
        sxxl = tl;
        sxxh = th;
        i_mul(ctx, yyl,yyh,xxl,xxh,&sl,&sh);
        i_add(ctx, sxyl,sxyh,sl,sh,&tl,&th);
        sxyl = tl;
        sxyh = th;
    }
    i_square(ctx, sxl,sxh,&sl,&sh);
    i_sub(ctx, nn * sxxl,nn * sxxh,sl,sh,&nl,&nh);
    i_mul(ctx, sxl,sxh,syl,syh,&sl,&sh);
    i_sub(ctx, nn * sxyl,nn * sxyh,sl,sh,&ul,&uh);

    /* variance must be positive */

    if (nl < ctx->IVRVARL) {
        nl = ctx->IVRVARL;
        if (nh < nl)
            nh = nl;
    } 
    i_div(ctx, ul,uh,nl,nh,rl,rh);          
/*       
tda_out(">>>>>>>>>> rl=%g rh=%g\n",*rl,*rh);
*/
    if (deriv) {

        if (nl <= 0.0 && nh >= 0.0) {
            for (i = 0; i < n; ++i) {
                gl[i] = -1.0;
                gh[i] =  1.0;
            }
        }
        else {
            for (i = 0; i < ctx->NOC; ++i) {

                if ((k = ctx->IVRYP[i]) >= 0) {   

                    if (ctx->IVRXP[i] >= 0) {
                        xxl = xl[ctx->IVRXP[i]];
                        xxh = xh[ctx->IVRXP[i]];
                    }
                    else
                        xxl = xxh = ctx->IVRX[i];

                    i_sub(ctx, nn * xxl,nn * xxh,sxl,sxh,&tl,&th);
                    gl[k] = tl;
                    gh[k] = th;

/*   tda_out("gradY k=%d gl=%g gh=%f\n",k,gl[k],gh[k]);    */


                }
                if ((k = ctx->IVRXP[i]) >= 0) {   

                    xxl = xl[k];
                    xxh = xh[k];
                       
                    if (ctx->IVRYP[i] >= 0) {
                        yyl = xl[ctx->IVRYP[i]];
                        yyh = xh[ctx->IVRYP[i]];
                    }
                    else
                        yyl = yyh = ctx->IVRY[i];

                    i_sub(ctx, nn * xxl,nn * xxh,sxl,sxh,&sl,&sh);
                    i_mul(ctx, 2.0 * *rl,2.0 * *rh,sl,sh,&tl,&th);
                    i_sub(ctx, nn * yyl,nn * yyh,syl,syh,&sl,&sh);
                    i_sub(ctx, sl,sh,tl,th,&ul,&uh);
                    gl[k] = ul;
                    gh[k] = uh;

/*  tda_out("gradX k=%d gl=%g gh=%f\n",k,gl[k],gh[k]);     */
                }
            }
        }
    }    
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  ivreg_h1(opt,n,xl,xh,rl,rh,deriv,gl,gh)                                 */
/*                                                                          */
/*  Preliminary substitute for ivreg_h.                                     */
/*                                                                          */
/*  Calculate function depending on opt.                                    */
/*                                                                          */
/*  If opt = 0 : standard function with arguments xl[i], i=0,...,n-1        */
/*               return function value in rl.                               */
/*  If opt = 1 : inclusion function at [xl[i],xh[i]], i = 0,...,n-1         */
/*               return interval in [rl,rh].                                */
/*                                                                          */
/*  If deriv != 0 and opt=1, calculate bounds for gradient in gl[i],gh[i].  */
/*                                                                          */
/*  Return 0 if OK, -1 if error                                             */

int ivreg_h1(TDAContext *ctx, int opt,int n,double *xl,double *xh, double *rl,double *rh,int deriv,double *gl,double *gh)
{
    (void)n;        /* unused: the signature is shared */
    register int i,k;   
    double nn,tmp,sxl,syl,sxxl,sxyl;
    double xxl,xxh,yyl,yyh,sl,sh,tl,th,ul,uh,nl,nh;
    double xml,xmh,yml,ymh;
       
/**    
tda_out("ivreg_h opt=%d\n",opt);

    for (i = 0; i < n; ++i) {
        printf1(ctx, "i=%d xl=%20.12f ",i,xl[i]);
        if (opt)
            printf1(ctx, "xh=%20.12f",xh[i]);
        newline(ctx);
    }     
**/
    nn = (double)ctx->NOC;     

    if (opt == 0) {
        sxl = syl = sxxl = sxyl = 0.0;
        for (i = 0; i < ctx->NOC; ++i) {
            if (ctx->IVRXP[i] >= 0)
                xxl = xl[ctx->IVRXP[i]];
            else
                xxl = ctx->IVRX[i];

            if (ctx->IVRYP[i] >= 0)
                yyl = xl[ctx->IVRYP[i]];
            else
                yyl = ctx->IVRY[i];

            sxl += xxl;
            syl += yyl;
            sxxl += xxl * xxl;
            sxyl += xxl * yyl;
        }
        tmp = nn * sxxl - sxl * sxl;
        if (tmp == 0.0) {
            printf1(ctx, "\nFATAL ERROR in ivreg_h: zero variance.\n");
            exit(0);
        }
        *rl = (nn * sxyl - sxl * syl) / tmp;
        return(0);
    }

    /* opt = 1 */

    xml = xmh = yml = ymh = 0.0;
    for (i = 0; i < ctx->NOC; ++i) {

        if (ctx->IVRXP[i] >= 0) {
            xxl = xl[ctx->IVRXP[i]];
            xxh = xh[ctx->IVRXP[i]];
        }
        else
            xxl = xxh = ctx->IVRX[i];

        if (ctx->IVRYP[i] >= 0) {
            yyl = xl[ctx->IVRYP[i]];
            yyh = xh[ctx->IVRYP[i]];
        }
        else
            yyl = yyh = ctx->IVRY[i];
      
        xml += xxl;
        xmh += xxh;
        yml += yyl;
        ymh += yyh;
    }

    *rl = *rh = 0.0;
    ul = uh = 0.0;

    for (i = 0; i < ctx->NOC; ++i) {

        if (ctx->IVRXP[i] >= 0) {
            xxl = xl[ctx->IVRXP[i]];
            xxh = xh[ctx->IVRXP[i]];
        }
        else
            xxl = xxh = ctx->IVRX[i];

        if (ctx->IVRYP[i] >= 0) {
            yyl = xl[ctx->IVRYP[i]];
            yyh = xh[ctx->IVRYP[i]];
        }
        else
            yyl = yyh = ctx->IVRY[i];
      
        i_sub(ctx, nn * yyl,nn * yyh,yml,ymh,&sl,&sh);
        i_mul(ctx, sl,sh,xxl,xxh,&tl,&th);
        i_add(ctx, ul,uh,tl,th,&sl,&sh);
        ul = sl;
        uh = sh;
    }
    if (i_div(ctx, ul,uh,nn * nn * ctx->IVRVARL,nn * nn * ctx->IVRVARU,rl,rh)) {
        printf1(ctx, "FATAL ERROR in ivreg_k.\n");
        exit(0);
    }

    if (deriv) {
        xml /= nn;
        xmh /= nn;
        yml /= nn;
        ymh /= nn;

        for (i = 0; i < ctx->NOC; ++i) {

            if (ctx->IVRXP[i] >= 0) {
                xxl = xl[ctx->IVRXP[i]];
                xxh = xh[ctx->IVRXP[i]];
            }
            else
                xxl = xxh = ctx->IVRX[i];

            if (ctx->IVRYP[i] >= 0) {
                yyl = xl[ctx->IVRYP[i]];
                yyh = xh[ctx->IVRYP[i]];
            }
            else
                yyl = yyh = ctx->IVRY[i];
          
            if ((k = ctx->IVRYP[i]) >= 0) {   

                i_sub(ctx, xxl,xxh,xml,xmh,&tl,&th);
                gl[k] = tl;
                gh[k] = th;
            }
            if ((k = ctx->IVRXP[i]) >= 0) {

                i_sub(ctx, yyl,yyh,yml,ymh,&sl,&sh);
                i_sub(ctx, xxl,xxh,xml,xmh,&tl,&th);
                i_mul(ctx, 2.0 * *rl,2.0 * *rh,tl,th,&nl,&nh);
                i_sub(ctx, sl,sh,nl,nh,&tl,&th);
                if (tl > 0.0)
                    gl[k] = gh[k] = 1.0;
                else if (th < 0.0)
                    gl[k] = gh[k] = -1.0;
                else {
                    gl[k] = -1.0;
                    gh[k] =  1.0;
                }
            }
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  ivar1f(n,xl,xh,ip,idx,v0,v1)                                            */
/*                                                                          */
/*  Calculate bounds for variance of (xl,xh). Return bounds in vmin,vmax.   */
/*  Return: 0 if both bounds can be calculated,                             */
/*          1 if only lower bound                                           */ 
/*          2 if only upper bound                                           */
/*         -1 if both cannot be calculated                                  */

int ivar1f(TDAContext *ctx, int n,double *xl,double *xh,double *ip,int *idx,double *v0,double *v1)
{
    register int i,j,np;
    int nn,fin,r;   
    double tmp,mmin,mmax,vmin,vmax,a,b,u,v,mmin1,mmax1;
         
    r = -1;
    np = 0;
    mmin = mmax = 0.0;

    for (i = 0; i < n; ++i) {
        ip[np++] = xl[i];
        ip[np++] = xh[i];
        mmin += xl[i];
        mmax += xh[i];
    }
    mmin /= (double)n;
    mmax /= (double)n;

    if (sortd(ctx, np,ip,0))
        return(-1);
   
    nn = 1;
    for (i = 1; i < np; ++i) {          /* ip contains induced partition */
        if (ip[i] > ip[i - 1])
            ip[nn++] = ip[i];
    }
    vmin = vmax = 0.0;
    for (i = 0; i < n; ++i) {
        u = mmin - xl[i];
        vmin += u * u;
        u = mmax - xh[i];
        vmax += u * u;
    }
    vmin /= (double)n;
    vmax /= (double)n;

    vmin = ctx->DBLMAX;
    fin = 0;
    for (i = 1; i < nn; ++i) {
        a = ip[i - 1];
        b = ip[i];
        if (b < mmin)
            continue;
        if (a > mmax)
            break;
            
        np = 0;
        u = 0.0;
        for (j = 0; j < n; ++j) {
            if (xh[j] <= a) {
                u += xh[j];
                np++;
            }
            else if (xl[j] >= b) {
                u += xl[j];
                np++;
            }
        }
        if (np > 0) {
            u /= (double)np;  
            tmp = ivarf(ctx, n,u,xl,xh);
            tmp /= (double)n;
            if (tmp < vmin) {
                vmin = tmp;
                fin = 1;
            }
        }
    }
    if (fin) {
        *v0 = vmin;
        r = 1;
    }
    for (i = 0; i < n; ++i)
        idx[i] = 0;

    mmin1 = mmin;
    mmax1 = mmax;
    fin = 0;
    for (i = 1; i < 10; ++i) {
        np = 0;
        mmin = mmax = 0.0;

        for (j = 0; j < n; ++j) {
            if (idx[j] == 0) {
                np++;

                tmp = (xl[j] + xh[j]) / 2.0;
                if (tmp <= mmin1) {
                    idx[j] = 1;
                    mmin += xl[j];
                    mmax += xl[j];
                }
                else if (tmp >= mmax1) {
                    xl[j] = xh[j];
                    idx[j] = 1;
                    mmin += xh[j];
                    mmax += xh[j];
                }
                else {
                    mmin += xl[j];
                    mmax += xh[j];
                }           
            }
            else {
                mmin += xl[j];
                mmax += xl[j];
            }
        }
        mmin /= (double)n;
        mmax /= (double)n;
        if (np == 0) {
            fin = 1;
            break;
        }
        mmin1 = mmin;
        mmax1 = mmax;
    }
    if (fin) {
        u = 0.0;
        for (j = 0; j < n; ++j)
            u += xl[j];
        u /= (double)n;

        v = 0.0;
        for (j = 0; j < n; ++j)  
            v += (xl[j] - u) * (xl[j] - u);
        v /= (double)n;

        *v1 = v;
        if (r == 1)
            r = 0;
        else
            r = 2;
    }
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  ivreg_z(np,rp,ns,opt)  Heuristic method to minimize (opt=0) or          */
/*                         maximize (opt = 1) IVRBETA. rp is number of      */
/*                         remaining coordinates, ns is block size.         */
/*                                                                          */  
/*  Return 0 if OK, -1 if error (insufficient memory)                       */

int ivreg_z(TDAContext *ctx, int np,int rp,int ns,int opt)
{
    register int i,j,k;
    int nx,first,iter,bf;
    double b;

    tda_err("\nHeuristic optimization.\n");
    printf1(ctx, "\nHeuristic optimization.\n");
    printf1(ctx, "Block  Size              Beta\n");

    iter = 0;

    while (1) {     
        iter++;
        if (alloc_aci(ctx, rp + 1))   
            return(-1);          

        /* AcI[.] encodes which RPar index is still open: -i for the
           even (center) slot, i+1 for the odd (radius) slot -- shifted
           by one so it is never 0, which used to collide with the even
           slot's own -i encoding at i=0 (both stored as plain 0, so
           the decode below always read it back as the even slot,
           silently never fixing the odd one at i=0; the outer while(1)
           loop then never saw k reach 0, looping forever whenever that
           was the last parameter left open). Decode below matches. */
        k = 0;
        for (i = 0; i < np; i += 2) {
            if (fabs(ctx->RParU[i] - ctx->RParL[i]) > ctx->TOLBW)
                ctx->AcI[k++] = -i;
            if (fabs(ctx->RParU[i + 1] - ctx->RParL[i + 1]) > ctx->TOLBW)
                ctx->AcI[k++] = i + 1;
        }
        if (k == 0)
            break;
        nx = imin(ctx, k,ns);

        if (alloc_acj(ctx, nx + 1))  
            return(-1);     
        if (alloc_acl(ctx, nx + 1))  
            return(-1);     

        if (alloc_acu(ctx, nx + 1))  
            return(-1);           
        if (alloc_acv(ctx, nx + 1))  
            return(-1);           

        for (i = 0; i < nx; ++i) {
            j = ctx->AcI[i];
            if (j <= 0) {
                ctx->AcU[i] = ctx->RParL[-j];
                ctx->AcV[i] = ctx->RParU[-j];
            }
            else {
                ctx->AcU[i] = ctx->RParL[j];
                ctx->AcV[i] = ctx->RParU[j];
            }
        }
        first = 1;
        while (comb_nm1(ctx, 2,nx,ctx->AcJ,first)) {
            first = 0;
            for (i = 0; i < nx; ++i) {
                j = ctx->AcI[i];
                if (j <= 0) {
                    if (ctx->AcJ[i] == 0)
                        ctx->RPar[-j] = ctx->RParL[-j] = ctx->RParU[-j] = ctx->AcU[i];
                    else
                        ctx->RPar[-j] = ctx->RParL[-j] = ctx->RParU[-j] = ctx->AcV[i];
                }
                else {
                    if (ctx->AcJ[i] == 0)
                        ctx->RPar[j] = ctx->RParL[j] = ctx->RParU[j] = ctx->AcU[i];
                    else
                        ctx->RPar[j] = ctx->RParL[j] = ctx->RParU[j] = ctx->AcV[i];
                }
            }
            b = ivreg_b0(ctx, np,ctx->RPar);   
            bf = 0;
            if (opt == 0 && b < ctx->IVRBETA) {  
                ctx->IVRBETA = b;
                bf = 1;
            }
            else if (opt != 0 && b > ctx->IVRBETA) {  
                ctx->IVRBETA = b;
                bf = 1;
            }
            if (bf) {
                tda_err("Block %4d  Beta %20.8f\n",iter,ctx->IVRBETA);
                for (i = 0; i < nx; ++i)
                    ctx->AcL[i] = ctx->AcJ[i];
            }
        }

        /* fix best values */

        for (i = 0; i < nx; ++i) {
            j = ctx->AcI[i];
            if (j <= 0) {  
                if (ctx->AcL[i] == 0)
                    ctx->RPar[-j] = ctx->RParL[-j] = ctx->RParU[-j] = ctx->AcU[i];
                else
                    ctx->RPar[-j] = ctx->RParL[-j] = ctx->RParU[-j] = ctx->AcV[i];
            }
            else {  
                if (ctx->AcL[i] == 0)
                    ctx->RPar[j] = ctx->RParL[j] = ctx->RParU[j] = ctx->AcU[i];
                else
                    ctx->RPar[j] = ctx->RParL[j] = ctx->RParU[j] = ctx->AcV[i];
            }
        }
        printf1(ctx, "%5d  %4d  %16.8f\n",iter,nx,ctx->IVRBETA);
#ifdef TDA_R_PACKAGE
        /* iteration, changed coordinates, current best slope value --
           the heuristic reports achieved points, not bounds */
        tda_export_cell(ctx, "iv.heuristic", (double)iter);
        tda_export_cell(ctx, "iv.heuristic", (double)nx);
        tda_export_cell(ctx, "iv.heuristic", ctx->IVRBETA);
        tda_export_endrow(ctx, "iv.heuristic");
#endif
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  ivreg1      Regression with interval data.                              */
/*                                                                          */
/*              ivreg1(                                                     */
/*                  nbox=...,       max number of boxes, def. 100           */
/*                  mxit=...,       max number of iterations, def. 100      */
/*                  tolbw=...,      tolerance box widght, def. 1.e-4        */
/*                  tolf=...,       tolerance for beta iteration, 1.e-6     */
/*                  tolfd=...,      tolerance function range, def. 1.e-10   */
/*                  tolbw=...,      tolerance global min/max, def. 1.e-10   */
/*                  fmt=...,        print format, def. 13.6                 */  
/*                  df=...,         output file containing data             */
/*                  prot=...,       protocol file                           */
/*              ) = YL,YH,XL,XH;                                            */
/*                                                                          */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int ivreg1(TDAContext *ctx)
{
    register int i,j,k;
    int err,n,npp,r;
    double xl,xh;
           
    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    ctx->TOLF  = 1.e-6;       /* def tolerance for beta iteration                 */
    ctx->TOLBW = 1.e-4;       /* def tolerance for box length                     */
    ctx->TOLFD = 1.e-10;      /* def tolerance for function range                 */
    ctx->TOLFE = 1.e-10;      /* def tolerance for global minimum                 */

    printf1(ctx, "Regression with interval data. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 6,4,1))     /* get parameters */
        goto IVREGFin;

    if (ctx->PMNV != 4) {   
        printf1(ctx, "Error: need four variables on right-hand side.\n");
        goto IVREGFin;
    }
    if (ctx->PMFmtF == 0)        /* default print format */
        pmfmt(ctx, 13,6);

    n = ctx->NOC;
    npp = 4;

    printf1(ctx, "\nNumber of data points: %d\n",n);
    if (n < 2) {
        printf1(ctx, "Error: need at least two data points.\n");
        goto IVREGFin;
    }

    if (alloc_acx(ctx, 2 * ctx->NOC + 1))     
        goto IVREGFin;
    if (alloc_acy(ctx, 2 * ctx->NOC + 1))     
        goto IVREGFin;

    k = 0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < 2; ++j) {
            xl = get_data(ctx, ctx->PMVIdx[j * 2],i);
            xh = get_data(ctx, ctx->PMVIdx[j * 2 + 1],i);
            if (xh < xl) {
                printf1(ctx, "Error: no valid interval in record %d.\n",i + 1);
                goto IVREGFin;
            }
            ctx->AcX[k] = xl;
            ctx->AcY[k] = xh;
            k++;
        }
    }
    /**
    tda_out("data\n"); 
    for (i = 0; i < np; ++i) {
        tda_out("%3d %f %f\n",i,AcX[i],AcY[i]);
    }
    **/ 
    /* Same scaled defaults as ivreg's opt=3, for the same reason: the
       flat 100/100 made even the shipped example give up uncertified.
       User-set values are respected exactly as before. */

    if (ctx->MxIter < 1 || ctx->MxItFlg == 0) {
        ctx->MxIter = 20000;
        for (i = 7; i < npp && ctx->MxIter < 1000000; ++i)
            ctx->MxIter *= 2;
        if (ctx->MxIter > 1000000)
            ctx->MxIter = 1000000;
    }
    if (ctx->PMNBOX < 1) {       /* max number of boxes */
        ctx->PMNBOX = ctx->MxIter / 10;
        if (ctx->PMNBOX < 1000)
            ctx->PMNBOX = 1000;
        if (ctx->PMNBOX > 100000)
            ctx->PMNBOX = 100000;
    }

    if (alloc_par(ctx, npp,1))
        goto IVREGFin;

    if (alloc_list(ctx, ctx->PMNBOX,npp)) 
        goto IVREGFin;

    printf1(ctx, "Number of parameters: %d\n",npp);              
    printf1(ctx, "Maximum number of iterations: %d\n",ctx->MxIter);              
    printf1(ctx, "Maximum number of boxes: %d\n",ctx->PMNBOX);
    printf1(ctx, "Tolerance for box width: %g\n",ctx->TOLBW);
    printf1(ctx, "Tolerance for function range: %g\n",ctx->TOLFD);
    printf1(ctx, "Tolerance for global minimum: %g\n\n",ctx->TOLFE);
      
    ctx->RParL[0] =  1.0;
    ctx->RParU[0] =  2.5;
    ctx->RParL[1] =   0.0;
    ctx->RParU[1] =   0.5;
    ctx->RParL[2] =   0.0;
    ctx->RParU[2] =   0.1;
    ctx->RParL[3] =   0.0;
    ctx->RParU[3] =   0.1;
    for (i = 0; i < npp; ++i)
        ctx->RPar[i] = (ctx->RParL[i] + ctx->RParU[i]) / 2.0;
      
    ctx->RPar[0] = 2.0889;      
    ctx->RPar[1] = 0.368047;
    ctx->RPar[2] = 0.0797882;
    ctx->RPar[3] = 0.05;        
#ifdef TDA_R_PACKAGE
    /* sbox= replaces the hard-coded search box above (see t_parm.c);
       the start point becomes the box's midpoint. */
    if (ctx->IVSBoxSet) {
        for (i = 0; i < 4 && i < npp; ++i) {
            ctx->RParL[i] = ctx->IVSBox[i];
            ctx->RParU[i] = ctx->IVSBox[4 + i];
            ctx->RPar[i] = (ctx->RParL[i] + ctx->RParU[i]) / 2.0;
        }
        ctx->IVSBoxSet = 0;
    }
#endif
             


    r = ivreg_min(ctx, 1,0,0,npp,ctx->PMNBOX,ctx->RPar,ctx->RParL,ctx->RParU,ctx->MxIter,ctx->TOLBW,ctx->TOLFD,ctx->TOLFE,0);
    if (r) {
        if (r == 1)
            printf1(ctx, "Exceeded maximum number of boxes.\n");
        printf1(ctx, "The value below is the best point found, not proven\n");
        printf1(ctx, "optimal; raise nbox= (and mxit=) to certify it.\n");
        printf1(ctx, "Current best value of beta:");
        rt_printf1_d(ctx, ctx->PMFmtS,ctx->GO_FMIN);
        newline(ctx);
        goto IVREGFin;
    }
    igmin_res(ctx, ctx->TOLFE,0);
    /***
    for (i = 0;i < npp; ++i)
        tda_out("i=%d par=%f\n",i,RPar[i]);
    ***/

          
    newline(ctx);
    err = 0;

IVREGFin:
    alloc_list(ctx, 0,0);
    alloc_par(ctx, 0,0);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ivreg1_h(mopt,opt,n,pl,ph,rl,rh,deriv,gl,gh,fminp)                      */
/*                                                                          */
/*  Calculate function depending on opt.                                    */
/*  If mopt == 0 absolute, else squared.                                    */
/*                                                                          */
/*  If opt = 0 : standard function with arguments pl[i], i=0,...,n-1        */
/*               return function value in rl.                               */
/*  If opt = 1 : inclusion function at [pl[i],ph[i]], i = 0,...,n-1         */
/*               return interval in [rl,rh].                                */
/*                                                                          */
/*  If deriv != 0 and opt=1, calculate bounds for gradient in gl[i],gh[i].  */
/*                                                                          */
/*  Return 0 if OK, -1 if error                                             */

int ivreg1_h(TDAContext *ctx, int mopt,int opt,int n,double *pl,double *ph, double *rl,double *rh,int deriv,double *gl,double *gh,double fminp)
{
    (void)deriv; (void)fminp; (void)gh; (void)gl; (void)mopt; (void)n;        /* unused: the signature is shared */
    register int i;   
    int np;
    double xc,xr,yc,yr,d1,d2;
    double acl,ach,arl,arh,bcl,bch,brl,brh,sl,sh,tl,th,ul,uh,vl,vh;
    double abl,abh,ax,d1l,d1h,d2l,d2h,ml,mh,ml1,mh1,ml2,mh2;

/**     
tda_out("ivreg1_h opt=%d mopt=%d\n",opt,mopt);

    for (i = 0; i < n; ++i) {
        printf1(ctx, "i=%d xl=%20.12f ",i,pl[i]);
        if (opt)
            printf1(ctx, "xh=%20.12f",ph[i]);
        newline(ctx);
    }     
**/       

    np = 2 * ctx->NOC;
        
    *rl = *rh = 0.0;

    if (opt == 0) {
   
        for (i = 0; i < np; i += 2) {
            yc = (ctx->AcX[i] + ctx->AcY[i]) / 2.0;
            yr = (ctx->AcY[i] - ctx->AcX[i]) / 2.0;
            xc = (ctx->AcX[i + 1] + ctx->AcY[i + 1]) / 2.0;
            xr = (ctx->AcY[i + 1] - ctx->AcX[i + 1]) / 2.0;

            d1 = fabs(yc - pl[0] - i_center(ctx, xc,xr,pl[2],pl[3]));
            d2 = fabs(yr - pl[1] - i_radius(ctx, xc,xr,pl[2],pl[3]));
            *rl += d1 * d1 + d2 * d2;
        }
/**
        tda_out("*rl=%g\n",*rl);          
exit(0);    
**/

        return(0);
    }

    acl = pl[0];
    ach = ph[0];
    arl = pl[1];
    arh = ph[1];
    bcl = pl[2];
    bch = ph[2];
    brl = pl[3];
    brh = ph[3];

    i_sub(ctx, bcl,bch,brl,brh,&sl,&sh);
    i_abs(ctx, sl,sh,&tl,&th);
    i_add(ctx, bcl,bch,brl,brh,&sl,&sh);
    i_abs(ctx, sl,sh,&ul,&uh);
    i_max(ctx, tl,th,ul,uh,&abl,&abh);



    for (i = 0; i < np; i += 2) {

        yc = (ctx->AcX[i] + ctx->AcY[i]) / 2.0;
        yr = (ctx->AcY[i] - ctx->AcX[i]) / 2.0;
        xc = (ctx->AcX[i + 1] + ctx->AcY[i + 1]) / 2.0;
        xr = (ctx->AcY[i + 1] - ctx->AcX[i + 1]) / 2.0;

        ax = dmax(ctx, fabs(xc - xr),fabs(xc + xr));

        i_mul(ctx, xr,xr,abl,abh,&ml1,&mh1);
        i_mul(ctx, brl,brh,ax,ax,&ml2,&mh2);

        i_max(ctx, ml1,mh1,ml2,mh2,&ml,&mh);

        i_abs(ctx, bcl,bch,&ul,&uh);
        i_mul(ctx, xr,xr,ul,uh,&sl,&sh);

        i_mul(ctx, brl,brh,fabs(xc),fabs(xc),&tl,&th);
        i_add(ctx, sl,sh,tl,th,&ml2,&mh2);
        i_max(ctx, ml,mh,ml2,mh2,&ml1,&mh1);

        i_add(ctx, arl,arh,ml1,mh1,&sl,&sh);
        i_sub(ctx, yr,yr,sl,sh,&ul,&uh);
        i_abs(ctx, ul,uh,&d2l,&d2h);
/*
        tda_out("d2l=%g %g\n",d2l,d2h);
*/
        i_mul(ctx, xc,xc,bcl,bch,&sl,&sh);

        i_abs(ctx, bcl,bch,&tl,&th);
        i_mul(ctx, xr,xr,tl,th,&ml1,&mh1);
        i_mul(ctx, brl,brh,fabs(xc),fabs(xc),&ml2,&mh2);
        i_min(ctx, ml1,mh1,ml2,mh2,&ml,&mh);
        i_mul(ctx, xr,xr,brl,brh,&ml2,&mh2);
        i_min(ctx, ml,mh,ml2,mh2,&ml1,&mh1);

        if (sl > 0.0) {
            i_add(ctx, sl,sh,ml1,mh1,&ul,&uh);
        }
        else if (sh < 0.0) {
            i_sub(ctx, sl,sh,ml1,mh1,&ul,&uh);
        }
        else {
            i_add(ctx, sl,sh,ml1,mh1,&tl,&th);
            i_sub(ctx, sl,sh,ml1,mh1,&vl,&vh);
            i_max(ctx, tl,th,vl,vh,&sl,&uh);
            i_min(ctx, tl,th,vl,vh,&ul,&sh);

/**
tda_out("tl=%g %g vl=%g %g sl=%g uh=%g ul=%g sh=%g\n", tl,th,vl,vh,sl,uh,ul,sh);
exit(0);
**/

        }
        i_add(ctx, acl,ach,ul,uh,&sl,&sh);
        i_sub(ctx, yc,yc,sl,sh,&tl,&th);
        i_abs(ctx, tl,th,&d1l,&d1h);
/*
        tda_out("d1l=%g %g\n",d1l,d1h);
*/

        i_square(ctx, d1l,d1h,&sl,&sh);
        i_square(ctx, d2l,d2h,&tl,&th);
        i_add(ctx, sl,sh,tl,th,&ul,&uh);
        i_add(ctx, *rl,*rh,ul,uh,&tl,&th);
        *rl = tl;
        *rh = th;

    }
/*
    tda_out("*rl=%g %g\n",*rl,*rh);             
**/         
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  i_center(ac,ar,bc,br)   Return center of a * b                          */

double i_center(TDAContext *ctx, double ac,double ar,double bc,double br)
{
    double c,d;             

    c = ac * bc;
    d = c + dsign(ctx, c) * dmin(ctx, ar * fabs(bc),dmin(ctx, br * fabs(ac),ar * br));
    return(d);
}

/* ------------------------------------------------------------------------ */
/*  i_radius(ac,ar,bc,br)   Return radius of a * b                          */

double i_radius(TDAContext *ctx, double ac,double ar,double bc,double br)
{
    double aa,ba,d;             

    aa = dmax(ctx, fabs(ac + ar),fabs(ac - ar));
    ba = dmax(ctx, fabs(bc + br),fabs(bc - br));
    d = dmax(ctx, ar * ba,dmax(ctx, br * aa,ar * fabs(bc) + br * fabs(ac)));
    return(d);
}

/* ------------------------------------------------------------------------ */
/*  i_bvar                                                                  */

void i_bvar(TDAContext *ctx, double xl,double xh,double bl,double bh,double *ul,double *uh)
{

    if (xl >= 0.0) {
        if (bl >= 0.0) {
            *ul = xl * bl;
            *uh = xh * bh;
        }
        else if (bh <= 0.0) {
            *ul = xh * bl;
            *uh = xl * bh;
        }
        else {
            *ul = xh * bl;
            *uh = xh * bh;
        }
    }
    else if (xh <= 0.0) {
        if (bl >= 0.0) {
            *ul = xl * bh;
            *uh = xh * bl;
        }
        else if (bh <= 0.0) {
            *ul = xh * bh;
            *uh = xl * bl;
        }
        else {
            *ul = xl * bh;
            *uh = xl * bl;
        }
    }
    else {
        if (bl >= 0.0) {
            *ul = xl * bh;
            *uh = xh * bh;
        }
        else if (bh <= 0.0) {
            *ul = xh * bl;
            *uh = xl * bl;
        }
        else {
            *ul = dmin(ctx, xl * bh,xh * bl);
            *uh = dmin(ctx, xl * bl,xh * bh);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  i_norm.   Return distance.                                              */

double i_norm(TDAContext *ctx, int opt,double xl,double xh,double yl,double yh)
{
    double d;             

    d = dmax(ctx, fabs(xl - yl),fabs(xh - yh));
    if (opt)
        d *= d;
    return(d);
}

/* ------------------------------------------------------------------------ */
/*  ivreg1_nrm.    Return distance.                                          */

double ivreg1_nrm(TDAContext *ctx, int opt,double xl,double xh,double yl,double yh, double al,double ah,double bl,double bh)
{
    double sl,sh,tl,th;

    i_mul(ctx, xl,xh,bl,bh,&sl,&sh);
    i_add(ctx, al,ah,sl,sh,&tl,&th);
    sl = dmax(ctx, fabs(yl - tl),fabs(yh - th));
    if (opt)
        sl *= sl;
    return(sl);
}

/* -##--------------------------------------------------------------------- */
/*  ivreg2      Regression with interval data.                              */
/*                                                                          */
/*              ivreg2(                                                     */
/*                  opt=...,        option, def. 1                          */
/*                                  1 = direct calculation                  */
/*                                  2 = direct search                       */
/*                                  3 = mina                                */
/*                                  4 = contour plot                        */
/*                                                                          */
/*                  slen=...        initial step length, def. 1             */
/*                  sred=...,       reduction factor, def. 0.5              */
/*                  fmt=...,        print format, def. 14.8                 */  
/*                                                                          */
/*                  n=              number of level (automatically created) */
/*                  a=bc,br,        required for n option                   */
/*                  x=...           explicitly defined level                */
/*                  nn=...,         grid specification, def. 10,10          */
/*                  lt=...,         line type, def. 1 (solid line)          */
/*                  lw=...,         line width, def. 0.2 (mm)               */  
/*                  gs=...,         grey scale value, def. 1 (white)        */
/*              ) = YL,YH,XL,XH;                                            */
/*                                                                          */
/*  There are two possibilities to create a contour plot. Both require a    */
/*  plot file set up.                                                       */  
/*  a)  with n and a parameters, it is created automatically.               */
/*  b)  function levels specified with the x parameter.                     */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int ivreg2(TDAContext *ctx)
{
    register int i;
    int err,r,rn,n,nw;
    double x,y,lsn,bc,br,ac,ar,f,d1,d2;
           
    err = -1;         
    if (check_cmd(ctx, 2))
        return(-1);

    set_mldef(ctx);        /* set defaults for minimization */
    ctx->SLen = 1.0;         /* default step length */
    ctx->SRed = 0.5;         /* default reduction factor */
    ctx->TOLS = 1.0e-8;      /* minimum of step length */

    printf1(ctx, "Regression with interval data. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 6,4,1))     /* get parameters */
        goto IVREGFin;

    if (ctx->PMOPT > 4)
        ctx->PMOPT = 4;

    if (ctx->PMNV != 4) {   
        printf1(ctx, "Error: need four variables on right-hand side.\n");
        goto IVREGFin;
    }
    if (ctx->PMFmtF == 0)        /* default print format */
        pmfmt(ctx, 14,8);

    printf1(ctx, "\nNumber of data intervals: %d\n",ctx->NOC);
    if (ctx->NOC < 2) {
        printf1(ctx, "Error: need at least two data intervals.\n");
        goto IVREGFin;
    }

    /* SET UP DATA ==================================================== */

    if (alloc_acx(ctx, ctx->NOC + 1))             /* X center */
        goto IVREGFin;
    if (alloc_acy(ctx, ctx->NOC + 1))             /* Y center */
        goto IVREGFin;
    if (alloc_acu(ctx, ctx->NOC + 1))             /* X radius */
        goto IVREGFin;
    if (alloc_acv(ctx, ctx->NOC + 1))             /* Y radius */
        goto IVREGFin;


    ctx->IVRXC = ctx->AcX;
    ctx->IVRXR = ctx->AcU;
    ctx->IVRYC = ctx->AcY;
    ctx->IVRYR = ctx->AcV;
    ctx->IVRXSC = 0.0;
    ctx->IVRXSCA = 0.0;
    ctx->IVRXSR = 0.0;
    ctx->IVRXSA = 0.0;
    ctx->IVRYSC = 0.0;
    ctx->IVRYSR = 0.0;

    nw = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        x = get_data(ctx, ctx->PMVIdx[0],i);
        y = get_data(ctx, ctx->PMVIdx[1],i);
        ctx->IVRYC[i] = (x + y) / 2.0;
        ctx->IVRYR[i] = (y - x) / 2.0;
        x = get_data(ctx, ctx->PMVIdx[2],i);
        y = get_data(ctx, ctx->PMVIdx[3],i);
        ctx->IVRXC[i] = (x + y) / 2.0;
        ctx->IVRXR[i] = (y - x) / 2.0;

        if (ctx->IVRXR[i] < -ctx->EPSI || ctx->IVRYR[i] < -ctx->EPSI) {
            printf1(ctx, "Error: no valid interval in record %d.\n",i + 1);
            goto IVREGFin;
        }
        if (fabs(ctx->IVRXC[i]) < ctx->IVRXR[i] - ctx->EPSI)
            nw++;

        ctx->IVRXSC += ctx->IVRXC[i];
        ctx->IVRXSCA += fabs(ctx->IVRXC[i]);
        ctx->IVRXSR += ctx->IVRXR[i];
        ctx->IVRXSA += fabs(ctx->IVRXC[i]) + ctx->IVRXR[i];
        ctx->IVRYSC += ctx->IVRYC[i];
        ctx->IVRYSR += ctx->IVRYR[i];
    }
    ctx->IVRXSC /= (double)ctx->NOC;
    ctx->IVRXSCA /= (double)ctx->NOC;
    ctx->IVRXSR /= (double)ctx->NOC;
    ctx->IVRXSA /= (double)ctx->NOC;
    ctx->IVRYSC /= (double)ctx->NOC;
    ctx->IVRYSR /= (double)ctx->NOC;
      
/** 
    tda_out("data\n"); 
    for (i = 0; i < NOC; ++i) {
        tda_out("%3d %f %f %f %f\n",i,IVRXC[i],IVRXR[i],IVRYC[i],IVRYR[i] );
    }
    tda_out("ivrysc=%g ivrysr=%g\n",IVRYSC,IVRYSR);
    tda_out("ivrxsc=%g ivrxsr=%g\n",IVRXSC,IVRXSR);
**/    

    /* CONTOUR PLOT =================================================== */

    if (ctx->PMOPT == 4) {

        if (ctx->PMNNFlg == 0 || ctx->PMNN1 < 1 || ctx->PMNN2 < 1)
            ctx->PMNN1 = ctx->PMNN2 = 10;

        printf1(ctx, "\nContour plot with %d,%d grid.\n",ctx->PMNN1,ctx->PMNN2);
        if (check_ps(ctx, 2))
            goto IVREGFin;

        if (ctx->PMN > 0) {
            if (ctx->PMAFlg == 0) {
                printf1(ctx, "Error: need a parameter (beta center and radius).\n");
                goto IVREGFin;
            }
            if (alloc_actmp(ctx, ctx->PMN + 2))                               
                goto IVREGFin;

            if (ctx->PMA2 < 0.0)
                ctx->PMA2 = 0.0;

            printf1(ctx, "Level generation with beta=%g,%g\n",ctx->PMA1,ctx->PMA2);
            if (ctx->PMA1 <= ctx->PA1[0] - ctx->EPSI || ctx->PMA1 >= ctx->PA2[0] + ctx->EPSI ||
                ctx->PMA2 <= ctx->PA1[1] - ctx->EPSI || ctx->PMA2 >= ctx->PA2[1] + ctx->EPSI) {  
                printf1(ctx, "Error: beta not inside grid.\n");
                goto IVREGFin;
            }
            d1 = (ctx->PMA1 - ctx->PA1[0]) / (double)ctx->PMN;
            d2 = (ctx->PA2[0] - ctx->PMA1) / (double)ctx->PMN;
            br = ctx->PMA2;
            i  = 0;

            if (d1 > d2) {
                bc = ctx->PA1[0] + d1 / 2.0;
                while (bc <= ctx->PMA1 + ctx->EPSI1) {
                    f = ivreg2f(ctx, bc,br);
                    printf1(ctx, "Level %2d  bc=%g Function value: %20.12f\n",i + 1,bc,f);
                    ctx->AcTmp[i++] = f;
                    if (i == ctx->PMN)
                        break;
                    bc += d1;
                }
            }   
            else {
                bc = ctx->PA2[0] - d2 / 2.0;
                while (bc >= ctx->PMA1 - ctx->EPSI1) {
                    f = ivreg2f(ctx, bc,br);
                    printf1(ctx, "Level %2d  bc=%g Function value: %20.12f\n",i + 1,bc,f);
                    ctx->AcTmp[i++] = f;
                    if (i == ctx->PMN)
                        break;
                    bc -= d2;
                }
            }   
        }
        else {
            if (ctx->PMNTP < 1) {
                printf1(ctx, "Error: need level definition with x parameter.\n");
                goto IVREGFin;
            }
        }        
        /***
        for (i = 0; i < PMN; ++i)
            tda_out("%d %g\n",i,AcTmp[i]);
        ***/

        if (ivreg2p(ctx, ctx->PMN,ctx->AcTmp))
            goto IVREGFin;
    
        goto IVREGFin;
    }
           
    /* OLS REGRESSION ================================================= */

    printf1(ctx, "\nOLS regression with center values.\n");

    if (alloc_acw(ctx, 3 * ctx->NOC + 1))     /* data for lsei */
        goto IVREGFin;
    if (alloc_acz(ctx, 3))               /* parameter */
        goto IVREGFin;

    for (i = 0; i < ctx->NOC; ++i) {
        ctx->AcW[i * 3 + 1] = 1.0;
        ctx->AcW[i * 3 + 2] = ctx->IVRXC[i];
        ctx->AcW[i * 3 + 3] = ctx->IVRYC[i];
    }
    r = lsei(ctx, ctx->AcW,0,ctx->NOC,0,2,ctx->AcZ,0,&x,&lsn,&rn,&n);
    if (r) {
        if (r == -4)  
            p_err(ctx, -2,1);
        else  
            printf1(ctx, "LSEI error return: %d\n",err);
        goto IVREGFin;  
    }
    printf1(ctx, "Rank of least squares matrix: %d\n",rn);
    printf1(ctx, "Norm of least squares residuals:");
    rt_printf1_d(ctx, ctx->PMFmtS,lsn);
    printf1(ctx, "\nParameter:");
    rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcZ[1]);
    rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcZ[2]);
    newline(ctx);

    /* DIRECT CALCULATION ============================================= */

    ivreg2d(ctx);

    if (nw > 0) {
        printf1(ctx, "Warning: Cases with x-radius greater than x-center: %d\n\n",nw);
    }

    /* FUNCTION MINIMIZATIOBN WITH DIRECT SEARCH  ===================== */

    if (ctx->PMOPT == 2) {

    printf1(ctx, "Direct search minimization.\n");
    printf1(ctx, "Step length: %g\n",ctx->SLen);
    printf1(ctx, "Reduction factor: %g\n",ctx->SRed);
    printf1(ctx, "Minimum of step length: %g\n\n",ctx->TOLS);

    bc = ctx->AcZ[2];
    br = 0.0;
    ctx->AcZ[1] = ctx->AcZ[2] = 0.0;

    if (get_dsv(ctx, 2,ctx->AcZ,0,ctx->AcZ,ctx->AcZ,0))   /* get starting values */ 
        goto IVREGFin;
       
    if (ctx->DSVFlg) {
        bc = ctx->AcZ[1];
        if (ctx->AcZ[2] >= 0.0)
            br = ctx->AcZ[2];
    }
    printf1(ctx, "Starting values:");
    rt_printf1_d(ctx, ctx->PMFmtS,bc);
    rt_printf1_d(ctx, ctx->PMFmtS,br);
    newline(ctx);

        f = ivreg2f(ctx, bc,br);

        printf1(ctx, "Function value: ");
        rt_printf1_d(ctx, ctx->PMFmtS,f);
        newline(ctx);
        
        f = ivreg2m(ctx, &bc,&br);

    }
    else if (ctx->PMOPT == 3) {

        ctx->NParm = 2;                  /* number of parameters */
        ctx->MINA = 2;                   /* always this algorithm */
        ctx->PMDOPT = 0;                 /* without derivatives */
        ctx->CCTyp = 0;                  /* no covariance matrix */
        set_mlopt(ctx);                /* adjust options */         

        if (ml_init(ctx, 1,1,0))         /* init function minimization */
            goto IVREGFin;

        ctx->Par[1] = ctx->AcZ[2];
        ctx->Par[2] = 0.0;

f = ivreg2f(ctx, ctx->Par[1],ctx->Par[2]);
tda_out("f=%g\n",f);


        if (ffmin(ctx, IVRMOD1,1,0,1)) {                               
            printf1(ctx, "ERROR\n");       
            exit(0);
        }
        prn_mlres(ctx, 1);   
        tda_out("FMIN=%g\n",ctx->FMin);
        tda_out("par=%20.12f %20.12f\n",ctx->Par[1],ctx->Par[2]);

        /* was a bare exit(0) here -- killed the whole process on
           success, discarding any command after this one in the
           current file (or any enclosing one) with no message.
           Its sibling debug exit(0) calls elsewhere in this file are
           commented out, as this one should have been; return to the
           caller instead. */
        err = 0;
        goto IVREGFin;
    }
    /* f, bc, br are only ever assigned above for opt=2 (direct
       search) via ivreg2m(); opt=1 (direct calculation, the default)
       falls through here with them still uninitialized, so this
       summary used to print a garbage "Best function value"
       (observed as -nan) and a garbage Beta -- the real opt=1 answer
       is already printed per-domain, with the valid domain(s) marked
       '*', by ivreg2d() just above. opt=3 (mina) prints its own
       FMIN=/par= via tda_out and returns before reaching here.
       Printed only where f/bc/br are real. */
    if (ctx->PMOPT == 2) {
        ar = ac = 0.0;
        for (i = 0; i < ctx->NOC; ++i) {
            ac += ctx->IVRYC[i] - i_center(ctx, ctx->IVRXC[i],ctx->IVRXR[i],bc,br);
            ar += ctx->IVRYR[i] - i_radius(ctx, ctx->IVRXC[i],ctx->IVRXR[i],bc,br);
        }
        ac /= (double)ctx->NOC;
        ar /= (double)ctx->NOC;

        printf1(ctx, "\nBest function value:  ");
        rt_printf1_d(ctx, ctx->PMFmtS,f);
        printf1(ctx, "\nAlpha (center-radius):");
        rt_printf1_d(ctx, ctx->PMFmtS,ac);
        rt_printf1_d(ctx, ctx->PMFmtS,ar);
        printf1(ctx, "\nBeta  (center-radius):");
        rt_printf1_d(ctx, ctx->PMFmtS,bc);
        rt_printf1_d(ctx, ctx->PMFmtS,br);
        newline(ctx);
    }

    newline(ctx);
    err = 0;

IVREGFin:
    con_free(ctx);
    ml_init(ctx, 0,0,0);      
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  ivreg2m(bc,br)  Function minimization for ivreg2. This function uses    */
/*                  a simple direct search method in 2-dim. space.          */
/*                  Starting values are given by bc and br, final values    */
/*                  are returned in these locations. Also return the        */
/*                  mimimal function value.                                 */
/*                                                                          */
/*                  Initial step size SLen, reduction factor: SRed.         */
/*                  Iterations are terminated when step size < TOLS.        */

double ivreg2m(TDAContext *ctx, double *bc,double *br)
{
    double sl,rf,f,f0,f1,f2,f3,f4,f5,f6,f7,f8;
    double bc0,bc1,bc2,bc3,bc4,bc5,bc6,bc7,bc8;
    double br0,br1,br2,br3,br4,br5,br6,br7,br8;

    sl = ctx->SLen;
    rf = ctx->SRed;
/**
    tda_out("hier sl=%g rf=%g tols=%g bc=%g br=%g\n",sl,rf,TOLS,*bc,*br  );
**/

    f  = ivreg2f(ctx, *bc,*br);

    if (ctx->SILENTFlg < 2) {
        printfe(ctx, "\nBeta center          Beta radius         Function value      Step size\n");
        printfe(ctx, "%19.12e %19.12e %19.12e %19.12e\n",*bc,*br,f,sl);
    }
        
    while (1) {
                 
        bc1 = *bc - sl; 
        bc5 = *bc + sl;
        bc3 = bc7 = *bc; 
        bc2 = bc8 = *bc - sl / 2.0;
        bc4 = bc6 = *bc + sl / 2.0;

        br1 = br5 = *br;
        br3 = *br + sl;
        br7 = *br - sl;
        br2 = br4 = *br + sl / 2.0;
        br6 = br8 = *br - sl / 2.0;

        f1 = ivreg2f(ctx, bc1,br1);
        f2 = ivreg2f(ctx, bc2,br2);
        f3 = ivreg2f(ctx, bc3,br3);
        f4 = ivreg2f(ctx, bc4,br4);
        f5 = ivreg2f(ctx, bc5,br5);
/**
        tda_out("\nslen=%14.12f\n",sl);
        tda_out("bc=%16.12f br=%16.12f f =%19.14f\n",*bc,*br,f);
        tda_out("bc=%16.12f br=%16.12f f1=%19.14f\n",bc1,br1,f1);
        tda_out("bc=%16.12f br=%16.12f f2=%19.14f\n",bc2,br2,f2);
        tda_out("bc=%16.12f br=%16.12f f3=%19.14f\n",bc3,br3,f3);
        tda_out("bc=%16.12f br=%16.12f f4=%19.14f\n",bc4,br4,f4);
        tda_out("bc=%16.12f br=%16.12f f5=%19.14f\n",bc5,br5,f5);
**/
        f0 = f1;
        bc0 = bc1;
        br0 = br1;
        if (f2 < f0) {
            f0 = f2;
            bc0 = bc2;
            br0 = br2;
        }
        if (f3 < f0) {
            f0 = f3;
            bc0 = bc3;
            br0 = br3;
        }
        if (f4 < f0) {
            f0 = f4;
            bc0 = bc4;
            br0 = br4;
        }
        if (f5 < f0) {
            f0 = f5;
            bc0 = bc5;
            br0 = br5;
        }
/**
        tda_out("k=%d f0=%19.14f\n",k,f0);
**/
        if (br6 >= -ctx->EPSI) {
            f6 = ivreg2f(ctx, bc6,br6);
            f8 = ivreg2f(ctx, bc8,br8);
  
/**
            tda_out("bc=%16.12f br=%16.12f f6=%19.14f\n",bc6,br6,f6);
            tda_out("bc=%16.12f br=%16.12f f8=%19.14f\n",bc8,br8,f8);
**/
            if (f6 < f0) {
                f0 = f6;
                bc0 = bc6;
                br0 = br6;

            }
            if (f8 < f0) {
                f0 = f8;
                bc0 = bc8;
                br0 = br8;
            }
        }
        if (br7 >= -ctx->EPSI) {
            f7 = ivreg2f(ctx, bc7,br7);
/*
            tda_out("bc=%16.12f br=%16.12f f7=%19.14f\n",bc7,br7,f7);
*/
            if (f7 < f0) {
                f0 = f7;
                bc0 = bc7;
                br0 = br7;
            }
        }
        if (sl <= ctx->TOLS) {
            break;
        }
        if (f0 >= f - ctx->EPSI) {
            sl *= rf; 
        }
        else {
            f = f0;
            *bc = bc0;
            *br = br0;
            if (ctx->SILENTFlg < 2)  
                printfe(ctx, "%19.12e %19.12e %19.12e %19.12e\n",*bc,*br,f,sl);
        }
    }
    if (ctx->SILENTFlg < 2)  
        printfe(ctx, "%19.12e %19.12e %19.12e %19.12e\n",*bc,*br,f,sl);

    return(f);
}

/* ------------------------------------------------------------------------ */
/*  ivreg2f(bc,br)      Return function value (using data as set up by      */
/*                      ivreg2().)                                          */

double ivreg2f(TDAContext *ctx, double bc,double br) 
{
    register int i;                    
    double xc,xr,yc,yr,d1,d2,xbc,xbr,xbcs,xbrs;
         
    d1 = d2 = xbcs = xbrs = 0.0;

    for (i = 0; i < ctx->NOC; ++i) {
        yc = ctx->IVRYC[i];
        yr = ctx->IVRYR[i];
        xc = ctx->IVRXC[i];
        xr = ctx->IVRXR[i];

        xbc = i_center(ctx, xc,xr,bc,br);
        xbr = i_radius(ctx, xc,xr,bc,br);
             
        d1 += xbc * (xbc - 2.0 * (yc - ctx->IVRYSC));
        d2 += xbr * (xbr - 2.0 * (yr - ctx->IVRYSR));

        xbcs += xbc;
        xbrs += xbr;
    }
    d1 -= xbcs * xbcs / (double)ctx->NOC;
    d2 -= xbrs * xbrs / (double)ctx->NOC;
    return(d1 + d2);
}

/*--------------------------------------------------------------------------*/
/*  ivreg2p(nlev,flev)   Contour plot for ivreg2(). If nlev > 0 the array   */
/*                       flev contains level.                               */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int ivreg2p(TDAContext *ctx, int nlev,double *flev)
{
    register int i,j;
    int nn,err,pcbita,pcfva;
    double bc,br,*fl;
           
    pcbita = pcfva = 0;
    err = -1;

    if (nlev > 0)  
        fl = flev;
    else {
        nlev = ctx->PMNTP;
        fl = ctx->PMTP;
    }


    if (ctx->PMNNFlg == 0 || ctx->PMNN1 < 1 || ctx->PMNN2 < 1)
        ctx->PMNN1 = ctx->PMNN2 = 10;
      
    ctx->PCDX = ctx->UXLen / (double) ctx->PMNN1;
    ctx->PCDY = ctx->UYLen / (double) ctx->PMNN2;
    ctx->PCNX = ctx->PMNN1 + 1;
    ctx->PCNY = ctx->PMNN2 + 1;

    if (!(ctx->PCFV = (float *) calloc((size_t)(ctx->PCNX) * (size_t)(ctx->PCNY) + 1,sizeof(float)))) {
        p_err(ctx, -2,1);
        goto IVCPFin;
    }
    memrq(ctx, ctx->PCNX * ctx->PCNY + 1,sizeof(float));
    pcfva = 1;

    nn = 2 * ctx->PCNX * ctx->PCNY * nlev / 8 + 1;

    if (!(ctx->PCBitM = (char *)calloc((size_t)(nn),sizeof(char)))) {
        p_err(ctx, -2,1);
        goto IVCPFin;
    }
    memrq(ctx, nn,sizeof(char));
    pcbita = 1;

    /* get function values in PCFV */
      
    for (i = 1; i <= ctx->PCNY; ++i) {
        for (j = 1; j <= ctx->PCNX; ++j) {
     
            bc = ctx->PA1[0] + (double)(j - 1) * ctx->PCDX;
            br = ctx->PA1[1] + (double)(i - 1) * ctx->PCDY;
            ctx->PCFV[(j - 1) * ctx->PCNX + i] = (float)ivreg2f(ctx, bc,br);
        }
    }
    fprintf(ctx->PSFd,"\n%%#%d: ivreg contour plot\n",++ctx->PSONUM);

    cont_plot(ctx, nlev,fl,0);
                    
    err = 0;

IVCPFin:
    if (pcfva) {
        free((char *)ctx->PCFV);
        memrq(ctx, -ctx->PCNX * ctx->PCNY - 1,sizeof(float));
    }
    if (pcbita) {
        free((char *)ctx->PCBitM);
        memrq(ctx, -nn,sizeof(char));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ivr1_fn()   Calculate function value for IVRMOD1.                       */
/*              Use parameter values in TPar[] (i=1,NParm).                 */
/*              Return function value in FTmp.                              */
/*                                                                          */
/*              VErsion for ivreg3                                          */


/*              Return 0 if OK, 1 if error in function                      */

int ivr1_fn(TDAContext *ctx)
{
    register int i;
    int err;
    double bc,br,xc,xr,yc,yr;

    bc = ctx->TPar[1];
    br = fabs(ctx->TPar[2]);
    ctx->FTmp = 0.0;

    for (i = 0; i < ctx->NOC; ++i) {
        yc = ctx->IVRYC[i];
        yr = ctx->IVRYR[i];
        xc = ctx->IVRXC[i];
        xr = ctx->IVRXR[i];

        ctx->FTmp += (xc * bc + xr * br) * (bc * (xc - ctx->IVRXSC) +
                                       br * (xr - ctx->IVRXSR) - 
                                      2.0 * (yc - ctx->IVRYSC)) +
                (xr * bc + xc * br) * (bc * (xr - ctx->IVRXSR) +
                                       br * (xc - ctx->IVRXSC) - 
                                      2.0 * (yr - ctx->IVRYSR)); 
    }
    err = checkov(ctx);        
    if (err) {
        printf1(ctx, "\nError in function evaluation.\n");
        printf1(ctx, "Numerical overflow.\n");
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ivreg2dp(n,na,bc,br)                                                    */

void ivreg2dp(TDAContext *ctx, int n,int na,double bc,double br) 
{
    register int i;
    double f,ac,ar; 

    printf1(ctx, "%4d ",n);
    if (na < 0) {
        printf1(ctx, "\n");
        return;
    }
    else if (na > 0)
        printf1(ctx, "* ");
    else
        printf1(ctx, "  ");

    prnchar(ctx, ' ',11 - ctx->PMFmt1,0);
    rt_printf1_d(ctx, ctx->PMFmtS,bc);
    prnchar(ctx, ' ',11 - ctx->PMFmt1,0);
    rt_printf1_d(ctx, ctx->PMFmtS,br);

    if (na > 0) {
        ar = ac = 0.0;
        for (i = 0; i < ctx->NOC; ++i) {
            ac += ctx->IVRYC[i] - i_center(ctx, ctx->IVRXC[i],ctx->IVRXR[i],bc,br);
            ar += ctx->IVRYR[i] - i_radius(ctx, ctx->IVRXC[i],ctx->IVRXR[i],bc,br);
        }
        ac /= (double)ctx->NOC;
        ar /= (double)ctx->NOC;

        prnchar(ctx, ' ',12 - ctx->PMFmt1,0);
        rt_printf1_d(ctx, ctx->PMFmtS,ac);
        prnchar(ctx, ' ',12 - ctx->PMFmt1,0);
        rt_printf1_d(ctx, ctx->PMFmtS,ar);
    }
    else
        prnchar(ctx, ' ',2 * imax(ctx, 13,ctx->PMFmt1 + 1),0);

    f = ivreg2f(ctx, bc,br);
    prnchar(ctx, ' ',16 - ctx->PMFmt1,0);
    rt_printf1_d(ctx, ctx->PMFmtS,f);
    newline(ctx);             
}

/* ------------------------------------------------------------------------ */
/*  ivreg2d()                                                               */

int ivreg2d(TDAContext *ctx)
{
    register int i;
    int na;
    double tmp,a,b,c,d,e,br,bc,xc,xr,xa,xca,yc,yr,vc,vr,va,vca,wc,wr;
    double vzs,vcs,vrs,vas,vcas,vcar,vcara,vss,wzs,wcs,wrs,urs,wcar,wcars;

    vzs = vcs = vrs = vas = vcas = wzs = wcs = wrs = urs = 0.0;
    wcars = wcar = vss = vcar = vcara = 0.0;

    for (i = 0; i < ctx->NOC; ++i) {
        yc = ctx->IVRYC[i];
        yr = ctx->IVRYR[i];
        xc = ctx->IVRXC[i];
        xr = ctx->IVRXR[i];
        xca = fabs(xc);
        xa = xca + xr;

        vc = (xc - ctx->IVRXSC);
        vca = (xca - ctx->IVRXSCA);
        vr = (xr - ctx->IVRXSR);
        va = (xa - ctx->IVRXSA);
        wc = (yc - ctx->IVRYSC);
        wr = (yr - ctx->IVRYSR);

        vcs += xc * vc;
        vrs += xr * vr;
        vas += xa * va;
        vcas += xca * vca;
        vcar += xca * vr;
        vcara += xr * vca;
        vzs += (xc + dsign(ctx, xc) * xr) * (vc + dsign(ctx, xc) * vr);
        vss += dsign(ctx, xc) * (xc * vr + xr * vc);

        wcs += xc * wc;
        wrs += xr * wr;
        wzs += (xc + dsign(ctx, xc) * xr) * wc;                     
        wcar += xca * wr;
        wcars += dsign(ctx, xc) * xr * wc;
        urs += xa * wr;


    }
    printf1(ctx, "\nDomain ");
    prnchar(ctx, ' ',ctx->PMFmt1 - 11,0);
    printf1(ctx, "Beta center ");
    prnchar(ctx, ' ',ctx->PMFmt1 - 11,0);
    printf1(ctx, "Beta radius ");
    prnchar(ctx, ' ',ctx->PMFmt1 - 12,0);
    printf1(ctx, "Alpha center ");
    prnchar(ctx, ' ',ctx->PMFmt1 - 12,0);
    printf1(ctx, "Alpha radius   ");
    prnchar(ctx, ' ',ctx->PMFmt1 - 14,0);
    printf1(ctx, "Function value\n");

    /* domain 1 */

    na = 0;
    bc = br = 0.0;
    tmp = vcs + vrs;
    if (tmp == 0.0)
        na = -1;
    else {
        bc = (wcs + wrs) / tmp;
        if (bc >= -ctx->EPSI1)
            na = 1;
    }
    ivreg2dp(ctx, 1,na,bc,br);

    /* domain 2 */

    na = 0;
    bc = br = 0.0;
    tmp = vcs + vrs;
    if (tmp == 0.0)
        na = -1;
    else {
        bc = (wcs - wrs) / tmp;
        if (bc <= ctx->EPSI1)
            na = 1;
    }
    ivreg2dp(ctx, 2,na,bc,br);

    /* domain 3 */

    na = 0;
    bc = br = 0.0;
    if (vas == 0.0 || vzs == 0.0)
        na = -1;
    else {
        bc = wzs / vzs;
        br = urs / vas;
        if (br >= -ctx->EPSI1 && fabs(bc) <= br + ctx->EPSI1)
            na = 1;
    }
    ivreg2dp(ctx, 3,na,bc,br);

    /* domain 4 */

    na = 0;
    bc = br = 0.0;

    a = vcs + vrs;
    b = vrs + vcas;
    c = vss + vcar + vcara;
    d = wcs + wrs;
    e = wcars + wcar;
    tmp = c * c - 4.0 * a * b;

    if (tmp == 0.0)
        na = -1;
    else {
        bc = (2.0 * c * e - 4.0 * b * d) / tmp;
        br = (2.0 * c * d - 4.0 * a * e) / tmp;
        if (br >= -ctx->EPSI1 && bc >= -ctx->EPSI1 && bc >= br - ctx->EPSI1)
            na = 1;
    }
    ivreg2dp(ctx, 4,na,bc,br);

    /* domain 5 */

    na = 0;
    bc = br = 0.0;

    a = vcs + vrs;
    b = vrs + vcas;
    c = -vss - vcar - vcara;
    d = wcs - wrs;
    e = -wcars + wcar;
    tmp = c * c - 4.0 * a * b;

    if (tmp == 0.0)
        na = -1;
    else {
        bc = (2.0 * c * e - 4.0 * b * d) / tmp;
        br = (2.0 * c * d - 4.0 * a * e) / tmp;
        if (br >= -ctx->EPSI1 && bc <= ctx->EPSI1 && -bc >= br - ctx->EPSI1)
            na = 1;
    }
    ivreg2dp(ctx, 5,na,bc,br);




    newline(ctx);
    return(0);
}


