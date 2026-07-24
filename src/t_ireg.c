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

/* ------------------------------------------------------------------------ */
/*  functions in t_ireg                                                     */

int ilsreg(void);
int imreg(void);

int ivls(void);
void mcross(int n,int m,double *x,double *a);

int ivreg(void);
void ivreg_m(int opt,int np,double *ml,double *mh);
void ivreg_y(int opt,int np,double xml,double xmh);
int ivreg_n(int opt,int np);
int ivreg_upd(int opt,int np,int i);
double ivreg_b0(int np,double *xy);              
int ivreg_bb(int n,double *xl,double *xh,double *bl,double *bh);

int ivreg_min(int typ,int opt,int mopt,int narg,int nbmax,double *par,double *lb,
    double *ub,int mxit,double tolbw,double tolfd,double tolfe,int gc);
int ivreg_h(int opt,int n,double *xl,double *xh,
    double *rl,double *rh,int deriv,double *gl,double *gh);

int ivar1f(int n,double *xl,double *xh,double *ip,int *idx,double *v0,double *v1);
int ivreg_z(int np,int rp,int ns,int opt);

int ivreg1(void);
int ivreg1_h(int mopt,int opt,int n,double *xl,double *xh,
    double *rl,double *rh,int deriv,double *gl,double *gh,double fminp);

double i_center(double ac,double ar,double bc,double br);
double i_radius(double ac,double ar,double bc,double br);
void i_bvar(double xl,double xh,double bl,double bh,double *ul,double *uh);
double i_norm(int opt,double xl,double xh,double yl,double yh);
double ivreg1_nrm(int opt,double xl,double xh,double yl,double yh,
    double al,double ah,double bl,double bh);

int ivreg2(void);
double ivreg2m(double *bc,double *br);
double ivreg2f(double bc,double br);
int ivreg2p(int nlev,double *flev);
int ivr1_fn(void);
int ivreg2d(void);
void ivreg2dp(int n,int na,double bc,double br);


/*  Global parameters and variables                                         */

double IVRXML = 0.0;        /* lower bound mean of X                        */
double IVRXMU = 0.0;        /* upper bound mean of X                        */
double IVRYML = 0.0;        /* lower bound mean of Y                        */
double IVRYMU = 0.0;        /* upper bound mean of Y                        */
double IVRBETA = 0.0;       /* current value of beta                        */
double IVRBL = 0.0;         /* lower bound for beta                         */
double IVRBU = 0.0;         /* upper bound for beta                         */
double IVRVARL = 0.0;       /* lower bound of variance                      */
double IVRVARU = 0.0;       /* upper bound of variance                      */
double *IVRX;               /* pointer to data used by ivreg_h              */
double *IVRY;
int *IVRXP;
int *IVRYP;

double *IVRXC;              /* center of X                                  */
double *IVRXR;              /* radius of X                                  */
double *IVRYC;              /* center of Y                                  */
double *IVRYR;              /* radius of Y                                  */
double IVRXSC = 0.0;        /* mean of X center                             */  
double IVRXSCA = 0.0;       /* mean of abs(X center)                        */  
double IVRXSR = 0.0;        /* mean of X radius                             */  
double IVRXSA = 0.0;        /* mean of abs(X)                               */  
double IVRYSC = 0.0;        /* mean of Y center                             */  
double IVRYSR = 0.0;        /* mean of Y radius                             */  

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

int ilsreg(void)
{
    register int i,j;
    int err,il,ih,ix;     
    double x,yl,yh,ylm,yhm,v,w,xm,xxm;
    double tmp,bmin,bmax,amin,amax,ni,nw;
         
    err = -1;         
    if (check_cmd(1))
        return(-1);

    printf1("LS regression with interval-valued dependent variable. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 6,4,1))     /* get parameters */
        goto ILSFin;

    if (PMOPT != 2)
        PMOPT = 1;

    if (PMNV != 3) {            /* need three variables */
        p_err(-1,1);
        goto ILSFin;
    }
    il = PMVIdx[0];
    ih = PMVIdx[1];
    ix = PMVIdx[2];

    if (alloc_acx(NOC + 1))     /* x variable */ 
        goto ILSFin;
    if (alloc_acy(NOC + 1))     /* lower bound */ 
        goto ILSFin;
    if (alloc_acw(NOC + 1))     /* width */
        goto ILSFin;

    tmp = v = xxm = xm = ylm = yhm = 0.0;
    for (i = 0; i < NOC; ++i) {
        x  = get_data(ix,i);
        yl = get_data(il,i);
        yh = get_data(ih,i);
        if (yh < yl) {
            printf1("Error in case %d, found: %lg,%lg\n",i+1,yl,yh);
            goto ILSFin;
        }
        AcX[i] = x;
        AcY[i] = yl;
        AcW[i] = yh - yl;
        ylm += yl;
        yhm += yh;
        xm += x;
        xxm += x * x;
        v += x * yl;

/*      printf1("i=%d x=%lg  yl=%lg w=%lg \n",i,x,  yl,AcW[i]); */      
    }
    w = (double)NOC * xxm - xm * xm;
    v = (double)NOC * v - xm * ylm;
    xm /= (double)NOC;
    ylm /= (double)NOC;
    yhm /= (double)NOC;

    if (w <= 0.0) {
        printf1("Error: indep. variable has zero variance.\n");
        goto ILSFin;
    }
    yl = yh = bmin = bmax = 0.0;
          
    for (i = 0; i < NOC; ++i) {
        tmp = AcX[i] - xm;   
        if (tmp > 0.0) {
            bmax += tmp * AcW[i];
            yl += AcY[i];
            yh += AcY[i] + AcW[i];
        }
        else if (tmp < 0.0) {
            bmin += tmp * AcW[i];
            yh += AcY[i];
            yl += AcY[i] + AcW[i];
        }
    }
    bmax = ((double)NOC * bmax + v) / w;
    bmin = ((double)NOC * bmin + v) / w;
    amax = yh / (double)NOC - xm * bmax;
    amin = yl / (double)NOC - xm * bmin;

    printf1("\nMean of independent variable: ");
    printf1(PMTFmtS,xm);
    printf1("\nRange of mean of dependent variable: ");
    printf1(PMTFmtS,ylm);
    printf1(PMTFmtS,yhm);
    printf1("\nBeta minimum: ");
    printf1(PMTFmtS,bmin);
    printf1("  corresponding alpha: ");
    printf1(PMTFmtS,amin);
    printf1("\nBeta maximum: ");
    printf1(PMTFmtS,bmax);
    printf1("  corresponding alpha: ");
    printf1(PMTFmtS,amax);
    newline();

    /* calculate bounds for given x values */

    if (PMNTP > 0) {

        newline();
        prnchar(' ',PMTFmt1 - 10,0);
        printf1("   x-value ");
        prnchar(' ',PMTFmt1 - 10,0);
        printf1(" y minimum ");
        prnchar(' ',PMTFmt1 - 10,0);
        printf1(" y maximum\n");
        prnchar('-',3 * PMTFmt1 + 2,1);

        for (j = 0; j < PMNTP; ++j) {
            x = PMTP[j];
            ni = 1.0 / (double)NOC;
            nw = (double)NOC * (x - xm) / w;

            yl = yh = bmin = bmax = 0.0;

            for (i = 0; i < NOC; ++i) {
                tmp = nw * (AcX[i] - xm) + ni;

                if (tmp > 0.0)  
                    bmax += tmp * AcW[i];
                else if (tmp < 0.0)  
                    bmin += tmp * AcW[i];
            }
            yh = ylm + (x - xm) * v / w + bmax;
            yl = ylm + (x - xm) * v / w + bmin;

            printf1(PMTFmtS,x);
            printf1(PMTFmtS,yl);
            printf1(PMTFmtS,yh);
            newline();

            if (PMF1Def) {
                fprintf(PMF1d,PMTFmtS,x);
                fprintf(PMF1d,PMTFmtS,yl);
                fprintf(PMF1d,PMTFmtS,yh);
                fprintf(PMF1d,"\n");
            }   
        }
        if (PMF1Def)  
            printf1("\n%d records written to: %s\n",PMNTP,PMF1dName);

    }
    newline();
    err = 0;

ILSFin:
    p_clean();
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

int imreg(void)
{
    register int i,j;
    int err,ilx,ihx,ily,ihy,n,nn,nrec;
    double xl,xh,yl,yh,ymin,ymax,a,b;
         
    err = -1;         
    if (check_cmd(1))
        return(-1);

    printf1("Conditional means of interval-valued data. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,4,1))     /* get parameters */
        goto IMREGFin;

    if (PMNV != 4) {            /* need four variables */
        p_err(-1,1);
        goto IMREGFin;
    }
    if (PMFmtF == 0)            /* default print format */
        pmfmt(10,4);

    ily = PMVIdx[0];
    ihy = PMVIdx[1];
    ilx = PMVIdx[2];
    ihx = PMVIdx[3];

    if (alloc_acx(NOC + 1))     /* x lower bound */ 
        goto IMREGFin;
    if (alloc_acu(NOC + 1))     /* x upper bound */ 
        goto IMREGFin;
    if (alloc_acy(NOC + 1))     /* y lower bound */
        goto IMREGFin;
    if (alloc_acv(NOC + 1))     /* y upper bound */
        goto IMREGFin;
    if (alloc_actmp(2 * NOC + 1))     /* induced partition */
        goto IMREGFin;

    n = 0;
    for (i = 0; i < NOC; ++i) {
        xl = get_data(ilx,i);
        xh = get_data(ihx,i);
        if (xh <= xl + EPSI1) {
            printf1("Error in case %d, found: %lg,%lg\n",i+1,xl,xh);
            goto IMREGFin;
        }
        yl = get_data(ily,i);
        yh = get_data(ihy,i);
        if (yh <= yl + EPSI1) {
            printf1("Error in case %d, found: %lg,%lg\n",i+1,yl,yh);
            goto IMREGFin;
        }
        AcX[i] = xl;
        AcU[i] = xh;
        AcY[i] = yl;
        AcV[i] = yh;
        AcTmp[n++] = xl;
        AcTmp[n++] = xh;
    }
    if (sortd(n,AcTmp,0))
        goto IMREGFin;
   
    nn = 1;
    for (i = 1; i < n; ++i) {           /* AcTmp contains induced partition */
        if (AcTmp[i] > AcTmp[i - 1])
            AcTmp[nn++] = AcTmp[i];
    }
    nrec = 0;

    for (i = 0; i < nn; ++i) {
        a = AcTmp[i];
        b = AcTmp[i + 1];
             
        n = 0;
        ymin = ymax = 0.0;
        for (j = 0; j < NOC; ++j) {
            if (AcX[j] <= a && AcU[j] >= b) {
                ymin += AcY[j];
                ymax += AcV[j];
                n++;
            }
        }
        if (n > 1) {
            ymin /= (double)n;
            ymax /= (double)n;
        }
        if (PMF1Def) {
            fprintf(PMF1d,"%6d ",i + 1);
            fprintf(PMF1d,PMFmtS,a);
            fprintf(PMF1d,PMFmtS,ymin);
            fprintf(PMF1d,PMFmtS,ymax);
            fprintf(PMF1d,"\n");

            fprintf(PMF1d,"%6d ",i + 1);
            fprintf(PMF1d,PMFmtS,b);
            fprintf(PMF1d,PMFmtS,ymin);
            fprintf(PMF1d,PMFmtS,ymax);
            fprintf(PMF1d,"\n");
            nrec += 2;
        }
        if (i >= nn - 2)
            break;
    }
    if (nrec > 0)
        printf1("%d records written to: %s\n",nrec,PMF1dName);

    newline();
    err = 0;

IMREGFin:
    p_clean();
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

int ivls(void)
{
    register int i,j,k;
    int err,n,nv,m,r,iter,ni,nf,m1,nlp;
    int *z,*p,*xlidx,*xhidx;
    double tol,xl,xh,x,d,v,c,c1;
    double *acmat,*admat,*bcvec,*bdvec,*rmat,*gmat,*x0vec,*gvec,*lpmat;
    double *xlvec,*xhvec,*dvec,*d1vec,*lpx,*lpy;
         
    err = -1;         
    if (check_cmd(1))
        return(-1);

    TOLP = 1.0e-8;

    printf1("Linear interval equations. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,4,1))     /* get parameters */
        goto IVLSFin;

    nv = PMNV / 2;
    if (2 * nv != PMNV) {   
        printf1("Error: need even number of variables on right-hand side.\n");
        goto IVLSFin;
    }
    if (nv < 2) {   
        printf1("Error: need at least two interval-valued variables.\n");
        goto IVLSFin;
    }
    n = NOC;
    m = nv - 1;

    if (n < m) {
        printf1("Error: number of cases must not be less than number of variables.\n");
        goto IVLSFin;
    }
    if (PMSC <= EPSI)  
        PMSC = 0.01;

    if (MxIter < 1 || MxItFlg == 0)
        MxIter = 1000;

    printf1("\nTolerance for preliminary bounds (sc): %lg\n",PMSC);
    printf1("Max number of iterations: %d\n",MxIter);
    printf1("Tolerance (tolp): %lg\n\n",TOLP);

    if (alloc_acx(2 * n * m + 4))   
        goto IVLSFin;
    acmat = AcX;
    admat = AcX + n * m + 2;

    if (alloc_acy(3 * n + 6))   
        goto IVLSFin;
    bcvec = AcY;
    bdvec = AcY + n + 2;
    dvec  = AcY + 2 * n + 4;

    if (alloc_acz(m * n + m * m + 4))   
        goto IVLSFin;
    rmat = AcZ;
    gmat = AcZ + m * n + 2;

    if (alloc_acw(5 * m + 10))       
        goto IVLSFin;
    x0vec = AcW;
    gvec  = AcW + m + 2;
    d1vec = AcW + 2 * m + 6;
    xlvec = AcW + 3 * m + 8;
    xhvec = AcW + 4 * m + 10;

    for (i = 0; i < n; ++i) {
        k = 0;
        for (j = 0; j < nv; ++j) {
            xl = get_data(PMVIdx[k++],i);
            xh = get_data(PMVIdx[k++],i);
            if (xh < xl) {
                printf1("Error: no valid interval in record %d.\n",i + 1);
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
        printf("yc=%lf ",bcvec[i]);
        printf("yd=%lf  xc: ",bdvec[i]);
        for (j = 1; j <= m; ++j)  
            printf("%lf ",acmat[i * m + j]);
        printf(" xd: ");
        for (j = 1; j <= m; ++j)  
            printf("%lf ",admat[i * m + j]);
        newline();
    }
***/

    mcross(n,m,acmat,gmat);
    r = ginv(m,m,gmat);
    if (r < 0) {
        printf1("Error in matrix inversion. Cannot continue.\n");
        goto IVLSFin;
    }
    else if (r != m) {
        printf1("Error: insufficient rank=%d of data matrix.\n",r);
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
      
    tol = DBLMAX / 1000.0;

    for (i = 0; i < m; ++i)
        d1vec[i] = 0.0;

    for (iter = 0; iter < MxIter; ++iter) {

        for (i = 0; i < m; ++i)
            dvec[i] = d1vec[i];

        for (i = 0; i < m; ++i) {
            x = 0.0;
            for (j = 1; j <= m; ++j)  
                x += gmat[i * m + j] * dvec[j - 1];
            d1vec[i] = x + gvec[i] + PMSC;
        }
        d = 0.0;
        for (i = 0; i < m; ++i)  
            d = dmax(d,fabs(dvec[i] - d1vec[i]));

        if (d < PMSC || d > tol)
            break;
    }
    if (d >= PMSC) {
        printf1("Error: no convergence (d=%lg).\n",d);
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
    printf1("Preliminary bounds\n");
    for (i = 0; i < m; ++i) {
        printf1(PMTFmtS,xlvec[i]);
        printf1(PMTFmtS,xhvec[i]);
        newline();
    }
    newline();

    /* try to find optimal enclosure */

    m1 = m + 1;

    if (alloc_acu((2 * n + 1) * m1 + 2))   
        goto IVLSFin;
    lpmat = AcU;

    if (alloc_acv(m + 2 * n + 4))   
        goto IVLSFin;
    lpx = AcV;
    lpy = AcV + m + 2;

    if (alloc_acn(2 * m + 4))   
        goto IVLSFin;
    z = AcN;
    p = AcN + m + 2;

    if (alloc_acm(2 * m + 4))   
        goto IVLSFin;
    xlidx = AcM;
    xhidx = AcM + m + 2;
             
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
    printf1("Trying to find optimal enclosure.\n");
    printf1("Need to solve %d lp problems.\n\n",nlp);
             
    for (i = 0; i < ni; ++i)
        z[p[i]] = -1;
                    
    while (1) {
        /***/
        printf("z: ");
        for (i = 0; i < m; ++i)
            printf("%2d ",z[i]);
        newline();
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
            printf("\nlpmat fuer min\n");
            for (i = 0; i <= 2 * n; ++i) {
                for (j = 1; j <= m1; ++j)
                    printf("%lf ",lpmat[i * m1 + j]);
                newline();
            }
            ***/
         
            r = lpf1(2 * n,m,0,lpmat,lpx,lpy,&v,TOLP);       

            /***/
            printf("min  r=%d\n",r);
            printf("lpx: ");
            for (i = 1; i <= m; ++i)
                printf("%g ",lpx[i]);
            newline();
            /***/

            if (r == 0) {
                if (xlidx[k] == 0) {
                    xlvec[k] = lpx[k + 1] * c1;
                    xlidx[k] = 1;
                }
                else
                    xlvec[k] = dmin(xlvec[k],lpx[k + 1] * c1);
            }
            lpmat[k + 1] = c1;

            /***
            printf("\nlpmat fuer max\n");
            for (i = 0; i <= 2 * n; ++i) {
                for (j = 1; j <= m1; ++j)
                    printf("%f ",lpmat[i * m1 + j]);
                newline();
            }
            ***/

            r = lpf1(2 * n,m,0,lpmat,lpx,lpy,&v,TOLP);       

            /***/
            printf("max r=%d\n",r);
            printf("lpx: ");
            for (i = 1; i <= m; ++i)
                printf("%g ",lpx[i]);
            newline();
            /***/

            if (r == 0) {
                if (xhidx[k] == 0) {
                    xhvec[k] = lpx[k + 1] * c1;
                    xhidx[k] = 1;
                }
                else
                    xhvec[k] = dmax(xhvec[k],lpx[k + 1] * c1);
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
        printf1("Cannot find a solution.\n");
    }   
    else {
        nf = 1;
        printf1("Final bounds\n");
        for (i = 0; i < m; ++i) {
            printf1(PMTFmtS,xlvec[i]);
            printf1(PMTFmtS,xhvec[i]);
            newline();
            if (xlvec[i] > xhvec[i] + EPSI)
                nf = 0;
        }
        if (nf == 0)
            printf1("\nThe solution set seems to be empty.\n");
    }
    newline();
    err = 0;

IVLSFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  mcross(n,m,x,a)  Return cross product x'x in a.                         */

void mcross(int n,int m,double *x,double *a)
{
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

int ivreg(void)
{
    register int i,j,k;
    int err,n,np,r,nx,ny,nx0,ny0,iter,rp;
    double xl,xh,beta;
           
    err = -1;         
    if (check_cmd(1))
        return(-1);

    TOLF  = 1.e-6;       /* def tolerance for beta iteration                 */
    TOLBW = 1.e-4;       /* def tolerance for box length                     */
    TOLFD = 1.e-10;      /* def tolerance for function range                 */
    TOLFE = 1.e-10;      /* def tolerance for global minimum                 */

    printf1("Regression with interval data. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,4,1))     /* get parameters */
        goto IVREGFin;

    if (PMOPT < 1 || PMOPT > 3)
        PMOPT = 1;

    if (PMN < 2)
        PMN = 20;

    if (PMNV != 4) {   
        printf1("Error: need four variables on right-hand side.\n");
        goto IVREGFin;
    }
    if (PMNS != 1)
        PMNS = 0;

    if (PMFmtF == 0)        /* default print format */
        pmfmt(13,6);

    n = NOC;
    np = 2 * n; 

    if (PMNS == 0)
        printf1("\nMinimization.\n");
    else
        printf1("\nMaximization.\n");
    printf1("Number of data points: %d\n",n);
    if (n < 2) {
        printf1("Error: need at least two data points.\n");
        goto IVREGFin;
    }
    printf1("Tolerance for box width: %g\n\n",TOLBW);

    /* allocate RPar, RParL, RParU (i = 0,...,np-1) */

    if (alloc_par(np,1))
        goto IVREGFin;

    /* temporarily used for variance calculation */

    if (alloc_acx(n + 1))              
        goto IVREGFin; 
    if (alloc_acy(n + 1))  
        goto IVREGFin; 
    if (alloc_actmp(2 * n + 1))  
        goto IVREGFin; 
    if (alloc_ack(n + 1))  
        goto IVREGFin; 

    /*  get data into RParL and RParU */

    k = 0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < 2; ++j) {
            xl = get_data(PMVIdx[j * 2],i);
            xh = get_data(PMVIdx[j * 2 + 1],i);
            if (xh < xl) {
                printf1("Error: no valid interval in record %d.\n",i + 1);
                goto IVREGFin;
            }
            RParL[k] = xl;
            RParU[k] = xh;
            RPar[k] = (xl + xh) / 2.0;
            k++;

            if (j == 1) {
                AcX[i] = xl;
                AcY[i] = xh;
            }
        }
    }
/**
    printf("data RPar, RParl and RParU\n"); 
    for (i = 0; i < np; ++i) {
        printf("%3d %f %f %f\n",i,RPar[i],RParL[i],RParU[i]);
    }
**/
    /* ============================================================ */

    ivreg_m(0,np,&IVRXML,&IVRXMU);      /* calculate mean of X */

    ivreg_y(PMNS,np,IVRXML,IVRXMU);     /* fix y values */

    ivreg_m(1,np,&IVRYML,&IVRYMU);      /* calculate mean of Y */

    ny = ivreg_n(1,np);                 /* number of fixed y values */

    printf1("First step: fixed %d values of Y.\n",ny);
    printf1("Mean Y = [ ");
    printf1(PMFmtS,IVRYML);  
    printf1(", ");
    printf1(PMFmtS,IVRYMU);  
    printf("]\nMean X = [ ");
    printf1(PMFmtS,IVRXML);  
    printf1(", ");
    printf1(PMFmtS,IVRXMU);  
    printf("]\n");

    /* calculate variance of X */

    r = ivar1f(n,AcX,AcY,AcTmp,AcK,&IVRVARL,&IVRVARU);
    if (r) {
        printf1("\nError: cannot find bounds for variance of X.\n");
        goto IVREGFin;
    }
    printf1("Var  X = [ ");
    printf1(PMFmtS,IVRVARL);  
    printf1(", ");
    printf1(PMFmtS,IVRVARU);  
    printf1("]\n\n");

/**
    printf("data RPar, RParl and RParU\n"); 
    for (i = 0; i < np; ++i) {
        printf("%3d %f %f %f\n",i,RPar[i],RParL[i],RParU[i]);
    }
**/
    /* ============================================================ */

    printf1("Second step: fixing values of X (and Y).\n");

    IVRBETA = ivreg_b0(np,RPar);                /* preliminary beta */
    ivreg_bb(np,RParL,RParU,&IVRBL,&IVRBU);     /* bounds for beta */

    if (PMNS == 0)                              /* minimization */
        IVRBU = IVRBETA;
    else    
        IVRBL = IVRBETA;

    printf1("Iteration     bounds of current beta          bounds of x-mean  x-fixed  y-fixed\n");
    printf1("%7d  %14.6f %12.6f %12.6f %12.6f %8d %8d\n",0,IVRBL,IVRBU,IVRXML,IVRXMU,0,ny);
           
    iter = nx0 = 0;
    ny0 = ny;
            
    while (1) {                                 

        for (i = 0; i < np; i += 2) {
            if (fabs(RParU[i + 1] - RParL[i + 1]) <= TOLBW)
                continue;

            ivreg_upd(PMNS,np,i);
        }
        if (nx0 < n) {
            ivreg_m(0,np,&IVRXML,&IVRXMU);      /* calculate new mean of X */
            nx = ivreg_n(0,np);                 /* number of fixed x values */
        }
        if (ny0 < n) {                          /* also update y values */
            ivreg_y(PMNS,np,IVRXML,IVRXMU);     /* fix y values */
            ivreg_m(1,np,&IVRYML,&IVRYMU);      /* calculate mean of Y */
            ny = ivreg_n(1,np);                 /* number of fixed y values */
        }
        if (nx == nx0 && ny == ny0)  
            break;
        nx0 = nx;
        ny0 = ny;

        /* try to update variance of X */

        k = 0;
        for (i = 0; i < np; i += 2) {
            AcX[k] = RParL[i + 1];
            AcY[k] = RParU[i + 1];
            k++;
        }
        ivar1f(n,AcX,AcY,AcTmp,AcK,&IVRVARL,&IVRVARU);            
/*
        printf("new var %g %g\n",IVRVARL,IVRVARU);
*/
        /* update beta */

        beta = IVRBETA;
        IVRBETA = ivreg_b0(np,RPar);                
        ivreg_bb(np,RParL,RParU,&IVRBL,&IVRBU);  
        if (PMNS == 0)              
            IVRBU = IVRBETA;
        else    
            IVRBL = IVRBETA;

        printf1("%7d  %14.6f %12.6f %12.6f %12.6f %8d %8d\n",++iter,IVRBL,IVRBU,IVRXML,IVRXMU,nx,ny);



    }
           
    /* calculate rp = remaining parameters */

    rp = 0; 
    for (i = 0; i < np; i += 2) {
        if (fabs(RParU[i] - RParL[i]) > TOLBW)  
            rp++;
        if (fabs(RParU[i + 1] - RParL[i + 1]) > TOLBW)  
            rp++;
    }
    printf1("\nRemaining coordinates: %d\n",rp);
              
    /* ============================================================ */

    if (PMOPT == 1 || rp == 0)
        goto IVREGCont; 
                   
    printf1("\nThird step: additional optimization (opt=%d).\n",PMOPT);

    if (PMOPT == 2) {
        ivreg_z(np,rp,PMN,PMNS);
        goto IVREGFin;

    }
    else if (PMOPT == 3) {

        if (MxIter < 1 || MxItFlg == 0)
            MxIter = 100;

        if (PMNBOX < 1)         /* max number of boxes */
            PMNBOX = 100;

        /* reorganize x- and y-values */

        if (alloc_acx(n + 1))  
            goto IVREGFin; 
        if (alloc_acy(n + 1))  
            goto IVREGFin; 
        if (alloc_aci(n + 1))  
            goto IVREGFin; 
        if (alloc_acj(n + 1))  
            goto IVREGFin; 

        IVRX  = AcX;
        IVRY  = AcY;
        IVRXP = AcI;
        IVRYP = AcJ;

        j = k = 0;
        for (i = 0; i < np; i += 2) {

            IVRY[k] = RPar[i];
            IVRX[k] = RPar[i + 1];

            IVRYP[k] = -1;
            if (fabs(RParU[i] - RParL[i]) > TOLBW) {
                IVRYP[k] = j;
                RParL[j] = RParL[i];
                RParU[j] = RParU[i];
                RPar[j]  = RPar[i];
                j++;
            }
            IVRXP[k] = -1;
            if (fabs(RParU[i + 1] - RParL[i + 1]) > TOLBW) {
                IVRXP[k] = j;
                RParL[j] = RParL[i + 1];
                RParU[j] = RParU[i + 1];
                RPar[j]  = RPar[i + 1];
                j++;
            }
            k++;
        }
        if (alloc_list(PMNBOX,rp)) 
            goto IVREGFin;

        printf1("Number of parameters: %d\n",rp);              
        printf1("Maximum number of iterations: %d\n",MxIter);              
        printf1("Maximum number of boxes: %d\n",PMNBOX);
        printf1("Tolerance for box width: %g\n",TOLBW);
        printf1("Tolerance for function range: %g\n",TOLFD);
        printf1("Tolerance for global minimum: %g\n\n",TOLFE);
      
        r = ivreg_min(0,PMNS,0,rp,PMNBOX,RPar,RParL,RParU,MxIter,TOLBW,TOLFD,TOLFE,1);
        if (r) {
            if (r == 1)
                printf1("Exceeded maximum number of boxes.\n");
            printf1("Current best value of beta:");
            printf1(PMFmtS,GO_FMIN);
            newline();
            goto IVREGFin;
        }
        igmin_res(TOLFE,0);
    }    

IVREGCont:
    if (PMF1Def) {          /* write to output file */

        if (PMOPT == 1 || rp == 0) {
            for (i = 0; i < np; i += 2) {
                fprintf(PMF1d,PMFmtS,RPar[i]);
                fprintf(PMF1d,PMFmtS,RPar[i + 1]);
                fprintf(PMF1d,"\n");
            }
        }
        else if (PMOPT == 3) {   

            for (i = 0; i < n; ++i) {
                if ((k = IVRYP[i]) >= 0)
                    fprintf(PMF1d,PMFmtS,RPar[k]);
                else
                    fprintf(PMF1d,PMFmtS,IVRY[i]);

                if ((k = IVRXP[i]) >= 0)
                    fprintf(PMF1d,PMFmtS,RPar[k]);
                else
                    fprintf(PMF1d,PMFmtS,IVRX[i]);
                fprintf(PMF1d,"\n");
            }                     
        }
        printf1("%d records written to: %s\n",n,PMF1dName);
    }
    newline();
    err = 0;

IVREGFin:
    alloc_list(0,0);
    alloc_par(0,0);
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ivreg_m(opt,np,ml,mh)                                                   */
/*                                                                          */
/*  Calculate bounds of mean value. If opt = 0 then X, otherwise Y.         */
/*  Return bounds in ml and mh.                                             */

void ivreg_m(int opt,int np,double *ml,double *mh)
{
    register int i;

    *ml = *mh = 0.0;
    for (i = 0; i < np; i += 2) {
        if (opt == 1) {
            *ml += RParL[i];
            *mh += RParU[i];
        }
        else {
            *ml += RParL[i + 1];
            *mh += RParU[i + 1];
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

void ivreg_y(int opt,int np,double xml,double xmh)
{
    register int i;

    for (i = 0; i < np; i += 2) {
        if (RParU[i + 1] < xml) {
            if (opt == 0)   
                RPar[i] = RParL[i] = RParU[i];
            else
                RPar[i] = RParU[i] = RParL[i];
        }
        else if (RParL[i + 1] > xmh) {
            if (opt == 0)   
                RPar[i] = RParU[i] = RParL[i];
            else
                RPar[i] = RParL[i] = RParU[i];
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  ivreg_n(opt,np)     If opt = 0 return number of fixed x values,         */
/*                      otherwise number of fixed y values.                 */

int ivreg_n(int opt,int np)
{
    register int i,n;

    n = 0;
    for (i = 0; i < np; i += 2) {
        if (opt == 1) {
            if (fabs(RParU[i] - RParL[i]) <= TOLBW)
                n++;          
        }
        else {
            if (fabs(RParU[i + 1] - RParL[i + 1]) <= TOLBW)
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

int ivreg_upd(int opt,int np,int i)
{
    int n;
    double xl,xh,yl,yh,tl,th;

    n = 0;
    i_sub(RParL[i],RParU[i],IVRYML,IVRYMU,&yl,&yh);
    i_sub(RParL[i + 1],RParU[i + 1],IVRXML,IVRXMU,&xl,&xh);
    i_mul(2.0 * IVRBL,2.0 * IVRBU,xl,xh,&tl,&th);
       
    if (yh < tl) {
        if (opt == 0)                      
            RPar[i + 1] = RParL[i + 1] = RParU[i + 1];
        else
            RPar[i + 1] = RParU[i + 1] = RParL[i + 1];
        return(1);
    }
    else if (yl > th) {
        if (opt == 0)                      
            RPar[i + 1] = RParU[i + 1] = RParL[i + 1];
        else
            RPar[i + 1] = RParL[i + 1] = RParU[i + 1];
        return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  ivreg_b0(np,xy)     Return beta based on interval centers.              */

double ivreg_b0(int np,double *xy)               
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
        printf1("ERROR in ivreg_b0. Zero variance.\n");
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

int ivreg_bb(int n,double *xl,double *xh,double *bl,double *bh)
{
    register int i;   
    double nn,sl,sh,tl,th,ul,uh;
       
    nn = (double)(n / 2);
    ul = uh = 0.0;

    for (i = 0; i < n; i += 2) {
        i_sub(nn * xl[i],nn * xh[i],nn * IVRYML,nn * IVRYMU,&sl,&sh);
        i_mul(sl,sh,xl[i + 1],xh[i + 1],&tl,&th);
        i_add(ul,uh,tl,th,&sl,&sh);
        ul = sl;
        uh = sh;
    }
    if (i_div(ul,uh,nn * nn * IVRVARL,nn * nn * IVRVARU,bl,bh)) {
        printf1("FATAL ERROR in ivreg_bb.\n");
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

int ivreg_min(int typ,int opt,int mopt,int narg,int nbmax,double *par,double *lb,
    double *ub,int mxit,double tolbw,double tolfd,double tolfe,int gc)
{
    register int i,j;
    int err,r,nb,is,iter,ja,jb,jj,first;
    double fl,fu,tmp,w,f,fminp;
    double *lptra,*lptrb,*uptra,*uptrb;

    GO_IT = 0;              /* number of iterations performed */
    GO_FN = 0;              /* number of function evaluations */
    GO_IFN = 0;             /* number of inclusion function evaluations */
    GO_NBU = 0;             /* number of boxes used */
    GO_NB = 0;              /* number of boxes in final list */

    err = -1;

    for (i = 0; i < nbmax; ++i)
        GO_FLG[i] = 0;

    if (alloc_actmp(narg + 1))  
        goto IGMINFin; 
    if (gc) {
        if (alloc_acu(narg + 1))    /* lower bounds of gradient */
            goto IGMINFin; 
        if (alloc_acv(narg + 1))    /* upper bounds of gradient */
            goto IGMINFin; 
    }
    err = 0;

    lptra = GO_LB;                      /* put initial box on list */
    uptra = GO_UB;

    for (i = 0; i < narg; ++i) {
        lptra[i] = lb[i];
        uptra[i] = ub[i];
    }
    GO_FLG[0] = -1;     /* not yet processed */

    /* calculate inclusion function with initial boxes */
   
    GO_IFN++;

    if (typ == 0)
        r = ivreg_h(1,narg,lptra,uptra,&fl,&fu,0,AcU,AcV);                                      
    else
        r = ivreg1_h(mopt,1,narg,lptra,uptra,&fl,&fu,0,AcU,AcV,fminp);                                      

    if (r) {
        printf1("Error in evaluating inclusion function.\n");
        err = -2;
        goto IGMINFin;
    }
    GO_LBF[0] = fl;
    GO_UBF[0] = fu;
    
    /* calculate function at starting values */

    GO_FN++;
    if (typ == 0)
        r = ivreg_h(0,narg,par,&tmp,&fminp,&tmp,0,AcU,AcV);                                      
    else
        r = ivreg1_h(mopt,0,narg,par,&tmp,&fminp,&tmp,0,AcU,AcV,fminp);                                      

    if (r) {
        printf1("Error in evaluating function.\n");
        err = -2;
        goto IGMINFin;
    }
    nb = 1;     /* number of boxes */
    iter = 0;                                     

    if (SILENTFlg < 2)
        printfe("\n  Iter    Function Value       NBox  FCall  IFCall\n");
    
    while (++iter <= mxit) {       

        if (SILENTFlg < 2)
            printfe("%5d  %20.13e  %6d %6d %7d\n",iter,fminp,nb,GO_FN,GO_IFN);

        /* find box with lowest lower, or largest upper, function bound,
           and check for boxes that can be dropped */

        ja = -1;
        first = 1;
        for (i = 0; i < nb; ++i) {

            if (GO_FLG[i]) {    /* check all boxes */
  
                if (opt == 0) {
                    if (GO_LBF[i] > fminp)  
                        GO_FLG[i] = 0;
 
                    else if (GO_FLG[i] < 0) {
                        if (first) {      
                            f = GO_LBF[i];
                            ja = i;
                            first = 0;    
                        }
                        else if (f > GO_LBF[i]) {
                            f = GO_LBF[i];
                            ja = i;
                        }
                    }
                }   
                else {

                    if (GO_UBF[i] < fminp)  
                        GO_FLG[i] = 0;
                    else if (GO_FLG[i] < 0) {
                        if (first) {      
                            f = GO_UBF[i];
                            ja = i;
                            first = 0;
                        }                 
                        else if (f < GO_UBF[i]) {
                            f = GO_UBF[i];
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
            if (GO_FLG[i] == 0) {
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
            if (GO_NBU < nb)
                GO_NBU = nb;
        }
        lptra = GO_LB + ja * narg;
        uptra = GO_UB + ja * narg;

        lptrb = GO_LB + jb * narg;
        uptrb = GO_UB + jb * narg;
   
        is = 0;                         /* coordinate with largest width */
        w = uptra[0] - lptra[0];

        for (i = 1; i < narg; ++i) {
            tmp = uptra[i] - lptra[i];
            if (w < tmp) {
                w = tmp;
                is = i;
            }
        }

        if (w <= TOLBW) {               /* temporarily accept this box */
            GO_FLG[ja] = 1;
            continue;
        }

        for (i = 0; i < narg; ++i) {    /* copy box ja to jb */
            lptrb[i] = lptra[i];
            uptrb[i] = uptra[i];
        }

        tmp = lptra[is] + w / 2.0;
        uptra[is] = lptrb[is] = tmp;
        GO_FLG[jb] = -1;

        for (jj = 0; jj < 2; ++jj) {           /* for both boxes */

            if (jj)   
                ja = jb;

            lptra = GO_LB + ja * narg;
            uptra = GO_UB + ja * narg;

            /* evaluate inclusion function */

            GO_IFN++;
            if (typ == 0)
                r = ivreg_h(1,narg,lptra,uptra,&fl,&fu,gc,AcU,AcV);                                      
            else
                r = ivreg1_h(mopt,1,narg,lptra,uptra,&fl,&fu,gc,AcU,AcV,fminp);                                      
            if (r) {
                printf1("Error in evaluating inclusion function.\n");
                err = -2;
                goto IGMINFin;
            }
            if (opt == 0) {
                if (fl > fminp) {
                    GO_FLG[ja] = 0;             /* free box */
                    if (ja == nb - 1)
                        nb--;    
                    continue;
                }
            }
            else {
                if (fu < fminp) {
                    GO_FLG[ja] = 0;             /* free box */
                    if (ja == nb - 1)
                        nb--;    
                    continue;
                }
            }
            GO_LBF[ja] = fl;
            GO_UBF[ja] = fu;

            if (opt == 0) {
                if (fminp > fu) {
                    fminp = fu;
                    igmin_cpar(narg,uptra,par);
                }
            }
            else {
                if (fminp < fl) {
                    fminp = fl;
                    igmin_cpar(narg,lptra,par);
                }
            }
            if (fabs(fu - fl) < tolfd)      /* accept for solution */
                GO_FLG[ja] = 1;
   
            if (gc) {                       /* gradient check */

                r = 0;
                for (j = 0; j < narg; ++j) {
                    if (opt == 0) {
                        if (AcU[j] > 0.0) {
                            if (lptra[j] > lb[j]) {
                                r = 1;
                                break;
                            }
                            uptra[j] = lptra[j];
                            r = -1;
                        }
                        else if (AcV[j] < 0.0) {
                            if (uptra[j] < ub[j]) {
                                r = 1;
                                break;
                            }
                            lptra[j] = uptra[j];
                            r = -1;
                        }
                    }
                    else {
                        if (AcU[j] > 0.0) {
                            if (uptra[j] < ub[j]) {
                                r = 1;
                                break;
                            }
                            lptra[j] = uptra[j];
                            r = -1;
                        }
                        else if (AcV[j] < 0.0) {
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
                    GO_FLG[ja] = 0;             /* free box */
                    if (ja == nb - 1)
                        nb--;    
                    continue;
                }
                else if (r == -1) {         /* update inclusion function */
                    GO_IFN++;
                    if (typ == 0)
                        r = ivreg_h(1,narg,lptra,uptra,&fl,&fu,0,AcU,AcV);                                      
                    else
                        r = ivreg1_h(mopt,1,narg,lptra,uptra,&fl,&fu,0,AcU,AcV,fminp);                                      
                    if (r) {
                        printf1("Error in evaluating inclusion function.\n");
                        err = -2;
                        goto IGMINFin;
                    }

                    if (opt == 0) {
                        if (fl > fminp) {
                            GO_FLG[ja] = 0;             /* free box */
                            if (ja == nb - 1)
                                nb--;    
                            continue;
                        }
                    }
                    else {
                        if (fu < fminp) {
                            GO_FLG[ja] = 0;             /* free box */
                            if (ja == nb - 1)
                                nb--;    
                            continue;
                        }
                    }
                    GO_LBF[ja] = fl;
                    GO_UBF[ja] = fu;
                }
            }

            /* update fminp with function value for midpoint of a box with
               least lower bound */
                   
            for (j = 0; j < narg; ++j)
                AcTmp[j] = (lptra[j] + uptra[j]) / 2.0;

            /* calculate function at midpoint */

            GO_FN++;
            if (typ == 0)
                r = ivreg_h(0,narg,AcTmp,&tmp,&f,&tmp,0,AcU,AcV);                                      
            else
                r = ivreg1_h(mopt,0,narg,AcTmp,&tmp,&f,&tmp,0,AcU,AcV,fminp);                                      

            if (r) {
                printf1("Error in evaluating function.\n");
                err = -2;
                goto IGMINFin;
            }
            if (opt == 0) {
                if (fminp > f) {
                    fminp = f;
                    igmin_cpar(narg,AcTmp,par);
                }
            }   
            else {
                if (fminp < f) {
                    fminp = f;
                    igmin_cpar(narg,AcTmp,par);
                }
            }   
        }
    }
    GO_NB = nb;
    GO_IT = iter;
    GO_FMIN = fminp;

    if (PMProtFDef) {
        fprintf(PMProtFd,"IGMIN optimization.\n");
        prval("Best function value",GO_FMIN);
        prvec("Parameters",narg,par - 1);
        fprintf(PMProtFd,"\n");
        
        if (opt == 0)
            tmp = GO_FMIN - tolfe;
        else
            tmp = GO_FMIN + tolfe;

        ja = 1;
        for (i = 0; i < nb; ++i) {
            if (GO_FLG[i] == 0)  
                continue;

            if (GO_FLG[i] >= 1) {

                if (opt == 0) {                 /* minimization */
                    if (GO_LBF[i] >= tmp)         
                        GO_FLG[i] = 2;
                }
                else {
                    if (GO_UBF[i] <= tmp)         
                        GO_FLG[i] = 2;
                }
            }
            lptra = GO_LB + i * narg;
            uptra = GO_UB + i * narg;

            w = uptra[0] - lptra[0];
            for (j = 1; j < narg; ++j)  
                w = dmax(w,uptra[j] - lptra[j]);

            fprintf(PMProtFd,"Box   Width                 Lower function bound  Upper function bound  Acceptance\n");
            fprintf(PMProtFd,"%3d  %21.14e %21.14e %21.14e  %6d\n",
                                 ja++,w,GO_LBF[i],GO_UBF[i],GO_FLG[i]);

            fprintf(PMProtFd,"Parameter vector\nLB: ");
            for (j = 0; j < narg; ++j)  
                fprintf(PMProtFd,PMPFmtS,lptra[j]);
            fprintf(PMProtFd,"\nUB: ");
            for (j = 0; j < narg; ++j)  
                fprintf(PMProtFd,PMPFmtS,uptra[j]);
            fprintf(PMProtFd,"\n\n");
        }
    }

IGMINFin:
    alloc_acu(0);           
    alloc_acv(0);           
    alloc_actmp(0);           
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ivreg_h(opt,n,xl,xh,rl,rh,deriv,gl,gh)                                  */
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

int ivreg_h(int opt,int n,double *xl,double *xh,
    double *rl,double *rh,int deriv,double *gl,double *gh)
{
    register int i,k;   
    double nn,tmp,sxl,sxh,syl,syh,sxxl,sxxh,sxyl,sxyh;
    double xxl,xxh,yyl,yyh,sl,sh,tl,th,ul,uh,nl,nh;
       
/**    
printf("ivreg_h opt=%d\n",opt);

    for (i = 0; i < n; ++i) {
        printf1("i=%d xl=%20.12f ",i,xl[i]);
        if (opt)
            printf1("xh=%20.12f",xh[i]);
        newline();
    }     
**/
    nn = (double)NOC;     

    if (opt == 0) {
        sxl = syl = sxxl = sxyl = 0.0;
        for (i = 0; i < NOC; ++i) {
            if (IVRXP[i] >= 0)
                xxl = xl[IVRXP[i]];
            else
                xxl = IVRX[i];

            if (IVRYP[i] >= 0)
                yyl = xl[IVRYP[i]];
            else
                yyl = IVRY[i];

            sxl += xxl;
            syl += yyl;
            sxxl += xxl * xxl;
            sxyl += xxl * yyl;
        }
        tmp = nn * sxxl - sxl * sxl;
        if (tmp == 0.0) {
            printf1("\nFATAL ERROR in ivreg_h: zero variance.\n");
            exit(0);
        }
        *rl = (nn * sxyl - sxl * syl) / tmp;
        return(0);
    }
    sxl = syl = sxxl = sxyl = 0.0;
    sxh = syh = sxxh = sxyh = 0.0;
    *rl = *rh = 0.0;

    for (i = 0; i < NOC; ++i) {

        if (IVRXP[i] >= 0) {
            xxl = xl[IVRXP[i]];
            xxh = xh[IVRXP[i]];
        }
        else
            xxl = xxh = IVRX[i];

        if (IVRYP[i] >= 0) {
            yyl = xl[IVRYP[i]];
            yyh = xh[IVRYP[i]];
        }
        else
            yyl = yyh = IVRY[i];

        i_add(sxl,sxh,xxl,xxh,&tl,&th);
        sxl = tl;
        sxh = th;
        i_add(syl,syh,yyl,yyh,&tl,&th);
        syl = tl;
        syh = th;
        i_square(xxl,xxh,&sl,&sh);
        i_add(sxxl,sxxh,sl,sh,&tl,&th);
        sxxl = tl;
        sxxh = th;
        i_mul(yyl,yyh,xxl,xxh,&sl,&sh);
        i_add(sxyl,sxyh,sl,sh,&tl,&th);
        sxyl = tl;
        sxyh = th;
    }
    i_square(sxl,sxh,&sl,&sh);
    i_sub(nn * sxxl,nn * sxxh,sl,sh,&nl,&nh);
    i_mul(sxl,sxh,syl,syh,&sl,&sh);
    i_sub(nn * sxyl,nn * sxyh,sl,sh,&ul,&uh);

    /* variance must be positive */

    if (nl < IVRVARL) {
        nl = IVRVARL;
        if (nh < nl)
            nh = nl;
    } 
    i_div(ul,uh,nl,nh,rl,rh);          
/*       
printf(">>>>>>>>>> rl=%g rh=%g\n",*rl,*rh);
*/
    if (deriv) {

        if (nl <= 0.0 && nh >= 0.0) {
            for (i = 0; i < n; ++i) {
                gl[i] = -1.0;
                gh[i] =  1.0;
            }
        }
        else {
            for (i = 0; i < NOC; ++i) {

                if ((k = IVRYP[i]) >= 0) {   

                    if (IVRXP[i] >= 0) {
                        xxl = xl[IVRXP[i]];
                        xxh = xh[IVRXP[i]];
                    }
                    else
                        xxl = xxh = IVRX[i];

                    i_sub(nn * xxl,nn * xxh,sxl,sxh,&tl,&th);
                    gl[k] = tl;
                    gh[k] = th;

/*   printf("gradY k=%d gl=%g gh=%f\n",k,gl[k],gh[k]);    */


                }
                if ((k = IVRXP[i]) >= 0) {   

                    xxl = xl[k];
                    xxh = xh[k];
                       
                    if (IVRYP[i] >= 0) {
                        yyl = xl[IVRYP[i]];
                        yyh = xh[IVRYP[i]];
                    }
                    else
                        yyl = yyh = IVRY[i];

                    i_sub(nn * xxl,nn * xxh,sxl,sxh,&sl,&sh);
                    i_mul(2.0 * *rl,2.0 * *rh,sl,sh,&tl,&th);
                    i_sub(nn * yyl,nn * yyh,syl,syh,&sl,&sh);
                    i_sub(sl,sh,tl,th,&ul,&uh);
                    gl[k] = ul;
                    gh[k] = uh;

/*  printf("gradX k=%d gl=%g gh=%f\n",k,gl[k],gh[k]);     */
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

int ivreg_h1(int opt,int n,double *xl,double *xh,
    double *rl,double *rh,int deriv,double *gl,double *gh)
{
    register int i,k;   
    double nn,tmp,sxl,syl,sxxl,sxyl;
    double xxl,xxh,yyl,yyh,sl,sh,tl,th,ul,uh,nl,nh;
    double xml,xmh,yml,ymh;
       
/**    
printf("ivreg_h opt=%d\n",opt);

    for (i = 0; i < n; ++i) {
        printf1("i=%d xl=%20.12f ",i,xl[i]);
        if (opt)
            printf1("xh=%20.12f",xh[i]);
        newline();
    }     
**/
    nn = (double)NOC;     

    if (opt == 0) {
        sxl = syl = sxxl = sxyl = 0.0;
        for (i = 0; i < NOC; ++i) {
            if (IVRXP[i] >= 0)
                xxl = xl[IVRXP[i]];
            else
                xxl = IVRX[i];

            if (IVRYP[i] >= 0)
                yyl = xl[IVRYP[i]];
            else
                yyl = IVRY[i];

            sxl += xxl;
            syl += yyl;
            sxxl += xxl * xxl;
            sxyl += xxl * yyl;
        }
        tmp = nn * sxxl - sxl * sxl;
        if (tmp == 0.0) {
            printf1("\nFATAL ERROR in ivreg_h: zero variance.\n");
            exit(0);
        }
        *rl = (nn * sxyl - sxl * syl) / tmp;
        return(0);
    }

    /* opt = 1 */

    xml = xmh = yml = ymh = 0.0;
    for (i = 0; i < NOC; ++i) {

        if (IVRXP[i] >= 0) {
            xxl = xl[IVRXP[i]];
            xxh = xh[IVRXP[i]];
        }
        else
            xxl = xxh = IVRX[i];

        if (IVRYP[i] >= 0) {
            yyl = xl[IVRYP[i]];
            yyh = xh[IVRYP[i]];
        }
        else
            yyl = yyh = IVRY[i];
      
        xml += xxl;
        xmh += xxh;
        yml += yyl;
        ymh += yyh;
    }

    *rl = *rh = 0.0;
    ul = uh = 0.0;

    for (i = 0; i < NOC; ++i) {

        if (IVRXP[i] >= 0) {
            xxl = xl[IVRXP[i]];
            xxh = xh[IVRXP[i]];
        }
        else
            xxl = xxh = IVRX[i];

        if (IVRYP[i] >= 0) {
            yyl = xl[IVRYP[i]];
            yyh = xh[IVRYP[i]];
        }
        else
            yyl = yyh = IVRY[i];
      
        i_sub(nn * yyl,nn * yyh,yml,ymh,&sl,&sh);
        i_mul(sl,sh,xxl,xxh,&tl,&th);
        i_add(ul,uh,tl,th,&sl,&sh);
        ul = sl;
        uh = sh;
    }
    if (i_div(ul,uh,nn * nn * IVRVARL,nn * nn * IVRVARU,rl,rh)) {
        printf1("FATAL ERROR in ivreg_k.\n");
        exit(0);
    }

    if (deriv) {
        xml /= nn;
        xmh /= nn;
        yml /= nn;
        ymh /= nn;

        for (i = 0; i < NOC; ++i) {

            if (IVRXP[i] >= 0) {
                xxl = xl[IVRXP[i]];
                xxh = xh[IVRXP[i]];
            }
            else
                xxl = xxh = IVRX[i];

            if (IVRYP[i] >= 0) {
                yyl = xl[IVRYP[i]];
                yyh = xh[IVRYP[i]];
            }
            else
                yyl = yyh = IVRY[i];
          
            if ((k = IVRYP[i]) >= 0) {   

                i_sub(xxl,xxh,xml,xmh,&tl,&th);
                gl[k] = tl;
                gh[k] = th;
            }
            if ((k = IVRXP[i]) >= 0) {

                i_sub(yyl,yyh,yml,ymh,&sl,&sh);
                i_sub(xxl,xxh,xml,xmh,&tl,&th);
                i_mul(2.0 * *rl,2.0 * *rh,tl,th,&nl,&nh);
                i_sub(sl,sh,nl,nh,&tl,&th);
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

int ivar1f(int n,double *xl,double *xh,double *ip,int *idx,double *v0,double *v1)
{
    register int i,j,np;
    int nn,fin,r;   
    double tmp,mmin,mmax,umin,vmin,vmax,a,b,u,v,mmin1,mmax1;
         
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

    if (sortd(np,ip,0))
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

    umin = vmin = DBLMAX;
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
            tmp = ivarf(n,u,xl,xh);
            tmp /= (double)n;
            if (tmp < vmin) {
                vmin = tmp;
                umin = u;             
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

int ivreg_z(int np,int rp,int ns,int opt)
{
    register int i,j,k;
    int nx,first,iter,bf;
    double b;

    fprintf(stderr,"\nHeuristic optimization.\n");
    printf1("\nHeuristic optimization.\n");
    printf1("Block  Size              Beta\n");

    iter = 0;

    while (1) {     
        iter++;
        if (alloc_aci(rp + 1))   
            return(-1);          

        k = 0;
        for (i = 0; i < np; i += 2) {
            if (fabs(RParU[i] - RParL[i]) > TOLBW)
                AcI[k++] = -i;
            if (fabs(RParU[i + 1] - RParL[i + 1]) > TOLBW)
                AcI[k++] = i;
        }
        if (k == 0)
            break;
        nx = imin(k,ns);

        if (alloc_acj(nx + 1))  
            return(-1);     
        if (alloc_acl(nx + 1))  
            return(-1);     

        if (alloc_acu(nx + 1))  
            return(-1);           
        if (alloc_acv(nx + 1))  
            return(-1);           

        for (i = 0; i < nx; ++i) {
            j = AcI[i];
            if (j <= 0) {
                AcU[i] = RParL[-j];
                AcV[i] = RParU[-j];
            }
            else {
                AcU[i] = RParL[j + 1];
                AcV[i] = RParU[j + 1];
            }
        }
        first = 1;
        while (comb_nm1(2,nx,AcJ,first)) {
            first = 0;
            for (i = 0; i < nx; ++i) {
                j = AcI[i];
                if (j <= 0) {
                    if (AcJ[i] == 0)
                        RPar[-j] = RParL[-j] = RParU[-j] = AcU[i];
                    else
                        RPar[-j] = RParL[-j] = RParU[-j] = AcV[i];
                }
                else {
                    if (AcJ[i] == 0)
                        RPar[j + 1] = RParL[j + 1] = RParU[j + 1] = AcU[i];
                    else
                        RPar[j + 1] = RParL[j + 1] = RParU[j + 1] = AcV[i];
                }
            }
            b = ivreg_b0(np,RPar);   
            bf = 0;
            if (opt == 0 && b < IVRBETA) {  
                IVRBETA = b;
                bf = 1;
            }
            else if (opt != 0 && b > IVRBETA) {  
                IVRBETA = b;
                bf = 1;
            }
            if (bf) {
                fprintf(stderr,"Block %4d  Beta %20.8f\n",iter,IVRBETA);
                for (i = 0; i < nx; ++i)
                    AcL[i] = AcJ[i];
            }
        }

        /* fix best values */

        for (i = 0; i < nx; ++i) {
            j = AcI[i];
            if (j <= 0) {  
                if (AcL[i] == 0)
                    RPar[-j] = RParL[-j] = RParU[-j] = AcU[i];
                else
                    RPar[-j] = RParL[-j] = RParU[-j] = AcV[i];
            }
            else {  
                if (AcL[i] == 0)
                    RPar[j + 1] = RParL[j + 1] = RParU[j + 1] = AcU[i];
                else
                    RPar[j + 1] = RParL[j + 1] = RParU[j + 1] = AcV[i];
            }
        }
        printf1("%5d  %4d  %16.8f\n",iter,nx,IVRBETA);
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

int ivreg1(void)
{
    register int i,j,k;
    int err,n,np,npp,r;
    double xl,xh;
           
    err = -1;         
    if (check_cmd(1))
        return(-1);

    TOLF  = 1.e-6;       /* def tolerance for beta iteration                 */
    TOLBW = 1.e-4;       /* def tolerance for box length                     */
    TOLFD = 1.e-10;      /* def tolerance for function range                 */
    TOLFE = 1.e-10;      /* def tolerance for global minimum                 */

    printf1("Regression with interval data. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 6,4,1))     /* get parameters */
        goto IVREGFin;

    if (PMNV != 4) {   
        printf1("Error: need four variables on right-hand side.\n");
        goto IVREGFin;
    }
    if (PMFmtF == 0)        /* default print format */
        pmfmt(13,6);

    n = NOC;
    np = 2 * n; 
    npp = 4;

    printf1("\nNumber of data points: %d\n",n);
    if (n < 2) {
        printf1("Error: need at least two data points.\n");
        goto IVREGFin;
    }

    if (alloc_acx(2 * NOC + 1))     
        goto IVREGFin;
    if (alloc_acy(2 * NOC + 1))     
        goto IVREGFin;

    k = 0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < 2; ++j) {
            xl = get_data(PMVIdx[j * 2],i);
            xh = get_data(PMVIdx[j * 2 + 1],i);
            if (xh < xl) {
                printf1("Error: no valid interval in record %d.\n",i + 1);
                goto IVREGFin;
            }
            AcX[k] = xl;
            AcY[k] = xh;
            k++;
        }
    }
    /**
    printf("data\n"); 
    for (i = 0; i < np; ++i) {
        printf("%3d %f %f\n",i,AcX[i],AcY[i]);
    }
    **/ 
    if (MxIter < 1 || MxItFlg == 0)
        MxIter = 100;

    if (PMNBOX < 1)         /* max number of boxes */
        PMNBOX = 100;

    if (alloc_par(npp,1))
        goto IVREGFin;

    if (alloc_list(PMNBOX,npp)) 
        goto IVREGFin;

    printf1("Number of parameters: %d\n",npp);              
    printf1("Maximum number of iterations: %d\n",MxIter);              
    printf1("Maximum number of boxes: %d\n",PMNBOX);
    printf1("Tolerance for box width: %g\n",TOLBW);
    printf1("Tolerance for function range: %g\n",TOLFD);
    printf1("Tolerance for global minimum: %g\n\n",TOLFE);
      
    RParL[0] =  1.0;
    RParU[0] =  2.5;
    RParL[1] =   0.0;
    RParU[1] =   0.5;
    RParL[2] =   0.0;
    RParU[2] =   0.1;
    RParL[3] =   0.0;
    RParU[3] =   0.1;
    for (i = 0; i < npp; ++i)
        RPar[i] = (RParL[i] + RParU[i]) / 2.0;
      
    RPar[0] = 2.0889;      
    RPar[1] = 0.368047;
    RPar[2] = 0.0797882;
    RPar[3] = 0.05;        
             


    r = ivreg_min(1,0,0,npp,PMNBOX,RPar,RParL,RParU,MxIter,TOLBW,TOLFD,TOLFE,0);
    if (r) {
        if (r == 1)
            printf1("Exceeded maximum number of boxes.\n");
        printf1("Current best value of beta:");
        printf1(PMFmtS,GO_FMIN);
        newline();
        goto IVREGFin;
    }
    igmin_res(TOLFE,0);
    /***
    for (i = 0;i < npp; ++i)
        printf("i=%d par=%f\n",i,RPar[i]);
    ***/

          
    newline();
    err = 0;

IVREGFin:
    alloc_list(0,0);
    alloc_par(0,0);
    p_clean();
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

int ivreg1_h(int mopt,int opt,int n,double *pl,double *ph,
    double *rl,double *rh,int deriv,double *gl,double *gh,double fminp)
{
    register int i;   
    int np;
    double xc,xr,yc,yr,d1,d2;
    double acl,ach,arl,arh,bcl,bch,brl,brh,sl,sh,tl,th,ul,uh,vl,vh;
    double abl,abh,ax,d1l,d1h,d2l,d2h,ml,mh,ml1,mh1,ml2,mh2;

/**     
printf("ivreg1_h opt=%d mopt=%d\n",opt,mopt);

    for (i = 0; i < n; ++i) {
        printf1("i=%d xl=%20.12f ",i,pl[i]);
        if (opt)
            printf1("xh=%20.12f",ph[i]);
        newline();
    }     
**/       

    np = 2 * NOC;
        
    *rl = *rh = 0.0;

    if (opt == 0) {
   
        for (i = 0; i < np; i += 2) {
            yc = (AcX[i] + AcY[i]) / 2.0;
            yr = (AcY[i] - AcX[i]) / 2.0;
            xc = (AcX[i + 1] + AcY[i + 1]) / 2.0;
            xr = (AcY[i + 1] - AcX[i + 1]) / 2.0;

            d1 = fabs(yc - pl[0] - i_center(xc,xr,pl[2],pl[3]));
            d2 = fabs(yr - pl[1] - i_radius(xc,xr,pl[2],pl[3]));
            *rl += d1 * d1 + d2 * d2;
        }
/**
        printf("*rl=%g\n",*rl);          
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

    i_sub(bcl,bch,brl,brh,&sl,&sh);
    i_abs(sl,sh,&tl,&th);
    i_add(bcl,bch,brl,brh,&sl,&sh);
    i_abs(sl,sh,&ul,&uh);
    i_max(tl,th,ul,uh,&abl,&abh);



    for (i = 0; i < np; i += 2) {

        yc = (AcX[i] + AcY[i]) / 2.0;
        yr = (AcY[i] - AcX[i]) / 2.0;
        xc = (AcX[i + 1] + AcY[i + 1]) / 2.0;
        xr = (AcY[i + 1] - AcX[i + 1]) / 2.0;

        ax = dmax(fabs(xc - xr),fabs(xc + xr));

        i_mul(xr,xr,abl,abh,&ml1,&mh1);
        i_mul(brl,brh,ax,ax,&ml2,&mh2);

        i_max(ml1,mh1,ml2,mh2,&ml,&mh);

        i_abs(bcl,bch,&ul,&uh);
        i_mul(xr,xr,ul,uh,&sl,&sh);

        i_mul(brl,brh,fabs(xc),fabs(xc),&tl,&th);
        i_add(sl,sh,tl,th,&ml2,&mh2);
        i_max(ml,mh,ml2,mh2,&ml1,&mh1);

        i_add(arl,arh,ml1,mh1,&sl,&sh);
        i_sub(yr,yr,sl,sh,&ul,&uh);
        i_abs(ul,uh,&d2l,&d2h);
/*
        printf("d2l=%g %g\n",d2l,d2h);
*/
        i_mul(xc,xc,bcl,bch,&sl,&sh);

        i_abs(bcl,bch,&tl,&th);
        i_mul(xr,xr,tl,th,&ml1,&mh1);
        i_mul(brl,brh,fabs(xc),fabs(xc),&ml2,&mh2);
        i_min(ml1,mh1,ml2,mh2,&ml,&mh);
        i_mul(xr,xr,brl,brh,&ml2,&mh2);
        i_min(ml,mh,ml2,mh2,&ml1,&mh1);

        if (sl > 0.0) {
            i_add(sl,sh,ml1,mh1,&ul,&uh);
        }
        else if (sh < 0.0) {
            i_sub(sl,sh,ml1,mh1,&ul,&uh);
        }
        else {
            i_add(sl,sh,ml1,mh1,&tl,&th);
            i_sub(sl,sh,ml1,mh1,&vl,&vh);
            i_max(tl,th,vl,vh,&sl,&uh);
            i_min(tl,th,vl,vh,&ul,&sh);

/**
printf("tl=%g %g vl=%g %g sl=%g uh=%g ul=%g sh=%g\n",
            tl,th,vl,vh,sl,uh,ul,sh);
exit(0);
**/

        }
        i_add(acl,ach,ul,uh,&sl,&sh);
        i_sub(yc,yc,sl,sh,&tl,&th);
        i_abs(tl,th,&d1l,&d1h);
/*
        printf("d1l=%g %g\n",d1l,d1h);
*/

        i_square(d1l,d1h,&sl,&sh);
        i_square(d2l,d2h,&tl,&th);
        i_add(sl,sh,tl,th,&ul,&uh);
        i_add(*rl,*rh,ul,uh,&tl,&th);
        *rl = tl;
        *rh = th;

    }
/*
    printf("*rl=%g %g\n",*rl,*rh);             
**/         
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  i_center(ac,ar,bc,br)   Return center of a * b                          */

double i_center(double ac,double ar,double bc,double br)
{
    double c,d;             

    c = ac * bc;
    d = c + dsign(c) * dmin(ar * fabs(bc),dmin(br * fabs(ac),ar * br));
    return(d);
}

/* ------------------------------------------------------------------------ */
/*  i_radius(ac,ar,bc,br)   Return radius of a * b                          */

double i_radius(double ac,double ar,double bc,double br)
{
    double aa,ba,d;             

    aa = dmax(fabs(ac + ar),fabs(ac - ar));
    ba = dmax(fabs(bc + br),fabs(bc - br));
    d = dmax(ar * ba,dmax(br * aa,ar * fabs(bc) + br * fabs(ac)));
    return(d);
}

/* ------------------------------------------------------------------------ */
/*  i_bvar                                                                  */

void i_bvar(double xl,double xh,double bl,double bh,double *ul,double *uh)
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
            *ul = dmin(xl * bh,xh * bl);
            *uh = dmin(xl * bl,xh * bh);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  i_norm.   Return distance.                                              */

double i_norm(int opt,double xl,double xh,double yl,double yh)
{
    double d;             

    d = dmax(fabs(xl - yl),fabs(xh - yh));
    if (opt)
        d *= d;
    return(d);
}

/* ------------------------------------------------------------------------ */
/*  ivreg1_nrm.    Return distance.                                          */

double ivreg1_nrm(int opt,double xl,double xh,double yl,double yh,
    double al,double ah,double bl,double bh)
{
    double sl,sh,tl,th;

    i_mul(xl,xh,bl,bh,&sl,&sh);
    i_add(al,ah,sl,sh,&tl,&th);
    sl = dmax(fabs(yl - tl),fabs(yh - th));
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

int ivreg2(void)
{
    register int i;
    int err,r,rn,n,nw;
    double x,y,lsn,bc,br,ac,ar,f,d1,d2;
           
    err = -1;         
    if (check_cmd(2))
        return(-1);

    set_mldef();        /* set defaults for minimization */
    SLen = 1.0;         /* default step length */
    SRed = 0.5;         /* default reduction factor */
    TOLS = 1.0e-8;      /* minimum of step length */

    printf1("Regression with interval data. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 6,4,1))     /* get parameters */
        goto IVREGFin;

    if (PMOPT > 4)
        PMOPT = 4;

    if (PMNV != 4) {   
        printf1("Error: need four variables on right-hand side.\n");
        goto IVREGFin;
    }
    if (PMFmtF == 0)        /* default print format */
        pmfmt(14,8);

    printf1("\nNumber of data intervals: %d\n",NOC);
    if (NOC < 2) {
        printf1("Error: need at least two data intervals.\n");
        goto IVREGFin;
    }

    /* SET UP DATA ==================================================== */

    if (alloc_acx(NOC + 1))             /* X center */
        goto IVREGFin;
    if (alloc_acy(NOC + 1))             /* Y center */
        goto IVREGFin;
    if (alloc_acu(NOC + 1))             /* X radius */
        goto IVREGFin;
    if (alloc_acv(NOC + 1))             /* Y radius */
        goto IVREGFin;


    IVRXC = AcX;
    IVRXR = AcU;
    IVRYC = AcY;
    IVRYR = AcV;
    IVRXSC = 0.0;
    IVRXSCA = 0.0;
    IVRXSR = 0.0;
    IVRXSA = 0.0;
    IVRYSC = 0.0;
    IVRYSR = 0.0;

    nw = 0;
    for (i = 0; i < NOC; ++i) {
        x = get_data(PMVIdx[0],i);
        y = get_data(PMVIdx[1],i);
        IVRYC[i] = (x + y) / 2.0;
        IVRYR[i] = (y - x) / 2.0;
        x = get_data(PMVIdx[2],i);
        y = get_data(PMVIdx[3],i);
        IVRXC[i] = (x + y) / 2.0;
        IVRXR[i] = (y - x) / 2.0;

        if (IVRXR[i] < -EPSI || IVRYR[i] < -EPSI) {
            printf1("Error: no valid interval in record %d.\n",i + 1);
            goto IVREGFin;
        }
        if (fabs(IVRXC[i]) < IVRXR[i] - EPSI)
            nw++;

        IVRXSC += IVRXC[i];
        IVRXSCA += fabs(IVRXC[i]);
        IVRXSR += IVRXR[i];
        IVRXSA += fabs(IVRXC[i]) + IVRXR[i];
        IVRYSC += IVRYC[i];
        IVRYSR += IVRYR[i];
    }
    IVRXSC /= (double)NOC;
    IVRXSCA /= (double)NOC;
    IVRXSR /= (double)NOC;
    IVRXSA /= (double)NOC;
    IVRYSC /= (double)NOC;
    IVRYSR /= (double)NOC;
      
/** 
    printf("data\n"); 
    for (i = 0; i < NOC; ++i) {
        printf("%3d %f %f %f %f\n",i,IVRXC[i],IVRXR[i],IVRYC[i],IVRYR[i] );
    }
    printf("ivrysc=%g ivrysr=%g\n",IVRYSC,IVRYSR);
    printf("ivrxsc=%g ivrxsr=%g\n",IVRXSC,IVRXSR);
**/    

    /* CONTOUR PLOT =================================================== */

    if (PMOPT == 4) {

        if (PMNNFlg == 0 || PMNN1 < 1 || PMNN2 < 1)
            PMNN1 = PMNN2 = 10;

        printf1("\nContour plot with %d,%d grid.\n",PMNN1,PMNN2);
        if (check_ps(2))
            goto IVREGFin;

        if (PMN > 0) {
            if (PMAFlg == 0) {
                printf1("Error: need a parameter (beta center and radius).\n");
                goto IVREGFin;
            }
            if (alloc_actmp(PMN + 2))                               
                goto IVREGFin;

            if (PMA2 < 0.0)
                PMA2 = 0.0;

            printf1("Level generation with beta=%g,%g\n",PMA1,PMA2);
            if (PMA1 <= PA1[0] - EPSI || PMA1 >= PA2[0] + EPSI ||
                PMA2 <= PA1[1] - EPSI || PMA2 >= PA2[1] + EPSI) {  
                printf1("Error: beta not inside grid.\n");
                goto IVREGFin;
            }
            d1 = (PMA1 - PA1[0]) / (double)PMN;
            d2 = (PA2[0] - PMA1) / (double)PMN;
            br = PMA2;
            i  = 0;

            if (d1 > d2) {
                bc = PA1[0] + d1 / 2.0;
                while (bc <= PMA1 + EPSI1) {
                    f = ivreg2f(bc,br);
                    printf1("Level %2d  bc=%g Function value: %20.12f\n",i + 1,bc,f);
                    AcTmp[i++] = f;
                    if (i == PMN)
                        break;
                    bc += d1;
                }
            }   
            else {
                bc = PA2[0] - d2 / 2.0;
                while (bc >= PMA1 - EPSI1) {
                    f = ivreg2f(bc,br);
                    printf1("Level %2d  bc=%g Function value: %20.12f\n",i + 1,bc,f);
                    AcTmp[i++] = f;
                    if (i == PMN)
                        break;
                    bc -= d2;
                }
            }   
        }
        else {
            if (PMNTP < 1) {
                printf1("Error: need level definition with x parameter.\n");
                goto IVREGFin;
            }
        }        
        /***
        for (i = 0; i < PMN; ++i)
            printf("%d %g\n",i,AcTmp[i]);
        ***/

        if (ivreg2p(PMN,AcTmp))
            goto IVREGFin;
    
        goto IVREGFin;
    }
           
    /* OLS REGRESSION ================================================= */

    printf1("\nOLS regression with center values.\n");

    if (alloc_acw(3 * NOC + 1))     /* data for lsei */
        goto IVREGFin;
    if (alloc_acz(3))               /* parameter */
        goto IVREGFin;

    for (i = 0; i < NOC; ++i) {
        AcW[i * 3 + 1] = 1.0;
        AcW[i * 3 + 2] = IVRXC[i];
        AcW[i * 3 + 3] = IVRYC[i];
    }
    r = lsei(AcW,0,NOC,0,2,AcZ,0,&x,&lsn,&rn,&n);
    if (r) {
        if (r == -4)  
            p_err(-2,1);
        else  
            printf1("LSEI error return: %d\n",err);
        goto IVREGFin;  
    }
    printf1("Rank of least squares matrix: %d\n",rn);
    printf1("Norm of least squares residuals:");
    printf1(PMFmtS,lsn);
    printf1("\nParameter:");
    printf1(PMFmtS,AcZ[1]);
    printf1(PMFmtS,AcZ[2]);
    newline();

    /* DIRECT CALCULATION ============================================= */

    ivreg2d();

    if (nw > 0) {
        printf1("Warning: Cases with x-radius greater than x-center: %d\n\n",nw);
    }

    /* FUNCTION MINIMIZATIOBN WITH DIRECT SEARCH  ===================== */

    if (PMOPT == 2) {

    printf1("Direct search minimization.\n");
    printf1("Step length: %g\n",SLen);
    printf1("Reduction factor: %g\n",SRed);
    printf1("Minimum of step length: %g\n\n",TOLS);

    bc = AcZ[2];
    br = 0.0;
    AcZ[1] = AcZ[2] = 0.0;

    if (get_dsv(2,AcZ,0,AcZ,AcZ,0))   /* get starting values */ 
        goto IVREGFin;
       
    if (DSVFlg) {
        bc = AcZ[1];
        if (AcZ[2] >= 0.0)
            br = AcZ[2];
    }
    printf1("Starting values:");
    printf1(PMFmtS,bc);
    printf1(PMFmtS,br);
    newline();

        f = ivreg2f(bc,br);

        printf1("Function value: ");
        printf1(PMFmtS,f);
        newline();
        
        f = ivreg2m(&bc,&br);

    }
    else if (PMOPT == 3) {

        NParm = 2;                  /* number of parameters */
        MINA = 2;                   /* always this algorithm */
        PMDOPT = 0;                 /* without derivatives */
        CCTyp = 0;                  /* no covariance matrix */
        set_mlopt();                /* adjust options */         

        if (ml_init(1,1,0))         /* init function minimization */
            goto IVREGFin;

        Par[1] = AcZ[2];
        Par[2] = 0.0;

f = ivreg2f(Par[1],Par[2]);
printf("f=%g\n",f);


        if (ffmin(IVRMOD1,1,0,1)) {                               
            printf1("ERROR\n");       
            exit(0);
        }
        prn_mlres(1);   
        printf("FMIN=%g\n",FMin);
        printf("par=%20.12f %20.12f\n",Par[1],Par[2]);

exit(0);



    }
    ar = ac = 0.0;
    for (i = 0; i < NOC; ++i) {
        ac += IVRYC[i] - i_center(IVRXC[i],IVRXR[i],bc,br);
        ar += IVRYR[i] - i_radius(IVRXC[i],IVRXR[i],bc,br);
    }
    ac /= (double)NOC;
    ar /= (double)NOC;
          
    printf1("\nBest function value:  ");
    printf1(PMFmtS,f);
    printf1("\nAlpha (center-radius):");
    printf1(PMFmtS,ac);
    printf1(PMFmtS,ar);
    printf1("\nBeta  (center-radius):");
    printf1(PMFmtS,bc);
    printf1(PMFmtS,br);
    newline();

    newline();
    err = 0;

IVREGFin:
    con_free();
    ml_init(0,0,0);      
    p_clean();
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

double ivreg2m(double *bc,double *br)
{
    int k;
    double sl,rf,f,f0,f1,f2,f3,f4,f5,f6,f7,f8;
    double bc0,bc1,bc2,bc3,bc4,bc5,bc6,bc7,bc8;
    double br0,br1,br2,br3,br4,br5,br6,br7,br8;

    sl = SLen;
    rf = SRed;
/**
    printf("hier sl=%g rf=%g tols=%g bc=%g br=%g\n",sl,rf,TOLS,*bc,*br  );
**/

    f  = ivreg2f(*bc,*br);

    if (SILENTFlg < 2) {
        printfe("\nBeta center          Beta radius         Function value      Step size\n");
        printfe("%19.12e %19.12e %19.12e %19.12e\n",*bc,*br,f,sl);
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

        f1 = ivreg2f(bc1,br1);
        f2 = ivreg2f(bc2,br2);
        f3 = ivreg2f(bc3,br3);
        f4 = ivreg2f(bc4,br4);
        f5 = ivreg2f(bc5,br5);
/**
        printf("\nslen=%14.12f\n",sl);
        printf("bc=%16.12f br=%16.12f f =%19.14f\n",*bc,*br,f);
        printf("bc=%16.12f br=%16.12f f1=%19.14f\n",bc1,br1,f1);
        printf("bc=%16.12f br=%16.12f f2=%19.14f\n",bc2,br2,f2);
        printf("bc=%16.12f br=%16.12f f3=%19.14f\n",bc3,br3,f3);
        printf("bc=%16.12f br=%16.12f f4=%19.14f\n",bc4,br4,f4);
        printf("bc=%16.12f br=%16.12f f5=%19.14f\n",bc5,br5,f5);
**/
        k = 1;
        f0 = f1;
        bc0 = bc1;
        br0 = br1;
        if (f2 < f0) {
            k = 2;
            f0 = f2;
            bc0 = bc2;
            br0 = br2;
        }
        if (f3 < f0) {
            k = 3;
            f0 = f3;
            bc0 = bc3;
            br0 = br3;
        }
        if (f4 < f0) {
            k = 4;
            f0 = f4;
            bc0 = bc4;
            br0 = br4;
        }
        if (f5 < f0) {
            k = 5;
            f0 = f5;
            bc0 = bc5;
            br0 = br5;
        }
/**
        printf("k=%d f0=%19.14f\n",k,f0);
**/
        if (br6 >= -EPSI) {
            f6 = ivreg2f(bc6,br6);
            f8 = ivreg2f(bc8,br8);
  
/**
            printf("bc=%16.12f br=%16.12f f6=%19.14f\n",bc6,br6,f6);
            printf("bc=%16.12f br=%16.12f f8=%19.14f\n",bc8,br8,f8);
**/
            if (f6 < f0) {
                k = 6;
                f0 = f6;
                bc0 = bc6;
                br0 = br6;

            }
            if (f8 < f0) {
                k = 8;
                f0 = f8;
                bc0 = bc8;
                br0 = br8;
            }
        }
        if (br7 >= -EPSI) {
            f7 = ivreg2f(bc7,br7);
/*
            printf("bc=%16.12f br=%16.12f f7=%19.14f\n",bc7,br7,f7);
*/
            if (f7 < f0) {
                k = 7;
                f0 = f7;
                bc0 = bc7;
                br0 = br7;
            }
        }
        if (sl <= TOLS) {
            break;
        }
        if (f0 >= f - EPSI) {
            sl *= rf; 
        }
        else {
            f = f0;
            *bc = bc0;
            *br = br0;
            if (SILENTFlg < 2)  
                printfe("%19.12e %19.12e %19.12e %19.12e\n",*bc,*br,f,sl);
        }
    }
    if (SILENTFlg < 2)  
        printfe("%19.12e %19.12e %19.12e %19.12e\n",*bc,*br,f,sl);

    return(f);
}

/* ------------------------------------------------------------------------ */
/*  ivreg2f(bc,br)      Return function value (using data as set up by      */
/*                      ivreg2().)                                          */

double ivreg2f(double bc,double br) 
{
    register int i;                    
    double xc,xr,yc,yr,d1,d2,xbc,xbr,xbcs,xbrs;
         
    d1 = d2 = xbcs = xbrs = 0.0;

    for (i = 0; i < NOC; ++i) {
        yc = IVRYC[i];
        yr = IVRYR[i];
        xc = IVRXC[i];
        xr = IVRXR[i];

        xbc = i_center(xc,xr,bc,br);
        xbr = i_radius(xc,xr,bc,br);
             
        d1 += xbc * (xbc - 2.0 * (yc - IVRYSC));
        d2 += xbr * (xbr - 2.0 * (yr - IVRYSR));

        xbcs += xbc;
        xbrs += xbr;
    }
    d1 -= xbcs * xbcs / (double)NOC;
    d2 -= xbrs * xbrs / (double)NOC;
    return(d1 + d2);
}

/*--------------------------------------------------------------------------*/
/*  ivreg2p(nlev,flev)   Contour plot for ivreg2(). If nlev > 0 the array   */
/*                       flev contains level.                               */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int ivreg2p(int nlev,double *flev)
{
    register int i,j;
    int nn,err,pcbita,pcfva;
    double bc,br,*fl;
           
    pcbita = pcfva = 0;
    err = -1;

    if (nlev > 0)  
        fl = flev;
    else {
        nlev = PMNTP;
        fl = PMTP;
    }


    if (PMNNFlg == 0 || PMNN1 < 1 || PMNN2 < 1)
        PMNN1 = PMNN2 = 10;
      
    PCDX = UXLen / (double) PMNN1;
    PCDY = UYLen / (double) PMNN2;
    PCNX = PMNN1 + 1;
    PCNY = PMNN2 + 1;

    if (!(PCFV = (float *) calloc(PCNX * PCNY + 1,sizeof(float)))) {
        p_err(-2,1);
        goto IVCPFin;
    }
    memrq(PCNX * PCNY + 1,sizeof(float));
    pcfva = 1;

    nn = 2 * PCNX * PCNY * nlev / 8 + 1;

    if (!(PCBitM = (char *)calloc(nn,sizeof(char)))) {
        p_err(-2,1);
        goto IVCPFin;
    }
    memrq(nn,sizeof(char));
    pcbita = 1;

    /* get function values in PCFV */
      
    for (i = 1; i <= PCNY; ++i) {
        for (j = 1; j <= PCNX; ++j) {
     
            bc = PA1[0] + (double)(j - 1) * PCDX;
            br = PA1[1] + (double)(i - 1) * PCDY;
            PCFV[(j - 1) * PCNX + i] = (float)ivreg2f(bc,br);
        }
    }
    fprintf(PSFd,"\n%%#%d: ivreg contour plot\n",++PSONUM);

    cont_plot(nlev,fl,0);
                    
    err = 0;

IVCPFin:
    if (pcfva) {
        free((char *)PCFV);
        memrq(-PCNX * PCNY - 1,sizeof(float));
    }
    if (pcbita) {
        free((char *)PCBitM);
        memrq(-nn,sizeof(char));
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

int ivr1_fn(void)
{
    register int i;
    int err;
    double bc,br,xc,xr,yc,yr;

    bc = TPar[1];
    br = fabs(TPar[2]);
    FTmp = 0.0;

    for (i = 0; i < NOC; ++i) {
        yc = IVRYC[i];
        yr = IVRYR[i];
        xc = IVRXC[i];
        xr = IVRXR[i];

        FTmp += (xc * bc + xr * br) * (bc * (xc - IVRXSC) +
                                       br * (xr - IVRXSR) - 
                                      2.0 * (yc - IVRYSC)) +
                (xr * bc + xc * br) * (bc * (xr - IVRXSR) +
                                       br * (xc - IVRXSC) - 
                                      2.0 * (yr - IVRYSR)); 
    }
    err = checkov();        
    if (err) {
        printf1("\nError in function evaluation.\n");
        printf1("Numerical overflow.\n");
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ivreg2dp(n,na,bc,br)                                                    */

void ivreg2dp(int n,int na,double bc,double br) 
{
    register int i;
    double f,ac,ar; 

    printf1("%4d ",n);
    if (na < 0) {
        printf1("\n");
        return;
    }
    else if (na > 0)
        printf1("* ");
    else
        printf1("  ");

    prnchar(' ',11 - PMFmt1,0);
    printf1(PMFmtS,bc);
    prnchar(' ',11 - PMFmt1,0);
    printf1(PMFmtS,br);

    if (na > 0) {
        ar = ac = 0.0;
        for (i = 0; i < NOC; ++i) {
            ac += IVRYC[i] - i_center(IVRXC[i],IVRXR[i],bc,br);
            ar += IVRYR[i] - i_radius(IVRXC[i],IVRXR[i],bc,br);
        }
        ac /= (double)NOC;
        ar /= (double)NOC;

        prnchar(' ',12 - PMFmt1,0);
        printf1(PMFmtS,ac);
        prnchar(' ',12 - PMFmt1,0);
        printf1(PMFmtS,ar);
    }
    else
        prnchar(' ',2 * imax(13,PMFmt1 + 1),0);

    f = ivreg2f(bc,br);
    prnchar(' ',16 - PMFmt1,0);
    printf1(PMFmtS,f);
    newline();             
}

/* ------------------------------------------------------------------------ */
/*  ivreg2d()                                                               */

int ivreg2d(void)
{
    register int i;
    int na;
    double tmp,a,b,c,d,e,br,bc,xc,xr,xa,xca,yc,yr,vc,vr,va,vca,wc,wr;
    double vzs,vcs,vrs,vas,vcas,vcar,vcara,vss,wzs,wcs,wrs,urs,wcar,wcars;

    vzs = vcs = vrs = vas = vcas = wzs = wcs = wrs = urs = 0.0;
    wcars = wcar = vss = vcar = vcara = 0.0;

    for (i = 0; i < NOC; ++i) {
        yc = IVRYC[i];
        yr = IVRYR[i];
        xc = IVRXC[i];
        xr = IVRXR[i];
        xca = fabs(xc);
        xa = xca + xr;

        vc = (xc - IVRXSC);
        vca = (xca - IVRXSCA);
        vr = (xr - IVRXSR);
        va = (xa - IVRXSA);
        wc = (yc - IVRYSC);
        wr = (yr - IVRYSR);

        vcs += xc * vc;
        vrs += xr * vr;
        vas += xa * va;
        vcas += xca * vca;
        vcar += xca * vr;
        vcara += xr * vca;
        vzs += (xc + dsign(xc) * xr) * (vc + dsign(xc) * vr);
        vss += dsign(xc) * (xc * vr + xr * vc);

        wcs += xc * wc;
        wrs += xr * wr;
        wzs += (xc + dsign(xc) * xr) * wc;                     
        wcar += xca * wr;
        wcars += dsign(xc) * xr * wc;
        urs += xa * wr;


    }
    printf1("\nDomain ");
    prnchar(' ',PMFmt1 - 11,0);
    printf1("Beta center ");
    prnchar(' ',PMFmt1 - 11,0);
    printf1("Beta radius ");
    prnchar(' ',PMFmt1 - 12,0);
    printf1("Alpha center ");
    prnchar(' ',PMFmt1 - 12,0);
    printf1("Alpha radius   ");
    prnchar(' ',PMFmt1 - 14,0);
    printf1("Function value\n");

    /* domain 1 */

    na = 0;
    bc = br = 0.0;
    tmp = vcs + vrs;
    if (tmp == 0.0)
        na = -1;
    else {
        bc = (wcs + wrs) / tmp;
        if (bc >= -EPSI1)
            na = 1;
    }
    ivreg2dp(1,na,bc,br);

    /* domain 2 */

    na = 0;
    bc = br = 0.0;
    tmp = vcs + vrs;
    if (tmp == 0.0)
        na = -1;
    else {
        bc = (wcs - wrs) / tmp;
        if (bc <= EPSI1)
            na = 1;
    }
    ivreg2dp(2,na,bc,br);

    /* domain 3 */

    na = 0;
    bc = br = 0.0;
    if (vas == 0.0 || vzs == 0.0)
        na = -1;
    else {
        bc = wzs / vzs;
        br = urs / vas;
        if (br >= -EPSI1 && fabs(bc) <= br + EPSI1)
            na = 1;
    }
    ivreg2dp(3,na,bc,br);

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
        if (br >= -EPSI1 && bc >= -EPSI1 && bc >= br - EPSI1)
            na = 1;
    }
    ivreg2dp(4,na,bc,br);

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
        if (br >= -EPSI1 && bc <= EPSI1 && -bc >= br - EPSI1)
            na = 1;
    }
    ivreg2dp(5,na,bc,br);




    newline();
    return(0);
}


