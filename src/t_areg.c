/****************************************************************************/
/*  t_areg                                                                  */
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
#include "t_sort.h"
#include "t_ds.h"
#include "t_smo.h"     
#include "t_gf.h"
#include "t_ml.h"

/*  functions in t_areg.c */

int npreg(void);
int npreg1(int opt,int kopt,int n,double *x,double *y,int nx,double *xx,double d);
void npreg2(int opt,int n,int ny,int *y);
double mkern(int opt,int n,double d,double a,double *x,double *y);
double mkern1(int opt,int n,double d,double a,double *x);
int lowess(double *x,double *y,int n,double f,int nsteps,double delta,
    double *ys,double *rw,double *res);
int lowest(double *x,double *y,int n,double xs,double *ys,int nleft,
    int nright,double *w,int userw,double *rw);
int rmidmean(int n,double *x,double *y,int r,double *xm,double *ym,
    float *ym1,float *ym2);
double midmean(int n,double *x,int opt,int *err);
int smidmean(int n,double *x,double *sm1,double *sm2);
int inpreg(void);
double imkern(int opt,int n,double d,double a,double *x,double *y);


/* ------------------------------------------------------------------------ */
/*  npreg()         non-parametric regression                               */
/*                                                                          */  
/*                  npreg(                                                  */
/*                      opt=...,    1 smoothed means                        */
/*                                  2 smoothed quantiles                    */
/*                                  3 smoothes frequencies                  */
/*                                  4 lowess                                */
/*                                  5 midmeans                              */
/*                      k=...,      kernel, def. 1                          */
/*                                  1 : uniform                             */
/*                                  2 : triangle                            */
/*                                  3 : quartic                             */
/*                                  4 : Epanechnikov                        */
/*                      x=...,      list of points for X axis               */
/*                      d=...,      widths of interval, def. 1              */
/*                      df=...,     output file                             */
/*                      fmt=...,    print format, def. 10.4                 */
/*                      dtda=...,   TDA description                         */
/*                      sel=...,    case selection                          */
/*                  ) = Y,X;                                                */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int npreg(void)
{
    register int i;
    int err,nn,ix,iy,nrec;

    err = -1;         
    if (check_cmd(1))
        return(-1);
    /*  prnchar('-',LLEN,1); */
    printf1("Nonparametric regression. Current memory: %d bytes.\n",MemReq);
         
    if (parm(CmdBuf + 5,4,1))     /* get parameters */
        goto NPRFin;

    if (PMNV != 2) {            /* need two variables */
        p_err(-1,1);
        goto NPRFin;
    }
                         /* PMNTP = number of arguments in PMTP[] (t=...) */
    if (PMNTP < 1) {
        printf1("Error: need points defined with x parameter.\n");
        goto NPRFin;
    }
    if (PMD < EPSI1)
        PMD = 1.0;

    if (!PMF1Def) {             /* need an output file */
        p_err(-27,1);
        goto NPRFin;
    }
    if (PMOPT < 1 || PMOPT > 5)
        PMOPT = 1;

    iy = PMVIdx[0];
    ix = PMVIdx[1];

    if (alloc_acx(NOC + 1))
        goto NPRFin;
    if (alloc_acy(NOC + 1))
        goto NPRFin;

    nn = getd2(AcX,AcY,ix,iy,1);
    if (nn == 0)  
        goto NPRFin;

    if (PMFmtF == 0)  
        pmfmt(10,4);

    printf1("Number of cases: %d\nMethod: ",nn);

    nrec = 0;
    if (PMOPT <= 3) {           /* means */
        printf1("smoothed ");
        if (PMOPT == 1)  
            printf1("means");
        else if (PMOPT == 2)
            printf1("quantiles");
        else            
            printf1("frequencies");

        printf1(". Interval width: %lg\n",PMD);
        if (PMOPT == 1) {
            printf1("Kernel: ");
            switch (PMK) {
                case 2: printf1("triangle\n"); break;
                case 3: printf1("quartic\n"); break;
                case 4: printf1("Epanechnikov\n"); break;
                default: PMK = 1; printf1("uniform\n"); break;
            }
        }
        nrec = npreg1(PMOPT,PMK,nn,AcX,AcY,PMNTP,PMTP,PMD);
    }
    else if (PMOPT == 4) {           /* lowess */

        if (alloc_acu(nn + 1))
            goto NPRFin;
        if (alloc_acv(nn + 1))
            goto NPRFin;
        if (alloc_actmp(nn + 1))
            goto NPRFin;

        if (PMNS < 0)               /* number of steps */
            PMNS = 2;
        if (PMSIG <= 0.0 || PMSIG >= 1.0)
            PMSIG = 0.5;
        if (PMD < 0.0)
            PMD = 0.0;
        printf1("Lowess: sig=%lg ns=%d d=%lg\n",PMSIG,PMNS,PMD);

        if (lowess(AcX - 1,AcY - 1,nn,PMSIG,PMNS,PMD,AcU,AcV,AcTmp))
            goto NPRFin;

        for (i = 1; i <= nn; ++i) {
            fprintf(PMF1d,"%6d ",i);
            fprintf(PMF1d,PMFmtS,AcX[i - 1]);
            fprintf(PMF1d,PMFmtS,AcY[i - 1]);
            fprintf(PMF1d,PMFmtS,AcU[i]);
            fprintf(PMF1d,"\n");
        }
        printf1("%d records written to: %s\n",nrec,PMF1dName);
    }
    else if (PMOPT == 5) {           /* midmeans */

        if (PMR < 2)
            PMR = 2;
        if (PMR >= nn) {
            printf1("Error: r must be less than number of cases.\n");
            goto NPRFin;
        }
        printf1("Midmeans: r=%d\n",PMR);
        /**
        for (i = 0; i < nn; ++i) {
            nneighbor(nn,AcX,PMR,i,&k1,&k2);
            AxU[i] = midmean(PMR,AcX + k1,0,&err);

            for (j = 0; j < PMR; ++j)
                AcTmp[j] = AcY[j + k1];

            AcV[i] = midmean(PMR,AcTmp,1,&err);
        }
        **/
    }
    err = 0;
NPRFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  npreg1(opt,kopt,n,x,y,nx,xx,d)                                          */
/*                                                                          */
/*  Nonparametric regression:                                               */
/*      opt = 1 : smoothed means.                                           */
/*      opt = 2 : smoothed quantiles.                                       */
/*      opt = 3 : smoothed frequencies.                                     */
/*                                                                          */
/*  Kernel defined by                                                       */
/*  kopt = 1 : uniform                                                      */
/*         2 : triangle                                                     */
/*         3 : quartic                                                      */
/*         4 : Epanechnikov                                                 */
/*                                                                          */
/*  n = number of data points                                               */
/*  data: x[i],y[i], i = 0,...,n-1                                          */
/*  xx[j], j = 0,...,nx-1 list of points on X axis.                         */
/*  d = interval width                                                      */
/*                                                                          */  
/*  return -1 if error (insuff. memory), otherwise number of records        */
/*  written to output file.                                                 */

int npreg1(int opt,int kopt,int n,double *x,double *y,int nx,double *xx,double d)
{
    register int i,j,k,m;
    int nm,nqx,nrec,nycat;
    double w,a,b,xm,ym,ys;
    double qx[11];
        
    nycat = 0;

    if (sortd2(n,x,y))     /* sort simultaneously according to x */
        return(-1);       

    if (alloc_acu(n + 1))
        return(-1);  
    if (alloc_acv(n + 1))
        return(-1);         
    if (alloc_acw(n + 1))
        return(-1);         

    if (opt == 2) {
        nqx = 11;
        qx[ 0] = 0.1 ;
        qx[ 1] = 0.2 ;
        qx[ 2] = 0.25;
        qx[ 3] = 0.3 ;
        qx[ 4] = 0.4 ;
        qx[ 5] = 0.5 ;
        qx[ 6] = 0.6 ;
        qx[ 7] = 0.7 ;
        qx[ 8] = 0.75;
        qx[ 9] = 0.8 ;
        qx[10] = 0.9 ;
    }
    else if (opt == 3) {    /* get categories of Y */

        if (alloc_acn(n + 1))
            return(-1);  

        for (i = 0; i < n; ++i) 
            AcU[i] = y[i];

        if (sortd(n,AcU,0))   
            return(-1);       

        nycat = 1;
        AcN[0] = (int)AcU[0];
        i = 0;
        while (++i < n) {
            if (AcU[i] != AcU[i - 1])
                AcN[nycat++] = (int)AcU[i];
        }
        printf1("Number of categories in dependent variable: %d\nCategories:",nycat);
        for (i = 0; i < nycat; ++i)
            printf1(" %d",AcN[i]);
        newline();     
         
        if (alloc_acm(nycat))
            return(-1);  
    }

    w = d / 2.0;
    for (i = 0; i < n; ++i)
        AcW[i] = 1.0;

    nrec = 0;
    for (j = 0; j < nx; ++j) {
  
        a = xx[j] - w;
        b = xx[j] + w;
        m = 0;
        for (i = 0; i < n; ++i) {
            if (x[i] < a)
                continue;
            if (x[i] > b)
                break;

            AcU[m] = x[i];
            AcV[m++] = y[i];
        }
        if (m == 0)
            continue;

        xm = xmean(m,AcU,AcW);

        fprintf(PMF1d,"%6d ",j + 1);
        fprintf(PMF1d,PMFmtS,xx[j]);
        fprintf(PMF1d,"%6d ",m);
        fprintf(PMF1d,PMFmtS,xm);   /* x mean */

        if (opt == 1) {
            ym = xstd(m,AcV,AcW,&ys);

            if (kopt > 1)
                ym = mkern(kopt,m,d,xx[j],AcU,AcV);

            fprintf(PMF1d,PMFmtS,ym);   /* y mean */
            fprintf(PMF1d,PMFmtS,ys);   /* y standard deviation */
        }
        else if (opt == 2) {
            if (m > 1) {
                if (sortd(m,AcV,0))
                    return(-1);    
            }
            for (i = 0; i < nqx; ++i) {
                ym = quantf(m,AcV,qx[i]);
                fprintf(PMF1d,PMFmtS,ym); 
            }
        }
        else if (opt == 3) {
            for (k = 0; k < nycat; ++k)
                AcM[k] = 0;
            for (i = 0; i < m; ++i) {
                nm = (int)AcV[i];
                for (k = 0; k < nycat; ++k) {
                    if (AcN[k] == nm) {
                        AcM[k] += 1;
                        break;
                    }
                }
            }
            for (k = 0; k < nycat; ++k)  
                fprintf(PMF1d,PMFmtS,(double)AcM[k] / (double)m);
        }
        fprintf(PMF1d,"\n");
        nrec++;
    }
    printf1("%d records written to: %s\n",nrec,PMF1dName);
    if (PMTDAFDef)  
        npreg2(opt,nrec,nycat,AcN);

    return(nrec);
}

/* ------------------------------------------------------------------------ */
/*  npreg2(opt,n,ny,y)                                                      */
/*                                                                          */
/*                  write TDA description file, n = # of cases              */
/*                  opt = 1 : smoothed means                                */
/*                        2 : smoothed quantiles                            */
/*                        3 : smoothed frquencies. In this case ny is the   */
/*                            number of Y categories in y[].                */
/*                        4 : lowess                                        */
/*                        5 : midmeans                                      */
/*                                                                          */

void npreg2(int opt,int n,int ny,int *y)
{
    register int i,k;
    int q[11];

    if (!PMTDAFDef)  
        return;

    fprintf(PMTDAFd,"# data written by npreg command.\n");
    fprintf(PMTDAFd,"nvar(\n");
    fprintf(PMTDAFd,"  dfile = %s,\n",PMF1dName);
    fprintf(PMTDAFd,"  noc = %d,\n",n);
    k = 0;
    fprintf(PMTDAFd,"  RECN [ 6.0] = c%-2d, # record number\n",++k);
    fprintf(PMTDAFd,"  X    [%d.%d] = c%-2d, # x value\n",PMFmt1,PMFmt2,++k);
    fprintf(PMTDAFd,"  NX   [ 6.0] = c%-2d, # number of cases\n",++k);
    fprintf(PMTDAFd,"  XM   [%d.%d] = c%-2d, # mean of X\n",PMFmt1,PMFmt2,++k);

    switch (opt) {
        case 1:
            fprintf(PMTDAFd,"  YM   [%d.%d] = c%-2d, # mean of Y\n",PMFmt1,PMFmt2,++k);
            fprintf(PMTDAFd,"  YSD  [%d.%d] = c%-2d, # standard deviation of Y\n",PMFmt1,PMFmt2,++k);
            break;

        case 2:
            q[ 0] = 10;
            q[ 1] = 20;
            q[ 2] = 25;
            q[ 3] = 30;
            q[ 4] = 40;
            q[ 5] = 50;
            q[ 6] = 60;
            q[ 7] = 70;
            q[ 8] = 75;
            q[ 9] = 80;
            q[10] = 90;

            for (i = 0; i < 11; ++i)  
                fprintf(PMTDAFd,"  Q%2d  [%d.%d] = c%-2d, # .%2d quantile of Y\n",q[i],PMFmt1,PMFmt2,++k,q[i]);
            break;

        case 3:
            for (i = 0; i < ny; ++i)  
                fprintf(PMTDAFd,"  Y%-4d[%d.%d] = c%-2d, # frequency of Y=%d\n",y[i],PMFmt1,PMFmt2,++k,y[i]);
            break;

        default:    break;
    }
    fprintf(PMTDAFd,");\n");
    printf1("TDA description written to: %s\n",PMTDAFName);
}

/* ------------------------------------------------------------------------ */
/*  mkern(opt,n,d,a,x,y)                                                    */
/*                                                                          */
/*  calculate mean of y at a, depending on type of kernel with bandwidth d. */
/*  opt = 2 :   triangle                                                    */
/*        3     quartic                                                     */
/*        4     Epanechnikov                                                */
/*                                                                          */
/*  return mean of y.                                                       */

double mkern(int opt,int n,double d,double a,double *x,double *y)
{
    register int i;
    double ym,xm,dd,w;

    w = xm = ym = 0.0;
    dd = d / 2.0;

    for (i = 0; i < n; ++i) {
        if (opt == 2)
            w = (1.0 - fabs(a - x[i]) / dd) / dd;
        else if (opt == 3) {
            w = (a - x[i]) / dd;
            w = 1.0 - w * w;  
            w = 15.0 * w * w / (8.0 * d);
        }
        else if (opt == 4) {
            w = (a - x[i]) / dd;
            w = 1.0 - w * w;  
            w = 3.0 * w / (2.0 * d);
        }
        xm += w;
        ym += y[i] * w;
    }
    if (xm > 0.0)
        ym /= xm;
    return(ym);
}

/* ------------------------------------------------------------------------ */
/*  mkern1(opt,n,d,a,x)                                                     */
/*  ##                                                                      */
/*  calculate mean of x at a, depending on type of kernel with bandwidth d. */
/*  opt = 1 :   uniform                                                     */
/*        2     triangle                                                    */
/*        3     quartic                                                     */
/*        4     Epanechnikov                                                */
/*                                                                          */
/*  return mean of x.                                                       */

double mkern1(int opt,int n,double d,double a,double *x)
{
    register int i;
    double xm,dd,w;

    w = xm = 0.0;
    dd = d / 2.0;

    for (i = 0; i < n; ++i) {
        if (opt == 1)
            w = 1.0 / d;
        else if (opt == 2)
            w = (1.0 - fabs(a - x[i]) / dd) / dd;
        else if (opt == 3) {
            w = (a - x[i]) / dd;
            w = 1.0 - w * w;  
            w = 15.0 * w * w / (8.0 * d);
        }
        else if (opt == 4) {
            w = (a - x[i]) / dd;
            w = 1.0 - w * w;  
            w = 3.0 * w / (2.0 * d);
        }
        xm += w;
    }
    return(xm);
}

/* ------------------------------------------------------------------------ */
/*                                                                          */
/*        lowess (x,y,n,f,nsteps,delta,ys,rw,res)                           */
/*                                                                          */
/*        LOWESS computes the smooth of a scatterplot of Y  against  X      */
/*        using  robust  locally  weighted regression.  Fitted values,      */
/*        YS, are computed at each of the  values  of  the  horizontal      */
/*        axis in X.                                                        */
/*                                                                          */
/*              X = Input; abscissas of the points on the                   */
/*                  scatterplot; the values in X must be ordered            */
/*                  from smallest to largest.                               */
/*              Y = Input; ordinates of the points on the                   */
/*                  scatterplot.                                            */
/*              N = Input; dimension of X,Y,YS,RW, and RES.                 */
/*              F = Input; specifies the amount of smoothing; F is          */
/*                  the fraction of points used to compute each             */
/*                  fitted value; as F increases the smoothed values        */
/*                  become smoother; choosing F in the range .2 to          */
/*                  idea which value to use, try F = .5.                    */
/*         NSTEPS = Input; the number of iterations in the robust           */
/*                  fit; if NSTEPS = 0, the nonrobust fit is                */
/*                  returned; setting NSTEPS equal to 2 should serve        */
/*                  most purposes.                                          */
/*          DELTA = input; nonnegative parameter which may be used          */
/*                  to save computations; if N is less than 100, set        */
/*                  DELTA equal to 0.0; if N is greater than 100 you        */
/*                  should find out how DELTA works by reading the          */
/*                  additional instructions section.                        */
/*             YS = Output; fitted values; YS(I) is the fitted value        */
/*                  at X(I); to summarize the scatterplot, YS(I)            */
/*                  should be plotted against X(I).                         */
/*             RW = Output; robustness weights; RW(I) is the weight         */
/*                  given to the point (X(I),Y(I)); if NSTEPS = 0,          */
/*                  RW is not used.                                         */
/*            RES = Output; residuals; RES(I) = Y(I)-YS(I).                 */
/*                                                                          */
/*      return 0 if OK, -1 if insufficient memory.                          */

int lowess(double *x,double *y,int n,double f,int nsteps,double delta,
            double *ys,double *rw,double *res) 
{
    register int i,j;
    int iter,ok,ns,nright,nleft,last,m1,m2,userw;
    double d1,d2,denom,alpha,cut,cmad,c1,c9,r;

    if (n < 2) {
        ys[1] = y[1];
        return(0);
    }

    /* at least two, at most n points */

    ns = imax(imin((int)(f * (double)n),n),2);
 
    for (iter = 1; iter <= nsteps + 1; iter++) {    /* robustness iterations */

        nleft = 1;
        nright = ns;
        last = 0;       /* index of prev estimated point */
        i = 1;          /* index of current point */
        while (1) {
            while(nright < n) { /* move nleft, nright to right if radius decreases */
                d1 = x[i] - x[nleft];
                d2 = x[nright + 1] - x[i];
                  /* if d1<=d2 with x[nright+1]==x[nright], lowest fixes */
                if (d1 <= d2)
                    break;
                                /* radius will not decrease by move right */
                nleft++;             
                nright++;               
            }
            if (iter > 1)
                userw = 1;
            else
                userw = 0;

            ok = lowest(x,y,n,x[i],&ys[i],nleft,nright,res,userw,rw);

            /* fitted value at x[i] */

            if (!ok)
                ys[i] = y[i];   /* all weights zero - copy over value (all rw==0) */

            if (last < i - 1) {     /* skipped points -- interpolate */

                denom = x[i] - x[last];         /* non-zero - proof? */

                for (j = last + 1; j < i; j++) {
                    alpha = (x[j] - x[last]) / denom;
                    ys[j] = alpha * ys[i] + (1.0 - alpha) * ys[last];
                }
            }
            last = i;                   /* last point actually estimated */
            cut = x[last] + delta;      /* x coord of close points */
            for(i = last + 1; i <= n; i++) {    /* find close points */
                if (x[i] > cut)                 /* i one beyond last pt within cut */
                    break;
                if(x[i] == x[last]) {   /* exact match in x */
                    ys[i] = ys[last];
                    last = i;
                }
            }
            i = imax(last + 1,i - 1);

            /* back 1 point so interpolation within delta, but always go forward */

            if (last >= n)
                break;
        }
        for (i = 1; i <= n; ++i)        /* residuals */
            res[i] = y[i] - ys[i];

        if (iter > nsteps)  /* compute robustness weights except last time */
            break;

        for (i = 1; i <= n; ++i)     
            rw[i] = fabs(res[i]);

        if (sortd(n,rw + 1,0))
            return(-1);

        m1 = 1 + n / 2;
        m2 = n - m1 + 1;

        cmad = 3.0 * (rw[m1] + rw[m2]);     /* 6 median abs resid */
        c9 = 0.999 * cmad;
        c1 = 0.001 * cmad;

        for (i = 1; i <= n; ++i) {
            r = fabs(res[i]);
            if (r <= c1)
                rw[i] = 1.0;    /* near 0, avoid underflow  */
            else if (r > c9)
                rw[i] = 0.0;    /* near 1, avoid underflow  */
            else
                rw[i] = pow(1.0 - pow(r / cmad,2.0),2.0);
        }
    }
    return(0);
}         
 
/* ------------------------------------------------------------------------ */
/*                                                                          */
/*      ok = lowest(x,y,n,xs,ys,nleft,nright,w,userw,rw)                    */
/*                                                                          */
/*      LOWEST is a support routine for LOWESS and  ordinarily  will        */
/*      not  be  called  by  the  user.   The  fitted  value, YS, is        */
/*      computed  at  the  value,  XS,  of  the   horizontal   axis.        */
/*      Robustness  weights,  RW,  can  be employed in computing the        */
/*      fit.                                                                */
/*                                                                          */
/*              X = Input; abscissas of the points on the                   */
/*                  scatterplot; the values in X must be ordered            */
/*                  from smallest to largest.                               */
/*              Y = Input; ordinates of the points on the                   */
/*                  scatterplot.                                            */
/*              N = Input; dimension of X,Y,W, and RW.                      */
/*             XS = Input; value of the horizontal axis at which the        */
/*                  smooth is computed.                                     */
/*             YS = Output; fitted value at XS.                             */
/*          NLEFT = Input; index of the first point which should be         */
/*                  considered in computing the fitted value.               */
/*         NRIGHT = Input; index of the last point which should be          */
/*                  considered in computing the fitted value.               */
/*              W = Output; W(I) is the weight for Y(I) used in the         */
/*                  expression for YS, which is the sum from                */
/*                  I = NLEFT to NRIGHT of W(I)*Y(I); W(I) is               */
/*                  defined only at locations NLEFT to NRIGHT.              */
/*          USERW = Input; logical variable; if USERW is .TRUE., a          */
/*                  robust fit is carried out using the weights in          */
/*                  RW; if USERW is .FALSE., the values in RW are           */
/*                  not used.                                               */
/*             RW = Input; robustness weights.                              */
/*             OK = Output; logical variable; if the weights for the        */
/*                  smooth are all 0.0, the fitted value, YS, is not        */
/*                  computed and OK is set equal to .FALSE.; if the         */
/*                  fitted value is computed OK is set equal to             */
/*                                                                          */
/*        Method                                                            */
/*                                                                          */
/*        The smooth at XS is computed using (robust) locally weighted      */
/*        regression of degree 1.  The tricube weight function is used      */
/*        with h equal to the maximum of XS-X(NLEFT) and X(NRIGHT)-XS.      */
/*        Two  cases  where  the  program  reverts to locally weighted      */
/*        regression of degree 0 are described  in  the  documentation      */
/*        for LOWESS.                                                       */
 
int lowest(double *x,double *y,int n,double xs,double *ys,int nleft,
           int nright,double *w,int userw,double *rw)
{
    register int j;
    int nrt,ok;
    double range,h,h1,h9,a,b,c,r;

    ok = 0;
    range = x[n] - x[1];
    h = dmax(xs - x[nleft],x[nright] - xs);
    h9 = 0.999 * h;
    h1 = 0.001 * h;
    a = 0.0;       /* sum of weights */

    for(j = nleft; j <= n; j++) {   /* compute weights (pick up all ties on right) */
        w[j] = 0.0;
        r = fabs(x[j] - xs);
        if (r <= h9) {              /* small enough for non-zero weight */
            if (r > h1)
                w[j] = pow(1.0 - pow(r / h,3.0),3.0);
            else    
                w[j] = 1.0;
            if (userw)
                w[j] = rw[j] * w[j];
            a = a + w[j];
        }
        else if (x[j] > xs)         /* get out at first zero wt on right */
            break;
    }
    nrt = j - 1;    /* rightmost pt (may be greater than nright because of ties) */
    if (a <= 0.0)
        ok = 0;   
    else {                          /* weighted least squares */
        ok = 1;    
        for (j = nleft; j <= nrt; ++j)  
            w[j] = w[j] / a;        /* make sum of w[j] == 1 */

        if (h > 0.0) {              /* use linear fit */
            a = 0.0;
            for (j = nleft; j <= nrt; ++j)
                a = a + w[j] * x[j];        /* weighted center of x values */
            b = xs - a;
            c = 0.0;

            for (j = nleft; j <= nrt; ++j)
                c = c + w[j] * pow(x[j] - a,2.0);
            if (sqrt(c) > 0.001 * range) {

                        /* points are spread out enough to compute slope */

                b = b / c;
                for (j = nleft; j <= nrt; ++j)
                    w[j] = w[j] * (1.0 + b * (x[j] - a));
            }
        }
        *ys = 0.0;
        for (j = nleft; j <= nrt; ++j)
            *ys = *ys + w[j] * y[j];
    }
    return(ok);
}         

/*--------------------------------------------------------------------------*/
/*  rmidmean(opt,n,x,y,r,xm,ym,ym1,ym2)                                     */
/*                                                                          */
/*  Return 0 if OK, -1 if error (insuff. memory)                            */

int rmidmean(int n,double *x,double *y,int r,double *xm,double *ym,
                                                    float *ym1,float *ym2)
{
    register int i,j;
    int k1,k2,err;
    double sm1,sm2;

    err = 0;
    if (alloc_actmp(n + 1)) {
        /* p_err(-2,1); */
        return(-1);  
    }
    for (i = 0; i < n; ++i) {
        nneighbor(n,x,r,i,&k1,&k2);
        xm[i] = midmean(r,x + k1,0,&err);
        if (err)
            break;

        for (j = 0; j < r; ++j)
            AcTmp[j] = y[j + k1];
        ym[i] = midmean(r,AcTmp,1,&err);
        if (err)
            break;

        for (j = 0; j < r; ++j)
            AcTmp[j] = y[j + k1] - ym[i];
        err = smidmean(r,AcTmp,&sm1,&sm2);
        if (err)
            break;
    
        ym1[i] = (float)(sm1 + ym[i]);
        ym2[i] = (float)(sm2 + ym[i]);
    }
    alloc_actmp(0);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  midmean(n,x,opt,err)                                                    */
/*                                                                          */
/*  return midmean of x[i], i=0,n-1. if opt == 0 the data are already       */
/*  sorted. err = 0 if OK, -1 if insuff. memory.                            */

double midmean(int n,double *x,int opt,int *err)
{
    register int i;
    int m;
    double qx1,qx2,tmp;

    *err = 0;
    if (opt) {      /* sort x[] */
        if (sortd(n,x,0)) { 
            /* p_err(-2,1); */
            *err = -1;
            return(0.0);
        }
    }
    qx1 = quantf(n,x,0.25);
    qx2 = quantf(n,x,0.75);

    tmp = 0.0;
    m = 0;
    for (i = 0; i < n; ++i) {
        if (x[i] >= qx1 && x[i] <= qx2) {
            tmp += x[i];
            m++;
        }
    }
    if (m > 0)
        tmp /= (double)m;

    return(tmp);
}

/*--------------------------------------------------------------------------*/
/*  smidmean(n,x,sm1,sm2)   x[i] i=0,n-1. return lower semi-midmean in sm1  */
/*                          and upper semi-midmean in sm2. the data are     */
/*                          not sorted.                                     */
/*                                                                          */
/*  Return 0 if OK, -1 if error (insuff. memory)                            */

int smidmean(int n,double *x,double *sm1,double *sm2)
{
    register int i;
    int m;
    double qx1,qx2,tmp;

    if (sortd(n,x,0)) { 
        /* p_err(-2,1); */
        return(-1);
    }
    qx1 = quantf(n,x,1.0 / 8.0);
    qx2 = quantf(n,x,3.0 / 8.0);

    tmp = 0.0;
    m = 0;
    for (i = 0; i < n; ++i) {
        if (x[i] >= qx1 && x[i] <= qx2) {
            tmp += x[i];
            m++;
        }
    }
    if (m > 0)
        tmp /= (double)m;

    *sm1 = tmp;

    qx1 = quantf(n,x,5.0 / 8.0);
    qx2 = quantf(n,x,7.0 / 8.0);

    tmp = 0.0;
    m = 0;
    for (i = 0; i < n; ++i) {
        if (x[i] >= qx1 && x[i] <= qx2) {
            tmp += x[i];
            m++;
        }
    }
    if (m > 0)
        tmp /= (double)m;

    *sm2 = tmp;
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  inpreg()        non-parametric regression with interval-valued          */
/*                  variables                                               */
/*                                                                          */  
/*                  npreg(                                                  */
/*                      k=...,      kernel, def. 1                          */
/*                                  1 : uniform                             */
/*                                  2 : triangle                            */
/*                                  3 : quartic                             */
/*                                  4 : Epanechnikov                        */
/*                      x=...,      list of points for X axis               */
/*                      d=...,      widths of interval, def. 1              */
/*                      df=...,     output file                             */
/*                      mxit=...,   max iterations, def. 10                 */
/*                      fmt=...,    print format, def. 10.4                 */
/*                  ) = YL,YH,XL,XH;                                        */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int inpreg(void)
{
    register int i,j;
    int err,nrec,ixl,ixh,iyl,iyh,m,iter,itmax;
    double xl,xh,yl,yh,xlm,xhm,ylm,yhm,w,x,a,b,ym;

    err = -1;         
    if (check_cmd(1))
        return(-1);

    printf1("Nonparametric regression with interval-valued variables.\n");
    printf1("Current memory: %d bytes.\n\n",MemReq);
         
    if (parm(CmdBuf + 6,4,1))     /* get parameters */
        goto INPFin;

    if (PMNV != 4) {            /* need four variables */
        p_err(-1,1);
        goto INPFin;
    }
                         /* PMNTP = number of arguments in PMTP[] (t=...) */
    if (PMNTP < 1) {
        printf1("Error: need points defined with x parameter.\n");
        goto INPFin;
    }
    if (PMD < EPSI1)
        PMD = 1.0;

    if (MxItFlg == 0)
        MxIter = 10;

    if (!PMF1Def) {             /* need an output file */
        p_err(-27,1);
        goto INPFin;
    }
    if (PMFmtF == 0)  
        pmfmt(10,4);

    if (PMOPT < 1 || PMOPT > 1)
        PMOPT = 1;

    iyl = PMVIdx[0];
    iyh = PMVIdx[1];
    ixl = PMVIdx[2];
    ixh = PMVIdx[3];

    if (alloc_acx(NOC + 1))     /* x lower bound */ 
        goto INPFin;
    if (alloc_acv(NOC + 1))     /* x upper bound */ 
        goto INPFin;
    if (alloc_acy(NOC + 1))     /* y lower bound */ 
        goto INPFin;
    if (alloc_acw(NOC + 1))     /* y upper bound */
        goto INPFin;

    if (alloc_actmp(NOC + 1))                       
        goto INPFin;
    if (alloc_actmp1(NOC + 1))                       
        goto INPFin;
    if (alloc_acz(NOC + 1))                       
        goto INPFin;
    if (alloc_act(NOC + 1))                       
        goto INPFin;
    if (alloc_acu(NOC + 1))                       
        goto INPFin;

    xlm = xhm = ylm = yhm = 0.0;
    for (i = 0; i < NOC; ++i) {
        xl = get_data(ixl,i);
        xh = get_data(ixh,i);
        yl = get_data(iyl,i);
        yh = get_data(iyh,i);
        if (yh < yl) {
            printf1("Error in case %d, found: %lg,%lg\n",i+1,yl,yh);
            goto INPFin;
        }
        AcX[i] = xl;
        AcV[i] = xh;
        AcY[i] = yl;
        AcW[i] = yh;
        xlm += xl;
        xhm += xh;
        ylm += yl;
        yhm += yh;
    }
    xlm /= (double)NOC;
    xhm /= (double)NOC;
    ylm /= (double)NOC;
    yhm /= (double)NOC;

    printf1("Smoothing means, kernel: ");
    switch (PMK) {
        case 2: printf1("triangle\n"); break;
        case 3: printf1("quartic\n"); break;
        case 4: printf1("Epanechnikov\n"); break;
        default: PMK = 1; printf1("uniform\n"); break;
    }
    printf1("Interval width: %lg\n",PMD);
    printf1("Max number of iterations: %d\n",MxIter);
    printf1("Number of cases: %d\n",NOC);
    printf1("Range of x means: ");
    printf1(PMFmtS,xlm);
    printf1(PMFmtS,xhm);
    printf1("\nRange of y means: ");
    printf1(PMFmtS,ylm);
    printf1(PMFmtS,yhm);
    newline();
/*
    for (i = 0; i < NOC; ++i)
        printf1("i=%d xl=%lf xh=%lf yl=%lf yh=%lf\n",i,AcX[i],AcV[i],AcY[i],AcW[i]);         
*/
    w = PMD / 2.0;

    itmax = nrec = 0;
    for (j = 0; j < PMNTP; ++j) {
  
        x = PMTP[j];
        a = x - w;
        b = x + w;

        m = 0;
        ylm = yhm = 0.0;
        for (i = 0; i < NOC; ++i) {
            if (AcX[i] <= b && AcV[i] >= a) {
                ylm += AcY[i];
                yhm += AcW[i];
                AcTmp[m] = AcY[i];
                AcTmp1[m] = AcW[i];
                AcT[m] = AcX[i];
                AcU[m] = AcV[i];
                m++;
            }
        }
/**
        printf("m=%d ",m);
        for (i = 0; i < m; ++i)
            printf(" (%lg,%lg) ",AcTmp[i],AcTmp1[i]);
        newline();

        printf("m=%d ",m);
        for (i = 0; i < m; ++i)
            printf(" (%lg,%lg) ",AcT[i],AcU[i]);
        newline();
**/
        if (m == 0)
            continue;

        ylm /= (double)m;
        yhm /= (double)m;


/*      printf("ylm=%lf yhm=%lf\n",ylm,yhm);*/
   
        /* first deal with minima */

        for (iter = 0; iter < MxIter; ++iter) {

            for (i = 0; i < m; ++i) {
                if (AcTmp[i] <= ylm) {
                    if (AcT[i] > x)
                        AcZ[i] = AcT[i];
                    else if (AcU[i] < x)
                        AcZ[i] = AcU[i];
                    else
                        AcZ[i] = x;
                }
                else {
                    if (fabs(x - AcT[i]) >= fabs(x - AcU[i]))
                        AcZ[i] = AcT[i];
                    else
                        AcZ[i] = AcU[i];
                }

            }
/**
            printf("Z::: ");
            for (i = 0; i < m; ++i)
                printf(" %lg ",AcZ[i]);
            newline();
**/
            ym = imkern(PMK,m,PMD,x,AcZ,AcTmp);

/*          printf("ym=%lf \n",ym); */

            if (ym >= ylm - EPSI1) {
                ylm = ym;
/*   printf("break ym=%lg ylm=%lg\n",ym,ylm);   */
                break;
            }
            ylm = ym;
        }
        itmax = imax(itmax,iter);

        /* now deal with maxima */

        for (iter = 0; iter < MxIter; ++iter) {

            for (i = 0; i < m; ++i) {
                if (AcTmp1[i] >= yhm) {
                    if (AcT[i] > x)
                        AcZ[i] = AcT[i];
                    else if (AcU[i] < x)
                        AcZ[i] = AcU[i];
                    else
                        AcZ[i] = x;
                }
                else {
                    if (fabs(x - AcT[i]) >= fabs(x - AcU[i]))
                        AcZ[i] = AcT[i];
                    else
                        AcZ[i] = AcU[i];
                }

            }
/*
            printf("Z::: ");
            for (i = 0; i < m; ++i)
                printf(" %lg ",AcZ[i]);
            newline();
*/
            ym = imkern(PMK,m,PMD,x,AcZ,AcTmp1);
  
/*          printf("ym=%lf \n",ym); */

            if (ym <= yhm + EPSI1) {
                yhm = ym;
                break;
            }
            yhm = ym;
        }
        itmax = imax(itmax,iter);

        fprintf(PMF1d,"%6d ",j + 1);
        fprintf(PMF1d,PMFmtS,x);
        fprintf(PMF1d,PMFmtS,ylm);                 
        fprintf(PMF1d,PMFmtS,yhm);                 
        fprintf(PMF1d,"\n");
        nrec++;
    }
    printf1("\n%d records written to: %s\n",nrec,PMF1dName);
    printf1("Max number of iterations used: %d\n",itmax);

    err = 0;
INPFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  imkern(opt,n,d,a,x,y)                                                   */
/*                                                                          */
/*  calculate mean of y at a, depending on type of kernel with bandwidth d. */
/*  opt = 1     rectangular                                                 */
/*        2     triangle                                                    */
/*        3     quartic                                                     */
/*        4     Epanechnikov                                                */
/*                                                                          */
/*  return mean of y.                                                       */

double imkern(int opt,int n,double d,double a,double *x,double *y)
{
    register int i;
    double ym,xm,dd,w;

    w = xm = ym = 0.0;
    dd = d / 2.0;

    for (i = 0; i < n; ++i) {
/*  printf("i=%d a=%lg xi=%lg yi=%lg\n",i,a,x[i],y[i]);     */
        if (fabs(x[i] - a) > dd)  
            continue;
  
        if (opt == 1)
            w = 1.0 / d;
        else if (opt == 2)
            w = (1.0 - fabs(a - x[i]) / dd) / dd;
        else if (opt == 3) {
            w = (a - x[i]) / dd;
            w = 1.0 - w * w;  
            w = 15.0 * w * w / (8.0 * d);
        }
        else if (opt == 4) {
            w = (a - x[i]) / dd;
            w = 1.0 - w * w;  
            w = 3.0 * w / (2.0 * d);
        }
        xm += w;
        ym += y[i] * w;
    }
    if (xm > 0.0)
        ym /= xm;
    return(ym);
}






