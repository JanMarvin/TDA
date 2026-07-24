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

/*  functions in t_smo.c */

int sma(void);
int s_sma(int n,double *x,int nw,double *w,int r,int opt);
int smd(void);
int s_smd(char *sm,int opt,int n,double *y);
int s_smdd(int n,double *x,double *y,int ns,int opt,int r,int dir);
double med3(double *x);
double med4(double *x);
double med5(double *x);
void nnsmooth(int n,double *x,double *y,int r);
void nneighbor(int n,double *x,int r,int k,int *k1,int *k2);
int spl(void);
int scplot(void);
void scplot1(int n,double *x,double *y,int ns);
int scplot2(int n);
void scplot2a(int n,double x,double y);

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

int sma(void)
{
    register int i;
    int err,n,iv,ivv,ix,nrec;
    double w;  

    err = -1;
    if (check_cmd(0))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    printf1("Moving averages. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 3,4,1))     /* get parameters */
        goto SMAFin;

    if (PMR < 1)
        PMR = 1;
    if (PMOPT != 2)
        PMOPT = 1;
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    if (PMNTP < 1) {
        printf1("Error: gss parameter must specify at least two weights.\n");
        goto SMAFin;
    }

    printf1("Weights:");
    for (i = PMNTP - 1; i >= 0; --i)
        printf1(" %g",PMTP[i]);
    for (i = 1; i < PMNTP; ++i)
        printf1(" %g",PMTP[i]);
    newline();

    n = PMNTP - 1;
    w = PMTP[0];
    for (i = 1; i <= n; ++i) {
        if (PMTP[i] < 0.0) {
            printf1("Error: negative weights are not allowed.\n");
            goto SMAFin;
        }
        w += 2.0 * PMTP[i];
    }
    printf1("Sum of weights: %g\n",w);

    printf1("Using: ");
    if (PMOPT == 1)
        printf1("copy-on end-value rule.\n");
    else            
        printf1("replicate end-value rule.\n");
    printf1("Repeat factor: %d\n",PMR);

    if (alloc_acx(NOC))         /* allocate AcX */
        goto SMAFin;

    if (alloc_acy(NOC))         /* allocate AcY */
        goto SMAFin;
     
    nrec = ix = 0;
    for (iv = 0; iv < PMNV; ++iv) {     /* for all variables */

        ivv = PMVIdx[iv];

        if (PMF1Def)  
            fprintf(PMF1d,"# Var: %s\n",VName[ivv]);
        else
            printf1("\nVariable: %s\n",VName[ivv]);

        n = getd1(AcX,ivv,-1,1);
        if (n == 0)  
            continue;    
             
        for (i = 0; i < n; ++i)
            AcY[i] = AcX[i];

        s_sma(n,AcY,PMNTP,PMTP,PMR,PMOPT);

        for (i = 0; i < n; ++i) {
            if (PMF1Def) {
                fprintf(PMF1d,"%3d %6d ",ix,i + 1);
                fprintf(PMF1d,PMFmtS,AcX[i]);
                fprintf(PMF1d,PMFmtS,AcY[i]);
                fprintf(PMF1d,"\n");
                nrec++;
            }
            else {
                printf1("%3d %6d ",ix,i + 1);
                printf1(PMFmtS,AcX[i]);
                printf1(PMFmtS,AcY[i]);
                printf1("\n");
            }
        }
        ix++;
    }
    if (PMF1Def)  
        printf1("%d records written to: %s\n",nrec,PMF1dName);

    err = 0;

SMAFin:
    p_clean();
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

int s_sma(int n,double *x,int nw,double *w,int r,int opt)
{
    register int i,j,k;
    int ii,nn,na;
    double tmp;

    if (n < 1 || nw < 1)
        return(0);
        
    if (alloc_actmp(n))         /* allocate AcTmp */
        return(-1); 

    for (i = 0; i < n; ++i)
        AcTmp[i] = x[i];

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
            AcTmp[i] = tmp;           
        }
        for (i = 0; i < n; ++i)
            x[i] = AcTmp[i];

        if (++ii >= r)
            break;
    }
    alloc_actmp(0);     /* free AcTmp */
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

int smd(void)
{
    register int i;
    int err,n,iv,ivv,ix,nrec;

    err = -1;
    if (check_cmd(0))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    printf1("Running Median Smoother. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 3,4,1))     /* get parameters */
        goto SMDFin;

    if (PMRHSTRA < 1 || PMOPT < 1 || PMOPT > 3) {
        p_err(-1,1);
        goto SMDFin;
    }
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    printf1("Using: sm=[%s]\nUsing: ",PMRHSTR);
    if (PMOPT == 1)
        printf1("copy-on end-value rule.\n");
    else if (PMOPT == 2)
        printf1("replicate end-value rule.\n");
    else
        printf1("end-value smoothing.\n");

    if (alloc_acx(NOC))         /* allocate AcX */
        goto SMDFin;

    if (alloc_acy(NOC))         /* allocate AcY */
        goto SMDFin;

    if (PMF1Def)  
        fprintf(PMF1d,"# Running medians [%s]\n",PMRHSTR);

    nrec = ix = 0;
    for (iv = 0; iv < PMNV; ++iv) {     /* for all variables */

        ivv = PMVIdx[iv];

        if (PMF1Def)  
            fprintf(PMF1d,"# Var: %s\n",VName[ivv]);
        else
            printf1("\nVariable: %s\n",VName[ivv]);

        n = getd1(AcX,ivv,-1,1);
        if (n == 0)  
            continue;    

        for (i = 0; i < n; ++i)
            AcY[i] = AcX[i];
           
        if (s_smd(PMRHSTR,PMOPT,n,AcY)) {
            printf1("Syntax error: [%s]\n",PMRHSTR);
            goto SMDFin;
        }
        for (i = 0; i < n; ++i) {
            if (PMF1Def) {
                fprintf(PMF1d,"%3d %6d ",ix,i + 1);
                fprintf(PMF1d,PMFmtS,AcX[i]);
                fprintf(PMF1d,PMFmtS,AcY[i]);
                fprintf(PMF1d,"\n");
                nrec++;
            }
            else {
                printf1("%3d %6d ",ix,i + 1);
                printf1(PMFmtS,AcX[i]);
                printf1(PMFmtS,AcY[i]);
                printf1("\n");
            }
        }
        ix++;
    }
    if (PMF1Def)  
        printf1("%d records written to: %s\n",nrec,PMF1dName);
    err = 0;

SMDFin:
    p_clean();
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

int s_smd(char *sm,int opt,int n,double *y)
{
    register int i;
    int r,dir,tflag,err;
    register char *p;   
    double tmp;

    if (n <= 1)
        return(0);

    err = -2;
    if (alloc_acu(n))         /* allocate AcU */
        goto SSMDFin;

    if (alloc_acv(n))         /* allocate AcV */
        goto SSMDFin;

    err = -1;
    for (i = 0; i < n; ++i)  
        AcV[i] = y[i];

    tflag = 0;
SMDRep:      
    dir = 0;
    p = sm;         
    while (*p) {
        r = 0;
        if (*p == '2') {
            p++;
            s_smdd(n,y,AcU,2,opt,r,dir);
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
            s_smdd(n,y,AcU,3,opt,r,0);
        }
        else if (*p == '4') {
            p++;
            s_smdd(n,y,AcU,4,opt,r,dir);
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
            s_smdd(n,y,AcU,5,opt,r,0);
        }
        else if (*p == 'H') {
            p++;
            s_smdd(n,y,AcU,0,opt,0,0);
        }
        else if (*p == 'S') {
            p++;
            s_smdd(n,y,AcU,1,opt,0,0);
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
            tmp = AcV[i];
            AcV[i] = y[i];
            y[i] = tmp - y[i];
        }
        goto SMDRep;
    }
    else if (tflag == 2) {
        for (i = 0; i < n; ++i)  
            y[i] += AcV[i];
    }
    err = 0;

SSMDFin:
    alloc_acu(0);       /* free AcU */
    alloc_acv(0);       /* free AcV */
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

int s_smdd(int n,double *x,double *y,int ns,int opt,int r,int dir)
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
                    y[i] = med3(s);

                    s[0] = x[i + 1];
                    s[1] = x[i + 2];
                    s[2] = 3.0 * x[i + 2] - 2.0 * x[i + 3];
                    y[i + 1] = med3(s);
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
                    y[i] = med3(x + i - 1);
            }
            else if (ns == 4) {
                if (i > 0 && i < nn) {
                    if (dir == 0) {
                        if (i == 1)
                            y[i] = med3(x);
                        else             
                            y[i] = med4(x + i - 2);
                    }
                    else {
                        if (i == nn - 1)
                            y[i] = med3(x + nn - 2);
                        else
                            y[i] = med4(x + i - 1);
                    }
                }
            }
            else if (ns == 5) {
                if (i > 0 && i < nn) {
                    if (i == 1)
                        y[i] = med3(x);
                    else if (i == nn - 1)
                        y[i] = med3(x + nn - 2);
                    else   
                        y[i] = med5(x + i - 2);
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
                y[0] = med3(s);
                s[0] = 3.0 * y[nn - 1] - 2.0 * y[nn - 2];
                s[1] = x[nn];
                s[2] = y[nn - 1];
                y[0] = med3(s);
            }
        }
        j = 0;
        for (i = 0; i < n; ++i) {
            if (fabs(x[i] - y[i]) > EPSI1)
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

double med3(double *x)
{
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

double med5(double *x)
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
    return(med3(s));
}

/* ------------------------------------------------------------------------ */
/*  med4(x)     return median of x[0],x[1],x[2],x[3]                        */

double med4(double *x)
{
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

void nnsmooth(int n,double *x,double *y,int r)
{
    register int i,j;
    int k1,k2;
    double tmp;

    for (i = 0; i < n; ++i) {
        y[i] = x[i];
        if (r < n) {
            nneighbor(n,x,r,i,&k1,&k2);
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

void nneighbor(int n,double *x,int r,int k,int *k1,int *k2)
{
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

int spl(void)
{
    register int i,j,l,k;
    int err,r,ix,iy,iw,nrec,nk;
    double xa,xb,x,y,eps;

    err = -1;
    if (check_cmd(0))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    printf1("Smoothing spline function. Current memory: %d bytes.\n",MemReq);
    MxIter = 0;

    if (parm(CmdBuf + 3,4,1))     /* get parameters */
        goto SPLFin;

    if (PMNV != 2 && PMNV != 3) {
        printf1("Error: need two or three variables on right-hand side.\n");
        goto SPLFin;
    }
    if (PMSIG < 0.0)
        PMSIG = 0.0;
    if (PMDEG < 2)
        PMDEG = 2;
    if (PMMax < 3 * PMDEG + 1)
        PMMax = imax(NOC / 2 + 1,3 * PMDEG + 1);
    if (MxIter < 1)
        MxIter = 10;
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    printf1("Degree: %d. Smoothing factor: %g\n",PMDEG,PMSIG);
    printf1("Max number of iterations: %d\n",MxIter);
    printf1("Max number of knots: %d\n",PMMax);

    if (alloc_acx(NOC + 1))         /* allocate AcX */
        goto SPLFin;

    if (alloc_acy(NOC + 1))         /* allocate AcY */
        goto SPLFin;
     
    if (alloc_acw(NOC + 1))         /* allocate AcW */
        goto SPLFin;
     
    if (alloc_acu(PMMax + 1))       /* allocate AcU */
        goto SPLFin;
     
    if (alloc_acv(PMMax + 1))       /* allocate AcV */
        goto SPLFin;
     
    if (alloc_actmp(PMDEG + 2))     /* allocate AcTmp */
        goto SPLFin;
     
    ix = PMVIdx[0];
    iy = PMVIdx[1];
    if (PMNV > 2)
        iw = PMVIdx[2];
    else
        iw = -1;

    for (i = 0; i < NOC; ++i) {
        AcX[i] = get_data(ix,i);
        AcY[i] = get_data(iy,i);
        if (iw >= 0) {
            AcW[i] = get_data(iw,i);
            if (AcW[i] <= 0.0) {
                printf1("Error: found zero or negative weight in case %d.\n",i + 1);
                goto SPLFin;
            }
        }
        else
            AcW[i] = 1.0;
    }
    if (sortd3(NOC,AcX,AcY,AcW))        /* sort in ascending order */
        goto SPLFin;

    for (i = 1; i < NOC; ++i) {
        if (AcX[i] <= AcX[i - 1]) {
            printf1("Error: cannot sort into strictly ascending order.\n");
            goto SPLFin;
        }   
    }
    nk = PMMax;
    xa = AcX[0];
    xb = AcX[NOC - 1];

    r = splf (NOC,AcX - 1,AcY - 1,AcW - 1,xa,xb,PMDEG,PMSIG,&nk,AcU,AcV,MxIter);

    if (r < 0) {
        printf1("Error (%d)",r);
        if (r == -2) 
            printf1(": exceeded max number of knots");
        else if (r == -3)
            printf1(": exceeded max number of iterations");
        printf1(".\n");
        goto SPLFin;
    }
    printf1("Successful calculation of ");
    if (r == 0)
        printf1("smoothing spline.\n");
    else if (r == 1)
        printf1("interpolating spline.\n");
    else
        printf1("OLS approximation.\n");
    printf1("Actual number of knots: %d\n",nk);

    xa = AcU[PMDEG + 1];
    xb = AcU[nk - PMDEG];

    printf1("Interval for interpolation: [%g,%g]\n",xa,xb);

    if (PMF1Def) {
        eps = 1.e-6;
        nrec = 0;
        if (PMRXFlg) {
            printf1("Request is interpolation in [%g,%g], increment: %g\n",
                                            PMRXA,PMRXB,PMRXD);

            r = (int)((PMRXB - PMRXA) / PMRXD) + 1; 
            x = PMRXA;
            interv(nk,AcU,x,1);           /* initialization */
            k = 1;
            for (i = 0; i <= r; ++i) {
                x = PMRXA + (double)i * PMRXD;
                if (x < xa)
                    continue;
                if (x > xb + eps)
                    break;

                fprintf(PMF1d,"%6d ",k++);
                fprintf(PMF1d,PMFmtS,x);

                if (fabs(x - xb) <= eps)
                    l = nk - PMDEG - 1;
                else
                    l = interv(nk,AcU,x,0);

                for (j = 0; j <= PMDEG; ++j) {
                    y = spld(nk,AcU,PMDEG,AcV,j,x,l,AcTmp);
                    fprintf(PMF1d,PMFmtS,y);
                }
                if (iw >= 0)
                    fprintf(PMF1d,PMFmtS,AcW[i]);

                fprintf(PMF1d,"\n");
                nrec++;
            }
        }
        else {
            interv(nk,AcU,x,1);           /* initialization */
            for (i = 0; i < NOC; ++i) {
                x = AcX[i];
                if (x < xa)
                    continue;
                if (x > xb)
                    break;

                fprintf(PMF1d,"%6d ",i + 1);
                fprintf(PMF1d,PMFmtS,x);
                fprintf(PMF1d,PMFmtS,AcY[i]);

                if (fabs(x - xb) <= eps)
                    l = nk - PMDEG - 1;
                else
                    l = interv(nk,AcU,x,0);

                for (j = 0; j <= PMDEG; ++j) {
                    y = spld(nk,AcU,PMDEG,AcV,j,x,l,AcTmp);
                    fprintf(PMF1d,PMFmtS,y);
                }
                if (iw >= 0)
                    fprintf(PMF1d,PMFmtS,AcW[i]);

                fprintf(PMF1d,"\n");
                nrec++;
            }

        }
        printf1("%d records written to: %s\n",nrec,PMF1dName);
    }

    err = 0;

SPLFin:
    p_clean();
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

int scplot(void)
{
    register int i;
    int nn,err,ix,iy,first,nrec;
    double x,y;
           
    err = - 1;
    nrec = 0;
    if (check_cmd(0))
        return(-1);

    printf1("Scatterplot. Current memory: %d bytes.\n",MemReq);

    if (check_ps(2))
        return(-1);

    if (parm(CmdBuf + 6,4,1))       /* get parameters */
        goto PLSCFin;

    if (PMOPT < 1 || PMOPT > 4)
        PMOPT = 1;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    if (NOC < 2)
        goto PLSCFin;

    if (PMNV != 2) {
        p_err(-1,1);
        goto PLSCFin;
    }
    ix = PMVIdx[0];
    iy = PMVIdx[1];

    if (alloc_acx(NOC + 1))
        goto PLSCFin;
    if (alloc_acy(NOC + 1))
        goto PLSCFin;

    nn = getd2(AcX,AcY,ix,iy,1);
    if (nn < 2) {
        p_err(-26,1);
        goto PLSCFin;
    }
    if (PMOPT > 1) {
        if (sortd2(nn,AcX,AcY))     /* sort simultaneously according to AcX */
            goto PLSCFin;
    }
    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    fprintf(PSFd,"gsave\n");
    if (PMNC != 1)
        set_clip();
    ps_lwidth(PMLW);
    ps_ltyp(PMLT);

    if (PMOPT != 2) {   /* standard symbols */

        if (PMS >= 1 && PMS <= 17 && PMFS > 0.0) {  /* plot symbols */

            fprintf(PSFd,"gsave\n");
            fprintf(PSFd,"/ssiz %4.2f def\n",PtMM * PMFS / 2.0);

            for (i = 0; i < nn; ++i) {
                x = ps_2dx(AcX[i]);       
                y = ps_2dy(AcY[i]);       
                ps_sym(PMS,x,y,PtMM * PMFS / 2.0);
            }
            fprintf(PSFd,"grestore\n");
        }
        /***
        ps_lwidth(PMLW);
        ps_ltyp(PMLT);
        ***/
    }

    if (PMOPT == 2) {               /* sunflower plot */
        err = scplot2(nn);
    }
    else if (PMOPT == 3) {          /* lowess */

        if (alloc_acu(nn + 1))
            goto PLSCFin;
        if (alloc_acv(nn + 1))
            goto PLSCFin;
        if (alloc_actmp(nn + 1))
            goto PLSCFin;

        if (PMNS < 0)               /* number of steps */
            PMNS = 2;
        if (PMSIG <= 0.0 || PMSIG >= 1.0)
            PMSIG = 0.5;
        if (PMD < 0.0)
            PMD = 0.0;
        printf1("Lowess: sig=%g ns=%d d=%g\n",PMSIG,PMNS,PMD);

        if (lowess(AcX - 1,AcY - 1,nn,PMSIG,PMNS,PMD,AcU,AcV,AcTmp))
            goto PLSCFin;

        first = 0;
        for (i = 0; i < nn; ++i) {
            ps_2dplot(AcX[i],AcU[i + 1],first);    
            first = 1;
        }
        fprintf(PSFd,"stroke\n");

        if (PMF1Def) {
            for (i = 0; i < nn; ++i) {
                fprintf(PMF1d,"%6d ",i + 1);
                fprintf(PMF1d,PMFmtS,AcX[i]);
                fprintf(PMF1d,PMFmtS,AcY[i]);
                fprintf(PMF1d,PMFmtS,AcU[i + 1]);
                fprintf(PMF1d,"\n");
                nrec++;
            }
        }
    }
    else if (PMOPT == 4) {                  /* midmeans */

        if (alloc_acu(nn + 1))
            goto PLSCFin;
        if (alloc_acv(nn + 1))
            goto PLSCFin;
        if (alloc_acxf(nn + 1))
            goto PLSCFin;
        if (alloc_acyf(nn + 1))
            goto PLSCFin;

        if (PMR < 2)
            PMR = 2;
        if (PMR >= nn) {
            printf1("Error: r must be less than number of cases.\n");
            goto PLSCFin;
        }
        printf1("Midmeans: r = %d\n",PMR);

        if (rmidmean(nn,AcX,AcY,PMR,AcU,AcV,AcXF,AcYF)) {
            p_err(-2,1);
            goto PLSCFin;
        }
        if (PMF1Def) {
            for (i = 0; i < nn; ++i) {
                fprintf(PMF1d,"%6d ",i + 1);
                fprintf(PMF1d,PMFmtS,AcX[i]);
                fprintf(PMF1d,PMFmtS,AcY[i]);
                fprintf(PMF1d,PMFmtS,AcU[i]);
                fprintf(PMF1d,PMFmtS,AcV[i]);
                fprintf(PMF1d,PMFmtS,(double)AcXF[i]);
                fprintf(PMF1d,PMFmtS,(double)AcYF[i]);
                fprintf(PMF1d,"\n");
                nrec++;
            }
        }

        if (PMNS <= 0 || PMNS >= nn)
            PMNS = 0;
        else {
            if (alloc_actmp(nn + 1))
                goto PLSCFin;
        }
        scplot1(nn,AcU,AcV,PMNS);

        for (i = 0; i < nn; ++i)
            AcV[i] = (double)AcXF[i];
        scplot1(nn,AcU,AcV,PMNS);

        for (i = 0; i < nn; ++i)
            AcV[i] = (double)AcYF[i];
        scplot1(nn,AcU,AcV,PMNS);

    }
    fprintf(PSFd,"grestore\n");

    if (PMF1Def && nrec > 0)
        printf1("%d records written to: %s\n",nrec,PMF1dName);

    err = 0;

PLSCFin:
    p_clean();       
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  scplot1(n,x,y,ns)   plot x,y. if ns > 0 then smooth.                    */

void scplot1(int n,double *x,double *y,int ns)
{
    register int i,first;

    if (ns)  
        nnsmooth(n,y,AcTmp,PMNS);
       
    first = 0;
    for (i = 0; i < n; ++i) {
        if (ns)
            ps_2dplot(x[i],AcTmp[i],first);    
        else
            ps_2dplot(x[i],y[i],first);    
        first = 1;
    }
    fprintf(PSFd,"stroke\n");
}

/*--------------------------------------------------------------------------*/
/*  scplot2     Sunflower plot                                              */
/*              n = number of points in AcX,AcY                             */
/*              return 0 if OK, -1 if insufficient memory.                  */

int scplot2(int n)
{
    register int i,j,ix,iy;
    int err,nx,ny,px,py;
    double ax,ay,x,y;

    err = -1;

    nx = PMNN1;
    ny = PMNN2;
    if (nx < 1)
        nx = 1;             
    if (ny < 1)
        ny = 1;           
    ax = UXLen / (double)nx;
    ay = UYLen / (double)ny;

    printf1("Sunflower plot: nx=%d [%g] ny=%d [%g]\n",nx,ax,ny,ay);

    if (alloc_acn(nx * ny))
        goto SCP2Fin;

    for (i = 0; i < n; ++i) {
        px = (int)((AcX[i] - PA1[0]) / ax);
        py = (int)((AcY[i] - PA1[1]) / ay);
        if (px < 0)
            px = 0;
        if (px >= nx)
            px = nx - 1;
        if (py < 0)
            py = 0;
        if (py >= ny)
            py = ny - 1;
        AcN[py * nx + px] += 1;
    }
    fprintf(PSFd,"/ssiz %4.2f def\n",PtMM * PMFS / 5.0);

    /*****
    xs = 0.35 * ax * PSXLen / UXLen;    
    ys = 0.35 * ay * PSYLen / UYLen;
    ******/

    y = PA1[1] + ay / 2.0;
    for (iy = 0; iy < ny; ++iy) {
        x = PA1[0] + ax / 2.0;
        for (ix = 0; ix < nx; ++ix) {
            j = AcN[iy * nx + ix];
            if (j > 0)                  /* plot */
                scplot2a(j,x,y);
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

void scplot2a(int n,double x,double y)
{
    int i;
    double a,b,siz;

    x = ps_2dx(x);
    y = ps_2dy(y);
    fprintf(PSFd,"gsave\n%5.2f %5.2f translate\n",x,y);

    if (n == 1)
        ps_sym(5,0.0,0.0,PtMM * PMFS / 5.0);      /* always symbol 5 */

    else if (n >= 2) {
        siz = PMFSS * PtMM / 2.0;    /* size in points */
        a = 2.0 * Pi / (double)n;
        b = Pi / 2.0;
        for (i = 0; i < n; ++i) {
            x = siz * cos(a * (double)i + b);
            y = siz * sin(a * (double)i + b);
            fprintf(PSFd,"%6.2f %6.2f m\n0 0 l\nstroke\n",x,y);             
        }
        /** fprintf(PSFd,"stroke\n"); **/         
    }
    fprintf(PSFd,"grestore\n");
}



