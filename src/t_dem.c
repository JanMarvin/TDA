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

/*  functions in t_dem.c */

int spmod(void);
double enorm(int n,double *x);
double r_alph(int n,double *x,double *y);
int etest(void);
int rap(void);
double rap_rate(int i,int j,int n,int na,int ny);

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

int spmod(void)
{
    register int i,j;
    int err,iter,i0,i1,i2,conv;
    double alpha,a;

    err = -1;
    if (check_cmd(1))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    printf1("Stationary solution of Leslie-matrix. Current memory: %d bytes.\n",MemReq);

    TOLF = 1.0e-6;

    if (parm(CmdBuf + 5,4,1))     /* get parameters */
        goto DEMFin;
   
    if (PMNV != 3) {
        printf1("Error: need three variables.\n");
        goto DEMFin;
    }
    if (PMFmtF == 0)
        pmfmt(8,6);

    i0 = PMVIdx[0];
    i1 = PMVIdx[1];
    i2 = PMVIdx[2];

    /* allocate memory */

    if (alloc_acx(NOC + 1))
        goto DEMFin;
    if (alloc_acy(NOC + 1))
        goto DEMFin;

    if (MxItFlg == 0)
        MxIter = 100;
    printf1("Max number of iterations: %d\n",MxIter);
    printf1("Tolerance for convergence: %lg\n\n",TOLF);

    for (i = 0; i < NOC; ++i)
        AcX[i] = get_data(i2,i);
    enorm(NOC,AcX);

    if (PMProtFDef) {
        fprintf(PMProtFd,"%6d ",iter);
        for (i = 0; i < NOC; ++i)
            fprintf(PMProtFd,PMFmtS,AcX[i]);
        fprintf(PMProtFd,"\n");
    }
    conv = 0;
    for (iter = 1; iter <= MxIter; ++iter) {

        AcY[0] = 0.0;
        for (i = 0; i < NOC; ++i)
            AcY[0] += AcX[i] * get_data(i0,i);

        for (i = 1; i < NOC; ++i)
            AcY[i] = AcX[i - 1] * get_data(i1,i - 1);

        alpha = r_alph(NOC,AcX,AcY);
        enorm(NOC,AcY);

        if (PMProtFDef) {
            fprintf(PMProtFd,"%6d ",iter);
            for (i = 0; i < NOC; ++i)
                fprintf(PMProtFd,PMFmtS,AcY[i]);
            fprintf(PMProtFd,PMFmtS,alpha);
            fprintf(PMProtFd,"\n");
        }

        if (iter > 0) {
            if (fabs(alpha - a) <= TOLF) {
                conv = 1;
                break;
            }
        }
        for (i = 0; i < NOC; ++i)  
            AcX[i] = AcY[i];
        a = alpha;
    }
    if (conv == 0) {        /* if no convergence */
        printf1("No convergence.\n");
        goto DEMFin;
    }
    a = 0.0;
    for (i = 0; i < NOC; ++i)    
        a = dmax(a,fabs(AcX[i] - AcY[i]));

    printf1("Convergence after %d iterations.\n",iter);
    printf1("Final growth factor: ");
    printf1(PMFmtS,alpha);
    printf1("\nMax difference in stationary vector: %lg\n",a);

    if (PMF1Def) {
        for (i = 0; i < NOC; ++i)
            AcX[i] = get_data(i2,i);
        enorm(NOC,AcX);

        for (i = 0; i < NOC; ++i) {
            fprintf(PMF1d,"%4d ",i + 1);
            fprintf(PMF1d,PMFmtS,AcX[i]);
            fprintf(PMF1d,PMFmtS,AcY[i]);
            fprintf(PMF1d,"\n");
        }
        printf1("Results written to: %s\n",PMF1dName);
    }   
    if (PMProtFd)
        printf1("Protocol written to: %s\n",PMProtFName);

    a = 0.0;
    for (i = 0; i < NOC; ++i) {
        a += AcY[i];
    }
    printf("Sum of Y = %g\n",a);


    for (i = 0; i < NOC; ++i) {
        a = 0.0;
        if (i == 0) {
            for (j = 0; j < NOC; ++j) 
                a += AcY[j] * get_data(i0,j);
        }
        else {
            a = get_data(i1,i - 1) * AcY[i - 1];
        }
        AcX[i] = a;
    }
    a = 0.0;
    for (i = 0; i < NOC; ++i) {
        a += AcX[i];
    }
    printf("Sum of next Y = %g\n",a);

    for (i = 0; i < NOC; ++i) {
        printf("%4d %16.8lf %16.8lf %16.8lf\n",i+1,AcY[i],AcX[i],AcY[i] * alpha);
    }
                




    err = 0;

DEMFin:
    p_clean();  
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  enorm(n,x)  return sum of components of x, scale x to                   */
/*              sum of components = 1.                                      */

double enorm(int n,double *x)
{
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

double r_alph(int n,double *x,double *y)
{
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

int etest(void)
{
    register int i,j,k;
    int err,r,row,col,ivflg;
    double dr,di,tmp,tmp1;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Test of eigenvalue/vector calculation. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,8,1))     /* get parameters */
        goto ETFin;
   
    if (PMFmtF == 0)
        pmfmt(-19,11);

    if (n_mexpr(PMRHSTR,1,0,&row,&col,&ivflg,NULL) == NULL)
        goto ETFin;
   
    if (row != col) {
        mat_err(1);
        goto ETFin;
    }
    if (PMALG < 1 || PMALG > 3)
        PMALG = 1;

    printf1("Algorithm: %d",PMALG);
    if (MxItFlg == 0)
        MxIter = 100;
    if (PMALG == 3)
        printf1(" [max iterations: %d]",MxIter);
    newline();

    if (alloc_acw(row * col + 1))
        goto ETFin;

    for (i = 1; i <= row * col; ++i)
        AcW[i] = MX[0][i];

    printf1("\n%s\n",PMRHSTR);
    for (i = 1; i <= row; ++i) {
        for (j = 1; j <= col; ++j)  
            printf1(PMFmtS,AcW[(i - 1) * col + j]);
        newline();
    }

    if (alloc_acn(col + 1))
        goto ETFin;
    if (alloc_acx(col + 1))
        goto ETFin;
    if (alloc_acy(col + 1))
        goto ETFin;
    if (alloc_acu(row * col + 1))
        goto ETFin;
    if (alloc_acv(row * col + 1))
        goto ETFin;

    if (PMALG == 1)  
        r = eigen(col,MX[0],AcX,AcY,AcU,1);
    else if (PMALG == 2)  
        r = eigen1(col,MX[0],AcX,AcY,AcU,AcV,AcN);
    else {
        r = eigen2(col,MX[0],AcX,AcY,AcU,1,MxIter);
        if (r > 0 && r < MxIter) {
            printf1("Required %d iterations.\n",r);
            r = 0;
        }
    }       
     
    if (r) {
        printf1("No success in calculations (err=%d).\n",r);
        goto ETFin;
    }
     
    if (PMALG != 2) {
        for (k = 1; k <= col; ++k) {
            if (fabs(AcY[k]) != 0.0) {
                for (i = 1; i <= row; ++i) {
                    AcV[(i - 1) * col + k] = AcU[(i - 1) * col + k + 1];
                    AcV[(i - 1) * col + k + 1] = -AcU[(i - 1) * col + k + 1];
                    AcU[(i - 1) * col + k + 1] = AcU[(i - 1) * col + k];
                }
                k++;
            }
        }
    }
    for (k = 1; k <= col; ++k) {
   
        printf1("\n(%3d) ",k);
        printf1(PMFmtS,AcX[k]);
        printf1(PMFmtS,AcY[k]);
        printf1("  [eigenvalue]\n\n");

        for (i = 1; i <= row; ++i) {
            printf1("      ");
            printf1(PMFmtS,AcU[(i - 1) * col + k]);
            printf1(PMFmtS,AcV[(i - 1) * col + k]);
            if (i == 1)
                printf1("  [eigenvector]");
            newline();
        }

        /* check */

        dr = di = 0.0;
        for (i = 1; i <= row; ++i) {
            tmp = 0.0;
            for (j = 1; j <= col; ++j)  
                tmp += AcW[(i - 1) * col + j] * AcU[(j - 1) * col + k];
            tmp1 = AcX[k] * AcU[(i - 1) * col + k] -
                               AcY[k] * AcV[(i - 1) * col + k];
            dr = dmax(dr,fabs(tmp - tmp1));
        }
        for (i = 1; i <= row; ++i) {
            tmp = 0.0;
            for (j = 1; j <= col; ++j)  
                tmp += AcW[(i - 1) * col + j] * AcV[(j - 1) * col + k];
            tmp1 = AcX[k] * AcV[(i - 1) * col + k] +
                               AcY[k] * AcU[(i - 1) * col + k];
            di = dmax(di,fabs(tmp - tmp1));
        }
        printf1("\n      ");
        printf1(PMFmtS,dr);
        printf1(PMFmtS,di);
        printf1("  [check]\n");
    }
    err = 0;

ETFin:
    mx_free();
    p_clean();  
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

int rap(void)
{
    register int i,j;
    int err,its,itc,itf,id,mints,maxts,d,t,ts,tc,tf,nn,ny,na,a,nr,ne;
    double r;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Rates in age-period form. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 3,4,1))     /* get parameters */
        goto RAPFin;
   
    if (PMNV != 4) {
        printf1("Error: need four variables on right-hand side.\n");
        goto RAPFin;
    }
    if (PMNC < 0)
        PMNC = 0;   

    its = PMVIdx[0];
    itc = PMVIdx[1];
    itf = PMVIdx[2];
    id  = PMVIdx[3];

    for (i = 0; i < NOC; ++i) {
        ts = get_data(its,i);
        tc = get_data(itc,i);
        tf = get_data(itf,i);
        if (ts > tc || tc > tf) {
            printf1("Inconsistent data in case %d\n",i + 1);
            goto RAPFin;
        }
    }
    mints = maxts = (int)get_data(its,0);

    for (i = 1; i < NOC; ++i) {
        t = get_data(its,i);
        mints = imin(mints,t);
        maxts = imax(maxts,t);
    }

    if (mints < 0) {
        printf1("Error: minimal starting time is negative: %d\n",mints);
        goto RAPFin;
    }
    nn = maxts - mints + 1;

    if (alloc_acn(nn + 1))
        goto RAPFin;
    if (alloc_aci(nn + 1))
        goto RAPFin;
    if (alloc_acr(nn + 1))
        goto RAPFin;
    if (alloc_acs(nn + 1))
        goto RAPFin;

    for (i = 0; i < nn; ++i)
        AcR[i] = AcS[i] = -1;

    for (i = 0; i < NOC; ++i) {
        ts = get_data(its,i);
        tc = get_data(itc,i);
        tf = get_data(itf,i);
        d  = get_data(id ,i);
        j = ts - mints;
        AcN[j] += 1;
        if (d != 0)
            AcI[j] += 1;

        if (AcR[j] < 0)
            AcR[j] = tc;
        else
            AcR[j] = imin(AcR[j],tc);

        if (AcS[j] < 0)
            AcS[j] = tf;
        else
            AcS[j] = imax(AcS[j],tf);
    }

    printf1("\n  Idx       TS        N      D=1      D=0   Min TC   Max TF\n");
    prnchar('-',59,1);

    j = 1;
    for (i = 0; i < nn; ++i) {
        if (AcN[i] < 1)
            continue;
        printf1("%5d %8d %8d %8d %8d %8d %8d\n",j++,mints + i,AcN[i],
            AcI[i],AcN[i] - AcI[i],AcR[i],AcS[i]);
    }
    newline();

    if (PMYEAR1 < 0 || PMAGE1 < 0) {
        err = 0;
        goto RAPFin;
    }
    printf1("Table construction for years %d-%d and ages %d-%d.\n\n",
        PMYEAR1,PMYEAR2,PMAGE1,PMAGE2);

    ny = PMYEAR2 - PMYEAR1 + 1;
    na = PMAGE2 - PMAGE1 + 1;

    if (alloc_acn(ny * na + 1))                      
        goto RAPFin;
    if (alloc_acm(ny * na + 1))
        goto RAPFin;

    nr = ne = 0;

    for (i = 0; i < NOC; ++i) {
        ts = get_data(its,i);
        tc = get_data(itc,i);
        tf = get_data(itf,i);
        d  = get_data(id ,i);

        for (j = tc; j <= tf; ++j) {
            if (j < PMYEAR1 || j > PMYEAR2)
                continue;
            a = j - ts;
            if (a < PMAGE1 || a > PMAGE2)
                continue;

            t = j - PMYEAR1;
            a -= PMAGE1;

            if (t < 0 || a < 0 || t >= ny || a >= na) {
                printf("t=%d a=%d\n",t,a);
                exit(0);
            }
            AcN[a * ny + t] += 1;
            nr++;
            if (d != 0 && j == tf) {
                AcM[a * ny + t] += 1;
                ne++;
            }
        }
    }

    printf1("                Entries       Cells       Empty\n");
    prnchar('-',47,1);
    t = d = 0;
    a = ny * na;
    for (i = 0; i < a; ++i) {
        if (AcN[i] == 0)
            t++;
        if (AcM[i] == 0)
            d++;
    }
    printf1("Risk table   %10d  %10d  %10d\n",nr,a,t);
    printf1("Event table  %10d  %10d  %10d\n",ne,a,d);
    newline();


    if (PMF1Def) {
        fprintf(PMF1d,"\nRisk table\n");
        fprnchar(PMF1d,' ',PMNFmt + 1,0);
        for (j = 0; j < ny; ++j)
            fprintf(PMF1d,PMNFmtS,PMYEAR1 + j);
        fprintf(PMF1d,"\n");
        for (i = na - 1; i >= 0; --i) {
            fprintf(PMF1d,PMNFmtS,PMAGE1 + i);
            for (j = 0; j < ny; ++j)
                fprintf(PMF1d,PMNFmtS,AcN[i * ny + j]);
            fprintf(PMF1d,"\n");
        }

        fprintf(PMF1d,"\nEvent table\n");
        fprnchar(PMF1d,' ',PMNFmt + 1,0);
        for (j = 0; j < ny; ++j)
            fprintf(PMF1d,PMNFmtS,PMYEAR1 + j);
        fprintf(PMF1d,"\n");
        for (i = na - 1; i >= 0; --i) {
            fprintf(PMF1d,PMNFmtS,PMAGE1 + i);
            for (j = 0; j < ny; ++j)
                fprintf(PMF1d,PMNFmtS,AcM[i * ny + j]);
            fprintf(PMF1d,"\n");
        }
        fprintf(PMF1d,"\nRates\n");
        fprnchar(PMF1d,' ',PMNFmt + 1,0);
        for (j = 0; j < ny; ++j)
            fprintf(PMF1d," %6d",PMYEAR1 + j);
        fprintf(PMF1d,"\n");
        for (i = na - 1; i >= 0; --i) {
            fprintf(PMF1d,PMNFmtS,PMAGE1 + i);
            for (j = 0; j < ny; ++j) {
                r = rap_rate(i,j,PMNC,na,ny);
                fprintf(PMF1d," %6.2f",r);
            }
            fprintf(PMF1d,"\n");
        }
        printf1("Tables written to: %s\n",PMF1dName);
    }

    err = 0;

RAPFin:
    p_clean();  
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  rap_rate(i,j,n)                                                         */

double rap_rate(int i,int j,int n,int na,int ny)
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

            nr += AcN[iii * ny + jjj];
            ne += AcM[iii * ny + jjj];
        }
    }
    if (nr == 0)
        r = -1.0;
    else 
        r = 100.0 * (double)ne / (double)nr;                  
    return(r);
}


