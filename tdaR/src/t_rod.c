/****************************************************************************/
/*  t_rod                                                                   */
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
#include "t_gdd.h"
#include "t_alloc.h"
#include "t_gdat.h"
#include "t_sort.h"
#include "t_gf.h"
#include "t_freq.h"
#include "t_com.h"
#include "tda_context.h"

/*  functions in t_rod.c */

int rod(TDAContext *ctx);
int rod_check(TDAContext *ctx, int n,int m,double *x,int *p);
int rod_icheck(TDAContext *ctx, int m,int *x,int *y);
int rod_dcheck(TDAContext *ctx, int m,double *x,double *y);
int rod_s(TDAContext *ctx, int m,double *x,int *p);
int rod_p(TDAContext *ctx, int opt,int n,int m,double *x,int *p,int *pt,int *pi,int *pj, int *pk,int *r,int rmax,int *dmin);
int rod_d(TDAContext *ctx, int m,double *x,int *p);
int rod_dd(TDAContext *ctx, int opt,int m,double *x,double *p);
int rod_di(TDAContext *ctx, int m,int *x,int *p);
int rod_median(TDAContext *ctx, int n,int m,double *x,int *b,int *p,int *pi,int *pj,   int *r,int max,int *dm);
void rod_ac(TDAContext *ctx, int m,int *r,int *a);  
int rod_acm(TDAContext *ctx, int m,int *ri,int *rj,int *a);  
int rod_acc(TDAContext *ctx, int m,int *a);  
int rod_ar(TDAContext *ctx, int m,int *a,int *r,int *p,int *rp);  
int rod_between(TDAContext *ctx, int m,int *a,int *b,int *r,int i,int j,int *p,int *pi, int *pj,int max,int *rn);

int cro(TDAContext *ctx);
int cro_gen(TDAContext *ctx, int m,int *p,int *q,int *s,int *c,int *d,int *c1,int *d1, int first,int *nt);
int ro_b(TDAContext *ctx, int m,int *a,int *b,int *c);
int ro_bd(TDAContext *ctx, int m,double *a,double *b,double *c);


/* ------------------------------------------------------------------------ */
/*  rod         Rank order data                                             */
/*                                                                          */
/*              rod(                                                        */
/*                  opt=...,        option, def. 1                          */
/*                                  1 = create standard format rank orders  */
/*                                  2 = create distance matrix              */
/*                                  3 = find median/mean rank orders        */
/*                                  4 = create graph                        */
/*                                  5 = create graph                        */
/*                                  6 = create graph                        */
/*                                  7 = find rank orders in between         */
/*                                  8 = tabulate rank orders                */
/*                  s=...,          option, def. 0                          */
/*                                  0 = median rank orders                  */
/*                                  1 = mean rank orders                    */
/*                  alg=...,        algorithm for option 3, def. 1          */
/*                                  1 = complete enumeration                */
/*                                  2 = for rank orders without ties        */
/*                                      ONLY median rank orders             */
/*                                  3 = experimental                        */
/*                  max=...,        max number of multiple rank orders,     */
/*                                  def. 1000                               */
/*                  w=...,          case weight variable for option 3       */
/*                                  (currently only with algorithm 1)       */
/*                  df=...,         output file                             */
/*                  fmt=...,        print format for df, def. 0.0           */
/*                  pcf=...,        plot file for neato (options 4 and 5)   */
/*              ) = X1,...,Xm;      variables.                              */
/*                                                                          */
/*  Values of the variables are interpreted as follows. If the values       */
/*  for case k are: (xk1,...,xkm), this means:                              */
/*  a)  alternative i is preferred to alternative j if xki > xkj            */
/*  b)  indifference if xki = xkj.                                          */
/*                                                                          */
/*  Option 1 is only done if there is an output file defined with the       */
/*      df option.                                                          */
/*                                                                          */
/*  Option 2 writes results into standard output or, alternatively, into    */
/*      the output file defined with the df option.                         */
/*                                                                          */
/*  Option 3 depends on the s parameter (also alg and max).                 */
/*      Alg 1 uses complete enumration of all rank orders including ties.   */
/*      Alg 2 is a two-step algorithm only for median rank orders that      */
/*            assumes that the input rank orders have no ties. It is also   */
/*            assumed that the input data are integer-valued.               */
/*      Alg 3 is experimental and should not be used.                       */
/*                                                                          */
/*  Option 4 creates a graph where two rank orders are connected when       */
/*      there is no other rank order in between.                            */
/*                                                                          */
/*  Option 5 creates a graph where two rank orders are connected when       */
/*      there is at least one rank order in between.                        */
/*                                                                          */  
/*  Option 6 creates a graph where each triple (r,r',r'') is connected      */
/*      when r' is between r and r''.                                       */
/*                                                                          */  
/*      For these options: if df is used then the edge list is written,     */
/*      if pcf is used then a plot file for neato is used.                  */
/*                                                                          */  
/*  Option 7 tries to find all rank orders that are in between of at least  */
/*      two rank orders in the input data. The maximum number of rank       */
/*      orders must be specified with the max parameter (default 1000).     */
/*      If the df option is used the rank orders are written into the       */
/*      output file.                                                        */
/*                                                                          */
/*  Option 8 creates an aggregates table where each rank order appears      */
/*      only once. An additional column contains the frequencies. If        */
/*      the df option is used the table is written into an output file.     */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int rod(TDAContext *ctx)
{
    register int i,j,k,l;
    int err,n,m,r,rn,d,di,dm,nt,first,ofs,mc,nw,max;
    double wt,*xi;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Rank order data. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 3,4,1))       /* get parameters */
        goto RODFin;

    if (ctx->PMOPT > 8)
        ctx->PMOPT = 8;
    if (ctx->PMALG < 1 || ctx->PMALG > 3)
        ctx->PMALG = 1;
    if (ctx->PMMax < 1)
        ctx->PMMax = 1000;
    if (ctx->PMS != 0)
        ctx->PMS = 1;

    n = ctx->NOC;
    m = ctx->PMNV;                       /* number of variables */

    if (m < 2) {
        printf1(ctx, "Error: need at least two variables.\n");
        goto RODFin;
    }
    if (alloc_acx(ctx, n * m + 1))           /* input data */
        goto RODFin;

    if (ctx->PMWVar >= 0) {
        if (alloc_acw(ctx, n + 1))           /* weights */
            goto RODFin;
    }

    for (i = 0; i < n; ++i) {           /* get data into AcX and AcW */
        for (j = 0; j < m; ++j)  
            ctx->AcX[i * m + j] = get_data(ctx, ctx->PMVIdx[j],i);

        if (ctx->PMWVar >= 0) {
            wt = get_data(ctx, ctx->PMWVar,i);
            if (wt < 0.0) {
                printf1(ctx, "Error: found negative weight in case %d\n",i + 1);
                goto RODFin;
            }
            ctx->AcW[i] = wt;
        }
    }
                     
       
    for (i = 0; i < n; ++i) {   
        for (j = 0; j < m; ++j)  
            tda_out("%g ",ctx->AcX[i * m + j]);

        if (ctx->PMWVar >= 0)  
            tda_out("%g ",ctx->AcW[i]);
        newline(ctx);
    }
    newline(ctx);

    if (ctx->PMOPT == 1) {             
        printf1(ctx, "Creating standard format rank orders.\n");
        if (ctx->PMF1Def == 0) {
            printf1(ctx, "Need output file (df option).\n");
            goto RODFin;
        }
        if (alloc_ack(ctx, m + 1))                                           
            goto RODFin;

        rn = 1;
        for (i = 0; i < n; ++i) {
            r = rod_s(ctx, m,ctx->AcX + i * m,ctx->AcK);
            rn = imax(ctx, r,rn);
        }
        for (i = 0; i < n; ++i) {
            for (j = 0; j < m; ++j) 
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[i * m + j]);
            fprintf(ctx->PMF1d,"\n");
        }
        printf1(ctx, "%d records written to: %s\n",n,ctx->PMF1dName);
        printf1(ctx, "Max size of tie blocks: %d\n",rn);
        goto RODFin1;
    }
    else if (ctx->PMOPT == 2) {      /* create distance matrix */

        printf1(ctx, "Distance matrix (s = %d).\n",ctx->PMS);

        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; ++j) {
                r = rod_dd(ctx, ctx->PMS,m,ctx->AcX + i * m,ctx->AcX + j * m);

                if (ctx->PMF1Def)  
                    fprintf(ctx->PMF1d,"%3d ",r);
                else
                    printf1(ctx, "%3d ",r);
            }
            if (ctx->PMF1Def)  
                fprintf(ctx->PMF1d,"\n");
            else
                newline(ctx);
        }
        if (ctx->PMF1Def)  
            printf1(ctx, "%d records written to: %s\n",n,ctx->PMF1dName);

        goto RODFin1;
    }
    else  if (ctx->PMOPT == 3) {     /* find mean or median */

        if (ctx->PMALG == 2)
            ctx->PMS = 0;

        printf1(ctx, "Searching for ");
        if (ctx->PMS == 0)
            printf1(ctx, "median");
        else
            printf1(ctx, "mean");
        printf1(ctx, " rank orders (alg=%d).\n",ctx->PMALG);

        if (ctx->PMWVar >= 0)  
            printf1(ctx, "Using weights defined by: %s\n",ctx->VName[ctx->PMWVar]);

        if (ctx->PMALG == 1) {
            mc = ofs = 0;

            if (alloc_acn(ctx, m + 1))   
                goto RODFin;
            if (alloc_acm(ctx, m + 1))   
                goto RODFin;
            if (alloc_ack(ctx, m + 1))   
                goto RODFin;
            if (alloc_aci(ctx, 2 * m + 3))   
                goto RODFin;
            if (alloc_acj(ctx, 2 * m + 3))   
                goto RODFin;
            if (alloc_acr(ctx, ctx->PMMax * m + 1))       
                goto RODFin;

            rn = 0;
            dm = 0;
            first = 1;
            while (cro_gen(ctx, m,ctx->AcN,ctx->AcM,ctx->AcK,ctx->AcI,ctx->AcJ,ctx->AcI + m + 1,ctx->AcJ + m + 1,first,&nt)) {
                mc++;
                di = 0;
                for (i = 0; i < n; ++i) {
                    d = rod_d(ctx, m,ctx->AcX + i * m,ctx->AcN);
                    if (ctx->PMS)
                        d *= d;
                    if (ctx->PMWVar >= 0)  
                        d *= (int)ctx->AcW[i];
                    di += d;
                }
                if (first) {
                    dm = di;
                    for (i = 0; i < m; ++i)
                        ctx->AcR[i] = ctx->AcN[i];
                    rn = 1;
                }
                else if (di < dm) {
                    dm = di;
                    for (i = 0; i < m; ++i)
                        ctx->AcR[i] = ctx->AcN[i];
                    rn = 1;
                }
                else if (di == dm) {
                    if (rn >= ctx->PMMax) {
                        rn = -1;
                        break;
                    }                
                    for (i = 0; i < m; ++i)
                        ctx->AcR[rn * m + i] = ctx->AcN[i];
                    rn++;
                }
                first = 0;
            }
        }
        else if (ctx->PMALG == 2) {

            mc = ofs = 0;

            if (alloc_acn(ctx, m * m + 1))       
                goto RODFin;
            if (alloc_aci(ctx, m + 1))       
                goto RODFin;
            if (alloc_acj(ctx, m + 1))       
                goto RODFin;
            if (alloc_ack(ctx, m + 1))       
                goto RODFin;
            if (alloc_acr(ctx, ctx->PMMax * m + 1))       
                goto RODFin;

            rn = rod_median(ctx, n,m,ctx->AcX,ctx->AcN,ctx->AcI,ctx->AcJ,ctx->AcK,ctx->AcR,ctx->PMMax,&dm);
        }
        else {                          /* PMALG == 3, clamped above */
            mc = 0;
            ofs = 1;

            if (alloc_ack(ctx, m + 1))       
                goto RODFin;
            if (alloc_aci(ctx, m + 1))       
                goto RODFin;
            if (alloc_acj(ctx, m + 1))       
                goto RODFin;
            if (alloc_acn(ctx, m + 1))       
                goto RODFin;
            if (alloc_acm(ctx, ctx->PMMax * m + 1))       
                goto RODFin;
            if (alloc_acr(ctx, ctx->PMMax * m + 1))       
                goto RODFin;

            rn = rod_p(ctx, ctx->PMS,n,m,ctx->AcX,ctx->AcM,ctx->AcK,ctx->AcI,ctx->AcJ,ctx->AcN,ctx->AcR,ctx->PMMax,&dm);

        }
        if (rn < 0) {
            printf1(ctx, "Error: exceeded max number of rank orders (max=%d).\n",ctx->PMMax);
            goto RODFin;
        }
        printf1(ctx, "Minimal sum of distances (s = %d, mc=%d): %d\n",ctx->PMS,mc,dm);
        printf1(ctx, "Number of rank orders: %d\n\n",rn);

        for (i = 0; i < rn; ++i) {
            printf1(ctx, "%3d : ",i + 1);
            for (j = 0; j < m; ++j)
                printf1(ctx, "%2d ",ctx->AcR[i * m + j] + ofs);
            newline(ctx);
        }
        goto RODFin1;
    }
    else  if (ctx->PMOPT >= 4 && ctx->PMOPT <= 6) {     /* create graph */

        if (alloc_acd(ctx, n * n + 1))   
            goto RODFin;

        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; ++j) {
                if (i == j)
                    continue;

                rn = 1;
                for (k = 0; k < n; ++k) {
                    if (k == i || k == j)
                        continue;
                    r = ro_bd(ctx, m,ctx->AcX + i * m,ctx->AcX + j * m,ctx->AcX + k * m);
                    if (r) {
                        rn = 0;
                        if (ctx->PMOPT <= 5)
                            break;

                        if (ctx->PMOPT == 6) {
                            ctx->AcD[i * n + k] = 1;
                            ctx->AcD[k * n + i] = 1;
                            ctx->AcD[k * n + j] = 1;
                            ctx->AcD[j * n + k] = 1;
                        }
                    }
                }
                if (rn && ctx->PMOPT <= 5) {
                    ctx->AcD[i * n + j] = 1;
                    ctx->AcD[j * n + i] = 1;
                }
            }
        }
        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; ++j)
                printf1(ctx, "%d ",(int)ctx->AcD[i * n + j]);
            newline(ctx);
        }

        if (ctx->PMPCFDef || ctx->PMF1Def) {
            k = 0;
            if (ctx->PMPCFDef)                 
                fprintf(ctx->PMPCFd,"graph g {\nsize=\"5,5\";\n");

            for (i = 1; i < n; ++i) {
                for (j = 0; j < i; ++j) {

                    if (((ctx->PMOPT == 4 || ctx->PMOPT == 6) && ctx->AcD[i * n + j] != 0) ||
                        (ctx->PMOPT == 5 && ctx->AcD[i * n + j] == 0)) {

                        if (ctx->PMPCFDef)   
                            fprintf(ctx->PMPCFd,"%2d -- %2d [len=2];\n",i + 1, j + 1);

                        if (ctx->PMF1Def)  
                            fprintf(ctx->PMF1d,"%2d %2d %2d\n",i + 1, j + 1,1);
                        k++;
                    }   
                }
            }
            if (ctx->PMF1Def)  
                printf1(ctx, "%d records written to: %s\n",k,ctx->PMF1dName);
            if (ctx->PMPCFDef) { 
                fprintf(ctx->PMPCFd,"}\n");
                printf1(ctx, "%d records written to: %s\n",k + 3,ctx->PMPCFName);
            }
        }
        goto RODFin1;
    }
    else  if (ctx->PMOPT == 7) {     /* find rank orders in between  */

        if (n < 2) {
            printf1(ctx, "Error: need at least two rank orders.\n");
            goto RODFin;
        }
        max = m * m * m * m;

        if (alloc_acr(ctx, m + 1))   
            goto RODFin;
        if (alloc_acs(ctx, m + 1))   
            goto RODFin;
        if (alloc_acn(ctx, m * m + 1))   
            goto RODFin;
        if (alloc_acm(ctx, m * m + 1))   
            goto RODFin;
        if (alloc_ack(ctx, m * m + 1))   
            goto RODFin;
        if (alloc_aci(ctx, 2 * m + 4))   
            goto RODFin;
        if (alloc_acj(ctx, ctx->PMMax * m + 1))   
            goto RODFin;
        if (alloc_acl(ctx, max * m + 1))   
            goto RODFin;

        rn = 0;
        for (i = 1; i < n; ++i) {
                 
            for (k = 0; k < m; ++k)  
                ctx->AcR[k] = (int)ctx->AcX[i * m + k];

            for (j = 0; j < i; ++j) {
                 
                for (k = 0; k < m; ++k)  
                    ctx->AcS[k] = (int)ctx->AcX[j * m + k];

                rod_ac(ctx, m,ctx->AcR,ctx->AcN);
                rod_ac(ctx, m,ctx->AcS,ctx->AcM);

                d = 0;
                if (rod_between(ctx, m,ctx->AcN,ctx->AcM,ctx->AcK,0,0,ctx->AcL,ctx->AcI,ctx->AcI + m + 1,max,&d)) {
                    printf1(ctx, "Error: exceeded max number of rank orders.\n");
                    goto RODFin;
                }
                for (k = 0; k < d; ++k) {
                    nt = 1;
                    for (l = 0; l < rn; ++l) {
                        if (rod_icheck(ctx, m,ctx->AcL + k * m,ctx->AcJ + l * m) == 0) {
                            nt = 0;
                            break;
                        }
                    }
                    if (nt) {
                        if (rn >= ctx->PMMax) {
                            printf1(ctx, "Error: exceeded max number of rank orders (max=%d).\n",ctx->PMMax);
                            goto RODFin;
                        }
                        for (l = 0; l < m; ++l)
                            ctx->AcJ[rn * m + l] = ctx->AcL[k * m + l];
                        rn++;
                    }
                }
            }
        }
        printf1(ctx, "Found %d rank orders.\n",rn);
        if (ctx->PMF1Def) {
            for (i = 0; i < rn; ++i) {
                for (j = 0; j < m; ++j) 
                    fprintf(ctx->PMF1d,"%2d ",ctx->AcJ[i * m + j]);
                fprintf(ctx->PMF1d,"\n");
            }
            printf1(ctx, "%d records written to: %s\n",rn,ctx->PMF1dName);
        }
        goto RODFin1;
    }
    else  if (ctx->PMOPT == 8) {     /* standard integer format */

        printf1(ctx, "Tabulate rank orders.\n");

        if (alloc_acy(ctx, n * m + 1))   
            goto RODFin;
        if (alloc_acn(ctx, n + 1))   
            goto RODFin;

        for (j = 0; j < m; ++j)
            ctx->AcY[j] = ctx->AcX[j];
        if (ctx->PMWVar >= 0)  
            ctx->AcN[0] = (int)ctx->AcW[0];
        else 
            ctx->AcN[0] = 1;

        rn = 1;
        for (i = 1; i < n; ++i) {
            r = -1;
            xi = ctx->AcX + i * m;
            for (j = 0; j < rn; ++j) {
                if (rod_dcheck(ctx, m,xi,ctx->AcY + j * m) == 0) {
                    r = j;
                    break;
                }   
            }   
            if (r >= 0) {
                if (ctx->PMWVar >= 0)  
                    ctx->AcN[r] += (int)ctx->AcW[i];
                else 
                    ctx->AcN[r] += 1;
            }
            else {
                for (j = 0; j < m; ++j)
                    ctx->AcY[rn * m + j] = xi[j];
                if (ctx->PMWVar >= 0)  
                    ctx->AcN[rn] = (int)ctx->AcW[i];
                else 
                    ctx->AcN[rn] = 1;
                rn++;
            }
        }
        nw = 0;
        for (i = 0; i < rn; ++i)
            nw += ctx->AcN[i];

        printf1(ctx, "Found %d different rank orders.\n",rn);
        printf1(ctx, "Sum of weights: %d\n",nw);

        if (ctx->PMF1Def) {
            for (i = 0; i < rn; ++i) {
                xi = ctx->AcY + i * m;
                for (j = 0; j < m; ++j)
                    fprintf(ctx->PMF1d,"%g ",xi[j]);
                fprintf(ctx->PMF1d,"%5d\n",ctx->AcN[i]);
            }
            printf1(ctx, "%d records written to: %s\n",rn,ctx->PMF1dName);
        }
        goto RODFin1;
    }

RODFin1:
    newline(ctx);
    err = 0;

RODFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  rod__check()    return 1 if p in x.                                     */
       
/* ------------------------------------------------------------------------ */
/*  rod__icheck()   return 1 if x and y are different rank orders.          */
       
int rod_icheck(TDAContext *ctx, int m,int *x,int *y)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j;

    for (i = 1; i < m; ++i) {
        for (j = 0; j < i; ++j) {
            if (x[j] < x[i]) {
                if (y[j] >= y[i])
                    return(1);
            }
            else if (x[j] > x[i]) {
                if (y[j] <= y[i])
                    return(1);
            }
            else if (x[j] == x[i]) {
                if (y[j] != y[i])
                    return(1);
            }
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  rod__dcheck()   return 1 if x and y are different rank orders.          */
       
int rod_dcheck(TDAContext *ctx, int m,double *x,double *y)
{
    register int i,j;

    for (i = 1; i < m; ++i) {
        for (j = 0; j < i; ++j) {
            if (x[j] < x[i] - ctx->EPSI1) {
                if (y[j] >= y[i] - ctx->EPSI1)
                    return(1);
            }
            else if (x[j] > x[i] + ctx->EPSI1) {
                if (y[j] <= y[i] + ctx->EPSI1)
                    return(1);
            }
            else if (fabs(x[j] - x[i]) <= ctx->EPSI1) {
                if (fabs(y[j] - y[i]) > ctx->EPSI1)
                    return(1);
            }
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  rod_s(m,x,p)    change x[i] (i=0,m-1) into standard rank order format.  */
/*                                                                          */
/*  Return max size of tie groups.                                          */
       
int rod_s(TDAContext *ctx, int m,double *x,int *p)
{
    register int i,j,n;       
    int r;
    double xa;

    r = 1;
    sortdp(ctx, m,x,p);
    i = 0;
    while (i < m) {
        xa = x[p[i]];
        n = 1;
        for (j = i + 1; j < m; ++j) {
            if (x[p[j]] != xa)
                break;
            n++;
        }
        xa = (double)i + (double)(n + 1) / 2.0;
        for (j = 1; j <= n; ++j) {
            x[p[i]] = xa;
            i++;
        }
        r = imax(ctx, r,n);
    }
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  rod_p(opt)  calculate median (opt = 0) or mean (opt = 1) rank order,    */
/*              returned in r[].                                            */
/*              min distance in dm. rmax is max number of multiple          */
/*              rank orders (rows in r[] and p[]).                          */
/*                                                                          */
/*  Return: number of rank orders in r[], or -1 if insufficient rmax        */
       
int rod_p(TDAContext *ctx, int opt,int n,int m,double *x,int *p,int *pt,int *pi,int *pj, int *pk,int *r,int rmax,int *dm)
{
    register int i,j,k,l,q;
    int first,d,di,npi,rn0,rn1,ii,a,b;

    /* first step: check with complete ties  */

    for (i = 0; i < m; ++i)
        pk[i] = 0;

    *dm = 0;
    for (i = 0; i < n; ++i) {
        di = rod_d(ctx, m,x + i * m,pk);
        if (opt)
            di *= di;
        *dm += di;
    }

/*  tda_out("FIRST dm=%d  \n",*dm);      */

    for (i = 0; i < m; ++i)     /* save in p[] */
        p[i] = pt[i];
    rn0 = 1;

    /* second step: check all permutations without ties, save optimal
       permutations in p[] (i = 0,...,rn0-1) */

    first = 1;
    while (perm(ctx, m,pt,pi,pj,first)) {

        d = 0;
        for (i = 0; i < n; ++i) {
            di = rod_d(ctx, m,x + i * m,pt);
            if (opt)
                di *= di;
            d += di;
            if (d > *dm)
                break;
        }
                            
        if (d < *dm) {
            *dm = d;
            for (i = 0; i < m; ++i)
                p[i] = pt[i];
            rn0 = 1;
        }
        else if (d == *dm) {
            if (rn0 >= rmax)  
                return(-1);

            for (i = 0; i < m; ++i)
                p[rn0 * m + i] = pt[i];
            rn0++;
        }
        first = 0;
    }
      
tda_out("\nAfter second step. rn0=%d dm=%d \n",rn0,*dm);
for (i = 0; i < rn0; ++i) {
    for (j = 0; j < m; ++j)
        tda_out("%d ",p[i * m + j]);
    newline(ctx);
}

tda_out("dm=%d\n",*dm);
       


    /*  third step: check all possible tie groups */

    for (i = 0; i < rn0 * m; ++i)   /* copy to r[] */
        r[i] = p[i];
    rn1 = rn0;

    for (ii = 0; ii < rn0; ++ii) {  /* for all permutations from second step */

        for (k = 2; k < m; ++k) {
tda_out("\nii=%d k=%d\n",ii,k);
            first = 1;
            while (partit1(ctx, m,k,pi,pj,first) == 0) {
      
                tda_out("pi: ");
                for (i = 1; i <= k; ++i)
                    tda_out("%d ",pi[i]);
                newline(ctx);
        


                l = 0;
                for (i = 1; i <= k; ++i) {
                    npi = pi[i];
                    q = l;
                    for (j = 0; j < npi; ++j)
                        pt[l++] = q;
                }

                       
                tda_out("pt: ");
                for (i = 0; i < m; ++i)
                    tda_out("%d ",pt[i]);
                newline(ctx);
                       

       
                for (i = 0; i < m; ++i)  
                    pk[i] = pt[p[ii * m + i]];
        
                         
                tda_out("pk: ");
                for (i = 0; i < m; ++i)
                    tda_out("%d ",pk[i]);
                newline(ctx);
        
                tda_out("xxx: ");
                for (i = 0; i < m; ++i)
                    tda_out("%g ",x[i]);
                newline(ctx);
                         

                d = 0;
                for (i = 0; i < n; ++i) {
                    di = rod_d(ctx, m,x + i * m,pk);
                    if (opt)
                        di *= di;
                    d += di;
                    if (d > *dm)
                        break;
                }


                tda_out("d========%d dm=%d  rn1=%d \n",d,*dm,rn1);   

                if (d < *dm) {
                    *dm = d;
                    for (i = 0; i < m; ++i)
                        r[i] = pk[i];
                    rn1 = 1;
                }
                else if (d == *dm) {
                    if (rn0 >= rmax)  
                        return(-1);

                    b = 1;
                    for (i = 0; i < rn1; ++i) {
                        a = 1;
                        for (j = 0; j < m; ++j) {
                            if (pk[j] != r[i * m + j]) {
                                a = 0;
                                break;
                            }
                        }
                        if (a == 1) {
                            b = 0;
                            break;
                        }
                    }
                    if (b) {
                        for (i = 0; i < m; ++i)
                            r[rn1 * m + i] = pk[i];
                        rn1++;
                    }
                }
                first = 0;
            }
        }
    }
    return(rn1);
}

/* ------------------------------------------------------------------------ */
/*  rod_d(m,x,p)    return distance between x[] and p[].                    */
/*                                                                          */
       
int rod_d(TDAContext *ctx, int m,double *x,int *p)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j;
    int d;
       
    d = 0;
    for (i = 1; i < m; ++i) {
        for (j = 0; j < i; ++j) {
            if (x[i] == x[j]) {
                if (p[i] != p[j])
                    d++;           
            }
            else if (x[i] < x[j]) {
                if (p[i] == p[j])
                    d++;
                else if (p[i] > p[j])
                    d += 2;
            }
            else {
                if (p[i] == p[j])
                    d++;
                else if (p[i] < p[j])
                    d += 2;
            }
        }
    }
    return(d);
}

/* ------------------------------------------------------------------------ */
/*  rod_dd(opt,m,x,p)   return distance between x[] and p[].                */
/*                      if opt=0 absolute, else square.                     */
       
int rod_dd(TDAContext *ctx, int opt,int m,double *x,double *p)
{
    register int i,j;
    int d,di;
       
    di = 0;
    for (i = 1; i < m; ++i) {
        for (j = 0; j < i; ++j) {
            d = 0;
            if (fabs(x[i] - x[j]) <= ctx->EPSI1) {
                if (fabs(p[i] - p[j]) > ctx->EPSI1)
                    d++;           
            }
            else if (x[i] < x[j] - ctx->EPSI1) {
                if (fabs(p[i] - p[j]) <= ctx->EPSI1)
                    d++;
                else if (p[i] > p[j] + ctx->EPSI1)
                    d += 2;
            }
            else {
                if (fabs(p[i] - p[j]) <= ctx->EPSI1)
                    d++;
                else if (p[i] < p[j] - ctx->EPSI1)
                    d += 2;
            }
            if (opt)
                d *= d;
            di += d;
        }
    }
    return(di);
}

/* ------------------------------------------------------------------------ */
/*  rod_di(m,x,p)    return distance between x[] and p[].                   */
/*                                                                          */
       
int rod_di(TDAContext *ctx, int m,int *x,int *p)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j;
    int d;
       
    d = 0;
    for (i = 1; i < m; ++i) {
        for (j = 0; j < i; ++j) {
            if (x[i] == x[j]) {
                if (p[i] != p[j])
                    d++;           
            }
            else if (x[i] < x[j]) {
                if (p[i] == p[j])
                    d++;
                else if (p[i] > p[j])
                    d += 2;
            }
            else {
                if (p[i] == p[j])
                    d++;
                else if (p[i] < p[j])
                    d += 2;
            }
        }
    }
    return(d);
}

/* ------------------------------------------------------------------------ */
/*  rod_median(n,m,x,b,p,pi,pj,r,max,dm)                                    */
/*                                                                          */
/*  n    = number of rank orders in x[]                                     */
/*  m    = size of rank orders                                              */
/*  x[]  = n,m matrix with rank orders                                      */
/*  b[]  = m,m matrix                                                       */
/*  p[]  = m vector for permutations                                        */
/*  pi[] = m vector                                                         */
/*  pj[] = m vector                                                         */
/*  r[]  = array for maximal max median rank orders                         */
/*  max  = max number of rank orders in r[]                                 */
/*  dm   = best value of criterion                                          */
/*                                                                          */
/*  Return: number of rank orders in r[], or -1 if max exceeded.            */

int rod_median(TDAContext *ctx, int n,int m,double *x,int *b,int *p,int *pi,int *pj, int *r,int max,int *dm)
{
    register int i,j,k,ii;
    int bjk,first,d,rn,rn0,*ri,*rii;
    double *xi;

    /* create b matrix */

    for (j = 0; j < m; ++j) {
        for (k = 0; k < m; ++k) {
            bjk = 0;
            if (j != k) {
                for (i = 0; i < n; ++i) {
                    xi = x + i * m;
                    if (xi[j] < xi[k] - ctx->EPSI1)
                        bjk++;
                    else if (xi[j] > xi[k] + ctx->EPSI1)
                        bjk--;
                }
            }
            b[j * m + k] = bjk;
        }
    }
    rn = 0;
    *dm = 0;
    first = 1;
    while (perm(ctx, m,p,pi,pj,first)) {
        first = 0;
        d = 0;
        for (j = 1; j < m; ++j) {
            for (k = 0; k < j; ++k) {
                if (p[j] < p[k])
                    d += b[j * m + k]; 
                else if (p[j] > p[k])
                    d -= b[j * m + k];
            }
        }
        if (d > *dm) {
            *dm = d;
            for (j = 0; j < m; ++j) 
                r[j] = p[j] + 1;
            rn = 1;
        }
        else if (d == *dm) {
            if (rn >= max)
                return(-1);
            for (j = 0; j < m; ++j) 
                r[rn * m + j] = p[j] + 1;
            rn++;
        }
    }

    /* find all mean rank orders */

    rn0 = rn;

    for (i = 1; i < rn; ++i) {
        ri = r + i * m;
        for (ii = 0; ii < i; ++ii) {
            rii = r + ii * m;
      
            if (rod_acm(ctx, m,ri,rii,b))  
                continue;
            if (rod_acc(ctx, m,b))  
                continue;
            if (rod_ar(ctx, m,b,p,pi,pj))             
                continue;

            /* check whether already there */

            first = 1;
            for (k = rn0; k < rn; ++k) {
                if (rod_icheck(ctx, m,p,r + k * m) == 0) {
                    first = 0;
                    break;
                }
            }
            if (first) {
                if (rn >= max)
                    return(-1);
                for (k = 0; k < m; ++k) 
                    r[rn * m + k] = p[k];
                rn++;
            }
        }
    }
    return(rn);
}

/* ------------------------------------------------------------------------ */
/*  rod_ac(m,r,a)                                                           */
/*                                                                          */
/*  For the rank order r[], create a rank order matrix in a[].              */

void rod_ac(TDAContext *ctx, int m,int *r,int *a)   
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j;
    int k;
       
    for (i = 1; i < m; ++i) {
        for (j = 0; j < i; ++j) {
            if (r[i] < r[j])
                k = 1;
            else if (r[j] < r[i])
                k = -1;
            else
                k = 0;

            a[i * m + j] = k;
            a[j * m + i] = -k;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  rod_acm(m,ri,rj,a)                                                      */
/*                                                                          */
/*  Create rank order matrix for the mean of ri[] and rj[] in a[]. If the   */
/*  mean cannot be generated return 1, otherwise 0.                         */

int rod_acm(TDAContext *ctx, int m,int *ri,int *rj,int *a)   
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j,k;

    for (i = 1; i < m; ++i) {
        for (j = 0; j < i; ++j) {
            if (ri[i] < ri[j])
                k = 1;
            else if (ri[j] < ri[i])
                k = -1;
            else
                k = 0;

            a[i * m + j] = k;
        }
    }

    for (i = 1; i < m; ++i) {
        for (j = 0; j < i; ++j) {
            k = a[i * m + j];
            if (rj[i] < rj[j])  
                k += 1;        
            else if (rj[j] < rj[i])
                k -= 1;      

            if (k != 0 && k != 2 && k != -2)
                return(1);
            k /= 2;
            a[i * m + j] = k;
            a[j * m + i] = -k;
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  rod_acc(m,a)                                                            */
/*                                                                          */
/*  a[] is a matrix containing -1, 0, +1. Return 0 if a[] is a rank order   */
/*  matrix, otherwise return 1.                                             */

int rod_acc(TDAContext *ctx, int m,int *a)   
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j,k;
    int aij,aik,akj;

    for (i = 0; i < m; ++i) {
        for (j = 0; j < m; ++j) {
            if (i == j) 
                continue;

            aij = a[i * m + j];
            if (aij != -a[j * m + i]) {
                tda_out("NOT0\n");
                return(1);
            }
            for (k = 0; k < m; ++k) {
                if (k == i || k == j)
                    continue;
                aik = a[i * m + k];
                akj = a[k * m + j];

                if (aik >= 0 && akj >= 0) {
                    if (aij < 0)  
                        return(1);
                    if (aik == 0 && akj == 0 && aij != 0)  
                        return(1);
                    if (aij == 0 && (aik != 0 || akj != 0))  
                        return(1);
                }
            }
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  rod_ar(m,a,r,p,rp)                                                      */
/*                                                                          */
/*  a[] is a (m,m) rank order matrix checked by rod_acc(). Try to create    */
/*  the corresponding rank order in r[]. If possible, return 0, otherwise 1 */

int rod_ar(TDAContext *ctx, int m,int *a,int *r,int *p,int *rp)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j,k,l;
    int pmax;
             
    for (i = 0; i < m; ++i)     /* record permutations */
        p[i] = i;

    pmax = m * m + 1;

AGAIN:          
    for (i = 0; i < m; ++i) {
        for (j = i + 1; j < m; ++j) {
            if (a[p[i] * m + p[j]] < 0) {
                k = p[i];
                p[i] = p[j];
                p[j] = k;
                if (--pmax < 0) {
                    tda_err("FATAL ERROR 1 in ROD_AR.\n");
                    exit(0);
                }
                goto AGAIN;
            }
        }
    }
    rp[0] = pmax = 1;            /* begin with first row */
    for (j = 1; j < m; ++j) {
        k = a[p[0] * m + p[j]];
        if (k > 0)
            rp[j] = ++pmax;
        else if (k < 0) {
            rp[j] = rp[j - 1];
            for (l = 0; l < j; ++l)
                rp[l] += 1;
            pmax++;
        }
        else
            rp[j] = rp[0];
    }

    /* now check all other rows */

    for (i = 1; i < m; ++i) {
        for (j = i + 1; j < m; ++j) {
            k = a[p[i] * m + p[j]];          
            if (k == 1) {
                if (rp[i] >= rp[j]) {
                    tda_err("FATAL ERROR 2 in ROD_AR.\n");
                    exit(0);
                }
            }
            else if (k == 0) {
                if (rp[i] != rp[j])  
                    rp[j] = rp[i];
            }
            else {
                tda_err("FATAL ERROR 3 in ROD_AR.\n");
                exit(0);
            }
        }
    }
    for (i = 0; i < m; ++i)
        r[p[i]] = rp[i];

    return(0);
}

/* ------------------------------------------------------------------------ */
/*  rod_between(m,a,b,r,i,j,p,pi,pj,n,max,rn)                               */
/*                                                                          */
/*  Find all rank orders in between a[] and b[]. a[] and b[] are assumed    */
/*  to be valid rank order matrices. r[] is m,m matrix used to create the   */
/*  new rank orders. pi[] and pj[] are m vectors. p is max times m matrix.  */
/*  max is maximum number of rank orders that can be stored. p[] will be    */
/*  used to store the new rank orders. rn is current number of rank         */
/*  orders in p[].                                                          */
/*                                                                          */
/*  The function returns -1 if max is insufficient, otherwise 0.            */

int rod_between(TDAContext *ctx, int m,int *a,int *b,int *r,int i,int j,int *p,int *pi, int *pj,int max,int *rn)
{
    int aij,bij,rij;

    if (++j >= m) {
        if (++i == m - 1) {
            for (i = 1; i < m; ++i) {
                for (j = 0; j < i; ++j)
                    r[i * m + j] = -r[j * m + i];
            }

            /* check for new rank order */

            if (rod_acc(ctx, m,r) == 0) {        /* if rank order matrix */
                if (*rn >= max)
                    return(-1);

                rod_ar(ctx, m,r,p + *rn * m,pi,pj);
                *rn += 1;
            }
            return(0);
        }
        j = i + 1;
    }

    aij = a[i * m + j];
    bij = b[i * m + j];

    if (aij <= bij) {
        for (rij = aij; rij <= bij; ++rij) {
            r[i * m + j] = rij;
            if (rod_between(ctx, m,a,b,r,i,j,p,pi,pj,max,rn))
                return(-1);
        }
    }
    else  {
        for (rij = aij; rij >= bij; --rij) {
            r[i * m + j] = rij;
            if (rod_between(ctx, m,a,b,r,i,j,p,pi,pj,max,rn))
                return(-1);
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  cro         Create all rank orders of size m.                           */
/*  ##                                                                      */
/*              cro(                                                        */
/*                  opt=...,    option, def. 1                              */
/*                              1 = create rank orders of size m            */
/*                              2 = find all rank orders between two        */
/*                                  extremes of size m.                     */
/*                              3 = create edge list for all rank orders    */
/*                                  with distance 1                         */
/*                  m=...,      def. 2                                      */
/*              ) = fname;                                                  */
/*                                                                          */
/*  The output file contains m + 2 columns. First column contains record    */  
/*  number, second column contains number of tie groups, then follows       */
/*  the rank order.                                                         */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int cro(TDAContext *ctx)
{
    register int i,j,k; 
    int err,n,nn,nt,first,d;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Creating rank orders. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto CROFin;

    if (ctx->PMOPT > 3)
        ctx->PMOPT = 3;
    if (ctx->PMM < 2)
        ctx->PMM = 2;

    if (alloc_acn(ctx, ctx->PMM + 1))   
        goto CROFin;
    if (alloc_acm(ctx, ctx->PMM + 1))   
        goto CROFin;
    if (alloc_ack(ctx, ctx->PMM + 1))   
        goto CROFin;
    if (alloc_aci(ctx, 2 * ctx->PMM + 3))   
        goto CROFin;
    if (alloc_acj(ctx, 2 * ctx->PMM + 3))   
        goto CROFin;

    n = 0;
    if (ctx->PMOPT == 1) {           /* create rank orders */
        first = 1;
        while (cro_gen(ctx, ctx->PMM,ctx->AcN,ctx->AcM,ctx->AcK,ctx->AcI,ctx->AcJ,ctx->AcI + ctx->PMM + 1,ctx->AcJ + ctx->PMM + 1,first,&nt)) {

            fprintf(ctx->PMFd,"%8d %2d  ",++n,nt);
            for (i = 0; i < ctx->PMM; ++i)
                fprintf(ctx->PMFd,"%2d ",ctx->AcN[i]);
            fprintf(ctx->PMFd,"\n");
            first = 0;
        }
    }
    else if (ctx->PMOPT == 2) {

        if (alloc_acr(ctx, ctx->PMM + 1))   
            goto CROFin;
        if (alloc_acs(ctx, ctx->PMM + 1))   
            goto CROFin;

        for (i = 0; i < ctx->PMM; ++i) {
            ctx->AcR[i] = i + 1;
            ctx->AcS[i] = ctx->PMM - i;
        }

        tda_out("R: ");
        for (i = 0; i < ctx->PMM; ++i)
            tda_out("%d ",ctx->AcR[i]);
        tda_out("\nS: ");
        for (i = 0; i < ctx->PMM; ++i)
            tda_out("%d ",ctx->AcS[i]);

        newline(ctx);
      
        first = 1;
        while (cro_gen(ctx, ctx->PMM,ctx->AcN,ctx->AcM,ctx->AcK,ctx->AcI,ctx->AcJ,ctx->AcI + ctx->PMM + 1,ctx->AcJ + ctx->PMM + 1,first,&nt)) {
            first = 0;

            if (ro_b(ctx, ctx->PMM,ctx->AcR,ctx->AcS,ctx->AcN)) {
                tda_out("N: ");
                for (i = 0; i < ctx->PMM; ++i)
                    tda_out("%d ",ctx->AcN[i]);
                newline(ctx);
            }
        }

    }
    else if (ctx->PMOPT == 3) {
        nn = 0;
        first = 1;
        while (cro_gen(ctx, ctx->PMM,ctx->AcN,ctx->AcM,ctx->AcK,ctx->AcI,ctx->AcJ,ctx->AcI + ctx->PMM + 1,ctx->AcJ + ctx->PMM + 1,first,&nt)) {
            nn++;
            first = 0;
        }
        tda_out("nn=%d\n",nn);
        if (alloc_acr(ctx, nn * ctx->PMM + 1))       
            goto CROFin;

        j = 0; 
        first = 1;
        while (cro_gen(ctx, ctx->PMM,ctx->AcN,ctx->AcM,ctx->AcK,ctx->AcI,ctx->AcJ,ctx->AcI + ctx->PMM + 1,ctx->AcJ + ctx->PMM + 1,first,&nt)) {
            first = 0;
            if (j < nn) {
                for (i = 0; i < ctx->PMM; ++i)
                    ctx->AcR[j * ctx->PMM + i] = ctx->AcN[i];
            }
            j++;
        }
        for (j = 0; j < nn; ++j) {
            tda_out("j = %5d ",j + 1);
            for (i = 0; i < ctx->PMM; ++i)
                tda_out("%3d ",ctx->AcR[j * ctx->PMM + i]);
            newline(ctx);
        }
         
        for (j = 0; j < nn; ++j) {
            for (k = j + 1; k < nn; ++k) {
                d = rod_di(ctx, ctx->PMM,ctx->AcR + j * ctx->PMM,ctx->AcR + k * ctx->PMM);
                if (d == 1) {
                    fprintf(ctx->PMFd,"%4d %4d %d\n",j+1,k+1,d);
                    n++;
                }
            }
        }
    }
    printf1(ctx, "%d records written to: %s\n",n,ctx->PMFdName);
    err = 0;

CROFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  cro_gen(m,p,q,s,c,d,c1,d1,first,nt)                                     */
/*                                                                          */
/*  create all rank orders of m alternatives, m >= 2.                       */  
/*  First call with first = 1, afterwards first = 0. The function returns   */
/*  1 if a new pattern found, otherwise 0. Number of ties groups returned   */
/*  in nt.                                                                  */

int cro_gen(TDAContext *ctx, int m,int *p,int *q,int *s,int *c,int *d,int *c1,int *d1, int first,int *nt)
{
    register int i;
    int r;

    if (first) {
        ctx->s_cro_gen_k = 2;           
        ctx->s_cro_gen_f = ctx->s_cro_gen_last = 1;
        ctx->s_cro_gen_nf = 0;

        *nt = 1;
        for (i = 0; i < m; ++i)
            p[i] = 1;
        return(1);
    }
    if (ctx->s_cro_gen_k == m)            
        goto CGP;

CGNXT:
    if (ctx->s_cro_gen_nf)
        goto CGNXT1;

    r = partit2(ctx, m,ctx->s_cro_gen_k,q,c,d,ctx->s_cro_gen_last);
    ctx->s_cro_gen_last = 0;

    if (r) {
        ctx->s_cro_gen_k++;
        if (ctx->s_cro_gen_k < m) {
            ctx->s_cro_gen_last = 1;
            goto CGNXT;
        }
        goto CGP;
    }

    /* all permutations for tie groups */

CGNXT1:
    if (ctx->s_cro_gen_nf == 0) {
        ctx->s_cro_gen_f1 = 1;
        ctx->s_cro_gen_nf = 1;
    }
    r = perm(ctx, ctx->s_cro_gen_k,s,c1,d1,ctx->s_cro_gen_f1);         
    ctx->s_cro_gen_f1 = 0;
    if (r) {  
        for (i = 0; i < m; ++i)  
            p[i] = s[q[i + 1] - 1] + 1;
        *nt = ctx->s_cro_gen_k;
        return(1);
    }
    ctx->s_cro_gen_nf = 0;
    goto CGNXT;

CGP:    /* all permutations */

    r = perm(ctx, m,q,c,d,ctx->s_cro_gen_f);         
    ctx->s_cro_gen_f = 0;
    if (r) {  
        for (i = 0; i < m; ++i)
            p[i] = q[i] + 1;
        *nt = m;
        return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  ro_b(m,a,b,c)   a[], b[], c[] are three rank orders of length m.        */
/*                                                                          */
/*  Return 1 if c[] is between a[] and b[], otherwise 0.                    */

int ro_b(TDAContext *ctx, int m,int *a,int *b,int *c)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j;

    for (i = 1; i < m; ++i) {
        for (j = 0; j < i; ++j) {
            if (a[i] < a[j]) {
                if (b[i] < b[j]) {
                    if (c[i] >= c[j])  
                        return(0);
                }
                else if (b[i] == b[j]) {
                    if (c[i] > c[j])  
                        return(0);
                }
            }
            else if (a[i] > a[j]) {
                if (b[i] == b[j]) {
                    if (c[i] < c[j])  
                        return(0);
                }
                else if (b[i] > b[j]) {
                    if (c[i] <= c[j])  
                        return(0);
                }
            }
            else {
                if (b[i] < b[j]) {
                    if (c[i] > c[j])  
                        return(0);
                }
                else if (b[i] == b[j]) {
                    if (c[i] != c[j])  
                        return(0);
                }
                else  {
                    if (c[i] < c[j])  
                        return(0);
                }
            }
        }
    }
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  ro_bd(m,a,b,c)   a[], b[], c[] are three rank orders of length m.       */
/*                                                                          */
/*  Return 1 if c[] is between a[] and b[], otherwise 0.                    */

int ro_bd(TDAContext *ctx, int m,double *a,double *b,double *c)
{
    register int i,j;

    for (i = 1; i < m; ++i) {
        for (j = 0; j < i; ++j) {

            if (fabs(a[i] - a[j]) <= ctx->EPSI1) {
                if (b[i] < b[j] - ctx->EPSI1) {
                    if (c[i] > c[j] + ctx->EPSI1)     
                        return(0);
                }
                else if (fabs(b[i] - b[j]) <= ctx->EPSI1) {
                    if (fabs(c[i] - c[j]) > ctx->EPSI1)            
                        return(0);
                }
                else  {
                    if (c[i] < c[j] - ctx->EPSI1)     
                        return(0);
                }
            }
            else if (a[i] < a[j] - ctx->EPSI1) {
                if (b[i] < b[j] - ctx->EPSI1) {
                    if (c[i] >= c[j] - ctx->EPSI1)  
                        return(0);
                }
                else if (fabs(b[i] - b[j]) <= ctx->EPSI1) {
                    if (c[i] > c[j] + ctx->EPSI1)  
                        return(0);
                }
            }
            else if (a[i] > a[j] + ctx->EPSI1) {
                if (fabs(b[i] - b[j]) <= ctx->EPSI1) {
                    if (c[i] < c[j] - ctx->EPSI1)  
                        return(0);
                }
                else if (b[i] > b[j] + ctx->EPSI1) {
                    if (c[i] <= c[j] + ctx->EPSI1)        
                        return(0);
                }
            }
        }
    }
    return(1);
}

