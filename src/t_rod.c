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

/*  functions in t_rod.c */

int rod(void);
int rod_check(int n,int m,double *x,int *p);
int rod_icheck(int m,int *x,int *y);
int rod_dcheck(int m,double *x,double *y);
int rod_s(int m,double *x,int *p);
int rod_p(int opt,int n,int m,double *x,int *p,int *pt,int *pi,int *pj,
    int *pk,int *r,int rmax,int *dmin);
int rod_d(int m,double *x,int *p);
int rod_dd(int opt,int m,double *x,double *p);
int rod_di(int m,int *x,int *p);
int rod_median(int n,int m,double *x,int *b,int *p,int *pi,int *pj,  
    int *r,int max,int *dm);
void rod_ac(int m,int *r,int *a);  
int rod_acm(int m,int *ri,int *rj,int *a);  
int rod_acc(int m,int *a);  
int rod_ar(int m,int *a,int *r,int *p,int *rp);  
int rod_between(int m,int *a,int *b,int *r,int i,int j,int *p,int *pi,
    int *pj,int max,int *rn);

int cro(void);
int cro_gen(int m,int *p,int *q,int *s,int *c,int *d,int *c1,int *d1,
    int first,int *nt);                      
int ro_b(int m,int *a,int *b,int *c);
int ro_bd(int m,double *a,double *b,double *c);


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

int rod(void)
{
    register int i,j,k,l;
    int err,n,m,r,rn,d,di,dm,opt,nt,first,ofs,mc,nw,max;
    double wt,*xi;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Rank order data. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 3,4,1))       /* get parameters */
        goto RODFin;

    if (PMOPT > 8)
        PMOPT = 8;
    if (PMALG < 1 || PMALG > 3)
        PMALG = 1;
    if (PMMax < 1)
        PMMax = 1000;
    if (PMS != 0)
        PMS = 1;

    n = NOC;
    m = PMNV;                       /* number of variables */

    if (m < 2) {
        printf1("Error: need at least two variables.\n");
        goto RODFin;
    }
    if (alloc_acx(n * m + 1))           /* input data */
        goto RODFin;

    if (PMWVar >= 0) {
        if (alloc_acw(n + 1))           /* weights */
            goto RODFin;
    }

    for (i = 0; i < n; ++i) {           /* get data into AcX and AcW */
        for (j = 0; j < m; ++j)  
            AcX[i * m + j] = get_data(PMVIdx[j],i);

        if (PMWVar >= 0) {
            wt = get_data(PMWVar,i);
            if (wt < 0.0) {
                printf1("Error: found negative weight in case %d\n",i + 1);
                goto RODFin;
            }
            AcW[i] = wt;
        }
    }
                     
       
    for (i = 0; i < n; ++i) {   
        for (j = 0; j < m; ++j)  
            printf("%g ",AcX[i * m + j]);

        if (PMWVar >= 0)  
            printf("%g ",AcW[i]);
        newline();
    }
    newline();

    if (PMOPT == 1) {             
        printf1("Creating standard format rank orders.\n");
        if (PMF1Def == 0) {
            printf1("Need output file (df option).\n");
            goto RODFin;
        }
        if (alloc_ack(m + 1))                                           
            goto RODFin;

        rn = 1;
        for (i = 0; i < n; ++i) {
            r = rod_s(m,AcX + i * m,AcK);
            rn = imax(r,rn);
        }
        for (i = 0; i < n; ++i) {
            for (j = 0; j < m; ++j) 
                fprintf(PMF1d,PMFmtS,AcX[i * m + j]);
            fprintf(PMF1d,"\n");
        }
        printf1("%d records written to: %s\n",n,PMF1dName);
        printf1("Max size of tie blocks: %d\n",rn);
        goto RODFin1;
    }
    else if (PMOPT == 2) {      /* create distance matrix */

        printf1("Distance matrix (s = %d).\n",PMS);

        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; ++j) {
                r = rod_dd(PMS,m,AcX + i * m,AcX + j * m);

                if (PMF1Def)  
                    fprintf(PMF1d,"%3d ",r);
                else
                    printf1("%3d ",r);
            }
            if (PMF1Def)  
                fprintf(PMF1d,"\n");
            else
                newline();
        }
        if (PMF1Def)  
            printf1("%d records written to: %s\n",n,PMF1dName);

        goto RODFin1;
    }
    else  if (PMOPT == 3) {     /* find mean or median */

        if (PMALG == 2)
            PMS = 0;

        printf1("Searching for ");
        if (PMS == 0)
            printf1("median");
        else
            printf1("mean");
        printf1(" rank orders (alg=%d).\n",PMALG);

        if (PMWVar >= 0)  
            printf1("Using weights defined by: %s\n",VName[PMWVar]);

        if (PMALG == 1) {
            mc = ofs = 0;

            if (alloc_acn(m + 1))   
                goto RODFin;
            if (alloc_acm(m + 1))   
                goto RODFin;
            if (alloc_ack(m + 1))   
                goto RODFin;
            if (alloc_aci(2 * m + 3))   
                goto RODFin;
            if (alloc_acj(2 * m + 3))   
                goto RODFin;
            if (alloc_acr(PMMax * m + 1))       
                goto RODFin;

            rn = 0;
            dm = 0;
            first = 1;
            while (cro_gen(m,AcN,AcM,AcK,AcI,AcJ,AcI + m + 1,AcJ + m + 1,first,&nt)) {
                mc++;
                di = 0;
                for (i = 0; i < n; ++i) {
                    d = rod_d(m,AcX + i * m,AcN);
                    if (PMS)
                        d *= d;
                    if (PMWVar >= 0)  
                        d *= (int)AcW[i];
                    di += d;
                }
                if (first) {
                    dm = di;
                    for (i = 0; i < m; ++i)
                        AcR[i] = AcN[i];
                    rn = 1;
                }
                else if (di < dm) {
                    dm = di;
                    for (i = 0; i < m; ++i)
                        AcR[i] = AcN[i];
                    rn = 1;
                }
                else if (di == dm) {
                    if (rn >= PMMax) {
                        rn = -1;
                        break;
                    }                
                    for (i = 0; i < m; ++i)
                        AcR[rn * m + i] = AcN[i];
                    rn++;
                }
                first = 0;
            }
        }
        else if (PMALG == 2) {

            mc = ofs = 0;

            if (alloc_acn(m * m + 1))       
                goto RODFin;
            if (alloc_aci(m + 1))       
                goto RODFin;
            if (alloc_acj(m + 1))       
                goto RODFin;
            if (alloc_ack(m + 1))       
                goto RODFin;
            if (alloc_acr(PMMax * m + 1))       
                goto RODFin;

            rn = rod_median(n,m,AcX,AcN,AcI,AcJ,AcK,AcR,PMMax,&dm);
        }
        else if (PMALG == 3) {
            mc = 0;
            ofs = 1;

            if (alloc_ack(m + 1))       
                goto RODFin;
            if (alloc_aci(m + 1))       
                goto RODFin;
            if (alloc_acj(m + 1))       
                goto RODFin;
            if (alloc_acn(m + 1))       
                goto RODFin;
            if (alloc_acm(PMMax * m + 1))       
                goto RODFin;
            if (alloc_acr(PMMax * m + 1))       
                goto RODFin;

            rn = rod_p(PMS,n,m,AcX,AcM,AcK,AcI,AcJ,AcN,AcR,PMMax,&dm);

        }
        if (rn < 0) {
            printf1("Error: exceeded max number of rank orders (max=%d).\n",PMMax);
            goto RODFin;
        }
        printf1("Minimal sum of distances (s = %d, mc=%d): %d\n",PMS,mc,dm);
        printf1("Number of rank orders: %d\n\n",rn);

        for (i = 0; i < rn; ++i) {
            printf1("%3d : ",i + 1);
            for (j = 0; j < m; ++j)
                printf1("%2d ",AcR[i * m + j] + ofs);
            newline();
        }
        goto RODFin1;
    }
    else  if (PMOPT >= 4 && PMOPT <= 6) {     /* create graph */

        if (alloc_acd(n * n + 1))   
            goto RODFin;

        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; ++j) {
                if (i == j)
                    continue;

                rn = 1;
                for (k = 0; k < n; ++k) {
                    if (k == i || k == j)
                        continue;
                    r = ro_bd(m,AcX + i * m,AcX + j * m,AcX + k * m);
                    if (r) {
                        rn = 0;
                        if (PMOPT <= 5)
                            break;

                        if (PMOPT == 6) {
                            AcD[i * n + k] = 1;
                            AcD[k * n + i] = 1;
                            AcD[k * n + j] = 1;
                            AcD[j * n + k] = 1;
                        }
                    }
                }
                if (rn && PMOPT <= 5) {
                    AcD[i * n + j] = 1;
                    AcD[j * n + i] = 1;
                }
            }
        }
        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; ++j)
                printf1("%d ",(int)AcD[i * n + j]);
            newline();
        }

        if (PMPCFDef || PMF1Def) {
            k = 0;
            if (PMPCFDef)                 
                fprintf(PMPCFd,"graph g {\nsize=\"5,5\";\n");

            for (i = 1; i < n; ++i) {
                for (j = 0; j < i; ++j) {

                    if (((PMOPT == 4 || PMOPT == 6) && AcD[i * n + j] != 0) ||
                        (PMOPT == 5 && AcD[i * n + j] == 0)) {

                        if (PMPCFDef)   
                            fprintf(PMPCFd,"%2d -- %2d [len=2];\n",i + 1, j + 1);

                        if (PMF1Def)  
                            fprintf(PMF1d,"%2d %2d %2d\n",i + 1, j + 1,1);
                        k++;
                    }   
                }
            }
            if (PMF1Def)  
                printf1("%d records written to: %s\n",k,PMF1dName);
            if (PMPCFDef) { 
                fprintf(PMPCFd,"}\n");
                printf1("%d records written to: %s\n",k + 3,PMPCFName);
            }
        }
        goto RODFin1;
    }
    else  if (PMOPT == 7) {     /* find rank orders in between  */

        if (n < 2) {
            printf1("Error: need at least two rank orders.\n");
            goto RODFin;
        }
        max = m * m * m * m;

        if (alloc_acr(m + 1))   
            goto RODFin;
        if (alloc_acs(m + 1))   
            goto RODFin;
        if (alloc_acn(m * m + 1))   
            goto RODFin;
        if (alloc_acm(m * m + 1))   
            goto RODFin;
        if (alloc_ack(m * m + 1))   
            goto RODFin;
        if (alloc_aci(2 * m + 4))   
            goto RODFin;
        if (alloc_acj(PMMax * m + 1))   
            goto RODFin;
        if (alloc_acl(max * m + 1))   
            goto RODFin;

        rn = 0;
        for (i = 1; i < n; ++i) {
                 
            for (k = 0; k < m; ++k)  
                AcR[k] = (int)AcX[i * m + k];

            for (j = 0; j < i; ++j) {
                 
                for (k = 0; k < m; ++k)  
                    AcS[k] = (int)AcX[j * m + k];

                rod_ac(m,AcR,AcN);
                rod_ac(m,AcS,AcM);

                d = 0;
                if (rod_between(m,AcN,AcM,AcK,0,0,AcL,AcI,AcI + m + 1,max,&d)) {
                    printf1("Error: exceeded max number of rank orders.\n");
                    goto RODFin;
                }
                for (k = 0; k < d; ++k) {
                    nt = 1;
                    for (l = 0; l < rn; ++l) {
                        if (rod_icheck(m,AcL + k * m,AcJ + l * m) == 0) {
                            nt = 0;
                            break;
                        }
                    }
                    if (nt) {
                        if (rn >= PMMax) {
                            printf1("Error: exceeded max number of rank orders (max=%d).\n",PMMax);
                            goto RODFin;
                        }
                        for (l = 0; l < m; ++l)
                            AcJ[rn * m + l] = AcL[k * m + l];
                        rn++;
                    }
                }
            }
        }
        printf1("Found %d rank orders.\n",rn);
        if (PMF1Def) {
            for (i = 0; i < rn; ++i) {
                for (j = 0; j < m; ++j) 
                    fprintf(PMF1d,"%2d ",AcJ[i * m + j]);
                fprintf(PMF1d,"\n");
            }
            printf1("%d records written to: %s\n",rn,PMF1dName);
        }
        goto RODFin1;
    }
    else  if (PMOPT == 8) {     /* standard integer format */

        printf1("Tabulate rank orders.\n");

        if (alloc_acy(n * m + 1))   
            goto RODFin;
        if (alloc_acn(n + 1))   
            goto RODFin;

        for (j = 0; j < m; ++j)
            AcY[j] = AcX[j];
        if (PMWVar >= 0)  
            AcN[0] = (int)AcW[0];
        else 
            AcN[0] = 1;

        rn = 1;
        for (i = 1; i < n; ++i) {
            r = -1;
            xi = AcX + i * m;
            for (j = 0; j < rn; ++j) {
                if (rod_dcheck(m,xi,AcY + j * m) == 0) {
                    r = j;
                    break;
                }   
            }   
            if (r >= 0) {
                if (PMWVar >= 0)  
                    AcN[r] += (int)AcW[i];
                else 
                    AcN[r] += 1;
            }
            else {
                for (j = 0; j < m; ++j)
                    AcY[rn * m + j] = xi[j];
                if (PMWVar >= 0)  
                    AcN[rn] = (int)AcW[i];
                else 
                    AcN[rn] = 1;
                rn++;
            }
        }
        nw = 0;
        for (i = 0; i < rn; ++i)
            nw += AcN[i];

        printf1("Found %d different rank orders.\n",rn);
        printf1("Sum of weights: %d\n",nw);

        if (PMF1Def) {
            for (i = 0; i < rn; ++i) {
                xi = AcY + i * m;
                for (j = 0; j < m; ++j)
                    fprintf(PMF1d,"%g ",xi[j]);
                fprintf(PMF1d,"%5d\n",AcN[i]);
            }
            printf1("%d records written to: %s\n",rn,PMF1dName);
        }
        goto RODFin1;
    }

RODFin1:
    newline();
    err = 0;

RODFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  rod__check()    return 1 if p in x.                                     */
       
int rod_check(int n,int m,double *x,int *p)
{
    register int i,j,fnd;
    double *xi;

    for (i = 0; i < n; ++i) {
        xi = x + i * m;
        fnd = 1;
        for (j = 0; j < m; ++j) {
            if (fabs(xi[j] - (double)p[j]) > EPSI1) {
                fnd = 0;
                break;
            }
        }
        if (fnd)
            return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  rod__icheck()   return 1 if x and y are different rank orders.          */
       
int rod_icheck(int m,int *x,int *y)
{
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
       
int rod_dcheck(int m,double *x,double *y)
{
    register int i,j;

    for (i = 1; i < m; ++i) {
        for (j = 0; j < i; ++j) {
            if (x[j] < x[i] - EPSI1) {
                if (y[j] >= y[i] - EPSI1)
                    return(1);
            }
            else if (x[j] > x[i] + EPSI1) {
                if (y[j] <= y[i] + EPSI1)
                    return(1);
            }
            else if (fabs(x[j] - x[i]) <= EPSI1) {
                if (fabs(y[j] - y[i]) > EPSI1)
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
       
int rod_s(int m,double *x,int *p)
{
    register int i,j,n;       
    int r;
    double xa;

    r = 1;
    sortdp(m,x,p);
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
        r = imax(r,n);
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
       
int rod_p(int opt,int n,int m,double *x,int *p,int *pt,int *pi,int *pj,
    int *pk,int *r,int rmax,int *dm)
{
    register int i,j,k,l,q;
    int first,d,di,npi,rn0,rn1,ii,a,b;

    /* first step: check with complete ties  */

    for (i = 0; i < m; ++i)
        pk[i] = 0;

    *dm = 0;
    for (i = 0; i < n; ++i) {
        di = rod_d(m,x + i * m,pk);
        if (opt)
            di *= di;
        *dm += di;
    }

/*  printf("FIRST dm=%d  \n",*dm);      */

    for (i = 0; i < m; ++i)     /* save in p[] */
        p[i] = pt[i];
    rn0 = 1;

    /* second step: check all permutations without ties, save optimal
       permutations in p[] (i = 0,...,rn0-1) */

    first = 1;
    while (perm(m,pt,pi,pj,first)) {

        d = 0;
        for (i = 0; i < n; ++i) {
            di = rod_d(m,x + i * m,pt);
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
      
printf("\nAfter second step. rn0=%d dm=%d \n",rn0,*dm);
for (i = 0; i < rn0; ++i) {
    for (j = 0; j < m; ++j)
        printf("%d ",p[i * m + j]);
    newline();
}

printf("dm=%d\n",*dm);
       


    /*  third step: check all possible tie groups */

    for (i = 0; i < rn0 * m; ++i)   /* copy to r[] */
        r[i] = p[i];
    rn1 = rn0;

    for (ii = 0; ii < rn0; ++ii) {  /* for all permutations from second step */

        for (k = 2; k < m; ++k) {
printf("\nii=%d k=%d\n",ii,k);
            first = 1;
            while (partit1(m,k,pi,pj,first) == 0) {
      
                printf("pi: ");
                for (i = 1; i <= k; ++i)
                    printf("%d ",pi[i]);
                newline();
        


                l = 0;
                for (i = 1; i <= k; ++i) {
                    npi = pi[i];
                    q = l;
                    for (j = 0; j < npi; ++j)
                        pt[l++] = q;
                }

                       
                printf("pt: ");
                for (i = 0; i < m; ++i)
                    printf("%d ",pt[i]);
                newline();
                       

       
                for (i = 0; i < m; ++i)  
                    pk[i] = pt[p[ii * m + i]];
        
                         
                printf("pk: ");
                for (i = 0; i < m; ++i)
                    printf("%d ",pk[i]);
                newline();
        
                printf("xxx: ");
                for (i = 0; i < m; ++i)
                    printf("%g ",x[i]);
                newline();
                         

                d = 0;
                for (i = 0; i < n; ++i) {
                    di = rod_d(m,x + i * m,pk);
                    if (opt)
                        di *= di;
                    d += di;
                    if (d > *dm)
                        break;
                }


                printf("d========%d dm=%d  rn1=%d \n",d,*dm,rn1);   

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
       
int rod_d(int m,double *x,int *p)
{
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
       
int rod_dd(int opt,int m,double *x,double *p)
{
    register int i,j;
    int d,di;
       
    di = 0;
    for (i = 1; i < m; ++i) {
        for (j = 0; j < i; ++j) {
            d = 0;
            if (fabs(x[i] - x[j]) <= EPSI1) {
                if (fabs(p[i] - p[j]) > EPSI1)
                    d++;           
            }
            else if (x[i] < x[j] - EPSI1) {
                if (fabs(p[i] - p[j]) <= EPSI1)
                    d++;
                else if (p[i] > p[j] + EPSI1)
                    d += 2;
            }
            else {
                if (fabs(p[i] - p[j]) <= EPSI1)
                    d++;
                else if (p[i] < p[j] - EPSI1)
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
       
int rod_di(int m,int *x,int *p)
{
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

int rod_median(int n,int m,double *x,int *b,int *p,int *pi,int *pj,
    int *r,int max,int *dm)
{
    register int i,j,k,l,ii;
    int bjk,first,d,rn,rn0,pmax,*ri,*rii;
    double *xi;

    /* create b matrix */

    for (j = 0; j < m; ++j) {
        for (k = 0; k < m; ++k) {
            bjk = 0;
            if (j != k) {
                for (i = 0; i < n; ++i) {
                    xi = x + i * m;
                    if (xi[j] < xi[k] - EPSI1)
                        bjk++;
                    else if (xi[j] > xi[k] + EPSI1)
                        bjk--;
                }
            }
            b[j * m + k] = bjk;
        }
    }
    rn = 0;
    *dm = 0;
    first = 1;
    while (perm(m,p,pi,pj,first)) {
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
      
            if (rod_acm(m,ri,rii,b))  
                continue;
            if (rod_acc(m,b))  
                continue;
            if (rod_ar(m,b,p,pi,pj))             
                continue;

            /* check whether already there */

            first = 1;
            for (k = rn0; k < rn; ++k) {
                if (rod_icheck(m,p,r + k * m) == 0) {
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

void rod_ac(int m,int *r,int *a)   
{
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

int rod_acm(int m,int *ri,int *rj,int *a)   
{
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

int rod_acc(int m,int *a)   
{
    register int i,j,k;
    int aij,aik,akj;

    for (i = 0; i < m; ++i) {
        for (j = 0; j < m; ++j) {
            if (i == j) 
                continue;

            aij = a[i * m + j];
            if (aij != -a[j * m + i]) {
                printf("NOT0\n");
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

int rod_ar(int m,int *a,int *r,int *p,int *rp)
{
    register int i,j,k,l;
    int pmax,aij,aik,akj;
             
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
                    fprintf(stderr,"FATAL ERROR 1 in ROD_AR.\n");
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
                    fprintf(stderr,"FATAL ERROR 2 in ROD_AR.\n");
                    exit(0);
                }
            }
            else if (k == 0) {
                if (rp[i] != rp[j])  
                    rp[j] = rp[i];
            }
            else {
                fprintf(stderr,"FATAL ERROR 3 in ROD_AR.\n");
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

int rod_between(int m,int *a,int *b,int *r,int i,int j,int *p,int *pi,
    int *pj,int max,int *rn)
{
    register int k,l;
    int aij,bij,rij;

    if (++j >= m) {
        if (++i == m - 1) {
            for (i = 1; i < m; ++i) {
                for (j = 0; j < i; ++j)
                    r[i * m + j] = -r[j * m + i];
            }

            /* check for new rank order */

            if (rod_acc(m,r) == 0) {        /* if rank order matrix */
                if (*rn >= max)
                    return(-1);

                rod_ar(m,r,p + *rn * m,pi,pj);
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
            if (rod_between(m,a,b,r,i,j,p,pi,pj,max,rn))
                return(-1);
        }
    }
    else  {
        for (rij = aij; rij >= bij; --rij) {
            r[i * m + j] = rij;
            if (rod_between(m,a,b,r,i,j,p,pi,pj,max,rn))
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

int cro(void)
{
    register int i,j,k; 
    int err,n,nn,nt,first,d;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Creating rank orders. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto CROFin;

    if (PMOPT > 3)
        PMOPT = 3;
    if (PMM < 2)
        PMM = 2;

    if (alloc_acn(PMM + 1))   
        goto CROFin;
    if (alloc_acm(PMM + 1))   
        goto CROFin;
    if (alloc_ack(PMM + 1))   
        goto CROFin;
    if (alloc_aci(2 * PMM + 3))   
        goto CROFin;
    if (alloc_acj(2 * PMM + 3))   
        goto CROFin;

    n = 0;
    if (PMOPT == 1) {           /* create rank orders */
        first = 1;
        while (cro_gen(PMM,AcN,AcM,AcK,AcI,AcJ,AcI + PMM + 1,AcJ + PMM + 1,first,&nt)) {

            fprintf(PMFd,"%8d %2d  ",++n,nt);
            for (i = 0; i < PMM; ++i)
                fprintf(PMFd,"%2d ",AcN[i]);
            fprintf(PMFd,"\n");
            first = 0;
        }
    }
    else if (PMOPT == 2) {

        if (alloc_acr(PMM + 1))   
            goto CROFin;
        if (alloc_acs(PMM + 1))   
            goto CROFin;

        for (i = 0; i < PMM; ++i) {
            AcR[i] = i + 1;
            AcS[i] = PMM - i;
        }

        printf("R: ");
        for (i = 0; i < PMM; ++i)
            printf("%d ",AcR[i]);
        printf("\nS: ");
        for (i = 0; i < PMM; ++i)
            printf("%d ",AcS[i]);

        newline();
      
        first = 1;
        while (cro_gen(PMM,AcN,AcM,AcK,AcI,AcJ,AcI + PMM + 1,AcJ + PMM + 1,first,&nt)) {
            first = 0;

            if (ro_b(PMM,AcR,AcS,AcN)) {
                printf("N: ");
                for (i = 0; i < PMM; ++i)
                    printf("%d ",AcN[i]);
                newline();
            }
        }

    }
    else if (PMOPT == 3) {
        nn = 0;
        first = 1;
        while (cro_gen(PMM,AcN,AcM,AcK,AcI,AcJ,AcI + PMM + 1,AcJ + PMM + 1,first,&nt)) {
            nn++;
            first = 0;
        }
        printf("nn=%d\n",nn);
        if (alloc_acr(nn * PMM + 1))       
            goto CROFin;

        j = 0; 
        first = 1;
        while (cro_gen(PMM,AcN,AcM,AcK,AcI,AcJ,AcI + PMM + 1,AcJ + PMM + 1,first,&nt)) {
            first = 0;
            if (j < nn) {
                for (i = 0; i < PMM; ++i)
                    AcR[j * PMM + i] = AcN[i];
            }
            j++;
        }
        for (j = 0; j < nn; ++j) {
            printf("j = %5d ",j + 1);
            for (i = 0; i < PMM; ++i)
                printf("%3d ",AcR[j * PMM + i]);
            newline();
        }
         
        for (j = 0; j < nn; ++j) {
            for (k = j + 1; k < nn; ++k) {
                d = rod_di(PMM,AcR + j * PMM,AcR + k * PMM);
                if (d == 1) {
                    fprintf(PMFd,"%4d %4d %d\n",j+1,k+1,d);
                    n++;
                }
            }
        }
    }
    printf1("%d records written to: %s\n",n,PMFdName);
    err = 0;

CROFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  cro_gen(m,p,q,s,c,d,c1,d1,first,nt)                                     */
/*                                                                          */
/*  create all rank orders of m alternatives, m >= 2.                       */  
/*  First call with first = 1, afterwards first = 0. The function returns   */
/*  1 if a new pattern found, otherwise 0. Number of ties groups returned   */
/*  in nt.                                                                  */

int cro_gen(int m,int *p,int *q,int *s,int *c,int *d,int *c1,int *d1,
    int first,int *nt)               
{
    register int i;
    int r;
    static int k,f,f1,last,nf;

    if (first) {
        k = 2;           
        f = last = 1;
        nf = 0;

        *nt = 1;
        for (i = 0; i < m; ++i)
            p[i] = 1;
        return(1);
    }
    if (k == m)            
        goto CGP;

CGNXT:
    if (nf)
        goto CGNXT1;

    r = partit2(m,k,q,c,d,last);
    last = 0;

    if (r) {
        k++;
        if (k < m) {
            last = 1;
            goto CGNXT;
        }
        goto CGP;
    }

    /* all permutations for tie groups */

CGNXT1:
    if (nf == 0) {
        f1 = 1;
        nf = 1;
    }
    r = perm(k,s,c1,d1,f1);         
    f1 = 0;
    if (r) {  
        for (i = 0; i < m; ++i)  
            p[i] = s[q[i + 1] - 1] + 1;
        *nt = k;
        return(1);
    }
    nf = 0;
    goto CGNXT;

CGP:    /* all permutations */

    r = perm(m,q,c,d,f);         
    f = 0;
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

int ro_b(int m,int *a,int *b,int *c)
{
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

int ro_bd(int m,double *a,double *b,double *c)
{
    register int i,j;

    for (i = 1; i < m; ++i) {
        for (j = 0; j < i; ++j) {

            if (fabs(a[i] - a[j]) <= EPSI1) {
                if (b[i] < b[j] - EPSI1) {
                    if (c[i] > c[j] + EPSI1)     
                        return(0);
                }
                else if (fabs(b[i] - b[j]) <= EPSI1) {
                    if (fabs(c[i] - c[j]) > EPSI1)            
                        return(0);
                }
                else  {
                    if (c[i] < c[j] - EPSI1)     
                        return(0);
                }
            }
            else if (a[i] < a[j] - EPSI1) {
                if (b[i] < b[j] - EPSI1) {
                    if (c[i] >= c[j] - EPSI1)  
                        return(0);
                }
                else if (fabs(b[i] - b[j]) <= EPSI1) {
                    if (c[i] > c[j] + EPSI1)  
                        return(0);
                }
            }
            else if (a[i] > a[j] + EPSI1) {
                if (fabs(b[i] - b[j]) <= EPSI1) {
                    if (c[i] < c[j] - EPSI1)  
                        return(0);
                }
                else if (b[i] > b[j] + EPSI1) {
                    if (c[i] <= c[j] + EPSI1)        
                        return(0);
                }
            }
        }
    }
    return(1);
}

