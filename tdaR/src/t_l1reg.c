/****************************************************************************/
/*  t_l1reg                                                                 */
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

#include "tda.h"
#include "t_gen.h"
#include "t_pgen.h"
#include "t_parm.h"
#include "t_alloc.h"
#include "t_var.h"
#include "t_gdat.h"
#include "t_gf.h"
#include "t_sort.h"
#include "t_lsei.h"
#include "t_mat.h"
#include "t_matf.h"
#include "tda_context.h"

/*  functions in t_l1reg */

int l1reg(TDAContext *ctx);
int l1regf(TDAContext *ctx, int m, int n, double *x, double *y, double *s, double *res, int *rank);
void l1reg_res(TDAContext *ctx, int n,int nx,double *x,double *y,double *res);

/* ------------------------------------------------------------------------ */
/*  l1reg()         L1 norm regression.                                     */
/*                                                                          */
/*                  l1reg(                                                  */
/*                      ni=...,     1 if without intercept, def. 0          */
/*                      sel=...,    optional case selection                 */
/*                      ppar=...,   print parameter to output file          */
/*                      tfmt=...,   print format for parameters, def. 10.4  */
/*                      pres=...,   print residuals                         */
/*                      pcov=...,   unscaled covariance matrix              */
/*                      mfmt=...,   printformat for pcov option, def. 12.4  */
/*                      df=...,     optional output file                    */
/*                      fmt=...,    print format for df option, def. 0.0    */  
/*                  ) = Y,X1,X2,...;                                        */
/*                                                                          */
/*                  Note: only with cross-sectional data.                   */
/*                                                                          */
/*                  Return 0 if OK, otherwise -1.                           */

int l1reg(TDAContext *ctx)
{
    register int i,j;
    int err,nv,nx,nx1,m,n,k1,k2;
    int rank,ierr,nn,df;
    double res,tmp,tmp1,e1,e2,tau,rnorm;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "L1-norm regression. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,4,1))     /* get parameters */
        goto L1Fin;

    if (ctx->PMNW != 1) {
        printf1(ctx, "Error: l1reg command requires cross-section data (nw=1).\n");
        goto L1Fin;
    }
    newline(ctx);
    nv = check_nvar(ctx, 0);         /* check variables */
    if (nv == 0)            
        goto L1Fin;        
       
    prn_nwvar(ctx, 0);               /* print variables */
    nx = nx1 = nv - 1;
    if (ctx->PMNI)  
        printf1(ctx, "Model without intercept.\n");
    else {
        ctx->PMNI = 0;
        nx1++;
        if (nx == 0)   
            printf1(ctx, "Model without independent variables.\n");
    }
    if (nx1 == 0)
        goto L1Fin;

    /* allocate data matrix AcX[], AcY[], and solution vector AcU[] */

    if (alloc_acx(ctx, nx1 * ctx->NOC + 1))
        goto L1Fin;
    if (alloc_acy(ctx, ctx->NOC + 1))
        goto L1Fin;
    if (alloc_acu(ctx, nx1 + 1))
        goto L1Fin;
    if (alloc_acv(ctx, nx1 + 1))
        goto L1Fin;

    nn = getdxy(ctx, ctx->AcX,ctx->AcY,nx1,ctx->PMNI,1);    /* get data */

    if (nn == 0)
        goto L1Fin;
  
    if (ctx->PMF1Def) {          
        prn1_data(ctx, ctx->AcY,ctx->AcX,nn,1,nx1,ctx->PMF1d,ctx->PMFmtS);
        p_wmsg(ctx, 1,ctx->PMF1dName,1);
    }
    if (nn < nx1) {
        p_err(ctx, -30,1);
        goto L1Fin;
    }
    newline(ctx);

    ierr = l1regf(ctx, nn,nx1,ctx->AcX,ctx->AcY,ctx->AcU,&res,&rank);  /* perform regression */

    if (ierr == -2) {   /* insuff memory */
        p_err(ctx, -2,1);
        goto L1Fin;
    }
    if (rank != nx1) {
        printf1(ctx, "Rank of data matrix  %d\n",rank); 
        printf1(ctx, "Less than full rank. No regression performed.\n");
        goto L1Fin;
    }
    switch (ierr) {
        case 0:  printf1(ctx, "Solution is unique\n"); 
                 break;
        case 1:  printf1(ctx, "Solution is possibly not unique\n"); 
                 break;
        case 2:  printf1(ctx, "Terminated prematurely due to rounding errors\n");
                 break;
        default: break;
    }
    printf1(ctx, "Rank of data matrix: %d\n",rank); 
    printf1(ctx, "Sum of absolute deviations: %lg\n",res);

    if (ctx->PMMPLogDef == 1)                /* create matrix for sum of resid */
        mp_putlog(ctx, res);

    if (ctx->PMMPParDef == 1)                /* create matrix for parameters */
        mp_putpar(ctx, nx1,ctx->AcU);

    /* calculate residuals in AcTmp */                

    if (alloc_actmp(ctx, nn + 1))
        goto L1Fin;

    for (i = 0; i < nn; ++i) {
        tmp = 0.0;
        for (j = 1; j <= nx1; ++j)
            tmp += ctx->AcX[i * nx1 + j] * ctx->AcU[j];
        ctx->AcTmp[i + 1] = ctx->AcY[i + 1] - tmp;
    }

    if (ctx->PMResFDef)                         /* write residuals */                 
        l1reg_res(ctx, nn,nx1,ctx->AcX,ctx->AcY,ctx->AcTmp);

    /* calculate covariance */

    df = 0;
    if (nn > nx1) {

        if (sortd(ctx, nn,ctx->AcTmp + 1,0))
            goto L1Fin;

        n = 0;
        for (i = 1; i <= nn; ++i) {
            if (fabs(ctx->AcTmp[i]) <= ctx->EPSI1)
                n++;
        }
        m = nn - n;
        printf1(ctx, "Number of non-zero residuals: %d\n",m);
        if (m >= 2) {

            tmp = ((double)m + 1.0) / 2.0;
            tmp1 = sqrt((double)m);
            k1 = (int)(tmp - tmp1 + 0.5);
            k2 = (int)(tmp + tmp1 + 0.5);
            if (k1 < 1)
                k1 = 1;
            if (k2 > m)
                k2 = m;

            e1 = ctx->AcTmp[k1];
            e2 = ctx->AcTmp[k2 + n];
            tau = tmp1 * (e2 - e1) / 4.0;

            if (tau > 0.0) {
                printf1(ctx, "\nTau: ");
                rt_printf1_d(ctx, ctx->PMTFmtS,tau);
                printf1(ctx, " [calculated from residuals: %lg,%lg]\n",e1,e2);

                if (alloc_aci(ctx, nx1 + 1))
                    goto L1Fin;
                if (alloc_acw(ctx, nx1 + 1))
                    goto L1Fin;

                rank = lhhfti(ctx, nn,nx1,nx1,1,ctx->AcX,ctx->AcY,&rnorm,2,ctx->AcV,ctx->AcW,ctx->AcI,ctx->EPSI1);

                printf1(ctx, "Pseudorank of data matrix: %d\n",rank); 
                printf1(ctx, "Norm of least squares residuals: "); 
                rt_printf1_d(ctx, ctx->PMTFmtS,rnorm);
                printf1(ctx, "\nDegrees of freedom: %d\n",nn - nx1);

                if (ctx->PMMPCovDef == 1)              /* create matrix for covariance */
                    mp_putcov(ctx, nx1,ctx->AcX);
#ifdef TDA_R_PACKAGE
                export_prn(ctx, "l1reg.vcov", nx1,nx1,nx1,ctx->AcX);
#endif

                if (ctx->PMCovFDef) {                          /* write cov matrix */
                    prn_data(ctx, nx1,nx1,nx1,ctx->AcX,ctx->PMCovFd,ctx->PMMFmtS);
                    printf1(ctx, "Unscaled covariance matrix written to: %s\n",ctx->PMCovFName);
                }
                if (rank == nx1) {
                    for (j = 1; j <= nx1; ++j)   
                        ctx->AcV[j] = tau * tau * ctx->AcX[(j - 1) * nx1 + j];
                    df = nn - nx1;
                }
            }
        }
    }
    prn1_coeff(ctx, nx1,ctx->AcU,ctx->AcV,ctx->PMNI,df,ctx->PMVIdx,1);     /* print estimated coeff. */

    mp_info(ctx);              /* print info about matrices */

    err = 0;

L1Fin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  l1regf      L1 norm regression.                                         */
/*              L1-norm (xb - y) ---> minimum                               */
/*                                                                          */
/*      Ref: Barrodale and Roberts, Solution of an Overdetermined System    */
/*      of Equations in the L1 Norm, Algorithm 478, Comm. ACM, Vol. 17,     */
/*      1974, p. 319 - 320.                                                 */
/*                                                                          */
/*      int l1regf(m,n,x,y,s,&res,&rank)                                    */
/*                                                                          */
/*      int m       input:  number of cases                                 */
/*      int n       input:  number of variables                             */
/*      double *x   input:  (m,n) matrix with data for x variables          */
/*                  output: not changed                                     */
/*      double *y   input:  (m) vector with y variables                     */
/*                  output: not changed                                     */
/*      double *s   output: solution vector                                 */
/*      double res  output: sum of absolute residuals                       */
/*      int rank    output: rank of x matrix                                */
/*                                                                          */
/*      Return:     0  if a unique optiomal solution was found              */
/*                  1  if an optimal, but not neccessarily unique,          */
/*                     solution was found                                   */
/*                  2  if premature end in consequence of rounding errors.  */
/*                 -2  if insufficient memory                               */
/*                                                                          */
/*      Note: the function works currently with a working array of          */
/*      dimension (m+1,n+2); the x matrix is in the first m rows and        */
/*      first n columns.                                                    */

#define toler 1.0e-8    /* tolerance for zero */

int l1regf(TDAContext *ctx, int m, int n, double *x, double *y, double *s, double *res, int *rank)
{
    register int i,j,k,l;
    int err,in = 0,out = 0,m1,n1,stage,test = 0,kount,kr,kl,m2,n2,*si;
    double bigdbl,sum,min,max,d,pivot,*a,*b,*e;

    bigdbl = ctx->DBLMAX / 100.0;

    err = -2;
    m1 = m + 1;
    m2 = m + 2;
    n1 = n + 1;
    n2 = n + 2;

    if (!(a = (double *) calloc((size_t)(m2) * (size_t)(n2) + 1,sizeof(double))))  
        goto L1FFin;
    memrq(ctx, m2 * n2 + 1,sizeof(double));
  
    if (!(b = (double *) calloc((size_t)(m2 + 1),sizeof(double))))  
        goto L1FFin1;
    memrq(ctx, m2 + 1,sizeof(double));
  
    if (!(e = (double *) calloc((size_t)(m  + 1),sizeof(double))))  
        goto L1FFin2;
    memrq(ctx, m + 1,sizeof(double));
  
    if (!(si = (int *) calloc((size_t)(m + 1),sizeof(int))))  
        goto L1FFin3;
    memrq(ctx, m + 1,sizeof(int));
        
    /* copy because of alternative matrix-format */

    for (i = 1; i <= m; ++i) {
        b[i] = y[i];
        for (j = 1; j <= n; ++j)
            a[(i - 1) * n2 + j] = x[(i - 1) * n + j];
    }

    for (j = 1; j <= n; ++j) {
        a[(m2 - 1) * n2 + j] = j;
        s[j] = 0.0;
    }
    for (i = 1; i <= m; ++i) {
        a[(i - 1) * n2 + n2] = n + i;
        a[(i - 1) * n2 + n1] = b[i];

        if (b[i] < 0.0) {
            for (j = 1; j <= n2; ++j)
                a[(i - 1) * n2 + j] = -a[(i - 1) * n2 + j];
        }
        e[i] = 0.0;
    }
    for (j = 1; j <= n1; ++j) {
        sum = 0.0;

        for (i = 1; i <= m; ++i)
            sum += a[(i - 1) * n2 + j];

        a[(m1 - 1) * n2 + j] = sum;
    }
    stage = 1;
    kount = 0;
    kr = kl = 1;

L1Reg70:
    max = -1.0;

    for (j = kr; j <= n; ++j) {

        if (fabs(a[(m2 - 1) * n2 + j]) <= n) {
            d = fabs(a[(m1 - 1) * n2 + j]);
            if (d > max) {
                max = d;
                in = j;
            }
        }
    }
    if (a[(m1 - 1) * n2 + in] < 0.0) {
        for (i = 1; i <= m2; ++i)
            a[(i - 1) * n2 + in] = -a[(i - 1) * n2 + in];
    }
L1Reg100:
    k = 0;

    for (i = kl; i <= m; ++i) {
        d = a[(i - 1) * n2 + in];
        if (d > toler) {
            k++;
            b[k] = a[(i - 1) * n2 + n1] / d;
            si[k] = i;
            test = 1;
        }
    }
    while (1) {

        if (k <= 0) {
            test = 0; 
        }
        else {
            min = bigdbl;
            for (i = 1; i <= k; ++i) {
                if (b[i] < min) {
                    j = i;
                    min = b[i];
                    out = si[i];
                }
            }
            b[j] = b[k];
            si[j] = si[k];
            k--;
        }
        if (test != 1 && stage != 0) {
            for (i = 1; i <= m2; ++i) {
                d = a[(i - 1) * n2 + kr];
                a[(i - 1) * n2 + kr] = a[(i - 1) * n2 + in];
                a[(i - 1) * n2 + in] = d;
            }
            kr++;
            goto L1Reg260;
        }
        else {
            if (test != 1) {
                a[(m2 - 1) * n2 + n1] = 2.0;
                goto L1Reg350;
            }
            pivot = a[(out - 1) * n2 + in];

            if ((a[(m1 - 1) * n2 + in] - pivot - pivot) <= toler) 
                break;

            for (j = kr; j <= n1; ++j) {
                d = a[(out - 1) * n2 + j];
                a[(m1 - 1) * n2 + j] -= d + d;
                a[(out - 1) * n2 + j] = -d;
            }
            a[(out - 1) * n2 + n2] = -a[(out - 1) * n2 + n2];
        }
    }
    for (j = kr; j <= n1; ++j) {
        if (j != in)  
            a[(out - 1) * n2 + j] /= pivot;
    }
    for (i = 1; i <= m1; ++i) {

        if (i != out) {
            d = a[(i - 1) * n2 + in];
            for (j = kr; j <= n1; ++j) {
                if (j != in)  
                    a[(i - 1) * n2 + j] -= d * a[(out - 1) * n2 + j];
            }
        }
    }
    for (i = 1; i <= m1; ++i) {
        if (i != out)
            a[(i - 1) * n2 + in] = -a[(i - 1) * n2 + in] / pivot;
    }
    a[(out - 1) * n2 + in] = 1.0 / pivot; 
    d = a[(out - 1) * n2 + n2]; 
    a[(out - 1) * n2 + n2] = a[(m2 - 1) * n2 + in];

    a[(m2 - 1) * n2 + in] = d;
    kount++;

    if (stage != 0) {
        kl++;
        for (j = kr; j <= n2; ++j) {
            d = a[(out - 1) * n2 + j]; 
            a[(out - 1) * n2 + j] = a[(kount - 1) * n2 + j];
            a[(kount - 1) * n2 + j] = d;
        }

L1Reg260:
        if ((kount + kr) != n1)
            goto L1Reg70;

        stage = 0;
    }
    max = -bigdbl;

    for (j = kr; j <= n; ++j) {
        d = a[(m1 - 1) * n2 + j];
        if (d < 0.0) {
            if (d > -2.0)
                goto L1Reg290;
            d = -d - 2.0;
        }
        if (d > max) {
            max = d; 
            in = j;
        }
L1Reg290: ;
    }
    if (max > toler) {
        if (a[(m1 - 1) * n2 + in] <= 0.0) {
            for (i = 1; i <= m2; ++i)
                a[(i - 1) * n2 + in] = -a[(i - 1) * n2 + in];

            a[(m1 - 1) * n2 + in] -= 2.0;
        }
        goto L1Reg100;
    }
    l = kl - 1;

    for (i = 1; i <= l; ++i) {
        if (a[(i - 1) * n2 + n1] < 0.0) {
            for (j = kr; j <= n2; ++j)
                a[(i - 1) * n2 + j] = -a[(i - 1) * n2 + j];
        }
    }
    a[(m2 - 1) * n2 + n1] = 0.0;

    if (kr == 1) {
        for (j = 1; j <= n; ++j) {
            d = fabs(a[(m1 - 1) * n2 + j]);
            if (d <= toler || (2.0 - d) <= toler) 
                goto L1Reg350;
        }
        a[(m2 - 1) * n2 + n1] = 1.0;
    }
L1Reg350:

    for (i = 1; i <= m; ++i) {
        k = (int)a[(i - 1) * n2 + n2]; 
        d = a[(i - 1) * n2 + n1];

        if (k <= 0) {
            k = -k;
            d = -d;
        }
        if (i < kl) {
            s[k] = d;
        }
        else {
            k -= n;
            e[k] = d;
        }
    }
    a[(m2 - 1) * n2 + n2] = kount; 
    a[(m1 - 1) * n2 + n2] = n1 - kr;
    sum = 0.0;

    for (i = kl; i <= m; ++i)
        sum += a[(i - 1) * n2 + n1];

    *res = sum;
    *rank = (int)a[(m1 - 1) * n2 + n2];

    if (a[(m2 - 1) * n2 + n1] == 1)
        err = 0;

    else if (a[(m2 - 1) * n2 + n1] == 0)
        err = 1;
    else           
        err = 2;
            
    free((char *)si);
    memrq(ctx, -m - 1,sizeof(int));
L1FFin3:
    free((char *)e);
    memrq(ctx, -m - 1,sizeof(double));
L1FFin2:
    free((char *)b);
    memrq(ctx, -m2 - 1,sizeof(double));
L1FFin1:
    free((char *)a);
    memrq(ctx, -m2 * n2 - 1,sizeof(double));
L1FFin:
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  l1reg_res(n,nx,x,y,res)    Write residuals to PMResFd.                  */
/*                                                                          */
/*  n  = number of cases in x,y                                             */
/*  nx = number of variables, columns in x                                  */
/*  x  = X matrix (n,nx)                                                    */
/*  y  = Y vector (n vector)                                                */
/*  res = residuals                                                         */
/*                                                                          */  

void l1reg_res(TDAContext *ctx, int n,int nx,double *x,double *y,double *res)
{
    register int i,j,l;

    l = 1;
    fprintf(ctx->PMResFd,"# residuals created by l1reg command.\n");
    fprintf(ctx->PMResFd,"# c%-3d : case number\n",l++);
    if (ctx->PMNI == 0)
        fprintf(ctx->PMResFd,"# c%-3d : constant one\n",l++);
    j = 1;     
    /* the header names the covariates, so it is bounded by the variable
       count, not by nx -- the caller passes nx1, which counts the
       intercept column too and walked one past PMVIdx (valgrind) */
    for (i = 0; i < ctx->PMNV - 1; ++i)   
        fprintf(ctx->PMResFd,"# c%-3d : %s\n",l++,ctx->VName[ctx->PMVIdx[j++]]);
       
    fprintf(ctx->PMResFd,"# c%-3d : %s (dependent)\n",l++,ctx->VName[ctx->PMVIdx[0]]);
    fprintf(ctx->PMResFd,"# c%-3d : predicted\n",l++);
    fprintf(ctx->PMResFd,"# c%-3d : residual\n",l++);
    fprintf(ctx->PMResFd,"\n");

    for (i = 1; i <= n; ++i) {                 /* number of cases */

        fprintf(ctx->PMResFd,"%6d ",i);
        for (j = 1; j <= nx; ++j)
            rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,x[(i - 1) * nx + j]);

        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,y[i]);
        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,y[i] - res[i]);
        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,res[i]);
        fprintf(ctx->PMResFd,"\n");
    }
    printf1(ctx, "\nResiduals written to: %s\n",ctx->PMResFName);
}


