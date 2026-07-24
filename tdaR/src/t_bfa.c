/****************************************************************************/
/*  t_bfa                                                                   */
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
#include "t_ml.h"
#include "t_gdat.h"
#include "t_alloc.h"
#include "t_var.h"
#include "t_gf.h"
#include "t_cdf.h"
#include "t_freq.h"
#include "t_con.h"
#include "tda_context.h"
#include "t_sort.h"


/*  functions in t_bfa.c */

int bfa(TDAContext *ctx);
int bfa_a(TDAContext *ctx, int mmax);
int bfa_free(TDAContext *ctx);
void bfa_p(TDAContext *ctx, unsigned int x,int n,int m);
void bfa_ps(TDAContext *ctx, unsigned int x,unsigned int y,int n,int m,int opt);
void bfa_ps_p(TDAContext *ctx, unsigned int x,unsigned int y,int n,int m);
int bfa_len(TDAContext *ctx, unsigned int x,int n,int m);
int bfa_find(TDAContext *ctx, int m,int n,int maximp);
int bfa_elim(TDAContext *ctx, int m,int n,int nip);
int bfa_chart(TDAContext *ctx, int m,int n,int nip,int maximp);
int bfa_ess(TDAContext *ctx, int nip,int ny1);
int bfa_red(TDAContext *ctx, int nip,int ny1,int *ni,int *nj);
int bfa_pm(TDAContext *ctx, int ni,int nj,int maximp);
int bfa_minc(TDAContext *ctx, int n,char *x, char *y);
int bfa_mince(TDAContext *ctx, int n,char *x, char *y);
int bfa_icheck(TDAContext *ctx, int n,char *x, char *y);
void bfa_copy(TDAContext *ctx, int n,char *x, char *y);
int bfa_fcheck(TDAContext *ctx, unsigned int x,int n,int m,int nei,int ni,char *sel);
void bfa_prn(TDAContext *ctx, int nei,int nf,int ni,int n,int m,int opt);
int bfa_lm(TDAContext *ctx, int ni,int nj);
int bfa_lm1(TDAContext *ctx, int ni,int nj);
int bfc(TDAContext *ctx);
void bfa_ps_n(TDAContext *ctx, unsigned int x,unsigned int y,int n,int m,int jv);


/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

#define IMPMax 50000
#define BMMax 15            /* maximum number of arguments                  */




/* ------------------------------------------------------------------------ */
/*  bfa     Boolean function analysis                                       */
/*                                                                          */
/*  bfa(                                                                    */
/*      opt=...,            treatment of undefined arguments, def. 1        */
/*                          1 treated as don't cares                        */
/*                          2 treated as logical 1                          */
/*                          3 treated as logical 0                          */
/*      alg=...,            algorithm for minimization, def. 0              */
/*                          0 not performed                                 */
/*                          1 = Lawler's algorithm                          */
/*                          2 = Petrick's method                            */
/*                          3 = Lawler's algorithm II                       */
/*      ptab=...,           print table with input data                     */
/*      ptab1=...,          print table with minimized boolean functions    */
/*      prot=...,           protocol file with additional information       */
/*      max=...,            max number of implicants, def. 10000            */
/*                                                                          */
/*  ) = Y,X1,...,Xm;        variables for Boolean function                  */
/*                          (Y is the dependent variable)                   */
/*                                                                          */
/*  Boolean function: Y = f(X1,...,Xm). All variables are interepreted      */
/*  as binary variables. A value not equal to zero will be interpreted      */
/*  as a logical 1.                                                         */
/*                                                                          */
/*  Return 0 if successful, -1 if error.                                    */

int bfa(TDAContext *ctx)
{
    register int i,j,k;
    int err,maximp,m,n,ii,ny,ny0,ny1,ny01,ma,fnd,nip,ni,nj;
    int nei,nni,nf,mlen,nsel;
    unsigned int ui;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Boolean function analysis. Current memory: %d bytes.\n\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 3,4,1))     /* get parameters */
        goto BFAFin;

    if (ctx->PMALG > 3)
        ctx->PMALG = 0;

    if (ctx->PMMax > 0)
        maximp = ctx->PMMax;
    else
        maximp = IMPMax;

    m = ctx->PMNV - 1;
    if (m < 1) {
        printf1(ctx, "Error: need at least two variables.\n");
        goto BFAFin;
    }
    else if (m > BMMax) {
        printf1(ctx, "Error: maximum number of arguments is %d.\n",BMMax);
        goto BFAFin;
    }
    n = (int)pow(2.0,(double)m);

    printf1(ctx, "Number of arguments: %d\n",m);
    printf1(ctx, "Rows of truth table: %d\n",n);

    if (alloc_acn(ctx, n + 1))                 /* if Y=1 */
        goto BFAFin;
    if (alloc_acm(ctx, n + 1))                 /* if Y=0 */
        goto BFAFin;

    for (i = 0; i < ctx->NOC; ++i) {
        ii = 0;
        k  = 1;
        for (j = ctx->PMNV - 1; j >= 1; --j) {
            if ((int)get_data(ctx, ctx->PMVIdx[j],i))
                ii += k;
            k *= 2;
        }
        if ((int)get_data(ctx, ctx->PMVIdx[0],i))
            ctx->AcN[ii] += 1;
        else
            ctx->AcM[ii] += 1;
    }
    ny = ny0 = ny1 = ny01 = 0;
    for (i = 0; i < n; ++i) {
        if (ctx->AcM[i] > 0 && ctx->AcN[i] > 0)
            ny01++;
        else if (ctx->AcM[i] > 0)
            ny0++;
        else if (ctx->AcN[i] > 0)
            ny1++;
        else
            ny++;
    }
    printf1(ctx, "- with y = 0 : %d\n",ny0);
    printf1(ctx, "- with y = 1 : %d\n",ny1);
    printf1(ctx, "- both       : %d\n",ny01);
    printf1(ctx, "- undefined  : %d\n\n",ny);

    if (ny > 0) {
        printf1(ctx, "Undefined truth table rows treated as ");
        if (ctx->PMOPT == 2)
            printf1(ctx, "logical 1.\n");
        else if (ctx->PMOPT == 3)
            printf1(ctx, "logical 0.\n");
        else {
            ctx->PMOPT = 1;
            printf1(ctx, "don't cares.\n");
        }
        if (ctx->PMOPT == 2 || ctx->PMOPT == 3) {
            for (i = 0; i < n; ++i) {
                if (ctx->AcM[i] == 0 && ctx->AcN[i] == 0) {
                    if (ctx->PMOPT == 2)
                        ctx->AcN[i] += 1;
                    else
                        ctx->AcM[i] += 1;
                }
            }
            if (ctx->PMOPT == 2)
                ny1 += ny;
            else
                ny0 += ny;
            ny = 0;
        }
    }

    /* -------------------------------------------------------------------- */
    /* If requested print the truth table to an output file.                */

    if (ctx->PMTabFDef) {            /* print table */
        for (i = 0; i < n; ++i) {
            fprintf(ctx->PMTabFd,"%10d  ",i);

            ui = (unsigned int)(n);
            for (j = 0; j < m; ++j) {
                ui = ui >> 1;
                if (ui & (unsigned int)i)
                    fprintf(ctx->PMTabFd,"1 ");
                else
                    fprintf(ctx->PMTabFd,"0 ");
            }
            fprintf(ctx->PMTabFd,"%3d %3d\n",ctx->AcM[i],ctx->AcN[i]);
        }
        printf1(ctx, "Table written to: %s\n",ctx->PMTabFName);
    }
    if (ctx->PMProtFDef)
        printf1(ctx, "Protocol file: %s\n",ctx->PMProtFName);
    newline(ctx);

    if (ny1 < 1) {
        printf1(ctx, "Number of true minterms is zero.\n");
        err = 0;
        goto BFAFin;
    }
    if (ny01 > 0) {
        printf1(ctx, "Minimization requires unique function values.\n");
        err = 0;
        goto BFAFin;
    }

    if (bfa_a(ctx, maximp))              /* allocate memory */
        goto BFAFin;

    nei = nni = ni = nj = nf = 0;

    /* -------------------------------------------------------------------- */
    /* Find all prime implicants and store in BMImplP.                      */
    /* nip is the number of prime implicants.                               */

    nip = bfa_find(ctx, m,n,maximp);
    if (nip < 0) {
        err = nip;
        goto BFAFin;
    }

    /* If table is incomplete and PMOPT = 1 eliminate redundant
       implicants */

    if (ny > 0 && ctx->PMOPT == 1)
        nip = bfa_elim(ctx, m,n,nip);

    /* print prime implicants */

    printf1(ctx, "Found %d prime implicants.\n",nip);
    for (i = 0; i < nip; ++i) {
        printf1(ctx, "%5d : ",i);
        bfa_ps(ctx, ctx->BMImplP[i],ctx->BMImplP1[i],n,m,0);
        printf1(ctx, "    ");
        bfa_ps(ctx, ctx->BMImplP[i],ctx->BMImplP1[i],n,m,1);
        newline(ctx);
    }
    newline(ctx);

    if (ctx->PMALG == 0) {   /* nothing more done if alg = 0 */
        err = 0;
        goto BFAFin;
    }

    /* -------------------------------------------------------------------- */
    /* Make prime implicant chart in BMImplC, size is nip x ny1             */

    ma = nip * ny1;
    ni = nip;
    nj = ny1;

    printf1(ctx, "Allocating prime implicant chart (%d x %d = %d bytes).\n",nip,ny1,ma);
    if (!(ctx->BMImplC = (char *)calloc((size_t)(ma),sizeof(char)))) {
        p_err(ctx, -2,1);
        goto BFAFin;
    }
    ctx->BMImplC_A = ma;
    memrq(ctx, ma,sizeof(char));

    if ((err = bfa_chart(ctx, m,n,nip,maximp)) < 0)
        goto BFAFin;

    /* -------------------------------------------------------------------- */
    /* Find essential prime implicants. Mark them by a 1 in BMImplX.        */
    /* The number of essential prime implicants is nei.                     */

    nei = bfa_ess(ctx, nip,ny1);

    printf1(ctx, "Found %d essential prime implicants.\n",nei);

    if (ctx->PMProtFDef) {
        fprintf(ctx->PMProtFd,"Prime implicants: %d. Essential: %d.\n",nip,nei);
        fprintf(ctx->PMProtFd,"Prime implicant chart (%d x %d).\n",nip,ny1);

        for (i = 0; i < nip; ++i) {
            fprintf(ctx->PMProtFd,"%5d : ",i);
            bfa_ps_p(ctx, ctx->BMImplP[i],ctx->BMImplP1[i],n,m);
            fprintf(ctx->PMProtFd," : %d : ",ctx->BMImplX[i]);
            for (j = 0; j < ny1; ++j)
                fprintf(ctx->PMProtFd,"%d ",ctx->BMImplC[i * ny1 + j]);
            fprintf(ctx->PMProtFd,"\n");
        }
    }

    /* -------------------------------------------------------------------- */
    /* Store essential prime implicants in BMImplE.                         */

    if (nei > 0) {
        if (!(ctx->BMImplE = (unsigned int *)calloc((size_t)(nei + 1),sizeof(int)))) {
            p_err(ctx, -2,1);
            goto BFAFin;
        }
        ctx->BMImplE_A = nei + 1;
        memrq(ctx, nei + 1,sizeof(int));

        if (!(ctx->BMImplE1 = (unsigned int *)calloc((size_t)(nei + 1),sizeof(int)))) {
            p_err(ctx, -2,1);
            goto BFAFin;
        }
        ctx->BMImplE1_A = nei + 1;
        memrq(ctx, nei + 1,sizeof(int));

        j = 0;
        for (i = 0; i < nip; ++i) {
            if (ctx->BMImplX[i]) {
                ctx->BMImplE[j] = ctx->BMImplP[i];
                ctx->BMImplE1[j] = ctx->BMImplP1[i];
                j++;
            }
        }
    }

    /* -------------------------------------------------------------------- */
    /* Store non-essential prime implicants in BMImplP.                     */
    /* The number of non-essential prime implicants is nni.                 */

    nni = 0;
    for (i = 0; i < nip; ++i) {
        if (ctx->BMImplX[i] == 0) {
            ctx->BMImplP[nni] = ctx->BMImplP[i];
            ctx->BMImplP1[nni] = ctx->BMImplP1[i];
            nni++;
        }
    }
    if (nni <= 0)
        goto BFAFin1;

    /* -------------------------------------------------------------------- */
    /* If nei > 0 make reduced table with only the non-essential prime      */
    /* implicants.                                                          */

    if (nei > 0) {
        if (bfa_red(ctx, nip,ny1,&ni,&nj))
            goto BFAFin;

        printf1(ctx, "Removing essential implicants; reduced chart: %d x %d\n",ni,nj);

        if (ctx->PMProtFDef) {
            fprintf(ctx->PMProtFd,"\nReduced chart (%d x %d).\n",ni,nj);
            if (ni > 0 && nj > 0) {
                for (i = 0; i < ni; ++i) {
                    fprintf(ctx->PMProtFd,"%5d : ",i);
                    bfa_ps_p(ctx, ctx->BMImplP[i],ctx->BMImplP1[i],n,m);
                    fprintf(ctx->PMProtFd," : ");

                    for (j = 0; j < nj; ++j)
                        fprintf(ctx->PMProtFd,"%d ",ctx->BMImplC[i * nj + j]);
                    fprintf(ctx->PMProtFd,"\n");
                }
            }
        }
    }
    if (ni != nni) {
        printf1(ctx, "Error. Can't continue.\n");
        goto BFAFin;
    }
    if (ni <= 0 || nj <= 0)
        goto BFAFin1;

    /* -------------------------------------------------------------------- */
    /* Find a minimal set of the non-essential implicants.                  */
    /* If alg == 0 use Lawler's method                                      */
    /*           1 use Petrick's method.                                    */
    /*                                                                      */
    /* Store the covers in the rows of BMImplS; length in BMImplL.          */
    /* If alg = 1, use maximp as an upper limit for the number of covers.   */
    /* If alg = 0, an upper limit is determined in bfa_lm().                */
    /*                                                                      */
    /* nf is the final number of covers in PMImplS and PMImplL.             */

    if (ctx->PMALG == 1) {
        printf1(ctx, "Find minimal covers with Lawler's method.\n");

        nf = bfa_lm(ctx, ni,nj);      /* Lawler's method */
        if (nf < 0) {
            err = nf;
            goto BFAFin;
        }
    }
    else if (ctx->PMALG == 3) {
        printf1(ctx, "Find minimal covers with Lawler's method II.\n");

        nf = bfa_lm1(ctx, ni,nj);      /* Lawler's method */
        if (nf < 0) {
            err = nf;
            goto BFAFin;
        }
    }
    else if (ctx->PMALG == 2) {
        printf1(ctx, "Find minimal covers with Petrick's method.\n");

        ma = ni * maximp + 1;
        if (!(ctx->BMImplS = (char *)calloc((size_t)(ma),sizeof(char)))) {
            p_err(ctx, -2,1);
            goto BFAFin;
        }
        ctx->BMImplS_A = ma;
        memrq(ctx, ma,sizeof(char));

        ma = maximp + 1;
        if (!(ctx->BMImplL = (short *)calloc((size_t)(ma),sizeof(short)))) {
            p_err(ctx, -2,1);
            goto BFAFin;
        }
        ctx->BMImplL_A = ma;
        memrq(ctx, ma,sizeof(short));

        nf = bfa_pm(ctx, ni,nj,maximp);      /* Petrick's method */
        if (nf < 0) {
            err = nf;
            goto BFAFin;
        }
    }
    if (alloc_ack(ctx, ni + 1))
        goto BFAFin;

    for (i = 0; i < ni; ++i)
        ctx->AcK[i] = bfa_len(ctx, ctx->BMImplP1[i],n,m);

    mlen = ctx->INTMAX;
    for (i = 0; i < nf; ++i) {
        k = 0;
        for (j = 0; j < ni; ++j)
            k += ctx->BMImplS[i * ni + j] * ctx->AcK[j];
        ctx->BMImplL[i] = (short)(k);
        mlen = imin(ctx, k,mlen);
    }
    nsel = 0;
    for (i = 0; i < nf; ++i) {
        if (ctx->BMImplL[i] == mlen) {
            ctx->BMImplL[i] = 1;
            nsel++;
        }
        else
            ctx->BMImplL[i] = 0;
    }
    printf1(ctx, "Number of covers for non-essential implicants: %d\n",nsel);

BFAFin1:

    /* -------------------------------------------------------------------- */
    /* Print final covers.                                                  */

    printf1(ctx, "\nFinal selection of prime implicants.\n");
    bfa_prn(ctx, nei,nf,ni,n,m,0);
    newline(ctx);
    bfa_prn(ctx, nei,nf,ni,n,m,1);


    /* -------------------------------------------------------------------- */
    /* Final check of selections.                                           */

    fnd = 0;
    for (i = 0; i < nf; ++i) {
        if (ctx->BMImplL[i]) {
            for (j = 0; j < n; ++j) {
                k = bfa_fcheck(ctx, (unsigned int)j,n,m,nei,ni,ctx->BMImplS + i * ni);
                if ((ctx->AcN[j] > 0 && k == 0) || (ctx->AcM[j] > 0 && k == 1)) {
                    printf1(ctx, "Final check - error in row %d - selection %d\n",j,i+1);
                    fnd++;
                }
            }
        }
    }
    printf1(ctx, "\nFinal check: %d errors.\n",fnd);
    err = 0;

BFAFin:
    if (err == -2) {
        printf1(ctx, "Exceeded maximum number of implicants (%d).\n",maximp);
        err = -1;
    }
    else if (err == -3) {
        printf1(ctx, "Exceeded maximum levels in bfa_lm().\n");
        err = -1;
    }
    bfa_free(ctx);
    if (ctx->BMImplC_A > 0) {
        free((char *)ctx->BMImplC);
        memrq(ctx, -ctx->BMImplC_A,sizeof(char));
        ctx->BMImplC_A = 0;
    }
    if (ctx->BMImplS_A > 0) {
        free((char *)ctx->BMImplS);
        memrq(ctx, -ctx->BMImplS_A,sizeof(char));
        ctx->BMImplS_A = 0;
    }
    if (ctx->BMImplL_A > 0) {
        free((char *)ctx->BMImplL);
        memrq(ctx, -ctx->BMImplL_A,sizeof(short));
        ctx->BMImplL_A = 0;
    }
    if (ctx->BMImplE_A > 0) {
        free((char *)ctx->BMImplE);
        memrq(ctx, -ctx->BMImplE_A,sizeof(int));
        ctx->BMImplE_A = 0;
    }
    if (ctx->BMImplE1_A > 0) {
        free((char *)ctx->BMImplE1);
        memrq(ctx, -ctx->BMImplE1_A,sizeof(int));
        ctx->BMImplE1_A = 0;
    }
    p_clean(ctx);
    return(err);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  bfa_a   allocated memory. Return 0 if successful, else -1.              */

int bfa_a(TDAContext *ctx, int mmax)
{
    if (!(ctx->BMImplA = (unsigned int *)calloc((size_t)(mmax),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);;
    }
    ctx->BMImplA_A = mmax;
    memrq(ctx, mmax,sizeof(int));

    if (!(ctx->BMImplA1 = (unsigned int *)calloc((size_t)(mmax),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);;
    }
    ctx->BMImplA1_A = mmax;
    memrq(ctx, mmax,sizeof(int));

    if (!(ctx->BMImplB = (unsigned int *)calloc((size_t)(mmax),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);;
    }
    ctx->BMImplB_A = mmax;
    memrq(ctx, mmax,sizeof(int));

    if (!(ctx->BMImplB1 = (unsigned int *)calloc((size_t)(mmax),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);;
    }
    ctx->BMImplB1_A = mmax;
    memrq(ctx, mmax,sizeof(int));

    if (!(ctx->BMImplP = (unsigned int *)calloc((size_t)(mmax),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);;
    }
    ctx->BMImplP_A = mmax;
    memrq(ctx, mmax,sizeof(int));

    if (!(ctx->BMImplP1 = (unsigned int *)calloc((size_t)(mmax),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);;
    }
    ctx->BMImplP1_A = mmax;
    memrq(ctx, mmax,sizeof(int));

    if (!(ctx->BMImplX = (char *)calloc((size_t)(mmax),sizeof(char)))) {
        p_err(ctx, -2,1);
        return(-1);;
    }
    ctx->BMImplX_A = mmax;
    memrq(ctx, mmax,sizeof(char));

    return(0);
}

/* ------------------------------------------------------------------------ */
/*  bfa_free()  free allocated memory.                                      */

int bfa_free(TDAContext *ctx)
{
    if (ctx->BMImplA_A > 0) {
        free((char *)ctx->BMImplA);
        memrq(ctx, -ctx->BMImplA_A,sizeof(int));
        ctx->BMImplA_A = 0;
    }
    if (ctx->BMImplA1_A > 0) {
        free((char *)ctx->BMImplA1);
        memrq(ctx, -ctx->BMImplA1_A,sizeof(int));
        ctx->BMImplA1_A = 0;
    }
    if (ctx->BMImplB_A > 0) {
        free((char *)ctx->BMImplB);
        memrq(ctx, -ctx->BMImplB_A,sizeof(int));
        ctx->BMImplB_A = 0;
    }
    if (ctx->BMImplB1_A > 0) {
        free((char *)ctx->BMImplB1);
        memrq(ctx, -ctx->BMImplB1_A,sizeof(int));
        ctx->BMImplB1_A = 0;
    }
    if (ctx->BMImplP_A > 0) {
        free((char *)ctx->BMImplP);
        memrq(ctx, -ctx->BMImplP_A,sizeof(int));
        ctx->BMImplP_A = 0;
    }
    if (ctx->BMImplP1_A > 0) {
        free((char *)ctx->BMImplP1);
        memrq(ctx, -ctx->BMImplP1_A,sizeof(int));
        ctx->BMImplP1_A = 0;
    }
    if (ctx->BMImplX_A > 0) {
        free((char *)ctx->BMImplX);
        memrq(ctx, -ctx->BMImplX_A,sizeof(char));
        ctx->BMImplX_A = 0;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  bfa_p(x) Print x as a bit pattern.                                      */

void bfa_p(TDAContext *ctx, unsigned int x,int n,int m)
{
    register int j;
    unsigned int ui;

    ui = (unsigned int)(n);
    for (j = 0; j < m; ++j) {
        ui = ui >> 1;
        if (ui & x)
            printf1(ctx, "1");
        else
            printf1(ctx, "0");
    }
}

/* ------------------------------------------------------------------------ */
/*  bfa_ps(x)   Print x as a bit pattern. If opt = 1 print variable names.  */

void bfa_ps(TDAContext *ctx, unsigned int x,unsigned int y,int n,int m,int opt)
{
    register int j,k,l;
    unsigned int ui;

    ui = (unsigned int)(n);
    if (opt == 0) {
        for (j = 0; j < m; ++j) {
            ui = ui >> 1;
            if (ui & (~y))
                printf1(ctx, "-");
            else if (ui & x)
                printf1(ctx, "1");
            else
                printf1(ctx, "0");
        }
    }
    else {
        l = k = 0;
        for (j = 0; j < m; ++j) {
            ui = ui >> 1;
            if (ui & (~y))
                ;
            else {
                l++;
                if (k)
                    printf1(ctx, ".");
                printf1(ctx, "%s",ctx->VName[ctx->PMVIdx[j + 1]]);
                k = 1;
                if (ui & x)
                    ;
                else
                    printf1(ctx, "'");
            }
        }
        if (l == 0)
            tda_out("[empty]");
    }
}

/* ------------------------------------------------------------------------ */
/*  bfa_ps_p(x)   Print x as a bit pattern to the protocol file.            */

void bfa_ps_p(TDAContext *ctx, unsigned int x,unsigned int y,int n,int m)
{
    register int j;
    unsigned int ui;

    if (ctx->PMProtFDef == 0)
        return;

    ui = (unsigned int)(n);
    for (j = 0; j < m; ++j) {
        ui = ui >> 1;
        if (ui & (~y))
            fprintf(ctx->PMProtFd,"-");
        else if (ui & x)
            fprintf(ctx->PMProtFd,"1");
        else
            fprintf(ctx->PMProtFd,"0");
    }
}

/* ------------------------------------------------------------------------ */
/*  bfa_len(x)  Return number of nonzero bits in x.                         */

int bfa_len(TDAContext *ctx, unsigned int x,int n,int m)
{
    (void)ctx;        /* unused: the signature is shared */
    register int j,r;
    unsigned int ui;

    r = 0;
    ui = (unsigned int)(n);
    for (j = 0; j < m; ++j) {
        ui = ui >> 1;
        if (ui & x)
            r++;
    }
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  bfa_find()  Find all implicants. The prime implicants are stored in     */
/*              BMImplP and BMImplP1. Return the number of implicants.      */
/*              Return -2 if exceeded storage for max number of implicants. */

int bfa_find(TDAContext *ctx, int m,int n,int maximp)
{
    register int i,j,k;
    int d,nia,nib,nip,fnd;
    unsigned int ai,aj,ax,xx,am;

    /* begin with implicants of order 0 */

    k = 0;
    for (i = 0; i < n; ++i) {
        if (ctx->AcM[i] == 0) {
            if (k >= maximp)
                return(-2);
            ctx->BMImplA[k] = (unsigned int)(i);
            k++;
        }
    }
    nia = k;
    for (i = 0; i < nia; ++i)
        ctx->BMImplA1[i] = 0xffffff;

    nip = 0;    /* number of prime implicants */

    for (d = 1; d <= m; ++d) {  /* d is order of implicant */

        tda_err("Level %d of %d (nia=%d, nip=%d)\n",d,m,nia,nip);
        tda_err_flush();

        if (nia == 0)
            break;

        nib = 0;
        for (i = 1; i < nia; ++i) {
            ai = ctx->BMImplA[i];
            ax = ctx->BMImplA1[i];
            ai = ai & ax;

            for (j = 0; j < i; ++j) {
                if (ax != ctx->BMImplA1[j])
                    continue;

                aj = ctx->BMImplA[j] & ax;
                xx = ai ^ aj;
                fnd = 0;
                for (k = 0; k < m; ++k) {
                    if (xx & 0x01) {
                        if (xx == 0x01)
                            fnd = 1;
                        break;
                    }
                    xx = xx >> 1;
                }
                if (fnd) {
                    ctx->BMImplX[i] = 1;     /* mark i and j as used */
                    ctx->BMImplX[j] = 1;

                    /* check if previously found */

                    am = ax & (~(ai ^ aj));

                    for (k = 0; k < nib; ++k) {
                        if (ctx->BMImplB1[k] == am && (ctx->BMImplB[k] & am) == (ai & am)) {
                            fnd = 0;
                            break;
                        }
                    }
                    if (fnd) {
                        if (nib >= maximp)
                           return(-2);
                        ctx->BMImplB[nib] = ai;
                        ctx->BMImplB1[nib] = am;
                        nib++;
                    }
                }
            }
        }

        /* save prime implicants  */

        for (i = 0; i < nia; ++i) {
            if (ctx->BMImplX[i] == 0) {
                if (nip >= maximp)
                    return(-2);
                ctx->BMImplP[nip] = ctx->BMImplA[i];
                ctx->BMImplP1[nip] = ctx->BMImplA1[i];
                nip++;
            }
        }
        for (i = 0; i < nib; ++i) {
            ctx->BMImplA[i] = ctx->BMImplB[i];
            ctx->BMImplA1[i] = ctx->BMImplB1[i];
            ctx->BMImplX[i] = 0;
        }
        nia = nib;
        nib = 0;
    }
    if (nia > 0) {
        for (i = 0; i < nia; ++i) {
            if (nip >= maximp)
                return(-2);
            ctx->BMImplP[nip] = ctx->BMImplA[i];
            ctx->BMImplP1[nip] = ctx->BMImplA1[i];
            nip++;
        }
    }
    return(nip);
}

/* ------------------------------------------------------------------------ */
/*  bfa_elim()  Eliminate redundant implicants from BMImplP, and BMImplP1.  */
/*              Return nip = new number of prime implicants.                */

int bfa_elim(TDAContext *ctx, int m,int n,int nip)
{
    (void)m;        /* unused: the signature is shared */
    register int i,j,k;
    unsigned int ax,ai;

    j = 0;
    for (i = 0; i < nip; ++i) {
        ai = ctx->BMImplP[i];
        ax = ctx->BMImplP1[i];
        ai = ai & ax;
        for (k = 0; k < n; ++k) {
            if (ctx->AcN[k] > 0 && ((unsigned int)k & ax) == ai) {
                ctx->BMImplP[j] = ctx->BMImplP[i];
                ctx->BMImplP1[j] = ctx->BMImplP1[i];
                j++;
                break;
            }
        }
    }
    return(j);
}

/* ------------------------------------------------------------------------ */
/*  bfa_chart() Make prime implicant chart in BMImplC (nip x ny1)           */
/*              Return  0 if successful.                                    */
/*                     -1 if insufficient memory.                           */
/*                     -2 if exceeded storage for max number of implicants. */

int bfa_chart(TDAContext *ctx, int m,int n,int nip,int maximp)
{
    (void)m;        /* unused: the signature is shared */
    register int i,j,k;
    int nia;
    unsigned int ai,ax,am;

    k = 0;
    for (i = 0; i < n; ++i) {
        if (ctx->AcM[i] == 0 && ctx->AcN[i] > 0) {
            if (k >= maximp)
                return(-2);
            ctx->BMImplA[k] = (unsigned int)(i);
            k++;
        }
    }
    nia = k;
    for (i = 0; i < nip; ++i) {
        ai = ctx->BMImplP[i];
        ax = ctx->BMImplP1[i];
        am = ai & ax;
        for (j = 0; j < nia; ++j) {
            if (am == (ctx->BMImplA[j] & ax))
                ctx->BMImplC[i * nia + j] = 1;
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  bfa_ess()   Find essential prime implicants form prime implicant chart. */
/*              Mark them by a 1 in BMImplX (i = 0,...,nip-1).              */
/*              Return number of essential implicants.                      */

int bfa_ess(TDAContext *ctx, int nip,int ny1)
{
    register int i,j,k,l;
    int nei;

    nei = 0;
    for (i = 0; i < nip; ++i)
        ctx->BMImplX[i] = 0;

    for (j = 0; j < ny1; ++j) {
        k = 0;
        l = -1;
        for (i = 0; i < nip; ++i) {
            if (ctx->BMImplC[i * ny1 + j]) {
                k++;
                if (l < 0)
                    l = i;
            }
        }
        if (k == 1 && l >= 0)
            ctx->BMImplX[l] = 1;
    }
    for (i = 0; i < nip; ++i) {
        if (ctx->BMImplX[i])
            nei++;
    }
    return(nei);
}

/* ------------------------------------------------------------------------ */
/*  bfa_red()   Reduce the prime implicant chart by deleting the essential  */
/*              implicants. It is assumed that BMImplX still marks the      */
/*              essential prime implicants.                                 */
/*                                                                          */
/*              Return 0 if successful (and ni and nj of reduced table)     */
/*                    -1 if insufficient memory.                            */

int bfa_red(TDAContext *ctx, int nip,int ny1,int *ni,int *nj)
{
    register int i,j,ii,jj;

    if (alloc_acj(ctx, ny1 + 1))
        return(-1);

    for (j = 0; j < ny1; ++j) {
        for (i = 0; i < nip; ++i) {
            if (ctx->BMImplX[i] && ctx->BMImplC[i * ny1 + j]) {
                ctx->AcJ[j] = 1;
                break;
            }
        }
    }
    *ni = *nj = 0;
    for (i = 0; i < nip; ++i) {
        if (ctx->BMImplX[i] == 0)
            *ni += 1;
    }
    for (j = 0; j < ny1; ++j) {
        if (ctx->AcJ[j] == 0)
            *nj += 1;
    }
    if (*ni == 0 || *nj == 0)
        return(0);

    ii = 0;
    for (i = 0; i < nip; ++i) {
        if (ctx->BMImplX[i] == 0) {
            jj = 0;
            for (j = 0; j < ny1; ++j) {
                if (ctx->AcJ[j] == 0) {
                    ctx->BMImplC[ii * *nj + jj] = ctx->BMImplC[i * ny1 + j];
                    jj++;
                }
            }
            ii++;
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  bfa_pm()    Enumeration (Petrick's method). The algorithm follows       */
/*              E.L. Lawler, Covering Problems: Duality Relations and a     */
/*              New Method of Solution, SIAM Journal on Applied             */
/*              Mathematics 14 (1966), 1115-32.                             */
/*                                                                          */
/*              Prime implicant chart is in BMImplC with ni rows and nj     */
/*              columns. The functions uses BMImplS with ni rows and        */
/*              maximp columns for storage of selections.                   */
/*                                                                          */
/*              Return -1 if insufficient memory.                           */
/*                     -2 if exceeded storage for max number of implicants. */

int bfa_pm(TDAContext *ctx, int ni,int nj,int maximp)
{
    register int i,j,k,l;
    int is,isel,ii,jj,fnd;
    char *xptr;

    if (alloc_acj(ctx, nj + 1))
        return(-1);
    if (alloc_ack(ctx, nj + 1))
        return(-1);

    for (i = 0; i < ni; ++i) {
        k = i * nj;
        for (j = 0; j < nj; ++j) {
            if (ctx->BMImplC[k + j])
                ctx->AcJ[j] += 1;
        }
    }
    if (sortdpi(ctx, nj,ctx->AcJ,ctx->AcK))
        return(-1);

    is = isel = 0;
    for (jj = 0; jj < nj; ++jj) {
        j = ctx->AcK[jj];
        /****************
        tda_err("column %d of %d\n",jj,nj);
        tda_err_flush();
        *******************/
        if (isel == 0) {
            for (i = 0; i < ni; ++i) {
                if (ctx->BMImplC[i * nj + j]) {
                    if (isel >= maximp)
                        return(-2);
                    for (k = 0; k < ni; ++k) {
                        if (k == i)
                            ctx->BMImplS[isel * ni + k] = 1;
                        else
                            ctx->BMImplS[isel * ni + k] = 0;
                    }
                    isel++;
                }
            }
            is = isel;
        }
        else {
            is = isel;
            for (i = 0; i < ni; ++i) {
                if (ctx->BMImplC[i * nj + j]) {
                    for (l = 0; l < isel; ++l) {
                        if (is >= maximp)
                            return(-2);

                        for (k = 0; k < ni; ++k) {
                            if (k == i)
                                ctx->BMImplS[is * ni + k] = 1;
                            else
                                ctx->BMImplS[is * ni + k] = ctx->BMImplS[l * ni + k];
                        }
                        is++;
                    }
                }
            }
        }

        /* find minimal elements */

        if (isel == is)
            continue;

        for (l = isel; l < is; ++l) {
            ctx->BMImplL[l] = 1;
            xptr = ctx->BMImplS + l * ni;
            for (k = isel; k < is; ++k) {
                if (k != l && bfa_minc(ctx, ni,ctx->BMImplS + k * ni,xptr)) {
                    ctx->BMImplL[l] = 0;
                    break;
                }
            }
        }
        ii = 0;
        for (l = isel; l < is; ++l) {
            if (ctx->BMImplL[l]) {
                fnd = 0;
                for (k = 0; k < ii; ++k) {  /* does vector already exist? */
                    if (bfa_icheck(ctx, ni,ctx->BMImplS + k * ni,ctx->BMImplS + l * ni)) {
                        fnd = 1;
                        break;
                    }
                }
                if (fnd == 0) {
                    for (k = 0; k < ni; ++k)
                        ctx->BMImplS[ii * ni + k] = ctx->BMImplS[l * ni + k];
                    ii++;
                }
            }
        }
        isel = ii;
    }
    fnd = ctx->INTMAX;
    for (l = 0; l < isel; ++l) {
        ii = 0;
        for (k = 0; k < ni; ++k) {
            if (ctx->BMImplS[l * ni + k])
                ii++;
        }
        ctx->BMImplL[l] = (short)(ii);
        fnd = imin(ctx, fnd,ii);
    }
    j = 0;
    for (i = 0; i < isel; ++i) {
        if (ctx->BMImplL[i] > fnd)
            continue;
        if (i > j) {
            bfa_copy(ctx, ni,ctx->BMImplS + j * ni,ctx->BMImplS + i * ni);
            ii++;
        }
        j++;
    }
    return(j);
}

/* ------------------------------------------------------------------------ */
/*  bfa_minc(n,x,y)     Return 1 if x is less than y, otherwise 0.          */

int bfa_minc(TDAContext *ctx, int n,char *x, char *y)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,m;

    m = 0;
    for (i = 0; i < n; ++i) {
        if (x[i] > y[i])
            return(0);
        else if (x[i] < y[i])
            m++;
    }
    if (m)
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  bfa_mince(n,x,y)    Return 1 if x is less or equal than y, otherwise 0. */

int bfa_mince(TDAContext *ctx, int n,char *x, char *y)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;

    for (i = 0; i < n; ++i) {
        if (x[i] > y[i])
            return(0);
    }
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  bfa_icheck(n,x,y)   Return 1 if x and y are identical, otherwise 0.     */

int bfa_icheck(TDAContext *ctx, int n,char *x, char *y)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;

    for (i = 0; i < n; ++i) {
        if (x[i] != y[i])
            return(0);
    }
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  bfa_copy(n,x,y)     Copy y into x.                                      */

void bfa_copy(TDAContext *ctx, int n,char *x, char *y)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;
    for (i = 0; i < n; ++i)
        x[i] = y[i];
}

/* ------------------------------------------------------------------------ */
/*  bfa_fcheck()    Check final selections. nei is the number of essential  */
/*                  implicants in BMImplE.                                  */
/*                                                                          */

int bfa_fcheck(TDAContext *ctx, unsigned int x,int n,int m,int nei,int ni,char *sel)
{
    (void)m; (void)n;        /* unused: the signature is shared */
    register int i;

    for (i = 0; i < nei; ++i) {
        if ((x & ctx->BMImplE1[i]) == (ctx->BMImplE[i] & ctx->BMImplE1[i]))
            return(1);
    }
    for (i = 0; i < ni; ++i) {
        if (sel[i]) {
            if ((x & ctx->BMImplP1[i]) == (ctx->BMImplP[i] & ctx->BMImplP1[i]))
                return(1);
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  bfa_prn()   Print final selections.                                     */

#ifdef TDA_R_PACKAGE
/* one row per term of each final selection: selection number, then a
   code per argument variable: 1, 0, or -1 for don't-care -- the same
   bit walk bfa_ps prints.  Emitted only for the opt==0 (pattern) call
   so the symbolic pass doesn't double it. */
static void bfa_export_term(TDAContext *ctx, int sel, unsigned int x,
                            unsigned int y, int n, int m)
{
    int j;
    unsigned int ui = (unsigned int)n;
    tda_export_cell(ctx, "bfa.selection", (double)sel);
    for (j = 0; j < m; ++j) {
        ui = ui >> 1;
        if (ui & (~y))
            tda_export_cell(ctx, "bfa.selection", -1.0);
        else if (ui & x)
            tda_export_cell(ctx, "bfa.selection", 1.0);
        else
            tda_export_cell(ctx, "bfa.selection", 0.0);
    }
    tda_export_endrow(ctx, "bfa.selection");
}
#endif

void bfa_prn(TDAContext *ctx, int nei,int nf,int ni,int n,int m,int opt)
{
    register int i,j,k,l;

    if (nf == 0) {
        for (j = 0; j < nei; ++j) {
            if (j)
                printf1(ctx, " + ");
            bfa_ps(ctx, ctx->BMImplE[j],ctx->BMImplE1[j],n,m,opt);
#ifdef TDA_R_PACKAGE
            if (opt == 0)
                bfa_export_term(ctx, 1, ctx->BMImplE[j],ctx->BMImplE1[j],n,m);
#endif
        }
        newline(ctx);
        return;
    }
    k = 0;
    for (i = 0; i < nf; ++i) {
        if (ctx->BMImplL[i]) {
            l = 0;
            printf1(ctx, "%3d  ",++k);
            for (j = 0; j < nei; ++j) {
                if (l)
                    printf1(ctx, " + ");
                bfa_ps(ctx, ctx->BMImplE[j],ctx->BMImplE1[j],n,m,opt);
#ifdef TDA_R_PACKAGE
                if (opt == 0)
                    bfa_export_term(ctx, k, ctx->BMImplE[j],ctx->BMImplE1[j],n,m);
#endif
                l = 1;
            }
            for (j = 0; j < ni; ++j) {
                if (ctx->BMImplS[i * ni + j]) {
                    if (l)
                        printf1(ctx, " + ");
                    bfa_ps(ctx, ctx->BMImplP[j],ctx->BMImplP1[j],n,m,opt);
#ifdef TDA_R_PACKAGE
                    if (opt == 0)
                        bfa_export_term(ctx, k, ctx->BMImplP[j],ctx->BMImplP1[j],n,m);
#endif
                    l = 1;
                }
            }
            newline(ctx);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  bfa_lm()    Find minimal covers. The algorithm is due to                */
/*              E.L. Lawler, Covering Problems: Duality Relations and a     */
/*              New Method of Solution, SIAM Journal on Applied             */
/*              Mathematics 14 (1966), 1115-32.                             */
/*                                                                          */
/*              Prime implicant chart is in BMImplC with ni rows and nj     */
/*              columns. The functions uses BMImplS with ni rows and        */
/*              maximp columns for storage of selections.                   */
/*                                                                          */
/*              Return -1 if insufficient memory.                           */
/*                     -2 if exceeded storage for max number of implicants. */
/*                     -3 if maxlev (or ma) insufficient.                   */

int bfa_lm(TDAContext *ctx, int ni,int nj)
{
    register int i,j,k,l;
    int ii,ik,jk,cmin,nej,ne,nlev,maxlev,ma,nf,jmin;

int ncnt;


    maxlev = nj + 10;
    nf = 0;

    if (alloc_acl(ctx, ni * maxlev + 1))
        return(-1);
    if (alloc_aci(ctx, maxlev + 1))
        return(-1);
    if (alloc_acns(ctx, ni * nj + 1))
        return(-1);
    if (alloc_acms(ctx, ni + 1))
        return(-1);
    if (alloc_acr(ctx, ni + 1))
        return(-1);
    if (alloc_acs(ctx, ni + 1))
        return(-1);
    if (alloc_acc(ctx, ni + 1))
        return(-1);
    if (alloc_acd(ctx, ni + 1))
        return(-1);
    if (alloc_ace(ctx, ni * nj + 1))
        return(-1);
    if (alloc_acf(ctx, ni * nj + 1))
        return(-1);

    k = ni * nj;                /* use AcF as local copy of BMImplC */
    for (i = 0; i < k; ++i)
        ctx->AcF[i] = ctx->BMImplC[i];
    nej = nj;
    nlev = 0;

    while (nej > 0) {
tda_out("nej=%d\n",nej);            /* select column from AcF */
/**
        for (i = 0; i < ni; ++i) {
            tda_out("AcF  ");
            for (j = 0; j < nej; ++j)
                tda_out("%d ",AcF[i * nj + j]);
            newline(ctx);
        }
**/



/*********

        if (alloc_acj(ctx, nej + 1))
            return(-1);
        if (alloc_ack(ctx, nej + 1))
            return(-1);

        for (i = 0; i < ni; ++i) {
            k = i * nj;
            for (j = 0; j < nej; ++j) {
                if (AcF[k + j])
                    AcJ[j] += 1;
            }
        }
        if (sortdpi(ctx, nej,AcJ,AcK))
            return(-1);
        j = AcK[0];
**********/

        jmin = ctx->INTMAX;
        for (j = 0; j < nej; ++j) {
            k = 0;
            for (i = 0; i < ni; ++i) {
                if (ctx->AcF[i * nj + j])
                    k++;
            }
            jmin = imin(ctx, jmin,k);
        }
        tda_out("JMIN=%d\n",jmin);

        for (j = 0; j < nej; ++j) {
            k = 0;
            for (i = 0; i < ni; ++i) {
                if (ctx->AcF[i * nj + j])
                    k++;
            }
            if (k == jmin)
                break;
        }

tda_out("Using j=%d nej=%d ni=%d nj=%d    \n",j,nej,ni,nj  );


        ik = 0;
        for (i = 0; i < ni; ++i) {
            ii = i * nj;
            if (ctx->AcF[ii + j] == 0)
                continue;

            if (nlev >= maxlev)
                return(-3);

            ctx->AcL[nlev * ni + ik] = i;
            ctx->AcI[nlev] += 1;

            jk = 0;
            for (k = 0; k < nej; ++k) {
                if (ctx->AcF[ii + k] == 0) {
                    ctx->AcNS[ik * nj + jk] = (short)(k);
                    jk++;
                }
            }
            ctx->AcMS[ik] = (short)(jk);
            ik++;
        }
        nlev++;
        if (nej <= 1)
            break;

        ne = 0;
        for (i = 0; i < ik; ++i)
            ctx->AcR[i] = 0;

ncnt = 0;

        while (1) {

            ncnt++;
                                      /* create or'ed row */
            for (i = 0; i < ni; ++i)
                ctx->AcC[i] = 0;
            for (k = 0; k < ik; ++k) {
                l = ctx->AcNS[k * nj + ctx->AcR[k]];
                for (i = 0; i < ni; ++i)
                    ctx->AcC[i] = ctx->AcC[i] | ctx->AcF[i * nj + l];
            }

            /* --------------------------------------------- */

            cmin = 1;       /* assume AcC is minimal */
            for (i = 0; i < ik; ++i)
                ctx->AcS[i] = 0;

            while (1) {
                for (i = 0; i < ni; ++i)
                    ctx->AcD[i] = 0;
                for (k = 0; k < ik; ++k) {
                    l = ctx->AcNS[k * nj + ctx->AcS[k]];
                    for (i = 0; i < ni; ++i)
                        ctx->AcD[i] = ctx->AcD[i] | ctx->AcF[i * nj + l];
                }
                if (bfa_minc(ctx, ni,ctx->AcD,ctx->AcC)) {
                    cmin = 0;
                    break;
                }
                for (k = 0; k < ik; ++k) {
                    ctx->AcS[k] += 1;
                    if (ctx->AcS[k] < ctx->AcMS[k])
                        break;
                    ctx->AcS[k] = 0;
                }
                if (k >= ik)
                    break;
            }

            /* --------------------------------------------- */

            if (cmin) {        /* save minimal column in AcE */
                for (i = 0; i < ni; ++i)
                    ctx->AcE[i * nj + ne] = ctx->AcC[i];
                ne++;
            }
            for (k = 0; k < ik; ++k) {
                ctx->AcR[k] += 1;
                if (ctx->AcR[k] < ctx->AcMS[k])
                    break;
                ctx->AcR[k] = 0;
            }
            if (k >= ik)
                break;
        }
tda_out("ncnt=%d\n",ncnt);


                                     /* copy AcE to AcF */
        for (i = 0; i < ni; ++i) {
            for (j = 0; j < ne; ++j)
                ctx->AcF[i * nj + j] = ctx->AcE[i * nj + j];
        }
        nej = ne;
    }

    tda_out("LEVEL nlev=%d\n",nlev  );
    for (i = 0; i < nlev; ++i) {
        tda_out("AcI: %d ",ctx->AcI[i]);
        for (j = 0; j < ctx->AcI[i]; ++j)
            tda_out("%d ",ctx->AcL[i * ni + j]);
        newline(ctx);
    }











    ma = 1;
    for (i = 0; i < nlev; ++i)
        ma *= ctx->AcI[i];

    if (!(ctx->BMImplS = (char *)calloc((size_t)(ma) * (size_t)(ni) + 1,sizeof(char)))) {
        p_err(ctx, -2,1);
        return(-1);
    }
    ctx->BMImplS_A = ma * ni + 1;
    memrq(ctx, ma * ni + 1,sizeof(char));

    if (!(ctx->BMImplL = (short *)calloc((size_t)(ma + 1),sizeof(short)))) {
        p_err(ctx, -2,1);
        return(-1);
    }
    ctx->BMImplL_A = ma + 1;
    memrq(ctx, ma + 1,sizeof(short));

    if (alloc_acr(ctx, nlev + 1))
        return(-1);

    while (1) {
        ii = 0;
        for (j = 0; j < nj; ++j) {
            ii = 0;
            for (i = 0; i < nlev; ++i) {
                k = ctx->AcL[i * ni + ctx->AcR[i]];
                ii += ctx->BMImplC[k * nj + j];
            }
            if (ii == 0)
                break;
        }
        if (ii) {
            if (nf >= ma)
                return(-3);

            for (i = 0; i < ni; ++i)
                ctx->BMImplS[nf * ni + i] = 0;

            for (i = 0; i < nlev; ++i) {
                k = ctx->AcL[i * ni + ctx->AcR[i]];
                ctx->BMImplS[nf * ni + k] = 1;
            }
            nf++;
        }
        for (i = 0; i < nlev; ++i) {
            ctx->AcR[i] += 1;
            if (ctx->AcR[i] < ctx->AcI[i])
                break;
            ctx->AcR[i] = 0;
        }
        if (i >= nlev)
            break;
    }
    return(nf);
}


/* ------------------------------------------------------------------------ */
/*  bfa_lm1()   Find minimal covers. The algorithm is due to                */
/*              E.L. Lawler, Covering Problems: Duality Relations and a     */
/*              New Method of Solution, SIAM Journal on Applied             */
/*              Mathematics 14 (1966), 1115-32.                             */
/*                                                                          */
/*              Prime implicant chart is in BMImplC with ni rows and nj     */
/*              columns. The functions uses BMImplS with ni rows and        */
/*              maximp columns for storage of selections.                   */
/*                                                                          */
/*              Return -1 if insufficient memory.                           */
/*                     -2 if exceeded storage for max number of implicants. */
/*                     -3 if maxlev (or ma) insufficient.                   */

/* ------------------------------------------------------------------------ */
/*  bfa_minc1(n,x,y)     Return 1 if x is less than y, otherwise 0.          */

int bfa_minc1(TDAContext *ctx, int n,char *x, char *y)
{
    register int i,m;

    tda_out("minc1 n=%d\n",n);
    tda_out("x ");
    for (i = 0; i < n; ++i)
        tda_out("%d ",x[i]);
    tda_out("\ny ");
    for (i = 0; i < n; ++i)
        tda_out("%d ",y[i]);
    newline(ctx);

    m = 0;
    for (i = 0; i < n; ++i) {
tda_out("%d (%d,%d)\n",i,x[i],y[i]);


        if (x[i] > y[i]) {
tda_out(" hier return 0\n");
            return(0);
        }
        else if (x[i] < y[i])
            m++;
    }
tda_out(" hier m=%d\n",m );
    if (m)
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  bfa_lm1()   Find minimal covers. The algorithm is due to                */
/*              E.L. Lawler, Covering Problems: Duality Relations and a     */
/*              New Method of Solution, SIAM Journal on Applied             */
/*              Mathematics 14 (1966), 1115-32.                             */
/*                                                                          */
/*              Prime implicant chart is in BMImplC with ni rows and nj     */
/*              columns. The functions uses BMImplS with ni rows and        */
/*              maximp columns for storage of selections.                   */
/*                                                                          */
/*              Return -1 if insufficient memory.                           */
/*                     -2 if exceeded storage for max number of implicants. */
/*                     -3 if maxlev (or ma) insufficient.                   */

int bfa_lm1(TDAContext *ctx, int ni,int nj)
{
    register int i,j,k,l;
    int ii,jj,ik,ll,jk,nej,ne,nlev,maxlev,ma,nf,nm,njmax,nec,maxnf;
    char x;

    njmax = 2 * nj;
    maxlev = nj + 10;
    maxnf = 200000;

    nf = 0;

    if (alloc_acl(ctx, ni * maxlev + 1))
        return(-1);
    if (alloc_aci(ctx, maxlev + 1))
        return(-1);
    if (alloc_ack(ctx, njmax + 1))
        return(-1);
    if (alloc_acns(ctx, ni * njmax + 1))
        return(-1);
    if (alloc_acs(ctx, ni + 1))
        return(-1);
    if (alloc_acr(ctx, ni + 1))
        return(-1);
    if (alloc_acc(ctx, ni + 1))
        return(-1);
    if (alloc_ace(ctx, ni * njmax + 1))
        return(-1);
    if (alloc_acf(ctx, ni * njmax + 1))
        return(-1);

                                /* use AcF as local copy of BMImplC */
    for (i = 0; i < ni; ++i) {
        for (j = 0; j < nj; ++j)
            ctx->AcF[i * njmax + j] = ctx->BMImplC[i * nj + j];
    }
    nej = nj;
    ne = nlev = 0;

    while (nej > 0) {
tda_out("nej=%d ni=%d \n",nej,ni  );            /* select column from AcF */

        for (i = 0; i < ni; ++i) {
            tda_out("AcF  ");
            for (j = 0; j < nej; ++j)
                tda_out("%d ",ctx->AcF[i * njmax + j]);
            newline(ctx);
        }

        ii = ctx->INTMAX;
        for (j = 0; j < nej; ++j) {
            k = 0;
            for (i = 0; i < ni; ++i) {
                if (ctx->AcF[i * njmax + j])
                    k++;
            }
            ii = imin(ctx, ii,k);
        }
        for (j = 0; j < nej; ++j) {
            k = 0;
            for (i = 0; i < ni; ++i) {
                if (ctx->AcF[i * njmax + j])
                    k++;
            }
            if (k == ii)
                break;
        }

tda_out("Using j=%d nej=%d ni=%d nj=%d ii=%d    \n",j,nej,ni,nj,ii);


        ik = 0;
        for (i = 0; i < ni; ++i) {
            ii = i * njmax;
            if (ctx->AcF[ii + j] == 0)
                continue;

            if (nlev >= maxlev)     {
tda_out("nlev=%d maxlev=%d\n",nlev,maxlev);

                return(-3);
}
            ctx->AcL[nlev * ni + ik] = i;
            ctx->AcI[nlev] += 1;

            jk = 0;
            for (k = 0; k < nej; ++k) {
                if (ctx->AcF[ii + k] == 0) {
                    ctx->AcNS[ik * njmax + jk] = (short)(k);
                    jk++;
                }
            }
            ctx->AcS[ik] = jk;
            ik++;
        }
        j = 0;
        for (i = 0; i < ik; ++i) {
            ii = ctx->AcS[i];
            if (ii > 0) {
                for (k = 0; k < ii; ++k)
                    ctx->AcNS[j * njmax + k] = ctx->AcNS[i * njmax + k];
                j++;
            }
        }
/****************

        if (j < ik) {
            ik = j;

            tda_out("Reduced AcNS\n");
            for (i = 0; i < ik; ++i) {
                tda_out("AcS=%3d : ",AcS[i]);
                for (k = 0; k < AcS[i]; ++k)
                    tda_out("%3d ",AcNS[i * njmax + k]);
                newline(ctx);
            }
        }
        ik = j;
***********/


        if (ik == 0)
            break;


        nlev++;
        if (nej <= 1)
            break;

        /* Create array with column numbers for or'ed columns */
        /* AcMS[nm x ik] for column numbers */
        /* AcD[nm] to index non-minimal elements */

        /* first find the number of equal column numbers: nec and save */
        /* in AcK */

        nec = 0;
        nm = 0;
        for (i = 0; i < ik; ++i) {
            ctx->AcR[i] = 0;
            if (ctx->AcS[i])
                nm = 1;
        }
        if (nm != 0) {

            while (1) {
                l = -1;
                ii = 0;
                for (k = 0; k < ik; ++k) {
                    if (ctx->AcS[k] > 0) {
                        if (l < 0) {
                            l = ctx->AcNS[k * njmax + ctx->AcR[k]];
                            ii = 1;
                        }
                        else if (l != ctx->AcNS[k * njmax + ctx->AcR[k]]) {
                            ii = 0;
                            break;
                        }
                    }
                }
                if (ii)
                    ctx->AcK[nec++] = l;

                for (k = 0; k < ik; ++k) {
                    ctx->AcR[k] += 1;
                    if (ctx->AcR[k] < ctx->AcS[k])
                        break;
                    ctx->AcR[k] = 0;
                }
                if (k >= ik)
                    break;
            }

            tda_out("nec=%d\nAcK: ",nec );
            for (k = 0; k < nec; ++k)
                tda_out("%d ",ctx->AcK[k]);
            newline(ctx);

            /* reorganize AcS and AcNS */

            for (i = 0; i < ik; ++i) {
                jj = ctx->AcS[i];
                ll = 0;
                for (k = 0; k < jj; ++k) {
                    ii = ctx->AcNS[i * njmax + k];
                    nf = 1;
                    for (l = 0; l < nec; ++l) {
                        if (ii == ctx->AcK[l]) {
                            nf = 0;
                            break;
                        }
                    }
                    if (nf) { /* not found */
                        ctx->AcNS[i * njmax + ll] = (short)(ii);
                        ll++;
                    }
                }
                ctx->AcS[i] = ll;
            }

            tda_out("NEW AcNS\n");
            for (i = 0; i < ik; ++i) {
                tda_out("AcS=%3d : ",ctx->AcS[i]);
                for (k = 0; k < ctx->AcS[i]; ++k)
                    tda_out("%3d ",ctx->AcNS[i * njmax + k]);
                newline(ctx);
            }
            nm = 0;
            for (i = 0; i < ik; ++i) {
                if (ctx->AcS[i])
                    nm = 1;
            }
        }

tda_out("nm==========%d ik=%d \n",nm,ik);

        /* Save minimal vectors from equal column numbers */

        ne = 0;
        for (j = 0; j < nec; ++j) {
            k = ctx->AcK[j];
            if (ne >= njmax) {
                printf1(ctx, "Fatal error in bfa_lm() [njmax=%d ne=%d].\n",njmax,ne);
                return(-1);
            }
            for (i = 0; i < ni; ++i)
                ctx->AcE[ne * ni + i] = ctx->AcF[i * njmax + k];
            ne++;
        }


tda_out("actual ne=%d\n",ne);

        /* If nm > 0 find remaining minimal vectors */

        if (nm != 0) {

            nm = 1;
            for (i = 0; i < ik; ++i) {
                if (ctx->AcS[i])
                    nm *= ctx->AcS[i];
            }

tda_out("Actual nm=%d\n",nm);


            if (alloc_acms(ctx, nm * ik + 1))
                return(-1);
            if (alloc_acd(ctx, nm + 1))
                return(-1);

            for (i = 0; i < ik; ++i)
                ctx->AcR[i] = 0;

            j = 0;
            while (1) {
                for (k = 0; k < ik; ++k) {
                    if (ctx->AcS[k])
                        ctx->AcMS[j * ik + k] = ctx->AcNS[k * njmax + ctx->AcR[k]];
                    else
                        ctx->AcMS[j * ik + k] = -1;
                }
                j++;
                for (k = 0; k < ik; ++k) {
                    ctx->AcR[k] += 1;
                    if (ctx->AcR[k] < ctx->AcS[k])
                        break;
                    ctx->AcR[k] = 0;
                }
                if (k >= ik)
                    break;
            }
            /* Find minimal elements. AcD[] = 1 if not minimal */

            for (j = 0; j < nm; ++j) {
                for (i = 0; i < ni; ++i)
                    ctx->AcC[i] = 0;

                for (k = 0; k < ik; ++k) {
                    if ((l = ctx->AcMS[j * ik + k]) >= 0) {
                        for (i = 0; i < ni; ++i)
                            ctx->AcC[i] |= ctx->AcF[i * njmax + l];
                    }
                }
                for (i = 0; i < ne; ++i) {
                    if (bfa_mince(ctx, ni,ctx->AcE + i * ni,ctx->AcC)) {
                        ctx->AcD[j] = 1;
                        break;
                    }
                }
                if (ctx->AcD[j])
                    continue;

                for (ll = j + 1; ll < nm; ++ll) {
                    ii = 0;
                    for (i = 0; i < ni; ++i) {
                        x = 0;
                        for (k = 0; k < ik; ++k) {
                            if ((l = ctx->AcMS[ll * ik + k]) >= 0)
                                x |= ctx->AcF[i * njmax + l];
                        }
                        if (ctx->AcC[i]) {
                            if (x == 0)
                                ii++;
                        }
                        else if (x) {
                            ii = 0;
                            break;
                        }
                    }
                    if (ii) {  /* l is less than j */
                        ctx->AcD[j] = 1;
                        break;
                    }
                }
                if (ctx->AcD[j] == 0) {      /* found new minimal vector */
                    if (ne >= njmax) {
                        printf1(ctx, "Fatal error in bfa_lm() [njmax=%d ne=%d].\n",njmax,ne);
                        return(-1);
                    }
                    for (i = 0; i < ni; ++i)
                        ctx->AcE[ne * ni + i] = ctx->AcC[i];
                    ne++;
                }
            }
            tda_out("NEUES ne=%d\n",ne);
        }

        for (i = 0; i < ni; ++i) {
            for (j = 0; j < ne; ++j)
                ctx->AcF[i * njmax + j] = ctx->AcE[j * ni + i];
        }
        nej = ne;


    }

    tda_out("LEVEL nlev=%d\n",nlev  );
    for (i = 0; i < nlev; ++i) {
        tda_out("AcI: %d ",ctx->AcI[i]);
        for (j = 0; j < ctx->AcI[i]; ++j)
            tda_out("%d ",ctx->AcL[i * ni + j]);
        newline(ctx);
    }

    /* Simplify levels */

    for (i = 1; i < nlev; ++i) {
        ii = ctx->AcI[i];
        for (j = 0; j < i; ++j) {
            jj = ctx->AcI[j];
            ll = 0;
            for (k = 0; k < jj; ++k) {
                ma = ctx->AcL[j * ni + k];
                for (l = 0; l < ii; ++l) {
                    if (ma == ctx->AcL[i * ni + l]) {
                        ll++;
                        break;
                    }
                }
            }
            if (ll == jj) {
                for (k = 0; k < jj; ++k) {
                    ma = ctx->AcL[j * ni + k];
                    for (l = 0; l < ii; ++l) {
                        if (ma == ctx->AcL[i * ni + l]) {
                            ctx->AcL[i * ni + l] = -1;
                            break;
                        }
                    }
                }
            }
        }
        j = 0;
        for (l = 0; l < ii; ++l) {
            if ((k = ctx->AcL[i * ni + l]) >= 0)
                ctx->AcL[i * ni + j++] = k;
        }
        ctx->AcI[i] = j;
    }
    tda_out("LEVEL nlev=%d\n",nlev  );
    for (i = 0; i < nlev; ++i) {
        tda_out("AcI: %d ",ctx->AcI[i]);
        for (j = 0; j < ctx->AcI[i]; ++j)
            tda_out("%d ",ctx->AcL[i * ni + j]);
        newline(ctx);
    }



    if (!(ctx->BMImplS = (char *)calloc((size_t)(maxnf) * (size_t)(ni) + 1,sizeof(char)))) {
        p_err(ctx, -2,1);
        return(-1);
    }
    ctx->BMImplS_A = maxnf * ni + 1;
    memrq(ctx, maxnf * ni + 1,sizeof(char));

    if (!(ctx->BMImplL = (short *)calloc((size_t)(maxnf + 1),sizeof(short)))) {
        p_err(ctx, -2,1);
        return(-1);
    }
    ctx->BMImplL_A = maxnf + 1;
    memrq(ctx, maxnf + 1,sizeof(short));

    if (alloc_acr(ctx, nlev + 1))
        return(-1);

    nf = 0;
    while (1) {
        ii = 0;
        for (j = 0; j < nj; ++j) {
            ii = 0;
            for (i = 0; i < nlev; ++i) {
                k = ctx->AcL[i * ni + ctx->AcR[i]];
                ii += ctx->BMImplC[k * nj + j];
            }
            if (ii == 0)
                break;
        }
        if (ii) {
            if (nf >= maxnf)
                return(-3);

            for (i = 0; i < ni; ++i)
                ctx->BMImplS[nf * ni + i] = 0;

            for (i = 0; i < nlev; ++i) {
                k = ctx->AcL[i * ni + ctx->AcR[i]];
                ctx->BMImplS[nf * ni + k] = 1;
            }
            ii = 0;
            for (k = 0; k < nf; ++k) {
                ii = 1;
                for (i = 0; i < ni; ++i) {
                    if (ctx->BMImplS[nf * ni + i] != ctx->BMImplS[k * ni + i]) {
                        ii = 0;
                        break;
                    }
                }
                if (ii)     /* already found */
                    break;
            }
            if (ii == 0)
                nf++;
        }
        for (i = 0; i < nlev; ++i) {
            ctx->AcR[i] += 1;
            if (ctx->AcR[i] < ctx->AcI[i])
                break;
            ctx->AcR[i] = 0;
        }
        if (i >= nlev)
            break;
    }

    tda_out("final nf=%d\n",nf);


    return(nf);
}

/* ------------------------------------------------------------------------ */
/*  bfc     Context-dependencies in Boolean functions                       */
/*                                                                          */
/*  bfc(                                                                    */
/*      opt=...,            treatment of undefined arguments, def. 1        */
/*                          1 treated as don't cares                        */
/*                          2 treated as logical 1                          */
/*                          3 treated as logical 0                          */
/*      prot=...,           protocol file with additional information       */
/*      max=...,            max number of implicants, def. 10000            */
/*                                                                          */
/*  ) = Y,X1,...,Xm;        variables for Boolean function                  */
/*                          (Y is the dependent variable)                   */
/*                                                                          */
/*  Boolean function: Y = f(X1,...,Xm). All variables are interepreted      */
/*  as binary variables. A value not equal to zero will be interpreted      */
/*  as a logical 1.                                                         */
/*                                                                          */
/*  Return 0 if successful, -1 if error.                                    */

int bfc(TDAContext *ctx)
{
    register int i,j,k;
    int err,maximp,m,n,m1,jv,ii,nip;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Context-dependencies in Boolean functions. Current memory: %d bytes.\n\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 3,4,1))     /* get parameters */
        goto BFCFin;

    if (ctx->PMMax > 0)
        maximp = ctx->PMMax;
    else
        maximp = IMPMax;

    m = ctx->PMNV - 1;
    if (m < 2) {
        printf1(ctx, "Error: need at least three variables.\n");
        goto BFCFin;
    }
    else if (m > BMMax) {
        printf1(ctx, "Error: maximum number of arguments is %d.\n",BMMax);
        goto BFCFin;
    }
    m1 = m - 1;
    n = (int)pow(2.0,(double)m1);

    printf1(ctx, "Number of arguments: %d\n",m);
    printf1(ctx, "Rows of reduced truth table: %d\n",n);
    if (ctx->PMProtFDef)
        printf1(ctx, "Protocol file: %s\n",ctx->PMProtFName);
    newline(ctx);

    if (bfa_a(ctx, maximp))              /* allocate memory */
        goto BFCFin;

    /* ---------------------------------------------------------------- */
    /* Do for all independent variables: ...                            */

    for (jv = 1; jv <= m; ++jv) {
        printf1(ctx, "Variable: %s\n",ctx->VName[ctx->PMVIdx[jv]]);

        /* AcN: X = 0 and Y = 0 */
        /* AcM: X = 0 and Y = 1 */
        /* AcI: X = 1 and Y = 0 */
        /* AcJ: X = 1 and Y = 1 */
        /* AcR = 0 if not known
                 1 if positive effect
                 2 if negative effect
                 3 if no effect
                 4 if inconsistent */

        if (alloc_acn(ctx, n + 1))
            goto BFCFin;
        if (alloc_acm(ctx, n + 1))
            goto BFCFin;
        if (alloc_aci(ctx, n + 1))
            goto BFCFin;
        if (alloc_acj(ctx, n + 1))
            goto BFCFin;
        if (alloc_acr(ctx, n + 1))
            goto BFCFin;

        /* create reduced truth table */

        for (i = 0; i < ctx->NOC; ++i) {
            ii = 0;
            k  = 1;
            for (j = ctx->PMNV - 1; j >= 1; --j) {
                if (j == jv)
                    continue;

                if ((int)get_data(ctx, ctx->PMVIdx[j],i))
                    ii += k;
                k *= 2;
            }
            if ((int)get_data(ctx, ctx->PMVIdx[0],i)) {
                if ((int)get_data(ctx, ctx->PMVIdx[jv],i))
                    ctx->AcJ[ii] += 1;
                else
                    ctx->AcM[ii] += 1;
            }
            else {
                if ((int)get_data(ctx, ctx->PMVIdx[jv],i))
                    ctx->AcI[ii] += 1;
                else
                    ctx->AcN[ii] += 1;
            }
        }
        ii = 0;

        for (i = 0; i < n; ++i) {
            if (ctx->AcN[i] > 0 && ctx->AcM[i] > 0) {
                ctx->AcR[i] = 4;
                ii = 1;
            }
            else if (ctx->AcI[i] > 0 && ctx->AcJ[i] > 0) {
                ctx->AcR[i] = 4;
                ii = 1;
            }
            else if (ctx->AcN[i] > 0 && ctx->AcJ[i] > 0)
                ctx->AcR[i] = 1;
            else if (ctx->AcM[i] > 0 && ctx->AcI[i] > 0)
                ctx->AcR[i] = 2;
            else if (ctx->AcN[i] > 0 && ctx->AcI[i] > 0)
                ctx->AcR[i] = 3;
            else if (ctx->AcM[i] > 0 && ctx->AcJ[i] > 0)
                ctx->AcR[i] = 3;
        }
        if (ctx->PMProtFDef) {
            fprintf(ctx->PMProtFd,"Variable: %s\n",ctx->VName[ctx->PMVIdx[jv]]);
            for (i = 0; i < n; ++i) {
                fprintf(ctx->PMProtFd,"%4d ",i);
                fprintf(ctx->PMProtFd,"%4d ",ctx->AcN[i]);
                fprintf(ctx->PMProtFd,"%4d ",ctx->AcM[i]);
                fprintf(ctx->PMProtFd,"%4d ",ctx->AcI[i]);
                fprintf(ctx->PMProtFd,"%4d  ",ctx->AcJ[i]);
                if (ctx->AcR[i] == 0)
                    fprintf(ctx->PMProtFd,"(not known)\n");
                else if (ctx->AcR[i] == 1)
                    fprintf(ctx->PMProtFd,"(positive effect)\n");
                else if (ctx->AcR[i] == 2)
                    fprintf(ctx->PMProtFd,"(negative effect)\n");
                else if (ctx->AcR[i] == 3)
                    fprintf(ctx->PMProtFd,"(no effect)\n");
                else
                    fprintf(ctx->PMProtFd,"(inconsistent)\n");
            }
        }
        if (ii) {
            printf1(ctx, "-- data are inconsistent.\n\n");
            continue;
        }

        for (ii = 0; ii <= 3; ++ii) {
            if (ii == 0)
                printf1(ctx, "Not known: ");
            else if (ii == 1)
                printf1(ctx, "Positive effect: ");
            else if (ii == 2)
                printf1(ctx, "Negative effect: ");
            else
                printf1(ctx, "No effect: ");

            for (i = 0; i < n; ++i) {
                if (ctx->AcR[i] == ii)
                    ctx->AcM[i] = 0;
                else
                    ctx->AcM[i] = 1;
            }
            nip = bfa_find(ctx, m1,n,maximp);
            if (nip < 0) {
                err = nip;
                goto BFCFin;
            }
            else if (nip == 0)
                printf1(ctx, "none\n");
            else {
                for (i = 0; i < nip; ++i) {
                    if (i)
                         printf1(ctx, ", ");
                    bfa_ps_n(ctx, ctx->BMImplP[i],ctx->BMImplP1[i],n,m1,jv-1);
                }
                newline(ctx);
            }
        }
        newline(ctx);
    }

    err = 0;

BFCFin:
    if (err == -2) {
        printf1(ctx, "Exceeded maximum number of implicants (%d).\n",maximp);
        err = -1;
    }
    else if (err == -3) {
        printf1(ctx, "Found inconsistent data (check with bfa).\n");
        err = -1;
    }
    bfa_free(ctx);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  bfa_ps_n(x) Print variable names according to x, skip jv.               */

void bfa_ps_n(TDAContext *ctx, unsigned int x,unsigned int y,int n,int m,int jv)
{
    register int j,k,l,j1;
    unsigned int ui;

    ui = (unsigned int)(n);
    l = k = 0;
    for (j = 0; j < m; ++j) {
        ui = ui >> 1;
        if (ui & (~y))
            ;
        else {
            l++;
            if (k)
                printf1(ctx, ".");

            if (j < jv)
                j1 = j;
            else
                j1 = j + 1;

            printf1(ctx, "%s",ctx->VName[ctx->PMVIdx[j1 + 1]]);
            k = 1;
            if (ui & x)
                ;
            else
                printf1(ctx, "'");
        }
    }
    if (l == 0)
        tda_out("always");
}
