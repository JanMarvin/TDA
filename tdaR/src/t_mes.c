/****************************************************************************/
/*  t_mes                                                                   */
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
#include "t_ml.h"
#include "t_lsei.h"
#include "t_gio.h"
#include "t_com.h"
#include "t_lp.h"
#include "tda_context.h"

/*  functions in t_mes.c */

int ghd(TDAContext *ctx);
void ghd_init(TDAContext *ctx);
void visit_ghd(TDAContext *ctx, int n,int k,int k1,int k2,double val,int first);
int ghd_d(TDAContext *ctx, int m,double *xi,double *xj);

int ghd1(TDAContext *ctx);
void ghd1_fnd(TDAContext *ctx, int n,int l,int nm,int ni,int *num,int *idx,int *lev);
int ghd1_h(TDAContext *ctx, int n,int ni,int *num,int *lev);

int tda_conj(TDAContext *ctx);
void conj_free(TDAContext *ctx);
int conj_r(TDAContext *ctx);
int conj_tab(TDAContext *ctx);
int conj_dmat(TDAContext *ctx, int opt);
int conj_com(void);
int conj_it(TDAContext *ctx);
int conj_eval(TDAContext *ctx);
void conj_prn(TDAContext *ctx);
void conj_res(TDAContext *ctx);
int conj_lp(TDAContext *ctx);

int mreg(TDAContext *ctx);
int mreg_it(TDAContext *ctx, int nx1);

int mreg_proj(TDAContext *ctx, int n,double *x);
int mreg_proj1(TDAContext *ctx, int n,int m,int *t,int *p,int *pp,double *x,double *y);
int mreg_proj2(TDAContext *ctx, int n,int m,int *t,double *x,double *y);
void mreg_ms(TDAContext *ctx, int n,double *y,double *m,double *s);
void mreg_std(TDAContext *ctx, int n,int nx,double *x,double *beta,double *z,double m,double s);
int mreg_ties(TDAContext *ctx, int n,double *y,int *t);

int nmca(TDAContext *ctx);
int nmca_com(TDAContext *ctx);






/* ------------------------------------------------------------------------ */
/*  ghd         Create Hasse-Diagram.                                       */
/*                                                                          */
/*              ghd(                                                        */
/*                  df=...,         edge list                               */
/*                  pcf=...,        dot file                                */
/*                  ptab=...,       file with interval-valued ranks         */
/*              ) = X1,...,Xm;                                              */
/*                                                                          */
/*  It is assumed that the (X1,...,Xm) patterns are unique.                 */
/*  The ptab option uses a version of gep().                                */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int ghd(TDAContext *ctx)
{
    register int i,j,k;
    int err,n,m,a,nn,nn1 = 0,ne,n2,na,nb,nna,nnb,la,lb,lmin,lmax;
    double x,xa,xb;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Hasse diagram. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 3,4,1))       /* get parameters */
        goto GHDFin;

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 8,4);

    n = ctx->NOC;
    n2 = n + 2;                         /* plus two additional nodes */
    m = ctx->PMNV;                           /* number of variables */

    if (alloc_acx(ctx, n2 * m + 1))          /* input data */
        goto GHDFin;
    if (alloc_ack(ctx, n + 1))               /* node list */
        goto GHDFin;

    for (i = 0; i < n; ++i) {           /* get data into AcX */
        for (j = 0; j < m; ++j)
            ctx->AcX[i * m + j] = get_data(ctx, ctx->PMVIdx[j],i);
    }

    /* create list of nodes with zero indegree */

    for (i = 0; i < n; ++i) {
        a = 1;
        for (k = 0; k < n; ++k) {
            if (k != i && ghd_d(ctx, m,ctx->AcX + k * m,ctx->AcX + i * m)) {
                a = 0;
                break;
            }
        }
        ctx->AcK[i] = a;
    }

    printf1(ctx, "\nNodes with zero indegree:");
    for (i = 0; i < n; ++i) {
        if (ctx->AcK[i])
            printf1(ctx, " %d",i + 1);
    }
    newline(ctx);

    if (ctx->PMPCFDef) {
        fprintf(ctx->PMPCFd,"digraph g {\n");
        for (i = 1; i <= n; ++i)
            fprintf(ctx->PMPCFd,"n%d [label=%d];\n",i,i);
        nn1 = n + 1;
    }
    ne = nn = 0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (j == i || ghd_d(ctx, m,ctx->AcX + i * m,ctx->AcX + j * m) == 0)
                continue;
            a = 1;
            for (k = 0; k < n; ++k) {
                if (k == i || k == j)
                    continue;
                if (ghd_d(ctx, m,ctx->AcX + i * m,ctx->AcX + k * m) &&
                    ghd_d(ctx, m,ctx->AcX + k * m,ctx->AcX + j * m)) {
                    a = 0;
                    break;
                }
            }
            if (a) {
                if (ctx->PMF1Def) {
                    fprintf(ctx->PMF1d,"%5d %5d\n",i + 1,j + 1);
                    nn++;
                }
                if (ctx->PMPCFDef) {
                    fprintf(ctx->PMPCFd,"n%d -> n%d;\n",i+1,j+1);
                    nn1++;
                }
                ne++;
            }
        }
    }
    printf1(ctx, "Hasse diagram: %d nodes, %d edges.\n",n,ne);

    if (ctx->PMF1Def)
        printf1(ctx, "%d records written to: %s\n",nn,ctx->PMF1dName);
    if (ctx->PMPCFDef) {
        fprintf(ctx->PMPCFd,"}\n");
        nn1++;
        printf1(ctx, "%d records written to: %s\n",nn1,ctx->PMPCFName);
    }
    newline(ctx);

    if (ctx->PMTabFDef == 0) {
        err = 0;
        goto GHDFin;
    }
    printf1(ctx, "Creating interval-valued ranks.\n");

    /***
    for (i = 0; i < n; ++i) {
        for (j = 0; j < m; ++j)
            tda_out("%g ",AcX[i * m + j]);
        newline(ctx);
    }
    ***/

    /*  create minimal and maximal elements in x[n] and x[n+1] */

    na = n;
    nb = n + 1;

    for (j = 0; j < m; ++j) {
        xa = xb = ctx->AcX[j];
        for (i = 1; i < n; ++i) {
            x = ctx->AcX[i * m + j];
            xa = dmin(ctx, xa,x);
            xb = dmax(ctx, xb,x);
        }
        ctx->AcX[na * m + j] = xa;
        ctx->AcX[nb * m + j] = xb;
    }
    nna = nnb = -1;

    for (i = 0; i < n; ++i) {
        a = 1;
        for (j = 0; j < m; ++j) {
            if (ctx->AcX[i * m + j] != ctx->AcX[na * m + j]) {
                a = 0;
                break;
            }
        }
        if (a) {
            nna = i;
            break;
        }
    }
    for (i = 0; i < n; ++i) {
        a = 1;
        for (j = 0; j < m; ++j) {
            if (ctx->AcX[i * m + j] != ctx->AcX[nb * m + j]) {
                a = 0;
                break;
            }
        }
        if (a) {
            nnb = i;
            break;
        }
    }
    nn = n;
    printf1(ctx, "Minimal element: %g",ctx->AcX[na * m]);
    for (j = 1; j < m; ++j)
        printf1(ctx, ",%g",ctx->AcX[na * m + j]);
    if (nna < 0) {
        nna = na;
        nn++;
        printf1(ctx, " (not contained, new node number: %d)\n",nna + 1);
    }
    else
        printf1(ctx, " (contained, node number: %d)\n",nna + 1);

    printf1(ctx, "Maximal element: %g",ctx->AcX[nb * m]);
    for (j = 1; j < m; ++j)
        printf1(ctx, ",%g",ctx->AcX[nb * m + j]);
    if (nnb < 0) {
        nnb = nn;
        if (nnb < nb) {
            for (j = 0; j < m; ++j)
                ctx->AcX[nnb * m + j] = ctx->AcX[nb * m + j];
        }
        nn++;
        printf1(ctx, " (not contained, new node number: %d)\n",nnb + 1);
    }
    else
        printf1(ctx, " (contained, node number: %d)\n",nnb + 1);

    /* create adjacency matrix in AcD[] */

    if (alloc_acd(ctx, nn * nn + 1))
        goto GHDFin;

    ne = 0;
    for (i = 0; i < nn; ++i) {
        for (j = 0; j < nn; ++j) {
            if (j == i || ghd_d(ctx, m,ctx->AcX + i * m,ctx->AcX + j * m) == 0)
                continue;
            a = 1;
            for (k = 0; k < nn; ++k) {
                if (k == i || k == j)
                    continue;
                if (ghd_d(ctx, m,ctx->AcX + i * m,ctx->AcX + k * m) &&
                    ghd_d(ctx, m,ctx->AcX + k * m,ctx->AcX + j * m)) {
                    a = 0;
                    break;
                }
            }
            if (a) {
                ne++;
                ctx->AcD[i * nn + j] = 1;
            }
        }
    }
    printf1(ctx, "New Hasse diagram: %d nodes, %d edges.\n",nn,ne);

    /***
    for (i = 0; i < nn; ++i) {
        for (j = 0; j < nn; ++j)
            tda_out("%d ",AcD[i * nn + j]);
        newline(ctx);
    }
    ***/

    /* find length of longest path from nna to nnb */

    if (alloc_acc(ctx, nn + 1))
        goto GHDFin;

    for (k = 0; k < nn; ++k)
        ctx->AcC[k] = 0;

    ghd_init(ctx);
    visit_ghd(ctx, nn,nna,nna,nnb,0.0,1);
    lmin = ctx->GSCON_MINLEN;
    lmax = ctx->GSCON_MAXLEN;

    /***
    tda_out("np=%d nsp=%d len=%d minlen=%d maxlen=%d val=%g  minval=%g maxval=%g cnt=%d\n",
        GSCON_NP,GSCON_NSP,GSCON_LEN,GSCON_MINLEN,GSCON_MAXLEN,
        GSCON_VAL,GSCON_MINVAL,GSCON_MAXVAL,GSCON_CNT);
    ***/

    printf1(ctx, "\nNumber of paths from node %d to %d: %d\n",nna+1,nnb+1,ctx->GSCON_NP);
    printf1(ctx, "Minimal length: %d, maximal length: %d\n",lmin,lmax);

    /** printf1("\nNode  MinLenght  MaxLength  Interval of ranks\n"); **/

    nn1 = 0;
    for (k = 0; k < nn; ++k) {
        if (k == nna || k == nnb)
            continue;

        ghd_init(ctx);
        visit_ghd(ctx, nn,nna,nna,k,0.0,1);
        la = ctx->GSCON_MINLEN;

        ghd_init(ctx);
        visit_ghd(ctx, nn,k,k,nnb,0.0,1);
        lb = ctx->GSCON_MINLEN;

        /** printf1("%4d  %9d  %9d  [%7d,%7d]\n",k+1,la,lb,la,lmax - lb); **/

        fprintf(ctx->PMTabFd,"%4d  %9d  %9d  %8d %8d\n",k+1,la,lb,la,lmax - lb);
        nn1++;
    }
    printf1(ctx, "%d records written to: %s\n\n",nn1,ctx->PMTabFName);

    err = 0;

GHDFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ghd_init()  Initialize values for visit_ghd().                          */

void ghd_init(TDAContext *ctx)
{
    ctx->GSCON_NP = 0;           /* # of paths */
    ctx->GSCON_NSP = 0;          /* # of shortest paths */
    ctx->GSCON_LEN = -1;         /* length of actual paths */
    ctx->GSCON_MINLEN = 0;       /* length of shortest paths */
    ctx->GSCON_MAXLEN = 0;       /* length of longest paths */
    ctx->GSCON_VAL = 0.0;        /* value of actual paths */
    ctx->GSCON_MINVAL = 0.0;     /* min value of paths */
    ctx->GSCON_MAXVAL = 0.0;     /* max value of paths */
    ctx->GSCON_CNT = 0;
}

/* ------------------------------------------------------------------------ */
/*  visit_ghd   called by ghd.                                              */

void visit_ghd(TDAContext *ctx, int n,int k,int k1,int k2,double val,int first)
{
    register int j;
    double tmp;

    /* printf1("visit_ghd k=%d k1=%d k2=%d first=%d\n",k,k1,k2,first); */

    if (first == 0 || k != k1 || k != k2)
        ctx->AcC[k] = 1;
    ctx->GSCON_CNT++;
    ctx->GSCON_LEN++;
    ctx->GSCON_VAL += val;

    if (k == k2 && first == 0) {

        if (k1 == k2 && ctx->GSCON_LEN <= 2)     /* for cycles min length 3 */
            goto V4CONT;

        ctx->GSCON_NP++;
        if (ctx->GSCON_MINVAL <= 0.0 || ctx->GSCON_MINVAL > ctx->GSCON_VAL) {
            ctx->GSCON_MINVAL = ctx->GSCON_VAL;
            ctx->GSCON_MINLEN = ctx->GSCON_LEN;
            ctx->GSCON_NSP = 1;
        }
        else if (fabs(ctx->GSCON_MINVAL - ctx->GSCON_VAL) < ctx->EPSI1) {
            ctx->GSCON_NSP++;
            if (ctx->GSCON_MINLEN > ctx->GSCON_LEN)
                ctx->GSCON_MINLEN = ctx->GSCON_LEN;
        }
        if (ctx->GSCON_MAXVAL <= 0.0 || ctx->GSCON_MAXVAL < ctx->GSCON_VAL) {
            ctx->GSCON_MAXVAL = ctx->GSCON_VAL;
            ctx->GSCON_MAXLEN = ctx->GSCON_LEN;
        }
        else if (fabs(ctx->GSCON_MAXVAL - ctx->GSCON_VAL) < ctx->EPSI1) {
            if (ctx->GSCON_MAXLEN < ctx->GSCON_LEN)
                ctx->GSCON_MAXLEN = ctx->GSCON_LEN;
        }

V4CONT:
        ctx->GSCON_LEN--;
        ctx->GSCON_VAL -= val;
        ctx->AcC[k] = 0;
        ctx->GSCON_CNT--;
        return;
    }
    for (j = 0; j < n; ++j) {
        if (ctx->AcD[k * n + j]) {          /* edge from k to j */
            if (ctx->AcC[j] == 0) {
                tmp = 1.0;
                if (tmp >= 0.0)
                    visit_ghd(ctx, n,j,k1,k2,tmp,0);
            }
        }
    }

/*******
    n0 = GD_FPN[k];
    if (n0 > 0) {
        for (l = 0; l < n0; ++l) {
            j = GD_FPI[k][l];
            if (AcC[j] == 0) {
                tmp = gdd_adj(ctx, k,j,PMGN);
                if (tmp >= 0.0)
                    visit_ghd(ctx, j,k1,k2,tmp,0);
            }
        }
    }
**/


    ctx->GSCON_LEN--;
    ctx->GSCON_VAL -= val;
    if (ctx->GSCON_CNT <= 0) {
        printfe(ctx, "ERROR in visit_ghd.\n");
        gerr_exit(ctx, 216);
    }
    ctx->AcC[k] = 0;
    ctx->GSCON_CNT--;
}

/* ------------------------------------------------------------------------ */
/*  ghd_d(n,m,i,j,x)                                                        */
/*                                                                          */
/*  Return: 1 if xi <= xj, otherwise 0.                                     */

int ghd_d(TDAContext *ctx, int m,double *xi,double *xj)
{
    (void)ctx;        /* unused: the signature is shared */
    register int k;

    for (k = 0; k < m; ++k) {
        if (xi[k] > xj[k])
            return(0);
    }
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  ghd1        Find maximal subgraph in Hasse-Diagram.                     */
/*              Experimental, should not be used.                           */
/*                                                                          */
/*              ghd1(                                                       */
/*                  df=...,         edge list                               */
/*                  pcf=...,        dot file                                */
/*              ) = X1,...,Xm;                                              */
/*                                                                          */
/*  It is assumed that the (X1,...,Xm) patterns are unique.                 */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int ghd1(TDAContext *ctx)
{
    register int i,j,k;
    int err,n,m,a,ne,nm;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Find maximal hierarchy. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,4,1))       /* get parameters */
        goto GHD1Fin;

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 8,4);

    n = ctx->NOC;
    m = ctx->PMNV;                           /* number of variables */

    if (alloc_acx(ctx, n * m + 1))           /* input data */
        goto GHD1Fin;
    if (alloc_acd(ctx, n * n + 1))           /* adjacency matrix */
        goto GHD1Fin;
    if (alloc_ack(ctx, n + 1))               /* flags node numbers */
        goto GHD1Fin;
    if (alloc_acn(ctx, n + 1))               /* node numbers */
        goto GHD1Fin;
    if (alloc_acs(ctx, n + 1))               /* level */
        goto GHD1Fin;
    if (alloc_acm(ctx, n + 1))               /* final selection */
        goto GHD1Fin;

    for (i = 0; i < n; ++i) {           /* get data into AcX */
        for (j = 0; j < m; ++j)
            ctx->AcX[i * m + j] = get_data(ctx, ctx->PMVIdx[j],i);
    }

    ne = 0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (j == i || ghd_d(ctx, m,ctx->AcX + i * m,ctx->AcX + j * m) == 0)
                continue;
            a = 1;
            for (k = 0; k < n; ++k) {
                if (k == i || k == j)
                    continue;
                if (ghd_d(ctx, m,ctx->AcX + i * m,ctx->AcX + k * m) &&
                    ghd_d(ctx, m,ctx->AcX + k * m,ctx->AcX + j * m)) {
                    a = 0;
                    break;
                }
            }
            if (a) {
                ctx->AcD[i * n + j] = 1;
                ne++;
            }
        }
    }
    printf1(ctx, "Hasse diagram: %d nodes, %d edges.\n",n,ne);

    ctx->GHD1Max = 1;
    nm = 0;
    for (k = 0; k < n; ++k) {
        ctx->AcN[0] = k;
        ctx->AcK[k] = 1;

        ghd1_fnd(ctx, n,1,nm,1,ctx->AcN,ctx->AcK,ctx->AcS);

        printf1(ctx, "\nGHD1Max=%d : ",ctx->GHD1Max);
        for (i = 0; i < ctx->GHD1Max; ++i)
            printf1(ctx, "%d ",ctx->AcM[i]);
        newline(ctx);

        break;

        ctx->AcK[k] = -1;
        nm++;
    }




    err = 0;

GHD1Fin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ghd1_fnd()                                                              */
/*                                                                          */

void ghd1_fnd(TDAContext *ctx, int n,int l,int nm,int ni,int *num,int *idx,int *lev)
{
    register int i,j;
    int l1;

    if (n - nm <= ctx->GHD1Max)
        return;

    /* add one more node */

    l1 = l + 1;

    for (i = 0; i < n; ++i) {
        if (idx[i])
            continue;

        num[ni] = i;
        idx[i] = l1;
        if (ghd1_h(ctx, n,ni + 1,num,lev) == 0) {
            idx[i] = -l1;
            nm++;
            continue;
        }
        if (ni + 1 > ctx->GHD1Max) {
            printf1(ctx, "GHD1Max=%d ni+1=%d\n",ctx->GHD1Max,ni + 1);
            ctx->GHD1Max = ni + 1;
            for (j = 0; j <= ni; ++j) {
                ctx->AcM[j] = num[j];
                printf1(ctx, "%d ",ctx->AcM[j]);
            }
            newline(ctx);
        }
        ghd1_fnd(ctx, n,l1,nm,ni + 1,num,idx,lev);
    }
    for (i = 0; i < n; ++i) {
        if (iabs(ctx, idx[i]) == l1)
            idx[i] = 0;
    }
}

/* ------------------------------------------------------------------------ */
/*  ghd1_h(n,ni,num,lev)    Return 1 if the nodes in num[] can be           */
/*                          hierarchically partitioned, otherwise 0.        */

int ghd1_h(TDAContext *ctx, int n,int ni,int *num,int *lev)
{
    register int i,j,in,jn,l;
    int a;

    if (ni <= 2)
        return(1);

    for (i = 0; i < ni; ++i)
        lev[i] = 0;

    for (j = 0; j < ni; ++j) {
        jn = num[j];
        a = 1;
        for (i = 0; i < ni; ++i) {
            in = num[i];
            if (in == jn)
                continue;
            if (ctx->AcD[in * n + jn]) {
                a = 0;
                break;
            }
        }
        if (a)
            lev[j] = 1;
    }
    for (l = 2; l <= ni; ++l) {
        for (i = 0; i < ni; ++i) {
            in = num[i];
            if (lev[i] == l - 1) {
                for (j = 0; j < ni; ++j) {
                    jn = num[j];
                    if (ctx->AcD[in * n + jn]) {
                        if (lev[j] > 0 && lev[j] != l)
                            return(0);
                        lev[j] = l;
                    }
                }
            }
        }
    }
    for (i = 0; i < ni; ++i) {
        if (lev[i] == 0)
            return(0);
    }
    return(1);
}


/* -##--------------------------------------------------------------------- */
/*  tda_conj        Conjoint analysis (alternating least squares).              */
/*                                                                          */
/*              tda_conj(                                                       */
/*                  alg=...,        algorithm, def. 1                       */
/*                                  1 = monotone regression                 */
/*                                  2 = linear programming                  */
/*                  opt=...,        1 or 2 or 3, def. 1                     */
/*                                  1 = ties ignored                        */
/*                                  2 = Kruskal's primary approach          */
/*                                  3 = Kruskal's secondary approach        */
/*                  xp=...,         starting values                         */
/*                  dsv=...,        input file with starting values         */
/*                  mxit=...,       max number of iterations, def. 50       */
/*                  tolf=...,       tolerance for stress, def. 1e-6         */
/*                  tolsp=...,      tolerance for par change, def. 1.e-6    */
/*                  pmat=...,       print design matrix                     */
/*                  nfmt = ...,     integer format, def. 2                  */
/*                  tfmt=...,       print format coefficients, 10.4         */
/*                  pres=...,       write rank order and estimates          */
/*                  fmt=...,        format for pres, def. 8.4               */
/*              ) = R,X1,...,Xm;                                            */
/*                                                                          */
/*  R contains rank ordering. Floating point values.                        */
/*  X1,...,Xm are integer valued.                                           */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int tda_conj(TDAContext *ctx)
{
    int err;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Conjoint analysis. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLF = 1.e-6;
    ctx->TOLSP = 1.e-6;

    if (parm(ctx, ctx->CmdBuf + 4,4,1))     /* get parameters */
        goto CONJFin;

    if (ctx->PMALG != 2)
        ctx->PMALG = 1;

    if (ctx->PMOPT > 3)
        ctx->PMOPT = 1;

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 8,4);

    if (ctx->PMNFmtF == 0) {
        ctx->PMNFmt = 2;
        makenfmt(ctx, &ctx->PMNFmt,ctx->PMNFmtS,sizeof(ctx->PMNFmtS),1,ctx->SEPC);
    }
    if (ctx->MxIter < 0 || ctx->MxItFlg == 0)
        ctx->MxIter = 50;

    if (ctx->PMNV < 3) {
        printf1(ctx, "Error: need at least two variables on right-hand side.\n");
        goto CONJFin;
    }
    ctx->CJDIM = ctx->PMNV - 1;               /* number of dimensions */

    printf1(ctx, "\nRank order variable: %s\n",ctx->VName[ctx->PMVIdx[0]]);

    if (conj_r(ctx))                   /* check rank order */
        goto CONJFin;

    if (conj_tab(ctx))                 /* check dimensions */
        goto CONJFin;

    /* check for ties */

    if (alloc_acm(ctx, ctx->NOC + 1))
        goto CONJFin;

    printf1(ctx, "Mean of rank order variable: %g\n",ctx->CJRMean);
    printf1(ctx, "Stand. dev. of rank order variable: %g\n\n",ctx->CJRSD);

    if (ctx->PMALG == 1) {
        printf1(ctx, "Algorithm: alternating least squares.\n\n");

        if (conj_dmat(ctx, 0))                   /* create design matrix */
            goto CONJFin;

        ctx->NTies = mreg_ties(ctx, ctx->NOC,ctx->CJR,ctx->AcM);     /* AcM reserved to record ties */
        printf1(ctx, "Number of tie blocks: %d ",ctx->NTies);
        if (ctx->NTies == ctx->NOC)
            printf1(ctx, "(no ties).");
        printf1(ctx, "\n\n");

        if (conj_it(ctx))
            goto CONJFin;
    }
    else if (ctx->PMALG == 2) {
        printf1(ctx, "Algorithm: linear programming.\n\n");

        if (conj_dmat(ctx, 1))               /* create design matrix */
            goto CONJFin;

        ctx->NTies = mreg_ties(ctx, ctx->NOC,ctx->CJR,ctx->AcM);     /* AcM reserved to record ties */
        printf1(ctx, "Number of tie blocks: %d ",ctx->NTies);
        if (ctx->NTies == ctx->NOC)
            printf1(ctx, "(no ties).");
        printf1(ctx, "\n\n");


        if (conj_lp(ctx))                  /* try to solve the problem */
            goto CONJFin;
    }
    err = 0;

CONJFin:
    conj_free(ctx);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  conj_free. Free allocated memory.                                       */

void conj_free(TDAContext *ctx)
{
    register int i;

    if (ctx->CJNCatA > 0) {
        free((char *)ctx->CJNCat);
        memrq(ctx, -ctx->CJNCatA,sizeof(int));
        ctx->CJNCatA = 0;
    }
    if (ctx->CJCatA > 0) {
        for (i = 0; i < ctx->CJDIM; ++i) {
            if (ctx->CJCatAI[i] > 0) {
                free((char *)ctx->CJCat[i]);
                memrq(ctx, -ctx->CJCatAI[i],sizeof(int));
            }
        }
        free((char *)ctx->CJCat);
        memrq(ctx, -ctx->CJCatA,sizeof(int *));
        ctx->CJCatA = 0;
    }
    if (ctx->CJCatAIA > 0) {
        free((char *)ctx->CJCatAI);
        memrq(ctx, -ctx->CJCatAIA,sizeof(int));
        ctx->CJCatAIA = 0;
    }
    if (ctx->CJDMatA > 0) {
        free((char *)ctx->CJDMat);
        memrq(ctx, -ctx->CJDMatA,sizeof(char));
        ctx->CJDMatA = 0;
    }
    if (ctx->CJRA > 0) {
        free((char *)ctx->CJR);
        memrq(ctx, -ctx->CJRA,sizeof(double));
        ctx->CJRA = 0;
    }
    if (ctx->CJRYA > 0) {
        free((char *)ctx->CJRY);
        memrq(ctx, -ctx->CJRYA,sizeof(double));
        ctx->CJRYA = 0;
    }
    if (ctx->CJEYA > 0) {
        free((char *)ctx->CJEY);
        memrq(ctx, -ctx->CJEYA,sizeof(double));
        ctx->CJEYA = 0;
    }
}

/* ------------------------------------------------------------------------ */
/*  conj_r      check rank order and calculate mean and std. dev.           */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int conj_r(TDAContext *ctx)
{
    register int i;
    double tmp;

    if (!(ctx->CJR = (double *) calloc((size_t)(ctx->NOC + 1),sizeof(double)))) {
        p_err(ctx, -2,1);
        return(-1);
    }
    ctx->CJRA = ctx->NOC + 1;
    memrq(ctx, ctx->CJRA,sizeof(double));

    if (!(ctx->CJRY = (double *) calloc((size_t)(ctx->NOC + 1),sizeof(double)))) {
        p_err(ctx, -2,1);
        return(-1);
    }
    ctx->CJRYA = ctx->NOC + 1;
    memrq(ctx, ctx->CJRYA,sizeof(double));

    if (!(ctx->CJEY = (double *) calloc((size_t)(ctx->NOC + 1),sizeof(double)))) {
        p_err(ctx, -2,1);
        return(-1);
    }
    ctx->CJEYA = ctx->NOC + 1;
    memrq(ctx, ctx->CJEYA,sizeof(double));

    for (i = 0; i < ctx->NOC; ++i)
        ctx->CJR[i] = get_data(ctx, ctx->PMVIdx[0],i);

    ctx->CJRMean = 0.0;
    for (i = 0; i < ctx->NOC; ++i)
        ctx->CJRMean += get_data(ctx, ctx->PMVIdx[0],i);
    ctx->CJRMean /= (double)ctx->NOC;

    ctx->CJRSD = 0.0;
    for (i = 0; i < ctx->NOC; ++i) {
        tmp = get_data(ctx, ctx->PMVIdx[0],i) - ctx->CJRMean;
        ctx->CJRSD += tmp * tmp;
    }
    ctx->CJRSD /= (double)ctx->NOC;
    if (ctx->CJRSD > 0.0)
        ctx->CJRSD = sqrt(ctx->CJRSD);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  conj_tab    Create frequency distritution.                              */
/*              Return: 0 if OK, -1 if error.                               */

int conj_tab(TDAContext *ctx)
{
    register int i,j,k;
    int err,n;
    double tmp;

    err = -1;

    /*  Calculation of categories of dimensions of the table:
        CJNCat[i] = number of categories in i.th dimension,
        CJCat[i][j] (j = 1,2,3,...,CJNCat[i]) is the value of the j.th
        categorie in the i.th dimension. */

    if (!(ctx->CJNCat = (int *) calloc((size_t)(ctx->CJDIM),sizeof(int)))) {
        p_err(ctx, -2,1);
        goto CONJTFin;
    }
    ctx->CJNCatA = ctx->CJDIM;
    memrq(ctx, ctx->CJNCatA,sizeof(int));

    if (!(ctx->CJCatAI = (int *) calloc((size_t)(ctx->CJDIM),sizeof(int)))) {
        p_err(ctx, -2,1);
        goto CONJTFin;
    }
    ctx->CJCatAIA = ctx->CJDIM;
    memrq(ctx, ctx->CJCatAIA,sizeof(int));

    if (!(ctx->CJCat = (int **) calloc((size_t)(ctx->CJDIM),sizeof(int *)))) {
        p_err(ctx, -2,1);
        goto CONJTFin;
    }
    ctx->CJCatA = ctx->CJDIM;
    memrq(ctx, ctx->CJCatA,sizeof(int *));

    if (alloc_acr(ctx, ctx->NOC + 1))             /* data */
        goto CONJTFin;
    if (alloc_acs(ctx, ctx->NOC + 1))             /* values of categories */
        goto CONJTFin;
    if (alloc_aci(ctx, ctx->NOC + 1))             /* frequencies */
        goto CONJTFin;
    if (alloc_acj(ctx, ctx->NOC + 1))             /* bf pointer */
        goto CONJTFin;

    printf1(ctx, "Dimension  Variable ");
    prnchar(ctx, ' ',ctx->VNameLen - 8,0);
    printf1(ctx, "  Categories\n");

    ctx->CJNNCat = 0;
    for (i = 1; i <= ctx->CJDIM; ++i) {

        for (j = 0; j < ctx->NOC; ++j) {
            ctx->AcR[j + 1] = (int)get_data(ctx, ctx->PMVIdx[i],j);
        }
        n = cfreq(ctx, 1,ctx->NOC,ctx->AcR,ctx->NOC,ctx->AcS,ctx->AcI,ctx->AcJ,0,&tmp,&tmp);

        if (n < 1) {
            printf1(ctx, "Error: can't sort dimension %d.\n",i);
            goto CONJTFin;
        }
        ctx->CJNCat[i - 1] = n;
        ctx->CJNNCat += n;

        if (!(ctx->CJCat[i - 1] = (int *) calloc((size_t)(n + 1),sizeof(int)))) {
            p_err(ctx, -2,1);
            goto CONJTFin;
        }
        ctx->CJCatAI[i - 1] = n + 1;
        memrq(ctx, n + 1,sizeof(int));

        printf1(ctx, "%6d      %s ",i,ctx->VName[ctx->PMVIdx[i]]);
        prnchar(ctx, ' ',ctx->VNameLen - (int)strlen(ctx->VName[ctx->PMVIdx[i]]),0);
        printf1(ctx, "  %3d : ",n);
        for (j = 1; j <= n; ++j) {
            k = ctx->AcS[ctx->AcJ[j]];
            ctx->CJCat[i - 1][j] = k;
            rt_printf1_i(ctx, ctx->PMNFmtS,k);
        }
        newline(ctx);
    }
    newline(ctx);

    err = 0;

CONJTFin:
    alloc_acr(ctx, 0);
    alloc_acs(ctx, 0);
    alloc_aci(ctx, 0);
    alloc_acj(ctx, 0);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  conj_dmat   Make design matrix in CJDMat. Sort according to ranks.      */
/*              AcN contains categories corresp to design matrix.           */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int conj_dmat(TDAContext *ctx, int opt)
{
    register int i,j,k;
    int err,n,s,scat,first;

    err = -1;

    if (alloc_aci(ctx, ctx->NOC + 1))             /* temp. sort pointer */
        goto CJDMFin;

    if (alloc_acj(ctx, ctx->NOC + 1))             /* inverse pointer */
        goto CJDMFin;

    if (sortdp(ctx, ctx->NOC,ctx->CJR,ctx->AcI))
        goto CJDMFin;

    for (i = 0; i < ctx->NOC; ++i) {
        ctx->CJRY[i] = ctx->CJR[ctx->AcI[i]];
        ctx->AcJ[ctx->AcI[i]] = i;
    }
    printf1(ctx, "Indicator variables for design matrix:\n\n");

    printf1(ctx, "Idx  Dimension  Variable  ");
    prnchar(ctx, ' ',ctx->VNameLen - 8,0);
    printf1(ctx, "  Category\n");

    ctx->CJDMCol = 0;            /* number of colums in design matrix */

    for (j = 0; j < ctx->CJDIM; ++j) {
        for (k = 1; k <= ctx->CJNCat[j]; ++k) {

            if (k == 1 && opt == 0)
                continue;

            scat = ctx->CJCat[j][k];

            for (i = 0; i < ctx->NOC; ++i) {
                s = (int)get_data(ctx, ctx->PMVIdx[j + 1],i);
                if (s == scat) {
                    ctx->CJDMCol++;
                    printf1(ctx, "%3d  %6d      %s ",ctx->CJDMCol,j + 1,ctx->VName[ctx->PMVIdx[j + 1]]);
                    prnchar(ctx, ' ',ctx->VNameLen - (int)strlen(ctx->VName[ctx->PMVIdx[j + 1]]),0);
                    printf1(ctx, "%7d\n",s);
                    break;
                }
            }
        }
    }
    newline(ctx);

    if (opt == 0)
        ctx->CJDMCol++;      /* include intercept */

    if (!(ctx->CJDMat = (char *) calloc((size_t)(ctx->NOC) * (size_t)(ctx->CJDMCol) + 1,sizeof(char)))) {
        p_err(ctx, -2,1);
        goto CJDMFin;
    }
    ctx->CJDMatA = ctx->NOC * ctx->CJDMCol + 1;
    memrq(ctx, ctx->CJDMatA,sizeof(char));

    /* create values of design matrix in CJDMat. Also create:
       AcR : reference to dimensions
       AcS : reference to categories */

    if (alloc_acr(ctx, ctx->CJDMCol + 1))
        goto CJDMFin;

    if (alloc_acs(ctx, ctx->CJDMCol + 1))
        goto CJDMFin;

    if (alloc_acn(ctx, ctx->NOC * ctx->CJDIM + 1))     /* categories */
        goto CJDMFin;

    n = 0;
    if (opt == 0) {
        n = 1;
        for (i = 0; i < ctx->NOC; ++i)
            ctx->CJDMat[i * ctx->CJDMCol + 1] = 1;
    }
    for (j = 0; j < ctx->CJDIM; ++j) {
        for (k = 1; k <= ctx->CJNCat[j]; ++k) {

            scat = ctx->CJCat[j][k];
            first = 1;
            for (i = 0; i < ctx->NOC; ++i) {
                s = (int)get_data(ctx, ctx->PMVIdx[j + 1],i);
                if (s == scat) {

                    ctx->AcN[ctx->AcJ[i] * ctx->CJDIM + j + 1] = s;

                    if (k == 1 && opt == 0)
                        continue;

                    if (first) {
                        n++;
                        ctx->AcR[n] = j + 1;
                        ctx->AcS[n] = scat;
                        first = 0;
                    }
                    ctx->CJDMat[ctx->AcJ[i] * ctx->CJDMCol + n] = 1;
                }
            }
        }
    }
    for (i = 0; i < ctx->NOC; ++i)
        ctx->CJR[i] = ctx->CJRY[i];

    if (ctx->PMTabFDef) {            /* print design matrix */

        for (i = 0; i < ctx->NOC; ++i) {
            for (j = 1; j <= ctx->CJDMCol; ++j)
                rt_fprintf_i(ctx, ctx->PMTabFd,ctx->PMNFmtS,(int)ctx->CJDMat[i * ctx->CJDMCol + j]);
            fprintf(ctx->PMTabFd,"  %g\n",ctx->CJR[i]);
        }
        printf1(ctx, "Design matrix written to: %s\n\n",ctx->PMTabFName);
    }
    err = 0;

CJDMFin:
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  conj_it.   Iterative procedure, alternating least squares.              */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int conj_it(TDAContext *ctx)
{
    register int i,j,k;
    int err,iter,ranka,nc,nl,r;
    double tmp,rnorml,sa = 0.0,sb = 0.0,mz,sz;

    err = -1;

    printf1(ctx, "Handling ties option: %d\n",ctx->PMOPT);
    printf1(ctx, "Max number of iterations: %d\n",ctx->MxIter);
    printf1(ctx, "Tolerance for stress: %g\n",ctx->TOLF);
    printf1(ctx, "Tolerance for parameter change: %g\n",ctx->TOLSP);
    printf1(ctx, "Number of parameters: %d\n\n",ctx->CJDMCol);

    if (alloc_acx(ctx, ctx->CJDMCol + 1))             /* coefficients */
        goto CJITFin;

    /* ##  try to get starting values */

    if (get_dsv(ctx, ctx->CJDMCol,ctx->AcX,0,ctx->AcX,ctx->AcX,1))
        goto CJITFin;


#if 0   /* leftover debug output, disabled -- see README, mreg/conj */
    for (i = 1; i <= ctx->CJDMCol; ++i)
        tda_out("i=%d AcX=%f\n",i,ctx->AcX[i]);
#endif  /* leftover debug output */


    nc = ctx->CJDMCol + 1;

    if (alloc_acy(ctx, ctx->CJDMCol + 1))             /* previous coefficients */
        goto CJITFin;

    if (alloc_acw(ctx, ctx->NOC * nc + 1))             /* lsei matrix */
        goto CJITFin;

    if (ctx->PMOPT != 1) {
        if (alloc_actmp(ctx, ctx->NOC + 1))
            goto CJITFin;

        if (ctx->PMOPT == 2) {
            if (alloc_aci(ctx, ctx->NOC + 1))
                goto CJITFin;
            if (alloc_acj(ctx, ctx->NOC + 1))
                goto CJITFin;
        }
    }


    if (ctx->SILENTFlg < 2) {
        fflushe(ctx);
        printfe(ctx, "\n  Iter    Function Value        Par Change    Norm of Residuals\n");
        fflushe(ctx);
    }

    for (iter = 1; iter <= ctx->MxIter; ++iter) {

        if (iter == 1 && ctx->DSVFlg)
            goto CJITNxt;

        /* create lsei matrix */

        i = k = 0;
        while (i < ctx->NOC) {
            for (j = 1; j <= ctx->CJDMCol; ++j)
                ctx->AcW[i * nc + j] = (double)ctx->CJDMat[k * ctx->CJDMCol + j];
            ctx->AcW[i * nc + nc] = ctx->CJRY[k];
            i++;
            k++;
        }
        /***
        printf1(ctx, "W\n");
        for (i = 0; i < NOC; ++i) {
            for (j = 1; j <= nc; ++j)
                printf1(ctx, "%g ",AcW[i * nc + j]);
            newline(ctx);
        }
        ***/

        err = lsei(ctx, ctx->AcW,0,ctx->NOC,0,ctx->CJDMCol,ctx->AcX,0,&tmp,&rnorml,&ranka,&nl);
        if (err) {
            if (err == -4)
                p_err(ctx, -2,1);
            else
                printf1(ctx, "LSEI error return: %d\n",err);
            goto CJITFin;
        }

        /* standardize estimated dep. variable in AcZ and AcW */

CJITNxt:

#if 0   /* leftover debug output, disabled -- see README, mreg/conj */
        for (i = 1; i <= ctx->CJDMCol; ++i)
            tda_out("vorher i=%d AcX=%f\n",i,ctx->AcX[i]);
#endif  /* leftover debug output */


        for (i = 0; i < ctx->NOC; ++i) {
            tmp = 0.0;
            for (j = 1; j <= ctx->CJDMCol; ++j) {
                if (ctx->CJDMat[i * ctx->CJDMCol + j])
                    tmp += ctx->AcX[j];
            }
            ctx->CJRY[i] = tmp;
        }
        mreg_ms(ctx, ctx->NOC,ctx->CJRY,&mz,&sz);


#if 0   /* leftover debug output, disabled -- see README, mreg/conj */
        tda_out("Mean of z = %f\n",mz);
        tda_out("STD  of z = %f\n",sz);
#endif  /* leftover debug output */


        if (sz > 0.0) {
            for (j = 1; j <= ctx->CJDMCol; ++j)
                ctx->AcX[j] *= ctx->CJRSD / sz;
        }
        mz = 0.0;
        for (i = 0; i < ctx->NOC; ++i) {
            tmp = 0.0;
            for (j = 1; j <= ctx->CJDMCol; ++j) {
                if (ctx->CJDMat[i * ctx->CJDMCol + j])
                    tmp += ctx->AcX[j];
            }
            ctx->CJRY[i] = tmp;
            mz += tmp;
        }
        mz /= (double)ctx->NOC;

#if 0   /* leftover debug output, disabled -- see README, mreg/conj */
        tda_out("New Mean of CJRY = %f\n",mz);
#endif  /* leftover debug output */

        ctx->AcX[1] = ctx->AcX[1] + ctx->CJRMean - mz;       /* adjust intercept */

        for (i = 0; i < ctx->NOC; ++i) {
            tmp = 0.0;
            for (j = 1; j <= ctx->CJDMCol; ++j) {
                if (ctx->CJDMat[i * ctx->CJDMCol + j])
                    tmp += ctx->AcX[j];
            }
            ctx->CJRY[i] = tmp;
        }
        mreg_ms(ctx, ctx->NOC,ctx->CJRY,&mz,&sz);


#if 0   /* leftover debug output, disabled -- see README, mreg/conj */
        tda_out("Mean of CJRY = %f\n",mz);
        tda_out("STD  of CJRY = %f\n",sz);
#endif  /* leftover debug output */

        for (i = 0; i < ctx->NOC; ++i)
            ctx->CJEY[i] = ctx->CJRY[i];

        /* calculate projection, depending on PMOPT */

        if (ctx->PMOPT == 1)
            mreg_proj(ctx, ctx->NOC,ctx->CJRY);

        else if (ctx->PMOPT == 2)
            mreg_proj1(ctx, ctx->NOC,ctx->NTies,ctx->AcM,ctx->AcI,ctx->AcJ,ctx->CJRY,ctx->AcTmp);

        else if (ctx->PMOPT == 3)
            mreg_proj2(ctx, ctx->NOC,ctx->NTies,ctx->AcM,ctx->CJRY,ctx->AcTmp);


#if 0   /* leftover debug output, disabled -- see README, mreg/conj */
        tda_out("\nNew AcZ/CJRY after projection\n");
        for (i = 0; i < ctx->NOC; ++i)
            tda_out("i=%3d cjry=%f cjey=%f\n",i,ctx->CJRY[i],ctx->CJEY[i]);
#endif  /* leftover debug output */


        /* calculate stress */

        tmp = 0.0;
        for (i = 0; i < ctx->NOC; ++i) {
            tmp += ctx->CJEY[i];
        }
        tmp /= (double)ctx->NOC;

        sa = sb = 0.0;
        for (i = 0; i < ctx->NOC; ++i) {
            sa += (ctx->CJRY[i] - ctx->CJEY[i]) * (ctx->CJRY[i] - ctx->CJEY[i]);
            sb += (ctx->CJEY[i] - tmp) * (ctx->CJEY[i] - tmp);
        }
        if (sb > 0.0)
            sa /= sb;
        if (sa > 0.0)
            sa = sqrt(sa);

        /* calculate par change */

        sb = 0.0;
        for (j = 1; j <= ctx->CJDMCol; ++j) {
            tmp = fabs(ctx->AcX[j] - ctx->AcY[j]) / dmax(ctx, fabs(ctx->AcX[j]),1.0);
            sb = dmax(ctx, sb,tmp);
            ctx->AcY[j] = ctx->AcX[j];
        }
        if (ctx->SILENTFlg < 2) {
            printfe(ctx, "  %3d  %20.13e %17.10e %17.10e\n",iter,sa,sb,rnorml);
            fflushe(ctx);
        }
        if (sa <= ctx->TOLF || sb <= ctx->TOLSP || iter >= ctx->MxIter)
            break;
    }
    printf1(ctx, "Performed %d iterations.\n",iter);
    printf1(ctx, "Final stress: %g\n",sa);
    printf1(ctx, "Final parameter change: %g\n",sb);
    printf1(ctx, "Final norm of least squares residuals: %g\n",rnorml);
    printf1(ctx, "Rank of design matrix: %d\n\n",ranka);

    conj_prn(ctx); /* print coefficients */

    r = conj_eval(ctx);
    printf1(ctx, "Type of equivalence: %d\n",r);

    if (ctx->PMResFDef)
        conj_res(ctx);

    err = 0;

CJITFin:
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  conj_eval()   Compare CJR and CJEY. Return                              */
/*                1  if strict monotone equivalent                          */
/*                2  if half-weak monotone equivalent                       */
/*                0  otherwise                                              */
/*                                                                          */
/*  It is assumed that CJR is sorted in increasing order.                   */

int conj_eval(TDAContext *ctx)
{
    register int i,j;
    int r;
    double tol = sqrt(ctx->EPSI1);

    r = 1;
    for (i = 1; i < ctx->NOC; ++i) {
        if (fabs(ctx->CJR[i] - ctx->CJR[i - 1]) <= tol) {
            if (fabs(ctx->CJEY[i] - ctx->CJEY[i - 1]) > tol) {
                r = 2;
                goto CJENxt;
            }
        }
        else if (ctx->CJEY[i] < ctx->CJEY[i - 1] + tol) {
            r = 2;
            goto CJENxt;
        }
    }

CJENxt:
    if (r == 1)
        return(1);

    for (i = 1; i < ctx->NOC; ++i) {
        if (ctx->CJR[i] > ctx->CJR[i - 1] + tol) {
            if (ctx->CJEY[i] < ctx->CJEY[i - 1] - tol)
                return(0);
        }
    }
    for (i = 1; i < ctx->NOC; ++i) {
        for (j = 0; j < i; ++j) {
            if (ctx->CJEY[j] < ctx->CJEY[i] - tol) {
                if (ctx->CJR[j] > ctx->CJR[i] + tol)
                    return(0);
            }
            else if (ctx->CJEY[j] > ctx->CJEY[i] + tol) {
                if (ctx->CJR[j] < ctx->CJR[i] - tol)
                    return(0);
            }
        }
    }
    return(2);
}

/* ------------------------------------------------------------------------ */
/*  conj_prn.  Print coefficients.                                          */
/*                                                                          */

void conj_prn(TDAContext *ctx)
{
    register int j,k;

    printf1(ctx, "Idx  Dimension  Variable  ");
    prnchar(ctx, ' ',ctx->VNameLen - 8,0);
    printf1(ctx, "  Category  Coefficient\n");

    for (j = 1; j <= ctx->CJDMCol; ++j) {
        k = ctx->AcR[j];
        if (j == 1) {
            printf1(ctx, "%3d       -      Const ",j - 1);
            prnchar(ctx, ' ',ctx->VNameLen - 5,0);
            printf1(ctx, "       -    ");
        }
        else {
            printf1(ctx, "%3d  %6d      %s ",j - 1,k,ctx->VName[ctx->PMVIdx[k]]);
            prnchar(ctx, ' ',ctx->VNameLen - (int)strlen(ctx->VName[ctx->PMVIdx[k]]),0);
            printf1(ctx, "  %6d    ",ctx->AcS[j]);
        }
        rt_printf1_d(ctx, ctx->PMTFmtS,ctx->AcX[j]);
        newline(ctx);
    }
    newline(ctx);
}

/* ------------------------------------------------------------------------ */
/*  conj_res.  Print estimation results                                     */
/*                                                                          */

void conj_res(TDAContext *ctx)
{
    register int i,j;

    for (i = 0; i < ctx->NOC; ++i) {
        fprintf(ctx->PMResFd,"%4d ",i + 1);
        for (j = 1; j <= ctx->CJDIM; ++j) {
            rt_fprintf_i(ctx, ctx->PMResFd,ctx->PMNFmtS,ctx->AcN[i * ctx->CJDIM + j]);
        }
        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMFmtS,ctx->CJR[i]);
        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMFmtS,ctx->CJRY[i]);
        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMFmtS,ctx->CJEY[i]);
        fprintf(ctx->PMResFd,"\n");
    }
    printf1(ctx, "Estimation results written to: %s\n\n",ctx->PMResFName);
}

/* ------------------------------------------------------------------------ */
/*  mreg()      Monotone regression.                                        */
/*                                                                          */
/*              mreg(                                                       */
/*                  opt=...,    1 or 2 or 3, def. 1                         */
/*                              1 = ties ignored                            */
/*                              2 = Kruskal's primary approach              */
/*                              3 = Kruskal's secondary approach            */
/*                  mxit=...,   max number of iterations, def. 50           */
/*                  tolf=...,   tolerance for stress, def. 1e-6             */
/*                  tolsp=...,  tolerance for par change, def. 1.e-6        */
/*                  ni=...,     1 if without intercept                      */
/*                  ppar=...,   print parameter to output file              */
/*                  tfmt=...,   print format for parameters, def. 10.4      */
/*                  df=...,     optional output file                        */
/*                  fmt=...,    print format for df option, def. 0.0        */
/*              ) = Y,X1,X2,...;                                            */
/*                                                                          */
/*              Note: only with cross-sectional data.                       */
/*              Return 0 if OK, otherwise -1.                               */

int mreg(TDAContext *ctx)
{
    register int i,j,k;
    int err,nv,nx1;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Monotone regression. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLF = 1.e-6;
    ctx->TOLSP = 1.e-6;

    if (parm(ctx, ctx->CmdBuf + 4,4,1))     /* get parameters */
        goto MRFin;

    if (ctx->PMNW != 1) {
        printf1(ctx, "Error: mreg command requires cross-section data (nw=1).\n");
        goto MRFin;
    }
    if (ctx->MxIter < 0 || ctx->MxItFlg == 0)
        ctx->MxIter = 50;

    newline(ctx);
    nv = check_nvar(ctx, 0);         /* check variables */
    if (nv == 0)
        goto MRFin;

    prn_nwvar(ctx, 0);               /* print variables */
    nx1 = nv - 1;
    nx1++;
    if (nx1 == 0)
        goto MRFin;

    /* allocate memory */

    if (alloc_acx(ctx, nx1 * ctx->NOC + 1))       /* x matrix */
        goto MRFin;
    if (alloc_acy(ctx, ctx->NOC + 1))             /* y vector */
        goto MRFin;
    if (alloc_ack(ctx, ctx->NOC + 1))             /* sort ptr for y */
        goto MRFin;


    /* get data for dependent variable and sort */

    for (i = 0; i < ctx->NOC; ++i)
        ctx->AcY[i] = get_data(ctx, ctx->PMVIdx[0],i);

    if (sortdp(ctx, ctx->NOC,ctx->AcY,ctx->AcK)) {
        p_err(ctx, -2,1);
        goto MRFin;
    }

    /* get data in sorted order */

    for (i = 0; i < ctx->NOC; ++i) {
        ctx->AcY[i] = get_data(ctx, ctx->PMVIdx[0],ctx->AcK[i]);
        k = 1;
        ctx->AcX[i * nx1 + 1] = 1.0;
        k = 2;

        for (j = 1; j < ctx->PMNV; ++j) {
            ctx->AcX[i * nx1 + k] = get_data(ctx, ctx->PMVIdx[j],ctx->AcK[i]);
            k++;
        }
    }

    /* get mean and std.dev. of AcY in CJRMean and CJRSD */

    mreg_ms(ctx, ctx->NOC,ctx->AcY,&ctx->CJRMean,&ctx->CJRSD);

    printf1(ctx, "Mean of dep. variable: %g\n",ctx->CJRMean);
    printf1(ctx, "Stand. dev. of dep. variable: %g\n",ctx->CJRSD);

    /* check for ties */

    if (alloc_acn(ctx, ctx->NOC + 1))
        goto MRFin;

    ctx->NTies = mreg_ties(ctx, ctx->NOC,ctx->AcY,ctx->AcN);

    printf1(ctx, "Number of tie blocks: %d ",ctx->NTies);
    if (ctx->NTies == ctx->NOC)
        printf1(ctx, "(no ties).");
    printf1(ctx, "\n\n");

    if (mreg_it(ctx, nx1))
        goto MRFin;
    err = 0;

MRFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  mreg_it.   Iterative procedure. AcX data matrix, AcY dep. variable.     */
/*             Both sorted according to AcY.                                */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int mreg_it(TDAContext *ctx, int nx1)
{
    register int i,j;
    int err,iter,ranka,n,nc;
    double tmp,rnorml,sa = 0.0,sb = 0.0;

    err = -1;

    printf1(ctx, "Handling ties option: %d\n",ctx->PMOPT);
    printf1(ctx, "Max number of iterations: %d\n",ctx->MxIter);
    printf1(ctx, "Tolerance for stress: %g\n",ctx->TOLF);
    printf1(ctx, "Tolerance for parameter change: %g\n",ctx->TOLSP);
    printf1(ctx, "Number of parameters: %d\n\n",nx1);

    nc = nx1 + 1;

    if (alloc_acz(ctx, ctx->NOC + 1))                 /* updates of AcY */
        goto MRITFin;

    if (alloc_acu(ctx, nx1 + 1))                 /* coefficients */
        goto MRITFin;

    if (alloc_acv(ctx, nx1 + 1))                 /* previous coefficients */
        goto MRITFin;

    if (alloc_acw(ctx, ctx->NOC * nc + 1))            /* lsei matrix */
        goto MRITFin;

    if (ctx->PMOPT != 1) {
        if (alloc_actmp(ctx, ctx->NOC + 1))
            goto MRITFin;

        if (ctx->PMOPT == 2) {
            if (alloc_aci(ctx, ctx->NOC + 1))
                goto MRITFin;
            if (alloc_acj(ctx, ctx->NOC + 1))
                goto MRITFin;
        }
    }

    if (ctx->SILENTFlg < 2) {
        fflushe(ctx);
        printfe(ctx, "\n  Iter         Stress           Par Change    Norm of Residuals\n");
        fflushe(ctx);
    }
    printf1(ctx, "  Iter         Stress           Par Change    Norm of Residuals\n");

    for (i = 0; i < ctx->NOC; ++i)
        ctx->AcZ[i] = ctx->AcY[i];

    for (iter = 1; iter <= ctx->MxIter; ++iter) {

        /* create lsei matrix */

#if 0   /* leftover debug output, disabled -- see README, mreg/conj */
newline(ctx);
for (i = 0; i < ctx->NOC; ++i)
tda_out("VREG   i=%3d  AcZ=%f\n",i,ctx->AcZ[i]);
#endif  /* leftover debug output */


        for (i = 0; i < ctx->NOC; ++i) {
            for (j = 1; j <= nx1; ++j)
                ctx->AcW[i * nc + j] = ctx->AcX[i * nx1 + j];
            ctx->AcW[i * nc + nc] = ctx->AcZ[i];
        }
        err = lsei(ctx, ctx->AcW,0,ctx->NOC,0,nx1,ctx->AcU,0,&tmp,&rnorml,&ranka,&n);
        if (err) {
            if (err == -4)
                p_err(ctx, -2,1);
            else
                printf1(ctx, "LSEI error return: %d\n",err);
            goto MRITFin;
        }

        /* standardize estimated dep. variable in AcZ and AcW */

        mreg_std(ctx, ctx->NOC,nx1,ctx->AcX,ctx->AcU,ctx->AcZ,ctx->CJRMean,ctx->CJRSD);

        for (i = 0; i < ctx->NOC; ++i)
            ctx->AcW[i] = ctx->AcZ[i];


        /* calculate projection, depending on PMOPT */

        if (ctx->PMOPT == 1)
            mreg_proj(ctx, ctx->NOC,ctx->AcZ);

        else if (ctx->PMOPT == 2)
            mreg_proj1(ctx, ctx->NOC,ctx->NTies,ctx->AcN,ctx->AcI,ctx->AcJ,ctx->AcZ,ctx->AcTmp);

        else if (ctx->PMOPT == 3)
            mreg_proj2(ctx, ctx->NOC,ctx->NTies,ctx->AcN,ctx->AcZ,ctx->AcTmp);



#if 0   /* leftover debug output, disabled -- see README, mreg/conj */
        tda_out("\nNew AcZ after projection\n");
        for (i = 0; i < ctx->NOC; ++i)
            tda_out("i=%3d acz=%f acw=%f\n",i,ctx->AcZ[i],ctx->AcW[i]);
#endif  /* leftover debug output */

        /* calculate stress */

        tmp = 0.0;
        for (i = 0; i < ctx->NOC; ++i) {
            tmp += ctx->AcW[i];
        }
        tmp /= (double)ctx->NOC;

        sa = sb = 0.0;
        for (i = 0; i < ctx->NOC; ++i) {
#if 0   /* leftover debug output, disabled -- see README, mreg/conj */
tda_out("Stress i=%d AcZ=%f AcW=%f\n",i,ctx->AcZ[i],ctx->AcW[i]);
#endif  /* leftover debug output */

            sa += (ctx->AcZ[i] - ctx->AcW[i]) * (ctx->AcZ[i] - ctx->AcW[i]);
            sb += (ctx->AcW[i] - tmp) * (ctx->AcW[i] - tmp);
        }
        if (sb > 0.0)
            sa /= sb;
        if (sa > 0.0)
            sa = sqrt(sa);
#if 0   /* leftover debug output, disabled -- see README, mreg/conj */
tda_out("Stress sa=%f\n\n",sa);
#endif  /* leftover debug output */

        /* calculate par change */

        sb = 0.0;
        for (j = 1; j <= nx1; ++j) {
            tmp = fabs(ctx->AcU[j] - ctx->AcV[j]) / dmax(ctx, fabs(ctx->AcU[j]),1.0);
            sb = dmax(ctx, sb,tmp);
            ctx->AcV[j] = ctx->AcU[j];
        }
        if (ctx->SILENTFlg < 2) {
            printfe(ctx, "  %3d  %20.13e %17.10e %17.10e\n",iter,sa,sb,rnorml);
            fflushe(ctx);
        }
        printf1(ctx, "  %3d  %20.13e %17.10e %17.10e\n",iter,sa,sb,rnorml);

        if (sa <= ctx->TOLF || sb <= ctx->TOLSP || iter >= ctx->MxIter)
            break;

    }
    newline(ctx);
    printf1(ctx, "Performed %d iterations.\n",iter);
    printf1(ctx, "Final stress: %g\n",sa);
    printf1(ctx, "Final parameter change: %g\n",sb);
    printf1(ctx, "Final norm of least squares residuals: %g\n",rnorml);
    printf1(ctx, "Rank of design matrix: %d\n\n",ranka);

    prn1_coeff(ctx, nx1,ctx->AcU,ctx->AcV,0,0,ctx->PMVIdx,0);   /* print parameters */

    if (ctx->PMF1Def) {          /* write data etc to output file */
        for (i = 0; i < ctx->NOC; ++i) {
            fprintf(ctx->PMF1d,"%4d ",i + 1);
            for (j = 1; j <= nx1; ++j)
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[i * nx1 + j]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[i]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcW[i]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcZ[i]);
            fprintf(ctx->PMF1d,"\n");
        }
        printf1(ctx, "\nData and results written to: %s\n\n",ctx->PMF1dName);
    }
    err = 0;

MRITFin:
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  mreg_proj(n,x)    Return monotone projection of x in x.                 */
/*                    Returns always 0.                                     */

int mreg_proj(TDAContext *ctx, int n, double *x)
{
    register int i,j,ia,nn;
    double r,rr,ry;

    while (1) {
        ia = i = 0;
        rr = ry = x[i];
        nn  = 1;
        while (i < n) {
            while (++i < n) {
                r  = x[i];
                if (r <= rr + ctx->EPSI1) {
                    ry += r;
                    nn++;
                }
                else
                    break;
                rr = ry / (double)nn;

            }
            for (j = ia; j < i; ++j)
                x[j] = rr;

            ia = i;
            rr = ry = x[i];
            nn = 1;
        }
        nn = 0;
        for (i = 1; i < n; ++i) {
            if (x[i] < x[i - 1] - ctx->EPSI1) {
                nn = 1;
                break;
            }
        }
        if (nn == 0)
            break;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mreg_proj1(n,m,t,p,pp,x,y)                                              */
/*                                                                          */
/*                         Return projection of x in x. If Ties, then       */
/*                         Kruskal's primary approach. t is integer array   */
/*                         containing the tie information.                  */
/*                                                                          */
/*                         p and pp are integer arrays for sorting.         */
/*                         Return 0 if OK, -1 if insuff. memory.            */

int mreg_proj1(TDAContext *ctx, int n,int m,int *t,int *p,int *pp,double *x,double *y)
{
    register int i,j,k;

    for (i = 0; i < n; ++i)
        p[i] = i;

    i = 0;
    for (j = 0; j < m; ++j) {
        if (t[j] == 2 && x[i] > x[i + 1]) {
            p[i] = i + 1;
            p[i + 1] = i;
        }
        else if (t[j] > 2) {
            if (sortdp(ctx, t[j],x + i,pp))
                return(-1);
            for (k = 0; k < t[j]; ++k)
                p[i + k] = pp[k] + i;
        }
        i += t[j];
    }
    for (i = 0; i < n; ++i)
        y[i] = x[p[i]];

    mreg_proj(ctx, n,y);

    for (i = 0; i < n; ++i)
        x[p[i]] = y[i];

    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mreg_proj2(n,m,t,x,y)  Return projection of x in x. If Ties, then       */
/*                         Kruskal's secondary approach. t is integer array */
/*                         containing the tie information. m is number of   */
/*                         entries in t. y auxilliary.                      */
/*                                                                          */
/*                         Returns always 0.                                */

int mreg_proj2(TDAContext *ctx, int n,int m,int *t,double *x,double *y)
{
    (void)n;        /* unused: the signature is shared */
    register int i,j,k;
    double tmp;

    i = 0;
    for (j = 0; j < m; ++j) {
        tmp = 0.0;
        for (k = 0; k < t[j]; ++k)
            tmp += x[i + k];
        tmp /= (double)t[j];
        y[j] = tmp;
        i += t[j];
    }
    mreg_proj(ctx, m,y);

    i = 0;
    for (j = 0; j < m; ++j) {
        for (k = 0; k < t[j]; ++k)
            x[i + k] = y[j];
        i += t[j];
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mreg_ms(n,y,m,s)   get mean and std.dev. of x.                          */
/*                                                                          */

void mreg_ms(TDAContext *ctx, int n,double *y,double *m,double *s)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;
    double r;

    *m = 0.0;
    for (i = 0; i < n; ++i)
        *m += y[i];
    *m /= (double)n;

    *s = 0.0;
    for (i = 0; i < n; ++i) {
        r = y[i] - *m;
        *s += r * r;
    }
    *s /= (double)n;
    if (*s > 0.0)
        *s = sqrt(*s);
}

/* ------------------------------------------------------------------------ */
/*  mreg_std(n,nx,x,beta,z,m,s)                                             */
/*                                                                          */
/*  Standardize estimated dependend variable with mean m and std.dev. s.    */
/*  x is n,nx matrix, beta are current parameters.                          */
/*  Return new parameters in beta and stand. dep. variable in z.            */
/*                                                                          */

void mreg_std(TDAContext *ctx, int n,int nx,double *x,double *beta,double *z,double m,double s)
{
    register int i,j;
    double tmp,mz,sz;

#if 0   /* leftover debug output, disabled -- see README, mreg/conj */
    for (i = 1; i <= nx; ++i)
        tda_out("vorher i=%d beta=%f\n",i,beta[i]);
#endif  /* leftover debug output */


    for (i = 0; i < n; ++i) {
        tmp = 0.0;
        for (j = 1; j <= nx; ++j)
            tmp += x[i * nx + j] * beta[j];
        z[i] = tmp;
    }
    mreg_ms(ctx, n,z,&mz,&sz);

#if 0   /* leftover debug output, disabled -- see README, mreg/conj */
tda_out("Mean of z = %f\n",mz);
tda_out("STD  of z = %f\n",sz);
#endif  /* leftover debug output */

    if (sz > 0.0) {
        for (j = 1; j <= nx; ++j)
            beta[j] *= s / sz;
    }
    mz = 0.0;
    for (i = 0; i < n; ++i) {
        tmp = 0.0;
        for (j = 1; j <= nx; ++j)
            tmp += x[i * nx + j] * beta[j];
        z[i] = tmp;
        mz += tmp;
    }
    mz /= (double)n;

#if 0   /* leftover debug output, disabled -- see README, mreg/conj */
tda_out("New Mean of z = %f\n",mz);
#endif  /* leftover debug output */

    beta[1] = beta[1] + m - mz;

    for (i = 0; i < n; ++i) {
        tmp = 0.0;
        for (j = 1; j <= nx; ++j)
            tmp += x[i * nx + j] * beta[j];
        z[i] = tmp;
    }
    mreg_ms(ctx, n,z,&mz,&sz);

#if 0   /* leftover debug output, disabled -- see README, mreg/conj */
tda_out("Mean of z = %f\n",mz);
tda_out("STD  of z = %f\n",sz);
#endif  /* leftover debug output */

}

/* ------------------------------------------------------------------------ */
/*  mreg_ties(n,y,t)                                                        */
/*                                                                          */
/*  Check n-vector y for ties. Return number of tie blocks. Sizes of        */
/*  tie block in t[].                                                       */

int mreg_ties(TDAContext *ctx, int n,double *y,int *t)
{
    register int i;
    int m;

    m = 0;
    t[m] = 1;

    for (i = 1; i < n; ++i) {
        if (fabs(y[i] - y[i - 1]) <= ctx->EPSI1)
            t[m] += 1;
        else {
            m++;
            t[m] = 1;
        }
    }
    m++;

    /***
    tda_out("\nNumber of tie blocks m=%d\n",m);
    for (i = 0; i < m; ++i)
        tda_out("i=%3d t[i]=%d\n",i,t[i]);
    **/
    return(m);
}

/* -##--------------------------------------------------------------------- */
/*  nmca        Non-metric conjoint analysis.                               */
/*                                                                          */
/*              nmca(                                                       */
/*                  df=...,         output file                             */
/*                  nfmt = ...,     integer format, def. 2                  */
/*                  tfmt=...,       print format coefficients, 10.4         */
/*                  fmt=...,        format for pres, def. 8.4               */
/*              ) = R,X1,...,Xm;                                            */
/*                                                                          */
/*  R contains rank ordering. Floating point values.                        */
/*  X1,...,Xm are integer valued.                                           */
/*                                                                          */
/*  The command finds all partial orders for the m dimensions such that     */
/*  the order R is a consistent continuation. The partial orders are        */
/*  written into the standard output. If the df option is used, the         */
/*  permuted input date are written into the output file, separately for    */
/*  each partial order.                                                     */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int nmca(TDAContext *ctx)
{
    int err,r;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Non-metric conjoint analysis. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,4,1))     /* get parameters */
        goto NMCAFin;

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 8,4);

    if (ctx->PMNFmtF == 0) {
        ctx->PMNFmt = 2;
        makenfmt(ctx, &ctx->PMNFmt,ctx->PMNFmtS,sizeof(ctx->PMNFmtS),1,ctx->SEPC);
    }
    if (ctx->PMNV < 3) {
        printf1(ctx, "Error: need at least two variables on right-hand side.\n");
        goto NMCAFin;
    }
    ctx->CJDIM = ctx->PMNV - 1;               /* number of dimensions */

    printf1(ctx, "\nRank order variable: %s\n",ctx->VName[ctx->PMVIdx[0]]);

    if (conj_r(ctx))                   /* check rank order */
        goto NMCAFin;

    if (conj_tab(ctx))                 /* check dimensions */
        goto NMCAFin;

    r = nmca_com(ctx);
    if (r < 0)
        goto NMCAFin;

    printf1(ctx, "\nFound %d partial orders.\n",r);
    if (ctx->PMF1Def)
        printf1(ctx, "%d records written to: %s\n",r * ctx->NOC,ctx->PMF1dName);

    err = 0;

NMCAFin:
    conj_free(ctx);
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  nmca_com().   Combinatorial procedure.                                  */
/*                                                                          */
/*  Return -1 if error, otherwise the number of partial orders that have    */
/*  been found.                                                             */

int nmca_com(TDAContext *ctx)
{
    register int i,j,k,ii;
    int err,rn,r,is,iis,np;

    err = -1;
    np = 0;

    /*  CJNCat[i] = number of categories in i.th dimension,
        CJCat[i][j] (j = 1,2,3,...,CJNCat[i]) is the value of the j.th
        categorie in the i.th dimension. */


    if (alloc_acy(ctx, ctx->NOC + 1))             /* rank order variable  */
        goto CONJCFin;
    if (alloc_acr(ctx, ctx->CJDIM * ctx->NOC + 1))     /* indep. variables  */
        goto CONJCFin;
    if (alloc_acs(ctx, ctx->NOC + 1))             /* temp. sort pointer */
        goto CONJCFin;

    for (i = 0; i < ctx->NOC; ++i) {
        ctx->AcY[i] = get_data(ctx, ctx->PMVIdx[0],i);
        for (j = 0; j < ctx->CJDIM; ++j) {
            k = (int)get_data(ctx, ctx->PMVIdx[j + 1],i) - 1;
            if (k < 0 || k >= ctx->CJNCat[j]) {
                printf1(ctx, "Error in record %d category %d (dim %d, cat %d).\n",
                                i + 1,k,j + 1,ctx->CJNCat[i]);
                return(-1);
            }
            ctx->AcR[i * ctx->CJDIM + j] = k;
        }
    }
    if (sortdp(ctx, ctx->NOC,ctx->AcY,ctx->AcS))
        goto CONJCFin;

    /***
    for (i = 0; i < NOC; ++i) {
        is = AcS[i];
        tda_out("%f ",AcY[is]);
        for (j = 0; j < CJDIM; ++j)
            tda_out("%2d ",AcR[is * CJDIM + j]);
        newline(ctx);
    }
    ***/

    if (alloc_acn(ctx, ctx->CJNNCat + 1))     /* record permutations */
        goto CONJCFin;
    if (alloc_aci(ctx, ctx->CJNNCat + 1))
        goto CONJCFin;
    if (alloc_acj(ctx, ctx->CJNNCat + 1))
        goto CONJCFin;
    if (alloc_ack(ctx, ctx->CJDIM + 1))       /* offsets */
        goto CONJCFin;

    k = 0;
    for (i = 0; i < ctx->CJDIM; ++i) {
        ctx->AcK[i] = k;
        k += ctx->CJNCat[i];
    }

    /***
    tda_out("offsets: ");
    for (i = 0; i < CJDIM; ++i)
        tda_out("%d ",AcK[i]);
    newline(ctx);
    ***/

    /* generate first permutation */

    for (i = 0; i < ctx->CJDIM; ++i) {
        k = ctx->AcK[i];
        perm(ctx, ctx->CJNCat[i],ctx->AcN + k,ctx->AcI + k,ctx->AcJ + k,1);
        ctx->AcR[i] = 0;
    }

    while (1) {

        /*********************
        newline(ctx);
        for (i = 0; i < n; ++i)
            tda_out("%d ",AcN[i]);
        newline(ctx);
        *********************/

        r = 1;
        for (i = 0; i < ctx->NOC; ++i) {
            is = ctx->AcS[i];
            for (ii = i + 1; ii < ctx->NOC; ++ii) {
                iis = ctx->AcS[ii];

                /******************************
                tda_out("    compare: ");
                for (j = 0; j < CJDIM; ++j)
                    tda_out("%2d ",AcR[is * CJDIM + j]);
                tda_out(" and: ");
                for (j = 0; j < CJDIM; ++j)
                    tda_out("%2d ",AcR[iis * CJDIM + j]);
                newline(ctx);
                tda_out("NEW compare: ");
                for (j = 0; j < CJDIM; ++j) {
                    k = AcK[j];
                    tda_out("%2d ",AcN[k + AcR[is * CJDIM + j]]);

                }
                tda_out(" and: ");
                for (j = 0; j < CJDIM; ++j) {
                    k = AcK[j];
                    tda_out("%2d ",AcN[k + AcR[iis * CJDIM + j]]);
                }
                newline(ctx);
                *******************/

                rn = 0;
                for (j = 0; j < ctx->CJDIM; ++j) {
                    k = ctx->AcK[j];
                    if (ctx->AcN[k + ctx->AcR[is * ctx->CJDIM + j]] >=
                        ctx->AcN[k + ctx->AcR[iis * ctx->CJDIM + j]]) {
                        rn++;
                    }
                }
                if (rn == ctx->CJDIM) {
                    r = 0;
                    break;
                }
            }
            if (r == 0)
                break;
        }
        if (r) {        /* found valid permutation */
            np++;
            printf1(ctx, "Found: ");
            for (i = 0; i < ctx->CJDIM; ++i) {
                printf1(ctx, "(");
                k = ctx->AcK[i];
                ii = ctx->CJNCat[i];
                for (j = 0; j < ii; ++j) {
                    rt_printf1_i(ctx, ctx->PMNFmtS,ctx->AcN[k + j] + 1);
                    if (j == ii - 1)
                        printf1(ctx, ") ");
                    else
                        printf1(ctx, ",");
                }
            }
            newline(ctx);

            if (ctx->PMF1Def) {      /* write to output file */

                for (i = 0; i < ctx->NOC; ++i) {
                    is = ctx->AcS[i];
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[is]);
                    for (j = 0; j < ctx->CJDIM; ++j) {
                        k = ctx->AcK[j];
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->AcN[k + ctx->AcR[is * ctx->CJDIM + j]] + 1);
                    }
                    fprintf(ctx->PMF1d,"\n");
                }
            }
        }

        /* generate next permutation */

        for (i = ctx->CJDIM - 1; i >= 0; --i) {
            k = ctx->AcK[i];
            r = perm(ctx, ctx->CJNCat[i],ctx->AcN + k,ctx->AcI + k,ctx->AcJ + k,0);
            if (r == 1)
                break;

            if (i == 0)
                goto CONJCFin1;

            perm(ctx, ctx->CJNCat[i],ctx->AcN + k,ctx->AcI + k,ctx->AcJ + k,1);
        }
    }

CONJCFin1:
    err = np;

CONJCFin:
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  conj_lp()   Conjoint analysis with linear programming.                  */
/*                                                                          */
/*  Return 0 if OK, 1 if no solution, -1 if insufficient memory.            */

int conj_lp(TDAContext *ctx)
{
    register int i,j,k;
    int err,n,m,m1,ai,bi,p;
    double z,tol;

    err = -1;
    tol = 1.0e-8;

    n = ctx->NOC;
    m = ctx->CJDMCol;
    m1 = m + 1;

    if (alloc_actmp(ctx, n * m1 + 1))
        goto CONJLPFin;

    if (alloc_acx(ctx, m + 1))
        goto CONJLPFin;

    if (alloc_acy(ctx, n + 1))
        goto CONJLPFin;

    for (j = 1; j <= m; ++j)
        ctx->AcTmp[j] = -1.0;

    p = 0;
    k = 0;
    for (i = 1; i < n; ++i) {
        if (fabs(ctx->CJR[i] - ctx->CJR[i - 1]) >= ctx->EPSI1) {       /* not tied */
            k++;
            for (j = 1; j <= m; ++j) {
                ai = ctx->CJDMat[(i - 1) * m + j];
                bi = ctx->CJDMat[i * m + j];
                ctx->AcTmp[k * m1 + j] = (double)(ai - bi);
            }
            ctx->AcTmp[k * m1 + m1] = -1.0;
        }
    }
    for (i = 1; i < n; ++i) {
        if (fabs(ctx->CJR[i] - ctx->CJR[i - 1]) < ctx->EPSI1) {       /* tied */
            p++;
            k++;
            for (j = 1; j <= m; ++j) {
                ai = ctx->CJDMat[(i - 1) * m + j];
                bi = ctx->CJDMat[i * m + j];
                ctx->AcTmp[k * m1 + j] = (double)(ai - bi);
            }
            ctx->AcTmp[k * m1 + m1] = 0.0;
        }
    }
    /***********
    tda_out("actmp p=%d  \n",p );
    for (i = 0; i < n; ++i) {
        for (j = 1; j <= m1; ++j)
            tda_out("%4.1f ",AcTmp[i * m1 + j]);
        newline(ctx);
    }
    **********/

    err = lpf1(ctx, n - 1,m,p,ctx->AcTmp,ctx->AcX,ctx->AcY,&z,tol);
    if (err) {
        if (err == -5) {
            err = -1;
            goto CONJLPFin;
        }
        else {
            err = 1;
            goto CONJLPFin;
        }
    }

    /* create estimated rank order in CJRY and CJEY */

    for (i = 0; i < n; ++i) {
        z = 0.0;
        for (j = 1; j <= m; ++j) {
            if (ctx->CJDMat[i * m + j])
                z += ctx->AcX[j];
        }
        ctx->CJRY[i] = ctx->CJEY[i] = z;
    }

    /* check for valid solution (strict monotone equivalence) */

    if (conj_eval(ctx) != 1) {
        printf1(ctx, "Cannot find a solution.\n");
        err = 1;
        goto CONJLPFin;
    }

    printf1(ctx, "Idx  Dimension  Variable  ");
    prnchar(ctx, ' ',ctx->VNameLen - 8,0);
    printf1(ctx, "  Category    Solution\n");

    for (j = 1; j <= ctx->CJDMCol; ++j) {
        k = ctx->AcR[j];
        printf1(ctx, "%3d  %6d      %s ",j - 1,k,ctx->VName[ctx->PMVIdx[k]]);
        prnchar(ctx, ' ',ctx->VNameLen - (int)strlen(ctx->VName[ctx->PMVIdx[k]]),0);
        printf1(ctx, "  %6d    ",ctx->AcS[j]);
        rt_printf1_d(ctx, ctx->PMTFmtS,ctx->AcX[j]);
        newline(ctx);
    }
    printf1(ctx, "\nType of equivalence: 1\n");
    if (ctx->PMResFDef)
        conj_res(ctx);

CONJLPFin:
    if (err == 1)
        printf1(ctx, "Cannot find a solution.\n");
    else if (err == -1)
        printf1(ctx, "Error: insufficient memory.\n");
    return(err);
}
