/****************************************************************************/
/*  t_ass                                                                   */
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
#include "t_gio.h"
#include "t_ml.h"
#include "tda_context.h"

/*  functions in t_ass.c */

int gap(TDAContext *ctx);
void g_aloc(TDAContext *ctx, int n,int *a,int *c,int *t,int *ch,int *lc,int *lr,int *lz,int *nz, int *slc,int *slr,int *rh,int *u);
int gqap(TDAContext *ctx);
int g_prnp(TDAContext *ctx, int n,int *perm);
int gqapd(TDAContext *ctx, int n,int niter,double alpha,double beta,int look4,int *seed, int *f,int *d,int *a,int *b,int *srtf,int *srtif,int *srtd,int *srtid, int *srtc,int *srtic,int *indexd,int *indexf,int *cost,int *fdind, int *opta,int *bestv,int *iter,int *ito);
void a_srtcst(TDAContext *ctx, int n,int n2,double beta,int *f,int *d,int *srtf,int *srtif, int *srtd,int *srtid,int *srtc,int *srtic,int *indexd,int *indexf, int *cost,int *fdind);
void a_stage1(TDAContext *ctx, int n,int n2,int *i,int *j,int *k,int *l,int *seed,double alpha, double beta,int *objv,int *indexd,int *indexf,int *fdind,int *cost, int *a,int *b);
void a_stage2(TDAContext *ctx, int n,int n2,int i,int j,int k,int l,int *seed,int *objv, double alpha,int *f,int *d,int *srtc,int *srtic,int *a,int *b);
void a_savsol(TDAContext *ctx, int n,int objv,int *bestv,int *a,int *opta);
void a_local(TDAContext *ctx, int n,int n2,int *objv,int *f,int *d,int *a,int *b);
void a_mkbseq(TDAContext *ctx, int n,int *a,int *b);
void a_insrtq(TDAContext *ctx, int n2,int v,int iv,int *sizeq,int *q,int *iq);
void a_removq(TDAContext *ctx, int n2,int *v,int *iv,int *sizeq,int *q,int *iq);
double a_randp(TDAContext *ctx, int *ix);
void a_evalij(TDAContext *ctx, int n,int n2,int i,int j,int* xgain,int *f,int *d,int *a);
int gloc(TDAContext *ctx);
int g_locate(TDAContext *ctx, int nvpt,int nfpt,int *vwt,int *vfwt,int *post);    
void g_netflo(TDAContext *ctx, int n,int e,int *endp1,int *endp2,int *endp3,int *endp4, int source,int sink,int *cut,int *cap,int *aux);
int qap_w(TDAContext *ctx, int n,double *c,double *f,double *d,int *loc3n,double *wtmp);


/* ------------------------------------------------------------------------ */
/*  gap         Column permutations. Graph must be integer-valued.          */
/*                                                                          */
/*              Uses CACM 548.                                              */ 
/*                                                                          */
/*              gap(                                                        */
/*                  gn = ...,       g1,g2, graph numbers                    */
/*                  nfmt = ...,     integer format, def. 4                  */
/*                  df=...,         write permuted matrix                   */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gap(TDAContext *ctx)
{
    register int i,j;
    int err,n,n1,nr,t;
    int *ch,*lc,*lr,*lz,*nz,*slc,*slr,*rh,*u;
    double tmp;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Column permutations. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto GAPFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GAPFin;
    if (gdd_tcheck(ctx, 0,0,0))      
        goto GAPFin;

    if (ctx->MxIter < 1 || ctx->MxItFlg == 0)
        ctx->MxIter = 100;
    if (ctx->PMIDFA <= 0.0 || ctx->PMIDFA > 1.0)
        ctx->PMIDFA = 0.25;
    if (ctx->PMIDFB <= 0.0 || ctx->PMIDFB > 1.0)
        ctx->PMIDFB = 0.5;

    if (ctx->PMALG != 1)
        ctx->PMALG = 1;
    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 10,4);
     
    n = ctx->GD_NP;
    n1 = n + 1;  

    if (alloc_acn(ctx, n * n1 + 1))      /* a matrix */
        goto GAPFin;
    if (alloc_acm(ctx, n + 1))          /* c vector */
        goto GAPFin;
    
    nr = 0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            tmp = gdd_adj(ctx, i,j,ctx->PMGN);
            if (tmp < 0.0) {
                tmp = 0.0;
                if (i != j)
                    nr++;
            }
            ctx->AcN[j * n1 + i + 1] = (int)tmp;
        }
    }
    if (nr > 0)
        printf1(ctx, "Substituted %d missing values.\n",nr);

    if (alloc_acr(ctx, 7 * n + 2 * n1 + 2))     /* working storage */
        goto GAPFin;

    ch = ctx->AcR;
    lc = ch + n;
    lr = lc + n;
    lz = lr + n;
    nz = lz + n;
    slc = nz + n;
    slr = slc + n;
    rh  = slr + n;
    u   = rh  + n1;
   
    g_aloc(ctx, n,ctx->AcN,ctx->AcM,&t,ch,lc,lr,lz,nz,slc,slr,rh,u); 

    printf1(ctx, "Minimal cost value: %d\n",t);
    for (i = 1; i <= n; ++i) {
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcM[i]);
        fprintf(ctx->PMFd,"\n");
    }
    printf1(ctx, "%d records written to: %s\n",n,ctx->PMFdName);

    if (ctx->PMF1Def) {                       /* write permuted matrix */

        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; ++j) {
                nr = (int)gdd_adj(ctx, i,ctx->AcM[j + 1] - 1,ctx->PMGN);
                if (nr < 0)  
                    nr = 0;
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,nr);
            }
            fprintf(ctx->PMF1d,"\n");
        }
        printf1(ctx, "%d records written to: %s\n",n,ctx->PMF1dName);
    }
    err = 0;

GAPFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_aloc()                                                                */
/*                                                                          */
/*  This function solves the simple assignment problem.                     */
/*  The code is adapted from G. Carpaneto, P. Toth, Algorithm 548 of ACM,   */
/*  Solution of the assignment problem. ACM Transactions on Mathematical    */
/*  Software 6 (1980), 104 - 111.                                           */
/*                                                                          */
/*  n = number of nodes                                                     */
/*  a(i,j) = (n,n+1) matrix                                                 */
/*  c(j) = row assigned to column j, j = 1,...,n                            */
/*  t = cost of optimal assignment                                          */
/*  ch[j], j = 1,...,n                                                      */
/*  lc[j], j = 1,...,n                                                      */
/*  lr[j], j = 1,...,n                                                      */
/*  lz[j], j = 1,...,n                                                      */
/*  nz[j], j = 1,...,n                                                      */
/*  slc[j],j = 1,...,n                                                      */
/*  slr[j],j = 1,...,n                                                      */
/*  rh[j], j = 1,...,n+1                                                    */
/*  u[j],  j = 1,...,n+1                                                    */
/*                                                                          */

void g_aloc(TDAContext *ctx, int n,int *a,int *c,int *t,int *ch,int *lc,int *lr,int *lz,int *nz, int *slc,int *slr,int *rh,int *u)
{
    register int i,j,k,l,r;
    int q,s,np1,m,nm,lj,lm,nl,h,kslc,kslr;


    /* initialization */

    np1 = n + 1;
    for (j = 1; j <= n; ++j) {
        c[j] = lz[j] = nz[j] =  u[j] = 0;
    }

    u[np1] = 0;
    *t = 0;

    /* reduction of the initial cost matrix */ 

    for (j = 1; j <= n; ++j) {
        s = a[j];  
        for (l = 2; l <= n; ++l) {
            if (a[(l-1) * np1 + j] < s)
                s = a[(l-1) * np1 + j];
        }
        t += s; 
        for (i = 1; i <= n; ++i)
          a[(i-1) * np1 + j] -= s;       
    }
    for (i = 1; i <= n; ++i) {
        q = a[(i-1) * np1 + 1];
        for (l = 2; l <= n; ++l) {
            if (a[(i-1) * np1 + l] < q)
                q = a[(i-1) * np1 + l];
        }
        *t += q; 
        l = np1;
        for (j = 1; j <= n; ++j) {
            a[(i-1) * np1 + j] -= q;            
            if (a[(i-1) * np1 + j] == 0 ) {
                a[(i-1) * np1 + l] = -j;
                l = j;
            }
        }
    }

    /* choice of the initial solution */

    k = np1;
    for (i = 1; i <= n; ++i) {
        lj = np1;
        j = -a[(i-1) * np1 + np1];
L80:
        if (c[j] != 0) {
            lj = j;
            j = -a[(i-1) * np1 + j];
            if (j != 0)
                goto L80; 
            lj = np1;
            j = -a[(i-1) * np1 + np1];
L90:
            r = c[j];
            lm = lz[r];
            m = nz[r];

            while (m != 0) {
                if (c[m] == 0)
                    goto L120;
                lm = m;
                m = -a[(r-1) * np1 + m];
            }
            lj = j;
            j = -a[(i-1) * np1 + j];
            if (j != 0)
                goto L90;
            u[k] = i;
            k = i;
            continue;  

L120:
            nz[r] = -a[(r-1) * np1 + m];
            lz[r] = j;
            a[(r-1) * np1 + lm] = -j;
            a[(r-1) * np1 + j] = a[(r-1) * np1 + m];
            a[(r-1) * np1 + m] = 0;
            c[m] = r;
        }
        c[j] = i;
        a[(i-1) * np1 + lj] = a[(i-1) * np1 + j];
        nz[i] = -a[(i-1) * np1 + j];
        lz[i] = lj;
        a[(i-1) * np1 + j] = 0;
    }

    /* research of a new assignment */

L150:
    if (u[np1] == 0)
        return;

    for (i = 1; i <= n; ++i)  
        ch[i] = lc[i] = lr[i] = rh[i] = 0;
     
    rh[np1] = -1;
    kslc = 0;
    kslr = 1;
    r = u[np1];
    lr[r] = -1;
    slr[1] = r;
    if (a[(r-1) * np1 + np1] == 0)
        goto L220;

L170:
    l = -a[(r-1) * np1 + np1];
    if (a[(r-1) * np1 + l] != 0) {
        if (rh[r] == 0) {
            rh[r] = rh[np1];
            ch[r] = -a[(r-1) * np1 + l];
            rh[np1] = r;
        }
    }

L180:
    if (lc[l] == 0)
        goto L200;
    if (rh[r] == 0)
        goto L210;

L190:
    l = ch[r];
    ch[r] = -a[(r-1) * np1 + l];
    if (a[(r-1) * np1 + l] != 0)
        goto L180;

    rh[np1] = rh[r];
    rh[r] = 0;
    goto L180;

L200:
    lc[l] = r;
    if (c[l] == 0) {   /* assignment of a new row */
      
        while (1) {
            c[l] = r;
            m = np1;
       
            while (1) {
                nm = -a[(r-1) * np1 + m];
                if (nm == l)
                    break;    
                m = nm;
            }
            a[(r-1) * np1 + m] = a[(r-1) * np1 + l];
            a[(r-1) * np1 + l] = 0;
            if (lr[r] < 0)
                break;       

            l = lr[r];
            a[(r-1) * np1 + l] = a[(r-1) * np1 + np1];
            a[(r-1) * np1 + np1] = -l;
            r = lc[l];
        }   
        u[np1] = u[r];
        u[r] = 0; 
        goto L150;
    }
    kslc += 1;    
    slc[kslc] = l;
    r = c[l];
    lr[r] = l;
    kslr += 1;     
    slr[kslr] = r;
    if (a[(r-1) * np1 + np1] != 0)
        goto L170;

L210:
    if (rh[np1] > 0)
        goto L350;

    /* reduction of the current cost matrix */

L220:
    h = ctx->INTMAX; 
    for (j = 1; j <= n; ++j) {

        if (lc[j] != 0)
            continue;  

        for (k = 1; k <= kslr; ++k) {
            i = slr[k];
            if (a[(i-1) * np1 + j] < h)
                h = a[(i-1) * np1 + j];
        }
    }

    *t += h; 
    for (j = 1; j <= n; ++j) {

        if (lc[j] != 0)
            continue;  

        for (k = 1; k <= kslr; ++k) {
            i = slr[k];
            a[(i-1) * np1 + j] -= h;        
            if (a[(i-1) * np1 + j] != 0)
                continue; 
            if (rh[i] == 0) {
                rh[i] = rh[np1];
                ch[i] = j;
                rh[np1] = i;
            }
            l = np1;
      
            while (1) {
                nl = -a[(i-1) * np1 + l];
                if (nl == 0)
                    break;     
                l = nl;
            }
            a[(i-1) * np1 + l] = -j;
        }
    }

    if (kslc != 0) {

        for (i = 1; i <= n; ++i) {

            if (lr[i] != 0)
                continue; 

            for (k = 1; k <= kslc; ++k) {
                j = slc[k];
                if (a[(i-1) * np1 + j] > 0) {
                    a[(i-1) * np1 + j] += h;       
                }
                else {
                    l = np1;

                    while (1) {
                        nl = - a[(i-1) * np1 + l];
                        if (nl == j)
                            break;     
                        l = nl;
                    }
                    a[(i-1) * np1 + l] = a[(i-1) * np1 + j];
                    a[(i-1) * np1 + j] = h;
                }          
            }
        }
    }
L350:
    r = rh[np1];
    goto L190;
}

/* -##--------------------------------------------------------------------- */
/*  gqap        Quadratic assignment. Graph is always undirected.           */
/*                                    and must be integer-valued.           */
/*                                                                          */
/*              gqap(                                                       */
/*                  alg = ...,      algorithm, def. 1                       */
/*                                  1 = GRASP for dense graphs, CACM 754    */
/*                  gn = ...,       g1,g2, graph numbers (F,D)              */
/*                  mxit=...,       max iterations, def. 100                */
/*                  idf=...,        alpha,beta, def. alpha=0.25,beta=0.5    */  
/*                  nfmt = ...,     integer format, def. 4                  */
/*                  df=...,         write graphs, F permuted                */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gqap(TDAContext *ctx)
{
    register int i,j,k,l;
    int err,n,nn,n1,n2,iter,ito,seed,bestv,look4;
    int *a,*b,*srtf,*srtif,*srtd,*srtid,*srtc,*srtic,*indexf,*indexd,*cost,*fdind;
    double tmp;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Quadratic assignment. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto GQAPFin;

    if (gdd_check2(ctx, ctx->PMGN,ctx->PMGN1,0))
        goto GQAPFin;
    /**
    if (gdd_tcheck(ctx, 0,1,0))  
        goto GQAPFin;
    **/

    if (ctx->MxIter < 1 || ctx->MxItFlg == 0)
        ctx->MxIter = 100;
    if (ctx->PMIDFA <= 0.0 || ctx->PMIDFA > 1.0)
        ctx->PMIDFA = 0.25;
    if (ctx->PMIDFB <= 0.0 || ctx->PMIDFB > 1.0)
        ctx->PMIDFB = 0.5;

    if (ctx->PMALG != 1)
        ctx->PMALG = 1;
    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 10,4);
     
    n = ctx->GD_NP;
    nn = n * n;

    if (alloc_acn(ctx, nn + 1))          /* f matrix */
        goto GQAPFin;
    if (alloc_acm(ctx, nn + 1))          /* d matrix */
        goto GQAPFin;
    
    n1 = n2 = k = l = 0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            tmp = gdd_adj(ctx, i,j,ctx->PMGN);
            if (tmp < 0.0) {
                tmp = 0.0;
                if (i != j)
                    n1++;
            }
            ctx->AcN[++k] =  (int)tmp;

            tmp = gdd_adj(ctx, i,j,ctx->PMGN1);
            if (tmp < 0.0) {
                tmp = 0.0;
                if (i != j)
                    n2++;
            }
            ctx->AcM[++l] = (int)tmp;
        }
    }
    printf1(ctx, "F matrix: graph %d.",ctx->PMGN);
    if (n1 > 0)
        printf1(ctx, " Substituted %d missing values.",n1);
    printf1(ctx, "\nD matrix: graph %d.",ctx->PMGN1);
    if (n2 > 0)
        printf1(ctx, " Substituted %d missing values.",n2);

    printf1(ctx, "\nMax number of iterations: %d (alpha %lg, beta %lg).\n\n",
        ctx->MxIter,ctx->PMIDFA,ctx->PMIDFB);

    /*******
    printf1(ctx, "F\n");
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) 
            printf1(ctx, "%2d ",AcN[(i - 1) * n + j]);
        newline(ctx);
    }

    printf1(ctx, "D\n");
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) 
            printf1(ctx, "%2d ",AcM[(i - 1) * n + j]);
        newline(ctx);
    }
    *****/

    if (alloc_acs(ctx, n + 1))           /* optimal permutation vector */
        goto GQAPFin;

    if (alloc_acr(ctx, 10 * nn + 2 * n + 2))     /* working storage */
        goto GQAPFin;

    a = ctx->AcR;
    b = a + n;
    srtf = b + n;
    srtif = srtf + nn;
    srtd  = srtif + nn;
    srtid = srtd + nn;
    srtc  = srtid + nn;
    srtic = srtc + nn;
    indexf = srtic + nn;
    indexd = indexf + nn;
    cost = indexd + nn;
    fdind = cost + nn;
    seed = 270001;
    look4 = -1;         /* NOT USED ! */ 
   
    gqapd(ctx, n,ctx->MxIter,ctx->PMIDFA,ctx->PMIDFB,look4,&seed,ctx->AcN,ctx->AcM,a,b,srtf,srtif,srtd,srtid,
          srtc,srtic,indexd,indexf,cost,fdind,ctx->AcS,&bestv,&iter,&ito);

    printf1(ctx, "Performed %d iterations.\n",iter);
    printf1(ctx, "Best cost value: %d  (found in %d iterations).\n",bestv,ito);

#ifdef TDA_R_PACKAGE
    /* direct exports (CONTRIBUTING.md): the optimum and the
       assignment permutation */
    tda_export_cell(ctx, "gqap.cost", (double)bestv);
    tda_export_endrow(ctx, "gqap.cost");
    for (i = 1; i <= n; ++i) {
        tda_export_cell(ctx, "gqap.assignment", (double)ctx->AcS[i]);
        tda_export_endrow(ctx, "gqap.assignment");
    }
    tda_export_flush(ctx, "gqap.cost");
    tda_export_flush(ctx, "gqap.assignment");
#endif
    n1 = 0;
    for (i = 1; i <= n; ++i) {
        k = ctx->AcS[i];
        for (j = 1; j <= n; ++j) {
            l = ctx->AcS[j];
            n1 += ctx->AcM[(i - 1) * n + j] * ctx->AcN[(k - 1) * n + l];
        }
    }
    printf1(ctx, "Checked cost value: %d\n",n1);

    for (i = 1; i <= n; ++i) {
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcS[i]);
        fprintf(ctx->PMFd,"\n");
    }
    printf1(ctx, "%d records written to: %s\n",n,ctx->PMFdName);

    if (ctx->PMF1Def) {                       /* write second output file */
        n1 = g_prnp(ctx, n,ctx->AcS + 1);
        printf1(ctx, "%d records written to: %s\n",n1,ctx->PMF1dName);
    }
    err = 0;

GQAPFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_prnp(n,perm)  write first two graphs to PMF1d. Permute first graph    */
/*                  according to perm[]. Return number of records written.  */

int g_prnp(TDAContext *ctx, int n,int *perm)
{
    register int i,j,k,l,v;

    for (i = 0; i < n; ++i) {
        k = perm[i] - 1;
        for (j = 0; j < n; ++j) {
            l = perm[j] - 1;

            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,i + 1);
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,j + 1);
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,gdd_node(ctx, i));
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,gdd_node(ctx, j));

            v = (int)gdd_adj(ctx, k,l,ctx->PMGN);
            if (v < 0)
                v = 0;
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,v);

            v = (int)gdd_adj(ctx, i,j,ctx->PMGN1);
            if (v < 0)
                v = 0;
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,v);
            fprintf(ctx->PMF1d,"\n");
        }
    }
    return(n * n);
}

/* ------------------------------------------------------------------------ */    
/*  A Greedy Randomized Adaptive Search Procedure (GRASP) for the           */
/*  Quadratic Assignment Problem (QAP)                                      */
/*                                                                          */
/*  Authors: M.G.C. Resende (AT&T Bell Laboratories)                        */
/*           [mgcr@research.att.com]                                        */
/*           Y. Li (Penn State University)                                  */
/*           [yong@cs.psu.edu]                                              */
/*           P.M. Pardalos (University of Florida)                          */
/*           [pardalos@ufl.edu]                                             */
/*                                                                          */
/*  gqapd(n,n2,niter,alpha,beta,look4,seed,f,d,a,b,                         */
/*          srtf,srtif,srtd,srtid,srtc,srtic,indexd,                        */
/*          indexf,cost,fdind,opta,bestv,iter)                              */
/*                                                                          */
/*  Subroutine for finding an approximate solution of a                     */
/*  dense symmetric quadratic assignment problem.                           */
/*                                                                          */
/*  infty - a large integer                                                 */
/*                                                                          */
/*  Passed input scalars:                                                   */
/*                                                                          */
/*         n      - dimension of qap problem                                */
/*         n2     - n * n                                                   */
/*         niter  - maximum number of GRASP iterations                      */
/*         alpha  - phase 1 parameter                                       */
/*         beta   - phase 1 parameter                                       */
/*         look4  - if permutation of cost look4 or less is found gqapd     */
/*                  returns that permutation                                */
/*                                                                          */
/*  Passed input/output scalar:                                             */
/*                                                                          */
/*  seed   - random number generator seed                                   */
/*                                                                          */
/*  Passed output scalars:                                                  */
/*                                                                          */
/*  bestv  - cost of best assignment found                                  */
/*  iter   - number of GRASP iterations taken                               */
/*  ito    - number of iterations used for best solution                    */
/*                                                                          */
/*  Passed input arrays:                                                    */
/*                                                                          */
/*  f      - flow matrix stored as a 1-dimensional array,                   */
/*           row by row (dim = n2).                                         */
/*  d      - distance matrix stored as a 1-dimensional array,               */
/*           row by row (dim = n2).                                         */
/*                                                                          */
/*  Passed work arrays:                                                     */
/*                                                                          */
/*  a      - permutation vector (dim = n).                                  */
/*  b      - permutation vector (dim = n).                                  */
/*  srtf   - sorted F values                                                */
/*  srtif  - sorted F values (indices)                                      */
/*  srtd   - sorted D values                                                */
/*  srtid  - sorted D values (indices)                                      */
/*  srtc   - sorted cost values                                             */
/*  srtic  - sorted cost values (indices)                                   */
/*  indexf - indices of facilities in unsorted cost matrix                  */
/*  indexd - indices of locations in unsorted cost matrix                   */
/*  cost   - sorted cost matrix                                             */
/*  fdind  - indices of sorted cost matrix                                  */
/*                                                                          */
/*  Passed output array:                                                    */
/*                                                                          */
/*  opta   - best permutation vector (dim = n).                             */
/*                                                                          */


int gqapd(TDAContext *ctx, int n,int niter,double alpha,double beta,int look4,int *seed, int *f,int *d,int *a,int *b,int *srtf,int *srtif,int *srtd,int *srtid, int *srtc,int *srtic,int *indexd,int *indexf,int *cost,int *fdind, int *opta,int *bestv,int *iter,int *ito)
{
    (void)look4;        /* unused: the signature is shared */
    int n2,objv,i,j,k,l;
              
    n2 = n * n;
    *bestv = ctx->INTMAX;

    /* Sort the cost = f(i,j) * d(k,l) in increasing order to be used 
       by the stage1 construction phase of GRASP. */

    a_srtcst(ctx, n,n2,beta,f,d,srtf,srtif,srtd,srtid,srtc,
             srtic,indexd,indexf,cost,fdind);
   
    /* Do GRASP iterations. */

    for (*iter = 1; *iter <= niter; ++*iter) {

        /* Stage 1 of GRASP construction phase. */

        a_stage1(ctx, n,n2,&i,&j,&k,&l,seed,alpha,beta,&objv,indexd,indexf,
                 fdind,cost,a,b);
   
        /* Stage 2 of GRASP construction phase. */

        a_stage2(ctx, n,n2,i,j,k,l,seed,&objv,alpha,f,d,srtc,srtic,a,b);
   
        /* Local search phase of GRASP. */

        a_local(ctx, n,n2,&objv,f,d,a,b);

        /* If cost assignment is best so far, save permutation and
           cost of assignment. */

        if (objv < *bestv) {        

            if (ctx->PMProtFDef)  
                fprintf(ctx->PMProtFd,"Iter %5d  value %d\n",*iter,objv);

            *ito = *iter;

            a_savsol(ctx, n,objv,bestv,a,opta);

            /* If cost of assignment is at least as good as requested,
               return best permutation found. */

            /* NOT USED ********
            if (*bestv <= look4)
                return(0);
            ********************/
        }        
    }
    *iter = niter;
    return(0);
}               

/* ------------------------------------------------------------------------ */    
/*  a_srtcst(n,n2,beta,f,d,srtf,srtif,srtd,srtid,srtc,                      */
/*         srtic,indexd,indexf,cost,fdind)                                  */
/*                                                                          */
/*  Sorts cost = f(i,j)*d(k,l) in increasing order.                         */
/*                                                                          */
/*  Passed input scalars:                                                   */
/*                                                                          */
/*  n      - qap dimension                                                  */
/*  n2     - n * n                                                          */
/*  beta   - construction phase parameter                                   */
/*                                                                          */
/*  Passed input arrays:                                                    */
/*                                                                          */
/*  f      - flow matrix (row major order)                                  */
/*  d      - distance matrix (row major order)                              */
/*                                                                          */
/*  Passed work arrays:                                                     */
/*                                                                          */
/*  srtf   - sorted flow matrix (values)                                    */
/*  srtif  - sorted flow matrix (indices)                                   */
/*  srtd   - sorted distance matrix (values)                                */
/*  srtid  - sorted distance matrix (indices)                               */
/*  srtc   - sorted cost matrix (values)                                    */
/*  srtic  - sorted cost matrix (indices)                                   */
/*                                                                          */
/*  Passed output arrays:                                                   */
/*                                                                          */
/*  indexd - indices of locations in unsorted cost matrix                   */
/*  indexf - indices of facilities in unsorted cost matrix                  */
/*  cost   - sorted cost matrix                                             */
/*  fdind  - indices of sorted cost matrix                                  */
/*                                                                          */
/*  Sort D in increasing order,                                             */
/*  F in decreasing order (-F in increasing order).                         */
/*  Keep only the (n*n-n)*beta best elements in each sorting.               */
 

void a_srtcst(TDAContext *ctx, int n,int n2,double beta,int *f,int *d,int *srtf,int *srtif, int *srtd,int *srtid,int *srtc,int *srtic,int *indexd,int *indexf, int *cost,int *fdind)
{
    int index,sizec,sized,sizef,dv,fv,dind,find,nbeta,i,j;

    /* Initialize cardinalities of sorted sets of elements of D, F,
       and cost */

    sized = sizef = sizec = 0;

    /* Insert all non-diagonal elements of D into D-priority heap
       and all non-diagonal elements of -F into F-priority heap. */

    index = 0;
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) {
            index++;       
            if (i != j) {     
                a_insrtq(ctx, n2,d[index],index,&sized,srtd,srtid);
                a_insrtq(ctx, n2,-f[index],index,&sizef,srtf,srtif);       

            }          
        }
    }

    /* Compute size of sorted sets. */

    nbeta = (int)(beta * (double)(n * n - n));
  
    /* Remove the nbeta smallest D elements from D-priority heap and
       the nbeta smallest -F elements from F-priority heap. */

    for (i = 1; i <= nbeta; ++i) {

        a_removq(ctx, n2,&dv,&dind,&sized,srtd,srtid);
        a_removq(ctx, n2,&fv,&find,&sizef,srtf,srtif);

        /* Cost is product of sorted flow and distance. */

        cost[i] = -dv * fv;
        indexd[i] = dind;
        indexf[i] = find;

        /* Insert cost into cost priority-heap. */

        a_insrtq(ctx, n2,cost[i],i,&sizec,srtc,srtic);
    }

    /* Remove nbeta sorted cost elements from cost priority-heap. */              

    for (i = 1; i <= nbeta; ++i) {
        a_removq(ctx, n2,&cost[i],&fdind[i],&sizec,srtc,srtic);
    }
}             

/* ------------------------------------------------------------------------ */    
/*  a_stage1(n,n2,i,j,k,l,seed,alpha,beta,objv,indexd,indexf,               */
/*         fdind,cost,a,b)                                                  */
/*                                                                          */
/*  stage1:  Builds the initial 2 assignments for the GRASP                 */
/*           construction phase (facility i to site k and                   */
/*           facility j to site l).                                         */
/*                                                                          */
/*  Passed input scalars:                                                   */
/*                                                                          */
/*  n      - qap dimension                                                  */
/*  n2     - n * n                                                          */
/*  alpha  - construction phase parameter                                   */
/*  beta   - construction phase parameter                                   */
/*                                                                          */
/*  Passed input/output scalar:                                             */
/*                                                                          */
/*  seed   - random number generator seed                                   */
/*                                                                          */
/*  Passed output scalars:                                                  */
/*                                                                          */
/*  i      - facility index                                                 */
/*  j      - facility index                                                 */
/*  k      - location index                                                 */
/*  l      - location index                                                 */
/*  objv   - cost of initial 2 assignments                                  */
/*                                                                          */
/*  Passed input arrays:                                                    */
/*                                                                          */
/*  indexd - indices of locations in unsorted cost matrix                   */
/*  indexf - indices of facilities in unsorted cost matrix                  */
/*  fdind  - indices of sorted cost matrix                                  */
/*  cost   - cost of assignment                                             */
/*                                                                          */
/*  Passed output arrays:                                                   */
/*                                                                          */
/*  a      - permutation array                                              */
/*  b      - permutation array                                              */
/*                                                                          */
                                                                             
void a_stage1(TDAContext *ctx, int n,int n2,int *i,int *j,int *k,int *l,int *seed,double alpha, double beta,int *objv,int *indexd,int *indexf,int *fdind,int *cost, int *a,int *b)
{
    int nselct,dind,find,high,ii,tmp;     

    /* Initialize permutations. */

    for (ii = 1; ii <= n; ++ii)  
        a[ii] = b[ii] = ii;

    /* Select element, at random, from the best (n*n-n)*alpha cost 
       elements. */

    a_randp(ctx, seed);
    high = (int)(alpha * beta * (n * n - n));
    if (high <= 0)  
        high = 1;
       
    nselct = 1 + *seed / (2147483647 / high);
    if (nselct > n2 - n)
        nselct = n2 - n;

    /* Initial assignment is facility i to location k
       facility j to location l. */

    dind = indexd[fdind[nselct]];
    find = indexf[fdind[nselct]];
    *i = (find - 1) / n + 1;
    *j = find - (*i - 1) * n;
    *k = (dind - 1) / n + 1;
    *l = dind - (*k - 1) * n;
    *objv = cost[nselct];

    /* Make initial assignments to permutation arrays:
       Assign facility i to location k. */

    a[1] = *i;
    a[*i] = 1;
    b[1] = *k;
    b[*k] = 1;

    /* Assign facility j to location l. */

    for (ii = 1; ii <= n; ++ii) {
        if (a[ii] == *j) {     
            tmp = a[2];
            a[2] = *j;
            a[ii] = tmp;
            break;    
        }        
    }
    for (ii = 1; ii <= n; ++ii) {
        if (b[ii] == *l) {      
            tmp = b[2];
            b[2] = *l;
            b[ii] = tmp;
            break;      
        }         
    }
}            

/* ------------------------------------------------------------------------ */    
/*  a_stage2(n,n2,i,j,k,l,seed,objv,alpha,f,d,srtc,srtic,a,b)               */
/*                                                                          */
/*  Builds a randomized greedy permutation starting from                    */
/*  the assignments made in stage1.                                         */
/*  Permutation is returned in array a(*).                                  */
/*                                                                          */
/*  Passed input scalars:                                                   */
/*                                                                          */
/*  n      - problem dimension                                              */
/*  n2     - n * n                                                          */
/*  i      - facility index                                                 */
/*  j      - location index                                                 */
/*  k      - facility index                                                 */
/*  l      - location index                                                 */
/*  alpha  - construction phase parameter                                   */
/*                                                                          */
/*  Passed input/output scalars:                                            */
/*                                                                          */
/*  seed   - random number generator seed                                   */
/*  objv   - cost of assignment                                             */
/*                                                                          */
/*  Passed input arrays:                                                    */
/*                                                                          */
/*  f      - flow matrix                                                    */
/*  d      - distance matrix                                                */
/*                                                                          */
/*  Passed work arrays:                                                     */
/*                                                                          */
/*  srtc   - sorted cost matrix (values)                                    */
/*  srtic  - sorted cost matrix (indices)                                   */
/*                                                                          */
/*  Passed input/output arrays:                                             */
/*                                                                          */
/*  a      - permutation array                                              */
/*  b      - permutation array                                              */
/*                                                                          */

void a_stage2(TDAContext *ctx, int n,int n2,int i,int j,int k,int l,int *seed,int *objv, double alpha,int *f,int *d,int *srtc,int *srtic,int *a,int *b)
{
    int high,assign,cost,sizec,nselct,tmp,kinv,linv,fdind;
    int akm1tn,blm1tn,anm1tn,bnm1tn;

    /* Main loop:  Assignments 3,4,..,n-1 are made. */

    for (assign = 3; assign < n; ++assign) {

        /* For all pairs not assigned yet, compute costs of all possible
           assignments, w.r.t. already-made assignments. */

        sizec = 0;

        for (k = assign; k <= n; ++k) {

            akm1tn = (a[k] - 1) * n;

            for (l = assign; l <= n; ++l) {

                blm1tn = (b[l] - 1) * n;
                cost = 0;

                for (i = 1; i < assign; ++i) {

                    /* Facility a(i) already assigned to location b(i):
                       Cost of assigning facility a(k) to location b(l)
                       relative to assignment of  facility a(i) to 
                       location b(i) */

                    cost += f[akm1tn + a[i]] * d[blm1tn + b[i]];
                }

                /* Insert cost element into cost-priority heap for 
                   sorting. */

                a_insrtq(ctx, n2,cost,akm1tn + b[l],&sizec,srtc,srtic);
            }
        }
 
        /* Select assignment, at random, from the best alpha*sizec 
           assignments. */

        a_randp(ctx, seed);
        high = (int)(alpha * sizec);
        if (high < 1)  
            high = 1;
  
        nselct = 1 + *seed / (2147483647 / high);

        for (i = 1; i <= nselct; ++i)  
            a_removq(ctx, n2,&cost,&fdind,&sizec,srtc,srtic);

        /* make assignment. */

        *objv += cost;
        kinv = (fdind - 1) / n + 1;
        linv = fdind - (kinv - 1) * n;

        for (i = assign; i <= n; ++i) {
            if (a[i] == kinv) {    
                k = i;
                break;          
            }           
        }
        for (j = assign; j <= n; ++j) {
            if (b[j] == linv) {      
                l = j;
                break;             
            }           
        }
        tmp = a[assign];
        a[assign] = a[k];
        a[k] = tmp;
        tmp = b[assign];
        b[assign] = b[l];
        b[l] = tmp;
    }
    anm1tn = (a[n] - 1) * n;
    bnm1tn = (b[n] - 1) * n;

    for (i = 1; i < n; ++i)  
        *objv += f[anm1tn + a[i]] * d[bnm1tn + b[i]];

    *objv *= 2;               
}

/* ------------------------------------------------------------------------ */    
/*  a_savsol(n,objv,bestv,a,opta)                                           */
/*                                                                          */
/*  Saves current best solution.                                            */
/*                                                                          */
/*  Passed input scalars:                                                   */
/*  n      - problem dimension                                              */
/*  objv   - objective function value                                       */
/*                                                                          */
/*  Passed output scalars:                                                  */
/*  bestv  - best objective function value so far                           */
/*                                                                          */
/*  Passed input array:                                                     */
/*  a      - permutation array                                              */
/*                                                                          */
/*  Passed output array:                                                    */
/*  opta   - array of best permutation so far                               */

void a_savsol(TDAContext *ctx, int n,int objv,int *bestv,int *a,int *opta)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;
    
    for (i = 1; i <= n; ++i)    
        opta[i] = a[i];
    *bestv = objv;
}

/* ------------------------------------------------------------------------ */    
/*  a_local(n,n2,objv,f,d,a,b)                                              */
/*                                                                          */
/*  Local 2-exchange on permutation array a.                                */
/*  Return improved permutation array a and objv.                           */
/*                                                                          */
/*  Passed input scalars:                                                   */
/*  n      - problem dimension                                              */
/*  n2     - n * n                                                          */
/*                                                                          */
/*  Passed input/output scalar:                                             */
/*  objv   - objective function value                                       */
/*                                                                          */
/*  Passed input arrays:                                                    */
/*  f      - flow matrix                                                    */
/*  d      - distance matrix                                                */
/*                                                                          */
/*  Passed input/output arrays:                                             */
/*  a      - permutation array                                              */
/*  b      - permutation array                                              */

void a_local(TDAContext *ctx, int n,int n2,int *objv,int *f,int *d,int *a,int *b)
{
    int i,j,temp,xgain,improv;

    /* Make array b(*) = (1,2,3,...,n) for local search. */

    a_mkbseq(ctx, n,a,b);

    /* Attempt to switch all pairs in permutation array a. */

LA10: 
    improv = 0;     

    for (i = 1; i < n; ++i) {
        for (j = i + 1; j <= n; ++j) {

            /* Evaluate cost difference by adopting switch of a(i) 
               and a(j). */
    
            a_evalij(ctx, n,n2,i,j,&xgain,f,d,a);
   
            /* If switch improves cost, adopt it. */

            if (xgain > 0) {     
                temp = a[i];
                a[i] = a[j];
                a[j] = temp;
                *objv -= xgain;
                improv = 1;   
            }                 
        }
    }

    /* If no switch improves cost (improv=.false.), return; else repeat. */

    if (improv)
        goto LA10;
}

/* ------------------------------------------------------------------------ */    
/*  a_mkbseq(n,a,b)                                                         */
/*                                                                          */
/*  Change permutation arrays a and b to make b = (1,2,...,n).              */
/*                                                                          */
/*  Passed input scalar:                                                    */
/*  n      - QAP dimension                                                  */
/*                                                                          */
/*  Passed input/output arrays:                                             */
/*  a      - permutation array                                              */
/*  b      - permutation array                                              */

void a_mkbseq(TDAContext *ctx, int n,int *a,int *b)
{
    (void)ctx;        /* unused: the signature is shared */
    int i,j,tmp;

    for (i = 1; i < n; ++i) {
        for (j = i + 1; j <= n; ++j) {
            if (b[j] == i) {     
                b[j] = b[i];
                b[i] = i;
                tmp = a[i];
                a[i] = a[j];
                a[j] = tmp;
                break;  
            }              
        }
    }
}

/* ------------------------------------------------------------------------ */    
/*  a_insrtq(n2,v,iv,sizeq,q,iq)                                            */
/*                                                                          */
/*  Insert an element (v,iv) into a queue (q,iq).                           */
/*                                                                          */
/*  Passed input scalars:                                                   */
/*  n2     - n * n                                                          */
/*  v      - heap element (value)                                           */
/*  iv     - heap element (index)                                           */
/*                                                                          */
/*  Passed input/output scalar:                                             */
/*  sizeq  - size of heap                                                   */
/*                                                                          */
/*  Passed input/output arrays:                                             */
/*  q      - heap (value)                                                   */
/*  iq     - heap (index)                                                   */
/*                                                                          */

void a_insrtq(TDAContext *ctx, int n2,int v,int iv,int *sizeq,int *q,int *iq)
{
    (void)ctx; (void)n2;        /* unused: the signature is shared */
    int sq,tsz;
       
    /* Insert element into heap. */

    *sizeq += 1;    
    q[*sizeq] = v;
    iq[*sizeq] = iv;

    /* Update heap to proper order. */

    sq = *sizeq;
    v = q[sq];
    iv = iq[sq];

L10:
    tsz = sq / 2;
    if (tsz != 0) {    
        if (q[tsz] > v) {     
            q[sq] = q[tsz];
            iq[sq] = iq[tsz];
            sq = tsz;
            goto L10;
        }        
    }
    q[sq] = v;
    iq[sq] = iv;
}

/* ------------------------------------------------------------------------ */    
/*  a_removq(n2,v,iv,sizeq,q,iq)                                            */
/*                                                                          */
/*  Remove smallest element (v,iv) from a priority                          */
/*  queue (q,iq).                                                           */
/*                                                                          */
/*  Passed input scalar:                                                    */
/*  n2     - n * n                                                          */
/*                                                                          */
/*  Passed input/output scalar:                                             */
/*  sizeq  - size of heap                                                   */
/*                                                                          */
/*  Passed output scalars:                                                  */
/*  v      - smallest element in heap (value)                               */
/*  iv     - smallest element in heap (index)                               */
/*                                                                          */
/*  Passed input/output arrays:                                             */
/*  q      - heap (value)                                                   */
/*  iq     - heap (index)                                                   */
/*                                                                          */

void a_removq(TDAContext *ctx, int n2,int *v,int *iv,int *sizeq,int *q,int *iq)
{
    (void)ctx; (void)n2;        /* unused: the signature is shared */
    int vtmp,ivtmp,k,j,szqd2;

    /* Remove element from heap. */

    *v = q[1];
    *iv = iq[1];
    q[1] = q[*sizeq];
    iq[1] = iq[*sizeq];
    *sizeq -= 1;     

    /* Update heap to proper order. */

    k = 1;
    vtmp = q[k];
    ivtmp = iq[k];
    szqd2 = *sizeq / 2;

LL10:
    if (k <= szqd2) {     
        j = k + k;
        if (j < *sizeq) {     
            if (q[j] > q[j + 1]) 
                j++;    
        }
        if (vtmp > q[j]) {     
            q[k] = q[j];
            iq[k] = iq[j];
            k = j;
            goto LL10;
        }         
    }        
    q[k] = vtmp;
    iq[k] = ivtmp;
}

/* ------------------------------------------------------------------------ */    
/*  (double)a_randp(ix)                                                     */
/*                                                                          */
/*  Portable pseudo-random number generator.                                */
/*  Reference: L. Schrage, "A More Portable Fortran                         */
/*  Random Number Generator", ACM Transactions on                           */
/*  Mathematical Software, Vol. 2, No. 2, (June, 1979).                     */
/*                                                                          */

double a_randp(TDAContext *ctx, int *ix)
{
    (void)ctx;        /* unused: the signature is shared */
    int xhi,xalo,leftlo,fhi,k;
    double randp;

    static int a = 16807;
    static int b15 = 32768;
    static int b16 = 65536;
    static int p = 2147483647;

    xhi = *ix / b16;
    xalo = (*ix - xhi * b16) * a;
    leftlo = xalo / b16;
    fhi = xhi * a + leftlo;
    k = fhi / b15;
    *ix = (((xalo - leftlo * b16) - p) + (fhi - k * b15) * b16) + k;
    if (*ix < 0)
        *ix += p;
    randp=(double)(*ix) * 4.656612875e-10;
    return(randp);
}

/* ------------------------------------------------------------------------ */    
/*  a_evalij(n,n2,i,j,xgain,f,d,a)                                          */
/*                                                                          */
/*  Computes the gain in objective function by switching                    */
/*  the locations of facilities i and j (i < j).                            */
/*                                                                          */
/*  Passed input scalars:                                                   */
/*  n      - QAP dimension                                                  */
/*  n2     - n * n                                                          */
/*  i      - permutation array index                                        */
/*  j      - permutation array index                                        */
/*                                                                          */
/*  Passed output scalar:                                                   */
/*  xgain  - gain achieved by swapping i and j in                           */
/*           permutation                                                    */
/*                                                                          */
/*  Passed input arrays:                                                    */
/*  f      - flow matrix                                                    */
/*  d      - distance matrix                                                */
/*  a      - permutation vector                                             */
/*                                                                          */ 

void a_evalij(TDAContext *ctx, int n,int n2,int i,int j,int* xgain,int *f,int *d,int *a)
{
    (void)ctx; (void)n2;        /* unused: the signature is shared */
    int k,aim1tn,ajm1tn,akm1tn,ai,aj,ak,im1tn,jm1tn,km1tn,dtmp1,dtmp2;

    *xgain = 0;
    ai = a[i];
    aj = a[j];
    aim1tn = (ai - 1) * n;
    ajm1tn = (aj - 1) * n;
    im1tn = (i - 1) * n;
    jm1tn = (j - 1) * n;

    km1tn = 0;
    for (k = 1; k <= n; ++k) {
        if (k != i && k != j) {       
            ak = a[k];
            akm1tn = (ak - 1) * n;
            dtmp1 = d[km1tn + i] - d[km1tn + j];
            dtmp2 = d[im1tn + k] - d[jm1tn + k];
            *xgain += dtmp1 * (f[akm1tn + ai] - f[akm1tn + aj]) +
                      dtmp2 * (f[aim1tn + ak] - f[ajm1tn + ak]);
        }        
        km1tn += n;
    }
    dtmp1 = d[im1tn + j] - d[jm1tn + i];
    *xgain += dtmp1 * (f[aim1tn + aj] - f[ajm1tn + ai]);
}

/* ------------------------------------------------------------------------ */
/*  gloc        Location with rectilinear distances.                        */
/*                                                                          */
/*              One-dimensional multifacility location (Picard/Ratliff/     */
/*              Cheung, see g_locate below).  The graph's nodes are split   */
/*              by n=: nodes 1..n are the variable points, the remaining    */
/*              nodes are the fixed points in increasing order of their     */
/*              positions on the line.                                      */
/*                                                                          */
/*              gloc(                                                       */
/*                  gn = ...,       g1,g2, graph numbers: g1's values       */
/*                                  among nodes 1..n are the weights        */
/*                                  between variable points, g2's values    */
/*                                  from nodes 1..n to the fixed nodes are  */
/*                                  the variable-fixed weights.  gn=g uses  */
/*                                  one graph for both.  Default 1,1.       */
/*                  n = ...,        number of variable points (required)    */
/*              );                                                          */
/*                                                                          */
/*              Weights must be non-negative integers; a missing edge       */
/*              counts as 0.                                                */
/*                                                                          */
/*              As shipped in TDA 6.4 this command parsed gn= and then      */
/*              solved one data set hard-coded into the source, whatever    */
/*              the input -- left mid-study, like ivreg1 and zreg1.  The    */
/*              solver g_locate() itself was complete; only this reader     */
/*              was missing.                                                */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gloc(TDAContext *ctx)
{
    register int i,j;
    int err,n,m,r,miss;
    int *v,*w,*post;
    double tmp;

    err = -1;
    v = w = post = NULL;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Location problem. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,0,0))       /* get parameters */
        goto GLOCFin;

    if (gdd_check2(ctx, ctx->PMGN,ctx->PMGN1,1))
        goto GLOCFin;
    if (ctx->PMGN1 < 1)
        ctx->PMGN1 = ctx->PMGN;

    n = ctx->PMN;                    /* variable points */
    if (n < 1) {
        printf1(ctx, "Error: n= (number of variable points) is required.\n");
        goto GLOCFin;
    }
    m = ctx->GD_NP - n;              /* fixed points */
    if (m < 1) {
        printf1(ctx, "Error: n=%d leaves no fixed points (graph has %d nodes).\n",
            n,ctx->GD_NP);
        goto GLOCFin;
    }
    printf1(ctx, "Variable points: %d (nodes 1..%d).\n",n,n);
    printf1(ctx, "Fixed points: %d (nodes %d..%d, in increasing order of their\n",
        m,n + 1,ctx->GD_NP);
    printf1(ctx, "positions on the line).\n");

    if (!(v = (int *)calloc((size_t)(n) * (size_t)(n) + 1,sizeof(int)))) {
        p_err(ctx, -2,1);
        goto GLOCFin;
    }
    if (!(w = (int *)calloc((size_t)(n) * (size_t)(m) + 1,sizeof(int)))) {
        p_err(ctx, -2,1);
        goto GLOCFin;
    }
    if (!(post = (int *)calloc((size_t)(n + 1),sizeof(int)))) {
        p_err(ctx, -2,1);
        goto GLOCFin;
    }

    miss = 0;
    for (i = 0; i < n; ++i) {         /* weights between variable points */
        for (j = 0; j < n; ++j) {
            tmp = gdd_adj(ctx, i,j,ctx->PMGN);
            if (tmp < 0.0) {
                tmp = 0.0;
                if (i != j)
                    miss++;
            }
            v[i * n + j + 1] = (int)tmp;
        }
    }
    for (i = 0; i < n; ++i) {         /* weights variable -> fixed */
        for (j = 0; j < m; ++j) {
            tmp = gdd_adj(ctx, i,n + j,ctx->PMGN1);
            if (tmp < 0.0) {
                tmp = 0.0;
                miss++;
            }
            w[i * m + j + 1] = (int)tmp;
        }
    }
    if (miss > 0)
        printf1(ctx, "Substituted %d missing values.\n",miss);

    r = g_locate(ctx, n,m,v,w,post);
    if (r) {
        printf1(ctx, "Error: location algorithm failed.\n");
        goto GLOCFin;
    }
    newline(ctx);
    printf1(ctx, "Optimal locations:\n");
    for (i = 1; i <= n; ++i) {
        printf1(ctx, "Variable point %d -> fixed point %d (node %d)\n",
            i,post[i],n + post[i]);
#ifdef TDA_R_PACKAGE
        /* variable point, fixed point, node -- beside the print */
        tda_export_cell(ctx, "gloc.assignment", (double)i);
        tda_export_cell(ctx, "gloc.assignment", (double)post[i]);
        tda_export_cell(ctx, "gloc.assignment", (double)(n + post[i]));
        tda_export_endrow(ctx, "gloc.assignment");
#endif
    }

    err = 0;

GLOCFin:
    if (v != NULL)
        free((char *)v);
    if (w != NULL)
        free((char *)w);
    if (post != NULL)
        free((char *)post);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */ 
/*  g_locate(nvpt,nfpt,vwt,vfwt,post)                                       */
/*                                                                          */
/*  THIS SUBROUTINE SOLVES THE FOLLOWING ONE-DIMENSIONAL                    */
/*  MULTIFACILITY LOCATION PROBLEM USING AN ALGORITHM DESIGNED              */
/*  BY PICARD, RATLIFF (SEE OPERATIONS RESEARCH 26 (1978),                  */
/*  PP.422-433) AND CHEUNG.                                                 */
/*    " GIVEN NFPT FIXED POINTS ON A LINE, DETERMINE THE                    */
/*      LOCATIONS OF NVPT VARIABLE POINTS SUCH THAT A WEIGHTED              */
/*      SUM OF THE DISTANCES BETWEEN THE POINTS IS MINIMUM."                */
/*  OUR ALGORITHM TRANSFORMS THE LOCATION PROBLEM TO A SEQUENCE             */
/*  OF MINIMUM-CUT PROBLEMS. A FORTRAN CODE (THE SUBROUTINE                 */
/*  NETFLO IN "COMBINATORIAL ALGORITHMS", WRITTEN BY NIJENHUIS              */
/*  AND WILF) OF DINIC'S ALGORITHM IS USED TO FIND THE MINIMUM              */
/*  CUTS.                                                                   */
/*  TO APPLY THE SUBROUTINE  'LOCATE', THE FIXED POINTS SHOULD              */
/*  BE SUBSCRIPTED ON THE LINE IN INCREASING ORDER OF THEIR                 */
/*  LOCATIONS. THEIR EXACT LOCATIONS ARE NOT REQUIRED.                      */
/*                                                                          */
/*  nvpt = number of variable points                                        */
/*  nfpt = number of fixed points                                           */
/*  vwt[i,j] = (nvpt,nvpt) square matrix, non-negative integers. The        */
/*      upper triangle contains the weights between the var. points.        */
/*  vfwt[i,j] = (nvpt,nfpt) matrix, non-negative integers with weights      */
/*      between fixed and variable points.                                  */
/*                                                                          */
/*  post[i], i = 1,...,nvpt, at output: indicating the optimal locations    */
/*      of the variable points. Var point i will coincide with the fixed    */
/*      point whose position is post[i].                                    */
/*                                                                          */
/*  working storage:                                                        */
/*  endpt[i], i = 1,...,4 * e, e = nvpt * (nvpt + 1).                       */
/*  cut[i], i = 1,...,nvpt + 2                                              */
/*                                                                          */
/*  Return 0 if OK, -1 if error (insufficient memory)                       */

int g_locate(TDAContext *ctx, int nvpt,int nfpt,int *vwt,int *vfwt,int *post)     
{
    register int i,j,k,l;
    int col,source,sum,nn,n2,nnode,il,jl,ea,ec,ex,iter;
    int *endp1,*endp2,*endp3,*endp4,*cap,*aux;
    int cut[1000],cvpt[1000];

    if (nvpt == 1) {

        /* IN CASE THE NUMBER OF VARIABLE POINTS IS 1, OUR PROBLEM IS
           THE SAME AS THE 1-MEDIAN LOCATION PROBLEM. OUR ALGORITHM
           IS ALSO EQUIVALENT TO THE FOLLOWING WELL-KNOWN SIMPLE
           PROCESS FOR SOLVING 1-MEDIAN LOCATION PROBLEMS. */

        sum = -vfwt[1];   
        for (i = 2; i <= nfpt; ++i)
            sum += vfwt[i];  

        i = 1;
L15:
        if (sum <= 0)
            goto L16;
        i++;   
        sum -=  2 * vfwt[i];     
        goto L15;

L16:
        post[1] = i;
        return(0);
    }         

    ea = nvpt * (nvpt + 1) + 1;
    ec = (nvpt + 2) * (nvpt + 2) + 1;
    ex = nvpt + 3;

    if (!(endp1 = (int *)calloc((size_t)(ea),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);
    }
    if (!(endp2 = (int *)calloc((size_t)(ea),sizeof(int)))) {
        free((char *)endp1);
        p_err(ctx, -2,1);
        return(-1);
    }
    if (!(endp3 = (int *)calloc((size_t)(ea),sizeof(int)))) {
        free((char *)endp1);
        free((char *)endp2);
        p_err(ctx, -2,1);
        return(-1);
    }
    if (!(endp4 = (int *)calloc((size_t)(ea),sizeof(int)))) {
        free((char *)endp1);
        free((char *)endp2);
        free((char *)endp3);
        p_err(ctx, -2,1);
        return(-1);
    }
    if (!(cap = (int *)calloc((size_t)(ec),sizeof(int)))) {
        free((char *)endp1);
        free((char *)endp2);
        free((char *)endp3);
        free((char *)endp4);
        p_err(ctx, -2,1);
        return(-1);
    }
    if (!(aux = (int *)calloc((size_t)(ex),sizeof(int)))) {
        free((char *)endp1);
        free((char *)endp2);
        free((char *)endp3);
        free((char *)endp4);
        free((char *)cap);
        p_err(ctx, -2,1);
        return(-1);
    }
    memrq(ctx, 4 * ea + ec + ex,sizeof(double));
 
    /* initialize endp for the first iteration. */

    source = nvpt + 1;
    col = 0;

    for (i = 1; i <= nvpt; ++i) {

        cvpt[i] = i;
        post[i] = 0;

        for (j = 1; j <= nvpt; ++j) {

            col++;    
            endp1[col] = i;
            endp2[col] = j;
            endp4[col] = 0;

            if (i < j)  
                endp3[col] = vwt[(i-1) * nvpt + j];
                
            else if (i == j) {
                endp1[col] = source;
                sum = 0;
                for (k = 1; k <= nfpt; ++k)  
                    sum += vfwt[(i-1) * nfpt + k];

                endp3[col] = sum;
            }
            else  
                 endp3[col] = vwt[(j-1) * nvpt + i];
        }
    }
    nn = nvpt;

#if 0   /* leftover debug output, disabled -- see README */
tda_out("nn=%d col=%d\n",nn,col);
#endif
 
    /* iterations */

    for (iter = 1; iter <= nfpt; ++iter) {
       
        nnode = nn + 2;
#if 0   /* leftover debug output, disabled -- see README */
tda_out("iter=%d nnode=%d nn=%d\n",iter,nnode,nn);
#endif
        for (i = 1; i <= nn; ++i) {
            col++;     
            endp1[col] = i;
            endp2[col] = nnode;
            l = cvpt[i];
            endp3[col] = 2 * vfwt[(l-1) * nfpt + iter];
            endp4[col] = 0;
        }

#if 0   /* leftover debug output, disabled -- see README */
tda_out("vor netflo nnode=%d col=%d\n",nnode,col);
#endif


        /* solve a min-cut problem by Dinic's algorithm. */

        g_netflo(ctx, nnode,col,endp1,endp2,endp3,endp4,source,nnode,cut,cap,aux);
 
        /* update post and initialize cvpt for the next iteration. */

        n2 = nn;
        nn = 0;
        for (i = 1; i <= n2; ++i) {
            if (cut[i] == 1) {
                nn++;     
                cvpt[nn] = cvpt[i];
            }
            else {         
                l = cvpt[i];
                post[l] = iter;
            }       
        }
 
        /* nn is the number of variable points in the next iteration. */

        if (nn == 0)  
            break;

        /* initialize endp for the next iteration. */

        source = nn + 1;
        col = il = 0;

        for (i = 1; i <= n2; ++i) {
            if (cut[i] == 0)          
                continue;
            il++;    
            jl = 0;
            for (j = 1; j <= n2; ++j) {
                if (cut[j] == 0)              
                    continue;
                col++;     
                jl++;    
                endp1[col] = il;
                if (i == j)
                    endp1[col] = source;
                endp2[col] = jl;
                l = (i-1) * n2 + j;
                endp3[col] = endp3[l] - endp4[l];
                endp4[col] = 0;
            }
        }
    }
    free((char *)endp1);
    free((char *)endp2);
    free((char *)endp3);
    free((char *)endp4);
    free((char *)cap);
    free((char *)aux);
    memrq(ctx, -4 * ea - ec - ex,sizeof(double));
    return(0);
}         

/* ------------------------------------------------------------------------ */ 
/*  g_netflo(n,e,endp1,endp2,endp3,endp4,source,sink,cut)                   */
/*                                                                          */
/*  THIS SUBROUTINE SOLVES A MAX-FLOW-MIN-CUT PROBLEM USING                 */
/*  DINIC'S ALGORITHM. THE CODE IS TAKEN FROM "COMBINATORIAL                */
/*  ALGORITHMS", ACADEMIC PRESS, WRITTEN BY A. NIJENHUIS                    */
/*  AND H.S. WILF. THERE ARE THREE MINOR CHANGES: (1) STATEMENT             */
/*  LABELS ARE RENUMBERED. (2) SOME DATA DECLARATION STATEMENTS             */
/*  AND THE FORMAL PARAMETER LIST HAVE BEEN MODIFIED. (3) THE               */
/*  WAY OF SUBSCRIPT CONSTRUCTION FOR SOME OF THE ARRAYS IS                 */
/*  MODIFIED. THE MODIFICATIONS ARE MAINLY FOR COMPLIANCE WITH              */
/*  ANSI FORTRAN STANDARD.                                                  */
/*  THANKS ARE DUE TO THE AUTHORS AND PUBLISHER OF                          */
/*  'COMBINATORIAL ALGORITHMS' FOR THEIR PERMISSION TO INCLUDE              */
/*  THIS SUBROUTINE.                                                        */
/*                                                                          */
/*  n = number of nodes                                                     */
/*  e = number of edges                                                     */
/*  endp1[i], i = 1,...,e       starting point of edge                      */
/*  endp2[i], i = 1,...,e       ending point of edge                        */
/*  endp3[i], i = 1,...,e       capacity of edge                            */
/*  endp4[i], i = 1,...,e       max flow through edge                       */
/*  cap[i,j], (n+2,n+2) working array                                       */
/*  aux[i], n + 2 working array                                             */

void g_netflo(TDAContext *ctx, int n,int e,int *endp1,int *endp2,int *endp3,int *endp4, int source,int sink,int *cut,int *cap,int *aux)
{
    register int i,j,k;
    int k0 = 0,m,c,p,q,rd,wr,delta,i1,i2,nmin,lblsnk,label;
    int vert[1000];


    /* initialization. */

    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) 
/*          CAP(I,J) = 0;           */
            cap[(i-1) * n + j] = 0;

    }
    for (i = 1; i <= e; ++i) {
        i1 = endp1[i];
        i2 = endp2[i];
/*      CAP(i1,i2) = endp3[i];      */
        cap[(i1 - 1) * n + i2] = endp3[i];

    }
    for (i = 1; i <= n; ++i) {
        k = 0;
        for (j = 1; j <= n; ++j) {
/*          if (CAP(i,j) + CAP(j,i) == 0)           */
            if (cap[(i-1) * n + j] + cap[(j-1) * n + i] == 0)
                continue;
            k++;   
/*          VERT(i,k) = j;      */
            vert[(i-1) * n + k] = j;
        }
/*      VERT(i,n) = k;      */

        vert[(i-1) * n + n] = k;
    }
    nmin = -n - 1;
 
    /* scanning and labeling */

L6:
    lblsnk = n;
    for (i = 1; i <= n; ++i)
        cut[i] = nmin;

    rd = wr = 0;
    p = source;
    label = -1;
L8:
/*  m = VERT(p,n);      */

    m = vert[(p-1) * n + n];

    i = 1;
L9:
    if (i > m)
        goto L13;

/*  q = VERT(p,i);      */

    q = vert[(p-1) * n + i];

/*  if(CAP(p,q) == 0)   */

    if (cap[(p - 1) * n + q] == 0)
        goto L12;

    if (q == sink)
        lblsnk = -label;

    if (cut[q] <= label) {
        if (cut[q] < label) {
            cut[q] = label;
            wr++;   
            aux[wr] = q;
        }
        i++;   
        goto L9;
    }
L12:         
/*      VERT(p,i) = VERT(p,m);
        VERT(p,m) = q;          */

        vert[(p-1) * n + i] = vert[(p-1) * n + m];
        vert[(p-1) * n + m] = q;


        m--;   
        goto L9;
       

L13:
    cut[p] = m;
    rd++;    
    if (rd > wr)
        goto L24;
    p = aux[rd];

    if (cut[p] + lblsnk == 0)
        goto L14;

    label = cut[p] - 1;
    goto L8;
 
    /* CONSTRUCTION OF PATH FROM SOURCE TO SINK */

L14:
    q = source;
    k = 0;
L15:
    k++;  
    aux[k] = q;
    if (k > lblsnk)
        goto L19;
L16:
    p = aux[k];
L17:
    m = cut[p];
    if (m == 0)
        goto L18;
/*  q = VERT(p,m);          */

    q = vert[(p-1) * n + m];


    goto L15;
L18:
    k--;  
    if (k == 0)
        goto L6;
    p = aux[k];
    cut[p] -= 1;      
    goto L17;
L19:
    if (q != sink)
        goto L18;
/*  delta = CAP(p,q);   */
    delta = cap[(p - 1) * n + q];

    for (i = 2; i <= k; ++i) {
        i1 = aux[i-1];
        i2 = aux[i];
/*      delta = imin(delta,CAP(I1,I2));     */

        delta = imin(ctx, delta,cap[(i1-1) * n + i2]);

    }
L21:
    k--;  
    if (k == 0)
        goto L23;
    p = aux[k];
/*  c = CAP(p,q) - delta;   */
    c = cap[(p - 1) * n + q] - delta;

    if (c > 0)
        goto L22;
    cut[p] -= 1;    
    k0 = k;
L22:
/*  CAP(P,Q) = c;
    CAP(Q,P) += delta;      */

    cap[(p - 1) * n + q] = c;
    cap[(q - 1) * n + p] += delta;

    q = p;
    goto L21;
L23:
    k = k0;
    goto L16;
 
    /* exit procedure. */

L24:
    for (i = 1; i <= n; ++i)  
        cut[i] = imin(ctx, 1,cut[i] - nmin);

    for (i = 1; i <= e; ++i) {
        i1 = endp1[i];
        i2 = endp2[i];
/*      endp4[i] = endp3[i] - CAP(I1,I2);       */
        endp4[i] = endp3[i] - cap[(i1 - 1) * n + i2];
    }
}

/* -####-------------------------------------------------------------------- */ 
/*  qap_w(n,c,f,d,loc3n,wtmp)                                               */  
/*                                                                          */  
/*  Approximate solution of the Quadratic Assignment Problem. Adapted       */
/*  from CACM algorithm 608 (D.H. West).                                    */
/*                                                                          */
/*  c, f, and d are (n,n)-matrices, n > 1. The main diagonal of f and d     */
/*  must contain zeros. The procedure tries to find a permutation p that    */
/*  minimizes                                                               */
/*                                                                          */
/*  sum_i ( c[i,p(i)] + sum_j f[i,j] d[p(i),p(j()] )                        */
/*                                                                          */
/*  loc3n must have length at least 3 * n, wtmp must have length at         */
/*  least 4 * n.                                                            */
/*                                                                          */
/*  The final function value is returned in wtmp[1], the permutation is     */
/*  return in loc3n[1,...,n].                                               */
/*                                                                          */
/*  Return 0 if OK, -1 if error (pnorm <= EPSI1 or n <= 1)                  */

int qap_w(TDAContext *ctx, int n,double *c,double *f,double *d,int *loc3n,double *wtmp)
{
    register int im,kp1,i,j,m;
    int n1,n2,k,is,ibar,jbar,jl,ibest = 0,jbest = 0,kp1s,kp1l,ibests,jbestl;
    int itemp,jtemp,il,ls1,kp11,lkp1,li,lj,lm,m1,lm1;
    int *ls,*ll;
    double tmp,obj,delbst,pnorm,scsl,sfs,sdl,cmx,fmx,dmx,fr,fc,dr,dc;
    double f1,f2,d1,d2,db1,dbp2,dbp5,rfnk,rfnk1,rfnk2;
    double bp[8],bpb[8],bpn[8],*w1,*w2,*w3,*w4;
 
    if (n <= 1) {
        printf1(ctx, "Error in qap_w: dimension is less than 2\n");
        return(-1);                             
    }
    n1 = n;
    n2 = n + n;

    ls = loc3n + n1;
    ll = loc3n + n2;
    w1 = wtmp;
    w2 = w1 + n;
    w3 = w2 + n;
    w4 = w3 + n;

    /* initialization */

    scsl = sfs = sdl = cmx = fmx = dmx = 0.0;
    rfnk = 1.0 / (double)n1; 
    rfnk1 = 0.0;
    rfnk1 = 1.0 / (double)(n1 - 1);
    rfnk2 = 0.0;
    if (n1 > 2)
        rfnk2 = 1.0 / (double)(n1 - 2);
 
    for (i = 1; i <= n1; ++i) {

        ls[i] = i;
        ll[i] = i;
        fr = fc = dr = dc = 0.0;
 
        for (j = 1; j <= n1; ++j) {
            tmp  = c[(i - 1) * n + j];
            scsl += tmp;                   
            cmx = dmax(ctx, cmx,fabs(tmp));
            tmp = f[(i - 1) * n + j];
            fr += tmp;                     
            fmx = dmax(ctx, fmx,fabs(tmp));
            fc += f[(j - 1) * n + i];
            tmp = d[(i - 1) * n + j];
            dr += tmp;                
            dmx = dmax(ctx, dmx,fabs(tmp));
            dc += d[(j - 1) * n + i];
        }
        sfs += fr;         
        sdl += dr;          
        w1[i] = fr;
        w3[i] = fc;
        w2[i] = dr;
        w4[i] = dc;
    }
    bp[1] = 0.0;
    bp[2] = 0.0;
    bp[3] = sfs * sdl;
    bp[4] = 0.0;
    bp[5] = scsl;
    bp[6] = sfs;
    bp[7] = sdl;
    pnorm = cmx + fmx * dmx;

    if (pnorm <= ctx->EPSI1) {
        printf1(ctx, "Error in qap_w: matrices are zero.\n");
        return(-1);
    }

    /* end of initialization */

    obj = (bp[3] * rfnk1 + bp[5]) * rfnk;
 
    for (kp1 = 1; kp1 < n; ++kp1) {

        k = kp1 - 1;
        delbst = pnorm;

        for (i = kp1; i <= n; ++i) {

            is = n1 + i;
            ibar = loc3n[is];
 
            for (j = kp1; j <= n; ++j) {

                jl = n2 + j;
                jbar = loc3n[jl];

                db1 = 0.0;
                dbp2 = w1[ibar] * w2[jbar] + w3[ibar] * w4[jbar];

                if (k > 0) {
                    for (im = 1; im <= k; ++im) {
                        m1 = ls[im];
                        lm1 = ll[im];
                        f1 = f[(ibar - 1) * n + m1];
                        f2 = f[(m1 - 1) * n + ibar];
                        d1 = d[(jbar - 1) * n + lm1];
                        d2 = d[(lm1 - 1) * n + jbar];
                        db1 = db1 + f1 * d1 + f2 * d2;
                        dbp2 = dbp2 - f2 * w2[lm1] - f1 * w4[lm1] -
                                           w1[m1] * d2 - w3[m1] * d1;
                    }
                    dbp2 += db1;
                }
                bpn[1] = bp[1] + db1;
                bpn[2] = bp[2] + dbp2;
                bpn[6] = bp[6] - w1[ibar] - w3[ibar];
                bpn[7] = bp[7] - w2[jbar] - w4[jbar];
                bpn[3] = bpn[6] * bpn[7];

                if (k >= n1 - 2)
                    bpn[3] = 0.0;

                bpn[4] = bp[4] + c[(ibar - 1) * n + jbar];
                dbp5 = c[(ibar - 1) * n + jbar];
 
                for (im = kp1; im <= n1; ++im) {
                    m1 = ls[im];
                    lm1 = ll[im];
                    dbp5 = dbp5 - c[(ibar - 1) * n + lm1] - c[(m1 - 1) * n + jbar];
                }
                bpn[5] = bp[5] + dbp5;
                tmp = db1 + c[(ibar - 1) * n + jbar] - rfnk * (bp[2] + bp[5]) +
                       rfnk1 * ((bpn[3] * rfnk2 - bp[3] * rfnk) + (bpn[2] + bpn[5]));

                if (tmp < delbst) {
                    ibest = i;
                    jbest = j;
                    delbst = tmp;
                    bpb[1] = bpn[1];
                    bpb[2] = bpn[2];
                    bpb[3] = bpn[3];
                    bpb[4] = bpn[4];
                    bpb[5] = bpn[5];
                    bpb[6] = bpn[6];
                    bpb[7] = bpn[7];
                }
            }
        }
        kp1s = n1 + kp1;
        kp1l = n2 + kp1;
        ibests = n1 + ibest;
        jbestl = n2 + jbest;
        itemp = loc3n[kp1s];
        loc3n[kp1s] = loc3n[ibests];
        loc3n[ibests] = itemp;

        jtemp = loc3n[kp1l];
        loc3n[kp1l] = loc3n[jbestl];
        loc3n[jbestl] = jtemp;
        obj += delbst;

        if (kp1 != n - 1) {

            kp11 = ls[kp1];
            lkp1 = ll[kp1];
 
            for (im = 1; im <= n1; ++im) {
                w1[im] -= f[(im - 1) * n + kp11];
                w3[im] -= f[(kp11 - 1) * n + im];
                w2[im] -= d[(im - 1) * n + lkp1];
                w4[im] -= d[(lkp1 - 1) * n + im];
            }
            bp[1] = bpb[1];
            bp[2] = bpb[2];
            bp[3] = bpb[3];
            bp[4] = bpb[4];
            bp[5] = bpb[5];
            bp[6] = bpb[6];
            bp[7] = bpb[7];

            rfnk = rfnk1;
            rfnk1 = rfnk2;
            rfnk2 = 0.0;

            if (kp1 < n1 - 2)
                rfnk2 = 1.0 / (double)(n1 - kp1 - 2);
        }
    }
    for (i = 1; i <= n; ++i) {
        il = n2 + i;
        is = n1 + i;
        ls1 = loc3n[is];
        loc3n[ls1] = loc3n[il];
    }
    for (m = 1; m <= n; ++m) {

        delbst = 0.0;
 
        for (i = 2; i <= n; ++i) {

            for (j = 1; j < i; ++j) {

                li = loc3n[i];
                lj = loc3n[j];
                tmp = c[(i - 1) * n + lj] - c[(i - 1) * n + li] + 
                      c[(j - 1) * n + li] - c[(j - 1) * n + lj] +
                      (f[(i - 1) * n + j] - f[(j - 1) * n + i]) *
                      (d[(lj - 1) * n + li] - d[(li - 1) * n + lj]);
 
                for (im = 1; im <= n1; ++im) {
                    if (im == i || im == j)
                        continue;

                    lm = loc3n[im];
                    tmp = tmp + (f[(i - 1) * n + im] -
                           f[(j - 1) * n + im]) *
                           (d[(lj - 1) * n + lm] - d[(li - 1) * n + lm]) +
                           (f[(im - 1) * n + i] - f[(im - 1) * n + j]) *
                           (d[(lm - 1) * n + lj] - d[(lm - 1) * n + li]);
                }
                if (tmp < delbst) {
                    ibest = i;
                    jbest = j;
                    delbst = tmp;
                }
            }
        }
        if (delbst >= 0.0)
/*          break;      */
            continue;

        li = loc3n[ibest];
        lj = loc3n[jbest];
        loc3n[ibest] = lj;
        loc3n[jbest] = li;
        obj += delbst;
    }
    w1[1] = obj;
    return(0);
}             










