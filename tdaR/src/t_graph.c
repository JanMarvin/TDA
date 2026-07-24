/****************************************************************************/
/*  t_graph                                                                 */
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
#include "t_ml.h"
#include "t_rand.h"
#include "t_gio.h"
#include "t_svd.h"
#include "t_com.h"
#include "tda_context.h"

/*  functions in t_graph.c */

int gni(TDAContext *ctx);
int gdln(TDAContext *ctx);
int gdp(TDAContext *ctx);
void p_gdp(TDAContext *ctx, int i,int j,int ii,int jj);
int gda(TDAContext *ctx);
int gdu(TDAContext *ctx);
int gde(TDAContext *ctx);
int gdcset(TDAContext *ctx);
int f_gdcset(TDAContext *ctx, int m,int *com,int *com1,int *nci,int *ncj);
int gsym(TDAContext *ctx);
int gdot(TDAContext *ctx);
int gtcl(TDAContext *ctx);
void g_tclos(TDAContext *ctx, int gn,int n,char *ai);
void g_spath(TDAContext *ctx, int gn,int n,float *a);
int gsort(TDAContext *ctx);
int gdcyc(TDAContext *ctx);
int g_dcyc(TDAContext *ctx, int gn,int n,int *p,int **h,int *hp,int opt,int ns,int cmax, int *minclen,int *maxclen,int dopt);
int gev(TDAContext *ctx); 
void g_op1(TDAContext *ctx, int n,double *z,double *w);
void g_evcheck(TDAContext *ctx, int r,int n,double *x,double *d);
int gsp(TDAContext *ctx);
int g_glist(TDAContext *ctx, int i,int *nodes);
void g_sp(TDAContext *ctx, int in,int n,int gn,double *d,int *ndec);
int tnet(TDAContext *ctx);


/* ------------------------------------------------------------------------ */
/*  gni()       Information about nodes.                                    */
/*                                                                          */
/*              gni(                                                        */
/*                  gn = ...,   graph number, def 1 or multigraph           */
/*                  nfmt = ..., integer format, def. 4                      */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gni(TDAContext *ctx)
{
    register int i,k;
    int err,n,nrec;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Degree of nodes. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto GNIFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto GNIFin;

    if (alloc_acr(ctx, ctx->GD_NG + 1))
        goto GNIFin;
    if (alloc_acs(ctx, ctx->GD_NG + 1))
        goto GNIFin;
    if (alloc_acn(ctx, ctx->GD_NG + 1))
        goto GNIFin;

    nrec = 0;
    for (i = 0; i < ctx->GD_NP; ++i) {
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);                
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));                
  
        for (k = 1; k <= ctx->GD_NG; ++k) {
            if (ctx->PMGN != 0 && k != ctx->PMGN)
                continue;
   
            n = g_pre(ctx, i,k);   
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,n);                
            if (ctx->AcR[k] < n)
                ctx->AcR[k] = n;

            if (ctx->GD_GT > 2)
                n = g_suc(ctx, i,k);   
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,n);                
            if (ctx->AcS[k] < n)
                ctx->AcS[k] = n;

            n = g_loop(ctx, i,k);   
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,n);                
            if (n > 0)
                ctx->AcN[k] += 1;
        }
        fprintf(ctx->PMFd,"\n");                
        nrec++;
    }
    printf1(ctx, "\nGraph  max in-degree  max out-degree  nodes with loops\n");
    for (k = 1; k <= ctx->GD_NG; ++k)  
        printf1(ctx, "%5d  %8d       %8d       %14d\n",k,ctx->AcR[k],ctx->AcS[k],ctx->AcN[k]);

    printf1(ctx, "\n%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

GNIFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdln()      Forward and backward links                                  */
/*              Requires a directed graph with GD_TYP = 1.                  */
/*                                                                          */
/*              gdln(                                                       */
/*                  opt = ...,      1 = forward links                       */
/*                                  2 = backward links                      */
/*                  nfmt=...,       integer format, def. 4                  */
/*                  gn  = graph number, def. 1                              */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gdln(TDAContext *ctx)
{
    register int i,j,k,l;
    int err,nrec,n,maxl;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Forward/backward links. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto GFLFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GFLFin;
    if (gdd_tcheck(ctx, 1,2,0))
        goto GFLFin;

    if (ctx->PMOPT == 1)
        printf1(ctx, "Option 1: forward links.\n");
    else {
        printf1(ctx, "Option 2: backward links.\n");
        ctx->PMOPT = 2;
    }
    if (alloc_acn(ctx, imax(ctx, ctx->GD_FPMax,ctx->GD_BPMax) + 2))
        goto GFLFin;

    maxl = nrec = 0;
    for (i = 0; i < ctx->GD_NP; ++i) {
        if (ctx->PMOPT == 1)
            n = ctx->GD_FPN[i];
        else
            n = ctx->GD_BPN[i];
        if (n <= 0)
            continue;

        l = 0;
        for (j = 0; j < n; ++j) {

            if (ctx->PMOPT == 1)
                k = ctx->GD_FPK[i][j];
            else
                k = ctx->GD_BPK[i][j];

            if (k >= 0 && gdd_ev(ctx, k,ctx->PMGN) >= 0.0) {
                if (ctx->PMOPT == 1)
                    ctx->AcN[l++] = gdd_node(ctx, ctx->GD_FPI[i][j]);        
                else
                    ctx->AcN[l++] = gdd_node(ctx, ctx->GD_BPI[i][j]);        
            }
        }
        if (l > 0) {
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);        
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));        
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,l);        
            for (j = 0; j < l; ++j)
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcN[j]);        
            fprintf(ctx->PMFd,"\n");                
            nrec++;
            if (maxl < l)   
                maxl = l;
        }
    }
    printf1(ctx, "Max number of links: %d\n",maxl);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

GFLFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdp()       Writing graph data                                          */
/*                                                                          */
/*              gdp(                                                        */
/*                  opt = ...,      1 = edge list                           */
/*                                  2 = lower triangle                      */
/*                                  3 = lower triangle with diagonal        */
/*                                  4 = full adj. matrix                    */
/*                                  5 = adj. matrix as square matrix        */
/*                  gn = ...,       gdd graph number, required for opt = 5  */
/*                                  default gn = 1                          */
/*                  sc=...,         substitute for missing values, def. -1  */
/*                  nfmt =...,      print format for integers, def. 4       */
/*                  fmt = ...,      print format for values, def. 10.4      */
/*                  if = ...,       input file with node numbers            */
/*              ) = fname;          output file name, required.             */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gdp(TDAContext *ctx)
{
    register int i,j,k;
    int err,nrec,m,ni,ii,jj;
    double tmp;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Writing graph data. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto GDPFin;
   
    if (gdd_check(ctx, ctx->PMGN,1))
        goto GDPFin;

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    if (ctx->PMOPT > 5)
        ctx->PMOPT = 1;

    if (ctx->PMSCFlg == 0)
        ctx->PMSC = -1.0;

    if (ctx->PMOPT >= 5 && ctx->PMGN == 0)  
        ctx->PMGN = 1;

    if (alloc_acn(ctx, ctx->GD_NP))                 
        goto GDPFin;

    ni = 0;
    if (ctx->PMF2Def) {
        if (alloc_ack(ctx, ctx->GD_NP))                 
            goto GDPFin;

        if ((ni = g_getnd(ctx, 0,ctx->AcK)) < 0)         /* get node list */
            goto GDPFin;

        printf1(ctx, "Number of input nodes: %d\n",ni);
        if (ni == 0)
            goto GDPFin;
    }
    if (alloc_actmp(ctx, ctx->GD_NG + 1))                 
        goto GDPFin;
       
    nrec = 0;

    if (ctx->PMOPT == 1) {
        for (i = 0; i < ctx->GD_NP; ++i) {
            n_info(ctx, i,100);
            if (ni > 0 && ctx->AcK[i] == 0)
                continue;
            ii = gdd_node(ctx, i);

            for (j = 0; j < ctx->GD_NP; ++j) {
                if (ni > 0 && ctx->AcK[j] == 0)
                    continue;

                jj = gdd_node(ctx, j);

                if (ctx->PMGN != 0) {
                    tmp = gdd_adj(ctx, i,j,ctx->PMGN);
                    if (tmp >= 0.0) {
                        p_gdp(ctx, i + 1,j + 1,ii,jj);
                        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
                        fprintf(ctx->PMFd,"\n");                
                        nrec++;    

                        ctx->AcN[i] = ctx->AcN[j] = 1;
                    }
                }
                else {
                    m = 0;
                    for (k = 1; k <= ctx->GD_NG; ++k) {
                        if (gdd_adj(ctx, i,j,k) >= 0.0) {
                            m = 1;
                            break;
                        }
                    }
                    if (m) {
                        p_gdp(ctx, i + 1,j + 1,ii,jj);
                        for (k = 1; k <= ctx->GD_NG; ++k) {
                            tmp = gdd_adj(ctx, i,j,k);
                            if (tmp < 0.0)
                                tmp = ctx->PMSC;
                            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
                        }
                        fprintf(ctx->PMFd,"\n");                
                        nrec++;    

                        ctx->AcN[i] = ctx->AcN[j] = 1;
                    }
                }
            }
        }
      
        /* add isolated nodes */

        tmp = ctx->PMSC;
        for (i = 0; i < ctx->GD_NP; ++i) {
            if (ctx->AcN[i] == 0) {
                if (ni > 0 && ctx->AcK[i] == 0)
                    continue;
                ii = gdd_node(ctx, i);

                if (ctx->PMGN != 0) {
                    p_gdp(ctx, i + 1,i + 1,ii,ii);
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
                }
                else {
                    p_gdp(ctx, i + 1,i + 1,ii,ii);
                    for (k = 1; k <= ctx->GD_NG; ++k)  
                        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
                }
                fprintf(ctx->PMFd,"\n");                
                nrec++;    
            }
        }
    }
    else if (ctx->PMOPT == 2 || ctx->PMOPT == 3 || ctx->PMOPT == 4) {
        for (i = 0; i < ctx->GD_NP; ++i) {
            n_info(ctx, i,100);
            if (ni > 0 && ctx->AcK[i] == 0)
                continue;

            ii = gdd_node(ctx, i);

            for (j = 0; j < ctx->GD_NP; ++j) {

                if ((ctx->PMOPT == 2 && j >= i) || (ctx->PMOPT == 3 && j > i))
                    continue;

                if (ni > 0 && ctx->AcK[j] == 0)
                    continue;

                jj = gdd_node(ctx, j);

                if (ctx->PMGN != 0) {
                    p_gdp(ctx, i + 1,j + 1,ii,jj);
                    tmp = gdd_adj(ctx, i,j,ctx->PMGN);
                    if (tmp < 0.0)
                        tmp = ctx->PMSC;
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
                    fprintf(ctx->PMFd,"\n");                
                    nrec++;    
                }
                else {
                    p_gdp(ctx, i + 1,j + 1,ii,jj);
                    for (k = 1; k <= ctx->GD_NG; ++k) {
                        tmp = gdd_adj(ctx, i,j,k);
                        if (tmp < 0.0)
                            tmp = ctx->PMSC;
                        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
                    }
                    fprintf(ctx->PMFd,"\n");                
                    nrec++;    
                }
            }
        }
    }
    else if (ctx->PMOPT == 5) {
        for (i = 0; i < ctx->GD_NP; ++i) {
            n_info(ctx, i,100);
            if (ni > 0 && ctx->AcK[i] == 0)
                continue;

            ii = gdd_node(ctx, i);

            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);                
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ii);                

            for (j = 0; j < ctx->GD_NP; ++j) {

                if (ni > 0 && ctx->AcK[j] == 0)
                    continue;

                tmp = gdd_adj(ctx, i,j,ctx->PMGN);
                if (tmp < 0.0)
                    tmp = ctx->PMSC;
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
            }
            fprintf(ctx->PMFd,"\n");                
            nrec++;    
        }
    }
    n_info_e(ctx);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

GDPFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  p_gdp(i,j,ii,jj)                                                        */

void p_gdp(TDAContext *ctx, int i,int j,int ii,int jj)
{
    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i);                
    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j);                
    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ii);                
    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,jj);                
}

/* ------------------------------------------------------------------------ */
/*  gda()       Aggregating nodes of a graph. Only with first graph!        */
/*                                                                          */
/*              gda(                                                        */
/*                  rcn=...,        recode information                      */
/*                  nfmt =...,      print format for integers, def. 4       */
/*                  fmt = ...,      print format for values, def. 10.4      */
/*              ) = fname;          output file name, required.             */
/*                                                                          */
/*  rcn = n1[k1,k2,...],...                                                 */
/*                                                                          */
/*  rcn = n1[k1,,k2,...],...    Then take all node numbers from k1 to k2.   */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int gda(TDAContext *ctx)
{
    register int i,j,k;
    int err,nrec,n,m,mm,nm,ii,jj,fnd;
    char *p;
    double tmp;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Aggregating nodes of a graph. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto GDAFin;
   
    if (gdd_check(ctx, ctx->PMGN,1))
        goto GDAFin;

    if (ctx->GD_GT != 2 && ctx->GD_GT != 4) {
        printf1(ctx, "Error: command requires valued graph.\n");
        goto GDAFin;
    }

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    if (ctx->PRCNAlloc == 0) {
        printf1(ctx, "No recode information. Nothing done.\n");
        err = 0;
        goto GDAFin;
    }
    printf1(ctx, "Recode: %s\n",ctx->PRCN);
               
    if (alloc_acn(ctx, ctx->GD_NP))                 
        goto GDAFin;
    if (alloc_acm(ctx, ctx->GD_NP))                 
        goto GDAFin;

    for (i = 0; i < ctx->GD_NP; ++i)  
        ctx->AcM[i] = ctx->AcN[i] = gdd_node(ctx, i);

    err = -2;
    p = ctx->PRCN;
    while (*p) {
        if (sscanf(p,"%d",&n) != 1)
            goto GDAFin;
        if (n < 1) {
            err = -4;
            goto GDAFin;
        }
        p = skip_int(ctx, p);
        if (*p++ != '[')
            goto GDAFin;

        while (*p) {
            if (sscanf(p,"%d",&m) != 1)
                goto GDAFin;

            fnd = 0;
            for (i = 0; i < ctx->GD_NP; ++i) {
                if (m == ctx->AcN[i]) {
                    ctx->AcM[i] = n;
                    fnd = 1;
                    break;
                }
            }
            if (fnd == 0) {
                err = -3;
                goto GDAFin; 
            }
            p = skip_int(ctx, p);
            if (*p == ']')
                break;

            if (*p++ != ',')   
                goto GDAFin;

            if (*p != ',')
                continue;

            if (sscanf(++p,"%d",&mm) != 1)
                goto GDAFin;

            fnd = 0;
            for (i = 0; i < ctx->GD_NP; ++i) {
                if (mm == ctx->AcN[i]) {
                    ctx->AcM[i] = n;
                    fnd = 1;
                    break;
                }
            }
            if (fnd == 0) {
                err = -3;
                goto GDAFin; 
            }
            for (j = m + 1; j < mm; ++j) {
                for (i = 0; i < ctx->GD_NP; ++i) {
                    if (j == ctx->AcN[i]) {
                        ctx->AcM[i] = n;
                        break;
                    }
                }
            }
            p = skip_int(ctx, p);
            if (*p == ']')
                break;

            if (*p++ != ',')   
                goto GDAFin;
        }
        if (*p++ != ']')
            goto GDAFin;

        if (*p && *p++ != ',')
            goto GDAFin;
    }
    for (i = 0; i < ctx->GD_NP; ++i)
        ctx->AcN[i] = ctx->AcM[i];

    if (sorti(ctx, ctx->GD_NP,ctx->AcN,0))
        goto GDAFin;

    mm = 1;
    for (i = 1; i < ctx->GD_NP; ++i) {
        if (ctx->AcN[i] != ctx->AcN[i - 1])
            ctx->AcN[mm++] = ctx->AcN[i];
    }
    nm = ctx->AcN[mm - 1];

    if (alloc_ack(ctx, nm + 1))                 
        goto GDAFin;

    for (i = 0; i < mm; ++i)  
        ctx->AcK[ctx->AcN[i]] = i;

    if (alloc_acx(ctx, mm * mm + 1))                 
        goto GDAFin;

    for (i = 0; i <= mm * mm; ++i)
        ctx->AcX[i] = -1.0;

    for (i = 0; i < ctx->GD_NP; ++i) {
        n_info(ctx, i,100);
        for (j = 0; j < ctx->GD_NP; ++j) {
            if (ctx->GD_GT == 2 && j > i)
                break;

            tmp = gdd_adj(ctx, i,j,ctx->PMGN);
            if (tmp >= 0.0) {
                ii = ctx->AcK[ctx->AcM[i]];
                jj = ctx->AcK[ctx->AcM[j]];
                       
                /** tda_out("i=%4d j=%4d ii=%5d jj=%5d : %d %d : v=%g      \n",i,j,ii,jj,AcN[ii],AcN[jj],tmp    ); **/
                       
                k = ii * mm + jj;           
                if (ctx->AcX[k] < 0.0)
                    ctx->AcX[k] = tmp;
                else
                    ctx->AcX[k] += tmp;
            }
        }
    }
    n_info_e(ctx);

    nrec = 0;
    for (i = 0; i < mm; ++i) {
        for (j = 0; j < mm; ++j) {
            k = i * mm + j;
            if (ctx->AcX[k] >= 0.0) {
                p_gdp(ctx, i + 1,j + 1,ctx->AcN[i],ctx->AcN[j]);
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcX[k]);
                fprintf(ctx->PMFd,"\n");                
                nrec++;    
            }
        }
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

GDAFin:       
    if (err == -2) {
        printf1(ctx, "Syntax error in recode string.\n");
        err = -1;
    }
    else if (err == -3) {
        printf1(ctx, "Error: unknown node number in recode string.\n");
        err = -1;
    }
    else if (err == -4) {
        printf1(ctx, "Error: node numbers must be positive.\n");
        err = -1;
    }
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdu         Union of two graphs.                                        */
/*                                                                          */
/*              gdu(                                                        */
/*                  df=...,         output file                             */
/*                  sc=...,         substitute for missing values, def. -1  */
/*                  nfmt=...,       integer print format, def. 4            */  
/*                  fmt=...,        print format for values, def. 10.4      */
/*              ) = I,J,V1,...;                                             */
/*                                                                          */  
/*  The command assumes that a graph (possibly a multigraph), called        */
/*  graph 1, is already defined with the perm option. In addition, the      */
/*  command expects on its right-hand side at least threee variables from   */
/*  the current data matrix that can be interpreted as an edge list         */
/*  defining a second graph, called graph 2. Interpretation is in the       */
/*  usual way: The first two variables provide node numbers that must be    */
/*  positive integers. Each additional variable specifies the edge values   */
/*  of a graph. In this way, also graph 2 can be a multigraph.              */
/*                                                                          */
/*  The command creates a new graph consisting of all nodes numbers and     */
/*  edges that are part of graph 1 or graph 2. The new graph is written as  */
/*  an edge list to the output file specified with the df parameter. If     */
/*  one of the graphs contains isolated nodes (without valid edge values)   */
/*  they are written at the end of the output file                          */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int gdu(TDAContext *ctx)
{
    register int i,j,k,l;        
    int err,n,nn,iv,jv,i1 = 0,i2,j1 = 0,j2,kp = 0,f,nrec,ne,ni;
    double tmp;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Union of two graphs. Current memory: %d bytes.\n",ctx->MemReq);

    if (gdd_check(ctx, 0,1))
        goto GDUFin;

    if (parm(ctx, ctx->CmdBuf + 3,4,1))       /* get parameters */
        goto GDUFin;
   
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    if (ctx->PMSCFlg == 0)
        ctx->PMSC = -1.0;

    if (ctx->PMNV < 3) {
        printf1(ctx, "Error: need at least three variables on right-hand side.\n");
        goto GDUFin;
    }
    nn = ctx->GD_NP + 2 * ctx->NOC;           
    if (alloc_acn(ctx, nn + 1))
        goto GDUFin;
    if (alloc_aci(ctx, ctx->NOC + 1))
        goto GDUFin;
    if (alloc_acj(ctx, ctx->NOC + 1))
        goto GDUFin;
    if (alloc_ack(ctx, ctx->NOC + 1))
        goto GDUFin;

    iv = (int)ctx->PMVIdx[0];
    jv = (int)ctx->PMVIdx[1];

    /* sort I and J in data matrix */

    for (i = 0; i < ctx->NOC; ++i) {
        ctx->AcI[i] = (int)get_data(ctx, iv,i);
        ctx->AcJ[i] = (int)get_data(ctx, jv,i);
        if (ctx->AcI[i] < 1 || ctx->AcJ[i] < 1) {
            printf1(ctx, "Error: node number in record %d is zero or negative.\n",i + 1);
            goto GDUFin;
        }
    }
    if (sortdpi2(ctx, ctx->NOC,ctx->AcI,ctx->AcJ,ctx->AcK))
        goto GDUFin;

    /* check for unique edges */

    i = ctx->AcI[ctx->AcK[0]];
    j = ctx->AcJ[ctx->AcK[0]];
    for (k = 1; k < ctx->NOC; ++k) {
        kp = ctx->AcK[k];
        if (ctx->AcI[kp] == i && ctx->AcJ[kp] == j) { 
            printf1(ctx, "Error: multiple occurrence of edge %d - %d.\n",i,j);
            goto GDUFin;
        }
        i = ctx->AcI[kp]; 
        j = ctx->AcJ[kp];
    }

    /* get list of all node numbers */

    k = 0;
    for (i = 0; i < ctx->GD_NP; ++i)
        ctx->AcN[k++] = gdd_node(ctx, i);

    for (i = 0; i < ctx->NOC; ++i) {
        ctx->AcN[k++] = ctx->AcI[i];
        ctx->AcN[k++] = ctx->AcJ[i];
    }
    if (sorti(ctx, nn,ctx->AcN,0))
        goto GDUFin;

    n = 1;
    for (i = 1; i < nn; ++i) {
        if (ctx->AcN[i] != ctx->AcN[i - 1])
            ctx->AcN[n++] = ctx->AcN[i];
    }
    printf1(ctx, "Union of graphs: %d nodes.\n",n);         

    /* create array to flag isolated nodes */

    nn = ctx->AcN[n - 1];     /* max node number */

    if (alloc_acns(ctx, nn + 1))
        goto GDUFin;

    for (i = 0; i < n; ++i)  
        ctx->AcNS[ctx->AcN[i]] = 1;

    /* write edges */

    ne = ni = nrec = 0;
    i = j = 0;

    for (k = 0; k <= ctx->NOC; ++k) {
        if (k < ctx->NOC) {
            kp = ctx->AcK[k];
            i2 = ctx->AcI[kp];
            j2 = ctx->AcJ[kp];
        }
        else  
            i2 = j2 = ctx->INTMAX;

        while (i < ctx->GD_NP && j < ctx->GD_NP) {
            i1 = gdd_node(ctx, i);
            j1 = gdd_node(ctx, j);
       
            f = 0;
            if (i1 < i2 || (i1 == i2 && j1 < j2)) {
                for (l = 1; l <= ctx->GD_NG; ++l) {
                    if (gdd_adj(ctx, i,j,l) >= 0.0) {
                        f = 1;
                        break;
                    }
                }
            }
            else if (i1 == i2 && j1 == j2)  
                f = 2;
            else
                break;

            if (f) {
                if (ctx->PMF1Def) {
                    rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,i1);
                    rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,j1);
                   
                    for (l = 1; l <= ctx->GD_NG; ++l) {
                        tmp = gdd_adj(ctx, i,j,l);
                        if (tmp < 0.0)
                            tmp = ctx->PMSC;
                        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,tmp);
                    }
                    for (l = 2; l < ctx->PMNV; ++l) {
                        tmp = ctx->PMSC;
                        if (f == 2 && k < ctx->NOC) {
                            tmp = get_data(ctx, (int)ctx->PMVIdx[l],kp);
                            if (tmp < 0.0)
                                tmp = ctx->PMSC;
                        }
                        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,tmp);
                    }
                    fprintf(ctx->PMF1d,"\n");                
                    nrec++;
                }
                ne++;
                ctx->AcNS[i1] = ctx->AcNS[j1] = 0;
            }
            if (++j >= ctx->GD_NP) {
                i++;
                j = 0;
            }
            if (f == 2)
                break;
        }
        if (k == ctx->NOC)
            break;

        if (i1 != i2 || j1 != j2) {

            f = 0;
            for (l = 2; l < ctx->PMNV; ++l) {
                if (get_data(ctx, (int)ctx->PMVIdx[l],kp) >= 0.0) {
                    f = 1;
                    break;
                }
            }
            if (f) {
                if (ctx->PMF1Def) {
                    rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,i2);
                    rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,j2);

                    for (l = 1; l <= ctx->GD_NG; ++l)  
                        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->PMSC);
                     
                    for (l = 2; l < ctx->PMNV; ++l) {
                        tmp = get_data(ctx, (int)ctx->PMVIdx[l],kp);
                        if (tmp < 0.0)
                            tmp = ctx->PMSC;
                        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,tmp);
                    }
                    fprintf(ctx->PMF1d,"\n");                
                    nrec++;
                }
                ne++;
                ctx->AcNS[i2] = ctx->AcNS[j2] = 0;
            }
        }
    }
    k = ctx->GD_NG + ctx->PMNV - 2;
    for (i = 1; i <= nn; ++i) {     /* check for isolated nodes */
        if (ctx->AcNS[i]) {
            if (ctx->PMF1Def) {
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,i);
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,i);
                for (l = 1; l <= k; ++l)  
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->PMSC);
                fprintf(ctx->PMF1d,"\n");                
                nrec++;
            }
            ni++;
        }
    }
    printf1(ctx, "Number of valid edges: %d\n",ne);
    printf1(ctx, "Number of isolated notes: %d\n",ni);
    if (ctx->PMF1Def)  
        printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);
    err = 0;

GDUFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gde         Create an undirected graph.                                 */
/*                                                                          */
/*              gde(                                                        */
/*                  df = ...,           output file                         */
/*                  prn = ...,          type of output, def. 0              */
/*                                      0 edge list                         */
/*                                      1 square matrix                     */
/*                  ni = ...,           1 without loops, def. 0             */
/*                  nfmt = ...,         integer print format, def. 4        */  
/*              ) = N,ID;                                                   */
/*                                                                          */
/*  The command requires two variables on the right-hand side. N must be    */
/*  positive integer and is interpreted as a list of node numbers. ID can   */
/*  be integer or floating point and is interpreted as the ID of an object  */
/*  that belongs to node N. The command constructs a undirected graph       */
/*  consisting of all different node numbers found in variable N. The       */
/*  value of an edge (i,j) is calculated as the number of times that an     */
/*  identical object (identified by ID) belongs both to node i and node j.  */
/*                                                                          */  
/*  Return: 0 if OK, -1 if error.                                           */

int gde(TDAContext *ctx)
{
    register int i,j,k = 0,l;          
    int err,n,ne,ii,iv,ni,nj,jp,lp,nmax;     
    double x;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Create an undirected graph. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 3,4,1))       /* get parameters */
        goto GDEFin;
   
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    if (ctx->PMNV != 2) {
        printf1(ctx, "Error: need exactly two variables on right-hand side.\n");
        goto GDEFin;
    }
    if (ctx->PMF1Def == 0) {
        printf1(ctx, "Error: need output file (df parameter).\n");
        goto GDEFin;
    }
    if (alloc_aci(ctx, ctx->NOC + 1))
        goto GDEFin;
    if (alloc_acj(ctx, ctx->NOC + 1))
        goto GDEFin;
    if (alloc_ack(ctx, ctx->NOC + 1))
        goto GDEFin;
    if (alloc_acx(ctx, ctx->NOC + 1))
        goto GDEFin;
    if (alloc_acy(ctx, ctx->NOC + 1))
        goto GDEFin;

    ii = (int)ctx->PMVIdx[0];
    iv = (int)ctx->PMVIdx[1];

    for (i = 0; i < ctx->NOC; ++i) {
        j = (int)get_data(ctx, ii,i);
        if (j < 1) {
            printf1(ctx, "Error: node number in record %d is zero or negative.\n",i + 1);
            goto GDEFin;
        }
        ctx->AcI[i] = j;                     
    }
    if (sortdpi(ctx, ctx->NOC,ctx->AcI,ctx->AcK))
        goto GDEFin;

    j = ctx->AcI[ctx->AcK[0]];
    ctx->AcJ[0] = j;
    n = 1;
    nj = 0;
    for (i = 1; i < ctx->NOC; ++i) {
        k = ctx->AcK[i];
        if (ctx->AcI[k] != j) {
            j = ctx->AcI[k];
            n++;
            ctx->AcJ[++nj] = j;
        }
    }
    nmax = ctx->AcJ[nj];
    nj++;
    printf1(ctx, "Number of nodes: %d\n",n);
   
    /**
    tda_out("AcJ nj=%d : ",nj);
    for (i = 0; i < nj; ++i)
        tda_out("%d ",AcJ[i]);
    newline(ctx);
    **/

    if (alloc_acr(ctx, nmax + 1))
        goto GDEFin;

    for (i = 0; i < nj; ++i)  
        ctx->AcR[ctx->AcJ[i]] = i;
    /**
    tda_out("AcR nm=%d : ",nmax);
    for (i = 0; i <= nmax; ++i)
        tda_out("%d ",AcR[i]);
    newline(ctx);
    **/

    if (alloc_acn(ctx, n * n + 1))
        goto GDEFin;

    for (i = 0; i < ctx->NOC; ++i) {
        ctx->AcX[i] = get_data(ctx, iv,i);
        ctx->AcY[i] = get_data(ctx, ii,i);
    }
    if (sortdp2(ctx, ctx->NOC,ctx->AcX,ctx->AcY,ctx->AcK))
        goto GDEFin;

    ne = 0;
    x = ctx->AcX[ctx->AcK[0]];
    ctx->AcI[0] = (int)ctx->AcY[ctx->AcK[0]];
    ni = 1;
    for (i = 1; i <= ctx->NOC; ++i) {
        if (i < ctx->NOC) 
            k = ctx->AcK[i];

        if (i == ctx->NOC || ctx->AcX[k] != x) {
            for (j = 1; j < ni; ++j) {
                jp = ctx->AcR[ctx->AcI[j]];
                for (l = 0; l < j; ++l) {
                    lp = ctx->AcR[ctx->AcI[l]];
                    ctx->AcN[jp * n + lp] += 1;
                    ctx->AcN[lp * n + jp] += 1;
                }
            }
            if (i == ctx->NOC)
                break;

            x = ctx->AcX[k];
            ctx->AcI[0] = (int)ctx->AcY[k];
            ni = 1;
        }
        else
            ctx->AcI[ni++] = (int)ctx->AcY[k];
    }
    /**
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j)
            tda_out("%4d ",AcN[i * n + j]);
        newline(ctx);
    }
    **/

    ne = ni = 0;
    if (ctx->PMPRNO == 0) {
        for (i = 0; i < n; ++i) {
            for (j = 0; j <= i; ++j) {
                if (j == i && ctx->PMNI == 1)
                    continue;

                k = ctx->AcN[i * n + j];
                if (k > 0) {
                    ne++;
                    if (ctx->PMF1Def) {
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->AcJ[j]);
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->AcJ[i]);
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,k);
                        fprintf(ctx->PMF1d,"\n");
                        ni++;
                    }
                }
            }
        }
    }
    else {
        for (i = 0; i < n; ++i) {
            if (ctx->PMF1Def)   
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->AcJ[i]);
            for (j = 0; j < n; ++j) {
                k = ctx->AcN[i * n + j];
                if (k > 0) {
                    if (ctx->PMNI == 0 || i != j)
                        ne++;
                }
                if (ctx->PMF1Def)   
                    rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,k);
            }
            fprintf(ctx->PMF1d,"\n");
            ni++;
        }
        if (ne > 1)
            ne /= 2;
    }
    printf1(ctx, "Number of edges: %d\n",ne);
    if (ctx->PMF1Def)  
        printf1(ctx, "%d records written to: %s\n",ni,ctx->PMF1dName);
    err = 0;

GDEFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdcset      Compact blocks for a directed and valued graph.             */
/*                                                                          */
/*              gdcset(                                                     */
/*                  gn=...,     graph number, def. 1                        */
/*                  nfmt=...,   integer print format, def. 4                */
/*                  fmt=...,    floating point print format, def. 10.4      */
/*              ) = fname;      output file required.                       */
/*                                                                          */
/*                                                                          */
/*  This command is in development and not finished yet.                    */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int gdcset(TDAContext *ctx)
{
    int err,m,r,first;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Find compact blocks. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 6,1,1))       /* get parameters */
        goto GDCSETFin;
   
    if (gdd_check(ctx, ctx->PMGN,1))
        goto GDCSETFin;

    if (ctx->GD_GT != 4) {
        printf1(ctx, "Error: command requires directed and valued graph.\n");
        goto GDCSETFin;
    }
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    if (alloc_acn(ctx, ctx->GD_NP + 1))
        goto GDCSETFin;
    if (alloc_acm(ctx, ctx->GD_NP + 1))
        goto GDCSETFin;
    if (alloc_aci(ctx, ctx->GD_NP + 1))
        goto GDCSETFin;
    if (alloc_acj(ctx, ctx->GD_NP + 1))
        goto GDCSETFin;

    for (m = 2; m < ctx->GD_NP; ++m) {

        first = 1;
        while (1) {
            if (comb_nm(ctx, ctx->GD_NP,m,ctx->AcN,first) == 0) 
                break;
            r = f_gdcset(ctx, m,ctx->AcN,ctx->AcM,ctx->AcI,ctx->AcJ);
            if (r < 0)
                goto GDCSETFin;
            first = 0;
        }
    }
    /**  tda_out("found nb=%d\n",nb); **/
    err = 0;

GDCSETFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  f_gdcset(m,com)                                                         */

int f_gdcset(TDAContext *ctx, int m,int *com,int *com1,int *nci,int *ncj)
{
    register int i,j,ii;
    int f,ne;
    double amin,aij;
/**
    tda_out("m=%3d : ",m);
    for (i = 0; i < m; ++i)
        tda_out("%2d ",com[i]);
    newline(ctx);
**/
    for (i = 0; i < m; ++i)
        com1[i] = com[i];

    ne = g_mst(ctx, ctx->PMGN,m,com1,nci,ncj,0);
    if (ne < m - 1)
        return(0);

    amin = ctx->DBLMAX;
    f = 0;
    for (i = 0; i < ne; ++i) {
        if ((aij = gdd_adj(ctx, nci[i],ncj[i],ctx->PMGN)) > 0.0) {
            amin = dmin(ctx, amin,aij);
            f = 1;
        }
        if ((aij = gdd_adj(ctx, ncj[i],nci[i],ctx->PMGN)) > 0.0) {
            amin = dmin(ctx, amin,aij);
            f = 1;
        }
    }
    if (f == 0) {
        tda_out("??????????\n");
        return(0);
    }
/*  tda_out("amin=%g\n",amin);
*/

    if (alloc_acr(ctx, ctx->GD_NP + 1))
        return(-1);       

    for (i = 0; i < m; ++i)
        ctx->AcR[com[i]] = 1;
     
    for (i = 0; i < m; ++i) {
        ii = com[i];
        for (j = 0; j < ctx->GD_NP; ++j) {
            if (ctx->AcR[j])
                continue;

            if ((aij = gdd_adj(ctx, ii,j,ctx->PMGN)) > 0.0) {
                if (aij >= amin)   
                    return(0);
            }
            if ((aij = gdd_adj(ctx, j,ii,ctx->PMGN)) > 0.0) {
                if (aij >= amin)  
                    return(0);
            }
        }
    }
    /**
    tda_out("found m=%3d : ",m);
    for (i = 0; i < m; ++i)
        tda_out("%2d ",com[i]);
    newline(ctx);
    **/

    return(1);
}

/* ------------------------------------------------------------------------ */
/*  gsym        Check whether a directed graph is symmetrical.              */
/*              gsym (gn=...);      gn = graph number, def. 1               */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int gsym(TDAContext *ctx)
{
    register int i,j;
    int err,n;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Check graph for symmetry. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,1,0))       /* get parameters */
        goto GSYMFin;
   
    if (gdd_check(ctx, ctx->PMGN,1))
        goto GSYMFin;

    if (ctx->GD_GT != 3 && ctx->GD_GT != 4) {
        printf1(ctx, "Error: command requires directed graph.\n");
        goto GSYMFin;
    }
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    n = 1;
    for (i = 1; i < ctx->GD_NP; ++i) {
        for (j = 0; j < i; ++j) {
            if (fabs(gdd_adj(ctx, i,j,ctx->PMGN) - gdd_adj(ctx, j,i,ctx->PMGN)) > 1.e-6) {
                n = 0;
                break;
            }   
        }
        if (n == 0)
            break;

    }
    printf1(ctx, "Graph is ");
    if (n == 0)
        printf1(ctx, "not ");
    printf1(ctx, "symmetrical.\n");
    if (n == 0) {
        printf1(ctx, "Nodes %5d %5d : %18.8f\n",i+1,j+1,gdd_adj(ctx, i,j,ctx->PMGN));
        printf1(ctx, "Nodes %5d %5d : %18.8f\n",j+1,i+1,gdd_adj(ctx, j,i,ctx->PMGN));
    }
    err = 0;

GSYMFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdot        Print directed graph into output file according to the      */
/*              dot syntax.                                                 */
/*                                                                          */
/*              gdot = fname;       output file name, required.             */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int gdot(TDAContext *ctx)
{
    register int i,j;
    int err,nrec;
    double tmp;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Print graph in dot format. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto GDOTFin;
   
    if (gdd_check(ctx, ctx->PMGN,1))
        goto GDOTFin;

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    nrec = 0;

    if (ctx->GD_GT == 1 || ctx->GD_GT == 2) {     /* undirected graph */
        fprintf(ctx->PMFd,"graph g {\n");
        for (i = 0; i < ctx->GD_NP; ++i) {
            fprintf(ctx->PMFd,"n%d [label=%d];\n",i + 1,gdd_node(ctx, i));
            nrec++;
        }
        for (i = 0; i < ctx->GD_NP; ++i) {
            /** ii = gdd_node(i); **/

            for (j = 0; j < i; ++j) {
                /** jj = gdd_node(j); **/

                tmp = gdd_adj(ctx, i,j,ctx->PMGN);
                if (tmp > 0.0) {
                    fprintf(ctx->PMFd,"n%d -- n%d;\n",i + 1,j + 1);
                    nrec++;
                }
            }
        }
    }
    else {          /* directed graph */
        fprintf(ctx->PMFd,"digraph g {\n");
        for (i = 0; i < ctx->GD_NP; ++i) {
            fprintf(ctx->PMFd,"n%d [label=%d];\n",i + 1,gdd_node(ctx, i));
            nrec++;
        }
        for (i = 0; i < ctx->GD_NP; ++i) {
            /** ii = gdd_node(i); **/

            for (j = 0; j < ctx->GD_NP; ++j) {
                /** jj = gdd_node(j); **/

                tmp = gdd_adj(ctx, i,j,ctx->PMGN);
                if (tmp > 0.0) {
                    fprintf(ctx->PMFd,"n%d -> n%d;\n",i + 1,j + 1);
                    nrec++;
                }
            }
        }
    }
    fprintf(ctx->PMFd,"}\n");
    nrec++;

    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

GDOTFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gtcl        Transitive closure                                          */
/*              If unvalued graph: Warshall's algorithm                     */
/*              If valued graph: Floyd's algorithm                          */
/*                                                                          */
/*              gtcl(                                                       */
/*                  gn = graph number, def. 1                               */
/*                  nfmt = ..., integer format, def. 4                      */
/*                  fmt = ...,  print format for values, def. 10.4          */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Output is always written as a square adj. matrix plus       */
/*              two leading columns containing node numbers.                */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gtcl(TDAContext *ctx)
{
    register int i,j,k;
    int err,typ;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Transitive closure. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto GTCLFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto GTCLFin;
    if (gdd_tcheck(ctx, 0,0,0))
        goto GTCLFin;

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    if (ctx->GD_GT == 1 || ctx->GD_GT == 3) {         /* unvalued graph */
        typ = 1;
        if (alloc_acc(ctx, ctx->GD_NP * ctx->GD_NP + 1))
            goto GTCLFin;

        g_tclos(ctx, ctx->PMGN,ctx->GD_NP,ctx->AcC);
    }
    else {                                  /* valued graph */
        typ = 0;
        if (alloc_acxf(ctx, ctx->GD_NP * ctx->GD_NP + 1))
            goto GTCLFin;

        g_spath(ctx, ctx->PMGN,ctx->GD_NP,ctx->AcXF);
    }
    for (i = 0; i < ctx->GD_NP; ++i) {

        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));

        k = i * ctx->GD_NP;
        for (j = 0; j < ctx->GD_NP; ++j) {
            if (typ)
                fprintf(ctx->PMFd,"%d ",(int)ctx->AcC[k++]);
            else
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,(double)ctx->AcXF[k++]);
        }
        fprintf(ctx->PMFd,"\n");
    }
    printf1(ctx, "%d records written to: %s\n",ctx->GD_NP,ctx->PMFdName);
    err = 0;

GTCLFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_tclos(gn,n,ai)                                                        */
/*                                                                          */
/*  Calculate transitive closure of graph gn with type gt. ai is an         */
/*  (n,n) character array, n = number of nodes.                             */

void g_tclos(TDAContext *ctx, int gn,int n,char *ai)
{
    register int i,j,k,ni,nk;

    for (i = 0; i < n; ++i) {
        ni = i * n;
        for (j = 0; j < n; ++j) {
            if (i == j || gdd_adj(ctx, i,j,gn) >= 0.0)
                ai[ni] = 1;
            else
                ai[ni] = 0;
            ni++;
        }
    }
    for (k = 0; k < n; ++k) {
        nk = k * n;
        for (i = 0; i < n; ++i) {
            ni = i * n;
            if (ai[ni + k]) {
                for (j = 0; j < n; ++j) {
                    if (ai[nk + j])
                        ai[ni + j] = 1;
                }
            }
        }
    }
    for (k = 0; k < n; ++k) {
        nk = k * n;
        if (gdd_adj(ctx, k,k,gn) < 0.0)
            ai[nk + k] = 0;
    }
}

/* ------------------------------------------------------------------------ */
/*  g_spath(gn,n,a)                                                         */
/*                                                                          */
/*  Calculate shortes patht matrix for a valued graph.                      */
/*                                                                          */
/*  Algorithm based on R.W. Floyd, Shortest Path.                           */
/*  Collected Algorithms from ACM 97.                                       */

void g_spath(TDAContext *ctx, int gn,int n,float *a)
{
    register int i,j,k,ni,nj;
    double s,max;

    max = (double)((float)(ctx->GD_EVMax[gn - 1])) * 10.0;

    for (i = 0; i < n; ++i) {
        ni = i * n;
        for (j = 0; j < n; ++j) {
            a[ni] = (float)(max);
            if (i != j) {
                s = (double)((float)gdd_adj(ctx, i,j,gn));          
                if (s >= 0.0)
                    a[ni] = (float)(s);
            }
            ni++;
        }
    }
    for (i = 0; i < n; ++i) {
        ni = i * n;
        for (j = 0; j < n; ++j) {
            nj = j * n;
            if ((double)(a[nj + i]) < max) {
                for (k = 0; k < n; ++k) {
                    if ((double)(a[ni + k]) < max) {
                        s = (double)(a[nj + i] + a[ni + k]);
                        if ((double)(s) < (double)(a[nj + k]))
                            a[nj + k] = (float)(s);
                    }   
                }
            }
        }
    }
    for (i = 0; i < n; ++i) {
        ni = i * n;
        for (j = 0; j < n; ++j) {
            if (i == j)
                a[ni] = 0.0;
            else if ((double)(a[ni]) >= max)  
                a[ni] = -1.0;
            ni++;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  gsort()     Topological sort. Needs GD_TYP = 1, directed.               */
/*                                                                          */
/*              gsort(                                                      */
/*                  opt = ...,  1 = list of new/old labels                  */
/*                              2 = graph with new labels                   */
/*                  gn = ...,   graph number, def 1                         */
/*                  nfmt = ..., integer print format, def. 4                */
/*                  fmt = ...,  print format, def. 10.4                     */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gsort(TDAContext *ctx)
{
    register int i,j,k,l;
    int err,n,nu,iia,iib,nr,nl;
    double tmp;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Topological sort. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,1,1))       /* get parameters */
        goto GSORTFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GSORTFin;
    if (gdd_tcheck(ctx, 1,2,0))
        goto GSORTFin;

    if (ctx->PMOPT > 2)
        ctx->PMOPT = 2;
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);
           
    if (alloc_acn(ctx, ctx->GD_NP))           /* in-degree */
        goto GSORTFin;
    if (alloc_acm(ctx, ctx->GD_NP))   
        goto GSORTFin;
    if (alloc_aci(ctx, ctx->GD_NP))       
        goto GSORTFin;
    if (alloc_acj(ctx, ctx->GD_NP))       
        goto GSORTFin;
    if (alloc_acr(ctx, ctx->GD_NP))       
        goto GSORTFin;
    if (alloc_acx(ctx, ctx->GD_NP))       
        goto GSORTFin;

    for (i = 0; i < ctx->GD_NP; ++i) {

        n = ctx->GD_FPN[i];
        if (n <= 0)
            continue;

        for (l = 0; l < n; ++l) {
            k = ctx->GD_FPK[i][l];
            if (k >= 0 && gdd_ev(ctx, k,ctx->PMGN) >= 0.0)   
                ctx->AcN[ctx->GD_FPI[i][l]] += 1;
        }
    }
    nl = iia = iib = nu = 0;
    for (i = 0; i < ctx->GD_NP; ++i) {
        if (ctx->AcN[i] == 0) {
            nu++;
            ctx->AcI[iib++] = i;
            ctx->AcR[i] = nl++;      /* new label */
        }
    }
    while (iia < iib) {
        i = ctx->AcI[iia++];
        n = ctx->GD_FPN[i];
        for (l = 0; l < n; ++l) {
            k = ctx->GD_FPK[i][l];
            if (k >= 0 && gdd_ev(ctx, k,ctx->PMGN) >= 0.0) { 
                j = ctx->GD_FPI[i][l];      /* edge from i to j */
                ctx->AcN[j] -= 1;
                if (ctx->AcN[j] == 0) {
                    ctx->AcI[iib++] = j;
                    ctx->AcR[j] = nl++;
                    nu++;
                }
            }
        }
    }   
    if (nu != ctx->GD_NP || iib != ctx->GD_NP) {
        printf1(ctx, "Sort not completed. There is at least one cycle or loop.\n");
    }    
    else {
        nr = 0;
        if (ctx->PMOPT == 1) {
            for (i = 0; i < ctx->GD_NP; ++i) {
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcR[i] + 1);
                fprintf(ctx->PMFd,"\n");
                nr++;
            }
        }
        else {
            for (i = 0; i < ctx->GD_NP; ++i) {
                nu = 0;
                iia = ctx->AcI[i];
                n = ctx->GD_FPN[iia];
                for (l = 0; l < n; ++l) {
                    k = ctx->GD_FPK[iia][l];
                    if (k >= 0 && (tmp = gdd_ev(ctx, k,ctx->PMGN)) >= 0.0) { 
                        j = ctx->GD_FPI[iia][l];      /* edge from iia to j */
                        ctx->AcN[nu] = j;
                        ctx->AcM[nu] = ctx->AcR[j];
                        ctx->AcX[nu++] = tmp;                    
                    }
                }
                if (nu > 0) {

                    if (sortdpi(ctx, nu,ctx->AcM,ctx->AcJ))
                        goto GSORTFin;
    
                    for (l = 0; l < nu; ++l) {
                        j = ctx->AcJ[l];        
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcR[iia] + 1);
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcR[ctx->AcN[j]] + 1);
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, iia));
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ctx->AcN[j]));
                        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcX[j]);
                        fprintf(ctx->PMFd,"\n");
                        nr++;
                    }
                }
            }
        }
        printf1(ctx, "%d records written to: %s\n",nr,ctx->PMFdName);
    }
    err = 0;

GSORTFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdcyc       All cycles in directed graphs. Requires a directed graph    */
/*                                             defined with gdd option 1.   */
/*              gdcyc(                                                      */
/*                  opt=...,    1 = one record for each cycle               */  
/*                              2 = one record for each node                */
/*                  max=...,    max length of cycles, def. any length.      */
/*                  ns=...,     max storage length, def. = # of nodes       */
/*                  gn=...,     graph number                                */
/*                  nfmt=...,   integer format, def. 4                      */
/*                  dopt=...,   break option, def. 0                        */
/*                              if dopt = 1 stop after the first cycle has  */
/*                              been found                                  */
/*              ) = fname;                                                  */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int gdcyc(TDAContext *ctx)
{
    int err,nc,ha,minclen,maxclen;
    int **h;

    ha = 0;
    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Cycles in digraphs. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 5,1,1))       /* get parameters */
        goto GDCYCFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GDCYCFin;

    if (gdd_tcheck(ctx, 1,2,0))      /* must be directed with gdd option 1 */
        goto GDCYCFin;
   
    if (ctx->PMOPT > 3)
        ctx->PMOPT = 3;
    if (ctx->PMDOPT != 1)
        ctx->PMDOPT = 0;

    if (alloc_acn(ctx, ctx->GD_NP))
        goto GDCYCFin;
    if (alloc_ack(ctx, ctx->GD_NP))
        goto GDCYCFin;

    if (!(h = (int **)calloc((size_t)(ctx->GD_NP),sizeof(int *)))) {
        p_err(ctx, -2,1);
        goto GDCYCFin; 
    }
    ha = ctx->GD_NP;
    memrq(ctx, ha,sizeof(int *));

    ctx->GSCON_CNT = 0;
    nc = g_dcyc(ctx, ctx->PMGN,ctx->GD_NP,ctx->AcN,h,ctx->AcK,ctx->PMOPT,ctx->PMNS,ctx->PMMax,&minclen,&maxclen,ctx->PMDOPT);
    if (nc < 0)
        goto GDCYCFin;

    printf1(ctx, "Number of cycles: %d\n",nc);
    if (nc > 0) {
        printf1(ctx, "Min length of cycles: %d\n",minclen);
        printf1(ctx, "Max length of cycles: %d\n",maxclen);
    }
    if (ctx->PMOPT <= 2)
        printf1(ctx, "%d records written to: %s\n",ctx->GSCON_CNT,ctx->PMFdName);
    err = 0;

GDCYCFin:       
    if (ha > 0) {
        free((char *)h);
        memrq(ctx, -ha,sizeof(int *));
    }
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_dcyc(gn,n,p,h,hp,opt,ns,cmax,minclen,maxclen,dopt)                    */
/*                                                                          */
/*  Enumeration of cycles in digraphs. Algorithm adapted from               */
/*  Tiernan,  An efficient search algorithm to find the elementary circuits */
/*  of a graph, Comm. ACM 13 (1970), 722 - 726.                             */
/*                                                                          */
/*  gn = graph number                                                       */
/*  n  = number of nodes                                                    */
/*  p  = array of size n                                                    */
/*  h  = array of size n with integer pointers                              */
/*  hp = array of size n                                                    */
/*  opt = 1 : write to PMFd, one record for each cycle                      */
/*  opt = 2 : write to PMFd, one record for each node in cycle              */
/*  In these cases report number of records written in GSCON_CNT.           */  
/*                                                                          */
/*  If dopt = 1 stop after the first cycle has been found.                  */
/*                                                                          */
/*  If ns > 0 this is max storage                                           */
/*  If cmax > 0 search only for cycles with max length cmax.                */
/*  Return min length of cycles found in minclen                            */
/*  Return max length of cycles found in maxclen                            */
/*                                                                          */
/*  Return number of cycles, or -1 if insufficient memory.                  */
/*  Return -2 if storage length insufficient.                               */

int g_dcyc(TDAContext *ctx, int gn,int n,int *p,int **h,int *hp,int opt,int ns,int cmax, int *minclen,int *maxclen,int dopt)
{
    register int k,i,j,l,pk;
    int ii,nc,ncc,nn,np,fnd,hmax;

    if (ns > 0)
        hmax = ns;
    else
        hmax = n;

    *minclen = n + 1;
    *maxclen = 0;
    for (i = 0; i < n; ++i)
        hp[i] = -1;                 /* if h[i] not allocated */

    nc = k = 0;
    p[0] = 0;
    np = 1;

    while (1) {
        pk = p[k];
        if (k == 0 && ctx->GD_BPN[pk] == 0)
            goto CCONT;

        fnd = 0;
        nn = ctx->GD_FPN[pk];
        for (i = 0; i < nn; ++i) {

            j = ctx->GD_FPI[pk][i];               
            if (ctx->GD_FPN[j] == 0)  
                continue;
                 
            ii = ctx->GD_FPK[pk][i];
            if (ii >= 0 && gdd_ev(ctx, ii,gn) >= 0.0) {  
         
                if (j <= p[0])
                    continue;

                fnd = 1;
                for (l = 0; l < np; ++l) {
                    if (p[l] == j) {
                        fnd = 0;
                        break;
                    }   
                }   
                if (fnd == 0)
                    continue;
 
                for (l = 0; l < hp[pk]; ++l) {
                    if (h[pk][l] == j) {
                        fnd = 0;
                        break;
                    }   
                }   
                if (fnd)
                    break;
            }
        }
        if (fnd) {
            p[++k] = j;
            np++;
            if (cmax > 0 && np > cmax)
                goto CCONT;
            continue;
        }

        fnd = 0;
        pk = p[k];

        nn = ctx->GD_FPN[pk];
        for (i = 0; i < nn; ++i) {
            j = ctx->GD_FPI[pk][i];               
            ii = ctx->GD_FPK[pk][i];
            if (ii >= 0 && gdd_ev(ctx, ii,gn) >= 0.0) {  
                if (p[0] == j) {
                    fnd = 1;  
                    break;
                }
            }
        }
        if (fnd) {
            nc++;
            if (*maxclen < np)
                *maxclen = np;
            if (*minclen > np)
                *minclen = np;

            if (opt == 1) {
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nc);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,np);
                for (j = 0; j < np; ++j)
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, p[j]));
                fprintf(ctx->PMFd,"\n");
                ctx->GSCON_CNT++;
            }
            else if (opt == 2) {
                for (j = 0; j < np; ++j) {
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nc);
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,np);
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, p[j]));
                    fprintf(ctx->PMFd,"\n");
                    ctx->GSCON_CNT++;
                }
            }
            if (dopt == 1)
                break;

            ncc = nc / 1000;
            if (ncc * 1000 == nc)
                tda_out("%10d\n",nc);


        }
CCONT:
        if (k == 0) {
            if (ctx->SILENTFlg < 2 && n >= 1000) {
                printfe(ctx, "Node: %7d     %c",p[0] + 1,CR);
                fflushe(ctx);
            }
            if (p[0] >= n - 1)  
                break;      

            p[0] += 1;
            k = 0;
            np = 1;

            for (j = 0; j < n; ++j) {
                if (hp[j] >= 0) {
                    free((char *)h[j]);
                    memrq(ctx, -hmax,sizeof(int));
                    hp[j] = -1;
                }
            }
        }
        else {
            pk = p[k];
            if (hp[pk] >= 0) {
                free((char *)h[pk]);
                memrq(ctx, -hmax,sizeof(int));
                hp[pk] = -1;
            }
            pk = p[k - 1];

            if (hp[pk] < 0) {
                if (!(h[pk] = (int *)calloc((size_t)(hmax),sizeof(int)))) {
                    p_err(ctx, -2,1);
                    nc = -1;       
                    break;
                }
                memrq(ctx, hmax,sizeof(int));
                hp[pk] = 0;
            }
            i = hp[pk];
            if (i >= hmax) {
                printf1(ctx, "Error: insufficient storage length.\n");
                nc = -2;
                break;
            }

            h[pk][i] = p[k];
            hp[pk] += 1;
            p[k] = 0;
            k--;
            np--;
        }
    }
    if (ctx->SILENTFlg < 2 && n >= 1000) {
        printfe(ctx, "\n");
        fflushe(ctx);
    }
    for (j = 0; j < n; ++j) {
        if (hp[j] >= 0) {
            free((char *)h[j]);
            memrq(ctx, -hmax,sizeof(int));
        }
    }
    return(nc);
}

/* ------------------------------------------------------------------------ */
/*  gev         Eigenvalues/vectors.  Requires an undirected graph.         */
/*                                                                          */
/*              gev(                                                        */
/*                  n=...,      number of eigenvalues/vectors, def. 1       */  
/*                  eps=...,    accuracy, def. 1.e-6                        */
/*                  mxit=...,   max number of iterations, def. 30           */
/*                  gn = graph number,                                      */
/*                  fmt = ...,  print format, def. 12.6                     */
/*                  prot=...,   protocol file with checks                   */
/*              ) = fname;                                                  */
/*                                                                          */
/*  If a protocol file is requested, then for each eigenvalue/vector,       */
/*  the maximal absolute difference between Cx and lambda * x is            */
/*  calculated and written into the protocol file.                          */
/*                                                                          */  
/*  Return: 0 if OK, -1 if error.                                           */

int gev(TDAContext *ctx)
{
    register int i,j;
    int err,r,n,nit;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Eigenvalues/vectors. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto GEVFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GEVFin;

    if (gdd_tcheck(ctx, 0,1,0))      /* must be undirected */
        goto GEVFin;
   
    if (ctx->PMOPT > 2)
        ctx->PMOPT = 2;
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 13,6);
    if (ctx->MxItFlg == 0)
        ctx->MxIter = 30;
    if (ctx->PMN < 1)
        ctx->PMN = 1;

    if (ctx->PMN >= ctx->GD_NP) { 
        printf1(ctx, "Error: n=%d exceeds maximum for current graph.\n",ctx->PMN);
        goto GEVFin;
    }
    n = ctx->PMN + 1;
    if (alloc_acx(ctx, ctx->GD_NP * n + 1))
        goto GEVFin;
    if (alloc_acy(ctx, ctx->GD_NP + 1))
        goto GEVFin;

    printf1(ctx, "Number of eigenvalues/vectors: %d\n",ctx->PMN);
    printf1(ctx, "Requested accuracy: %g\n",ctx->PMEPS);
    printf1(ctx, "Max number of iterations: %d\n",ctx->MxIter);

    r = simitz(ctx, ctx->GD_NP,n,ctx->MxIter,ctx->PMEPS,ctx->PMN,ctx->AcX,ctx->AcY,&nit);
    if (r < 0)
        goto GEVFin;

    printf1(ctx, "Iterations performed: %d\n",nit);
    if (r < 1)
        printf1(ctx, "No success.\n");
    else {
        printf1(ctx, "\nEigenvalue(s)\n");
        for (i = 1; i <= r; ++i) {
            printf1(ctx, "%3d  ",i);
            rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcY[i]); 
            newline(ctx);
        }
        newline(ctx);

        if (ctx->PMProtFDef)   
            g_evcheck(ctx, r,ctx->GD_NP,ctx->AcX,ctx->AcY);

        for (i = 1; i <= ctx->GD_NP; ++i) {
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i);
            for (j = 0; j < r; ++j) 
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcX[j * ctx->GD_NP + i]);
            fprintf(ctx->PMFd,"\n");
        }
        printf1(ctx, "%d records written to: %s\n",ctx->GD_NP,ctx->PMFdName);
    }
    err = 0;

GEVFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_op1(n,z,w) returns w = Cz where C is the n x n matrix to be used      */
/*               for the calculation of eigenvalues/vectors.                */
/*               called by simitz().                                        */

void g_op1(TDAContext *ctx, int n,double *z,double *w)
{
    register int i,j,l;
    int m;
    double a,tmp;

    if (ctx->GD_TYP == 1) {                  /* edge list */

        for (i = 0; i < n; ++i) {
            tmp = 0.0;
            m = ctx->GD_FPN[i];
            if (m > 0) {
                for (l = 0; l < m; ++l) {
                    j = ctx->GD_FPI[i][l];            /* edge from i to j */
                    a = gdd_adj(ctx, i,j,ctx->PMGN);
                    if (a > 0.0)
                        tmp += a * z[j + 1];
                }
            }
            w[i + 1] = tmp;
        }
    }
    else {
        for (i = 0; i < n; ++i) {
            tmp = 0.0;
            for (j = 0; j < n; ++j) {
                a = gdd_adj(ctx, i,j,ctx->PMGN);
                if (a > 0.0)
                    tmp += a * z[j + 1];
            }
            w[i + 1] = tmp;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  g_evcheck(r,n,x,d)  check calculated eigenvalues/vectors                */
/*                      and print to protocol file.                         */

void g_evcheck(TDAContext *ctx, int r,int n,double *x,double *d)
{
    register int i,j;
    double tmp,*xp;

    if (alloc_acu(ctx, n + 1))
        return;   

    for (i = 1; i <= r; ++i) {
        xp = x + (i - 1) * n;
        g_op1(ctx, n,xp,ctx->AcU);
        tmp = 0.0;
        for (j = 1; j <= n; ++j)  
            tmp = dmax(ctx, tmp,fabs(ctx->AcU[j] - d[i] * xp[j]));

        fprintf(ctx->PMProtFd,"Check %3d : ",i);
        rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,tmp);
        fprintf(ctx->PMProtFd,"\n");
    }
}

/* -##--------------------------------------------------------------------- */
/*  gsp         Shortest paths.                                             */
/*                                                                          */
/*              gsp(                                                        */
/*                  opt = ...,  1 = only reachable nodes                    */
/*                              2 = square matrix plus two leading cols     */
/*                              3 = without leading columns.                */
/*                  gn = graph number, def. 1                               */
/*                  nfmt = ..., integer format, def. 4                      */
/*                  fmt = ...,  print format for values, def. 10.4          */
/*              ) = fname;                                                  */
/*                                                                          */
/*  NOTE: this command requires that the graph is defined by an edge list.  */
/*                                                                          */  
/*  Return: 0 if OK, -1 if error.                                           */

int gsp(TDAContext *ctx)
{
    register int i,j;
    int err,n;   
    double tmp;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Shortest paths. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto GSPFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto GSPFin;
    if (gdd_tcheck(ctx, 1,0,0))
        goto GSPFin;

    if (ctx->PMOPT > 3)
        ctx->PMOPT = 3;
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    if (alloc_acx(ctx, ctx->GD_NP))
        goto GSPFin;
    if (alloc_ack(ctx, ctx->GD_NP + 1))
        goto GSPFin;

    n = 0;
    for (i = 0; i < ctx->GD_NP; ++i) {

        g_sp(ctx, i,ctx->GD_NP,ctx->PMGN,ctx->AcX,ctx->AcK);      

        if (ctx->PMOPT <= 2) {
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);                
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));                
        }

        for (j = 0; j < ctx->GD_NP; ++j) {
            tmp = ctx->AcX[j];
            if (ctx->PMOPT == 1) {
                if (j == i)
                    continue;
                if (tmp < ctx->DBLMAX) {
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, j));                
                    /** fprintf(PMFd," ["); **/            
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);                
                    /** fprintf(PMFd,"] "); **/            
                }
            }
            else {
                if (tmp >= ctx->DBLMAX)
                    tmp = -1.0;
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);                
            }
        }
        fprintf(ctx->PMFd,"\n");                
        n++;
    }
    printf1(ctx, "%d records written to: %s\n",n,ctx->PMFdName);
    err = 0;

GSPFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_glist(i,nodes)    Get list of nodes that are reachable from node i.   */
/*                      Return number of nodes in nodes[]. Skip loops.      */

int g_glist(TDAContext *ctx, int i,int *nodes)
{
    register int j,l,n,m;

    m = 0;
    n = ctx->GD_FPN[i];  
    for (l = 0; l < n; ++l) {
        j = ctx->GD_FPI[i][l];   
        if (j == i)
            continue;
        if (gdd_adj(ctx, i,j,ctx->PMGN) >= 0.0)
            nodes[m++] = j;
    }
    return(m);
}

/* ------------------------------------------------------------------------ */
/*  g_sp(in,n,gn,d,w,nflg,                                                  */
/*                                                                          */
/*  Find shortest path from node in to all other nodes. Algorithm adapted   */
/*  from U. Pape, Shortest Path Lengths, CALGO 562.                         */  
/*                                                                          */
/*  in = node number                                                        */
/*  n  = number of nodes                                                    */
/*  gn = graph number                                                       */
/*  d[i] (i = 0,...,GD_NP-1)    list of shortest path values                */
/*  ndec[i] array used for double-sided deque.                              */
/*                                                                          */

void g_sp(TDAContext *ctx, int in,int n,int gn,double *d,int *ndec)
{
    (void)gn;        /* unused: the signature is shared */
    register int i,j,l,np;   
    int m,max;
    double dj,tmp;   

    for (i = 0; i < n; ++i) {
        d[i] = ctx->DBLMAX;
        ndec[i] = 0;
    }
    max = n + 2;
    d[in] = 0.0;
    ndec[in] = max;  
    np = i = in;

    while (1) {
        m = ctx->GD_FPN[i];
        for (l = 0; l < m; ++l) {
            j = ctx->GD_FPI[i][l];                   

            tmp = gdd_adj(ctx, i,j,ctx->PMGN);
            if (tmp >= 0.0) {
                dj = d[i] + tmp;
                if (dj < d[j]) {
                    d[j] = dj;
                    /** w[j] = i; **/
                    if (ndec[j] == 0) {
                        ndec[np] = j + 1;
                        np = j;
                        ndec[j] = max;
                    }
                    else if (ndec[j] < 0) {
                        ndec[j] = ndec[i];
                        ndec[i] = j + 1;
                        if (np == i)
                            np = j;
                    }
                }        
            }
        }
        j = ndec[i];
        if (j >= max)
            break;
        ndec[i] = -j;        
        i = j - 1;
        if (i < 0)
            gerr_exit(ctx, 104);
    }
}

/* ------------------------------------------------------------------------ */
/*  tnet   ###  Temporal networks.                                          */
/*                                                                          */
/*              tnet(                                                       */
/*                  opt = ...,  1 = basic information.                      */
/*                              2 = creates temporal subgraphs              */
/*                              3 = distribution of indegrees               */
/*                              4 = distribution of outdegrees              */
/*                              5 = layers in historical network            */
/*                                  (requires directed graph)               */
/*                  dopt=...,   0 = temporal (default)                      */
/*                              1 = historical                              */
/*                  nfmt = ..., integer format, def. 4                      */
/*                  fmt = ...,  print format for values, def. 10.4          */
/*              ) [= fname];                                                */
/*                                                                          */
/*  Uses always graph number 1,2,3.                                         */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int tnet(TDAContext *ctx)
{
    register int i,j,k,l;
    int err,m,n,ne,iv,is,it,ii,jj,nrec,nrec1,im,n0,nl,ln;
    double tmp,s,t,tmin,tmax;

    err = -1;
    nrec1 = nrec = 0;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Temporal networks. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 4,1,0))       /* get parameters */
        goto TNETFin;

    if (ctx->GD_NG < 3) {
        printf1(ctx, "Error: need multigraph with at least three subgraphs.\n");
        goto TNETFin;
    }
    if (ctx->GD_TYP != 1) {
        printf1(ctx, "Error: need multigraph defined as an edge list.\n");
        goto TNETFin;
    }
    if (ctx->PMOPT > 5)
        ctx->PMOPT = 1;
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);


    iv = ctx->GD_EV[0];
    is = ctx->GD_EV[1];
    it = ctx->GD_EV[2];

    tmin = ctx->DBLMAX;
    tmax = 0.0;
    ne = 0;
    for (i = 0; i < ctx->NOC; ++i) {

        tmp = get_data(ctx, iv,i);

        if (tmp >= 0.0) {
            ii  = (int)get_data(ctx, ctx->GD_EI,i);
            jj  = (int)get_data(ctx, ctx->GD_EJ,i);
            s   = get_data(ctx, is,i);
            t   = get_data(ctx, it,i);

            if (s < 0.0 || t < s) {
                printf1(ctx, "Error: invalid dates in edge %d (%d,%d)\n",i+1,ii,jj);
                goto TNETFin;
            }
            tmin = dmin(ctx, tmin,s);
            tmax = dmax(ctx, tmax,t);
            ne++;
        }
    }
    printf1(ctx, "\nValid temporal network.\n");
    printf1(ctx, "Number of nodes: %d\n",ctx->GD_NP);
    printf1(ctx, "Number of edges: %d\n",ne);
    printf1(ctx, "Time axis: %g to %g\n",tmin,tmax);
    printf1(ctx, "Interpretation: ");
    if (ctx->PMDOPT != 1) {
        ctx->PMDOPT = 0;
        printf1(ctx, "temporal\n\n");
    }
    else
        printf1(ctx, "historical\n\n");


    if (ctx->PMOPT == 1) {
        err = 0;
        goto TNETFin;
    }
    if (ctx->PMNTP < 1) {
        printf1(ctx, "Error: need time points defined with tp option.\n");
        goto TNETFin;
    }
    if (alloc_acn(ctx, ctx->GD_NP + 1))
        goto TNETFin;
    if (alloc_acm(ctx, ctx->GD_NP + 1))
        goto TNETFin;
    if (alloc_ack(ctx, ctx->GD_NP + 1))
        goto TNETFin;

    if (ctx->PMOPT == 2) {           /* temporal subgraphs */

        printf1(ctx, "Option 2: creating temporal/historical subgraphs\n\n");

        prnchar(ctx, ' ',ctx->PMFmt1 - 4,0);
        printf1(ctx, "Time ");
        printf1(ctx, "   nodes    edges\n");

        for (j = 0; j < ctx->PMNTP; ++j) {
            tmp = ctx->PMTP[j];

            for (i = 0; i < ctx->GD_NP; ++i)
                ctx->AcN[i] = 0;

            m = n = 0; 
            for (i = 0; i < ctx->NOC; ++i) {
                s   = get_data(ctx, is,i);
                t   = get_data(ctx, it,i);
                if (s <= tmp && (ctx->PMDOPT == 1 || t >= tmp)) {
                    m++;
                    ii  = (int)get_data(ctx, ctx->GD_EI,i);
                    k = gdd_fndni(ctx, ii);
                    if (k >= 0 && k < ctx->GD_NP)
                        ctx->AcN[k] += 1;

                    ii  = (int)get_data(ctx, ctx->GD_EJ,i);
                    k = gdd_fndni(ctx, ii);
                    if (k >= 0 && k < ctx->GD_NP)
                        ctx->AcN[k] += 1;
                }      
            }
            for (k = 0; k < ctx->GD_NP;++k) {
                if (ctx->AcN[k] > 0)
                    n++;
            }
            rt_printf1_d(ctx, ctx->PMFmtS,tmp);
            printf1(ctx, "%8d %8d\n",n,m);

            if (ctx->PMFDef) {
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);                
                fprintf(ctx->PMFd,"%8d %8d\n",n,m);
                nrec++;
            }
        }
        newline(ctx);
        if (ctx->PMFDef)
            printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);

        if (ctx->PMF1Def) {

            for (i = 0; i < ctx->NOC; ++i) {
                ii  = (int)get_data(ctx, ctx->GD_EI,i);
                jj  = (int)get_data(ctx, ctx->GD_EJ,i);
                s   = get_data(ctx, is,i);
                t   = get_data(ctx, it,i);
                fprintf(ctx->PMF1d,"%8d %8d ",ii,jj);

                for (j = 0; j < ctx->PMNTP; ++j) {
                    tmp = ctx->PMTP[j];
                    if (s <= tmp && (ctx->PMDOPT == 1 || t >= tmp)) 
                        fprintf(ctx->PMF1d,"  1");
                    else
                        fprintf(ctx->PMF1d," -1");
                }
                fprintf(ctx->PMF1d,"\n");
                nrec1++;
            }
            printf1(ctx, "%d records written to: %s\n",nrec1,ctx->PMF1dName);
        }
        err = 0;
        goto TNETFin;
    }

    if (ctx->PMOPT == 3 || ctx->PMOPT == 4) {           /* in/outdegrees */

        printf1(ctx, "Option %d: distribution of ",ctx->PMOPT);
        if (ctx->PMOPT == 3)                   
            printf1(ctx, "in");
        else
            printf1(ctx, "out");

        printf1(ctx, "degrees\nInterpretation: ");
        if (ctx->GD_GT <= 2)
            printf1(ctx, "un");
        printf1(ctx, "directed\n\n");

        prnchar(ctx, ' ',ctx->PMFmt1 - 4,0);
        printf1(ctx, "Time ");
        printf1(ctx, "   nodes    edges   frequency of in- or outdegrees (0,1,2,...)\n");

        for (j = 0; j < ctx->PMNTP; ++j) {
            tmp = ctx->PMTP[j];

            im = 0;
            for (i = 0; i < ctx->GD_NP; ++i)
                ctx->AcK[i] = ctx->AcM[i] = ctx->AcN[i] = 0;

            m = n = 0; 
            for (i = 0; i < ctx->NOC; ++i) {
                ii  = (int)get_data(ctx, ctx->GD_EI,i);
                jj  = (int)get_data(ctx, ctx->GD_EJ,i);
                s   = get_data(ctx, is,i);
                t   = get_data(ctx, it,i);
                if (s <= tmp && (ctx->PMDOPT == 1 || t >= tmp)) {
                    m++;
                    k = gdd_fndni(ctx, ii);
                    if (k >= 0 && k < ctx->GD_NP) {
                        ctx->AcN[k] += 1;
                        if (ctx->GD_GT <= 2 || ctx->PMOPT == 4)
                            ctx->AcM[k] += 1;
                    }
                    k = gdd_fndni(ctx, jj);
                    if (k >= 0 && k < ctx->GD_NP) {
                        ctx->AcN[k] += 1;
                        if (ctx->GD_GT <= 2 || ctx->PMOPT == 3)
                            ctx->AcM[k] += 1;
                    }
                }      
            }
            for (k = 0; k < ctx->GD_NP;++k) {
                if (ctx->AcN[k] > 0) {
                    n++;
                    l = ctx->AcM[k];
                    im = imax(ctx, im,l);
                    ctx->AcK[l] += 1;
                }
            }
            if (im < 10)
                im = imin(ctx, 10,ctx->GD_NP);

            rt_printf1_d(ctx, ctx->PMFmtS,tmp);
            printf1(ctx, "%8d %8d ",n,m);
            for (k = 0; k <= im; ++k)
                rt_printf1_i(ctx, ctx->PMNFmtS,ctx->AcK[k]);
            newline(ctx);
       
            if (ctx->PMFDef) {
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);                
                fprintf(ctx->PMFd,"%8d %8d ",n,m);
                for (k = 0; k <= im; ++k)
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcK[k]);
                fprintf(ctx->PMFd,"\n");
                nrec++;
            }
        }
        newline(ctx);
        if (ctx->PMFDef)
            printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);

        err = 0;
        goto TNETFin;
    }

    if (ctx->PMOPT == 5) {           /* layers */

        printf1(ctx, "Option 5: layers in historical network.\n\n");
        if (ctx->GD_GT < 3) {                 
            printf1(ctx, "Error: requires directed graph.\n");
            goto TNETFin;
        }
        prnchar(ctx, ' ',ctx->PMFmt1 - 4,0);
        printf1(ctx, "Time ");
        printf1(ctx, " layer  members  node numbers\n");

        for (j = 0; j < ctx->PMNTP; ++j) {
            tmp = ctx->PMTP[j];

            for (i = 0; i < ctx->GD_NP; ++i) {
                ctx->AcK[i] = 0;
                ctx->AcN[i] = -1;
            }

            /* find members of layer 0 in AcN, number is n0 */

            ln = n0 = 0;
            for (i = 0; i < ctx->NOC; ++i) {
                ii  = (int)get_data(ctx, ctx->GD_EI,i);
                jj  = (int)get_data(ctx, ctx->GD_EJ,i);
                s   = get_data(ctx, is,i);
                if (s <= tmp) {
                    ii  = (int)get_data(ctx, ctx->GD_EI,i);
                    k = gdd_fndni(ctx, ii);
                    if (ctx->AcK[k] == 0)
                        ctx->AcK[k] = 1;  
                    jj  = (int)get_data(ctx, ctx->GD_EJ,i);
                    k = gdd_fndni(ctx, jj);
                    ctx->AcK[k] = 2;  
                }
            }
            for (k = 0; k < ctx->GD_NP; ++k) {
                if (ctx->AcK[k] == 1) {
                    ctx->AcN[k] = 0;
                    n0++;
                }
            }
            if (n0 == 0)
                continue;

            nl = 0;
            for (k = 0; k < ctx->GD_NP; ++k) {
                if (ctx->AcN[k] == 0)
                    nl++;
            }
            rt_printf1_d(ctx, ctx->PMFmtS,tmp);
            printf1(ctx, "%5d %8d ",ln,nl);

            if (ctx->PMFDef) {
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);                
                fprintf(ctx->PMFd,"%5d %8d ",ln,nl);
            }
            for (k = 0; k < ctx->GD_NP; ++k) {
                if (ctx->AcN[k] == 0) {
                    rt_printf1_i(ctx, ctx->PMNFmtS,ctx->GD_ND[k]);
                    if (ctx->PMFDef)  
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GD_ND[k]);
                }

            }
            newline(ctx);
            if (ctx->PMFDef) {
                fprintf(ctx->PMFd,"\n");
                nrec++;
            }
  
            while (1) {
                n0 = 0;
                for (i = 0; i < ctx->NOC; ++i) {
                    ii  = (int)get_data(ctx, ctx->GD_EI,i);
                    jj  = (int)get_data(ctx, ctx->GD_EJ,i);
                    s   = get_data(ctx, is,i);

                    if (s <= tmp) {
                        l = gdd_fndni(ctx, ii);
                        if (ctx->AcN[l] == ln) {
                            k = gdd_fndni(ctx, jj);
                            if (ctx->AcN[k] < 0) {
                                ctx->AcN[k] = ln + 1;
                                n0++;
                            }
                        }   
                    }
                }
                if (n0 == 0)
                    break;

                ln++;
                nl = 0;
                for (k = 0; k < ctx->GD_NP; ++k) {
                    if (ctx->AcN[k] == ln)
                        nl++;
                }
                rt_printf1_d(ctx, ctx->PMFmtS,tmp);
                printf1(ctx, "%5d %8d ",ln,nl);
                if (ctx->PMFDef) {
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);                
                    fprintf(ctx->PMFd,"%5d %8d ",ln,nl);
                }
                for (k = 0; k < ctx->GD_NP; ++k) {
                    if (ctx->AcN[k] == ln) {
                        rt_printf1_i(ctx, ctx->PMNFmtS,ctx->GD_ND[k]);
                        if (ctx->PMFDef)  
                            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GD_ND[k]);
                    }
                }
                newline(ctx);
                if (ctx->PMFDef) {
                    fprintf(ctx->PMFd,"\n");
                    nrec++;
                }
            }
        }
        newline(ctx);
        if (ctx->PMFDef)
            printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);

        err = 0;
        goto TNETFin;
    }


TNETFin:       
    p_clean(ctx);
    return(err);
}













