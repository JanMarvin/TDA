/****************************************************************************/
/*  t_gio                                                                   */
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
#include "t_ml.h"
#include "t_gf.h"
#include "t_gg.h"
#include "tda_context.h"

/*  functions in t_gio.c */

int gnc(TDAContext *ctx);
int gcon(TDAContext *ctx);
void u_visit(TDAContext *ctx, int k,int gn);
void u_visit1(TDAContext *ctx, int k,int gn);
int gdcon(TDAContext *ctx);
void d_visit(TDAContext *ctx, int k,int gn);
void d_visit1(TDAContext *ctx, int k,int gn); 
int gst(TDAContext *ctx);
void u1_visit(TDAContext *ctx, int k,int gn);
void u1_visit1(TDAContext *ctx, int k,int gn);
int gmst(TDAContext *ctx);
void u2_visit(TDAContext *ctx, int k,int gn);
void u2_visit1(TDAContext *ctx, int k,int gn);
int g_mst(TDAContext *ctx, int gn,int n,int *nodes,int *aci,int *acj,int mflag);
int g_mst1(TDAContext *ctx, int gn,int gtyp,int n,int *nodes,int *aci,int *acj,int mflag);
int g_ufind(TDAContext *ctx, int *dad,int x,int y,int doit);
int gcut(TDAContext *ctx);
int m_visit(TDAContext *ctx, int k,int gn);
int m_visit1(TDAContext *ctx, int k,int gn);
int gnst(TDAContext *ctx);
int gnst_grow(TDAContext *ctx, int n,char *g,short *s,char *b,char *a,int ns,int opt);
int gnst_prn(TDAContext *ctx, int n,char *a,int opt);
void put_bit(TDAContext *ctx, int i,char *bf,int v);
int get_bit(TDAContext *ctx, int i,char *bf);
void put_bit2(TDAContext *ctx, int i,int j,int n,char *bf,int v);
int get_bit2(TDAContext *ctx, int i,int j,int n,char *bf);
void cpy_bit2(TDAContext *ctx, int i,char *aci,int j,char *acj,int n);
int gcyc(TDAContext *ctx);
int g_cyc(TDAContext *ctx, int nc,int ne,int n,int *nodes,int opt);
int g_cycprn(TDAContext *ctx, int n,int *nodes,int i,int j,int *nptr,int *tptr,int *jptr, int ccnt,int nc,int opt,int **cptr,int *cptrn);
void g_cycprn1(TDAContext *ctx, int n,int *nodes,int nb,int **cptr,int *cptrn,int nc);
int g_cycprn2(TDAContext *ctx, int n,int *nodes,int ne,int nb,int **cptr,int *cptrn,int nc,int opt);
int gio(TDAContext *ctx);
void visit2(TDAContext *ctx, int k,double val,double val1);
int gfcf(TDAContext *ctx);
int gbcf(TDAContext *ctx);
void visit3(TDAContext *ctx, int k);
int gep(TDAContext *ctx);
void visit4(TDAContext *ctx, int k,int k1,int k2,double val,int first,int opt,int opt1);
int gflow(TDAContext *ctx);
double g_flow(TDAContext *ctx, int n,int i0,int i1,int gn,int *q,int *pre,char *vis,int *fptr,   double *flow,double *fd);
double g_flowbfs(TDAContext *ctx, int n,int i0,int i1,int gn,int *q,int *pre,char *vis,int *fptr,   double *flow,double *fd);
int gfc(TDAContext *ctx);
double g_cflow(TDAContext *ctx, int n,int i0,int i1,int gn,int *q,int *pre,char *vis,int *fptr,   double *flow,double *fd,int k0);
double g_cflowbfs(TDAContext *ctx, int n,int i0,int i1,int gn,int *q,int *pre,char *vis,int *fptr,   double *flow,double *fd,int k0);
int fc_visit(TDAContext *ctx, int k,int k0,int gn);




/* ------------------------------------------------------------------------ */
/*  gnc         node-centered networks.                                     */
/*                                                                          */
/*              gnc(                                                        */
/*                  gn = ...,       graph number,                           */
/*                  nfmt = ...,     integer format, def. 4                  */
/*                  fmt=...,        print format for values, def. 10.4      */
/*              ) = fname;                                                  */
/*                                                                          */
/*  The command assumes an undirected graph. The output file will contain   */
/*  one record for each node of the graph with following entries:           */
/*                                                                          */
/*  Column 1: internal node number                                          */
/*  Column 2: external node number                                          */
/*  Column 3: ni = degree of node i plus 1 (= number of nodes in the        */
/*            node-centered network of node i).                             */
/*  Column 4: ki = number of edges in the node-centered network             */
/*  Column 5: density, calculated as di = (2 ki)/(ni (ni - 1))              */
/*                                                                          */  
/*  Return: 0 if OK, -1 if error.                                           */

int gnc(TDAContext *ctx)
{
    register int i,j,k;
    int err,nrec,n,nn;
    double tmp;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Node-centered networks. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto GNCFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GNCFin;

    if (gdd_tcheck(ctx, 0,1,0))      /* must be undirected */
        goto GNCFin;         

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 10,4);

    if (alloc_acn(ctx, ctx->GD_NP))                 
        goto GNCFin;

    nrec = 0;
    for (i = 0; i < ctx->GD_NP; ++i) {
        n = 0;
        for (j = 0; j < ctx->GD_NP; ++j) {
            if (j == i)
                continue;
            tmp = gdd_adj(ctx, i,j,ctx->PMGN);
            if (tmp >= 0.0)  
                ctx->AcN[n++] = j;
        }   
        nn = n;
        for (j = 0; j < n; ++j) {
            for (k = j + 1; k < n; ++k) {
                if (gdd_adj(ctx, ctx->AcN[j],ctx->AcN[k],ctx->PMGN) >= 0.0)
                    nn++;
            }
        }
        n++;
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,n);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nn);
        if (n <= 1)
            tmp = 0.0;
        else
            tmp = 2.0 * (double)nn / ((double)n * (double)(n - 1));
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
        fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
        /* the density column is a ratio, and the file format rounds it
           to four decimals; without this the reader returned 0.6667
           for a value of 2/3 */
        {
            double erow[5];
            erow[0] = (double)(i + 1);
            erow[1] = (double)gdd_node(ctx, i);
            erow[2] = (double)n;
            erow[3] = (double)nn;
            erow[4] = tmp;
            tda_export_row(ctx, "gnc.table", erow, 5);
        }
#endif
        nrec++;
    }






    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

GNCFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gcon        Connected components. Graph is always undirected.           */
/*                                                                          */
/*              gcon(                                                       */
/*                  opt = ...,      1 = node list                           */
/*                                  2 = edge list                           */
/*                  gn = ...,       graph number,                           */
/*                  nfmt = ...,     integer format, def. 4                  */
/*                  fmt=...,        print format for values, def. 10.4      */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gcon(TDAContext *ctx)
{
    register int i,j,k;
    int err,nrec;
    double tmp;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Connected components. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto GCONFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GCONFin;

    if (gdd_tcheck(ctx, 0,1,0))      /* must be undirected */
        goto GCONFin;         

    if (ctx->PMOPT > 2)
        ctx->PMOPT = 2;
    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 10,4);

    if (alloc_acn(ctx, ctx->GD_NP))                 
        goto GCONFin;

    ctx->GSCON_ID = 1;
    i = nrec = 0;

    while (i >= 0) {

        if (ctx->SILENTFlg < 2 && ctx->GD_NP > 1000) {
            printfe(ctx, "Component: %7d     %c",ctx->GSCON_ID,CR);
            fflushe(ctx);
        }
        ctx->GSCON_CNT = 0;
        if (ctx->GD_TYP == 1)            /* edge list */
            u_visit(ctx, i,ctx->PMGN);
        else
            u_visit1(ctx, i,ctx->PMGN);

        i = -1;
        for (j = 0; j < ctx->GD_NP; ++j) {
            if (ctx->AcN[j] == ctx->GSCON_ID) {
                if (ctx->PMOPT == 1) {
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_ID);
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_CNT);
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j + 1);
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, j));
                    fprintf(ctx->PMFd,"\n");
                    nrec++;
                }
                else if (ctx->GSCON_CNT == 1) {
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_ID);
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_CNT);
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j + 1);
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j + 1);
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, j));
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, j));
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,gdd_adj(ctx, j,j,ctx->PMGN));
                    fprintf(ctx->PMFd,"\n");
                    nrec++;
                }
                else {
                    for (k = j + 1; k < ctx->GD_NP; ++k) {
                        if (ctx->AcN[k] == ctx->GSCON_ID) {
                            tmp = gdd_adj(ctx, j,k,ctx->PMGN);
                            if (tmp >= 0.0) {
                                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_ID);
                                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_CNT);
                                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j + 1);
                                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,k + 1);
                                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, j));
                                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, k));
                                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
                                fprintf(ctx->PMFd,"\n");
                                nrec++;
                            }
                        }
                    }
                }
            }
            else if (i < 0 && ctx->AcN[j] == 0)
                i = j;
        }
        ctx->GSCON_ID++;  
    }
    if (ctx->SILENTFlg < 2 && ctx->GD_NP > 1000) {
        printfe(ctx, "\n");
        fflushe(ctx);
    }
    printf1(ctx, "Number of components: %d\n",ctx->GSCON_ID - 1);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

GCONFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  u_visit(k,gn)   called by gcon. k is node number, gn is graph number.   */

void u_visit(TDAContext *ctx, int k,int gn)
{
    register int j,l,k1,k2;
    int n;

    ctx->AcN[k] = ctx->GSCON_ID;
    ctx->GSCON_CNT++;

    n = ctx->GD_FPN[k];
    if (n <= 0)
        return;          

    for (l = 0; l < n; ++l) {
        
        j = ctx->GD_FPI[k][l];                    /* edge from k to j */
        if (ctx->AcN[j] == 0) {
            k1 = ctx->GD_FPK[k][l];
            k2 = ctx->GD_FPK1[k][l];
            if ((k1 >= 0 && gdd_ev(ctx, k1,gn) >= 0.0) ||
                (k2 >= 0 && gdd_ev(ctx, k2,gn) >= 0.0))    
                u_visit(ctx, j,gn);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  u_visit1(k,gn)  called by gcon. k is node number, gn is graph number.   */

void u_visit1(TDAContext *ctx, int k,int gn)
{
    register int j;

    ctx->AcN[k] = ctx->GSCON_ID;
    ctx->GSCON_CNT++;

    for (j = 0; j < ctx->GD_NP; ++j) {

        if (ctx->AcN[j] == 0) {
            if (gdd_adj(ctx, k,j,gn) >= 0.0)
                u_visit1(ctx, j,gn);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  gdcon       Search for reachable nodes. Requires directed graph.        */
/*                                                                          */
/*              gdcon(                                                      */
/*                  opt=...,    1 = only number of reachable nodes          */
/*                              2 = additionally, the node numbers          */
/*                              3 = strong components, node list            */
/*                              4 = strong components, edge list            */
/*                  nfmt=...,   integer print format, def. 4                */
/*                  fmt=...,    print format for values, def. 10.4          */
/*                  gn=...,     graph number, def. 1                        */
/*                  if=...,     input file with node numbers,               */
/*              ) = fname;                                                  */
/*                                                                          */
/*  For options 3 and 4 the if parameter is ignored.                        */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int gdcon(TDAContext *ctx)
{
    register int i,j,k;
    int err,nrec,n,ni,maxl,nc;
    double v;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Search for reachable nodes. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 5,1,1))       /* get parameters */
        goto GDCONFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GDCONFin;

    if (gdd_tcheck(ctx, 0,2,0))      /* must be directed */
        goto GDCONFin;

    if (ctx->PMOPT > 4)
        ctx->PMOPT = 4;

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 10,4);

    ni = 0;
    if (ctx->PMOPT <= 2 && ctx->PMF2Def) {
        if (alloc_ack(ctx, ctx->GD_NP))                 
            goto GDCONFin;

        if ((ni = g_getnd(ctx, 0,ctx->AcK)) < 0)         /* get node list */
            goto GDCONFin;

        if (ni > 0)
            printf1(ctx, "Number of input nodes: %d\n",ni);

    }
    if (alloc_acn(ctx, ctx->GD_NP))                 
        goto GDCONFin;

    nrec = 0;

    if (ctx->PMOPT <= 2) {

        ctx->GSCON_NN = 0;
        ctx->GSCON_ID = 1;
        maxl = 0;
        for (i = 0; i < ctx->GD_NP; ++i) {

            if (ni > 0 && ctx->AcK[i] == 0)
                continue;

            for (j = 0; j < ctx->GD_NP; ++j)
                ctx->AcN[j] = 0;

            if (ctx->GD_TYP == 1)
                d_visit(ctx, i,ctx->PMGN);
            else
                d_visit1(ctx, i,ctx->PMGN);

            n = -1;
            for (j = 0; j < ctx->GD_NP; ++j) {
                if (ctx->AcN[j])
                    n++;
            }
            if (n > 0) {
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,n);
                if (ctx->PMOPT == 2) {
                    for (j = 0; j < ctx->GD_NP; ++j) {
                        if (ctx->AcN[j] && j != i)
                            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, j));
                    }
                }
                fprintf(ctx->PMFd,"\n");
                nrec++;
                ctx->GSCON_ID += 1;
                if (maxl < n)
                    maxl = n;
            }
            n_info(ctx, i + 1,100);
        }
        n_info_e(ctx);
        printf1(ctx, "Max number of reachable nodes: %d\n",maxl);
        printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    }
    else {                          /* directed graph */

        if (alloc_acm(ctx, ctx->GD_NP))                 
            goto GDCONFin;
        if (alloc_ack(ctx, ctx->GD_NP))                 
            goto GDCONFin;

        ctx->GSCON_NN = 0;
        ctx->GSCON_ID = 1;
        nc = 0;
                      
        for (i = 0; i < ctx->GD_NP; ++i) {

            if (ctx->AcK[i])
                continue;

            for (j = 0; j < ctx->GD_NP; ++j)            
                ctx->AcN[j] = 0;            
             
            if (ctx->GD_TYP == 1)
                d_visit(ctx, i,ctx->PMGN);
            else
                d_visit1(ctx, i,ctx->PMGN);
    
            /* save all node numbers that can be reached from i in AcM */
    
            for (j = 0; j < ctx->GD_NP; ++j)            
                ctx->AcM[j] = ctx->AcN[j];
    
            nc++;
            ctx->AcK[i] = nc;

            for (j = 0; j < ctx->GD_NP; ++j) {
                if (j != i && ctx->AcM[j]) {
                    for (k = 0; k < ctx->GD_NP; ++k)
                        ctx->AcN[k] = 0;

                    if (ctx->GD_TYP == 1)
                        d_visit(ctx, j,ctx->PMGN);
                    else
                        d_visit1(ctx, j,ctx->PMGN);

                    if (ctx->AcN[i])
                        ctx->AcK[j] = nc;
                }
            }
            n_info(ctx, i + 1,100);
        }
        n_info_e(ctx);
        printf1(ctx, "Number of strong components: %d\n",nc);

        if (ctx->PMOPT == 3) {
            for (i = 0; i < ctx->GD_NP; ++i) {
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcK[i]);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));
                fprintf(ctx->PMFd,"\n");
                nrec++;
            }
        }
        else {
            for (k = 1; k <= nc; ++k) {
                for (i = 0; i < ctx->GD_NP; ++i) {
                    if (ctx->AcK[i] != k)
                        continue;
                    n = 0;
                    for (j = 0; j < ctx->GD_NP; ++j) {
                        if (ctx->AcK[j] == k && (v = gdd_adj(ctx, i,j,ctx->PMGN)) >= 0.0) {
                            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,k);
                            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
                            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j + 1);
                            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));
                            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, j));
                            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,v);
                            fprintf(ctx->PMFd,"\n");
                            nrec++;
                            n++;
                        }
                    }
                    if (n == 0) {
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,k);
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));
                        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,-1.0);
                        fprintf(ctx->PMFd,"\n");
                        nrec++;
                    }
                }
            }
        }
        printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    }
    err = 0;

GDCONFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  d_visit(k,gn)   directed visit, called by gdcon. k is current node,     */
/*                  gn is graph number.                                     */

void d_visit(TDAContext *ctx, int k,int gn)
{
    register int j,l,kp;
    int n;

    ctx->AcN[k] = ctx->GSCON_ID;

    n = ctx->GD_FPN[k];
    if (n <= 0)
        return;          

    for (l = 0; l < n; ++l) {

        j = ctx->GD_FPI[k][l];               
        if (ctx->AcN[j] == 0) {
            kp = ctx->GD_FPK[k][l];
            if (kp >= 0 && gdd_ev(ctx, kp,gn) >= 0.0)    
                d_visit(ctx, j,gn);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  d_visit1(k,gn)  directed visit, called by gdcon. k is current node,     */
/*                  gn is graph number.                                     */

void d_visit1(TDAContext *ctx, int k,int gn)
{
    register int j;

    ctx->AcN[k] = ctx->GSCON_ID;

    for (j = 0; j < ctx->GD_NP; ++j) {

        if (ctx->AcN[j] == 0) {
            if (gdd_adj(ctx, k,j,gn) >= 0.0)
                d_visit1(ctx, j,gn);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  gst         Depth-first spanning tree. Graph is always undirected.      */
/*                                                                          */
/*              gst(                                                        */
/*                  gn = graph number,                                      */
/*                  nfmt = ..., integer format, def. 4                      */
/*                  fmt = ...,  print format for values, def. 10.4          */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gst(TDAContext *ctx)
{
    register int i,j;
    int err,nrec;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Spanning tree. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto GSTFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GSTFin;

    if (gdd_tcheck(ctx, 0,1,0))      /* must be undirected */
        goto GSTFin;

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 10,4);

    if (alloc_acn(ctx, ctx->GD_NP))                 
        goto GSTFin;
    if (alloc_aci(ctx, ctx->GD_NP))                 
        goto GSTFin;
    if (alloc_acj(ctx, ctx->GD_NP))                 
        goto GSTFin;

    ctx->GSCON_ID = 1;
    i = nrec = 0;

    while (i >= 0) {

        if (ctx->SILENTFlg < 2 && ctx->GD_NP > 500) {
            printfe(ctx, "Tree: %7d     %c",ctx->GSCON_ID,CR);
            fflushe(ctx);
        }
        ctx->GSCON_NN = 0;
        if (ctx->GD_TYP == 1)
            u1_visit(ctx, i,ctx->PMGN);
        else
            u1_visit1(ctx, i,ctx->PMGN);

        for (j = 0; j < ctx->GSCON_NN; ++j) {
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_ID);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ctx->AcI[j]));
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ctx->AcJ[j]));
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,gdd_adj(ctx, ctx->AcI[j],ctx->AcJ[j],ctx->PMGN)); 
            fprintf(ctx->PMFd,"\n");
            nrec++;
        }
        i = -1;
        for (j = 0; j < ctx->GD_NP; ++j) {
            if (ctx->AcN[j] == 0) {
                i = j;
                break;
            }
        }
        ctx->GSCON_ID++;  
    }
    if (ctx->SILENTFlg < 2 && ctx->GD_NP > 500) {
        printfe(ctx, "\n");
        fflushe(ctx);
    }
    printf1(ctx, "Number of trees: %d\n",ctx->GSCON_ID - 1);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

GSTFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  u1_visit(k,gn)   called by gcon. k is node number, gn is graph number.   */

void u1_visit(TDAContext *ctx, int k,int gn)
{
    register int j,l,k1,k2;
    int n;

    ctx->AcN[k] = ctx->GSCON_ID;

    n = ctx->GD_FPN[k];
    if (n <= 0)
        return;          

    for (l = 0; l < n; ++l) {

        j = ctx->GD_FPI[k][l];                    /* edge from k to j */
        if (ctx->AcN[j] == 0) {
            k1 = ctx->GD_FPK[k][l];
            k2 = ctx->GD_FPK1[k][l];
            if ((k1 >= 0 && gdd_ev(ctx, k1,gn) >= 0.0) ||
                (k2 >= 0 && gdd_ev(ctx, k2,gn) >= 0.0)) {  
                ctx->AcI[ctx->GSCON_NN] = k;
                ctx->AcJ[ctx->GSCON_NN++] = j;
                u1_visit(ctx, j,gn);
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  u1_visit1(k,gn)  called by gcon. k is node number, gn is graph number.   */

void u1_visit1(TDAContext *ctx, int k,int gn)
{
    register int j;

    ctx->AcN[k] = ctx->GSCON_ID;

    for (j = 0; j < ctx->GD_NP; ++j) {

        if (ctx->AcN[j] == 0) {
            if (gdd_adj(ctx, k,j,gn) >= 0.0) {
                ctx->AcI[ctx->GSCON_NN] = k;
                ctx->AcJ[ctx->GSCON_NN++] = j;
                u1_visit1(ctx, j,gn);
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  gmst        Minimum spanning tree. Graph is always undirected/valued.   */
/*                                     Loops are not used.                  */
/*              gmst(                                                       */
/*                  alg=...,        1 = Kruskal                             */
/*                                  2 = Prim                                */
/*                  max = ...,      1 for maximum spanning tree, def. 0     */
/*                  gn = graph number,                                      */
/*                  sort,       sort wrt to i,j                             */
/*                  nfmt = ..., integer format, def. 4                      */
/*                  fmt = ...,  print format for values, def. 10.4          */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gmst(TDAContext *ctx)
{
    register int i,j,k;
    int err,nrec,ne,nt,nte;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Minimum (maximum) spanning tree. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto GMSTFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GMSTFin;

    if (gdd_tcheck(ctx, 0,1,2))      /* must be undirected */
        goto GMSTFin;

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 10,4);

    if (ctx->PMMax != 1)
        ctx->PMMax = 0;

    if (ctx->PMALG != 2)
        ctx->PMALG = 1;

    printf1(ctx, "Algorithm: ");
    if (ctx->PMALG == 1)
        printf1(ctx, "Kruskal.\n");
    else  
        printf1(ctx, "Prim.\n");
         
    if (alloc_acn(ctx, ctx->GD_NP))                 
        goto GMSTFin;
    if (alloc_acm(ctx, ctx->GD_NP))                 
        goto GMSTFin;
    if (alloc_aci(ctx, ctx->GD_NP))                 
        goto GMSTFin;
    if (alloc_acj(ctx, ctx->GD_NP))                 
        goto GMSTFin;

    ctx->GSCON_ID = 1;
    nte = nt = i = nrec = 0;

    while (i >= 0) {

        if (ctx->SILENTFlg < 2 && ctx->GD_NP > 500) {
            printfe(ctx, "Tree: %7d     %c",ctx->GSCON_ID,CR);
            fflushe(ctx);
        }
        ctx->GSCON_NN = 0;
        if (ctx->GD_TYP == 1)
            u2_visit(ctx, i,ctx->PMGN);
        else
            u2_visit1(ctx, i,ctx->PMGN);
     
        if (ctx->GSCON_NN == 0)
            break;

        if (ctx->GSCON_NN == 1) {
            ctx->AcI[0] = ctx->AcJ[0] = ctx->AcM[0];
            ne = 1;
            nte++;
        }
        else {
            if (ctx->PMALG == 1)
                ne = g_mst1(ctx, ctx->PMGN,ctx->GD_TYP,ctx->GSCON_NN,ctx->AcM,ctx->AcI,ctx->AcJ,ctx->PMMax);
            else
                ne = g_mst(ctx, ctx->PMGN,ctx->GSCON_NN,ctx->AcM,ctx->AcI,ctx->AcJ,ctx->PMMax);         

            if (ne <= 0) {
                if (ne == 0)  
                    printf1(ctx, "ERROR: gmst.\n");
                goto GMSTFin;
            }   
        }
        if (ctx->PMSORTFlg) {
            if (alloc_ack(ctx, ne))
                goto GMSTFin;
            if (sortdpi2(ctx, ne,ctx->AcI,ctx->AcJ,ctx->AcK))
                goto GMSTFin;
        }
        for (j = 0; j < ne; ++j) {
            if (ctx->PMSORTFlg)
                k = ctx->AcK[j];
            else
                k = j;

            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_ID);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ctx->AcI[k]));
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ctx->AcJ[k]));
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,gdd_adj(ctx, ctx->AcI[k],ctx->AcJ[k],ctx->PMGN)); 
            fprintf(ctx->PMFd,"\n");
            nrec++;
        }
        nt++;
        ctx->GSCON_ID++;  

        i = -1;
        for (j = 0; j < ctx->GD_NP; ++j) {
            if (ctx->AcN[j] == 0) {
                i = j;
                break;
            }
        }
    }
    if (ctx->SILENTFlg < 2 && ctx->GD_NP > 500) {
        printfe(ctx, "\n");
        fflushe(ctx);
    }
    printf1(ctx, "Number of trees: %d",nt);
    if (nte > 0)
        printf1(ctx, " [includes %d isolated nodes]",nte);
    printf1(ctx, "\n%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

GMSTFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  u2_visit(k,gn)   called by gcon. k is node number, gn is graph number.   */

void u2_visit(TDAContext *ctx, int k,int gn)
{
    register int j,l,k1,k2;
    int n;

    ctx->AcN[k] = ctx->GSCON_ID;
    ctx->AcM[ctx->GSCON_NN++] = k;

    n = ctx->GD_FPN[k];
    if (n <= 0)
        return;          

    for (l = 0; l < n; ++l) {

        j = ctx->GD_FPI[k][l];                    /* edge from k to j */
        if (ctx->AcN[j] == 0) {
            k1 = ctx->GD_FPK[k][l];
            k2 = ctx->GD_FPK1[k][l];
            if ((k1 >= 0 && gdd_ev(ctx, k1,gn) >= 0.0) ||
                (k2 >= 0 && gdd_ev(ctx, k2,gn) >= 0.0)) {  
                u2_visit(ctx, j,gn);
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  u2_visit1(k,gn)  called by gcon. k is node number, gn is graph number.   */

void u2_visit1(TDAContext *ctx, int k,int gn)
{
    register int j;

    ctx->AcN[k] = ctx->GSCON_ID;
    ctx->AcM[ctx->GSCON_NN++] = k;

    for (j = 0; j < ctx->GD_NP; ++j) {

        if (ctx->AcN[j] == 0) {
            if (gdd_adj(ctx, k,j,gn) >= 0.0)  
                u2_visit1(ctx, j,gn);
        }
    }
}
    
/* ------------------------------------------------------------------------ */
/*  g_mst(gn,n,nodes,aci,acj,mflag)       MST with Prim's algorithm.        */
/*                                                                          */
/*  gn is graph number                                                      */
/*  nodes[i] (i = 0,...,n-1) contains node numbers for current component.   */
/*  at return, aci[j], acj[j] (j = 0,...,ne-1) contain the edges.           */
/*  By default: minimum spanning tree, if mflag != 0 then max spanning tree */
/*                                                                          */
/*  Return: ne, or -1 if error.                                             */

int g_mst(TDAContext *ctx, int gn,int n,int *nodes,int *aci,int *acj,int mflag)
{
    register int i,im,u,w;
    int ne;
    double tmp,mval,cmin;

    if (alloc_acx(ctx, n))
        return(-1);
    if (alloc_acr(ctx, n))
        return(-1);

    if (mflag == 0)
        mval = ctx->GD_EVMax[gn - 1] + 1.0;
    else
        mval = -1.0;

    u = nodes[0];
    nodes[0] = -1;
    im = -1;
    cmin = mval;
    for (i = 1; i < n; ++i) {
        tmp = gdd_adj(ctx, u,nodes[i],gn);
        if (tmp < 0.0)
            tmp = mval;
        ctx->AcX[i] = tmp;   
        ctx->AcR[i] = u;
        if (mflag == 0) {
            if (cmin > tmp) {
                cmin = tmp;
                im = i;
            }
        }
        else {
            if (cmin < tmp) {
                cmin = tmp;
                im = i;
            }
        }
    }
    ne = 0;
    while (im >= 0) {
      
        /* new edge from AcR[im] to nodes[im] */

        w = nodes[im];
        nodes[im] = -1;
        aci[ne] = imin(ctx, ctx->AcR[im],w);
        acj[ne++] = imax(ctx, ctx->AcR[im],w);          

        im = -1;
        cmin = mval;
        for (i = 0; i < n; ++i) {
            if (nodes[i] >= 0) {
                u = nodes[i];
                tmp = gdd_adj(ctx, w,u,gn);
                if (tmp >= 0.0) {
                    if (mflag == 0 && tmp < ctx->AcX[i]) {
                        ctx->AcR[i] = w;
                        ctx->AcX[i] = tmp;
                    }
                    else if (mflag == 1 && tmp > ctx->AcX[i]) {
                        ctx->AcR[i] = w;
                        ctx->AcX[i] = tmp;
                    }
                }
                if (mflag == 0) {
                    if (cmin > ctx->AcX[i]) {
                        cmin = ctx->AcX[i];
                        im = i;
                    }
                }
                else {
                    if (cmin < ctx->AcX[i]) {
                        cmin = ctx->AcX[i];
                        im = i;
                    }
                }
            }
        }
    }
    return(ne);
}

/* ------------------------------------------------------------------------ */
/*  g_mst1(gn,gtyp,n,nodes,aci,acj,mflag)   MST with Kruskal's algorithm.   */
/*                                                                          */
/*  gn is graph number.                                                     */
/*  gtyp is the gdd option.                                                 */
/*  nodes[i] (i = 0,...,n-1) contains node numbers for current component.   */
/*  at return, aci[j], acj[j] (j = 0,...,ne-1) contain the edges.           */
/*  By default: minimum spanning tree, if mflag != 0 then max spanning tree */
/*                                                                          */
/*  Return: ne, or -1 if error.                                             */

int g_mst1(TDAContext *ctx, int gn,int gtyp,int n,int *nodes,int *aci,int *acj,int mflag)
{
    register int i,j,l,ip;
    int err,ne,nemax,nn,da;
    int *dad;
    double tmp;

    err = -1;   

    nemax = ctx->GD_NE[gn - 1];
    da = 0;
    if (!(dad = (int *)calloc((size_t)(ctx->GD_NP + 2),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);
    }
    da = ctx->GD_NP + 2;
    memrq(ctx, da,sizeof(int));

    if (alloc_acr(ctx, nemax))
        goto MST1Fin;
    if (alloc_acs(ctx, nemax))
        goto MST1Fin;
    if (alloc_acx(ctx, nemax))
        goto MST1Fin;

    ne = 0;
    if (gtyp == 1) {                        /* edge list */
        for (i = 0; i < n; ++i) {
            ip = nodes[i];
            nn = ctx->GD_FPN[ip];
            for (l = 0; l < nn; ++l) {
                j = ctx->GD_FPI[ip][l];           /* edge from ip to j */
                if (j < ip) {
                    tmp = gdd_adj(ctx, ip,j,gn);
                    if (tmp >= 0.0) {
                        ctx->AcR[ne] = ip;
                        ctx->AcS[ne] = j;
                        ctx->AcX[ne++] = tmp;
                    }
                }
            }
        }
    }
    else {
        for (i = 1; i < n; ++i) {
            for (j = 0; j < i; ++j) {
                tmp = gdd_adj(ctx, nodes[i],nodes[j],gn);
                if (tmp >= 0.0) {
                    ctx->AcR[ne] = nodes[i];
                    ctx->AcS[ne] = nodes[j];
                    ctx->AcX[ne++] = tmp;
                }
            }
        }
    }
    if (alloc_ack(ctx, ne + 1))
        goto MST1Fin;

    if (sortdp(ctx, ne,ctx->AcX,ctx->AcK))
        goto MST1Fin;

    nn = 0;
    if (mflag == 0) {
        for (i = 0; i < ne; ++i) {
            ip = ctx->AcK[i];
            j = ctx->AcR[ip];
            l = ctx->AcS[ip];
                          
            if (g_ufind(ctx, dad,j + 1,l + 1,1)) {
                aci[nn] = imin(ctx, j,l);
                acj[nn++] = imax(ctx, j,l);
            }
        }
    }
    else {
        for (i = ne - 1; i >= 0; i--) {
            ip = ctx->AcK[i];
            j = ctx->AcR[ip];
            l = ctx->AcS[ip];
                          
            if (g_ufind(ctx, dad,j + 1,l + 1,1)) {
                aci[nn] = imin(ctx, j,l);
                acj[nn++] = imax(ctx, j,l);
            }
        }
    }
    err = nn;

MST1Fin:
    if (da > 0) {
        free((char *)dad);
        memrq(ctx, -da,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_ufind(dad,x,y,doit)   union find function, see Sedgewick 1990, p.446  */

int g_ufind(TDAContext *ctx, int *dad,int x,int y,int doit)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j,t;

    i = x;
    j = y;
    while (dad[i] > 0)  
        i = dad[i];
    while (dad[j] > 0)
        j = dad[j];

    while (dad[x] > 0) {
        t = x;
        x = dad[x];
        dad[t] = i;
    }
    while (dad[y] > 0) {
        t = y;
        y = dad[y]; 
        dad[t] = j;
    }
    if (doit && i != j) {
        if (dad[j] < dad[i]) {
            dad[j] += dad[i] - 1;
            dad[i] = j;
        }
        else {
            dad[i] += dad[j] - 1;
            dad[j] = i;
        }
    }
    return(i != j);
}

/* ------------------------------------------------------------------------ */
/*  gcut        Cut nodes and blocks. Graph is undirected.                  */
/*                                                                          */
/*              gcut(                                                       */
/*                  opt=...,    1 = one record for each component           */  
/*                              2 = one record for each node                */
/*                  gn = graph number,                                      */
/*                  nfmt = ..., integer format, def. 4                      */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gcut(TDAContext *ctx)
{
    register int i,j;
    int err,nrec,nb,nc,ncp,n,nid;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Cut nodes and blocks. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto GCUTFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GCUTFin;

    if (gdd_tcheck(ctx, 0,1,0))      /* must be undirected */
        goto GCUTFin;
   
    if (ctx->PMOPT > 2)
        ctx->PMOPT = 2;

    if (alloc_acn(ctx, ctx->GD_NP))                 
        goto GCUTFin;
    if (alloc_acm(ctx, ctx->GD_NP))                 
        goto GCUTFin;
    if (alloc_acns(ctx, ctx->GD_NP))                 
        goto GCUTFin;

    ctx->GSCON_ID = 0;
    nb = ncp = nc = i = nrec = 0;

    while (i >= 0) {

        if (ctx->SILENTFlg < 2 && ctx->GD_NP > 1000) {
            printfe(ctx, "Component: %7d     %c",nc+1,CR);
            fflushe(ctx);
        }

        ctx->GSCON_NSP = i;
        ctx->GSCON_CNT = 0;              /* number of paths from root */
        ctx->GSCON_NN  = 0;              /* number of elements in component */
        nid = ctx->GSCON_ID;
        if (ctx->GD_TYP == 1)
            m_visit(ctx, i,ctx->PMGN);
        else
            m_visit1(ctx, i,ctx->PMGN);
    
        nc++;
        if (ctx->GSCON_CNT < 2)
            ctx->AcNS[i] = 0;

        n = 0;
        for (j = 0; j < ctx->GD_NP; ++j) {
            if (ctx->AcNS[j]) {
                if (ctx->PMOPT == 1) {
                    ctx->AcM[n] = j;
                    ctx->AcNS[j] = 0;
                }
                n++;
                ncp++;
            }
        }
        if (n == 0)
            nb++;

        if (ctx->PMOPT == 1) {
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nc);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_NN);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,n);
            for (j = 0; j < n; ++j)
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ctx->AcM[j]));
            fprintf(ctx->PMFd,"\n");
            nrec++;
        }
        i = -1;
        for (j = 0; j < ctx->GD_NP; ++j) {
            if (ctx->AcN[j] > nid) {
                if (ctx->PMOPT == 2) {
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nc);
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j + 1);
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, j));
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,n);
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcNS[j]);
                    fprintf(ctx->PMFd,"\n");
                    nrec++;

                    if (ctx->AcNS[j] > 0)  
                        ctx->AcNS[j] = 0;
                }
            }
            else if (i < 0 && ctx->AcN[j] == 0)
                i = j;
        }
        ctx->GSCON_ID++;  
    }
    if (ctx->SILENTFlg < 2 && ctx->GD_NP > 1000) {
        printfe(ctx, "\n");
        fflushe(ctx);
    }
    printf1(ctx, "Number of components: %d\n",nc);
    printf1(ctx, "Number of cut nodes: %d\n",ncp);
    printf1(ctx, "Number of blocks: %d\n",nb);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

GCUTFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  m_visit(k,gn)   called by gcut. k is node number, gn is graph number.   */
/*                  GSCON_NSP = k in first call.                            */

int m_visit(TDAContext *ctx, int k,int gn)
{
    register int j,l,k1,k2;
    int n,m,min;

    ctx->AcN[k] = ++ctx->GSCON_ID;
    min = ctx->GSCON_ID;
    ctx->GSCON_NN++;
   
    n = ctx->GD_FPN[k];
    for (l = 0; l < n; ++l) {

        j = ctx->GD_FPI[k][l];                    /* edge from k to j */
        k1 = ctx->GD_FPK[k][l];
        k2 = ctx->GD_FPK1[k][l];
        if ((k1 >= 0 && gdd_ev(ctx, k1,gn) >= 0.0) ||
                      (k2 >= 0 && gdd_ev(ctx, k2,gn) >= 0.0)) {  

            if (ctx->AcN[j] == 0) {
               if (k == ctx->GSCON_NSP)
                   ctx->GSCON_CNT++;

                m = m_visit(ctx, j,gn);
                if (m < min)     
                    min = m;
     
                if (m >= ctx->AcN[k])  
                    ctx->AcNS[k] = 1;
            }
            else if (ctx->AcN[j] < min)   
                min = ctx->AcN[j];
        }
    }
    return(min);
}

/* ------------------------------------------------------------------------ */
/*  m_visit1(k,gn)  called by gcut. k is node number, gn is graph number.   */
/*                  GSCON_NSP = k in first call.                            */

int m_visit1(TDAContext *ctx, int k,int gn)
{
    register int j;
    int m,min;

    ctx->AcN[k] = ++ctx->GSCON_ID;
    min = ctx->GSCON_ID;
    ctx->GSCON_NN++;

    for (j = 0; j < ctx->GD_NP; ++j) {
        if (gdd_adj(ctx, k,j,gn) >= 0.0) {
            if (ctx->AcN[j] == 0) {
                if (k == ctx->GSCON_NSP)
                    ctx->GSCON_CNT++;

                m = m_visit1(ctx, j,gn);

                if (m < min)     
                    min = m;
     
                if (m >= ctx->AcN[k])    
                    ctx->AcNS[k] = 1;
            }
            else if (ctx->AcN[j] < min)   
                min = ctx->AcN[j];
        }
    }
    return(min);
}

/* ------------------------------------------------------------------------ */
/*  gnst        Enumeration of spanning trees. Requires an undirected graph */
/*              Uses CACM 354 separately for each component.                */
/*                                                                          */
/*              gnst(                                                       */
/*                  opt=...,    1 = one record for each tree                */  
/*                              2 = edge list                               */
/*                  gn = graph number,                                      */
/*                  nfmt = ..., integer format, def. 4                      */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gnst(TDAContext *ctx)
{
    register int i,j,ii;
    int err,n,ni,nc,in,jn;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Enumeration of spanning trees. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto GNSTFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GNSTFin;

    if (gdd_tcheck(ctx, 0,1,0))      /* must be undirected */
        goto GNSTFin;
   
    if (ctx->PMOPT > 2)
        ctx->PMOPT = 2;

    if (alloc_acn(ctx, ctx->GD_NP))                 
        goto GNSTFin;
    if (alloc_acm(ctx, ctx->GD_NP))                 
        goto GNSTFin;

    ctx->GSCON_NN = 0;
    ctx->GSCON_ID = 1;
    nc = ni = ii = 0;

    printf1(ctx, "\nComponent  nodes    trees\n");

    while (ii >= 0) {
        ctx->GSCON_CNT = 0;
        if (ctx->GD_TYP == 1)
            u_visit(ctx, ii,ctx->PMGN);
        else
            u_visit1(ctx, ii,ctx->PMGN);

        ii = -1;
        n = 0;
        for (j = 0; j < ctx->GD_NP; ++j) {
            if (ctx->AcN[j] == ctx->GSCON_ID)  
                ctx->AcM[n++] = j;
            else if (ii < 0 && ctx->AcN[j] == 0)
                ii = j;
        }
        if (n == 1)             /* isolated node */
            ni++;

        if (n > 0) {            /* new component */

            ctx->GSCON_CNT = 0;
            if (ctx->SILENTFlg < 2 && ctx->GD_NP > 10) {
                printfe(ctx, "Component: %7d     %c",ctx->GSCON_ID,CR);
                fflushe(ctx);
            }
            if (n > 1) {
                if (alloc_acc(ctx, n * n / 8 + 1))
                    goto GNSTFin;
                if (alloc_ace(ctx, n * n / 8 + 1))
                    goto GNSTFin;
                if (alloc_acd(ctx, n / 8 + 1))
                    goto GNSTFin;
                if (alloc_acns(ctx, n))
                    goto GNSTFin;

                for (i = 1; i < n; ++i) {
                    in = ctx->AcM[i];
                    for (j = 0; j < i; ++j) {
                        jn = ctx->AcM[j];
                        if (gdd_adj(ctx, in,jn,ctx->PMGN) >= 0.0) {
                            put_bit2(ctx, i,j,n,ctx->AcC,1);
                            put_bit2(ctx, j,i,n,ctx->AcC,1);
                        }
                    }
                }
                ctx->AcNS[0] = 1;
                for (j = 0; j < n; ++j) {
                    if (get_bit2(ctx, 0,j,n,ctx->AcC))
                        put_bit(ctx, j,ctx->AcD,1);
                }   
                if (gnst_grow(ctx, n,ctx->AcC,ctx->AcNS,ctx->AcD,ctx->AcE,1,ctx->PMOPT))
                    goto GNSTFin;      
            }
            printf1(ctx, "%6d %9d %8d\n",++nc,n,ctx->GSCON_CNT);
        }
        ctx->GSCON_ID++;  
    }
    if (ctx->SILENTFlg < 2 && ctx->GD_NP > 10) {
        printfe(ctx, "\n");
        fflushe(ctx);
    }
    printf1(ctx, "\n%d records written to: %s\n",ctx->GSCON_NN,ctx->PMFdName);
    err = 0;

GNSTFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gnst_grow   find all spanning trees for an undirected graph.            */
/*              the algorithm is adapted from:                              */
/*              M.D. McIlroy, Generator of Spanning Trees, CACM 354, 1969.  */
/*                                                                          */
/*  Return 0 if OK, -1 if error (insuff. memory)                            */

int gnst_grow(TDAContext *ctx, int n,char *g,short *s,char *b,char *a,int ns,int opt)
{
    register int i,j,k;
    int err,t,nfin,gfa,bfa;
    char *gf,*bf;

    err = gfa = bfa = 0;
     
    if (ns == n) {
        if (gnst_prn(ctx, n,a,opt))      /* print */
            return(-1);
    }
    else {
        if (!(gf = (char *)calloc((size_t)n * (size_t)n / 8 + 1,sizeof(char)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        gfa = n * n / 8 + 1;
        memrq(ctx, gfa,sizeof(char));
        memcpy(gf,g,(size_t)(gfa));

        for (i = 0; i < n; ++i) {
            if (get_bit(ctx, i,b)) {
                nfin = 0;
                for (j = 0; j < n; ++j) {
                    k = get_bit2(ctx, i,j,n,gf);
                    put_bit2(ctx, i,j,n,a,k * s[j]);

                    k *= (1 - s[j]);
                    put_bit2(ctx, i,j,n,gf,k);
                    nfin += k;
                }
                t = s[i]; 
                if (t == 0) {
                    s[i] = 1;
                    ns++;
                }
                if (!(bf = (char *)calloc((size_t)(n / 8 + 1),sizeof(char)))) {
                    p_err(ctx, -2,1);
                    err = -1;
                    break;           
                }
                bfa = n / 8 + 1;  
                memrq(ctx, bfa,sizeof(char));
                     
                for (j = 0; j < n; ++j) {
                    if (get_bit2(ctx, i,j,n,gf) || (j > i && get_bit(ctx, j,b)))
                        put_bit(ctx, j,bf,1);
                }
                err = gnst_grow(ctx, n,gf,s,bf,a,ns,opt);

                free((char *)bf);
                memrq(ctx, -bfa,sizeof(char));

                s[i] = (short)(t);
                if (t == 0)
                    ns--;
                if (err || nfin == 0)
                    break;         
            }
        }
    }
    if (gfa > 0) {
        free((char *)gf);
        memrq(ctx, -gfa,sizeof(char));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gnst_prn    Print trees created by gnst_grow. Print to PMFd.            */
/*              Count number of trees in GSCON_CNT, records in GSCON_NN.    */
/*              Use AcM[] for translation of node numbers.                  */

int gnst_prn(TDAContext *ctx, int n,char *a,int opt)
{
    register int i,j;
    int fin,m;

    if (alloc_acms(ctx, n))
        return(-1);    
    if (alloc_ack(ctx, n))
        return(-1);    

    m = 0;
    fin = 1;
    for (i = 0; i < n; ++i) {
        ctx->AcK[i] = -1;
        for (j = 0; j < n; ++j) {
            if (get_bit2(ctx, i,j,n,a)) {
                ctx->AcMS[i] = (short)(ctx->AcK[i] = j);
                m++;
                fin = 0;
                break;
            }
        }
    }
    while (fin == 0) {
        ctx->GSCON_CNT++;
        if (opt == 1) {
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_ID);                
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_CNT);                
            for (j = 0; j < n; ++j) {
                i = ctx->AcMS[j];
                if (ctx->AcK[j] >= 0)
                    i = gdd_node(ctx, ctx->AcM[i]);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i);
            }
            fprintf(ctx->PMFd,"\n");                
            ctx->GSCON_NN++;
        }
        else {
            i = 1;
            for (j = 0; j < n; ++j) {
                if (ctx->AcK[j] >= 0) {
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_ID);                
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_CNT);                
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,m);                
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i++);                
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ctx->AcM[j]));
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ctx->AcM[(int)ctx->AcMS[j]]));
                    fprintf(ctx->PMFd,"\n");                
                    ctx->GSCON_NN++;
                }
            }
        }
        fin = 1;
        for (i = n - 1; i >= 0; --i) {
            if (ctx->AcK[i] < 0)
                continue;

            for (j = ctx->AcMS[i] + 1; j < n; ++j) {
                if (get_bit2(ctx, i,j,n,a)) {
                    ctx->AcMS[i] = (short)(j);
                    fin = 0;
                    break;
                }
            }
            if (fin == 0)
                break;
            ctx->AcMS[i] = (short)(ctx->AcK[i]);
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  put_bit(i,bf,v)                                                         */

void put_bit(TDAContext *ctx, int i,char *bf,int v)
{
    register char *p;

    p = bf + i / 8;
    if (v)
        *p |= ctx->BMsk[i % 8];
    else
        *p &= ~ctx->BMsk[i % 8];
}

/* ------------------------------------------------------------------------ */
/*  get_bit(i,bf)                                                           */

int get_bit(TDAContext *ctx, int i,char *bf)
{
    register char *p;

    p = bf + i / 8;
    if (*p & ctx->BMsk[i % 8])
        return(1);
    else
        return(0);
}

/* ------------------------------------------------------------------------ */
/*  put_bit2(i,j,n,bf,v)                                                    */

void put_bit2(TDAContext *ctx, int i,int j,int n,char *bf,int v)
{
    register char *p;

    i = i * n + j;
    p = bf + i / 8;
    if (v)
        *p |= ctx->BMsk[i % 8];
    else
        *p &= ~ctx->BMsk[i % 8];
}

/* ------------------------------------------------------------------------ */
/*  get_bit2(i,j,n,bf)                                                      */

int get_bit2(TDAContext *ctx, int i,int j,int n,char *bf)
{
    register char *p;

    i = i * n + j;
    p = bf + i / 8;
    if (*p & ctx->BMsk[i % 8])
        return(1);
    else
        return(0);
}

/* ------------------------------------------------------------------------ */
/*  cpy_bit2(i,aci,j,acj,n)                                                 */
/*                                                                          */
/*  Copy row j in acj to row i in aci, n is number of columns               */

void cpy_bit2(TDAContext *ctx, int i,char *aci,int j,char *acj,int n)
{
    register int k,ii,jj;
    register char *pi,*pj;

    ii = i * n;
    jj = j * n;
    for (k = 0; k < n; ++k) {
        pi = aci + ii / 8;
        pj = acj + jj / 8;
        *pi &= ~ctx->BMsk[ii % 8];
        if (*pj & ctx->BMsk[jj % 8])
            *pi |= ctx->BMsk[ii % 8];
        ii++;
        jj++;
    }
}                     

/* ------------------------------------------------------------------------ */
/*  gcyc        Cycles in undirected graphs. Only with gdd option 1.        */
/*                                                                          */
/*              The algorithm to find a fundamental set of cycles is        */
/*              adapted from K. Paton, An Algorithm for Finding a           */
/*              Fundamental Set of Cycles of a Graph, Comm ACM 12, 1969,    */
/*              pp. 514 - 518. The algorithm to generate all cycles from    */
/*              the basis is taken from: N.E. Gibbs, A Cycle Generation     */
/*              Algorithm for Finite Undirected Linear Graphs, Journal      */
/*              of the ACM 16 (1969), pp. 564 - 568.                        */
/*                                                                          */
/*              gcyc                                                        */
/*                  opt=...,    1 = fundamental set of cycles, node list    */  
/*                              2 = fundamental set of cycles, edge list    */  
/*                              3 = all cycles, print version 1             */
/*                              4 = all cycles, print version 2             */
/*                                                                          */
/*                  gn = graph number,                                      */
/*                  nfmt = ..., integer format, def. 4                      */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gcyc(TDAContext *ctx)
{
    register int i,j,ii;
    int err,n,nc,nb,ne;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Cycles in undirected graphs. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto GCYCFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GCYCFin;

    if (gdd_tcheck(ctx, 1,1,0))      /* must be undirected and gdd option 1 */
        goto GCYCFin;
   
    if (ctx->PMOPT > 4)
        ctx->PMOPT = 4;

    if (alloc_acn(ctx, ctx->GD_NP))                 
        goto GCYCFin;
    if (alloc_acm(ctx, ctx->GD_NP))                 
        goto GCYCFin;

    printf1(ctx, "\nComponent  nodes  edges  fund cycles  all cycles\n");

    ctx->GSCON_NN = 0;
    ctx->GSCON_ID = 1;
    nc = ii = 0;

    while (ii >= 0) {

        if (ctx->SILENTFlg < 2 && ctx->GD_NP > 50) {
            printfe(ctx, "Component: %7d     %c",ctx->GSCON_ID,CR);
            fflushe(ctx);
        }

        ctx->GSCON_CNT = 0;
        if (ctx->GD_TYP == 1)
            u_visit(ctx, ii,ctx->PMGN);
        else
            u_visit1(ctx, ii,ctx->PMGN);

        ii = -1;
        n = 0;
        for (j = 0; j < ctx->GD_NP; ++j) {
            if (ctx->AcN[j] == ctx->GSCON_ID)  
                ctx->AcM[n++] = j;
            else if (ii < 0 && ctx->AcN[j] == 0)
                ii = j;
        }
        nc++;                   /* new component */
        ne = 0;                 /* number of edges */
        if (n > 0) {
            for (i = 0; i < n; ++i) {
                for (j = i + 1; j < n; ++j) {
                    if (gdd_adj(ctx, ctx->AcM[i],ctx->AcM[j],ctx->PMGN) >= 0.0)
                        ne++;
                }
            }
        }
        nb = ne - n + 1;        /* number of basic cycles, if not negative */
        ctx->GSCON_NP = 0;           /* number of all cycles */

        if (nb > 0) {
            nb = g_cyc(ctx, nc,ne,n,ctx->AcM,ctx->PMOPT);
            if (nb < 0)
                goto GCYCFin;
        }
        else
            nb = 0;

        printf1(ctx, "%6d %9d %6d  %8d   %10d\n",nc,n,ne,nb,ctx->GSCON_NP);

        ctx->GSCON_ID++;  
    }
    if (ctx->SILENTFlg < 2 && ctx->GD_NP > 50) {
        printfe(ctx, "\n");
        fflushe(ctx);
    }
    printf1(ctx, "\n%d records written to: %s\n",ctx->GSCON_NN,ctx->PMFdName);
    err = 0;

GCYCFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_cyc       Fundamental set of cycles. Requires an undirected graph.    */
/*                                                                          */
/*                                                                          */
/*              Return number of basic cycles, or -1 if error.              */

int g_cyc(TDAContext *ctx, int nc,int ne,int n,int *nodes,int opt)
{
    register int i,j,l,k1,k2;
    int err,it,ip,in,jn,nj,nb,cptra;
    int **cptr = NULL;

    err = -1;
    cptra = 0; 
    nb = ne - n + 1;        /* number of basic cycles, if not negative */
    if (nb <= 0)
        return(0);

    if (opt >= 2) {

        if (!(cptr = (int **)calloc((size_t)(nb),sizeof(int *)))) {
            p_err(ctx, -2,1);
            goto G_CYCFin;  
        }
        cptra = nb;    
        memrq(ctx, cptra,sizeof(int *));

        if (alloc_acs(ctx, nb))
            goto G_CYCFin; 
    }

    if (alloc_acr(ctx, ctx->GD_NP + 1))
        goto G_CYCFin; 
    if (alloc_acc(ctx, n))
        goto G_CYCFin; 
    if (alloc_acd(ctx, n))
        goto G_CYCFin; 
    if (alloc_aci(ctx, n))
        goto G_CYCFin; 
    if (alloc_acj(ctx, n))
        goto G_CYCFin; 
    if (alloc_ack(ctx, n))
        goto G_CYCFin; 

    for (i = 0; i < n; ++i)  
        ctx->AcR[nodes[i]] = i;

    nb = i = it = 0;
    ctx->AcD[0] = 1;
    ctx->AcI[0] = 0;
    ctx->AcJ[0] = -1;
    ip = 1;
     
    while (i >= 0) {
        in = nodes[i];
        nj = ctx->GD_FPN[in];
        for (l = 0; l < nj; ++l) {
            jn = ctx->GD_FPI[in][l];
            j = ctx->AcR[jn];

            if (ctx->AcD[j])  
                continue;

            k1 = ctx->GD_FPK[in][l];
            k2 = ctx->GD_FPK1[in][l];
            if ((k1 >= 0 && gdd_ev(ctx, k1,ctx->PMGN) >= 0.0) ||
                             (k2 >= 0 && gdd_ev(ctx, k2,ctx->PMGN) >= 0.0)) {  

                if (ctx->AcC[j]) {
                    if (g_cycprn(ctx, n,nodes,i,j,ctx->AcK,ctx->AcJ,ctx->AcI,nb,nc,opt,cptr,ctx->AcS))
                        goto G_CYCFin;
                    nb++;
                }
                else {
                    ctx->AcC[j] = 1;
                    ctx->AcK[j] = ip;
                    ctx->AcI[ip] = j;
                    ctx->AcJ[ip++] = ctx->AcK[i];
                    if (it < j)
                        it = j;
                }
            }
        }
        ctx->AcD[i] = 1;
        i = -1;
        for (j = it; j >= 0; --j) {
            if (ctx->AcC[j] != 0 && ctx->AcD[j] == 0) {
                i = j;
                it = j - 1;
                break;
            }
        }
    }
    if (opt == 2)  
        g_cycprn1(ctx, n,nodes,nb,cptr,ctx->AcS,nc);
    else if (opt >= 3)  {
        if (g_cycprn2(ctx, n,nodes,ne,nb,cptr,ctx->AcS,nc,opt))
            goto G_CYCFin;
    }
    err = nb;

G_CYCFin:
    if (cptra > 0) {
        for (i = 0; i < nb; ++i) {
            if (ctx->AcS[i] > 0) {
                free((char *)cptr[i]);
                memrq(ctx, -ctx->AcS[i],sizeof(int));
            }
        }
        free((char *)cptr);
        memrq(ctx, -cptra,sizeof(int *));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_cycprn    If opt == 1 print current cycle, if opt >= 2 save in        */
/*              cptr, cptrn.                                                */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int g_cycprn(TDAContext *ctx, int n,int *nodes,int i,int j,int *nptr,int *tptr,int *jptr, int ccnt,int nc,int opt,int **cptr,int *cptrn)
{
    register int k,ni,nj;
    int nn;

    if (alloc_acns(ctx, n))
        return(-1);    
    if (alloc_acms(ctx, n))
        return(-1);    

    ctx->AcNS[0] = (short)(i);
    ctx->AcMS[0] = (short)(j);
    ni = nj = 1;

    k = nptr[i];
    while ((k = tptr[k]) >= 0)  
        ctx->AcNS[ni++] = (short)(jptr[k]);

    k = nptr[j];
    while ((k = tptr[k]) >= 0)  
        ctx->AcMS[nj++] = (short)(jptr[k]);

    nn = 0;
    while (ni > 0 && nj > 0 && ctx->AcNS[--ni] == ctx->AcMS[--nj])  
        nn++; 
    if (nn == 0)  
        gerr_exit(ctx, 101);
    ni += 2;
    while (nj >= 0)
       ctx->AcNS[ni++] = ctx->AcMS[nj--];

    if (opt == 1) {
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nc);                
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ccnt + 1);                
        for (k = 0; k < ni; ++k)
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, nodes[ctx->AcNS[k]]));
        fprintf(ctx->PMFd,"\n");                
        ctx->GSCON_NN++;
    }
    else if (opt >= 2) {

        if (!(cptr[ccnt] = (int *)calloc((size_t)(ni),sizeof(int)))) {
            p_err(ctx, -2,1);
            return(-1);     
        }
        cptrn[ccnt] = ni;
        memrq(ctx, ni,sizeof(int));
        for (k = 0; k < ni; ++k)
            cptr[ccnt][k] = ctx->AcNS[k];
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  g_cycprn1   print fundamental cycles for option 2.                      */

void g_cycprn1(TDAContext *ctx, int n,int *nodes,int nb,int **cptr,int *cptrn,int nc)
{
    register int i,j,l,k,k1,k2;
    int nn,u,it,in,jn;

    nn = 0;
    for (i = 0; i < n; ++i) {
        in = nodes[i];
        for (j = i + 1; j < n; ++j) {
            jn = nodes[j];
            if (gdd_adj(ctx, in,jn,ctx->PMGN) >= 0.0) {

                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nc);                
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,++nn);                
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, nodes[i]));                
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, nodes[j]));                
                
                for (l = 0; l < nb; ++l) {
                    u = 0;
                    it = cptrn[l];
                    k1 = cptr[l][it - 1];
                    for (k = 0; k < it; ++k) {
                        k2 = cptr[l][k];
                        if (i == imin(ctx, k1,k2) && j == imax(ctx, k1,k2)) {
                            u = 1;
                            break;
                        }
                        k1 = k2;
                    }
                    fprintf(ctx->PMFd,"%d ",u);
                }
                fprintf(ctx->PMFd,"\n");
                ctx->GSCON_NN++;
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  g_cycprn2   print all cycles (option 3).                                */
/*              count number of cycles in GSCON_NP.                         */

int g_cycprn2(TDAContext *ctx, int n,int *nodes,int ne,int nb,int **cptr,int *cptrn,int nc,int opt) 
{
    register int i,j,l,k,k1,k2,iq;
    int u,it,in,jn,ie,qmax,nq,nq0,ns;

    /* first put edge list for fund cycles into AcC */

    if (alloc_acc(ctx, nb * ne / 8 + 1))
        return(-1);     

    ie = 0;
    for (i = 0; i < n; ++i) {
        in = nodes[i];
        for (j = i + 1; j < n; ++j) {
            jn = nodes[j];
            if (gdd_adj(ctx, in,jn,ctx->PMGN) >= 0.0) {

                for (l = 0; l < nb; ++l) {
                    u = 0;
                    it = cptrn[l];
                    k1 = cptr[l][it - 1];
                    for (k = 0; k < it; ++k) {
                        k2 = cptr[l][k];
                        if (i == imin(ctx, k1,k2) && j == imax(ctx, k1,k2)) {
                            u = 1;
                            break;
                        }
                        k1 = k2;
                    }
                    put_bit2(ctx, l,ie,ne,ctx->AcC,u);
                }
                ie++;
            }
        }
    }
    if (ie != ne)
        gerr_exit(ctx, 102);

    qmax = (int)pow(2.0,(double)nb) + 1;
    if (qmax < 1)
        gerr_exit(ctx, 103);

    if (alloc_acd(ctx, qmax * ne / 8 + 1))
        return(-1);     
    if (alloc_acns(ctx, qmax))
        return(-1);    

    cpy_bit2(ctx, 0,ctx->AcD,0,ctx->AcC,ne);
    ns = nq = 1;

    for (it = 1; it < nb; ++it) {

        if (alloc_ace(ctx, qmax / 8 + 1))
            return(-1);     

        nq0 = nq;    
        for (iq = 0; iq < nq0; ++iq) {
            u = 0;
            for (j = 0; j < ne; ++j) {
                k = 0;
                k1 = get_bit2(ctx, iq,j,ne,ctx->AcD);
                k2 = get_bit2(ctx, it,j,ne,ctx->AcC);
                if ((k1 != 0 && k2 == 0) || (k1 == 0 && k2 != 0))  
                    k = 1;
                put_bit2(ctx, nq,j,ne,ctx->AcD,k);
                if (k1 != 0 && k2 != 0)
                    u = 1;
            }
            if (u)
                put_bit(ctx, nq,ctx->AcE,1);
            else
                put_bit(ctx, nq,ctx->AcE,0);
            nq++;
        }
        for (iq = 0; iq < nq; ++iq) {
            if (get_bit(ctx, iq,ctx->AcE) == 0)
                continue;

            for (i = 0; i < iq; ++i) {
                if (get_bit(ctx, i,ctx->AcE) == 0)
                    continue;

                u = 1;
                for (j = 0; j < ne; ++j) {
                    if (get_bit2(ctx, iq,j,ne,ctx->AcD) != 0 && get_bit2(ctx, i,j,ne,ctx->AcD) == 0) {
                        u = 0;
                        break;
                    }   
                }
                if (u)  
                    put_bit(ctx, i,ctx->AcE,0);
                else {
                    u = 1;
                    for (j = 0; j < ne; ++j) {
                        if (get_bit2(ctx, i,j,ne,ctx->AcD) != 0 && get_bit2(ctx, iq,j,ne,ctx->AcD) == 0) {
                            u = 0;
                            break;
                        }   
                    }
                    if (u)  
                        put_bit(ctx, iq,ctx->AcE,0);
                }   
            }
        }
        cpy_bit2(ctx, nq,ctx->AcD,it,ctx->AcC,ne);
        for (iq = 0; iq < nq; ++iq) {
            if (get_bit(ctx, iq,ctx->AcE))
                ctx->AcNS[ns++] = (short)(iq);
        }
        ctx->AcNS[ns++] = (short)(nq);   
        nq++;
    }
    if (opt >= 3) {
        ie = 0;
        for (i = 0; i < n; ++i) {
            in = nodes[i];
            for (j = i + 1; j < n; ++j) {
                jn = nodes[j];
                if (gdd_adj(ctx, in,jn,ctx->PMGN) >= 0.0) {

                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nc);                
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ie + 1);                
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, nodes[i]));                
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, nodes[j]));                
                    if (opt == 3) {
                        for (k = 0; k < ns; ++k) {
                            l = ctx->AcNS[k];
                            fprintf(ctx->PMFd,"%d ",(int)get_bit2(ctx, l,ie,ne,ctx->AcD));
                        }
                    }
                    fprintf(ctx->PMFd,"\n");
                    ctx->GSCON_NN++;
                    ie++;
                }
            }
        }
    }
    if (opt == 4) {

        for (k = 0; k < ns; ++k) {
            l = ctx->AcNS[k];
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nc);                
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,k + 1);                

            for (i = 0; i < ne; ++i)
                fprintf(ctx->PMFd,"%d ",(int)get_bit2(ctx, l,i,ne,ctx->AcD));
            fprintf(ctx->PMFd,"\n");
            ctx->GSCON_NN++;
        }
    }
    ctx->GSCON_NP = ns;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  gio         Integrated Ownership. Graph is always directed, valued.     */
/*                                    gdd option 1 required                 */
/*              gio(                                                        */
/*                  mxit=...,   max number of iterations, def. 20           */
/*                  eps=...,    tolerance for concergence, def. 1.e-4       */
/*                  if=...,     input file with node numbers                */
/*                  gn = ...,   graph number, def. 1                        */
/*                  nfmt = ..., integer format, def. 4                      */
/*                  fmt=...,    print format for values, def. 7.4           */
/*                  opt=...,    output option                               */
/*                              1 = write all records                       */  
/*                              2 = only records with non-zero entries      */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gio(TDAContext *ctx)
{
    register int i,j,l,nn;
    int err,nrec,ii,iter,n,nni,ni,ne;
    double val,tmp;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Integrated ownership. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto GIOFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GIOFin;
    if (gdd_tcheck(ctx, 1,2,2))      /* gdd option 1 required */
        goto GIOFin;            /* and directed, valued */

    if (ctx->PMOPT > 2)
        ctx->PMOPT = 2;

    printf1(ctx, "Maximum of edge values: %lg\n",ctx->GD_EVMax[ctx->PMGN - 1]);
    if (ctx->GD_EVMax[ctx->PMGN - 1] > 1.0) {
        printf1(ctx, "Must not exceed 1.\n");
        goto GIOFin;
    }
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 7,4);

    if (ctx->MxItFlg == 0 || ctx->MxIter < 1)
        ctx->MxIter = 20;
    if (ctx->PMEPSFlg == 0)
        ctx->PMEPS = 1.e-4;

    printf1(ctx, "Max number of iterations: %d\n",ctx->MxIter);
    printf1(ctx, "Tolerance for convergence: %lg\n",ctx->PMEPS);

    ni = 0;
    if (ctx->PMF2Def) {
        if (alloc_ack(ctx, ctx->GD_NP))                 
            goto GIOFin;

        if ((ni = g_getnd(ctx, 0,ctx->AcK)) < 0)         /* get node list */
            goto GIOFin;

        if (ni > 0)
            printf1(ctx, "Number of input nodes: %d\n",ni);
    }
    if (alloc_acx(ctx, ctx->GD_NP))                 
        goto GIOFin;

    if (alloc_acy(ctx, ctx->GD_NP))                 
        goto GIOFin;

    if (alloc_acn(ctx, ctx->GD_NP))                 
        goto GIOFin;

    nni = nrec = 0;

    for (ii = 0; ii < ctx->GD_NP; ++ii) {

        if (ni > 0 && ctx->AcK[ii] == 0)
            continue;

        ctx->GSCON_ID = 0;
        for (i = 0; i < ctx->GD_NP; ++i) {
            ctx->AcN[i] = 0;
            ctx->AcX[i] = 0.0;
        }

        visit2(ctx, ii,1.0,1.0);

        n = 0;
        for (i = 0; i < ctx->GD_NP; ++i) {
            if (ctx->AcN[i])    
                n++;
        }
        if (n == 0) {
            nni++;
            continue;
        }

        for (iter = 0; iter < ctx->MxIter; ++iter) {

            for (i = 0; i < ctx->GD_NP; ++i) {
                if (ctx->AcX[i] > 0.0) {
    
                    nn = ctx->GD_BPN[i];
                    if (nn > 0) {
                        for (l = 0; l < nn; ++l) {
                            j = ctx->GD_BPI[i][l];          /* edge from j to i */
                            if (j != i) {
                                ne = ctx->GD_BPK[i][l];
                                if (ne >= 0 && (val = gdd_ev(ctx, ne,ctx->PMGN)) > 0.0) {  
                                    if (j == ii)
                                        ctx->AcY[i] += val;
                                    else
                                        ctx->AcY[i] += ctx->AcX[j] * val;
                                }
                            }
                        }
                    }   
                }
            }
            val = 0.0;
            for (i = 0; i < ctx->GD_NP; ++i) {
                if (ctx->AcN[i]) {
                    tmp = fabs(ctx->AcX[i] - ctx->AcY[i]);
                    if (val < tmp)
                        val = tmp;
                      
                    ctx->AcX[i] = ctx->AcY[i];
                    ctx->AcY[i] = 0.0;
                }
            }
            if (val <= ctx->PMEPS)
                break;
        }
        for (i = 0; i < ctx->GD_NP; ++i) {
            if (ctx->AcN[i]) {
                tmp = gdd_adj(ctx, ii,i,ctx->PMGN);
                if (tmp < 0.0)
                    tmp = 0.0;

                if (ctx->PMOPT == 2 && tmp < ctx->PMSC && ctx->AcX[i] < ctx->PMSC)
                    continue;

                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ii + 1);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ii));
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcX[i]);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,n);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,iter);
                fprintf(ctx->PMFd,"\n");
                nrec++;
            }
        }
        n_info(ctx, ii + 1,10);
    }
    n_info_e(ctx);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    if (nni > 0)
        printf1(ctx, "skipped %d isolated nodes.\n",nni);
    err = 0;

GIOFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  visit2      called by gio(). Graph always directed, valued.             */

void visit2(TDAContext *ctx, int k,double val,double val1)
{
    register int j,l,nn,kk;
    double tmp;

    ctx->AcN[k] = ++ctx->GSCON_ID;
    ctx->AcX[k] += val * val1;

    nn = ctx->GD_FPN[k];
    if (nn > 0) {
        for (l = 0; l < nn; ++l) {
            j = ctx->GD_FPI[k][l];                    /* edge from k to j */
            if (ctx->AcN[j] == 0) {
                kk = ctx->GD_FPK[k][l];
                if (kk >= 0) {
                    tmp = gdd_ev(ctx, kk,ctx->PMGN);
                    if (tmp >= 0.0)
                        visit2(ctx, j,tmp,val);
                }
            }
        }
    }
}
   
/* -##--------------------------------------------------------------------- */
/*  gfcf        Forward control flows.                                      */
/*              Graph type is always directed (unvalued).                   */
/*                                   gdd option 1 required.                 */
/*              gfcf(                                                       */
/*                  gn = graph number,                                      */
/*                  if = input file with list of nodes (optional)           */
/*                  sc = ...,           def. 0.5                            */
/*                  opt = 1 or 2,       def. 1                              */
/*                  nfmt = ...,         def. 4                              */
/*                  fmt = ...,          def. 10.4                           */
/*                  df = ...,           second output file                  */
/*                  tab = ...,          third output file                   */
/*              ) = fname;                                                  */
/*                                                                          */
/*  The second output file will contain a list of all input nodes and, for  */
/*  each node i, the number of nodes controlled by i, and their numbers.    */
/*                                                                          */
/*  The third output file will contain a list of all nodes and, for each    */
/*  of these nodes: the internal and external node number, and the number   */
/*  of input nodes that can control the node.                               */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int gfcf(TDAContext *ctx)
{
    register int i,j,k,l,nn;
    int err,nrec,nrec1,nrec2,n,nr,nb,level,ni,nii,ne,iflag;
    double sc,tmp,tmp1;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Forward control flows. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto GFCFFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GFCFFin;
    if (gdd_tcheck(ctx, 1,2,0))      /* directed, gdd option 1 required */
        goto GFCFFin;   

    if (ctx->PMSCFlg && ctx->PMSC >= 0.0)
        sc = ctx->PMSC;
    else
        sc = 0.5;   

    printf1(ctx, "Minimum level of influence: %lg\n",sc);

    if (ctx->PMOPT > 2)
        ctx->PMOPT = 2;
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    if (alloc_ack(ctx, ctx->GD_NP))                 
        goto GFCFFin;

    iflag = 0;
    if (ctx->PMF2Def) {                      /* get node list from input file */
        if ((ni = g_getnd(ctx, 0,ctx->AcK)) < 0)         
            goto GFCFFin;
        printf1(ctx, "Number of input nodes: %d\n",ni);
        iflag = 1;
    }
    if (alloc_acn(ctx, ctx->GD_NP))                 
        goto GFCFFin;

    if (alloc_acm(ctx, ctx->GD_NP))                 
        goto GFCFFin;

    if (alloc_acx(ctx, ctx->GD_NP))                 
        goto GFCFFin;

    if (ctx->PMTabFDef) {
        if (alloc_acr(ctx, ctx->GD_NP))                 
            goto GFCFFin;
    }
    nii = nrec = nrec1 = nrec2 = 0;
    ctx->GSCON_ID = 1;

    for (i = 0; i < ctx->GD_NP; ++i) {       /* loop throuh all requested nodes */

        if (iflag && ctx->AcK[i] == 0)
            continue;

        n_info(ctx, i + 1,1);
   
        for (j = 0; j < ctx->GD_NP; ++j) {
            ctx->AcN[j] = 0;
            ctx->AcM[j] = -1;
        }
        visit3(ctx, i);

        nr = -1;
        for (j = 0; j < ctx->GD_NP; ++j) {
            if (ctx->AcN[j])
                nr++;
        }

        ctx->AcM[i] = 0;
        level = 1;
        nb = 0;
        while (1) {
            n = 0;
            for (j = 0; j < ctx->GD_NP; ++j) {
                if (ctx->AcN[j] && ctx->AcM[j] < 0) {
                    tmp = 0.0;
                    nn = ctx->GD_BPN[j];
                    if (nn > 0) {
                        for (l = 0; l < nn; ++l) {
                            k = ctx->GD_BPI[j][l];          /* edge from k to j */
                            if (ctx->AcM[k] >= 0) {
                                ne = ctx->GD_BPK[j][l];
                                if (ne >= 0 && (tmp1 = gdd_ev(ctx, ne,ctx->PMGN)) > 0.0)    
                                    tmp += tmp1;
                            }
                        }
                    }
                    ctx->AcX[j] = tmp;
                    if (tmp > sc) {
                        ctx->AcM[j] = -2;
                        n++;
                    }
                }
            }
            if (n == 0)
                break;

            for (j = 0; j < ctx->GD_NP; ++j) {
                if (ctx->AcM[j] == -2) {
                    ctx->AcM[j] = level;
                    nb++;
                }
            }
            level++;
        }

        if (ctx->PMF1Def) {
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,i + 1);           
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,gdd_node(ctx, i));           
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,nb);           
        }
        if (nb > 0) {
            if (ctx->PMOPT == 2) {
                nr = 0;
                for (j = 0; j < ctx->GD_NP; ++j) {
                    if (ctx->AcM[j] > 0)
                        nr++;
                }
            }
            k = 0;
            for (j = 0; j < ctx->GD_NP; ++j) {
                if (ctx->AcN[j] && j != i) {
                    if (ctx->PMOPT == 2 && ctx->AcM[j] < 1)
                        continue;

                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);           
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));           
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, j)); 
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcM[j]);           
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcX[j]);           
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nb);           
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nr);           
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,++k);           
                    fprintf(ctx->PMFd,"\n");           
                    nrec++;

                    if (ctx->AcM[j] >= 0) {
                        if (ctx->PMF1Def)  
                            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,gdd_node(ctx, j));           

                        if (ctx->PMTabFDef)  
                            ctx->AcR[j] += 1;
                    }
                }
            }
            if (k)
                nii++;
        }
        if (ctx->PMF1Def) {
            fprintf(ctx->PMF1d,"\n");           
            nrec1++;
        }
    }
    n_info_e(ctx);

    if (ctx->PMTabFDef) {
        for (i = 0; i < ctx->GD_NP; ++i) {
            rt_fprintf_i(ctx, ctx->PMTabFd,ctx->PMNFmtS,i + 1);           
            rt_fprintf_i(ctx, ctx->PMTabFd,ctx->PMNFmtS,gdd_node(ctx, i));           
            rt_fprintf_i(ctx, ctx->PMTabFd,ctx->PMNFmtS,ctx->AcR[i]);           
            fprintf(ctx->PMTabFd,"\n");           
            nrec2++;
        }
    }
    printf1(ctx, "Number of nodes controlling at least one other node: %d\n",nii);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    if (nrec1 > 0)
        printf1(ctx, "%d records written to: %s\n",nrec1,ctx->PMF1dName);
    if (nrec2 > 0)
        printf1(ctx, "%d records written to: %s\n",nrec2,ctx->PMTabFName);
    err = 0;

GFCFFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gbcf        Backward control flows.                                     */
/*              Graph type is always directed (unvalued), gdd option 1.     */
/*                                                                          */
/*              gbcf(                                                       */
/*                  gn = graph number,                                      */
/*                  if = input file with list of nodes,                     */
/*                  sc = ...,           def. 0.5                            */
/*                  opt = 1 or 2,       def. 1                              */
/*                  nfmt = ...,         def. 4                              */
/*                  fmt = ...,          def. 10.4                           */
/*                  df = ...,           second output file                  */
/*                  ns = ...,           max number of controlled nodes,     */
/*                                      def. 2 * NP                         */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gbcf(TDAContext *ctx)
{
    register int i,j,k,l,nn;
    int err,nrec,n,nb,level,ni,nni,nnic,fnd,nd,ndd,cptr,ne;   
    double sc,tmp,tmp1;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Backward control flows. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto GBCFFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GBCFFin;
    if (gdd_tcheck(ctx, 1,2,0))      /* directed, gdd option 1 required */
        goto GBCFFin;   

    if (ctx->PMOPT > 2)
        ctx->PMOPT = 2;
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    if (ctx->PMSCFlg && ctx->PMSC >= 0.0)
        sc = ctx->PMSC;
    else
        sc = 0.5;  

    if (ctx->PMNS < 2 * ctx->GD_NP)
        ctx->PMNS = 2 * ctx->GD_NP;

    printf1(ctx, "Max number of controlled nodes: %d\n",ctx->PMNS);
    printf1(ctx, "Minimum level of influence: %lg\n",sc);

    ni = 0;
    if (ctx->PMF2Def) {
        if (alloc_ack(ctx, ctx->GD_NP))                 
            goto GBCFFin;
    
        if ((ni = g_getnd(ctx, 0,ctx->AcK)) <= 0)         /* get node list */
            goto GBCFFin;

        if (ni > 0)
            printf1(ctx, "Number of input nodes: %d\n",ni);
    }
    if (alloc_acj(ctx, ctx->PMNS))        /* list of nodes controlled by i */
        goto GBCFFin;

    if (alloc_aci(ctx, ctx->PMNS))        /* corresponding controllers */
        goto GBCFFin;

    if (alloc_acns(ctx, ctx->GD_NP))      /* number of nodes controlled by i */
        goto GBCFFin;

    if (alloc_acn(ctx, ctx->GD_NP))                 
        goto GBCFFin;

    if (alloc_acm(ctx, ctx->GD_NP))                 
        goto GBCFFin;

    /* first make list of all nodes controlled by an ultimate node */

    ctx->GSCON_ID = 1;
    nrec = nnic = nni = cptr = 0;

    for (i = 0; i < ctx->GD_NP; ++i) {       /* loop throuh all nodes */
   
        n = g_pre(ctx, i,ctx->PMGN);
        if (n > 0)                      /* use only nodes with zero in-degree */
            continue;

        nni++;
        n_info(ctx, i + 1,1);

        for (j = 0; j < ctx->GD_NP; ++j) {
            ctx->AcN[j] = 0;
            ctx->AcM[j] = -1;
        }
        visit3(ctx, i);

        ctx->AcM[i] = 0;
        level = 1;
        nb = 0;
        while (1) {
            n = 0;
            for (j = 0; j < ctx->GD_NP; ++j) {
                if (ctx->AcN[j] && ctx->AcM[j] < 0) {
                    tmp = 0.0;

                    nn = ctx->GD_BPN[j];
                    if (nn > 0) {
                        for (l = 0; l < nn; ++l) {
                            k = ctx->GD_BPI[j][l];          /* edge from k to j */
                            if (ctx->AcM[k] >= 0) {
                                ne = ctx->GD_BPK[j][l];
                                if (ne >= 0 && (tmp1 = gdd_ev(ctx, ne,ctx->PMGN)) > 0.0)    
                                    tmp += tmp1;
                            }
                        }
                    }

/**
                    k = GD_JPtr[j];
                    if (k >= 0) {
                        nn = GD_JEN[j];
                        while (nn-- > 0) {
                            if (AcM[GD_JEI[k]] >= 0) {
                                  tmp1 = get_data(ctx, GD_EV[0],GD_JE[k]);
                                  if (tmp1 > 0)   
                                      tmp += tmp1;
                            }
                            k++;
                        }
                    }
**/

                    if (tmp > sc) {
                        ctx->AcM[j] = -2;
                        n++;
                    }
                }
            }
            if (n == 0)
                break;

            for (j = 0; j < ctx->GD_NP; ++j) {
                if (ctx->AcM[j] == -2) {
                    ctx->AcM[j] = level;
                    nb++;
                }
            }
            level++;
        }
        ctx->AcNS[i] = (short)(nb);
        if (nb > 0) {
            k = 0;
            for (j = 0; j < ctx->GD_NP; ++j) {
                if (ctx->AcN[j] && ctx->AcM[j] > 0 && j != i) {
                    if (cptr >= ctx->PMNS) {
                        printf1(ctx, "Error: exceeded max number of controlled nodes (adjust ns).\n");
                        goto GBCFFin;
                    }
                    ctx->AcI[cptr] = i;
                    ctx->AcJ[cptr++] = j;
                }
            }   
        }
        if (ctx->PMF1Def) {
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,i + 1);        
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,gdd_node(ctx, i));     
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,nb);        
            fprintf(ctx->PMF1d,"\n");
            nrec++;
        }
        if (nb > 0)
            nnic++;
    }       
    n_info_e(ctx);

    printf1(ctx, "\nNumber of ultimate nodes: %d\n",nni);
    printf1(ctx, "Controlling at least one other node: %d\n",nnic);
    if (ctx->PMF1Def)  
        printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);
    if (nnic == 0)  
        goto GBCFFin;

    if (alloc_acx(ctx, ctx->GD_BPMax))                 
        goto GBCFFin;

    if (alloc_acy(ctx, ctx->GD_BPMax))                 
        goto GBCFFin;

    if (alloc_acm(ctx, ctx->GD_BPMax))                 
        goto GBCFFin;

    if (alloc_acn(ctx, ctx->GD_BPMax))                 
        goto GBCFFin;

    if (alloc_acr(ctx, ctx->GD_BPMax))                 
        goto GBCFFin;

    nni = nrec = 0;
    for (i = 0; i < ctx->GD_NP; ++i) {       /* loop through all requested nodes */

        if (ni > 0 && ctx->AcK[i] == 0)
            continue;
    
        n_info(ctx, i + 1,1);

        /* find direct backward links in AcY */

        nd = 0;
        for (j = 0; j < ctx->GD_BPMax; ++j) {
            ctx->AcX[j] = ctx->AcY[j] = 0.0;
            ctx->AcN[j] = ctx->AcM[j] = -1;
        }

        nn = ctx->GD_BPN[i];
        if (nn > 0) {
            for (l = 0; l < nn; ++l) {
                k = ctx->GD_BPI[i][l];          /* edge from k to i */
                ne = ctx->GD_BPK[i][l];
                if (ne >= 0 && (tmp = gdd_ev(ctx, ne,ctx->PMGN)) > 0.0) {  
                    ctx->AcM[nd] = k;
                    ctx->AcY[nd++] = tmp;
                }
            }
        }

/**
        k = GD_JPtr[i];
        if (k >= 0) {
            for (j = 0; j < GD_JEN[i]; ++j) {
                tmp = gdd_adj(ctx, GD_JEI[k],i,PMGN);
                if (tmp > 0.0) {
                    AcM[nd] = GD_JEI[k];
                    AcY[nd++] = tmp;
                }
                k++;
            }
        }
**/

        if (nd == 0) {
            nni++;
            continue;
        }
        for (j = 0; j < nd; ++j) {
            k = ctx->AcM[j];
            fnd = 0;
            for (l = 0; l < cptr; ++l) {
                if (ctx->AcJ[l] == k) {
                    ctx->AcN[j] = ctx->AcI[l];
                    fnd = 1;
                    break;
                }
            }
            if (fnd == 0)
                ctx->AcN[j] = k;
        }
        ndd = 0;
        for (j = 0; j < nd; ++j) {
            if (ctx->AcN[j] >= 0) {
                ctx->AcR[ndd] = ctx->AcN[j];
                for (k = j; k < nd; ++k) {
                    if (ctx->AcN[k] == ctx->AcR[ndd]) {
                        ctx->AcX[ndd] += ctx->AcY[k];
                        ctx->AcN[k] = -1;
                    }
                }
                ndd++;
            }
        }
        if (ctx->PMOPT == 1) {
            for (j = 0; j < nd; ++j) {

                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);        
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nd);        
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ctx->AcM[j]));        
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcY[j]);        
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ndd);        
                if (j < ndd) {
                    k = gdd_node(ctx, ctx->AcR[j]);
                    tmp = ctx->AcX[j];
                }
                else {
                    k = -1;
                    tmp = 0.0;
                }
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,k);        
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);        
                fprintf(ctx->PMFd,"\n");
                nrec++;
            }
        }
        else {
            
            if (sortd(ctx, nd,ctx->AcY,1))
                break;
            
            if (sortd(ctx, ndd,ctx->AcX,1))
                break;
            
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);        
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GD_BPMax);        
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nd);        
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ndd);        
            for (j = 0; j < ctx->GD_BPMax; ++j)
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcY[j]);                
            for (j = 0; j < ctx->GD_BPMax; ++j)
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcX[j]);                
            fprintf(ctx->PMFd,"\n");                
            nrec++;
        }
    }
    n_info_e(ctx);

    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    if (nni > 0)
        printf1(ctx, "Skipped %d input nodes having zero in-degree.\n",nni);

    err = 0;

GBCFFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  visit3      called by gfcf(), gbcf()                                    */

void visit3(TDAContext *ctx, int k)
{
    register int j,l,kk,nn;
   
    ctx->AcN[k] = ctx->GSCON_ID;

    nn = ctx->GD_FPN[k];
    if (nn > 0) {
        for (l = 0; l < nn; ++l) {
            j = ctx->GD_FPI[k][l];                    /* edge from k to j */
            if (ctx->AcN[j] == 0) {
                kk = ctx->GD_FPK[k][l];
                if (kk >= 0 && gdd_ev(ctx, kk,ctx->PMGN) >= 0.0)
                    visit3(ctx, j);
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  gep         Enumeration of paths.   gdd option 1 required.              */
/*                                                                          */
/*              gep(                                                        */
/*                  opt = ...,  1 = i,j, # of path, min length              */
/*                              2 = i,j, one rec for each path              */
/*                              3 = adj matrix with min values              */
/*                              4 = adj matrix with max values              */
/*                              5 = same as 1 plus location of nodes        */
/*                              6 = all plus location of nodes              */
/*                  max=...,    max number of paths, def. 100               */
/*                  gn = graph number,                                      */
/*                  nfmt = ..., integer print format, def. 4                */
/*                  fmt = ...,  print format, def. 10.4                     */
/*                  if = ...,   input file with edges                       */
/*                  ns = ...,   max number of edges, def. 100               */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gep(TDAContext *ctx)
{
    register int i,j,k,l,kk;
    int err,opt,opt1,ne,ii,maxp,mpi = 0,mpj = 0;
    double tmp;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Enumeration of paths. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto GEPFin;

    if (ctx->PMOPT > 6)
        ctx->PMOPT = 6;
    if (ctx->PMGN < 1)
        ctx->PMGN = 1;
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);
    if (ctx->PMMax < 1)
        ctx->PMMax = 100;
    if (ctx->PMNS < 1)
        ctx->PMNS = 100;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GEPFin;
    if (gdd_tcheck(ctx, 1,0,0))      /* gdd option 1 required */
        goto GEPFin;

    if (alloc_acn(ctx, ctx->GD_NP + 2))                 
        goto GEPFin;
    if (alloc_acc(ctx, ctx->GD_NP))                 
        goto GEPFin;

    ne = 0;
    if (ctx->PMF2Def) {
        if (alloc_acr(ctx, ctx->PMNS))                 
            goto GEPFin;
        if (alloc_acs(ctx, ctx->PMNS))                 
            goto GEPFin;

        if ((ne = g_getne(ctx, ctx->PMNS,ctx->AcR,ctx->AcS)) < 0)         /* get edge list */
            goto GEPFin;

        printf1(ctx, "Number of edges read from %s: %d\n",ctx->PMF2dName,ne);
        if (ne < 1)
            goto GEPFin;
    }
    opt1 = opt = 0;
    if (ctx->PMOPT == 2)
        opt = 1;
    else if (ctx->PMOPT == 5) {
        opt1 = 1;
        if (alloc_acm(ctx, ctx->GD_NP))                 
            goto GEPFin;
    }
    else if (ctx->PMOPT == 6) {
        printf1(ctx, "Max number of paths: %d\n",ctx->PMMax);
        opt1 = 2;
        if (alloc_acm(ctx, ctx->PMMax * ctx->GD_NP))       /* for internal nodes */
            goto GEPFin;
        if (alloc_acx(ctx, ctx->PMMax))               /* for values */
            goto GEPFin;
        if (alloc_acj(ctx, ctx->PMMax))               /* for length */
            goto GEPFin;
        if (alloc_aci(ctx, ctx->PMMax))               /* pointer for sorting */
            goto GEPFin;
    }
    maxp = 0;
    ctx->GSCON_NN = 0;
    ii = i = 0;
GEPCONTI:
    if (ne > 0) {
        i = ctx->AcR[ii];
        j = ctx->AcS[ii];
    }     

    if (ctx->PMOPT == 3 || ctx->PMOPT == 4) {
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));    
    }
    if (ne <= 0)
        j = 0;

GEPCONTJ:
    if (ctx->PMOPT == 3 || ctx->PMOPT == 4) {
        if (i == j) {
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,0.0);
            goto GEPCONT;        
        }
    }
    else if (j == i)  
        goto GEPCONT;
         
    else if (ne <= 0 && j < i && ctx->GD_GT <= 2)            
        goto GEPCONT;
           

    for (k = 0; k < ctx->GD_NP; ++k) {
        ctx->AcC[k] = 0;
        if (opt1 == 1)
            ctx->AcM[k] = 0;
    }

    ctx->GSCON_NP = 0;           /* # of paths */
    ctx->GSCON_NSP = 0;          /* # of shortest paths */
    ctx->GSCON_LEN = -1;         /* length of actual paths */
    ctx->GSCON_MINLEN = 0;       /* length of shortest paths */
    ctx->GSCON_MAXLEN = 0;       /* length of longest paths */
    ctx->GSCON_VAL = 0.0;        /* value of actual paths */
    ctx->GSCON_MINVAL = 0.0;     /* min value of paths */
    ctx->GSCON_MAXVAL = 0.0;     /* max value of paths */
    ctx->GSCON_CNT = 0;
    visit4(ctx, i,i,j,0.0,1,opt,opt1);

    if (maxp < ctx->GSCON_NP) {
        maxp = ctx->GSCON_NP;
        mpi  = i;
        mpj  = j;
    }
    tmp = -1.0;
    if (ctx->PMOPT == 1 || opt1 == 1) {
        if (ctx->GSCON_NP > 0) {
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j + 1);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, j));
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_NP);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_NSP);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_MINLEN);
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->GSCON_MINVAL);
            if (opt1 == 0) {
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_MAXLEN);
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->GSCON_MAXVAL);
            }   
            else {
                for (k = 0; k < ctx->GD_NP; ++k)
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcM[k]);
            }
            fprintf(ctx->PMFd,"\n");
            ctx->GSCON_NN++;
        }
    }
    else if (ctx->PMOPT == 3) {
        if (ctx->GSCON_NP > 0)
            tmp = ctx->GSCON_MINVAL;
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
    }
    else if (ctx->PMOPT == 4) {
        if (ctx->GSCON_NP > 0)
             tmp = ctx->GSCON_MAXVAL;
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
    }
    else if (ctx->PMOPT == 6) {
        if (ctx->GSCON_NP > ctx->PMMax) {
            printf1(ctx, "Error: exceeded max number of paths (%d,%d).\n",
                                                   gdd_node(ctx, i),gdd_node(ctx, j));
            goto GEPFin;
        }
        if (sortdp(ctx, ctx->GSCON_NP,ctx->AcX,ctx->AcI))
            goto GEPFin;

        for (k = 0; k <  ctx->GSCON_NP; ++k) {
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j + 1);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, j));
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_NP);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,k + 1);

            kk = ctx->AcI[k];
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcJ[kk]);
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcX[kk]);
            kk *= ctx->GD_NP;
            for (l = 0; l < ctx->GD_NP; ++l)
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcM[kk + l]);
            fprintf(ctx->PMFd,"\n");
            ctx->GSCON_NN++;
        }
    }

GEPCONT:
    if (ne > 0) {
        if (ctx->PMOPT == 3 || ctx->PMOPT == 4) {
            fprintf(ctx->PMFd,"\n");
            ctx->GSCON_NN++;
        }
        if (++ii < ne)
            goto GEPCONTI;
    }
    else {
        if (++j < ctx->GD_NP)
            goto GEPCONTJ;
        if (ctx->PMOPT == 3 || ctx->PMOPT == 4) {
            fprintf(ctx->PMFd,"\n");
            ctx->GSCON_NN++;
        }
        n_info(ctx, i + 1,1);
        if (++i < ctx->GD_NP)
            goto GEPCONTI;
    }
    n_info_e(ctx); 
    printf1(ctx, "Max number of paths (nodes: %d,%d): %d\n",
                                         gdd_node(ctx, mpi),gdd_node(ctx, mpj),maxp);
    printf1(ctx, "%d records written to: %s\n",ctx->GSCON_NN,ctx->PMFdName);
    err = 0;

GEPFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  visit4      called by gep                                               */
/*              opt  = 1  print for (i,j) pairs (paths)                     */
/*              opt  = 2  print for nodes (cycles)                          */
/*              opt1 = 1  count how often nodes occur in paths, in AcM[]    */
/*              opt1 = 2  save up to PMMax paths in AcM, values in AcX      */
/*                        length in AcJ.                                    */
/*                        count number of paths in GSCON_NP.                */

void visit4(TDAContext *ctx, int k,int k1,int k2,double val,int first,int opt,int opt1)
{
    register int j,l,nn;
    double tmp;

    /* printf1("visit4 k=%d k1=%d k2=%d first=%d\n",k,k1,k2,first); */

    if (opt1 == 2 && ctx->GSCON_NP >= ctx->PMMax) {
        ctx->GSCON_NP = ctx->PMMax + 1;
        ctx->AcC[k] = 1;
        return;
    }
    if (first == 0 || k != k1 || k != k2)  
        ctx->AcC[k] = 1;         
    ctx->AcN[ctx->GSCON_CNT++] = k;
    ctx->GSCON_LEN++;
    ctx->GSCON_VAL += val;
       
    if (k == k2 && first == 0) {

        if (k1 == k2 && ctx->GSCON_LEN <= 2)     /* for cycles min length 3 */
            goto V4CONT;

        if (opt1 == 2) {
            ctx->AcX[ctx->GSCON_NP] = ctx->GSCON_VAL;
            ctx->AcJ[ctx->GSCON_NP] = ctx->GSCON_LEN;
            nn = ctx->GSCON_NP * ctx->GD_NP;
            for (j = 0; j < ctx->GD_NP; ++j)  
                ctx->AcM[nn + j] = 0;
            for (j = 1; j < ctx->GSCON_CNT - 1; ++j)  
                ctx->AcM[nn + ctx->AcN[j]] = 1;
        }

        ctx->GSCON_NP++;
        if (ctx->GSCON_MINVAL <= 0.0 || ctx->GSCON_MINVAL > ctx->GSCON_VAL) {
            ctx->GSCON_MINVAL = ctx->GSCON_VAL;
            ctx->GSCON_MINLEN = ctx->GSCON_LEN;
            ctx->GSCON_NSP = 1;

            if (opt1 == 1) {          
                for (j = 0; j < ctx->GD_NP; ++j)  
                    ctx->AcM[j] = 0;
                for (j = 1; j < ctx->GSCON_CNT - 1; ++j)  
                    ctx->AcM[ctx->AcN[j]] = 1;
            }      
        }
        else if (fabs(ctx->GSCON_MINVAL - ctx->GSCON_VAL) < ctx->EPSI1) {
            ctx->GSCON_NSP++;
            if (ctx->GSCON_MINLEN > ctx->GSCON_LEN)
                ctx->GSCON_MINLEN = ctx->GSCON_LEN;

            if (opt1 == 1) {          
                for (j = 1; j < ctx->GSCON_CNT - 1; ++j)  
                    ctx->AcM[ctx->AcN[j]] += 1;
            }      
        }
        if (ctx->GSCON_MAXVAL <= 0.0 || ctx->GSCON_MAXVAL < ctx->GSCON_VAL) {
            ctx->GSCON_MAXVAL = ctx->GSCON_VAL;
            ctx->GSCON_MAXLEN = ctx->GSCON_LEN;
        }
        else if (fabs(ctx->GSCON_MAXVAL - ctx->GSCON_VAL) < ctx->EPSI1) {
            if (ctx->GSCON_MAXLEN < ctx->GSCON_LEN)
                ctx->GSCON_MAXLEN = ctx->GSCON_LEN;
        }

        if (opt) {
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,k1 + 1);
            if (opt == 1)
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,k2 + 1);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, k1));
            if (opt == 1)
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, k2));

            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_NP);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_LEN);
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->GSCON_VAL);
            for (j = 0; j < ctx->GSCON_CNT; ++j)  
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ctx->AcN[j]));
            fprintf(ctx->PMFd,"\n");
            ctx->GSCON_NN++;
        }   
        
V4CONT:
        ctx->GSCON_LEN--;      
        ctx->GSCON_VAL -= val;
        ctx->AcC[k] = 0;
        ctx->GSCON_CNT--;             
        return;
    }
    nn = ctx->GD_FPN[k];
    if (nn > 0) {  
        for (l = 0; l < nn; ++l) {
            j = ctx->GD_FPI[k][l];                    /* edge from k to j */
            if (ctx->AcC[j] == 0) {
                tmp = gdd_adj(ctx, k,j,ctx->PMGN);
                if (tmp >= 0.0)  
                    visit4(ctx, j,k1,k2,tmp,0,opt,opt1);
            }
        }
    }                  
    ctx->GSCON_LEN--;      
    ctx->GSCON_VAL -= val;
    if (ctx->GSCON_CNT <= 0) {
        printfe(ctx, "ERROR in visit4.\n");
        gerr_exit(ctx, 216);
    }
    ctx->AcC[k] = 0;
    ctx->GSCON_CNT--;
}

/* ------------------------------------------------------------------------ */
/*  gflow       Maximal flows.                                              */
/*                                                                          */
/*              gflow(                                                      */
/*                  opt = ...,      1 = one record for each (i,j)           */
/*                                  2 = square matrix                       */
/*                  gn = ...,       graph number,                           */
/*                  nfmt = ...,     integer format, def. 4                  */
/*                  fmt=...,        print format for values, def. 10.4      */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gflow(TDAContext *ctx)
{
    register int i,j;
    int err,n;
    double f;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Maximal flows. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 5,1,1))       /* get parameters */
        goto GFLOWFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GFLOWFin;

    if (gdd_tcheck(ctx, 1,2,0))      /* must be directed, gdd option 1 */
        goto GFLOWFin;

    if (ctx->PMOPT > 2)
        ctx->PMOPT = 2;
    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 10,4);

    if (alloc_acn(ctx, ctx->GD_NP))       /* used as queue in g_flowbfs(ctx) */
        goto GFLOWFin;
    if (alloc_aci(ctx, ctx->GD_NP))       /* used for pre[] in g_flowbfs */ 
        goto GFLOWFin;
    if (alloc_acc(ctx, ctx->GD_NP))       /* used for vis[] in g_flowbfs */
        goto GFLOWFin;
    if (alloc_acj(ctx, ctx->GD_NP))       /* used for pointers to flow[] */
        goto GFLOWFin;
    if (alloc_acz(ctx, ctx->GD_NP))       /* used for fd in g_flowbfs */
        goto GFLOWFin;

    n = 0;
    for (i = 0; i < ctx->GD_NP; ++i) {

        if (ctx->SILENTFlg < 2 && ctx->GD_NP >= 500) {
            printfe(ctx, "Node: %7d     %c",i + 1,CR);
            fflushe(ctx);
        }
        /* the row prefix (index, node) belongs to the matrix form,
           opt=2, where each source node has one line.  With opt=1
           each pair is its own complete record, and the prefix put
           two stray numbers in front of a source node's first record
           (and left them dangling for a node with no record at all);
           the manual's Box 1 of 7.2.8.1 shows the six-column record
           without them. */
        if (ctx->PMOPT == 2) {
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));
        }

        for (j = 0; j < ctx->GD_NP; ++j) {

            if (j == i) {
                if (ctx->PMOPT == 2)
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,0.0);
                continue;
            }

            if (alloc_acx(ctx, ctx->GD_NOC))      /* used for flow[] */
                goto GFLOWFin;

            f = g_flow(ctx, ctx->GD_NP,i,j,ctx->PMGN,ctx->AcN,ctx->AcI,ctx->AcC,ctx->AcJ,ctx->AcX,ctx->AcZ);

            if (ctx->PMOPT == 1) {
                if (f < 0.0)
                    continue;

                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,++n);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j + 1);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, j));
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,f);
                fprintf(ctx->PMFd,"\n");
            }
            else  
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,f);

        }
        if (ctx->PMOPT == 2) {
            fprintf(ctx->PMFd,"\n");
            n++;
        }
    }
    if (ctx->SILENTFlg < 2 && ctx->GD_NP >= 500) {
        printfe(ctx, "\n");
        fflushe(ctx);
    }
    printf1(ctx, "%d records written to: %s\n",n,ctx->PMFdName);
    err = 0;

GFLOWFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_flow.     Find maximal flow between nodes i0 and i1.                  */
/*              Uses the Edmonds and Karp algorithm as described in         */
/*              V. Turau, Algorithmische Graphentheorie, Addison-Wesley,    */  
/*              1996, pp. 171ff.                                            */
/*                                                                          */
/*  Return value of max flow, or -1 if there is no path from i0 to i1.      */

double g_flow(TDAContext *ctx, int n,int i0,int i1,int gn,int *q,int *pre,char *vis,int *fptr,   double *flow,double *fd)
{
    register int i,k,first;
    double fdd,fmax;        

    first = 1;
    fmax = 0.0;
    while (1) {
        fdd = g_flowbfs(ctx, n,i0,i1,gn,q,pre,vis,fptr,flow,fd);               
        if (first && fdd < 0.0)
            return(-1.0);
        else if (fdd <= 0.0)
            break;
        first = 0;
        fmax += fdd;
        i = i1;                     /* modify current flow */
        while (pre[i] >= 0) {
            k = fptr[i];
            if (k > 0)  
                flow[k - 1] += fdd;
            else  
                flow[-k -1] -= fdd;
            i = pre[i];
        }
    }
    return(fmax);
}

/* ------------------------------------------------------------------------ */
/*  g_flowbfs()     Used by g_flow. Tries to find a new flow from i0 to     */
/*                  i1 with higher value. Uses breadth-first search.        */
/*                  Return increment of current flow, or -1 if there        */
/*                  is no path from i0 to i1.                               */

double g_flowbfs(TDAContext *ctx, int n,int i0,int i1,int gn,int *q,int *pre,char *vis,int *fptr,   double *flow,double *fd)
{
    register int i,j,l,k;
    int aptr,bptr,m,nf;
    double tmp;

    for (i = 0; i < n; ++i)
        vis[i] = 0;

    nf = 1;
    pre[i0] = -1;  
    vis[i0] = 1;
    fd[i0] = ctx->DBLMAX;
    fd[i1] = 0.0;
    q[0] = i0;
    aptr = 0;
    bptr = 1;

    while (aptr < bptr) {
        i = q[aptr++];
        m = ctx->GD_FPN[i];
        for (l = 0; l < m; ++l) {
            j = ctx->GD_FPI[i][l];           /* edge from i to j */
            if (vis[j] == 0) {
                k = ctx->GD_FPK[i][l];           /* pointer to edge and flow */
                tmp = gdd_ev(ctx, k,gn);   
                if (tmp >= 0.0) {
                    if (j == i1)
                        nf = 0;
                    if (flow[k] < tmp) {
                        vis[j] = 1;     
                        pre[j] = i;
                        fptr[j] = k + 1;
                        q[bptr++] = j;
                        fd[j] = dmin(ctx, fd[i],tmp - flow[k]);
                    }
                }
            }
        }
        m = ctx->GD_BPN[i];
        for (l = 0; l < m; ++l) {
            j = ctx->GD_BPI[i][l];           /* edge from j to i */
            if (vis[j] == 0) {
                k = ctx->GD_BPK[i][l];       /* pointer to edge and flow */
                tmp = gdd_ev(ctx, k,gn);   
                if (tmp >= 0.0) {
                    if (flow[k] > 0.0) {
                        vis[j] = 1;     
                        pre[j] = i;
                        fptr[j] = -(k + 1);
                        q[bptr++] = j;
                        fd[j] = dmin(ctx, fd[i],flow[k]);
                    }
                }
            }
        }
        if (fd[i1] > 0.0)
            break;
    }
    if (nf)
        return(-1.0);
    return(fd[i1]);
}

/* ------------------------------------------------------------------------ */
/*  gfc         Flow control in directed valued graphs.                     */
/*                                                                          */
/*              gfc(                                                        */
/*                  gn = ...,       graph number,                           */
/*                  nfmt = ...,     integer format, def. 4                  */
/*                  fmt=...,        print format for values, def. 10.4      */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gfc(TDAContext *ctx)
{
    register int i,j,k;
    int err,n;
    double f,f0;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Flow control. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto GFCFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GFCFin;

    if (gdd_tcheck(ctx, 1,2,0))      /* must be directed, gdd option 1 */
        goto GFCFin;

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 10,4);

    if (alloc_acn(ctx, ctx->GD_NP))       /* used as queue in g_flowbfs(ctx) */
        goto GFCFin;
    if (alloc_aci(ctx, ctx->GD_NP))       /* used for pre[] in g_flowbfs */ 
        goto GFCFin;
    if (alloc_acc(ctx, ctx->GD_NP))       /* used for vis[] in g_flowbfs */
        goto GFCFin;
    if (alloc_acj(ctx, ctx->GD_NP))       /* used for pointers to flow[] */
        goto GFCFin;
    if (alloc_acz(ctx, ctx->GD_NP))       /* used for fd in g_flowbfs */
        goto GFCFin;

    n = 0;
    for (i = 0; i < ctx->GD_NP; ++i) {

        if (ctx->SILENTFlg < 2 && ctx->GD_NP >= 500) {
            printfe(ctx, "Node: %7d     %c",i + 1,CR);
            fflushe(ctx);
        }
        for (j = 0; j < ctx->GD_NP; ++j) {
            if (j == i)  
                continue;

            /* first calculate max flow from i to j */

            if (alloc_acx(ctx, ctx->GD_NOC))      /* used for flow[] */
                goto GFCFin;

            f0 = g_flow(ctx, ctx->GD_NP,i,j,ctx->PMGN,ctx->AcN,ctx->AcI,ctx->AcC,ctx->AcJ,ctx->AcX,ctx->AcZ);

            if (f0 < 0.0)
                continue;

            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j + 1);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, j));
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,f0);

            for (k = 0; k < ctx->GD_NP; ++k) {

                f = -1.0;
                if (k != i && k != j) {
                                    
                    if (alloc_ace(ctx, ctx->GD_NOC))  
                        goto GFCFin;
                                    
                    if (fc_visit(ctx, i,k,ctx->PMGN)) {      

                        if (alloc_ace(ctx, ctx->GD_NOC))  
                            goto GFCFin;
                                    
                        if (fc_visit(ctx, k,j,ctx->PMGN)) {      

                            if (alloc_acx(ctx, ctx->GD_NOC))      /* used for flow[] */
                                goto GFCFin;

                            f = g_cflow(ctx, ctx->GD_NP,i,j,ctx->PMGN,ctx->AcN,ctx->AcI,ctx->AcC,ctx->AcJ,ctx->AcX,ctx->AcZ,k);
                            if (f < 0.0)
                                f = 0.0;
                        }   
                    }
                }
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,f);
            }
            fprintf(ctx->PMFd,"\n");
            n++;
        }
    }
    if (ctx->SILENTFlg < 2 && ctx->GD_NP >= 500) {
        printfe(ctx, "\n");
        fflushe(ctx);
    }
    printf1(ctx, "%d records written to: %s\n",n,ctx->PMFdName);
    err = 0;

GFCFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_cflow.    Find maximal flow between nodes i0 and i1.                  */
/*              Refers to a reduced graph with node k0 deleted.             */
/*              Uses the Edmonds and Karp algorithm as described in         */
/*              V. Turau, Algorithmische Graphentheorie, Addison-Wesley,    */  
/*              1996, pp. 171ff.                                            */
/*                                                                          */
/*  Return value of max flow, or -1 if there is no path from i0 to i1.      */

double g_cflow(TDAContext *ctx, int n,int i0,int i1,int gn,int *q,int *pre,char *vis,int *fptr,   double *flow,double *fd,int k0)
{
    register int i,k,first;
    double fdd,fmax;        

    first = 1;
    fmax = 0.0;
    while (1) {
        fdd = g_cflowbfs(ctx, n,i0,i1,gn,q,pre,vis,fptr,flow,fd,k0);               
        if (first && fdd < 0.0)
            return(-1.0);
        else if (fdd <= 0.0)
            break;
        first = 0;
        fmax += fdd;
        i = i1;                     /* modify current flow */
        while (pre[i] >= 0) {
            k = fptr[i];
            if (k > 0)  
                flow[k - 1] += fdd;
            else  
                flow[-k -1] -= fdd;
            i = pre[i];
        }
    }
    return(fmax);
}

/* ------------------------------------------------------------------------ */
/*  g_cflowbfs()    Used by g_flow. Tries to find a new flow from i0 to     */
/*                  i1 with higher value. Uses breadth-first search.        */
/*                  Refers to reduced graph where node k0 is deleted.       */
/*                  Return increment of current flow, or -1 if there        */
/*                  is no path from i0 to i1.                               */

double g_cflowbfs(TDAContext *ctx, int n,int i0,int i1,int gn,int *q,int *pre,char *vis,int *fptr,   double *flow,double *fd,int k0)
{
    register int i,j,l,k;
    int aptr,bptr,m,nf;
    double tmp;

    for (i = 0; i < n; ++i)
        vis[i] = 0;

    vis[k0] = 1;
    pre[k0] = -1;

    nf = 1;
    pre[i0] = -1;  
    vis[i0] = 1;
    fd[i0] = ctx->DBLMAX;
    fd[i1] = 0.0;
    q[0] = i0;
    aptr = 0;
    bptr = 1;

    while (aptr < bptr) {
        i = q[aptr++];
        m = ctx->GD_FPN[i];
        for (l = 0; l < m; ++l) {
            j = ctx->GD_FPI[i][l];           /* edge from i to j */
            if (vis[j] == 0) {
                k = ctx->GD_FPK[i][l];           /* pointer to edge and flow */
                tmp = gdd_ev(ctx, k,gn);   
                if (tmp >= 0.0) {
                    if (j == i1)
                        nf = 0;
                    if (flow[k] < tmp) {
                        vis[j] = 1;     
                        pre[j] = i;
                        fptr[j] = k + 1;
                        q[bptr++] = j;
                        fd[j] = dmin(ctx, fd[i],tmp - flow[k]);
                    }
                }
            }
        }
        m = ctx->GD_BPN[i];
        for (l = 0; l < m; ++l) {
            j = ctx->GD_BPI[i][l];           /* edge from j to i */
            if (vis[j] == 0) {
                k = ctx->GD_BPK[i][l];       /* pointer to edge and flow */
                tmp = gdd_ev(ctx, k,gn);   
                if (tmp >= 0.0) {
                    if (flow[k] > 0.0) {
                        vis[j] = 1;     
                        pre[j] = i;
                        fptr[j] = -(k + 1);
                        q[bptr++] = j;
                        fd[j] = dmin(ctx, fd[i],flow[k]);
                    }
                }
            }
        }
        if (fd[i1] > 0.0)
            break;
    }
    if (nf)
        return(-1.0);
    return(fd[i1]);
}

/* ------------------------------------------------------------------------ */
/*  fc_visit(k,k0,gn)   Return 1 if k0 can be reached from k, otherwise 0   */

int fc_visit(TDAContext *ctx, int k,int k0,int gn)
{
    register int j,l;
    int n;

    if (k == k0)
        return(1);

    ctx->AcE[k] = 1;         

    n = ctx->GD_FPN[k];
    if (n <= 0)
        return(0);       

    for (l = 0; l < n; ++l) {

        j = ctx->GD_FPI[k][l];                    /* edge from k to j */
        if (ctx->AcE[j] == 0) {
            if (gdd_adj(ctx, k,j,gn) >= 0.0) {
                if (fc_visit(ctx, j,k0,gn))
                    return(1);
            }
        }
    }
    return(0);
}

