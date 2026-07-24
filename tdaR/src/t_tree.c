/****************************************************************************/
/*  t_tree                                                                  */
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
#include "t_gdd.h"
#include "t_plot.h"
#include "t_alloc.h"
#include "tda_context.h"

/*  functions in t_tree.c */

int gt_alloc(TDAContext *ctx, int n);
int g_tree(TDAContext *ctx, int r);
int get_root(TDAContext *ctx);
int ptree(TDAContext *ctx);
int pltree(TDAContext *ctx);
int pl_tree(TDAContext *ctx, int ptyp,int r,int lt,double lw,double fs);

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */



/*  additional tree structure */


/* ------------------------------------------------------------------------ */
/*  gt_alloc(n)     If n > 0 allocate memory for tree data structure,       */
/*                  with n nodes. If n = 0 free previously allocated        */
/*                  memory.                                                 */
/*                                                                          */
/*  Return: 0 if OK, -1 if insufficient memory                              */

int gt_alloc(TDAContext *ctx, int n)
{
    int err = 0;

    if (n == 0)  
        goto GTALFin;
       
    err = -1;

    if (!(ctx->GT_EI = (int *)calloc((size_t)(n + 1),sizeof(int))))  
        goto GTALFin; 
    memrq(ctx, n + 1,sizeof(int));
    ctx->GT_EIA = n + 1;

    if (!(ctx->GT_EJ = (int *)calloc((size_t)(n + 1),sizeof(int))))  
        goto GTALFin; 
    memrq(ctx, n + 1,sizeof(int));
    ctx->GT_EJA = n + 1;

    if (!(ctx->GT_EV = (float *)calloc((size_t)(n + 1),sizeof(float))))  
        goto GTALFin; 
    memrq(ctx, n + 1,sizeof(float));
    ctx->GT_EVA = n + 1;

    if (!(ctx->GT_PR = (int *)calloc((size_t)(n + 1),sizeof(int))))  
        goto GTALFin; 
    memrq(ctx, n + 1,sizeof(int));
    ctx->GT_PRA = n + 1;

    if (!(ctx->GT_SU = (int *)calloc((size_t)(n + 1),sizeof(int))))  
        goto GTALFin; 
    memrq(ctx, n + 1,sizeof(int));
    ctx->GT_SUA = n + 1;

    if (!(ctx->GT_FR = (int *)calloc((size_t)(n + 1),sizeof(int))))  
        goto GTALFin; 
    memrq(ctx, n + 1,sizeof(int));
    ctx->GT_FRA = n + 1;

    return(0);

GTALFin:
    if (ctx->GT_EIA > 0) {
        free((char *)ctx->GT_EI);
        memrq(ctx, -ctx->GT_EIA,sizeof(int));
        ctx->GT_EIA = 0;
    }
    if (ctx->GT_EJA > 0) {
        free((char *)ctx->GT_EJ);
        memrq(ctx, -ctx->GT_EJA,sizeof(int));
        ctx->GT_EJA = 0;
    }
    if (ctx->GT_EVA > 0) {
        free((char *)ctx->GT_EV);
        memrq(ctx, -ctx->GT_EVA,sizeof(float));
        ctx->GT_EVA = 0;
    }
    if (ctx->GT_PRA > 0) {
        free((char *)ctx->GT_PR);
        memrq(ctx, -ctx->GT_PRA,sizeof(int));
        ctx->GT_PRA = 0;
    }
    if (ctx->GT_SUA > 0) {
        free((char *)ctx->GT_SU);
        memrq(ctx, -ctx->GT_SUA,sizeof(int));
        ctx->GT_SUA = 0;
    }
    if (ctx->GT_FRA > 0) {
        free((char *)ctx->GT_FR);
        memrq(ctx, -ctx->GT_FRA,sizeof(int));
        ctx->GT_FRA = 0;
    }
    ctx->GT_RT = 0;
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_tree(r)     Create tree data structure, as proposed by                */  
/*                Barthelemy and Guenoche, Trees and Proximity              */
/*                Representations, 1988, p. 7.                              */
/*                                                                          */
/*                Based on current graph data (PMGN).                       */
/*                r = root; if r = 0 choose an arbitrary root.              */
/*                                                                          */
/*                Note: tree data structure uses internal node numbers      */
/*                1,...,GD_NP, instead of 0,...,GD_NP - 1.                  */
/*                                                                          */
/*  Return: 0 if OK, -2 if not a tree, -1 if insufficient memory            */

int g_tree(TDAContext *ctx, int r)
{
    register int i,j,u,v = 0,k;
    int err,m,n,ne,sdeg;                   
    double tmp;
               
    n = ctx->GD_NP;
            
    if (ctx->GD_NP != ctx->GD_NE[ctx->PMGN - 1] + 1) {
        err = -2;
        goto GTREEFin;
    }
              
    if (gt_alloc(ctx, n)) {
        p_err(ctx, -2,1);
        return(-1);
    }
    err = -1;
    if (alloc_acn(ctx, n + 1))
        goto GTREEFin;

    sdeg = m = 0;
    for (i = 2; i <= n; ++i) {              /* count edges and degrees */
        for (j = 1; j < i; ++j) {           /* and save edges */

            tmp = (double)((float)gdd_adj(ctx, i - 1,j - 1,ctx->PMGN));         
            if (tmp >= 0.0) {
                m++;
                if (m < n) {
                    ctx->GT_EI[m] = i;
                    ctx->GT_EJ[m] = j;
                    ctx->GT_EV[m] = (float)(tmp);
                }
                ctx->AcN[i] += 1;
                ctx->AcN[j] += 1;
                sdeg += 2;
            }
        }
    }
    if (m + 1 != n) {                   /* if no tree */
        err = -2;
        goto GTREEFin;
    }
    if (alloc_aci(ctx, n + 1))
        goto GTREEFin;

    if (alloc_acj(ctx, n + 1))
        goto GTREEFin;

    if (alloc_acxf(ctx, n + 1))
        goto GTREEFin;

    /* find a representation for any root */

    while (sdeg > 0) {

        for (i = 1; i <= n; ++i)  
            ctx->AcI[i] = 0;

        ne = 0;
        for (i = 1; i <=  n; ++i) {

            if (ctx->AcJ[i])   
                continue;

            if (ctx->AcN[ctx->GT_EI[i]] == 1) {
                u = ctx->GT_EI[i];
                v = ctx->GT_EJ[i];
            }
            else if (ctx->AcN[ctx->GT_EJ[i]] == 1) {
                u = ctx->GT_EJ[i];
                v = ctx->GT_EI[i];
            }
            else
                continue;

            ctx->GT_PR[u] = v;
            ctx->AcJ[i] = 1;
            ne++;
            ctx->AcI[u] += 1;
            ctx->AcI[v] += 1;
            ctx->AcXF[u] = ctx->GT_EV[i];
            k = ctx->GT_SU[v];
            if (k == 0) 
                ctx->GT_SU[v] = u;
            else {
                while (ctx->GT_FR[k] > 0)  
                    k = ctx->GT_FR[k];
                ctx->GT_FR[k] = u;
            }
        }
        if (ne == 0) {
            err = -2;
            goto GTREEFin;
        }
        for (i = 1; i <= n; ++i) {
            ctx->AcN[i] -= ctx->AcI[i];
            sdeg -= ctx->AcI[i];
        }
    }
    for (i = 1; i <= n; ++i)
        ctx->GT_EV[i] = ctx->AcXF[i];

    /* current root is v. if v != r find new representation for root r */

    ctx->GT_RT = v;
    if (r == 0 || v == r) {
        err = 0;
        goto GTREEFin;
    }

    v = ctx->GT_RT;  /* current root */
    u = r;      /* required root */                 
    i = 0;
    ctx->GT_EV[v] = ctx->GT_EV[u];

    while (u > 0) {
        v = ctx->GT_PR[u];
        ctx->GT_PR[u] = i;
        i = u;
        if (v == 0)
            break;    

        k = ctx->GT_SU[u];
        if (k == 0)  
            ctx->GT_SU[u] = v;
        else {
            j = ctx->GT_FR[k];
            while (j > 0) {
                k = j;
                j = ctx->GT_FR[k];
            }
            ctx->GT_FR[k] = v;
        }
        k = ctx->GT_SU[v];
        if (k == u) {
            ctx->GT_SU[v] = ctx->GT_FR[u];
            ctx->GT_FR[u] = 0;
        }
        else {
            j = ctx->GT_FR[k];
            while (j > 0 && j != u) {
                k = j;
                j = ctx->GT_FR[k];
            }
            ctx->GT_FR[k] = ctx->GT_FR[u];
            ctx->GT_FR[u] = 0;
        }
        u = v;
        if (u <= 0)
            break;

        tmp = (double)(ctx->GT_EV[r]);
        ctx->GT_EV[r] = ctx->GT_EV[v];
        ctx->GT_EV[v] = (float)(tmp);
    }
    ctx->GT_RT = r;   /* new root */
    ctx->GT_EV[r] = 0.0;
    err = 0;

GTREEFin:
    alloc_acn(ctx, 0);
    alloc_aci(ctx, 0);
    alloc_acj(ctx, 0);
    alloc_acxf(ctx, 0);
    if (err) {
        if (err == -2)
            printf1(ctx, "Error: current graph is not a tree.\n");
        gt_alloc(ctx, 0);
    }
    else
        printf1(ctx, "Selected root: %d\n",ctx->GD_ND[ctx->GT_RT - 1]);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_root()  return internal root number based on PMRT                   */   
/*              if PMRT  not specified return 0.                            */
/*              if error, message and return -1.                            */

int get_root(TDAContext *ctx)
{
    register int i;
    int rt;

    if (ctx->PMRT < 1)  
        rt = 0;
    else {
        rt = -1;
        for (i = 0; i < ctx->GD_NP; ++i) {
            if (ctx->GD_ND[i] == ctx->PMRT) {
                rt = i + 1;
                break;
            }
        }
    }
    if (rt < 0) {
        printf1(ctx, "Error: root node %d not available.\n",ctx->PMRT);
        return(-1);     
    }
    return(rt);
}

/* ------------------------------------------------------------------------ */
/*  ptree()     Print tree data structure.                                  */
/*                                                                          */
/*              ptree(                                                      */
/*                  rt=...,     root                                        */
/*                  nfmt=...,   integer format, def. 4                      */
/*                  fmt=...,    format for values, def. 10.4                */
/*                  gn=...,     graph number, def. 1                        */
/*              ) = fname                                                   */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int ptree(TDAContext *ctx)
{
    register int i;
    int err,rt;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Print tree data structure. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,1,1))       /* get parameters */
        goto PTREEFin;

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 10,4);
/**
    if (gd_setup(PMGN,PMGT,PMGTT,0,0,0,1,2))     always undirected       
        goto PTREEFin;
**/
    /*  There was no check that a graph exists -- the gd_setup() call
        that would have done it is commented out just above.  With no
        graph, g_tree() indexed GD_NE[PMGN - 1] on an unallocated array
        and read wild memory (Valgrind, ptree.cf, which runs ptree() on
        a bare session).  gdd_check reports it the way every other graph
        command does.  */
    if (gdd_check(ctx, ctx->PMGN,0))
        goto PTREEFin;

    if ((rt = get_root(ctx)) < 0)
        goto PTREEFin;
       
    if (g_tree(ctx, rt))             /* setup tree data structure */
        goto PTREEFin;

    fprintf(ctx->PMFd,"# tree data structure with root %d\n",ctx->GD_ND[ctx->GT_RT - 1]);
   
    for (i = 1; i <= ctx->GD_NP; ++i) {
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GD_ND[i - 1]);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GT_EI[i]);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GT_EJ[i]);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,(double)(ctx->GT_EV[i]));
#ifdef TDA_R_PACKAGE
        {
            double erow[8];
            erow[0] = (double)i;
            erow[1] = (double)ctx->GD_ND[i - 1];
            erow[2] = (double)ctx->GT_EI[i];
            erow[3] = (double)ctx->GT_EJ[i];
            erow[4] = (double)(ctx->GT_EV[i]);
            erow[5] = (double)ctx->GT_PR[i];
            erow[6] = (double)ctx->GT_SU[i];
            erow[7] = (double)ctx->GT_FR[i];
            tda_export_row(ctx, "ptree.table", erow, 8);
        }
#endif
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GT_PR[i]);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GT_SU[i]);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GT_FR[i]);
        fprintf(ctx->PMFd,"\n");
    }    
    printf1(ctx, "%d records written to: %s\n",ctx->GD_NP,ctx->PMFdName);
    err = 0;

PTREEFin:
    gt_alloc(ctx, 0);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  pltree()    Plot tree. Graph must be undrected.                         */
/*  ##                                                                      */  
/*              pltree(                                                     */
/*                  gn=...,     graph number, def. 1.                       */
/*                  pl=...,     type of drawing, def. 3                     */
/*                  rt=...,     root                                        */
/*                  nc=...,     label nodes up to nc                        */
/*                  ni=1,       labeL only leaves, def. 0                   */
/*                  lt=...,     line type, def. 1                           */
/*                  lw=...,     line width, def. 0.2 mm                     */
/*                  fs=...,     font size, def. 2 mm                        */
/*              );                                                          */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int pltree(TDAContext *ctx)
{
    int err,rt;

    err = -1;
    if (check_pcmd(ctx, 1,2))  /* check plot command for PostScript set up */ 
        return(-1);

    printf1(ctx, "Plot tree. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 6 ,1,0))       /* get parameters */
        goto PLTREEFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto PLTREEFin;
    if (gdd_tcheck(ctx, 0,1,2))          /* need undirected valued graph */
        goto PLTREEFin;

    if ((rt = get_root(ctx)) < 0)
        goto PLTREEFin;
     
    if (g_tree(ctx, rt))             /* setup tree data structure */
        goto PLTREEFin;

    if (ctx->PMPLFlg < 1 || ctx->PMPLFlg > 3)
        ctx->PMPLFlg = 3;

    printf1(ctx, "Type of drawing: %d ",ctx->PMPLFlg);

    if (ctx->PMPLFlg == 1) 
        printf1(ctx, "(radial)\n");
    else if (ctx->PMPLFlg == 2) 
        printf1(ctx, "(hierarchical)\n");
    else               
        printf1(ctx, "(axial)\n");

    err = pl_tree(ctx, ctx->PMPLFlg,ctx->GT_RT,ctx->PMLT,ctx->PMLW,ctx->PMFS);
    if (err) {
        if (err == -2) {
            printf1(ctx, "Unable to plot this tree, try another root.\n");
            err = 0;
        }
    }
    else  
        printf1(ctx, "PostScript output written to: %s\n",ctx->PSFName);

PLTREEFin:
    gt_alloc(ctx, 0);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  pl_tree(ptyp,r,lt,lw,fs)                                                */  
/*                                                                          */
/*                Plot tree, as proposed by                                 */  
/*                Barthelemy and Guenoche, Trees and Proximity              */
/*                Representations, 1988, p. 23ff.                           */
/*                                                                          */  
/*  r = root, lt = line type, lw = line width, fs = font size               */
/*                                                                          */
/*  ptyp = 1:   radial drawing                                              */ 
/*         2:   hierarchical                                                */
/*         3:   axial                                                       */
/*                                                                          */
/*  The function assumes that a tree structure has been created with        */
/*  g_tree().                                                               */  
/*                                                                          */
/*  Return: 0 if OK, -1 if insufficient memory, -2 if unable to plot this   */
/*  tree (try another root).                                                */  

int pl_tree(TDAContext *ctx, int ptyp,int r,int lt,double lw,double fs) 
{
    register int i,j,k,mi;
    int err,n,m,np,s = 0,p,q,max,lc,rr,dr,ga,f;
    double xmin,xmax,ymin,ymax,a,ai,lg,gx,gy,dl,dm,lag,lad,ad,ag,zd,zg;
    double pi = 3.14159;

    n = ctx->GD_NP;
    m = n - 1;

    err = -1;
    if (alloc_aci(ctx, n + 1))
        goto PLTFin;

    if (alloc_acj(ctx, n + 1))
        goto PLTFin;

    if (alloc_ack(ctx, n + 1))
        goto PLTFin;

    if (alloc_acu(ctx, n + 1))
        goto PLTFin;

    if (alloc_acxf(ctx, n + 1))
        goto PLTFin;

    if (alloc_acyf(ctx, n + 1))
        goto PLTFin;

    if (alloc_acn(ctx, n + 1))
        goto PLTFin;

    err = -2;

    for (i = 1; i <= n; ++i)  
        ctx->AcN[i] = 0;

    for (i = 1; i <= m; ++i) {
        ctx->AcN[ctx->GT_EI[i]] += 1;
        ctx->AcN[ctx->GT_EJ[i]] += 1;
    }
    if (ptyp == 3)  
        goto PLTCont;

                        /* if radial or hierarchical drawing */
    mi = 0;
    k = r;
    ctx->AcJ[0] = r;
       
    while (k) {
        i = ctx->GT_SU[k];
        if (i > 0) {
            mi++;
            ctx->AcJ[mi] = k = i;
        }
        else {
            while (k) {
                i = ctx->GT_FR[k];
                if (i > 0) {
                    mi++;
                    ctx->AcJ[mi] = k = i;
                    break;
                }
                else  
                    k = ctx->GT_PR[k];
            }
        }
    }

PLTCont:
    switch (ptyp) {

        case 1:         /* radial drawing */
       
            np = 0;
            for (i = 1; i <= n; ++i) {
                ctx->AcU[i] = 0.0;
                if (ctx->AcN[i] == 1)
                    np++;
            }
            ai = 2.0 * pi / (double)np;
           
            k = 0;
            for (i = m; i >= 1; --i) {
                s = ctx->AcJ[i];
                if (ctx->AcN[s] > 1)  
                    ctx->AcU[s] /=  (double)(ctx->AcN[s] - 1);
                else {
                    ctx->AcU[s] = ai * (double)k;
                    k++;
                }
                p = ctx->GT_PR[s];
                ctx->AcU[p] += ctx->AcU[s];
            }
            ctx->AcXF[r] = ctx->AcYF[r] = 0.0;
            for (i = 1; i <= m; ++i) {
                s = ctx->AcJ[i];
                p = ctx->GT_PR[s];
                ai = ctx->AcU[s];
                lg = (double)(ctx->GT_EV[s]);
                ctx->AcXF[s] = (float)((double)(ctx->AcXF[p]) + lg * cos(ai));
                ctx->AcYF[s] = (float)((double)(ctx->AcYF[p]) - lg * sin(ai));
            }
            break;
 
        case 2:                 /* hierarchical drawing */

            ctx->AcXF[r] = ctx->AcYF[r] = 0.0;
            for (i = 1; i <= m; ++i) {
                s = ctx->AcJ[i];
                p = ctx->GT_PR[s];
                ctx->AcYF[s] = ctx->AcYF[p] + ctx->GT_EV[s];
                ctx->AcXF[s] = 0.0;
            }
            k = -1;
            for (i = m; i >= 1; --i) {
                s = ctx->AcJ[i];
                if (ctx->AcN[s] <= 1) {
                    k++;
                    ctx->AcXF[s] = (float)(50.0 * (double)k);
                }
                else
                    ctx->AcXF[s] = (float)((double)(ctx->AcXF[s]) / ((double)(ctx->AcN[s] - 1)));
                p = ctx->GT_PR[s];
                ctx->AcXF[p] += ctx->AcXF[s];
            }
            if (ctx->AcN[r] > 1)
                ctx->AcXF[r] = (float)((double)ctx->AcXF[r] / (double)ctx->AcN[r]);
            break;
           
        case 3:                         /* axial drawing */
            ctx->AcI[0] = r;
            mi = 0;
            i = -1;

            while (1) {
                i++;
                k = ctx->AcI[i];
                if (i > mi)   
                    break;
                p = ctx->GT_SU[k];
                if (p != 0) {
                    mi++;
                    ctx->AcI[mi] = p;

                    while (1) {
                        p = ctx->GT_FR[p];
                        if (p == 0)
                            break;
                        mi++;
                        ctx->AcI[mi] = p;
                    }
                }
            }
            ctx->AcK[r] = 0;
            for (i = 1; i <= mi; ++i) {
                j = ctx->AcI[i];
                ctx->AcK[j] = ctx->AcK[ctx->GT_PR[j]] + 1;
            }
            max = 0;
            for (i = 1; i <= n; ++i) {
                if (ctx->AcK[i] > max) {
                    s = i;
                    max = ctx->AcK[i];
                }
            }
            lc = 1;
            ctx->AcJ[lc] = s;

            while (1) {
                s = ctx->GT_PR[s];
                lc++;
                ctx->AcJ[lc] = s;
                if (s == r)
                    break;       
            }
            rr = lc;
            ctx->AcI[1] = ctx->AcJ[lc - 1];
            mi = 1;
    
            i = 0;
            while (1) {

                i++;
                k = ctx->AcI[i];
                if (i > mi)  
                    break;

                p = ctx->GT_SU[k];
                if (p != 0) {
                    mi++;
                    ctx->AcI[mi] = p;
                    while (1) {
                        p = ctx->GT_FR[p];
                        if (p == 0)
                            break;       
                        mi++;
                        ctx->AcI[mi] = p;
                    }
                }
            }
            for (i = 1; i <= mi; ++i) {
                ctx->AcK[ctx->AcI[i]] = 0;
            }
            max = 0;
            for (i = 1; i <= n; ++i) {
                if (ctx->AcK[i] > max) {
                    s = i;
                    max = ctx->AcK[i];
                }
            }
            if (max == 0)  
                goto PLTFin;
               
            lc += max;
            k = lc;
            ctx->AcJ[k] = s;

            while (1) {
                s = ctx->GT_PR[s];
                if (s == r)
                    break;         
                k--;
                ctx->AcJ[k] = s;
            }

            /* coordinates of diameter and vertices adjacent to endpoints */

            s = ctx->AcJ[2];
            ctx->AcXF[s] = ctx->AcYF[s] = 0.0;
            dr = ctx->AcJ[3];
            ai = pi / (double)ctx->AcN[s];
            p = ctx->GT_SU[s];
            k = 0;
            if (p == dr)
                goto L1250;
L1230:
            k++;
            a = pi / 2.0 + (double)k * ai;
            ctx->AcXF[p] = (float)((double)(ctx->AcXF[s]) + (double)(ctx->GT_EV[p]) * cos(a));
            ctx->AcYF[p] = (float)((double)(ctx->AcYF[s]) + (double)(ctx->GT_EV[p]) * sin(a));

L1250:
            p = ctx->GT_FR[p];
            if (p == dr)
                goto L1250;
            else if (p > 0)
                goto L1230;
           
            for (i = 3; i < lc; ++i) {
                s = ctx->AcJ[i];
                dr = ctx->AcJ[i + 1];
                ga = ctx->AcJ[i - 1];
                ctx->AcYF[s] = 0.0;
                if (i <= rr) {                          
                    ctx->AcXF[s] = ctx->AcXF[ga] + ctx->GT_EV[ga];
                }
                else
                    ctx->AcXF[s] = ctx->AcXF[ga] + ctx->GT_EV[s];
            }
            s = ctx->AcJ[lc - 1];
            ga = ctx->AcJ[lc - 2];
            ai = pi / (double)ctx->AcN[s];
            p = ctx->GT_SU[s];
            k = 0;
            if (p == ga)
                goto L1330;
L1310:
            k++;
            a = -pi / 2.0 + (double)k * ai;
            ctx->AcXF[p] = (float)((double)(ctx->AcXF[s]) + (double)(ctx->GT_EV[p]) * cos(a));
            ctx->AcYF[p] = (float)((double)(ctx->AcYF[s]) + (double)(ctx->GT_EV[p]) * sin(a));
L1330:
            p = ctx->GT_FR[p];
            if (p == ga)
                goto L1330;
            else if (p > 0)
                goto L1310;

            f = 1;                          /* location of other vertices */
            for (j = 3; j < lc - 1; ++j) {
                s = ctx->AcJ[j];
                dr = ctx->AcJ[j + 1];
                ga = ctx->AcJ[j - 1];
L1430:
                mi = 0;
                p = ctx->GT_SU[s];
                if (p == 0)
                    goto L1780;

                while (1) {
                    if (p != dr && p != ga) {
                        mi++;
                        ctx->AcI[mi] = p;
                    }
                    p = ctx->GT_FR[p];
                    if (p <= 0)
                        break;     
                }
                if (mi == 0)
                    goto L1780;
                else if (mi == 1) {
                    p = ctx->AcI[1];          /* only one successor */
                    ctx->AcXF[p] = ctx->AcXF[s];
                    ctx->AcYF[p] = (float)((double)ctx->AcYF[s] +
                               (double)f * (double)ctx->GT_EV[p]);
                    s = p;
                    goto L1430;
                }
                i = 0;
                while (1) {
                    i++;
                    k = ctx->AcI[i];
                    if (i > mi)  
                        break;

                    p = ctx->GT_SU[k];
                    if (p != 0) {
                        mi++;
                        ctx->AcI[mi] = p;
                        while (1) {
                            p = ctx->GT_FR[p];
                            if (p == 0)
                                break;       
                            mi++;
                            ctx->AcI[mi] = p;
                        }
                    }
                }
                np = 0;
                for (i = 1; i <= mi; ++i) {
                    if (ctx->AcN[ctx->AcI[i]] == 1)
                        np++;
                }
                dm = 0;
                for (i = 1; i <= mi; ++i) {
                    p = ctx->AcI[i];
                    ctx->AcU[p] = 0.0;
                    if (ctx->AcN[p] <= 1) {
                        dl = (double)(ctx->GT_EV[p]);

                        while (1) {
                            p = ctx->GT_PR[p];
                            if (p == s)  
                                break;
                            dl += (double)(ctx->GT_EV[p]);
                        }
                        if (dl > dm)
                            dm = dl;
                    }
                }
                lag = (double)(ctx->AcXF[s] - ctx->AcXF[ga]);
                lad = (double)(ctx->AcXF[dr] - ctx->AcXF[s]);

                /* calculation of angles at extreme left and right */

                if (dm <= lag)  
                    ag = pi / 4.0;
                else {
                    zg = sqrt(dm * dm - lag * lag);
                    ag = 0.0;
                    if (lag != 0.0)
                        ag = atan2(zg,lag);
                }
                if (dm <= lad)  
                    ad = pi / 4.0;
                else {
                    zd = sqrt(dm * dm - lad * lad);
                    ad = 0.0;
                    if (lad != 0.0)
                        ad = atan2(zd,lad);
                }
                if (f == 1)
                    ag = pi - ag;
                else {
                    ag = -pi + ag; 
                    ad = -ad;
                }
                ai = (ag - ad) / (np - 1);

                            /* coordinates for vertices in subtrees */
                k = 0;
                for (i = mi; i >= 1; --i) {
                    p = ctx->AcI[i];
                    if (ctx->AcN[p] > 1)  
                        ctx->AcU[p] /= (ctx->AcN[p] - 1);
                    else {
                        ctx->AcU[p] = ad + k * ai;
                        k++;
                    }
                    q = ctx->GT_PR[p];
                    if (q != s) 
                        ctx->AcU[q] += ctx->AcU[p];
                }
                for (i = 1; i <= mi; ++i) {
                    s = ctx->AcI[i];
                    p = ctx->GT_PR[s];
                    ctx->AcXF[s] = (float)((double)(ctx->AcXF[p]) + (double)(ctx->GT_EV[s]) * cos(ctx->AcU[s]));
                    ctx->AcYF[s] = (float)((double)(ctx->AcYF[p]) + (double)(ctx->GT_EV[s]) * sin(ctx->AcU[s]));
                }
L1780:
                f = -f;
            }
            break;
    }

    /* scale the plot for given coordinate system */

    xmax = (double)ctx->AcXF[1];


    xmin = (double)xmax;
    ymax = (double)ctx->AcYF[1];

    ymin = (double)ymax;
    for (i = 2; i <= n; ++i) {
        if ((double)(xmin) > (double)(ctx->AcXF[i])) 
            xmin = (double)(ctx->AcXF[i]);
        if ((double)(xmax) < (double)(ctx->AcXF[i])) 
            xmax = (double)(ctx->AcXF[i]);
        if ((double)(ymin) > (double)(ctx->AcYF[i])) 
            ymin = (double)(ctx->AcYF[i]);
        if ((double)(ymax) < (double)(ctx->AcYF[i])) 
            ymax = (double)(ctx->AcYF[i]);
    }
    zd = zg = 0.5;
    gx = gy = 0.0;
    if (xmax > xmin) {
        gx = 0.86 * ctx->UXLen / (xmax - xmin);
        zd = 0.07;
    } 
    if (ymax > ymin) {
        gy = 0.86 * ctx->UYLen / (ymax - ymin);
        zg = 0.07;
    }
    if (ptyp != 2 && gx > 0.0 && gy > 0.0) {
        a = gx * ctx->PXLen * ctx->UYLen / (ctx->PYLen * ctx->UXLen);
        if (gy >= a)
            gy = a;
        else
            gx *= gy / a;
    }
    for (i = 1; i <= n; ++i) {
        ctx->AcXF[i] = (float)(((double)(ctx->AcXF[i]) - xmin) * gx + zd * ctx->UXLen + ctx->PA1[0]);
        ctx->AcYF[i] = (float)(((double)(ctx->AcYF[i]) - ymin) * gy + zg * ctx->UYLen + ctx->PA1[1]);
    }
    if (ptyp == 2) {        /* change direction for hierarchical drawing */
        for (i = 1; i <= n; ++i)  
            ctx->AcYF[i] = (float)(ctx->PA2[1] - ((double)(ctx->AcYF[i]) - ctx->PA1[1]));
    }   

    /* drawing */
    
    fprintf(ctx->PSFd,"\n%%#%d: tree-plot\n",++ctx->PSONUM);

    fprintf(ctx->PSFd,"gsave\n");

    ps_ltyp(ctx, lt);     /* set line type */
    ps_lwidth(ctx, lw);   /* and line width */

    for (i = 1; i <= n; ++i) {
        j = ctx->GT_PR[i];
        if (j) {
            if (ptyp == 2) {
                ps_2dplot(ctx, (double)ctx->AcXF[i],(double)ctx->AcYF[j],0);    
                ps_2dplot(ctx, (double)ctx->AcXF[j],(double)ctx->AcYF[j],1);    
                fprintf(ctx->PSFd,"stroke\n");
                ps_2dplot(ctx, (double)ctx->AcXF[i],(double)ctx->AcYF[i],0);    
                ps_2dplot(ctx, (double)ctx->AcXF[i],(double)ctx->AcYF[j],1);    
            }
            else {
                ps_2dplot(ctx, (double)ctx->AcXF[i],(double)ctx->AcYF[i],0);    
                ps_2dplot(ctx, (double)ctx->AcXF[j],(double)ctx->AcYF[j],1);    
            }
            fprintf(ctx->PSFd,"stroke\n");
        }
    }
    if (fs > 0.0) {                 /* print labels */

        /* ## if PMNI != 0 only for leaves */
       
        for (i = 1; i <= ctx->PMNC; ++i) {
            if (i > n)
                break;

            if (ctx->PMNI == 0 || ctx->AcN[i] == 1)  
                ps_nlab(ctx, ctx->GD_ND[i - 1],(double)ctx->AcXF[i],(double)ctx->AcYF[i],fs,1);
        }
    }
    fprintf(ctx->PSFd,"grestore\n");
    err = 0;

PLTFin:
    alloc_acxf(ctx, 0);
    alloc_acyf(ctx, 0);
    alloc_acu(ctx, 0);
    alloc_aci(ctx, 0);
    alloc_acj(ctx, 0);
    alloc_ack(ctx, 0);
    alloc_acn(ctx, 0);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_mtree   Minimum spanning tree.                                        */
/*                                                                          */
/*      The code is adapted from AS13, Applied Statistics 18 (1969), 103f   */
/*                                                                          */
/*      int g_mtree(idx,dist,opt)                                           */
/*                                                                          */
/*      The function assumes a graph with GD_NP data points and uses        */
/*      gdd_adj(i,j,gn) to get the edge lengths for an undirected graph.    */
/*                                                                          */
/*      If opt == 0 the algorithm interprets zero edge values as valid,     */
/*      otherwise the algorithm assumes only edges with positive values.    */
/*                                                                          */
/*      Output: idx[i] = number of point connected with point i (i >= 1)    */ 
/*              dist[i] = distance between i and idx[i].                    */ 
/*                                                                          */
/*      Return: number of entries in idx[] if successful,                   */
/*      otherwise -1 if insufficient memory.                                */

int g_mtree(TDAContext *ctx, int gn,int *idx,double *dist,int opt)
{
    register int i,j,k,next;
    int n,nc;
    double big,min,d;

    n = ctx->GD_NP;
    big = ctx->DBLMAX;           
 
    if (alloc_acc(ctx, n))
        return(-1);

    for (i = 0; i < n; ++i) {
        idx[i] = -1;
        dist[i] = big;
    }
    next = j = 0;
    for (i = 1; i < n; ++i) {
        min = big;
        for (k = 1; k < n; ++k) {
            if (ctx->AcC[k] == 0) {
                d = gdd_adj(ctx, j,k,gn);   /* assume undirected graph */

                if (opt && d <= 0.0)
                    d = big;

                if (d < dist[k]) {
                    dist[k] = d;
                    idx[k] = j;
                }
                if (min > dist[k]) {
                    min = dist[k];
                    next = k;
                }
            }
        }
        j = next;
        ctx->AcC[j] = 1;
    }
    alloc_acc(ctx, 0);
    nc = 0;
    for (i = 1; i < n; ++i) {
        if (idx[i] >= 0)
            nc++;
    }
    return(nc);
}


