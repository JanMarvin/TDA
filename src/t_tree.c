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

/*  functions in t_tree.c */

int gt_alloc(int n);
int g_tree(int r);
int get_root(void);
int ptree(void);
int pltree(void);
int pl_tree(int ptyp,int r,int lt,double lw,double fs);

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */

int GT_RT = 0;      /* tree structure: root                                 */
int *GT_EI;         /* tree structure: starting point of edge               */
int *GT_EJ;         /* tree structure: end point of edge                    */
float *GT_EV;       /* tree structure: value of edge                        */
int *GT_PR;         /* tree structure: predecessor                          */
int *GT_SU;         /* tree structure: successor                            */
int *GT_FR;         /* tree structure: pointer to successor                 */

int GT_EIA = 0;     /* set if GT_EI allocated                               */
int GT_EJA = 0;     /* set if GT_EJ allocated                               */
int GT_EVA = 0;     /* set if GT_EV allocated                               */
int GT_PRA = 0;     /* set if GT_PR allocated                               */
int GT_SUA = 0;     /* set if GT_SU allocated                               */
int GT_FRA = 0;     /* set if GT_FR allocated                               */

/*  additional tree structure */

int TR_NE = 0;      /* number of edges                                      */
int *TR_PR;         /* predecessor                                          */
float *TR_LG;       /* length of edges                                      */
int TR_PRA = 0;     /* set if TR_PR allocated                               */
int TR_LGA = 0;     /* set if TR_LG allocated                               */

/* ------------------------------------------------------------------------ */
/*  gt_alloc(n)     If n > 0 allocate memory for tree data structure,       */
/*                  with n nodes. If n = 0 free previously allocated        */
/*                  memory.                                                 */
/*                                                                          */
/*  Return: 0 if OK, -1 if insufficient memory                              */

int gt_alloc(int n)
{
    int err = 0;

    if (n == 0)  
        goto GTALFin;
       
    err = -1;

    if (!(GT_EI = (int *)calloc(n + 1,sizeof(int))))  
        goto GTALFin; 
    memrq(n + 1,sizeof(int));
    GT_EIA = n + 1;

    if (!(GT_EJ = (int *)calloc(n + 1,sizeof(int))))  
        goto GTALFin; 
    memrq(n + 1,sizeof(int));
    GT_EJA = n + 1;

    if (!(GT_EV = (float *)calloc(n + 1,sizeof(float))))  
        goto GTALFin; 
    memrq(n + 1,sizeof(float));
    GT_EVA = n + 1;

    if (!(GT_PR = (int *)calloc(n + 1,sizeof(int))))  
        goto GTALFin; 
    memrq(n + 1,sizeof(int));
    GT_PRA = n + 1;

    if (!(GT_SU = (int *)calloc(n + 1,sizeof(int))))  
        goto GTALFin; 
    memrq(n + 1,sizeof(int));
    GT_SUA = n + 1;

    if (!(GT_FR = (int *)calloc(n + 1,sizeof(int))))  
        goto GTALFin; 
    memrq(n + 1,sizeof(int));
    GT_FRA = n + 1;

    return(0);

GTALFin:
    if (GT_EIA > 0) {
        free((char *)GT_EI);
        memrq(-GT_EIA,sizeof(int));
        GT_EIA = 0;
    }
    if (GT_EJA > 0) {
        free((char *)GT_EJ);
        memrq(-GT_EJA,sizeof(int));
        GT_EJA = 0;
    }
    if (GT_EVA > 0) {
        free((char *)GT_EV);
        memrq(-GT_EVA,sizeof(float));
        GT_EVA = 0;
    }
    if (GT_PRA > 0) {
        free((char *)GT_PR);
        memrq(-GT_PRA,sizeof(int));
        GT_PRA = 0;
    }
    if (GT_SUA > 0) {
        free((char *)GT_SU);
        memrq(-GT_SUA,sizeof(int));
        GT_SUA = 0;
    }
    if (GT_FRA > 0) {
        free((char *)GT_FR);
        memrq(-GT_FRA,sizeof(int));
        GT_FRA = 0;
    }
    GT_RT = 0;
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

int g_tree(int r)
{
    register int i,j,u,v,k;
    int err,m,n,ne,sdeg;                   
    float tmp;
               
    n = GD_NP;
            
    if (GD_NP != GD_NE[PMGN - 1] + 1) {
        err = -2;
        goto GTREEFin;
    }
              
    if (gt_alloc(n)) {
        p_err(-2,1);
        return(-1);
    }
    err = -1;
    if (alloc_acn(n + 1))
        goto GTREEFin;

    sdeg = m = 0;
    for (i = 2; i <= n; ++i) {              /* count edges and degrees */
        for (j = 1; j < i; ++j) {           /* and save edges */

            tmp = (float)gdd_adj(i - 1,j - 1,PMGN);         
            if (tmp >= 0.0) {
                m++;
                if (m < n) {
                    GT_EI[m] = i;
                    GT_EJ[m] = j;
                    GT_EV[m] = tmp;
                }
                AcN[i] += 1;
                AcN[j] += 1;
                sdeg += 2;
            }
        }
    }
    if (m + 1 != n) {                   /* if no tree */
        err = -2;
        goto GTREEFin;
    }
    if (alloc_aci(n + 1))
        goto GTREEFin;

    if (alloc_acj(n + 1))
        goto GTREEFin;

    if (alloc_acxf(n + 1))
        goto GTREEFin;

    /* find a representation for any root */

    while (sdeg > 0) {

        for (i = 1; i <= n; ++i)  
            AcI[i] = 0;

        ne = 0;
        for (i = 1; i <=  n; ++i) {

            if (AcJ[i])   
                continue;

            if (AcN[GT_EI[i]] == 1) {
                u = GT_EI[i];
                v = GT_EJ[i];
            }
            else if (AcN[GT_EJ[i]] == 1) {
                u = GT_EJ[i];
                v = GT_EI[i];
            }
            else
                continue;

            GT_PR[u] = v;
            AcJ[i] = 1;
            ne++;
            AcI[u] += 1;
            AcI[v] += 1;
            AcXF[u] = GT_EV[i];
            k = GT_SU[v];
            if (k == 0) 
                GT_SU[v] = u;
            else {
                while (GT_FR[k] > 0)  
                    k = GT_FR[k];
                GT_FR[k] = u;
            }
        }
        if (ne == 0) {
            err = -2;
            goto GTREEFin;
        }
        for (i = 1; i <= n; ++i) {
            AcN[i] -= AcI[i];
            sdeg -= AcI[i];
        }
    }
    for (i = 1; i <= n; ++i)
        GT_EV[i] = AcXF[i];

    /* current root is v. if v != r find new representation for root r */

    GT_RT = v;
    if (r == 0 || v == r) {
        err = 0;
        goto GTREEFin;
    }

    v = GT_RT;  /* current root */
    u = r;      /* required root */                 
    i = 0;
    GT_EV[v] = GT_EV[u];

    while (u > 0) {
        v = GT_PR[u];
        GT_PR[u] = i;
        i = u;
        if (v == 0)
            break;    

        k = GT_SU[u];
        if (k == 0)  
            GT_SU[u] = v;
        else {
            j = GT_FR[k];
            while (j > 0) {
                k = j;
                j = GT_FR[k];
            }
            GT_FR[k] = v;
        }
        k = GT_SU[v];
        if (k == u) {
            GT_SU[v] = GT_FR[u];
            GT_FR[u] = 0;
        }
        else {
            j = GT_FR[k];
            while (j > 0 && j != u) {
                k = j;
                j = GT_FR[k];
            }
            GT_FR[k] = GT_FR[u];
            GT_FR[u] = 0;
        }
        u = v;
        if (u <= 0)
            break;

        tmp = GT_EV[r];
        GT_EV[r] = GT_EV[v];
        GT_EV[v] = tmp;
    }
    GT_RT = r;   /* new root */
    GT_EV[r] = 0.0;
    err = 0;

GTREEFin:
    alloc_acn(0);
    alloc_aci(0);
    alloc_acj(0);
    alloc_acxf(0);
    if (err) {
        if (err == -2)
            printf1("Error: current graph is not a tree.\n");
        gt_alloc(0);
    }
    else
        printf1("Selected root: %d\n",GD_ND[GT_RT - 1]);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_root()  return internal root number based on PMRT                   */   
/*              if PMRT  not specified return 0.                            */
/*              if error, message and return -1.                            */

int get_root(void)
{
    register int i;
    int rt;

    if (PMRT < 1)  
        rt = 0;
    else {
        rt = -1;
        for (i = 0; i < GD_NP; ++i) {
            if (GD_ND[i] == PMRT) {
                rt = i + 1;
                break;
            }
        }
    }
    if (rt < 0) {
        printf1("Error: root node %d not available.\n",PMRT);
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

int ptree(void)
{
    register int i;
    int err,rt;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Print tree data structure. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,1,1))       /* get parameters */
        goto PTREEFin;

    if (PMFmtF == 0)            /* default print format */
        pmfmt(10,4);
/**
    if (gd_setup(PMGN,PMGT,PMGTT,0,0,0,1,2))     always undirected       
        goto PTREEFin;
**/
    if ((rt = get_root()) < 0)
        goto PTREEFin;
       
    if (g_tree(rt))             /* setup tree data structure */
        goto PTREEFin;

    fprintf(PMFd,"# tree data structure with root %d\n",GD_ND[GT_RT - 1]);
   
    for (i = 1; i <= GD_NP; ++i) {
        fprintf(PMFd,PMNFmtS,i);
        fprintf(PMFd,PMNFmtS,GD_ND[i - 1]);
        fprintf(PMFd,PMNFmtS,GT_EI[i]);
        fprintf(PMFd,PMNFmtS,GT_EJ[i]);
        fprintf(PMFd,PMFmtS,GT_EV[i]);
        fprintf(PMFd,PMNFmtS,GT_PR[i]);
        fprintf(PMFd,PMNFmtS,GT_SU[i]);
        fprintf(PMFd,PMNFmtS,GT_FR[i]);
        fprintf(PMFd,"\n");
    }    
    printf1("%d records written to: %s\n",GD_NP,PMFdName);
    err = 0;

PTREEFin:
    gt_alloc(0);
    p_clean();
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

int pltree(void)
{
    int err,rt;

    err = -1;
    if (check_pcmd(1,2))  /* check plot command for PostScript set up */ 
        return(-1);

    printf1("Plot tree. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 6 ,1,0))       /* get parameters */
        goto PLTREEFin;

    if (gdd_check(PMGN,1))
        goto PLTREEFin;
    if (gdd_tcheck(0,1,2))          /* need undirected valued graph */
        goto PLTREEFin;

    if ((rt = get_root()) < 0)
        goto PLTREEFin;
     
    if (g_tree(rt))             /* setup tree data structure */
        goto PLTREEFin;

    if (PMPLFlg < 1 || PMPLFlg > 3)
        PMPLFlg = 3;

    printf1("Type of drawing: %d ",PMPLFlg);

    if (PMPLFlg == 1) 
        printf1("(radial)\n");
    else if (PMPLFlg == 2) 
        printf1("(hierarchical)\n");
    else               
        printf1("(axial)\n");

    err = pl_tree(PMPLFlg,GT_RT,PMLT,PMLW,PMFS);
    if (err) {
        if (err == -2) {
            printf1("Unable to plot this tree, try another root.\n");
            err = 0;
        }
    }
    else  
        printf1("PostScript output written to: %s\n",PSFName);

PLTREEFin:
    gt_alloc(0);
    p_clean();
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

int pl_tree(int ptyp,int r,int lt,double lw,double fs) 
{
    register int i,j,k,mi;
    int err,n,m,np,s,p,q,max,lc,rr,dr,ga,f;
    float xmin,xmax,ymin,ymax,a,ai,lg,gx,gy,dl,dm,lag,lad,ad,ag,zd,zg;
    double pi = 3.14159;

    n = GD_NP;
    m = n - 1;

    err = -1;
    if (alloc_aci(n + 1))
        goto PLTFin;

    if (alloc_acj(n + 1))
        goto PLTFin;

    if (alloc_ack(n + 1))
        goto PLTFin;

    if (alloc_acu(n + 1))
        goto PLTFin;

    if (alloc_acxf(n + 1))
        goto PLTFin;

    if (alloc_acyf(n + 1))
        goto PLTFin;

    if (alloc_acn(n + 1))
        goto PLTFin;

    err = -2;

    for (i = 1; i <= n; ++i)  
        AcN[i] = 0;

    for (i = 1; i <= m; ++i) {
        AcN[GT_EI[i]] += 1;
        AcN[GT_EJ[i]] += 1;
    }
    if (ptyp == 3)  
        goto PLTCont;

                        /* if radial or hierarchical drawing */
    mi = 0;
    k = r;
    AcJ[0] = r;
       
    while (k) {
        i = GT_SU[k];
        if (i > 0) {
            mi++;
            AcJ[mi] = k = i;
        }
        else {
            while (k) {
                i = GT_FR[k];
                if (i > 0) {
                    mi++;
                    AcJ[mi] = k = i;
                    break;
                }
                else  
                    k = GT_PR[k];
            }
        }
    }

PLTCont:
    switch (ptyp) {

        case 1:         /* radial drawing */
       
            np = 0;
            for (i = 1; i <= n; ++i) {
                AcU[i] = 0.0;
                if (AcN[i] == 1)
                    np++;
            }
            ai = 2.0 * pi / (double)np;
           
            k = 0;
            for (i = m; i >= 1; --i) {
                s = AcJ[i];
                if (AcN[s] > 1)  
                    AcU[s] /=  (double)(AcN[s] - 1);
                else {
                    AcU[s] = ai * (double)k;
                    k++;
                }
                p = GT_PR[s];
                AcU[p] += AcU[s];
            }
            AcXF[r] = AcYF[r] = 0.0;
            for (i = 1; i <= m; ++i) {
                s = AcJ[i];
                p = GT_PR[s];
                ai = AcU[s];
                lg = GT_EV[s];
                AcXF[s] = AcXF[p] + lg * cos(ai);
                AcYF[s] = AcYF[p] - lg * sin(ai);
            }
            break;
 
        case 2:                 /* hierarchical drawing */

            AcXF[r] = AcYF[r] = 0.0;
            for (i = 1; i <= m; ++i) {
                s = AcJ[i];
                p = GT_PR[s];
                AcYF[s] = AcYF[p] + GT_EV[s];
                AcXF[s] = 0.0;
            }
            k = -1;
            for (i = m; i >= 1; --i) {
                s = AcJ[i];
                if (AcN[s] <= 1) {
                    k++;
                    AcXF[s] = 50.0 * (double)k;
                }
                else
                    AcXF[s] /= (double)(AcN[s] - 1);
                p = GT_PR[s];
                AcXF[p] += AcXF[s];
            }
            if (AcN[r] > 1)
                AcXF[r] /= AcN[r];
            break;
           
        case 3:                         /* axial drawing */
            AcI[0] = r;
            mi = 0;
            i = -1;

            while (1) {
                i++;
                k = AcI[i];
                if (i > mi)   
                    break;
                p = GT_SU[k];
                if (p != 0) {
                    mi++;
                    AcI[mi] = p;

                    while (1) {
                        p = GT_FR[p];
                        if (p == 0)
                            break;
                        mi++;
                        AcI[mi] = p;
                    }
                }
            }
            AcK[r] = 0;
            for (i = 1; i <= mi; ++i) {
                j = AcI[i];
                AcK[j] = AcK[GT_PR[j]] + 1;
            }
            max = 0;
            for (i = 1; i <= n; ++i) {
                if (AcK[i] > max) {
                    s = i;
                    max = AcK[i];
                }
            }
            lc = 1;
            AcJ[lc] = s;

            while (1) {
                s = GT_PR[s];
                lc++;
                AcJ[lc] = s;
                if (s == r)
                    break;       
            }
            rr = lc;
            AcI[1] = AcJ[lc - 1];
            mi = 1;
    
            i = 0;
            while (1) {

                i++;
                k = AcI[i];
                if (i > mi)  
                    break;

                p = GT_SU[k];
                if (p != 0) {
                    mi++;
                    AcI[mi] = p;
                    while (1) {
                        p = GT_FR[p];
                        if (p == 0)
                            break;       
                        mi++;
                        AcI[mi] = p;
                    }
                }
            }
            for (i = 1; i <= mi; ++i) {
                AcK[AcI[i]] = 0;
            }
            max = 0;
            for (i = 1; i <= n; ++i) {
                if (AcK[i] > max) {
                    s = i;
                    max = AcK[i];
                }
            }
            if (max == 0)  
                goto PLTFin;
               
            lc += max;
            k = lc;
            AcJ[k] = s;

            while (1) {
                s = GT_PR[s];
                if (s == r)
                    break;         
                k--;
                AcJ[k] = s;
            }

            /* coordinates of diameter and vertices adjacent to endpoints */

            s = AcJ[2];
            AcXF[s] = AcYF[s] = 0.0;
            dr = AcJ[3];
            ai = pi / (double)AcN[s];
            p = GT_SU[s];
            k = 0;
            if (p == dr)
                goto L1250;
L1230:
            k++;
            a = pi / 2.0 + (double)k * ai;
            AcXF[p] = AcXF[s] + GT_EV[p] * cos(a);
            AcYF[p] = AcYF[s] + GT_EV[p] * sin(a);

L1250:
            p = GT_FR[p];
            if (p == dr)
                goto L1250;
            else if (p > 0)
                goto L1230;
           
            for (i = 3; i < lc; ++i) {
                s = AcJ[i];
                dr = AcJ[i + 1];
                ga = AcJ[i - 1];
                AcYF[s] = 0.0;
                if (i <= rr) {                          
                    AcXF[s] = AcXF[ga] + GT_EV[ga];
                }
                else
                    AcXF[s] = AcXF[ga] + GT_EV[s];
            }
            s = AcJ[lc - 1];
            ga = AcJ[lc - 2];
            ai = pi / (double)AcN[s];
            p = GT_SU[s];
            k = 0;
            if (p == ga)
                goto L1330;
L1310:
            k++;
            a = -pi / 2.0 + (double)k * ai;
            AcXF[p] = AcXF[s] + GT_EV[p] * cos(a);
            AcYF[p] = AcYF[s] + GT_EV[p] * sin(a);
L1330:
            p = GT_FR[p];
            if (p == ga)
                goto L1330;
            else if (p > 0)
                goto L1310;

            f = 1;                          /* location of other vertices */
            for (j = 3; j < lc - 1; ++j) {
                s = AcJ[j];
                dr = AcJ[j + 1];
                ga = AcJ[j - 1];
L1430:
                mi = 0;
                p = GT_SU[s];
                if (p == 0)
                    goto L1780;

                while (1) {
                    if (p != dr && p != ga) {
                        mi++;
                        AcI[mi] = p;
                    }
                    p = GT_FR[p];
                    if (p <= 0)
                        break;     
                }
                if (mi == 0)
                    goto L1780;
                else if (mi == 1) {
                    p = AcI[1];          /* only one successor */
                    AcXF[p] = AcXF[s];
                    AcYF[p] = AcYF[s] + f * GT_EV[p];
                    s = p;
                    goto L1430;
                }
                i = 0;
                while (1) {
                    i++;
                    k = AcI[i];
                    if (i > mi)  
                        break;

                    p = GT_SU[k];
                    if (p != 0) {
                        mi++;
                        AcI[mi] = p;
                        while (1) {
                            p = GT_FR[p];
                            if (p == 0)
                                break;       
                            mi++;
                            AcI[mi] = p;
                        }
                    }
                }
                np = 0;
                for (i = 1; i <= mi; ++i) {
                    if (AcN[AcI[i]] == 1)
                        np++;
                }
                dm = 0;
                for (i = 1; i <= mi; ++i) {
                    p = AcI[i];
                    AcU[p] = 0.0;
                    if (AcN[p] <= 1) {
                        dl = GT_EV[p];

                        while (1) {
                            p = GT_PR[p];
                            if (p == s)  
                                break;
                            dl += GT_EV[p];
                        }
                        if (dl > dm)
                            dm = dl;
                    }
                }
                lag = AcXF[s] - AcXF[ga];
                lad = AcXF[dr] - AcXF[s];

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
                    p = AcI[i];
                    if (AcN[p] > 1)  
                        AcU[p] /= (AcN[p] - 1);
                    else {
                        AcU[p] = ad + k * ai;
                        k++;
                    }
                    q = GT_PR[p];
                    if (q != s) 
                        AcU[q] += AcU[p];
                }
                for (i = 1; i <= mi; ++i) {
                    s = AcI[i];
                    p = GT_PR[s];
                    AcXF[s] = AcXF[p] + GT_EV[s] * cos(AcU[s]);
                    AcYF[s] = AcYF[p] + GT_EV[s] * sin(AcU[s]);
                }
L1780:
                f = -f;
            }
            break;
    }

    /* scale the plot for given coordinate system */

    xmin = xmax = AcXF[1];
    ymin = ymax = AcYF[1];
    for (i = 2; i <= n; ++i) {
        if (xmin > AcXF[i]) 
            xmin = AcXF[i];
        if (xmax < AcXF[i]) 
            xmax = AcXF[i];
        if (ymin > AcYF[i]) 
            ymin = AcYF[i];
        if (ymax < AcYF[i]) 
            ymax = AcYF[i];
    }
    zd = zg = 0.5;
    gx = gy = 0.0;
    if (xmax > xmin) {
        gx = 0.86 * UXLen / (xmax - xmin);
        zd = 0.07;
    } 
    if (ymax > ymin) {
        gy = 0.86 * UYLen / (ymax - ymin);
        zg = 0.07;
    }
    if (ptyp != 2 && gx > 0.0 && gy > 0.0) {
        a = gx * PXLen * UYLen / (PYLen * UXLen);
        if (gy >= a)
            gy = a;
        else
            gx *= gy / a;
    }
    for (i = 1; i <= n; ++i) {
        AcXF[i] = (AcXF[i] - xmin) * gx + zd * UXLen + PA1[0];
        AcYF[i] = (AcYF[i] - ymin) * gy + zg * UYLen + PA1[1];
    }
    if (ptyp == 2) {        /* change direction for hierarchical drawing */
        for (i = 1; i <= n; ++i)  
            AcYF[i] = PA2[1] - (AcYF[i] - PA1[1]);
    }   

    /* drawing */
    
    fprintf(PSFd,"\n%%#%d: tree-plot\n",++PSONUM);

    fprintf(PSFd,"gsave\n");

    ps_ltyp(lt);     /* set line type */
    ps_lwidth(lw);   /* and line width */

    for (i = 1; i <= n; ++i) {
        j = GT_PR[i];
        if (j) {
            if (ptyp == 2) {
                ps_2dplot((double)AcXF[i],(double)AcYF[j],0);    
                ps_2dplot((double)AcXF[j],(double)AcYF[j],1);    
                fprintf(PSFd,"stroke\n");
                ps_2dplot((double)AcXF[i],(double)AcYF[i],0);    
                ps_2dplot((double)AcXF[i],(double)AcYF[j],1);    
            }
            else {
                ps_2dplot((double)AcXF[i],(double)AcYF[i],0);    
                ps_2dplot((double)AcXF[j],(double)AcYF[j],1);    
            }
            fprintf(PSFd,"stroke\n");
        }
    }
    if (fs > 0.0) {                 /* print labels */

        /* ## if PMNI != 0 only for leaves */
       
        for (i = 1; i <= PMNC; ++i) {
            if (i > n)
                break;

            if (PMNI == 0 || AcN[i] == 1)  
                ps_nlab(GD_ND[i - 1],(double)AcXF[i],(double)AcYF[i],fs,1);
        }
    }
    fprintf(PSFd,"grestore\n");
    err = 0;

PLTFin:
    alloc_acxf(0);
    alloc_acyf(0);
    alloc_acu(0);
    alloc_aci(0);
    alloc_acj(0);
    alloc_ack(0);
    alloc_acn(0);
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

int g_mtree(int gn,int *idx,double *dist,int opt)
{
    register int i,j,k,next;
    int n,nc;
    double big,min,d;

    n = GD_NP;
    big = DBLMAX;           
 
    if (alloc_acc(n))
        return(-1);

    for (i = 0; i < n; ++i) {
        idx[i] = -1;
        dist[i] = big;
    }
    next = j = 0;
    for (i = 1; i < n; ++i) {
        min = big;
        for (k = 1; k < n; ++k) {
            if (AcC[k] == 0) {
                d = gdd_adj(j,k,gn);   /* assume undirected graph */

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
        AcC[j] = 1;
    }
    alloc_acc(0);
    nc = 0;
    for (i = 1; i < n; ++i) {
        if (idx[i] >= 0)
            nc++;
    }
    return(nc);
}


