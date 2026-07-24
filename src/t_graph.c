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

/*  functions in t_graph.c */

int gni(void);
int gdln(void);
int gdp(void);
void p_gdp(int i,int j,int ii,int jj);
int gda(void);
int gdu(void);
int gde(void);
int gdcset(void);
int f_gdcset(int m,int *com,int *com1,int *nci,int *ncj);
int gsym(void);
int gdot(void);
int gtcl(void);
void g_tclos(int gn,int n,char *ai);
void g_spath(int gn,int n,float *a);
int gsort(void);
int gdcyc(void);
int g_dcyc(int gn,int n,int *p,int **h,int *hp,int opt,int ns,int cmax,
    int *minclen,int *maxclen,int dopt);
int gev(void); 
void g_op1(int n,double *z,double *w);
void g_evcheck(int r,int n,double *x,double *d);
int gsp(void);
int g_glist(int i,int *nodes);
void g_sp(int in,int n,int gn,double *d,int *ndec);
int tnet(void);


/* ------------------------------------------------------------------------ */
/*  gni()       Information about nodes.                                    */
/*                                                                          */
/*              gni(                                                        */
/*                  gn = ...,   graph number, def 1 or multigraph           */
/*                  nfmt = ..., integer format, def. 4                      */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gni(void)
{
    register int i,k;
    int err,n,nrec;

    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("Degree of nodes. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto GNIFin;

    if (gdd_check(PMGN,1))
        goto GNIFin;

    if (alloc_acr(GD_NG + 1))
        goto GNIFin;
    if (alloc_acs(GD_NG + 1))
        goto GNIFin;
    if (alloc_acn(GD_NG + 1))
        goto GNIFin;

    nrec = 0;
    for (i = 0; i < GD_NP; ++i) {
        fprintf(PMFd,PMNFmtS,i + 1);                
        fprintf(PMFd,PMNFmtS,gdd_node(i));                
  
        for (k = 1; k <= GD_NG; ++k) {
            if (PMGN != 0 && k != PMGN)
                continue;
   
            n = g_pre(i,k);   
            fprintf(PMFd,PMNFmtS,n);                
            if (AcR[k] < n)
                AcR[k] = n;

            if (GD_GT > 2)
                n = g_suc(i,k);   
            fprintf(PMFd,PMNFmtS,n);                
            if (AcS[k] < n)
                AcS[k] = n;

            n = g_loop(i,k);   
            fprintf(PMFd,PMNFmtS,n);                
            if (n > 0)
                AcN[k] += 1;
        }
        fprintf(PMFd,"\n");                
        nrec++;
    }
    printf1("\nGraph  max in-degree  max out-degree  nodes with loops\n");
    for (k = 1; k <= GD_NG; ++k)  
        printf1("%5d  %8d       %8d       %14d\n",k,AcR[k],AcS[k],AcN[k]);

    printf1("\n%d records written to: %s\n",nrec,PMFdName);
    err = 0;

GNIFin:       
    p_clean();
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

int gdln(void)
{
    register int i,j,k,l;
    int err,nrec,n,maxl;

    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("Forward/backward links. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto GFLFin;

    if (gdd_check(PMGN,0))
        goto GFLFin;
    if (gdd_tcheck(1,2,0))
        goto GFLFin;

    if (PMOPT == 1)
        printf1("Option 1: forward links.\n");
    else {
        printf1("Option 2: backward links.\n");
        PMOPT = 2;
    }
    if (alloc_acn(imax(GD_FPMax,GD_BPMax) + 2))
        goto GFLFin;

    maxl = nrec = 0;
    for (i = 0; i < GD_NP; ++i) {
        if (PMOPT == 1)
            n = GD_FPN[i];
        else
            n = GD_BPN[i];
        if (n <= 0)
            continue;

        l = 0;
        for (j = 0; j < n; ++j) {

            if (PMOPT == 1)
                k = GD_FPK[i][j];
            else
                k = GD_BPK[i][j];

            if (k >= 0 && gdd_ev(k,PMGN) >= 0.0) {
                if (PMOPT == 1)
                    AcN[l++] = gdd_node(GD_FPI[i][j]);        
                else
                    AcN[l++] = gdd_node(GD_BPI[i][j]);        
            }
        }
        if (l > 0) {
            fprintf(PMFd,PMNFmtS,i + 1);        
            fprintf(PMFd,PMNFmtS,gdd_node(i));        
            fprintf(PMFd,PMNFmtS,l);        
            for (j = 0; j < l; ++j)
                fprintf(PMFd,PMNFmtS,AcN[j]);        
            fprintf(PMFd,"\n");                
            nrec++;
            if (maxl < l)   
                maxl = l;
        }
    }
    printf1("Max number of links: %d\n",maxl);
    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

GFLFin:       
    p_clean();
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

int gdp(void)
{
    register int i,j,k;
    int err,nrec,m,ni,ii,jj;
    double tmp;

    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("Writing graph data. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto GDPFin;
   
    if (gdd_check(PMGN,1))
        goto GDPFin;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    if (PMOPT > 5)
        PMOPT = 1;

    if (PMSCFlg == 0)
        PMSC = -1.0;

    if (PMOPT >= 5 && PMGN == 0)  
        PMGN = 1;

    if (alloc_acn(GD_NP))                 
        goto GDPFin;

    ni = 0;
    if (PMF2Def) {
        if (alloc_ack(GD_NP))                 
            goto GDPFin;

        if ((ni = g_getnd(0,AcK)) < 0)         /* get node list */
            goto GDPFin;

        printf1("Number of input nodes: %d\n",ni);
        if (ni == 0)
            goto GDPFin;
    }
    if (alloc_actmp(GD_NG + 1))                 
        goto GDPFin;
       
    nrec = 0;

    if (PMOPT == 1) {
        for (i = 0; i < GD_NP; ++i) {
            n_info(i,100);
            if (ni > 0 && AcK[i] == 0)
                continue;
            ii = gdd_node(i);

            for (j = 0; j < GD_NP; ++j) {
                if (ni > 0 && AcK[j] == 0)
                    continue;

                jj = gdd_node(j);

                if (PMGN != 0) {
                    tmp = gdd_adj(i,j,PMGN);
                    if (tmp >= 0.0) {
                        p_gdp(i + 1,j + 1,ii,jj);
                        fprintf(PMFd,PMFmtS,tmp);
                        fprintf(PMFd,"\n");                
                        nrec++;    

                        AcN[i] = AcN[j] = 1;
                    }
                }
                else {
                    m = 0;
                    for (k = 1; k <= GD_NG; ++k) {
                        if (gdd_adj(i,j,k) >= 0.0) {
                            m = 1;
                            break;
                        }
                    }
                    if (m) {
                        p_gdp(i + 1,j + 1,ii,jj);
                        for (k = 1; k <= GD_NG; ++k) {
                            tmp = gdd_adj(i,j,k);
                            if (tmp < 0.0)
                                tmp = PMSC;
                            fprintf(PMFd,PMFmtS,tmp);
                        }
                        fprintf(PMFd,"\n");                
                        nrec++;    

                        AcN[i] = AcN[j] = 1;
                    }
                }
            }
        }
      
        /* add isolated nodes */

        tmp = PMSC;
        for (i = 0; i < GD_NP; ++i) {
            if (AcN[i] == 0) {
                if (ni > 0 && AcK[i] == 0)
                    continue;
                ii = gdd_node(i);

                if (PMGN != 0) {
                    p_gdp(i + 1,i + 1,ii,ii);
                    fprintf(PMFd,PMFmtS,tmp);
                }
                else {
                    p_gdp(i + 1,i + 1,ii,ii);
                    for (k = 1; k <= GD_NG; ++k)  
                        fprintf(PMFd,PMFmtS,tmp);
                }
                fprintf(PMFd,"\n");                
                nrec++;    
            }
        }
    }
    else if (PMOPT == 2 || PMOPT == 3 || PMOPT == 4) {
        for (i = 0; i < GD_NP; ++i) {
            n_info(i,100);
            if (ni > 0 && AcK[i] == 0)
                continue;

            ii = gdd_node(i);

            for (j = 0; j < GD_NP; ++j) {

                if ((PMOPT == 2 && j >= i) || (PMOPT == 3 && j > i))
                    continue;

                if (ni > 0 && AcK[j] == 0)
                    continue;

                jj = gdd_node(j);

                if (PMGN != 0) {
                    p_gdp(i + 1,j + 1,ii,jj);
                    tmp = gdd_adj(i,j,PMGN);
                    if (tmp < 0.0)
                        tmp = PMSC;
                    fprintf(PMFd,PMFmtS,tmp);
                    fprintf(PMFd,"\n");                
                    nrec++;    
                }
                else {
                    p_gdp(i + 1,j + 1,ii,jj);
                    for (k = 1; k <= GD_NG; ++k) {
                        tmp = gdd_adj(i,j,k);
                        if (tmp < 0.0)
                            tmp = PMSC;
                        fprintf(PMFd,PMFmtS,tmp);
                    }
                    fprintf(PMFd,"\n");                
                    nrec++;    
                }
            }
        }
    }
    else if (PMOPT == 5) {
        for (i = 0; i < GD_NP; ++i) {
            n_info(i,100);
            if (ni > 0 && AcK[i] == 0)
                continue;

            ii = gdd_node(i);

            fprintf(PMFd,PMNFmtS,i + 1);                
            fprintf(PMFd,PMNFmtS,ii);                

            for (j = 0; j < GD_NP; ++j) {

                if (ni > 0 && AcK[j] == 0)
                    continue;

                tmp = gdd_adj(i,j,PMGN);
                if (tmp < 0.0)
                    tmp = PMSC;
                fprintf(PMFd,PMFmtS,tmp);
            }
            fprintf(PMFd,"\n");                
            nrec++;    
        }
    }
    n_info_e();
    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

GDPFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  p_gdp(i,j,ii,jj)                                                        */

void p_gdp(int i,int j,int ii,int jj)
{
    fprintf(PMFd,PMNFmtS,i);                
    fprintf(PMFd,PMNFmtS,j);                
    fprintf(PMFd,PMNFmtS,ii);                
    fprintf(PMFd,PMNFmtS,jj);                
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

int gda(void)
{
    register int i,j,k;
    int err,nrec,n,m,mm,nm,ii,jj,fnd,con; 
    char *p;
    double tmp;

    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("Aggregating nodes of a graph. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto GDAFin;
   
    if (gdd_check(PMGN,1))
        goto GDAFin;

    if (GD_GT != 2 && GD_GT != 4) {
        printf1("Error: command requires valued graph.\n");
        goto GDAFin;
    }

    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    if (PRCNAlloc == 0) {
        printf1("No recode information. Nothing done.\n");
        err = 0;
        goto GDAFin;
    }
    printf1("Recode: %s\n",PRCN);
               
    if (alloc_acn(GD_NP))                 
        goto GDAFin;
    if (alloc_acm(GD_NP))                 
        goto GDAFin;

    for (i = 0; i < GD_NP; ++i)  
        AcM[i] = AcN[i] = gdd_node(i);

    err = -2;
    p = PRCN;
    while (*p) {
        if (sscanf(p,"%d",&n) != 1)
            goto GDAFin;
        if (n < 1) {
            err = -4;
            goto GDAFin;
        }
        p = skip_int(p);
        if (*p++ != '[')
            goto GDAFin;

        while (*p) {
            if (sscanf(p,"%d",&m) != 1)
                goto GDAFin;

            fnd = 0;
            for (i = 0; i < GD_NP; ++i) {
                if (m == AcN[i]) {
                    AcM[i] = n;
                    fnd = 1;
                    break;
                }
            }
            if (fnd == 0) {
                err = -3;
                goto GDAFin; 
            }
            p = skip_int(p);
            if (*p == ']')
                break;

            if (*p++ != ',')   
                goto GDAFin;

            if (*p != ',')
                continue;

            if (sscanf(++p,"%d",&mm) != 1)
                goto GDAFin;

            fnd = 0;
            for (i = 0; i < GD_NP; ++i) {
                if (mm == AcN[i]) {
                    AcM[i] = n;
                    fnd = 1;
                    break;
                }
            }
            if (fnd == 0) {
                err = -3;
                goto GDAFin; 
            }
            for (j = m + 1; j < mm; ++j) {
                for (i = 0; i < GD_NP; ++i) {
                    if (j == AcN[i]) {
                        AcM[i] = n;
                        break;
                    }
                }
            }
            p = skip_int(p);
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
    for (i = 0; i < GD_NP; ++i)
        AcN[i] = AcM[i];

    if (sorti(GD_NP,AcN,0))
        goto GDAFin;

    mm = 1;
    for (i = 1; i < GD_NP; ++i) {
        if (AcN[i] != AcN[i - 1])
            AcN[mm++] = AcN[i];
    }
    nm = AcN[mm - 1];

    if (alloc_ack(nm + 1))                 
        goto GDAFin;

    for (i = 0; i < mm; ++i)  
        AcK[AcN[i]] = i;

    if (alloc_acx(mm * mm + 1))                 
        goto GDAFin;

    for (i = 0; i <= mm * mm; ++i)
        AcX[i] = -1.0;

    for (i = 0; i < GD_NP; ++i) {
        n_info(i,100);
        for (j = 0; j < GD_NP; ++j) {
            if (GD_GT == 2 && j > i)
                break;

            tmp = gdd_adj(i,j,PMGN);
            if (tmp >= 0.0) {
                ii = AcK[AcM[i]];
                jj = AcK[AcM[j]];
                       
                /** printf("i=%4d j=%4d ii=%5d jj=%5d : %d %d : v=%g      \n",i,j,ii,jj,AcN[ii],AcN[jj],tmp    ); **/
                       
                k = ii * mm + jj;           
                if (AcX[k] < 0.0)
                    AcX[k] = tmp;
                else
                    AcX[k] += tmp;
            }
        }
    }
    n_info_e();

    nrec = 0;
    for (i = 0; i < mm; ++i) {
        for (j = 0; j < mm; ++j) {
            k = i * mm + j;
            if (AcX[k] >= 0.0) {
                p_gdp(i + 1,j + 1,AcN[i],AcN[j]);
                fprintf(PMFd,PMFmtS,AcX[k]);
                fprintf(PMFd,"\n");                
                nrec++;    
            }
        }
    }
    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

GDAFin:       
    if (err == -2) {
        printf1("Syntax error in recode string.\n");
        err = -1;
    }
    else if (err == -3) {
        printf1("Error: unknown node number in recode string.\n");
        err = -1;
    }
    else if (err == -4) {
        printf1("Error: node numbers must be positive.\n");
        err = -1;
    }
    p_clean();
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

int gdu(void)
{
    register int i,j,k,l;        
    int err,n,nn,iv,jv,i1,i2,j1,j2,kp,f,nrec,ne,ni;
    double tmp;

    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("Union of two graphs. Current memory: %d bytes.\n",MemReq);

    if (gdd_check(0,1))
        goto GDUFin;

    if (parm(CmdBuf + 3,4,1))       /* get parameters */
        goto GDUFin;
   
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    if (PMSCFlg == 0)
        PMSC = -1.0;

    if (PMNV < 3) {
        printf1("Error: need at least three variables on right-hand side.\n");
        goto GDUFin;
    }
    nn = GD_NP + 2 * NOC;           
    if (alloc_acn(nn + 1))
        goto GDUFin;
    if (alloc_aci(NOC + 1))
        goto GDUFin;
    if (alloc_acj(NOC + 1))
        goto GDUFin;
    if (alloc_ack(NOC + 1))
        goto GDUFin;

    iv = (int)PMVIdx[0];
    jv = (int)PMVIdx[1];

    /* sort I and J in data matrix */

    for (i = 0; i < NOC; ++i) {
        AcI[i] = (int)get_data(iv,i);
        AcJ[i] = (int)get_data(jv,i);
        if (AcI[i] < 1 || AcJ[i] < 1) {
            printf1("Error: node number in record %d is zero or negative.\n",i + 1);
            goto GDUFin;
        }
    }
    if (sortdpi2(NOC,AcI,AcJ,AcK))
        goto GDUFin;

    /* check for unique edges */

    i = AcI[AcK[0]];
    j = AcJ[AcK[0]];
    for (k = 1; k < NOC; ++k) {
        kp = AcK[k];
        if (AcI[kp] == i && AcJ[kp] == j) { 
            printf1("Error: multiple occurrence of edge %d - %d.\n",i,j);
            goto GDUFin;
        }
        i = AcI[kp]; 
        j = AcJ[kp];
    }

    /* get list of all node numbers */

    k = 0;
    for (i = 0; i < GD_NP; ++i)
        AcN[k++] = gdd_node(i);

    for (i = 0; i < NOC; ++i) {
        AcN[k++] = AcI[i];
        AcN[k++] = AcJ[i];
    }
    if (sorti(nn,AcN,0))
        goto GDUFin;

    n = 1;
    for (i = 1; i < nn; ++i) {
        if (AcN[i] != AcN[i - 1])
            AcN[n++] = AcN[i];
    }
    printf1("Union of graphs: %d nodes.\n",n);         

    /* create array to flag isolated nodes */

    nn = AcN[n - 1];     /* max node number */

    if (alloc_acns(nn + 1))
        goto GDUFin;

    for (i = 0; i < n; ++i)  
        AcNS[AcN[i]] = 1;

    /* write edges */

    ne = ni = nrec = 0;
    i = j = 0;

    for (k = 0; k <= NOC; ++k) {
        if (k < NOC) {
            kp = AcK[k];
            i2 = AcI[kp];
            j2 = AcJ[kp];
        }
        else  
            i2 = j2 = INTMAX;

        while (i < GD_NP && j < GD_NP) {
            i1 = gdd_node(i);
            j1 = gdd_node(j);
       
            f = 0;
            if (i1 < i2 || (i1 == i2 && j1 < j2)) {
                for (l = 1; l <= GD_NG; ++l) {
                    if (gdd_adj(i,j,l) >= 0.0) {
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
                if (PMF1Def) {
                    fprintf(PMF1d,PMNFmtS,i1);
                    fprintf(PMF1d,PMNFmtS,j1);
                   
                    for (l = 1; l <= GD_NG; ++l) {
                        tmp = gdd_adj(i,j,l);
                        if (tmp < 0.0)
                            tmp = PMSC;
                        fprintf(PMF1d,PMFmtS,tmp);
                    }
                    for (l = 2; l < PMNV; ++l) {
                        tmp = PMSC;
                        if (f == 2 && k < NOC) {
                            tmp = get_data((int)PMVIdx[l],kp);
                            if (tmp < 0.0)
                                tmp = PMSC;
                        }
                        fprintf(PMF1d,PMFmtS,tmp);
                    }
                    fprintf(PMF1d,"\n");                
                    nrec++;
                }
                ne++;
                AcNS[i1] = AcNS[j1] = 0;
            }
            if (++j >= GD_NP) {
                i++;
                j = 0;
            }
            if (f == 2)
                break;
        }
        if (k == NOC)
            break;

        if (i1 != i2 || j1 != j2) {

            f = 0;
            for (l = 2; l < PMNV; ++l) {
                if (get_data((int)PMVIdx[l],kp) >= 0.0) {
                    f = 1;
                    break;
                }
            }
            if (f) {
                if (PMF1Def) {
                    fprintf(PMF1d,PMNFmtS,i2);
                    fprintf(PMF1d,PMNFmtS,j2);

                    for (l = 1; l <= GD_NG; ++l)  
                        fprintf(PMF1d,PMFmtS,PMSC);
                     
                    for (l = 2; l < PMNV; ++l) {
                        tmp = get_data((int)PMVIdx[l],kp);
                        if (tmp < 0.0)
                            tmp = PMSC;
                        fprintf(PMF1d,PMFmtS,tmp);
                    }
                    fprintf(PMF1d,"\n");                
                    nrec++;
                }
                ne++;
                AcNS[i2] = AcNS[j2] = 0;
            }
        }
    }
    k = GD_NG + PMNV - 2;
    for (i = 1; i <= nn; ++i) {     /* check for isolated nodes */
        if (AcNS[i]) {
            if (PMF1Def) {
                fprintf(PMF1d,PMNFmtS,i);
                fprintf(PMF1d,PMNFmtS,i);
                for (l = 1; l <= k; ++l)  
                    fprintf(PMF1d,PMFmtS,PMSC);
                fprintf(PMF1d,"\n");                
                nrec++;
            }
            ni++;
        }
    }
    printf1("Number of valid edges: %d\n",ne);
    printf1("Number of isolated notes: %d\n",ni);
    if (PMF1Def)  
        printf1("%d records written to: %s\n",nrec,PMF1dName);
    err = 0;

GDUFin:       
    p_clean();
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

int gde(void)
{
    register int i,j,k,l;          
    int err,n,ne,ii,iv,ni,nj,jp,lp,nmax;     
    double x;

    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("Create an undirected graph. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 3,4,1))       /* get parameters */
        goto GDEFin;
   
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    if (PMNV != 2) {
        printf1("Error: need exactly two variables on right-hand side.\n");
        goto GDEFin;
    }
    if (PMF1Def == 0) {
        printf1("Error: need output file (df parameter).\n");
        goto GDEFin;
    }
    if (alloc_aci(NOC + 1))
        goto GDEFin;
    if (alloc_acj(NOC + 1))
        goto GDEFin;
    if (alloc_ack(NOC + 1))
        goto GDEFin;
    if (alloc_acx(NOC + 1))
        goto GDEFin;
    if (alloc_acy(NOC + 1))
        goto GDEFin;

    ii = (int)PMVIdx[0];
    iv = (int)PMVIdx[1];

    for (i = 0; i < NOC; ++i) {
        j = (int)get_data(ii,i);
        if (j < 1) {
            printf1("Error: node number in record %d is zero or negative.\n",i + 1);
            goto GDEFin;
        }
        AcI[i] = j;                     
    }
    if (sortdpi(NOC,AcI,AcK))
        goto GDEFin;

    j = AcI[AcK[0]];
    AcJ[0] = j;
    n = 1;
    nj = 0;
    for (i = 1; i < NOC; ++i) {
        k = AcK[i];
        if (AcI[k] != j) {
            j = AcI[k];
            n++;
            AcJ[++nj] = j;
        }
    }
    nmax = AcJ[nj];
    nj++;
    printf1("Number of nodes: %d\n",n);
   
    /**
    printf("AcJ nj=%d : ",nj);
    for (i = 0; i < nj; ++i)
        printf("%d ",AcJ[i]);
    newline();
    **/

    if (alloc_acr(nmax + 1))
        goto GDEFin;

    for (i = 0; i < nj; ++i)  
        AcR[AcJ[i]] = i;
    /**
    printf("AcR nm=%d : ",nmax);
    for (i = 0; i <= nmax; ++i)
        printf("%d ",AcR[i]);
    newline();
    **/

    if (alloc_acn(n * n + 1))
        goto GDEFin;

    for (i = 0; i < NOC; ++i) {
        AcX[i] = get_data(iv,i);
        AcY[i] = get_data(ii,i);
    }
    if (sortdp2(NOC,AcX,AcY,AcK))
        goto GDEFin;

    ne = 0;
    x = AcX[AcK[0]];
    AcI[0] = (int)AcY[AcK[0]];
    ni = 1;
    for (i = 1; i <= NOC; ++i) {
        if (i < NOC) 
            k = AcK[i];

        if (i == NOC || AcX[k] != x) {
            for (j = 1; j < ni; ++j) {
                jp = AcR[AcI[j]];
                for (l = 0; l < j; ++l) {
                    lp = AcR[AcI[l]];
                    AcN[jp * n + lp] += 1;
                    AcN[lp * n + jp] += 1;
                }
            }
            if (i == NOC)
                break;

            x = AcX[k];
            AcI[0] = (int)AcY[k];
            ni = 1;
        }
        else
            AcI[ni++] = (int)AcY[k];
    }
    /**
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j)
            printf("%4d ",AcN[i * n + j]);
        newline();
    }
    **/

    ne = ni = 0;
    if (PMPRNO == 0) {
        for (i = 0; i < n; ++i) {
            for (j = 0; j <= i; ++j) {
                if (j == i && PMNI == 1)
                    continue;

                k = AcN[i * n + j];
                if (k > 0) {
                    ne++;
                    if (PMF1Def) {
                        fprintf(PMF1d,PMNFmtS,AcJ[j]);
                        fprintf(PMF1d,PMNFmtS,AcJ[i]);
                        fprintf(PMF1d,PMNFmtS,k);
                        fprintf(PMF1d,"\n");
                        ni++;
                    }
                }
            }
        }
    }
    else {
        for (i = 0; i < n; ++i) {
            if (PMF1Def)   
                fprintf(PMF1d,PMNFmtS,AcJ[i]);
            for (j = 0; j < n; ++j) {
                k = AcN[i * n + j];
                if (k > 0) {
                    if (PMNI == 0 || i != j)
                        ne++;
                }
                if (PMF1Def)   
                    fprintf(PMF1d,PMNFmtS,k);
            }
            fprintf(PMF1d,"\n");
            ni++;
        }
        if (ne > 1)
            ne /= 2;
    }
    printf1("Number of edges: %d\n",ne);
    if (PMF1Def)  
        printf1("%d records written to: %s\n",ni,PMF1dName);
    err = 0;

GDEFin:       
    p_clean();
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

int gdcset(void)
{
    int err,m,r,first,nb;

    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("Find compact blocks. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 6,1,1))       /* get parameters */
        goto GDCSETFin;
   
    if (gdd_check(PMGN,1))
        goto GDCSETFin;

    if (GD_GT != 4) {
        printf1("Error: command requires directed and valued graph.\n");
        goto GDCSETFin;
    }
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    if (alloc_acn(GD_NP + 1))
        goto GDCSETFin;
    if (alloc_acm(GD_NP + 1))
        goto GDCSETFin;
    if (alloc_aci(GD_NP + 1))
        goto GDCSETFin;
    if (alloc_acj(GD_NP + 1))
        goto GDCSETFin;

    nb = 0;
    for (m = 2; m < GD_NP; ++m) {

        first = 1;
        while (1) {
            if (comb_nm(GD_NP,m,AcN,first) == 0) 
                break;
            r = f_gdcset(m,AcN,AcM,AcI,AcJ);
            if (r < 0)
                goto GDCSETFin;
            if (r > 0)
                nb++;
            first = 0;
        }
    }
    /**  printf("found nb=%d\n",nb); **/
    err = 0;

GDCSETFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  f_gdcset(m,com)                                                         */

int f_gdcset(int m,int *com,int *com1,int *nci,int *ncj)
{
    register int i,j,ii;
    int f,ne;
    double amin,aij;
/**
    printf("m=%3d : ",m);
    for (i = 0; i < m; ++i)
        printf("%2d ",com[i]);
    newline();
**/
    for (i = 0; i < m; ++i)
        com1[i] = com[i];

    ne = g_mst(PMGN,m,com1,nci,ncj,0);
    if (ne < m - 1)
        return(0);

    amin = DBLMAX;
    f = 0;
    for (i = 0; i < ne; ++i) {
        if ((aij = gdd_adj(nci[i],ncj[i],PMGN)) > 0.0) {
            amin = dmin(amin,aij);
            f = 1;
        }
        if ((aij = gdd_adj(ncj[i],nci[i],PMGN)) > 0.0) {
            amin = dmin(amin,aij);
            f = 1;
        }
    }
    if (f == 0) {
        printf("??????????\n");
        return(0);
    }
/*  printf("amin=%g\n",amin);
*/

    if (alloc_acr(GD_NP + 1))
        return(-1);       

    for (i = 0; i < m; ++i)
        AcR[com[i]] = 1;
     
    for (i = 0; i < m; ++i) {
        ii = com[i];
        for (j = 0; j < GD_NP; ++j) {
            if (AcR[j])
                continue;

            if ((aij = gdd_adj(ii,j,PMGN)) > 0.0) {
                if (aij >= amin)   
                    return(0);
            }
            if ((aij = gdd_adj(j,ii,PMGN)) > 0.0) {
                if (aij >= amin)  
                    return(0);
            }
        }
    }
    /**
    printf("found m=%3d : ",m);
    for (i = 0; i < m; ++i)
        printf("%2d ",com[i]);
    newline();
    **/

    return(1);
}

/* ------------------------------------------------------------------------ */
/*  gsym        Check whether a directed graph is symmetrical.              */
/*              gsym (gn=...);      gn = graph number, def. 1               */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int gsym(void)
{
    register int i,j;
    int err,n;

    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("Check graph for symmetry. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,1,0))       /* get parameters */
        goto GSYMFin;
   
    if (gdd_check(PMGN,1))
        goto GSYMFin;

    if (GD_GT != 3 && GD_GT != 4) {
        printf1("Error: command requires directed graph.\n");
        goto GSYMFin;
    }
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    n = 1;
    for (i = 1; i < GD_NP; ++i) {
        for (j = 0; j < i; ++j) {
            if (fabs(gdd_adj(i,j,PMGN) - gdd_adj(j,i,PMGN)) > 1.e-6) {
                n = 0;
                break;
            }   
        }
        if (n == 0)
            break;

    }
    printf1("Graph is ");
    if (n == 0)
        printf1("not ");
    printf1("symmetrical.\n");
    if (n == 0) {
        printf1("Nodes %5d %5d : %18.8f\n",i+1,j+1,gdd_adj(i,j,PMGN));
        printf1("Nodes %5d %5d : %18.8f\n",j+1,i+1,gdd_adj(j,i,PMGN));
    }
    err = 0;

GSYMFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdot        Print directed graph into output file according to the      */
/*              dot syntax.                                                 */
/*                                                                          */
/*              gdot = fname;       output file name, required.             */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int gdot(void)
{
    register int i,j;
    int err,nrec;
    double tmp;

    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("Print graph in dot format. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto GDOTFin;
   
    if (gdd_check(PMGN,1))
        goto GDOTFin;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    nrec = 0;

    if (GD_GT == 1 || GD_GT == 2) {     /* undirected graph */
        fprintf(PMFd,"graph g {\n");
        for (i = 0; i < GD_NP; ++i) {
            fprintf(PMFd,"n%d [label=%d];\n",i + 1,gdd_node(i));
            nrec++;
        }
        for (i = 0; i < GD_NP; ++i) {
            /** ii = gdd_node(i); **/

            for (j = 0; j < i; ++j) {
                /** jj = gdd_node(j); **/

                tmp = gdd_adj(i,j,PMGN);
                if (tmp > 0.0) {
                    fprintf(PMFd,"n%d -- n%d;\n",i + 1,j + 1);
                    nrec++;
                }
            }
        }
    }
    else {          /* directed graph */
        fprintf(PMFd,"digraph g {\n");
        for (i = 0; i < GD_NP; ++i) {
            fprintf(PMFd,"n%d [label=%d];\n",i + 1,gdd_node(i));
            nrec++;
        }
        for (i = 0; i < GD_NP; ++i) {
            /** ii = gdd_node(i); **/

            for (j = 0; j < GD_NP; ++j) {
                /** jj = gdd_node(j); **/

                tmp = gdd_adj(i,j,PMGN);
                if (tmp > 0.0) {
                    fprintf(PMFd,"n%d -> n%d;\n",i + 1,j + 1);
                    nrec++;
                }
            }
        }
    }
    fprintf(PMFd,"}\n");
    nrec++;

    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

GDOTFin:       
    p_clean();
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

int gtcl(void)
{
    register int i,j,k;
    int err,typ;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Transitive closure. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto GTCLFin;

    if (gdd_check(PMGN,1))
        goto GTCLFin;
    if (gdd_tcheck(0,0,0))
        goto GTCLFin;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    if (GD_GT == 1 || GD_GT == 3) {         /* unvalued graph */
        typ = 1;
        if (alloc_acc(GD_NP * GD_NP + 1))
            goto GTCLFin;

        g_tclos(PMGN,GD_NP,AcC);
    }
    else {                                  /* valued graph */
        typ = 0;
        if (alloc_acxf(GD_NP * GD_NP + 1))
            goto GTCLFin;

        g_spath(PMGN,GD_NP,AcXF);
    }
    for (i = 0; i < GD_NP; ++i) {

        fprintf(PMFd,PMNFmtS,i + 1);
        fprintf(PMFd,PMNFmtS,gdd_node(i));

        k = i * GD_NP;
        for (j = 0; j < GD_NP; ++j) {
            if (typ)
                fprintf(PMFd,"%d ",(int)AcC[k++]);
            else
                fprintf(PMFd,PMFmtS,(double)AcXF[k++]);
        }
        fprintf(PMFd,"\n");
    }
    printf1("%d records written to: %s\n",GD_NP,PMFdName);
    err = 0;

GTCLFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_tclos(gn,n,ai)                                                        */
/*                                                                          */
/*  Calculate transitive closure of graph gn with type gt. ai is an         */
/*  (n,n) character array, n = number of nodes.                             */

void g_tclos(int gn,int n,char *ai)
{
    register int i,j,k,ni,nk;

    for (i = 0; i < n; ++i) {
        ni = i * n;
        for (j = 0; j < n; ++j) {
            if (i == j || gdd_adj(i,j,gn) >= 0.0)
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
        if (gdd_adj(k,k,gn) < 0.0)
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

void g_spath(int gn,int n,float *a)
{
    register int i,j,k,ni,nj;
    float s,max;

    max = (float)GD_EVMax[gn - 1] * 10.0;

    for (i = 0; i < n; ++i) {
        ni = i * n;
        for (j = 0; j < n; ++j) {
            a[ni] = max;
            if (i != j) {
                s = (float)gdd_adj(i,j,gn);          
                if (s >= 0.0)
                    a[ni] = s;
            }
            ni++;
        }
    }
    for (i = 0; i < n; ++i) {
        ni = i * n;
        for (j = 0; j < n; ++j) {
            nj = j * n;
            if (a[nj + i] < max) {
                for (k = 0; k < n; ++k) {
                    if (a[ni + k] < max) {
                        s = a[nj + i] + a[ni + k];
                        if (s < a[nj + k])
                            a[nj + k] = s;
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
            else if (a[ni] >= max)  
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

int gsort(void)
{
    register int i,j,k,l;
    int err,n,nu,iia,iib,nr,nl;
    double tmp;

    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("Topological sort. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,1,1))       /* get parameters */
        goto GSORTFin;

    if (gdd_check(PMGN,0))
        goto GSORTFin;
    if (gdd_tcheck(1,2,0))
        goto GSORTFin;

    if (PMOPT > 2)
        PMOPT = 2;
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);
           
    if (alloc_acn(GD_NP))           /* in-degree */
        goto GSORTFin;
    if (alloc_acm(GD_NP))   
        goto GSORTFin;
    if (alloc_aci(GD_NP))       
        goto GSORTFin;
    if (alloc_acj(GD_NP))       
        goto GSORTFin;
    if (alloc_acr(GD_NP))       
        goto GSORTFin;
    if (alloc_acx(GD_NP))       
        goto GSORTFin;

    for (i = 0; i < GD_NP; ++i) {

        n = GD_FPN[i];
        if (n <= 0)
            continue;

        for (l = 0; l < n; ++l) {
            k = GD_FPK[i][l];
            if (k >= 0 && gdd_ev(k,PMGN) >= 0.0)   
                AcN[GD_FPI[i][l]] += 1;
        }
    }
    nl = iia = iib = nu = 0;
    for (i = 0; i < GD_NP; ++i) {
        if (AcN[i] == 0) {
            nu++;
            AcI[iib++] = i;
            AcR[i] = nl++;      /* new label */
        }
    }
    while (iia < iib) {
        i = AcI[iia++];
        n = GD_FPN[i];
        for (l = 0; l < n; ++l) {
            k = GD_FPK[i][l];
            if (k >= 0 && gdd_ev(k,PMGN) >= 0.0) { 
                j = GD_FPI[i][l];      /* edge from i to j */
                AcN[j] -= 1;
                if (AcN[j] == 0) {
                    AcI[iib++] = j;
                    AcR[j] = nl++;
                    nu++;
                }
            }
        }
    }   
    if (nu != GD_NP || iib != GD_NP) {
        printf1("Sort not completed. There is at least one cycle or loop.\n");
    }    
    else {
        nr = 0;
        if (PMOPT == 1) {
            for (i = 0; i < GD_NP; ++i) {
                fprintf(PMFd,PMNFmtS,i + 1);
                fprintf(PMFd,PMNFmtS,gdd_node(i));
                fprintf(PMFd,PMNFmtS,AcR[i] + 1);
                fprintf(PMFd,"\n");
                nr++;
            }
        }
        else {
            for (i = 0; i < GD_NP; ++i) {
                nu = 0;
                iia = AcI[i];
                n = GD_FPN[iia];
                for (l = 0; l < n; ++l) {
                    k = GD_FPK[iia][l];
                    if (k >= 0 && (tmp = gdd_ev(k,PMGN)) >= 0.0) { 
                        j = GD_FPI[iia][l];      /* edge from iia to j */
                        AcN[nu] = j;
                        AcM[nu] = AcR[j];
                        AcX[nu++] = tmp;                    
                    }
                }
                if (nu > 0) {

                    if (sortdpi(nu,AcM,AcJ))
                        goto GSORTFin;
    
                    for (l = 0; l < nu; ++l) {
                        j = AcJ[l];        
                        fprintf(PMFd,PMNFmtS,AcR[iia] + 1);
                        fprintf(PMFd,PMNFmtS,AcR[AcN[j]] + 1);
                        fprintf(PMFd,PMNFmtS,gdd_node(iia));
                        fprintf(PMFd,PMNFmtS,gdd_node(AcN[j]));
                        fprintf(PMFd,PMFmtS,AcX[j]);
                        fprintf(PMFd,"\n");
                        nr++;
                    }
                }
            }
        }
        printf1("%d records written to: %s\n",nr,PMFdName);
    }
    err = 0;

GSORTFin:       
    p_clean();
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

int gdcyc(void)
{
    int err,nc,ha,minclen,maxclen;
    int **h;

    ha = 0;
    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Cycles in digraphs. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 5,1,1))       /* get parameters */
        goto GDCYCFin;

    if (gdd_check(PMGN,0))
        goto GDCYCFin;

    if (gdd_tcheck(1,2,0))      /* must be directed with gdd option 1 */
        goto GDCYCFin;
   
    if (PMOPT > 3)
        PMOPT = 3;
    if (PMDOPT != 1)
        PMDOPT = 0;

    if (alloc_acn(GD_NP))
        goto GDCYCFin;
    if (alloc_ack(GD_NP))
        goto GDCYCFin;

    if (!(h = (int **)calloc(GD_NP,sizeof(int *)))) {
        p_err(-2,1);
        goto GDCYCFin; 
    }
    ha = GD_NP;
    memrq(ha,sizeof(int *));

    GSCON_CNT = 0;
    nc = g_dcyc(PMGN,GD_NP,AcN,h,AcK,PMOPT,PMNS,PMMax,&minclen,&maxclen,PMDOPT);
    if (nc < 0)
        goto GDCYCFin;

    printf1("Number of cycles: %d\n",nc);
    if (nc > 0) {
        printf1("Min length of cycles: %d\n",minclen);
        printf1("Max length of cycles: %d\n",maxclen);
    }
    if (PMOPT <= 2)
        printf1("%d records written to: %s\n",GSCON_CNT,PMFdName);
    err = 0;

GDCYCFin:       
    if (ha > 0) {
        free((char *)h);
        memrq(-ha,sizeof(int *));
    }
    p_clean();
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

int g_dcyc(int gn,int n,int *p,int **h,int *hp,int opt,int ns,int cmax,
    int *minclen,int *maxclen,int dopt)
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
        if (k == 0 && GD_BPN[pk] == 0)
            goto CCONT;

        fnd = 0;
        nn = GD_FPN[pk];
        for (i = 0; i < nn; ++i) {

            j = GD_FPI[pk][i];               
            if (GD_FPN[j] == 0)  
                continue;
                 
            ii = GD_FPK[pk][i];
            if (ii >= 0 && gdd_ev(ii,gn) >= 0.0) {  
         
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

        nn = GD_FPN[pk];
        for (i = 0; i < nn; ++i) {
            j = GD_FPI[pk][i];               
            ii = GD_FPK[pk][i];
            if (ii >= 0 && gdd_ev(ii,gn) >= 0.0) {  
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
                fprintf(PMFd,PMNFmtS,nc);
                fprintf(PMFd,PMNFmtS,np);
                for (j = 0; j < np; ++j)
                    fprintf(PMFd,PMNFmtS,gdd_node(p[j]));
                fprintf(PMFd,"\n");
                GSCON_CNT++;
            }
            else if (opt == 2) {
                for (j = 0; j < np; ++j) {
                    fprintf(PMFd,PMNFmtS,nc);
                    fprintf(PMFd,PMNFmtS,np);
                    fprintf(PMFd,PMNFmtS,gdd_node(p[j]));
                    fprintf(PMFd,"\n");
                    GSCON_CNT++;
                }
            }
            if (dopt == 1)
                break;

            ncc = nc / 1000;
            if (ncc * 1000 == nc)
                printf("%10d\n",nc);


        }
CCONT:
        if (k == 0) {
            if (SILENTFlg < 2 && n >= 1000) {
                printfe("Node: %7d     %c",p[0] + 1,CR);
                fflushe();
            }
            if (p[0] >= n - 1)  
                break;      

            p[0] += 1;
            k = 0;
            np = 1;

            for (j = 0; j < n; ++j) {
                if (hp[j] >= 0) {
                    free((char *)h[j]);
                    memrq(-hmax,sizeof(int));
                    hp[j] = -1;
                }
            }
        }
        else {
            pk = p[k];
            if (hp[pk] >= 0) {
                free((char *)h[pk]);
                memrq(-hmax,sizeof(int));
                hp[pk] = -1;
            }
            pk = p[k - 1];

            if (hp[pk] < 0) {
                if (!(h[pk] = (int *)calloc(hmax,sizeof(int)))) {
                    p_err(-2,1);
                    nc = -1;       
                    break;
                }
                memrq(hmax,sizeof(int));
                hp[pk] = 0;
            }
            i = hp[pk];
            if (i >= hmax) {
                printf1("Error: insufficient storage length.\n");
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
    if (SILENTFlg < 2 && n >= 1000) {
        printfe("\n");
        fflushe();
    }
    for (j = 0; j < n; ++j) {
        if (hp[j] >= 0) {
            free((char *)h[j]);
            memrq(-hmax,sizeof(int));
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

int gev(void)
{
    register int i,j;
    int err,r,n,nit;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Eigenvalues/vectors. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto GEVFin;

    if (gdd_check(PMGN,0))
        goto GEVFin;

    if (gdd_tcheck(0,1,0))      /* must be undirected */
        goto GEVFin;
   
    if (PMOPT > 2)
        PMOPT = 2;
    if (PMFmtF == 0)                /* default print format */
        pmfmt(13,6);
    if (MxItFlg == 0)
        MxIter = 30;
    if (PMN < 1)
        PMN = 1;

    if (PMN >= GD_NP) { 
        printf1("Error: n=%d exceeds maximum for current graph.\n",PMN);
        goto GEVFin;
    }
    n = PMN + 1;
    if (alloc_acx(GD_NP * n + 1))
        goto GEVFin;
    if (alloc_acy(GD_NP + 1))
        goto GEVFin;

    printf1("Number of eigenvalues/vectors: %d\n",PMN);
    printf1("Requested accuracy: %g\n",PMEPS);
    printf1("Max number of iterations: %d\n",MxIter);

    r = simitz(GD_NP,n,MxIter,PMEPS,PMN,AcX,AcY,&nit);
    if (r < 0)
        goto GEVFin;

    printf1("Iterations performed: %d\n",nit);
    if (r < 1)
        printf1("No success.\n");
    else {
        printf1("\nEigenvalue(s)\n");
        for (i = 1; i <= r; ++i) {
            printf1("%3d  ",i);
            printf1(PMFmtS,AcY[i]); 
            newline();
        }
        newline();

        if (PMProtFDef)   
            g_evcheck(r,GD_NP,AcX,AcY);

        for (i = 1; i <= GD_NP; ++i) {
            fprintf(PMFd,PMNFmtS,i);
            for (j = 0; j < r; ++j) 
                fprintf(PMFd,PMFmtS,AcX[j * GD_NP + i]);
            fprintf(PMFd,"\n");
        }
        printf1("%d records written to: %s\n",GD_NP,PMFdName);
    }
    err = 0;

GEVFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_op1(n,z,w) returns w = Cz where C is the n x n matrix to be used      */
/*               for the calculation of eigenvalues/vectors.                */
/*               called by simitz().                                        */

void g_op1(int n,double *z,double *w)
{
    register int i,j,l;
    int m;
    double a,tmp;

    if (GD_TYP == 1) {                  /* edge list */

        for (i = 0; i < n; ++i) {
            tmp = 0.0;
            m = GD_FPN[i];
            if (m > 0) {
                for (l = 0; l < m; ++l) {
                    j = GD_FPI[i][l];            /* edge from i to j */
                    a = gdd_adj(i,j,PMGN);
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
                a = gdd_adj(i,j,PMGN);
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

void g_evcheck(int r,int n,double *x,double *d)
{
    register int i,j;
    double tmp,*xp;

    if (alloc_acu(n + 1))
        return;   

    for (i = 1; i <= r; ++i) {
        xp = x + (i - 1) * n;
        g_op1(n,xp,AcU);
        tmp = 0.0;
        for (j = 1; j <= n; ++j)  
            tmp = dmax(tmp,fabs(AcU[j] - d[i] * xp[j]));

        fprintf(PMProtFd,"Check %3d : ",i);
        fprintf(PMProtFd,PMPFmtS,tmp);
        fprintf(PMProtFd,"\n");
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

int gsp(void)
{
    register int i,j;
    int err,n;   
    double tmp;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Shortest paths. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto GSPFin;

    if (gdd_check(PMGN,1))
        goto GSPFin;
    if (gdd_tcheck(1,0,0))
        goto GSPFin;

    if (PMOPT > 3)
        PMOPT = 3;
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    if (alloc_acx(GD_NP))
        goto GSPFin;
    if (alloc_ack(GD_NP + 1))
        goto GSPFin;

    n = 0;
    for (i = 0; i < GD_NP; ++i) {

        g_sp(i,GD_NP,PMGN,AcX,AcK);      

        if (PMOPT <= 2) {
            fprintf(PMFd,PMNFmtS,i + 1);                
            fprintf(PMFd,PMNFmtS,gdd_node(i));                
        }

        for (j = 0; j < GD_NP; ++j) {
            tmp = AcX[j];
            if (PMOPT == 1) {
                if (j == i)
                    continue;
                if (tmp < DBLMAX) {
                    fprintf(PMFd,PMNFmtS,gdd_node(j));                
                    /** fprintf(PMFd," ["); **/            
                    fprintf(PMFd,PMFmtS,tmp);                
                    /** fprintf(PMFd,"] "); **/            
                }
            }
            else {
                if (tmp >= DBLMAX)
                    tmp = -1.0;
                fprintf(PMFd,PMFmtS,tmp);                
            }
        }
        fprintf(PMFd,"\n");                
        n++;
    }
    printf1("%d records written to: %s\n",n,PMFdName);
    err = 0;

GSPFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_glist(i,nodes)    Get list of nodes that are reachable from node i.   */
/*                      Return number of nodes in nodes[]. Skip loops.      */

int g_glist(int i,int *nodes)
{
    register int j,l,n,m;

    m = 0;
    n = GD_FPN[i];  
    for (l = 0; l < n; ++l) {
        j = GD_FPI[i][l];   
        if (j == i)
            continue;
        if (gdd_adj(i,j,PMGN) >= 0.0)
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

void g_sp(int in,int n,int gn,double *d,int *ndec)
{
    register int i,j,l,np;   
    int m,max;
    double dj,tmp;   

    for (i = 0; i < n; ++i) {
        d[i] = DBLMAX;
        ndec[i] = 0;
    }
    max = n + 2;
    d[in] = 0.0;
    ndec[in] = max;  
    np = i = in;

    while (1) {
        m = GD_FPN[i];
        for (l = 0; l < m; ++l) {
            j = GD_FPI[i][l];                   

            tmp = gdd_adj(i,j,PMGN);
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
            gerr_exit(104);
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

int tnet(void)
{
    register int i,j,k,l;
    int err,m,n,ne,iv,is,it,ii,jj,nrec,nrec1,im,n0,nl,ln;
    double tmp,s,t,tmin,tmax;

    err = -1;
    nrec1 = nrec = 0;
    if (check_cmd(1))
        return(-1);

    printf1("Temporal networks. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 4,1,0))       /* get parameters */
        goto TNETFin;

    if (GD_NG < 3) {
        printf1("Error: need multigraph with at least three subgraphs.\n");
        goto TNETFin;
    }
    if (GD_TYP != 1) {
        printf1("Error: need multigraph defined as an edge list.\n");
        goto TNETFin;
    }
    if (PMOPT > 5)
        PMOPT = 1;
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);


    iv = GD_EV[0];
    is = GD_EV[1];
    it = GD_EV[2];

    tmin = DBLMAX;
    tmax = 0.0;
    ne = 0;
    for (i = 0; i < NOC; ++i) {

        tmp = get_data(iv,i);

        if (tmp >= 0.0) {
            ii  = (int)get_data(GD_EI,i);
            jj  = (int)get_data(GD_EJ,i);
            s   = get_data(is,i);
            t   = get_data(it,i);

            if (s < 0.0 || t < s) {
                printf1("Error: invalid dates in edge %d (%d,%d)\n",i+1,ii,jj);
                goto TNETFin;
            }
            tmin = dmin(tmin,s);
            tmax = dmax(tmax,t);
            ne++;
        }
    }
    printf1("\nValid temporal network.\n");
    printf1("Number of nodes: %d\n",GD_NP);
    printf1("Number of edges: %d\n",ne);
    printf1("Time axis: %g to %g\n",tmin,tmax);
    printf1("Interpretation: ");
    if (PMDOPT != 1) {
        PMDOPT = 0;
        printf1("temporal\n\n");
    }
    else
        printf1("historical\n\n");


    if (PMOPT == 1) {
        err = 0;
        goto TNETFin;
    }
    if (PMNTP < 1) {
        printf1("Error: need time points defined with tp option.\n");
        goto TNETFin;
    }
    if (alloc_acn(GD_NP + 1))
        goto TNETFin;
    if (alloc_acm(GD_NP + 1))
        goto TNETFin;
    if (alloc_ack(GD_NP + 1))
        goto TNETFin;

    if (PMOPT == 2) {           /* temporal subgraphs */

        printf1("Option 2: creating temporal/historical subgraphs\n\n");

        prnchar(' ',PMFmt1 - 4,0);
        printf1("Time ");
        printf1("   nodes    edges\n");

        for (j = 0; j < PMNTP; ++j) {
            tmp = PMTP[j];

            for (i = 0; i < GD_NP; ++i)
                AcN[i] = 0;

            m = n = 0; 
            for (i = 0; i < NOC; ++i) {
                s   = get_data(is,i);
                t   = get_data(it,i);
                if (s <= tmp && (PMDOPT == 1 || t >= tmp)) {
                    m++;
                    ii  = (int)get_data(GD_EI,i);
                    k = gdd_fndni(ii);
                    if (k >= 0 && k < GD_NP)
                        AcN[k] += 1;

                    ii  = (int)get_data(GD_EJ,i);
                    k = gdd_fndni(ii);
                    if (k >= 0 && k < GD_NP)
                        AcN[k] += 1;
                }      
            }
            for (k = 0; k < GD_NP;++k) {
                if (AcN[k] > 0)
                    n++;
            }
            printf1(PMFmtS,tmp);
            printf1("%8d %8d\n",n,m);

            if (PMFDef) {
                fprintf(PMFd,PMFmtS,tmp);                
                fprintf(PMFd,"%8d %8d\n",n,m);
                nrec++;
            }
        }
        newline();
        if (PMFDef)
            printf1("%d records written to: %s\n",nrec,PMFdName);

        if (PMF1Def) {

            for (i = 0; i < NOC; ++i) {
                ii  = (int)get_data(GD_EI,i);
                jj  = (int)get_data(GD_EJ,i);
                s   = get_data(is,i);
                t   = get_data(it,i);
                fprintf(PMF1d,"%8d %8d ",ii,jj);

                for (j = 0; j < PMNTP; ++j) {
                    tmp = PMTP[j];
                    if (s <= tmp && (PMDOPT == 1 || t >= tmp)) 
                        fprintf(PMF1d,"  1");
                    else
                        fprintf(PMF1d," -1");
                }
                fprintf(PMF1d,"\n");
                nrec1++;
            }
            printf1("%d records written to: %s\n",nrec1,PMF1dName);
        }
        err = 0;
        goto TNETFin;
    }

    if (PMOPT == 3 || PMOPT == 4) {           /* in/outdegrees */

        printf1("Option %d: distribution of ",PMOPT);
        if (PMOPT == 3)                   
            printf1("in");
        else
            printf1("out");

        printf1("degrees\nInterpretation: ");
        if (GD_GT <= 2)
            printf1("un");
        printf1("directed\n\n");

        prnchar(' ',PMFmt1 - 4,0);
        printf1("Time ");
        printf1("   nodes    edges   frequency of in- or outdegrees (0,1,2,...)\n");

        for (j = 0; j < PMNTP; ++j) {
            tmp = PMTP[j];

            im = 0;
            for (i = 0; i < GD_NP; ++i)
                AcK[i] = AcM[i] = AcN[i] = 0;

            m = n = 0; 
            for (i = 0; i < NOC; ++i) {
                ii  = (int)get_data(GD_EI,i);
                jj  = (int)get_data(GD_EJ,i);
                s   = get_data(is,i);
                t   = get_data(it,i);
                if (s <= tmp && (PMDOPT == 1 || t >= tmp)) {
                    m++;
                    k = gdd_fndni(ii);
                    if (k >= 0 && k < GD_NP) {
                        AcN[k] += 1;
                        if (GD_GT <= 2 || PMOPT == 4)
                            AcM[k] += 1;
                    }
                    k = gdd_fndni(jj);
                    if (k >= 0 && k < GD_NP) {
                        AcN[k] += 1;
                        if (GD_GT <= 2 || PMOPT == 3)
                            AcM[k] += 1;
                    }
                }      
            }
            for (k = 0; k < GD_NP;++k) {
                if (AcN[k] > 0) {
                    n++;
                    l = AcM[k];
                    im = imax(im,l);
                    AcK[l] += 1;
                }
            }
            if (im < 10)
                im = imin(10,GD_NP);

            printf1(PMFmtS,tmp);
            printf1("%8d %8d ",n,m);
            for (k = 0; k <= im; ++k)
                printf1(PMNFmtS,AcK[k]);
            newline();
       
            if (PMFDef) {
                fprintf(PMFd,PMFmtS,tmp);                
                fprintf(PMFd,"%8d %8d ",n,m);
                for (k = 0; k <= im; ++k)
                    fprintf(PMFd,PMNFmtS,AcK[k]);
                fprintf(PMFd,"\n");
                nrec++;
            }
        }
        newline();
        if (PMFDef)
            printf1("%d records written to: %s\n",nrec,PMFdName);

        err = 0;
        goto TNETFin;
    }

    if (PMOPT == 5) {           /* layers */

        printf1("Option 5: layers in historical network.\n\n");
        if (GD_GT < 3) {                 
            printf1("Error: requires directed graph.\n");
            goto TNETFin;
        }
        prnchar(' ',PMFmt1 - 4,0);
        printf1("Time ");
        printf1(" layer  members  node numbers\n");

        for (j = 0; j < PMNTP; ++j) {
            tmp = PMTP[j];

            for (i = 0; i < GD_NP; ++i) {
                AcK[i] = 0;
                AcN[i] = -1;
            }

            /* find members of layer 0 in AcN, number is n0 */

            ln = n0 = 0;
            for (i = 0; i < NOC; ++i) {
                ii  = (int)get_data(GD_EI,i);
                jj  = (int)get_data(GD_EJ,i);
                s   = get_data(is,i);
                if (s <= tmp) {
                    ii  = (int)get_data(GD_EI,i);
                    k = gdd_fndni(ii);
                    if (AcK[k] == 0)
                        AcK[k] = 1;  
                    jj  = (int)get_data(GD_EJ,i);
                    k = gdd_fndni(jj);
                    AcK[k] = 2;  
                }
            }
            for (k = 0; k < GD_NP; ++k) {
                if (AcK[k] == 1) {
                    AcN[k] = 0;
                    n0++;
                }
            }
            if (n0 == 0)
                continue;

            nl = 0;
            for (k = 0; k < GD_NP; ++k) {
                if (AcN[k] == 0)
                    nl++;
            }
            printf1(PMFmtS,tmp);
            printf1("%5d %8d ",ln,nl);

            if (PMFDef) {
                fprintf(PMFd,PMFmtS,tmp);                
                fprintf(PMFd,"%5d %8d ",ln,nl);
            }
            for (k = 0; k < GD_NP; ++k) {
                if (AcN[k] == 0) {
                    printf1(PMNFmtS,GD_ND[k]);
                    if (PMFDef)  
                        fprintf(PMFd,PMNFmtS,GD_ND[k]);
                }

            }
            newline();
            if (PMFDef) {
                fprintf(PMFd,"\n");
                nrec++;
            }
  
            while (1) {
                n0 = 0;
                for (i = 0; i < NOC; ++i) {
                    ii  = (int)get_data(GD_EI,i);
                    jj  = (int)get_data(GD_EJ,i);
                    s   = get_data(is,i);

                    if (s <= tmp) {
                        l = gdd_fndni(ii);
                        if (AcN[l] == ln) {
                            k = gdd_fndni(jj);
                            if (AcN[k] < 0) {
                                AcN[k] = ln + 1;
                                n0++;
                            }
                        }   
                    }
                }
                if (n0 == 0)
                    break;

                ln++;
                nl = 0;
                for (k = 0; k < GD_NP; ++k) {
                    if (AcN[k] == ln)
                        nl++;
                }
                printf1(PMFmtS,tmp);
                printf1("%5d %8d ",ln,nl);
                if (PMFDef) {
                    fprintf(PMFd,PMFmtS,tmp);                
                    fprintf(PMFd,"%5d %8d ",ln,nl);
                }
                for (k = 0; k < GD_NP; ++k) {
                    if (AcN[k] == ln) {
                        printf1(PMNFmtS,GD_ND[k]);
                        if (PMFDef)  
                            fprintf(PMFd,PMNFmtS,GD_ND[k]);
                    }
                }
                newline();
                if (PMFDef) {
                    fprintf(PMFd,"\n");
                    nrec++;
                }
            }
        }
        newline();
        if (PMFDef)
            printf1("%d records written to: %s\n",nrec,PMFdName);

        err = 0;
        goto TNETFin;
    }


TNETFin:       
    p_clean();
    return(err);
}













