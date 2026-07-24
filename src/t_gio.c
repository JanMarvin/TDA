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

/*  functions in t_gio.c */

int gnc(void);
int gcon(void);
void u_visit(int k,int gn);
void u_visit1(int k,int gn);
int gdcon(void);
void d_visit(int k,int gn);
void d_visit1(int k,int gn); 
int gst(void);
void u1_visit(int k,int gn);
void u1_visit1(int k,int gn);
int gmst(void);
void u2_visit(int k,int gn);
void u2_visit1(int k,int gn);
int g_mst(int gn,int n,int *nodes,int *aci,int *acj,int mflag);
int g_mst1(int gn,int gtyp,int n,int *nodes,int *aci,int *acj,int mflag);
int g_ufind(int *dad,int x,int y,int doit);
int gcut(void);
int m_visit(int k,int gn);
int m_visit1(int k,int gn);
int gnst(void);
int gnst_grow(int n,char *g,short *s,char *b,char *a,int ns,int opt);
int gnst_prn(int n,char *a,int opt);
void put_bit(int i,char *bf,int v);
int get_bit(int i,char *bf);
void put_bit2(int i,int j,int n,char *bf,int v);
int get_bit2(int i,int j,int n,char *bf);
void cpy_bit2(int i,char *aci,int j,char *acj,int n);
int gcyc(void);
int g_cyc(int nc,int ne,int n,int *nodes,int opt);
int g_cycprn(int n,int *nodes,int i,int j,int *nptr,int *tptr,int *jptr,
    int ccnt,int nc,int opt,int **cptr,int *cptrn);
void g_cycprn1(int n,int *nodes,int nb,int **cptr,int *cptrn,int nc);
int g_cycprn2(int n,int *nodes,int ne,int nb,int **cptr,int *cptrn,int nc,int opt);
int gio(void);
void visit2(int k,double val,double val1);
int gfcf(void);
int gbcf(void);
void visit3(int k);
int gep(void);
void visit4(int k,int k1,int k2,double val,int first,int opt,int opt1);
int gflow(void);
double g_flow(int n,int i0,int i1,int gn,int *q,int *pre,char *vis,int *fptr,  
    double *flow,double *fd);
double g_flowbfs(int n,int i0,int i1,int gn,int *q,int *pre,char *vis,int *fptr,  
    double *flow,double *fd);
int gfc(void);
double g_cflow(int n,int i0,int i1,int gn,int *q,int *pre,char *vis,int *fptr,  
    double *flow,double *fd,int k0);
double g_cflowbfs(int n,int i0,int i1,int gn,int *q,int *pre,char *vis,int *fptr,  
    double *flow,double *fd,int k0);
int fc_visit(int k,int k0,int gn);


int GSCON_ID = 0;
int GSCON_NN = 0;
int GSCON_NP = 0;
int GSCON_NSP = 0;
int GSCON_CNT = 0;
int GSCON_LEN = 0;
int GSCON_MINLEN = 0;
int GSCON_MAXLEN = 0;
double GSCON_VAL = 0.0;
double GSCON_MINVAL = 0.0;
double GSCON_MAXVAL = 0.0;


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

int gnc(void)
{
    register int i,j,k;
    int err,nrec,n,nn;
    double tmp;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Node-centered networks. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto GNCFin;

    if (gdd_check(PMGN,0))
        goto GNCFin;

    if (gdd_tcheck(0,1,0))      /* must be undirected */
        goto GNCFin;         

    if (PMFmtF == 0)            /* default print format */
        pmfmt(10,4);

    if (alloc_acn(GD_NP))                 
        goto GNCFin;

    nrec = 0;
    for (i = 0; i < GD_NP; ++i) {
        n = 0;
        for (j = 0; j < GD_NP; ++j) {
            if (j == i)
                continue;
            tmp = gdd_adj(i,j,PMGN);
            if (tmp >= 0.0)  
                AcN[n++] = j;
        }   
        nn = n;
        for (j = 0; j < n; ++j) {
            for (k = j + 1; k < n; ++k) {
                if (gdd_adj(AcN[j],AcN[k],PMGN) >= 0.0)
                    nn++;
            }
        }
        n++;
        fprintf(PMFd,PMNFmtS,i + 1);
        fprintf(PMFd,PMNFmtS,gdd_node(i));
        fprintf(PMFd,PMNFmtS,n);
        fprintf(PMFd,PMNFmtS,nn);
        if (n <= 1)
            tmp = 0.0;
        else
            tmp = 2.0 * (double)nn / ((double)n * (double)(n - 1));
        fprintf(PMFd,PMFmtS,tmp);
        fprintf(PMFd,"\n");
        nrec++;
    }






    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

GNCFin:       
    p_clean();
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

int gcon(void)
{
    register int i,j,k;
    int err,nrec;
    double tmp;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Connected components. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto GCONFin;

    if (gdd_check(PMGN,0))
        goto GCONFin;

    if (gdd_tcheck(0,1,0))      /* must be undirected */
        goto GCONFin;         

    if (PMOPT > 2)
        PMOPT = 2;
    if (PMFmtF == 0)            /* default print format */
        pmfmt(10,4);

    if (alloc_acn(GD_NP))                 
        goto GCONFin;

    GSCON_ID = 1;
    i = nrec = 0;

    while (i >= 0) {

        if (SILENTFlg < 2 && GD_NP > 1000) {
            printfe("Component: %7d     %c",GSCON_ID,CR);
            fflushe();
        }
        GSCON_CNT = 0;
        if (GD_TYP == 1)            /* edge list */
            u_visit(i,PMGN);
        else
            u_visit1(i,PMGN);

        i = -1;
        for (j = 0; j < GD_NP; ++j) {
            if (AcN[j] == GSCON_ID) {
                if (PMOPT == 1) {
                    fprintf(PMFd,PMNFmtS,GSCON_ID);
                    fprintf(PMFd,PMNFmtS,GSCON_CNT);
                    fprintf(PMFd,PMNFmtS,j + 1);
                    fprintf(PMFd,PMNFmtS,gdd_node(j));
                    fprintf(PMFd,"\n");
                    nrec++;
                }
                else if (GSCON_CNT == 1) {
                    fprintf(PMFd,PMNFmtS,GSCON_ID);
                    fprintf(PMFd,PMNFmtS,GSCON_CNT);
                    fprintf(PMFd,PMNFmtS,j + 1);
                    fprintf(PMFd,PMNFmtS,j + 1);
                    fprintf(PMFd,PMNFmtS,gdd_node(j));
                    fprintf(PMFd,PMNFmtS,gdd_node(j));
                    fprintf(PMFd,PMFmtS,gdd_adj(j,j,PMGN));
                    fprintf(PMFd,"\n");
                    nrec++;
                }
                else {
                    for (k = j + 1; k < GD_NP; ++k) {
                        if (AcN[k] == GSCON_ID) {
                            tmp = gdd_adj(j,k,PMGN);
                            if (tmp >= 0.0) {
                                fprintf(PMFd,PMNFmtS,GSCON_ID);
                                fprintf(PMFd,PMNFmtS,GSCON_CNT);
                                fprintf(PMFd,PMNFmtS,j + 1);
                                fprintf(PMFd,PMNFmtS,k + 1);
                                fprintf(PMFd,PMNFmtS,gdd_node(j));
                                fprintf(PMFd,PMNFmtS,gdd_node(k));
                                fprintf(PMFd,PMFmtS,tmp);
                                fprintf(PMFd,"\n");
                                nrec++;
                            }
                        }
                    }
                }
            }
            else if (i < 0 && AcN[j] == 0)
                i = j;
        }
        GSCON_ID++;  
    }
    if (SILENTFlg < 2 && GD_NP > 1000) {
        printfe("\n");
        fflushe();
    }
    printf1("Number of components: %d\n",GSCON_ID - 1);
    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

GCONFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  u_visit(k,gn)   called by gcon. k is node number, gn is graph number.   */

void u_visit(int k,int gn)
{
    register int j,l,k1,k2;
    int n;

    AcN[k] = GSCON_ID;
    GSCON_CNT++;

    n = GD_FPN[k];
    if (n <= 0)
        return;          

    for (l = 0; l < n; ++l) {
        
        j = GD_FPI[k][l];                    /* edge from k to j */
        if (AcN[j] == 0) {
            k1 = GD_FPK[k][l];
            k2 = GD_FPK1[k][l];
            if ((k1 >= 0 && gdd_ev(k1,gn) >= 0.0) ||
                (k2 >= 0 && gdd_ev(k2,gn) >= 0.0))    
                u_visit(j,gn);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  u_visit1(k,gn)  called by gcon. k is node number, gn is graph number.   */

void u_visit1(int k,int gn)
{
    register int j;

    AcN[k] = GSCON_ID;
    GSCON_CNT++;

    for (j = 0; j < GD_NP; ++j) {

        if (AcN[j] == 0) {
            if (gdd_adj(k,j,gn) >= 0.0)
                u_visit1(j,gn);
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

int gdcon(void)
{
    register int i,j,k;
    int err,nrec,n,ni,maxl,nc;
    double v;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Search for reachable nodes. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 5,1,1))       /* get parameters */
        goto GDCONFin;

    if (gdd_check(PMGN,0))
        goto GDCONFin;

    if (gdd_tcheck(0,2,0))      /* must be directed */
        goto GDCONFin;

    if (PMOPT > 4)
        PMOPT = 4;

    if (PMFmtF == 0)            /* default print format */
        pmfmt(10,4);

    ni = 0;
    if (PMOPT <= 2 && PMF2Def) {
        if (alloc_ack(GD_NP))                 
            goto GDCONFin;

        if ((ni = g_getnd(0,AcK)) < 0)         /* get node list */
            goto GDCONFin;

        if (ni > 0)
            printf1("Number of input nodes: %d\n",ni);

    }
    if (alloc_acn(GD_NP))                 
        goto GDCONFin;

    nrec = 0;

    if (PMOPT <= 2) {

        GSCON_NN = 0;
        GSCON_ID = 1;
        maxl = 0;
        for (i = 0; i < GD_NP; ++i) {

            if (ni > 0 && AcK[i] == 0)
                continue;

            for (j = 0; j < GD_NP; ++j)
                AcN[j] = 0;

            if (GD_TYP == 1)
                d_visit(i,PMGN);
            else
                d_visit1(i,PMGN);

            n = -1;
            for (j = 0; j < GD_NP; ++j) {
                if (AcN[j])
                    n++;
            }
            if (n > 0) {
                fprintf(PMFd,PMNFmtS,i + 1);
                fprintf(PMFd,PMNFmtS,gdd_node(i));
                fprintf(PMFd,PMNFmtS,n);
                if (PMOPT == 2) {
                    for (j = 0; j < GD_NP; ++j) {
                        if (AcN[j] && j != i)
                            fprintf(PMFd,PMNFmtS,gdd_node(j));
                    }
                }
                fprintf(PMFd,"\n");
                nrec++;
                GSCON_ID += 1;
                if (maxl < n)
                    maxl = n;
            }
            n_info(i + 1,100);
        }
        n_info_e();
        printf1("Max number of reachable nodes: %d\n",maxl);
        printf1("%d records written to: %s\n",nrec,PMFdName);
    }
    else {                          /* directed graph */

        if (alloc_acm(GD_NP))                 
            goto GDCONFin;
        if (alloc_ack(GD_NP))                 
            goto GDCONFin;

        GSCON_NN = 0;
        GSCON_ID = 1;
        nc = 0;
                      
        for (i = 0; i < GD_NP; ++i) {

            if (AcK[i])
                continue;

            for (j = 0; j < GD_NP; ++j)            
                AcN[j] = 0;            
             
            if (GD_TYP == 1)
                d_visit(i,PMGN);
            else
                d_visit1(i,PMGN);
    
            /* save all node numbers that can be reached from i in AcM */
    
            for (j = 0; j < GD_NP; ++j)            
                AcM[j] = AcN[j];
    
            nc++;
            AcK[i] = nc;

            for (j = 0; j < GD_NP; ++j) {
                if (j != i && AcM[j]) {
                    for (k = 0; k < GD_NP; ++k)
                        AcN[k] = 0;

                    if (GD_TYP == 1)
                        d_visit(j,PMGN);
                    else
                        d_visit1(j,PMGN);

                    if (AcN[i])
                        AcK[j] = nc;
                }
            }
            n_info(i + 1,100);
        }
        n_info_e();
        printf1("Number of strong components: %d\n",nc);

        if (PMOPT == 3) {
            for (i = 0; i < GD_NP; ++i) {
                fprintf(PMFd,PMNFmtS,AcK[i]);
                fprintf(PMFd,PMNFmtS,i + 1);
                fprintf(PMFd,PMNFmtS,gdd_node(i));
                fprintf(PMFd,"\n");
                nrec++;
            }
        }
        else {
            for (k = 1; k <= nc; ++k) {
                for (i = 0; i < GD_NP; ++i) {
                    if (AcK[i] != k)
                        continue;
                    n = 0;
                    for (j = 0; j < GD_NP; ++j) {
                        if (AcK[j] == k && (v = gdd_adj(i,j,PMGN)) >= 0.0) {
                            fprintf(PMFd,PMNFmtS,k);
                            fprintf(PMFd,PMNFmtS,i + 1);
                            fprintf(PMFd,PMNFmtS,j + 1);
                            fprintf(PMFd,PMNFmtS,gdd_node(i));
                            fprintf(PMFd,PMNFmtS,gdd_node(j));
                            fprintf(PMFd,PMFmtS,v);
                            fprintf(PMFd,"\n");
                            nrec++;
                            n++;
                        }
                    }
                    if (n == 0) {
                        fprintf(PMFd,PMNFmtS,k);
                        fprintf(PMFd,PMNFmtS,i + 1);
                        fprintf(PMFd,PMNFmtS,i + 1);
                        fprintf(PMFd,PMNFmtS,gdd_node(i));
                        fprintf(PMFd,PMNFmtS,gdd_node(i));
                        fprintf(PMFd,PMFmtS,-1.0);
                        fprintf(PMFd,"\n");
                        nrec++;
                    }
                }
            }
        }
        printf1("%d records written to: %s\n",nrec,PMFdName);
    }
    err = 0;

GDCONFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  d_visit(k,gn)   directed visit, called by gdcon. k is current node,     */
/*                  gn is graph number.                                     */

void d_visit(int k,int gn)
{
    register int j,l,kp;
    int n;

    AcN[k] = GSCON_ID;

    n = GD_FPN[k];
    if (n <= 0)
        return;          

    for (l = 0; l < n; ++l) {

        j = GD_FPI[k][l];               
        if (AcN[j] == 0) {
            kp = GD_FPK[k][l];
            if (kp >= 0 && gdd_ev(kp,gn) >= 0.0)    
                d_visit(j,gn);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  d_visit1(k,gn)  directed visit, called by gdcon. k is current node,     */
/*                  gn is graph number.                                     */

void d_visit1(int k,int gn)
{
    register int j;

    AcN[k] = GSCON_ID;

    for (j = 0; j < GD_NP; ++j) {

        if (AcN[j] == 0) {
            if (gdd_adj(k,j,gn) >= 0.0)
                d_visit1(j,gn);
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

int gst(void)
{
    register int i,j;
    int err,nrec;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Spanning tree. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto GSTFin;

    if (gdd_check(PMGN,0))
        goto GSTFin;

    if (gdd_tcheck(0,1,0))      /* must be undirected */
        goto GSTFin;

    if (PMFmtF == 0)            /* default print format */
        pmfmt(10,4);

    if (alloc_acn(GD_NP))                 
        goto GSTFin;
    if (alloc_aci(GD_NP))                 
        goto GSTFin;
    if (alloc_acj(GD_NP))                 
        goto GSTFin;

    GSCON_ID = 1;
    i = nrec = 0;

    while (i >= 0) {

        if (SILENTFlg < 2 && GD_NP > 500) {
            printfe("Tree: %7d     %c",GSCON_ID,CR);
            fflushe();
        }
        GSCON_NN = 0;
        if (GD_TYP == 1)
            u1_visit(i,PMGN);
        else
            u1_visit1(i,PMGN);

        for (j = 0; j < GSCON_NN; ++j) {
            fprintf(PMFd,PMNFmtS,GSCON_ID);
            fprintf(PMFd,PMNFmtS,gdd_node(AcI[j]));
            fprintf(PMFd,PMNFmtS,gdd_node(AcJ[j]));
            fprintf(PMFd,PMFmtS,gdd_adj(AcI[j],AcJ[j],PMGN)); 
            fprintf(PMFd,"\n");
            nrec++;
        }
        i = -1;
        for (j = 0; j < GD_NP; ++j) {
            if (AcN[j] == 0) {
                i = j;
                break;
            }
        }
        GSCON_ID++;  
    }
    if (SILENTFlg < 2 && GD_NP > 500) {
        printfe("\n");
        fflushe();
    }
    printf1("Number of trees: %d\n",GSCON_ID - 1);
    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

GSTFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  u1_visit(k,gn)   called by gcon. k is node number, gn is graph number.   */

void u1_visit(int k,int gn)
{
    register int j,l,k1,k2;
    int n;

    AcN[k] = GSCON_ID;

    n = GD_FPN[k];
    if (n <= 0)
        return;          

    for (l = 0; l < n; ++l) {

        j = GD_FPI[k][l];                    /* edge from k to j */
        if (AcN[j] == 0) {
            k1 = GD_FPK[k][l];
            k2 = GD_FPK1[k][l];
            if ((k1 >= 0 && gdd_ev(k1,gn) >= 0.0) ||
                (k2 >= 0 && gdd_ev(k2,gn) >= 0.0)) {  
                AcI[GSCON_NN] = k;
                AcJ[GSCON_NN++] = j;
                u1_visit(j,gn);
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  u1_visit1(k,gn)  called by gcon. k is node number, gn is graph number.   */

void u1_visit1(int k,int gn)
{
    register int j;

    AcN[k] = GSCON_ID;

    for (j = 0; j < GD_NP; ++j) {

        if (AcN[j] == 0) {
            if (gdd_adj(k,j,gn) >= 0.0) {
                AcI[GSCON_NN] = k;
                AcJ[GSCON_NN++] = j;
                u1_visit1(j,gn);
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

int gmst(void)
{
    register int i,j,k;
    int err,nrec,ne,nt,nte;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Minimum (maximum) spanning tree. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto GMSTFin;

    if (gdd_check(PMGN,0))
        goto GMSTFin;

    if (gdd_tcheck(0,1,2))      /* must be undirected */
        goto GMSTFin;

    if (PMFmtF == 0)            /* default print format */
        pmfmt(10,4);

    if (PMMax != 1)
        PMMax = 0;

    if (PMALG != 2)
        PMALG = 1;

    printf1("Algorithm: ");
    if (PMALG == 1)
        printf1("Kruskal.\n");
    else  
        printf1("Prim.\n");
         
    if (alloc_acn(GD_NP))                 
        goto GMSTFin;
    if (alloc_acm(GD_NP))                 
        goto GMSTFin;
    if (alloc_aci(GD_NP))                 
        goto GMSTFin;
    if (alloc_acj(GD_NP))                 
        goto GMSTFin;

    GSCON_ID = 1;
    nte = nt = i = nrec = 0;

    while (i >= 0) {

        if (SILENTFlg < 2 && GD_NP > 500) {
            printfe("Tree: %7d     %c",GSCON_ID,CR);
            fflushe();
        }
        GSCON_NN = 0;
        if (GD_TYP == 1)
            u2_visit(i,PMGN);
        else
            u2_visit1(i,PMGN);
     
        if (GSCON_NN == 0)
            break;

        if (GSCON_NN == 1) {
            AcI[0] = AcJ[0] = AcM[0];
            ne = 1;
            nte++;
        }
        else {
            if (PMALG == 1)
                ne = g_mst1(PMGN,GD_TYP,GSCON_NN,AcM,AcI,AcJ,PMMax);
            else
                ne = g_mst(PMGN,GSCON_NN,AcM,AcI,AcJ,PMMax);         

            if (ne <= 0) {
                if (ne == 0)  
                    printf1("ERROR: gmst.\n");
                goto GMSTFin;
            }   
        }
        if (PMSORTFlg) {
            if (alloc_ack(ne))
                goto GMSTFin;
            if (sortdpi2(ne,AcI,AcJ,AcK))
                goto GMSTFin;
        }
        for (j = 0; j < ne; ++j) {
            if (PMSORTFlg)
                k = AcK[j];
            else
                k = j;

            fprintf(PMFd,PMNFmtS,GSCON_ID);
            fprintf(PMFd,PMNFmtS,gdd_node(AcI[k]));
            fprintf(PMFd,PMNFmtS,gdd_node(AcJ[k]));
            fprintf(PMFd,PMFmtS,gdd_adj(AcI[k],AcJ[k],PMGN)); 
            fprintf(PMFd,"\n");
            nrec++;
        }
        nt++;
        GSCON_ID++;  

        i = -1;
        for (j = 0; j < GD_NP; ++j) {
            if (AcN[j] == 0) {
                i = j;
                break;
            }
        }
    }
    if (SILENTFlg < 2 && GD_NP > 500) {
        printfe("\n");
        fflushe();
    }
    printf1("Number of trees: %d",nt);
    if (nte > 0)
        printf1(" [includes %d isolated nodes]",nte);
    printf1("\n%d records written to: %s\n",nrec,PMFdName);
    err = 0;

GMSTFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  u2_visit(k,gn)   called by gcon. k is node number, gn is graph number.   */

void u2_visit(int k,int gn)
{
    register int j,l,k1,k2;
    int n;

    AcN[k] = GSCON_ID;
    AcM[GSCON_NN++] = k;

    n = GD_FPN[k];
    if (n <= 0)
        return;          

    for (l = 0; l < n; ++l) {

        j = GD_FPI[k][l];                    /* edge from k to j */
        if (AcN[j] == 0) {
            k1 = GD_FPK[k][l];
            k2 = GD_FPK1[k][l];
            if ((k1 >= 0 && gdd_ev(k1,gn) >= 0.0) ||
                (k2 >= 0 && gdd_ev(k2,gn) >= 0.0)) {  
                u2_visit(j,gn);
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  u2_visit1(k,gn)  called by gcon. k is node number, gn is graph number.   */

void u2_visit1(int k,int gn)
{
    register int j;

    AcN[k] = GSCON_ID;
    AcM[GSCON_NN++] = k;

    for (j = 0; j < GD_NP; ++j) {

        if (AcN[j] == 0) {
            if (gdd_adj(k,j,gn) >= 0.0)  
                u2_visit1(j,gn);
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

int g_mst(int gn,int n,int *nodes,int *aci,int *acj,int mflag)
{
    register int i,im,u,w;
    int ne;
    double tmp,mval,cmin;

    if (alloc_acx(n))
        return(-1);
    if (alloc_acr(n))
        return(-1);

    if (mflag == 0)
        mval = GD_EVMax[gn - 1] + 1.0;
    else
        mval = -1.0;

    u = nodes[0];
    nodes[0] = -1;
    im = -1;
    cmin = mval;
    for (i = 1; i < n; ++i) {
        tmp = gdd_adj(u,nodes[i],gn);
        if (tmp < 0.0)
            tmp = mval;
        AcX[i] = tmp;   
        AcR[i] = u;
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
        aci[ne] = imin(AcR[im],w);
        acj[ne++] = imax(AcR[im],w);          

        im = -1;
        cmin = mval;
        for (i = 0; i < n; ++i) {
            if (nodes[i] >= 0) {
                u = nodes[i];
                tmp = gdd_adj(w,u,gn);
                if (tmp >= 0.0) {
                    if (mflag == 0 && tmp < AcX[i]) {
                        AcR[i] = w;
                        AcX[i] = tmp;
                    }
                    else if (mflag == 1 && tmp > AcX[i]) {
                        AcR[i] = w;
                        AcX[i] = tmp;
                    }
                }
                if (mflag == 0) {
                    if (cmin > AcX[i]) {
                        cmin = AcX[i];
                        im = i;
                    }
                }
                else {
                    if (cmin < AcX[i]) {
                        cmin = AcX[i];
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

int g_mst1(int gn,int gtyp,int n,int *nodes,int *aci,int *acj,int mflag)
{
    register int i,j,l,ip;
    int err,ne,nemax,nn,da;
    int *dad;
    double tmp;

    err = -1;   

    nemax = GD_NE[gn - 1];
    da = 0;
    if (!(dad = (int *)calloc(GD_NP + 2,sizeof(int)))) {
        p_err(-2,1);
        return(-1);
    }
    da = GD_NP + 2;
    memrq(da,sizeof(int));

    if (alloc_acr(nemax))
        goto MST1Fin;
    if (alloc_acs(nemax))
        goto MST1Fin;
    if (alloc_acx(nemax))
        goto MST1Fin;

    ne = 0;
    if (gtyp == 1) {                        /* edge list */
        for (i = 0; i < n; ++i) {
            ip = nodes[i];
            nn = GD_FPN[ip];
            for (l = 0; l < nn; ++l) {
                j = GD_FPI[ip][l];           /* edge from ip to j */
                if (j < ip) {
                    tmp = gdd_adj(ip,j,gn);
                    if (tmp >= 0.0) {
                        AcR[ne] = ip;
                        AcS[ne] = j;
                        AcX[ne++] = tmp;
                    }
                }
            }
        }
    }
    else {
        for (i = 1; i < n; ++i) {
            for (j = 0; j < i; ++j) {
                tmp = gdd_adj(nodes[i],nodes[j],gn);
                if (tmp >= 0.0) {
                    AcR[ne] = nodes[i];
                    AcS[ne] = nodes[j];
                    AcX[ne++] = tmp;
                }
            }
        }
    }
    if (alloc_ack(ne + 1))
        goto MST1Fin;

    if (sortdp(ne,AcX,AcK))
        goto MST1Fin;

    nn = 0;
    if (mflag == 0) {
        for (i = 0; i < ne; ++i) {
            ip = AcK[i];
            j = AcR[ip];
            l = AcS[ip];
                          
            if (g_ufind(dad,j + 1,l + 1,1)) {
                aci[nn] = imin(j,l);
                acj[nn++] = imax(j,l);
            }
        }
    }
    else {
        for (i = ne - 1; i >= 0; i--) {
            ip = AcK[i];
            j = AcR[ip];
            l = AcS[ip];
                          
            if (g_ufind(dad,j + 1,l + 1,1)) {
                aci[nn] = imin(j,l);
                acj[nn++] = imax(j,l);
            }
        }
    }
    err = nn;

MST1Fin:
    if (da > 0) {
        free((char *)dad);
        memrq(-da,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_ufind(dad,x,y,doit)   union find function, see Sedgewick 1990, p.446  */

int g_ufind(int *dad,int x,int y,int doit)
{
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

int gcut(void)
{
    register int i,j;
    int err,nrec,nb,nc,ncp,n,nid;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Cut nodes and blocks. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto GCUTFin;

    if (gdd_check(PMGN,0))
        goto GCUTFin;

    if (gdd_tcheck(0,1,0))      /* must be undirected */
        goto GCUTFin;
   
    if (PMOPT > 2)
        PMOPT = 2;

    if (alloc_acn(GD_NP))                 
        goto GCUTFin;
    if (alloc_acm(GD_NP))                 
        goto GCUTFin;
    if (alloc_acns(GD_NP))                 
        goto GCUTFin;

    GSCON_ID = 0;
    nb = ncp = nc = i = nrec = 0;

    while (i >= 0) {

        if (SILENTFlg < 2 && GD_NP > 1000) {
            printfe("Component: %7d     %c",nc+1,CR);
            fflushe();
        }

        GSCON_NSP = i;
        GSCON_CNT = 0;              /* number of paths from root */
        GSCON_NN  = 0;              /* number of elements in component */
        nid = GSCON_ID;
        if (GD_TYP == 1)
            m_visit(i,PMGN);
        else
            m_visit1(i,PMGN);
    
        nc++;
        if (GSCON_CNT < 2)
            AcNS[i] = 0;

        n = 0;
        for (j = 0; j < GD_NP; ++j) {
            if (AcNS[j]) {
                if (PMOPT == 1) {
                    AcM[n] = j;
                    AcNS[j] = 0;
                }
                n++;
                ncp++;
            }
        }
        if (n == 0)
            nb++;

        if (PMOPT == 1) {
            fprintf(PMFd,PMNFmtS,nc);
            fprintf(PMFd,PMNFmtS,GSCON_NN);
            fprintf(PMFd,PMNFmtS,n);
            for (j = 0; j < n; ++j)
                fprintf(PMFd,PMNFmtS,gdd_node(AcM[j]));
            fprintf(PMFd,"\n");
            nrec++;
        }
        i = -1;
        for (j = 0; j < GD_NP; ++j) {
            if (AcN[j] > nid) {
                if (PMOPT == 2) {
                    fprintf(PMFd,PMNFmtS,nc);
                    fprintf(PMFd,PMNFmtS,j + 1);
                    fprintf(PMFd,PMNFmtS,gdd_node(j));
                    fprintf(PMFd,PMNFmtS,n);
                    fprintf(PMFd,PMNFmtS,AcNS[j]);
                    fprintf(PMFd,"\n");
                    nrec++;

                    if (AcNS[j] > 0)  
                        AcNS[j] = 0;
                }
            }
            else if (i < 0 && AcN[j] == 0)
                i = j;
        }
        GSCON_ID++;  
    }
    if (SILENTFlg < 2 && GD_NP > 1000) {
        printfe("\n");
        fflushe();
    }
    printf1("Number of components: %d\n",nc);
    printf1("Number of cut nodes: %d\n",ncp);
    printf1("Number of blocks: %d\n",nb);
    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

GCUTFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  m_visit(k,gn)   called by gcut. k is node number, gn is graph number.   */
/*                  GSCON_NSP = k in first call.                            */

int m_visit(int k,int gn)
{
    register int j,l,k1,k2;
    int n,m,min;

    AcN[k] = ++GSCON_ID;
    min = GSCON_ID;
    GSCON_NN++;
   
    n = GD_FPN[k];
    for (l = 0; l < n; ++l) {

        j = GD_FPI[k][l];                    /* edge from k to j */
        k1 = GD_FPK[k][l];
        k2 = GD_FPK1[k][l];
        if ((k1 >= 0 && gdd_ev(k1,gn) >= 0.0) ||
                      (k2 >= 0 && gdd_ev(k2,gn) >= 0.0)) {  

            if (AcN[j] == 0) {
               if (k == GSCON_NSP)
                   GSCON_CNT++;

                m = m_visit(j,gn);
                if (m < min)     
                    min = m;
     
                if (m >= AcN[k])  
                    AcNS[k] = 1;
            }
            else if (AcN[j] < min)   
                min = AcN[j];
        }
    }
    return(min);
}

/* ------------------------------------------------------------------------ */
/*  m_visit1(k,gn)  called by gcut. k is node number, gn is graph number.   */
/*                  GSCON_NSP = k in first call.                            */

int m_visit1(int k,int gn)
{
    register int j;
    int m,min;

    AcN[k] = ++GSCON_ID;
    min = GSCON_ID;
    GSCON_NN++;

    for (j = 0; j < GD_NP; ++j) {
        if (gdd_adj(k,j,gn) >= 0.0) {
            if (AcN[j] == 0) {
                if (k == GSCON_NSP)
                    GSCON_CNT++;

                m = m_visit1(j,gn);

                if (m < min)     
                    min = m;
     
                if (m >= AcN[k])    
                    AcNS[k] = 1;
            }
            else if (AcN[j] < min)   
                min = AcN[j];
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

int gnst(void)
{
    register int i,j,ii;
    int err,n,ni,nc,in,jn;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Enumeration of spanning trees. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto GNSTFin;

    if (gdd_check(PMGN,0))
        goto GNSTFin;

    if (gdd_tcheck(0,1,0))      /* must be undirected */
        goto GNSTFin;
   
    if (PMOPT > 2)
        PMOPT = 2;

    if (alloc_acn(GD_NP))                 
        goto GNSTFin;
    if (alloc_acm(GD_NP))                 
        goto GNSTFin;

    GSCON_NN = 0;
    GSCON_ID = 1;
    nc = ni = ii = 0;

    printf1("\nComponent  nodes    trees\n");

    while (ii >= 0) {
        GSCON_CNT = 0;
        if (GD_TYP == 1)
            u_visit(ii,PMGN);
        else
            u_visit1(ii,PMGN);

        ii = -1;
        n = 0;
        for (j = 0; j < GD_NP; ++j) {
            if (AcN[j] == GSCON_ID)  
                AcM[n++] = j;
            else if (ii < 0 && AcN[j] == 0)
                ii = j;
        }
        if (n == 1)             /* isolated node */
            ni++;

        if (n > 0) {            /* new component */

            GSCON_CNT = 0;
            if (SILENTFlg < 2 && GD_NP > 10) {
                printfe("Component: %7d     %c",GSCON_ID,CR);
                fflushe();
            }
            if (n > 1) {
                if (alloc_acc(n * n / 8 + 1))
                    goto GNSTFin;
                if (alloc_ace(n * n / 8 + 1))
                    goto GNSTFin;
                if (alloc_acd(n / 8 + 1))
                    goto GNSTFin;
                if (alloc_acns(n))
                    goto GNSTFin;

                for (i = 1; i < n; ++i) {
                    in = AcM[i];
                    for (j = 0; j < i; ++j) {
                        jn = AcM[j];
                        if (gdd_adj(in,jn,PMGN) >= 0.0) {
                            put_bit2(i,j,n,AcC,1);
                            put_bit2(j,i,n,AcC,1);
                        }
                    }
                }
                AcNS[0] = 1;
                for (j = 0; j < n; ++j) {
                    if (get_bit2(0,j,n,AcC))
                        put_bit(j,AcD,1);
                }   
                if (gnst_grow(n,AcC,AcNS,AcD,AcE,1,PMOPT))
                    goto GNSTFin;      
            }
            printf1("%6d %9d %8d\n",++nc,n,GSCON_CNT);
        }
        GSCON_ID++;  
    }
    if (SILENTFlg < 2 && GD_NP > 10) {
        printfe("\n");
        fflushe();
    }
    printf1("\n%d records written to: %s\n",GSCON_NN,PMFdName);
    err = 0;

GNSTFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gnst_grow   find all spanning trees for an undirected graph.            */
/*              the algorithm is adapted from:                              */
/*              M.D. McIlroy, Generator of Spanning Trees, CACM 354, 1969.  */
/*                                                                          */
/*  Return 0 if OK, -1 if error (insuff. memory)                            */

int gnst_grow(int n,char *g,short *s,char *b,char *a,int ns,int opt)
{
    register int i,j,k;
    int err,t,nfin,gfa,bfa;
    char *gf,*bf;

    err = gfa = bfa = 0;
     
    if (ns == n) {
        if (gnst_prn(n,a,opt))      /* print */
            return(-1);
    }
    else {
        if (!(gf = (char *)calloc(n * n / 8 + 1,sizeof(char)))) {
            p_err(-2,1);
            return(-1);
        }
        gfa = n * n / 8 + 1;
        memrq(gfa,sizeof(char));
        memcpy(gf,g,gfa);

        for (i = 0; i < n; ++i) {
            if (get_bit(i,b)) {
                nfin = 0;
                for (j = 0; j < n; ++j) {
                    k = get_bit2(i,j,n,gf);
                    put_bit2(i,j,n,a,k * s[j]);

                    k *= (1 - s[j]);
                    put_bit2(i,j,n,gf,k);
                    nfin += k;
                }
                t = s[i]; 
                if (t == 0) {
                    s[i] = 1;
                    ns++;
                }
                if (!(bf = (char *)calloc(n / 8 + 1,sizeof(char)))) {
                    p_err(-2,1);
                    err = -1;
                    break;           
                }
                bfa = n / 8 + 1;  
                memrq(bfa,sizeof(char));
                     
                for (j = 0; j < n; ++j) {
                    if (get_bit2(i,j,n,gf) || (j > i && get_bit(j,b)))
                        put_bit(j,bf,1);
                }
                err = gnst_grow(n,gf,s,bf,a,ns,opt);

                free((char *)bf);
                memrq(-bfa,sizeof(char));

                s[i] = t;
                if (t == 0)
                    ns--;
                if (err || nfin == 0)
                    break;         
            }
        }
    }
    if (gfa > 0) {
        free((char *)gf);
        memrq(-gfa,sizeof(char));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gnst_prn    Print trees created by gnst_grow. Print to PMFd.            */
/*              Count number of trees in GSCON_CNT, records in GSCON_NN.    */
/*              Use AcM[] for translation of node numbers.                  */

int gnst_prn(int n,char *a,int opt)
{
    register int i,j;
    int fin,m;

    if (alloc_acms(n))
        return(-1);    
    if (alloc_ack(n))
        return(-1);    

    m = 0;
    fin = 1;
    for (i = 0; i < n; ++i) {
        AcK[i] = -1;
        for (j = 0; j < n; ++j) {
            if (get_bit2(i,j,n,a)) {
                AcMS[i] = AcK[i] = j;
                m++;
                fin = 0;
                break;
            }
        }
    }
    while (fin == 0) {
        GSCON_CNT++;
        if (opt == 1) {
            fprintf(PMFd,PMNFmtS,GSCON_ID);                
            fprintf(PMFd,PMNFmtS,GSCON_CNT);                
            for (j = 0; j < n; ++j) {
                i = AcMS[j];
                if (AcK[j] >= 0)
                    i = gdd_node(AcM[i]);
                fprintf(PMFd,PMNFmtS,i);
            }
            fprintf(PMFd,"\n");                
            GSCON_NN++;
        }
        else {
            i = 1;
            for (j = 0; j < n; ++j) {
                if (AcK[j] >= 0) {
                    fprintf(PMFd,PMNFmtS,GSCON_ID);                
                    fprintf(PMFd,PMNFmtS,GSCON_CNT);                
                    fprintf(PMFd,PMNFmtS,m);                
                    fprintf(PMFd,PMNFmtS,i++);                
                    fprintf(PMFd,PMNFmtS,gdd_node(AcM[j]));
                    fprintf(PMFd,PMNFmtS,gdd_node(AcM[(int)AcMS[j]]));
                    fprintf(PMFd,"\n");                
                    GSCON_NN++;
                }
            }
        }
        fin = 1;
        for (i = n - 1; i >= 0; --i) {
            if (AcK[i] < 0)
                continue;

            for (j = AcMS[i] + 1; j < n; ++j) {
                if (get_bit2(i,j,n,a)) {
                    AcMS[i] = j;
                    fin = 0;
                    break;
                }
            }
            if (fin == 0)
                break;
            AcMS[i] = AcK[i];
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  put_bit(i,bf,v)                                                         */

void put_bit(int i,char *bf,int v)
{
    register char *p;

    p = bf + i / 8;
    if (v)
        *p |= BMsk[i % 8];
    else
        *p &= ~BMsk[i % 8];
}

/* ------------------------------------------------------------------------ */
/*  get_bit(i,bf)                                                           */

int get_bit(int i,char *bf)
{
    register char *p;

    p = bf + i / 8;
    if (*p & BMsk[i % 8])
        return(1);
    else
        return(0);
}

/* ------------------------------------------------------------------------ */
/*  put_bit2(i,j,n,bf,v)                                                    */

void put_bit2(int i,int j,int n,char *bf,int v)
{
    register char *p;

    i = i * n + j;
    p = bf + i / 8;
    if (v)
        *p |= BMsk[i % 8];
    else
        *p &= ~BMsk[i % 8];
}

/* ------------------------------------------------------------------------ */
/*  get_bit2(i,j,n,bf)                                                      */

int get_bit2(int i,int j,int n,char *bf)
{
    register char *p;

    i = i * n + j;
    p = bf + i / 8;
    if (*p & BMsk[i % 8])
        return(1);
    else
        return(0);
}

/* ------------------------------------------------------------------------ */
/*  cpy_bit2(i,aci,j,acj,n)                                                 */
/*                                                                          */
/*  Copy row j in acj to row i in aci, n is number of columns               */

void cpy_bit2(int i,char *aci,int j,char *acj,int n)
{
    register int k,ii,jj;
    register char *pi,*pj;

    ii = i * n;
    jj = j * n;
    for (k = 0; k < n; ++k) {
        pi = aci + ii / 8;
        pj = acj + jj / 8;
        *pi &= ~BMsk[ii % 8];
        if (*pj & BMsk[jj % 8])
            *pi |= BMsk[ii % 8];
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

int gcyc(void)
{
    register int i,j,ii;
    int err,n,nc,nb,ne;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Cycles in undirected graphs. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto GCYCFin;

    if (gdd_check(PMGN,0))
        goto GCYCFin;

    if (gdd_tcheck(1,1,0))      /* must be undirected and gdd option 1 */
        goto GCYCFin;
   
    if (PMOPT > 4)
        PMOPT = 4;

    if (alloc_acn(GD_NP))                 
        goto GCYCFin;
    if (alloc_acm(GD_NP))                 
        goto GCYCFin;

    printf1("\nComponent  nodes  edges  fund cycles  all cycles\n");

    GSCON_NN = 0;
    GSCON_ID = 1;
    nc = ii = 0;

    while (ii >= 0) {

        if (SILENTFlg < 2 && GD_NP > 50) {
            printfe("Component: %7d     %c",GSCON_ID,CR);
            fflushe();
        }

        GSCON_CNT = 0;
        if (GD_TYP == 1)
            u_visit(ii,PMGN);
        else
            u_visit1(ii,PMGN);

        ii = -1;
        n = 0;
        for (j = 0; j < GD_NP; ++j) {
            if (AcN[j] == GSCON_ID)  
                AcM[n++] = j;
            else if (ii < 0 && AcN[j] == 0)
                ii = j;
        }
        nc++;                   /* new component */
        ne = 0;                 /* number of edges */
        if (n > 0) {
            for (i = 0; i < n; ++i) {
                for (j = i + 1; j < n; ++j) {
                    if (gdd_adj(AcM[i],AcM[j],PMGN) >= 0.0)
                        ne++;
                }
            }
        }
        nb = ne - n + 1;        /* number of basic cycles, if not negative */
        GSCON_NP = 0;           /* number of all cycles */

        if (nb > 0) {
            nb = g_cyc(nc,ne,n,AcM,PMOPT);
            if (nb < 0)
                goto GCYCFin;
        }
        else
            nb = 0;

        printf1("%6d %9d %6d  %8d   %10d\n",nc,n,ne,nb,GSCON_NP);

        GSCON_ID++;  
    }
    if (SILENTFlg < 2 && GD_NP > 50) {
        printfe("\n");
        fflushe();
    }
    printf1("\n%d records written to: %s\n",GSCON_NN,PMFdName);
    err = 0;

GCYCFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_cyc       Fundamental set of cycles. Requires an undirected graph.    */
/*                                                                          */
/*                                                                          */
/*              Return number of basic cycles, or -1 if error.              */

int g_cyc(int nc,int ne,int n,int *nodes,int opt)
{
    register int i,j,l,k1,k2;
    int err,it,ip,in,jn,nj,nb,cptra;
    int **cptr;

    err = -1;
    cptra = 0; 
    nb = ne - n + 1;        /* number of basic cycles, if not negative */
    if (nb <= 0)
        return(0);

    if (opt >= 2) {

        if (!(cptr = (int **)calloc(nb,sizeof(int *)))) {
            p_err(-2,1);
            goto G_CYCFin;  
        }
        cptra = nb;    
        memrq(cptra,sizeof(int *));

        if (alloc_acs(nb))
            goto G_CYCFin; 
    }

    if (alloc_acr(GD_NP + 1))
        goto G_CYCFin; 
    if (alloc_acc(n))
        goto G_CYCFin; 
    if (alloc_acd(n))
        goto G_CYCFin; 
    if (alloc_aci(n))
        goto G_CYCFin; 
    if (alloc_acj(n))
        goto G_CYCFin; 
    if (alloc_ack(n))
        goto G_CYCFin; 

    for (i = 0; i < n; ++i)  
        AcR[nodes[i]] = i;

    nb = i = it = 0;
    AcD[0] = 1;
    AcI[0] = 0;
    AcJ[0] = -1;
    ip = 1;
     
    while (i >= 0) {
        in = nodes[i];
        nj = GD_FPN[in];
        for (l = 0; l < nj; ++l) {
            jn = GD_FPI[in][l];
            j = AcR[jn];

            if (AcD[j])  
                continue;

            k1 = GD_FPK[in][l];
            k2 = GD_FPK1[in][l];
            if ((k1 >= 0 && gdd_ev(k1,PMGN) >= 0.0) ||
                             (k2 >= 0 && gdd_ev(k2,PMGN) >= 0.0)) {  

                if (AcC[j]) {
                    if (g_cycprn(n,nodes,i,j,AcK,AcJ,AcI,nb,nc,opt,cptr,AcS))
                        goto G_CYCFin;
                    nb++;
                }
                else {
                    AcC[j] = 1;
                    AcK[j] = ip;
                    AcI[ip] = j;
                    AcJ[ip++] = AcK[i];
                    if (it < j)
                        it = j;
                }
            }
        }
        AcD[i] = 1;
        i = -1;
        for (j = it; j >= 0; --j) {
            if (AcC[j] != 0 && AcD[j] == 0) {
                i = j;
                it = j - 1;
                break;
            }
        }
    }
    if (opt == 2)  
        g_cycprn1(n,nodes,nb,cptr,AcS,nc);
    else if (opt >= 3)  {
        if (g_cycprn2(n,nodes,ne,nb,cptr,AcS,nc,opt))
            goto G_CYCFin;
    }
    err = nb;

G_CYCFin:
    if (cptra > 0) {
        for (i = 0; i < nb; ++i) {
            if (AcS[i] > 0) {
                free((char *)cptr[i]);
                memrq(-AcS[i],sizeof(int));
            }
        }
        free((char *)cptr);
        memrq(-cptra,sizeof(int *));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_cycprn    If opt == 1 print current cycle, if opt >= 2 save in        */
/*              cptr, cptrn.                                                */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int g_cycprn(int n,int *nodes,int i,int j,int *nptr,int *tptr,int *jptr,
    int ccnt,int nc,int opt,int **cptr,int *cptrn)
{
    register int k,ni,nj;
    int nn;

    if (alloc_acns(n))
        return(-1);    
    if (alloc_acms(n))
        return(-1);    

    AcNS[0] = i;
    AcMS[0] = j;
    ni = nj = 1;

    k = nptr[i];
    while ((k = tptr[k]) >= 0)  
        AcNS[ni++] = jptr[k];

    k = nptr[j];
    while ((k = tptr[k]) >= 0)  
        AcMS[nj++] = jptr[k];

    nn = 0;
    while (ni > 0 && nj > 0 && AcNS[--ni] == AcMS[--nj])  
        nn++; 
    if (nn == 0)  
        gerr_exit(101);
    ni += 2;
    while (nj >= 0)
       AcNS[ni++] = AcMS[nj--];

    if (opt == 1) {
        fprintf(PMFd,PMNFmtS,nc);                
        fprintf(PMFd,PMNFmtS,ccnt + 1);                
        for (k = 0; k < ni; ++k)
            fprintf(PMFd,PMNFmtS,gdd_node(nodes[AcNS[k]]));
        fprintf(PMFd,"\n");                
        GSCON_NN++;
    }
    else if (opt >= 2) {

        if (!(cptr[ccnt] = (int *)calloc(ni,sizeof(int)))) {
            p_err(-2,1);
            return(-1);     
        }
        cptrn[ccnt] = ni;
        memrq(ni,sizeof(int));
        for (k = 0; k < ni; ++k)
            cptr[ccnt][k] = AcNS[k];
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  g_cycprn1   print fundamental cycles for option 2.                      */

void g_cycprn1(int n,int *nodes,int nb,int **cptr,int *cptrn,int nc)
{
    register int i,j,l,k,k1,k2;
    int nn,u,it,in,jn;

    nn = 0;
    for (i = 0; i < n; ++i) {
        in = nodes[i];
        for (j = i + 1; j < n; ++j) {
            jn = nodes[j];
            if (gdd_adj(in,jn,PMGN) >= 0.0) {

                fprintf(PMFd,PMNFmtS,nc);                
                fprintf(PMFd,PMNFmtS,++nn);                
                fprintf(PMFd,PMNFmtS,gdd_node(nodes[i]));                
                fprintf(PMFd,PMNFmtS,gdd_node(nodes[j]));                
                
                for (l = 0; l < nb; ++l) {
                    u = 0;
                    it = cptrn[l];
                    k1 = cptr[l][it - 1];
                    for (k = 0; k < it; ++k) {
                        k2 = cptr[l][k];
                        if (i == imin(k1,k2) && j == imax(k1,k2)) {
                            u = 1;
                            break;
                        }
                        k1 = k2;
                    }
                    fprintf(PMFd,"%d ",u);
                }
                fprintf(PMFd,"\n");
                GSCON_NN++;
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  g_cycprn2   print all cycles (option 3).                                */
/*              count number of cycles in GSCON_NP.                         */

int g_cycprn2(int n,int *nodes,int ne,int nb,int **cptr,int *cptrn,int nc,int opt) 
{
    register int i,j,l,k,k1,k2,iq;
    int u,it,in,jn,ie,qmax,nq,nq0,ns;

    /* first put edge list for fund cycles into AcC */

    if (alloc_acc(nb * ne / 8 + 1))
        return(-1);     

    ie = 0;
    for (i = 0; i < n; ++i) {
        in = nodes[i];
        for (j = i + 1; j < n; ++j) {
            jn = nodes[j];
            if (gdd_adj(in,jn,PMGN) >= 0.0) {

                for (l = 0; l < nb; ++l) {
                    u = 0;
                    it = cptrn[l];
                    k1 = cptr[l][it - 1];
                    for (k = 0; k < it; ++k) {
                        k2 = cptr[l][k];
                        if (i == imin(k1,k2) && j == imax(k1,k2)) {
                            u = 1;
                            break;
                        }
                        k1 = k2;
                    }
                    put_bit2(l,ie,ne,AcC,u);
                }
                ie++;
            }
        }
    }
    if (ie != ne)
        gerr_exit(102);

    qmax = (int)pow(2.0,(double)nb) + 1;
    if (qmax < 1)
        gerr_exit(103);

    if (alloc_acd(qmax * ne / 8 + 1))
        return(-1);     
    if (alloc_acns(qmax))
        return(-1);    

    cpy_bit2(0,AcD,0,AcC,ne);
    ns = nq = 1;

    for (it = 1; it < nb; ++it) {

        if (alloc_ace(qmax / 8 + 1))
            return(-1);     

        nq0 = nq;    
        for (iq = 0; iq < nq0; ++iq) {
            u = 0;
            for (j = 0; j < ne; ++j) {
                k = 0;
                k1 = get_bit2(iq,j,ne,AcD);
                k2 = get_bit2(it,j,ne,AcC);
                if ((k1 != 0 && k2 == 0) || (k1 == 0 && k2 != 0))  
                    k = 1;
                put_bit2(nq,j,ne,AcD,k);
                if (k1 != 0 && k2 != 0)
                    u = 1;
            }
            if (u)
                put_bit(nq,AcE,1);
            else
                put_bit(nq,AcE,0);
            nq++;
        }
        for (iq = 0; iq < nq; ++iq) {
            if (get_bit(iq,AcE) == 0)
                continue;

            for (i = 0; i < iq; ++i) {
                if (get_bit(i,AcE) == 0)
                    continue;

                u = 1;
                for (j = 0; j < ne; ++j) {
                    if (get_bit2(iq,j,ne,AcD) != 0 && get_bit2(i,j,ne,AcD) == 0) {
                        u = 0;
                        break;
                    }   
                }
                if (u)  
                    put_bit(i,AcE,0);
                else {
                    u = 1;
                    for (j = 0; j < ne; ++j) {
                        if (get_bit2(i,j,ne,AcD) != 0 && get_bit2(iq,j,ne,AcD) == 0) {
                            u = 0;
                            break;
                        }   
                    }
                    if (u)  
                        put_bit(iq,AcE,0);
                }   
            }
        }
        cpy_bit2(nq,AcD,it,AcC,ne);
        for (iq = 0; iq < nq; ++iq) {
            if (get_bit(iq,AcE))
                AcNS[ns++] = iq;
        }
        AcNS[ns++] = nq;   
        nq++;
    }
    if (opt >= 3) {
        ie = 0;
        for (i = 0; i < n; ++i) {
            in = nodes[i];
            for (j = i + 1; j < n; ++j) {
                jn = nodes[j];
                if (gdd_adj(in,jn,PMGN) >= 0.0) {

                    fprintf(PMFd,PMNFmtS,nc);                
                    fprintf(PMFd,PMNFmtS,ie + 1);                
                    fprintf(PMFd,PMNFmtS,gdd_node(nodes[i]));                
                    fprintf(PMFd,PMNFmtS,gdd_node(nodes[j]));                
                    if (opt == 3) {
                        for (k = 0; k < ns; ++k) {
                            l = AcNS[k];
                            fprintf(PMFd,"%d ",(int)get_bit2(l,ie,ne,AcD));
                        }
                    }
                    fprintf(PMFd,"\n");
                    GSCON_NN++;
                    ie++;
                }
            }
        }
    }
    if (opt == 4) {

        for (k = 0; k < ns; ++k) {
            l = AcNS[k];
            fprintf(PMFd,PMNFmtS,nc);                
            fprintf(PMFd,PMNFmtS,k + 1);                

            for (i = 0; i < ne; ++i)
                fprintf(PMFd,"%d ",(int)get_bit2(l,i,ne,AcD));
            fprintf(PMFd,"\n");
            GSCON_NN++;
        }
    }
    GSCON_NP = ns;
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

int gio(void)
{
    register int i,j,l,nn;
    int err,nrec,nr,ii,iter,n,nni,ni,ne;
    double val,tmp;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Integrated ownership. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto GIOFin;

    if (gdd_check(PMGN,0))
        goto GIOFin;
    if (gdd_tcheck(1,2,2))      /* gdd option 1 required */
        goto GIOFin;            /* and directed, valued */

    if (PMOPT > 2)
        PMOPT = 2;

    printf1("Maximum of edge values: %lg\n",GD_EVMax[PMGN - 1]);
    if (GD_EVMax[PMGN - 1] > 1.0) {
        printf1("Must not exceed 1.\n");
        goto GIOFin;
    }
    if (PMFmtF == 0)                /* default print format */
        pmfmt(7,4);

    if (MxItFlg == 0 || MxIter < 1)
        MxIter = 20;
    if (PMEPSFlg == 0)
        PMEPS = 1.e-4;

    printf1("Max number of iterations: %d\n",MxIter);
    printf1("Tolerance for convergence: %lg\n",PMEPS);

    ni = 0;
    if (PMF2Def) {
        if (alloc_ack(GD_NP))                 
            goto GIOFin;

        if ((ni = g_getnd(0,AcK)) < 0)         /* get node list */
            goto GIOFin;

        if (ni > 0)
            printf1("Number of input nodes: %d\n",ni);
    }
    if (alloc_acx(GD_NP))                 
        goto GIOFin;

    if (alloc_acy(GD_NP))                 
        goto GIOFin;

    if (alloc_acn(GD_NP))                 
        goto GIOFin;

    nr = nni = nrec = 0;

    for (ii = 0; ii < GD_NP; ++ii) {

        if (ni > 0 && AcK[ii] == 0)
            continue;

        GSCON_ID = 0;
        for (i = 0; i < GD_NP; ++i) {
            AcN[i] = 0;
            AcX[i] = 0.0;
        }

        visit2(ii,1.0,1.0);

        n = 0;
        for (i = 0; i < GD_NP; ++i) {
            if (AcN[i])    
                n++;
        }
        if (n == 0) {
            nni++;
            continue;
        }

        for (iter = 0; iter < MxIter; ++iter) {

            for (i = 0; i < GD_NP; ++i) {
                if (AcX[i] > 0.0) {
    
                    nn = GD_BPN[i];
                    if (nn > 0) {
                        for (l = 0; l < nn; ++l) {
                            j = GD_BPI[i][l];          /* edge from j to i */
                            if (j != i) {
                                ne = GD_BPK[i][l];
                                if (ne >= 0 && (val = gdd_ev(ne,PMGN)) > 0.0) {  
                                    if (j == ii)
                                        AcY[i] += val;
                                    else
                                        AcY[i] += AcX[j] * val;
                                }
                            }
                        }
                    }   
                }
            }
            val = 0.0;
            for (i = 0; i < GD_NP; ++i) {
                if (AcN[i]) {
                    tmp = fabs(AcX[i] - AcY[i]);
                    if (val < tmp)
                        val = tmp;
                      
                    AcX[i] = AcY[i];
                    AcY[i] = 0.0;
                }
            }
            if (val <= PMEPS)
                break;
        }
        for (i = 0; i < GD_NP; ++i) {
            if (AcN[i]) {
                tmp = gdd_adj(ii,i,PMGN);
                if (tmp < 0.0)
                    tmp = 0.0;

                if (PMOPT == 2 && tmp < PMSC && AcX[i] < PMSC)
                    continue;

                fprintf(PMFd,PMNFmtS,ii + 1);
                fprintf(PMFd,PMNFmtS,gdd_node(ii));
                fprintf(PMFd,PMNFmtS,gdd_node(i));
                fprintf(PMFd,PMFmtS,tmp);
                fprintf(PMFd,PMFmtS,AcX[i]);
                fprintf(PMFd,PMNFmtS,n);
                fprintf(PMFd,PMNFmtS,iter);
                fprintf(PMFd,"\n");
                nrec++;
            }
        }
        nr++;
        n_info(ii + 1,10);
    }
    n_info_e();
    printf1("%d records written to: %s\n",nrec,PMFdName);
    if (nni > 0)
        printf1("skipped %d isolated nodes.\n",nni);
    err = 0;

GIOFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  visit2      called by gio(). Graph always directed, valued.             */

void visit2(int k,double val,double val1)
{
    register int j,l,nn,kk;
    double tmp;

    AcN[k] = ++GSCON_ID;
    AcX[k] += val * val1;

    nn = GD_FPN[k];
    if (nn > 0) {
        for (l = 0; l < nn; ++l) {
            j = GD_FPI[k][l];                    /* edge from k to j */
            if (AcN[j] == 0) {
                kk = GD_FPK[k][l];
                if (kk >= 0) {
                    tmp = gdd_ev(kk,PMGN);
                    if (tmp >= 0.0)
                        visit2(j,tmp,val);
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

int gfcf(void)
{
    register int i,j,k,l,nn;
    int err,nrec,nrec1,nrec2,n,nr,nb,level,ni,nii,ne,iflag;
    double sc,tmp,tmp1;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Forward control flows. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto GFCFFin;

    if (gdd_check(PMGN,0))
        goto GFCFFin;
    if (gdd_tcheck(1,2,0))      /* directed, gdd option 1 required */
        goto GFCFFin;   

    if (PMSCFlg && PMSC >= 0.0)
        sc = PMSC;
    else
        sc = 0.5;   

    printf1("Minimum level of influence: %lg\n",sc);

    if (PMOPT > 2)
        PMOPT = 2;
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    if (alloc_ack(GD_NP))                 
        goto GFCFFin;

    iflag = 0;
    if (PMF2Def) {                      /* get node list from input file */
        if ((ni = g_getnd(0,AcK)) < 0)         
            goto GFCFFin;
        printf1("Number of input nodes: %d\n",ni);
        iflag = 1;
    }
    if (alloc_acn(GD_NP))                 
        goto GFCFFin;

    if (alloc_acm(GD_NP))                 
        goto GFCFFin;

    if (alloc_acx(GD_NP))                 
        goto GFCFFin;

    if (PMTabFDef) {
        if (alloc_acr(GD_NP))                 
            goto GFCFFin;
    }
    nii = nrec = nrec1 = nrec2 = 0;
    GSCON_ID = 1;

    for (i = 0; i < GD_NP; ++i) {       /* loop throuh all requested nodes */

        if (iflag && AcK[i] == 0)
            continue;

        n_info(i + 1,1);
   
        for (j = 0; j < GD_NP; ++j) {
            AcN[j] = 0;
            AcM[j] = -1;
        }
        visit3(i);

        nr = -1;
        for (j = 0; j < GD_NP; ++j) {
            if (AcN[j])
                nr++;
        }

        AcM[i] = 0;
        level = 1;
        nb = 0;
        while (1) {
            n = 0;
            for (j = 0; j < GD_NP; ++j) {
                if (AcN[j] && AcM[j] < 0) {
                    tmp = 0.0;
                    nn = GD_BPN[j];
                    if (nn > 0) {
                        for (l = 0; l < nn; ++l) {
                            k = GD_BPI[j][l];          /* edge from k to j */
                            if (AcM[k] >= 0) {
                                ne = GD_BPK[j][l];
                                if (ne >= 0 && (tmp1 = gdd_ev(ne,PMGN)) > 0.0)    
                                    tmp += tmp1;
                            }
                        }
                    }
                    AcX[j] = tmp;
                    if (tmp > sc) {
                        AcM[j] = -2;
                        n++;
                    }
                }
            }
            if (n == 0)
                break;

            for (j = 0; j < GD_NP; ++j) {
                if (AcM[j] == -2) {
                    AcM[j] = level;
                    nb++;
                }
            }
            level++;
        }

        if (PMF1Def) {
            fprintf(PMF1d,PMNFmtS,i + 1);           
            fprintf(PMF1d,PMNFmtS,gdd_node(i));           
            fprintf(PMF1d,PMNFmtS,nb);           
        }
        if (nb > 0) {
            if (PMOPT == 2) {
                nr = 0;
                for (j = 0; j < GD_NP; ++j) {
                    if (AcM[j] > 0)
                        nr++;
                }
            }
            k = 0;
            for (j = 0; j < GD_NP; ++j) {
                if (AcN[j] && j != i) {
                    if (PMOPT == 2 && AcM[j] < 1)
                        continue;

                    fprintf(PMFd,PMNFmtS,i + 1);           
                    fprintf(PMFd,PMNFmtS,gdd_node(i));           
                    fprintf(PMFd,PMNFmtS,gdd_node(j)); 
                    fprintf(PMFd,PMNFmtS,AcM[j]);           
                    fprintf(PMFd,PMFmtS,AcX[j]);           
                    fprintf(PMFd,PMNFmtS,nb);           
                    fprintf(PMFd,PMNFmtS,nr);           
                    fprintf(PMFd,PMNFmtS,++k);           
                    fprintf(PMFd,"\n");           
                    nrec++;

                    if (AcM[j] >= 0) {
                        if (PMF1Def)  
                            fprintf(PMF1d,PMNFmtS,gdd_node(j));           

                        if (PMTabFDef)  
                            AcR[j] += 1;
                    }
                }
            }
            if (k)
                nii++;
        }
        if (PMF1Def) {
            fprintf(PMF1d,"\n");           
            nrec1++;
        }
    }
    n_info_e();

    if (PMTabFDef) {
        for (i = 0; i < GD_NP; ++i) {
            fprintf(PMTabFd,PMNFmtS,i + 1);           
            fprintf(PMTabFd,PMNFmtS,gdd_node(i));           
            fprintf(PMTabFd,PMNFmtS,AcR[i]);           
            fprintf(PMTabFd,"\n");           
            nrec2++;
        }
    }
    printf1("Number of nodes controlling at least one other node: %d\n",nii);
    printf1("%d records written to: %s\n",nrec,PMFdName);
    if (nrec1 > 0)
        printf1("%d records written to: %s\n",nrec1,PMF1dName);
    if (nrec2 > 0)
        printf1("%d records written to: %s\n",nrec2,PMTabFName);
    err = 0;

GFCFFin:       
    p_clean();
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

int gbcf(void)
{
    register int i,j,k,l,nn;
    int err,nrec,n,nb,level,ni,nni,nnic,fnd,nd,ndd,cptr,ne;   
    double sc,tmp,tmp1;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Backward control flows. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto GBCFFin;

    if (gdd_check(PMGN,0))
        goto GBCFFin;
    if (gdd_tcheck(1,2,0))      /* directed, gdd option 1 required */
        goto GBCFFin;   

    if (PMOPT > 2)
        PMOPT = 2;
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    if (PMSCFlg && PMSC >= 0.0)
        sc = PMSC;
    else
        sc = 0.5;  

    if (PMNS < 2 * GD_NP)
        PMNS = 2 * GD_NP;

    printf1("Max number of controlled nodes: %d\n",PMNS);
    printf1("Minimum level of influence: %lg\n",sc);

    ni = 0;
    if (PMF2Def) {
        if (alloc_ack(GD_NP))                 
            goto GBCFFin;
    
        if ((ni = g_getnd(0,AcK)) <= 0)         /* get node list */
            goto GBCFFin;

        if (ni > 0)
            printf1("Number of input nodes: %d\n",ni);
    }
    if (alloc_acj(PMNS))        /* list of nodes controlled by i */
        goto GBCFFin;

    if (alloc_aci(PMNS))        /* corresponding controllers */
        goto GBCFFin;

    if (alloc_acns(GD_NP))      /* number of nodes controlled by i */
        goto GBCFFin;

    if (alloc_acn(GD_NP))                 
        goto GBCFFin;

    if (alloc_acm(GD_NP))                 
        goto GBCFFin;

    /* first make list of all nodes controlled by an ultimate node */

    GSCON_ID = 1;
    nrec = nnic = nni = cptr = 0;

    for (i = 0; i < GD_NP; ++i) {       /* loop throuh all nodes */
   
        n = g_pre(i,PMGN);
        if (n > 0)                      /* use only nodes with zero in-degree */
            continue;

        nni++;
        n_info(i + 1,1);

        for (j = 0; j < GD_NP; ++j) {
            AcN[j] = 0;
            AcM[j] = -1;
        }
        visit3(i);

        AcM[i] = 0;
        level = 1;
        nb = 0;
        while (1) {
            n = 0;
            for (j = 0; j < GD_NP; ++j) {
                if (AcN[j] && AcM[j] < 0) {
                    tmp = 0.0;

                    nn = GD_BPN[j];
                    if (nn > 0) {
                        for (l = 0; l < nn; ++l) {
                            k = GD_BPI[j][l];          /* edge from k to j */
                            if (AcM[k] >= 0) {
                                ne = GD_BPK[j][l];
                                if (ne >= 0 && (tmp1 = gdd_ev(ne,PMGN)) > 0.0)    
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
                                  tmp1 = get_data(GD_EV[0],GD_JE[k]);
                                  if (tmp1 > 0)   
                                      tmp += tmp1;
                            }
                            k++;
                        }
                    }
**/

                    if (tmp > sc) {
                        AcM[j] = -2;
                        n++;
                    }
                }
            }
            if (n == 0)
                break;

            for (j = 0; j < GD_NP; ++j) {
                if (AcM[j] == -2) {
                    AcM[j] = level;
                    nb++;
                }
            }
            level++;
        }
        AcNS[i] = nb;
        if (nb > 0) {
            k = 0;
            for (j = 0; j < GD_NP; ++j) {
                if (AcN[j] && AcM[j] > 0 && j != i) {
                    if (cptr >= PMNS) {
                        printf1("Error: exceeded max number of controlled nodes (adjust ns).\n");
                        goto GBCFFin;
                    }
                    AcI[cptr] = i;
                    AcJ[cptr++] = j;
                }
            }   
        }
        if (PMF1Def) {
            fprintf(PMF1d,PMNFmtS,i + 1);        
            fprintf(PMF1d,PMNFmtS,gdd_node(i));     
            fprintf(PMF1d,PMNFmtS,nb);        
            fprintf(PMF1d,"\n");
            nrec++;
        }
        if (nb > 0)
            nnic++;
    }       
    n_info_e();

    printf1("\nNumber of ultimate nodes: %d\n",nni);
    printf1("Controlling at least one other node: %d\n",nnic);
    if (PMF1Def)  
        printf1("%d records written to: %s\n",nrec,PMF1dName);
    if (nnic == 0)  
        goto GBCFFin;

    if (alloc_acx(GD_BPMax))                 
        goto GBCFFin;

    if (alloc_acy(GD_BPMax))                 
        goto GBCFFin;

    if (alloc_acm(GD_BPMax))                 
        goto GBCFFin;

    if (alloc_acn(GD_BPMax))                 
        goto GBCFFin;

    if (alloc_acr(GD_BPMax))                 
        goto GBCFFin;

    nni = nrec = 0;
    for (i = 0; i < GD_NP; ++i) {       /* loop through all requested nodes */

        if (ni > 0 && AcK[i] == 0)
            continue;
    
        n_info(i + 1,1);

        /* find direct backward links in AcY */

        nd = 0;
        for (j = 0; j < GD_BPMax; ++j) {
            AcX[j] = AcY[j] = 0.0;
            AcN[j] = AcM[j] = -1;
        }

        nn = GD_BPN[i];
        if (nn > 0) {
            for (l = 0; l < nn; ++l) {
                k = GD_BPI[i][l];          /* edge from k to i */
                ne = GD_BPK[i][l];
                if (ne >= 0 && (tmp = gdd_ev(ne,PMGN)) > 0.0) {  
                    AcM[nd] = k;
                    AcY[nd++] = tmp;
                }
            }
        }

/**
        k = GD_JPtr[i];
        if (k >= 0) {
            for (j = 0; j < GD_JEN[i]; ++j) {
                tmp = gdd_adj(GD_JEI[k],i,PMGN);
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
            k = AcM[j];
            fnd = 0;
            for (l = 0; l < cptr; ++l) {
                if (AcJ[l] == k) {
                    AcN[j] = AcI[l];
                    fnd = 1;
                    break;
                }
            }
            if (fnd == 0)
                AcN[j] = k;
        }
        ndd = 0;
        for (j = 0; j < nd; ++j) {
            if (AcN[j] >= 0) {
                AcR[ndd] = AcN[j];
                for (k = j; k < nd; ++k) {
                    if (AcN[k] == AcR[ndd]) {
                        AcX[ndd] += AcY[k];
                        AcN[k] = -1;
                    }
                }
                ndd++;
            }
        }
        if (PMOPT == 1) {
            for (j = 0; j < nd; ++j) {

                fprintf(PMFd,PMNFmtS,i + 1);        
                fprintf(PMFd,PMNFmtS,gdd_node(i));
                fprintf(PMFd,PMNFmtS,nd);        
                fprintf(PMFd,PMNFmtS,gdd_node(AcM[j]));        
                fprintf(PMFd,PMFmtS,AcY[j]);        
                fprintf(PMFd,PMNFmtS,ndd);        
                if (j < ndd) {
                    k = gdd_node(AcR[j]);
                    tmp = AcX[j];
                }
                else {
                    k = -1;
                    tmp = 0.0;
                }
                fprintf(PMFd,PMNFmtS,k);        
                fprintf(PMFd,PMFmtS,tmp);        
                fprintf(PMFd,"\n");
                nrec++;
            }
        }
        else {
            
            if (sortd(nd,AcY,1))
                break;
            
            if (sortd(ndd,AcX,1))
                break;
            
            fprintf(PMFd,PMNFmtS,i + 1);        
            fprintf(PMFd,PMNFmtS,gdd_node(i));
            fprintf(PMFd,PMNFmtS,GD_BPMax);        
            fprintf(PMFd,PMNFmtS,nd);        
            fprintf(PMFd,PMNFmtS,ndd);        
            for (j = 0; j < GD_BPMax; ++j)
                fprintf(PMFd,PMFmtS,AcY[j]);                
            for (j = 0; j < GD_BPMax; ++j)
                fprintf(PMFd,PMFmtS,AcX[j]);                
            fprintf(PMFd,"\n");                
            nrec++;
        }
    }
    n_info_e();

    printf1("%d records written to: %s\n",nrec,PMFdName);
    if (nni > 0)
        printf1("Skipped %d input nodes having zero in-degree.\n",nni);

    err = 0;

GBCFFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  visit3      called by gfcf(), gbcf()                                    */

void visit3(int k)
{
    register int j,l,kk,nn;
   
    AcN[k] = GSCON_ID;

    nn = GD_FPN[k];
    if (nn > 0) {
        for (l = 0; l < nn; ++l) {
            j = GD_FPI[k][l];                    /* edge from k to j */
            if (AcN[j] == 0) {
                kk = GD_FPK[k][l];
                if (kk >= 0 && gdd_ev(kk,PMGN) >= 0.0)
                    visit3(j);
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

int gep(void)
{
    register int i,j,k,l,kk;
    int err,opt,opt1,ne,ii,maxp,mpi,mpj;
    double tmp;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Enumeration of paths. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto GEPFin;

    if (PMOPT > 6)
        PMOPT = 6;
    if (PMGN < 1)
        PMGN = 1;
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);
    if (PMMax < 1)
        PMMax = 100;
    if (PMNS < 1)
        PMNS = 100;

    if (gdd_check(PMGN,0))
        goto GEPFin;
    if (gdd_tcheck(1,0,0))      /* gdd option 1 required */
        goto GEPFin;

    if (alloc_acn(GD_NP + 2))                 
        goto GEPFin;
    if (alloc_acc(GD_NP))                 
        goto GEPFin;

    ne = 0;
    if (PMF2Def) {
        if (alloc_acr(PMNS))                 
            goto GEPFin;
        if (alloc_acs(PMNS))                 
            goto GEPFin;

        if ((ne = g_getne(PMNS,AcR,AcS)) < 0)         /* get edge list */
            goto GEPFin;

        printf1("Number of edges read from %s: %d\n",PMF2dName,ne);
        if (ne < 1)
            goto GEPFin;
    }
    opt1 = opt = 0;
    if (PMOPT == 2)
        opt = 1;
    else if (PMOPT == 5) {
        opt1 = 1;
        if (alloc_acm(GD_NP))                 
            goto GEPFin;
    }
    else if (PMOPT == 6) {
        printf1("Max number of paths: %d\n",PMMax);
        opt1 = 2;
        if (alloc_acm(PMMax * GD_NP))       /* for internal nodes */
            goto GEPFin;
        if (alloc_acx(PMMax))               /* for values */
            goto GEPFin;
        if (alloc_acj(PMMax))               /* for length */
            goto GEPFin;
        if (alloc_aci(PMMax))               /* pointer for sorting */
            goto GEPFin;
    }
    maxp = 0;
    GSCON_NN = 0;
    ii = i = 0;
GEPCONTI:
    if (ne > 0) {
        i = AcR[ii];
        j = AcS[ii];
    }     

    if (PMOPT == 3 || PMOPT == 4) {
        fprintf(PMFd,PMNFmtS,i + 1);
        fprintf(PMFd,PMNFmtS,gdd_node(i));    
    }
    if (ne <= 0)
        j = 0;

GEPCONTJ:
    if (PMOPT == 3 || PMOPT == 4) {
        if (i == j) {
            fprintf(PMFd,PMFmtS,0.0);
            goto GEPCONT;        
        }
    }
    else if (j == i)  
        goto GEPCONT;
         
    else if (ne <= 0 && j < i && GD_GT <= 2)            
        goto GEPCONT;
           

    for (k = 0; k < GD_NP; ++k) {
        AcC[k] = 0;
        if (opt1 == 1)
            AcM[k] = 0;
    }

    GSCON_NP = 0;           /* # of paths */
    GSCON_NSP = 0;          /* # of shortest paths */
    GSCON_LEN = -1;         /* length of actual paths */
    GSCON_MINLEN = 0;       /* length of shortest paths */
    GSCON_MAXLEN = 0;       /* length of longest paths */
    GSCON_VAL = 0.0;        /* value of actual paths */
    GSCON_MINVAL = 0.0;     /* min value of paths */
    GSCON_MAXVAL = 0.0;     /* max value of paths */
    GSCON_CNT = 0;
    visit4(i,i,j,0.0,1,opt,opt1);

    if (maxp < GSCON_NP) {
        maxp = GSCON_NP;
        mpi  = i;
        mpj  = j;
    }
    tmp = -1.0;
    if (PMOPT == 1 || opt1 == 1) {
        if (GSCON_NP > 0) {
            fprintf(PMFd,PMNFmtS,i + 1);
            fprintf(PMFd,PMNFmtS,j + 1);
            fprintf(PMFd,PMNFmtS,gdd_node(i));
            fprintf(PMFd,PMNFmtS,gdd_node(j));
            fprintf(PMFd,PMNFmtS,GSCON_NP);
            fprintf(PMFd,PMNFmtS,GSCON_NSP);
            fprintf(PMFd,PMNFmtS,GSCON_MINLEN);
            fprintf(PMFd,PMFmtS,GSCON_MINVAL);
            if (opt1 == 0) {
                fprintf(PMFd,PMNFmtS,GSCON_MAXLEN);
                fprintf(PMFd,PMFmtS,GSCON_MAXVAL);
            }   
            else {
                for (k = 0; k < GD_NP; ++k)
                    fprintf(PMFd,PMNFmtS,AcM[k]);
            }
            fprintf(PMFd,"\n");
            GSCON_NN++;
        }
    }
    else if (PMOPT == 3) {
        if (GSCON_NP > 0)
            tmp = GSCON_MINVAL;
        fprintf(PMFd,PMFmtS,tmp);
    }
    else if (PMOPT == 4) {
        if (GSCON_NP > 0)
             tmp = GSCON_MAXVAL;
        fprintf(PMFd,PMFmtS,tmp);
    }
    else if (PMOPT == 6) {
        if (GSCON_NP > PMMax) {
            printf1("Error: exceeded max number of paths (%d,%d).\n",
                                                   gdd_node(i),gdd_node(j));
            goto GEPFin;
        }
        if (sortdp(GSCON_NP,AcX,AcI))
            goto GEPFin;

        for (k = 0; k <  GSCON_NP; ++k) {
            fprintf(PMFd,PMNFmtS,i + 1);
            fprintf(PMFd,PMNFmtS,j + 1);
            fprintf(PMFd,PMNFmtS,gdd_node(i));
            fprintf(PMFd,PMNFmtS,gdd_node(j));
            fprintf(PMFd,PMNFmtS,GSCON_NP);
            fprintf(PMFd,PMNFmtS,k + 1);

            kk = AcI[k];
            fprintf(PMFd,PMNFmtS,AcJ[kk]);
            fprintf(PMFd,PMFmtS,AcX[kk]);
            kk *= GD_NP;
            for (l = 0; l < GD_NP; ++l)
                fprintf(PMFd,PMNFmtS,AcM[kk + l]);
            fprintf(PMFd,"\n");
            GSCON_NN++;
        }
    }

GEPCONT:
    if (ne > 0) {
        if (PMOPT == 3 || PMOPT == 4) {
            fprintf(PMFd,"\n");
            GSCON_NN++;
        }
        if (++ii < ne)
            goto GEPCONTI;
    }
    else {
        if (++j < GD_NP)
            goto GEPCONTJ;
        if (PMOPT == 3 || PMOPT == 4) {
            fprintf(PMFd,"\n");
            GSCON_NN++;
        }
        n_info(i + 1,1);
        if (++i < GD_NP)
            goto GEPCONTI;
    }
    n_info_e(); 
    printf1("Max number of paths (nodes: %d,%d): %d\n",
                                         gdd_node(mpi),gdd_node(mpj),maxp);
    printf1("%d records written to: %s\n",GSCON_NN,PMFdName);
    err = 0;

GEPFin:       
    p_clean();
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

void visit4(int k,int k1,int k2,double val,int first,int opt,int opt1)
{
    register int j,l,nn;
    double tmp;

    /* printf1("visit4 k=%d k1=%d k2=%d first=%d\n",k,k1,k2,first); */

    if (opt1 == 2 && GSCON_NP >= PMMax) {
        GSCON_NP = PMMax + 1;
        AcC[k] = 1;
        return;
    }
    if (first == 0 || k != k1 || k != k2)  
        AcC[k] = 1;         
    AcN[GSCON_CNT++] = k;
    GSCON_LEN++;
    GSCON_VAL += val;
       
    if (k == k2 && first == 0) {

        if (k1 == k2 && GSCON_LEN <= 2)     /* for cycles min length 3 */
            goto V4CONT;

        if (opt1 == 2) {
            AcX[GSCON_NP] = GSCON_VAL;
            AcJ[GSCON_NP] = GSCON_LEN;
            nn = GSCON_NP * GD_NP;
            for (j = 0; j < GD_NP; ++j)  
                AcM[nn + j] = 0;
            for (j = 1; j < GSCON_CNT - 1; ++j)  
                AcM[nn + AcN[j]] = 1;
        }

        GSCON_NP++;
        if (GSCON_MINVAL <= 0.0 || GSCON_MINVAL > GSCON_VAL) {
            GSCON_MINVAL = GSCON_VAL;
            GSCON_MINLEN = GSCON_LEN;
            GSCON_NSP = 1;

            if (opt1 == 1) {          
                for (j = 0; j < GD_NP; ++j)  
                    AcM[j] = 0;
                for (j = 1; j < GSCON_CNT - 1; ++j)  
                    AcM[AcN[j]] = 1;
            }      
        }
        else if (fabs(GSCON_MINVAL - GSCON_VAL) < EPSI1) {
            GSCON_NSP++;
            if (GSCON_MINLEN > GSCON_LEN)
                GSCON_MINLEN = GSCON_LEN;

            if (opt1 == 1) {          
                for (j = 1; j < GSCON_CNT - 1; ++j)  
                    AcM[AcN[j]] += 1;
            }      
        }
        if (GSCON_MAXVAL <= 0.0 || GSCON_MAXVAL < GSCON_VAL) {
            GSCON_MAXVAL = GSCON_VAL;
            GSCON_MAXLEN = GSCON_LEN;
        }
        else if (fabs(GSCON_MAXVAL - GSCON_VAL) < EPSI1) {
            if (GSCON_MAXLEN < GSCON_LEN)
                GSCON_MAXLEN = GSCON_LEN;
        }

        if (opt) {
            fprintf(PMFd,PMNFmtS,k1 + 1);
            if (opt == 1)
                fprintf(PMFd,PMNFmtS,k2 + 1);
            fprintf(PMFd,PMNFmtS,gdd_node(k1));
            if (opt == 1)
                fprintf(PMFd,PMNFmtS,gdd_node(k2));

            fprintf(PMFd,PMNFmtS,GSCON_NP);
            fprintf(PMFd,PMNFmtS,GSCON_LEN);
            fprintf(PMFd,PMFmtS,GSCON_VAL);
            for (j = 0; j < GSCON_CNT; ++j)  
                fprintf(PMFd,PMNFmtS,gdd_node(AcN[j]));
            fprintf(PMFd,"\n");
            GSCON_NN++;
        }   
        
V4CONT:
        GSCON_LEN--;      
        GSCON_VAL -= val;
        AcC[k] = 0;
        GSCON_CNT--;             
        return;
    }
    nn = GD_FPN[k];
    if (nn > 0) {  
        for (l = 0; l < nn; ++l) {
            j = GD_FPI[k][l];                    /* edge from k to j */
            if (AcC[j] == 0) {
                tmp = gdd_adj(k,j,PMGN);
                if (tmp >= 0.0)  
                    visit4(j,k1,k2,tmp,0,opt,opt1);
            }
        }
    }                  
    GSCON_LEN--;      
    GSCON_VAL -= val;
    if (GSCON_CNT <= 0) {
        printfe("ERROR in visit4.\n");
        gerr_exit(216);
    }
    AcC[k] = 0;
    GSCON_CNT--;
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

int gflow(void)
{
    register int i,j;
    int err,n;
    double f;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Maximal flows. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 5,1,1))       /* get parameters */
        goto GFLOWFin;

    if (gdd_check(PMGN,0))
        goto GFLOWFin;

    if (gdd_tcheck(1,2,0))      /* must be directed, gdd option 1 */
        goto GFLOWFin;

    if (PMOPT > 2)
        PMOPT = 2;
    if (PMFmtF == 0)            /* default print format */
        pmfmt(10,4);

    if (alloc_acn(GD_NP))       /* used as queue in g_flowbfs() */
        goto GFLOWFin;
    if (alloc_aci(GD_NP))       /* used for pre[] in g_flowbfs */ 
        goto GFLOWFin;
    if (alloc_acc(GD_NP))       /* used for vis[] in g_flowbfs */
        goto GFLOWFin;
    if (alloc_acj(GD_NP))       /* used for pointers to flow[] */
        goto GFLOWFin;
    if (alloc_acz(GD_NP))       /* used for fd in g_flowbfs */
        goto GFLOWFin;

    n = 0;
    for (i = 0; i < GD_NP; ++i) {

        if (SILENTFlg < 2 && GD_NP >= 500) {
            printfe("Node: %7d     %c",i + 1,CR);
            fflushe();
        }
        fprintf(PMFd,PMNFmtS,i + 1);
        fprintf(PMFd,PMNFmtS,gdd_node(i));

        for (j = 0; j < GD_NP; ++j) {

            if (j == i) {
                if (PMOPT == 2)
                    fprintf(PMFd,PMFmtS,0.0);
                continue;
            }

            if (alloc_acx(GD_NOC))      /* used for flow[] */
                goto GFLOWFin;

            f = g_flow(GD_NP,i,j,PMGN,AcN,AcI,AcC,AcJ,AcX,AcZ);

            if (PMOPT == 1) {
                if (f < 0.0)
                    continue;

                fprintf(PMFd,PMNFmtS,++n);
                fprintf(PMFd,PMNFmtS,i + 1);
                fprintf(PMFd,PMNFmtS,j + 1);
                fprintf(PMFd,PMNFmtS,gdd_node(i));
                fprintf(PMFd,PMNFmtS,gdd_node(j));
                fprintf(PMFd,PMFmtS,f);
                fprintf(PMFd,"\n");
            }
            else  
                fprintf(PMFd,PMFmtS,f);

        }
        if (PMOPT == 2) {
            fprintf(PMFd,"\n");
            n++;
        }
    }
    if (SILENTFlg < 2 && GD_NP >= 500) {
        printfe("\n");
        fflushe();
    }
    printf1("%d records written to: %s\n",n,PMFdName);
    err = 0;

GFLOWFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_flow.     Find maximal flow between nodes i0 and i1.                  */
/*              Uses the Edmonds and Karp algorithm as described in         */
/*              V. Turau, Algorithmische Graphentheorie, Addison-Wesley,    */  
/*              1996, pp. 171ff.                                            */
/*                                                                          */
/*  Return value of max flow, or -1 if there is no path from i0 to i1.      */

double g_flow(int n,int i0,int i1,int gn,int *q,int *pre,char *vis,int *fptr,  
    double *flow,double *fd)
{
    register int i,k,first;
    double fdd,fmax;        

    first = 1;
    fmax = 0.0;
    while (1) {
        fdd = g_flowbfs(n,i0,i1,gn,q,pre,vis,fptr,flow,fd);               
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

double g_flowbfs(int n,int i0,int i1,int gn,int *q,int *pre,char *vis,int *fptr,  
    double *flow,double *fd)
{
    register int i,j,l,k;
    int aptr,bptr,m,nf;
    double tmp;

    for (i = 0; i < n; ++i)
        vis[i] = 0;

    nf = 1;
    pre[i0] = -1;  
    vis[i0] = 1;
    fd[i0] = DBLMAX;
    fd[i1] = 0.0;
    q[0] = i0;
    aptr = 0;
    bptr = 1;

    while (aptr < bptr) {
        i = q[aptr++];
        m = GD_FPN[i];
        for (l = 0; l < m; ++l) {
            j = GD_FPI[i][l];           /* edge from i to j */
            if (vis[j] == 0) {
                k = GD_FPK[i][l];           /* pointer to edge and flow */
                tmp = gdd_ev(k,gn);   
                if (tmp >= 0.0) {
                    if (j == i1)
                        nf = 0;
                    if (flow[k] < tmp) {
                        vis[j] = 1;     
                        pre[j] = i;
                        fptr[j] = k + 1;
                        q[bptr++] = j;
                        fd[j] = dmin(fd[i],tmp - flow[k]);
                    }
                }
            }
        }
        m = GD_BPN[i];
        for (l = 0; l < m; ++l) {
            j = GD_BPI[i][l];           /* edge from j to i */
            if (vis[j] == 0) {
                k = GD_BPK[i][l];       /* pointer to edge and flow */
                tmp = gdd_ev(k,gn);   
                if (tmp >= 0.0) {
                    if (flow[k] > 0.0) {
                        vis[j] = 1;     
                        pre[j] = i;
                        fptr[j] = -(k + 1);
                        q[bptr++] = j;
                        fd[j] = dmin(fd[i],flow[k]);
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

int gfc(void)
{
    register int i,j,k;
    int err,n;
    double f,f0;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Flow control. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto GFCFin;

    if (gdd_check(PMGN,0))
        goto GFCFin;

    if (gdd_tcheck(1,2,0))      /* must be directed, gdd option 1 */
        goto GFCFin;

    if (PMFmtF == 0)            /* default print format */
        pmfmt(10,4);

    if (alloc_acn(GD_NP))       /* used as queue in g_flowbfs() */
        goto GFCFin;
    if (alloc_aci(GD_NP))       /* used for pre[] in g_flowbfs */ 
        goto GFCFin;
    if (alloc_acc(GD_NP))       /* used for vis[] in g_flowbfs */
        goto GFCFin;
    if (alloc_acj(GD_NP))       /* used for pointers to flow[] */
        goto GFCFin;
    if (alloc_acz(GD_NP))       /* used for fd in g_flowbfs */
        goto GFCFin;

    n = 0;
    for (i = 0; i < GD_NP; ++i) {

        if (SILENTFlg < 2 && GD_NP >= 500) {
            printfe("Node: %7d     %c",i + 1,CR);
            fflushe();
        }
        for (j = 0; j < GD_NP; ++j) {
            if (j == i)  
                continue;

            /* first calculate max flow from i to j */

            if (alloc_acx(GD_NOC))      /* used for flow[] */
                goto GFCFin;

            f0 = g_flow(GD_NP,i,j,PMGN,AcN,AcI,AcC,AcJ,AcX,AcZ);

            if (f0 < 0.0)
                continue;

            fprintf(PMFd,PMNFmtS,i + 1);
            fprintf(PMFd,PMNFmtS,j + 1);
            fprintf(PMFd,PMNFmtS,gdd_node(i));
            fprintf(PMFd,PMNFmtS,gdd_node(j));
            fprintf(PMFd,PMFmtS,f0);

            for (k = 0; k < GD_NP; ++k) {

                f = -1.0;
                if (k != i && k != j) {
                                    
                    if (alloc_ace(GD_NOC))  
                        goto GFCFin;
                                    
                    if (fc_visit(i,k,PMGN)) {      

                        if (alloc_ace(GD_NOC))  
                            goto GFCFin;
                                    
                        if (fc_visit(k,j,PMGN)) {      

                            if (alloc_acx(GD_NOC))      /* used for flow[] */
                                goto GFCFin;

                            f = g_cflow(GD_NP,i,j,PMGN,AcN,AcI,AcC,AcJ,AcX,AcZ,k);
                            if (f < 0.0)
                                f = 0.0;
                        }   
                    }
                }
                fprintf(PMFd,PMFmtS,f);
            }
            fprintf(PMFd,"\n");
            n++;
        }
    }
    if (SILENTFlg < 2 && GD_NP >= 500) {
        printfe("\n");
        fflushe();
    }
    printf1("%d records written to: %s\n",n,PMFdName);
    err = 0;

GFCFin:       
    p_clean();
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

double g_cflow(int n,int i0,int i1,int gn,int *q,int *pre,char *vis,int *fptr,  
    double *flow,double *fd,int k0)
{
    register int i,k,first;
    double fdd,fmax;        

    first = 1;
    fmax = 0.0;
    while (1) {
        fdd = g_cflowbfs(n,i0,i1,gn,q,pre,vis,fptr,flow,fd,k0);               
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

double g_cflowbfs(int n,int i0,int i1,int gn,int *q,int *pre,char *vis,int *fptr,  
    double *flow,double *fd,int k0)
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
    fd[i0] = DBLMAX;
    fd[i1] = 0.0;
    q[0] = i0;
    aptr = 0;
    bptr = 1;

    while (aptr < bptr) {
        i = q[aptr++];
        m = GD_FPN[i];
        for (l = 0; l < m; ++l) {
            j = GD_FPI[i][l];           /* edge from i to j */
            if (vis[j] == 0) {
                k = GD_FPK[i][l];           /* pointer to edge and flow */
                tmp = gdd_ev(k,gn);   
                if (tmp >= 0.0) {
                    if (j == i1)
                        nf = 0;
                    if (flow[k] < tmp) {
                        vis[j] = 1;     
                        pre[j] = i;
                        fptr[j] = k + 1;
                        q[bptr++] = j;
                        fd[j] = dmin(fd[i],tmp - flow[k]);
                    }
                }
            }
        }
        m = GD_BPN[i];
        for (l = 0; l < m; ++l) {
            j = GD_BPI[i][l];           /* edge from j to i */
            if (vis[j] == 0) {
                k = GD_BPK[i][l];       /* pointer to edge and flow */
                tmp = gdd_ev(k,gn);   
                if (tmp >= 0.0) {
                    if (flow[k] > 0.0) {
                        vis[j] = 1;     
                        pre[j] = i;
                        fptr[j] = -(k + 1);
                        q[bptr++] = j;
                        fd[j] = dmin(fd[i],flow[k]);
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

int fc_visit(int k,int k0,int gn)
{
    register int j,l;
    int n;

    if (k == k0)
        return(1);

    AcE[k] = 1;         

    n = GD_FPN[k];
    if (n <= 0)
        return(0);       

    for (l = 0; l < n; ++l) {

        j = GD_FPI[k][l];                    /* edge from k to j */
        if (AcE[j] == 0) {
            if (gdd_adj(k,j,gn) >= 0.0) {
                if (fc_visit(j,k0,gn))
                    return(1);
            }
        }
    }
    return(0);
}

