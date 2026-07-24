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

/*  functions in t_mes.c */

int ghd(void);
void ghd_init(void);
void visit_ghd(int n,int k,int k1,int k2,double val,int first);
int ghd_d(int m,double *xi,double *xj);

int ghd1(void);
void ghd1_fnd(int n,int l,int nm,int ni,int *num,int *idx,int *lev);
int ghd1_h(int n,int ni,int *num,int *lev);

int conj(void);
void conj_free(void);
int conj_r(void);
int conj_tab(void); 
int conj_dmat(int opt);
int conj_com(void);
int conj_it(void); 
int conj_eval(void);              
void conj_prn(void);              
void conj_res(void);              
int conj_lp(void);

int mreg(void);
int mreg_it(int nx1);

int mreg_proj(int n,double *x);
int mreg_proj1(int n,int m,int *t,int *p,int *pp,double *x,double *y);
int mreg_proj2(int n,int m,int *t,double *x,double *y);
void mreg_ms(int n,double *y,double *m,double *s);
void mreg_std(int n,int nx,double *x,double *beta,double *z,double m,double s);
int mreg_ties(int n,double *y,int *t);

int nmca(void);
int nmca_com(void);


int GHD1Max = 0;            /* max number of nodes                          */

int CJDIM = 0;              /* number of dimensions                         */
int *CJNCat;                /* number of categories in i.th dimension       */
int CJNCatA = 0;            /* if allocated                                 */
int **CJCat;                /* value of categories in the i.th dimension    */
int CJCatA = 0;             /* if allocated                                 */
int *CJCatAI;               /* if CTCat[i] is allocated                     */
int CJCatAIA = 0;
char *CJDMat;               /* design matrix                                */
int CJDMatA = 0;
int CJDMCol = 0;            /* columns in design matrix                     */
int CJNNCat = 0;            /* total number of categories                   */  

double *CJR;                /* original rank order                          */
int CJRA = 0; 
double *CJRY;               /* estimated rank order                         */
int CJRYA = 0; 
double *CJEY;               /* estimated dep variable                       */
int CJEYA = 0; 
double CJRMean = 0.0;       /* mean of original rank order                  */
double CJRSD   = 0.0;       /* stand. deviation of original rank order      */
int NTies = 0;              /* number of tie blocks                         */


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

int ghd(void)
{
    register int i,j,k;
    int err,n,m,a,nn,nn1,ne,n2,na,nb,nna,nnb,la,lb,lmin,lmax;
    double x,xa,xb;

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("Hasse diagram. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 3,4,1))       /* get parameters */
        goto GHDFin;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(8,4);

    n = NOC;
    n2 = n + 2;                         /* plus two additional nodes */
    m = PMNV;                           /* number of variables */

    if (alloc_acx(n2 * m + 1))          /* input data */
        goto GHDFin;
    if (alloc_ack(n + 1))               /* node list */
        goto GHDFin;

    for (i = 0; i < n; ++i) {           /* get data into AcX */
        for (j = 0; j < m; ++j)  
            AcX[i * m + j] = get_data(PMVIdx[j],i);
    }

    /* create list of nodes with zero indegree */

    for (i = 0; i < n; ++i) {
        a = 1;
        for (k = 0; k < n; ++k) {
            if (k != i && ghd_d(m,AcX + k * m,AcX + i * m)) {
                a = 0;
                break;
            }
        }
        AcK[i] = a;
    }

    printf1("\nNodes with zero indegree:");
    for (i = 0; i < n; ++i) {
        if (AcK[i])
            printf1(" %d",i + 1);
    }
    newline();          

    if (PMPCFDef) {
        fprintf(PMPCFd,"digraph g {\n");
        for (i = 1; i <= n; ++i) 
            fprintf(PMPCFd,"n%d [label=%d];\n",i,i);
        nn1 = n + 1;
    }
    ne = nn = 0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (j == i || ghd_d(m,AcX + i * m,AcX + j * m) == 0)   
                continue;
            a = 1;
            for (k = 0; k < n; ++k) {
                if (k == i || k == j)
                    continue;
                if (ghd_d(m,AcX + i * m,AcX + k * m) && 
                    ghd_d(m,AcX + k * m,AcX + j * m)) { 
                    a = 0;
                    break;
                }
            }
            if (a) {
                if (PMF1Def) {
                    fprintf(PMF1d,"%5d %5d\n",i + 1,j + 1);
                    nn++;
                }
                if (PMPCFDef) {
                    fprintf(PMPCFd,"n%d -> n%d;\n",i+1,j+1);
                    nn1++;
                }
                ne++;
            }
        }
    }
    printf1("Hasse diagram: %d nodes, %d edges.\n",n,ne);

    if (PMF1Def)  
        printf1("%d records written to: %s\n",nn,PMF1dName);
    if (PMPCFDef) {
        fprintf(PMPCFd,"}\n");
        nn1++;
        printf1("%d records written to: %s\n",nn1,PMPCFName);
    }
    newline();

    if (PMTabFDef == 0) {
        err = 0;
        goto GHDFin;
    }
    printf1("Creating interval-valued ranks.\n");

    /***
    for (i = 0; i < n; ++i) {   
        for (j = 0; j < m; ++j)  
            printf("%g ",AcX[i * m + j]);
        newline();
    }
    ***/
           
    /*  create minimal and maximal elements in x[n] and x[n+1] */

    na = n;
    nb = n + 1;
                   
    for (j = 0; j < m; ++j) {
        xa = xb = AcX[j];
        for (i = 1; i < n; ++i) {
            x = AcX[i * m + j];
            xa = dmin(xa,x);  
            xb = dmax(xb,x);
        }
        AcX[na * m + j] = xa;
        AcX[nb * m + j] = xb;
    }
    nna = nnb = -1;

    for (i = 0; i < n; ++i) {
        a = 1;
        for (j = 0; j < m; ++j) {
            if (AcX[i * m + j] != AcX[na * m + j]) {
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
            if (AcX[i * m + j] != AcX[nb * m + j]) {
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
    printf1("Minimal element: %g",AcX[na * m]); 
    for (j = 1; j < m; ++j)
        printf1(",%g",AcX[na * m + j]);
    if (nna < 0) {
        nna = na;
        nn++;
        printf1(" (not contained, new node number: %d)\n",nna + 1);
    }
    else           
        printf1(" (contained, node number: %d)\n",nna + 1);

    printf1("Maximal element: %g",AcX[nb * m]); 
    for (j = 1; j < m; ++j)
        printf1(",%g",AcX[nb * m + j]);
    if (nnb < 0) {
        nnb = nn;
        if (nnb < nb) {
            for (j = 0; j < m; ++j)
                AcX[nnb * m + j] = AcX[nb * m + j];
        }
        nn++;
        printf1(" (not contained, new node number: %d)\n",nnb + 1);
    }
    else           
        printf1(" (contained, node number: %d)\n",nnb + 1);

    /* create adjacency matrix in AcD[] */

    if (alloc_acd(nn * nn + 1))     
        goto GHDFin;

    ne = 0;
    for (i = 0; i < nn; ++i) {
        for (j = 0; j < nn; ++j) {
            if (j == i || ghd_d(m,AcX + i * m,AcX + j * m) == 0)   
                continue;
            a = 1;
            for (k = 0; k < nn; ++k) {
                if (k == i || k == j)
                    continue;
                if (ghd_d(m,AcX + i * m,AcX + k * m) && 
                    ghd_d(m,AcX + k * m,AcX + j * m)) { 
                    a = 0;
                    break;
                }
            }
            if (a) {
                ne++;
                AcD[i * nn + j] = 1;
            }
        }
    }
    printf1("New Hasse diagram: %d nodes, %d edges.\n",nn,ne);

    /***
    for (i = 0; i < nn; ++i) {
        for (j = 0; j < nn; ++j)
            printf("%d ",AcD[i * nn + j]);
        newline(); 
    }
    ***/

    /* find length of longest path from nna to nnb */

    if (alloc_acc(nn + 1))     
        goto GHDFin;

    for (k = 0; k < nn; ++k)
        AcC[k] = 0;

    ghd_init();
    visit_ghd(nn,nna,nna,nnb,0.0,1);
    lmin = GSCON_MINLEN;
    lmax = GSCON_MAXLEN;

    /***
    printf("np=%d nsp=%d len=%d minlen=%d maxlen=%d val=%g  minval=%g maxval=%g cnt=%d\n",
        GSCON_NP,GSCON_NSP,GSCON_LEN,GSCON_MINLEN,GSCON_MAXLEN,
        GSCON_VAL,GSCON_MINVAL,GSCON_MAXVAL,GSCON_CNT);
    ***/
        
    printf1("\nNumber of paths from node %d to %d: %d\n",nna+1,nnb+1,GSCON_NP);
    printf1("Minimal length: %d, maximal length: %d\n",lmin,lmax);

    /** printf1("\nNode  MinLenght  MaxLength  Interval of ranks\n"); **/

    nn1 = 0;
    for (k = 0; k < nn; ++k) {
        if (k == nna || k == nnb)
            continue;

        ghd_init();
        visit_ghd(nn,nna,nna,k,0.0,1);
        la = GSCON_MINLEN;

        ghd_init();
        visit_ghd(nn,k,k,nnb,0.0,1);
        lb = GSCON_MINLEN;

        /** printf1("%4d  %9d  %9d  [%7d,%7d]\n",k+1,la,lb,la,lmax - lb); **/

        fprintf(PMTabFd,"%4d  %9d  %9d  %8d %8d\n",k+1,la,lb,la,lmax - lb);
        nn1++;
    }
    printf1("%d records written to: %s\n\n",nn1,PMTabFName);
                        
    err = 0;

GHDFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ghd_init()  Initialize values for visit_ghd().                          */

void ghd_init(void)
{
    GSCON_NP = 0;           /* # of paths */
    GSCON_NSP = 0;          /* # of shortest paths */
    GSCON_LEN = -1;         /* length of actual paths */
    GSCON_MINLEN = 0;       /* length of shortest paths */
    GSCON_MAXLEN = 0;       /* length of longest paths */
    GSCON_VAL = 0.0;        /* value of actual paths */
    GSCON_MINVAL = 0.0;     /* min value of paths */
    GSCON_MAXVAL = 0.0;     /* max value of paths */
    GSCON_CNT = 0;
}

/* ------------------------------------------------------------------------ */
/*  visit_ghd   called by ghd.                                              */

void visit_ghd(int n,int k,int k1,int k2,double val,int first)
{
    register int j;
    double tmp;

    /* printf1("visit_ghd k=%d k1=%d k2=%d first=%d\n",k,k1,k2,first); */

    if (first == 0 || k != k1 || k != k2)  
        AcC[k] = 1;         
    GSCON_CNT++;
    GSCON_LEN++;
    GSCON_VAL += val;
       
    if (k == k2 && first == 0) {

        if (k1 == k2 && GSCON_LEN <= 2)     /* for cycles min length 3 */
            goto V4CONT;

        GSCON_NP++;
        if (GSCON_MINVAL <= 0.0 || GSCON_MINVAL > GSCON_VAL) {
            GSCON_MINVAL = GSCON_VAL;
            GSCON_MINLEN = GSCON_LEN;
            GSCON_NSP = 1;
        }
        else if (fabs(GSCON_MINVAL - GSCON_VAL) < EPSI1) {
            GSCON_NSP++;
            if (GSCON_MINLEN > GSCON_LEN)
                GSCON_MINLEN = GSCON_LEN;
        }
        if (GSCON_MAXVAL <= 0.0 || GSCON_MAXVAL < GSCON_VAL) {
            GSCON_MAXVAL = GSCON_VAL;
            GSCON_MAXLEN = GSCON_LEN;
        }
        else if (fabs(GSCON_MAXVAL - GSCON_VAL) < EPSI1) {
            if (GSCON_MAXLEN < GSCON_LEN)
                GSCON_MAXLEN = GSCON_LEN;
        }
        
V4CONT:
        GSCON_LEN--;      
        GSCON_VAL -= val;
        AcC[k] = 0;
        GSCON_CNT--;             
        return;
    }
    for (j = 0; j < n; ++j) {
        if (AcD[k * n + j]) {          /* edge from k to j */
            if (AcC[j] == 0) {
                tmp = 1.0;                
                if (tmp >= 0.0)  
                    visit_ghd(n,j,k1,k2,tmp,0);
            }
        }
    }

/*******
    n0 = GD_FPN[k];
    if (n0 > 0) {  
        for (l = 0; l < n0; ++l) {
            j = GD_FPI[k][l];               
            if (AcC[j] == 0) {
                tmp = gdd_adj(k,j,PMGN);
                if (tmp >= 0.0)  
                    visit_ghd(j,k1,k2,tmp,0);
            }
        }
    }                  
**/


    GSCON_LEN--;      
    GSCON_VAL -= val;
    if (GSCON_CNT <= 0) {
        printfe("ERROR in visit_ghd.\n");
        gerr_exit(216);
    }
    AcC[k] = 0;
    GSCON_CNT--;
}

/* ------------------------------------------------------------------------ */
/*  ghd_d(n,m,i,j,x)                                                        */
/*                                                                          */
/*  Return: 1 if xi <= xj, otherwise 0.                                     */

int ghd_d(int m,double *xi,double *xj)
{
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

int ghd1(void)
{
    register int i,j,k;
    int err,n,m,a,ne,nm;

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("Find maximal hierarchy. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,4,1))       /* get parameters */
        goto GHD1Fin;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(8,4);

    n = NOC;
    m = PMNV;                           /* number of variables */

    if (alloc_acx(n * m + 1))           /* input data */
        goto GHD1Fin;
    if (alloc_acd(n * n + 1))           /* adjacency matrix */
        goto GHD1Fin;
    if (alloc_ack(n + 1))               /* flags node numbers */
        goto GHD1Fin;
    if (alloc_acn(n + 1))               /* node numbers */
        goto GHD1Fin;
    if (alloc_acs(n + 1))               /* level */
        goto GHD1Fin;
    if (alloc_acm(n + 1))               /* final selection */
        goto GHD1Fin;

    for (i = 0; i < n; ++i) {           /* get data into AcX */
        for (j = 0; j < m; ++j)  
            AcX[i * m + j] = get_data(PMVIdx[j],i);
    }

    ne = 0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (j == i || ghd_d(m,AcX + i * m,AcX + j * m) == 0)   
                continue;
            a = 1;
            for (k = 0; k < n; ++k) {
                if (k == i || k == j)
                    continue;
                if (ghd_d(m,AcX + i * m,AcX + k * m) && 
                    ghd_d(m,AcX + k * m,AcX + j * m)) { 
                    a = 0;
                    break;
                }
            }
            if (a) {
                AcD[i * n + j] = 1;
                ne++;
            }
        }
    }
    printf1("Hasse diagram: %d nodes, %d edges.\n",n,ne);

    GHD1Max = 1;
    nm = 0;
    for (k = 0; k < n; ++k) {
        AcN[0] = k;
        AcK[k] = 1;
          
        ghd1_fnd(n,1,nm,1,AcN,AcK,AcS);

        printf1("\nGHD1Max=%d : ",GHD1Max);
        for (i = 0; i < GHD1Max; ++i)
            printf1("%d ",AcM[i]);
        newline();

        break;

        AcK[k] = -1;
        nm++;
    }




    err = 0;

GHD1Fin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ghd1_fnd()                                                              */
/*                                                                          */

void ghd1_fnd(int n,int l,int nm,int ni,int *num,int *idx,int *lev)
{
    register int i,j;
    int l1;

    if (n - nm <= GHD1Max)  
        return;

    /* add one more node */

    l1 = l + 1;
                
    for (i = 0; i < n; ++i) {
        if (idx[i])       
            continue;

        num[ni] = i;
        idx[i] = l1;
        if (ghd1_h(n,ni + 1,num,lev) == 0) {
            idx[i] = -l1;
            nm++;
            continue;
        }
        if (ni + 1 > GHD1Max) {
            printf1("GHD1Max=%d ni+1=%d\n",GHD1Max,ni + 1);         
            GHD1Max = ni + 1;
            for (j = 0; j <= ni; ++j) {
                AcM[j] = num[j];
                printf1("%d ",AcM[j]);
            }
            newline();
        }
        ghd1_fnd(n,l1,nm,ni + 1,num,idx,lev);
    }
    for (i = 0; i < n; ++i) {
        if (iabs(idx[i]) == l1)
            idx[i] = 0;
    }
}

/* ------------------------------------------------------------------------ */
/*  ghd1_h(n,ni,num,lev)    Return 1 if the nodes in num[] can be           */
/*                          hierarchically partitioned, otherwise 0.        */

int ghd1_h(int n,int ni,int *num,int *lev)
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
            if (AcD[in * n + jn]) {
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
                    if (AcD[in * n + jn]) {
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
/*  conj        Conjoint analysis (alternating least squares).              */
/*                                                                          */
/*              conj(                                                       */
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

int conj(void)
{
    int err;

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("Conjoint analysis. Current memory: %d bytes.\n",MemReq);

    TOLF = 1.e-6;
    TOLSP = 1.e-6;
        
    if (parm(CmdBuf + 4,4,1))     /* get parameters */
        goto CONJFin;

    if (PMALG != 2)
        PMALG = 1;

    if (PMOPT > 3)
        PMOPT = 1;

    if (PMFmtF == 0)            /* default print format */
        pmfmt(8,4);

    if (PMNFmtF == 0) {
        PMNFmt = 2;
        makenfmt(&PMNFmt,PMNFmtS,1,SEPC);
    }
    if (MxIter < 0 || MxItFlg == 0)
        MxIter = 50;

    if (PMNV < 3) {
        printf1("Error: need at least two variables on right-hand side.\n");
        goto CONJFin;
    }
    CJDIM = PMNV - 1;               /* number of dimensions */

    printf1("\nRank order variable: %s\n",VName[PMVIdx[0]]);

    if (conj_r())                   /* check rank order */
        goto CONJFin;

    if (conj_tab())                 /* check dimensions */
        goto CONJFin;

    /* check for ties */

    if (alloc_acm(NOC + 1))            
        goto CONJFin;

    printf1("Mean of rank order variable: %g\n",CJRMean);
    printf1("Stand. dev. of rank order variable: %g\n\n",CJRSD);

    if (PMALG == 1) {
        printf1("Algorithm: alternating least squares.\n\n");

        if (conj_dmat(0))                   /* create design matrix */
            goto CONJFin;

        NTies = mreg_ties(NOC,CJR,AcM);     /* AcM reserved to record ties */
        printf1("Number of tie blocks: %d ",NTies);
        if (NTies == NOC)
            printf1("(no ties).");
        printf1("\n\n");

        if (conj_it())  
            goto CONJFin;
    }
    else if (PMALG == 2) {
        printf1("Algorithm: linear programming.\n\n");

        if (conj_dmat(1))               /* create design matrix */
            goto CONJFin;

        NTies = mreg_ties(NOC,CJR,AcM);     /* AcM reserved to record ties */
        printf1("Number of tie blocks: %d ",NTies);
        if (NTies == NOC)
            printf1("(no ties).");
        printf1("\n\n");


        if (conj_lp())                  /* try to solve the problem */
            goto CONJFin;
    }
    err = 0;

CONJFin:       
    conj_free();
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  conj_free. Free allocated memory.                                       */

void conj_free(void)
{
    register int i;

    if (CJNCatA > 0) {
        free((char *)CJNCat);
        memrq(-CJNCatA,sizeof(int));
        CJNCatA = 0;    
    }
    if (CJCatA > 0) {
        for (i = 0; i < CJDIM; ++i) {
            if (CJCatAI[i] > 0) {
                free((char *)CJCat[i]);
                memrq(-CJCatAI[i],sizeof(int));
            }
        }
        free((char *)CJCat);
        memrq(-CJCatA,sizeof(int *));
        CJCatA = 0;    
    }
    if (CJCatAIA > 0) {
        free((char *)CJCatAI);
        memrq(-CJCatAIA,sizeof(int));
        CJCatAIA = 0;    
    }
    if (CJDMatA > 0) {
        free((char *)CJDMat);
        memrq(-CJDMatA,sizeof(char));
        CJDMatA = 0;    
    }
    if (CJRA > 0) {
        free((char *)CJR);
        memrq(-CJRA,sizeof(double));
        CJRA = 0;    
    }
    if (CJRYA > 0) {
        free((char *)CJRY);
        memrq(-CJRYA,sizeof(double));
        CJRYA = 0;    
    }
    if (CJEYA > 0) {
        free((char *)CJEY);
        memrq(-CJEYA,sizeof(double));
        CJEYA = 0;    
    }
}

/* ------------------------------------------------------------------------ */
/*  conj_r      check rank order and calculate mean and std. dev.           */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int conj_r(void)
{
    register int i;
    double tmp;

    if (!(CJR = (double *) calloc(NOC + 1,sizeof(double)))) {   
        p_err(-2,1);
        return(-1);   
    }
    CJRA = NOC + 1;
    memrq(CJRA,sizeof(double));

    if (!(CJRY = (double *) calloc(NOC + 1,sizeof(double)))) {   
        p_err(-2,1);
        return(-1);   
    }
    CJRYA = NOC + 1;
    memrq(CJRYA,sizeof(double));

    if (!(CJEY = (double *) calloc(NOC + 1,sizeof(double)))) {   
        p_err(-2,1);
        return(-1);   
    }
    CJEYA = NOC + 1;
    memrq(CJEYA,sizeof(double));

    for (i = 0; i < NOC; ++i)  
        CJR[i] = get_data(PMVIdx[0],i);

    CJRMean = 0.0;
    for (i = 0; i < NOC; ++i)  
        CJRMean += get_data(PMVIdx[0],i);
    CJRMean /= (double)NOC;

    CJRSD = 0.0;
    for (i = 0; i < NOC; ++i) {
        tmp = get_data(PMVIdx[0],i) - CJRMean;
        CJRSD += tmp * tmp;
    }
    CJRSD /= (double)NOC;
    if (CJRSD > 0.0)
        CJRSD = sqrt(CJRSD);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  conj_tab    Create frequency distritution.                              */
/*              Return: 0 if OK, -1 if error.                               */

int conj_tab(void)
{
    register int i,j,k;
    int err,n;
    double tmp;

    err = -1;

    /*  Calculation of categories of dimensions of the table:
        CJNCat[i] = number of categories in i.th dimension,
        CJCat[i][j] (j = 1,2,3,...,CJNCat[i]) is the value of the j.th
        categorie in the i.th dimension. */

    if (!(CJNCat = (int *) calloc(CJDIM,sizeof(int)))) {   
        p_err(-2,1);
        goto CONJTFin;
    }
    CJNCatA = CJDIM;
    memrq(CJNCatA,sizeof(int));

    if (!(CJCatAI = (int *) calloc(CJDIM,sizeof(int)))) {   
        p_err(-2,1);
        goto CONJTFin;
    }
    CJCatAIA = CJDIM;
    memrq(CJCatAIA,sizeof(int));

    if (!(CJCat = (int **) calloc(CJDIM,sizeof(int *)))) {   
        p_err(-2,1);
        goto CONJTFin;
    }
    CJCatA = CJDIM;
    memrq(CJCatA,sizeof(int *));

    if (alloc_acr(NOC + 1))             /* data */
        goto CONJTFin;
    if (alloc_acs(NOC + 1))             /* values of categories */
        goto CONJTFin;
    if (alloc_aci(NOC + 1))             /* frequencies */
        goto CONJTFin;
    if (alloc_acj(NOC + 1))             /* bf pointer */
        goto CONJTFin;

    printf1("Dimension  Variable ");
    prnchar(' ',VNameLen - 8,0);
    printf1("  Categories\n");

    CJNNCat = 0;
    for (i = 1; i <= CJDIM; ++i) {

        for (j = 0; j < NOC; ++j) {
            AcR[j + 1] = (int)get_data(PMVIdx[i],j);
        }
        n = cfreq(1,NOC,AcR,NOC,AcS,AcI,AcJ,0,&tmp,&tmp);

        if (n < 1) {
            printf1("Error: can't sort dimension %d.\n",i);
            goto CONJTFin;
        }
        CJNCat[i - 1] = n;
        CJNNCat += n;

        if (!(CJCat[i - 1] = (int *) calloc(n + 1,sizeof(int)))) {   
            p_err(-2,1);
            goto CONJTFin;
        }
        CJCatAI[i - 1] = n + 1;
        memrq(n + 1,sizeof(int));

        printf1("%6d      %s ",i,VName[PMVIdx[i]]);
        prnchar(' ',VNameLen - strlen(VName[PMVIdx[i]]),0);
        printf1("  %3d : ",n);
        for (j = 1; j <= n; ++j) {
            k = AcS[AcJ[j]];
            CJCat[i - 1][j] = k;
            printf1(PMNFmtS,k);
        }
        newline();
    }
    newline();

    err = 0;

CONJTFin:
    alloc_acr(0);
    alloc_acs(0);
    alloc_aci(0);
    alloc_acj(0);
    return(err);
}  

/* -##--------------------------------------------------------------------- */
/*  conj_dmat   Make design matrix in CJDMat. Sort according to ranks.      */
/*              AcN contains categories corresp to design matrix.           */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int conj_dmat(int opt)
{
    register int i,j,k;
    int err,n,s,scat,first;

    err = -1;

    if (alloc_aci(NOC + 1))             /* temp. sort pointer */
        goto CJDMFin;

    if (alloc_acj(NOC + 1))             /* inverse pointer */
        goto CJDMFin;

    if (sortdp(NOC,CJR,AcI))
        goto CJDMFin;

    for (i = 0; i < NOC; ++i) {
        CJRY[i] = CJR[AcI[i]];
        AcJ[AcI[i]] = i;
    }
    printf1("Indicator variables for design matrix:\n\n");

    printf1("Idx  Dimension  Variable  ");
    prnchar(' ',VNameLen - 8,0);
    printf1("  Category\n");

    CJDMCol = 0;            /* number of colums in design matrix */

    for (j = 0; j < CJDIM; ++j) {
        for (k = 1; k <= CJNCat[j]; ++k) {

            if (k == 1 && opt == 0)
                continue;

            scat = CJCat[j][k];

            for (i = 0; i < NOC; ++i) {
                s = (int)get_data(PMVIdx[j + 1],i);
                if (s == scat) { 
                    CJDMCol++;
                    printf1("%3d  %6d      %s ",CJDMCol,j + 1,VName[PMVIdx[j + 1]]);
                    prnchar(' ',VNameLen - strlen(VName[PMVIdx[j + 1]]),0);
                    printf1("%7d\n",s);
                    break;  
                }
            }
        }
    }
    newline();                    

    if (opt == 0)
        CJDMCol++;      /* include intercept */

    if (!(CJDMat = (char *) calloc(NOC * CJDMCol + 1,sizeof(char)))) {   
        p_err(-2,1);
        goto CJDMFin;
    }
    CJDMatA = NOC * CJDMCol + 1;
    memrq(CJDMatA,sizeof(char));

    /* create values of design matrix in CJDMat. Also create:
       AcR : reference to dimensions
       AcS : reference to categories */

    if (alloc_acr(CJDMCol + 1))     
        goto CJDMFin;

    if (alloc_acs(CJDMCol + 1)) 
        goto CJDMFin;

    if (alloc_acn(NOC * CJDIM + 1))     /* categories */
        goto CJDMFin;

    n = 0;
    if (opt == 0) {
        n = 1;
        for (i = 0; i < NOC; ++i)  
            CJDMat[i * CJDMCol + 1] = 1;
    }      
    for (j = 0; j < CJDIM; ++j) {
        for (k = 1; k <= CJNCat[j]; ++k) {

            scat = CJCat[j][k];
            first = 1;
            for (i = 0; i < NOC; ++i) {
                s = (int)get_data(PMVIdx[j + 1],i);
                if (s == scat) { 

                    AcN[AcJ[i] * CJDIM + j + 1] = s;

                    if (k == 1 && opt == 0)
                        continue;

                    if (first) {
                        n++;
                        AcR[n] = j + 1;
                        AcS[n] = scat;
                        first = 0;
                    }
                    CJDMat[AcJ[i] * CJDMCol + n] = 1;
                }
            }
        }
    }
    for (i = 0; i < NOC; ++i)
        CJR[i] = CJRY[i];

    if (PMTabFDef) {            /* print design matrix */ 

        for (i = 0; i < NOC; ++i) {
            for (j = 1; j <= CJDMCol; ++j) 
                fprintf(PMTabFd,PMNFmtS,(int)CJDMat[i * CJDMCol + j]);
            fprintf(PMTabFd,"  %g\n",CJR[i]);
        }
        printf1("Design matrix written to: %s\n\n",PMTabFName);
    }
    err = 0;

CJDMFin:
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  conj_it.   Iterative procedure, alternating least squares.              */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int conj_it(void)
{
    register int i,j,k;
    int err,iter,ranka,nc,nl,r;
    double tmp,rnorml,sa,sb,mz,sz;

    err = -1;

    printf1("Handling ties option: %d\n",PMOPT);
    printf1("Max number of iterations: %d\n",MxIter);
    printf1("Tolerance for stress: %g\n",TOLF);
    printf1("Tolerance for parameter change: %g\n",TOLSP);
    printf1("Number of parameters: %d\n\n",CJDMCol);

    if (alloc_acx(CJDMCol + 1))             /* coefficients */
        goto CJITFin;

    /* ##  try to get starting values */

    if (get_dsv(CJDMCol,AcX,0,AcX,AcX,1))
        goto CJITFin;

           
    for (i = 1; i <= CJDMCol; ++i)
        printf("i=%d AcX=%f\n",i,AcX[i]);
               

    nc = CJDMCol + 1;

    if (alloc_acy(CJDMCol + 1))             /* previous coefficients */
        goto CJITFin;

    if (alloc_acw(NOC * nc + 1))             /* lsei matrix */
        goto CJITFin;

    if (PMOPT != 1) {
        if (alloc_actmp(NOC + 1))   
            goto CJITFin;

        if (PMOPT == 2) {
            if (alloc_aci(NOC + 1))   
                goto CJITFin;
            if (alloc_acj(NOC + 1))   
                goto CJITFin;
        }
    }


    if (SILENTFlg < 2) {
        fflushe();
        printfe("\n  Iter    Function Value        Par Change    Norm of Residuals\n");
        fflushe();
    }

    for (iter = 1; iter <= MxIter; ++iter) {

        if (iter == 1 && DSVFlg)
            goto CJITNxt;

        /* create lsei matrix */

        i = k = 0;
        while (i < NOC) {
            for (j = 1; j <= CJDMCol; ++j)
                AcW[i * nc + j] = (double)CJDMat[k * CJDMCol + j];
            AcW[i * nc + nc] = CJRY[k];
            i++;
            k++;
        }
        /*** 
        printf1("W\n");
        for (i = 0; i < NOC; ++i) {
            for (j = 1; j <= nc; ++j)
                printf1("%g ",AcW[i * nc + j]);
            newline();
        }
        ***/

        err = lsei(AcW,0,NOC,0,CJDMCol,AcX,0,&tmp,&rnorml,&ranka,&nl);
        if (err) {
            if (err == -4)  
                p_err(-2,1);
            else  
                printf1("LSEI error return: %d\n",err);
            goto CJITFin;  
        }

        /* standardize estimated dep. variable in AcZ and AcW */
           
CJITNxt:
              
        for (i = 1; i <= CJDMCol; ++i)
            printf("vorher i=%d AcX=%f\n",i,AcX[i]);
                

        for (i = 0; i < NOC; ++i) {
            tmp = 0.0;
            for (j = 1; j <= CJDMCol; ++j) {
                if (CJDMat[i * CJDMCol + j])
                    tmp += AcX[j];
            }   
            CJRY[i] = tmp;
        }
        mreg_ms(NOC,CJRY,&mz,&sz);

                
        printf("Mean of z = %f\n",mz);
        printf("STD  of z = %f\n",sz);
                 
 
        if (sz > 0.0) {
            for (j = 1; j <= CJDMCol; ++j)
                AcX[j] *= CJRSD / sz;
        }
        mz = 0.0;
        for (i = 0; i < NOC; ++i) {
            tmp = 0.0;
            for (j = 1; j <= CJDMCol; ++j) {
                if (CJDMat[i * CJDMCol + j])
                    tmp += AcX[j];
            }   
            CJRY[i] = tmp;
            mz += tmp;
        }
        mz /= (double)NOC; 
             
        printf("New Mean of CJRY = %f\n",mz);
              
        AcX[1] = AcX[1] + CJRMean - mz;       /* adjust intercept */

        for (i = 0; i < NOC; ++i) {
            tmp = 0.0;
            for (j = 1; j <= CJDMCol; ++j) {
                if (CJDMat[i * CJDMCol + j])
                    tmp += AcX[j];
            }   
            CJRY[i] = tmp;
        }
        mreg_ms(NOC,CJRY,&mz,&sz);

             
        printf("Mean of CJRY = %f\n",mz);
        printf("STD  of CJRY = %f\n",sz);
              
        for (i = 0; i < NOC; ++i)  
            CJEY[i] = CJRY[i];        

        /* calculate projection, depending on PMOPT */
             
        if (PMOPT == 1)
            mreg_proj(NOC,CJRY);
          
        else if (PMOPT == 2)
            mreg_proj1(NOC,NTies,AcM,AcI,AcJ,CJRY,AcTmp);

        else if (PMOPT == 3)
            mreg_proj2(NOC,NTies,AcM,CJRY,AcTmp);
        
             
        printf("\nNew AcZ/CJRY after projection\n");
        for (i = 0; i < NOC; ++i) 
            printf("i=%3d cjry=%f cjey=%f\n",i,CJRY[i],CJEY[i]);
              

        /* calculate stress */

        tmp = 0.0;
        for (i = 0; i < NOC; ++i) {
            tmp += CJEY[i];
        }
        tmp /= (double)NOC;

        sa = sb = 0.0;
        for (i = 0; i < NOC; ++i) {
            sa += (CJRY[i] - CJEY[i]) * (CJRY[i] - CJEY[i]);
            sb += (CJEY[i] - tmp) * (CJEY[i] - tmp);
        }
        if (sb > 0.0)
            sa /= sb;
        if (sa > 0.0)
            sa = sqrt(sa);

        /* calculate par change */

        sb = 0.0;
        for (j = 1; j <= CJDMCol; ++j) {
            tmp = fabs(AcX[j] - AcY[j]) / dmax(fabs(AcX[j]),1.0);
            sb = dmax(sb,tmp);
            AcY[j] = AcX[j];
        }
        if (SILENTFlg < 2) {
            printfe("  %3d  %20.13e %17.10e %17.10e\n",iter,sa,sb,rnorml);
            fflushe();
        }
        if (sa <= TOLF || sb <= TOLSP || iter >= MxIter)
            break;
    }
    printf1("Performed %d iterations.\n",iter);
    printf1("Final stress: %g\n",sa);
    printf1("Final parameter change: %g\n",sb);
    printf1("Final norm of least squares residuals: %g\n",rnorml);
    printf1("Rank of design matrix: %d\n\n",ranka);

    conj_prn(); /* print coefficients */

    r = conj_eval();              
    printf1("Type of equivalence: %d\n",r);
       
    if (PMResFDef)  
        conj_res();

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

int conj_eval(void)               
{
    register int i,j;
    int r;
    double tol = sqrt(EPSI1);

    r = 1;
    for (i = 1; i < NOC; ++i) {
        if (fabs(CJR[i] - CJR[i - 1]) <= tol) {
            if (fabs(CJEY[i] - CJEY[i - 1]) > tol) {
                r = 2;
                goto CJENxt;
            }
        }
        else if (CJEY[i] < CJEY[i - 1] + tol) {
            r = 2;
            goto CJENxt;
        }
    }

CJENxt:
    if (r == 1)
        return(1);

    for (i = 1; i < NOC; ++i) {
        if (CJR[i] > CJR[i - 1] + tol) {
            if (CJEY[i] < CJEY[i - 1] - tol)  
                return(0);
        }
    }
    for (i = 1; i < NOC; ++i) {
        for (j = 0; j < i; ++j) {
            if (CJEY[j] < CJEY[i] - tol) {
                if (CJR[j] > CJR[i] + tol)  
                    return(0);
            }
            else if (CJEY[j] > CJEY[i] + tol) {
                if (CJR[j] < CJR[i] - tol)  
                    return(0);
            }
        }
    }
    return(2);
}

/* ------------------------------------------------------------------------ */
/*  conj_prn.  Print coefficients.                                          */
/*                                                                          */

void conj_prn(void)               
{
    register int j,k;

    printf1("Idx  Dimension  Variable  ");
    prnchar(' ',VNameLen - 8,0);
    printf1("  Category  Coefficient\n");

    for (j = 1; j <= CJDMCol; ++j) {
        k = AcR[j];
        if (j == 1) {
            printf1("%3d       -      Const ",j - 1);
            prnchar(' ',VNameLen - 5,0);
            printf1("       -    ");
        }
        else {
            printf1("%3d  %6d      %s ",j - 1,k,VName[PMVIdx[k]]);
            prnchar(' ',VNameLen - strlen(VName[PMVIdx[k]]),0);
            printf1("  %6d    ",AcS[j]);
        }
        printf1(PMTFmtS,AcX[j]);
        newline();
    }
    newline();                    
}

/* ------------------------------------------------------------------------ */
/*  conj_res.  Print estimation results                                     */
/*                                                                          */

void conj_res(void)               
{
    register int i,j;

    for (i = 0; i < NOC; ++i) {
        fprintf(PMResFd,"%4d ",i + 1);
        for (j = 1; j <= CJDIM; ++j) {
            fprintf(PMResFd,PMNFmtS,AcN[i * CJDIM + j]);
        }
        fprintf(PMResFd,PMFmtS,CJR[i]);
        fprintf(PMResFd,PMFmtS,CJRY[i]);
        fprintf(PMResFd,PMFmtS,CJEY[i]);
        fprintf(PMResFd,"\n");
    }
    printf1("Estimation results written to: %s\n\n",PMResFName);
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

int mreg(void)
{
    register int i,j,k;
    int err,nv,nx,nx1;            

    err = -1;
    if (check_cmd(0))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    printf1("Monotone regression. Current memory: %d bytes.\n",MemReq);

    TOLF = 1.e-6;
    TOLSP = 1.e-6;

    if (parm(CmdBuf + 4,4,1))     /* get parameters */
        goto MRFin;

    if (PMNW != 1) {
        printf1("Error: mreg command requires cross-section data (nw=1).\n");
        goto MRFin;
    }
    if (MxIter < 0 || MxItFlg == 0)
        MxIter = 50;

    newline();
    nv = check_nvar(0);         /* check variables */
    if (nv == 0)            
        goto MRFin;        
       
    prn_nwvar(0);               /* print variables */
    nx = nx1 = nv - 1;
    nx1++;
    if (nx1 == 0)
        goto MRFin;

    /* allocate memory */

    if (alloc_acx(nx1 * NOC + 1))       /* x matrix */
        goto MRFin;
    if (alloc_acy(NOC + 1))             /* y vector */
        goto MRFin;
    if (alloc_ack(NOC + 1))             /* sort ptr for y */
        goto MRFin;


    /* get data for dependent variable and sort */

    for (i = 0; i < NOC; ++i)  
        AcY[i] = get_data(PMVIdx[0],i);

    if (sortdp(NOC,AcY,AcK)) {
        p_err(-2,1);
        goto MRFin;
    }

    /* get data in sorted order */

    for (i = 0; i < NOC; ++i) {
        AcY[i] = get_data(PMVIdx[0],AcK[i]);
        k = 1;
        AcX[i * nx1 + 1] = 1.0;
        k = 2;
                   
        for (j = 1; j < PMNV; ++j) {
            AcX[i * nx1 + k] = get_data(PMVIdx[j],AcK[i]);
            k++;
        }
    }

    /* get mean and std.dev. of AcY in CJRMean and CJRSD */

    mreg_ms(NOC,AcY,&CJRMean,&CJRSD);

    printf1("Mean of dep. variable: %g\n",CJRMean);
    printf1("Stand. dev. of dep. variable: %g\n",CJRSD);

    /* check for ties */

    if (alloc_acn(NOC + 1))            
        goto MRFin;

    NTies = mreg_ties(NOC,AcY,AcN);
            
    printf1("Number of tie blocks: %d ",NTies);
    if (NTies == NOC)
        printf1("(no ties).");
    printf1("\n\n");

    if (mreg_it(nx1))  
        goto MRFin;
    err = 0;

MRFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  mreg_it.   Iterative procedure. AcX data matrix, AcY dep. variable.     */
/*             Both sorted according to AcY.                                */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int mreg_it(int nx1)
{
    register int i,j,k;
    int err,iter,ranka,n,nc;
    double tmp,tmp1,rnorml,sa,sb,mz,sz;

    err = -1;

    printf1("Handling ties option: %d\n",PMOPT);
    printf1("Max number of iterations: %d\n",MxIter);
    printf1("Tolerance for stress: %g\n",TOLF);
    printf1("Tolerance for parameter change: %g\n",TOLSP);
    printf1("Number of parameters: %d\n\n",nx1);

    nc = nx1 + 1;

    if (alloc_acz(NOC + 1))                 /* updates of AcY */
        goto MRITFin;

    if (alloc_acu(nx1 + 1))                 /* coefficients */
        goto MRITFin;

    if (alloc_acv(nx1 + 1))                 /* previous coefficients */
        goto MRITFin;

    if (alloc_acw(NOC * nc + 1))            /* lsei matrix */
        goto MRITFin;

    if (PMOPT != 1) {
        if (alloc_actmp(NOC + 1))            
            goto MRITFin;
       
        if (PMOPT == 2) {
            if (alloc_aci(NOC + 1))            
                goto MRITFin;
            if (alloc_acj(NOC + 1))            
                goto MRITFin;
        }
    }

    if (SILENTFlg < 2) {
        fflushe();
        printfe("\n  Iter         Stress           Par Change    Norm of Residuals\n");
        fflushe();
    }
    printf1("  Iter         Stress           Par Change    Norm of Residuals\n");

    for (i = 0; i < NOC; ++i)
        AcZ[i] = AcY[i];

    for (iter = 1; iter <= MxIter; ++iter) {

        /* create lsei matrix */

newline();
for (i = 0; i < NOC; ++i)
printf("VREG   i=%3d  AcZ=%f\n",i,AcZ[i]);


        for (i = 0; i < NOC; ++i) {
            for (j = 1; j <= nx1; ++j)  
                AcW[i * nc + j] = AcX[i * nx1 + j];
            AcW[i * nc + nc] = AcZ[i];
        }
        err = lsei(AcW,0,NOC,0,nx1,AcU,0,&tmp,&rnorml,&ranka,&n);
        if (err) {
            if (err == -4)  
                p_err(-2,1);
            else  
                printf1("LSEI error return: %d\n",err);
            goto MRITFin;  
        }
                    
        /* standardize estimated dep. variable in AcZ and AcW */

        mreg_std(NOC,nx1,AcX,AcU,AcZ,CJRMean,CJRSD); 

        for (i = 0; i < NOC; ++i)  
            AcW[i] = AcZ[i];        


        /* calculate projection, depending on PMOPT */

        if (PMOPT == 1)
            mreg_proj(NOC,AcZ);

        else if (PMOPT == 2)
            mreg_proj1(NOC,NTies,AcN,AcI,AcJ,AcZ,AcTmp);

        else if (PMOPT == 3)
            mreg_proj2(NOC,NTies,AcN,AcZ,AcTmp);



        printf("\nNew AcZ after projection\n");
        for (i = 0; i < NOC; ++i) 
            printf("i=%3d acz=%f acw=%f\n",i,AcZ[i],AcW[i]);

        /* calculate stress */
      
        tmp = 0.0;
        for (i = 0; i < NOC; ++i) {
            tmp += AcW[i];
        }
        tmp /= (double)NOC;

        sa = sb = 0.0;
        for (i = 0; i < NOC; ++i) {
printf("Stress i=%d AcZ=%f AcW=%f\n",i,AcZ[i],AcW[i]);

            sa += (AcZ[i] - AcW[i]) * (AcZ[i] - AcW[i]);
            sb += (AcW[i] - tmp) * (AcW[i] - tmp);
        }
        if (sb > 0.0)
            sa /= sb;
        if (sa > 0.0)
            sa = sqrt(sa);
printf("Stress sa=%f\n\n",sa);

        /* calculate par change */
       
        sb = 0.0;
        for (j = 1; j <= nx1; ++j) {
            tmp = fabs(AcU[j] - AcV[j]) / dmax(fabs(AcU[j]),1.0);
            sb = dmax(sb,tmp);
            AcV[j] = AcU[j];
        }
        if (SILENTFlg < 2) {
            printfe("  %3d  %20.13e %17.10e %17.10e\n",iter,sa,sb,rnorml);
            fflushe();
        }
        printf1("  %3d  %20.13e %17.10e %17.10e\n",iter,sa,sb,rnorml);

        if (sa <= TOLF || sb <= TOLSP || iter >= MxIter)
            break;
       
    }
    newline();
    printf1("Performed %d iterations.\n",iter);
    printf1("Final stress: %g\n",sa);
    printf1("Final parameter change: %g\n",sb);
    printf1("Final norm of least squares residuals: %g\n",rnorml);
    printf1("Rank of design matrix: %d\n\n",ranka);
      
    prn1_coeff(nx1,AcU,AcV,0,0,PMVIdx,0);   /* print parameters */

    if (PMF1Def) {          /* write data etc to output file */
        for (i = 0; i < NOC; ++i) {
            fprintf(PMF1d,"%4d ",i + 1);
            for (j = 1; j <= nx1; ++j)  
                fprintf(PMF1d,PMFmtS,AcX[i * nx1 + j]);
            fprintf(PMF1d,PMFmtS,AcY[i]);
            fprintf(PMF1d,PMFmtS,AcW[i]);
            fprintf(PMF1d,PMFmtS,AcZ[i]);
            fprintf(PMF1d,"\n");
        }
        printf1("\nData and results written to: %s\n\n",PMF1dName);
    }
    err = 0;

MRITFin:
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  mreg_proj(n,x)    Return monotone projection of x in x.                 */
/*                    Returns always 0.                                     */

int mreg_proj(int n, double *x)
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
                if (r <= rr + EPSI1) {
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
            if (x[i] < x[i - 1] - EPSI1) {
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

int mreg_proj1(int n,int m,int *t,int *p,int *pp,double *x,double *y) 
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
            if (sortdp(t[j],x + i,pp))
                return(-1);
            for (k = 0; k < t[j]; ++k)  
                p[i + k] = pp[k] + i;
        }
        i += t[j];
    }
    for (i = 0; i < n; ++i)
        y[i] = x[p[i]];

    mreg_proj(n,y);

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

int mreg_proj2(int n,int m,int *t,double *x,double *y) 
{
    register int i,j,k,l;
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
    mreg_proj(m,y);

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

void mreg_ms(int n,double *y,double *m,double *s)
{
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

void mreg_std(int n,int nx,double *x,double *beta,double *z,double m,double s)
{
    register int i,j; 
    double tmp,mz,sz;

    for (i = 1; i <= nx; ++i)
        printf("vorher i=%d beta=%f\n",i,beta[i]);


    for (i = 0; i < n; ++i) {
        tmp = 0.0;
        for (j = 1; j <= nx; ++j)  
            tmp += x[i * nx + j] * beta[j];
        z[i] = tmp;
    }
    mreg_ms(n,z,&mz,&sz);

printf("Mean of z = %f\n",mz);
printf("STD  of z = %f\n",sz);
 
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

printf("New Mean of z = %f\n",mz);

    beta[1] = beta[1] + m - mz;

    for (i = 0; i < n; ++i) {
        tmp = 0.0;
        for (j = 1; j <= nx; ++j)  
            tmp += x[i * nx + j] * beta[j];
        z[i] = tmp;
    }
    mreg_ms(n,z,&mz,&sz);

printf("Mean of z = %f\n",mz);
printf("STD  of z = %f\n",sz);
 
}

/* ------------------------------------------------------------------------ */
/*  mreg_ties(n,y,t)                                                        */
/*                                                                          */
/*  Check n-vector y for ties. Return number of tie blocks. Sizes of        */
/*  tie block in t[].                                                       */

int mreg_ties(int n,double *y,int *t)
{
    register int i; 
    int m;                

    m = 0;
    t[m] = 1;
                      
    for (i = 1; i < n; ++i) {
        if (fabs(y[i] - y[i - 1]) <= EPSI1)  
            t[m] += 1;            
        else {
            m++;
            t[m] = 1;
        }
    }
    m++;

    /***
    printf("\nNumber of tie blocks m=%d\n",m);
    for (i = 0; i < m; ++i)
        printf("i=%3d t[i]=%d\n",i,t[i]);
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

int nmca(void)
{
    int err,r;

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("Non-metric conjoint analysis. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,4,1))     /* get parameters */
        goto NMCAFin;

    if (PMFmtF == 0)            /* default print format */
        pmfmt(8,4);

    if (PMNFmtF == 0) {
        PMNFmt = 2;
        makenfmt(&PMNFmt,PMNFmtS,1,SEPC);
    }
    if (PMNV < 3) {
        printf1("Error: need at least two variables on right-hand side.\n");
        goto NMCAFin;
    }
    CJDIM = PMNV - 1;               /* number of dimensions */

    printf1("\nRank order variable: %s\n",VName[PMVIdx[0]]);

    if (conj_r())                   /* check rank order */
        goto NMCAFin;

    if (conj_tab())                 /* check dimensions */
        goto NMCAFin;

    r = nmca_com();  
    if (r < 0)
        goto NMCAFin;

    printf1("\nFound %d partial orders.\n",r);
    if (PMF1Def)  
        printf1("%d records written to: %s\n",r * NOC,PMF1dName);

    err = 0;

NMCAFin:       
    conj_free();
    p_clean();
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  nmca_com().   Combinatorial procedure.                                  */
/*                                                                          */
/*  Return -1 if error, otherwise the number of partial orders that have    */
/*  been found.                                                             */

int nmca_com(void)
{
    register int i,j,k,ii;
    int err,n,rn,r,is,iis,np;

    err = -1;
    np = 0;

    /*  CJNCat[i] = number of categories in i.th dimension,
        CJCat[i][j] (j = 1,2,3,...,CJNCat[i]) is the value of the j.th
        categorie in the i.th dimension. */

    n = CJNNCat;        /* total number of categories */

    if (alloc_acy(NOC + 1))             /* rank order variable  */
        goto CONJCFin;
    if (alloc_acr(CJDIM * NOC + 1))     /* indep. variables  */
        goto CONJCFin;
    if (alloc_acs(NOC + 1))             /* temp. sort pointer */
        goto CONJCFin;

    for (i = 0; i < NOC; ++i) {
        AcY[i] = get_data(PMVIdx[0],i);
        for (j = 0; j < CJDIM; ++j) {
            k = (int)get_data(PMVIdx[j + 1],i) - 1;
            if (k < 0 || k >= CJNCat[j]) {
                printf1("Error in record %d category %d (dim %d, cat %d).\n",
                                i + 1,k,j + 1,CJNCat[i]);
                return(-1);
            }
            AcR[i * CJDIM + j] = k;
        }
    }
    if (sortdp(NOC,AcY,AcS))
        goto CONJCFin;

    /***
    for (i = 0; i < NOC; ++i) {
        is = AcS[i];
        printf("%f ",AcY[is]);
        for (j = 0; j < CJDIM; ++j)
            printf("%2d ",AcR[is * CJDIM + j]);
        newline();
    }
    ***/

    if (alloc_acn(CJNNCat + 1))     /* record permutations */
        goto CONJCFin;
    if (alloc_aci(CJNNCat + 1))     
        goto CONJCFin;
    if (alloc_acj(CJNNCat + 1))     
        goto CONJCFin;
    if (alloc_ack(CJDIM + 1))       /* offsets */
        goto CONJCFin;

    k = 0;
    for (i = 0; i < CJDIM; ++i) {
        AcK[i] = k;
        k += CJNCat[i];
    }

    /***
    printf("offsets: ");
    for (i = 0; i < CJDIM; ++i)
        printf("%d ",AcK[i]);
    newline();
    ***/

    /* generate first permutation */

    for (i = 0; i < CJDIM; ++i) {
        k = AcK[i];
        perm(CJNCat[i],AcN + k,AcI + k,AcJ + k,1);
        AcR[i] = 0;
    }

    while (1) {

        /*********************
        newline();
        for (i = 0; i < n; ++i)
            printf("%d ",AcN[i]);
        newline();
        *********************/

        r = 1;
        for (i = 0; i < NOC; ++i) {
            is = AcS[i];
            for (ii = i + 1; ii < NOC; ++ii) {
                iis = AcS[ii];

                /******************************
                printf("    compare: ");
                for (j = 0; j < CJDIM; ++j)
                    printf("%2d ",AcR[is * CJDIM + j]);
                printf(" and: ");
                for (j = 0; j < CJDIM; ++j)
                    printf("%2d ",AcR[iis * CJDIM + j]);
                newline();
                printf("NEW compare: ");
                for (j = 0; j < CJDIM; ++j) {
                    k = AcK[j];
                    printf("%2d ",AcN[k + AcR[is * CJDIM + j]]);

                }
                printf(" and: ");
                for (j = 0; j < CJDIM; ++j) {
                    k = AcK[j];
                    printf("%2d ",AcN[k + AcR[iis * CJDIM + j]]);
                }
                newline();
                *******************/

                rn = 0;  
                for (j = 0; j < CJDIM; ++j) {
                    k = AcK[j];
                    if (AcN[k + AcR[is * CJDIM + j]] >=
                        AcN[k + AcR[iis * CJDIM + j]]) {
                        rn++;   
                    }
                }
                if (rn == CJDIM) {
                    r = 0;
                    break; 
                }
            }
            if (r == 0)
                break;
        }
        if (r) {        /* found valid permutation */
            np++;
            printf1("Found: "); 
            for (i = 0; i < CJDIM; ++i) {
                printf1("(");
                k = AcK[i];
                ii = CJNCat[i];
                for (j = 0; j < ii; ++j) {
                    printf1(PMNFmtS,AcN[k + j] + 1);
                    if (j == ii - 1)
                        printf1(") ");
                    else
                        printf1(",");
                }
            }
            newline();

            if (PMF1Def) {      /* write to output file */

                for (i = 0; i < NOC; ++i) {
                    is = AcS[i];
                    fprintf(PMF1d,PMFmtS,AcY[is]);
                    for (j = 0; j < CJDIM; ++j) {
                        k = AcK[j];
                        fprintf(PMF1d,PMNFmtS,AcN[k + AcR[is * CJDIM + j]] + 1);
                    }
                    fprintf(PMF1d,"\n");
                }
            }
        }

        /* generate next permutation */

        for (i = CJDIM - 1; i >= 0; --i) {
            k = AcK[i];
            r = perm(CJNCat[i],AcN + k,AcI + k,AcJ + k,0);
            if (r == 1)
                break;

            if (i == 0)
                goto CONJCFin1;

            perm(CJNCat[i],AcN + k,AcI + k,AcJ + k,1);
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

int conj_lp(void)
{
    register int i,j,k;
    int err,n,m,m1,ai,bi,p;
    double z,tol;

    err = -1;
    tol = 1.0e-8;

    n = NOC;
    m = CJDMCol;
    m1 = m + 1;

    if (alloc_actmp(n * m1 + 1))       
        goto CONJLPFin;

    if (alloc_acx(m + 1))       
        goto CONJLPFin;

    if (alloc_acy(n + 1))       
        goto CONJLPFin;
      
    for (j = 1; j <= m; ++j)
        AcTmp[j] = -1.0;
      
    p = 0;
    k = 0;
    for (i = 1; i < n; ++i) {
        if (fabs(CJR[i] - CJR[i - 1]) >= EPSI1) {       /* not tied */
            k++;
            for (j = 1; j <= m; ++j) {
                ai = CJDMat[(i - 1) * m + j];
                bi = CJDMat[i * m + j];
                AcTmp[k * m1 + j] = (double)(ai - bi);
            }
            AcTmp[k * m1 + m1] = -1.0;
        }
    }
    for (i = 1; i < n; ++i) {                       
        if (fabs(CJR[i] - CJR[i - 1]) < EPSI1) {       /* tied */
            p++;
            k++;
            for (j = 1; j <= m; ++j) {
                ai = CJDMat[(i - 1) * m + j];
                bi = CJDMat[i * m + j];
                AcTmp[k * m1 + j] = (double)(ai - bi);
            }
            AcTmp[k * m1 + m1] = 0.0;
        }
    }
    /***********
    printf("actmp p=%d  \n",p );
    for (i = 0; i < n; ++i) {
        for (j = 1; j <= m1; ++j)
            printf("%4.1f ",AcTmp[i * m1 + j]);
        newline();
    }
    **********/

    err = lpf1(n - 1,m,p,AcTmp,AcX,AcY,&z,tol);
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
            if (CJDMat[i * m + j])
                z += AcX[j];   
        }
        CJRY[i] = CJEY[i] = z;
    }

    /* check for valid solution (strict monotone equivalence) */

    if (conj_eval() != 1) {       
        printf1("Cannot find a solution.\n");
        err = 1;
        goto CONJLPFin;
    }

    printf1("Idx  Dimension  Variable  ");
    prnchar(' ',VNameLen - 8,0);
    printf1("  Category    Solution\n");

    for (j = 1; j <= CJDMCol; ++j) {
        k = AcR[j];
        printf1("%3d  %6d      %s ",j - 1,k,VName[PMVIdx[k]]);
        prnchar(' ',VNameLen - strlen(VName[PMVIdx[k]]),0);
        printf1("  %6d    ",AcS[j]);
        printf1(PMTFmtS,AcX[j]);
        newline();
    }
    printf1("\nType of equivalence: 1\n");
    if (PMResFDef)                   
        conj_res();

CONJLPFin:
    if (err == 1)
        printf1("Cannot find a solution.\n");
    else if (err == -1)
        printf1("Error: insufficient memory.\n");
    return(err);
}

