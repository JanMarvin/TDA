/****************************************************************************/
/*  t_gdd                                                                   */
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
#include "t_alloc.h"
#include "t_gdat.h"
#include "t_gf.h"
#include "t_sort.h"
#include "t_rand.h"
#include "t_mat.h"

/*  functions in t_gdd.c */

int gdd(void);
void gdd_info(void);
int gdd_palloc(int n,int m);
int gdd_alloc(int n); 
void gdd_free(int opt);
void gdd_setgt(int typ,int gt,int gtt);
int gdd_setnd(void);
int gdd_setptr(int opt);
int gdd_setap(void);
int gdd_setne(int opt);
void gdd_prot(void);
int gdd_fndni(int k);
int check_gdj(int j);
double gdd_adj(int i,int j,int gn);
int g_suc(int i,int gn);
int g_pre(int i,int gn);
int g_loop(int i,int gn);
int g_getnd(int n,int *nd);
int g_getne(int n,int *ni,int *nj);
int gdd_check(int gn,int opt);
int gdd_check2(int gn,int gn1,int opt);
int gdd_tcheck(int typ,int dir,int val);
int gdd_ei(int k);             
int gdd_ej(int k);             
double gdd_ev(int k,int gn);      
int gdd_node(int i);             
void n_info(int i,int n);
void n_info_e(void);
int gcd(void);


int GD_TYP = 0;         /* type of data structure (= opt in gdd)            */
int GD_GT = 0;          /* type of graph                                    */
int GD_GTT = 0;         /* handling of bidirectional edges                  */
int GD_NG = 0;          /* number of graphs                                 */
int GD_NP = 0;          /* number of nodes                                  */
int GD_NDMAX = 0;       /* max node number                                  */
int *GD_ND;             /* list of node numbers                             */
int GD_NDA = 0;         /* allocated in GD_ND                               */
int GD_EI = -1;         /* index of variable for starting node              */
int GD_EJ = -1;         /* index of variable for ending node                */
int *GD_NE;             /* number of edges                                  */
int GD_NEA = 0;         /* allocated for GD_NE                              */
int *GD_IN;             /* number of isolated nodes                         */
int GD_INA = 0;         /* allocated for GD_IN                              */
int *GD_NL;             /* number of loops                                  */
int GD_NLA = 0;         /* allocated for GD_NL                              */
int *GD_EV;             /* list of variable numbers for edge values         */
int GD_EVA = 0;         /* allocated for GD_EV                              */
double GD_EVMin = 0.0;  /* min value of valid edges                         */
double *GD_EVMax;       /* max value of edges                               */
int GD_EVMaxA = 0;      /* allocated for GD_EVMax                           */
int *GD_APtr;           /* pointer to adjacency matrix (option 2)           */
int GD_APtrA = 0;       /* allocated for GD_APtr                            */
int *GD_FPN;            /* number of forward links                          */
int GD_FPNA = 0;        /* allocated GD_FPN                                 */
int **GD_FPI;           /* forward pointer: internal node numbers           */
int GD_FPIA = 0;
int **GD_FPK;           /* forward pointer: data matrix rows                */
int GD_FPKA = 0;
int **GD_FPK1;          /* backward pointer: undirected graph               */
int GD_FPK1A = 0;

int *GD_BPN;            /* number of backward links                         */
int GD_BPNA = 0;        /* allocated GD_BPN                                 */
int **GD_BPI;           /* backward pointer: internal node numbers          */
int GD_BPIA = 0;
int **GD_BPK;           /* backward pointer: data matrix rows               */
int GD_BPKA = 0;
int GD_FPMax = 0;       /* max value of GD_FPN[]                            */   
int GD_BPMax = 0;       /* max value of GD_BPN[]                            */   

int GD_PERM = 0;        /* set if gdd1 command is used                      */
int GD_NOC = 0;         /* number of cases for graph data                   */
int *GD_EIP;            /* values of EI                                     */
int GD_EIPA = 0;        /* allocated                                        */
int *GD_EJP;            /* values of EJ                                     */
int GD_EJPA = 0;        /* allocated                                        */
float *GD_EP;           /* edge values                                      */
int GD_EPA = 0;         /* allocated                                        */

/* ------------------------------------------------------------------------ */
/*  gdd()   Setup graph data.                                               */
/*                                                                          */  
/*          gdd(                                                            */
/*              opt=...,    1   edge list                                   */
/*                          2   pointer to adj. matrix                      */
/*                          3   lower triangle                              */
/*                          4   upper triangle                              */
/*                          5   lower triangle plus diagonal                */
/*                          6   upper triangle plus diagonal                */
/*                          7   full adj. matrix (variables)                */
/*                          8   full adj. matrix (matrices)                 */
/*              gt=...,     graph type                                      */
/*              perm=...,   perm=1 to make permanent, def. perm=0           */
/*              sc=...,     min value for valid edges, def. 0.0             */
/*                          (must be >= 0.0)                                */
/*          ) = varlist.                                                    */
/*                                                                          */
/*  opt 1: gdd = I,J,V1,...,    gt default 4                                */
/*  opt 2: gdd = I,J,V1,...,    gt default 4                                */
/*  opt 3: gdd = V1,...,        gt default 2                                */
/*  opt 4: gdd = V1,...,        gt default 2                                */
/*  opt 5: gdd = V1,...,        gt default 2                                */
/*  opt 6: gdd = V1,...,        gt default 2                                */
/*  opt 7: gdd = V1,...,        gt default 4                                */
/*  opt 8: gdd = A1,...,        gt default 4                                */
/*                                                                          */
/*  gt = 1  undirected, unvalued                                            */
/*  gt = 2  undirected, valued                                              */
/*  gt = 3  directed, unvalued                                              */
/*  gt = 4  directed, valued                                                */
/*                                                                          */
/*  If gt = 2, then gt(n) = 2 possible.                                     */
/*                                                                          */
/*  n = 1   minimum of bidirectional edges                                  */
/*  n = 2   maximum of bidirectional edges                                  */
/*  n = 3   sum of bidirectional edges                                      */
/*                                                                          */
/*  If perm = 1 save basic data in                                          */
/*  GD_NOC = number of cases                                                */
/*  GD_EIP = values of EI                                                   */
/*  GD_EJP = values of EJ                                                   */
/*  GD_EP  = edge values, dimension is GD_NOC rows, GD_NG columns           */
/*  and set GD_PERM = 1.                                                    */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int gdd(void)
{
    register int i,j,k;
    int err,n;     
    double tmp;

    err = -1;

    if (check_cmd(1))
        return(-1);

    if (!strcmp(CmdBuf,"gdd")) {
        gdd_info();
        return(0);
    }

    if (!strcmp(CmdBuf,"gdd=off")) {
        if (GD_TYP) {
            gdd_free(0);
            printf1("Removed graph data structure. Current memory: %d bytes.\n",MemReq);
        }
        err = 0;
        goto GDDFin;
    }
    gdd_free(0);         /* free previously defined graph data */

    printf1("New graph data. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 3,14,1))     /* get parameters */
        goto GDDFin;

    if (PMOPT < 1 || PMOPT > 8) {
        printf1("Error: possible options are: 1 - 8.\n");
        goto GDDFin;
    }
    if (PMOPT != 8 && PMNVTyp == 1) {
        printf1("Error: matrix names only with option 8.\n");
        goto GDDFin;
    }
    GD_EVMin = PMSC;
    if (GD_EVMin < 0.0) {
        printf1("Error: minimal value for valid edges must not be negative.\n");
        goto GDDFin;
    }

    gdd_setgt(PMOPT,PMGT,PMGTT);

    printf1("Option %d: ",PMOPT);

    if (PMOPT <= 2) {               /* setup node list: GD_ND */

        if (PMOPT == 1)
            printf1("edge list.\n");
        else
            printf1("pointer to adjacency matrix.\n");

        if (PMNV < 3) {
            printf1("Error: this option requires at least three variables.\n");
            goto GDDFin;
        }
        GD_EI = PMVIdx[0];          /* index of variable for starting node */
        GD_EJ = PMVIdx[1];          /* index of variable for ending node */
        GD_NG = PMNV - 2;           /* number of graphs */

        if (gdd_alloc(GD_NG)) 
            goto GDDFin;

        for (i = 0; i < GD_NG; ++i)
            GD_EV[i] = PMVIdx[i + 2];

        if (gdd_setnd())            /* setup node list */
            goto GDDFin;

        if (PMOPT == 1) {           /* setup pointer */
            n = 0;
            if (GD_GT <= 2)
                n = 1;
            if (gdd_setptr(n))     
                goto GDDFin;
        }
        else {
            if (gdd_setap())     
                goto GDDFin;
        }
        GD_NOC = NOC;
    }
    else if (PMOPT <= 7) {

        if (PMOPT == 3 || PMOPT == 4) {           /* lower/upper triangle */
            if (PMOPT == 3)
                printf1("lower triangle.\n");
            else
                printf1("upper triangle.\n");

            n = (int)((sqrt(8.0 * (double)NOC + 1.0) + 1.0) / 2.0);
            if (n < 1 || n * (n - 1) / 2 != NOC)  
                err = -2; 
        }
        else if (PMOPT == 5 || PMOPT == 6) {    /* lower/upper triangle plus diagonal */
            if (PMOPT == 5)
                printf1("lower");
            else
                printf1("upper");
            printf1(" triangle plus diagonal.\n");

            n = (int)((sqrt(8.0 * (double)NOC - 1.0) + 1.0) / 2.0);
            if (n < 1 || n * (n + 1) / 2 != NOC)  
                err = -2;
        }
        else if (PMOPT == 7) {          /* full adj. matrix */

            printf1("full adjacency matrix (variables).\n");

            n = (int)sqrt((double)NOC);
            if (n * n != NOC)
                err = -2;
        }
        if (err == -2) {
            printf1("Error in number of cases.\n");
            goto GDDFin;
        }
        GD_NDMAX = GD_NP = n;
        GD_NG = PMNV;

        if (gdd_alloc(GD_NG)) 
            goto GDDFin;

        for (i = 0; i < GD_NG; ++i)
            GD_EV[i] = PMVIdx[i];

        GD_TYP = PMOPT;
        if (gdd_setne(GD_TYP))
            goto GDDFin;

        GD_NOC = NOC;
    }
    else {                          /* ## list of matrices */
        
        printf1("full adjacency matrix (matrices).\n");
        if (PMNVTyp != 1) {
            printf1("Invalid list of matrix names.\n");
            goto GDDFin;
        }

        GD_NG = 0;         
        for (i = 0; i < PMNV; ++i) {
            j = PMVIdx[i];

            if (j < 0 || j >= MaxMat || MatAlloc[j] == 0) {
                printf1("Error: undefined matrix on right-hand side.\n");
                goto GDDFin;
            }
            if (MatRow[j] != MatCol[j]) {
                printf1("Error: can only use square matrices.\n");
                goto GDDFin;
            }
            if (i == 0)  
                n = MatRow[j];
            else if (MatRow[j] != n) {
                printf1("Error: all matrices must have the same order.\n");
                goto GDDFin;
            }
            GD_NG++;
        }
        if (GD_NG == 0) {
            printf1("Error: need at least one matrix.\n");
            goto GDDFin;
        }
        GD_NDMAX = GD_NP = n;

        if (gdd_alloc(GD_NG)) 
            goto GDDFin;

        for (i = 0; i < GD_NG; ++i)
            GD_EV[i] = PMVIdx[i];

        GD_TYP = PMOPT;
        if (gdd_setne(GD_TYP))
            goto GDDFin;

        GD_NOC = 0;
    }
    GD_TYP = PMOPT;

    printf1("\nMinimal level for valid edges: %g\n",GD_EVMin);
    printf1("Number of nodes: %d. Largest node number: %d\n",GD_NP,GD_NDMAX);
    printf1("Graph type %d ",GD_GT);
    if (GD_GT == 1)
        printf1("(undirected, unvalued).\n");
    else if (GD_GT == 2) {
        printf1("(undirected, valued, bidirectional edges: ");
        if (GD_GTT == 1)
            printf1("minimum).\n");
        else if (GD_GTT == 2)
            printf1("maximum).\n");
        else if (GD_GTT == 3)
            printf1("sum).\n");
    }
    else if (GD_GT == 3)
        printf1("(directed, unvalued).\n");
    else if (GD_GT == 4)
        printf1("(directed, valued).\n");

    printf1("\nGraph   edges   loops  max edge value  isolated nodes\n");
    for (i = 0; i < GD_NG; ++i)  
        printf1("%4d %8d %7d  %14.6f  %12d\n",i + 1,GD_NE[i],GD_NL[i],GD_EVMax[i],GD_IN[i]);

    printf1("\nSuccessfully created new graph data structure. ");
    printf1("Current memory: %d bytes.\n",MemReq);

    if (PMProtFDef)   
        gdd_prot();

    if (PMPERM == 1) {     /* permanently save the input data */

        if (GD_TYP != 1) {
            printf1("Warning: perm parameter requires option 1 (edge list).\n");
            printf1("Graph data not permanently saved.\n");
        }
        else {
            if (gdd_palloc(NOC,GD_NG))
                goto GDDFin;

            for (i = 0; i < NOC; ++i) {
                GD_EIP[i] = (int)get_data(GD_EI,i);
                GD_EJP[i] = (int)get_data(GD_EJ,i);
                k = i * GD_NG;
                for (j = 0; j < GD_NG; ++j) {
                    tmp = get_data(GD_EV[j],i);
                    if (tmp <= GD_EVMin)
                        tmp = -1.0;
                    GD_EP[k++] = (float)tmp;                       
                }
            }
            GD_PERM = 1;
            printf1("Saved data for permanent use. ");
            printf1("Current memory: %d bytes.\n",MemReq);
        }
    }
    err = 0;

GDDFin:
    if (err)  
        gdd_free(0);
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdd_info()  Info about current graph data                               */

void gdd_info(void)
{
    register int i;

    printf1("Current relational (graph) data.\n");
    if (GD_TYP == 0) {
        printf1("not defined.\n");
        return;
    }
    printf1("Number of nodes: %d. Largest node number: %d\n",GD_NP,GD_NDMAX);
    printf1("Graph type %d ",GD_GT);
    if (GD_GT == 1)
        printf1("(undirected, unvalued).\n");
    else if (GD_GT == 2) {
        printf1("(undirected, valued, bidirectional edges: ");
        if (GD_GTT == 1)
            printf1("minimum).\n");
        else if (GD_GTT == 2)
            printf1("maximum).\n");
        else if (GD_GTT == 3)
            printf1("sum).\n");
    }
    else if (GD_GT == 3)
        printf1("(directed, unvalued).\n");
    else if (GD_GT == 4)
        printf1("(directed, valued).\n");

    printf1("\nGraph   edges   loops  max edge value  isolated nodes\n");
    for (i = 0; i < GD_NG; ++i)  
        printf1("%4d %8d %7d  %14.6f  %12d\n",i + 1,GD_NE[i],GD_NL[i],GD_EVMax[i],GD_IN[i]);
}

/* ------------------------------------------------------------------------ */
/*  gdd_palloc(n,m)  allocate arrays for gdd1 option. n is number of cases  */
/*                   If n = 0 free previously allocated arrays.             */
/*                   m is number of graphs.                                 */
/*                   Return 0 if OK, -1 if error.                           */

int gdd_palloc(int n,int m)
{
    int err = -1;

    if (n == 0) {
        err = 0;
        goto GDD_PAFin;
    }
    if (n > 0) {
        if (!(GD_EIP = (int *)calloc(n,sizeof(int))))   
            goto GDD_PAFin;
        memrq(n,sizeof(int));
        GD_EIPA = n;

        if (!(GD_EJP = (int *)calloc(n,sizeof(int))))   
            goto GDD_PAFin;
        memrq(n,sizeof(int));
        GD_EJPA = n;

        if (!(GD_EP = (float *)calloc(n * m,sizeof(float))))   
            goto GDD_PAFin;
        memrq(n * m,sizeof(float));
        GD_EPA = n * m;
        return(0);
    }

GDD_PAFin:
    if (err)  
        printf1("Error: insufficient memory for saving graph data.\n");

    if (GD_EIPA > 0) {
        free((char *)GD_EIP);
        memrq(-GD_EIPA,sizeof(int));
        GD_EIPA = 0;
    }
    if (GD_EJPA > 0) {
        free((char *)GD_EJP);
        memrq(-GD_EJPA,sizeof(int));
        GD_EJPA = 0;
    }
    if (GD_EPA > 0) {
        free((char *)GD_EP);
        memrq(-GD_EPA,sizeof(float));
        GD_EPA = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdd_alloc(n) allocate GD_NE, GD_IN, GD_NL, GD_EV, GD_EVMax for n graphs */
/*               if n = 0 free previously allocated memory.                 */
/*               Return 0 if OK, -1 if error.                               */

int gdd_alloc(int n)  
{
    if (GD_NEA > 0) {
        free((char *)GD_NE);
        memrq(-GD_NEA,sizeof(int));
        GD_NEA = 0;
    }
    if (GD_INA > 0) {
        free((char *)GD_IN);
        memrq(-GD_INA,sizeof(int));
        GD_INA = 0;
    }
    if (GD_NLA > 0) {
        free((char *)GD_NL);
        memrq(-GD_NLA,sizeof(int));
        GD_NLA = 0;
    }
    if (GD_EVA > 0) {
        free((char *)GD_EV);
        memrq(-GD_EVA,sizeof(int));
        GD_EVA = 0;
    }
    if (GD_EVMaxA > 0) {
        free((char *)GD_EVMax);
        memrq(-GD_EVMaxA,sizeof(double));
        GD_EVMaxA = 0;
    }
    if (n == 0)
        return(0);

    if (!(GD_NE = (int *)calloc(n,sizeof(int)))) { 
        p_err(-2,1);
        return(-1);
    }
    memrq(n,sizeof(int));
    GD_NEA = n;

    if (!(GD_IN = (int *)calloc(n,sizeof(int)))) { 
        p_err(-2,1);
        return(-1);
    }
    memrq(n,sizeof(int));
    GD_INA = n;

    if (!(GD_NL = (int *)calloc(n,sizeof(int)))) { 
        p_err(-2,1);
        return(-1);
    }
    memrq(n,sizeof(int));
    GD_NLA = n;

    if (!(GD_EV = (int *)calloc(n,sizeof(int)))) { 
        p_err(-2,1);
        return(-1);
    }
    memrq(n,sizeof(int));
    GD_EVA = n;

    if (!(GD_EVMax = (double *)calloc(n,sizeof(double)))) { 
        p_err(-2,1);
        return(-1);
    }
    memrq(n,sizeof(double));
    GD_EVMaxA = n;

    return(0);
}

/* ------------------------------------------------------------------------ */
/*  gdd_free(opt)       free graph data structures, if opt !== message      */

void gdd_free(int opt)
{
    register int i;
    int n;

    gdd_alloc(0); 

    if (GD_NDA > 0) {
        free((char *)GD_ND);
        memrq(-GD_NDA,sizeof(int));
        GD_NDA = 0;
    }
    if (GD_APtrA > 0) {
        free((char *)GD_APtr);
        memrq(-GD_APtrA,sizeof(int));
        GD_APtrA = 0;
    }
    if (GD_FPNA > 0) {
        for (i = 0; i < GD_NP; ++i) {
            n = GD_FPN[i];
            if (n > 0) {
                if (GD_FPIA > 0) {
                    free((char *)GD_FPI[i]);
                    memrq(-n,sizeof(int));
                }
                if (GD_FPKA > 0) {
                    free((char *)GD_FPK[i]);
                    memrq(-n,sizeof(int));
                }
                if (GD_FPK1A > 0) {
                    free((char *)GD_FPK1[i]);
                    memrq(-n,sizeof(int));
                }
            }
        }
        free((char *)GD_FPN);
        memrq(-GD_FPNA,sizeof(int));
        GD_FPNA = 0;

        if (GD_FPIA > 0) {
            free((char *)GD_FPI);
            memrq(-GD_FPIA,sizeof(int *));
            GD_FPIA = 0;
        }
        if (GD_FPKA > 0) {
            free((char *)GD_FPK);
            memrq(-GD_FPKA,sizeof(int *));
            GD_FPKA = 0;
        }
        if (GD_FPK1A > 0) {
            free((char *)GD_FPK1);
            memrq(-GD_FPK1A,sizeof(int *));
            GD_FPK1A = 0;
        }
    }
    if (GD_BPNA > 0) {
        for (i = 0; i < GD_NP; ++i) {
            n = GD_BPN[i];
            if (n > 0) {
                if (GD_BPIA > 0) {
                    free((char *)GD_BPI[i]);
                    memrq(-n,sizeof(int));
                }
                if (GD_BPKA > 0) {
                    free((char *)GD_BPK[i]);
                    memrq(-n,sizeof(int));
                }
            }
        }
        free((char *)GD_BPN);
        memrq(-GD_BPNA,sizeof(int));
        GD_BPNA = 0;

        if (GD_BPIA > 0) {
            free((char *)GD_BPI);
            memrq(-GD_BPIA,sizeof(int *));
            GD_BPIA = 0;
        }
        if (GD_BPKA > 0) {
            free((char *)GD_BPK);
            memrq(-GD_BPKA,sizeof(int *));
            GD_BPKA = 0;
        }
    }
    gdd_palloc(0,0);

    if (GD_TYP && opt)
        printf1("Removed graph data.\n");

    GD_EVMin = 0.0;
    GD_PERM = GD_NG = GD_TYP = 0;
    GD_EI = GD_EJ = -1;
}

/* ------------------------------------------------------------------------ */
/*  gdd_setgt(typ,gt,gtt)   Setup GD_GT, GD_GTT for typ.                    */   

void gdd_setgt(int typ,int gt,int gtt)
{
    GD_GT = gt;
    GD_GTT = gtt;

    if (gt < 1 || gt > 4) {
        if (typ == 3 || typ == 4)
            GD_GT = 2;
        else
            GD_GT = 4;
    }
    if (typ >= 3 && typ <= 6 && GD_GT > 2)
        GD_GT = 2;
}

/* -###-------------------------------------------------------------------- */
/*  gdd_setnd()         Setup nodes list GD_ND, number of nodes GD_NP,      */
/*                      and maximal number of nodes, GD_NDMAX.              */
/*                      This function is called for opt 1 and 2.            */
/*                      Return 0 if OK, -1 if error.                        */

int gdd_setnd(void)
{
    register int i,j,k,l;
    int err,ei,ej,ip;
    double tmp;
                       
    err = -1;
    ei = PMVIdx[0];         /* index of variable for starting node */
    ej = PMVIdx[1];         /* index of variable for ending node */

    if (alloc_acn(2 * NOC))                 
        goto GDDNDFin;

    GD_NDMAX = ip = 0;
    for (k = 0; k < NOC; ++k) {
        i = (int)get_data(ei,k);
        j = (int)get_data(ej,k);
        if (i < 1 || j < 1) {
            printf1("Error: found negative, or zero, node number in case %d.\n",k + 1);
            goto GDDNDFin;
        }
        AcN[ip++] = i;
        AcN[ip++] = j;

        if (GD_NDMAX < i)
            GD_NDMAX = i;
        if (GD_NDMAX < j)
            GD_NDMAX = j;

        for (l = 0; l < GD_NG; ++l) {
            if ((tmp = get_data(GD_EV[l],k)) >= GD_EVMin) {
                GD_NE[l] += 1;
                if (i == j)
                    GD_NL[l] += 1;

                if (GD_EVMax[l] < tmp)
                    GD_EVMax[l] = tmp;
            }
        }
    }
    if (sorti(2 * NOC,AcN,0))
        goto GDDNDFin;

    k = 1;
    for (i = 1; i < 2 * NOC; ++i) {
        if (AcN[i] != AcN[k - 1])
            AcN[k++] = AcN[i]; 
    }
    GD_NP = k;

    if (!(GD_ND = (int *)calloc(GD_NP,sizeof(int)))) { 
        p_err(-2,1);
        goto GDDNDFin;
    }
    memrq(GD_NP,sizeof(int));
    GD_NDA = GD_NP;       
    for (k = 0; k < GD_NP; ++k)  
        GD_ND[k] = AcN[k];

    /* calculate number of isolated nodes */

    for (l = 0; l < GD_NG; ++l) {
           
        if (alloc_acn(GD_NDMAX + 1))                 
            goto GDDNDFin;

        for (k = 0; k < NOC; ++k) {
            if (get_data(GD_EV[l],k) >= GD_EVMin) {
                i = (int)get_data(ei,k);
                j = (int)get_data(ej,k);
                if (i != j)
                    AcN[i] = AcN[j] = 1;
            }
        }
        GD_IN[l] = 0;
        for (k = 0; k < GD_NP; ++k) {
            if (AcN[GD_ND[k]] == 0)
                GD_IN[l] += 1;
        }
    }
    err = 0;

GDDNDFin:
    alloc_acn(0);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  gdd_setptr(opt)     Setup pointer for gdd option 1.                     */
/*                      If opt=0 directed, opt != 0 undirected graph.       */
/*                      Return 0 if OK, -1 if error.                        */

int gdd_setptr(int opt)
{
    register int i,j,k,l;
    int err,n,m,ei,ej,ii,jj,ki,kj,n1,n2,nn,i1,i2,k1,k2;
    short vidx[2];


    err = -1;
    ei = PMVIdx[0];         /* index of variable for starting node */
    ej = PMVIdx[1];         /* index of variable for ending node */
    GD_FPMax = GD_BPMax = 0;

    /* forward pointers */

    if (!(GD_FPN = (int *)calloc(GD_NP,sizeof(int)))) { 
        p_err(-2,1);
        goto GDDPFin;
    }
    memrq(GD_NP,sizeof(int));
    GD_FPNA = GD_NP;       

    if (!(GD_FPI = (int **)calloc(GD_NP,sizeof(int *)))) { 
        p_err(-2,1);
        goto GDDPFin;
    }
    memrq(GD_NP,sizeof(int *));
    GD_FPIA = GD_NP;       

    if (!(GD_FPK = (int **)calloc(GD_NP,sizeof(int *)))) { 
        p_err(-2,1);
        goto GDDPFin;
    }
    memrq(GD_NP,sizeof(int *));
    GD_FPKA = GD_NP;       

    vidx[0] = ei;
    vidx[1] = ej;
    if (vsort(2,vidx,1,0,0))
        goto GDDPFin;
        

    k = VSORTPtr[0];
    i = (int)get_data(ei,k);
    j = (int)get_data(ej,k);
    for (l = 1; l < NOC; ++l) {
        k = VSORTPtr[l];
        ki = (int)get_data(ei,k);
        kj = (int)get_data(ej,k);
        if (ki == i && kj == j) {
            printf1("Error: edge %d - %d occurs two times.\n",i,j);
            goto GDDPFin;
        }
        i = ki;
        j = kj;
    }

    if (alloc_ack(NOC))
        goto GDDPFin;
    if (alloc_aci(NOC))
        goto GDDPFin;
   
    l = 0;
    k = VSORTPtr[0];
    ki = (int)get_data(ei,k);
    for (i = 0; i < GD_NP; ++i) {
        ii = GD_ND[i];
        n = 0;
        while (ki < ii && ++l < NOC) {
            k = VSORTPtr[l];
            ki = (int)get_data(ei,k);
        }
        if (ki == ii) {
            kj = (int)get_data(ej,k);
            AcI[n] = gdd_fndni(kj);
            AcK[n++] = k;
            while (++l < NOC) {
                k = VSORTPtr[l];
                ki = (int)get_data(ei,k);
                if (ki == ii) {
                    kj = (int)get_data(ej,k);
                    AcI[n] = gdd_fndni(kj);
                    AcK[n++] = k;
                }
                else
                    break;
            }
            if (!(GD_FPI[i] = (int *)calloc(n,sizeof(int)))) { 
                p_err(-2,1);
                goto GDDPFin;
            }
            if (!(GD_FPK[i] = (int *)calloc(n,sizeof(int)))) { 
                free((char *)GD_FPI[i]);
                p_err(-2,1);
                goto GDDPFin;
            }
            memrq(2 * n,sizeof(int));
            GD_FPN[i] = n;
            for (j = 0; j < n; ++j) {
                GD_FPI[i][j] = AcI[j];
                GD_FPK[i][j] = AcK[j];
            }
            if (GD_FPMax < n)
                GD_FPMax = n;
        }
    }
    vsort(2,vidx,0,0,0);      /* free */

    /* backward pointers */

    if (opt == 0) {

        if (!(GD_BPN = (int *)calloc(GD_NP,sizeof(int)))) { 
            p_err(-2,1);
            goto GDDPFin;
        }
        memrq(GD_NP,sizeof(int));
        GD_BPNA = GD_NP;       

        if (!(GD_BPI = (int **)calloc(GD_NP,sizeof(int *)))) { 
            p_err(-2,1);
            goto GDDPFin;
        }
        memrq(GD_NP,sizeof(int *));
        GD_BPIA = GD_NP;       

        if (!(GD_BPK = (int **)calloc(GD_NP,sizeof(int *)))) { 
            p_err(-2,1);
            goto GDDPFin;
        }
        memrq(GD_NP,sizeof(int *));
        GD_BPKA = GD_NP;       
    }
    else {
        if (!(GD_FPK1 = (int **)calloc(GD_NP,sizeof(int *)))) { 
            p_err(-2,1);
            goto GDDPFin;
        }
        memrq(GD_NP,sizeof(int *));
        GD_FPK1A = GD_NP;       

        if (alloc_acn(NOC))
            goto GDDPFin;
        if (alloc_acr(NOC))
            goto GDDPFin;
        if (alloc_acs(NOC))
            goto GDDPFin;
    }
    vidx[0] = ej;
    vidx[1] = ei;
    if (vsort(2,vidx,1,0,0))
        goto GDDPFin;
    
    l = 0;
    k = VSORTPtr[0];
    kj = (int)get_data(ej,k);
    for (j = 0; j < GD_NP; ++j) {
        jj = GD_ND[j];
        n = 0;
        while (kj < jj && ++l < NOC) {
            k = VSORTPtr[l];
            kj = (int)get_data(ej,k);
        }
        if (kj == jj) {
            ki = (int)get_data(ei,k);
            AcI[n] = gdd_fndni(ki);
            AcK[n++] = k;
            while (++l < NOC) {
                k = VSORTPtr[l];
                kj = (int)get_data(ej,k);
                if (kj == jj) {
                    ki = (int)get_data(ei,k);
                    AcI[n] = gdd_fndni(ki);
                    AcK[n++] = k;
                }
                else
                    break;
            }
            if (opt == 0) {

                if (!(GD_BPI[j] = (int *)calloc(n,sizeof(int)))) { 
                    p_err(-2,1);
                    goto GDDPFin;
                }
                if (!(GD_BPK[j] = (int *)calloc(n,sizeof(int)))) { 
                    free((char *)GD_BPI[j]);
                    p_err(-2,1);
                    goto GDDPFin;
                }
                memrq(2 * n,sizeof(int));
                GD_BPN[j] = n;
                for (i = 0; i < n; ++i) {
                    GD_BPI[j][i] = AcI[i];
                    GD_BPK[j][i] = AcK[i];
                }
                if (GD_BPMax < n)
                    GD_BPMax = n;
            }
        }
/*##*/  if (opt != 0) {         /* undirected graph */
            m = GD_FPN[j];

            if (n > 0 || m > 0) {
      
                i1 = i2 = nn = 0;
                while (1) {
                    if (i1 < m) {
                        n1 = GD_FPI[j][i1];
                        k1 = GD_FPK[j][i1++];
                    }
                    else
                        k1 = n1 = -1;

                    if (i2 < n) {
                        n2 = AcI[i2];
                        k2 = AcK[i2++];
                    }
                    else
                        k2 = n2 = -1;

                    if (n1 < 0 && n2 < 0)
                        break;
                    else if (n1 < 0) {
                        AcN[nn] = n2;
                        AcR[nn] = -1;
                        AcS[nn++] = k2;
                    }
                    else if (n2 < 0) {
                        AcN[nn] = n1;
                        AcR[nn] = k1;
                        AcS[nn++] = -1;
                    }
                    else if (n1 < n2) {
                        AcN[nn] = n1;
                        AcR[nn] = k1;
                        AcS[nn++] = -1;

                        while (i1 < m) {
                            n1 = GD_FPI[j][i1];              
                            k1 = GD_FPK[j][i1];              
                            if (n1 >= n2) {
                                if (n1 == n2) {
                                    AcN[nn] = n1;
                                    AcR[nn] = k1;
                                    AcS[nn++] = k2;
                                    i1++;
                                }
                                break;
                            }
                            AcN[nn] = n1;
                            AcR[nn] = k1;
                            AcS[nn++] = -1;
                            i1++;
                        }
                        if (n2 < n1 || i1 >= m) {
                            AcN[nn] = n2;
                            AcR[nn] = -1;
                            AcS[nn++] = k2;
                        }
                    }
                    else if (n2 < n1) {
                        AcN[nn] = n2;
                        AcR[nn] = -1;
                        AcS[nn++] = k2;

                        while (i2 < n) {
                            n2 = AcI[i2];           
                            k2 = AcK[i2];           
                            if (n2 >= n1) {
                                if (n2 == n1) {
                                    AcN[nn] = n2;
                                    AcR[nn] = k1;
                                    AcS[nn++] = k2;
                                    i2++;
                                }
                                break;
                            }
                            AcN[nn] = n2;
                            AcS[nn] = k2;
                            AcR[nn++] = -1;
                            i2++;
                        }
                        if (n1 < n2 || i2 >= n) {
                            AcN[nn] = n1;
                            AcR[nn] = k1;
                            AcS[nn++] = -1;
                        }
                    }
                    else {
                        AcN[nn] = n1;
                        AcR[nn] = k1;
                        AcS[nn++] = k2;
                    }
                }
                if (m > 0) {
                    free((char *)GD_FPI[j]);
                    free((char *)GD_FPK[j]);
                    memrq(-2 * m,sizeof(int));
                }
                if (!(GD_FPI[j] = (int *)calloc(nn,sizeof(int)))) { 
                    p_err(-2,1);
                    goto GDDPFin;
                }
                if (!(GD_FPK[j] = (int *)calloc(nn,sizeof(int)))) { 
                    free((char *)GD_FPI[j]);
                    p_err(-2,1);
                    goto GDDPFin;
                }
                if (!(GD_FPK1[j] = (int *)calloc(nn,sizeof(int)))) { 
                    free((char *)GD_FPI[j]);
                    free((char *)GD_FPK[j]);
                    p_err(-2,1);
                    goto GDDPFin;
                }
                memrq(3 * nn,sizeof(int));
                GD_FPN[j] = nn;
                for (i = 0; i < nn; ++i) {
                    GD_FPI[j][i] = AcN[i];
                    GD_FPK[j][i] = AcR[i];
                    GD_FPK1[j][i] = AcS[i];
                }
                if (GD_FPMax < nn)
                    GD_FPMax = nn;
            }
        }
    }
    err = 0;

GDDPFin:
    alloc_acs(0);   
    alloc_acr(0);   
    alloc_acn(0);   
    alloc_aci(0);   
    alloc_ack(0);   
    vsort(2,vidx,0,0,0);      /* free */
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdd_setap()         Setup GD_APtr for gdd option 2.                     */
/*                      Return 0 if OK, -1 if error.                        */

int gdd_setap(void)
{
    register int i,j,k,ii,jj;
    int err,ei,ej;


    err = -1;

    ei = PMVIdx[0];         /* index of variable for starting node */
    ej = PMVIdx[1];         /* index of variable for ending node */

    if (!(GD_APtr = (int *)calloc(GD_NP * GD_NP,sizeof(int)))) { 
        p_err(-2,1);
        goto GDDAPFin;
    }
    GD_APtrA = GD_NP * GD_NP;
    memrq(GD_APtrA,sizeof(int));

    for (k = 0; k < NOC; ++k) {
        i = (int)get_data(ei,k);
        j = (int)get_data(ej,k);
        ii = gdd_fndni(i);
        jj = gdd_fndni(j);
        if (ii < 0 || jj < 0) {
            printfe("ERROR in gdd [i=%d j=%d ii=%d jj=%d]\n",i,j,ii,jj);
            gerr_exit(215);
        }
        if (GD_APtr[ii * GD_NP + jj] > 0) {
            printf1("Error: edge (%d,%d) occurs twice.\n",i,j);
            goto GDDAPFin;
        }
        GD_APtr[ii * GD_NP + jj] = k + 1;
    }
    err = 0;

GDDAPFin:
    return(err);
}

/* -###-------------------------------------------------------------------- */
/*  gdd_setne(opt)      Setup GD_EVMax, GD_NE, GD_NL for otions 3 - 7 or 8. */
/*                      Return 0 if OK, -1 if error.                        */

int gdd_setne(int opt)
{
    register int i,j,k;
    int ne,nl;
    double tmp,emax;

    for (k = 1; k <= GD_NG; ++k) {

        if (alloc_acn(GD_NDMAX + 1))                 
            return(-1);        

        ne = nl = 0;
        emax = 0.0;
        for (i = 0; i < GD_NP; ++i) {
            for (j = 0; j < GD_NP; ++j) {
                if (opt == 3) {
                    if (j >= i)
                        continue;
                }
                else if (opt == 4) {
                    if (j <= i)
                        continue;
                }
                else if (opt == 5) {
                    if (j > i)
                        continue;
                }
                else if (opt == 6) {
                    if (j < i)
                        continue;
                }
                tmp = gdd_adj(i,j,k);
                if (tmp >= GD_EVMin) {
                    ne++;
                    if (i == j)
                        nl++;
                    emax = dmax(emax,tmp);
                    if (i != j)
                        AcN[i] = AcN[j] = 1;
                }
            }
        }
        GD_NE[k - 1] = ne;
        GD_NL[k - 1] = nl;
        GD_EVMax[k - 1] = emax;

        GD_IN[k - 1] = 0;                   /* number of isolated nodes */
        for (i = 0; i < GD_NP; ++i) {
            if (AcN[i] == 0)
                GD_IN[k - 1] += 1;
        }
    }
    alloc_acn(0);                
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  gdd_prot()  print current gdd data structure to protocol file.          */

void gdd_prot(void)
{     
    register int i,j,n;

    fprintf(PMProtFd,"Current gdd data structure\n\n");
    fprintf(PMProtFd,"GD_TYP   = %d\n",GD_TYP);
    fprintf(PMProtFd,"GD_GT    = %d\n",GD_GT);
    fprintf(PMProtFd,"GD_GTT   = %d\n",GD_GTT);
    fprintf(PMProtFd,"GD_NOC   = %d\n",GD_NOC);
    fprintf(PMProtFd,"GD_NG    = %d\n",GD_NG);
    fprintf(PMProtFd,"GD_NP    = %d\n",GD_NP);
    fprintf(PMProtFd,"GD_NDMAX = %d\n",GD_NDMAX);

    fprintf(PMProtFd,"\nGD_EV\n");
    for (i = 0; i < GD_NG; ++i)
        fprintf(PMProtFd,"%5d %5d\n",i,GD_EV[i]);

    fprintf(PMProtFd,"\nGD_EVMax\n");
    for (i = 0; i < GD_NG; ++i)
        fprintf(PMProtFd,"%5d     %g\n",i,GD_EVMax[i]);

    fprintf(PMProtFd,"\nGD_NE\n");
    for (i = 0; i < GD_NG; ++i)
        fprintf(PMProtFd,"%5d %5d\n",i,GD_NE[i]);

    fprintf(PMProtFd,"\nGD_NL\n");
    for (i = 0; i < GD_NG; ++i)
        fprintf(PMProtFd,"%5d %5d\n",i,GD_NL[i]);

    fprintf(PMProtFd,"\nGD_IN\n");
    for (i = 0; i < GD_NG; ++i)
        fprintf(PMProtFd,"%5d %5d\n",i,GD_IN[i]);

    if (GD_TYP <= 2) {
        fprintf(PMProtFd,"\nGD_ND\n");
        for (i = 0; i < GD_NP; ++i) {
            fprintf(PMProtFd,"%5d ",i);
            fprintf(PMProtFd,"%5d ",GD_ND[i]);
            fprintf(PMProtFd,"\n");
        }
    }
    if (GD_TYP == 1) {
        fprintf(PMProtFd,"\nForward pointer (GD_FPMax=%d)\n",GD_FPMax);
        for (i = 0; i < GD_NP; ++i) {
            n = GD_FPN[i];
            fprintf(PMProtFd,"%5d %5d : ",i,n);
            for (j = 0; j < n; ++j) {
                fprintf(PMProtFd,"%5d [%5d",GD_FPI[i][j],GD_FPK[i][j]);
                if (GD_GT <= 2)
                    fprintf(PMProtFd,",%5d] ",GD_FPK1[i][j]);
                else
                    fprintf(PMProtFd,"] ");
            }
            fprintf(PMProtFd,"\n");
        }
        if (GD_GT > 2) {
            fprintf(PMProtFd,"\nBackward pointer (GD_BPMax=%d)\n",GD_BPMax);
            for (i = 0; i < GD_NP; ++i) {
                n = GD_BPN[i];
                fprintf(PMProtFd,"%5d %5d : ",i,n);
                for (j = 0; j < n; ++j)
                    fprintf(PMProtFd,"%5d [%5d] ",GD_BPI[i][j],GD_BPK[i][j]);
                fprintf(PMProtFd,"\n");
            }
        }
    }
    if (GD_TYP == 2) {
        fprintf(PMProtFd,"\nGD_APtr\n");
        for (i = 0; i < GD_NP; ++i) {
            fprintf(PMProtFd,"%5d : ",i);
            for (j = 0; j < GD_NP; ++j)  
                fprintf(PMProtFd,"%4d ",GD_APtr[i * GD_NP + j]);
            fprintf(PMProtFd,"\n");
        }
    }
    fprintf(PMProtFd,"\n");
}

/* ------------------------------------------------------------------------ */
/*  gdd_fndni(k)    Find internal node number for external node number k.   */
/*                  return -1 if not found.                                 */

int gdd_fndni(int k)
{
    register int i,j,ij;

    if (GD_TYP > 2) {
        if (k < 1 || k > GD_NP)
            return(-1);
        return(k - 1);
    }
    i = 0;
    j = GD_NP - 1;
    while (i <= j) {
        ij = (i + j) / 2;
        if (k < GD_ND[ij])
            j = ij - 1;
        else if (k > GD_ND[ij])
            i = ij + 1;
        else  
            break;
    }
    if (k != GD_ND[ij])
        return(-1);
    return(ij);
}

/* ------------------------------------------------------------------------ */
/*  check_gdj(j)    Return 1 if variable j is required for current gdd      */
/*                  data structure, otherwise return 0.                     */

int check_gdj(int j)
{
    register int i;

    if (GD_TYP == 0 || GD_PERM != 0)
        return(0);
    if (j == GD_EI || j == GD_EJ)
        return(1);
    for (i = 0; i < GD_NG; ++i) {
        if (j == GD_EV[i])
            return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  gdd_adj(i,j,gn)     return value of edge (i,j) for graph number gn.     */
/*                                                                          */
/*  Return -1.0 if not available. Use GD_GT and GD_GTT.                     */
/*                                                                          */
/*  NOTE: this function uses internal node number: 0,...,GD_NP - 1.         */
/*  and requires that gd_setup() has been run correctly.                    */

double gdd_adj(int i,int j,int gn)
{
    register int k,l,n;
    double tmp1,tmp2;

    if (GD_TYP == 0)
        return(-1.0);

    if (i < 0 || j < 0 || i >= GD_NP || j >= GD_NP || gn < 1 || gn > GD_NG)
        return(-1.0);

    tmp1 = tmp2 = -1.0;

    switch (GD_TYP) {
        case 1:     n = GD_FPN[i];
                    for (k = 0; k < n; ++k) {
                        if (j == GD_FPI[i][k]) {
                            l = GD_FPK[i][k];
                            if (l >= 0)
                                tmp1 = gdd_ev(l,gn);

                            if (GD_GT <= 2) {
                                if (i == j)
                                    tmp2 = tmp1;
                                else  
                                    l = GD_FPK1[i][k];
                                if (l >= 0)
                                    tmp2 = gdd_ev(l,gn);
                            }
                            break;
                        }
                    }
                    break;

        case 2:     l = GD_APtr[i * GD_NP + j];
                    if (l > 0)
                        tmp1 = gdd_ev(l - 1,gn);

                    if (GD_GT <= 2) {
                        if (i == j)
                            tmp2 = tmp1;
                        else {
                            l = GD_APtr[j * GD_NP + i];
                            if (l > 0)
                                tmp2 = gdd_ev(l - 1,gn);
                        }
                    }
                    break;

        case 3:     if (j != i) {
                        if (j < i)
                            l = i * (i - 1) / 2 + j;
                        else               
                            l = j * (j - 1) / 2 + i;

                        tmp1 = tmp2 = gdd_ev(l,gn);
                    }
                    break;

        case 4:     if (j != i) {
                        if (i < j)
                            l = i * GD_NP - i * (i + 1) / 2 + j - i - 1;
                        else
                            l = j * GD_NP - j * (j + 1) / 2 + i - j - 1;

                        tmp1 = tmp2 = gdd_ev(l,gn);
                    }
                    break;


        case 5:     if (j <= i)  
                        l = i * (i + 1) / 2 + j;
                    else
                        l = j * (j + 1) / 2 + i;
                    tmp1 = tmp2 = gdd_ev(l,gn);
                    break;

        case 6:     if (j >= i)  
                        l = i * GD_NP - i * (i - 1) / 2 + j - i;
                    else
                        l = j * GD_NP - j * (j - 1) / 2 + i - j;
                    tmp1 = tmp2 = gdd_ev(l,gn);
                    break;

        case 7:     l = i * GD_NP + j;
                    tmp1 = gdd_ev(l,gn);

                    if (GD_GT <= 2) {
                        l = j * GD_NP + i;
                        tmp2 = gdd_ev(l,gn);
                    }
                    break;
                
        case 8:     l = i * GD_NP + j + 1;
                    tmp1 = gdd_ev(l,gn);

                    if (GD_GT <= 2) {
                        l = j * GD_NP + i + 1;
                        tmp2 = gdd_ev(l,gn);
                    }
                    break;

    }

    if (GD_GT == 1) {
        if (tmp1 >= 0.0 || tmp2 >= 0.0)
            return(1.0);
        else
            return(-1.0);
    }
    else if (GD_GT == 2) {
        if (tmp1 < 0.0)
            return(tmp2);
        else if (tmp2 < 0.0)
            return(tmp1);
        else if (GD_GTT == 1)
            return(dmin(tmp1,tmp2));
        else if (GD_GTT == 2)
            return(dmax(tmp1,tmp2));
        else if (i == j)
            return(tmp1);
        else
            return(tmp1 + tmp2);
    }
    else if (GD_GT == 3) {
        if (tmp1 < 0.0)
            return(-1.0);
        else
            return(1.0);
    }
    else
        return(tmp1);
}

/* ------------------------------------------------------------------------ */
/*  g_suc(i)    return number of successors of node i, without loops        */

int g_suc(int i,int gn)
{
    register int j,dp;
    int nn,n;

    if (GD_TYP == 0)
        return(0);
           
    if (i < 0 || i >= GD_NP || gn < 1 || gn > GD_NG)
        return(0);

    nn = 0;
    if (GD_TYP == 1) {
        n = GD_FPN[i];
        for (j = 0; j < n; ++j) {
            if (GD_FPI[i][j] == i)   
                continue;
     
            dp = GD_FPK[i][j];
            if (dp >= 0 && gdd_ev(dp,gn) >= 0.0)  
                    nn++;
            else if (GD_GT <= 2) {
                dp = GD_FPK1[i][j];
                if (dp >= 0 && gdd_ev(dp,gn) >= 0.0)
                    nn++;
            }
        }
    }
    else {
        for (j = 0; j < GD_NP; ++j) {
            if (i != j && gdd_adj(i,j,gn) >= 0.0)
                nn++;
        }
    }
    return(nn);
}

/* ------------------------------------------------------------------------ */
/*  g_pre(i)    return number of predecessors of node i, graph gn           */

int g_pre(int i,int gn)
{
    register int j,dp;
    int nn,n;

    if (GD_TYP == 0)
        return(0);

    if (i < 0 || i >= GD_NP || gn < 1 || gn > GD_NG)
        return(0);

    nn = 0;
    if (GD_TYP == 1) {
        if (GD_GT <= 2)
            nn = g_suc(i,gn);
        else {
            n = GD_BPN[i];
            for (j = 0; j < n; ++j) {
                if (GD_BPI[i][j] == i)
                    continue;
                dp = GD_BPK[i][j];
                if (dp >= 0 && gdd_ev(dp,gn) >= 0.0)
                    nn++;
            }
        }
    }
    else {
        for (j = 0; j < GD_NP; ++j) {
            if (j != i && gdd_adj(j,i,gn) >= 0.0)
                nn++;
        }
    }
    return(nn);
}

/* ------------------------------------------------------------------------ */
/*  g_loop(i)   return 1 if node i has a loop, graph number gn.             */

int g_loop(int i,int gn)
{
    if (gdd_adj(i,i,gn) >= 0.0)
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  g_getnd(n,nd)   Get list of input nodes from PMF2d into nd[], max       */
/*                  number of nodes is n. If n = 0 expect an array          */
/*                  nd[0,...,GDU_NP-1] and sets flags.                      */
/*                                                                          */
/*                  Return:  -1 if error, otherwise number of nodes         */

int g_getnd(int n,int *nd)
{
    int nn,k;
    char buf[101],*p;

    if (PMF2Def == 0) {
        printf1("Error: need if parameter for list of nodes.\n");
        return(-1);
    }
    nn = 0;
    while (fgets(buf,100,PMF2d)) {
        p = skip_b(buf);
        if (sscanf(p,"%d",&k) == 1 && k >= 1) {

            if (GD_TYP <= 2) {
                if (k <= GD_NDMAX) {                         
                    k = gdd_fndni(k);
                    if (k < 0)
                        continue;
                }
            }
            else {      
                if (k > GD_NP)
                    continue;
                k--;
            }

            if (n > 0) {
                if (nn >= n) {
                    printf1("Error: exceeded maximal number of nodes for input (%d).\n",n);
                    return(-1);
                }
                nd[nn++] = k; 
            }
            else {
                if (nd[k] == 0) {
                    nd[k] = 1;
                    nn++;
                }
            }
        }
    }
    return(nn);
}

/* ------------------------------------------------------------------------ */
/*  g_getne(n,ni,nj)   Get list of edges from PMF2d into ni[], nj[].        */
/*                     Max number of edges is n.                            */
/*                                                                          */
/*                  Return:  -1 if error, otherwise number of nodes         */

int g_getne(int n,int *ni,int *nj)
{
    int nn,k1,k2;
    char buf[101],*p;

    if (PMF2Def == 0) {
        printf1("Error: need if parameter for list of edges.\n");
        return(-1);
    }
    nn = 0;
    while (fgets(buf,100,PMF2d)) {
        p = skip_b(buf);
        if (sscanf(p,"%d",&k1) != 1 || k1 < 1 || k1 > GD_NDMAX)
            continue;
        p = skip_int(p);
        p = skip_b(p);
        if (sscanf(p,"%d",&k2) != 1 || k2 < 1 || k2 > GD_NDMAX)
            continue;

        if (nn >= n) {
            printf1("Error: exceeded maximal number of edges for input (%d).\n",n);
            return(-1);
        }
        if (GD_TYP > 2) {
            k1--;
            k2--;
        }
        else {
            k1 = gdd_fndni(k1);
            if (k1 < 0)
                continue;

            k2 = gdd_fndni(k2);
            if (k2 < 0)
                continue;
        }
        ni[nn] = k1;
        nj[nn++] = k2;
    }
    return(nn);
}

/* ------------------------------------------------------------------------ */
/*  gdd_check(gn,opt)       Check for gdd data.                             */
/*                          If opt = 0 and gn = 0 select gn = 1, if         */
/*                          opt = 1 and gn = 0, don't change.               */
/*                          Set PMGN.                                       */
/*                          Return 0 if OK, -1 if error.                    */

int gdd_check(int gn,int opt)
{
    if (GD_TYP == 0) {
        printf1("Error: no graph data defined.\n");
        return(-1);
    }
    if (gn > GD_NG) {
        printf1("Error: graph number %d not defined.\n",gn);
        return(-1);
    }
    PMGN = gn;  
    if (gn == 0) {
        if (opt == 0 || GD_NG == 1)
            PMGN = 1;
    }
    if (PMGN == 0)
        printf1("Using multigraph ");
    else
        printf1("Using graph number %d ",PMGN);

    if (GD_GT == 1)
        printf1("(undirected, unvalued).\n");
    else if (GD_GT == 2) {
        printf1("(undirected, valued, bidirectional edges: ");
        if (GD_GTT == 1)
            printf1("minimum).\n");
        else if (GD_GTT == 2)
            printf1("maximum).\n");
        else if (GD_GTT == 3)
            printf1("sum).\n");
    }
    else if (GD_GT == 3)
        printf1("(directed, unvalued).\n");
    else if (GD_GT == 4)
        printf1("(directed, valued).\n");
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  gdd_check2(gn,gn1,opt)  Check for gdd data. Two graphs, gn and gn1.     */
/*                          If opt = 0 and gn = 0 and gn1 == 0              */
/*                          select gn=1, gn1=2.                             */
/*                          Set PMGN and PMGN1.                             */
/*                          Return 0 if OK, -1 if error.                    */

int gdd_check2(int gn,int gn1,int opt)
{
    if (GD_TYP == 0) {
        printf1("Error: no graph data defined.\n");
        return(-1);
    }
    if (gn < 0)
        gn = 0;
    if (gn1 < 0)
        gn1 = 0;

    if (gn == 0 && gn1 == 0 && opt == 0) {
        gn = 1;
        gn1 = 2;
    }
    if (gn > GD_NG || gn1 > GD_NG) {
        printf1("Error: at least one graph number not defined.\n");
        return(-1);
    }
    PMGN = gn;  
    PMGN1 = gn1;
    printf1("Using graphs %d and %d ",PMGN,PMGN1);

    if (GD_GT == 1)
        printf1("(undirected, unvalued).\n");
    else if (GD_GT == 2) {
        printf1("(undirected, valued, bidirectional edges: ");
        if (GD_GTT == 1)
            printf1("minimum).\n");
        else if (GD_GTT == 2)
            printf1("maximum).\n");
        else if (GD_GTT == 3)
            printf1("sum).\n");
    }
    else if (GD_GT == 3)
        printf1("(directed, unvalued).\n");
    else if (GD_GT == 4)
        printf1("(directed, valued).\n");
    return(0);
}
   
/* ------------------------------------------------------------------------ */
/*  gdd_tcheck(typ,dir,val)   Check graph type.                             */
/*                                                                          */
/*                        typ = 1: graph must have GD_TYP = 1               */
/*                        dir = 1: graph must be undirected                 */
/*                        dir = 2: graph must be directed                   */
/*                        val = 1: graph must be unvalued                   */
/*                        val = 2: graph must be valued                     */
/*                        Return 0 if OK, -1 if error.                      */

int gdd_tcheck(int typ,int dir,int val)
{
    if (typ == 1) {
        if (GD_TYP != 1) {
            printf1("Error: need a graph defined with gdd option 1 (edge list).\n");
            return(-1);
        }
    }
    if (dir == 1) {
        if (GD_GT > 2) {
            printf1("Error: graph must be undirected.\n");
            return(-1);
        }
    }
    else if (dir == 2) {
        if (GD_GT <= 2) {
            printf1("Error: graph must be directed.\n");
            return(-1);
        }
    }
    if (val == 1) {
        if (GD_GT != 1 && GD_GT != 3) {
            printf1("Error: graph must be unvalued.\n");
            return(-1);
        }
    }
    else if (val == 2) {
        if (GD_GT != 2 && GD_GT != 4) {
            printf1("Error: graph must be valued.\n");
            return(-1);
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  gdd_ei(k)   Return node number in k.th data matrix row.                 */

int gdd_ei(int k)              
{
    if (GD_PERM)
        return(GD_EIP[k]);
    return((int)get_data(GD_EI,k));
}

/* ------------------------------------------------------------------------ */
/*  gdd_ej(k)   Return node number in k.th data matrix row.                 */

int gdd_ej(int k)              
{
    if (GD_PERM)
        return(GD_EJP[k]);
    return((int)get_data(GD_EJ,k));
}

/* ------------------------------------------------------------------------ */
/*  gdd_ev(k,gn)  Return edge value in k.th data matrix row, graph gn.      */
/*                Graphs are counted 1,...,GD_GN - 1.                       */
               
double gdd_ev(int k,int gn)       
{
    double tmp;

    if (k < 0)
        return(-1.0);
              
    if (GD_PERM)
        return((double)GD_EP[k * GD_NG + gn - 1]);

    if (GD_TYP == 8)  
        tmp = MatVal[GD_EV[gn - 1]][k];
    else
        tmp = get_data(GD_EV[gn - 1],k); 

    if (tmp < GD_EVMin)
        tmp = -1.0;
    return(tmp); 
}

/* ------------------------------------------------------------------------ */
/*  gdd_node(i)     Return external node number corresponding to i.         */

int gdd_node(int i)              
{
    if (GD_TYP == 1 || GD_TYP == 2)
        return(GD_ND[i]);
    return(i + 1);
}

/* ------------------------------------------------------------------------ */
/*  n_info      Print info about processed nodes to stderr.                 */
  
void n_info(int i,int n)
{
    register int j;

    if (i > 0 && SILENTFlg < 2 && GD_NP > 500) {
        j = i / n;
        if (j * n == i) {        
            printfe("Node: %7d     %c",i,CR);
            fflushe();
        }
    }
}          
void n_info_e(void)
{
    printfe("\n");
    fflushe();
}          

/* ------------------------------------------------------------------------ */
/*  gcd         Create graphs with randomly selected edges.                 */
/*                                                                          */
/*              gcd(                                                        */
/*                  opt=...,    1 = complete graph without loops            */
/*                              2 = complete graph with loops               */
/*                              3 = randomly select m edges                 */
/*                  n = ...,    number of nodes, def. 10                    */
/*                  m = ...,    number of edges, def. 1                     */
/*                  nfmt=...,   print format for nodes numbers, def. 4      */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gcd(void)
{
    register int i,j;
    int err,nrec,ne,ni;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Create unvalued graph. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto GCDFin;

    if (PMN < 1)
        PMN = 10;
    if (PMM < 0)
        PMM = 0;
    else if (PMM > PMN * PMN)
        PMM = PMN * PMN;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);
    if (PMOPT > 3)
        PMOPT = 3;

    if (PMOPT == 3 && PMM == PMN * PMN)  
        PMOPT = 2;

    printf1("Graph with %d nodes, ",PMN);
    if (PMOPT == 1)
        printf1("complete without loops.\n");
    else if (PMOPT == 2)
        printf1("complete, including loops.\n");
    else  
        printf1("%d edges (randomly selected).\n",PMM);

    ni = nrec = 0;
    if (PMOPT <= 2) {
        for (i = 1; i <= PMN; ++i) {
            for (j = 1; j <= PMN; ++j) {
                if (i == j && PMOPT ==1)
                    continue;

                fprintf(PMFd,PMNFmtS,i);
                fprintf(PMFd,PMNFmtS,j);
                fprintf(PMFd,PMNFmtS,1);
                fprintf(PMFd,"\n");
                nrec++;
            }
        }
    }
    else if (PMM == 0) {
        for (i = 1; i <= PMN; ++i) {
            fprintf(PMFd,PMNFmtS,i);
            fprintf(PMFd,PMNFmtS,i);
            fprintf(PMFd,PMNFmtS,-1);
            fprintf(PMFd,"\n");
            nrec++;
        }
    }
    else {                          /* random selection */
         if (alloc_aci(PMN))
            goto GCDFin;

        ne = 0;
        while (ne < PMM) {
            j = (int)((double)PMN * random1());
            if (j >= PMN)
                j = PMN - 1;
            if (AcI[j] < PMN) {
                AcI[j] += 1;
                ne++;
            }

        }
      
        for (i = 0; i < PMN; ++i) {
            if (AcI[i] == 0) {
                fprintf(PMFd,PMNFmtS,i + 1);
                fprintf(PMFd,PMNFmtS,i + 1);
                fprintf(PMFd,PMNFmtS,-1);
                fprintf(PMFd,"\n");
                nrec++;
                ni++;
            }
            else if (AcI[i] == PMN) {
                for (j = 1; j <= PMN; ++j) {
                    fprintf(PMFd,PMNFmtS,i + 1);
                    fprintf(PMFd,PMNFmtS,j);
                    fprintf(PMFd,PMNFmtS,1);
                    fprintf(PMFd,"\n");
                    nrec++;
                }
            }
            else {
                if (alloc_acc(PMN))
                    goto GCDFin;
                ne = 0;
                while (ne < AcI[i]) {
                    j = (int)((double)PMN * random1());
                    if (j >= PMN)
                        j = PMN - 1;
                    if (AcC[j] == 0) {
                        AcC[j] = 1;
                        ne++;
                    }
                }
                for (j = 0; j < PMN; ++j) {
                    if (AcC[j]) {
                        fprintf(PMFd,PMNFmtS,i + 1);
                        fprintf(PMFd,PMNFmtS,j + 1);
                        fprintf(PMFd,PMNFmtS,1);
                        fprintf(PMFd,"\n");
                        nrec++;
                    }
                }
            }
            if (i > 0 && SILENTFlg < 2 && PMN > 100) {
                j = i / 100;
                if (j * 100 == i) {        
                    printfe("Node: %7d     %c",i,CR);
                    fflushe();
                }
            }
        }
        if (SILENTFlg < 2 && PMN > 100) {
            printfe("\n");
            fflushe();
        }
    }
    printf1("%d records written to: %s\n",nrec,PMFdName);
    if (PMOPT == 3 && ni > 0)
        printf1("Graph contains %d isolated nodes.\n",ni);    

    err = 0;

GCDFin:       
    p_clean();
    return(err);
}

