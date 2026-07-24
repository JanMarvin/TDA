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
#include "tda_context.h"

/*  functions in t_gdd.c */

int gdd(TDAContext *ctx);
void gdd_info(TDAContext *ctx);
int gdd_palloc(TDAContext *ctx, int n,int m);
int gdd_alloc(TDAContext *ctx, int n); 
void gdd_free(TDAContext *ctx, int opt);
void gdd_setgt(TDAContext *ctx, int typ,int gt,int gtt);
int gdd_setnd(TDAContext *ctx);
int gdd_setptr(TDAContext *ctx, int opt);
int gdd_setap(TDAContext *ctx);
int gdd_setne(TDAContext *ctx, int opt);
void gdd_prot(TDAContext *ctx);
int gdd_fndni(TDAContext *ctx, int k);
int check_gdj(TDAContext *ctx, int j);
double gdd_adj(TDAContext *ctx, int i,int j,int gn);
int g_suc(TDAContext *ctx, int i,int gn);
int g_pre(TDAContext *ctx, int i,int gn);
int g_loop(TDAContext *ctx, int i,int gn);
int g_getnd(TDAContext *ctx, int n,int *nd);
int g_getne(TDAContext *ctx, int n,int *ni,int *nj);
int gdd_check(TDAContext *ctx, int gn,int opt);
int gdd_check2(TDAContext *ctx, int gn,int gn1,int opt);
int gdd_tcheck(TDAContext *ctx, int typ,int dir,int val);
int gdd_ei(TDAContext *ctx, int k);             
int gdd_ej(TDAContext *ctx, int k);             
double gdd_ev(TDAContext *ctx, int k,int gn);      
int gdd_node(TDAContext *ctx, int i);             
void n_info(TDAContext *ctx, int i,int n);
void n_info_e(TDAContext *ctx);
int gcd(TDAContext *ctx);





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

int gdd(TDAContext *ctx)
{
    register int i,j,k;
    int err,n = 0;     
    double tmp;

    err = -1;

    if (check_cmd(ctx, 1))
        return(-1);

    if (!strcmp(ctx->CmdBuf,"gdd")) {
        gdd_info(ctx);
        return(0);
    }

    if (!strcmp(ctx->CmdBuf,"gdd=off")) {
        if (ctx->GD_TYP) {
            gdd_free(ctx, 0);
            printf1(ctx, "Removed graph data structure. Current memory: %d bytes.\n",ctx->MemReq);
        }
        err = 0;
        goto GDDFin;
    }
    gdd_free(ctx, 0);         /* free previously defined graph data */

    printf1(ctx, "New graph data. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 3,14,1))     /* get parameters */
        goto GDDFin;

    if (ctx->PMOPT < 1 || ctx->PMOPT > 8) {
        printf1(ctx, "Error: possible options are: 1 - 8.\n");
        goto GDDFin;
    }
    if (ctx->PMOPT != 8 && ctx->PMNVTyp == 1) {
        printf1(ctx, "Error: matrix names only with option 8.\n");
        goto GDDFin;
    }
    ctx->GD_EVMin = ctx->PMSC;
    if (ctx->GD_EVMin < 0.0) {
        printf1(ctx, "Error: minimal value for valid edges must not be negative.\n");
        goto GDDFin;
    }

    gdd_setgt(ctx, ctx->PMOPT,ctx->PMGT,ctx->PMGTT);

    printf1(ctx, "Option %d: ",ctx->PMOPT);

    if (ctx->PMOPT <= 2) {               /* setup node list: GD_ND */

        if (ctx->PMOPT == 1)
            printf1(ctx, "edge list.\n");
        else
            printf1(ctx, "pointer to adjacency matrix.\n");

        if (ctx->PMNV < 3) {
            printf1(ctx, "Error: this option requires at least three variables.\n");
            goto GDDFin;
        }
        ctx->GD_EI = ctx->PMVIdx[0];          /* index of variable for starting node */
        ctx->GD_EJ = ctx->PMVIdx[1];          /* index of variable for ending node */
        ctx->GD_NG = ctx->PMNV - 2;           /* number of graphs */

        if (gdd_alloc(ctx, ctx->GD_NG)) 
            goto GDDFin;

        for (i = 0; i < ctx->GD_NG; ++i)
            ctx->GD_EV[i] = ctx->PMVIdx[i + 2];

        if (gdd_setnd(ctx))            /* setup node list */
            goto GDDFin;

        if (ctx->PMOPT == 1) {           /* setup pointer */
            n = 0;
            if (ctx->GD_GT <= 2)
                n = 1;
            if (gdd_setptr(ctx, n))     
                goto GDDFin;
        }
        else {
            if (gdd_setap(ctx))     
                goto GDDFin;
        }
        ctx->GD_NOC = ctx->NOC;
    }
    else if (ctx->PMOPT <= 7) {

        if (ctx->PMOPT == 3 || ctx->PMOPT == 4) {           /* lower/upper triangle */
            if (ctx->PMOPT == 3)
                printf1(ctx, "lower triangle.\n");
            else
                printf1(ctx, "upper triangle.\n");

            n = (int)((sqrt(8.0 * (double)ctx->NOC + 1.0) + 1.0) / 2.0);
            if (n < 1 || n * (n - 1) / 2 != ctx->NOC)  
                err = -2; 
        }
        else if (ctx->PMOPT == 5 || ctx->PMOPT == 6) {    /* lower/upper triangle plus diagonal */
            if (ctx->PMOPT == 5)
                printf1(ctx, "lower");
            else
                printf1(ctx, "upper");
            printf1(ctx, " triangle plus diagonal.\n");

            n = (int)((sqrt(8.0 * (double)ctx->NOC - 1.0) + 1.0) / 2.0);
            if (n < 1 || n * (n + 1) / 2 != ctx->NOC)  
                err = -2;
        }
        else if (ctx->PMOPT == 7) {          /* full adj. matrix */

            printf1(ctx, "full adjacency matrix (variables).\n");

            n = (int)sqrt((double)ctx->NOC);
            if (n * n != ctx->NOC)
                err = -2;
        }
        if (err == -2) {
            printf1(ctx, "Error in number of cases.\n");
            goto GDDFin;
        }
        ctx->GD_NDMAX = ctx->GD_NP = n;
        ctx->GD_NG = ctx->PMNV;

        if (gdd_alloc(ctx, ctx->GD_NG)) 
            goto GDDFin;

        for (i = 0; i < ctx->GD_NG; ++i)
            ctx->GD_EV[i] = ctx->PMVIdx[i];

        ctx->GD_TYP = ctx->PMOPT;
        if (gdd_setne(ctx, ctx->GD_TYP))
            goto GDDFin;

        ctx->GD_NOC = ctx->NOC;
    }
    else {                          /* ## list of matrices */
        
        printf1(ctx, "full adjacency matrix (matrices).\n");
        if (ctx->PMNVTyp != 1) {
            printf1(ctx, "Invalid list of matrix names.\n");
            goto GDDFin;
        }

        ctx->GD_NG = 0;         
        for (i = 0; i < ctx->PMNV; ++i) {
            j = ctx->PMVIdx[i];

            if (j < 0 || j >= ctx->MaxMat || ctx->MatAlloc[j] == 0) {
                printf1(ctx, "Error: undefined matrix on right-hand side.\n");
                goto GDDFin;
            }
            if (ctx->MatRow[j] != ctx->MatCol[j]) {
                printf1(ctx, "Error: can only use square matrices.\n");
                goto GDDFin;
            }
            if (i == 0)  
                n = ctx->MatRow[j];
            else if (ctx->MatRow[j] != n) {
                printf1(ctx, "Error: all matrices must have the same order.\n");
                goto GDDFin;
            }
            ctx->GD_NG++;
        }
        if (ctx->GD_NG == 0) {
            printf1(ctx, "Error: need at least one matrix.\n");
            goto GDDFin;
        }
        ctx->GD_NDMAX = ctx->GD_NP = n;

        if (gdd_alloc(ctx, ctx->GD_NG)) 
            goto GDDFin;

        for (i = 0; i < ctx->GD_NG; ++i)
            ctx->GD_EV[i] = ctx->PMVIdx[i];

        ctx->GD_TYP = ctx->PMOPT;
        if (gdd_setne(ctx, ctx->GD_TYP))
            goto GDDFin;

        ctx->GD_NOC = 0;
    }
    ctx->GD_TYP = ctx->PMOPT;

    printf1(ctx, "\nMinimal level for valid edges: %g\n",ctx->GD_EVMin);
    printf1(ctx, "Number of nodes: %d. Largest node number: %d\n",ctx->GD_NP,ctx->GD_NDMAX);
    printf1(ctx, "Graph type %d ",ctx->GD_GT);
    if (ctx->GD_GT == 1)
        printf1(ctx, "(undirected, unvalued).\n");
    else if (ctx->GD_GT == 2) {
        printf1(ctx, "(undirected, valued, bidirectional edges: ");
        if (ctx->GD_GTT == 1)
            printf1(ctx, "minimum).\n");
        else if (ctx->GD_GTT == 2)
            printf1(ctx, "maximum).\n");
        else if (ctx->GD_GTT == 3)
            printf1(ctx, "sum).\n");
    }
    else if (ctx->GD_GT == 3)
        printf1(ctx, "(directed, unvalued).\n");
    else if (ctx->GD_GT == 4)
        printf1(ctx, "(directed, valued).\n");

    printf1(ctx, "\nGraph   edges   loops  max edge value  isolated nodes\n");
    for (i = 0; i < ctx->GD_NG; ++i) { 
        printf1(ctx, "%4d %8d %7d  %14.6f  %12d\n",i + 1,ctx->GD_NE[i],ctx->GD_NL[i],ctx->GD_EVMax[i],ctx->GD_IN[i]);
#ifdef TDA_R_PACKAGE
        {
            double erow[5];
            erow[0] = (double)(i + 1);
            erow[1] = (double)ctx->GD_NE[i];
            erow[2] = (double)ctx->GD_NL[i];
            erow[3] = ctx->GD_EVMax[i];
            erow[4] = (double)ctx->GD_IN[i];
            tda_export_row(ctx, "graph.summary", erow, 5);
        }
#endif
    }
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "graph.summary");
#endif

    printf1(ctx, "\nSuccessfully created new graph data structure. ");
    printf1(ctx, "Current memory: %d bytes.\n",ctx->MemReq);

    if (ctx->PMProtFDef)   
        gdd_prot(ctx);

    if (ctx->PMPERM == 1) {     /* permanently save the input data */

        if (ctx->GD_TYP != 1) {
            printf1(ctx, "Warning: perm parameter requires option 1 (edge list).\n");
            printf1(ctx, "Graph data not permanently saved.\n");
        }
        else {
            if (gdd_palloc(ctx, ctx->NOC,ctx->GD_NG))
                goto GDDFin;

            for (i = 0; i < ctx->NOC; ++i) {
                ctx->GD_EIP[i] = (int)get_data(ctx, ctx->GD_EI,i);
                ctx->GD_EJP[i] = (int)get_data(ctx, ctx->GD_EJ,i);
                k = i * ctx->GD_NG;
                for (j = 0; j < ctx->GD_NG; ++j) {
                    tmp = get_data(ctx, ctx->GD_EV[j],i);
                    if (tmp <= ctx->GD_EVMin)
                        tmp = -1.0;
                    ctx->GD_EP[k++] = (float)tmp;                       
                }
            }
            ctx->GD_PERM = 1;
            printf1(ctx, "Saved data for permanent use. ");
            printf1(ctx, "Current memory: %d bytes.\n",ctx->MemReq);
        }
    }
    err = 0;

GDDFin:
    if (err)  
        gdd_free(ctx, 0);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdd_info()  Info about current graph data                               */

void gdd_info(TDAContext *ctx)
{
    register int i;

    printf1(ctx, "Current relational (graph) data.\n");
    if (ctx->GD_TYP == 0) {
        printf1(ctx, "not defined.\n");
        return;
    }
    printf1(ctx, "Number of nodes: %d. Largest node number: %d\n",ctx->GD_NP,ctx->GD_NDMAX);
    printf1(ctx, "Graph type %d ",ctx->GD_GT);
    if (ctx->GD_GT == 1)
        printf1(ctx, "(undirected, unvalued).\n");
    else if (ctx->GD_GT == 2) {
        printf1(ctx, "(undirected, valued, bidirectional edges: ");
        if (ctx->GD_GTT == 1)
            printf1(ctx, "minimum).\n");
        else if (ctx->GD_GTT == 2)
            printf1(ctx, "maximum).\n");
        else if (ctx->GD_GTT == 3)
            printf1(ctx, "sum).\n");
    }
    else if (ctx->GD_GT == 3)
        printf1(ctx, "(directed, unvalued).\n");
    else if (ctx->GD_GT == 4)
        printf1(ctx, "(directed, valued).\n");

    printf1(ctx, "\nGraph   edges   loops  max edge value  isolated nodes\n");
    for (i = 0; i < ctx->GD_NG; ++i) { 
        printf1(ctx, "%4d %8d %7d  %14.6f  %12d\n",i + 1,ctx->GD_NE[i],ctx->GD_NL[i],ctx->GD_EVMax[i],ctx->GD_IN[i]);
#ifdef TDA_R_PACKAGE
        {
            double erow[5];
            erow[0] = (double)(i + 1);
            erow[1] = (double)ctx->GD_NE[i];
            erow[2] = (double)ctx->GD_NL[i];
            erow[3] = ctx->GD_EVMax[i];
            erow[4] = (double)ctx->GD_IN[i];
            tda_export_row(ctx, "graph.summary", erow, 5);
        }
#endif
    }
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "graph.summary");
#endif
}

/* ------------------------------------------------------------------------ */
/*  gdd_palloc(n,m)  allocate arrays for gdd1 option. n is number of cases  */
/*                   If n = 0 free previously allocated arrays.             */
/*                   m is number of graphs.                                 */
/*                   Return 0 if OK, -1 if error.                           */

int gdd_palloc(TDAContext *ctx, int n,int m)
{
    int err = -1;

    if (n == 0) {
        err = 0;
        goto GDD_PAFin;
    }
    if (n > 0) {
        if (!(ctx->GD_EIP = (int *)calloc((size_t)(n),sizeof(int))))   
            goto GDD_PAFin;
        memrq(ctx, n,sizeof(int));
        ctx->GD_EIPA = n;

        if (!(ctx->GD_EJP = (int *)calloc((size_t)(n),sizeof(int))))   
            goto GDD_PAFin;
        memrq(ctx, n,sizeof(int));
        ctx->GD_EJPA = n;

        if (!(ctx->GD_EP = (float *)calloc((size_t)(n) * (size_t)(m),sizeof(float))))   
            goto GDD_PAFin;
        memrq(ctx, n * m,sizeof(float));
        ctx->GD_EPA = n * m;
        return(0);
    }

GDD_PAFin:
    if (err)  
        printf1(ctx, "Error: insufficient memory for saving graph data.\n");

    if (ctx->GD_EIPA > 0) {
        free((char *)ctx->GD_EIP);
        memrq(ctx, -ctx->GD_EIPA,sizeof(int));
        ctx->GD_EIPA = 0;
    }
    if (ctx->GD_EJPA > 0) {
        free((char *)ctx->GD_EJP);
        memrq(ctx, -ctx->GD_EJPA,sizeof(int));
        ctx->GD_EJPA = 0;
    }
    if (ctx->GD_EPA > 0) {
        free((char *)ctx->GD_EP);
        memrq(ctx, -ctx->GD_EPA,sizeof(float));
        ctx->GD_EPA = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdd_alloc(n) allocate GD_NE, GD_IN, GD_NL, GD_EV, GD_EVMax for n graphs */
/*               if n = 0 free previously allocated memory.                 */
/*               Return 0 if OK, -1 if error.                               */

int gdd_alloc(TDAContext *ctx, int n)  
{
    if (ctx->GD_NEA > 0) {
        free((char *)ctx->GD_NE);
        memrq(ctx, -ctx->GD_NEA,sizeof(int));
        ctx->GD_NEA = 0;
    }
    if (ctx->GD_INA > 0) {
        free((char *)ctx->GD_IN);
        memrq(ctx, -ctx->GD_INA,sizeof(int));
        ctx->GD_INA = 0;
    }
    if (ctx->GD_NLA > 0) {
        free((char *)ctx->GD_NL);
        memrq(ctx, -ctx->GD_NLA,sizeof(int));
        ctx->GD_NLA = 0;
    }
    if (ctx->GD_EVA > 0) {
        free((char *)ctx->GD_EV);
        memrq(ctx, -ctx->GD_EVA,sizeof(int));
        ctx->GD_EVA = 0;
    }
    if (ctx->GD_EVMaxA > 0) {
        free((char *)ctx->GD_EVMax);
        memrq(ctx, -ctx->GD_EVMaxA,sizeof(double));
        ctx->GD_EVMaxA = 0;
    }
    if (n == 0)
        return(0);

    if (!(ctx->GD_NE = (int *)calloc((size_t)(n),sizeof(int)))) { 
        p_err(ctx, -2,1);
        return(-1);
    }
    memrq(ctx, n,sizeof(int));
    ctx->GD_NEA = n;

    if (!(ctx->GD_IN = (int *)calloc((size_t)(n),sizeof(int)))) { 
        p_err(ctx, -2,1);
        return(-1);
    }
    memrq(ctx, n,sizeof(int));
    ctx->GD_INA = n;

    if (!(ctx->GD_NL = (int *)calloc((size_t)(n),sizeof(int)))) { 
        p_err(ctx, -2,1);
        return(-1);
    }
    memrq(ctx, n,sizeof(int));
    ctx->GD_NLA = n;

    if (!(ctx->GD_EV = (int *)calloc((size_t)(n),sizeof(int)))) { 
        p_err(ctx, -2,1);
        return(-1);
    }
    memrq(ctx, n,sizeof(int));
    ctx->GD_EVA = n;

    if (!(ctx->GD_EVMax = (double *)calloc((size_t)(n),sizeof(double)))) { 
        p_err(ctx, -2,1);
        return(-1);
    }
    memrq(ctx, n,sizeof(double));
    ctx->GD_EVMaxA = n;

    return(0);
}

/* ------------------------------------------------------------------------ */
/*  gdd_free(opt)       free graph data structures, if opt !== message      */

void gdd_free(TDAContext *ctx, int opt)
{
    register int i;
    int n;

    gdd_alloc(ctx, 0); 

    if (ctx->GD_NDA > 0) {
        free((char *)ctx->GD_ND);
        memrq(ctx, -ctx->GD_NDA,sizeof(int));
        ctx->GD_NDA = 0;
    }
    if (ctx->GD_APtrA > 0) {
        free((char *)ctx->GD_APtr);
        memrq(ctx, -ctx->GD_APtrA,sizeof(int));
        ctx->GD_APtrA = 0;
    }
    if (ctx->GD_FPNA > 0) {
        for (i = 0; i < ctx->GD_NP; ++i) {
            n = ctx->GD_FPN[i];
            if (n > 0) {
                if (ctx->GD_FPIA > 0) {
                    free((char *)ctx->GD_FPI[i]);
                    memrq(ctx, -n,sizeof(int));
                }
                if (ctx->GD_FPKA > 0) {
                    free((char *)ctx->GD_FPK[i]);
                    memrq(ctx, -n,sizeof(int));
                }
                if (ctx->GD_FPK1A > 0) {
                    free((char *)ctx->GD_FPK1[i]);
                    memrq(ctx, -n,sizeof(int));
                }
            }
        }
        free((char *)ctx->GD_FPN);
        memrq(ctx, -ctx->GD_FPNA,sizeof(int));
        ctx->GD_FPNA = 0;

        if (ctx->GD_FPIA > 0) {
            free((char *)ctx->GD_FPI);
            memrq(ctx, -ctx->GD_FPIA,sizeof(int *));
            ctx->GD_FPIA = 0;
        }
        if (ctx->GD_FPKA > 0) {
            free((char *)ctx->GD_FPK);
            memrq(ctx, -ctx->GD_FPKA,sizeof(int *));
            ctx->GD_FPKA = 0;
        }
        if (ctx->GD_FPK1A > 0) {
            free((char *)ctx->GD_FPK1);
            memrq(ctx, -ctx->GD_FPK1A,sizeof(int *));
            ctx->GD_FPK1A = 0;
        }
    }
    if (ctx->GD_BPNA > 0) {
        for (i = 0; i < ctx->GD_NP; ++i) {
            n = ctx->GD_BPN[i];
            if (n > 0) {
                if (ctx->GD_BPIA > 0) {
                    free((char *)ctx->GD_BPI[i]);
                    memrq(ctx, -n,sizeof(int));
                }
                if (ctx->GD_BPKA > 0) {
                    free((char *)ctx->GD_BPK[i]);
                    memrq(ctx, -n,sizeof(int));
                }
            }
        }
        free((char *)ctx->GD_BPN);
        memrq(ctx, -ctx->GD_BPNA,sizeof(int));
        ctx->GD_BPNA = 0;

        if (ctx->GD_BPIA > 0) {
            free((char *)ctx->GD_BPI);
            memrq(ctx, -ctx->GD_BPIA,sizeof(int *));
            ctx->GD_BPIA = 0;
        }
        if (ctx->GD_BPKA > 0) {
            free((char *)ctx->GD_BPK);
            memrq(ctx, -ctx->GD_BPKA,sizeof(int *));
            ctx->GD_BPKA = 0;
        }
    }
    gdd_palloc(ctx, 0,0);

    if (ctx->GD_TYP && opt)
        printf1(ctx, "Removed graph data.\n");

    ctx->GD_EVMin = 0.0;
    ctx->GD_PERM = ctx->GD_NG = ctx->GD_TYP = 0;
    ctx->GD_EI = ctx->GD_EJ = -1;
}

/* ------------------------------------------------------------------------ */
/*  gdd_setgt(typ,gt,gtt)   Setup GD_GT, GD_GTT for typ.                    */   

void gdd_setgt(TDAContext *ctx, int typ,int gt,int gtt)
{
    ctx->GD_GT = gt;
    ctx->GD_GTT = gtt;

    if (gt < 1 || gt > 4) {
        if (typ == 3 || typ == 4)
            ctx->GD_GT = 2;
        else
            ctx->GD_GT = 4;
    }
    if (typ >= 3 && typ <= 6 && ctx->GD_GT > 2)
        ctx->GD_GT = 2;
}

/* -###-------------------------------------------------------------------- */
/*  gdd_setnd()         Setup nodes list GD_ND, number of nodes GD_NP,      */
/*                      and maximal number of nodes, GD_NDMAX.              */
/*                      This function is called for opt 1 and 2.            */
/*                      Return 0 if OK, -1 if error.                        */

int gdd_setnd(TDAContext *ctx)
{
    register int i,j,k,l;
    int err,ei,ej,ip;
    double tmp;
                       
    err = -1;
    ei = ctx->PMVIdx[0];         /* index of variable for starting node */
    ej = ctx->PMVIdx[1];         /* index of variable for ending node */

    if (alloc_acn(ctx, 2 * ctx->NOC))                 
        goto GDDNDFin;

    ctx->GD_NDMAX = ip = 0;
    for (k = 0; k < ctx->NOC; ++k) {
        i = (int)get_data(ctx, ei,k);
        j = (int)get_data(ctx, ej,k);
        if (i < 1 || j < 1) {
            printf1(ctx, "Error: found negative, or zero, node number in case %d.\n",k + 1);
            goto GDDNDFin;
        }
        ctx->AcN[ip++] = i;
        ctx->AcN[ip++] = j;

        if (ctx->GD_NDMAX < i)
            ctx->GD_NDMAX = i;
        if (ctx->GD_NDMAX < j)
            ctx->GD_NDMAX = j;

        for (l = 0; l < ctx->GD_NG; ++l) {
            if ((tmp = get_data(ctx, ctx->GD_EV[l],k)) >= ctx->GD_EVMin) {
                ctx->GD_NE[l] += 1;
                if (i == j)
                    ctx->GD_NL[l] += 1;

                if (ctx->GD_EVMax[l] < tmp)
                    ctx->GD_EVMax[l] = tmp;
            }
        }
    }
    if (sorti(ctx, 2 * ctx->NOC,ctx->AcN,0))
        goto GDDNDFin;

    k = 1;
    for (i = 1; i < 2 * ctx->NOC; ++i) {
        if (ctx->AcN[i] != ctx->AcN[k - 1])
            ctx->AcN[k++] = ctx->AcN[i]; 
    }
    ctx->GD_NP = k;

    if (!(ctx->GD_ND = (int *)calloc((size_t)(ctx->GD_NP),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto GDDNDFin;
    }
    memrq(ctx, ctx->GD_NP,sizeof(int));
    ctx->GD_NDA = ctx->GD_NP;       
    for (k = 0; k < ctx->GD_NP; ++k)  
        ctx->GD_ND[k] = ctx->AcN[k];

    /* calculate number of isolated nodes */

    for (l = 0; l < ctx->GD_NG; ++l) {
           
        if (alloc_acn(ctx, ctx->GD_NDMAX + 1))                 
            goto GDDNDFin;

        for (k = 0; k < ctx->NOC; ++k) {
            if (get_data(ctx, ctx->GD_EV[l],k) >= ctx->GD_EVMin) {
                i = (int)get_data(ctx, ei,k);
                j = (int)get_data(ctx, ej,k);
                if (i != j)
                    ctx->AcN[i] = ctx->AcN[j] = 1;
            }
        }
        ctx->GD_IN[l] = 0;
        for (k = 0; k < ctx->GD_NP; ++k) {
            if (ctx->AcN[ctx->GD_ND[k]] == 0)
                ctx->GD_IN[l] += 1;
        }
    }
    err = 0;

GDDNDFin:
    alloc_acn(ctx, 0);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  gdd_setptr(opt)     Setup pointer for gdd option 1.                     */
/*                      If opt=0 directed, opt != 0 undirected graph.       */
/*                      Return 0 if OK, -1 if error.                        */

int gdd_setptr(TDAContext *ctx, int opt)
{
    register int i,j,k,l;
    int err,n,m,ei,ej,ii,jj,ki,kj,n1,n2,nn,i1,i2,k1,k2;
    short vidx[2];


    err = -1;
    ei = ctx->PMVIdx[0];         /* index of variable for starting node */
    ej = ctx->PMVIdx[1];         /* index of variable for ending node */
    ctx->GD_FPMax = ctx->GD_BPMax = 0;

    /* forward pointers */

    if (!(ctx->GD_FPN = (int *)calloc((size_t)(ctx->GD_NP),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto GDDPFin;
    }
    memrq(ctx, ctx->GD_NP,sizeof(int));
    ctx->GD_FPNA = ctx->GD_NP;       

    if (!(ctx->GD_FPI = (int **)calloc((size_t)(ctx->GD_NP),sizeof(int *)))) { 
        p_err(ctx, -2,1);
        goto GDDPFin;
    }
    memrq(ctx, ctx->GD_NP,sizeof(int *));
    ctx->GD_FPIA = ctx->GD_NP;       

    if (!(ctx->GD_FPK = (int **)calloc((size_t)(ctx->GD_NP),sizeof(int *)))) { 
        p_err(ctx, -2,1);
        goto GDDPFin;
    }
    memrq(ctx, ctx->GD_NP,sizeof(int *));
    ctx->GD_FPKA = ctx->GD_NP;       

    vidx[0] = (short)(ei);
    vidx[1] = (short)(ej);
    if (vsort(ctx, 2,vidx,1,0,0))
        goto GDDPFin;
        

    k = ctx->VSORTPtr[0];
    i = (int)get_data(ctx, ei,k);
    j = (int)get_data(ctx, ej,k);
    for (l = 1; l < ctx->NOC; ++l) {
        k = ctx->VSORTPtr[l];
        ki = (int)get_data(ctx, ei,k);
        kj = (int)get_data(ctx, ej,k);
        if (ki == i && kj == j) {
            printf1(ctx, "Error: edge %d - %d occurs two times.\n",i,j);
            goto GDDPFin;
        }
        i = ki;
        j = kj;
    }

    if (alloc_ack(ctx, ctx->NOC))
        goto GDDPFin;
    if (alloc_aci(ctx, ctx->NOC))
        goto GDDPFin;
   
    l = 0;
    k = ctx->VSORTPtr[0];
    ki = (int)get_data(ctx, ei,k);
    for (i = 0; i < ctx->GD_NP; ++i) {
        ii = ctx->GD_ND[i];
        n = 0;
        while (ki < ii && ++l < ctx->NOC) {
            k = ctx->VSORTPtr[l];
            ki = (int)get_data(ctx, ei,k);
        }
        if (ki == ii) {
            kj = (int)get_data(ctx, ej,k);
            ctx->AcI[n] = gdd_fndni(ctx, kj);
            ctx->AcK[n++] = k;
            while (++l < ctx->NOC) {
                k = ctx->VSORTPtr[l];
                ki = (int)get_data(ctx, ei,k);
                if (ki == ii) {
                    kj = (int)get_data(ctx, ej,k);
                    ctx->AcI[n] = gdd_fndni(ctx, kj);
                    ctx->AcK[n++] = k;
                }
                else
                    break;
            }
            if (!(ctx->GD_FPI[i] = (int *)calloc((size_t)(n),sizeof(int)))) { 
                p_err(ctx, -2,1);
                goto GDDPFin;
            }
            if (!(ctx->GD_FPK[i] = (int *)calloc((size_t)(n),sizeof(int)))) { 
                free((char *)ctx->GD_FPI[i]);
                p_err(ctx, -2,1);
                goto GDDPFin;
            }
            memrq(ctx, 2 * n,sizeof(int));
            ctx->GD_FPN[i] = n;
            for (j = 0; j < n; ++j) {
                ctx->GD_FPI[i][j] = ctx->AcI[j];
                ctx->GD_FPK[i][j] = ctx->AcK[j];
            }
            if (ctx->GD_FPMax < n)
                ctx->GD_FPMax = n;
        }
    }
    vsort(ctx, 2,vidx,0,0,0);      /* free */

    /* backward pointers */

    if (opt == 0) {

        if (!(ctx->GD_BPN = (int *)calloc((size_t)(ctx->GD_NP),sizeof(int)))) { 
            p_err(ctx, -2,1);
            goto GDDPFin;
        }
        memrq(ctx, ctx->GD_NP,sizeof(int));
        ctx->GD_BPNA = ctx->GD_NP;       

        if (!(ctx->GD_BPI = (int **)calloc((size_t)(ctx->GD_NP),sizeof(int *)))) { 
            p_err(ctx, -2,1);
            goto GDDPFin;
        }
        memrq(ctx, ctx->GD_NP,sizeof(int *));
        ctx->GD_BPIA = ctx->GD_NP;       

        if (!(ctx->GD_BPK = (int **)calloc((size_t)(ctx->GD_NP),sizeof(int *)))) { 
            p_err(ctx, -2,1);
            goto GDDPFin;
        }
        memrq(ctx, ctx->GD_NP,sizeof(int *));
        ctx->GD_BPKA = ctx->GD_NP;       
    }
    else {
        if (!(ctx->GD_FPK1 = (int **)calloc((size_t)(ctx->GD_NP),sizeof(int *)))) { 
            p_err(ctx, -2,1);
            goto GDDPFin;
        }
        memrq(ctx, ctx->GD_NP,sizeof(int *));
        ctx->GD_FPK1A = ctx->GD_NP;       

        if (alloc_acn(ctx, ctx->NOC))
            goto GDDPFin;
        if (alloc_acr(ctx, ctx->NOC))
            goto GDDPFin;
        if (alloc_acs(ctx, ctx->NOC))
            goto GDDPFin;
    }
    vidx[0] = (short)(ej);
    vidx[1] = (short)(ei);
    if (vsort(ctx, 2,vidx,1,0,0))
        goto GDDPFin;
    
    l = 0;
    k = ctx->VSORTPtr[0];
    kj = (int)get_data(ctx, ej,k);
    for (j = 0; j < ctx->GD_NP; ++j) {
        jj = ctx->GD_ND[j];
        n = 0;
        while (kj < jj && ++l < ctx->NOC) {
            k = ctx->VSORTPtr[l];
            kj = (int)get_data(ctx, ej,k);
        }
        if (kj == jj) {
            ki = (int)get_data(ctx, ei,k);
            ctx->AcI[n] = gdd_fndni(ctx, ki);
            ctx->AcK[n++] = k;
            while (++l < ctx->NOC) {
                k = ctx->VSORTPtr[l];
                kj = (int)get_data(ctx, ej,k);
                if (kj == jj) {
                    ki = (int)get_data(ctx, ei,k);
                    ctx->AcI[n] = gdd_fndni(ctx, ki);
                    ctx->AcK[n++] = k;
                }
                else
                    break;
            }
            if (opt == 0) {

                if (!(ctx->GD_BPI[j] = (int *)calloc((size_t)(n),sizeof(int)))) { 
                    p_err(ctx, -2,1);
                    goto GDDPFin;
                }
                if (!(ctx->GD_BPK[j] = (int *)calloc((size_t)(n),sizeof(int)))) { 
                    free((char *)ctx->GD_BPI[j]);
                    p_err(ctx, -2,1);
                    goto GDDPFin;
                }
                memrq(ctx, 2 * n,sizeof(int));
                ctx->GD_BPN[j] = n;
                for (i = 0; i < n; ++i) {
                    ctx->GD_BPI[j][i] = ctx->AcI[i];
                    ctx->GD_BPK[j][i] = ctx->AcK[i];
                }
                if (ctx->GD_BPMax < n)
                    ctx->GD_BPMax = n;
            }
        }
/*##*/  if (opt != 0) {         /* undirected graph */
            m = ctx->GD_FPN[j];

            if (n > 0 || m > 0) {
      
                i1 = i2 = nn = 0;
                while (1) {
                    if (i1 < m) {
                        n1 = ctx->GD_FPI[j][i1];
                        k1 = ctx->GD_FPK[j][i1++];
                    }
                    else
                        k1 = n1 = -1;

                    if (i2 < n) {
                        n2 = ctx->AcI[i2];
                        k2 = ctx->AcK[i2++];
                    }
                    else
                        k2 = n2 = -1;

                    if (n1 < 0 && n2 < 0)
                        break;
                    else if (n1 < 0) {
                        ctx->AcN[nn] = n2;
                        ctx->AcR[nn] = -1;
                        ctx->AcS[nn++] = k2;
                    }
                    else if (n2 < 0) {
                        ctx->AcN[nn] = n1;
                        ctx->AcR[nn] = k1;
                        ctx->AcS[nn++] = -1;
                    }
                    else if (n1 < n2) {
                        ctx->AcN[nn] = n1;
                        ctx->AcR[nn] = k1;
                        ctx->AcS[nn++] = -1;

                        while (i1 < m) {
                            n1 = ctx->GD_FPI[j][i1];              
                            k1 = ctx->GD_FPK[j][i1];              
                            if (n1 >= n2) {
                                if (n1 == n2) {
                                    ctx->AcN[nn] = n1;
                                    ctx->AcR[nn] = k1;
                                    ctx->AcS[nn++] = k2;
                                    i1++;
                                }
                                break;
                            }
                            ctx->AcN[nn] = n1;
                            ctx->AcR[nn] = k1;
                            ctx->AcS[nn++] = -1;
                            i1++;
                        }
                        if (n2 < n1 || i1 >= m) {
                            ctx->AcN[nn] = n2;
                            ctx->AcR[nn] = -1;
                            ctx->AcS[nn++] = k2;
                        }
                    }
                    else if (n2 < n1) {
                        ctx->AcN[nn] = n2;
                        ctx->AcR[nn] = -1;
                        ctx->AcS[nn++] = k2;

                        while (i2 < n) {
                            n2 = ctx->AcI[i2];           
                            k2 = ctx->AcK[i2];           
                            if (n2 >= n1) {
                                if (n2 == n1) {
                                    ctx->AcN[nn] = n2;
                                    ctx->AcR[nn] = k1;
                                    ctx->AcS[nn++] = k2;
                                    i2++;
                                }
                                break;
                            }
                            ctx->AcN[nn] = n2;
                            ctx->AcS[nn] = k2;
                            ctx->AcR[nn++] = -1;
                            i2++;
                        }
                        if (n1 < n2 || i2 >= n) {
                            ctx->AcN[nn] = n1;
                            ctx->AcR[nn] = k1;
                            ctx->AcS[nn++] = -1;
                        }
                    }
                    else {
                        ctx->AcN[nn] = n1;
                        ctx->AcR[nn] = k1;
                        ctx->AcS[nn++] = k2;
                    }
                }
                if (m > 0) {
                    free((char *)ctx->GD_FPI[j]);
                    free((char *)ctx->GD_FPK[j]);
                    memrq(ctx, -2 * m,sizeof(int));
                }
                if (!(ctx->GD_FPI[j] = (int *)calloc((size_t)(nn),sizeof(int)))) { 
                    p_err(ctx, -2,1);
                    goto GDDPFin;
                }
                if (!(ctx->GD_FPK[j] = (int *)calloc((size_t)(nn),sizeof(int)))) { 
                    free((char *)ctx->GD_FPI[j]);
                    p_err(ctx, -2,1);
                    goto GDDPFin;
                }
                if (!(ctx->GD_FPK1[j] = (int *)calloc((size_t)(nn),sizeof(int)))) { 
                    free((char *)ctx->GD_FPI[j]);
                    free((char *)ctx->GD_FPK[j]);
                    p_err(ctx, -2,1);
                    goto GDDPFin;
                }
                memrq(ctx, 3 * nn,sizeof(int));
                ctx->GD_FPN[j] = nn;
                for (i = 0; i < nn; ++i) {
                    ctx->GD_FPI[j][i] = ctx->AcN[i];
                    ctx->GD_FPK[j][i] = ctx->AcR[i];
                    ctx->GD_FPK1[j][i] = ctx->AcS[i];
                }
                if (ctx->GD_FPMax < nn)
                    ctx->GD_FPMax = nn;
            }
        }
    }
    err = 0;

GDDPFin:
    alloc_acs(ctx, 0);   
    alloc_acr(ctx, 0);   
    alloc_acn(ctx, 0);   
    alloc_aci(ctx, 0);   
    alloc_ack(ctx, 0);   
    vsort(ctx, 2,vidx,0,0,0);      /* free */
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdd_setap()         Setup GD_APtr for gdd option 2.                     */
/*                      Return 0 if OK, -1 if error.                        */

int gdd_setap(TDAContext *ctx)
{
    register int i,j,k,ii,jj;
    int err,ei,ej;


    err = -1;

    ei = ctx->PMVIdx[0];         /* index of variable for starting node */
    ej = ctx->PMVIdx[1];         /* index of variable for ending node */

    if (!(ctx->GD_APtr = (int *)calloc((size_t)(ctx->GD_NP) * (size_t)(ctx->GD_NP),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto GDDAPFin;
    }
    ctx->GD_APtrA = ctx->GD_NP * ctx->GD_NP;
    memrq(ctx, ctx->GD_APtrA,sizeof(int));

    for (k = 0; k < ctx->NOC; ++k) {
        i = (int)get_data(ctx, ei,k);
        j = (int)get_data(ctx, ej,k);
        ii = gdd_fndni(ctx, i);
        jj = gdd_fndni(ctx, j);
        if (ii < 0 || jj < 0) {
            printfe(ctx, "ERROR in gdd [i=%d j=%d ii=%d jj=%d]\n",i,j,ii,jj);
            gerr_exit(ctx, 215);
        }
        if (ctx->GD_APtr[ii * ctx->GD_NP + jj] > 0) {
            printf1(ctx, "Error: edge (%d,%d) occurs twice.\n",i,j);
            goto GDDAPFin;
        }
        ctx->GD_APtr[ii * ctx->GD_NP + jj] = k + 1;
    }
    err = 0;

GDDAPFin:
    return(err);
}

/* -###-------------------------------------------------------------------- */
/*  gdd_setne(opt)      Setup GD_EVMax, GD_NE, GD_NL for otions 3 - 7 or 8. */
/*                      Return 0 if OK, -1 if error.                        */

int gdd_setne(TDAContext *ctx, int opt)
{
    register int i,j,k;
    int ne,nl;
    double tmp,emax;

    for (k = 1; k <= ctx->GD_NG; ++k) {

        if (alloc_acn(ctx, ctx->GD_NDMAX + 1))                 
            return(-1);        

        ne = nl = 0;
        emax = 0.0;
        for (i = 0; i < ctx->GD_NP; ++i) {
            for (j = 0; j < ctx->GD_NP; ++j) {
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
                tmp = gdd_adj(ctx, i,j,k);
                if (tmp >= ctx->GD_EVMin) {
                    ne++;
                    if (i == j)
                        nl++;
                    emax = dmax(ctx, emax,tmp);
                    if (i != j)
                        ctx->AcN[i] = ctx->AcN[j] = 1;
                }
            }
        }
        ctx->GD_NE[k - 1] = ne;
        ctx->GD_NL[k - 1] = nl;
        ctx->GD_EVMax[k - 1] = emax;

        ctx->GD_IN[k - 1] = 0;                   /* number of isolated nodes */
        for (i = 0; i < ctx->GD_NP; ++i) {
            if (ctx->AcN[i] == 0)
                ctx->GD_IN[k - 1] += 1;
        }
    }
    alloc_acn(ctx, 0);                
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  gdd_prot()  print current gdd data structure to protocol file.          */

void gdd_prot(TDAContext *ctx)
{     
    register int i,j,n;

    fprintf(ctx->PMProtFd,"Current gdd data structure\n\n");
    fprintf(ctx->PMProtFd,"GD_TYP   = %d\n",ctx->GD_TYP);
    fprintf(ctx->PMProtFd,"GD_GT    = %d\n",ctx->GD_GT);
    fprintf(ctx->PMProtFd,"GD_GTT   = %d\n",ctx->GD_GTT);
    fprintf(ctx->PMProtFd,"GD_NOC   = %d\n",ctx->GD_NOC);
    fprintf(ctx->PMProtFd,"GD_NG    = %d\n",ctx->GD_NG);
    fprintf(ctx->PMProtFd,"GD_NP    = %d\n",ctx->GD_NP);
    fprintf(ctx->PMProtFd,"GD_NDMAX = %d\n",ctx->GD_NDMAX);

    fprintf(ctx->PMProtFd,"\nGD_EV\n");
    for (i = 0; i < ctx->GD_NG; ++i)
        fprintf(ctx->PMProtFd,"%5d %5d\n",i,ctx->GD_EV[i]);

    fprintf(ctx->PMProtFd,"\nGD_EVMax\n");
    for (i = 0; i < ctx->GD_NG; ++i)
        fprintf(ctx->PMProtFd,"%5d     %g\n",i,ctx->GD_EVMax[i]);

    fprintf(ctx->PMProtFd,"\nGD_NE\n");
    for (i = 0; i < ctx->GD_NG; ++i)
        fprintf(ctx->PMProtFd,"%5d %5d\n",i,ctx->GD_NE[i]);

    fprintf(ctx->PMProtFd,"\nGD_NL\n");
    for (i = 0; i < ctx->GD_NG; ++i)
        fprintf(ctx->PMProtFd,"%5d %5d\n",i,ctx->GD_NL[i]);

    fprintf(ctx->PMProtFd,"\nGD_IN\n");
    for (i = 0; i < ctx->GD_NG; ++i)
        fprintf(ctx->PMProtFd,"%5d %5d\n",i,ctx->GD_IN[i]);

    if (ctx->GD_TYP <= 2) {
        fprintf(ctx->PMProtFd,"\nGD_ND\n");
        for (i = 0; i < ctx->GD_NP; ++i) {
            fprintf(ctx->PMProtFd,"%5d ",i);
            fprintf(ctx->PMProtFd,"%5d ",ctx->GD_ND[i]);
            fprintf(ctx->PMProtFd,"\n");
        }
    }
    if (ctx->GD_TYP == 1) {
        fprintf(ctx->PMProtFd,"\nForward pointer (GD_FPMax=%d)\n",ctx->GD_FPMax);
        for (i = 0; i < ctx->GD_NP; ++i) {
            n = ctx->GD_FPN[i];
            fprintf(ctx->PMProtFd,"%5d %5d : ",i,n);
            for (j = 0; j < n; ++j) {
                fprintf(ctx->PMProtFd,"%5d [%5d",ctx->GD_FPI[i][j],ctx->GD_FPK[i][j]);
                if (ctx->GD_GT <= 2)
                    fprintf(ctx->PMProtFd,",%5d] ",ctx->GD_FPK1[i][j]);
                else
                    fprintf(ctx->PMProtFd,"] ");
            }
            fprintf(ctx->PMProtFd,"\n");
        }
        if (ctx->GD_GT > 2) {
            fprintf(ctx->PMProtFd,"\nBackward pointer (GD_BPMax=%d)\n",ctx->GD_BPMax);
            for (i = 0; i < ctx->GD_NP; ++i) {
                n = ctx->GD_BPN[i];
                fprintf(ctx->PMProtFd,"%5d %5d : ",i,n);
                for (j = 0; j < n; ++j)
                    fprintf(ctx->PMProtFd,"%5d [%5d] ",ctx->GD_BPI[i][j],ctx->GD_BPK[i][j]);
                fprintf(ctx->PMProtFd,"\n");
            }
        }
    }
    if (ctx->GD_TYP == 2) {
        fprintf(ctx->PMProtFd,"\nGD_APtr\n");
        for (i = 0; i < ctx->GD_NP; ++i) {
            fprintf(ctx->PMProtFd,"%5d : ",i);
            for (j = 0; j < ctx->GD_NP; ++j)  
                fprintf(ctx->PMProtFd,"%4d ",ctx->GD_APtr[i * ctx->GD_NP + j]);
            fprintf(ctx->PMProtFd,"\n");
        }
    }
    fprintf(ctx->PMProtFd,"\n");
}

/* ------------------------------------------------------------------------ */
/*  gdd_fndni(k)    Find internal node number for external node number k.   */
/*                  return -1 if not found.                                 */

int gdd_fndni(TDAContext *ctx, int k)
{
    register int i,j,ij = 0;

    if (ctx->GD_TYP > 2) {
        if (k < 1 || k > ctx->GD_NP)
            return(-1);
        return(k - 1);
    }
    i = 0;
    j = ctx->GD_NP - 1;
    while (i <= j) {
        ij = (i + j) / 2;
        if (k < ctx->GD_ND[ij])
            j = ij - 1;
        else if (k > ctx->GD_ND[ij])
            i = ij + 1;
        else  
            break;
    }
    if (k != ctx->GD_ND[ij])
        return(-1);
    return(ij);
}

/* ------------------------------------------------------------------------ */
/*  check_gdj(j)    Return 1 if variable j is required for current gdd      */
/*                  data structure, otherwise return 0.                     */

int check_gdj(TDAContext *ctx, int j)
{
    register int i;

    if (ctx->GD_TYP == 0 || ctx->GD_PERM != 0)
        return(0);
    if (j == ctx->GD_EI || j == ctx->GD_EJ)
        return(1);
    for (i = 0; i < ctx->GD_NG; ++i) {
        if (j == ctx->GD_EV[i])
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

double gdd_adj(TDAContext *ctx, int i,int j,int gn)
{
    register int k,l,n;
    double tmp1,tmp2;

    if (ctx->GD_TYP == 0)
        return(-1.0);

    if (i < 0 || j < 0 || i >= ctx->GD_NP || j >= ctx->GD_NP || gn < 1 || gn > ctx->GD_NG)
        return(-1.0);

    tmp1 = tmp2 = -1.0;

    switch (ctx->GD_TYP) {
        case 1:     n = ctx->GD_FPN[i];
                    for (k = 0; k < n; ++k) {
                        if (j == ctx->GD_FPI[i][k]) {
                            l = ctx->GD_FPK[i][k];
                            if (l >= 0)
                                tmp1 = gdd_ev(ctx, l,gn);

                            if (ctx->GD_GT <= 2) {
                                if (i == j)
                                    tmp2 = tmp1;
                                else  
                                    l = ctx->GD_FPK1[i][k];
                                if (l >= 0)
                                    tmp2 = gdd_ev(ctx, l,gn);
                            }
                            break;
                        }
                    }
                    break;

        case 2:     l = ctx->GD_APtr[i * ctx->GD_NP + j];
                    if (l > 0)
                        tmp1 = gdd_ev(ctx, l - 1,gn);

                    if (ctx->GD_GT <= 2) {
                        if (i == j)
                            tmp2 = tmp1;
                        else {
                            l = ctx->GD_APtr[j * ctx->GD_NP + i];
                            if (l > 0)
                                tmp2 = gdd_ev(ctx, l - 1,gn);
                        }
                    }
                    break;

        case 3:     if (j != i) {
                        if (j < i)
                            l = i * (i - 1) / 2 + j;
                        else               
                            l = j * (j - 1) / 2 + i;

                        tmp1 = tmp2 = gdd_ev(ctx, l,gn);
                    }
                    break;

        case 4:     if (j != i) {
                        if (i < j)
                            l = i * ctx->GD_NP - i * (i + 1) / 2 + j - i - 1;
                        else
                            l = j * ctx->GD_NP - j * (j + 1) / 2 + i - j - 1;

                        tmp1 = tmp2 = gdd_ev(ctx, l,gn);
                    }
                    break;


        case 5:     if (j <= i)  
                        l = i * (i + 1) / 2 + j;
                    else
                        l = j * (j + 1) / 2 + i;
                    tmp1 = tmp2 = gdd_ev(ctx, l,gn);
                    break;

        case 6:     if (j >= i)  
                        l = i * ctx->GD_NP - i * (i - 1) / 2 + j - i;
                    else
                        l = j * ctx->GD_NP - j * (j - 1) / 2 + i - j;
                    tmp1 = tmp2 = gdd_ev(ctx, l,gn);
                    break;

        case 7:     l = i * ctx->GD_NP + j;
                    tmp1 = gdd_ev(ctx, l,gn);

                    if (ctx->GD_GT <= 2) {
                        l = j * ctx->GD_NP + i;
                        tmp2 = gdd_ev(ctx, l,gn);
                    }
                    break;
                
        case 8:     l = i * ctx->GD_NP + j + 1;
                    tmp1 = gdd_ev(ctx, l,gn);

                    if (ctx->GD_GT <= 2) {
                        l = j * ctx->GD_NP + i + 1;
                        tmp2 = gdd_ev(ctx, l,gn);
                    }
                    break;

    }

    if (ctx->GD_GT == 1) {
        if (tmp1 >= 0.0 || tmp2 >= 0.0)
            return(1.0);
        else
            return(-1.0);
    }
    else if (ctx->GD_GT == 2) {
        if (tmp1 < 0.0)
            return(tmp2);
        else if (tmp2 < 0.0)
            return(tmp1);
        else if (ctx->GD_GTT == 1)
            return(dmin(ctx, tmp1,tmp2));
        else if (ctx->GD_GTT == 2)
            return(dmax(ctx, tmp1,tmp2));
        else if (i == j)
            return(tmp1);
        else
            return(tmp1 + tmp2);
    }
    else if (ctx->GD_GT == 3) {
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

int g_suc(TDAContext *ctx, int i,int gn)
{
    register int j,dp;
    int nn,n;

    if (ctx->GD_TYP == 0)
        return(0);
           
    if (i < 0 || i >= ctx->GD_NP || gn < 1 || gn > ctx->GD_NG)
        return(0);

    nn = 0;
    if (ctx->GD_TYP == 1) {
        n = ctx->GD_FPN[i];
        for (j = 0; j < n; ++j) {
            if (ctx->GD_FPI[i][j] == i)   
                continue;
     
            dp = ctx->GD_FPK[i][j];
            if (dp >= 0 && gdd_ev(ctx, dp,gn) >= 0.0)  
                    nn++;
            else if (ctx->GD_GT <= 2) {
                dp = ctx->GD_FPK1[i][j];
                if (dp >= 0 && gdd_ev(ctx, dp,gn) >= 0.0)
                    nn++;
            }
        }
    }
    else {
        for (j = 0; j < ctx->GD_NP; ++j) {
            if (i != j && gdd_adj(ctx, i,j,gn) >= 0.0)
                nn++;
        }
    }
    return(nn);
}

/* ------------------------------------------------------------------------ */
/*  g_pre(i)    return number of predecessors of node i, graph gn           */

int g_pre(TDAContext *ctx, int i,int gn)
{
    register int j,dp;
    int nn,n;

    if (ctx->GD_TYP == 0)
        return(0);

    if (i < 0 || i >= ctx->GD_NP || gn < 1 || gn > ctx->GD_NG)
        return(0);

    nn = 0;
    if (ctx->GD_TYP == 1) {
        if (ctx->GD_GT <= 2)
            nn = g_suc(ctx, i,gn);
        else {
            n = ctx->GD_BPN[i];
            for (j = 0; j < n; ++j) {
                if (ctx->GD_BPI[i][j] == i)
                    continue;
                dp = ctx->GD_BPK[i][j];
                if (dp >= 0 && gdd_ev(ctx, dp,gn) >= 0.0)
                    nn++;
            }
        }
    }
    else {
        for (j = 0; j < ctx->GD_NP; ++j) {
            if (j != i && gdd_adj(ctx, j,i,gn) >= 0.0)
                nn++;
        }
    }
    return(nn);
}

/* ------------------------------------------------------------------------ */
/*  g_loop(i)   return 1 if node i has a loop, graph number gn.             */

int g_loop(TDAContext *ctx, int i,int gn)
{
    if (gdd_adj(ctx, i,i,gn) >= 0.0)
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  g_getnd(n,nd)   Get list of input nodes from PMF2d into nd[], max       */
/*                  number of nodes is n. If n = 0 expect an array          */
/*                  nd[0,...,GDU_NP-1] and sets flags.                      */
/*                                                                          */
/*                  Return:  -1 if error, otherwise number of nodes         */

int g_getnd(TDAContext *ctx, int n,int *nd)
{
    int nn,k;
    char buf[101],*p;

    if (ctx->PMF2Def == 0) {
        printf1(ctx, "Error: need if parameter for list of nodes.\n");
        return(-1);
    }
    nn = 0;
    while (fgets(buf,100,ctx->PMF2d)) {
        p = skip_b(ctx, buf);
        if (sscanf(p,"%d",&k) == 1 && k >= 1) {

            if (ctx->GD_TYP <= 2) {
                if (k <= ctx->GD_NDMAX) {                         
                    k = gdd_fndni(ctx, k);
                    if (k < 0)
                        continue;
                }
            }
            else {      
                if (k > ctx->GD_NP)
                    continue;
                k--;
            }

            if (n > 0) {
                if (nn >= n) {
                    printf1(ctx, "Error: exceeded maximal number of nodes for input (%d).\n",n);
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

int g_getne(TDAContext *ctx, int n,int *ni,int *nj)
{
    int nn,k1,k2;
    char buf[101],*p;

    if (ctx->PMF2Def == 0) {
        printf1(ctx, "Error: need if parameter for list of edges.\n");
        return(-1);
    }
    nn = 0;
    while (fgets(buf,100,ctx->PMF2d)) {
        p = skip_b(ctx, buf);
        if (sscanf(p,"%d",&k1) != 1 || k1 < 1 || k1 > ctx->GD_NDMAX)
            continue;
        p = skip_int(ctx, p);
        p = skip_b(ctx, p);
        if (sscanf(p,"%d",&k2) != 1 || k2 < 1 || k2 > ctx->GD_NDMAX)
            continue;

        if (nn >= n) {
            printf1(ctx, "Error: exceeded maximal number of edges for input (%d).\n",n);
            return(-1);
        }
        if (ctx->GD_TYP > 2) {
            k1--;
            k2--;
        }
        else {
            k1 = gdd_fndni(ctx, k1);
            if (k1 < 0)
                continue;

            k2 = gdd_fndni(ctx, k2);
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

int gdd_check(TDAContext *ctx, int gn,int opt)
{
    if (ctx->GD_TYP == 0) {
        printf1(ctx, "Error: no graph data defined.\n");
        return(-1);
    }
    if (gn > ctx->GD_NG) {
        printf1(ctx, "Error: graph number %d not defined.\n",gn);
        return(-1);
    }
    ctx->PMGN = gn;  
    if (gn == 0) {
        if (opt == 0 || ctx->GD_NG == 1)
            ctx->PMGN = 1;
    }
    if (ctx->PMGN == 0)
        printf1(ctx, "Using multigraph ");
    else
        printf1(ctx, "Using graph number %d ",ctx->PMGN);

    if (ctx->GD_GT == 1)
        printf1(ctx, "(undirected, unvalued).\n");
    else if (ctx->GD_GT == 2) {
        printf1(ctx, "(undirected, valued, bidirectional edges: ");
        if (ctx->GD_GTT == 1)
            printf1(ctx, "minimum).\n");
        else if (ctx->GD_GTT == 2)
            printf1(ctx, "maximum).\n");
        else if (ctx->GD_GTT == 3)
            printf1(ctx, "sum).\n");
    }
    else if (ctx->GD_GT == 3)
        printf1(ctx, "(directed, unvalued).\n");
    else if (ctx->GD_GT == 4)
        printf1(ctx, "(directed, valued).\n");
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  gdd_check2(gn,gn1,opt)  Check for gdd data. Two graphs, gn and gn1.     */
/*                          If opt = 0 and gn = 0 and gn1 == 0              */
/*                          select gn=1, gn1=2.                             */
/*                          Set PMGN and PMGN1.                             */
/*                          Return 0 if OK, -1 if error.                    */

int gdd_check2(TDAContext *ctx, int gn,int gn1,int opt)
{
    if (ctx->GD_TYP == 0) {
        printf1(ctx, "Error: no graph data defined.\n");
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
    if (gn > ctx->GD_NG || gn1 > ctx->GD_NG) {
        printf1(ctx, "Error: at least one graph number not defined.\n");
        return(-1);
    }
    ctx->PMGN = gn;  
    ctx->PMGN1 = gn1;
    printf1(ctx, "Using graphs %d and %d ",ctx->PMGN,ctx->PMGN1);

    if (ctx->GD_GT == 1)
        printf1(ctx, "(undirected, unvalued).\n");
    else if (ctx->GD_GT == 2) {
        printf1(ctx, "(undirected, valued, bidirectional edges: ");
        if (ctx->GD_GTT == 1)
            printf1(ctx, "minimum).\n");
        else if (ctx->GD_GTT == 2)
            printf1(ctx, "maximum).\n");
        else if (ctx->GD_GTT == 3)
            printf1(ctx, "sum).\n");
    }
    else if (ctx->GD_GT == 3)
        printf1(ctx, "(directed, unvalued).\n");
    else if (ctx->GD_GT == 4)
        printf1(ctx, "(directed, valued).\n");
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

int gdd_tcheck(TDAContext *ctx, int typ,int dir,int val)
{
    if (typ == 1) {
        if (ctx->GD_TYP != 1) {
            printf1(ctx, "Error: need a graph defined with gdd option 1 (edge list).\n");
            return(-1);
        }
    }
    if (dir == 1) {
        if (ctx->GD_GT > 2) {
            printf1(ctx, "Error: graph must be undirected.\n");
            return(-1);
        }
    }
    else if (dir == 2) {
        if (ctx->GD_GT <= 2) {
            printf1(ctx, "Error: graph must be directed.\n");
            return(-1);
        }
    }
    if (val == 1) {
        if (ctx->GD_GT != 1 && ctx->GD_GT != 3) {
            printf1(ctx, "Error: graph must be unvalued.\n");
            return(-1);
        }
    }
    else if (val == 2) {
        if (ctx->GD_GT != 2 && ctx->GD_GT != 4) {
            printf1(ctx, "Error: graph must be valued.\n");
            return(-1);
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  gdd_ei(k)   Return node number in k.th data matrix row.                 */

int gdd_ei(TDAContext *ctx, int k)              
{
    if (ctx->GD_PERM)
        return(ctx->GD_EIP[k]);
    return((int)get_data(ctx, ctx->GD_EI,k));
}

/* ------------------------------------------------------------------------ */
/*  gdd_ej(k)   Return node number in k.th data matrix row.                 */

int gdd_ej(TDAContext *ctx, int k)              
{
    if (ctx->GD_PERM)
        return(ctx->GD_EJP[k]);
    return((int)get_data(ctx, ctx->GD_EJ,k));
}

/* ------------------------------------------------------------------------ */
/*  gdd_ev(k,gn)  Return edge value in k.th data matrix row, graph gn.      */
/*                Graphs are counted 1,...,GD_GN - 1.                       */
               
double gdd_ev(TDAContext *ctx, int k,int gn)       
{
    double tmp;

    if (k < 0)
        return(-1.0);
              
    if (ctx->GD_PERM)
        return((double)ctx->GD_EP[k * ctx->GD_NG + gn - 1]);

    if (ctx->GD_TYP == 8)  
        tmp = ctx->MatVal[ctx->GD_EV[gn - 1]][k];
    else
        tmp = get_data(ctx, ctx->GD_EV[gn - 1],k); 

    if (tmp < ctx->GD_EVMin)
        tmp = -1.0;
    return(tmp); 
}

/* ------------------------------------------------------------------------ */
/*  gdd_node(i)     Return external node number corresponding to i.         */

int gdd_node(TDAContext *ctx, int i)              
{
    if (ctx->GD_TYP == 1 || ctx->GD_TYP == 2)
        return(ctx->GD_ND[i]);
    return(i + 1);
}

/* ------------------------------------------------------------------------ */
/*  n_info      Print info about processed nodes to stderr.                 */
  
void n_info(TDAContext *ctx, int i,int n)
{
    register int j;

    if (i > 0 && ctx->SILENTFlg < 2 && ctx->GD_NP > 500) {
        j = i / n;
        if (j * n == i) {        
            printfe(ctx, "Node: %7d     %c",i,CR);
            fflushe(ctx);
        }
    }
}          
void n_info_e(TDAContext *ctx)
{
    printfe(ctx, "\n");
    fflushe(ctx);
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

int gcd(TDAContext *ctx)
{
    register int i,j;
    int err,nrec,ne,ni;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);



    printf1(ctx, "Create unvalued graph. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto GCDFin;

    if (ctx->PMN < 1)
        ctx->PMN = 10;
    if (ctx->PMM < 0)
        ctx->PMM = 0;
    else if (ctx->PMM > ctx->PMN * ctx->PMN)
        ctx->PMM = ctx->PMN * ctx->PMN;

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);
    if (ctx->PMOPT > 3)
        ctx->PMOPT = 3;

    if (ctx->PMOPT == 3 && ctx->PMM == ctx->PMN * ctx->PMN)  
        ctx->PMOPT = 2;

    printf1(ctx, "Graph with %d nodes, ",ctx->PMN);
    if (ctx->PMOPT == 1)
        printf1(ctx, "complete without loops.\n");
    else if (ctx->PMOPT == 2)
        printf1(ctx, "complete, including loops.\n");
    else  
        printf1(ctx, "%d edges (randomly selected).\n",ctx->PMM);

    ni = nrec = 0;
    if (ctx->PMOPT <= 2) {
        for (i = 1; i <= ctx->PMN; ++i) {
            for (j = 1; j <= ctx->PMN; ++j) {
                if (i == j && ctx->PMOPT ==1)
                    continue;

                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,1);
                fprintf(ctx->PMFd,"\n");
                nrec++;
            }
        }
    }
    else if (ctx->PMM == 0) {
        for (i = 1; i <= ctx->PMN; ++i) {
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,-1);
            fprintf(ctx->PMFd,"\n");
            nrec++;
        }
    }
    else {                          /* random selection */
         if (alloc_aci(ctx, ctx->PMN))
            goto GCDFin;

        ne = 0;
        while (ne < ctx->PMM) {
            j = (int)((double)ctx->PMN * random1(ctx));
            if (j >= ctx->PMN)
                j = ctx->PMN - 1;
            if (ctx->AcI[j] < ctx->PMN) {
                ctx->AcI[j] += 1;
                ne++;
            }

        }
      
        for (i = 0; i < ctx->PMN; ++i) {
            if (ctx->AcI[i] == 0) {
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,-1);
                fprintf(ctx->PMFd,"\n");
                nrec++;
                ni++;
            }
            else if (ctx->AcI[i] == ctx->PMN) {
                for (j = 1; j <= ctx->PMN; ++j) {
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j);
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,1);
                    fprintf(ctx->PMFd,"\n");
                    nrec++;
                }
            }
            else {
                if (alloc_acc(ctx, ctx->PMN))
                    goto GCDFin;
                ne = 0;
                while (ne < ctx->AcI[i]) {
                    j = (int)((double)ctx->PMN * random1(ctx));
                    if (j >= ctx->PMN)
                        j = ctx->PMN - 1;
                    if (ctx->AcC[j] == 0) {
                        ctx->AcC[j] = 1;
                        ne++;
                    }
                }
                for (j = 0; j < ctx->PMN; ++j) {
                    if (ctx->AcC[j]) {
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j + 1);
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,1);
                        fprintf(ctx->PMFd,"\n");
                        nrec++;
                    }
                }
            }
            if (i > 0 && ctx->SILENTFlg < 2 && ctx->PMN > 100) {
                j = i / 100;
                if (j * 100 == i) {        
                    printfe(ctx, "Node: %7d     %c",i,CR);
                    fflushe(ctx);
                }
            }
        }
        if (ctx->SILENTFlg < 2 && ctx->PMN > 100) {
            printfe(ctx, "\n");
            fflushe(ctx);
        }
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    if (ctx->PMOPT == 3 && ni > 0)
        printf1(ctx, "Graph contains %d isolated nodes.\n",ni);    

    err = 0;

GCDFin:       
    p_clean(ctx);
    return(err);
}

