/****************************************************************************/
/*  t_gdf                                                                   */
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
#include "t_con.h"
#include "t_cdf.h"
#include "t_lsei.h"
#include "t_mat.h"
#include "t_sort.h"
#include "t_svd.h"
#include "t_ml.h"
#include "t_gf.h"
#include "t_freq.h"
#include "tda_context.h"

/*  functions in t_gdf.c */

int gdf(TDAContext *ctx);
void gdf_pmet(TDAContext *ctx, int met,int opt);
int gdf_dcheck(TDAContext *ctx);
int alloc_gdfytyp(TDAContext *ctx, int n);
int alloc_gdfgen(TDAContext *ctx, int nu,int ndim);
int gdf_gcheck(TDAContext *ctx);
void gdf_pmtyp(TDAContext *ctx);
int gdf_marg(TDAContext *ctx, int opt);
void gdf_pdat1(TDAContext *ctx, int d,int n,double *y,double *f);
void gdf_pdat1a(TDAContext *ctx, int d,int n,double *y,double *f,short *cen,int opt);
int gdf_edf1(TDAContext *ctx, int n,double *y,double *f,int opt);
int gdf_edf2(TDAContext *ctx, int n,double *y,short *cen,double *y1,double *h,int opt);
int gdf_joint(TDAContext *ctx, int met,int opt);
void gdf_pdat2(TDAContext *ctx, int t,int m,int ndim,double *val,double *f);
void gdf_pdat2a(TDAContext *ctx, int t,int m,int ndim,double *val,int *ptr,int opt);
void gdf_cptr(TDAContext *ctx, int t,int nu,int ndim,int *ptr);
int gdf_domain(TDAContext *ctx, int t,int nu,int ndim,int *ptr,double *yval,int opt,double d);
int gdf_edf1m(TDAContext *ctx, int t,int nu,int ndim,int *ptr,double *val,double *f,int opt);
int gdf_edf3m(TDAContext *ctx, int t,int nu,int ndim,int *ptr,double *yval,double *val, double *f,int opt,int prn,int mxit,double tol);
int gdf_edf4m(TDAContext *ctx, int t,int nu,int ndim,int *ptr,double *yval,double *val, double *f,double d,int opt,int prn);
void gdf_getb(TDAContext *ctx, int n,int ndim,int opt);
int gdf_boxptr(TDAContext *ctx, int ndim,int *ptr);
void gdf_getmval(TDAContext *ctx, int ne,int ndim,int np,int nred,int *acm,double *mval);

int lsreg1(TDAContext *ctx);
void ls1_xx(TDAContext *ctx, int nx,int iflag);
void ls1_gpar(TDAContext *ctx, int n,double *y,int nx,int iflag,double *xx,double *b,double *h);
void ls1_res(TDAContext *ctx, int n,int nx,int iflag,double *b,double *yp,double *res);
void ls1_ppar(TDAContext *ctx, int nx,double *beta);
void ls1_pres(TDAContext *ctx, int nx,int iflag,double *yp,double *res);
int ls1_marg(TDAContext *ctx, int n,double *res,short *ytyp);
int ls1_joint(TDAContext *ctx, int n,double *res,short *ytyp,int opt);




/* ------------------------------------------------------------------------ */
/*  gdf()           General distribution functions                          */
/*                                                                          */
/*                  gdf(                                                    */
/*                      opt=...,        1 (default) marginal calculations   */
/*                                      2 joint calculations, method 1      */
/*                                      2 joint calculations, method 2      */
/*                      prn=...,        0 (default) distribution function   */
/*                                      1 survivor function                 */
/*                                      2 expected values                   */
/*                      yl=...,         variable (lower bound)              */
/*                      cen=...,        censoring indicator                 */
/*                      n=...,          number of boxes, def. 100           */
/*                      sc=...,         offset for domain, def. 0.0         */
/*                      d=...,          proportion of delta, def. 0.1       */
/*                      mxit=...,       max iterations, def. 20             */
/*                      tolf=...,       tolerance for convergence, 0.001    */
/*                      grp=ID,L1,      defines hierarchical structure      */
/*                      fmt = ...,      print format, def. 10.4             */
/*                      prot=...,       additional output file with         */
/*                                      diagnostic information              */
/*                  ) = dfile;                                              */
/*                                                                          */
/*  grp = ID,L1,                                                            */
/*                                                                          */
/*  If this parameter is not used, each data matrix row is treated as a     */
/*  separate unit. Otherwise, ID is an ID variable that defines units,      */
/*  and L1 is a variables that defines the dimension (level 1) to which     */
/*  the corresponding data matrix row belongs. Values of L1 must be         */
/*  positive integers. It is not required that these values are contiguous. */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int gdf(TDAContext *ctx)
{
    int err;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "General distribution functions. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->MxItFlg = 0;
    ctx->TOLF = 0.001; 
    if (parm(ctx, ctx->CmdBuf + 3,1,1))     /* get parameters */
        goto GDFFin;

    if (ctx->PMYH >= 0) {
        printf1(ctx, "Error: yh parameter cannot be used with gdf command.\n");
        goto GDFFin;
    }
    if (ctx->MxItFlg == 0)
        ctx->MxIter = 20;
    else if (ctx->MxIter < 0)
        ctx->MxIter = 0;
    if (ctx->PMOPT > 3)
        ctx->PMOPT = 3;
    if (ctx->PMPRNO > 2)
        ctx->PMPRNO = 2;
    if (ctx->PMD <= 0.0)
        ctx->PMD = 0.1;
    else if (ctx->PMD > 1.0)
        ctx->PMD = 1.0;
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    newline(ctx);
    if (gdf_dcheck(ctx))           /* check dependent variable */
        goto GDFFin;

    if (gdf_gcheck(ctx))           /* check group structure */
        goto GDFFin;

    ctx->GDFWNOC = 0;                /* records written to output file */

    if (ctx->PMOPT == 1) {           /* marginal distributions */
        if (gdf_marg(ctx, ctx->PMPRNO))
            goto GDFFin;
    }
    else {                      /* joint distributions */
        if (gdf_joint(ctx, ctx->PMOPT,ctx->PMPRNO))
            goto GDFFin;
    }
    if (ctx->GDFWNOC > 0) 
        printf1(ctx, "\n%d records written to: %s\n",ctx->GDFWNOC,ctx->PMFdName);
    err = 0;

GDFFin:
    alloc_gdfytyp(ctx, 0);
    alloc_gdfgen(ctx, 0,0);
    if (ctx->PM1NV > 0)                     
        vsort(ctx, 0,ctx->PM1VIdx,0,0,0);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdf_pmet(met,opt)      Print method and option.                         */

void gdf_pmet(TDAContext *ctx, int met,int opt)
{
    if (met == 1)
        printf1(ctx, "Marginal calculation:");
    else if (met == 2)
        printf1(ctx, "Joint calculation (method 1):");
    else                   
        printf1(ctx, "Joint calculation (method 2):");
    if (opt == 0)
        printf1(ctx, " distribution functions.\n\n");
    else if (opt == 1)
        printf1(ctx, " survivor functions.\n\n");
    else
        printf1(ctx, " expected values.\n\n");
}

/* ------------------------------------------------------------------------ */
/*  gdf_dcheck()     Check dependent variables.                             */
/*                                                                          */
/*  Set: GDFYTyp =   1 if exact                                             */
/*                   2 if interval censored                                 */
/*                   3 if right censored                                    */
/*                                                                          */
/*  GDFNEX  = number of exact cases                                         */
/*  GDFNINT = number of interval censored cases                             */
/*  GDFNCEN = number of right censored cases                                */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int gdf_dcheck(TDAContext *ctx)
{
    register int i;
    double yl,yh;    

    ctx->GDFNEX = ctx->GDFNINT = ctx->GDFNCEN = 0;

    if (ctx->PMYL < 0) {
        printf1(ctx, "Error: yl parameter is required.\n");
        return(-1);
    }
    printf1(ctx, "yl: %s",ctx->VName[ctx->PMYL]);
    if (ctx->PMYH >= 0)
        printf1(ctx, "  yh: %s",ctx->VName[ctx->PMYH]);
    if (ctx->PMCEN >= 0)
        printf1(ctx, "  cen: %s",ctx->VName[ctx->PMCEN]);
    newline(ctx);

    if (alloc_gdfytyp(ctx, ctx->NOC))
        return(-1);   

    for (i = 0; i < ctx->NOC; ++i) {
        ctx->GDFYTyp[i] = 1;
        yl = get_data(ctx, ctx->PMYL,i);
        if (ctx->PMYH >= 0) {
            yh = get_data(ctx, ctx->PMYH,i);
            if (yh > yl)
                ctx->GDFYTyp[i] = 2;
            else if (yh < yl) {
                printf1(ctx, "Error: no valid interval in case %d.\n",i + 1);
                return(-1);   
            }
        }
        if (ctx->PMCEN >= 0) {
            if (ctx->GDFYTyp[i] == 1 && fabs(get_data(ctx, ctx->PMCEN,i)) <= ctx->EPSI1)
                ctx->GDFYTyp[i] = 3;
        }
        if (ctx->GDFYTyp[i] == 1)
            ctx->GDFNEX++;
        else if (ctx->GDFYTyp[i] == 2)
            ctx->GDFNINT++;
        else
            ctx->GDFNCEN++;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_gdfytyp(n)    If n > 0 allocate GDFYTyp, otherwise free.          */
/*                      Return 0 if OK, -1 if error.                        */

int alloc_gdfytyp(TDAContext *ctx, int n)
{
    if (ctx->GDFYTypA > 0) {
        free((char *)ctx->GDFYTyp);
        memrq(ctx, -ctx->GDFYTypA,sizeof(short));
        ctx->GDFYTypA = 0;
    }
    if (n > 0) {
        if (!(ctx->GDFYTyp = (short *)calloc((size_t)(n),sizeof(short)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(short));
        ctx->GDFYTypA = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_gdfgen(nu,ndim)   if nu and ndim > 0 allocate basic data          */
/*                          structures, otherwise free.                     */
/*                          Return 0 if OK, -1 if error.                    */

int alloc_gdfgen(TDAContext *ctx, int nu,int ndim)
{
    int err = -1;

    if (nu == 0 && ndim == 0) {
        err = 0;
        goto GDFGFin;
    }
    if (ndim > 0) {
        if (!(ctx->GDFDim = (int *)calloc((size_t)(ndim),sizeof(int))))  
            goto GDFGFin; 
        memrq(ctx, ndim,sizeof(int));
        ctx->GDFDimA = ndim;

        if (!(ctx->GDFDCnt = (int *)calloc((size_t)(ndim),sizeof(int))))  
            goto GDFGFin; 
        memrq(ctx, ndim,sizeof(int));
        ctx->GDFDCntA = ndim;

        if (!(ctx->GDFDInt = (int *)calloc((size_t)(ndim),sizeof(int))))  
            goto GDFGFin; 
        memrq(ctx, ndim,sizeof(int));
        ctx->GDFDIntA = ndim;

        if (!(ctx->GDFDCen = (int *)calloc((size_t)(ndim),sizeof(int))))  
            goto GDFGFin; 
        memrq(ctx, ndim,sizeof(int));
        ctx->GDFDCenA = ndim;

        if (!(ctx->GDFDL = (double *)calloc((size_t)(ndim),sizeof(double))))  
            goto GDFGFin; 
        memrq(ctx, ndim,sizeof(double));
        ctx->GDFDLA = ndim;

        if (!(ctx->GDFDH = (double *)calloc((size_t)(ndim),sizeof(double))))  
            goto GDFGFin; 
        memrq(ctx, ndim,sizeof(double));
        ctx->GDFDHA = ndim;

        if (!(ctx->GDFBLen = (double *)calloc((size_t)(ndim),sizeof(double))))  
            goto GDFGFin; 
        memrq(ctx, ndim,sizeof(double));
        ctx->GDFBLenA = ndim;
    }
    if (nu > 0) {
        if (!(ctx->GDFPtr = (int *)calloc((size_t)(nu),sizeof(int))))  
            goto GDFGFin;
        memrq(ctx, nu,sizeof(int));
        ctx->GDFPtrA = nu;

        if (!(ctx->GDFULen = (short *)calloc((size_t)(nu),sizeof(short))))  
            goto GDFGFin;
        memrq(ctx, nu,sizeof(short));
        ctx->GDFULenA = nu;
    }
    if (nu > 0 && ndim > 0) {
        if (!(ctx->GDFMarg = (short *)calloc((size_t)(nu) * (size_t)(ndim) + 1,sizeof(short))))  
            goto GDFGFin;
        memrq(ctx, nu * ndim + 1,sizeof(short));
        ctx->GDFMargA = nu * ndim + 1;

        if (!(ctx->GDFMLen = (short *)calloc((size_t)(nu),sizeof(short))))  
            goto GDFGFin;
        memrq(ctx, nu,sizeof(short));
        ctx->GDFMLenA = nu;

        if (!(ctx->GDFMUCnt = (int *)calloc((size_t)(nu),sizeof(int))))  
            goto GDFGFin;
        memrq(ctx, nu,sizeof(int));
        ctx->GDFMUCntA = nu;

        if (!(ctx->GDFMTyp = (int *)calloc((size_t)(nu),sizeof(int))))  
            goto GDFGFin;
        memrq(ctx, nu,sizeof(int));
        ctx->GDFMTypA = nu;
    }
    return(0);

GDFGFin:
    if (err)
        p_err(ctx, -2,1);

    if (ctx->GDFDimA > 0) {
        free((char *)ctx->GDFDim);
        memrq(ctx, -ctx->GDFDimA,sizeof(int));
        ctx->GDFDimA = 0;
    }
    if (ctx->GDFDCntA > 0) {
        free((char *)ctx->GDFDCnt);
        memrq(ctx, -ctx->GDFDCntA,sizeof(int));
        ctx->GDFDCntA = 0;
    }
    if (ctx->GDFDIntA > 0) {
        free((char *)ctx->GDFDInt);
        memrq(ctx, -ctx->GDFDIntA,sizeof(int));
        ctx->GDFDIntA = 0;
    }
    if (ctx->GDFDCenA > 0) {
        free((char *)ctx->GDFDCen);
        memrq(ctx, -ctx->GDFDCenA,sizeof(int));
        ctx->GDFDCenA = 0;
    }
    if (ctx->GDFDLA > 0) {
        free((char *)ctx->GDFDL);
        memrq(ctx, -ctx->GDFDLA,sizeof(double));
        ctx->GDFDLA = 0;
    }
    if (ctx->GDFDHA > 0) {
        free((char *)ctx->GDFDH);
        memrq(ctx, -ctx->GDFDHA,sizeof(double));
        ctx->GDFDHA = 0;
    }
    if (ctx->GDFBLenA > 0) {
        free((char *)ctx->GDFBLen);
        memrq(ctx, -ctx->GDFBLenA,sizeof(double));
        ctx->GDFBLenA = 0;
    }
    if (ctx->GDFPtrA > 0) {
        free((char *)ctx->GDFPtr);
        memrq(ctx, -ctx->GDFPtrA,sizeof(int));
        ctx->GDFPtrA = 0;
    }
    if (ctx->GDFULenA > 0) {
        free((char *)ctx->GDFULen);
        memrq(ctx, -ctx->GDFULenA,sizeof(short));
        ctx->GDFULenA = 0;
    }
    if (ctx->GDFMargA > 0) {
        free((char *)ctx->GDFMarg);
        memrq(ctx, -ctx->GDFMargA,sizeof(short));
        ctx->GDFMargA = 0;
    }
    if (ctx->GDFMLenA > 0) {
        free((char *)ctx->GDFMLen);
        memrq(ctx, -ctx->GDFMLenA,sizeof(short));
        ctx->GDFMLenA = 0;
    }
    if (ctx->GDFMUCntA > 0) {
        free((char *)ctx->GDFMUCnt);
        memrq(ctx, -ctx->GDFMUCntA,sizeof(int));
        ctx->GDFMUCntA = 0;
    }
    if (ctx->GDFMTypA > 0) {
        free((char *)ctx->GDFMTyp);
        memrq(ctx, -ctx->GDFMTypA,sizeof(int));
        ctx->GDFMTypA = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdf_gcheck()    Check group structure and create data structures.       */
/*                  GDFNU = number of units                                 */
/*                  GDFNDim = number of dimensions                          */
/*                  GDFGRP = 1 if with groups, otherwise 0                  */
/*                                                                          */
/*                  GDFMTyp = 1 exact                                       */
/*                            2 i-censored                                  */
/*                            3 r-censored                                  */
/*                            4 i-censored and r-censored                   */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int gdf_gcheck(TDAContext *ctx)
{
    register int i,j,k,l,ii;
    int err,n,i0,l1 = 0,iid,il1,t,tt,fnd,cen2,cen3;
    double id = 0.0,id0;

    err = -1;
    if (ctx->NOC < 1) {
        printf1(ctx, "Error: no cases in the data matrix.\n");
        return(-1);
    }
    if (ctx->PM1NV == 0) {           /* no group structure */

        ctx->GDFGRP = 0;
        ctx->GDFNU = ctx->NOC;
        ctx->GDFNDim = 1;

        if (!(ctx->GDFDim = (int *)calloc(1,sizeof(int)))) {
            p_err(ctx, -2,1);
            goto GDFCFin; 
        }
        memrq(ctx, 1,sizeof(int));
        ctx->GDFDimA = 1;     
        ctx->GDFDim[0] = 1;

        err = 0;
        goto GDFCON;
    }
    if (ctx->PM1NV < 2) {
        printf1(ctx, "Error: grp parameter requires at least two variables.\n");
        goto GDFCFin;
    }
    if (vsort(ctx, ctx->PM1NV,ctx->PM1VIdx,1,1,0))     /* sort, create ptr VSORTPtr[] */
        goto GDFCFin;

    /* check units and observations */

    if (alloc_acm(ctx, ctx->NOC + 1))     /* used to sort L1 values */
        goto GDFCFin;

    ii = l = 0;
    ctx->GDFNU = 0;
    j = ctx->PM1VIdx[0];
    k = ctx->PM1VIdx[1];
    id0 = get_data(ctx, j,ctx->VSORTPtr[0]) - 1.0;
    for (i = 0; i < ctx->NOC; ++i) {
        n = (int)get_data(ctx, k,i);
        if (n < 1) {
            printf1(ctx, "Error: level variables must be positive integers.\n");
            goto GDFCFin;
        }
        ctx->AcM[ii++] = n;
        id = get_data(ctx, j,ctx->VSORTPtr[i]);
        if (id != id0) {
            ctx->GDFNU++;
            id0 = id;   
            l = 0;
        }
        l++;
    }
    if (sorti(ctx, ii,ctx->AcM,0))        /* sort and calculate GDFNDim */
        goto GDFCFin;

    j = ctx->AcM[0];
    ctx->GDFNDim = 1;
    for (i = 1; i < ii; ++i) {
        l = ctx->AcM[i];
        if (l != j) {
            ctx->GDFNDim++;
            j = l;
        }
    }
    if (alloc_gdfgen(ctx, ctx->GDFNU,ctx->GDFNDim))
        goto GDFCFin;

    ctx->GDFDim[0] = ctx->AcM[0];                 /* create GDFDim[] */
    k = 0;
    for (i = 1; i < ii; ++i) {
        l = ctx->AcM[i];
        if (l != ctx->GDFDim[k])  
            ctx->GDFDim[++k] = l;
    }
    if (alloc_acm(ctx, ctx->GDFNDim + 1))
        goto GDFCFin;

    ctx->GDFNTyp = 0;
    iid = ctx->PM1VIdx[0];
    il1 = ctx->PM1VIdx[1];
    i0 = 0;
    k = l = i = 0;
    id0 = get_data(ctx, iid,ctx->VSORTPtr[0]);
    cen2 = cen3 = 0;
    while (i <= ctx->NOC) {
        if (i < ctx->NOC) {
            ii = ctx->VSORTPtr[i];
            id = get_data(ctx, iid,ii);
            l1 = (int)get_data(ctx, il1,ii);
            for (t = 0; t < ctx->GDFNDim; ++t) {
                if (l1 == ctx->GDFDim[t]) {
                    ctx->GDFDCnt[t] += 1;
                    if (ctx->GDFYTyp[ii] == 2)
                        ctx->GDFDInt[t] += 1;
                    else if (ctx->GDFYTyp[ii] == 3)
                        ctx->GDFDCen[t] += 1;
                    break;
                }
            }
        }
        if (id != id0 || i == ctx->NOC) {

            tt = -1;
            for (t = 0; t < ctx->GDFNTyp; ++t) {
                if (ctx->GDFMLen[t] != l)
                    continue;
                fnd = 1;
                for (j = 0; j < l; ++j) {
                    if (ctx->GDFMarg[t * ctx->GDFNDim + j] != ctx->AcM[j + 1]) {
                        fnd = 0;
                        break;
                    }
                }
                if (fnd) {
                    tt = t;
                    break;
                }
            }
            if (tt < 0) {
                tt = ctx->GDFNTyp;
                for (j = 0; j < l; ++j) {
                    ctx->GDFMarg[tt * ctx->GDFNDim + j] = (short)(ctx->AcM[j + 1]);
                    if (j > 0 && ctx->AcM[j] == ctx->AcM[j + 1]) {
                        printf1(ctx, "Error in unit %g. Observations in the same unit\n",id0);
                        printf1(ctx, "must belong to different dimensions.\n");
                        goto GDFCFin;
                    }
                }
                ctx->GDFMLen[tt] = (short)(l);
                ctx->GDFNTyp++;
            }
            ctx->GDFPtr[k] = i0;
            ctx->GDFULen[k] = (short)(l);
            ctx->GDFMUCnt[tt] += 1;
            ctx->GDFMTyp[tt] = imax(ctx, ctx->GDFMTyp[tt],cen2 + 2 * cen3 + 1);
            k++;
            i0 = i;
            cen2 = cen3 = l = 0;
            id0 = id;
        }
        if (i == ctx->NOC)
            break;

        if (ctx->GDFYTyp[ii] == 2)
            cen2 = 1;
        else if (ctx->GDFYTyp[ii] == 3)
            cen3 = 1;
        l++;
        ctx->AcM[l] = l1;
        i++;
    }
    ctx->GDFGRP = 1;
    err = 0;

GDFCON:
    printf1(ctx, "Number of cases: %d\n",ctx->NOC);
    printf1(ctx, "Number of units: %d\n",ctx->GDFNU);
    printf1(ctx, "Number of dimensions (level 1): %d\n",ctx->GDFNDim);
    printf1(ctx, "Number of exact cases: %d\n",ctx->GDFNEX);
    printf1(ctx, "Number of interval censored cases: %d\n",ctx->GDFNINT);
    printf1(ctx, "Number of right censored cases: %d\n",ctx->GDFNCEN);

    if (ctx->GDFGRP) {
        printf1(ctx, "Number of marginal patterns: %d\n\n",ctx->GDFNTyp);
        gdf_pmtyp(ctx);

        printf1(ctx, "Dimension  level-1  observations  i-censored  r-censored\n");
        for (i = 0; i < ctx->GDFNDim; ++i)  
            printf1(ctx, "%6d     %5d    %12d  %10d  %10d\n",i + 1,
                              ctx->GDFDim[i],ctx->GDFDCnt[i],ctx->GDFDInt[i],ctx->GDFDCen[i]);
    }
    newline(ctx);

GDFCFin:
    alloc_acm(ctx, 0);
    alloc_acn(ctx, 0);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdf_pmtyp()     Print info about marginal patterns.                     */

void gdf_pmtyp(TDAContext *ctx)
{
    register int j,t;
    int l,d;
    
    printf1(ctx, "Pattern  type     units  dimensions\n");

    for (t = 0; t < ctx->GDFNTyp; ++t) {
        printf1(ctx, "%5d  %5d %10d  ",t + 1,ctx->GDFMTyp[t],ctx->GDFMUCnt[t]);
        l = ctx->GDFMLen[t];
        for (j = 0; j < l; ++j) {
            d = ctx->GDFMarg[t * ctx->GDFNDim + j];
            printf1(ctx, "%2d ",d);
        }
        newline(ctx);
    } 
    newline(ctx);
}   

/* ------------------------------------------------------------------------ */
/*  gdf_marg(opt)   create marginal distr. functions, write to output file. */
/*                  opt = 0 : distribution function                         */
/*                  opt = 1 : survivor function                             */
/*                  opt = 2 : expected values                               */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int gdf_marg(TDAContext *ctx, int opt)
{
    register int i,n;
    int err,d,dl,m,nn,nc,il1 = 0;

    err = -1;
    gdf_pmet(ctx, 1,opt);

    for (d = 0; d < ctx->GDFNDim; ++d) {

        dl = ctx->GDFDim[d];
        if (ctx->GDFGRP) {
            il1 = ctx->PM1VIdx[1];           /* index of L1 variable */
            nn = ctx->GDFDCnt[d];
            nc = ctx->GDFDCen[d];
        }
        else {
            nn = ctx->NOC;
            nc = ctx->GDFNCEN;
        }
        if (alloc_acy(ctx, nn))
            goto GDFMARGFin;

        if (nc > 0) {                   /* AcNS used for censoring indicator */
            if (alloc_acns(ctx, nn))
                goto GDFMARGFin;
            if (alloc_acu(ctx, nn))          /* used for new y values */
                goto GDFMARGFin;
        }
        if (alloc_acw(ctx, nn))
            goto GDFMARGFin;

        n = 0;
        for (i = 0; i < ctx->NOC; ++i) {
            if (ctx->GDFGRP == 0 || (int)get_data(ctx, il1,i) == dl) {
                ctx->AcY[n] = get_data(ctx, ctx->PMYL,i);
                if (nc > 0) {
                    if (ctx->GDFYTyp[i] == 3)
                        ctx->AcNS[n] = 0;
                    else
                        ctx->AcNS[n] = 1;
                }
                n++;
            }
        }
        if (nc == 0) {                      /* only exact observations */

            printf1(ctx, "Dimension %3d: EDF with exact data.\n",dl);

            m = gdf_edf1(ctx, n,ctx->AcY,ctx->AcW,opt);    /* get edf */         
            if (m < 1)
                goto GDFMARGFin;
            if (opt == 2)
                gdf_pdat1a(ctx, dl,m,ctx->AcY,ctx->AcW,ctx->AcNS,1);     /* print to output file */
            else
                gdf_pdat1(ctx, dl,m,ctx->AcY,ctx->AcW);        
        }
        else {                              /* exact and r-censored obs. */ 

            printf1(ctx, "Dimension %3d: Kaplan-Meier.\n",dl);

            m = gdf_edf2(ctx, n,ctx->AcY,ctx->AcNS,ctx->AcU,ctx->AcW,opt);   /* get edf */         
            if (m < 1)
                goto GDFMARGFin;
            if (opt == 2)
                gdf_pdat1a(ctx, dl,m,ctx->AcU,ctx->AcW,ctx->AcNS,0);     /* print to output file */
            else
                gdf_pdat1(ctx, dl,m,ctx->AcU,ctx->AcW);        
        }
    }
    err = 0;

GDFMARGFin:
    alloc_acns(ctx, 0);
    alloc_acu(ctx, 0);
    alloc_acw(ctx, 0);
    alloc_acy(ctx, 0);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdf_pdat1(d,n,y,f)    Print distribution function to output file.       */
   
void gdf_pdat1(TDAContext *ctx, int d,int n,double *y,double *f)
{
    register int i;

    for (i = 0; i < n; ++i) {
        fprintf(ctx->PMFd,"%4d ",d);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,y[i]);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,f[i]);
        fprintf(ctx->PMFd,"\n");
        ctx->GDFWNOC++;
#ifdef TDA_R_PACKAGE
        {
            double erow[3];
            erow[0] = (double)d;
            erow[1] = y[i];
            erow[2] = f[i];
            tda_export_row(ctx, "gdf.table", erow, 3);
        }
#endif
    }
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "gdf.table");
#endif
}

/* ------------------------------------------------------------------------ */
/*  gdf_pdat1a(d,n,y,f,cen,opt)   Print distribution function to output     */
/*                                file. If opt != 0 do not use cen but      */
/*                                always print 1.                           */
   
void gdf_pdat1a(TDAContext *ctx, int d,int n,double *y,double *f,short *cen,int opt)
{
    register int i,m;

    m = 1;
    for (i = 0; i < n; ++i) {
        fprintf(ctx->PMFd,"%4d ",d);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,y[i]);
        if (opt == 0)
            m = cen[i];                
        fprintf(ctx->PMFd,"%2d ",m);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,f[i]);
        fprintf(ctx->PMFd,"\n");
        ctx->GDFWNOC++;
    }
}

/* ------------------------------------------------------------------------ */
/*  gdf_edf1(n,y,f,opt)                                                     */
/*                                                                          */
/*                      calculate empirical distr function with exact data. */
/*                      y[i] (i=0,...,n-1) contain the input data. Return   */
/*                      y[j],f[j] for j = 0,...,m-1 where y[j] are the      */
/*                      unique values and f[j] the corresponding value of   */
/*                      the edf.                                            */
/*                                                                          */
/*  opt = 0 : distr. function                                               */
/*        1 : survivor function                                             */
/*        2 : expected = observed values (only for completeness)            */
/*            in this case, do not sort the data.                           */
/*                                                                          */
/*  Return m > 0 if OK, or -1 if error.                                     */
   
int gdf_edf1(TDAContext *ctx, int n,double *y,double *f,int opt)
{
    register int i,k,m;
    double y0,y1 = 0.0,f0;

    if (opt == 2) {
        for (i = 0; i < n; ++i)  
            f[i] = y[i];
        return(n);
    }
    if (alloc_actmp(ctx, n))
        return(-1);

    for (i = 0; i < n; ++i)
        ctx->AcTmp[i] = y[i];

    if (sortd(ctx, n,ctx->AcTmp,0))
        return(-1);    

    y0 = ctx->AcTmp[0];
    m = 0;
    k = 1;
    f0 = 0.0;
    for (i = 1; i <= n; ++i) {
        if (i < n)
            y1 = ctx->AcTmp[i];

        if (i == n || y1 > y0) {
            f0 += (double)k / (double)n;
            y[m] = y0;
            f[m++] = f0;
            y0 = y1;
            k = 0;
        }
        k++;
    }
    if (opt == 1) {
        for (i = 0; i < m; ++i)
            f[i] = 1.0 - f[i];
    }
    alloc_actmp(ctx, 0);
    return(m);
}

/* ------------------------------------------------------------------------ */
/*  gdf_edf2(n,y,cen,y1,h,opt)                                              */
/*                                                                          */
/*  y[i], cen[i] (i = 0,...,n-1) are the input data. y[i] is right          */
/*  censored if cen[i] = 0.                                                 */
/*                                                                          */
/*  If opt=0 or opt=1:                                                      */
/*  the function creates y1[j] (j=0,...,m-1) containing the different       */
/*  values in y[i]. And for each y1[j] it calculates the corresponding      */
/*  Kaplan-Meyer distribution function (if opt=0), or the survivor          */
/*  function (if opt=1).                                                    */
/*  The largest observation is always assumed to be uncensored.             */
/*                                                                          */
/*  If opt=2, the function returns in y1[0,...,n-1] the input               */  
/*  values and in h[i] the corresponding expected values. In this case the  */
/*  values will not be sorted!                                              */
/*                                                                          */
/*  Return m > 0 if OK, or -1 if error.                                     */

int gdf_edf2(TDAContext *ctx, int n,double *y,short *cen,double *y1,double *h,int opt)
{
    register int i,j,k,l;
    int m;
    double tmp,tmp1;

    if (alloc_ack(ctx, n))
        return(-1);
    if (alloc_actmp(ctx, n))
        return(-1);

    if (sortdp2a(ctx, n,y,cen,ctx->AcK,0))      /* sort */
        return(-1);    

    tmp = 1.0 / (double)n;
    m = n - 1;
    for (i = 0; i < n; ++i) {
        j = ctx->AcK[i];
        ctx->AcTmp[j] += tmp;
        if (cen[j] == 0 && i < n - 1) {
            if (m > 0) {
                tmp1 = ctx->AcTmp[j] / (double)m;
                for (k = i + 1; k < n; ++k)  
                    ctx->AcTmp[ctx->AcK[k]] += tmp1;
            }
            ctx->AcTmp[j] = 0.0;
        }
        m--;
    }
    if (opt == 2) {         /* expected values */
        for (i = 0; i < n; ++i) {
            j = ctx->AcK[i];
            y1[j] = y[j];
            if (cen[j] || i == n - 1)  
                h[j] = y[j];
            else {
                tmp1 = tmp = 0.0;
                for (k = i + 1; k < n; ++k) {
                    l = ctx->AcK[k];
                    tmp1 += y[l] * ctx->AcTmp[l];
                    tmp += ctx->AcTmp[l];
                }
                if (tmp > 0.0)
                    h[j] = tmp1 / tmp;
                else
                    h[j] = 0.0;
            }
        }
        m = n;
    }
    else {
        m = 0;
        tmp = 0.0;
        j = ctx->AcK[0];
        y1[m] = y[j];
        h[m] = ctx->AcTmp[j];

        for (i = 1; i < n; ++i) {
            j = ctx->AcK[i];
            if (ctx->AcTmp[j] < ctx->EPSI1)
                continue;

            if (y[j] == y1[m])
                h[m] += ctx->AcTmp[j];
            else {
                m++;
                y1[m] = y[j];
                h[m] = h[m - 1] + ctx->AcTmp[j];
            }
        }
        m++;
        if (opt == 1) {                 /* surv function */
            for (i = 0; i < m; ++i)
                h[i] = 1.0 - h[i];
        }
    }
    alloc_ack(ctx, 0);
    alloc_actmp(ctx, 0);
    return(m);
}

/* ------------------------------------------------------------------------ */
/*  gdf_joint(met,opt)    joint estimation.                                 */
/*                                                                          */
/*  met = 2 : iterative estimation based on grid                            */
/*        3 : local Kaplan-Meier estimation                                 */
/*                                                                          */
/*  opt = 0 : distr. function                                               */
/*        1 : survivor function                                             */
/*        2 : expected values                                               */
/*                                                                          */
/*  If iterative procedure, then if PMProtFDef, write diagnostic            */
/*  information to PMProtFd.                                                */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int gdf_joint(TDAContext *ctx, int met,int opt)
{
    int err,t,m,nu,ndim;

    err = -1;
    
    gdf_pmet(ctx, met,opt);
    if (ctx->GDFGRP == 0) {
        printf1(ctx, "Error: grp parameter required.\n");
        return(-1);
    }

    for (t = 0; t < ctx->GDFNTyp; ++t) {     /* for all patterns */

        nu = ctx->GDFMUCnt[t];
        ndim = ctx->GDFMLen[t];

        if (ctx->PMProtFDef) {
            fprintf(ctx->PMProtFd,"Current marginal pattern: %d  Typ: %d\n",t+1,ctx->GDFMTyp[t]);
            fprintf(ctx->PMProtFd,"Number of units: %d  dimensions: %d\n",nu,ndim);
        }   
        if (alloc_aci(ctx, nu * ndim + 1))   /* pointer for cases */
            goto GDFJFin;

        gdf_cptr(ctx, t,nu,ndim,ctx->AcI);        /* create pointer to data matrix cases */

        if (ctx->GDFMTyp[t] == 1) {          /* only exact observations */

            printf1(ctx, "\nPattern %d: EDF with exact data.\n",t + 1);
            if (gdf_domain(ctx, t,nu,ndim,ctx->AcI,NULL,1,0.0))
                goto GDFJFin;

            if (alloc_acx(ctx, nu * ndim + 1))   /* values */
                goto GDFJFin;
            if (alloc_acw(ctx, nu + 1))          /* distr function */
                goto GDFJFin;

            m = gdf_edf1m(ctx, t,nu,ndim,ctx->AcI,ctx->AcX,ctx->AcW,opt);
            if (m < 1)
                goto GDFJFin;

            if (opt == 2)
                gdf_pdat2a(ctx, t,m,ndim,ctx->AcX,ctx->AcI,1);  
            else
                gdf_pdat2(ctx, t,m,ndim,ctx->AcX,ctx->AcW);    

        }
        else if (ctx->GDFMTyp[t] == 3) {     /* exact and r-censored observations */

            printf1(ctx, "\nPattern %d. Units: %d, dimensions: %d\n",t + 1,nu,ndim);

            if (alloc_acv(ctx, nu * ndim + 1))           /* values/lower bounds */
                goto GDFJFin;
            if (alloc_acw(ctx, nu + 1))                  /* distr function */
                goto GDFJFin;

            if (met == 2) {
                printf1(ctx, "Iterative procedure.\n");
                printf1(ctx, "Max iterations: %d  tolerance: %g\n",ctx->MxIter,ctx->TOLF);

                if (gdf_domain(ctx, t,nu,ndim,ctx->AcI,NULL,1,0.0))
                    goto GDFJFin;

                gdf_getb(ctx, ctx->PMN,ndim,1);                   /* get number of boxes */
    
                if (gdf_edf3m(ctx, t,nu,ndim,ctx->AcI,NULL,ctx->AcV,ctx->AcW,opt,1,ctx->MxIter,ctx->TOLF))
                    goto GDFJFin;
            }
            else {
                printf1(ctx, "Local Kaplan-Meier.\n");
                if (gdf_domain(ctx, t,nu,ndim,ctx->AcI,NULL,1,ctx->PMD))
                    goto GDFJFin;

                if (gdf_edf4m(ctx, t,nu,ndim,ctx->AcI,NULL,ctx->AcV,ctx->AcW,ctx->PMD,opt,1))
                    goto GDFJFin;
            }
            if (opt == 2)
                gdf_pdat2a(ctx, t,nu,ndim,ctx->AcV,ctx->AcI,0);  
            else
                gdf_pdat2(ctx, t,nu,ndim,ctx->AcV,ctx->AcW);    
        }
    }
    err = 0;

GDFJFin:
    alloc_aci(ctx, 0);
    alloc_acx(ctx, 0);
    alloc_acw(ctx, 0);
    alloc_acv(ctx, 0);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdf_pdat2(t,m,ndim,val,f)                                               */
/*                                                                          */
/*  write joint distr function to output file. t is number of pattern.      */
/*  val is (m,ndim), f is m vector with masses.                             */
   
void gdf_pdat2(TDAContext *ctx, int t,int m,int ndim,double *val,double *f)
{
    register int i,j;

    for (i = 0; i < m; ++i) {
        fprintf(ctx->PMFd,"%3d  ",t + 1);
        for (j = 0; j < ctx->GDFNDim; ++j)
            fprintf(ctx->PMFd,"%2d ",ctx->GDFMarg[t * ctx->GDFNDim + j]);

        for (j = 0; j < ndim; ++j)
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,val[i * ndim + j]);

        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,f[i]);
        fprintf(ctx->PMFd,"\n");
        ctx->GDFWNOC++;
#ifdef TDA_R_PACKAGE
        {
            double *erow = (double *)malloc((size_t)(2 + ctx->GDFNDim +
                                                     ndim) *
                                            sizeof(double));
            if (erow != NULL) {
                int ec = 0;
                erow[ec++] = (double)(t + 1);
                for (j = 0; j < ctx->GDFNDim; ++j)
                    erow[ec++] =
                        (double)ctx->GDFMarg[t * ctx->GDFNDim + j];
                for (j = 0; j < ndim; ++j)
                    erow[ec++] = val[i * ndim + j];
                erow[ec++] = f[i];
                tda_export_row(ctx, "gdf.table", erow, ec);
                free(erow);
            }
        }
#endif
    }
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "gdf.table");
#endif
}

/* ------------------------------------------------------------------------ */
/*  gdf_pdat2a(t,m,ndim,val,ptr,opt)                                        */
/*                                                                          */
/*  write expected values. If opt = 1 val contains the exact values and     */
/*  ptr is not used. If opt = 0, val contains the expected values and       */
/*  ptr provides pointer to original observations. In this case, m = nu.    */
   
void gdf_pdat2a(TDAContext *ctx, int t,int m,int ndim,double *val,int *ptr,int opt)
{
    register int i,j,ii;
    int n;
    double tmp;

    n = 1;
    for (i = 0; i < m; ++i) {
        ii = i * ndim;
        fprintf(ctx->PMFd,"%3d  ",t + 1);
        for (j = 0; j < ctx->GDFNDim; ++j)
            fprintf(ctx->PMFd,"%2d ",ctx->GDFMarg[t * ctx->GDFNDim + j]);

        if (opt) {
            for (j = 0; j < ndim; ++j)
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,val[ii + j]);
            for (j = 0; j < ndim; ++j)
                fprintf(ctx->PMFd,"%d ",n);
            for (j = 0; j < ndim; ++j)
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,val[ii + j]);
        }
        else {
            for (j = 0; j < ndim; ++j) {
                tmp = get_data(ctx, ctx->PMYL,ptr[ii + j]);
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
            }
            for (j = 0; j < ndim; ++j) {
                n = 1;
                if (ctx->GDFYTyp[ptr[ii + j]] == 3)
                    n = 0;
                fprintf(ctx->PMFd,"%d ",n);
            }
            for (j = 0; j < ndim; ++j)
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,val[ii + j]);
        }
        fprintf(ctx->PMFd,"\n");
        ctx->GDFWNOC++;
    }
}

/* ------------------------------------------------------------------------ */
/*  gdf_cptr(t,nu,ndim,ptr)                                                 */
/*                                                                          */
/*  Create pointer ptr[i,j] (i=0,...,nu-1, j=0,...,ndim-1)  such that       */
/*  ptr[i,j] is pointer to data matrix row for the corresponding case.      */
/*                                                                          */
   
void gdf_cptr(TDAContext *ctx, int t,int nu,int ndim,int *ptr)
{
    register int i,j,k,l;
    int l1,il1,len,t1,d1,iu;

    il1 = ctx->PM1VIdx[1];
    iu = 0;
    for (i = 0; i < ctx->GDFNU; ++i) {
        if (iu >= nu)
            break;

        t1 = t * ctx->GDFNDim;
        d1 = ctx->GDFMarg[t1]; 
        len = ctx->GDFULen[i];

        if (len != ndim)     
            continue;

        k = ctx->GDFPtr[i];   
        for (l = 0; l < len; ++l) {
            d1 = ctx->GDFMarg[t1++]; 
            j = ctx->VSORTPtr[k + l];

            get_data(ctx, ctx->PMYL,j);
            l1 = (int)get_data(ctx, il1,j);
            if (l1 != d1) {
                d1 = -1;
                break;
            }
            ptr[iu * ndim + l] = j;
        }
        if (d1 < 0)
            continue;
        iu++;
    }
    if (iu != nu) {
        printfe(ctx, "ERROR in gdf_cptr (iu=%d, nu=%d)\n",iu,nu);
        gerr_exit(ctx, 213);
    }
}

/* ------------------------------------------------------------------------ */
/*  gdf_domain(t,nu,ndim,ptr,yval,opt,d)                                    */
/*                                                                          */
/*  Also create bounds for domain in GDFDL and GDFDH.                       */
/*  If opt = 1 print to standard output, if opt = 2 print to prot file.     */
/*                                                                          */
/*  If yval != NULL, use these values: yval[i,j], i=0,nu-1,j=0,ndim-1       */
/*                                                                          */  
/*  Return 0 if OK, -1 if error (empty domain)                              */
   
int gdf_domain(TDAContext *ctx, int t,int nu,int ndim,int *ptr,double *yval,int opt,double d)
{
    register int i,j,k;
    int err;
    double y,dl,dh;

    err = 0;
    for (j = 0; j < ndim; ++j) {
        k = ptr[j];
        if (yval == NULL)
            dl = dh = get_data(ctx, ctx->PMYL,k);
        else                 
            dl = dh = yval[k];

        if (ctx->PMYH >= 0 && ctx->GDFYTyp[k] < 3) {
            if (yval == NULL)
                dh = get_data(ctx, ctx->PMYH,k);
        }

        for (i = 1; i < nu; ++i) {
            k = ptr[i * ndim + j];
            if (yval == NULL)
                y = get_data(ctx, ctx->PMYL,k);
            else
                y = yval[k];

            dl = dmin(ctx, dl,y);
            if (ctx->PMYH >= 0 && ctx->GDFYTyp[k] < 3) {
                if (yval == NULL)
                    y = get_data(ctx, ctx->PMYH,k);
            }
            dh = dmax(ctx, dh,y);
        }
        if (ctx->PMSC > 0.0) {
            dl -= ctx->PMSC;
            dh += ctx->PMSC;
        }
        ctx->GDFDL[j] = dl;
        ctx->GDFDH[j] = dh;
        if (dh <= dl + ctx->EPSI1)
            err = -1;
    }
    if (opt == 1) {  
        printf1(ctx, "Domain");
        if (ctx->PMSC > 0.0)
            printf1(ctx, " (sc=%g)",ctx->PMSC);
        newline(ctx);
        for (j = 0; j < ndim; ++j) {
            printf1(ctx, "%3d ",ctx->GDFMarg[t * ctx->GDFNDim + j]);
            rt_printf1_d(ctx, ctx->PMTFmtS,ctx->GDFDL[j]);
            rt_printf1_d(ctx, ctx->PMTFmtS,ctx->GDFDH[j]);
            if (d > 0.0)
                printf1(ctx, "  delta: %g",d * (ctx->GDFDH[j] - ctx->GDFDL[j]));
            newline(ctx);
        }
        newline(ctx);
    }
    else if (opt == 2 && ctx->PMProtFDef) {
        fprintf(ctx->PMProtFd,"Domain");
        if (ctx->PMSC > 0.0)
            fprintf(ctx->PMProtFd," (sc=%g)",ctx->PMSC);
        fprintf(ctx->PMProtFd,"\n");
        for (j = 0; j < ndim; ++j) {
            fprintf(ctx->PMProtFd,"%3d ",ctx->GDFMarg[t * ctx->GDFNDim + j]);
            rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMTFmtS,ctx->GDFDL[j]);
            rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMTFmtS,ctx->GDFDH[j]);
            if (d > 0.0)
                fprintf(ctx->PMProtFd,"  delta: %g",d * (ctx->GDFDH[j] - ctx->GDFDL[j]));
            fprintf(ctx->PMProtFd,"\n");
        }
        fprintf(ctx->PMProtFd,"\n");
    }
    if (err)
        printf1(ctx, "Error: domain is empty.\n");
    return(err);
}
       
/* ------------------------------------------------------------------------ */
/*  gdf_edf1m(t,nu,ndim,ptr,val,f,opt)                                      */
/*                                                                          */
/*  Calculate joint distribution (opt = 0) or survivor (opt != 0) function. */
/*  Assume exact observations only.                                         */
/*  For pattern t, ptr created by gdf_cptr. Return values in val[] and      */
/*  value of distr. function in f.                                          */
/*                                                                          */
/*  If opt=2 calculate expected (= observed) values. Simply copy sorted     */
/*  values into val[] and let f = 0.                                        */
/*                                                                          */  
/*  Return m > 0 if OK, or -1 if error.                                     */
   
int gdf_edf1m(TDAContext *ctx, int t,int nu,int ndim,int *ptr,double *val,double *f,int opt)
{
    (void)t;        /* unused: the signature is shared */
    register int i,j,k,l,ii,kk;
    int sflg,m,mm;
    double tmp;

    if (alloc_actmp(ctx, nu * ndim + 1))
        return(-1);
    if (alloc_acj(ctx, ndim + 1))
        return(-1);
    if (alloc_ack(ctx, nu + 1))
        return(-1);

    for (j = 1; j <= ndim; ++j)
        ctx->AcJ[j] = j;

    for (i = 0; i < nu; ++i) {
        l = i * ndim;
        for (j = 0; j < ndim; ++j) {
            k = ptr[l];
            l++;
            ctx->AcTmp[l] = get_data(ctx, ctx->PMYL,k);
        }
    }
    if (sortdpn(ctx, nu,ndim,ctx->AcTmp,ndim,ctx->AcJ,ctx->AcK))
        return(-1);

    if (opt == 2) {
        for (i = 0; i < nu; ++i) {
            k = (ctx->AcK[i + 1] - 1) * ndim + 1;
            ii = i * ndim;
            for (j = 0; j < ndim; ++j) 
                val[ii + j] = ctx->AcTmp[k + j];
        }
        m = nu;
        goto EDF1MFin;
    }

    tmp = 1.0 / (double)nu;

    if (opt == 0) {                     /* distribution function */
        m = 0;
        for (i = 0; i < nu; ++i) {
            k = (ctx->AcK[i + 1] - 1) * ndim + 1;
            if (i == 0) {
                for (j = 0; j < ndim; ++j) 
                    val[j] = ctx->AcTmp[k + j];
                f[0] = tmp;
            }
            else {
                sflg = 1;
                mm = m * ndim;
                for (j = 0; j < ndim; ++j) {
                    if (val[mm + j] != ctx->AcTmp[k + j]) {
                        sflg = 0;
                        break;
                    }
                }
                if (sflg)  
                    f[m] += tmp;
                else {
                    m++;
                    mm = m * ndim;
                    for (j = 0; j < ndim; ++j) 
                        val[mm + j] = ctx->AcTmp[k + j];
                    f[m] = tmp;
                    for (ii = 0; ii < i; ++ii) {
                        kk = (ctx->AcK[ii + 1] - 1) * ndim + 1;
                        sflg = 1;
                        for (j = 0; j < ndim; ++j) {
                            if (val[mm + j] < ctx->AcTmp[kk + j]) {
                                sflg = 0;
                                break;
                            }
                        }
                        if (sflg)
                            f[m] += tmp;
                    }   
                }
            }
        }
        m++;
    }
    else {                     /* survivor function */
        m = 0;
        for (i = nu - 1; i >= 0; --i) {
            k = (ctx->AcK[i + 1] - 1) * ndim + 1;
            if (i == nu - 1) {
                for (j = 0; j < ndim; ++j) 
                    val[j] = ctx->AcTmp[k + j];
                f[0] = 0.0;
            }
            else {
                sflg = 1;
                mm = m * ndim;
                for (j = 0; j < ndim; ++j) {
                    if (val[mm + j] != ctx->AcTmp[k + j]) {
                        sflg = 0;
                        break;
                    }
                }
                if (sflg == 0) {
                    m++;
                    mm = m * ndim;
                    for (j = 0; j < ndim; ++j) 
                        val[mm + j] = ctx->AcTmp[k + j];
                    f[m] = 0.0;
                    for (ii = i + 1; ii < nu; ++ii) {
                        kk = (ctx->AcK[ii + 1] - 1) * ndim + 1;
                        sflg = 1;
                        for (j = 0; j < ndim; ++j) {
                            if (val[mm + j] >= ctx->AcTmp[kk + j]) {
                                sflg = 0;
                                break;
                            }
                        }
                        if (sflg)
                            f[m] += tmp;
                    }   
                }
            }
        }
        m++;
    }
EDF1MFin:
    alloc_actmp(ctx, 0);
    alloc_acj(ctx, 0);
    alloc_ack(ctx, 0);
    return(m);
}

/* ------------------------------------------------------------------------ */
/*  gdf_edf3m(t,nu,ndim,ptr,yval,val,f,opt,prn,mxit,tol)                    */
/*                                                                          */
/*  Iterative estimation of distribution function with r-censored data.     */
/*                                                                          */
/*  Use pattern t, ptr created by gdf_cptr.                                 */
/*                                                                          */
/*  t = number of current marginal pattern.                                 */
/*  nu = number of units                                                    */
/*  ndim = number of dimensions in current marg pattern                     */
/*  ptr = ptr to observations, created by gdf_cptr().                       */
/*                                                                          */
/*  If yval == NULL use data matrix, otherwise use yval[0,...,NOC-1].       */
/*                                                                          */
/*  If prn != 0 print table with number of i/r-censored observations        */
/*  in the current pattern.                                                 */
/*                                                                          */
/*  Perform maximal mxit iterations. Tolerance is given by tol.             */  
/*                                                                          */  
/*  If opt = 0  return values in val[] and distribution function in f[]     */
/*  If opt = 1  return values in val[] and survivor function in f[]         */
/*  If opt = 2  return expected values in val[] and ignore f[].             */
/*                                                                          */
/*  Return 0 if OK, or -1 if error.                                         */
   
int gdf_edf3m(TDAContext *ctx, int t,int nu,int ndim,int *ptr,double *yval,double *val, double *f,int opt,int prn,int mxit,double tol)
{
    register int i,j,k,l;
    int fin,cen,n,np,iter,nred,ne,kp,nep;
    double y,tmp0,tmp,tmax,tol1;

    tmp0 = 1.0 / (double)nu;

    if (alloc_ack(ctx, nu))              /* number of boxes for unit */
        return(-1);
    if (alloc_acu(ctx, ndim))            /* temporary for y values */
        return(-1);
    if (alloc_acm(ctx, ndim))            /* for different use */
        return(-1);
    if (alloc_acj(ctx, ndim))            /* for different use */
        return(-1);
    if (alloc_acns(ctx, nu))             /* count censored dimensions */
        return(-1);
    if (alloc_acn(ctx, nu * ndim + 1))   /* pointer to grid */
        return(-1);
    if (alloc_acc(ctx, nu * ndim + 1))   /* set with status of observation */
        return(-1);

    /* determine pointer to the location of points in grid for
       all units, save in AcN[]. Set AcC to status of observation
       (1 - exact, 2 = icensored, 3 = r-censored).
       Set AcNS = -1 if unit falls outside of domain.       
                   0 if exact in all dimensions
                   1 otherwise     
    */

    if (prn)            /* print table */
        printf1(ctx, "Dimension  i-censored  r-censored  observations\n");

    for (i = 0; i < nu; ++i) {
        l = i * ndim;
        cen = 0;
        for (j = 0; j < ndim; ++j) {
            k = ptr[l + j];
            if (yval == NULL)
                y = get_data(ctx, ctx->PMYL,k);
            else
                y = yval[k];

            ctx->AcC[l + j] = (char)ctx->GDFYTyp[k];
            if (ctx->GDFYTyp[k] == 2) { 
                ctx->AcM[j] += 1;
                cen = 1;
            }   
            else if (ctx->GDFYTyp[k] == 3) {
                ctx->AcJ[j] += 1;
                cen = 1;
                y += ctx->EPSI1;             /* add small offset */
            }
            ctx->AcU[j] = y;
        }
        ctx->AcNS[i] = (short)(cen);

        for (j = ndim - 1; j >= 0; --j) {
            n = (int)ceil((ctx->AcU[j] - ctx->GDFDL[j]) / ctx->GDFBLen[j]);
            if (n < 1 || n > ctx->GDFNBOX) {
                ctx->AcNS[i] = -1;
                break;
            }
            ctx->AcN[l + j] = n;
        }
    }
    if (prn) {      /* print table with i/r-censored observations */
        for (j = 0; j < ndim; ++j)  
            printf1(ctx, "%6d    %11d %11d  %12d\n",
                             ctx->GDFMarg[t * ctx->GDFNDim + j],ctx->AcM[j],ctx->AcJ[j],nu);
        newline(ctx);
    }

    /*  create pointers for grid in AcJ */
    /*  first set AcJ > 0 for all boxes in grid that have an observation */
    /*  AcJ = 1 for exact observations, AcJ = 2 otherwise */
    /*  Also set AcK = number of boxes for unit */
    /*  count number of boxes that have exact observations in ne */

    ne = 0;

    if (alloc_acj(ctx, ctx->GDFTBOX + 1))  
        return(-1);

    if (ctx->PMProtFDef)  
        fprintf(ctx->PMProtFd,"\nLocation of observations in grid.\n");

    /***
    printf1(ctx, "ACN\n");
    for (i=0;i< nu; ++i) {
        for (j=0;j<ndim; ++j)
            printf1(ctx, "%d ",AcN[i * ndim + j]);
        newline(ctx);
    }
    ****/

    for (i = 0; i < nu; ++i) {      /* for all units */

        if (ctx->AcNS[i] < 0)            /* skip observations outside of grid */
            continue;   

        l = i * ndim;
        for (j = 0; j < ndim; ++j)   
            ctx->AcM[j] = ctx->AcN[l + j];

        np = gdf_boxptr(ctx, ndim,ctx->AcM);  /* pointer to first box of observation */

        if (ctx->PMProtFDef)  
            fprintf(ctx->PMProtFd,"%6d %3d %5d\n",i,ctx->AcNS[i],np);

        ctx->AcK[i] = 1;                 /* count boxes */
        if (ctx->AcNS[i] == 0) {         /* if exact */    
            ctx->AcJ[np] = 1;
            continue;
        }

        fin = 0;
        while (fin == 0) {
        
            if (ctx->AcJ[np] == 0)   
                ctx->AcJ[np] = 2;
        
        
            fin = 1;
            for (j = 0; j < ndim; ++j) {
                if (ctx->AcC[l + j] == 3 && ctx->AcM[j] < ctx->GDFNBOX) {
                    ctx->AcM[j] += 1;
                    for (k = 0; k < j; ++k)     /* reset */
                        ctx->AcM[k] = ctx->AcN[l + k];
                    fin = 0;
                    break;
                }
            }
            if (fin == 0) {
                np = gdf_boxptr(ctx, ndim,ctx->AcM);
                ctx->AcK[i] += 1;
                if (ctx->PMProtFDef)  
                    fprintf(ctx->PMProtFd,"%6d %3d %5d\n",i,ctx->AcNS[i],np);
            }
        }
    }
    if (ctx->PMProtFDef)  
        fprintf(ctx->PMProtFd,"\nPointers to reduced grid\n");

    /* AcJ[np] is index of np in reduced grid */
    /* nred is the number of locations in the reduced grid */

    


    nred = 0;
    for (i = 0; i < ctx->GDFTBOX; ++i) {
        if (ctx->AcJ[i] > 0) {
            if (ctx->AcJ[i] == 1)
                ne++;
            ctx->AcJ[i] = nred++;
            if (ctx->PMProtFDef)  
                fprintf(ctx->PMProtFd,"%4d %4d\n",i,ctx->AcJ[i]);
        }
        else
            ctx->AcJ[i] = -1;
    }
    if (ctx->PMProtFDef) {
        fprintf(ctx->PMProtFd,"Number of boxes, total: %d  reduced: %d\n",ctx->GDFTBOX,nred);
        fprintf(ctx->PMProtFd,"Number of boxes with exact observations: %d\n",ne);
    }

    /* create two arrays for reduced grid.
       AcXF for densities, AcYF for redistribution weights */

    if (alloc_acxf(ctx, nred))          
        return(-1);
    if (alloc_acyf(ctx, nred))          
        return(-1);

    /*  create initial weights in AcXF */

    for (i = 0; i < nu; ++i) {      /* for all units */

        if (ctx->AcNS[i] < 0)            /* skip observations outside of grid */
            continue;   

        l = i * ndim;
        for (j = 0; j < ndim; ++j) 
            ctx->AcM[j] = ctx->AcN[l + j];

        np = gdf_boxptr(ctx, ndim,ctx->AcM);  /* pointer to first box of observation */

        if (ctx->AcNS[i] == 0) {         /* if exact then unit mass */
            ctx->AcXF[ctx->AcJ[np]] = (float)((double)(ctx->AcXF[ctx->AcJ[np]]) + (tmp0));  
            continue;
        }
        fin = 0;
        while (fin == 0) {

            ctx->AcXF[ctx->AcJ[np]] = (float)((double)(ctx->AcXF[ctx->AcJ[np]]) + (tmp0 / (double)ctx->AcK[i]));

            fin = 1;
            for (j = 0; j < ndim; ++j) {
                if (ctx->AcC[l + j] == 3 && ctx->AcM[j] < ctx->GDFNBOX) {
                    ctx->AcM[j] += 1;
                    for (k = 0; k < j; ++k)     /* reset */
                        ctx->AcM[k] = ctx->AcN[l + k];
                    fin = 0;
                    break;
                }
            }
            if (fin == 0)  
                np = gdf_boxptr(ctx, ndim,ctx->AcM);
        }
    }
    if (ctx->PMProtFDef) {
        fprintf(ctx->PMProtFd,"\nInitial distribution\n");
        tmp = 0.0;
        for (j = 0; j < nred; ++j)  {
            tmp+=(double)(ctx->AcXF[j]);
            fprintf(ctx->PMProtFd,"%3d  %8.6f\n",j,(double)(ctx->AcXF[j]));
        }
        fprintf(ctx->PMProtFd,"Sum of values: %g\n",tmp);
    }

    /* perform iterations */

    tmax = 0.0;
    for (iter = 1; iter <= mxit; ++iter) {

        /* update distribution. Use weights in AcXF, create new
           densities in AcYF */

        for (i = 0; i < nu; ++i) {      /* for all units */

            if (ctx->AcNS[i] < 0)            /* skip observations outside of grid */
                continue;   

            l = i * ndim;
            for (j = 0; j < ndim; ++j) 
                ctx->AcM[j] = ctx->AcN[l + j];

            np = gdf_boxptr(ctx, ndim,ctx->AcM);  /* pointer to first box of observation */

            if (ctx->AcNS[i] == 0) {         /* if exact then unit mass */
                ctx->AcYF[ctx->AcJ[np]] = (float)((double)(ctx->AcYF[ctx->AcJ[np]]) + (tmp0));  
                continue;
            }
            tmp = 0.0;                  /* first find sum of weights */
            fin = 0;
            while (fin == 0) {

                tmp += (double)(ctx->AcXF[ctx->AcJ[np]]);
                fin = 1;
                for (j = 0; j < ndim; ++j) {
                    if (ctx->AcC[l + j] == 3 && ctx->AcM[j] < ctx->GDFNBOX) {
                        ctx->AcM[j] += 1;
                        for (k = 0; k < j; ++k)     /* reset */
                            ctx->AcM[k] = ctx->AcN[l + k];
                        fin = 0;
                        break;
                    }
                }
                if (fin == 0)  
                    np = gdf_boxptr(ctx, ndim,ctx->AcM);
            }

            /* now add new weights */

            if (tmp < ctx->EPSI1)
                continue;
            tmp = tmp0 / tmp;

            for (j = 0; j < ndim; ++j) 
                ctx->AcM[j] = ctx->AcN[l + j];

            np = gdf_boxptr(ctx, ndim,ctx->AcM);  /* pointer to first box of observation */

            fin = 0;                /* add new weights */
            while (fin == 0) {

                ctx->AcYF[ctx->AcJ[np]] = (float)((double)(ctx->AcYF[ctx->AcJ[np]]) + ((double)(ctx->AcXF[ctx->AcJ[np]]) * tmp));

                fin = 1;
                for (j = 0; j < ndim; ++j) {
                    if (ctx->AcC[l + j] == 3 && ctx->AcM[j] < ctx->GDFNBOX) {
                        ctx->AcM[j] += 1;
                        for (k = 0; k < j; ++k)     /* reset */
                            ctx->AcM[k] = ctx->AcN[l + k];
                        fin = 0;
                        break;
                    }
                }
                if (fin == 0)  
                    np = gdf_boxptr(ctx, ndim,ctx->AcM);
            }
        }
        tmax = 0.0;
        for (j = 0; j < nred; ++j) {
            if ((double)(ctx->AcYF[j]) < ctx->EPSI)
                ctx->AcYF[j] = 0.0;
            tmax = (double)(dmax(ctx, tmax,fabs((double)ctx->AcXF[j] - (double)ctx->AcYF[j])));
            ctx->AcXF[j] = ctx->AcYF[j];
            ctx->AcYF[j] = 0.0;
        }
        printfe(ctx, "Iter %2d  Crit %15.8e\n",iter,tmax);
        if (iter >= mxit || tmax <= tol)
            break;
    }
    printf1(ctx, "Finished after %d iterations. Max change: %g\n",iter,tmax);

    if (ctx->PMProtFDef) {
        fprintf(ctx->PMProtFd,"\nFinal densities after %d iterations.\n",iter);
        tmp = 0.0;
        for (j = 0; j < nred; ++j)  {
            fprintf(ctx->PMProtFd,"%3d  %8.6f\n",j,(double)(ctx->AcXF[j]));
            tmp += (double)(ctx->AcXF[j]);
        }
        fprintf(ctx->PMProtFd,"Sum: %g  Maximal change: %g\n\n",tmp,tmax);
    }

    if (opt <= 1) {

        /*  calculate distribution or survivor function, depending on opt */
        /*  values in val[nu,ndim], function values in f[] */

        if (alloc_acr(ctx, ndim))    
            return(-1);
        if (alloc_acs(ctx, ndim))    
            return(-1);

        for (i = 0; i < nu; ++i) {

            l = i * ndim;
            for (j = 0; j < ndim; ++j) {
                k = ptr[l + j];
                if (yval == NULL)
                    val[l + j] = get_data(ctx, ctx->PMYL,k);
                else
                    val[l + j] = yval[k];

                if (opt == 0) {         /* distribution function */
                    ctx->AcR[j] = ctx->AcM[j] = 1;
                    ctx->AcS[j] = ctx->AcN[l + j];
                }
                else {                  /* survivor function */
                    ctx->AcS[j] = ctx->GDFNBOX;
                    ctx->AcM[j] = ctx->AcN[l + j];
                    /*******************
                    if (AcM[j] < GDFNBOX)
                        AcM[j] += 1;
                    ********************/
                    ctx->AcR[j] = ctx->AcM[j];
                }
            }
            if (ctx->AcNS[i] < 0) {
                f[i] = -1.0;
                continue;
            }
            np = gdf_boxptr(ctx, ndim,ctx->AcM);  /* pointer to first box of observation */

            tmp = 0.0;
            fin = 0;    
            while (fin == 0) {
                if (np >= 0 && np < ctx->GDFTBOX) {
                    k = ctx->AcJ[np];
                    if (k >= 0 && k < nred)
                        tmp += (double)(ctx->AcXF[k]);
                }
                fin = 1;
                for (j = 0; j < ndim; ++j) {
                    if (ctx->AcM[j] < ctx->AcS[j]) {
                        ctx->AcM[j] += 1;
                        for (k = 0; k < j; ++k)     /* reset */
                            ctx->AcM[k] = ctx->AcR[k];
                        fin = 0;
                        break;
                    }
                }
                if (fin == 0)  
                    np = gdf_boxptr(ctx, ndim,ctx->AcM);
            }
            f[i] = tmp;
        }
    }
    else {      /* return expected values in val[], ignore f[] */

        /* first calculate array AcYF[ne,ndim] with means of exact 
           observations in boxes. Also create pointer AcR[0,...,nred-1]
           such that AcR[AcJ[np]] points to the proper row of AcYF[].     
           AcS[ne] counts number of exact observations. */

        if (ne > 0) {           /* number of boxes with exact observations */

            if (alloc_acyf(ctx, ne * ndim))          
                return(-1);
            if (alloc_acs(ctx, ne))          
                return(-1);
            if (alloc_acr(ctx, nred))          
                return(-1);

            for (i = 0; i < nred; ++i)
                ctx->AcR[i] = -1;

            nep = 0;
            for (i = 0; i < nu; ++i) {      /* for all units */

                if (ctx->AcNS[i] != 0)           /* only exact observations inside of grid */
                    continue;   

                l = i * ndim;
                for (j = 0; j < ndim; ++j)   
                    ctx->AcM[j] = ctx->AcN[l + j];

                np = gdf_boxptr(ctx, ndim,ctx->AcM);  /* pointer to first box of observation */

                k = ctx->AcJ[np];
                if (ctx->AcR[k] < 0)  
                    ctx->AcR[k] = nep++;
                kp = ctx->AcR[k];
                for (j = 0; j < ndim; ++j) {
                    if (yval == NULL)
                        ctx->AcYF[kp * ndim + j] = (float)((double)ctx->AcYF[kp * ndim + j] +
                                   get_data(ctx, ctx->PMYL,ptr[l + j]));
                    else
                        ctx->AcYF[kp * ndim + j] = (float)((double)ctx->AcYF[kp * ndim + j] +
                                   (double)yval[ptr[l + j]]);
                }
                ctx->AcS[kp] += 1;
            }
            if (ne != nep) {          
                printfe(ctx, "ERROR in gdf_edf3m: ne=%d nep=%d\n",ne,nep);
                gerr_exit(ctx, 214);
            }
            for (i = 0; i < ne; ++i) {
                if (ctx->AcS[i] > 1) {
                    k = i * ndim;
                    /* Rohwer writes this as AcYF[k++] /= (double)AcS[i],
                       so k advances once per dimension.  Rewriting the
                       compound assignment for the explicit casts moved
                       the k++ out of the subscript and past the loop,
                       which divided AcYF[k] ndim times and AcYF[k+1] not
                       at all: with two dimensions the first box mean came
                       out over the count squared and the second stayed a
                       sum. */
                    for (j = 0; j < ndim; ++j) {
                        ctx->AcYF[k] = (float)((double)ctx->AcYF[k] /
                                               (double)ctx->AcS[i]);
                        k++;
                    }
                }
            }

            if (ctx->PMProtFDef) {
                fprintf(ctx->PMProtFd,"Pointer for mapping exact observations\n");
                for (i = 0; i < nred; ++i)  
                    fprintf(ctx->PMProtFd,"%6d %6d\n",i,ctx->AcR[i]);
                fprintf(ctx->PMProtFd,"Array with mean exact observations\n");
                for (i = 0; i < ne; ++i) {
                    fprintf(ctx->PMProtFd,"%6d %6d ",i,ctx->AcS[i]);
                    for (j = 0; j < ndim; ++j)
                        rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMTFmtS,(double)(ctx->AcYF[i * ndim + j]));
                    fprintf(ctx->PMProtFd,"\n");
                }
            }
        }
        if (ctx->PMProtFDef)
            fprintf(ctx->PMProtFd,"\nUnit  status  values  expectation\n");

        if (alloc_acy1(ctx, ndim))    
            return(-1);
        if (alloc_actmp(ctx, ndim))    
            return(-1);

        tol1 = tol / 10.0;

        for (i = 0; i < nu; ++i) {

            l = i * ndim;

            if (ctx->PMProtFDef) {
                fprintf(ctx->PMProtFd,"%4d  ",i);
                for (j = 0; j < ndim; ++j)
                    fprintf(ctx->PMProtFd,"%2d ",(int)ctx->AcC[l + j]);
            }
            for (j = 0; j < ndim; ++j) {
                k = ptr[l + j];
                if (yval == NULL)
                    val[l + j] = get_data(ctx, ctx->PMYL,k);
                else
                    val[l + j] = yval[k];

                if (ctx->PMProtFDef)   
                    rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMFmtS,val[l + j]);
            }
            if (ctx->AcNS[i] <= 0) {
                if (ctx->PMProtFDef)   
                    fprintf(ctx->PMProtFd,"\n");
                continue;
            }
            for (j = 0; j < ndim; ++j) {               
                ctx->AcM[j] = ctx->AcN[l + j];
                /*************************************** 
                if (AcC[l + j] == 3 && AcM[j] < GDFNBOX)  
                    AcM[j] += 1;
                **************************************/
                ctx->AcTmp[j] = ctx->AcY1[j] = 0.0;
            }
            np = gdf_boxptr(ctx, ndim,ctx->AcM);  /* pointer to first box of observation */

            fin = 0;
            while (fin == 0) {
                if (np >= 0 && np < ctx->GDFTBOX) {
                    k = ctx->AcJ[np];
                    if (k >= 0 && k < nred) {
                        tmp = (double)(ctx->AcXF[k]);

                        if (tmp > tol1) {
                            gdf_getmval(ctx, ne,ndim,np,nred,ctx->AcM,ctx->AcU);

                            for (j = 0; j < ndim; ++j) {
                                if (ctx->AcC[l + j] == 3) {
                                    ctx->AcY1[j] += ctx->AcU[j] * tmp;
                                    ctx->AcTmp[j] += tmp;
                                }   
                            }
                        }
                    }
                }
                fin = 1;
                for (j = 0; j < ndim; ++j) {
                    if (ctx->AcC[l + j] == 3 && ctx->AcM[j] < ctx->GDFNBOX) {
                        ctx->AcM[j] += 1;
                        for (k = 0; k < j; ++k) {    /* reset */
                            ctx->AcM[k] = ctx->AcN[l + k];
                            /*************************************** 
                            if (AcC[l + k] == 3 && AcM[k] < GDFNBOX)  
                                AcM[k] += 1;
                            ***************************************/
                        }
                        fin = 0;
                        break;
                    }
                }
                if (fin == 0)  
                    np = gdf_boxptr(ctx, ndim,ctx->AcM);
            }
            for (j = 0; j < ndim; ++j) {
                if (ctx->AcC[l + j] == 3 && ctx->AcTmp[j] > 0.0)  
                    val[l + j] = dmax(ctx, val[l + j],ctx->AcY1[j] / ctx->AcTmp[j]);
                if (ctx->PMProtFDef)   
                    rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMFmtS,val[l + j]);
            }
            if (ctx->PMProtFDef)   
                fprintf(ctx->PMProtFd,"\n");
        }
        if (ctx->PMProtFDef)   
            fprintf(ctx->PMProtFd,"\n");
    }
    alloc_actmp(ctx, 0);
    alloc_acy1(ctx, 0);
    alloc_acs(ctx, 0);
    alloc_acr(ctx, 0);
    alloc_acyf(ctx, 0);
    alloc_acxf(ctx, 0);
    alloc_acc(ctx, 0);
    alloc_acn(ctx, 0);
    alloc_acns(ctx, 0);
    alloc_acj(ctx, 0);
    alloc_acm(ctx, 0);
    alloc_acu(ctx, 0);
    alloc_ack(ctx, 0);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  gdf_getb(n,ndim,opt)   Calculate number of boxes and set in GDFNBOX     */
/*                         Also calculate GDFBLen.                          */
/*  If opt=1 print to stdout, if opt=2 print to protocol file.              */
   
void gdf_getb(TDAContext *ctx, int n,int ndim,int opt)
{
    register int j; 

    if (n < 1)  
        n = 100;     

    ctx->GDFNBOX = (int)(rexp(ctx, rlog(ctx, (double)n) / (double)ndim) + 0.5);
    if (ctx->GDFNBOX < 1)
        ctx->GDFNBOX = 1;

    ctx->GDFTBOX = (int)pow((double)ctx->GDFNBOX,(double)ndim);

    if (opt == 1)
        printf1(ctx, "Number of boxes in each dimension: %d. Total: %d\n",ctx->GDFNBOX,ctx->GDFTBOX);
    else if (opt == 2 && ctx->PMProtFDef)
        fprintf(ctx->PMProtFd,"Number of boxes in each dimension: %d. Total: %d\n",ctx->GDFNBOX,ctx->GDFTBOX);

    for (j = 0; j < ndim; ++j)    
        ctx->GDFBLen[j] = (ctx->GDFDH[j] - ctx->GDFDL[j]) / (double)ctx->GDFNBOX;
}

/* ------------------------------------------------------------------------ */
/*  gdf_boxptr(ndim,ptr)                                                    */
   
int gdf_boxptr(TDAContext *ctx, int ndim,int *ptr)
{
    register int j,np;

    np = 0;
    for (j = ndim - 1; j >= 0; --j) {
        np *= ctx->GDFNBOX;
        np += (ptr[j] - 1);    
    }
    return(np);
}

/* ------------------------------------------------------------------------ */
/*  gdf_getmval(ne,ndim,np,nred,acm,mval)                                   */
/*                                                                          */
/*  find mean value for box np and return in mval[].                        */
   
void gdf_getmval(TDAContext *ctx, int ne,int ndim,int np,int nred,int *acm,double *mval)
{
    register int k,j; 

    if (ne > 0 && np >= 0 && np < ctx->GDFTBOX) {

        if ((k = ctx->AcJ[np]) >= 0 && k < nred) {
            if ((k = ctx->AcR[k]) >= 0) {
                for (j = 0; j < ndim; ++j)       
                    mval[j] = (double)(ctx->AcYF[k * ndim + j]);
                return;
            }
        }
    }
                                  /* use box means */
    for (j = 0; j < ndim; ++j)  
        mval[j] = ctx->GDFDL[j] + ((double)acm[j] - 0.5) * ctx->GDFBLen[j];
}

/* ------------------------------------------------------------------------ */
/*  gdf_edf4m(t,nu,ndim,ptr,yval,val,f,d,opt,prn)                           */
/*                                                                          */
/*  Local Kaplan-Meier procedure.                                           */
/*                                                                          */
/*  Use pattern t, ptr created by gdf_cptr.                                 */
/*                                                                          */
/*  t = number of current marginal pattern.                                 */
/*  nu = number of units                                                    */
/*  ndim = number of dimensions in current marg pattern                     */
/*  ptr = ptr to observations, created by gdf_cptr().                       */
/*                                                                          */
/*  If yval == NULL use data matrix, otherwise use yval[0,...,NOC-1].       */
/*                                                                          */
/*  If prn != 0 print table with number of i/r-censored observations        */
/*  in the current pattern.                                                 */
/*                                                                          */
/*  If opt = 0  return values in val[] and distribution function in f[]     */
/*  If opt = 1  return values in val[] and survivor function in f[]         */
/*  If opt = 2  return expected values in val[] and ignore f[].             */
/*                                                                          */
/*  Return 0 if OK, or -1 if error.                                         */
   
int gdf_edf4m(TDAContext *ctx, int t,int nu,int ndim,int *ptr,double *yval,double *val, double *f,double d,int opt,int prn)
{
    register int i,j,k,l,ii;
    int err,cen,n,ll,kk,sflg;
    double tmp,tmp0;

    err = -1;
    if (alloc_acns(ctx, nu))             /* index of censored dimension */
        goto EDF4MFin;              /* -1 if not censored */
    if (alloc_acj(ctx, ndim))            /* number of censored observations */
        goto EDF4MFin;              /* in dimension j */

    if (prn)                        /* print table */
        printf1(ctx, "Dimension  i-censored  r-censored  observations\n");

    n = 0;
    for (i = 0; i < nu; ++i) {
        l = i * ndim;
        cen = 0;
        ctx->AcNS[i] = -1;
        for (j = 0; j < ndim; ++j) {
            k = ptr[l + j];
            if (yval == NULL)
                val[l + j] = get_data(ctx, ctx->PMYL,k);
            else
                val[l + j] = yval[k];

            if (ctx->GDFYTyp[k] == 3) {
                ctx->AcNS[i] = (short)(j);  
                ctx->AcJ[j] += 1;
                cen++;   
            }
        }
        if (cen > 1)
            n++;      
    }
    if (prn) {      /* print table with i/r-censored observations */
        for (j = 0; j < ndim; ++j)  
            printf1(ctx, "%6d    %11d %11d  %12d\n",
                             ctx->GDFMarg[t * ctx->GDFNDim + j],0,ctx->AcJ[j],nu);
        newline(ctx);
    }
    if (n > 0) {
        printf1(ctx, "Error: found %d observation(s) having more than a single censored dimension.\n",n);
        goto EDF4MFin;
    }

    /* sort observations */

    if (alloc_acj(ctx, ndim + 1))
        goto EDF4MFin;
    if (alloc_ack(ctx, nu + 1))
        goto EDF4MFin;

    for (j = 1; j <= ndim; ++j)
        ctx->AcJ[j] = j;

    if (sortdpn1(ctx, nu,ndim,val - 1,ndim,ctx->AcJ,ctx->AcK - 1,ctx->AcNS - 1))
        goto EDF4MFin;

    /* calculate densities in AcXF */

    if (alloc_acxf(ctx, nu))
        goto EDF4MFin;
    if (alloc_acu(ctx, nu))              /* used for delta / 2 */
        goto EDF4MFin;
    if (alloc_acr(ctx, nu))      
        goto EDF4MFin;

    for (j = 0; j < ndim; ++j)   
        ctx->AcU[j] = d * (ctx->GDFDH[j] - ctx->GDFDL[j]) / 2.0;

    tmp0 = 1.0 / (double)nu;
    for (i = 0; i < nu; ++i)  
        ctx->AcXF[i] = (float)(tmp0);

    for (i = 0; i < nu; ++i) {
        l = ctx->AcK[i] - 1;
        k = l * ndim;
        cen = ctx->AcNS[l];
        if (cen < 0)
            continue;

        /* find observations for distributing the mass of the current one */

        n = 0;
        for (ii = i + 1; ii < nu; ++ii) {
            ll = ctx->AcK[ii] - 1;
            kk = ll * ndim;
            sflg = 1;
            for (j = 0; j < ndim; ++j) {
                if (j == cen)
                    continue;
                if (ctx->AcNS[ll] >= 0 && ctx->AcNS[ll] != cen) {
                    sflg = 0;
                    break;
                }
                if (fabs(val[kk + j] - val[k + j]) > ctx->AcU[j]) {
                    sflg = 0;
                    break;
                }
            }
            if (sflg)  
                ctx->AcR[n++] = ll;
        }
        if (n == 0)
            continue;

        tmp = (double)(ctx->AcXF[l]) / (double)n;                       
        for (j = 0; j < n; ++j)  
            ctx->AcXF[ctx->AcR[j]] = (float)((double)(ctx->AcXF[ctx->AcR[j]]) + (tmp));
        ctx->AcXF[l] = 0.0;
    }
    if (ctx->PMProtFDef) {
        fprintf(ctx->PMProtFd,"\nDensity after redistribution.\n");
        tmp = 0.0;
        for (i = 0; i < nu; ++i) {
            fprintf(ctx->PMProtFd,"%6d ",i);
            for (j = 0; j < ndim; ++j) 
                rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMTFmtS,val[i * ndim + j]);
            rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMTFmtS,(double)(ctx->AcXF[i]));
            fprintf(ctx->PMProtFd,"\n");
            tmp += (double)(ctx->AcXF[i]);
        }
        fprintf(ctx->PMProtFd,"Sum: %g\n",tmp);
    }

    if (opt == 0) {                         /* distribution function */
        for (i = 0; i < nu; ++i) {
            l = ctx->AcK[i] - 1;
            k = l * ndim;
            tmp = 0.0;
            for (ii = 0; ii < nu; ++ii) {
                ll = ctx->AcK[ii] - 1;
                kk = ll * ndim;
                sflg = 1;
                for (j = 0; j < ndim; ++j) {
                    if (val[kk + j] > val[k + j]) {
                        sflg = 0;
                        break;
                    }   
                }
                if (sflg)   
                    tmp += (double)(ctx->AcXF[ll]);
            }
            f[l] = tmp;
        }
    }
    else if (opt == 1) {                    /* survivor function */
        for (i = 0; i < nu; ++i) {
            l = ctx->AcK[i] - 1;
            k = l * ndim;
            tmp = 0.0;
            for (ii = i + 1; ii < nu; ++ii) {
                ll = ctx->AcK[ii] - 1;
                kk = ll * ndim;
                sflg = 1;
                for (j = 0; j < ndim; ++j) {
                    if (val[kk + j] <= val[k + j]) {
                        sflg = 0;
                        break;
                    }   
                }
                if (sflg)   
                    tmp += (double)(ctx->AcXF[ll]);
            }
            f[l] = tmp;
        }
    }
    else {                                  /* expected values */

        for (i = 0; i < nu; ++i) {
            l = ctx->AcK[i] - 1;
            k = l * ndim;
            cen = ctx->AcNS[l];
            if (cen < 0)
                continue;

            n = 0;
            for (ii = i + 1; ii < nu; ++ii) {
                ll = ctx->AcK[ii] - 1;
                kk = ll * ndim;
                sflg = 1;
                for (j = 0; j < ndim; ++j) {
                    if (j == cen)
                        continue;
                    if (ctx->AcNS[ll] >= 0 && ctx->AcNS[ll] != cen) {
                        sflg = 0;
                        break;
                    }
                    if (fabs(val[kk + j] - val[k + j]) > ctx->AcU[j]) {
                        sflg = 0;
                        break;
                    }
                }
                if (sflg)  
                    ctx->AcR[n++] = ll;
            }
            if (n == 0)
                continue;

            tmp0 = tmp = 0.0;
            for (j = 0; j < n; ++j) {
                ll = ctx->AcR[j];
                tmp += (double)((val[ll * ndim + cen])) * (double)(ctx->AcXF[ll]);
                tmp0 += (double)(ctx->AcXF[ll]);
            }
            if (tmp0 > 0.0)
                val[k + cen] = tmp / tmp0;
        }   
    } 
    err = 0;

EDF4MFin:
    alloc_acr(ctx, 0);
    alloc_acu(ctx, 0);
    alloc_acxf(ctx, 0);
    alloc_ack(ctx, 0);
    alloc_acj(ctx, 0);
    alloc_acns(ctx, 0);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  lsreg1()        Regression with censored data.                          */
/*  ##                                                                      */
/*                                                                          */
/*                  lsreg1(                                                 */
/*                      opt=...,        option, default 1                   */
/*                                      1 = marginal df                     */
/*                                      2 = joint df, method 1              */
/*                                      3 = joint df, method 2              */
/*                      yl=...,         dep. variable                       */
/*                      cen=...,        censoring indicator                 */
/*                      grp=ID,L1,      defines hierarchical structure      */
/*                      ni=...,         ni=1 if without intercept, def. 0   */
/*                      mxit=...,       max iterations, def. 20             */
/*                      tolp=...,       tolerance for convergence, 0.001    */
/*                      n=...,          number of boxes, def. 100           */
/*                      sc=...,         offset for domain, def. 0.0         */
/*                      mxitl=...,      max iterations for survivor         */
/*                                      function estimation, def. 10        */     
/*                      tolf=...,       tolerance for convergence, 0.001    */
/*                      d=...,          specification of delta, def. 0.1    */
/*                      prot=...,       diagnostic information              */
/*                      tfmt=...,       print format for parameters         */
/*                      ppar=...,       print estimated parameters          */
/*                      pres=...,       write data and residuals            */
/*                      fmt = ...,      print format, def. 10.4             */
/*                      dtda =...,      write TDA description               */
/*                                                                          */
/*                  ) [= varlist]; (independent variables, optional)        */
/*                                                                          */
/*  grp = ID,L1,                                                            */
/*                                                                          */
/*  If this parameter is not used, each data matrix row is treated as a     */
/*  separate unit. Otherwise, ID is an ID variable that defines units,      */
/*  and L1 is a variables that defines the dimension (level 1) to which     */
/*  the corresponding data matrix row belongs. Values of L1 must be         */
/*  positive integers. It is not required that these values are contiguous. */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int lsreg1(TDAContext *ctx)
{
    register int i;
    int err,nx,r,iter;
    double tmp = 0.0;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Regression with censored data. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->MxItFlg = 0;
    ctx->MxIt1 = 10;
    ctx->TOLP = 0.001;       /* tolerance for regression parameters */
    ctx->TOLF = 0.001;       /* tolerance for surv function estimation */

    if (parm(ctx, ctx->CmdBuf + 6,4,0))     /* get parameters */
        goto LS1Fin;

    if (ctx->MxItFlg == 0)
        ctx->MxIter = 20;
    else if (ctx->MxIter < 1)
        ctx->MxIter = 1;
    if (ctx->MxIt1 < 1)
        ctx->MxIt1 = 1;
    if (ctx->PMD <= 0.0)
        ctx->PMD = 0.1;
    else if (ctx->PMD > 1.0)
        ctx->PMD = 1.0;
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    newline(ctx);
    if (ctx->PMNI != 0) {
        if (ctx->PMNV == 0) {
            printf1(ctx, "Error: need intercept or independent variables.\n");
            goto LS1Fin;
        }
        ctx->PMNI = 1;
    }
    if (ctx->PMYL < 0 || ctx->PMCEN < 0) {
        printf1(ctx, "Error: need both, yl and cen parameters.\n");     
        goto LS1Fin;
    }

    if (gdf_dcheck(ctx))           /* check dependent variable */
        goto LS1Fin;            /* creates GDFYTyp: 1 exact,
                                                    2 i-censored,
                                                    3 r-censored  */
    if (ctx->PMNV > 0) {
        printf1(ctx, "Indep. variables: %s",ctx->VName[ctx->PMVIdx[0]]);
        for (i = 1; i < ctx->PMNV; ++i)  
            printf1(ctx, ",%s",ctx->VName[ctx->PMVIdx[i]]);
        newline(ctx);
    }
    else
        printf1(ctx, "No independent variables.\n");
    newline(ctx);
        
    if (gdf_gcheck(ctx))           /* check group structure */
        goto LS1Fin;

    if (ctx->PMOPT == 1)  
        printf1(ctx, "Option 1: marginal estimation");
    else {
        if (ctx->PMOPT > 3)
            ctx->PMOPT = 3;
        printf1(ctx, "Option %d: joint estimation (method %d)",ctx->PMOPT,ctx->PMOPT - 1);
    }
    printf1(ctx, " of conditional expectations.\n");
    if (ctx->PMOPT > 1 && ctx->GDFGRP == 0) {
        printf1(ctx, "Error: grp parameter required.\n");
        goto LS1Fin;
    }

    printf1(ctx, "Max number of iterations: %d\n",ctx->MxIter);
    printf1(ctx, "Tolerance for convergence: %g\n",ctx->TOLP);

    nx = ctx->PMNV;
    if (ctx->PMNI == 0)
        nx++;

    if (alloc_acy(ctx, ctx->NOC))                 /* Y vector */
        goto LS1Fin;
    if (alloc_acz(ctx, imax(ctx, ctx->NOC,nx)))        /* residuals */
        goto LS1Fin;
    if (alloc_acx(ctx, nx * nx + 1))         /* X'X */
        goto LS1Fin;
    if (alloc_acv(ctx, nx))                  /* beta */
        goto LS1Fin;
    if (alloc_acw(ctx, nx))                  /* new beta */
        goto LS1Fin;

    ls1_xx(ctx, nx,1 - ctx->PMNI);        /* create X'X */
    r = ginv(ctx, nx,nx,ctx->AcX);        /* generalized inverse */
    if (r < nx) {
        if (r < 0)
            p_err(ctx, -2,1);
        else  
            printf1(ctx, "Error: X matrix is rank deficient, rank = %d\n",r);
        goto LS1Fin;
    }

    /* start with YL */

    for (i = 0; i < ctx->NOC; ++i)
        ctx->AcY[i] = get_data(ctx, ctx->PMYL,i);

    /*  get initial parameters into AcV */

    ls1_gpar(ctx, ctx->NOC,ctx->AcY,nx,1 - ctx->PMNI,ctx->AcX,ctx->AcV,ctx->AcZ);  

    if (ctx->PMProtFDef) {
        fprintf(ctx->PMProtFd,"Regression with censored data.\n");
        prvec(ctx, "Initial parameters",nx,ctx->AcV - 1);
    }
    for (iter = 1; iter <= ctx->MxIter; ++iter) {

        /* calculate predicted values in AcY, residuals in AcZ */

        ls1_res(ctx, ctx->NOC,nx,1 - ctx->PMNI,ctx->AcV,ctx->AcY,ctx->AcZ);  

        if (ctx->PMOPT == 1) {                       /* marginal estimates */
            if (ls1_marg(ctx, ctx->NOC,ctx->AcZ,ctx->GDFYTyp))
                goto LS1Fin;
        }
        else {                                  /* joint estimates */
            if (ls1_joint(ctx, ctx->NOC,ctx->AcZ,ctx->GDFYTyp,ctx->PMOPT))
                goto LS1Fin;
        }

        /*  calculate new Y vector */

        for (i = 0; i < ctx->NOC; ++i) {
            if (ctx->GDFYTyp[i] == 3)
                ctx->AcY[i] += ctx->AcZ[i];              
            else
                ctx->AcY[i] = get_data(ctx, ctx->PMYL,i);
        }

        /*  get new parameters into AcW */

        ls1_gpar(ctx, ctx->NOC,ctx->AcY,nx,1 - ctx->PMNI,ctx->AcX,ctx->AcW,ctx->AcZ);  
        /* check convergence and copy new parameter estimates into AcV */

        tmp = 0.0;  
        for (i = 0; i < nx; ++i) {
            tmp = dmax(ctx, tmp,fabs(ctx->AcW[i] - ctx->AcV[i]) / dmax(ctx, ctx->AcW[i],1.0));
            ctx->AcV[i] = ctx->AcW[i];
        }
        if (ctx->PMProtFDef) {
            fprintf(ctx->PMProtFd,"\nIteration: %d\n",iter);
            prvec(ctx, "Current parameters",nx,ctx->AcW - 1);
            prval(ctx, "Criterion",tmp);
        }
        if (ctx->SILENTFlg < 2) {
            printfe(ctx, "Iteration%3d  Crit: %20.13e\n",iter,tmp);
            fflushe(ctx);
        }   
        if (tmp < ctx->TOLP)
            break;
    }
    printf1(ctx, "\nFinished after %d iterations.\n",iter);
    printf1(ctx, "Final maximal parameter change: %g\n\n",tmp);     


    if (ctx->PMPPFDef)               /* write parameter to output file */
        ls1_ppar(ctx, nx,ctx->AcV);
      
    if (ctx->PMResFDef) {            /* write data and residuals to output file */
        ls1_res(ctx, ctx->NOC,nx,1 - ctx->PMNI,ctx->AcV,ctx->AcY,ctx->AcZ);  
        ls1_pres(ctx, nx,1 - ctx->PMNI,ctx->AcY,ctx->AcZ);
    }
    err = 0;

LS1Fin:
    alloc_gdfytyp(ctx, 0);
    alloc_gdfgen(ctx, 0,0);
    if (ctx->PM1NV > 0)                     
        vsort(ctx, 0,ctx->PM1VIdx,0,0,0);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ls1_xx(nx,iflag)   Create X'X in AcX                                    */
/*  iflag = 1 if with intercept, otherwise iflag = 0.                       */ 

void ls1_xx(TDAContext *ctx, int nx,int iflag)
{
    register int i,j,k;
    double tmp,tmp1,tmp2;

    for (i = 0; i < nx; ++i) {
        for (j = 0; j < nx; ++j) {

            tmp = 0.0;
            for (k = 0; k < ctx->NOC; ++k) {

                if (i == 0 && iflag)
                    tmp1 = 1.0;
                else
                    tmp1 = get_data(ctx, ctx->PMVIdx[i - iflag],k);

                if (j == 0 && iflag)
                    tmp2 = 1.0;
                else
                    tmp2 = get_data(ctx, ctx->PMVIdx[j - iflag],k);
                tmp += tmp1 * tmp2;
            }
            ctx->AcX[i * nx + j + 1] = tmp;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  ls1_gpar(n,y,nx,iflag,xx,b,h)                                           */
/*                                                                          */
/*  calculate: b = xx X' y  (new parameter vector: b{0,...,nx-1].           */
/*                                                                          */
/*  It is assumed that n = NOC, iflag = 1 if with intercept, otherwise      */
/*  iflag = 0. xx is (nx,nx) matrix containing (X'X)(-1).                   */
/*  Y and X are taken directly from data matrix. h[] is used for            */
/*  intermediate calculations.                                              */

void ls1_gpar(TDAContext *ctx, int n,double *y,int nx,int iflag,double *xx,double *b,double *h)
{
    register int i,j;
    double tmp,tmp1;

    for (j = 0; j < nx; ++j) {
        tmp = 0.0;
        for (i = 0; i < n; ++i) {
            if (j == 0 && iflag)
                tmp1 = 1.0;
            else
                tmp1 = get_data(ctx, ctx->PMVIdx[j - iflag],i);
            tmp += y[i] * tmp1;
        }
        h[j] = tmp;
    }
    for (i = 0; i < nx; ++i) {
        tmp = 0.0;
        for (j = 0; j < nx; ++j)   
            tmp += xx[i * nx + j + 1] * h[j];
        b[i] = tmp;
    }
}

/* ------------------------------------------------------------------------ */
/*  ls1_res(n,nx,iflag,b,yp,res)                                            */
/*                                                                          */
/*  calculate residuals in res[], predicted values in yp[].                 */

void ls1_res(TDAContext *ctx, int n,int nx,int iflag,double *b,double *yp,double *res)
{
    register int i,j;
    double tmp,tmp1;

    for (i = 0; i < n; ++i) {
        tmp = 0.0;
        for (j = 0; j < nx; ++j) {
            if (j == 0 && iflag)
                tmp1 = 1.0;
            else 
                tmp1 = get_data(ctx, ctx->PMVIdx[j - iflag],i);
            tmp += tmp1 * b[j];
        }
        yp[i] = tmp;
        res[i] = get_data(ctx, ctx->PMYL,i) - tmp;
    }
}

/* ------------------------------------------------------------------------ */
/*  ls1_ppar(nx,beta)       Write parameter to PMPPFd.                      */

void ls1_ppar(TDAContext *ctx, int nx,double *beta)
{
    register int i;

    for (i = 0; i < nx; ++i) {
        rt_fprintf_d(ctx, ctx->PMPPFd,ctx->PMTFmtS,beta[i]);
        fprintf(ctx->PMPPFd,"\n");
    }
    printf1(ctx, "Parameter estimates written to: %s\n",ctx->PMPPFName);
}

/* ------------------------------------------------------------------------ */
/*  ls1_pres(nx,iflag,yp,res)                                               */
/*  Write data and residuals to PMResFd.                                    */

void ls1_pres(TDAContext *ctx, int nx,int iflag,double *yp,double *res)
{
    register int i,j;
    double tmp;

    for (i = 0; i < ctx->NOC; ++i) {
        fprintf(ctx->PMResFd,"%6d ",i + 1);
        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMFmtS,get_data(ctx, ctx->PMYL,i));
        if (ctx->PMYH >= 0)
            rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMFmtS,get_data(ctx, ctx->PMYH,i));
        if (ctx->PMCEN >= 0)
            rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMFmtS,get_data(ctx, ctx->PMCEN,i));
        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMFmtS,yp[i]);
        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMFmtS,res[i]);

        for (j = 0; j < nx; ++j) {
            if (j == 0 && iflag)
                tmp = 1.0;
            else
                tmp = get_data(ctx, ctx->PMVIdx[j - iflag],i);
            rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMFmtS,tmp);
        }
        fprintf(ctx->PMResFd,"\n");
    }
    printf1(ctx, "Data and residuals written to: %s\n",ctx->PMResFName);

    if (ctx->PMTDAFDef) {
        fprintf(ctx->PMTDAFd,"# data written by lsreg1 command.\n");
        fprintf(ctx->PMTDAFd,"nvar(\n");
        fprintf(ctx->PMTDAFd,"  dfile = %s,\n",ctx->PMResFName);
        fprintf(ctx->PMTDAFd,"  noc = %d,\n",ctx->NOC);
        j = 0;
        fprintf(ctx->PMTDAFd,"  CASE [6.0] = c%-2d,\n",++j);
        fprintf(ctx->PMTDAFd,"  YL   [%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,++j);
        if (ctx->PMYH >= 0)
            fprintf(ctx->PMTDAFd,"  YH   [%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,++j);
        if (ctx->PMCEN >= 0)
            fprintf(ctx->PMTDAFd,"  CEN  [%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,++j);
        fprintf(ctx->PMTDAFd,"  YPRED[%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,++j);
        fprintf(ctx->PMTDAFd,"  RES  [%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,++j);

        if (iflag)
            fprintf(ctx->PMTDAFd,"  INT  [%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,++j);
        for (i = 0; i < ctx->PMNV; ++i)  
            fprintf(ctx->PMTDAFd,"  %s [%d.%d] = c%-2d,\n",ctx->VName[ctx->PMVIdx[i]],ctx->PMFmt1,ctx->PMFmt2,++j);
        fprintf(ctx->PMTDAFd,");\n");
        printf1(ctx, "TDA description written to: %s\n",ctx->PMTDAFName);
    }
}

/* ------------------------------------------------------------------------ */
/*  ls1_marg(n,res,ytyp)                                                    */
/*                                                                          */
/*  res[i], ytyp[i] (i = 0,...,n-1) are the input data.                     */
/*  ytyp[i] = 1  if exact, 2 if i-censored, 3 if r-censored.                */
/*                                                                          */
/*  The function returns res[i] = 0 if not censored, otherwise the          */
/*  conditional expectation.                                                */
/*  The largest observation is always assumed to be uncensored.             */
/*                                                                          */
/*  Return 0 if OK, or -1 if error.                                         */

int ls1_marg(TDAContext *ctx, int n,double *res,short *ytyp)
{
    register int i,j,k;
    int m;
    double tmp,tmp1,fsum,wsum;

/**
tda_out("res vor marg\n");
for (i = 0; i < n; ++i)
tda_out("i=%4d res=%g\n",i,res[i]);
**/


    if (alloc_ack(ctx, n))
        return(-1);
    if (alloc_actmp(ctx, n))
        return(-1);
           
    if (sortdp2a(ctx, n,res,ytyp,ctx->AcK,1))      /* sort */
        return(-1);    

    tmp = 1.0 / (double)n;
    m = n - 1;
    for (i = 0; i < n; ++i) {
        j = ctx->AcK[i];
        ctx->AcTmp[j] += tmp;
        if (ytyp[j] == 3 && i < n - 1) {
            tmp1 = ctx->AcTmp[j] / (double)m;
            for (k = i + 1; k < n; ++k)  
                ctx->AcTmp[ctx->AcK[k]] += tmp1;
            ctx->AcTmp[j] = 0.0;
        }
        m--;
    }
    fsum = wsum = 0.0;
    for (i = n - 1; i >= 0; --i) {
        j = ctx->AcK[i];
        if (ytyp[j] == 1 || i == n - 1) {       /* if not censored */
            fsum += res[j] * ctx->AcTmp[j];
            wsum += ctx->AcTmp[j];
            if (ytyp[j] == 1)
                res[j] = 0.0;
        }
        else if (ytyp[j] == 3) {        /* if r-censored */
            if (wsum > 0.0)
                res[j] = fsum / wsum;
        }
    }
/**
tda_out("res in marg\n");
for (i = 0; i < n; ++i)
tda_out("i=%4d res=%g\n",i,res[i]);
**/
    alloc_ack(ctx, 0);
    alloc_actmp(ctx, 0);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  ls1_joint(n,res,ytyp,opt)                                               */
/*                                                                          */
/*  res[i], ytyp[i] (i = 0,...,n-1) are the input data.                     */
/*  ytyp[i] = 1  if exact, 3 if r-censored.                                 */
/*                                                                          */
/*  The function returns res[i] = 0 if not censored, otherwise the          */
/*  conditional expectation.                                                */
/*  The largest observation is always assumed to be uncensored.             */
/*                                                                          */
/*  Return 0 if OK, or -1 if error.                                         */

int ls1_joint(TDAContext *ctx, int n,double *res,short *ytyp,int opt)
{
    register int i,j,k,ii;
    int err,t,nu,ndim;
    double tmp;

    err = -1;
/**
    printf1(ctx, "noc=%d res at begin\n",n);
    for (i = 0; i < n; ++i) {
        printf1(ctx, "i=%3d res=%g\n",i,res[i] );
    }
**/


    for (t = 0; t < ctx->GDFNTyp; ++t) {     /* for all patterns */

        nu = ctx->GDFMUCnt[t];
        ndim = ctx->GDFMLen[t];

        /**  printf1("Pattern t=%d nu=%d ndim=%d\n",t,nu,ndim);  **/


        if (ctx->PMProtFDef) {
            fprintf(ctx->PMProtFd,"\nCurrent marginal pattern: %d  Typ: %d\n",t+1,ctx->GDFMTyp[t]);
            fprintf(ctx->PMProtFd,"Number of units: %d  dimensions: %d\n",nu,ndim);
        }   
        if (alloc_aci(ctx, nu * ndim + 1))   /* pointer for cases */
            goto LS1JFin;

        gdf_cptr(ctx, t,nu,ndim,ctx->AcI);        /* create pointer to data matrix cases */

        if (ctx->GDFMTyp[t] != 1) {          /* exact and censored observations */

            if (gdf_domain(ctx, t,nu,ndim,ctx->AcI,res,2,0.0))
                goto LS1JFin;

            gdf_getb(ctx, ctx->PMN,ndim,2);                   /* get number of boxes */
    
            if (alloc_act(ctx, nu * ndim + 1))           /* values/lower bounds */
                goto LS1JFin;
   
            if (opt == 2) {
                if (gdf_edf3m(ctx, t,nu,ndim,ctx->AcI,res,ctx->AcT,&tmp,2,1,ctx->MxIt1,ctx->TOLF))
                    goto LS1JFin;
            }
            else  {
                if (gdf_edf4m(ctx, t,nu,ndim,ctx->AcI,res,ctx->AcT,&tmp,ctx->PMD,2,1))
                    goto LS1JFin;
            }

            /* replace values in res[] with cond. expectations */

            for (i = 0; i < nu; ++i) {
                ii = i * ndim;
                for (j = 0; j < ndim; ++j)  {
                    k = ctx->AcI[ii + j]; 
                    res[k] = ctx->AcT[ii + j];
                }
            }
        }
    }

    /* set res[i] = 0 for exact observations */

    for (i = 0; i < n; ++i) {
        if (ytyp[i] == 1)
            res[i] = 0.0;          
    }
/**
    printf1(ctx, "noc=%d final\n",n);
    for (i = 0; i < n; ++i) {
        printf1(ctx, "i=%3d res=%g\n",i,res[i] );
    }
**/

    err = 0; 

LS1JFin:
    alloc_act(ctx, 0);
    alloc_aci(ctx, 0);
    return(err);
}

