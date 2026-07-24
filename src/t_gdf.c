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

/*  functions in t_gdf.c */

int gdf(void);
void gdf_pmet(int met,int opt);
int gdf_dcheck(void);
int alloc_gdfytyp(int n);
int alloc_gdfgen(int nu,int ndim);
int gdf_gcheck(void);
void gdf_pmtyp(void);
int gdf_marg(int opt);
void gdf_pdat1(int d,int n,double *y,double *f);
void gdf_pdat1a(int d,int n,double *y,double *f,short *cen,int opt);
int gdf_edf1(int n,double *y,double *f,int opt);
int gdf_edf2(int n,double *y,short *cen,double *y1,double *h,int opt);
int gdf_joint(int met,int opt);
void gdf_pdat2(int t,int m,int ndim,double *val,double *f);
void gdf_pdat2a(int t,int m,int ndim,double *val,int *ptr,int opt);
void gdf_cptr(int t,int nu,int ndim,int *ptr);
int gdf_domain(int t,int nu,int ndim,int *ptr,double *yval,int opt,double d);
int gdf_edf1m(int t,int nu,int ndim,int *ptr,double *val,double *f,int opt);
int gdf_edf3m(int t,int nu,int ndim,int *ptr,double *yval,double *val,
    double *f,int opt,int prn,int mxit,double tol);
int gdf_edf4m(int t,int nu,int ndim,int *ptr,double *yval,double *val,
    double *f,double d,int opt,int prn);                     
void gdf_getb(int n,int ndim,int opt);
int gdf_boxptr(int ndim,int *ptr);
void gdf_getmval(int ne,int ndim,int np,int nred,int *acm,double *mval);

int lsreg1(void);
void ls1_xx(int nx,int iflag);
void ls1_gpar(int n,double *y,int nx,int iflag,double *xx,double *b,double *h);
void ls1_res(int n,int nx,int iflag,double *b,double *yp,double *res);
void ls1_ppar(int nx,double *beta);
void ls1_pres(int nx,int iflag,double *yp,double *res);
int ls1_marg(int n,double *res,short *ytyp);
int ls1_joint(int n,double *res,short *ytyp,int opt);


int GDFTyp = 0;
int GDFWNOC = 0;            /* records written to output file               */
int GDFNU = 0;              /* number of units                              */
int GDFNDim = 0;            /* number of dimensions                         */
int GDFGRP = 0;             /* set if with groups                           */
int GDFNEX = 0;             /* number of exact cases                        */
int GDFNINT = 0;            /* number of interval censored cases            */
int GDFNCEN = 0;            /* number of right censored cases               */
int GDFNBOX = 0;            /* number of boxes in each dimension            */
int GDFTBOX = 0;            /* total number of boxes                        */

int GDFNTyp = 0;            /* number of m-types                            */
short *GDFYTyp;             /* type of dependent variable                   */
int GDFYTypA = 0;           /* if allocated                                 */
int *GDFDim;                /* maps L1-values to dimensions                 */
int GDFDimA = 0;            /* if allocated                                 */
int *GDFDCnt;               /* number of observations in dimension          */
int GDFDCntA = 0;           /* if allocated                                 */
int *GDFDInt;               /* number of interval observations in dim.      */
int GDFDIntA = 0;           /* if allocated                                 */
int *GDFDCen;               /* number of censored observations in dim.      */
int GDFDCenA = 0;           /* if allocated                                 */
int *GDFPtr;                /* pointer for units                            */
int GDFPtrA = 0;            /* if allocated                                 */
short *GDFULen;             /* number of observations for unit              */
int GDFULenA = 0;           /* if allocated                                 */
short *GDFMarg;             /* list of m-types                              */
int GDFMargA = 0;           /* if allocated                                 */
short *GDFMLen;             /* length of m-type                             */
int GDFMLenA = 0;           /* if allocated                                 */
int *GDFMUCnt;              /* number of units for m-type                   */
int GDFMUCntA = 0;          /* if allocated                                 */
int *GDFMTyp;               /* type of pattern                              */
int GDFMTypA = 0;           /* if allocated                                 */
double *GDFDL;              /* domain, lower bound                          */
int GDFDLA = 0;             /* if allocated                                 */
double *GDFDH;              /* domain, upper bound                          */
int GDFDHA = 0;             /* if allocated                                 */
double *GDFBLen;            /* length of boxes                              */
int GDFBLenA = 0;           /* if allocated                                 */

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

int gdf(void)
{
    int err,opt;

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("General distribution functions. Current memory: %d bytes.\n",MemReq);

    MxItFlg = 0;
    TOLF = 0.001; 
    if (parm(CmdBuf + 3,1,1))     /* get parameters */
        goto GDFFin;

    if (PMYH >= 0) {
        printf1("Error: yh parameter cannot be used with gdf command.\n");
        goto GDFFin;
    }
    if (MxItFlg == 0)
        MxIter = 20;
    else if (MxIter < 0)
        MxIter = 0;
    if (PMOPT > 3)
        PMOPT = 3;
    if (PMPRNO > 2)
        PMPRNO = 2;
    if (PMD <= 0.0)
        PMD = 0.1;
    else if (PMD > 1.0)
        PMD = 1.0;
    if (PMFmtF == 0)
        pmfmt(10,4);

    newline();
    if (gdf_dcheck())           /* check dependent variable */
        goto GDFFin;

    if (gdf_gcheck())           /* check group structure */
        goto GDFFin;

    GDFWNOC = 0;                /* records written to output file */
    opt = 0;

    if (PMOPT == 1) {           /* marginal distributions */
        if (gdf_marg(PMPRNO))
            goto GDFFin;
    }
    else {                      /* joint distributions */
        if (gdf_joint(PMOPT,PMPRNO))
            goto GDFFin;
    }
    if (GDFWNOC > 0) 
        printf1("\n%d records written to: %s\n",GDFWNOC,PMFdName);
    err = 0;

GDFFin:
    alloc_gdfytyp(0);
    alloc_gdfgen(0,0);
    if (PM1NV > 0)                     
        vsort(0,PM1VIdx,0,0,0);
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdf_pmet(met,opt)      Print method and option.                         */

void gdf_pmet(int met,int opt)
{
    if (met == 1)
        printf1("Marginal calculation:");
    else if (met == 2)
        printf1("Joint calculation (method 1):");
    else                   
        printf1("Joint calculation (method 2):");
    if (opt == 0)
        printf1(" distribution functions.\n\n");
    else if (opt == 1)
        printf1(" survivor functions.\n\n");
    else
        printf1(" expected values.\n\n");
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

int gdf_dcheck(void)
{
    register int i;
    double yl,yh;    

    GDFNEX = GDFNINT = GDFNCEN = 0;

    if (PMYL < 0) {
        printf1("Error: yl parameter is required.\n");
        return(-1);
    }
    printf1("yl: %s",VName[PMYL]);
    if (PMYH >= 0)
        printf1("  yh: %s",VName[PMYH]);
    if (PMCEN >= 0)
        printf1("  cen: %s",VName[PMCEN]);
    newline();

    if (alloc_gdfytyp(NOC))
        return(-1);   

    for (i = 0; i < NOC; ++i) {
        GDFYTyp[i] = 1;
        yl = get_data(PMYL,i);
        if (PMYH >= 0) {
            yh = get_data(PMYH,i);
            if (yh > yl)
                GDFYTyp[i] = 2;
            else if (yh < yl) {
                printf1("Error: no valid interval in case %d.\n",i + 1);
                return(-1);   
            }
        }
        if (PMCEN >= 0) {
            if (GDFYTyp[i] == 1 && fabs(get_data(PMCEN,i)) <= EPSI1)
                GDFYTyp[i] = 3;
        }
        if (GDFYTyp[i] == 1)
            GDFNEX++;
        else if (GDFYTyp[i] == 2)
            GDFNINT++;
        else
            GDFNCEN++;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_gdfytyp(n)    If n > 0 allocate GDFYTyp, otherwise free.          */
/*                      Return 0 if OK, -1 if error.                        */

int alloc_gdfytyp(int n)
{
    if (GDFYTypA > 0) {
        free((char *)GDFYTyp);
        memrq(-GDFYTypA,sizeof(short));
        GDFYTypA = 0;
    }
    if (n > 0) {
        if (!(GDFYTyp = (short *)calloc(n,sizeof(short)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(short));
        GDFYTypA = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_gdfgen(nu,ndim)   if nu and ndim > 0 allocate basic data          */
/*                          structures, otherwise free.                     */
/*                          Return 0 if OK, -1 if error.                    */

int alloc_gdfgen(int nu,int ndim)
{
    int err = -1;

    if (nu == 0 && ndim == 0) {
        err = 0;
        goto GDFGFin;
    }
    if (ndim > 0) {
        if (!(GDFDim = (int *)calloc(ndim,sizeof(int))))  
            goto GDFGFin; 
        memrq(ndim,sizeof(int));
        GDFDimA = ndim;

        if (!(GDFDCnt = (int *)calloc(ndim,sizeof(int))))  
            goto GDFGFin; 
        memrq(ndim,sizeof(int));
        GDFDCntA = ndim;

        if (!(GDFDInt = (int *)calloc(ndim,sizeof(int))))  
            goto GDFGFin; 
        memrq(ndim,sizeof(int));
        GDFDIntA = ndim;

        if (!(GDFDCen = (int *)calloc(ndim,sizeof(int))))  
            goto GDFGFin; 
        memrq(ndim,sizeof(int));
        GDFDCenA = ndim;

        if (!(GDFDL = (double *)calloc(ndim,sizeof(double))))  
            goto GDFGFin; 
        memrq(ndim,sizeof(double));
        GDFDLA = ndim;

        if (!(GDFDH = (double *)calloc(ndim,sizeof(double))))  
            goto GDFGFin; 
        memrq(ndim,sizeof(double));
        GDFDHA = ndim;

        if (!(GDFBLen = (double *)calloc(ndim,sizeof(double))))  
            goto GDFGFin; 
        memrq(ndim,sizeof(double));
        GDFBLenA = ndim;
    }
    if (nu > 0) {
        if (!(GDFPtr = (int *)calloc(nu,sizeof(int))))  
            goto GDFGFin;
        memrq(nu,sizeof(int));
        GDFPtrA = nu;

        if (!(GDFULen = (short *)calloc(nu,sizeof(short))))  
            goto GDFGFin;
        memrq(nu,sizeof(short));
        GDFULenA = nu;
    }
    if (nu > 0 && ndim > 0) {
        if (!(GDFMarg = (short *)calloc(nu * ndim + 1,sizeof(short))))  
            goto GDFGFin;
        memrq(nu * ndim + 1,sizeof(short));
        GDFMargA = nu * ndim + 1;

        if (!(GDFMLen = (short *)calloc(nu,sizeof(short))))  
            goto GDFGFin;
        memrq(nu,sizeof(short));
        GDFMLenA = nu;

        if (!(GDFMUCnt = (int *)calloc(nu,sizeof(int))))  
            goto GDFGFin;
        memrq(nu,sizeof(int));
        GDFMUCntA = nu;

        if (!(GDFMTyp = (int *)calloc(nu,sizeof(int))))  
            goto GDFGFin;
        memrq(nu,sizeof(int));
        GDFMTypA = nu;
    }
    return(0);

GDFGFin:
    if (err)
        p_err(-2,1);

    if (GDFDimA > 0) {
        free((char *)GDFDim);
        memrq(-GDFDimA,sizeof(int));
        GDFDimA = 0;
    }
    if (GDFDCntA > 0) {
        free((char *)GDFDCnt);
        memrq(-GDFDCntA,sizeof(int));
        GDFDCntA = 0;
    }
    if (GDFDIntA > 0) {
        free((char *)GDFDInt);
        memrq(-GDFDIntA,sizeof(int));
        GDFDIntA = 0;
    }
    if (GDFDCenA > 0) {
        free((char *)GDFDCen);
        memrq(-GDFDCenA,sizeof(int));
        GDFDCenA = 0;
    }
    if (GDFDLA > 0) {
        free((char *)GDFDL);
        memrq(-GDFDLA,sizeof(double));
        GDFDLA = 0;
    }
    if (GDFDHA > 0) {
        free((char *)GDFDH);
        memrq(-GDFDHA,sizeof(double));
        GDFDHA = 0;
    }
    if (GDFBLenA > 0) {
        free((char *)GDFBLen);
        memrq(-GDFBLenA,sizeof(double));
        GDFBLenA = 0;
    }
    if (GDFPtrA > 0) {
        free((char *)GDFPtr);
        memrq(-GDFPtrA,sizeof(int));
        GDFPtrA = 0;
    }
    if (GDFULenA > 0) {
        free((char *)GDFULen);
        memrq(-GDFULenA,sizeof(short));
        GDFULenA = 0;
    }
    if (GDFMargA > 0) {
        free((char *)GDFMarg);
        memrq(-GDFMargA,sizeof(short));
        GDFMargA = 0;
    }
    if (GDFMLenA > 0) {
        free((char *)GDFMLen);
        memrq(-GDFMLenA,sizeof(short));
        GDFMLenA = 0;
    }
    if (GDFMUCntA > 0) {
        free((char *)GDFMUCnt);
        memrq(-GDFMUCntA,sizeof(int));
        GDFMUCntA = 0;
    }
    if (GDFMTypA > 0) {
        free((char *)GDFMTyp);
        memrq(-GDFMTypA,sizeof(int));
        GDFMTypA = 0;
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

int gdf_gcheck(void)
{
    register int i,j,k,l,ii;
    int err,n,i0,l1,iid,il1,t,tt,fnd,cen2,cen3;
    double id,id0;

    err = -1;
    if (PM1NV == 0) {           /* no group structure */

        GDFGRP = 0;
        GDFNU = NOC;
        GDFNDim = 1;

        if (!(GDFDim = (int *)calloc(1,sizeof(int)))) {
            p_err(-2,1);
            goto GDFCFin; 
        }
        memrq(1,sizeof(int));
        GDFDimA = 1;     
        GDFDim[0] = 1;

        err = 0;
        goto GDFCON;
    }
    if (PM1NV < 2) {
        printf1("Error: grp parameter requires at least two variables.\n");
        goto GDFCFin;
    }
    if (vsort(PM1NV,PM1VIdx,1,1,0))     /* sort, create ptr VSORTPtr[] */
        goto GDFCFin;

    /* check units and observations */

    if (alloc_acm(NOC + 1))     /* used to sort L1 values */
        goto GDFCFin;

    ii = l = 0;
    GDFNU = 0;
    j = PM1VIdx[0];
    k = PM1VIdx[1];
    id0 = get_data(j,VSORTPtr[0]) - 1.0;
    for (i = 0; i < NOC; ++i) {
        n = (int)get_data(k,i);
        if (n < 1) {
            printf1("Error: level variables must be positive integers.\n");
            goto GDFCFin;
        }
        AcM[ii++] = n;
        id = get_data(j,VSORTPtr[i]);
        if (id != id0) {
            GDFNU++;
            id0 = id;   
            l = 0;
        }
        l++;
    }
    if (sorti(ii,AcM,0))        /* sort and calculate GDFNDim */
        goto GDFCFin;

    j = AcM[0];
    GDFNDim = 1;
    for (i = 1; i < ii; ++i) {
        l = AcM[i];
        if (l != j) {
            GDFNDim++;
            j = l;
        }
    }
    if (alloc_gdfgen(GDFNU,GDFNDim))
        goto GDFCFin;

    GDFDim[0] = AcM[0];                 /* create GDFDim[] */
    k = 0;
    for (i = 1; i < ii; ++i) {
        l = AcM[i];
        if (l != GDFDim[k])  
            GDFDim[++k] = l;
    }
    if (alloc_acm(GDFNDim + 1))
        goto GDFCFin;

    GDFNTyp = 0;
    iid = PM1VIdx[0];
    il1 = PM1VIdx[1];
    i0 = 0;
    k = l = i = 0;
    id0 = get_data(iid,VSORTPtr[0]);
    cen2 = cen3 = 0;
    while (i <= NOC) {
        if (i < NOC) {
            ii = VSORTPtr[i];
            id = get_data(iid,ii);
            l1 = (int)get_data(il1,ii);
            for (t = 0; t < GDFNDim; ++t) {
                if (l1 == GDFDim[t]) {
                    GDFDCnt[t] += 1;
                    if (GDFYTyp[ii] == 2)
                        GDFDInt[t] += 1;
                    else if (GDFYTyp[ii] == 3)
                        GDFDCen[t] += 1;
                    break;
                }
            }
        }
        if (id != id0 || i == NOC) {

            tt = -1;
            for (t = 0; t < GDFNTyp; ++t) {
                if (GDFMLen[t] != l)
                    continue;
                fnd = 1;
                for (j = 0; j < l; ++j) {
                    if (GDFMarg[t * GDFNDim + j] != AcM[j + 1]) {
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
                tt = GDFNTyp;
                for (j = 0; j < l; ++j) {
                    GDFMarg[tt * GDFNDim + j] = AcM[j + 1];
                    if (j > 0 && AcM[j] == AcM[j + 1]) {
                        printf1("Error in unit %g. Observations in the same unit\n",id0);
                        printf1("must belong to different dimensions.\n");
                        goto GDFCFin;
                    }
                }
                GDFMLen[tt] = l;
                GDFNTyp++;
            }
            GDFPtr[k] = i0;
            GDFULen[k] = l;
            GDFMUCnt[tt] += 1;
            GDFMTyp[tt] = imax(GDFMTyp[tt],cen2 + 2 * cen3 + 1);
            k++;
            i0 = i;
            cen2 = cen3 = l = 0;
            id0 = id;
        }
        if (i == NOC)
            break;

        if (GDFYTyp[ii] == 2)
            cen2 = 1;
        else if (GDFYTyp[ii] == 3)
            cen3 = 1;
        l++;
        AcM[l] = l1;
        i++;
    }
    GDFGRP = 1;
    err = 0;

GDFCON:
    printf1("Number of cases: %d\n",NOC);
    printf1("Number of units: %d\n",GDFNU);
    printf1("Number of dimensions (level 1): %d\n",GDFNDim);
    printf1("Number of exact cases: %d\n",GDFNEX);
    printf1("Number of interval censored cases: %d\n",GDFNINT);
    printf1("Number of right censored cases: %d\n",GDFNCEN);

    if (GDFGRP) {
        printf1("Number of marginal patterns: %d\n\n",GDFNTyp);
        gdf_pmtyp();

        printf1("Dimension  level-1  observations  i-censored  r-censored\n");
        for (i = 0; i < GDFNDim; ++i)  
            printf1("%6d     %5d    %12d  %10d  %10d\n",i + 1,
                              GDFDim[i],GDFDCnt[i],GDFDInt[i],GDFDCen[i]);
    }
    newline();

GDFCFin:
    alloc_acm(0);
    alloc_acn(0);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdf_pmtyp()     Print info about marginal patterns.                     */

void gdf_pmtyp(void)
{
    register int j,t;
    int l,d;
    
    printf1("Pattern  type     units  dimensions\n");

    for (t = 0; t < GDFNTyp; ++t) {
        printf1("%5d  %5d %10d  ",t + 1,GDFMTyp[t],GDFMUCnt[t]);
        l = GDFMLen[t];
        for (j = 0; j < l; ++j) {
            d = GDFMarg[t * GDFNDim + j];
            printf1("%2d ",d);
        }
        newline();
    } 
    newline();
}   

/* ------------------------------------------------------------------------ */
/*  gdf_marg(opt)   create marginal distr. functions, write to output file. */
/*                  opt = 0 : distribution function                         */
/*                  opt = 1 : survivor function                             */
/*                  opt = 2 : expected values                               */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int gdf_marg(int opt)
{
    register int i,n;
    int err,d,dl,m,nn,nc,il1;

    err = -1;
    gdf_pmet(1,opt);

    for (d = 0; d < GDFNDim; ++d) {

        dl = GDFDim[d];
        if (GDFGRP) {
            il1 = PM1VIdx[1];           /* index of L1 variable */
            nn = GDFDCnt[d];
            nc = GDFDCen[d];
        }
        else {
            nn = NOC;
            nc = GDFNCEN;
        }
        if (alloc_acy(nn))
            goto GDFMARGFin;

        if (nc > 0) {                   /* AcNS used for censoring indicator */
            if (alloc_acns(nn))
                goto GDFMARGFin;
            if (alloc_acu(nn))          /* used for new y values */
                goto GDFMARGFin;
        }
        if (alloc_acw(nn))
            goto GDFMARGFin;

        n = 0;
        for (i = 0; i < NOC; ++i) {
            if (GDFGRP == 0 || (int)get_data(il1,i) == dl) {
                AcY[n] = get_data(PMYL,i);
                if (nc > 0) {
                    if (GDFYTyp[i] == 3)
                        AcNS[n] = 0;
                    else
                        AcNS[n] = 1;
                }
                n++;
            }
        }
        if (nc == 0) {                      /* only exact observations */

            printf1("Dimension %3d: EDF with exact data.\n",dl);

            m = gdf_edf1(n,AcY,AcW,opt);    /* get edf */         
            if (m < 1)
                goto GDFMARGFin;
            if (opt == 2)
                gdf_pdat1a(dl,m,AcY,AcW,AcNS,1);     /* print to output file */
            else
                gdf_pdat1(dl,m,AcY,AcW);        
        }
        else {                              /* exact and r-censored obs. */ 

            printf1("Dimension %3d: Kaplan-Meier.\n",dl);

            m = gdf_edf2(n,AcY,AcNS,AcU,AcW,opt);   /* get edf */         
            if (m < 1)
                goto GDFMARGFin;
            if (opt == 2)
                gdf_pdat1a(dl,m,AcU,AcW,AcNS,0);     /* print to output file */
            else
                gdf_pdat1(dl,m,AcU,AcW);        
        }
    }
    err = 0;

GDFMARGFin:
    alloc_acns(0);
    alloc_acu(0);
    alloc_acw(0);
    alloc_acy(0);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdf_pdat1(d,n,y,f)    Print distribution function to output file.       */
   
void gdf_pdat1(int d,int n,double *y,double *f)
{
    register int i;

    for (i = 0; i < n; ++i) {
        fprintf(PMFd,"%4d ",d);
        fprintf(PMFd,PMFmtS,y[i]);
        fprintf(PMFd,PMFmtS,f[i]);
        fprintf(PMFd,"\n");
        GDFWNOC++;
    }
}

/* ------------------------------------------------------------------------ */
/*  gdf_pdat1a(d,n,y,f,cen,opt)   Print distribution function to output     */
/*                                file. If opt != 0 do not use cen but      */
/*                                always print 1.                           */
   
void gdf_pdat1a(int d,int n,double *y,double *f,short *cen,int opt)
{
    register int i,m;

    m = 1;
    for (i = 0; i < n; ++i) {
        fprintf(PMFd,"%4d ",d);
        fprintf(PMFd,PMFmtS,y[i]);
        if (opt == 0)
            m = cen[i];                
        fprintf(PMFd,"%2d ",m);
        fprintf(PMFd,PMFmtS,f[i]);
        fprintf(PMFd,"\n");
        GDFWNOC++;
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
   
int gdf_edf1(int n,double *y,double *f,int opt)
{
    register int i,k,m;
    double y0,y1,f0;

    if (opt == 2) {
        for (i = 0; i < n; ++i)  
            f[i] = y[i];
        return(n);
    }
    if (alloc_actmp(n))
        return(-1);

    for (i = 0; i < n; ++i)
        AcTmp[i] = y[i];

    if (sortd(n,AcTmp,0))
        return(-1);    

    y0 = AcTmp[0];
    m = 0;
    k = 1;
    f0 = 0.0;
    for (i = 1; i <= n; ++i) {
        if (i < n)
            y1 = AcTmp[i];

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
    alloc_actmp(0);
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

int gdf_edf2(int n,double *y,short *cen,double *y1,double *h,int opt)
{
    register int i,j,k,l;
    int m;
    double tmp,tmp1;

    if (alloc_ack(n))
        return(-1);
    if (alloc_actmp(n))
        return(-1);

    if (sortdp2a(n,y,cen,AcK,0))      /* sort */
        return(-1);    

    tmp = 1.0 / (double)n;
    m = n - 1;
    for (i = 0; i < n; ++i) {
        j = AcK[i];
        AcTmp[j] += tmp;
        if (cen[j] == 0 && i < n - 1) {
            if (m > 0) {
                tmp1 = AcTmp[j] / (double)m;
                for (k = i + 1; k < n; ++k)  
                    AcTmp[AcK[k]] += tmp1;
            }
            AcTmp[j] = 0.0;
        }
        m--;
    }
    if (opt == 2) {         /* expected values */
        for (i = 0; i < n; ++i) {
            j = AcK[i];
            y1[j] = y[j];
            if (cen[j] || i == n - 1)  
                h[j] = y[j];
            else {
                tmp1 = tmp = 0.0;
                for (k = i + 1; k < n; ++k) {
                    l = AcK[k];
                    tmp1 += y[l] * AcTmp[l];
                    tmp += AcTmp[l];
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
        j = AcK[0];
        y1[m] = y[j];
        h[m] = AcTmp[j];

        for (i = 1; i < n; ++i) {
            j = AcK[i];
            if (AcTmp[j] < EPSI1)
                continue;

            if (y[j] == y1[m])
                h[m] += AcTmp[j];
            else {
                m++;
                y1[m] = y[j];
                h[m] = h[m - 1] + AcTmp[j];
            }
        }
        m++;
        if (opt == 1) {                 /* surv function */
            for (i = 0; i < m; ++i)
                h[i] = 1.0 - h[i];
        }
    }
    alloc_ack(0);
    alloc_actmp(0);
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

int gdf_joint(int met,int opt)
{
    int err,t,m,nu,ndim;

    err = -1;
    
    gdf_pmet(met,opt);
    if (GDFGRP == 0) {
        printf1("Error: grp parameter required.\n");
        return(-1);
    }

    for (t = 0; t < GDFNTyp; ++t) {     /* for all patterns */

        nu = GDFMUCnt[t];
        ndim = GDFMLen[t];

        if (PMProtFDef) {
            fprintf(PMProtFd,"Current marginal pattern: %d  Typ: %d\n",t+1,GDFMTyp[t]);
            fprintf(PMProtFd,"Number of units: %d  dimensions: %d\n",nu,ndim);
        }   
        if (alloc_aci(nu * ndim + 1))   /* pointer for cases */
            goto GDFJFin;

        gdf_cptr(t,nu,ndim,AcI);        /* create pointer to data matrix cases */

        if (GDFMTyp[t] == 1) {          /* only exact observations */

            printf1("\nPattern %d: EDF with exact data.\n",t + 1);
            if (gdf_domain(t,nu,ndim,AcI,NULL,1,0.0))
                goto GDFJFin;

            if (alloc_acx(nu * ndim + 1))   /* values */
                goto GDFJFin;
            if (alloc_acw(nu + 1))          /* distr function */
                goto GDFJFin;

            m = gdf_edf1m(t,nu,ndim,AcI,AcX,AcW,opt);
            if (m < 1)
                goto GDFJFin;

            if (opt == 2)
                gdf_pdat2a(t,m,ndim,AcX,AcI,1);  
            else
                gdf_pdat2(t,m,ndim,AcX,AcW);    

        }
        else if (GDFMTyp[t] == 3) {     /* exact and r-censored observations */

            printf1("\nPattern %d. Units: %d, dimensions: %d\n",t + 1,nu,ndim);

            if (alloc_acv(nu * ndim + 1))           /* values/lower bounds */
                goto GDFJFin;
            if (alloc_acw(nu + 1))                  /* distr function */
                goto GDFJFin;

            if (met == 2) {
                printf1("Iterative procedure.\n");
                printf1("Max iterations: %d  tolerance: %g\n",MxIter,TOLF);

                if (gdf_domain(t,nu,ndim,AcI,NULL,1,0.0))
                    goto GDFJFin;

                gdf_getb(PMN,ndim,1);                   /* get number of boxes */
    
                if (gdf_edf3m(t,nu,ndim,AcI,NULL,AcV,AcW,opt,1,MxIter,TOLF))
                    goto GDFJFin;
            }
            else {
                printf1("Local Kaplan-Meier.\n");
                if (gdf_domain(t,nu,ndim,AcI,NULL,1,PMD))
                    goto GDFJFin;

                if (gdf_edf4m(t,nu,ndim,AcI,NULL,AcV,AcW,PMD,opt,1))
                    goto GDFJFin;
            }
            if (opt == 2)
                gdf_pdat2a(t,nu,ndim,AcV,AcI,0);  
            else
                gdf_pdat2(t,nu,ndim,AcV,AcW);    
        }
    }
    err = 0;

GDFJFin:
    alloc_aci(0);
    alloc_acx(0);
    alloc_acw(0);
    alloc_acv(0);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gdf_pdat2(t,m,ndim,val,f)                                               */
/*                                                                          */
/*  write joint distr function to output file. t is number of pattern.      */
/*  val is (m,ndim), f is m vector with masses.                             */
   
void gdf_pdat2(int t,int m,int ndim,double *val,double *f)
{
    register int i,j;

    for (i = 0; i < m; ++i) {
        fprintf(PMFd,"%3d  ",t + 1);
        for (j = 0; j < GDFNDim; ++j)
            fprintf(PMFd,"%2d ",GDFMarg[t * GDFNDim + j]);

        for (j = 0; j < ndim; ++j)
            fprintf(PMFd,PMFmtS,val[i * ndim + j]);

        fprintf(PMFd,PMFmtS,f[i]);
        fprintf(PMFd,"\n");
        GDFWNOC++;
    }
}

/* ------------------------------------------------------------------------ */
/*  gdf_pdat2a(t,m,ndim,val,ptr,opt)                                        */
/*                                                                          */
/*  write expected values. If opt = 1 val contains the exact values and     */
/*  ptr is not used. If opt = 0, val contains the expected values and       */
/*  ptr provides pointer to original observations. In this case, m = nu.    */
   
void gdf_pdat2a(int t,int m,int ndim,double *val,int *ptr,int opt)
{
    register int i,j,ii;
    int n;
    double tmp;

    n = 1;
    for (i = 0; i < m; ++i) {
        ii = i * ndim;
        fprintf(PMFd,"%3d  ",t + 1);
        for (j = 0; j < GDFNDim; ++j)
            fprintf(PMFd,"%2d ",GDFMarg[t * GDFNDim + j]);

        if (opt) {
            for (j = 0; j < ndim; ++j)
                fprintf(PMFd,PMFmtS,val[ii + j]);
            for (j = 0; j < ndim; ++j)
                fprintf(PMFd,"%d ",n);
            for (j = 0; j < ndim; ++j)
                fprintf(PMFd,PMFmtS,val[ii + j]);
        }
        else {
            for (j = 0; j < ndim; ++j) {
                tmp = get_data(PMYL,ptr[ii + j]);
                fprintf(PMFd,PMFmtS,tmp);
            }
            for (j = 0; j < ndim; ++j) {
                n = 1;
                if (GDFYTyp[ptr[ii + j]] == 3)
                    n = 0;
                fprintf(PMFd,"%d ",n);
            }
            for (j = 0; j < ndim; ++j)
                fprintf(PMFd,PMFmtS,val[ii + j]);
        }
        fprintf(PMFd,"\n");
        GDFWNOC++;
    }
}

/* ------------------------------------------------------------------------ */
/*  gdf_cptr(t,nu,ndim,ptr)                                                 */
/*                                                                          */
/*  Create pointer ptr[i,j] (i=0,...,nu-1, j=0,...,ndim-1)  such that       */
/*  ptr[i,j] is pointer to data matrix row for the corresponding case.      */
/*                                                                          */
   
void gdf_cptr(int t,int nu,int ndim,int *ptr)
{
    register int i,j,k,l;
    int l1,il1,len,t1,d1,iu;
    double y;

    il1 = PM1VIdx[1];
    iu = 0;
    for (i = 0; i < GDFNU; ++i) {
        if (iu >= nu)
            break;

        t1 = t * GDFNDim;
        d1 = GDFMarg[t1]; 
        len = GDFULen[i];

        if (len != ndim)     
            continue;

        k = GDFPtr[i];   
        for (l = 0; l < len; ++l) {
            d1 = GDFMarg[t1++]; 
            j = VSORTPtr[k + l];

            y = get_data(PMYL,j);
            l1 = (int)get_data(il1,j);
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
        printfe("ERROR in gdf_cptr (iu=%d, nu=%d)\n",iu,nu);
        gerr_exit(213);
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
   
int gdf_domain(int t,int nu,int ndim,int *ptr,double *yval,int opt,double d)
{
    register int i,j,k;
    int err;
    double y,dl,dh;

    err = 0;
    for (j = 0; j < ndim; ++j) {
        k = ptr[j];
        if (yval == NULL)
            dl = dh = get_data(PMYL,k);
        else                 
            dl = dh = yval[k];

        if (PMYH >= 0 && GDFYTyp[k] < 3) {
            if (yval == NULL)
                dh = get_data(PMYH,k);
        }

        for (i = 1; i < nu; ++i) {
            k = ptr[i * ndim + j];
            if (yval == NULL)
                y = get_data(PMYL,k);
            else
                y = yval[k];

            dl = dmin(dl,y);
            if (PMYH >= 0 && GDFYTyp[k] < 3) {
                if (yval == NULL)
                    y = get_data(PMYH,k);
            }
            dh = dmax(dh,y);
        }
        if (PMSC > 0.0) {
            dl -= PMSC;
            dh += PMSC;
        }
        GDFDL[j] = dl;
        GDFDH[j] = dh;
        if (dh <= dl + EPSI1)
            err = -1;
    }
    if (opt == 1) {  
        printf1("Domain");
        if (PMSC > 0.0)
            printf1(" (sc=%g)",PMSC);
        newline();
        for (j = 0; j < ndim; ++j) {
            printf1("%3d ",GDFMarg[t * GDFNDim + j]);
            printf1(PMTFmtS,GDFDL[j]);
            printf1(PMTFmtS,GDFDH[j]);
            if (d > 0.0)
                printf1("  delta: %g",d * (GDFDH[j] - GDFDL[j]));
            newline();
        }
        newline();
    }
    else if (opt == 2 && PMProtFDef) {
        fprintf(PMProtFd,"Domain");
        if (PMSC > 0.0)
            fprintf(PMProtFd," (sc=%g)",PMSC);
        fprintf(PMProtFd,"\n");
        for (j = 0; j < ndim; ++j) {
            fprintf(PMProtFd,"%3d ",GDFMarg[t * GDFNDim + j]);
            fprintf(PMProtFd,PMTFmtS,GDFDL[j]);
            fprintf(PMProtFd,PMTFmtS,GDFDH[j]);
            if (d > 0.0)
                fprintf(PMProtFd,"  delta: %g",d * (GDFDH[j] - GDFDL[j]));
            fprintf(PMProtFd,"\n");
        }
        fprintf(PMProtFd,"\n");
    }
    if (err)
        printf1("Error: domain is empty.\n");
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
   
int gdf_edf1m(int t,int nu,int ndim,int *ptr,double *val,double *f,int opt)
{
    register int i,j,k,l,ii,kk;
    int sflg,m,mm;
    double tmp;

    if (alloc_actmp(nu * ndim + 1))
        return(-1);
    if (alloc_acj(ndim + 1))
        return(-1);
    if (alloc_ack(nu + 1))
        return(-1);

    for (j = 1; j <= ndim; ++j)
        AcJ[j] = j;

    for (i = 0; i < nu; ++i) {
        l = i * ndim;
        for (j = 0; j < ndim; ++j) {
            k = ptr[l];
            l++;
            AcTmp[l] = get_data(PMYL,k);
        }
    }
    if (sortdpn(nu,ndim,AcTmp,ndim,AcJ,AcK))
        return(-1);

    if (opt == 2) {
        for (i = 0; i < nu; ++i) {
            k = (AcK[i + 1] - 1) * ndim + 1;
            ii = i * ndim;
            for (j = 0; j < ndim; ++j) 
                val[ii + j] = AcTmp[k + j];
        }
        m = nu;
        goto EDF1MFin;
    }

    tmp = 1.0 / (double)nu;

    if (opt == 0) {                     /* distribution function */
        m = 0;
        for (i = 0; i < nu; ++i) {
            k = (AcK[i + 1] - 1) * ndim + 1;
            if (i == 0) {
                for (j = 0; j < ndim; ++j) 
                    val[j] = AcTmp[k + j];
                f[0] = tmp;
            }
            else {
                sflg = 1;
                mm = m * ndim;
                for (j = 0; j < ndim; ++j) {
                    if (val[mm + j] != AcTmp[k + j]) {
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
                        val[mm + j] = AcTmp[k + j];
                    f[m] = tmp;
                    for (ii = 0; ii < i; ++ii) {
                        kk = (AcK[ii + 1] - 1) * ndim + 1;
                        sflg = 1;
                        for (j = 0; j < ndim; ++j) {
                            if (val[mm + j] < AcTmp[kk + j]) {
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
            k = (AcK[i + 1] - 1) * ndim + 1;
            if (i == nu - 1) {
                for (j = 0; j < ndim; ++j) 
                    val[j] = AcTmp[k + j];
                f[0] = 0.0;
            }
            else {
                sflg = 1;
                mm = m * ndim;
                for (j = 0; j < ndim; ++j) {
                    if (val[mm + j] != AcTmp[k + j]) {
                        sflg = 0;
                        break;
                    }
                }
                if (sflg == 0) {
                    m++;
                    mm = m * ndim;
                    for (j = 0; j < ndim; ++j) 
                        val[mm + j] = AcTmp[k + j];
                    f[m] = 0.0;
                    for (ii = i + 1; ii < nu; ++ii) {
                        kk = (AcK[ii + 1] - 1) * ndim + 1;
                        sflg = 1;
                        for (j = 0; j < ndim; ++j) {
                            if (val[mm + j] >= AcTmp[kk + j]) {
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
    alloc_actmp(0);
    alloc_acj(0);
    alloc_ack(0);
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
   
int gdf_edf3m(int t,int nu,int ndim,int *ptr,double *yval,double *val,
    double *f,int opt,int prn,int mxit,double tol)
{
    register int i,j,k,l;
    int fin,cen,n,np,iter,nred,ne,kp,nep;
    double y,tmp0,tmp,tmax,tol1;

    tmp0 = 1.0 / (double)nu;

    if (alloc_ack(nu))              /* number of boxes for unit */
        return(-1);
    if (alloc_acu(ndim))            /* temporary for y values */
        return(-1);
    if (alloc_acm(ndim))            /* for different use */
        return(-1);
    if (alloc_acj(ndim))            /* for different use */
        return(-1);
    if (alloc_acns(nu))             /* count censored dimensions */
        return(-1);
    if (alloc_acn(nu * ndim + 1))   /* pointer to grid */
        return(-1);
    if (alloc_acc(nu * ndim + 1))   /* set with status of observation */
        return(-1);

    /* determine pointer to the location of points in grid for
       all units, save in AcN[]. Set AcC to status of observation
       (1 - exact, 2 = icensored, 3 = r-censored).
       Set AcNS = -1 if unit falls outside of domain.       
                   0 if exact in all dimensions
                   1 otherwise     
    */

    if (prn)            /* print table */
        printf1("Dimension  i-censored  r-censored  observations\n");

    for (i = 0; i < nu; ++i) {
        l = i * ndim;
        cen = 0;
        for (j = 0; j < ndim; ++j) {
            k = ptr[l + j];
            if (yval == NULL)
                y = get_data(PMYL,k);
            else
                y = yval[k];

            AcC[l + j] = (char)GDFYTyp[k];
            if (GDFYTyp[k] == 2) { 
                AcM[j] += 1;
                cen = 1;
            }   
            else if (GDFYTyp[k] == 3) {
                AcJ[j] += 1;
                cen = 1;
                y += EPSI1;             /* add small offset */
            }
            AcU[j] = y;
        }
        AcNS[i] = cen;

        for (j = ndim - 1; j >= 0; --j) {
            n = (int)ceil((AcU[j] - GDFDL[j]) / GDFBLen[j]);
            if (n < 1 || n > GDFNBOX) {
                AcNS[i] = -1;
                break;
            }
            AcN[l + j] = n;
        }
    }
    if (prn) {      /* print table with i/r-censored observations */
        for (j = 0; j < ndim; ++j)  
            printf1("%6d    %11d %11d  %12d\n",
                             GDFMarg[t * GDFNDim + j],AcM[j],AcJ[j],nu);
        newline();
    }

    /*  create pointers for grid in AcJ */
    /*  first set AcJ > 0 for all boxes in grid that have an observation */
    /*  AcJ = 1 for exact observations, AcJ = 2 otherwise */
    /*  Also set AcK = number of boxes for unit */
    /*  count number of boxes that have exact observations in ne */

    ne = 0;

    if (alloc_acj(GDFTBOX + 1))  
        return(-1);

    if (PMProtFDef)  
        fprintf(PMProtFd,"\nLocation of observations in grid.\n");

    /***
    printf1("ACN\n");
    for (i=0;i< nu; ++i) {
        for (j=0;j<ndim; ++j)
            printf1("%d ",AcN[i * ndim + j]);
        newline();
    }
    ****/

    for (i = 0; i < nu; ++i) {      /* for all units */

        if (AcNS[i] < 0)            /* skip observations outside of grid */
            continue;   

        l = i * ndim;
        for (j = 0; j < ndim; ++j)   
            AcM[j] = AcN[l + j];

        np = gdf_boxptr(ndim,AcM);  /* pointer to first box of observation */

        if (PMProtFDef)  
            fprintf(PMProtFd,"%6d %3d %5d\n",i,AcNS[i],np);

        AcK[i] = 1;                 /* count boxes */
        if (AcNS[i] == 0) {         /* if exact */    
            AcJ[np] = 1;
            continue;
        }

        fin = 0;
        while (fin == 0) {
        
            if (AcJ[np] == 0)   
                AcJ[np] = 2;
        
        
            fin = 1;
            for (j = 0; j < ndim; ++j) {
                if (AcC[l + j] == 3 && AcM[j] < GDFNBOX) {
                    AcM[j] += 1;
                    for (k = 0; k < j; ++k)     /* reset */
                        AcM[k] = AcN[l + k];
                    fin = 0;
                    break;
                }
            }
            if (fin == 0) {
                np = gdf_boxptr(ndim,AcM);
                AcK[i] += 1;
                if (PMProtFDef)  
                    fprintf(PMProtFd,"%6d %3d %5d\n",i,AcNS[i],np);
            }
        }
    }
    if (PMProtFDef)  
        fprintf(PMProtFd,"\nPointers to reduced grid\n");

    /* AcJ[np] is index of np in reduced grid */
    /* nred is the number of locations in the reduced grid */

    


    nred = 0;
    for (i = 0; i < GDFTBOX; ++i) {
        if (AcJ[i] > 0) {
            if (AcJ[i] == 1)
                ne++;
            AcJ[i] = nred++;
            if (PMProtFDef)  
                fprintf(PMProtFd,"%4d %4d\n",i,AcJ[i]);
        }
        else
            AcJ[i] = -1;
    }
    if (PMProtFDef) {
        fprintf(PMProtFd,"Number of boxes, total: %d  reduced: %d\n",GDFTBOX,nred);
        fprintf(PMProtFd,"Number of boxes with exact observations: %d\n",ne);
    }

    /* create two arrays for reduced grid.
       AcXF for densities, AcYF for redistribution weights */

    if (alloc_acxf(nred))          
        return(-1);
    if (alloc_acyf(nred))          
        return(-1);

    /*  create initial weights in AcXF */

    for (i = 0; i < nu; ++i) {      /* for all units */

        if (AcNS[i] < 0)            /* skip observations outside of grid */
            continue;   

        l = i * ndim;
        for (j = 0; j < ndim; ++j) 
            AcM[j] = AcN[l + j];

        np = gdf_boxptr(ndim,AcM);  /* pointer to first box of observation */

        if (AcNS[i] == 0) {         /* if exact then unit mass */
            AcXF[AcJ[np]] += tmp0;  
            continue;
        }
        fin = 0;
        while (fin == 0) {

            AcXF[AcJ[np]] += tmp0 / (double)AcK[i];

            fin = 1;
            for (j = 0; j < ndim; ++j) {
                if (AcC[l + j] == 3 && AcM[j] < GDFNBOX) {
                    AcM[j] += 1;
                    for (k = 0; k < j; ++k)     /* reset */
                        AcM[k] = AcN[l + k];
                    fin = 0;
                    break;
                }
            }
            if (fin == 0)  
                np = gdf_boxptr(ndim,AcM);
        }
    }
    if (PMProtFDef) {
        fprintf(PMProtFd,"\nInitial distribution\n");
        tmp = 0.0;
        for (j = 0; j < nred; ++j)  {
            tmp+=AcXF[j];
            fprintf(PMProtFd,"%3d  %8.6f\n",j,AcXF[j]);
        }
        fprintf(PMProtFd,"Sum of values: %g\n",tmp);
    }

    /* perform iterations */

    tmax = 0.0;
    for (iter = 1; iter <= mxit; ++iter) {

        /* update distribution. Use weights in AcXF, create new
           densities in AcYF */

        for (i = 0; i < nu; ++i) {      /* for all units */

            if (AcNS[i] < 0)            /* skip observations outside of grid */
                continue;   

            l = i * ndim;
            for (j = 0; j < ndim; ++j) 
                AcM[j] = AcN[l + j];

            np = gdf_boxptr(ndim,AcM);  /* pointer to first box of observation */

            if (AcNS[i] == 0) {         /* if exact then unit mass */
                AcYF[AcJ[np]] += tmp0;  
                continue;
            }
            tmp = 0.0;                  /* first find sum of weights */
            fin = 0;
            while (fin == 0) {

                tmp += AcXF[AcJ[np]];
                fin = 1;
                for (j = 0; j < ndim; ++j) {
                    if (AcC[l + j] == 3 && AcM[j] < GDFNBOX) {
                        AcM[j] += 1;
                        for (k = 0; k < j; ++k)     /* reset */
                            AcM[k] = AcN[l + k];
                        fin = 0;
                        break;
                    }
                }
                if (fin == 0)  
                    np = gdf_boxptr(ndim,AcM);
            }

            /* now add new weights */

            if (tmp < EPSI1)
                continue;
            tmp = tmp0 / tmp;

            for (j = 0; j < ndim; ++j) 
                AcM[j] = AcN[l + j];

            np = gdf_boxptr(ndim,AcM);  /* pointer to first box of observation */

            fin = 0;                /* add new weights */
            while (fin == 0) {

                AcYF[AcJ[np]] += AcXF[AcJ[np]] * tmp;

                fin = 1;
                for (j = 0; j < ndim; ++j) {
                    if (AcC[l + j] == 3 && AcM[j] < GDFNBOX) {
                        AcM[j] += 1;
                        for (k = 0; k < j; ++k)     /* reset */
                            AcM[k] = AcN[l + k];
                        fin = 0;
                        break;
                    }
                }
                if (fin == 0)  
                    np = gdf_boxptr(ndim,AcM);
            }
        }
        tmax = 0.0;
        for (j = 0; j < nred; ++j) {
            if (AcYF[j] < EPSI)
                AcYF[j] = 0.0;
            tmax = dmax(tmax,fabs(AcXF[j] - AcYF[j]));
            AcXF[j] = AcYF[j];
            AcYF[j] = 0.0;
        }
        printfe("Iter %2d  Crit %15.8e\n",iter,tmax);
        if (iter >= mxit || tmax <= tol)
            break;
    }
    printf1("Finished after %d iterations. Max change: %g\n",iter,tmax);

    if (PMProtFDef) {
        fprintf(PMProtFd,"\nFinal densities after %d iterations.\n",iter);
        tmp = 0.0;
        for (j = 0; j < nred; ++j)  {
            fprintf(PMProtFd,"%3d  %8.6f\n",j,AcXF[j]);
            tmp += AcXF[j];
        }
        fprintf(PMProtFd,"Sum: %g  Maximal change: %g\n\n",tmp,tmax);
    }

    if (opt <= 1) {

        /*  calculate distribution or survivor function, depending on opt */
        /*  values in val[nu,ndim], function values in f[] */

        if (alloc_acr(ndim))    
            return(-1);
        if (alloc_acs(ndim))    
            return(-1);

        for (i = 0; i < nu; ++i) {

            l = i * ndim;
            for (j = 0; j < ndim; ++j) {
                k = ptr[l + j];
                if (yval == NULL)
                    val[l + j] = get_data(PMYL,k);
                else
                    val[l + j] = yval[k];

                if (opt == 0) {         /* distribution function */
                    AcR[j] = AcM[j] = 1;
                    AcS[j] = AcN[l + j];
                }
                else {                  /* survivor function */
                    AcS[j] = GDFNBOX;
                    AcM[j] = AcN[l + j];
                    /*******************
                    if (AcM[j] < GDFNBOX)
                        AcM[j] += 1;
                    ********************/
                    AcR[j] = AcM[j];
                }
            }
            if (AcNS[i] < 0) {
                f[i] = -1.0;
                continue;
            }
            np = gdf_boxptr(ndim,AcM);  /* pointer to first box of observation */

            tmp = 0.0;
            fin = 0;    
            while (fin == 0) {
                if (np >= 0 && np < GDFTBOX) {
                    k = AcJ[np];
                    if (k >= 0 && k < nred)
                        tmp += AcXF[k];
                }
                fin = 1;
                for (j = 0; j < ndim; ++j) {
                    if (AcM[j] < AcS[j]) {
                        AcM[j] += 1;
                        for (k = 0; k < j; ++k)     /* reset */
                            AcM[k] = AcR[k];
                        fin = 0;
                        break;
                    }
                }
                if (fin == 0)  
                    np = gdf_boxptr(ndim,AcM);
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

            if (alloc_acyf(ne * ndim))          
                return(-1);
            if (alloc_acs(ne))          
                return(-1);
            if (alloc_acr(nred))          
                return(-1);

            for (i = 0; i < nred; ++i)
                AcR[i] = -1;

            nep = 0;
            for (i = 0; i < nu; ++i) {      /* for all units */

                if (AcNS[i] != 0)           /* only exact observations inside of grid */
                    continue;   

                l = i * ndim;
                for (j = 0; j < ndim; ++j)   
                    AcM[j] = AcN[l + j];

                np = gdf_boxptr(ndim,AcM);  /* pointer to first box of observation */

                k = AcJ[np];
                if (AcR[k] < 0)  
                    AcR[k] = nep++;
                kp = AcR[k];
                for (j = 0; j < ndim; ++j) {
                    if (yval == NULL)
                        AcYF[kp * ndim + j] += get_data(PMYL,ptr[l + j]);
                    else
                        AcYF[kp * ndim + j] += yval[ptr[l + j]];
                }
                AcS[kp] += 1;
            }
            if (ne != nep) {          
                printfe("ERROR in gdf_edf3m: ne=%d nep=%d\n",ne,nep);
                gerr_exit(214);
            }
            for (i = 0; i < ne; ++i) {
                if (AcS[i] > 1) {
                    k = i * ndim;
                    for (j = 0; j < ndim; ++j) 
                        AcYF[k++] /= (double)AcS[i];
                }
            }

            if (PMProtFDef) {
                fprintf(PMProtFd,"Pointer for mapping exact observations\n");
                for (i = 0; i < nred; ++i)  
                    fprintf(PMProtFd,"%6d %6d\n",i,AcR[i]);
                fprintf(PMProtFd,"Array with mean exact observations\n");
                for (i = 0; i < ne; ++i) {
                    fprintf(PMProtFd,"%6d %6d ",i,AcS[i]);
                    for (j = 0; j < ndim; ++j)
                        fprintf(PMProtFd,PMTFmtS,AcYF[i * ndim + j]);
                    fprintf(PMProtFd,"\n");
                }
            }
        }
        if (PMProtFDef)
            fprintf(PMProtFd,"\nUnit  status  values  expectation\n");

        if (alloc_acy1(ndim))    
            return(-1);
        if (alloc_actmp(ndim))    
            return(-1);

        tol1 = tol / 10.0;

        for (i = 0; i < nu; ++i) {

            l = i * ndim;

            if (PMProtFDef) {
                fprintf(PMProtFd,"%4d  ",i);
                for (j = 0; j < ndim; ++j)
                    fprintf(PMProtFd,"%2d ",(int)AcC[l + j]);
            }
            for (j = 0; j < ndim; ++j) {
                k = ptr[l + j];
                if (yval == NULL)
                    val[l + j] = get_data(PMYL,k);
                else
                    val[l + j] = yval[k];

                if (PMProtFDef)   
                    fprintf(PMProtFd,PMFmtS,val[l + j]);
            }
            if (AcNS[i] <= 0) {
                if (PMProtFDef)   
                    fprintf(PMProtFd,"\n");
                continue;
            }
            for (j = 0; j < ndim; ++j) {               
                AcM[j] = AcN[l + j];
                /*************************************** 
                if (AcC[l + j] == 3 && AcM[j] < GDFNBOX)  
                    AcM[j] += 1;
                **************************************/
                AcTmp[j] = AcY1[j] = 0.0;
            }
            np = gdf_boxptr(ndim,AcM);  /* pointer to first box of observation */

            fin = 0;
            while (fin == 0) {
                if (np >= 0 && np < GDFTBOX) {
                    k = AcJ[np];
                    if (k >= 0 && k < nred) {
                        tmp = AcXF[k];

                        if (tmp > tol1) {
                            gdf_getmval(ne,ndim,np,nred,AcM,AcU);

                            for (j = 0; j < ndim; ++j) {
                                if (AcC[l + j] == 3) {
                                    AcY1[j] += AcU[j] * tmp;
                                    AcTmp[j] += tmp;
                                }   
                            }
                        }
                    }
                }
                fin = 1;
                for (j = 0; j < ndim; ++j) {
                    if (AcC[l + j] == 3 && AcM[j] < GDFNBOX) {
                        AcM[j] += 1;
                        for (k = 0; k < j; ++k) {    /* reset */
                            AcM[k] = AcN[l + k];
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
                    np = gdf_boxptr(ndim,AcM);
            }
            for (j = 0; j < ndim; ++j) {
                if (AcC[l + j] == 3 && AcTmp[j] > 0.0)  
                    val[l + j] = dmax(val[l + j],AcY1[j] / AcTmp[j]);
                if (PMProtFDef)   
                    fprintf(PMProtFd,PMFmtS,val[l + j]);
            }
            if (PMProtFDef)   
                fprintf(PMProtFd,"\n");
        }
        if (PMProtFDef)   
            fprintf(PMProtFd,"\n");
    }
    alloc_actmp(0);
    alloc_acy1(0);
    alloc_acs(0);
    alloc_acr(0);
    alloc_acyf(0);
    alloc_acxf(0);
    alloc_acc(0);
    alloc_acn(0);
    alloc_acns(0);
    alloc_acj(0);
    alloc_acm(0);
    alloc_acu(0);
    alloc_ack(0);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  gdf_getb(n,ndim,opt)   Calculate number of boxes and set in GDFNBOX     */
/*                         Also calculate GDFBLen.                          */
/*  If opt=1 print to stdout, if opt=2 print to protocol file.              */
   
void gdf_getb(int n,int ndim,int opt)
{
    register int j; 

    if (n < 1)  
        n = 100;     

    GDFNBOX = (int)(rexp(rlog((double)n) / (double)ndim) + 0.5);
    if (GDFNBOX < 1)
        GDFNBOX = 1;

    GDFTBOX = (int)pow((double)GDFNBOX,(double)ndim);

    if (opt == 1)
        printf1("Number of boxes in each dimension: %d. Total: %d\n",GDFNBOX,GDFTBOX);
    else if (opt == 2 && PMProtFDef)
        fprintf(PMProtFd,"Number of boxes in each dimension: %d. Total: %d\n",GDFNBOX,GDFTBOX);

    for (j = 0; j < ndim; ++j)    
        GDFBLen[j] = (GDFDH[j] - GDFDL[j]) / (double)GDFNBOX;
}

/* ------------------------------------------------------------------------ */
/*  gdf_boxptr(ndim,ptr)                                                    */
   
int gdf_boxptr(int ndim,int *ptr)
{
    register int j,np;

    np = 0;
    for (j = ndim - 1; j >= 0; --j) {
        np *= GDFNBOX;
        np += (ptr[j] - 1);    
    }
    return(np);
}

/* ------------------------------------------------------------------------ */
/*  gdf_getmval(ne,ndim,np,nred,acm,mval)                                   */
/*                                                                          */
/*  find mean value for box np and return in mval[].                        */
   
void gdf_getmval(int ne,int ndim,int np,int nred,int *acm,double *mval)
{
    register int k,j; 

    if (ne > 0 && np >= 0 && np < GDFTBOX) {

        if ((k = AcJ[np]) >= 0 && k < nred) {
            if ((k = AcR[k]) >= 0) {
                for (j = 0; j < ndim; ++j)       
                    mval[j] = AcYF[k * ndim + j];
                return;
            }
        }
    }
                                  /* use box means */
    for (j = 0; j < ndim; ++j)  
        mval[j] = GDFDL[j] + ((double)acm[j] - 0.5) * GDFBLen[j];
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
   
int gdf_edf4m(int t,int nu,int ndim,int *ptr,double *yval,double *val,
    double *f,double d,int opt,int prn)                      
{
    register int i,j,k,l,ii;
    int err,cen,n,ll,kk,sflg;
    double tmp,tmp0;

    err = -1;
    if (alloc_acns(nu))             /* index of censored dimension */
        goto EDF4MFin;              /* -1 if not censored */
    if (alloc_acj(ndim))            /* number of censored observations */
        goto EDF4MFin;              /* in dimension j */

    if (prn)                        /* print table */
        printf1("Dimension  i-censored  r-censored  observations\n");

    n = 0;
    for (i = 0; i < nu; ++i) {
        l = i * ndim;
        cen = 0;
        AcNS[i] = -1;
        for (j = 0; j < ndim; ++j) {
            k = ptr[l + j];
            if (yval == NULL)
                val[l + j] = get_data(PMYL,k);
            else
                val[l + j] = yval[k];

            if (GDFYTyp[k] == 3) {
                AcNS[i] = j;  
                AcJ[j] += 1;
                cen++;   
            }
        }
        if (cen > 1)
            n++;      
    }
    if (prn) {      /* print table with i/r-censored observations */
        for (j = 0; j < ndim; ++j)  
            printf1("%6d    %11d %11d  %12d\n",
                             GDFMarg[t * GDFNDim + j],0,AcJ[j],nu);
        newline();
    }
    if (n > 0) {
        printf1("Error: found %d observation(s) having more than a single censored dimension.\n",n);
        goto EDF4MFin;
    }

    /* sort observations */

    if (alloc_acj(ndim + 1))
        goto EDF4MFin;
    if (alloc_ack(nu + 1))
        goto EDF4MFin;

    for (j = 1; j <= ndim; ++j)
        AcJ[j] = j;

    if (sortdpn1(nu,ndim,val - 1,ndim,AcJ,AcK - 1,AcNS - 1))
        goto EDF4MFin;

    /* calculate densities in AcXF */

    if (alloc_acxf(nu))
        goto EDF4MFin;
    if (alloc_acu(nu))              /* used for delta / 2 */
        goto EDF4MFin;
    if (alloc_acr(nu))      
        goto EDF4MFin;

    for (j = 0; j < ndim; ++j)   
        AcU[j] = d * (GDFDH[j] - GDFDL[j]) / 2.0;

    tmp0 = 1.0 / (double)nu;
    for (i = 0; i < nu; ++i)  
        AcXF[i] = tmp0;

    for (i = 0; i < nu; ++i) {
        l = AcK[i] - 1;
        k = l * ndim;
        cen = AcNS[l];
        if (cen < 0)
            continue;

        /* find observations for distributing the mass of the current one */

        n = 0;
        for (ii = i + 1; ii < nu; ++ii) {
            ll = AcK[ii] - 1;
            kk = ll * ndim;
            sflg = 1;
            for (j = 0; j < ndim; ++j) {
                if (j == cen)
                    continue;
                if (AcNS[ll] >= 0 && AcNS[ll] != cen) {
                    sflg = 0;
                    break;
                }
                if (fabs(val[kk + j] - val[k + j]) > AcU[j]) {
                    sflg = 0;
                    break;
                }
            }
            if (sflg)  
                AcR[n++] = ll;
        }
        if (n == 0)
            continue;

        tmp = AcXF[l] / (double)n;                       
        for (j = 0; j < n; ++j)  
            AcXF[AcR[j]] += tmp;
        AcXF[l] = 0.0;
    }
    if (PMProtFDef) {
        fprintf(PMProtFd,"\nDensity after redistribution.\n");
        tmp = 0.0;
        for (i = 0; i < nu; ++i) {
            fprintf(PMProtFd,"%6d ",i);
            for (j = 0; j < ndim; ++j) 
                fprintf(PMProtFd,PMTFmtS,val[i * ndim + j]);
            fprintf(PMProtFd,PMTFmtS,AcXF[i]);
            fprintf(PMProtFd,"\n");
            tmp += AcXF[i];
        }
        fprintf(PMProtFd,"Sum: %g\n",tmp);
    }

    if (opt == 0) {                         /* distribution function */
        for (i = 0; i < nu; ++i) {
            l = AcK[i] - 1;
            k = l * ndim;
            tmp = 0.0;
            for (ii = 0; ii < nu; ++ii) {
                ll = AcK[ii] - 1;
                kk = ll * ndim;
                sflg = 1;
                for (j = 0; j < ndim; ++j) {
                    if (val[kk + j] > val[k + j]) {
                        sflg = 0;
                        break;
                    }   
                }
                if (sflg)   
                    tmp += AcXF[ll];
            }
            f[l] = tmp;
        }
    }
    else if (opt == 1) {                    /* survivor function */
        for (i = 0; i < nu; ++i) {
            l = AcK[i] - 1;
            k = l * ndim;
            tmp = 0.0;
            for (ii = i + 1; ii < nu; ++ii) {
                ll = AcK[ii] - 1;
                kk = ll * ndim;
                sflg = 1;
                for (j = 0; j < ndim; ++j) {
                    if (val[kk + j] <= val[k + j]) {
                        sflg = 0;
                        break;
                    }   
                }
                if (sflg)   
                    tmp += AcXF[ll];
            }
            f[l] = tmp;
        }
    }
    else {                                  /* expected values */

        for (i = 0; i < nu; ++i) {
            l = AcK[i] - 1;
            k = l * ndim;
            cen = AcNS[l];
            if (cen < 0)
                continue;

            n = 0;
            for (ii = i + 1; ii < nu; ++ii) {
                ll = AcK[ii] - 1;
                kk = ll * ndim;
                sflg = 1;
                for (j = 0; j < ndim; ++j) {
                    if (j == cen)
                        continue;
                    if (AcNS[ll] >= 0 && AcNS[ll] != cen) {
                        sflg = 0;
                        break;
                    }
                    if (fabs(val[kk + j] - val[k + j]) > AcU[j]) {
                        sflg = 0;
                        break;
                    }
                }
                if (sflg)  
                    AcR[n++] = ll;
            }
            if (n == 0)
                continue;

            tmp0 = tmp = 0.0;
            for (j = 0; j < n; ++j) {
                ll = AcR[j];
                tmp += val[ll * ndim + cen] * AcXF[ll];
                tmp0 += AcXF[ll];
            }
            if (tmp0 > 0.0)
                val[k + cen] = tmp / tmp0;
        }   
    } 
    err = 0;

EDF4MFin:
    alloc_acr(0);
    alloc_acu(0);
    alloc_acxf(0);
    alloc_ack(0);
    alloc_acj(0);
    alloc_acns(0);
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

int lsreg1(void)
{
    register int i;
    int err,nx,r,iter;
    double tmp;

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("Regression with censored data. Current memory: %d bytes.\n",MemReq);

    MxItFlg = 0;
    MxIt1 = 10;
    TOLP = 0.001;       /* tolerance for regression parameters */
    TOLF = 0.001;       /* tolerance for surv function estimation */

    if (parm(CmdBuf + 6,4,0))     /* get parameters */
        goto LS1Fin;

    if (MxItFlg == 0)
        MxIter = 20;
    else if (MxIter < 1)
        MxIter = 1;
    if (MxIt1 < 1)
        MxIt1 = 1;
    if (PMD <= 0.0)
        PMD = 0.1;
    else if (PMD > 1.0)
        PMD = 1.0;
    if (PMFmtF == 0)
        pmfmt(10,4);

    newline();
    if (PMNI != 0) {
        if (PMNV == 0) {
            printf1("Error: need intercept or independent variables.\n");
            goto LS1Fin;
        }
        PMNI = 1;
    }
    if (PMYL < 0 || PMCEN < 0) {
        printf1("Error: need both, yl and cen parameters.\n");     
        goto LS1Fin;
    }

    if (gdf_dcheck())           /* check dependent variable */
        goto LS1Fin;            /* creates GDFYTyp: 1 exact,
                                                    2 i-censored,
                                                    3 r-censored  */
    if (PMNV > 0) {
        printf1("Indep. variables: %s",VName[PMVIdx[0]]);
        for (i = 1; i < PMNV; ++i)  
            printf1(",%s",VName[PMVIdx[i]]);
        newline();
    }
    else
        printf1("No independent variables.\n");
    newline();
        
    if (gdf_gcheck())           /* check group structure */
        goto LS1Fin;

    if (PMOPT == 1)  
        printf1("Option 1: marginal estimation");
    else {
        if (PMOPT > 3)
            PMOPT = 3;
        printf1("Option %d: joint estimation (method %d)",PMOPT,PMOPT - 1);
    }
    printf1(" of conditional expectations.\n");
    if (PMOPT > 1 && GDFGRP == 0) {
        printf1("Error: grp parameter required.\n");
        goto LS1Fin;
    }

    printf1("Max number of iterations: %d\n",MxIter);
    printf1("Tolerance for convergence: %g\n",TOLP);

    nx = PMNV;
    if (PMNI == 0)
        nx++;

    if (alloc_acy(NOC))                 /* Y vector */
        goto LS1Fin;
    if (alloc_acz(imax(NOC,nx)))        /* residuals */
        goto LS1Fin;
    if (alloc_acx(nx * nx + 1))         /* X'X */
        goto LS1Fin;
    if (alloc_acv(nx))                  /* beta */
        goto LS1Fin;
    if (alloc_acw(nx))                  /* new beta */
        goto LS1Fin;

    ls1_xx(nx,1 - PMNI);        /* create X'X */
    r = ginv(nx,nx,AcX);        /* generalized inverse */
    if (r < nx) {
        if (r < 0)
            p_err(-2,1);
        else  
            printf1("Error: X matrix is rank deficient, rank = %d\n",r);
        goto LS1Fin;
    }

    /* start with YL */

    for (i = 0; i < NOC; ++i)
        AcY[i] = get_data(PMYL,i);

    /*  get initial parameters into AcV */

    ls1_gpar(NOC,AcY,nx,1 - PMNI,AcX,AcV,AcZ);  

    if (PMProtFDef) {
        fprintf(PMProtFd,"Regression with censored data.\n");
        prvec("Initial parameters",nx,AcV - 1);
    }
    for (iter = 1; iter <= MxIter; ++iter) {

        /* calculate predicted values in AcY, residuals in AcZ */

        ls1_res(NOC,nx,1 - PMNI,AcV,AcY,AcZ);  

        if (PMOPT == 1) {                       /* marginal estimates */
            if (ls1_marg(NOC,AcZ,GDFYTyp))
                goto LS1Fin;
        }
        else {                                  /* joint estimates */
            if (ls1_joint(NOC,AcZ,GDFYTyp,PMOPT))
                goto LS1Fin;
        }

        /*  calculate new Y vector */

        for (i = 0; i < NOC; ++i) {
            if (GDFYTyp[i] == 3)
                AcY[i] += AcZ[i];              
            else
                AcY[i] = get_data(PMYL,i);
        }

        /*  get new parameters into AcW */

        ls1_gpar(NOC,AcY,nx,1 - PMNI,AcX,AcW,AcZ);  
        /* check convergence and copy new parameter estimates into AcV */

        tmp = 0.0;  
        for (i = 0; i < nx; ++i) {
            tmp = dmax(tmp,fabs(AcW[i] - AcV[i]) / dmax(AcW[i],1.0));
            AcV[i] = AcW[i];
        }
        if (PMProtFDef) {
            fprintf(PMProtFd,"\nIteration: %d\n",iter);
            prvec("Current parameters",nx,AcW - 1);
            prval("Criterion",tmp);
        }
        if (SILENTFlg < 2) {
            printfe("Iteration%3d  Crit: %20.13e\n",iter,tmp);
            fflushe();
        }   
        if (tmp < TOLP)
            break;
    }
    printf1("\nFinished after %d iterations.\n",iter);
    printf1("Final maximal parameter change: %g\n\n",tmp);     


    if (PMPPFDef)               /* write parameter to output file */
        ls1_ppar(nx,AcV);
      
    if (PMResFDef) {            /* write data and residuals to output file */
        ls1_res(NOC,nx,1 - PMNI,AcV,AcY,AcZ);  
        ls1_pres(nx,1 - PMNI,AcY,AcZ);
    }
    err = 0;

LS1Fin:
    alloc_gdfytyp(0);
    alloc_gdfgen(0,0);
    if (PM1NV > 0)                     
        vsort(0,PM1VIdx,0,0,0);
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ls1_xx(nx,iflag)   Create X'X in AcX                                    */
/*  iflag = 1 if with intercept, otherwise iflag = 0.                       */ 

void ls1_xx(int nx,int iflag)
{
    register int i,j,k;
    double tmp,tmp1,tmp2;

    for (i = 0; i < nx; ++i) {
        for (j = 0; j < nx; ++j) {

            tmp = 0.0;
            for (k = 0; k < NOC; ++k) {

                if (i == 0 && iflag)
                    tmp1 = 1.0;
                else
                    tmp1 = get_data(PMVIdx[i - iflag],k);

                if (j == 0 && iflag)
                    tmp2 = 1.0;
                else
                    tmp2 = get_data(PMVIdx[j - iflag],k);
                tmp += tmp1 * tmp2;
            }
            AcX[i * nx + j + 1] = tmp;
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

void ls1_gpar(int n,double *y,int nx,int iflag,double *xx,double *b,double *h)
{
    register int i,j;
    double tmp,tmp1;

    for (j = 0; j < nx; ++j) {
        tmp = 0.0;
        for (i = 0; i < n; ++i) {
            if (j == 0 && iflag)
                tmp1 = 1.0;
            else
                tmp1 = get_data(PMVIdx[j - iflag],i);
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

void ls1_res(int n,int nx,int iflag,double *b,double *yp,double *res)
{
    register int i,j;
    double tmp,tmp1;

    for (i = 0; i < n; ++i) {
        tmp = 0.0;
        for (j = 0; j < nx; ++j) {
            if (j == 0 && iflag)
                tmp1 = 1.0;
            else 
                tmp1 = get_data(PMVIdx[j - iflag],i);
            tmp += tmp1 * b[j];
        }
        yp[i] = tmp;
        res[i] = get_data(PMYL,i) - tmp;
    }
}

/* ------------------------------------------------------------------------ */
/*  ls1_ppar(nx,beta)       Write parameter to PMPPFd.                      */

void ls1_ppar(int nx,double *beta)
{
    register int i;

    for (i = 0; i < nx; ++i) {
        fprintf(PMPPFd,PMTFmtS,beta[i]);
        fprintf(PMPPFd,"\n");
    }
    printf1("Parameter estimates written to: %s\n",PMPPFName);
}

/* ------------------------------------------------------------------------ */
/*  ls1_pres(nx,iflag,yp,res)                                               */
/*  Write data and residuals to PMResFd.                                    */

void ls1_pres(int nx,int iflag,double *yp,double *res)
{
    register int i,j;
    double tmp;

    for (i = 0; i < NOC; ++i) {
        fprintf(PMResFd,"%6d ",i + 1);
        fprintf(PMResFd,PMFmtS,get_data(PMYL,i));
        if (PMYH >= 0)
            fprintf(PMResFd,PMFmtS,get_data(PMYH,i));
        if (PMCEN >= 0)
            fprintf(PMResFd,PMFmtS,get_data(PMCEN,i));
        fprintf(PMResFd,PMFmtS,yp[i]);
        fprintf(PMResFd,PMFmtS,res[i]);

        for (j = 0; j < nx; ++j) {
            if (j == 0 && iflag)
                tmp = 1.0;
            else
                tmp = get_data(PMVIdx[j - iflag],i);
            fprintf(PMResFd,PMFmtS,tmp);
        }
        fprintf(PMResFd,"\n");
    }
    printf1("Data and residuals written to: %s\n",PMResFName);

    if (PMTDAFDef) {
        fprintf(PMTDAFd,"# data written by lsreg1 command.\n");
        fprintf(PMTDAFd,"nvar(\n");
        fprintf(PMTDAFd,"  dfile = %s,\n",PMResFName);
        fprintf(PMTDAFd,"  noc = %d,\n",NOC);
        j = 0;
        fprintf(PMTDAFd,"  CASE [6.0] = c%-2d,\n",++j);
        fprintf(PMTDAFd,"  YL   [%d.%d] = c%-2d,\n",PMFmt1,PMFmt2,++j);
        if (PMYH >= 0)
            fprintf(PMTDAFd,"  YH   [%d.%d] = c%-2d,\n",PMFmt1,PMFmt2,++j);
        if (PMCEN >= 0)
            fprintf(PMTDAFd,"  CEN  [%d.%d] = c%-2d,\n",PMFmt1,PMFmt2,++j);
        fprintf(PMTDAFd,"  YPRED[%d.%d] = c%-2d,\n",PMFmt1,PMFmt2,++j);
        fprintf(PMTDAFd,"  RES  [%d.%d] = c%-2d,\n",PMFmt1,PMFmt2,++j);

        if (iflag)
            fprintf(PMTDAFd,"  INT  [%d.%d] = c%-2d,\n",PMFmt1,PMFmt2,++j);
        for (i = 0; i < PMNV; ++i)  
            fprintf(PMTDAFd,"  %s [%d.%d] = c%-2d,\n",VName[PMVIdx[i]],PMFmt1,PMFmt2,++j);
        fprintf(PMTDAFd,");\n");
        printf1("TDA description written to: %s\n",PMTDAFName);
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

int ls1_marg(int n,double *res,short *ytyp)
{
    register int i,j,k;
    int m;
    double tmp,tmp1,fsum,wsum;

/**
printf("res vor marg\n");
for (i = 0; i < n; ++i)
printf("i=%4d res=%g\n",i,res[i]);
**/


    if (alloc_ack(n))
        return(-1);
    if (alloc_actmp(n))
        return(-1);
           
    if (sortdp2a(n,res,ytyp,AcK,1))      /* sort */
        return(-1);    

    tmp = 1.0 / (double)n;
    m = n - 1;
    for (i = 0; i < n; ++i) {
        j = AcK[i];
        AcTmp[j] += tmp;
        if (ytyp[j] == 3 && i < n - 1) {
            tmp1 = AcTmp[j] / (double)m;
            for (k = i + 1; k < n; ++k)  
                AcTmp[AcK[k]] += tmp1;
            AcTmp[j] = 0.0;
        }
        m--;
    }
    fsum = wsum = 0.0;
    for (i = n - 1; i >= 0; --i) {
        j = AcK[i];
        if (ytyp[j] == 1 || i == n - 1) {       /* if not censored */
            fsum += res[j] * AcTmp[j];
            wsum += AcTmp[j];
            if (ytyp[j] == 1)
                res[j] = 0.0;
        }
        else if (ytyp[j] == 3) {        /* if r-censored */
            if (wsum > 0.0)
                res[j] = fsum / wsum;
        }
    }
/**
printf("res in marg\n");
for (i = 0; i < n; ++i)
printf("i=%4d res=%g\n",i,res[i]);
**/
    alloc_ack(0);
    alloc_actmp(0);
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

int ls1_joint(int n,double *res,short *ytyp,int opt)
{
    register int i,j,k,ii;
    int err,t,nu,ndim;
    double tmp;

    err = -1;
/**
    printf1("noc=%d res at begin\n",n);
    for (i = 0; i < n; ++i) {
        printf1("i=%3d res=%g\n",i,res[i] );
    }
**/


    for (t = 0; t < GDFNTyp; ++t) {     /* for all patterns */

        nu = GDFMUCnt[t];
        ndim = GDFMLen[t];

        /**  printf1("Pattern t=%d nu=%d ndim=%d\n",t,nu,ndim);  **/


        if (PMProtFDef) {
            fprintf(PMProtFd,"\nCurrent marginal pattern: %d  Typ: %d\n",t+1,GDFMTyp[t]);
            fprintf(PMProtFd,"Number of units: %d  dimensions: %d\n",nu,ndim);
        }   
        if (alloc_aci(nu * ndim + 1))   /* pointer for cases */
            goto LS1JFin;

        gdf_cptr(t,nu,ndim,AcI);        /* create pointer to data matrix cases */

        if (GDFMTyp[t] != 1) {          /* exact and censored observations */

            if (gdf_domain(t,nu,ndim,AcI,res,2,0.0))
                goto LS1JFin;

            gdf_getb(PMN,ndim,2);                   /* get number of boxes */
    
            if (alloc_act(nu * ndim + 1))           /* values/lower bounds */
                goto LS1JFin;
   
            if (opt == 2) {
                if (gdf_edf3m(t,nu,ndim,AcI,res,AcT,&tmp,2,1,MxIt1,TOLF))
                    goto LS1JFin;
            }
            else  {
                if (gdf_edf4m(t,nu,ndim,AcI,res,AcT,&tmp,PMD,2,1))
                    goto LS1JFin;
            }

            /* replace values in res[] with cond. expectations */

            for (i = 0; i < nu; ++i) {
                ii = i * ndim;
                for (j = 0; j < ndim; ++j)  {
                    k = AcI[ii + j]; 
                    res[k] = AcT[ii + j];
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
    printf1("noc=%d final\n",n);
    for (i = 0; i < n; ++i) {
        printf1("i=%3d res=%g\n",i,res[i] );
    }
**/

    err = 0; 

LS1JFin:
    alloc_act(0);
    alloc_aci(0);
    return(err);
}

