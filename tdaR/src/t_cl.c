/****************************************************************************/
/*  t_cl                                                                    */
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
#include "t_gio.h"
#include "t_rand.h"
#include "t_com.h"
#include "t_ml.h"
#include "tda_context.h"
#include "t_min.h"

/*  functions in t_cl.c */

int hcls(TDAContext *ctx);
int h_cls(TDAContext *ctx, int ctyp,int gn,int n,float *d,int *a,int *b,int *p,int *q,float *h);
float cl_getdu(TDAContext *ctx, int i,int j,float *d);
void cl_putdu(TDAContext *ctx, float x,int i,int j,float *d);
void h_cls_den(TDAContext *ctx, int n,int *a,int *b,float *h,int *p,int *q,float *x,float *y);
void h_cls_el(TDAContext *ctx, int nmax,int n,int *a,int *b,float *h,int *new,float *lev);
void h_cls_dm(TDAContext *ctx, int n,int *a,int *b,float *h,int *new,int *ll,float *d);
double h_cls_mu(TDAContext *ctx, int n,float *d);

int nncl(TDAContext *ctx);
int nn_cl(TDAContext *ctx, int *cn,double t);
int mn_cl(TDAContext *ctx, int k,int *idx,int *mptr,int *iptr,int *jptr,double *val);
int mn_cl1(TDAContext *ctx, int i,int k,int *nptr,double *val,int *sptr,int *idx);
int hcld(TDAContext *ctx);
int h_cld1(TDAContext *ctx, int n,int min,int *cidx,int *ncn,int *ncnt,int opt);
int h_cld2(TDAContext *ctx, int n,int nsmax,int *cidx,double *dia,int *nodes,             int *in,int *jn,int opt);
void h_cld_prn(TDAContext *ctx, int n,int *cidx,int cn,int k,int nn,double *dia);
double h_cld_dia(TDAContext *ctx, int n,int *cidx,int k);
int becl(TDAContext *ctx);
int be_clprn(TDAContext *ctx, int opt,int n,int *idx,double mval);
double be_cl1(TDAContext *ctx, int opt,int mflag,int i0,int n,int *pos,int *idx,int *tidx);
double be_ccheck(TDAContext *ctx, int n,int *idx);                         
double be_rcheck(TDAContext *ctx, int n,int *idx);                         
int acl(TDAContext *ctx);
float acl1(TDAContext *ctx, int n,int m,float *d,int *p,int *q,float mev,float cmin);

int udcl(TDAContext *ctx, int n,int nc,float *d,int *p,int l,int u,float *c,double *pc,int *qa,int *fp);
void udcl1(TDAContext *ctx, int n,float *d,int *p,float *c,int u);
int ucl(TDAContext *ctx);
int hclsp(TDAContext *ctx);
void hclsp_visit(TDAContext *ctx, int i,int gn);
void hclsp_prn(TDAContext *ctx, int i,int gn,int nlev,int level,int mlen);
void hclsp_prn1(TDAContext *ctx, int i,int gn,int m);
int scla(TDAContext *ctx);
int scla_check(TDAContext *ctx, int i,int n,int *acn,int *p,double *dp);
int scla_find(TDAContext *ctx, int n,int *acn,int np,int *p,int *com);
int clp(TDAContext *ctx);
int clp_f(TDAContext *ctx, int meth,int m,int n,int itmax, double *t,int *z,int *mj,double *e,double *c);
int clp_p(TDAContext *ctx);
void clp_rep(TDAContext *ctx, int m,int r,int a,int b);

int clu(TDAContext *ctx);
double clu_ld(TDAContext *ctx, double *d);
double clu_pd(TDAContext *ctx, double *d);
int clu_fn(TDAContext *ctx);
double clu_check(TDAContext *ctx);

int clpyr(TDAContext *ctx);
double clpyr_dis(TDAContext *ctx, int i,int j,float *d);
int clpyr_con(TDAContext *ctx, int i,int j,int *ack,int *acs,int *acl,int *acr);
int clpyr_ord(TDAContext *ctx, int i,int j,int n,int *idx_left,int *idx_right,int *idx_min,int *idx_max);
int clpyr_left(TDAContext *ctx, int i,int j,int *acl);
int clpyr_right(TDAContext *ctx, int i,int j,int *idx_right);
int clpyr_mcheck(TDAContext *ctx, int i,int *idx_e1,int *idx_e2,char *idx_mark);
void clpyr_mark(TDAContext *ctx, int i,int *idx_e1,int *idx_e2,char *idx_mark,int n);
int clpyr_cont(TDAContext *ctx, int i,int j,int *idx_e1,int *idx_e2);
int clpyr_ccl(TDAContext *ctx, int i,int j,int n,int nc,int *idx_e1,int *idx_e2);
void clpyr_den(TDAContext *ctx, int n,int nc,int *iord,int *idx_e1,int *idx_e2,float *acy,float *tmp);


/* global parameters */


/* -##--------------------------------------------------------------------- */
/*  hcls        Hierarchical clustering, SAHN algorithms.                   */
/*              Requires an undirected valued graph.                        */
/*                                                                          */
/*              hcls(                                                       */
/*                  opt=...,    option, def. 1                              */
/*                              1:   single-link                            */  
/*                              2:   complete-link                          */  
/*                              3:   weighted average                       */  
/*                              4:   weighted centroid                      */  
/*                              5:   group average                          */  
/*                              6:   unweighted centroid                    */  
/*                              7:   Ward's min variance                    */  
/*                  gn=...,     graph number, def. 1                        */
/*                  ptab=...,   print implied distance matrix               */
/*                  pcf=...,    create command file for dendrogram          */
/*                  df=...,     print dendrogram as edge list               */
/*                  np=...,     if np=1, add external node numbers          */
/*                  nfmt=...,   integer format, def. 4                      */
/*                  fmt=...,    print format for values, def. 10.4          */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int hcls(TDAContext *ctx)
{
    register int i,j;
    int err,nrec,nmax,n,nn,m;                 
    double d,t;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "SAHN clustering. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto HCLSFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto HCLSFin;
    if (gdd_tcheck(ctx, 0,1,2))          /* need undirected valued graph */
        goto HCLSFin;

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    if (ctx->PMNP == 1)  
        nmax = ctx->GD_NDMAX;
    else
        nmax = 0;

    if (ctx->PMOPT > 7)
        ctx->PMOPT = 7;
    printf1(ctx, "\nOption %d: ",ctx->PMOPT);
    if (ctx->PMOPT == 1)
        printf1(ctx, "single-link.\n");
    else if (ctx->PMOPT == 2)
        printf1(ctx, "complete-link.\n");
    else if (ctx->PMOPT == 3)
        printf1(ctx, "weighted average.\n");
    else if (ctx->PMOPT == 4)
        printf1(ctx, "weighted centroid.\n");
    else if (ctx->PMOPT == 5)
        printf1(ctx, "group average.\n");
    else if (ctx->PMOPT == 6)
        printf1(ctx, "unweighted centroid.\n");
    else                 
        printf1(ctx, "Ward's minimum variance.\n");
    newline(ctx);

    n  = ctx->GD_NP;
    nn = n * (n - 1) / 2;

    if (alloc_acn(ctx, n + 1))       /* a */   
        goto HCLSFin;
    if (alloc_acm(ctx, n + 1))       /* b */   
        goto HCLSFin;
    if (alloc_aci(ctx, n + 1))       /* p */   
        goto HCLSFin;
    if (alloc_acj(ctx, n + 1))       /* q */   
        goto HCLSFin;
    if (alloc_acyf(ctx, n + 1))      /* h */   
        goto HCLSFin;
    if (alloc_acxf(ctx, nn + 1))     /* d */   
        goto HCLSFin;

    if (h_cls(ctx, ctx->PMOPT,ctx->PMGN,n,ctx->AcXF,ctx->AcN,ctx->AcM,ctx->AcI,ctx->AcJ,ctx->AcYF)) {
        printf1(ctx, "Error: found negative or missing edge value.\n");
        goto HCLSFin;
    }

    /* check for monotone level function */

    m = 0;
    for (i = 2; i < n; ++i) {
        if (ctx->AcYF[i] < ctx->AcYF[i-1]) {
            m = 1;
            break;
        }
    }
    printf1(ctx, "Index function is ");
    if (m)
        printf1(ctx, "NOT ");
    printf1(ctx, "monotone.\n");

    /* create new distance matrix in AcUF */

    if (alloc_acr(ctx, 2 * n + 2))                 
        goto HCLSFin;
    if (alloc_acuf(ctx, nn + 1))                 
        goto HCLSFin;

    h_cls_dm(ctx, n,ctx->AcN,ctx->AcM,ctx->AcYF,ctx->AcI,ctx->AcR,ctx->AcUF);

    d = h_cls_mu(ctx, n,ctx->AcUF); 

    printf1(ctx, "Ultrametric condition: %lg\n",d);

    /* calculate fit of distance matrices */

    for (i = 1; i < n; ++i) {
        for (j = i + 1; j <= n; ++j) {
            t = (double)((float)gdd_adj(ctx, i - 1,j - 1,ctx->PMGN));
            cl_putdu(ctx,(float)(t),i,j,ctx->AcXF);    
        }
    }
    d = 0.0;
    for (i = 1; i <= nn; ++i) {
        t = (double)(ctx->AcXF[i] - ctx->AcUF[i]);
        d += t * t;
    }
    t = sqrt(2.0 * d);
    printf1(ctx, "Fit of distances: %lg\n\n",t);


    nrec = 0;
    for (i = 1; i < n; ++i) {
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,(double)ctx->AcYF[i]);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcN[i]);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcM[i]);
#ifdef TDA_R_PACKAGE
        tda_export_cell(ctx, "hcls.dend", (double)i);
        tda_export_cell(ctx, "hcls.dend", (double)ctx->AcYF[i]);
        tda_export_cell(ctx, "hcls.dend", (double)ctx->AcN[i]);
        tda_export_cell(ctx, "hcls.dend", (double)ctx->AcM[i]);
#endif
        if (ctx->PMNP == 1) {
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ctx->AcN[i] - 1));
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ctx->AcM[i] - 1));
#ifdef TDA_R_PACKAGE
            tda_export_cell(ctx, "hcls.dend",
                            (double)gdd_node(ctx, ctx->AcN[i] - 1));
            tda_export_cell(ctx, "hcls.dend",
                            (double)gdd_node(ctx, ctx->AcM[i] - 1));
#endif
        }
        fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
        tda_export_endrow(ctx, "hcls.dend");
#endif
        nrec++;
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);

    if (ctx->PMTabFDef) {        /* write the distance matrix */
        for (i = 1; i <= n; ++i) {
            for (j = 1; j<= n; ++j) {
                if (j > i)
                    d = (double)cl_getdu(ctx, i,j,ctx->AcUF);
                else if (j < i)
                    d = (double)cl_getdu(ctx, j,i,ctx->AcUF);
                else 
                    d = 0.0;
                rt_fprintf_d(ctx, ctx->PMTabFd,ctx->PMFmtS,d);
            }
            fprintf(ctx->PMTabFd,"\n");
        }
        printf1(ctx, "%d records written to: %s\n",n,ctx->PMTabFName);
    }

    if (ctx->PMF1Def) {                       /* write dendrogram as edge list */
        if (alloc_aczf(ctx, n + 1))  
            goto HCLSFin;
        h_cls_el(ctx, nmax,n,ctx->AcN,ctx->AcM,ctx->AcYF,ctx->AcI,ctx->AcZF);
    }    
    if (ctx->PMPCFDef) {                      /* plot dendrogram */
        if (alloc_aczf(ctx, n + 1))  
            goto HCLSFin;
        if (alloc_acxf(ctx, n + 1))  
            goto HCLSFin;
        h_cls_den(ctx, ctx->GD_NP,ctx->AcN,ctx->AcM,ctx->AcYF,ctx->AcI,ctx->AcJ,ctx->AcZF,ctx->AcXF);
    }
    err = 0;

HCLSFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  h_cls(ctyp,gn,n,d,a,b,p,q,h)                                            */
/*                                                                          */
/*              Hierarchical clustering, SAHN algorithms.                   */
/*              Adopted from H. Spaeth, Cluster-Analyse-Algorithmen,        */
/*              Muenchen 1975, S. 172.                                      */
/*                                                                          */
/*  ctyp = 1:   single-link                                                 */  
/*         2:   complete-link                                               */  
/*         3:   weighted average                                            */  
/*         4:   weighted centroid                                           */  
/*         5:   group average                                               */  
/*         6:   unweighted centroid                                         */  
/*         7:   Ward's min variance                                         */  
/*                                                                          */
/*  n = number of objects                                                   */
/*  gn = graph number                                                       */
/*  d = array of float for upper triangle of dissimilarity matrix           */
/*  a,b,p,q integer arrays of length n + 1                                  */
/*  h = array of float, length n + 1                                        */
/*                                                                          */
/*  Return 0 if OK, -1 if a negative distance found.                        */

int h_cls(TDAContext *ctx, int ctyp,int gn,int n,float *d,int *a,int *b,int *p,int *q,float *h)
{
    register int i,j,k,l;
    int ic,jc,k1,k2,qi = 0,qic = 0,qjc = 0;
    double dm,dmx,t = 0.0,f = 0.0,dj,dk;

    dm = 0.0;
    for (i = 1; i <= n; ++i) {
        p[i] = 0;
        q[i] = 1;
        for (j = i + 1; j <= n; ++j) {
            t = (double)((float)gdd_adj(ctx, i - 1,j - 1,gn));
            if (t < 0.0)  
                return(-1);
            cl_putdu(ctx,(float)(t),i,j,d);    
            dm = dmax(ctx, t,dm);
        }
    }

    ic = jc = k = 0;
    while (++k < n) {

        dmx = dm + 10.0;
        for (i = 1; i < n; ++i) {
            if (p[i] == 0) {
                for (j = i + 1; j <= n; ++j) {
                    if (p[j] == 0) {
                        t = (double)(cl_getdu(ctx, i,j,d));
                        if (t <= dmx) {
                            ic = i;
                            jc = j;
                            dmx = t;
                        }
                    }
                }
            }
        }
        p[jc] = 1;
        a[k] = ic;
        b[k] = jc;
        h[k] = (float)(dmx);

        if (ctyp >= 5) {
            qic = q[ic];
            qjc = q[jc];

            if (ctyp != 7)
                f = 1.0 / (double)((float)(qic + qjc));
        }

        for (i = 1; i <= n; ++i) {
            if (i != ic && p[i] == 0) {
                if (ctyp == 7) {
                    qi = q[i];
                    f = 1.0 / (double)((float)(qi + qic + qjc));
                }
                j  = imin(ctx, ic,i);
                l  = imax(ctx, ic,i);
                k1 = imin(ctx, jc,i);
                k2 = imax(ctx, jc,i);
                dj = (double)(cl_getdu(ctx, j,l,d));
                dk = (double)(cl_getdu(ctx, k1,k2,d));

                switch (ctyp) {
                    case 1:     t = (double)((float)dmin(ctx, (double)dj,(double)dk));
                                break; 
                    case 2:     t = (double)((float)dmax(ctx, (double)dj,(double)dk));
                                break; 
                    case 3:     t = 0.5 * (dj + dk);
                                break; 
                    case 4:     t = 0.5 * (dj + dk) - 0.25 * dmx;
                                break; 
                    case 5:     t = f * ((double)((float)(qic)) * (double)(dj) + (double)((float)(qjc)) * (double)(dk));
                                break; 
                    case 6:     t = f * ((double)((float)(qic)) * (double)(dj) + (double)((float)(qjc)) * (double)(dk) -
                                        (double)((float)qic) * (double)((float)(qjc)) * (double)(f) * dmx);
                                break; 
                    case 7:     t = f * ((double)(((float)qic + (float)qi)) * dj +
                                        (double)(((float)qjc + (float)qi)) * dk - 
                                         (double)((float)(qi)) * (double)(dmx));
                                break; 
                }
                cl_putdu(ctx,(float)(t),j,l,d);
                dm = dmax(ctx, dm,t);
            }
        }
        q[ic] += q[jc];

    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  cl_getdu(i,j,d)     return d(i,j) from upper triangle of d[].           */

float cl_getdu(TDAContext *ctx, int i,int j,float *d)
{
    (void)ctx;        /* unused: the signature is shared */
    return(d[(j - 1) * (j - 2) / 2 + i]);
}

/* ------------------------------------------------------------------------ */
/*  cl_putdu(x,i,j,d)   put x into d(i,j), upper triangle.                  */

void cl_putdu(TDAContext *ctx, float x,int i,int j,float *d)
{
    (void)ctx;        /* unused: the signature is shared */
    d[(j - 1) * (j - 2) / 2 + i] = x;
}

/* ------------------------------------------------------------------------ */
/*  h_cls_den   Create a command file for plotting a dendrogram.            */
/*              Assume data created by hcls().                              */

void h_cls_den(TDAContext *ctx, int n,int *a,int *b,float *h,int *p,int *q,float *x,float *y)
{
    register int i,j,k,l;
    int ic,jc;
    double max,d,t;

    for (i = 1; i <= n; ++i) {
        q[i] = p[i] = 0;
        x[i] = 0.0;
    }
    k = 0;
    for (i = n - 1; i >= 1; --i) {
        ic = imin(ctx, a[i],b[i]);
        jc = imax(ctx, a[i],b[i]);
        if (p[ic] == 0) {
            p[ic] = ++k;
            q[k] = ic;
        }
        else {
            l = p[ic] + 1;
            for (j = 1; j <= n; ++j) {
                if (p[j] >= l) {
                    p[j] += 1;
                    q[p[j]] = j;
                }
            }
            p[jc] = l;
            q[l] = jc;
        }
        if (p[jc] == 0) {
            p[jc] = ++k;
            q[k] = jc;
        }
    }
    max = (double)(h[n - 1]);
    d = max / 15.0;
    fprintf(ctx->PMPCFd,"psfile = %s.ps;\n",ctx->PMPCFName);
    fprintf(ctx->PMPCFd,"psetup(pxa=%g,%g,pya=0,%d);\n",-2.0 * d,max + d,n+1);
    for (i = 1; i <= n; ++i) {
        if (ctx->PMNP == 1)
            fprintf(ctx->PMPCFd,"pltext(xy=%g,%d)=%d;\n",-d,i,gdd_node(ctx, q[i] - 1));
        else
            fprintf(ctx->PMPCFd,"pltext(xy=%g,%d)=%d;\n",-d,i,q[i]);
    }
    for (i = 1; i <= n; ++i)  
        y[q[i]] = (float)i;

    for (i = 1; i < n; ++i) {
        ic = a[i];
        jc = b[i];
        t = (double)(h[i]);
        fprintf(ctx->PMPCFd,"plotp=%g,%g,%g,%g,%g,%g,%g,%g;\n",    
                          (double)(x[ic]),(double)(y[ic]),t,(double)(y[ic]),t,(double)(y[jc]),(double)(x[jc]),(double)(y[jc]));
        x[ic] = x[jc] = (float)t;
        y[ic] = y[jc] = (float)(((double)y[ic] + (double)y[jc]) / 2.0);
    }
    printf1(ctx, "Output file for plot of dendrogram: %s\n",ctx->PMPCFName);
}

/* ------------------------------------------------------------------------ */
/*  h_cls_el    Write dendrogram as an edge list to PMF1d.                  */
/*              Assume data created by hcls().                              */
/*                                                                          */
/*  If nmax = 0 begin counting the internal nodes with GD_NP + 1, otherwise */
/*  begin with nmax + 1.                                                    */

void h_cls_el(TDAContext *ctx, int nmax,int n,int *a,int *b,float *h,int *new,float *lev)
{
    register int i,ic,jc;
    int nrec,nn,n1,n2,m1,m2;
    double l1,l2;

    for (i = 1; i <= n; ++i) {
        new[i] = i;
        lev[i] = 0.0;
    }
    if (nmax == 0)
        nn = ctx->GD_NP + 1;
    else
        nn = nmax + 1;

    nrec = 0;
    for (i = 1; i < n; ++i) {
        ic = a[i];
        jc = b[i];
        if (new[ic] <= new[jc]) {
            n1 = new[ic];
            n2 = new[jc];
            l1 = (double)(h[i] - lev[ic]);
            l2 = (double)(h[i] - lev[jc]);
        }
        else {
            n1 = new[jc];
            n2 = new[ic];
            l1 = (double)(h[i] - lev[jc]);
            l2 = (double)(h[i] - lev[ic]);
        }
        if (nmax == 0)
            m1 = n1;
        else {
            if (n1 <= ctx->GD_NP)   
                m1 = gdd_node(ctx, n1 - 1);
            else
                m1 = n1;
        }
        if (nmax == 0)
            m2 = n2;
        else {
            if (n2 <= ctx->GD_NP)   
                m2 = gdd_node(ctx, n2 - 1);
            else
                m2 = n2;
        }
        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,m1);
        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,nn);
        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,l1);
        fprintf(ctx->PMF1d,"\n");
        nrec++;
        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,m2);
        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,nn);
        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,l2);
        fprintf(ctx->PMF1d,"\n");
        nrec++;
        new[ic] = new[jc] = nn;
        lev[ic] = lev[jc] = h[i];
        nn++;
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);
}

/* -##--------------------------------------------------------------------- */
/*  h_cls_dm    Create the distance matrix implied by the dendrogram        */
/*              in the array d[]; only upper triangle.                      */
/*                                                                          */

void h_cls_dm(TDAContext *ctx, int n,int *a,int *b,float *h,int *new,int *ll,float *d)
{
    register int i,j,ic,jc;
    int nn,n1,n2;          

    for (i = 1; i <= n; ++i)   
        new[i] = i;
    nn = n + 1;
    for (i = 1; i < n; ++i) {
        ic = a[i];
        jc = b[i];
        if (new[ic] <= new[jc]) {
            n1 = new[ic];
            n2 = new[jc];
        }
        else {
            n1 = new[jc];
            n2 = new[ic];
        }
        ll[n1] = nn;
        ll[n2] = nn;
        new[ic] = new[jc] = nn;
        nn++;
    }
    for (i = 1; i < n; ++i) {
        for (j = i + 1; j <= n; ++j) {
            ic = i;
            jc = j;
            while (1) {
                while (ic < jc)
                    ic = ll[ic];
                while (jc < ic)
                    jc = ll[jc];
                if (ic == jc)
                    break;
            }
            cl_putdu(ctx, h[ic - n],i,j,d);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  h_cls_mu    Return measure of ultrametric condition.                    */
/*              Assume upper triangle of distance matrix in d.              */
/*                                                                          */

double h_cls_mu(TDAContext *ctx, int n,float *d)
{
    register int i,j,k;               
    int ii,jj,kk;
    double t,dij,dik,djk;
         
    t = 0.0;
    for (i = 1; i < n; ++i) {
        ii = (i - 1) * (i - 2) / 2;
        for (j = i + 1; j <= n; ++j) {
            jj = (j - 1) * (j - 2) / 2;
            dij = (double)(d[jj + i]);
            for (k = 1; k <= n; ++k) {
                if (k == i || k == j)
                    continue;

                kk = (k - 1) * (k - 2) / 2;

                if (k < i)
                    dik = (double)(d[ii + k]);
                else
                    dik = (double)(d[kk + i]);
                                           
                if (dij > dik)
                    continue;

                if (k < j)
                    djk = (double)(d[jj + k]);
                else
                    djk = (double)(d[kk + j]);

                if (dij <= djk && fabs(dik - djk) > ctx->EPSI1)
                    t += (dik - djk) * (dik - djk);
            }
        }
    }
    return((double)t);
}

/* ------------------------------------------------------------------------ */
/*  nncl        Nearest-neighbor clustering.                                */
/*              See Jain and Dubes, p. 128.                                 */
/*                                                                          */
/*              Requires an undirected valued graph.                        */
/*                                                                          */
/*              nncl(                                                       */
/*                  alg = ...,  algorithm, def. 1                           */
/*                              1 = nearest neighbors                       */
/*                              2 = mutual neighborhood clustering          */
/*                  gn = ...,   graph number, def. 1                        */
/*                  sc = ...,   threshold, def. 1                           */
/*                  ns=...,     size of neighborhood, def. 2                */
/*                  pcf=...,    create command file for dendrogram          */
/*                  df=...,     print dendrogram as edge list               */
/*                  nfmt = ..., integer format, def. 4                      */
/*                  fmt = ...,  print format for values, def. 10.4          */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int nncl(TDAContext *ctx)
{
    register int i;
    int err,nrec,n;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Nearest-neighbor clustering. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto NNCLFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto NNCLFin;
    if (gdd_tcheck(ctx, 0,1,2))          /* need undirected valued graph */
        goto NNCLFin;

    if (ctx->PMALG != 2)
        ctx->PMALG = 1;

    printf1(ctx, "Algorithm %d: ",ctx->PMALG);
    if (ctx->PMALG == 1) {
        if (ctx->PMSC <= 0.0)
            ctx->PMSC = 1.0;
        printf1(ctx, "nearest neighbors with threshold %g\n",ctx->PMSC);
    }
    else {
        if (ctx->PMNS < 2)
            ctx->PMNS = 2;
        printf1(ctx, "mutual neighborhoods of size %d\n",ctx->PMNS);
    }
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    nrec = 0;
    if (ctx->PMALG == 1) {
        if (alloc_acn(ctx, ctx->GD_NP))           /* cluster numbers */
            goto NNCLFin;

        n = nn_cl(ctx, ctx->AcN,ctx->PMSC);

        printf1(ctx, "Found %d clusters.\n",n);
        for (i = 0; i < ctx->GD_NP; ++i) {
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcN[i]);
            fprintf(ctx->PMFd,"\n");
            nrec++;
        }
    }
    else {
        if (alloc_acn(ctx, ctx->GD_NP * ctx->PMNS + 1)) 
            goto NNCLFin;
        if (alloc_acr(ctx, ctx->GD_NP * ctx->PMNS + 1)) 
            goto NNCLFin;
        if (alloc_aci(ctx, ctx->GD_NP * ctx->PMNS + 1)) 
            goto NNCLFin;
        if (alloc_acj(ctx, ctx->GD_NP * ctx->PMNS + 1)) 
            goto NNCLFin;
        if (alloc_acx(ctx, ctx->GD_NP + 1)) 
            goto NNCLFin;

        nrec = mn_cl(ctx, ctx->PMNS,ctx->AcN,ctx->AcR,ctx->AcI,ctx->AcJ,ctx->AcX);
        if (nrec < 0)
            goto NNCLFin;
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

NNCLFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  nn_cl       Nearest-neighbor clustering.                                */
/*              Return: 0 if OK, -1 if error.                               */

int nn_cl(TDAContext *ctx, int *cn,double t)
{
    register int i,j,k;
    int m;
    double d;

    m = 1;                     /* number of clusters */
    cn[0] = 1;

    for (i = 1; i < ctx->GD_NP; ++i) {
        if (ctx->SILENTFlg < 2 && ctx->GD_NP >= 1000) {
            printfe(ctx, "Node: %7d     %c",i + 1,CR);
            fflushe(ctx);
        }
        k = -1;
        for (j = 0; j < i; ++j) {
            d = gdd_adj(ctx, i,j,ctx->PMGN);
            if (d >= 0.0 && d <= t) {
                k = cn[j];
            }
        }
        if (k < 0)  
            cn[i] = ++m;
        else   
            cn[i] = k;
    }
    if (ctx->SILENTFlg < 2 && ctx->GD_NP >= 1000) {
        printfe(ctx, "\n");
        fflushe(ctx);
    }
    return(m);
}

/* ------------------------------------------------------------------------ */
/*  mn_cl       Mutual neighborhood clustering.                             */
/*              Return number of records written to output file, or -1 if   */
/*              error (= insufficient memory).                              */  

int mn_cl(TDAContext *ctx, int k,int *idx,int *mptr,int *iptr,int *jptr,double *val)
{
    register int i,j,l;
    int m,n,ip,jp,n1,ip1,jp1,nrec;

    for (i = 0; i < ctx->GD_NP; ++i) {
        if (ctx->SILENTFlg < 2 && ctx->GD_NP >= 1000) {
            printfe(ctx, "Node: %7d     %c",i + 1,CR);
            fflushe(ctx);
        }
        n = mn_cl1(ctx, i,k,mptr,val,iptr,jptr); 
        if (n < 0)
            return(-1);

        for (j = 0; j < k; ++j) {
            if (j < n)
                idx[i * k + j] = mptr[j];
            else
                idx[i * k + j] = ctx->INTMAX;
        }
    }
    if (ctx->SILENTFlg < 2 && ctx->GD_NP >= 1000) {
        printfe(ctx, "\n");
        fflushe(ctx);
    }
    m = 0;
    for (i = 0; i < ctx->GD_NP; ++i) {
        for (j = 0; j < k; ++j) {
            ip = idx[i * k + j];
            if (ip >= ctx->GD_NP)
                continue;
            n = 0;
            for (l = 0; l < k; ++l) {
                if (idx[ip * k + l] == i) {
                    n = 2 + j + l;
                    break;
                }
            }
            if (n) {
                mptr[m] = n;
                iptr[m] = imin(ctx, i,ip);
                jptr[m++] = imax(ctx, i,ip);
            }
        }
    }
    if (sortdpi2(ctx, m,mptr,iptr,idx))
        return(-1);

    nrec = 0;
    n1 = ip1 = jp1 = -1;
    for (i = 0; i < m; ++i) {
        j = idx[i];
        n  = mptr[j];
        ip = iptr[j];
        jp = jptr[j];
        if (n != n1 || ip != ip1 || jp != jp1) {
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nrec + 1);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,n);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ip + 1);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,jp + 1);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ip));
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, jp));
            fprintf(ctx->PMFd,"\n");
            nrec++;
        }
        n1 = n;
        ip1 = ip;
        jp1 = jp;
    }
    return(nrec);
}

/* ------------------------------------------------------------------------ */
/*  mn_cl1.         Get nearest neighbors into nptr, for node i.            */
/*                  Return number of neighbors, or -1 if error.             */

int mn_cl1(TDAContext *ctx, int i,int k,int *nptr,double *val,int *sptr,int *idx)
{
    register int j;
    int n;
    double tmp;

    n = 0;
    for (j = 0; j < ctx->GD_NP; ++j) {
        if (j != i) {
            tmp = gdd_adj(ctx, i,j,ctx->PMGN);
            if (tmp >= 0.0) {
                idx[n] = j;
                val[n++] = tmp;
            }
        }
    }
    if (n < 1)
        return(0);

    if (sortdp(ctx, n,val,sptr))
        return(-1);

    n = imin(ctx, n,k);
    for (j = 0; j < n; ++j)
        nptr[j] = idx[sptr[j]];

    return(n);
}

/* ------------------------------------------------------------------------ */
/*  hcld        Hierarchical divisive clustering.                           */
/*              Requires an undirected valued graph.                        */
/*              Graph must be connected.                                    */
/*                                                                          */
/*              hcld(                                                       */
/*                  alg = ...,  algorithm, def. 1                           */
/*                              1:  simple binary splits                    */  
/*                              2:  min diameter clustering (Prim)          */
/*                  opt=...,    output option, def. 1                       */
/*                              1 = write all clusters                      */
/*                              2 = write only final clusters               */
/*                  gn = ...,   graph number, def. 1                        */
/*                  min=...,    min size of clusters, def. 2    alg 1       */
/*                  max=...,    max number of splits, def. 1    alg 2       */
/*                  df=...,     print dendrogram as edge list               */
/*                  nfmt = ..., integer format, def. 4                      */
/*                  fmt = ...,  print format for values, def. 10.4          */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int hcld(TDAContext *ctx)
{
    int err,r,cmax;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Hierarchical divisive clustering. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto HCLDFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto HCLDFin;
    if (gdd_tcheck(ctx, 0,1,2))          /* need undirected valued graph */
        goto HCLDFin;

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);
    if (ctx->PMMin < 2)
        ctx->PMMin = 2;
    if (ctx->PMMax < 1)
        ctx->PMMax = 1;

    if (ctx->PMOPT > 2)
        ctx->PMOPT = 2;

    if (ctx->PMALG != 2)
        ctx->PMALG = 1;

    printf1(ctx, "Algorithm %d: ",ctx->PMALG);
    if (ctx->PMALG == 1)
        printf1(ctx, "maximal distance splits.\n");
    else                   
        printf1(ctx, "minimal diameter splits.\n");

    if (alloc_acn(ctx, ctx->GD_NP + 1))      
        goto HCLDFin;
    if (alloc_acm(ctx, ctx->GD_NP + 1))      
        goto HCLDFin;
    if (alloc_aci(ctx, ctx->GD_NP + 1))      
        goto HCLDFin;

    if (ctx->PMALG == 2) {
        cmax = 2 * ctx->PMMax + 1;           /* max number of clusters */

        if (alloc_acu(ctx, cmax + 2))      
            goto HCLDFin;
        if (alloc_acj(ctx, ctx->GD_NP + 1))      
            goto HCLDFin;
    }

    if (ctx->PMALG == 1)
        r = h_cld1(ctx, ctx->GD_NP,ctx->PMMin,ctx->AcN,ctx->AcM,ctx->AcI,ctx->PMOPT);
    else
        r = h_cld2(ctx, ctx->GD_NP,ctx->PMMax,ctx->AcN,ctx->AcU,ctx->AcM,ctx->AcI,ctx->AcJ,ctx->PMOPT);

    if (r < 0)
        goto HCLDFin;

    err = 0;

HCLDFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  h_cld1      Binary splits according to largest distance.                */
/*                                                                          */
/*  n = number of nodes                                                     */
/*  min = min size of clusters                                              */
/*  cidx[] number of cluster assigned to node i                             */
/*  ncn[] current cluster size                                              */
/*  ncnt[] new cluster size.                                                */
/*  if opt == 1 write all clusters, if opt = 2 only final clusters.         */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int h_cld1(TDAContext *ctx, int n,int min,int *cidx,int *ncn,int *ncnt,int opt)
{
    register int i,j,k,l;
    int fin,nc,ncl,nch,nnc,i0,j0,kn,nr,nr1;
    double dmax,d,di,dj;
       
    nr = nr1 = 0;
    nc = 1;    
    ncl = 1;
    nch = 2;
    ncn[0] = n;
    for (k = 0; k < n; ++k)
        cidx[k] = 1;
       
    fin = 0;
    while (fin == 0) {
        fin = 1;
        nnc = 0;

        for (k = 0; k < n; ++k)
            ncnt[k] = 0;

        for (k = 0; k < nc; ++k) {      /* for all clusters */

            kn = ncl + k;

            if (ncn[k] <= min || opt == 1) {

                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,kn);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ncn[k]);
                for (l = 0; l < n; ++l) {
                    if (cidx[l] == kn)
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,l + 1);
                }
                fprintf(ctx->PMFd,"\n");
                nr++;
            }
            if (ncn[k] <= min)   
                continue;

            /* split cluster k */

            dmax = 0.0;
            i0 = j0 = -1;
            for (i = 0; i < n; ++i) {
                if (cidx[i] != kn)
                    continue;
                for (j = i + 1; j < n; ++j) {
                    if (cidx[j] != kn)
                        continue;
                    d = gdd_adj(ctx, i,j,ctx->PMGN);
                    if (d >= dmax) {
                        dmax = d;
                        i0 = i;
                        j0 = j;
                    }
                }
            }
            if (i0 < 0)  
                continue;

            fin = 0;

            /* create two new clusters */

            for (i = 0; i < n; ++i) {
                if (cidx[i] != kn)
                    continue;
                di = gdd_adj(ctx, i0,i,ctx->PMGN);
                dj = gdd_adj(ctx, j0,i,ctx->PMGN);
                if (di <= dj) {
                    cidx[i] = nch + nnc;
                    ncnt[nnc] += 1;
                }
                else {
                    cidx[i] = nch + nnc + 1;
                    ncnt[nnc + 1] += 1;
                }
            }
            if (ctx->PMF1Def) {   /* write split tree as edge list */

                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,nch + nnc);
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,kn);
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ncnt[nnc]);
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ncn[k]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,dmax / 2.0);
                fprintf(ctx->PMF1d,"\n");
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,nch + nnc + 1);
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,kn);
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ncnt[nnc + 1]);
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ncn[k]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,dmax / 2.0);
                fprintf(ctx->PMF1d,"\n");
                nr1 += 2;
            }
            nnc += 2;
        }
        for (k = 0; k < nnc; ++k)  
            ncn[k] = ncnt[k];
        ncl += nc;      
        nch += nnc;
        nc = nnc;
    }
    printf1(ctx, "%d records written to: %s\n",nr,ctx->PMFdName);
    if (nr1 > 0)
        printf1(ctx, "%d records written to: %s\n",nr1,ctx->PMF1dName);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  h_cld2      Binary splits according to minimal diameter.                */
/*                                                                          */
/*  n = number of nodes                                                     */
/*  nsmax = max number of splits, max number of clusters is 2 nsmax + 1     */
/*  cidx[] number of cluster assigned to node i                             */
/*  cptr[k], k=0,...,ncmax-1 pointer to cluster indices                     */
/*  dia[k],  k=0,...,ncmax-1 diameters                                      */
/*                                                                          */
/*  if opt == 1 write all clusters, if opt = 2 only final clusters.         */
/*  MST always with Prim's algorithm.                                       */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int h_cld2(TDAContext *ctx, int n,int nsmax,int *cidx,double *dia,int *nodes,              int *in,int *jn,int opt)
{
    register int i,j,k;
    int err,nr,nr1,nc,cn,ne,m,ns,ns1,ns2,n1,n2;
    double tmp;
       
    nr = nr1 = 0;
    err = -1;

    for (i = 0; i < n; ++i)
        cidx[i] = 1;

    nc = 1;
    dia[1] = h_cld_dia(ctx, n,cidx,1);

    for (ns = 1; ns <= nsmax; ++ns) {

        /* select cluster cn with maximal diameter */

        tmp = -1.0;
        cn = -1;
        for (k = 1; k <= nc; ++k) {
            if (tmp < dia[k]) {
                cn = k;
                tmp = dia[k];
            }   
        }
        if (cn < 0)
            break;

        /* select nodes for cluster cn */

        m = 0;
        for (i = 0; i < n; ++i) {
            if (cidx[i] == cn)  
                nodes[m++] = i;
        }
        if (m <= 1)  
            break;          

        if (m == 2) {
            in[0] = nodes[0];
            jn[0] = nodes[1];
            ne = 1;
        }
        else {              /* create max spanning tree */

            /* note: the mst function uses AcX and AcR          
               and destroys nodes[] */

            ne = g_mst(ctx, ctx->PMGN,m,nodes,in,jn,1);
            if (ne <= 0) {
                if (ne == 0)  
                    printf1(ctx, "ERROR: hcld2.\n");
                goto HCLD2Fin;
            }   
        }
        n1 = n2 = 1;        /* number of cluster elements */
        i = in[0];
        j = jn[0];
        ns1 = 2 * ns;       /* bicoloring */
        ns2 = ns1 + 1;
        cidx[i] = ns1;
        cidx[j] = ns2;

        for (k = 1; k < ne; ++k) {
            i = in[k];
            j = jn[k];
            if (cidx[i] == ns1) {
                cidx[j] = ns2;
                n2++;
            }
            else if (cidx[i] == ns2) {
                cidx[j] = ns1;
                n1++;
            }
            else if (cidx[j] == ns1) {
                cidx[i] = ns2;
                n2++;
            }
            else if (cidx[j] == ns2) {
                cidx[i] = ns1;
                n1++;
            }
            else {
                printf1(ctx, "ERROR1: hcld2.\n");
                goto HCLD2Fin;
            }
        }

        /* find diameters of new clusters */

        dia[2 * ns] = h_cld_dia(ctx, n,cidx,2 * ns);
        dia[2 * ns + 1] = h_cld_dia(ctx, n,cidx,2 * ns + 1);

        /* print */

        if (opt == 1) {
            h_cld_prn(ctx, n,cidx,cn,ns1,n1,dia);
            h_cld_prn(ctx, n,cidx,cn,ns2,n2,dia);
            nr += 2;
        }
        if (ctx->PMF1Def) {   /* write split tree as edge list */

            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ns1);
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,cn);
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,n1);
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,m);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,dia[ns1]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,dia[cn]);
            fprintf(ctx->PMF1d,"\n");

            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ns2);
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,cn);
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,n2);
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,m);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,dia[ns2]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,dia[cn]);
            fprintf(ctx->PMF1d,"\n");
            nr1 += 2;
        }
        dia[cn] = -1.0;
        nc += 2;
    }
    printf1(ctx, "%d records written to: %s\n",nr,ctx->PMFdName);
    if (nr1 > 0)
        printf1(ctx, "%d records written to: %s\n",nr1,ctx->PMF1dName);

    err = 0;
         
HCLD2Fin:
    return(err);
}
   
/* ------------------------------------------------------------------------ */
/*  h_cld_prn()             print cluster.                                  */

void h_cld_prn(TDAContext *ctx, int n,int *cidx,int cn,int k,int nn,double *dia)
{
    register int i;

    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,k);
    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,cn);
    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,dia[k]);
    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nn);
    for (i = 0; i < n; ++i) {
        if (cidx[i] == k)
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
    }
    fprintf(ctx->PMFd,"\n");
}

/* ------------------------------------------------------------------------ */
/*  h_cld_dia(n,cidx,k)     Return diameter of subgraph with nodes in       */
/*                          cidx[i] = k.                                    */

double h_cld_dia(TDAContext *ctx, int n,int *cidx,int k)
{
    register int i,j;
    double dia,tmp;

    dia = 0.0;

    for (i = 0; i < n; ++i) {
        if (cidx[i] != k)
            continue;

        for (j = 0; j < n; ++j) {
            if (j == i || cidx[j] != k)
                continue;

            tmp = gdd_adj(ctx, i,j,ctx->PMGN);
            if (tmp >= 0.0)  
                dia = dmax(ctx, dia,tmp);
        }
    }
    return(dia);
}

/* ------------------------------------------------------------------------ */
/*  becl        Clustering with bond energy algorithm.                      */
/*                                                                          */
/*              becl(                                                       */
/*                  alg = ...,  option, def. 1                              */
/*                              1 column permutations                       */
/*                              2 row permutations                          */
/*                              3 row and column permutations               */
/*                  min = ...,  1 minimize criterion, def. maximize         */
/*                  opt = ...,  output option, def. 1                       */
/*                              1 only node list                            */
/*                              2 edge list, dep. on directed/undir.        */
/*                              3 lower triangle                            */
/*                              4 square matrix                             */
/*                              5 square matrix plus leading columns        */
/*                  gn = ...,   graph number, def. 1                        */
/*                  cn=...,     list of starting nodes, def. 1              */
/*                  sc=...,     substitute for miss values, def. -1         */
/*                  nfmt = ..., integer format, def. 4                      */
/*                  fmt = ...,  print format for values, def. 10.4          */
/*                  sc=..,      substitute for missings if opt = 4,5        */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int becl(TDAContext *ctx)
{
    register int i,j;
    int err,nrec,i0,cp,isel;
    double tmp,dmax;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "BE clustering. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto BECLFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto BECLFin;
    if (gdd_tcheck(ctx, 0,0,0))
        goto BECLFin;

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);
    if (ctx->PMSCFlg == 0)
        ctx->PMSC = -1.0;
    if (ctx->PMOPT < 1 || ctx->PMOPT > 5)
        ctx->PMOPT = 1;

    if (ctx->PMMin == 1) {
        printf1(ctx, "Minimize with ");
        dmax = ctx->DBLMAX;
    }
    else {
        ctx->PMMin = 0;
        printf1(ctx, "Maximize with ");
        dmax = -1.0;
    }
    if (ctx->PMALG < 1 || ctx->PMALG > 3)
        ctx->PMALG = 1;

    if (ctx->PMALG == 1) 
        printf1(ctx, "column");
    else if (ctx->PMALG == 2) 
        printf1(ctx, "row");
    else              
        printf1(ctx, "row and column");
    printf1(ctx, " permutations.\n");

    if (alloc_acn(ctx, ctx->GD_NP))      
        goto BECLFin;
    if (alloc_aci(ctx, ctx->GD_NP))      
        goto BECLFin;
    if (alloc_acj(ctx, ctx->GD_NP))      
        goto BECLFin;
    if (alloc_acs(ctx, ctx->GD_NP))      
        goto BECLFin;
   
    printf1(ctx, "\nStarting node  final BE criterion\n");

    isel = -1;
    i0 = -1;
    cp = 0;
    while (cp < ctx->PMNCN) {
        j = ctx->PMCN[cp++] - 1;
        if (j >= 0 && j < ctx->GD_NP) {
            i0 = j;
            break;
        }
    }
    if (i0 < 0)
        i0 = 0;

    while (1) {
        if (ctx->PMALG == 1 || ctx->PMALG == 3) {
            tmp = be_cl1(ctx, 1,ctx->PMMin,i0,ctx->GD_NP,ctx->AcN,ctx->AcI,ctx->AcJ);
            if ((ctx->PMMin == 0 && dmax < tmp) || (ctx->PMMin && dmax > tmp)) {
                dmax = tmp;
                isel = i0;
                for (i = 0; i < ctx->GD_NP; ++i)
                    ctx->AcS[i] = ctx->AcI[i];
            }
        }
        if (ctx->PMALG == 2 || ctx->PMALG == 3) {
            tmp = be_cl1(ctx, 2,ctx->PMMin,i0,ctx->GD_NP,ctx->AcN,ctx->AcI,ctx->AcJ);
            if ((ctx->PMMin == 0 && dmax < tmp) || (ctx->PMMin && dmax > tmp)) {
                dmax = tmp;
                isel = i0;
                for (i = 0; i < ctx->GD_NP; ++i)
                    ctx->AcS[i] = ctx->AcI[i];
            }
        }
        printf1(ctx, "%11d  %18.4f\n",i0 + 1,dmax);

        i0 = -1;
        while (cp < ctx->PMNCN) {
            j = ctx->PMCN[cp++] - 1;
            if (j >= 0 && j < ctx->GD_NP) {
                i0 = j;
                break;
            }
        }
        if (i0 < 0)
            break;
    }
    printf1(ctx, "\nSelected for output: starting node %d.\n",isel + 1);
    nrec = be_clprn(ctx, ctx->PMOPT,ctx->GD_NP,ctx->AcS,ctx->PMSC);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

BECLFin:       
    p_clean(ctx);
    return(err);
}
   
/* ------------------------------------------------------------------------ */
/*  be_clprn(opt,n,idx)     print to output file, idx[] contains permut.    */
/*                          opt = 1 only node list                          */
/*                                2 edge list, dep. on directed/undir.      */
/*                                3 lower triangle                          */
/*                                4 square matrix                           */
/*                                5 square matrix plus leading columns      */
/*                                                                          */
/*  Return number of records written to output file.                        */  

int be_clprn(TDAContext *ctx, int opt,int n,int *idx,double mval)
{
    register int i,j;
    int ip,jp,m = 0;
    double tmp;

    if (opt == 1) {
        for (i = 0; i < n; ++i) {
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,idx[i] + 1);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, idx[i]));
            fprintf(ctx->PMFd,"\n");
            m++;
        }
    }
    else if (opt == 2) {
        if (ctx->GD_GT <= 2) {               /* undirected */
            for (i = 0; i < n; ++i) {
                ip = idx[i];
                for (j = i; j < n; ++j) {
                    jp = idx[j];
                    tmp = gdd_adj(ctx, ip,jp,ctx->PMGN);
                    if (tmp >= 0.0) {
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j + 1);
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ip + 1);
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,jp + 1);
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ip));
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, jp));
                        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
                        fprintf(ctx->PMFd,"\n");
                        m++;
                    }
                }
            }
        }   
        else {                             /* directed */
            for (i = 0; i < n; ++i) {
                ip = idx[i];
                for (j = 0; j < n; ++j) {
                    jp = idx[j];
                    tmp = gdd_adj(ctx, ip,jp,ctx->PMGN);
                    if (tmp >= 0.0) {
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j + 1);
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ip + 1);
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,jp + 1);
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ip));
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, jp));
                        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
                        fprintf(ctx->PMFd,"\n");
                        m++;
                    }
                }
            }
        }   
    }
    else if (opt == 3) {            /* lower triangle */

        for (i = 1; i < n; ++i) {
            ip = idx[i];
            for (j = 0; j < i; ++j) {
                jp = idx[j];
                tmp = gdd_adj(ctx, ip,jp,ctx->PMGN);
                if (tmp < 0.0)
                    tmp = mval;

                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j + 1);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ip + 1);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,jp + 1);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ip));
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, jp));
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
                fprintf(ctx->PMFd,"\n");
                m++;
            }
        }
    }   
    else {                      /* square matrix */

        for (i = 0; i < n; ++i) {
            ip = idx[i];

            if (opt == 5) {
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ip + 1);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ip));
            }
            for (j = 0; j < n; ++j) {
                jp = idx[j];
                tmp = gdd_adj(ctx, ip,jp,ctx->PMGN);
                if (tmp < 0.0)
                    tmp = mval;
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
            }
            fprintf(ctx->PMFd,"\n");
            m++;
        }   
    }
    return(m);
}

/* ------------------------------------------------------------------------ */
/*  be_cl1(opt,i0,n,pos,idx,tidx)                                           */
/*                                                                          */
/*  Calculates BE measure. Return optimal permutation in idx[].             */
/*  Begin with starting node i0.                                            */
/*  opt = 1 : column permutations                                           */
/*  opt = 2 : row permutations                                              */
/*  if mflag != 0 minimize the criterion                                    */

double be_cl1(TDAContext *ctx, int opt,int mflag,int i0,int n,int *pos,int *idx,int *tidx)
{
    register int i = 0,j = 0,k = 0;
    int nidx = 0,knum = 0,kdir = 0,kidx = 0,ktdir = 0,ktidx = 0,pos0 = 0,nt = 0,fin = 0,nn = 0;
    double bmax = 0.0,be = 0.0;

    for (i = 0; i < n; ++i)  
        pos[i] = -2;

    bmax = 0.0;
    pos0 = i0;   
    pos[pos0] = -1;
    idx[0] = i0;
    nidx = 1;

    while (nidx < n) {
     
        if (mflag == 0)
            bmax = -1.0;
        else
            bmax = ctx->DBLMAX;

        for (k = 0; k < n; ++k) {
            if (pos[k] >= -1)
                continue;

            /* try column k */

            j = pos0;
            fin = i = 0;
            while (fin == 0) {
                nt = 0;
                j = pos0;
                nn = 0;
                while (j >= 0) {
                    if (i == nn) {
                        if (nt == 0) {
                            ktdir = 0;
                            ktidx = j;
                        }
                        else {          
                            ktdir = 1;
                            ktidx = tidx[nt - 1];
                        }
                        tidx[nt++] = k;
                    }
                    tidx[nt++] = j;
                    j = pos[j];
                    nn++;
                }
                if (i == nn) {
                    ktdir = 1;
                    ktidx = tidx[nt - 1];
                    tidx[nt++] = k;
                    fin = 1;
                }     
                i++;

                if (opt == 1)
                    be = be_ccheck(ctx, nt,tidx);
                else
                    be = be_rcheck(ctx, nt,tidx);

                if ((mflag == 0 && bmax < be) || (mflag && bmax > be)) {
                    knum = k;
                    kdir = ktdir;
                    kidx = ktidx;  
                    bmax = be;
                }
            }
        }
        if (kdir == 0) {        /* insert at beginning */
            pos[knum] = pos0;
            pos0 = knum;
        }
        else {                  /* insert behind kidx */
            j = pos[kidx];
            pos[kidx] = knum;
            pos[knum] = j;
        }
        idx[nidx++] = knum;    
    }
    i = 0;                  /* return permutation in idx */
    j = pos0;
    while (j >= 0) {
        idx[i++] = j;
        j = pos[j];
    }
    /**
    printf1(ctx, "IDX: ");
    for (i = 0; i < n; ++i)
        printf1(ctx, "%d ",idx[i]);
    newline(ctx);
    **/
    return(bmax);
}

/* ------------------------------------------------------------------------ */
/*  be_ccheck(n,idx)                                                        */
/*                                                                          */
/*  Calculate BE measure for columns with nodes in idx[i], i=0,...,n-1.     */

double be_ccheck(TDAContext *ctx, int n,int *idx)                          
{
    register int i,j,k;
    double be,tmp,tmp1;

    be = 0.0;
    for (i = 0; i < ctx->GD_NP; ++i) {
        for (j = 0; j < n; ++j) {
            k = idx[j];
            tmp = gdd_adj(ctx, i,k,ctx->PMGN);
            if (tmp <= 0.0)
                continue;

            if (j > 0 && (tmp1 = gdd_adj(ctx, i,idx[j - 1],ctx->PMGN)) > 0.0)
                be += tmp * tmp1;
            if (j < n - 1 && (tmp1 = gdd_adj(ctx, i,idx[j + 1],ctx->PMGN)) > 0.0)
                be += tmp * tmp1;
        }
    }
    return(be);
}

/* ------------------------------------------------------------------------ */
/*  be_rcheck(n,idx)                                                        */
/*                                                                          */
/*  Calculate BE measure for rows with nodes in idx[i], i=0,...,n-1.        */

double be_rcheck(TDAContext *ctx, int n,int *idx)                          
{
    register int i,j,k;
    double be,tmp,tmp1;

    be = 0.0;
    for (i = 0; i < ctx->GD_NP; ++i) {
        for (j = 0; j < n; ++j) {
            k = idx[j];
            tmp = gdd_adj(ctx, i,k,ctx->PMGN);
            if (tmp <= 0.0)
                continue;

            if (j > 0 && (tmp1 = gdd_adj(ctx, idx[j - 1],i,ctx->PMGN)) > 0.0)
                be += tmp * tmp1;
            if (j < n - 1 && (tmp1 = gdd_adj(ctx, idx[j + 1],i,ctx->PMGN)) > 0.0)
                be += tmp * tmp1;
        }
    }
    return(be);
}

/* ------------------------------------------------------------------------ */
/*  acl         Clustering with assignment.                                 */
/*                                                                          */
/*              Requires an undirected valued graph.                        */
/*                                                                          */
/*              acl(                                                        */
/*                  cn=...,     cn=k1,k2,...  (1 <= k1,k2,.. <= n/2)        */
/*                  gn = ...,   graph number, def. 1                        */
/*                  df = ...,   additional output file                      */  
/*                  nfmt = ..., integer print format, def. 4                */
/*                  fmt = ...,  print format, def. 8.2                      */
/*              ) = fname;                                                  */
/*                                                                          */
/*  fname will contain indicators for the assignment of objects to          */
/*  prototypes.                                                             */
/*                                                                          */  
/*  If an additional output file is requested with the df parameter, the    */
/*  file will contain weights for the assignment of objects to prototypes.  */
/*  Only if all distances are different from zero.                          */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int acl(TDAContext *ctx)
{
    register int i,j,k; 
    int err,np,nk,nc,ncmax,nn,first,kk,pflag;
    double r;
    double c,mev;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Clustering with assignment. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto ACLFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto ACLFin;
    if (gdd_tcheck(ctx, 0,1,2))          /* need undirected valued graph */
        goto ACLFin;

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 8,2);

    np = ctx->GD_NP;                     /* number of objects */

    nk = ctx->PMNCN;
    if (nk == 0) {
        printf1(ctx, "Need value for cn parameter.\n");
        goto ACLFin;
    }
    ncmax = 0;
    for (k = 0; k < nk; ++k) {
        nc = ctx->PMCN[k];
        if (nc < 1 || nc > np / 2) {
            printf1(ctx, "\nError: range of cn parameter is %d to %d\n",1,np/2);
            goto ACLFin;
        }
        ncmax = imax(ctx, ncmax,nc);
    }
    newline(ctx);

    if (alloc_acxf(ctx, np * np + 1))    /* distance matrix */
        goto ACLFin;
    if (alloc_acm(ctx, ncmax))           /* prototype numbers */
        goto ACLFin;
    if (alloc_acn(ctx, nk * ncmax))      /* final prototype numbers */
        goto ACLFin;
    if (alloc_acr(ctx, np))              /* cluster membership indicators */
        goto ACLFin;
    if (alloc_acns(ctx, np * nk))        /* final cluster membership indicators */
        goto ACLFin;
    if (alloc_acyf(ctx, nk))             /* final assignment value */
        goto ACLFin;

    mev = ctx->GD_EVMax[ctx->PMGN-1];

    pflag = 1;                      /* set to 0 if zero distances */
    k = 0;
    for (i = 0; i < np; ++i) {
        for (j = 0; j < np; ++j) {
            if (j == i)
                ctx->AcXF[k++] = 0.0;
            else {
                ctx->AcXF[k] = (float)(gdd_adj(ctx, i,j,ctx->PMGN));
                if ((double)(ctx->AcXF[k]) == 0.0)
                    pflag = 0;
                k++;
            }
        }
    }
    /*******
    for (i = 0; i < np; ++i) {
        for (j = 0; j < np; ++j) {
            tda_out("%2.0f ",AcXF[i * np + j]);
        }
        newline(ctx);
    }
    ***/

    for (kk = 0; kk < nk; ++kk) {       /* for all cluster sizes */
        nc = ctx->PMCN[kk];
        ctx->AcYF[kk] = (float)(np * mev);

        first = 1;
        while (comb_nm(ctx, np,nc,ctx->AcM,first)) {
            for (i = 0; i < np; ++i)
                ctx->AcR[i] = -1;
            for (j = 0; j < nc; ++j)
                ctx->AcR[ctx->AcM[j]] = j;

            c = (double)((float)(acl1(ctx, np,nc,ctx->AcXF,ctx->AcM,ctx->AcR,(float)(mev),ctx->AcYF[kk])));

            if (first || (double)(c) < (double)(ctx->AcYF[kk])) {
                ctx->AcYF[kk] = (float)(c);
                j = kk * ncmax;
                for (k = 0; k < nc; ++k)  
                    ctx->AcN[j++] = ctx->AcM[k];
                j = kk * np;
                for (k = 0; k < np; ++k)
                    ctx->AcNS[j++] = (short)(ctx->AcR[k]);
            }
            first = 0;
        }
    }
    for (kk = 0; kk < nk; ++kk) {
        nc = ctx->PMCN[kk];
        printf1(ctx, "nc=%3d  cmin=%12.2f  prototypes: ",nc,(double)(ctx->AcYF[kk]));
        j = kk * ncmax;
        for (k = 0; k < nc; ++k)
            printf1(ctx, "%3d ",ctx->AcN[j++] + 1);
        newline(ctx);
    }
    newline(ctx);

    /* print to output file */

    for (i = 0; i < np; ++i) {
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
        for (k = 0; k < nk; ++k)  
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcNS[k * np + i] + 1);
        fprintf(ctx->PMFd,"\n");
    }
    printf1(ctx, "%d records written to: %s\n",np,ctx->PMFdName);

    if (ctx->PMF1Def) {          /* additional output file */
        if (pflag == 0) {
            printf1(ctx, "df option ignored.\n");
            goto ACLFin;
        }
        if (alloc_actmp(ctx, ncmax))             
            goto ACLFin;

        nn = 0;
        for (kk = 0; kk < nk; ++kk) {
            nc = ctx->PMCN[kk];
            fprintf(ctx->PMF1d,"Assignment weights for nc = %d [%d",nc,ctx->AcN[kk * ncmax] + 1);
            for (k = 1; k < nc; ++k)
                fprintf(ctx->PMF1d,",%d",ctx->AcN[kk * ncmax + k] + 1);
            fprintf(ctx->PMF1d,"]\n");
            nn++;    

            for (i = 0; i < np; ++i) {
                r = 0.0;    
                for (k = 0; k < nc; ++k) {
                    j = ctx->AcN[kk * ncmax + k];
                    if (i == j) {
                        for (j = 0; j < nc; ++j)
                            ctx->AcTmp[j] = 0.0;
                        ctx->AcTmp[k] = 1.0;
                        r = 0.0;    
                        break;
                    }
                    ctx->AcTmp[k] = 1.0 / (double)(ctx->AcXF[i * np + j]);
                    r += ctx->AcTmp[k]; 
                }
                if (r > 0.0) {
                    for (k = 0; k < nc; ++k)
                        ctx->AcTmp[k] /= r;
                }
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,i + 1);
                for (k = 0; k < nc; ++k)
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcTmp[k]);
                fprintf(ctx->PMF1d,"\n");
                nn++;
            }
            fprintf(ctx->PMF1d,"\n");
            nn++;
        }
        printf1(ctx, "%d records written to: %s\n",nn,ctx->PMF1dName);
    }
    err = 0;

ACLFin:       
    p_clean(ctx);
    return(err);
}


/* ------------------------------------------------------------------------ */
/*  acl1(n,m,d,p,q,mev)                                                     */
/*                                                                          */
/*  n is number of objects, m number of clusters, d[] distance matrix,      */
/*  p[] indices of prototypes. q[] used to index cluster membership.        */
/*  mev is max edge value. Set indices in q[] and return min. assignment    */
/*  value.                                                                  */

float acl1(TDAContext *ctx, int n,int m,float *d,int *p,int *q,float mev,float cmin)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j,k = 0;
    double dm,t,c;
         
    c = 0.0;
    for (i = 0; i < n; ++i) {
        if (q[i] >= 0)
            continue;

        dm = (double)(mev);
        for (j = 0; j < m; ++j) {
            t = (double)(d[i * n + p[j]]);
            if (t < dm) {
                dm = t;
                k = j;
            }
        }
        q[i] = k;
        c += dm;
        if ((double)(c) > (double)(cmin))  
            break;
    }
    return ((float)(c));
}

/* ------------------------------------------------------------------------ */
/*  udcl(n,nc,d,p,l,u,c,pc,qa,fp)                                           */
/*                                                                          */
/*  Find best partition of an ordering into nc contiguous groups.           */
/*  Dynamic programming algorithm as described by C.J. Alpert, A.B. Kahng,  */
/*  Splitting an Ordering into a Partition to Minimize Diameter.            */
/*  J. of Classification 14 (1997), p. 64.                                  */
/*                                                                          */
/*  Note: this function always minimizes the maximum of the diameters of    */
/*  the partitions.                                                         */
/*                                                                          */
/*  n       number of data points in distance matrix d(i,j)                 */
/*  nc      number of clusters                                              */
/*  d[]     distance matrix (i,j=0,n-1)                                     */
/*  p[]     permutation of objects (i=0,n-1)                                */
/*  l       min size of partition                                           */
/*  u       max size of partition                                           */
/*  c[]     cost matrix (i,j=0,n-1)                                         */
/*  pc[]    cost matrix (i=0,nc-1,j=0,n-1), DBLMAX marks infeasible states  */
/*  qa[]    pointer for partitions (i=0,nc-1,j=0,n-1)                       */
/*  fp[]    final partition (i=0,nc-1) [starting indices]                   */
/*                                                                          */
/*  The final cost is in pc[nc-1,n-1].                                      */
/*                                                                          */
/*  Return 0 if OK, -1 if no feasible partition exists.                     */

int udcl(TDAContext *ctx, int n,int nc,float *d,int *p,int l,int u,float *c,double *pc,int *qa,int *fp)
{
    register int j,k,m;
    int km1,mlo,mhi,a;
    double fb,f;

    udcl1(ctx, n,d,p,c,u);           /* initial cost */

    m = n * nc;
    for (j = 0; j < m; ++j)
        pc[j] = ctx->DBLMAX;

    for (j = l - 1; j < u && j < n; ++j)
        pc[j] = (double)c[j];

    for (k = 1; k < nc; ++k) {
        km1 = k - 1;
        for (j = 0; j < n; ++j) {
            fb = ctx->DBLMAX;
            mlo = imax(ctx, j - u,k * l - 1);
            mhi = imin(ctx, j - l,k * u - 1);
            for (m = mlo; m <= mhi; ++m) {
                if (m < 0)
                    continue;
                if (pc[km1 * n + m] >= ctx->DBLMAX)
                    continue;
                f = dmax(ctx, pc[km1 * n + m],(double)c[(m + 1) * n + j]);
                if (fb > f) {
                    fb = f;
                    pc[k * n + j] = f;
                    qa[k * n + j] = m;
                }
            }
        }
    }
    if (pc[(nc - 1) * n + n - 1] >= ctx->DBLMAX)
        return(-1);

    m = n - 1;
    for (k = nc - 1; k >= 1; --k) {
        a = qa[k * n + m];
        fp[k] = a + 1;
        m = a;
    }
    fp[0] = 0;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  udcl1(n,d,p,c,u)            calculate cluster costs c[i,j]              */
/*                                                                          */
/*  c[i,j] is the diameter of the segment p[i],...,p[j], computed for       */
/*  segment lengths up to u.                                                */
/*                                                                          */
/*  n       number of data points in distance matrix d(i,j)                 */
/*  d[]     distance matrix (i,j = 0,...,n-1)                               */
/*  p[]     permutation of objects (i = 0,...,n-1)                          */
/*  c[]     cluster costs (i,j = 0,...,n-1)                                 */
/*  u       max cluster size                                                */

void udcl1(TDAContext *ctx, int n,float *d,int *p,float *c,int u)
{
    register int i,j,l;

    for (i = 0; i < n * n; ++i)
        c[i] = 0.0;

    for (l = 1; l < u; ++l) {
        for (i = 0; i + l < n; ++i) {
            j = i + l;
            c[i * n + j] =
            (float)dmax(ctx, dmax(ctx, (double)c[i * n + j - 1],(double)c[(i + 1) * n + j]),
                      (double)d[p[i] * n + p[j]]);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  ucl         Clustering of an ordering.                                  */
/*              Splits an ordering of the nodes into contiguous clusters    */
/*              so that the maximal cluster diameter is minimized (Alpert   */
/*              and Kahng, J. of Classification 14, 1997, p. 64).           */
/*                                                                          */
/*              Requires an undirected valued graph.                        */
/*                                                                          */
/*              ucl(                                                        */
/*                  nc=...,     number of clusters, def. 2                  */
/*                  min=...,    min cluster size, def. 1                    */
/*                  max=...,    max cluster size, def. n                    */
/*                  cn=...,     ordering of nodes (internal indices),       */
/*                              def. 1,2,...,n                              */
/*                  gn=...,     graph number, def. 1                        */
/*                  nfmt=...,   integer print format, def. 4                */
/*                  fmt=...,    print format, def. 10.4                     */
/*              ) = fname;                                                  */
/*                                                                          */
/*  fname contains one record per node: position in the ordering, internal  */
/*  node index, cluster number, external node number.                       */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int ucl(TDAContext *ctx)
{
    register int i,j,k; 
    int err,np,nc,l,u,cl,nrec;
    double t;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Clustering of an ordering. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto UCLFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto UCLFin;

    if (ctx->PMFDef == 0) {
        printf1(ctx, "Error: need an output file.\n");
        goto UCLFin;
    }

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    np = ctx->GD_NP;                     /* number of objects */

    nc = ctx->PMNC;
    if (nc < 2)
        nc = 2;
    if (nc > np) {
        printf1(ctx, "Error: number of clusters must be at most %d.\n",np);
        goto UCLFin;
    }
    l = ctx->PMMin;
    if (l < 1)
        l = 1;
    u = ctx->PMMax;
    if (u < 1)
        u = np;
    if (u < l)
        u = l;
    if (nc * l > np || nc * u < np) {
        printf1(ctx, "Error: no partition of %d nodes into %d clusters of size %d to %d.\n",np,nc,l,u);
        goto UCLFin;
    }
    printf1(ctx, "Number of clusters: %d\n",nc);
    printf1(ctx, "Cluster sizes: %d to %d\n",l,u);

    if (alloc_acxf(ctx, np * np + 1))    /* distance matrix */
        goto UCLFin;
    if (alloc_acyf(ctx, np * np + 1))    /* cost matrix */
        goto UCLFin;
    if (alloc_acu(ctx, nc * np + 1))     /* cost matrix for partitions */
        goto UCLFin;
    if (alloc_ack(ctx, np))              /* ordering */
        goto UCLFin;
    if (alloc_acn(ctx, nc * np + 1))     /* pointer to partitions */
        goto UCLFin;
    if (alloc_aci(ctx, nc + 1))          /* final partition */
        goto UCLFin;

    if (ctx->PMNCN > 0) {                /* ordering given with cn */
        if (ctx->PMNCN != np) {
            printf1(ctx, "Error: cn parameter needs %d node indices.\n",np);
            goto UCLFin;
        }
        if (alloc_acj(ctx, np))
            goto UCLFin;
        for (i = 0; i < np; ++i)
            ctx->AcJ[i] = 0;
        for (i = 0; i < np; ++i) {
            j = ctx->PMCN[i] - 1;
            if (j < 0 || j >= np || ctx->AcJ[j]) {
                printf1(ctx, "Error: cn parameter is not a permutation of 1,...,%d.\n",np);
                goto UCLFin;
            }
            ctx->AcJ[j] = 1;
            ctx->AcK[i] = j;
        }
    }
    else {
        for (i = 0; i < np; ++i)
            ctx->AcK[i] = i;
    }

    k = 0;
    for (i = 0; i < np; ++i) {
        for (j = 0; j < np; ++j) {
            if (j == i)
                ctx->AcXF[k++] = 0.0;
            else {
                t = gdd_adj(ctx, i,j,ctx->PMGN);
                if (t < 0.0) {
                    printf1(ctx, "Error: found negative or missing edge value.\n");
                    goto UCLFin;
                }
                ctx->AcXF[k++] = (float)t;
            }
        }
    }

    if (udcl(ctx, np,nc,ctx->AcXF,ctx->AcK,l,u,ctx->AcYF,ctx->AcU,ctx->AcN,ctx->AcI)) {
        printf1(ctx, "Error: no feasible partition found.\n");
        goto UCLFin;
    }

    printf1(ctx, "Maximal cluster diameter: %lg\n",ctx->AcU[(nc - 1) * np + np - 1]);

    nrec = 0;
    cl = 0;
    for (j = 0; j < np; ++j) {
        if (cl < nc - 1 && j == ctx->AcI[cl + 1])
            cl++;
        i = ctx->AcK[j];
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j + 1);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,cl + 1);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, i));
        fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
        tda_export_cell(ctx, "ucl.part", (double)(j + 1));
        tda_export_cell(ctx, "ucl.part", (double)(i + 1));
        tda_export_cell(ctx, "ucl.part", (double)(cl + 1));
        tda_export_cell(ctx, "ucl.part", (double)gdd_node(ctx, i));
        tda_export_endrow(ctx, "ucl.part");
#endif
        nrec++;
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

UCLFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  hclsp       Partitioning the leafs of a tree.                           */
/*                                                                          */
/*              hclsp(                                                      */
/*                  nlev=...,   number of levels, def. 0                    */
/*                  cn=...,     sequence of node numbers                    */
/*                  nfmt = ..., integer print format, def. 4                */
/*              ) = output_file;                                            */
/*                                                                          */
/*  Required: directed graph defined with gdd option 1 (edge list). The     */
/*  graph number is always 1. The graph must be a tree with a single node   */
/*  that has outdegree 0. This is used as the root node. In particular,     */
/*  the tree can be created with the df option in the hcls command.         */
/*                                                                          */
/*  There are two options (that can be used simultaneously).                */
/*                                                                          */
/*  (1) If nlev > 0, the command determines the root node, that is, the     */
/*      single node with outdegree 0 and reconstructs the tree for nlev     */  
/*      levels, beginning from the root node. The truncated tree is written */
/*      to the output file and shows, in addition to the node numbers, the  */
/*      number of leafs that have a directed path to the node.              */
/*                                                                          */
/*  (2) The second option requires the specification of a sequence of       */
/*      nodes (external node numbers) with the parameter                    */
/*                                                                          */  
/*          cn = n1,n2,...,                                                 */
/*                                                                          */
/*      The output file will then contain one record for each leaf of the   */
/*      tree that has a path to one of the nodes n1,n2,... There will be    */
/*      two columns containing one of the node numbers from n1,n2,..,       */
/*      followed by the node number of the leaf.                            */
/*                                                                          */
/*      Note that is is not required that the nodes n1,n2,... define a      */
/*      complete partition of the leafs. If a leaf has two or more nodes    */  
/*      from the set n1,n2,... as followers, the assignment is arbitrary.   */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int hclsp(TDAContext *ctx)
{
    register int i; 
    int err,gn,n,np,ne,rn,nl,level,len;
    char buf[60];

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Partitioning the leafs of a tree. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 5,1,1))           /* get parameters */
        goto HCLSPFin;

    gn = 1;                             /* always graph number 1 */
    if (gdd_check(ctx, gn,0))        
        goto HCLSPFin;
    if (gdd_tcheck(ctx, 1,2,0))              /* need undirected valued graph */
        goto HCLSPFin;

    if (ctx->PMFmtF == 0)                    /* default print format */
        pmfmt(ctx, 8,2);

    np = ctx->GD_NP;                         /* number of nodes */
    ne = ctx->GD_NE[0];                      /* number of edges */
    printf1(ctx, "Number of nodes: %d\n",np);
    printf1(ctx, "Number of edges: %d\n",ne);
    if (ne != np - 1) {
        printf1(ctx, "Error: graph is not a tree.\n");
        goto HCLSPFin;
    }
    nl = 0;
    rn = -1;
    for (i = 0; i < np; ++i) {
        if (ctx->GD_BPN[i] == 0)
            nl++;
        if (ctx->GD_FPN[i] == 0) {
            if (rn < 0)
                rn = i;
            else {
                rn = -1;
                break;
            }
        }
    }
    if (nl == 0 || rn < 0) {
        printf1(ctx, "Error: graph is not a tree or has no single root.\n");
        goto HCLSPFin;
    }
    printf1(ctx, "Number of leafs: %d\n",nl);
    printf1(ctx, "Root node number: %d\n",ctx->GD_ND[rn]);

    len = 1;
    for (i = 0; i < np; ++i) {
        snprintf(buf,sizeof(buf),"(%d,%d) ",ctx->GD_ND[i],nl);
        len = (int)(imax(ctx, len,(int)(strlen(buf))));
    }
    if (alloc_acn(ctx, np + 1))
        goto HCLSPFin;

    for (i = 0; i < np; ++i) {
        if (ctx->GD_BPN[i] == 0)
            hclsp_visit(ctx, i,gn);
    }
    if (ctx->PMNLEV > 0) {
        fprintf(ctx->PMFd,"Hierarchical structure of internal nodes (number of leafs in brackets).\n\n");
        level = 0;
        hclsp_prn(ctx, rn,gn,ctx->PMNLEV,level,len);
        fprintf(ctx->PMFd,"\n");   
    }
    if (ctx->PMNCN > 0) {
        if (alloc_acd(ctx, np + 1))
            goto HCLSPFin;

        for (i = 0; i < ctx->PMNCN; ++i) {
            n = gdd_fndni(ctx, ctx->PMCN[i]);
            if (n < 0) {
                printf1(ctx, "Error: node number %d not defined.\n",ctx->PMCN[i]);
                goto HCLSPFin;
            }
            hclsp_prn1(ctx, n,gn,ctx->PMCN[i]);
        }
        n = 0; 
        for (i = 0; i < np; ++i) {
            if (ctx->AcD[i] && ctx->GD_BPN[i] == 0)
                n++;
        }
        printf1(ctx, "Number of leafs in the output file: %d\n",n);
    }
    err = 0;

HCLSPFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  hclsp_visit(i,gn)                                                       */

void hclsp_visit(TDAContext *ctx, int i,int gn)
{
    register int j,l,k;
    int n;

    ctx->AcN[i] += 1;               
    n = ctx->GD_FPN[i];
    if (n <= 0)
        return;          

    for (l = 0; l < n; ++l) {
        j = ctx->GD_FPI[i][l];                    /* edge from k to j */
        k = ctx->GD_FPK[i][l];
        if (k >= 0 && gdd_ev(ctx, k,gn) >= 0.0)          
            hclsp_visit(ctx, j,gn);
    }
}

/* ------------------------------------------------------------------------ */
/*  hclsp_prn(i,gn,nlev,level,mlen)                                         */

void hclsp_prn(TDAContext *ctx, int i,int gn,int nlev,int level,int mlen)
{
    register int j,l,k,ll;
    int n,nline,len;
    char buf[60];

    if (level >= nlev)  
        return;
           
    if (level != 0)  
        fprintf(ctx->PMFd,"<- ");
    snprintf(buf,sizeof(buf),"(%d,%d) ",ctx->GD_ND[i],ctx->AcN[i]);
    len = (int)(strlen(buf));
    fprintf(ctx->PMFd,"%s",buf);
    fprnchar(ctx, ctx->PMFd,' ',mlen - len,0);

    if (level + 1 == nlev || (n = ctx->GD_BPN[i]) <= 0) {
        fprintf(ctx->PMFd,"\n");   
        return;
    }
    nline = 0;
    for (l = 0; l < n; ++l) {
        j = ctx->GD_BPI[i][l];                    /* edge from k to j */
        k = ctx->GD_BPK[i][l];
        if (k >= 0 && gdd_ev(ctx, k,gn) >= 0.0) {        
            if (nline) {
                fprnchar(ctx, ctx->PMFd,' ',mlen,0);
                for (ll = 0; ll < level; ++ll) 
                    fprnchar(ctx, ctx->PMFd,' ',mlen + 3,0);
            }
            hclsp_prn(ctx, j,gn,nlev,level + 1,mlen);
            nline = 1;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  hclsp_prn1(i,gn)                                                       */

void hclsp_prn1(TDAContext *ctx, int i,int gn,int m)
{
    register int j,l,k;
    int n;

    if (ctx->AcD[i])
        return;

    ctx->AcD[i] = 1;
    if (ctx->GD_BPN[i] == 0) {       /* print only for leafs */
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,m);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GD_ND[i]);
        fprintf(ctx->PMFd,"\n");
    }                
    if ((n = ctx->GD_BPN[i]) <= 0)
        return;
                               
    for (l = 0; l < n; ++l) {
        j = ctx->GD_BPI[i][l];                    /* edge from k to j */
        if (ctx->AcD[j] == 0) {
            k = ctx->GD_BPK[i][l];
            if (k >= 0 && gdd_ev(ctx, k,gn) >= 0.0)          
                hclsp_prn1(ctx, j,gn,m);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  scla        Find separable clusters.                                    */
/*              Requires an undirected valued graph.                        */
/*                                                                          */
/*              scla(                                                       */
/*                  gn = ...,   graph number, def. 1                        */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int scla(TDAContext *ctx)
{
    register int i,j;
    int err,fin,nnsp,n;


    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Find separable clusters. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto SCLAFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto SCLAFin;
    if (gdd_tcheck(ctx, 0,1,2))          /* need undirected valued graph */
        goto SCLAFin;

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    printf1(ctx, "Number of nodes: %d\n",ctx->GD_NP);

    /* For i=0,GD_NP-1: AcN[i] = 1,2,3,... cluster number, or 0 if node i
       cannot be used for an s-cluster. */

    if (alloc_acn(ctx, ctx->GD_NP))                 
        goto SCLAFin;
    if (alloc_acm(ctx, ctx->GD_NP))                 
        goto SCLAFin;
    if (alloc_acx(ctx, ctx->GD_NP))                 
        goto SCLAFin;

    for (i = 0; i < ctx->GD_NP; ++i)
        ctx->AcN[i] = -1;

    nnsp = 0;       /* number of point that cannot be userd for s-clusters */
    fin = 0;
    while (fin == 0) {
        fin = 1;
        for (i = 0; i < ctx->GD_NP; ++i) {
            if (ctx->AcN[i] >= 0)
                continue;
            if (scla_check(ctx, i,ctx->GD_NP,ctx->AcN,ctx->AcM,ctx->AcX)) {
                ctx->AcN[i] = 0;
                nnsp++;
                fin = 0;
                break;
            }
        }
    }
    printf1(ctx, "\n%d isolated points",nnsp);
    if (nnsp == 0)
        printf1(ctx, ".\n");
    else {
        j = 0;
        for (i = 0; i < ctx->GD_NP; ++i) {
            if (ctx->AcN[i] == 0) {
                if (j && j < nnsp)
                    printf1(ctx, ",");
                else
                    printf1(ctx, ": ");
                printf1(ctx, "%d",i + 1);
                j++;
            }
        }
        printf1(ctx, "\n");
    }
    n = 0;
    for (i = 0; i < ctx->GD_NP; ++i) {
        if (ctx->AcN[i] < 0)  
            ctx->AcM[n++] = i;
    }
    printf1(ctx, "Remaining number of points: %d\n\n",n);
    if (n < 4) {
        err = 0;
        goto SCLAFin;
    }
    if (alloc_ack(ctx, ctx->GD_NP))                 
        goto SCLAFin;
    scla_find(ctx, ctx->GD_NP,ctx->AcN,n,ctx->AcM,ctx->AcK);
    err = 0;

SCLAFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  scla_check()                                                            */
/*                                                                          */
/*  Return  1 if i cannot be used for an s-cluster                          */
/*          0 otherwise                                                     */
        
int scla_check(TDAContext *ctx, int i,int n,int *acn,int *p,double *dp)
{
    register int j,k,jx;
    int np,nn,fin;  
    double x,dmin,dm;

    p[0] = i;
    dp[0] = 0.0;
    np = 1;
    acn[i] = -2;

    /* first add nearest node to p[] */

    jx = -1;
    dmin = ctx->DBLMAX;
    for (j = 0; j < n; ++j) {
        if (acn[j] == -1) {
            x = gdd_adj(ctx, i,j,ctx->PMGN);            
            if (x < dmin) {
                dmin = x;
                jx = j;
            }
        }
    }
    if (jx < 0)  
        return(1);
    p[np++] = jx;
    acn[jx] = -2;
    dp[0] = dp[1] = dmin;
      
    fin = 0; 
    while (fin == 0) {
        fin = 1;
        jx = -1;
        for (j = 0; j < n; ++j) {
            if (acn[j] == -1) {
                jx = -1;
                for (k = 0; k < np; ++k) {
                    if (gdd_adj(ctx, p[k],j,ctx->PMGN) <= dp[k]) {
                        jx = j;
                        break;
                    }
                }
                if (jx >= 0) {
                    dm = 0.0;
                    for (k = 0; k < np; ++k) {
                        x = gdd_adj(ctx, p[k],jx,ctx->PMGN);
                        if (dp[k] < x)
                            dp[k] = x;
                        dm = dmax(ctx, dm,x);
                    }   
                    dp[np] = dm;       
                    p[np++] = jx;
                    fin = 0;
                    acn[jx] = -2;
                }
            }
        }
    }
    nn = 0;
    for (j = 0; j < n; ++j) {
        if (acn[j] == -1)
            nn++;
        else if (acn[j] == -2)
            acn[j] = -1;
    }
    if (nn > 0)
        return(0);

    return(1);
}

/* ------------------------------------------------------------------------ */
/*  scla_find()  Find largest s-clusters.                                   */
/*                                                                          */
/*  use all points with acn[0,...,n-1] = -1.                                */
/*  Return -1 if error                                                      */
/*          number of clusters found.                                       */
        
int scla_find(TDAContext *ctx, int n,int *acn,int np,int *p,int *com)
{
    register int i,j,k,l;
    int first,nc,nu,nf;
    double dm;

    for (nc = 2; nc <= np - 2; ++nc) {
        printf1(ctx, "Size %d\n",nc);
        nf = 0;
        first = 1;
        while (comb_nm(ctx, np,nc,com,first)) {
            first = 0;
            for (i = 0; i < nc; ++i) {
                j = p[com[i]];
                acn[j] = 1;
            }
            nu = 0;
            for (i = 0; i < nc; ++i) {
                j = p[com[i]];
                dm = 0.0;
                for (k = 0; k < nc; ++k) {
                    if (k != i)
                        dm = dmax(ctx, dm,gdd_adj(ctx, j,p[com[k]],ctx->PMGN));            
                }
                for (l = 0; l < n; ++l) {
                    if (acn[l] < 0 && gdd_adj(ctx, l,j,ctx->PMGN) <= dm) {
                        nu = 1;
                        break;
                    }
                }
                if (nu)
                    break;
            }
            for (i = 0; i < nc; ++i) {
                j = p[com[i]];
                acn[j] = -1;
            }
            if (nu==0) {
                printf1(ctx, "%4d -: ",++nf);
                for (i = 0; i < nc; ++i)
                    printf1(ctx, "%d ",p[com[i]] + 1);
                newline(ctx);
            }
        }
    }
    return(0);
}


/* ------------------------------------------------------------------------ */
/*  clp         Partitional clustering.                                     */
/*                                                                          */
/*              clp(                                                        */
/*                  meth=...,   method (1,2,3), def. 2                      */
/*                  gn=...,     graph number, def. 1                        */
/*                  nc=...,     number of clusters, def. 2                  */
/*                  ns=...,     number of random restarts, def. 1           */
/*                  mxit=...,   max number of iterations, def. 15           */
/*                  fmt=...,    print format for values, def. 10.4          */
/*              ) = fname;                                                  */
/*                                                                          */
/*                                                                          */
/*              Requires an undirected valued graph. Uses the exchange      */
/*              method described in H. Spaeth, Cluster-Formation und        */
/*              -Analyse, 1983, p. 143.                                     */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int clp(TDAContext *ctx)
{
    register int i,j,k;
    int err,n,nn,r,kmin;    
    double t,cmin;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Partitional clustering. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto CLPFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto CLPFin;
    if (gdd_tcheck(ctx, 0,0,0))
        goto CLPFin;

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    if (ctx->PMMETH < 1 || ctx->PMMETH > 4)
        ctx->PMMETH = 2;

    if (ctx->PMNC < 2)
        ctx->PMNC = 2;
    if (ctx->PMNC >= ctx->GD_NP) {
        printf1(ctx, "Number of clusters must be less than %d\n",ctx->GD_NP);
        goto CLPFin;
    }        

    if (ctx->PMNS < 1)
        ctx->PMNS = 1;

    if (ctx->MxItFlg == 0)
        ctx->MxIter = 15;

    printf1(ctx, "\nMethod: %d\n",ctx->PMMETH);

    if (ctx->PMMETH == 4) {
        printf1(ctx, "Find %d points with maximal min-distance.\n",ctx->PMNC);
        err = clp_p(ctx);
        goto CLPFin;
    }

    printf1(ctx, "Number of clusters: %d\n",ctx->PMNC);
    printf1(ctx, "Number of random restarts: %d\n",ctx->PMNS);
    printf1(ctx, "Maximal number of iterations: %d\n",ctx->MxIter);

    n = ctx->GD_NP;
    nn = n * (n - 1) / 2;

    if (alloc_act(ctx, nn + 1))      
        goto CLPFin;
    if (alloc_acu(ctx, ctx->PMNC + 1))      
        goto CLPFin;
    if (alloc_acv(ctx, ctx->PMNC + 1))      
        goto CLPFin;
    if (alloc_acn(ctx, ctx->PMNC + 1))      
        goto CLPFin;
    if (alloc_ack(ctx, n + 1))      
        goto CLPFin;
    if (alloc_acl(ctx, n + 1))      
        goto CLPFin;
   
    k = 0;
    for (i = 1; i < n; ++i) {
        for (j = 0; j < i; ++j)  
            ctx->AcT[++k] = gdd_adj(ctx, i,j,ctx->PMGN);
    }
    newline(ctx);
                
    cmin = ctx->DBLMAX;
    kmin = -1;

    for (k = 1; k <= ctx->PMNS; ++k) {           /* repetitions */
        for (i = 1; i <= n; ++i) {          /* create initial partition */
            t = random1(ctx);
            t *= (double)ctx->PMNC;
            j = (int)(t + 1.0);
            if (j < 1)
                j = 1;
            else if (j > ctx->PMNC)
                j = ctx->PMNC;
            ctx->AcK[i] = j;
        }
        r = clp_f(ctx, ctx->PMMETH,n,ctx->PMNC,ctx->MxIter,ctx->AcT,ctx->AcK,ctx->AcN,ctx->AcU,ctx->AcV);

        if (r >= 0) {
            t = 0.0;
            for (j = 1; j <= ctx->PMNC; ++j)
                t += ctx->AcU[j];
            if (t < cmin) {
                kmin = k;
                cmin = t;
                for (i = 1; i <= n; ++i)
                    ctx->AcL[i] = ctx->AcK[i];
            }
        }
        else
            t = -1.0;

        printf1(ctx, "%4d %2d ",k,r);
        rt_printf1_d(ctx, ctx->PMFmtS,t);
        for (j = 1; j <= ctx->PMNC; ++j) {
            printf1(ctx, "%3d ",ctx->AcN[j]);
            rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcU[j]);
        }
        newline(ctx);
    }
    printf1(ctx, "\nBest solution: %d\n",kmin);

    if (kmin > 0) {
        for (i = 1; i <= n; ++i)  
            fprintf(ctx->PMFd,"%4d %3d\n",i,ctx->AcL[i]);
        printf1(ctx, "%d records written to: %s\n",n,ctx->PMFdName);
    }
    err = 0;

CLPFin:       
    p_clean(ctx);
    return(err);
}
   

/* ------------------------------------------------------------------------ */
/*  clp_f       Partitional clustering.                                     */
/*                                                                          */
/*  meth = 1 or 2 or 3                                                      */
/*  m = number of objects. m12 = m * (m-1) / 2                              */
/*  t[k], k=1,...,m12 should contain the distances.                         */
/*  distance(i,h) corresponds to t[k] with k = (i-1)(i-2)/2 + h             */
/*                                                                          */
/*  n = number of clusters                                                  */  
/*  z[i] (i = 1,...,m) contains the initial partition. z[i] = j if          */
/*  object i belongs to cluster j.                                          */
/*  mj[j] is number of elements in cluster j                                */
/*  e[j] (j=1,...,n)                                                        */
/*  c[j] (j=1,...,n)                                                        */
/*                                                                          */
/*  Return  0 if successful                                                 */
/*          1 if maximal number of iterations                               */
/*         -1 if a cluster contains less than minimal number of objects     */


int clp_f(TDAContext *ctx, int meth,int m,int n,int itmax, double *t,int *z,int *mj,double *e,double *c)
{
    register int i,j,h,l,k;
    int it,mo,p,q = 0,im,ki,kh,is;
    double d,f,v,v1,u,u1,w = 0.0,w1 = 0.0,ej,ep = 0.0,eq,bj,cj;
    double r = 0.999;

    it = 0;
    d = 0.0;
    for (j = 1; j <= n; ++j) {
        mj[j] = 0; 
        e[j] = 0.0;
    }

    mo = 1;                         /* minimal cluster size */
    if (meth == 3)
        mo = 2;

    for (i = 1; i <= m; ++i)        /* determine cluster sizes */
        mj[z[i]] += 1;
        
    for (j = 1; j <= n; ++j) {      /* check min cluster size */
        if (mj[j] < mo)
            return(-1);
    }

    k = 0;
    for (i = 2; i <= m; ++i) {
        p = z[i];
        l = i - 1;
        for (h = 1; h <= l; ++h) {
            q = z[h];
            k++;
            if (q == p)
                e[q] += t[k];
        }
    }
    for (j = 1; j <= n; ++j) {
        f = e[j];
        if (meth != 1) {
            v = (double)mj[j];
            if (meth == 2)
                f /= v;
            else if (meth == 3)
                f /= (v * (v - 1.0));
            e[j] = f;
        }
        d += f;
    }
    if (n <= 1)
        return(0);

    i = is = 0;

L9:
    is++;
    if (is > m)
        return(0);

L10:
    i++;
    if (i > m) {
        it++;
        if (it > itmax)
            return(1);
        i = 1;
    }
  
    p = z[i];
    l = mj[p];
    if (l <= mo)
        goto L9;

    v = (double)l;
    v1 = v - 1.0;
    im = ((i - 1) * (i - 2)) / 2;
    for (j = 1; j <= n; ++j)
        c[j] = 0.0;

    for (h = 1; h <= m; ++h) {
        if (h == i)
            continue;
        j = z[h];
        kh = h;
        ki = im;
        if (i <= h) {
            ki = i;
            kh = ((h - 1) * (h - 2)) / 2;
        }
        k = kh + ki;
        c[j] += t[k];
    }
    eq = ctx->DBLMAX;

    for (j = 1; j <= n; ++j) {
        bj = e[j];
        cj = c[j];
        u = (double)mj[j];
        u1 = u + 1.0;
        if (j == p) {
            if (meth == 1) 
                ep = cj;
            else if (meth == 2)
                ep = (cj - bj) /v1;
            else
                ep = (cj - 2.0 * v1 * bj) / (v1 * (v - 2.0));
                         
        }
        else {
            if (meth == 1)
                ej = cj;
            else if (meth == 2)
                ej = (cj - bj) / u1;
            else
                ej = (cj - 2.0 * u * bj) / (u * u1);

            if (ej < eq) {
                eq = ej;
                q = j;
                w = u;
                w1 = u1;
            }
        }
    }

    if (eq >= ep * r)
        goto L9;

    is = 0;
    d = d - ep + eq;
    if (meth == 1) {
        e[p] -= c[p];
        e[q] += c[q];
    }
    else if (meth == 2) {
        e[p] = (v * e[p] - c[p]) / v1;
        e[q] = (w * e[q] + c[q]) / w1;
    }
    else {
        e[p] = (v * v1 * e[p] - c[p]) / (v1 * (v - 2.0));
        e[q] = (w * (w - 1.0) * e[q] + c[q]) / (w * w1);
    }
    mj[p] = l - 1;
    mj[q] += 1;
    z[i] = q;
    goto L10;
}


/* ------------------------------------------------------------------------ */
/*  clp_p       Partitional clustering.                                     */
/*              Method 4: find PMNC points with maximal min-distance        */
/*              by complete enumeration.                                    */  
/*                                                                          */
/*  Return 0 if successful, -1 if insuff. memory.                           */
  

int clp_p(TDAContext *ctx)
{
    register int i,j;
    int n,m;

    n = ctx->GD_NP;
    m = ctx->PMNC;
    if (alloc_acn(ctx, m))                     
        return(-1);          
    if (alloc_acm(ctx, m))                     
        return(-1);          

    for (i = 0; i < m; ++i)     /* current best config */
        ctx->AcM[i] = i;

    ctx->CLP_DMIN = ctx->DBLMAX;
    for (i = 1; i < ctx->PMNC; ++i) {
        for (j = 0; j < i; ++j)  
            ctx->CLP_DMIN = dmin(ctx, ctx->CLP_DMIN,gdd_adj(ctx, i,j,ctx->PMGN));
    }
    printf1(ctx, "CLP_DMIN=%lg\n",ctx->CLP_DMIN);

    clp_rep(ctx, ctx->PMNC,ctx->PMNC,0,n - m);

    printf1(ctx, "Best value: %lg\n",ctx->CLP_DMIN);
    printf1(ctx, "Points: ");
    for (i = 0; i < m; ++i)
        printf1(ctx, "%d ",ctx->AcM[i] + 1);
    newline(ctx);

    return(0);
}

void clp_rep(TDAContext *ctx, int m,int r,int a,int b)
{
    register int i,j,k;
    int cflag;
    double vmin;

    if (a > b)
        return;
    for (k = a; k <= b; ++k) {
           
        cflag = 0;
        for (j = 0; j < m - r; ++j) {
            if (gdd_adj(ctx, j,k,ctx->PMGN) < ctx->CLP_DMIN) {
                cflag = 1;
/*              tda_out(" BREAK bei r=%d k==%d a=%d b=%d\n",r, k,a,b  ); */
                break;
            }
        }
        if (cflag)
            continue;
          
        ctx->AcN[m - r] = k;

/*      tda_out(" r=%d k==%d a=%d b=%d\n",r, k,a,b  );   */
        if (r > 1)
            clp_rep(ctx, m,r-1,k+1,b+1);
        else {
/***************
            tda_out("ACN: ");
            for (j = 0; j < m; ++j)
                tda_out("%d ",AcN[j]);
            newline(ctx);

            tda_out("CURRENT dmin=%lg\n",CLP_DMIN);
*******/
            /* try to find better value */
                  
            vmin = ctx->DBLMAX;
            for (i = 1; i < m; ++i) {
                for (j = 0; j < i; ++j)  
                    vmin = dmin(ctx, vmin,gdd_adj(ctx, ctx->AcN[i],ctx->AcN[j],ctx->PMGN));
            }
/*          printf1("vmin=%lg\n",vmin);     */
            if (vmin > ctx->CLP_DMIN) {
                tda_out("NEW CLP_DMIN=%lg\n",vmin);           
                ctx->CLP_DMIN = vmin;

                for (i = 0; i < m; ++i)
                    ctx->AcM[i] = ctx->AcN[i];
            }
                  
        }
    }

}

/* -##--------------------------------------------------------------------- */
/*  clu     LS fit of ultrametric distances.                                */
/*                                                                          */
/*          clu(                                                            */
/*              gn=...,    graph number, def. 1                             */
/*              mina=...,  minimization algorithm, def. 7                   */
/*                         (possible are: 1, 2, 3, 4, 7)                    */
/*              mxit1=..., maximal number of iteratinos in main loop,       */
/*                         def. 10                                          */
/*              tol=...,   tolerance for convergence, def. 1e-4             */
/*              mxit=...,  maximal number of iterations, default            */
/*                         depends on algorithm                             */
/*              xp=...,    starting values                                  */
/*              dvs=...,   starting values                                  */
/*              df=...,    additional output file containing the            */
/*                         estimated distance matrix                        */
/*              fmt=...,   print format, def. 10.4                          */
/*          ) = fname;                                                      */
/*                                                                          */
/*  The procedure uses: CLMOD1 with clu_fn().                               */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */
       
int clu(TDAContext *ctx)
{
    register int i,j,k;
    int err,n,nn,it;
    double d,t,v;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "LS fit of ultrametric distances). Current memory: %d bytes.\n",ctx->MemReq);

    ctx->MxIt1 = 10;

    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto CLUFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto CLUFin;
    if (gdd_tcheck(ctx, 0,1,2))          /* need undirected valued graph */
        goto CLUFin;

    if (ctx->GD_NP < 4) {
        printf1(ctx, "Error: need at least four nodes.\n");
        goto CLUFin;
    }
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    ctx->CLN = n = ctx->GD_NP;                /* number of points */

    /* setup MDSDIS vector with distances */

    nn = n * (n - 1) / 2;
    if (alloc_acu(ctx, nn + 1))          /* distances */      
        goto CLUFin; 
    if (alloc_act(ctx, nn + 1))          /* previous parameter vector */      
        goto CLUFin; 

    ctx->CLDis = ctx->AcU;

    k = 0;
    for (i = 2; i <= n; ++i) {      /* get distances */
        for (j = 1; j < i; ++j)  
            ctx->CLDis[++k] = gdd_adj(ctx, i-1,j-1,ctx->PMGN);
    }

    tda_out("CLDis ");
    for (i = 1; i <= nn; ++i)
        tda_out("%lg\n",ctx->CLDis[i]);
    newline(ctx);

    /* prepare for min algorithm, always MINA = 7 */

    if (ctx->MINA == 5 || ctx->MINA == 6 || ctx->MINA == 8)
        ctx->MINA = 7;

    ctx->NParm = nn;                  
    ctx->PMDOPT = 1;                     /* analytical gradient */
    ctx->CCTyp = 0;                      /* no covariance matrix */

    set_mlopt(ctx);                    /* adjust options */         

    if (ml_init(ctx, 1,1,0))             /* init function minimization */
        goto CLUFin;

    if (get_dsv(ctx, ctx->NParm,ctx->Par,0,ctx->ParLB,ctx->ParUB,1)) {    /* try to get starting values */
        printf1(ctx, "Error in starting values.\n");
        goto CLUFin;
    }  

    if (ctx->DSVFlg == 0) {              /* create starting values */
        t = d = 0.0;
        for (i = 1; i <= nn; ++i) {
            d += ctx->CLDis[i];
            t += ctx->CLDis[i] * ctx->CLDis[i];
        }
        v = (t / (double)nn) - (d * d / (double)(nn * nn));
        printf1(ctx, "Variance of distances: %lg\n",v);
        v = sqrt(v /= 3.0);

        for (i = 1; i <= nn; ++i)    
            ctx->Par[i] = ctx->CLDis[i] + v * normal(ctx);
    }
    for (i = 1; i <= nn; ++i)   
        tda_out("%10.4lf %10.4lf\n",ctx->CLDis[i],ctx->Par[i]);



    d = clu_ld(ctx, ctx->Par);
    t = clu_pd(ctx, ctx->Par);

    if (fabs(t) <= ctx->EPSI1) {
        printf1(ctx, "\nPenalty term in starting values: %lg\n",t);
        err = 0;
        goto CLUFin;
    }
    ctx->CLRho = d / t;
    if (ctx->CLRho > 1.0)
        ctx->CLRho = 1.0;

    tda_out("ld=%lg pd=%lg rho=%lg\n",d,t,ctx->CLRho);
                  
               
    for (it = 1; it <= ctx->MxIt1; ++it) {   /* main loop */

        printf1(ctx, "Iteration %d  CLRho=%lg\n",it,ctx->CLRho);

        for (i = 1; i <= nn; ++i)       /* save parameter vector */
            ctx->AcT[i] = ctx->Par[i];

        if (ffmin(ctx, CLMOD1,0,1,1))        /* minimization */
            goto CLUFin;                /* insuff memory or fn error */

        prn_mlres(ctx, 1);                   /* print info about minimization */

        /* print result into output file */

        if (ctx->PMFDef) {
            fprintf(ctx->PMFd,"%4d %3d %12.8lf ",it,ctx->LConv,ctx->FMin);
            d = clu_ld(ctx, ctx->Par);
            fprintf(ctx->PMFd,"%12.8lf ",d);
            d = clu_pd(ctx, ctx->Par);
            fprintf(ctx->PMFd,"%12.8lf  ",d);
            for (i = 1; i <= ctx->NParm; ++i)
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->Par[i]);
            fprintf(ctx->PMFd,"\n");
        }

        /* check convergence */

        t = 0.0;
        for (i = 1; i <= nn; ++i) {
            d = ctx->Par[i] - ctx->AcT[i];
            t += d * d;
        }
        t = sqrt(t);
        printf1(ctx, "Change in parameter vector: %lg\n",t);
        if (t <= ctx->PMTol)  
            break;
        ctx->CLRho *= 10.0;
    }
    for (i = 1; i <= nn; ++i) {
        tda_out("%10.4lf %10.4lf %10.4lf\n",ctx->CLDis[i],ctx->AcT[i],ctx->Par[i]);
    }

    d = clu_ld(ctx, ctx->Par);
    printf1(ctx, "Final value of L(D): %lg\n",d);
    d = clu_pd(ctx, ctx->Par);
    printf1(ctx, "Final value of P(D): %lg\n",d);
    d = clu_check(ctx);
    printf1(ctx, "Ultrametric condition: %lg\n",d);

    if (ctx->PMF1Def) {      /* write distance matrix */
        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= n; ++j) {
                if (j < i) 
                    d = ctx->Par[(i - 1) * (i - 2) / 2 + j];
                else if (i < j) 
                    d = ctx->Par[(j - 1) * (j - 2) / 2 + i];
                else
                    d = 0.0;
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,d);
            }
            fprintf(ctx->PMF1d,"\n");
        }
        printf1(ctx, "%d records written to: %s\n",n,ctx->PMF1dName);
    }
    err = 0;

CLUFin:       
    ml_init(ctx, 0,0,0);      
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  clu_ld()    Return L(D).                                                */

double clu_ld(TDAContext *ctx, double *d)
{
    register int i;
    double t,f;
         
    f = 0.0;
    for (i = 1; i <= ctx->NParm; ++i) {
        t = ctx->CLDis[i] - d[i];
        f += t * t;
    }   
    return(f);
}

/* ------------------------------------------------------------------------ */
/*  clu_pd()    Return P(D).                                                */

double clu_pd(TDAContext *ctx, double *d)
{
    register int i,j,k;               
    int ii,jj,kk;
    double f,dij,dik,djk;

    f = 0.0;
    for (i = 2; i <= ctx->CLN; ++i) {
        ii = (i - 1) * (i - 2) / 2;
        for (j = 1; j < i; ++j) {
            jj = (j - 1) * (j - 2) / 2;
            dij = d[ii + j];
            for (k = 1; k <= ctx->CLN; ++k) {
                if (k == i || k == j)
                    continue;

                kk = (k - 1) * (k - 2) / 2;

                if (k < i)
                    dik = d[ii + k];
                else
                    dik = d[kk + i];

                if (dij > dik)
                    continue;

                if (k < j)
                    djk = d[jj + k];
                else
                    djk = d[kk + j];

                if (dij <= djk && fabs(dik - djk) > ctx->EPSI1)
                    f += (dik - djk) * (dik - djk);
            }
        }
    }
    return(f);
}

/* ------------------------------------------------------------------------ */
/*  clu_fn()    Calculate function value for CLMOD1                         */
/*              Use parameter values in TPar[] (i=1,NParm).                 */
/*              Return function value in FTmp, gradient in Grad[].          */
/*                                                                          */
/*              Return 0 if OK, 1 if error in function                      */

int clu_fn(TDAContext *ctx)
{
    register int i,j,k,l;             
    int err,ii,jj,kk;
    double d,dij,dik,djk,fa,fb;
         
    fa = clu_ld(ctx, ctx->TPar);
    fb = clu_pd(ctx, ctx->TPar);
    if (ctx->LFunc)
        ctx->FTmp = fa + ctx->CLRho * fb;                 

    if (ctx->LGrad) {
        for (l = 1; l <= ctx->NParm; ++l)  
            ctx->Grad[l] = 2.0 * (ctx->TPar[l] - ctx->CLDis[l]);
           
        for (l = 1; l <= ctx->NParm; ++l) {
            for (i = 2; i <= ctx->CLN; ++i) {
                ii = (i - 1) * (i - 2) / 2;
                for (j = 1; j < i; ++j) {
                    jj = (j - 1) * (j - 2) / 2;
                    dij = ctx->TPar[ii + j];

                    for (k = 1; k <= ctx->CLN; ++k) {
                        if (k == i || k == j)
                            continue;

                        kk = (k - 1) * (k - 2) / 2;

                        if (k < i)
                            dik = ctx->TPar[ii + k];
                        else
                            dik = ctx->TPar[kk + i];

                        if (dij > dik)
                            continue;

                        if (k < j)
                            djk = ctx->TPar[jj + k];
                        else
                            djk = ctx->TPar[kk + j];

                        if (dij <= djk && fabs(dik - djk) > ctx->EPSI1) {

                            d = ctx->CLRho * 2.0 * (dik - djk);

                            if (k < i && l == (ii + k))               
                                ctx->Grad[l] += d;
                            else if (k > i && l == (kk + i))               
                                ctx->Grad[l] += d;
                            if (k < j && l == (jj + k))               
                                ctx->Grad[l] -= d;
                            else if (k > j && l == (kk + j))               
                                ctx->Grad[l] -= d;
                        }
                    }
                }
            }
        }
    }
    err = checkov(ctx);          /* check overflow */
    if (err) {
        printf1(ctx, "\nError in function evaluation.\n");
        printf1(ctx, "Numerical overflow.\n");
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  clu_check()     Check for ultrametric condition.                        */

double clu_check(TDAContext *ctx)
{
    register int i,j,k;               
    int ii,jj,kk;
    double d,dij,dik,djk;
         
    d = 0.0;
    for (i = 2; i <= ctx->CLN; ++i) {
        ii = (i - 1) * (i - 2) / 2;
        for (j = 1; j < i; ++j) {
            jj = (j - 1) * (j - 2) / 2;
            dij = ctx->Par[ii + j];
            for (k = 1; k <= ctx->CLN; ++k) {
                if (k == i || k == j)
                    continue;

                kk = (k - 1) * (k - 2) / 2;

                if (k < i)
                    dik = ctx->Par[ii + k];
                else
                    dik = ctx->Par[kk + i];
                                           
                if (dij > dik)
                    continue;

                if (k < j)
                    djk = ctx->Par[jj + k];
                else
                    djk = ctx->Par[kk + j];

                if (dij <= djk && fabs(dik - djk) > ctx->EPSI1)
                    d += (dik - djk) * (dik - djk);
            }
        }
    }
    return(d);
}

/* -###-------------------------------------------------------------------- */
/*  clpyr       Pyramidal clustering.                                       */
/*                                                                          */
/*              clpyr(                                                      */
/*                  opt=...,    option, def. 1                              */
/*                              1:   single-link                            */  
/*                              2:   complete-link                          */  
/*      #                       3:   weighted average                       */  
/*      #                       4:   weighted centroid                      */  
/*      #                       5:   group average                          */  
/*      #                       6:   unweighted centroid                    */  
/*      #                       7:   Ward's min variance                    */  
/*                  gn=...,     graph number, def. 1                        */
/*                  df=...,     printf new distance matrix                  */
/*                  pcf=...,    create command file for dendrogram          */
/*                  nfmt=...,   integer print format, def. 4                */
/*                  fmt=...,    print format for values, def. 10.4          */
/*              ) = fname;                                                  */
/*                                                                          */
/*                                                                          */
/*              Requires an undirected valued graph.                        */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int clpyr(TDAContext *ctx)
{
    register int i,j,k,l;
    int err,n,nmax,nn,nc,r = 0,isel,jsel,ip,sflag;
    int *idx_min,*idx_max,*idx_left,*idx_right,*idx_e1,*idx_e2;
    char *idx_mark;
    double d,dm;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Pyramidal clustering. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 5,1,1))       /* get parameters */
        goto CLPYRFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto CLPYRFin;
    if (gdd_tcheck(ctx, 0,0,0))
        goto CLPYRFin;

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    n = ctx->GD_NP;
    nmax = (n * (n + 1)) / 2 + 1;
    nn = nmax * (nmax - 1) / 2 + 2;

    if (alloc_acxf(ctx, nn + 1))         /* distances */   
        goto CLPYRFin;
    if (alloc_aci(ctx, nmax + 1))     
        goto CLPYRFin;
    if (alloc_acj(ctx, nmax + 1))     
        goto CLPYRFin;
    if (alloc_acn(ctx, nmax + 1))        /* number of elements */
        goto CLPYRFin;
    if (alloc_acns(ctx, nmax + 1))       /* times used */
        goto CLPYRFin;
    if (alloc_ack(ctx, nmax + 1))    
        goto CLPYRFin;
    if (alloc_acs(ctx, nmax + 1))    
        goto CLPYRFin;
    if (alloc_acl(ctx, nmax + 1))     
        goto CLPYRFin;
    if (alloc_acr(ctx, nmax + 1))     
        goto CLPYRFin;
    if (alloc_acyf(ctx, nmax + 1))       /* level */
        goto CLPYRFin;
    if (alloc_acc(ctx, n + 1))           /* mark */
        goto CLPYRFin;

    idx_e1 = ctx->AcI;                   /* index of first element */
    idx_e2 = ctx->AcJ;                   /* index of second element */
    idx_min = ctx->AcK;                  /* index of minimal element */
    idx_max = ctx->AcS;                  /* index of maximal element */
    idx_left = ctx->AcL;                 /* index of left neighbour */
    idx_right = ctx->AcR;                /* index of right neighbour */
    idx_mark = ctx->AcC;                 /* temporary membership marks */
   
    k = 1;
    for (i = 1; i < n; ++i) {
        for (j = 0; j < i; ++j)   
            ctx->AcXF[k++] = (float)gdd_adj(ctx, i,j,ctx->PMGN);                       
    }
    for (i = 1; i <= n; ++i) { 
        ctx->AcN[i] = 1;
        idx_min[i] = idx_max[i] = i;
    }


tda_out("n=%d nmax=%d nn=%d\n",n,nmax,nn);

    for (i = 2; i <= n; ++i) {
        ip = (i - 1) * (i - 2) / 2;
        for (j = 1; j < i; ++j)
            tda_out("%g ",(double)(ctx->AcXF[ip + j]));
        newline(ctx);
    }



    nc = n;         /* number of nodes */

    while (1) {
                                 
        dm = ctx->DBLMAX;
        isel = jsel = 0;

        for (i = 2; i <= nc; ++i) {
            if (ctx->AcNS[i] >= 2)
                continue;
            if (idx_left[i] && idx_right[i])
                continue;

            ip = ((i - 1) * (i - 2)) / 2;  

            for (j = 1; j < i; ++j) {
                if (ctx->AcNS[j] >= 2)
                    continue;

                d = (double)ctx->AcXF[ip + j];
                if (d >= dm)  
                    continue;

                /* check order requirements; only if i and/or j are
                   elementary objects */

                if (i <= n || j <= n) {
                    r = clpyr_ord(ctx, i,j,n,idx_left,idx_right,idx_min,idx_max);
                    if (r == 0)
                        continue;
                }
 
                /* check whether i and j are contained in the same cluster */

                r = imax(ctx, ctx->AcN[i],ctx->AcN[j]);

                sflag = 1;
                for (k = 1; k <= nc; ++k) {
                    if (ctx->AcN[k] < r)
                        continue;

                    for (l = 1; l <= n; ++l)
                        idx_mark[l] = 0;

                    clpyr_mark(ctx, k,idx_e1,idx_e2,idx_mark,n);
                    if (clpyr_mcheck(ctx, i,idx_e1,idx_e2,idx_mark)) 
                        continue;
                    if (clpyr_mcheck(ctx, j,idx_e1,idx_e2,idx_mark) == 0) {
                        sflag = 0;
                        break;
                    }
                }
                if (sflag == 0) {
                    tda_out("sflag = 0. continue\n");
                    continue;
                }

                /* check if connected */

                if (clpyr_con(ctx, i,j,idx_min,idx_max,idx_left,idx_right) == 0)
                    continue;

                /* check for additional condition */

                sflag = 1;
                for (l = 1; l <= n; ++l)
                    idx_mark[l] = 0;
                clpyr_mark(ctx, i,idx_e1,idx_e2,idx_mark,n);
                clpyr_mark(ctx, j,idx_e1,idx_e2,idx_mark,n);
                r = 0;
                for (l = 1; l <= n; ++l) {
                    if (idx_mark[l])
                        r++;
                }
                if (r >= 3) {
                    for (k = n + 1; k <= nc; ++k) {
                        if (k == i || k == j || ctx->AcN[k] < 2)
                            continue;

                        r = clpyr_cont(ctx, k,i,idx_e1,idx_e2);
                        tda_out("clpyr_cont k=%d i=%d r=%d\n",k,i,r);

                        if (r)
                            continue;

                        r = clpyr_cont(ctx, k,j,idx_e1,idx_e2);
                        tda_out("clpyr_cont k=%d j=%d r=%d\n",k,j,r);

                        if (r)
                            continue;

                        if (clpyr_mcheck(ctx, k,idx_e1,idx_e2,idx_mark) == 0) {
                            sflag = 0;
                            break;
                        }
                    }
                 
                    if (sflag == 0) {
                        tda_out("NOT OK i=%d j=%d\n",i,j   );
                        continue;
                    }
                }
 

tda_out("vorlaeufig ok: i=%d j=%d d=%g\n",i,j,d  );
                isel = i; 
                jsel = j;
                dm = d;
            }
        }

        if (isel == 0 || jsel == 0)     /* nor more fusions */
            break;

        /* isel and jsel have been selected for fusion */

        if (nc >= nmax) {
            printf1(ctx, "Error: exceeded maximal storage (nmax=%d).\n",nmax);
            goto CLPYRFin;
        }
        nc++;
        idx_e1[nc] = isel;
        idx_e2[nc] = jsel;
        ctx->AcNS[isel] += 1;
        ctx->AcNS[jsel] += 1;

        if (isel <= n || jsel <= n) {
            r = clpyr_ord(ctx, isel,jsel,n,idx_left,idx_right,idx_min,idx_max);
            if (r == 1) {
                if (isel <= n && idx_right[isel] == 0) {
                    k = idx_min[jsel];
                    idx_right[isel] = k;               
                    if (idx_left[k] == 0)
                        idx_left[k] = isel;
                }
                if (jsel <= n && idx_left[jsel] == 0) {
                    k = idx_max[isel];
                    idx_left[jsel] = k;               
                    if (idx_right[k] == 0)
                        idx_right[k] = jsel;
                }
            }
            else if (r == 2) {
                if (isel <= n && idx_left[isel] == 0) {
                    k = idx_max[jsel];
                    idx_left[isel] = k;                   
                    if (idx_right[k] == 0)
                        idx_right[k] = isel;
                }
                if (jsel <= n && idx_right[jsel] == 0) {
                    k = idx_min[isel];
                    idx_right[jsel] = k;                 
                    if (idx_left[k] == 0)
                        idx_left[k] = jsel;
                }
            }
        }
        ctx->AcYF[nc] = (float)dm;

        /* number of elements */

        for (l = 1; l <= n; ++l)
            idx_mark[l] = 0;
        clpyr_mark(ctx, nc,idx_e1,idx_e2,idx_mark,n);
        for (l = 1; l <= n; ++l) {
            if (idx_mark[l])
                ctx->AcN[nc] += 1;
        }


tda_out("nc=%d : ",nc);
for (l = 1; l <= n; ++l)
    tda_out("%d ",idx_mark[l]);
newline(ctx);


        /* update indices of minimal and maximal elements of group */

        for (l = 1; l <= n; ++l) {
            if (idx_mark[l]) {
                k = idx_left[l];
                if (k == 0 || idx_mark[k] == 0)
                    idx_min[nc] = l;
                k = idx_right[l];
                if (k == 0 || idx_mark[k] == 0)
                    idx_max[nc] = l;
            }
        }
        

        tda_out("nc=%2d i=%2d j=%2d dm=%lg acn[nc]=%2d k[nc]=%2d s[nc]=%2d l(i)=%2d r(i)=%2d l(j)=%2d r(j)=%2d r=%d y=%g\n",
                                   nc,isel,jsel,dm,ctx->AcN[nc],idx_min[nc],idx_max[nc],
                                   idx_left[isel],idx_right[isel],idx_left[jsel],idx_right[jsel],r,(double)(ctx->AcYF[nc]));     


        /* update dissimilarity matrix */

        ip = (nc - 1) * (nc - 2) / 2;
        for (k = 1; k <= nc; ++k) {
            ctx->AcXF[++ip] = (float)dmax(ctx, clpyr_dis(ctx, k,isel,ctx->AcXF),clpyr_dis(ctx, k,jsel,ctx->AcXF));
        }
    }

    fprintf(ctx->PMFd,"# 1. cluster number\n");
    fprintf(ctx->PMFd,"# 2. first element\n");
    fprintf(ctx->PMFd,"# 3. second element\n");
    fprintf(ctx->PMFd,"# 4. number of objects\n");
    fprintf(ctx->PMFd,"# 5. minimal element\n");
    fprintf(ctx->PMFd,"# 6. maximal element\n");
    fprintf(ctx->PMFd,"# 7. next on left side\n");
    fprintf(ctx->PMFd,"# 8. next on right side\n");
    fprintf(ctx->PMFd,"# 9. index\n\n");
             
    for (i = 1; i <= nc; ++i) {
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,idx_e1[i]);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,idx_e2[i]);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcN[i]);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,idx_min[i]);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,idx_max[i]);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,idx_left[i]);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,idx_right[i]);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,(double)(ctx->AcYF[i]));
        fprintf(ctx->PMFd,"\n");
    }
    printf1(ctx, "%d records written to: %s\n",nc + 10,ctx->PMFdName);

    /* create order of objects in AcK */

    sflag = 1;
    l = 0;
    for (i = 1; i <= n; ++i) {
        if (idx_left[i] == 0) {
            l = i;
            break;
        }
    }
    if (l == 0)   
        goto CLPYRCon;

    ctx->AcK[1] = l;
    for (i = 2; i <= n; ++i) {
        l = idx_right[l];
        if (l < 1 || l > n)
            goto CLPYRCon;
        ctx->AcK[i] = l;
    }
    sflag = 0;

    printf1(ctx, "Order: ");
    for (i = 1; i <= n; ++i) 
        printf1(ctx, "%d ",ctx->AcK[i]);
    newline(ctx);

    if (ctx->PMF1Def) {   /* print new distance matrix */

        if (alloc_aczf(ctx, n * (n - 1) / 2 + 1))                     
            goto CLPYRFin;

        for (i = 2; i <= n; ++i) {
            for (j = 1; j < i; ++j) {
                k = clpyr_ccl(ctx, i,j,n,nc,idx_e1,idx_e2);
                if (k == 0) {
                    tda_out("Error in data structure clpyr_ccl (i=%d, j=%d).\n",i,j);
                    goto CLPYRFin;
                }
tda_out("i=%d j=%d k=%d d=%lg \n",i,j,k,(double)(ctx->AcYF[k]));
                ctx->AcZF[(i - 1) * (i - 2) / 2 + j] = ctx->AcYF[k];
            }
        }
        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= n; ++j) {
                if (i == j)
                    d = 0.0;
                else if (j < i)
                    d = (double)(ctx->AcZF[(i - 1) * (i - 2) / 2 + j]);
                else             
                    d = (double)(ctx->AcZF[(j - 1) * (j - 2) / 2 + i]);

                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,d);
            }
            fprintf(ctx->PMF1d,"\n");
        }
        printf1(ctx, "%d records written to: %s\n",n,ctx->PMF1dName);
    }


CLPYRCon:
    if (sflag)
        printf1(ctx, "Error: cannot find a valid ordering.\n");

    /* check for compatible ordering */




    if (ctx->PMPCFDef) {                      /* plot dendrogram */
        if (alloc_aczf(ctx, nc + 1))                                
            goto CLPYRFin;
        clpyr_den(ctx, n,nc,ctx->AcK,idx_e1,idx_e2,ctx->AcYF,ctx->AcZF);
    }


    err = 0;

CLPYRFin:       
    p_clean(ctx);
    return(err);
}
   
/* return distance between node i and j */

double clpyr_dis(TDAContext *ctx, int i,int j,float *d)
{
    (void)ctx;        /* unused: the signature is shared */
    int ip;

    if (i == j)
        return(0.0);

    if (j < i)
        ip = ((i - 1) * (i - 2)) / 2 + j;
    else               
        ip = ((j - 1) * (j - 2)) / 2 + i;
    return((double)d[ip]);
}

/* return 1 if the fusion of i and j is connected */

int clpyr_con(TDAContext *ctx, int i,int j,int *idx_min,int *idx_max,int *idx_left,int *idx_right)
{
    int il,ir,jl,jr;

    il = idx_min[i];
    jl = idx_min[j];
    ir = idx_max[i];
    jr = idx_max[j];

tda_out("clpyr_con: i=%d j=%d il=%d ir=%d jl=%d jr=%d ",i,j,il,ir,jl,jr);


    if (clpyr_left(ctx, ir,jl,idx_left)) {
        tda_out("ir left from jl\n");
        tda_out("idx_right(ir) = %d\n",idx_right[ir]);
                             
        if (idx_right[ir] != jl) {
tda_out(" return 0\n");
            return(0);
        }
    }
    if (clpyr_left(ctx, jr,il,idx_left)) {
        tda_out("jr left from il\n");
        tda_out("idx_right(jr) = %d\n",idx_right[jr]);
                             
        if (idx_right[jr] != il) {
tda_out(" return 0\n");
              return(0);
        }
    }
tda_out(" return 1\n");
    return(1);
}

/* return 1 if i << j possible, 2 if j << i possible, 0 if nothing possible */

int clpyr_ord(TDAContext *ctx, int i,int j,int n,int *idx_left,int *idx_right,int *idx_min,int *idx_max)
{
    int il,ir,jl,jr;

    if (i <= n)
        il = ir = i;
    else {
        il = idx_min[i];
        ir = idx_max[i];
    }
    if (j <= n)
        jl = jr = j;
    else {
        jl = idx_min[j];
        jr = idx_max[j];
    }

    if (i <= n) {
        if (idx_left[i] && idx_right[i]) {
            return(0);        
        }
        else if (idx_left[i]) {
            if (j <= n) {
                if (idx_left[j] || clpyr_left(ctx, j,i,idx_left))
                    return(0);
                else
                    return(1);
            }
            else {
                if (idx_left[jl] || clpyr_left(ctx, jl,i,idx_left))
                    return(0);
                else
                    return(1);
            }
        }
        else if (idx_right[i]) {
            if (j <= n) {
                if (idx_right[j] || clpyr_right(ctx, j,i,idx_right))
                    return(0);
                else
                    return(2);
            }
            else {
                if (idx_right[jr] || clpyr_right(ctx, jr,i,idx_right))
                    return(0);
                else
                    return(2);
            }
        }
        else {
            if (j <= n) {
                if (idx_left[j] == 0 && clpyr_left(ctx, i,j,idx_left) == 0)
                    return(1);
                else if (idx_right[j] == 0 && clpyr_right(ctx, i,j,idx_right) == 0)
                    return(2);
                else
                    return(0);
            }
            else {
                if (idx_left[jl] == 0 && clpyr_left(ctx, i,jl,idx_left) == 0)
                    return(1);
                else if (idx_right[jr] == 0 && clpyr_right(ctx, i,jr,idx_right) == 0)
                    return(2);
                else
                    return(0);
            }
        }
    }
    else {      /* i > n */ 
        if (j > n)
            return(1);

        /* j <= n */

        if (idx_left[j]) {
            if (idx_left[il] || clpyr_left(ctx, j,il,idx_left))
                return(0);
            else
                return(1);
        }
        else if (idx_right[j]) {
            if (idx_right[ir] || clpyr_right(ctx, j,ir,idx_right))
                return(0);
            else
                return(2);
        }
        else {
            if (idx_left[il] == 0 && clpyr_left(ctx, j,il,idx_left) == 0)
                return(2);
            else if (idx_right[ir] == 0 && clpyr_right(ctx, j,ir,idx_right) == 0)
                return(1);
            else
                return(0);
        }
    }
}

/* return 1 if i is placed left of j */                            

int clpyr_left(TDAContext *ctx, int i,int j,int *idx_left)
{
    (void)ctx;        /* unused: the signature is shared */
    while (j) {
        j = idx_left[j];
        if (j == i)
            return(1);
    }
    return(0);
}

/* return 1 if i is placed right of j */                            

int clpyr_right(TDAContext *ctx, int i,int j,int *idx_right)
{
    (void)ctx;        /* unused: the signature is shared */
    while (j) {
        j = idx_right[j];
        if (j == i)
            return(1);
    }
    return(0);
}

/* return 1 if i is not contained in idx_mark */

int clpyr_mcheck(TDAContext *ctx, int i,int *idx_e1,int *idx_e2,char *idx_mark)
{
    register int ii,jj;

    if (i == 0)
        return(0);
    ii = idx_e1[i];
    jj = idx_e2[i];
    if (ii) {
        if (clpyr_mcheck(ctx, ii,idx_e1,idx_e2,idx_mark))
            return(1);
        if (clpyr_mcheck(ctx, jj,idx_e1,idx_e2,idx_mark))
            return(1);
    }
    else if (idx_mark[i] != 1)
        return(1);
    return(0);
}

/* mark idx_mark[l] = 1 if l is contained in i */

void clpyr_mark(TDAContext *ctx, int i,int *idx_e1,int *idx_e2,char *idx_mark,int n)
{
    register int ii,jj;

    if (i == 0)
        return;
    ii = idx_e1[i];
    jj = idx_e2[i];
    if (ii) {
        clpyr_mark(ctx, ii,idx_e1,idx_e2,idx_mark,n);
        clpyr_mark(ctx, jj,idx_e1,idx_e2,idx_mark,n);
    }
    else  
        idx_mark[i] = 1;            
    return;
}

/* return 1 if i is contained in j */

int clpyr_cont(TDAContext *ctx, int i,int j,int *idx_e1,int *idx_e2)
{
    register int j1,j2;

    if (i == j)
        return(1);
    if (j == 0)
        return(0);
    j1 = idx_e1[j];
    j2 = idx_e2[j];
    if (j1 == 0)
        return(0);

    if (clpyr_cont(ctx, i,j1,idx_e1,idx_e2))
        return(1);
    if (clpyr_cont(ctx, i,j2,idx_e1,idx_e2))
        return(1);
    return(0);                  
}

/* return index of lowest cluster which contains i and j */

int clpyr_ccl(TDAContext *ctx, int i,int j,int n,int nc,int *idx_e1,int *idx_e2)
{
    register int k;
         
    for (k = n + 1; k <= nc; ++k) {
        if (clpyr_cont(ctx, i,k,idx_e1,idx_e2) && clpyr_cont(ctx, j,k,idx_e1,idx_e2))
            return(k);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  clpyr_den   Create a command file for plotting the pyramid.             */

void clpyr_den(TDAContext *ctx, int n,int nc,int *iord,int *idx_e1,int *idx_e2,float *acy,float *tmp)
{
    register int i,j,k;
    double d,hy,h,hl,hr,xl,xr; 

    hy = 0.0;
    for (i = 1; i <= nc; ++i) {
        if ((double)(hy) < (double)(acy[i]))
            hy = (double)(acy[i]);
    }
    fprintf(ctx->PMPCFd,"psfile = %s.ps;\n",ctx->PMPCFName);
    fprintf(ctx->PMPCFd,"psetup(pxa=0,%d,pya=0,%g);\n",n + 1,hy);
    fprintf(ctx->PMPCFd,"plxa;\n");
    fprintf(ctx->PMPCFd,"plya;\n");
    d = -hy / 20;
    for (i = 1; i <= n; ++i)  
        fprintf(ctx->PMPCFd,"pltext(xy=%d,%lg,fs=1.5)=%d;\n",i,d,iord[i]);

    for (i = 1; i <= n; ++i)
        tmp[iord[i]] = (float)i;

    for (k = n + 1; k <= nc; ++k) {
        i = idx_e1[k];
        j = idx_e2[k];

        h = (double)(acy[k]);
        xl = (double)tmp[i];
        hl = (double)acy[i];
        xr = (double)tmp[j];
        hr = (double)acy[j];
        d = (xr - xl) / 4.0;

        fprintf(ctx->PMPCFd,"plotp=%g,%g,%g,%g,%g,%g,%g,%g;\n",    
            xl,hl,xl + d,h,xr - d,h,xr,hr);

        tmp[k] = (float)((xl + xr) / 2.0);

    }

    printf1(ctx, "Output file for plot of dendrogram: %s\n",ctx->PMPCFName);
}

