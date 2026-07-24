/****************************************************************************/
/*  t_mds                                                                   */
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
#include "t_ml.h"
#include "t_rand.h"
#include "t_imat.h"
#include "t_freq.h"
#include "t_com.h"
#include "t_plot.h"
#include "t_xplot.h"
#include "t_graph.h"
#include "t_min.h"
#include "t_svd.h"
#include "t_mes.h"
#include "t_ass.h"
#include "t_lp.h"
#include "tda_context.h"

/*  functions in t_mds.c */

int dmet(TDAContext *ctx); 
int dmet1(TDAContext *ctx); 
double dmet_c(TDAContext *ctx, int n,float *d);
void m_spath(TDAContext *ctx, int n,float *d,float *a);

int gpro(TDAContext *ctx);
void gpro_pcf(TDAContext *ctx, int n,double *x);
int gpro_sval(TDAContext *ctx, int n,double *x,int *idx);
double gpro_fun(TDAContext *ctx, int n,double *x);
void gpro_df(TDAContext *ctx, int n,double *x);
int gpro_min(TDAContext *ctx, int n,double *par,double *grad,double *step);

int uds(TDAContext *ctx);
void uds_check(TDAContext *ctx, int n,double *a,double *b);
void uds_ncheck(TDAContext *ctx, int n,int *a,int *b);
int uds_1(TDAContext *ctx);  
int uds_0(TDAContext *ctx, int n,int k);
double uds_get_dip(TDAContext *ctx, int n,int i,int *p);
double uds_get_dpi(TDAContext *ctx, int n,int i,int *p);
double uds_get_f(TDAContext *ctx, int n,int *p);
int uds_r(TDAContext *ctx, int n,int k);
int uds_r0(TDAContext *ctx, int n,int k,double ck,double ckk0,double cl);
int uds_r1(TDAContext *ctx, int n,int k,int l,double ck,double cl,double cll0);
int uds_qa1(TDAContext *ctx);  
int uds_qa2(TDAContext *ctx);  

int sga(TDAContext *ctx);
int sga1(TDAContext *ctx, int m,int *r);
int hpo(TDAContext *ctx, int n,int m,int *x);
int hpo1(TDAContext *ctx, int m,int *x,int *y);

int unf(TDAContext *ctx);
void unf1(TDAContext *ctx, int n,int m,int *p,int *r,double *x,int *s,int *rr,int prn, short *pp,int wf,double *w);
int unf2(TDAContext *ctx, int m,int *r,int *s,int *rr);

int mdsc(TDAContext *ctx);
int mdsm(TDAContext *ctx);
int mds4_fn(TDAContext *ctx);


int mdsn(TDAContext *ctx);
void mdsn_dis(TDAContext *ctx, double *dis);
int mdsn_norm(TDAContext *ctx);
int mdsn_grad(TDAContext *ctx, double s,double t,double stress);   

int mdsn1(TDAContext *ctx);
void mdsn1_dis(TDAContext *ctx, int n,double *par,double *dis);
void mdsn1_proj(TDAContext *ctx, int opt,int nn,double *dis,double *tmp,double *proj);
int mdsn1_fn(TDAContext *ctx, int n,int nn,double *par,double *fval,double *grad,double *dis,double *disp);

double mds_edis(TDAContext *ctx, double x1,double y1,double x2,double y2);
double mds_elen(TDAContext *ctx, double x,double y);

int mdsx(TDAContext *ctx);
int mds_m(TDAContext *ctx, int ip,int ip0,int ip1);
int mds3_fn(TDAContext *ctx);
int mds_r(TDAContext *ctx, int n,float *d,double *x,double *y,double *r,int *p,double *ds, int *is,int *js,int *ptr);
void mds_p(TDAContext *ctx, int n,double *x,double *y,double *r);
void mds_p1(TDAContext *ctx, int n,double *x,double *y);
int mdsxr_m(TDAContext *ctx, int narg,int nbmax,double *par,double *lb,double *ub,double *par0);
int mdsxr_f(TDAContext *ctx, int opt,int n,double *xl,double *xh,double *f);
double mdsxr_fm(TDAContext *ctx);
int mdsxr_f1(TDAContext *ctx, int opt,int n,double *xl,double *xh,double *f,double fminp);
double mdsxr_dm(TDAContext *ctx, double ai,double bi,double aj,double bj);
int mdsxr_c2(TDAContext *ctx, int n,double *lb,double *ub,double *x,double *y);
int mdsxr_c3(TDAContext *ctx, int n,double *lb,double *ub,double *x,double *y);

int mdsr(TDAContext *ctx);
double mdsr_f(TDAContext *ctx, int n,double *x,double *y);

int rfit(TDAContext *ctx);
int rfit1(TDAContext *ctx);

/* ------------------------------------------------------------------------ */



/* ------------------------------------------------------------------------ */
/*  dmet        Change distance matrix intro metric.                        */
/*              Requires undirected valued graph.                           */
/*                                                                          */
/*              dmet(                                                       */
/*                  gn = ...,   graph number, def. 1                        */
/*                  fmt=...,    printf format, def. 10.4                    */
/*              ) = fname;                                                  */
/*                                                                          */
/*  The command calculates                                                  */
/*  c = max{0, max_i,j,k {d_i,j - d_i,k - d_kj}}                            */
/*  and changes the distances d_ij into d_ij + c.                           */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */
       
int dmet(TDAContext *ctx)
{
    register int i,j,k;
    int err,n;
    double a,dij,dm;              

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Change distance matrix into metric. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLFD = 1.0e-4;                 /* default box width */

    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto DMETFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto DMETFin;
    if (gdd_tcheck(ctx, 0,1,2))          /* need undirected valued graph */
        goto DMETFin;

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    n = ctx->GD_NP;         
    if (n < 3) {
        printf1(ctx, "Nothing done for less than 3 points.\n");
        err = 0;
        goto DMETFin;
    }
    if (alloc_acx(ctx, n * n + 1))  
        goto DMETFin;

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (i == j)
                ctx->AcX[i * n + j] = 0.0;
            else
                ctx->AcX[i * n + j] = gdd_adj(ctx, i,j,ctx->PMGN);
        }
    }
    fprintf(ctx->PMFd,"Original distance matrix\n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcX[i * n + j]);
#ifdef TDA_R_PACKAGE
            tda_export_cell(ctx, "dmet.original", ctx->AcX[i * n + j]);
#endif
        }
        fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
        tda_export_endrow(ctx, "dmet.original");
#endif
    }
    dm = 0.0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (i == j)
                continue;
            dij = ctx->AcX[i * n + j];
            for (k = 0; k < n; ++k) {
                if (k == i || k == j)
                    continue;
                a = dij - ctx->AcX[i * n + k] - ctx->AcX[k * n + j];
                if (dm < a) {
tda_out("i=%d j=%d k=%d a=%lg dm=%lg ",i+1,j+1,k+1,a,dm);
                    dm = a;
tda_out(" new dm=%lg\n",dm); 

                }
            }
        }
    }
    printf1(ctx, "Additive constant: %lg\n",dm);

    fprintf(ctx->PMFd,"Modified distance matrix\n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            a = ctx->AcX[i * n + j];
            if (i != j)       
                a += dm;
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,a);
#ifdef TDA_R_PACKAGE
            tda_export_cell(ctx, "dmet.modified", a);
#endif
        }
        fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
        tda_export_endrow(ctx, "dmet.modified");
#endif
    }
#ifdef TDA_R_PACKAGE
    /* the additive constant the console reports, so the reader need
       not regex it back out of the output */
    tda_export_row(ctx, "dmet.constant", &dm, 1);
#endif
    err = 0;

DMETFin:
    p_clean(ctx);
    return(err);
}

       
/* ------------------------------------------------------------------------ */
/*  dmet1       Change distance matrix intro metric.                        */
/*              Requires undirected valued graph.                           */
/*                                                                          */
/*              dmet1(                                                      */
/*                  gn = ...,   graph number, def. 1                        */
/*                  tolfd=...,  tolerance for max deviation, def. 1.e-4     */
/*                  mxit=...,   max number of iterations, def. 20           */  
/*                  fmt=...,    printf format, def. 10.4                    */
/*              ) = fname;                                                  */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */
       
int dmet1(TDAContext *ctx)
{
    register int i,j,k;
    int err,n,iter;
    double dev;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Change distance matrix into metric. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLFD = 1.0e-4;                 /* default box width */

    if (parm(ctx, ctx->CmdBuf + 5,1,1))       /* get parameters */
        goto DMET1Fin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto DMET1Fin;
    if (gdd_tcheck(ctx, 0,1,2))          /* need undirected valued graph */
        goto DMET1Fin;

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);
    if (ctx->MxItFlg == 0)
        ctx->MxIter = 20;

    printf1(ctx, "Maximal number of iterations: %d\n",ctx->MxIter);
    printf1(ctx, "Tolerance for maximal deviation: %g\n\n",ctx->TOLFD);

    n = ctx->GD_NP;         
    if (n < 3) {
        printf1(ctx, "Nothing done for less than 3 points.\n");
        err = 0;
        goto DMET1Fin;
    }
    if (alloc_acxf(ctx, n * n + 1))  
        goto DMET1Fin;
    if (alloc_acyf(ctx, n * n + 1))  
        goto DMET1Fin;

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (i == j)
                ctx->AcXF[i * n + j] = 0.0;
            else
                ctx->AcXF[i * n + j] = (float)gdd_adj(ctx, i,j,ctx->PMGN);
        }
    }
    fprintf(ctx->PMFd,"Original distance matrix\n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j)
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,(double)(ctx->AcXF[i * n + j]));
        fprintf(ctx->PMFd,"\n");
    }
    dev = dmet_c(ctx, n,ctx->AcXF);     

    printf1(ctx, "Maximal deviation in input data: %g\n",dev);
    if (dev <= ctx->TOLFD) {
        err = 0;                 /* already metric: success, not an error */
        goto DMET1Fin;
    }

    for (iter = 1; iter <= ctx->MxIter; ++iter) {

        m_spath(ctx, n,ctx->AcXF,ctx->AcYF);       /* find shortest paths */

        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; ++j) {
                if (i != j) {
                    k = i * n + j;
                    ctx->AcXF[k] = (float)dmin(ctx, (double)ctx->AcXF[k],(double)ctx->AcYF[k]);
                }
            }
        }
        dev = dmet_c(ctx, n,ctx->AcXF);
/*      printf1("dev for modified X = %g\n",dev);   */

        if (dev <= ctx->TOLFD)
            break;
    }
    if (iter > ctx->MxIter)
        iter--;

    printf1(ctx, "Number of iterations performed: %d\n",iter);
    printf1(ctx, "Remaining deviation: %g\n\n",dev);

    fprintf(ctx->PMFd,"Modified distance matrix\n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j)
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,(double)(ctx->AcXF[i * n + j]));
        fprintf(ctx->PMFd,"\n");
    }
    err = 0;

DMET1Fin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  dmet_c(n,d)     Check d[] for metric. Return max deviation.             */
       
double dmet_c(TDAContext *ctx, int n,float *d)
{
    register int i,j,k;
    double dij,dik,dkj,dev;
                  
    dev = 0.0;
    for (i = 1; i < n; ++i) {
        for (j = 0; j < i; ++j) {
            dij = (double)(d[i * n + j]);
            for (k = 0; k < n; ++k) {
                if (k == i || k == j)
                    continue;
                dik = (double)(d[i * n + k]);
                dkj = (double)(d[k * n + j]);
                dev = dmax(ctx, dev,dij - dik - dkj);
            }
        }
    }
    return(dev);
}

/* ------------------------------------------------------------------------ */
/*  m_spath(n,d,a)                                                          */
/*                                                                          */
/*  Calculate shortes patht matrix for a valued graph.                      */
/*                                                                          */
/*  Algorithm based on R.W. Floyd, Shortest Path.                           */
/*  Collected Algorithms from ACM 97.                                       */

void m_spath(TDAContext *ctx, int n,float *d,float *a)
{
    register int i,j,k,ni,nj;
    double s,max;

    max = 0.0; 
    for (i = 1; i < n; ++i) {
        for (j = 0; j < i; ++j) 
            max = (double)(dmax(ctx, max,(double)d[i * n + j]));
    }
    max *= 10.0;

    for (i = 0; i < n; ++i) {
        ni = i * n;
        for (j = 0; j < n; ++j) {
            a[ni] = (float)(max);
            if (i != j) {
                s = (double)(d[i * n + j]);                                        
                if (s >= 0.0)
                    a[ni] = (float)(s);
            }
            ni++;
        }
    }
    for (i = 0; i < n; ++i) {
        ni = i * n;
        for (j = 0; j < n; ++j) {
            nj = j * n;
            if ((double)(a[nj + i]) < max) {
                for (k = 0; k < n; ++k) {
                    if ((double)(a[ni + k]) < max) {
                        s = (double)(a[nj + i] + a[ni + k]);
                        if ((double)(s) < (double)(a[nj + k]))
                            a[nj + k] = (float)(s);
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
            else if ((double)(a[ni]) >= max)  
                a[ni] = -1.0;
            ni++;
        }
    }
}


/* ------------------------------------------------------------------------ */
/*  gpro        General projection of proximities.                          */
/*              Requires an undirected valued graph.                        */
/*                                                                          */
/*              gpro(                                                       */
/*                  alg = ...,  algorithm, def. 1                           */
/*                              1 = direct search                           */  
/*                  gn = ...,   graph number, def. 1                        */
/*                  df=...,     print distances                             */
/*                  pcf=...,    cf file for plot of configuration           */
/*                  nfmt = ..., integer format, def. 4                      */
/*                  fmt = ...,  print format for values, def. 10.4          */
/*                  mxit=...,   max number of iterations, def. 50           */
/*                  tols=...,   min of reduction factor, def. 10e-4         */
/*                  slen=...,   step length, def. 1                         */
/*                  sred=...,   reduction factor, def. 0.5                  */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gpro(TDAContext *ctx)
{
    register int i;
    int err,n,nn,nrec;
    double x,y;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Projection of proximities. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto GPROFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto GPROFin;
    if (gdd_tcheck(ctx, 0,1,2))          /* need undirected valued graph */
        goto GPROFin;

    if (ctx->GD_NP < 3) {
        printf1(ctx, "Error: need at least three nodes.\n");
        goto GPROFin;
    }

    if (ctx->PMALG != 1)
        ctx->PMALG = 1;

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);
    if (ctx->MxItFlg == 0)
        ctx->MxIter = 50;

    printf1(ctx, "\nAlgorithm %d: ",ctx->PMALG);
    if (ctx->PMALG == 1) {
        printf1(ctx, "direct search.\n");
        printf1(ctx, "Starting value of step length: %g\n",ctx->SLen);
        printf1(ctx, "Step length reduction factor: %g\n",ctx->SRed);
        printf1(ctx, "Minimum value of step length: %g\n",ctx->TOLS);
    }
    printf1(ctx, "Max number of iterations: %d\n\n",ctx->MxIter);

    if (alloc_acx(ctx, 2 * ctx->GD_NP + 1))       /* coordinates */
        goto GPROFin;
    if (alloc_acy(ctx, 2 * ctx->GD_NP + 1))      
        goto GPROFin;
    if (alloc_actmp(ctx, 2 * ctx->GD_NP + 1))   
        goto GPROFin;
    if (alloc_acn(ctx, ctx->GD_NP + 1))   
        goto GPROFin;

    ctx->MDSIter = ctx->MDSFCalls = 0;
    n = ctx->GD_NP - 1;

    /* calculate starting values */

    if (gpro_sval(ctx, n,ctx->AcX,ctx->AcN)) {
        printf1(ctx, "Error: graph not connected.\n");
        goto GPROFin;
    }
    nn = 2 * n;
    gpro_min(ctx, nn,ctx->AcX,ctx->AcY,ctx->AcTmp);

    printf1(ctx, "Performed iterations: %d (%d function calls)\n",ctx->MDSIter,ctx->MDSFCalls);
    printf1(ctx, "Final stress: %g\n\n",ctx->MDSFN);

    nrec = 0;
    for (i = 0; i <= n; ++i) {
        if (i == 0)
            x = y = 0.0;
        else {
            x = ctx->AcX[2 * i - 1];
            y = ctx->AcX[2 * i];
        }
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,x);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,y);
        fprintf(ctx->PMFd,"\n");
        nrec++;
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);

    if (ctx->PMF1Def)            /* write estimated distances to PMF1d */
        gpro_df(ctx, n,ctx->AcX);

    if (ctx->PMPCFDef)           /* create cf file for plot of configuration */
        gpro_pcf(ctx, n,ctx->AcX);

    err = 0;

GPROFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gpro_pcf    Create a command file for plotting the configuration.       */

void gpro_pcf(TDAContext *ctx, int n,double *x)
{
    register int i;
    double xmin,xmax,ymin,ymax,xd,yd;

    xmin = ymin = xmax = ymax = 0.0;
    for (i = 1; i <= n; ++i) {
        xmin = dmin(ctx, xmin,x[2 * i - 1]);
        xmax = dmax(ctx, xmax,x[2 * i - 1]);
        ymin = dmin(ctx, ymin,x[2 * i]);
        ymax = dmax(ctx, ymax,x[2 * i]);
    }
    xd = (xmax - xmin) / 10.0;
    yd = (ymax - ymin) / 10.0;

    fprintf(ctx->PMPCFd,"psfile = %s.ps;\n",ctx->PMPCFName);
    fprintf(ctx->PMPCFd,"psetup(pxa=%g,%g,pya=%g,%g);\n",
                      xmin - xd,xmax + xd,ymin - yd,ymax + yd);

    fprintf(ctx->PMPCFd,"plg(\n  nmax=%d,\n",n + 1);
    for (i = 0; i <= n; ++i) {
        if (i == 0)  
            xd = yd = 0.0;
        else {
            xd = x[2 * i - 1];
            yd = x[2 * i];
        }
        fprintf(ctx->PMPCFd,"  node(n=%d)=%g,%g,\n",i + 1,xd,yd);
    }
    fprintf(ctx->PMPCFd,");\n");
    printf1(ctx, "Output file for plot of configuration: %s\n",ctx->PMPCFName);
}

/* ------------------------------------------------------------------------ */
/*  gpro_sval(n,x,idx)  Calcuate starting values for function minimization. */
/*                      Return 0 if ok, -1 if graph not connected.          */

int gpro_sval(TDAContext *ctx, int n,double *x,int *idx)
{
    register int i,j;
    int fnd;
    double dij,xi,xj,yi,b,c,tmp;


    for (i = 1; i <= n; ++i)   
        idx[i] = 0;

    x[0] = 0.0;
    idx[0] = 1;
    fnd = 0;
    for (j = 1; j <= n; ++j) {
        dij = gdd_adj(ctx, 0,j,ctx->PMGN);
        if (dij >= 0.0) {
            x[2 * j - 1] =  dij;
            x[2 * j] = 0.0;
            idx[j] = 1;
            fnd = 1;
            break;
        }
    }
    if (fnd == 0)
        return(-1);

    /*  x is indexed 1-based here -- x[2 * i - 1] and x[2 * i] -- so i
        starting at 0 read one double BEFORE the array, which
        read one double before the array on gpro.cf.

        NOTE the idx[0] test below: it does not depend on i, so it either
        skips the whole loop or none of it.  That looks like a typo for
        idx[i], but changing it would change results, and there is no
        evidence here for what it should be -- left as Rohwer wrote it,
        recorded in CONTRIBUTING.md.  */
    for (i = 1; i <= n; ++i) {
        if (idx[0] == 0)
            continue;
        for (j = i + 1; j <= n; ++j) {
            if (idx[j])
                continue;
            dij = gdd_adj(ctx, i,j,ctx->PMGN);
            if (dij < 0.0)
                continue;
            xi = x[2 * i - 1];
            yi = x[2 * i];
            c  = 0.5 * (dij * dij - xi * xi - yi * yi);
            b  = xi + yi;
            tmp = b * b + 4.0 * c;
            if (tmp >= 0.0)  
                xj = (b + sqrt(tmp)) / 2.0;
            else
                xj = random1(ctx);
            x[2 * j - 1] = x[2 * j] = xj;
            idx[j] = 1;
            fnd++;
        }
    }
    if (fnd != n)
        return(-1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  gpro_fun(n,x)       Calcuates the function value and gradient for       */
/*                      gpro. n is 2 * number of nodes. x is the current    */
/*                      configuration. The function returns the function    */
/*                      value calculated at given configuration x, and      */
/*                      the gradient in g[], only if g != NULL              */

double gpro_fun(TDAContext *ctx, int n,double *x)
{
    register int i,j;
    int m;
    double dij,aij,xij,yij,f,f1;    

    ctx->MDSFCalls++;
    f1 = f = 0.0;
    m = n / 2;
    for (i = 0; i < m; ++i) {
        for (j = i + 1; j <= m; ++j) {
            dij = gdd_adj(ctx, i,j,ctx->PMGN);
            if (dij < 0.0)
                continue;

            if (i == 0) {
                xij = -x[2 * j - 1];
                yij = -x[2 * j];
            }
            else {
                xij = x[2 * i - 1] - x[2 * j - 1];
                yij = x[2 * i    ] - x[2 * j    ];
            }
            aij = sqrt(xij * xij + yij * yij);
            f += (aij - dij) * (aij - dij);
            f1 += aij * aij;
        }
    }
    if (f1 > 0.0)
        f /= f1;
    if (f > 0.0)
        f = sqrt(f);
    return(f);
}

/* ------------------------------------------------------------------------ */
/*  gpro_df(n,par)  write estimated distances to PMF1d.                     */

void gpro_df(TDAContext *ctx, int n,double *x)
{
    register int i,j,nrec;
    double dij,aij,xij,yij;

    nrec = 0;
    for (i = 0; i < n; ++i) {
        for (j = i + 1; j <= n; ++j) {
            dij = gdd_adj(ctx, i,j,ctx->PMGN);
            if (dij < 0.0)
                continue;

            if (i == 0) {
                xij = -x[2 * j - 1];
                yij = -x[2 * j];
            }
            else {
                xij = x[2 * i - 1] - x[2 * j - 1];
                yij = x[2 * i    ] - x[2 * j    ];
            }
            aij = sqrt(xij * xij + yij * yij);

            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,i + 1);
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,j + 1);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,dij);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,aij);
            fprintf(ctx->PMF1d,"\n");
            nrec++;
        }
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);
}

/* ------------------------------------------------------------------------ */
/*  gpro_min(n,par,par1,step)                                               */
/*                                                                          */
/*      Function minimization with a direct search method, used by gpro.    */
/*                                                                          */
/*      Return:     0 if successful.                                        */
/*                 -1 if exceeded max number of iterations.                 */

int gpro_min(TDAContext *ctx, int n,double *par,double *par1,double *step)
{
    register int k;
    double delta,fm,fm1,theta,tmp;

    if (ctx->SILENTFlg < 2)
        printfe(ctx, "\n  Iter    Function Value       Step Length     Par Change   FCall\n");

    delta = ctx->SLen;
    for (k = 1; k <= n; ++k) {
        tmp = fabs(par[k]);
        if (tmp < ctx->EPSI1)
            step[k] = delta;
        else
            step[k] = delta * tmp;
    }
    ctx->MDSFN = gpro_fun(ctx, n,par);

FM1L1:  
    fm = ctx->MDSFN;
    for (k = 1; k <= n; ++k)
        par1[k] = par[k];

    for (k = 1; k <= n; ++k) {
        par1[k] += step[k];
        fm1 = gpro_fun(ctx, n,par1);

        if (fm1 < fm) 
            fm = fm1;
        else {
            step[k] = -step[k];
            par1[k] += 2.0 * step[k];

            fm1 = gpro_fun(ctx, n,par1);

            if (fm1 < fm)
                fm = fm1;
            else
                par1[k] -= step[k];
        }
    }
    if (fm < ctx->MDSFN) {

FM1L2:
        for (k = 1; k <= n; ++k) {
            if ((par1[k] > par[k] && step[k] <  0.0) ||
                                            (par1[k] <= par[k] && step[k] >= 0.0))   
                step[k] = -step[k];
            theta = par[k];
            par[k] = par1[k];
            par1[k] = 2.0 * par1[k] - theta;
        }
        ctx->MDSFN = fm;
        fm = fm1 = gpro_fun(ctx, n,par1);

        for (k = 1; k <= n; ++k) {
            par1[k] += step[k];

            fm1 = gpro_fun(ctx, n,par1);

            if (fm1 < fm) 
                fm = fm1;
            else {
                step[k] = -step[k];
                par1[k] += 2.0 * step[k];

                fm1 = gpro_fun(ctx, n,par1);
                if (fm1 < fm)
                    fm = fm1;
                else
                    par1[k] -= step[k];
            }
        }
        if (fm >= ctx->MDSFN)
            goto FM1L1;

        for (k = 1; k <= n; ++k) {
            if (fabs(par1[k] - par[k]) > 0.5 * fabs(step[k])) 
                goto FM1L2;
        }
    }
    ctx->MDSIter++;

    if (ctx->SILENTFlg < 2) {
        printfe(ctx, "  %3d  %20.13e %17.10e ",ctx->MDSIter,ctx->MDSFN,delta);
        printfe(ctx, "        --   %6d\n",ctx->MDSFCalls);
    }
    if (delta > ctx->TOLS) {
        if (ctx->MDSIter >= ctx->MxIter) {
            ctx->LConv = 2;
            goto FM1Fin;
        }
        delta *= ctx->SRed;
        for (k = 1; k <= n; ++k)
            step[k] *= ctx->SRed;
        goto FM1L1;
    }
    else  
        ctx->LConv = 0;

FM1Fin:
    return(0);
}  

/* -##--------------------------------------------------------------------- */
/*  uds         Uni-dimensional scaling.                                    */
/*              Requires an undirected valued graph.                        */
/*                                                                          */
/*              uds(                                                        */
/*                  alg=...,    algorithm, def. 1                           */
/*                              1 = branch and bound (Defays 1978)          */
/*                              2 = approx. quadr. assignment (CACM 608)    */
/*                              3 = approx. quadr. assignment (CACM 754)    */
/*                                  (presupposes integer-valued distance    */
/*                                   matix)                                 */
/*                  gn = ...,   graph number, def. 1                        */
/*                  mxit=...,   alg 3: max iterations, def. 100             */
/*                  idf=...,    alg 3: alpha,beta, def. alpha=0.25,beta=0.5 */  
/*                  prot=...,   protocol file                               */
/*                  fmt=...,    print format, def. 10.4                     */
/*                  nfmt=...,   integer print format, def. 4                */
/*              ) = fname;                                                  */
/*                                                                          */
/*              The original and the permuted distance matrices are         */
/*              written into the output file fname.                         */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int uds(TDAContext *ctx)
{
    int err,n;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Uni-dimensional scaling. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto UDSFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto UDSFin;
    if (gdd_tcheck(ctx, 0,1,2))          /* need undirected valued graph */
        goto UDSFin;

    if (ctx->GD_NP < 3) {
        printf1(ctx, "Error: need at least three nodes.\n");
        goto UDSFin;
    }
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    n = ctx->GD_NP;     /* number of points */
    printf1(ctx, "Number of points: %d\n",n);              

    if (ctx->PMALG < 1 || ctx->PMALG > 3)
        ctx->PMALG = 1;

    if (ctx->PMALG == 1) {
        printf1(ctx, "\nAlgorithm 1: Branch and bound (Defays 1978).\n");
        if (n >= 100) {
            printf1(ctx, "Current maximum: 100 points.\n");
            goto UDSFin;
        }
        uds_1(ctx);   
    } 
    else if (ctx->PMALG == 2) {
        printf1(ctx, "\nAlgorithm 2: approx. quadratic assignment.\n");
        if (uds_qa1(ctx))  
            goto UDSFin;

    }
    else {
        printf1(ctx, "\nAlgorithm 3: approx. quadratic assignment.\n");

        if (ctx->MxIter < 1 || ctx->MxItFlg == 0)
            ctx->MxIter = 100;
        if (ctx->PMIDFA <= 0.0 || ctx->PMIDFA > 1.0)
            ctx->PMIDFA = 0.25;
        if (ctx->PMIDFB <= 0.0 || ctx->PMIDFB > 1.0)
            ctx->PMIDFB = 0.5;
        printf1(ctx, "Max number of iterations: %d (alpha %lg, beta %lg).\n\n",
            ctx->MxIter,ctx->PMIDFA,ctx->PMIDFB);

        if (uds_qa2(ctx))  
            goto UDSFin;
    }
    err = 0;

UDSFin:       
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  uds_check(n,a,b)                                                        */
/*                                                                          */
/*  Check whether distances in a and b have same order.                     */

void uds_check(TDAContext *ctx, int n,double *a,double *b)
{
    register int i,j,k,l;
    int d0,d1;
    double aij,bij,akl,bkl;

    d0 = d1 = 0;
    for (i = 1; i < n; ++i) {
        for (j = 0; j < i; ++j) {
            aij = a[i * n + j + 1];
            bij = b[i * n + j + 1];

            for (k = i; k < n; ++k) {
                if (k == i)
                    l = j + 1;
                else
                    l = 0;
                for ( ; l < k; ++l) {
                    d0++;
                    akl = a[k * n + l + 1];
                    bkl = b[k * n + l + 1];
                          
                    if (aij < akl) {
                        if (bij > bkl)
                            d1++;
                    }
                    else if (aij > akl) {
                        if (bij < bkl)
                            d1++;
                    }
                    /*******
                    tda_out("i=%3d j=%3d k=%d l=%d aij=%6.1lf akl=%6.1lf bij=%6.1lf bkl=%6.1lf  %6d %6d\n",
                    i,j,k,l,aij,akl,bij,bkl,d0,d1);
                    ***************************/
                }
            }
        }
    }
    printf1(ctx, "%d of %d compares are inconsistent.\n",d1,d0);

}

/* -##--------------------------------------------------------------------- */
/*  uds_ncheck(n,a,b)                                                       */
/*                                                                          */
/*  Check whether distances in a and b have same order.                     */

void uds_ncheck(TDAContext *ctx, int n,int *a,int *b)
{
    register int i,j,k,l;
    int d0,d1,aij,bij,akl,bkl;

    d0 = d1 = 0;
    for (i = 1; i < n; ++i) {
        for (j = 0; j < i; ++j) {
            aij = a[i * n + j + 1];
            bij = b[i * n + j + 1];

            for (k = i; k < n; ++k) {
                if (k == i)
                    l = j + 1;
                else
                    l = 0;
                for ( ; l < k; ++l) {
                    d0++;
                    akl = a[k * n + l + 1];
                    bkl = b[k * n + l + 1];
                          
                    if (aij < akl) {
                        if (bij > bkl)
                            d1++;
                    }
                    else if (aij > akl) {
                        if (bij < bkl)
                            d1++;
                    }
                    /*******
                    tda_out("i=%3d j=%3d k=%d l=%d aij=%6.1lf akl=%6.1lf bij=%6.1lf bkl=%6.1lf  %6d %6d\n",
                    i,j,k,l,aij,akl,bij,bkl,d0,d1);
                    ***************************/
                }
            }
        }
    }
    printf1(ctx, "%d of %d compares are inconsistent.\n",d1,d0);

}

/* -##--------------------------------------------------------------------- */
/*  uds_1()     uni-dimensional scaling with branch and bound.              */
/*                                                                          */
/*  return  0 if successful,                                                */
/*         -1 if insuff. memory.                                            */
/*         -2 if undef. elements in distance matrix.                        */


int uds_1(TDAContext *ctx)   
{
    register int i,j;
    int n;
    double d,f,tmp;

    n = ctx->GD_NP;     /* number of points */

    if (alloc_aci(ctx, n + 1))          /* inverse permutation */
        return(-1);       

    fprintf(ctx->PMFd,"Original distance matrix\n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            ctx->UDS_D[i * n + j + 1] = gdd_adj(ctx, i,j,ctx->PMGN);
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->UDS_D[i * n + j + 1]);           
        }
        fprintf(ctx->PMFd,"\n");           
    }
    ctx->UDSDPP = 0.0;
    for (i = 0; i < n; ++i) {
        d = 0.0;
        for (j = 0; j < n; ++j)
            d += gdd_adj(ctx, i,j,ctx->PMGN);
        if (ctx->UDSDPP < d)
            ctx->UDSDPP = d;
    }
    printf1(ctx, "UDSDPP = %g\n",ctx->UDSDPP);
    uds_0(ctx, n,1);    /* find first feasible permutation */
    uds_r(ctx, n,0); 

    for (i = 1; i <= n; ++i) {
        d = uds_get_dip(ctx, n,i,ctx->UDS_E0) - uds_get_dpi(ctx, n,i,ctx->UDS_E0);
        d /= (double)n;
        ctx->UDS_X[i] = d;                  
        ctx->AcI[ctx->UDS_E0[i]] = i;     /* inverse permutation */
    }
    /***********************
    tda_out("UDS_X: ");
    for (i=1; i <= n; ++i)
        tda_out("%6.2lf ",UDS_X[i]);
    newline(ctx);
    tda_out("UDS_E0: ");
    for (i=1; i <= n; ++i)
        tda_out("%d ",UDS_E0[i]);
    newline(ctx);
    tda_out("ACI   : ");
    for (i=1; i <= n; ++i)
        tda_out("%d ",AcI[i]);
    newline(ctx);
    ***************/

    printf1(ctx, "\nBest solution: %lg\n",ctx->UDS_F0);

    printf1(ctx, "\nIndex  X-value     Object  Permutation\n");
    for (i = 1; i <= n; ++i) {
        printf1(ctx, "%4d  ",i);
        rt_printf1_d(ctx, ctx->PMFmtS,ctx->UDS_X[i]);
        printf1(ctx, "%7d   %8d\n",ctx->UDS_E0[i],ctx->AcI[i]);
    }
    f = 0.0;
    fprintf(ctx->PMFd,"\nFitted distance matrix\n");
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) {
            d = fabs(ctx->UDS_X[ctx->AcI[i]] - ctx->UDS_X[ctx->AcI[j]]);
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,d);
            if (j != i) {
                tmp = ctx->UDS_D[(i - 1) * n + j] - d;
                f += tmp * tmp;
            }
        }
        fprintf(ctx->PMFd,"\n");
    }
    if (f > 0.0)
        f = sqrt(f);

    printf1(ctx, "\nDistance between original and fitted matrix: %lg\n",f);

    f = 0.0;
    for (i = 2; i <= n; ++i) {
        for (j = 1; j < i; ++j) {
            d = fabs(ctx->UDS_X[ctx->AcI[i]] - ctx->UDS_X[ctx->AcI[j]]);
            f += fabs(ctx->UDS_D[(i - 1) * n + j] - d);
        }
    }
    f /= (double)(n * (n - 1) / 2);
    printf1(ctx, "Mean difference between original and fitted distances: %lg\n",f);
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  uds_0(n,k)    find first feasible permutation.                          */

int uds_0(TDAContext *ctx, int n,int k)
{
    (void)k;        /* unused: the signature is shared */
    register int r,i,j,jj;
    int imax = 0,first;
    double b,bmax;   

    for (i = 1; i <= n; ++i) {
        ctx->UDS_E0[i] = ctx->UDS_PI[i] = i;
        ctx->UDS_IF[i] = 0;
    }
    for (r = 1; r <= n; ++r) {
        for (i = 1; i <= n; ++i)  
            ctx->UDS_PI[i] = ctx->UDS_E0[i];
        first = 1;
        for (j = 1; j <= n; ++j) {
            if (ctx->UDS_IF[j])
                continue;
            jj = 0;
            for (i = r; i <= n; ++i) {
                if (ctx->UDS_PI[i] == j) {
                    jj = i;
                    break;
                }
            }   
            ctx->UDS_PI[jj] = ctx->UDS_PI[r];
            ctx->UDS_PI[r] = j;

            b = 0.0;
            for (i = 1; i <= r; ++i) 
                b += (uds_get_dip(ctx, n,i,ctx->UDS_PI) - uds_get_dpi(ctx, n,i,ctx->UDS_PI));

            if (first || b > bmax) {
                bmax = b;
                imax = j;
                for (i = r; i <= n; ++i)
                    ctx->UDS_E0[i] = ctx->UDS_PI[i];
                first = 0;
            }
        }
        ctx->UDS_IF[imax] = 1;
    }
    ctx->UDS_F0 = uds_get_f(ctx, n,ctx->UDS_E0);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  uds_get_dip(n,i,p)    return d(i,.) for permutation p.                  */

double uds_get_dip(TDAContext *ctx, int n,int i,int *p)
{
    register int j,k;
    double d;      

    d = 0.0;
    k = (p[i] - 1) * n;
    for (j = i + 1; j <= n; ++j)     
        d += ctx->UDS_D[k + p[j]];
    return(d);
}

/* ------------------------------------------------------------------------ */
/*  uds_get_dpi(n,i,p)    return d(.,i) for permutation p.                  */

double uds_get_dpi(TDAContext *ctx, int n,int i,int *p)
{
    register int j,k;
    double d;      

    d = 0.0;
    k = (p[i] - 1) * n;
    for (j = 1; j < i; ++j)  
        d += ctx->UDS_D[k + p[j]];
    return(d);
}

/* ------------------------------------------------------------------------ */
/*  uds_get_f(n,p)    return function value for permutation p.              */

double uds_get_f(TDAContext *ctx, int n,int *p)
{
    register int i;   
    double d,tmp;  

    d = 0.0;
    for (i = 1; i <= n; ++i) {
        tmp =  uds_get_dip(ctx, n,i,p) - uds_get_dpi(ctx, n,i,p);
        d += tmp * tmp;
    }
    return(d);
}

/* ------------------------------------------------------------------------ */
/*  uds_r(n,k)      recursively check all permutations.                     */

int uds_r(TDAContext *ctx, int n,int k)
{
    register int i,j,jj;
    double b,ck,cl,ckk0;

    printf1(ctx, "Enter uds_r n=%d k=%d E0: ",n,k);
    for (i = 1; i <= n; ++i)
        printf1(ctx, "%d ",ctx->UDS_E0[i]);
    printf1(ctx, "   F0 = %g\n",ctx->UDS_F0);

    for (i = 1; i <= n; ++i) {
        ctx->UDS_PI[i] = i;
        ctx->UDS_IF[i] = 0;
    }
    for (j = 1; j <= n; ++j) {

        if (ctx->UDS_PI[1] != j) {
            jj = 0;
            for (i = 2; i <= n; ++i) {
                if (ctx->UDS_PI[i] == j) {
                    jj = i;
                    break;
                }
            }   
            i = ctx->UDS_PI[1];
            ctx->UDS_IF[i] = 0;
            ctx->UDS_PI[jj] = i;
            ctx->UDS_PI[1] = j;
        }
        ctx->UDS_IF[j] = 1;

        ckk0 = uds_get_dip(ctx, n,1,ctx->UDS_PI) - uds_get_dpi(ctx, n,1,ctx->UDS_PI);
        ck = ckk0 * ckk0;
        b = ck + (double)(n - 1) * ctx->UDSDPP * ctx->UDSDPP;

        if (b >= ctx->UDS_F0) {
            cl = 0.0;
            uds_r0(ctx, n,1,ck,ckk0,cl);        /* branch to left */
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  uds_r0(n,k,ck,ckk0,cl)                                                  */
/*                                                                          */
/*  branch to node R(r(1),...,r(k) | r(n-k+1),...,r(n))                     */

int uds_r0(TDAContext *ctx, int n,int k,double ck,double ckk0,double cl)
{
    register int i,j,l,jj,ll;
    int lk1;
    double b,ckk,cll,cll0,cll1,tmp;

    ckk = ckk0 * ckk0;
    l = n - k + 1;
    lk1 = l - k - 1;

    for (i = 1; i <= n; ++i) {
        if (k == 1 && i <= k)
            continue;

        if (ctx->UDS_IF[i])
            continue;

        if (ctx->UDS_PI[l] != i) {
            jj = 0;
            for (j = k + 1; j < l; ++j) {
                if (ctx->UDS_PI[j] == i) {
                    jj = j;
                    break;
                }
            }   
            ll = ctx->UDS_PI[l];
            ctx->UDS_IF[ll] = 0;
            ctx->UDS_PI[jj] = ll;
            ctx->UDS_PI[l] = i;
        }
        ctx->UDS_IF[i] = 1;
        cll0 = uds_get_dip(ctx, n,l,ctx->UDS_PI) - uds_get_dpi(ctx, n,l,ctx->UDS_PI);
        cll  = cll0 * cll0;
                                      
        if (cll0 - ctx->EPSI1 > ckk0) {          /* not feasible */ 
            ctx->UDS_IF[i] = 0;
            continue;
        }
        
        if (l < n) {                /* second check */
         
            cll1 = uds_get_dip(ctx, n,l+1,ctx->UDS_PI) - uds_get_dpi(ctx, n,l+1,ctx->UDS_PI);
            if (cll0 < cll1 + 2.0 * ctx->UDS_D[(ctx->UDS_PI[l] - 1) * n + ctx->UDS_PI[l+1]]) {
                ctx->UDS_IF[i] = 0;
                continue;
            }
        }
        if (lk1 > 1) {    
            b = ck + cl + cll + (double)lk1 * dmax(ctx, ckk,cll);
            if (b + ctx->EPSI1 < ctx->UDS_F0) {
                ctx->UDS_IF[i] = 0;
                continue;
            }
            uds_r1(ctx, n,k,l,ck,cl + cll,cll0);     /* branch to right */
        }
        else {
            tmp = uds_get_dip(ctx, n,k+1,ctx->UDS_PI) - uds_get_dpi(ctx, n,k+1,ctx->UDS_PI);
            b = ck + cl + cll + tmp * tmp;
            if (b - ctx->EPSI1 > ctx->UDS_F0) {
                ctx->UDS_F0 = b;
                for (j = 1; j <=n; ++j)
                    ctx->UDS_E0[j] = ctx->UDS_PI[j];
            }
        }
        ctx->UDS_IF[i] = 0;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  uds_r1(n,k,l,ck,cl,cll0)                                                */
/*                                                                          */
/*  branch to node R(r(1),...,r(k+1) | r(l),...,r(n))                       */

int uds_r1(TDAContext *ctx, int n,int k,int l,double ck,double cl,double cll0)
{
    register int i,j,jj,kk,k1;
    int lk2;
    double b,ckk,ckk0,ckk1,cll,tmp;
  
    k1 = k + 1;
    lk2 = l - k1 - 1;
    cll = cll0 * cll0;

    for (i = 1; i <= n; ++i) {
        if (ctx->UDS_IF[i])
            continue;

        if (ctx->UDS_PI[k1] != i) {
            jj = 0;
            for (j = k1 + 1; j < l; ++j) {
                if (ctx->UDS_PI[j] == i) {
                    jj = j;
                    break;
                }
            }   
            kk = ctx->UDS_PI[k1];
            ctx->UDS_IF[kk] = 0;
            ctx->UDS_PI[jj] = kk;
            ctx->UDS_PI[k1] = i;
        }
        ctx->UDS_IF[i] = 1;
 
        ckk0 = uds_get_dip(ctx, n,k1,ctx->UDS_PI) - uds_get_dpi(ctx, n,k1,ctx->UDS_PI);
       
        if (cll0 - ctx->EPSI1 > ckk0) {                  /* not feasible */
            ctx->UDS_IF[i] = 0;
            continue;
        }

                /* second check */
          
        ckk1 = uds_get_dip(ctx, n,k,ctx->UDS_PI) - uds_get_dpi(ctx, n,k,ctx->UDS_PI);
        if (ckk1 - 2.0 * ctx->UDS_D[(ctx->UDS_PI[k] - 1) * n + ctx->UDS_PI[k1]] + ctx->EPSI1 < ckk0) {
            ctx->UDS_IF[i] = 0;
            continue;
        }
       
        ckk = ckk0 * ckk0;  

        if (lk2 > 1) {
            b = ck + ckk + cl + (double)lk2 * dmax(ctx, ckk,cll);
            if (b + ctx->EPSI1 < ctx->UDS_F0) {
                ctx->UDS_IF[i] = 0;
                continue;
            }
            uds_r0(ctx, n,k1,ck + ckk,ckk0,cl);        /* branch to left */
        }
        else {
            tmp = uds_get_dip(ctx, n,k1+1,ctx->UDS_PI) - uds_get_dpi(ctx, n,k1+1,ctx->UDS_PI);
            b = ck + ckk + cl + tmp * tmp;
            if (b - ctx->EPSI1 > ctx->UDS_F0) {
                ctx->UDS_F0 = b;
                for (j = 1; j <=n; ++j)
                    ctx->UDS_E0[j] = ctx->UDS_PI[j];
            }
        }
        ctx->UDS_IF[i] = 0;
    }
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  uds_qa1()   Uni-dimensional scaling with approx. quadratic assignment.  */
/*              (uses                                                       */
/*                                                                          */
/*  return  0 if successful,                                                */
/*         -1 if insuff. memory.                                            */
/*         -2 if undef. elements in distance matrix.                        */
/*         -3 if error in qap_w().                                          */

int uds_qa1(TDAContext *ctx)   
{
    register int i,j,k;
    int n,nn,err;
    double tmp;

    n = ctx->GD_NP;
    nn = n * n;

    if (alloc_acu(ctx, nn + 1))          /* C matrix */
        return(-1);       
    if (alloc_acv(ctx, nn + 1))          /* F matrix */
        return(-1);   
    if (alloc_acz(ctx, nn + 1))          /* D matrix */
        return(-1);   
    if (alloc_acn(ctx, 3 * n + 1))       /* loc3n */
        return(-1);   
    if (alloc_actmp(ctx, 4 * n + 1))       /* wtmp */
        return(-1);   
    
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            tmp = gdd_adj(ctx, i,j,ctx->PMGN);
            if (tmp < 0.0) {
                tmp = 0.0;
                if (i != j) {
                    printf1(ctx, "Error: undefined distance (%d,%d).\n",i+1,j+1);
                    return(-2);
                }
            }
            ctx->AcZ[i * n + j + 1] = tmp;
        }
    }

    fprintf(ctx->PMFd,"Original distance matrix\n");
    for (i = 0; i < n; ++i) {
        for (j = 1; j <= n; ++j)
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcZ[i * n + j]);
        fprintf(ctx->PMFd,"\n");
    }

    for (i = 0; i < n; ++i) {
        k = i;
        for (j = 0; j <= i; ++j) {
            ctx->AcV[i * n + j + 1] = (double)(-k);
            k--;
        }
        k = 1;
        for (j = i + 1; j < n; ++j) {
            ctx->AcV[i * n + j + 1] = (double)(-k);
            k++;
        }
    }

    /******* 
    printf1(ctx, "D\n");
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) 
            printf1(ctx, "%6.3lf ",AcZ[(i - 1) * n + j]);
        newline(ctx);
    }

    printf1(ctx, "F\n");
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) 
            printf1(ctx, "%6.3lf ",AcV[(i - 1) * n + j]);
        newline(ctx);
    }
    ****/

    err = qap_w(ctx, n,ctx->AcU,ctx->AcV,ctx->AcZ,ctx->AcN,ctx->AcTmp);

    if (err) {
        printf1(ctx, "Error in function qap_w.\n");
        return(-3);
    }         

    /* make inverse permutation */

    for (i = 1; i <= n; ++i)
        ctx->AcN[ctx->AcN[i] + n] = i;

    printf1(ctx, "Best function value: %lg\n",-ctx->AcTmp[1]);
    printf1(ctx, "Optimal ordering\n");
    for (i = 1; i <= n; ++i)
        printf1(ctx, "%4d %4d\n",i,ctx->AcN[i]);
    newline(ctx);
    printf1(ctx, "Permutation\n");
    for (i = 1; i <= n; ++i)
        printf1(ctx, "%4d %4d\n",i,ctx->AcN[i + n]);
    newline(ctx);

    fprintf(ctx->PMFd,"\nFitted distance matrix\n");
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) {
            ctx->AcU[(i - 1) * n + j] = tmp = (double)iabs(ctx, ctx->AcN[i + n] - ctx->AcN[j + n]);
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
        }
        fprintf(ctx->PMFd,"\n");
    }

    /* check numer of inconstencies */
          
    uds_check(ctx, n,ctx->AcZ,ctx->AcU);
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  uds_qa2()   uni-dimensional scaling with approx. quadratic assignment.  */
/*                                                                          */
/*  return  0 if successful,                                                */
/*         -1 if insuff. memory.                                            */
/*         -2 if undef. elements in distance matrix.                        */
/*                                                                          */
/*  NOTE: This functions presupposes an intervalued distance matrix!        */
/*                                                                          */

int uds_qa2(TDAContext *ctx)   
{
    register int i,j,k,l;
    int n,nn,n1,iter,ito,seed,bestv,look4;
    int *a,*b,*srtf,*srtif,*srtd,*srtid,*srtc,*srtic,*indexf,*indexd,*cost,*fdind;
    double tmp;

    n = ctx->GD_NP;
    nn = n * n;

    if (alloc_acn(ctx, nn + 1))          /* f matrix */
        return(-1);       
    if (alloc_acm(ctx, nn + 1))          /* d matrix */
        return(-1);   
    
    n1 = k = l = 0;
    for (i = 0; i < n; ++i) {
        l = i;
        for (j = 0; j < n; ++j) {
            tmp = gdd_adj(ctx, i,j,ctx->PMGN);
            if (tmp < 0.0) {
                tmp = 0.0;
                if (i != j) {
                    printf1(ctx, "Error: undefined distance (%d,%d).\n",i+1,j+1);
                    return(-2);
                }
            }
            ctx->AcN[++k] = (int)tmp;
            ctx->AcM[k] = -l;
            if (j >= i)
                l++;
            else
                l--;
        }
    }

    fprintf(ctx->PMFd,"Original distance matrix\n");
    for (i = 0; i < n; ++i) {
        for (j = 1; j <= n; ++j)
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcN[i * n + j]);
        fprintf(ctx->PMFd,"\n");
    }


    /******* 
    printf1(ctx, "F\n");
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) 
            printf1(ctx, "%2d ",AcN[(i - 1) * n + j]);
        newline(ctx);
    }

    printf1(ctx, "D\n");
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) 
            printf1(ctx, "%2d ",AcM[(i - 1) * n + j]);
        newline(ctx);
    }
    ****/
              
    if (alloc_acs(ctx, n + 1))           /* optimal permutation vector */
        return(-1);      

    if (alloc_acr(ctx, 10 * nn + 2 * n + 2))     /* working storage */
        return(-1);            

    a = ctx->AcR;
    b = a + n;
    srtf = b + n;
    srtif = srtf + nn;
    srtd  = srtif + nn;
    srtid = srtd + nn;
    srtc  = srtid + nn;
    srtic = srtc + nn;
    indexf = srtic + nn;
    indexd = indexf + nn;
    cost = indexd + nn;
    fdind = cost + nn;
    seed = 270001;
    look4 = -1;
   
    gqapd(ctx, n,ctx->MxIter,ctx->PMIDFA,ctx->PMIDFB,look4,&seed,ctx->AcN,ctx->AcM,a,b,srtf,srtif,srtd,srtid,
          srtc,srtic,indexd,indexf,cost,fdind,ctx->AcS,&bestv,&iter,&ito);

    printf1(ctx, "Performed %d iterations.\n",iter);
    printf1(ctx, "Best cost value: %d  (found in %d iterations).\n",-bestv,ito);
    n1 = 0;
    for (i = 1; i <= n; ++i) {
        k = ctx->AcS[i];
        for (j = 1; j <= n; ++j) {
            l = ctx->AcS[j];
            n1 += ctx->AcM[(i - 1) * n + j] * ctx->AcN[(k - 1) * n + l];
        }
    }
    printf1(ctx, "Checked cost value: %d\n",-n1);
        
    /* make inverse permutation */

    for (i = 1; i <= n; ++i)
        ctx->AcR[ctx->AcS[i]] = i;

    printf1(ctx, "Optimal ordering\n");
    for (i = 1; i <= n; ++i)
        printf1(ctx, "%4d %4d\n",i,ctx->AcS[i]);
    newline(ctx);
    printf1(ctx, "Permutation\n");
    for (i = 1; i <= n; ++i)
        printf1(ctx, "%4d %4d\n",i,ctx->AcR[i]);
    newline(ctx);

    fprintf(ctx->PMFd,"\nFitted distance matrix\n");
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) {
            ctx->AcM[(i - 1) * n + j] = k = iabs(ctx, ctx->AcR[i] - ctx->AcR[j]);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,k);
        }
        fprintf(ctx->PMFd,"\n");
    }

    /* check numer of inconstencies */
          
    uds_ncheck(ctx, n,ctx->AcN,ctx->AcM);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  sga         Scalogram analysis.                                         */
/*                                                                          */
/*              sga(                                                        */
/*                                                                          */
/*              ) = X1,...,Xm;      variables.                              */
/*                                                                          */
/*  Variables are translated into binary indicators:                        */
/*  X(i,j) >  0  ->  1  (individual i solves task j successfully)           */
/*         <= 0  ->  0  otherwise                                           */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int sga(TDAContext *ctx)
{
    register int i,j,k;
    int err,n,m,s,nn,nf,len,err0,err1,err2;
    double d;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Scalogram analysis. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 3,4,1))       /* get parameters */
        goto SGAFin;

    n = ctx->NOC;
    m = ctx->PMNV;                       /* number of variables */

    if (alloc_acn(ctx, m * n + 1))
        goto SGAFin;
    if (alloc_acs(ctx, n + 1))
        goto SGAFin;
    if (alloc_aci(ctx, n + 1))
        goto SGAFin;
    if (alloc_acr(ctx, m + 1))
        goto SGAFin;
    if (alloc_ack(ctx, m + 1))
        goto SGAFin;
    if (alloc_acm(ctx, n + 1))
        goto SGAFin;

    for (i = 0; i < ctx->NOC; ++i) {                 /* number of cases */
        for (j = 0; j < m; ++j) {
            d = get_data(ctx, ctx->PMVIdx[j],i);
            if (d > 0.0) {
                ctx->AcS[i] += 1;
                ctx->AcR[j] += 1;
            }
        }
    }
    if (sortdpi(ctx, n,ctx->AcS,ctx->AcI))     /* sort scores */
        goto SGAFin;

    if (sortdpi(ctx, m,ctx->AcR,ctx->AcK))     /* sort difficulties */
        goto SGAFin;

    /* put reduced and sorted data matrix into AcN */

    k = ctx->AcI[0];
    for (j = 0; j < m; ++j) {
        d = get_data(ctx, ctx->PMVIdx[ctx->AcK[j]],k);
        if (d > 0.0)
            ctx->AcN[j] = 1;
    }
    ctx->AcM[0] = 1;
    nn = 1;
    for (i = 1; i < n; ++i) {
        k = ctx->AcI[i];
        nf = 1;
        for (j = 0; j < m; ++j) {
            if (get_data(ctx, ctx->PMVIdx[ctx->AcK[j]],k) > 0.0)
                s = 1;
            else
                s = 0;
            ctx->AcN[nn * m + j] = s;
            if (s != ctx->AcN[(nn - 1) * m + j])
                nf = 0;
        }
        if (nf == 1)
            ctx->AcM[nn - 1] += 1;
        else {
            ctx->AcM[nn] = 1;
            nn++;
        }
    }
    if (alloc_acr(ctx, m))
        goto SGAFin;
    if (alloc_acs(ctx, nn))
        goto SGAFin;

    for (i = 0; i < nn; ++i) {
        for (j = 0; j < m; ++j) {
            if (ctx->AcN[i * m + j] == 1) {
                ctx->AcS[i] += 1;
                ctx->AcR[j] += ctx->AcM[i];
            }
        }
    }
    if (alloc_aci(ctx, nn))      /* error1 (Goodenough) */
        goto SGAFin;
    if (alloc_acj(ctx, nn))      /* error2 (Guttman) */
        goto SGAFin;

    len = 4;
    for (j = 0; j < m; ++j)
        len = (int)(imax(ctx, len,(int)(strlen(ctx->VName[ctx->PMVIdx[j]]))));

    printf1(ctx, "\n Idx ");
    for (j = 0; j < m; ++j) {
        k = ctx->AcK[j];
        prnchar(ctx, ' ',len - (int)strlen(ctx->VName[ctx->PMVIdx[k]]),0);
        printf1(ctx, "%s ",ctx->VName[ctx->PMVIdx[k]]);
    }
    printf1(ctx, "   Freq   Value   Err0   Err1   Err2\n");
    prnchar(ctx, '-',41 + m * (len + 1),1);

    err0 = err1 = err2 = 0;
    for (i = 0; i < nn; ++i) {
        printf1(ctx, "%4d ",i + 1);

        s = m - ctx->AcS[i];
        for (j = 0; j < m; ++j) {
            k = ctx->AcN[i * m + j];                         
            prnchar(ctx, ' ',len - 4,0);
            printf1(ctx, "%4d ",k);                         
            
            if (k == 1) {
                if (j < s)
                    ctx->AcI[i] += 1;
            }
            else {
                if (j >= s)
                    ctx->AcI[i] += 1;
            }
        }
        s = imin(ctx, ctx->AcI[i],1);
        err0 += s * ctx->AcM[i];

        err1 += ctx->AcI[i] * ctx->AcM[i];

        ctx->AcJ[i] = sga1(ctx, m,ctx->AcN + i * m);
        err2 += ctx->AcJ[i] * ctx->AcM[i];

        printf1(ctx, "%7d %7d %6d %6d %6d\n",ctx->AcM[i],ctx->AcS[i],s,ctx->AcI[i],ctx->AcJ[i]);                  
    }
    prnchar(ctx, '-',41 + m * (len + 1),1);
    prnchar(ctx, ' ',5,0);

    for (j = 0; j < m; ++j) {
        prnchar(ctx, ' ',len - 4,0);
        printf1(ctx, "%4d ",ctx->AcR[j]);
    }
    printf1(ctx, "%7d     %10d %6d %6d\n\n",n,err0,err1,err2);

    printf1(ctx, "Reproducibility\n");
    printf1(ctx, "Err0: %6.4f\n",1.0 - (double)err0 / (double)n);
    printf1(ctx, "Err1: %6.4f\n",1.0 - (double)err1 / (double)(n * m));
    printf1(ctx, "Err2: %6.4f\n\n",1.0 - (double)err2 / (double)(n * m));


    hpo(ctx, nn,m,ctx->AcN);     

    err = 0;

SGAFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sga1(m,r)   given profile (r1,...,rm) return minimal distance to        */
/*              perfect profile (Guttman method).                           */

int sga1(TDAContext *ctx, int m,int *r)
{
    register int i,j,mm,e;
    int err;

    err = ctx->INTMAX;
    for (i = 0; i <= m; ++i) {
        mm = m - i;
        e = 0;
        for (j = 0; j < m; ++j) {
            if (j < mm) {
                if (r[j] != 0)
                    e++;
            }
            else {
                if (r[j] == 0)
                    e++;
            }
        }   
        err = imin(ctx, err,e);
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  hpo(n,m,x)                                                              */

int hpo(TDAContext *ctx, int n,int m,int *x)      
{
    register int i,j,k,im;
    int err,d;

    err = -1;
    if (alloc_acd(ctx, n * n))   
        goto HPO1Fin;

    for (i = 0; i < n; ++i) {
        im = i * n;
        for (j = 0; j < n; ++j) {
            if (j == i)
                continue;

            d = hpo1(ctx, m,x + i * m,x + j * m);
            if (d != 1)
                continue;
        
            for (k = 0; k < n; ++k) {
                if (k == i || k == j) {
                    continue;
                }
                if (hpo1(ctx, m,x + i * m,x + k * m) == 1 &&
                    hpo1(ctx, m,x + k * m,x + j * m) == 1) {
                    d = 0;
                    break;
                }
       
            }
            if (d == 1)
                ctx->AcD[im + j] = (char)d;
        }
    }

    for (i = 0; i < n; ++i) {
        im = i * n;
        for (j = 0; j < n; ++j)  
            printf1(ctx, "%d ",(int)ctx->AcD[im + j]);
        newline(ctx);
    }

    printf1(ctx, "digraph g {\n");
    for (i = 0; i < n; ++i) {
        printf1(ctx, "n%d [label=%d];\n",i+1,i+1);
    }

    for (i = 0; i < n; ++i) {
        im = i * n;
        for (j = 0; j < n; ++j) {
            if (ctx->AcD[im + j])
                printf1(ctx, "n%d -> n%d;\n",i+1,j+1);
        }
    }



    err = 0;
HPO1Fin:
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  hpo1(m,x,y)     return 1 if x <= y, otherwise 0.                        */

int hpo1(TDAContext *ctx, int m,int *x,int *y)
{
    (void)ctx;        /* unused: the signature is shared */
    register int j;

    for (j = 0; j < m; ++j) {
        if (x[j] > y[j])
            return(0);
    }
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  unf         Unfolding                                                   */
/*                                                                          */
/*              unf(                                                        */
/*                  w=...,          case weight variable                    */
/*                  df=...,         create output file                      */
/*                  fmt=...,        print format, def. 10.4                 */
/*                  x=...,          rank order                              */
/*                                                                          */
/*              ) = X1,...,Xm;      variables.                              */
/*                                                                          */
/*  Values of the variables are interpreted as follows. If the values       */
/*  for case i are: (x1,...,xm), this means:                                */
/*  a)  alternative i is preferred to alternative j if xi > xj              */
/*  b)  indifference if xi = xj.                                            */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

double UNFMax = 0.0;

int unf(TDAContext *ctx)
{
    register int i,j,k;
    int err,n,m,first,wf;
    double wt,nn;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Unfolding. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 3,4,1))       /* get parameters */
        goto UNFFin;

    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    n = ctx->NOC;
    m = ctx->PMNV;                       /* number of variables */

    if (m < 2) {
        printf1(ctx, "Error: need at least two variables.\n");
        goto UNFFin;
    }
    if (alloc_acx(ctx, n * m + 1))           /* input data */
        goto UNFFin;

    wf = 0;
    if (ctx->PMWVar >= 0) {
        printf1(ctx, "Using case weights defined by: %s\n",ctx->VName[ctx->PMWVar]);
        if (alloc_acw(ctx, n + 1))   
            goto UNFFin;
        wf = 1;
    }
    if (alloc_acm(ctx, m + 1))               /* best permutation */
        goto UNFFin;

    for (i = 0; i < n; ++i) {           /* get data into AcX */
        for (j = 0; j < m; ++j)  
            ctx->AcX[i * m + j] = get_data(ctx, ctx->PMVIdx[j],i);

        if (ctx->PMWVar >= 0) {
            wt = get_data(ctx, ctx->PMWVar,i);
            if (wt < 0.0) {
                printf1(ctx, "Error: found negative weight in case %d\n",i + 1);
                goto UNFFin;
            }
            ctx->AcW[i] = wt;
        }
    }
                     
/***
    for (i = 0; i < n; ++i) {   
        for (j = 0; j < m; ++j)  
            printf1(ctx, "%g ",AcX[i * m + j]);
        newline(ctx);
    }
***/

    /*  do for all permutations of m alternatives */

    if (alloc_ack(ctx, m))
        goto UNFFin;
    if (alloc_aci(ctx, m))
        goto UNFFin;
    if (alloc_acj(ctx, m))
        goto UNFFin;
    if (alloc_acr(ctx, m))
        goto UNFFin;
    if (alloc_acs(ctx, m * m + 1))
        goto UNFFin;
    if (alloc_acn(ctx, m + 1))
        goto UNFFin;
    if (alloc_acns(ctx, m + 1))
        goto UNFFin;

    UNFMax = -(ctx->DBLMAX - 100);    
    nn = (double)(m * (m - 1) / 2);

    if (ctx->PMNTP > 0) {
        printf1(ctx, "\nCompare with rank order: %d",(int)ctx->PMTP[0]);
        for (j = 1; j < ctx->PMNTP; ++j)
            printf1(ctx, ", %g",ctx->PMTP[j]);
        newline(ctx);
        if (ctx->PMNTP != m) {
            printf1(ctx, "Error in number of alternatives.\n");
            goto UNFFin;
        }
        for (j = 0; j < m; ++j) {
            ctx->AcK[j] = k = (int)ctx->PMTP[j] - 1;
            if (k < 0 || k >= m) {
                printf1(ctx, "Error: invalid permutation.\n");
                goto UNFFin;
            }
        }
        unf1(ctx, n,m,ctx->AcK,ctx->AcR,ctx->AcX,ctx->AcS,ctx->AcN,0,ctx->AcNS,wf,ctx->AcW);  

        printf1(ctx, "Value (sum of taus): %g\n\n",UNFMax / nn);
        /***
        printf1(ctx, "Best permutation: ");
        for (j = 0; j < m; ++j)
            printf1(ctx, "%d ",AcM[j]);
        newline(ctx);
        **/
        err = 0;
        goto UNFFin;
    }











    first = 1;
    while (perm(ctx, m,ctx->AcK,ctx->AcI,ctx->AcJ,first)) {

        if (ctx->AcK[0] > ctx->AcK[m - 1])        /* half of permutations suffices */
            continue;

        unf1(ctx, n,m,ctx->AcK,ctx->AcR,ctx->AcX,ctx->AcS,ctx->AcN,0,ctx->AcNS,wf,ctx->AcW);  /* check sub-permutations */

printf1(ctx, "UNFMax=%g\n",UNFMax);

        first = 0;
    }
    printf1(ctx, "\nBest value (sum of taus): %g\n",UNFMax / nn);
    printf1(ctx, "Best permutation: ");
    for (j = 0; j < m; ++j)
        printf1(ctx, "%d ",ctx->AcM[j]);
    newline(ctx);

    if (ctx->PMF1Def) {
        for (j = 0; j < m; ++j)
            ctx->AcK[j] = ctx->AcM[j] - 1;

        unf1(ctx, n,m,ctx->AcK,ctx->AcR,ctx->AcX,ctx->AcS,ctx->AcN,1,ctx->AcNS,wf,ctx->AcW);

        printf1(ctx, "\n%d records written to: %s\n",n,ctx->PMF1dName);
    }
    err = 0;

UNFFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  unf1(n,m,p,r,x,s,rr,prn,pp,wf,w)                                        */
/*                                                                          */
/*                      check all sub-permutations of p[]. x[] contains     */
/*                      input data (rank orders), r[] is used to build the  */
/*                      sub-permutations. If a new maximum is found its     */
/*                      value is saved in UNFMax, the corresponding         */  
/*                      permutation is saved in AcM[].                      */
/*                                                                          */
/*                      If wf != 0 the vector w[] contains case weights.    */
/*                                                                          */
/*                      If prn != 0 print to output file.                   */

void unf1(TDAContext *ctx, int n,int m,int *p,int *r,double *x,int *s,int *rr,int prn, short *pp,int wf,double *w)
{
    register int i,j,k,l,t;  
    int ii,c,d;
    double xij,cmax,nn = 0.0;
        
    if (prn)
        nn = (double)(m * (m - 1) / 2);

    cmax = 0.0;                   /* max over all cases */
    for (ii = 0; ii < n; ++ii) {

        c = -(ctx->INTMAX - 100);    
                         
        l = ii * m;
        for (j = 0; j < m; ++j) {   
            r[j]= p[j] + 1;
            t = j * m;
            for (k = 0; k < j; ++k)  
                s[t + k] = (int)dsign(ctx, x[l + j] - x[l + k]);
        }
        d = unf2(ctx, m,r,s,rr);

        if (d > c) {
            c = d;
            if (prn) {
                for (j = 0; j < m; ++j)
                    pp[r[j] - 1] = (short)(m - j);              
            }
        }
    
        /* create all sub-permutations */
            
        for (i = 0; i < m - 1; ++i) {
            k = i % 2;
            l = k;
            while (l <= i) {
                t = r[l];
                r[l] = r[l + 1];
                r[l + 1] = t;
                l += 2;
            }
            d = unf2(ctx, m,r,s,rr);
            if (d > c) {
                c = d;
                if (prn) {
                    for (j = 0; j < m; ++j)
                        pp[r[j] - 1] = (short)(m - j);              
                }
            }
        }
        for (i = m - 3; i >= 0; --i) {
            k = i % 2;
            l = k;
            while (l <= i) {
                t = r[l];
                r[l] = r[l + 1];
                r[l + 1] = t;
                l += 2;
            }
            d = unf2(ctx, m,r,s,rr);
            if (d > c) {
                c = d;
                if (prn) {
                    for (j = 0; j < m; ++j)
                        pp[r[j] - 1] = (short)(m - j);              
                }
            }
        }
        if (wf)
            cmax += (double)c * w[ii];
        else
            cmax += (double)c;

        if (prn) {
            fprintf(ctx->PMF1d,"%5d  ",ii + 1);
            for (j = 0; j < m; ++j)
                fprintf(ctx->PMF1d,"%2d ",pp[j]);
            fprintf(ctx->PMF1d," ");
            t = 0;
            for (j = 0; j < m; ++j) {
                xij = x[ii * m + j];
                fprintf(ctx->PMF1d,"%2lg ",xij);
                if (xij != (double)pp[j])
                    t = 1;
            }
            fprintf(ctx->PMF1d,"  %d ",t);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,(double)c / nn);
            if (wf)
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,w[ii]);
            fprintf(ctx->PMF1d,"\n");
        }
    }
    if (prn)
        return;

    if (cmax > UNFMax) {        /* save new max and permutation */
        UNFMax = cmax;
        for (j = 0; j < m; ++j) 
            ctx->AcM[j] = p[j] + 1;
    }
}

/* ------------------------------------------------------------------------ */
/*  unf2(m,r,s)   return rank correlation (tau) between rr and x.           */
/*                r is set up by unf1. to get rr create rank orders:        */
/*                rr[r[j]-1] = m-j.                                         */
/*                                                                          */
/*                s[] is a m x m matrix containing sign(xi - xj).           */

int unf2(TDAContext *ctx, int m,int *r,int *s,int *rr)
{
    register int i,j;
    int t;
       
    for (j = 0; j < m; ++j)
        rr[r[j] - 1] = m - j;              

    t = 0;
    for (i = 1; i < m; ++i) {
        for (j = 0; j < i; ++j)  
            t += isign(ctx, rr[i] - rr[j]) * s[i * m + j];
    }

tda_out("rr: ");
for (j = 0; j < m; ++j)
tda_out("%d ",rr[j]);
newline(ctx);

         
    return(t);
}   

/* ------------------------------------------------------------------------ */
/*  mdsc        MDS with principal coordinates.                             */
/*              Requires undirected valued graph.                           */
/*                                                                          */
/*              mdsc(                                                       */
/*                  gn=...,     graph number, def. 1                        */
/*                  opt=        method if negative eigenvalues, def. 1      */
/*                               1 = do nothing                             */
/*                               2 = add constant to distance matrix        */
/*                  df=...,     additional output file                      */
/*                  fmt=...,    print format, def. 10.4                     */
/*              ) = fname;                                                  */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */
       
int mdsc(TDAContext *ctx)
{
    register int i,j,k;
    int err,n,r,nn,rflag,n1rec;
    double d,di,t,xi1,xi2,xj1,xj2,evmin,c;

    n1rec = rflag = 0;
    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "MDS with principal coordinates. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto MDSCFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto MDSCFin;
    if (gdd_tcheck(ctx, 0,1,2))          /* need undirected valued graph */
        goto MDSCFin;

    if (ctx->GD_NP < 3) {
        printf1(ctx, "Error: need at least three nodes.\n");
        goto MDSCFin;
    }
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    n = ctx->GD_NP;                      /* number of points */

    /* create A matrix in AcU */

    if (alloc_acu(ctx, n * n + 2))                                    
        goto MDSCFin; 
    if (alloc_acx(ctx, n + 1))                                    
        goto MDSCFin; 


MDSCRep:

    d = 0.0;
    for (i = 0; i < n; ++i) {
        di = 0.0;
        for (j = 0; j < n; ++j) {
            if (j != i) {
                t = gdd_adj(ctx, i,j,ctx->PMGN);
                if (rflag)
                    t += c;

                di += t * t;
                d  += t * t;
            }
        }
        ctx->AcX[i] = di;
    }   
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (i != j) {
                t = gdd_adj(ctx, i,j,ctx->PMGN);
                if (rflag)
                    t += c;
            }
            else 
                t = 0.0;

            ctx->AcU[i * n + j + 1] = -0.5 * (t * t - ctx->AcX[i] / (double)n
                                               - ctx->AcX[j] / (double)n
                                               + d / (double)(n * n));
        }
    }

    /* if requested write A-matrix */

    if (ctx->PMF1Def) {
        if (rflag == 0)
            fprintf(ctx->PMF1d,"A-matrix\n");
        else
            fprintf(ctx->PMF1d,"A-matrix (of modified distance matrix)\n");
        n1rec++;
        for (i = 0; i < n; ++i) {
            for (j = 1; j <= n; ++j)  
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcU[i * n + j]);
            fprintf(ctx->PMF1d,"\n");
            n1rec++;
        }
    }

    r = evecf(ctx, n,ctx->AcX,ctx->AcU);
    if (r) {
        if (r == -1) {
            printf1(ctx, "No success in eigenvalue calculation.\n");
            printf1(ctx, "Exceeded maximum number of iterations.\n");
        }
        else
            printf1(ctx, "Insufficient memory.\n");
        goto MDSCFin;
    }
            
    if (ctx->PMF1Def) {
        fprintf(ctx->PMF1d,"Eigenvalues\n");
        n1rec++;
        for (j = 1; j <= n; ++j)   
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[j]);
        fprintf(ctx->PMF1d,"\n");
        n1rec++;
    }
    printf1(ctx, "\nEigenvalue         per cent\n");
    r = 0;
    t = 0.0;
    for (j = 1; j <= n; ++j) {
        if (fabs(ctx->AcX[j]) <= ctx->EPSI1)
            ctx->AcX[j] = 0.0;
        if (ctx->AcX[j] < 0.0)  
            r++;
        t += ctx->AcX[j];
    }
    for (j = 1; j <= n; ++j) {
        prnchar(ctx, ' ',16 - ctx->PMFmt1,0);
        rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcX[j]);
        if (t > 0.0 && r == 0)
            printf1(ctx, "  %8.2f",100.0 * ctx->AcX[j] / t);
        newline(ctx);
#ifdef TDA_R_PACKAGE
        {
            double erow[2];
            erow[0] = ctx->AcX[j];
            erow[1] = (t > 0.0 && r == 0) ? 100.0 * ctx->AcX[j] / t : (double)(NAN);
            tda_export_row(ctx, "mds.eigenvalues", erow, 2);
        }
#endif
    }
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "mds.eigenvalues");
#endif
    newline(ctx);

    if (rflag == 0 && ctx->PMOPT == 2) {       /* find additive constant */
        evmin = ctx->AcX[n];
        c = 4.0 * evmin * evmin - 2.0 * evmin;
        if (c > 0.0)
            c = sqrt(c);
        else
            c = 0.0;
        c -= 2.0 * evmin;
        c += ctx->EPSI1;
        printf1(ctx, "Additive constant: %lg\n",c);
        if (c <= 0.0)
            goto MDSCFin;
        rflag = 1;
        goto MDSCRep;
    }
                 
    if (ctx->PMF1Def) {
        fprintf(ctx->PMF1d,"Eigenvectors\n");
        n1rec++;
        for (i = 0; i < n; ++i) {
            for (j = 1; j <= n; ++j)  
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcU[i * n + j]);
            fprintf(ctx->PMF1d,"\n");
            n1rec++;
        }
    }

    r = 0;
    for (j = 1; j <= n; ++j) {
        if (ctx->AcX[j] > 0.0)
            r++;
        else
            break;
    }
    if (r > 0) {
        for (i = 0; i < n; ++i) {
            for (j = 1; j <= r; ++j)  
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,sqrt(ctx->AcX[j]) * ctx->AcU[i * n + j]);
            fprintf(ctx->PMFd,"\n");
        }
        printf1(ctx, "%d records written to: %s\n",n,ctx->PMFdName);
    }
    if (n1rec > 0)
        printf1(ctx, "%d records written to: %s\n",n1rec,ctx->PMF1dName);


    /*******************************/
    if (ctx->PMTabFDef && r >= 2) {
        nn = n * (n - 1) / 2;
        if (alloc_acv(ctx, nn + 1))                                    
            goto MDSCFin; 
        if (alloc_acw(ctx, nn + 1))                                    
            goto MDSCFin; 
        if (alloc_ack(ctx, nn + 1))                                    
            goto MDSCFin; 

        k = 0;
        for (i = 0; i < n; ++i) {
            xi1 = sqrt(ctx->AcX[1]) * ctx->AcU[i * n + 1];
            xi2 = sqrt(ctx->AcX[2]) * ctx->AcU[i * n + 2];

            for (j = 0; j < i; ++j) {
                ctx->AcV[k] = gdd_adj(ctx, i,j,ctx->PMGN);
                xj1 = sqrt(ctx->AcX[1]) * ctx->AcU[j * n + 1];
                xj2 = sqrt(ctx->AcX[2]) * ctx->AcU[j * n + 2];

                ctx->AcW[k] = sqrt(pow(xi1 - xj1,2) + pow(xi2 - xj2,2));
                k++;
            }
        }
        if (sortdp2(ctx, k,ctx->AcV,ctx->AcW,ctx->AcK))
            goto MDSCFin;

        for (i = 0; i < k; ++i) {
            j = ctx->AcK[i];
            fprintf(ctx->PMTabFd,"%3d ",i);
            rt_fprintf_d(ctx, ctx->PMTabFd,ctx->PMFmtS,ctx->AcV[j]);
            rt_fprintf_d(ctx, ctx->PMTabFd,ctx->PMFmtS,ctx->AcW[j]);
            fprintf(ctx->PMTabFd,"\n");
        }
        printf1(ctx, "%d records written to: %s\n",k,ctx->PMTabFName);
    }
    /*****************/

    err = 0;

MDSCFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  mdsm        Metric MDS                                                  */
/*              Requires undirected valued graph.                           */
/*                                                                          */
/*              mdsm(                                                       */
/*                  gn=...,    graph number, def. 1                         */
/*                  mina=...,  minimization algorithm, def. 7               */
/*                             (possible are: 1, 2, 3, 4, 7)                */
/*                  mxit=...,  maximal number of iterations, default        */
/*                             depends on algorithm                         */
/*                  ns=...,    number of random restarts, def. 1            */
/*                  xp=...,    starting values                              */
/*                  dvs=...,   starting values                              */
/*                  df=...,    additional output file containing the        */
/*                             estimated distance matrix                    */
/*                  fmt=...,   print format, def. 10.4                      */
/*              ) = fname;                                                  */
/*                                                                          */
/*  Note: the iteration protocol can be suppressed with a separate          */
/*  silent=2 command.                                                       */
/*                                                                          */
/*  The procedure uses: MDSMOD4 mit mds4_fn().                              */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */
       
int mdsm(TDAContext *ctx)
{
    register int i,j,k;
    int err,nn,ii;
    double d,a,b,aa,bb,dm,fm = 0.0,ds;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Metric MDS. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto MDSMFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto MDSMFin;
    if (gdd_tcheck(ctx, 0,1,2))          /* need undirected valued graph */
        goto MDSMFin;

    if (ctx->GD_NP < 3) {
        printf1(ctx, "Error: need at least three nodes.\n");
        goto MDSMFin;
    }
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    ctx->MDSRN = ctx->GD_NP;                  /* number of points */

    /* setup MDSDIS vector with distances */

    nn = ctx->MDSRN * (ctx->MDSRN - 1) / 2;
    if (alloc_acv(ctx, nn + 1))                                        
        goto MDSMFin; 
    if (alloc_acw(ctx, nn + 1))                                        
        goto MDSMFin; 
    ctx->MDSDIS = ctx->AcV;
    ctx->MDSDIS1 = ctx->AcW;

    ds = dm = 0.0;
    k = 0;
    for (i = 2; i <= ctx->MDSRN; ++i) {
        for (j = 1; j < i; ++j) {
            d = gdd_adj(ctx, i-1,j-1,ctx->PMGN);
            ctx->MDSDIS[k++] = d;
            dm = dmax(ctx, dm,d);
            ds += d * d;
        }
    }

    /* prepare for min algorithm, always MINA = 7 */

    if (ctx->MINA == 5 || ctx->MINA == 6 || ctx->MINA == 8)
        ctx->MINA = 7;

    ctx->NParm = 2 * (ctx->MDSRN - 1);
    ctx->PMDOPT = 1;                     /* analytical gradient */
    ctx->CCTyp = 0;                      /* no covariance matrix */

    set_mlopt(ctx);                    /* adjust options */         

    if (ml_init(ctx, 1,1,0))             /* init function minimization */
        goto MDSMFin;

    if (get_dsv(ctx, ctx->NParm,ctx->Par,0,ctx->ParLB,ctx->ParUB,1))   /* try to get starting values */
        goto MDSMFin;
       
    if (ctx->DSVFlg)  
        ctx->PMNS = 1;
    else if (ctx->PMNS < 1)
        ctx->PMNS = 1;

    if (alloc_act(ctx, nn + 1))      /* used for best coordinates */   
        goto MDSMFin; 

    for (ii = 1; ii <= ctx->PMNS; ++ii) {    /* repeat */

        if (ctx->DSVFlg == 0) {              /* get initial configuration */
            for (i = 1; i <= ctx->NParm; ++i) 
                ctx->Par[i] = dm * random1(ctx);
        }
        if (ffmin(ctx, MDSMOD4,0,1,1))       /* minimization */
            goto MDSMFin;               /* insuff memory or fn error */

        prn_mlres(ctx, 1);                   /* print info about minimization */

        /* print result into output file */

        fprintf(ctx->PMFd,"%4d %3d %12.8lf  ",ii,ctx->LConv,ctx->FMin);
        for (i = 1; i <= ctx->NParm; ++i)
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->Par[i]);
        fprintf(ctx->PMFd,"\n");

        if (ii == 1) {
            fm = ctx->FMin;
            for (i = 1; i <= ctx->NParm; ++i)
                ctx->AcT[i] = ctx->Par[i];
        }
        else if (ctx->FMin < fm) {
            fm = ctx->FMin;     
            for (i = 1; i <= ctx->NParm; ++i)
                ctx->AcT[i] = ctx->Par[i];
        }
        if (ctx->DSVFlg)
            break;
    }
    printf1(ctx, "Best function value: %14.6lf (stress: %8.6lf)\n",fm,sqrt(fm/ds));
    printf1(ctx, "Coordinates\n");
    for (i = 1; i <= ctx->MDSRN; ++i) {
        printf1(ctx, "%3d ",i);
        if (i == 1)
            a = b = 0.0;
        else {
            a = ctx->AcT[i - 1];
            b = ctx->AcT[ctx->MDSRN + i - 2];
        }
        rt_printf1_d(ctx, ctx->PMFmtS,a);
        rt_printf1_d(ctx, ctx->PMFmtS,b);
        newline(ctx);
    }
    newline(ctx);

    if (ctx->PMF1Def) {      /* write estimated distance matrix */
        for (i = 1; i <= ctx->MDSRN; ++i) {
            if (i == 1)
                a = b = 0.0;
            else {
                a = ctx->AcT[i - 1];
                b = ctx->AcT[ctx->MDSRN + i - 2];
            }
            for (j = 1; j <= ctx->MDSRN; ++j) {
                if (j == 1)
                    aa = bb = 0.0;
                else {
                    aa = ctx->AcT[j - 1];
                    bb = ctx->AcT[ctx->MDSRN + j - 2];
                }
                d = mds_edis(ctx, a,b,aa,bb);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,d);
            }
            fprintf(ctx->PMF1d,"\n");
        }
        printf1(ctx, "%d records written to %s\n",ctx->MDSRN,ctx->PMF1dName);
    }
    printf1(ctx, "%d records written to %s\n",ctx->PMNS,ctx->PMFdName);

    err = 0;

MDSMFin:       
    ml_init(ctx, 0,0,0);      
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  mds4_fn()   Calculate function value for MDSMOD4.                       */
/*              Use parameter values in TPar[] (i=1,NParm).                 */
/*              Return function value in FTmp, gradient in Grad[].          */
/*                                                                          */
/*              Return 0 if OK, 1 if error in function                      */

int mds4_fn(TDAContext *ctx)
{
    register int i,j,k,l,il,lj;
    int err;
    double ai,bi,aj,bj,al,bl,d,t;
         
    /* get estimated distances from parameters */

    k = 0;
    for (i = 2; i <= ctx->MDSRN; ++i) {
        ai = ctx->TPar[i - 1];
        bi = ctx->TPar[i + ctx->MDSRN - 2];                       
        for (j = 1; j < i; ++j) {
            if (j == 1)  
                aj = bj = 0.0;
            else {
                aj = ctx->TPar[j - 1];
                bj = ctx->TPar[j + ctx->MDSRN - 2];                       
            }
            d = mds_edis(ctx, ai,bi,aj,bj);
            ctx->MDSDIS1[k] = d;
            d = ctx->MDSDIS[k] - d;
            if (ctx->LFunc)
                ctx->FTmp += d * d;
            k++;
        }
    }
    if (ctx->LGrad) {
        for (l = 2; l <= ctx->MDSRN; ++l) {
            al = ctx->TPar[l - 1];
            bl = ctx->TPar[l + ctx->MDSRN - 2];                       
            lj = ((l - 1) * (l - 2)) / 2;
            t = ctx->MDSDIS[lj];
            if ((d = mds_elen(ctx, al,bl)) != 0.0)
                t /= d;
            ctx->Grad[l-1] -= 2.0 * al * (t - 1.0);

            for (j = 2; j < l; ++j) {
                aj = ctx->TPar[j - 1];
                bj = ctx->TPar[j + ctx->MDSRN - 2];                       
                lj = ((l - 1) * (l - 2)) / 2 + j - 1;
                t = ctx->MDSDIS[lj];
                if ((d = mds_edis(ctx, al,bl,aj,bj)) != 0.0)
                    t /= d;
                ctx->Grad[l-1] -= 2.0 * (al - aj) * (t - 1.0);
            }
            for (i = l + 1; i <= ctx->MDSRN; ++i) {
                ai = ctx->TPar[i - 1];
                bi = ctx->TPar[i + ctx->MDSRN - 2];                       
                il = ((i - 1) * (i - 2)) / 2 + l - 1;
                t = ctx->MDSDIS[il];
                if ((d = mds_edis(ctx, ai,bi,al,bl)) != 0.0)
                    t /= d;
                ctx->Grad[l-1] += 2.0 * (ai - al) * (t - 1.0);
            }
        }
        for (l = 2; l <= ctx->MDSRN; ++l) {
            al = ctx->TPar[l - 1];
            bl = ctx->TPar[l + ctx->MDSRN - 2];                       
            lj = ((l - 1) * (l - 2)) / 2;
            t = ctx->MDSDIS[lj];
            if ((d = mds_elen(ctx, al,bl)) != 0.0)
                t /= d;
            ctx->Grad[ctx->MDSRN + l - 2] -= 2.0 * bl * (t - 1.0);

            for (j = 2; j < l; ++j) {
                aj = ctx->TPar[j - 1];
                bj = ctx->TPar[j + ctx->MDSRN - 2];                       
                lj = ((l - 1) * (l - 2)) / 2 + j - 1;
                t = ctx->MDSDIS[lj];
                if ((d = mds_edis(ctx, al,bl,aj,bj)) != 0.0)
                    t /= d;
                ctx->Grad[ctx->MDSRN + l - 2] -= 2.0 * (bl - bj) * (t - 1.0);
            }
            for (i = l + 1; i <= ctx->MDSRN; ++i) {
                ai = ctx->TPar[i - 1];
                bi = ctx->TPar[i + ctx->MDSRN - 2];                       
                il = ((i - 1) * (i - 2)) / 2 + l - 1;
                t = ctx->MDSDIS[il];
                if ((d = mds_edis(ctx, ai,bi,al,bl)) != 0.0)
                    t /= d;
                ctx->Grad[ctx->MDSRN + l - 2] += 2.0 * (bi - bl) * (t - 1.0);
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
/*  mdsn        Nonmetric MDS                                               */
/*              Requires undirected valued graph.                           */
/*                                                                          */
/*              mdsn(                                                       */
/*                  gn=...,    graph number, def. 1                         */
/*                  ndim=...   number of dimensions (1,2,3), def. 2         */
/*                  opt=...,   1 or 2 or 3, def. 1                          */
/*                             1 = ties ignored                             */
/*                             2 = Kruskal's primary approach               */
/*                             3 = Kruskal's secondary approach             */
/*                  mxit=...,  maximal number of iteration, default 100     */
/*                  tolg=...,  tolerance for gradient, def. 1.e-6           */ 
/*                  slen=...,  initial step length, def. 0.2                */
/*                  tols=...,  tolerance for stress, def. 1.e-4             */
/*                  ns=...,    number of random restarts, def. 1            */
/*                  xp=...,    starting values                              */
/*                  dvs=...,   starting values                              */
/*                  df=...,    additional output file containing the        */
/*                             estimated distance matrix                    */
/*                  pcf=...,   additional output file with a version of     */
/*                             Shepard's diagram                            */
/*                  fmt=...,   print format, def. 10.4                      */
/*              ) = fname;                                                  */
/*                                                                          */
/*  a) If the stress value of a configuration is not greater than tols,     */ 
/*     a perfect match will be assumed.                                     */
/*                                                                          */
/*  b) No repeats are performed if the user supplies starting values.       */ 
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */
       
int mdsn(TDAContext *ctx)
{
    register int i,j,k;
    int err,ii,nties,iter;
    double d,a,b,dm,st = 0.0,s,t,stmin = 0.0,ngrad;
    double mag,slen,stp = 0.0,ngradp = 0.0,af,gf,rf,cosphi,st5;
    double sthis[6];

    ctx->SLen = 0.2;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Nonmetric MDS. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto MDSNFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto MDSNFin;
    if (gdd_tcheck(ctx, 0,1,2))          /* need undirected valued graph */
        goto MDSNFin;

    if (ctx->GD_NP < 3) {
        printf1(ctx, "Error: need at least three nodes.\n");
        goto MDSNFin;
    }
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    if (ctx->MxItFlg == 0)
        ctx->MxIter = 100;

    ctx->MDSRN = ctx->GD_NP;                       /* number of points */
    ctx->MDSRNN = ctx->MDSRN * (ctx->MDSRN - 1) / 2;
    ctx->MDSND = ctx->PMNDIM;
    if (ctx->MDSND < 0)
        ctx->MDSND = 2;
    ctx->MDSNP = ctx->MDSRN * ctx->MDSND;

    printf1(ctx, "Number of points: %d\n",ctx->MDSRN);
    printf1(ctx, "Number of dimensions: %d\n",ctx->MDSND);
    if (ctx->MDSND < 1 || ctx->MDSND > 3) {
        tda_out("Error in number of dimensions.\n");
        goto MDSNFin;
    }

    /* allocate memory */

    if (alloc_acv(ctx, ctx->MDSRNN + 1))      /* used for MDSDIS (sorted) */             
        goto MDSNFin; 
    if (alloc_acw(ctx, ctx->MDSRNN + 1))      /* used for MDSDIS1 */            
        goto MDSNFin; 
    if (alloc_acu(ctx, ctx->MDSRNN + 1))      /* used for MDSDIS2 */            
        goto MDSNFin; 
    if (alloc_act(ctx, ctx->MDSRNN + 1))      /* used for MDSDIS3 */            
        goto MDSNFin; 
    if (alloc_acs(ctx, ctx->MDSRNN + 1))      /* used for MDSPTR (ptr for MDSDIS) */             
        goto MDSNFin; 
    if (alloc_acr(ctx, ctx->MDSRNN + 1))      /* inverse pointer */             
        goto MDSNFin; 
    if (alloc_acm(ctx, ctx->MDSRNN + 1))      /* used for tie blocks */             
        goto MDSNFin; 
    if (alloc_aci(ctx, ctx->MDSRNN + 1))                                            
        goto MDSNFin; 
    if (alloc_acj(ctx, ctx->MDSRNN + 1))                                            
        goto MDSNFin; 
    if (alloc_acx(ctx, ctx->MDSNP + 1))       /* used for configuration (= parameters) */
        goto MDSNFin; 
    if (alloc_acz(ctx, ctx->MDSNP + 1))       /* used for best configuration */
        goto MDSNFin; 
    if (alloc_acy(ctx, ctx->MDSNP + 1))       /* used for gradient */
        goto MDSNFin; 
    if (alloc_actmp(ctx, ctx->MDSNP + 1))     /* used for previous gradient */
        goto MDSNFin; 

    ctx->MDSDIS = ctx->AcV;
    /* private storage instead of aliasing the shared Ac* scratch:
       a later alloc_* resize frees those buffers underneath the
       aliases (a use-after-free at the projection write-back).  Freed in the
       function's Fin path. */
    ctx->MDSDIS1 = (double *)calloc((size_t)ctx->MDSRNN + 1,sizeof(double));
    ctx->MDSDIS2 = (double *)calloc((size_t)ctx->MDSRNN + 1,sizeof(double));
    ctx->MDSDIS3 = (double *)calloc((size_t)ctx->MDSRNN + 1,sizeof(double));
    if (!ctx->MDSDIS1 || !ctx->MDSDIS2 || !ctx->MDSDIS3) {
        p_err(ctx, -2,1);
        goto MDSNFin;
    }
    ctx->MDSPTR = ctx->AcS;
    ctx->MDSPTRI = ctx->AcR;
    ctx->MDSX = ctx->AcX;
    ctx->MDSG = ctx->AcY;
    ctx->MDSGP = ctx->AcTmp;
    ctx->MDSXB = ctx->AcZ;        

    /* put lower triangle of distance matrix into MDSDIS */

    dm = 0.0;
    k = 0;
    for (i = 2; i <= ctx->MDSRN; ++i) {
        for (j = 1; j < i; ++j) {
            d = gdd_adj(ctx, i-1,j-1,ctx->PMGN);
            ctx->MDSDIS[k++] = d;
            dm = dmax(ctx, dm,d);
        }
    }

    printf1(ctx, "Largest distance value: %lg\n",dm);

    /* sort MDSDIS in ascending order. pointer: MDSPTR */

    if (sortdp(ctx, ctx->MDSRNN,ctx->MDSDIS,ctx->MDSPTR))
        goto MDSNFin;

    /* make inverse pointer in MDSPTRI */

    for (i = 0; i < ctx->MDSRNN; ++i)
        ctx->MDSPTRI[ctx->MDSPTR[i]] = i;

    /* rearrange MDSDIS into sorted order */

    for (i = 0; i < ctx->MDSRNN; ++i)
        ctx->MDSDIS1[i] = ctx->MDSDIS[ctx->MDSPTR[i]];

    for (i = 0; i < ctx->MDSRNN; ++i)
        ctx->MDSDIS[i] = ctx->MDSDIS1[i];





    nties = 0;
    a = -1.0;
    for (i = 0; i < ctx->MDSRNN; ++i) {
        b = ctx->MDSDIS[i];
        if (b <= a + ctx->EPSI1)
            nties++;
        else
            a = b;
    }
    printf1(ctx, "Number of ties: %d\n",nties);
    if (nties > 0) {
        printf1(ctx, "Handling of ties: ");
        if (ctx->PMOPT == 2)
            printf1(ctx, "Kruskal's primary approach.\n");
        else if (ctx->PMOPT == 3)
            printf1(ctx, "Kruskal's secondary approach.\n");
        else {
            ctx->PMOPT = 1;
            printf1(ctx, "ignored.\n");
        }
    }
    nties = mreg_ties(ctx, ctx->MDSRNN,ctx->MDSDIS,ctx->AcM);


    /* PMNS defaults to -1 (unset) until an ns= option sets it; the
       repeat loop below is bounded by PMNS, so leaving it unclamped
       skips the loop entirely, never computing a configuration and
       printing stmin/MDSXB uninitialized. Clamped exactly as the
       sibling mds functions (mdsc, mdsm) clamp the same field. */
    if (ctx->PMNS < 1)
        ctx->PMNS = 1;
      

    printf1(ctx, "\nParameters for iterative minimization.\n");
    printf1(ctx, "Maximal number of iterations: %d\n",ctx->MxIter);
    printf1(ctx, "Initial step length: %lg\n",ctx->SLen);
    printf1(ctx, "Tolerance for norm of gradient: %lg\n",ctx->TOLG);
    printf1(ctx, "Tolerance for stress: %lg\n\n",ctx->TOLS);

    printf1(ctx, "Number of random repeats: %d\n\n",ctx->PMNS);

    for (ii = 1; ii <= ctx->PMNS; ++ii) {    /* repeat */

        printf1(ctx, "\nRepeat %d\n",ii);

        if (ctx->DSVFlg == 0) {              /* get initial configuration */
            for (i = 1; i <= ctx->MDSNP; ++i) 
                ctx->MDSX[i] = dm * random1(ctx);
        }
/**
        tda_out("CONFIG: ");
        for (i = 1; i <= MDSNP; ++i)
            tda_out("%6.3lf ",MDSX[i]);
        newline(ctx);
**/
        slen = ctx->SLen;

        for (iter = 1; iter <= ctx->MxIter; ++iter) {

            /* normalize the configuration and create distances in MDSDIS1 */

            if (mdsn_norm(ctx))
                goto MDSNFin;

/**
            tda_out("MDSDIS1 nach norm: ");
            for (i = 0; i < MDSRNN; ++i)
                tda_out("%6.3lf ",MDSDIS1[i]);
            newline(ctx);
**/

            /* find monotone projection in MDSDIS2 */

            for (i = 0; i < ctx->MDSRNN; ++i)
                ctx->MDSDIS2[i] = ctx->MDSDIS1[i];
     
            if (ctx->PMOPT == 1)
                mreg_proj(ctx, ctx->MDSRNN,ctx->MDSDIS2);
    
            else if (ctx->PMOPT == 2)
                mreg_proj1(ctx, ctx->MDSRNN,nties,ctx->AcM,ctx->AcI,ctx->AcJ,ctx->MDSDIS2,ctx->MDSDIS3);
      
            else if (ctx->PMOPT == 3)
                mreg_proj2(ctx, ctx->MDSRNN,nties,ctx->AcM,ctx->MDSDIS2,ctx->MDSDIS3);

/**
            tda_out("MDSDIS2 nach proj: ");
            for (i = 0; i < MDSRNN; ++i)
                tda_out("%6.3lf ",MDSDIS2[i]);
            newline(ctx);
**/


            /* calculate stress */ 
    
            s = t = 0.0;
            for (i = 0; i < ctx->MDSRNN; ++i) {
                a = ctx->MDSDIS1[i];
                b = ctx->MDSDIS2[i];
                s += (a - b) * (a - b);
                t += a * a;
            }
            if (fabs(t) <= ctx->EPSI1) {
                printf1(ctx, "Error in stress calculation.\n");
                goto MDSNFin;
            }
            st = sqrt(s/t);

            if (st <= ctx->TOLS) {
                printf1(ctx, "Found configuration with stress: %lg\n",st);
                break;
            }
            mdsn_grad(ctx, s,t,st);   

            ngrad = 0.0;
            for (i = 1; i <= ctx->MDSNP; ++i)  
                ngrad += ctx->MDSG[i] * ctx->MDSG[i];
            if (ngrad > 0.0)
                ngrad = sqrt(ngrad);

            printf1(ctx, "Iteration %3d  stress %lg  norm of grad: %lg\n",iter,st,ngrad);
            if (ngrad <= ctx->TOLG) {
                printf1(ctx, "Convergence reached in %d iterations.\n",iter);
                printf1(ctx, "Final norm of gradient: %lg\n",ngrad);
                break; 
            }

            /* update configuration */

            a = 0.0;
            for (i = 1; i <= ctx->MDSNP; ++i)   
                a += ctx->MDSX[i] * ctx->MDSX[i];

            if (a <= ctx->EPSI1) {
                printf1(ctx, "Error in distances. Can't continue.\n");
                goto MDSNFin;
            }
            mag = ngrad / sqrt(a);

            /* Kruskal's (1964) step length: the step is slen times the
               configuration's own scale along the unit gradient, and
               slen is adapted from the angle between successive
               gradients (angle factor), the stress ratio over the last
               five steps (relaxation factor) and over the last step
               (good-luck factor).  The original left this block
               commented out with a fixed step, which stalls far from
               the minimum; see doc/changes-from-tda.md. */
            if (iter == 1)
                af = gf = 1.0;
            else {
                gf = dmin(ctx, 1.0,st / stp);
                a = 0.0;
                for (i = 1; i <= ctx->MDSNP; ++i)
                    a += ctx->MDSG[i] * ctx->MDSGP[i];
                cosphi = (ngrad > 0.0 && ngradp > 0.0) ?
                         a / (ngrad * ngradp) : 0.0;
                af = pow(4.0,pow(cosphi,3.0));
            }
            if (iter <= 5) {
                sthis[iter] = st;
                rf = 0.65;
            }
            else {
                st5 = sthis[1];
                for (i = 1; i <= 4; ++i)
                    sthis[i] = sthis[i + 1];
                sthis[5] = st;
                rf = 1.3 / (1.0 + pow(dmin(ctx, 1.0,st / st5),5.0));
            }
            slen *= af * rf * gf;

            for (i = 1; i <= ctx->MDSNP; ++i)
                ctx->MDSX[i] -= slen * ctx->MDSG[i] / mag;

            for (i = 1; i <= ctx->MDSNP; ++i)    /* save gradient */
                ctx->MDSGP[i] = ctx->MDSG[i];
            stp = st;
            ngradp = ngrad;

        }
           
        printf1(ctx, "Best stress value (in repeat %d): %lg\n",ii,st);
                        
        /* print result into output file */
                     
        fprintf(ctx->PMFd,"%4d %12.8lf  ",ii,st);
        for (i = 1; i <= ctx->MDSNP; ++i)
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->MDSX[i]);
        fprintf(ctx->PMFd,"\n");

        /* save best configuration */

        if (ii == 1) {
            stmin = st;
            for (i = 1; i <= ctx->MDSNP; ++i)
                ctx->MDSXB[i] = ctx->MDSX[i];
        }
        else if (st < stmin) {
            stmin = st;     
            for (i = 1; i <= ctx->MDSNP; ++i)
                ctx->MDSXB[i] = ctx->MDSX[i];
        }
    }

    printf1(ctx, "\nFinal best stress value: %lg\n",stmin);
    printf1(ctx, "Coordinates\n");
    for (i = 1; i <= ctx->MDSRN; ++i) {
        printf1(ctx, "%3d ",i);
        for (j = 1; j <= ctx->MDSND; ++j)  
            rt_printf1_d(ctx, ctx->PMFmtS,ctx->MDSXB[(i - 1) * ctx->MDSND + j]);
        newline(ctx);
    }
    newline(ctx);

    /* the distances and the monotone projection of the best
       configuration, for the two optional files */
    for (i = 1; i <= ctx->MDSNP; ++i)
        ctx->MDSX[i] = ctx->MDSXB[i];
    mdsn_dis(ctx, ctx->MDSDIS1);
    for (i = 0; i < ctx->MDSRNN; ++i)
        ctx->MDSDIS2[i] = ctx->MDSDIS1[i];
    if (ctx->PMOPT == 1)
        mreg_proj(ctx, ctx->MDSRNN,ctx->MDSDIS2);
    else if (ctx->PMOPT == 2)
        mreg_proj1(ctx, ctx->MDSRNN,nties,ctx->AcM,ctx->AcI,ctx->AcJ,ctx->MDSDIS2,ctx->MDSDIS3);
    else if (ctx->PMOPT == 3)
        mreg_proj2(ctx, ctx->MDSRNN,nties,ctx->AcM,ctx->MDSDIS2,ctx->MDSDIS3);

    if (ctx->PMF1Def) {      /* write estimated distance matrix */
        for (i = 1; i <= ctx->MDSRN; ++i) {
            for (j = 1; j <= ctx->MDSRN; ++j) {
                d = 0.0;
                for (k = 1; k <= ctx->MDSND; ++k) {
                    a = ctx->MDSXB[(i - 1) * ctx->MDSND + k] -
                        ctx->MDSXB[(j - 1) * ctx->MDSND + k];
                    d += a * a;
                }
                if (d > 0.0)
                    d = sqrt(d);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,d);
            }
            fprintf(ctx->PMF1d,"\n");
        }
        printf1(ctx, "%d records written to %s\n",ctx->MDSRN,ctx->PMF1dName);
    }
    if (ctx->PMPCFDef) {     /* write Shepard's diagram: dissimilarity,
                                fitted distance, monotone projection, in
                                ascending order of the dissimilarity */
        for (i = 0; i < ctx->MDSRNN; ++i) {
            rt_fprintf_d(ctx, ctx->PMPCFd,ctx->PMFmtS,ctx->MDSDIS[i]);
            rt_fprintf_d(ctx, ctx->PMPCFd,ctx->PMFmtS,ctx->MDSDIS1[i]);
            rt_fprintf_d(ctx, ctx->PMPCFd,ctx->PMFmtS,ctx->MDSDIS2[i]);
            fprintf(ctx->PMPCFd,"\n");
        }
        printf1(ctx, "%d records written to %s\n",ctx->MDSRNN,ctx->PMPCFName);
    }
    printf1(ctx, "%d records written to %s\n",ctx->PMNS,ctx->PMFdName);

    err = 0;

MDSNFin:       
    free(ctx->MDSDIS1); ctx->MDSDIS1 = NULL;
    free(ctx->MDSDIS2); ctx->MDSDIS2 = NULL;
    free(ctx->MDSDIS3); ctx->MDSDIS3 = NULL;
    ml_init(ctx, 0,0,0);      
    p_clean(ctx);
    return(err);
}



/* ------------------------------------------------------------------------ */
/*  mdsn_dis(dis)       Create distances from configuration MDSX in the     */
/*                      array dis[] in sorted order.                        */

void mdsn_dis(TDAContext *ctx, double *dis)
{
    register int i,j,k,l;
    double d,t;

    l = 0;
    for (i = 2; i <= ctx->MDSRN; ++i) {
        for (j = 1; j < i; ++j) {
            d = 0.0;
            for (k = 1; k <= ctx->MDSND; ++k) {
                t = ctx->MDSX[(i - 1) * ctx->MDSND + k] - ctx->MDSX[(j - 1) * ctx->MDSND + k];  
                d += t * t;
            }
            if (d > 0.0)
                d = sqrt(d);

            dis[ctx->MDSPTRI[l++]] = d;
        }
    }
}


/* ------------------------------------------------------------------------ */
/*  mdsn_norm()     Normalize the configuration in MDSX and update the      */
/*                  distances in MDSDIS1.                                   */
/*                                                                          */
/*  Return 0 if OK, 1 if error.                                             */

int mdsn_norm(TDAContext *ctx)
{
    register int i,j;
    double s;             
/**
    tda_out("X vorher: ");
    for (i = 1; i <= MDSNP; ++i)
        tda_out("%6.3lf ",MDSX[i]);
    newline(ctx);
**/

    for (j = 1; j <= ctx->MDSND; ++j) {
        s = 0.0;
        for (i = 1; i <= ctx->MDSRN; ++i)   
            s += ctx->MDSX[(i - 1) * ctx->MDSND + j];
        s /= (double)ctx->MDSRN;
        for (i = 1; i <= ctx->MDSRN; ++i)   
            ctx->MDSX[(i - 1) * ctx->MDSND + j] -= s;
    }
/********
    tda_out("X nach 1: ");
    for (i = 1; i <= MDSNP; ++i)
        tda_out("%6.3lf ",MDSX[i]);
    newline(ctx);
**/
    mdsn_dis(ctx, ctx->MDSDIS1);    /* calculate distances */


/**
    tda_out("MDSDIS1 nach 1 sorted: ");
    for (i = 0; i < MDSRNN; ++i)
        tda_out("%6.3lf ",MDSDIS1[i]);
    newline(ctx);

    tda_out("MDSDIS1 nach 1 unsorted: ");
    for (i = 0; i < MDSRNN; ++i)
        tda_out("%6.3lf ",MDSDIS1[MDSPTRI[i]]);
    newline(ctx);
**/


    s = 0.0;
    for (i = 0; i < ctx->MDSRNN; ++i)  
        s += ctx->MDSDIS1[i] * ctx->MDSDIS1[i];

    if (s <= ctx->EPSI1) {
        printf1(ctx, "Error in distances. Can't continue.\n");
        return(1);          
    }
    s = sqrt((double)ctx->MDSRNN / s);

/**
    printf1(ctx, "s=======%lg\n",s);
**/
    for (i = 1; i <= ctx->MDSNP; ++i)
        ctx->MDSX[i] *= s;

/**
    tda_out("X nach 2: ");
    for (i = 1; i <= MDSNP; ++i)
        tda_out("%6.3lf ",MDSX[i]);
    newline(ctx);
**/


    mdsn_dis(ctx, ctx->MDSDIS1);    /* calculate new distances */

/**

    s = 0.0;
    for (i = 0; i < MDSRNN; ++i)  
        s += MDSDIS1[i] * MDSDIS1[i];

    s = sqrt((double)MDSRNN / s);
    printf1(ctx, "new s=======%lg\n",s);

    tda_out("MDSDIS1 nach 2 sorted: ");
    for (i = 0; i < MDSRNN; ++i)
        tda_out("%6.3lf ",MDSDIS1[i]);
    newline(ctx);

    tda_out("MDSDIS1 nach 2 unsorted: ");
    for (i = 0; i < MDSRNN; ++i)
        tda_out("%6.3lf ",MDSDIS1[MDSPTRI[i]]);
    newline(ctx);


**/


    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mdsn_grad() Calculate gradient of stress function. Assume distances     */
/*              in MDSDIS1, projection values in MDSDIS2. Return gradient   */
/*              in MDSG. Assume stress = sqrt(s/t).                         */
/*              Return 0 if OK, 1 if error.                                 */

int mdsn_grad(TDAContext *ctx, double s,double t,double stress)    
{
    register int i,j,k,l,ij,kl;
    double dij,dijp,xil,xjl,g,dki,dkj;

    kl = 0;
    for (k = 1; k <= ctx->MDSRN; ++k) {
        for (l = 1; l <= ctx->MDSND; ++l) {

            g = 0.0;
            ij = 0;
            for (i = 2; i <= ctx->MDSRN; ++i) {
                if (i == k)
                    dki = 1.0;
                else
                    dki = 0.0;
                for (j = 1; j < i; ++j) {
                    if (j == k)
                        dkj = 1.0;
                    else
                        dkj = 0.0;

                    dij = ctx->MDSDIS1[ctx->MDSPTRI[ij]];
                    dijp = ctx->MDSDIS2[ctx->MDSPTRI[ij]];

/**
tda_out("k=%2d l=%2d ij=%3d i=%2d j=%2d dkj=%lg dki=%lg dij=%6.3lf dijp=%6.3lf\n", k,l,ij,i,j,dkj,dki,dij,dijp);
**/
                    ij++;

                    if (dij <= 0.0 || dijp <= 0.0)  {
                        continue;
                    }
                    xil = ctx->MDSX[(i - 1) * ctx->MDSND + l];
                    xjl = ctx->MDSX[(j - 1) * ctx->MDSND + l];
                                
                    g += (dki - dkj) * ((dij - dijp) / s - dij / t) *
                            (xil - xjl) / dij;
                }
            }
            ctx->MDSG[++kl] = g * stress;
        }
    }
    return(0);
}


/* ------------------------------------------------------------------------ */
/*  mdsn1       Nonmetric MDS                                               */
/*              Requires undirected valued graph.                           */
/*                                                                          */
/*              mdsn(                                                       */
/*                  gn=...,    graph number, def. 1                         */
/*                  opt=...,   1 or 2 or 3, def. 1                          */
/*                             1 = ties ignored                             */
/*                             2 = Kruskal's primary approach               */
/*                             3 = Kruskal's secondary approach             */
/*                  mxit=...,  maximal number of iteration, default 100     */
/*                  tolg=...,  tolerance for gradient, def. 1.e-6           */ 
/*                  slen=...,  initial step length, def. 0.2                */
/*                  tols=...,  tolerance for stress, def. 1.e-4             */
/*                  ns=...,    number of random restarts, def. 1            */
/*                  xp=...,    starting values                              */
/*                  dvs=...,   starting values                              */
/*                  df=...,    additional output file containing the        */
/*                             estimated distance matrix                    */
/*                  pcf=...,   additional output file with a version of     */
/*                             Shepard's diagram                            */
/*                  fmt=...,   print format, def. 10.4                      */
/*              ) = fname;                                                  */
/*                                                                          */
/*  a) If the stress value of a configuration is not greater than tols,     */ 
/*     a perfect match will be assumed.                                     */
/*                                                                          */
/*  b) No repeats are performed if the user supplies starting values.       */ 
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */
       
int mdsn1(TDAContext *ctx)
{
    register int i,j,k;
    int err,ii,nties,np,iter;
    double d,a,b,aa,bb,dm,st = 0.0,s,t,stmin = 0.0,fval,ngrad;
    double mag,slen,stp = 0.0,ngradp = 0.0,af,gf,rf,cosphi,st5;
    double sthis[6];

    ctx->SLen = 0.2;

/*  TOLG = 1.e-4;   */

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Nonmetric MDS. Current memory: %d bytes.\n",ctx->MemReq);

    /*  The command name is "mdsn1", five characters, so its option
        block starts at +5.  This read from +4, leaving the trailing
        "1" in front of the options, so every call was a syntax
        error -- written when the function was expected to be
        reached as "mdsn".  */
    if (parm(ctx, ctx->CmdBuf + 5,1,1))       /* get parameters */
        goto MDSNFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto MDSNFin;
    if (gdd_tcheck(ctx, 0,1,2))          /* need undirected valued graph */
        goto MDSNFin;

    if (ctx->GD_NP < 3) {
        printf1(ctx, "Error: need at least three nodes.\n");
        goto MDSNFin;
    }
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    if (ctx->MxItFlg == 0)
        ctx->MxIter = 100;

    ctx->MDSRN = ctx->GD_NP;                       /* number of points */
    ctx->MDSRNN = ctx->MDSRN * (ctx->MDSRN - 1) / 2;

    /* setup MDSDIS vector with distances */

    if (alloc_acv(ctx, ctx->MDSRNN + 1))      /* used for MDSDIS */             
        goto MDSNFin; 
    if (alloc_acw(ctx, ctx->MDSRNN + 1))      /* used for MDSDIS1 */            
        goto MDSNFin; 
    if (alloc_acu(ctx, ctx->MDSRNN + 1))      /* used for MDSDIS2 */            
        goto MDSNFin; 
    if (alloc_act(ctx, ctx->MDSRNN + 1))      /* used for MDSDIS3 */            
        goto MDSNFin; 
    if (alloc_acz(ctx, ctx->MDSRNN + 1))      /* used to save best configuration */
        goto MDSNFin; 
    if (alloc_acs(ctx, ctx->MDSRNN + 1))      /* used for MDSPTR */             
        goto MDSNFin; 

    ctx->MDSDIS = ctx->AcV;
    /* private storage instead of aliasing the shared Ac* scratch:
       a later alloc_* resize frees those buffers underneath the
       aliases (a use-after-free at the projection write-back).  Freed in the
       function's Fin path. */
    ctx->MDSDIS1 = (double *)calloc((size_t)ctx->MDSRNN + 1,sizeof(double));
    ctx->MDSDIS2 = (double *)calloc((size_t)ctx->MDSRNN + 1,sizeof(double));
    ctx->MDSDIS3 = (double *)calloc((size_t)ctx->MDSRNN + 1,sizeof(double));
    if (!ctx->MDSDIS1 || !ctx->MDSDIS2 || !ctx->MDSDIS3) {
        p_err(ctx, -2,1);
        goto MDSNFin;
    }
    ctx->MDSPTR = ctx->AcS;

    np = 2 * (ctx->MDSRN - 1);       /* number of parameters */

    if (alloc_acx(ctx, np + 1))      /* used for configuration (= parameters) */
        goto MDSNFin; 
    if (alloc_acy(ctx, np + 1))      /* used for gradient */           
        goto MDSNFin; 
    if (alloc_actmp(ctx, np + 1))    /* used for previous gradient */    
        goto MDSNFin; 

    /* put lower triangle of distance matrix into MDSDIS */

    dm = 0.0;
    k = 0;
    for (i = 2; i <= ctx->MDSRN; ++i) {
        for (j = 1; j < i; ++j) {
            d = gdd_adj(ctx, i-1,j-1,ctx->PMGN);
            ctx->MDSDIS[k++] = d;
            dm = dmax(ctx, dm,d);
        }
    }

    printf1(ctx, "Largest distance value: %lg\n",dm);

    /* sort MDSDIS in ascending order. pointer: MDSPTR */

    if (sortdp(ctx, ctx->MDSRNN,ctx->MDSDIS,ctx->MDSPTR))
        goto MDSNFin;


    printf1(ctx, "MDSPTR: ");
    for (i = 0; i < ctx->MDSRNN; ++i)
        tda_out("%d ",ctx->MDSPTR[i]);
    newline(ctx);


    nties = 0;
    a = -1.0;
    for (i = 0; i < ctx->MDSRNN; ++i) {
        b = ctx->MDSDIS[ctx->MDSPTR[i]];
        if (b <= a + ctx->EPSI1)
            nties++;
        else
            a = b;
    }
    printf1(ctx, "Number of ties: %d\n",nties);
    if (nties > 0) {
        printf1(ctx, "Handling of ties: ");
        if (ctx->PMOPT == 2)
            printf1(ctx, "Kruskal's primary approach.\n");
        else if (ctx->PMOPT == 3)
            printf1(ctx, "Kruskal's secondary approach.\n");
        else {
            ctx->PMOPT = 1;
            printf1(ctx, "ignored.\n");
        }
    }

/**
    nties = mreg_ties(ctx, MDSRNN,double *y,int *t);
**/


    /***
    tda_out("\nNumber of tie blocks m=%d\n",m);
    for (i = 0; i < m; ++i)
        tda_out("i=%3d t[i]=%d\n",i,t[i]);
    **/

    /* PMNS defaults to -1 (unset) until an ns= option sets it; the
       repeat loop below is bounded by PMNS, so leaving it unclamped
       skips the loop entirely, never computing a configuration and
       printing stmin/MDSXB uninitialized. Clamped exactly as the
       sibling mds functions (mdsc, mdsm) clamp the same field. */
    if (ctx->PMNS < 1)
        ctx->PMNS = 1;

    if (alloc_act(ctx, ctx->MDSRNN + 1))      /* used for best coordinates */   
        goto MDSNFin; 

    printf1(ctx, "\nParameters for iterative minimization.\n");
    printf1(ctx, "Maximal number of iterations: %d\n",ctx->MxIter);
    printf1(ctx, "Initial step length: %lg\n",ctx->SLen);
    printf1(ctx, "Tolerance for norm of gradient: %lg\n",ctx->TOLG);
    printf1(ctx, "Tolerance for stress: %lg\n\n",ctx->TOLS);

    printf1(ctx, "Number of random repeats: %d\n\n",ctx->PMNS);

    for (ii = 1; ii <= ctx->PMNS; ++ii) {    /*  repeat */

        printf1(ctx, "\nRepeat %d\n",ii);

        if (ctx->DSVFlg == 0) {              /* get initial configuration */
            for (i = 1; i <= np; ++i) 
                ctx->AcX[i] = dm * random1(ctx);
        }

        slen = ctx->SLen;

        for (iter = 1; iter <= ctx->MxIter; ++iter) {

            /* create distances of configuration in MDSDIS1 */

            mdsn1_dis(ctx, ctx->MDSRN,ctx->AcX,ctx->MDSDIS1);

            /* normalize the configuration */

            a = 0.0;
            for (i = 0; i < ctx->MDSRNN; ++i)  
                a += ctx->MDSDIS1[i] * ctx->MDSDIS1[i];

            if (a <= ctx->EPSI1) {
                printf1(ctx, "Error in distances. Can't continue.\n");
                goto MDSNFin;
            }
            b = sqrt((double)ctx->MDSRNN / a);

            for (i = 1; i <= np; ++i)
                ctx->AcX[i] *= b;

            mdsn1_dis(ctx, ctx->MDSRN,ctx->AcX,ctx->MDSDIS1);


            a = 0.0;
            for (i = 0; i < ctx->MDSRNN; ++i)  
                a += ctx->MDSDIS1[i] * ctx->MDSDIS1[i];

            if (a <= ctx->EPSI1) {
                printf1(ctx, "Error in distances. Can't continue.\n");
                goto MDSNFin;
            }
            b = sqrt((double)ctx->MDSRNN / a);




            /* find monotone projection in MDSDIS3 */

            mdsn1_proj(ctx, ctx->PMOPT,ctx->MDSRNN,ctx->MDSDIS1,ctx->MDSDIS2,ctx->MDSDIS3);

            /* calculate stress */ 
    
            s = t = 0.0;
            for (i = 0; i < ctx->MDSRNN; ++i) {
                a = ctx->MDSDIS1[i];
                b = ctx->MDSDIS3[i];
                s += (a - b) * (a - b);
                t += a * a;
            }
            if (fabs(t) <= ctx->EPSI1) {
                printf1(ctx, "Error in stress calculation.\n");
                goto MDSNFin;
            }
            st = sqrt(s/t);
            printf1(ctx, "Iteration %3d  stress %lg\n",iter,st);
      
            if (st <= ctx->TOLS) {
                printf1(ctx, "Found configuration with stress: %lg\n",st);
                break;
            }
    
            if (mdsn1_fn(ctx, ctx->MDSRN,ctx->MDSRNN,ctx->AcX,&fval,ctx->AcY,ctx->MDSDIS1,ctx->MDSDIS3)) {
                printf1(ctx, "Error in function evaluation. Can't continue.\n");
                goto MDSNFin;
            }
            ngrad = 0.0;
            for (i = 1; i <= np; ++i)
                ngrad += ctx->AcY[i] * ctx->AcY[i];

            if (ngrad > 0.0)
                ngrad = sqrt(ngrad);

            if (ngrad <= ctx->TOLG) {
                printf1(ctx, "Convergence reached with %d iterations.\n",iter);
                printf1(ctx, "Final norm of gradient: %lg\n",ngrad);
                break; 
            }

            /* update configuration */

            a = 0.0;
            for (i = 1; i <= np; ++i)   
                a += ctx->AcX[i] * ctx->AcX[i];

            if (a <= ctx->EPSI1) {
                printf1(ctx, "Error in distances. Can't continue.\n");
                goto MDSNFin;
            }
            mag = ngrad / sqrt(a);

            /* Kruskal's step length adaptation, as in mdsn() above */
            if (iter == 1)
                af = gf = 1.0;
            else {
                gf = dmin(ctx, 1.0,st / stp);
                a = 0.0;
                for (i = 1; i <= np; ++i)
                    a += ctx->AcY[i] * ctx->AcTmp[i];
                cosphi = (ngrad > 0.0 && ngradp > 0.0) ?
                         a / (ngrad * ngradp) : 0.0;
                af = pow(4.0,pow(cosphi,3.0));
            }
            if (iter <= 5) {
                sthis[iter] = st;
                rf = 0.65;
            }
            else {
                st5 = sthis[1];
                for (i = 1; i <= 4; ++i)
                    sthis[i] = sthis[i + 1];
                sthis[5] = st;
                rf = 1.3 / (1.0 + pow(dmin(ctx, 1.0,st / st5),5.0));
            }
            slen *= af * rf * gf;

            for (i = 1; i <= np; ++i)
                ctx->AcX[i] -= slen * ctx->AcY[i] / mag;

            for (i = 1; i <= np; ++i)       /* save gradient */
                ctx->AcTmp[i] = ctx->AcY[i];
            stp = st;
            ngradp = ngrad;

        }

        printf1(ctx, "Best stress value (in repeat %d): %lg\n",ii,st);
                        
        /* print result into output file */
                     
        fprintf(ctx->PMFd,"%4d %12.8lf  ",ii,st);
        for (i = 1; i <= np; ++i)
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcX[i]);
        fprintf(ctx->PMFd,"\n");

        /* save best configuration */

        if (ii == 1) {
            stmin = st;
            for (i = 1; i <= np; ++i)
                ctx->AcZ[i] = ctx->AcX[i];
        }
        else if (st < stmin) {
            stmin = st;     
            for (i = 1; i <= ctx->NParm; ++i)
                ctx->AcZ[i] = ctx->AcX[i];
        }
        /*******************
        if (DSVFlg)        no repeats if starting values 
            break;
        *******************/
    }
    printf1(ctx, "\nFinal best stress value: %lg\n",stmin);
    printf1(ctx, "Coordinates\n");
    for (i = 1; i <= ctx->MDSRN; ++i) {
        printf1(ctx, "%3d ",i);
        if (i == 1)
            a = b = 0.0;
        else {
            a = ctx->AcZ[i - 1];
            b = ctx->AcZ[ctx->MDSRN + i - 2];
        }
        rt_printf1_d(ctx, ctx->PMFmtS,a);
        rt_printf1_d(ctx, ctx->PMFmtS,b);
        newline(ctx);
    }
    newline(ctx);

    if (ctx->PMF1Def) {      /* write estimated distance matrix */
        for (i = 1; i <= ctx->MDSRN; ++i) {
            if (i == 1)
                a = b = 0.0;
            else {
                a = ctx->AcZ[i - 1];
                b = ctx->AcZ[ctx->MDSRN + i - 2];
            }
            for (j = 1; j <= ctx->MDSRN; ++j) {
                if (j == 1)
                    aa = bb = 0.0;
                else {
                    aa = ctx->AcZ[j - 1];
                    bb = ctx->AcZ[ctx->MDSRN + j - 2];
                }
                d = mds_edis(ctx, a,b,aa,bb);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,d);
            }
            fprintf(ctx->PMF1d,"\n");
        }
        printf1(ctx, "%d records written to %s\n",ctx->MDSRN,ctx->PMF1dName);
    }
    if (ctx->PMPCFDef) {     /* write Shepard's diagram */
        k = 0;
        for (i = 2; i <= ctx->MDSRN; ++i) {
            a = ctx->AcX[i - 1];
            b = ctx->AcX[ctx->MDSRN + i - 2];
            for (j = 1; j < i; ++j) {
                if (j == 1)
                    aa = bb = 0.0;
                else {
                    aa = ctx->AcX[j - 1];
                    bb = ctx->AcX[ctx->MDSRN + j - 2];
                }
                ctx->MDSDIS1[k++] = mds_edis(ctx, a,b,aa,bb);
            }
        }
        for (i = 0; i < ctx->MDSRNN; ++i) {
            rt_fprintf_d(ctx, ctx->PMPCFd,ctx->PMFmtS,ctx->MDSDIS[ctx->MDSPTR[i]]);
            rt_fprintf_d(ctx, ctx->PMPCFd,ctx->PMFmtS,ctx->MDSDIS1[ctx->MDSPTR[i]]);
            fprintf(ctx->PMPCFd,"\n");
        }
        printf1(ctx, "%d records written to %s\n",ctx->MDSRNN,ctx->PMPCFName);
    }
    printf1(ctx, "%d records written to %s\n",ctx->PMNS,ctx->PMFdName);

    err = 0;

MDSNFin:       
    free(ctx->MDSDIS1); ctx->MDSDIS1 = NULL;
    free(ctx->MDSDIS2); ctx->MDSDIS2 = NULL;
    free(ctx->MDSDIS3); ctx->MDSDIS3 = NULL;
    ml_init(ctx, 0,0,0);      
    p_clean(ctx);
    return(err);
}



/* ------------------------------------------------------------------------ */
/*  mdsn1_dis(n,par,dis) Create distances from configuration in par[]        */
/*                      in the array dis[].                                 */

void mdsn1_dis(TDAContext *ctx, int n,double *par,double *dis)
{
    register int i,j,k;
    double a,b,aa,bb;

    k = 0;
    for (i = 2; i <= n; ++i) {
        a = par[i - 1];
        b = par[n + i - 2];

        for (j = 1; j < i; ++j) {
            if (j == 1)
                aa = bb = 0.0;
            else {
                aa = par[j - 1];
                bb = par[n + j - 2];
            }
            dis[k++] = mds_edis(ctx, a,b,aa,bb);
        }
    }


}

/* ------------------------------------------------------------------------ */
/*  mdsn1_proj()     Find projection, depending on opt.                      */
/*                  nn = n (n-1) / 2 (n = number of points)                 */

void mdsn1_proj(TDAContext *ctx, int opt,int nn,double *dis,double *tmp,double *proj)
{
    register int i;


    /* find monotone projection */
    
    for (i = 0; i < nn; ++i)
        tmp[i] = dis[ctx->MDSPTR[i]];
    
    if (opt == 1)
        mreg_proj(ctx, nn,tmp);
    
    /**********
    else if (PMOPT == 2)
        mreg_proj1(ctx, NOC,NTies,AcM,AcI,AcJ,CJRY,AcTmp);
    
    else if (PMOPT == 3)
        mreg_proj2(ctx, NOC,NTies,AcM,CJRY,AcTmp);
    **********/
            
    /* replace according to original order */
    
    for (i = 0; i < nn; ++i)
        proj[ctx->MDSPTR[i]] = tmp[i];
    
}   

/* ------------------------------------------------------------------------ */
/*  mdsn1_fn()  Calculate function value and gradient. n is number of       */
/*              points, nn = n(n-1)/2. par[] contains 2(n-1) parameters.    */  
/*              Return function value in fval, gradient in grad[].          */
/*                                                                          */
/*              Return 0 if OK, 1 if error in function                      */

int mdsn1_fn(TDAContext *ctx, int n,int nn,double *par,double *fval,double *grad,double *dis,double *disp)
{
    (void)dis; (void)disp;        /* unused: the signature is shared */
    register int i,j,k,l,il,lj;
    double ai,bi,aj,bj,al,bl,d,s,t,st,sal,sbl,tal,tbl;

    /* get estimated distances from parameters, use MDSDIS1 */

    k = 0;
    for (i = 2; i <= n; ++i) {
        ai = par[i - 1];
        bi = par[i + n - 2];                       
        for (j = 1; j < i; ++j) {
            if (j == 1)  
                aj = bj = 0.0;
            else {
                aj = par[j - 1];
                bj = par[j + n - 2];                       
            }
            d = mds_edis(ctx, ai,bi,aj,bj);
            ctx->MDSDIS1[k] = d;
            d = ctx->MDSDIS1[k] - d;
            k++;
        }
    }

    /* calculate stress terms */

    s = t = 0.0;
    for (i = 0; i < nn; ++i) {
        ai = ctx->MDSDIS1[i];
        bi = ctx->MDSDIS3[i];
        s += (ai - bi) * (ai - bi);
        t += ai * ai;
    }
    if (fabs(s) <= ctx->EPSI1 || fabs(t) <= ctx->EPSI1) {
        printf1(ctx, "\nError in function mdsn1_fn().\n");
        return(1);
    }
    st = sqrt(s / t);

    *fval = st;

    for (l = 2; l <= n; ++l) {
        al = par[l - 1];
        bl = par[l + n - 2];                       
        tal = tbl = sal = sbl = 0.0;
        for (j = 1; j < l; ++j) {
            lj = ((l - 1) * (l - 2)) / 2 + j - 1;
            if (j == 1)
                aj = bj = 0.0;
            else {
                aj = par[j - 1];
                bj = par[j + n - 2];                       
            }
            tal += 2.0 * (al - aj);
            tbl += 2.0 * (bl - bj);
            sal += 2.0 * (al - aj) * (1.0 - ctx->MDSDIS3[lj] / ctx->MDSDIS1[lj]);
            sbl += 2.0 * (bl - bj) * (1.0 - ctx->MDSDIS3[lj] / ctx->MDSDIS1[lj]);
        }
        for (i = l + 1; i <= n; ++i) {
            il = ((i - 1) * (i - 2)) / 2 + l - 1;
            ai = par[i - 1];
            bi = par[i + n - 2];                       
            tal -= 2.0 * (ai - al);
            tbl -= 2.0 * (bi - bl);
            sal += 2.0 * (ai - al) * (ctx->MDSDIS3[il] / ctx->MDSDIS1[il] - 1.0);
            sbl += 2.0 * (bi - bl) * (ctx->MDSDIS3[il] / ctx->MDSDIS1[il] - 1.0);
        }
        /* assigned, not accumulated: the caller never clears grad[],
           and with += the array carried the sum of every gradient
           since the command started (doc/changes-from-tda.md) */
        grad[l - 1] = 0.5 * st * (sal / s - tal / t);
        grad[l + n - 2] = 0.5 * st * (sbl / s - tbl / t);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mds_edis(x1,y1,x2,y2) Return euclidean distance of (x1,y1) and (x2,y2). */

double mds_edis(TDAContext *ctx, double x1,double y1,double x2,double y2)
{
    (void)ctx;        /* unused: the signature is shared */
    return(sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)));
}

/* ------------------------------------------------------------------------ */
/*  mds_elen(x,y)    Return euclidean length of (x,y).                      */

double mds_elen(TDAContext *ctx, double x,double y)
{
    (void)ctx;        /* unused: the signature is shared */
    return(sqrt(x * x + y * y));
}

/* ------------------------------------------------------------------------ */
/*  mdsx        Multidimensional scaling.                                   */
/*              Requires undirected valued graph.                           */
/*                                                                          */
/*              mdsx(                                                       */
/*                  opt=...,    option, def. 1                              */
/*                              1 = minmax approach                         */
/*                              2 = least squares                           */
/*                  gn = ...,   graph number, def. 1                        */
/*                  mxit=...,   max number of iterations, def. 1000         */
/*                  tolv=...,   tolerance for variance, def. 1.e-10         */
/*                  mxitl=..,   max iterations for step 0, def. 10000       */
/*                  nbox=...,   max number of boxes, def. 10000             */
/*                  tolbw=...,  tolerance box width, def. 0.01              */
/*                                                                          */
/*                  sc=...,     box width in verify step, def. 1.0          */
/*                  tolfd=...,  final min box width, def. 1e-3              */
/*                  fmt=...,    print format, def. 10.4                     */
/*                  pcf=...,    create plot command file                    */
/*              ) = fname;                                                  */
/*                                                                          */
/*  If successful, fname will contain the final coordinates and radii.      */
/*  1. index                                                                */
/*  2. first coordinate                                                     */
/*  3. second coordinate                                                    */
/*  4. radius                                                               */
/*                                                                          */
/*  If sc > 0 the command executes an additional algorithm ...              */


/*  Return: 0 if OK, -1 if error.                                           */
       
int mdsx(TDAContext *ctx)
{
    register int i,j;
    int err,rn,ip,ip0,ip1,nn,sflag;
    double r,dm,dmk,fmin0 = 0.0;

double a,b,d;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Multidimensional scaling. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLBW = 0.01;                   /* default box width */
    ctx->TOLFD = 0.1;       

    /* "mdsx" is four characters, so the parser has to start at +4.  At +3 it
       saw the trailing 'x' and rejected every invocation that carried options
       with "Error: x(...)".  mdsm, mdsn, mdsc and mdsr all use +4; this one
       did not, in the 6.4 release. */
    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto MDSXFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto MDSXFin;
    if (gdd_tcheck(ctx, 0,1,2))          /* need undirected valued graph */
        goto MDSXFin;

    if (ctx->GD_NP < 3) {
        printf1(ctx, "Error: need at least three nodes.\n");
        goto MDSXFin;
    }
    if (ctx->PMFmtF == 0)                /* default print format */
        pmfmt(ctx, 10,4);

    if (ctx->PMOPT == 1) {
        printf1(ctx, "Minmax approach.\n\n");
        ctx->MDSTYP = 0;
    }
    else {
        ctx->MDSTYP = 1;
        printf1(ctx, "Least squares approach.\n\n");
    }

    ctx->MDSRN = ctx->GD_NP;                  /* number of points */

    if (ctx->MxItFlg == 0)               /* max number of iterations */
        ctx->MxIter = 1000;

    /* setup parameters for simplex algorithm */

    ctx->NParm = 2 * ctx->MDSRN - 3;          /* number of parameters */
    ctx->MINA = 2;                       /* always this algorithm */
    ctx->PMDOPT = 0;                     /* without derivatives */
    ctx->CCTyp = 0;                      /* no covariance matrix */
    set_mlopt(ctx);                    /* adjust options */         

    if (ml_init(ctx, 1,1,0))             /* init function minimization */
        goto MDSXFin;


    /* setup parameters for intervall algorithm */

    if (ctx->MxIt1 <= 50)        /* def. max number of iterations */
        ctx->MxIt1 = 10000;
    if (ctx->PMNBOX < 1)         /* def. max number of boxes */
        ctx->PMNBOX = 10000;

    printf1(ctx, "Parameters for step 0 algorithm.\n");
    printf1(ctx, "Max number of iterations: %d\n",ctx->MxIt1);
    printf1(ctx, "Maximum number of boxes: %d\n",ctx->PMNBOX);
    printf1(ctx, "Tolerance for box width: %g\n\n",ctx->TOLBW);
            
    if (alloc_par(ctx, 2 * ctx->MDSRN,1))         /* allocate arrays for parameters */
        goto MDSXFin;

    if (alloc_list(ctx, ctx->PMNBOX,2 * ctx->MDSRN))   /* allocate boxes */
        goto MDSXFin;

    nn = ctx->MDSRN * (ctx->MDSRN - 1) / 2;

    if (alloc_acx(ctx, ctx->MDSRN + 1))       /* first coordinate (a) */
        goto MDSXFin; 
    if (alloc_acy(ctx, ctx->MDSRN + 1))       /* second coordinate (b) */
        goto MDSXFin; 
    if (alloc_acu(ctx, ctx->MDSRN + 1))     
        goto MDSXFin; 
    if (alloc_actmp(ctx, ctx->MDSRN + 1))     
        goto MDSXFin; 
    if (alloc_acn(ctx, ctx->MDSRN + 1))     
        goto MDSXFin; 
    if (alloc_acm(ctx, ctx->MDSRN + 1))     
        goto MDSXFin; 
    if (alloc_acxf(ctx, ctx->MDSRN * ctx->MDSRN + 1))     
        goto MDSXFin; 
    if (alloc_aci(ctx, nn + 1))     
        goto MDSXFin; 
    if (alloc_acj(ctx, nn + 1))     
        goto MDSXFin; 
    if (alloc_acw(ctx, nn + 1))     
        goto MDSXFin; 
    if (alloc_ack(ctx, nn + 1))     
        goto MDSXFin; 

    ctx->MDSRI = ctx->AcN;
    ctx->MDSRF = ctx->AcM;
    ctx->MDSRD = ctx->AcXF;
                     
    /* find points ip0 and ip1 that have maximal distance */

    dm = 0.0;
    ip0 = ip1 = -1;

    for (i = 0; i < ctx->MDSRN; ++i) {
        for (j = 0; j < ctx->MDSRN; ++j) {
            if (i != j) {
                r = gdd_adj(ctx, i,j,ctx->PMGN);
                if (r >= 0.0) {
                    ctx->MDSRD[i * ctx->MDSRN + j] = (float)(r);                  
                    if (r > dm) {
                        dm = r;
                        ip0 = i;
                        ip1 = j;
                    }
                }
            }
        }
    }
    if (ip0 > ip1) {
        i = ip0;
        ip0 = ip1;
        ip1 = i;
    }
    ctx->MDSRIP0 = ip0;
    ctx->MDSRIP1 = ip1;

    tda_out("dm=%g ip0=%d ip1=%d  \n",dm,ip0,ip1);
    for (i = 0; i < ctx->MDSRN; ++i) {
        for (j = 0; j < ctx->MDSRN; ++j)
            tda_out("%f ",(double)(ctx->MDSRD[i * ctx->MDSRN + j]));
        newline(ctx);
    }
    newline(ctx);

    if (ip0 < 0) {
        printf1(ctx, "Error: no positive distances.\n");
        goto MDSXFin;
    }
    ctx->MDSRNI = 2;
    ctx->MDSRI[0] = ip0;
    ctx->MDSRI[1] = ip1;
    ctx->MDSRF[ip0] = 1;
    ctx->MDSRF[ip1] = 1;

    ctx->AcX[ip0] = 1.0;             /* first point */
    ctx->AcY[ip0] = 1.0;

    ctx->AcX[ip1] = ctx->AcX[ip0] + dm;   /* second point */
    ctx->AcY[ip1] = 1.0;
        
    sflag = ctx->SILENTFlg;
    ctx->SILENTFlg = 2;

    printf1(ctx, "Step  NPnt  Iterations  NBoxes  Accepted  Variance   Min Function\n");
    if (sflag < 2)
        printfe(ctx, "\nStep  NPnt  Iterations  NBoxes  Accepted  Variance   Min Function\n");

    /*  do for all further points. fmin0 is the current best minimum */
    /*  found in the first step. */

    while (ctx->MDSRNI < ctx->MDSRN) {

        ip = 0;              /* find next point */
        dmk = -1.0;
        for (i = 0; i < ctx->MDSRN; ++i) {                       
            if (ctx->MDSRF[i])
                continue;

            r = 0.0;
            for (j = 0; j < ctx->MDSRNI; ++j)  
                r = (double)(dmax(ctx, r,(double)ctx->MDSRD[i * ctx->MDSRN + j]));                 
            if (r > dmk) {
                ip = i;
                dmk = r;
            }
        }
        ctx->MDSRI[ctx->MDSRNI] = ctx->MDSRIP = ip;  
        ctx->MDSRF[ip] = 1;
        ctx->MDSRNI++;

        /* find starting values with interval algorithm */

        ctx->RParL[0] = ctx->AcX[ip0];
        ctx->RParU[0] = ctx->AcX[ip1];
        ctx->RParL[1] = ctx->AcY[ip0];                 
        ctx->RParU[1] = ctx->AcY[ip0] + dm;                  
        ctx->RPar[0] = (ctx->RParL[0] + ctx->RParU[0]) / 2.0;
        ctx->RPar[1] = (ctx->RParL[1] + ctx->RParU[1]) / 2.0;

        rn = mdsxr_m(ctx, 2,ctx->PMNBOX,ctx->RPar,ctx->RParL,ctx->RParU,ctx->AcTmp);             
       
        if (rn) {
            printf1(ctx, "Error %d in mdsxr_m.\n",rn);
            goto MDSXFin;
        }
       



        ctx->AcX[ip] = ctx->RPar[0];
        ctx->AcY[ip] = ctx->RPar[1];
        fmin0 = mdsxr_fm(ctx);

        printf1(ctx, "%4d  %4d  %10d %7d %9d           %14.7e\n",0,ctx->MDSRNI,ctx->GO_IT,ctx->GO_NB,ctx->GO_NA,fmin0);
        if (sflag < 2)
            printfe(ctx, "%4d  %4d  %10d %7d %9d           %14.7e\n",0,ctx->MDSRNI,ctx->GO_IT,ctx->GO_NB,ctx->GO_NA,fmin0);


tda_out("\nNEU nach step 0 a,b,r  \n");
        for (i = 0; i < ctx->MDSRN; ++i)
            tda_out("%f %f %f\n",ctx->AcX[i],ctx->AcY[i],ctx->AcU[i]);


        rn = mds_m(ctx, ip,ip0,ip1);         /* call minimization */
        if (rn) {
            printf1(ctx, "Error %d in mdsxr_m.\n",rn);
            goto MDSXFin;
        }
        fmin0 = ctx->FMin;

        printf1(ctx, "%4d  %4d  %10d                   %9.2e %14.7e\n",1,ctx->MDSRNI,ctx->Iter,ctx->CValV,fmin0);
        if (sflag < 2)
            printfe(ctx, "%4d  %4d  %10d                   %9.2e %14.7e\n",1,ctx->MDSRNI,ctx->Iter,ctx->CValV,fmin0);
    }
    ctx->SILENTFlg = sflag;
    printf1(ctx, "\nOptimal value of criterion: %g\n",fmin0);
    if (ctx->LConv)  
        printf1(ctx, "Warning: final step didn't converge.\n");
    newline(ctx);

    /*  find optimal radii */

    mds_r(ctx, ctx->MDSRN,ctx->MDSRD,ctx->AcX,ctx->AcY,ctx->AcU,ctx->MDSRI,ctx->AcW,ctx->AcI,ctx->AcJ,ctx->AcK);

tda_out("\nNEU Final values a,b,r  \n");
    for (i = 0; i < ctx->MDSRN; ++i)
        tda_out("%f %f %f\n",ctx->AcX[i],ctx->AcY[i],ctx->AcU[i]);

    newline(ctx);

    printf1(ctx, "Neue Abstandsmatrix\n");
    for (i = 0; i < ctx->MDSRN; ++i) {
        for (j = 0; j < ctx->MDSRN; ++j) {
            a = ctx->AcX[i] - ctx->AcX[j];
            b = ctx->AcY[i] - ctx->AcY[j];
            d = sqrt(a * a + b * b);
            printf1(ctx, "%8.4lf ",d);
        }
        newline(ctx);
    }


    for (i = 0; i < ctx->MDSRN; ++i) {       /* write to output file */
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcX[i]);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcY[i]);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcU[i]);
        fprintf(ctx->PMFd,"\n");
    }
    printf1(ctx, "%d records written to: %s\n",ctx->MDSRN,ctx->PMFdName);

    /* create plot command file */

    if (ctx->PMPCFDef)       
        mds_p(ctx, ctx->MDSRN,ctx->AcX,ctx->AcY,ctx->AcU);


    /**********  *********/

    if (ctx->PMSC > 0.0) {
        printf1(ctx, "\nStarting additional algorithm.\n");
        printf1(ctx, "Initial box width (sc): %g\n",ctx->PMSC);
        printf1(ctx, "Minimal box width (tolfd): %g\n",ctx->TOLFD);

        if (alloc_acu(ctx, 2 * ctx->MDSRN + 1))      
            goto MDSXFin; 
        if (alloc_acv(ctx, 2 * ctx->MDSRN + 1))   
            goto MDSXFin; 

        for (i = 0; i < ctx->MDSRN; ++i) {
            ctx->RPar[i] = ctx->AcX[i];
            if (i == ip0)  
                ctx->RParL[i] = ctx->RParU[i] = ctx->RPar[i];
            else {
                ctx->RParL[i] = ctx->RPar[i] - ctx->PMSC / 2.0;
                ctx->RParU[i] = ctx->RPar[i] + ctx->PMSC / 2.0 + ctx->EPSI1;
            }     
            j = i + ctx->MDSRN;
            ctx->RPar[j] = ctx->AcY[i];
            if (i == ip0 || i == ip1) 
                ctx->RParL[j] = ctx->RParU[j] = ctx->RPar[j];
            else {
                ctx->RParL[j] = ctx->RPar[j] - ctx->PMSC / 2.0;
                ctx->RParU[j] = ctx->RPar[j] + ctx->PMSC / 2.0 + ctx->EPSI1;
            }       
        }
        for (i = 0; i < 2 * ctx->MDSRN; ++i)
            tda_out("i=%3d %20.16f %20.16f %20.16f\n",i,ctx->RPar[i],ctx->RParL[i],ctx->RParU[i]);

/**
        MDSFMIN = fmin0;
        rn = mdsxr_c2(ctx, 2 * MDSRN,RParL,RParU,AcU,AcV);
**/


        ctx->MDSFMIN = ctx->PMSC; 
        rn = mdsxr_c3(ctx, 2 * ctx->MDSRN,ctx->RParL,ctx->RParU,ctx->AcU,ctx->AcV);
        tda_out("rn=%d\n",rn);

        for (i = 0; i < 2 * ctx->MDSRN; ++i)
            tda_out("i=%3d %20.16f %20.16f\n",i,ctx->AcU[i],ctx->AcV[i]);


        tda_out("MDSNCALL: %d\n",ctx->MDSNCALL);
        tda_out("MDSNSKIP: %d\n",ctx->MDSNSKIP);

        if (ctx->PMPCFDef)       
            mds_p1(ctx, 2 * ctx->MDSRN,ctx->AcU,ctx->AcV);

    }
    err = 0;

MDSXFin:       
    alloc_list(ctx, 0,0);
    alloc_par(ctx, 0,0);
    ml_init(ctx, 0,0,0);      
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  mds_m       Function minimization with simplex method, used by mds().   */
/*              Two steps:                                                  */
/*              1) find optimal position for new point given the position   */
/*                 of previous points.                                      */
/*              2) optimize positions of all points.                        */

int mds_m(TDAContext *ctx, int ip,int ip0,int ip1)
{
    (void)ip0; (void)ip1;        /* unused: the signature is shared */
    register int i,ii,k;
    int fflag;
    double a = 0.0,b = 0.0;


/**     
    tda_out("new point is %d\n",ip);
    tda_out("distance ip0 to new: %f\n",d0n);
    tda_out("distance ip1 to new: %f\n",d1n);
    tda_out("distance ip0 to ip1: %f\n",d01);

tda_out("current XY vor startwert \n");
    for (i = 0; i < MDSRN; ++i)
        tda_out("%d %f %f\n",i,AcX[i],AcY[i]);
newline(ctx);

**/     
         
fflag = 0;
    if (fflag) {
        ctx->NParm1 = ctx->NParm = 2;
        ctx->Par[1] = a;
        ctx->Par[2] = b;

        ctx->MDSMSTEP = 0;

        if (ffmin(ctx, MDSMOD3,1,0,1))
            return(-1);
        ctx->AcX[ip] = ctx->TPar[1];
        ctx->AcY[ip] = ctx->TPar[2];

        tda_out("NEWxPar: %g %g FMIN=%g  \n",ctx->Par[1],ctx->Par[2],ctx->FMin );
    }


    ctx->MDSMSTEP = 1;
    ctx->NParm1 = ctx->NParm = ctx->MDSRNI * 2 - 3;

    k = 0;
    for (i = 0; i < ctx->MDSRNI; ++i) {
        ii = ctx->MDSRI[i];
        if (ii == ctx->MDSRIP0)
            continue;
        else if (ii == ctx->MDSRIP1)
            ctx->Par[++k] = ctx->AcX[ii];
        else {
            ctx->Par[++k] = ctx->AcX[ii];
            ctx->Par[++k] = ctx->AcY[ii];
        }
    }

/**       
tda_out("k=%d NParm=%d\n",k,NParm);

    tda_out("PAR\n");
    for (i = 1; i <= NParm; ++i)
        tda_out("%d par=%f\n",i,Par[i]);

tda_out("current XY nach min  \n");
    for (i = 0; i < MDSRN; ++i)
        tda_out("%d %f %f\n",i,AcX[i],AcY[i]);
newline(ctx);

***/     
           
    if (ffmin(ctx, MDSMOD3,1,0,1))
        return(-1);

    /* update values in AcX and AcY with new parameters */

    k = 0;
    for (i = 0; i < ctx->MDSRNI; ++i) {
        ii = ctx->MDSRI[i];
        if (ii == ctx->MDSRIP0)
            continue;
        else if (ii == ctx->MDSRIP1)
            ctx->AcX[ii] = ctx->Par[++k];                 
        else {
            ctx->AcX[ii] = ctx->Par[++k];
            ctx->AcY[ii] = ctx->Par[++k];
        }
    }
       
tda_out("NEW XY nach step 1 fmin=%g  \n",ctx->FMin );
    for (i = 0; i < ctx->MDSRN; ++i)
        tda_out("%d %f %f\n",i,ctx->AcX[i],ctx->AcY[i]);
newline(ctx);

       

    return(0);
                
}

/* ------------------------------------------------------------------------ */
/*  mds3_fn()   Calculate function value for MDSMOD2, depends on MDSMSTEP.  */
/*              Use parameter values in TPar[] (i=1,NParm).                 */
/*              Return function value in FTmp.                              */
/*                                                                          */
/*              MDSTYP = 0 : minmax                                         */
/*                       1 : least squares                                  */
/*                                                                          */
/*              Return 0 if OK, 1 if error in function                      */

int mds3_fn(TDAContext *ctx)
{
    register int i,j,ii,jj,ip,jp;
    int err;
    double tmp,dij,ai,aj,bi,bj,aij,bij;
      
/**      
    tda_out("mds3_fn NParm=%d %d \n",NParm,NParm1   );
    tda_out("PAR: ");
    for (i = 1; i <= NParm; ++i)
        tda_out("%g ",TPar[i]);
    newline(ctx);
**/      

    ctx->FTmp = 0.0;

/**
    for (i = 1; i <= NParm; ++i) {
        if (TPar[i] < 10.0 * EPSI1) {
            FTmp = sqrt(DBLMAX);                
tda_out("HIER\n");
           
            return(0);
        }
    }
**/

    if (ctx->MDSMSTEP == 0) {
                      
        for (i = 0; i < ctx->MDSRNI; ++i) {
            ip = ctx->MDSRI[i];
            if (ip == ctx->MDSRIP)
                continue;
            aij = ctx->AcX[ip] - ctx->TPar[1];
            bij = ctx->AcY[ip] - ctx->TPar[2];
            dij = (double)(ctx->MDSRD[ip * ctx->MDSRN + ctx->MDSRIP]);
            tmp = fabs(dij - sqrt(aij * aij + bij * bij));
/**
tda_out("i=%d ip=%d aij=%g bij=%g ai=%g bi=%g aip=%g bip=%g dij=%g tmp=%g\n", i,ip,aij,bij,AcX[ip],AcY[ip],TPar[1],TPar[2],dij ,tmp);
**/
            if (ctx->MDSTYP == 0)
                ctx->FTmp = dmax(ctx, ctx->FTmp,tmp);
            else
                ctx->FTmp += tmp * tmp;            
        }
    }
    else {
        ii = 1;
        for (i = 1; i < ctx->MDSRNI; ++i) {
            ip = ctx->MDSRI[i];
            if (ip == ctx->MDSRIP0) {
                ai = 1.0;
                bi = 1.0;
            }
            else if (ip == ctx->MDSRIP1) {
                ai = ctx->TPar[ii++];
                bi = 1.0;
            }
            else {
                ai = ctx->TPar[ii++];
                bi = ctx->TPar[ii++];
            }
            jj = 1;
            for (j = 0; j < i; ++j) {
                jp = ctx->MDSRI[j];
                if (jp == ctx->MDSRIP0) {
                    aj = 1.0;
                    bj = 1.0;
                }
                else if (jp == ctx->MDSRIP1) {
                    aj = ctx->TPar[jj++];
                    bj = 1.0;
                }
                else {
                    aj = ctx->TPar[jj++];
                    bj = ctx->TPar[jj++];
                }
                dij = (double)(ctx->MDSRD[ip * ctx->MDSRN + jp]);
                aij = ai - aj;          
                bij = bi - bj;                   
                tmp = fabs(dij - sqrt(aij * aij + bij * bij));
                if (ctx->MDSTYP == 0)
                    ctx->FTmp = dmax(ctx, ctx->FTmp,tmp);
                else
                    ctx->FTmp += tmp * tmp;            
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
/*  mds_r(n,d,x,y,r,p,ds)                                                   */
/*                                                                          */
/*  Calculate radii.                                                        */
/*                                                                          */
/*  Return 0 if OK, -1 if error                                             */

int mds_r(TDAContext *ctx, int n,float *d,double *x,double *y,double *r,int *p,double *ds, int *is,int *js,int *ptr)
{
    (void)p;        /* unused: the signature is shared */
    register int i,j,k,  ii,jj,i0;
    int nn;
    double a,b,tmp,ri,rj,rd;

    nn = 0;
    for (i = 1; i < n; ++i) {
        for (j = 0; j < i; ++j) {

            a = x[i] - x[j];
            b = y[i] - y[j];
            tmp = a * a + b * b;
            if (tmp > 0.0)
                tmp = sqrt(tmp);
            ds[nn] = fabs((double)d[i * n + j] - tmp);
            is[nn] = i;
            js[nn] = j;
            nn++;
        }
    }

    tda_out("vor sort\n");
    for (i = 0; i < nn; ++i)
        tda_out("i=%3d ds=%f is=%d js=%d\n",i,ds[i],is[i],js[i]);

newline(ctx);

    if (sortdp(ctx, nn,ds,ptr))
        return(-1);

    tda_out("nach sort\n");
    for (i = 0; i < nn; ++i) {
        k = ptr[i];
        tda_out("i=%3d ds=%f is=%d js=%d\n",i,ds[k],is[k],js[k]);
    }
newline(ctx);

              

    for (i = 0; i < n; ++i)
        r[i] = 0.0;

    for (i = nn - 1; i >= 0; --i) {
        k = ptr[i];
        ii = is[k];
        jj = js[k];
        tmp = ds[k];

        ri = r[ii];
        rj = r[jj];

        if (tmp > ri + rj) {
            tmp -= (ri + rj);

            if (ri != rj) {
                if (ri < rj) {
                    rd = rj - ri;
                    i0 = ii;
                }
                else {
                    rd = ri - rj;
                    i0 = jj;
                }
                rd = dmin(ctx, rd,tmp);
                r[i0] += rd;
                tmp -= rd;
            }
            if (tmp > 0.0) {
                r[ii] += tmp / 2.0;
                r[jj] += tmp / 2.0;
            }
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mds_p   Create a command file for plotting the configuration.           */

void mds_p(TDAContext *ctx, int n,double *x,double *y,double *r)
{
    register int i,m;
    int nx,ny,icx,icy;
    double xmin,xmax,ymin,ymax,xa,xb,ya,yb,rmin,rmax,tx,ty,scx,scy;  

    rmax = r[0];
    for (i = 1; i < n; ++i)
        rmax = dmax(ctx, rmax,r[i]);

    xmin = xmax = x[0];
    ymin = ymax = y[0];
                   
    for (i = 1; i < n; ++i) {
        xmin = dmin(ctx, xmin,x[i]);
        xmax = dmax(ctx, xmax,x[i]);
        ymin = dmin(ctx, ymin,y[i]);
        ymax = dmax(ctx, ymax,y[i]);
    }
    rmin = 0.01 * xmax;

    tx = 0.1 * (xmax + xmin + 2.0 * rmax) / 2.0;
    ty = 0.1 * (ymax + ymin + 2.0 * rmax) / 2.0;

    nx = x_getaxval(ctx, xmin-tx,xmax+tx,&xa,&xb,0);     /* get values for axes */
    ny = x_getaxval(ctx, ymin-ty,ymax+ty,&ya,&yb,0);

    icx = 0;
    scx = pow(10.0,(double)nx);
    if ((xb - xa) / scx > 20.0) {    
        scx = pow(10.0,(double)(nx+1));
        icx = 10.0;
    }
    icy = 0;
    scy = pow(10.0,(double)ny);
    if ((yb - ya) / scy > 20.0) {    
        scy = pow(10.0,(double)(ny+1));
        icy = 10.0;
    }

    fprintf(ctx->PMPCFd,"psfile = %s.ps;\n",ctx->PMPCFName);
    fprintf(ctx->PMPCFd,"psetup(pxlen=100,pylen=100,pxa=%g,%g,pya=%g,%g);\n",xa,xb,ya,yb);
    fprintf(ctx->PMPCFd,"plxa(sc=%g,ic=%d);\n",scx,icx);      
    fprintf(ctx->PMPCFd,"plya(sc=%g,ic=%d);\n",scy,icy);

    m = 0;
    for (i = 0; i < n; ++i) {
        if (r[i] >= rmin) {
            fprintf(ctx->PMPCFd,"ploto(xy=%g,%g,lt=1,nc=1)=%g;\n",x[i],y[i],r[i]);
            m++;
        }
    }
    if (m < n) {
        nx = 0;
        for (i = 0; i < n; ++i) {
            if (r[i] < rmin) {
                if (nx == 0) {
                    fprintf(ctx->PMPCFd,"plotp(s=5,fs=1,nc=1,lw=0)=%g,%g",x[i],y[i]);
                    nx = 1;
                }
                else
                    fprintf(ctx->PMPCFd,",%g,%g",x[i],y[i]);
            }
        }
        fprintf(ctx->PMPCFd,";\n");
    }
    printf1(ctx, "Output file for plot: %s\n",ctx->PMPCFName);
}

/* ------------------------------------------------------------------------ */
/*  mds_p1  Create a command file for plotting the configuration.           */

void mds_p1(TDAContext *ctx, int n,double *a,double *b)
{
    register int i,m;
    double x,y,xw,yw;

    m = n / 2;
    for (i = 0; i < m; ++i) {
        x = a[i];
        y = a[i + m];
        xw = b[i] - a[i];
        yw = b[i + m] - a[i + m];
        if (xw <= ctx->EPSI1)
            xw = 1000.0 * ctx->EPSI1;
        if (yw <= ctx->EPSI1)
            yw = 1000.0 * ctx->EPSI1;

        fprintf(ctx->PMPCFd,"plrec()=%g,%g,%g,%g;\n",x,y,xw,yw);
    }
}

/* ------------------------------------------------------------------------ */
/*  mdsxr_m(narg,nbmax,par,lb,ub,par0)                                       */
/*                                                                          */
/*  The parameters for the best function value are returned in par[].       */
/*  Bounds are return in lb[] and ub[].                                     */
/*                                                                          */
/*  Return  0 if successful                                                 */
/*         -1 error in function evaluation                                  */ 
/*          1 if exceeded max number of boxes                               */
/*          2 if number of accepted boxes is empty                          */
/*          3 if exceeded max number of iterations                          */

int mdsxr_m(TDAContext *ctx, int narg,int nbmax,double *par,double *lb,double *ub,double *par0)
{
    (void)par0;        /* unused: the signature is shared */
    register int i,j;
    int err,r,nb,is,iter,ja,jb,jj,first,na;
    double fl,tmp,w,f,fminp;
    double *lptra,*lptrb,*uptra,*uptrb;

    ctx->GO_IT = 0;              /* number of iterations performed */
    ctx->GO_FN = 0;              /* number of function evaluations */
    ctx->GO_IFN = 0;             /* number of inclusion function evaluations */
    ctx->GO_NBU = 0;             /* number of boxes used */
    ctx->GO_NB = 0;              /* number of boxes in final list */

    err = 0;
    for (i = 0; i < nbmax; ++i)
        ctx->GO_FLG[i] = 0;

    lptra = ctx->GO_LB;                      /* put initial box on list */
    uptra = ctx->GO_UB;

    for (i = 0; i < narg; ++i) {
        lptra[i] = lb[i];
        uptra[i] = ub[i];
    }
    ctx->GO_FLG[0] = -1;     /* not yet processed */

    /* calculate inclusion function with initial boxes */
   
    ctx->GO_IFN++;

    r = mdsxr_f(ctx, 1,narg,lptra,uptra,&fl);
    if (r) {
        printf1(ctx, "Error in evaluating inclusion function.\n");
        err = -1;
        goto MDSRMFin;
    }
    ctx->GO_LBF[0] = fl;
    
    /* calculate function at starting values */

   
    ctx->GO_FN++;
    r = mdsxr_f(ctx, 0,narg,par,par,&fminp);
    if (r) {
        printf1(ctx, "Error in evaluating function.\n");
        err = -1;
        goto MDSRMFin;
    }
    nb = 1;     /* number of boxes */
    na = 0;     /* number of accepted boxes */
    iter = 0;                                     

    /***
    if (SILENTFlg < 2)
        printfe(ctx, "\n  Iter    Function Value       NBox  FCall  IFCall\n");
    ***/

    
    while (++iter <= ctx->MxIt1) {       

        j = iter / 100;

        /**
        if (SILENTFlg < 2 && 100 * j == iter)
            printfe(ctx, "%5d  %20.13e  %6d %6d %7d\n",iter,fminp,nb,GO_FN,GO_IFN);
        **/

        /* find box with lowest lower function bound,
           and check for boxes that can be dropped */

        ja = -1;
        first = 1;
        for (i = 0; i < nb; ++i) {

            if (ctx->GO_FLG[i]) {    /* check all boxes */
  
                if (ctx->GO_LBF[i] > fminp) {
                    ctx->GO_FLG[i] = 0;
                }
                else if (ctx->GO_FLG[i] < 0) {
                    if (first) {      
                        f = ctx->GO_LBF[i];
                        ja = i;
                        first = 0;    
                    }
                    else if (f > ctx->GO_LBF[i]) {
                        f = ctx->GO_LBF[i];
                        ja = i;
                    }
                }   
            }
        }
        if (ja < 0)      /* here we have no more unprocessed boxes */
            break;
                  
        /* bisect box ja and enter subboxes into list */

        jb = -1;        /* find free box */

        for (i = 0; i < nbmax; ++i) {
            if (ctx->GO_FLG[i] == 0) {
                jb = i;
                break;
            }
        }
        if (jb < 0) {   
            printf1(ctx, "Exceeded maximum number of boxes.\n");
            err = 1;
            goto MDSRM1Fin; 
        }
        if (jb >= nb) {
            nb = jb + 1;
            if (ctx->GO_NBU < nb)
                ctx->GO_NBU = nb;
        }
        lptra = ctx->GO_LB + ja * narg;
        uptra = ctx->GO_UB + ja * narg;

        lptrb = ctx->GO_LB + jb * narg;
        uptrb = ctx->GO_UB + jb * narg;
   
        is = 0;                         /* coordinate with largest width */
        w = uptra[0] - lptra[0];

        for (i = 1; i < narg; ++i) {
            tmp = uptra[i] - lptra[i];
            if (w < tmp) {
                w = tmp;
                is = i;
            }
        }
        if (w <= ctx->TOLBW) {                   /* accept this box */
            ctx->GO_FLG[ja] = 1;
            na++;
            continue;
        }

        for (i = 0; i < narg; ++i) {      /* copy box ja to jb */
            lptrb[i] = lptra[i];
            uptrb[i] = uptra[i];
        }

        tmp = lptra[is] + w / 2.0;
        uptra[is] = lptrb[is] = tmp;

        ctx->GO_FLG[ja] = -1;                   
        ctx->GO_FLG[jb] = -1;

        for (jj = 0; jj < 2; ++jj) {           /* for both boxes */

            if (jj)   
                ja = jb;

            lptra = ctx->GO_LB + ja * narg;
            uptra = ctx->GO_UB + ja * narg;

            /* evaluate inclusion function */

            ctx->GO_IFN++;
            r = mdsxr_f(ctx, 1,narg,lptra,uptra,&fl);
            if (r) {
                printf1(ctx, "Error in evaluating inclusion function.\n");
                err = -1;
                goto MDSRMFin;
            }
            if (fl > fminp) {
                ctx->GO_FLG[ja] = 0;             /* free box */
                if (ja == nb - 1)
                    nb--;    
                continue;
            }
            ctx->GO_LBF[ja] = fl;

            /* update fminp with function value for midpoint */             
                   
            for (j = 0; j < narg; ++j)
                ctx->AcTmp[j] = (lptra[j] + uptra[j]) / 2.0;

            ctx->GO_FN++;
            r = mdsxr_f(ctx, 0,narg,ctx->AcTmp,ctx->AcTmp,&f);
            if (r) {
                printf1(ctx, "Error in evaluating function.\n");
                err = -1;
                goto MDSRMFin;
            }
            if (fminp > f) {
                fminp = f;
                igmin_cpar(ctx, narg,ctx->AcTmp,par);
            }   
        }
    }
MDSRM1Fin:
    ctx->GO_NB = nb;
    ctx->GO_NA = na;
    ctx->GO_IT = iter;
    ctx->GO_FMIN = fminp;

    if (err == 0) {
        if (ctx->Iter >= ctx->MxIt1)
            err = 3;
        else if (na == 0)  
            err = 2;
    }


MDSRMFin:
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  mdsxr_f(opt,n,xl,xh,f)                                                   */
/*                                                                          */
/*  MDSTYP = 0 : minmax                                                     */
/*           1 : least squares                                              */
/*                                                                          */
/*  Calculate function depending on opt.                                    */
/*  If opt = 0 : standard function with arguments xl[i], i=0,...,n-1        */
/*               return function value in f.                                */
/*  If opt = 1 : inclusion function at [xl[i],xh[i]], i = 0,...,n-1         */
/*               return lower bound  in f.                                  */
/*                                                                          */
/*  Return 0 if OK, -1 if error                                             */

int mdsxr_f(TDAContext *ctx, int opt,int n,double *xl,double *xh,double *f)
{
    (void)n;        /* unused: the signature is shared */
    register int j,jj;
    double dij,ai,bi,aj,bj,aij,bij,sl,sh,tl,th,ul,uh,tmp;
    double aih,bih,ajh,bjh;

/***************** 
    tda_out("MDSR_F n=%d MDSRIP=%d    \n",n,MDSRIP       );    
    for (i = 0; i < n; ++i)
        printf1(ctx, "i=%d xl=%f xh=%f\n",i,xl[i],xh[i]);
**************/

    *f = 0.0;
    if (opt == 0) {
        ai = xl[0];
        bi = xl[1];
        for (j = 0; j < ctx->MDSRNI; ++j) {
            jj = ctx->MDSRI[j];
            if (jj == ctx->MDSRIP)  
                continue;

            dij = (double)(ctx->MDSRD[ctx->MDSRIP * ctx->MDSRN + jj]);
            if (dij < 0.0)
                continue;
            aj = ctx->AcX[jj];
            bj = ctx->AcY[jj];
            aij = ai - aj;
            bij = bi - bj;
            tmp = fabs(dij - sqrt(aij * aij + bij * bij));
            if (ctx->MDSTYP == 0)
                *f = dmax(ctx, *f,tmp);
            else
                *f += tmp * tmp;
        }
        return(0);
    }

    /* calculate inclusion function */

    ai = xl[0];
    bi = xl[1];
    aih = xh[0];
    bih = xh[1];

    for (j = 0; j < ctx->MDSRNI; ++j) {
        jj = ctx->MDSRI[j];
        if (jj == ctx->MDSRIP)
            continue;

        dij = (double)(ctx->MDSRD[ctx->MDSRIP * ctx->MDSRN + jj]);
        if (dij < 0.0)
            continue;

        aj = ajh = ctx->AcX[jj];
        bj = bjh = ctx->AcY[jj];

        i_sub(ctx, ai,aih,aj,ajh,&sl,&sh);       
        i_square(ctx, sl,sh,&tl,&th);

        i_sub(ctx, bi,bih,bj,bjh,&sl,&sh);       
        i_square(ctx, sl,sh,&ul,&uh);

        i_add(ctx, tl,th,ul,uh,&sl,&sh);       
        i_sqrt(ctx, sl,sh,&tl,&th);
        i_sub(ctx, dij,dij,tl,th,&sl,&sh);
        i_abs(ctx, sl,sh,&tl,&th);

        if (ctx->MDSTYP == 0)
            *f = dmax(ctx, *f,tl);
        else {
            i_square(ctx, tl,th,&ul,&uh);
            *f += ul;
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mdsxr_fm()       Return current value of objective function.             */
/*                  MDSTYP = 0 : minmax                                     */
/*                           1 : least squares                              */
/*                                                                          */

double mdsxr_fm(TDAContext *ctx)
{
    register int i,j,ii,jj;
    double f,dij,ai,bi,aj,bj,aij,bij,tmp;

    f = 0.0;
    for (i = 1; i < ctx->MDSRNI; ++i) {
        ii = ctx->MDSRI[i];
        ai = ctx->AcX[ii];
        bi = ctx->AcY[ii];
        for (j = 0; j < i; ++j) {
            jj = ctx->MDSRI[j];
            dij = (double)(ctx->MDSRD[ii * ctx->MDSRN + jj]);
            if (dij < 0.0)
                continue;

            aj = ctx->AcX[jj];
            bj = ctx->AcY[jj];
            aij = ai - aj;
            bij = bi - bj;
            tmp = fabs(dij - sqrt(aij * aij + bij * bij));
            if (ctx->MDSTYP == 0)
                f = dmax(ctx, f,tmp);
            else
                f += tmp * tmp;
        }
    }
    return(f);
}

/* ------------------------------------------------------------------------ */
/*  mdsxr_f1(opt,n,xl,xh,f,fminp)     always minmax.                         */
/*                                                                          */
/*  Calculate function depending on opt.                                    */
/*  If opt = 0 : standard function with arguments xl[i], i=0,...,n-1        */
/*               return function value in f.                                */
/*  If opt = 1 : inclusion function at [xl[i],xh[i]], i = 0,...,n-1         */
/*               return lower bound  in f.                                  */
/*                                                                          */
/*  Return 0 if OK, -1 if error                                             */

int mdsxr_f1(TDAContext *ctx, int opt,int n,double *xl,double *xh,double *f,double fminp)
{
    (void)n;        /* unused: the signature is shared */
    register int i,j;
    double aij,bij,dij,tmp,ail,aih,bil,bih,ajl,ajh,bjl,bjh,d,d1,d2,d3,d4;

/** 
        tda_out("MDSR_F n=%d MDSRIP=%d    \n",n,MDSRIP       );    
        for (i = 0; i < n; ++i)
            printf1(ctx, "i=%d xl=%f xh=%f\n",i,xl[i],xh[i]);
**/ 

    *f = 0.0;
    if (opt == 0) {
        for (i = 1; i < ctx->MDSRN; ++i) {
            for (j = 0; j < i; ++j) {
                dij = (double)(ctx->MDSRD[i * ctx->MDSRN + j]);
                if (dij < 0.0)
                    continue;
                aij = xl[i] - xl[j];
                bij = xl[i + ctx->MDSRN] - xl[j + ctx->MDSRN];
                tmp = fabs(dij - sqrt(aij * aij + bij * bij));
                *f = dmax(ctx, *f,tmp);
            }
        }
        return(0);
    }

    for (i = 1; i < ctx->MDSRN; ++i) {
        ail = xl[i];
        aih = xh[i];
        bil = xl[i + ctx->MDSRN];
        bih = xh[i + ctx->MDSRN];

        for (j = 0; j < i; ++j) {
            ajl = xl[j];
            ajh = xh[j];
            bjl = xl[j + ctx->MDSRN];
            bjh = xh[j + ctx->MDSRN];

            dij = (double)(ctx->MDSRD[i * ctx->MDSRN + j]);
            if (dij < 0.0)
                continue;

            tmp = 0.0;
            if (aih <= ajl) {
                if (fabs(bih - bjl) >= fabs(bil - bjh))  
                    d = mdsxr_dm(ctx, ail,bih,ajh,bjl);
                else
                    d = mdsxr_dm(ctx, ail,bil,ajh,bjh);

                if (dij > d) {
                    tmp = dij - d;
                }
                else {
                    if (bil > bjh || bih < bjl) {
                        if (fabs(bil - bjh) <= fabs(bih - bjl))  
                            d = mdsxr_dm(ctx, aih,bil,ajl,bjh);
                        else
                            d = mdsxr_dm(ctx, aih,bih,ajl,bjl);
                    }   
                    else
                        d = mdsxr_dm(ctx, aih,0.0,ajl,0.0);

                    if (dij < d)
                        tmp = d - dij;
                }
            }
            else if (ajh <= ail) {
                if (fabs(bih - bjl) >= fabs(bil - bjh))  
                    d = mdsxr_dm(ctx, aih,bih,ajl,bjl);
                else
                    d = mdsxr_dm(ctx, aih,bil,ajl,bjh);

                if (dij > d) {
                    tmp = dij - d;
                }
                else {
                    if (bil > bjh || bih < bjl) {
                        if (fabs(bil - bjh) <= fabs(bih - bjl))  
                            d = mdsxr_dm(ctx, ail,bil,ajh,bjh);
                        else
                            d = mdsxr_dm(ctx, ail,bih,ajh,bjl);
                    }   
                    else
                        d = mdsxr_dm(ctx, ail,0.0,ajh,0.0);

                    if (dij < d)
                        tmp = d - dij;
                }
            }
            else if (bih <= bjl) {
                if (fabs(aih - ajl) > fabs(ail - ajh))  
                    d = mdsxr_dm(ctx, aih,bil,ajl,bjh);
                else
                    d = mdsxr_dm(ctx, ail,bil,ajh,bjh);

                if (dij > d) {
                    tmp = dij - d;
                }
                else {
                    d = mdsxr_dm(ctx, 0.0,bih,0.0,bjl);
                    if (dij < d)
                        tmp = d - dij;
                }
            }
            else if (bil >= bjh) {
                if (fabs(aih - ajl) > fabs(ail - ajh))  
                    d = mdsxr_dm(ctx, aih,bih,ajl,bjl);
                else
                    d = mdsxr_dm(ctx, ail,bih,ajh,bjl);

                if (dij > d) {
                    tmp = dij - d;
                }
                else {
                    d = mdsxr_dm(ctx, 0.0,bil,0.0,bjh);
                    if (dij < d)
                        tmp = d - dij;
                }
            }
            else {
                d1 = mdsxr_dm(ctx, ail,bil,ajh,bjh);
                d2 = mdsxr_dm(ctx, aih,bil,ajl,bjh);
                d3 = mdsxr_dm(ctx, aih,bih,ajl,bjl);
                d4 = mdsxr_dm(ctx, ail,bih,ajh,bjl);
                d = dmax(ctx, dmax(ctx, d1,d2),dmax(ctx, d3,d4));
                if (dij > d) {
                    tmp = dij - d;
                }
                else {
                    ;
                }
            }
            *f = dmax(ctx, *f,tmp);
            if (*f > fminp)
                return(0);
        }
    }
    return(0);
}

double mdsxr_dm(TDAContext *ctx, double ai,double bi,double aj,double bj)
{
    (void)ctx;        /* unused: the signature is shared */
    double aij,bij;

    aij = ai - aj;
    bij = bi - bj;
    return(sqrt(aij * aij + bij * bij));
}

/* ------------------------------------------------------------------------ */
/*  mdsxr_c2()                                                               */
/*                                                                          */

int mdsxr_c2(TDAContext *ctx, int n,double *lb,double *ub,double *x,double *y)
{
    register int i;
    int is;
    double a,b,c,f,w,tmp;

    ctx->MDSNCALL++;

    mdsxr_f1(ctx, 1,n,lb,ub,&f,ctx->MDSFMIN);
    if (f > ctx->MDSFMIN) {
        ctx->MDSNSKIP++;
/*      tda_out("skipped\n");        */     
        return(1);
    }

    w = ub[0] - lb[0];
    is = 0;
    for (i = 1; i < n; ++i) {
        tmp = ub[i] - lb[i];
        if (tmp > w) {
            w = tmp;
            is = i;
        }
    }
/*   tda_out("w=%g is=%d  \n",w,is);        */
          
    if (w > ctx->TOLFD) {
        a = lb[is];
        b = ub[is];
        c = (a + b) / 2.0;

        ub[is] = c;
        mdsxr_c2(ctx, n,lb,ub,x,y);

        ub[is] = b;
        lb[is] = c;

        mdsxr_c2(ctx, n,lb,ub,x,y);

        lb[is] = a;
    }
    else {
                              
        if (f < ctx->MDSFMIN) {
            tda_out("1-MDSFMIN=%g <- f=%g TOLBW=%g w=%g \n",ctx->MDSFMIN,f,ctx->TOLFD,w  );
            ctx->MDSFMIN = f;

            for (i = 0; i < n; ++i) {
                x[i] = lb[i];
                y[i] = ub[i];

                tda_out("i=%3d %20.16f %20.16f\n",i,lb[i],ub[i]);
            }
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mdsxr_c3()                                                               */
/*                                                                          */

int mdsxr_c3(TDAContext *ctx, int n,double *lb,double *ub,double *x,double *y)
{
    register int i;
    int is;
    double a,b,c,f,w,tmp;

    ctx->MDSNCALL++;
         
    mdsxr_f1(ctx, 1,n,lb,ub,&f,ctx->TOLFE);

tda_out("f=%g TOLFE=%g\n",f,ctx->TOLFE);
    if (f > ctx->TOLFE) {
tda_out("skipped\n");
        ctx->MDSNSKIP++;
        return(1);
    }
          


    w = ub[0] - lb[0];
    is = 0;
    for (i = 1; i < n; ++i) {
        tmp = ub[i] - lb[i];
        if (tmp > w) {
            w = tmp;
            is = i;
        }
    }


    tda_out("w=%g is=%d MDSFMIN=%g \n",w,is,ctx->MDSFMIN);  
          
    if (w >= ctx->MDSFMIN) {
        a = lb[is];
        b = ub[is];
        c = (a + b) / 2.0;

        ub[is] = c;
        mdsxr_c3(ctx, n,lb,ub,x,y);

        ub[is] = b;
        lb[is] = c;

        mdsxr_c3(ctx, n,lb,ub,x,y);

        lb[is] = a;
    }
    else if (w < ctx->MDSFMIN) {
tda_out("hier\n");             
        if (f < ctx->TOLFE) {
            tda_out("1-MDSFMIN=%g <- w=%g TOLBW=%g f=%g TOLFE=%g\n",ctx->MDSFMIN,w,ctx->TOLFD,f,ctx->TOLFE);
            ctx->MDSFMIN = w;

            for (i = 0; i < n; ++i) {
                x[i] = lb[i];
                y[i] = ub[i];

                tda_out("i=%3d %20.16f %20.16f\n",i,lb[i],ub[i]);
            }
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mdsr    Find axes in configuration.                                     */
/*                                                                          */
/*  mdsr(                                                                   */
/*      xv=...,         variables                                           */
/*      fmt=...,        print format for df option, def. 10.4               */
/*  ) = X,y;            configuration                                       */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int mdsr(TDAContext *ctx)
{
    register int i;
    int err,ix,iy,iv,nn,vz,vzmax = 0;
    double phi,sinphi,cosphi,r,phimax,rmax,xmin,xmax,ymin,ymax,xmean,ymean;
    double d,ya,yb;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Find axes in configuration. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,4,1))     /* get parameters */
        goto MDSRFin;

    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    if (ctx->PMNV != 2) {
        printf1(ctx, "Error: need two variables on right-hand side.\n");
        goto MDSRFin;
    }
    if (ctx->PM1NV < 1) {
        printf1(ctx, "Error: no xv variables.\n");
        goto MDSRFin;
    }
    ix = ctx->PMVIdx[0];
    iy = ctx->PMVIdx[1];
    printf1(ctx, "\nConfiguration: %s, %s\n",ctx->VName[ix],ctx->VName[iy]);

    if (alloc_acx(ctx, ctx->NOC + 1))  
        goto MDSRFin;
    if (alloc_acy(ctx, ctx->NOC + 1))  
        goto MDSRFin;

    nn = 1.0;
    if (ctx->NOC > 1)
        nn = (int)(ctx->NOC * (ctx->NOC - 1) / 2.0);

    xmin = xmax = get_data(ctx, ix,0);
    ymin = ymax = get_data(ctx, ix,1);
    xmean = xmin;       
    ymean = ymin;

    for (i = 1; i < ctx->NOC; ++i) {
        r = get_data(ctx, ix,i);
        xmean += r;
        xmin = dmin(ctx, xmin,r);
        xmax = dmax(ctx, xmax,r);
        r = get_data(ctx, iy,i);
        ymean += r;
        ymin = dmin(ctx, ymin,r);
        ymax = dmax(ctx, ymax,r);
    }
    xmean /= (double)ctx->NOC;
    ymean /= (double)ctx->NOC;

    printf1(ctx, "X minimum: ");
    rt_printf1_d(ctx, ctx->PMFmtS,xmin);
    printf1(ctx, "  maximum: ");
    rt_printf1_d(ctx, ctx->PMFmtS,xmax);
    printf1(ctx, "  mean: ");
    rt_printf1_d(ctx, ctx->PMFmtS,xmean);
    printf1(ctx, "\nY minimum: ");
    rt_printf1_d(ctx, ctx->PMFmtS,ymin);
    printf1(ctx, "  maximum: ");
    rt_printf1_d(ctx, ctx->PMFmtS,ymax);
    printf1(ctx, "  mean: ");
    rt_printf1_d(ctx, ctx->PMFmtS,ymean);
    newline(ctx);

    if (fabs(xmin - xmax) <= ctx->EPSI1)
        goto MDSRFin; 

    for (iv = 0; iv < ctx->PM1NV; ++iv) {
        printf1(ctx, "\nVariable: %s\n",ctx->VName[ctx->PM1VIdx[iv]]);

        for (i = 0; i < ctx->NOC; ++i)  
            ctx->AcY[i] = get_data(ctx, ctx->PM1VIdx[iv],i);

        rmax = phimax = 0.0;
        for (phi = 0.0; phi < Pi; phi += 0.01) {
            sinphi = sin(phi);
            cosphi = cos(phi);
            for (i = 0; i < ctx->NOC; ++i)  
                ctx->AcX[i] = cosphi * get_data(ctx, ix,i) + sinphi * get_data(ctx, iy,i);

            r = mdsr_f(ctx, ctx->NOC,ctx->AcX,ctx->AcY);
            if (r < 0.0) {
                vz = -1;
                r = -1.0;
            }
            else
                vz = 1;

            if (rmax < r) {
                rmax = r;
                phimax = phi;
                vzmax = vz;
            }
        }
        rmax /= nn; 
        printf1(ctx, "Optimal phi: %lg\n",phimax);
        printf1(ctx, "Rank correlation: %lg\n",rmax * (double)vzmax);

        phi = phimax;
        sinphi = sin(phi);
        cosphi = cos(phi);

        /***************************
        for (i = 0; i < NOC; ++i) {
            printf1(ctx, PMFmtS,cosphi * get_data(ctx, ix,i) + sinphi * get_data(ctx, iy,i));
            printf1(ctx, PMFmtS,-sinphi * get_data(ctx, ix,i) + cosphi * get_data(ctx, iy,i));
            newline(ctx);
        }
        *********************************/

        d = sinphi * (xmean - xmin);
        ya = ymean - d;
        d = sinphi * (xmax - xmean);
        yb = ymean + d;

        printf1(ctx, "Axis: (%lg,%lg) to (%lg,%lg)\n",xmin,ya,xmax,yb);
    }
    err = 0;

MDSRFin:
    p_clean(ctx);
    return(err);
}


/* ------------------------------------------------------------------------ */
/*  mdsr_f  return rank correlation                                         */
/*                                                                          */

double mdsr_f(TDAContext *ctx, int n,double *x,double *y)
{
    register int i,j;
    double r;                  

    r = 0.0;
    for (i = 1; i < n; ++i) {
        for (j = 0; j < i; ++j)  
            r += dsign(ctx, x[i] - x[j]) * dsign(ctx, y[i] - y[j]);
    }
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  rfit        Fit of binary relation(s).                                  */
/*                                                                          */
/*              rfit(                                                       */
/*                  gn=...,   graph number, def. 1                          */
/*                  rel=...,  properties of relation                        */
/*                  max=...,  upper limit to number of solutions, def.10    */
/*                  df=...,   additional output file                        */
/*                                                                          */
/*              ) = fname;                                                  */
/*                                                                          */
/*                                                                          */
/*  rel=... is a comma-separated list of integers:                          */
/*   1 reflexive                                                            */
/*   2 symmetric                                                            */
/*   3 antisymmetric                                                        */
/*   4 transitive                                                           */
/*   5 complete                                                             */
/*                                                                          */
/*  Solutions are written into the output file specified on the             */
/*  right-hand side.                                                        */  
/*                                                                          */
/*  If requested, the data for the linear programming problem is written    */
/*  into the file specified with the df parameter.                          */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */
       
int rfit(TDAContext *ctx)
{
    register int i,j,k,ii;
    int err,n,nh,nt,m,r,vmax,ncon;
    /* nn bounds the cleanup loop at the exit label, which is reached by
       goto before nn is assigned when parm() or a graph check fails.
       Uninitialised it walked AcI/AcIPtr out of bounds and freed
       whatever it found there.  Present in the 6.4 release. */
    int nn = 0;
    int *ap;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Fit of binary relations. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto RFITFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto RFITFin;
    if (gdd_tcheck(ctx, 0,0,2))          /* need valued graph */
        goto RFITFin;

    if (ctx->PMMax < 1)
        ctx->PMMax = 10;

    n = ctx->GD_NP;         
    nn = n * n;
    nh = (n * (n - 1)) / 2;

    printf1(ctx, "\nNumber of nodes: %d\n",n);
    printf1(ctx, "Entries in adjacency matrix: %d\n",nn);

    ncon = 0;                            /* number of constraints */
    printf1(ctx, "\nRequested properties:\n");
    if (ctx->PMREL[1] == 1) {
        printf1(ctx, "- reflexive\n");
        ncon += n;
    }
    if (ctx->PMREL[2] == 1) {
        printf1(ctx, "- symmetric\n");
        ncon += 2 * nh;
    }
    if (ctx->PMREL[3] == 1) {
        printf1(ctx, "- anti-symmetric\n");
        ncon += nh;
    }
    if (ctx->PMREL[4] == 1) {
        printf1(ctx, "- transitive\n");
        nt = 0;
        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= n; ++j) {
                if (j == i)
                    continue;
                for (k = 1; k <= n; ++k) { 
                    if (k == i || k == j)
                        continue;
                    nt++;
                }
            }
        }
        ncon += nt;
    }
    if (ctx->PMREL[5] == 1) {
        printf1(ctx, "- complete\n");
        ncon += nh;
    }
    printf1(ctx, "\nNumber of constraints: %d\n",ncon);

    m = ncon + 1;

    if (alloc_acr(ctx, m * nn + 1))  
        goto RFITFin;

    if (alloc_acs(ctx, m + 1))  
        goto RFITFin;

    k = 0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (i == j)
                r = 1;
            else
                r = (int)gdd_adj(ctx, i,j,ctx->PMGN);                   
            ctx->AcR[++k] = 2 * r - 1;                                    
        }
    }

    /* add constraints */

    ii = 2;
    if (ctx->PMREL[1] == 1) {               /* reflexive */
        for (i = 1; i <= n; ++i) {                                         
            ctx->AcR[(ii - 1) * nn + (i - 1) * n + i] = 1;
            ctx->AcS[ii] = 1;
            ii++;
        }
    }
    if (ctx->PMREL[2] == 1) {               /* symmetric */
        for (i = 1; i < n; ++i) {                                          
            for (j = i + 1; j <= n; ++j) {
                ctx->AcR[(ii - 1) * nn + (i - 1) * n + j] =  1;
                ctx->AcR[(ii - 1) * nn + (j - 1) * n + i] = -1;
                ctx->AcS[ii] = 0;
                ii++;
                ctx->AcR[(ii - 1) * nn + (i - 1) * n + j] = -1;
                ctx->AcR[(ii - 1) * nn + (j - 1) * n + i] =  1;
                ctx->AcS[ii] = 0;
                ii++;
            }
        }
    }
    if (ctx->PMREL[3] == 1) {               /* anti-symmetric */
        for (i = 1; i < n; ++i) {                                          
            for (j = i + 1; j <= n; ++j) {
                ctx->AcR[(ii - 1) * nn + (i - 1) * n + j] = -1;
                ctx->AcR[(ii - 1) * nn + (j - 1) * n + i] = -1;
                ctx->AcS[ii] = -1;
                ii++;
            }
        }
    }
    if (ctx->PMREL[4] == 1) {               /* transitive */
        for (i = 1; i <= n; ++i) {                                    
            for (j = 1; j <= n; ++j) {
                if (j == i)
                    continue;
                for (k = 1; k <= n; ++k) { 
                    if (k == i || k == j)
                        continue;
                    ctx->AcR[(ii - 1) * nn + (i - 1) * n + k] =  1;
                    ctx->AcR[(ii - 1) * nn + (i - 1) * n + j] = -1;
                    ctx->AcR[(ii - 1) * nn + (j - 1) * n + k] = -1;
                    ctx->AcS[ii] = -1;
                    ii++;
                }
            }
        }
    } 
    if (ctx->PMREL[5] == 1) {               /* complete */
        for (i = 1; i < n; ++i) {                                        
            for (j = i + 1; j <= n; ++j) {
                ctx->AcR[(ii - 1) * nn + (i - 1) * n + j] = 1;
                ctx->AcR[(ii - 1) * nn + (j - 1) * n + i] = 1;
                ctx->AcS[ii] = 1;
                ii++;
            }
        }
    }
    if (ctx->PMF1Def) {
        for (i = 0; i < m; ++i) {
            for (j = 1; j <= nn; ++j) 
                fprintf(ctx->PMF1d,"%2d ",ctx->AcR[i * nn + j]);
            fprintf(ctx->PMF1d,"%3d\n",ctx->AcS[i + 1]);
        }
        printf1(ctx, "%d records written to: %s\n",m,ctx->PMF1dName);
    }
    newline(ctx);

    if (alloc_acn(ctx, ctx->PMMax * nn + 1))
        goto RFITFin;

    if (alloc_acm(ctx, 4 * m + 6 * nn + m * nn + 2))
        goto RFITFin;

    r = lpi(ctx, nn,m,ctx->PMMax,ctx->AcR,ctx->AcS,ctx->AcN,ctx->AcM,&vmax);

    if (r < 1) {
        switch (r) {
            case -1:    printf1(ctx, "Error: inconsistent constraints.\n");
                        break;
            case -2:    printf1(ctx, "Error: exceeded maximal number of solutions.\n");
                        break;
            default:    printf1(ctx, "Error (%d).\n",r);
                        break;
        }   
        goto RFITFin;
    }
    printf1(ctx, "Value: %d\nNumber of solutions: %d\n",vmax,r);

    nt = 0;
    for (ii = 0; ii < r; ++ii) {

        for (j = 1; j <= nn; ++j)
            fprintf(ctx->PMFd,"%d ",ctx->AcN[ii * nn + j]);
        fprintf(ctx->PMFd,"\n");
        nt++;
        ap = ctx->AcN + ii * nn;
        k = 1;
        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= n; ++j)
                fprintf(ctx->PMFd,"%d ",ap[k++]);
            fprintf(ctx->PMFd,"\n");
            nt++;
        }
        fprintf(ctx->PMFd,"\n");
        nt++;
    }
    printf1(ctx, "%d records written to: %s\n",nt,ctx->PMFdName);
    err = 0;

RFITFin:
    p_clean(ctx);
    return(err);
}

/* -###-------------------------------------------------------------------- */
/*  rfit1       Fit of binary relation(s).                                  */
/*                                                                          */
/*              rfit(                                                       */
/*                  gn=...,   graph number, def. 1                          */
/*                  rel=...,  properties of relation                        */
/*                  max=...,  upper limit to number of solutions, def.10    */
/*                  df=...,   additional output file                        */
/*                                                                          */
/*              ) = fname;                                                  */
/*                                                                          */
/*                                                                          */
/*  rel=... is a comma-separated list of integers:                          */
/*   1 reflexive                                                            */
/*   2 symmetric                                                            */
/*   3 antisymmetric                                                        */
/*   4 transitive                                                           */
/*   5 complete                                                             */
/*                                                                          */
/*  Solutions are written into the output file specified on the             */
/*  right-hand side.                                                        */  
/*                                                                          */
/*  If requested, the data for the linear programming problem is written    */
/*  into the file specified with the df parameter.                          */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */
       
int rfit1(TDAContext *ctx)
{
    register int i,j,k,ii;
    int err,n,nh,nt,m,r,vmax,ncon;
    /* nn bounds the cleanup loop at the exit label, which is reached by
       goto before nn is assigned when parm() or a graph check fails.
       Uninitialised it walked AcI/AcIPtr out of bounds and freed
       whatever it found there.  Present in the 6.4 release. */
    int nn = 0;
    int *ap;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Fit of binary relations. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,1,1))       /* get parameters */
        goto RFITFin;

    if (gdd_check(ctx, ctx->PMGN,1))
        goto RFITFin;
    if (gdd_tcheck(ctx, 0,0,2))          /* need valued graph */
        goto RFITFin;

    if (ctx->PMMax < 1)
        ctx->PMMax = 10;

    n = ctx->GD_NP;         
    nn = n * n;
    nh = (n * (n - 1)) / 2;

    printf1(ctx, "\nNumber of nodes: %d\n",n);
    printf1(ctx, "Entries in adjacency matrix: %d\n",nn);

    ncon = 0;                            /* number of constraints */
    printf1(ctx, "\nRequested properties:\n");
    if (ctx->PMREL[1] == 1) {
        printf1(ctx, "- reflexive\n");
        ncon += n;
    }
    if (ctx->PMREL[2] == 1) {
        printf1(ctx, "- symmetric\n");
        ncon += 2 * nh;
    }
    if (ctx->PMREL[3] == 1) {
        printf1(ctx, "- anti-symmetric\n");
        ncon += nh;
    }
    if (ctx->PMREL[4] == 1) {
        printf1(ctx, "- transitive\n");
        nt = 0;
        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= n; ++j) {
                if (j == i)
                    continue;
                for (k = 1; k <= n; ++k) { 
                    if (k == i || k == j)
                        continue;
                    nt++;
                }
            }
        }
        ncon += nt;
    }
    if (ctx->PMREL[5] == 1) {
        printf1(ctx, "- complete\n");
        ncon += nh;
    }
    printf1(ctx, "\nNumber of constraints: %d\n",ncon);

    m = ncon + 1;

    if (alloc_acr(ctx, m * nn + 1))  
        goto RFITFin;

    if (alloc_acs(ctx, m + 1))  
        goto RFITFin;

    k = 0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (i == j)
                r = 1;
            else
                r = (int)gdd_adj(ctx, i,j,ctx->PMGN);                   
            ctx->AcR[++k] = 2 * r - 1;                                    
        }
    }

    /* add constraints */

    ii = 2;
    if (ctx->PMREL[1] == 1) {               /* reflexive */
        for (i = 1; i <= n; ++i) {                                         
            ctx->AcR[(ii - 1) * nn + (i - 1) * n + i] = 1;
            ctx->AcS[ii] = 1;
            ii++;
        }
    }
    if (ctx->PMREL[2] == 1) {               /* symmetric */
        for (i = 1; i < n; ++i) {                                          
            for (j = i + 1; j <= n; ++j) {
                ctx->AcR[(ii - 1) * nn + (i - 1) * n + j] =  1;
                ctx->AcR[(ii - 1) * nn + (j - 1) * n + i] = -1;
                ctx->AcS[ii] = 0;
                ii++;
                ctx->AcR[(ii - 1) * nn + (i - 1) * n + j] = -1;
                ctx->AcR[(ii - 1) * nn + (j - 1) * n + i] =  1;
                ctx->AcS[ii] = 0;
                ii++;
            }
        }
    }
    if (ctx->PMREL[3] == 1) {               /* anti-symmetric */
        for (i = 1; i < n; ++i) {                                          
            for (j = i + 1; j <= n; ++j) {
                ctx->AcR[(ii - 1) * nn + (i - 1) * n + j] = -1;
                ctx->AcR[(ii - 1) * nn + (j - 1) * n + i] = -1;
                ctx->AcS[ii] = -1;
                ii++;
            }
        }
    }
    if (ctx->PMREL[4] == 1) {               /* transitive */
        for (i = 1; i <= n; ++i) {                                    
            for (j = 1; j <= n; ++j) {
                if (j == i)
                    continue;
                for (k = 1; k <= n; ++k) { 
                    if (k == i || k == j)
                        continue;
                    ctx->AcR[(ii - 1) * nn + (i - 1) * n + k] =  1;
                    ctx->AcR[(ii - 1) * nn + (i - 1) * n + j] = -1;
                    ctx->AcR[(ii - 1) * nn + (j - 1) * n + k] = -1;
                    ctx->AcS[ii] = -1;
                    ii++;
                }
            }
        }
    } 
    if (ctx->PMREL[5] == 1) {               /* complete */
        for (i = 1; i < n; ++i) {                                        
            for (j = i + 1; j <= n; ++j) {
                ctx->AcR[(ii - 1) * nn + (i - 1) * n + j] = 1;
                ctx->AcR[(ii - 1) * nn + (j - 1) * n + i] = 1;
                ctx->AcS[ii] = 1;
                ii++;
            }
        }
    }
    if (ctx->PMF1Def) {
        for (i = 0; i < m; ++i) {
            for (j = 1; j <= nn; ++j) 
                fprintf(ctx->PMF1d,"%2d ",ctx->AcR[i * nn + j]);
            fprintf(ctx->PMF1d,"%3d\n",ctx->AcS[i + 1]);
        }
        printf1(ctx, "%d records written to: %s\n",m,ctx->PMF1dName);
    }
    newline(ctx);

    if (alloc_ack(ctx, nn + 1))
        goto RFITFin;
    if (alloc_acl(ctx, nn + 1))
        goto RFITFin;
    if (alloc_aci(ctx, nn + 1))
        goto RFITFin;
    if (alloc_acj(ctx, nn + 1))
        goto RFITFin;
    if (alloc_aciptr(ctx, nn + 1))
        goto RFITFin;
    if (alloc_acjptr(ctx, nn + 1))
        goto RFITFin;

    for (j = 1; j <= nn; ++j) {
        tda_out("%5d : ",j);
        nt = 0;
        for (i = 0; i < m; ++i) {
            if (ctx->AcR[i * nn + j] != 0)
                nt++;
        }
        tda_out("%d\n",nt);

        ctx->AcK[j] = ctx->AcL[j] = nt;

        if (!(ctx->AcIPtr[j] = (int *)calloc((size_t)(2) * (size_t)(nt) + 2,sizeof(int)))) {
            p_err(ctx, -2,1);
            goto RFITFin;
        }
        ctx->AcI[j] = 2 * nt + 2;
        memrq(ctx, 2 * nt + 2,sizeof(int));
        if (!(ctx->AcJPtr[j] = (int *)calloc((size_t)(2) * (size_t)(nt) + 2,sizeof(int)))) {
            p_err(ctx, -2,1);
            goto RFITFin;
        }
        ctx->AcJ[j] = 2 * nt + 2;
        memrq(ctx, 2 * nt + 2,sizeof(int));

        ii = nt + 1;
        for (i = 0; i < m; ++i) {
            if ((k = ctx->AcR[i * nn + j]) != 0) {
                ctx->AcIPtr[j][ii] = k;
                ctx->AcJPtr[j][ii] = i + 1;
                ii++;
            }
        }
    }
    alloc_acr(ctx, 0);       /* no longer required */

    if (alloc_acn(ctx, ctx->PMMax * nn + 1))
        goto RFITFin;

    if (alloc_acm(ctx, 4 * m + 6 * nn + 2000            ))
        goto RFITFin;

/**
    r = lpia(ctx, nn,m,PMMax,AcR,AcS,AcN,AcM,&vmax);               
    r = 0;

int lpia(TDAContext *ctx, int n,int m,int nest,int *an,int *an0,int **av,int **ar, int *b0,int *opts,int *itmp,int *vmax);


**/
  
    r = lpia(ctx, nn,m,ctx->PMMax,ctx->AcK,ctx->AcL,ctx->AcIPtr,ctx->AcJPtr,ctx->AcS,ctx->AcN,ctx->AcM,&vmax);
  
     

    if (r < 1) {
        switch (r) {
            case -1:    printf1(ctx, "Error: inconsistent constraints.\n");
                        break;
            case -2:    printf1(ctx, "Error: exceeded maximal number of solutions.\n");
                        break;
            default:    printf1(ctx, "Error (%d).\n",r);
                        break;
        }   
        goto RFITFin;
    }
    printf1(ctx, "Value: %d\nNumber of solutions: %d\n",vmax,r);

    nt = 0;
    for (ii = 0; ii < r; ++ii) {

        for (j = 1; j <= nn; ++j)
            fprintf(ctx->PMFd,"%d ",ctx->AcN[ii * nn + j]);
        fprintf(ctx->PMFd,"\n");
        nt++;
        ap = ctx->AcN + ii * nn;
        k = 1;
        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= n; ++j)
                fprintf(ctx->PMFd,"%d ",ap[k++]);
            fprintf(ctx->PMFd,"\n");
            nt++;
        }
        fprintf(ctx->PMFd,"\n");
        nt++;
    }
    printf1(ctx, "%d records written to: %s\n",nt,ctx->PMFdName);
    err = 0;

RFITFin:
    for (j = 1; j <= nn; ++j) {
        if (ctx->AcI[j] > 0) {
            free((char *)ctx->AcIPtr[j]);
            memrq(ctx, -ctx->AcI[j],sizeof(int));
            ctx->AcI[j] = 0;        
        }
        if (ctx->AcJ[j] > 0) {
            free((char *)ctx->AcJPtr[j]);
            memrq(ctx, -ctx->AcJ[j],sizeof(int));
            ctx->AcI[j] = 0;        
        }
    }
    p_clean(ctx);
    return(err);
}








