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

/*  functions in t_mds.c */

int dmet(void); 
int dmet1(void); 
double dmet_c(int n,float *d);
void m_spath(int n,float *d,float *a);

int gpro(void);
void gpro_pcf(int n,double *x);
int gpro_sval(int n,double *x,int *idx);
double gpro_fun(int n,double *x);
void gpro_df(int n,double *x);
int gpro_min(int n,double *par,double *grad,double *step);

int uds(void);
void uds_check(int n,double *a,double *b);
void uds_ncheck(int n,int *a,int *b);
int uds_1(void);  
int uds_0(int n,int k);
double uds_get_dip(int n,int i,int *p);
double uds_get_dpi(int n,int i,int *p);
double uds_get_f(int n,int *p);
int uds_r(int n,int k);
int uds_r0(int n,int k,double ck,double ckk0,double cl);
int uds_r1(int n,int k,int l,double ck,double cl,double cll0);
int uds_qa1(void);  
int uds_qa2(void);  

int sga(void);
int sga1(int m,int *r);
int hpo(int n,int m,int *x);
int hpo1(int m,int *x,int *y);

int unf(void);
void unf1(int n,int m,int *p,int *r,double *x,int *s,int *rr,int prn,
    short *pp,int wf,double *w);
int unf2(int m,int *r,int *s,int *rr);

int mdsc(void);
int mdsm(void);
int mds4_fn(void);


int mdsn(void);
void mdsn_dis(double *dis);
int mdsn_norm(void);
int mdsn_grad(double s,double t,double stress);   

int mdsn1(void);
void mdsn1_dis(int n,double *par,double *dis);
void mdsn1_proj(int opt,int nn,double *dis,double *tmp,double *proj);
int mdsn1_fn(int n,int nn,double *par,double *fval,double *grad,double *dis,double *disp);

double mds_edis(double x1,double y1,double x2,double y2);
double mds_elen(double x,double y);

int mdsx(void);
int mds_m(int ip,int ip0,int ip1);
int mds3_fn(void);
int mds_r(int n,float *d,double *x,double *y,double *r,int *p,double *ds,
    int *is,int *js,int *ptr);
void mds_p(int n,double *x,double *y,double *r);
void mds_p1(int n,double *x,double *y);
int mdsxr_m(int narg,int nbmax,double *par,double *lb,double *ub,double *par0);
int mdsxr_f(int opt,int n,double *xl,double *xh,double *f);
double mdsxr_fm(void);
int mdsxr_f1(int opt,int n,double *xl,double *xh,double *f,double fminp);
double mdsxr_dm(double ai,double bi,double aj,double bj);
int mdsxr_c2(int n,double *lb,double *ub,double *x,double *y);
int mdsxr_c3(int n,double *lb,double *ub,double *x,double *y);

int mdsr(void);
double mdsr_f(int n,double *x,double *y);

int rfit(void);
int rfit1(void);

/* ------------------------------------------------------------------------ */

double MDSFN = 0.0;
int MDSIter = 0;
int MDSFCalls = 0;
double MDSTOLS = 1.e-6;
double MDSTOLG = 1.e-4;

int MDSTYP = 0;                 /* 0 = minmax, 1 = least squares            */
int MDSRN = 0;                  /* number of points                         */
int MDSND = 2;                  /* number of dimensions                     */
int MDSNP = 0;                  /* number of parameters                     */
int MDSRNN = 0;                 /* number of distances                      */
double *MDSX;                   /* configuration, parameters                */
double *MDSXB;                  /* best configuration                       */
double *MDSG;                   /* gradient                                 */
double *MDSGP;                  /* previous gradient                        */
double *MDSDIS;                 /* distances                                */
double *MDSDIS1;                /* distances                                */
double *MDSDIS2;                /* distances                                */
double *MDSDIS3;                /* distances                                */
int *MDSPTR;                    /* pointer for sort of MDSDIS               */
int *MDSPTRI;                   /* inverse pointer for sort of MDSDIS       */
float *MDSRD;                   /* distance matrix                          */
int MDSRNI = 0;                 /* number of points in MDSRI[]              */
int *MDSRI;                     /* list of points                           */
int *MDSRF;                     /* flags already used points                */
int MDSRIP = 0;                 /* current point                            */
int MDSRIP0 = 0;                /* first fixed point                        */
int MDSRIP1 = 0;                /* second fixed point                       */
int MDSMSTEP = 0;               /* 0/1 step in optimization                 */
int MDSNBOX = 0;                /* count critical boxes in verify step      */
int MDSNCALL = 0;               /* number of function calls in verify step  */
int MDSNSKIP = 0;           
double MDSFMIN = 0.0;           /* global minimum for verify step           */

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
       
int dmet(void)
{
    register int i,j,k;
    int err,n;
    double a,dij,dm;              

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Change distance matrix into metric. Current memory: %d bytes.\n",MemReq);

    TOLFD = 1.0e-4;                 /* default box width */

    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto DMETFin;

    if (gdd_check(PMGN,1))
        goto DMETFin;
    if (gdd_tcheck(0,1,2))          /* need undirected valued graph */
        goto DMETFin;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    n = GD_NP;         
    if (n < 3) {
        printf1("Nothing done for less than 3 points.\n");
        err = 0;
        goto DMETFin;
    }
    if (alloc_acx(n * n + 1))  
        goto DMETFin;

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (i == j)
                AcX[i * n + j] = 0.0;
            else
                AcX[i * n + j] = gdd_adj(i,j,PMGN);
        }
    }
    fprintf(PMFd,"Original distance matrix\n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j)
            fprintf(PMFd,PMFmtS,AcX[i * n + j]);
        fprintf(PMFd,"\n");
    }
    dm = 0.0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (i == j)
                continue;
            dij = AcX[i * n + j];
            for (k = 0; k < n; ++k) {
                if (k == i || k == j)
                    continue;
                a = dij - AcX[i * n + k] - AcX[k * n + j];
                if (dm < a) {
printf("i=%d j=%d k=%d a=%lg dm=%lg ",i+1,j+1,k+1,a,dm);
                    dm = a;
printf(" new dm=%lg\n",dm); 

                }
            }
        }
    }
    printf1("Additive constant: %lg\n",dm);

    fprintf(PMFd,"Modified distance matrix\n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            a = AcX[i * n + j];
            if (i != j)       
                a += dm;
            fprintf(PMFd,PMFmtS,a);
        }
        fprintf(PMFd,"\n");
    }
    err = 0;

DMETFin:
    p_clean();
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
       
int dmet1(void)
{
    register int i,j,k;
    int err,n,iter;
    double d,dij,dik,djk,dev; 

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Change distance matrix into metric. Current memory: %d bytes.\n",MemReq);

    TOLFD = 1.0e-4;                 /* default box width */

    if (parm(CmdBuf + 5,1,1))       /* get parameters */
        goto DMET1Fin;

    if (gdd_check(PMGN,1))
        goto DMET1Fin;
    if (gdd_tcheck(0,1,2))          /* need undirected valued graph */
        goto DMET1Fin;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);
    if (MxItFlg == 0)
        MxIter = 20;

    printf1("Maximal number of iterations: %d\n",MxIter);
    printf1("Tolerance for maximal deviation: %g\n\n",TOLFD);

    n = GD_NP;         
    if (n < 3) {
        printf1("Nothing done for less than 3 points.\n");
        err = 0;
        goto DMET1Fin;
    }
    if (alloc_acxf(n * n + 1))  
        goto DMET1Fin;
    if (alloc_acyf(n * n + 1))  
        goto DMET1Fin;

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (i == j)
                AcXF[i * n + j] = 0.0;
            else
                AcXF[i * n + j] = (float)gdd_adj(i,j,PMGN);
        }
    }
    fprintf(PMFd,"Original distance matrix\n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j)
            fprintf(PMFd,PMFmtS,AcXF[i * n + j]);
        fprintf(PMFd,"\n");
    }
    dev = dmet_c(n,AcXF);     

    printf1("Maximal deviation in input data: %g\n",dev);
    if (dev <= TOLFD) 
        goto DMET1Fin;

    for (iter = 1; iter <= MxIter; ++iter) {

        m_spath(n,AcXF,AcYF);       /* find shortest paths */

        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; ++j) {
                if (i != j) {
                    k = i * n + j;
                    AcXF[k] = dmin(AcXF[k],AcYF[k]);
                }
            }
        }
        dev = dmet_c(n,AcXF);
/*      printf1("dev for modified X = %g\n",dev);   */

        if (dev <= TOLFD)
            break;
    }
    if (iter > MxIter)
        iter--;

    printf1("Number of iterations performed: %d\n",iter);
    printf1("Remaining deviation: %g\n\n",dev);

    fprintf(PMFd,"Modified distance matrix\n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j)
            fprintf(PMFd,PMFmtS,AcXF[i * n + j]);
        fprintf(PMFd,"\n");
    }
    err = 0;

DMET1Fin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  dmet_c(n,d)     Check d[] for metric. Return max deviation.             */
       
double dmet_c(int n,float *d)
{
    register int i,j,k;
    double dij,dik,dkj,dev;
                  
    dev = 0.0;
    for (i = 1; i < n; ++i) {
        for (j = 0; j < i; ++j) {
            dij = d[i * n + j];
            for (k = 0; k < n; ++k) {
                if (k == i || k == j)
                    continue;
                dik = d[i * n + k];
                dkj = d[k * n + j];
                dev = dmax(dev,dij - dik - dkj);
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

void m_spath(int n,float *d,float *a)
{
    register int i,j,k,ni,nj;
    float s,max;

    max = 0.0; 
    for (i = 1; i < n; ++i) {
        for (j = 0; j < i; ++j) 
            max = dmax(max,d[i * n + j]);
    }
    max *= 10.0;

    for (i = 0; i < n; ++i) {
        ni = i * n;
        for (j = 0; j < n; ++j) {
            a[ni] = max;
            if (i != j) {
                s = d[i * n + j];                                        
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

int gpro(void)
{
    register int i;
    int err,r,n,nn,nrec;
    double x,y;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Projection of proximities. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto GPROFin;

    if (gdd_check(PMGN,1))
        goto GPROFin;
    if (gdd_tcheck(0,1,2))          /* need undirected valued graph */
        goto GPROFin;

    if (GD_NP < 3) {
        printf1("Error: need at least three nodes.\n");
        goto GPROFin;
    }

    if (PMALG != 1)
        PMALG = 1;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);
    if (MxItFlg == 0)
        MxIter = 50;

    printf1("\nAlgorithm %d: ",PMALG);
    if (PMALG == 1) {
        printf1("direct search.\n");
        printf1("Starting value of step length: %g\n",SLen);
        printf1("Step length reduction factor: %g\n",SRed);
        printf1("Minimum value of step length: %g\n",TOLS);
    }
    printf1("Max number of iterations: %d\n\n",MxIter);

    if (alloc_acx(2 * GD_NP + 1))       /* coordinates */
        goto GPROFin;
    if (alloc_acy(2 * GD_NP + 1))      
        goto GPROFin;
    if (alloc_actmp(2 * GD_NP + 1))   
        goto GPROFin;
    if (alloc_acn(GD_NP + 1))   
        goto GPROFin;

    MDSIter = MDSFCalls = 0;
    n = GD_NP - 1;

    /* calculate starting values */

    if (gpro_sval(n,AcX,AcN)) {
        printf1("Error: graph not connected.\n");
        goto GPROFin;
    }
    nn = 2 * n;
    r = gpro_min(nn,AcX,AcY,AcTmp);

    printf1("Performed iterations: %d (%d function calls)\n",MDSIter,MDSFCalls);
    printf1("Final stress: %g\n\n",MDSFN);

    nrec = 0;
    for (i = 0; i <= n; ++i) {
        if (i == 0)
            x = y = 0.0;
        else {
            x = AcX[2 * i - 1];
            y = AcX[2 * i];
        }
        fprintf(PMFd,PMNFmtS,i + 1);
        fprintf(PMFd,PMFmtS,x);
        fprintf(PMFd,PMFmtS,y);
        fprintf(PMFd,"\n");
        nrec++;
    }
    printf1("%d records written to: %s\n",nrec,PMFdName);

    if (PMF1Def)            /* write estimated distances to PMF1d */
        gpro_df(n,AcX);

    if (PMPCFDef)           /* create cf file for plot of configuration */
        gpro_pcf(n,AcX);

    err = 0;

GPROFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gpro_pcf    Create a command file for plotting the configuration.       */

void gpro_pcf(int n,double *x)
{
    register int i;
    double xmin,xmax,ymin,ymax,xd,yd;

    xmin = ymin = xmax = ymax = 0.0;
    for (i = 1; i <= n; ++i) {
        xmin = dmin(xmin,x[2 * i - 1]);
        xmax = dmax(xmax,x[2 * i - 1]);
        ymin = dmin(ymin,x[2 * i]);
        ymax = dmax(ymax,x[2 * i]);
    }
    xd = (xmax - xmin) / 10.0;
    yd = (ymax - ymin) / 10.0;

    fprintf(PMPCFd,"psfile = %s.ps;\n",PMPCFName);
    fprintf(PMPCFd,"psetup(pxa=%g,%g,pya=%g,%g);\n",
                      xmin - xd,xmax + xd,ymin - yd,ymax + yd);

    fprintf(PMPCFd,"plg(\n  nmax=%d,\n",n + 1);
    for (i = 0; i <= n; ++i) {
        if (i == 0)  
            xd = yd = 0.0;
        else {
            xd = x[2 * i - 1];
            yd = x[2 * i];
        }
        fprintf(PMPCFd,"  node(n=%d)=%g,%g,\n",i + 1,xd,yd);
    }
    fprintf(PMPCFd,");\n");
    printf1("Output file for plot of configuration: %s\n",PMPCFName);
}

/* ------------------------------------------------------------------------ */
/*  gpro_sval(n,x,idx)  Calcuate starting values for function minimization. */
/*                      Return 0 if ok, -1 if graph not connected.          */

int gpro_sval(int n,double *x,int *idx)
{
    register int i,j;
    int err,fnd;
    double dij,xi,xj,yi,b,c,tmp;

    err = -1;

    for (i = 1; i <= n; ++i)   
        idx[i] = 0;

    x[0] = 0.0;
    idx[0] = 1;
    fnd = 0;
    for (j = 1; j <= n; ++j) {
        dij = gdd_adj(0,j,PMGN);
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

    for (i = 0; i < n; ++i) {
        if (idx[0] == 0)
            continue;
        for (j = i + 1; j <= n; ++j) {
            if (idx[j])
                continue;
            dij = gdd_adj(i,j,PMGN);
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
                xj = random1();
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

double gpro_fun(int n,double *x)
{
    register int i,j,k;
    int m;
    double dij,aij,xij,yij,f,f1;    

    MDSFCalls++;
    f1 = f = 0.0;
    m = n / 2;
    for (i = 0; i < m; ++i) {
        for (j = i + 1; j <= m; ++j) {
            dij = gdd_adj(i,j,PMGN);
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

void gpro_df(int n,double *x)
{
    register int i,j,nrec;
    double dij,aij,xij,yij;

    nrec = 0;
    for (i = 0; i < n; ++i) {
        for (j = i + 1; j <= n; ++j) {
            dij = gdd_adj(i,j,PMGN);
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

            fprintf(PMF1d,PMNFmtS,i + 1);
            fprintf(PMF1d,PMNFmtS,j + 1);
            fprintf(PMF1d,PMFmtS,dij);
            fprintf(PMF1d,PMFmtS,aij);
            fprintf(PMF1d,"\n");
            nrec++;
        }
    }
    printf1("%d records written to: %s\n",nrec,PMF1dName);
}

/* ------------------------------------------------------------------------ */
/*  gpro_min(n,par,par1,step)                                               */
/*                                                                          */
/*      Function minimization with a direct search method, used by gpro.    */
/*                                                                          */
/*      Return:     0 if successful.                                        */
/*                 -1 if exceeded max number of iterations.                 */

int gpro_min(int n,double *par,double *par1,double *step)
{
    register int k;
    double delta,fm,fm1,theta,tmp;

    if (SILENTFlg < 2)
        printfe("\n  Iter    Function Value       Step Length     Par Change   FCall\n");

    delta = SLen;
    for (k = 1; k <= n; ++k) {
        tmp = fabs(par[k]);
        if (tmp < EPSI1)
            step[k] = delta;
        else
            step[k] = delta * tmp;
    }
    MDSFN = gpro_fun(n,par);

FM1L1:  
    fm = MDSFN;
    for (k = 1; k <= n; ++k)
        par1[k] = par[k];

    for (k = 1; k <= n; ++k) {
        par1[k] += step[k];
        fm1 = gpro_fun(n,par1);

        if (fm1 < fm) 
            fm = fm1;
        else {
            step[k] = -step[k];
            par1[k] += 2.0 * step[k];

            fm1 = gpro_fun(n,par1);

            if (fm1 < fm)
                fm = fm1;
            else
                par1[k] -= step[k];
        }
    }
    if (fm < MDSFN) {

FM1L2:
        for (k = 1; k <= n; ++k) {
            if (par1[k] > par[k] && step[k] <  0.0 ||
                                            par1[k] <= par[k] && step[k] >= 0.0)   
                step[k] = -step[k];
            theta = par[k];
            par[k] = par1[k];
            par1[k] = 2.0 * par1[k] - theta;
        }
        MDSFN = fm;
        fm = fm1 = gpro_fun(n,par1);

        for (k = 1; k <= n; ++k) {
            par1[k] += step[k];

            fm1 = gpro_fun(n,par1);

            if (fm1 < fm) 
                fm = fm1;
            else {
                step[k] = -step[k];
                par1[k] += 2.0 * step[k];

                fm1 = gpro_fun(n,par1);
                if (fm1 < fm)
                    fm = fm1;
                else
                    par1[k] -= step[k];
            }
        }
        if (fm >= MDSFN)
            goto FM1L1;

        for (k = 1; k <= n; ++k) {
            if (fabs(par1[k] - par[k]) > 0.5 * fabs(step[k])) 
                goto FM1L2;
        }
    }
    MDSIter++;

    if (SILENTFlg < 2) {
        printfe("  %3d  %20.13e %17.10e ",MDSIter,MDSFN,delta);
        printfe("        --   %6d\n",MDSFCalls);
    }
    if (delta > TOLS) {
        if (MDSIter >= MxIter) {
            LConv = 2;
            goto FM1Fin;
        }
        delta *= SRed;
        for (k = 1; k <= n; ++k)
            step[k] *= SRed;
        goto FM1L1;
    }
    else  
        LConv = 0;

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

int uds(void)
{
    register int i,j;
    int err,n;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Uni-dimensional scaling. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto UDSFin;

    if (gdd_check(PMGN,1))
        goto UDSFin;
    if (gdd_tcheck(0,1,2))          /* need undirected valued graph */
        goto UDSFin;

    if (GD_NP < 3) {
        printf1("Error: need at least three nodes.\n");
        goto UDSFin;
    }
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    n = GD_NP;     /* number of points */
    printf1("Number of points: %d\n",n);              

    if (PMALG < 1 || PMALG > 3)
        PMALG = 1;

    if (PMALG == 1) {
        printf1("\nAlgorithm 1: Branch and bound (Defays 1978).\n");
        if (n >= 100) {
            printf1("Current maximum: 100 points.\n");
            goto UDSFin;
        }
        uds_1();   
    } 
    else if (PMALG == 2) {
        printf1("\nAlgorithm 2: approx. quadratic assignment.\n");
        if (uds_qa1())  
            goto UDSFin;

    }
    else {
        printf1("\nAlgorithm 3: approx. quadratic assignment.\n");

        if (MxIter < 1 || MxItFlg == 0)
            MxIter = 100;
        if (PMIDFA <= 0.0 || PMIDFA > 1.0)
            PMIDFA = 0.25;
        if (PMIDFB <= 0.0 || PMIDFB > 1.0)
            PMIDFB = 0.5;
        printf1("Max number of iterations: %d (alpha %lg, beta %lg).\n\n",
            MxIter,PMIDFA,PMIDFB);

        if (uds_qa2())  
            goto UDSFin;
    }
    err = 0;

UDSFin:       
    p_clean();
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  uds_check(n,a,b)                                                        */
/*                                                                          */
/*  Check whether distances in a and b have same order.                     */

void uds_check(int n,double *a,double *b)
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
                    printf("i=%3d j=%3d k=%d l=%d aij=%6.1lf akl=%6.1lf bij=%6.1lf bkl=%6.1lf  %6d %6d\n",
                    i,j,k,l,aij,akl,bij,bkl,d0,d1);
                    ***************************/
                }
            }
        }
    }
    printf1("%d of %d compares are inconsistent.\n",d1,d0);

}

/* -##--------------------------------------------------------------------- */
/*  uds_ncheck(n,a,b)                                                       */
/*                                                                          */
/*  Check whether distances in a and b have same order.                     */

void uds_ncheck(int n,int *a,int *b)
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
                    printf("i=%3d j=%3d k=%d l=%d aij=%6.1lf akl=%6.1lf bij=%6.1lf bkl=%6.1lf  %6d %6d\n",
                    i,j,k,l,aij,akl,bij,bkl,d0,d1);
                    ***************************/
                }
            }
        }
    }
    printf1("%d of %d compares are inconsistent.\n",d1,d0);

}

/* -##--------------------------------------------------------------------- */
/*  uds_1()     uni-dimensional scaling with branch and bound.              */
/*                                                                          */
/*  return  0 if successful,                                                */
/*         -1 if insuff. memory.                                            */
/*         -2 if undef. elements in distance matrix.                        */

double UDSDPP = 0.0;
double UDS_F0 = 0.0;        /* function value */
int UDS_E0[100];            /* currently best permutation */
int UDS_PI[100];            /* working permutation */
int UDS_IF[100];
double UDS_X[100];          /* x values */  
double UDS_D[10000];        /* distance matrix */

int uds_1(void)   
{
    register int i,j;
    int err,n;
    double d,f,tmp;

    n = GD_NP;     /* number of points */

    if (alloc_aci(n + 1))          /* inverse permutation */
        return(-1);       

    fprintf(PMFd,"Original distance matrix\n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            UDS_D[i * n + j + 1] = gdd_adj(i,j,PMGN);
            fprintf(PMFd,PMFmtS,UDS_D[i * n + j + 1]);           
        }
        fprintf(PMFd,"\n");           
    }
    UDSDPP = 0.0;
    for (i = 0; i < n; ++i) {
        d = 0.0;
        for (j = 0; j < n; ++j)
            d += gdd_adj(i,j,PMGN);
        if (UDSDPP < d)
            UDSDPP = d;
    }
    printf1("UDSDPP = %g\n",UDSDPP);
    uds_0(n,1);    /* find first feasible permutation */
    uds_r(n,0); 

    for (i = 1; i <= n; ++i) {
        d = uds_get_dip(n,i,UDS_E0) - uds_get_dpi(n,i,UDS_E0);
        d /= (double)n;
        UDS_X[i] = d;                  
        AcI[UDS_E0[i]] = i;     /* inverse permutation */
    }
    /***********************
    printf("UDS_X: ");
    for (i=1; i <= n; ++i)
        printf("%6.2lf ",UDS_X[i]);
    newline();
    printf("UDS_E0: ");
    for (i=1; i <= n; ++i)
        printf("%d ",UDS_E0[i]);
    newline();
    printf("ACI   : ");
    for (i=1; i <= n; ++i)
        printf("%d ",AcI[i]);
    newline();
    ***************/

    printf1("\nBest solution: %lg\n",UDS_F0);

    printf1("\nIndex  X-value     Object  Permutation\n");
    for (i = 1; i <= n; ++i) {
        printf1("%4d  ",i);
        printf1(PMFmtS,UDS_X[i]);
        printf1("%7d   %8d\n",UDS_E0[i],AcI[i]);
    }
    fprintf(PMFd,"\nFitted distance matrix\n");
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) {
            d = fabs(UDS_X[AcI[i]] - UDS_X[AcI[j]]);
            fprintf(PMFd,PMFmtS,d);
            if (j != i) {
                tmp = UDS_D[(i - 1) * n + j] - d;
                f += tmp * tmp;
            }
        }
        fprintf(PMFd,"\n");
    }
    if (f > 0.0)
        f = sqrt(f);

    printf1("\nDistance between original and fitted matrix: %lg\n",f);

    f = 0.0;
    for (i = 2; i <= n; ++i) {
        for (j = 1; j < i; ++j) {
            d = fabs(UDS_X[AcI[i]] - UDS_X[AcI[j]]);
            f += fabs(UDS_D[(i - 1) * n + j] - d);
        }
    }
    f /= (double)(n * (n - 1) / 2);
    printf1("Mean difference between original and fitted distances: %lg\n",f);
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  uds_0(n,k)    find first feasible permutation.                          */

int uds_0(int n,int k)
{
    register int r,i,j,jj;
    int imax,first;
    double b,bmax;   

    for (i = 1; i <= n; ++i) {
        UDS_E0[i] = UDS_PI[i] = i;
        UDS_IF[i] = 0;
    }
    for (r = 1; r <= n; ++r) {
        for (i = 1; i <= n; ++i)  
            UDS_PI[i] = UDS_E0[i];
        first = 1;
        for (j = 1; j <= n; ++j) {
            if (UDS_IF[j])
                continue;
            jj = 0;
            for (i = r; i <= n; ++i) {
                if (UDS_PI[i] == j) {
                    jj = i;
                    break;
                }
            }   
            UDS_PI[jj] = UDS_PI[r];
            UDS_PI[r] = j;

            b = 0.0;
            for (i = 1; i <= r; ++i) 
                b += (uds_get_dip(n,i,UDS_PI) - uds_get_dpi(n,i,UDS_PI));

            if (first || b > bmax) {
                bmax = b;
                imax = j;
                for (i = r; i <= n; ++i)
                    UDS_E0[i] = UDS_PI[i];
                first = 0;
            }
        }
        UDS_IF[imax] = 1;
    }
    UDS_F0 = uds_get_f(n,UDS_E0);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  uds_get_dip(n,i,p)    return d(i,.) for permutation p.                  */

double uds_get_dip(int n,int i,int *p)
{
    register int j,k;
    double d;      

    d = 0.0;
    k = (p[i] - 1) * n;
    for (j = i + 1; j <= n; ++j)     
        d += UDS_D[k + p[j]];
    return(d);
}

/* ------------------------------------------------------------------------ */
/*  uds_get_dpi(n,i,p)    return d(.,i) for permutation p.                  */

double uds_get_dpi(int n,int i,int *p)
{
    register int j,k;
    double d;      

    d = 0.0;
    k = (p[i] - 1) * n;
    for (j = 1; j < i; ++j)  
        d += UDS_D[k + p[j]];
    return(d);
}

/* ------------------------------------------------------------------------ */
/*  uds_get_f(n,p)    return function value for permutation p.              */

double uds_get_f(int n,int *p)
{
    register int i;   
    double d,tmp;  

    d = 0.0;
    for (i = 1; i <= n; ++i) {
        tmp =  uds_get_dip(n,i,p) - uds_get_dpi(n,i,p);
        d += tmp * tmp;
    }
    return(d);
}

/* ------------------------------------------------------------------------ */
/*  uds_r(n,k)      recursively check all permutations.                     */

int uds_r(int n,int k)
{
    register int i,j,l,jj;
    double b,ck,cl,ckk0;

    printf1("Enter uds_r n=%d k=%d E0: ",n,k);
    for (i = 1; i <= n; ++i)
        printf1("%d ",UDS_E0[i]);
    printf1("   F0 = %g\n",UDS_F0);

    for (i = 1; i <= n; ++i) {
        UDS_PI[i] = i;
        UDS_IF[i] = 0;
    }
    for (j = 1; j <= n; ++j) {

        if (UDS_PI[1] != j) {
            jj = 0;
            for (i = 2; i <= n; ++i) {
                if (UDS_PI[i] == j) {
                    jj = i;
                    break;
                }
            }   
            i = UDS_PI[1];
            UDS_IF[i] = 0;
            UDS_PI[jj] = i;
            UDS_PI[1] = j;
        }
        UDS_IF[j] = 1;

        ckk0 = uds_get_dip(n,1,UDS_PI) - uds_get_dpi(n,1,UDS_PI);
        ck = ckk0 * ckk0;
        b = ck + (double)(n - 1) * UDSDPP * UDSDPP;

        if (b >= UDS_F0) {
            cl = 0.0;
            uds_r0(n,1,ck,ckk0,cl);        /* branch to left */
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  uds_r0(n,k,ck,ckk0,cl)                                                  */
/*                                                                          */
/*  branch to node R(r(1),...,r(k) | r(n-k+1),...,r(n))                     */

int uds_r0(int n,int k,double ck,double ckk0,double cl)
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

        if (UDS_IF[i])
            continue;

        if (UDS_PI[l] != i) {
            jj = 0;
            for (j = k + 1; j < l; ++j) {
                if (UDS_PI[j] == i) {
                    jj = j;
                    break;
                }
            }   
            ll = UDS_PI[l];
            UDS_IF[ll] = 0;
            UDS_PI[jj] = ll;
            UDS_PI[l] = i;
        }
        UDS_IF[i] = 1;
        cll0 = uds_get_dip(n,l,UDS_PI) - uds_get_dpi(n,l,UDS_PI);
        cll  = cll0 * cll0;
                                      
        if (cll0 - EPSI1 > ckk0) {          /* not feasible */ 
            UDS_IF[i] = 0;
            continue;
        }
        
        if (l < n) {                /* second check */
         
            cll1 = uds_get_dip(n,l+1,UDS_PI) - uds_get_dpi(n,l+1,UDS_PI);
            if (cll0 < cll1 + 2.0 * UDS_D[(UDS_PI[l] - 1) * n + UDS_PI[l+1]]) {
                UDS_IF[i] = 0;
                continue;
            }
        }
        if (lk1 > 1) {    
            b = ck + cl + cll + (double)lk1 * dmax(ckk,cll);
            if (b + EPSI1 < UDS_F0) {
                UDS_IF[i] = 0;
                continue;
            }
            uds_r1(n,k,l,ck,cl + cll,cll0);     /* branch to right */
        }
        else {
            tmp = uds_get_dip(n,k+1,UDS_PI) - uds_get_dpi(n,k+1,UDS_PI);
            b = ck + cl + cll + tmp * tmp;
            if (b - EPSI1 > UDS_F0) {
                UDS_F0 = b;
                for (j = 1; j <=n; ++j)
                    UDS_E0[j] = UDS_PI[j];
            }
        }
        UDS_IF[i] = 0;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  uds_r1(n,k,l,ck,cl,cll0)                                                */
/*                                                                          */
/*  branch to node R(r(1),...,r(k+1) | r(l),...,r(n))                       */

int uds_r1(int n,int k,int l,double ck,double cl,double cll0)
{
    register int i,j,jj,kk,k1;
    int lk2;
    double b,ckk,ckk0,ckk1,cll,tmp;
  
    k1 = k + 1;
    lk2 = l - k1 - 1;
    cll = cll0 * cll0;

    for (i = 1; i <= n; ++i) {
        if (UDS_IF[i])
            continue;

        if (UDS_PI[k1] != i) {
            jj = 0;
            for (j = k1 + 1; j < l; ++j) {
                if (UDS_PI[j] == i) {
                    jj = j;
                    break;
                }
            }   
            kk = UDS_PI[k1];
            UDS_IF[kk] = 0;
            UDS_PI[jj] = kk;
            UDS_PI[k1] = i;
        }
        UDS_IF[i] = 1;
 
        ckk0 = uds_get_dip(n,k1,UDS_PI) - uds_get_dpi(n,k1,UDS_PI);
       
        if (cll0 - EPSI1 > ckk0) {                  /* not feasible */
            UDS_IF[i] = 0;
            continue;
        }

                /* second check */
          
        ckk1 = uds_get_dip(n,k,UDS_PI) - uds_get_dpi(n,k,UDS_PI);
        if (ckk1 - 2.0 * UDS_D[(UDS_PI[k] - 1) * n + UDS_PI[k1]] + EPSI1 < ckk0) {
            UDS_IF[i] = 0;
            continue;
        }
       
        ckk = ckk0 * ckk0;  

        if (lk2 > 1) {
            b = ck + ckk + cl + (double)lk2 * dmax(ckk,cll);
            if (b + EPSI1 < UDS_F0) {
                UDS_IF[i] = 0;
                continue;
            }
            uds_r0(n,k1,ck + ckk,ckk0,cl);        /* branch to left */
        }
        else {
            tmp = uds_get_dip(n,k1+1,UDS_PI) - uds_get_dpi(n,k1+1,UDS_PI);
            b = ck + ckk + cl + tmp * tmp;
            if (b - EPSI1 > UDS_F0) {
                UDS_F0 = b;
                for (j = 1; j <=n; ++j)
                    UDS_E0[j] = UDS_PI[j];
            }
        }
        UDS_IF[i] = 0;
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

int uds_qa1(void)   
{
    register int i,j,k;
    int n,nn,err;
    double tmp;

    n = GD_NP;
    nn = n * n;

    if (alloc_acu(nn + 1))          /* C matrix */
        return(-1);       
    if (alloc_acv(nn + 1))          /* F matrix */
        return(-1);   
    if (alloc_acz(nn + 1))          /* D matrix */
        return(-1);   
    if (alloc_acn(3 * n + 1))       /* loc3n */
        return(-1);   
    if (alloc_actmp(4 * n + 1))       /* wtmp */
        return(-1);   
    
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            tmp = gdd_adj(i,j,PMGN);
            if (tmp < 0.0) {
                tmp = 0.0;
                if (i != j) {
                    printf1("Error: undefined distance (%d,%d).\n",i+1,j+1);
                    return(-2);
                }
            }
            AcZ[i * n + j + 1] = tmp;
        }
    }

    fprintf(PMFd,"Original distance matrix\n");
    for (i = 0; i < n; ++i) {
        for (j = 1; j <= n; ++j)
            fprintf(PMFd,PMFmtS,AcZ[i * n + j]);
        fprintf(PMFd,"\n");
    }

    for (i = 0; i < n; ++i) {
        k = i;
        for (j = 0; j <= i; ++j) {
            AcV[i * n + j + 1] = (double)(-k);
            k--;
        }
        k = 1;
        for (j = i + 1; j < n; ++j) {
            AcV[i * n + j + 1] = (double)(-k);
            k++;
        }
    }

    /******* 
    printf1("D\n");
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) 
            printf1("%6.3lf ",AcZ[(i - 1) * n + j]);
        newline();
    }

    printf1("F\n");
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) 
            printf1("%6.3lf ",AcV[(i - 1) * n + j]);
        newline();
    }
    ****/

    err = qap_w(n,AcU,AcV,AcZ,AcN,AcTmp);

    if (err) {
        printf1("Error in function qap_w.\n");
        return(-3);
    }         

    /* make inverse permutation */

    for (i = 1; i <= n; ++i)
        AcN[AcN[i] + n] = i;

    printf1("Best function value: %lg\n",-AcTmp[1]);
    printf1("Optimal ordering\n");
    for (i = 1; i <= n; ++i)
        printf1("%4d %4d\n",i,AcN[i]);
    newline();
    printf1("Permutation\n");
    for (i = 1; i <= n; ++i)
        printf1("%4d %4d\n",i,AcN[i + n]);
    newline();

    fprintf(PMFd,"\nFitted distance matrix\n");
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) {
            AcU[(i - 1) * n + j] = tmp = (double)iabs(AcN[i + n] - AcN[j + n]);
            fprintf(PMFd,PMFmtS,tmp);
        }
        fprintf(PMFd,"\n");
    }

    /* check numer of inconstencies */
          
    uds_check(n,AcZ,AcU);
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

int uds_qa2(void)   
{
    register int i,j,k,l;
    int err,n,nn,n1,iter,ito,seed,bestv,look4,d,d1,x,x1;
    int *a,*b,*srtf,*srtif,*srtd,*srtid,*srtc,*srtic,*indexf,*indexd,*cost,*fdind;
    double tmp;

    n = GD_NP;
    nn = n * n;

    if (alloc_acn(nn + 1))          /* f matrix */
        return(-1);       
    if (alloc_acm(nn + 1))          /* d matrix */
        return(-1);   
    
    n1 = k = l = 0;
    for (i = 0; i < n; ++i) {
        l = i;
        for (j = 0; j < n; ++j) {
            tmp = gdd_adj(i,j,PMGN);
            if (tmp < 0.0) {
                tmp = 0.0;
                if (i != j) {
                    printf1("Error: undefined distance (%d,%d).\n",i+1,j+1);
                    return(-2);
                }
            }
            AcN[++k] = (int)tmp;
            AcM[k] = -l;
            if (j >= i)
                l++;
            else
                l--;
        }
    }

    fprintf(PMFd,"Original distance matrix\n");
    for (i = 0; i < n; ++i) {
        for (j = 1; j <= n; ++j)
            fprintf(PMFd,PMNFmtS,AcN[i * n + j]);
        fprintf(PMFd,"\n");
    }


    /******* 
    printf1("F\n");
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) 
            printf1("%2d ",AcN[(i - 1) * n + j]);
        newline();
    }

    printf1("D\n");
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) 
            printf1("%2d ",AcM[(i - 1) * n + j]);
        newline();
    }
    ****/
              
    if (alloc_acs(n + 1))           /* optimal permutation vector */
        return(-1);      

    if (alloc_acr(10 * nn + 2 * n + 2))     /* working storage */
        return(-1);            

    a = AcR;
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
   
    gqapd(n,MxIter,PMIDFA,PMIDFB,look4,&seed,AcN,AcM,a,b,srtf,srtif,srtd,srtid,
          srtc,srtic,indexd,indexf,cost,fdind,AcS,&bestv,&iter,&ito);

    printf1("Performed %d iterations.\n",iter);
    printf1("Best cost value: %d  (found in %d iterations).\n",-bestv,ito);
    n1 = 0;
    for (i = 1; i <= n; ++i) {
        k = AcS[i];
        for (j = 1; j <= n; ++j) {
            l = AcS[j];
            n1 += AcM[(i - 1) * n + j] * AcN[(k - 1) * n + l];
        }
    }
    printf1("Checked cost value: %d\n",-n1);
        
    /* make inverse permutation */

    for (i = 1; i <= n; ++i)
        AcR[AcS[i]] = i;

    printf1("Optimal ordering\n");
    for (i = 1; i <= n; ++i)
        printf1("%4d %4d\n",i,AcS[i]);
    newline();
    printf1("Permutation\n");
    for (i = 1; i <= n; ++i)
        printf1("%4d %4d\n",i,AcR[i]);
    newline();

    fprintf(PMFd,"\nFitted distance matrix\n");
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) {
            AcM[(i - 1) * n + j] = k = iabs(AcR[i] - AcR[j]);
            fprintf(PMFd,PMNFmtS,k);
        }
        fprintf(PMFd,"\n");
    }

    /* check numer of inconstencies */
          
    uds_ncheck(n,AcN,AcM);
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

int sga(void)
{
    register int i,j,k;
    int err,n,m,s,nn,nf,len,err0,err1,err2;
    double d;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Scalogram analysis. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 3,4,1))       /* get parameters */
        goto SGAFin;

    n = NOC;
    m = PMNV;                       /* number of variables */

    if (alloc_acn(m * n + 1))
        goto SGAFin;
    if (alloc_acs(n + 1))
        goto SGAFin;
    if (alloc_aci(n + 1))
        goto SGAFin;
    if (alloc_acr(m + 1))
        goto SGAFin;
    if (alloc_ack(m + 1))
        goto SGAFin;
    if (alloc_acm(n + 1))
        goto SGAFin;

    for (i = 0; i < NOC; ++i) {                 /* number of cases */
        for (j = 0; j < m; ++j) {
            d = get_data(PMVIdx[j],i);
            if (d > 0.0) {
                AcS[i] += 1;
                AcR[j] += 1;
            }
        }
    }
    if (sortdpi(n,AcS,AcI))     /* sort scores */
        goto SGAFin;

    if (sortdpi(m,AcR,AcK))     /* sort difficulties */
        goto SGAFin;

    /* put reduced and sorted data matrix into AcN */

    k = AcI[0];
    for (j = 0; j < m; ++j) {
        d = get_data(PMVIdx[AcK[j]],k);
        if (d > 0.0)
            AcN[j] = 1;
    }
    AcM[0] = 1;
    nn = 1;
    for (i = 1; i < n; ++i) {
        k = AcI[i];
        nf = 1;
        for (j = 0; j < m; ++j) {
            if (get_data(PMVIdx[AcK[j]],k) > 0.0)
                s = 1;
            else
                s = 0;
            AcN[nn * m + j] = s;
            if (s != AcN[(nn - 1) * m + j])
                nf = 0;
        }
        if (nf == 1)
            AcM[nn - 1] += 1;
        else {
            AcM[nn] = 1;
            nn++;
        }
    }
    if (alloc_acr(m))
        goto SGAFin;
    if (alloc_acs(nn))
        goto SGAFin;

    for (i = 0; i < nn; ++i) {
        for (j = 0; j < m; ++j) {
            if (AcN[i * m + j] == 1) {
                AcS[i] += 1;
                AcR[j] += AcM[i];
            }
        }
    }
    if (alloc_aci(nn))      /* error1 (Goodenough) */
        goto SGAFin;
    if (alloc_acj(nn))      /* error2 (Guttman) */
        goto SGAFin;

    len = 4;
    for (j = 0; j < m; ++j)
        len = imax(len,strlen(VName[PMVIdx[j]]));

    printf1("\n Idx ");
    for (j = 0; j < m; ++j) {
        k = AcK[j];
        prnchar(' ',len - strlen(VName[PMVIdx[k]]),0);
        printf1("%s ",VName[PMVIdx[k]]);
    }
    printf1("   Freq   Value   Err0   Err1   Err2\n");
    prnchar('-',41 + m * (len + 1),1);

    err0 = err1 = err2 = 0;
    for (i = 0; i < nn; ++i) {
        printf1("%4d ",i + 1);

        s = m - AcS[i];
        for (j = 0; j < m; ++j) {
            k = AcN[i * m + j];                         
            prnchar(' ',len - 4,0);
            printf1("%4d ",k);                         
            
            if (k == 1) {
                if (j < s)
                    AcI[i] += 1;
            }
            else {
                if (j >= s)
                    AcI[i] += 1;
            }
        }
        s = imin(AcI[i],1);
        err0 += s * AcM[i];

        err1 += AcI[i] * AcM[i];

        AcJ[i] = sga1(m,AcN + i * m);
        err2 += AcJ[i] * AcM[i];

        printf1("%7d %7d %6d %6d %6d\n",AcM[i],AcS[i],s,AcI[i],AcJ[i]);                  
    }
    prnchar('-',41 + m * (len + 1),1);
    prnchar(' ',5,0);

    for (j = 0; j < m; ++j) {
        prnchar(' ',len - 4,0);
        printf1("%4d ",AcR[j]);
    }
    printf1("%7d     %10d %6d %6d\n\n",n,err0,err1,err2);

    printf1("Reproducibility\n");
    printf1("Err0: %6.4f\n",1.0 - (double)err0 / (double)n);
    printf1("Err1: %6.4f\n",1.0 - (double)err1 / (double)(n * m));
    printf1("Err2: %6.4f\n\n",1.0 - (double)err2 / (double)(n * m));


    hpo(nn,m,AcN);     

    err = 0;

SGAFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sga1(m,r)   given profile (r1,...,rm) return minimal distance to        */
/*              perfect profile (Guttman method).                           */

int sga1(int m,int *r)
{
    register int i,j,mm,e;
    int err;

    err = INTMAX;
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
        err = imin(err,e);
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  hpo(n,m,x)                                                              */

int hpo(int n,int m,int *x)      
{
    register int i,j,k,im;
    int err,d;

    err = -1;
    if (alloc_acd(n * n))   
        goto HPO1Fin;

    for (i = 0; i < n; ++i) {
        im = i * n;
        for (j = 0; j < n; ++j) {
            if (j == i)
                continue;

            d = hpo1(m,x + i * m,x + j * m);
            if (d != 1)
                continue;
        
            for (k = 0; k < n; ++k) {
                if (k == i || k == j) {
                    continue;
                }
                if (hpo1(m,x + i * m,x + k * m) == 1 &&
                    hpo1(m,x + k * m,x + j * m) == 1) {
                    d = 0;
                    break;
                }
       
            }
            if (d == 1)
                AcD[im + j] = (char)d;
        }
    }

    for (i = 0; i < n; ++i) {
        im = i * n;
        for (j = 0; j < n; ++j)  
            printf1("%d ",(int)AcD[im + j]);
        newline();
    }

    printf1("digraph g {\n");
    for (i = 0; i < n; ++i) {
        printf1("n%d [label=%d];\n",i+1,i+1);
    }

    for (i = 0; i < n; ++i) {
        im = i * n;
        for (j = 0; j < n; ++j) {
            if (AcD[im + j])
                printf1("n%d -> n%d;\n",i+1,j+1);
        }
    }



    err = 0;
HPO1Fin:
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  hpo1(m,x,y)     return 1 if x <= y, otherwise 0.                        */

int hpo1(int m,int *x,int *y)
{
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

int unf(void)
{
    register int i,j,k;
    int err,n,m,first,wf;
    double wt,nn;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Unfolding. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 3,4,1))       /* get parameters */
        goto UNFFin;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    n = NOC;
    m = PMNV;                       /* number of variables */

    if (m < 2) {
        printf1("Error: need at least two variables.\n");
        goto UNFFin;
    }
    if (alloc_acx(n * m + 1))           /* input data */
        goto UNFFin;

    wf = 0;
    if (PMWVar >= 0) {
        printf1("Using case weights defined by: %s\n",VName[PMWVar]);
        if (alloc_acw(n + 1))   
            goto UNFFin;
        wf = 1;
    }
    if (alloc_acm(m + 1))               /* best permutation */
        goto UNFFin;

    for (i = 0; i < n; ++i) {           /* get data into AcX */
        for (j = 0; j < m; ++j)  
            AcX[i * m + j] = get_data(PMVIdx[j],i);

        if (PMWVar >= 0) {
            wt = get_data(PMWVar,i);
            if (wt < 0.0) {
                printf1("Error: found negative weight in case %d\n",i + 1);
                goto UNFFin;
            }
            AcW[i] = wt;
        }
    }
                     
/***
    for (i = 0; i < n; ++i) {   
        for (j = 0; j < m; ++j)  
            printf1("%g ",AcX[i * m + j]);
        newline();
    }
***/

    /*  do for all permutations of m alternatives */

    if (alloc_ack(m))
        goto UNFFin;
    if (alloc_aci(m))
        goto UNFFin;
    if (alloc_acj(m))
        goto UNFFin;
    if (alloc_acr(m))
        goto UNFFin;
    if (alloc_acs(m * m + 1))
        goto UNFFin;
    if (alloc_acn(m + 1))
        goto UNFFin;
    if (alloc_acns(m + 1))
        goto UNFFin;

    UNFMax = -(DBLMAX - 100);    
    nn = (double)(m * (m - 1) / 2);

    if (PMNTP > 0) {
        printf1("\nCompare with rank order: %d",(int)PMTP[0]);
        for (j = 1; j < PMNTP; ++j)
            printf1(", %g",PMTP[j]);
        newline();
        if (PMNTP != m) {
            printf1("Error in number of alternatives.\n");
            goto UNFFin;
        }
        for (j = 0; j < m; ++j) {
            AcK[j] = k = (int)PMTP[j] - 1;
            if (k < 0 || k >= m) {
                printf1("Error: invalid permutation.\n");
                goto UNFFin;
            }
        }
        unf1(n,m,AcK,AcR,AcX,AcS,AcN,0,AcNS,wf,AcW);  

        printf1("Value (sum of taus): %g\n\n",UNFMax / nn);
        /***
        printf1("Best permutation: ");
        for (j = 0; j < m; ++j)
            printf1("%d ",AcM[j]);
        newline();
        **/
        err = 0;
        goto UNFFin;
    }











    first = 1;
    while (perm(m,AcK,AcI,AcJ,first)) {

        if (AcK[0] > AcK[m - 1])        /* half of permutations suffices */
            continue;

        unf1(n,m,AcK,AcR,AcX,AcS,AcN,0,AcNS,wf,AcW);  /* check sub-permutations */

printf1("UNFMax=%g\n",UNFMax);

        first = 0;
    }
    printf1("\nBest value (sum of taus): %g\n",UNFMax / nn);
    printf1("Best permutation: ");
    for (j = 0; j < m; ++j)
        printf1("%d ",AcM[j]);
    newline();

    if (PMF1Def) {
        for (j = 0; j < m; ++j)
            AcK[j] = AcM[j] - 1;

        unf1(n,m,AcK,AcR,AcX,AcS,AcN,1,AcNS,wf,AcW);

        printf1("\n%d records written to: %s\n",n,PMF1dName);
    }
    err = 0;

UNFFin:       
    p_clean();
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

void unf1(int n,int m,int *p,int *r,double *x,int *s,int *rr,int prn,
    short *pp,int wf,double *w)
{
    register int i,j,k,l,t;  
    int ii,c,d;
    double xij,cmax,nn;
        
    if (prn)
        nn = (double)(m * (m - 1) / 2);

    cmax = 0.0;                   /* max over all cases */
    for (ii = 0; ii < n; ++ii) {

        c = -(INTMAX - 100);    
                         
        l = ii * m;
        for (j = 0; j < m; ++j) {   
            r[j]= p[j] + 1;
            t = j * m;
            for (k = 0; k < j; ++k)  
                s[t + k] = (int)dsign(x[l + j] - x[l + k]);
        }
printf("vor unf2\n");
printf("r: ");
for (j = 0; j < m; ++j)
    printf("%d ",r[j]);
newline();
printf("x: ");
for (j = 0; j < m; ++j)
    printf("%g ",x[l + j]);
newline();

        d = unf2(m,r,s,rr);
printf("   corr=%d\n",d);
exit(0);

        if (d > c) {
            c = d;
            if (prn) {
                for (j = 0; j < m; ++j)
                    pp[r[j] - 1] = m - j;              
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
            d = unf2(m,r,s,rr);
            if (d > c) {
                c = d;
                if (prn) {
                    for (j = 0; j < m; ++j)
                        pp[r[j] - 1] = m - j;              
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
            d = unf2(m,r,s,rr);
            if (d > c) {
                c = d;
                if (prn) {
                    for (j = 0; j < m; ++j)
                        pp[r[j] - 1] = m - j;              
                }
            }
        }
        if (wf)
            cmax += (double)c * w[ii];
        else
            cmax += (double)c;

        if (prn) {
            fprintf(PMF1d,"%5d  ",ii + 1);
            for (j = 0; j < m; ++j)
                fprintf(PMF1d,"%2d ",pp[j]);
            fprintf(PMF1d," ");
            t = 0;
            for (j = 0; j < m; ++j) {
                xij = x[ii * m + j];
                fprintf(PMF1d,"%2lg ",xij);
                if (xij != (double)pp[j])
                    t = 1;
            }
            fprintf(PMF1d,"  %d ",t);
            fprintf(PMF1d,PMFmtS,(double)c / nn);
            if (wf)
                fprintf(PMF1d,PMFmtS,w[ii]);
            fprintf(PMF1d,"\n");
        }
    }
    if (prn)
        return;

    if (cmax > UNFMax) {        /* save new max and permutation */
        UNFMax = cmax;
        for (j = 0; j < m; ++j) 
            AcM[j] = p[j] + 1;
    }
}

/* ------------------------------------------------------------------------ */
/*  unf2(m,r,s)   return rank correlation (tau) between rr and x.           */
/*                r is set up by unf1. to get rr create rank orders:        */
/*                rr[r[j]-1] = m-j.                                         */
/*                                                                          */
/*                s[] is a m x m matrix containing sign(xi - xj).           */

int unf2(int m,int *r,int *s,int *rr)
{
    register int i,j;
    int t;
       
    for (j = 0; j < m; ++j)
        rr[r[j] - 1] = m - j;              

    t = 0;
    for (i = 1; i < m; ++i) {
        for (j = 0; j < i; ++j)  
            t += isign(rr[i] - rr[j]) * s[i * m + j];
    }

printf("rr: ");
for (j = 0; j < m; ++j)
printf("%d ",rr[j]);
newline();

         
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
       
int mdsc(void)
{
    register int i,j,k;
    int err,n,r,nn,rflag,n1rec;
    double d,di,t,xi1,xi2,xj1,xj2,evmin,c;

    n1rec = rflag = 0;
    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("MDS with principal coordinates. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto MDSCFin;

    if (gdd_check(PMGN,1))
        goto MDSCFin;
    if (gdd_tcheck(0,1,2))          /* need undirected valued graph */
        goto MDSCFin;

    if (GD_NP < 3) {
        printf1("Error: need at least three nodes.\n");
        goto MDSCFin;
    }
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    n = GD_NP;                      /* number of points */

    /* create A matrix in AcU */

    if (alloc_acu(n * n + 2))                                    
        goto MDSCFin; 
    if (alloc_acx(n + 1))                                    
        goto MDSCFin; 


MDSCRep:

    d = 0.0;
    for (i = 0; i < n; ++i) {
        di = 0.0;
        for (j = 0; j < n; ++j) {
            if (j != i) {
                t = gdd_adj(i,j,PMGN);
                if (rflag)
                    t += c;

                di += t * t;
                d  += t * t;
            }
        }
        AcX[i] = di;
    }   
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (i != j) {
                t = gdd_adj(i,j,PMGN);
                if (rflag)
                    t += c;
            }
            else 
                t = 0.0;

            AcU[i * n + j + 1] = -0.5 * (t * t - AcX[i] / (double)n
                                               - AcX[j] / (double)n
                                               + d / (double)(n * n));
        }
    }

    /* if requested write A-matrix */

    if (PMF1Def) {
        if (rflag == 0)
            fprintf(PMF1d,"A-matrix\n");
        else
            fprintf(PMF1d,"A-matrix (of modified distance matrix)\n");
        n1rec++;
        for (i = 0; i < n; ++i) {
            for (j = 1; j <= n; ++j)  
                fprintf(PMF1d,PMFmtS,AcU[i * n + j]);
            fprintf(PMF1d,"\n");
            n1rec++;
        }
    }

    r = evecf(n,AcX,AcU);
    if (r) {
        if (r == -1) {
            printf1("No success in eigenvalue calculation.\n");
            printf1("Exceeded maximum number of iterations.\n");
        }
        else
            printf1("Insufficient memory.\n");
        goto MDSCFin;
    }
            
    if (PMF1Def) {
        fprintf(PMF1d,"Eigenvalues\n");
        n1rec++;
        for (j = 1; j <= n; ++j)   
            fprintf(PMF1d,PMFmtS,AcX[j]);
        fprintf(PMF1d,"\n");
        n1rec++;
    }
    printf1("\nEigenvalue         per cent\n");
    r = 0;
    t = 0.0;
    for (j = 1; j <= n; ++j) {
        if (fabs(AcX[j]) <= EPSI1)
            AcX[j] = 0.0;
        if (AcX[j] < 0.0)  
            r++;
        t += AcX[j];
    }
    for (j = 1; j <= n; ++j) {
        prnchar(' ',16 - PMFmt1,0);
        printf1(PMFmtS,AcX[j]);
        if (t > 0.0 && r == 0)
            printf1("  %8.2f",100.0 * AcX[j] / t);
        newline();
    }
    newline();

    if (rflag == 0 && PMOPT == 2) {       /* find additive constant */
        evmin = AcX[n];
        c = 4.0 * evmin * evmin - 2.0 * evmin;
        if (c > 0.0)
            c = sqrt(c);
        else
            c = 0.0;
        c -= 2.0 * evmin;
        c += EPSI1;
        printf1("Additive constant: %lg\n",c);
        if (c <= 0.0)
            goto MDSCFin;
        rflag = 1;
        goto MDSCRep;
    }
                 
    if (PMF1Def) {
        fprintf(PMF1d,"Eigenvectors\n");
        n1rec++;
        for (i = 0; i < n; ++i) {
            for (j = 1; j <= n; ++j)  
                fprintf(PMF1d,PMFmtS,AcU[i * n + j]);
            fprintf(PMF1d,"\n");
            n1rec++;
        }
    }

    r = 0;
    for (j = 1; j <= n; ++j) {
        if (AcX[j] > 0.0)
            r++;
        else
            break;
    }
    if (r > 0) {
        for (i = 0; i < n; ++i) {
            for (j = 1; j <= r; ++j)  
                fprintf(PMFd,PMFmtS,sqrt(AcX[j]) * AcU[i * n + j]);
            fprintf(PMFd,"\n");
        }
        printf1("%d records written to: %s\n",n,PMFdName);
    }
    if (n1rec > 0)
        printf1("%d records written to: %s\n",n1rec,PMF1dName);


    /*******************************/
    if (PMTabFDef && r >= 2) {
        nn = n * (n - 1) / 2;
        if (alloc_acv(nn + 1))                                    
            goto MDSCFin; 
        if (alloc_acw(nn + 1))                                    
            goto MDSCFin; 
        if (alloc_ack(nn + 1))                                    
            goto MDSCFin; 

        k = 0;
        for (i = 0; i < n; ++i) {
            xi1 = sqrt(AcX[1]) * AcU[i * n + 1];
            xi2 = sqrt(AcX[2]) * AcU[i * n + 2];

            for (j = 0; j < i; ++j) {
                AcV[k] = gdd_adj(i,j,PMGN);
                xj1 = sqrt(AcX[1]) * AcU[j * n + 1];
                xj2 = sqrt(AcX[2]) * AcU[j * n + 2];

                AcW[k] = sqrt(pow(xi1 - xj1,2) + pow(xi2 - xj2,2));
                k++;
            }
        }
        if (sortdp2(k,AcV,AcW,AcK))
            goto MDSCFin;

        for (i = 0; i < k; ++i) {
            j = AcK[i];
            fprintf(PMTabFd,"%3d ",i);
            fprintf(PMTabFd,PMFmtS,AcV[j]);
            fprintf(PMTabFd,PMFmtS,AcW[j]);
            fprintf(PMTabFd,"\n");
        }
        printf1("%d records written to: %s\n",k,PMTabFName);
    }
    /*****************/

    err = 0;

MDSCFin:       
    p_clean();
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
       
int mdsm(void)
{
    register int i,j,k;
    int err,nn,ii;
    double d,a,b,aa,bb,dm,fm,ds;

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("Metric MDS (experimental). Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto MDSMFin;

    if (gdd_check(PMGN,1))
        goto MDSMFin;
    if (gdd_tcheck(0,1,2))          /* need undirected valued graph */
        goto MDSMFin;

    if (GD_NP < 3) {
        printf1("Error: need at least three nodes.\n");
        goto MDSMFin;
    }
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    MDSRN = GD_NP;                  /* number of points */

    /* setup MDSDIS vector with distances */

    nn = MDSRN * (MDSRN - 1) / 2;
    if (alloc_acv(nn + 1))                                        
        goto MDSMFin; 
    if (alloc_acw(nn + 1))                                        
        goto MDSMFin; 
    MDSDIS = AcV;
    MDSDIS1 = AcW;

    ds = dm = 0.0;
    k = 0;
    for (i = 2; i <= MDSRN; ++i) {
        for (j = 1; j < i; ++j) {
            d = gdd_adj(i-1,j-1,PMGN);
            MDSDIS[k++] = d;
            dm = dmax(dm,d);
            ds += d * d;
        }
    }

    /* prepare for min algorithm, always MINA = 7 */

    if (MINA == 5 || MINA == 6 || MINA == 8)
        MINA = 7;

    NParm = 2 * (MDSRN - 1);
    PMDOPT = 1;                     /* analytical gradient */
    CCTyp = 0;                      /* no covariance matrix */

    set_mlopt();                    /* adjust options */         

    if (ml_init(1,1,0))             /* init function minimization */
        goto MDSMFin;

    if (get_dsv(NParm,Par,0,ParLB,ParUB,1))   /* try to get starting values */
        goto MDSMFin;
       
    if (DSVFlg)  
        PMNS = 1;
    else if (PMNS < 1)
        PMNS = 1;

    if (alloc_act(nn + 1))      /* used for best coordinates */   
        goto MDSMFin; 

    for (ii = 1; ii <= PMNS; ++ii) {    /* repeat */

        if (DSVFlg == 0) {              /* get initial configuration */
            for (i = 1; i <= NParm; ++i) 
                Par[i] = dm * random1();
        }
        if (ffmin(MDSMOD4,0,1,1))       /* minimization */
            goto MDSMFin;               /* insuff memory or fn error */

        prn_mlres(1);                   /* print info about minimization */

        /* print result into output file */

        fprintf(PMFd,"%4d %3d %12.8lf  ",ii,LConv,FMin);
        for (i = 1; i <= NParm; ++i)
            fprintf(PMFd,PMFmtS,Par[i]);
        fprintf(PMFd,"\n");

        if (ii == 1) {
            fm = FMin;
            for (i = 1; i <= NParm; ++i)
                AcT[i] = Par[i];
        }
        else if (FMin < fm) {
            fm = FMin;     
            for (i = 1; i <= NParm; ++i)
                AcT[i] = Par[i];
        }
        if (DSVFlg)
            break;
    }
    printf1("Best function value: %14.6lf (stress: %8.6lf)\n",fm,sqrt(fm/ds));
    printf1("Coordinates\n");
    for (i = 1; i <= MDSRN; ++i) {
        printf1("%3d ",i);
        if (i == 1)
            a = b = 0.0;
        else {
            a = AcT[i - 1];
            b = AcT[MDSRN + i - 2];
        }
        printf1(PMFmtS,a);
        printf1(PMFmtS,b);
        newline();
    }
    newline();

    if (PMF1Def) {      /* write estimated distance matrix */
        for (i = 1; i <= MDSRN; ++i) {
            if (i == 1)
                a = b = 0.0;
            else {
                a = AcT[i - 1];
                b = AcT[MDSRN + i - 2];
            }
            for (j = 1; j <= MDSRN; ++j) {
                if (j == 1)
                    aa = bb = 0.0;
                else {
                    aa = AcT[j - 1];
                    bb = AcT[MDSRN + j - 2];
                }
                d = mds_edis(a,b,aa,bb);
                fprintf(PMF1d,PMFmtS,d);
            }
            fprintf(PMF1d,"\n");
        }
        printf1("%d records written to %s\n",MDSRN,PMF1dName);
    }
    printf1("%d records written to %s\n",PMNS,PMFdName);

    err = 0;

MDSMFin:       
    ml_init(0,0,0);      
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  mds4_fn()   Calculate function value for MDSMOD4.                       */
/*              Use parameter values in TPar[] (i=1,NParm).                 */
/*              Return function value in FTmp, gradient in Grad[].          */
/*                                                                          */
/*              Return 0 if OK, 1 if error in function                      */

int mds4_fn(void)
{
    register int i,j,k,l,il,lj;
    int err;
    double ai,bi,aj,bj,al,bl,d,t;
         
    /* get estimated distances from parameters */

    k = 0;
    for (i = 2; i <= MDSRN; ++i) {
        ai = TPar[i - 1];
        bi = TPar[i + MDSRN - 2];                       
        for (j = 1; j < i; ++j) {
            if (j == 1)  
                aj = bj = 0.0;
            else {
                aj = TPar[j - 1];
                bj = TPar[j + MDSRN - 2];                       
            }
            d = mds_edis(ai,bi,aj,bj);
            MDSDIS1[k] = d;
            d = MDSDIS[k] - d;
            if (LFunc)
                FTmp += d * d;
            k++;
        }
    }
    if (LGrad) {
        for (l = 2; l <= MDSRN; ++l) {
            al = TPar[l - 1];
            bl = TPar[l + MDSRN - 2];                       
            lj = ((l - 1) * (l - 2)) / 2;
            t = MDSDIS[lj];
            if ((d = mds_elen(al,bl)) != 0.0)
                t /= d;
            Grad[l-1] -= 2.0 * al * (t - 1.0);

            for (j = 2; j < l; ++j) {
                aj = TPar[j - 1];
                bj = TPar[j + MDSRN - 2];                       
                lj = ((l - 1) * (l - 2)) / 2 + j - 1;
                t = MDSDIS[lj];
                if ((d = mds_edis(al,bl,aj,bj)) != 0.0)
                    t /= d;
                Grad[l-1] -= 2.0 * (al - aj) * (t - 1.0);
            }
            for (i = l + 1; i <= MDSRN; ++i) {
                ai = TPar[i - 1];
                bi = TPar[i + MDSRN - 2];                       
                il = ((i - 1) * (i - 2)) / 2 + l - 1;
                t = MDSDIS[il];
                if ((d = mds_edis(ai,bi,al,bl)) != 0.0)
                    t /= d;
                Grad[l-1] += 2.0 * (ai - al) * (t - 1.0);
            }
        }
        for (l = 2; l <= MDSRN; ++l) {
            al = TPar[l - 1];
            bl = TPar[l + MDSRN - 2];                       
            lj = ((l - 1) * (l - 2)) / 2;
            t = MDSDIS[lj];
            if ((d = mds_elen(al,bl)) != 0.0)
                t /= d;
            Grad[MDSRN + l - 2] -= 2.0 * bl * (t - 1.0);

            for (j = 2; j < l; ++j) {
                aj = TPar[j - 1];
                bj = TPar[j + MDSRN - 2];                       
                lj = ((l - 1) * (l - 2)) / 2 + j - 1;
                t = MDSDIS[lj];
                if ((d = mds_edis(al,bl,aj,bj)) != 0.0)
                    t /= d;
                Grad[MDSRN + l - 2] -= 2.0 * (bl - bj) * (t - 1.0);
            }
            for (i = l + 1; i <= MDSRN; ++i) {
                ai = TPar[i - 1];
                bi = TPar[i + MDSRN - 2];                       
                il = ((i - 1) * (i - 2)) / 2 + l - 1;
                t = MDSDIS[il];
                if ((d = mds_edis(ai,bi,al,bl)) != 0.0)
                    t /= d;
                Grad[MDSRN + l - 2] += 2.0 * (bi - bl) * (t - 1.0);
            }
        }
    }
    err = checkov();          /* check overflow */
    if (err) {
        printf1("\nError in function evaluation.\n");
        printf1("Numerical overflow.\n");
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
       
int mdsn(void)
{
    register int i,j,k;
    int err,ii,nties,iter;
    double d,a,b,aa,bb,dm,fm,ds,st,stp,s,t,stmin,fval,ngrad,ngradp;
    double mag,slen,af,rf,gf,cosphi,st5,sthis[6];

    SLen = 0.2;

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("Nonmetric MDS (experimental). Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto MDSNFin;

    if (gdd_check(PMGN,1))
        goto MDSNFin;
    if (gdd_tcheck(0,1,2))          /* need undirected valued graph */
        goto MDSNFin;

    if (GD_NP < 3) {
        printf1("Error: need at least three nodes.\n");
        goto MDSNFin;
    }
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    if (MxItFlg == 0)
        MxIter = 100;

    MDSRN = GD_NP;                       /* number of points */
    MDSRNN = MDSRN * (MDSRN - 1) / 2;
    MDSND = PMNDIM;
    if (MDSND < 0)
        MDSND = 2;
    MDSNP = MDSRN * MDSND;

    printf1("Number of points: %d\n",MDSRN);
    printf1("Number of dimensions: %d\n",MDSND);
    if (MDSND < 1 || MDSND > 3) {
        printf("Error in number of dimensions.\n");
        goto MDSNFin;
    }

    /* allocate memory */

    if (alloc_acv(MDSRNN + 1))      /* used for MDSDIS (sorted) */             
        goto MDSNFin; 
    if (alloc_acw(MDSRNN + 1))      /* used for MDSDIS1 */            
        goto MDSNFin; 
    if (alloc_acu(MDSRNN + 1))      /* used for MDSDIS2 */            
        goto MDSNFin; 
    if (alloc_act(MDSRNN + 1))      /* used for MDSDIS3 */            
        goto MDSNFin; 
    if (alloc_acs(MDSRNN + 1))      /* used for MDSPTR (ptr for MDSDIS) */             
        goto MDSNFin; 
    if (alloc_acr(MDSRNN + 1))      /* inverse pointer */             
        goto MDSNFin; 
    if (alloc_acm(MDSRNN + 1))      /* used for tie blocks */             
        goto MDSNFin; 
    if (alloc_aci(MDSRNN + 1))                                            
        goto MDSNFin; 
    if (alloc_acj(MDSRNN + 1))                                            
        goto MDSNFin; 
    if (alloc_acx(MDSNP + 1))       /* used for configuration (= parameters) */
        goto MDSNFin; 
    if (alloc_acz(MDSNP + 1))       /* used for best configuration */
        goto MDSNFin; 
    if (alloc_acy(MDSNP + 1))       /* used for gradient */
        goto MDSNFin; 
    if (alloc_actmp(MDSNP + 1))     /* used for previous gradient */
        goto MDSNFin; 

    MDSDIS = AcV;
    MDSDIS1 = AcW;
    MDSDIS2 = AcU;
    MDSDIS3 = AcT;
    MDSPTR = AcS;
    MDSPTRI = AcR;
    MDSX = AcX;
    MDSG = AcY;
    MDSGP = AcTmp;
    MDSXB = AcZ;        

    /* put lower triangle of distance matrix into MDSDIS */

    ds = dm = 0.0;
    k = 0;
    for (i = 2; i <= MDSRN; ++i) {
        for (j = 1; j < i; ++j) {
            d = gdd_adj(i-1,j-1,PMGN);
            MDSDIS[k++] = d;
            dm = dmax(dm,d);
            ds += d * d;
        }
    }

    printf("MDSDIS unsorted: ");
    for (i = 0; i < MDSRNN; ++i)
        printf("%6.3lf ",MDSDIS[i]);
    newline();

    printf1("Largest distance value: %lg\n",dm);

    /* sort MDSDIS in ascending order. pointer: MDSPTR */

    if (sortdp(MDSRNN,MDSDIS,MDSPTR))
        goto MDSNFin;

    printf("PTR: ");
    for (i = 0; i < MDSRNN; ++i)
        printf("%d (%d) ",MDSPTR[i],i);
    newline();

    /* make inverse pointer in MDSPTRI */

    for (i = 0; i < MDSRNN; ++i)
        MDSPTRI[MDSPTR[i]] = i;

    /* rearrange MDSDIS into sorted order */

    for (i = 0; i < MDSRNN; ++i)
        MDSDIS1[i] = MDSDIS[MDSPTR[i]];

    for (i = 0; i < MDSRNN; ++i)
        MDSDIS[i] = MDSDIS1[i];

    printf("MDSDIS   sorted: ");
    for (i = 0; i < MDSRNN; ++i)
        printf("%6.3lf ",MDSDIS[i]);
    newline();

    printf("MDSDIS unsorted: ");
    for (i = 0; i < MDSRNN; ++i)
        printf("%6.3lf ",MDSDIS[MDSPTRI[i]]);
    newline();




    nties = 0;
    a = -1.0;
    for (i = 0; i < MDSRNN; ++i) {
        b = MDSDIS[i];
        if (b <= a + EPSI1)
            nties++;
        else
            a = b;
    }
    printf1("Number of ties: %d\n",nties);
    if (nties > 0) {
        printf1("Handling of ties: ");
        if (PMOPT == 2)
            printf1("Kruskal's primary approach.\n");
        else if (PMOPT == 3)
            printf1("Kruskal's secondary approach.\n");
        else {
            PMOPT = 1;
            printf1("ignored.\n");
        }
    }
    nties = mreg_ties(MDSRNN,MDSDIS,AcM);
         
    printf("\nNumber of tie blocks nties=%d\n",nties);
    for (i = 0; i < nties; ++i)
        printf("i=%3d AcM[i]=%d\n",i,AcM[i]);
                           

    /*************************************************************************
    if (get_dsv(NParm,Par,0,ParLB,ParUB,1))      try to get starting values   
        goto MDSNFin;
    if (DSVFlg)  
        PMNS = 1;
    else if (PMNS < 1)
        PMNS = 1;
    *****************************/
      

    printf1("\nParameters for iterative minimization.\n");
    printf1("Maximal number of iterations: %d\n",MxIter);
    printf1("Initial step length: %lg\n",SLen);
    printf1("Tolerance for norm of gradient: %lg\n",TOLG);
    printf1("Tolerance for stress: %lg\n\n",TOLS);

    printf1("Number of random repeats: %d\n\n",PMNS);

    for (ii = 1; ii <= PMNS; ++ii) {    /* repeat */

        printf1("\nRepeat %d\n",ii);

        if (DSVFlg == 0) {              /* get initial configuration */
            for (i = 1; i <= MDSNP; ++i) 
                MDSX[i] = dm * random1();
        }
/**
        printf("CONFIG: ");
        for (i = 1; i <= MDSNP; ++i)
            printf("%6.3lf ",MDSX[i]);
        newline();
**/
        slen = SLen;

        for (iter = 1; iter <= MxIter; ++iter) {

            /* normalize the configuration and create distances in MDSDIS1 */

            if (mdsn_norm())
                goto MDSNFin;

/**
            printf("MDSDIS1 nach norm: ");
            for (i = 0; i < MDSRNN; ++i)
                printf("%6.3lf ",MDSDIS1[i]);
            newline();
**/

            /* find monotone projection in MDSDIS2 */

            for (i = 0; i < MDSRNN; ++i)
                MDSDIS2[i] = MDSDIS1[i];
     
            if (PMOPT == 1)
                mreg_proj(MDSRNN,MDSDIS2);
    
            else if (PMOPT == 2)
                mreg_proj1(MDSRNN,nties,AcM,AcI,AcJ,MDSDIS2,MDSDIS3);
      
            else if (PMOPT == 3)
                mreg_proj2(MDSRNN,nties,AcM,MDSDIS2,MDSDIS3);

/**
            printf("MDSDIS2 nach proj: ");
            for (i = 0; i < MDSRNN; ++i)
                printf("%6.3lf ",MDSDIS2[i]);
            newline();
**/


            /* calculate stress */ 
    
            s = t = 0.0;
            for (i = 0; i < MDSRNN; ++i) {
                a = MDSDIS1[i];
                b = MDSDIS2[i];
                s += (a - b) * (a - b);
                t += a * a;
            }
            if (fabs(t) <= EPSI1) {
                printf1("Error in stress calculation.\n");
                goto MDSNFin;
            }
            st = sqrt(s/t);
/** 
            printf("s=%lg t=%lg st=%lg\n",s,t,st);
**/
         
            if (st <= TOLS) {
                printf1("Found configuration with stress: %lg\n",st);
                break;
            }
            mdsn_grad(s,t,st);   

            ngrad = 0.0;
            for (i = 1; i <= MDSNP; ++i)  
                ngrad += MDSG[i] * MDSG[i];
            if (ngrad > 0.0)
                ngrad = sqrt(ngrad);

            printf1("Iteration %3d  stress %lg  norm of grad: %lg\n",iter,st,ngrad);
            if (ngrad <= TOLG) {
                printf1("Convergence reached in %d iterations.\n",iter);
                printf1("Final norm of gradient: %lg\n",ngrad);
                break; 
            }

            /* update configuration */

            a = 0.0;
            for (i = 1; i <= MDSNP; ++i)   
                a += MDSX[i] * MDSX[i];

            if (a <= EPSI1) {
                printf1("Error in distances. Can't continue.\n");
                goto MDSNFin;
            }
            mag = ngrad / sqrt(a);

            printf("mag=%lg\n",mag);

            /***************************************
            if (iter == 1)
                af = gf = 1.0; 
            else {
                gf = dmin(1.0,st/stp);

                a = 0.0;
                for (i = 1; i <= MDSNP; ++i)
                    a += AcY[i] * AcTmp[i];
                cosphi = a / (ngrad * ngradp);
                af = pow(4.0,pow(cosphi,3.0));
            }
            if (iter <= 5) {
                sthis[iter] = st;
                rf = 0.65;
            }
            else {
                st5 = sthis[1];
                for (i = 1; i <= 4; ++i) 
                    sthis[i] = sthis[i+1];
                sthis[5] = st;
                rf = 1.3 / (1.0 + pow(dmin(1.0,st / st5),5.0));
            }
        
            printf("gf=%lg rf=%lg af=%lg cosphi=%lg a=%lg\n",
                                            gf,rf,af,cosphi,a );

            slen *= af * rf * gf;          
            ****************************************************/

            slen = SLen;  
        
/*          printf("SLEN=%lg\n",slen);          */

/**************
            printf("CUR CONF: ");
            for (i = 1; i <= MDSNP; ++i)
                printf("%6.3lf ",MDSX[i]);
            newline();
            printf("NEW CONF: ");
***********/
mag=1.0;
              
            for (i = 1; i <= MDSNP; ++i)  
                MDSX[i] -= slen * MDSG[i] / mag;

            stp = st;                       /* save stress value */
            ngradp = ngrad;                 /* save norm of gradient */
            for (i = 1; i <= MDSNP; ++i)    /* save gradient */
                MDSGP[i] = MDSG[i];                                      

        }
           
        printf1("Best stress value (in repeat %d): %lg\n",ii,st);
                        
        /* print result into output file */
                     
        fprintf(PMFd,"%4d %12.8lf  ",ii,st);
        for (i = 1; i <= MDSNP; ++i)
            fprintf(PMFd,PMFmtS,MDSX[i]);
        fprintf(PMFd,"\n");

        /* save best configuration */

        if (ii == 1) {
            stmin = st;
            for (i = 1; i <= MDSNP; ++i)
                MDSXB[i] = MDSX[i];
        }
        else if (st < stmin) {
            stmin = st;     
            for (i = 1; i <= MDSNP; ++i)
                MDSXB[i] = MDSX[i];
        }
printf("HIER MDSXB: ");
for (i = 1; i <= MDSNP; ++i)
    printf("%lg ",MDSXB[i]);
newline();


        /*******************
        if (DSVFlg)        no repeats if starting values 
            break;
        *******************/
    }

printf("MDSRN=%d MDSND=%d MDSNP=%d\n",MDSRN,MDSND,MDSNP);

    printf1("\nFinal best stress value: %lg\n",stmin);
    printf1("Coordinates\n");
    for (i = 1; i <= MDSRN; ++i) {
        printf1("%3d ",i);
        for (j = 1; j <= MDSND; ++j)  
            printf1(PMFmtS,MDSXB[(i - 1) * MDSND + j]);
        newline();
    }
    newline();

    if (PMF1Def) {      /* write estimated distance matrix */
        /******************
        for (i = 1; i <= MDSRN; ++i) {
            if (i == 1)
                a = b = 0.0;
            else {
                a = AcZ[i - 1];
                b = AcZ[MDSRN + i - 2];
            }
            for (j = 1; j <= MDSRN; ++j) {
                if (j == 1)
                    aa = bb = 0.0;
                else {
                    aa = AcZ[j - 1];
                    bb = AcZ[MDSRN + j - 2];
                }
                d = mds_edis(a,b,aa,bb);
                fprintf(PMF1d,PMFmtS,d);
            }
            fprintf(PMF1d,"\n");
        }
        ******************/
        printf1("%d records written to %s\n",MDSRN,PMF1dName);
    }
    if (PMPCFDef) {     /* write Shepard's diagram */
        /*************************************************
        for (i = 0; i < MDSRNN; ++i) {
            fprintf(PMPCFd,PMFmtS,MDSDIS[MDSPTR[i]]);
            fprintf(PMPCFd,PMFmtS,MDSDIS1[MDSPTR[i]]);
            fprintf(PMPCFd,"\n");
        }
        **********************/
        printf1("%d records written to %s\n",MDSRNN,PMPCFName);
    }
    printf1("%d records written to %s\n",PMNS,PMFdName);

    err = 0;

MDSNFin:       
    ml_init(0,0,0);      
    p_clean();
    return(err);
}



/* ------------------------------------------------------------------------ */
/*  mdsn_dis(dis)       Create distances from configuration MDSX in the     */
/*                      array dis[] in sorted order.                        */

void mdsn_dis(double *dis)
{
    register int i,j,k,l;
    double d,t;

    l = 0;
    for (i = 2; i <= MDSRN; ++i) {
        for (j = 1; j < i; ++j) {
            d = 0.0;
            for (k = 1; k <= MDSND; ++k) {
                t = MDSX[(i - 1) * MDSND + k] - MDSX[(j - 1) * MDSND + k];  
                d += t * t;
            }
            if (d > 0.0)
                d = sqrt(d);

            dis[MDSPTRI[l++]] = d;
        }
    }
}


/* ------------------------------------------------------------------------ */
/*  mdsn_norm()     Normalize the configuration in MDSX and update the      */
/*                  distances in MDSDIS1.                                   */
/*                                                                          */
/*  Return 0 if OK, 1 if error.                                             */

int mdsn_norm(void)
{
    register int i,j;
    double s;             
/**
    printf("X vorher: ");
    for (i = 1; i <= MDSNP; ++i)
        printf("%6.3lf ",MDSX[i]);
    newline();
**/

    for (j = 1; j <= MDSND; ++j) {
        s = 0.0;
        for (i = 1; i <= MDSRN; ++i)   
            s += MDSX[(i - 1) * MDSND + j];
        s /= (double)MDSRN;
        for (i = 1; i <= MDSRN; ++i)   
            MDSX[(i - 1) * MDSND + j] -= s;
    }
/********
    printf("X nach 1: ");
    for (i = 1; i <= MDSNP; ++i)
        printf("%6.3lf ",MDSX[i]);
    newline();
**/
    mdsn_dis(MDSDIS1);    /* calculate distances */


/**
    printf("MDSDIS1 nach 1 sorted: ");
    for (i = 0; i < MDSRNN; ++i)
        printf("%6.3lf ",MDSDIS1[i]);
    newline();

    printf("MDSDIS1 nach 1 unsorted: ");
    for (i = 0; i < MDSRNN; ++i)
        printf("%6.3lf ",MDSDIS1[MDSPTRI[i]]);
    newline();
**/


    s = 0.0;
    for (i = 0; i < MDSRNN; ++i)  
        s += MDSDIS1[i] * MDSDIS1[i];

    if (s <= EPSI1) {
        printf1("Error in distances. Can't continue.\n");
        return(1);          
    }
    s = sqrt((double)MDSRNN / s);

/**
    printf1("s=======%lg\n",s);
**/
    for (i = 1; i <= MDSNP; ++i)
        MDSX[i] *= s;

/**
    printf("X nach 2: ");
    for (i = 1; i <= MDSNP; ++i)
        printf("%6.3lf ",MDSX[i]);
    newline();
**/


    mdsn_dis(MDSDIS1);    /* calculate new distances */

/**

    s = 0.0;
    for (i = 0; i < MDSRNN; ++i)  
        s += MDSDIS1[i] * MDSDIS1[i];

    s = sqrt((double)MDSRNN / s);
    printf1("new s=======%lg\n",s);

    printf("MDSDIS1 nach 2 sorted: ");
    for (i = 0; i < MDSRNN; ++i)
        printf("%6.3lf ",MDSDIS1[i]);
    newline();

    printf("MDSDIS1 nach 2 unsorted: ");
    for (i = 0; i < MDSRNN; ++i)
        printf("%6.3lf ",MDSDIS1[MDSPTRI[i]]);
    newline();


**/


    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mdsn_grad() Calculate gradient of stress function. Assume distances     */
/*              in MDSDIS1, projection values in MDSDIS2. Return gradient   */
/*              in MDSG. Assume stress = sqrt(s/t).                         */
/*              Return 0 if OK, 1 if error.                                 */

int mdsn_grad(double s,double t,double stress)    
{
    register int i,j,k,l,ij,kl;
    double dij,dijp,xil,xjl,g,dki,dkj;

    kl = 0;
    for (k = 1; k <= MDSRN; ++k) {
        for (l = 1; l <= MDSND; ++l) {

            g = 0.0;
            ij = 0;
            for (i = 2; i <= MDSRN; ++i) {
                if (i == k)
                    dki = 1.0;
                else
                    dki = 0.0;
                for (j = 1; j < i; ++j) {
                    if (j == k)
                        dkj = 1.0;
                    else
                        dkj = 0.0;

                    dij = MDSDIS1[MDSPTRI[ij]];
                    dijp = MDSDIS2[MDSPTRI[ij]];

/**
printf("k=%2d l=%2d ij=%3d i=%2d j=%2d dkj=%lg dki=%lg dij=%6.3lf dijp=%6.3lf\n",
            k,l,ij,i,j,dkj,dki,dij,dijp);
**/
                    ij++;

                    if (dij <= 0.0 || dijp <= 0.0)  {
printf("HIER\n");
exit(0);
                        continue;
                    }
                    xil = MDSX[(i - 1) * MDSND + l];
                    xjl = MDSX[(j - 1) * MDSND + l];
                                
                    g += (dki - dkj) * ((dij - dijp) / s - dij / t) *
                            (xil - xjl) / dij;
                }
            }
            MDSG[++kl] = g * stress;
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
       
int mdsn1(void)
{
    register int i,j,k;
    int err,ii,nties,np,iter;
    double d,a,b,aa,bb,dm,fm,ds,st,stp,s,t,stmin,fval,ngrad,ngradp;
    double mag,slen,af,rf,gf,cosphi,st5,sthis[6];

    SLen = 0.2;

/*  TOLG = 1.e-4;   */

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("Nonmetric MDS (experimental). Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto MDSNFin;

    if (gdd_check(PMGN,1))
        goto MDSNFin;
    if (gdd_tcheck(0,1,2))          /* need undirected valued graph */
        goto MDSNFin;

    if (GD_NP < 3) {
        printf1("Error: need at least three nodes.\n");
        goto MDSNFin;
    }
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    if (MxItFlg == 0)
        MxIter = 100;

    MDSRN = GD_NP;                       /* number of points */
    MDSRNN = MDSRN * (MDSRN - 1) / 2;

    /* setup MDSDIS vector with distances */

    if (alloc_acv(MDSRNN + 1))      /* used for MDSDIS */             
        goto MDSNFin; 
    if (alloc_acw(MDSRNN + 1))      /* used for MDSDIS1 */            
        goto MDSNFin; 
    if (alloc_acu(MDSRNN + 1))      /* used for MDSDIS2 */            
        goto MDSNFin; 
    if (alloc_act(MDSRNN + 1))      /* used for MDSDIS3 */            
        goto MDSNFin; 
    if (alloc_acz(MDSRNN + 1))      /* used to save best configuration */
        goto MDSNFin; 
    if (alloc_acs(MDSRNN + 1))      /* used for MDSPTR */             
        goto MDSNFin; 

    MDSDIS = AcV;
    MDSDIS1 = AcW;
    MDSDIS2 = AcU;
    MDSDIS3 = AcT;
    MDSPTR = AcS;

    np = 2 * (MDSRN - 1);       /* number of parameters */

    if (alloc_acx(np + 1))      /* used for configuration (= parameters) */
        goto MDSNFin; 
    if (alloc_acy(np + 1))      /* used for gradient */           
        goto MDSNFin; 
    if (alloc_actmp(np + 1))    /* used for previous gradient */    
        goto MDSNFin; 

    /* put lower triangle of distance matrix into MDSDIS */

    ds = dm = 0.0;
    k = 0;
    for (i = 2; i <= MDSRN; ++i) {
        for (j = 1; j < i; ++j) {
            d = gdd_adj(i-1,j-1,PMGN);
            MDSDIS[k++] = d;
            dm = dmax(dm,d);
            ds += d * d;
        }
    }

    printf1("Largest distance value: %lg\n",dm);

    /* sort MDSDIS in ascending order. pointer: MDSPTR */

    if (sortdp(MDSRNN,MDSDIS,MDSPTR))
        goto MDSNFin;

    printf1("d not sorted: ");
    for (i = 0; i < MDSRNN; ++i)
        printf("%lg ",MDSDIS[i]);
    newline();

    printf1("d sorted: ");
    for (i = 0; i < MDSRNN; ++i)
        printf("%lg ",MDSDIS[MDSPTR[i]]);
    newline();

    printf1("MDSPTR: ");
    for (i = 0; i < MDSRNN; ++i)
        printf("%d ",MDSPTR[i]);
    newline();


    nties = 0;
    a = -1.0;
    for (i = 0; i < MDSRNN; ++i) {
        b = MDSDIS[MDSPTR[i]];
        if (b <= a + EPSI1)
            nties++;
        else
            a = b;
    }
    printf1("Number of ties: %d\n",nties);
    if (nties > 0) {
        printf1("Handling of ties: ");
        if (PMOPT == 2)
            printf1("Kruskal's primary approach.\n");
        else if (PMOPT == 3)
            printf1("Kruskal's secondary approach.\n");
        else {
            PMOPT = 1;
            printf1("ignored.\n");
        }
    }

/**
    nties = mreg_ties(MDSRNN,double *y,int *t);
**/


    /***
    printf("\nNumber of tie blocks m=%d\n",m);
    for (i = 0; i < m; ++i)
        printf("i=%3d t[i]=%d\n",i,t[i]);
    **/

    /*************************************************************************
    if (get_dsv(NParm,Par,0,ParLB,ParUB,1))      try to get starting values   
        goto MDSNFin;
    if (DSVFlg)  
        PMNS = 1;
    else if (PMNS < 1)
        PMNS = 1;
    *****************************/

    if (alloc_act(MDSRNN + 1))      /* used for best coordinates */   
        goto MDSNFin; 

    printf1("\nParameters for iterative minimization.\n");
    printf1("Maximal number of iterations: %d\n",MxIter);
    printf1("Initial step length: %lg\n",SLen);
    printf1("Tolerance for norm of gradient: %lg\n",TOLG);
    printf1("Tolerance for stress: %lg\n\n",TOLS);

    printf1("Number of random repeats: %d\n\n",PMNS);

    for (ii = 1; ii <= PMNS; ++ii) {    /*  repeat */

        printf1("\nRepeat %d\n",ii);

        if (DSVFlg == 0) {              /* get initial configuration */
            for (i = 1; i <= np; ++i) 
                AcX[i] = dm * random1();
        }

        printf("CONFIG: ");
        for (i = 1; i <= np; ++i)
            printf("%lg ",AcX[i]);
        newline();

        slen = SLen;

        for (iter = 1; iter <= MxIter; ++iter) {

            /* create distances of configuration in MDSDIS1 */

            mdsn1_dis(MDSRN,AcX,MDSDIS1);

            /* normalize the configuration */

            a = 0.0;
            for (i = 0; i < MDSRNN; ++i)  
                a += MDSDIS1[i] * MDSDIS1[i];

            if (a <= EPSI1) {
                printf1("Error in distances. Can't continue.\n");
                goto MDSNFin;
            }
            b = sqrt((double)MDSRNN / a);
            printf1("b=======%lg\n",b);

            for (i = 1; i <= np; ++i)
                AcX[i] *= b;

            mdsn1_dis(MDSRN,AcX,MDSDIS1);


            a = 0.0;
            for (i = 0; i < MDSRNN; ++i)  
                a += MDSDIS1[i] * MDSDIS1[i];

            if (a <= EPSI1) {
                printf1("Error in distances. Can't continue.\n");
                goto MDSNFin;
            }
            b = sqrt((double)MDSRNN / a);
            printf1("new b=======%lg\n",b);




            /* find monotone projection in MDSDIS3 */

            mdsn1_proj(PMOPT,MDSRNN,MDSDIS1,MDSDIS2,MDSDIS3);

            /* calculate stress */ 
    
            s = t = 0.0;
            for (i = 0; i < MDSRNN; ++i) {
                a = MDSDIS1[i];
                b = MDSDIS3[i];
                s += (a - b) * (a - b);
                t += a * a;
            }
            if (fabs(t) <= EPSI1) {
                printf1("Error in stress calculation.\n");
                goto MDSNFin;
            }
            st = sqrt(s/t);
    
            printf("s=%lg t=%lg st=%lg\n",s,t,st);

            printf1("Iteration %3d  stress %lg\n",iter,st);
      
            if (st <= TOLS) {
                printf1("Found configuration with stress: %lg\n",st);
                break;
            }
    
            if (mdsn1_fn(MDSRN,MDSRNN,AcX,&fval,AcY,MDSDIS1,MDSDIS3)) {
                printf1("Error in function evaluation. Can't continue.\n");
                goto MDSNFin;
            }
            printf("fval=%lg\n",fval);
          
            printf("GRAD: ");

            ngrad = 0.0;
            for (i = 1; i <= np; ++i) {
                ngrad += AcY[i] * AcY[i];
                printf("%lg ",AcY[i]);
            }
            newline();

            if (ngrad > 0.0)
                ngrad = sqrt(ngrad);

            printf("Norm of gradient: %lg\n",ngrad);

            if (ngrad <= TOLG) {
                printf1("Convergence reached with %d iterations.\n",iter);
                printf1("Final norm of gradient: %lg\n",ngrad);
                break; 
            }

            /* update configuration */

            a = 0.0;
            for (i = 1; i <= np; ++i)   
                a += AcX[i] * AcX[i];

            if (a <= EPSI1) {
                printf1("Error in distances. Can't continue.\n");
                goto MDSNFin;
            }
            mag = ngrad / sqrt(a);

            printf("mag=%lg\n",mag);

            /***************************************
            if (iter == 1)
                af = gf = 1.0; 
            else {
                gf = dmin(1.0,st/stp);

                a = 0.0;
                for (i = 1; i <= np; ++i)
                    a += AcY[i] * AcTmp[i];
                cosphi = a / (ngrad * ngradp);
                af = pow(4.0,pow(cosphi,3.0));
            }
            if (iter <= 5) {
                sthis[iter] = st;
                rf = 0.65;
            }
            else {
                st5 = sthis[1];
                for (i = 1; i <= 4; ++i) 
                    sthis[i] = sthis[i+1];
                sthis[5] = st;
                rf = 1.3 / (1.0 + pow(dmin(1.0,st / st5),5.0));
            }
            printf("gf=%lg rf=%lg af=%lg cosphi=%lg a=%lg\n",
                                            gf,rf,af,cosphi,a );

            slen *= af * rf * gf;          
            ****************************************************/

            slen = SLen;  

            printf("SLEN=%lg\n",slen);

            printf("CUR CONF: ");
            for (i = 1; i <= np; ++i)
                printf("%6.3lf ",AcX[i]);
            newline();
            printf("NEW CONF: ");

            for (i = 1; i <= np; ++i) {
                AcX[i] -= slen * AcY[i] / mag;
                printf("%6.3lf ",AcX[i]);
            }
            newline();

            stp = st;                       /* save stress value */
            ngradp = ngrad;                 /* save norm of gradient */
            for (i = 1; i <= np; ++i)       /* save gradient */
                AcTmp[i] = AcY[i];                                      

        }

        printf1("Best stress value (in repeat %d): %lg\n",ii,st);
                        
        /* print result into output file */
                     
        fprintf(PMFd,"%4d %12.8lf  ",ii,st);
        for (i = 1; i <= np; ++i)
            fprintf(PMFd,PMFmtS,AcX[i]);
        fprintf(PMFd,"\n");

        /* save best configuration */

        if (ii == 1) {
            stmin = st;
            for (i = 1; i <= np; ++i)
                AcZ[i] = AcX[i];
        }
        else if (st < stmin) {
            stmin = st;     
            for (i = 1; i <= NParm; ++i)
                AcZ[i] = AcX[i];
        }
        /*******************
        if (DSVFlg)        no repeats if starting values 
            break;
        *******************/
    }
    printf1("\nFinal best stress value: %lg\n",stmin);
    printf1("Coordinates\n");
    for (i = 1; i <= MDSRN; ++i) {
        printf1("%3d ",i);
        if (i == 1)
            a = b = 0.0;
        else {
            a = AcZ[i - 1];
            b = AcZ[MDSRN + i - 2];
        }
        printf1(PMFmtS,a);
        printf1(PMFmtS,b);
        newline();
    }
    newline();

    if (PMF1Def) {      /* write estimated distance matrix */
        for (i = 1; i <= MDSRN; ++i) {
            if (i == 1)
                a = b = 0.0;
            else {
                a = AcZ[i - 1];
                b = AcZ[MDSRN + i - 2];
            }
            for (j = 1; j <= MDSRN; ++j) {
                if (j == 1)
                    aa = bb = 0.0;
                else {
                    aa = AcZ[j - 1];
                    bb = AcZ[MDSRN + j - 2];
                }
                d = mds_edis(a,b,aa,bb);
                fprintf(PMF1d,PMFmtS,d);
            }
            fprintf(PMF1d,"\n");
        }
        printf1("%d records written to %s\n",MDSRN,PMF1dName);
    }
    if (PMPCFDef) {     /* write Shepard's diagram */
        k = 0;
        for (i = 2; i <= MDSRN; ++i) {
            a = AcX[i - 1];
            b = AcX[MDSRN + i - 2];
            for (j = 1; j < i; ++j) {
                if (j == 1)
                    aa = bb = 0.0;
                else {
                    aa = AcX[j - 1];
                    bb = AcX[MDSRN + j - 2];
                }
                MDSDIS1[k++] = mds_edis(a,b,aa,bb);
            }
        }
        for (i = 0; i < MDSRNN; ++i) {
            fprintf(PMPCFd,PMFmtS,MDSDIS[MDSPTR[i]]);
            fprintf(PMPCFd,PMFmtS,MDSDIS1[MDSPTR[i]]);
            fprintf(PMPCFd,"\n");
        }
        printf1("%d records written to %s\n",MDSRNN,PMPCFName);
    }
    printf1("%d records written to %s\n",PMNS,PMFdName);

    err = 0;

MDSNFin:       
    ml_init(0,0,0);      
    p_clean();
    return(err);
}



/* ------------------------------------------------------------------------ */
/*  mdsn1_dis(n,par,dis) Create distances from configuration in par[]        */
/*                      in the array dis[].                                 */

void mdsn1_dis(int n,double *par,double *dis)
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
            dis[k++] = mds_edis(a,b,aa,bb);
        }
    }

        printf1("x dis not sorted: ");
        for (i = 0; i < k; ++i)
            printf("%lg ",dis[i]);
        newline();
            
        printf1("x dis sorted: ");
        for (i = 0; i < k; ++i)
            printf("%lg ",dis[MDSPTR[i]]);
        newline();

}

/* ------------------------------------------------------------------------ */
/*  mdsn1_proj()     Find projection, depending on opt.                      */
/*                  nn = n (n-1) / 2 (n = number of points)                 */

void mdsn1_proj(int opt,int nn,double *dis,double *tmp,double *proj)
{
    register int i,j,k;


    /* find monotone projection */
    
    for (i = 0; i < nn; ++i)
        tmp[i] = dis[MDSPTR[i]];
    
    if (opt == 1)
        mreg_proj(nn,tmp);
    
    /**********
    else if (PMOPT == 2)
        mreg_proj1(NOC,NTies,AcM,AcI,AcJ,CJRY,AcTmp);
    
    else if (PMOPT == 3)
        mreg_proj2(NOC,NTies,AcM,CJRY,AcTmp);
    **********/
            
    printf1("projection: ");
    for (i = 0; i < nn; ++i)
        printf("%lg ",tmp[i]);
    newline();
    
    /* replace according to original order */
    
    for (i = 0; i < nn; ++i)
        proj[MDSPTR[i]] = tmp[i];
    
    printf1("projection (ORG): ");
    for (i = 0; i < nn; ++i)
        printf("%lg ",proj[i]);
    newline();
}   

/* ------------------------------------------------------------------------ */
/*  mdsn1_fn()  Calculate function value and gradient. n is number of       */
/*              points, nn = n(n-1)/2. par[] contains 2(n-1) parameters.    */  
/*              Return function value in fval, gradient in grad[].          */
/*                                                                          */
/*              Return 0 if OK, 1 if error in function                      */

int mdsn1_fn(int n,int nn,double *par,double *fval,double *grad,double *dis,double *disp)
{
    register int i,j,k,l,il,lj;
    int err;
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
            d = mds_edis(ai,bi,aj,bj);
            MDSDIS1[k] = d;
            d = MDSDIS1[k] - d;
            k++;
        }
    }

    /* calculate stress terms */

    s = t = 0.0;
    for (i = 0; i < nn; ++i) {
        ai = MDSDIS1[i];
        bi = MDSDIS3[i];
        s += (ai - bi) * (ai - bi);
        t += ai * ai;
    }
    if (fabs(s) <= EPSI1 || fabs(t) <= EPSI1) {
        printf1("\nError in function mdsn1_fn().\n");
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
            sal += 2.0 * (al - aj) * (1.0 - MDSDIS3[lj] / MDSDIS1[lj]);
            sbl += 2.0 * (bl - bj) * (1.0 - MDSDIS3[lj] / MDSDIS1[lj]);
        }
        for (i = l + 1; i <= n; ++i) {
            il = ((i - 1) * (i - 2)) / 2 + l - 1;
            ai = par[i - 1];
            bi = par[i + n - 2];                       
            tal -= 2.0 * (ai - al);
            tbl -= 2.0 * (bi - bl);
            sal += 2.0 * (ai - al) * (MDSDIS3[il] / MDSDIS1[il] - 1.0);
            sbl += 2.0 * (bi - bl) * (MDSDIS3[il] / MDSDIS1[il] - 1.0);
        }
        grad[l - 1] += 0.5 * st * (sal / s - tal / t);
        grad[l + n - 2] += 0.5 * st * (sbl / s - tbl / t);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mds_edis(x1,y1,x2,y2) Return euclidean distance of (x1,y1) and (x2,y2). */

double mds_edis(double x1,double y1,double x2,double y2)
{
    return(sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)));
}

/* ------------------------------------------------------------------------ */
/*  mds_elen(x,y)    Return euclidean length of (x,y).                      */

double mds_elen(double x,double y)
{
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
       
int mdsx(void)
{
    register int i,j,ii,jj;
    int err,rn,ip,ip0,ip1,np,i0,nn,sflag;
    double r,dm,dmk,fminp,fmin0,fmin1;

double a,b,d;

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("Multidimensional scaling. Current memory: %d bytes.\n",MemReq);

    TOLBW = 0.01;                   /* default box width */
    TOLFD = 0.1;       

    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto MDSXFin;

    if (gdd_check(PMGN,1))
        goto MDSXFin;
    if (gdd_tcheck(0,1,2))          /* need undirected valued graph */
        goto MDSXFin;

    if (GD_NP < 3) {
        printf1("Error: need at least three nodes.\n");
        goto MDSXFin;
    }
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    if (PMOPT == 1) {
        printf1("Minmax approach.\n\n");
        MDSTYP = 0;
    }
    else {
        MDSTYP = 1;
        printf1("Least squares approach.\n\n");
    }

    MDSRN = GD_NP;                  /* number of points */

    if (MxItFlg == 0)               /* max number of iterations */
        MxIter = 1000;

    /* setup parameters for simplex algorithm */

    NParm = 2 * MDSRN - 3;          /* number of parameters */
    MINA = 2;                       /* always this algorithm */
    PMDOPT = 0;                     /* without derivatives */
    CCTyp = 0;                      /* no covariance matrix */
    set_mlopt();                    /* adjust options */         

    if (ml_init(1,1,0))             /* init function minimization */
        goto MDSXFin;


    /* setup parameters for intervall algorithm */

    if (MxIt1 <= 50)        /* def. max number of iterations */
        MxIt1 = 10000;
    if (PMNBOX < 1)         /* def. max number of boxes */
        PMNBOX = 10000;

    printf1("Parameters for step 0 algorithm.\n");
    printf1("Max number of iterations: %d\n",MxIt1);
    printf1("Maximum number of boxes: %d\n",PMNBOX);
    printf1("Tolerance for box width: %g\n\n",TOLBW);
            
    if (alloc_par(2 * MDSRN,1))         /* allocate arrays for parameters */
        goto MDSXFin;

    if (alloc_list(PMNBOX,2 * MDSRN))   /* allocate boxes */
        goto MDSXFin;

    nn = MDSRN * (MDSRN - 1) / 2;

    if (alloc_acx(MDSRN + 1))       /* first coordinate (a) */
        goto MDSXFin; 
    if (alloc_acy(MDSRN + 1))       /* second coordinate (b) */
        goto MDSXFin; 
    if (alloc_acu(MDSRN + 1))     
        goto MDSXFin; 
    if (alloc_actmp(MDSRN + 1))     
        goto MDSXFin; 
    if (alloc_acn(MDSRN + 1))     
        goto MDSXFin; 
    if (alloc_acm(MDSRN + 1))     
        goto MDSXFin; 
    if (alloc_acxf(MDSRN * MDSRN + 1))     
        goto MDSXFin; 
    if (alloc_aci(nn + 1))     
        goto MDSXFin; 
    if (alloc_acj(nn + 1))     
        goto MDSXFin; 
    if (alloc_acw(nn + 1))     
        goto MDSXFin; 
    if (alloc_ack(nn + 1))     
        goto MDSXFin; 

    MDSRI = AcN;
    MDSRF = AcM;
    MDSRD = AcXF;
                     
    /* find points ip0 and ip1 that have maximal distance */

    dm = 0.0;
    ip0 = ip1 = -1;

    for (i = 0; i < MDSRN; ++i) {
        for (j = 0; j < MDSRN; ++j) {
            if (i != j) {
                r = gdd_adj(i,j,PMGN);
                if (r >= 0.0) {
                    MDSRD[i * MDSRN + j] = r;                  
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
    MDSRIP0 = ip0;
    MDSRIP1 = ip1;

    printf("dm=%g ip0=%d ip1=%d  \n",dm,ip0,ip1);
    for (i = 0; i < MDSRN; ++i) {
        for (j = 0; j < MDSRN; ++j)
            printf("%f ",MDSRD[i * MDSRN + j]);
        newline();
    }
    newline();

    if (ip0 < 0) {
        printf1("Error: no positive distances.\n");
        goto MDSXFin;
    }
    MDSRNI = 2;
    MDSRI[0] = ip0;
    MDSRI[1] = ip1;
    MDSRF[ip0] = 1;
    MDSRF[ip1] = 1;

    AcX[ip0] = 1.0;             /* first point */
    AcY[ip0] = 1.0;

    AcX[ip1] = AcX[ip0] + dm;   /* second point */
    AcY[ip1] = 1.0;
        
    sflag = SILENTFlg;
    SILENTFlg = 2;

    printf1("Step  NPnt  Iterations  NBoxes  Accepted  Variance   Min Function\n");
    if (sflag < 2)
        printfe("\nStep  NPnt  Iterations  NBoxes  Accepted  Variance   Min Function\n");

    /*  do for all further points. fmin0 is the current best minimum */
    /*  found in the first step. */

    while (MDSRNI < MDSRN) {

        ip = 0;              /* find next point */
        dmk = -1.0;
        for (i = 0; i < MDSRN; ++i) {                       
            if (MDSRF[i])
                continue;

            r = 0.0;
            for (j = 0; j < MDSRNI; ++j)  
                r = dmax(r,MDSRD[i * MDSRN + j]);                 
            if (r > dmk) {
                ip = i;
                dmk = r;
            }
        }
        MDSRI[MDSRNI] = MDSRIP = ip;  
        MDSRF[ip] = 1;
        MDSRNI++;

        /* find starting values with interval algorithm */

        RParL[0] = AcX[ip0];
        RParU[0] = AcX[ip1];
        RParL[1] = AcY[ip0];                 
        RParU[1] = AcY[ip0] + dm;                  
        RPar[0] = (RParL[0] + RParU[0]) / 2.0;
        RPar[1] = (RParL[1] + RParU[1]) / 2.0;

        rn = mdsxr_m(2,PMNBOX,RPar,RParL,RParU,AcTmp);             
       
        if (rn) {
            printf1("Error %d in mdsxr_m.\n",rn);
            goto MDSXFin;
        }
       
printf("HIER nach step 0 RPar %g %g\n",RPar[0],RPar[1]);

printf("IT=%d NB=%d NA=%d\n",GO_IT,GO_NB,GO_NA);



        AcX[ip] = RPar[0];
        AcY[ip] = RPar[1];
        fmin0 = mdsxr_fm();

        printf1("%4d  %4d  %10d %7d %9d           %14.7e\n",0,MDSRNI,GO_IT,GO_NB,GO_NA,fmin0);
        if (sflag < 2)
            printfe("%4d  %4d  %10d %7d %9d           %14.7e\n",0,MDSRNI,GO_IT,GO_NB,GO_NA,fmin0);


printf("\nNEU nach step 0 a,b,r  \n");
        for (i = 0; i < MDSRN; ++i)
            printf("%f %f %f\n",AcX[i],AcY[i],AcU[i]);


        rn = mds_m(ip,ip0,ip1);         /* call minimization */
        if (rn) {
            printf1("Error %d in mdsxr_m.\n",rn);
            goto MDSXFin;
        }
        fmin0 = FMin;

        printf1("%4d  %4d  %10d                   %9.2e %14.7e\n",1,MDSRNI,Iter,CValV,fmin0);
        if (sflag < 2)
            printfe("%4d  %4d  %10d                   %9.2e %14.7e\n",1,MDSRNI,Iter,CValV,fmin0);
    }
    SILENTFlg = sflag;
    printf1("\nOptimal value of criterion: %g\n",fmin0);
    if (LConv)  
        printf1("Warning: final step didn't converge.\n");
    newline();

    /*  find optimal radii */

    mds_r(MDSRN,MDSRD,AcX,AcY,AcU,MDSRI,AcW,AcI,AcJ,AcK);

printf("\nNEU Final values a,b,r  \n");
    for (i = 0; i < MDSRN; ++i)
        printf("%f %f %f\n",AcX[i],AcY[i],AcU[i]);

    newline();

    printf1("Neue Abstandsmatrix\n");
    for (i = 0; i < MDSRN; ++i) {
        for (j = 0; j < MDSRN; ++j) {
            a = AcX[i] - AcX[j];
            b = AcY[i] - AcY[j];
            d = sqrt(a * a + b * b);
            printf1("%8.4lf ",d);
        }
        newline();
    }


    for (i = 0; i < MDSRN; ++i) {       /* write to output file */
        fprintf(PMFd,PMNFmtS,i + 1);
        fprintf(PMFd,PMFmtS,AcX[i]);
        fprintf(PMFd,PMFmtS,AcY[i]);
        fprintf(PMFd,PMFmtS,AcU[i]);
        fprintf(PMFd,"\n");
    }
    printf1("%d records written to: %s\n",MDSRN,PMFdName);

    /* create plot command file */

    if (PMPCFDef)       
        mds_p(MDSRN,AcX,AcY,AcU);


    /**********  *********/

    if (PMSC > 0.0) {
        printf1("\nStarting additional algorithm.\n");
        printf1("Initial box width (sc): %g\n",PMSC);
        printf1("Minimal box width (tolfd): %g\n",TOLFD);

        if (alloc_acu(2 * MDSRN + 1))      
            goto MDSXFin; 
        if (alloc_acv(2 * MDSRN + 1))   
            goto MDSXFin; 

        for (i = 0; i < MDSRN; ++i) {
            RPar[i] = AcX[i];
            if (i == ip0)  
                RParL[i] = RParU[i] = RPar[i];
            else {
                RParL[i] = RPar[i] - PMSC / 2.0;
                RParU[i] = RPar[i] + PMSC / 2.0 + EPSI1;
            }     
            j = i + MDSRN;
            RPar[j] = AcY[i];
            if (i == ip0 || i == ip1) 
                RParL[j] = RParU[j] = RPar[j];
            else {
                RParL[j] = RPar[j] - PMSC / 2.0;
                RParU[j] = RPar[j] + PMSC / 2.0 + EPSI1;
            }       
        }
        for (i = 0; i < 2 * MDSRN; ++i)
            printf("i=%3d %20.16f %20.16f %20.16f\n",i,RPar[i],RParL[i],RParU[i]);

/**
        MDSFMIN = fmin0;
        rn = mdsxr_c2(2 * MDSRN,RParL,RParU,AcU,AcV);
**/


        MDSFMIN = PMSC; 
        rn = mdsxr_c3(2 * MDSRN,RParL,RParU,AcU,AcV);
        printf("rn=%d\n",rn);

        for (i = 0; i < 2 * MDSRN; ++i)
            printf("i=%3d %20.16f %20.16f\n",i,AcU[i],AcV[i]);


        printf("MDSNCALL: %d\n",MDSNCALL);
        printf("MDSNSKIP: %d\n",MDSNSKIP);

        if (PMPCFDef)       
            mds_p1(2 * MDSRN,AcU,AcV);

    }
    err = 0;

MDSXFin:       
    alloc_list(0,0);
    alloc_par(0,0);
    ml_init(0,0,0);      
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  mds_m       Function minimization with simplex method, used by mds().   */
/*              Two steps:                                                  */
/*              1) find optimal position for new point given the position   */
/*                 of previous points.                                      */
/*              2) optimize positions of all points.                        */

int mds_m(int ip,int ip0,int ip1)
{
    register int i,ii,k;
    int fflag;
    double d0n,d1n,d01,a,b,a0,a1;

    d0n = MDSRD[ip0 * MDSRN + ip];
    d1n = MDSRD[ip1 * MDSRN + ip];
    d01 = MDSRD[ip0 * MDSRN + ip1];

/**     
    printf("new point is %d\n",ip);
    printf("distance ip0 to new: %f\n",d0n);
    printf("distance ip1 to new: %f\n",d1n);
    printf("distance ip0 to ip1: %f\n",d01);

printf("current XY vor startwert \n");
    for (i = 0; i < MDSRN; ++i)
        printf("%d %f %f\n",i,AcX[i],AcY[i]);
newline();

**/     
         
printf("current XY nach startwert\n");
    for (i = 0; i < MDSRN; ++i)
        printf("%d %f %f\n",i,AcX[i],AcY[i]);
newline();
fflag = 0;
    if (fflag) {
        NParm1 = NParm = 2;
        Par[1] = a;
        Par[2] = b;

        MDSMSTEP = 0;

        if (ffmin(MDSMOD3,1,0,1)) {                               
            printf("ERROR\n");       
            exit(0);    
        }
        AcX[ip] = TPar[1];
        AcY[ip] = TPar[2];

        printf("NEWxPar: %g %g FMIN=%g  \n",Par[1],Par[2],FMin );
    }


    MDSMSTEP = 1;
    NParm1 = NParm = MDSRNI * 2 - 3;

    k = 0;
    for (i = 0; i < MDSRNI; ++i) {
        ii = MDSRI[i];
        if (ii == MDSRIP0)
            continue;
        else if (ii == MDSRIP1)
            Par[++k] = AcX[ii];
        else {
            Par[++k] = AcX[ii];
            Par[++k] = AcY[ii];
        }
    }

/**       
printf("k=%d NParm=%d\n",k,NParm);

    printf("PAR\n");
    for (i = 1; i <= NParm; ++i)
        printf("%d par=%f\n",i,Par[i]);

printf("current XY nach min  \n");
    for (i = 0; i < MDSRN; ++i)
        printf("%d %f %f\n",i,AcX[i],AcY[i]);
newline();

***/     
           
    if (ffmin(MDSMOD3,1,0,1)) {                               
        printf("ERROR\n");       
        exit(0);    
    }

    /* update values in AcX and AcY with new parameters */

    k = 0;
    for (i = 0; i < MDSRNI; ++i) {
        ii = MDSRI[i];
        if (ii == MDSRIP0)
            continue;
        else if (ii == MDSRIP1)
            AcX[ii] = Par[++k];                 
        else {
            AcX[ii] = Par[++k];
            AcY[ii] = Par[++k];
        }
    }
       
printf("NEW XY nach step 1 fmin=%g  \n",FMin );
    for (i = 0; i < MDSRN; ++i)
        printf("%d %f %f\n",i,AcX[i],AcY[i]);
newline();

       

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

int mds3_fn(void)
{
    register int i,j,ii,jj,ip,jp;
    int err;
    double tmp,dij,ai,aj,bi,bj,aij,bij;
      
/**      
    printf("mds3_fn NParm=%d %d \n",NParm,NParm1   );
    printf("PAR: ");
    for (i = 1; i <= NParm; ++i)
        printf("%g ",TPar[i]);
    newline();
**/      

    FTmp = 0.0;

/**
    for (i = 1; i <= NParm; ++i) {
        if (TPar[i] < 10.0 * EPSI1) {
            FTmp = sqrt(DBLMAX);                
printf("HIER\n");
           
            return(0);
        }
    }
**/

    if (MDSMSTEP == 0) {
                      
        for (i = 0; i < MDSRNI; ++i) {
            ip = MDSRI[i];
            if (ip == MDSRIP)
                continue;
            aij = AcX[ip] - TPar[1];
            bij = AcY[ip] - TPar[2];
            dij = MDSRD[ip * MDSRN + MDSRIP];
            tmp = fabs(dij - sqrt(aij * aij + bij * bij));
/**
printf("i=%d ip=%d aij=%g bij=%g ai=%g bi=%g aip=%g bip=%g dij=%g tmp=%g\n",
        i,ip,aij,bij,AcX[ip],AcY[ip],TPar[1],TPar[2],dij ,tmp);
**/
            if (MDSTYP == 0)
                FTmp = dmax(FTmp,tmp);
            else
                FTmp += tmp * tmp;            
        }
    }
    else {
        ii = 1;
        for (i = 1; i < MDSRNI; ++i) {
            ip = MDSRI[i];
            if (ip == MDSRIP0) {
                ai = 1.0;
                bi = 1.0;
            }
            else if (ip == MDSRIP1) {
                ai = TPar[ii++];
                bi = 1.0;
            }
            else {
                ai = TPar[ii++];
                bi = TPar[ii++];
            }
            jj = 1;
            for (j = 0; j < i; ++j) {
                jp = MDSRI[j];
                if (jp == MDSRIP0) {
                    aj = 1.0;
                    bj = 1.0;
                }
                else if (jp == MDSRIP1) {
                    aj = TPar[jj++];
                    bj = 1.0;
                }
                else {
                    aj = TPar[jj++];
                    bj = TPar[jj++];
                }
                dij = MDSRD[ip * MDSRN + jp];
                aij = ai - aj;          
                bij = bi - bj;                   
                tmp = fabs(dij - sqrt(aij * aij + bij * bij));
                if (MDSTYP == 0)
                    FTmp = dmax(FTmp,tmp);
                else
                    FTmp += tmp * tmp;            
            }
        }
    }
    err = checkov();          /* check overflow */
    if (err) {
        printf1("\nError in function evaluation.\n");
        printf1("Numerical overflow.\n");
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  mds_r(n,d,x,y,r,p,ds)                                                   */
/*                                                                          */
/*  Calculate radii.                                                        */
/*                                                                          */
/*  Return 0 if OK, -1 if error                                             */

int mds_r(int n,float *d,double *x,double *y,double *r,int *p,double *ds,
    int *is,int *js,int *ptr)
{
    register int i,j,k,  ii,jj,i0;
    int nn;
    double dij,a,b,tmp,ri,rj,rd;

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

    printf("vor sort\n");
    for (i = 0; i < nn; ++i)
        printf("i=%3d ds=%f is=%d js=%d\n",i,ds[i],is[i],js[i]);

newline();

    if (sortdp(nn,ds,ptr))
        return(-1);

    printf("nach sort\n");
    for (i = 0; i < nn; ++i) {
        k = ptr[i];
        printf("i=%3d ds=%f is=%d js=%d\n",i,ds[k],is[k],js[k]);
    }
newline();

              

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
                rd = dmin(rd,tmp);
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

void mds_p(int n,double *x,double *y,double *r)
{
    register int i,m;
    int nx,ny,icx,icy;
    double xmin,xmax,ymin,ymax,xa,xb,ya,yb,rmin,rmax,tx,ty,scx,scy;  

    rmax = r[0];
    for (i = 1; i < n; ++i)
        rmax = dmax(rmax,r[i]);

    xmin = xmax = x[0];
    ymin = ymax = y[0];
                   
    for (i = 1; i < n; ++i) {
        xmin = dmin(xmin,x[i]);
        xmax = dmax(xmax,x[i]);
        ymin = dmin(ymin,y[i]);
        ymax = dmax(ymax,y[i]);
    }
    rmin = 0.01 * xmax;

    tx = 0.1 * (xmax + xmin + 2.0 * rmax) / 2.0;
    ty = 0.1 * (ymax + ymin + 2.0 * rmax) / 2.0;

    nx = x_getaxval(xmin-tx,xmax+tx,&xa,&xb,0);     /* get values for axes */
    ny = x_getaxval(ymin-ty,ymax+ty,&ya,&yb,0);

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

    fprintf(PMPCFd,"psfile = %s.ps;\n",PMPCFName);
    fprintf(PMPCFd,"psetup(pxlen=100,pylen=100,pxa=%g,%g,pya=%g,%g);\n",xa,xb,ya,yb);
    fprintf(PMPCFd,"plxa(sc=%g,ic=%d);\n",scx,icx);      
    fprintf(PMPCFd,"plya(sc=%g,ic=%d);\n",scy,icy);

    m = 0;
    for (i = 0; i < n; ++i) {
        if (r[i] >= rmin) {
            fprintf(PMPCFd,"ploto(xy=%g,%g,lt=1,nc=1)=%g;\n",x[i],y[i],r[i]);
            m++;
        }
    }
    if (m < n) {
        nx = 0;
        for (i = 0; i < n; ++i) {
            if (r[i] < rmin) {
                if (nx == 0) {
                    fprintf(PMPCFd,"plotp(s=5,fs=1,nc=1,lw=0)=%g,%g",x[i],y[i]);
                    nx = 1;
                }
                else
                    fprintf(PMPCFd,",%g,%g",x[i],y[i]);
            }
        }
        fprintf(PMPCFd,";\n");
    }
    printf1("Output file for plot: %s\n",PMPCFName);
}

/* ------------------------------------------------------------------------ */
/*  mds_p1  Create a command file for plotting the configuration.           */

void mds_p1(int n,double *a,double *b)
{
    register int i,m;
    double x,y,xw,yw;

    m = n / 2;
    for (i = 0; i < m; ++i) {
        x = a[i];
        y = a[i + m];
        xw = b[i] - a[i];
        yw = b[i + m] - a[i + m];
        if (xw <= EPSI1)
            xw = 1000.0 * EPSI1;
        if (yw <= EPSI1)
            yw = 1000.0 * EPSI1;

        fprintf(PMPCFd,"plrec()=%g,%g,%g,%g;\n",x,y,xw,yw);
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

int mdsxr_m(int narg,int nbmax,double *par,double *lb,double *ub,double *par0)
{
    register int i,j;
    int err,r,nb,is,iter,ja,jb,jj,first,nf,na;
    double fl,tmp,w,f,fminp;
    double *lptra,*lptrb,*uptra,*uptrb;

    GO_IT = 0;              /* number of iterations performed */
    GO_FN = 0;              /* number of function evaluations */
    GO_IFN = 0;             /* number of inclusion function evaluations */
    GO_NBU = 0;             /* number of boxes used */
    GO_NB = 0;              /* number of boxes in final list */

    err = 0;
    for (i = 0; i < nbmax; ++i)
        GO_FLG[i] = 0;

    lptra = GO_LB;                      /* put initial box on list */
    uptra = GO_UB;

    for (i = 0; i < narg; ++i) {
        lptra[i] = lb[i];
        uptra[i] = ub[i];
    }
    GO_FLG[0] = -1;     /* not yet processed */

    /* calculate inclusion function with initial boxes */
   
    GO_IFN++;

    r = mdsxr_f(1,narg,lptra,uptra,&fl);
    if (r) {
        printf1("Error in evaluating inclusion function.\n");
        err = -1;
        goto MDSRMFin;
    }
    GO_LBF[0] = fl;
    
    /* calculate function at starting values */

   
    GO_FN++;
    r = mdsxr_f(0,narg,par,par,&fminp);
    if (r) {
        printf1("Error in evaluating function.\n");
        err = -1;
        goto MDSRMFin;
    }
    nb = 1;     /* number of boxes */
    na = 0;     /* number of accepted boxes */
    iter = 0;                                     

    /***
    if (SILENTFlg < 2)
        printfe("\n  Iter    Function Value       NBox  FCall  IFCall\n");
    ***/

    
    while (++iter <= MxIt1) {       

        j = iter / 100;

        /**
        if (SILENTFlg < 2 && 100 * j == iter)
            printfe("%5d  %20.13e  %6d %6d %7d\n",iter,fminp,nb,GO_FN,GO_IFN);
        **/

        /* find box with lowest lower function bound,
           and check for boxes that can be dropped */

        ja = -1;
        first = 1;
        for (i = 0; i < nb; ++i) {

            if (GO_FLG[i]) {    /* check all boxes */
  
                if (GO_LBF[i] > fminp) {
                    GO_FLG[i] = 0;
                }
                else if (GO_FLG[i] < 0) {
                    if (first) {      
                        f = GO_LBF[i];
                        ja = i;
                        first = 0;    
                    }
                    else if (f > GO_LBF[i]) {
                        f = GO_LBF[i];
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
            if (GO_FLG[i] == 0) {
                jb = i;
                break;
            }
        }
        if (jb < 0) {   
            printf1("Exceeded maximum number of boxes.\n");
            err = 1;
            goto MDSRM1Fin; 
        }
        if (jb >= nb) {
            nb = jb + 1;
            if (GO_NBU < nb)
                GO_NBU = nb;
        }
        lptra = GO_LB + ja * narg;
        uptra = GO_UB + ja * narg;

        lptrb = GO_LB + jb * narg;
        uptrb = GO_UB + jb * narg;
   
        is = 0;                         /* coordinate with largest width */
        w = uptra[0] - lptra[0];

        for (i = 1; i < narg; ++i) {
            tmp = uptra[i] - lptra[i];
            if (w < tmp) {
                w = tmp;
                is = i;
            }
        }
        if (w <= TOLBW) {                   /* accept this box */
            GO_FLG[ja] = 1;
            na++;
            continue;
        }

        for (i = 0; i < narg; ++i) {      /* copy box ja to jb */
            lptrb[i] = lptra[i];
            uptrb[i] = uptra[i];
        }

        tmp = lptra[is] + w / 2.0;
        uptra[is] = lptrb[is] = tmp;

        GO_FLG[ja] = -1;                   
        GO_FLG[jb] = -1;

        for (jj = 0; jj < 2; ++jj) {           /* for both boxes */

            if (jj)   
                ja = jb;

            lptra = GO_LB + ja * narg;
            uptra = GO_UB + ja * narg;

            /* evaluate inclusion function */

            GO_IFN++;
            r = mdsxr_f(1,narg,lptra,uptra,&fl);
            if (r) {
                printf1("Error in evaluating inclusion function.\n");
                err = -1;
                goto MDSRMFin;
            }
            if (fl > fminp) {
                GO_FLG[ja] = 0;             /* free box */
                if (ja == nb - 1)
                    nb--;    
                continue;
            }
            GO_LBF[ja] = fl;

            /* update fminp with function value for midpoint */             
                   
            for (j = 0; j < narg; ++j)
                AcTmp[j] = (lptra[j] + uptra[j]) / 2.0;

            GO_FN++;
            r = mdsxr_f(0,narg,AcTmp,AcTmp,&f);
            if (r) {
                printf1("Error in evaluating function.\n");
                err = -1;
                goto MDSRMFin;
            }
            if (fminp > f) {
                fminp = f;
                igmin_cpar(narg,AcTmp,par);
            }   
        }
    }
MDSRM1Fin:
    GO_NB = nb;
    GO_NA = na;
    GO_IT = iter;
    GO_FMIN = fminp;

    if (err == 0) {
        if (Iter >= MxIt1)
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

int mdsxr_f(int opt,int n,double *xl,double *xh,double *f)
{
    register int i,j,ii,jj;
    int first;
    double dij,ai,bi,aj,bj,aij,bij,al,ah,bl,bh,sl,sh,tl,th,ul,uh,tmp;
    double aih,bih,ajh,bjh,ri;

/***************** 
    printf("MDSR_F n=%d MDSRIP=%d    \n",n,MDSRIP       );    
    for (i = 0; i < n; ++i)
        printf1("i=%d xl=%f xh=%f\n",i,xl[i],xh[i]);
**************/

    *f = 0.0;
    if (opt == 0) {
        ai = xl[0];
        bi = xl[1];
        for (j = 0; j < MDSRNI; ++j) {
            jj = MDSRI[j];
            if (jj == MDSRIP)  
                continue;

            dij = MDSRD[MDSRIP * MDSRN + jj];
            if (dij < 0.0)
                continue;
            aj = AcX[jj];
            bj = AcY[jj];
            aij = ai - aj;
            bij = bi - bj;
            tmp = fabs(dij - sqrt(aij * aij + bij * bij));
            if (MDSTYP == 0)
                *f = dmax(*f,tmp);
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

    for (j = 0; j < MDSRNI; ++j) {
        jj = MDSRI[j];
        if (jj == MDSRIP)
            continue;

        dij = MDSRD[MDSRIP * MDSRN + jj];
        if (dij < 0.0)
            continue;

        aj = ajh = AcX[jj];
        bj = bjh = AcY[jj];

        i_sub(ai,aih,aj,ajh,&sl,&sh);       
        i_square(sl,sh,&tl,&th);

        i_sub(bi,bih,bj,bjh,&sl,&sh);       
        i_square(sl,sh,&ul,&uh);

        i_add(tl,th,ul,uh,&sl,&sh);       
        i_sqrt(sl,sh,&tl,&th);
        i_sub(dij,dij,tl,th,&sl,&sh);
        i_abs(sl,sh,&tl,&th);

        if (MDSTYP == 0)
            *f = dmax(*f,tl);
        else {
            i_square(tl,th,&ul,&uh);
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

double mdsxr_fm(void)
{
    register int i,j,ii,jj;
    double f,dij,ai,bi,aj,bj,aij,bij,tmp;

    f = 0.0;
    for (i = 1; i < MDSRNI; ++i) {
        ii = MDSRI[i];
        ai = AcX[ii];
        bi = AcY[ii];
        for (j = 0; j < i; ++j) {
            jj = MDSRI[j];
            dij = MDSRD[ii * MDSRN + jj];
            if (dij < 0.0)
                continue;

            aj = AcX[jj];
            bj = AcY[jj];
            aij = ai - aj;
            bij = bi - bj;
            tmp = fabs(dij - sqrt(aij * aij + bij * bij));
            if (MDSTYP == 0)
                f = dmax(f,tmp);
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

int mdsxr_f1(int opt,int n,double *xl,double *xh,double *f,double fminp)
{
    register int i,j;
    double aij,bij,dij,tmp,ail,aih,bil,bih,ajl,ajh,bjl,bjh,d,d1,d2,d3,d4;

/** 
        printf("MDSR_F n=%d MDSRIP=%d    \n",n,MDSRIP       );    
        for (i = 0; i < n; ++i)
            printf1("i=%d xl=%f xh=%f\n",i,xl[i],xh[i]);
**/ 

    *f = 0.0;
    if (opt == 0) {
        for (i = 1; i < MDSRN; ++i) {
            for (j = 0; j < i; ++j) {
                dij = MDSRD[i * MDSRN + j];
                if (dij < 0.0)
                    continue;
                aij = xl[i] - xl[j];
                bij = xl[i + MDSRN] - xl[j + MDSRN];
                tmp = fabs(dij - sqrt(aij * aij + bij * bij));
                *f = dmax(*f,tmp);
            }
        }
        return(0);
    }

    for (i = 1; i < MDSRN; ++i) {
        ail = xl[i];
        aih = xh[i];
        bil = xl[i + MDSRN];
        bih = xh[i + MDSRN];

        for (j = 0; j < i; ++j) {
            ajl = xl[j];
            ajh = xh[j];
            bjl = xl[j + MDSRN];
            bjh = xh[j + MDSRN];

            dij = MDSRD[i * MDSRN + j];
            if (dij < 0.0)
                continue;

            tmp = 0.0;
            if (aih <= ajl) {
                if (fabs(bih - bjl) >= fabs(bil - bjh))  
                    d = mdsxr_dm(ail,bih,ajh,bjl);
                else
                    d = mdsxr_dm(ail,bil,ajh,bjh);

                if (dij > d) {
                    tmp = dij - d;
                }
                else {
                    if (bil > bjh || bih < bjl) {
                        if (fabs(bil - bjh) <= fabs(bih - bjl))  
                            d = mdsxr_dm(aih,bil,ajl,bjh);
                        else
                            d = mdsxr_dm(aih,bih,ajl,bjl);
                    }   
                    else
                        d = mdsxr_dm(aih,0.0,ajl,0.0);

                    if (dij < d)
                        tmp = d - dij;
                }
            }
            else if (ajh <= ail) {
                if (fabs(bih - bjl) >= fabs(bil - bjh))  
                    d = mdsxr_dm(aih,bih,ajl,bjl);
                else
                    d = mdsxr_dm(aih,bil,ajl,bjh);

                if (dij > d) {
                    tmp = dij - d;
                }
                else {
                    if (bil > bjh || bih < bjl) {
                        if (fabs(bil - bjh) <= fabs(bih - bjl))  
                            d = mdsxr_dm(ail,bil,ajh,bjh);
                        else
                            d = mdsxr_dm(ail,bih,ajh,bjl);
                    }   
                    else
                        d = mdsxr_dm(ail,0.0,ajh,0.0);

                    if (dij < d)
                        tmp = d - dij;
                }
            }
            else if (bih <= bjl) {
                if (fabs(aih - ajl) > fabs(ail - ajh))  
                    d = mdsxr_dm(aih,bil,ajl,bjh);
                else
                    d = mdsxr_dm(ail,bil,ajh,bjh);

                if (dij > d) {
                    tmp = dij - d;
                }
                else {
                    d = mdsxr_dm(0.0,bih,0.0,bjl);
                    if (dij < d)
                        tmp = d - dij;
                }
            }
            else if (bil >= bjh) {
                if (fabs(aih - ajl) > fabs(ail - ajh))  
                    d = mdsxr_dm(aih,bih,ajl,bjl);
                else
                    d = mdsxr_dm(ail,bih,ajh,bjl);

                if (dij > d) {
                    tmp = dij - d;
                }
                else {
                    d = mdsxr_dm(0.0,bil,0.0,bjh);
                    if (dij < d)
                        tmp = d - dij;
                }
            }
            else {
                d1 = mdsxr_dm(ail,bil,ajh,bjh);
                d2 = mdsxr_dm(aih,bil,ajl,bjh);
                d3 = mdsxr_dm(aih,bih,ajl,bjl);
                d4 = mdsxr_dm(ail,bih,ajh,bjl);
                d = dmax(dmax(d1,d2),dmax(d3,d4));
                if (dij > d) {
                    tmp = dij - d;
                }
                else {
                    ;
                }
            }
            *f = dmax(*f,tmp);
            if (*f > fminp)
                return(0);
        }
    }
    return(0);
}

double mdsxr_dm(double ai,double bi,double aj,double bj)
{
    double aij,bij;

    aij = ai - aj;
    bij = bi - bj;
    return(sqrt(aij * aij + bij * bij));
}

/* ------------------------------------------------------------------------ */
/*  mdsxr_c2()                                                               */
/*                                                                          */

int mdsxr_c2(int n,double *lb,double *ub,double *x,double *y)
{
    register int i;
    int is;
    double a,b,c,f,w,tmp;

    MDSNCALL++;

    mdsxr_f1(1,n,lb,ub,&f,MDSFMIN);
    if (f > MDSFMIN) {
        MDSNSKIP++;
/*      printf("skipped\n");        */     
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
/*   printf("w=%g is=%d  \n",w,is);        */
          
    if (w > TOLFD) {
        a = lb[is];
        b = ub[is];
        c = (a + b) / 2.0;

        ub[is] = c;
        mdsxr_c2(n,lb,ub,x,y);

        ub[is] = b;
        lb[is] = c;

        mdsxr_c2(n,lb,ub,x,y);

        lb[is] = a;
    }
    else {
                              
        if (f < MDSFMIN) {
            printf("1-MDSFMIN=%g <- f=%g TOLBW=%g w=%g \n",MDSFMIN,f,TOLFD,w  );
            MDSFMIN = f;

            for (i = 0; i < n; ++i) {
                x[i] = lb[i];
                y[i] = ub[i];

                printf("i=%3d %20.16f %20.16f\n",i,lb[i],ub[i]);
            }
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mdsxr_c3()                                                               */
/*                                                                          */

int mdsxr_c3(int n,double *lb,double *ub,double *x,double *y)
{
    register int i;
    int is;
    double a,b,c,f,w,tmp;

    MDSNCALL++;
         
    mdsxr_f1(1,n,lb,ub,&f,TOLFE);

printf("f=%g TOLFE=%g\n",f,TOLFE);
    if (f > TOLFE) {
printf("skipped\n");
        MDSNSKIP++;
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


    printf("w=%g is=%d MDSFMIN=%g \n",w,is,MDSFMIN);  
          
    if (w >= MDSFMIN) {
        a = lb[is];
        b = ub[is];
        c = (a + b) / 2.0;

        ub[is] = c;
        mdsxr_c3(n,lb,ub,x,y);

        ub[is] = b;
        lb[is] = c;

        mdsxr_c3(n,lb,ub,x,y);

        lb[is] = a;
    }
    else if (w < MDSFMIN) {
printf("hier\n");             
        if (f < TOLFE) {
            printf("1-MDSFMIN=%g <- w=%g TOLBW=%g f=%g TOLFE=%g\n",MDSFMIN,w,TOLFD,f,TOLFE);
            MDSFMIN = w;

            for (i = 0; i < n; ++i) {
                x[i] = lb[i];
                y[i] = ub[i];

                printf("i=%3d %20.16f %20.16f\n",i,lb[i],ub[i]);
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

int mdsr(void)
{
    register int i;
    int err,ix,iy,iv,nn,vz,vzmax;
    double phi,sinphi,cosphi,r,phimax,rmax,xmin,xmax,ymin,ymax,xmean,ymean;
    double d,ya,yb;

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("Find axes in configuration. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,4,1))     /* get parameters */
        goto MDSRFin;

    if (PMFmtF == 0)
        pmfmt(10,4);

    if (PMNV != 2) {
        printf1("Error: need two variables on right-hand side.\n");
        goto MDSRFin;
    }
    if (PM1NV < 1) {
        printf1("Error: no xv variables.\n");
        goto MDSRFin;
    }
    ix = PMVIdx[0];
    iy = PMVIdx[1];
    printf1("\nConfiguration: %s, %s\n",VName[ix],VName[iy]);

    if (alloc_acx(NOC + 1))  
        goto MDSRFin;
    if (alloc_acy(NOC + 1))  
        goto MDSRFin;

    nn = 1.0;
    if (NOC > 1)
        nn = NOC * (NOC - 1) / 2.0;

    xmin = xmax = get_data(ix,0);
    ymin = ymax = get_data(ix,1);
    xmean = xmin;       
    ymean = ymin;

    for (i = 1; i < NOC; ++i) {
        r = get_data(ix,i);
        xmean += r;
        xmin = dmin(xmin,r);
        xmax = dmax(xmax,r);
        r = get_data(iy,i);
        ymean += r;
        ymin = dmin(ymin,r);
        ymax = dmax(ymax,r);
    }
    xmean /= (double)NOC;
    ymean /= (double)NOC;

    printf1("X minimum: ");
    printf1(PMFmtS,xmin);
    printf1("  maximum: ");
    printf1(PMFmtS,xmax);
    printf1("  mean: ");
    printf1(PMFmtS,xmean);
    printf1("\nY minimum: ");
    printf1(PMFmtS,ymin);
    printf1("  maximum: ");
    printf1(PMFmtS,ymax);
    printf1("  mean: ");
    printf1(PMFmtS,ymean);
    newline();

    if (fabs(xmin - xmax) <= EPSI1)
        goto MDSRFin; 

    for (iv = 0; iv < PM1NV; ++iv) {
        printf1("\nVariable: %s\n",VName[PM1VIdx[iv]]);

        for (i = 0; i < NOC; ++i)  
            AcY[i] = get_data(PM1VIdx[iv],i);

        rmax = phimax = 0.0;
        for (phi = 0.0; phi < Pi; phi += 0.01) {
            sinphi = sin(phi);
            cosphi = cos(phi);
            for (i = 0; i < NOC; ++i)  
                AcX[i] = cosphi * get_data(ix,i) + sinphi * get_data(iy,i);

            r = mdsr_f(NOC,AcX,AcY);
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
        printf1("Optimal phi: %lg\n",phimax);
        printf1("Rank correlation: %lg\n",rmax * (double)vzmax);

        phi = phimax;
        sinphi = sin(phi);
        cosphi = cos(phi);

        /***************************
        for (i = 0; i < NOC; ++i) {
            printf1(PMFmtS,cosphi * get_data(ix,i) + sinphi * get_data(iy,i));
            printf1(PMFmtS,-sinphi * get_data(ix,i) + cosphi * get_data(iy,i));
            newline();
        }
        *********************************/

        d = sinphi * (xmean - xmin);
        ya = ymean - d;
        d = sinphi * (xmax - xmean);
        yb = ymean + d;

        printf1("Axis: (%lg,%lg) to (%lg,%lg)\n",xmin,ya,xmax,yb);
    }
    err = 0;

MDSRFin:
    p_clean();
    return(err);
}


/* ------------------------------------------------------------------------ */
/*  mdsr_f  return rank correlation                                         */
/*                                                                          */

double mdsr_f(int n,double *x,double *y)
{
    register int i,j;
    double r;                  

    r = 0.0;
    for (i = 1; i < n; ++i) {
        for (j = 0; j < i; ++j)  
            r += dsign(x[i] - x[j]) * dsign(y[i] - y[j]);
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
       
int rfit(void)
{
    register int i,j,k,ii;
    int err,n,nn,nh,nt,m,r,vmax,ncon;
    int *ap;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Fit of binary relations. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto RFITFin;

    if (gdd_check(PMGN,1))
        goto RFITFin;
    if (gdd_tcheck(0,0,2))          /* need valued graph */
        goto RFITFin;

    if (PMMax < 1)
        PMMax = 10;

    n = GD_NP;         
    nn = n * n;
    nh = (n * (n - 1)) / 2;

    printf1("\nNumber of nodes: %d\n",n);
    printf1("Entries in adjacency matrix: %d\n",nn);

    ncon = 0;                            /* number of constraints */
    printf1("\nRequested properties:\n");
    if (PMREL[1] == 1) {
        printf1("- reflexive\n");
        ncon += n;
    }
    if (PMREL[2] == 1) {
        printf1("- symmetric\n");
        ncon += 2 * nh;
    }
    if (PMREL[3] == 1) {
        printf1("- anti-symmetric\n");
        ncon += nh;
    }
    if (PMREL[4] == 1) {
        printf1("- transitive\n");
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
    if (PMREL[5] == 1) {
        printf1("- complete\n");
        ncon += nh;
    }
    printf1("\nNumber of constraints: %d\n",ncon);

    m = ncon + 1;

    if (alloc_acr(m * nn + 1))  
        goto RFITFin;

    if (alloc_acs(m + 1))  
        goto RFITFin;

    k = 0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (i == j)
                r = 1;
            else
                r = (int)gdd_adj(i,j,PMGN);                   
            AcR[++k] = 2 * r - 1;                                    
        }
    }

    /* add constraints */

    ii = 2;
    if (PMREL[1] == 1) {               /* reflexive */
        for (i = 1; i <= n; ++i) {                                         
            AcR[(ii - 1) * nn + (i - 1) * n + i] = 1;
            AcS[ii] = 1;
            ii++;
        }
    }
    if (PMREL[2] == 1) {               /* symmetric */
        for (i = 1; i < n; ++i) {                                          
            for (j = i + 1; j <= n; ++j) {
                AcR[(ii - 1) * nn + (i - 1) * n + j] =  1;
                AcR[(ii - 1) * nn + (j - 1) * n + i] = -1;
                AcS[ii] = 0;
                ii++;
                AcR[(ii - 1) * nn + (i - 1) * n + j] = -1;
                AcR[(ii - 1) * nn + (j - 1) * n + i] =  1;
                AcS[ii] = 0;
                ii++;
            }
        }
    }
    if (PMREL[3] == 1) {               /* anti-symmetric */
        for (i = 1; i < n; ++i) {                                          
            for (j = i + 1; j <= n; ++j) {
                AcR[(ii - 1) * nn + (i - 1) * n + j] = -1;
                AcR[(ii - 1) * nn + (j - 1) * n + i] = -1;
                AcS[ii] = -1;
                ii++;
            }
        }
    }
    if (PMREL[4] == 1) {               /* transitive */
        for (i = 1; i <= n; ++i) {                                    
            for (j = 1; j <= n; ++j) {
                if (j == i)
                    continue;
                for (k = 1; k <= n; ++k) { 
                    if (k == i || k == j)
                        continue;
                    AcR[(ii - 1) * nn + (i - 1) * n + k] =  1;
                    AcR[(ii - 1) * nn + (i - 1) * n + j] = -1;
                    AcR[(ii - 1) * nn + (j - 1) * n + k] = -1;
                    AcS[ii] = -1;
                    ii++;
                }
            }
        }
    } 
    if (PMREL[5] == 1) {               /* complete */
        for (i = 1; i < n; ++i) {                                        
            for (j = i + 1; j <= n; ++j) {
                AcR[(ii - 1) * nn + (i - 1) * n + j] = 1;
                AcR[(ii - 1) * nn + (j - 1) * n + i] = 1;
                AcS[ii] = 1;
                ii++;
            }
        }
    }
    if (PMF1Def) {
        for (i = 0; i < m; ++i) {
            for (j = 1; j <= nn; ++j) 
                fprintf(PMF1d,"%2d ",AcR[i * nn + j]);
            fprintf(PMF1d,"%3d\n",AcS[i + 1]);
        }
        printf1("%d records written to: %s\n",m,PMF1dName);
    }
    newline();

    if (alloc_acn(PMMax * nn + 1))
        goto RFITFin;

    if (alloc_acm(4 * m + 6 * nn + m * nn + 2))
        goto RFITFin;

    r = lpi(nn,m,PMMax,AcR,AcS,AcN,AcM,&vmax);

    if (r < 1) {
        switch (r) {
            case -1:    printf1("Error: inconsistent constraints.\n");
                        break;
            case -2:    printf1("Error: exceeded maximal number of solutions.\n");
                        break;
            default:    printf1("Error (%d).\n",r);
                        break;
        }   
        goto RFITFin;
    }
    printf1("Value: %d\nNumber of solutions: %d\n",vmax,r);

    nt = 0;
    for (ii = 0; ii < r; ++ii) {

        for (j = 1; j <= nn; ++j)
            fprintf(PMFd,"%d ",AcN[ii * nn + j]);
        fprintf(PMFd,"\n");
        nt++;
        ap = AcN + ii * nn;
        k = 1;
        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= n; ++j)
                fprintf(PMFd,"%d ",ap[k++]);
            fprintf(PMFd,"\n");
            nt++;
        }
        fprintf(PMFd,"\n");
        nt++;
    }
    printf1("%d records written to: %s\n",nt,PMFdName);
    err = 0;

RFITFin:
    p_clean();
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
       
int rfit1(void)
{
    register int i,j,k,ii;
    int err,n,nn,nh,nt,m,r,vmax,ncon;
    int *ap;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Fit of binary relations. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,1,1))       /* get parameters */
        goto RFITFin;

    if (gdd_check(PMGN,1))
        goto RFITFin;
    if (gdd_tcheck(0,0,2))          /* need valued graph */
        goto RFITFin;

    if (PMMax < 1)
        PMMax = 10;

    n = GD_NP;         
    nn = n * n;
    nh = (n * (n - 1)) / 2;

    printf1("\nNumber of nodes: %d\n",n);
    printf1("Entries in adjacency matrix: %d\n",nn);

    ncon = 0;                            /* number of constraints */
    printf1("\nRequested properties:\n");
    if (PMREL[1] == 1) {
        printf1("- reflexive\n");
        ncon += n;
    }
    if (PMREL[2] == 1) {
        printf1("- symmetric\n");
        ncon += 2 * nh;
    }
    if (PMREL[3] == 1) {
        printf1("- anti-symmetric\n");
        ncon += nh;
    }
    if (PMREL[4] == 1) {
        printf1("- transitive\n");
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
    if (PMREL[5] == 1) {
        printf1("- complete\n");
        ncon += nh;
    }
    printf1("\nNumber of constraints: %d\n",ncon);

    m = ncon + 1;

    if (alloc_acr(m * nn + 1))  
        goto RFITFin;

    if (alloc_acs(m + 1))  
        goto RFITFin;

    k = 0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (i == j)
                r = 1;
            else
                r = (int)gdd_adj(i,j,PMGN);                   
            AcR[++k] = 2 * r - 1;                                    
        }
    }

    /* add constraints */

    ii = 2;
    if (PMREL[1] == 1) {               /* reflexive */
        for (i = 1; i <= n; ++i) {                                         
            AcR[(ii - 1) * nn + (i - 1) * n + i] = 1;
            AcS[ii] = 1;
            ii++;
        }
    }
    if (PMREL[2] == 1) {               /* symmetric */
        for (i = 1; i < n; ++i) {                                          
            for (j = i + 1; j <= n; ++j) {
                AcR[(ii - 1) * nn + (i - 1) * n + j] =  1;
                AcR[(ii - 1) * nn + (j - 1) * n + i] = -1;
                AcS[ii] = 0;
                ii++;
                AcR[(ii - 1) * nn + (i - 1) * n + j] = -1;
                AcR[(ii - 1) * nn + (j - 1) * n + i] =  1;
                AcS[ii] = 0;
                ii++;
            }
        }
    }
    if (PMREL[3] == 1) {               /* anti-symmetric */
        for (i = 1; i < n; ++i) {                                          
            for (j = i + 1; j <= n; ++j) {
                AcR[(ii - 1) * nn + (i - 1) * n + j] = -1;
                AcR[(ii - 1) * nn + (j - 1) * n + i] = -1;
                AcS[ii] = -1;
                ii++;
            }
        }
    }
    if (PMREL[4] == 1) {               /* transitive */
        for (i = 1; i <= n; ++i) {                                    
            for (j = 1; j <= n; ++j) {
                if (j == i)
                    continue;
                for (k = 1; k <= n; ++k) { 
                    if (k == i || k == j)
                        continue;
                    AcR[(ii - 1) * nn + (i - 1) * n + k] =  1;
                    AcR[(ii - 1) * nn + (i - 1) * n + j] = -1;
                    AcR[(ii - 1) * nn + (j - 1) * n + k] = -1;
                    AcS[ii] = -1;
                    ii++;
                }
            }
        }
    } 
    if (PMREL[5] == 1) {               /* complete */
        for (i = 1; i < n; ++i) {                                        
            for (j = i + 1; j <= n; ++j) {
                AcR[(ii - 1) * nn + (i - 1) * n + j] = 1;
                AcR[(ii - 1) * nn + (j - 1) * n + i] = 1;
                AcS[ii] = 1;
                ii++;
            }
        }
    }
    if (PMF1Def) {
        for (i = 0; i < m; ++i) {
            for (j = 1; j <= nn; ++j) 
                fprintf(PMF1d,"%2d ",AcR[i * nn + j]);
            fprintf(PMF1d,"%3d\n",AcS[i + 1]);
        }
        printf1("%d records written to: %s\n",m,PMF1dName);
    }
    newline();

    if (alloc_ack(nn + 1))
        goto RFITFin;
    if (alloc_acl(nn + 1))
        goto RFITFin;
    if (alloc_aci(nn + 1))
        goto RFITFin;
    if (alloc_acj(nn + 1))
        goto RFITFin;
    if (alloc_aciptr(nn + 1))
        goto RFITFin;
    if (alloc_acjptr(nn + 1))
        goto RFITFin;

    for (j = 1; j <= nn; ++j) {
        printf("%5d : ",j);
        nt = 0;
        for (i = 0; i < m; ++i) {
            if (AcR[i * nn + j] != 0)
                nt++;
        }
        printf("%d\n",nt);

        AcK[j] = AcL[j] = nt;

        if (!(AcIPtr[j] = (int *)calloc(2 * nt + 2,sizeof(int)))) {
            p_err(-2,1);
            goto RFITFin;
        }
        AcI[j] = 2 * nt + 2;
        memrq(2 * nt + 2,sizeof(int));
        if (!(AcJPtr[j] = (int *)calloc(2 * nt + 2,sizeof(int)))) {
            p_err(-2,1);
            goto RFITFin;
        }
        AcJ[j] = 2 * nt + 2;
        memrq(2 * nt + 2,sizeof(int));

        ii = nt + 1;
        for (i = 0; i < m; ++i) {
            if ((k = AcR[i * nn + j]) != 0) {
                AcIPtr[j][ii] = k;
                AcJPtr[j][ii] = i + 1;
                ii++;
            }
        }
    }
    alloc_acr(0);       /* no longer required */

    if (alloc_acn(PMMax * nn + 1))
        goto RFITFin;

    if (alloc_acm(4 * m + 6 * nn + 2000            ))
        goto RFITFin;

/**
    r = lpia(nn,m,PMMax,AcR,AcS,AcN,AcM,&vmax);               
    r = 0;

int lpia(int n,int m,int nest,int *an,int *an0,int **av,int **ar,
         int *b0,int *opts,int *itmp,int *vmax);


**/
  
    r = lpia(nn,m,PMMax,AcK,AcL,AcIPtr,AcJPtr,AcS,AcN,AcM,&vmax);
  
     

    if (r < 1) {
        switch (r) {
            case -1:    printf1("Error: inconsistent constraints.\n");
                        break;
            case -2:    printf1("Error: exceeded maximal number of solutions.\n");
                        break;
            default:    printf1("Error (%d).\n",r);
                        break;
        }   
        goto RFITFin;
    }
    printf1("Value: %d\nNumber of solutions: %d\n",vmax,r);

    nt = 0;
    for (ii = 0; ii < r; ++ii) {

        for (j = 1; j <= nn; ++j)
            fprintf(PMFd,"%d ",AcN[ii * nn + j]);
        fprintf(PMFd,"\n");
        nt++;
        ap = AcN + ii * nn;
        k = 1;
        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= n; ++j)
                fprintf(PMFd,"%d ",ap[k++]);
            fprintf(PMFd,"\n");
            nt++;
        }
        fprintf(PMFd,"\n");
        nt++;
    }
    printf1("%d records written to: %s\n",nt,PMFdName);
    err = 0;

RFITFin:
    for (j = 1; j <= nn; ++j) {
        if (AcI[j] > 0) {
            free((char *)AcIPtr[j]);
            memrq(-AcI[j],sizeof(int));
            AcI[j] = 0;        
        }
        if (AcJ[j] > 0) {
            free((char *)AcJPtr[j]);
            memrq(-AcJ[j],sizeof(int));
            AcI[j] = 0;        
        }
    }
    p_clean();
    return(err);
}








