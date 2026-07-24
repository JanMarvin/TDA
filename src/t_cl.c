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

/*  functions in t_cl.c */

int hcls(void);
int h_cls(int ctyp,int gn,int n,float *d,int *a,int *b,int *p,int *q,float *h);
float cl_getdu(int i,int j,float *d);
void cl_putdu(float x,int i,int j,float *d);
void h_cls_den(int n,int *a,int *b,float *h,int *p,int *q,float *x,float *y);
void h_cls_el(int nmax,int n,int *a,int *b,float *h,int *new,float *lev);
void h_cls_dm(int n,int *a,int *b,float *h,int *new,int *ll,float *d);
double h_cls_mu(int n,float *d);

int nncl(void);
int nn_cl(int *cn,double t);
int mn_cl(int k,int *idx,int *mptr,int *iptr,int *jptr,double *val);
int mn_cl1(int i,int k,int *nptr,double *val,int *sptr,int *idx);
int hcld(void);
int h_cld1(int n,int min,int *cidx,int *ncn,int *ncnt,int opt);
int h_cld2(int n,int nsmax,int *cidx,double *dia,int *nodes,            
    int *in,int *jn,int opt);
void h_cld_prn(int n,int *cidx,int cn,int k,int nn,double *dia);
double h_cld_dia(int n,int *cidx,int k);
int becl(void);
int be_clprn(int opt,int n,int *idx,double mval);
double be_cl1(int opt,int mflag,int i0,int n,int *pos,int *idx,int *tidx);
double be_ccheck(int n,int *idx);                         
double be_rcheck(int n,int *idx);                         
int acl(void);
float acl1(int n,int m,float *d,int *p,int *q,float mev,float cmin);

void udcl(int n,int nc,float *d,int *p,int l,int u,float *c,double *pc,
    int *qa,int *qb,int *fp);
void udcl1(int n,float *d,int *p,float *c,int u);
int ucl(void);
int hclsp(void);
void hclsp_visit(int i,int gn);
void hclsp_prn(int i,int gn,int nlev,int level,int mlen);
void hclsp_prn1(int i,int gn,int m);
int scla(void);
int scla_check(int i,int n,int *acn,int *p,double *dp);
int scla_find(int n,int *acn,int np,int *p,int *com);
int clp(void);
int clp_f(int meth,int m,int n,int itmax,
          double *t,int *z,int *mj,double *e,double *c);
int clp_p(void);
void clp_rep(int m,int r,int a,int b);

int clu(void);
double clu_ld(double *d);
double clu_pd(double *d);
int clu_fn(void);
double clu_check(void);

int clpyr(void);
double clpyr_dis(int i,int j,float *d);
int clpyr_con(int i,int j,int *ack,int *acs,int *acl,int *acr);
int clpyr_ord(int i,int j,int n,int *idx_left,int *idx_right,int *idx_min,int *idx_max);
int clpyr_left(int i,int j,int *acl);
int clpyr_right(int i,int j,int *idx_right);
int clpyr_mcheck(int i,int *idx_e1,int *idx_e2,char *idx_mark);
void clpyr_mark(int i,int *idx_e1,int *idx_e2,char *idx_mark,int n);
int clpyr_cont(int i,int j,int *idx_e1,int *idx_e2);
int clpyr_ccl(int i,int j,int n,int nc,int *idx_e1,int *idx_e2);
void clpyr_den(int n,int nc,int *iord,int *idx_e1,int *idx_e2,float *acy,float *tmp);


/* global parameters */

int CLN = 0.0;
double *CLDis;
double CLRho = 0.0;

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

int hcls(void)
{
    register int i,j;
    int err,nrec,nmax,n,nn,m;                 
    double d,t;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("SAHN clustering. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto HCLSFin;

    if (gdd_check(PMGN,1))
        goto HCLSFin;
    if (gdd_tcheck(0,1,2))          /* need undirected valued graph */
        goto HCLSFin;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    if (PMNP == 1)  
        nmax = GD_NDMAX;
    else
        nmax = 0;

    if (PMOPT > 7)
        PMOPT = 7;
    printf1("\nOption %d: ",PMOPT);
    if (PMOPT == 1)
        printf1("single-link.\n");
    else if (PMOPT == 2)
        printf1("complete-link.\n");
    else if (PMOPT == 3)
        printf1("weighted average.\n");
    else if (PMOPT == 4)
        printf1("weighted centroid.\n");
    else if (PMOPT == 5)
        printf1("group average.\n");
    else if (PMOPT == 6)
        printf1("unweighted centroid.\n");
    else                 
        printf1("Ward's minimum variance.\n");
    newline();

    n  = GD_NP;
    nn = n * (n - 1) / 2;

    if (alloc_acn(n + 1))       /* a */   
        goto HCLSFin;
    if (alloc_acm(n + 1))       /* b */   
        goto HCLSFin;
    if (alloc_aci(n + 1))       /* p */   
        goto HCLSFin;
    if (alloc_acj(n + 1))       /* q */   
        goto HCLSFin;
    if (alloc_acyf(n + 1))      /* h */   
        goto HCLSFin;
    if (alloc_acxf(nn + 1))     /* d */   
        goto HCLSFin;

    if (h_cls(PMOPT,PMGN,n,AcXF,AcN,AcM,AcI,AcJ,AcYF)) {
        printf1("Error: found negative or missing edge value.\n");
        goto HCLSFin;
    }

    /* check for monotone level function */

    m = 0;
    for (i = 2; i < n; ++i) {
        if (AcYF[i] < AcYF[i-1]) {
            m = 1;
            break;
        }
    }
    printf1("Index function is ");
    if (m)
        printf1("NOT ");
    printf1("monotone.\n");

    /* create new distance matrix in AcUF */

    if (alloc_acr(2 * n + 2))                 
        goto HCLSFin;
    if (alloc_acuf(nn + 1))                 
        goto HCLSFin;

    h_cls_dm(n,AcN,AcM,AcYF,AcI,AcR,AcUF);

    d = h_cls_mu(n,AcUF); 

    printf1("Ultrametric condition: %lg\n",d);

    /* calculate fit of distance matrices */

    for (i = 1; i < n; ++i) {
        for (j = i + 1; j <= n; ++j) {
            t = (float)gdd_adj(i - 1,j - 1,PMGN);
            cl_putdu(t,i,j,AcXF);    
        }
    }
    d = 0.0;
    for (i = 1; i <= nn; ++i) {
        t = AcXF[i] - AcUF[i];
        d += t * t;
    }
    t = sqrt(2.0 * d);
    printf1("Fit of distances: %lg\n\n",t);


    nrec = 0;
    for (i = 1; i < n; ++i) {
        fprintf(PMFd,PMNFmtS,i);
        fprintf(PMFd,PMFmtS,(double)AcYF[i]);
        fprintf(PMFd,PMNFmtS,AcN[i]);
        fprintf(PMFd,PMNFmtS,AcM[i]);
        if (PMNP == 1) {
            fprintf(PMFd,PMNFmtS,gdd_node(AcN[i] - 1));
            fprintf(PMFd,PMNFmtS,gdd_node(AcM[i] - 1));
        }
        fprintf(PMFd,"\n");
        nrec++;
    }
    printf1("%d records written to: %s\n",nrec,PMFdName);

    if (PMTabFDef) {        /* write the distance matrix */
        for (i = 1; i <= n; ++i) {
            for (j = 1; j<= n; ++j) {
                if (j > i)
                    d = (double)cl_getdu(i,j,AcUF);
                else if (j < i)
                    d = (double)cl_getdu(j,i,AcUF);
                else 
                    d = 0.0;
                fprintf(PMTabFd,PMFmtS,d);
            }
            fprintf(PMTabFd,"\n");
        }
        printf1("%d records written to: %s\n",n,PMTabFName);
    }

    if (PMF1Def) {                       /* write dendrogram as edge list */
        if (alloc_aczf(n + 1))  
            goto HCLSFin;
        h_cls_el(nmax,n,AcN,AcM,AcYF,AcI,AcZF);
    }    
    if (PMPCFDef) {                      /* plot dendrogram */
        if (alloc_aczf(n + 1))  
            goto HCLSFin;
        if (alloc_acxf(n + 1))  
            goto HCLSFin;
        h_cls_den(GD_NP,AcN,AcM,AcYF,AcI,AcJ,AcZF,AcXF);
    }
    err = 0;

HCLSFin:       
    p_clean();
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

int h_cls(int ctyp,int gn,int n,float *d,int *a,int *b,int *p,int *q,float *h)
{
    register int i,j,k,l;
    int ic,jc,k1,k2,qi,qic,qjc;
    float dm,dmx,t,f,dj,dk;

    dm = 0.0;
    for (i = 1; i <= n; ++i) {
        p[i] = 0;
        q[i] = 1;
        for (j = i + 1; j <= n; ++j) {
            t = (float)gdd_adj(i - 1,j - 1,gn);
            if (t < 0.0)  
                return(-1);
            cl_putdu(t,i,j,d);    
            dm = dmax(t,dm);
        }
    }

    ic = jc = k = 0;
    while (++k < n) {

        dmx = dm + 10.0;
        for (i = 1; i < n; ++i) {
            if (p[i] == 0) {
                for (j = i + 1; j <= n; ++j) {
                    if (p[j] == 0) {
                        t = cl_getdu(i,j,d);
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
        h[k] = dmx;

        if (ctyp >= 5) {
            qic = q[ic];
            qjc = q[jc];

            if (ctyp != 7)
                f = 1.0 / (float)(qic + qjc);
        }

        for (i = 1; i <= n; ++i) {
            if (i != ic && p[i] == 0) {
                if (ctyp == 7) {
                    qi = q[i];
                    f = 1.0 / (float)(qi + qic + qjc);
                }
                j  = imin(ic,i);
                l  = imax(ic,i);
                k1 = imin(jc,i);
                k2 = imax(jc,i);
                dj = cl_getdu(j,l,d);
                dk = cl_getdu(k1,k2,d);

                switch (ctyp) {
                    case 1:     t = (float)dmin((double)dj,(double)dk);
                                break; 
                    case 2:     t = (float)dmax((double)dj,(double)dk);
                                break; 
                    case 3:     t = 0.5 * (dj + dk);
                                break; 
                    case 4:     t = 0.5 * (dj + dk) - 0.25 * dmx;
                                break; 
                    case 5:     t = f * ((float)qic * dj + (float)qjc * dk);
                                break; 
                    case 6:     t = f * ((float)qic * dj + (float)qjc * dk -
                                        (float)qic * (float)qjc * f * dmx);
                                break; 
                    case 7:     t = f * (((float)qic + (float)qi) * dj +
                                        ((float)qjc + (float)qi) * dk - 
                                         (float)qi * dmx);
                                break; 
                }
                cl_putdu(t,j,l,d);
                dm = dmax(dm,t);
            }
        }
        q[ic] += q[jc];

    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  cl_getdu(i,j,d)     return d(i,j) from upper triangle of d[].           */

float cl_getdu(int i,int j,float *d)
{
    return(d[(j - 1) * (j - 2) / 2 + i]);
}

/* ------------------------------------------------------------------------ */
/*  cl_putdu(x,i,j,d)   put x into d(i,j), upper triangle.                  */

void cl_putdu(float x,int i,int j,float *d)
{
    d[(j - 1) * (j - 2) / 2 + i] = x;
}

/* ------------------------------------------------------------------------ */
/*  h_cls_den   Create a command file for plotting a dendrogram.            */
/*              Assume data created by hcls().                              */

void h_cls_den(int n,int *a,int *b,float *h,int *p,int *q,float *x,float *y)
{
    register int i,j,k,l;
    int ic,jc;
    float max,d,t;

    for (i = 1; i <= n; ++i) {
        q[i] = p[i] = 0;
        x[i] = 0.0;
    }
    k = 0;
    for (i = n - 1; i >= 1; --i) {
        ic = imin(a[i],b[i]);
        jc = imax(a[i],b[i]);
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
    max = h[n - 1];
    d = max / 15.0;
    fprintf(PMPCFd,"psfile = %s.ps;\n",PMPCFName);
    fprintf(PMPCFd,"psetup(pxa=%g,%g,pya=0,%d);\n",-2.0 * d,max + d,n+1);
    for (i = 1; i <= n; ++i) {
        if (PMNP == 1)
            fprintf(PMPCFd,"pltext(xy=%g,%d)=%d;\n",-d,i,gdd_node(q[i] - 1));
        else
            fprintf(PMPCFd,"pltext(xy=%g,%d)=%d;\n",-d,i,q[i]);
    }
    for (i = 1; i <= n; ++i)  
        y[q[i]] = (float)i;

    for (i = 1; i < n; ++i) {
        ic = a[i];
        jc = b[i];
        t  = h[i];
        fprintf(PMPCFd,"plotp=%g,%g,%g,%g,%g,%g,%g,%g;\n",    
                          x[ic],y[ic],t,y[ic],t,y[jc],x[jc],y[jc]);
        x[ic] = x[jc] = t;
        y[ic] = y[jc] = (y[ic] + y[jc]) / 2.0;
    }
    printf1("Output file for plot of dendrogram: %s\n",PMPCFName);
}

/* ------------------------------------------------------------------------ */
/*  h_cls_el    Write dendrogram as an edge list to PMF1d.                  */
/*              Assume data created by hcls().                              */
/*                                                                          */
/*  If nmax = 0 begin counting the internal nodes with GD_NP + 1, otherwise */
/*  begin with nmax + 1.                                                    */

void h_cls_el(int nmax,int n,int *a,int *b,float *h,int *new,float *lev)
{
    register int i,ic,jc;
    int nrec,nn,n1,n2,m1,m2;
    float l1,l2;

    for (i = 1; i <= n; ++i) {
        new[i] = i;
        lev[i] = 0.0;
    }
    if (nmax == 0)
        nn = GD_NP + 1;
    else
        nn = nmax + 1;

    nrec = 0;
    for (i = 1; i < n; ++i) {
        ic = a[i];
        jc = b[i];
        if (new[ic] <= new[jc]) {
            n1 = new[ic];
            n2 = new[jc];
            l1 = h[i] - lev[ic];
            l2 = h[i] - lev[jc];
        }
        else {
            n1 = new[jc];
            n2 = new[ic];
            l1 = h[i] - lev[jc];
            l2 = h[i] - lev[ic];
        }
        if (nmax == 0)
            m1 = n1;
        else {
            if (n1 <= GD_NP)   
                m1 = gdd_node(n1 - 1);
            else
                m1 = n1;
        }
        if (nmax == 0)
            m2 = n2;
        else {
            if (n2 <= GD_NP)   
                m2 = gdd_node(n2 - 1);
            else
                m2 = n2;
        }
        fprintf(PMF1d,PMNFmtS,m1);
        fprintf(PMF1d,PMNFmtS,nn);
        fprintf(PMF1d,PMFmtS,l1);
        fprintf(PMF1d,"\n");
        nrec++;
        fprintf(PMF1d,PMNFmtS,m2);
        fprintf(PMF1d,PMNFmtS,nn);
        fprintf(PMF1d,PMFmtS,l2);
        fprintf(PMF1d,"\n");
        nrec++;
        new[ic] = new[jc] = nn;
        lev[ic] = lev[jc] = h[i];
        nn++;
    }
    printf1("%d records written to: %s\n",nrec,PMF1dName);
}

/* -##--------------------------------------------------------------------- */
/*  h_cls_dm    Create the distance matrix implied by the dendrogram        */
/*              in the array d[]; only upper triangle.                      */
/*                                                                          */

void h_cls_dm(int n,int *a,int *b,float *h,int *new,int *ll,float *d)
{
    register int i,j,k,ic,jc;
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
    k = 0;
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
            cl_putdu(h[ic - n],i,j,d);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  h_cls_mu    Return measure of ultrametric condition.                    */
/*              Assume upper triangle of distance matrix in d.              */
/*                                                                          */

double h_cls_mu(int n,float *d)
{
    register int i,j,k;               
    int ii,jj,kk;
    float t,dij,dik,djk;
         
    t = 0.0;
    for (i = 1; i < n; ++i) {
        ii = (i - 1) * (i - 2) / 2;
        for (j = i + 1; j <= n; ++j) {
            jj = (j - 1) * (j - 2) / 2;
            dij = d[jj + i];
            for (k = 1; k <= n; ++k) {
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

                if (dij <= djk && fabs(dik - djk) > EPSI1)
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

int nncl(void)
{
    register int i;
    int err,nrec,n;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Nearest-neighbor clustering. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto NNCLFin;

    if (gdd_check(PMGN,1))
        goto NNCLFin;
    if (gdd_tcheck(0,1,2))          /* need undirected valued graph */
        goto NNCLFin;

    if (PMALG != 2)
        PMALG = 1;

    printf1("Algorithm %d: ",PMALG);
    if (PMALG == 1) {
        if (PMSC <= 0.0)
            PMSC = 1.0;
        printf1("nearest neighbors with threshold %g\n",PMSC);
    }
    else {
        if (PMNS < 2)
            PMNS = 2;
        printf1("mutual neighborhoods of size %d\n",PMNS);
    }
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    nrec = 0;
    if (PMALG == 1) {
        if (alloc_acn(GD_NP))           /* cluster numbers */
            goto NNCLFin;

        n = nn_cl(AcN,PMSC);

        printf1("Found %d clusters.\n",n);
        for (i = 0; i < GD_NP; ++i) {
            fprintf(PMFd,PMNFmtS,i + 1);
            fprintf(PMFd,PMNFmtS,gdd_node(i));
            fprintf(PMFd,PMNFmtS,AcN[i]);
            fprintf(PMFd,"\n");
            nrec++;
        }
    }
    else {
        if (alloc_acn(GD_NP * PMNS + 1)) 
            goto NNCLFin;
        if (alloc_acr(GD_NP * PMNS + 1)) 
            goto NNCLFin;
        if (alloc_aci(GD_NP * PMNS + 1)) 
            goto NNCLFin;
        if (alloc_acj(GD_NP * PMNS + 1)) 
            goto NNCLFin;
        if (alloc_acx(GD_NP + 1)) 
            goto NNCLFin;

        nrec = mn_cl(PMNS,AcN,AcR,AcI,AcJ,AcX);
        if (nrec < 0)
            goto NNCLFin;
    }
    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

NNCLFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  nn_cl       Nearest-neighbor clustering.                                */
/*              Return: 0 if OK, -1 if error.                               */

int nn_cl(int *cn,double t)
{
    register int i,j,k;
    int m,nm;
    double d,dmin;

    nm = m = 1;                     /* number of clusters */
    cn[0] = 1;

    for (i = 1; i < GD_NP; ++i) {
        if (SILENTFlg < 2 && GD_NP >= 1000) {
            printfe("Node: %7d     %c",i + 1,CR);
            fflushe();
        }
        k = -1;
        dmin = 0.0;
        for (j = 0; j < i; ++j) {
            d = gdd_adj(i,j,PMGN);
            if (d >= 0.0 && d <= t) {
                dmin = d;
                k = cn[j];
            }
        }
        if (k < 0)  
            cn[i] = ++m;
        else   
            cn[i] = k;
    }
    if (SILENTFlg < 2 && GD_NP >= 1000) {
        printfe("\n");
        fflushe();
    }
    return(m);
}

/* ------------------------------------------------------------------------ */
/*  mn_cl       Mutual neighborhood clustering.                             */
/*              Return number of records written to output file, or -1 if   */
/*              error (= insufficient memory).                              */  

int mn_cl(int k,int *idx,int *mptr,int *iptr,int *jptr,double *val)
{
    register int i,j,l;
    int m,n,ip,jp,n1,ip1,jp1,nrec;

    for (i = 0; i < GD_NP; ++i) {
        if (SILENTFlg < 2 && GD_NP >= 1000) {
            printfe("Node: %7d     %c",i + 1,CR);
            fflushe();
        }
        n = mn_cl1(i,k,mptr,val,iptr,jptr); 
        if (n < 0)
            return(-1);

        for (j = 0; j < k; ++j) {
            if (j < n)
                idx[i * k + j] = mptr[j];
            else
                idx[i * k + j] = INTMAX;
        }
    }
    if (SILENTFlg < 2 && GD_NP >= 1000) {
        printfe("\n");
        fflushe();
    }
    m = 0;
    for (i = 0; i < GD_NP; ++i) {
        for (j = 0; j < k; ++j) {
            ip = idx[i * k + j];
            if (ip >= GD_NP)
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
                iptr[m] = imin(i,ip);
                jptr[m++] = imax(i,ip);
            }
        }
    }
    if (sortdpi2(m,mptr,iptr,idx))
        return(-1);

    nrec = 0;
    n1 = ip1 = jp1 = -1;
    for (i = 0; i < m; ++i) {
        j = idx[i];
        n  = mptr[j];
        ip = iptr[j];
        jp = jptr[j];
        if (n != n1 || ip != ip1 || jp != jp1) {
            fprintf(PMFd,PMNFmtS,nrec + 1);
            fprintf(PMFd,PMNFmtS,n);
            fprintf(PMFd,PMNFmtS,ip + 1);
            fprintf(PMFd,PMNFmtS,jp + 1);
            fprintf(PMFd,PMNFmtS,gdd_node(ip));
            fprintf(PMFd,PMNFmtS,gdd_node(jp));
            fprintf(PMFd,"\n");
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

int mn_cl1(int i,int k,int *nptr,double *val,int *sptr,int *idx)
{
    register int j;
    int n;
    double tmp;

    n = 0;
    for (j = 0; j < GD_NP; ++j) {
        if (j != i) {
            tmp = gdd_adj(i,j,PMGN);
            if (tmp >= 0.0) {
                idx[n] = j;
                val[n++] = tmp;
            }
        }
    }
    if (n < 1)
        return(0);

    if (sortdp(n,val,sptr))
        return(-1);

    n = imin(n,k);
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

int hcld(void)
{
    int err,r,cmax;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Hierarchical divisive clustering. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto HCLDFin;

    if (gdd_check(PMGN,1))
        goto HCLDFin;
    if (gdd_tcheck(0,1,2))          /* need undirected valued graph */
        goto HCLDFin;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);
    if (PMMin < 2)
        PMMin = 2;
    if (PMMax < 1)
        PMMax = 1;

    if (PMOPT > 2)
        PMOPT = 2;

    if (PMALG != 2)
        PMALG = 1;

    printf1("Algorithm %d: ",PMALG);
    if (PMALG == 1)
        printf1("maximal distance splits.\n");
    else                   
        printf1("minimal diameter splits.\n");

    if (alloc_acn(GD_NP + 1))      
        goto HCLDFin;
    if (alloc_acm(GD_NP + 1))      
        goto HCLDFin;
    if (alloc_aci(GD_NP + 1))      
        goto HCLDFin;

    if (PMALG == 2) {
        cmax = 2 * PMMax + 1;           /* max number of clusters */

        if (alloc_acu(cmax + 2))      
            goto HCLDFin;
        if (alloc_acj(GD_NP + 1))      
            goto HCLDFin;
    }

    if (PMALG == 1)
        r = h_cld1(GD_NP,PMMin,AcN,AcM,AcI,PMOPT);
    else
        r = h_cld2(GD_NP,PMMax,AcN,AcU,AcM,AcI,AcJ,PMOPT);

    if (r < 0)
        goto HCLDFin;

    err = 0;

HCLDFin:       
    p_clean();
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

int h_cld1(int n,int min,int *cidx,int *ncn,int *ncnt,int opt)
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

                fprintf(PMFd,PMNFmtS,kn);
                fprintf(PMFd,PMNFmtS,ncn[k]);
                for (l = 0; l < n; ++l) {
                    if (cidx[l] == kn)
                        fprintf(PMFd,PMNFmtS,l + 1);
                }
                fprintf(PMFd,"\n");
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
                    d = gdd_adj(i,j,PMGN);
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
                di = gdd_adj(i0,i,PMGN);
                dj = gdd_adj(j0,i,PMGN);
                if (di <= dj) {
                    cidx[i] = nch + nnc;
                    ncnt[nnc] += 1;
                }
                else {
                    cidx[i] = nch + nnc + 1;
                    ncnt[nnc + 1] += 1;
                }
            }
            if (PMF1Def) {   /* write split tree as edge list */

                fprintf(PMF1d,PMNFmtS,nch + nnc);
                fprintf(PMF1d,PMNFmtS,kn);
                fprintf(PMF1d,PMNFmtS,ncnt[nnc]);
                fprintf(PMF1d,PMNFmtS,ncn[k]);
                fprintf(PMF1d,PMFmtS,dmax / 2.0);
                fprintf(PMF1d,"\n");
                fprintf(PMF1d,PMNFmtS,nch + nnc + 1);
                fprintf(PMF1d,PMNFmtS,kn);
                fprintf(PMF1d,PMNFmtS,ncnt[nnc + 1]);
                fprintf(PMF1d,PMNFmtS,ncn[k]);
                fprintf(PMF1d,PMFmtS,dmax / 2.0);
                fprintf(PMF1d,"\n");
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
    printf1("%d records written to: %s\n",nr,PMFdName);
    if (nr1 > 0)
        printf1("%d records written to: %s\n",nr1,PMF1dName);
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

int h_cld2(int n,int nsmax,int *cidx,double *dia,int *nodes,             
    int *in,int *jn,int opt)
{
    register int i,j,k;
    int err,nr,nr1,nc,cn,ne,m,ns,ns1,ns2,n1,n2;
    double tmp;
       
    nr = nr1 = 0;
    err = -1;

    for (i = 0; i < n; ++i)
        cidx[i] = 1;

    nc = 1;
    dia[1] = h_cld_dia(n,cidx,1);

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

            ne = g_mst(PMGN,m,nodes,in,jn,1);
            if (ne <= 0) {
                if (ne == 0)  
                    printf1("ERROR: hcld2.\n");
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
                printf1("ERROR1: hcld2.\n");
                goto HCLD2Fin;
            }
        }

        /* find diameters of new clusters */

        dia[2 * ns] = h_cld_dia(n,cidx,2 * ns);
        dia[2 * ns + 1] = h_cld_dia(n,cidx,2 * ns + 1);

        /* print */

        if (opt == 1) {
            h_cld_prn(n,cidx,cn,ns1,n1,dia);
            h_cld_prn(n,cidx,cn,ns2,n2,dia);
            nr += 2;
        }
        if (PMF1Def) {   /* write split tree as edge list */

            fprintf(PMF1d,PMNFmtS,ns1);
            fprintf(PMF1d,PMNFmtS,cn);
            fprintf(PMF1d,PMNFmtS,n1);
            fprintf(PMF1d,PMNFmtS,m);
            fprintf(PMF1d,PMFmtS,dia[ns1]);
            fprintf(PMF1d,PMFmtS,dia[cn]);
            fprintf(PMF1d,"\n");

            fprintf(PMF1d,PMNFmtS,ns2);
            fprintf(PMF1d,PMNFmtS,cn);
            fprintf(PMF1d,PMNFmtS,n2);
            fprintf(PMF1d,PMNFmtS,m);
            fprintf(PMF1d,PMFmtS,dia[ns2]);
            fprintf(PMF1d,PMFmtS,dia[cn]);
            fprintf(PMF1d,"\n");
            nr1 += 2;
        }
        dia[cn] = -1.0;
        nc += 2;
    }
    printf1("%d records written to: %s\n",nr,PMFdName);
    if (nr1 > 0)
        printf1("%d records written to: %s\n",nr1,PMF1dName);

    err = 0;
         
HCLD2Fin:
    return(err);
}
   
/* ------------------------------------------------------------------------ */
/*  h_cld_prn()             print cluster.                                  */

void h_cld_prn(int n,int *cidx,int cn,int k,int nn,double *dia)
{
    register int i;

    fprintf(PMFd,PMNFmtS,k);
    fprintf(PMFd,PMNFmtS,cn);
    fprintf(PMFd,PMFmtS,dia[k]);
    fprintf(PMFd,PMNFmtS,nn);
    for (i = 0; i < n; ++i) {
        if (cidx[i] == k)
            fprintf(PMFd,PMNFmtS,i + 1);
    }
    fprintf(PMFd,"\n");
}

/* ------------------------------------------------------------------------ */
/*  h_cld_dia(n,cidx,k)     Return diameter of subgraph with nodes in       */
/*                          cidx[i] = k.                                    */

double h_cld_dia(int n,int *cidx,int k)
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

            tmp = gdd_adj(i,j,PMGN);
            if (tmp >= 0.0)  
                dia = dmax(dia,tmp);
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

int becl(void)
{
    register int i,j;
    int err,nrec,i0,cp,isel;
    double tmp,dmax;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("BE clustering. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto BECLFin;

    if (gdd_check(PMGN,1))
        goto BECLFin;
    if (gdd_tcheck(0,0,0))
        goto BECLFin;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);
    if (PMSCFlg == 0)
        PMSC = -1.0;
    if (PMOPT < 1 || PMOPT > 5)
        PMOPT = 1;

    if (PMMin == 1) {
        printf1("Minimize with ");
        dmax = DBLMAX;
    }
    else {
        PMMin = 0;
        printf1("Maximize with ");
        dmax = -1.0;
    }
    if (PMALG < 1 || PMALG > 3)
        PMALG = 1;

    if (PMALG == 1) 
        printf1("column");
    else if (PMALG == 2) 
        printf1("row");
    else              
        printf1("row and column");
    printf1(" permutations.\n");

    if (alloc_acn(GD_NP))      
        goto BECLFin;
    if (alloc_aci(GD_NP))      
        goto BECLFin;
    if (alloc_acj(GD_NP))      
        goto BECLFin;
    if (alloc_acs(GD_NP))      
        goto BECLFin;
   
    printf1("\nStarting node  final BE criterion\n");

    isel = -1;
    i0 = -1;
    cp = 0;
    while (cp < PMNCN) {
        j = PMCN[cp++] - 1;
        if (j >= 0 && j < GD_NP) {
            i0 = j;
            break;
        }
    }
    if (i0 < 0)
        i0 = 0;

    while (1) {
        if (PMALG == 1 || PMALG == 3) {
            tmp = be_cl1(1,PMMin,i0,GD_NP,AcN,AcI,AcJ);
            if ((PMMin == 0 && dmax < tmp) || (PMMin && dmax > tmp)) {
                dmax = tmp;
                isel = i0;
                for (i = 0; i < GD_NP; ++i)
                    AcS[i] = AcI[i];
            }
        }
        if (PMALG == 2 || PMALG == 3) {
            tmp = be_cl1(2,PMMin,i0,GD_NP,AcN,AcI,AcJ);
            if ((PMMin == 0 && dmax < tmp) || (PMMin && dmax > tmp)) {
                dmax = tmp;
                isel = i0;
                for (i = 0; i < GD_NP; ++i)
                    AcS[i] = AcI[i];
            }
        }
        printf1("%11d  %18.4f\n",i0 + 1,dmax);

        i0 = -1;
        while (cp < PMNCN) {
            j = PMCN[cp++] - 1;
            if (j >= 0 && j < GD_NP) {
                i0 = j;
                break;
            }
        }
        if (i0 < 0)
            break;
    }
    printf1("\nSelected for output: starting node %d.\n",isel + 1);
    nrec = be_clprn(PMOPT,GD_NP,AcS,PMSC);
    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

BECLFin:       
    p_clean();
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

int be_clprn(int opt,int n,int *idx,double mval)
{
    register int i,j;
    int ip,jp,m = 0;
    double tmp;

    if (opt == 1) {
        for (i = 0; i < n; ++i) {
            fprintf(PMFd,PMNFmtS,i + 1);
            fprintf(PMFd,PMNFmtS,idx[i] + 1);
            fprintf(PMFd,PMNFmtS,gdd_node(i));
            fprintf(PMFd,PMNFmtS,gdd_node(idx[i]));
            fprintf(PMFd,"\n");
            m++;
        }
    }
    else if (opt == 2) {
        if (GD_GT <= 2) {               /* undirected */
            for (i = 0; i < n; ++i) {
                ip = idx[i];
                for (j = i; j < n; ++j) {
                    jp = idx[j];
                    tmp = gdd_adj(ip,jp,PMGN);
                    if (tmp >= 0.0) {
                        fprintf(PMFd,PMNFmtS,i + 1);
                        fprintf(PMFd,PMNFmtS,j + 1);
                        fprintf(PMFd,PMNFmtS,ip + 1);
                        fprintf(PMFd,PMNFmtS,jp + 1);
                        fprintf(PMFd,PMNFmtS,gdd_node(ip));
                        fprintf(PMFd,PMNFmtS,gdd_node(jp));
                        fprintf(PMFd,PMFmtS,tmp);
                        fprintf(PMFd,"\n");
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
                    tmp = gdd_adj(ip,jp,PMGN);
                    if (tmp >= 0.0) {
                        fprintf(PMFd,PMNFmtS,i + 1);
                        fprintf(PMFd,PMNFmtS,j + 1);
                        fprintf(PMFd,PMNFmtS,ip + 1);
                        fprintf(PMFd,PMNFmtS,jp + 1);
                        fprintf(PMFd,PMNFmtS,gdd_node(ip));
                        fprintf(PMFd,PMNFmtS,gdd_node(jp));
                        fprintf(PMFd,PMFmtS,tmp);
                        fprintf(PMFd,"\n");
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
                tmp = gdd_adj(ip,jp,PMGN);
                if (tmp < 0.0)
                    tmp = mval;

                fprintf(PMFd,PMNFmtS,i + 1);
                fprintf(PMFd,PMNFmtS,j + 1);
                fprintf(PMFd,PMNFmtS,ip + 1);
                fprintf(PMFd,PMNFmtS,jp + 1);
                fprintf(PMFd,PMNFmtS,gdd_node(ip));
                fprintf(PMFd,PMNFmtS,gdd_node(jp));
                fprintf(PMFd,PMFmtS,tmp);
                fprintf(PMFd,"\n");
                m++;
            }
        }
    }   
    else {                      /* square matrix */

        for (i = 0; i < n; ++i) {
            ip = idx[i];

            if (opt == 5) {
                fprintf(PMFd,PMNFmtS,i + 1);
                fprintf(PMFd,PMNFmtS,ip + 1);
                fprintf(PMFd,PMNFmtS,gdd_node(ip));
            }
            for (j = 0; j < n; ++j) {
                jp = idx[j];
                tmp = gdd_adj(ip,jp,PMGN);
                if (tmp < 0.0)
                    tmp = mval;
                fprintf(PMFd,PMFmtS,tmp);
            }
            fprintf(PMFd,"\n");
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

double be_cl1(int opt,int mflag,int i0,int n,int *pos,int *idx,int *tidx)
{
    register int i,j,k;
    int nidx,knum,kdir,kidx,ktdir,ktidx,pos0,nt,fin,nn;
    double bmax,be;

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
            bmax = DBLMAX;

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
                    be = be_ccheck(nt,tidx);
                else
                    be = be_rcheck(nt,tidx);

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
    printf1("IDX: ");
    for (i = 0; i < n; ++i)
        printf1("%d ",idx[i]);
    newline();
    **/
    return(bmax);
}

/* ------------------------------------------------------------------------ */
/*  be_ccheck(n,idx)                                                        */
/*                                                                          */
/*  Calculate BE measure for columns with nodes in idx[i], i=0,...,n-1.     */

double be_ccheck(int n,int *idx)                          
{
    register int i,j,k;
    double be,tmp,tmp1;

    be = 0.0;
    for (i = 0; i < GD_NP; ++i) {
        for (j = 0; j < n; ++j) {
            k = idx[j];
            tmp = gdd_adj(i,k,PMGN);
            if (tmp <= 0.0)
                continue;

            if (j > 0 && (tmp1 = gdd_adj(i,idx[j - 1],PMGN)) > 0.0)
                be += tmp * tmp1;
            if (j < n - 1 && (tmp1 = gdd_adj(i,idx[j + 1],PMGN)) > 0.0)
                be += tmp * tmp1;
        }
    }
    return(be);
}

/* ------------------------------------------------------------------------ */
/*  be_rcheck(n,idx)                                                        */
/*                                                                          */
/*  Calculate BE measure for rows with nodes in idx[i], i=0,...,n-1.        */

double be_rcheck(int n,int *idx)                          
{
    register int i,j,k;
    double be,tmp,tmp1;

    be = 0.0;
    for (i = 0; i < GD_NP; ++i) {
        for (j = 0; j < n; ++j) {
            k = idx[j];
            tmp = gdd_adj(i,k,PMGN);
            if (tmp <= 0.0)
                continue;

            if (j > 0 && (tmp1 = gdd_adj(idx[j - 1],i,PMGN)) > 0.0)
                be += tmp * tmp1;
            if (j < n - 1 && (tmp1 = gdd_adj(idx[j + 1],i,PMGN)) > 0.0)
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

int acl(void)
{
    register int i,j,k; 
    int err,np,nk,nc,ncmax,nn,first,kk,pflag;
    double r;
    float c,mev;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Clustering with assignment. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto ACLFin;

    if (gdd_check(PMGN,1))
        goto ACLFin;
    if (gdd_tcheck(0,1,2))          /* need undirected valued graph */
        goto ACLFin;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(8,2);

    np = GD_NP;                     /* number of objects */

    nk = PMNCN;
    if (nk == 0) {
        printf1("Need value for cn parameter.\n");
        goto ACLFin;
    }
    ncmax = 0;
    for (k = 0; k < nk; ++k) {
        nc = PMCN[k];
        if (nc < 1 || nc > np / 2) {
            printf1("\nError: range of cn parameter is %d to %d\n",1,np/2);
            goto ACLFin;
        }
        ncmax = imax(ncmax,nc);
    }
    newline();

    if (alloc_acxf(np * np + 1))    /* distance matrix */
        goto ACLFin;
    if (alloc_acm(ncmax))           /* prototype numbers */
        goto ACLFin;
    if (alloc_acn(nk * ncmax))      /* final prototype numbers */
        goto ACLFin;
    if (alloc_acr(np))              /* cluster membership indicators */
        goto ACLFin;
    if (alloc_acns(np * nk))        /* final cluster membership indicators */
        goto ACLFin;
    if (alloc_acyf(nk))             /* final assignment value */
        goto ACLFin;

    mev = GD_EVMax[PMGN-1];

    pflag = 1;                      /* set to 0 if zero distances */
    k = 0;
    for (i = 0; i < np; ++i) {
        for (j = 0; j < np; ++j) {
            if (j == i)
                AcXF[k++] = 0.0;
            else {
                AcXF[k] = gdd_adj(i,j,PMGN);
                if (AcXF[k] == 0.0)
                    pflag = 0;
                k++;
            }
        }
    }
    /*******
    for (i = 0; i < np; ++i) {
        for (j = 0; j < np; ++j) {
            printf("%2.0f ",AcXF[i * np + j]);
        }
        newline();
    }
    ***/

    for (kk = 0; kk < nk; ++kk) {       /* for all cluster sizes */
        nc = PMCN[kk];
        AcYF[kk] = np * mev;

        first = 1;
        while (comb_nm(np,nc,AcM,first)) {
            for (i = 0; i < np; ++i)
                AcR[i] = -1;
            for (j = 0; j < nc; ++j)
                AcR[AcM[j]] = j;

            c = acl1(np,nc,AcXF,AcM,AcR,mev,AcYF[kk]);

            if (first || c < AcYF[kk]) {
                AcYF[kk] = c;
                j = kk * ncmax;
                for (k = 0; k < nc; ++k)  
                    AcN[j++] = AcM[k];
                j = kk * np;
                for (k = 0; k < np; ++k)
                    AcNS[j++] = AcR[k];
            }
            first = 0;
        }
    }
    for (kk = 0; kk < nk; ++kk) {
        nc = PMCN[kk];
        printf1("nc=%3d  cmin=%12.2f  prototypes: ",nc,AcYF[kk]);
        j = kk * ncmax;
        for (k = 0; k < nc; ++k)
            printf1("%3d ",AcN[j++] + 1);
        newline();
    }
    newline();

    /* print to output file */

    for (i = 0; i < np; ++i) {
        fprintf(PMFd,PMNFmtS,i + 1);
        for (k = 0; k < nk; ++k)  
            fprintf(PMFd,PMNFmtS,AcNS[k * np + i] + 1);
        fprintf(PMFd,"\n");
    }
    printf1("%d records written to: %s\n",np,PMFdName);

    if (PMF1Def) {          /* additional output file */
        if (pflag == 0) {
            printf1("df option ignored.\n");
            goto ACLFin;
        }
        if (alloc_actmp(ncmax))             
            goto ACLFin;

        nn = 0;
        for (kk = 0; kk < nk; ++kk) {
            nc = PMCN[kk];
            fprintf(PMF1d,"Assignment weights for nc = %d [%d",nc,AcN[kk * ncmax] + 1);
            for (k = 1; k < nc; ++k)
                fprintf(PMF1d,",%d",AcN[kk * ncmax + k] + 1);
            fprintf(PMF1d,"]\n");
            nn++;    

            for (i = 0; i < np; ++i) {
                r = 0.0;    
                for (k = 0; k < nc; ++k) {
                    j = AcN[kk * ncmax + k];
                    if (i == j) {
                        for (j = 0; j < nc; ++j)
                            AcTmp[j] = 0.0;
                        AcTmp[k] = 1.0;
                        r = 0.0;    
                        break;
                    }
                    AcTmp[k] = 1.0 / AcXF[i * np + j];
                    r += AcTmp[k]; 
                }
                if (r > 0.0) {
                    for (k = 0; k < nc; ++k)
                        AcTmp[k] /= r;
                }
                fprintf(PMF1d,PMNFmtS,i + 1);
                for (k = 0; k < nc; ++k)
                    fprintf(PMF1d,PMFmtS,AcTmp[k]);
                fprintf(PMF1d,"\n");
                nn++;
            }
            fprintf(PMF1d,"\n");
            nn++;
        }
        printf1("%d records written to: %s\n",nn,PMF1dName);
    }
    err = 0;

ACLFin:       
    p_clean();
    return(err);
}


/* ------------------------------------------------------------------------ */
/*  acl1(n,m,d,p,q,mev)                                                     */
/*                                                                          */
/*  n is number of objects, m number of clusters, d[] distance matrix,      */
/*  p[] indices of prototypes. q[] used to index cluster membership.        */
/*  mev is max edge value. Set indices in q[] and return min. assignment    */
/*  value.                                                                  */

float acl1(int n,int m,float *d,int *p,int *q,float mev,float cmin)
{
    register int i,j,k;
    float dm,t,c;
         
    c = 0.0;
    for (i = 0; i < n; ++i) {
        if (q[i] >= 0)
            continue;

        dm = mev;
        for (j = 0; j < m; ++j) {
            t = d[i * n + p[j]];
            if (t < dm) {
                dm = t;
                k = j;
            }
        }
        q[i] = k;
        c += dm;
        if (c > cmin)  
            break;
    }
    return(c);
}

/* ------------------------------------------------------------------------ */
/*  udcl(n,nc,d,p,l,u,c,pc,qa,qb,fp)                                        */
/*                                                                          */
/*  Find best partition into nc groups. Dynamic programming algorithm       */
/*  as described by C.J. Alpert, A.B. Kahng, Splitting an Ordering into a   */
/*  Partition to Minimize Diameter. J. of Classification 14 (1997), p. 64.  */
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
/*  pc[]    temporary cost matrix (i=0,nc-1,j=0,n-1)                        */
/*  qa[]    pointer for partitions (i=0,nc-1,j=0,n-1)                       */
/*  qb[]    pointer for partitions (i=0,nc-1,j=0,n-1)                       */
/*  fp[]    final partition (i=0,nc-1) [starting indices]                   */
/*                                                                          */
/*  The final cost is returned in pc[nc-1,n-1]                              */
/*                                                                          */

void udcl(int n,int nc,float *d,int *p,int l,int u,float *c,double *pc,
    int *qa,int *qb,int *fp)
{
    register int j,k,m;
    int km1,mp1,a,b;
    double fb,f;

    udcl1(n,d,p,c,u);           /* initial cost */
 
/********
printf("c matrix\n");
for (k = 0; k < n; ++k) {
    for (j=0; j < n; ++j)
        printf("%g ",c[k * n + j]);
    newline();
}
**/
    m = n * nc;
    for (j = 0; j < m; ++j)
        pc[j] = 0.0;

    for (j = 0; j < n; ++j)
        pc[j] = (double)c[j];

/**
printf("pc matrix\n");
for (k = 0; k < nc; ++k) {
    for (j=0; j < n; ++j)
        printf("%g ",pc[k * n + j]);
    newline();
}
**/




    for (k = 1; k < nc; ++k) {
        km1 = k - 1;

/* printf("\nk=%d km1=%d   l=%d u=%d    \n",k,km1,l,u );    */

        for (j = 0; j < n; ++j) {
            fb = DBLMAX;
/*  printf("jjjjj=%d fb=%g\n",j,fb);       */

            for (m = l - 1; m < u; ++m) {
/*  printf("m=%d l-1=%d u=%d\n",m,l-1,u);       */
                    if (m < j - u)  
                        continue;
                    if (m > j - l)  
                        break;
                      
                mp1 = m + 1;
/**
printf("pc[km1*n+m]=%g km1=%d m=%d\n",pc[km1*n+m],km1,m  );
printf(" c[mp1*n+j]=%g mp1=%d j=%d\n", c[mp1*n+j],mp1,j  );
**/
                f = dmax(pc[km1 * n + m],c[mp1 * n + j]);

/*  printf("fb=%g  f==%g\n\n",fb, f);     */

           
            
                if (fb > f) {
                    fb = f;

/*  printf("put f=%g into pc[k*n+j] k=%d j=%d\n\n",f,k,j); */


                    pc[k * n + j] = f;               
                    qa[k * n + j] = m;
                    qb[k * n + j] = j;
                }
            } 
        }
    }

/***********
printf("pc matrix\n");
for (k = 0; k < nc; ++k) {
    for (j=0; j < n; ++j)
        printf("%g ",pc[k * n + j]);
    newline();
}
printf("qa matrix\n");
for (k = 0; k < nc; ++k) {
    for (j=0; j < n; ++j)
        printf("%d ",qa[k * n + j]);
    newline();
}
printf("qb matrix\n");
for (k = 0; k < nc; ++k) {
    for (j=0; j < n; ++j)
        printf("%d ",qb[k * n + j]);
    newline();
}
*********/

    m = n - 1;
    for (k = nc - 1; k >= 0; --k) {
        a = qa[k * n + m];
        b = qb[k * n + m];

        printf("k=%d  qa=%d  qb=%d m=%d    \n",k,a,b,m );          

        fp[k] = a + 1;

        m = a;
    }
    fp[0] = 0;
/**
printf("FP: ");
for (j = 0; j < nc; ++j)
    printf("%d ",fp[j]);
newline();

**/

}

/* ------------------------------------------------------------------------ */
/*  udcl1(n,d,p,c,u)            calculate cluster costs c[i,j]              */
/*                                                                          */
/*  n       number of data points in distance matrix d(i,j)                 */
/*  d[]     distance matrix (i,j = 0,...,n-1)                               */
/*  p[]     permutation of objects (i = 0,...,n-1)                          */
/*  c[]     cluster costs (i,j = 0,...,n-1)                                 */
/*  u       max cluster size                                                */

void udcl1(int n,float *d,int *p,float *c,int u)
{
    register int i,j,l;

    for (i = 0; i < n; ++i)
        c[i * n + i] = 0.0;

    for (l = 1; l < u; ++l) {
        for (i = 0; i < n; ++i) {
            j = (i + l) % n;
            c[i * n + j] =
            dmax(dmax((double)c[i * n + j - 1],(double)c[(i+1) * n + j]),
                      (double)d[p[i] * n + p[j]]);
        }
    }
}



/* ------------------------------------------------------------------------ */
/*  ucl         test                                                        */
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

int ucl(void)
{
    register int i,j,k; 
    int err,np,nc,pflag,l,u;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("ucl. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto UCLFin;

    if (gdd_check(PMGN,1))
        goto UCLFin;
    if (gdd_tcheck(0,1,2))          /* need undirected valued graph */
        goto UCLFin;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(8,2);

    np = GD_NP;                     /* number of objects */

    nc = 3;
    l = 1;
    u = 3;

    if (alloc_acxf(np * np + 1))    /* distance matrix */
        goto UCLFin;
    if (alloc_acyf(np * np + 1))    /* cost matrix */
        goto UCLFin;
    if (alloc_acu(nc * np + 1))     /* cost matrix for partitions */
        goto UCLFin;
    if (alloc_ack(np))              /* permutation */
        goto UCLFin;
    if (alloc_acn(nc * np + 1))     /* pointer to partitions */
        goto UCLFin;
    if (alloc_acm(nc * np + 1))     /* pointer to partitions */
        goto UCLFin;
    if (alloc_aci(nc + 1))          /* final partitions */
        goto UCLFin;


    for (i = 0; i < np; ++i)
        AcK[i] = i;

    AcK[0] = 1;
    AcK[1] = 0;

    pflag = 1;                      /* set to 0 if zero distances */
    k = 0;
    for (i = 0; i < np; ++i) {
        for (j = 0; j < np; ++j) {
            if (j == i)
                AcXF[k++] = 0.0;
            else {
                AcXF[k] = gdd_adj(i,j,PMGN);
                if (AcXF[k] == 0.0)
                    pflag = 0;
                k++;
            }
        }
    }
                    
    for (i = 0; i < np; ++i) {
        for (j = 0; j < np; ++j) {
            printf1("%2.0f ",AcXF[i * np + j]);
        }
        newline();
    }
              
               
    udcl(np,nc,AcXF,AcK,l,u,AcYF,AcU,AcN,AcM,AcI);


    printf1("Min objective function: %g\n",AcU[(nc - 1) * np + np - 1]);
    printf1("Final partition: ");
    for (j = 0; j < nc; ++j)
        printf1("%d ",AcI[j]);
    newline();

    err = 0;

UCLFin:       
    p_clean();
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

int hclsp(void)
{
    register int i; 
    int err,gn,n,np,ne,rn,nl,level,len;
    char buf[60];

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Partitioning the leafs of a tree. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 5,1,1))           /* get parameters */
        goto HCLSPFin;

    gn = 1;                             /* always graph number 1 */
    if (gdd_check(gn,0))        
        goto HCLSPFin;
    if (gdd_tcheck(1,2,0))              /* need undirected valued graph */
        goto HCLSPFin;

    if (PMFmtF == 0)                    /* default print format */
        pmfmt(8,2);

    np = GD_NP;                         /* number of nodes */
    ne = GD_NE[0];                      /* number of edges */
    printf1("Number of nodes: %d\n",np);
    printf1("Number of edges: %d\n",ne);
    if (ne != np - 1) {
        printf1("Error: graph is not a tree.\n");
        goto HCLSPFin;
    }
    nl = 0;
    rn = -1;
    for (i = 0; i < np; ++i) {
        if (GD_BPN[i] == 0)
            nl++;
        if (GD_FPN[i] == 0) {
            if (rn < 0)
                rn = i;
            else {
                rn = -1;
                break;
            }
        }
    }
    if (nl == 0 || rn < 0) {
        printf1("Error: graph is not a tree or has no single root.\n");
        goto HCLSPFin;
    }
    printf1("Number of leafs: %d\n",nl);
    printf1("Root node number: %d\n",GD_ND[rn]);

    len = 1;
    for (i = 0; i < np; ++i) {
        sprintf(buf,"(%d,%d) ",GD_ND[i],nl);
        len = imax(len,strlen(buf));
    }
    if (alloc_acn(np + 1))
        goto HCLSPFin;

    for (i = 0; i < np; ++i) {
        if (GD_BPN[i] == 0)
            hclsp_visit(i,gn);
    }
    if (PMNLEV > 0) {
        fprintf(PMFd,"Hierarchical structure of internal nodes (number of leafs in brackets).\n\n");
        level = 0;
        hclsp_prn(rn,gn,PMNLEV,level,len);
        fprintf(PMFd,"\n");   
    }
    if (PMNCN > 0) {
        if (alloc_acd(np + 1))
            goto HCLSPFin;

        for (i = 0; i < PMNCN; ++i) {
            n = gdd_fndni(PMCN[i]);
            if (n < 0) {
                printf1("Error: node number %d not defined.\n",PMCN[i]);
                goto HCLSPFin;
            }
            hclsp_prn1(n,gn,PMCN[i]);
        }
        n = 0; 
        for (i = 0; i < np; ++i) {
            if (AcD[i] && GD_BPN[i] == 0)
                n++;
        }
        printf1("Number of leafs in the output file: %d\n",n);
    }
    err = 0;

HCLSPFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  hclsp_visit(i,gn)                                                       */

void hclsp_visit(int i,int gn)
{
    register int j,l,k;
    int n;

    AcN[i] += 1;               
    n = GD_FPN[i];
    if (n <= 0)
        return;          

    for (l = 0; l < n; ++l) {
        j = GD_FPI[i][l];                    /* edge from k to j */
        k = GD_FPK[i][l];
        if (k >= 0 && gdd_ev(k,gn) >= 0.0)          
            hclsp_visit(j,gn);
    }
}

/* ------------------------------------------------------------------------ */
/*  hclsp_prn(i,gn,nlev,level,mlen)                                         */

void hclsp_prn(int i,int gn,int nlev,int level,int mlen)
{
    register int j,l,k,ll;
    int n,nline,len;
    char buf[60];

    if (level >= nlev)  
        return;
           
    if (level != 0)  
        fprintf(PMFd,"<- ");
    sprintf(buf,"(%d,%d) ",GD_ND[i],AcN[i]);
    len = strlen(buf);
    fprintf(PMFd,"%s",buf);
    fprnchar(PMFd,' ',mlen - len,0);

    if (level + 1 == nlev || (n = GD_BPN[i]) <= 0) {
        fprintf(PMFd,"\n");   
        return;
    }
    nline = 0;
    for (l = 0; l < n; ++l) {
        j = GD_BPI[i][l];                    /* edge from k to j */
        k = GD_BPK[i][l];
        if (k >= 0 && gdd_ev(k,gn) >= 0.0) {        
            if (nline) {
                fprnchar(PMFd,' ',mlen,0);
                for (ll = 0; ll < level; ++ll) 
                    fprnchar(PMFd,' ',mlen + 3,0);
            }
            hclsp_prn(j,gn,nlev,level + 1,mlen);
            nline = 1;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  hclsp_prn1(i,gn)                                                       */

void hclsp_prn1(int i,int gn,int m)
{
    register int j,l,k;
    int n;

    if (AcD[i])
        return;

    AcD[i] = 1;
    if (GD_BPN[i] == 0) {       /* print only for leafs */
        fprintf(PMFd,PMNFmtS,m);
        fprintf(PMFd,PMNFmtS,GD_ND[i]);
        fprintf(PMFd,"\n");
    }                
    if ((n = GD_BPN[i]) <= 0)
        return;
                               
    for (l = 0; l < n; ++l) {
        j = GD_BPI[i][l];                    /* edge from k to j */
        if (AcD[j] == 0) {
            k = GD_BPK[i][l];
            if (k >= 0 && gdd_ev(k,gn) >= 0.0)          
                hclsp_prn1(j,gn,m);
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

int scla(void)
{
    register int i,j;
    int err,fin,nnsp,nc,n;
    double x;


    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Find separable clusters. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto SCLAFin;

    if (gdd_check(PMGN,1))
        goto SCLAFin;
    if (gdd_tcheck(0,1,2))          /* need undirected valued graph */
        goto SCLAFin;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    printf1("Number of nodes: %d\n",GD_NP);

    /* For i=0,GD_NP-1: AcN[i] = 1,2,3,... cluster number, or 0 if node i
       cannot be used for an s-cluster. */

    if (alloc_acn(GD_NP))                 
        goto SCLAFin;
    if (alloc_acm(GD_NP))                 
        goto SCLAFin;
    if (alloc_acx(GD_NP))                 
        goto SCLAFin;

    for (i = 0; i < GD_NP; ++i)
        AcN[i] = -1;

    nnsp = 0;       /* number of point that cannot be userd for s-clusters */
    fin = 0;
    while (fin == 0) {
        fin = 1;
        for (i = 0; i < GD_NP; ++i) {
            if (AcN[i] >= 0)
                continue;
            if (scla_check(i,GD_NP,AcN,AcM,AcX)) {
                AcN[i] = 0;
                nnsp++;
                fin = 0;
                break;
            }
        }
    }
    printf1("\n%d isolated points",nnsp);
    if (nnsp == 0)
        printf1(".\n");
    else {
        j = 0;
        for (i = 0; i < GD_NP; ++i) {
            if (AcN[i] == 0) {
                if (j && j < nnsp)
                    printf1(",");
                else
                    printf1(": ");
                printf1("%d",i + 1);
                j++;
            }
        }
        printf1("\n");
    }
    n = 0;
    for (i = 0; i < GD_NP; ++i) {
        if (AcN[i] < 0)  
            AcM[n++] = i;
    }
    printf1("Remaining number of points: %d\n\n",n);
    if (n < 4) {
        err = 0;
        goto SCLAFin;
    }
    if (alloc_ack(GD_NP))                 
        goto SCLAFin;
    nc = scla_find(GD_NP,AcN,n,AcM,AcK);       
    err = 0;

SCLAFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  scla_check()                                                            */
/*                                                                          */
/*  Return  1 if i cannot be used for an s-cluster                          */
/*          0 otherwise                                                     */
        
int scla_check(int i,int n,int *acn,int *p,double *dp)
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
    dmin = DBLMAX;
    for (j = 0; j < n; ++j) {
        if (acn[j] == -1) {
            x = gdd_adj(i,j,PMGN);            
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
                    if (gdd_adj(p[k],j,PMGN) <= dp[k]) {
                        jx = j;
                        break;
                    }
                }
                if (jx >= 0) {
                    dm = 0.0;
                    for (k = 0; k < np; ++k) {
                        x = gdd_adj(p[k],jx,PMGN);
                        if (dp[k] < x)
                            dp[k] = x;
                        dm = dmax(dm,x);
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
        
int scla_find(int n,int *acn,int np,int *p,int *com)
{
    register int i,j,k,l;
    int first,nc,nu,nf;
    double dm;

    for (nc = 2; nc <= np - 2; ++nc) {
        printf1("Size %d\n",nc);
        nf = 0;
        first = 1;
        while (comb_nm(np,nc,com,first)) {
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
                        dm = dmax(dm,gdd_adj(j,p[com[k]],PMGN));            
                }
                for (l = 0; l < n; ++l) {
                    if (acn[l] < 0 && gdd_adj(l,j,PMGN) <= dm) {
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
                printf1("%4d -: ",++nf);
                for (i = 0; i < nc; ++i)
                    printf1("%d ",p[com[i]] + 1);
                newline();
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

int clp(void)
{
    register int i,j,k;
    int err,n,nn,r,kmin;    
    double t,cmin;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Partitional clustering. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto CLPFin;

    if (gdd_check(PMGN,1))
        goto CLPFin;
    if (gdd_tcheck(0,0,0))
        goto CLPFin;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    if (PMMETH < 1 || PMMETH > 4)
        PMMETH = 2;

    if (PMNC < 2)
        PMNC = 2;
    if (PMNC >= GD_NP) {
        printf1("Number of clusters must be less than %d\n",GD_NP);
        goto CLPFin;
    }        

    if (PMNS < 1)
        PMNS = 1;

    if (MxItFlg == 0)
        MxIter = 15;

    printf1("\nMethod: %d\n",PMMETH);

    if (PMMETH == 4) {
        printf1("Find %d points with maximal min-distance.\n",PMNC);
        err = clp_p();
        goto CLPFin;
    }

    printf1("Number of clusters: %d\n",PMNC);
    printf1("Number of random restarts: %d\n",PMNS);
    printf1("Maximal number of iterations: %d\n",MxIter);

    n = GD_NP;
    nn = n * (n - 1) / 2;

    if (alloc_act(nn + 1))      
        goto CLPFin;
    if (alloc_acu(PMNC + 1))      
        goto CLPFin;
    if (alloc_acv(PMNC + 1))      
        goto CLPFin;
    if (alloc_acn(PMNC + 1))      
        goto CLPFin;
    if (alloc_ack(n + 1))      
        goto CLPFin;
    if (alloc_acl(n + 1))      
        goto CLPFin;
   
    k = 0;
    for (i = 1; i < n; ++i) {
        for (j = 0; j < i; ++j)  
            AcT[++k] = gdd_adj(i,j,PMGN);
    }
    newline();
                
    cmin = DBLMAX;
    kmin = -1;

    for (k = 1; k <= PMNS; ++k) {           /* repetitions */
        for (i = 1; i <= n; ++i) {          /* create initial partition */
            t = random1();
            t *= (double)PMNC;
            j = (int)(t + 1.0);
            if (j < 1)
                j = 1;
            else if (j > PMNC)
                j = PMNC;
            AcK[i] = j;
        }
        r = clp_f(PMMETH,n,PMNC,MxIter,AcT,AcK,AcN,AcU,AcV);

        if (r >= 0) {
            t = 0.0;
            for (j = 1; j <= PMNC; ++j)
                t += AcU[j];
            if (t < cmin) {
                kmin = k;
                cmin = t;
                for (i = 1; i <= n; ++i)
                    AcL[i] = AcK[i];
            }
        }
        else
            t = -1.0;

        printf1("%4d %2d ",k,r);
        printf1(PMFmtS,t);
        for (j = 1; j <= PMNC; ++j) {
            printf1("%3d ",AcN[j]);
            printf1(PMFmtS,AcU[j]);
        }
        newline();
    }
    printf1("\nBest solution: %d\n",kmin);

    if (kmin > 0) {
        for (i = 1; i <= n; ++i)  
            fprintf(PMFd,"%4d %3d\n",i,AcL[i]);
        printf1("%d records written to: %s\n",n,PMFdName);
    }
    err = 0;

CLPFin:       
    p_clean();
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


int clp_f(int meth,int m,int n,int itmax,
          double *t,int *z,int *mj,double *e,double *c)
{
    register int i,j,h,l,k;
    int it,mo,p,q,im,ki,kh,is;
    double d,f,v,v1,u,u1,w,w1,ej,ep,eq,bj,cj;
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
    eq = DBLMAX;

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
  
double CLP_DMIN = 0.0;

int clp_p(void)
{
    register int i,j;
    int n,m;

    n = GD_NP;
    m = PMNC;
    if (alloc_acn(m))                     
        return(-1);          
    if (alloc_acm(m))                     
        return(-1);          

    for (i = 0; i < m; ++i)     /* current best config */
        AcM[i] = i;

    CLP_DMIN = DBLMAX;
    for (i = 1; i < PMNC; ++i) {
        for (j = 0; j < i; ++j)  
            CLP_DMIN = dmin(CLP_DMIN,gdd_adj(i,j,PMGN));
    }
    printf1("CLP_DMIN=%lg\n",CLP_DMIN);

    clp_rep(PMNC,PMNC,0,n - m);

    printf1("Best value: %lg\n",CLP_DMIN);
    printf1("Points: ");
    for (i = 0; i < m; ++i)
        printf1("%d ",AcM[i] + 1);
    newline();

    return(0);
}

void clp_rep(int m,int r,int a,int b)
{
    register int i,j,k;
    int cflag;
    double vmin;

    if (a > b)
        return;
    for (k = a; k <= b; ++k) {
           
        cflag = 0;
        for (j = 0; j < m - r; ++j) {
            if (gdd_adj(j,k,PMGN) < CLP_DMIN) {
                cflag = 1;
/*              printf(" BREAK bei r=%d k==%d a=%d b=%d\n",r, k,a,b  ); */
                break;
            }
        }
        if (cflag)
            continue;
          
        AcN[m - r] = k;

/*      printf(" r=%d k==%d a=%d b=%d\n",r, k,a,b  );   */
        if (r > 1)
            clp_rep(m,r-1,k+1,b+1);
        else {
/***************
            printf("ACN: ");
            for (j = 0; j < m; ++j)
                printf("%d ",AcN[j]);
            newline();

            printf("CURRENT dmin=%lg\n",CLP_DMIN);
*******/
            /* try to find better value */
                  
            vmin = DBLMAX;
            for (i = 1; i < m; ++i) {
                for (j = 0; j < i; ++j)  
                    vmin = dmin(vmin,gdd_adj(AcN[i],AcN[j],PMGN));
            }
/*          printf1("vmin=%lg\n",vmin);     */
            if (vmin > CLP_DMIN) {
                printf("NEW CLP_DMIN=%lg\n",vmin);           
                CLP_DMIN = vmin;

                for (i = 0; i < m; ++i)
                    AcM[i] = AcN[i];
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
       
int clu(void)
{
    register int i,j,k;
    int err,n,nn,it;
    double d,t,v;

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("LS fit of ultrametric distances). Current memory: %d bytes.\n",MemReq);

    MxIt1 = 10;

    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto CLUFin;

    if (gdd_check(PMGN,1))
        goto CLUFin;
    if (gdd_tcheck(0,1,2))          /* need undirected valued graph */
        goto CLUFin;

    if (GD_NP < 4) {
        printf1("Error: need at least four nodes.\n");
        goto CLUFin;
    }
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    CLN = n = GD_NP;                /* number of points */

    /* setup MDSDIS vector with distances */

    nn = n * (n - 1) / 2;
    if (alloc_acu(nn + 1))          /* distances */      
        goto CLUFin; 
    if (alloc_act(nn + 1))          /* previous parameter vector */      
        goto CLUFin; 

    CLDis = AcU;

    k = 0;
    for (i = 2; i <= n; ++i) {      /* get distances */
        for (j = 1; j < i; ++j)  
            CLDis[++k] = gdd_adj(i-1,j-1,PMGN);
    }

    printf("CLDis ");
    for (i = 1; i <= nn; ++i)
        printf("%lg\n",CLDis[i]);
    newline();

    /* prepare for min algorithm, always MINA = 7 */

    if (MINA == 5 || MINA == 6 || MINA == 8)
        MINA = 7;

    NParm = nn;                  
    PMDOPT = 1;                     /* analytical gradient */
    CCTyp = 0;                      /* no covariance matrix */

    set_mlopt();                    /* adjust options */         

    if (ml_init(1,1,0))             /* init function minimization */
        goto CLUFin;

    if (get_dsv(NParm,Par,0,ParLB,ParUB,1)) {    /* try to get starting values */
        printf1("Error in starting values.\n");
        goto CLUFin;
    }  

    if (DSVFlg == 0) {              /* create starting values */
        t = d = 0.0;
        for (i = 1; i <= nn; ++i) {
            d += CLDis[i];
            t += CLDis[i] * CLDis[i];
        }
        v = (t / (double)nn) - (d * d / (double)(nn * nn));
        printf1("Variance of distances: %lg\n",v);
        v = sqrt(v /= 3.0);

        for (i = 1; i <= nn; ++i)    
            Par[i] = CLDis[i] + v * normal();
    }
    for (i = 1; i <= nn; ++i)   
        printf("%10.4lf %10.4lf\n",CLDis[i],Par[i]);



    d = clu_ld(Par);
    t = clu_pd(Par);

    if (fabs(t) <= EPSI1) {
        printf1("\nPenalty term in starting values: %lg\n",t);
        err = 0;
        goto CLUFin;
    }
    CLRho = d / t;
    if (CLRho > 1.0)
        CLRho = 1.0;

    printf("ld=%lg pd=%lg rho=%lg\n",d,t,CLRho);
                  
               
    for (it = 1; it <= MxIt1; ++it) {   /* main loop */

        printf1("Iteration %d  CLRho=%lg\n",it,CLRho);

        for (i = 1; i <= nn; ++i)       /* save parameter vector */
            AcT[i] = Par[i];

        if (ffmin(CLMOD1,0,1,1))        /* minimization */
            goto CLUFin;                /* insuff memory or fn error */

        prn_mlres(1);                   /* print info about minimization */

        /* print result into output file */

        if (PMFDef) {
            fprintf(PMFd,"%4d %3d %12.8lf ",it,LConv,FMin);
            d = clu_ld(Par);
            fprintf(PMFd,"%12.8lf ",d);
            d = clu_pd(Par);
            fprintf(PMFd,"%12.8lf  ",d);
            for (i = 1; i <= NParm; ++i)
                fprintf(PMFd,PMFmtS,Par[i]);
            fprintf(PMFd,"\n");
        }

        /* check convergence */

        t = 0.0;
        for (i = 1; i <= nn; ++i) {
            d = Par[i] - AcT[i];
            t += d * d;
        }
        t = sqrt(t);
        printf1("Change in parameter vector: %lg\n",t);
        if (t <= PMTol)  
            break;
        CLRho *= 10.0;
    }
    for (i = 1; i <= nn; ++i) {
        printf("%10.4lf %10.4lf %10.4lf\n",CLDis[i],AcT[i],Par[i]);
    }

    d = clu_ld(Par);
    printf1("Final value of L(D): %lg\n",d);
    d = clu_pd(Par);
    printf1("Final value of P(D): %lg\n",d);
    d = clu_check();
    printf1("Ultrametric condition: %lg\n",d);

    if (PMF1Def) {      /* write distance matrix */
        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= n; ++j) {
                if (j < i) 
                    d = Par[(i - 1) * (i - 2) / 2 + j];
                else if (i < j) 
                    d = Par[(j - 1) * (j - 2) / 2 + i];
                else
                    d = 0.0;
                fprintf(PMF1d,PMFmtS,d);
            }
            fprintf(PMF1d,"\n");
        }
        printf1("%d records written to: %s\n",n,PMF1dName);
    }
    err = 0;

CLUFin:       
    ml_init(0,0,0);      
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  clu_ld()    Return L(D).                                                */

double clu_ld(double *d)
{
    register int i,j,k;               
    int err,ii,jj,kk;
    double t,f;
         
    f = 0.0;
    for (i = 1; i <= NParm; ++i) {
        t = CLDis[i] - d[i];
        f += t * t;
    }   
    return(f);
}

/* ------------------------------------------------------------------------ */
/*  clu_pd()    Return P(D).                                                */

double clu_pd(double *d)
{
    register int i,j,k;               
    int err,ii,jj,kk;
    double f,dij,dik,djk;

    f = 0.0;
    for (i = 2; i <= CLN; ++i) {
        ii = (i - 1) * (i - 2) / 2;
        for (j = 1; j < i; ++j) {
            jj = (j - 1) * (j - 2) / 2;
            dij = d[ii + j];
            for (k = 1; k <= CLN; ++k) {
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

                if (dij <= djk && fabs(dik - djk) > EPSI1)
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

int clu_fn(void)
{
    register int i,j,k,l;             
    int err,ii,jj,kk;
    double d,dij,dik,djk,fa,fb;
         
    fa = clu_ld(TPar);
    fb = clu_pd(TPar);
    if (LFunc)
        FTmp = fa + CLRho * fb;                 

    if (LGrad) {
        for (l = 1; l <= NParm; ++l)  
            Grad[l] = 2.0 * (TPar[l] - CLDis[l]);
           
        for (l = 1; l <= NParm; ++l) {
            for (i = 2; i <= CLN; ++i) {
                ii = (i - 1) * (i - 2) / 2;
                for (j = 1; j < i; ++j) {
                    jj = (j - 1) * (j - 2) / 2;
                    dij = TPar[ii + j];

                    for (k = 1; k <= CLN; ++k) {
                        if (k == i || k == j)
                            continue;

                        kk = (k - 1) * (k - 2) / 2;

                        if (k < i)
                            dik = TPar[ii + k];
                        else
                            dik = TPar[kk + i];

                        if (dij > dik)
                            continue;

                        if (k < j)
                            djk = TPar[jj + k];
                        else
                            djk = TPar[kk + j];

                        if (dij <= djk && fabs(dik - djk) > EPSI1) {

                            d = CLRho * 2.0 * (dik - djk);

                            if (k < i && l == (ii + k))               
                                Grad[l] += d;
                            else if (k > i && l == (kk + i))               
                                Grad[l] += d;
                            if (k < j && l == (jj + k))               
                                Grad[l] -= d;
                            else if (k > j && l == (kk + j))               
                                Grad[l] -= d;
                        }
                    }
                }
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
/*  clu_check()     Check for ultrametric condition.                        */

double clu_check(void)
{
    register int i,j,k;               
    int ii,jj,kk;
    double d,dij,dik,djk;
         
    d = 0.0;
    for (i = 2; i <= CLN; ++i) {
        ii = (i - 1) * (i - 2) / 2;
        for (j = 1; j < i; ++j) {
            jj = (j - 1) * (j - 2) / 2;
            dij = Par[ii + j];
            for (k = 1; k <= CLN; ++k) {
                if (k == i || k == j)
                    continue;

                kk = (k - 1) * (k - 2) / 2;

                if (k < i)
                    dik = Par[ii + k];
                else
                    dik = Par[kk + i];
                                           
                if (dij > dik)
                    continue;

                if (k < j)
                    djk = Par[jj + k];
                else
                    djk = Par[kk + j];

                if (dij <= djk && fabs(dik - djk) > EPSI1)
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

int clpyr(void)
{
    register int i,j,k,l;
    int err,n,nmax,nn,nc,r,isel,jsel,ip,sflag,il,ir;
    int *idx_min,*idx_max,*idx_left,*idx_right,*idx_e1,*idx_e2;
    char *idx_mark;
    double d,dm;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Pyramidal clustering. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 5,1,1))       /* get parameters */
        goto CLPYRFin;

    if (gdd_check(PMGN,1))
        goto CLPYRFin;
    if (gdd_tcheck(0,0,0))
        goto CLPYRFin;

    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    n = GD_NP;
    nmax = (n * (n + 1)) / 2 + 1;
    nn = nmax * (nmax - 1) / 2 + 2;

    if (alloc_acxf(nn + 1))         /* distances */   
        goto CLPYRFin;
    if (alloc_aci(nmax + 1))     
        goto CLPYRFin;
    if (alloc_acj(nmax + 1))     
        goto CLPYRFin;
    if (alloc_acn(nmax + 1))        /* number of elements */
        goto CLPYRFin;
    if (alloc_acns(nmax + 1))       /* times used */
        goto CLPYRFin;
    if (alloc_ack(nmax + 1))    
        goto CLPYRFin;
    if (alloc_acs(nmax + 1))    
        goto CLPYRFin;
    if (alloc_acl(nmax + 1))     
        goto CLPYRFin;
    if (alloc_acr(nmax + 1))     
        goto CLPYRFin;
    if (alloc_acyf(nmax + 1))       /* level */
        goto CLPYRFin;
    if (alloc_acc(n + 1))           /* mark */
        goto CLPYRFin;

    idx_e1 = AcI;                   /* index of first element */
    idx_e2 = AcJ;                   /* index of second element */
    idx_min = AcK;                  /* index of minimal element */
    idx_max = AcS;                  /* index of maximal element */
    idx_left = AcL;                 /* index of left neighbour */
    idx_right = AcR;                /* index of right neighbour */
    idx_mark = AcC;                 /* temporary membership marks */
   
    k = 1;
    for (i = 1; i < n; ++i) {
        for (j = 0; j < i; ++j)   
            AcXF[k++] = (float)gdd_adj(i,j,PMGN);                       
    }
    for (i = 1; i <= n; ++i) { 
        AcN[i] = 1;
        idx_min[i] = idx_max[i] = i;
    }


printf("n=%d nmax=%d nn=%d\n",n,nmax,nn);

    for (i = 2; i <= n; ++i) {
        ip = (i - 1) * (i - 2) / 2;
        for (j = 1; j < i; ++j)
            printf("%g ",AcXF[ip + j]);
        newline();
    }



    nc = n;         /* number of nodes */

    while (1) {
                                 
        dm = DBLMAX;
        isel = jsel = 0;

        for (i = 2; i <= nc; ++i) {
            if (AcNS[i] >= 2)
                continue;
            if (idx_left[i] && idx_right[i])
                continue;

            ip = ((i - 1) * (i - 2)) / 2;  

            for (j = 1; j < i; ++j) {
                if (AcNS[j] >= 2)
                    continue;

                d = (double)AcXF[ip + j];
                if (d >= dm)  
                    continue;

                /* check order requirements; only if i and/or j are
                   elementary objects */

                if (i <= n || j <= n) {
                    r = clpyr_ord(i,j,n,idx_left,idx_right,idx_min,idx_max);
                    if (r == 0)
                        continue;
                }
 
                /* check whether i and j are contained in the same cluster */

                r = imax(AcN[i],AcN[j]);

                sflag = 1;
                for (k = 1; k <= nc; ++k) {
                    if (AcN[k] < r)
                        continue;

                    for (l = 1; l <= n; ++l)
                        idx_mark[l] = 0;

                    clpyr_mark(k,idx_e1,idx_e2,idx_mark,n);
                    if (clpyr_mcheck(i,idx_e1,idx_e2,idx_mark)) 
                        continue;
                    if (clpyr_mcheck(j,idx_e1,idx_e2,idx_mark) == 0) {
                        sflag = 0;
                        break;
                    }
                }
                if (sflag == 0) {
                    printf("sflag = 0. continue\n");
                    continue;
                }

                /* check if connected */

                if (clpyr_con(i,j,idx_min,idx_max,idx_left,idx_right) == 0)
                    continue;

                /* check for additional condition */

                sflag = 1;
                for (l = 1; l <= n; ++l)
                    idx_mark[l] = 0;
                clpyr_mark(i,idx_e1,idx_e2,idx_mark,n);
                clpyr_mark(j,idx_e1,idx_e2,idx_mark,n);
                r = 0;
                for (l = 1; l <= n; ++l) {
                    if (idx_mark[l])
                        r++;
                }
                if (r >= 3) {
                    for (k = n + 1; k <= nc; ++k) {
                        if (k == i || k == j || AcN[k] < 2)
                            continue;

                        r = clpyr_cont(k,i,idx_e1,idx_e2);
                        printf("clpyr_cont k=%d i=%d r=%d\n",k,i,r);

                        if (r)
                            continue;

                        r = clpyr_cont(k,j,idx_e1,idx_e2);
                        printf("clpyr_cont k=%d j=%d r=%d\n",k,j,r);

                        if (r)
                            continue;

                        if (clpyr_mcheck(k,idx_e1,idx_e2,idx_mark) == 0) {
                            sflag = 0;
                            break;
                        }
                    }
                 
                    if (sflag == 0) {
                        printf("NOT OK i=%d j=%d\n",i,j   );
                        continue;
                    }
                }
 

printf("vorlaeufig ok: i=%d j=%d d=%g\n",i,j,d  );
                isel = i; 
                jsel = j;
                dm = d;
            }
        }

        if (isel == 0 || jsel == 0)     /* nor more fusions */
            break;

        /* isel and jsel have been selected for fusion */

        if (nc >= nmax) {
            printf1("Error: exceeded maximal storage (nmax=%d).\n",nmax);
            goto CLPYRFin;
        }
        nc++;
        idx_e1[nc] = isel;
        idx_e2[nc] = jsel;
        AcNS[isel] += 1;
        AcNS[jsel] += 1;

        if (isel <= n || jsel <= n) {
            r = clpyr_ord(isel,jsel,n,idx_left,idx_right,idx_min,idx_max);
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
        AcYF[nc] = (float)dm;

        /* number of elements */

        for (l = 1; l <= n; ++l)
            idx_mark[l] = 0;
        clpyr_mark(nc,idx_e1,idx_e2,idx_mark,n);
        for (l = 1; l <= n; ++l) {
            if (idx_mark[l])
                AcN[nc] += 1;
        }


printf("nc=%d : ",nc);
for (l = 1; l <= n; ++l)
    printf("%d ",idx_mark[l]);
newline();


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
        

        printf("nc=%2d i=%2d j=%2d dm=%lg acn[nc]=%2d k[nc]=%2d s[nc]=%2d l(i)=%2d r(i)=%2d l(j)=%2d r(j)=%2d r=%d y=%g\n",
                                   nc,isel,jsel,dm,AcN[nc],idx_min[nc],idx_max[nc],
                                   idx_left[isel],idx_right[isel],idx_left[jsel],idx_right[jsel],r,AcYF[nc] );     


        /* update dissimilarity matrix */

        ip = (nc - 1) * (nc - 2) / 2;
        for (k = 1; k <= nc; ++k) {
            AcXF[++ip] = (float)dmax(clpyr_dis(k,isel,AcXF),clpyr_dis(k,jsel,AcXF));
        }
    }

    printf("Final dis matrix\n");
    for (i = 1; i <= nc; ++i) {
        printf("%4d : ",i);
        for (j = 1; j <= nc; ++j) {
            printf("%lg ",clpyr_dis(i,j,AcXF));
        }
        newline();
    }

    fprintf(PMFd,"# 1. cluster number\n");
    fprintf(PMFd,"# 2. first element\n");
    fprintf(PMFd,"# 3. second element\n");
    fprintf(PMFd,"# 4. number of objects\n");
    fprintf(PMFd,"# 5. minimal element\n");
    fprintf(PMFd,"# 6. maximal element\n");
    fprintf(PMFd,"# 7. next on left side\n");
    fprintf(PMFd,"# 8. next on right side\n");
    fprintf(PMFd,"# 9. index\n\n");
             
    for (i = 1; i <= nc; ++i) {
        fprintf(PMFd,PMNFmtS,i);
        fprintf(PMFd,PMNFmtS,idx_e1[i]);
        fprintf(PMFd,PMNFmtS,idx_e2[i]);
        fprintf(PMFd,PMNFmtS,AcN[i]);
        fprintf(PMFd,PMNFmtS,idx_min[i]);
        fprintf(PMFd,PMNFmtS,idx_max[i]);
        fprintf(PMFd,PMNFmtS,idx_left[i]);
        fprintf(PMFd,PMNFmtS,idx_right[i]);
        fprintf(PMFd,PMFmtS,AcYF[i]);
        fprintf(PMFd,"\n");
    }
    printf1("%d records written to: %s\n",nc + 10,PMFdName);

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

    AcK[1] = l;
    for (i = 2; i <= n; ++i) {
        l = idx_right[l];
        if (l < 1 || l > n)
            goto CLPYRCon;
        AcK[i] = l;
    }
    sflag = 0;

    printf1("Order: ");
    for (i = 1; i <= n; ++i) 
        printf1("%d ",AcK[i]);
    newline();

    if (PMF1Def) {   /* print new distance matrix */

        if (alloc_aczf(n * (n - 1) / 2 + 1))                     
            goto CLPYRFin;

        for (i = 2; i <= n; ++i) {
            for (j = 1; j < i; ++j) {
                k = clpyr_ccl(i,j,n,nc,idx_e1,idx_e2);
                if (k == 0) {
                    printf("Error in data structure clpyr_ccl (i=%d, j=%d).\n",i,j);
                    goto CLPYRFin;
                }
printf("i=%d j=%d k=%d d=%lg \n",i,j,k,AcYF[k] );
                AcZF[(i - 1) * (i - 2) / 2 + j] = AcYF[k];
            }
        }
        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= n; ++j) {
                if (i == j)
                    d = 0.0;
                else if (j < i)
                    d = AcZF[(i - 1) * (i - 2) / 2 + j];
                else             
                    d = AcZF[(j - 1) * (j - 2) / 2 + i];

                fprintf(PMF1d,PMFmtS,d);
            }
            fprintf(PMF1d,"\n");
        }
        printf1("%d records written to: %s\n",n,PMF1dName);
    }


CLPYRCon:
    if (sflag)
        printf1("Error: cannot find a valid ordering.\n");

    /* check for compatible ordering */




    if (PMPCFDef) {                      /* plot dendrogram */
        if (alloc_aczf(nc + 1))                                
            goto CLPYRFin;
        clpyr_den(n,nc,AcK,idx_e1,idx_e2,AcYF,AcZF);
    }


    err = 0;

CLPYRFin:       
    p_clean();
    return(err);
}
   
/* return distance between node i and j */

double clpyr_dis(int i,int j,float *d)
{
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

int clpyr_con(int i,int j,int *idx_min,int *idx_max,int *idx_left,int *idx_right)
{
    int il,ir,jl,jr;

    il = idx_min[i];
    jl = idx_min[j];
    ir = idx_max[i];
    jr = idx_max[j];

printf("clpyr_con: i=%d j=%d il=%d ir=%d jl=%d jr=%d ",i,j,il,ir,jl,jr);


    if (clpyr_left(ir,jl,idx_left)) {
        printf("ir left from jl\n");
        printf("idx_right(ir) = %d\n",idx_right[ir]);
                             
        if (idx_right[ir] != jl) {
printf(" return 0\n");
            return(0);
        }
    }
    if (clpyr_left(jr,il,idx_left)) {
        printf("jr left from il\n");
        printf("idx_right(jr) = %d\n",idx_right[jr]);
                             
        if (idx_right[jr] != il) {
printf(" return 0\n");
              return(0);
        }
    }
printf(" return 1\n");
    return(1);
}

/* return 1 if i << j possible, 2 if j << i possible, 0 if nothing possible */

int clpyr_ord(int i,int j,int n,int *idx_left,int *idx_right,int *idx_min,int *idx_max)
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
                if (idx_left[j] || clpyr_left(j,i,idx_left))
                    return(0);
                else
                    return(1);
            }
            else {
                if (idx_left[jl] || clpyr_left(jl,i,idx_left))
                    return(0);
                else
                    return(1);
            }
        }
        else if (idx_right[i]) {
            if (j <= n) {
                if (idx_right[j] || clpyr_right(j,i,idx_right))
                    return(0);
                else
                    return(2);
            }
            else {
                if (idx_right[jr] || clpyr_right(jr,i,idx_right))
                    return(0);
                else
                    return(2);
            }
        }
        else {
            if (j <= n) {
                if (idx_left[j] == 0 && clpyr_left(i,j,idx_left) == 0)
                    return(1);
                else if (idx_right[j] == 0 && clpyr_right(i,j,idx_right) == 0)
                    return(2);
                else
                    return(0);
            }
            else {
                if (idx_left[jl] == 0 && clpyr_left(i,jl,idx_left) == 0)
                    return(1);
                else if (idx_right[jr] == 0 && clpyr_right(i,jr,idx_right) == 0)
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
            if (idx_left[il] || clpyr_left(j,il,idx_left))
                return(0);
            else
                return(1);
        }
        else if (idx_right[j]) {
            if (idx_right[ir] || clpyr_right(j,ir,idx_right))
                return(0);
            else
                return(2);
        }
        else {
            if (idx_left[il] == 0 && clpyr_left(j,il,idx_left) == 0)
                return(2);
            else if (idx_right[ir] == 0 && clpyr_right(j,ir,idx_right) == 0)
                return(1);
            else
                return(0);
        }
    }
}

/* return 1 if i is placed left of j */                            

int clpyr_left(int i,int j,int *idx_left)
{
    while (j) {
        j = idx_left[j];
        if (j == i)
            return(1);
    }
    return(0);
}

/* return 1 if i is placed right of j */                            

int clpyr_right(int i,int j,int *idx_right)
{
    while (j) {
        j = idx_right[j];
        if (j == i)
            return(1);
    }
    return(0);
}

/* return 1 if i is not contained in idx_mark */

int clpyr_mcheck(int i,int *idx_e1,int *idx_e2,char *idx_mark)
{
    register int ii,jj;

    if (i == 0)
        return(0);
    ii = idx_e1[i];
    jj = idx_e2[i];
    if (ii) {
        if (clpyr_mcheck(ii,idx_e1,idx_e2,idx_mark))
            return(1);
        if (clpyr_mcheck(jj,idx_e1,idx_e2,idx_mark))
            return(1);
    }
    else if (idx_mark[i] != 1)
        return(1);
    return(0);
}

/* mark idx_mark[l] = 1 if l is contained in i */

void clpyr_mark(int i,int *idx_e1,int *idx_e2,char *idx_mark,int n)
{
    register int ii,jj;

    if (i == 0)
        return;
    ii = idx_e1[i];
    jj = idx_e2[i];
    if (ii) {
        clpyr_mark(ii,idx_e1,idx_e2,idx_mark,n);
        clpyr_mark(jj,idx_e1,idx_e2,idx_mark,n);
    }
    else  
        idx_mark[i] = 1;            
    return;
}

/* return 1 if i is contained in j */

int clpyr_cont(int i,int j,int *idx_e1,int *idx_e2)
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

    if (clpyr_cont(i,j1,idx_e1,idx_e2))
        return(1);
    if (clpyr_cont(i,j2,idx_e1,idx_e2))
        return(1);
    return(0);                  
}

/* return index of lowest cluster which contains i and j */

int clpyr_ccl(int i,int j,int n,int nc,int *idx_e1,int *idx_e2)
{
    register int k;
         
    for (k = n + 1; k <= nc; ++k) {
        if (clpyr_cont(i,k,idx_e1,idx_e2) && clpyr_cont(j,k,idx_e1,idx_e2))
            return(k);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  clpyr_den   Create a command file for plotting the pyramid.             */

void clpyr_den(int n,int nc,int *iord,int *idx_e1,int *idx_e2,float *acy,float *tmp)
{
    register int i,j,k,l;
    int ic,jc;
    double d,hy,h,hl,hr,xl,xr; 

    hy = 0.0;
    for (i = 1; i <= nc; ++i) {
        if (hy < acy[i])
            hy = acy[i];
    }
    fprintf(PMPCFd,"psfile = %s.ps;\n",PMPCFName);
    fprintf(PMPCFd,"psetup(pxa=0,%d,pya=0,%g);\n",n + 1,hy);
    fprintf(PMPCFd,"plxa;\n");
    fprintf(PMPCFd,"plya;\n");
    d = -hy / 20;
    for (i = 1; i <= n; ++i)  
        fprintf(PMPCFd,"pltext(xy=%d,%lg,fs=1.5)=%d;\n",i,d,iord[i]);

    for (i = 1; i <= n; ++i)
        tmp[iord[i]] = (float)i;

    for (k = n + 1; k <= nc; ++k) {
        i = idx_e1[k];
        j = idx_e2[k];

        h = acy[k];
        xl = (double)tmp[i];
        hl = (double)acy[i];
        xr = (double)tmp[j];
        hr = (double)acy[j];
        d = (xr - xl) / 4.0;

        fprintf(PMPCFd,"plotp=%g,%g,%g,%g,%g,%g,%g,%g;\n",    
            xl,hl,xl + d,h,xr - d,h,xr,hr);

        tmp[k] = (xl + xr) / 2.0;

    }

    printf1("Output file for plot of dendrogram: %s\n",PMPCFName);
}

