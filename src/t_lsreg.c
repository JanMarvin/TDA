/****************************************************************************/
/*  t_lsreg                                                                 */
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
#include "t_matf.h"
#include "t_sort.h"
#include "t_svd.h"
#include "t_ml.h"
#include "t_gf.h"
#include "t_freq.h"

/*  functions in t_lsreg.c */

int lsreg(void);
void lsreg_res(int nx,double *b,int cov,int nw,double *c,double sig);
void lsreg_pdat(int mw,int nw,int ne,int ni,int nif,int nv,int opt);
int lsreg_dgrp(int opt);
void lsreg_d1grp(double *beta,double *cov,int nc,int ni);
int lsreg_ncov(int nw,int nx1,double *b,double sig);
int zreg(void);
int zreg1(void);


double *DGRPW[MaxDGRP];     /* weights */
short DGRPWA[MaxDGRP];      /* if allocated */
short DGRPNV[MaxDGRP];      /* number of variables in DGRPVIdx[] */     
short *DGRPVIdx[MaxDGRP];   /* variables */
char *DGRPDef[MaxDGRP];     /* definition */
short DGRPDefA[MaxDGRP];    /* allocated */
int DGRPDLen = 0;           /* length of definition */
int NDGRP = 0;              /* number of groups */

/* ------------------------------------------------------------------------ */
/*  lsreg()         Least squares regression.                               */
/*                                                                          */
/*  lsreg(                                                                  */
/*      w=...,          case weight variable, def. no weights               */
/*      ni=1,           if without intercept, def. ni=0                     */
/*      lsecon=...,     equality constraints                                */
/*      lsicon=...,     inequality constraints                              */
/*      dgrp=...,       estimate dummy variables with constraints           */
/*      tfmt=...,       print format for results, def. tfmt=10.4            */
/*      s=1,            use robust covariance matrix, def. s=0              */
/*      df=...,         print data to an output file                        */
/*      fmt=...,        print format for df option, def. 10.4               */
/*      dtda=...,       TDA description file for df option                  */
/*      ppar=...,       print estimated coefficients to output file         */
/*      pcov=...,       print covariance matrix to output file              */
/*      pres=...,       print residuals to output file                      */
/*      mfmt=...,       print format for pcov and pres option               */
/*      mplog=...,      write norm of residuals into matrix                 */
/*      mppar=...,      write parameters into matrix                        */
/*      mpcov=...,      write covariance matrix into matrix                 */
/*  ) = varlist;                                                            */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int lsreg(void)
{
    register int i,j,k,l;
    int ll,err,nv,ncon,nx,nx1,nw,mw,ne,ni;
    int cov,ierr,nn,ranka,ranke,rdef,df;
    double tmp,rnorme,rnorml,sse,ssq,rr,fs,sy,syy,sig,wt;

    err = -1;
    if (check_cmd(2))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    printf1("Least squares regression. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,4,1))     /* get parameters */
        goto LSRFin;

    if (PMFmtF == 0)
        pmfmt(10,4);

    newline();
    if (PMNW != 1) {
        printf1("Error: lsreg command requires cross-section data (nw=1).\n");
        goto LSRFin;
    }
    nv = check_nvar(0);         /* check variables */
    if (nv == 0)                /* PMNV is number of variables */       
        goto LSRFin;        
       
    prn_nwvar(0);               /* print variables */
    nx = nx1 = nv - 1;
    if (PMNI)  
        printf1("Model without intercept.\n");
         
    else {
        PMNI = 0;
        nx1++;
        if (nx == 0)  
            printf1("Model without independent variables.\n");
    }
    if (DGRPFlg) {              /* process dgrp option */
        if (PMWVar >= 0) {
            printf1("Error: dgrp option not compatible with case weights.\n");
            goto LSRFin;
        }
        if (lsreg_dgrp(1))
            goto LSRFin;
    }
    if (PMS != 0 && PMWVar >= 0) {
        printf1("Error: s=1 option not compatible with case weights.\n");
        goto LSRFin;   
    }

    /* allocate data matrix w[] and solution vector x[] */

    ne = ni = ncon = 0;                 /* number of constraints */
    nw = nx1 + 1;                       /* number of columns */
    mw = NOC + NCONSTR + NDGRP;         /* number of rows */

    if (alloc_acw(mw * nw + 1))
        goto LSRFin;

    if (alloc_acx(nx1 + 1))
        goto LSRFin;

    if (alloc_acy(nx1 + 1))
        goto LSRFin;

    if (alloc_acv(nx1 * nx1 + 1))   /* used for cov matrix */
        goto LSRFin;

    if (NCONSTR > 0) {
        if (p_con(CmdBuf,nx,PMNI,mw,nw,AcW,&ne,&ni))      /* get constraints */
            goto LSRFin;
    }
    if (NDGRP > 0) {        /* add constraints to data matrix */

        for (i = 0; i < NDGRP; ++i) {
            k = (ne + i) * nw + 1;
            if (PMNI == 0)
                k++;

            for (j = 1; j <= nx; ++j)  
                AcW[k++] = DGRPW[i][j];
        }
        ne += NDGRP;
    }
    ncon = ne + ni;

    printf1("Equality constraints: %d\n",ne);        
    printf1("Inequality constraints: %d\n\n",ni);          

    printf1("Reading data. Cases: %d\n",NOC);
               
    if (SVEFlg)
        p_warn(-1,1);

    if (PMWVar >= 0)  
        printf1("Using weights defined by: %s\n",VName[PMWVar]);

    /* get data into AcW, also calculate sy (sum of Y) and syy (sum of YY) */

    sy = syy = 0.0;
    k = ne;
    wt = 1.0;
    for (i = 0; i < NOC; ++i) {                 /* number of cases */

        if (PMWVar >= 0) {
            wt = get_data(PMWVar,i);
            if (wt < 0.0) {
                printf1("Error: found negative weight in case %d\n",i + 1);
                goto LSRFin;
            }
            if (wt > 0.0)
                wt = sqrt(wt);
        }

        if (PMNI == 0 && k >= ne && k < mw - ni)
            AcW[k * nw + 1] = wt;

        tmp = get_data(PMVIdx[0],i) * wt;  
        AcW[k * nw + nw] = tmp;                         
        sy += tmp;
        syy += tmp * tmp;

        j = 1;           
        if (PMNI)
            ll = 1;
        else
            ll = 2;
        for (l = 0; l < nx; ++l) {
            AcW[k * nw + ll] = get_data(PMVIdx[j],i) * wt;  
            j++;         
            ll++;
        }
        k++;
    }
    if (PMF1Def)                               /* write data matrix */
        lsreg_pdat(mw,nw,ne,ni,PMNI,nv,0);

    newline();
               
    /* the covariance matrix is only calculated when: there are no 
       inequality constraints and the number of cases
       is greater than the number of variables */

    nn = mw - ncon;     /* number of least squares cases */

    if (ni > 0 || nn <= nx1) {
        cov = 0;    
        printf1("Covariance matrix not calculated.\n");
    }
    else
        cov = 1;

    /*  solve least squares problem */

    ierr = lsei(AcW,ne,nn,ni,nx1,AcX,cov,&rnorme,&rnorml,&ranka,&ranke);

    if (ierr == -1) {
        printf1("Equality constraints are contradictory.\n");
        printf1("Calculated a least squares solution.\n");
    }
    else if (ierr) {

        if (ierr == -2)  
            printf1("Inequality constraints are incompatible.\n");
        else if (ierr == -3)  
            printf1("Equality and inequality constraints are contradictory.\n");
        else if (ierr == -4)  
            p_err(-2,1);
        else  
            printf1("LSEI error return: %d\n",ierr);
        goto LSRFin;  
    }
    rdef = 0;
    if (ranka != nx1 - ranke) {    /* set if rank deficient */
        rdef = 1;
        cov = df = 0;
    }
    else {
        df = nn - ranka;         
        if (df <= 0)
            cov = df = 0;
    }

    if (ne < nx1) {
        printf1("Rank of "); 
        if (ne > 0)
            printf1("reduced ");
        printf1("least squares data matrix: %d\n",ranka); 
    }
    if (rdef)
        printf1("Warning: data matrix is rank deficient.\n");

    if (ne < nx1)
        printf1("Norm of least squares residuals: %lg\n",rnorml); 

    if (PMMPLogDef == 1)                /* create matrix for norm of resid */
        mp_putlog(rnorml);
        
    if (PMMPParDef == 1)                /* create matrix for parameters */
        mp_putpar(nx1,AcX);

    if (ne > 0) {
        printf1("Rank of equality constraints: %d\n",ranke); 
        if (ne != ranke)
            printf1("Warning: equality constraints not linear independent.\n");

        printf1("Norm of residuals of equality constraints: %lg\n",rnorme);
    }
    newline();

    /*  If residuals are almost zero ignore standard errors */

    if (fabs(rnorme + rnorml) < EPSI1)  
        cov = 0;
       
    if (cov == 0)
        df = 0;
    else  
        printf1("Degrees of freedom: %d\n",df);

    /*  calculate sse = sum of squared residuals, also variance of
        residuals, squared multiple correlation, etc. */

    sig = 0.0;

    if (nn > 0) {       /* do only if there are least squares equations */

        sse  = rnorml * rnorml;
        printf1("Sum of squared residuals: %lg\n",sse); 

        tmp = 0.0;
        if (df > 0) {
            sig = sse / (double)df;
            printf1("Variance of residuals: %lg\n",sig);
        }
        ssq = syy - sy * sy / (double)nn;

        if (ssq > 0.0 && ni == 0 && df > 0) {     /* only without inequality constraints */
            rr = (ssq - sse) / ssq;
            if (rr >= 0.0)  
                printf1("Squared multiple correlation: %lg\n",rr);
               
            if (PMNI == 0 && df > 0) {
                tmp = 1.0 - (double)(nn - 1) * (1.0 - rr) / (double)df;
                if (tmp >= 0.0)  
                    printf1("Adjusted: %lg\n",tmp);
            }
        }
        if (ncon == 0 && rdef == 0) {   /* currently only without constraints */

            /*  Calculate Log likelihood */
            /*****************************
            tmp = (double)(-nn) * (L2PI + rlog(sse / (double)nn) + 1.0) / 2.0;
            printf1("\nLog likelihood: %lg\n,tmp);
            *********************/

            /*  Calculate F statistics */

            if (nx1 > 1 && df >= 1 && rr > 0.0 && rr < 1.0) {
                fs = (rr / (double)(nx1 - 1)) / ((1.0 - rr) / (double)df);  
                if (fs > 0.0) {
                    tmp = cdff(fs,nx1 - 1,df);
                    printf1("F-statistic: %lg\n",fs); 
                    printf1("Level of significance: %lg\n",tmp); 
                }
            }
        }
    }
    else
        cov = df = 0;

    if (cov) {      /* create new cov matrix in AcV, depending on PMS */

        if (PMResFDef)      /* write residuals */                 
            lsreg_res(nx,AcX,cov,nw,AcW,sig);

        if (lsreg_ncov(nw,nx1,AcX,sig))
            goto LSRFin;
   
        if (PMMPCovDef == 1)              /* create matrix for covariance */
            mp_putcov(nx1,AcV);

        if (PMCovFDef) {                          /* write cov matrix */
            prn_data(nx1,nx1,nx1,AcV,PMCovFd,PMMFmtS);
            p_wmsg(2,PMCovFName,1);
        }
      
        for (i = 1; i <= nx1; ++i)
            AcY[i] = AcV[(i - 1) * nx1 + i];
    }
    prn1_coeff(nx1,AcX,AcY,PMNI,df,PMVIdx,1); /* print coefficients */

    if (NDGRP > 0 && cov && df > 0)         /*  print adj st.dev. of coeff. */
        lsreg_d1grp(AcX,AcV,nx1,PMNI);

    mp_info();              /* print info about matrices */

    newline();
    err = 0;

LSRFin:
    lsreg_dgrp(0);
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  lsreg_res(nx,b)           Write residuals to PMResFd.                   */

void lsreg_res(int nx,double *b,int cov,int nw,double *c,double sig)
{
    register int i,j,l;
    int ll,nx1;
    double h,y,ye,se,se1,tmp;

    if (alloc_actmp(nx + 2))
        return;       

    if (sig < EPSI1)
        cov = 0;

    if (PMNI == 0) {
        AcTmp[1] = 1.0;
        ll = 1;
        nx1 = nx + 1;
    }
    else {
        ll = 0;
        nx1 = nx;
    }
    l = 1;
    fprintf(PMResFd,"# residuals created by lsreg command.\n");
    fprintf(PMResFd,"# c%-3d : case number\n",l++);
    l = 2;
    if (PMNI == 0)
        fprintf(PMResFd,"# c%-3d : constant one\n",l++);
    j = 1;     
    for (i = 0; i < nx; ++i)   
        fprintf(PMResFd,"# c%-3d : %s\n",l++,VName[PMVIdx[j++]]);
       
    fprintf(PMResFd,"# c%-3d : %s (dependent)\n",l++,VName[PMVIdx[0]]);
    fprintf(PMResFd,"# c%-3d : predicted\n",l++);
    fprintf(PMResFd,"# c%-3d : residual\n",l++);
    fprintf(PMResFd,"# c%-3d : leverage\n",l++);
    fprintf(PMResFd,"# c%-3d : std. dev. residual\n",l++);
    fprintf(PMResFd,"# c%-3d : stand. residual\n",l++);
    if (PMWVar >= 0)
        fprintf(PMResFd,"# c%-3d : case weights (%s)\n",l++,VName[PMWVar]);
    fprintf(PMResFd,"\n");

    for (i = 0; i < NOC; ++i) {                 /* number of cases */

        j = 1;           
        for (l = 1; l <= nx; ++l) {
            AcTmp[ll + l] = get_data(PMVIdx[j],i);          
            j++;         
        }
        fprintf(PMResFd,"%6d ",i + 1);
        for (l = 1; l <= nx1; ++l)
            fprintf(PMResFd,PMMFmtS,AcTmp[l]);

        y = get_data(PMVIdx[0],i);      /* dep variable */
        fprintf(PMResFd,PMMFmtS,y);

        ye = 0.0;
        for (l = 1; l <= nx1; ++l)  
            ye += AcTmp[l] * b[l];

        fprintf(PMResFd,PMMFmtS,ye);
        fprintf(PMResFd,PMMFmtS,y - ye);

        h = se1 = se = -1.0;
        if (cov) {              /* stand dev */
            tmp = 0.0;
            for (l = 1; l <= nx1; ++l) {
                for (j = 1; j <= nx1; ++j)  
                    tmp += c[(l - 1) * nw + j] * AcTmp[l] * AcTmp[j];
            }
            if (tmp > 0.0) {
                h = tmp / sig;
                tmp = sig - tmp;
                if (tmp > 0.0)  
                    se = sqrt(tmp);
                se1 = (y - ye) / se;
            }
        }
        fprintf(PMResFd,PMMFmtS,h);
        fprintf(PMResFd,PMMFmtS,se);
        fprintf(PMResFd,PMMFmtS,se1);

        if (PMWVar >= 0) {
            tmp = get_data(PMWVar,i);
            fprintf(PMResFd,PMMFmtS,tmp);
        }
        fprintf(PMResFd,"\n");
    }
    printf1("\nResiduals written to: %s\n",PMResFName);
}

/* ------------------------------------------------------------------------ */
/*  lsreg_pdat(mw,nw,ne,ni,nif,nv,opt)                                      */
/*                                                                          */
/*  Write data to output file. Optionally also dtda description file.       */

void lsreg_pdat(int mw,int nw,int ne,int ni,int nif,int nv,int opt)
{
    register int i,j,k,l;
    int nn;

    nn = mw - ne - ni;

    prn_data(nn,nw,nw,AcW + ne * nw,PMF1d,PMFmtS);
    p_wmsg(1,PMF1dName,1);
   
    if (PMTDAFDef) {
        fprintf(PMTDAFd,"# data written by lsreg command.\n");
        fprintf(PMTDAFd,"nvar(\n");
        fprintf(PMTDAFd,"  dfile = %s,\n",PMF1dName);
        fprintf(PMTDAFd,"  noc = %d,\n",nn);
        k = 0;
        if (nif == 0)  
            fprintf(PMTDAFd,"  INTZ [%d.%d] = c%-2d,\n",PMFmt1,PMFmt2,++k);
           
        l = 1;   
        for (i = 0; i < nv; ++i) {
            if (i == nv - 1)
                j = PMVIdx[0];
            else
                j = PMVIdx[opt + l++];

            fprintf(PMTDAFd,"  %s [%d.%d]",VName[j],PMFmt1,PMFmt2);
            if (VLabel[j] != NULL)                                                      
                fprintf(PMTDAFd,"(%s)",VLabel[j]);                                           
            fprintf(PMTDAFd," = c%d,\n",++k);
        }
        fprintf(PMTDAFd,");\n");
    
        printf1("TDA description written to: %s\n",PMTDAFName);
    }
}

/* ------------------------------------------------------------------------ */
/*  lsreg_dgrp(opt)     Process dgrp option. Weights are calculated in      */
/*                      DGRP[i][], i = 0,...,NDGRP - 1. If opt == 0         */
/*                      free previously allocated memory.                   */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int lsreg_dgrp(int opt)
{
    register int i,j,k,kk;
    int err,fin,nv,nx,nb;
    double tmp,w,sum;
    register char *p,*q,*s;

    if (opt == 0) {
        err = 0;
        goto LSDGRPFin;
    }   
    printf1("Processing dgrp option.\n");
    if (PMRHSTRA == 0) {
        p_err(-1,1);                 
        return(-1);
    }
    err = -1;
    p = PMRHSTR;
    fin = 0;
    nx = PMNV - 1;

    while (*p) {
        if (NDGRP >= MaxDGRP) {
            printf1("Error: exceeded max number of groups.\n");
            goto LSDGRPFin;
        }
        q = skip_nc(p);
        if (*p != '[' || *(q - 1) != ']') {
            p_err(-1,1);
            goto LSDGRPFin;
        }   
        if (!*q)
            fin = 1;
        else if (*q != ',') {
            p_err(-1,1);
            goto LSDGRPFin;
        }
        *q = '\0';

        printf1("Group selection: %s\n",p);
        s = p;
        p = get_nvia(p + 1,&nv,1,&nb);
        if (*p != ']') {
            p_err(-1,1);
            goto LSDGRPFin;
        }
        if (nb) {
            p_err(-42,1);
            goto LSDGRPFin;
        }   
        if (nv < 2) {
            if (nv > 0)
                printf1("Error: at least two variables from regression varlist must be specified.\n");
            goto LSDGRPFin;
        }
        k = strlen(s) + 1;
        if (DGRPDLen < k)
            DGRPDLen = k;

        if (!(DGRPDef[NDGRP] = (char *)calloc(k,sizeof(char)))) {
            p_err(-2,1);
            goto LSDGRPFin;
        }
        memrq(k,sizeof(char));
        DGRPDefA[NDGRP] = k;       
        strcpy(DGRPDef[NDGRP],s);

        if (!(DGRPVIdx[NDGRP] = (short *)calloc(nv,sizeof(short)))) {
            p_err(-2,1);
            goto LSDGRPFin;
        }
        memrq(nv,sizeof(short));
        DGRPNV[NDGRP] = nv;       

        kk = 0;
        for (i = 0; i < nv; ++i) {
            k = VLVIdx[i];
            for (j = 0; j < i; ++j) {
                if (k == VLVIdx[j]) {
                    printf1("Error: variables in constraints should be unique.\n");
                    goto LSDGRPFin;
                }
            }
            for (j = 1; j < PMNV; ++j) {
                if (k == PMVIdx[j]) {
                    DGRPVIdx[NDGRP][kk++] = j;
                    k = -1;
                    break;
                }
            }
            if (k >= 0) {
                printf1("Error: at least one variable is not a regressor variable.\n");
                goto LSDGRPFin;
            }
        }
        if (!(DGRPW[NDGRP] = (double *)calloc(nx + 1,sizeof(double)))) {
            p_err(-2,1);
            goto LSDGRPFin;
        }
        memrq(nx + 1,sizeof(double));
        DGRPWA[NDGRP] = nx + 1;

        sum = 0.0;
        for (i = 0; i < NOC; ++i) {                 /* number of cases */
            for (k = 0; k < nv; ++k) {
                kk = DGRPVIdx[NDGRP][k];
                tmp = get_data(PMVIdx[kk],i);       
                if (fabs(tmp) < EPSI1)
                    ;                
                else if (fabs(tmp - 1.0) < EPSI1) {
                    DGRPW[NDGRP][kk] += 1.0;
                    sum += 1.0;
                }
                else {
                    printf1("Error: variable %s in case %d isn't a dummy variable.\n",VName[PMVIdx[kk]],i + 1);
                    goto LSDGRPFin;
                }
            }
        }
        j = 10; 
        if (j < VNameLen)
            j = VNameLen;

        printf1("\nVariable  ");
        prnchar(' ',j - 10,0);
        printf1("   cases       weight\n");
        prnchar('-',21 + j,1);
        for (k = 0; k < nv; ++k) {
            kk = DGRPVIdx[NDGRP][k];
            i = PMVIdx[kk];
            printf1("%s ",VName[i]);
            prnchar(' ',j - strlen(VName[i]),0);
            w = tmp = DGRPW[NDGRP][kk];
            if (sum > 0.0)
                w /= sum;
            printf1("%7d %12.4lf\n",(int)tmp,w);
            DGRPW[NDGRP][kk] = w;
        }
        newline();
        if (sum < EPSI1) {
            printf1("Error: sum of weights is zero.\n");
            goto LSDGRPFin;
        }
        NDGRP++;

        if (fin)
            break;

        p = q + 1;
    }
    err = 0;

LSDGRPFin:
    if (err || opt == 0) {
        for (i = 0; i < MaxDGRP; ++i) {
            if (DGRPDefA[i] > 0) {    
                free((char *)DGRPDef[i]);            
                memrq(-DGRPDefA[i],sizeof(char));
                DGRPDefA[i] = 0;
            }
            if (DGRPWA[i] > 0) {    
                free((char *)DGRPW[i]);            
                memrq(-DGRPWA[i],sizeof(double));
                DGRPWA[i] = 0;
            }
            if (DGRPNV[i] > 0) {    
                free((char *)DGRPVIdx[i]);            
                memrq(-DGRPNV[i],sizeof(short));
                DGRPNV[i] = 0;
            }
        }
        NDGRP = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  lsreg_d1grp()       Print adj stand dev. of coefficients.               */
/*                                                                          */

void lsreg_d1grp(double *beta,double *cov,int nc,int ni)
{
    register int i,j,k;
    int iflag,ib,nv,len;
    double sig;                 

    iflag = 1;
    if (ni)
        iflag = 0;

    len = DGRPDLen;
    if (len < 8)
        len = 8;

    printf1("\nDGroup");
    prnchar(' ',len - 6,0);
    printf1("weighted adj. stand. dev.\n");
    prnchar('-',25 + len,1);
    for (i = 0; i < NDGRP; ++i) {
        nv = DGRPNV[i];
        sig = 0.0;
        for (j = 0; j < nv; ++j) {
            k = DGRPVIdx[i][j];
            ib = k + iflag;
            sig += DGRPW[i][k] * (beta[ib] * beta[ib] - cov[(ib - 1) * nc + ib]);
        }
        printf1("%s",DGRPDef[i]);
        prnchar(' ',10 + len - strlen(DGRPDef[i]),0);
        if (sig >= 0.0) {
            sig = sqrt(sig);
            printf1(PMTFmtS,sig);
        }  
        else
            printf1("  ---");
        newline();
    }
}

/* ------------------------------------------------------------------------ */
/*  lsreg_ncov(nw,nx1,b,sig)                                                */
/*                                                                          */
/*  nw = column dimension of AcW (containing stand cov matrix)              */
/*  nx1 = number of variables, including intercept.                         */
/*  note: PMNI == 1 if there is no intercept.                               */
/*  b[] contains parameters.                                                */  
/*  sig = variance of residuals                                             */  
/*                                                                          */
/*  Create covariance matrix in AcV. If PMS != 0, create a robust           */
/*  cov matrix based on White (1980).                                       */
/*                                                                          */
/*  Return 0 if OK, -1 if error (insufficient memory).                      */

int lsreg_ncov(int nw,int nx1,double *b,double sig)
{
    register int i,j,k,l,i1,j1;
    int nxx,ll;
    double tmp,tmp1,tmp2,ye;

    nxx = nx1 * nx1;

    if (PMS == 0) {
        k = 1;
        for (i = 0; i < nx1; ++i) {
            for (j = 1; j <= nx1; ++j) 
                AcV[k++] = AcW[i * nw + j];
        }
        return(0);
    }

    printf1("Calculating White's covariance matrix.\n");
    if (sig < EPSI1) {
        printf1("Error: variance of residuals almost zero.\n");
        return(-1);
    }
    sig *= sig;

    if (PMNI == 0)
        ll = 1;
    else
        ll = 0;

    /* first create AcU[] containing residuals */

    if (alloc_acu(NOC + 1))
        return(-1);   

    if (alloc_actmp(nx1 + 1))
        return(-1);   

    for (i = 0; i < NOC; ++i) {                 /* number of cases */

        ye = 0.0;
        i1 = 1;
        for (l = 1; l <= nx1; ++l) {
            if (PMNI == 0 && l == 1)
                ye += b[i1++];
            else  
                ye += b[i1++] * get_data(PMVIdx[l - ll],i);
        }
        AcU[i + 1] = get_data(PMVIdx[0],i) - ye;                     
    }

    /* create X' diag(res2) X matrix in AcTmp[] */

    if (alloc_actmp(nxx + 1))
        return(-1);   

    for (i1 = 1; i1 <= nx1; ++i1) {
        for (j1 = 1; j1 <= nx1; ++j1) {
            tmp = 0.0;
                     
            for (i = 0; i < NOC; ++i) {                 /* number of cases */

                tmp1 = tmp2 = 1.0;
                if (PMNI || i1 > 1)
                    tmp1 = get_data(PMVIdx[i1 - ll],i);

                if (PMNI || j1 > 1)
                    tmp2 = get_data(PMVIdx[j1 - ll],i);

                tmp += tmp1 * tmp2 * AcU[i + 1] * AcU[i + 1];
            }
            AcTmp[(i1 - 1) * nx1 + j1] = tmp;
        }
    }

    /* create new cov matrix in AcV */

    for (i1 = 1; i1 <= nx1; ++i1) {
        for (j1 = 1; j1 <= nx1; ++j1) {
            tmp = 0.0;
            for (k = 1; k <= nx1; ++k) {
                tmp1 = 0.0;
                for (l = 1; l <= nx1; ++l)  
                    tmp1 += AcTmp[(k - 1) * nx1 + l] * AcW[(l - 1) * nw + j1];
                tmp += tmp1 * AcW[(i1 - 1) * nw + k];
            }
            AcV[(i1 - 1) * nx1 + j1] = tmp / sig;
        }
    }
    alloc_acu(0);
    alloc_actmp(0);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  zreg()          Least squares regression with censored data             */
/*                  (Buckley-James).                                        */
/*                                                                          */
/*  zreg(                                                                   */
/*      yw=...,         name of censoring indicator (variable)              */
/*      ni=1,           if without intercept, def. ni=0                     */
/*      mxit=...,       max number of iterations, def. 20                   */
/*      tolp=...,       tolerance for convergence, def. 1e-6                */
/*      tfmt=...,       print format for results, def. tfmt=10.4            */
/*      ppar=...,       print estimated coefficients to output file         */
/*      df=...,         print dep. var. to an output file                   */
/*      fmt=...,        print format for df option, def. 10.4               */
/*      prot=...,       protocol file                                       */
/*      pfmt=...,       print format, -19.11                                */
/*  ) = Y,X1,...;       varlist, first dep. variable.                       */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int zreg(void)
{
    register int i,j,k,l;
    int ll,err,nv,nx,nx1,nw,iter;
    int conv,cov,ierr,nn,ranka,ranke,m;
    double rnorme,rnorml,delta,tmp,tmp1;

    err = -1;
    if (check_cmd(2))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    printf1("Regression with censored data. Current memory: %d bytes.\n",MemReq);

    TOLP = 1.e-6;

    if (parm(CmdBuf + 4,4,1))     /* get parameters */
        goto ZRFin;

    if (MxIter < 1 || MxItFlg == 0)
        MxIter = 20;

    if (PMFmtF == 0)
        pmfmt(10,4);

    newline();
    if (PMNW != 1) {
        printf1("Error: zreg command requires cross-section data (nw=1).\n");
        goto ZRFin;
    }
    nv = check_nvar(0);         /* check variables */
    if (nv == 0)                /* PMNV is number of variables */       
        goto ZRFin;        
       
    prn_nwvar(0);               /* print variables */

    if (PMYWVar < 0) {
        printf1("Error: need censoring information.\n");
        goto ZRFin;
    }
    printf1("Number of cases: %d\n",NOC);
    printf1("Censoring indicator: %s\n",VName[PMYWVar]);
    printf1("Max number of iterations: %d\n",MxIter);
    printf1("Tolerance for convergence: %lg\n\n",TOLG);

    nx = nx1 = nv - 1;
    if (PMNI)  
        printf1("Model without intercept.\n");
         
    else {
        PMNI = 0;
        nx1++;
        if (nx == 0)  
            printf1("Model without independent variables.\n");
    }

    /* allocate data matrix w[] and solution vector x[] */

    nw = nx1 + 1;                   /* number of columns */
    nn = NOC;                       /* number of rows */

    if (alloc_acw(nn * nw + 1))     /* least squares matrix */
        goto ZRFin;
    if (alloc_acx(nx1 + 1))         /* used for parameters */
        goto ZRFin;
    if (alloc_acy1(nx1 + 1))        /* used for previous parameters */
        goto ZRFin;
    if (alloc_acy(nn + 1))          /* used for dependent variable */
        goto ZRFin;
    if (alloc_acz(nn + 1))          /* used for residuals */
        goto ZRFin;
    if (alloc_acns(nn + 1))         /* used for censoring indicator */
        goto ZRFin;
    if (alloc_acn(nn + 1))          /* used for sorting pointer */
        goto ZRFin;
    if (alloc_acu(nn + 1))          /* used for jumps of K-M */
        goto ZRFin;

    for (i = 0; i < NOC; ++i) {
        AcY[i] = get_data(PMVIdx[0],i);  
        if ((int)get_data(PMYWVar,i) != 0)  
            AcNS[i] = 1;
    }

    /* begin iterations */

    if (SILENTFlg < 2)
        printfe("\n  Iter   Norm of Residuals     Par Change\n");

    if (PMProtFDef)  
        fprintf(PMProtFd,"Iteration  Parameter vector\n");

    cov = 0;    /* do not calculate covariance matrix */
    conv = 0;
    delta = 0.0;

    for (iter = 0; iter <= MxIter; ++iter) {
 
        /* get data into AcW  */
       
        k = 0;
        for (i = 0; i < NOC; ++i) {         
       
            if (PMNI == 0)
                AcW[k * nw + 1] = 1.0;

            AcW[k * nw + nw] = AcY[i];                      

            j = 1;           
            if (PMNI)
                ll = 1;
            else
                ll = 2;
            for (l = 0; l < nx; ++l) {
                AcW[k * nw + ll] = get_data(PMVIdx[j],i);  
                j++;         
                ll++;
            }
            k++;
        }

        /*  solve least squares problem */

        ierr = lsei(AcW,0,nn,0,nx1,AcX,cov,&rnorme,&rnorml,&ranka,&ranke);

        if (ierr) {
            if (ierr == -4)  
                p_err(-2,1);
            else  
                printf1("LSEI error return: %d\n",ierr);
            goto ZRFin;  
        }
        if (ranka != nx1 - ranke) {    /* set if rank deficient */
            printf1("Rank of least squares matrix: %d\n",ranka);
            printf1("Will not continue.\n");
            goto ZRFin;
        }
        if (iter == 0) {
            printf1("Regression results with original values.\n");
            printf1("Rank of least squares data matrix: %d\n",ranka); 
            printf1("Norm of least squares residuals: %lg\n",rnorml); 

            prn1_coeff(nx1,AcX,AcX,PMNI,0,PMVIdx,1); /* print coefficients */
            newline();
        }

        /* print parameters to protocol file */
   
        if (PMProtFDef) {
            fprintf(PMProtFd,"%7d    ",iter);
            for (j = 1; j <= nx1; ++j)
                fprintf(PMProtFd,PMPFmtS,AcX[j]);
            fprintf(PMProtFd,"\n");
        }
   
        /* check convergence */

        if (iter > 0) {
            delta = 0.0;
            for (j = 1; j <= nx1; ++j) {
                tmp = fabs(AcX[j] - AcY1[j]);
                if (AcX[j])
                    tmp /= fabs(AcX[j]);
                delta = dmax(delta,tmp);
            }
        }
        if (SILENTFlg < 2)  
            printfe("  %3d  %20.13e %17.10e\n",iter,rnorml,delta);

        if (iter > 0 && delta <= TOLP) {
            conv = 1;
            break;
        }
        for (j = 1; j <= nx1; ++j)  
            AcY1[j] = AcX[j];           

        /* update dependent variable. first calculate residuals in AcZ[] */

/**   
printf("VORHER\n");
        for (i = 0; i < NOC; ++i)  
printf("i=%d d=%d y=%lg\n",i,AcNS[i], AcY[i]);
newline();
**/        


        for (i = 0; i < NOC; ++i) {
            tmp = 0.0;
            j = 1;
            if (PMNI == 0) {
                tmp = AcX[j];
                j++;   
            }
            k = 1;
            while (j <= nx1) {
                tmp += AcX[j] * get_data(PMVIdx[k],i);  
                j++;
                k++;
            }
            AcZ[i] = get_data(PMVIdx[0],i) - tmp;
            if (AcNS[i] == 0)
                AcY[i] = tmp;

        }
        if (sortdp2a(NOC,AcZ,AcNS,AcN,0))      /* sort */
            goto ZRFin;    

/**      
printf("Erster Schritt\n");
        for (i = 0; i < NOC; ++i)  
printf("i=%d d=%d y=%lg\n",i,AcNS[i], AcY[i]);
newline();
**/    
        for (i = 0; i < NOC; ++i)  
            AcU[i] = 0.0;

        tmp = 1.0 / (double)NOC;
        m = NOC - 1;
        for (i = 0; i < NOC; ++i) {
            j = AcN[i];
            AcU[j] += tmp;
            if (AcNS[j] == 0 && i < NOC - 1) {
                tmp1 = AcU[j] / (double)m;
                for (k = i + 1; k < NOC; ++k)  
                    AcU[AcN[k]] += tmp1;
                AcU[j] = 0.0;
            }
            m--;
        }
        tmp = 0.0;
        for (i = 0; i < NOC; ++i) {
            k = AcN[i];
            tmp += AcU[k];
/*      
printf("i=%d AcNS[i]=%d z=%10.4lf  F=%10.4lf  U=%10.4lf\n",i,AcNS[i],
                                            AcZ[i],tmp,AcU[i]);
*/     
        }

        for (i = 0; i < NOC; ++i) {
            j = AcN[i];
            if (AcNS[j] == 0) {         /* if censored */
                if (i == NOC - 1) {
                    AcY[j] += AcZ[j];
                }
                else {
                    tmp = tmp1 = 0.0;
                    for (k = i + 1; k < NOC; ++k) {
                        l = AcN[k];
/*                      printf(".... l=%d d=%d z=%lg\n",l,AcNS[l],AcZ[l]); */
                        if (AcNS[l] != 0 || k == NOC - 1) {
                            tmp += AcZ[l] * AcU[l];
                            tmp1 += AcU[l];
                        }
                    }
                    if (tmp1 != 0.0)
                        AcY[j] += tmp / tmp1;
                }
            }
        }
/*    
printf("ENDE\n");
        for (i = 0; i < NOC; ++i)  
printf("i=%d d=%d y=%lg\n",i,AcNS[i], AcY[i]);
newline();
*/       

    }
    printf1("Convergence ");   
    if (conv == 0)
        printf1("not ");
    printf1("reached in %d iterations.\n",iter);

    printf1("Rank of least squares data matrix: %d\n",ranka); 
    printf1("Norm of least squares residuals: %lg\n",rnorml); 

    prn1_coeff(nx1,AcX,AcX,PMNI,0,PMVIdx,1); /* print coefficients */
    newline();

    if (PMF1Def) {
        for (i = 0; i < NOC; ++i) {
            tmp = get_data(PMVIdx[0],i);
            fprintf(PMF1d,"%6d ",i + 1);
            fprintf(PMF1d,PMFmtS,tmp);  
            fprintf(PMF1d,"%2d ",AcNS[i]);
            fprintf(PMF1d,PMFmtS,AcY[i]);
            fprintf(PMF1d,"\n");        
        }
        printf1("%d records written to: %s\n",NOC,PMF1dName);
    }
    err = 0;

ZRFin:
    p_clean();
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  zreg1()         Least squares regression with censored data             */
/*                  (Buckley-James).                                        */
/*                                                                          */
/*  zreg(                                                                   */
/*      yw=...,         name of censoring indicator (variable)              */
/*      xv=...,         time-varying covariates (dates)                     */
/*      mxit=...,       max number of iterations, def. 20                   */
/*      tolp=...,       tolerance for convergence, def. 1e-6                */
/*      tfmt=...,       print format for results, def. tfmt=10.4            */
/*      ppar=...,       print estimated coefficients to output file         */
/*      df=...,         print dep. var. to an output file                   */
/*      fmt=...,        print format for df option, def. 10.4               */
/*      prot=...,       protocol file                                       */
/*      pfmt=...,       print format, -19.11                                */
/*  ) = Y,X1,...;       varlist, first dep. variable.                       */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int zreg1(void)
{
    register int i,j,k,l;
    int err,nv,nx,nxu,nw,nxv,nk,iter,m,n,nn;
    int conv,ierr,ranka,ranke,ncen;
    double rnorme,rnorml,delta,tmp,tmp1,time;

    err = -1;
    if (check_cmd(2))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    printf1("Regression with censored data. Current memory: %d bytes.\n",MemReq);

    TOLP = 1.e-6;

    if (parm(CmdBuf + 5,4,1))     /* get parameters */
        goto ZR1Fin;

    if (MxIter < 1 || MxItFlg == 0)
        MxIter = 20;

    if (PMFmtF == 0)
        pmfmt(10,4);

    newline();
    if (PMNW != 1) {
        printf1("Error: zreg command requires cross-section data (nw=1).\n");
        goto ZR1Fin;
    }
    nv = check_nvar(0);         /* check variables */
    if (nv == 0)                /* PMNV is number of variables */       
        goto ZR1Fin;        
       
    prn_nwvar(0);               /* print variables */

    if (PM1NV > 0) {
        printf1("Plus time-varying indicator variables: %s",VName[PM1VIdx[0]]);
        for (i = 1; i < PM1NV; ++i)
            printf1(", %s",VName[PM1VIdx[i]]);
        newline();
        newline();
    }

    if (PMYWVar < 0) {
        printf1("Error: need censoring information.\n");
        goto ZR1Fin;
    }
/***###*************/

    for (i = 0; i < NOC; ++i) {
        k = (int)get_data(PMVIdx[0],i);        
        l = (int)get_data(PMVIdx[1],i);        
        printf("%4d %d %6d %6d",i + 1,
                (int)get_data(PMYWVar,i),k,l);   

        for (j = 0; j <= 428; ++j) {
            if (l > j || j >= k    )
                printf(" 0");
            else
                printf(" 1");
        }
        newline();
    }
    goto ZR1Fin;
/********************/








    nxu = nv - 1;
    nxv = PM1NV;
    nx = 1 + nxu + nxv;

    printf1("Number of cases: %d\n",NOC);
    printf1("Number of time-independent covariates: %d\n",nxu);
    printf1("Number of time-varying covariates: %d\n",nxv);
    printf1("Censoring indicator: %s\n",VName[PMYWVar]);
    printf1("Max number of iterations: %d\n",MxIter);
    printf1("Tolerance for convergence: %lg\n\n",TOLG);

    /* allocate data matrix w[] and solution vector x[] */

    nw = nx + 1;                    /* number of columns */

    if (alloc_acw(NOC * nw + 1))    /* least squares matrix */
        goto ZR1Fin;
    if (alloc_acx(nx + 1))          /* used for parameters */
        goto ZR1Fin;
    if (alloc_acy1(nx + 1))         /* used for previous parameters */
        goto ZR1Fin;
    if (alloc_acy(NOC + 1))         /* used for dependent variable */
        goto ZR1Fin;
    if (alloc_acz(NOC + 1))         /* used for residuals */
        goto ZR1Fin;
    if (alloc_acns(NOC + 1))        /* used for censoring indicator */
        goto ZR1Fin;
    if (alloc_acms(NOC + 1))        /* used for censoring indicator */
        goto ZR1Fin;
    if (alloc_acn(NOC + 1))         /* used for sorting pointer */
        goto ZR1Fin;
    if (alloc_acj(NOC + 1))         /* used for sorting pointer */
        goto ZR1Fin;
    if (alloc_acm(nxv + 1))         /* used for counting covariate values */
        goto ZR1Fin;
    if (alloc_acu(NOC + 1))         /* used for jumps of K-M */
        goto ZR1Fin;

    ncen = 0;
    for (i = 0; i < NOC; ++i) {
        AcY[i] = get_data(PMVIdx[0],i);  
        if ((int)get_data(PMYWVar,i) != 0)  
            AcNS[i] = 1;
        else
            ncen++;

        if (AcY[i] <= 0.0) {
            printf1("Error: found non-positive value in record %d\n",i + 1);
            goto ZR1Fin;
        }
    }

    if (sortdp(NOC,AcY,AcN))      /* sort */
        goto ZR1Fin;   

    printf1("Time      N    Cen ");
    l = 42;
    for (j = 0; j < nxv; ++j) {              
        k = PM1VIdx[j];
        prnchar(' ',6 - strlen(VName[k]),0);
        printf1("%s ",VName[k]);
        l += imax(7,strlen(VName[k]));
    }
    printf1(" Iter  Model parameters\n");
    prnchar('-',l,1);

    time = 0.0;
    nk = 0;
    n = NOC;
    while (n > 10) {
        while (nk < NOC) {
            j = AcN[nk];
            if (AcY[j] > time)
                break;
            if (AcNS[j] == 0)
                ncen--;
            n--;    
            nk++;
        }       
        printf1("%4d %6d %6d ",(int)time,n,ncen);


        /* create least squares matrix */

        for (l = 0; l < nxv; ++l)
            AcM[l] = 0;

        nn = 0;
        for (i = nk; i < NOC; ++i) {
            j = AcN[i];
            AcW[nn * nw + nw] = AcY[j] - time;               
            AcW[nn * nw + 1] = 1.0;
       
            for (l = 1; l <= nxu; ++l)  
                AcW[nn * nw + 1 + l] = get_data(PMVIdx[l],j);  
     
            for (l = 0; l < nxv; ++l) {
                tmp = get_data(PM1VIdx[l],j);          
                if (tmp <= time) {
                    AcM[l] += 1;
                    tmp = 1.0;
                }
                else
                    tmp = 0.0;

                AcW[nn * nw + 2 + nxu + l] = tmp;                       
            }
            nn++;
        }
        if (n != nn) {
            printf1("\nFatal error: n=%d nn=%d\n",n,nn);
            goto ZR1Fin;
        }
        for (l = 0; l < nxv; ++l)
            printf1("%6.4lf ",(double)AcM[l]/(double)n);
/**
printf("lsmatrix\n");
for (i = 0; i < nn; ++i) {
    for (j = 1; j <= nw; ++j)
        printf("%lf ",AcW[i * nw + j]);
    newline();
}
**/
        conv = 0;
        delta = 0.0;

        for (iter = 0; iter <= MxIter; ++iter) {
/**
printf("\niter=%d\n",iter);
printf("lsmatrix\n");
for (i = 0; i < nn; ++i) {
    for (j = 1; j <= nw; ++j)
        printf("%lf ",AcW[i * nw + j]);
    newline();
}
**/
            ierr = lsei(AcW,0,nn,0,nx,AcX,0,&rnorme,&rnorml,&ranka,&ranke);
            if (ierr) {
                newline();
                if (ierr == -4)  
                    p_err(-2,1);
                else  
                    printf1("LSEI error return: %d\n",ierr);
                goto ZR1Fin;  
            }
            if (ranka != nx - ranke) {    /* set if rank deficient */
                printf1("Rank of least squares matrix: %d\n",ranka);
                printf1("Will not continue.\n");
                goto ZR1Fin;
            }

            /* check convergence */

            if (iter > 0) {
                delta = 0.0;
                for (j = 1; j <= nx; ++j) {
                    tmp = fabs(AcX[j] - AcY1[j]);
                    if (AcX[j])
                        tmp /= fabs(AcX[j]);
                    delta = dmax(delta,tmp);
                }
                if (delta <= TOLP) {
                    conv = 1;
                    break;
                }
            }
            for (j = 1; j <= nx; ++j)  
                AcY1[j] = AcX[j];           

            /* update dependent variable  */
            /* first recreate least squares matrix and calculate residuals
               in AcZ[] */

            nn = 0;
            for (i = nk; i < NOC; ++i) {
                j = AcN[i];
                AcW[nn * nw + nw] = AcY[j] - time;               
                AcW[nn * nw + 1] = 1.0;
                tmp1 = AcX[1];

                for (l = 1; l <= nxu; ++l) {
                    tmp  = get_data(PMVIdx[l],j);  
                    AcW[nn * nw + 1 + l] = tmp;                     
                    tmp1 += tmp * AcX[l + 1];
                }

                for (l = 0; l < nxv; ++l) {
                    tmp = get_data(PM1VIdx[l],j);          
                    if (tmp <= time)   
                        tmp = 1.0;
                    else
                        tmp = 0.0;
                    AcW[nn * nw + 2 + nxu + l] = tmp;                       
                    tmp1 += tmp * AcX[2 + nxu + l];
                }
                AcZ[nn] = AcY[j] - time - tmp1;
                AcMS[nn] = AcNS[j];

                if (AcNS[j] == 0)
                    AcW[nn * nw + nw] = tmp1;                      

                nn++;
            }
/**
printf("##  2-lsmatrix\n");
for (i = 0; i < nn; ++i) {
    for (j = 1; j <= nw; ++j)
        printf("%lf ",AcW[i * nw + j]);
printf("z=%lg AcMS=%d ",AcZ[i],AcMS[i]   );
    newline();
}

**/





            if (sortdp2a(nn,AcZ,AcMS,AcJ,0))      /* sort */
                goto ZR1Fin;   
       
            for (i = 0; i < nn; ++i)  
                AcU[i] = 0.0;

            tmp = 1.0 / (double)nn;
            m = nn - 1;
            for (i = 0; i < nn; ++i) {
                j = AcJ[i];
                AcU[j] += tmp;
                if (AcMS[j] == 0 && i < nn - 1) {
                    tmp1 = AcU[j] / (double)m;
                    for (k = i + 1; k < nn; ++k)  
                        AcU[AcJ[k]] += tmp1;
                    AcU[j] = 0.0;
                }
                m--;
            }
            tmp = 0.0;
            for (i = 0; i < nn; ++i) {
                k = AcJ[i];
                tmp += AcU[k];
            }
            for (i = 0; i < nn; ++i) {
                j = AcJ[i];
                if (AcMS[j] == 0) {         /* if censored */
                    if (i == nn - 1) {
                        AcW[j * nw + nw] += AcZ[j];                      
                    }
                    else {
                        tmp = tmp1 = 0.0;
                        for (k = i + 1; k < nn; ++k) {
                            l = AcJ[k];
                            if (AcMS[l] != 0 || k == nn - 1) {
                                tmp += AcZ[l] * AcU[l];
                                tmp1 += AcU[l];
                            }
                        }
                        if (tmp1 != 0.0) {
                            AcW[j * nw + nw] += tmp / tmp1;                  
                        }
                    }
                }
            }



        } /* end of iteration */


/**

newline();
    printf1("Convergence ");   
    if (conv == 0)
        printf1("not ");
    printf1("reached in %d iterations.\n",iter);

    printf1("Rank of least squares data matrix: %d\n",ranka); 
    printf1("Norm of least squares residuals: %lg\n",rnorml); 

    prn1_coeff(nx,AcX,AcX,PMNI,0,PMVIdx,1); 
    newline();
**/

    printf1("%5d ",iter);
    for (j = 1; j <= nx; ++j)
        printf1(PMTFmtS,AcX[j]);
    newline();
/*
    if (time >= 12.0)
        break;

*/
        time += 1.0;
    }

    err = 0;

ZR1Fin:
    p_clean();
    return(err);
}






