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
#include "tda_context.h"

/*  functions in t_lsreg.c */

int lsreg(TDAContext *ctx);
void lsreg_res(TDAContext *ctx, int nx,double *b,int cov,int nw,double *c,double sig);
void lsreg_pdat(TDAContext *ctx, int mw,int nw,int ne,int ni,int nif,int nv,int opt);
int lsreg_dgrp(TDAContext *ctx, int opt);
void lsreg_d1grp(TDAContext *ctx, double *beta,double *cov,int nc,int ni);
int lsreg_ncov(TDAContext *ctx, int nw,int nx1,double *b,double sig);
int zreg(TDAContext *ctx);
int zreg1(TDAContext *ctx);



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

int lsreg(TDAContext *ctx)
{
    register int i,j,k,l;
    int ll,err,nv,ncon,nx,nx1,nw,mw,ne,ni;
    int cov,ierr,nn,ranka,ranke,rdef,df;
    double tmp,rnorme,rnorml,sse,ssq,rr = 0.0,fs,sy,syy,sig,wt;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Least squares regression. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,4,1))     /* get parameters */
        goto LSRFin;

    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    newline(ctx);
    if (ctx->PMNW != 1) {
        printf1(ctx, "Error: lsreg command requires cross-section data (nw=1).\n");
        goto LSRFin;
    }
    nv = check_nvar(ctx, 0);         /* check variables */
    if (nv == 0)                /* PMNV is number of variables */       
        goto LSRFin;        
       
    prn_nwvar(ctx, 0);               /* print variables */
    nx = nx1 = nv - 1;
    if (ctx->PMNI)  
        printf1(ctx, "Model without intercept.\n");
         
    else {
        ctx->PMNI = 0;
        nx1++;
        if (nx == 0)  
            printf1(ctx, "Model without independent variables.\n");
    }
    if (ctx->DGRPFlg) {              /* process dgrp option */
        if (ctx->PMWVar >= 0) {
            printf1(ctx, "Error: dgrp option not compatible with case weights.\n");
            goto LSRFin;
        }
        if (lsreg_dgrp(ctx, 1))
            goto LSRFin;
    }
    if (ctx->PMS != 0 && ctx->PMWVar >= 0) {
        printf1(ctx, "Error: s=1 option not compatible with case weights.\n");
        goto LSRFin;   
    }

    /* allocate data matrix w[] and solution vector x[] */

    ne = ni = ncon = 0;                 /* number of constraints */
    nw = nx1 + 1;                       /* number of columns */
    mw = ctx->NOC + ctx->NCONSTR + ctx->NDGRP;         /* number of rows */

    if (alloc_acw(ctx, mw * nw + 1))
        goto LSRFin;

    if (alloc_acx(ctx, nx1 + 1))
        goto LSRFin;

    if (alloc_acy(ctx, nx1 + 1))
        goto LSRFin;

    if (alloc_acv(ctx, nx1 * nx1 + 1))   /* used for cov matrix */
        goto LSRFin;

    if (ctx->NCONSTR > 0) {
        if (p_con(ctx, ctx->CmdBuf,nx,ctx->PMNI,mw,nw,ctx->AcW,&ne,&ni))      /* get constraints */
            goto LSRFin;
    }
    if (ctx->NDGRP > 0) {        /* add constraints to data matrix */

        for (i = 0; i < ctx->NDGRP; ++i) {
            k = (ne + i) * nw + 1;
            if (ctx->PMNI == 0)
                k++;

            for (j = 1; j <= nx; ++j)  
                ctx->AcW[k++] = ctx->DGRPW[i][j];
        }
        ne += ctx->NDGRP;
    }
    ncon = ne + ni;

    printf1(ctx, "Equality constraints: %d\n",ne);        
    printf1(ctx, "Inequality constraints: %d\n\n",ni);          

    printf1(ctx, "Reading data. Cases: %d\n",ctx->NOC);
               
    if (ctx->SVEFlg)
        p_warn(ctx, -1,1);

    if (ctx->PMWVar >= 0)  
        printf1(ctx, "Using weights defined by: %s\n",ctx->VName[ctx->PMWVar]);

    /* get data into AcW, also calculate sy (sum of Y) and syy (sum of YY) */

    sy = syy = 0.0;
    k = ne;
    wt = 1.0;
    for (i = 0; i < ctx->NOC; ++i) {                 /* number of cases */

        if (ctx->PMWVar >= 0) {
            wt = get_data(ctx, ctx->PMWVar,i);
            if (wt < 0.0) {
                printf1(ctx, "Error: found negative weight in case %d\n",i + 1);
                goto LSRFin;
            }
            if (wt > 0.0)
                wt = sqrt(wt);
        }

        if (ctx->PMNI == 0 && k >= ne && k < mw - ni)
            ctx->AcW[k * nw + 1] = wt;

        tmp = get_data(ctx, ctx->PMVIdx[0],i) * wt;  
        ctx->AcW[k * nw + nw] = tmp;                         
        sy += tmp;
        syy += tmp * tmp;

        j = 1;           
        if (ctx->PMNI)
            ll = 1;
        else
            ll = 2;
        for (l = 0; l < nx; ++l) {
            ctx->AcW[k * nw + ll] = get_data(ctx, ctx->PMVIdx[j],i) * wt;  
            j++;         
            ll++;
        }
        k++;
    }
    if (ctx->PMF1Def)                               /* write data matrix */
        lsreg_pdat(ctx, mw,nw,ne,ni,ctx->PMNI,nv,0);

    newline(ctx);
               
    /* the covariance matrix is only calculated when: there are no 
       inequality constraints and the number of cases
       is greater than the number of variables */

    nn = mw - ncon;     /* number of least squares cases */

    if (ni > 0 || nn <= nx1) {
        cov = 0;    
        printf1(ctx, "Covariance matrix not calculated.\n");
    }
    else
        cov = 1;

    /*  solve least squares problem */

    ierr = lsei(ctx, ctx->AcW,ne,nn,ni,nx1,ctx->AcX,cov,&rnorme,&rnorml,&ranka,&ranke);

    if (ierr == -1) {
        printf1(ctx, "Equality constraints are contradictory.\n");
        printf1(ctx, "Calculated a least squares solution.\n");
    }
    else if (ierr) {

        if (ierr == -2)  
            printf1(ctx, "Inequality constraints are incompatible.\n");
        else if (ierr == -3)  
            printf1(ctx, "Equality and inequality constraints are contradictory.\n");
        else if (ierr == -4)  
            p_err(ctx, -2,1);
        else  
            printf1(ctx, "LSEI error return: %d\n",ierr);
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
        printf1(ctx, "Rank of "); 
        if (ne > 0)
            printf1(ctx, "reduced ");
        printf1(ctx, "least squares data matrix: %d\n",ranka); 
    }
    if (rdef)
        printf1(ctx, "Warning: data matrix is rank deficient.\n");

    if (ne < nx1)
        printf1(ctx, "Norm of least squares residuals: %lg\n",rnorml); 
#ifdef TDA_R_PACKAGE
        tda_export_mat(ctx, "lsreg.rnorm", &rnorml, 1, 1);
#endif

    if (ctx->PMMPLogDef == 1)                /* create matrix for norm of resid */
        mp_putlog(ctx, rnorml);
        
    if (ctx->PMMPParDef == 1)                /* create matrix for parameters */
        mp_putpar(ctx, nx1,ctx->AcX);

    if (ne > 0) {
        printf1(ctx, "Rank of equality constraints: %d\n",ranke); 
        if (ne != ranke)
            printf1(ctx, "Warning: equality constraints not linear independent.\n");

        printf1(ctx, "Norm of residuals of equality constraints: %lg\n",rnorme);
#ifdef TDA_R_PACKAGE
        tda_export_mat(ctx, "lsreg.rnorme", &rnorme, 1, 1);
#endif
    }
    newline(ctx);

    /*  If residuals are almost zero ignore standard errors */

    if (fabs(rnorme + rnorml) < ctx->EPSI1)  
        cov = 0;
       
    if (cov == 0)
        df = 0;
    else  
        printf1(ctx, "Degrees of freedom: %d\n",df);

    /*  calculate sse = sum of squared residuals, also variance of
        residuals, squared multiple correlation, etc. */

    sig = 0.0;

    if (nn > 0) {       /* do only if there are least squares equations */

        sse  = rnorml * rnorml;
        printf1(ctx, "Sum of squared residuals: %lg\n",sse); 
#ifdef TDA_R_PACKAGE
        tda_export_mat(ctx, "lsreg.sse", &sse, 1, 1);
#endif

        tmp = 0.0;
        if (df > 0) {
            sig = sse / (double)df;
            printf1(ctx, "Variance of residuals: %lg\n",sig);
#ifdef TDA_R_PACKAGE
            tda_export_mat(ctx, "lsreg.sigma2", &sig, 1, 1);
#endif
        }
        ssq = syy - sy * sy / (double)nn;

        if (ssq > 0.0 && ni == 0 && df > 0) {     /* only without inequality constraints */
            rr = (ssq - sse) / ssq;
            if (rr >= 0.0)  
                printf1(ctx, "Squared multiple correlation: %lg\n",rr);
#ifdef TDA_R_PACKAGE
                tda_export_mat(ctx, "lsreg.r2", &rr, 1, 1);
#endif
               
            if (ctx->PMNI == 0 && df > 0) {
                tmp = 1.0 - (double)(nn - 1) * (1.0 - rr) / (double)df;
                if (tmp >= 0.0)  
                    printf1(ctx, "Adjusted: %lg\n",tmp);
#ifdef TDA_R_PACKAGE
                    tda_export_mat(ctx, "lsreg.adj", &tmp, 1, 1);
#endif
            }
        }
        if (ncon == 0 && rdef == 0) {   /* currently only without constraints */

            /*  Calculate Log likelihood */
            /*****************************
            tmp = (double)(-nn) * (L2PI + rlog(ctx, sse / (double)nn) + 1.0) / 2.0;
            printf1(ctx, "\nLog likelihood: %lg\n,tmp);
            *********************/

            /*  Calculate F statistics */

            if (nx1 > 1 && df >= 1 && rr > 0.0 && rr < 1.0) {
                fs = (rr / (double)(nx1 - 1)) / ((1.0 - rr) / (double)df);  
                if (fs > 0.0) {
                    tmp = cdff(ctx, fs,nx1 - 1,df);
                    printf1(ctx, "F-statistic: %lg\n",fs); 
                    printf1(ctx, "Level of significance: %lg\n",tmp); 
#ifdef TDA_R_PACKAGE
                    tda_export_mat(ctx, "lsreg.f", &fs, 1, 1);
                    tda_export_mat(ctx, "lsreg.f.signif", &tmp, 1, 1);
#endif
                }
            }
        }
    }
    else
        cov = df = 0;

    if (cov) {      /* create new cov matrix in AcV, depending on PMS */

        if (ctx->PMResFDef)      /* write residuals */                 
            lsreg_res(ctx, nx,ctx->AcX,cov,nw,ctx->AcW,sig);

        if (lsreg_ncov(ctx, nw,nx1,ctx->AcX,sig))
            goto LSRFin;
   
        if (ctx->PMMPCovDef == 1)              /* create matrix for covariance */
            mp_putcov(ctx, nx1,ctx->AcV);
#ifdef TDA_R_PACKAGE
        export_prn(ctx, "lsreg.vcov", nx1,nx1,nx1,ctx->AcV);
#endif

        if (ctx->PMCovFDef) {                          /* write cov matrix */
            prn_data(ctx, nx1,nx1,nx1,ctx->AcV,ctx->PMCovFd,ctx->PMMFmtS);
            p_wmsg(ctx, 2,ctx->PMCovFName,1);
        }
      
        for (i = 1; i <= nx1; ++i)
            ctx->AcY[i] = ctx->AcV[(i - 1) * nx1 + i];
    }
    prn1_coeff(ctx, nx1,ctx->AcX,ctx->AcY,ctx->PMNI,df,ctx->PMVIdx,1); /* print coefficients */

    if (ctx->NDGRP > 0 && cov && df > 0)         /*  print adj st.dev. of coeff. */
        lsreg_d1grp(ctx, ctx->AcX,ctx->AcV,nx1,ctx->PMNI);

    mp_info(ctx);              /* print info about matrices */

    newline(ctx);
    err = 0;

LSRFin:
    lsreg_dgrp(ctx, 0);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  lsreg_res(nx,b)           Write residuals to PMResFd.                   */

void lsreg_res(TDAContext *ctx, int nx,double *b,int cov,int nw,double *c,double sig)
{
    register int i,j,l;
    int ll,nx1;
    double h,y,ye,se,se1,tmp;

    if (alloc_actmp(ctx, nx + 2))
        return;       

    if (sig < ctx->EPSI1)
        cov = 0;

    if (ctx->PMNI == 0) {
        ctx->AcTmp[1] = 1.0;
        ll = 1;
        nx1 = nx + 1;
    }
    else {
        ll = 0;
        nx1 = nx;
    }
    l = 1;
    fprintf(ctx->PMResFd,"# residuals created by lsreg command.\n");
    fprintf(ctx->PMResFd,"# c%-3d : case number\n",l++);
    l = 2;
    if (ctx->PMNI == 0)
        fprintf(ctx->PMResFd,"# c%-3d : constant one\n",l++);
    j = 1;     
    for (i = 0; i < nx; ++i)   
        fprintf(ctx->PMResFd,"# c%-3d : %s\n",l++,ctx->VName[ctx->PMVIdx[j++]]);
       
    fprintf(ctx->PMResFd,"# c%-3d : %s (dependent)\n",l++,ctx->VName[ctx->PMVIdx[0]]);
    fprintf(ctx->PMResFd,"# c%-3d : predicted\n",l++);
    fprintf(ctx->PMResFd,"# c%-3d : residual\n",l++);
    fprintf(ctx->PMResFd,"# c%-3d : leverage\n",l++);
    fprintf(ctx->PMResFd,"# c%-3d : std. dev. residual\n",l++);
    fprintf(ctx->PMResFd,"# c%-3d : stand. residual\n",l++);
    if (ctx->PMWVar >= 0)
        fprintf(ctx->PMResFd,"# c%-3d : case weights (%s)\n",l++,ctx->VName[ctx->PMWVar]);
    fprintf(ctx->PMResFd,"\n");

    for (i = 0; i < ctx->NOC; ++i) {                 /* number of cases */

        j = 1;           
        for (l = 1; l <= nx; ++l) {
            ctx->AcTmp[ll + l] = get_data(ctx, ctx->PMVIdx[j],i);          
            j++;         
        }
        fprintf(ctx->PMResFd,"%6d ",i + 1);
#ifdef TDA_R_PACKAGE
        /* the residual row, cell by cell: its width depends on how many
           covariates the fit has and whether case weights were given */
        tda_export_cell(ctx, "lsreg.residuals", (double)(i + 1));
#endif
        /* braces are required here: the loop body was a single
           statement, so an unbraced export call landed OUTSIDE the loop
           and read AcTmp[nx1 + 1] -- one past the end -- once per row
           instead of every column in order. */
        for (l = 1; l <= nx1; ++l) {
            rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,ctx->AcTmp[l]);
#ifdef TDA_R_PACKAGE
            tda_export_cell(ctx, "lsreg.residuals", ctx->AcTmp[l]);
#endif
        }

        y = get_data(ctx, ctx->PMVIdx[0],i);      /* dep variable */
        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,y);
#ifdef TDA_R_PACKAGE
        tda_export_cell(ctx, "lsreg.residuals", y);
#endif

        ye = 0.0;
        for (l = 1; l <= nx1; ++l)  
            ye += ctx->AcTmp[l] * b[l];

        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,ye);
        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,y - ye);
#ifdef TDA_R_PACKAGE
        tda_export_cell(ctx, "lsreg.residuals", ye);
        tda_export_cell(ctx, "lsreg.residuals", y - ye);
#endif

        h = se1 = se = -1.0;
        if (cov) {              /* stand dev */
            tmp = 0.0;
            for (l = 1; l <= nx1; ++l) {
                for (j = 1; j <= nx1; ++j)  
                    tmp += c[(l - 1) * nw + j] * ctx->AcTmp[l] * ctx->AcTmp[j];
            }
            if (tmp > 0.0) {
                h = tmp / sig;
                tmp = sig - tmp;
                if (tmp > 0.0)  
                    se = sqrt(tmp);
                se1 = (y - ye) / se;
            }
        }
        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,h);
        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,se);
        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,se1);
#ifdef TDA_R_PACKAGE
        tda_export_cell(ctx, "lsreg.residuals", h);
        tda_export_cell(ctx, "lsreg.residuals", se);
        tda_export_cell(ctx, "lsreg.residuals", se1);
#endif

        if (ctx->PMWVar >= 0) {
            tmp = get_data(ctx, ctx->PMWVar,i);
            rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,tmp);
#ifdef TDA_R_PACKAGE
            tda_export_cell(ctx, "lsreg.residuals", tmp);
#endif
        }
        fprintf(ctx->PMResFd,"\n");
#ifdef TDA_R_PACKAGE
        tda_export_endrow(ctx, "lsreg.residuals");
#endif
    }
    printf1(ctx, "\nResiduals written to: %s\n",ctx->PMResFName);
}

/* ------------------------------------------------------------------------ */
/*  lsreg_pdat(mw,nw,ne,ni,nif,nv,opt)                                      */
/*                                                                          */
/*  Write data to output file. Optionally also dtda description file.       */

void lsreg_pdat(TDAContext *ctx, int mw,int nw,int ne,int ni,int nif,int nv,int opt)
{
    register int i,j,k,l;
    int nn;

    nn = mw - ne - ni;

    prn_data(ctx, nn,nw,nw,ctx->AcW + ne * nw,ctx->PMF1d,ctx->PMFmtS);
    p_wmsg(ctx, 1,ctx->PMF1dName,1);
   
    if (ctx->PMTDAFDef) {
        fprintf(ctx->PMTDAFd,"# data written by lsreg command.\n");
        fprintf(ctx->PMTDAFd,"nvar(\n");
        fprintf(ctx->PMTDAFd,"  dfile = %s,\n",ctx->PMF1dName);
        fprintf(ctx->PMTDAFd,"  noc = %d,\n",nn);
        k = 0;
        if (nif == 0)  
            fprintf(ctx->PMTDAFd,"  INTZ [%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,++k);
           
        l = 1;   
        for (i = 0; i < nv; ++i) {
            if (i == nv - 1)
                j = ctx->PMVIdx[0];
            else
                j = ctx->PMVIdx[opt + l++];

            fprintf(ctx->PMTDAFd,"  %s [%d.%d]",ctx->VName[j],ctx->PMFmt1,ctx->PMFmt2);
            if (ctx->VLabel[j] != NULL)                                                      
                fprintf(ctx->PMTDAFd,"(%s)",ctx->VLabel[j]);                                           
            fprintf(ctx->PMTDAFd," = c%d,\n",++k);
        }
        fprintf(ctx->PMTDAFd,");\n");
    
        printf1(ctx, "TDA description written to: %s\n",ctx->PMTDAFName);
    }
}

/* ------------------------------------------------------------------------ */
/*  lsreg_dgrp(opt)     Process dgrp option. Weights are calculated in      */
/*                      DGRP[i][], i = 0,...,NDGRP - 1. If opt == 0         */
/*                      free previously allocated memory.                   */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int lsreg_dgrp(TDAContext *ctx, int opt)
{
    register int i,j,k,kk;
    int err,fin,nv,nx,nb;
    double tmp,w,sum;
    register char *p,*q,*s;

    if (opt == 0) {
        err = 0;
        goto LSDGRPFin;
    }   
    printf1(ctx, "Processing dgrp option.\n");
    if (ctx->PMRHSTRA == 0) {
        p_err(ctx, -1,1);                 
        return(-1);
    }
    err = -1;
    p = ctx->PMRHSTR;
    fin = 0;
    nx = ctx->PMNV - 1;

    while (*p) {
        if (ctx->NDGRP >= MaxDGRP) {
            printf1(ctx, "Error: exceeded max number of groups.\n");
            goto LSDGRPFin;
        }
        q = skip_nc(ctx, p);
        if (*p != '[' || *(q - 1) != ']') {
            p_err(ctx, -1,1);
            goto LSDGRPFin;
        }   
        if (!*q)
            fin = 1;
        else if (*q != ',') {
            p_err(ctx, -1,1);
            goto LSDGRPFin;
        }
        *q = '\0';

        printf1(ctx, "Group selection: %s\n",p);
        s = p;
        p = get_nvia(ctx, p + 1,&nv,1,&nb);
        if (*p != ']') {
            p_err(ctx, -1,1);
            goto LSDGRPFin;
        }
        if (nb) {
            p_err(ctx, -42,1);
            goto LSDGRPFin;
        }   
        if (nv < 2) {
            if (nv > 0)
                printf1(ctx, "Error: at least two variables from regression varlist must be specified.\n");
            goto LSDGRPFin;
        }
        k = (int)(strlen(s) + 1);
        if (ctx->DGRPDLen < k)
            ctx->DGRPDLen = k;

        if (!(ctx->DGRPDef[ctx->NDGRP] = (char *)calloc((size_t)(k),sizeof(char)))) {
            p_err(ctx, -2,1);
            goto LSDGRPFin;
        }
        memrq(ctx, k,sizeof(char));
        ctx->DGRPDefA[ctx->NDGRP] = (short)(k);       
        strcpy(ctx->DGRPDef[ctx->NDGRP],s);

        if (!(ctx->DGRPVIdx[ctx->NDGRP] = (short *)calloc((size_t)(nv),sizeof(short)))) {
            p_err(ctx, -2,1);
            goto LSDGRPFin;
        }
        memrq(ctx, nv,sizeof(short));
        ctx->DGRPNV[ctx->NDGRP] = (short)(nv);       

        kk = 0;
        for (i = 0; i < nv; ++i) {
            k = ctx->VLVIdx[i];
            for (j = 0; j < i; ++j) {
                if (k == ctx->VLVIdx[j]) {
                    printf1(ctx, "Error: variables in constraints should be unique.\n");
                    goto LSDGRPFin;
                }
            }
            for (j = 1; j < ctx->PMNV; ++j) {
                if (k == ctx->PMVIdx[j]) {
                    ctx->DGRPVIdx[ctx->NDGRP][kk++] = (short)(j);
                    k = -1;
                    break;
                }
            }
            if (k >= 0) {
                printf1(ctx, "Error: at least one variable is not a regressor variable.\n");
                goto LSDGRPFin;
            }
        }
        if (!(ctx->DGRPW[ctx->NDGRP] = (double *)calloc((size_t)(nx + 1),sizeof(double)))) {
            p_err(ctx, -2,1);
            goto LSDGRPFin;
        }
        memrq(ctx, nx + 1,sizeof(double));
        ctx->DGRPWA[ctx->NDGRP] = (short)(nx + 1);

        sum = 0.0;
        for (i = 0; i < ctx->NOC; ++i) {                 /* number of cases */
            for (k = 0; k < nv; ++k) {
                kk = ctx->DGRPVIdx[ctx->NDGRP][k];
                tmp = get_data(ctx, ctx->PMVIdx[kk],i);       
                if (fabs(tmp) < ctx->EPSI1)
                    ;                
                else if (fabs(tmp - 1.0) < ctx->EPSI1) {
                    ctx->DGRPW[ctx->NDGRP][kk] += 1.0;
                    sum += 1.0;
                }
                else {
                    printf1(ctx, "Error: variable %s in case %d isn't a dummy variable.\n",ctx->VName[ctx->PMVIdx[kk]],i + 1);
                    goto LSDGRPFin;
                }
            }
        }
        j = 10; 
        if (j < ctx->VNameLen)
            j = ctx->VNameLen;

        printf1(ctx, "\nVariable  ");
        prnchar(ctx, ' ',j - 10,0);
        printf1(ctx, "   cases       weight\n");
        prnchar(ctx, '-',21 + j,1);
        for (k = 0; k < nv; ++k) {
            kk = ctx->DGRPVIdx[ctx->NDGRP][k];
            i = ctx->PMVIdx[kk];
            printf1(ctx, "%s ",ctx->VName[i]);
            prnchar(ctx, ' ',j - (int)strlen(ctx->VName[i]),0);
            w = tmp = ctx->DGRPW[ctx->NDGRP][kk];
            if (sum > 0.0)
                w /= sum;
            printf1(ctx, "%7d %12.4lf\n",(int)tmp,w);
            ctx->DGRPW[ctx->NDGRP][kk] = w;
        }
        newline(ctx);
        if (sum < ctx->EPSI1) {
            printf1(ctx, "Error: sum of weights is zero.\n");
            goto LSDGRPFin;
        }
        ctx->NDGRP++;

        if (fin)
            break;

        p = q + 1;
    }
    err = 0;

LSDGRPFin:
    if (err || opt == 0) {
        for (i = 0; i < MaxDGRP; ++i) {
            if (ctx->DGRPDefA[i] > 0) {    
                free((char *)ctx->DGRPDef[i]);            
                memrq(ctx, -ctx->DGRPDefA[i],sizeof(char));
                ctx->DGRPDefA[i] = 0;
            }
            if (ctx->DGRPWA[i] > 0) {    
                free((char *)ctx->DGRPW[i]);            
                memrq(ctx, -ctx->DGRPWA[i],sizeof(double));
                ctx->DGRPWA[i] = 0;
            }
            if (ctx->DGRPNV[i] > 0) {    
                free((char *)ctx->DGRPVIdx[i]);            
                memrq(ctx, -ctx->DGRPNV[i],sizeof(short));
                ctx->DGRPNV[i] = 0;
            }
        }
        ctx->NDGRP = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  lsreg_d1grp()       Print adj stand dev. of coefficients.               */
/*                                                                          */

void lsreg_d1grp(TDAContext *ctx, double *beta,double *cov,int nc,int ni)
{
    register int i,j,k;
    int iflag,ib,nv,len;
    double sig;                 

    iflag = 1;
    if (ni)
        iflag = 0;

    len = ctx->DGRPDLen;
    if (len < 8)
        len = 8;

    printf1(ctx, "\nDGroup");
    prnchar(ctx, ' ',len - 6,0);
    printf1(ctx, "weighted adj. stand. dev.\n");
    prnchar(ctx, '-',25 + len,1);
    for (i = 0; i < ctx->NDGRP; ++i) {
        nv = ctx->DGRPNV[i];
        sig = 0.0;
        for (j = 0; j < nv; ++j) {
            k = ctx->DGRPVIdx[i][j];
            ib = k + iflag;
            sig += ctx->DGRPW[i][k] * (beta[ib] * beta[ib] - cov[(ib - 1) * nc + ib]);
        }
        printf1(ctx, "%s",ctx->DGRPDef[i]);
        prnchar(ctx, ' ',10 + len - (int)strlen(ctx->DGRPDef[i]),0);
        if (sig >= 0.0) {
            sig = sqrt(sig);
            rt_printf1_d(ctx, ctx->PMTFmtS,sig);
        }  
        else
            printf1(ctx, "  ---");
        newline(ctx);
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

int lsreg_ncov(TDAContext *ctx, int nw,int nx1,double *b,double sig)
{
    register int i,j,k,l,i1,j1;
    int nxx,ll;
    double tmp,tmp1,tmp2,ye;

    nxx = nx1 * nx1;

    if (ctx->PMS == 0) {
        k = 1;
        for (i = 0; i < nx1; ++i) {
            for (j = 1; j <= nx1; ++j) 
                ctx->AcV[k++] = ctx->AcW[i * nw + j];
        }
        return(0);
    }

    printf1(ctx, "Calculating White's covariance matrix.\n");
    if (sig < ctx->EPSI1) {
        printf1(ctx, "Error: variance of residuals almost zero.\n");
        return(-1);
    }
    sig *= sig;

    if (ctx->PMNI == 0)
        ll = 1;
    else
        ll = 0;

    /* first create AcU[] containing residuals */

    if (alloc_acu(ctx, ctx->NOC + 1))
        return(-1);   

    if (alloc_actmp(ctx, nx1 + 1))
        return(-1);   

    for (i = 0; i < ctx->NOC; ++i) {                 /* number of cases */

        ye = 0.0;
        i1 = 1;
        for (l = 1; l <= nx1; ++l) {
            if (ctx->PMNI == 0 && l == 1)
                ye += b[i1++];
            else  
                ye += b[i1++] * get_data(ctx, ctx->PMVIdx[l - ll],i);
        }
        ctx->AcU[i + 1] = get_data(ctx, ctx->PMVIdx[0],i) - ye;                     
    }

    /* create X' diag(res2) X matrix in AcTmp[] */

    if (alloc_actmp(ctx, nxx + 1))
        return(-1);   

    for (i1 = 1; i1 <= nx1; ++i1) {
        for (j1 = 1; j1 <= nx1; ++j1) {
            tmp = 0.0;
                     
            for (i = 0; i < ctx->NOC; ++i) {                 /* number of cases */

                tmp1 = tmp2 = 1.0;
                if (ctx->PMNI || i1 > 1)
                    tmp1 = get_data(ctx, ctx->PMVIdx[i1 - ll],i);

                if (ctx->PMNI || j1 > 1)
                    tmp2 = get_data(ctx, ctx->PMVIdx[j1 - ll],i);

                tmp += tmp1 * tmp2 * ctx->AcU[i + 1] * ctx->AcU[i + 1];
            }
            ctx->AcTmp[(i1 - 1) * nx1 + j1] = tmp;
        }
    }

    /* create new cov matrix in AcV */

    for (i1 = 1; i1 <= nx1; ++i1) {
        for (j1 = 1; j1 <= nx1; ++j1) {
            tmp = 0.0;
            for (k = 1; k <= nx1; ++k) {
                tmp1 = 0.0;
                for (l = 1; l <= nx1; ++l)  
                    tmp1 += ctx->AcTmp[(k - 1) * nx1 + l] * ctx->AcW[(l - 1) * nw + j1];
                tmp += tmp1 * ctx->AcW[(i1 - 1) * nw + k];
            }
            ctx->AcV[(i1 - 1) * nx1 + j1] = tmp / sig;
        }
    }
    alloc_acu(ctx, 0);
    alloc_actmp(ctx, 0);
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

int zreg(TDAContext *ctx)
{
    register int i,j,k,l;
    int ll,err,nv,nx,nx1,nw,iter;
    int conv,cov,ierr,nn,ranka,ranke,m;
    double rnorme,rnorml,delta,tmp,tmp1;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Regression with censored data. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLP = 1.e-6;

    if (parm(ctx, ctx->CmdBuf + 4,4,1))     /* get parameters */
        goto ZRFin;

    if (ctx->MxIter < 1 || ctx->MxItFlg == 0)
        ctx->MxIter = 20;

    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    newline(ctx);
    if (ctx->PMNW != 1) {
        printf1(ctx, "Error: zreg command requires cross-section data (nw=1).\n");
        goto ZRFin;
    }
    nv = check_nvar(ctx, 0);         /* check variables */
    if (nv == 0)                /* PMNV is number of variables */       
        goto ZRFin;        
       
    prn_nwvar(ctx, 0);               /* print variables */

    if (ctx->PMYWVar < 0) {
        printf1(ctx, "Error: need censoring information.\n");
        goto ZRFin;
    }
    printf1(ctx, "Number of cases: %d\n",ctx->NOC);
    printf1(ctx, "Censoring indicator: %s\n",ctx->VName[ctx->PMYWVar]);
    printf1(ctx, "Max number of iterations: %d\n",ctx->MxIter);
    printf1(ctx, "Tolerance for convergence: %lg\n\n",ctx->TOLG);

    nx = nx1 = nv - 1;
    if (ctx->PMNI)  
        printf1(ctx, "Model without intercept.\n");
         
    else {
        ctx->PMNI = 0;
        nx1++;
        if (nx == 0)  
            printf1(ctx, "Model without independent variables.\n");
    }

    /* allocate data matrix w[] and solution vector x[] */

    nw = nx1 + 1;                   /* number of columns */
    nn = ctx->NOC;                       /* number of rows */

    if (alloc_acw(ctx, nn * nw + 1))     /* least squares matrix */
        goto ZRFin;
    if (alloc_acx(ctx, nx1 + 1))         /* used for parameters */
        goto ZRFin;
    if (alloc_acy1(ctx, nx1 + 1))        /* used for previous parameters */
        goto ZRFin;
    if (alloc_acy(ctx, nn + 1))          /* used for dependent variable */
        goto ZRFin;
    if (alloc_acz(ctx, nn + 1))          /* used for residuals */
        goto ZRFin;
    if (alloc_acns(ctx, nn + 1))         /* used for censoring indicator */
        goto ZRFin;
    if (alloc_acn(ctx, nn + 1))          /* used for sorting pointer */
        goto ZRFin;
    if (alloc_acu(ctx, nn + 1))          /* used for jumps of K-M */
        goto ZRFin;

    for (i = 0; i < ctx->NOC; ++i) {
        ctx->AcY[i] = get_data(ctx, ctx->PMVIdx[0],i);  
        if ((int)get_data(ctx, ctx->PMYWVar,i) != 0)  
            ctx->AcNS[i] = 1;
    }

    /* begin iterations */

    if (ctx->SILENTFlg < 2)
        printfe(ctx, "\n  Iter   Norm of Residuals     Par Change\n");

    if (ctx->PMProtFDef)  
        fprintf(ctx->PMProtFd,"Iteration  Parameter vector\n");

    cov = 0;    /* do not calculate covariance matrix */
    conv = 0;
    delta = 0.0;

    for (iter = 0; iter <= ctx->MxIter; ++iter) {
 
        /* get data into AcW  */
       
        k = 0;
        for (i = 0; i < ctx->NOC; ++i) {         
       
            if (ctx->PMNI == 0)
                ctx->AcW[k * nw + 1] = 1.0;

            ctx->AcW[k * nw + nw] = ctx->AcY[i];                      

            j = 1;           
            if (ctx->PMNI)
                ll = 1;
            else
                ll = 2;
            for (l = 0; l < nx; ++l) {
                ctx->AcW[k * nw + ll] = get_data(ctx, ctx->PMVIdx[j],i);  
                j++;         
                ll++;
            }
            k++;
        }

        /*  solve least squares problem */

        ierr = lsei(ctx, ctx->AcW,0,nn,0,nx1,ctx->AcX,cov,&rnorme,&rnorml,&ranka,&ranke);

        if (ierr) {
            if (ierr == -4)  
                p_err(ctx, -2,1);
            else  
                printf1(ctx, "LSEI error return: %d\n",ierr);
            goto ZRFin;  
        }
        if (ranka != nx1 - ranke) {    /* set if rank deficient */
            printf1(ctx, "Rank of least squares matrix: %d\n",ranka);
            printf1(ctx, "Will not continue.\n");
            goto ZRFin;
        }
        if (iter == 0) {
            printf1(ctx, "Regression results with original values.\n");
            printf1(ctx, "Rank of least squares data matrix: %d\n",ranka); 
            printf1(ctx, "Norm of least squares residuals: %lg\n",rnorml); 
#ifdef TDA_R_PACKAGE
            tda_export_mat(ctx, "lsreg.rnorm", &rnorml, 1, 1);
#endif

            prn1_coeff(ctx, nx1,ctx->AcX,ctx->AcX,ctx->PMNI,0,ctx->PMVIdx,1); /* print coefficients */
            newline(ctx);
        }

        /* print parameters to protocol file */
   
        if (ctx->PMProtFDef) {
            fprintf(ctx->PMProtFd,"%7d    ",iter);
            for (j = 1; j <= nx1; ++j)
                rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,ctx->AcX[j]);
            fprintf(ctx->PMProtFd,"\n");
        }
   
        /* check convergence */

        if (iter > 0) {
            delta = 0.0;
            for (j = 1; j <= nx1; ++j) {
                tmp = fabs(ctx->AcX[j] - ctx->AcY1[j]);
                if ((ctx->AcX[j]) != 0.0)
                    tmp /= fabs(ctx->AcX[j]);
                delta = dmax(ctx, delta,tmp);
            }
        }
        if (ctx->SILENTFlg < 2)  
            printfe(ctx, "  %3d  %20.13e %17.10e\n",iter,rnorml,delta);

        if (iter > 0 && delta <= ctx->TOLP) {
            conv = 1;
            break;
        }
        for (j = 1; j <= nx1; ++j)  
            ctx->AcY1[j] = ctx->AcX[j];           

        /* update dependent variable. first calculate residuals in AcZ[] */

/**   
tda_out("VORHER\n");
        for (i = 0; i < NOC; ++i)  
tda_out("i=%d d=%d y=%lg\n",i,AcNS[i], AcY[i]);
newline(ctx);
**/        


        for (i = 0; i < ctx->NOC; ++i) {
            tmp = 0.0;
            j = 1;
            if (ctx->PMNI == 0) {
                tmp = ctx->AcX[j];
                j++;   
            }
            k = 1;
            while (j <= nx1) {
                tmp += ctx->AcX[j] * get_data(ctx, ctx->PMVIdx[k],i);  
                j++;
                k++;
            }
            ctx->AcZ[i] = get_data(ctx, ctx->PMVIdx[0],i) - tmp;
            if (ctx->AcNS[i] == 0)
                ctx->AcY[i] = tmp;

        }
        if (sortdp2a(ctx, ctx->NOC,ctx->AcZ,ctx->AcNS,ctx->AcN,0))      /* sort */
            goto ZRFin;    

/**      
tda_out("Erster Schritt\n");
        for (i = 0; i < NOC; ++i)  
tda_out("i=%d d=%d y=%lg\n",i,AcNS[i], AcY[i]);
newline(ctx);
**/    
        for (i = 0; i < ctx->NOC; ++i)  
            ctx->AcU[i] = 0.0;

        tmp = 1.0 / (double)ctx->NOC;
        m = ctx->NOC - 1;
        for (i = 0; i < ctx->NOC; ++i) {
            j = ctx->AcN[i];
            ctx->AcU[j] += tmp;
            if (ctx->AcNS[j] == 0 && i < ctx->NOC - 1) {
                tmp1 = ctx->AcU[j] / (double)m;
                for (k = i + 1; k < ctx->NOC; ++k)  
                    ctx->AcU[ctx->AcN[k]] += tmp1;
                ctx->AcU[j] = 0.0;
            }
            m--;
        }
        tmp = 0.0;
        for (i = 0; i < ctx->NOC; ++i) {
            k = ctx->AcN[i];
            tmp += ctx->AcU[k];
/*      
tda_out("i=%d AcNS[i]=%d z=%10.4lf  F=%10.4lf  U=%10.4lf\n",i,AcNS[i], AcZ[i],tmp,AcU[i]);
*/     
        }

        for (i = 0; i < ctx->NOC; ++i) {
            j = ctx->AcN[i];
            if (ctx->AcNS[j] == 0) {         /* if censored */
                if (i == ctx->NOC - 1) {
                    ctx->AcY[j] += ctx->AcZ[j];
                }
                else {
                    tmp = tmp1 = 0.0;
                    for (k = i + 1; k < ctx->NOC; ++k) {
                        l = ctx->AcN[k];
/*                      tda_out(".... l=%d d=%d z=%lg\n",l,AcNS[l],AcZ[l]); */
                        if (ctx->AcNS[l] != 0 || k == ctx->NOC - 1) {
                            tmp += ctx->AcZ[l] * ctx->AcU[l];
                            tmp1 += ctx->AcU[l];
                        }
                    }
                    if (tmp1 != 0.0)
                        ctx->AcY[j] += tmp / tmp1;
                }
            }
        }
/*    
tda_out("ENDE\n");
        for (i = 0; i < NOC; ++i)  
tda_out("i=%d d=%d y=%lg\n",i,AcNS[i], AcY[i]);
newline(ctx);
*/       

    }
    printf1(ctx, "Convergence ");   
    if (conv == 0)
        printf1(ctx, "not ");
    printf1(ctx, "reached in %d iterations.\n",iter);

    printf1(ctx, "Rank of least squares data matrix: %d\n",ranka); 
    printf1(ctx, "Norm of least squares residuals: %lg\n",rnorml); 
#ifdef TDA_R_PACKAGE
    tda_export_mat(ctx, "lsreg.rnorm", &rnorml, 1, 1);
#endif

    prn1_coeff(ctx, nx1,ctx->AcX,ctx->AcX,ctx->PMNI,0,ctx->PMVIdx,1); /* print coefficients */
    newline(ctx);

    if (ctx->PMF1Def) {
        for (i = 0; i < ctx->NOC; ++i) {
            tmp = get_data(ctx, ctx->PMVIdx[0],i);
            fprintf(ctx->PMF1d,"%6d ",i + 1);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,tmp);  
            fprintf(ctx->PMF1d,"%2d ",ctx->AcNS[i]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[i]);
            fprintf(ctx->PMF1d,"\n");        
        }
        printf1(ctx, "%d records written to: %s\n",ctx->NOC,ctx->PMF1dName);
    }
    err = 0;

ZRFin:
    p_clean(ctx);
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

int zreg1(TDAContext *ctx)
{
    register int i,j,k,l;
    int err,nv,nx,nxu,nw,nxv,nk,iter,m,n,nn;
    int ierr,ranka,ranke,ncen;
    double rnorme,rnorml,delta,tmp,tmp1,time;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Regression with censored data. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLP = 1.e-6;

    if (parm(ctx, ctx->CmdBuf + 5,4,1))     /* get parameters */
        goto ZR1Fin;

    if (ctx->MxIter < 1 || ctx->MxItFlg == 0)
        ctx->MxIter = 20;

    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    newline(ctx);
    if (ctx->PMNW != 1) {
        printf1(ctx, "Error: zreg command requires cross-section data (nw=1).\n");
        goto ZR1Fin;
    }
    nv = check_nvar(ctx, 0);         /* check variables */
    if (nv == 0)                /* PMNV is number of variables */       
        goto ZR1Fin;        
       
    prn_nwvar(ctx, 0);               /* print variables */

    if (ctx->PM1NV > 0) {
        printf1(ctx, "Plus time-varying indicator variables: %s",ctx->VName[ctx->PM1VIdx[0]]);
        for (i = 1; i < ctx->PM1NV; ++i)
            printf1(ctx, ", %s",ctx->VName[ctx->PM1VIdx[i]]);
        newline(ctx);
        newline(ctx);
    }

    if (ctx->PMYWVar < 0) {
        printf1(ctx, "Error: need censoring information.\n");
        goto ZR1Fin;
    }
/***###*************/
#if 0
    /* Debug scaffolding left in the 6.4 sources: this block dumps, for
       every case, a 429-column risk/indicator matrix built from the
       first covariate and the response, and then jumps straight to
       ZR1Fin -- so as shipped, zreg1 printed this dump and never
       reached the estimator below at all.  Disabled (not deleted) for
       the 6.4q build so the complete Buckley-James residual-life
       estimator that follows runs; see doc/changes-from-tda.md. */
    for (i = 0; i < ctx->NOC; ++i) {
        k = (int)get_data(ctx, ctx->PMVIdx[0],i);        
        l = (int)get_data(ctx, ctx->PMVIdx[1],i);        
        tda_out("%4d %d %6d %6d",i + 1,
                (int)get_data(ctx, ctx->PMYWVar,i),k,l);   

        for (j = 0; j <= 428; ++j) {
            if (l > j || j >= k    )
                tda_out(" 0");
            else
                tda_out(" 1");
        }
        newline(ctx);
    }
    goto ZR1Fin;
#endif
/********************/








    nxu = nv - 1;
    nxv = ctx->PM1NV;
    nx = 1 + nxu + nxv;

    printf1(ctx, "Number of cases: %d\n",ctx->NOC);
    printf1(ctx, "Number of time-independent covariates: %d\n",nxu);
    printf1(ctx, "Number of time-varying covariates: %d\n",nxv);
    printf1(ctx, "Censoring indicator: %s\n",ctx->VName[ctx->PMYWVar]);
    printf1(ctx, "Max number of iterations: %d\n",ctx->MxIter);
    printf1(ctx, "Tolerance for convergence: %lg\n\n",ctx->TOLG);

    /* allocate data matrix w[] and solution vector x[] */

    nw = nx + 1;                    /* number of columns */

    if (alloc_acw(ctx, ctx->NOC * nw + 1))    /* least squares matrix */
        goto ZR1Fin;
    if (alloc_acx(ctx, nx + 1))          /* used for parameters */
        goto ZR1Fin;
    if (alloc_acy1(ctx, nx + 1))         /* used for previous parameters */
        goto ZR1Fin;
    if (alloc_acy(ctx, ctx->NOC + 1))         /* used for dependent variable */
        goto ZR1Fin;
    if (alloc_acz(ctx, ctx->NOC + 1))         /* used for residuals */
        goto ZR1Fin;
    if (alloc_acns(ctx, ctx->NOC + 1))        /* used for censoring indicator */
        goto ZR1Fin;
    if (alloc_acms(ctx, ctx->NOC + 1))        /* used for censoring indicator */
        goto ZR1Fin;
    if (alloc_acn(ctx, ctx->NOC + 1))         /* used for sorting pointer */
        goto ZR1Fin;
    if (alloc_acj(ctx, ctx->NOC + 1))         /* used for sorting pointer */
        goto ZR1Fin;
    if (alloc_acm(ctx, nxv + 1))         /* used for counting covariate values */
        goto ZR1Fin;
    if (alloc_acu(ctx, ctx->NOC + 1))         /* used for jumps of K-M */
        goto ZR1Fin;

    ncen = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        ctx->AcY[i] = get_data(ctx, ctx->PMVIdx[0],i);  
        if ((int)get_data(ctx, ctx->PMYWVar,i) != 0)  
            ctx->AcNS[i] = 1;
        else
            ncen++;

        if (ctx->AcY[i] <= 0.0) {
            printf1(ctx, "Error: found non-positive value in record %d\n",i + 1);
            goto ZR1Fin;
        }
    }

    if (sortdp(ctx, ctx->NOC,ctx->AcY,ctx->AcN))      /* sort */
        goto ZR1Fin;   

    printf1(ctx, "Time      N    Cen ");
    l = 42;
    for (j = 0; j < nxv; ++j) {              
        k = ctx->PM1VIdx[j];
        prnchar(ctx, ' ',(int)(6 - strlen(ctx->VName[k])),0);
        printf1(ctx, "%s ",ctx->VName[k]);
        l = (int)(l + (imax(ctx, 7,(int)(strlen(ctx->VName[k])))));
    }
    printf1(ctx, " Iter  Model parameters\n");
    prnchar(ctx, '-',l,1);

    time = 0.0;
    nk = 0;
    n = ctx->NOC;
    while (n > 10) {
        while (nk < ctx->NOC) {
            j = ctx->AcN[nk];
            if (ctx->AcY[j] > time)
                break;
            if (ctx->AcNS[j] == 0)
                ncen--;
            n--;    
            nk++;
        }       
        printf1(ctx, "%4d %6d %6d ",(int)time,n,ncen);


        /* create least squares matrix */

        for (l = 0; l < nxv; ++l)
            ctx->AcM[l] = 0;

        nn = 0;
        for (i = nk; i < ctx->NOC; ++i) {
            j = ctx->AcN[i];
            ctx->AcW[nn * nw + nw] = ctx->AcY[j] - time;               
            ctx->AcW[nn * nw + 1] = 1.0;
       
            for (l = 1; l <= nxu; ++l)  
                ctx->AcW[nn * nw + 1 + l] = get_data(ctx, ctx->PMVIdx[l],j);  
     
            for (l = 0; l < nxv; ++l) {
                tmp = get_data(ctx, ctx->PM1VIdx[l],j);          
                if (tmp <= time) {
                    ctx->AcM[l] += 1;
                    tmp = 1.0;
                }
                else
                    tmp = 0.0;

                ctx->AcW[nn * nw + 2 + nxu + l] = tmp;                       
            }
            nn++;
        }
        if (n != nn) {
            printf1(ctx, "\nFatal error: n=%d nn=%d\n",n,nn);
            goto ZR1Fin;
        }
        for (l = 0; l < nxv; ++l)
            printf1(ctx, "%6.4lf ",(double)ctx->AcM[l]/(double)n);
/**
tda_out("lsmatrix\n");
for (i = 0; i < nn; ++i) {
    for (j = 1; j <= nw; ++j)
        tda_out("%lf ",AcW[i * nw + j]);
    newline(ctx);
}
**/
        delta = 0.0;

        for (iter = 0; iter <= ctx->MxIter; ++iter) {
/**
tda_out("\niter=%d\n",iter);
tda_out("lsmatrix\n");
for (i = 0; i < nn; ++i) {
    for (j = 1; j <= nw; ++j)
        tda_out("%lf ",AcW[i * nw + j]);
    newline(ctx);
}
**/
            ierr = lsei(ctx, ctx->AcW,0,nn,0,nx,ctx->AcX,0,&rnorme,&rnorml,&ranka,&ranke);
            if (ierr) {
                newline(ctx);
                if (ierr == -4)  
                    p_err(ctx, -2,1);
                else  
                    printf1(ctx, "LSEI error return: %d\n",ierr);
                goto ZR1Fin;  
            }
            if (ranka != nx - ranke) {    /* set if rank deficient */
                printf1(ctx, "Rank of least squares matrix: %d\n",ranka);
                printf1(ctx, "Will not continue.\n");
                goto ZR1Fin;
            }

            /* check convergence */

            if (iter > 0) {
                delta = 0.0;
                for (j = 1; j <= nx; ++j) {
                    tmp = fabs(ctx->AcX[j] - ctx->AcY1[j]);
                    if ((ctx->AcX[j]) != 0.0)
                        tmp /= fabs(ctx->AcX[j]);
                    delta = dmax(ctx, delta,tmp);
                }
                if (delta <= ctx->TOLP) {
                    break;
                }
            }
            for (j = 1; j <= nx; ++j)  
                ctx->AcY1[j] = ctx->AcX[j];           

            /* update dependent variable  */
            /* first recreate least squares matrix and calculate residuals
               in AcZ[] */

            nn = 0;
            for (i = nk; i < ctx->NOC; ++i) {
                j = ctx->AcN[i];
                ctx->AcW[nn * nw + nw] = ctx->AcY[j] - time;               
                ctx->AcW[nn * nw + 1] = 1.0;
                tmp1 = ctx->AcX[1];

                for (l = 1; l <= nxu; ++l) {
                    tmp  = get_data(ctx, ctx->PMVIdx[l],j);  
                    ctx->AcW[nn * nw + 1 + l] = tmp;                     
                    tmp1 += tmp * ctx->AcX[l + 1];
                }

                for (l = 0; l < nxv; ++l) {
                    tmp = get_data(ctx, ctx->PM1VIdx[l],j);          
                    if (tmp <= time)   
                        tmp = 1.0;
                    else
                        tmp = 0.0;
                    ctx->AcW[nn * nw + 2 + nxu + l] = tmp;                       
                    tmp1 += tmp * ctx->AcX[2 + nxu + l];
                }
                ctx->AcZ[nn] = ctx->AcY[j] - time - tmp1;
                ctx->AcMS[nn] = ctx->AcNS[j];

                if (ctx->AcNS[j] == 0)
                    ctx->AcW[nn * nw + nw] = tmp1;                      

                nn++;
            }
/**
tda_out("##  2-lsmatrix\n");
for (i = 0; i < nn; ++i) {
    for (j = 1; j <= nw; ++j)
        tda_out("%lf ",AcW[i * nw + j]);
tda_out("z=%lg AcMS=%d ",AcZ[i],AcMS[i]   );
    newline(ctx);
}

**/





            if (sortdp2a(ctx, nn,ctx->AcZ,ctx->AcMS,ctx->AcJ,0))      /* sort */
                goto ZR1Fin;   
       
            for (i = 0; i < nn; ++i)  
                ctx->AcU[i] = 0.0;

            tmp = 1.0 / (double)nn;
            m = nn - 1;
            for (i = 0; i < nn; ++i) {
                j = ctx->AcJ[i];
                ctx->AcU[j] += tmp;
                if (ctx->AcMS[j] == 0 && i < nn - 1) {
                    tmp1 = ctx->AcU[j] / (double)m;
                    for (k = i + 1; k < nn; ++k)  
                        ctx->AcU[ctx->AcJ[k]] += tmp1;
                    ctx->AcU[j] = 0.0;
                }
                m--;
            }
            tmp = 0.0;
            for (i = 0; i < nn; ++i) {
                k = ctx->AcJ[i];
                tmp += ctx->AcU[k];
            }
            for (i = 0; i < nn; ++i) {
                j = ctx->AcJ[i];
                if (ctx->AcMS[j] == 0) {         /* if censored */
                    if (i == nn - 1) {
                        ctx->AcW[j * nw + nw] += ctx->AcZ[j];                      
                    }
                    else {
                        tmp = tmp1 = 0.0;
                        for (k = i + 1; k < nn; ++k) {
                            l = ctx->AcJ[k];
                            if (ctx->AcMS[l] != 0 || k == nn - 1) {
                                tmp += ctx->AcZ[l] * ctx->AcU[l];
                                tmp1 += ctx->AcU[l];
                            }
                        }
                        if (tmp1 != 0.0) {
                            ctx->AcW[j * nw + nw] += tmp / tmp1;                  
                        }
                    }
                }
            }



        } /* end of iteration */


/**

newline(ctx);
    printf1(ctx, "Convergence ");   
    if (conv == 0)
        printf1(ctx, "not ");
    printf1(ctx, "reached in %d iterations.\n",iter);

    printf1(ctx, "Rank of least squares data matrix: %d\n",ranka); 
    printf1(ctx, "Norm of least squares residuals: %lg\n",rnorml); 
#ifdef TDA_R_PACKAGE
    tda_export_mat(ctx, "lsreg.rnorm", &rnorml, 1, 1);
#endif

    prn1_coeff(ctx, nx,AcX,AcX,PMNI,0,PMVIdx,1); 
    newline(ctx);
**/

    printf1(ctx, "%5d ",iter);
    for (j = 1; j <= nx; ++j)
        rt_printf1_d(ctx, ctx->PMTFmtS,ctx->AcX[j]);
    newline(ctx);
#ifdef TDA_R_PACKAGE
    /* one row of the residual-life trajectory, exactly as printed:
       time, cases at risk, censored among them, each date's share
       (already written by the loop above via acm), the BJ iteration
       count, then the parameter vector (CONTRIBUTING.md).  The
       shares are recomputed here from the same acm counters the print
       used; nxv is small. */
    {
        double *erow = (double *)malloc((size_t)(3 + nxv + 1 + nx) *
                                        sizeof(double));
        if (erow != NULL) {
            int ec = 0;
            erow[ec++] = time;
            erow[ec++] = (double)n;
            erow[ec++] = (double)ncen;
            for (j = 0; j < nxv; ++j)
                erow[ec++] = (double)ctx->AcM[j] / (double)n;
            erow[ec++] = (double)iter;
            for (j = 1; j <= nx; ++j)
                erow[ec++] = ctx->AcX[j];
            tda_export_row(ctx, "zreg1.path", erow, ec);
            free(erow);
        }
    }
#endif
/*
    if (time >= 12.0)
        break;

*/
        time += 1.0;
    }

    err = 0;

ZR1Fin:
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "zreg1.path");
#endif
    p_clean(ctx);
    return(err);
}






