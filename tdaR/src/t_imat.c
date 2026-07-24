/****************************************************************************/
/*  t_imat                                                                  */
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
#include "t_gdat.h"
#include "t_gf.h"
#include "t_alloc.h"
#include "t_mat.h"
#include "t_matc.h"
#include "t_matf.h"
#include "t_sort.h"
#include "t_eval.h"  
#include "t_eval1.h"  
#include "t_eval3.h"  
#include "t_ml.h"  
#include "t_gmin.h"  
#include "t_gdd.h"  
#include "t_svd.h"  
#include "t_lp.h"  
#include "t_lsei.h"  
#include "t_com.h"  
#include "tda_context.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_imat                                                     */

int m_midf(TDAContext *ctx, char *cmd);
int m_midf1(TDAContext *ctx, char *cmd);
int m_midf2(TDAContext *ctx, char *cmd);
int get_ipart(TDAContext *ctx, int n,double *xl,double *xu,double *ip);
int m_midf3(TDAContext *ctx, char *cmd);
int gmin(TDAContext *ctx);
int range(TDAContext *ctx);
int alloc_par(TDAContext *ctx, int n,int opt);
int igmin(TDAContext *ctx, int opt,int narg,int nbmax,double *par,double *lb,double *ub,int mxit, double tolbw,double tolfd,double tolfe,int gc,int typ);
void igmin_cpar(TDAContext *ctx, int n,double *x,double *par);
int alloc_list(TDAContext *ctx, int n,int m);
int igmin_res(TDAContext *ctx, double tolfe,int opt);
void i_add(TDAContext *ctx, double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_sub(TDAContext *ctx, double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_mul(TDAContext *ctx, double xl,double xh,double yl,double yh,double *rl,double *rh);
int  i_div(TDAContext *ctx, double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_abs(TDAContext *ctx, double xl,double xh,double *rl,double *rh);
void i_max(TDAContext *ctx, double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_min(TDAContext *ctx, double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_square(TDAContext *ctx, double xl,double xh,double *rl,double *rh);
void i_sqrt(TDAContext *ctx, double xl,double xh,double *rl,double *rh);

int igmin_fun(TDAContext *ctx, int opt,int typ,int n,double *xl,double *xh, double *rl,double *rh,int deriv,double *gl,double *gh);
int idf(TDAContext *ctx);
int sddf(TDAContext *ctx);
int iddf(TDAContext *ctx);
void prn_ddf(TDAContext *ctx, int al,int m,int *low,int *high,double *mdf,int noc);
void prn_ddf1(TDAContext *ctx, int al,int m,double *df);
void prn_ddf2(TDAContext *ctx, int al,int m,int *low,int *high,double *mdf,double *cdf,int noc); 
int imean(TDAContext *ctx);
int f_range(TDAContext *ctx, int fn,int na,int nb,int mxit,double tolbw,double tolfd, double tolfe,int il,int ih,double *fminp,double *fmaxp);
int ivar(TDAContext *ctx);
int icov2(TDAContext *ctx, int typ);
int igini(TDAContext *ctx);
int ivar1(TDAContext *ctx);
double ivarf(TDAContext *ctx, int n,double x,double *xl,double *xh);


/*  Global parameters and variables                                         */





/*--------------------------------------------------------------------------*/
/*  m_midf(cmd)     midf(XL,XU,DL,DU,DM)        matrix command              */
/*                                                                          */
/*                  Calculates distribution function DM, and lower (DL)     */
/*                  and upper (DU) bounds for an interval vector (XL,XU).   */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_midf(TDAContext *ctx, char *cmd)
{
    register int i,j;
    int err,n,nn,row,col,idx,idx1,idx2,idx3,ivflg;
    register char *p;
    double xl,xu,xx,tmp,f,fl,fu;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    if ((p = n_mexpr(ctx, cmd + 5,2,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MIDFFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MIDFFin;
    if (col != 1) {
        mat_err(ctx, 7);
        goto MIDFFin;
    }
    if (alloc_actmp(ctx, 2 * row))
        goto MIDFFin;

    n = 0;
    for (i = 1; i <= row; ++i) {
        xl = ctx->MX[0][i];
        xu = ctx->MX[1][i];
        if (fabs(xl - xu) <= ctx->EPSI1) {
            mat_err(ctx, 11);
            goto MIDFFin;
        }
        if (xl > xu) {
            tmp = xl;
            xl = xu;
            xu = tmp;
            ctx->MX[0][i] = xl;
            ctx->MX[1][i] = xu;
        }
        ctx->AcTmp[n++] = xl;
        ctx->AcTmp[n++] = xu;
    }
    if (sortd(ctx, n,ctx->AcTmp,0))
        goto MIDFFin;

    nn = 1;
    for (i = 1; i < n; ++i) {
        if (ctx->AcTmp[i] > ctx->AcTmp[i - 1])
            ctx->AcTmp[nn++] = ctx->AcTmp[i];
    }
    if ((p = m_getmat(ctx, p,nn,1,&idx,cmd,0)) == NULL)
        goto MIDFFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MIDFFin;
    if ((p = m_getmat(ctx, p,nn,1,&idx1,cmd,0)) == NULL)
        goto MIDFFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MIDFFin;
    if ((p = m_getmat(ctx, p,nn,1,&idx2,cmd,0)) == NULL)
        goto MIDFFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MIDFFin;
    if ((p = m_getmat(ctx, p,nn,1,&idx3,cmd,1)) == NULL)
        goto MIDFFin;

    /* save values of induced partition */

    for (i = 1; i <= nn; ++i)
        ctx->MatVal[idx][i] = ctx->AcTmp[i - 1];

    if (sortd2(ctx, row,ctx->MX[0] + 1,ctx->MX[1] + 1))
        goto MIDFFin;

    /* calculate distribution functions */

    for (i = 1; i <= nn; ++i) {
        tmp = ctx->AcTmp[i - 1];
        fl = fu = f = 0.0;
        for (j = 1; j <= row; ++j) {
            xl = ctx->MX[0][j];
            xu = ctx->MX[1][j];
            if (xu <= tmp)
                fl += 1.0;
            if (xl <= tmp) {
                fu += 1.0;
                if (xl < tmp) {
                    xx = dmin(ctx, xu,tmp);
                    f += (xx - xl) / (xu - xl);
                }
            }
            else
                break;
        }
        ctx->MatVal[idx1][i] = fl / (double)row;
        ctx->MatVal[idx2][i] = fu / (double)row;
        ctx->MatVal[idx3][i] = f  / (double)row;
    }
    err = 0;

MIDFFin:
    alloc_actmp(ctx, 0);
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_midf1(cmd)    midf(XL,XU,F)           matrix command                  */
/*                                                                          */
/*                  Calculates mean distribution function F for lower end   */
/*                  points of [XL,XU].                                      */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_midf1(TDAContext *ctx, char *cmd)
{
    register int i,j,jj;
    int err,row,col,idx,ivflg;
    register char *p;
    double xl,xu,xx,tmp,f;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    if ((p = n_mexpr(ctx, cmd + 6,2,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MIDF1Fin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MIDF1Fin;
    if (col != 1) {
        mat_err(ctx, 7);
        goto MIDF1Fin;
    }
    if ((p = m_getmat(ctx, p,row,1,&idx,cmd,1)) == NULL)
        goto MIDF1Fin;

    if (alloc_acn(ctx, row))
        goto MIDF1Fin;

    if (sortdp(ctx, row,ctx->MX[0] + 1,ctx->AcN))
        goto MIDF1Fin;

    for (i = 1; i <= row; ++i) {
        tmp = ctx->MX[0][i];
        f = 0.0;
        for (j = 0; j < row; ++j) {
            jj = ctx->AcN[j] + 1;
            xl = ctx->MX[0][jj];
            xu = ctx->MX[1][jj];
            if (xl < tmp) {
                xx = dmin(ctx, xu,tmp);
                f += (xx - xl) / (xu - xl);
            }
            else
                break;
        }
        ctx->MatVal[idx][i] = f  / (double)row;
    }
    err = 0;

MIDF1Fin:
    alloc_acn(ctx, 0);
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_midf2(cmd)    midf(XL,XU,F)                 matrix command            */
/*                                                                          */
/*                  Calculates mean distribution function F for lower end   */
/*                  points of [XL,XU].                                      */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_midf2(TDAContext *ctx, char *cmd)
{
    register int i,j,jj;
    int err,nn,row,col,idx,ivflg;
    register char *p;
    double xl,xu,xx,tmp,f;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    if ((p = n_mexpr(ctx, cmd + 6,2,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MIDF2Fin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MIDF2Fin;
    if (col != 1) {
        mat_err(ctx, 7);
        goto MIDF2Fin;
    }
    if ((p = m_getmat(ctx, p,row,1,&idx,cmd,1)) == NULL)
        goto MIDF2Fin;

    if (alloc_actmp(ctx, 2 * row))
        goto MIDF2Fin;

    if ((nn = get_ipart(ctx, row,ctx->MX[0],ctx->MX[1],ctx->AcTmp)) < 2)
        goto MIDF2Fin;
   
    if (alloc_acn(ctx, row))
        goto MIDF2Fin;

    if (sortdp(ctx, row,ctx->MX[0] + 1,ctx->AcN))
        goto MIDF2Fin;
     
    if (alloc_acu(ctx, nn))          /* used for distr function */
        goto MIDF2Fin;

    for (i = 0; i < nn; ++i) {
        tmp = ctx->AcTmp[i];              
        f = 0.0;
        for (j = 0; j < row; ++j) {
            jj = ctx->AcN[j] + 1;
            xl = ctx->MX[0][jj];
            xu = ctx->MX[1][jj];
            if (xl < tmp) {
                xx = dmin(ctx, xu,tmp);
                f += (xx - xl) / (xu - xl);
            }
            else
                break;
        }
        ctx->AcU[i] = f  / (double)row;
    }
    for (i = 1; i <= row; ++i) {
        tmp = ctx->MX[0][i];
        for (j = 0; j < nn; ++j) {
            if (fabs(tmp - ctx->AcTmp[j]) < ctx->EPSI1)
                break;
        }
        tmp = 0.0;
        for (jj = j; jj < nn - 1; ++jj) {
            xl = ctx->AcTmp[jj];
            xu = ctx->AcTmp[jj + 1];
            f  = (ctx->AcU[jj + 1] - ctx->AcU[jj]) / (xu - xl);
            tmp += f * (xu * xu - xl * xl);
        }
        tmp /= 2.0;
        if (1.0 - ctx->AcU[j] <= ctx->EPSI) 
            goto MIDF2Fin;

        ctx->MatVal[idx][i] = tmp / (1.0 - ctx->AcU[j]);
    }
    err = 0;

MIDF2Fin:
    alloc_acn(ctx, 0);
    alloc_actmp(ctx, 0);
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  get_ipart(n,xl,xu,ip)                                                   */
/*                                                                          */
/*  Calculate induced partition for [xl,xu] (i=1,...,n). Return in          */
/*  ip[j], j = 0,...,m-1. Return m, or -1 if error.                         */

int get_ipart(TDAContext *ctx, int n,double *xl,double *xu,double *ip)
{
    register int i;
    int m,nn;
    double xxl,xxu;

    m = 0;
    for (i = 1; i <= n; ++i) {
        xxl = xl[i];
        xxu = xu[i];
        if (xxl >= xxu - ctx->EPSI1) {
            mat_err(ctx, 11);
            return(-1);     
        }
        ip[m++] = xxl;
        ip[m++] = xxu;
    }
    if (sortd(ctx, m,ip,0))
        return(-1);    

    nn = 1;
    for (i = 1; i < m; ++i) {
        if (ip[i] > ip[i - 1])
            ip[nn++] = ip[i];
    }
    return(nn);
}

/*--------------------------------------------------------------------------*/
/*  m_midf3(cmd)    midf(XL,XU,XL1,XU1)                                     */
/*                                                                          */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_midf3(TDAContext *ctx, char *cmd)
{
    register int i,j;
    int err,n1,n2,row,col,idx,idx1,ivflg;
    register char *p;
    double xl,xu,xl1,xu1,tmp1,tmp2;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    if ((p = n_mexpr(ctx, cmd + 6,2,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MIDF3Fin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MIDF3Fin;
    if (col != 1) {
        mat_err(ctx, 7);
        goto MIDF3Fin;
    }
    if ((p = m_getmat(ctx, p,row,1,&idx,cmd,0)) == NULL)
        goto MIDF3Fin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MIDF3Fin;
    if ((p = m_getmat(ctx, p,row,1,&idx1,cmd,1)) == NULL)
        goto MIDF3Fin;

    for (i = 1; i <= row; ++i) {
        xl = ctx->MX[0][i];
        xu = ctx->MX[1][i];
        n1 = n2 = 0;
        tmp1 = tmp2 = 0.0;
        for (j = 1; j <= row; ++j) {
            xl1 = ctx->MX[0][j];
            xu1 = ctx->MX[1][j];
            if (xl1 >= xl && xl1 <= xu) {
                tmp1 += xl1;
                n1++;
            }
            if (xu1 >= xl && xu1 <= xu) {
                tmp2 += xu1;
                n2++;
            }
            if (n1 > 0)
                ctx->MatVal[idx][i] = tmp1 / (double)n1;
            if (n2 > 0)
                ctx->MatVal[idx1][i] = tmp2 / (double)n2;
        }
    }
    err = 0;

MIDF3Fin:
    mx_free(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gmin()          Global minimization.                                    */
/*                                                                          */
/*                  gmin(                                                   */
/*                      ns=...,         0 without derivatives, def. 0       */
/*                                      1 use derivatives                   */
/*                      xp=...,         starting values plus boxes          */
/*                      dsv=...,                                            */
/*                      mxit=...,       max number of iterations, def. 100  */
/*                      nbox=...,       max number of boxes, def. 100       */
/*                      tolbw=...,      tolerance box width, def. 1.e-4     */
/*                      tolfd=...,      tolerance for function range, 1e-10 */
/*                      tolfe=...,      for global min/max, def. 1.e-10     */
/*                      fmt=...,        print format, def. 13.6             */
/*                      prot=...,       protocol file                       */
/*                  ) = function;                                           */
/*                                                                          */
/*  It must be possible to evaluate the function with interval arithmetic.  */
/*  Function may not contain interval operators.                            */
/*                                                                          */  
/*  Return 0 if OK, otherwise -1.                                           */

int gmin(TDAContext *ctx)
{
    register int i,j,k;
    int err,r,n,l,gc;          
    double w,*uptr,*lptr;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Global minimization. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLBW = 1.e-4;       /* def tolerance for box length, branch and bound   */
    ctx->TOLFD = 1.e-10;      /* def tolerance for function range                 */
    ctx->TOLFE = 1.e-10;      /* def tolerance for global minimum                 */

    if (parm(ctx, ctx->CmdBuf + 4,12,1))   /* get parameters */
        goto GMINFin;

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 13,6);
           
    newline(ctx);

    if (ctx->FNArgN < 1) {
        printf1(ctx, "Error: function should contain at least one argument.\n");
        goto GMINFin;
    }
    if (alloc_par(ctx, ctx->FNArgN,1))
        goto GMINFin;

    if (get_dsv(ctx, ctx->FNArgN,ctx->RPar-1,1,ctx->RParL-1,ctx->RParU-1,1))   /* get starting values */ 
        goto GMINFin;

    /* print parameters and starting values */
     
    if (prn_fsval(ctx, ctx->FNArgN,ctx->RPar,1,ctx->RParL,ctx->RParU,0))
        goto GMINFin;

    if (ctx->MxIter < 1 || ctx->MxItFlg == 0)
        ctx->MxIter = 100;

    if (ctx->PMNS == 1)
        gc = 1;
    else
        gc = 0;

    if (ctx->PMNBOX < 1)         /* max number of boxes */
        ctx->PMNBOX = 100;

    /*  allocate box list for PMNBOX entries 

        GO_LB       lower bounds
        GO_UB       upper bounds
        GO_LBF      lower function values
        GO_UBF      upper function values
        GO_FLG      -1 if not yet processed
                     0 if dropped
                     1 if temporarily accepted
                     2 if finally accepted
    */

    if (alloc_list(ctx, ctx->PMNBOX,ctx->FNArgN)) 
        goto GMINFin;

    printf1(ctx, "Starting branch and bound algorithm.\n");
    printf1(ctx, "Maximum number of iterations: %d\n",ctx->MxIter);              
    printf1(ctx, "Maximum number of boxes: %d\n",ctx->PMNBOX);
    printf1(ctx, "Tolerance for box width: %g\n",ctx->TOLBW);
    printf1(ctx, "Tolerance for function range: %g\n",ctx->TOLFD);
    printf1(ctx, "Tolerance for global minimum: %g\n",ctx->TOLFE);
    if (gc) {
        printf1(ctx, "Using derivatives for monotonicity test.\n");

        if (fnd_alloc(ctx, 1,1,ctx->FNArgN,ctx->FNPN,1))
            goto GMINFin;
    }
    newline(ctx);

    r = igmin(ctx, 0,ctx->FNArgN,ctx->PMNBOX,ctx->RPar,ctx->RParL,ctx->RParU,ctx->MxIter,ctx->TOLBW,ctx->TOLFD,ctx->TOLFE,gc,0);
    if (r) {
        if (r == 1)
            printf1(ctx, "Exceeded maximum number of boxes.\n");
        goto GMINFin;
    }
    n = igmin_res(ctx, ctx->TOLFE,0);      
    newline(ctx);
              
    if (n > 0) {

        printf1(ctx, "Box  Acc  Width             Lower function bound  Upper function bound\n");
        r = 1;
        for (i = 0; i < ctx->GO_NB; ++i) {
            if (ctx->GO_FLG[i] == 2) {

                lptr = ctx->GO_LB + i * ctx->FNArgN;
                uptr = ctx->GO_UB + i * ctx->FNArgN;

                w = uptr[0] - lptr[0];
                for (j = 1; j < ctx->FNArgN; ++j)  
                    w = dmax(ctx, w,uptr[j] - lptr[j]);

                printf1(ctx, "%3d  %3d %17.10e %21.14e %21.14e\n",
                                    r++,ctx->GO_FLG[i],w,ctx->GO_LBF[i],ctx->GO_UBF[i]);
            }
        }

        l = 10;
        if (l < ctx->FNMLen)
            l = ctx->FNMLen;

        r = 1;
        for (i = 0; i < ctx->GO_NB; ++i) {
            if (ctx->GO_FLG[i] == 2) {
                printf1(ctx, "\nBox  Idx  Parameter\n");

                lptr = ctx->GO_LB + i * ctx->FNArgN;
                uptr = ctx->GO_UB + i * ctx->FNArgN;

                for (j = 0; j < ctx->FNArgN; ++j) {
                    k = ctx->FNArgSP[j];
                    printf1(ctx, "%3d  %3d  %s   ",r,j + 1,ctx->FNArgDef[k]);
                    prnchar(ctx, ' ',l - (int)strlen(ctx->FNArgDef[k]),0);
                    rt_printf1_d(ctx, ctx->PMFmtS,lptr[j]);
                    rt_printf1_d(ctx, ctx->PMFmtS,uptr[j]);
                    newline(ctx);
                }
                r++;
            }
        }
    }
    newline(ctx);
    err = 0;

GMINFin:
    alloc_list(ctx, 0,0);
    alloc_par(ctx, 0,0);
    fnd_alloc(ctx, 0,0,0,0,0);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  range()         Range of interval valued functions.                     */
/*                                                                          */
/*                  range(                                                  */
/*                      ns=...,         0 without derivatives, def. 0       */
/*                                      1 use derivatives                   */
/*                      xp=...,         starting values plus boxes          */
/*                      dsv=...,                                            */
/*                      mxit=...,       max number of iterations, def. 100  */
/*                      nbox=...,       max number of boxes, def. 100       */
/*                      tolbw=...,      tolerance box width, def. 1.e-6     */
/*                      tolfd=...,      tolerance for function range, 1e-10 */
/*                      tolfe=...,      for global min/max, def. 1.e-10     */
/*                      fmt=...,        print format, def. 13.6             */
/*                      prot=...,       protocol file                       */
/*                  ) = function;                                           */
/*                                                                          */
/*  It must be possible to evaluate the function with interval arithmetic.  */
/*  Function may contain interval operators (not completed yet).            */
/*                                                                          */  
/*  Return 0 if OK, otherwise -1.                                           */

int range(TDAContext *ctx)
{
    int err,r,ii,gc;              

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Range of interval-valued function. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLBW = 1.e-4;       /* def tolerance for box length, branch and bound   */
    ctx->TOLFD = 1.e-10;      /* def tolerance for function range                 */
    ctx->TOLFE = 1.e-10;      /* def tolerance for global minimum                 */

    if (parm(ctx, ctx->CmdBuf + 5,12,1))   /* get parameters */
        goto RNGFin;

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 13,6);
           
    newline(ctx);

    if (ctx->FNArgN < 1) {
        printf1(ctx, "Error: function should contain at least one argument.\n");
        goto RNGFin;
    }
    if (alloc_par(ctx, ctx->FNArgN,1))
        goto RNGFin;

    if (get_dsv(ctx, ctx->FNArgN,ctx->RPar-1,1,ctx->RParL-1,ctx->RParU-1,1))   /* get starting values */ 
        goto RNGFin;
    
    /* print parameters and starting values */
     
    if (prn_fsval(ctx, ctx->FNArgN,ctx->RPar,1,ctx->RParL,ctx->RParU,0))
        goto RNGFin;

    if (ctx->PMNS == 1) {
        gc = 1;
        if (fnd_alloc(ctx, 1,1,ctx->FNArgN,ctx->FNPN,1))
            goto RNGFin;
    }
    else
        gc = 0;

    if (ctx->MxIter < 1 || ctx->MxItFlg == 0)
        ctx->MxIter = 100;

    if (ctx->PMNBOX < 1)         /* max number of boxes */
        ctx->PMNBOX = 100;

    /*  allocate box list for PMNBOX entries 

        GO_LB       lower bounds
        GO_UB       upper bounds
        GO_LBF      lower function values
        GO_UBF      upper function values
        GO_FLG      -1 if not yet processed
                     0 if dropped
                     1 if temporarily accepted
                     2 if finally accepted
    */

    if (alloc_list(ctx, ctx->PMNBOX,ctx->FNArgN)) 
        goto RNGFin;

    printf1(ctx, "Starting branch and bound algorithm.\n");
    printf1(ctx, "Maximum number of iterations: %d\n",ctx->MxIter);              
    printf1(ctx, "Maximum number of boxes: %d\n",ctx->PMNBOX);
    printf1(ctx, "Tolerance for box width: %g\n",ctx->TOLBW);
    printf1(ctx, "Tolerance for function range: %g\n",ctx->TOLFD);
    printf1(ctx, "Tolerance for global minimum: %g\n",ctx->TOLFE);
    if (gc)
        printf1(ctx, "Using derivatives for monotonicity test.\n");

    for (ii = 0; ii <= 1; ++ii) {
        printf1(ctx, "\nFunction ");
        if (ii == 0)
            printf1(ctx, "minimization.\n");
        else
            printf1(ctx, "maximization.\n");

        r = igmin(ctx, ii,ctx->FNArgN,ctx->PMNBOX,ctx->RPar,ctx->RParL,ctx->RParU,ctx->MxIter,ctx->TOLBW,ctx->TOLFD,ctx->TOLFE,gc,0);

        if (r) {
            if (r == 1)
                printf1(ctx, "Exceeded maximum number of boxes.\n");
            goto RNGFin;
        }
        igmin_res(ctx, ctx->TOLFE,ii);
    }
    newline(ctx);
    err = 0;

RNGFin:
    alloc_list(ctx, 0,0);
    alloc_par(ctx, 0,0);
    fnd_alloc(ctx, 0,0,0,0,0);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  alloc_par(n,opt)    If opt != 0 allocate memory for n parameters and    */
/*                      bounds. Otherwise free previously allocated memory  */
/*                                                                          */
/*          Return 0 if OK, otherwise -1.                                   */

int alloc_par(TDAContext *ctx, int n,int opt)
{
    int err = 0;

    if (opt == 0)
        goto APARFin;

    err = -1;
    if (!(ctx->RPar = (double *)calloc((size_t)(n),sizeof(double)))) { 
        p_err(ctx, -2,1);
        goto APARFin;
    }   
    memrq(ctx, n,sizeof(double));
    ctx->RParA = n;             
              
    if (!(ctx->RParL = (double *)calloc((size_t)(n),sizeof(double)))) { 
        p_err(ctx, -2,1);
        goto APARFin;
    }   
    memrq(ctx, n,sizeof(double));
    ctx->RParLA = n;             
              
    if (!(ctx->RParU = (double *)calloc((size_t)(n),sizeof(double)))) { 
        p_err(ctx, -2,1);
        goto APARFin;
    }   
    memrq(ctx, n,sizeof(double));
    ctx->RParUA = n;             
    return(0);

APARFin:
    if (ctx->RParA > 0) {
        free((char *)ctx->RPar);
        memrq(ctx, -ctx->RParA,sizeof(double));
        ctx->RParA = 0;
    }
    if (ctx->RParLA > 0) {
        free((char *)ctx->RParL);
        memrq(ctx, -ctx->RParLA,sizeof(double));
        ctx->RParLA = 0;
    }
    if (ctx->RParUA > 0) {
        free((char *)ctx->RParU);
        memrq(ctx, -ctx->RParUA,sizeof(double));
        ctx->RParUA = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  igmin(opt,narg,nbmax,par,lb,ub,mxit,tolbw,tolfd,tolfe,gc,typ)           */
/*                                                                          */
/*              Global optimization with branch and bound method.           */
/*                                                                          */
/*              opt         0 if minimization                               */
/*                          1 if maximization                               */
/*              narg        number of arguments                             */
/*              nbmax       max number of boxes                             */
/*              par         initial parameters, par[i], i = 1,...,narg      */
/*              lb          initial lower bounds                            */
/*              ub          initial upper bounds                            */
/*              mxit        max number of iterations                        */
/*              tolbw       tolerance for box width                         */
/*              tolfd       tolerance for function range                    */
/*              tolfe       tolerance for global minimum                    */
/*              gc          0 don't use gradients, 1 use gradients          */
/*                                                                          */
/*  It is asssumed that a list for nbmax boxes is already allocated:        */
/*                                                                          */
/*      GO_LB       lower bounds                                            */
/*      GO_UB       upper bounds                                            */
/*      GO_LBF      lower function values                                   */
/*      GO_UBF      upper function values                                   */
/*      GO_FLG      -1 if not yet processed                                 */
/*                   0 if dropped                                           */
/*                   1 if temporarily accepted                              */
/*                   2 if finally accepted                                  */
/*                                                                          */
/*  If typ >= 1 use                                                         */
/*                                                                          */
/*      igmin_fun(opt,typ,narg,xl,xh,rl,rh,deriv,gl,gh)                     */
/*                                                                          */
/*  for function evaluation.                                                */
/*                                                                          */
/*  The igmin function sets the following global variables                  */
/*                                                                          */
/*  GO_FN       number of function evaluations                              */
/*  GO_IFN      number of inclusion function evalutaions                    */
/*  GO_IT       number of iterations performed                              */
/*  GO_NB       number of boxes in final list                               */
/*  GO_NBU      number of boxes used                                        */
/*  GO_FMIN     best function value                                         */
/*                                                                          */  
/*  The parameters for the best function value are returned in par[].       */
/*                                                                          */
/*  Return  0 if successful                                                 */
/*         -1 if insuff memory                                              */
/*         -2 error in function evaluation                                  */ 
/*          1 if exceeded max number of boxes                               */

int igmin(TDAContext *ctx, int opt,int narg,int nbmax,double *par,double *lb,double *ub,int mxit, double tolbw,double tolfd,double tolfe,int gc,int typ)
{
    (void)tolbw;        /* unused: the signature is shared */
    register int i,j;
    int err,r,nb,is,iter,ja,jb,jj,first;
    double fl,fu,tmp,w,f,fminp;
    double *lptra,*lptrb,*uptra,*uptrb;

    ctx->GO_IT = 0;              /* number of iterations performed */
    ctx->GO_FN = 0;              /* number of function evaluations */
    ctx->GO_IFN = 0;             /* number of inclusion function evaluations */
    ctx->GO_NBU = 0;             /* number of boxes used */
    ctx->GO_NB = 0;              /* number of boxes in final list */

    err = -1;

    for (i = 0; i < nbmax; ++i)
        ctx->GO_FLG[i] = 0;

    if (alloc_actmp(ctx, narg + 1))  
        goto IGMINFin; 
    if (gc) {
        if (alloc_acu(ctx, narg + 1))    /* lower bounds of gradient */
            goto IGMINFin; 
        if (alloc_acv(ctx, narg + 1))    /* upper bounds of gradient */
            goto IGMINFin; 
    }
    err = 0;

    lptra = ctx->GO_LB;                      /* put initial box on list */
    uptra = ctx->GO_UB;

    for (i = 0; i < narg; ++i) {
        lptra[i] = lb[i];
        uptra[i] = ub[i];
    }
    ctx->GO_FLG[0] = -1;     /* not yet processed */

    /* calculate inclusion function with initial boxes */
   
    ctx->GO_IFN++;

    if (typ == 0) {
        r = get_ifval(ctx, &fl,&fu,narg,lptra,uptra,0,ctx->AcU,ctx->AcV);
        if (r) {
            printf1(ctx, "Error in evaluating inclusion function.\n");
            prn_emsg2(ctx, r);              
            err = -2;
            goto IGMINFin;
        }
    }
    else {
        r = igmin_fun(ctx, 1,typ,narg,lptra,uptra,&fl,&fu,0,ctx->AcU,ctx->AcV);                                      
        if (r) {
            printf1(ctx, "Error in evaluating inclusion function (typ %d).\n",typ);
            err = -2;
            goto IGMINFin;
        }
    }
    ctx->GO_LBF[0] = fl;
    ctx->GO_UBF[0] = fu;
    
    /* calculate function at starting values */

    ctx->GO_FN++;
    if (typ == 0) {
        r = get_flval(ctx, &fminp,narg,par,0,0,&tmp,&tmp,&tmp);
        if (r) {      /* r from v_eval1(ctx) */
            printf1(ctx, "Error in function evaluation with starting values.\n");
            prn_emsg2(ctx, r);              
            err = -2;
            goto IGMINFin;
        }
    }
    else {
        r = igmin_fun(ctx, 0,typ,narg,par,&tmp,&fminp,&tmp,0,ctx->AcU,ctx->AcV);                                      
        if (r) {
            printf1(ctx, "Error in evaluating function (typ %d).\n",typ);
            err = -2;
            goto IGMINFin;
        }
    }
    nb = 1;     /* number of boxes */
    iter = 0;                                     

    if (ctx->SILENTFlg < 2)
        printfe(ctx, "\n  Iter    Function Value       NBox  FCall  IFCall\n");
    
    while (++iter <= mxit) {       
#ifdef TDA_R_PACKAGE
        /* an interrupt behaves exactly like exhausting the iteration
           limit: the loop ends, the ordinary bookkeeping runs, and the
           result is reported as uncertified with this line saying why */
        if (tda_check_interrupt()) {
            printf1(ctx, "Interrupted by user.\n");
            break;
        }
#endif

        if (ctx->SILENTFlg < 2)
            printfe(ctx, "%5d  %20.13e  %6d %6d %7d\n",iter,fminp,nb,ctx->GO_FN,ctx->GO_IFN);

        /* find box with lowest lower, or largest upper, function bound,
           and check for boxes that can be dropped */

        ja = -1;
        first = 1;
        for (i = 0; i < nb; ++i) {

            if (ctx->GO_FLG[i]) {    /* check all boxes */
  
                if (opt == 0) {
                    if (ctx->GO_LBF[i] > fminp)  
                        ctx->GO_FLG[i] = 0;
 
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
                else {

                    if (ctx->GO_UBF[i] < fminp)  
                        ctx->GO_FLG[i] = 0;
                    else if (ctx->GO_FLG[i] < 0) {
                        if (first) {      
                            f = ctx->GO_UBF[i];
                            ja = i;
                            first = 0;
                        }                 
                        else if (f < ctx->GO_UBF[i]) {
                            f = ctx->GO_UBF[i];
                            ja = i;
                        }   
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
        if (jb < 0) {       /* exceeded max number of boxes */
            err = 1;
            break;          
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

        if (w <= ctx->TOLBW) {                   /* temporarily accept this box */
            ctx->GO_FLG[ja] = 1;
            continue;
        }

        for (i = 0; i < narg; ++i) {      /* copy box ja to jb */
            lptrb[i] = lptra[i];
            uptrb[i] = uptra[i];
        }

        tmp = lptra[is] + w / 2.0;
        uptra[is] = lptrb[is] = tmp;
        ctx->GO_FLG[jb] = -1;

        for (jj = 0; jj < 2; ++jj) {           /* for both boxes */

            if (jj)   
                ja = jb;

            lptra = ctx->GO_LB + ja * narg;
            uptra = ctx->GO_UB + ja * narg;

            /* evaluate inclusion function */

            ctx->GO_IFN++;
            if (typ == 0) {
                r = get_ifval(ctx, &fl,&fu,narg,lptra,uptra,gc,ctx->AcU,ctx->AcV);
                if (r) {
                    printf1(ctx, "Error in evaluating inclusion function.\n");
                    prn_emsg2(ctx, r);              
                    err = -2;
                    goto IGMINFin;
                }
            }
            else {
                r = igmin_fun(ctx, 1,typ,narg,lptra,uptra,&fl,&fu,gc,ctx->AcU,ctx->AcV);                                      
                if (r) {
                    printf1(ctx, "Error in evaluating inclusion function (typ %d).\n",typ);
                    err = -2;
                    goto IGMINFin;
                }
            }
            if (opt == 0) {
                if (fl > fminp) {
                    ctx->GO_FLG[ja] = 0;             /* free box */
                    if (ja == nb - 1)
                        nb--;    
                    continue;
                }
            }
            else {
                if (fu < fminp) {
                    ctx->GO_FLG[ja] = 0;             /* free box */
                    if (ja == nb - 1)
                        nb--;    
                    continue;
                }
            }
            ctx->GO_LBF[ja] = fl;
            ctx->GO_UBF[ja] = fu;

            if (opt == 0) {
                if (fminp > fu) {
                    fminp = fu;
                    igmin_cpar(ctx, narg,uptra,par);
                }
            }
            else {
                if (fminp < fl) {
                    fminp = fl;
                    igmin_cpar(ctx, narg,lptra,par);
                }
            }
            if (fabs(fu - fl) < tolfd)      /* accept for solution */
                ctx->GO_FLG[ja] = 1;
   
            if (gc) {                       /* gradient check */

                r = 0;
                for (j = 0; j < narg; ++j) {
                    if (opt == 0) {
                        if (ctx->AcU[j] > 0.0) {
                            if (lptra[j] > lb[j]) {
                                r = 1;
                                break;
                            }
                            uptra[j] = lptra[j];
                            r = -1;
                        }
                        else if (ctx->AcV[j] < 0.0) {
                            if (uptra[j] < ub[j]) {
                                r = 1;
                                break;
                            }
                            lptra[j] = uptra[j];
                            r = -1;
                        }
                    }
                    else {
                        if (ctx->AcU[j] > 0.0) {
                            if (uptra[j] < ub[j]) {
                                r = 1;
                                break;
                            }
                            lptra[j] = uptra[j];
                            r = -1;
                        }
                        else if (ctx->AcV[j] < 0.0) {
                            if (lptra[j] > lb[j]) {
                                r = 1;
                                break;
                            }
                            uptra[j] = lptra[j];
                            r = -1;
                        }
                    }
                }
                if (r == 1) {
                    ctx->GO_FLG[ja] = 0;             /* free box */
                    if (ja == nb - 1)
                        nb--;    
                    continue;
                }
                else if (r == -1) {         /* update inclusion function */
                    ctx->GO_IFN++;
                    if (typ == 0) {
                        r = get_ifval(ctx, &fl,&fu,narg,lptra,uptra,0,&tmp,&tmp);
                        if (r) {
                            printf1(ctx, "Error in evaluating inclusion function.\n");
                            prn_emsg2(ctx, r);              
                            err = -2;
                            goto IGMINFin;
                        }
                    }
                    else {
                        r = igmin_fun(ctx, 1,typ,narg,lptra,uptra,&fl,&fu,0,ctx->AcU,ctx->AcV);                                      
                        if (r) {
                            printf1(ctx, "Error in evaluating inclusion function (typ %d).\n",typ);
                            err = -2;
                            goto IGMINFin;
                        }
                    }

                    if (opt == 0) {
                        if (fl > fminp) {
                            ctx->GO_FLG[ja] = 0;             /* free box */
                            if (ja == nb - 1)
                                nb--;    
                            continue;
                        }
                    }
                    else {
                        if (fu < fminp) {
                            ctx->GO_FLG[ja] = 0;             /* free box */
                            if (ja == nb - 1)
                                nb--;    
                            continue;
                        }
                    }
                    ctx->GO_LBF[ja] = fl;
                    ctx->GO_UBF[ja] = fu;
                }
            }

            /* update fminp with function value for midpoint of a box with
               least lower bound */
                   
            for (j = 0; j < narg; ++j)
                ctx->AcTmp[j] = (lptra[j] + uptra[j]) / 2.0;

            /* calculate function at midpoint */

            ctx->GO_FN++;
            if (typ == 0) {
                r = get_flval(ctx, &f,narg,ctx->AcTmp,0,0,&tmp,&tmp,&tmp);
                if (r) {      /* r from v_eval1(ctx) */
                    printf1(ctx, "Error in function evaluation.\n");
                    prn_emsg2(ctx, r);              
                    err = -2;
                    goto IGMINFin;
                }
            }
            else {
                r = igmin_fun(ctx, 0,typ,narg,ctx->AcTmp,&tmp,&f,&tmp,0,ctx->AcU,ctx->AcV);                                      
                if (r) {
                    printf1(ctx, "Error in evaluating function (typ %d).\n",typ);
                    err = -2;
                    goto IGMINFin;
                }
            }

            if (opt == 0) {
                if (fminp > f) {
                    fminp = f;
                    igmin_cpar(ctx, narg,ctx->AcTmp,par);
                }
            }   
            else {
                if (fminp < f) {
                    fminp = f;
                    igmin_cpar(ctx, narg,ctx->AcTmp,par);
                }
            }   
        }
    }
    ctx->GO_NB = nb;
    ctx->GO_IT = iter;
    ctx->GO_FMIN = fminp;

    if (ctx->PMProtFDef) {
        fprintf(ctx->PMProtFd,"IGMIN optimization.\n");
        prval(ctx, "Best function value",ctx->GO_FMIN);
        prvec(ctx, "Parameters",narg,par - 1);
        fprintf(ctx->PMProtFd,"\n");
        
        if (opt == 0)
            tmp = ctx->GO_FMIN - tolfe;
        else
            tmp = ctx->GO_FMIN + tolfe;

        ja = 1;
        for (i = 0; i < nb; ++i) {
            if (ctx->GO_FLG[i] == 0)  
                continue;

            if (ctx->GO_FLG[i] >= 1) {

                if (opt == 0) {                 /* minimization */
                    if (ctx->GO_LBF[i] >= tmp)         
                        ctx->GO_FLG[i] = 2;
                }
                else {
                    if (ctx->GO_UBF[i] <= tmp)         
                        ctx->GO_FLG[i] = 2;
                }
            }
            lptra = ctx->GO_LB + i * narg;
            uptra = ctx->GO_UB + i * narg;

            w = uptra[0] - lptra[0];
            for (j = 1; j < narg; ++j)  
                w = dmax(ctx, w,uptra[j] - lptra[j]);

            fprintf(ctx->PMProtFd,"Box   Width                 Lower function bound  Upper function bound  Acceptance\n");
            fprintf(ctx->PMProtFd,"%3d  %21.14e %21.14e %21.14e  %6d\n",
                                 ja++,w,ctx->GO_LBF[i],ctx->GO_UBF[i],ctx->GO_FLG[i]);

            fprintf(ctx->PMProtFd,"Parameter vector\nLB: ");
            for (j = 0; j < narg; ++j)  
                rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,lptra[j]);
            fprintf(ctx->PMProtFd,"\nUB: ");
            for (j = 0; j < narg; ++j)  
                rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,uptra[j]);
            fprintf(ctx->PMProtFd,"\n\n");
        }
    }

IGMINFin:
    alloc_acu(ctx, 0);           
    alloc_acv(ctx, 0);           
    alloc_actmp(ctx, 0);           
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  igmin_cpar(n,x,par) copy parameters from x[] into par[].                */

void igmin_cpar(TDAContext *ctx, int n,double *x,double *par)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;

    for (i = 0; i < n; ++i)
        par[i] = x[i];
}

/* ------------------------------------------------------------------------ */
/*  alloc_list(n,m)     Allocate box list with n entries, m parameters.     */
/*                      If n == 0 free previously allocated memory.         */
/*                      Return 0 if OK, -1 if error.                        */

int alloc_list(TDAContext *ctx, int n,int m)
{
    int err,nm;

    if (n <= 0) {
        err = 0;          
        goto AListFin;
    }
    err = -1;
    nm = n * m;
    if (!(ctx->GO_LB = (double *)calloc((size_t)(nm),sizeof(double)))) { 
        p_err(ctx, -2,1);
        goto AListFin;
    }
    memrq(ctx, nm,sizeof(double));
    ctx->GO_LBA = nm;

    if (!(ctx->GO_UB = (double *)calloc((size_t)(nm),sizeof(double)))) { 
        p_err(ctx, -2,1);
        goto AListFin;
    }
    memrq(ctx, nm,sizeof(double));
    ctx->GO_UBA = nm;

    if (!(ctx->GO_LBF = (double *)calloc((size_t)(n),sizeof(double)))) { 
        p_err(ctx, -2,1);
        goto AListFin;
    }
    memrq(ctx, n,sizeof(double));
    ctx->GO_LBFA = n;

    if (!(ctx->GO_UBF = (double *)calloc((size_t)(n),sizeof(double)))) { 
        p_err(ctx, -2,1);
        goto AListFin;
    }
    memrq(ctx, n,sizeof(double));
    ctx->GO_UBFA = n;

    if (!(ctx->GO_FLG = (short *)calloc((size_t)(n),sizeof(short)))) { 
        p_err(ctx, -2,1);
        goto AListFin;
    }
    memrq(ctx, n,sizeof(short));
    ctx->GO_FLGA = n;

    return(0);

AListFin:
    if (ctx->GO_LBA > 0) {
        free((char *)ctx->GO_LB);
        memrq(ctx, -ctx->GO_LBA,sizeof(double));
        ctx->GO_LBA = 0;
    }
    if (ctx->GO_UBA > 0) {
        free((char *)ctx->GO_UB);
        memrq(ctx, -ctx->GO_UBA,sizeof(double));
        ctx->GO_UBA = 0;
    }
    if (ctx->GO_LBFA > 0) {
        free((char *)ctx->GO_LBF);
        memrq(ctx, -ctx->GO_LBFA,sizeof(double));
        ctx->GO_LBFA = 0;
    }
    if (ctx->GO_UBFA > 0) {
        free((char *)ctx->GO_UBF);
        memrq(ctx, -ctx->GO_UBFA,sizeof(double));
        ctx->GO_UBFA = 0;
    }
    if (ctx->GO_FLGA > 0) {
        free((char *)ctx->GO_FLG);
        memrq(ctx, -ctx->GO_FLGA,sizeof(short));
        ctx->GO_FLGA = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  igmin_res(tolfe,opt)                                                    */
/*                                                                          */
/*  Get results from igmin(). If opt = 0 minimum, otherwise maximum         */
/*  Return number of finally accepted boxes.                                */

int igmin_res(TDAContext *ctx, double tolfe,int opt)
{
    register int i;
    int n0,n1,n2;
    double tmp,a;

    n0 = n1 = n2 = 0;         
    if (opt == 0)
        tmp = ctx->GO_FMIN - tolfe;
    else
        tmp = ctx->GO_FMIN + tolfe;
    a = tmp;

    for (i = 0; i < ctx->GO_NB; ++i) {
        if (ctx->GO_FLG[i] >= 1) {
            n1++;
            if (opt == 0) {
                if (ctx->GO_LBF[i] >= tmp) {       
                    n2++;
                    ctx->GO_FLG[i] = 2;
                }
                a = dmin(ctx, a,ctx->GO_LBF[i]);
            }
            else {
                if (ctx->GO_UBF[i] <= tmp) {       
                    n2++;
                    ctx->GO_FLG[i] = 2;
                }
                a = dmax(ctx, a,ctx->GO_UBF[i]);
            }
        }
        else if (ctx->GO_FLG[i] < 0)
            n0++;
    }

    printf1(ctx, "Number of iterations performed: %d\n",ctx->GO_IT);
    printf1(ctx, "Number of function evaluations: %d\n",ctx->GO_FN);
    printf1(ctx, "Number of inclusion function evaluations: %d\n",ctx->GO_IFN);
    printf1(ctx, "Number of boxes used: %d\n",ctx->GO_NBU);
    printf1(ctx, "Number of temporarily accepted boxes: %d\n",n1);
    printf1(ctx, "Number of finally accepted boxes: %d\n",n2);
    if (n0 > 0)
        printf1(ctx, "Warning: %d boxes have not been processed.\n",n0);

    if (opt == 0)  
        printf1(ctx, "\nBest minimal function value: ");
    else
        printf1(ctx, "\nBest maximal function value: ");
    rt_printf1_d(ctx, ctx->PMFmtS,ctx->GO_FMIN);
    if (opt == 0)  
        printf1(ctx, "  best lower bound: ");
    else
        printf1(ctx, "  best upper bound: ");
    if (n1 == 0 && n0 == 0) {
        /* the queue is exhausted and no candidate region survives:
           every box was excluded by an inclusion bound strictly above
           the then-current best, the best only improved afterwards,
           and it is attained -- so the best value IS the certified
           optimum.  Discard tests compare against the running best
           (see the GO_LBF[i] > fminp sites), which makes this sound.
           Whether any box also passes the final width/tolerance
           acceptance is a platform-marginal question (a few ulps
           around tolbw decide it), and the certificate must not
           depend on that. */
        rt_printf1_d(ctx, ctx->PMFmtS,ctx->GO_FMIN);
        newline(ctx);
        printf1(ctx, "(certified by exhaustion: every region was excluded\n");
        printf1(ctx, "against the best value.)\n");
    }
    else if (n1 == 0) {
        printf1(ctx, "***");
        newline(ctx);
        printf1(ctx, "No certified bound: the search stopped at its limits\n");
        printf1(ctx, "(mxit=, nbox=).  The value above is the best point found,\n");
        printf1(ctx, "not proven optimal.  Raise the limits to certify it.\n");
    }
    else
        rt_printf1_d(ctx, ctx->PMFmtS,a);
    newline(ctx);
    return(n2);
}

/* ------------------------------------------------------------------------ */
/*  i_add(xl,xh,yl,yh,rl,rh)                                                */
/*                                                                          */
/*  [rl,rh] = [xl,xh] + [yl,yh]                                             */

void i_add(TDAContext *ctx, double xl,double xh,double yl,double yh,double *rl,double *rh)
{
    (void)ctx;        /* unused: the signature is shared */
    *rl = xl + yl;
    *rh = xh + yh;
}

/* ------------------------------------------------------------------------ */
/*  i_sub(xl,xh,yl,yh,rl,rh)                                                */
/*                                                                          */
/*  [rl,rh] = [xl,xh] - [yl,yh]                                             */

void i_sub(TDAContext *ctx, double xl,double xh,double yl,double yh,double *rl,double *rh)
{
    (void)ctx;        /* unused: the signature is shared */
    *rl = xl - yh;
    *rh = xh - yl;
}

/* ------------------------------------------------------------------------ */
/*  i_mul(xl,xh,yl,yh,rl,rh)                                                */
/*                                                                          */
/*  [rl,rh] = [xl,xh] * [yl,yh]                                             */

void i_mul(TDAContext *ctx, double xl,double xh,double yl,double yh,double *rl,double *rh)
{
    if (0.0 <= xl) {
        if (0.0 <= yl) {
            *rl = xl * yl;
            *rh = xh * yh;
        }
        else if (yh <= 0.0) {
            *rl = xh * yl;
            *rh = xl * yh;
        }
        else {
            *rl = xh * yl;
            *rh = xh * yh;
        }
    }
    else if (xh <= 0.0) {
        if (0.0 <= yl) {
            *rl = xl * yh;
            *rh = xh * yl;
        }
        else if (yh <= 0.0) {
            *rl = xh * yh;
            *rh = xl * yl;
        }
        else {
            *rl = xl * yh;
            *rh = xl * yl;
        }
    }
    else {
        if (0.0 <= yl) {
            *rl = xl * yh;
            *rh = xh * yh;
        }
        else if (yh <= 0.0) {
            *rl = xh * yl;
            *rh = xl * yl;
        }
        else {
            *rl = dmin(ctx, xl * yh,xh * yl);
            *rh = dmax(ctx, xl * yl,xh * yh);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  i_div(xl,xh,yl,yh,rl,rh)                                                */
/*                                                                          */
/*  [rl,rh] = [xl,xh] / [yl,yh]                                             */
/*                                                                          */
/*  Return 0 if OK, -1 if error                                             */

int i_div(TDAContext *ctx, double xl,double xh,double yl,double yh,double *rl,double *rh)
{
    (void)ctx;        /* unused: the signature is shared */
    if (0.0 <= xl) {
        if (0.0 < yl) {
            *rl = xl / yh;
            *rh = xh / yl;
        }
        else if (yh < 0.0) {
            *rl = xh / yh;
            *rh = xl / yl;
        }
        else  
            return(-1);    
    }
    else if (xh <= 0.0) {
        if (0.0 < yl) {
            *rl = xl / yl;
            *rh = xh / yh;
        }
        else if (yh < 0.0) {
            *rl = xh / yl;
            *rh = xl / yh;
        }
        else  
            return(-1);     
    }
    else {
        if (0.0 < yl) {
            *rl = xl / yl;
            *rh = xh / yl;
        }
        else if (yh < 0.0) {
            *rl = xh / yh;
            *rh = xl / yh;
        }
        else  
            return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  i_abs(xl,xh,rl,rh)                                                      */
/*                                                                          */
/*  [rl,rh] = abs([xl,xh])                                                  */

void i_abs(TDAContext *ctx, double xl,double xh,double *rl,double *rh)
{
    if (xl >= 0.0) {
        *rl = xl;
        *rh = xh;
    }
    else if (xh <= 0.0) {
        *rl = -xh;
        *rh = -xl;
    }
    else {
        *rl = 0.0;
        *rh = dmax(ctx, -xl,xh);
    }
}

/* ------------------------------------------------------------------------ */
/*  i_max(xl,xh,yl,yh,rl,rh)                                                */
/*                                                                          */

void i_max(TDAContext *ctx, double xl,double xh,double yl,double yh,double *rl,double *rh)
{

    *rl = dmax(ctx, xl,yl);
    *rh = dmax(ctx, xh,yh);
}

/* ------------------------------------------------------------------------ */
/*  i_min(xl,xh,yl,yh,rl,rh)                                                */
/*                                                                          */

void i_min(TDAContext *ctx, double xl,double xh,double yl,double yh,double *rl,double *rh)
{

    *rl = dmin(ctx, xl,yl);
    *rh = dmin(ctx, xh,yh);
}

/* ------------------------------------------------------------------------ */
/*  i_square(xl,xh,rl,rh)                                                   */
/*                                                                          */

void i_square(TDAContext *ctx, double xl,double xh,double *rl,double *rh)
{
    if (xl >= 0.0) {
        *rl = xl * xl;
        *rh = xh * xh;
    }
    else if (xh <= 0.0) {
        *rl = xh * xh;
        *rh = xl * xl;
    }
    else {
        *rl = 0.0;          
        *rh = dmax(ctx, xl * xl,xh * xh);
    }
}

/* ------------------------------------------------------------------------ */
/*  i_sqrt(xl,xh,rl,rh)                                                     */
/*                                                                          */
/*  [rl,rh] = [sqrt(xl),sqrt(xh)]                                           */

void i_sqrt(TDAContext *ctx, double xl,double xh,double *rl,double *rh)
{
    (void)ctx;        /* unused: the signature is shared */
    if (xl < 0.0) {
        tda_err("ERROR in i_sqrt.\n");
        exit(0);
    }
    *rl = sqrt(xl);
    *rh = sqrt(xh);
}

/* ------------------------------------------------------------------------ */
/*  igmin_fun(opt,typ,n,xl,xh,rl,rh,deriv,gl,gh)                            */
/*                                                                          */
/*  Calculate function depending on typ and opt.                            */
/*  If opt = 0 : standard function with arguments xl[i], i=0,...,n-1        */
/*               return function value in rl.                               */
/*  If opt = 1 : inclusion function at [xl[i],xh[i]], i = 0,...,n-1         */
/*               return interval in [rl,rh].                                */
/*                                                                          */
/*  If deriv != 0 and opt=1, calculate bounds for gradient in gl[i],gh[i].  */
/*                                                                          */
/*  typ = 1  -- variance                                                    */
/*        2  -- Gini coefficient (uses IGMINFUNVal as argument)             */
/*        3  -- multi-dimensional scaling (mds1_fn  in t_mds)               */
/*                                                                          */  
/*  Return 0 if OK, -1 if error                                             */

int igmin_fun(TDAContext *ctx, int opt,int typ,int n,double *xl,double *xh, double *rl,double *rh,int deriv,double *gl,double *gh)
{
    register int i,j,k;
    double x,y,x1,y1,tl,th,sl,sh,tmp;
/**   
    for (i = 0; i < n; ++i)
        printf1(ctx, "i=%d xl=%lf xh=%lf\n",i,xl[i],xh[i]);
**/      

    switch (typ) {
        case 1:                 /* variance */

            if (opt == 0) {
                x = 0.0;
                for (i = 0; i < n; ++i)
                    x += xl[i];
                x /= (double)n;
                *rl = 0.0;
                for (i = 0; i < n; ++i) {
                    y = (xl[i] - x);
                    *rl += y * y;
                }
                *rl /= (double)n;
                return(0);
            }
            *rl = *rh = x = y = 0.0;
            for (i = 0; i < n; ++i) {
                x += xl[i];
                y += xh[i];
            }
            x /= (double)n;
            y /= (double)n;
            for (i = 0; i < n; ++i) {
                i_sub(ctx, xl[i],xh[i],x,y,&tl,&th);
                i_mul(ctx, tl,th,tl,th,&sl,&sh);
                i_add(ctx, *rl,*rh,sl,sh,&tl,&th);
                *rl = tl;
                *rh = th;
            }
            *rl /= (double)n;
            *rh /= (double)n;

            if (deriv) {
                for (i = 0; i < n; ++i) {
                    gl[i] = 2.0 * (xl[i] - y) / (double)n;
                    gh[i] = 2.0 * (xh[i] - x) / (double)n;
                }
            }
            return(0);
          
        case 2:                 /* Gini */

            if (opt == 0) {
                y = x = 0.0;
                for (i = 0; i < n; ++i) {
                    x += xl[i];
                    if (xl[i] <= ctx->IGMINFUNVal)
                        y += xl[i];
                }
                *rl = y / x;
                return(0);
            }
            *rl = *rh = x = y = x1 = y1 = 0.0;
            for (i = 0; i < n; ++i) {
                x += xl[i];
                y += xh[i];
                if (xl[i] <= ctx->IGMINFUNVal)
                    x1 += xl[i];
                if (xh[i] >= ctx->IGMINFUNVal)
                    y1 += xh[i];
            }
            i_div(ctx, x1,y1,x,y,rl,rh);

            if (deriv) {
                for (i = 0; i < n; ++i) {
                    gl[i] = gh[i] = -1.0;
                }
            }
            return(0);

        case 3:                 /* multi-dimensional scaling */

            if (opt == 0) {
                *rl = 0.0;
                for (i = 1; i < ctx->GD_NP; ++i) {
                    for (j = 0; j < i; ++j) {
                        x = gdd_adj(ctx, i,j,ctx->PMGN);
                        if (x >= 0.0) {
                            tmp = x - fabs(xl[i] - xl[j]);
                            *rl += tmp * tmp;
                        }
                    }
                }
                return(0);
            }
            *rl = *rh = 0.0;

            for (i = 1; i < ctx->GD_NP; ++i) {
                for (j = 0; j < i; ++j) {
                    x = gdd_adj(ctx, i,j,ctx->PMGN);
                    if (x >= 0.0) {
                        i_sub(ctx, xl[i],xh[i],xl[j],xh[j],&tl,&th);
                        i_abs(ctx, tl,th,&sl,&sh);
                        i_sub(ctx, x,x,sl,sh,&tl,&th);
                        i_mul(ctx, tl,th,tl,th,&sl,&sh);
                        i_add(ctx, *rl,*rh,sl,sh,&tl,&th);
                        *rl = tl;
                        *rh = th;
                    }
                }
            }
            if (deriv) {
                for (k = 0; k < n; ++k) {
                    gl[k] = gh[k] = 0;
                    for (i = 1; i < ctx->GD_NP; ++i) {
                        for (j = 0; j < i; ++j) {
                            if (i == k || j == k) {
                                x = gdd_adj(ctx, i,j,ctx->PMGN);
                                if (x >= 0.0) {
                                    i_sub(ctx, xl[i],xh[i],xl[j],xh[j],&tl,&th);
                                    i_abs(ctx, tl,th,&sl,&sh);
                                    i_sub(ctx, x,x,sl,sh,&tl,&th);

                                    if (k == i) {
                                        if (xh[i] + ctx->EPSI1 < xl[j])        
                                            i_mul(ctx, tl,th,2,2,&sl,&sh);
                                        else if (xh[j] + ctx->EPSI1 < xl[i])   
                                            i_mul(ctx, tl,th,-2,-2,&sl,&sh);
                                        else
                                            i_mul(ctx, tl,th,-2,2,&sl,&sh);
                                    }
                                    else {
                                        if (xh[j] + ctx->EPSI1 < xl[i])        
                                            i_mul(ctx, tl,th,2,2,&sl,&sh);
                                        else if (xh[i] + ctx->EPSI1 < xl[j])   
                                            i_mul(ctx, tl,th,-2,-2,&sl,&sh);
                                        else
                                            i_mul(ctx, tl,th,-2,2,&sl,&sh);
                                    }
                                    tl = gl[k];
                                    th = gh[k];
                                    i_add(ctx, tl,th,sl,sh,&x,&y);
                                    gl[k] = x;
                                    gh[k] = y;
                                }
                            }
                        }
                    }
                }
            }
            return(0);
          
        case 4:                 /* covariance: first n/2 coordinates are
                                   the x observations, last n/2 the y */
        case 5: {               /* correlation = cov / sqrt(varx*vary) */
            int m = n / 2;
            double mxl,mxh,myl,myh,cl,ch,vxl,vxh,vyl,vyh,dl,dh;
            if (m < 2 || 2 * m != n)
                return(-1);

            if (opt == 0) {     /* point evaluation at xl[] */
                x = y = 0.0;
                for (i = 0; i < m; ++i) {
                    x += xl[i];
                    y += xl[m + i];
                }
                x /= (double)m;
                y /= (double)m;
                tl = sl = sh = 0.0;
                for (i = 0; i < m; ++i) {
                    tl += (xl[i] - x) * (xl[m + i] - y);
                    sl += (xl[i] - x) * (xl[i] - x);
                    sh += (xl[m + i] - y) * (xl[m + i] - y);
                }
                tl /= (double)m;
                if (typ == 4) {
                    *rl = tl;
                    return(0);
                }
                sl /= (double)m;
                sh /= (double)m;
                if (sl <= 0.0 || sh <= 0.0)
                    return(-1);
                *rl = tl / sqrt(sl * sh);
                return(0);
            }

            /* interval inclusion: means of the halves, then the sums */
            mxl = mxh = myl = myh = 0.0;
            for (i = 0; i < m; ++i) {
                mxl += xl[i];      mxh += xh[i];
                myl += xl[m + i];  myh += xh[m + i];
            }
            mxl /= (double)m; mxh /= (double)m;
            myl /= (double)m; myh /= (double)m;

            cl = ch = 0.0;
            vxl = vxh = vyl = vyh = 0.0;
            for (i = 0; i < m; ++i) {
                i_sub(ctx, xl[i],xh[i],mxl,mxh,&tl,&th);
                i_sub(ctx, xl[m + i],xh[m + i],myl,myh,&sl,&sh);
                i_mul(ctx, tl,th,sl,sh,&dl,&dh);
                cl += dl; ch += dh;
                i_mul(ctx, tl,th,tl,th,&dl,&dh);
                vxl += dl; vxh += dh;
                i_mul(ctx, sl,sh,sl,sh,&dl,&dh);
                vyl += dl; vyh += dh;
            }
            cl /= (double)m; ch /= (double)m;

            if (typ == 4) {
                *rl = cl;
                *rh = ch;
                if (deriv) {
                    /* d cov/d x_i = (y_i - my)/m, d cov/d y_i = (x_i - mx)/m */
                    for (i = 0; i < m; ++i) {
                        gl[i] = (xl[m + i] - myh) / (double)m;
                        gh[i] = (xh[m + i] - myl) / (double)m;
                        gl[m + i] = (xl[i] - mxh) / (double)m;
                        gh[m + i] = (xh[i] - mxl) / (double)m;
                    }
                }
                return(0);
            }

            /* correlation: variances can touch 0 inside the box, in
               which case the quotient is unbounded -- but correlation
               itself never leaves [-1,1], so intersecting with it is
               both valid and the only honest answer there */
            vxl /= (double)m; vxh /= (double)m;
            vyl /= (double)m; vyh /= (double)m;
            if (vxl < 0.0) vxl = 0.0;
            if (vyl < 0.0) vyl = 0.0;
            i_mul(ctx, vxl,vxh,vyl,vyh,&dl,&dh);
            if (dl < 0.0) dl = 0.0;
            i_sqrt(ctx, dl,dh,&sl,&sh);
            if (sl <= 0.0) {
                *rl = -1.0;
                *rh = 1.0;
            }
            else {
                i_div(ctx, cl,ch,sl,sh,rl,rh);
                if (*rl < -1.0) *rl = -1.0;
                if (*rh > 1.0) *rh = 1.0;
            }
            if (deriv) {
                /* quotient rule, interval-evaluated:
                   d corr/d x_i = (y_i - my)/(m s) - c (x_i - mx) vy/(m s^3)
                   d corr/d y_i symmetric; s = sqrt(vx vy).  Crude but
                   correct bounds are all the monotonicity test needs. */
                double s3l,s3h,t1l,t1h,t2l,t2h,u1l,u1h;
                i_mul(ctx, sl,sh,sl,sh,&t1l,&t1h);
                i_mul(ctx, t1l,t1h,sl,sh,&s3l,&s3h);
                for (i = 0; i < m; ++i) {
                    i_sub(ctx, xl[m + i],xh[m + i],myl,myh,&t1l,&t1h);
                    i_div(ctx, t1l,t1h,(double)m * sl,(double)m * sh,&t1l,&t1h);
                    i_sub(ctx, xl[i],xh[i],mxl,mxh,&t2l,&t2h);
                    i_mul(ctx, cl,ch,t2l,t2h,&t2l,&t2h);
                    i_mul(ctx, t2l,t2h,vyl,vyh,&t2l,&t2h);
                    i_div(ctx, t2l,t2h,(double)m * s3l,(double)m * s3h,&t2l,&t2h);
                    i_sub(ctx, t1l,t1h,t2l,t2h,&u1l,&u1h);
                    gl[i] = u1l; gh[i] = u1h;

                    i_sub(ctx, xl[i],xh[i],mxl,mxh,&t1l,&t1h);
                    i_div(ctx, t1l,t1h,(double)m * sl,(double)m * sh,&t1l,&t1h);
                    i_sub(ctx, xl[m + i],xh[m + i],myl,myh,&t2l,&t2h);
                    i_mul(ctx, cl,ch,t2l,t2h,&t2l,&t2h);
                    i_mul(ctx, t2l,t2h,vxl,vxh,&t2l,&t2h);
                    i_div(ctx, t2l,t2h,(double)m * s3l,(double)m * s3h,&t2l,&t2h);
                    i_sub(ctx, t1l,t1h,t2l,t2h,&u1l,&u1h);
                    gl[m + i] = u1l; gh[m + i] = u1h;
                }
            }
            return(0);
        }

        default:                 
            break;
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  idf()           Distribution and density function for interval-valued   */
/*                  variables.                                              */
/*                                                                          */  
/*                  idf(                                                    */
/*                      opt=...,    1 = min,max,mean                        */
/*                                  2 = self-cons. distribution             */
/*                      mxit=...,   max number of iterations, def. 50       */
/*                                                                          */  
/*                  ) = XL,XU;                                              */
/*                                                                          */
/*  There must be exactly two variables on the right-hand side to be        */  
/*  interpreted as lower and upper bounds of an interval-valued variable.   */  
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int idf(TDAContext *ctx)
{
    register int i,j,k;
    int err,il,ih,n,nn,ai,bi,iter,suc = 0;
    double xx,xl,xh,tmp,f,fl,fh,c;
         
    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Interval-valued distribution. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLF = 0.001;                /* as iddf(); tolf= overrides */

    if (parm(ctx, ctx->CmdBuf + 3,4,1))     /* get parameters */
        goto IDFFin;

    if (ctx->PMOPT != 2)
        ctx->PMOPT = 1;

    if (ctx->PMNV != 2) {            /* need two variables */
        p_err(ctx, -1,1);
        goto IDFFin;
    }
    il = ctx->PMVIdx[0];
    ih = ctx->PMVIdx[1];

    if (alloc_acx(ctx, ctx->NOC + 1))     /* lower bound */ 
        goto IDFFin;
    if (alloc_acy(ctx, ctx->NOC + 1))     /* upper bound */
        goto IDFFin;
    if (alloc_actmp(ctx, 2 * ctx->NOC + 1))     /* induced partition */
        goto IDFFin;

    n = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        xl = get_data(ctx, il,i);
        xh = get_data(ctx, ih,i);
        if (xh <= xl + ctx->EPSI1) {
            printf1(ctx, "Error in case %d, found: %g,%g\n",i+1,xl,xh);
            goto IDFFin;
        }
        ctx->AcX[i] = xl;
        ctx->AcY[i] = xh;

        /* printf1("i=%d x=%g y=%g \n",i,AcX[i],AcY[i]); */
         
        ctx->AcTmp[n++] = xl;
        ctx->AcTmp[n++] = xh;
    }
    if (sortd(ctx, n,ctx->AcTmp,0))
        goto IDFFin;
   
    nn = 1;
    for (i = 1; i < n; ++i) {           /* AcTmp contains induced partition */
        if (ctx->AcTmp[i] > ctx->AcTmp[i - 1])
            ctx->AcTmp[nn++] = ctx->AcTmp[i];
    }
    if (sortd2(ctx, ctx->NOC,ctx->AcX,ctx->AcY))            /* sort observations */
        goto IDFFin;

    /* calculate distribution functions */

    if (alloc_acz(ctx, nn))     /* mean df */
        goto IDFFin;

    printf1(ctx, "\n  Idx      Partition  Lower Bound  Upper Bound      Mean DF\n");

    for (i = 0; i < nn; ++i) {
        tmp = ctx->AcTmp[i];
        fl = fh = f = 0.0;
        for (j = 0; j < ctx->NOC; ++j) {
            xl = ctx->AcX[j];
            xh = ctx->AcY[j];
            if (xh <= tmp)
                fl += 1.0;
            if (xl <= tmp) {
                fh += 1.0;
                if (xl < tmp) {
                    xx = dmin(ctx, xh,tmp);
                    f += (xx - xl) / (xh - xl);
                }
            }
            else
                break;
        }
        ctx->AcZ[i] = f / (double)ctx->NOC;

        printf1(ctx, "%5d ",i + 1);
        printf1(ctx, "%14.6lf ",tmp);
        printf1(ctx, "%12.6lf ",fl / (double)ctx->NOC);
        printf1(ctx, "%12.6lf ",fh / (double)ctx->NOC);
        printf1(ctx, "%12.6lf\n",f  / (double)ctx->NOC);
#ifdef TDA_R_PACKAGE
        /* idf prints its table inline (no prn_* function); same row,
           full precision, for the idf.table export */
        {
            double erow[5];
            erow[0] = (double)(i + 1);
            erow[1] = tmp;
            erow[2] = fl / (double)ctx->NOC;
            erow[3] = fh / (double)ctx->NOC;
            erow[4] = f / (double)ctx->NOC;
            tda_export_row(ctx, "idf.table", erow, 5);
        }
#endif

    }
    newline(ctx);
#ifdef TDA_R_PACKAGE
    /* flush per printed table: the self-consistent iterations below
       reprint it, and each print becomes its own export (.2, ...) */
    tda_export_flush(ctx, "idf.table");
#endif

    if (ctx->PMOPT != 2) {
        err = 0;
        goto IDFFin;
    }

    /* iterations for self-consistent df */

    if (ctx->MxIter < 1 || ctx->MxItFlg == 0)
        ctx->MxIter = 50;

    printf1(ctx, "Iterative calculation of self-consistent distribution.\n");
    printf1(ctx, "Max number of iterations: %d\n",ctx->MxIter);
    printf1(ctx, "Tolerance for convergence: %g\n",ctx->TOLF);

    if (alloc_acu(ctx, nn)) 
        goto IDFFin;

    for (iter = 0; iter < ctx->MxIter; ++iter) {


        for (i = 0; i < nn; ++i) {

            xx = ctx->AcTmp[i];

/*   printf1("i=%d xx=%g\n",i,xx);
*/

            tmp = 0.0;

            for (j = 0; j < ctx->NOC; ++j) {
                xl = ctx->AcX[j];
                xh = ctx->AcY[j];
                ai = bi = -1;

                if (xl > xx)  
                    ;
                else if (xh < xx)
                    tmp += 1.0;
                else {
                    ai = bi = -1;
                    for (k = 0; k < nn; ++k) {
                        if (fabs(xl - ctx->AcTmp[k]) < 0.01)  
                            ai = k;
                        if (fabs(xh - ctx->AcTmp[k]) < 0.01)  
                            bi = k;
                    }
                    if (ai < 0 || bi < 0) {
                        printf2(ctx, "ERROR\n");
                        gerr_exit(ctx, 223);
                    }   
                    tmp += (ctx->AcZ[i] - ctx->AcZ[ai]) / (ctx->AcZ[bi] - ctx->AcZ[ai]);
                }
                /****
                printf1(ctx, "j=%d xl=%lf xh=%lf tmp=%lf ai=%2d bi=%2d\n",
                            j,xl,xh,tmp,ai,bi);
                ***/

            }
            tmp /= (double)ctx->NOC;
            ctx->AcU[i] = tmp;
        }
        /* The iterate becomes the current distribution, and the
           iteration stops when the largest change falls under tolf --
           the same rule iddf() applies to its discrete counterpart.  The
           original had the update and any report of the result inside a
           comment, so opt=2 announced the iteration and showed nothing;
           see doc/changes-from-tda.md. */
        c = 0.0;
        for (i = 0; i < nn; ++i) {
            c = dmax(ctx, c,fabs(ctx->AcU[i] - ctx->AcZ[i]));
            ctx->AcZ[i] = ctx->AcU[i];
        }
        if (c <= ctx->TOLF) {
            suc = 1;
            break;
        }
    }
    printf1(ctx, "Convergence ");
    if (suc == 0)
        printf1(ctx, "not ");
    printf1(ctx, "reached after %d iterations.\n",iter);

    printf1(ctx, "\n  Idx      Partition           DF\n");
    for (i = 0; i < nn; ++i) {
        printf1(ctx, "%5d ",i + 1);
        printf1(ctx, "%14.6lf ",ctx->AcTmp[i]);
        printf1(ctx, "%12.6lf\n",ctx->AcZ[i]);
#ifdef TDA_R_PACKAGE
        {
            double erow[3];
            erow[0] = (double)(i + 1);
            erow[1] = ctx->AcTmp[i];
            erow[2] = ctx->AcZ[i];
            tda_export_row(ctx, "idf.scdf", erow, 3);
        }
#endif
    }
    newline(ctx);
    err = 0;

IDFFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sddf()          Distribution of set-valued discrete variable.           */
/*                                                                          */  
/*                  sddf(                                                   */
/*                      opt=...,        1 only min/max/mean df              */
/*                                      2 plus self-consistent df           */
/*                      mxit=...,       max iterations, def. 50             */
/*                      tolf=...,       tolerance, def. 1e-3                */
/*                      df=...,         print to output file                */
/*                      prot=...,       protocol file                       */
/*                  ) = X1,...,Xm;                                          */
/*                                                                          */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int sddf(TDAContext *ctx)
{
    register int i,j,k;
    int err,m,n,iter,suc;
    double c,p,tmp;

    err = -1;         
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Distribution of set-valued discrete variable. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLF = 0.001;

    if (parm(ctx, ctx->CmdBuf + 4,4,1))     /* get parameters */
        goto SDDFFin;

    if (ctx->PMOPT != 2)
        ctx->PMOPT = 1;

    if (ctx->MxIter < 1 || ctx->MxItFlg == 0)
        ctx->MxIter = 50;

    m = ctx->PMNV;                   /* number of variables = categories */

    if (alloc_acn(ctx, m + 1))       /* lower bound */ 
        goto SDDFFin;
    if (alloc_acm(ctx, m + 1))       /* upper bound */
        goto SDDFFin;
    if (alloc_acx(ctx, m + 1))       /* mean distribution */
        goto SDDFFin;
    if (alloc_ack(ctx, m + 1))                           
        goto SDDFFin;


    for (i = 0; i < ctx->NOC; ++i) {
        k = -1;
        n = 0;
        for (j = 1; j <= m; ++j) {
            if ((get_data(ctx, ctx->PMVIdx[j-1],i)) != 0.0) {
                k = j;
                n++;
                ctx->AcM[k] += 1;
                ctx->AcK[k] = 1;
            }
            else
                ctx->AcK[j] = 0;
        }
        if (n == 1)  
            ctx->AcN[k] += 1;
        for (j = 1; j <= m; ++j)
            ctx->AcX[j] += (double)ctx->AcK[j] / (double)n;
    }
    if (ctx->PMOPT != 2 || ctx->PMF1Def == 0) {
        prn_ddf(ctx, 0,m,ctx->AcN,ctx->AcM,ctx->AcX,ctx->NOC);
        if (ctx->PMOPT != 2) {
            err = 0;
            goto SDDFFin;
        }
    }
    if (ctx->PMF1Def) {
        if (alloc_acu(ctx, m + 1))   
            goto SDDFFin;

        for (j = 1; j <= m; ++j)
            ctx->AcU[j] = ctx->AcX[j];
    }

    /* iterations for self-consistent df */

    printf1(ctx, "Iterative calculation of self-consistent distribution.\n");
    printf1(ctx, "Max number of iterations: %d\n",ctx->MxIter);
    printf1(ctx, "Tolerance for convergence: %g\n",ctx->TOLF);

    if (alloc_acy(ctx, m + 1))   
        goto SDDFFin;
    if (alloc_acz(ctx, m + 1))   
        goto SDDFFin;
    
    suc = 0;
    for (iter = 0; iter < ctx->MxIter; ++iter) {
    
        for (i = 0; i < ctx->NOC; ++i) {
            p = 0.0;
            for (j = 1; j <= m; ++j) {
                if ((get_data(ctx, ctx->PMVIdx[j-1],i)) != 0.0) {
                    p += ctx->AcX[j];
                    ctx->AcZ[j] = ctx->AcX[j];
                }
                else
                    ctx->AcZ[j] = 0.0;
            }
            for (j = 1; j <= m; ++j) {
                if (p > 0.0)
                    ctx->AcY[j] += ctx->AcZ[j] / p;
            }
        }
        if (ctx->PMProtFDef)  
            fprintf(ctx->PMProtFd,"Iter%3d ",iter);

        c = 0.0;
        for (j = 1; j <= m; ++j) {
            tmp = ctx->AcX[j];
            ctx->AcX[j] = ctx->AcY[j] / (double)ctx->NOC;
            c = dmax(ctx, c,fabs(ctx->AcX[j] - tmp));
            ctx->AcY[j] = 0.0;
            if (ctx->PMProtFDef)  
                rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,ctx->AcX[j]);
        }
        if (ctx->PMProtFDef)  
            fprintf(ctx->PMProtFd,"\n");
    
        if (c <= ctx->TOLF) {
            suc = 1;
            break;
        }
    }
    printf1(ctx, "Convergence ");
    if (suc == 0)
        printf1(ctx, "not ");
    printf1(ctx, "reached after %d iterations.\n",iter);

    if (ctx->PMF1Def)
        prn_ddf2(ctx, 0,m,ctx->AcN,ctx->AcM,ctx->AcU,ctx->AcX,ctx->NOC);
    else
        prn_ddf1(ctx, 0,m,ctx->AcX);

    err = 0;

SDDFFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  iddf()          Distribution of interval-valued discrete variable.      */
/*                                                                          */  
/*                  iddf(                                                   */
/*                      opt=...,        1 only min/max/mean df              */
/*                                      2 plus self-consistent df           */
/*                      mxit=...,       max iterations, def. 50             */
/*                      tolf=...,       tolerance, def. 1e-3                */
/*                      df=...,         print to output file                */
/*                      prot=...,       protocol file                       */
/*                  ) = XL,XH;                                              */
/*                                                                          */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int iddf(TDAContext *ctx)
{
    register int i,j;
    int err,m,n,iter,suc,il,ih,al,ah,l,h;
    double c,p,tmp;

    err = -1;         
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Distribution of interval-valued discrete variable. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLF = 0.001;

    if (parm(ctx, ctx->CmdBuf + 4,4,1))     /* get parameters */
        goto IDDFFin;

    if (ctx->PMNV != 2) {            /* need two variables */
        p_err(ctx, -1,1);
        goto IDDFFin;
    }
    il = ctx->PMVIdx[0];
    ih = ctx->PMVIdx[1];

    if (ctx->PMOPT != 2)
        ctx->PMOPT = 1;

    if (ctx->MxIter < 1 || ctx->MxItFlg == 0)
        ctx->MxIter = 50;

    al = (int)get_data(ctx, il,0);
    ah = (int)get_data(ctx, ih,0);

    for (i = 1; i < ctx->NOC; ++i) {
        al = imin(ctx, al,(int)get_data(ctx, il,i));  
        ah = imax(ctx, ah,(int)get_data(ctx, ih,i));  
    }
    printf1(ctx, "Range of variable: [%d,%d]\n",al,ah);

    m = ah - al + 1;

    if (alloc_acn(ctx, m + 1))       /* lower bound */ 
        goto IDDFFin;
    if (alloc_acm(ctx, m + 1))       /* upper bound */
        goto IDDFFin;
    if (alloc_acx(ctx, m + 1))       /* mean distribution */
        goto IDDFFin;
    
    for (i = 0; i < ctx->NOC; ++i) {

        l = (int)get_data(ctx, il,i) - al + 1;
        h = (int)get_data(ctx, ih,i) - al + 1;
        n = h - l + 1;

        if (l == h)
            ctx->AcN[l] += 1;

        for (j = 1; j <= m; ++j) {
            if (j >= l && j <= h) {
                ctx->AcM[j] += 1;
                ctx->AcX[j] += 1.0 / (double)n;
            }
        }
    }
    if (ctx->PMOPT != 2 || ctx->PMF1Def == 0) {
        prn_ddf(ctx, al-1,m,ctx->AcN,ctx->AcM,ctx->AcX,ctx->NOC);
        if (ctx->PMOPT != 2) {
            err = 0;
            goto IDDFFin;
        }
    }
    if (ctx->PMF1Def) {
        if (alloc_acu(ctx, m + 1))   
            goto IDDFFin;

        for (j = 1; j <= m; ++j)
            ctx->AcU[j] = ctx->AcX[j];
    }

    /* iterations for self-consistent df */

    printf1(ctx, "Iterative calculation of self-consistent distribution.\n");
    printf1(ctx, "Max number of iterations: %d\n",ctx->MxIter);
    printf1(ctx, "Tolerance for convergence: %g\n",ctx->TOLF);

    if (alloc_acy(ctx, m + 1))   
        goto IDDFFin;
    if (alloc_acz(ctx, m + 1))   
        goto IDDFFin;
    
    suc = 0;
    for (iter = 0; iter < ctx->MxIter; ++iter) {
    
        for (i = 0; i < ctx->NOC; ++i) {

            l = (int)get_data(ctx, il,i) - al + 1;
            h = (int)get_data(ctx, ih,i) - al + 1;
            n = h - l + 1;

            p = 0.0;
            for (j = l; j <= h; ++j)
                p += ctx->AcX[j];

            for (j = l; j <= h; ++j) {
                if (p > 0.0)
                    ctx->AcY[j] += ctx->AcX[j] / p;
            }
        }
        if (ctx->PMProtFDef)  
            fprintf(ctx->PMProtFd,"Iter%3d ",iter);

        c = 0.0;
        for (j = 1; j <= m; ++j) {
            tmp = ctx->AcX[j];
            ctx->AcX[j] = ctx->AcY[j] / (double)ctx->NOC;
            c = dmax(ctx, c,fabs(ctx->AcX[j] - tmp));
            ctx->AcY[j] = 0.0;
            if (ctx->PMProtFDef)  
                rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,ctx->AcX[j]);
        }
        if (ctx->PMProtFDef)  
            fprintf(ctx->PMProtFd,"\n");
    
        if (c <= ctx->TOLF) {
            suc = 1;
            break;
        }
    }
    printf1(ctx, "Convergence ");
    if (suc == 0)
        printf1(ctx, "not ");
    printf1(ctx, "reached after %d iterations.\n",iter);

    if (ctx->PMF1Def)
        prn_ddf2(ctx, al-1,m,ctx->AcN,ctx->AcM,ctx->AcU,ctx->AcX,ctx->NOC);
    else
        prn_ddf1(ctx, al-1,m,ctx->AcX);

    err = 0;

IDDFFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_ddf(al,m,low,high,mdf,noc)                                          */
/*                                                                          */  
/*  Print distributions. If PMF1Def print to PMF1d.                         */

void prn_ddf(TDAContext *ctx, int al,int m,int *low,int *high,double *mdf,int noc)
{
    FILE * fd;
    register int j;

    if (ctx->PMF1Def)
        fd = ctx->PMF1d;
    else {
        fd = TDA_CONSOLE;
        fprintf(fd,"\nValue  Lower Bound  Upper Bound      Mean DF\n");
    }
    for (j = 1; j <= m; ++j) {
        fprintf(fd,"%5d ",j + al);
        fprintf(fd,"%12.6f ",(double)low[j] / (double)noc);
        fprintf(fd,"%12.6f ",(double)high[j] / (double)noc);
        fprintf(fd,"%12.6f\n",mdf[j] / (double)noc);  
#ifdef TDA_R_PACKAGE
        {
            double erow[4];
            erow[0] = (double)(j + al);
            erow[1] = (double)low[j] / (double)noc;
            erow[2] = (double)high[j] / (double)noc;
            erow[3] = mdf[j] / (double)noc;
            tda_export_row(ctx, "imat.df", erow, 4);
        }
#endif
    }
    fprintf(fd,"\n");
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "imat.df");
#endif
    if (ctx->PMF1Def)
        printf1(ctx, "%d records written to: %s\n",m,ctx->PMF1dName);
}

/* ------------------------------------------------------------------------ */
/*  prn_ddf1(al,n,df)                                                       */
/*                                                                          */  
/*  Print distributions. If PMF1Def print to PMF1d.                         */

void prn_ddf1(TDAContext *ctx, int al,int m,double *df)
{
    register int j;

    printf1(ctx, "\nValue  Distribution\n");

    for (j = 1; j <= m; ++j) {
        printf1(ctx, "%5d  ",j + al);
        printf1(ctx, "%12.6f\n",df[j]);  
#ifdef TDA_R_PACKAGE
        {
            double erow[2];
            erow[0] = (double)(j + al);
            erow[1] = df[j];
            tda_export_row(ctx, "imat.df1", erow, 2);
        }
#endif
    }
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "imat.df1");
#endif
    newline(ctx);
}

/* ------------------------------------------------------------------------ */
/*  prn_ddf2(al,m,low,high,mdf,cdf,noc)                                     */
/*                                                                          */  
/*  Print distributions to PMF1d.                                           */

void prn_ddf2(TDAContext *ctx, int al,int m,int *low,int *high,double *mdf,double *cdf,int noc)
{
    register int j;

    for (j = 1; j <= m; ++j) {
        fprintf(ctx->PMF1d,"%5d ",j + al);
        fprintf(ctx->PMF1d,"%12.6f ",(double)low[j] / (double)noc);
        fprintf(ctx->PMF1d,"%12.6f ",(double)high[j] / (double)noc);
        fprintf(ctx->PMF1d,"%12.6f ",mdf[j] / (double)noc);  
        fprintf(ctx->PMF1d,"%12.6f\n",cdf[j]);  
#ifdef TDA_R_PACKAGE
        {
            double erow[5];
            erow[0] = (double)(j + al);
            erow[1] = (double)low[j] / (double)noc;
            erow[2] = (double)high[j] / (double)noc;
            erow[3] = mdf[j] / (double)noc;
            erow[4] = cdf[j];
            tda_export_row(ctx, "imat.df2", erow, 5);
        }
#endif
    }
    fprintf(ctx->PMF1d,"\n");
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "imat.df2");
#endif
    printf1(ctx, "%d records written to: %s\n",m,ctx->PMF1dName);
}

/* ------------------------------------------------------------------------ */
/*  imean()         Mean of interval-valued variable.                       */
/*                                                                          */  
/*                  imean(                                                  */
/*                      fmt=...,    print format 10.4                       */
/*                  ) = XL,XH;                                              */
/*                                                                          */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int imean(TDAContext *ctx)
{
    register int i;
    int err,il,ih;
    double al,ah,xl,xh,mmin,mmax;

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Mean of interval-valued variable. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,4,1))     /* get parameters */
        goto IMEANFin;

    if (ctx->PMNV != 2) {            /* need two variables */
        p_err(ctx, -1,1);
        goto IMEANFin;
    }
    il = ctx->PMVIdx[0];
    ih = ctx->PMVIdx[1];

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 10,4);

    al = get_data(ctx, il,0);
    ah = get_data(ctx, ih,0);

    mmin = mmax = 0.0;

    for (i = 0; i < ctx->NOC; ++i) {
        xl = get_data(ctx, il,i);  
        xh = get_data(ctx, ih,i);  
        al = dmin(ctx, al,xl);
        ah = dmax(ctx, ah,xh);
        mmin += xl;
        mmax += xh;
    }
    mmin /= (double)ctx->NOC;
    mmax /= (double)ctx->NOC;

    printf1(ctx, "Range of variable: [%g,%g]\n",al,ah);
    printf1(ctx, "Minimum mean value: ");
    rt_printf1_d(ctx, ctx->PMFmtS,mmin);
    newline(ctx);
    printf1(ctx, "Maximum mean value: ");
    rt_printf1_d(ctx, ctx->PMFmtS,mmax);
    newline(ctx);
    err = 0;

IMEANFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  f_range(fn,na,nb,mxit,tolbw,tolfd,tolfe,il,ih,fminp,fmaxp)              */
/*                                                                          */
/*  Calculate range of predefined functions.                                */
/*                                                                          */
/*  fn      function nubmer                                                 */
/*           1 = variance                                                   */ 
/*  na      number of function arguments                                    */
/*  nb      number of boxes                                                 */
/*  mxit    max number of iterations                                        */
/*  tolbw   tolerance for box width                                         */
/*  tolfd   tolerance for function                                          */
/*  tolfe   tolerance for final check                                       */
/*  il      index of lower bound variable                                   */
/*  ih      index of upper bound variable                                   */
/*  fminp   min function value                                              */
/*  fmaxp   max function value                                              */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int f_range(TDAContext *ctx, int fn,int na,int nb,int mxit,double tolbw,double tolfd, double tolfe,int il,int ih,double *fminp,double *fmaxp)
{
    (void)fmaxp; (void)fminp;        /* unused: the signature is shared */
    register int i;
    int ii,r,err;
    err = -1;

    /*  allocate box list for PMNBOX entries 

        GO_LB       lower bounds
        GO_UB       upper bounds
        GO_LBF      lower function values
        GO_UBF      upper function values
        GO_FLG      -1 if not yet processed
                     0 if dropped
                     1 if temporarily accepted
                     2 if finally accepted
    */

    if (alloc_list(ctx, nb,na)) 
        goto RNGFFin;

    if (alloc_par(ctx, na,1))
        goto RNGFFin;

    for (i = 0; i < na; ++i) {
        ctx->RParL[i] = get_data(ctx, il,i);  
        ctx->RParU[i] = get_data(ctx, ih,i);  
        ctx->RPar[i] = (ctx->RParL[i] + ctx->RParU[i]) / 2.0;
    }

    printf1(ctx, "Starting branch and bound algorithm.\n");
    printf1(ctx, "Maximum number of iterations: %d\n",ctx->MxIter);              
    printf1(ctx, "Maximum number of boxes: %d\n",nb);
    printf1(ctx, "Tolerance for box width: %g\n",tolbw);
    printf1(ctx, "Tolerance for function range: %g\n",tolfd);
    printf1(ctx, "Tolerance for global minimum: %g\n",tolfe);
    printf1(ctx, "Using derivatives for monotonicity test.\n");

    for (ii = 0; ii <= 1; ++ii) {
        printf1(ctx, "\nFunction ");
        if (ii == 0)
            printf1(ctx, "minimization.\n");
        else
            printf1(ctx, "maximization.\n");
  
        r = igmin(ctx, ii,na,nb,ctx->RPar,ctx->RParL,ctx->RParU,mxit,tolbw,tolfd,tolfe,1,fn);

        if (r == 1) {
            /* same honesty as ivreg: stopped at its limits -- report
               the best point found, uncertified, and still run the
               other direction */
            printf1(ctx, "Exceeded maximum number of boxes.\n");
            igmin_res(ctx, tolfe,ii);
            continue;
        }
        if (r)
            goto RNGFFin;
        igmin_res(ctx, tolfe,ii);
    }
    newline(ctx);
    err = 0;

RNGFFin:
    alloc_list(ctx, 0,0);
    alloc_par(ctx, 0,0);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  icov()          Covariance of two interval-valued variables.            */
/*  icorr()         Correlation of two interval-valued variables.          */
/*                                                                          */
/*                  icov(mxit=, nbox=, tolbw=, tolfd=, tolfe=, fmt=,        */
/*                       prot=) = XL,XU,YL,YU;                              */
/*                  icorr( ... same options ... ) = XL,XU,YL,YU;            */
/*                                                                          */
/*  Same estimand as ivar: the interval of all values the statistic can     */
/*  take when each observation varies in its interval (1/n convention).     */
/*  Not in Rohwer's TDA; added on request next to ivar, reusing its         */
/*  branch-and-bound (igmin typ 4 and 5).  Both searches use interval       */
/*  derivative monotonicity tests (quotient rule for the correlation).      */
/*                                                                          */
/*  Practical envelope, measured: covariance certifies at real sizes        */
/*  (n=11 exact against 4 million corner combinations); correlation         */
/*  certifies small problems (n=3 instantly, exact against a dense grid)    */
/*  but its clamped direct-form inclusion is too loose near the flat        */
/*  optimum to certify much beyond that -- larger runs terminate at their   */
/*  limits and say so.  A centered-form inclusion would be the next step.   */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int icov2(TDAContext *ctx, int typ)
{
    register int i;
    int err,na,r,ii;
    int iv[4];
    double fminp;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    if (typ == 4)
        printf1(ctx, "Covariance of interval-valued variables. Current memory: %d bytes.\n",ctx->MemReq);
    else
        printf1(ctx, "Correlation of interval-valued variables. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLBW = 1.e-4;
    ctx->TOLFD = 1.e-6;
    ctx->TOLFE = 1.e-6;

    if (parm(ctx, ctx->CmdBuf + (typ == 4 ? 4 : 5),4,1))
        goto ICOVFin;

    if (ctx->PMNV != 4) {            /* need four variables */
        p_err(ctx, -1,1);
        goto ICOVFin;
    }
    for (i = 0; i < 4; ++i)
        iv[i] = ctx->PMVIdx[i];

    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    if (ctx->MxIter < 1 || ctx->MxItFlg == 0)
        ctx->MxIter = 100;
    if (ctx->PMNBOX < 1)
        ctx->PMNBOX = 100;

    na = 2 * ctx->NOC;
    if (alloc_list(ctx, ctx->PMNBOX,na))
        goto ICOVFin;
    if (alloc_par(ctx, na,1))
        goto ICOVFin;

    for (i = 0; i < ctx->NOC; ++i) {
        ctx->RParL[i] = get_data(ctx, iv[0],i);
        ctx->RParU[i] = get_data(ctx, iv[1],i);
        ctx->RParL[ctx->NOC + i] = get_data(ctx, iv[2],i);
        ctx->RParU[ctx->NOC + i] = get_data(ctx, iv[3],i);
        if (ctx->RParU[i] < ctx->RParL[i] ||
            ctx->RParU[ctx->NOC + i] < ctx->RParL[ctx->NOC + i]) {
            printf1(ctx, "Error: upper bound below lower bound in case %d.\n",i + 1);
            goto ICOVFin;
        }
        ctx->RPar[i] = (ctx->RParL[i] + ctx->RParU[i]) / 2.0;
        ctx->RPar[ctx->NOC + i] =
            (ctx->RParL[ctx->NOC + i] + ctx->RParU[ctx->NOC + i]) / 2.0;
    }

    newline(ctx);
    printf1(ctx, "Starting branch and bound algorithm.\n");
    printf1(ctx, "Maximum number of iterations: %d\n",ctx->MxIter);
    printf1(ctx, "Maximum number of boxes: %d\n",ctx->PMNBOX);
    printf1(ctx, "Tolerance for box width: %g\n",ctx->TOLBW);
    printf1(ctx, "Tolerance for function range: %g\n",ctx->TOLFD);
    printf1(ctx, "Tolerance for global minimum: %g\n",ctx->TOLFE);

    for (ii = 0; ii <= 1; ++ii) {
        printf1(ctx, "\nFunction ");
        if (ii == 0)
            printf1(ctx, "minimization.\n");
        else
            printf1(ctx, "maximization.\n");

        r = igmin(ctx, ii,na,ctx->PMNBOX,ctx->RPar,ctx->RParL,ctx->RParU,ctx->MxIter,ctx->TOLBW,ctx->TOLFD,ctx->TOLFE,1,typ);
        if (r == 1) {
            /* like ivreg: the search stopped at its limits -- still
               report the best point found, honestly uncertified, and
               keep going so the other direction runs too */
            printf1(ctx, "Exceeded maximum number of boxes.\n");
            igmin_res(ctx, ctx->TOLFE,ii);
            continue;
        }
        if (r)
            goto ICOVFin;
        igmin_res(ctx, ctx->TOLFE,ii);
    }
    newline(ctx);
    err = 0;

ICOVFin:
    alloc_list(ctx, 0,0);
    alloc_par(ctx, 0,0);
    (void)fminp;
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ivar()          Variance of interval-valued variable.                   */
/*                                                                          */
/*                  ivar(                                                   */
/*                      mxit=...,       max number of iterations, def. 100  */
/*                      nbox=...,       max number of boxes, def. 100       */
/*                      tolbw=...,      tolerance box width, def. 1.e-4     */
/*                      tolfd=...,      tolerance for function range, 1e-6  */
/*                      tolfe=...,      for global min/max, def. 1.e-6      */
/*                      fmt=...,        print format, def. 10.4             */
/*                      prot=...,       protocol file                       */
/*                  ) = function;                                           */
/*                                                                          */
/*  It must be possible to evaluate the function with interval arithmetic.  */
/*  Function may contain interval operators (not completed yet).            */
/*                                                                          */  
/*  Return 0 if OK, otherwise -1.                                           */

int ivar(TDAContext *ctx)
{
    int err,il,ih;                
    double fminp,fmaxp;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Variance of interval-valued variable. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLBW = 1.e-4;       /* def tolerance for box length, branch and bound   */
    ctx->TOLFD = 1.e-6 ;      /* def tolerance for function range                 */
    ctx->TOLFE = 1.e-6 ;      /* def tolerance for global minimum                 */

    if (parm(ctx, ctx->CmdBuf + 4,4,1))     /* get parameters */
        goto IVARFin;

    if (ctx->PMNV != 2) {            /* need two variables */
        p_err(ctx, -1,1);
        goto IVARFin;
    }
    il = ctx->PMVIdx[0];
    ih = ctx->PMVIdx[1];

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 10,4);

    if (ctx->MxIter < 1 || ctx->MxItFlg == 0)
        ctx->MxIter = 100;

    if (ctx->PMNBOX < 1)         /* max number of boxes */
        ctx->PMNBOX = 100;

    newline(ctx);

    err = f_range(ctx, 1,ctx->NOC,ctx->PMNBOX,ctx->MxIter,ctx->TOLBW,ctx->TOLFD,ctx->TOLFE,il,ih,&fminp,&fmaxp);

IVARFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  igini()         Gini coefficient of interval valued functions.          */
/*                                                                          */
/*                  igini(                                                  */
/*                      x=...,          points for evaluation               */
/*                      df=...,         output file                         */
/*                      mxit=...,       max number of iterations, def. 100  */
/*                      nbox=...,       max number of boxes, def. 100       */
/*                      tolbw=...,      tolerance box width, def. 1.e-4     */
/*                      tolfd=...,      tolerance for function range, 1e-6  */
/*                      tolfe=...,      for global min/max, def. 1.e-6      */
/*                      fmt=...,        print format, def. 10.4             */
/*                      prot=...,       protocol file                       */
/*                  ) = function;                                           */
/*                                                                          */
/*  It must be possible to evaluate the function with interval arithmetic.  */
/*  Function may contain interval operators (not completed yet).            */
/*                                                                          */  
/*  Return 0 if OK, otherwise -1.                                           */

int igini(TDAContext *ctx)
{
    int err,r,il,ih,i;              
    double al,ah,x,fminp = 0.0,fmaxp = 0.0;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Gini coefficient of interval-valued variable. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLBW = 1.e-4;       /* def tolerance for box length, branch and bound   */
    ctx->TOLFD = 1.e-6 ;      /* def tolerance for function range                 */
    ctx->TOLFE = 1.e-6 ;      /* def tolerance for global minimum                 */

    if (parm(ctx, ctx->CmdBuf + 5,4,1))     /* get parameters */
        goto IGINIFin;

    if (ctx->PMNV != 2) {            /* need two variables */
        p_err(ctx, -1,1);
        goto IGINIFin;
    }
    il = ctx->PMVIdx[0];
    ih = ctx->PMVIdx[1];

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 10,4);

    if (ctx->MxIter < 1 || ctx->MxItFlg == 0)
        ctx->MxIter = 100;

    if (ctx->PMNBOX < 1)         /* max number of boxes */
        ctx->PMNBOX = 100;

    if (ctx->PMNTP < 1) {
        printf1(ctx, "Error: need x parameter.\n");
        goto IGINIFin;
    }
    newline(ctx);

    al = get_data(ctx, il,0);
    ah = get_data(ctx, ih,0);

    for (i = 0; i < ctx->NOC; ++i) {
        al = dmin(ctx, al,get_data(ctx, il,i));  
        ah = dmin(ctx, ah,get_data(ctx, ih,i));  
    }
    printf1(ctx, "Range of variable: [%g,%g]\n",al,ah);
     

    for (i = 0; i < ctx->PMNTP; ++i) {

        x = ctx->PMTP[i];

        /** tda_out(" %g\n ",x); **/
             
        if (x > al && x < ah) {

            ctx->IGMINFUNVal = x;

            r = f_range(ctx, 2,ctx->NOC,ctx->PMNBOX,ctx->MxIter,ctx->TOLBW,ctx->TOLFD,ctx->TOLFE,il,ih,&fminp,&fmaxp);
            if (r)  
                fminp = fmaxp = -1.0;
        }
        else if (x <= al)  
            fminp = fmaxp = 0.0;
        else               
            fminp = fmaxp = 1.0;
            
        if (ctx->PMF1Def) {
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,x);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,fminp);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,fmaxp);
            fprintf(ctx->PMF1d,"\n");
        }
    }
    err = 0;

IGINIFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ivar1()     New algorithm for calculating the variance of an interval-  */
/*              valued variable.                                            */
/*                                                                          */  
/*                  ivar1(                                                  */
/*                      fmt=...,        print format, def. 0.0              */  
/*                  ) = XL,XU;                                              */
/*                                                                          */
/*  There must be exactly two variables on the right-hand side to be        */  
/*  interpreted as lower and upper bounds of an interval-valued variable.   */  
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int ivar1(TDAContext *ctx)
{
    register int i,j;
    int err,il,ih,n,nn,fin;     
    double xl,xh,tmp,mmin,mmax,umin,vmin,vmax,a,b,u,v,mmin1,mmax1;
         
    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Variance of interval-valued variable. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,4,1))     /* get parameters */
        goto IVFin;

    if (ctx->PMNV != 2) {            /* need two variables */
        p_err(ctx, -1,1);
        goto IVFin;
    }
    il = ctx->PMVIdx[0];
    ih = ctx->PMVIdx[1];

    if (alloc_acx(ctx, ctx->NOC + 1))     /* lower bound */ 
        goto IVFin;
    if (alloc_acy(ctx, ctx->NOC + 1))     /* upper bound */
        goto IVFin;
    if (alloc_aci(ctx, ctx->NOC + 1))     /* index */
        goto IVFin;
    if (alloc_actmp(ctx, 2 * ctx->NOC + 1))     /* induced partition */
        goto IVFin;

    n = 0;
    mmin = mmax = 0.0;

    for (i = 0; i < ctx->NOC; ++i) {
        xl = get_data(ctx, il,i);
        xh = get_data(ctx, ih,i);
        if (xh <= xl + ctx->EPSI1) {
            printf1(ctx, "Error in case %d, found: %g,%g\n",i+1,xl,xh);
            goto IVFin;
        }
        ctx->AcX[i] = xl;
        ctx->AcY[i] = xh;

        /* printf1("i=%d x=%g y=%g \n",i,AcX[i],AcY[i]); */
         
        ctx->AcTmp[n++] = xl;
        ctx->AcTmp[n++] = xh;

        mmin += xl;
        mmax += xh;
    }
    mmin /= (double)ctx->NOC;
    mmax /= (double)ctx->NOC;

    if (sortd(ctx, n,ctx->AcTmp,0))
        goto IVFin;
   
    nn = 1;
    for (i = 1; i < n; ++i) {           /* AcTmp contains induced partition */
        if (ctx->AcTmp[i] > ctx->AcTmp[i - 1])
            ctx->AcTmp[nn++] = ctx->AcTmp[i];
    }

    vmin = vmax = 0.0;
    for (i = 0; i < ctx->NOC; ++i) {
        u = mmin - ctx->AcX[i];
        vmin += u * u;
        u = mmax - ctx->AcY[i];
        vmax += u * u;
    }
    vmin /= (double)ctx->NOC;
    vmax /= (double)ctx->NOC;

    printf1(ctx, "\nRange of mean values:      ");
    rt_printf1_d(ctx, ctx->PMFmtS,mmin);
    printf1(ctx, "\n                           ");
    rt_printf1_d(ctx, ctx->PMFmtS,mmax);
    printf1(ctx, "\nVariance at lower bounds:  ");
    rt_printf1_d(ctx, ctx->PMFmtS,vmin);
    printf1(ctx, "\nVariance at upper bounds:  ");
    rt_printf1_d(ctx, ctx->PMFmtS,vmax);
    newline(ctx);

    umin = vmin = ctx->DBLMAX;
    fin = 0;
    for (i = 1; i < nn; ++i) {
        a = ctx->AcTmp[i - 1];
        b = ctx->AcTmp[i];
        if (b < mmin)
            continue;
        if (a > mmax)
            break;
            
        n = 0;
        u = 0.0;
        for (j = 0; j < ctx->NOC; ++j) {
            xl = ctx->AcX[j];
            xh = ctx->AcY[j];
            if (xh <= a) {
                u += xh;
                n++;
            }
            else if (xl >= b) {
                u += xl;
                n++;
            }
        }
        if (n > 0) {
            u /= (double)n;  
            tmp = ivarf(ctx, ctx->NOC,u,ctx->AcX,ctx->AcY);
            tmp /= (double)ctx->NOC;
            if (tmp < vmin) {
                vmin = tmp;
                umin = u;             
                fin = 1;
            }
        }
    }
    if (fin == 0) {
        printf1(ctx, "Cannot find unique minimum.\n");
    }   
    else {
        printf1(ctx, "\nMinimum value of variance: ");
        rt_printf1_d(ctx, ctx->PMFmtS,vmin);
        printf1(ctx, "\nCorresponding mean value:  ");
        rt_printf1_d(ctx, ctx->PMFmtS,umin);
        newline(ctx);
    }

  
    /** tda_out("mmin=%f mmax=%f\n",mmin,mmax); **/
  

    mmin1 = mmin;
    mmax1 = mmax;
    fin = 0;

    for (i = 1; i < 10; ++i) {
        n = 0;
        mmin = mmax = 0.0;

        for (j = 0; j < ctx->NOC; ++j) {
            if (ctx->AcI[j] == 0) {
                n++;

                xl = ctx->AcX[j];
                xh = ctx->AcY[j];
                tmp = (xl + xh) / 2.0;
                if (tmp <= mmin1) {
                    ctx->AcI[j] = 1;
                    mmin += xl;
                    mmax += xl;
                }
                else if (tmp >= mmax1) {
                    ctx->AcX[j] = xh;
                    ctx->AcI[j] = 1;
                    mmin += xh;
                    mmax += xh;
                }
                else {
                    mmin += xl;
                    mmax += xh;
                }           
            }
            else {
                mmin += ctx->AcX[j];
                mmax += ctx->AcX[j];
            }
        }
        mmin /= (double)ctx->NOC;
        mmax /= (double)ctx->NOC;
  
        /** tda_out("mmin=%f mmax=%f n=%d \n",mmin,mmax,n ); **/
  
        if (n == 0) {
            fin = 1;
            break;
        }
        mmin1 = mmin;
        mmax1 = mmax;
    }
    if (fin == 0) {
        printf1(ctx, "\nCannot calculate maximum of variance.\n");
    }
    else {
        u = 0.0;
        for (j = 0; j < ctx->NOC; ++j)
            u += ctx->AcX[j];
        u /= (double)ctx->NOC;

        v = 0.0;
        for (j = 0; j < ctx->NOC; ++j)  
            v += (ctx->AcX[j] - u) * (ctx->AcX[j] - u);
        v /= (double)ctx->NOC;

        printf1(ctx, "\nMaximum value of variance: ");
        rt_printf1_d(ctx, ctx->PMFmtS,v);
        printf1(ctx, "\nCorresponding mean value:  ");
        rt_printf1_d(ctx, ctx->PMFmtS,u);
        newline(ctx);
    }
    newline(ctx);
    err = 0;

IVFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ivarf()                                                                 */
/*                                                                          */

double ivarf(TDAContext *ctx, int n,double x,double *xl,double *xh)
{
    (void)ctx;        /* unused: the signature is shared */
    register int j;
    double d,v;

    v = 0.0;
    for (j = 0; j < n; ++j) {
        if (xl[j] > x)
            d = xl[j] - x;
        else if (xh[j] < x)
            d = x - xh[j];
        else  
            d = 0;
        v += d * d;
    }
    return(v);
}

