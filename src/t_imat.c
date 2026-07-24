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

/* ------------------------------------------------------------------------ */
/*  functions in t_imat                                                     */

int m_midf(char *cmd);
int m_midf1(char *cmd);
int m_midf2(char *cmd);
int get_ipart(int n,double *xl,double *xu,double *ip);
int m_midf3(char *cmd);
int gmin(void);
int range(void);
int alloc_par(int n,int opt);
int igmin(int opt,int narg,int nbmax,double *par,double *lb,double *ub,int mxit,
    double tolbw,double tolfd,double tolfe,int gc,int typ);
void igmin_cpar(int n,double *x,double *par);
int alloc_list(int n,int m);
int igmin_res(double tolfe,int opt);
void i_add(double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_sub(double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_mul(double xl,double xh,double yl,double yh,double *rl,double *rh);
int  i_div(double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_abs(double xl,double xh,double *rl,double *rh);
void i_max(double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_min(double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_square(double xl,double xh,double *rl,double *rh);
void i_sqrt(double xl,double xh,double *rl,double *rh);

int igmin_fun(int opt,int typ,int n,double *xl,double *xh,
    double *rl,double *rh,int deriv,double *gl,double *gh);
int idf(void);
int sddf(void);
int iddf(void);
void prn_ddf(int al,int m,int *low,int *high,double *mdf,int noc);
void prn_ddf1(int al,int m,double *df);
void prn_ddf2(int al,int m,int *low,int *high,double *mdf,double *cdf,int noc); 
int imean(void);
int f_range(int fn,int na,int nb,int mxit,double tolbw,double tolfd,
    double tolfe,int il,int ih,double *fminp,double *fmaxp);
int ivar(void);
int igini(void);
int ivar1(void);
double ivarf(int n,double x,double *xl,double *xh);


/*  Global parameters and variables                                         */

double TOLBW = 1.e-4;       /* tolerance for box length, branch and bound   */
double TOLBC = 1.e-4;       /* tolerance for constraint function boxes      */
double TOLFD = 1.e-10;      /* tolerance for function range                 */
double TOLFE = 1.e-10;      /* tolerance for global minimum                 */

double *RPar;               /* parameters                                   */
double *RParL;              /* lower bounds                                 */
double *RParU;              /* upper bounds                                 */
int RParA = 0;              /* allocation                                   */
int RParLA = 0;             /* allocation                                   */
int RParUA = 0;             /* allocation                                   */

double *GO_LB;              /* box list                                     */
double *GO_UB; 
double *GO_LBF; 
double *GO_UBF; 
short  *GO_FLG;
int GO_LBA = 0;
int GO_UBA = 0;
int GO_LBFA = 0;
int GO_UBFA = 0;
int GO_FLGA = 0;

int GO_IT = 0;              /* number of iterations performed               */
int GO_FN = 0;              /* number of function evaluations               */
int GO_IFN = 0;             /* number of inclusion function evaluations     */
int GO_NB = 0;              /* number of boxes in final list                */
int GO_NA = 0;              /* number of accepted boxes                     */
int GO_NBU = 0;             /* number of boxes used                         */
double GO_FMIN = 0.0;       /* currently best global function value         */
double IGMINFUNVal = 0.0;   /* argument in igmin_fun                        */

/*--------------------------------------------------------------------------*/
/*  m_midf(cmd)     midf(XL,XU,DL,DU,DM)        matrix command              */
/*                                                                          */
/*                  Calculates distribution function DM, and lower (DL)     */
/*                  and upper (DU) bounds for an interval vector (XL,XU).   */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_midf(char *cmd)
{
    register int i,j;
    int err,n,nn,row,col,idx,idx1,idx2,idx3,ivflg;
    register char *p;
    double xl,xu,xx,tmp,f,fl,fu;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    if ((p = n_mexpr(cmd + 5,2,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MIDFFin;
    if ((p = m_check1(p)) == NULL)
        goto MIDFFin;
    if (col != 1) {
        mat_err(7);
        goto MIDFFin;
    }
    if (alloc_actmp(2 * row))
        goto MIDFFin;

    n = 0;
    for (i = 1; i <= row; ++i) {
        xl = MX[0][i];
        xu = MX[1][i];
        if (fabs(xl - xu) <= EPSI1) {
            mat_err(11);
            goto MIDFFin;
        }
        if (xl > xu) {
            tmp = xl;
            xl = xu;
            xu = tmp;
            MX[0][i] = xl;
            MX[1][i] = xu;
        }
        AcTmp[n++] = xl;
        AcTmp[n++] = xu;
    }
    if (sortd(n,AcTmp,0))
        goto MIDFFin;

    nn = 1;
    for (i = 1; i < n; ++i) {
        if (AcTmp[i] > AcTmp[i - 1])
            AcTmp[nn++] = AcTmp[i];
    }
    if ((p = m_getmat(p,nn,1,&idx,cmd,0)) == NULL)
        goto MIDFFin;
    if ((p = m_check1(p)) == NULL)
        goto MIDFFin;
    if ((p = m_getmat(p,nn,1,&idx1,cmd,0)) == NULL)
        goto MIDFFin;
    if ((p = m_check1(p)) == NULL)
        goto MIDFFin;
    if ((p = m_getmat(p,nn,1,&idx2,cmd,0)) == NULL)
        goto MIDFFin;
    if ((p = m_check1(p)) == NULL)
        goto MIDFFin;
    if ((p = m_getmat(p,nn,1,&idx3,cmd,1)) == NULL)
        goto MIDFFin;

    /* save values of induced partition */

    for (i = 1; i <= nn; ++i)
        MatVal[idx][i] = AcTmp[i - 1];

    if (sortd2(row,MX[0] + 1,MX[1] + 1))
        goto MIDFFin;

    /* calculate distribution functions */

    for (i = 1; i <= nn; ++i) {
        tmp = AcTmp[i - 1];
        fl = fu = f = 0.0;
        for (j = 1; j <= row; ++j) {
            xl = MX[0][j];
            xu = MX[1][j];
            if (xu <= tmp)
                fl += 1.0;
            if (xl <= tmp) {
                fu += 1.0;
                if (xl < tmp) {
                    xx = dmin(xu,tmp);
                    f += (xx - xl) / (xu - xl);
                }
            }
            else
                break;
        }
        MatVal[idx1][i] = fl / (double)row;
        MatVal[idx2][i] = fu / (double)row;
        MatVal[idx3][i] = f  / (double)row;
    }
    err = 0;

MIDFFin:
    alloc_actmp(0);
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_midf1(cmd)    midf(XL,XU,F)           matrix command                  */
/*                                                                          */
/*                  Calculates mean distribution function F for lower end   */
/*                  points of [XL,XU].                                      */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_midf1(char *cmd)
{
    register int i,j,jj;
    int err,row,col,idx,ivflg;
    register char *p;
    double xl,xu,xx,tmp,f;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    if ((p = n_mexpr(cmd + 6,2,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MIDF1Fin;
    if ((p = m_check1(p)) == NULL)
        goto MIDF1Fin;
    if (col != 1) {
        mat_err(7);
        goto MIDF1Fin;
    }
    if ((p = m_getmat(p,row,1,&idx,cmd,1)) == NULL)
        goto MIDF1Fin;

    if (alloc_acn(row))
        goto MIDF1Fin;

    if (sortdp(row,MX[0] + 1,AcN))
        goto MIDF1Fin;

    for (i = 1; i <= row; ++i) {
        tmp = MX[0][i];
        f = 0.0;
        for (j = 0; j < row; ++j) {
            jj = AcN[j] + 1;
            xl = MX[0][jj];
            xu = MX[1][jj];
            if (xl < tmp) {
                xx = dmin(xu,tmp);
                f += (xx - xl) / (xu - xl);
            }
            else
                break;
        }
        MatVal[idx][i] = f  / (double)row;
    }
    err = 0;

MIDF1Fin:
    alloc_acn(0);
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_midf2(cmd)    midf(XL,XU,F)                 matrix command            */
/*                                                                          */
/*                  Calculates mean distribution function F for lower end   */
/*                  points of [XL,XU].                                      */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_midf2(char *cmd)
{
    register int i,j,jj;
    int err,nn,row,col,idx,ivflg;
    register char *p;
    double xl,xu,xx,tmp,f;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    if ((p = n_mexpr(cmd + 6,2,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MIDF2Fin;
    if ((p = m_check1(p)) == NULL)
        goto MIDF2Fin;
    if (col != 1) {
        mat_err(7);
        goto MIDF2Fin;
    }
    if ((p = m_getmat(p,row,1,&idx,cmd,1)) == NULL)
        goto MIDF2Fin;

    if (alloc_actmp(2 * row))
        goto MIDF2Fin;

    if ((nn = get_ipart(row,MX[0],MX[1],AcTmp)) < 2)
        goto MIDF2Fin;
   
    if (alloc_acn(row))
        goto MIDF2Fin;

    if (sortdp(row,MX[0] + 1,AcN))
        goto MIDF2Fin;
     
    if (alloc_acu(nn))          /* used for distr function */
        goto MIDF2Fin;

    for (i = 0; i < nn; ++i) {
        tmp = AcTmp[i];              
        f = 0.0;
        for (j = 0; j < row; ++j) {
            jj = AcN[j] + 1;
            xl = MX[0][jj];
            xu = MX[1][jj];
            if (xl < tmp) {
                xx = dmin(xu,tmp);
                f += (xx - xl) / (xu - xl);
            }
            else
                break;
        }
        AcU[i] = f  / (double)row;
    }
    for (i = 1; i <= row; ++i) {
        tmp = MX[0][i];
        for (j = 0; j < nn; ++j) {
            if (fabs(tmp - AcTmp[j]) < EPSI1)
                break;
        }
        tmp = 0.0;
        for (jj = j; jj < nn - 1; ++jj) {
            xl = AcTmp[jj];
            xu = AcTmp[jj + 1];
            f  = (AcU[jj + 1] - AcU[jj]) / (xu - xl);
            tmp += f * (xu * xu - xl * xl);
        }
        tmp /= 2.0;
        if (1.0 - AcU[j] <= EPSI) 
            goto MIDF2Fin;

        MatVal[idx][i] = tmp / (1.0 - AcU[j]);
    }
    err = 0;

MIDF2Fin:
    alloc_acn(0);
    alloc_actmp(0);
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  get_ipart(n,xl,xu,ip)                                                   */
/*                                                                          */
/*  Calculate induced partition for [xl,xu] (i=1,...,n). Return in          */
/*  ip[j], j = 0,...,m-1. Return m, or -1 if error.                         */

int get_ipart(int n,double *xl,double *xu,double *ip)
{
    register int i;
    int m,nn;
    double xxl,xxu;

    m = 0;
    for (i = 1; i <= n; ++i) {
        xxl = xl[i];
        xxu = xu[i];
        if (xxl >= xxu - EPSI1) {
            mat_err(11);
            return(-1);     
        }
        ip[m++] = xxl;
        ip[m++] = xxu;
    }
    if (sortd(m,ip,0))
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

int m_midf3(char *cmd)
{
    register int i,j;
    int err,n1,n2,row,col,idx,idx1,ivflg;
    register char *p;
    double xl,xu,xl1,xu1,tmp1,tmp2;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    if ((p = n_mexpr(cmd + 6,2,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MIDF3Fin;
    if ((p = m_check1(p)) == NULL)
        goto MIDF3Fin;
    if (col != 1) {
        mat_err(7);
        goto MIDF3Fin;
    }
    if ((p = m_getmat(p,row,1,&idx,cmd,0)) == NULL)
        goto MIDF3Fin;
    if ((p = m_check1(p)) == NULL)
        goto MIDF3Fin;
    if ((p = m_getmat(p,row,1,&idx1,cmd,1)) == NULL)
        goto MIDF3Fin;

    for (i = 1; i <= row; ++i) {
        xl = MX[0][i];
        xu = MX[1][i];
        n1 = n2 = 0;
        tmp1 = tmp2 = 0.0;
        for (j = 1; j <= row; ++j) {
            xl1 = MX[0][j];
            xu1 = MX[1][j];
            if (xl1 >= xl && xl1 <= xu) {
                tmp1 += xl1;
                n1++;
            }
            if (xu1 >= xl && xu1 <= xu) {
                tmp2 += xu1;
                n2++;
            }
            if (n1 > 0)
                MatVal[idx][i] = tmp1 / (double)n1;
            if (n2 > 0)
                MatVal[idx1][i] = tmp2 / (double)n2;
        }
    }
    err = 0;

MIDF3Fin:
    mx_free();
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

int gmin(void)
{
    register int i,j,k;
    int err,r,n,l,gc;          
    double w,*uptr,*lptr;

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("Global minimization. Current memory: %d bytes.\n",MemReq);

    TOLBW = 1.e-4;       /* def tolerance for box length, branch and bound   */
    TOLFD = 1.e-10;      /* def tolerance for function range                 */
    TOLFE = 1.e-10;      /* def tolerance for global minimum                 */

    if (parm(CmdBuf + 4,12,1))   /* get parameters */
        goto GMINFin;

    if (PMFmtF == 0)            /* default print format */
        pmfmt(13,6);
           
    newline();

    if (FNArgN < 1) {
        printf1("Error: function should contain at least one argument.\n");
        goto GMINFin;
    }
    if (alloc_par(FNArgN,1))
        goto GMINFin;

    if (get_dsv(FNArgN,RPar-1,1,RParL-1,RParU-1,1))   /* get starting values */ 
        goto GMINFin;

    /* print parameters and starting values */
     
    if (prn_fsval(FNArgN,RPar,1,RParL,RParU,0))
        goto GMINFin;

    if (MxIter < 1 || MxItFlg == 0)
        MxIter = 100;

    if (PMNS == 1)
        gc = 1;
    else
        gc = 0;

    if (PMNBOX < 1)         /* max number of boxes */
        PMNBOX = 100;

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

    if (alloc_list(PMNBOX,FNArgN)) 
        goto GMINFin;

    printf1("Starting branch and bound algorithm.\n");
    printf1("Maximum number of iterations: %d\n",MxIter);              
    printf1("Maximum number of boxes: %d\n",PMNBOX);
    printf1("Tolerance for box width: %g\n",TOLBW);
    printf1("Tolerance for function range: %g\n",TOLFD);
    printf1("Tolerance for global minimum: %g\n",TOLFE);
    if (gc) {
        printf1("Using derivatives for monotonicity test.\n");

        if (fnd_alloc(1,1,FNArgN,FNPN,1))
            goto GMINFin;
    }
    newline();

    r = igmin(0,FNArgN,PMNBOX,RPar,RParL,RParU,MxIter,TOLBW,TOLFD,TOLFE,gc,0);
    if (r) {
        if (r == 1)
            printf1("Exceeded maximum number of boxes.\n");
        goto GMINFin;
    }
    n = igmin_res(TOLFE,0);      
    newline();
              
    if (n > 0) {

        printf1("Box  Acc  Width             Lower function bound  Upper function bound\n");
        r = 1;
        for (i = 0; i < GO_NB; ++i) {
            if (GO_FLG[i] == 2) {

                lptr = GO_LB + i * FNArgN;
                uptr = GO_UB + i * FNArgN;

                w = uptr[0] - lptr[0];
                for (j = 1; j < FNArgN; ++j)  
                    w = dmax(w,uptr[j] - lptr[j]);

                printf1("%3d  %3d %17.10e %21.14e %21.14e\n",
                                    r++,GO_FLG[i],w,GO_LBF[i],GO_UBF[i]);
            }
        }

        l = 10;
        if (l < FNMLen)
            l = FNMLen;

        r = 1;
        for (i = 0; i < GO_NB; ++i) {
            if (GO_FLG[i] == 2) {
                printf1("\nBox  Idx  Parameter\n");

                lptr = GO_LB + i * FNArgN;
                uptr = GO_UB + i * FNArgN;

                for (j = 0; j < FNArgN; ++j) {
                    k = FNArgSP[j];
                    printf1("%3d  %3d  %s   ",r,j + 1,FNArgDef[k]);
                    prnchar(' ',l - strlen(FNArgDef[k]),0);
                    printf1(PMFmtS,lptr[j]);
                    printf1(PMFmtS,uptr[j]);
                    newline();
                }
                r++;
            }
        }
    }
    newline();
    err = 0;

GMINFin:
    alloc_list(0,0);
    alloc_par(0,0);
    fnd_alloc(0,0,0,0,0);
    p_clean();
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

int range(void)
{
    int err,r,ii,gc;              

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("Range of interval-valued function. Current memory: %d bytes.\n",MemReq);

    TOLBW = 1.e-4;       /* def tolerance for box length, branch and bound   */
    TOLFD = 1.e-10;      /* def tolerance for function range                 */
    TOLFE = 1.e-10;      /* def tolerance for global minimum                 */

    if (parm(CmdBuf + 5,12,1))   /* get parameters */
        goto RNGFin;

    if (PMFmtF == 0)            /* default print format */
        pmfmt(13,6);
           
    newline();

    if (FNArgN < 1) {
        printf1("Error: function should contain at least one argument.\n");
        goto RNGFin;
    }
    if (alloc_par(FNArgN,1))
        goto RNGFin;

    if (get_dsv(FNArgN,RPar-1,1,RParL-1,RParU-1,1))   /* get starting values */ 
        goto RNGFin;
    
    /* print parameters and starting values */
     
    if (prn_fsval(FNArgN,RPar,1,RParL,RParU,0))
        goto RNGFin;

    if (PMNS == 1) {
        gc = 1;
        if (fnd_alloc(1,1,FNArgN,FNPN,1))
            goto RNGFin;
    }
    else
        gc = 0;

    if (MxIter < 1 || MxItFlg == 0)
        MxIter = 100;

    if (PMNBOX < 1)         /* max number of boxes */
        PMNBOX = 100;

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

    if (alloc_list(PMNBOX,FNArgN)) 
        goto RNGFin;

    printf1("Starting branch and bound algorithm.\n");
    printf1("Maximum number of iterations: %d\n",MxIter);              
    printf1("Maximum number of boxes: %d\n",PMNBOX);
    printf1("Tolerance for box width: %g\n",TOLBW);
    printf1("Tolerance for function range: %g\n",TOLFD);
    printf1("Tolerance for global minimum: %g\n",TOLFE);
    if (gc)
        printf1("Using derivatives for monotonicity test.\n");

    for (ii = 0; ii <= 1; ++ii) {
        printf1("\nFunction ");
        if (ii == 0)
            printf1("minimization.\n");
        else
            printf1("maximization.\n");

        r = igmin(ii,FNArgN,PMNBOX,RPar,RParL,RParU,MxIter,TOLBW,TOLFD,TOLFE,gc,0);

        if (r) {
            if (r == 1)
                printf1("Exceeded maximum number of boxes.\n");
            goto RNGFin;
        }
        igmin_res(TOLFE,ii);
    }
    newline();
    err = 0;

RNGFin:
    alloc_list(0,0);
    alloc_par(0,0);
    fnd_alloc(0,0,0,0,0);
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  alloc_par(n,opt)    If opt != 0 allocate memory for n parameters and    */
/*                      bounds. Otherwise free previously allocated memory  */
/*                                                                          */
/*          Return 0 if OK, otherwise -1.                                   */

int alloc_par(int n,int opt)
{
    int err = 0;

    if (opt == 0)
        goto APARFin;

    err = -1;
    if (!(RPar = (double *)calloc(n,sizeof(double)))) { 
        p_err(-2,1);
        goto APARFin;
    }   
    memrq(n,sizeof(double));
    RParA = n;             
              
    if (!(RParL = (double *)calloc(n,sizeof(double)))) { 
        p_err(-2,1);
        goto APARFin;
    }   
    memrq(n,sizeof(double));
    RParLA = n;             
              
    if (!(RParU = (double *)calloc(n,sizeof(double)))) { 
        p_err(-2,1);
        goto APARFin;
    }   
    memrq(n,sizeof(double));
    RParUA = n;             
    return(0);

APARFin:
    if (RParA > 0) {
        free((char *)RPar);
        memrq(-RParA,sizeof(double));
        RParA = 0;
    }
    if (RParLA > 0) {
        free((char *)RParL);
        memrq(-RParLA,sizeof(double));
        RParLA = 0;
    }
    if (RParUA > 0) {
        free((char *)RParU);
        memrq(-RParUA,sizeof(double));
        RParUA = 0;
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

int igmin(int opt,int narg,int nbmax,double *par,double *lb,double *ub,int mxit,
    double tolbw,double tolfd,double tolfe,int gc,int typ)
{
    register int i,j;
    int err,r,nb,is,iter,ja,jb,jj,first;
    double fl,fu,tmp,w,f,fminp;
    double *lptra,*lptrb,*uptra,*uptrb;

    GO_IT = 0;              /* number of iterations performed */
    GO_FN = 0;              /* number of function evaluations */
    GO_IFN = 0;             /* number of inclusion function evaluations */
    GO_NBU = 0;             /* number of boxes used */
    GO_NB = 0;              /* number of boxes in final list */

    err = -1;

    for (i = 0; i < nbmax; ++i)
        GO_FLG[i] = 0;

    if (alloc_actmp(narg + 1))  
        goto IGMINFin; 
    if (gc) {
        if (alloc_acu(narg + 1))    /* lower bounds of gradient */
            goto IGMINFin; 
        if (alloc_acv(narg + 1))    /* upper bounds of gradient */
            goto IGMINFin; 
    }
    err = 0;

    lptra = GO_LB;                      /* put initial box on list */
    uptra = GO_UB;

    for (i = 0; i < narg; ++i) {
        lptra[i] = lb[i];
        uptra[i] = ub[i];
    }
    GO_FLG[0] = -1;     /* not yet processed */

    /* calculate inclusion function with initial boxes */
   
    GO_IFN++;

    if (typ == 0) {
        r = get_ifval(&fl,&fu,narg,lptra,uptra,0,AcU,AcV);
        if (r) {
            printf1("Error in evaluating inclusion function.\n");
            prn_emsg2(r);              
            err = -2;
            goto IGMINFin;
        }
    }
    else {
        r = igmin_fun(1,typ,narg,lptra,uptra,&fl,&fu,0,AcU,AcV);                                      
        if (r) {
            printf1("Error in evaluating inclusion function (typ %d).\n",typ);
            err = -2;
            goto IGMINFin;
        }
    }
    GO_LBF[0] = fl;
    GO_UBF[0] = fu;
    
    /* calculate function at starting values */

    GO_FN++;
    if (typ == 0) {
        r = get_flval(&fminp,narg,par,0,0,&tmp,&tmp,&tmp);
        if (r) {      /* r from v_eval1() */
            printf1("Error in function evaluation with starting values.\n");
            prn_emsg2(r);              
            err = -2;
            goto IGMINFin;
        }
    }
    else {
        r = igmin_fun(0,typ,narg,par,&tmp,&fminp,&tmp,0,AcU,AcV);                                      
        if (r) {
            printf1("Error in evaluating function (typ %d).\n",typ);
            err = -2;
            goto IGMINFin;
        }
    }
    nb = 1;     /* number of boxes */
    iter = 0;                                     

    if (SILENTFlg < 2)
        printfe("\n  Iter    Function Value       NBox  FCall  IFCall\n");
    
    while (++iter <= mxit) {       

        if (SILENTFlg < 2)
            printfe("%5d  %20.13e  %6d %6d %7d\n",iter,fminp,nb,GO_FN,GO_IFN);

        /* find box with lowest lower, or largest upper, function bound,
           and check for boxes that can be dropped */

        ja = -1;
        first = 1;
        for (i = 0; i < nb; ++i) {

            if (GO_FLG[i]) {    /* check all boxes */
  
                if (opt == 0) {
                    if (GO_LBF[i] > fminp)  
                        GO_FLG[i] = 0;
 
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
                else {

                    if (GO_UBF[i] < fminp)  
                        GO_FLG[i] = 0;
                    else if (GO_FLG[i] < 0) {
                        if (first) {      
                            f = GO_UBF[i];
                            ja = i;
                            first = 0;
                        }                 
                        else if (f < GO_UBF[i]) {
                            f = GO_UBF[i];
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
            if (GO_FLG[i] == 0) {
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

        if (w <= TOLBW) {                   /* temporarily accept this box */
            GO_FLG[ja] = 1;
            continue;
        }

        for (i = 0; i < narg; ++i) {      /* copy box ja to jb */
            lptrb[i] = lptra[i];
            uptrb[i] = uptra[i];
        }

        tmp = lptra[is] + w / 2.0;
        uptra[is] = lptrb[is] = tmp;
        GO_FLG[jb] = -1;

        for (jj = 0; jj < 2; ++jj) {           /* for both boxes */

            if (jj)   
                ja = jb;

            lptra = GO_LB + ja * narg;
            uptra = GO_UB + ja * narg;

            /* evaluate inclusion function */

            GO_IFN++;
            if (typ == 0) {
                r = get_ifval(&fl,&fu,narg,lptra,uptra,gc,AcU,AcV);
                if (r) {
                    printf1("Error in evaluating inclusion function.\n");
                    prn_emsg2(r);              
                    err = -2;
                    goto IGMINFin;
                }
            }
            else {
                r = igmin_fun(1,typ,narg,lptra,uptra,&fl,&fu,gc,AcU,AcV);                                      
                if (r) {
                    printf1("Error in evaluating inclusion function (typ %d).\n",typ);
                    err = -2;
                    goto IGMINFin;
                }
            }
            if (opt == 0) {
                if (fl > fminp) {
                    GO_FLG[ja] = 0;             /* free box */
                    if (ja == nb - 1)
                        nb--;    
                    continue;
                }
            }
            else {
                if (fu < fminp) {
                    GO_FLG[ja] = 0;             /* free box */
                    if (ja == nb - 1)
                        nb--;    
                    continue;
                }
            }
            GO_LBF[ja] = fl;
            GO_UBF[ja] = fu;

            if (opt == 0) {
                if (fminp > fu) {
                    fminp = fu;
                    igmin_cpar(narg,uptra,par);
                }
            }
            else {
                if (fminp < fl) {
                    fminp = fl;
                    igmin_cpar(narg,lptra,par);
                }
            }
            if (fabs(fu - fl) < tolfd)      /* accept for solution */
                GO_FLG[ja] = 1;
   
            if (gc) {                       /* gradient check */

                r = 0;
                for (j = 0; j < narg; ++j) {
                    if (opt == 0) {
                        if (AcU[j] > 0.0) {
                            if (lptra[j] > lb[j]) {
                                r = 1;
                                break;
                            }
                            uptra[j] = lptra[j];
                            r = -1;
                        }
                        else if (AcV[j] < 0.0) {
                            if (uptra[j] < ub[j]) {
                                r = 1;
                                break;
                            }
                            lptra[j] = uptra[j];
                            r = -1;
                        }
                    }
                    else {
                        if (AcU[j] > 0.0) {
                            if (uptra[j] < ub[j]) {
                                r = 1;
                                break;
                            }
                            lptra[j] = uptra[j];
                            r = -1;
                        }
                        else if (AcV[j] < 0.0) {
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
                    GO_FLG[ja] = 0;             /* free box */
                    if (ja == nb - 1)
                        nb--;    
                    continue;
                }
                else if (r == -1) {         /* update inclusion function */
                    GO_IFN++;
                    if (typ == 0) {
                        r = get_ifval(&fl,&fu,narg,lptra,uptra,0,&tmp,&tmp);
                        if (r) {
                            printf1("Error in evaluating inclusion function.\n");
                            prn_emsg2(r);              
                            err = -2;
                            goto IGMINFin;
                        }
                    }
                    else {
                        r = igmin_fun(1,typ,narg,lptra,uptra,&fl,&fu,0,AcU,AcV);                                      
                        if (r) {
                            printf1("Error in evaluating inclusion function (typ %d).\n",typ);
                            err = -2;
                            goto IGMINFin;
                        }
                    }

                    if (opt == 0) {
                        if (fl > fminp) {
                            GO_FLG[ja] = 0;             /* free box */
                            if (ja == nb - 1)
                                nb--;    
                            continue;
                        }
                    }
                    else {
                        if (fu < fminp) {
                            GO_FLG[ja] = 0;             /* free box */
                            if (ja == nb - 1)
                                nb--;    
                            continue;
                        }
                    }
                    GO_LBF[ja] = fl;
                    GO_UBF[ja] = fu;
                }
            }

            /* update fminp with function value for midpoint of a box with
               least lower bound */
                   
            for (j = 0; j < narg; ++j)
                AcTmp[j] = (lptra[j] + uptra[j]) / 2.0;

            /* calculate function at midpoint */

            GO_FN++;
            if (typ == 0) {
                r = get_flval(&f,narg,AcTmp,0,0,&tmp,&tmp,&tmp);
                if (r) {      /* r from v_eval1() */
                    printf1("Error in function evaluation.\n");
                    prn_emsg2(r);              
                    err = -2;
                    goto IGMINFin;
                }
            }
            else {
                r = igmin_fun(0,typ,narg,AcTmp,&tmp,&f,&tmp,0,AcU,AcV);                                      
                if (r) {
                    printf1("Error in evaluating function (typ %d).\n",typ);
                    err = -2;
                    goto IGMINFin;
                }
            }

            if (opt == 0) {
                if (fminp > f) {
                    fminp = f;
                    igmin_cpar(narg,AcTmp,par);
                }
            }   
            else {
                if (fminp < f) {
                    fminp = f;
                    igmin_cpar(narg,AcTmp,par);
                }
            }   
        }
    }
    GO_NB = nb;
    GO_IT = iter;
    GO_FMIN = fminp;

    if (PMProtFDef) {
        fprintf(PMProtFd,"IGMIN optimization.\n");
        prval("Best function value",GO_FMIN);
        prvec("Parameters",narg,par - 1);
        fprintf(PMProtFd,"\n");
        
        if (opt == 0)
            tmp = GO_FMIN - tolfe;
        else
            tmp = GO_FMIN + tolfe;

        ja = 1;
        for (i = 0; i < nb; ++i) {
            if (GO_FLG[i] == 0)  
                continue;

            if (GO_FLG[i] >= 1) {

                if (opt == 0) {                 /* minimization */
                    if (GO_LBF[i] >= tmp)         
                        GO_FLG[i] = 2;
                }
                else {
                    if (GO_UBF[i] <= tmp)         
                        GO_FLG[i] = 2;
                }
            }
            lptra = GO_LB + i * narg;
            uptra = GO_UB + i * narg;

            w = uptra[0] - lptra[0];
            for (j = 1; j < narg; ++j)  
                w = dmax(w,uptra[j] - lptra[j]);

            fprintf(PMProtFd,"Box   Width                 Lower function bound  Upper function bound  Acceptance\n");
            fprintf(PMProtFd,"%3d  %21.14e %21.14e %21.14e  %6d\n",
                                 ja++,w,GO_LBF[i],GO_UBF[i],GO_FLG[i]);

            fprintf(PMProtFd,"Parameter vector\nLB: ");
            for (j = 0; j < narg; ++j)  
                fprintf(PMProtFd,PMPFmtS,lptra[j]);
            fprintf(PMProtFd,"\nUB: ");
            for (j = 0; j < narg; ++j)  
                fprintf(PMProtFd,PMPFmtS,uptra[j]);
            fprintf(PMProtFd,"\n\n");
        }
    }

IGMINFin:
    alloc_acu(0);           
    alloc_acv(0);           
    alloc_actmp(0);           
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  igmin_cpar(n,x,par) copy parameters from x[] into par[].                */

void igmin_cpar(int n,double *x,double *par)
{
    register int i;

    for (i = 0; i < n; ++i)
        par[i] = x[i];
}

/* ------------------------------------------------------------------------ */
/*  alloc_list(n,m)     Allocate box list with n entries, m parameters.     */
/*                      If n == 0 free previously allocated memory.         */
/*                      Return 0 if OK, -1 if error.                        */

int alloc_list(int n,int m)
{
    int err,nm;

    if (n <= 0) {
        err = 0;          
        goto AListFin;
    }
    err = -1;
    nm = n * m;
    if (!(GO_LB = (double *)calloc(nm,sizeof(double)))) { 
        p_err(-2,1);
        goto AListFin;
    }
    memrq(nm,sizeof(double));
    GO_LBA = nm;

    if (!(GO_UB = (double *)calloc(nm,sizeof(double)))) { 
        p_err(-2,1);
        goto AListFin;
    }
    memrq(nm,sizeof(double));
    GO_UBA = nm;

    if (!(GO_LBF = (double *)calloc(n,sizeof(double)))) { 
        p_err(-2,1);
        goto AListFin;
    }
    memrq(n,sizeof(double));
    GO_LBFA = n;

    if (!(GO_UBF = (double *)calloc(n,sizeof(double)))) { 
        p_err(-2,1);
        goto AListFin;
    }
    memrq(n,sizeof(double));
    GO_UBFA = n;

    if (!(GO_FLG = (short *)calloc(n,sizeof(short)))) { 
        p_err(-2,1);
        goto AListFin;
    }
    memrq(n,sizeof(short));
    GO_FLGA = n;

    return(0);

AListFin:
    if (GO_LBA > 0) {
        free((char *)GO_LB);
        memrq(-GO_LBA,sizeof(double));
        GO_LBA = 0;
    }
    if (GO_UBA > 0) {
        free((char *)GO_UB);
        memrq(-GO_UBA,sizeof(double));
        GO_UBA = 0;
    }
    if (GO_LBFA > 0) {
        free((char *)GO_LBF);
        memrq(-GO_LBFA,sizeof(double));
        GO_LBFA = 0;
    }
    if (GO_UBFA > 0) {
        free((char *)GO_UBF);
        memrq(-GO_UBFA,sizeof(double));
        GO_UBFA = 0;
    }
    if (GO_FLGA > 0) {
        free((char *)GO_FLG);
        memrq(-GO_FLGA,sizeof(short));
        GO_FLGA = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  igmin_res(tolfe,opt)                                                    */
/*                                                                          */
/*  Get results from igmin(). If opt = 0 minimum, otherwise maximum         */
/*  Return number of finally accepted boxes.                                */

int igmin_res(double tolfe,int opt)
{
    register int i;
    int n0,n1,n2;
    double tmp,a;

    n0 = n1 = n2 = 0;         
    if (opt == 0)
        tmp = GO_FMIN - tolfe;
    else
        tmp = GO_FMIN + tolfe;
    a = tmp;

    for (i = 0; i < GO_NB; ++i) {
        if (GO_FLG[i] >= 1) {
            n1++;
            if (opt == 0) {
                if (GO_LBF[i] >= tmp) {       
                    n2++;
                    GO_FLG[i] = 2;
                }
                a = dmin(a,GO_LBF[i]);
            }
            else {
                if (GO_UBF[i] <= tmp) {       
                    n2++;
                    GO_FLG[i] = 2;
                }
                a = dmax(a,GO_UBF[i]);
            }
        }
        else if (GO_FLG[i] < 0)
            n0++;
    }

    printf1("Number of iterations performed: %d\n",GO_IT);
    printf1("Number of function evaluations: %d\n",GO_FN);
    printf1("Number of inclusion function evaluations: %d\n",GO_IFN);
    printf1("Number of boxes used: %d\n",GO_NBU);
    printf1("Number of temporarily accepted boxes: %d\n",n1);
    printf1("Number of finally accepted boxes: %d\n",n2);
    if (n0 > 0)
        printf1("Warning: %d boxes have not been processed.\n",n0);

    if (opt == 0)  
        printf1("\nBest minimal function value: ");
    else
        printf1("\nBest maximal function value: ");
    printf1(PMFmtS,GO_FMIN);
    if (opt == 0)  
        printf1("  best lower bound: ");
    else
        printf1("  best upper bound: ");
    if (n1 == 0)
        printf1("***");
    else
        printf1(PMFmtS,a);
    newline();
    return(n2);
}

/* ------------------------------------------------------------------------ */
/*  i_add(xl,xh,yl,yh,rl,rh)                                                */
/*                                                                          */
/*  [rl,rh] = [xl,xh] + [yl,yh]                                             */

void i_add(double xl,double xh,double yl,double yh,double *rl,double *rh)
{
    *rl = xl + yl;
    *rh = xh + yh;
}

/* ------------------------------------------------------------------------ */
/*  i_sub(xl,xh,yl,yh,rl,rh)                                                */
/*                                                                          */
/*  [rl,rh] = [xl,xh] - [yl,yh]                                             */

void i_sub(double xl,double xh,double yl,double yh,double *rl,double *rh)
{
    *rl = xl - yh;
    *rh = xh - yl;
}

/* ------------------------------------------------------------------------ */
/*  i_mul(xl,xh,yl,yh,rl,rh)                                                */
/*                                                                          */
/*  [rl,rh] = [xl,xh] * [yl,yh]                                             */

void i_mul(double xl,double xh,double yl,double yh,double *rl,double *rh)
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
            *rl = dmin(xl * yh,xh * yl);
            *rh = dmax(xl * yl,xh * yh);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  i_div(xl,xh,yl,yh,rl,rh)                                                */
/*                                                                          */
/*  [rl,rh] = [xl,xh] / [yl,yh]                                             */
/*                                                                          */
/*  Return 0 if OK, -1 if error                                             */

int i_div(double xl,double xh,double yl,double yh,double *rl,double *rh)
{
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

void i_abs(double xl,double xh,double *rl,double *rh)
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
        *rh = dmax(-xl,xh);
    }
}

/* ------------------------------------------------------------------------ */
/*  i_max(xl,xh,yl,yh,rl,rh)                                                */
/*                                                                          */

void i_max(double xl,double xh,double yl,double yh,double *rl,double *rh)
{

    *rl = dmax(xl,yl);
    *rh = dmax(xh,yh);
}

/* ------------------------------------------------------------------------ */
/*  i_min(xl,xh,yl,yh,rl,rh)                                                */
/*                                                                          */

void i_min(double xl,double xh,double yl,double yh,double *rl,double *rh)
{

    *rl = dmin(xl,yl);
    *rh = dmin(xh,yh);
}

/* ------------------------------------------------------------------------ */
/*  i_square(xl,xh,rl,rh)                                                   */
/*                                                                          */

void i_square(double xl,double xh,double *rl,double *rh)
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
        *rh = dmax(xl * xl,xh * xh);
    }
}

/* ------------------------------------------------------------------------ */
/*  i_sqrt(xl,xh,rl,rh)                                                     */
/*                                                                          */
/*  [rl,rh] = [sqrt(xl),sqrt(xh)]                                           */

void i_sqrt(double xl,double xh,double *rl,double *rh)
{
    if (xl < 0.0) {
        fprintf(stderr,"ERROR in i_sqrt.\n");
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

int igmin_fun(int opt,int typ,int n,double *xl,double *xh,
    double *rl,double *rh,int deriv,double *gl,double *gh)
{
    register int i,j,k;
    double x,y,x1,y1,tl,th,sl,sh,tmp;
/**   
    for (i = 0; i < n; ++i)
        printf1("i=%d xl=%lf xh=%lf\n",i,xl[i],xh[i]);
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
                i_sub(xl[i],xh[i],x,y,&tl,&th);
                i_mul(tl,th,tl,th,&sl,&sh);
                i_add(*rl,*rh,sl,sh,&tl,&th);
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
                    if (xl[i] <= IGMINFUNVal)
                        y += xl[i];
                }
                *rl = y / x;
                return(0);
            }
            *rl = *rh = x = y = x1 = y1 = 0.0;
            for (i = 0; i < n; ++i) {
                x += xl[i];
                y += xh[i];
                if (xl[i] <= IGMINFUNVal)
                    x1 += xl[i];
                if (xh[i] >= IGMINFUNVal)
                    y1 += xh[i];
            }
            i_div(x1,y1,x,y,rl,rh);

            if (deriv) {
                for (i = 0; i < n; ++i) {
                    gl[i] = gh[i] = -1.0;
                }
            }
            return(0);

        case 3:                 /* multi-dimensional scaling */

            if (opt == 0) {
                *rl = 0.0;
                for (i = 1; i < GD_NP; ++i) {
                    for (j = 0; j < i; ++j) {
                        x = gdd_adj(i,j,PMGN);
                        if (x >= 0.0) {
                            tmp = x - fabs(xl[i] - xl[j]);
                            *rl += tmp * tmp;
                        }
                    }
                }
                return(0);
            }
            *rl = *rh = 0.0;

            for (i = 1; i < GD_NP; ++i) {
                for (j = 0; j < i; ++j) {
                    x = gdd_adj(i,j,PMGN);
                    if (x >= 0.0) {
                        i_sub(xl[i],xh[i],xl[j],xh[j],&tl,&th);
                        i_abs(tl,th,&sl,&sh);
                        i_sub(x,x,sl,sh,&tl,&th);
                        i_mul(tl,th,tl,th,&sl,&sh);
                        i_add(*rl,*rh,sl,sh,&tl,&th);
                        *rl = tl;
                        *rh = th;
                    }
                }
            }
            if (deriv) {
                for (k = 0; k < n; ++k) {
                    gl[k] = gh[k] = 0;
                    for (i = 1; i < GD_NP; ++i) {
                        for (j = 0; j < i; ++j) {
                            if (i == k || j == k) {
                                x = gdd_adj(i,j,PMGN);
                                if (x >= 0.0) {
                                    i_sub(xl[i],xh[i],xl[j],xh[j],&tl,&th);
                                    i_abs(tl,th,&sl,&sh);
                                    i_sub(x,x,sl,sh,&tl,&th);

                                    if (k == i) {
                                        if (xh[i] + EPSI1 < xl[j])        
                                            i_mul(tl,th,2,2,&sl,&sh);
                                        else if (xh[j] + EPSI1 < xl[i])   
                                            i_mul(tl,th,-2,-2,&sl,&sh);
                                        else
                                            i_mul(tl,th,-2,2,&sl,&sh);
                                    }
                                    else {
                                        if (xh[j] + EPSI1 < xl[i])        
                                            i_mul(tl,th,2,2,&sl,&sh);
                                        else if (xh[i] + EPSI1 < xl[j])   
                                            i_mul(tl,th,-2,-2,&sl,&sh);
                                        else
                                            i_mul(tl,th,-2,2,&sl,&sh);
                                    }
                                    tl = gl[k];
                                    th = gh[k];
                                    i_add(tl,th,sl,sh,&x,&y);
                                    gl[k] = x;
                                    gh[k] = y;
                                }
                            }
                        }
                    }
                }
            }
            return(0);
          
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

int idf(void)
{
    register int i,j,k;
    int err,il,ih,n,nn,ai,bi,iter;
    double xx,xl,xh,tmp,f,fl,fh;
         
    err = -1;         
    if (check_cmd(1))
        return(-1);

    printf1("Interval-valued distribution. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 3,4,1))     /* get parameters */
        goto IDFFin;

    if (PMOPT != 2)
        PMOPT = 1;

    if (PMNV != 2) {            /* need two variables */
        p_err(-1,1);
        goto IDFFin;
    }
    il = PMVIdx[0];
    ih = PMVIdx[1];

    if (alloc_acx(NOC + 1))     /* lower bound */ 
        goto IDFFin;
    if (alloc_acy(NOC + 1))     /* upper bound */
        goto IDFFin;
    if (alloc_actmp(2 * NOC + 1))     /* induced partition */
        goto IDFFin;

    n = 0;
    for (i = 0; i < NOC; ++i) {
        xl = get_data(il,i);
        xh = get_data(ih,i);
        if (xh <= xl + EPSI1) {
            printf1("Error in case %d, found: %g,%g\n",i+1,xl,xh);
            goto IDFFin;
        }
        AcX[i] = xl;
        AcY[i] = xh;

        /* printf1("i=%d x=%g y=%g \n",i,AcX[i],AcY[i]); */
         
        AcTmp[n++] = xl;
        AcTmp[n++] = xh;
    }
    if (sortd(n,AcTmp,0))
        goto IDFFin;
   
    nn = 1;
    for (i = 1; i < n; ++i) {           /* AcTmp contains induced partition */
        if (AcTmp[i] > AcTmp[i - 1])
            AcTmp[nn++] = AcTmp[i];
    }
    if (sortd2(NOC,AcX,AcY))            /* sort observations */
        goto IDFFin;

    /* calculate distribution functions */

    if (alloc_acz(nn))     /* mean df */
        goto IDFFin;

    printf1("\n  Idx      Partition  Lower Bound  Upper Bound      Mean DF\n");

    for (i = 0; i < nn; ++i) {
        tmp = AcTmp[i];
        fl = fh = f = 0.0;
        for (j = 0; j < NOC; ++j) {
            xl = AcX[j];
            xh = AcY[j];
            if (xh <= tmp)
                fl += 1.0;
            if (xl <= tmp) {
                fh += 1.0;
                if (xl < tmp) {
                    xx = dmin(xh,tmp);
                    f += (xx - xl) / (xh - xl);
                }
            }
            else
                break;
        }
        AcZ[i] = f / (double)NOC;

        printf1("%5d ",i + 1);
        printf1("%14.6lf ",tmp);
        printf1("%12.6lf ",fl / (double)NOC);
        printf1("%12.6lf ",fh / (double)NOC);
        printf1("%12.6lf\n",f  / (double)NOC);

    }
    newline();

    if (PMOPT != 2) {
        err = 0;
        goto IDFFin;
    }

    /* iterations for self-consistent df */

    if (MxIter < 1 || MxItFlg == 0)
        MxIter = 50;

    printf1("Iterative calculation of self-consistent distribution.\n");
    printf1("Max number of iterations: %d\n",MxIter);
    printf1("Tolerance for convergence: %g\n",TOLF);

    if (alloc_acu(nn)) 
        goto IDFFin;

    for (iter = 0; iter < MxIter; ++iter) {


        for (i = 0; i < nn; ++i) {

            xx = AcTmp[i];

/*   printf1("i=%d xx=%g\n",i,xx);
*/

            tmp = 0.0;

            for (j = 0; j < NOC; ++j) {
                xl = AcX[j];
                xh = AcY[j];
                ai = bi = -1;

                if (xl > xx)  
                    ;
                else if (xh < xx)
                    tmp += 1.0;
                else {
                    ai = bi = -1;
                    for (k = 0; k < nn; ++k) {
                        if (fabs(xl - AcTmp[k]) < 0.01)  
                            ai = k;
                        if (fabs(xh - AcTmp[k]) < 0.01)  
                            bi = k;
                    }
                    if (ai < 0 || bi < 0) {
                        printf2("ERROR\n");
                        gerr_exit(223);
                    }   
                    tmp += (AcZ[i] - AcZ[ai]) / (AcZ[bi] - AcZ[ai]);
                }
                /****
                printf1("j=%d xl=%lf xh=%lf tmp=%lf ai=%2d bi=%2d\n",
                            j,xl,xh,tmp,ai,bi);
                ***/

            }
            tmp /= (double)NOC;
            AcU[i] = tmp;
            /***  printf1("tmp=%lf\n",tmp); ***/
        }
        /***
        printf1("acz: ");
        for (i = 0; i < nn; ++i)
            printf1("%lf ",AcZ[i]);
        newline();

        printf1("acu: ");
        for (i = 0; i < nn; ++i)    {
            printf1("%lf ",AcU[i]);

            AcZ[i] = AcU[i];
        }
        newline();
        ***/
    }
    err = 0;

IDFFin:
    p_clean();
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

int sddf(void)
{
    register int i,j,k;
    int err,m,n,iter,suc;
    double c,p,tmp;

    err = -1;         
    if (check_cmd(2))
        return(-1);

    printf1("Distribution of set-valued discrete variable. Current memory: %d bytes.\n",MemReq);

    TOLF = 0.001;

    if (parm(CmdBuf + 4,4,1))     /* get parameters */
        goto SDDFFin;

    if (PMOPT != 2)
        PMOPT = 1;

    if (MxIter < 1 || MxItFlg == 0)
        MxIter = 50;

    m = PMNV;                   /* number of variables = categories */

    if (alloc_acn(m + 1))       /* lower bound */ 
        goto SDDFFin;
    if (alloc_acm(m + 1))       /* upper bound */
        goto SDDFFin;
    if (alloc_acx(m + 1))       /* mean distribution */
        goto SDDFFin;
    if (alloc_ack(m + 1))                           
        goto SDDFFin;


    for (i = 0; i < NOC; ++i) {
        k = -1;
        n = 0;
        for (j = 1; j <= m; ++j) {
            if (get_data(PMVIdx[j-1],i)) {
                k = j;
                n++;
                AcM[k] += 1;
                AcK[k] = 1;
            }
            else
                AcK[j] = 0;
        }
        if (n == 1)  
            AcN[k] += 1;
        for (j = 1; j <= m; ++j)
            AcX[j] += (double)AcK[j] / (double)n;
    }
    if (PMOPT != 2 || PMF1Def == 0) {
        prn_ddf(0,m,AcN,AcM,AcX,NOC);
        if (PMOPT != 2) {
            err = 0;
            goto SDDFFin;
        }
    }
    if (PMF1Def) {
        if (alloc_acu(m + 1))   
            goto SDDFFin;

        for (j = 1; j <= m; ++j)
            AcU[j] = AcX[j];
    }

    /* iterations for self-consistent df */

    printf1("Iterative calculation of self-consistent distribution.\n");
    printf1("Max number of iterations: %d\n",MxIter);
    printf1("Tolerance for convergence: %g\n",TOLF);

    if (alloc_acy(m + 1))   
        goto SDDFFin;
    if (alloc_acz(m + 1))   
        goto SDDFFin;
    
    suc = 0;
    for (iter = 0; iter < MxIter; ++iter) {
    
        for (i = 0; i < NOC; ++i) {
            p = 0.0;
            for (j = 1; j <= m; ++j) {
                if (get_data(PMVIdx[j-1],i)) {
                    p += AcX[j];
                    AcZ[j] = AcX[j];
                }
                else
                    AcZ[j] = 0.0;
            }
            for (j = 1; j <= m; ++j) {
                if (p > 0.0)
                    AcY[j] += AcZ[j] / p;
            }
        }
        if (PMProtFDef)  
            fprintf(PMProtFd,"Iter%3d ",iter);

        c = 0.0;
        for (j = 1; j <= m; ++j) {
            tmp = AcX[j];
            AcX[j] = AcY[j] / (double)NOC;
            c = dmax(c,fabs(AcX[j] - tmp));
            AcY[j] = 0.0;
            if (PMProtFDef)  
                fprintf(PMProtFd,PMPFmtS,AcX[j]);
        }
        if (PMProtFDef)  
            fprintf(PMProtFd,"\n");
    
        if (c <= TOLF) {
            suc = 1;
            break;
        }
    }
    printf1("Convergence ");
    if (suc == 0)
        printf1("not ");
    printf1("reached after %d iterations.\n",iter);

    if (PMF1Def)
        prn_ddf2(0,m,AcN,AcM,AcU,AcX,NOC);
    else
        prn_ddf1(0,m,AcX);

    err = 0;

SDDFFin:
    p_clean();
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

int iddf(void)
{
    register int i,j;
    int err,m,n,iter,suc,il,ih,al,ah,l,h;
    double c,p,tmp;

    err = -1;         
    if (check_cmd(2))
        return(-1);

    printf1("Distribution of interval-valued discrete variable. Current memory: %d bytes.\n",MemReq);

    TOLF = 0.001;

    if (parm(CmdBuf + 4,4,1))     /* get parameters */
        goto IDDFFin;

    if (PMNV != 2) {            /* need two variables */
        p_err(-1,1);
        goto IDDFFin;
    }
    il = PMVIdx[0];
    ih = PMVIdx[1];

    if (PMOPT != 2)
        PMOPT = 1;

    if (MxIter < 1 || MxItFlg == 0)
        MxIter = 50;

    al = (int)get_data(il,0);
    ah = (int)get_data(ih,0);

    for (i = 1; i < NOC; ++i) {
        al = imin(al,(int)get_data(il,i));  
        ah = imax(ah,(int)get_data(ih,i));  
    }
    printf1("Range of variable: [%d,%d]\n",al,ah);

    m = ah - al + 1;

    if (alloc_acn(m + 1))       /* lower bound */ 
        goto IDDFFin;
    if (alloc_acm(m + 1))       /* upper bound */
        goto IDDFFin;
    if (alloc_acx(m + 1))       /* mean distribution */
        goto IDDFFin;
    
    for (i = 0; i < NOC; ++i) {

        l = (int)get_data(il,i) - al + 1;
        h = (int)get_data(ih,i) - al + 1;
        n = h - l + 1;

        if (l == h)
            AcN[l] += 1;

        for (j = 1; j <= m; ++j) {
            if (j >= l && j <= h) {
                AcM[j] += 1;
                AcX[j] += 1.0 / (double)n;
            }
        }
    }
    if (PMOPT != 2 || PMF1Def == 0) {
        prn_ddf(al-1,m,AcN,AcM,AcX,NOC);
        if (PMOPT != 2) {
            err = 0;
            goto IDDFFin;
        }
    }
    if (PMF1Def) {
        if (alloc_acu(m + 1))   
            goto IDDFFin;

        for (j = 1; j <= m; ++j)
            AcU[j] = AcX[j];
    }

    /* iterations for self-consistent df */

    printf1("Iterative calculation of self-consistent distribution.\n");
    printf1("Max number of iterations: %d\n",MxIter);
    printf1("Tolerance for convergence: %g\n",TOLF);

    if (alloc_acy(m + 1))   
        goto IDDFFin;
    if (alloc_acz(m + 1))   
        goto IDDFFin;
    
    suc = 0;
    for (iter = 0; iter < MxIter; ++iter) {
    
        for (i = 0; i < NOC; ++i) {

            l = (int)get_data(il,i) - al + 1;
            h = (int)get_data(ih,i) - al + 1;
            n = h - l + 1;

            p = 0.0;
            for (j = l; j <= h; ++j)
                p += AcX[j];

            for (j = l; j <= h; ++j) {
                if (p > 0.0)
                    AcY[j] += AcX[j] / p;
            }
        }
        if (PMProtFDef)  
            fprintf(PMProtFd,"Iter%3d ",iter);

        c = 0.0;
        for (j = 1; j <= m; ++j) {
            tmp = AcX[j];
            AcX[j] = AcY[j] / (double)NOC;
            c = dmax(c,fabs(AcX[j] - tmp));
            AcY[j] = 0.0;
            if (PMProtFDef)  
                fprintf(PMProtFd,PMPFmtS,AcX[j]);
        }
        if (PMProtFDef)  
            fprintf(PMProtFd,"\n");
    
        if (c <= TOLF) {
            suc = 1;
            break;
        }
    }
    printf1("Convergence ");
    if (suc == 0)
        printf1("not ");
    printf1("reached after %d iterations.\n",iter);

    if (PMF1Def)
        prn_ddf2(al-1,m,AcN,AcM,AcU,AcX,NOC);
    else
        prn_ddf1(al-1,m,AcX);

    err = 0;

IDDFFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_ddf(al,m,low,high,mdf,noc)                                          */
/*                                                                          */  
/*  Print distributions. If PMF1Def print to PMF1d.                         */

void prn_ddf(int al,int m,int *low,int *high,double *mdf,int noc)
{
    FILE * fd;
    register int j;

    if (PMF1Def)
        fd = PMF1d;
    else {
        fd = stdout;
        fprintf(fd,"\nValue  Lower Bound  Upper Bound      Mean DF\n");
    }
    for (j = 1; j <= m; ++j) {
        fprintf(fd,"%5d ",j + al);
        fprintf(fd,"%12.6f ",(double)low[j] / (double)noc);
        fprintf(fd,"%12.6f ",(double)high[j] / (double)noc);
        fprintf(fd,"%12.6f\n",mdf[j] / (double)noc);  
    }
    fprintf(fd,"\n");
    if (PMF1Def)
        printf1("%d records written to: %s\n",m,PMF1dName);
}

/* ------------------------------------------------------------------------ */
/*  prn_ddf1(al,n,df)                                                       */
/*                                                                          */  
/*  Print distributions. If PMF1Def print to PMF1d.                         */

void prn_ddf1(int al,int m,double *df)
{
    register int j;

    printf1("\nValue  Distribution\n");

    for (j = 1; j <= m; ++j) {
        printf1("%5d  ",j + al);
        printf1("%12.6f\n",df[j]);  
    }
    newline();
}

/* ------------------------------------------------------------------------ */
/*  prn_ddf2(al,m,low,high,mdf,cdf,noc)                                     */
/*                                                                          */  
/*  Print distributions to PMF1d.                                           */

void prn_ddf2(int al,int m,int *low,int *high,double *mdf,double *cdf,int noc)
{
    register int j;

    for (j = 1; j <= m; ++j) {
        fprintf(PMF1d,"%5d ",j + al);
        fprintf(PMF1d,"%12.6f ",(double)low[j] / (double)noc);
        fprintf(PMF1d,"%12.6f ",(double)high[j] / (double)noc);
        fprintf(PMF1d,"%12.6f ",mdf[j] / (double)noc);  
        fprintf(PMF1d,"%12.6f\n",cdf[j]);  
    }
    fprintf(PMF1d,"\n");
    printf1("%d records written to: %s\n",m,PMF1dName);
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

int imean(void)
{
    register int i;
    int err,il,ih;
    double al,ah,xl,xh,mmin,mmax;

    err = -1;         
    if (check_cmd(1))
        return(-1);

    printf1("Mean of interval-valued variable. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,4,1))     /* get parameters */
        goto IMEANFin;

    if (PMNV != 2) {            /* need two variables */
        p_err(-1,1);
        goto IMEANFin;
    }
    il = PMVIdx[0];
    ih = PMVIdx[1];

    if (PMFmtF == 0)            /* default print format */
        pmfmt(10,4);

    al = get_data(il,0);
    ah = get_data(ih,0);

    mmin = mmax = 0.0;

    for (i = 0; i < NOC; ++i) {
        xl = get_data(il,i);  
        xh = get_data(ih,i);  
        al = dmin(al,xl);
        ah = dmax(ah,xh);
        mmin += xl;
        mmax += xh;
    }
    mmin /= (double)NOC;
    mmax /= (double)NOC;

    printf1("Range of variable: [%g,%g]\n",al,ah);
    printf1("Minimum mean value: ");
    printf1(PMFmtS,mmin);
    newline();
    printf1("Maximum mean value: ");
    printf1(PMFmtS,mmax);
    newline();
    err = 0;

IMEANFin:
    p_clean();
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

int f_range(int fn,int na,int nb,int mxit,double tolbw,double tolfd,
    double tolfe,int il,int ih,double *fminp,double *fmaxp)
{
    register int i;
    int ii,r,err;

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

    if (alloc_list(nb,na)) 
        goto RNGFFin;

    if (alloc_par(na,1))
        goto RNGFFin;

    for (i = 0; i < na; ++i) {
        RParL[i] = get_data(il,i);  
        RParU[i] = get_data(ih,i);  
        RPar[i] = (RParL[i] + RParU[i]) / 2.0;
    }

    printf1("Starting branch and bound algorithm.\n");
    printf1("Maximum number of iterations: %d\n",MxIter);              
    printf1("Maximum number of boxes: %d\n",nb);
    printf1("Tolerance for box width: %g\n",tolbw);
    printf1("Tolerance for function range: %g\n",tolfd);
    printf1("Tolerance for global minimum: %g\n",tolfe);
    printf1("Using derivatives for monotonicity test.\n");

    for (ii = 0; ii <= 1; ++ii) {
        printf1("\nFunction ");
        if (ii == 0)
            printf1("minimization.\n");
        else
            printf1("maximization.\n");
  
        r = igmin(ii,na,nb,RPar,RParL,RParU,mxit,tolbw,tolfd,tolfe,1,fn);

        if (r) {
            if (r == 1)
                printf1("Exceeded maximum number of boxes.\n");
            goto RNGFFin;
        }
        igmin_res(tolfe,ii);
    }
    newline();
    err = 0;

RNGFFin:
    alloc_list(0,0);
    alloc_par(0,0);
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

int ivar(void)
{
    int err,il,ih;                
    double fminp,fmaxp;

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("Variance of interval-valued variable. Current memory: %d bytes.\n",MemReq);

    TOLBW = 1.e-4;       /* def tolerance for box length, branch and bound   */
    TOLFD = 1.e-6 ;      /* def tolerance for function range                 */
    TOLFE = 1.e-6 ;      /* def tolerance for global minimum                 */

    if (parm(CmdBuf + 4,4,1))     /* get parameters */
        goto IVARFin;

    if (PMNV != 2) {            /* need two variables */
        p_err(-1,1);
        goto IVARFin;
    }
    il = PMVIdx[0];
    ih = PMVIdx[1];

    if (PMFmtF == 0)            /* default print format */
        pmfmt(10,4);

    if (MxIter < 1 || MxItFlg == 0)
        MxIter = 100;

    if (PMNBOX < 1)         /* max number of boxes */
        PMNBOX = 100;

    newline();

    err = f_range(1,NOC,PMNBOX,MxIter,TOLBW,TOLFD,TOLFE,il,ih,&fminp,&fmaxp);

IVARFin:
    p_clean();
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

int igini(void)
{
    int err,r,il,ih,i;              
    double al,ah,x,fminp,fmaxp;

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("Gini coefficient of interval-valued variable. Current memory: %d bytes.\n",MemReq);

    TOLBW = 1.e-4;       /* def tolerance for box length, branch and bound   */
    TOLFD = 1.e-6 ;      /* def tolerance for function range                 */
    TOLFE = 1.e-6 ;      /* def tolerance for global minimum                 */

    if (parm(CmdBuf + 5,4,1))     /* get parameters */
        goto IGINIFin;

    if (PMNV != 2) {            /* need two variables */
        p_err(-1,1);
        goto IGINIFin;
    }
    il = PMVIdx[0];
    ih = PMVIdx[1];

    if (PMFmtF == 0)            /* default print format */
        pmfmt(10,4);

    if (MxIter < 1 || MxItFlg == 0)
        MxIter = 100;

    if (PMNBOX < 1)         /* max number of boxes */
        PMNBOX = 100;

    if (PMNTP < 1) {
        printf1("Error: need x parameter.\n");
        goto IGINIFin;
    }
    newline();

    al = get_data(il,0);
    ah = get_data(ih,0);

    for (i = 0; i < NOC; ++i) {
        al = dmin(al,get_data(il,i));  
        ah = dmin(ah,get_data(ih,i));  
    }
    printf1("Range of variable: [%g,%g]\n",al,ah);
     

    for (i = 0; i < PMNTP; ++i) {

        x = PMTP[i];

        /** printf(" %g\n ",x); **/
             
        if (x > al && x < ah) {

            IGMINFUNVal = x;

            r = f_range(2,NOC,PMNBOX,MxIter,TOLBW,TOLFD,TOLFE,il,ih,&fminp,&fmaxp);
            if (r)  
                fminp = fmaxp = -1.0;
        }
        else if (x <= al)  
            fminp = fmaxp = 0.0;
        else               
            fminp = fmaxp = 1.0;
            
        if (PMF1Def) {
            fprintf(PMF1d,PMFmtS,x);
            fprintf(PMF1d,PMFmtS,fminp);
            fprintf(PMF1d,PMFmtS,fmaxp);
            fprintf(PMF1d,"\n");
        }
    }
    err = 0;

IGINIFin:
    p_clean();
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

int ivar1(void)
{
    register int i,j;
    int err,il,ih,n,nn,fin;     
    double xl,xh,tmp,mmin,mmax,umin,vmin,vmax,a,b,u,v,mmin1,mmax1;
         
    err = -1;         
    if (check_cmd(1))
        return(-1);

    printf1("Variance of interval-valued variable. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,4,1))     /* get parameters */
        goto IVFin;

    if (PMNV != 2) {            /* need two variables */
        p_err(-1,1);
        goto IVFin;
    }
    il = PMVIdx[0];
    ih = PMVIdx[1];

    if (alloc_acx(NOC + 1))     /* lower bound */ 
        goto IVFin;
    if (alloc_acy(NOC + 1))     /* upper bound */
        goto IVFin;
    if (alloc_aci(NOC + 1))     /* index */
        goto IVFin;
    if (alloc_actmp(2 * NOC + 1))     /* induced partition */
        goto IVFin;

    n = 0;
    mmin = mmax = 0.0;

    for (i = 0; i < NOC; ++i) {
        xl = get_data(il,i);
        xh = get_data(ih,i);
        if (xh <= xl + EPSI1) {
            printf1("Error in case %d, found: %g,%g\n",i+1,xl,xh);
            goto IVFin;
        }
        AcX[i] = xl;
        AcY[i] = xh;

        /* printf1("i=%d x=%g y=%g \n",i,AcX[i],AcY[i]); */
         
        AcTmp[n++] = xl;
        AcTmp[n++] = xh;

        mmin += xl;
        mmax += xh;
    }
    mmin /= (double)NOC;
    mmax /= (double)NOC;

    if (sortd(n,AcTmp,0))
        goto IVFin;
   
    nn = 1;
    for (i = 1; i < n; ++i) {           /* AcTmp contains induced partition */
        if (AcTmp[i] > AcTmp[i - 1])
            AcTmp[nn++] = AcTmp[i];
    }

    vmin = vmax = 0.0;
    for (i = 0; i < NOC; ++i) {
        u = mmin - AcX[i];
        vmin += u * u;
        u = mmax - AcY[i];
        vmax += u * u;
    }
    vmin /= (double)NOC;
    vmax /= (double)NOC;

    printf1("\nRange of mean values:      ");
    printf1(PMFmtS,mmin);
    printf1("\n                           ");
    printf1(PMFmtS,mmax);
    printf1("\nVariance at lower bounds:  ");
    printf1(PMFmtS,vmin);
    printf1("\nVariance at upper bounds:  ");
    printf1(PMFmtS,vmax);
    newline();

    umin = vmin = DBLMAX;
    fin = 0;
    for (i = 1; i < nn; ++i) {
        a = AcTmp[i - 1];
        b = AcTmp[i];
        if (b < mmin)
            continue;
        if (a > mmax)
            break;
            
        n = 0;
        u = 0.0;
        for (j = 0; j < NOC; ++j) {
            xl = AcX[j];
            xh = AcY[j];
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
            tmp = ivarf(NOC,u,AcX,AcY);
            tmp /= (double)NOC;
            if (tmp < vmin) {
                vmin = tmp;
                umin = u;             
                fin = 1;
            }
        }
    }
    if (fin == 0) {
        printf1("Cannot find unique minimum.\n");
    }   
    else {
        printf1("\nMinimum value of variance: ");
        printf1(PMFmtS,vmin);
        printf1("\nCorresponding mean value:  ");
        printf1(PMFmtS,umin);
        newline();
    }

  
    /** printf("mmin=%f mmax=%f\n",mmin,mmax); **/
  

    mmin1 = mmin;
    mmax1 = mmax;
    fin = 0;

    for (i = 1; i < 10; ++i) {
        n = 0;
        mmin = mmax = 0.0;

        for (j = 0; j < NOC; ++j) {
            if (AcI[j] == 0) {
                n++;

                xl = AcX[j];
                xh = AcY[j];
                tmp = (xl + xh) / 2.0;
                if (tmp <= mmin1) {
                    AcI[j] = 1;
                    mmin += xl;
                    mmax += xl;
                }
                else if (tmp >= mmax1) {
                    AcX[j] = xh;
                    AcI[j] = 1;
                    mmin += xh;
                    mmax += xh;
                }
                else {
                    mmin += xl;
                    mmax += xh;
                }           
            }
            else {
                mmin += AcX[j];
                mmax += AcX[j];
            }
        }
        mmin /= (double)NOC;
        mmax /= (double)NOC;
  
        /** printf("mmin=%f mmax=%f n=%d \n",mmin,mmax,n ); **/
  
        if (n == 0) {
            fin = 1;
            break;
        }
        mmin1 = mmin;
        mmax1 = mmax;
    }
    if (fin == 0) {
        printf1("\nCannot calculate maximum of variance.\n");
    }
    else {
        u = 0.0;
        for (j = 0; j < NOC; ++j)
            u += AcX[j];
        u /= (double)NOC;

        v = 0.0;
        for (j = 0; j < NOC; ++j)  
            v += (AcX[j] - u) * (AcX[j] - u);
        v /= (double)NOC;

        printf1("\nMaximum value of variance: ");
        printf1(PMFmtS,v);
        printf1("\nCorresponding mean value:  ");
        printf1(PMFmtS,u);
        newline();
    }
    newline();
    err = 0;

IVFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ivarf()                                                                 */
/*                                                                          */

double ivarf(int n,double x,double *xl,double *xh)
{
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

