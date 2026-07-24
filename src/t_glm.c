/****************************************************************************/
/*  t_glm                                                                   */
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
#include "t_gdat.h"   
#include "t_eval.h"   
#include "t_eval3.h"
#include "t_alloc.h"   
#include "t_ml.h"   
#include "t_gmin.h"   
#include "t_lsei.h"   
#include "t_gf.h"   
#include "t_var.h"   
#include "t_con.h"   
#include "t_cdf.h"   

/*  functions in t_glm.c */

int glm(void);
int prn_glmd(int d);
void glm_prot(int nx);
int glm_it(int nx,int intflg,int ne,int ni);
int get_mue(double eta,double *mue);
int rtbis(double x1,double x2,double *x,double eta,double *dx);
int glm_lfd(double mue,double *d);
int glm_dev(int nx,int intflg,int df);
void glm_dtda(int intflg,int nx);

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */

int GLMD = 1;           /* type of distribution                             */
int GLMLF = 0;          /* type of link function                            */
int GLMIT = 0;          /* number of iterations                             */
int GLMRK = 0;          /* rank                                             */  
double GLMPC = 0.0;     /* scaled parameter change                          */
double GLMDEV = 0.0;    /* deviance                                         */
double GLMX = 0.0;      /* Pearson statistic                                */
double PHI_ML = 0.0;    /* ML-based scale parameter                         */
double PHI_DB = 0.0;    /* deviance-based scale parameter                   */

/* ------------------------------------------------------------------------ */
/*  glm()   Generalized linear models.                                      */  
/*                                                                          */
/*          glm(                                                            */
/*              d =...,     1 = normal                                      */
/*                          2 = binomial                                    */
/*                          3 = poisson                                     */
/*                          4 = gamma                                       */
/*                          5 = inverse gaussian                            */
/*              link=...,   1 = identity                                    */
/*                          2 = log                                         */
/*                          3 = logit                                       */
/*                          4 = reciprocal                                  */
/*                          5 = probit                                      */
/*                          6 = comp. log-log                               */
/*                          7 = square root                                 */
/*                          8 = quadratic inverse                           */
/*              v=...,      varlist: Y,X1,X2,...                            */
/*              yw=...,     variable for proportions                        */
/*              ni=...,     1 if without intercept                          */
/*              lsecon=     equality constraints                            */
/*              lsicon=     inequality constraints                          */
/*              mxit=...,   max number of iterations, def. 20               */
/*              tolsp=...,  tolerance for convergence, def. 1e-6            */
/*              xp=...,     starting values                                 */
/*              dsv=...,    starting values                                 */
/*              tfmt=...,   print format for parameters, def. 10.4          */
/*              ppar=...,   print parameter                                 */
/*              pcov=...,   print covariance matrix                         */
/*              mfmt=...,   print format pcov                               */
/*              pres=...,   print residuals                                 */
/*              fmt=...,    print format pres                               */
/*              prot=...,   protocol file                                   */
/*              pfmt=...,   print format for protocol file, def. -19.11     */
/*              ab=...,     domain for link function, def. (0,1)            */
/*              tolf=...,   tolerance for inverse link function, 10-12      */
/*                                                                          */
/*          ) = link_function;                                              */
/*                                                                          */
/*          Return 0 if OK, otherwise -1.                                   */

int glm(void)
{
    register int i;
    int err,r,nv,nx,nx1,nw,mw,intflg,df,ncon,ne,ni;

    err = -1;
    if (check_cmd(2))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    printf1("Generalized linear models. Current memory: %d bytes.\n",MemReq);

    TOLSP = 1.e-6;
    TOLF  = 1.e-12;

    if (parm(CmdBuf + 3,9,0))   /* get parameters */
        goto GLMFin;

    if (PMNW != 1) {
        printf1("Error: glm command requires cross-section data (nw=1).\n");
        goto GLMFin;
    }
    if (PMFmtF == 0)                /* default print format */
        pmfmt(10,4);

    newline();
    if (prn_glmd((int)PMD))
        goto GLMFin;

    nv = check_nvar(0);         /* check variables */
    if (nv == 0) {              /* PMNV is number of variables */       
        printf1("Error: no variables.\n");
        goto GLMFin;        
    }
    prn_nwvar(0);               /* print variables */

    nx = nx1 = nv - 1;
    if (PMNI) {
        printf1("Model without intercept.\n");
        PMNI = 1;
        intflg = 0;
    }
    else {
        PMNI = 0;
        intflg = 1;
        nx1++;
        if (nx == 0)  
            printf1("Model without independent variables.\n");
    }
    if (nx1 == 0) {
        printf1("Error: no model parameters.\n");
        goto GLMFin;
    }
    if (PMYWVar >= 0) {
        printf1("Variable used to define proportions: %s",VName[PMYWVar]);
        if (PMD != 2) {
            printf1(" (will be ignored)");
            PMYWVar = -1;
        }
        printf1("\n\n");
    }   

    /* allocate memory for derivatives of link function */

    if (fnd_alloc(1,1,FNArgN,FNPN,0))
        goto GLMFin;

    /* data matrix: AcW, parameter vector: AcX */

    ne = ni = ncon = 0;         /* number of constraints */
    nw = nx1 + 1;
    r = mw = NOC + NCONSTR;     
    if (r < nw)
        r = nw;

    if (alloc_acw(r * nw + 1))
        goto GLMFin;

    if (alloc_acx(nx1 + 1))
        goto GLMFin;

    if (alloc_acy(nx1 + 1))
        goto GLMFin;

    if (MxItFlg == 0)
        MxIter = 20;

    printf1("Estimation with iteratively re-weighted least squares.\n");
    if (WIVar >= 0)  
        printf1("Using weights defined by: %s\n",VName[WIVar]);

    printf1("Number of model parameters: %d\n",nx1);
         
    if (NCONSTR > 0) {
        newline();
        if (p_con(CmdBuf,nx,PMNI,mw,nw,AcW,&ne,&ni))      /* get constraints */
            goto GLMFin;
        ncon = ne + ni;
        printf1("Equality constraints: %d\n",ne);        
        printf1("Inequality constraints: %d\n\n",ni);          

        /* need some additional memory */

        if (alloc_acu(ne * nw + 1))
            goto GLMFin;

        if (alloc_acv(ne * nw + 1))
            goto GLMFin;

        r = ne * nw;
        for (i = 1; i <= r; ++i)
            AcU[i] = AcW[i];

        r = ni * nw;
        df = (ne + NOC) * nw;
        for (i = 1; i <= r; ++i)
            AcV[i] = AcW[df + i];
    }
    if (NOC + ne < nx1 + 1) {
        printf1("Error: need at least %d cases.\n",nx1 + 1);
        goto GLMFin;
    }
    printf1("Maximum number of iterations: %d\n",MxIter);
    printf1("Tolerance for scaled parameter change: %g\n",TOLSP);

    if (get_dsv(nx1,AcX,0,ParLB,ParUB,1))   /* try to get starting values */
        goto GLMFin;                        /* ParLB,ParUB not used */

    glm_prot(nx1);      /* init protocol file */
    newline();

    r = glm_it(nx,intflg,ne,ni);    /* perform iterations */
    if (r < 0) {
        if (r == -2)  
            printf1("Convergence not reached.\nExceeded maximum number of iterations.\n");
        goto GLMFin;
    }
    printf1("Convergence reached in %d iterations.\n",GLMIT);
    printf1("Final scaled parameter change: %g\n",GLMPC);

    df = NOC - GLMRK;
    if (df < 0)
        df = 0;

    if (glm_dev(nx,intflg,df))      /* calculate and print deviance etc */
        goto GLMFin;

    if (r > 0)      /* no covariance matrix calculated */
        df = 0;
    else {
        for (i = 1; i <= nx1; ++i)  
            AcY[i] = AcW[(i - 1) * nw + i];
    }  
    prn1_coeff(nx1,AcX,AcY,PMNI,df,PMVIdx,1);

    newline();
    if (r == 0 && PMCovFDef) {                          /* write cov matrix */
        prn_data(nx1,nw,nx1,AcW,PMCovFd,PMMFmtS);
        p_wmsg(2,PMCovFName,1);
    }
    err = 0;

GLMFin:
    fnd_alloc(0,0,0,0,0);
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_glmd(d)     Type type of distribution and set GLMD and GLMLF        */

int prn_glmd(int d)
{
    int nl;

    printf1("Distribution: ");
    switch (d) {
        case 2: printf1("binomial.\n");
                nl = 3;
                break;
        case 3: printf1("Poisson.\n");
                nl = 2;
                break;
        case 4: printf1("gamma.\n");
                nl = 4;
                break;
        case 5: printf1("inverse Gaussian.\n");
                nl = 8;
                break;
        default: printf1("normal.\n");
                d = 1;
                nl = 1;
                break;
    }
    GLMD = d;

    printf1("Link function: ");

    if (FNFlg) {        /* user-defined link function */
        GLMLF = 0;
        printf1("user-defined.\n");
        if (FNArgN != 1) {
            printf1("Error: link function must have exactly one argument.\n");
            return(-1);  
        }
        if (FVFlg) {
            printf1("Error: link function must not refer to data matrix variables.\n");
            return(-1);   
        }
        if (PMXYFlg == 0) {
            PMX = 0.0;
            PMY = 1.0;
        }
        if (PMX >= PMY) {
            printf1("Error in definition of range for mue.\n");
            return(-1);  
        }
        printf1("Domain of link function: (%g,%g)\n",PMX,PMY);
    }
    else {
        if (PMLINK >= 1 && PMLINK <= 8)
            GLMLF = PMLINK;
        else
            GLMLF = nl;

        switch (GLMLF) {
            case  1:    printf1("identity.\n"); break;
            case  2:    printf1("log.\n"); break;
            case  3:    printf1("logit.\n"); break;
            case  4:    printf1("reciprocal.\n"); break;
            case  5:    printf1("probit.\n"); break;
            case  6:    printf1("complementary log-log.\n"); break;
            case  7:    printf1("square root.\n"); break;
            case  8:    printf1("quadratic inverse.\n"); break;
        }
    }
    newline();
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  glm_prot()  protocol file.                                              */

void glm_prot(int nx)
{
    register int i;

    if (PMProtFDef == 0)
        return;

    printf1("Protocol will be written to: %s\n",PMProtFName);

    fprintf(PMProtFd,"Generalized linear model.\n");
    fprintf(PMProtFd,"Number of model parameters: %d\n",nx);
    fprintf(PMProtFd,"\nStarting values.\n");
    for (i = 1; i <= nx; ++i)
        fprintf(PMProtFd,PMPFmtS,AcX[i]);
    fprintf(PMProtFd,"\n\n");
}

/* ------------------------------------------------------------------------ */
/*  glm_it(nx,intflg,ne,ni) Iterated weighted regression for glm.           */
/*                          nx = number of x variables, intflg = 1 if       */
/*                          with intercept. ne and ni are the numbers of    */
/*                          equality and inequality constraints, resp.      */
/*                                                                          */
/*                          AcX contains current parameter vector.          */
/*                          Function requires AcY and AcW as working        */
/*                          space.                                          */
/*                                                                          */
/*  If PMResFDef write data and estimated Y values into output file, and    */  
/*  optionally an TDA description file.                                     */
/*                                                                          */
/*  Return:  0  if successful.                                              */
/*          -1  if error                                                    */
/*          -2  if max number of iterations reached.                        */
/*           1  if covariance matrix not calculated.                        */

int glm_it(int nx,int intflg,int ne,int ni)
{
    register int i,j,k,iw;
    int iter,err,nw,nx1,ranka,ranke,cov;
    double y,yn,eta,tmp,mue,d,w,rnorme,rnorml,wt,wsum;

    rnorme = rnorml = 0.0;
    ranka = ranke = GLMIT = GLMRK = 0;
    nx1 = nx;
    if (intflg)
        nx1++;
    nw = nx1 + 1;
    cov = 0;
    yn = wt = 1.0;
               
    if (SILENTFlg < 2)
        printfe("\n  Iter    Function Value\n");

    for (iter = 0; iter <= MxIter; ++iter) {

        for (j = 1; j <= nx1; ++j)      /* save current parameters */
            AcY[j] = AcX[j];

        if (PMProtFDef)
            fprintf(PMProtFd,"\nIteration %d\n",iter);
 
        /* set up matrix for weighted regression */

        iw = ne;
        wsum = 0.0;

        for (i = 0; i < NOC; ++i) {

            y = get_data(PMVIdx[0],i);

            eta = 0.0;                  /* eta = xi * beta */
            j = 1;
            if (intflg)
                eta += AcX[j++];
            for (k = 1; k <= nx; ++k) {
                tmp = get_data(PMVIdx[k],i);
                eta += tmp * AcX[j++];
            }
            if (get_mue(eta,&mue))      /* g(mue) = eta */
                return(-1);
   
            if (PMYWVar >= 0 && PMD == 2) {
                yn = get_data(PMYWVar,i);
                if (y < 0.0 || yn < 1.0 || y > yn) {
                    printf1("Error in case %d: y = %g, count = %g\n",i+1,y,yn);
                    return(-1);
                }
                mue *= yn;
            }

            /* ### evaluate link function with mue, and get first derivative in d */
      
            /* printf("iter=%d i=%d eta=%g mue=%g yn=%g\n",iter,i,eta,mue,yn); */
   
            tmp = mue / yn;
            if ((err = glm_lfd(tmp,&d))) {
                printf1("Error: cannot evaluate link function, or derivative, with mue = %g (case %d)\n",mue / yn,i + 1);
                if (err > 0)
                    prn_emsg2(err);
                return(-1);
            }

            /* printf("yn=%g d=%g\n",yn,d); */

            if (PMYWVar >= 0 && PMD == 2)
                d /= yn;
     
            /* calculate working weights */

            if (WIVar >= 0) {
                wt = get_data(WIVar,i);
                if (wt < 0.0) {
                    printf1("Error: found negative weight in case %d.\n",i + 1);
                    return(-1);
                }
            }
            wsum += wt;

            tmp = d * d;     
            switch (GLMD) {
                case 2:     tmp *= mue * (1.0 - mue / yn);
                            break;
                case 3:     tmp *= mue;                 
                            break;
                case 4:     tmp *= mue * mue;
                            break;
                case 5:     tmp *= mue * mue * mue;
                            break;
                default:    break;
            }
            if (wt > 0.0) {
                if (tmp <= 0.0) {
                    printf1("Error: cannot calculate weight in case %d\n",i + 1);
                    return(-1);
                }
                w = sqrt(wt / tmp);
            }
            else     
                w = 0.0;
                 
            j = 1;
            if (intflg) {
                AcW[iw * nw + j] = w;
                j++;
            }
            for (k = 1; k <= nx; ++k) {
                tmp = get_data(PMVIdx[k],i);
                AcW[iw * nw + j] = tmp * w;
                j++;
            }
            if (w > 0.0)   
                AcW[iw * nw + j] = w * (eta  + (y - mue) * d);
            else
                AcW[iw * nw + j] = 0.0;

            if (cov && PMResFDef) {

                fprintf(PMResFd,"%6d ",i + 1);
                tmp = get_data(PMVIdx[0],i);
                fprintf(PMResFd,PMFmtS,y);
                fprintf(PMResFd,PMFmtS,mue);
                fprintf(PMResFd,PMFmtS,eta);
                for (k = 1; k <= nx1; ++k) {
                    tmp = AcW[iw * nw + k];           
                    if (w > 0.0)
                        tmp /= w;
                    fprintf(PMResFd,PMFmtS,tmp);           
                }
                fprintf(PMResFd,PMFmtS,wt);
                fprintf(PMResFd,PMFmtS,w * w);
                fprintf(PMResFd,"\n");
            }
            iw++;
        }
        if (cov && PMResFDef) {
            printf1("Data and estimated values written to: %s\n",PMResFName);
            if (PMTDAFDef)  
                glm_dtda(intflg,nx);
            newline();
            if (cov == 2)
                return(1);
        }
        if (wsum < EPSI1) {
            printf1("Error: sum of weights is almost zero.\n");
            return(-1);
        }
        if (ne > 0) {       /* add equality constraints */
            j = ne * nw;
            for (i = 1; i <= j; ++i)
                AcW[i] = AcU[i];
        }
        if (ni > 0) {       /* add inequality constraints */
            j = ni * nw;
            k = (ne + NOC) * nw;
            for (i = 1; i <= j; ++i)
                AcW[k + i] = AcV[i];
        }

        /**** 
        for (i = 0; i < NOC + ne + ni; ++i) {
            for (j = 1; j <= nw; ++j) 
                printf1("%f ",AcW[i * nw + j]);
            printf1("\n");
        }
        newline();
        **********/
         
        /*  solve least squares problem */

        err = lsei(AcW,ne,NOC,ni,nx1,AcX,cov,&rnorme,&rnorml,&ranka,&ranke);

        if (err == -1) {
            printf1("Equality constraints are contradictory.\n");
            printf1("Calculated a least squares solution.\n");
        }
        else if (err) {
            printf1("Error: cannot solve least squares problem (%d)\n",err);
            if (err == -2)  
                printf1("Inequality constraints are incompatible.\n");
            else if (err == -3)  
                printf1("Equality and inequality constraints are contradictory.\n");
            else if (err == -4)  
                p_err(-2,1);
            return(-1);
        }
        if (ranka != nx1 - ranke) {
            printf1("Error: rank of least squares data matrix: %d\n",ranka);
            return(-1);
        }
        GLMPC = 0.0;
        for (j = 1; j <= nx1; ++j) {
            tmp = fabs(AcX[j] - AcY[j]) / dmax(fabs(AcX[j]),1.0);                  
            GLMPC = dmax(GLMPC,tmp);                                            
        }
        if (SILENTFlg < 2)
            printfe("  %3d  %20.13e\n",iter,GLMPC);

        if (PMProtFDef) {
            fprintf(PMProtFd,"Rank of least squares data matrix: %d\n",ranka);
            fprintf(PMProtFd,"Norm of least squares residuals: %g\n",rnorml); 
            if (ne > 0) {
                fprintf(PMProtFd,"Rank of equality constraints: %d\n",ranke);
                fprintf(PMProtFd,"Norm of residuals of equality constraints: %g\n",rnorme); 
            }
            prvec("\nParameter vector",nx1,AcX);
            prval("Scaled parameter change",GLMPC);
        }
        if (GLMPC <= TOLSP) {
            GLMIT = iter;
            GLMRK = ranka;

            /* do not calculate cov matrix if ... */
            if (ni > 0 || NOC <= ranka || fabs(rnorme + rnorml) < EPSI1) {
                if (PMResFDef) {
                    if (cov)
                        return(1);
                    else
                        cov = 2;
                }
                else
                    return(1);
            }
            else if (cov == 1) {    /* convergence reached and covariance
                                       matrix successfully calculated */

                /* if not normal distribution, we need to rescale the
                   covariance matrix */

                if (NOC <= GLMRK || rnorml < EPSI)
                    return(1);

                if (GLMD == 1)
                    return(0);

                tmp = (double)(NOC - GLMRK) / (rnorml * rnorml);
                for (i = 1; i <= nx1; ++i) {
                    for (j = 1; j <= nx1; ++j)
                        AcW[(i - 1) * nw + j] *= tmp;
                }
                return(0);
            }
            else if (cov == 2)
                return(1);
            else  
                cov = 1;
        }
    }
    return(-2);
}

/* ------------------------------------------------------------------------ */
/*  get_mue(eta,mue)    get mue such that g(mue) = eta.                     */
/*                      return 0 if successful, -1 if error.                */

int get_mue(double eta,double *mue)
{
    int err;
    double dx;

    dx = 0.0;
    err = 0;
    switch (GLMLF) {  
        case 0:     err = rtbis(PMX,PMY,mue,eta,&dx);   /* user-defined */
                    break;
        case 1:     *mue = eta;                         /* identity */
                    break;
        case 2:     *mue = rexp(eta);                   /* log */
                    break;
        case 3:     dx = rexp(eta);                     /* logit */
                    *mue = dx / (1.0 + dx);
                    break;
        case 4:     if (eta != 0.0)                     /* reciprocal */
                        *mue = 1.0 / eta;
                    else
                        err = -1;
                    break;
        case 5:     *mue = cdnf(eta);                   /* probit */
                    break;
        case 6:     dx = rexp(eta);                     /* compl. log-log */
                    *mue = 1.0 - rexp(-dx);
                    break;
        case 7:     *mue = eta * eta;                   /* square root */
                    break;
        case 8:     if (eta <= 0.0)                     /* quadratic inverse */
                        err = -1;
                    *mue = sqrt(1.0 / eta);
                    break;
    }
    if (err) {
        printf1("Error: cannot find mue corresponding to eta = %g [dx=%g]\n",eta,dx);
        if (GLMLF == 0)
            printf1("Check link function and range for mue.\n");
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  rtbis(x1,x2,x,eta)      find x such that g(x) = eta.                    */
/*  ###                     g(x) is the link function                       */
/*                          x1,x2 should be brackets such that              */
/*                          g(x1) <= eta <= g(x2) or                        */
/*                          g(x2) <= eta <= g(x1)                           */
/*                                                                          */
/*  Return 0 if successful, otherwise -1.                                   */

#define RTJMAX 500

int rtbis(double x1,double x2,double *x,double eta,double *dx)
{
    int j,err;
    double f1,f2,fmid,xmid;
    double tol,eps;

    eps = 10.e-15;
    tol = TOLF;     
    *dx = 0.0;
    FNArgVal[0] = x1 + eps;
    err = get_fival(0,&f1,0,0);
    /* printf1("x1=%g f1=%g err=%d\n",x1,f1,err); */ 
    if (err)
        return(err);
    f1 -= eta;

    FNArgVal[0] = x2 - eps;
    err = get_fival(0,&f2,0,0);
    /* printf1("x2=%g f2=%g err=%d\n",x2,f2,err); */
    if (err)
        return(err);
    f2 -= eta;

    if (f1 * f2 >= 0.0)  
        return(-1);
 
    if (f1 < 0.0) {
        *x = x1;
        *dx = x2 - x1;
    }
    else {
        *x = x2;
        *dx = x1 - x2;
    }
    for (j = 0; j < RTJMAX; ++j) {
        *dx *= 0.5;
        xmid = *x + *dx;
        FNArgVal[0] = xmid;
        err = get_fival(0,&fmid,0,0);
        if (err)
            return(err);
        fmid -= eta;
        if (fmid <= 0.0)
            *x = xmid;
        if (fabs(*dx) < tol || fmid == 0.0) {
            /* printf1("eta=%20.16e x=%20.16e\n",eta,*x); */ 
            return(0);
        }
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  glm_lfd(mue,d)  Evaluate link function with mue, and get first          */
/*                  derivative, return in d.                                */
/*                  If OK return 0, else -1.                                */

int glm_lfd(double mue,double *d)
{
    int err;
    double tmp;

    *d = 0.0;
    err = 0;
    switch (GLMLF) {  
        case 0:     FNArgVal[0] = mue;                  /* user-defined */
                    err = get_fival(0,&tmp,1,0);
                    if (err == 0)
                        *d = FNGrad[0][0];
                    break;
        case 1:     *d = 1.0;                           /* identity */
                    break;
        case 2:     if (mue == 0.0)                     /* log */
                        err = -1;
                    else
                        *d = 1.0 / mue;
                    break;
        case 3:     tmp = mue * (1.0 - mue);            /* logit */
                    if (tmp == 0.0)          
                        err = -1;
                    else
                        *d = 1.0 / tmp;
                    break;
        case 4:     if (mue == 0.0)                     /* reciprocal */
                        err = -1;         
                    else
                        *d = -1.0 / (mue * mue);
                    break;
        case 5:     if (mue <= 0.0 || mue >= 1.0)
                        err = -1;
                    else {
                        tmp = cdnif1(mue);
                        *d = 1.0 / dnf(tmp);            /* probit */
                    }
                    break;
        case 6:     tmp = 1.0 - mue;
                    if (tmp <= 0.0)                     /* compl. log-log */
                        err = -1;            
                    else  
                        *d = -1.0 / (tmp * rlog(tmp));
                    break;
        case 7:     if (mue <= 0.0)                     /* square root */
                        err = -1;
                    else
                        *d = 1.0 / (2.0 * sqrt(mue));
                    break;
        case 8:     tmp = mue * mue * mue;              /* quadratic inverse */
                    if (tmp == 0.0)
                        err = -1;
                    else
                        *d = -2.0 / tmp;
                    break;
    }
    return(err);
}
  
/* ------------------------------------------------------------------------ */
/*  glm_dev(nx,intflg,df)   This function calculates:                       */
/*                                                                          */
/*  GLMDEV = deviance                                                       */
/*  GLMX   = Pearson statistic                                              */
/*  PHI_ML = ML-based scale parameter                                       */
/*  PHI_DB = deviance-based scale parameter                                 */
/*                                                                          */  
/*  Return:  0  if successful, -1 if error.                                 */

int glm_dev(int nx,int intflg,int df)
{
    register int i,j,k;
    int nw,nx1;
    double eta,y,yn,mue,wt,wsum,tmp,tmp1,dev;

    yn = 1.0;
    nx1 = nx;
    if (intflg)
        nx1++;
    nw = nx1 + 1;
    wt = 1.0;
    wsum = 0.0;

    GLMDEV = 0.0;   /* deviance */
    GLMX = 0.0;     /* Pearson statistic */
    PHI_ML = 0.0;   /* ML-based scale parameter */
    PHI_DB = 0.0;   /* deviance-based scale parameter */

    for (i = 0; i < NOC; ++i) {

        y = get_data(PMVIdx[0],i);
        eta = 0.0;
        j = 1;
        if (intflg)
            eta += AcX[j++];
        for (k = 1; k <= nx; ++k) {
            tmp = get_data(PMVIdx[k],i);
            eta += tmp * AcX[j++];
        }
        if (get_mue(eta,&mue))  
            return(-1);

        if (PMYWVar >= 0 && PMD == 2) {
            yn = get_data(PMYWVar,i);
            if (y < 0.0 || yn < 1.0 || y > yn) {
                printf1("Error in case %d: y = %g, count = %g\n",i+1,y,yn);
                return(-1);
            }
            mue *= yn;
        }
        if (WIVar >= 0) {
            wt = get_data(WIVar,i);
            if (wt < 0.0) {
                printf1("Error: found negative weight in case %d.\n",i + 1);
                return(-1);
            }
        }
        wsum += wt;

        dev = 0.0;
        switch (GLMD) {
            case 2:     if (y > 0.0)                /* binomial */
                            dev += 2.0 * y * rlog(y / mue);
                        tmp = yn - y;
                        tmp1 = yn - mue;
                        if (tmp > 0.0 && tmp1 > 0.0)
                            dev += 2.0 * tmp * rlog(tmp / tmp1);
                        break;
            case 3:     if (y > 0.0)                /* Poisson */
                            dev = 2.0 * y * rlog(y / mue);
                        dev -= 2.0 * (y - mue);
                        break;

            case 4:     if (y > 0.0 && mue > 0.0)              /* Gamma */
                            dev = -2.0 * rlog(y / mue);
                        if (mue > 0.)
                            dev += 2.0 * (y - mue) / mue;
                        break;

            case 5:     if (y > 0.0 && mue > 0.0)   /* Inverse Gaussian */
                            dev = (y - mue) * (y - mue) / (y * mue * mue);
                        break;

            default:    tmp = y - mue;              /* normal */
                        tmp1 = tmp * tmp;
                        dev = tmp1;
                        GLMX += tmp1;
                        PHI_ML += tmp1;
                        break;
        }
        GLMDEV += wt * dev;
    }
    if (wsum < EPSI1) {
        printf1("Error: sum of weights is almost zero.\n");
        return(-1);
    }
    switch (GLMD) {
        case  2:    break;
        case  3:    break;
        case  4:    break;
        case  5:    break;
        default:    PHI_ML /= wsum;
                    break;
    }
    if (df > 0)
        PHI_DB = GLMDEV / (double)df;

    printf1("\nRank of data matrix: %d\n",GLMRK);
    printf1("Degrees of freedom: %d\n",df);
    prn_sfmt("Deviance",31,PMTFmtS,GLMDEV);

    if (GLMD == 1) {
        prn_sfmt("Pearson statistic",31,PMTFmtS,GLMX);
        prn_sfmt("ML-based scaling factor",31,PMTFmtS,PHI_ML);
        prn_sfmt("Deviance-based scaling factor",31,PMTFmtS,PHI_DB);
    }
    else if (GLMD == 2)
        prn_sfmt("Fixed scaling factor",31,PMTFmtS,1.0);


    return(0);
}

/* ------------------------------------------------------------------------ */
/*  glm_dtda()      write TDA description file                              */

void glm_dtda(int intflg,int nx)
{
    register int i,j,k;
   
    if (PMTDAFDef) {
        fprintf(PMTDAFd,"# data written by glm command.\n");
        fprintf(PMTDAFd,"nvar(\n");
        fprintf(PMTDAFd,"  dfile = %s,\n",PMResFName);
        fprintf(PMTDAFd,"  noc = %d,\n",NOC);
        k = 1;
        fprintf(PMTDAFd,"  Case");
        fprnchar(PMTDAFd,' ',VNameLen - 4,0);
        fprintf(PMTDAFd,"[6.0] = c%-2d,\n",k++);
        j = PMVIdx[0];
        fprintf(PMTDAFd,"  %s",VName[j]);
        fprnchar(PMTDAFd,' ',VNameLen - strlen(VName[j]),0);
        fprintf(PMTDAFd,"[%d.%d] = c%-2d,\n",PMFmt1,PMFmt2,k++);

        fprintf(PMTDAFd,"  Mue");
        fprnchar(PMTDAFd,' ',VNameLen - 3,0);
        fprintf(PMTDAFd,"[%d.%d] = c%-2d,\n",PMFmt1,PMFmt2,k++);

        fprintf(PMTDAFd,"  Eta");
        fprnchar(PMTDAFd,' ',VNameLen - 3,0);
        fprintf(PMTDAFd,"[%d.%d] = c%-2d,\n",PMFmt1,PMFmt2,k++);

        if (intflg) {
            fprintf(PMTDAFd,"  Int");
            fprnchar(PMTDAFd,' ',VNameLen - 3,0);
            fprintf(PMTDAFd,"[%d.%d] = c%-2d,\n",PMFmt1,PMFmt2,k++);
        }
        for (i = 1; i <= nx; ++i) {
            j = PMVIdx[i];
            fprintf(PMTDAFd,"  %s",VName[j]);
            fprnchar(PMTDAFd,' ',VNameLen - strlen(VName[j]),0);
            fprintf(PMTDAFd,"[%d.%d] = c%-2d,\n",PMFmt1,PMFmt2,k++);
        }
        fprintf(PMTDAFd,"  Weight");
        fprnchar(PMTDAFd,' ',VNameLen - 6,0);
        fprintf(PMTDAFd,"[%d.%d] = c%-2d,\n",PMFmt1,PMFmt2,k++);

        fprintf(PMTDAFd,"  WWeight");
        fprnchar(PMTDAFd,' ',VNameLen - 7,0);
        fprintf(PMTDAFd,"[%d.%d] = c%-2d,\n",PMFmt1,PMFmt2,k++);
        fprintf(PMTDAFd,");\n");
    
        printf1("TDA description written to: %s\n",PMTDAFName);
    }
}


