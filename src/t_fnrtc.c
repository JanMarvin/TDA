/****************************************************************************/
/*  t_fnrtc                                                                 */
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
#include "t_gf.h"
#include "t_ml.h"
#include "t_gdat.h"
#include "t_edat.h"
#include "t_rate.h"
#include "t_parm.h"
#include "t_mat.h"
#include "t_matf.h"

/*  functions in t_fnrtc.c */

int fn_dis(int mod);

/****************************************************************************/
/*  fn_dis  Calculation of log-likelihood, gradient, hessian for            */
/*          discrete time models.                                           */
/*                                                                          */
/*      DLR: Logistic Regression                                            */
/*      CLL: Complementary Log-Log                                          */
/*                                                                          */

int fn_dis(int mod)
{
    register int j,j1,k,jv,jk,kv;
    int err,icase,sn,ip,ip1,ip2,spl,nspl,tsi,tfi,tl,next,nexta,org,des,cen;
    double ts,tf,xa,xb,xj,xk,wt,a,b,a1,b1,rj,rj1,rjk,rjk1,dl,dl1;
    double tmp,gtmp,htmp,htmp1,htmpx;

    wt = 1.0;   NOCUsed = 0;

    get_spell(1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
    while (get_spell(0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

        /* ip points to begin of transitions for this spell */

        ip = TranPtr[(sn - 1) * MaxOrg1 + org];
        if (ip < 0)
            continue;

        if (WIVar >= 0)                             /* get weights */
            wt = get_data(WIVar,icase) * WNorm;

        tsi = (int)ts;
        tfi = (int)tf;

        /*  the  main loop goes over all time points for this episode       */

        for (tl = tsi + 1; tl <= tfi; ++tl) {

            /*  calculate the overall rate rj for origin state org and      */
            /*  time point tl.                                              */

            rj = 0.0;

            /*  make a loop over all transitions (ip1) starting with org.    */
            /*  j is start index in the overall parameter vector.           */

            j = PIdxPtr[ip]; 
            ip1 = ip - 1;

            while (sn == SnTran1[++ip1] && org == OrgTran1[ip1]) {

                if (des == DesTran1[ip1])
                    cen = 0;
                else
                    cen = 1;

                nexta = j;      /* nexta is the beginning index for this    */
                next = j - 1;   /* transition.                              */

                /* get parameters and covariate values for this transition  */

                xa = TPar[j];                                                   
                while ((jv = PIdx[++j]))                                                 
                    xa += TPar[j] * get_data(jv - 1,icase);
                a = rexp(xa);

                if (PMDEG) {
                    dl = (double)tl;
                    xb = TPar[j] * dl;                                              
                    while ((jv = PIdx[++j])) {                                               
                        dl *= (double)tl;
                        xb += TPar[j] * dl;
                    }
                    b = rexp(xb);
                }
                else
                    b = 1.0;

                /* calculate the rate, rjk, for this transition and for     */
                /* the point of time l.                                     */

                if (mod == DLR)
                    rjk = a * b / (1.0 + a * b);
                else
                    rjk = 1.0 - rexp(-a * b);

                /* if tl = tf, and if the episode is not censored and has   */
                /* a destination state according to this transition, then   */
                /* calculate contribution to log-likelihood and derivatives */

                if (tl == tfi && cen == 0) {

                    if (LFunc)  
                        FTmp += wt * rlog(rjk);                        

                    /* if necessary calculate contributions to first and    */
                    /* second derivatives.                                  */
        
                    if (LGrad || LSec) {
                                                
                        /* calculate contribution to gradient in gtmp, and  */
                        /* contribution to hessian in htmp, weighted.       */

                        j1 = nexta;   
                        jk = (j1 * (j1 - 1)) / 2 - next;
                        
                        if (mod == DLR) {
                            gtmp = (1.0 - rjk) * wt;
                            htmp = -rjk * gtmp;
                        }
                        else {
                            gtmp = a * b * (1.0 / rjk - 1.0) * wt;
                            htmp = gtmp * (1.0 - a * b / rjk);
                        }
                        dl = 1.0;

                        /* derivatives with respect to alpha, first */
                
                        xj = 1.0;

                        while (1) {
                            if (LGrad) {
                                Grad[j1] += gtmp * xj;
                                if (CGradFlg && NOCUsed < MPGradRow)
                                    MatVal[MPGradIdx][NOCUsed * MPGradCol + j1] += gtmp * xj;
                            }
                            if (LSec) {
                                k = nexta;
                                dl = 1.0;
                                xk = 1.0;
                                jk += next;
                                
                                while (1) {     /* alpha, alpha */
                                    if (k < j1)  
                                        Hess[++jk] += htmp * xj * xk;
                                    else if (k == j1) {
                                        Diag[k] += htmp * xj * xk;
                                        break;
                                    }
                                    if (!(kv = PIdx[++k]))                                                 
                                        break;
                                    xk = get_data(kv - 1,icase);
                                }
                            }
                            if (!(jv = PIdx[++j1]))                                                 
                                break;
                            xj = get_data(jv - 1,icase);
                        }
                        
                        /* derivatives with respect to beta */

                        dl = tf;
                        
                        while (PMDEG) {

                            if (LGrad) {
                                Grad[j1] += gtmp * dl;
                                if (CGradFlg && NOCUsed < MPGradRow)
                                    MatVal[MPGradIdx][NOCUsed * MPGradCol + j1] += gtmp * dl;
                            }
                            if (LSec) {
                                k = nexta;
                                jk += next;
                                xk = 1.0;

                                while (1) {     /* beta, alpha */
                                    if (k < j1)  
                                        Hess[++jk] += htmp * dl * xk;
                                    else if (k == j1) {
                                        Diag[k] += htmp * dl * xk;
                                        break;
                                    }
                                    if (!(kv = PIdx[++k]))                                                 
                                        break;
                                    xk = get_data(kv - 1,icase);
                                }
                                
                                dl1 = dl * tf;
                            
                                while (1) {     /* beta, beta */
                                    if (k < j1)  
                                        Hess[++jk] += htmp * dl1;
                                    else if (k == j1) {
                                        Diag[k] += htmp * dl1;
                                        break;
                                    }
                                    if (!(kv = PIdx[++k]))                                                 
                                        break;
                                    dl1 *= tf;
                                }
                            }   
                            if (!(jv = PIdx[++j1]))                                                 
                                break;
                            dl *= tf;
                        }
                    }
                }

                /* add up the overall transition rate. And save rjk for */
                /* later use with derivatives. If CLL we need in        */
                /* addition transition-specific a and b values.         */

                rj += rjk;
                DRJKW[ip1] = rjk;
                if (mod == CLL) {
                    DAJKW[ip1] = a;
                    DBJKW[ip1] = b;
                }

            }   /* end of loop over all transitions */

            /*  calculate contributions to log-likelihood and derivatives.  */
            /*  omit episodes which are not censored for time point tf.     */
            
            if (tl != tfi || des == org) {
            
                if (LFunc)  
                    FTmp += wt * rlog(1.0 - rj);                   

                if (LGrad || LSec) {
                                                
                    rj1 = 1.0 - rj;

                    /*  we need derivatives for all parameters with         */
                    /*  origin state org, so we need a new loop over all    */
                    /*  transitions starting with org.                      */

                    j = PIdxPtr[ip]; 
                    nexta = j;          /* nexta is the beginning index     */
                    next = j - 1;       /* for this transition. jk is       */
                                        /* the index for the hessian.       */
                    jk = (j * (j - 1)) / 2 - next;

                    ip1 = ip - 1;

                    while (sn == SnTran1[++ip1] && org == OrgTran1[ip1]) {

                        /* calculate contribution to gradient in gtmp and   */
                        /* something to simplify calculation of hessian.    */

                        rjk = DRJKW[ip1];
                        if (mod == DLR) {
                            tmp = -rjk * (1.0 - rjk) / rj1;
                            htmpx = tmp * (1.0 - 2.0 * rjk - tmp);
                            htmp1 = -rjk * (1.0 - rjk) / (rj1 * rj1);
                        }
                        else {
                            tmp = a * b * (rjk - 1.0) / rj1;
                            htmpx = tmp * (1.0 - a * b * (rjk - rj) / rj1);
                            htmp1 = -a * b * (1.0 - rjk) / (rj1 * rj1);
                        }
                        gtmp = wt * tmp;
                        htmpx *= wt;
                        htmp1 *= wt;

                        /* derivatives with respect to alpha, first */
                
                        xj = 1.0;

                        while (1) {

                            if (LGrad) {
                                Grad[j] += gtmp * xj;
                                if (CGradFlg && NOCUsed < MPGradRow)
                                    MatVal[MPGradIdx][NOCUsed * MPGradCol + j] += gtmp * xj;
                            }
                            if (LSec) {

                                /* if second derivatives we need a second   */
                                /* loop over transitions.                   */

                                k = nexta;
                                jk += next;

                                ip2 = ip - 1;
                                while (sn == SnTran1[++ip2] && org == OrgTran1[ip2]) {

                                    if (ip2 > ip1)
                                        break;

                                    if (ip2 < ip1) {
                                        rjk1 = DRJKW[ip2];
                                        if (mod == DLR)
                                            htmp = htmp1 * rjk1 * (1.0 - rjk1);
                                        else {
                                            a1 = DAJKW[ip2];
                                            b1 = DBJKW[ip2];
                                            htmp = htmp1 * a1 * b1 * (1.0 - rjk1);
                                        }
                                    }
                                    else  
                                        htmp = htmpx;

                                    xk = 1.0;

                                    while (k <= j) {    /* alpha, alpha */
                                        if (k < j)  
                                            Hess[++jk] += htmp * xj * xk;
                                        else if (k == j) {
                                            Diag[k++] += htmp * xj * xk;
                                            break;
                                        }
                                        if (!(kv = PIdx[++k]))                                                 
                                            break;
                                        xk = get_data(kv - 1,icase);
                                    }
                                    dl1 = (double)tl;

                                    while (PMDEG && k <= j) {  /* alpha, beta */
                                        if (k < j)   
                                            Hess[++jk] += htmp * xj * dl1;
                                        else if (k == j) {
                                            Diag[k++] += htmp * xj * dl1;
                                            break;
                                        }
                                        if (!(kv = PIdx[++k]))                                                 
                                            break;
                                        dl1 *= (double)tl;
                                    }
                                }
                            }
                            if (!(jv = PIdx[++j]))                                                 
                                break;
                            xj = get_data(jv - 1,icase);
                        }
                        
                        /* derivatives with respect to beta */

                        dl = (double)tl;
                        
                        while (PMDEG) {

                            if (LGrad) {
                                Grad[j] += gtmp * dl;
                                if (CGradFlg && NOCUsed < MPGradRow)
                                    MatVal[MPGradIdx][NOCUsed * MPGradCol + j] += gtmp * dl;
                            }
                            if (LSec) {
                            
                                /* if second derivatives we need a second   */
                                /* loop over transitions.                   */

                                k = nexta;
                                jk += next;

                                ip2 = ip - 1;
                                while (sn == SnTran1[++ip2] && org == OrgTran1[ip2]) {

                                    if (ip2 > ip1)
                                        break;

                                    if (ip2 < ip1) {
                                        rjk1 = DRJKW[ip2];
                                        if (mod == DLR)
                                            htmp = htmp1 * rjk1 * (1.0 - rjk1);
                                        else {
                                            a1 = DAJKW[ip2];
                                            b1 = DBJKW[ip2];
                                            htmp = htmp1 * a1 * b1 * (1.0 - rjk1);
                                        }
                                    }
                                    else  
                                        htmp = htmpx;

                                    xk = 1.0;

                                    while (k <= j) {    /* beta, alpha */
                                        if (k < j)   
                                            Hess[++jk] += htmp * dl * xk;
                                        else if (k == j) {
                                            Diag[k++] += htmp * dl * xk;
                                            break;
                                        }
                                        if (!(kv = PIdx[++k]))                                                 
                                            break;
                                        xk = get_data(kv - 1,icase);
                                    }
                                
                                    dl1 = dl * (double)tl;
                            
                                    while (k <= j) {    /* beta, beta */
                                        if (k < j)   
                                            Hess[++jk] += htmp * dl1;
                                        else if (k == j) {
                                            Diag[k++] += htmp * dl1;
                                            break;
                                        }
                                        if (!(kv = PIdx[++k]))                                                 
                                            break;
                                        dl1 *= (double)tl;
                                    }
                                }
                            }   
                            if (!(jv = PIdx[++j]))                                                 
                                break;
                            dl *= (double)tl;
                        }
                    }
                }
            }

        }   /* end of loop over all time points */

        if (CGradFlg && NOCUsed < MPGradRow && PM2NV > 0)
            mp_putvar(NOCUsed,MPGradCol,MatVal[MPGradIdx],PM2NV,PM2VIdx,icase);
        NOCUsed++;
    }
    err = checkov();
    return(err);
}



