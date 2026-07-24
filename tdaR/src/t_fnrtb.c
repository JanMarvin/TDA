/****************************************************************************/
/*  t_fnrtb                                                                 */
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
#include "t_int.h"
#include "t_cdf.h"
#include "t_mat.h"
#include "t_matf.h"
#include "tda_context.h"

/*  functions in t_fnrtb.c */

int fn_exp(TDAContext *ctx, int resid);
int fn_exp1(TDAContext *ctx, int resid);
int fn_exp2(TDAContext *ctx, int resid);
int fn_pol(TDAContext *ctx, int resid);
int fn_pol1(TDAContext *ctx, int resid);
double get_mpol1(TDAContext *ctx, double t,int *err); 
double get_mpol2(TDAContext *ctx, double t,int *err); 
double get_mpol3(TDAContext *ctx, double t,int *err); 
int fn_gm(TDAContext *ctx, int resid);
int fn_wei(TDAContext *ctx, int resid);
int fn_sic(TDAContext *ctx, int resid);
int fn_ll(TDAContext *ctx, int resid);
int fn_ll2(TDAContext *ctx, int resid);
int fn_ln(TDAContext *ctx, int resid);
int fn_ig(TDAContext *ctx, int resid);
int fn_gam(TDAContext *ctx, int resid); 


/* ------------------------------------------------------------------------ */
/*  fn_exp(resid)                                                           */
/*      Exponential model with optional gamma mixture.                      */
/*      Calculation of log-likelihood, gradient, hessian.                   */
/*      Episode splitting possible if no mixture.                           */
/*      If resid != 0 print residuals.                                      */

int fn_exp(TDAContext *ctx, int resid)
{
    register int j,k,jv,jk,kv,ip;
    int err,icase,sn,next,org,des,cen,fin,fin1,spl,nspl;
    double ts,tf,xa,xd,xj,xk,wt,a,d = 0.0,hc,ldh = 0.0,dh,dh1 = 0.0,dh2 = 0.0,dh3 = 0.0,ui = 0.0;
    double tmp1,tmp2,tmp3 = 0.0,tmp4,tmp5 = 0.0;

    err = 0;    ctx->NOCUsed = 0;
    wt = 1.0;
    get_spell(ctx, 1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
    while (get_spell(ctx, 0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

        /***********
        printf1(ctx, "ic=%d sn=%d  org=%d des=%d ts=%lg tf=%lg\n",icase,sn,org,des,ts,tf);
        *******/
        if (ctx->WIVar >= 0)                             /* get weights */
            wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;

        ip = ctx->TranPtr[(sn - 1) * ctx->MaxOrg1 + org];
        if (ip < 0)  
            continue;
 
        j = ctx->PIdxPtr[ip];
        /*****************
        printf1(ctx, "icase=%d ip=%d sn=%d org=%d des=%d\n",icase,ip,sn,org,des);
        ***********/

        while (sn == ctx->SnTran1[ip] && org == ctx->OrgTran1[ip]) {
            /*****************
            printf1(ctx, "in while: ip=%d sn(ip)=%d org(ip)=%d des(ip)=%d\n",
            ip,SnTran1[ip],OrgTran1[ip],DesTran1[ip]);
            *************/
            if (des == ctx->DesTran1[ip])
                cen = 0;
            else
                cen = 1;
            ip++;

            next = j;                                                       

            xa = ctx->TPar[j];                                                   
            while ((jv = ctx->PIdx[++j]))                                                 
                xa += ctx->TPar[j] * get_data(ctx, jv - 1,icase);

            a = rexp(ctx, xa);
            hc = a * (tf - ts);

            /* print residuals */
   
            if (resid) {
                tmp1 = rexp(ctx, -hc);
                prn_resid1(ctx, icase + 1,org,des,ts,tf,a,tmp1,hc,wt);            
            }
            if (ctx->MixTyp == 1) {
                xd = ctx->TPar[j];                                                   
                while ((jv = ctx->PIdx[++j]))                                                 
                    xd += ctx->TPar[j] * get_data(ctx, jv - 1,icase);

                d = rexp(ctx, xd);
                ui = -1.0 / d;
                if (!cen)
                    ui -= 1.0;
                ldh = rlog(ctx, 1.0 + d * hc);
                dh = d / (1.0 + d * hc);
                dh1 = dh * hc;
                dh2 = 1.0 - dh1;
                dh3 = dh1 * dh2;
            }
        
            if (ctx->LFunc) {

                if (ctx->MixTyp == 1) 
                    tmp1 = ui * ldh;
                else
                    tmp1 = -hc;

                if (!cen)  
                    tmp1 += xa;
                tmp1 *= wt;
                if (ctx->LFunc)  
                    ctx->FTmp += tmp1;                                     
            }
            if (ctx->LGrad || ctx->LSec) {
                j = next;
                jk = j * (j - 1) / 2;
  
                if (ctx->MixTyp == 1) {
                    tmp1 = ui * dh1;
                    tmp2 = tmp1 * dh2;
                    tmp3 = dh1 / d + ui * dh3;
                    tmp4 = ldh / d + ui * dh1;
                    tmp5 = (2.0 * dh1 - ldh) / d + ui * dh3;
                    tmp3 *= wt;
                    tmp4 *= wt;
                    tmp5 *= wt;
                }
                else
                    tmp1 = tmp2 = -hc;

                if (!cen)  
                    tmp1 += 1.0;
                tmp1 *= wt;
                tmp2 *= wt;

                jv = 1;
                xj = 1.0;
                fin = 0;

                while (1) {
                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp1 * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp1 * xj;
                    }
                    if (ctx->LSec)  
                        ctx->Diag[j] += tmp2 * xj * xj;

                    if (!(jv = ctx->PIdx[++j])) {                                               
                        if (ctx->MixTyp == 0 || fin)
                            break;

                        tmp1 = tmp4;
                        tmp2 = tmp5;
                        xj = 1.0;
                        fin = 1;
                    }
                    else  
                        xj = get_data(ctx, jv - 1,icase);

                    if (ctx->LSec) {
                        k = next;                                                       
                        jk += next;                                                   
                        xk = 1.0;
                        fin1 = 0;
                        while (1) {
                            if (!fin || fin1)
                                ctx->Hess[jk] += tmp2 * xj * xk;
                            else
                                ctx->Hess[jk] += tmp3 * xj * xk;
                            k++;                                                          
                            if (k == j)
                                break;
                            if (!(kv = ctx->PIdx[k])) {                                               
                                fin1 = 1;
                                xk = 1.0;
                            }
                            else  
                                xk = get_data(ctx, kv - 1,icase);

                            jk++;                                                         
                        }
                    }
                }
            }
            if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
                mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,icase);
            ctx->NOCUsed++;
        }
    }
    err = checkov(ctx);
    return(err);
}

/****************************************************************************/
/*  fn_exp1(resid)                                                          */
/*      Exponential Model with Time Periods.                                */
/*      Calculation of log-likelihood, gradient, hessian.                   */
/*      Episode splitting possible.                                         */
/*      If resid != 0 print residuals.                                      */
/*                                                                          */

int fn_exp1(TDAContext *ctx, int resid)
{
    register int j,k,jv,jk,kv,ip,tp,tps,tpf;
    int err,icase,sn,next,next1,org,des,cen,spl,nspl;
    double ts,tf,dt,xa,xj,xk,wt,extp,tmp = 0.0,res,rate = 0.0;

    wt = 1.0;  ctx->NOCUsed = 0;
           
    get_spell(ctx, 1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
    while (get_spell(ctx, 0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

        if (ctx->WIVar >= 0)
            wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;

        ip = ctx->TranPtr[(sn - 1) * ctx->MaxOrg + org];
        if (ip < 0)
            continue;
        next1 = ctx->PIdxPtr[ip] - 1;

        while (sn == ctx->SnTran1[ip] && org == ctx->OrgTran1[ip]) {

            if (des == ctx->DesTran1[ip])
                cen = 0;
            else
                cen = 1;
            ip++;
            next = next1;
            j = next + ctx->PMNTP;

            xa = 0.0;
            while ((jv = ctx->PIdx[++j]))                                                 
                xa += ctx->TPar[j] * get_data(ctx, jv - 1,icase);

            next1 = j - 1;
            res = 0.0;              /* accumulate residuals */

            /* get time periods */

            for (tps = 1; tps < ctx->PMNTP; ++tps) {
                if (ts < ctx->PMTP[tps])  
                    break;
            }
            for (tpf = 1; tpf < ctx->PMNTP; ++tpf) {
                if (tf < ctx->PMTP[tpf])  
                    break;
            }
            for (tp = tps; tp <= tpf; ++tp) {

                if (tps == tpf)
                    dt = tf - ts;
                else if (tp == tps)
                    dt = ctx->PMTP[tp] - ts;
                else if (tp == tpf)
                    dt = tf - ctx->PMTP[tp - 1];
                else
                    dt = ctx->PMTP[tp] - ctx->PMTP[tp - 1];

                extp = -dt * rexp(ctx, ctx->TPar[next + tp] + xa) * wt;

                if (resid) {
                    res -= extp;
                    if (tp == tpf)
                        rate = rexp(ctx, ctx->TPar[next + tp] + xa);
                }
                if (ctx->LFunc) {
                    tmp = extp;
                    if (tp == tpf && !cen)
                        tmp += (ctx->TPar[next + tp] + xa) * wt;

                    ctx->FTmp += tmp;
                }
                if (ctx->LGrad || ctx->LSec) {

                    j = next + ctx->PMNTP;
                    jk  = (next * (next - 1)) / 2;
                    jk += (ctx->PMNTP * (ctx->PMNTP - 1)) / 2;
                    jk += ctx->PMNTP * next;

                    if (ctx->LGrad) {
                        tmp = extp;
                        if (tp == tpf && !cen)
                            tmp += wt;
                        ctx->Grad[next + tp] += tmp;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + next + tp] = tmp;
                    }
                    if (ctx->LSec)  
                        ctx->Diag[next + tp] += extp;

                    while ((jv = ctx->PIdx[++j])) {

                        xj = get_data(ctx, jv - 1,icase);

                        if (ctx->LGrad) {
                            ctx->Grad[j] += tmp * xj;
                            if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                                ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp * xj;
                        }
                        if (ctx->LSec) {      

                            ctx->Diag[j] += extp * xj * xj;

                            k = next + ctx->PMNTP;
                            jk += next; 

                            ctx->Hess[jk + tp] += extp * xj;
                            jk += ctx->PMNTP;

                            while ((kv = ctx->PIdx[++k])) {
                                if (k == j)
                                    break;
                                xk = get_data(ctx, kv - 1,icase);
                                ctx->Hess[++jk] += extp * xj * xk;
                            }
                        }
                    }
                }
            }
            if (resid) {
                tmp = rexp(ctx, -res);
                prn_resid1(ctx, icase + 1,org,des,ts,tf,rate,tmp,res,wt);       
            }
            if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
                mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,icase);
            ctx->NOCUsed++;
        }
    }
    err = checkov(ctx);
    return(err);
}

/****************************************************************************/
/*  fn_exp2(resid)                                                          */
/*      Exponential Model with Time Periods, type II                        */
/*      Calculation of log-likelihood, gradient, hessian.                   */
/*      Episode splitting possible.                                         */
/*      If resid != 0 print residuals.                                      */
/*                                                                          */

int fn_exp2(TDAContext *ctx, int resid)
{
    register int j,k,jv,jk,kv,tp,tps,tpf;
    int err,ip,idx,icase,sn,spl,nspl,next,next1,org,des,cen;
    double ts,tf,dt,xa,xj,xk,wt,extp,tmp = 0.0,res,rate = 0.0;

    wt = 1.0; ctx->NOCUsed = 0;
    err = 0;
    get_spell(ctx, 1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
    while (get_spell(ctx, 0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

        if (ctx->WIVar >= 0)
            wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;

        ip = ctx->TranPtr[(sn - 1) * ctx->MaxOrg + org];
        if (ip < 0)
            continue;
        next1 = ctx->PIdxPtr[ip] - 1;

        while (sn == ctx->SnTran1[ip] && org == ctx->OrgTran1[ip]) {

            if (des == ctx->DesTran1[ip])
                cen = 0;
            else
                cen = 1;
            ip++;
            next = next1;
            j = idx = next + ctx->PMNTP;
            while (ctx->PIdx[++j]) ;
            next1 = j - 1;

            res = 0.0;              /* accumulate residuals */

            /* get time periods */

            for (tps = 1; tps < ctx->PMNTP; ++tps) {
                if (ts < ctx->PMTP[tps])  
                    break;
            }
            for (tpf = 1; tpf < ctx->PMNTP; ++tpf) {
                if (tf < ctx->PMTP[tpf])  
                    break;
            }

            for (tp = tps; tp <= tpf; ++tp) {

                if (tps == tpf)
                    dt = tf - ts;
                else if (tp == tps)
                    dt = ctx->PMTP[tp] - ts;
                else if (tp == tpf)
                    dt = tf - ctx->PMTP[tp - 1];
                else
                    dt = ctx->PMTP[tp] - ctx->PMTP[tp - 1];

                /* get X.alpha for time period tp */

                xa = 0.0;
                j  = idx + tp;
                while (j <= next1 && (jv = ctx->PIdx[j])) {                                               
                    xa += ctx->TPar[j] * get_data(ctx, jv - 1,icase);
                    j += ctx->PMNTP;
                }

                extp = -dt * rexp(ctx, ctx->TPar[next + tp] + xa) * wt;

                if (resid) {
                    res -= extp;
                    if (tp == tpf)
                        rate = rexp(ctx, ctx->TPar[next + tp] + xa);
                }

                if (ctx->LFunc) {
                    tmp = extp;
                    if (tp == tpf && !cen)
                        tmp += (ctx->TPar[next + tp] + xa) * wt;

                    ctx->FTmp += tmp;
                }
                    
                if (ctx->LGrad || ctx->LSec) {

                    if (ctx->LGrad) {
                        tmp = extp;
                        if (tp == tpf && !cen)
                            tmp += wt;
                        ctx->Grad[next + tp] += tmp;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + next + tp] = tmp;
                    }
                    if (ctx->LSec)  
                        ctx->Diag[next + tp] += extp;

                    j  = idx + tp;
                    while (j <= next1 && (jv = ctx->PIdx[j])) {

                        xj = get_data(ctx, jv - 1,icase);

                        if (ctx->LGrad) {
                            ctx->Grad[j] += tmp * xj;
                            if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                                ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp * xj;
                        }   
                        if (ctx->LSec) {      

                            ctx->Diag[j] += extp * xj * xj;

                            jk = next + tp + (j - 1) * (j - 2) / 2;
                                        
                            ctx->Hess[jk] += extp * xj;

                            k = idx + tp;
                            while (k <= next1 && (kv = ctx->PIdx[k])) {
                                if (k == j)
                                    break;
                                xk = get_data(ctx, kv - 1,icase);
                                jk += ctx->PMNTP;
                                ctx->Hess[jk] += extp * xj * xk;
                                k += ctx->PMNTP;
                            }
                        }
                        j += ctx->PMNTP;
                    }
                }
            }
            if (resid) {
                tmp = rexp(ctx, -res);
                prn_resid1(ctx, icase + 1,org,des,ts,tf,rate,tmp,res,wt);
            }
            if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
                mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,icase);
            ctx->NOCUsed++;
        }
    }
    err = checkov(ctx);
    return(err);
}

/****************************************************************************/
/*  fn_pol(resid)                                                           */
/*      Model with polynomial rates (I)                                     */
/*      Calculation of log-likelihood, gradient, hessian.                   */
/*      Episode splitting possible.                                         */
/*      If resid != 0 print residuals.                                      */

int fn_pol(TDAContext *ctx, int resid)
{
    register int j,k,jv,jk,kv,l,l1;
    int err,icase,sn,ip,next,org,des,cen,fin,fin1,spl,nspl;
    double ts,tf,xj,xk,wt,xa,dl,ds = 0.0,df,a,b,bb,ab,ab1,ab2,aab,ast,tmp1,tmp2,tmp3;
    double res;

    wt = 1.0;   ctx->NOCUsed = 0;
    err = 0;
    get_spell(ctx, 1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
    while (get_spell(ctx, 0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

        if (ctx->WIVar >= 0)                             /* get weights */
            wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;

        ip = ctx->TranPtr[(sn - 1) * ctx->MaxOrg1 + org];
        if (ip < 0)  
            continue;
 
        j = ctx->PIdxPtr[ip];

        while (sn == ctx->SnTran1[ip] && org == ctx->OrgTran1[ip]) {
            if (des == ctx->DesTran1[ip])
                cen = 0;
            else
                cen = 1;
            ip++;
            next = j;                                                       

            xa = ctx->TPar[j];                                                   
            while ((jv = ctx->PIdx[++j]))                                                 
                xa += ctx->TPar[j] * get_data(ctx, jv - 1,icase);
            a = rexp(ctx, xa);

            b = bb = 0.0;

            if (ctx->PMDEG) {
                l = 1;
                dl = 2.0;
                df = tf;
                ctx->POLWrk1[l] = df;
                b  += ctx->TPar[j] * ctx->POLWrk1[l];
                bb -= ctx->TPar[j] * ctx->POLWrk1[l] * tf / dl;

                if (ts > 0.0) {
                    ds = ts;
                    ctx->POLWrk2[l] = ds;
                    bb += ctx->TPar[j] * ctx->POLWrk2[l] * ts / dl;
                    ctx->POLWrk2[l] /= ctx->POLWrk1[l];
                }
                while ((jv = ctx->PIdx[++j])) {                                               
                    l++;
                    dl += 1.0;
                    df *= tf;
                    ctx->POLWrk1[l] = df;
                    b  += ctx->TPar[j] * ctx->POLWrk1[l];
                    bb -= ctx->TPar[j] * ctx->POLWrk1[l] * tf / dl;

                    if (ts > 0.0) {
                        ds *= ts;
                        ctx->POLWrk2[l] = ds;
                        bb += ctx->TPar[j] * ctx->POLWrk2[l] * ts / dl;
                        ctx->POLWrk2[l] /= ctx->POLWrk1[l];
                    }
                }
            }
            if (ts > 0.0)
                ast = a * (ts - tf);
            else
                ast = -a * tf;

            ab = a + b;
            ab1 = 1.0 / ab;
            ab2 = ab1 / ab;
            aab = a / ab;

            if (resid) {
                res = -ast - bb;
                tmp1 = rexp(ctx, -res);
                prn_resid1(ctx, icase + 1,org,des,ts,tf,ab,tmp1,res,wt);
            }

            if (ctx->LFunc) {

                tmp1 = ast + bb;
                if (!cen)
                    tmp1 += rlog(ctx, ab);
                tmp1 *= wt;

                ctx->FTmp += tmp1;                                   
            }
            if (ctx->LGrad || ctx->LSec) {
                j = next;
                jk = j * (j - 1) / 2;

                tmp2 = tmp1 = ast;
                if (!cen) {
                    tmp1 += aab;
                    tmp2 += aab * (1.0 - aab);
                }
                tmp1 *= wt;
                tmp2 *= wt;

                jv = 1;
                xj = 1.0;
                fin = 0;

                while (1) {
                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp1 * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp1 * xj;
                    }
                    if (ctx->LSec)  
                        ctx->Diag[j] += tmp2 * xj * xj;

                    if (!(jv = ctx->PIdx[++j])) {                                               
                        if (fin || !ctx->PMDEG)
                            break;
                        fin = 1;
                        l = 0;
                    }
                    if (fin) {
                        l++;
                        tmp1 = -tf;
                        tmp2 = tmp3 = 0.0;
                        if (ts > 0.0)  
                            tmp1 += ts * ctx->POLWrk2[l];

                        tmp1 /= (1.0 + (double)l);
                        if (!cen) {
                            tmp1 += ab1;
                            tmp2 -= ab2;
                            tmp3 -= a * ab2;
                        }
                        tmp1 *= wt;
                        tmp2 *= wt;
                        tmp3 *= wt;
                        xj = ctx->POLWrk1[l];
                    }
                    else
                        xj = get_data(ctx, jv - 1,icase);

                    if (ctx->LSec) {
                        k = next;                                                       
                        jk += next;                                                   
                        xk = 1.0;
                        fin1 = 0;
                        while (1) {
                            if (!fin || fin1)
                                ctx->Hess[jk] += tmp2 * xj * xk;
                            else
                                ctx->Hess[jk] += tmp3 * xj * xk;
                            k++;                                                          
                            if (k == j)
                                break;
                            if (!(kv = ctx->PIdx[k])) {                                               
                                fin1 = 1;
                                l1 = 0;
                            }
                            if (fin1)  
                                xk = ctx->POLWrk1[++l1];
                            else
                                xk = get_data(ctx, kv - 1,icase);
                            jk++;                                                         
                        }
                    }
                }
            }
            if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
                mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,icase);
            ctx->NOCUsed++;
        }
    }
    err = checkov(ctx);
    return(err);
}

/****************************************************************************/
/*  fn_pol1(resid)                                                          */
/*      Model with polynomial rates in an exponential function.             */
/*      Calculation of log-likelihood, gradient, hessematrix.               */
/*      Episode splitting possible.                                         */
/*      If resid != 0 print residuals.                                      */
/*                                                                          */
/*      Numerical integration is done with ni_gen() depending on PMM.       */
/*                                                                          */

int fn_pol1(TDAContext *ctx, int resid)
{
    register int i,ii,j,k,jv,jk,kv;
    int err,icase,sn,ip,spl,nspl,next,nexta,nextb,org,des,cen;
    double wt,ts,tf,xa,a,b,res,rate;
    double dl,tmp,tmp1,tmp2,xj,xk,gint,tgsum;

    wt = 1.0; ctx->NOCUsed = 0;
    err = 0;
    get_spell(ctx, 1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
    while (get_spell(ctx, 0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

        if (ctx->WIVar >= 0)                             /* get weights */
            wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;

        ip = ctx->TranPtr[(sn - 1) * ctx->MaxOrg1 + org];
        if (ip < 0)  
            continue;
 
        j = ctx->PIdxPtr[ip];

        while (sn == ctx->SnTran1[ip] && org == ctx->OrgTran1[ip]) {

            if (des == ctx->DesTran1[ip])
                cen = 0;
            else
                cen = 1;
            ip++;

            nexta = j;
            next = j - 1;                                                   
            
            xa = ctx->TPar[j];                                                   
            while ((jv = ctx->PIdx[++j]))                                                 
                xa += ctx->TPar[j] * get_data(ctx, jv - 1,icase);
            nextb = j;
            while ((jv = ctx->PIdx[++j])) ;     /* now we have j for the next    */       
                                           /* transition.                   */
            a = -rexp(ctx, xa);
            tgsum = 0.0;
            if (!ctx->PMDEG) {                   /* in this case we have a simple */
                gint = tf - ts;            /* exponential model.            */
            }
            else {
                dl = 1.0;
                for (i = 0; i < ctx->PMDEG; ++i) {
                    b = ctx->TPar[nextb + i];
                    ctx->BCBeta[i] = b;
                    ctx->BCTG[i] = pow(tf,dl);
                    tgsum += b * ctx->BCTG[i];                  
                    dl += 1.0;
                }
                gint = ni_gen(ctx, ts,tf,2,0,0, &err);
            }

            if (resid) {
                res = -a * gint;
                tmp1 = rexp(ctx, -res);
                rate = rexp(ctx, xa + tgsum);
                prn_resid1(ctx, icase + 1,org,des,ts,tf,rate,tmp1,res,wt);
            }

            if (ctx->LFunc) {

                tmp = a * gint; 
                if (!cen)
                    tmp += xa + tgsum;
                tmp *= wt;

                ctx->FTmp += tmp;                                     
            }
            if (ctx->LGrad || ctx->LSec) {
                                            /* next beginnt mit 0 */
                j = nexta;                  /* nexta mit 1        */

                jk = (j * (j - 1)) / 2 - next;
                xj = 1.0;   

                tmp1 = tmp = a * gint;
                if (!cen) 
                    tmp += 1.0;

                tmp  *= wt;
                tmp1 *= wt;

                while (1) {         /* alpha j */
                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;
     
                        while (1) {     /* alpha k */
                            if (k < j)  
                                ctx->Hess[++jk] += tmp1 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp1 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
                i = 0;
                while (ctx->PMDEG) {        /* beta j */

                    ctx->MPOL1I = i + 1;
                    tmp1 = a * ni_gen(ctx, ts,tf,3,0,0,&err);

                    if (ctx->LGrad) {
                        tmp = tmp1;
                        if (!cen)  
                            tmp += ctx->BCTG[i];              

                        ctx->Grad[j] += wt * tmp;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = wt * tmp;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (1) {         /* alpha k */

                            ctx->Hess[++jk] += wt * tmp1 * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextb;
                        ii = 0;
                        while (1) {         /* beta k */

                            ctx->MPOL1J = ii + 1;
                            tmp2 = ni_gen(ctx, ts,tf,4,0,0,&err);
                            tmp2 *= wt * a;

                            if (k < j)  
                                ctx->Hess[++jk] += tmp2;
                            else if (k == j) {
                                ctx->Diag[k] += tmp2;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            ii++;
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    i++;
                }
            }
            if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
                mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,icase);
            ctx->NOCUsed++;
        }
    }
    if (err == 0)
        err = checkov(ctx);
    return(err);
}

/****************************************************************************/
/*  get_mpol1(t,err)                                                        */
/*  Get function value for numerical integration.                           */

double get_mpol1(TDAContext *ctx, double t,int *err)
{
    register int i;
    double f;

    *err = 0;
    f = 0.0;
    for (i = 0; i < ctx->PMDEG; ++i) {
        f += ctx->BCBeta[i] * t;
        t *= t;
    }
    f = rexp(ctx, f);
    return(f);
}

/****************************************************************************/
/*  get_mpol2(t,err)                                                        */
/*  Get function value for numerical integration.                           */
/*  First derivatives according to MPOL1I                                   */

double get_mpol2(TDAContext *ctx, double t,int *err)
{
    double f;

    f = get_mpol1(ctx, t,err);
    f *= pow(t,(double)ctx->MPOL1I);
    return(f);
}

/****************************************************************************/
/*  get_mpol3(t,err)                                                        */
/*  Get function value for numerical integration.                           */
/*  Second derivatives according to MPOL1I and MPOL1J.                      */

double get_mpol3(TDAContext *ctx, double t,int *err)
{
    double f;

    f = get_mpol1(ctx, t,err);
    f *= pow(t,(double)(ctx->MPOL1I + ctx->MPOL1J));
    return(f);
}

/****************************************************************************/
/*  fn_gm(resid)                                                            */
/*      Gompertz-Makeham                                                    */
/*      Calculation of log-likelihood, gradient, hessian                    */
/*      Episode splitting possible.                                         */
/*      If resid != 0 print residuals.                                      */

int fn_gm(TDAContext *ctx, int resid)
{
    register int j,k,jv,jk,kv;
    int err,icase,sn,ip,spl,nspl,next,nexta,nextb,nextc,org,des,cen;
    double wt,ts,tf,a,b,c,c1,bc,ast,ect,ecs,ecst1,ecst2,bcecs;
    double abect,aabect,babect,tmp,tmp1,tmp2,xj,xk;

    wt = 1.0;   ctx->NOCUsed = 0;
    err = 0;    
    get_spell(ctx, 1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
    while (get_spell(ctx, 0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

        if (ctx->WIVar >= 0)                             /* get weights */
            wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;

        ip = ctx->TranPtr[(sn - 1) * ctx->MaxOrg1 + org];
        if (ip < 0)  
            continue;
 
        j = ctx->PIdxPtr[ip];

        while (sn == ctx->SnTran1[ip] && org == ctx->OrgTran1[ip]) {

            if (des == ctx->DesTran1[ip])
                cen = 0;
            else
                cen = 1;
            ip++;

            nexta = j;
            next = j - 1;                                                   

            a = b = c = 0.0;

            if (ctx->XMFlg[1]) {
                a = ctx->TPar[j];                                                   
                while ((jv = ctx->PIdx[++j]))                                                 
                    a += ctx->TPar[j] * get_data(ctx, jv - 1,icase);
                a = rexp(ctx, a);
            }
            if (ctx->XMFlg[2]) {
                b = ctx->TPar[j];                                                   
                while ((jv = ctx->PIdx[++j]))                                                 
                    b += ctx->TPar[j] * get_data(ctx, jv - 1,icase);
                b = rexp(ctx, b);
                c = ctx->TPar[j];                                                   
                while ((jv = ctx->PIdx[++j]))                                                 
                    c += ctx->TPar[j] * get_data(ctx, jv - 1,icase);
            }
            ast = a * (ts - tf);

            if (fabs(c) > ctx->EPSI) {
                ect = rexp(ctx, c * tf);
                ecs = rexp(ctx, c * ts);
                bc = b / c;
                c1 = 1.0 / c;
                tmp = c1 / c;
                tmp1 = ts - c1;
                tmp2 = tf - c1;
                bcecs = bc * (ecs - ect);
                ecst1 = ecs * tmp1 - ect * tmp2;
                ecst2 = ecs * (tmp1 * tmp1 + tmp) - ect * (tmp2 * tmp2 + tmp);
            }
            else {
                ecs = ect = 1.0;
                ecst2 = ecst1 = bcecs = bc = 0.0;
            }
            abect = a + b * ect;
            aabect = a / abect;
            babect = b * ect / abect;

            if (resid) {
                tmp1 = rexp(ctx, ast + bcecs);
                tmp2 = -ast - bcecs;
                prn_resid1(ctx, icase + 1,org,des,ts,tf,abect,tmp1,tmp2,wt);
            }

            if (ctx->LFunc) {
                tmp = ast + bcecs;    
                if (!cen)
                    tmp += rlog(ctx, abect);
                tmp *= wt;

                ctx->FTmp += tmp;                                     
            }
            if (ctx->LGrad || ctx->LSec) {
                                            /* next beginnt mit 0 */
                j = nexta;                  /* nexta mit 1        */

                jk = (j * (j - 1)) / 2 - next;
                xj = 1.0;   

                while (ctx->XMFlg[1]) {         /*  alpha j  */
                    if (ctx->LGrad) {
                        tmp = ast;
                        if (!cen)
                            tmp += aabect;
                        tmp *= wt;
                        ctx->Grad[j] += tmp * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;
     
                        while (1) {     /* alpha k */

                            tmp = ast;
                            if (!cen)
                                tmp += aabect * (1.0 - aabect);
                            tmp *= wt;

                            if (k < j)  
                                ctx->Hess[++jk] += tmp * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
                nextb = j;
                xj = 1.0;   

                while (ctx->XMFlg[2]) {     /* beta j */
                    if (ctx->LGrad) {
                        tmp = bcecs;
                        if (!cen)
                            tmp += babect;
                        tmp *= wt;
                        ctx->Grad[j] += tmp * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (ctx->XMFlg[1]) {         /* alpha k */

                            tmp = 0.0;
                            if (!cen)
                                tmp -= aabect * babect;
                            tmp *= wt;

                            ctx->Hess[++jk] += tmp * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextb;
                        xk = 1.0;

                        while (ctx->XMFlg[2]) {         /* beta k  */

                            tmp = bcecs;
                            if (!cen)
                                tmp += babect * (1.0 - babect);
                            tmp *= wt;
    
                            if (k < j)  
                                ctx->Hess[++jk] += tmp * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
                nextc = j;
                xj = 1.0;   

                while (ctx->XMFlg[3]) {             /* gamma j */
                    if (ctx->LGrad) {
                        tmp = bc * ecst1;
                        if (!cen)
                            tmp += tf * babect;
                        tmp *= wt;
                        ctx->Grad[j] += tmp * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (ctx->XMFlg[1]) {         /* alpha k */

                            tmp = 0.0;
                            if (!cen)
                                tmp -= tf * aabect * babect;
                            tmp *= wt;

                            ctx->Hess[++jk] += tmp * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextb;
                        xk = 1.0;

                        while (ctx->XMFlg[2]) {         /* beta k */

                            tmp = bc * ecst1;
                            if (!cen)
                                tmp += tf * babect * (1.0 - babect);
                            tmp *= wt;

                            ctx->Hess[++jk] += tmp * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextc;
                        xk = 1.0;

                        while (ctx->XMFlg[3]) {         /* gamma k */

                            tmp = bc * ecst2;
                            if (!cen)
                                tmp += tf * tf * babect * (1.0 - babect);
                            tmp *= wt;

                            if (k < j)  
                                ctx->Hess[++jk] += tmp * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
            }
            if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
                mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,icase);
            ctx->NOCUsed++;
        }
    }
    err = checkov(ctx);
    return(err);
}

/****************************************************************************/
/*  fn_wei(resid)                                                           */
/*      Weibull model, optional gamma mixture                               */
/*      Calculation of log-likelihood, gradient, hessian.                   */
/*      Episode splitting possible if no mixture.                           */
/*      If resid != 0 print residuals.                                      */

int fn_wei(TDAContext *ctx, int resid)
{
    register int j,k,jv,jk,kv;
    int err,icase,sn,ip,spl,nspl,next,nexta,nextb,nextc,org,des,cen;
    double ts,tf,xa,xb,xd,xj,xk,wt,b,d = 0.0,las = 0.0,lat,las1 = 0.0,lat1,hc,hcs = 0.0,ldh = 0.0;
    double tmp,tmp1,tmp2,tmp3,tmp4,ha,haa,hb,hbb,hab,lts = 0.0,ltf,ui = 0.0,dh = 0.0,dh1 = 0.0;

    wt = 1.0;   ctx->NOCUsed = 0;
    err = 0;
    get_spell(ctx, 1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
    while (get_spell(ctx, 0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

        if (ctx->WIVar >= 0)                             /* get weights */
            wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;

        ip = ctx->TranPtr[(sn - 1) * ctx->MaxOrg1 + org];
        if (ip < 0)  
            continue;
 
        j = ctx->PIdxPtr[ip];

        ltf = rlog(ctx, tf);
        if (ts > 0.0)
            lts = rlog(ctx, ts);

        while (sn == ctx->SnTran1[ip] && org == ctx->OrgTran1[ip]) {

            if (des == ctx->DesTran1[ip])
                cen = 0;
            else
                cen = 1;
            ip++;

            nexta = j;
            next = j - 1;                                                   

            xa = ctx->TPar[j];                                                   
            while ((jv = ctx->PIdx[++j]))                                                 
                xa += ctx->TPar[j] * get_data(ctx, jv - 1,icase);
            xb = ctx->TPar[j];                                                   
            while ((jv = ctx->PIdx[++j]))                                                 
                xb += ctx->TPar[j] * get_data(ctx, jv - 1,icase);

            b = rexp(ctx, xb);
            lat = xa + ltf;
            lat1 = 1.0 + b * lat;
            hc = rexp(ctx, b * lat);

            if (ts > 0.0) {
                las = xa + lts;
                las1 = 1.0 + b * las;
                hcs = rexp(ctx, b * las);
            }
            if (resid) {
                               
                tmp = xb + b * xa + (b - 1.0) * ltf;
                tmp1 = rexp(ctx, tmp);
                tmp = hc;
                if (ts > 0.0)
                    tmp -= hcs;
                tmp2 = rexp(ctx, -tmp);
                prn_resid1(ctx, icase + 1,org,des,ts,tf,tmp1,tmp2,tmp,wt);
            }

            if (ctx->MixTyp) {

                xd = ctx->TPar[j];                                                   
                while ((jv = ctx->PIdx[++j]))                                                 
                    xd += ctx->TPar[j] * get_data(ctx, jv - 1,icase);

                d = rexp(ctx, xd);
                ui = -1.0 / d;
                if (!cen)
                    ui -= 1.0;

                if (hc > 1.e50)
                    ldh = xd + b * lat;
                else
                    ldh = rlog(ctx, 1.0 + d * hc);

                dh = d / (1.0 + d * hc);
                dh1 = dh * hc;
            }
            if (ctx->LFunc) {

                if (ctx->MixTyp) 
                    tmp = ui * ldh;
                else {
                    tmp = -hc;
                    if (ts > 0.0)
                        tmp += hcs;
                }
                if (!cen)  
                    tmp += xb + b * xa + (b - 1.0) * ltf;

                tmp *= wt;

                ctx->FTmp += tmp;                                     
            }
            if (ctx->LGrad || ctx->LSec) {
                                            /* next beginnt mit 0 */
                j = nexta;                  /* nexta mit 1        */

                jk = (j * (j - 1)) / 2 - next;
                xj = 1.0;   

                ha = hc;
                hb = hc * lat;
                hab = hc * lat1;
                hbb = lat * hab;

                if (ts > 0.0) {
                    ha -= hcs;
                    hb -= hcs * las;
                    hab -= hcs * las1;
                    hbb -= las * hcs * las1;
                }
                ha *= b;
                hb *= b;
                hab *= b;
                hbb *= b;
                haa = ha * b;

                if (ctx->MixTyp) {
                    tmp1 = ui * dh * ha;
                    tmp2 = ui * dh * (haa - dh * ha * ha);
                }
                else {
                    tmp1 = -ha;
                    tmp2 = -haa;
                }
                if (!cen)
                    tmp1 += b;

                tmp1 *= wt;
                tmp2 *= wt;

                while (1) {         /* Alpha j */

                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp1 * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp1 * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (1) {     /* Alpha k */

                            if (k < j)  
                                ctx->Hess[++jk] += tmp2 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp2 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
                nextb = j;
                xj = 1.0;   

                if (ctx->MixTyp) {
                    tmp1 = ui * dh * hb;
                    tmp2 = ui * dh * (hab - dh * ha * hb);
                    tmp3 = ui * dh * (hbb - dh * hb * hb);
                }
                else {
                    tmp1 = -hb;
                    tmp2 = -hab;
                    tmp3 = -hbb;
                }
                if (!cen) {
                    tmp1 += lat1;
                    tmp2 += b;
                    tmp3 += b * lat;
                }
                tmp1 *= wt;
                tmp2 *= wt;
                tmp3 *= wt;

                while (1) {     /* Beta j */

                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp1 * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp1 * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (1) {         /* Alpha k */

                            ctx->Hess[++jk] += tmp2 * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextb;
                        xk = 1.0;

                        while (1) {         /* Beta k */
                            if (k < j)  
                                ctx->Hess[++jk] += tmp3 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp3 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
                nextc = j;
                xj = 1.0;   

                if (ctx->MixTyp) {
                    tmp1 = wt * (ldh / d + ui * dh1);
                    tmp2 = tmp3 = wt * dh * (1.0 / d + ui * (1.0 - dh1));
                    tmp2 *= ha;
                    tmp3 *= hb;
                    tmp4 = wt * ((2.0 * dh1 - ldh) / d + ui * dh1 * (1.0 - dh1));
                }

                while (ctx->MixTyp) {             /* Delta j */

                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp1 * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp1 * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (1) {         /* Alpha k */

                            ctx->Hess[++jk] += tmp2 * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextb;
                        xk = 1.0;

                        while (1) {         /* Beta k */

                            ctx->Hess[++jk] += tmp3 * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextc;
                        xk = 1.0;

                        while (1) {         /* Delta k */
                            if (k < j)  
                                ctx->Hess[++jk] += tmp4 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp4 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
            }
            if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
                mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,icase);
            ctx->NOCUsed++;
        }
    }
    err = checkov(ctx);
    return(err);
}

/****************************************************************************/
/*  fn_sic(resid)                                                           */
/*      Sickle model, optional gamma mixture                                */
/*      Calculation of log-likelihood, gradient, hessian.                   */
/*      Episode splitting possible if no mixture.                           */
/*      If resid != 0 print residuals.                                      */

int fn_sic(TDAContext *ctx, int resid)
{
    register int j,k,jv,jk,kv;
    int err,icase,sn,ip,spl,nspl,next,nexta,nextb,nextc,org,des,cen;
    double ts,tf,xa,xb,xd,xj,xk,wt,a,b,d = 0.0,hc,hcs,ldh = 0.0,ab,sb,tb,tsb,tfb,etb,esb;
    double tmp,tmp1,tmp2,tmp3,tmp4,ha,haa,hb,hbb,hab,ltf,ui = 0.0,dh = 0.0,dh1 = 0.0;

    wt = 1.0;   ctx->NOCUsed = 0;
    err = 0;
    get_spell(ctx, 1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
    while (get_spell(ctx, 0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

        if (ctx->WIVar >= 0)                             /* get weights */
            wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;

        ip = ctx->TranPtr[(sn - 1) * ctx->MaxOrg1 + org];
        if (ip < 0)  
            continue;
 
        ltf = rlog(ctx, tf);

        j = ctx->PIdxPtr[ip];

        while (sn == ctx->SnTran1[ip] && org == ctx->OrgTran1[ip]) {

            if (des == ctx->DesTran1[ip])
                cen = 0;
            else
                cen = 1;
            ip++;

            nexta = j;
            next = j - 1;                                                   

            xa = ctx->TPar[j];                                                   
            while ((jv = ctx->PIdx[++j]))                                                 
                xa += ctx->TPar[j] * get_data(ctx, jv - 1,icase);
            xb = ctx->TPar[j];                                                   
            while ((jv = ctx->PIdx[++j]))                                                 
                xb += ctx->TPar[j] * get_data(ctx, jv - 1,icase);

            a = rexp(ctx, xa);
            b = rexp(ctx, xb);

            ab = a * b;
            tb = tf + b;
            tfb = tf / b;
            etb = rexp(ctx, -tfb);
            hc  = ab * (b - tb * etb);
            if (ts > 0.0) {
                sb = ts + b;
                tsb = ts / b;
                esb = rexp(ctx, -tsb);
                hcs = ab * (b - sb * esb);
            }
            else {
                sb = b;
                tsb = 0.0;
                esb = 1.0;
                hcs = 0.0;
            }
            if (resid) {
                               
                tmp1 = rexp(ctx, ltf + xa - tfb);
                tmp = hc - hcs;
                tmp2 = rexp(ctx, -tmp);
                prn_resid1(ctx, icase + 1,org,des,ts,tf,tmp1,tmp2,tmp,wt);
            }

            if (ctx->MixTyp) {

                xd = ctx->TPar[j];                                                   
                while ((jv = ctx->PIdx[++j]))                                                 
                    xd += ctx->TPar[j] * get_data(ctx, jv - 1,icase);

                d = rexp(ctx, xd);
                ui = -1.0 / d;
                if (!cen)
                    ui -= 1.0;

                ldh = rlog(ctx, 1.0 + d * hc);
                dh = d / (1.0 + d * hc);
                dh1 = dh * hc;
            }
            if (ctx->LFunc) {

                if (ctx->MixTyp) 
                    tmp = ui * ldh;
                else  
                    tmp = hcs - hc;
                if (!cen)  
                    tmp += ltf + xa - tfb;

                tmp *= wt;

                ctx->FTmp += tmp;                                     
            }
            if (ctx->LGrad || ctx->LSec) {
                                            /* next starts with 0 */
                j = nexta;                  /* nexta wiht 1       */

                jk = (j * (j - 1)) / 2 - next;
                xj = 1.0;   

                ha = hc - hcs;
                hb = ab * (esb * (2.0 * sb + ts * tsb) - 
                           etb * (2.0 * tb + tf * tfb)); 
                haa = ha;
                hbb = ab * (esb * (4.0 * sb + ts * ts * (2.0 + tsb) / b) -
                            etb * (4.0 * tb + tf * tf * (2.0 + tfb) / b));

                hab = hb;

                if (ctx->MixTyp) {
                    tmp1 = ui * dh * ha;
                    tmp2 = ui * dh * (haa - dh * ha * ha);
                }
                else {
                    tmp1 = -ha;
                    tmp2 = -haa;
                }
                if (!cen)  
                    tmp1 += 1.0;

                tmp1 *= wt;
                tmp2 *= wt;

                while (1) {         /* Alpha j */

                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp1 * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp1 * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (1) {     /* Alpha k */

                            if (k < j)  
                                ctx->Hess[++jk] += tmp2 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp2 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
                nextb = j;
                xj = 1.0;   

                if (ctx->MixTyp) {
                    tmp1 = ui * dh * hb;
                    tmp2 = ui * dh * (hab - dh * ha * hb);
                    tmp3 = ui * dh * (hbb - dh * hb * hb);
                }
                else {
                    tmp1 = -hb;
                    tmp2 = -hab;
                    tmp3 = -hbb;
                }
                if (!cen) {
                    tmp1 += tfb;
                    tmp3 -= tfb;
                }
                tmp1 *= wt;
                tmp2 *= wt;
                tmp3 *= wt;

                while (1) {     /* Beta j */

                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp1 * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp1 * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (1) {         /* Alpha k */

                            ctx->Hess[++jk] += tmp2 * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextb;
                        xk = 1.0;

                        while (1) {         /* Beta k */
                            if (k < j)  
                                ctx->Hess[++jk] += tmp3 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp3 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
                nextc = j;
                xj = 1.0;   

                if (ctx->MixTyp) {
                    tmp1 = wt * (ldh / d + ui * dh1);
                    tmp2 = tmp3 = wt * dh * (1.0 / d + ui * (1.0 - dh1));
                    tmp2 *= ha;
                    tmp3 *= hb;
                    tmp4 = wt * ((2.0 * dh1 - ldh) / d + ui * dh1 * (1.0 - dh1));
                }

                while (ctx->MixTyp) {             /* Delta j */

                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp1 * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp1 * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (1) {         /* Alpha k */

                            ctx->Hess[++jk] += tmp2 * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextb;
                        xk = 1.0;

                        while (1) {         /* Beta k */

                            ctx->Hess[++jk] += tmp3 * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextc;
                        xk = 1.0;

                        while (1) {         /* Delta k */
                            if (k < j)  
                                ctx->Hess[++jk] += tmp4 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp4 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
            }
            if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
                mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,icase);
            ctx->NOCUsed++;
        }
    }
    err = checkov(ctx);
    return(err);
}

/****************************************************************************/
/*  fn_ll(resid)                                                            */
/*      Log-logistic model type I, optional gamma mixture                   */
/*      Calculation of log-likelihood, gradient, hessian.                   */
/*      Episode splitting possible if no mixture.                           */
/*      If resid != 0 print residuals.                                      */

int fn_ll(TDAContext *ctx, int resid)
{
    register int j,k,jv,jk,kv;
    int err,icase,sn,ip,spl,nspl,next,nexta,nextb,nextc,org,des,cen;
    double ts,tf,xa,xb,xd,xj,xk,wt,b,d = 0.0,hc,hcs = 0.0,ldh = 0.0;
    double atb,atb1,asb = 0.0,asb1 = 0.0,lat,las = 0.0,atd,asd;
    double tmp,tmp1,tmp2,tmp3,tmp4,ha,haa,hb,hbb,hab,lts = 0.0,ltf,ui = 0.0,dh = 0.0,dh1 = 0.0;

    wt = 1.0;   ctx->NOCUsed++;
    err = 0; 
    get_spell(ctx, 1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
    while (get_spell(ctx, 0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

        if (ctx->WIVar >= 0)                             /* get weights */
            wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;

        ip = ctx->TranPtr[(sn - 1) * ctx->MaxOrg1 + org];
        if (ip < 0)  
            continue;
 
        ltf = rlog(ctx, tf);
        if (ts > 0.0)
            lts = rlog(ctx, ts);

        j = ctx->PIdxPtr[ip];

        while (sn == ctx->SnTran1[ip] && org == ctx->OrgTran1[ip]) {

            if (des == ctx->DesTran1[ip])
                cen = 0;
            else
                cen = 1;
            ip++;

            nexta = j;
            next = j - 1;                                                   

            xa = ctx->TPar[j];                                                   
            while ((jv = ctx->PIdx[++j]))                                                 
                xa += ctx->TPar[j] * get_data(ctx, jv - 1,icase);
            xb = ctx->TPar[j];                                                   
            while ((jv = ctx->PIdx[++j]))                                                 
                xb += ctx->TPar[j] * get_data(ctx, jv - 1,icase);

            b = rexp(ctx, xb);

            lat = xa + ltf;
            /*******************  
            atb = pow(a * tf,b);
            ********************/
            atb = rexp(ctx, b * (xa + ltf));
            atb1 = atb + 1.0;
            if (atb1 > 1.e50)
                hc = b * (xa + ltf);
            else
                hc = rlog(ctx, atb1);

            if (ts > 0.0) {
                las = xa + lts;
                /*******************  
                asb = pow(a * ts,b);
                ********************/
                asb = rexp(ctx, b * (xa + lts));
                asb1 = asb + 1.0;
                if (asb1 > 1.e50)
                    hcs = b * (xa + lts);
                else
                    hcs = rlog(ctx, asb1);
            }
            if (ctx->MixTyp) {

                xd = ctx->TPar[j];                                                   
                while ((jv = ctx->PIdx[++j]))                                                 
                    xd += ctx->TPar[j] * get_data(ctx, jv - 1,icase);

                d = rexp(ctx, xd);
                ui = -1.0 / d;
                if (!cen)
                    ui -= 1.0;

                ldh = rlog(ctx, 1.0 + d * hc);
                dh = d / (1.0 + d * hc);
                dh1 = dh * hc;
            }
            if (resid) {
                                
                tmp = xb + b * xa + (b - 1.0) * ltf - hc;
                tmp1 = rexp(ctx, tmp);

                tmp = hc;
                if (ts > 0.0)
                    tmp -= hcs;
                tmp2 = rexp(ctx, -tmp);

                prn_resid1(ctx, icase + 1,org,des,ts,tf,tmp1,tmp2,tmp,wt);
            }

            if (ctx->LFunc) {

                if (ctx->MixTyp) 
                    tmp = ui * ldh;
                else {
                    tmp = -hc;
                    if (ts > 0.0)
                        tmp += hcs;
                }
                if (!cen)  
                    tmp += xb + b * xa + (b - 1.0) * ltf - hc;

                tmp *= wt;

                ctx->FTmp += tmp;                                     
            }
            if (ctx->LGrad || ctx->LSec) {
                                            /* next starts with 0 */
                j = nexta;                  /* nexta wiht 1       */

                jk = (j * (j - 1)) / 2 - next;
                xj = 1.0;   

                atd = atb / atb1;
                if (ts > 0.0)
                    asd = asb / asb1;

                ha = atd;
                hb = lat * atd;
                haa = atd * (1.0 - atd);
                hbb = lat * atd * (1.0 + b * lat / atb1);
                hab = atd * (1.0 + b * lat / atb1);

                if (ts > 0.0) {
                    ha -= asd;
                    hb -= las * asd;
                    haa -= asd * (1.0 - asd);
                    hbb -= las * asd * (1.0 + b * las / asb1);
                    hab -= asd * (1.0 + b * las / asb1);
                }
                ha *= b;
                hb *= b;
                haa *= b * b;
                hbb *= b;
                hab *= b;

                if (ctx->MixTyp) {
                    tmp1 = ui * dh * ha;
                    tmp2 = ui * dh * (haa - dh * ha * ha);
                }
                else {
                    tmp1 = -ha;
                    tmp2 = -haa;
                }
                if (!cen) {
                    tmp1 += b * (1.0 - atd);
                    tmp2 += b * b * atd * (atd - 1.0);
                }
                tmp1 *= wt;
                tmp2 *= wt;

                while (1) {         /* Alpha j */

                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp1 * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp1 * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (1) {     /* Alpha k */

                            if (k < j)  
                                ctx->Hess[++jk] += tmp2 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp2 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
                nextb = j;
                xj = 1.0;   

                if (ctx->MixTyp) {
                    tmp1 = ui * dh * hb;
                    tmp2 = ui * dh * (hab - dh * ha * hb);
                    tmp3 = ui * dh * (hbb - dh * hb * hb);
                }
                else {
                    tmp1 = -hb;
                    tmp2 = -hab;
                    tmp3 = -hbb;
                }
                if (!cen) {
                    tmp1 += 1.0 + b * lat / atb1;
                    tmp2 += b * (1.0 - atd * (1.0 + b * lat / atb1));
                    tmp3 += b * (lat / atb1) * (1.0 - b * lat * atd);
                }
                tmp1 *= wt;
                tmp2 *= wt;
                tmp3 *= wt;

                while (1) {     /* Beta j */

                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp1 * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp1 * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (1) {         /* Alpha k */

                            ctx->Hess[++jk] += tmp2 * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextb;
                        xk = 1.0;

                        while (1) {         /* Beta k */
                            if (k < j)  
                                ctx->Hess[++jk] += tmp3 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp3 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
                nextc = j;
                xj = 1.0;   

                if (ctx->MixTyp) {
                    tmp1 = wt * (ldh / d + ui * dh1);
                    tmp2 = tmp3 = wt * dh * (1.0 / d + ui * (1.0 - dh1));
                    tmp2 *= ha;
                    tmp3 *= hb;
                    tmp4 = wt * ((2.0 * dh1 - ldh) / d + ui * dh1 * (1.0 - dh1));
                }

                while (ctx->MixTyp) {             /* Delta j */

                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp1 * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp1 * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (1) {         /* Alpha k */

                            ctx->Hess[++jk] += tmp2 * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextb;
                        xk = 1.0;

                        while (1) {         /* Beta k */

                            ctx->Hess[++jk] += tmp3 * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextc;
                        xk = 1.0;

                        while (1) {         /* Delta k */
                            if (k < j)  
                                ctx->Hess[++jk] += tmp4 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp4 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
            }
            if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
                mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,icase);
            ctx->NOCUsed++;
        }
    }
    err = checkov(ctx);
    return(err);
}

/****************************************************************************/
/*  fn_ll2(resid)                                                           */
/*      Three parameter log-logistic model, type II.                        */
/*      Calculation of log-likelihood, gradient, hessian                    */
/*      If resid != 0 print residuals.                                      */

int fn_ll2(TDAContext *ctx, int resid)
{
    register int j,k,jv,jk,kv;
    int err,icase,sn,ip,spl,nspl,next,nexta,nextb,nextc,org,des,cen;
    double wt,ts,tf,tslog,tflog,xa,xb,xc,lam,del,gam,glam,zs = 0.0,zf,ezs = 0.0,ezf,lts,ltf;
    double ezs1 = 0.0,ezf1,ezf2,zs1 = 0.0,zf1,llts = 0.0,lltf,tmp,tmp1,tmp2,tmp3,xj,xk;

    wt = 1.0;   ctx->NOCUsed++;
    err = 0;
    get_spell(ctx, 1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
    while (get_spell(ctx, 0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

        if (ctx->WIVar >= 0)                             /* get weights */
            wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;

        ip = ctx->TranPtr[(sn - 1) * ctx->MaxOrg1 + org];
        if (ip < 0)  
            continue;
 
        j = ctx->PIdxPtr[ip];

        while (sn == ctx->SnTran1[ip] && org == ctx->OrgTran1[ip]) {

            if (des == ctx->DesTran1[ip])
                cen = 0;
            else
                cen = 1;
            ip++;

            nexta = j;
            next = j - 1;                                                   
            xa = ctx->TPar[j];                                                   
            while ((jv = ctx->PIdx[++j]))                                                 
                xa += ctx->TPar[j] * get_data(ctx, jv - 1,icase);
            xb = ctx->TPar[j];                                                   
            while ((jv = ctx->PIdx[++j]))                                                 
                xb += ctx->TPar[j] * get_data(ctx, jv - 1,icase);
            xc = ctx->TPar[j];                                                   
            while ((jv = ctx->PIdx[++j]))                                                 
                xc += ctx->TPar[j] * get_data(ctx, jv - 1,icase);

            lam = rexp(ctx, xa);
            del = rexp(ctx, xb);
            gam = rexp(ctx, xc);
            glam = gam / lam;

            tflog = rlog(ctx, tf);
            zf = del * (tflog + xa);
            ezf = rexp(ctx, zf);
            ezf1 = ezf / (1.0 + ezf);
            ezf2 = ezf1 / (1.0 + ezf);
            zf1  = zf / (1.0 + ezf);
            ltf  = 1.0 + pow(tf * lam,del);
            lltf = rlog(ctx, ltf);

            if (ts > 0.0) {
                tslog = rlog(ctx, ts);
                zs = del * (tslog + xa);
                ezs = rexp(ctx, zs);
                ezs1 = ezs / (1.0 + ezs);
                zs1  = zs / (1.0 + ezs);
                lts = 1.0 + pow(ts * lam,del);
                llts = rlog(ctx, lts);
            }

            if (resid) {
                tmp  = xc + xb - xa + zf - tflog - lltf;
                tmp1 = rexp(ctx, tmp);

                tmp = lltf;
                if (ts > 0.0)
                    tmp -= llts;
                tmp *= glam;
                tmp2 = rexp(ctx, -tmp);

                prn_resid1(ctx, icase + 1,org,des,ts,tf,tmp1,tmp2,tmp,wt);
            }

            if (ctx->LFunc) {

                tmp = -lltf;
                if (ts > 0.0)
                    tmp += llts;
                tmp *= glam;
                if (!cen)
                    tmp += xc + xb - xa + zf - tflog - lltf;
                tmp *= wt;

                ctx->FTmp += tmp;                                     
            }
            if (ctx->LGrad || ctx->LSec) {
                                            /* next beginnt mit 0 */
                j = nexta;                  /* nexta mit 1        */

                jk = (j * (j - 1)) / 2 - next;
                xj = 1.0;   

                while (1) {         /* alpha j */
                    if (ctx->LGrad) {
                        tmp = (lltf - del * ezf1) * glam;
                        if (ts > 0.0)
                            tmp -= (llts - del * ezs1) * glam;
                        if (!cen)
                            tmp += del / (1.0 + ezf) - 1.0;
                        tmp *= wt;
                        ctx->Grad[j] += tmp * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;
     
                        while (1) {     /* alpha k */

                            tmp1 = -glam * (lltf - del * ezf1 +
                                    del * ezf1 *
                                   (del / (1.0 + ezf) - 1.0));

                            if (ts > 0.0)
                                tmp1 += glam * (llts - del * ezs1 +
                                        del * ezs1 *
                                       (del / (1.0 + ezs) - 1.0));

                            if (!cen)
                                tmp1 -= del * del * ezf2;

                            tmp1 *= wt;

                            if (k < j)  
                                ctx->Hess[++jk] += tmp1 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp1 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
                nextb = j;
                xj = 1.0;   

                while (1) {                 /* beta j */
                    if (ctx->LGrad) {
                        tmp = -zf * ezf1;
                        if (ts > 0.0)
                            tmp += zs * ezs1;
                        tmp *= glam;
                        if (!cen)
                            tmp += zf1 + 1.0;
                        tmp *= wt;
                        ctx->Grad[j] += tmp * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (1) {         /* alpha k */

                            tmp1 = ezf1 * (zf - del - del * zf1);
                            if (ts > 0.0)
                                tmp1 -= ezs1 * (zs - del - del * zs1);
                            tmp1 *= glam;
                            if (!cen)
                                tmp1 += del * (1.0 - zf * ezf1) / (1.0 + ezf);
                            tmp1 *= wt;

                            ctx->Hess[++jk] += tmp1 * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextb;
                        xk = 1.0;

                        while (1) {         /* beta k */

                            tmp2 = -zf * ezf1 * (zf1 + 1.0);
                            if (ts > 0.0)
                                tmp2 += zs * ezs1 * (zs1 + 1.0);
                            tmp2 *= glam;
                            if (!cen)
                                tmp2 += zf1 * (1.0 - zf * ezf1);
                            tmp2 *= wt;
    
                            if (k < j)  
                                ctx->Hess[++jk] += tmp2 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp2 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
                nextc = j;
                xj = 1.0;   

                while (1) {             /* gamma j */
                    if (ctx->LGrad) {
                        tmp = -lltf;
                        if (ts > 0.0)
                            tmp += llts;
                        tmp *= glam;
                        if (!cen)
                            tmp += 1.0;
                        tmp *= wt;
                        ctx->Grad[j] += tmp * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (1) {         /* alpha k */

                            tmp1 = (lltf - del * ezf1);
                            if (ts > 0.0)
                                tmp1 -= (llts - del * ezs1);
                            tmp1 *= wt * glam;

                            ctx->Hess[++jk] += tmp1 * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextb;
                        xk = 1.0;

                        while (1) {         /* beta k */

                            tmp2 = -zf * ezf1;
                            if (ts > 0.0)
                                tmp2 += zs * ezs1;
                            tmp2 *= wt * glam;

                            ctx->Hess[++jk] += tmp2 * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextc;
                        xk = 1.0;

                        while (1) {         /* gamma k */

                            tmp3 = -lltf;
                            if (ts > 0.0)
                                tmp3 += llts;
                            tmp3 *= wt * glam;

                            if (k < j)  
                                ctx->Hess[++jk] += tmp3 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp3 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
            }
            if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
                mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,icase);
            ctx->NOCUsed++;
        }
    }
    err = checkov(ctx);
    return(err);
}

/****************************************************************************/
/*  fn_ln(resid)                                                            */
/*      Log-normal model (I and II), optional gamma mixture with model I.   */
/*      Calculation of log-likelihood, gradient, hessian.                   */
/*      Episode splitting possible if no mixture.                           */
/*      If resid != 0 print residuals.                                      */

int fn_ln(TDAContext *ctx, int resid)
{
    register int j,k,jv,jk,kv;
    int err,icase,sn,ip,spl,nspl,next,nexta,nextb,nextc,nextd,org,des,cen;
    double bt,bs,xb,xc,xd,xj,xk,wt,a,b,bb,c = 0.0,cbt = 0.0,cbs = 0.0,d = 0.0,hc,ldh = 0.0,zt,zs = 0.0,vzt,vzs = 0.0;
    double ts,tf,tmp,tmp1,tmp2,tmp3,tmp4,ha,haa,hb,hbb,hab,ui = 0.0,dh = 0.0,dh1 = 0.0;

    wt = 1.0;   ctx->NOCUsed = 0;
    err = 0;
    get_spell(ctx, 1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
    while (get_spell(ctx, 0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

        if (ctx->WIVar >= 0)                             /* get weights */
            wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;

        ip = ctx->TranPtr[(sn - 1) * ctx->MaxOrg1 + org];
        if (ip < 0)  
            continue;
 
        j = ctx->PIdxPtr[ip];

        while (sn == ctx->SnTran1[ip] && org == ctx->OrgTran1[ip]) {

            if (des == ctx->DesTran1[ip])
                cen = 0;
            else
                cen = 1;
            ip++;

            nexta = j;
            next = j - 1;                                                   

            a = ctx->TPar[j];                                                   
            while ((jv = ctx->PIdx[++j]))                                                 
                a += ctx->TPar[j] * get_data(ctx, jv - 1,icase);
            xb = ctx->TPar[j];                                                   
            while ((jv = ctx->PIdx[++j]))                                                 
                xb += ctx->TPar[j] * get_data(ctx, jv - 1,icase);

            if (ctx->XMFlg[3]) {     /* if three parameter log-normal model */
                xc = ctx->TPar[j];                                                   
                while ((jv = ctx->PIdx[++j]))                                                 
                    xc += ctx->TPar[j] * get_data(ctx, jv - 1,icase);

                c = rexp(ctx, xc);
                if ((tf -= c) <= 0.0)
                    tf = 0.01;
                if (ts > 0.0 && (ts -= c) < 0.0)
                    ts = 0.0;
            }
            b = rexp(ctx, xb);
            bb = b * b;

            bt = b * tf;
            if (ctx->XMFlg[3]) {
                cbt = c / bt;
                if (ts > 0.0) {
                    bs = b * ts;
                    cbs = c / bs;
                }
            }
            zt = (rlog(ctx, tf) - a) / b;
            vzt = cdnm(ctx, zt);

            if (ts > 0.0) {
                zs = (rlog(ctx, ts) - a) / b;
                vzs = cdnm(ctx, zs);
            }
            tmp = 1.0 - cdnf(ctx, zt);
            if (ts > 0.0)
                tmp /= (1.0 - cdnf(ctx, zs));
            hc = -rlog(ctx, tmp);

            if (resid) {            /* print residuals */
                tmp1 = vzt / bt;
                prn_resid1(ctx, icase + 1,org,des,ts,tf,tmp1,tmp,hc,wt);
            }

            if (ctx->MixTyp) {

                xd = ctx->TPar[j];                                                   
                while ((jv = ctx->PIdx[++j]))                                                 
                    xd += ctx->TPar[j] * get_data(ctx, jv - 1,icase);

                d = rexp(ctx, xd);
                ui = -1.0 / d;
                if (!cen)
                    ui -= 1.0;

                ldh = rlog(ctx, 1.0 + d * hc);
                dh = d / (1.0 + d * hc);
                dh1 = dh * hc;
            }
            if (ctx->LFunc) {

                if (ctx->MixTyp) 
                    tmp = ui * ldh;
                else  
                    tmp = -hc;
                if (!cen)  
                    tmp += rlog(ctx, vzt / bt);

                tmp *= wt;

                ctx->FTmp += tmp;                                     
            }
            if (ctx->LGrad || ctx->LSec) {
                                            /* next starts with 0 */
                j = nexta;                  /* nexta wiht 1       */

                jk = (j * (j - 1)) / 2 - next;
                xj = 1.0;   

                ha = -vzt;
                hb = -zt * vzt;
                haa = -vzt * (zt - vzt);
                hbb = -zt * vzt * (zt * (zt - vzt) - 1.0);
                hab = -vzt * (zt * (zt - vzt) - 1.0);

                if (ts > 0.0) {
                    ha += vzs;
                    hb += zs * vzs;
                    haa += vzs * (zs - vzs);
                    hbb += zs * vzs * (zs * (zs - vzs) - 1.0);
                    hab += vzs * (zs * (zs - vzs) - 1.0);
                }
                ha /= b;
                haa /= bb;
                hab /= b;

                if (ctx->MixTyp) {
                    tmp1 = ui * dh * ha;
                    tmp2 = ui * dh * (haa - dh * ha * ha);
                }
                else {
                    tmp1 = -ha;
                    tmp2 = -haa;
                }
                if (!cen) {
                    tmp1 += (zt - vzt) / b;
                    tmp2 += (vzt * (vzt - zt) - 1.0) / bb;
                }
                tmp1 *= wt;
                tmp2 *= wt;

                while (1) {         /* alpha j */

                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp1 * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp1 * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (1) {     /* alpha k */

                            if (k < j)  
                                ctx->Hess[++jk] += tmp2 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp2 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
                nextb = j;
                xj = 1.0;   

                if (ctx->MixTyp) {
                    tmp1 = ui * dh * hb;
                    tmp2 = ui * dh * (hab - dh * ha * hb);
                    tmp3 = ui * dh * (hbb - dh * hb * hb);
                }
                else {
                    tmp1 = -hb;
                    tmp2 = -hab;
                    tmp3 = -hbb;
                }
                if (!cen) {
                    tmp1 += zt * (zt - vzt) - 1.0;
                    tmp2 += ((1.0 + zt * vzt) * (vzt - zt) - zt) / b;
                    tmp3 += (zt * vzt * (1.0 + zt * (vzt - zt)) -
                                                     2.0 * zt * zt);
                }
                tmp1 *= wt;
                tmp2 *= wt;
                tmp3 *= wt;

                while (1) {     /* beta j */

                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp1 * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp1 * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (1) {         /* alpha k */

                            ctx->Hess[++jk] += tmp2 * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextb;
                        xk = 1.0;

                        while (1) {         /* beta k */
                            if (k < j)  
                                ctx->Hess[++jk] += tmp3 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp3 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
                nextc = j;
                xj = 1.0;   

                if (ctx->XMFlg[3]) {     /* if three-parameter model */

                    tmp1 = cbt * vzt;
                    tmp2 = tmp1 * (zt - vzt) / b;
                    tmp3 = tmp1 * (zt * (zt - vzt) - 1.0);  
                    tmp4 = tmp1 * (tf + c - c * (vzt - zt) / b) / tf;

                    if (ts > 0.0) {
                        tmp = cbs * vzs;
                        tmp1 -= tmp;
                        tmp2 -= tmp * (zs - vzs) / b;
                        tmp3 -= tmp * (zs * (zs - vzs) - 1.0);  
                        tmp4 -= tmp * (ts + c - c * (vzs - zs) / b) / ts;
                    }
                    if (!cen) {
                        tmp1 += cbt * (b + zt - vzt);
                        tmp2 += cbt * (vzt * (vzt - zt) - 1.0) / b;
                        tmp3 += cbt * ((vzt - zt) * (1.0 + zt * vzt) - zt);
                        tmp4 += cbt * ((tf + c) * (b + zt - vzt) / tf -
                                cbt * (1.0 - vzt * (vzt - zt)));
                    }
                    tmp1 *= wt;
                    tmp2 *= wt;
                    tmp3 *= wt;
                    tmp4 *= wt;
                }
                while (ctx->XMFlg[3]) {     /* gamma j */

                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp1 * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp1 * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (1) {         /* alpha k */

                            ctx->Hess[++jk] += tmp2 * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextb;
                        xk = 1.0;

                        while (1) {         /* beta k */

                            ctx->Hess[++jk] += tmp3 * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextc;
                        xk = 1.0;

                        while (1) {         /* gamma k */
                            if (k < j)  
                                ctx->Hess[++jk] += tmp4 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp4 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
                nextd = j;
                xj = 1.0;   

                if (ctx->MixTyp) {
                    tmp1 = wt * (ldh / d + ui * dh1);
                    tmp2 = tmp3 = wt * dh * (1.0 / d + ui * (1.0 - dh1));
                    tmp2 *= ha;
                    tmp3 *= hb;
                    tmp4 = wt * ((2.0 * dh1 - ldh) / d + ui * dh1 * (1.0 - dh1));
                }

                while (ctx->MixTyp) {             /* Delta j */

                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp1 * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp1 * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (1) {         /* Alpha k */

                            ctx->Hess[++jk] += tmp2 * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextb;
                        xk = 1.0;

                        while (1) {         /* Beta k */

                            ctx->Hess[++jk] += tmp3 * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextd;
                        xk = 1.0;

                        while (1) {         /* Delta k */
                            if (k < j)  
                                ctx->Hess[++jk] += tmp4 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp4 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
            }
            if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
                mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,icase);
            ctx->NOCUsed++;
        }
    }
    err = checkov(ctx);
    return(err);
}

/****************************************************************************/
/*  fn_ig(int)                                                              */
/*      Inverse gaussian model.                                             */
/*      Calculation of log-likelihood, gradient, hessian.                   */
/*      Episode splitting possible.                                         */
/*      If resid != 0 print residuals.                                      */

int fn_ig(TDAContext *ctx, int resid)
{
    register int j,k,jv,jk,kv;
    int err,icase,sn,ip,spl,nspl,next,nexta,nextb,org,des,cen;
    double ts,tf,xb,xj,xk,wt,a,b,bb,tmp,tmp1,tmp2,tmp3;
    double t12,s12 = 0.0,z1t,z1tt,z1tt1,z1s = 0.0,z2t,z2s = 0.0,ee,phiz1t,phiz1s = 0.0,phiz2t,phiz2s = 0.0;
    double psiz1t,psiz1s,psiz2t,psiz2s = 0.0,ft,gt,gs = 0.0,gta,gsa = 0.0,gtb,gsb;
    double gtaa,gsaa,gtbb,gsbb,gtab,gsab,res;

    wt = 1.0;   ctx->NOCUsed = 0;
    err = 0;
    get_spell(ctx, 1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
    while (get_spell(ctx, 0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

        if (ctx->WIVar >= 0)                             /* get weights */
            wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;

        ip = ctx->TranPtr[(sn - 1) * ctx->MaxOrg1 + org];
        if (ip < 0)  
            continue;
 
        j = ctx->PIdxPtr[ip];

        while (sn == ctx->SnTran1[ip] && org == ctx->OrgTran1[ip]) {

            if (des == ctx->DesTran1[ip])
                cen = 0;
            else
                cen = 1;
            ip++;

            nexta = j;
            next = j - 1;                                                   

            a = ctx->TPar[j];                                                   
            while ((jv = ctx->PIdx[++j]))                                                 
                a += ctx->TPar[j] * get_data(ctx, jv - 1,icase);
            xb = ctx->TPar[j];                                                   
            while ((jv = ctx->PIdx[++j]))                                                 
                xb += ctx->TPar[j] * get_data(ctx, jv - 1,icase);

            b = rexp(ctx, xb);
            bb = b * b;

            t12 = sqrt(tf);
            z1t = (a * tf - 1.0) / (b * t12);
            z2t = (a * tf + 1.0) / (b * t12);
            phiz1t = dnf(ctx, z1t);
            psiz1t = cdnf(ctx, -z1t);
            psiz2t = cdnf(ctx, -z2t);
            ee  = rexp(ctx, 2.0 * a / bb);

            ft = phiz1t / (b * tf * t12);
            gt = psiz1t - ee * psiz2t;
            if (gt <= 0.0)
                gt = 1.e-12;

            if (ts > 0.0) {
                s12 = sqrt(ts);
                z1s = (a * ts - 1.0) / (b * s12);
                z2s = (a * ts + 1.0) / (b * s12);
                psiz1s = cdnf(ctx, -z1s);
                psiz2s = cdnf(ctx, -z2s);
                gs = psiz1s - ee * psiz2s;
                if (gs <= 0.0)
                    gs = 1.e-12;
            }
            if (resid) {            /* print residuals */
                res = -rlog(ctx, gt);
                tmp2 = gt;
                if (ts > 0.0) {
                    res += rlog(ctx, gs);
                    tmp2 /= gs;
                }
                tmp1 = ft / gt;
                prn_resid1(ctx, icase + 1,org,des,ts,tf,tmp1,tmp2,res,wt);
            }


            if (ctx->LFunc) {

                if (!cen)  
                    tmp = rlog(ctx, ft);
                else
                    tmp = rlog(ctx, gt);

                if (ts > 0.0)   
                    tmp -= rlog(ctx, gs);

                tmp *= wt;

                ctx->FTmp += tmp;                                     
            }
            if (ctx->LGrad || ctx->LSec) {
                                            /* next starts with 0 */
                j = nexta;                  /* nexta wiht 1       */

                jk = (j * (j - 1)) / 2 - next;
                xj = 1.0;   

                z1tt = z1t * z1t;
                z1tt1 = z1tt - 1.0;

                phiz2t = dnf(ctx, z2t);

                if (!cen) {
                    tmp1 = -z1t * t12 / b;  
                    tmp2 = -tf / bb;
                }
                else {
                    gta = (t12 * (ee * phiz2t - phiz1t) -
                                                 ee * 2.0 * psiz2t / b) / b;
                    tmp1 = gta / gt;
                    gtaa = (4.0 * ee * (t12 * phiz2t / b - psiz2t / bb) +
                              tf * (z1t * phiz1t - z2t * phiz2t * ee)) / bb;
                    tmp2 = gtaa / gt - tmp1 * tmp1;
                }
                if (ts > 0.0) {
                    phiz1s = dnf(ctx, z1s);
                    phiz2s = dnf(ctx, z2s);

                    gsa = (s12 * (ee * phiz2s - phiz1s) -
                                                 ee * 2.0 * psiz2s / b) / b;
                    tmp   = gsa / gs;
                    tmp1 -= tmp;
                    gsaa = (4.0 * ee * (s12 * phiz2s / b - psiz2s / bb) +
                              ts * (z1s * phiz1s - z2s * phiz2s * ee)) / bb;
                    tmp2 -= gsaa / gs - tmp * tmp;
                }
                tmp1 *= wt;
                tmp2 *= wt;

                while (1) {         /* Alpha j */

                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp1 * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp1 * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (1) {     /* Alpha k */

                            if (k < j)  
                                ctx->Hess[++jk] += tmp2 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp2 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
                nextb = j;
                xj = 1.0;   

                if (!cen) {
                    tmp1 = z1tt1;
                    tmp2 = 2.0 * z1t * t12 / b; 
                    tmp3 = - 2.0 * z1tt;
                }
                else {
                    gtb = z1t * phiz1t + ee * (4.0 * a * psiz2t / bb - z2t *
                                                               phiz2t);
                    tmp1 = gtb / gt;

                    gtab = t12 * (phiz1t * (1.0 - z1t * z1t) - ee * phiz2t *
                            (1.0 + 4.0 * a / bb - z2t * z2t)) / b +
                            ee * 2.0 * (psiz2t * (4.0 * a / bb + 2.0) -
                            z2t * phiz2t) / bb;

                    tmp2 = gtab / gt - tmp1 * gta / gt; 

                    gtbb = z1t * phiz1t * (z1t * z1t - 1.0) + ee * z2t *
                            phiz2t * (1.0 + 8.0 * a / bb - z2t * z2t) -
                            ee * psiz2t * 8.0 * a * (1.0 + 2.0 * a / bb) / bb;

                    tmp3 = gtbb / gt - tmp1 * tmp1; 
                }
                if (ts > 0.0) {
                    gsb = z1s * phiz1s + ee * (4.0 * a * psiz2s / bb - z2s *
                                                               phiz2s);
                    tmp = gsb / gs;
                    tmp1 -= tmp;

                    gsab = s12 * (phiz1s * (1.0 - z1s * z1s) - ee * phiz2s *
                            (1.0 + 4.0 * a / bb - z2s * z2s)) / b +
                            ee * 2.0 * (psiz2s * (4.0 * a / bb + 2.0) -
                            z2s * phiz2s) / bb;

                    tmp2 -= gsab / gs - tmp * gsa / gs; 

                    gsbb = z1s * phiz1s * (z1s * z1s - 1.0) + ee * z2s *
                            phiz2s * (1.0 + 8.0 * a / bb - z2s * z2s) -
                            ee * psiz2s * 8.0 * a * (1.0 + 2.0 * a / bb) / bb;

                    tmp3 -= gsbb / gs - tmp * tmp; 
                }
                tmp1 *= wt;
                tmp2 *= wt;
                tmp3 *= wt;

                while (1) {     /* Beta j */

                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp1 * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp1 * xj;
                    }
                    if (ctx->LSec) {
                        k = nexta;
                        xk = 1.0;
                        jk += next;

                        while (1) {         /* Alpha k */

                            ctx->Hess[++jk] += tmp2 * xj * xk;

                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                        k = nextb;
                        xk = 1.0;

                        while (1) {         /* Beta k */
                            if (k < j)  
                                ctx->Hess[++jk] += tmp3 * xj * xk;
                            else if (k == j) {
                                ctx->Diag[k] += tmp3 * xj * xk;
                                break;
                            }
                            if (!(kv = ctx->PIdx[++k]))                                                 
                                break;
                            xk = get_data(ctx, kv - 1,icase);
                        }
                    }
                    if (!(jv = ctx->PIdx[++j]))                                                 
                        break;
                    xj = get_data(ctx, jv - 1,icase);
                }
            }
            if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
                mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,icase);
            ctx->NOCUsed++;
        }
    }
    err = checkov(ctx);
    return(err);
}

/****************************************************************************/
/*  fn_gam(resid)                                                           */
/*      Generalized gamma model.                                            */
/*      Calculation of log-likelihood, gradient, hessian                    */
/*      Episode splitting possible.                                         */
/*      If resid != 0 print residuals.                                      */

int fn_gam(TDAContext *ctx, int resid)
{
    register int j,k,jv,jk,kv;
    int err,icase,sn,ip,spl,nspl,next,org,des,cen,fin,fin1;
    double ts,tf,xb,xj,xk,wt,a,b,bb,zs = 0.0,zt,qs,qt,qs1 = 0.0,qt1 = 0.0,k1,k2;
    double lgam,git = 0.0,gis,qsk = 0.0,qtk,qskz = 0.0,qtkz = 0.0,ztk1,tmp1,tmp2,tmp3 = 0.0;

    ctx->NOCUsed = 0;

    k1 = sqrt(ctx->PMKGam);
    k2 = (ctx->PMKGam - 0.5) * rlog(ctx, ctx->PMKGam);
    lgam = loggam(ctx, ctx->PMKGam,&err);           /* Log of gamma function */
    if (err)
        return(err);
         
    wt = 1.0;
    err = 0;
    get_spell(ctx, 1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
    while (get_spell(ctx, 0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

        if (ctx->WIVar >= 0)                             /* get weights */
            wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;

        ip = ctx->TranPtr[(sn - 1) * ctx->MaxOrg1 + org];
        if (ip < 0)  
            continue;
 
        j = ctx->PIdxPtr[ip];

        while (sn == ctx->SnTran1[ip] && org == ctx->OrgTran1[ip]) {

            if (des == ctx->DesTran1[ip])
                cen = 0;
            else
                cen = 1;
            ip++;

            next = j;                                                       
            a = ctx->TPar[j];                                                   
            while ((jv = ctx->PIdx[++j]))                                                 
                a += ctx->TPar[j] * get_data(ctx, jv - 1,icase);
            xb = ctx->TPar[j];                                                   
            while ((jv = ctx->PIdx[++j]))                                                 
                xb += ctx->TPar[j] * get_data(ctx, jv - 1,icase);

            b = rexp(ctx, xb);
            bb = b * b;
 
            zt = (rlog(ctx, tf) - a) / b;
            qt = ctx->PMKGam * rexp(ctx, zt / k1);
            qtk = qt / k1 - k1;
            ztk1 = k1 - (zt / k1 + 1.0) * qt / k1;
            gis = 1.0;

            if (ts > 0.0) {
                zs = (rlog(ctx, ts) - a) / b;
                qs = ctx->PMKGam * rexp(ctx, zs / k1);
                qsk = qs / k1 - k1;
                gis = 1.0 - icgam(ctx, qs,ctx->PMKGam,&err);                  
                if (err)
                    return(err);

                /*****************
                err = icdgam(ctx, qs,PMKGam,dg);
                if (err)
                     printf1(ctx, "ERROR k=%lg qs=%lg\n",PMKGam,qs);
                gis = 1.0 - dg[6];
                *********************/

                qs1 = rexp(ctx, k2 + k1 * zs - qs - lgam) / gis;
                qskz = ((qsk - qs1) * zs - 1.0) * qs1;
            }
            if (cen || resid) {
                git = 1.0 - icgam(ctx, qt,ctx->PMKGam,&err);     
                if (err)
                    return(err);

                /*************************
                err = icdgam(ctx, qt,PMKGam,dg);
                if (err)
                     printf1(ctx, "ERROR k=%lg qs=%lg\n",PMKGam,qt);
                git = 1.0 - dg[6];
                *************************/
            }
            if (cen) {
                qt1 = rexp(ctx, k2 + k1 * zt - qt - lgam) / git;
                qtkz = ((qtk - qt1) * zt - 1.0) * qt1;
            }
            if (resid) {
                tmp1 = git;
                if (ts > 0.0 && gis > 0.0 && gis < 1.0)
                    tmp1 /= gis;
                tmp2 = -rlog(ctx, tmp1);
                tmp3 = rexp(ctx, k2 - xb - rlog(ctx, tf) - lgam + k1 * zt - qt);
                if (git > 0.0)
                    tmp3 /= git;
                else
                    tmp3 = 0.0;

                prn_resid1(ctx, icase + 1,org,des,ts,tf,tmp3,tmp1,tmp2,wt);
            }

            if (ctx->LFunc) {
                if (!cen)  
                    tmp1 = k2 - xb - rlog(ctx, tf) - lgam + k1 * zt - qt;
                else  
                    tmp1 = rlog(ctx, git);
                if (ts > 0.0)
                    tmp1 -= rlog(ctx, gis);
                tmp1 *= wt;
                ctx->FTmp += tmp1;                                     
            }
            if (ctx->LGrad || ctx->LSec) {
                j = next;
                jk = j * (j - 1) / 2;
                if (!cen) {
                    tmp1 = qtk / b;
                    tmp2 = -qt / (ctx->PMKGam * bb);
                }
                else {
                    tmp1 = qt1 / b;
                    tmp2 = (qtk - qt1) * qt1 / bb;
                }
                if (ts > 0.0) {
                    tmp1 -= qs1 / b;
                    tmp2 -= (qsk - qs1) * qs1 / bb;
                }
                tmp1 *= wt;
                tmp2 *= wt;

                jv = 1;
                xj = 1.0;
                fin = 0;
                while (1) {
                    if (ctx->LGrad) {
                        ctx->Grad[j] += tmp1 * xj;
                        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                            ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = tmp1 * xj;
                    }
                    if (ctx->LSec)  
                        ctx->Diag[j] += tmp2 * xj * xj;

                    if (!(jv = ctx->PIdx[++j])) {                                               
                        if (fin)
                            break;
                        if (!cen) {
                            tmp1 = qtk * zt - 1.0;
                            tmp2 = ztk1 * zt;
                            tmp3 = ztk1 / b;
                        }
                        else {
                            tmp1 = qt1 * zt;
                            tmp2 = qtkz * zt;
                            tmp3 = qtkz / b;
                        }
                        if (ts > 0.0) {
                            tmp1 -= qs1 * zs;
                            tmp2 -= qskz * zs;
                            tmp3 -= qskz / b;
                        }
                        tmp1 *= wt;
                        tmp2 *= wt;
                        tmp3 *= wt;

                        xj = 1.0;
                        fin = 1;
                    }
                    else
                        xj = get_data(ctx, jv - 1,icase);

                    if (ctx->LSec) {
                        k = next;                                                       
                        jk += next;                                                   
                        xk = 1.0;
                        fin1 = 0;
                        while (1) {
                            if (!fin || fin1)
                                ctx->Hess[jk] += tmp2 * xj * xk;
                            else
                                ctx->Hess[jk] += tmp3 * xj * xk;
                            k++;                                                          
                            if (k == j)
                                break;
                            if (!(kv = ctx->PIdx[k])) {                                               
                                fin1 = 1;
                                xk = 1.0;
                            }
                            else
                                xk = get_data(ctx, kv - 1,icase);
                            jk++;                                                         
                        }
                    }
                }
            }
            if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
                mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,icase);
            ctx->NOCUsed++;
        }
    }
    err = checkov(ctx);
    return(err);
}






