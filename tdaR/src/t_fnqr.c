/****************************************************************************/
/*  t_fnqr                                                                  */
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
#include "t_gdat.h"
#include "t_qrmod.h"
#include "t_ml.h"
#include "t_gf.h"
#include "t_cdf.h"
#include "t_lin.h"
#include "t_int.h"
#include "t_parm.h"
#include "t_mat.h"
#include "t_matf.h"
#include "tda_context.h"

/*  functions in t_fnqr.c */

int fn_logit1(TDAContext *ctx);
int fn_logit2(TDAContext *ctx);
int fn_logit3(TDAContext *ctx);
int fn_logit4(TDAContext *ctx);
int combo1(TDAContext *ctx, int n, int m, int *iset, int last);
int fn_probit1(TDAContext *ctx);
int fn_probit2(TDAContext *ctx);
int fn_probit3(TDAContext *ctx);
int fn_prob3(TDAContext *ctx, int icase,int wave,int cat,double *par,double *prob);
int fn_probit4(TDAContext *ctx);
int fn_prob4(TDAContext *ctx, double *par,double *t,double *arg,double *prob);
int fn_rmod1(TDAContext *ctx);


/* ------------------------------------------------------------------------ */
/*  fn_logit1                                                               */
/*      Likelihood calculation for a logit model with a binary dependent    */
/*      variable. Multiple waves are pooled.                                */

int fn_logit1(TDAContext *ctx)
{
    register int j,jj,k,jk,wave,icase;
    int cat,yv;
    double wt,xb,exb,exb1,tmp,tmp1,tmp2,xj,xk;

    wt = 1.0; ctx->NOCUsed = 0; 

    for (wave = 0; wave < ctx->NWave; ++wave) {

        yv = (int) ctx->PYVar[wave];      /* dependent variable */

        for (icase = 0; icase < ctx->NOC; ++icase) {
   
            if (!ctx->PDMCFlg || ctx->PDMCIdx[icase * ctx->NWave + wave]) {

                cat = ctx->YCatI[(int)get_data(ctx, yv,icase)];   /* internal category */

                /* cat == 0 is the reference category */

                if (ctx->WIVar >= 0)
                    wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;

                if (ctx->IntFlg)
                    xb = ctx->TPar[1];
                else
                    xb = 0.0;
                for (j = ctx->IntFlg; j < ctx->NParm; ++j) {
                    jj = j - ctx->IntFlg;
                    xb += ctx->TPar[j + 1] * get_data(ctx, ctx->PXVar[jj][wave],icase);
                }
                exb = rexp(ctx, xb);
                exb1 = exb + 1.0;

                if (ctx->LFunc) {
                    tmp = -rlog(ctx, exb1);
                    if (cat)
                        tmp += xb;
                    ctx->FTmp += wt * tmp;                                 
                }
                if (ctx->LGrad || ctx->LSec) {

                    tmp1 = -exb / exb1;
                    if (ctx->CCovFlg != 1)
                        tmp2 = tmp1 / exb1;
                    if (cat)
                        tmp1 += 1.0;
                    if (ctx->CCovFlg == 1)
                        tmp2 = tmp1 * tmp1;
  
                    tmp1 *= wt;
                    tmp2 *= wt;

                    jk = 1;
                    for (j = 0; j < ctx->NParm; ++j) {
                        if (j == 0 && ctx->IntFlg)
                            xj = 1.0;
                        else {
                            jj = j - ctx->IntFlg;
                            xj = get_data(ctx, ctx->PXVar[jj][wave],icase);
                        }
                        if (ctx->LGrad) {
                            ctx->Grad[j + 1] += tmp1 * xj;
                            if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                                ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j + 1] = tmp1 * xj;
                        }

                        if (ctx->LSec) {
                            ctx->Diag[j + 1] += tmp2 * xj * xj;
                            for (k = 0; k < j; ++k) {
                                if (k == 0 && ctx->IntFlg)
                                    xk = 1.0;
                                else {
                                    jj = k - ctx->IntFlg;
                                    xk = get_data(ctx, ctx->PXVar[jj][wave],icase);
                                }
                                ctx->Hess[jk++] += tmp2 * xj * xk;
                            }
                        }
                    }
                }
                if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
                    mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,icase);
                ctx->NOCUsed++;
            }
        }
    }
    ctx->QR_INIT = 0;
    checkov(ctx);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  fn_logit2                                                               */
/*      Likelihood calculation for a logit model with ordered categories    */
/*      (cumulative logits). Multiple waves are pooled.                     */

int fn_logit2(TDAContext *ctx)
{
    register int k,k1,u,v,uv,j,ji,icase,wave,yv;
    double tmp,tmp1,tmp2,tmp3;
    double wt,xu,xv,gtmp,phi,pk,pk1,ck,ck1,dk,dk1;

    wt = 1.0;   ctx->NOCUsed = 0;

    for (wave = 0; wave < ctx->NWave; ++wave) {

        yv = (int) ctx->PYVar[wave];

        for (icase = 0; icase < ctx->NOC; ++icase) {
        
            if (!ctx->PDMCFlg || ctx->PDMCIdx[icase * ctx->NWave + wave]) {

                k1 = ctx->YCatI[(int)get_data(ctx, yv,icase)];   /* internal category */
                k = k1 + 1;

                if (ctx->WIVar >= 0)
                    wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;

                j = 0;
                tmp = 0.0;
                for (u = ctx->NYCat; u <= ctx->NParm; ++u)   
                    tmp += ctx->TPar[u] * get_data(ctx, ctx->PXVar[j++][wave],icase);

                if (k < ctx->NYCat) {
                    tmp1 = rexp(ctx, ctx->TPar[k] + tmp); 
                    tmp2 = 1.0 + tmp1;
                    pk = tmp1 / tmp2;
                    ck = pk / tmp2;
                    dk = ck / tmp2;
                }
                else {
                    pk = ck = dk = 0.0;
                }
                if (k1 >= 1) {
                    tmp1 = rexp(ctx, ctx->TPar[k1] + tmp); 
                    tmp2 = 1.0 + tmp1;
                    pk1 = tmp1 / tmp2;
                    ck1 = pk1 / tmp2;
                    dk1 = ck1 / tmp2;
                }
                else {
                    pk1 = 1.0;
                    ck1 = dk1 = 0.0;
                }
                phi = pk1 - pk;

                if (phi <= 0.0) {
                    phi = ctx->EPSI1;
                    ctx->NPFlgs[6]++;
                }
                if (ctx->LFunc) {
                    ctx->FTmp  += wt * rlog(ctx, phi);  
                }
                if (ctx->LGrad || ctx->LSec) {

                    j = 0;
                    uv = 1;

                    for (u = 1; u <= ctx->NParm; ++u) {

                        if (u < ctx->NYCat) {   /* alpha(u) */

                            /* Gradient alpha(u)) */

                            if (ctx->LGrad) {
                                if (u == k) {
                                    ctx->Grad[u] -= wt * ck / phi;
                                    if (ctx->CCovFlg == 1)
                                        ctx->GTmp[u] = -ck / phi;

                                    if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                                        ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + u] = -wt * ck / phi;
                                }
                                else if (u == k1) {
                                    ctx->Grad[u] += wt * ck1 / phi;
                                    if (ctx->CCovFlg == 1)
                                        ctx->GTmp[u] = ck1 / phi;

                                    if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                                        ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + u] = wt * ck1 / phi;
                                }
                            }
                            if (ctx->LSec) {    /* Hessian alpha(u),alpha(v) */
    
                                for (v = 1; v <= u; ++v) {

                                    if (ctx->CCovFlg == 1) {
                                        tmp = ctx->GTmp[u] * ctx->GTmp[v];
                                    }
                                    else {
                                        tmp1 = tmp2 = tmp3 = 0.0;
                                        if (u == k1) {
                                            tmp2 += ck1;
                                            if (v == k1)
                                                tmp1 += dk1 - ck1 * pk1; 
                                        }
                                        else if (u == k) {
                                            tmp2 -= ck;
                                            if (v == k)
                                                tmp1 -= dk - ck * pk; 
                                        }
                                        if (v == k1)
                                            tmp3 += ck1;
                                        else if (v == k)
                                            tmp3 -= ck;

                                        tmp = (tmp1 - tmp2 * tmp3 / phi) / phi;
                                    }
                                    if (v == u)
                                        ctx->Diag[v] += wt * tmp;
                                    else
                                        ctx->Hess[uv++] += wt * tmp;
                                }
                            }
                        }
                        else {  /* beta */

                            xu = get_data(ctx, ctx->PXVar[j++][wave],icase);
                            gtmp = (ck1 - ck) / phi;

                            if (ctx->LGrad) {        /* Gradient beta */

                                ctx->Grad[u] += wt * xu * gtmp;
                                if (ctx->CCovFlg == 1)
                                    ctx->GTmp[u] = xu * gtmp;

                                if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                                    ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + u] = wt * xu * gtmp;

                            }
                            if (ctx->LSec) {   /* Hessian beta(ctx, u),alpha(v) */

                                ji = 0;
                                for (v = 1; v <= u; ++v) {

                                    if (ctx->CCovFlg == 1) {
                                        tmp = ctx->GTmp[u] * ctx->GTmp[v];
                                    }
                                    else if (v < ctx->NYCat) {

                                              /* Hessian beta(u),alpha(v) */
    
                                        tmp1 = tmp2 = 0.0;
                                        if (v == k1) {
                                            tmp1 += dk1 - ck1 * pk1;
                                            tmp2 += ck1;
                                        }
                                        else if (v == k) {
                                            tmp1 -= dk - ck * pk;
                                            tmp2 -= ck;
                                        }
                                        tmp = xu * (tmp1 - tmp2 * gtmp) / phi; 
                                    }
                                    else { 
                                        /* Hessian beta(u),beta(v) */
  
                                        xv = get_data(ctx, ctx->PXVar[ji++][wave],icase);
                                        tmp1 = dk1 - dk - ck1 * pk1 + ck * pk;
                                        tmp = xu * xv * (tmp1 / phi - gtmp * gtmp);
                                    }
                                    if (v == u)
                                        ctx->Diag[v] += tmp * wt;
                                    else
                                        ctx->Hess[uv++] += tmp * wt;
                                }
                            }
                        }
                    }
                }
                if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
                    mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,icase);
                ctx->NOCUsed++;
            }
        }
    }
    ctx->QR_INIT = 0;
    checkov(ctx);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  fn_logit3                                                               */
/*      Likelihood calculation for a multinomial logit model.               */
/*      Multiple waves are pooled.                                          */

int fn_logit3(TDAContext *ctx)
{
    register int i,j,u,v,uv,ju,jv,iu,iv,zu,zv;
    int k,k1,k2,l,wave,icase,cat,yv;
    double d,wt,tmp,tmp1,tmp2,tmp3,xu,xv;

    wt = 1.0; ctx->NOCUsed = 0;

    for (wave = 0; wave < ctx->NWave; ++wave) {

        yv = (int) ctx->PYVar[wave];

        for (icase = 0; icase < ctx->NOC; ++icase) {

            if (!ctx->PDMCFlg || ctx->PDMCIdx[icase * ctx->NWave + wave]) {

                cat = ctx->YCatI[(int)get_data(ctx, yv,icase)];   /* internal category */

                /* cat == 0 is the reference category */
    
                cat++;
    
                if (ctx->WIVar >= 0)
                    wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;
    
                if (ctx->QR_INIT && ctx->DSVFlg == 0) {
                    for (j = 1; j <= ctx->NYCat; ++j)  
                        ctx->QRProb[j] = 1.0 / (double)ctx->NYCat;
                }
                else {
                    d = 1.0;        /* d is the inclusive value */
                    u = 1;
                    for (j = 2; j <= ctx->NYCat; ++j) {   
    
                        if (ctx->IntFlg)
                            tmp = ctx->TPar[u++];
                        else
                            tmp = 0.0;
    
                        for (i = 0; i < ctx->NPX; ++i) {
                            k = ctx->PXVar[i][wave];
                            tmp += ctx->TPar[u++] * get_data(ctx, k,icase);
                        }
                        l = (ctx->NYCat - 1) * ctx->NPX1 + 1;
    
                        for (i = 0; i < ctx->NPZ; ++i) {
                            k1 = ctx->PZVar[i][wave * ctx->NYCat + j - 1];
                            k2 = ctx->PZVar[i][wave * ctx->NYCat];
                            tmp += ctx->TPar[l++] * (get_data(ctx, k1,icase) -
                                                get_data(ctx, k2,icase));
                        }
   
                        if (tmp < ExpMin || tmp > ExpMax)
                            ctx->NPFlgs[6]++;
      
                        ctx->QRProb[j] = rexp(ctx, tmp);
                        d += ctx->QRProb[j];
                    }
                    ctx->QRProb[1] = 1.0;
                    for (i = 2; i <= ctx->NYCat; ++i) {   
                        ctx->QRProb[i] /= d;
                        ctx->QRProb[1] -= ctx->QRProb[i];
                    }
                }
                if (ctx->LFunc)  
                    ctx->FTmp += wt * rlog(ctx, ctx->QRProb[cat]);                        
    
                if (ctx->LGrad || ctx->LSec) {
    
                    uv = u = 1;
    
                    if (ctx->NPX1) {
    
                        for (ju = 2; ju <= ctx->NYCat; ++ju) {
    
                            tmp = -ctx->QRProb[ju];
                            if (ju == cat)
                                tmp += 1.0;
                        
                            tmp1 = ctx->QRProb[ju] * (1.0 - ctx->QRProb[ju]);
    
                            for (iu = 1; iu <= ctx->NPX1; ++iu) {
    
                                if (iu == 1 && ctx->IntFlg)  
                                    xu = 1.0;
                                else {
                                    k = ctx->PXVar[iu - ctx->IntFlg - 1][wave];
                                    xu = get_data(ctx, k,icase);
                                }
                                if (ctx->LGrad) {
                                    ctx->Grad[u] += wt * tmp * xu;
                                    if (ctx->CCovFlg == 1)
                                        ctx->GTmp[u] = tmp * xu;

                                    if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                                        ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + u] = wt * tmp * xu;

                                }
                                if (ctx->LSec) {
    
                                    tmp2 = tmp1 * xu;
                                    v = 1;
    
                                    for (jv = 2; jv <= ctx->NYCat; ++jv) {
                            
                                        tmp3 = ctx->QRProb[ju] * ctx->QRProb[jv] * xu;
    
                                        for (iv = 1; iv <= ctx->NPX1; ++iv) {
    
                                            if (iv == 1 && ctx->IntFlg)
                                                xv = 1.0;
                                            else {
                                                k = ctx->PXVar[iv - ctx->IntFlg - 1][wave];
                                                xv = get_data(ctx, k,icase);
                                            }
                                            if (jv == ju) {
                                                if (v == u) {
                                                    if (ctx->CCovFlg == 1)
                                                        ctx->Diag[u] += wt *
                                                                ctx->GTmp[u] * ctx->GTmp[v];
                                                    else
                                                        ctx->Diag[u] -= wt * tmp2 * xv;
                                                    goto L3C;
                                                }
                                                else {
                                                    if (ctx->CCovFlg == 1)
                                                        ctx->Hess[uv++] += wt *
                                                                ctx->GTmp[u] * ctx->GTmp[v];
    
                                                    else
                                                        ctx->Hess[uv++] -= wt * tmp2 * xv;
                                                }
                                            }
                                            else {
                                                if (ctx->CCovFlg == 1)
                                                    ctx->Hess[uv++] += wt *
                                                                  ctx->GTmp[u] * ctx->GTmp[v];
                                                else
                                                    ctx->Hess[uv++] += wt * tmp3 * xv;
                                            }
                                            v++;
                                        }
                                    }
                                }
L3C:                            u++;
                            }
                        }
                    }
  
                    if (ctx->NPZ) {      /* generic variables */ 
    
                        for (zu = 0; zu < ctx->NPZ; ++zu) {
    
                            tmp = 0.0;
                            for (ju = 2; ju <= ctx->NYCat; ++ju) {
                                tmp1 = -ctx->QRProb[ju];
                                if (ju == cat)
                                    tmp1 += 1.0;
    
                                k1 = ctx->PZVar[zu][wave * ctx->NYCat + ju - 1];
                                k2 = ctx->PZVar[zu][wave * ctx->NYCat];
                                tmp += tmp1 * (get_data(ctx, k1,icase) -
                                               get_data(ctx, k2,icase));
                            }
                            if (ctx->LGrad) {
                                ctx->Grad[u] += tmp * wt;
                                if (ctx->CCovFlg == 1)
                                    ctx->GTmp[u] = tmp;

                                if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                                    ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + u] = wt * tmp;
                            }
                            if (ctx->LSec) {
    
                                v = 1;
    
                                if (ctx->NPX1) {
                                    for (ju = 2; ju <= ctx->NYCat; ++ju) {
                                        for (iu = 1; iu <= ctx->NPX1; ++iu) {
                                            tmp = 0.0;
                                            for (jv = 2; jv <= ctx->NYCat; ++jv) {
                                                k1 = ctx->PZVar[zu][wave * ctx->NYCat + jv - 1];
                                                k2 = ctx->PZVar[zu][wave * ctx->NYCat];
                                                tmp2 = (get_data(ctx, k1,icase) -
                                                        get_data(ctx, k2,icase));
                                                if (jv == ju) 
                                                    tmp += tmp2 *
                                                       ctx->QRProb[jv] * (1.0 - ctx->QRProb[jv]);
                                                else
                                                    tmp -= tmp2 *
                                                           ctx->QRProb[jv] * ctx->QRProb[ju];
                                            }
                                            if (iu == 1 && ctx->IntFlg)  
                                                xu = 1.0;
                                            else {
                                                k = ctx->PXVar[iu - ctx->IntFlg - 1][wave];
                                                xu = get_data(ctx, k,icase);
                                            }
                                            if (ctx->CCovFlg == 1)
                                                ctx->Hess[uv++] += wt * ctx->GTmp[u] * ctx->GTmp[v];
                                            else
                                                ctx->Hess[uv++] -= tmp * xu * wt;
                                            v++;
                                        }
                                    }
                                }
                                for (zv = 0; zv < ctx->NPZ; ++zv) {
                                    tmp = 0.0;
                                    for (ju = 2; ju <= ctx->NYCat; ++ju) {
                                        tmp1 = 0.0;
                                        for (jv = 2; jv <= ctx->NYCat; ++jv) {
    
                                            k1 = ctx->PZVar[zu][wave * ctx->NYCat + jv - 1];
                                            k2 = ctx->PZVar[zu][wave * ctx->NYCat];
                                            tmp2 = (get_data(ctx, k1,icase) -
                                                    get_data(ctx, k2,icase));
                                            if (jv == ju) 
                                                tmp1 += tmp2 *
                                                        ctx->QRProb[jv] * (1.0 - ctx->QRProb[jv]);
                                            else
                                                tmp1 -= tmp2 *
                                                        ctx->QRProb[jv] * ctx->QRProb[ju];
                                        }
                                        k1 = ctx->PZVar[zv][wave * ctx->NYCat + ju - 1];
                                        k2 = ctx->PZVar[zv][wave * ctx->NYCat];
                                        tmp2 = (get_data(ctx, k1,icase) -
                                                get_data(ctx, k2,icase));
                                        tmp += tmp1 * tmp2;
                                    }
                                    if (zu == zv) {
                                        if (ctx->CCovFlg == 1)
                                            ctx->Diag[u] += wt * ctx->GTmp[u] * ctx->GTmp[u];
                                        else
                                            ctx->Diag[u] -= wt * tmp;
                                        break;
                                    }
                                    else {
                                        if (ctx->CCovFlg == 1)
                                            ctx->Hess[uv++] += wt * ctx->GTmp[u] * ctx->GTmp[v];
                                        else
                                            ctx->Hess[uv++] -= wt * tmp;
                                    }
                                    v++;
                                }
                            }
                            u++;
                        }
                    }
                }
                if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
                    mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,icase);
                ctx->NOCUsed++;
            }
        }
    }
    ctx->QR_INIT = 0;
    checkov(ctx);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  fn_logit4                                                               */
/*  ##  Conditional logistic regression for binary response variable.       */
/*      WrkD, WrkE, WrkH are used for temporary storage.                    */

int fn_logit4(TDAContext *ctx)
{
    register int j,k,l,jk,wave,icase,jj;
    int cat,yv,ny,last,uwave = 0;
    double wt,xby,xbd,exbd,sxbd,xjd,xkd,tmp,tmp1;
       
    ctx->NOCUsed = 0;

    /*  prepare L4_IPTR to map L4_ISET to original wave numbers */
    /*  if PDMCFlg == 1 this must be done inside the icase loop */

    if (ctx->PDMCFlg == 0) {
        uwave = ctx->NWave;
        for (j = 1; j <= ctx->NWave; ++j)  
            ctx->L4_IPTR[j] = j;
    }

    wt = 1.0;             
    ctx->NFLUsed = 0;        /* count individuals used for likelihood */

    for (icase = 0; icase < ctx->NOC; ++icase) {
   
        if (ctx->PDMCFlg) {
            uwave = 0;
            k = 1;
            jj = icase * ctx->NWave;    
            for (j = 0; j < ctx->NWave; ++j) {
                if (ctx->PDMCIdx[jj++]) {
                    uwave++;
                    ctx->L4_IPTR[k++] = j + 1;
                }
            }
        }
        if (ctx->WIVar >= 0)
            wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;

        ny = 0;
        xby = 0.0;
        for (wave = 0; wave < ctx->NWave; ++wave) {

            if (!ctx->PDMCFlg || ctx->PDMCIdx[icase * ctx->NWave + wave]) {
       
                yv = (int) ctx->PYVar[wave];

                /* get internal category of dependent variable with possible
                   values: 0 and 1. 0 is taken as the reference category. */

                cat = ctx->YCatI[(int)get_data(ctx, yv,icase)];

                tmp = 0.0;
                for (j = 1; j <= ctx->NParm; ++j)  
                    tmp += ctx->TPar[j] * get_data(ctx, ctx->PXVar[j - 1][wave],icase);
                ctx->L4_XB[wave] = tmp;

                if (cat) {
                    ny++;
                    xby += tmp;
                }
            }
        }

        if (ny > 0 && ny < uwave) {

            ctx->NFLUsed++;

            if (ctx->LGrad || ctx->LSec) {
                dclear(ctx, ctx->NParm,ctx->WrkD);
                dclear(ctx, ctx->NParm,ctx->WrkE);
                dclear(ctx, ctx->HSiz,ctx->WrkH);
            }
            sxbd = 0.0;
            last = 1;

            while (!(last = combo1(ctx, uwave,ny,ctx->L4_ISET,last))) {

                /************************ 
                printf1(ctx, "\nny=%d\n",ny);
                for (l = 1; l <= ny; ++l)
                    printf1(ctx, "%2d (%2d) ",L4_ISET[l],L4_IPTR[L4_ISET[l]]);
                *************************/

                xbd = 0.0;
                for (l = 1; l <= ny; ++l)  
                    xbd += ctx->L4_XB[ctx->L4_IPTR[ctx->L4_ISET[l]] - 1];
                exbd = rexp(ctx, xbd);
                sxbd += exbd;

                if (ctx->LGrad || ctx->LSec) {
                    jk = 1;
                    for (j = 1; j <= ctx->NParm; ++j) {
                        xjd = 0.0;
                        for (l = 1; l <= ny; ++l)  
                            xjd += get_data(ctx, ctx->PXVar[j - 1][ctx->L4_IPTR[ctx->L4_ISET[l]]-1],icase);

                        if (ctx->LGrad)
                            ctx->WrkD[j] += exbd * xjd;

                        if (ctx->LSec && ctx->CCovFlg != 1) {

                            ctx->WrkE[j] += exbd * xjd * xjd;
                            for (k = 1; k < j; ++k) {
                                xkd = 0.0;
                                for (l = 1; l <= ny; ++l)  
                                    xkd += get_data(ctx, ctx->PXVar[k - 1][ctx->L4_IPTR[ctx->L4_ISET[l]]-1],icase);
                                ctx->WrkH[jk] += exbd * xjd * xkd;
                                jk++;
                            }
                        }
                    }
                }
            }
            /********************************  
            printf1(ctx, "icase=%d\n",icase);
            printfe(ctx, "icase=%d\n",icase);
            ********************************/
            if (ctx->LFunc)  
                ctx->FTmp += wt * (xby - rlog(ctx, sxbd));                   

            if (ctx->LGrad || ctx->LSec) {

                jk = 1;
                for (j = 1; j <= ctx->NParm; ++j) {

                    if (ctx->LGrad || (ctx->LSec && ctx->CCovFlg == 1)) {

                        tmp = 0.0;
                        for (wave = 0; wave < ctx->NWave; ++wave) {

                            if (!ctx->PDMCFlg || ctx->PDMCIdx[icase * ctx->NWave + wave]) {
                                yv = (int) ctx->PYVar[wave];
                                cat = ctx->YCatI[(int)get_data(ctx, yv,icase)];

                                if (cat)
                                    tmp += get_data(ctx, ctx->PXVar[j - 1][wave],icase);
                            }
                        }
                        tmp1 = tmp - ctx->WrkD[j] / sxbd;
                        if (ctx->LGrad) {
                            ctx->Grad[j] += wt * tmp1;

                            if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                                ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j] = wt * tmp1;
                        }
                        if (ctx->CCovFlg == 1)
                            ctx->WrkE[j] = tmp1;
                    }
                    if (ctx->LSec) {

                        if (ctx->CCovFlg == 1)
                            ctx->Diag[j] += wt * ctx->WrkE[j] * ctx->WrkE[j];
                        else
                            ctx->Diag[j] += wt * (ctx->WrkD[j] * ctx->WrkD[j] / sxbd - ctx->WrkE[j]) / sxbd;

                        for (k = 1; k < j; ++k) {
                            if (ctx->CCovFlg == 1)
                                ctx->Hess[jk] += wt * ctx->WrkE[j] * ctx->WrkE[k];
                            else
                                ctx->Hess[jk] += wt * (ctx->WrkD[j] * ctx->WrkD[k] / sxbd - ctx->WrkH[jk]) / sxbd;
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
    ctx->QR_INIT = 0;
    checkov(ctx);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  combo1. Generate combinations used by fn_logit4().                      */

int combo1(TDAContext *ctx, int n, int m, int *iset, int last)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,k,l; 

    if (last) {
        for (i = 1; i <= m; ++i)
            iset[i] = i;
    }
    else {
        k = m;

        while (k > 0) {
            l = iset[k] + 1;
            if ((l + m - k) <= n) {
                for (i = k; i <= m; ++i)  
                    iset[i] = l++;
                return(last);
            }
            k--;
        }
    }
    if (last)
        last = 0;
    else
        last = 1;
    return(last);
}

/* ------------------------------------------------------------------------ */
/*  fn_probit1                                                              */
/*      Likelihood calculation for a probit model with a binary dependent   */
/*      variable.                                                           */

int fn_probit1(TDAContext *ctx)
{
    register int j,jj,k,jk,wave,icase;
    int cat,yv;
    double wt,xb,phi,phi1,v,tmp1,tmp2,xj,xk;

    wt = 1.0; ctx->NOCUsed = 0;  

    for (wave = 0; wave < ctx->NWave; ++wave) {

        yv = (int) ctx->PYVar[wave];

        for (icase = 0; icase < ctx->NOC; ++icase) {
   
            if (!ctx->PDMCFlg || ctx->PDMCIdx[icase * ctx->NWave + wave]) {

                cat = ctx->YCatI[(int)get_data(ctx, yv,icase)];   /* internal category */
    
                /* cat == 0 is the reference category */
    
                if (ctx->WIVar >= 0)
                    wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;
    
                if (ctx->IntFlg)
                    xb = ctx->TPar[1];
                else
                    xb = 0.0;
                for (j = ctx->IntFlg; j < ctx->NParm; ++j) {
                    jj = j - ctx->IntFlg;
                    xb += ctx->TPar[j + 1] * get_data(ctx, ctx->PXVar[jj][wave],icase);
                }
                phi = cdnf(ctx, xb);
                if (cat)
                    phi1 = phi;
                else
                    phi1 = 1.0 - phi;
    
                if (ctx->LFunc)  
                    ctx->FTmp += wt * rlog(ctx, phi1);                          
    
                if (ctx->LGrad || ctx->LSec) {
    
                    v = dnf(ctx, xb) / phi1;
                    if (cat) {
                        tmp1 = v;
                        tmp2 = -v * (xb + v);
                    }
                    else {
                        tmp1 = -v;
                        tmp2 = v * (xb - v);
                    }
                    if (ctx->CCovFlg == 1)               /* outer product of Grad */
                        tmp2 = wt * tmp1 * tmp1;
                    else
                        tmp2 *= wt;
                    tmp1 *= wt;
    
                    jk = 1;
                    for (j = 0; j < ctx->NParm; ++j) {
                        if (j == 0 && ctx->IntFlg)
                            xj = 1.0;
                        else {
                            jj = j - ctx->IntFlg;
                            xj = get_data(ctx, ctx->PXVar[jj][wave],icase);
                        }
                        if (ctx->LGrad) {
                            ctx->Grad[j + 1] += tmp1 * xj;

                            if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                                ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + j + 1] = tmp1 * xj;
                        }
                        if (ctx->LSec) {
                            ctx->Diag[j + 1] += tmp2 * xj * xj;
                            for (k = 0; k < j; ++k) {
                                if (k == 0 && ctx->IntFlg)
                                    xk = 1.0;
                                else {
                                    jj = k - ctx->IntFlg;
                                    xk = get_data(ctx, ctx->PXVar[jj][wave],icase);
                                }
                                ctx->Hess[jk++] += tmp2 * xj * xk;
                            }
                        }
                    }
                }
                if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
                    mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,icase);
                ctx->NOCUsed++;
            }
        }
    }
    ctx->QR_INIT = 0;
    checkov(ctx);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  fn_probit2                                                              */
/*      Likelihood calculation for a probit model with ordered categories   */
/*      Multiple waves are pooled.                                          */

int fn_probit2(TDAContext *ctx)
{
    register int k,k1,u,v,uv,j,ji,icase,wave,yv;
    double tmp,tmp1,tmp2,tmp3;
    double wt,xu,xv,gtmp,phi,pk,pk1,ck,ck1,dk,dk1;

    wt = 1.0; ctx->NOCUsed = 0;

    for (wave = 0; wave < ctx->NWave; ++wave) {

        yv = (int) ctx->PYVar[wave];

        for (icase = 0; icase < ctx->NOC; ++icase) {
   
            if (!ctx->PDMCFlg || ctx->PDMCIdx[icase * ctx->NWave + wave]) {
 
                k1 = ctx->YCatI[(int)get_data(ctx, yv,icase)];   /* internal category */
                k = k1 + 1;
    
                if (ctx->WIVar >= 0)
                    wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;
    
                j = 0;
                tmp = 0.0;
                for (u = ctx->NYCat; u <= ctx->NParm; ++u)   
                    tmp += ctx->TPar[u] * get_data(ctx, ctx->PXVar[j++][wave],icase);
    
                if (k < ctx->NYCat) {
                    ck = tmp + ctx->TPar[k];
                    pk = cdnf(ctx, ck);
                    dk = dnf(ctx, ck);
                }
                else {
                    pk = ck = dk = 0.0;
                }
                if (k1 >= 1) {
                    ck1 = tmp + ctx->TPar[k1];
                    pk1 = cdnf(ctx, ck1);
                    dk1 = dnf(ctx, ck1);
                }
                else {
                    pk1 = 1.0;
                    ck1 = dk1 = 0.0;
                }
                phi = pk1 - pk;
    
                if (phi <= 0.0) {
                    phi = ctx->EPSI1;
                    ctx->NPFlgs[6]++;
                }
                if (ctx->LFunc) {
                    ctx->FTmp  += wt * rlog(ctx, phi);  
                }
                if (ctx->LGrad || ctx->LSec) {
    
                    j = 0;
                    uv = 1;
    
                    for (u = 1; u <= ctx->NParm; ++u) {
    
                        if (u < ctx->NYCat) {   /* alpha(u) */
    
                            /* Gradient alpha(u)) */
    
                            if (ctx->LGrad) {
                                if (u == k) {
                                    ctx->Grad[u] -= wt * dk / phi;
                                    if (ctx->CCovFlg == 1)
                                        ctx->GTmp[u] = -dk / phi;

                                    if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                                        ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + u] = -wt * dk / phi;
                                }
                                else if (u == k1) {
                                    ctx->Grad[u] += wt * dk1 / phi;
                                    if (ctx->CCovFlg == 1)
                                        ctx->GTmp[u] = dk1 / phi;

                                    if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                                        ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + u] = wt * dk1 / phi;
                                }
                            }
                            if (ctx->LSec) {    /* Hessian alpha(u),alpha(v) */
    
                                for (v = 1; v <= u; ++v) {
    
                                    if (ctx->CCovFlg == 1) {
                                        tmp = ctx->GTmp[u] * ctx->GTmp[v];
                                    }
                                    else {
                                        tmp1 = tmp2 = tmp3 = 0.0;
                                        if (u == k1) {
                                            tmp2 += dk1;
                                            if (v == k1)
                                                tmp1 -= dk1 * ck1; 
                                        }
                                        else if (u == k) {
                                            tmp2 -= dk;
                                            if (v == k)
                                                tmp1 += dk * ck; 
                                        }
                                        if (v == k1)
                                            tmp3 += dk1;
                                        else if (v == k)
                                            tmp3 -= dk;
    
                                        tmp = (tmp1 - tmp2 * tmp3 / phi) / phi;
                                    }
                                    if (v == u)
                                        ctx->Diag[v] += wt * tmp;
                                    else
                                        ctx->Hess[uv++] += wt * tmp;
                                }
                            }
                        }
                        else {  /* beta */
    
                            xu = get_data(ctx, ctx->PXVar[j++][wave],icase);
                            gtmp = (dk1 - dk) / phi;
    
                            if (ctx->LGrad) {        /* Gradient beta */
    
                                ctx->Grad[u] += xu * wt * gtmp;
                                if (ctx->CCovFlg == 1)
                                    ctx->GTmp[u] = xu * gtmp;

                                if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                                    ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + u] = wt * xu * gtmp;
                            }
                            if (ctx->LSec) {   /* Hessian beta(ctx, u),alpha(v) */
    
                                ji = 0;
                                for (v = 1; v <= u; ++v) {
    
                                    if (ctx->CCovFlg == 1) {
                                        tmp = ctx->GTmp[u] * ctx->GTmp[v];
                                    }
                                    else if (v < ctx->NYCat) {
                                                    
                                             /* Hessian beta(u),alpha(v) */
    
                                        tmp1 = tmp2 = 0.0;
                                        if (v == k1) {
                                            tmp1 -= dk1 * ck1;
                                            tmp2 += dk1;
                                        }
                                        else if (v == k) {
                                            tmp1 += dk * ck;
                                            tmp2 -= dk;
                                        }
                                        tmp = xu * (tmp1 - tmp2 * gtmp) / phi; 
                                    }
                                    else {
                                             /* Hessian beta(u),beta(v) */
    
                                        xv = get_data(ctx, ctx->PXVar[ji++][wave],icase);
                                        tmp1 = dk * ck - dk1 * ck1;
                                        tmp = xu * xv * (tmp1 / phi - gtmp * gtmp);
                                    }
                                    if (v == u)
                                        ctx->Diag[v] += tmp * wt;
                                    else
                                        ctx->Hess[uv++] += tmp * wt;
                                }
                            }
                        }
                    }
                }
                if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
                    mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,icase);
                ctx->NOCUsed++;
            }
        }
    }
    ctx->QR_INIT = 0;
    checkov(ctx);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  fn_probit3                                                              */
/*      Likelihood calculation for a multivariate probit model.             */
/*      Multiple waves are pooled.                                          */
/*                                                                          */
/*      Parameterization depends on PMOPT (1,2,3).                          */
/*                                                                          */  

int fn_probit3(TDAContext *ctx)
{
    int wave,icase,err;                     
    double wt,prob;                  

    wt = 1.0;   ctx->NOCUsed = 0;

    for (wave = 0; wave < ctx->NWave; ++wave) {

        for (icase = 0; icase < ctx->NOC; ++icase) {

            if (!ctx->PDMCFlg || ctx->PDMCIdx[icase * ctx->NWave + wave]) {

                /* fn_prob3's degenerate-covariance returns leave prob
                   unwritten while this loop continues past them (only
                   err < 0 aborts); seed it so the rlog clamp turns
                   those steps into a deterministic penalty instead of
                   reading an uninitialised value (valgrind, s. 41) */
                prob = 0.0;
                err = fn_prob3(ctx, icase,wave,-1,ctx->TPar,&prob);
                if (err) {
                    if (err < 0) {
                        ctx->NPFlgs[9] += 1;
                        return(1);
                    }
                    ctx->NPFlgs[10] += 1;
                }
                if (ctx->WIVar >= 0)
                    wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;

                if (ctx->LFunc)   
                    ctx->FTmp += wt * rlog(ctx, prob);                        

                ctx->NOCUsed++;
            }
        }
    }
    ctx->QR_INIT = 0;
    checkov(ctx);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  fn_prob3    Get probability for icase and wave. Return error.           */
/*              If cat >= 0 probability for this category, otherwise        */
/*              for observed category.                                      */
/*                                                                          */
/*  Parameterization depends on PMOPT.                                      */
/*  PMOPT = 1: no transformation                                            */
/*  PMOPT = 2: sigma = exp(par) / (1 + exp(par))                            */
/*  PMOPT = 3: assumes that par are given as Cholesky factors               */

int fn_prob3(TDAContext *ctx, int icase,int wave,int cat,double *par,double *prob)
{
    register int i,j,k,l,kk;
    int cat1,yv,jptr,ny1,m,k1,l1;
    double xi,s,t,tmp,tmp1,tmp2;
    double ss[DMVMax + 1];

    ny1 = ctx->NYCat - 1;
    if (cat < 0) {
        yv = (int) ctx->PYVar[wave];
        cat = ctx->YCatI[(int)get_data(ctx, yv,icase)];   /* internal category */
    }
    cat1 = cat + 1;

    for (l = 1; l <= ny1; ++l)
        ctx->DMVArg[l] = 0.0;

    jptr = (cat - 1) * ctx->NPX1 + 1;

    if (ctx->IntFlg)
        i = -1;
    else
        i = 0;

    kk = 1;
    for ( ; i < ctx->NPX; ++i) {
        if (i < 0)
            xi = 1.0;
        else {
            k = ctx->PXVar[i][wave];
            xi = get_data(ctx, k,icase);
        }
        tmp2 = 0.0;
        if (cat > 0)
            tmp2 = par[jptr];

        l = 1;
        for (k = 0; k < ctx->NYCat; ++k) {
            if (k != cat) {
                tmp1 = 0.0;
                if (k > 0)
                    tmp1 = par[(k - 1) * ctx->NPX1 + kk];

                ctx->DMVArg[l] += xi * (tmp2 - tmp1);
                l++;
            }
        }
        jptr++;
        kk++;
    }
    if (ctx->NPZ > 0) {

        j = ny1 * ctx->NPX1 + 1;

        for (i = 0; i < ctx->NPZ; ++i) {

            l = 1;
            for (k = 0; k < ctx->NYCat; ++k) {

                if (k != cat) {
                    tmp1 = get_data(ctx, ctx->PZVar[i][wave * ctx->NYCat + cat],icase);
                    tmp2 = get_data(ctx, ctx->PZVar[i][wave * ctx->NYCat + k],icase);

                    ctx->DMVArg[l] += (tmp1 - tmp2) * par[j];
                    l++;
                }
            }
            j++;
        }
    }

    /* set up covariance matrix in QRProb */

    jptr = ny1 * ctx->NPX1 + ctx->NPZ;    /* points to sigma in par */

    for (i = 1; i <= ny1; ++i) {
        for (j = 1; j <= ny1; ++j) {
            tmp = 0.0;
            for (k = 1; k <= ctx->NYCat; ++k) {
                tmp1 = 0.0;
                for (l = 1; l <= ctx->NYCat; ++l) {

                    if (ctx->PMOPT == 3) {
                        par[0] = 1.0;
                        if (k == 1)  
                            k1 = 0;
                        else  
                            k1 = jptr + ((k - 1) * (k - 2) / 2) + 1;
                        if (l == 1)  
                            l1 = 0;
                        else  
                            l1 = jptr + ((l - 1) * (l - 2) / 2) + 1;

                        m = imin(ctx, k,l);
                        s = 0.0;
                        for (kk = 1; kk < m; ++kk)  
                            s += par[l1++] * par[k1++];
                        if (k > m)
                            s += par[k1];
                        else if (l > m)
                            s += par[l1];
                        else
                            s += 1.0;

                    }
                    else {
                        if (k == l)
                            s = 1.0;
                        else {
                            if (l < k)      
                                s = par[jptr + ((k - 1) * (k - 2)) / 2 + l];
                            else               
                                s = par[jptr + ((l - 1) * (l - 2)) / 2 + k];

                            if (ctx->PMOPT == 2) {
                                s = rexp(ctx, s);
                                s = 2.0 * s / (1.0 + s) - 1.0;
                            }
                        }
                    }

                    if (l == cat1)
                        t = -1.0;
                    else if (l < cat1 && l == j)
                        t = 1.0;
                    else if (l > cat1 && l == j + 1)
                        t = 1.0;
                    else
                        t = 0.0;

                    tmp1 += s * t;
                }
                if (k == cat1)
                    t = -1.0;
                else if (k < cat1 && k == i)
                    t = 1.0;
                else if (k > cat1 && k == i + 1)
                    t = 1.0;
                else
                    t = 0.0;

                tmp += tmp1 * t;     
            }
            ctx->QRProb[(i - 1) * ny1 + j] = tmp;
        }   
    }
    /***
    printf1(ctx, "QRProb\n");
    for (i = 1; i <= ny1; ++i) {
        for (j = 1; j <= ny1; ++j)
            printf1(ctx, "%lg ",QRProb[(i - 1) * ny1 + j]);
        newline(ctx);
    }
    ***/

    /* change into correlation matrix */

    for (i = 1; i <= ny1; ++i) {
        tmp = ctx->QRProb[(i - 1) * ny1 + i];
        if (tmp <= 0.0) {
            ctx->NPFlgs[9] += 1;
            return(1);
        }
        ss[i] = sqrt(tmp);
    }
    for (i = 1; i <= ny1; ++i) {
        for (j = 1; j <= ny1; ++j)  
            ctx->QRProb[(i - 1) * ny1 + j] /= (ss[i] * ss[j]);
    }
    /***
    printf1(ctx, "CORR \n");
    for (i = 1; i <= ny1; ++i) {
        for (j = 1; j <= ny1; ++j)
            printf1(ctx, "%lg ",QRProb[(i - 1) * ny1 + j]);
        newline(ctx);
    }
    ***/

    k = dmv(ctx, ny1,ctx->PMNHP,ctx->DMVArg,ctx->QRProb,prob,ctx->PMEPS,&tmp1);
    return(k);
}

/* ------------------------------------------------------------------------ */
/*  fn_probit4                                                              */
/*      Likelihood calculation for simultaneous binary probit model.        */
/*                                                                          */
/*      Parameterization depends on PMOPT (1,2,3).                          */
/*                                                                          */  

int fn_probit4(TDAContext *ctx)
{
    register int i,j;
    int wave,icase,err,cat;                 
    double wt,xi,tmp;           

    wt = 1.0;   ctx->NOCUsed = 0;

    for (icase = 0; icase < ctx->NOC; ++icase) {

        for (wave = 0; wave < ctx->NWave; ++wave) {

            ctx->QRTmp[wave + 1] = 1.0;

            if (!ctx->PDMCFlg || ctx->PDMCIdx[icase * ctx->NWave + wave]) {

                tmp = 0.0;

                if (ctx->IntFlg)
                    i = -1;
                else
                    i = 0;

                j = wave * ctx->NPX1;
                for ( ; i < ctx->NPX; ++i) {
                    if (i < 0)
                        xi = 1.0;
                    else  
                        xi = get_data(ctx, ctx->PXVar[i][wave],icase);
                    tmp += xi * ctx->TPar[++j];
                }
                cat = ctx->YCatI[(int)get_data(ctx, ctx->PYVar[wave],icase)];   /* internal category */
                if (cat == 0) {
                    tmp *= -1.0;
                    ctx->QRTmp[wave + 1] = -1.0;
                }
                ctx->NOCUsed++;
            }
            else
                tmp = 12.0;         /* some high value */

            ctx->DMVArg[wave + 1] = tmp;
        }
        /****
        printf1(ctx, "DMVArg: "); 
        for (i = 1; i <= NWave; ++i)
            printf1(ctx, "%lg ",DMVArg[i]);
        printf1(ctx, "\nQRTmp: "); 
        for (i = 1; i <= NWave; ++i)
            printf1(ctx, "%lg ",QRTmp[i]);
        newline(ctx);
        ***/

        err = fn_prob4(ctx, ctx->TPar,ctx->QRTmp,ctx->DMVArg,&tmp);

        if (err) {
            if (err < 0) {
                ctx->NPFlgs[9] += 1;
                return(1);
            }
            ctx->NPFlgs[10] += 1;
        }
        if (ctx->WIVar >= 0)
            wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;

        if (ctx->LFunc)   
            ctx->FTmp += wt * rlog(ctx, tmp);                        
    }
    ctx->QR_INIT = 0;
    checkov(ctx);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  fn_prob4    Get probability for icase and wave. Return error.           */
/*              If cat >= 0 probability for this category, otherwise        */
/*              for observed category.                                      */
/*                                                                          */
/*  Parameterization depends on PMOPT.                                      */
/*  PMOPT = 1: no transformation                                            */
/*  PMOPT = 2: sigma = exp(par) / (1 + exp(par))                            */
/*  PMOPT = 3: assumes that par are given as Cholesky factors               */

int fn_prob4(TDAContext *ctx, double *par,double *t,double *arg,double *prob)
{
    register int k,l,kk;
    int jptr,m,k1,l1;
    double s,tmp;
    double ss[DMVMax + 1];

    /* set up covariance matrix in QRProb */

    jptr = ctx->NWave * ctx->NPX1;    /* points to sigma in par */

    for (k = 1; k <= ctx->NWave; ++k) {
        for (l = 1; l <= ctx->NWave; ++l) {

            if (ctx->PMOPT == 3) {
                par[0] = 1.0;
                if (k == 1)  
                    k1 = 0;
                else  
                    k1 = jptr + ((k - 1) * (k - 2) / 2) + 1;
                if (l == 1)  
                    l1 = 0;
                else  
                    l1 = jptr + ((l - 1) * (l - 2) / 2) + 1;

                m = imin(ctx, k,l);
                s = 0.0;
                for (kk = 1; kk < m; ++kk)  
                    s += par[l1++] * par[k1++];
                if (k > m)
                    s += par[k1];
                else if (l > m)
                    s += par[l1];
                else
                    s += 1.0;

            }
            else {
                if (k == l)
                    s = 1.0;
                else {
                    if (l < k)      
                        s = par[jptr + ((k - 1) * (k - 2)) / 2 + l];
                    else               
                        s = par[jptr + ((l - 1) * (l - 2)) / 2 + k];

                    if (ctx->PMOPT == 2) {
                        s = rexp(ctx, s);
                        s = 2.0 * s / (1.0 + s) - 1.0;
                    }
                }
            }
            s *= t[k] * t[l];
            ctx->QRProb[(k - 1) * ctx->NWave + l] = s;   
        }   
    }
    /******* 
    printf1(ctx, "QRProb\n");
    for (k = 1; k <= NWave; ++k) {
        for (l = 1; l <= NWave; ++l)
            printf1(ctx, "%lg ",QRProb[(k - 1) * NWave + l]);
        newline(ctx);
    }
    ****/
  
    /* change into correlation matrix */

    for (k = 1; k <= ctx->NWave; ++k) {
        tmp = ctx->QRProb[(k - 1) * ctx->NWave + k];
        if (tmp <= 0.0) {
            ctx->NPFlgs[9] += 1;
            return(1);
        }
        ss[k] = sqrt(tmp);
    }
    for (k = 1; k <= ctx->NWave; ++k) {
        for (l = 1; l <= ctx->NWave; ++l)  
            ctx->QRProb[(k - 1) * ctx->NWave + l] /= (ss[k] * ss[l]);
    }
    /****
    printf1(ctx, "CORR \n");
    for (k = 1; k <= NWave; ++k) {
        for (l = 1; l <= NWave; ++l)
            printf1(ctx, "%lg ",QRProb[(k - 1) * NWave + l]);
        newline(ctx);
    }
    ****/
  
    k = dmv(ctx, ctx->NWave,ctx->PMNHP,arg,ctx->QRProb,prob,ctx->PMEPS,&tmp);
    return(k);
}

/* ------------------------------------------------------------------------ */
/*  fn_rmod1                                                                */
/*      Likelihood calculation for simple Rasch model.                      */

int fn_rmod1(TDAContext *ctx)
{
    register int i,j;
    int xij;
    double tmp,p;

    for (i = 1; i <= ctx->RMPN; ++i) {
        tmp = 0.0;
        for (j = 1; j <= ctx->RMPV; ++j) {
            xij = ctx->RMPAT[(i - 1) * ctx->RMPV + j];

            /* the second term is the pattern (person) parameter, so it
               is indexed by the pattern i; indexing it by the item j
               left the person parameters out of the likelihood
               entirely and confounded the item parameters pairwise,
               which made the hessian singular for every data set.
               The parameter count RMPV + RMPN allocated for this
               command matches the i-indexed reading. */
            p = ctx->TPar[j] + ctx->TPar[ctx->RMPV + i];
            tmp += (double)xij * p - rlog(ctx, 1.0 + rexp(ctx, p));       

            /* the minimization requests the gradient (LGrad); the
               dispatcher clears Grad and the model function fills it.
               Without this the cleared gradient reads as zero and the
               BFGS algorithm stops after zero iterations.  The score
               of the Rasch log likelihood is (x_ij - p_ij), summed
               into the item and the pattern parameter. */
            if (ctx->LGrad) {
                double pij = rexp(ctx, p) / (1.0 + rexp(ctx, p));
                double sc = ((double)xij - pij) * ctx->RMPF[i];
                ctx->Grad[j] += sc;
                ctx->Grad[ctx->RMPV + i] += sc;
            }
        }
        if (ctx->LFunc) {
            ctx->FTmp += tmp * ctx->RMPF[i];                       
        }
    }
    checkov(ctx);
    return(0);
}
