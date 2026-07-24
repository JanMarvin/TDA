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

/*  functions in t_fnqr.c */

int fn_logit1(void);
int fn_logit2(void);
int fn_logit3(void);
int fn_logit4(void);
int combo1(int n, int m, int *iset, int last);
int fn_probit1(void);
int fn_probit2(void);
int fn_probit3(void);
int fn_prob3(int icase,int wave,int cat,double *par,double *prob);
int fn_probit4(void);
int fn_prob4(double *par,double *t,double *arg,double *prob);
int fn_rmod1(void);

int NFLUsed = 0;    /* number of cases used in fn_logit4() */

/* ------------------------------------------------------------------------ */
/*  fn_logit1                                                               */
/*      Likelihood calculation for a logit model with a binary dependent    */
/*      variable. Multiple waves are pooled.                                */

int fn_logit1(void)
{
    register int j,jj,k,jk,wave,icase;
    int cat,yv;
    double wt,xb,exb,exb1,tmp,tmp1,tmp2,xj,xk;

    wt = 1.0; NOCUsed = 0; 

    for (wave = 0; wave < NWave; ++wave) {

        yv = (int) PYVar[wave];      /* dependent variable */

        for (icase = 0; icase < NOC; ++icase) {
   
            if (!PDMCFlg || PDMCIdx[icase * NWave + wave]) {

                cat = YCatI[(int)get_data(yv,icase)];   /* internal category */

                /* cat == 0 is the reference category */

                if (WIVar >= 0)
                    wt = get_data(WIVar,icase) * WNorm;

                if (IntFlg)
                    xb = TPar[1];
                else
                    xb = 0.0;
                for (j = IntFlg; j < NParm; ++j) {
                    jj = j - IntFlg;
                    xb += TPar[j + 1] * get_data(PXVar[jj][wave],icase);
                }
                exb = rexp(xb);
                exb1 = exb + 1.0;

                if (LFunc) {
                    tmp = -rlog(exb1);
                    if (cat)
                        tmp += xb;
                    FTmp += wt * tmp;                                 
                }
                if (LGrad || LSec) {

                    tmp1 = -exb / exb1;
                    if (CCovFlg != 1)
                        tmp2 = tmp1 / exb1;
                    if (cat)
                        tmp1 += 1.0;
                    if (CCovFlg == 1)
                        tmp2 = tmp1 * tmp1;
  
                    tmp1 *= wt;
                    tmp2 *= wt;

                    jk = 1;
                    for (j = 0; j < NParm; ++j) {
                        if (j == 0 && IntFlg)
                            xj = 1.0;
                        else {
                            jj = j - IntFlg;
                            xj = get_data(PXVar[jj][wave],icase);
                        }
                        if (LGrad) {
                            Grad[j + 1] += tmp1 * xj;
                            if (CGradFlg && NOCUsed < MPGradRow)
                                MatVal[MPGradIdx][NOCUsed * MPGradCol + j + 1] = tmp1 * xj;
                        }

                        if (LSec) {
                            Diag[j + 1] += tmp2 * xj * xj;
                            for (k = 0; k < j; ++k) {
                                if (k == 0 && IntFlg)
                                    xk = 1.0;
                                else {
                                    jj = k - IntFlg;
                                    xk = get_data(PXVar[jj][wave],icase);
                                }
                                Hess[jk++] += tmp2 * xj * xk;
                            }
                        }
                    }
                }
                if (CGradFlg && NOCUsed < MPGradRow && PM2NV > 0)
                    mp_putvar(NOCUsed,MPGradCol,MatVal[MPGradIdx],PM2NV,PM2VIdx,icase);
                NOCUsed++;
            }
        }
    }
    QR_INIT = 0;
    checkov();
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  fn_logit2                                                               */
/*      Likelihood calculation for a logit model with ordered categories    */
/*      (cumulative logits). Multiple waves are pooled.                     */

int fn_logit2(void)
{
    register int k,k1,u,v,uv,j,ji,icase,wave,yv;
    double tmp,tmp1,tmp2,tmp3;
    double wt,xu,xv,gtmp,phi,pk,pk1,ck,ck1,dk,dk1;

    wt = 1.0;   NOCUsed = 0;

    for (wave = 0; wave < NWave; ++wave) {

        yv = (int) PYVar[wave];

        for (icase = 0; icase < NOC; ++icase) {
        
            if (!PDMCFlg || PDMCIdx[icase * NWave + wave]) {

                k1 = YCatI[(int)get_data(yv,icase)];   /* internal category */
                k = k1 + 1;

                if (WIVar >= 0)
                    wt = get_data(WIVar,icase) * WNorm;

                j = 0;
                tmp = 0.0;
                for (u = NYCat; u <= NParm; ++u)   
                    tmp += TPar[u] * get_data(PXVar[j++][wave],icase);

                if (k < NYCat) {
                    tmp1 = rexp(TPar[k] + tmp); 
                    tmp2 = 1.0 + tmp1;
                    pk = tmp1 / tmp2;
                    ck = pk / tmp2;
                    dk = ck / tmp2;
                }
                else {
                    pk = ck = dk = 0.0;
                }
                if (k1 >= 1) {
                    tmp1 = rexp(TPar[k1] + tmp); 
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
                    phi = EPSI1;
                    NPFlgs[6]++;
                }
                if (LFunc) {
                    FTmp  += wt * rlog(phi);  
                }
                if (LGrad || LSec) {

                    j = 0;
                    uv = 1;

                    for (u = 1; u <= NParm; ++u) {

                        if (u < NYCat) {   /* alpha(u) */

                            /* Gradient alpha(u)) */

                            if (LGrad) {
                                if (u == k) {
                                    Grad[u] -= wt * ck / phi;
                                    if (CCovFlg == 1)
                                        GTmp[u] = -ck / phi;

                                    if (CGradFlg && NOCUsed < MPGradRow)
                                        MatVal[MPGradIdx][NOCUsed * MPGradCol + u] = -wt * ck / phi;
                                }
                                else if (u == k1) {
                                    Grad[u] += wt * ck1 / phi;
                                    if (CCovFlg == 1)
                                        GTmp[u] = ck1 / phi;

                                    if (CGradFlg && NOCUsed < MPGradRow)
                                        MatVal[MPGradIdx][NOCUsed * MPGradCol + u] = wt * ck1 / phi;
                                }
                            }
                            if (LSec) {    /* Hessian alpha(u),alpha(v) */
    
                                for (v = 1; v <= u; ++v) {

                                    if (CCovFlg == 1) {
                                        tmp = GTmp[u] * GTmp[v];
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
                                        Diag[v] += wt * tmp;
                                    else
                                        Hess[uv++] += wt * tmp;
                                }
                            }
                        }
                        else {  /* beta */

                            xu = get_data(PXVar[j++][wave],icase);
                            gtmp = (ck1 - ck) / phi;

                            if (LGrad) {        /* Gradient beta */

                                Grad[u] += wt * xu * gtmp;
                                if (CCovFlg == 1)
                                    GTmp[u] = xu * gtmp;

                                if (CGradFlg && NOCUsed < MPGradRow)
                                    MatVal[MPGradIdx][NOCUsed * MPGradCol + u] = wt * xu * gtmp;

                            }
                            if (LSec) {   /* Hessian beta(u),alpha(v) */

                                ji = 0;
                                for (v = 1; v <= u; ++v) {

                                    if (CCovFlg == 1) {
                                        tmp = GTmp[u] * GTmp[v];
                                    }
                                    else if (v < NYCat) {

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
  
                                        xv = get_data(PXVar[ji++][wave],icase);
                                        tmp1 = dk1 - dk - ck1 * pk1 + ck * pk;
                                        tmp = xu * xv * (tmp1 / phi - gtmp * gtmp);
                                    }
                                    if (v == u)
                                        Diag[v] += tmp * wt;
                                    else
                                        Hess[uv++] += tmp * wt;
                                }
                            }
                        }
                    }
                }
                if (CGradFlg && NOCUsed < MPGradRow && PM2NV > 0)
                    mp_putvar(NOCUsed,MPGradCol,MatVal[MPGradIdx],PM2NV,PM2VIdx,icase);
                NOCUsed++;
            }
        }
    }
    QR_INIT = 0;
    checkov();
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  fn_logit3                                                               */
/*      Likelihood calculation for a multinomial logit model.               */
/*      Multiple waves are pooled.                                          */

int fn_logit3(void)
{
    register int i,j,u,v,uv,ju,jv,iu,iv,zu,zv;
    int k,k1,k2,l,wave,icase,cat,yv;
    double d,wt,tmp,tmp1,tmp2,tmp3,xu,xv;

    wt = 1.0; NOCUsed = 0;

    for (wave = 0; wave < NWave; ++wave) {

        yv = (int) PYVar[wave];

        for (icase = 0; icase < NOC; ++icase) {

            if (!PDMCFlg || PDMCIdx[icase * NWave + wave]) {

                cat = YCatI[(int)get_data(yv,icase)];   /* internal category */

                /* cat == 0 is the reference category */
    
                cat++;
    
                if (WIVar >= 0)
                    wt = get_data(WIVar,icase) * WNorm;
    
                if (QR_INIT && DSVFlg == 0) {
                    for (j = 1; j <= NYCat; ++j)  
                        QRProb[j] = 1.0 / (double)NYCat;
                }
                else {
                    d = 1.0;        /* d is the inclusive value */
                    u = 1;
                    for (j = 2; j <= NYCat; ++j) {   
    
                        if (IntFlg)
                            tmp = TPar[u++];
                        else
                            tmp = 0.0;
    
                        for (i = 0; i < NPX; ++i) {
                            k = PXVar[i][wave];
                            tmp += TPar[u++] * get_data(k,icase);
                        }
                        l = (NYCat - 1) * NPX1 + 1;
    
                        for (i = 0; i < NPZ; ++i) {
                            k1 = PZVar[i][wave * NYCat + j - 1];
                            k2 = PZVar[i][wave * NYCat];
                            tmp += TPar[l++] * (get_data(k1,icase) -
                                                get_data(k2,icase));
                        }
   
                        if (tmp < ExpMin || tmp > ExpMax)
                            NPFlgs[6]++;
      
                        QRProb[j] = rexp(tmp);
                        d += QRProb[j];
                    }
                    QRProb[1] = 1.0;
                    for (i = 2; i <= NYCat; ++i) {   
                        QRProb[i] /= d;
                        QRProb[1] -= QRProb[i];
                    }
                }
                if (LFunc)  
                    FTmp += wt * rlog(QRProb[cat]);                        
    
                if (LGrad || LSec) {
    
                    uv = u = 1;
    
                    if (NPX1) {
    
                        for (ju = 2; ju <= NYCat; ++ju) {
    
                            tmp = -QRProb[ju];
                            if (ju == cat)
                                tmp += 1.0;
                        
                            tmp1 = QRProb[ju] * (1.0 - QRProb[ju]);
    
                            for (iu = 1; iu <= NPX1; ++iu) {
    
                                if (iu == 1 && IntFlg)  
                                    xu = 1.0;
                                else {
                                    k = PXVar[iu - IntFlg - 1][wave];
                                    xu = get_data(k,icase);
                                }
                                if (LGrad) {
                                    Grad[u] += wt * tmp * xu;
                                    if (CCovFlg == 1)
                                        GTmp[u] = tmp * xu;

                                    if (CGradFlg && NOCUsed < MPGradRow)
                                        MatVal[MPGradIdx][NOCUsed * MPGradCol + u] = wt * tmp * xu;

                                }
                                if (LSec) {
    
                                    tmp2 = tmp1 * xu;
                                    v = 1;
    
                                    for (jv = 2; jv <= NYCat; ++jv) {
                            
                                        tmp3 = QRProb[ju] * QRProb[jv] * xu;
    
                                        for (iv = 1; iv <= NPX1; ++iv) {
    
                                            if (iv == 1 && IntFlg)
                                                xv = 1.0;
                                            else {
                                                k = PXVar[iv - IntFlg - 1][wave];
                                                xv = get_data(k,icase);
                                            }
                                            if (jv == ju) {
                                                if (v == u) {
                                                    if (CCovFlg == 1)
                                                        Diag[u] += wt *
                                                                GTmp[u] * GTmp[v];
                                                    else
                                                        Diag[u] -= wt * tmp2 * xv;
                                                    goto L3C;
                                                }
                                                else {
                                                    if (CCovFlg == 1)
                                                        Hess[uv++] += wt *
                                                                GTmp[u] * GTmp[v];
    
                                                    else
                                                        Hess[uv++] -= wt * tmp2 * xv;
                                                }
                                            }
                                            else {
                                                if (CCovFlg == 1)
                                                    Hess[uv++] += wt *
                                                                  GTmp[u] * GTmp[v];
                                                else
                                                    Hess[uv++] += wt * tmp3 * xv;
                                            }
                                            v++;
                                        }
                                    }
                                }
L3C:                            u++;
                            }
                        }
                    }
  
                    if (NPZ) {      /* generic variables */ 
    
                        for (zu = 0; zu < NPZ; ++zu) {
    
                            tmp = 0.0;
                            for (ju = 2; ju <= NYCat; ++ju) {
                                tmp1 = -QRProb[ju];
                                if (ju == cat)
                                    tmp1 += 1.0;
    
                                k1 = PZVar[zu][wave * NYCat + ju - 1];
                                k2 = PZVar[zu][wave * NYCat];
                                tmp += tmp1 * (get_data(k1,icase) -
                                               get_data(k2,icase));
                            }
                            if (LGrad) {
                                Grad[u] += tmp * wt;
                                if (CCovFlg == 1)
                                    GTmp[u] = tmp;

                                if (CGradFlg && NOCUsed < MPGradRow)
                                    MatVal[MPGradIdx][NOCUsed * MPGradCol + u] = wt * tmp;
                            }
                            if (LSec) {
    
                                v = 1;
    
                                if (NPX1) {
                                    for (ju = 2; ju <= NYCat; ++ju) {
                                        for (iu = 1; iu <= NPX1; ++iu) {
                                            tmp = 0.0;
                                            for (jv = 2; jv <= NYCat; ++jv) {
                                                k1 = PZVar[zu][wave * NYCat + jv - 1];
                                                k2 = PZVar[zu][wave * NYCat];
                                                tmp2 = (get_data(k1,icase) -
                                                        get_data(k2,icase));
                                                if (jv == ju) 
                                                    tmp += tmp2 *
                                                       QRProb[jv] * (1.0 - QRProb[jv]);
                                                else
                                                    tmp -= tmp2 *
                                                           QRProb[jv] * QRProb[ju];
                                            }
                                            if (iu == 1 && IntFlg)  
                                                xu = 1.0;
                                            else {
                                                k = PXVar[iu - IntFlg - 1][wave];
                                                xu = get_data(k,icase);
                                            }
                                            if (CCovFlg == 1)
                                                Hess[uv++] += wt * GTmp[u] * GTmp[v];
                                            else
                                                Hess[uv++] -= tmp * xu * wt;
                                            v++;
                                        }
                                    }
                                }
                                for (zv = 0; zv < NPZ; ++zv) {
                                    tmp = 0.0;
                                    for (ju = 2; ju <= NYCat; ++ju) {
                                        tmp1 = 0.0;
                                        for (jv = 2; jv <= NYCat; ++jv) {
    
                                            k1 = PZVar[zu][wave * NYCat + jv - 1];
                                            k2 = PZVar[zu][wave * NYCat];
                                            tmp2 = (get_data(k1,icase) -
                                                    get_data(k2,icase));
                                            if (jv == ju) 
                                                tmp1 += tmp2 *
                                                        QRProb[jv] * (1.0 - QRProb[jv]);
                                            else
                                                tmp1 -= tmp2 *
                                                        QRProb[jv] * QRProb[ju];
                                        }
                                        k1 = PZVar[zv][wave * NYCat + ju - 1];
                                        k2 = PZVar[zv][wave * NYCat];
                                        tmp2 = (get_data(k1,icase) -
                                                get_data(k2,icase));
                                        tmp += tmp1 * tmp2;
                                    }
                                    if (zu == zv) {
                                        if (CCovFlg == 1)
                                            Diag[u] += wt * GTmp[u] * GTmp[u];
                                        else
                                            Diag[u] -= wt * tmp;
                                        break;
                                    }
                                    else {
                                        if (CCovFlg == 1)
                                            Hess[uv++] += wt * GTmp[u] * GTmp[v];
                                        else
                                            Hess[uv++] -= wt * tmp;
                                    }
                                    v++;
                                }
                            }
                            u++;
                        }
                    }
                }
                if (CGradFlg && NOCUsed < MPGradRow && PM2NV > 0)
                    mp_putvar(NOCUsed,MPGradCol,MatVal[MPGradIdx],PM2NV,PM2VIdx,icase);
                NOCUsed++;
            }
        }
    }
    QR_INIT = 0;
    checkov();
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  fn_logit4                                                               */
/*  ##  Conditional logistic regression for binary response variable.       */
/*      WrkD, WrkE, WrkH are used for temporary storage.                    */

int fn_logit4(void)
{
    register int j,k,l,jk,wave,icase,jj;
    int cat,yv,ny,last,uwave;
    double wt,xby,xbd,exbd,sxbd,xjd,xkd,tmp,tmp1;
       
    NOCUsed = 0;

    /*  prepare L4_IPTR to map L4_ISET to original wave numbers */
    /*  if PDMCFlg == 1 this must be done inside the icase loop */

    if (PDMCFlg == 0) {
        uwave = NWave;
        for (j = 1; j <= NWave; ++j)  
            L4_IPTR[j] = j;
    }

    wt = 1.0;             
    NFLUsed = 0;        /* count individuals used for likelihood */

    for (icase = 0; icase < NOC; ++icase) {
   
        if (PDMCFlg) {
            uwave = 0;
            k = 1;
            jj = icase * NWave;    
            for (j = 0; j < NWave; ++j) {
                if (PDMCIdx[jj++]) {
                    uwave++;
                    L4_IPTR[k++] = j + 1;
                }
            }
        }
        if (WIVar >= 0)
            wt = get_data(WIVar,icase) * WNorm;

        ny = 0;
        xby = 0.0;
        for (wave = 0; wave < NWave; ++wave) {

            if (!PDMCFlg || PDMCIdx[icase * NWave + wave]) {
       
                yv = (int) PYVar[wave];

                /* get internal category of dependent variable with possible
                   values: 0 and 1. 0 is taken as the reference category. */

                cat = YCatI[(int)get_data(yv,icase)];

                tmp = 0.0;
                for (j = 1; j <= NParm; ++j)  
                    tmp += TPar[j] * get_data(PXVar[j - 1][wave],icase);
                L4_XB[wave] = tmp;

                if (cat) {
                    ny++;
                    xby += tmp;
                }
            }
        }

        if (ny > 0 && ny < uwave) {

            NFLUsed++;

            if (LGrad || LSec) {
                dclear(NParm,WrkD);
                dclear(NParm,WrkE);
                dclear(HSiz,WrkH);
            }
            sxbd = 0.0;
            last = 1;

            while (!(last = combo1(uwave,ny,L4_ISET,last))) {

                /************************ 
                printf1("\nny=%d\n",ny);
                for (l = 1; l <= ny; ++l)
                    printf1("%2d (%2d) ",L4_ISET[l],L4_IPTR[L4_ISET[l]]);
                *************************/

                xbd = 0.0;
                for (l = 1; l <= ny; ++l)  
                    xbd += L4_XB[L4_IPTR[L4_ISET[l]] - 1];
                exbd = rexp(xbd);
                sxbd += exbd;

                if (LGrad || LSec) {
                    jk = 1;
                    for (j = 1; j <= NParm; ++j) {
                        xjd = 0.0;
                        for (l = 1; l <= ny; ++l)  
                            xjd += get_data(PXVar[j - 1][L4_IPTR[L4_ISET[l]]-1],icase);

                        if (LGrad)
                            WrkD[j] += exbd * xjd;

                        if (LSec && CCovFlg != 1) {

                            WrkE[j] += exbd * xjd * xjd;
                            for (k = 1; k < j; ++k) {
                                xkd = 0.0;
                                for (l = 1; l <= ny; ++l)  
                                    xkd += get_data(PXVar[k - 1][L4_IPTR[L4_ISET[l]]-1],icase);
                                WrkH[jk] += exbd * xjd * xkd;
                                jk++;
                            }
                        }
                    }
                }
            }
            /********************************  
            printf1("icase=%d\n",icase);
            printfe("icase=%d\n",icase);
            ********************************/
            if (LFunc)  
                FTmp += wt * (xby - rlog(sxbd));                   

            if (LGrad || LSec) {

                jk = 1;
                for (j = 1; j <= NParm; ++j) {

                    if (LGrad || (LSec && CCovFlg == 1)) {

                        tmp = 0.0;
                        for (wave = 0; wave < NWave; ++wave) {

                            if (!PDMCFlg || PDMCIdx[icase * NWave + wave]) {
                                yv = (int) PYVar[wave];
                                cat = YCatI[(int)get_data(yv,icase)];

                                if (cat)
                                    tmp += get_data(PXVar[j - 1][wave],icase);
                            }
                        }
                        tmp1 = tmp - WrkD[j] / sxbd;
                        if (LGrad) {
                            Grad[j] += wt * tmp1;

                            if (CGradFlg && NOCUsed < MPGradRow)
                                MatVal[MPGradIdx][NOCUsed * MPGradCol + j] = wt * tmp1;
                        }
                        if (CCovFlg == 1)
                            WrkE[j] = tmp1;
                    }
                    if (LSec) {

                        if (CCovFlg == 1)
                            Diag[j] += wt * WrkE[j] * WrkE[j];
                        else
                            Diag[j] += wt * (WrkD[j] * WrkD[j] / sxbd - WrkE[j]) / sxbd;

                        for (k = 1; k < j; ++k) {
                            if (CCovFlg == 1)
                                Hess[jk] += wt * WrkE[j] * WrkE[k];
                            else
                                Hess[jk] += wt * (WrkD[j] * WrkD[k] / sxbd - WrkH[jk]) / sxbd;
                            jk++;
                        }
                    }
                }
            }
            if (CGradFlg && NOCUsed < MPGradRow && PM2NV > 0)
                mp_putvar(NOCUsed,MPGradCol,MatVal[MPGradIdx],PM2NV,PM2VIdx,icase);
            NOCUsed++; 
        }
    }
    QR_INIT = 0;
    checkov();
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  combo1. Generate combinations used by fn_logit4().                      */

int combo1(int n, int m, int *iset, int last)
{
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

int fn_probit1(void)
{
    register int j,jj,k,jk,wave,icase;
    int cat,yv;
    double wt,xb,phi,phi1,v,tmp1,tmp2,xj,xk;

    wt = 1.0; NOCUsed = 0;  

    for (wave = 0; wave < NWave; ++wave) {

        yv = (int) PYVar[wave];

        for (icase = 0; icase < NOC; ++icase) {
   
            if (!PDMCFlg || PDMCIdx[icase * NWave + wave]) {

                cat = YCatI[(int)get_data(yv,icase)];   /* internal category */
    
                /* cat == 0 is the reference category */
    
                if (WIVar >= 0)
                    wt = get_data(WIVar,icase) * WNorm;
    
                if (IntFlg)
                    xb = TPar[1];
                else
                    xb = 0.0;
                for (j = IntFlg; j < NParm; ++j) {
                    jj = j - IntFlg;
                    xb += TPar[j + 1] * get_data(PXVar[jj][wave],icase);
                }
                phi = cdnf(xb);
                if (cat)
                    phi1 = phi;
                else
                    phi1 = 1.0 - phi;
    
                if (LFunc)  
                    FTmp += wt * rlog(phi1);                          
    
                if (LGrad || LSec) {
    
                    v = dnf(xb) / phi1;
                    if (cat) {
                        tmp1 = v;
                        tmp2 = -v * (xb + v);
                    }
                    else {
                        tmp1 = -v;
                        tmp2 = v * (xb - v);
                    }
                    if (CCovFlg == 1)               /* outer product of Grad */
                        tmp2 = wt * tmp1 * tmp1;
                    else
                        tmp2 *= wt;
                    tmp1 *= wt;
    
                    jk = 1;
                    for (j = 0; j < NParm; ++j) {
                        if (j == 0 && IntFlg)
                            xj = 1.0;
                        else {
                            jj = j - IntFlg;
                            xj = get_data(PXVar[jj][wave],icase);
                        }
                        if (LGrad) {
                            Grad[j + 1] += tmp1 * xj;

                            if (CGradFlg && NOCUsed < MPGradRow)
                                MatVal[MPGradIdx][NOCUsed * MPGradCol + j + 1] = tmp1 * xj;
                        }
                        if (LSec) {
                            Diag[j + 1] += tmp2 * xj * xj;
                            for (k = 0; k < j; ++k) {
                                if (k == 0 && IntFlg)
                                    xk = 1.0;
                                else {
                                    jj = k - IntFlg;
                                    xk = get_data(PXVar[jj][wave],icase);
                                }
                                Hess[jk++] += tmp2 * xj * xk;
                            }
                        }
                    }
                }
                if (CGradFlg && NOCUsed < MPGradRow && PM2NV > 0)
                    mp_putvar(NOCUsed,MPGradCol,MatVal[MPGradIdx],PM2NV,PM2VIdx,icase);
                NOCUsed++;
            }
        }
    }
    QR_INIT = 0;
    checkov();
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  fn_probit2                                                              */
/*      Likelihood calculation for a probit model with ordered categories   */
/*      Multiple waves are pooled.                                          */

int fn_probit2(void)
{
    register int k,k1,u,v,uv,j,ji,icase,wave,yv;
    double tmp,tmp1,tmp2,tmp3;
    double wt,xu,xv,gtmp,phi,pk,pk1,ck,ck1,dk,dk1;

    wt = 1.0; NOCUsed = 0;

    for (wave = 0; wave < NWave; ++wave) {

        yv = (int) PYVar[wave];

        for (icase = 0; icase < NOC; ++icase) {
   
            if (!PDMCFlg || PDMCIdx[icase * NWave + wave]) {
 
                k1 = YCatI[(int)get_data(yv,icase)];   /* internal category */
                k = k1 + 1;
    
                if (WIVar >= 0)
                    wt = get_data(WIVar,icase) * WNorm;
    
                j = 0;
                tmp = 0.0;
                for (u = NYCat; u <= NParm; ++u)   
                    tmp += TPar[u] * get_data(PXVar[j++][wave],icase);
    
                if (k < NYCat) {
                    ck = tmp + TPar[k];
                    pk = cdnf(ck);
                    dk = dnf(ck);
                }
                else {
                    pk = ck = dk = 0.0;
                }
                if (k1 >= 1) {
                    ck1 = tmp + TPar[k1];
                    pk1 = cdnf(ck1);
                    dk1 = dnf(ck1);
                }
                else {
                    pk1 = 1.0;
                    ck1 = dk1 = 0.0;
                }
                phi = pk1 - pk;
    
                if (phi <= 0.0) {
                    phi = EPSI1;
                    NPFlgs[6]++;
                }
                if (LFunc) {
                    FTmp  += wt * rlog(phi);  
                }
                if (LGrad || LSec) {
    
                    j = 0;
                    uv = 1;
    
                    for (u = 1; u <= NParm; ++u) {
    
                        if (u < NYCat) {   /* alpha(u) */
    
                            /* Gradient alpha(u)) */
    
                            if (LGrad) {
                                if (u == k) {
                                    Grad[u] -= wt * dk / phi;
                                    if (CCovFlg == 1)
                                        GTmp[u] = -dk / phi;

                                    if (CGradFlg && NOCUsed < MPGradRow)
                                        MatVal[MPGradIdx][NOCUsed * MPGradCol + u] = -wt * dk / phi;
                                }
                                else if (u == k1) {
                                    Grad[u] += wt * dk1 / phi;
                                    if (CCovFlg == 1)
                                        GTmp[u] = dk1 / phi;

                                    if (CGradFlg && NOCUsed < MPGradRow)
                                        MatVal[MPGradIdx][NOCUsed * MPGradCol + u] = wt * dk1 / phi;
                                }
                            }
                            if (LSec) {    /* Hessian alpha(u),alpha(v) */
    
                                for (v = 1; v <= u; ++v) {
    
                                    if (CCovFlg == 1) {
                                        tmp = GTmp[u] * GTmp[v];
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
                                        Diag[v] += wt * tmp;
                                    else
                                        Hess[uv++] += wt * tmp;
                                }
                            }
                        }
                        else {  /* beta */
    
                            xu = get_data(PXVar[j++][wave],icase);
                            gtmp = (dk1 - dk) / phi;
    
                            if (LGrad) {        /* Gradient beta */
    
                                Grad[u] += xu * wt * gtmp;
                                if (CCovFlg == 1)
                                    GTmp[u] = xu * gtmp;

                                if (CGradFlg && NOCUsed < MPGradRow)
                                    MatVal[MPGradIdx][NOCUsed * MPGradCol + u] = wt * xu * gtmp;
                            }
                            if (LSec) {   /* Hessian beta(u),alpha(v) */
    
                                ji = 0;
                                for (v = 1; v <= u; ++v) {
    
                                    if (CCovFlg == 1) {
                                        tmp = GTmp[u] * GTmp[v];
                                    }
                                    else if (v < NYCat) {
                                                    
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
    
                                        xv = get_data(PXVar[ji++][wave],icase);
                                        tmp1 = dk * ck - dk1 * ck1;
                                        tmp = xu * xv * (tmp1 / phi - gtmp * gtmp);
                                    }
                                    if (v == u)
                                        Diag[v] += tmp * wt;
                                    else
                                        Hess[uv++] += tmp * wt;
                                }
                            }
                        }
                    }
                }
                if (CGradFlg && NOCUsed < MPGradRow && PM2NV > 0)
                    mp_putvar(NOCUsed,MPGradCol,MatVal[MPGradIdx],PM2NV,PM2VIdx,icase);
                NOCUsed++;
            }
        }
    }
    QR_INIT = 0;
    checkov();
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  fn_probit3                                                              */
/*      Likelihood calculation for a multivariate probit model.             */
/*      Multiple waves are pooled.                                          */
/*                                                                          */
/*      Parameterization depends on PMOPT (1,2,3).                          */
/*                                                                          */  

int fn_probit3(void)
{
    int wave,icase,err;                     
    double wt,prob;                  

    wt = 1.0;   NOCUsed = 0;

    for (wave = 0; wave < NWave; ++wave) {

        for (icase = 0; icase < NOC; ++icase) {

            if (!PDMCFlg || PDMCIdx[icase * NWave + wave]) {

                err = fn_prob3(icase,wave,-1,TPar,&prob);
                if (err) {
                    if (err < 0) {
                        NPFlgs[9] += 1;
                        return(1);
                    }
                    NPFlgs[10] += 1;
                }
                if (WIVar >= 0)
                    wt = get_data(WIVar,icase) * WNorm;

                if (LFunc)   
                    FTmp += wt * rlog(prob);                        

                NOCUsed++;
            }
        }
    }
    QR_INIT = 0;
    checkov();
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

int fn_prob3(int icase,int wave,int cat,double *par,double *prob)
{
    register int i,j,k,l,kk;
    int cat1,yv,jptr,ny1,m,k1,l1;
    double xi,s,t,tmp,tmp1,tmp2;
    double ss[DMVMax + 1];

    ny1 = NYCat - 1;
    if (cat < 0) {
        yv = (int) PYVar[wave];
        cat = YCatI[(int)get_data(yv,icase)];   /* internal category */
    }
    cat1 = cat + 1;

    for (l = 1; l <= ny1; ++l)
        DMVArg[l] = 0.0;

    jptr = (cat - 1) * NPX1 + 1;

    if (IntFlg)
        i = -1;
    else
        i = 0;

    kk = 1;
    for ( ; i < NPX; ++i) {
        if (i < 0)
            xi = 1.0;
        else {
            k = PXVar[i][wave];
            xi = get_data(k,icase);
        }
        tmp2 = 0.0;
        if (cat > 0)
            tmp2 = par[jptr];

        l = 1;
        for (k = 0; k < NYCat; ++k) {
            if (k != cat) {
                tmp1 = 0.0;
                if (k > 0)
                    tmp1 = par[(k - 1) * NPX1 + kk];

                DMVArg[l] += xi * (tmp2 - tmp1);
                l++;
            }
        }
        jptr++;
        kk++;
    }
    if (NPZ > 0) {

        j = ny1 * NPX1 + 1;

        for (i = 0; i < NPZ; ++i) {

            l = 1;
            for (k = 0; k < NYCat; ++k) {

                if (k != cat) {
                    tmp1 = get_data(PZVar[i][wave * NYCat + cat],icase);
                    tmp2 = get_data(PZVar[i][wave * NYCat + k],icase);

                    DMVArg[l] += (tmp1 - tmp2) * par[j];
                    l++;
                }
            }
            j++;
        }
    }

    /* set up covariance matrix in QRProb */

    jptr = ny1 * NPX1 + NPZ;    /* points to sigma in par */

    for (i = 1; i <= ny1; ++i) {
        for (j = 1; j <= ny1; ++j) {
            tmp = 0.0;
            for (k = 1; k <= NYCat; ++k) {
                tmp1 = 0.0;
                for (l = 1; l <= NYCat; ++l) {

                    if (PMOPT == 3) {
                        par[0] = 1.0;
                        if (k == 1)  
                            k1 = 0;
                        else  
                            k1 = jptr + ((k - 1) * (k - 2) / 2) + 1;
                        if (l == 1)  
                            l1 = 0;
                        else  
                            l1 = jptr + ((l - 1) * (l - 2) / 2) + 1;

                        m = imin(k,l);
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

                            if (PMOPT == 2) {
                                s = rexp(s);
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
            QRProb[(i - 1) * ny1 + j] = tmp;
        }   
    }
    /***
    printf1("QRProb\n");
    for (i = 1; i <= ny1; ++i) {
        for (j = 1; j <= ny1; ++j)
            printf1("%lg ",QRProb[(i - 1) * ny1 + j]);
        newline();
    }
    ***/

    /* change into correlation matrix */

    for (i = 1; i <= ny1; ++i) {
        tmp = QRProb[(i - 1) * ny1 + i];
        if (tmp <= 0.0) {
            NPFlgs[9] += 1;
            return(1);
        }
        ss[i] = sqrt(tmp);
    }
    for (i = 1; i <= ny1; ++i) {
        for (j = 1; j <= ny1; ++j)  
            QRProb[(i - 1) * ny1 + j] /= (ss[i] * ss[j]);
    }
    /***
    printf1("CORR \n");
    for (i = 1; i <= ny1; ++i) {
        for (j = 1; j <= ny1; ++j)
            printf1("%lg ",QRProb[(i - 1) * ny1 + j]);
        newline();
    }
    ***/

    k = dmv(ny1,PMNHP,DMVArg,QRProb,prob,PMEPS,&tmp1);
    return(k);
}

/* ------------------------------------------------------------------------ */
/*  fn_probit4                                                              */
/*      Likelihood calculation for simultaneous binary probit model.        */
/*                                                                          */
/*      Parameterization depends on PMOPT (1,2,3).                          */
/*                                                                          */  

int fn_probit4(void)
{
    register int i,j;
    int wave,icase,err,cat;                 
    double wt,xi,tmp;           

    wt = 1.0;   NOCUsed = 0;

    for (icase = 0; icase < NOC; ++icase) {

        for (wave = 0; wave < NWave; ++wave) {

            QRTmp[wave + 1] = 1.0;

            if (!PDMCFlg || PDMCIdx[icase * NWave + wave]) {

                tmp = 0.0;

                if (IntFlg)
                    i = -1;
                else
                    i = 0;

                j = wave * NPX1;
                for ( ; i < NPX; ++i) {
                    if (i < 0)
                        xi = 1.0;
                    else  
                        xi = get_data(PXVar[i][wave],icase);
                    tmp += xi * TPar[++j];
                }
                cat = YCatI[(int)get_data(PYVar[wave],icase)];   /* internal category */
                if (cat == 0) {
                    tmp *= -1.0;
                    QRTmp[wave + 1] = -1.0;
                }
                NOCUsed++;
            }
            else
                tmp = 12.0;         /* some high value */

            DMVArg[wave + 1] = tmp;
        }
        /****
        printf1("DMVArg: "); 
        for (i = 1; i <= NWave; ++i)
            printf1("%lg ",DMVArg[i]);
        printf1("\nQRTmp: "); 
        for (i = 1; i <= NWave; ++i)
            printf1("%lg ",QRTmp[i]);
        newline();
        ***/

        err = fn_prob4(TPar,QRTmp,DMVArg,&tmp);

        if (err) {
            if (err < 0) {
                NPFlgs[9] += 1;
                return(1);
            }
            NPFlgs[10] += 1;
        }
        if (WIVar >= 0)
            wt = get_data(WIVar,icase) * WNorm;

        if (LFunc)   
            FTmp += wt * rlog(tmp);                        
    }
    QR_INIT = 0;
    checkov();
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

int fn_prob4(double *par,double *t,double *arg,double *prob)
{
    register int k,l,kk;
    int jptr,m,k1,l1;
    double s,tmp;
    double ss[DMVMax + 1];

    /* set up covariance matrix in QRProb */

    jptr = NWave * NPX1;    /* points to sigma in par */

    for (k = 1; k <= NWave; ++k) {
        for (l = 1; l <= NWave; ++l) {

            if (PMOPT == 3) {
                par[0] = 1.0;
                if (k == 1)  
                    k1 = 0;
                else  
                    k1 = jptr + ((k - 1) * (k - 2) / 2) + 1;
                if (l == 1)  
                    l1 = 0;
                else  
                    l1 = jptr + ((l - 1) * (l - 2) / 2) + 1;

                m = imin(k,l);
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

                    if (PMOPT == 2) {
                        s = rexp(s);
                        s = 2.0 * s / (1.0 + s) - 1.0;
                    }
                }
            }
            s *= t[k] * t[l];
            QRProb[(k - 1) * NWave + l] = s;   
        }   
    }
    /******* 
    printf1("QRProb\n");
    for (k = 1; k <= NWave; ++k) {
        for (l = 1; l <= NWave; ++l)
            printf1("%lg ",QRProb[(k - 1) * NWave + l]);
        newline();
    }
    ****/
  
    /* change into correlation matrix */

    for (k = 1; k <= NWave; ++k) {
        tmp = QRProb[(k - 1) * NWave + k];
        if (tmp <= 0.0) {
            NPFlgs[9] += 1;
            return(1);
        }
        ss[k] = sqrt(tmp);
    }
    for (k = 1; k <= NWave; ++k) {
        for (l = 1; l <= NWave; ++l)  
            QRProb[(k - 1) * NWave + l] /= (ss[k] * ss[l]);
    }
    /****
    printf1("CORR \n");
    for (k = 1; k <= NWave; ++k) {
        for (l = 1; l <= NWave; ++l)
            printf1("%lg ",QRProb[(k - 1) * NWave + l]);
        newline();
    }
    ****/
  
    k = dmv(NWave,PMNHP,arg,QRProb,prob,PMEPS,&tmp);
    return(k);
}

/* ------------------------------------------------------------------------ */
/*  fn_rmod1                                                                */
/*      Likelihood calculation for simple Rasch model.                      */

int fn_rmod1(void)
{
    register int i,j;
    int xij;
    double tmp,p;

    for (i = 1; i <= RMPN; ++i) {
        tmp = 0.0;
        for (j = 1; j <= RMPV; ++j) {
            xij = RMPAT[(i - 1) * RMPV + j];

            p = TPar[j] + TPar[RMPV + j];
            tmp += (double)xij * p - rlog(1.0 + rexp(p));       
        }
        if (LFunc) {
            FTmp += tmp * RMPF[i];                       
        }
    }
    checkov();
    return(0);
}
