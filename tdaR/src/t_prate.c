/****************************************************************************/
/*  t_prate                                                                 */
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
#include "t_edat.h"
#include "t_rate.h"
#include "t_ml.h"
#include "t_gf.h"
#include "t_cdf.h"
#include "t_int.h"
#include "t_fnrta.h"
#include "tda_context.h"

/*  functions in t_prate.c */

int prate(TDAContext *ctx, char *pcmd,int mod);
void prn_prs(TDAContext *ctx, int mod);
void prstab(TDAContext *ctx, int mod,int idx);
int get_rate(TDAContext *ctx, int mod,double time,double *mt, double *r,double *s,double *d,int idx);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */


/* ------------------------------------------------------------------------ */
/*  prate(pcmd,mod)     Print estimated transition rates etc.               */
/*                      Return 0 if OK, otherwise -1.                       */

int prate(TDAContext *ctx, char *pcmd,int mod)
{
    register int j;
    int err,k,l,vi,xaflg,tabflg,fflag;
    register char *p,*q,*pp;
    char c,*fn;
    char fname[FNMaxLen + 1],vname[VNLMax + 1];
    double x;

    err = -1;
    printf1(ctx, "Calculating estimated rates.\n");

    /**********
    for (i = 1; i <= NParm; ++i) {
        j = PIdx[i];
        printf1(ctx, "i=%d PIdx=%d  ",i,j);
        if (j > 0)
            printf1(ctx, " %s",VName[j - 1]);
        printf1(ctx, "\n");
    }
    *********************/

    if (mod == MCOX) {
        fflag = 0;
        if (alloc_acw(ctx, ctx->NTran1 + 1))
            goto PRTFin;
    }
    fflag = k = xaflg = tabflg = 0;
    fname[0] = '\0';
    p = pcmd;
    while (*p) {
        if (!strncmp(p,"prate",5)) {
            tabflg = 0;
            pp = q = p + 5;
            if (*q == '(')  
                pp = q = skip_blev(ctx, p + 5);

            if (*q++ != '=') {
                p_err(ctx, -1,1);
                goto PRTFin;
            }
            fn = q;
            q = skip_com(ctx, q);
            c = *q;
            *q = '\0';
            printf1(ctx, "[%2d] %s\n",++k,p);

            if (strcmp(fn,fname)) {     /* new filename */

                if (fflag) {
                    fclose(ctx->PRTFd);
                    fflag = 0;
                }
                if (!(ctx->PRTFd = fopen(fn,OPEN_WR))) {    
                    printf1(ctx, "Error: can't open: %s\n",fn);
                    goto PRTCONT;
                }
                fflag = 1;
                strcpy(fname,fn);
                ctx->PRTIdx = 0;
                if (mod == MCOX)
                    fprintf(ctx->PRTFd,"# Cox Model. Baseline Rate Calculation.\n");
                else
                    fprintf(ctx->PRTFd,"# Estimated Rates.\n");
            }
            xaflg = 0;

            p += 5;
            while (++p < pp) {
                if (sscanf(p,"tab=%lf(%lf)%lf",&ctx->PRTTA,&ctx->PRTTD,&ctx->PRTTB) == 3 &&
                            ctx->PRTTA >= 0.0 && ctx->PRTTD > 0.0 && ctx->PRTTB > ctx->PRTTA) {
                    p = skip_dbl(ctx, p + 4);
                    p = skip_dbl(ctx, p + 1);
                    p = skip_dbl(ctx, p + 1);
                    tabflg = 1;
                }
                else {
                    if ((l = get_vnlen(ctx, p)) > 0 && (vi = get_vidx(ctx, p,vname)) >= 0 && *(p += l) == '=') {
                        if (xaflg == 0) {
                            xaflg = ctx->NParm + 1;
                            if (alloc_acx(ctx, xaflg))
                                goto PRTFin;
                        }
                        if (sscanf(++p,"%lg",&x) == 1) {
                            for (j = 1; j <= ctx->NParm; ++j) {
                                if (ctx->PIdx[j] == vi + 1 || (mod == MCOX && -ctx->PIdx[j] == vi + 1))
                                    ctx->AcX[j] = x;
                            }
                            p = skip_dbl(ctx, p);
                        }
                        else {
                            p_err(ctx, -1,1);
                            goto PRTCONT;
                        }
                    }
                    else if (*p != ')') {
                        p_err(ctx, -1,1);
                        goto PRTCONT;
                    }
                }
            }

            if (xaflg == 0) {
                xaflg = ctx->NParm + 1;
                if (alloc_acx(ctx, xaflg))
                    goto PRTFin;
            }
            if (alloc_actmp(ctx, ctx->MTerm + 1))
                goto PRTFin;

            if (tabflg == 0 && mod != MCOX) {
                printf1(ctx, "Error: need tab parameter for time axis.\n");
                goto PRTCONT;
            }

            prn_prs(ctx, mod);       /* print tables to output file */

PRTCONT:
            *q = c;
            p = q;
        }
        p++;
    }
    err = 0;

PRTFin:
    if (fflag)
        fclose(ctx->PRTFd);
    alloc_actmp(ctx, 0);
    alloc_acw(ctx, 0);
    alloc_acx(ctx, 0);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_prs(mod)    This function calculates the model terms for each       */
/*                  transition, and then calls prstab for printing.         */
/*                                                                          */  
/*  AcX[j] contains covariate values, j = 1,...,NParm.                      */
/*  AcTmp[k] is used to accumulate period-specific model terms.             */
/*  k = 1,...,MTerm.                                                        */
/*                                                                          */
/*  Note: different procedure for Cox models. Then calculate transition-    */
/*  specific xa terms in AcW and call fn_cox() with COXBLFlg.               */

void prn_prs(TDAContext *ctx, int mod)
{
    register int j,k;
    int sn,org,des,t,idx,iv;
    double tmp;

    if (mod != MCOX)
        fprintf(ctx->PRTFd,"\n# Time axis: %g (%g) %g\n",ctx->PRTTA,ctx->PRTTD,ctx->PRTTB);

    k = 1; idx = 0;

    for (t = 0; t < ctx->NTran1; ++t) {

        sn = ctx->SnTran1[t];
        org = ctx->OrgTran1[t];
        des = ctx->DesTran1[t];

        if (t > 0 && mod == MCOX)
            fprintf(ctx->PRTFd,"\n");

        fprintf(ctx->PRTFd,"\n# Idx SN Org Des MT Variable ");
        fprnchar(ctx, ctx->PRTFd,' ',ctx->VNameLen - 8,0);
        fprnchar(ctx, ctx->PRTFd,' ',ctx->PMTFmt1 - 10,0); fprintf(ctx->PRTFd,"     Coeff");
        fprnchar(ctx, ctx->PRTFd,' ',ctx->PMTFmt1 -  9,0); fprintf(ctx->PRTFd," Covariate\n# ");

        fprnchar(ctx, ctx->PRTFd,'-',18 + 2 * (ctx->PMTFmt1 + 1) + ctx->VNameLen,0);

        for (j = 1; j <= ctx->MTerm; ++j) {

            ctx->AcTmp[j] = 0.0;

            if (ctx->XMFlg[j]) {

                while (1) {

                    tmp = 0.0;  /* covariate value */

                    fprintf(ctx->PRTFd,"\n# %3d %2d %3d %3d  %c ",
                                              k,sn,org,des,(char)('A'+j-1));
                    
                    iv = ctx->PIdx[k];
                    if (iv > 0 || (iv < 0 && mod == MCOX)) {
                        if (mod == MCOX && iv < 0)
                            iv = -iv;

                        fprintf(ctx->PRTFd,"%s",ctx->VName[iv - 1]);
                        fprnchar(ctx, ctx->PRTFd,' ',ctx->VNameLen + 1 - (int)strlen(ctx->VName[iv - 1]),0);
                        tmp = ctx->AcX[k];
                    }
                    else {

                        if (j == 2 && !ctx->PIdx[k])
                            idx = k;

                        if (mod == MEXP1) {
                            fprintf(ctx->PRTFd,"Period-%-2d",-ctx->PIdx[k] + 1);
                            fprnchar(ctx, ctx->PRTFd,' ',ctx->VNameLen - 8,0);
                            if (!ctx->PIdx[k])
                                idx = k;
                        }
                        else if (j == 2 && (mod == MPOL || mod == MPOL1 ||
                                            mod == DLR  || mod == CLL)) {
 
                            fprintf(ctx->PRTFd,"Beta-%-2d",-ctx->PIdx[k] + 1);
                            fprnchar(ctx, ctx->PRTFd,' ',ctx->VNameLen - 6,0);
                            if (!ctx->PIdx[k])
                                idx = k;
                        }
                        else {
                            fprintf(ctx->PRTFd,"Constant");
                            fprnchar(ctx, ctx->PRTFd,' ',ctx->VNameLen - 7,0);
                        }
                        tmp = 1.0;
                    }
                    rt_fprintf_d(ctx, ctx->PRTFd,ctx->PMTFmtS,ctx->Par[k]);
                    rt_fprintf_d(ctx, ctx->PRTFd,ctx->PMTFmtS,tmp);
   
                    /* omit time-period constants in the    */
                    /* case of an MEXP1 model.              */
       
                    if (mod != MEXP1 || ctx->PIdx[k] > 0)
                        ctx->AcTmp[j] += tmp * ctx->Par[k];

                    if (!ctx->PIdx[++k] || (mod == MCOX && ctx->PIdx[k] < 0))
                        break;
                }
            }
        }
        fprintf(ctx->PRTFd,"\n");
            
        /* now the actual par * covariate values are in AcTmp[], */
        /* and the estimated rates etc. can be calculated        */
        /* and printed by prstab.                                */

        if (mod != MCOX) {
            prstab(ctx, mod,idx);               
            fprintf(ctx->PRTFd,"\n");
        }
        else {                              /* save exp(xa) */
            ctx->AcW[t] = rexp(ctx, ctx->AcTmp[1]);
        }
    }
    if (mod == MCOX) {
        ctx->COXBLFlg = 1;
        fn_cox(ctx);
        ctx->COXBLFlg = 0;
    }
}

/* ------------------------------------------------------------------------ */
/*  prstab(mod)                                                             */
/*                                                                          */
/*          If MPOL or MPOL1, idx is index of B-term parameters.            */
/*          If MEXP1, idx is index of A-term parameters.                    */

void prstab(TDAContext *ctx, int mod,int idx)
{
    double r,s,d,time;
       
    /*************
    fprintf(PRTFd,"print ... \n");
    fprintf(PRTFd,"time = %g (%g) %g\n",PRTTA,PRTTD,PRTTB);
    fprintf(PRTFd,"\nModel terms\n");
    for (j = 1; j <= MTerm; ++j) 
        fprintf(PRTFd,"%g ",AcTmp[j]);
    fprintf(PRTFd,"\n");
    ***********/
    
    fprintf(ctx->PRTFd,"\n# ID ");
    fprnchar(ctx, ctx->PRTFd,' ',ctx->PMMFmt1 - 4,0); fprintf(ctx->PRTFd,"Time");
    fprnchar(ctx, ctx->PRTFd,' ',ctx->PMMFmt1 - 5,0); fprintf(ctx->PRTFd,"Surv.F ");
    fprnchar(ctx, ctx->PRTFd,' ',ctx->PMMFmt1 - 7,0); fprintf(ctx->PRTFd,"Density");
    fprnchar(ctx, ctx->PRTFd,' ',ctx->PMMFmt1 - 3,0); fprintf(ctx->PRTFd,"Rate\n# ");
    fprnchar(ctx, ctx->PRTFd,'-',6 + 4 * ctx->PMMFmt1,0);
    fprintf(ctx->PRTFd,"\n");

    for (time = ctx->PRTTA; time <= ctx->PRTTB + ctx->EPSI1; time += ctx->PRTTD) {

        if (get_rate(ctx, mod,time,ctx->AcTmp,&r,&s,&d,idx)) {
            fprintf(ctx->PRTFd,"# Calculation not possible.\n");
            break;
        }
        fprintf(ctx->PRTFd,"%4d ",ctx->PRTIdx);
        rt_fprintf_d(ctx, ctx->PRTFd,ctx->PMMFmtS,time);
        rt_fprintf_d(ctx, ctx->PRTFd,ctx->PMMFmtS,s);
        rt_fprintf_d(ctx, ctx->PRTFd,ctx->PMMFmtS,d);
        rt_fprintf_d(ctx, ctx->PRTFd,ctx->PMMFmtS,r);
        fprintf(ctx->PRTFd,"\n");
#ifdef TDA_R_PACKAGE
        {
            double erow[5];
            erow[0] = (double)ctx->PRTIdx;
            erow[1] = time;
            erow[2] = s;
            erow[3] = d;
            erow[4] = r;
            tda_export_row(ctx, "prate.table", erow, 5);
        }
#endif
    }
    ctx->PRTIdx++;
#ifdef TDA_R_PACKAGE
    /* one export per printed profile (per idx), suffixed */
    tda_export_flush(ctx, "prate.table");
#endif
}

/* ------------------------------------------------------------------------ */
/*  get_rate(mod,time,sp,mt,&r,&s,&d,idx)                                   */
/*                                                                          */

int get_rate(TDAContext *ctx, int mod,double time,double *mt, double *r,double *s,double *d,int idx)
{
    register int i;
    int err;
    double xa,xb,xc,rate,surv,dens,a,b,c,z,tmp,ab,ca,tmp1,tmp2;
    double k1,k2,g1,ez;

    *r = 0.0;
    *s = 0.0;
    *d = 0.0;
      
    xa = mt[1];
    xb = mt[2];
    xc = mt[3];
    err = 0;

    switch (mod) {

        case  MEXP:     rate = rexp(ctx, xa);    
                        surv = rexp(ctx, -time * rate);
                        dens = surv * rate;
                        break;
                     
        case  MEXP1:    tmp = 0.0;
                        for (i = 1; i < ctx->PMNTP; ++i) {
                            rate = rexp(ctx, xa + ctx->TPar[idx + i - 1]);
                            if (time < ctx->PMTP[i]) {
                                tmp += rate * (time - ctx->PMTP[i - 1]);
                                break;
                            }
                            tmp += rate * (ctx->PMTP[i] - ctx->PMTP[i - 1]);
                        }
                        if (i >= ctx->PMNTP) {
                            rate = rexp(ctx, xa + ctx->TPar[idx + ctx->PMNTP - 1]);
                            tmp += rate * (time - ctx->PMTP[ctx->PMNTP - 1]);
                        }
                        surv = rexp(ctx, -tmp);
                        dens = surv * rate;
                        break;
        
        case  MPOL:     a = rexp(ctx, xa);
                        c = z = 0.0;
                        for (i = ctx->PMDEG - 1; i >= 0; --i) {
                            z += ctx->Par[idx + i];
                            z *= time;
                            c += ctx->Par[idx + i] / (double)(i + 2);
                            c *= time;
                        }
                        rate = a + z;
                        surv = rexp(ctx, -(a + c) * time);
                        dens = rate * surv;
                        break;
         
        case  MPOL1:    z = 0.0;
                        for (i = ctx->PMDEG - 1; i >= 0; --i) {
                            z += ctx->Par[idx + i];
                            z *= time;
                        }
                        rate = rexp(ctx, xa + z);
                        if (ctx->PMDEG > 0)
                            tmp = ni_gen(ctx, 0.0,time,2,0,0,&err);
                        else
                            tmp = time;

                        surv = rexp(ctx, -rexp(ctx, xa) * tmp);
                        dens = rate * surv;
                        break;
     
        case  MGM:      if (ctx->XMFlg[1])
                            a = rexp(ctx, xa);
                        else
                            a = 0.0;
                        
                        if (ctx->XMFlg[2]) {
                            b = rexp(ctx, xb);
                            c = xc;
                            tmp = b * rexp(ctx, c * time);
                            rate = a + tmp;
                            z = -a * time;
                            if (fabs(c) > ctx->EPSI)
                                z -= (tmp - b) / c;
                            surv = rexp(ctx, z);
                        }
                        else {
                            rate = a;
                            surv = rexp(ctx, -a * time);
                        }
                        dens = rate * surv;
                        break;
   
        case  MWEI:     if (time < ctx->EPSI1) {
                            rate = dens = 0.0;
                            surv = 1.0;
                            break;
                        }
                        a = rexp(ctx, xa);
                        b = rexp(ctx, xb);        
                        tmp = pow(a * time,b);
                        surv = rexp(ctx, -tmp);
                        rate = b * tmp / time;
                        dens = rate * surv;
                        break;
   
        case  MSIC:     a = rexp(ctx, xa);    
                        b = rexp(ctx, xb);
                        z = rexp(ctx, -time / b);
                        rate = a * time * z;
                        surv = rexp(ctx, -a * b * (b - (time + b) * z));
                        dens = rate * surv;
                        break;
       
        case  MLL:      if (time < ctx->EPSI1) {
                            rate = dens = 0.0;
                            surv = 1.0;
                            break;
                        }
                        a = rexp(ctx, xa);
                        b = rexp(ctx, xb);    
                        ab = pow(a,b);
                        surv = 1.0 / (1.0 + ab * pow(time,b));
                        rate = b * ab * pow(time,b) * surv / time;
                        dens = rate * surv;
                        break;

        case  MLL2:     if (time < ctx->EPSI1) {
                            rate = dens = 0.0;
                            surv = 1.0;
                            break;
                        }
                        a = rexp(ctx, xa);
                        b = rexp(ctx, xb);
                        c = rexp(ctx, xc);
                        ca = c / a;
                        ab = pow(a,b);
                        tmp = 1.0 / (1.0 + ab * pow(time,b));
                        surv = pow(tmp,ca);
                        rate = ca * b * ab * pow(time,b) * tmp / time;
                        dens = rate * surv;
                        break;
             
        case  MLN:      if (ctx->XMFlg[3])       /* if three parameter model */
                            c = rexp(ctx, xc);
                        else
                            c = 0.0;
                        
                        if (time > c) {
                            tmp = time - c;
                            a = xa;
                            b = rexp(ctx, xb);    
                            z = (rlog(ctx, tmp) - a) / b;
                            b *= tmp;
                            surv = 1.0 - cdnf(ctx, z);
                            dens = dnf(ctx, z) / b;
                            rate = dens / surv;
                        }
                        else {
                            rate = dens = 0.0;
                            surv = 1.0;
                        }
                        break;
               
        case  MGAM:     rate = dens = 0.0;
                        surv = 1.0;
                        if (time < ctx->EPSI1)
                            break;
                               
                        k1 = sqrt(ctx->PMKGam);
                        k2 = (ctx->PMKGam - 0.5) * rlog(ctx, ctx->PMKGam);
                        g1 = loggam(ctx, ctx->PMKGam,&err); 
                        if (err)
                            break;
                        b = rexp(ctx, xb);    
                        z = (rlog(ctx, time) - xa) / b;
                        b *= time;
                        ez = rexp(ctx, z / k1);
                        dens = rexp(ctx, k2 - xb + k1 * z - ctx->PMKGam * ez - g1) / time;
                        surv = 1.0 - icgam(ctx, ctx->PMKGam * ez,ctx->PMKGam,&err);
                        if (err)
                            break;
                        if (surv > 0.0)
                            rate = dens / surv;
                        break;
       
        case  MIG:      if (time < ctx->EPSI1) {
                            rate = dens = 0.0;
                            surv = 1.0;
                            break;
                        }
                        rate = 0.0;
                        a = xa;
                        b = rexp(ctx, xb);
                        tmp = sqrt(time);
                        z = (a * time - 1.0) / (b * tmp);
                        dens = dnf(ctx, z) / (b * time * tmp);
                        tmp1 = 1.0 - cdnf(ctx, z);
                        z = (a * time + 1.0) / (b * tmp);
                        tmp2 = 1.0 - cdnf(ctx, z);
                        surv = tmp1 - tmp2 * rexp(ctx, 2.0 * a / (b * b));
                        if (surv > 0.0)
                            rate = dens / surv;
                        break;
/**     

        case  DLR:  
        case  CLL:      if (init) {
                            for (k = 1; k < (int)time; ++k) {
                                ab = xa;
                                if (PMDEG) {
                                    j = n;
                                    tmp = (double)k; 
                                    ab += TPar[j] * tmp;                                              
                                    while (PIdx[++j]) {                                               
                                        tmp *= (double)k;
                                        ab += TPar[j] * tmp;
                                    }
                                }
                                tmp1 = rexp(ctx, ab);
                                if (mod == DLR)
                                    rate = tmp1 / (1.0 + tmp1);
                                else
                                    rate = 1.0 - rexp(ctx, -tmp1);
                                surv1 *= (1.0 - rate);
                            }
                            init = 0;
                        }
                        ab = xa;
                        if (PMDEG) {
                            j = n;
                            tmp = (double)time; 
                            ab += TPar[j] * tmp;                                              
                            while (PIdx[++j]) {                                               
                                tmp *= (double)time;
                                ab += TPar[j] * tmp;
                            }
                        }
                        tmp1 = rexp(ctx, ab);
                        if (mod == DLR)
                            rate = tmp1 / (1.0 + tmp1);
                        else
                            rate = 1.0 - rexp(ctx, -tmp1);

                        dens = rate * surv1;
                        surv1 *= (1.0 - rate);
                        surv = surv1;                    
                        break;
*******/
        default:        err = -1;   
                        break;  
    }

    if (err == 0) {
        *r = rate;
        *s = surv;
        *d = dens;
    }
    return(err);
}





