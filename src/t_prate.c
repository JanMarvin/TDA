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

/*  functions in t_prate.c */

int prate(char *pcmd,int mod);
void prn_prs(int mod);
void prstab(int mod,int idx);
int get_rate(int mod,double time,double *mt,
    double *r,double *s,double *d,int idx);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

FILE *PRTFd;            /* output file name                                 */
int PRTIdx = 0;         /* index for table                                  */
double PRTTA = 0.0;     /* time axis                                        */
double PRTTB = 0.0;
double PRTTD = 0.0;

/* ------------------------------------------------------------------------ */
/*  prate(pcmd,mod)     Print estimated transition rates etc.               */
/*                      Return 0 if OK, otherwise -1.                       */

int prate(char *pcmd,int mod)
{
    register int j;
    int err,k,l,vi,xaflg,tabflg,fflag;
    register char *p,*q,*pp;
    char c,*fn;
    char fname[FNMaxLen + 1],vname[VNLMax + 1];
    double x;

    err = -1;
    printf1("Calculating estimated rates.\n");

    /**********
    for (i = 1; i <= NParm; ++i) {
        j = PIdx[i];
        printf1("i=%d PIdx=%d  ",i,j);
        if (j > 0)
            printf1(" %s",VName[j - 1]);
        printf1("\n");
    }
    *********************/

    if (mod == MCOX) {
        if (alloc_acw(NTran1 + 1))
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
                pp = q = skip_blev(p + 5);

            if (*q++ != '=') {
                p_err(-1,1);
                goto PRTFin;
            }
            fn = q;
            q = skip_com(q);
            c = *q;
            *q = '\0';
            printf1("[%2d] %s\n",++k,p);

            if (strcmp(fn,fname)) {     /* new filename */

                if (fflag) {
                    fclose(PRTFd);
                    fflag = 0;
                }
                if (!(PRTFd = fopen(fn,OPEN_WR))) {    
                    printf1("Error: can't open: %s\n",fn);
                    goto PRTCONT;
                }
                fflag = 1;
                strcpy(fname,fn);
                PRTIdx = 0;
                if (mod == MCOX)
                    fprintf(PRTFd,"# Cox Model. Baseline Rate Calculation.\n");
                else
                    fprintf(PRTFd,"# Estimated Rates.\n");
            }
            xaflg = 0;

            p += 5;
            while (++p < pp) {
                if (sscanf(p,"tab=%lf(%lf)%lf",&PRTTA,&PRTTD,&PRTTB) == 3 &&
                            PRTTA >= 0.0 && PRTTD > 0.0 && PRTTB > PRTTA) {
                    p = skip_dbl(p + 4);
                    p = skip_dbl(p + 1);
                    p = skip_dbl(p + 1);
                    tabflg = 1;
                }
                else {
                    if ((l = get_vnlen(p)) > 0 && (vi = get_vidx(p,vname)) >= 0 && *(p += l) == '=') {
                        if (xaflg == 0) {
                            xaflg = NParm + 1;
                            if (alloc_acx(xaflg))
                                goto PRTFin;
                        }
                        if (sscanf(++p,"%lg",&x) == 1) {
                            for (j = 1; j <= NParm; ++j) {
                                if (PIdx[j] == vi + 1 || (mod == MCOX && -PIdx[j] == vi + 1))
                                    AcX[j] = x;
                            }
                            p = skip_dbl(p);
                        }
                        else {
                            p_err(-1,1);
                            goto PRTCONT;
                        }
                    }
                    else if (*p != ')') {
                        p_err(-1,1);
                        goto PRTCONT;
                    }
                }
            }

            if (xaflg == 0) {
                xaflg = NParm + 1;
                if (alloc_acx(xaflg))
                    goto PRTFin;
            }
            if (alloc_actmp(MTerm + 1))
                goto PRTFin;

            if (tabflg == 0 && mod != MCOX) {
                printf1("Error: need tab parameter for time axis.\n");
                goto PRTCONT;
            }

            prn_prs(mod);       /* print tables to output file */

PRTCONT:
            *q = c;
            p = q;
        }
        p++;
    }
    err = 0;

PRTFin:
    if (fflag)
        fclose(PRTFd);
    alloc_actmp(0);
    alloc_acw(0);
    alloc_acx(0);
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

void prn_prs(int mod)
{
    register int j,k;
    int sn,org,des,t,idx,iv;
    double tmp;

    if (mod != MCOX)
        fprintf(PRTFd,"\n# Time axis: %g (%g) %g\n",PRTTA,PRTTD,PRTTB);

    k = 1; idx = 0;

    for (t = 0; t < NTran1; ++t) {

        sn = SnTran1[t];
        org = OrgTran1[t];
        des = DesTran1[t];

        if (t > 0 && mod == MCOX)
            fprintf(PRTFd,"\n");

        fprintf(PRTFd,"\n# Idx SN Org Des MT Variable ");
        fprnchar(PRTFd,' ',VNameLen - 8,0);
        fprnchar(PRTFd,' ',PMTFmt1 - 10,0); fprintf(PRTFd,"     Coeff");
        fprnchar(PRTFd,' ',PMTFmt1 -  9,0); fprintf(PRTFd," Covariate\n# ");

        fprnchar(PRTFd,'-',18 + 2 * (PMTFmt1 + 1) + VNameLen,0);

        for (j = 1; j <= MTerm; ++j) {

            AcTmp[j] = 0.0;

            if (XMFlg[j]) {

                while (1) {

                    tmp = 0.0;  /* covariate value */

                    fprintf(PRTFd,"\n# %3d %2d %3d %3d  %c ",
                                              k,sn,org,des,(char)('A'+j-1));
                    
                    iv = PIdx[k];
                    if (iv > 0 || (iv < 0 && mod == MCOX)) {
                        if (mod == MCOX && iv < 0)
                            iv = -iv;

                        fprintf(PRTFd,"%s",VName[iv - 1]);
                        fprnchar(PRTFd,' ',VNameLen + 1 - strlen(VName[iv - 1]),0);
                        tmp = AcX[k];
                    }
                    else {

                        if (j == 2 && !PIdx[k])
                            idx = k;

                        if (mod == MEXP1) {
                            fprintf(PRTFd,"Period-%-2d",-PIdx[k] + 1);
                            fprnchar(PRTFd,' ',VNameLen - 8,0);
                            if (!PIdx[k])
                                idx = k;
                        }
                        else if (j == 2 && (mod == MPOL || mod == MPOL1 ||
                                            mod == DLR  || mod == CLL)) {
 
                            fprintf(PRTFd,"Beta-%-2d",-PIdx[k] + 1);
                            fprnchar(PRTFd,' ',VNameLen - 6,0);
                            if (!PIdx[k])
                                idx = k;
                        }
                        else {
                            fprintf(PRTFd,"Constant");
                            fprnchar(PRTFd,' ',VNameLen - 7,0);
                        }
                        tmp = 1.0;
                    }
                    fprintf(PRTFd,PMTFmtS,Par[k]);
                    fprintf(PRTFd,PMTFmtS,tmp);
   
                    /* omit time-period constants in the    */
                    /* case of an MEXP1 model.              */
       
                    if (mod != MEXP1 || PIdx[k] > 0)
                        AcTmp[j] += tmp * Par[k];

                    if (!PIdx[++k] || (mod == MCOX && PIdx[k] < 0))
                        break;
                }
            }
        }
        fprintf(PRTFd,"\n");
            
        /* now the actual par * covariate values are in AcTmp[], */
        /* and the estimated rates etc. can be calculated        */
        /* and printed by prstab.                                */

        if (mod != MCOX) {
            prstab(mod,idx);               
            fprintf(PRTFd,"\n");
        }
        else {                              /* save exp(xa) */
            AcW[t] = rexp(AcTmp[1]);
        }
    }
    if (mod == MCOX) {
        COXBLFlg = 1;
        fn_cox();
        COXBLFlg = 0;
    }
}

/* ------------------------------------------------------------------------ */
/*  prstab(mod)                                                             */
/*                                                                          */
/*          If MPOL or MPOL1, idx is index of B-term parameters.            */
/*          If MEXP1, idx is index of A-term parameters.                    */

void prstab(int mod,int idx)
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
    
    fprintf(PRTFd,"\n# ID ");
    fprnchar(PRTFd,' ',PMMFmt1 - 4,0); fprintf(PRTFd,"Time");
    fprnchar(PRTFd,' ',PMMFmt1 - 5,0); fprintf(PRTFd,"Surv.F ");
    fprnchar(PRTFd,' ',PMMFmt1 - 7,0); fprintf(PRTFd,"Density");
    fprnchar(PRTFd,' ',PMMFmt1 - 3,0); fprintf(PRTFd,"Rate\n# ");
    fprnchar(PRTFd,'-',6 + 4 * PMMFmt1,0);
    fprintf(PRTFd,"\n");

    for (time = PRTTA; time <= PRTTB + EPSI1; time += PRTTD) {

        if (get_rate(mod,time,AcTmp,&r,&s,&d,idx)) {
            fprintf(PRTFd,"# Calculation not possible.\n");
            break;
        }
        fprintf(PRTFd,"%4d ",PRTIdx);
        fprintf(PRTFd,PMMFmtS,time);
        fprintf(PRTFd,PMMFmtS,s);
        fprintf(PRTFd,PMMFmtS,d);
        fprintf(PRTFd,PMMFmtS,r);
        fprintf(PRTFd,"\n");
    }
    PRTIdx++;
}

/* ------------------------------------------------------------------------ */
/*  get_rate(mod,time,sp,mt,&r,&s,&d,idx)                                   */
/*                                                                          */

int get_rate(int mod,double time,double *mt,
    double *r,double *s,double *d,int idx)
{
    register int i;
    int err;
    double xa,xb,xc,xd,rate,surv,dens,a,b,c,z,tmp,ab,ca,tmp1,tmp2;
    double k1,k2,g1,ez;

    *r = 0.0;
    *s = 0.0;
    *d = 0.0;
      
    xa = mt[1];
    xb = mt[2];
    xc = mt[3];
    xd = mt[4];
    err = 0;

    switch (mod) {

        case  MEXP:     rate = rexp(xa);    
                        surv = rexp(-time * rate);
                        dens = surv * rate;
                        break;
                     
        case  MEXP1:    tmp = 0.0;
                        for (i = 1; i < PMNTP; ++i) {
                            rate = rexp(xa + TPar[idx + i - 1]);
                            if (time < PMTP[i]) {
                                tmp += rate * (time - PMTP[i - 1]);
                                break;
                            }
                            tmp += rate * (PMTP[i] - PMTP[i - 1]);
                        }
                        if (i >= PMNTP) {
                            rate = rexp(xa + TPar[idx + PMNTP - 1]);
                            tmp += rate * (time - PMTP[PMNTP - 1]);
                        }
                        surv = rexp(-tmp);
                        dens = surv * rate;
                        break;
        
        case  MPOL:     a = rexp(xa);
                        c = z = 0.0;
                        for (i = PMDEG - 1; i >= 0; --i) {
                            z += Par[idx + i];
                            z *= time;
                            c += Par[idx + i] / (double)(i + 2);
                            c *= time;
                        }
                        rate = a + z;
                        surv = rexp(-(a + c) * time);
                        dens = rate * surv;
                        break;
         
        case  MPOL1:    z = 0.0;
                        for (i = PMDEG - 1; i >= 0; --i) {
                            z += Par[idx + i];
                            z *= time;
                        }
                        rate = rexp(xa + z);
                        if (PMDEG > 0)
                            tmp = ni_gen(0.0,time,2,0,0,&err);
                        else
                            tmp = time;

                        surv = rexp(-rexp(xa) * tmp);
                        dens = rate * surv;
                        break;
     
        case  MGM:      if (XMFlg[1])
                            a = rexp(xa);
                        else
                            a = 0.0;
                        
                        if (XMFlg[2]) {
                            b = rexp(xb);
                            c = xc;
                            tmp = b * rexp(c * time);
                            rate = a + tmp;
                            z = -a * time;
                            if (fabs(c) > EPSI)
                                z -= (tmp - b) / c;
                            surv = rexp(z);
                        }
                        else {
                            rate = a;
                            surv = rexp(-a * time);
                        }
                        dens = rate * surv;
                        break;
   
        case  MWEI:     if (time < EPSI1) {
                            rate = dens = 0.0;
                            surv = 1.0;
                            break;
                        }
                        a = rexp(xa);
                        b = rexp(xb);        
                        tmp = pow(a * time,b);
                        surv = rexp(-tmp);
                        rate = b * tmp / time;
                        dens = rate * surv;
                        break;
   
        case  MSIC:     a = rexp(xa);    
                        b = rexp(xb);
                        z = rexp(-time / b);
                        rate = a * time * z;
                        surv = rexp(-a * b * (b - (time + b) * z));
                        dens = rate * surv;
                        break;
       
        case  MLL:      if (time < EPSI1) {
                            rate = dens = 0.0;
                            surv = 1.0;
                            break;
                        }
                        a = rexp(xa);
                        b = rexp(xb);    
                        ab = pow(a,b);
                        surv = 1.0 / (1.0 + ab * pow(time,b));
                        rate = b * ab * pow(time,b) * surv / time;
                        dens = rate * surv;
                        break;

        case  MLL2:     if (time < EPSI1) {
                            rate = dens = 0.0;
                            surv = 1.0;
                            break;
                        }
                        a = rexp(xa);
                        b = rexp(xb);
                        c = rexp(xc);
                        ca = c / a;
                        ab = pow(a,b);
                        tmp = 1.0 / (1.0 + ab * pow(time,b));
                        surv = pow(tmp,ca);
                        rate = ca * b * ab * pow(time,b) * tmp / time;
                        dens = rate * surv;
                        break;
             
        case  MLN:      if (XMFlg[3])       /* if three parameter model */
                            c = rexp(xc);
                        else
                            c = 0.0;
                        
                        if (time > c) {
                            tmp = time - c;
                            a = xa;
                            b = rexp(xb);    
                            z = (rlog(tmp) - a) / b;
                            b *= tmp;
                            surv = 1.0 - cdnf(z);
                            dens = dnf(z) / b;
                            rate = dens / surv;
                        }
                        else {
                            rate = dens = 0.0;
                            surv = 1.0;
                        }
                        break;
               
        case  MGAM:     rate = dens = 0.0;
                        surv = 1.0;
                        if (time < EPSI1)
                            break;
                               
                        k1 = sqrt(PMKGam);
                        k2 = (PMKGam - 0.5) * rlog(PMKGam);
                        g1 = loggam(PMKGam,&err); 
                        if (err)
                            break;
                        b = rexp(xb);    
                        z = (rlog(time) - xa) / b;
                        b *= time;
                        ez = rexp(z / k1);
                        dens = rexp(k2 - xb + k1 * z - PMKGam * ez - g1) / time;
                        surv = 1.0 - icgam(PMKGam * ez,PMKGam,&err);
                        if (err)
                            break;
                        if (surv > 0.0)
                            rate = dens / surv;
                        break;
       
        case  MIG:      if (time < EPSI1) {
                            rate = dens = 0.0;
                            surv = 1.0;
                            break;
                        }
                        rate = 0.0;
                        a = xa;
                        b = rexp(xb);
                        tmp = sqrt(time);
                        z = (a * time - 1.0) / (b * tmp);
                        dens = dnf(z) / (b * time * tmp);
                        tmp1 = 1.0 - cdnf(z);
                        z = (a * time + 1.0) / (b * tmp);
                        tmp2 = 1.0 - cdnf(z);
                        surv = tmp1 - tmp2 * rexp(2.0 * a / (b * b));
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
                                tmp1 = rexp(ab);
                                if (mod == DLR)
                                    rate = tmp1 / (1.0 + tmp1);
                                else
                                    rate = 1.0 - rexp(-tmp1);
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
                        tmp1 = rexp(ab);
                        if (mod == DLR)
                            rate = tmp1 / (1.0 + tmp1);
                        else
                            rate = 1.0 - rexp(-tmp1);

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





