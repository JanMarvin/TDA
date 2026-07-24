/****************************************************************************/
/*  t_ltb                                                                   */
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

#include "tda.h"
#include "t_gen.h"
#include "t_pgen.h"
#include "t_parm.h"
#include "t_edat.h"
#include "t_gdat.h"
#include "t_var.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_ltb.c                                                    */

int ltb(void);
int lifetab(void);
void prltb1_head(int sn,int org,int grp,int n,double w,char *df);
void prltb1 (int l,double ne,double zz,double rs,double *ee,char *df);
void prltb2_head(int ndes,char *df);
void prltb2 (int l, double surv, double sterr, double *dens, double derr,
    double *rate, double rerr, int ndes, char *df);
void prstar(void);

/* ------------------------------------------------------------------------ */
/*  ltb()           Life table calculation. ltb command in CmdBuf.          */
/*                  Return 0 if OK, otherwise -1.                           */

int ltb(void)
{
    int err = -1;

    if (check_cmd(1))
        return(-1);

    printf1("Life table estimation. Current memory: %d bytes.\n",MemReq);

    if (EDAvail == 0) {
        p_err(-15,1);
        return(0);
    }
    if (parm(CmdBuf + 3,1,1))     /* get parameters */
        goto LTBFin;

    if (PMNTP == 0) {
        p_err(-17,1);
        goto LTBFin;
    }
    if (PMFmtF == 0)            /* default print format */
        pmfmt(7,5);

    if (NSP > 0)  
        p_warn(-2,1);

    newline();
    err = lifetab();
    if (err) {
        if (err < 0)
            p_err(-2,1);
        else {
            printf1("\nWarning: at least one starting time is not zero.\n");
            printf1("Results are probably wrong.\n");
            err = -1;
        }
    }

LTBFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  lifetab()       Lifetable calculation.                                  */
/*                  Return 0 if OK, -1 if insufficient memory               */
/*                  Return 1 if not all starting times are zero.            */

int lifetab(void)
{
    register int k,l,grp;
    int err,err1,nt1,isn,iorg,ip,icase,sn,org,des,ndes,n0,ntab,spl,nspl;
    int a_ne,a_zz,a_rs,a_es,a_ee,a_subd,a_dens,a_rate,a_df;
    double tf,ts,wt,tmp,tmp1,pl,ql,med,surv,surv1,rate1,sterr,derr,rerr,sum,tt;
    double *es,*ee,*zz,*ne,*rs,*subd,*dens,*rate;
    char *df;

    a_df = a_rate = a_dens = a_subd = a_ee = a_es = a_rs = a_zz = a_ne = 0;
    err1 = ntab = 0;

    /*  memory allocation:                                          
        ne = number entering an interval                            
        zz = number of censored episodes                            
        rs = risk set                                               
        es = number of events                                       
        ee = number of events, diff. according to destination state 
        df = flags for existing destination states                  
        subd = subdistribution function                             
        dens = subdensity function                                  
        rate = transition rates. */                                 

    err = -1;
    nt1 = PMNTP + 1;

    if (!(ne = (double *)calloc(nt1,sizeof(double))))  
        goto LTBFin;
    a_ne = 1;
    memrq(nt1,sizeof(double));

    if (!(zz = (double *)calloc(nt1,sizeof(double))))  
        goto LTBFin;
    a_zz = 1;
    memrq(nt1,sizeof(double));

    if (!(rs = (double *)calloc(nt1,sizeof(double))))  
        goto LTBFin;
    a_rs = 1;
    memrq(nt1,sizeof(double));

    if (!(es = (double *)calloc(nt1,sizeof(double))))  
        goto LTBFin;
    a_es = 1;
    memrq(nt1,sizeof(double));
   
    if (!(ee = (double *)calloc(PMNTP * MaxDes1,sizeof(double))))  
        goto LTBFin;
    a_ee = 1;
    memrq(PMNTP * MaxDes1,sizeof(double));
  
    if (!(subd = (double *)calloc(MaxDes1 + 1,sizeof(double))))  
        goto LTBFin;
    a_subd = 1;
    memrq(MaxDes1 + 1,sizeof(double));
  
    if (!(dens = (double *)calloc(MaxDes1 + 1,sizeof(double))))  
        goto LTBFin;
    a_dens = 1;
    memrq(MaxDes1 + 1,sizeof(double));
  
    if (!(rate = (double *)calloc(MaxDes1 + 1,sizeof(double))))  
        goto LTBFin;
    a_rate = 1;
    memrq(MaxDes1 + 1,sizeof(double));
  
    if (!(df = (char *)calloc(MaxDes1,sizeof(char))))  
        goto LTBFin;
    a_df = 1;
    memrq(MaxDes1,sizeof(char));

    /* loop over all sn,org combinations */

    printf1("SN  Org  Group ");
    if (PM1NV > 0)  
        prnchar(' ',VNameLen + 1,0);
    printf1("  Median  Dest.States  Episodes  Weighted\n");
    if (PM1NV > 0)
        prnchar('-',VNameLen + 1,0);
    prnchar('-',56,1);

    wt = 1.0;

    for (isn = 1; isn <= MaxSnn; ++isn) {      /* loop over all sn */

        for (iorg = 0; iorg <= MaxOrg; ++iorg) {    /* loop over all origins */
       
            ip = TranPtr[(isn - 1) * MaxOrg1 + iorg];
            if (ip < 0)
                continue;

            if (PM1NV) grp = 0; else grp = -1;

            for (; grp < PM1NV; ++grp) {               /* loop for all groups */
    
                n0 = 0;
                ne[0] = 0.0;
                for (l = 0; l < PMNTP; ++l)
                    es[l] = zz[l] = 0.0;
                k = MaxDes1 * PMNTP;
                for (l = 0; l < k; ++l)
                    ee[l] = 0.0;
                for (k = 0; k <= MaxDes; ++k)
                    df[k] = '\0';
    
                /* Loop over all episodes. Calculation of the number of     */
                /* episodes with events in es, differentiated with respect  */
                /* to destination states in ee.                             */

                get_spell(1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
                while (get_spell(0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

                    if (fabs(ts) > EPSI1)
                        err1 = 1;

                    if (sn != isn || org != iorg)               
                        continue;

                    if (grp >= 0) {
                        if (fabs(get_data(PM1VIdx[grp],icase)) <= EPSI1)
                            continue;
                    }

                    if (WIVar >= 0)             /* get weights */
                        wt = get_data(WIVar,icase) * WNorm;
   
                    n0++;
                    ne[0] += wt;
                    tf -= ts;       /* always durations !! */

                    for (l = 1; l < PMNTP; ++l) {
                        if (tf < PMTP[l])  
                            break;
                    }
                    l--;    /* index of time period */

                    if (des == org)     /* censored cases */
                        zz[l] += wt;
                    else {              /* episodes with an event */
                        df[des] = 1;
                        es[l] += wt;
                        ee[l * MaxDes1 + des] += wt;
                    }
                }

                /* Calculate number of destination states in ndes */

                ndes = 0;
                for (k = 0; k <= MaxDes; ++k) {
                    if (df[k])
                        ndes++;
                }
                if (n0 == 0 || ndes == 0)
                    continue;

                /* Calculate number entering the intervals in ne[] */

                for (l = 1; l < PMNTP; ++l)  
                    ne[l] = ne[l - 1] - es[l - 1] - zz[l - 1];
                
                /*  Print headers of requested tables.                      */
   
                prltb1_head(sn,iorg,grp,n0,ne[0],df);
    
                /*  Calculate the risk set and print body of first part of  */
                /*  the life table.                                         */

                for (l = 0; l < PMNTP; ++l) {
                    if (ne[l] <= 0.0)
                        break;
                    rs[l] = ne[l] - PMCFrac * zz[l];
                    if (rs[l] < 0.0)
                        break;
   
                    prltb1(l,ne[l],zz[l],rs[l],ee + l * MaxDes1,df);
                }

                /* Print header of second life table part */
    
                prltb2_head(ndes,df);
   
                /*  Calculate overall survivor function surv, with standard */
                /*  error sterr. Destination-specific densities and rates   */
                /*  calculated in the arrays dens and rate. If only one     */
                /*  destination state, standard errors of the density and   */
                /*  rate function are calculated in derr and rerr.          */

                /*  The median is always calculated by linear interpolation */
                /*  of the overall survivor function.                       */

                surv = 1.0;
                med = sum = 0.0;
                for (k = 0; k <= MaxDes; ++k)
                    subd[k] = 0.0;

                for (l = 0; l < PMNTP; ++l) {

                    if (rs[l] <= 0.0)
                        break;

                    sterr = derr = rerr = -1.0;

                    ql = es[l] / rs[l];
                    pl = 1.0 - ql;

                    /* stand error of survivor function */

                    if (sum >= 0.0)
                        sterr = surv * sqrt(sum);

                    if (l < PMNTP - 1) {

                        surv1 = surv * pl;      /* next survivor function */

                        tt = PMTP[l + 1] - PMTP[l]; /* length of time interval */

                        for (k = 0; k <= MaxDes; ++k) {

                            if (df[k]) {
                                tmp = subd[k] + surv * ee[l * MaxDes1 + k] / rs[l];
                                dens[k] = (tmp - subd[k]) / tt;
                                tmp1 = surv + surv1;
                                if (tmp1 > 0.0)
                                    rate1 = rate[k] = 2.0 * dens[k] / tmp1;
                                else
                                    rate1 = rate[k] = 0.0;
                                subd[k] = tmp;
                            }
                        }

                        /* calculate the median by linear interpolation */
   
                        tmp = surv - surv1;
                        if (surv > 0.5 && surv1 <= 0.5 && tmp > 0.0)
                            med = PMTP[l] + tt * (surv - 0.5) / tmp;
                     
                        /* calculation of stand. errors for density and rate
                           only if there is only one destination state */
     
                        if (ndes == 1 && ql > 0.0) {

                            tmp = pl / (ql * rs[l]);
                            derr = ql * surv * sqrt(sum + tmp) / tt;

                            tmp = rate1 * tt / 2.0;
                            tmp1 = (1.0 - tmp * tmp) / (ql * rs[l]);
                            if (tmp1 > 0.0)
                                rerr = rate1 * sqrt(tmp1);
                        }
                    }
  
                    /* sum is needed for surv function stand. error */

                    if (sum >= 0.0 && pl > 0.0)
                        sum += ql / (pl * rs[l]);
                    else
                        sum = -1.0;
     
                    prltb2(l,surv,sterr,dens,derr,rate,rerr,ndes,df);

                    surv = surv1;
                }       

                printf1("%2d %4d  ",isn,iorg); 
                if (PM1NV == 0)  
                    printf1("  --  ");
                else {
                    printf1("%4d  ",grp + 1);
                    prn_vname(PM1VIdx[grp]);
                }

                if (med > 0.0) {
                    fprintf(PMFd,"\n# Median duration: %4.2f",med);
                    printf1("%8.2f ",med);
                }
                else
                    printf1("     --- ");
                printf1("  %5d   %12d %9.2f\n",ndes,n0,ne[0]);
                ntab++;        
                fprintf(PMFd,"\n\n");

            }   /* end of loop for all groups */
        }       /* end of loop for all origins */
    }           /* end of loop for all serial numbers */

    if (PM1NV > 0)
        prnchar('-',VNameLen + 1,0);
    prnchar('-',56,1);
    printf1("%d table(s) written to: %s\n",ntab,PMFdName);
    err = err1;

LTBFin:
    if (a_ne) {
        free((char *)ne);
        memrq(-nt1,sizeof(double));
    }
    if (a_zz) {
        free((char *)zz);
        memrq(-nt1,sizeof(double));
    }
    if (a_rs) {
        free((char *)rs);
        memrq(-nt1,sizeof(double));
    }
    if (a_es) {
        free((char *)es);
        memrq(-nt1,sizeof(double));
    }
    if (a_ee) {
        free((char *)ee);
        memrq(-PMNTP * MaxDes1,sizeof(double));
    }
    if (a_subd) {
        free((char *)subd);
        memrq(-MaxDes1 - 1,sizeof(double));
    }
    if (a_dens) {
        free((char *)dens);
        memrq(-MaxDes1 - 1,sizeof(double));
    }
    if (a_rate) {
        free((char *)rate);
        memrq(-MaxDes1 - 1,sizeof(double));
    }
    if (a_df) {
        free((char *)df);
        memrq(-MaxDes1,sizeof(char));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prltb1_head.  Print header of first part of life table.                 */

void prltb1_head(int sn,int org,int grp,int n,double w,char *df)
{
    register int k;

    fprintf(PMFd,"# Life table. SN %d. Origin state %d. ",sn,org);
    if (grp >= 0)
        fprintf(PMFd,"\n# Group: %s",VName[PM1VIdx[grp]]);
    fprintf(PMFd,"\n# Cases: %d  weighted: %g\n",n,w);
    fprintf(PMFd,"\n# Start of          Number   Number   Exposed   ");
    for (k = 0; k <= MaxDes; ++k) {
        if (df[k]) {
            fprintf(PMFd,"  D-State %-3d",k);
            fprnchar(PMFd,' ',PMFmt1 - 5,0);
        }
    }
    fprintf(PMFd,"\n# Interval Midpoint Entering Censored   to Risk ");
    for (k = 0; k <= MaxDes; ++k) {
        if (df[k]) {
            fprintf(PMFd,"Events ");
            fprnchar(PMFd,' ',PMFmt1 - 5,0);
            fprintf(PMFd,"Prob  ");
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  prltb1  Print one line of first life table part.                        */

void prltb1(int l,double ne,double zz,double rs,double *ee,char *df)
{
    register int k;
     
    fprintf(PMFd,"\n%10.2f",PMTP[l]);
    if (l < PMNTP - 1)
        fprintf(PMFd," %8.2f",(PMTP[l] + PMTP[l + 1]) / 2.0);
    else
        fprintf(PMFd,"        *");
    fprintf(PMFd," %8d %8d %9.1f ",(int)ne,(int)zz,rs);
    for (k = 0; k <= MaxDes; ++k) {
        if (df[k]) {
            fprintf(PMFd,"%6d ",(int)ee[k]);
            if (rs > 0.0)
                fprintf(PMFd,PMFmtS,ee[k] / rs);
            else
                prstar();
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  prltb2_head().  Print header of second part of life table.              */

void prltb2_head(int ndes,char *df)
{
    register int k;

    fprintf(PMFd,"\n\n# Start of          ");
    fprnchar(PMFd,' ',PMFmt1 - 8,0); fprintf(PMFd,"     Survivor");
    fprnchar(PMFd,' ',PMFmt1 - 3,0);
  
    for (k = 0; k <= MaxDes; ++k) {
        if (df[k]) {
            fprnchar(PMFd,' ',2 * PMFmt1 - 13,0);
            fprintf(PMFd,"  D-State %-3d  ",k);
            if (ndes == 1) {
                fprnchar(PMFd,' ',2 * PMFmt1 - 13,0);
                fprintf(PMFd,"  D-State %-3d  ",k);
            }
        }
    }
    fprintf(PMFd,"\n# Interval Midpoint ");
    fprnchar(PMFd,' ',PMFmt1 - 7,0); fprintf(PMFd,"Function");
    fprnchar(PMFd,' ',PMFmt1 - 5,0); fprintf(PMFd,"Error ");
    for (k = 0; k <= MaxDes; ++k) {
        if (df[k]) {
            fprnchar(PMFd,' ',PMFmt1 - 7,0); fprintf(PMFd,"Density ");
            fprnchar(PMFd,' ',PMFmt1 - 5,0);
            if (ndes == 1) {
                fprintf(PMFd,"Error "); fprnchar(PMFd,' ',PMFmt1 - 7,0); 
                fprintf(PMFd,"   Rate "); fprnchar(PMFd,' ',PMFmt1 - 5,0);
                fprintf(PMFd,"Error");
            }
            else
                fprintf(PMFd," Rate ");
        }
    }
}    

/* ------------------------------------------------------------------------ */
/*  prltb2  Print one line of second life table part.                       */

void prltb2 (int l, double surv, double sterr, double *dens, double derr,
    double *rate, double rerr, int ndes, char *df)
{
    register int k;

    fprintf(PMFd,"\n%10.2f",PMTP[l]);
    if (l < PMNTP - 1)
        fprintf(PMFd," %8.2f ",(PMTP[l] + PMTP[l + 1]) / 2.0);
    else
        fprintf(PMFd,"        * ");

    fprintf(PMFd,PMFmtS,surv);
    if (sterr >= 0.0)
        fprintf(PMFd,PMFmtS,sterr);
    else  
        prstar();

    for (k = 0; k <= MaxDes; ++k) {
        if (df[k]) {
            if (dens[k] >= 0.0)
                fprintf(PMFd,PMFmtS,dens[k]);
            else
                prstar();
  
            if (ndes == 1) {
                if (derr > 0.0)
                    fprintf(PMFd,PMFmtS,derr);
                else
                    prstar();

                if (rate[k] >= 0.0)
                    fprintf(PMFd,PMFmtS,rate[k]);
                else
                    prstar();

                if (rerr > 0.0)
                    fprintf(PMFd,PMFmtS,rerr);
                else
                    prstar();
            }
            else {
                if (rate[k] >= 0.0)
                    fprintf(PMFd,PMFmtS,rate[k]);
                else
                    prstar();
            }
        }
    }
}

void prstar(void)
{
    fprnchar(PMFd,' ',PMFmt1 - 2,0);
    fprintf(PMFd," * ");
}



