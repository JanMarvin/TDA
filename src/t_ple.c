/****************************************************************************/
/*  t_ple                                                                   */
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
#include "t_edat.h"
#include "t_gdat.h"
#include "t_var.h"
#include "t_gf.h"
#include "t_cdf.h"
#include "t_lsei.h"
#include "t_parm.h"
#include "t_pgen.h"
#include "t_alloc.h"
#include "t_ml.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_ple.c                                                    */

int ple(void); 
int plest(int qtyp);
void prple_head(int sn, int org, int des, int grp);
void prple(int id, int idx, double time, double evnts,
    double cens, double rs, double surv, double sterr);
void prstar1(void);
void prquant_head(void);
void prquant(int qtyp,int sn, int org, int des, int grp);
int csf_stat(void);
int plestat(void);
int dple(void);
int dltb(void);
int diple(void);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

int *RStat;         /* degrees of freedom and                               */  
double *TStat;      /* test statistics for comparing survivor functions     */

/* ------------------------------------------------------------------------ */
/*  ple()           Product-limit estimation. ple command in CmdBuf.        */
/*                  Return 0 if OK, otherwise -1.                           */

int ple(void)
{
    register int i;
    int err,qtyp;

    err = -1;         
    if (check_cmd(1))
        return(-1);

    printf1("Product-limit estimation. Current memory: %d bytes.\n",MemReq);

    if (EDAvail == 0) {
        p_err(-15,1);
        return(0);
    }
    if (parm(CmdBuf + 3,1,0))     /* get parameters */
        goto PLEFin;

    if (PMFmtF == 0)            /* default print format */
        pmfmt(7,5);

    prn_cwt();
    newline();

    if (ESortFlg == 0) {        /* sort data */
        if (sort_ed(1))      
            goto PLEFin;
    }
    qtyp = 0;
    if (PMQOFlg || PMQTFlg) {           /* quantiles or orders */
        if (alloc_acx(PMNTP + 1))       /* allocate AcX */
            goto PLEFin;

        if (PMQOFlg) {
            qtyp = 1;
            for (i = 0; i < PMNTP; ++i) {
                if (PMTP[i] >= 1.0 || PMTP[i] <= 0.0 ||
                       (i > 0 && PMTP[i] >= PMTP[i - 1])) {
                    printf1("Error in qo parameter.\n");
                    goto PLEFin;
                }
            }
        }
        else {
            qtyp = 2;
            for (i = 0; i < PMNTP; ++i) {
                if (PMTP[i] <= 0.0 || 
                       (i > 0 && PMTP[i] < PMTP[i - 1])) {
                    printf1("Error in qt parameter.");
                    goto PLEFin;
                }
            }
        }
    }
    else
        free_tp();      /* set PMNTP = 0 */

    if (PMFDef || qtyp) {
        printf1("\nProduct-limit estimation.\n");
        plest(qtyp);
    }
    if (PMCSF) {    
        printf1("\nComparing survivor functions.\n");
        if (PM1NV < 2)  
            p_err(-31,1);
        else
            csf_stat();
    }
    err = 0;

PLEFin:
    if (ESortFlg)
        sort_ed(0);     /* free memory */

    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  plest(qtyp1,qytp2)  Product-limit estimation.                           */
/*                                                                          */
/*                      If PMFDef print tables to PMFd.                     */
/*                      If qtyp=1 calcualte quantiles in AcX[], orders are  */
/*                      given in PMTP[] (0,...,PMNTP-1).                    */
/*                      If qtyp=2 calcualte orders in AcX[], quantiles are  */
/*                      given in PMTP[] (0,...,PMNTP-1).                    */
/*                                                                          */
/*                      Return 0 if OK, -1 if error.                        */

int plest(int qtyp)
{
    register int grp,icase,jcase,j;
    int isn,iorg,ides,sn;
    int org,des,des1,tran,id,idx,io,jo,nn,evflg,nrpar,ntab;
    double q,tmp,surv,surv1,evnts,cens,cens1,cens2,rs,ws,wt,wt1,ts,tf,time,time1;
    double med,stsum,sterr,ord,ts1,tf1;

    if (PMFDef == 0 && qtyp == 0)
        return(0);

    if (qtyp)
        prquant_head();     /* print header for quantile table */

    id = -1;    /* id for tables */
    nrpar = 0;  /* total number of parameters for plestr() */
    ntab = 0;   /* number of tables in output file */

    /* loop over all sn,org combinations */

    wt1 = wt = 1.0;

    for (tran = 0; tran < NTran1; ++tran) { /* loop over all transitions */

        isn  = SnTran1[tran];
        iorg = OrgTran1[tran];
        ides = DesTran1[tran];

        if (PM1NV) grp = 0; else grp = -1;

        for (; grp < PM1NV; ++grp) {               /* loop for all groups */

            /* print header for ple table */
                
            if (PMFDef) {
                prple_head(isn,iorg,ides,grp);
                ntab++;
            }
  
            /* Calculate the risk set for time 0 in rs. Take all        */
            /* episodes with origin state org (and contained in group   */
            /* grp) with starting times zero. The number of episodes is */
            /* calculated in nn, weighted in ws. It is assumed that     */
            /* episodes are sorted according to their starting times.   */

            rs = ws = 0.0;
            jcase = nn = 0;

            while (++jcase <= NOC) {

                jo = TSIdx[jcase - 1];
           
                get_edat(jo,&sn,&org,&des,&ts,&tf);

                if (sn == isn) {

                    if (ts > 0.0) {
                        jcase--;
                        break;
                    }
                    if (org != iorg)   
                        continue;

                    if (grp >= 0) {
                        if (fabs(get_data(PM1VIdx[grp],jo)) <= EPSI1)
                            continue;
                    }
                    if (WIVar >= 0)             /* get weights */
                        wt1 = get_data(WIVar,jo) * WNorm;
   
                    rs += wt1;
                    ws += wt1;
                    nn++;
                }
            }   
            id++;
            idx = 0;
            time = time1 = 0.0;     
            surv  = 1.0;            
            evnts = cens = cens1 = cens2 = 0.0;
            med  = 0.0;                
            stsum = 0.0;
        
            for (j = 0; j < PMNTP; ++j)  
                AcX[j] = 0.0;
           
            evflg = j = 0;

            /* Now start the main loop for all episodes with origin     */
            /* state (and contained in grp), according to their ending  */
            /* times.                                                   */

            for (icase = 0; icase < NOC; ++icase) {

                io = TFIdx[icase];

                get_edat(io,&sn,&org,&des1,&ts,&tf);

                if (sn != isn || org != iorg)
                    continue;

                if (grp >= 0) {
                    if (fabs(get_data(PM1VIdx[grp],io)) <= EPSI1)
                        continue;
                }

                if (WIVar >= 0)             /* get weights */
                    wt = get_data(WIVar,io) * WNorm;

                if (tf > time) {    

                    if (evnts || time == 0.0) {

                        /* Update the risk set. Add all episodes    */
                        /* with starting time less than time.       */

                        while (jcase < NOC && ++jcase <= NOC) {

                            jo = TSIdx[jcase - 1];
                            get_edat(jo,&sn,&org,&des,&ts1,&tf1);
        
                            if (sn == isn) {

                                if (ts1 >= time) {
                                    jcase--;
                                    break;
                                }
                                if (org != iorg)   
                                    continue;

                                if (grp >= 0) {
                                    if (fabs(get_data(PM1VIdx[grp],jo)) <= EPSI1)
                                        continue;
                                }
                                if (WIVar >= 0)             /* get weights */
                                    wt1 = get_data(WIVar,jo) * WNorm;

                                rs += wt1;
                                ws += wt1;
                                nn++;
                            }
                        }             
                        if (rs <= 0.0 && jcase >= NOC)
                            break;

                        /* calculate actual value of survivor       */
                        /* function in surv1.                       */

                        if (rs > 0.0)
                            surv1 = surv * (1.0 - evnts / rs);
                        else
                            surv1 = surv;

                        /* calculate the median */

                        if (surv1 <= 0.5 && surv > 0.5)  
                            med = time1 + (surv - 0.5) * (time - time1) /
                                                             (surv - surv1);

                        /* If quantils are requested, then:         */
                        /* if qtyp == 1, calculate quantiles.       */
                        /* if qtyp == 2, calculate orders.          */
                         
                        while (j < PMNTP) {
                            if (qtyp == 1) {
                                ord = PMTP[j];
                                if (surv1 <= ord && surv > ord) {
                                        tmp = time1 + (surv - ord) * 
                                            (time - time1) / (surv - surv1);
                                        AcX[j++] = tmp;
                                }
                                else
                                    break;
                            }
                            else if (qtyp == 2) {
                                q = PMTP[j];
                                if (surv1 < surv &&
                                              q > time1 && q <= time) {
                                    tmp = surv - (q - time1) * 
                                        (surv - surv1) / (time - time1);
                                    AcX[j++] = tmp;
                                }
                                else
                                    break;
                            }
                        }
                 
                        /* update value of the survivor function */

                        surv = surv1;

                        /* calculate standard errors in sterr */

                        sterr = -1.0;
                        if (rs > evnts) { 
                            stsum += evnts / (rs * (rs - evnts));
                            sterr = surv * sqrt(stsum);
                        }

                        /* now print one line of the ple table */
   
                        if (PMFDef)
                            prple(id,idx,time,evnts,cens1,rs,surv,sterr);
           
                        idx++;
                        cens1 = cens2;
                        cens2 = 0.0;
                        evflg = 0;

                        time1 = time;  /* save actual time in time1 */
                    }
                    rs -= (evnts + cens);    /* update the risk set */
                    evnts = cens = 0.0;      /* and actual time.    */
                    time = tf;
                }
                if (des1 != ides) {
                    cens += wt;  
                    if (!evflg)
                        cens1 += wt;
                    else
                        cens2 += wt;
                }
                else {           
                    evnts += wt;
                    evflg = 1;
                }
            }
        
            /* one more line if events or censored episodes are left */

            cens1 += cens2;

            if (evnts || cens1) {
                sterr = -1.0;
                if (rs > 0.0) {
                    surv *= (1.0 - evnts / rs);
                    if (surv > 0.0) {
                        stsum += evnts / (rs * (rs - evnts));
                        sterr = surv * sqrt(stsum);
                    }
                }
                else
                    surv = -1.0;
  
                if (PMFDef)
                    prple(id,idx,time,evnts,cens1,rs,surv,sterr);
            }
       
            /* print the median, and number of episodes */
       
            if (PMFDef) {
                if (med > 0.0) {
                    fprintf(PMFd,"\n# Median Duration: %4.2f",med);
                    if (cens1)
                        fprintf(PMFd,"\n# Duration times limited to: %g",time1);
                }
                fprintf(PMFd,"\n# Cases: %d  weighted: %g\n\n",nn,ws);
            }

            /* print table of quantiles and count parameters */
  
            prquant(qtyp,isn,iorg,ides,grp);
            nrpar += PMNTP;         
  
        }   /* end of group selection */
    }       /* end of loop for transitions */

    if (PMFDef) {
        if (qtyp)  
            printf1("\n");           
        printf1("%d table(s) written to: %s\n",ntab,PMFdName);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prple_head.  Print header of product-limit estimation table             */

void prple_head(int sn, int org, int des, int grp)
{
    fprintf(PMFd,"# SN %d. Transition: %d,%d - Product-Limit Estimation\n",
                                                              sn,org,des);
    if (grp >= 0)
        fprintf(PMFd,"\n# Group: %s",VName[PM1VIdx[grp]]);
    fprintf(PMFd,"\n#                      Number  Number    Exposed ");
    fprnchar(PMFd,' ',PMFmt1 - 7,0); fprintf(PMFd," Survivor");
    fprnchar(PMFd,' ',PMFmt1 - 5,0); fprintf(PMFd," Std. ");
    fprnchar(PMFd,' ',PMFmt1 - 5,0); fprintf(PMFd," Cum. ");
    fprintf(PMFd,"\n# ID Index      Time   Events  Censored  to Risk ");
    fprnchar(PMFd,' ',PMFmt1 - 7,0); fprintf(PMFd," Function");
    fprnchar(PMFd,' ',PMFmt1 - 5,0); fprintf(PMFd,"Error ");
    fprnchar(PMFd,' ',PMFmt1 - 5,0); fprintf(PMFd," Rate ");
}

/* ------------------------------------------------------------------------ */
/*  prple.  Print one line of the product-limit estimation table.           */

void prple(int id, int idx, double time, double evnts,
    double cens, double rs, double surv, double sterr) 
{
    double tmp;

    if (evnts < 0.0 || rs < 0.0 || surv < 0.0 || cens < 0.0)
        return;

    if (evnts > 0 || surv == 1.0)
        fprintf(PMFd,"\n%4d %5d",id,idx);
    else
        fprintf(PMFd,"\n# %2d %5d",id,idx);
    fprintf(PMFd," %9.2f",time);
    fprintf(PMFd," %8d %8d",(int)evnts,(int)cens);
    if (evnts > 0 || surv == 1.0) {
        fprintf(PMFd," %9d  ",(int)rs);
        fprintf(PMFd,PMFmtS,surv);
        if (sterr >= 0.0)
            fprintf(PMFd,PMFmtS,sterr);
        else
            prstar1();
        if (surv > 0.0) {
            if (surv == 1.0)
                tmp = 0.0;
            else
                tmp = -rlog(surv);
            fprintf(PMFd,PMFmtS,tmp);
        }
        else  
            prstar1();
    }
}

void prstar1(void)
{
    fprnchar(PMFd,' ',PMFmt1 - 2,0);
    fprintf(PMFd," * ");
}

/* ------------------------------------------------------------------------ */
/*  prquant_head.  Print header of table for quantiles.                     */

void prquant_head(void)
{
    if (PMNTP) {
        printf1("\n");
        prnchar(' ',14,0);
        if (PM1NV > 0)   
            prnchar(' ',VNameLen + 7,0);
  
        printf1("  Survivor      Time\n");
        printf1("SN  Org  Des  ");
  
        if (PM1NV > 0) {            
            printf1("Group  ");
            prnchar(' ',VNameLen,0);
        }
        printf1("  Function    Quantile\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  prquant.  Print table of quantiles.                                     */

void prquant(int qtyp,int sn, int org, int des, int grp)
{
    register int j;

    if (PMNTP <= 0)
        return;

    prnchar('-',36,0);
  
    if (PM1NV > 0)   
        prnchar('-',VNameLen + 7,0);

    for (j = 0; j < PMNTP; ++j) {

        printf1("\n");           

        if (!j) {
            printf1("%2d %4d %4d",sn,org,des);
            if (PM1NV > 0) {
                printf1("%4d  ",grp + 1);
                prn_vname(PM1VIdx[grp]);
            }
        }
        else {
            prnchar(' ',12,0);
            if (PM1NV > 0)   
                prnchar(' ',VNameLen + 7,0);
        }
        if (AcX[j] <= 0.0 || PMTP[j] <= 0.0)  
            break;
            /*  printf1("        ***          ***");       */
        else if (qtyp == 1)
            printf1(" %10.4f %12.4f",PMTP[j],AcX[j]);
        else
            printf1(" %10.4f %12.4f",AcX[j],PMTP[j]);
    }
    printf1("\n");
}

/* ------------------------------------------------------------------------ */
/*  csf_stat   Test statistics for comparing survivor functions.            */
/*                                                                          */

#define PLENTest 4      /* number of test statistics */

int csf_stat(void)
{
    register int j,tran;
    int err,n,nn,sn,org,des,tsa,rsa;
    double tmp;

    nn = NTran1 * PLENTest;
    tsa = rsa = err = 0;

    if (!(TStat = (double *)calloc(nn,sizeof(double)))) { 
        err = -2;
        goto CSFFin;
    }
    tsa = 1;
    memrq(nn,sizeof(double));

    if (!(RStat = (int *)calloc(nn,sizeof(int)))) { 
        err = -2;
        goto CSFFin;
    }
    rsa = 1;
    memrq(nn,sizeof(int));

    if (plestat())      /* perform calculations */
        goto CSFFin;

    printf1("\nSN  Org Des   Test Statistic                 T-Stat   DF   Signif");

    for (tran = 0; tran < NTran1; ++tran) { /* loop over all transitions */

        sn  = SnTran1[tran];
        org = OrgTran1[tran];
        des = DesTran1[tran];

        printf1("\n");
        prnchar('-',65,0);

        for (j = 0; j < PLENTest; ++j) {
            printf1("\n%2d  %3d %3d    ",sn,org,des);
            if (!j)  
                printf1("Log-Rank (Savage)     ");
            else {
                printf1("Wilcoxon (");
                if (j == 1)  
                    printf1("Breslow)    ");          
                else if (j == 2)  
                    printf1("Tarone-Ware)");          
                else if (j == 3)  
                    printf1("Prentice)   ");          
            }
            tmp = TStat[tran * PLENTest + j];
            printf1("%14.4f ",tmp);
            n = RStat[tran * PLENTest + j];
            printf1("%4d ",n);
            if (n > 0 && tmp > 0.0)  
                printf1("%8.4f",cdchif(tmp,n));
            else
                printf1("       *");
        }
    }
    printf1("\n");         

CSFFin:
    if (err == -2)  
        p_err(-2,1);
       
    if (tsa) {
        free((char *)TStat);
        memrq(-nn,sizeof(double));
    }
    if (rsa) {
        free((char *)RStat);
        memrq(-nn,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  plestat     Calculation of test statistics for comparison of            */
/*              survivor functions.                                         */
/*              Return 0 if OK, -1 if error (insuff. memory)                */

int plestat(void)
{
    register int l,i,j,icase,jcase;
    int err,rsa,eva,cena,va[PLENTest],ua[PLENTest];
    int sn,isn,iorg,ides,ii,nt,nn,g,org,des,des1,tran,io,jo;
    double tmp,tmp1,tmp2,wt,wt1,ts,tf,ts1,tf1,time,r0,e0,w[PLENTest];
    double *rs,*ev,*cen,*v[PLENTest],*u[PLENTest];

    rsa = eva = cena = err = 0;
    nt = PLENTest;     
    nn = PM1NV - 1;

    if (!(rs = (double *)calloc(PM1NV + 1,sizeof(double)))) {
        err = -1;
        goto PLESTATFin;
    }
    rsa = 1;
    memrq(PM1NV + 1,sizeof(double));

    if (!(ev = (double *)calloc(PM1NV + 1,sizeof(double)))) {
        err = -1;
        goto PLESTATFin;
    }
    eva = 1;
    memrq(PM1NV + 1,sizeof(double));

    if (!(cen = (double *)calloc(PM1NV + 1,sizeof(double)))) {
        err = -1;
        goto PLESTATFin;
    }
    cena = 1;
    memrq(PM1NV + 1,sizeof(double));

    for (i = 0; i < PLENTest; ++i) {
        if (!(v[i] = (double *)calloc(nn * nn + 1,sizeof(double)))) {
            err = -1;
            goto PLESTATFin;
        }
        va[i] = 1;
        memrq(nn * nn + 1,sizeof(double));

        if (!(u[i] = (double *)calloc(nn + 1,sizeof(double)))) {
            err = -1;
            goto PLESTATFin;
        }
        ua[i] = 1;
        memrq(nn + 1,sizeof(double));
    }
     
    wt = wt1 = 1.0;             

    for (tran = 0; tran < NTran1; ++tran) { /* loop over all transitions */

        isn  = SnTran1[tran];
        iorg = OrgTran1[tran];
        ides = DesTran1[tran];

        for (i = 0; i < PLENTest; ++i) {
            w[i] = 1.0;
            l = 1;
            for (j = 1; j <= nn; ++j) {
                u[i][j] = 0.0;
                for (l = 1; l <= nn; ++l)  
                    v[i][(j - 1) * nn + l] = 0.0;
            }
        }

        /* build risk set for time = 0 */

        for (i = 0; i <= PM1NV; ++i)
            rs[i] = 0.0;

        jcase = 0;
        while (++jcase <= NOC) {

            jo = TSIdx[jcase - 1];

            get_edat(jo,&sn,&org,&des,&ts,&tf);

            if (sn == isn) {

                if (ts > 0.0) {
                    jcase--;
                    break;
                }
                if (org != iorg)   
                    continue;
  
                for (l = 0; l < PM1NV; ++l) {
                    if (fabs(get_data(PM1VIdx[l],jo)) > EPSI1)
                        break;
                }
                if (l < PM1NV) {

                    if (WIVar >= 0)             /* get weights */
                        wt1 = get_data(WIVar,jo) * WNorm;

                    rs[0] += wt1;
                    rs[l + 1] += wt1;
                }
            }   
        }
        time = 0.0; 
        for (i = 0; i <= PM1NV; ++i)
            ev[i] = cen[i] = 0.0;

        for (icase = 0; icase < NOC; ++icase) {

            io = TFIdx[icase];

            get_edat(io,&sn,&org,&des1,&ts,&tf);

            if (sn != isn || org != iorg)
                continue;

            for (l = 0; l < PM1NV; ++l) {
                if (fabs(get_data(PM1VIdx[l],io)) > EPSI1)
                    break;
            }
            if (l >= PM1NV)  
                continue;

            g = l + 1;              /* group of episode */

            if (WIVar >= 0)             /* get weights */
                wt = get_data(WIVar,io) * WNorm;

            if (tf > time) {    

                if (ev[0]) {        /* update the risk set */

                    while (jcase < NOC && ++jcase <= NOC) {

                        jo = TSIdx[jcase - 1];
                        get_edat(jo,&sn,&org,&des,&ts1,&tf1);
        
                        if (sn == isn) {

                            if (ts1 >= time) {
                                jcase--;
                                break;
                            }
                            if (org != iorg)   
                                continue;

                            for (l = 0; l < PM1NV; ++l) {
                                if (fabs(get_data(PM1VIdx[l],jo)) > EPSI1)
                                    break;
                            }
                            if (l < PM1NV) {

                                if (WIVar >= 0)             /* get weights */
                                    wt1 = get_data(WIVar,jo) * WNorm;

                                rs[0] += wt1;
                                rs[l + 1] += wt1;
                            }
                        }   
                    }
                                /* update U and V */

                    r0 = rs[0];
                    e0 = ev[0];
                    w[1] = r0 / WSum;                 
                    w[2] = sqrt(w[1]);
                    w[3] *= (r0 - e0 + 1.0) / (r0 + 1.0);
                    tmp = e0 * (r0 - e0) / (r0 * (r0 - 1.0));

                    for (i = 1; i <= nn; ++i) {
                        tmp1 = ev[i] - e0 * rs[i] / r0;
                        for (l = 0; l < PLENTest; ++l)  
                            u[l][i] += tmp1 * w[l];

                        tmp1 = tmp * rs[i];
                        for (j = 1; j <= nn; ++j) {
                            if (i == j)
                                tmp2 = tmp1 * (1.0 - rs[j] / r0);
                            else
                                tmp2 = tmp1 * (-rs[j] / r0);

                            for (l = 0; l < PLENTest; ++l)
                                v[l][(i - 1) * nn + j] += tmp2 * w[l] * w[l];
                        }
                    }
                    if (rs[0] <= 0.0 && jcase >= NOC)
                        break;
                }
                for (i = 0; i <= PM1NV; ++i) {
                    rs[i] -= (ev[i] + cen[i]);
                    ev[i] = cen[i] = 0.0;
                }
                time = tf;
            }
            if (des1 != ides) {    /* Censored Cases */
                cen[0] += wt;
                cen[g] += wt;
            }
            else {
                ev[0] += wt;
                ev[g] += wt;
            }
        }
        if (ev[0] > 0.0 && rs[0] > 1.0) {
                                                       /* update U and V */
            r0 = rs[0];
            e0 = ev[0];
            w[1] = r0 / WSum;   
            w[2] = sqrt(w[1]);
            w[3] *= (r0 - e0 + 1.0) / (r0 + 1.0);
            tmp = e0 * (r0 - e0) / (r0 * (r0 - 1.0));

            for (i = 1; i <= nn; ++i) {
                tmp1 = ev[i] - e0 * rs[i] / r0;
                for (l = 0; l < PLENTest; ++l)  
                    u[l][i] += tmp1 * w[l];

                tmp1 = tmp * rs[i];
                for (j = 1; j <= nn; ++j) {
                    if (i == j)
                        tmp2 = tmp1 * (1.0 - rs[j] / r0);
                    else
                        tmp2 = tmp1 * (-rs[j] / r0);

                    for (l = 0; l < PLENTest; ++l)
                        v[l][(i - 1) * nn + j] += tmp2 * w[l] * w[l];
                }
            }
        }
        if (PMProtFDef) {
            fprintf(PMProtFd,"SN %d Transition: %d,%d\n",isn,iorg,ides);
            for (l = 0; l < PLENTest; ++l) {
                fprintf(PMProtFd,"\nTest statistic %d",l + 1);
                fprintf(PMProtFd,"\n U: ");
                for (i = 1; i <= nn; ++i)
                    fprintf(PMProtFd,PMPFmtS,u[l][i]);
    
                for (i = 1; i <= nn; ++i) {
                    fprintf(PMProtFd,"\n V: ");
                    for (j = 1; j <= nn; ++j)
                        fprintf(PMProtFd,PMPFmtS,v[l][(i - 1) * nn + j]);
                }
                fprintf(PMProtFd,"\n");
            }
            fprintf(PMProtFd,"\n");
        }

        ii = tran * PLENTest;
        for (l = 0; l < PLENTest; ++l) {
            RStat[ii + l] = syminv(nn,v[l]);
            TStat[ii + l] = 0.0;
            for (i = 1; i <= nn; ++i) {
                tmp = 0.0;
                for (j = 1; j <= nn; ++j)  
                    tmp += v[l][(i - 1) * nn + j] * u[l][j];
                TStat[ii + l] += tmp * u[l][i];
            }
        }
    }       /* end of tran loop */

PLESTATFin:
    if (err == -1)  
        p_err(-2,1);
    if (rsa) {
        free((char *)rs);
        memrq(-PM1NV - 1,sizeof(double));
    }
    if (eva) {
        free((char *)ev);
        memrq(-PM1NV - 1,sizeof(double));
    }
    if (cena) {
        free((char *)cen);
        memrq(-PM1NV - 1,sizeof(double));
    }
    for (i = 0; i < PLENTest; ++i) {
        if (va[i]) {
            free((char *)v[i]);
            memrq(-nn * nn - 1,sizeof(double));
        }
        if (ua[i]) {
            free((char *)u[i]);
            memrq(-nn - 1,sizeof(double));
        }
    }
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  dple()      Product-limit estimation with discrete times.               */
/*                                                                          */  
/*              dple(                                                       */
/*                  df=...,         output file (required)                  */
/*                  fmt=...,        print format, def. 12.4                 */
/*                  opt=...,        output option, def. 1                   */
/*                                  1 : all time points                     */ 
/*                                  2 : only time points occurring in T     */
/*                  grp=...,        list of group variables                 */
/*              ) = [S,] T, D;                                              */
/*                                                                          */
/*  There must be two or three variable names on the right-hand side        */
/*  defining durations on a process time axis beginning at zero. S is the   */
/*  beginning of the observation period (zero if not specified), and T is   */
/*  the duration, if D is not equal to zero, and otherwise is the end of    */
/*  the observation period. It is assumed that values of S and T are        */
/*  integers, furthermore: 0 <= S <= T.                                     */
/*                                                                          */
/*  For each time point t (min(S) <= t <= max(T)), the command counts:      */  
/*      e(t) = number of events at t, that is, the number of cases i with   */
/*             T(i) = t and D(i) != 0; and                                  */
/*      r(t) = number of elements in the risk set at t, that is, the number */
/*             of cases i with S(i) <= t and T(i) >= t                      */
/*                                                                          */  
/*  Finally, for each time point t, one record is written to the output     */
/*  file specified by the df parameter, containing the following columns:   */
/*      (1)  index of table                                                 */
/*      (2)  time                                                           */  
/*      (3)  r(t)                                                           */
/*      (4)  e(t)                                                           */
/*      (5)  number of censored cases                                       */
/*      (6)  the rate function e(t) / r(t)                                  */
/*      (7)  the survivor function calculated as the product of             */
/*           (1 - e(j)/r(j)) for j < t.                                     */
/*                                                                          */
/*  If the grp parameter is used to specify groups, the calculations are    */  
/*  done separately for each group. If opt = 1, the output tables contain   */
/*  one row for each time point beginning at min(S) and ending at max(T).   */
/*  If opt = 2, only rows with at least one observation (value of T) are    */  
/*  printed.                                                                */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int dple(void)
{
    register int i,j;
    int err,grp,idx,iis,iif,iid,t,ta,tmin,tmax,nt,nrec;
    double r,g,m,rr,gg,w,n,wt,ncen,nev;

    err = -1;         
    if (check_cmd(1))
        return(-1);

    printf1("PLE with discrete times. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,4,1))     /* get parameters */
        goto DPLEFin;
   
    if (PMF1Def == 0) {
        printf1("Error: need an output file.\n");
        goto DPLEFin;
    }
    if (PMFmtF == 0)
        pmfmt(12,4);

    if (PMNV == 2) {
        iis = -1;         
        iif = PMVIdx[0];
        iid = PMVIdx[1];
    }
    else if (PMNV == 3) {
        iis = PMVIdx[0];
        iif = PMVIdx[1];
        iid = PMVIdx[2];
    }
    else {
        printf1("Error: need two or three variables on right-hand side.\n");
        goto DPLEFin;
    }
    newline();
    nrec = 0;

    if (PM1NV) {
        grp = 0;
        printf1("Group  defined by      ");
        prnchar(' ',PMFmt1 - 5,0);
        printf1("Cases ");
        printf1(" Minimum  Maximum\n");
    }
    else  
        grp = -1;
    idx = -1;

    for (; grp < PM1NV; ++grp) {               /* loop for all groups */
        idx++;
        tmax = 0;                        
        if (iis == -1)
            tmin = 0;
        else
            tmin = INTMAX;

        wt = 1.0;
        n = 0.0;
        for (i = 0; i < NOC; ++i) {

            if (grp >= 0) {
                if (fabs(get_data(PM1VIdx[grp],i)) <= EPSI1)
                    continue;
            }
            if (WIVar >= 0) {                           /* get weights */
                wt = get_data(WIVar,i) * WNorm;
                n += wt;
            }
            else
                n += 1.0;

            ta = 0;
            if (iis >= 0) {
                ta = (int)get_data(iis,i);
                tmin = imin(ta,tmin);
                tmax = imax(ta,tmax);
            }
            t = (int)get_data(iif,i);
            tmin = imin(t,tmin);
            tmax = imax(t,tmax);
           
            if (ta < 0 || t < 0 || t < ta) {
                printf1("Error: found inconsistent data in case %d.\n",i + 1);
                goto DPLEFin;
            }
        }
        if (grp >= 0) {
            printf1("%5d  %s",grp,VName[PM1VIdx[grp]]);
            prnchar(' ',16 - strlen(VName[PM1VIdx[grp]]),0);
            printf1(PMFmtS,n);
            if (n > 0.0)
                printf1("%8d %8d",tmin,tmax);
            printf1("\n");
        }
        else {
            printf1("Number of cases (weighted): %g\n",n);
            printf1("Minimal time value: %d\n",tmin);
            printf1("Maximal time value: %d\n",tmax);
        }
        if (n <= EPSI1)
            continue;

        nt = tmax - tmin + 1;

        if (alloc_acx(nt + 1))
            goto DPLEFin;
        if (alloc_acy(nt + 1))
            goto DPLEFin;
        if (alloc_acz(nt + 1))
            goto DPLEFin;

        wt = 1.0;
        for (i = 0; i < NOC; ++i) {

            if (grp >= 0) {
                if (fabs(get_data(PM1VIdx[grp],i)) <= EPSI1)
                    continue;
            }
            if (WIVar >= 0)                             /* get weights */
                wt = get_data(WIVar,i) * WNorm;

            t = (int)get_data(iif,i);
            if (get_data(iid,i) != 0.0)
                AcX[t - tmin] += wt;
            else
                AcY[t - tmin] += wt;

            if (iis >= 0)
                ta = (int)get_data(iis,i);
            else
                ta = 0;

            for (j = ta; j <= t; ++j) 
                AcZ[j - tmin] += wt;

        }
        if (grp >= 0) {
            if (idx > 0) {
                fprintf(PMF1d,"\n");
                nrec++;
            }
            fprintf(PMF1d,"# Group %d defined by: %s\n",grp,VName[PM1VIdx[grp]]);
            nrec++;       
        }
        fprintf(PMF1d,"#");
        fprnchar(PMF1d,' ',11 + 5 * (PMFmt1 + 1) - 8,0);
        fprintf(PMF1d,"survivor\n");
        fprintf(PMF1d,"# Idx   Time ");
        fprnchar(PMF1d,' ',PMFmt1 - 7,0);
        fprintf(PMF1d,"at risk ");
        fprnchar(PMF1d,' ',PMFmt1 - 6,0);
        fprintf(PMF1d,"events ");
        fprnchar(PMF1d,' ',PMFmt1 - 8,0);
        fprintf(PMF1d,"censored ");
        fprnchar(PMF1d,' ',PMFmt1 - 4,0);
        fprintf(PMF1d,"rate ");
        fprnchar(PMF1d,' ',PMFmt1 - 8,0);
        fprintf(PMF1d,"function\n");
        nrec += 2;

        g = 1.0;
        ncen = nev = 0.0;
        for (i = 0; i < nt; ++i) {
            if (PMOPT == 2) {
                if (AcX[i] == 0.0 && AcY[i] == 0.0)
                    continue;
            }
            nev += AcX[i];
            ncen += AcY[i];

            fprintf(PMF1d,"%5d ",idx);
            fprintf(PMF1d,"%6d ",tmin + i);
            fprintf(PMF1d,PMFmtS,AcZ[i]);
            fprintf(PMF1d,PMFmtS,AcX[i]);
            fprintf(PMF1d,PMFmtS,AcY[i]);

            r = 0.0;
            if (AcZ[i] > 0.0)  
                r = AcX[i] / AcZ[i];

            fprintf(PMF1d,PMFmtS,r);
            fprintf(PMF1d,PMFmtS,g);

            if (PMOPT == 3) {           /* calculate mean values */

                gg = 1.0;
                w = m = 0.0;
                for (j = i; j < nt; ++j) {

                    rr = 0.0;
                    if (AcZ[j] > 0.0)  
                        rr = AcX[j] / AcZ[j];
                    if (j == nt - 1)
                        rr = 1.0;
                    m += (double)(tmin + j) * gg * rr;
                    w += gg * rr;

                    gg *= (1.0 - rr);
                }
                fprintf(PMF1d,"%12.4f ",m);
                fprintf(PMF1d,"%12.4f ",m - (double)(tmin + i));
            }
            fprintf(PMF1d,"\n");
            nrec++;
            g *= (1.0 - r);
        }
        fprintf(PMF1d,"# Number of events: %g\n",nev);
        fprintf(PMF1d,"# Number of censored cases: %g\n",ncen);
    }
    printf1("\n%d records written to: %s\n",nrec,PMF1dName);
    err = 0;

DPLEFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  dltb()      Life-table with discrete times.                             */
/*                                                                          */  
/*              dltb(                                                       */
/*                  df=...,         output file (required)                  */
/*                  fmt=...,        print format, def. 8.6                  */
/*                  grp=...,        list of group variables                 */
/*              ) = AGE,D;                                                  */
/*                                                                          */
/*                                                                          */
/*  If the grp parameter is used to specify groups, the calculations are    */  
/*  done separately for each group.                                         */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int dltb(void)
{
    register int i,j,k;
    int err,grp,idx,iia,iid,a,n,ne,nrec,amax,amax1;
    double r,g,m,w,tmp;

    err = -1;         
    if (check_cmd(1))
        return(-1);

    printf1("Life-table with discrete times. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,4,1))     /* get parameters */
        goto DLTBFin;
   
    if (PMF1Def == 0) {
        printf1("Error: need an output file.\n");
        goto DLTBFin;
    }
    if (PMFmtF == 0)
        pmfmt(10,6);

    if (PMNV != 2) {
        printf1("Error: need two variables on right-hand side.\n");
        goto DLTBFin;
    }
    newline();

    iia = PMVIdx[0];
    iid = PMVIdx[1];
    nrec = 0;

    amax = 0;
    for (i = 0; i < NOC; ++i) {
        a = (int)get_data(iia,i);
        amax = imax(a,amax);
    }
    if (amax < 1) {
        printf1("Error: maximal age value is %d\n",amax);
        goto DLTBFin;
    }
    if (alloc_acn(amax + 2))
        goto DLTBFin;
    if (alloc_acm(amax + 2))
        goto DLTBFin;
    if (alloc_acx(amax + 2))
        goto DLTBFin;
    if (alloc_acy(amax + 2))
        goto DLTBFin;
    if (alloc_acz(amax + 2))
        goto DLTBFin;

    if (PM1NV) {
        grp = 0;
        printf1("Group  defined by          Cases   Events\n");
    }
    else  
        grp = -1;
    idx = -1;

    for (; grp < PM1NV; ++grp) {               /* loop for all groups */
        idx++;

        for (i = 0; i <= amax; ++i)
            AcN[i] = AcM[i] = 0;

        amax1 = n = ne = 0;
        for (i = 0; i < NOC; ++i) {

            if (grp >= 0) {
                if (fabs(get_data(PM1VIdx[grp],i)) <= EPSI1)
                    continue;
            }
            n++;
            a = (int)get_data(iia,i);
            AcN[a] += 1;
            if ((int)get_data(iid,i)) {
                AcM[a] += 1;
                ne++;
            }
            amax1 = imax(amax1,a);
        }
        if (grp >= 0) {
            printf1("%5d  %s",grp,VName[PM1VIdx[grp]]);
            prnchar(' ',16 - strlen(VName[PM1VIdx[grp]]),0);
            printf1(" %8d  %8d\n",n,ne);
        }
        else {
            printf1("Number of cases: %d\n",n);
            printf1("Number of events: %d\n",ne);
        }
        if (n == 0)
            continue;


        g = 1.0;
        for (j = 0; j <= amax1; ++j) {
            if (AcN[j] > 0)
                r = (double)AcM[j] / (double)AcN[j];
            else
                r = 0.0;
            AcX[j] = r;
            AcY[j] = g;
            g *= (1.0 - r);
        }
        AcY[amax1 + 1] = 0.0;

        for (j = 0; j <= amax1; ++j) {
            w = m = 0.0;
            for (k = j; k <= amax1; ++k) {
                tmp = AcY[k] - AcY[k + 1];
                m += (double)k * tmp;                        
                w += tmp;                         
            }
            if (w > 0.0)
                AcZ[j] = m / w;
            else
                AcZ[j] = 0.0;
        }

        if (grp >= 0) {
            if (idx > 0) {
                fprintf(PMF1d,"\n");
                nrec++;
            }
            fprintf(PMF1d,"# Group %d defined by: %s\n",grp,VName[PM1VIdx[grp]]);
            nrec++;       
        }
        fprintf(PMF1d,"# Cases: %d  Events: %d\n",n,ne);
        fprintf(PMF1d,"# Idx    Age    Cases   Events ");
        fprnchar(PMF1d,' ',PMFmt1 - 4,0);
        fprintf(PMF1d,"Rate ");
        fprnchar(PMF1d,' ',PMFmt1 - 6,0);
        fprintf(PMF1d,"1-Rate ");
        fprnchar(PMF1d,' ',PMFmt1 - 2,0);
        fprintf(PMF1d,"SF\n");
        nrec += 2;

        for (j = 0; j <= amax1; ++j) {
            fprintf(PMF1d,"%5d ",idx);
            fprintf(PMF1d,"%6d ",j);
            fprintf(PMF1d,"%8d ",AcN[j]);
            fprintf(PMF1d,"%8d ",AcM[j]);
            fprintf(PMF1d,PMFmtS,AcX[j]);
            fprintf(PMF1d,PMFmtS,1.0 - AcX[j]);
            fprintf(PMF1d,PMFmtS,AcY[j]);
            fprintf(PMF1d,PMFmtS,AcZ[j]);
            tmp = AcZ[j] - (double)j;
            if (tmp < 0.0)
                tmp = 0.0;
            fprintf(PMF1d,PMFmtS,tmp);
            fprintf(PMF1d,"\n");
            nrec++;
        }
    }
    printf1("\n%d records written to: %s\n",nrec,PMF1dName);
    err = 0;

DLTBFin:
    p_clean();
    return(err);
}

/* -###-------------------------------------------------------------------- */
/*  diple()     Product-limit estimation with discrete data. Data may be    */
/*              interval-valued. Values are then equally distributed in     */
/*              the intervals.                                              */
/*                                                                          */  
/*              diple(                                                      */
/*                  df=...,         output file (required)                  */
/*                  fmt=...,        print format, def. 12.4                 */
/*              ) = S,TL,TU,D;                                              */
/*                                                                          */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int diple(void)
{
    register int i,t;
    int err,iis,iitl,iitu,iid,d,s,tl,tu,smin,tmax,tmaxe,nrec,ilen;
    double nlt,nrc,nex,nic,r,w,ws,wt,tmp;

    err = -1;         
    if (check_cmd(1))
        return(-1);

    printf1("PLE with discrete times. Current memory: %d bytes.\n",MemReq);

    TOLF = 0.001;

    if (parm(CmdBuf + 5,4,1))     /* get parameters */
        goto DIPLEFin;
   
    if (PMF1Def == 0) {
        printf1("Error: need an output file.\n");
        goto DIPLEFin;
    }
    if (PMFmtF == 0)
        pmfmt(12,4);

    if (PMNV == 4) {
        iis  = PMVIdx[0];
        iitl = PMVIdx[1];
        iitu = PMVIdx[2];
        iid  = PMVIdx[3];
    }
    else {
        printf1("Error: need four variables on right-hand side.\n");
        goto DIPLEFin;
    }
    newline();

    smin = INTMAX;
    tmaxe = tmax  = 0; 

    nlt = nrc = nic = nex = ws = 0.0;
    wt = 1.0;

    for (i = 0; i < NOC; ++i) {

        if (WIVar >= 0)                             /* get weights */
            wt = get_data(WIVar,i) * WNorm;
        if (wt <= 0.0)
            continue;

        ws += wt;

        s  = (int)get_data(iis,i);
        tl = (int)get_data(iitl,i);
        tu = (int)get_data(iitu,i);
        d  = (int)get_data(iid,i);

        if (s < 0 || tl < 0 || tu < tl || s > tl) {
            printf1("Error: inconsistent data in case %d\n",i + 1);
            goto DIPLEFin;
        }
        smin = imin(smin,s);
        if (d) {
            tmax = imax(tmax,tu);
            tmaxe = imax(tmaxe,tu);
            if (tu > tl)
                nic += wt;
            else
                nex += wt;
        }
        else {
            tmax = imax(tmax,tl);
            nrc += wt;
        }
        if (s > 0)
            nlt += wt;

    }
    printf1("Number of cases (weighted): %g\n",ws);
    printf1("Right censored: %g\n",nrc);
    printf1("Left truncated: %g\n",nlt);
    printf1("Exact duration: %g\n",nex);
    printf1("Interval-valued duration: %g\n\n",nic);
    printf1("Minimal time value of left truncation: %d\n",smin);
    printf1("Maximal time value of events: %d\n",tmaxe);
    printf1("Maximal time value for calculations: %d\n",tmax);

    if (ws <= EPSI1) {
        printf1("Number of cases is almost zero.\n");
        err = 0;
        goto DIPLEFin;
    }

    if (alloc_acx(tmax + 2))   
        goto DIPLEFin;
    if (alloc_acy(tmax + 2))   
        goto DIPLEFin;
    
    /* find mean distribution */

    wt = 1.0;
    for (i = 0; i < NOC; ++i) {

        if (WIVar >= 0)                     
            wt = get_data(WIVar,i) * WNorm;
        if (wt <= 0.0)
            continue;

        s  = (int)get_data(iis,i);
        tl = (int)get_data(iitl,i);
        tu = (int)get_data(iitu,i);
        d  = (int)get_data(iid,i);

        if (d == 0) {
            for (t = s; t <= tu; ++t)
                AcX[t] += wt;                
        }
        else if (tl == tu) {
            for (t = s; t <= tu; ++t)
                AcX[t] += wt;                
            AcY[tu] += wt;
        }
        else if (tl < tu) {

            w = wt;
            ilen = tu - tl + 1;
            tmp = wt / (double)ilen;

            for (t = s; t < tl; ++t)
                AcX[t] += wt;                

            for (t = tl; t <= tu; ++t) {
                AcY[t] += tmp;
                AcX[t] += w;
                w -= tmp;
            }
        }
    }
    nrec = 0;
    tmp = 1.0;
    for (t = smin; t <= tmax; ++t) {
        if (fabs(AcX[t]) <= EPSI1)   
            r = 0.0;
        else
            r = AcY[t] / AcX[t];

        fprintf(PMF1d,"%6d ",t);
        fprintf(PMF1d,PMFmtS,AcX[t]);
        fprintf(PMF1d,PMFmtS,AcY[t]);
        fprintf(PMF1d,PMFmtS,r);
        fprintf(PMF1d,PMFmtS,tmp);
        fprintf(PMF1d,"\n");
        nrec++;

        tmp *= (1.0 - r);
    }
    printf1("\n%d records written to: %s\n",nrec,PMF1dName);
    err = 0;

DIPLEFin:
    p_clean();
    return(err);
}


