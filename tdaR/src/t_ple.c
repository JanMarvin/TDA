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
#include "tda_context.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_ple.c                                                    */

int ple(TDAContext *ctx); 
int plest(TDAContext *ctx, int qtyp);
void prple_head(TDAContext *ctx, int sn, int org, int des, int grp);
void prple(TDAContext *ctx, int id, int idx, double time, double evnts, double cens, double rs, double surv, double sterr);
void prstar1(TDAContext *ctx);
void prquant_head(TDAContext *ctx);
void prquant(TDAContext *ctx, int qtyp,int sn, int org, int des, int grp);
int csf_stat(TDAContext *ctx);
int plestat(TDAContext *ctx);
int dple(TDAContext *ctx);
int dltb(TDAContext *ctx);
int diple(TDAContext *ctx);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */


/* ------------------------------------------------------------------------ */
/*  ple()           Product-limit estimation. ple command in CmdBuf.        */
/*                  Return 0 if OK, otherwise -1.                           */

int ple(TDAContext *ctx)
{
    register int i;
    int err,qtyp;

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Product-limit estimation. Current memory: %d bytes.\n",ctx->MemReq);

    if (ctx->EDAvail == 0) {
        p_err(ctx, -15,1);
        return(0);
    }
    if (parm(ctx, ctx->CmdBuf + 3,1,0))     /* get parameters */
        goto PLEFin;

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 7,5);

    prn_cwt(ctx);
    newline(ctx);

    if (ctx->ESortFlg == 0) {        /* sort data */
        if (sort_ed(ctx, 1))      
            goto PLEFin;
    }
    qtyp = 0;
    if (ctx->PMQOFlg || ctx->PMQTFlg) {           /* quantiles or orders */
        if (alloc_acx(ctx, ctx->PMNTP + 1))       /* allocate AcX */
            goto PLEFin;

        if (ctx->PMQOFlg) {
            qtyp = 1;
            for (i = 0; i < ctx->PMNTP; ++i) {
                if (ctx->PMTP[i] >= 1.0 || ctx->PMTP[i] <= 0.0 ||
                       (i > 0 && ctx->PMTP[i] >= ctx->PMTP[i - 1])) {
                    printf1(ctx, "Error in qo parameter.\n");
                    goto PLEFin;
                }
            }
        }
        else {
            qtyp = 2;
            for (i = 0; i < ctx->PMNTP; ++i) {
                if (ctx->PMTP[i] <= 0.0 || 
                       (i > 0 && ctx->PMTP[i] < ctx->PMTP[i - 1])) {
                    printf1(ctx, "Error in qt parameter.");
                    goto PLEFin;
                }
            }
        }
    }
    else
        free_tp(ctx);      /* set PMNTP = 0 */

    if (ctx->PMFDef || qtyp) {
        printf1(ctx, "\nProduct-limit estimation.\n");
        plest(ctx, qtyp);
    }
    if (ctx->PMCSF) {    
        printf1(ctx, "\nComparing survivor functions.\n");
        if (ctx->PM1NV < 2)  
            p_err(ctx, -31,1);
        else
            csf_stat(ctx);
    }
    err = 0;

PLEFin:
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "ple.table");
#endif
    if (ctx->ESortFlg)
        sort_ed(ctx, 0);     /* free memory */

    p_clean(ctx);
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

int plest(TDAContext *ctx, int qtyp)
{
    register int grp,icase,jcase,j;
    int isn,iorg,ides,sn;
    int org,des,des1,tran,id,idx,io,jo,nn,evflg,ntab;
    double q,tmp,surv,surv1,evnts,cens,cens1,cens2,rs,ws,wt,wt1,ts,tf,time,time1;
    double med,stsum,sterr,ord,ts1,tf1;

    if (ctx->PMFDef == 0 && qtyp == 0)
        return(0);

    if (qtyp)
        prquant_head(ctx);     /* print header for quantile table */

    id = -1;    /* id for tables */
    ntab = 0;   /* number of tables in output file */

    /* loop over all sn,org combinations */

    wt1 = wt = 1.0;

    for (tran = 0; tran < ctx->NTran1; ++tran) { /* loop over all transitions */

        isn  = ctx->SnTran1[tran];
        iorg = ctx->OrgTran1[tran];
        ides = ctx->DesTran1[tran];

        if (ctx->PM1NV) grp = 0; else grp = -1;

        for (; grp < ctx->PM1NV; ++grp) {               /* loop for all groups */

            /* print header for ple table */
                
            if (ctx->PMFDef) {
                prple_head(ctx, isn,iorg,ides,grp);
                ntab++;
            }
  
            /* Calculate the risk set for time 0 in rs. Take all        */
            /* episodes with origin state org (and contained in group   */
            /* grp) with starting times zero. The number of episodes is */
            /* calculated in nn, weighted in ws. It is assumed that     */
            /* episodes are sorted according to their starting times.   */

            rs = ws = 0.0;
            jcase = nn = 0;

            while (++jcase <= ctx->NOC) {

                jo = ctx->TSIdx[jcase - 1];
           
                get_edat(ctx, jo,&sn,&org,&des,&ts,&tf);

                if (sn == isn) {

                    if (ts > 0.0) {
                        jcase--;
                        break;
                    }
                    if (org != iorg)   
                        continue;

                    if (grp >= 0) {
                        if (fabs(get_data(ctx, ctx->PM1VIdx[grp],jo)) <= ctx->EPSI1)
                            continue;
                    }
                    if (ctx->WIVar >= 0)             /* get weights */
                        wt1 = get_data(ctx, ctx->WIVar,jo) * ctx->WNorm;
   
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
        
            for (j = 0; j < ctx->PMNTP; ++j)  
                ctx->AcX[j] = 0.0;
           
            evflg = j = 0;

            /* Now start the main loop for all episodes with origin     */
            /* state (and contained in grp), according to their ending  */
            /* times.                                                   */

            for (icase = 0; icase < ctx->NOC; ++icase) {

                io = ctx->TFIdx[icase];

                get_edat(ctx, io,&sn,&org,&des1,&ts,&tf);

                if (sn != isn || org != iorg)
                    continue;

                if (grp >= 0) {
                    if (fabs(get_data(ctx, ctx->PM1VIdx[grp],io)) <= ctx->EPSI1)
                        continue;
                }

                if (ctx->WIVar >= 0)             /* get weights */
                    wt = get_data(ctx, ctx->WIVar,io) * ctx->WNorm;

                if (tf > time) {    

                    if (evnts != 0 || time == 0.0) {

                        /* Update the risk set. Add all episodes    */
                        /* with starting time less than time.       */

                        while (jcase < ctx->NOC && ++jcase <= ctx->NOC) {

                            jo = ctx->TSIdx[jcase - 1];
                            get_edat(ctx, jo,&sn,&org,&des,&ts1,&tf1);
        
                            if (sn == isn) {

                                if (ts1 >= time) {
                                    jcase--;
                                    break;
                                }
                                if (org != iorg)   
                                    continue;

                                if (grp >= 0) {
                                    if (fabs(get_data(ctx, ctx->PM1VIdx[grp],jo)) <= ctx->EPSI1)
                                        continue;
                                }
                                if (ctx->WIVar >= 0)             /* get weights */
                                    wt1 = get_data(ctx, ctx->WIVar,jo) * ctx->WNorm;

                                rs += wt1;
                                ws += wt1;
                                nn++;
                            }
                        }             
                        if (rs <= 0.0 && jcase >= ctx->NOC)
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
                         
                        while (j < ctx->PMNTP) {
                            if (qtyp == 1) {
                                ord = ctx->PMTP[j];
                                if (surv1 <= ord && surv > ord) {
                                        tmp = time1 + (surv - ord) * 
                                            (time - time1) / (surv - surv1);
                                        ctx->AcX[j++] = tmp;
                                }
                                else
                                    break;
                            }
                            else if (qtyp == 2) {
                                q = ctx->PMTP[j];
                                if (surv1 < surv &&
                                              q > time1 && q <= time) {
                                    tmp = surv - (q - time1) * 
                                        (surv - surv1) / (time - time1);
                                    ctx->AcX[j++] = tmp;
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
   
                        if (ctx->PMFDef)
                            prple(ctx, id,idx,time,evnts,cens1,rs,surv,sterr);
           
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

            if (evnts != 0 || cens1 != 0.0) {
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
  
                if (ctx->PMFDef)
                    prple(ctx, id,idx,time,evnts,cens1,rs,surv,sterr);
            }
       
            /* print the median, and number of episodes */
       
            if (ctx->PMFDef) {
                if (med > 0.0) {
                    fprintf(ctx->PMFd,"\n# Median Duration: %4.2f",med);
                    if ((cens1) != 0.0)
                        fprintf(ctx->PMFd,"\n# Duration times limited to: %g",time1);
                }
                fprintf(ctx->PMFd,"\n# Cases: %d  weighted: %g\n\n",nn,ws);
#ifdef TDA_R_PACKAGE
                {
                    double srow[4];
                    srow[0] = med > 0.0 ? med : (double)(NAN);
                    srow[1] = (med > 0.0 && cens1 != 0.0) ? time1 : (double)NAN;
                    srow[2] = (double)nn;
                    srow[3] = ws;
                    tda_export_row(ctx, "ple.summary", srow, 4);
                }
#endif
            }

            /* print table of quantiles and count parameters */
  
            prquant(ctx, qtyp,isn,iorg,ides,grp);
  
        }   /* end of group selection */
    }       /* end of loop for transitions */

    if (ctx->PMFDef) {
        if (qtyp)  
            printf1(ctx, "\n");           
        printf1(ctx, "%d table(s) written to: %s\n",ntab,ctx->PMFdName);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prple_head.  Print header of product-limit estimation table             */

void prple_head(TDAContext *ctx, int sn, int org, int des, int grp)
{
#ifdef TDA_R_PACKAGE
    /* One row per printed block, in block order: the transition it is
       for, and how many rows of ple.table precede it.  ple.table is one
       concatenated export, so without this the reader has to take the
       block boundaries from the parsed text to split it. */
    {
        double erow[5];
        erow[0] = (double)sn;
        erow[1] = (double)org;
        erow[2] = (double)des;
        erow[3] = (double)grp;
        erow[4] = (double)ctx->RExportPleRows;
        tda_export_row(ctx, "ple.blocks", erow, 5);
    }
#endif
    fprintf(ctx->PMFd,"# SN %d. Transition: %d,%d - Product-Limit Estimation\n",
                                                              sn,org,des);
    if (grp >= 0)
        fprintf(ctx->PMFd,"\n# Group: %s",ctx->VName[ctx->PM1VIdx[grp]]);
    fprintf(ctx->PMFd,"\n#                      Number  Number    Exposed ");
    fprnchar(ctx, ctx->PMFd,' ',ctx->PMFmt1 - 7,0); fprintf(ctx->PMFd," Survivor");
    fprnchar(ctx, ctx->PMFd,' ',ctx->PMFmt1 - 5,0); fprintf(ctx->PMFd," Std. ");
    fprnchar(ctx, ctx->PMFd,' ',ctx->PMFmt1 - 5,0); fprintf(ctx->PMFd," Cum. ");
    fprintf(ctx->PMFd,"\n# ID Index      Time   Events  Censored  to Risk ");
    fprnchar(ctx, ctx->PMFd,' ',ctx->PMFmt1 - 7,0); fprintf(ctx->PMFd," Function");
    fprnchar(ctx, ctx->PMFd,' ',ctx->PMFmt1 - 5,0); fprintf(ctx->PMFd,"Error ");
    fprnchar(ctx, ctx->PMFd,' ',ctx->PMFmt1 - 5,0); fprintf(ctx->PMFd," Rate ");
}

/* ------------------------------------------------------------------------ */
/*  prple.  Print one line of the product-limit estimation table.           */

void prple(TDAContext *ctx, int id, int idx, double time, double evnts, double cens, double rs, double surv, double sterr)
{
    double tmp;

    if (evnts < 0.0 || rs < 0.0 || surv < 0.0 || cens < 0.0)
        return;

    if (evnts > 0 || surv == 1.0)
        fprintf(ctx->PMFd,"\n%4d %5d",id,idx);
    else
        fprintf(ctx->PMFd,"\n# %2d %5d",id,idx);
    fprintf(ctx->PMFd," %9.2f",time);
    fprintf(ctx->PMFd," %8d %8d",(int)evnts,(int)cens);
    if (evnts > 0 || surv == 1.0) {
#ifdef TDA_R_PACKAGE
        /* Exactly the rows the file's own reader keeps.  The
           censoring-only rows below are printed with a leading "#" and
           skipped by that reader, so exporting them here would make the
           two paths disagree -- tried, and it put NA rows into $blocks,
           $table and tda_survivor(), which broke plotting outright.
           They remain reachable through the generic file tap
           (file.PMFd.values); nothing is lost, it just does not belong
           in this table. */
        {
            double erow[9];
            erow[0] = (double)id;
            erow[1] = (double)idx;
            erow[2] = time;
            erow[3] = evnts;
            erow[4] = cens;
            erow[5] = rs;
            erow[6] = surv;
            erow[7] = sterr >= 0.0 ? sterr : (double)NAN;
            erow[8] = surv > 0.0 ?
                (surv == 1.0 ? 0.0 : -rlog(ctx, surv)) : (double)NAN;
            tda_export_row(ctx, "ple.table", erow, 9);
            /* counted HERE, not at the top of prple(): counting every
               PRINTED row put ple.blocks' offsets past the end of each
               block, since the censoring-only rows print but are not
               exported. */
            ctx->RExportPleRows++;
        }
#endif
        fprintf(ctx->PMFd," %9d  ",(int)rs);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,surv);
        if (sterr >= 0.0)
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,sterr);
        else
            prstar1(ctx);
        if (surv > 0.0) {
            if (surv == 1.0)
                tmp = 0.0;
            else
                tmp = -rlog(ctx, surv);
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
        }
        else  
            prstar1(ctx);
    }
}

void prstar1(TDAContext *ctx)
{
    fprnchar(ctx, ctx->PMFd,' ',ctx->PMFmt1 - 2,0);
    fprintf(ctx->PMFd," * ");
}

/* ------------------------------------------------------------------------ */
/*  prquant_head.  Print header of table for quantiles.                     */

void prquant_head(TDAContext *ctx)
{
    if (ctx->PMNTP) {
        printf1(ctx, "\n");
        prnchar(ctx, ' ',14,0);
        if (ctx->PM1NV > 0)   
            prnchar(ctx, ' ',ctx->VNameLen + 7,0);
  
        printf1(ctx, "  Survivor      Time\n");
        printf1(ctx, "SN  Org  Des  ");
  
        if (ctx->PM1NV > 0) {            
            printf1(ctx, "Group  ");
            prnchar(ctx, ' ',ctx->VNameLen,0);
        }
        printf1(ctx, "  Function    Quantile\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  prquant.  Print table of quantiles.                                     */

void prquant(TDAContext *ctx, int qtyp,int sn, int org, int des, int grp)
{
    register int j;

    if (ctx->PMNTP <= 0)
        return;

    prnchar(ctx, '-',36,0);
  
    if (ctx->PM1NV > 0)   
        prnchar(ctx, '-',ctx->VNameLen + 7,0);

    for (j = 0; j < ctx->PMNTP; ++j) {

        printf1(ctx, "\n");           

        if (!j) {
            printf1(ctx, "%2d %4d %4d",sn,org,des);
            if (ctx->PM1NV > 0) {
                printf1(ctx, "%4d  ",grp + 1);
                prn_vname(ctx, ctx->PM1VIdx[grp]);
            }
        }
        else {
            prnchar(ctx, ' ',12,0);
            if (ctx->PM1NV > 0)   
                prnchar(ctx, ' ',ctx->VNameLen + 7,0);
        }
        if (ctx->AcX[j] <= 0.0 || ctx->PMTP[j] <= 0.0)  
            break;
            /*  printf1("        ***          ***");       */
        else if (qtyp == 1)
            printf1(ctx, " %10.4f %12.4f",ctx->PMTP[j],ctx->AcX[j]);
        else
            printf1(ctx, " %10.4f %12.4f",ctx->AcX[j],ctx->PMTP[j]);
#ifdef TDA_R_PACKAGE
        /* the printed pair in printed order (qtyp swaps the columns),
           with the row's sn/org/des/group lead; full precision where
           the text has four decimals */
        {
            double erow[6];
            erow[0] = (double)sn;
            erow[1] = (double)org;
            erow[2] = (double)des;
            erow[3] = (double)(ctx->PM1NV > 0 ? grp + 1 : 0);
            if (qtyp == 1) {
                erow[4] = ctx->PMTP[j];
                erow[5] = ctx->AcX[j];
            }
            else {
                erow[4] = ctx->AcX[j];
                erow[5] = ctx->PMTP[j];
            }
            tda_export_row(ctx, "ple.quantiles", erow, 6);
        }
#endif
    }
    printf1(ctx, "\n");
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "ple.quantiles");
#endif
}

/* ------------------------------------------------------------------------ */
/*  csf_stat   Test statistics for comparing survivor functions.            */
/*                                                                          */

#define PLENTest 4      /* number of test statistics */

int csf_stat(TDAContext *ctx)
{
    register int j,tran;
    int err,n,nn,sn,org,des,tsa,rsa;
    double tmp;

    nn = ctx->NTran1 * PLENTest;
    tsa = rsa = err = 0;

    if (!(ctx->TStat = (double *)calloc((size_t)(nn),sizeof(double)))) { 
        err = -2;
        goto CSFFin;
    }
    tsa = 1;
    memrq(ctx, nn,sizeof(double));

    if (!(ctx->RStat = (int *)calloc((size_t)(nn),sizeof(int)))) { 
        err = -2;
        goto CSFFin;
    }
    rsa = 1;
    memrq(ctx, nn,sizeof(int));

    if (plestat(ctx))      /* perform calculations */
        goto CSFFin;

    printf1(ctx, "\nSN  Org Des   Test Statistic                 T-Stat   DF   Signif");

    for (tran = 0; tran < ctx->NTran1; ++tran) { /* loop over all transitions */

        sn  = ctx->SnTran1[tran];
        org = ctx->OrgTran1[tran];
        des = ctx->DesTran1[tran];

        printf1(ctx, "\n");
        prnchar(ctx, '-',65,0);

        for (j = 0; j < PLENTest; ++j) {
            printf1(ctx, "\n%2d  %3d %3d    ",sn,org,des);
            if (!j)  
                printf1(ctx, "Log-Rank (Savage)     ");
            else {
                printf1(ctx, "Wilcoxon (");
                if (j == 1)  
                    printf1(ctx, "Breslow)    ");          
                else if (j == 2)  
                    printf1(ctx, "Tarone-Ware)");          
                else if (j == 3)  
                    printf1(ctx, "Prentice)   ");          
            }
#ifdef TDA_R_PACKAGE
            /* the same four names the block prints, so the R side can
               rebuild the table from the export alone: a row whose
               statistic is nan prints "-nan" and "*", which no numeric
               parser can read back */
            tda_export_str_row(ctx, "ple.comparison.names",
                               j == 0 ? "Log-Rank (Savage)" :
                               j == 1 ? "Wilcoxon (Breslow)" :
                               j == 2 ? "Wilcoxon (Tarone-Ware)" :
                                        "Wilcoxon (Prentice)");
#endif
            tmp = ctx->TStat[tran * PLENTest + j];
            printf1(ctx, "%14.4f ",tmp);
            n = ctx->RStat[tran * PLENTest + j];
            printf1(ctx, "%4d ",n);
            if (n > 0 && tmp > 0.0)  
                printf1(ctx, "%8.4f",cdchif(ctx, tmp,n));
            else
                printf1(ctx, "       *");
#ifdef TDA_R_PACKAGE
            {
                double erow[6];
                erow[0] = (double)sn;
                erow[1] = (double)org;
                erow[2] = (double)des;
                erow[3] = tmp;
                erow[4] = (double)n;
                erow[5] = (n > 0 && tmp > 0.0) ? cdchif(ctx, tmp,n) : (double)(NAN);
                tda_export_row(ctx, "ple.comparison", erow, 6);
            }
#endif
        }
    }
    printf1(ctx, "\n");         
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "ple.comparison");
    tda_export_str_flush(ctx, "ple.comparison.names");
#endif

CSFFin:
    if (err == -2)  
        p_err(ctx, -2,1);
       
    if (tsa) {
        free((char *)ctx->TStat);
        memrq(ctx, -nn,sizeof(double));
    }
    if (rsa) {
        free((char *)ctx->RStat);
        memrq(ctx, -nn,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  plestat     Calculation of test statistics for comparison of            */
/*              survivor functions.                                         */
/*              Return 0 if OK, -1 if error (insuff. memory)                */

int plestat(TDAContext *ctx)
{
    register int l,i,j,icase,jcase;
    int err,rsa,eva,cena,va[PLENTest],ua[PLENTest];
    int sn,isn,iorg,ides,ii,nn,g,org,des,des1,tran,io,jo;
    double tmp,tmp1,tmp2,wt,wt1,ts,tf,ts1,tf1,time,r0,e0,w[PLENTest];
    double *rs,*ev,*cen,*v[PLENTest],*u[PLENTest];

    rsa = eva = cena = err = 0;
    nn = ctx->PM1NV - 1;

    if (!(rs = (double *)calloc((size_t)(ctx->PM1NV + 1),sizeof(double)))) {
        err = -1;
        goto PLESTATFin;
    }
    rsa = 1;
    memrq(ctx, ctx->PM1NV + 1,sizeof(double));

    if (!(ev = (double *)calloc((size_t)(ctx->PM1NV + 1),sizeof(double)))) {
        err = -1;
        goto PLESTATFin;
    }
    eva = 1;
    memrq(ctx, ctx->PM1NV + 1,sizeof(double));

    if (!(cen = (double *)calloc((size_t)(ctx->PM1NV + 1),sizeof(double)))) {
        err = -1;
        goto PLESTATFin;
    }
    cena = 1;
    memrq(ctx, ctx->PM1NV + 1,sizeof(double));

    for (i = 0; i < PLENTest; ++i) {
        if (!(v[i] = (double *)calloc((size_t)(nn) * (size_t)(nn) + 1,sizeof(double)))) {
            err = -1;
            goto PLESTATFin;
        }
        va[i] = 1;
        memrq(ctx, nn * nn + 1,sizeof(double));

        if (!(u[i] = (double *)calloc((size_t)(nn + 1),sizeof(double)))) {
            err = -1;
            goto PLESTATFin;
        }
        ua[i] = 1;
        memrq(ctx, nn + 1,sizeof(double));
    }
     
    wt = wt1 = 1.0;             

    for (tran = 0; tran < ctx->NTran1; ++tran) { /* loop over all transitions */

        isn  = ctx->SnTran1[tran];
        iorg = ctx->OrgTran1[tran];
        ides = ctx->DesTran1[tran];

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

        for (i = 0; i <= ctx->PM1NV; ++i)
            rs[i] = 0.0;

        jcase = 0;
        while (++jcase <= ctx->NOC) {

            jo = ctx->TSIdx[jcase - 1];

            get_edat(ctx, jo,&sn,&org,&des,&ts,&tf);

            if (sn == isn) {

                if (ts > 0.0) {
                    jcase--;
                    break;
                }
                if (org != iorg)   
                    continue;
  
                for (l = 0; l < ctx->PM1NV; ++l) {
                    if (fabs(get_data(ctx, ctx->PM1VIdx[l],jo)) > ctx->EPSI1)
                        break;
                }
                if (l < ctx->PM1NV) {

                    if (ctx->WIVar >= 0)             /* get weights */
                        wt1 = get_data(ctx, ctx->WIVar,jo) * ctx->WNorm;

                    rs[0] += wt1;
                    rs[l + 1] += wt1;
                }
            }   
        }
        time = 0.0; 
        for (i = 0; i <= ctx->PM1NV; ++i)
            ev[i] = cen[i] = 0.0;

        for (icase = 0; icase < ctx->NOC; ++icase) {

            io = ctx->TFIdx[icase];

            get_edat(ctx, io,&sn,&org,&des1,&ts,&tf);

            if (sn != isn || org != iorg)
                continue;

            for (l = 0; l < ctx->PM1NV; ++l) {
                if (fabs(get_data(ctx, ctx->PM1VIdx[l],io)) > ctx->EPSI1)
                    break;
            }
            if (l >= ctx->PM1NV)  
                continue;

            g = l + 1;              /* group of episode */

            if (ctx->WIVar >= 0)             /* get weights */
                wt = get_data(ctx, ctx->WIVar,io) * ctx->WNorm;

            if (tf > time) {    

                if ((ev[0]) != 0.0) {        /* update the risk set */

                    while (jcase < ctx->NOC && ++jcase <= ctx->NOC) {

                        jo = ctx->TSIdx[jcase - 1];
                        get_edat(ctx, jo,&sn,&org,&des,&ts1,&tf1);
        
                        if (sn == isn) {

                            if (ts1 >= time) {
                                jcase--;
                                break;
                            }
                            if (org != iorg)   
                                continue;

                            for (l = 0; l < ctx->PM1NV; ++l) {
                                if (fabs(get_data(ctx, ctx->PM1VIdx[l],jo)) > ctx->EPSI1)
                                    break;
                            }
                            if (l < ctx->PM1NV) {

                                if (ctx->WIVar >= 0)             /* get weights */
                                    wt1 = get_data(ctx, ctx->WIVar,jo) * ctx->WNorm;

                                rs[0] += wt1;
                                rs[l + 1] += wt1;
                            }
                        }   
                    }
                                /* update U and V */

                    r0 = rs[0];
                    e0 = ev[0];
                    w[1] = r0 / ctx->WSum;                 
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
                    if (rs[0] <= 0.0 && jcase >= ctx->NOC)
                        break;
                }
                for (i = 0; i <= ctx->PM1NV; ++i) {
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
            w[1] = r0 / ctx->WSum;   
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
        if (ctx->PMProtFDef) {
            fprintf(ctx->PMProtFd,"SN %d Transition: %d,%d\n",isn,iorg,ides);
            for (l = 0; l < PLENTest; ++l) {
                fprintf(ctx->PMProtFd,"\nTest statistic %d",l + 1);
                fprintf(ctx->PMProtFd,"\n U: ");
                for (i = 1; i <= nn; ++i)
                    rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,u[l][i]);
    
                for (i = 1; i <= nn; ++i) {
                    fprintf(ctx->PMProtFd,"\n V: ");
                    for (j = 1; j <= nn; ++j)
                        rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,v[l][(i - 1) * nn + j]);
                }
                fprintf(ctx->PMProtFd,"\n");
            }
            fprintf(ctx->PMProtFd,"\n");
        }

        ii = tran * PLENTest;
        for (l = 0; l < PLENTest; ++l) {
            ctx->RStat[ii + l] = syminv(ctx, nn,v[l]);
            ctx->TStat[ii + l] = 0.0;
            for (i = 1; i <= nn; ++i) {
                tmp = 0.0;
                for (j = 1; j <= nn; ++j)  
                    tmp += v[l][(i - 1) * nn + j] * u[l][j];
                ctx->TStat[ii + l] += tmp * u[l][i];
            }
        }
    }       /* end of tran loop */

PLESTATFin:
    if (err == -1)  
        p_err(ctx, -2,1);
    if (rsa) {
        free((char *)rs);
        memrq(ctx, -ctx->PM1NV - 1,sizeof(double));
    }
    if (eva) {
        free((char *)ev);
        memrq(ctx, -ctx->PM1NV - 1,sizeof(double));
    }
    if (cena) {
        free((char *)cen);
        memrq(ctx, -ctx->PM1NV - 1,sizeof(double));
    }
    for (i = 0; i < PLENTest; ++i) {
        if (va[i]) {
            free((char *)v[i]);
            memrq(ctx, -nn * nn - 1,sizeof(double));
        }
        if (ua[i]) {
            free((char *)u[i]);
            memrq(ctx, -nn - 1,sizeof(double));
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

int dple(TDAContext *ctx)
{
    register int i,j;
    int err,grp,idx,iis,iif,iid,t,ta,tmin,tmax,nt,nrec;
    double r,g,m,rr,gg,n,wt,ncen,nev;

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "PLE with discrete times. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,4,1))     /* get parameters */
        goto DPLEFin;
   
    if (ctx->PMF1Def == 0) {
        printf1(ctx, "Error: need an output file.\n");
        goto DPLEFin;
    }
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 12,4);

    if (ctx->PMNV == 2) {
        iis = -1;         
        iif = ctx->PMVIdx[0];
        iid = ctx->PMVIdx[1];
    }
    else if (ctx->PMNV == 3) {
        iis = ctx->PMVIdx[0];
        iif = ctx->PMVIdx[1];
        iid = ctx->PMVIdx[2];
    }
    else {
        printf1(ctx, "Error: need two or three variables on right-hand side.\n");
        goto DPLEFin;
    }
    newline(ctx);
    nrec = 0;

    if (ctx->PM1NV) {
        grp = 0;
        printf1(ctx, "Group  defined by      ");
        prnchar(ctx, ' ',ctx->PMFmt1 - 5,0);
        printf1(ctx, "Cases ");
        printf1(ctx, " Minimum  Maximum\n");
    }
    else  
        grp = -1;
    idx = -1;

    for (; grp < ctx->PM1NV; ++grp) {               /* loop for all groups */
        idx++;
        tmax = 0;                        
        if (iis == -1)
            tmin = 0;
        else
            tmin = ctx->INTMAX;

        wt = 1.0;
        n = 0.0;
        for (i = 0; i < ctx->NOC; ++i) {

            if (grp >= 0) {
                if (fabs(get_data(ctx, ctx->PM1VIdx[grp],i)) <= ctx->EPSI1)
                    continue;
            }
            if (ctx->WIVar >= 0) {                           /* get weights */
                wt = get_data(ctx, ctx->WIVar,i) * ctx->WNorm;
                n += wt;
            }
            else
                n += 1.0;

            ta = 0;
            if (iis >= 0) {
                ta = (int)get_data(ctx, iis,i);
                tmin = imin(ctx, ta,tmin);
                tmax = imax(ctx, ta,tmax);
            }
            t = (int)get_data(ctx, iif,i);
            tmin = imin(ctx, t,tmin);
            tmax = imax(ctx, t,tmax);
           
            if (ta < 0 || t < 0 || t < ta) {
                printf1(ctx, "Error: found inconsistent data in case %d.\n",i + 1);
                goto DPLEFin;
            }
        }
        if (grp >= 0) {
            printf1(ctx, "%5d  %s",grp,ctx->VName[ctx->PM1VIdx[grp]]);
            prnchar(ctx, ' ',(int)(16 - strlen(ctx->VName[ctx->PM1VIdx[grp]])),0);
            rt_printf1_d(ctx, ctx->PMFmtS,n);
            if (n > 0.0)
                printf1(ctx, "%8d %8d",tmin,tmax);
            printf1(ctx, "\n");
        }
        else {
            printf1(ctx, "Number of cases (weighted): %g\n",n);
            printf1(ctx, "Minimal time value: %d\n",tmin);
            printf1(ctx, "Maximal time value: %d\n",tmax);
        }
        if (n <= ctx->EPSI1)
            continue;

        nt = tmax - tmin + 1;

        if (alloc_acx(ctx, nt + 1))
            goto DPLEFin;
        if (alloc_acy(ctx, nt + 1))
            goto DPLEFin;
        if (alloc_acz(ctx, nt + 1))
            goto DPLEFin;

        wt = 1.0;
        for (i = 0; i < ctx->NOC; ++i) {

            if (grp >= 0) {
                if (fabs(get_data(ctx, ctx->PM1VIdx[grp],i)) <= ctx->EPSI1)
                    continue;
            }
            if (ctx->WIVar >= 0)                             /* get weights */
                wt = get_data(ctx, ctx->WIVar,i) * ctx->WNorm;

            t = (int)get_data(ctx, iif,i);
            if (get_data(ctx, iid,i) != 0.0)
                ctx->AcX[t - tmin] += wt;
            else
                ctx->AcY[t - tmin] += wt;

            if (iis >= 0)
                ta = (int)get_data(ctx, iis,i);
            else
                ta = 0;

            for (j = ta; j <= t; ++j) 
                ctx->AcZ[j - tmin] += wt;

        }
        if (grp >= 0) {
            if (idx > 0) {
                fprintf(ctx->PMF1d,"\n");
                nrec++;
            }
            fprintf(ctx->PMF1d,"# Group %d defined by: %s\n",grp,ctx->VName[ctx->PM1VIdx[grp]]);
            nrec++;       
        }
        fprintf(ctx->PMF1d,"#");
        fprnchar(ctx, ctx->PMF1d,' ',11 + 5 * (ctx->PMFmt1 + 1) - 8,0);
        fprintf(ctx->PMF1d,"survivor\n");
        fprintf(ctx->PMF1d,"# Idx   Time ");
        fprnchar(ctx, ctx->PMF1d,' ',ctx->PMFmt1 - 7,0);
        fprintf(ctx->PMF1d,"at risk ");
        fprnchar(ctx, ctx->PMF1d,' ',ctx->PMFmt1 - 6,0);
        fprintf(ctx->PMF1d,"events ");
        fprnchar(ctx, ctx->PMF1d,' ',ctx->PMFmt1 - 8,0);
        fprintf(ctx->PMF1d,"censored ");
        fprnchar(ctx, ctx->PMF1d,' ',ctx->PMFmt1 - 4,0);
        fprintf(ctx->PMF1d,"rate ");
        fprnchar(ctx, ctx->PMF1d,' ',ctx->PMFmt1 - 8,0);
        fprintf(ctx->PMF1d,"function\n");
        nrec += 2;

        g = 1.0;
        ncen = nev = 0.0;
        for (i = 0; i < nt; ++i) {
            if (ctx->PMOPT == 2) {
                if (ctx->AcX[i] == 0.0 && ctx->AcY[i] == 0.0)
                    continue;
            }
            nev += ctx->AcX[i];
            ncen += ctx->AcY[i];

            fprintf(ctx->PMF1d,"%5d ",idx);
            fprintf(ctx->PMF1d,"%6d ",tmin + i);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcZ[i]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[i]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[i]);
#ifdef TDA_R_PACKAGE
            /* opt=3 appends two mean columns, so the row's width is
               only known as it is written */
            tda_export_cell(ctx, "dple.table", (double)idx);
            tda_export_cell(ctx, "dple.table", (double)(tmin + i));
            tda_export_cell(ctx, "dple.table", ctx->AcZ[i]);
            tda_export_cell(ctx, "dple.table", ctx->AcX[i]);
            tda_export_cell(ctx, "dple.table", ctx->AcY[i]);
#endif

            r = 0.0;
            if (ctx->AcZ[i] > 0.0)  
                r = ctx->AcX[i] / ctx->AcZ[i];

            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,r);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,g);
#ifdef TDA_R_PACKAGE
            tda_export_cell(ctx, "dple.table", r);
            tda_export_cell(ctx, "dple.table", g);
#endif

            if (ctx->PMOPT == 3) {           /* calculate mean values */

                gg = 1.0;
                m = 0.0;
                for (j = i; j < nt; ++j) {

                    rr = 0.0;
                    if (ctx->AcZ[j] > 0.0)  
                        rr = ctx->AcX[j] / ctx->AcZ[j];
                    if (j == nt - 1)
                        rr = 1.0;
                    m += (double)(tmin + j) * gg * rr;

                    gg *= (1.0 - rr);
                }
                fprintf(ctx->PMF1d,"%12.4f ",m);
                fprintf(ctx->PMF1d,"%12.4f ",m - (double)(tmin + i));
#ifdef TDA_R_PACKAGE
                tda_export_cell(ctx, "dple.table", m);
                tda_export_cell(ctx, "dple.table", m - (double)(tmin + i));
#endif
            }
            fprintf(ctx->PMF1d,"\n");
#ifdef TDA_R_PACKAGE
            tda_export_endrow(ctx, "dple.table");
#endif
            nrec++;
            g *= (1.0 - r);
        }
        fprintf(ctx->PMF1d,"# Number of events: %g\n",nev);
        fprintf(ctx->PMF1d,"# Number of censored cases: %g\n",ncen);
    }
    printf1(ctx, "\n%d records written to: %s\n",nrec,ctx->PMF1dName);
    err = 0;

DPLEFin:
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "dple.table");
#endif
    p_clean(ctx);
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

int dltb(TDAContext *ctx)
{
    register int i,j,k;
    int err,grp,idx,iia,iid,a,n,ne,nrec,amax,amax1;
    double r,g,m,w,tmp;

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Life-table with discrete times. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,4,1))     /* get parameters */
        goto DLTBFin;
   
    if (ctx->PMF1Def == 0) {
        printf1(ctx, "Error: need an output file.\n");
        goto DLTBFin;
    }
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,6);

    if (ctx->PMNV != 2) {
        printf1(ctx, "Error: need two variables on right-hand side.\n");
        goto DLTBFin;
    }
    newline(ctx);

    iia = ctx->PMVIdx[0];
    iid = ctx->PMVIdx[1];
    nrec = 0;

    amax = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        a = (int)get_data(ctx, iia,i);
        amax = imax(ctx, a,amax);
    }
    if (amax < 1) {
        printf1(ctx, "Error: maximal age value is %d\n",amax);
        goto DLTBFin;
    }
    if (alloc_acn(ctx, amax + 2))
        goto DLTBFin;
    if (alloc_acm(ctx, amax + 2))
        goto DLTBFin;
    if (alloc_acx(ctx, amax + 2))
        goto DLTBFin;
    if (alloc_acy(ctx, amax + 2))
        goto DLTBFin;
    if (alloc_acz(ctx, amax + 2))
        goto DLTBFin;

    if (ctx->PM1NV) {
        grp = 0;
        printf1(ctx, "Group  defined by          Cases   Events\n");
    }
    else  
        grp = -1;
    idx = -1;

    for (; grp < ctx->PM1NV; ++grp) {               /* loop for all groups */
        idx++;

        for (i = 0; i <= amax; ++i)
            ctx->AcN[i] = ctx->AcM[i] = 0;

        amax1 = n = ne = 0;
        for (i = 0; i < ctx->NOC; ++i) {

            if (grp >= 0) {
                if (fabs(get_data(ctx, ctx->PM1VIdx[grp],i)) <= ctx->EPSI1)
                    continue;
            }
            n++;
            a = (int)get_data(ctx, iia,i);
            ctx->AcN[a] += 1;
            if ((int)get_data(ctx, iid,i)) {
                ctx->AcM[a] += 1;
                ne++;
            }
            amax1 = imax(ctx, amax1,a);
        }
        if (grp >= 0) {
            printf1(ctx, "%5d  %s",grp,ctx->VName[ctx->PM1VIdx[grp]]);
            prnchar(ctx, ' ',(int)(16 - strlen(ctx->VName[ctx->PM1VIdx[grp]])),0);
            printf1(ctx, " %8d  %8d\n",n,ne);
        }
        else {
            printf1(ctx, "Number of cases: %d\n",n);
            printf1(ctx, "Number of events: %d\n",ne);
        }
        if (n == 0)
            continue;


        g = 1.0;
        for (j = 0; j <= amax1; ++j) {
            if (ctx->AcN[j] > 0)
                r = (double)ctx->AcM[j] / (double)ctx->AcN[j];
            else
                r = 0.0;
            ctx->AcX[j] = r;
            ctx->AcY[j] = g;
            g *= (1.0 - r);
        }
        ctx->AcY[amax1 + 1] = 0.0;

        for (j = 0; j <= amax1; ++j) {
            w = m = 0.0;
            for (k = j; k <= amax1; ++k) {
                tmp = ctx->AcY[k] - ctx->AcY[k + 1];
                m += (double)k * tmp;                        
                w += tmp;                         
            }
            if (w > 0.0)
                ctx->AcZ[j] = m / w;
            else
                ctx->AcZ[j] = 0.0;
        }

        if (grp >= 0) {
            if (idx > 0) {
                fprintf(ctx->PMF1d,"\n");
                nrec++;
            }
            fprintf(ctx->PMF1d,"# Group %d defined by: %s\n",grp,ctx->VName[ctx->PM1VIdx[grp]]);
            nrec++;       
        }
        fprintf(ctx->PMF1d,"# Cases: %d  Events: %d\n",n,ne);
        fprintf(ctx->PMF1d,"# Idx    Age    Cases   Events ");
        fprnchar(ctx, ctx->PMF1d,' ',ctx->PMFmt1 - 4,0);
        fprintf(ctx->PMF1d,"Rate ");
        fprnchar(ctx, ctx->PMF1d,' ',ctx->PMFmt1 - 6,0);
        fprintf(ctx->PMF1d,"1-Rate ");
        fprnchar(ctx, ctx->PMF1d,' ',ctx->PMFmt1 - 2,0);
        fprintf(ctx->PMF1d,"SF\n");
        nrec += 2;

        for (j = 0; j <= amax1; ++j) {
            fprintf(ctx->PMF1d,"%5d ",idx);
            fprintf(ctx->PMF1d,"%6d ",j);
            fprintf(ctx->PMF1d,"%8d ",ctx->AcN[j]);
            fprintf(ctx->PMF1d,"%8d ",ctx->AcM[j]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[j]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,1.0 - ctx->AcX[j]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[j]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcZ[j]);
            tmp = ctx->AcZ[j] - (double)j;
            if (tmp < 0.0)
                tmp = 0.0;
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,tmp);
            fprintf(ctx->PMF1d,"\n");
#ifdef TDA_R_PACKAGE
            {
                double erow[9];
                erow[0] = (double)idx;
                erow[1] = (double)j;
                erow[2] = (double)ctx->AcN[j];
                erow[3] = (double)ctx->AcM[j];
                erow[4] = ctx->AcX[j];
                erow[5] = 1.0 - ctx->AcX[j];
                erow[6] = ctx->AcY[j];
                erow[7] = ctx->AcZ[j];
                erow[8] = tmp;
                tda_export_row(ctx, "dltb.table", erow, 9);
            }
#endif
            nrec++;
        }
    }
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "dltb.table");
#endif
    printf1(ctx, "\n%d records written to: %s\n",nrec,ctx->PMF1dName);
    err = 0;

DLTBFin:
    p_clean(ctx);
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

int diple(TDAContext *ctx)
{
    register int i,t;
    int err,iis,iitl,iitu,iid,d,s,tl,tu,smin,tmax,tmaxe,nrec,ilen;
    double nlt,nrc,nex,nic,r,w,ws,wt,tmp;

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "PLE with discrete times. Current memory: %d bytes.\n",ctx->MemReq);

    ctx->TOLF = 0.001;

    if (parm(ctx, ctx->CmdBuf + 5,4,1))     /* get parameters */
        goto DIPLEFin;
   
    if (ctx->PMF1Def == 0) {
        printf1(ctx, "Error: need an output file.\n");
        goto DIPLEFin;
    }
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 12,4);

    if (ctx->PMNV == 4) {
        iis  = ctx->PMVIdx[0];
        iitl = ctx->PMVIdx[1];
        iitu = ctx->PMVIdx[2];
        iid  = ctx->PMVIdx[3];
    }
    else {
        printf1(ctx, "Error: need four variables on right-hand side.\n");
        goto DIPLEFin;
    }
    newline(ctx);

    smin = ctx->INTMAX;
    tmaxe = tmax  = 0; 

    nlt = nrc = nic = nex = ws = 0.0;
    wt = 1.0;

    for (i = 0; i < ctx->NOC; ++i) {

        if (ctx->WIVar >= 0)                             /* get weights */
            wt = get_data(ctx, ctx->WIVar,i) * ctx->WNorm;
        if (wt <= 0.0)
            continue;

        ws += wt;

        s  = (int)get_data(ctx, iis,i);
        tl = (int)get_data(ctx, iitl,i);
        tu = (int)get_data(ctx, iitu,i);
        d  = (int)get_data(ctx, iid,i);

        if (s < 0 || tl < 0 || tu < tl || s > tl) {
            printf1(ctx, "Error: inconsistent data in case %d\n",i + 1);
            goto DIPLEFin;
        }
        smin = imin(ctx, smin,s);
        if (d) {
            tmax = imax(ctx, tmax,tu);
            tmaxe = imax(ctx, tmaxe,tu);
            if (tu > tl)
                nic += wt;
            else
                nex += wt;
        }
        else {
            tmax = imax(ctx, tmax,tl);
            nrc += wt;
        }
        if (s > 0)
            nlt += wt;

    }
    printf1(ctx, "Number of cases (weighted): %g\n",ws);
    printf1(ctx, "Right censored: %g\n",nrc);
    printf1(ctx, "Left truncated: %g\n",nlt);
    printf1(ctx, "Exact duration: %g\n",nex);
    printf1(ctx, "Interval-valued duration: %g\n\n",nic);
    printf1(ctx, "Minimal time value of left truncation: %d\n",smin);
    printf1(ctx, "Maximal time value of events: %d\n",tmaxe);
    printf1(ctx, "Maximal time value for calculations: %d\n",tmax);

    if (ws <= ctx->EPSI1) {
        printf1(ctx, "Number of cases is almost zero.\n");
        err = 0;
        goto DIPLEFin;
    }

    if (alloc_acx(ctx, tmax + 2))   
        goto DIPLEFin;
    if (alloc_acy(ctx, tmax + 2))   
        goto DIPLEFin;
    
    /* find mean distribution */

    wt = 1.0;
    for (i = 0; i < ctx->NOC; ++i) {

        if (ctx->WIVar >= 0)                     
            wt = get_data(ctx, ctx->WIVar,i) * ctx->WNorm;
        if (wt <= 0.0)
            continue;

        s  = (int)get_data(ctx, iis,i);
        tl = (int)get_data(ctx, iitl,i);
        tu = (int)get_data(ctx, iitu,i);
        d  = (int)get_data(ctx, iid,i);

        if (d == 0) {
            for (t = s; t <= tu; ++t)
                ctx->AcX[t] += wt;                
        }
        else if (tl == tu) {
            for (t = s; t <= tu; ++t)
                ctx->AcX[t] += wt;                
            ctx->AcY[tu] += wt;
        }
        else if (tl < tu) {

            w = wt;
            ilen = tu - tl + 1;
            tmp = wt / (double)ilen;

            for (t = s; t < tl; ++t)
                ctx->AcX[t] += wt;                

            for (t = tl; t <= tu; ++t) {
                ctx->AcY[t] += tmp;
                ctx->AcX[t] += w;
                w -= tmp;
            }
        }
    }
    nrec = 0;
    tmp = 1.0;
    for (t = smin; t <= tmax; ++t) {
        if (fabs(ctx->AcX[t]) <= ctx->EPSI1)   
            r = 0.0;
        else
            r = ctx->AcY[t] / ctx->AcX[t];

        fprintf(ctx->PMF1d,"%6d ",t);
        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[t]);
        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[t]);
        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,r);
        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,tmp);
        fprintf(ctx->PMF1d,"\n");
#ifdef TDA_R_PACKAGE
        {
            double erow[5];
            erow[0] = (double)t;
            erow[1] = ctx->AcX[t];
            erow[2] = ctx->AcY[t];
            erow[3] = r;
            erow[4] = tmp;
            tda_export_row(ctx, "diple.table", erow, 5);
        }
#endif
        nrec++;

        tmp *= (1.0 - r);
    }
    printf1(ctx, "\n%d records written to: %s\n",nrec,ctx->PMF1dName);
    err = 0;

DIPLEFin:
    p_clean(ctx);
    return(err);
}


