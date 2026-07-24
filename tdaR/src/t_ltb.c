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
#include "tda_context.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_ltb.c                                                    */

int ltb(TDAContext *ctx);
int lifetab(TDAContext *ctx);
void prltb1_head(TDAContext *ctx, int sn,int org,int grp,int n,double w,char *df);
void prltb1(TDAContext *ctx, int l,double ne,double zz,double rs,double *ee,char *df);
void prltb2_head(TDAContext *ctx, int ndes,char *df);
void prltb2(TDAContext *ctx, int l, double surv, double sterr, double *dens, double derr, double *rate, double rerr, int ndes, char *df);
void prstar(TDAContext *ctx);

/* ------------------------------------------------------------------------ */
/*  ltb()           Life table calculation. ltb command in CmdBuf.          */
/*                  Return 0 if OK, otherwise -1.                           */

int ltb(TDAContext *ctx)
{
    int err = -1;

    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Life table estimation. Current memory: %d bytes.\n",ctx->MemReq);

    if (ctx->EDAvail == 0) {
        p_err(ctx, -15,1);
        return(0);
    }
    if (parm(ctx, ctx->CmdBuf + 3,1,1))     /* get parameters */
        goto LTBFin;

    if (ctx->PMNTP == 0) {
        p_err(ctx, -17,1);
        goto LTBFin;
    }
    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 7,5);

    if (ctx->NSP > 0)  
        p_warn(ctx, -2,1);

    newline(ctx);
    err = lifetab(ctx);
    if (err) {
        if (err < 0)
            p_err(ctx, -2,1);
        else {
            printf1(ctx, "\nWarning: at least one starting time is not zero.\n");
            printf1(ctx, "Results are probably wrong.\n");
            err = -1;
        }
    }

LTBFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  lifetab()       Lifetable calculation.                                  */
/*                  Return 0 if OK, -1 if insufficient memory               */
/*                  Return 1 if not all starting times are zero.            */

int lifetab(TDAContext *ctx)
{
    register int k,l,grp;
    int err,err1,nt1,isn,iorg,ip,icase,sn,org,des,ndes,n0,ntab,spl,nspl;
    int a_ne,a_zz,a_rs,a_es,a_ee,a_subd,a_dens,a_rate,a_df;
    double tf,ts,wt,tmp,tmp1,pl,ql,med,surv,surv1,rate1 = 0.0,sterr,derr,rerr,sum,tt;
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
    nt1 = ctx->PMNTP + 1;

    if (!(ne = (double *)calloc((size_t)(nt1),sizeof(double))))  
        goto LTBFin;
    a_ne = 1;
    memrq(ctx, nt1,sizeof(double));

    if (!(zz = (double *)calloc((size_t)(nt1),sizeof(double))))  
        goto LTBFin;
    a_zz = 1;
    memrq(ctx, nt1,sizeof(double));

    if (!(rs = (double *)calloc((size_t)(nt1),sizeof(double))))  
        goto LTBFin;
    a_rs = 1;
    memrq(ctx, nt1,sizeof(double));

    if (!(es = (double *)calloc((size_t)(nt1),sizeof(double))))  
        goto LTBFin;
    a_es = 1;
    memrq(ctx, nt1,sizeof(double));
   
    if (!(ee = (double *)calloc((size_t)(ctx->PMNTP) * (size_t)(ctx->MaxDes1),sizeof(double))))  
        goto LTBFin;
    a_ee = 1;
    memrq(ctx, ctx->PMNTP * ctx->MaxDes1,sizeof(double));
  
    if (!(subd = (double *)calloc((size_t)(ctx->MaxDes1 + 1),sizeof(double))))  
        goto LTBFin;
    a_subd = 1;
    memrq(ctx, ctx->MaxDes1 + 1,sizeof(double));
  
    if (!(dens = (double *)calloc((size_t)(ctx->MaxDes1 + 1),sizeof(double))))  
        goto LTBFin;
    a_dens = 1;
    memrq(ctx, ctx->MaxDes1 + 1,sizeof(double));
  
    if (!(rate = (double *)calloc((size_t)(ctx->MaxDes1 + 1),sizeof(double))))  
        goto LTBFin;
    a_rate = 1;
    memrq(ctx, ctx->MaxDes1 + 1,sizeof(double));
  
    if (!(df = (char *)calloc((size_t)(ctx->MaxDes1),sizeof(char))))  
        goto LTBFin;
    a_df = 1;
    memrq(ctx, ctx->MaxDes1,sizeof(char));

    /* loop over all sn,org combinations */

    printf1(ctx, "SN  Org  Group ");
    if (ctx->PM1NV > 0)  
        prnchar(ctx, ' ',ctx->VNameLen + 1,0);
    printf1(ctx, "  Median  Dest.States  Episodes  Weighted\n");
    if (ctx->PM1NV > 0)
        prnchar(ctx, '-',ctx->VNameLen + 1,0);
    prnchar(ctx, '-',56,1);

    wt = 1.0;

    for (isn = 1; isn <= ctx->MaxSnn; ++isn) {      /* loop over all sn */

        for (iorg = 0; iorg <= ctx->MaxOrg; ++iorg) {    /* loop over all origins */
       
            ip = ctx->TranPtr[(isn - 1) * ctx->MaxOrg1 + iorg];
            if (ip < 0)
                continue;

            if (ctx->PM1NV) grp = 0; else grp = -1;

            for (; grp < ctx->PM1NV; ++grp) {               /* loop for all groups */
    
                n0 = 0;
                ne[0] = 0.0;
                for (l = 0; l < ctx->PMNTP; ++l)
                    es[l] = zz[l] = 0.0;
                k = ctx->MaxDes1 * ctx->PMNTP;
                for (l = 0; l < k; ++l)
                    ee[l] = 0.0;
                for (k = 0; k <= ctx->MaxDes; ++k)
                    df[k] = '\0';
    
                /* Loop over all episodes. Calculation of the number of     */
                /* episodes with events in es, differentiated with respect  */
                /* to destination states in ee.                             */

                get_spell(ctx, 1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
                while (get_spell(ctx, 0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

                    if (fabs(ts) > ctx->EPSI1)
                        err1 = 1;

                    if (sn != isn || org != iorg)               
                        continue;

                    if (grp >= 0) {
                        if (fabs(get_data(ctx, ctx->PM1VIdx[grp],icase)) <= ctx->EPSI1)
                            continue;
                    }

                    if (ctx->WIVar >= 0)             /* get weights */
                        wt = get_data(ctx, ctx->WIVar,icase) * ctx->WNorm;
   
                    n0++;
                    ne[0] += wt;
                    tf -= ts;       /* always durations !! */

                    for (l = 1; l < ctx->PMNTP; ++l) {
                        if (tf < ctx->PMTP[l])  
                            break;
                    }
                    l--;    /* index of time period */

                    if (des == org)     /* censored cases */
                        zz[l] += wt;
                    else {              /* episodes with an event */
                        df[des] = 1;
                        es[l] += wt;
                        ee[l * ctx->MaxDes1 + des] += wt;
                    }
                }

                /* Calculate number of destination states in ndes */

                ndes = 0;
                for (k = 0; k <= ctx->MaxDes; ++k) {
                    if (df[k])
                        ndes++;
                }
                if (n0 == 0 || ndes == 0)
                    continue;

                /* Calculate number entering the intervals in ne[] */

                for (l = 1; l < ctx->PMNTP; ++l)  
                    ne[l] = ne[l - 1] - es[l - 1] - zz[l - 1];
                
                /*  Print headers of requested tables.                      */
   
                prltb1_head(ctx, sn,iorg,grp,n0,ne[0],df);
    
                /*  Calculate the risk set and print body of first part of  */
                /*  the life table.                                         */

                for (l = 0; l < ctx->PMNTP; ++l) {
                    if (ne[l] <= 0.0)
                        break;
                    rs[l] = ne[l] - ctx->PMCFrac * zz[l];
                    if (rs[l] < 0.0)
                        break;
   
                    prltb1(ctx, l,ne[l],zz[l],rs[l],ee + l * ctx->MaxDes1,df);
                }

                /* Print header of second life table part */
    
                prltb2_head(ctx, ndes,df);
   
                /*  Calculate overall survivor function surv, with standard */
                /*  error sterr. Destination-specific densities and rates   */
                /*  calculated in the arrays dens and rate. If only one     */
                /*  destination state, standard errors of the density and   */
                /*  rate function are calculated in derr and rerr.          */

                /*  The median is always calculated by linear interpolation */
                /*  of the overall survivor function.                       */

                surv = 1.0;
                med = sum = 0.0;
                for (k = 0; k <= ctx->MaxDes; ++k)
                    subd[k] = 0.0;

                for (l = 0; l < ctx->PMNTP; ++l) {

                    if (rs[l] <= 0.0)
                        break;

                    sterr = derr = rerr = -1.0;

                    ql = es[l] / rs[l];
                    pl = 1.0 - ql;

                    /* stand error of survivor function */

                    if (sum >= 0.0)
                        sterr = surv * sqrt(sum);

                    if (l < ctx->PMNTP - 1) {

                        surv1 = surv * pl;      /* next survivor function */

                        tt = ctx->PMTP[l + 1] - ctx->PMTP[l]; /* length of time interval */

                        for (k = 0; k <= ctx->MaxDes; ++k) {

                            if (df[k]) {
                                tmp = subd[k] + surv * ee[l * ctx->MaxDes1 + k] / rs[l];
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
                            med = ctx->PMTP[l] + tt * (surv - 0.5) / tmp;
                     
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
     
                    prltb2(ctx, l,surv,sterr,dens,derr,rate,rerr,ndes,df);

                    surv = surv1;
                }       

                printf1(ctx, "%2d %4d  ",isn,iorg); 
                if (ctx->PM1NV == 0)  
                    printf1(ctx, "  --  ");
                else {
                    printf1(ctx, "%4d  ",grp + 1);
                    prn_vname(ctx, ctx->PM1VIdx[grp]);
                }

                if (med > 0.0) {
                    fprintf(ctx->PMFd,"\n# Median duration: %4.2f",med);
                    printf1(ctx, "%8.2f ",med);
                }
                else
                    printf1(ctx, "     --- ");
                printf1(ctx, "  %5d   %12d %9.2f\n",ndes,n0,ne[0]);
#ifdef TDA_R_PACKAGE
                /* the summary line the ltb table header announces;
                   accumulated across groups and flushed by
                   tda_export_flush_all() when the run ends */
                {
                    double srow[5];
                    srow[0] = (double)isn;
                    srow[1] = (double)iorg;
                    srow[2] = med > 0.0 ? med : (double)(NAN);
                    srow[3] = (double)n0;
                    srow[4] = ne[0];
                    tda_export_row(ctx, "ltb.summary", srow, 5);
                }
#endif
                ntab++;        
                fprintf(ctx->PMFd,"\n\n");
#ifdef TDA_R_PACKAGE
                /* one export per printed table, so groups come back as
                   ltb.risk / ltb.risk.2 / ... beside their blocks */
                tda_export_flush(ctx, "ltb.risk");
                tda_export_flush(ctx, "ltb.est");
#endif

            }   /* end of loop for all groups */
        }       /* end of loop for all origins */
    }           /* end of loop for all serial numbers */

    if (ctx->PM1NV > 0)
        prnchar(ctx, '-',ctx->VNameLen + 1,0);
    prnchar(ctx, '-',56,1);
    printf1(ctx, "%d table(s) written to: %s\n",ntab,ctx->PMFdName);
    err = err1;

LTBFin:
    if (a_ne) {
        free((char *)ne);
        memrq(ctx, -nt1,sizeof(double));
    }
    if (a_zz) {
        free((char *)zz);
        memrq(ctx, -nt1,sizeof(double));
    }
    if (a_rs) {
        free((char *)rs);
        memrq(ctx, -nt1,sizeof(double));
    }
    if (a_es) {
        free((char *)es);
        memrq(ctx, -nt1,sizeof(double));
    }
    if (a_ee) {
        free((char *)ee);
        memrq(ctx, -ctx->PMNTP * ctx->MaxDes1,sizeof(double));
    }
    if (a_subd) {
        free((char *)subd);
        memrq(ctx, -ctx->MaxDes1 - 1,sizeof(double));
    }
    if (a_dens) {
        free((char *)dens);
        memrq(ctx, -ctx->MaxDes1 - 1,sizeof(double));
    }
    if (a_rate) {
        free((char *)rate);
        memrq(ctx, -ctx->MaxDes1 - 1,sizeof(double));
    }
    if (a_df) {
        free((char *)df);
        memrq(ctx, -ctx->MaxDes1,sizeof(char));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prltb1_head.  Print header of first part of life table.                 */

void prltb1_head(TDAContext *ctx, int sn,int org,int grp,int n,double w,char *df)
{
    register int k;

    fprintf(ctx->PMFd,"# Life table. SN %d. Origin state %d. ",sn,org);
    if (grp >= 0)
        fprintf(ctx->PMFd,"\n# Group: %s",ctx->VName[ctx->PM1VIdx[grp]]);
    fprintf(ctx->PMFd,"\n# Cases: %d  weighted: %g\n",n,w);
    fprintf(ctx->PMFd,"\n# Start of          Number   Number   Exposed   ");
    for (k = 0; k <= ctx->MaxDes; ++k) {
        if (df[k]) {
            fprintf(ctx->PMFd,"  D-State %-3d",k);
            fprnchar(ctx, ctx->PMFd,' ',ctx->PMFmt1 - 5,0);
        }
    }
    fprintf(ctx->PMFd,"\n# Interval Midpoint Entering Censored   to Risk ");
    for (k = 0; k <= ctx->MaxDes; ++k) {
        if (df[k]) {
            fprintf(ctx->PMFd,"Events ");
            fprnchar(ctx, ctx->PMFd,' ',ctx->PMFmt1 - 5,0);
            fprintf(ctx->PMFd,"Prob  ");
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  prltb1  Print one line of first life table part.                        */

void prltb1(TDAContext *ctx, int l,double ne,double zz,double rs,double *ee,char *df)
{
    register int k;
#ifdef TDA_R_PACKAGE
    /* the same row the prt= file gets, staged for ltb.risk
       (CONTRIBUTING.md); the last interval's open midpoint and a
       zero risk set print as * and export as NaN */
    {
        double *erow = (double *)malloc((size_t)(5 + 2 * (ctx->MaxDes + 1))
                                        * sizeof(double));
        int ec = 0;
        if (erow == NULL)
            goto ltb1_noexp;
        erow[ec++] = ctx->PMTP[l];
        erow[ec++] = l < ctx->PMNTP - 1 ?
            (ctx->PMTP[l] + ctx->PMTP[l + 1]) / 2.0 : (double)NAN;
        erow[ec++] = ne;
        erow[ec++] = zz;
        erow[ec++] = rs;
        for (k = 0; k <= ctx->MaxDes; ++k) {
            if (df[k]) {
                erow[ec++] = ee[k];
                erow[ec++] = rs > 0.0 ? ee[k] / rs : (double)NAN;
            }
        }
        tda_export_row(ctx, "ltb.risk", erow, ec);
        free(erow);
ltb1_noexp: ;
    }
#endif
     
    fprintf(ctx->PMFd,"\n%10.2f",ctx->PMTP[l]);
    if (l < ctx->PMNTP - 1)
        fprintf(ctx->PMFd," %8.2f",(ctx->PMTP[l] + ctx->PMTP[l + 1]) / 2.0);
    else
        fprintf(ctx->PMFd,"        *");
    fprintf(ctx->PMFd," %8d %8d %9.1f ",(int)ne,(int)zz,rs);
    for (k = 0; k <= ctx->MaxDes; ++k) {
        if (df[k]) {
            fprintf(ctx->PMFd,"%6d ",(int)ee[k]);
            if (rs > 0.0)
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ee[k] / rs);
            else
                prstar(ctx);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  prltb2_head().  Print header of second part of life table.              */

void prltb2_head(TDAContext *ctx, int ndes,char *df)
{
    register int k;

    fprintf(ctx->PMFd,"\n\n# Start of          ");
    fprnchar(ctx, ctx->PMFd,' ',ctx->PMFmt1 - 8,0); fprintf(ctx->PMFd,"     Survivor");
    fprnchar(ctx, ctx->PMFd,' ',ctx->PMFmt1 - 3,0);
  
    for (k = 0; k <= ctx->MaxDes; ++k) {
        if (df[k]) {
            fprnchar(ctx, ctx->PMFd,' ',2 * ctx->PMFmt1 - 13,0);
            fprintf(ctx->PMFd,"  D-State %-3d  ",k);
            if (ndes == 1) {
                fprnchar(ctx, ctx->PMFd,' ',2 * ctx->PMFmt1 - 13,0);
                fprintf(ctx->PMFd,"  D-State %-3d  ",k);
            }
        }
    }
    fprintf(ctx->PMFd,"\n# Interval Midpoint ");
    fprnchar(ctx, ctx->PMFd,' ',ctx->PMFmt1 - 7,0); fprintf(ctx->PMFd,"Function");
    fprnchar(ctx, ctx->PMFd,' ',ctx->PMFmt1 - 5,0); fprintf(ctx->PMFd,"Error ");
    for (k = 0; k <= ctx->MaxDes; ++k) {
        if (df[k]) {
            fprnchar(ctx, ctx->PMFd,' ',ctx->PMFmt1 - 7,0); fprintf(ctx->PMFd,"Density ");
            fprnchar(ctx, ctx->PMFd,' ',ctx->PMFmt1 - 5,0);
            if (ndes == 1) {
                fprintf(ctx->PMFd,"Error "); fprnchar(ctx, ctx->PMFd,' ',ctx->PMFmt1 - 7,0); 
                fprintf(ctx->PMFd,"   Rate "); fprnchar(ctx, ctx->PMFd,' ',ctx->PMFmt1 - 5,0);
                fprintf(ctx->PMFd,"Error");
            }
            else
                fprintf(ctx->PMFd," Rate ");
        }
    }
}    

/* ------------------------------------------------------------------------ */
/*  prltb2  Print one line of second life table part.                       */

void prltb2(TDAContext *ctx, int l, double surv, double sterr, double *dens, double derr, double *rate, double rerr, int ndes, char *df)
{
    register int k;
#ifdef TDA_R_PACKAGE
    /* ltb.est: time, midpoint, survivor, its error, then per printed
       destination the same value-or-* sequence the file shows (density,
       and with a single destination also its error, the rate and the
       rate's error); * exports as NaN */
    {
        double *erow = (double *)malloc((size_t)(4 + 4 * (ctx->MaxDes + 1))
                                        * sizeof(double));
        int ec = 0;
        if (erow == NULL)
            goto ltb2_noexp;
        erow[ec++] = ctx->PMTP[l];
        erow[ec++] = l < ctx->PMNTP - 1 ?
            (ctx->PMTP[l] + ctx->PMTP[l + 1]) / 2.0 : (double)NAN;
        erow[ec++] = surv;
        erow[ec++] = sterr >= 0.0 ? sterr : (double)NAN;
        for (k = 0; k <= ctx->MaxDes; ++k) {
            if (df[k]) {
                erow[ec++] = dens[k] >= 0.0 ? dens[k] : (double)NAN;
                if (ndes == 1) {
                    erow[ec++] = derr > 0.0 ? derr : (double)NAN;
                    erow[ec++] = rate[k] >= 0.0 ? rate[k] : (double)NAN;
                    erow[ec++] = rerr > 0.0 ? rerr : (double)NAN;
                }
            }
        }
        tda_export_row(ctx, "ltb.est", erow, ec);
        free(erow);
ltb2_noexp: ;
    }
#endif

    fprintf(ctx->PMFd,"\n%10.2f",ctx->PMTP[l]);
    if (l < ctx->PMNTP - 1)
        fprintf(ctx->PMFd," %8.2f ",(ctx->PMTP[l] + ctx->PMTP[l + 1]) / 2.0);
    else
        fprintf(ctx->PMFd,"        * ");

    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,surv);
    if (sterr >= 0.0)
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,sterr);
    else  
        prstar(ctx);

    for (k = 0; k <= ctx->MaxDes; ++k) {
        if (df[k]) {
            if (dens[k] >= 0.0)
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,dens[k]);
            else
                prstar(ctx);
  
            if (ndes == 1) {
                if (derr > 0.0)
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,derr);
                else
                    prstar(ctx);

                if (rate[k] >= 0.0)
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,rate[k]);
                else
                    prstar(ctx);

                if (rerr > 0.0)
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,rerr);
                else
                    prstar(ctx);
            }
            else {
                if (rate[k] >= 0.0)
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,rate[k]);
                else
                    prstar(ctx);
            }
        }
    }
}

void prstar(TDAContext *ctx)
{
    fprnchar(ctx, ctx->PMFd,' ',ctx->PMFmt1 - 2,0);
    fprintf(ctx->PMFd," * ");
}



