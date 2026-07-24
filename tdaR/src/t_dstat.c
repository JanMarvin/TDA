/****************************************************************************/
/*  t_dstat                                                                 */
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
#include "t_ds.h"
#include "t_sort.h"
#include "t_freq.h"
#include "t_gf.h"
#include "t_cdf.h"
#include "t_matf.h"
#include "t_svd.h"
#include "tda_context.h"

/*  functions in t_dstat.c */

int dstat(TDAContext *ctx);
int quant(TDAContext *ctx);
int mfreq(TDAContext *ctx, int typ);
int ct_stat(TDAContext *ctx, int opt,int m,int *value,int *freq,int *bf,int widx,double *wfreq, int ridx,int cidx);
int pcov(TDAContext *ctx, int typ);
int atab(TDAContext *ctx);
int dma(TDAContext *ctx);
void dma_prn1(TDAContext *ctx, int m,int n,double *x,int p);
void dma_prn2(TDAContext *ctx, int m,int n,double *x,int *vidx);
void dma_prn3(TDAContext *ctx, FILE *fd,int m,int n,double *x,int p);
void dma_scal(TDAContext *ctx, int n,double *d,double *x);
void dma_stand(TDAContext *ctx, int n,int m,double *x,int opt);
void dma_corr(TDAContext *ctx, int n,int m,double *x,double *a);
void dma_pcf(TDAContext *ctx, int nr,double *rx,int nc,double *cx);

/* ------------------------------------------------------------------------ */
/*  dstat()     descriptive statistics.                                     */
/*                                                                          */
/*              dstat(                                                      */
/*                  grp=...,                    group variables             */
/*                  fmt=...,                    print format, def. 10.4     */
/*                  df=...,                     output file                 */
/*                  mppar=...,                  create matrix               */
/*              ) [=varlist];                   optional varlist            */  
/*                                                                          */
/*              Return 0 if OK, otherwise -1.                               */

int dstat(TDAContext *ctx)
{
    register int i = 0,j = 0;
    int err = 0,ig = 0,ign = 0,ix = 0,n = 0,l = 0,w = 0,nn = 0,nrow = 0;
    double xmin = 0.0,xmax = 0.0,xm = 0.0,xsd = 0.0,xsum = 0.0,wsum = 0.0,wt = 0.0,tmp = 0.0;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Descriptive statistics. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,4,0))     /* get parameters */
        goto DSFin;

    if (ctx->PMNV == 0) {        /* use all variables */

        if (pm_valloc(ctx, ctx->NVAR))  
            goto DSFin;
             
        w = j = 0;
        i = ctx->VIFirst;
        while (i >= 0) {
            ctx->PMVIdx[j++] = (short)(i);
            if (ctx->VTyp[i] == 5)
                w = 1;
            i = ctx->VNxt[i];
        }
        if (w)
            p_warn(ctx, -3,1);
    }   
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    if (alloc_acx(ctx, ctx->NOC + 1))
        goto DSFin;

    if (ctx->SVEFlg) {           /* info about case selection */
        prn_sve(ctx);
        if (ctx->WIVar >= 0) {
            printf1(ctx, "\nCannot use case weights with sel option.\n");
            goto DSFin;
        }
    }
    prn_cwt(ctx);              /* info about case weights */

    if (ctx->PMMPParDef == 1) {      /* use AcTmp for mppar */
        n = nn = ctx->PMNV * 6;
        nrow = ctx->PMNV;
        if (ctx->PM1NV > 0) {
            n *= ctx->PM1NV;
            nrow *= ctx->PM1NV;
        }
        if (alloc_actmp(ctx, n + 1))
            goto DSFin;
    }
    ig = -1;
    ign = 0;

DSNXT:
    if (ctx->PM1NV > 0) {
        ig = ctx->PM1VIdx[ign];
        printf1(ctx, "\nGroup: %s\n",ctx->VName[ig]);

        if (ctx->PMF1Def)  
            fprintf(ctx->PMF1d,"# Group: %s\n",ctx->VName[ig]);

    }         
    newline(ctx);
    prn_hvar(ctx);
    prn_hlabel(ctx);
    prnchar(ctx, ' ',ctx->PMFmt1 -  8,0); printf1(ctx, " Minimum"); 
    prnchar(ctx, ' ',ctx->PMFmt1 -  8,0); printf1(ctx, "  Maximum");
    prnchar(ctx, ' ',ctx->PMFmt1 -  8,0); printf1(ctx, "     Mean");
    prnchar(ctx, ' ',ctx->PMFmt1 -  8,0); printf1(ctx, " Std.Dev.");
    printf1(ctx, "    Sum of values\n");
       
    l = ctx->VNameLen + 18;
    if (ctx->VLabelLen > 0)
        l += ctx->VLabelLen + 1;
    prnchar(ctx, '-',l + 4 * (ctx->PMFmt1 + 1) - 1,1);

    for (j = 0; j < ctx->PMNV; ++j) {

        ix = ctx->PMVIdx[j];
        n = getd1(ctx, ctx->AcX,ix,ig,0);

        xmin = xmax = xm = xsd = xsum = 0.0;
        if (n > 0) {
            xmin = ctx->DBLMAX;
            xmax = -xmin;  

            wt = 1.0;
            if (ctx->WIVar >= 0)
                wsum = 0.0;
            else
                wsum = (double)n;

            for (i = 0; i < n; ++i) {
                tmp = ctx->AcX[i];
                if (ctx->WIVar >= 0) {       
                    wt = get_data(ctx, ctx->WIVar,i) * ctx->WNorm;
                    wsum += wt;
                }
                xsum += tmp * wt;
                if (xmax < tmp)
                    xmax = tmp;
                if (xmin > tmp)
                    xmin = tmp;
            }
            if (wsum > 1.0) {
                xm = xsum / wsum;      

                for (i = 0; i < n; ++i) {
                    tmp = ctx->AcX[i];

                    if (ctx->WIVar >= 0)         
                        wt = get_data(ctx, ctx->WIVar,i) * ctx->WNorm;
                    tmp -= xm;
                    xsd += tmp * tmp * wt;
                }
                xsd /= (wsum - 1.0);
                if (xsd > 0.0)
                    xsd = sqrt(xsd);
                else
                    xsd = 0.0;
            }
        }
        prn_vname(ctx, ix);
        prn_vlabel(ctx, ix);
        rt_printf1_d(ctx, ctx->PMFmtS,xmin);
        rt_printf1_d(ctx, ctx->PMFmtS,xmax);
        rt_printf1_d(ctx, ctx->PMFmtS,xm);
        rt_printf1_d(ctx, ctx->PMFmtS,xsd);
        printf1(ctx, "%16.4f\n",xsum);
#ifdef TDA_R_PACKAGE
        /* the same row the table prints, and its variable name, staged
           for dstat.stats / dstat.names (CONTRIBUTING.md) */
        {
            double erow[5];
            erow[0] = xmin; erow[1] = xmax; erow[2] = xm;
            erow[3] = xsd; erow[4] = xsum;
            tda_export_row(ctx, "dstat.stats", erow, 5);
            tda_export_str_row(ctx, "dstat.names", ctx->VName[ix]);
        }
#endif

        if (ctx->PMF1Def) {
            fprintf(ctx->PMF1d,"%4d ",j + 1);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,xmin);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,xmax);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,xm);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,xsd);
            fprintf(ctx->PMF1d,"%16.4f\n",xsum);
        }
        if (ctx->PMMPParDef == 1) {
            ctx->AcTmp[ign * nn + j * 6 + 1] = (double)ign;
            ctx->AcTmp[ign * nn + j * 6 + 2] = xmin;
            ctx->AcTmp[ign * nn + j * 6 + 3] = xmax;
            ctx->AcTmp[ign * nn + j * 6 + 4] = xm;
            ctx->AcTmp[ign * nn + j * 6 + 5] = xsd;
            ctx->AcTmp[ign * nn + j * 6 + 6] = xsum;
        }
    }
    if (++ign < ctx->PM1NV)  
        goto DSNXT;
       
    if (ctx->PMMPParDef == 1) {
        mp_putmpar(ctx, nrow,6,ctx->AcTmp);
        newline(ctx);
        mp_info(ctx);
    }
    err = 0;

DSFin:
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "dstat.stats");
    tda_export_str_flush(ctx, "dstat.names");
#endif
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  quant()     Quantiles.                                                  */
/*                                                                          */
/*              quant(                                                      */
/*                  fmt=...,                    print format, def. 10.4     */
/*                  df=...,                     output file                 */
/*                  mppar=...,                  create matrix               */
/*              ) =varlist;                     varlist required            */  
/*                                                                          */
/*              Return 0 if OK, otherwise -1.                               */

int quant(TDAContext *ctx)
{
    register int i;
    int err,n,iv,ivv,nqx,len;
    double tmp,qx[11];
#ifdef TDA_R_PACKAGE
    double equant[11];
#endif

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Quantiles. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,4,1))     /* get parameters */
        goto QUANTFin;

    if (alloc_acx(ctx, ctx->NOC))     /* allocate AcX */
        goto QUANTFin;

    if (ctx->SVEFlg)             /* info about case selection */
        prn_sve(ctx);

    if (ctx->PMMPParDef == 1) {      /* use AcTmp for mppar */
        if (alloc_actmp(ctx, ctx->PMNV * 11 + 1))
            goto QUANTFin;
    }

    nqx = 11;
    qx[ 0] = 0.1 ;
    qx[ 1] = 0.2 ;
    qx[ 2] = 0.25;
    qx[ 3] = 0.3 ;
    qx[ 4] = 0.4 ;
    qx[ 5] = 0.5 ;
    qx[ 6] = 0.6 ;
    qx[ 7] = 0.7 ;
    qx[ 8] = 0.75;
    qx[ 9] = 0.8 ;
    qx[10] = 0.9 ;

    if (ctx->PMFmtF == 0)  
        pmfmt(ctx, 7,2);
      
    len = ctx->VNameLen;
    if (len < ctx->PMFmt1)
        len = ctx->PMFmt1;

    printf1(ctx, "\nVariable ");            
    prnchar(ctx, ' ',len - 8,0);
    for (i = 0; i < nqx; ++i)  
        rt_printf1_d(ctx, ctx->PMFmtS,qx[i]);
    newline(ctx);
    prnchar(ctx, '-',len + nqx * (ctx->PMFmt1 + 1),1);

    for (iv = 0; iv < ctx->PMNV; ++iv) {         /* for all variables */ 

        ivv = ctx->PMVIdx[iv];

        n = getd1(ctx, ctx->AcX,ivv,-1,0);

        printf1(ctx, "%s ",ctx->VName[ivv]);
        prnchar(ctx, ' ',len - (int)strlen(ctx->VName[ivv]),0);

        if (n == 0) {
            printf1(ctx, " number of cases is zero.\n");
            continue;       
        }
        if (sortd(ctx, n,ctx->AcX,0))
            goto QUANTFin;

        for (i = 0; i < nqx; ++i) {
            tmp = quantf(ctx, n,ctx->AcX,qx[i]);
            rt_printf1_d(ctx, ctx->PMFmtS,tmp);
            if (ctx->PMF1Def)
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,tmp);

            if (ctx->PMMPParDef == 1)  
                ctx->AcTmp[iv * 11 + i + 1] = tmp;          
#ifdef TDA_R_PACKAGE
            equant[i] = tmp;
#endif
        }
#ifdef TDA_R_PACKAGE
        tda_export_row(ctx, "quant.table", equant, nqx);
        tda_export_str_row(ctx, "quant.names", ctx->VName[ivv]);
#endif
        newline(ctx);
        if (ctx->PMF1Def)
            fprintf(ctx->PMF1d,"\n");
    }
    if (ctx->PMMPParDef == 1) {
        mp_putmpar(ctx, ctx->PMNV,11,ctx->AcTmp);
        newline(ctx);
        mp_info(ctx);
    }
    err = 0;

QUANTFin:
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "quant.table");
    tda_export_str_flush(ctx, "quant.names");
#endif
    p_clean(ctx); 
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  mfreq(typ)  Process CmdDef[idx] with syntax                             */
/*              typ = 0 : freq  (joint freq distribution)                   */
/*              typ = 1 : freq1 (separate freq distributions)               */
/*              typ = 2 : freq2 (2-dim table, optional with cont. meas.)    */
/*                                                                          */
/*              Return 0 if successful, otherwise -1.                       */

int mfreq(TDAContext *ctx, int typ)
{
    register int i,j,k,l;
    int err,nvv,vl,len,iv,kv,r;
    double f,sf,tmp,stmp,sum;
#ifdef TDA_R_PACKAGE
    double *erow = NULL;
    int ecol = 0;
#endif

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Frequency tables. Current memory: %d bytes.\n",ctx->MemReq);
    i = 0;
    if (typ > 0)
        i = 1;

    if (parm(ctx, ctx->CmdBuf + 4 + i,4,1))     /* get parameters */
        goto MFFin;
  
    if (typ == 2 && ctx->PMNV != 2) {
        printf1(ctx, "Error: freq2 command needs exactly two variables.\n");
        goto MFFin;
    }

    printf1(ctx, "Maximum number of categories: %d\n",ctx->PMMaxCat);

    /* allocate memory */
   
    if (alloc_ack(ctx, ctx->PMMaxCat * ctx->PMNV + 1))     /* value */
        goto MFFin;

    if (alloc_acn(ctx, ctx->PMMaxCat + 1))            /* freq */
        goto MFFin;

    if (alloc_aci(ctx, ctx->PMMaxCat + 1))            /* bf pointer */
        goto MFFin;

    if (ctx->WIVar >= 0) {   /* weight command active */

        if (alloc_acw(ctx, ctx->PMMaxCat + 1))    
            goto MFFin;
        prn_cwt(ctx);                  /* case weight info */
    }
             
    /* prepare print format */

    if (ctx->PMFmtF == 0)
        ctx->PMFmt1 = 3;

    vl = get_vnl(ctx, ctx->PMNV,ctx->PMVIdx);
    if (vl < ctx->PMFmt1)
        vl = ctx->PMFmt1;

    if (vl > ctx->PMFmt1)
        pmfmt(ctx, vl,0);
    else if (ctx->PMFmtF == 0)
        pmfmt(ctx, ctx->PMFmt1,0);

    for (i = 0; i < ctx->PMNV; ++i) {      /* make all distributions */

        printf1(ctx, "\nFrequency distribution for variable(s): %s",ctx->VName[ctx->PMVIdx[i]]);

        if (typ == 1) {
            r = cfreq1(ctx, 1,ctx->PMVIdx + i,ctx->PMMaxCat,ctx->AcK,ctx->AcN,ctx->AcI,ctx->WIVar,ctx->AcW);
        }
        else {
            for (j = 1; j < ctx->PMNV; ++j)  
                printf1(ctx, ",%s",ctx->VName[ctx->PMVIdx[j]]);
                     
            r = cfreq1(ctx, ctx->PMNV,ctx->PMVIdx,ctx->PMMaxCat,ctx->AcK,ctx->AcN,ctx->AcI,ctx->WIVar,ctx->AcW);

        }
        newline(ctx);

        if (r <= 0) {
            if (r == -1)
                printf1(ctx, "Max number of categories is too small.\n");
            else if (r == -2)
                printf1(ctx, "Stack size is too small.\n");
            else if (r == -3)
                printf1(ctx, "Insufficient memory for frequency distribution.\n");
            else 
                printf1(ctx, "Error %d in frequency distribution.\n",r);
            goto MFFin;
        }
        printf1(ctx, "Number of categories: %d\n",r);
        sum = 0.0;
        for (j = 1; j <= r; ++j) {
            if (ctx->WIVar >= 0)
                sum += ctx->AcW[j];
            else 
                sum += (double)ctx->AcN[j];
        }
        if (sum == ctx->EPSI)  
            printf1(ctx, "Sum of (weighted) frequencies is zero.\n");
        else if (typ < 2) {

            len = 42;
            printf1(ctx, "\nIndex ");
            if (typ == 1) {
                nvv = 1;
                iv = ctx->PMVIdx[i];
                prnchar(ctx, ' ',vl - (int)strlen(ctx->VName[iv]),0);
                printf1(ctx, " %s",ctx->VName[iv]);
                len += vl + 1;
            }
            else {
                nvv = ctx->PMNV;
                for (j = 0; j < ctx->PMNV; ++j) {
                    iv = ctx->PMVIdx[j];
                    prnchar(ctx, ' ',vl - (int)strlen(ctx->VName[iv]),0);
                    printf1(ctx, " %s",ctx->VName[iv]);
                    len += vl + 1;
                }
            }
            printf1(ctx, "  Frequency   Pct   Cumulated   Pct\n");
            prnchar(ctx, '-',len,1);

            sf = stmp = 0.0;
#ifdef TDA_R_PACKAGE
            ecol = nvv + 5;
            free(erow);
            erow = (double *)malloc((size_t)ecol * sizeof(double));
#endif

            for (j = 1; j <= r; ++j) {
                printf1(ctx, "%5d  ",j);
                if (ctx->PMF1Def)
                    fprintf(ctx->PMF1d,"%5d  ",j);

                k = ctx->AcI[j];
                for (l = 1; l <= nvv; ++l) {
                    kv = ctx->AcK[(k - 1) * nvv + l];
                    rt_printf1_d(ctx, ctx->PMFmtS,(double)kv);
                    if (ctx->PMF1Def)
                        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,(double)kv);
#ifdef TDA_R_PACKAGE
                    if (erow != NULL)
                        erow[l] = (double)kv;
#endif
                }
                if (ctx->WIVar >= 0)
                    f = ctx->AcW[k];
                else 
                    f = (double)ctx->AcN[k];

                tmp = 100.0 * f / sum;
                sf += f;
                stmp += tmp;

                printf1(ctx, "%10.2f %6.2f %10.2f %6.2f\n",f,tmp,sf,stmp);
                if (ctx->PMF1Def)
                    fprintf(ctx->PMF1d,"%10.2f %6.2f %10.2f %6.2f\n",f,tmp,sf,stmp);
#ifdef TDA_R_PACKAGE
                if (erow != NULL) {
                    erow[0] = (double)j;
                    erow[nvv + 1] = f;
                    erow[nvv + 2] = tmp;
                    erow[nvv + 3] = sf;
                    erow[nvv + 4] = stmp;
                    tda_export_row(ctx, "freq.table", erow, ecol);
                }
#endif
            }
#ifdef TDA_R_PACKAGE
            tda_export_flush(ctx, "freq.table");
#endif
            prnchar(ctx, '-',len,1);
            printf1(ctx, "Sum   ");
            prnchar(ctx, ' ',len - 42,0);
            printf1(ctx, " %10.2f %6.2f\n",sum,100.0);
            if (ctx->PMF1Def)  
                printf1(ctx, "\nTable written to: %s\n",ctx->PMF1dName);
        }
        else {          /* freq2 command */

            ct_stat(ctx, (int)ctx->PMSC,r,ctx->AcK,ctx->AcN,ctx->AcI,ctx->WIVar,ctx->AcW,(int)ctx->PMVIdx[0],(int)ctx->PMVIdx[1]);
        }
        if (typ != 1)
            break;
    }
    err = 0;

MFFin:
#ifdef TDA_R_PACKAGE
    free(erow);
    tda_export_flush(ctx, "freq.table");
#endif
    p_clean(ctx);  
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ct_stat(opt,m,value,freq,bf,widx,wfreq,ridx,cidx)                       */
/*                                                                          */
/*  Calculate and print a two-dimensional table.                            */
/*  stored in value[], with m lines. bf is a pointer to sort value.         */
/*  If widx >= 0 the frequencies are in wfreq (double), otherwise the       */
/*  frequencies are in freq. ridx and cidx are the internal variable        */
/*  number of the row and column variables, respectively.                   */
/*                                                                          */
/*  If opt != 0, calculate contingency measures.                            */
/*  Return 0 if OK, -1 if error.                                            */

int ct_stat(TDAContext *ctx, int opt,int m,int *value,int *freq,int *bf,int widx,double *wfreq, int ridx,int cidx)
{
    register int i,j,k,jj,kk;
    int r,c,n,n1,n2,df,rc,err,nflag;
    double lr,s,e,a,a1,a2,b1,b2,sum,t,tmp;

    err = -1;

    if (alloc_acr(ctx, m + 1))
        goto CTFin;

    if (alloc_acs(ctx, m + 1))
        goto CTFin;

    if (alloc_acm(ctx, m + 1))
        goto CTFin;

    if (alloc_acj(ctx, m + 1))
        goto CTFin;


    /* sort categories for columns */

    for (i = 1; i <= m; ++i)
        ctx->AcS[i] = value[i * 2];  

    c = cfreq(ctx, 1,m,ctx->AcS,m,ctx->AcR,ctx->AcM,ctx->AcJ,0,&a,&tmp);
    if (c < 1) {
        printf1(ctx, "Error: can't sort categories.\n");
        goto CTFin;
    }
    for (i = 1; i <= c; ++i)
        ctx->AcS[i - 1] = ctx->AcR[ctx->AcJ[i]];

    /* sort categories for rows */

    r = 0;
    k = bf[1];
    n = value[(k - 1) * 2 + 1] - 1;
    for (i = 1; i <= m; ++i) {
        k = bf[i];
        j = value[(k - 1) * 2 + 1];
        if (j > n) {
            ctx->AcR[r] = j;
            n = j;
            r++;
        }
    }

    rc = r * c;

    if (alloc_actmp(ctx, rc))        /* used to store the table */
        goto CTFin;

    if (alloc_acx(ctx, r))           /* row sum */
        goto CTFin;

    if (alloc_acy(ctx, c))           /* column sum */
        goto CTFin;

    nflag = 0;
    sum = 0.0;
    for (i = 1; i <= m; ++i) {
        n = (bf[i] - 1) * 2;
        jj = value[++n];
        kk = value[++n];
        for (j = 0; j < r; ++j) {
            if (ctx->AcR[j] == jj)
                break;
        }
        for (k = 0; k < c; ++k) {
            if (ctx->AcS[k] == kk)
                break;
        }
        if (widx >= 0)
            tmp = wfreq[bf[i]];
        else  
            tmp = (double)freq[bf[i]];

        ctx->AcTmp[j * c + k] += tmp;                   
        sum += tmp;
        ctx->AcX[j] += tmp;
        ctx->AcY[k] += tmp;
        if (tmp < 0.0)
            nflag++;
    }

    /* print the table */

    printf1(ctx, "\nFrequency|\n");
    printf1(ctx, "Percent  |\n");
    printf1(ctx, "Row Pct  |\n");
    printf1(ctx, "Col Pct  |");
    for (k = 0; k < c; ++k)  
        printf1(ctx, "%8d |",ctx->AcS[k]);
    printf1(ctx, "   Total\n");
    prnchar(ctx, '-',(c + 2) * 10 - 2,1);
    for (j = 0; j < r; ++j) {
        printf1(ctx, "%8d |",ctx->AcR[j]);
        for (k = 0; k < c; ++k) {
            printf1(ctx, "%8.2f |",ctx->AcTmp[j * c + k]);
            if (ctx->PMF1Def)
                fprintf(ctx->PMF1d,"%8.2f ",ctx->AcTmp[j * c + k]);
        }
        if (ctx->PMF1Def)
            fprintf(ctx->PMF1d,"\n");
#ifdef TDA_R_PACKAGE
        tda_export_row(ctx, "freq2.table", ctx->AcTmp + j * c, c);
#endif

        printf1(ctx, "%8.2f\n         |",ctx->AcX[j]);
        for (k = 0; k < c; ++k) {
            if (sum != 0.0)
                tmp = 100.0 * ctx->AcTmp[j * c + k] / sum;
            else
                tmp = 0.0;
            printf1(ctx, "%8.2f |",tmp);
        }
        if (sum != 0.0)
            tmp = 100.0 * ctx->AcX[j] / sum;
        else
            tmp = 0.0;
        printf1(ctx, "%8.2f\n         |",tmp);
        for (k = 0; k < c; ++k) {
            if (ctx->AcX[j] != 0.0) 
                tmp = 100.0 * ctx->AcTmp[j * c + k] / ctx->AcX[j];
            else
                tmp = 0.0;
            printf1(ctx, "%8.2f |",tmp);
        }
        printf1(ctx, "\n         |");
        for (k = 0; k < c; ++k) {
            if (ctx->AcY[k] != 0.0)
                tmp = 100.0 * ctx->AcTmp[j * c + k] / ctx->AcY[k];
            else
                tmp = 0.0;
            printf1(ctx, "%8.2f |",tmp);
        }
        newline(ctx);
        prnchar(ctx, '-',(c + 2) * 10 - 2,1);
    }
    printf1(ctx, "Total     ");
    for (k = 0; k < c; ++k)  
        printf1(ctx, "%8.2f  ",ctx->AcY[k]);
    printf1(ctx, "%8.2f\n          ",sum);
    for (k = 0; k < c; ++k) {
        if (sum != 0.0)
            tmp = 100.0 * ctx->AcY[k] / sum;
        else
            tmp = 0.0;
        printf1(ctx, "%8.2f  ",tmp);
    }
    printf1(ctx, "  100.00\n");
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "freq2.table");
#endif
    if (opt == 0)
        goto CTFin;

    /* ----------------------------------------------------------------- */

    newline(ctx);
    prnchar(ctx, '-',LLEN,1);
    printf1(ctx, "Contingency measures for table: X (%s) x Y (%s)\n\n",
                                                   ctx->VName[ridx],ctx->VName[cidx]);
    if (r < 2) {
        printf1(ctx, "Table has less than two categories in: %s\n",ctx->VName[ridx]);
        goto CTFin;
    }
    if (c < 2) {
        printf1(ctx, "Table has less than two categories in: %s\n",ctx->VName[cidx]);
        goto CTFin;
    }
    if (nflag > 0) {
        printf1(ctx, "Error: table has %d negative elements.\n",nflag);
        goto CTFin;
    }
    if (sum <= ctx->EPSI1) {
        printf1(ctx, "Error: sum of table elements is zero or less.\n");
        goto CTFin;
    }

    /*  chi square */

    err = n1 = n2 = 0;
    lr = s = 0.0;
    for (j = 0; j < r; ++j) {
        for (k = 0; k < c; ++k) {
            t = ctx->AcTmp[j * c + k];
            e = ctx->AcX[j] * ctx->AcY[k] / sum;
            if (t < 5.0)
                n1++;
            if (e < 5.0)
                n2++;

            if (e <= ctx->EPSI1)
                err++;
            else {
                s += (t - e) * (t - e) / e;
                lr += t * rlog(ctx, t / e);
            }
        }
    }
    df = (r - 1) * (c - 1);
    a = 100.0 * (double)n1 / (double)rc;
    printf1(ctx, "Number of cells with less than 5 elements: %3d (%4.1f %%)\n",
                                                                       n1,a);
    a = 100.0 * (double)n2 / (double)rc;
    printf1(ctx, "Cells with expected frequency less than 5: %3d (%4.1f %%)\n",     
                                                                       n2,a);
    printf1(ctx, "Degrees of freedom: %d\n\n",df);
    if (err > 0) {
        printf1(ctx, "Can't calculate Chi-square.\n");
        err = -1;
        goto CTFin;
    }
    prn_sfmt(ctx, "Chi-square (Pearson)",35,ctx->PMTFmtS,s);
    a = 1.0 - cdchif(ctx, s,df);
    prn_sfmt(ctx, "Prob:",35,ctx->PMTFmtS,a);
    lr *= 2.0;
    prn_sfmt(ctx, "Chi-square (likelihood ratio)",35,ctx->PMTFmtS,lr);
    a =  1.0 - cdchif(ctx, lr,df);
    prn_sfmt(ctx, "Prob:",35,ctx->PMTFmtS,a);

    if (s > ctx->EPSI1) {
        a = sqrt(s / sum);
        newline(ctx);
        prn_sfmt(ctx, "Phi",35,ctx->PMTFmtS,a);
        n = r;
        if (n > c)
            n = c;
        if (n > 1) {
            a = sqrt(s / (sum * (double)(n - 1)));
            prn_sfmt(ctx, "Cramer's V",35,ctx->PMTFmtS,a);
        }
        if (s + sum != 0.0 && (tmp = s / (s + sum)) > ctx->EPSI1) {
            a = sqrt(tmp);
            prn_sfmt(ctx, "Contingency coefficient",35,ctx->PMTFmtS,a);
        }

        /* lamda's */

        a1 = 0.0;
        for (k = 0; k < c; ++k) {
            tmp = ctx->AcTmp[k];
            for (j = 1; j < r; ++j) {
                t = ctx->AcTmp[j * c + k];
                if (tmp < t)
                    tmp = t;
            }
            a1 += tmp;
        }
        b1 = ctx->AcX[0];
        for (j = 1; j < r; ++j) {
            if (b1 < ctx->AcX[j])
                b1 = ctx->AcX[j];
        }
        if (sum > b1) {
            s = (a1 - b1) / (sum - b1);
            if (s > 0.0)  
                prn_sfmt(ctx, "Lambda (X dependent)",35,ctx->PMTFmtS,s);
        }

        a2 = 0.0;
        for (j = 0; j < r; ++j) {
            tmp = ctx->AcTmp[j * c];
            for (k = 1; k < c; ++k) {
                t = ctx->AcTmp[j * c + k];
                if (tmp < t)
                    tmp = t;
            }
            a2 += tmp;
        }
        b2 = ctx->AcY[0];
        for (j = 1; j < c; ++j) {
            if (b2 < ctx->AcY[j])
                b2 = ctx->AcY[j];
        }
        if (sum > b2) {
            s = (a2 - b2) / (sum - b2);
            if (s > 0.0)  
                prn_sfmt(ctx, "Lambda (Y dependent)",35,ctx->PMTFmtS,s);
        }
        tmp = b1 + b2;
        if (2.0 * sum > tmp) {
            s = (a1 + a2 - tmp) / (2.0 * sum - tmp);        
            if (s > 0.0)  
                prn_sfmt(ctx, "Lambda (symmetric)",35,ctx->PMTFmtS,s);
        }

        /* Goodman - Kruskal Tau */

        a1 = 0.0;
        for (k = 0; k < c; ++k) {
            if (ctx->AcY[k] < ctx->EPSI1) {
                a1 = -1.0;
                break;
            }
            for (j = 0; j < r; ++j) {
                tmp = ctx->AcTmp[j * c + k];
                a1 += tmp * tmp / ctx->AcY[k];
            }
        }
        if (a1 > 0.0) {
            b1 = 0.0;
            for (j = 0; j < r; ++j)
                b1 += ctx->AcX[j] * ctx->AcX[j];

            tmp = sum * sum - b1;
            if (tmp > 0.0) {
                tmp = (sum * a1 - b1) / tmp;
                if (tmp >= 0.0)  
                    prn_sfmt(ctx, "Goodman/Kruskal Tau (X dependent)",35,ctx->PMTFmtS,tmp);
            }
        }
        a1 = 0.0;
        for (j = 0; j < r; ++j) {
            if (ctx->AcX[j] < ctx->EPSI1) {
                a1 = -1.0;
                break;
            }
            for (k = 0; k < c; ++k) {
                tmp = ctx->AcTmp[j * c + k];
                a1 += tmp * tmp / ctx->AcX[j];
            }
        }
        if (a1 > 0.0) {
            b1 = 0.0;
            for (k = 0; k < c; ++k)
                b1 += ctx->AcY[k] * ctx->AcY[k];

            tmp = sum * sum - b1;
            if (tmp > 0.0) {
                tmp = (sum * a1 - b1) / tmp;
                if (tmp >= 0.0)  
                    prn_sfmt(ctx, "Goodman/Kruskal Tau (Y dependent)",35,ctx->PMTFmtS,tmp);
            }
        }

        /* Uncertainty coefficient */

        a1 = 0.0;
        for (j = 0; j < r; ++j) {
            if (ctx->AcX[j] > 0.0) {
                tmp = ctx->AcX[j] / sum;
                a1 -= tmp * rlog(ctx, tmp);  
            }
        }
        b1 = 0.0;
        for (k = 0; k < c; ++k) {
            if (ctx->AcY[k] > 0.0) {
                tmp = ctx->AcY[k] / sum;
                b1 -= tmp * rlog(ctx, tmp);  
            }
        }
        s = 0.0;
        for (j = 0; j < r; ++j) {
            for (k = 0; k < c; ++k) {
                tmp = ctx->AcTmp[j * c + k];
                if (tmp > 0.0) {
                    tmp /= sum;
                    s -= tmp * rlog(ctx, tmp);  
                }
            }
        }
        if (fabs(a1) > ctx->EPSI1) {
            tmp = (a1 + b1 - s) / a1;
            if (tmp > 0.0)  
                prn_sfmt(ctx, "Uncertainty coeff (X dependent)",35,ctx->PMTFmtS,tmp);
        }
        if (fabs(b1) > ctx->EPSI1) {
            tmp = (a1 + b1 - s) / b1;
            if (tmp > 0.0)  
                prn_sfmt(ctx, "Uncertainty coeff (Y dependent)",35,ctx->PMTFmtS,tmp);
        }

        /* Kendall's tau, etc */

        a1 = b1 = 0.0;
        for (j = 0; j < r; ++j) {
            for (k = 0; k < c; ++k) {
                tmp = ctx->AcTmp[j * c + k];
                if (tmp > 0.0) {
                    t = 0.0;
                    for (jj = 0; jj < j; ++jj) {
                        for (kk = 0; kk < k; ++kk)  
                            t += ctx->AcTmp[jj * c + kk];
                    }
                    for (jj = j + 1; jj < r; ++jj) {
                        for (kk = k + 1; kk < c; ++kk)  
                            t += ctx->AcTmp[jj * c + kk];
                    }
                    a1 += tmp * t;

                    t = 0.0;
                    for (jj = 0; jj < j; ++jj) {
                        for (kk = k + 1; kk < c; ++kk)  
                            t += ctx->AcTmp[jj * c + kk];
                    }
                    for (jj = j + 1; jj < r; ++jj) {
                        for (kk = 0; kk < k; ++kk)  
                            t += ctx->AcTmp[jj * c + kk];
                    }
                    b1 += tmp * t;
                }
            }
        }
        a2 = sum * sum;
        for (j = 0; j < r; ++j)
            a2 -= ctx->AcX[j] * ctx->AcX[j];

        b2 = sum * sum;
        for (k = 0; k < c; ++k)
            b2 -= ctx->AcY[k] * ctx->AcY[k];

        tmp = a2 * b2;
        if (tmp > ctx->EPSI1) {
            s = (a1 - b1) / sqrt(tmp);
            if (s >= 0.0)  
                prn_sfmt(ctx, "Kendall's tau-b",35,ctx->PMTFmtS,s);
        }
        if (n > 1) {            /* n = min(r,c) */
            tmp = sum * sum * ((double)(n - 1));
            if (tmp > ctx->EPSI1) {
                s = (double)n * (a1 - b1) / tmp;
                if (s >= 0.0)  
                    prn_sfmt(ctx, "Kendall's tau-c",35,ctx->PMTFmtS,s);
            }
        }
        tmp = a1 + b1;
        if (tmp > ctx->EPSI1) {
            s = (a1 - b1) / tmp;
            if (s >= 0.0)  
                prn_sfmt(ctx, "Gamma",35,ctx->PMTFmtS,s);
        }
        if (b2 > ctx->EPSI1) {
            s = (a1 - b1) / b2;
            if (s >= 0.0)  
                prn_sfmt(ctx, "Somers' D (X dependent)",35,ctx->PMTFmtS,s);
        }
        if (a2 > ctx->EPSI1) {
            s = (a1 - b1) / a2;
            if (s >= 0.0)  
                prn_sfmt(ctx, "Somers' D (Y dependent)",35,ctx->PMTFmtS,s);
        }
    }
    err = 0;

CTFin:
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  pcov(typ)   typ 0 : cov(...)=...   print covariance matrix.             */
/*              typ 1 : corr(...)=...  print correlation matrix.            */
/*              typ 2 : rcorr(...)=... print rank correlation matrix.       */
/*                                                                          */
/*              cov/corr(                                                   */
/*                  fmt=...,                    print format, def. 10.4     */
/*                  df=...,                     output file                 */
/*                  prn=...,                    print option                */
/*                  mpcov=...,                  create matrix               */
/*              ) [=varlist];                   optional varlist            */  
/*                                                                          */
/*              prn= 0 (def) : lower triangle, including main diagonal      */
/*              prn= 1 : full square matrix                                 */
/*              prn= 2 : column vector: lower triangle                      */
/*              prn= 3 : column vector: lower triangle incl. diag.          */
/*              prn= 4 : column vector: full matrix                         */
/*                                                                          */
/*              Return 0 if successful, otherwise -1.                       */

int pcov(TDAContext *ctx, int typ)
{
    register int i,j,k;
    int nn,vl,err,w;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    if (typ == 0)
        printf1(ctx, "Covariance");
    else if (typ == 1)
        printf1(ctx, "Correlation");
    else if (typ == 2)
        printf1(ctx, "Rank correlation");

    printf1(ctx, " matrix. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 3 + typ,4,0))     /* get parameters */
        goto COVFin;

    if (ctx->PMNV == 0) {        /* use all variables */

        if (pm_valloc(ctx, ctx->NVAR))  
            goto COVFin;
           
        w = j = 0;
        i = ctx->VIFirst;
        while (i >= 0) {
            ctx->PMVIdx[j++] = (short)(i);
            if (ctx->VTyp[i] == 5)
                w = 1;
            i = ctx->VNxt[i];
        }
        if (w)
            p_warn(ctx, -3,1);
    }   
    if (ctx->PMNV < 2) {
        p_err(ctx, -8,1);
        goto COVFin;
    }
    if (typ <= 1)
        prn_cwt(ctx);          /* case weight information */

    /* allocate memory */

    nn = ctx->PMNV * (ctx->PMNV + 1) / 2;
    if (alloc_acx(ctx, nn + 1))
        goto COVFin;

    if (alloc_acy(ctx, ctx->PMNV + 1))
        goto COVFin;

    if (ctx->PMMPCovDef == 1) {      /* use AcTmp for mpcov */
        if (alloc_actmp(ctx, ctx->PMNV * ctx->PMNV + 1))
            goto COVFin;
    }

    if (typ == 0) {
        if (cov(ctx, ctx->PMNV,ctx->PMVIdx,ctx->AcY,ctx->AcX))       /* get covariance matrix */
            goto COVFin;
    }
    else if (typ == 1) {
        if (alloc_acu(ctx, ctx->PMNV + 1))
            goto COVFin;

        if (corr(ctx, ctx->PMNV,ctx->PMVIdx,ctx->AcY,ctx->AcU,ctx->AcX))  /* get correlation matrix */
            goto COVFin;
    }
    else if (typ == 2) {
        rcorr(ctx, ctx->PMNV,ctx->PMVIdx,ctx->AcX);             /* get rank correlation matrix */
    }
   
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 7,4);

    vl = get_vnl(ctx, ctx->PMNV,ctx->PMVIdx);
    if (vl < ctx->PMFmt1)
        vl = ctx->PMFmt1;

    newline(ctx);
    prnchar(ctx, ' ',vl,0);
    k = vl;

    for (j = 0; j < ctx->PMNV; ++j) {
        prnchar(ctx, ' ',vl - (int)strlen(ctx->VName[ctx->PMVIdx[j]]),0);
        printf1(ctx, "%s ",ctx->VName[ctx->PMVIdx[j]]);
        k += vl + 1;
    }
    newline(ctx);    
    prnchar(ctx, '-',k - 1,1);
    k = 0;
    for (i = 0; i < ctx->PMNV; ++i) {

        printf1(ctx, "%s",ctx->VName[ctx->PMVIdx[i]]);
        prnchar(ctx, ' ',vl - (int)strlen(ctx->VName[ctx->PMVIdx[i]]),0);

        for (j = 0; j <= i; ++j) {
            rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcX[k]);

            if (ctx->PMMPCovDef == 1) {
                ctx->AcTmp[i * ctx->PMNV + j + 1] =               
                ctx->AcTmp[j * ctx->PMNV + i + 1] = ctx->AcX[k];       
            }
            k++;
        }
        newline(ctx);
    }
#ifdef TDA_R_PACKAGE
    /* the full symmetric matrix the triangle above prints, under the
       command's own name -- cov, corr (typ 1) or rcorr (typ 2) */
    {
        int ei, ej, ek = 0;
        double *ev = (double *)malloc((size_t)ctx->PMNV *
                                      (size_t)ctx->PMNV * sizeof(double));
        if (ev != NULL) {
            for (ei = 0; ei < ctx->PMNV; ++ei)
                for (ej = 0; ej <= ei; ++ej) {
                    ev[ej * ctx->PMNV + ei] =
                        ev[ei * ctx->PMNV + ej] = ctx->AcX[ek];
                    ek++;
                }
            tda_export_mat(ctx,
                           typ == 1 ? "corr.matrix" :
                           typ == 2 ? "rcorr.matrix" : "cov.matrix",
                           ev, ctx->PMNV, ctx->PMNV);
            free(ev);
            for (ei = 0; ei < ctx->PMNV; ++ei)
                tda_export_str_row(ctx, "covcorr.names",
                                   ctx->VName[ctx->PMVIdx[ei]]);
            tda_export_str_flush(ctx, "covcorr.names");
        }
    }
#endif
    if (ctx->PMF1Def) {      /* write to output file */

        prn_f1mat(ctx, ctx->PMF1d,ctx->PMPRNO,ctx->PMFmtS,ctx->PMNV,ctx->AcX);
        printf1(ctx, "\nMatrix written to: %s\n",ctx->PMF1dName);
    }
    if (ctx->PMMPCovDef == 1) {
        mp_putcov(ctx, ctx->PMNV,ctx->AcTmp);
        newline(ctx);
        mp_info(ctx);
    }
    err = 0;

COVFin:
    p_clean(ctx);      
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  atab()      Aggregate tabulation                                        */
/*                                                                          */
/*              atab(                                                       */
/*                  x=...,                  definition of classes           */
/*                  fmt=...,                print format, def. 10.4         */
/*                  s = 1,                  classes left open, def. s=0     */
/*                  r = 1,                  print only nonempty classes     */
/*                  mppar=...,              create matrix                   */
/*              ) = X [,Y];                 optional varlist                */  
/*                                                                          */
/*              Classify X, and optionally,                                 */
/*              tabulate means of Y according to classifying X              */
/*                                                                          */
/*              x=x1,x2,...     classes for X                               */
/*              s=1 classes open on the left, default open on right         */
/*              r=1 only print non-empty classes.                           */

int atab(TDAContext *ctx)
{
    register int i,j,k;
    int n,l,k1,k2,ix,iy,err,nrow,ncol = 0;
    double wt,tmp,xmin,xmax,ymin = 0.0,ymax = 0.0;
#ifdef TDA_R_PACKAGE
    double erow[6];
    int ecol;
#endif

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Aggregated table. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,4,1))     /* get parameters */
        goto ATFin;
   
    n = ctx->PMNTP;      /* number of arguments in PMTP[] (t=...) */
    if (n < 1) {
        printf1(ctx, "Error: need classes defined with x parameter.\n");
        goto ATFin;
    }
    if (ctx->PMNV == 1) {
        ix = ctx->PMVIdx[0];
        iy = -1;
    }
    else if (ctx->PMNV == 2) {
        ix = ctx->PMVIdx[0];
        iy = ctx->PMVIdx[1];
    }
    else {
        printf1(ctx, "Error: need one or two variables.\n");
        goto ATFin;
    }

    /* allocate memory */

    if (alloc_acn(ctx, n + 1))
        goto ATFin;

    if (alloc_acw(ctx, n + 1))
        goto ATFin;

    if (alloc_acx(ctx, n + 1))
        goto ATFin;

    if (ctx->PMNV == 2) {
        if (alloc_acy(ctx, n + 1))
            goto ATFin;
    }
    prn_cwt(ctx);          /* case weight information */

    xmin = xmax = get_data(ctx, ix,0);
    if (iy >= 0)
        ymin = ymax = get_data(ctx, iy,0);

    wt = 1.0;
    for (i = 0; i < ctx->NOC; ++i) {

        if (ctx->WIVar >= 0)  
            wt = get_data(ctx, ctx->WIVar,i) * ctx->WNorm;
          
        tmp = get_data(ctx, ix,i);

        if (xmin > tmp)
            xmin = tmp;
        if (xmax < tmp)
            xmax = tmp;

        k = n;

        for (j = 0; j < n; ++j) {
            if (tmp < ctx->PMTP[j]) {
                k = j;
                break;
            }
        }
        ctx->AcN[k] += 1;
        ctx->AcW[k] += wt;
        ctx->AcX[k] += tmp * wt;
        if (iy >= 0) {
            tmp = get_data(ctx, iy,i);
            ctx->AcY[k] += tmp * wt;
            if (ymin > tmp)
                ymin = tmp;
            if (ymax < tmp)
                ymax = tmp;
        }
    }
    for (k = 0; k <= n; ++k) {
        tmp = ctx->AcW[k];
        if (fabs(tmp) >= ctx->EPSI) {
            ctx->AcX[k] /= tmp;
            if (iy >= 0)
                ctx->AcY[k] /= tmp;
        }
    }
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);
  
    printf1(ctx, "\nX: ");            
    prn_vname(ctx, ix);
    printf1(ctx, " Minimum: %12.4f  Maximum: %12.4f\n",xmin,xmax);
    if (iy >= 0) {
        printf1(ctx, "Y: ");         
        prn_vname(ctx, iy);
        printf1(ctx, " Minimum: %12.4f  Maximum: %12.4f\n",ymin,ymax);
    }
    newline(ctx);

    l = 4;
    prnchar(ctx, ' ',ctx->PMFmt1 -  5,0); printf1(ctx, "Lower "); 
    prnchar(ctx, ' ',ctx->PMFmt1 -  5,0); printf1(ctx, "Upper "); 
                                printf1(ctx, " Frequency"); 
    prnchar(ctx, ' ',ctx->PMFmt1 -  8,0); printf1(ctx, " Weighted"); 
    prnchar(ctx, ' ',ctx->PMFmt1 -  8,0); printf1(ctx, "  Mean(X)"); 
    if (iy >= 0) {
        prnchar(ctx, ' ',ctx->PMFmt1 -  8,0); printf1(ctx, "  Mean(Y)");             
        l++;
    }
    newline(ctx);
    prnchar(ctx, '-',11 + l * (ctx->PMFmt1 + 1) - 1,1);

    k1 = 0;
    for (k = 0; k <= n; ++k) {
        if (ctx->AcN[k]) {
            k1 = k;
            break;
        }
    }
    k2 = n;
    for (k = n; k >= 0; --k) {
        if (ctx->AcN[k]) {
            k2 = k;
            break;
        }
    }
   
    if (ctx->PMMPParDef == 1) {      /* use AcTmp for mppar */
        nrow = k2 - k1 + 1;
        ncol = 5;
        if (iy >= 0)
            ncol++;
        if (alloc_actmp(ctx, nrow * ncol + 1))
            goto ATFin;
    }
    nrow = 0;
    for (k = k1; k <= k2; ++k) {

        if (ctx->PMR != 0 && ctx->AcN[k] == 0)
            continue;

        if (k > 0) { 
            rt_printf1_d(ctx, ctx->PMFmtS,ctx->PMTP[k - 1]);
            if (ctx->PMMPParDef == 1)  
                ctx->AcTmp[nrow * ncol + 1] = ctx->PMTP[k - 1];  
#ifdef TDA_R_PACKAGE
            erow[0] = ctx->PMTP[k - 1];
#endif
        }
        else {
            prnchar(ctx, ' ',ctx->PMFmt1 + 1,0);
#ifdef TDA_R_PACKAGE
            erow[0] = (double)(NAN);
#endif
        }

        if (k < n) {
            rt_printf1_d(ctx, ctx->PMFmtS,ctx->PMTP[k]);
            if (ctx->PMMPParDef == 1)  
                ctx->AcTmp[nrow * ncol + 2] = ctx->PMTP[k];  
#ifdef TDA_R_PACKAGE
            erow[1] = ctx->PMTP[k];
#endif
        }
        else {
            prnchar(ctx, ' ',ctx->PMFmt1 + 1,0);
#ifdef TDA_R_PACKAGE
            erow[1] = (double)(NAN);
#endif
        }

        printf1(ctx, "%10d ",ctx->AcN[k]);
        rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcW[k]);
        rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcX[k]);

        if (ctx->PMMPParDef == 1) {
            ctx->AcTmp[nrow * ncol + 3] = (double)ctx->AcN[k];
            ctx->AcTmp[nrow * ncol + 4] = ctx->AcW[k];
            ctx->AcTmp[nrow * ncol + 5] = ctx->AcX[k];
        }   
        if (iy >= 0) {
            rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcY[k]);
            if (ctx->PMMPParDef == 1)  
                ctx->AcTmp[nrow * ncol + 6] = ctx->AcY[k];
#ifdef TDA_R_PACKAGE
            erow[5] = ctx->AcY[k];
#endif
        }
#ifdef TDA_R_PACKAGE
        erow[2] = (double)ctx->AcN[k];
        erow[3] = ctx->AcW[k];
        erow[4] = ctx->AcX[k];
        ecol = iy >= 0 ? 6 : 5;
        tda_export_row(ctx, "atab.table", erow, ecol);
#endif
        nrow++;
        newline(ctx);
    }
    if (ctx->PMMPParDef == 1) {
        mp_putmpar(ctx, nrow,ncol,ctx->AcTmp);
        newline(ctx);
        mp_info(ctx);
    }
    err = 0;

ATFin:
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "atab.table");
#endif
    p_clean(ctx);  
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  dma         Data matrix analysis                                        */
/*                                                                          */
/*              dma(                                                        */
/*                  v=...,      optional varlist, def. all var              */
/*                  alg=...,    algorithm, def. 1                           */
/*                              1  princ. components (data matrix)          */
/*                              2  princ. components (cov./corr)            */
/*                              3  factor analysis (data matrix)            */
/*                              4  factor analysis (corr. matrix)           */
/*                              5  dual scaling                             */
/*                              6  correspondence analysis                  */
/*                              7  direct svd-based projection              */
/*                  opt=...,    preprocessing (only for algorithms 1 and 7) */
/*                              def. 1                                      */ 
/*                              for alg = 1:                                */
/*                                1  do nothing                             */
/*                                2  mean centering of variables            */
/*                                3  standardization of variables           */
/*                              for alg = 5 and 7:                          */
/*                                1  do nothing                             */
/*                                2  relative frequencies of rows           */
/*                                3  relative frequencies of table          */
/*                  df=...,     write eigenvectors to output file           */
/*                  prn=1,      additional stand. output, def 0             */
/*                  ns=...,     print only ... vectors, def all             */
/*                  fmt=...,    print format, def. 10.4                     */
/*                  pcf=...,       command file for ca plot                 */
/*              ) [= output file]; optional output file                     */  
/*                                                                          */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int dma(TDAContext *ctx)
{
    register int i,j,k;
    int err,nc,nv,nnc,nnv,nv1,r,vlen,wrec,tflag;
    double s,tmp,tmp1,tmp2,tsum;   

    wrec = 0;
    tsum = 0.0;
    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Data matrix analysis. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 3,1,0))     /* get parameters */
        goto DMAFin;

    if (ctx->PMFmtF == 0 || ctx->PMFmt1 < 8)
        pmfmt(ctx, 10,4);
   
    if (ctx->PMNV > 0)
        nv = ctx->PMNV;
    else
        nv = ctx->NVAR;

    if (ctx->PMALG < 1 || ctx->PMALG > 7)
        ctx->PMALG = 1;
    printf1(ctx, "\nAlgorithm %d: ",ctx->PMALG);
    if (ctx->PMALG == 1)
        printf1(ctx, "principal components (based on data matrix).\n");
    else if (ctx->PMALG == 2)
        printf1(ctx, "principal components (based on cov./corr. matrix).\n");
    else if (ctx->PMALG == 3)
        printf1(ctx, "factor analysis (based on data matrix).\n");
    else if (ctx->PMALG == 4)
        printf1(ctx, "factor analysis (based on correlation matrix).\n");
    if (ctx->PMALG == 5)
        printf1(ctx, "dual scaling (based on frequency table).\n");
    else if (ctx->PMALG == 6)
        printf1(ctx, "correspondence analysis (based on frequency table).\n");
    else if (ctx->PMALG == 7)
        printf1(ctx, "direct svd-based projection.\n");
    
    printf1(ctx, "Number of rows (cases): %d\n",ctx->NOC);
    printf1(ctx, "Number of columns (variables): %d\n",nv);

    if (ctx->PMALG == 2 || ctx->PMALG == 4) {
        if (ctx->NOC != nv) {
            printf1(ctx, "Error: need a symmetric input matrix.\n");
            goto DMAFin;
        }
    }
    else if (ctx->PMALG == 7) {
        if (ctx->NOC < nv) {
            printf1(ctx, "Error: number of rows is less than number of columns.\n");
            goto DMAFin;
        }
    }
    newline(ctx);

    if (alloc_aci(ctx, nv))
        goto DMAFin;

    if (ctx->PMNV > 0) {
        for (j = 0; j < ctx->PMNV; ++j)
            ctx->AcI[j] = ctx->PMVIdx[j];
    }
    else { 
        j = 0;
        i = ctx->VIFirst;
        while (i >= 0) {
            ctx->AcI[j++] = i;
            i = ctx->VNxt[i];
        }
    }
    vlen = 8;
    for (j = 0; j < nv; ++j) 
        vlen = (int)(imax(ctx, vlen,(int)(strlen(ctx->VName[ctx->AcI[j]]))));

    nc = ctx->NOC;
    if (alloc_acx(ctx, nc * nv + 1))
        goto DMAFin;

    for (i = 0; i < nc; ++i) {
        for (j = 0; j < nv; ++j) 
            ctx->AcX[i * nv + j + 1] = get_data(ctx, ctx->AcI[j],i);
    }

    /* ------------------------ alg 1 and 2 ------------------------------- */

    if (ctx->PMALG == 1 || ctx->PMALG == 2) {

        if (alloc_acy(ctx, nv + 1))            /* used for eigenvalues */
            goto DMAFin;

        if (alloc_acv(ctx, nv * nv + 1))       /* used for eigenvectors */
            goto DMAFin;

        if (ctx->PMALG == 1) {

            if (ctx->PMOPT == 2) {
                printf1(ctx, "Preprocessing: mean centering of variables.\n");
                dma_stand(ctx, nc,nv,ctx->AcX,0);
            }
            else if (ctx->PMOPT == 3) {
                printf1(ctx, "Preprocessing: standardization of variables.\n");
                dma_stand(ctx, nc,nv,ctx->AcX,1);
            }
            if (ctx->PMPRNO) {
                printf1(ctx, "Input data (after preprocessing [opt=%d])\n",ctx->PMOPT);
                dma_prn1(ctx, nc,nv,ctx->AcX,-1);
                newline(ctx);
            }
            dma_corr(ctx, nc,nv,ctx->AcX,ctx->AcV);
        }
        else if (ctx->PMALG == 2) {
            for (i = 1; i <= nv * nv; ++i)
                ctx->AcV[i] = ctx->AcX[i];
        }
        if (ctx->PMPRNO) {
            printf1(ctx, "Cross-product matrix\n");
            dma_prn1(ctx, nv,nv,ctx->AcV,-1);
            newline(ctx);
        }
        r = evecf(ctx, nv,ctx->AcY,ctx->AcV);
        if (r) {
            printf1(ctx, "No success in eigenvalue calculation.\n");
            goto DMAFin;
        }
        printf1(ctx, "Eigenvalue         per cent\n");
        tmp = 0.0;
        for (j = 1; j <= nv; ++j)
            tmp += ctx->AcY[j];

        for (j = 1; j <= nv; ++j) {
            prnchar(ctx, ' ',16 - ctx->PMFmt1,0);
            rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcY[j]);
            if (tmp > 0.0)
                printf1(ctx, "  %8.2f",100.0 * ctx->AcY[j] / tmp);
            newline(ctx);
#ifdef TDA_R_PACKAGE
            {
                double erow[2];
                erow[0] = ctx->AcY[j];
                erow[1] = tmp > 0.0 ? 100.0 * ctx->AcY[j] / tmp : (double)(NAN);
                tda_export_row(ctx, "mds.eigenvalues", erow, 2);
            }
#endif
        }
#ifdef TDA_R_PACKAGE
        tda_export_flush(ctx, "mds.eigenvalues");
#endif
        newline(ctx);
        if (ctx->PMPRNO) {
            printf1(ctx, "Eigenvectors\n");
            dma_prn1(ctx, nv,nv,ctx->AcV,ctx->PMNS);
        }
        newline(ctx);

        if (ctx->PMF1Def) {
            dma_prn3(ctx, ctx->PMF1d,nv,nv,ctx->AcV,ctx->PMNS);
            printf1(ctx, "Eigenvectors (%d records) written to: %s\n",nv,ctx->PMF1dName);
        }

        if (ctx->PMALG == 1 && ctx->PMFDef) {  /* princ. components */

            if (ctx->PMNS >= 1 && ctx->PMNS <= nv)
                r = ctx->PMNS;
            else
                r = nv;

            for (i = 0; i < nc; ++i) {
                for (j = 1; j <= r; ++j) {
                    tmp = 0.0;
                    for (k = 1; k <= r; ++k) 
                        tmp += ctx->AcX[i * nv + k] * ctx->AcV[(k - 1) * nv + j];
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
#ifdef TDA_R_PACKAGE
                    tda_export_cell(ctx, "dma.scores", tmp);
#endif
                }
                fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
                tda_export_endrow(ctx, "dma.scores");
#endif
            }
            printf1(ctx, "Principal components (%d records) written to: %s\n",nc,ctx->PMFdName);
        }
    }

    /* ------------------------ alg 3 and 4 ------------------------------- */

    else if (ctx->PMALG == 3 || ctx->PMALG == 4) {

        if (alloc_acy(ctx, nv + 1))            /* used for eigenvalues */
            goto DMAFin;

        if (alloc_acv(ctx, nv * nv + 1))       /* used for eigenvectors */
            goto DMAFin;

        if (ctx->PMALG == 3) {                 /* create correlation matrix */
            dma_stand(ctx, nc,nv,ctx->AcX,1);
            dma_corr(ctx, nc,nv,ctx->AcX,ctx->AcV);
        }
        else if (ctx->PMALG == 4) {
            for (i = 1; i <= nv * nv; ++i)
                ctx->AcV[i] = ctx->AcX[i];
        }
        if (ctx->PMPRNO) {
            printf1(ctx, "Correlation matrix\n");
            dma_prn1(ctx, nv,nv,ctx->AcV,-1);
            newline(ctx);
        }
        r = evecf(ctx, nv,ctx->AcY,ctx->AcV);
        if (r) {
            printf1(ctx, "No success in eigenvalue calculation.\n");
            goto DMAFin;
        }
        printf1(ctx, "Eigenvalue         per cent\n");
        tmp = 0.0;
        for (j = 1; j <= nv; ++j)
            tmp += ctx->AcY[j];

        for (j = 1; j <= nv; ++j) {
            prnchar(ctx, ' ',16 - ctx->PMFmt1,0);
            rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcY[j]);
            if (tmp > 0.0)
                printf1(ctx, "  %8.2f",100.0 * ctx->AcY[j] / tmp);
            newline(ctx);
#ifdef TDA_R_PACKAGE
            {
                double erow[2];
                erow[0] = ctx->AcY[j];
                erow[1] = tmp > 0.0 ? 100.0 * ctx->AcY[j] / tmp : (double)(NAN);
                tda_export_row(ctx, "mds.eigenvalues", erow, 2);
            }
#endif
        }
#ifdef TDA_R_PACKAGE
        tda_export_flush(ctx, "mds.eigenvalues");
#endif
        newline(ctx);

        if (ctx->PMALG == 3 && ctx->PMFDef) {  /* factors */

            if (ctx->PMNS >= 1 && ctx->PMNS <= nv)
                r = ctx->PMNS;
            else
                r = nv;

            for (i = 0; i < nc; ++i) {
                for (j = 1; j <= r; ++j) {
                    tmp = 0.0;
                    for (k = 1; k <= r; ++k) 
                        tmp += ctx->AcX[i * nv + k] * ctx->AcV[(k - 1) * nv + j];
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
#ifdef TDA_R_PACKAGE
                    tda_export_cell(ctx, "dma.scores", tmp);
#endif
                }
                fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
                tda_export_endrow(ctx, "dma.scores");
#endif
            }
            wrec = nc;
        }
        dma_scal(ctx, nv,ctx->AcY,ctx->AcV);
        if (ctx->PMPRNO) {
            printf1(ctx, "Rescaled eigenvectors/factor loadings\n");
            dma_prn1(ctx, nv,nv,ctx->AcV,ctx->PMNS);
        }
        newline(ctx);
        if (ctx->PMF1Def) {
            dma_prn3(ctx, ctx->PMF1d,nv,nv,ctx->AcV,ctx->PMNS);
            printf1(ctx, "Eigenvectors/factor loadings (%d records) written to: %s\n",nv,ctx->PMF1dName);
        }
        if (wrec > 0)
            printf1(ctx, "Factors (%d records) written to: %s\n",nc,ctx->PMFdName);
    }
 
    /* -##--------------------- alg 5 ------------------------------- */

    else if (ctx->PMALG == 5) {
       
        if (ctx->PMPRNO) {
            printf1(ctx, "Frequency table\n");
            dma_prn2(ctx, nc,nv,ctx->AcX,ctx->AcI);
            newline(ctx);
        }

        if (ctx->PMOPT == 2 || ctx->PMOPT == 3) {
            if (ctx->PMOPT == 2) {
                for (i = 0; i < nc; ++i) {
                    tmp = 0.0;
                    for (j = 1; j <= nv; ++j)
                        tmp += ctx->AcX[i * nv + j];
                    if (tmp != 0.0) {
                        for (j = 1; j <= nv; ++j)
                            ctx->AcX[i * nv + j] /= tmp;
                    }
                }
            }
            else if (ctx->PMOPT == 3) {
                tmp = 0.0;
                for (i = 1; i <= nc * nv; ++i)  
                    tmp += ctx->AcX[i];
                if (tmp != 0.0) {
                    for (i = 1; i <= nc * nv; ++i)
                        ctx->AcX[i] /= tmp;
                }
            }

            if (ctx->PMPRNO) {
                printf1(ctx, "Input data after preprocessing (opt=%d)\n",ctx->PMOPT);
                dma_prn1(ctx, nc,nv,ctx->AcX,-1);
                newline(ctx);
            }
        }

        tflag = 0;
        nnc = nc;
        nnv = nv;
        if (nc < nv) {
            tflag = 1;
            nnc = nv;
            nnv = nc;
            for (i = 0; i < nnc; ++i) {
                for (j = 0; j < nnv; ++j) 
                    ctx->AcX[i * nnv + j + 1] = get_data(ctx, ctx->AcI[i],j);
            }
        }
        if (alloc_acz(ctx, nnc + 1))           /* used for row sums */
            goto DMAFin;
        if (alloc_acw(ctx, nnv + 1))           /* used for column sums */
            goto DMAFin;

        for (i = 1; i <= nnc; ++i) {
            for (j = 1; j <= nnv; ++j) {
                tmp = ctx->AcX[(i - 1) * nnv + j];
                ctx->AcZ[i] += tmp;
                ctx->AcW[j] += tmp;
            }
        }
        for (i = 1; i <= nnc; ++i) {
            for (j = 1; j <= nnv; ++j) {
                tmp = ctx->AcZ[i] * ctx->AcW[j];
                if (tmp <= 0.0) {
                    printf1(ctx, "Error: table has zero row or column.\n");
                    goto DMAFin;
                }
                ctx->AcX[(i - 1) * nnv + j] /= sqrt(tmp);
            }
        }
        if (ctx->PMF1Def) {
            for (i = 0; i < nnc; ++i) {
                for (j = 1; j <= nnv; ++j)  
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[i * nnv + j]);
                fprintf(ctx->PMF1d,"\n");
            }
        }

        if (alloc_acy(ctx, nnv + 1))   
            goto DMAFin;
        if (alloc_acu(ctx, nnc * nnv + 1))  
            goto DMAFin;
        if (alloc_acv(ctx, nnv * nnv + 1))  
            goto DMAFin;

        r = svdecomp(ctx, nnc,nnv,ctx->AcX,ctx->AcY,ctx->AcU,ctx->AcV,3);
        if (r) {
            p_err(ctx, -2,1);
            goto DMAFin;
        }                        
        printf1(ctx, "Singular value         eigenvalue  per cent\n");
        tmp = 0.0;
        for (j = 1; j <= nnv; ++j)
            tmp += ctx->AcY[j] * ctx->AcY[j];

        for (j = 1; j <= nnv; ++j) {
            s = ctx->AcY[j];
            prnchar(ctx, ' ',16 - ctx->PMFmt1,0);
            rt_printf1_d(ctx, ctx->PMFmtS,s);
            s *= s;
            prnchar(ctx, ' ',16 - ctx->PMFmt1,0);
            rt_printf1_d(ctx, ctx->PMFmtS,s);
            if (tmp > 0.0)
                printf1(ctx, " %8.2f",100.0 * s / tmp);
            newline(ctx);
        }

        printf1(ctx, "\nRow ");
        for (j = 1; j <= nnv; ++j) {
            prnchar(ctx, ' ',ctx->PMFmt1 - 8,0);
            printf1(ctx, "score%3d ",j);
        }
        newline(ctx);
        for (i = 1; i <= nc; ++i) {
            printf1(ctx, "%3d ",i);
            for (j = 1; j <= nnv; ++j) {
                if (tflag == 0)
                    rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcU[(i - 1) * nnv + j] / sqrt(ctx->AcZ[i]));
                else
                    rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcV[(i - 1) * nnv + j] / sqrt(ctx->AcW[i]));
            }
            newline(ctx);
        }

        printf1(ctx, "\nCol ");
        for (j = 1; j <= nnv; ++j) {
            prnchar(ctx, ' ',ctx->PMFmt1 - 8,0);
            printf1(ctx, "score%3d ",j);
        }
        newline(ctx);
        for (i = 1; i <= nv; ++i) {
            printf1(ctx, "%3d ",i);
            for (j = 1; j <= nnv; ++j) {
                if (tflag == 0)
                    rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcV[(i - 1) * nnv + j] / sqrt(ctx->AcW[i]));
                else
                    rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcU[(i - 1) * nnv + j] / sqrt(ctx->AcZ[i]));
            }
            newline(ctx);
        }
        if (ctx->PMF1Def)
            printf1(ctx, "\nA-matrix (%d records) written to: %s\n\n",nnc,ctx->PMF1dName);
    }

    /* ------------------------ alg 6 ------------------------------- */

    else if (ctx->PMALG == 6) {

        if (ctx->PMPRNO) {
            printf1(ctx, "Frequency table\n");
            dma_prn2(ctx, nc,nv,ctx->AcX,ctx->AcI);
            newline(ctx);
        }
        tflag = 0;      /* used for transposition */
        nnc = nc;
        nnv = nv;
        if (nc < nv) {
            tflag = 1;
            nnc = nv;
            nnv = nc;
            for (i = 0; i < nnc; ++i) {
                for (j = 0; j < nnv; ++j) 
                    ctx->AcX[i * nnv + j + 1] = get_data(ctx, ctx->AcI[i],j);
            }
        }
        tsum = 0.0;
        for (i = 1; i <= nc * nv; ++i)
            tsum += ctx->AcX[i];

        if (tsum <= ctx->EPSI1 || nc < 2 || nv < 2) {
            printf1(ctx, "Error: invalid table.\n");
            goto DMAFin;
        }
        for (i = 1; i <= nc * nv; ++i)
            ctx->AcX[i] /= tsum;

        if (alloc_acz(ctx, nnc + 1))           /* used for row sums */
            goto DMAFin;
        if (alloc_acw(ctx, nnv + 1))           /* used for column sums */
            goto DMAFin;

        for (i = 1; i <= nnc; ++i) {
            for (j = 1; j <= nnv; ++j) {
                tmp = ctx->AcX[(i - 1) * nnv + j];
                ctx->AcZ[i] += tmp;
                ctx->AcW[j] += tmp;
            }
        }
        for (i = 1; i <= nnc; ++i) {
            for (j = 1; j <= nnv; ++j) {
                tmp = ctx->AcZ[i] * ctx->AcW[j];
                if (tmp <= 0.0) {
                    printf1(ctx, "Error: table has zero row or column.\n");
                    goto DMAFin;
                }
                ctx->AcX[(i - 1) * nnv + j] =
                      (ctx->AcX[(i - 1) * nnv + j] - ctx->AcZ[i] * ctx->AcW[j]) / sqrt(tmp);
            }
        }
        if (ctx->PMPRNO) {
            printf1(ctx, "Standardized residuals\n");
            dma_prn1(ctx, nc,nv,ctx->AcX,-1);
            newline(ctx);
        }

        if (alloc_acy(ctx, nnv + 1))   
            goto DMAFin;
        if (alloc_acu(ctx, nnc * nnv + 1))  
            goto DMAFin;
        if (alloc_acv(ctx, nnv * nnv + 1))  
            goto DMAFin;

        r = svdecomp(ctx, nnc,nnv,ctx->AcX,ctx->AcY,ctx->AcU,ctx->AcV,3);
        if (r) {
            p_err(ctx, -2,1);
            goto DMAFin;
        }                        
        printf1(ctx, "Singular value         eigenvalue  per cent\n");
        tmp = 0.0;
        for (j = 1; j <= nnv; ++j)
            tmp += ctx->AcY[j] * ctx->AcY[j];

        for (j = 1; j <= nnv; ++j) {
            s = ctx->AcY[j];
            prnchar(ctx, ' ',16 - ctx->PMFmt1,0);
            rt_printf1_d(ctx, ctx->PMFmtS,s);
            s *= s;
            prnchar(ctx, ' ',16 - ctx->PMFmt1,0);
            rt_printf1_d(ctx, ctx->PMFmtS,s);
            if (tmp > 0.0)
                printf1(ctx, " %8.2f",100.0 * s / tmp);
            newline(ctx);
        }
        if (ctx->PMPCFDef) {
            if (alloc_actmp(ctx, 2 * nc + 1))   
                goto DMAFin;
            if (alloc_actmp1(ctx, 2 * nv + 1))   
                goto DMAFin;
        }

        printf1(ctx, "\nRow ");
        prnchar(ctx, ' ',ctx->PMFmt1 - 4,0); printf1(ctx, "Mass ");
        prnchar(ctx, ' ',ctx->PMFmt1 - 3,0); printf1(ctx, "PC1 ");
        prnchar(ctx, ' ',ctx->PMFmt1 - 3,0); printf1(ctx, "PC2 ");
        prnchar(ctx, ' ',ctx->PMFmt1 - 3,0); printf1(ctx, "SC1 ");
        prnchar(ctx, ' ',ctx->PMFmt1 - 3,0); printf1(ctx, "SC2\n");

        for (i = 1; i <= nc; ++i) {
            printf1(ctx, "%3d ",i);
            if (tflag == 0) {
                tmp = ctx->AcZ[i];
                tmp1 = ctx->AcU[(i - 1) * nnv + 1] / sqrt(ctx->AcZ[i]);
                tmp2 = ctx->AcU[(i - 1) * nnv + 2] / sqrt(ctx->AcZ[i]);
            }
            else {
                tmp = ctx->AcW[i];
                tmp1 = ctx->AcV[(i - 1) * nnv + 1] / sqrt(ctx->AcW[i]);
                tmp2 = ctx->AcV[(i - 1) * nnv + 2] / sqrt(ctx->AcW[i]);
            }
            rt_printf1_d(ctx, ctx->PMFmtS,tmp);
            rt_printf1_d(ctx, ctx->PMFmtS,tmp1 * ctx->AcY[1]);
            rt_printf1_d(ctx, ctx->PMFmtS,tmp2 * ctx->AcY[2]);
            rt_printf1_d(ctx, ctx->PMFmtS,tmp1);
            rt_printf1_d(ctx, ctx->PMFmtS,tmp2);
            newline(ctx);

            if (ctx->PMPCFDef) {
                ctx->AcTmp[(i - 1) * 2 + 1] = tmp1 * ctx->AcY[1];
                ctx->AcTmp[(i - 1) * 2 + 2] = tmp2 * ctx->AcY[2];
            }
            if (ctx->PMF1Def) {
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,tmp);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,tmp1 * ctx->AcY[1]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,tmp2 * ctx->AcY[2]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,tmp1);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,tmp2);
                fprintf(ctx->PMF1d,"\n");
            }
        }

        printf1(ctx, "\nCol ");
        prnchar(ctx, ' ',ctx->PMFmt1 - 4,0); printf1(ctx, "Mass ");
        prnchar(ctx, ' ',ctx->PMFmt1 - 3,0); printf1(ctx, "PC1 ");
        prnchar(ctx, ' ',ctx->PMFmt1 - 3,0); printf1(ctx, "PC2 ");
        prnchar(ctx, ' ',ctx->PMFmt1 - 3,0); printf1(ctx, "SC1 ");
        prnchar(ctx, ' ',ctx->PMFmt1 - 3,0); printf1(ctx, "SC2\n");

        for (i = 1; i <= nv; ++i) {
            printf1(ctx, "%3d ",i);
            if (tflag == 0) {
                tmp = ctx->AcW[i];
                tmp1 = ctx->AcV[(i - 1) * nnv + 1] / sqrt(ctx->AcW[i]);
                tmp2 = ctx->AcV[(i - 1) * nnv + 2] / sqrt(ctx->AcW[i]);
            }
            else {
                tmp = ctx->AcZ[i];
                tmp1 = ctx->AcU[(i - 1) * nnv + 1] / sqrt(ctx->AcZ[i]);
                tmp2 = ctx->AcU[(i - 1) * nnv + 2] / sqrt(ctx->AcZ[i]);
            }
            rt_printf1_d(ctx, ctx->PMFmtS,tmp);
            rt_printf1_d(ctx, ctx->PMFmtS,tmp1 * ctx->AcY[1]);
            rt_printf1_d(ctx, ctx->PMFmtS,tmp2 * ctx->AcY[2]);
            rt_printf1_d(ctx, ctx->PMFmtS,tmp1);
            rt_printf1_d(ctx, ctx->PMFmtS,tmp2);
            newline(ctx);

            if (ctx->PMPCFDef) {
                ctx->AcTmp1[(i - 1) * 2 + 1] = tmp1 * ctx->AcY[1];
                ctx->AcTmp1[(i - 1) * 2 + 2] = tmp2 * ctx->AcY[2];
            }

            if (ctx->PMFDef) {
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp1 * ctx->AcY[1]);
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp2 * ctx->AcY[2]);
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp1);
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp2);
                fprintf(ctx->PMFd,"\n");
            }
        }
        newline(ctx);
        if (ctx->PMF1Def)
            printf1(ctx, "Row coordinates (%d records) written to: %s\n",nc,ctx->PMF1dName);
        if (ctx->PMFDef)
            printf1(ctx, "Column coordinates (%d records) written to: %s\n",nv,ctx->PMFdName);
        if (ctx->PMPCFDef)
            dma_pcf(ctx, nc,ctx->AcTmp,nv,ctx->AcTmp1);

    }

    /* --###------------------- alg 7 ------------------------------- */

    else if (ctx->PMALG == 7) {

        if (ctx->PMPRNO) {
            printf1(ctx, "Input data\n");
            dma_prn1(ctx, nc,nv,ctx->AcX,-1);
            newline(ctx);
        }
                              
        if (ctx->PMOPT == 2 || ctx->PMOPT == 3) {
            if (ctx->PMOPT == 2) {
                for (i = 0; i < nc; ++i) {
                    tmp = 0.0;
                    for (j = 1; j <= nv; ++j)
                        tmp += ctx->AcX[i * nv + j];
                    if (tmp != 0.0) {
                        for (j = 1; j <= nv; ++j)
                            ctx->AcX[i * nv + j] /= tmp;
                    }
                }
            }
            else if (ctx->PMOPT == 3) {
                tmp = 0.0;
                for (i = 1; i <= nc * nv; ++i)  
                    tmp += ctx->AcX[i];
                if (tmp != 0.0) {
                    for (i = 1; i <= nc * nv; ++i)
                        ctx->AcX[i] /= tmp;
                }
            }
        }
        if (ctx->PMPRNO) {
            printf1(ctx, "Input data after preprocessing (opt=%d)\n",ctx->PMOPT);
            dma_prn1(ctx, nc,nv,ctx->AcX,-1);
            newline(ctx);
        }

        if (alloc_acy(ctx, nv + 1))   
            goto DMAFin;
        if (alloc_acu(ctx, nc * nv + 1))  
            goto DMAFin;
        if (alloc_acv(ctx, nv * nv + 1))  
            goto DMAFin;

        r = svdecomp(ctx, nc,nv,ctx->AcX,ctx->AcY,ctx->AcU,ctx->AcV,3);
        if (r) {
            p_err(ctx, -2,1);
            goto DMAFin;
        }                        
        printf1(ctx, "Singular value         eigenvalue  per cent\n");
        tmp = 0.0;
        for (j = 1; j <= nv; ++j)
            tmp += ctx->AcY[j] * ctx->AcY[j];

        for (j = 1; j <= nv; ++j) {
            s = ctx->AcY[j];
            prnchar(ctx, ' ',16 - ctx->PMFmt1,0);
            rt_printf1_d(ctx, ctx->PMFmtS,s);
            s *= s;
            prnchar(ctx, ' ',16 - ctx->PMFmt1,0);
            rt_printf1_d(ctx, ctx->PMFmtS,s);
            if (tmp > 0.0)
                printf1(ctx, " %8.2f",100.0 * s / tmp);
            newline(ctx);
        }
        nv1 = imin(ctx, nv,2);
        printf1(ctx, "\nCoordinates\n");
        for (i = 1; i <= nc; ++i) {
            for (j = 1; j <= nv1; ++j) {
                tmp = 0.0;
                for (k = 1; k <= nv; ++k) 
                    tmp += ctx->AcX[(i - 1) * nv + k] * ctx->AcV[(k - 1) * nv + j];
                rt_printf1_d(ctx, ctx->PMFmtS,tmp);

                if (ctx->PMFDef)  
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
            }   
            newline(ctx);
            if (ctx->PMFDef)  
                fprintf(ctx->PMFd,"\n");
        }
        newline(ctx);
        if (ctx->PMFDef)
            printf1(ctx, "Coordinates (%d records) written to: %s\n",nc,ctx->PMFdName);
    }
    err = 0;

DMAFin:
    p_clean(ctx);  
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  dma_prn1(m,n,x,p)     print (m,n)-matrix to stdout                      */

void dma_prn1(TDAContext *ctx, int m,int n,double *x,int p)
{
    register int i,j;

    if (p < 1 || p > n)
        p = n;

    for (i = 0; i < m; ++i) {
        for (j = 1; j <= p; ++j)
            rt_printf1_d(ctx, ctx->PMFmtS,x[i * n + j]);
        newline(ctx);
    }
}

/* ------------------------------------------------------------------------ */
/*  dma_prn2(m,n,x,vidx)   print table to stdout                            */

void dma_prn2(TDAContext *ctx, int m,int n,double *x,int *vidx)
{
    register int i,j,l;
    double tmp,sum;

    if (alloc_actmp(ctx, n)) {
        p_err(ctx, -2,1);
        return;
    }
    printf1(ctx, " Row ");
    for (j = 0; j < n; ++j) {
        l = (int)(strlen(ctx->VName[vidx[j]]));
        prnchar(ctx, ' ',ctx->PMFmt1 - l,0);
        printf1(ctx, "%s ",ctx->VName[vidx[j]]);
    }
    prnchar(ctx, ' ',ctx->PMFmt1 - 3,0);
    printf1(ctx, "Sum\n");

    for (i = 0; i < m; ++i) {
        printf1(ctx, "%4d ",i + 1);
        sum = 0.0;
        for (j = 0; j < n; ++j) {
            l = (int)(strlen(ctx->VName[vidx[j]]));
            prnchar(ctx, ' ',l - ctx->PMFmt1,0);
            tmp = x[i * n + j + 1];
            rt_printf1_d(ctx, ctx->PMFmtS,tmp);
            sum += tmp;
            ctx->AcTmp[j] += tmp;
        }
        rt_printf1_d(ctx, ctx->PMFmtS,sum);
        newline(ctx);
    }
    printf1(ctx, " Sum ");
    sum = 0.0;
    for (j = 0; j < n; ++j) {
        l = (int)(strlen(ctx->VName[vidx[j]]));
        prnchar(ctx, ' ',l - ctx->PMFmt1,0);
        tmp = ctx->AcTmp[j];
        rt_printf1_d(ctx, ctx->PMFmtS,tmp);
        sum += tmp;
    }
    rt_printf1_d(ctx, ctx->PMFmtS,sum);
    newline(ctx);
}

/* ------------------------------------------------------------------------ */
/*  dma_prn3(fd,m,n,x)   print (m,n)-matrix to fd                           */

void dma_prn3(TDAContext *ctx, FILE *fd,int m,int n,double *x,int p)
{
    register int i,j;

    if (p < 1 || p > n)
        p = n;

    for (i = 0; i < m; ++i) {
        for (j = 1; j <= p; ++j) {
            rt_fprintf_d(ctx, fd,ctx->PMFmtS,x[i * n + j]);
#ifdef TDA_R_PACKAGE
            tda_export_cell(ctx, "dma.vectors", x[i * n + j]);
#endif
        }
        fprintf(fd,"\n");
#ifdef TDA_R_PACKAGE
        tda_export_endrow(ctx, "dma.vectors");
#endif
    }
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "dma.vectors");
#endif
}

/* ------------------------------------------------------------------------ */
/*  dma_scal(n,d,x)   rescale x with d                                      */

void dma_scal(TDAContext *ctx, int n,double *d,double *x)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j;
    double tmp;

    for (j = 1; j <= n; ++j) {
        tmp = 0.0;
        for (i = 0; i < n; ++i)
            tmp += x[i * n + j] * x[i * n + j];
     
        if (tmp > 0.0 && d[j] > 0.0) {
            tmp = sqrt(d[j] / tmp);
            for (i = 0; i < n; ++i)
                x[i * n + j] *= tmp;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  dma_stand(n,m,x,opt)   If opt = 0 mean-centering, if opt = 1 stand. of  */
/*                         columns in (n,m)-matrix x.                       */

void dma_stand(TDAContext *ctx, int n,int m,double *x,int opt)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j;
    double tmp,s;

    for (j = 1; j <= m; ++j) {
        tmp = 0.0;
        for (i = 0; i < n; ++i) 
            tmp += x[i * m + j];
        tmp /= (double)n;
        for (i = 0; i < n; ++i)
            x[i * m + j] -= tmp;
    }
    if (opt) {
        for (j = 1; j <= m; ++j) {
            s = 0.0;
            for (i = 0; i < n; ++i)
                s += x[i * m + j] * x[i * m + j];
            if (n > 1)
                s /= (double)(n - 1);
            if (s > 0.0) {
                s = sqrt(s);
                for (i = 0; i < n; ++i)
                    x[i * m + j] /= s;
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  dma_corr(n,m,x,a)  create correlation matrix of x in a.                 */

void dma_corr(TDAContext *ctx, int n,int m,double *x,double *a)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j,k;
    double tmp,s;

    s = (double)n;
    if (s >= 2.0)
        s -= 1.0;

    for (i = 1; i <= m; ++i) {
        for (j = 1; j <= m; ++j) {
            tmp = 0.0;
            for (k = 0; k < n; ++k)
                tmp += x[ k * m + i] * x[ k * m + j];
            a[(i - 1) * m + j] = tmp / s;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  dma_pcf(nr,rx,nc,cx)                                                    */

void dma_pcf(TDAContext *ctx, int nr,double *rx,int nc,double *cx)
{
    register int i;
    double xmin,xmax,ymin,ymax,dx,dy,d;

    xmin = xmax = rx[1];
    for (i = 2; i <= nr; ++i) {
        xmin = dmin(ctx, xmin,rx[(i - 1) * 2 + 1]);
        xmax = dmax(ctx, xmax,rx[(i - 1) * 2 + 1]);
    }
    for (i = 1; i <= nc; ++i) {
        xmin = dmin(ctx, xmin,cx[(i - 1) * 2 + 1]);
        xmax = dmax(ctx, xmax,cx[(i - 1) * 2 + 1]);
    }

    ymin = ymax = rx[2];
    for (i = 2; i <= nr; ++i) {
        ymin = dmin(ctx, ymin,rx[(i - 1) * 2 + 2]);
        ymax = dmax(ctx, ymax,rx[(i - 1) * 2 + 2]);
    }
    for (i = 1; i <= nc; ++i) {
        ymin = dmin(ctx, ymin,cx[(i - 1) * 2 + 2]);
        ymax = dmax(ctx, ymax,cx[(i - 1) * 2 + 2]);
    }
    dx = xmax - xmin;
    dy = ymax - ymin;
    xmin -= 0.1 * dx;
    xmax += 0.1 * dx;
    ymin -= 0.1 * dy;
    ymax += 0.1 * dy;
    d = dx / 40.0;

    fprintf(ctx->PMPCFd,"psfile = %s.ps;\n",ctx->PMPCFName);
    fprintf(ctx->PMPCFd,"psetup(\n");
    fprintf(ctx->PMPCFd,"    pxa=%8.4f,%8.4f,\n",xmin,xmax);
    fprintf(ctx->PMPCFd,"    pya=%8.4f,%8.4f,\n",ymin,ymax);
    fprintf(ctx->PMPCFd,"    pxlen=100,\n");
    fprintf(ctx->PMPCFd,"    pylen=100);\n");
    
    fprintf(ctx->PMPCFd,"plotp=%8.4f,%8.4f,%8.4f,%8.4f;\n",0.0,ymin,0.0,ymax);
    fprintf(ctx->PMPCFd,"plotp=%8.4f,%8.4f,%8.4f,%8.4f;\n",xmin,0.0,xmax,0.0);
  
    fprintf(ctx->PMPCFd,"plotp(s=4,fs=2,lt=0)=\n");
    for (i = 1; i <= nr; ++i) {
        fprintf(ctx->PMPCFd,"    %8.4f,%8.4f",rx[(i - 1) * 2 + 1],rx[(i - 1) * 2 + 2]);
        if (i < nr)
            fprintf(ctx->PMPCFd,",\n");
        else
            fprintf(ctx->PMPCFd,";\n");
    }
    fprintf(ctx->PMPCFd,"plotp(s=5,fs=2,lt=0)=\n");
    for (i = 1; i <= nc; ++i) {
        fprintf(ctx->PMPCFd,"    %8.4f,%8.4f",cx[(i - 1) * 2 + 1],cx[(i - 1) * 2 + 2]);
        if (i < nc)
            fprintf(ctx->PMPCFd,",\n");
        else
            fprintf(ctx->PMPCFd,";\n");
    }
    for (i = 1; i <= nr; ++i) {
        fprintf(ctx->PMPCFd,"pltext(fs=2,xy=%8.4f,%8.4f)=R%d;\n",
                        rx[(i - 1) * 2 + 1] + d,rx[(i - 1) * 2 + 2],i);
    }
    for (i = 1; i <= nc; ++i) {
        fprintf(ctx->PMPCFd,"pltext(fs=2,xy=%8.4f,%8.4f)=C%d;\n",
                        cx[(i - 1) * 2 + 1] + d,cx[(i - 1) * 2 + 2],i);
    }

    printf1(ctx, "Plot command file written to: %s\n",ctx->PMPCFName);
}

