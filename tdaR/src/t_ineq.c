/****************************************************************************/
/*  t_ineq                                                                  */
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
#include "t_gdat.h"
#include "t_var.h"
#include "t_ds.h"
#include "t_alloc.h"
#include "t_matf.h"
#include "tda_context.h"
#include "tda_compat.h"

/*  functions in t_ineq.c */

int ineq(TDAContext *ctx);
int segr(TDAContext *ctx);
double xgini(TDAContext *ctx, int n,double *x,double *wt);
int ds_comp(const void *, const void *, void *);

/*--------------------------------------------------------------------------*/
/*  xgini(n,x,wt)   Calculate Gini coefficient for x[i], i = 0,...,n-1.     */
/*                  Return: Gini.                                           */

double xgini(TDAContext *ctx, int n,double *x,double *wt)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j;
    double g,gw,tmp,ws,mean;

    if (n < 2)
        return(0.0);

    gw = g = 0.0;
    ws = wt[0];
    mean = x[0] * wt[0];

    for (i = 1; i < n; ++i) {
        for (j = 0; j < i; ++j) {
            tmp = wt[i] * wt[j];
            g += fabs(x[i] - x[j]) * tmp;
            gw += tmp; 
        }
        mean += x[i] * wt[i];
        ws += wt[i];
    }
    if (ws > 0.0 && gw > 0.0) {
        tmp = 2.0 * gw * mean / ws;
        if (tmp != 0.0)
            g /= tmp;
        else
            g = 0.0;
    }
    else
        g = 0.0;
    return(g);
}

/* ------------------------------------------------------------------------ */
/*  ineq(idx)   Inequality measures.                                        */
/*                                                                          */
/*              ineq(                                                       */
/*                  fmt=...,                    print format, def. 10.4     */
/*                  df=...,                     output file                 */
/*                  mppar=...,                  create matrix               */
/*              ) = varlist;                    required varlist            */  
/*                                                                          */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int ineq(TDAContext *ctx)
{
    register int i,j,jj;
    int err,l,n;           
    double tmp,a = 0.0,b = 0.0;
#ifdef TDA_R_PACKAGE
    double erow[7];
#endif

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Inequality measures. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,4,1))     /* get parameters */
        goto INEQFin;

    prn_cwt(ctx);      /* case weight information */
       
    if (alloc_acx(ctx, ctx->NOC))
        goto INEQFin;

    if (alloc_acw(ctx, ctx->NOC))
        goto INEQFin;

    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);
   
    if (ctx->PMMPParDef == 1) {      /* use AcTmp for mppar */
        if (alloc_actmp(ctx, ctx->PMNV * 7 + 1))
            goto INEQFin;
    }

    newline(ctx);
    prn_hvar(ctx);
    prn_hlabel(ctx);
    printf1(ctx, "  Cases "); 
    prnchar(ctx, ' ',ctx->PMFmt1 -  8,0); printf1(ctx, " Minimum"); 
    prnchar(ctx, ' ',ctx->PMFmt1 -  8,0); printf1(ctx, "  Maximum");
    prnchar(ctx, ' ',ctx->PMFmt1 -  8,0); printf1(ctx, "     Mean");
    prnchar(ctx, ' ',ctx->PMFmt1 -  6,0); printf1(ctx, " StdDev"); 
    prnchar(ctx, ' ',ctx->PMFmt1 -  6,0); printf1(ctx, " VCoeff"); 
    prnchar(ctx, ' ',ctx->PMFmt1 -  4,0); printf1(ctx, " Gini"); 
    printf1(ctx, "\n");
       
    l = ctx->VNameLen + 9; 
    if (ctx->VLabelLen > 0)
        l += ctx->VLabelLen + 1;
    prnchar(ctx, '-',l + 6 * (ctx->PMFmt1 + 1) - 1,1);
              
    for (j = 0; j < ctx->PMNV; ++j) {
        jj = ctx->PMVIdx[j];
        prn_vname(ctx, jj);
        prn_vlabel(ctx, jj);

        l = n = 0;
        for (i = 0; i < ctx->NOC; ++i) {
            tmp = get_data(ctx, jj,i);
            if (tmp >= 0.0) {
                if (l++ == 0)  
                    a = b = tmp;
                else {
                    if (a > tmp)
                        a = tmp;
                    if (b < tmp)
                        b = tmp;
                }
                ctx->AcX[n] = tmp;
                if (ctx->WIVar >= 0)
                    ctx->AcW[n] = get_data(ctx, ctx->WIVar,i);
                else
                    ctx->AcW[n] = 1.0;
                n++;
            }
        }
        printf1(ctx, "%7d ",n);
        rt_printf1_d(ctx, ctx->PMFmtS,a);
        rt_printf1_d(ctx, ctx->PMFmtS,b);
#ifdef TDA_R_PACKAGE
        erow[0] = (double)n;
        erow[1] = a;
        erow[2] = b;
#endif

        if (ctx->PMMPParDef == 1) {
            ctx->AcTmp[j * 7 + 1] = (double)n;    
            ctx->AcTmp[j * 7 + 2] = a;    
            ctx->AcTmp[j * 7 + 3] = b;    
        }

        if (ctx->PMF1Def) {
            fprintf(ctx->PMF1d,"%5d %7d ",j + 1,n);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,a);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,b);
        }
        b = xstd(ctx, n,ctx->AcX,ctx->AcW,&tmp);      /* mean and std. deviation */
        if (b > 0.0)                   /* var coefficient */
            a = tmp / b;
        else
            a = 0.0;

        rt_printf1_d(ctx, ctx->PMFmtS,b);
        rt_printf1_d(ctx, ctx->PMFmtS,tmp);
        rt_printf1_d(ctx, ctx->PMFmtS,a);
#ifdef TDA_R_PACKAGE
        erow[3] = b;
        erow[4] = tmp;
        erow[5] = a;
#endif

        if (ctx->PMMPParDef == 1) {
            ctx->AcTmp[j * 7 + 4] = b;    
            ctx->AcTmp[j * 7 + 5] = tmp;    
            ctx->AcTmp[j * 7 + 6] = a;    
        }

        if (ctx->PMF1Def) {
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,b);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,tmp);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,a);
        }
        tmp = xgini(ctx, n,ctx->AcX,ctx->AcW);        /* gini coefficient */
        rt_printf1_d(ctx, ctx->PMFmtS,tmp);
#ifdef TDA_R_PACKAGE
        erow[6] = tmp;
        tda_export_row(ctx, "ineq.table", erow, 7);
#endif

        if (ctx->PMMPParDef == 1)  
            ctx->AcTmp[j * 7 + 7] = tmp;    

        if (ctx->PMF1Def)  
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,tmp);

        printf1(ctx, "\n");
        if (ctx->PMF1Def)  
            fprintf(ctx->PMF1d,"\n");
    }
    if (ctx->PMF1Def)  
        printf1(ctx, "\nData written to: %s\n",ctx->PMF1dName);

    if (ctx->PMMPParDef == 1) {
        mp_putmpar(ctx, ctx->PMNV,7,ctx->AcTmp);
        newline(ctx);
        mp_info(ctx);
    }
    err = 0;

INEQFin:
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "ineq.table");
#endif
    p_clean(ctx);      
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  segr()      Segregation measures.                                       */
/*  ##          segr (g=, v= [,df=,fmt=] ) [=fname]                         */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int segr(TDAContext *ctx)
{
    register int i,j,k,jj;
    int err,l,n,nc;
    double tmp,tmp1,d,g,gs0,gs1;
#ifdef TDA_R_PACKAGE
    double srow[7];
#endif

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Segregation measures. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,4,1))     /* get parameters */
        goto SEGRFin;

    if (ctx->PMGIdx < 0) {
        printf1(ctx, "Error: no grouping variable.\n");
        goto SEGRFin;
    }
    printf1(ctx, "Grouping variable: %s\n",ctx->VName[ctx->PMGIdx]);

    prn_cwt(ctx);      /* case weight information */
    
    if (alloc_ack(ctx, ctx->NOC)) 
        goto SEGRFin;

    if (alloc_acn(ctx, ctx->NOC)) 
        goto SEGRFin;

    if (alloc_acm(ctx, ctx->NOC)) 
        goto SEGRFin;

    if (alloc_acc(ctx, ctx->NOC)) 
        goto SEGRFin;

    if (alloc_acw(ctx, ctx->NOC)) 
        goto SEGRFin;

    if (alloc_acx(ctx, ctx->NOC)) 
        goto SEGRFin;

    if (alloc_acy(ctx, ctx->NOC)) 
        goto SEGRFin;

    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    if (ctx->PMTFmtF == 0)
        pmtfmt(ctx, 10,4);

    newline(ctx);
    prn_hvar(ctx);
    prn_hlabel(ctx);
    printf1(ctx, " Classes      Cases    Group_0    Group_1");
    prnchar(ctx, ' ',ctx->PMFmt1 -  6,0); printf1(ctx, "D-Index"); 
    prnchar(ctx, ' ',ctx->PMFmt1 -  6,0); printf1(ctx, "V-Ratio"); 
    prnchar(ctx, ' ',ctx->PMFmt1 -  3,0); printf1(ctx, "Gini"); 
    printf1(ctx, "\n");
       
    l = ctx->VNameLen + 43; 
    if (ctx->VLabelLen > 0)
        l += ctx->VLabelLen + 1;
    prnchar(ctx, '-',l + 3 * (ctx->PMFmt1 + 1) - 1,1);
              
    for (j = 0; j < ctx->PMNV; ++j) {
        jj = ctx->PMVIdx[j];
        prn_vname(ctx, jj);
        prn_vlabel(ctx, jj);

        n = 0;
        for (i = 0; i < ctx->NOC; ++i) {

            g = (double)((float)get_data(ctx, ctx->PMGIdx,i));
            l = (int)get_data(ctx, jj,i);

            if (g >= 0.0 && l >= 0) {
                ctx->AcK[n] = l;                     
                if (g > 0)
                    ctx->AcC[n] = 1;
                else
                    ctx->AcC[n] = 0;

                if (ctx->WIVar >= 0)
                    ctx->AcW[n] = get_data(ctx, ctx->WIVar,i);
                else
                    ctx->AcW[n] = 1.0;
                n++;
            }   
        }
        if (n < 1) {
            printf1(ctx, "number of cases is zero.\n");
            continue;
        }

        /* get ptr to sort vidx values */

        for (i = 0; i < n; ++i)  
            ctx->AcN[i] = i;
   
        tda_qsort_r((char *)ctx->AcN,(size_t)(n),sizeof(int), ds_comp, ctx);

        for (i = 0; i < n; ++i)
            ctx->AcX[i] = ctx->AcY[i] = 0.0;

        k = ctx->AcN[0];
        l = ctx->AcK[k];
        ctx->AcM[k] = l;
        if (ctx->AcC[k])
            ctx->AcY[0] = ctx->AcW[k];
        else
            ctx->AcX[0] = ctx->AcW[k];
        nc = 0;

        for (i = 1; i < n; ++i) {
            k = ctx->AcN[i];
            if (ctx->AcK[k] != l) {
                l = ctx->AcK[k];
                nc++;
                ctx->AcM[nc] = l;
                if (ctx->AcC[k])
                    ctx->AcY[nc] = ctx->AcW[k];
                else
                    ctx->AcX[nc] = ctx->AcW[k];
            }
            else {
                if (ctx->AcC[k])
                    ctx->AcY[nc] += ctx->AcW[k];
                else
                    ctx->AcX[nc] += ctx->AcW[k];
            }
        }
        nc++;
        gs0 = gs1 = 0.0;
        for (i = 0; i < nc; ++i) {
            gs0 += ctx->AcX[i];
            gs1 += ctx->AcY[i];
        }
        g = gs0 + gs1;

        /* ## print to second output file if requested */

        if (ctx->PMTabFDef) {
            fprintf(ctx->PMTabFd,"Variable: %s\n",ctx->VName[jj]);
            fprintf(ctx->PMTabFd,"Class      Cases    Group_0    Group_1    D-Index\n");
            d = 0.0;                                                          
            for (i = 0; i < nc; ++i) {
                fprintf(ctx->PMTabFd,"%5d %10.2f %10.2f %10.2f ",ctx->AcM[i],
                                        ctx->AcX[i] + ctx->AcY[i],ctx->AcX[i],ctx->AcY[i]);

                if (gs0 > 0.0 && gs1 > 0.0) {
                    tmp = 0.5 * fabs((ctx->AcX[i] / gs0) - (ctx->AcY[i] / gs1));
                    fprintf(ctx->PMTabFd,"%10.4f",tmp);
                    d += tmp;
                }
                fprintf(ctx->PMTabFd,"\n");
            }
            fprintf(ctx->PMTabFd,"Total %10.2f %10.2f %10.2f ",g,gs0,gs1);
            if (gs0 > 0.0 && gs1 > 0.0)  
                fprintf(ctx->PMTabFd,"%10.4f",d);
            fprintf(ctx->PMTabFd,"\n");
        }

        printf1(ctx, " %7d %10.2f %10.2f %10.2f ",nc,g,gs0,gs1);

        if (ctx->PMF1Def)  
            fprintf(ctx->PMF1d,"%4d %7d %10.2f %10.2f %10.2f ",j + 1,nc,g,gs0,gs1);
#ifdef TDA_R_PACKAGE
        srow[0] = (double)nc;
        srow[1] = g;
        srow[2] = gs0;
        srow[3] = gs1;
        srow[6] = (double)NAN;
    srow[5] = srow[6];
    srow[4] = srow[5];
#endif

        if (gs0 <= 0.0 || gs1 <= 0) {
            newline(ctx);         
            if (ctx->PMF1Def)
                fprintf(ctx->PMF1d,"\n");
#ifdef TDA_R_PACKAGE
            tda_export_row(ctx, "segr.table", srow, 7);
#endif
            continue;
        }

        /* calculate dissimilarity index */

        d = 0.0;
        for (i = 0; i < nc; ++i)  
            d += fabs((ctx->AcX[i] / gs0) - (ctx->AcY[i] / gs1));
        d *= 0.5;

        rt_printf1_d(ctx, ctx->PMFmtS,d);
        if (ctx->PMF1Def)
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,d);
#ifdef TDA_R_PACKAGE
        srow[4] = d;
#endif

        /* calculate variance ratio index */

        d = 0.0;
        for (i = 0; i < nc; ++i) {
            tmp1 = ctx->AcX[i] + ctx->AcY[i];
            if (tmp1 > 0.0) {
                tmp = gs1 / g;
                if (ctx->AcY[i] > 0.0)
                    tmp -= ctx->AcY[i] / tmp1;
                d += tmp * tmp * tmp1;
            }
        }
        d /= (gs0 * gs1 / g);
        rt_printf1_d(ctx, ctx->PMFmtS,d);
        if (ctx->PMF1Def)
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,d);
#ifdef TDA_R_PACKAGE
        srow[5] = d;
#endif

        /* calculate gini index */

        d = 0.0;
        for (i = 0; i < nc; ++i) {
            tmp = ctx->AcX[i] + ctx->AcY[i];           
            if (tmp > 0.0) {
                for (k = 0; k < nc; ++k) {
                    if (k != i) {
                        tmp1 = ctx->AcX[k] + ctx->AcY[k];           
                        if (tmp1 > 0.0)  
                            d += tmp * tmp1 *
                                       fabs((ctx->AcY[i] / tmp) - (ctx->AcY[k] / tmp1));
                    }
                }
            }
        }
        d = 0.5 * d / (gs0 * gs1);
        rt_printf1_d(ctx, ctx->PMFmtS,d);
        if (ctx->PMF1Def)
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,d);
#ifdef TDA_R_PACKAGE
        srow[6] = d;
        tda_export_row(ctx, "segr.table", srow, 7);
#endif

        newline(ctx);
        if (ctx->PMF1Def)
            fprintf(ctx->PMF1d,"\n");
    }
    if (ctx->PMF1Def || ctx->PMTabFDef)
        newline(ctx);
    if (ctx->PMF1Def)  
        printf1(ctx, "Data written to: %s\n",ctx->PMF1dName);
    if (ctx->PMTabFDef)  
        printf1(ctx, "Additional data written to: %s\n",ctx->PMTabFName);

    err = 0;

SEGRFin:
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "segr.table");
#endif
    p_clean(ctx); 
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ds_comp()   compare function for AcK                                    */

int ds_comp(const void *arg1, const void *arg2, void *_ctx)
{
    TDAContext *ctx = (TDAContext *)_ctx;     
    int n; 

    n = ctx->AcK[*(int *)arg1] - ctx->AcK[*(int *)arg2];   
    if (n > 0)
        return(1);
    else if (n < 0)
        return(-1);
    return(0);
}

