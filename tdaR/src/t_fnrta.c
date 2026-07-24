/****************************************************************************/
/*  t_fnrta                                                                 */
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
#include "t_gf.h"
#include "t_ml.h"
#include "t_gdat.h"
#include "t_edat.h"
#include "t_rate.h"
#include "t_parm.h"
#include "t_eval.h"
#include "t_eval1.h"
#include "t_var.h"
#include "t_cdf.h"
#include "t_pgen.h"
#include "t_lsei.h"
#include "t_alloc.h"
#include "t_prate.h"
#include "t_lin.h"
#include "t_mat.h"
#include "t_matf.h"
#include "tda_context.h"

/*  functions in t_fnrta.c */

int fn_cox(TDAContext *ctx);
void fn_cox_s(TDAContext *ctx, int next,int jo,double  rate);
void ts_cox(TDAContext *ctx, int tran,int dim,int tpi,int nep);
int test_cox(TDAContext *ctx);
void prscox_hd(TDAContext *ctx, int sn,int org,int des,int grp);
void prscox(TDAContext *ctx, int tranfirst,int tran,double time,double evnt,double cens, double rs,double crate);



/****************************************************************************/
/*  fn_cox                                                                  */
/*      Cox model. Calculation of log-likelihood, gradient, hessian.        */
/*      In addition: Calculation of baseline rates if COXBLFlg is set.      */
/*                                                                          */
/*      If COXTCFlg set, calculation of global goodness-of-fit              */
/*      statistics in the array TCStat.                                     */
/*                                                                          */
/*      Working Areas: WrkD, WrkE, WrkH.                                    */

int fn_cox(TDAContext *ctx)
{
    register int ii,j,k,jv,jk,sn;
    int err,sn1,icase1,icase2,io,jo,tran;
    int org,des,org1,des1,org2,des2;
    int icase,jcase,grp,next,l,l1,tpi = 0,dim = 0,sflg,nep = 0,tranfirst;
    double x,xa,xj,wt,wt1,tmp,rate,srate = 0.0,time,ts,tf,ts2,tf2;
    double chrate,rs,evnt,cens,rs1,cens1;

    ctx->NOCUsed = err = 0;
    sflg = ctx->LGrad + ctx->LSec + ctx->COXTCFlg;
     
    if (ctx->COXBLFlg)   
        ctx->LFunc = ctx->LGrad = ctx->LSec = ctx->COXTCFlg = sflg = 0;
     
    wt = wt1 = 1.0;           
    next = 0;            /* pointer to first parameter of sn and transition */

    for (tran = 0; tran < ctx->NTran1; ++tran) {     /* loop for transitions */

        tranfirst = 1;

        sn  = ctx->SnTran1[tran];
        org = ctx->OrgTran1[tran];
        des = ctx->DesTran1[tran];

        next = ctx->PIdxPtr[tran];       /* start of parameters */

        /*  If COXTCFlg initialize the calculation of goodness-of-fit */
        /*  statistic. dim is number of parameters. tpi is a counter */
        /*  for the time period, nep is the number of episodes. */
     
        if (ctx->COXTCFlg) {
            dim = ctx->VTNum[tran];
            if (ctx->PMProtFDef) {
                fprintf(ctx->PMProtFd,"\nCox Model Goodness-of-fit Statistic.\n");
                fprintf(ctx->PMProtFd,"SN=%d  Org=%d  Des=%d  Dim=%d\n",sn,org,des,dim);
            }
            dclear(ctx, ctx->VTMax,ctx->TCU);
            dclear(ctx, ctx->VTMax * ctx->VTMax,ctx->TCV);
            tpi = 1;
            nep = 0;
        }
       
        if (ctx->PM1NV) grp = 0; else grp = -1;

        for (; grp < ctx->PM1NV; ++grp) {               /* loop for all groups */
   
            /* Print header if baseline rate calculation, and make      */
            /* actual covariate constellations in PEVal.                */
    
            if (ctx->COXBLFlg)  
                prscox_hd(ctx, sn,org,des,grp);
    
            chrate = rs = rs1 = evnt = cens = cens1 = 0.0;
            time = 0.0;
            jcase  = 0;
            icase1 = 1;
       
            if (ctx->TranTVar[tran] == 0) {
                srate = 0.0;
                if (sflg) {
                    dclear(ctx, ctx->NParm,ctx->WrkD);                                           
                    dclear(ctx, ctx->NParm,ctx->WrkE);                                           
                    dclear(ctx, ctx->HSiz,ctx->WrkH);                                          
                }
            }
      
            /* Main loop over all episodes according to their ending    */
            /* times. Only episodes with sn and org.                    */

            for (icase = 1; icase <= ctx->NOC; ++icase) {

                io = ctx->TFIdx[icase - 1];
           
                get_edat(ctx, io,&sn1,&org1,&des1,&ts,&tf);
                if (sn1 != sn || org1 != org)
                    continue;     

                if (grp >= 0) {
                    if (fabs(get_data(ctx, ctx->PM1VIdx[grp],io)) <= ctx->EPSI1)
                        continue;
                }
                if (ctx->WIVar >= 0)             /* get weights */
                    wt = get_data(ctx, ctx->WIVar,io) * ctx->WNorm;
   
                if (des1 != des) {
                    cens += wt;
     
                    /* additional variables for censored cases must be written here */

                    if (ctx->CGradFlg && io < ctx->MPGradRow && ctx->PM2NV > 0)
                        mp_putvar(ctx, io,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,io);
                    continue;       
                }
                nep++;  /* count uncensored episodes for tcstat */

                /* Check for time-varying variables or end of cases */
       
                if (jcase < ctx->NOC && ctx->TranTVar[tran] == 0) {

                    /* Update of the risk set if no time-varying variables */
                    /* Use sorted starting times of the episodes: add      */
                    /* all episodes with starting times less than tf.      */
                    /* rs is the number of elements in the risk set.       */

                    while (++jcase <= ctx->NOC) {

                        jo = ctx->TSIdx[jcase - 1];

                        get_edat(ctx, jo,&sn1,&org2,&des2,&ts2,&tf2);
                        if (sn1 != sn || org2 != org)
                            continue;     

                        if (ts2 >= tf) {
                            jcase--;
                            break;
                        }
                        if (grp >= 0) {
                            if (fabs(get_data(ctx, ctx->PM1VIdx[grp],jo)) <= ctx->EPSI1)
                                continue;
                        }
                        if (ctx->WIVar >= 0)             /* get weights */
                            wt1 = get_data(ctx, ctx->WIVar,jo) * ctx->WNorm;

                        rs += wt1;
   
                        j = next;                                                     
                        jv = -ctx->PIdx[j];
                        xa = ctx->TPar[j] * get_data(ctx, jv - 1,jo);
                        while ((jv = ctx->PIdx[++j]) > 0)  
                            xa += ctx->TPar[j] * get_data(ctx, jv - 1,jo);

                        rate = wt1 * rexp(ctx, xa);                                            
                        srate += rate;                                 

                        /* if we need first and second derivatives */
                                
                        if (sflg)   
                            fn_cox_s(ctx, next,jo,rate);
                    }
                }               /* end of update of the risk set if there are */
                                /* no time varying variables.                 */

                if (tf <= time)      
                    goto FN1CONT;

                /* If time-varying variables continue at FN1B */
                    
                if (ctx->TranTVar[tran])
                    goto FN1B;
                
                /* Reduction of the risk set in the case of no      */
                /* time-varying variables. Take all episodes of the */
                /* icase-loop which have tf2 < tf.                  */
                    
                for (ii = icase1; ii < icase; ++ii) {

                    jo = ctx->TFIdx[ii - 1];

                    get_edat(ctx, jo,&sn1,&org2,&des2,&ts2,&tf2);
                    if (sn1 != sn || org2 != org)
                        continue;     

                    if (tf2 >= tf)
                        break;

                    if (grp >= 0) {
                        if (fabs(get_data(ctx, ctx->PM1VIdx[grp],jo)) <= ctx->EPSI1)
                            continue;
                    }
                    if (ctx->WIVar >= 0)             /* get weights */
                        wt1 = get_data(ctx, ctx->WIVar,jo) * ctx->WNorm;

                    rs -= wt1;

                    j = next;                                                     
                    jv = -ctx->PIdx[j];
                    xa = ctx->TPar[j] * get_data(ctx, jv - 1,jo);
                    while ((jv = ctx->PIdx[++j]) > 0)  
                        xa += ctx->TPar[j] * get_data(ctx, jv - 1,jo);
                            
                    rate = -wt1 * rexp(ctx, xa);                                            
                    srate += rate;                                 

                    if (sflg)  
                        fn_cox_s(ctx, next,jo,rate);
                }
                icase1 = ii;    /* end of reduction of risk set if  */
                goto FN1C;      /* no time-varying covariates.      */

FN1B: ;         /* Build risk set at tf if time-varying variables.  */

                rs = srate = 0.0;
       
                if (sflg) {
                    dclear(ctx, ctx->NParm,ctx->WrkD);                                           
                    dclear(ctx, ctx->NParm,ctx->WrkE);                                           
                    dclear(ctx, ctx->HSiz,ctx->WrkH);                                          
                }
         
                /* Loop over all episodes according to ending times. */
                /* But there are possibly some episodes before icase */
                /* with the same ending time but different destinat. */
                /* state; they have to be included in the risk set.  */

                icase2 = icase;
                while (--icase2 > 0) {

                    jo = ctx->TFIdx[icase2 - 1];

                    get_edat(ctx, jo,&sn1,&org2,&des2,&ts2,&tf2);
                    /****************************/
                    if (sn1 != sn || org2 != org)
                        continue;     
                    /******************/
                    if (tf2 < tf)
                        break;
                }
                icase2++;
                for (jcase = icase2; jcase <= ctx->NOC; ++jcase) {

                    jo = ctx->TFIdx[jcase - 1];

                    get_edat(ctx, jo,&sn1,&org2,&des2,&ts2,&tf2);
                    if (sn1 != sn || org2 != org)
                        continue;     

                    /* starting time must be less than tf */      
                        
                    if (ts2 < tf) {  

                        if (grp >= 0) {
                            if (fabs(get_data(ctx, ctx->PM1VIdx[grp],jo)) <= ctx->EPSI1)
                                continue;
                        }
                        if (ctx->WIVar >= 0)             /* get weights */
                            wt1 = get_data(ctx, ctx->WIVar,jo) * ctx->WNorm;

                        rs += wt1;
                        j = next;                                                     
                        xa = 0.0;                                                     
                        jv = -ctx->PIdx[j];
                        while (1) {

                            /* If VTyp is 5 the variable */
                            /* is time-varying.          */
                         
                            if (ctx->VTyp[jv - 1] == 5) {

                                ctx->EDVALSn = sn1;      
                                ctx->EDVALOrg = org2;
                                ctx->EDVALDes = des2;
                                ctx->EDVALTs = ts2;
                                ctx->EDVALTf = tf2;
                                ctx->EDVALTime = tf;

                                err = v_eval1(ctx, jo,ctx->VESCnt[jv-1],ctx->VESTyp[jv-1],ctx->VESVal[jv-1],ctx->ESIdx,&x,0,0,0,0,0);
                                if (err) {
                                    ctx->NPFlgs[5] += 1;
                                    x = 0.0;     
                                }
                                put_data(ctx, x,jv-1,jo);
                            }
                            xa += ctx->TPar[j] * get_data(ctx, jv - 1,jo);
                            if ((jv = ctx->PIdx[++j]) < 0)
                                break;
                        }
                        rate = wt1 * rexp(ctx, xa);                                            
                        srate += rate;                                 
                        if (sflg)   
                            fn_cox_s(ctx, next,jo,rate);
                    }
                }               /* end of risk set building in the case of  */
                                /* time-varying variables.                  */
FN1C:
                /* If baseline rate calculation printing for time   */
                /* is done here.                                    */
     
                if (ctx->COXBLFlg) {
                    if (time == 0.0)
                        rs1 = rs;
                    prscox(ctx, tranfirst,tran,time,evnt,cens1,rs1,chrate);
                    tranfirst = 0;
                    cens1 = cens;
                    evnt  = cens = 0.0;
                    rs1   = rs;
                }
       
                /* Set time to tf, so if there are more episodes    */
                /* with this ending time it is not necessary to re- */
                /* calculate the risk set.                          */

                time = tf;

                /* If COXTCFlg calculate test statistics */
         
                if (ctx->COXTCFlg && tpi < ctx->PMNTP && time >= ctx->PMTP[tpi]) {
                    ts_cox(ctx, tran,dim,tpi,nep - 1);
                    nep = 1;
                    while (++tpi < ctx->PMNTP && time >= ctx->PMTP[tpi])
                        ts_cox(ctx, tran,dim,tpi,nep - 1);
                }
           
FN1CONT:        if (ctx->COXBLFlg) { /* calculation of baseline rate */
                    evnt += wt;
                    if (srate > 0.0)  
                        chrate += wt / srate;
                    goto FN1Skip;
                }
          
                /* Now calculation of function value and gradient */
                /* and hessian if requested.                      */

                j = next;
                jv = -ctx->PIdx[j];
                xa = ctx->TPar[j] * get_data(ctx, jv - 1,io);
                while ((jv = ctx->PIdx[++j]) > 0)  
                    xa += ctx->TPar[j] * get_data(ctx, jv - 1,io);

                if (ctx->LFunc)  
                    ctx->FTmp += wt * (xa - rlog(ctx, srate));                             

                if (ctx->LGrad || ctx->COXTCFlg) {
                    j = next;
                    jv = -ctx->PIdx[j];                                                 
                    l = 1;
                    while (1) {
                        xj = get_data(ctx, jv - 1,io);
                        if (ctx->LGrad) { 
                            ctx->Grad[j] += wt * (xj - ctx->WrkD[j] / srate);                            
                            if (ctx->CGradFlg && io < ctx->MPGradRow) {
                                ctx->MatVal[ctx->MPGradIdx][io * ctx->MPGradCol + j] = 
                                                wt * (xj - ctx->WrkD[j] / srate);                            

                                if (ctx->PM2NV > 0)
                                    mp_putvar(ctx, io,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],
                                                             ctx->PM2NV,ctx->PM2VIdx,io);
                            }
                        }
                        if (ctx->COXTCFlg)
                            ctx->TCU[l++] += wt * (xj - ctx->WrkD[j] / srate);
         
                        if ((jv = ctx->PIdx[++j]) < 0 )
                            break;
                    }
                }
                if (ctx->LSec || ctx->COXTCFlg) {
                    j = next;
                    jk = (j * (j - 1)) / 2 + 1;                                   
                    jv = -ctx->PIdx[j];                                                 
                    tmp = ctx->WrkD[j] / srate;                                      
                    l = 1;
                    while (1) {
                        if (ctx->LSec)
                            ctx->Diag[j] += wt * (tmp * ctx->WrkD[j] - ctx->WrkE[j]) / srate;       
  
                        if (ctx->COXTCFlg)
                            ctx->TCV[(l - 1) * dim + l] -=
                                wt * (tmp * ctx->WrkD[j] - ctx->WrkE[j]) / srate;       
    
                        if ((jv = ctx->PIdx[++j]) < 0) 
                            break;
                        
                        k = next;
                        jk += next - 1; 
                        tmp = ctx->WrkD[j] / srate;                                      
                        l++;
                        l1 = 1;
                        while (k < j) {
                            if (ctx->LSec)
                                ctx->Hess[jk] +=
                                     wt * (tmp * ctx->WrkD[k] - ctx->WrkH[jk]) / srate;   
   
                            if (ctx->COXTCFlg)
                                ctx->TCV[(l - 1) * dim + l1] -=
                                  wt * (tmp * ctx->WrkD[k] - ctx->WrkH[jk]) / srate;   
    
                            k++;
                            jk++;
                            if (k == j)
                                break;
                            l1++;
                        }
                    } 
                }

FN1Skip:   ;
            }     /* end of loop over all episodes acc. to ending times */
      
            /* If baseline rate calculation */
  
            if (ctx->COXBLFlg) {
                prscox(ctx, tranfirst,tran,time,evnt,cens1,rs1,chrate);
                tranfirst = 0;
                ctx->PRTIdx++;
            }
        }           /* end of current group */

        /* If COXTCFlg calculate test statistics for last time periods */
  
        if (ctx->COXTCFlg) {
            ts_cox(ctx, tran,dim,tpi,nep);
            nep = 0;
            while (++tpi <= ctx->PMNTP)
                ts_cox(ctx, tran,dim,tpi,nep);

            ctx->TCDf[tran] = (short)(ctx->TCDf[tran] - (dim));  /* adjust degrees of freedom */
        }
    }               /* end of main loop for all transitions     */
    ctx->NOCUsed = ctx->NOC;
    if (err == 0)
        err = checkov(ctx);
    return(err);
}

void fn_cox_s(TDAContext *ctx, int next,int jo,double  rate)
{
    register int j,k,jv,kv,jk;
    double xj,xk,xjrate;

    j = next;
    jk = (j * (j - 1)) / 2 + 1;                                       
    jv = -ctx->PIdx[j];                                                 
    xj = get_data(ctx, jv - 1,jo);
    xjrate = xj * rate;                                             

    while (1) {
    
        ctx->WrkD[j] += xjrate;                                     
    
        if (ctx->LSec || ctx->COXTCFlg)  
            ctx->WrkE[j] += xj * xjrate;                  
    
        if ((jv = ctx->PIdx[++j]) < 0)     
            break;
    
        xj = get_data(ctx, jv - 1,jo);
        xjrate = xj * rate;                                             
    
        if (ctx->LSec || ctx->COXTCFlg) {
            k = next;
            jk += next - 1; 
            kv = -ctx->PIdx[k];
            while (k < j) {
                xk = get_data(ctx, kv - 1,jo);
                ctx->WrkH[jk++] += xk * xjrate;                            
                kv = ctx->PIdx[++k];
            }
        }
    }
}

/****************************************************************************/
/*  ts_cox     Goodness-of-fit statistic for Cox models.                    */

void ts_cox(TDAContext *ctx, int tran,int dim,int tpi,int nep)
{
    register int i,j,r;
    double tmp,tmp1;

    if (ctx->PMProtFDef) {
        fprintf(ctx->PMProtFd,"\nTime period: %d\nUncensored episodes: %d\n",tpi,nep);
        prvec(ctx, "TCU vector",dim,ctx->TCU);
        prmat(ctx, "TCV matrix",dim,dim,ctx->TCV);
    }
    if (!nep) {         /* skip periods without episodes */ 
        ctx->TCDf[tran] = 0;
        return;
    }  
    ctx->TCDf[tran] = (short)(ctx->TCDf[tran] + (dim));      /* degrees of freedom */

    r = syminv(ctx, dim,ctx->TCV);    /* matrix inversion */
    if (ctx->PMProtFDef)  
        fprintf(ctx->PMProtFd,"Rank of V: %d\n",r);

    if (r == dim) {
        tmp = 0.0;
        for (i = 1; i <= dim; ++i) {
            tmp1 = 0.0;
            for (j = 1; j <= dim; ++j)
                tmp1 += ctx->TCU[j] * ctx->TCV[(j - 1) * dim + i];
            tmp += tmp1 * ctx->TCU[i];
        }
        ctx->TCStat[tran] += tmp;
        if (ctx->PMProtFDef)  
            fprintf(ctx->PMProtFd,"TStat added: %8.4lf\n",tmp);
    }
    else {
        ctx->TCStat[tran] = -1.0;
        if (ctx->PMProtFDef)  
            fprintf(ctx->PMProtFd,"TStat not calculated\n");
    }
    dclear(ctx, ctx->VTMax,ctx->TCU);
    dclear(ctx, ctx->VTMax * ctx->VTMax,ctx->TCV);
}

/****************************************************************************/
/*  test_cox    Goodness-of-fit statistic for Cox model.                    */
/*              Return 0 if OK, -1 if error.                                */

int test_cox(TDAContext *ctx)
{
    int err,n,tran,sn,org,des;
    double tmp;
   
    printf1(ctx, "\nGlobal Goodness-of-fit.\n\n");

    err = -1;
    if (!(ctx->TCDf = (short *)calloc((size_t)(ctx->MaxSnn) * (size_t)(ctx->NTran1) + 1,sizeof(short))))  
        goto TSTCFin;
    ctx->TCDfA = ctx->MaxSnn * ctx->NTran1 + 1;
    memrq(ctx, ctx->TCDfA,sizeof(short));

    if (!(ctx->TCStat = (double *)calloc((size_t)(ctx->MaxSnn) * (size_t)(ctx->NTran1) + 1,sizeof(double))))  
        goto TSTCFin;
    ctx->TCStatA = ctx->MaxSnn * ctx->NTran1 + 1;
    memrq(ctx, ctx->TCStatA,sizeof(double));

    if (!(ctx->TCU = (double *)calloc((size_t)(ctx->VTMax + 1),sizeof(double))))  
        goto TSTCFin;
    ctx->TCUA = ctx->VTMax + 1;
    memrq(ctx, ctx->TCUA,sizeof(double));

    if (!(ctx->TCV = (double *)calloc((size_t)(ctx->VTMax) * (size_t)(ctx->VTMax) + 1,sizeof(double))))  
        goto TSTCFin;
    ctx->TCVA = ctx->VTMax * ctx->VTMax + 1;
    memrq(ctx, ctx->TCVA,sizeof(double));

    ctx->COXTCFlg = 1;
    ctx->LFunc = ctx->LGrad = ctx->LSec = 0; 

    fn_cox(ctx);   /*  calculations are done directly by fn_cox    */

    printf1(ctx, "SN  Org  Des        TStat    DF   Signif\n");        
    prnchar(ctx, '-',40,1);  

    for (tran = 0; tran < ctx->NTran1; ++tran) {

        sn  = ctx->SnTran1[tran];
        org = ctx->OrgTran1[tran];
        des = ctx->DesTran1[tran];

        printf1(ctx, "%2d %4d %4d ",sn,org,des);          

        tmp = ctx->TCStat[tran];
        if (tmp > 0.0)
            printf1(ctx, "%12.4lf ",tmp);
        else                          
            printf1(ctx, "         --- ");
        n = ctx->TCDf[tran];
        printf1(ctx, "%5d ",n);
        if (tmp > 0.0 && n > 0)
            printf1(ctx, "%8.4lf\n",cdchif(ctx, tmp,n));
        else
            printf1(ctx, "     ---\n");
#ifdef TDA_R_PACKAGE
        {
            double erow[6];
            erow[0] = (double)sn;
            erow[1] = (double)org;
            erow[2] = (double)des;
            erow[3] = tmp > 0.0 ? tmp : (double)(NAN);
            erow[4] = (double)n;
            erow[5] = (tmp > 0.0 && n > 0) ? cdchif(ctx, tmp,n) : (double)(NAN);
            tda_export_row(ctx, "gof.table", erow, 6);
        }
#endif
    }
    newline(ctx);    
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "gof.table");
#endif
    err = 0;

TSTCFin:
    if (err)
        p_err(ctx, -2,1);

    if (ctx->TCDfA) {
        free((char *)ctx->TCDf);
        memrq(ctx, -ctx->TCDfA,sizeof(short));
        ctx->TCDfA = 0;
    }
    if (ctx->TCStatA) {
        free((char *)ctx->TCStat);
        memrq(ctx, -ctx->TCStatA,sizeof(double));
        ctx->TCStatA = 0;
    }
    if (ctx->TCUA) {
        free((char *)ctx->TCU);
        memrq(ctx, -ctx->TCUA,sizeof(double));
        ctx->TCUA = 0;
    }
    if (ctx->TCVA) {
        free((char *)ctx->TCV);
        memrq(ctx, -ctx->TCVA,sizeof(double));
        ctx->TCVA = 0;
    }
    ctx->COXTCFlg = 0;
    return(err);
}

/****************************************************************************/
/*  prscox_hd    Print header of prs table for Cox models.                  */
 
void prscox_hd(TDAContext *ctx, int sn,int org,int des,int grp) 
{
    fprintf(ctx->PRTFd,"\n# Transition: SN %d  Org %d  Des %d\n",sn,org,des);

    if (grp >= 0)
        fprintf(ctx->PRTFd,"# Group: %s\n",ctx->VName[ctx->PM1VIdx[grp]]);

    fprintf(ctx->PRTFd,"\n# ID ");
    fprnchar(ctx, ctx->PRTFd,' ',ctx->PMMFmt1 - 4,0); fprintf(ctx->PRTFd,"Time ");
    fprnchar(ctx, ctx->PRTFd,' ',ctx->PMMFmt1 - 6,0); fprintf(ctx->PRTFd,"Events ");
    fprnchar(ctx, ctx->PRTFd,' ',ctx->PMMFmt1 - 8,0); fprintf(ctx->PRTFd,"Censored ");
    fprnchar(ctx, ctx->PRTFd,' ',ctx->PMMFmt1 - 8,0); fprintf(ctx->PRTFd,"Risk Set ");
    fprnchar(ctx, ctx->PRTFd,' ',ctx->PMMFmt1 - 7,0); fprintf(ctx->PRTFd,"Surv.F. ");
    fprnchar(ctx, ctx->PRTFd,' ',ctx->PMMFmt1 - 8,0); fprintf(ctx->PRTFd,"Cum.Rate");
    fprintf(ctx->PRTFd,"   Baseline-Rate\n# ");
    fprnchar(ctx, ctx->PRTFd,'-',2 + 6 * (ctx->PMMFmt1 + 1) + 16,0);
    fprintf(ctx->PRTFd,"\n");
}

/****************************************************************************/
/*  prscox             print one line of prate table for Cox models.        */

void prscox(TDAContext *ctx, int tranfirst,int tran, double time, double evnt, double cens, double rs, double crate)
{
    double tmp,surv;
     
    tmp = crate * ctx->AcW[tran];
    surv = rexp(ctx, -tmp);

    fprintf(ctx->PRTFd,"%4d ",ctx->PRTIdx);
    rt_fprintf_d(ctx, ctx->PRTFd,ctx->PMMFmtS,time);
    rt_fprintf_d(ctx, ctx->PRTFd,ctx->PMMFmtS,evnt);
    rt_fprintf_d(ctx, ctx->PRTFd,ctx->PMMFmtS,cens);
    rt_fprintf_d(ctx, ctx->PRTFd,ctx->PMMFmtS,rs);
    rt_fprintf_d(ctx, ctx->PRTFd,ctx->PMMFmtS,surv);
    rt_fprintf_d(ctx, ctx->PRTFd,ctx->PMMFmtS,tmp);

    tmp = 0.0;
    if (tranfirst)     
        ctx->s_prscox_prev_time = ctx->s_prscox_prev_crate = -1.0;
    else if (time > ctx->s_prscox_prev_time)
        tmp = (crate - ctx->s_prscox_prev_crate) / (time - ctx->s_prscox_prev_time);
    fprintf(ctx->PRTFd,"%15.10f ",tmp);
    fprintf(ctx->PRTFd,"\n");
   
    ctx->s_prscox_prev_time = time;
    ctx->s_prscox_prev_crate = crate;
}




