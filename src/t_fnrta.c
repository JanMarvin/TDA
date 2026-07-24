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

/*  functions in t_fnrta.c */

int fn_cox(void);
void fn_cox_s(int next,int jo,double  rate);
void ts_cox(int tran,int dim,int tpi,int nep);
int test_cox(void);
void prscox_hd(int sn,int org,int des,int grp);
void prscox(int tranfirst,int tran,double time,double evnt,double cens,
    double rs,double crate);

int COXBLFlg = 0;           /* set during base line rate calculation        */  
int COXTCFlg = 0;           /* set during goodness of fit calculation.      */  

short *TCDf;
int TCDfA = 0;                     
double *TCStat;
int TCStatA = 0;                     
double *TCU;
int TCUA = 0;               
double *TCV;
int TCVA = 0;               

/****************************************************************************/
/*  fn_cox                                                                  */
/*      Cox model. Calculation of log-likelihood, gradient, hessian.        */
/*      In addition: Calculation of baseline rates if COXBLFlg is set.      */
/*                                                                          */
/*      If COXTCFlg set, calculation of global goodness-of-fit              */
/*      statistics in the array TCStat.                                     */
/*                                                                          */
/*      Working Areas: WrkD, WrkE, WrkH.                                    */

int fn_cox(void)
{
    register int ii,j,k,jv,jk,sn;
    int err,sn1,icase1,icase2,io,jo,tran;
    int org,des,org1,des1,org2,des2;
    int icase,jcase,grp,next,l,l1,tpi,dim,sflg,nep,tranfirst;
    double x,xa,xj,wt,wt1,tmp,rate,srate,time,ts,tf,ts2,tf2;
    double chrate,rs,evnt,cens,rs1,cens1;

    NOCUsed = err = 0;
    sflg = LGrad + LSec + COXTCFlg;
     
    if (COXBLFlg)   
        LFunc = LGrad = LSec = COXTCFlg = sflg = 0;
     
    wt = wt1 = 1.0;           
    next = 0;            /* pointer to first parameter of sn and transition */

    for (tran = 0; tran < NTran1; ++tran) {     /* loop for transitions */

        tranfirst = 1;

        sn  = SnTran1[tran];
        org = OrgTran1[tran];
        des = DesTran1[tran];

        next = PIdxPtr[tran];       /* start of parameters */

        /*  If COXTCFlg initialize the calculation of goodness-of-fit */
        /*  statistic. dim is number of parameters. tpi is a counter */
        /*  for the time period, nep is the number of episodes. */
     
        if (COXTCFlg) {
            dim = VTNum[tran];
            if (PMProtFDef) {
                fprintf(PMProtFd,"\nCox Model Goodness-of-fit Statistic.\n");
                fprintf(PMProtFd,"SN=%d  Org=%d  Des=%d  Dim=%d\n",sn,org,des,dim);
            }
            dclear(VTMax,TCU);
            dclear(VTMax * VTMax,TCV);
            tpi = 1;
            nep = 0;
        }
       
        if (PM1NV) grp = 0; else grp = -1;

        for (; grp < PM1NV; ++grp) {               /* loop for all groups */
   
            /* Print header if baseline rate calculation, and make      */
            /* actual covariate constellations in PEVal.                */
    
            if (COXBLFlg)  
                prscox_hd(sn,org,des,grp);
    
            chrate = rs = rs1 = evnt = cens = cens1 = 0.0;
            time = 0.0;
            jcase  = 0;
            icase1 = 1;
       
            if (TranTVar[tran] == 0) {
                srate = 0.0;
                if (sflg) {
                    dclear(NParm,WrkD);                                           
                    dclear(NParm,WrkE);                                           
                    dclear(HSiz,WrkH);                                          
                }
            }
      
            /* Main loop over all episodes according to their ending    */
            /* times. Only episodes with sn and org.                    */

            for (icase = 1; icase <= NOC; ++icase) {

                io = TFIdx[icase - 1];
           
                get_edat(io,&sn1,&org1,&des1,&ts,&tf);
                if (sn1 != sn || org1 != org)
                    continue;     

                if (grp >= 0) {
                    if (fabs(get_data(PM1VIdx[grp],io)) <= EPSI1)
                        continue;
                }
                if (WIVar >= 0)             /* get weights */
                    wt = get_data(WIVar,io) * WNorm;
   
                if (des1 != des) {
                    cens += wt;
     
                    /* additional variables for censored cases must be written here */

                    if (CGradFlg && io < MPGradRow && PM2NV > 0)
                        mp_putvar(io,MPGradCol,MatVal[MPGradIdx],PM2NV,PM2VIdx,io);
                    continue;       
                }
                nep++;  /* count uncensored episodes for tcstat */

                /* Check for time-varying variables or end of cases */
       
                if (jcase < NOC && TranTVar[tran] == 0) {

                    /* Update of the risk set if no time-varying variables */
                    /* Use sorted starting times of the episodes: add      */
                    /* all episodes with starting times less than tf.      */
                    /* rs is the number of elements in the risk set.       */

                    while (++jcase <= NOC) {

                        jo = TSIdx[jcase - 1];

                        get_edat(jo,&sn1,&org2,&des2,&ts2,&tf2);
                        if (sn1 != sn || org2 != org)
                            continue;     

                        if (ts2 >= tf) {
                            jcase--;
                            break;
                        }
                        if (grp >= 0) {
                            if (fabs(get_data(PM1VIdx[grp],jo)) <= EPSI1)
                                continue;
                        }
                        if (WIVar >= 0)             /* get weights */
                            wt1 = get_data(WIVar,jo) * WNorm;

                        rs += wt1;
   
                        j = next;                                                     
                        jv = -PIdx[j];
                        xa = TPar[j] * get_data(jv - 1,jo);
                        while ((jv = PIdx[++j]) > 0)  
                            xa += TPar[j] * get_data(jv - 1,jo);

                        rate = wt1 * rexp(xa);                                            
                        srate += rate;                                 

                        /* if we need first and second derivatives */
                                
                        if (sflg)   
                            fn_cox_s(next,jo,rate);
                    }
                }               /* end of update of the risk set if there are */
                                /* no time varying variables.                 */

                if (tf <= time)      
                    goto FN1CONT;

                /* If time-varying variables continue at FN1B */
                    
                if (TranTVar[tran])
                    goto FN1B;
                
                /* Reduction of the risk set in the case of no      */
                /* time-varying variables. Take all episodes of the */
                /* icase-loop which have tf2 < tf.                  */
                    
                for (ii = icase1; ii < icase; ++ii) {

                    jo = TFIdx[ii - 1];

                    get_edat(jo,&sn1,&org2,&des2,&ts2,&tf2);
                    if (sn1 != sn || org2 != org)
                        continue;     

                    if (tf2 >= tf)
                        break;

                    if (grp >= 0) {
                        if (fabs(get_data(PM1VIdx[grp],jo)) <= EPSI1)
                            continue;
                    }
                    if (WIVar >= 0)             /* get weights */
                        wt1 = get_data(WIVar,jo) * WNorm;

                    rs -= wt1;

                    j = next;                                                     
                    jv = -PIdx[j];
                    xa = TPar[j] * get_data(jv - 1,jo);
                    while ((jv = PIdx[++j]) > 0)  
                        xa += TPar[j] * get_data(jv - 1,jo);
                            
                    rate = -wt1 * rexp(xa);                                            
                    srate += rate;                                 

                    if (sflg)  
                        fn_cox_s(next,jo,rate);
                }
                icase1 = ii;    /* end of reduction of risk set if  */
                goto FN1C;      /* no time-varying covariates.      */

FN1B: ;         /* Build risk set at tf if time-varying variables.  */

                rs = srate = 0.0;
       
                if (sflg) {
                    dclear(NParm,WrkD);                                           
                    dclear(NParm,WrkE);                                           
                    dclear(HSiz,WrkH);                                          
                }
         
                /* Loop over all episodes according to ending times. */
                /* But there are possibly some episodes before icase */
                /* with the same ending time but different destinat. */
                /* state; they have to be included in the risk set.  */

                icase2 = icase;
                while (--icase2 > 0) {

                    jo = TFIdx[icase2 - 1];

                    get_edat(jo,&sn1,&org2,&des2,&ts2,&tf2);
                    /****************************/
                    if (sn1 != sn || org2 != org)
                        continue;     
                    /******************/
                    if (tf2 < tf)
                        break;
                }
                icase2++;
                for (jcase = icase2; jcase <= NOC; ++jcase) {

                    jo = TFIdx[jcase - 1];

                    get_edat(jo,&sn1,&org2,&des2,&ts2,&tf2);
                    if (sn1 != sn || org2 != org)
                        continue;     

                    /* starting time must be less than tf */      
                        
                    if (ts2 < tf) {  

                        if (grp >= 0) {
                            if (fabs(get_data(PM1VIdx[grp],jo)) <= EPSI1)
                                continue;
                        }
                        if (WIVar >= 0)             /* get weights */
                            wt1 = get_data(WIVar,jo) * WNorm;

                        rs += wt1;
                        j = next;                                                     
                        xa = 0.0;                                                     
                        jv = -PIdx[j];
                        while (1) {

                            /* If VTyp is 5 the variable */
                            /* is time-varying.          */
                         
                            if (VTyp[jv - 1] == 5) {

                                EDVALSn = sn1;      
                                EDVALOrg = org2;
                                EDVALDes = des2;
                                EDVALTs = ts2;
                                EDVALTf = tf2;
                                EDVALTime = tf;

                                err = v_eval1(jo,VESCnt[jv-1],VESTyp[jv-1],VESVal[jv-1],ESIdx,&x,0,0,0,0,0);
                                if (err) {
                                    NPFlgs[5] += 1;
                                    x = 0.0;     
                                }
                                put_data(x,jv-1,jo);
                            }
                            xa += TPar[j] * get_data(jv - 1,jo);
                            if ((jv = PIdx[++j]) < 0)
                                break;
                        }
                        rate = wt1 * rexp(xa);                                            
                        srate += rate;                                 
                        if (sflg)   
                            fn_cox_s(next,jo,rate);
                    }
                }               /* end of risk set building in the case of  */
                                /* time-varying variables.                  */
FN1C:
                /* If baseline rate calculation printing for time   */
                /* is done here.                                    */
     
                if (COXBLFlg) {
                    if (time == 0.0)
                        rs1 = rs;
                    prscox(tranfirst,tran,time,evnt,cens1,rs1,chrate);
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
         
                if (COXTCFlg && tpi < PMNTP && time >= PMTP[tpi]) {
                    ts_cox(tran,dim,tpi,nep - 1);
                    nep = 1;
                    while (++tpi < PMNTP && time >= PMTP[tpi])
                        ts_cox(tran,dim,tpi,nep - 1);
                }
           
FN1CONT:        if (COXBLFlg) { /* calculation of baseline rate */
                    evnt += wt;
                    if (srate > 0.0)  
                        chrate += wt / srate;
                    goto FN1Skip;
                }
          
                /* Now calculation of function value and gradient */
                /* and hessian if requested.                      */

                j = next;
                jv = -PIdx[j];
                xa = TPar[j] * get_data(jv - 1,io);
                while ((jv = PIdx[++j]) > 0)  
                    xa += TPar[j] * get_data(jv - 1,io);

                if (LFunc)  
                    FTmp += wt * (xa - rlog(srate));                             

                if (LGrad || COXTCFlg) {
                    j = next;
                    jv = -PIdx[j];                                                 
                    l = 1;
                    while (1) {
                        xj = get_data(jv - 1,io);
                        if (LGrad) { 
                            Grad[j] += wt * (xj - WrkD[j] / srate);                            
                            if (CGradFlg && io < MPGradRow) {
                                MatVal[MPGradIdx][io * MPGradCol + j] = 
                                                wt * (xj - WrkD[j] / srate);                            

                                if (PM2NV > 0)
                                    mp_putvar(io,MPGradCol,MatVal[MPGradIdx],
                                                             PM2NV,PM2VIdx,io);
                            }
                        }
                        if (COXTCFlg)
                            TCU[l++] += wt * (xj - WrkD[j] / srate);
         
                        if ((jv = PIdx[++j]) < 0 )
                            break;
                    }
                }
                if (LSec || COXTCFlg) {
                    j = next;
                    jk = (j * (j - 1)) / 2 + 1;                                   
                    jv = -PIdx[j];                                                 
                    tmp = WrkD[j] / srate;                                      
                    l = 1;
                    while (1) {
                        if (LSec)
                            Diag[j] += wt * (tmp * WrkD[j] - WrkE[j]) / srate;       
  
                        if (COXTCFlg)
                            TCV[(l - 1) * dim + l] -=
                                wt * (tmp * WrkD[j] - WrkE[j]) / srate;       
    
                        if ((jv = PIdx[++j]) < 0) 
                            break;
                        
                        k = next;
                        jk += next - 1; 
                        tmp = WrkD[j] / srate;                                      
                        l++;
                        l1 = 1;
                        while (k < j) {
                            if (LSec)
                                Hess[jk] +=
                                     wt * (tmp * WrkD[k] - WrkH[jk]) / srate;   
   
                            if (COXTCFlg)
                                TCV[(l - 1) * dim + l1] -=
                                  wt * (tmp * WrkD[k] - WrkH[jk]) / srate;   
    
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
  
            if (COXBLFlg) {
                prscox(tranfirst,tran,time,evnt,cens1,rs1,chrate);
                tranfirst = 0;
                PRTIdx++;
            }
        }           /* end of current group */

        /* If COXTCFlg calculate test statistics for last time periods */
  
        if (COXTCFlg) {
            ts_cox(tran,dim,tpi,nep);
            nep = 0;
            while (++tpi <= PMNTP)
                ts_cox(tran,dim,tpi,nep);

            TCDf[tran] -= dim;  /* adjust degrees of freedom */
        }
    }               /* end of main loop for all transitions     */
    NOCUsed = NOC;
    if (err == 0)
        err = checkov();
    return(err);
}

void fn_cox_s(int next,int jo,double  rate)
{
    register int j,k,jv,kv,jk;
    double xj,xk,xjrate;

    j = next;
    jk = (j * (j - 1)) / 2 + 1;                                       
    jv = -PIdx[j];                                                 
    xj = get_data(jv - 1,jo);
    xjrate = xj * rate;                                             

    while (1) {
    
        WrkD[j] += xjrate;                                     
    
        if (LSec || COXTCFlg)  
            WrkE[j] += xj * xjrate;                  
    
        if ((jv = PIdx[++j]) < 0)     
            break;
    
        xj = get_data(jv - 1,jo);
        xjrate = xj * rate;                                             
    
        if (LSec || COXTCFlg) {
            k = next;
            jk += next - 1; 
            kv = -PIdx[k];
            while (k < j) {
                xk = get_data(kv - 1,jo);
                WrkH[jk++] += xk * xjrate;                            
                kv = PIdx[++k];
            }
        }
    }
}

/****************************************************************************/
/*  ts_cox     Goodness-of-fit statistic for Cox models.                    */

void ts_cox(int tran,int dim,int tpi,int nep)
{
    register int i,j,r;
    double tmp,tmp1;

    if (PMProtFDef) {
        fprintf(PMProtFd,"\nTime period: %d\nUncensored episodes: %d\n",tpi,nep);
        prvec("TCU vector",dim,TCU);
        prmat("TCV matrix",dim,dim,TCV);
    }
    if (!nep) {         /* skip periods without episodes */ 
        TCDf[tran] = 0;
        return;
    }  
    TCDf[tran] += dim;      /* degrees of freedom */

    r = syminv(dim,TCV);    /* matrix inversion */
    if (PMProtFDef)  
        fprintf(PMProtFd,"Rank of V: %d\n",r);

    if (r == dim) {
        tmp = 0.0;
        for (i = 1; i <= dim; ++i) {
            tmp1 = 0.0;
            for (j = 1; j <= dim; ++j)
                tmp1 += TCU[j] * TCV[(j - 1) * dim + i];
            tmp += tmp1 * TCU[i];
        }
        TCStat[tran] += tmp;
        if (PMProtFDef)  
            fprintf(PMProtFd,"TStat added: %8.4lf\n",tmp);
    }
    else {
        TCStat[tran] = -1.0;
        if (PMProtFDef)  
            fprintf(PMProtFd,"TStat not calculated\n");
    }
    dclear(VTMax,TCU);
    dclear(VTMax * VTMax,TCV);
}

/****************************************************************************/
/*  test_cox    Goodness-of-fit statistic for Cox model.                    */
/*              Return 0 if OK, -1 if error.                                */

int test_cox(void)
{
    int err,n,tran,sn,org,des;
    double tmp;
   
    printf1("\nGlobal Goodness-of-fit.\n\n");

    err = -1;
    if (!(TCDf = (short *)calloc(MaxSnn * NTran1 + 1,sizeof(short))))  
        goto TSTCFin;
    TCDfA = MaxSnn * NTran1 + 1;
    memrq(TCDfA,sizeof(short));

    if (!(TCStat = (double *)calloc(MaxSnn * NTran1 + 1,sizeof(double))))  
        goto TSTCFin;
    TCStatA = MaxSnn * NTran1 + 1;
    memrq(TCStatA,sizeof(double));

    if (!(TCU = (double *)calloc(VTMax + 1,sizeof(double))))  
        goto TSTCFin;
    TCUA = VTMax + 1;
    memrq(TCUA,sizeof(double));

    if (!(TCV = (double *)calloc(VTMax * VTMax + 1,sizeof(double))))  
        goto TSTCFin;
    TCVA = VTMax * VTMax + 1;
    memrq(TCVA,sizeof(double));

    COXTCFlg = 1;
    LFunc = LGrad = LSec = 0; 

    fn_cox();   /*  calculations are done directly by fn_cox    */

    printf1("SN  Org  Des        TStat    DF   Signif\n");        
    prnchar('-',40,1);  

    for (tran = 0; tran < NTran1; ++tran) {

        sn  = SnTran1[tran];
        org = OrgTran1[tran];
        des = DesTran1[tran];

        printf1("%2d %4d %4d ",sn,org,des);          

        tmp = TCStat[tran];
        if (tmp > 0.0)
            printf1("%12.4lf ",tmp);
        else                          
            printf1("         --- ");
        n = TCDf[tran];
        printf1("%5d ",n);
        if (tmp > 0.0 && n > 0)
            printf1("%8.4lf\n",cdchif(tmp,n));
        else
            printf1("     ---\n");
    }
    newline();    
    err = 0;

TSTCFin:
    if (err)
        p_err(-2,1);

    if (TCDfA) {
        free((char *)TCDf);
        memrq(-TCDfA,sizeof(short));
        TCDfA = 0;
    }
    if (TCStatA) {
        free((char *)TCStat);
        memrq(-TCStatA,sizeof(double));
        TCStatA = 0;
    }
    if (TCUA) {
        free((char *)TCU);
        memrq(-TCUA,sizeof(double));
        TCUA = 0;
    }
    if (TCVA) {
        free((char *)TCV);
        memrq(-TCVA,sizeof(double));
        TCVA = 0;
    }
    COXTCFlg = 0;
    return(err);
}

/****************************************************************************/
/*  prscox_hd    Print header of prs table for Cox models.                  */
 
void prscox_hd(int sn,int org,int des,int grp) 
{
    fprintf(PRTFd,"\n# Transition: SN %d  Org %d  Des %d\n",sn,org,des);

    if (grp >= 0)
        fprintf(PRTFd,"# Group: %s\n",VName[PM1VIdx[grp]]);

    fprintf(PRTFd,"\n# ID ");
    fprnchar(PRTFd,' ',PMMFmt1 - 4,0); fprintf(PRTFd,"Time ");
    fprnchar(PRTFd,' ',PMMFmt1 - 6,0); fprintf(PRTFd,"Events ");
    fprnchar(PRTFd,' ',PMMFmt1 - 8,0); fprintf(PRTFd,"Censored ");
    fprnchar(PRTFd,' ',PMMFmt1 - 8,0); fprintf(PRTFd,"Risk Set ");
    fprnchar(PRTFd,' ',PMMFmt1 - 7,0); fprintf(PRTFd,"Surv.F. ");
    fprnchar(PRTFd,' ',PMMFmt1 - 8,0); fprintf(PRTFd,"Cum.Rate");
    fprintf(PRTFd,"   Baseline-Rate\n# ");
    fprnchar(PRTFd,'-',2 + 6 * (PMMFmt1 + 1) + 16,0);
    fprintf(PRTFd,"\n");
}

/****************************************************************************/
/*  prscox             print one line of prate table for Cox models.        */

void prscox(int tranfirst,int tran, double time, double evnt, double cens,
    double rs, double crate)
{
    double tmp,surv;
    static double prev_time = -1.0,prev_crate = -1.0;
     
    tmp = crate * AcW[tran];
    surv = rexp(-tmp);

    fprintf(PRTFd,"%4d ",PRTIdx);
    fprintf(PRTFd,PMMFmtS,time);
    fprintf(PRTFd,PMMFmtS,evnt);
    fprintf(PRTFd,PMMFmtS,cens);
    fprintf(PRTFd,PMMFmtS,rs);
    fprintf(PRTFd,PMMFmtS,surv);
    fprintf(PRTFd,PMMFmtS,tmp);

    tmp = 0.0;
    if (tranfirst)     
        prev_time = prev_crate = -1.0;
    else if (time > prev_time)
        tmp = (crate - prev_crate) / (time - prev_time);
    fprintf(PRTFd,"%15.10f ",tmp);
    fprintf(PRTFd,"\n");
   
    prev_time = time;
    prev_crate = crate;
}




