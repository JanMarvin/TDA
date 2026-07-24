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

/*  functions in t_ineq.c */

int ineq(void);
int segr(void);
double xgini(int n,double *x,double *wt);
int ds_comp(const void *arg1,const void *arg2);

/*--------------------------------------------------------------------------*/
/*  xgini(n,x,wt)   Calculate Gini coefficient for x[i], i = 0,...,n-1.     */
/*                  Return: Gini.                                           */

double xgini(int n,double *x,double *wt)
{
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

int ineq(void)
{
    register int i,j,jj;
    int err,l,n;           
    double tmp,a,b;

    err = -1;
    if (check_cmd(1))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    printf1("Inequality measures. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,4,1))     /* get parameters */
        goto INEQFin;

    prn_cwt();      /* case weight information */
       
    if (alloc_acx(NOC))
        goto INEQFin;

    if (alloc_acw(NOC))
        goto INEQFin;

    if (PMFmtF == 0)
        pmfmt(10,4);
   
    if (PMMPParDef == 1) {      /* use AcTmp for mppar */
        if (alloc_actmp(PMNV * 7 + 1))
            goto INEQFin;
    }

    newline();
    prn_hvar();
    prn_hlabel();
    printf1("  Cases "); 
    prnchar(' ',PMFmt1 -  8,0); printf1(" Minimum"); 
    prnchar(' ',PMFmt1 -  8,0); printf1("  Maximum");
    prnchar(' ',PMFmt1 -  8,0); printf1("     Mean");
    prnchar(' ',PMFmt1 -  6,0); printf1(" StdDev"); 
    prnchar(' ',PMFmt1 -  6,0); printf1(" VCoeff"); 
    prnchar(' ',PMFmt1 -  4,0); printf1(" Gini"); 
    printf1("\n");
       
    l = VNameLen + 9; 
    if (VLabelLen > 0)
        l += VLabelLen + 1;
    prnchar('-',l + 6 * (PMFmt1 + 1) - 1,1);
              
    for (j = 0; j < PMNV; ++j) {
        jj = PMVIdx[j];
        prn_vname(jj);
        prn_vlabel(jj);

        l = n = 0;
        for (i = 0; i < NOC; ++i) {
            tmp = get_data(jj,i);
            if (tmp >= 0.0) {
                if (l++ == 0)  
                    a = b = tmp;
                else {
                    if (a > tmp)
                        a = tmp;
                    if (b < tmp)
                        b = tmp;
                }
                AcX[n] = tmp;
                if (WIVar >= 0)
                    AcW[n] = get_data(WIVar,i);
                else
                    AcW[n] = 1.0;
                n++;
            }
        }
        printf1("%7d ",n);
        printf1(PMFmtS,a);
        printf1(PMFmtS,b);

        if (PMMPParDef == 1) {
            AcTmp[j * 7 + 1] = (double)n;    
            AcTmp[j * 7 + 2] = a;    
            AcTmp[j * 7 + 3] = b;    
        }

        if (PMF1Def) {
            fprintf(PMF1d,"%5d %7d ",j + 1,n);
            fprintf(PMF1d,PMFmtS,a);
            fprintf(PMF1d,PMFmtS,b);
        }
        b = xstd(n,AcX,AcW,&tmp);      /* mean and std. deviation */
        if (b > 0.0)                   /* var coefficient */
            a = tmp / b;
        else
            a = 0.0;

        printf1(PMFmtS,b);
        printf1(PMFmtS,tmp);
        printf1(PMFmtS,a);

        if (PMMPParDef == 1) {
            AcTmp[j * 7 + 4] = b;    
            AcTmp[j * 7 + 5] = tmp;    
            AcTmp[j * 7 + 6] = a;    
        }

        if (PMF1Def) {
            fprintf(PMF1d,PMFmtS,b);
            fprintf(PMF1d,PMFmtS,tmp);
            fprintf(PMF1d,PMFmtS,a);
        }
        tmp = xgini(n,AcX,AcW);        /* gini coefficient */
        printf1(PMFmtS,tmp);

        if (PMMPParDef == 1)  
            AcTmp[j * 7 + 7] = tmp;    

        if (PMF1Def)  
            fprintf(PMF1d,PMFmtS,tmp);

        printf1("\n");
        if (PMF1Def)  
            fprintf(PMF1d,"\n");
    }
    if (PMF1Def)  
        printf1("\nData written to: %s\n",PMF1dName);

    if (PMMPParDef == 1) {
        mp_putmpar(PMNV,7,AcTmp);
        newline();
        mp_info();
    }
    err = 0;

INEQFin:
    p_clean();      
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  segr()      Segregation measures.                                       */
/*  ##          segr (g=, v= [,df=,fmt=] ) [=fname]                         */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int segr(void)
{
    register int i,j,k,jj;
    int err,l,n,nc;
    double tmp,tmp1,d,g,gs0,gs1;

    err = -1;
    if (check_cmd(1))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    printf1("Segregation measures. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,4,1))     /* get parameters */
        goto SEGRFin;

    if (PMGIdx < 0) {
        printf1("Error: no grouping variable.\n");
        goto SEGRFin;
    }
    printf1("Grouping variable: %s\n",VName[PMGIdx]);

    prn_cwt();      /* case weight information */
    
    if (alloc_ack(NOC)) 
        goto SEGRFin;

    if (alloc_acn(NOC)) 
        goto SEGRFin;

    if (alloc_acm(NOC)) 
        goto SEGRFin;

    if (alloc_acc(NOC)) 
        goto SEGRFin;

    if (alloc_acw(NOC)) 
        goto SEGRFin;

    if (alloc_acx(NOC)) 
        goto SEGRFin;

    if (alloc_acy(NOC)) 
        goto SEGRFin;

    if (PMFmtF == 0)
        pmfmt(10,4);

    if (PMTFmtF == 0)
        pmtfmt(10,4);

    newline();
    prn_hvar();
    prn_hlabel();
    printf1(" Classes      Cases    Group_0    Group_1");
    prnchar(' ',PMFmt1 -  6,0); printf1("D-Index"); 
    prnchar(' ',PMFmt1 -  6,0); printf1("V-Ratio"); 
    prnchar(' ',PMFmt1 -  3,0); printf1("Gini"); 
    printf1("\n");
       
    l = VNameLen + 43; 
    if (VLabelLen > 0)
        l += VLabelLen + 1;
    prnchar('-',l + 3 * (PMFmt1 + 1) - 1,1);
              
    for (j = 0; j < PMNV; ++j) {
        jj = PMVIdx[j];
        prn_vname(jj);
        prn_vlabel(jj);

        n = 0;
        for (i = 0; i < NOC; ++i) {

            g = (float)get_data(PMGIdx,i);
            l = (int)get_data(jj,i);

            if (g >= 0.0 && l >= 0) {
                AcK[n] = l;                     
                if (g > 0)
                    AcC[n] = 1;
                else
                    AcC[n] = 0;

                if (WIVar >= 0)
                    AcW[n] = get_data(WIVar,i);
                else
                    AcW[n] = 1.0;
                n++;
            }   
        }
        if (n < 1) {
            printf1("number of cases is zero.\n");
            continue;
        }

        /* get ptr to sort vidx values */

        for (i = 0; i < n; ++i)  
            AcN[i] = i;
   
        qsort((char *)AcN,n,sizeof(int),ds_comp);

        for (i = 0; i < n; ++i)
            AcX[i] = AcY[i] = 0.0;

        k = AcN[0];
        l = AcK[k];
        AcM[k] = l;
        if (AcC[k])
            AcY[0] = AcW[k];
        else
            AcX[0] = AcW[k];
        nc = 0;

        for (i = 1; i < n; ++i) {
            k = AcN[i];
            if (AcK[k] != l) {
                l = AcK[k];
                nc++;
                AcM[nc] = l;
                if (AcC[k])
                    AcY[nc] = AcW[k];
                else
                    AcX[nc] = AcW[k];
            }
            else {
                if (AcC[k])
                    AcY[nc] += AcW[k];
                else
                    AcX[nc] += AcW[k];
            }
        }
        nc++;
        gs0 = gs1 = 0.0;
        for (i = 0; i < nc; ++i) {
            gs0 += AcX[i];
            gs1 += AcY[i];
        }
        g = gs0 + gs1;

        /* ## print to second output file if requested */

        if (PMTabFDef) {
            fprintf(PMTabFd,"Variable: %s\n",VName[jj]);
            fprintf(PMTabFd,"Class      Cases    Group_0    Group_1    D-Index\n");
            d = 0.0;                                                          
            for (i = 0; i < nc; ++i) {
                fprintf(PMTabFd,"%5d %10.2f %10.2f %10.2f ",AcM[i],
                                        AcX[i] + AcY[i],AcX[i],AcY[i]);

                if (gs0 > 0.0 && gs1 > 0.0) {
                    tmp = 0.5 * fabs((AcX[i] / gs0) - (AcY[i] / gs1));
                    fprintf(PMTabFd,"%10.4f",tmp);
                    d += tmp;
                }
                fprintf(PMTabFd,"\n");
            }
            fprintf(PMTabFd,"Total %10.2f %10.2f %10.2f ",g,gs0,gs1);
            if (gs0 > 0.0 && gs1 > 0.0)  
                fprintf(PMTabFd,"%10.4f",d);
            fprintf(PMTabFd,"\n");
        }

        printf1(" %7d %10.2f %10.2f %10.2f ",nc,g,gs0,gs1);

        if (PMF1Def)  
            fprintf(PMF1d,"%4d %7d %10.2f %10.2f %10.2f ",j + 1,nc,g,gs0,gs1);

        if (gs0 <= 0.0 || gs1 <= 0) {
            newline();         
            if (PMF1Def)
                fprintf(PMF1d,"\n");
            continue;
        }

        /* calculate dissimilarity index */

        d = 0.0;
        for (i = 0; i < nc; ++i)  
            d += fabs((AcX[i] / gs0) - (AcY[i] / gs1));
        d *= 0.5;

        printf1(PMFmtS,d);
        if (PMF1Def)
            fprintf(PMF1d,PMFmtS,d);

        /* calculate variance ratio index */

        d = 0.0;
        for (i = 0; i < nc; ++i) {
            tmp1 = AcX[i] + AcY[i];
            if (tmp1 > 0.0) {
                tmp = gs1 / g;
                if (AcY[i] > 0.0)
                    tmp -= AcY[i] / tmp1;
                d += tmp * tmp * tmp1;
            }
        }
        d /= (gs0 * gs1 / g);
        printf1(PMFmtS,d);
        if (PMF1Def)
            fprintf(PMF1d,PMFmtS,d);

        /* calculate gini index */

        d = 0.0;
        for (i = 0; i < nc; ++i) {
            tmp = AcX[i] + AcY[i];           
            if (tmp > 0.0) {
                for (k = 0; k < nc; ++k) {
                    if (k != i) {
                        tmp1 = AcX[k] + AcY[k];           
                        if (tmp1 > 0.0)  
                            d += tmp * tmp1 *
                                       fabs((AcY[i] / tmp) - (AcY[k] / tmp1));
                    }
                }
            }
        }
        d = 0.5 * d / (gs0 * gs1);
        printf1(PMFmtS,d);
        if (PMF1Def)
            fprintf(PMF1d,PMFmtS,d);

        newline();
        if (PMF1Def)
            fprintf(PMF1d,"\n");
    }
    if (PMF1Def || PMTabFDef)
        newline();
    if (PMF1Def)  
        printf1("Data written to: %s\n",PMF1dName);
    if (PMTabFDef)  
        printf1("Additional data written to: %s\n",PMTabFName);

    err = 0;

SEGRFin:
    p_clean(); 
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ds_comp()   compare function for AcK                                    */

int ds_comp(const void *arg1,const void *arg2)
{     
    int n; 

    n = AcK[*(int *)arg1] - AcK[*(int *)arg2];   
    if (n > 0)
        return(1);
    else if (n < 0)
        return(-1);
    return(0);
}

