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

/*  functions in t_dstat.c */

int dstat(void);
int quant(void);
int mfreq(int typ);
int ct_stat(int opt,int m,int *value,int *freq,int *bf,int widx,double *wfreq,
    int ridx,int cidx);
int pcov(int typ);
int atab(void);
int dma(void);
void dma_prn1(int m,int n,double *x,int p);
void dma_prn2(int m,int n,double *x,int *vidx);
void dma_prn3(FILE *fd,int m,int n,double *x,int p);
void dma_scal(int n,double *d,double *x);
void dma_stand(int n,int m,double *x,int opt);
void dma_corr(int n,int m,double *x,double *a);
void dma_pcf(int nr,double *rx,int nc,double *cx);

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

int dstat(void)
{
    register int i,j;
    int err,ig,ign,ix,n,l,w,nn,nrow;
    double xmin,xmax,xm,xsd,xsum,wsum,wt,tmp;

    err = -1;
    if (check_cmd(1))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    printf1("Descriptive statistics. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,4,0))     /* get parameters */
        goto DSFin;

    if (PMNV == 0) {        /* use all variables */

        if (pm_valloc(NVAR))  
            goto DSFin;
             
        w = j = 0;
        i = VIFirst;
        while (i >= 0) {
            PMVIdx[j++] = i;
            if (VTyp[i] == 5)
                w = 1;
            i = VNxt[i];
        }
        if (w)
            p_warn(-3,1);
    }   
    if (PMFmtF == 0)
        pmfmt(10,4);

    if (alloc_acx(NOC + 1))
        goto DSFin;

    if (SVEFlg) {           /* info about case selection */
        prn_sve();
        if (WIVar >= 0) {
            printf1("\nCannot use case weights with sel option.\n");
            goto DSFin;
        }
    }
    prn_cwt();              /* info about case weights */

    if (PMMPParDef == 1) {      /* use AcTmp for mppar */
        n = nn = PMNV * 6;
        nrow = PMNV;
        if (PM1NV > 0) {
            n *= PM1NV;
            nrow *= PM1NV;
        }
        if (alloc_actmp(n + 1))
            goto DSFin;
    }
    ig = -1;
    ign = 0;

DSNXT:
    if (PM1NV > 0) {
        ig = PM1VIdx[ign];
        printf1("\nGroup: %s\n",VName[ig]);

        if (PMF1Def)  
            fprintf(PMF1d,"# Group: %s\n",VName[ig]);

    }         
    newline();
    prn_hvar();
    prn_hlabel();
    prnchar(' ',PMFmt1 -  8,0); printf1(" Minimum"); 
    prnchar(' ',PMFmt1 -  8,0); printf1("  Maximum");
    prnchar(' ',PMFmt1 -  8,0); printf1("     Mean");
    prnchar(' ',PMFmt1 -  8,0); printf1(" Std.Dev.");
    printf1("    Sum of values\n");
       
    l = VNameLen + 18;
    if (VLabelLen > 0)
        l += VLabelLen + 1;
    prnchar('-',l + 4 * (PMFmt1 + 1) - 1,1);

    for (j = 0; j < PMNV; ++j) {

        ix = PMVIdx[j];
        n = getd1(AcX,ix,ig,0);

        xmin = xmax = xm = xsd = xsum = 0.0;
        if (n > 0) {
            xmin = DBLMAX;
            xmax = -xmin;  

            wt = 1.0;
            if (WIVar >= 0)
                wsum = 0.0;
            else
                wsum = (double)n;

            for (i = 0; i < n; ++i) {
                tmp = AcX[i];
                if (WIVar >= 0) {       
                    wt = get_data(WIVar,i) * WNorm;
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
                    tmp = AcX[i];

                    if (WIVar >= 0)         
                        wt = get_data(WIVar,i) * WNorm;
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
        prn_vname(ix);
        prn_vlabel(ix);
        printf1(PMFmtS,xmin);
        printf1(PMFmtS,xmax);
        printf1(PMFmtS,xm);
        printf1(PMFmtS,xsd);
        printf1("%16.4f\n",xsum);

        if (PMF1Def) {
            fprintf(PMF1d,"%4d ",j + 1);
            fprintf(PMF1d,PMFmtS,xmin);
            fprintf(PMF1d,PMFmtS,xmax);
            fprintf(PMF1d,PMFmtS,xm);
            fprintf(PMF1d,PMFmtS,xsd);
            fprintf(PMF1d,"%16.4f\n",xsum);
        }
        if (PMMPParDef == 1) {
            AcTmp[ign * nn + j * 6 + 1] = (double)ign;
            AcTmp[ign * nn + j * 6 + 2] = xmin;
            AcTmp[ign * nn + j * 6 + 3] = xmax;
            AcTmp[ign * nn + j * 6 + 4] = xm;
            AcTmp[ign * nn + j * 6 + 5] = xsd;
            AcTmp[ign * nn + j * 6 + 6] = xsum;
        }
    }
    if (++ign < PM1NV)  
        goto DSNXT;
       
    if (PMMPParDef == 1) {
        mp_putmpar(nrow,6,AcTmp);
        newline();
        mp_info();
    }
    err = 0;

DSFin:
    p_clean();
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

int quant(void)
{
    register int i;
    int err,n,iv,ivv,nqx,len;
    double tmp,qx[11];

    err = -1;         
    if (check_cmd(1))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    printf1("Quantiles. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,4,1))     /* get parameters */
        goto QUANTFin;

    if (alloc_acx(NOC))     /* allocate AcX */
        goto QUANTFin;

    if (SVEFlg)             /* info about case selection */
        prn_sve();

    if (PMMPParDef == 1) {      /* use AcTmp for mppar */
        if (alloc_actmp(PMNV * 11 + 1))
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

    if (PMFmtF == 0)  
        pmfmt(7,2);
      
    len = VNameLen;
    if (len < PMFmt1)
        len = PMFmt1;

    printf1("\nVariable ");            
    prnchar(' ',len - 8,0);
    for (i = 0; i < nqx; ++i)  
        printf1(PMFmtS,qx[i]);
    newline();
    prnchar('-',len + nqx * (PMFmt1 + 1),1);

    for (iv = 0; iv < PMNV; ++iv) {         /* for all variables */ 

        ivv = PMVIdx[iv];

        n = getd1(AcX,ivv,-1,0);

        printf1("%s ",VName[ivv]);
        prnchar(' ',len - strlen(VName[ivv]),0);

        if (n == 0) {
            printf1(" number of cases is zero.\n");
            continue;       
        }
        if (sortd(n,AcX,0))
            goto QUANTFin;

        for (i = 0; i < nqx; ++i) {
            tmp = quantf(n,AcX,qx[i]);
            printf1(PMFmtS,tmp);
            if (PMF1Def)
                fprintf(PMF1d,PMFmtS,tmp);

            if (PMMPParDef == 1)  
                AcTmp[iv * 11 + i + 1] = tmp;          
        }
        newline();
        if (PMF1Def)
            fprintf(PMF1d,"\n");
    }
    if (PMMPParDef == 1) {
        mp_putmpar(PMNV,11,AcTmp);
        newline();
        mp_info();
    }
    err = 0;

QUANTFin:
    p_clean(); 
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  mfreq(typ)  Process CmdDef[idx] with syntax                             */
/*              typ = 0 : freq  (joint freq distribution)                   */
/*              typ = 1 : freq1 (separate freq distributions)               */
/*              typ = 2 : freq2 (2-dim table, optional with cont. meas.)    */
/*                                                                          */
/*              Return 0 if successful, otherwise -1.                       */

int mfreq(int typ)
{
    register int i,j,k,l;
    int err,nvv,vl,len,iv,kv,r;
    double f,sf,tmp,stmp,sum;

    err = -1;         
    if (check_cmd(1))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    printf1("Frequency tables. Current memory: %d bytes.\n",MemReq);
    i = 0;
    if (typ > 0)
        i = 1;

    if (parm(CmdBuf + 4 + i,4,1))     /* get parameters */
        goto MFFin;
  
    if (typ == 2 && PMNV != 2) {
        printf1("Error: freq2 command needs exactly two variables.\n");
        goto MFFin;
    }

    printf1("Maximum number of categories: %d\n",PMMaxCat);

    /* allocate memory */
   
    if (alloc_ack(PMMaxCat * PMNV + 1))     /* value */
        goto MFFin;

    if (alloc_acn(PMMaxCat + 1))            /* freq */
        goto MFFin;

    if (alloc_aci(PMMaxCat + 1))            /* bf pointer */
        goto MFFin;

    if (WIVar >= 0) {   /* weight command active */

        if (alloc_acw(PMMaxCat + 1))    
            goto MFFin;
        prn_cwt();                  /* case weight info */
    }
             
    /* prepare print format */

    if (PMFmtF == 0)
        PMFmt1 = 3;

    vl = get_vnl(PMNV,PMVIdx);
    if (vl < PMFmt1)
        vl = PMFmt1;

    if (vl > PMFmt1)
        pmfmt(vl,0);
    else if (PMFmtF == 0)
        pmfmt(PMFmt1,0);

    for (i = 0; i < PMNV; ++i) {      /* make all distributions */

        printf1("\nFrequency distribution for variable(s): %s",VName[PMVIdx[i]]);

        if (typ == 1) {
            r = cfreq1(1,PMVIdx + i,PMMaxCat,AcK,AcN,AcI,WIVar,AcW);
        }
        else {
            for (j = 1; j < PMNV; ++j)  
                printf1(",%s",VName[PMVIdx[j]]);
                     
            r = cfreq1(PMNV,PMVIdx,PMMaxCat,AcK,AcN,AcI,WIVar,AcW);

        }
        newline();

        if (r <= 0) {
            if (r == -1)
                printf1("Max number of categories is too small.\n");
            else if (r == -2)
                printf1("Stack size is too small.\n");
            else if (r == -3)
                printf1("Insufficient memory for frequency distribution.\n");
            else 
                printf1("Error %d in frequency distribution.\n",r);
            goto MFFin;
        }
        printf1("Number of categories: %d\n",r);
        sum = 0.0;
        for (j = 1; j <= r; ++j) {
            if (WIVar >= 0)
                sum += AcW[j];
            else 
                sum += (double)AcN[j];
        }
        if (sum == EPSI)  
            printf1("Sum of (weighted) frequencies is zero.\n");
        else if (typ < 2) {

            len = 42;
            printf1("\nIndex ");
            if (typ == 1) {
                nvv = 1;
                iv = PMVIdx[i];
                prnchar(' ',vl - strlen(VName[iv]),0);
                printf1(" %s",VName[iv]);
                len += vl + 1;
            }
            else {
                nvv = PMNV;
                for (j = 0; j < PMNV; ++j) {
                    iv = PMVIdx[j];
                    prnchar(' ',vl - strlen(VName[iv]),0);
                    printf1(" %s",VName[iv]);
                    len += vl + 1;
                }
            }
            printf1("  Frequency   Pct   Cumulated   Pct\n");
            prnchar('-',len,1);

            sf = stmp = 0.0;

            for (j = 1; j <= r; ++j) {
                printf1("%5d  ",j);
                if (PMF1Def)
                    fprintf(PMF1d,"%5d  ",j);

                k = AcI[j];
                for (l = 1; l <= nvv; ++l) {
                    kv = AcK[(k - 1) * nvv + l];
                    printf1(PMFmtS,(double)kv);
                    if (PMF1Def)
                        fprintf(PMF1d,PMFmtS,(double)kv);
                }
                if (WIVar >= 0)
                    f = AcW[k];
                else 
                    f = (double)AcN[k];

                tmp = 100.0 * f / sum;
                sf += f;
                stmp += tmp;

                printf1("%10.2f %6.2f %10.2f %6.2f\n",f,tmp,sf,stmp);
                if (PMF1Def)
                    fprintf(PMF1d,"%10.2f %6.2f %10.2f %6.2f\n",f,tmp,sf,stmp);
            }
            prnchar('-',len,1);
            printf1("Sum   ");
            prnchar(' ',len - 42,0);
            printf1(" %10.2f %6.2f\n",sum,100.0);
            if (PMF1Def)  
                printf1("\nTable written to: %s\n",PMF1dName);
        }
        else {          /* freq2 command */

            ct_stat((int)PMSC,r,AcK,AcN,AcI,WIVar,AcW,(int)PMVIdx[0],(int)PMVIdx[1]);
        }
        if (typ != 1)
            break;
    }
    err = 0;

MFFin:
    p_clean();  
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

int ct_stat(int opt,int m,int *value,int *freq,int *bf,int widx,double *wfreq,
    int ridx,int cidx)
{
    register int i,j,k,jj,kk;
    int r,c,n,n1,n2,df,rc,err,nflag;
    double lr,s,e,a,a1,a2,b1,b2,sum,t,tmp;

    err = -1;

    if (alloc_acr(m + 1))
        goto CTFin;

    if (alloc_acs(m + 1))
        goto CTFin;

    if (alloc_acm(m + 1))
        goto CTFin;

    if (alloc_acj(m + 1))
        goto CTFin;


    /* sort categories for columns */

    for (i = 1; i <= m; ++i)
        AcS[i] = value[i * 2];  

    c = cfreq(1,m,AcS,m,AcR,AcM,AcJ,0,&a,&tmp);
    if (c < 1) {
        printf1("Error: can't sort categories.\n");
        goto CTFin;
    }
    for (i = 1; i <= c; ++i)
        AcS[i - 1] = AcR[AcJ[i]];

    /* sort categories for rows */

    r = 0;
    k = bf[1];
    n = value[(k - 1) * 2 + 1] - 1;
    for (i = 1; i <= m; ++i) {
        k = bf[i];
        j = value[(k - 1) * 2 + 1];
        if (j > n) {
            AcR[r] = j;
            n = j;
            r++;
        }
    }

    rc = r * c;

    if (alloc_actmp(rc))        /* used to store the table */
        goto CTFin;

    if (alloc_acx(r))           /* row sum */
        goto CTFin;

    if (alloc_acy(c))           /* column sum */
        goto CTFin;

    nflag = 0;
    sum = 0.0;
    for (i = 1; i <= m; ++i) {
        n = (bf[i] - 1) * 2;
        jj = value[++n];
        kk = value[++n];
        for (j = 0; j < r; ++j) {
            if (AcR[j] == jj)
                break;
        }
        for (k = 0; k < c; ++k) {
            if (AcS[k] == kk)
                break;
        }
        if (widx >= 0)
            tmp = wfreq[bf[i]];
        else  
            tmp = (double)freq[bf[i]];

        AcTmp[j * c + k] += tmp;                   
        sum += tmp;
        AcX[j] += tmp;
        AcY[k] += tmp;
        if (tmp < 0.0)
            nflag++;
    }

    /* print the table */

    printf1("\nFrequency|\n");
    printf1("Percent  |\n");
    printf1("Row Pct  |\n");
    printf1("Col Pct  |");
    for (k = 0; k < c; ++k)  
        printf1("%8d |",AcS[k]);
    printf1("   Total\n");
    prnchar('-',(c + 2) * 10 - 2,1);
    for (j = 0; j < r; ++j) {
        printf1("%8d |",AcR[j]);
        for (k = 0; k < c; ++k) {
            printf1("%8.2f |",AcTmp[j * c + k]);
            if (PMF1Def)
                fprintf(PMF1d,"%8.2f ",AcTmp[j * c + k]);
        }
        if (PMF1Def)
            fprintf(PMF1d,"\n");

        printf1("%8.2f\n         |",AcX[j]);
        for (k = 0; k < c; ++k) {
            if (sum != 0.0)
                tmp = 100.0 * AcTmp[j * c + k] / sum;
            else
                tmp = 0.0;
            printf1("%8.2f |",tmp);
        }
        if (sum != 0.0)
            tmp = 100.0 * AcX[j] / sum;
        else
            tmp = 0.0;
        printf1("%8.2f\n         |",tmp);
        for (k = 0; k < c; ++k) {
            if (AcX[j] != 0.0) 
                tmp = 100.0 * AcTmp[j * c + k] / AcX[j];
            else
                tmp = 0.0;
            printf1("%8.2f |",tmp);
        }
        printf1("\n         |");
        for (k = 0; k < c; ++k) {
            if (AcY[k] != 0.0)
                tmp = 100.0 * AcTmp[j * c + k] / AcY[k];
            else
                tmp = 0.0;
            printf1("%8.2f |",tmp);
        }
        newline();
        prnchar('-',(c + 2) * 10 - 2,1);
    }
    printf1("Total     ");
    for (k = 0; k < c; ++k)  
        printf1("%8.2f  ",AcY[k]);
    printf1("%8.2f\n          ",sum);
    for (k = 0; k < c; ++k) {
        if (sum != 0.0)
            tmp = 100.0 * AcY[k] / sum;
        else
            tmp = 0.0;
        printf1("%8.2f  ",tmp);
    }
    printf1("  100.00\n");
    if (opt == 0)
        goto CTFin;

    /* ----------------------------------------------------------------- */

    newline();
    prnchar('-',LLEN,1);
    printf1("Contingency measures for table: X (%s) x Y (%s)\n\n",
                                                   VName[ridx],VName[cidx]);
    if (r < 2) {
        printf1("Table has less than two categories in: %s\n",VName[ridx]);
        goto CTFin;
    }
    if (c < 2) {
        printf1("Table has less than two categories in: %s\n",VName[cidx]);
        goto CTFin;
    }
    if (nflag > 0) {
        printf1("Error: table has %d negative elements.\n",nflag);
        goto CTFin;
    }
    if (sum <= EPSI1) {
        printf1("Error: sum of table elements is zero or less.\n");
        goto CTFin;
    }

    /*  chi square */

    err = n1 = n2 = 0;
    lr = s = 0.0;
    for (j = 0; j < r; ++j) {
        for (k = 0; k < c; ++k) {
            t = AcTmp[j * c + k];
            e = AcX[j] * AcY[k] / sum;
            if (t < 5.0)
                n1++;
            if (e < 5.0)
                n2++;

            if (e <= EPSI1)
                err++;
            else {
                s += (t - e) * (t - e) / e;
                lr += t * rlog(t / e);
            }
        }
    }
    df = (r - 1) * (c - 1);
    a = 100.0 * (double)n1 / (double)rc;
    printf1("Number of cells with less than 5 elements: %3d (%4.1f %%)\n",
                                                                       n1,a);
    a = 100.0 * (double)n2 / (double)rc;
    printf1("Cells with expected frequency less than 5: %3d (%4.1f %%)\n",     
                                                                       n2,a);
    printf1("Degrees of freedom: %d\n\n",df);
    if (err > 0) {
        printf1("Can't calculate Chi-square.\n");
        err = -1;
        goto CTFin;
    }
    prn_sfmt("Chi-square (Pearson)",35,PMTFmtS,s);
    a = 1.0 - cdchif(s,df);
    prn_sfmt("Prob:",35,PMTFmtS,a);
    lr *= 2.0;
    prn_sfmt("Chi-square (likelihood ratio)",35,PMTFmtS,lr);
    a =  1.0 - cdchif(lr,df);
    prn_sfmt("Prob:",35,PMTFmtS,a);

    if (s > EPSI1) {
        a = sqrt(s / sum);
        newline();
        prn_sfmt("Phi",35,PMTFmtS,a);
        n = r;
        if (n > c)
            n = c;
        if (n > 1) {
            a = sqrt(s / (sum * (double)(n - 1)));
            prn_sfmt("Cramer's V",35,PMTFmtS,a);
        }
        if (s + sum != 0.0 && (tmp = s / (s + sum)) > EPSI1) {
            a = sqrt(tmp);
            prn_sfmt("Contingency coefficient",35,PMTFmtS,a);
        }

        /* lamda's */

        a1 = 0.0;
        for (k = 0; k < c; ++k) {
            tmp = AcTmp[k];
            for (j = 1; j < r; ++j) {
                t = AcTmp[j * c + k];
                if (tmp < t)
                    tmp = t;
            }
            a1 += tmp;
        }
        b1 = AcX[0];
        for (j = 1; j < r; ++j) {
            if (b1 < AcX[j])
                b1 = AcX[j];
        }
        if (sum > b1) {
            s = (a1 - b1) / (sum - b1);
            if (s > 0.0)  
                prn_sfmt("Lambda (X dependent)",35,PMTFmtS,s);
        }

        a2 = 0.0;
        for (j = 0; j < r; ++j) {
            tmp = AcTmp[j * c];
            for (k = 1; k < c; ++k) {
                t = AcTmp[j * c + k];
                if (tmp < t)
                    tmp = t;
            }
            a2 += tmp;
        }
        b2 = AcY[0];
        for (j = 1; j < c; ++j) {
            if (b2 < AcY[j])
                b2 = AcY[j];
        }
        if (sum > b2) {
            s = (a2 - b2) / (sum - b2);
            if (s > 0.0)  
                prn_sfmt("Lambda (Y dependent)",35,PMTFmtS,s);
        }
        tmp = b1 + b2;
        if (2.0 * sum > tmp) {
            s = (a1 + a2 - tmp) / (2.0 * sum - tmp);        
            if (s > 0.0)  
                prn_sfmt("Lambda (symmetric)",35,PMTFmtS,s);
        }

        /* Goodman - Kruskal Tau */

        a1 = 0.0;
        for (k = 0; k < c; ++k) {
            if (AcY[k] < EPSI1) {
                a1 = -1.0;
                break;
            }
            for (j = 0; j < r; ++j) {
                tmp = AcTmp[j * c + k];
                a1 += tmp * tmp / AcY[k];
            }
        }
        if (a1 > 0.0) {
            b1 = 0.0;
            for (j = 0; j < r; ++j)
                b1 += AcX[j] * AcX[j];

            tmp = sum * sum - b1;
            if (tmp > 0.0) {
                tmp = (sum * a1 - b1) / tmp;
                if (tmp >= 0.0)  
                    prn_sfmt("Goodman/Kruskal Tau (X dependent)",35,PMTFmtS,tmp);
            }
        }
        a1 = 0.0;
        for (j = 0; j < r; ++j) {
            if (AcX[j] < EPSI1) {
                a1 = -1.0;
                break;
            }
            for (k = 0; k < c; ++k) {
                tmp = AcTmp[j * c + k];
                a1 += tmp * tmp / AcX[j];
            }
        }
        if (a1 > 0.0) {
            b1 = 0.0;
            for (k = 0; k < c; ++k)
                b1 += AcY[k] * AcY[k];

            tmp = sum * sum - b1;
            if (tmp > 0.0) {
                tmp = (sum * a1 - b1) / tmp;
                if (tmp >= 0.0)  
                    prn_sfmt("Goodman/Kruskal Tau (Y dependent)",35,PMTFmtS,tmp);
            }
        }

        /* Uncertainty coefficient */

        a1 = 0.0;
        for (j = 0; j < r; ++j) {
            if (AcX[j] > 0.0) {
                tmp = AcX[j] / sum;
                a1 -= tmp * rlog(tmp);  
            }
        }
        b1 = 0.0;
        for (k = 0; k < c; ++k) {
            if (AcY[k] > 0.0) {
                tmp = AcY[k] / sum;
                b1 -= tmp * rlog(tmp);  
            }
        }
        s = 0.0;
        for (j = 0; j < r; ++j) {
            for (k = 0; k < c; ++k) {
                tmp = AcTmp[j * c + k];
                if (tmp > 0.0) {
                    tmp /= sum;
                    s -= tmp * rlog(tmp);  
                }
            }
        }
        if (fabs(a1) > EPSI1) {
            tmp = (a1 + b1 - s) / a1;
            if (tmp > 0.0)  
                prn_sfmt("Uncertainty coeff (X dependent)",35,PMTFmtS,tmp);
        }
        if (fabs(b1) > EPSI1) {
            tmp = (a1 + b1 - s) / b1;
            if (tmp > 0.0)  
                prn_sfmt("Uncertainty coeff (Y dependent)",35,PMTFmtS,tmp);
        }

        /* Kendall's tau, etc */

        a1 = b1 = 0.0;
        for (j = 0; j < r; ++j) {
            for (k = 0; k < c; ++k) {
                tmp = AcTmp[j * c + k];
                if (tmp > 0.0) {
                    t = 0.0;
                    for (jj = 0; jj < j; ++jj) {
                        for (kk = 0; kk < k; ++kk)  
                            t += AcTmp[jj * c + kk];
                    }
                    for (jj = j + 1; jj < r; ++jj) {
                        for (kk = k + 1; kk < c; ++kk)  
                            t += AcTmp[jj * c + kk];
                    }
                    a1 += tmp * t;

                    t = 0.0;
                    for (jj = 0; jj < j; ++jj) {
                        for (kk = k + 1; kk < c; ++kk)  
                            t += AcTmp[jj * c + kk];
                    }
                    for (jj = j + 1; jj < r; ++jj) {
                        for (kk = 0; kk < k; ++kk)  
                            t += AcTmp[jj * c + kk];
                    }
                    b1 += tmp * t;
                }
            }
        }
        a2 = sum * sum;
        for (j = 0; j < r; ++j)
            a2 -= AcX[j] * AcX[j];

        b2 = sum * sum;
        for (k = 0; k < c; ++k)
            b2 -= AcY[k] * AcY[k];

        tmp = a2 * b2;
        if (tmp > EPSI1) {
            s = (a1 - b1) / sqrt(tmp);
            if (s >= 0.0)  
                prn_sfmt("Kendall's tau-b",35,PMTFmtS,s);
        }
        if (n > 1) {            /* n = min(r,c) */
            tmp = sum * sum * ((double)(n - 1));
            if (tmp > EPSI1) {
                s = (double)n * (a1 - b1) / tmp;
                if (s >= 0.0)  
                    prn_sfmt("Kendall's tau-c",35,PMTFmtS,s);
            }
        }
        tmp = a1 + b1;
        if (tmp > EPSI1) {
            s = (a1 - b1) / tmp;
            if (s >= 0.0)  
                prn_sfmt("Gamma",35,PMTFmtS,s);
        }
        if (b2 > EPSI1) {
            s = (a1 - b1) / b2;
            if (s >= 0.0)  
                prn_sfmt("Somers' D (X dependent)",35,PMTFmtS,s);
        }
        if (a2 > EPSI1) {
            s = (a1 - b1) / a2;
            if (s >= 0.0)  
                prn_sfmt("Somers' D (Y dependent)",35,PMTFmtS,s);
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

int pcov(int typ)
{
    register int i,j,k;
    int nn,vl,err,w;

    err = -1;
    if (check_cmd(1))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    if (typ == 0)
        printf1("Covariance");
    else if (typ == 1)
        printf1("Correlation");
    else if (typ == 2)
        printf1("Rank correlation");

    printf1(" matrix. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 3 + typ,4,0))     /* get parameters */
        goto COVFin;

    if (PMNV == 0) {        /* use all variables */

        if (pm_valloc(NVAR))  
            goto COVFin;
           
        w = j = 0;
        i = VIFirst;
        while (i >= 0) {
            PMVIdx[j++] = i;
            if (VTyp[i] == 5)
                w = 1;
            i = VNxt[i];
        }
        if (w)
            p_warn(-3,1);
    }   
    if (PMNV < 2) {
        p_err(-8,1);
        goto COVFin;
    }
    if (typ <= 1)
        prn_cwt();          /* case weight information */

    /* allocate memory */

    nn = PMNV * (PMNV + 1) / 2;
    if (alloc_acx(nn + 1))
        goto COVFin;

    if (alloc_acy(PMNV + 1))
        goto COVFin;

    if (PMMPCovDef == 1) {      /* use AcTmp for mpcov */
        if (alloc_actmp(PMNV * PMNV + 1))
            goto COVFin;
    }

    if (typ == 0) {
        if (cov(PMNV,PMVIdx,AcY,AcX))       /* get covariance matrix */
            goto COVFin;
    }
    else if (typ == 1) {
        if (alloc_acu(PMNV + 1))
            goto COVFin;

        if (corr(PMNV,PMVIdx,AcY,AcU,AcX))  /* get correlation matrix */
            goto COVFin;
    }
    else if (typ == 2) {
        rcorr(PMNV,PMVIdx,AcX);             /* get rank correlation matrix */
    }
   
    if (PMFmtF == 0)
        pmfmt(7,4);

    vl = get_vnl(PMNV,PMVIdx);
    if (vl < PMFmt1)
        vl = PMFmt1;

    newline();
    prnchar(' ',vl,0);
    k = vl;

    for (j = 0; j < PMNV; ++j) {
        prnchar(' ',vl - strlen(VName[PMVIdx[j]]),0);
        printf1("%s ",VName[PMVIdx[j]]);
        k += vl + 1;
    }
    newline();    
    prnchar('-',k - 1,1);
    k = 0;
    for (i = 0; i < PMNV; ++i) {

        printf1("%s",VName[PMVIdx[i]]);
        prnchar(' ',vl - strlen(VName[PMVIdx[i]]),0);

        for (j = 0; j <= i; ++j) {
            printf1(PMFmtS,AcX[k]);

            if (PMMPCovDef == 1) {
                AcTmp[i * PMNV + j + 1] =               
                AcTmp[j * PMNV + i + 1] = AcX[k];       
            }
            k++;
        }
        newline();
    }
    if (PMF1Def) {      /* write to output file */

        prn_f1mat(PMF1d,PMPRNO,PMFmtS,PMNV,AcX);
        printf1("\nMatrix written to: %s\n",PMF1dName);
    }
    if (PMMPCovDef == 1) {
        mp_putcov(PMNV,AcTmp);
        newline();
        mp_info();
    }
    err = 0;

COVFin:
    p_clean();      
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

int atab(void)
{
    register int i,j,k;
    int n,l,k1,k2,ix,iy,err,nrow,ncol;
    double wt,tmp,xmin,xmax,ymin,ymax;

    err = -1;
    if (check_cmd(1))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    printf1("Aggregated table. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,4,1))     /* get parameters */
        goto ATFin;
   
    n = PMNTP;      /* number of arguments in PMTP[] (t=...) */
    if (n < 1) {
        printf1("Error: need classes defined with x parameter.\n");
        goto ATFin;
    }
    if (PMNV == 1) {
        ix = PMVIdx[0];
        iy = -1;
    }
    else if (PMNV == 2) {
        ix = PMVIdx[0];
        iy = PMVIdx[1];
    }
    else {
        printf1("Error: need one or two variables.\n");
        goto ATFin;
    }

    /* allocate memory */

    if (alloc_acn(n + 1))
        goto ATFin;

    if (alloc_acw(n + 1))
        goto ATFin;

    if (alloc_acx(n + 1))
        goto ATFin;

    if (PMNV == 2) {
        if (alloc_acy(n + 1))
            goto ATFin;
    }
    prn_cwt();          /* case weight information */

    xmin = xmax = get_data(ix,0);
    if (iy >= 0)
        ymin = ymax = get_data(iy,0);

    wt = 1.0;
    for (i = 0; i < NOC; ++i) {

        if (WIVar >= 0)  
            wt = get_data(WIVar,i) * WNorm;
          
        tmp = get_data(ix,i);

        if (xmin > tmp)
            xmin = tmp;
        if (xmax < tmp)
            xmax = tmp;

        k = n;

        for (j = 0; j < n; ++j) {
            if (tmp < PMTP[j]) {
                k = j;
                break;
            }
        }
        AcN[k] += 1;
        AcW[k] += wt;
        AcX[k] += tmp * wt;
        if (iy >= 0) {
            tmp = get_data(iy,i);
            AcY[k] += tmp * wt;
            if (ymin > tmp)
                ymin = tmp;
            if (ymax < tmp)
                ymax = tmp;
        }
    }
    for (k = 0; k <= n; ++k) {
        tmp = AcW[k];
        if (fabs(tmp) >= EPSI) {
            AcX[k] /= tmp;
            if (iy >= 0)
                AcY[k] /= tmp;
        }
    }
    if (PMFmtF == 0)
        pmfmt(10,4);
  
    printf1("\nX: ");            
    prn_vname(ix);
    printf1(" Minimum: %12.4f  Maximum: %12.4f\n",xmin,xmax);
    if (iy >= 0) {
        printf1("Y: ");         
        prn_vname(iy);
        printf1(" Minimum: %12.4f  Maximum: %12.4f\n",ymin,ymax);
    }
    newline();

    l = 4;
    prnchar(' ',PMFmt1 -  5,0); printf1("Lower "); 
    prnchar(' ',PMFmt1 -  5,0); printf1("Upper "); 
                                printf1(" Frequency"); 
    prnchar(' ',PMFmt1 -  8,0); printf1(" Weighted"); 
    prnchar(' ',PMFmt1 -  8,0); printf1("  Mean(X)"); 
    if (iy >= 0) {
        prnchar(' ',PMFmt1 -  8,0); printf1("  Mean(Y)");             
        l++;
    }
    newline();
    prnchar('-',11 + l * (PMFmt1 + 1) - 1,1);

    k1 = 0;
    for (k = 0; k <= n; ++k) {
        if (AcN[k]) {
            k1 = k;
            break;
        }
    }
    k2 = n;
    for (k = n; k >= 0; --k) {
        if (AcN[k]) {
            k2 = k;
            break;
        }
    }
   
    if (PMMPParDef == 1) {      /* use AcTmp for mppar */
        nrow = k2 - k1 + 1;
        ncol = 5;
        if (iy >= 0)
            ncol++;
        if (alloc_actmp(nrow * ncol + 1))
            goto ATFin;
    }
    nrow = 0;
    for (k = k1; k <= k2; ++k) {

        if (PMR != 0 && AcN[k] == 0)
            continue;

        if (k > 0) { 
            printf1(PMFmtS,PMTP[k - 1]);
            if (PMMPParDef == 1)  
                AcTmp[nrow * ncol + 1] = PMTP[k - 1];  
        }
        else
            prnchar(' ',PMFmt1 + 1,0);

        if (k < n) {
            printf1(PMFmtS,PMTP[k]);
            if (PMMPParDef == 1)  
                AcTmp[nrow * ncol + 2] = PMTP[k];  
        }
        else
            prnchar(' ',PMFmt1 + 1,0);

        printf1("%10d ",AcN[k]);
        printf1(PMFmtS,AcW[k]);
        printf1(PMFmtS,AcX[k]);

        if (PMMPParDef == 1) {
            AcTmp[nrow * ncol + 3] = (double)AcN[k];
            AcTmp[nrow * ncol + 4] = AcW[k];
            AcTmp[nrow * ncol + 5] = AcX[k];
        }   
        if (iy >= 0) {
            printf1(PMFmtS,AcY[k]);
            if (PMMPParDef == 1)  
                AcTmp[nrow * ncol + 6] = AcY[k];
        }
        nrow++;
        newline();
    }
    if (PMMPParDef == 1) {
        mp_putmpar(nrow,ncol,AcTmp);
        newline();
        mp_info();
    }
    err = 0;

ATFin:
    p_clean();  
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

int dma(void)
{
    register int i,j,k;
    int err,nc,nv,nnc,nnv,nv1,r,vlen,wrec,tflag;
    double s,tmp,tmp1,tmp2,tsum;   

    wrec = 0;
    tsum = 0.0;
    err = -1;
    if (check_cmd(2))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    printf1("Data matrix analysis. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 3,1,0))     /* get parameters */
        goto DMAFin;

    if (PMFmtF == 0 || PMFmt1 < 8)
        pmfmt(10,4);
   
    if (PMNV > 0)
        nv = PMNV;
    else
        nv = NVAR;

    if (PMALG < 1 || PMALG > 7)
        PMALG = 1;
    printf1("\nAlgorithm %d: ",PMALG);
    if (PMALG == 1)
        printf1("principal components (based on data matrix).\n");
    else if (PMALG == 2)
        printf1("principal components (based on cov./corr. matrix).\n");
    else if (PMALG == 3)
        printf1("factor analysis (based on data matrix).\n");
    else if (PMALG == 4)
        printf1("factor analysis (based on correlation matrix).\n");
    if (PMALG == 5)
        printf1("dual scaling (based on frequency table).\n");
    else if (PMALG == 6)
        printf1("correspondence analysis (based on frequency table).\n");
    else if (PMALG == 7)
        printf1("direct svd-based projection.\n");
    
    printf1("Number of rows (cases): %d\n",NOC);
    printf1("Number of columns (variables): %d\n",nv);

    if (PMALG == 2 || PMALG == 4) {
        if (NOC != nv) {
            printf1("Error: need a symmetric input matrix.\n");
            goto DMAFin;
        }
    }
    else if (PMALG == 7) {
        if (NOC < nv) {
            printf1("Error: number of rows is less than number of columns.\n");
            goto DMAFin;
        }
    }
    newline();

    if (alloc_aci(nv))
        goto DMAFin;

    if (PMNV > 0) {
        for (j = 0; j < PMNV; ++j)
            AcI[j] = PMVIdx[j];
    }
    else { 
        j = 0;
        i = VIFirst;
        while (i >= 0) {
            AcI[j++] = i;
            i = VNxt[i];
        }
    }
    vlen = 8;
    for (j = 0; j < nv; ++j) 
        vlen = imax(vlen,strlen(VName[AcI[j]]));

    nc = NOC;
    if (alloc_acx(nc * nv + 1))
        goto DMAFin;

    for (i = 0; i < nc; ++i) {
        for (j = 0; j < nv; ++j) 
            AcX[i * nv + j + 1] = get_data(AcI[j],i);
    }

    /* ------------------------ alg 1 and 2 ------------------------------- */

    if (PMALG == 1 || PMALG == 2) {

        if (alloc_acy(nv + 1))            /* used for eigenvalues */
            goto DMAFin;

        if (alloc_acv(nv * nv + 1))       /* used for eigenvectors */
            goto DMAFin;

        if (PMALG == 1) {

            if (PMOPT == 2) {
                printf1("Preprocessing: mean centering of variables.\n");
                dma_stand(nc,nv,AcX,0);
            }
            else if (PMOPT == 3) {
                printf1("Preprocessing: standardization of variables.\n");
                dma_stand(nc,nv,AcX,1);
            }
            if (PMPRNO) {
                printf1("Input data (after preprocessing [opt=%d])\n",PMOPT);
                dma_prn1(nc,nv,AcX,-1);
                newline();
            }
            dma_corr(nc,nv,AcX,AcV);
        }
        else if (PMALG == 2) {
            for (i = 1; i <= nv * nv; ++i)
                AcV[i] = AcX[i];
        }
        if (PMPRNO) {
            printf1("Cross-product matrix\n");
            dma_prn1(nv,nv,AcV,-1);
            newline();
        }
        r = evecf(nv,AcY,AcV);
        if (r) {
            printf1("No success in eigenvalue calculation.\n");
            goto DMAFin;
        }
        printf1("Eigenvalue         per cent\n");
        tmp = 0.0;
        for (j = 1; j <= nv; ++j)
            tmp += AcY[j];

        for (j = 1; j <= nv; ++j) {
            prnchar(' ',16 - PMFmt1,0);
            printf1(PMFmtS,AcY[j]);
            if (tmp > 0.0)
                printf1("  %8.2f",100.0 * AcY[j] / tmp);
            newline();
        }
        newline();
        if (PMPRNO) {
            printf1("Eigenvectors\n");
            dma_prn1(nv,nv,AcV,PMNS);
        }
        newline();

        if (PMF1Def) {
            dma_prn3(PMF1d,nv,nv,AcV,PMNS);
            printf1("Eigenvectors (%d records) written to: %s\n",nv,PMF1dName);
        }

        if (PMALG == 1 && PMFDef) {  /* princ. components */

            if (PMNS >= 1 && PMNS <= nv)
                r = PMNS;
            else
                r = nv;

            for (i = 0; i < nc; ++i) {
                for (j = 1; j <= r; ++j) {
                    tmp = 0.0;
                    for (k = 1; k <= r; ++k) 
                        tmp += AcX[i * nv + k] * AcV[(k - 1) * nv + j];
                    fprintf(PMFd,PMFmtS,tmp);
                }
                fprintf(PMFd,"\n");
            }
            printf1("Principal components (%d records) written to: %s\n",nc,PMFdName);
        }
    }

    /* ------------------------ alg 3 and 4 ------------------------------- */

    else if (PMALG == 3 || PMALG == 4) {

        if (alloc_acy(nv + 1))            /* used for eigenvalues */
            goto DMAFin;

        if (alloc_acv(nv * nv + 1))       /* used for eigenvectors */
            goto DMAFin;

        if (PMALG == 3) {                 /* create correlation matrix */
            dma_stand(nc,nv,AcX,1);
            dma_corr(nc,nv,AcX,AcV);
        }
        else if (PMALG == 4) {
            for (i = 1; i <= nv * nv; ++i)
                AcV[i] = AcX[i];
        }
        if (PMPRNO) {
            printf1("Correlation matrix\n");
            dma_prn1(nv,nv,AcV,-1);
            newline();
        }
        r = evecf(nv,AcY,AcV);
        if (r) {
            printf1("No success in eigenvalue calculation.\n");
            goto DMAFin;
        }
        printf1("Eigenvalue         per cent\n");
        tmp = 0.0;
        for (j = 1; j <= nv; ++j)
            tmp += AcY[j];

        for (j = 1; j <= nv; ++j) {
            prnchar(' ',16 - PMFmt1,0);
            printf1(PMFmtS,AcY[j]);
            if (tmp > 0.0)
                printf1("  %8.2f",100.0 * AcY[j] / tmp);
            newline();
        }
        newline();

        if (PMALG == 3 && PMFDef) {  /* factors */

            if (PMNS >= 1 && PMNS <= nv)
                r = PMNS;
            else
                r = nv;

            for (i = 0; i < nc; ++i) {
                for (j = 1; j <= r; ++j) {
                    tmp = 0.0;
                    for (k = 1; k <= r; ++k) 
                        tmp += AcX[i * nv + k] * AcV[(k - 1) * nv + j];
                    fprintf(PMFd,PMFmtS,tmp);
                }
                fprintf(PMFd,"\n");
            }
            wrec = nc;
        }
        dma_scal(nv,AcY,AcV);
        if (PMPRNO) {
            printf1("Rescaled eigenvectors/factor loadings\n");
            dma_prn1(nv,nv,AcV,PMNS);
        }
        newline();
        if (PMF1Def) {
            dma_prn3(PMF1d,nv,nv,AcV,PMNS);
            printf1("Eigenvectors/factor loadings (%d records) written to: %s\n",nv,PMF1dName);
        }
        if (wrec > 0)
            printf1("Factors (%d records) written to: %s\n",nc,PMFdName);
    }
 
    /* -##--------------------- alg 5 ------------------------------- */

    else if (PMALG == 5) {
       
        if (PMPRNO) {
            printf1("Frequency table\n");
            dma_prn2(nc,nv,AcX,AcI);
            newline();
        }

        if (PMOPT == 2 || PMOPT == 3) {
            if (PMOPT == 2) {
                for (i = 0; i < nc; ++i) {
                    tmp = 0.0;
                    for (j = 1; j <= nv; ++j)
                        tmp += AcX[i * nv + j];
                    if (tmp != 0.0) {
                        for (j = 1; j <= nv; ++j)
                            AcX[i * nv + j] /= tmp;
                    }
                }
            }
            else if (PMOPT == 3) {
                tmp = 0.0;
                for (i = 1; i <= nc * nv; ++i)  
                    tmp += AcX[i];
                if (tmp != 0.0) {
                    for (i = 1; i <= nc * nv; ++i)
                        AcX[i] /= tmp;
                }
            }

            if (PMPRNO) {
                printf1("Input data after preprocessing (opt=%d)\n",PMOPT);
                dma_prn1(nc,nv,AcX,-1);
                newline();
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
                    AcX[i * nnv + j + 1] = get_data(AcI[i],j);
            }
        }
        if (alloc_acz(nnc + 1))           /* used for row sums */
            goto DMAFin;
        if (alloc_acw(nnv + 1))           /* used for column sums */
            goto DMAFin;

        for (i = 1; i <= nnc; ++i) {
            for (j = 1; j <= nnv; ++j) {
                tmp = AcX[(i - 1) * nnv + j];
                AcZ[i] += tmp;
                AcW[j] += tmp;
            }
        }
        for (i = 1; i <= nnc; ++i) {
            for (j = 1; j <= nnv; ++j) {
                tmp = AcZ[i] * AcW[j];
                if (tmp <= 0.0) {
                    printf1("Error: table has zero row or column.\n");
                    goto DMAFin;
                }
                AcX[(i - 1) * nnv + j] /= sqrt(tmp);
            }
        }
        if (PMF1Def) {
            for (i = 0; i < nnc; ++i) {
                for (j = 1; j <= nnv; ++j)  
                    fprintf(PMF1d,PMFmtS,AcX[i * nnv + j]);
                fprintf(PMF1d,"\n");
            }
        }

        if (alloc_acy(nnv + 1))   
            goto DMAFin;
        if (alloc_acu(nnc * nnv + 1))  
            goto DMAFin;
        if (alloc_acv(nnv * nnv + 1))  
            goto DMAFin;

        r = svdecomp(nnc,nnv,AcX,AcY,AcU,AcV,3);
        if (r) {
            p_err(-2,1);
            goto DMAFin;
        }                        
        printf1("Singular value         eigenvalue  per cent\n");
        tmp = 0.0;
        for (j = 1; j <= nnv; ++j)
            tmp += AcY[j] * AcY[j];

        for (j = 1; j <= nnv; ++j) {
            s = AcY[j];
            prnchar(' ',16 - PMFmt1,0);
            printf1(PMFmtS,s);
            s *= s;
            prnchar(' ',16 - PMFmt1,0);
            printf1(PMFmtS,s);
            if (tmp > 0.0)
                printf1(" %8.2f",100.0 * s / tmp);
            newline();
        }

        printf1("\nRow ");
        for (j = 1; j <= nnv; ++j) {
            prnchar(' ',PMFmt1 - 8,0);
            printf1("score%3d ",j);
        }
        newline();
        for (i = 1; i <= nc; ++i) {
            printf1("%3d ",i);
            for (j = 1; j <= nnv; ++j) {
                if (tflag == 0)
                    printf1(PMFmtS,AcU[(i - 1) * nnv + j] / sqrt(AcZ[i]));
                else
                    printf1(PMFmtS,AcV[(i - 1) * nnv + j] / sqrt(AcW[i]));
            }
            newline();
        }

        printf1("\nCol ");
        for (j = 1; j <= nnv; ++j) {
            prnchar(' ',PMFmt1 - 8,0);
            printf1("score%3d ",j);
        }
        newline();
        for (i = 1; i <= nv; ++i) {
            printf1("%3d ",i);
            for (j = 1; j <= nnv; ++j) {
                if (tflag == 0)
                    printf1(PMFmtS,AcV[(i - 1) * nnv + j] / sqrt(AcW[i]));
                else
                    printf1(PMFmtS,AcU[(i - 1) * nnv + j] / sqrt(AcZ[i]));
            }
            newline();
        }
        if (PMF1Def)
            printf1("\nA-matrix (%d records) written to: %s\n\n",nnc,PMF1dName);
    }

    /* ------------------------ alg 6 ------------------------------- */

    else if (PMALG == 6) {

        if (PMPRNO) {
            printf1("Frequency table\n");
            dma_prn2(nc,nv,AcX,AcI);
            newline();
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
                    AcX[i * nnv + j + 1] = get_data(AcI[i],j);
            }
        }
        tsum = 0.0;
        for (i = 1; i <= nc * nv; ++i)
            tsum += AcX[i];

        if (tsum <= EPSI1 || nc < 2 || nv < 2) {
            printf1("Error: invalid table.\n");
            goto DMAFin;
        }
        for (i = 1; i <= nc * nv; ++i)
            AcX[i] /= tsum;

        if (alloc_acz(nnc + 1))           /* used for row sums */
            goto DMAFin;
        if (alloc_acw(nnv + 1))           /* used for column sums */
            goto DMAFin;

        for (i = 1; i <= nnc; ++i) {
            for (j = 1; j <= nnv; ++j) {
                tmp = AcX[(i - 1) * nnv + j];
                AcZ[i] += tmp;
                AcW[j] += tmp;
            }
        }
        for (i = 1; i <= nnc; ++i) {
            for (j = 1; j <= nnv; ++j) {
                tmp = AcZ[i] * AcW[j];
                if (tmp <= 0.0) {
                    printf1("Error: table has zero row or column.\n");
                    goto DMAFin;
                }
                AcX[(i - 1) * nnv + j] =
                      (AcX[(i - 1) * nnv + j] - AcZ[i] * AcW[j]) / sqrt(tmp);
            }
        }
        if (PMPRNO) {
            printf1("Standardized residuals\n");
            dma_prn1(nc,nv,AcX,-1);
            newline();
        }

        if (alloc_acy(nnv + 1))   
            goto DMAFin;
        if (alloc_acu(nnc * nnv + 1))  
            goto DMAFin;
        if (alloc_acv(nnv * nnv + 1))  
            goto DMAFin;

        r = svdecomp(nnc,nnv,AcX,AcY,AcU,AcV,3);
        if (r) {
            p_err(-2,1);
            goto DMAFin;
        }                        
        printf1("Singular value         eigenvalue  per cent\n");
        tmp = 0.0;
        for (j = 1; j <= nnv; ++j)
            tmp += AcY[j] * AcY[j];

        for (j = 1; j <= nnv; ++j) {
            s = AcY[j];
            prnchar(' ',16 - PMFmt1,0);
            printf1(PMFmtS,s);
            s *= s;
            prnchar(' ',16 - PMFmt1,0);
            printf1(PMFmtS,s);
            if (tmp > 0.0)
                printf1(" %8.2f",100.0 * s / tmp);
            newline();
        }
        if (PMPCFDef) {
            if (alloc_actmp(2 * nc + 1))   
                goto DMAFin;
            if (alloc_actmp1(2 * nv + 1))   
                goto DMAFin;
        }

        printf1("\nRow ");
        prnchar(' ',PMFmt1 - 4,0); printf1("Mass ");
        prnchar(' ',PMFmt1 - 3,0); printf1("PC1 ");
        prnchar(' ',PMFmt1 - 3,0); printf1("PC2 ");
        prnchar(' ',PMFmt1 - 3,0); printf1("SC1 ");
        prnchar(' ',PMFmt1 - 3,0); printf1("SC2\n");

        for (i = 1; i <= nc; ++i) {
            printf1("%3d ",i);
            if (tflag == 0) {
                tmp = AcZ[i];
                tmp1 = AcU[(i - 1) * nnv + 1] / sqrt(AcZ[i]);
                tmp2 = AcU[(i - 1) * nnv + 2] / sqrt(AcZ[i]);
            }
            else {
                tmp = AcW[i];
                tmp1 = AcV[(i - 1) * nnv + 1] / sqrt(AcW[i]);
                tmp2 = AcV[(i - 1) * nnv + 2] / sqrt(AcW[i]);
            }
            printf1(PMFmtS,tmp);
            printf1(PMFmtS,tmp1 * AcY[1]);
            printf1(PMFmtS,tmp2 * AcY[2]);
            printf1(PMFmtS,tmp1);
            printf1(PMFmtS,tmp2);
            newline();

            if (PMPCFDef) {
                AcTmp[(i - 1) * 2 + 1] = tmp1 * AcY[1];
                AcTmp[(i - 1) * 2 + 2] = tmp2 * AcY[2];
            }
            if (PMF1Def) {
                fprintf(PMF1d,PMFmtS,tmp);
                fprintf(PMF1d,PMFmtS,tmp1 * AcY[1]);
                fprintf(PMF1d,PMFmtS,tmp2 * AcY[2]);
                fprintf(PMF1d,PMFmtS,tmp1);
                fprintf(PMF1d,PMFmtS,tmp2);
                fprintf(PMF1d,"\n");
            }
        }

        printf1("\nCol ");
        prnchar(' ',PMFmt1 - 4,0); printf1("Mass ");
        prnchar(' ',PMFmt1 - 3,0); printf1("PC1 ");
        prnchar(' ',PMFmt1 - 3,0); printf1("PC2 ");
        prnchar(' ',PMFmt1 - 3,0); printf1("SC1 ");
        prnchar(' ',PMFmt1 - 3,0); printf1("SC2\n");

        for (i = 1; i <= nv; ++i) {
            printf1("%3d ",i);
            if (tflag == 0) {
                tmp = AcW[i];
                tmp1 = AcV[(i - 1) * nnv + 1] / sqrt(AcW[i]);
                tmp2 = AcV[(i - 1) * nnv + 2] / sqrt(AcW[i]);
            }
            else {
                tmp = AcZ[i];
                tmp1 = AcU[(i - 1) * nnv + 1] / sqrt(AcZ[i]);
                tmp2 = AcU[(i - 1) * nnv + 2] / sqrt(AcZ[i]);
            }
            printf1(PMFmtS,tmp);
            printf1(PMFmtS,tmp1 * AcY[1]);
            printf1(PMFmtS,tmp2 * AcY[2]);
            printf1(PMFmtS,tmp1);
            printf1(PMFmtS,tmp2);
            newline();

            if (PMPCFDef) {
                AcTmp1[(i - 1) * 2 + 1] = tmp1 * AcY[1];
                AcTmp1[(i - 1) * 2 + 2] = tmp2 * AcY[2];
            }

            if (PMFDef) {
                fprintf(PMFd,PMFmtS,tmp);
                fprintf(PMFd,PMFmtS,tmp1 * AcY[1]);
                fprintf(PMFd,PMFmtS,tmp2 * AcY[2]);
                fprintf(PMFd,PMFmtS,tmp1);
                fprintf(PMFd,PMFmtS,tmp2);
                fprintf(PMFd,"\n");
            }
        }
        newline();
        if (PMF1Def)
            printf1("Row coordinates (%d records) written to: %s\n",nc,PMF1dName);
        if (PMFDef)
            printf1("Column coordinates (%d records) written to: %s\n",nv,PMFdName);
        if (PMPCFDef)
            dma_pcf(nc,AcTmp,nv,AcTmp1);

    }

    /* --###------------------- alg 7 ------------------------------- */

    else if (PMALG == 7) {

        if (PMPRNO) {
            printf1("Input data\n");
            dma_prn1(nc,nv,AcX,-1);
            newline();
        }
                              
        if (PMOPT == 2 || PMOPT == 3) {
            if (PMOPT == 2) {
                for (i = 0; i < nc; ++i) {
                    tmp = 0.0;
                    for (j = 1; j <= nv; ++j)
                        tmp += AcX[i * nv + j];
                    if (tmp != 0.0) {
                        for (j = 1; j <= nv; ++j)
                            AcX[i * nv + j] /= tmp;
                    }
                }
            }
            else if (PMOPT == 3) {
                tmp = 0.0;
                for (i = 1; i <= nc * nv; ++i)  
                    tmp += AcX[i];
                if (tmp != 0.0) {
                    for (i = 1; i <= nc * nv; ++i)
                        AcX[i] /= tmp;
                }
            }
        }
        if (PMPRNO) {
            printf1("Input data after preprocessing (opt=%d)\n",PMOPT);
            dma_prn1(nc,nv,AcX,-1);
            newline();
        }

        if (alloc_acy(nv + 1))   
            goto DMAFin;
        if (alloc_acu(nc * nv + 1))  
            goto DMAFin;
        if (alloc_acv(nv * nv + 1))  
            goto DMAFin;

        r = svdecomp(nc,nv,AcX,AcY,AcU,AcV,3);
        if (r) {
            p_err(-2,1);
            goto DMAFin;
        }                        
        printf1("Singular value         eigenvalue  per cent\n");
        tmp = 0.0;
        for (j = 1; j <= nv; ++j)
            tmp += AcY[j] * AcY[j];

        for (j = 1; j <= nv; ++j) {
            s = AcY[j];
            prnchar(' ',16 - PMFmt1,0);
            printf1(PMFmtS,s);
            s *= s;
            prnchar(' ',16 - PMFmt1,0);
            printf1(PMFmtS,s);
            if (tmp > 0.0)
                printf1(" %8.2f",100.0 * s / tmp);
            newline();
        }
        nv1 = imin(nv,2);
        printf1("\nCoordinates\n");
        for (i = 1; i <= nc; ++i) {
            for (j = 1; j <= nv1; ++j) {
                tmp = 0.0;
                for (k = 1; k <= nv; ++k) 
                    tmp += AcX[(i - 1) * nv + k] * AcV[(k - 1) * nv + j];
                printf1(PMFmtS,tmp);

                if (PMFDef)  
                    fprintf(PMFd,PMFmtS,tmp);
            }   
            newline();
            if (PMFDef)  
                fprintf(PMFd,"\n");
        }
        newline();
        if (PMFDef)
            printf1("Coordinates (%d records) written to: %s\n",nc,PMFdName);
    }
    err = 0;

DMAFin:
    p_clean();  
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  dma_prn1(m,n,x,p)     print (m,n)-matrix to stdout                      */

void dma_prn1(int m,int n,double *x,int p)
{
    register int i,j;

    if (p < 1 || p > n)
        p = n;

    for (i = 0; i < m; ++i) {
        for (j = 1; j <= p; ++j)
            printf1(PMFmtS,x[i * n + j]);
        newline();
    }
}

/* ------------------------------------------------------------------------ */
/*  dma_prn2(m,n,x,vidx)   print table to stdout                            */

void dma_prn2(int m,int n,double *x,int *vidx)
{
    register int i,j,l;
    double tmp,sum;

    if (alloc_actmp(n)) {
        p_err(-2,1);
        return;
    }
    printf1(" Row ");
    for (j = 0; j < n; ++j) {
        l = strlen(VName[vidx[j]]);
        prnchar(' ',PMFmt1 - l,0);
        printf1("%s ",VName[vidx[j]]);
    }
    prnchar(' ',PMFmt1 - 3,0);
    printf1("Sum\n");

    for (i = 0; i < m; ++i) {
        printf1("%4d ",i + 1);
        sum = 0.0;
        for (j = 0; j < n; ++j) {
            l = strlen(VName[vidx[j]]);
            prnchar(' ',l - PMFmt1,0);
            tmp = x[i * n + j + 1];
            printf1(PMFmtS,tmp);
            sum += tmp;
            AcTmp[j] += tmp;
        }
        printf1(PMFmtS,sum);
        newline();
    }
    printf1(" Sum ");
    sum = 0.0;
    for (j = 0; j < n; ++j) {
        l = strlen(VName[vidx[j]]);
        prnchar(' ',l - PMFmt1,0);
        tmp = AcTmp[j];
        printf1(PMFmtS,tmp);
        sum += tmp;
    }
    printf1(PMFmtS,sum);
    newline();
}

/* ------------------------------------------------------------------------ */
/*  dma_prn3(fd,m,n,x)   print (m,n)-matrix to fd                           */

void dma_prn3(FILE *fd,int m,int n,double *x,int p)
{
    register int i,j;

    if (p < 1 || p > n)
        p = n;

    for (i = 0; i < m; ++i) {
        for (j = 1; j <= p; ++j)
            fprintf(fd,PMFmtS,x[i * n + j]);
        fprintf(fd,"\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  dma_scal(n,d,x)   rescale x with d                                      */

void dma_scal(int n,double *d,double *x)
{
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

void dma_stand(int n,int m,double *x,int opt)
{
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

void dma_corr(int n,int m,double *x,double *a)
{
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

void dma_pcf(int nr,double *rx,int nc,double *cx)
{
    register int i;
    double xmin,xmax,ymin,ymax,dx,dy,d;

    xmin = xmax = rx[1];
    for (i = 2; i <= nr; ++i) {
        xmin = dmin(xmin,rx[(i - 1) * 2 + 1]);
        xmax = dmax(xmax,rx[(i - 1) * 2 + 1]);
    }
    for (i = 1; i <= nc; ++i) {
        xmin = dmin(xmin,cx[(i - 1) * 2 + 1]);
        xmax = dmax(xmax,cx[(i - 1) * 2 + 1]);
    }

    ymin = ymax = rx[2];
    for (i = 2; i <= nr; ++i) {
        ymin = dmin(ymin,rx[(i - 1) * 2 + 2]);
        ymax = dmax(ymax,rx[(i - 1) * 2 + 2]);
    }
    for (i = 1; i <= nc; ++i) {
        ymin = dmin(ymin,cx[(i - 1) * 2 + 2]);
        ymax = dmax(ymax,cx[(i - 1) * 2 + 2]);
    }
    dx = xmax - xmin;
    dy = ymax - ymin;
    xmin -= 0.1 * dx;
    xmax += 0.1 * dx;
    ymin -= 0.1 * dy;
    ymax += 0.1 * dy;
    d = dx / 40.0;

    fprintf(PMPCFd,"psfile = %s.ps;\n",PMPCFName);
    fprintf(PMPCFd,"psetup(\n");
    fprintf(PMPCFd,"    pxa=%8.4f,%8.4f,\n",xmin,xmax);
    fprintf(PMPCFd,"    pya=%8.4f,%8.4f,\n",ymin,ymax);
    fprintf(PMPCFd,"    pxlen=100,\n");
    fprintf(PMPCFd,"    pylen=100);\n");
    
    fprintf(PMPCFd,"plotp=%8.4f,%8.4f,%8.4f,%8.4f;\n",0.0,ymin,0.0,ymax);
    fprintf(PMPCFd,"plotp=%8.4f,%8.4f,%8.4f,%8.4f;\n",xmin,0.0,xmax,0.0);
  
    fprintf(PMPCFd,"plotp(s=4,fs=2,lt=0)=\n");
    for (i = 1; i <= nr; ++i) {
        fprintf(PMPCFd,"    %8.4f,%8.4f",rx[(i - 1) * 2 + 1],rx[(i - 1) * 2 + 2]);
        if (i < nr)
            fprintf(PMPCFd,",\n");
        else
            fprintf(PMPCFd,";\n");
    }
    fprintf(PMPCFd,"plotp(s=5,fs=2,lt=0)=\n");
    for (i = 1; i <= nc; ++i) {
        fprintf(PMPCFd,"    %8.4f,%8.4f",cx[(i - 1) * 2 + 1],cx[(i - 1) * 2 + 2]);
        if (i < nc)
            fprintf(PMPCFd,",\n");
        else
            fprintf(PMPCFd,";\n");
    }
    for (i = 1; i <= nr; ++i) {
        fprintf(PMPCFd,"pltext(fs=2,xy=%8.4f,%8.4f)=R%d;\n",
                        rx[(i - 1) * 2 + 1] + d,rx[(i - 1) * 2 + 2],i);
    }
    for (i = 1; i <= nc; ++i) {
        fprintf(PMPCFd,"pltext(fs=2,xy=%8.4f,%8.4f)=C%d;\n",
                        cx[(i - 1) * 2 + 1] + d,cx[(i - 1) * 2 + 2],i);
    }

    printf1("Plot command file written to: %s\n",PMPCFName);
}

