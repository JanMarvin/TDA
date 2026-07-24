/****************************************************************************/
/*  t_seqe                                                                  */
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
#include "t_gdat.h"
#include "t_var.h"
#include "t_gf.h"
#include "t_parm.h"
#include "t_alloc.h"
#include "t_seq.h"
#include "tda_context.h"

/*  functions in t_seqe.c */

int seqev(TDAContext *ctx);
int seqevd(TDAContext *ctx);
int get_ev(TDAContext *ctx, int sn,int ns,int *nev);
int seqmd(TDAContext *ctx);
int get_nxe(TDAContext *ctx, int n,int *na);

/*--------------------------------------------------------------------------*/
/*  seqev()         Info about events.                                      */
/*                  seqev(sn=,sel=) = fname (fname optional)                */
/*                  Return 0 if OK, -1 if error.                            */

int seqev(TDAContext *ctx)
{
    register int i,j;
    int err,n,nn,ns,sn;

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);
         
    printf1(ctx, "Counting events.\n");
         
    if (parm(ctx, ctx->CmdBuf + 5,1,0))    /* get parameters */
        goto SEQEVFin;
        
    sn = seq_getsn(ctx, ctx->PMSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQEVFin;

    ns = ctx->SeqSTN[sn];            /* number of states */

    if (alloc_acn(ctx, ns * ns + 1))
        goto SEQEVFin;

    nn = get_ev(ctx, sn,ns,ctx->AcN);
    if (nn < 1) {
        printf1(ctx, "Error: number of cases is zero.\n");
        goto SEQEVFin;
    }
    if (ctx->PMFDef == 0) {
        printf1(ctx, "\nEvent type    Number\n");
        prnchar(ctx, '-',20,1);
    }
    nn = 0;
    for (i = 0; i < ns; ++i) {
        for (j = 0; j < ns; ++j) {
            if (j != i) {
                n = ctx->AcN[i * ns + j];
                if (n > 0) {
                    if (ctx->PMFDef)
                        fprintf(ctx->PMFd,"%4d %4d %10d\n",ctx->SeqSTNI[sn][i],ctx->SeqSTNI[sn][j],n);
                    else
                        printf1(ctx, "%4d %4d %10d\n",ctx->SeqSTNI[sn][i],ctx->SeqSTNI[sn][j],n);
#ifdef TDA_R_PACKAGE
                    {
                        double erow[3];
                        erow[0] = (double)ctx->SeqSTNI[sn][i];
                        erow[1] = (double)ctx->SeqSTNI[sn][j];
                        erow[2] = (double)n;
                        tda_export_row(ctx, "seqev.table", erow, 3);
                    }
#endif
                    nn++;
                }
            }
        }
    }   
    if (ctx->PMFDef) 
        printf1(ctx, "%d records written to: %s\n",nn,ctx->PMFdName);

    err = 0;

SEQEVFin:
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "seqev.table");
#endif
    p_clean(ctx);
    return(err);
}

/*--##----------------------------------------------------------------------*/
/*  seqevd()        data file with all events and time points.              */
/*                  seqevd(sn=,sel=) = fname                                */
/*                  Return 0 if OK, -1 if error.                            */

int seqevd(TDAContext *ctx)
{
    register int i,j,k;
    int err,n,n1,ns,nss,sn,a,b,t,nn,nrec;

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);
         
    printf1(ctx, "Counting events separately for each time point.\n");
         
    if (parm(ctx, ctx->CmdBuf + 6,1,1))    /* get parameters */
        goto SEQEVDFin;
        
    sn = seq_getsn(ctx, ctx->PMSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQEVDFin;

    ns = ctx->SeqSTN[sn];            /* number of states */
    nss = ns * ns;

    if (alloc_acn(ctx, nss + 1))
        goto SEQEVDFin;

    if (alloc_acm(ctx, nss + 1))
        goto SEQEVDFin;

    nn = get_ev(ctx, sn,ns,ctx->AcN);
    if (nn < 1) {
        printf1(ctx, "Error: number of cases is zero.\n");
        goto SEQEVDFin;
    }
    a = ctx->SeqTMin[sn];
    b = ctx->SeqTMax[sn];

    nrec = 0;
    for (t = a; t < b; ++t) {
   
        for (i = 0; i < nss; ++i)
            ctx->AcM[i] = 0;         
    
        nn = 0;
        for (i = 0; i < ctx->NOC; ++i) {

            if (eval_sve(ctx, i) == 0)   
                continue;

            n = seq_sget(ctx, i,t,sn);
            if (n >= 0) {
                n = ctx->SeqSTNII[sn][n];
                n1 = seq_sget(ctx, i,t + 1,sn);
                if (n1 >= 0) {
                    n1 = ctx->SeqSTNII[sn][n1];
                    ctx->AcM[n * ns + n1] += 1;
                    nn++;
                }
            }   
        }
        fprintf(ctx->PMFd,"%6d %6d ",t + 1,nn);
#ifdef TDA_R_PACKAGE
        tda_export_cell(ctx, "seqevd.table", (double)(t + 1));
        tda_export_cell(ctx, "seqevd.table", (double)nn);
#endif
        nn = 0;
        for (i = 0; i < ns; ++i) {
            for (j = 0; j < ns; ++j) {
                if (j != i && ctx->AcN[i * ns + j] > 0) {
                    fprintf(ctx->PMFd,"%4d ",ctx->AcM[i * ns + j]);
#ifdef TDA_R_PACKAGE
                    tda_export_cell(ctx, "seqevd.table",
                                    (double)ctx->AcM[i * ns + j]);
#endif
                    nn += ctx->AcM[i * ns + j];
                }
            }
        }
        fprintf(ctx->PMFd,"%6d\n",nn);
#ifdef TDA_R_PACKAGE
        tda_export_cell(ctx, "seqevd.table", (double)nn);
        tda_export_endrow(ctx, "seqevd.table");
#endif
        nrec++;
    }   
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);

    if (ctx->PMTDAFDef && nrec > 0) {                    /* write TDA description */

        fprintf(ctx->PMTDAFd,"nvar(\n");
        fprintf(ctx->PMTDAFd,"  dfile = %s,\n",ctx->PMFdName);
        fprintf(ctx->PMTDAFd,"  noc = %d,\n",nrec);

        fprintf(ctx->PMTDAFd,"  TIME <5>[6.0] = c1 , # time\n");
        fprintf(ctx->PMTDAFd,"  NCAS <5>[6.0] = c2 , # number of cases\n");
           
        k = 3;
        for (i = 0; i < ns; ++i) {
            for (j = 0; j < ns; ++j) {
                if (j != i && ctx->AcN[i * ns + j] > 0) {
                    n  = ctx->SeqSTNI[sn][i];
                    n1 = ctx->SeqSTNI[sn][j];
                    fprintf(ctx->PMTDAFd,"  EV%d_%d <5>[4.0] = c%d , # number of events (%d,%d)\n",  
                                            n,n1,k++,n,n1);
                }
            }
        }
        fprintf(ctx->PMTDAFd,"  NEV <5>[6.0] = c%d , # total number of events\n);\n",k);
        printf1(ctx, "TDA description written to: %s\n",ctx->PMTDAFName);
    }
    err = 0;

SEQEVDFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  get_ev(sn,ns,nve)   Create matrix nev[] with number of events for       */
/*                      sequence data structure sn. ns is number of states. */
/*                      Return number of cases.                             */

int get_ev(TDAContext *ctx, int sn,int ns,int *nev)
{
    register int i,t;
    int a,b,n,n1,nn;

    a = ctx->SeqTMin[sn];
    b = ctx->SeqTMax[sn];
    nn = 0;
    for (t = a; t < b; ++t) {
       
        for (i = 0; i < ctx->NOC; ++i) {

            if (eval_sve(ctx, i) == 0)   
                continue;

            n = seq_sget(ctx, i,t,sn);
            if (n >= 0) {
                n = ctx->SeqSTNII[sn][n];
                n1 = seq_sget(ctx, i,t + 1,sn);
                if (n1 >= 0) {
                    n1 = ctx->SeqSTNII[sn][n1];
                    nev[n * ns + n1] += 1;
                    nn++;
                }
            }   
        }
    }
    return(nn);
}

/*--------------------------------------------------------------------------*/
/*  seqmd()     Creating data for regression models for events.             */
/*              Parameter:                                                  */
/*              sn =        number of sequence data structure.              */
/*              ev = [j,k]  selected event.                                 */
/*                  Return 0 if OK, -1 if error.                            */

int seqmd(TDAContext *ctx)
{
    register int i,j,k,l;
    int err,n,n1,n2,ns,nss,sn,ta,tb,tua,tub,t,t1,nn,nnu,nt;
    int ej,ek,ne,nxe,kk,jj,s1,s2;

    err = -1;         
    if (check_cmd(ctx, 0))
        return(-1);
         
    printf1(ctx, "Creating data for model estimation. Current memory: %d bytes.\n",ctx->MemReq);
         
    if (parm(ctx, ctx->CmdBuf + 5,1,0))    /* get parameters */
        goto SEQMDFin;
        
    if (ctx->PMEVSN < 1 || ctx->PMEV1 < 0 || ctx->PMEV2 < 0) {
        printf1(ctx, "Error: need the definition of an event.\n");
        goto SEQMDFin;
    }
    printf1(ctx, "Selected target event: [%d,%d,%d]\n",ctx->PMEVSN,ctx->PMEV1,ctx->PMEV2);

    sn = seq_getsn(ctx, ctx->PMEVSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQMDFin;

    ns = ctx->SeqSTN[sn];            /* number of states */
    nss = ns * ns;

    if (ctx->PMEV1 == ctx->PMEV2) {
        printf1(ctx, "Error: not an event.\n");
        goto SEQMDFin;
    }
    j = ctx->SeqSTNH[sn] - 1;
                       
    if (ctx->PMEV1 > j || ctx->PMEV2 > j) {
        printf1(ctx, "Error: maximal state number is %d.\n",j);
        goto SEQMDFin;
    }
    ej = ctx->SeqSTNII[sn][ctx->PMEV1];
    ek = ctx->SeqSTNII[sn][ctx->PMEV2];
    if (ej < 0 || ek < 0) {
        printf1(ctx, "Error: event does not occur in the selected sequence data structure.\n");
        goto SEQMDFin;
    }
    if (alloc_acn(ctx, nss + 1))
        goto SEQMDFin;

    nn = get_ev(ctx, sn,ns,ctx->AcN);         /* get transition matrix */
    if (nn < 1) {
        printf1(ctx, "Error: number of selected cases is zero.\n");
        goto SEQMDFin;
    }
    ne = ctx->AcN[ej * ns + ek];
    if (ne < 1) {
        printf1(ctx, "Error: number of events is zero.\n");
        goto SEQMDFin;
    }
    ta = ctx->SeqTMin[sn];
    tb = ctx->SeqTMax[sn];

    if (ctx->PMNTP > 0) {
        tua = (int)ctx->PMTP[0];
        tub = (int)ctx->PMTP[ctx->PMNTP - 1];
        nt = ctx->PMNTP;
    }
    else {                          /* default time axis */
        tua = ta + 1;
        tub = tb;
        nt = tub - tua + 1;
    }
    printf1(ctx, "Selected time axis for events: %d -- %d\n",tua,tub);
    if (tua <= ta || tub > tb || tub < tua) {
        printf1(ctx, "Error: cannot use this time axis. Max range is: %d -- %d\n",ta + 1,tb);
        goto SEQMDFin;
    }

    /* create risk set and number of events   
       AcS[] = time points
       AcR[] = risk set
       AcN[] = number of events */

    if (alloc_acs(ctx, nt))
        goto SEQMDFin;

    if (alloc_acr(ctx, nt))
        goto SEQMDFin;

    if (alloc_acn(ctx, nt))
        goto SEQMDFin;

    for (k = 0; k < nt; ++k) {
        if (ctx->PMNTP > 0)  
            ctx->AcS[k] = (int)ctx->PMTP[k];
        else
            ctx->AcS[k] = tua + k;
    }

    nn = 0;     /* records in data matrix */
    nnu = 0;    /* units in data matrix */

    for (i = 0; i < ctx->NOC; ++i) {

        if (eval_sve(ctx, i) == 0)   
            continue;

        j = 0;
        for (k = 0; k < nt; ++k) {

            t = ctx->AcS[k];
            n = seq_sget(ctx, i,t - 1,sn);
            if (n == ctx->PMEV1) {
                n1 = seq_sget(ctx, i,t,sn);
                if (n1 >= 0) {
                    nn++;
                    j++;
                    ctx->AcR[k] += 1;
                    if (n1 == ctx->PMEV2)
                        ctx->AcN[k] += 1;
                }
            }   
        }
        if (j)
            nnu++;
    }   

    printf1(ctx, "\nTime  Risk Set  Events\n");
    prnchar(ctx, '-',22,1);
             
    for (k = 0; k < nt; ++k) { 
        printf1(ctx, "%4d  %8d %7d\n",ctx->AcS[k],ctx->AcR[k],ctx->AcN[k]);
#ifdef TDA_R_PACKAGE
        {
            double erow[3];
            erow[0] = (double)ctx->AcS[k];
            erow[1] = (double)ctx->AcR[k];
            erow[2] = (double)ctx->AcN[k];
            tda_export_row(ctx, "seqmd.risk", erow, 3);
        }
#endif
    }
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "seqmd.risk");
#endif
       
    if (ctx->PMFDef == 0) {
        err = 0;
        goto SEQMDFin;
    }
    printf1(ctx, "\nCreating data matrix: %d records, %d units.\n",nn,nnu);
    printf1(ctx, "Number of time-independent covariates: %d\n",ctx->PMNV);
    nxe = ctx->XEFlg;
    printf1(ctx, "Number of event variables: %d\n",nxe);

    /* first check of event variables */

    if (nxe > 0) {
        if (alloc_ack(ctx, nxe * 4))
            goto SEQMDFin;

        if (get_nxe(ctx, nxe,ctx->AcK))
            goto SEQMDFin;

        if (alloc_acm(ctx, nxe))
            goto SEQMDFin;

        if (alloc_aci(ctx, nxe))
            goto SEQMDFin;

        if (alloc_acj(ctx, nxe))
            goto SEQMDFin;

        newline(ctx);
    }
    nn = 0;
    for (i = 0; i < ctx->NOC; ++i) {

        if (eval_sve(ctx, i) == 0)   
            continue;

        if (nxe > 0) {      /* initialize count of events */

            for (l = 0; l < nxe; ++l) {
                kk = ctx->AcK[l * 4];         
                if (kk >= 0) {
                    ctx->AcM[l] = 0;
                    ctx->AcI[l] = ctx->SeqTMin[kk] + 1;
                    ctx->AcJ[l] = -1;        
                }
            }
        }

        for (k = 0; k < nt; ++k) {

            t = ctx->AcS[k];
            n = seq_sget(ctx, i,t - 1,sn);
            if (n == ctx->PMEV1) {
                n1 = seq_sget(ctx, i,t,sn);
                if (n1 >= 0) {
                    if (n1 == ctx->PMEV2)
                        j = 1;
                    else
                        j = 0;

                    /* basic variables: id, time, indicator for event */

                    fprintf(ctx->PMFd,"%6d %4d %d ",i + 1,t,j);
#ifdef TDA_R_PACKAGE
                    tda_export_cell(ctx, "seqmd.table", (double)(i + 1));
                    tda_export_cell(ctx, "seqmd.table", (double)t);
                    tda_export_cell(ctx, "seqmd.table", (double)j);
#endif

                    /* add period dummies */

                    for (l = 0; l < nt; ++l) {
                        t1 = ctx->AcS[l];
                        j = 0;
                        if (t1 == t)
                            j = 1;
                        fprintf(ctx->PMFd,"%d ",j);
#ifdef TDA_R_PACKAGE
                        tda_export_cell(ctx, "seqmd.table", (double)j);
#endif
                    }

                    /* add time-independent variables */

                    if (ctx->PMNV > 0) {
                        for (l = 0; l < ctx->PMNV; ++l) {
                            kk = ctx->PMVIdx[l];
                            rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[kk],get_data(ctx, kk,i));
#ifdef TDA_R_PACKAGE
                            tda_export_cell(ctx, "seqmd.table",
                                            get_data(ctx, kk,i));
#endif
                        }
                    }

                    /* add event-specific variables */

                    if (nxe > 0) {
                        for (l = 0; l < nxe; ++l) {
                            kk = ctx->AcK[l * 4];         
                            if (kk < 0) {
                                t1 = (int)get_data(ctx, ctx->AcK[l * 4 + 1],i);
                                fprintf(ctx->PMFd,"%4d ",t - t1);
#ifdef TDA_R_PACKAGE
                                tda_export_cell(ctx, "seqmd.table",
                                                (double)(t - t1));
#endif
                            }
                            else {

                                /* first update number of events */

                                j  = ctx->AcI[l];
                                jj = imin(ctx, t - 1,ctx->SeqTMax[kk]);

                                if (j <= jj) {
                                    s1 = ctx->AcK[l * 4 + 1];
                                    s2 = ctx->AcK[l * 4 + 2];
                                    n1 = seq_sget(ctx, i,j - 1,kk);

                                    while (j <= jj) {
                                        n2 = seq_sget(ctx, i,j,kk);
                       
                                        if (n1 == s1 && n2 == s2) {
                                            ctx->AcM[l] += 1;
                                            ctx->AcJ[l] = j;
                                        }
                                        n1 = n2;
                                        j++;
                                    }
                                    ctx->AcI[l] = jj + 1;
                                }
                                fprintf(ctx->PMFd,"%4d %4d ",ctx->AcM[l],ctx->AcJ[l]);
#ifdef TDA_R_PACKAGE
                                tda_export_cell(ctx, "seqmd.table",
                                                (double)ctx->AcM[l]);
                                tda_export_cell(ctx, "seqmd.table",
                                                (double)ctx->AcJ[l]);
#endif

                                /* add effect shape variables */

                                n1 = ctx->AcJ[l];
                                for (j = 0; j <= ctx->AcK[l * 4 + 3]; ++j) {
                                    n2 = 0;
                                    if (n1 >= 0 && n1 + j == t)
                                        n2 = 1;
                                    fprintf(ctx->PMFd,"%d ",n2);
#ifdef TDA_R_PACKAGE
                                    tda_export_cell(ctx, "seqmd.table",
                                                    (double)n2);
#endif
                                }
                            }
                        }
                    }
                    fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
                    tda_export_endrow(ctx, "seqmd.table");
#endif
                    nn++;
                }
            }   
        }
    }   
    printf1(ctx, "%d records written to: %s\n",nn,ctx->PMFdName);

    if (ctx->PMTDAFDef && nn > 0) {                    /* write TDA description */

        fprintf(ctx->PMTDAFd,"nvar(\n");
        fprintf(ctx->PMTDAFd,"  dfile = %s,\n",ctx->PMFdName);
        fprintf(ctx->PMTDAFd,"  noc = %d,\n",nn);

        fprintf(ctx->PMTDAFd,"  ID   <5>[6.0] = c1 , # ID\n");
        fprintf(ctx->PMTDAFd,"  TIME <5>[4.0] = c2 , # time\n");
        fprintf(ctx->PMTDAFd,"  Z    <1>[1.0] = c3 , # target event\n");
        k = 4;
        for (l = 0; l < nt; ++l) {
            t = ctx->AcS[l];
            fprintf(ctx->PMTDAFd,"  P%-4d<1>[1.0] = c%d , # period %d\n",t,k++,t);
        }
        for (l = 0; l < ctx->PMNV; ++l) {
            j = ctx->PMVIdx[l];
            fprintf(ctx->PMTDAFd,"  %s <%d>[%d.%d]",
                                    ctx->VName[j],ctx->VSLen[j],ctx->VPFmt1[j],ctx->VPFmt2[j]);
            if (ctx->VLabel[j] != NULL)                                                      
                fprintf(ctx->PMTDAFd,"(%s)",ctx->VLabel[j]);                                           
            fprintf(ctx->PMTDAFd," = c%d,\n",k++);
        }
        if (nxe > 0) {
            for (l = 0; l < nxe; ++l) {
                kk = ctx->AcK[l * 4];         
                if (kk < 0) {
                    j = ctx->AcK[l * 4 + 1];
                    fprintf(ctx->PMTDAFd,"  %s_D <5>[4.0] = c%d,\n",ctx->VName[j],k++);
                }
                else {
                    s1 = ctx->AcK[l * 4 + 1];
                    s2 = ctx->AcK[l * 4 + 2];

                    fprintf(ctx->PMTDAFd,"  E_%d_%d  <2>[4.0] = c%d,\n",s1,s2,k++);
                    fprintf(ctx->PMTDAFd,"  ED_%d_%d <2>[4.0] = c%d,\n",s1,s2,k++);

                    for (j = 0; j <= ctx->AcK[l * 4 + 3]; ++j)  
                        fprintf(ctx->PMTDAFd,"  E%d_%d_%d <2>[4.0] = c%d,\n",j,s1,s2,k++);
                }
            }
        }
        fprintf(ctx->PMTDAFd,");\n");
        printf1(ctx, "TDA description written to: %s\n",ctx->PMTDAFName);
    }
    err = 0;

SEQMDFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  get_nxe(n,na)   Get parameters for event variables into na[]            */
/*                  na[i + 4 + 0] = sequence number, or -1 if a standard    */
/*                                  variable                                */  
/*                  na[i + 4 + 1] = origin state, or number of standard     */
/*                                  variable                                */
/*                  na[i + 4 + 2] = destination state                       */
/*                  na[i + 4 + 3] = dur                                     */
/*                  i = 0,...,n - 1.                                        */
/*                  String to be scanned is in PMRHSTR                      */
/*                  Return 0 if OK, -1 if error.                            */

int get_nxe(TDAContext *ctx, int n,int *na)
{
    register int i,l;
    int err,sn,j,k,d;
    register char *p,*q;
    char vname[VNLMax + 1];

    err = -1;
    p = ctx->PMRHSTR;
         printf1(ctx, "p:%s:\n",p);      

    printf1(ctx, "\nDefinition       Variables\n");
    prnchar(ctx, '-',26,1);
    for (i = 0; i < n; ++i) {
   
        if (*p != '[') {
         printf1(ctx, "1p:%s:\n",p);      
            err = -2;
            goto GETNXEFin;
        }
        q = p + 1;
        while (*q && *q != ']')
            q++;
        if (*q++ != ']') {
         printf1(ctx, "2p:%s:\n",p);      
            err = -2;
            goto GETNXEFin;
        }
        *q = '\0';
        printf1(ctx, "%s",p);
        prnchar(ctx, ' ',(int)(17 - strlen(p)),0);

        if ((j = get_vnlen(ctx, p + 1)) > 1) {
            k = get_vidx1(ctx, p + 1,vname);
            if (k < 0) {
                printf1(ctx, "error: variable does not exist.\n");
                goto GETNXEFin;
            }
            /*************
            if (j > 14) {
                printf1(ctx, "error: exceeded max length of variable name.\n");
                goto GETNXEFin;
            }
            *******************/
            na[i * 4] = -1;
            na[i * 4 + 1] = k;                 

            printf1(ctx, "%s_D\n",vname);
        }
        else if (sscanf(p,"[%d,%d,%d,%d]",&sn,&j,&k,&d) == 4) {
            sn--;
            if (sn < 0 || sn > ctx->SeqDNH || ctx->SeqDT[sn] == 0) {
                printf1(ctx, "error: sequence %d not defined.\n",sn + 1);
                goto GETNXEFin;
            }
            if (j < 0 || k < 0 || j >= ctx->SeqSTNH[sn] || k >= ctx->SeqSTNH[sn]) {
                printf1(ctx, "error: event [%d,%d] not defined.\n",j,k);
                goto GETNXEFin;
            }
            na[i * 4] = sn;
            na[i * 4 + 1] = j;
            na[i * 4 + 2] = k;
            na[i * 4 + 3] = d;                 

            printf1(ctx, "E_%d_%d  ",j,k);
            printf1(ctx, "ED_%d_%d  ",j,k);
            for (l = 0; l <= d; ++l)  
                printf1(ctx, "E%d_%d_%d  ",l,j,k);
            printf1(ctx, "\n");
        }
        else {
            err = -2;
            goto GETNXEFin;
        }
        p = q + 1;
    }
    err = 0;

GETNXEFin:
    if (err == -2)
        printf1(ctx, "Syntax error in definition of event variables.\n");
    return(err);
}

