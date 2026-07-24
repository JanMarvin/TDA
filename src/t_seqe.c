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

/*  functions in t_seqe.c */

int seqev(void);
int seqevd(void);
int get_ev(int sn,int ns,int *nev);
int seqmd(void);
int get_nxe(int n,int *na);

/*--------------------------------------------------------------------------*/
/*  seqev()         Info about events.                                      */
/*                  seqev(sn=,sel=) = fname (fname optional)                */
/*                  Return 0 if OK, -1 if error.                            */

int seqev(void)
{
    register int i,j;
    int err,n,nn,ns,sn;

    err = -1;         
    if (check_cmd(1))
        return(-1);
         
    printf1("Counting events.\n");
         
    if (parm(CmdBuf + 5,1,0))    /* get parameters */
        goto SEQEVFin;
        
    sn = seq_getsn(PMSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQEVFin;

    ns = SeqSTN[sn];            /* number of states */

    if (alloc_acn(ns * ns + 1))
        goto SEQEVFin;

    nn = get_ev(sn,ns,AcN);
    if (nn < 1) {
        printf1("Error: number of cases is zero.\n");
        goto SEQEVFin;
    }
    if (PMFDef == 0) {
        printf1("\nEvent type    Number\n");
        prnchar('-',20,1);
    }
    nn = 0;
    for (i = 0; i < ns; ++i) {
        for (j = 0; j < ns; ++j) {
            if (j != i) {
                n = AcN[i * ns + j];
                if (n > 0) {
                    if (PMFDef)
                        fprintf(PMFd,"%4d %4d %10d\n",SeqSTNI[sn][i],SeqSTNI[sn][j],n);
                    else
                        printf1("%4d %4d %10d\n",SeqSTNI[sn][i],SeqSTNI[sn][j],n);
                    nn++;
                }
            }
        }
    }   
    if (PMFDef) 
        printf1("%d records written to: %s\n",nn,PMFdName);

    err = 0;

SEQEVFin:
    p_clean();
    return(err);
}

/*--##----------------------------------------------------------------------*/
/*  seqevd()        data file with all events and time points.              */
/*                  seqevd(sn=,sel=) = fname                                */
/*                  Return 0 if OK, -1 if error.                            */

int seqevd(void)
{
    register int i,j,k;
    int err,n,n1,ns,nss,sn,a,b,t,nn,nrec;

    err = -1;         
    if (check_cmd(1))
        return(-1);
         
    printf1("Counting events separately for each time point.\n");
         
    if (parm(CmdBuf + 6,1,1))    /* get parameters */
        goto SEQEVDFin;
        
    sn = seq_getsn(PMSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQEVDFin;

    ns = SeqSTN[sn];            /* number of states */
    nss = ns * ns;

    if (alloc_acn(nss + 1))
        goto SEQEVDFin;

    if (alloc_acm(nss + 1))
        goto SEQEVDFin;

    nn = get_ev(sn,ns,AcN);
    if (nn < 1) {
        printf1("Error: number of cases is zero.\n");
        goto SEQEVDFin;
    }
    a = SeqTMin[sn];
    b = SeqTMax[sn];

    nrec = 0;
    for (t = a; t < b; ++t) {
   
        for (i = 0; i < nss; ++i)
            AcM[i] = 0;         
    
        nn = 0;
        for (i = 0; i < NOC; ++i) {

            if (eval_sve(i) == 0)   
                continue;

            n = seq_sget(i,t,sn);
            if (n >= 0) {
                n = SeqSTNII[sn][n];
                n1 = seq_sget(i,t + 1,sn);
                if (n1 >= 0) {
                    n1 = SeqSTNII[sn][n1];
                    AcM[n * ns + n1] += 1;
                    nn++;
                }
            }   
        }
        fprintf(PMFd,"%6d %6d ",t + 1,nn);
        nn = 0;
        for (i = 0; i < ns; ++i) {
            for (j = 0; j < ns; ++j) {
                if (j != i && AcN[i * ns + j] > 0) {
                    fprintf(PMFd,"%4d ",AcM[i * ns + j]);
                    nn += AcM[i * ns + j];
                }
            }
        }
        fprintf(PMFd,"%6d\n",nn);
        nrec++;
    }   
    printf1("%d records written to: %s\n",nrec,PMFdName);

    if (PMTDAFDef && nrec > 0) {                    /* write TDA description */

        fprintf(PMTDAFd,"nvar(\n");
        fprintf(PMTDAFd,"  dfile = %s,\n",PMFdName);
        fprintf(PMTDAFd,"  noc = %d,\n",nrec);

        fprintf(PMTDAFd,"  TIME <5>[6.0] = c1 , # time\n");
        fprintf(PMTDAFd,"  NCAS <5>[6.0] = c2 , # number of cases\n");
           
        k = 3;
        for (i = 0; i < ns; ++i) {
            for (j = 0; j < ns; ++j) {
                if (j != i && AcN[i * ns + j] > 0) {
                    n  = SeqSTNI[sn][i];
                    n1 = SeqSTNI[sn][j];
                    fprintf(PMTDAFd,"  EV%d_%d <5>[4.0] = c%d , # number of events (%d,%d)\n",  
                                            n,n1,k++,n,n1);
                }
            }
        }
        fprintf(PMTDAFd,"  NEV <5>[6.0] = c%d , # total number of events\n);\n",k);
        printf1("TDA description written to: %s\n",PMTDAFName);
    }
    err = 0;

SEQEVDFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  get_ev(sn,ns,nve)   Create matrix nev[] with number of events for       */
/*                      sequence data structure sn. ns is number of states. */
/*                      Return number of cases.                             */

int get_ev(int sn,int ns,int *nev)
{
    register int i,t;
    int a,b,n,n1,nn;

    a = SeqTMin[sn];
    b = SeqTMax[sn];
    nn = 0;
    for (t = a; t < b; ++t) {
       
        for (i = 0; i < NOC; ++i) {

            if (eval_sve(i) == 0)   
                continue;

            n = seq_sget(i,t,sn);
            if (n >= 0) {
                n = SeqSTNII[sn][n];
                n1 = seq_sget(i,t + 1,sn);
                if (n1 >= 0) {
                    n1 = SeqSTNII[sn][n1];
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

int seqmd(void)
{
    register int i,j,k,l;
    int err,n,n1,n2,ns,nss,sn,ta,tb,tua,tub,t,t1,nn,nnu,nt;
    int ej,ek,ne,nxe,kk,jj,s1,s2;

    err = -1;         
    if (check_cmd(0))
        return(-1);
         
    printf1("Creating data for model estimation. Current memory: %d bytes.\n",MemReq);
         
    if (parm(CmdBuf + 5,1,0))    /* get parameters */
        goto SEQMDFin;
        
    if (PMEVSN < 1 || PMEV1 < 0 || PMEV2 < 0) {
        printf1("Error: need the definition of an event.\n");
        goto SEQMDFin;
    }
    printf1("Selected target event: [%d,%d,%d]\n",PMEVSN,PMEV1,PMEV2);

    sn = seq_getsn(PMEVSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQMDFin;

    ns = SeqSTN[sn];            /* number of states */
    nss = ns * ns;

    if (PMEV1 == PMEV2) {
        printf1("Error: not an event.\n");
        goto SEQMDFin;
    }
    j = SeqSTNH[sn] - 1;
                       
    if (PMEV1 > j || PMEV2 > j) {
        printf1("Error: maximal state number is %d.\n",j);
        goto SEQMDFin;
    }
    ej = SeqSTNII[sn][PMEV1];
    ek = SeqSTNII[sn][PMEV2];
    if (ej < 0 || ek < 0) {
        printf1("Error: event does not occur in the selected sequence data structure.\n");
        goto SEQMDFin;
    }
    if (alloc_acn(nss + 1))
        goto SEQMDFin;

    nn = get_ev(sn,ns,AcN);         /* get transition matrix */
    if (nn < 1) {
        printf1("Error: number of selected cases is zero.\n");
        goto SEQMDFin;
    }
    ne = AcN[ej * ns + ek];
    if (ne < 1) {
        printf1("Error: number of events is zero.\n");
        goto SEQMDFin;
    }
    ta = SeqTMin[sn];
    tb = SeqTMax[sn];

    if (PMNTP > 0) {
        tua = (int)PMTP[0];
        tub = (int)PMTP[PMNTP - 1];
        nt = PMNTP;
    }
    else {                          /* default time axis */
        tua = ta + 1;
        tub = tb;
        nt = tub - tua + 1;
    }
    printf1("Selected time axis for events: %d -- %d\n",tua,tub);
    if (tua <= ta || tub > tb || tub < tua) {
        printf1("Error: cannot use this time axis. Max range is: %d -- %d\n",ta + 1,tb);
        goto SEQMDFin;
    }

    /* create risk set and number of events   
       AcS[] = time points
       AcR[] = risk set
       AcN[] = number of events */

    if (alloc_acs(nt))
        goto SEQMDFin;

    if (alloc_acr(nt))
        goto SEQMDFin;

    if (alloc_acn(nt))
        goto SEQMDFin;

    for (k = 0; k < nt; ++k) {
        if (PMNTP > 0)  
            AcS[k] = (int)PMTP[k];
        else
            AcS[k] = tua + k;
    }

    nn = 0;     /* records in data matrix */
    nnu = 0;    /* units in data matrix */

    for (i = 0; i < NOC; ++i) {

        if (eval_sve(i) == 0)   
            continue;

        j = 0;
        for (k = 0; k < nt; ++k) {

            t = AcS[k];
            n = seq_sget(i,t - 1,sn);
            if (n == PMEV1) {
                n1 = seq_sget(i,t,sn);
                if (n1 >= 0) {
                    nn++;
                    j++;
                    AcR[k] += 1;
                    if (n1 == PMEV2)
                        AcN[k] += 1;
                }
            }   
        }
        if (j)
            nnu++;
    }   

    printf1("\nTime  Risk Set  Events\n");
    prnchar('-',22,1);
             
    for (k = 0; k < nt; ++k)  
        printf1("%4d  %8d %7d\n",AcS[k],AcR[k],AcN[k]);
       
    if (PMFDef == 0) {
        err = 0;
        goto SEQMDFin;
    }
    printf1("\nCreating data matrix: %d records, %d units.\n",nn,nnu);
    printf1("Number of time-independent covariates: %d\n",PMNV);
    nxe = XEFlg;
    printf1("Number of event variables: %d\n",nxe);

    /* first check of event variables */

    if (nxe > 0) {
        if (alloc_ack(nxe * 4))
            goto SEQMDFin;

        if (get_nxe(nxe,AcK))
            goto SEQMDFin;

        if (alloc_acm(nxe))
            goto SEQMDFin;

        if (alloc_aci(nxe))
            goto SEQMDFin;

        if (alloc_acj(nxe))
            goto SEQMDFin;

        newline();
    }
    nn = 0;
    for (i = 0; i < NOC; ++i) {

        if (eval_sve(i) == 0)   
            continue;

        if (nxe > 0) {      /* initialize count of events */

            for (l = 0; l < nxe; ++l) {
                kk = AcK[l * 4];         
                if (kk >= 0) {
                    AcM[l] = 0;
                    AcI[l] = SeqTMin[kk] + 1;
                    AcJ[l] = -1;        
                }
            }
        }

        for (k = 0; k < nt; ++k) {

            t = AcS[k];
            n = seq_sget(i,t - 1,sn);
            if (n == PMEV1) {
                n1 = seq_sget(i,t,sn);
                if (n1 >= 0) {
                    if (n1 == PMEV2)
                        j = 1;
                    else
                        j = 0;

                    /* basic variables: id, time, indicator for event */

                    fprintf(PMFd,"%6d %4d %d ",i + 1,t,j);

                    /* add period dummies */

                    for (l = 0; l < nt; ++l) {
                        t1 = AcS[l];
                        j = 0;
                        if (t1 == t)
                            j = 1;
                        fprintf(PMFd,"%d ",j);
                    }

                    /* add time-independent variables */

                    if (PMNV > 0) {
                        for (l = 0; l < PMNV; ++l) {
                            kk = PMVIdx[l];
                            fprintf(PMFd,VPFmtS[kk],get_data(kk,i));
                        }
                    }

                    /* add event-specific variables */

                    if (nxe > 0) {
                        for (l = 0; l < nxe; ++l) {
                            kk = AcK[l * 4];         
                            if (kk < 0) {
                                t1 = (int)get_data(AcK[l * 4 + 1],i);
                                fprintf(PMFd,"%4d ",t - t1);
                            }
                            else {

                                /* first update number of events */

                                j  = AcI[l];
                                jj = imin(t - 1,SeqTMax[kk]);

                                if (j <= jj) {
                                    s1 = AcK[l * 4 + 1];
                                    s2 = AcK[l * 4 + 2];
                                    n1 = seq_sget(i,j - 1,kk);

                                    while (j <= jj) {
                                        n2 = seq_sget(i,j,kk);
                       
                                        if (n1 == s1 && n2 == s2) {
                                            AcM[l] += 1;
                                            AcJ[l] = j;
                                        }
                                        n1 = n2;
                                        j++;
                                    }
                                    AcI[l] = jj + 1;
                                }
                                fprintf(PMFd,"%4d %4d ",AcM[l],AcJ[l]);

                                /* add effect shape variables */

                                n1 = AcJ[l];
                                for (j = 0; j <= AcK[l * 4 + 3]; ++j) {
                                    n2 = 0;
                                    if (n1 >= 0 && n1 + j == t)
                                        n2 = 1;
                                    fprintf(PMFd,"%d ",n2);
                                }
                            }
                        }
                    }
                    fprintf(PMFd,"\n");
                    nn++;
                }
            }   
        }
    }   
    printf1("%d records written to: %s\n",nn,PMFdName);

    if (PMTDAFDef && nn > 0) {                    /* write TDA description */

        fprintf(PMTDAFd,"nvar(\n");
        fprintf(PMTDAFd,"  dfile = %s,\n",PMFdName);
        fprintf(PMTDAFd,"  noc = %d,\n",nn);

        fprintf(PMTDAFd,"  ID   <5>[6.0] = c1 , # ID\n");
        fprintf(PMTDAFd,"  TIME <5>[4.0] = c2 , # time\n");
        fprintf(PMTDAFd,"  Z    <1>[1.0] = c3 , # target event\n");
        k = 4;
        for (l = 0; l < nt; ++l) {
            t = AcS[l];
            fprintf(PMTDAFd,"  P%-4d<1>[1.0] = c%d , # period %d\n",t,k++,t);
        }
        for (l = 0; l < PMNV; ++l) {
            j = PMVIdx[l];
            fprintf(PMTDAFd,"  %s <%d>[%d.%d]",
                                    VName[j],VSLen[j],VPFmt1[j],VPFmt2[j]);
            if (VLabel[j] != NULL)                                                      
                fprintf(PMTDAFd,"(%s)",VLabel[j]);                                           
            fprintf(PMTDAFd," = c%d,\n",k++);
        }
        if (nxe > 0) {
            for (l = 0; l < nxe; ++l) {
                kk = AcK[l * 4];         
                if (kk < 0) {
                    j = AcK[l * 4 + 1];
                    fprintf(PMTDAFd,"  %s_D <5>[4.0] = c%d,\n",VName[j],k++);
                }
                else {
                    s1 = AcK[l * 4 + 1];
                    s2 = AcK[l * 4 + 2];

                    fprintf(PMTDAFd,"  E_%d_%d  <2>[4.0] = c%d,\n",s1,s2,k++);
                    fprintf(PMTDAFd,"  ED_%d_%d <2>[4.0] = c%d,\n",s1,s2,k++);

                    for (j = 0; j <= AcK[l * 4 + 3]; ++j)  
                        fprintf(PMTDAFd,"  E%d_%d_%d <2>[4.0] = c%d,\n",j,s1,s2,k++);
                }
            }
        }
        fprintf(PMTDAFd,");\n");
        printf1("TDA description written to: %s\n",PMTDAFName);
    }
    err = 0;

SEQMDFin:
    p_clean();
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

int get_nxe(int n,int *na)
{
    register int i,l;
    int err,sn,j,k,d;
    register char *p,*q;
    char vname[VNLMax + 1];

    err = -1;
    p = PMRHSTR;
         printf1("p:%s:\n",p);      

    printf1("\nDefinition       Variables\n");
    prnchar('-',26,1);
    for (i = 0; i < n; ++i) {
   
        if (*p != '[') {
         printf1("1p:%s:\n",p);      
            err = -2;
            goto GETNXEFin;
        }
        q = p + 1;
        while (*q && *q != ']')
            q++;
        if (*q++ != ']') {
         printf1("2p:%s:\n",p);      
            err = -2;
            goto GETNXEFin;
        }
        *q = '\0';
        printf1("%s",p);
        prnchar(' ',17 - strlen(p),0);

        if ((j = get_vnlen(p + 1)) > 1) {
            k = get_vidx1(p + 1,vname);
            if (k < 0) {
                printf1("error: variable does not exist.\n");
                goto GETNXEFin;
            }
            /*************
            if (j > 14) {
                printf1("error: exceeded max length of variable name.\n");
                goto GETNXEFin;
            }
            *******************/
            na[i * 4] = -1;
            na[i * 4 + 1] = k;                 

            printf1("%s_D\n",vname);
        }
        else if (sscanf(p,"[%d,%d,%d,%d]",&sn,&j,&k,&d) == 4) {
            sn--;
            if (sn < 0 || sn > SeqDNH || SeqDT[sn] == 0) {
                printf1("error: sequence %d not defined.\n",sn + 1);
                goto GETNXEFin;
            }
            if (j < 0 || k < 0 || j >= SeqSTNH[sn] || k >= SeqSTNH[sn]) {
                printf1("error: event [%d,%d] not defined.\n",j,k);
                goto GETNXEFin;
            }
            na[i * 4] = sn;
            na[i * 4 + 1] = j;
            na[i * 4 + 2] = k;
            na[i * 4 + 3] = d;                 

            printf1("E_%d_%d  ",j,k);
            printf1("ED_%d_%d  ",j,k);
            for (l = 0; l <= d; ++l)  
                printf1("E%d_%d_%d  ",l,j,k);
            printf1("\n");
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
        printf1("Syntax error in definition of event variables.\n");
    return(err);
}

