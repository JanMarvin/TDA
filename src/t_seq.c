/****************************************************************************/
/*  t_seq                                                                   */
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
#include "t_rand.h"
#include "t_gf.h"
#include "t_parm.h"
#include "t_alloc.h"
#include "t_seqe.h"
#include "t_seqm.h"

/*  functions in t_seq.c */

int t_seq(void);
int check_sdj(int j);
int get_sdata(int i,int j,int s);
int seq_sget(int i,int t,int n);
int seq_tsget(int i,int n);
int seq_tfget(int i,int n);
int seq_getsn(int sn,int opt);
int seq_scheck(int i,int j);
void seq_free(int i);                   
void seq_afree(int opt);                   
int seqdel(void);                   
int seqdef(void);
int seq_rc(int s,int smax);
void seqadj(void);
void prnsd(void);
int prnsdi(void);                   
int seqlg(void);
int seqld(void);
void prn_ssel(int n);
int seq_csel(void);
void seq_dtda(int typ,int sn,char *fname,int noc,int nv,short *vidx);
int seqgc(void);
int seqsd(void);
int seqen(void);
int seqtp(void);
int seqrd(void);
int seqsi(void);
int seqpe(void);
void seqpe1(int i,int id,int *st);
int seqpd(void);
int seqsn(void);
int seqgm(void);
int seqgsp(void);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

int SeqDN = 0;              /* number of sequences                          */  
int SeqDNH = 0;             /* highest sequence number                      */
char SeqDT[SEQ_Max];        /* type of sequence                             */
short SeqDNV[SEQ_Max];      /* number of variables                          */
short *SeqDV[SEQ_Max];      /* indices of variables                         */
int SeqGL[SEQ_Max];         /* max length of internal gaps                  */

short SeqSTN[SEQ_Max];      /* number of different states                   */
short *SeqSTNI[SEQ_Max];    /* array with SeqSTN[] state numbers            */
short *SeqSTNII[SEQ_Max];   /* array with SeqSTNH[] inverse state numbers   */
short SeqSTNH[SEQ_Max];     /* highest state number plus 1                  */

int SeqSTNHMax = 0;         /* highest state number in all sequences        */
int SeqTMin[SEQ_Max];       /* minimum of time points                       */
int SeqTMax[SEQ_Max];       /* maximum of time points                       */
int SeqTCMin = 0;           /* minimum of time points for all sequences     */
int SeqTCMax = 0;           /* maximum of time points for all sequences     */

short *SeqRC[SEQ_Max];      /* array with recode information                */
short SeqRCA[SEQ_Max];      /* length of SeqRC if allocated                 */

/*--------------------------------------------------------------------------*/
/*  t_seq()     Entry point for sequence commands.                          */
/*              Return 0 if OK, -1 if error, 1 if command cannot be         */
/*              interpreted.                                                */

int t_seq(void)
{
    char *p;
    int err = 0;

    p = CmdBuf;

    if (!strncmp(p,"seqdef",6))         /* new sequence data */
        err = seqdef();                 
    else if (!strncmp(p,"seqdel",6))    /* delete sequence data */
        err = seqdel();                 
    else if (!strncmp(p,"seqlg",5))     /* sequence length */
        err = seqlg();                 
    else if (!strncmp(p,"seqld",5))     /* last occurrence */
        err = seqld();                 
    else if (!strncmp(p,"seqgc",5))     /* sequence characteristics */
        err = seqgc();                 
    else if (!strncmp(p,"seqsd",5))     /* state distributions */
        err = seqsd();                 
    else if (!strncmp(p,"seqen",5))     /* entropy */
        err = seqen();                 
    else if (!strncmp(p,"seqtp",5))     /* transition probabilities */
        err = seqtp();                 
    else if (!strncmp(p,"seqrd",5))     /* random sequences */
        err = seqrd();                 
    else if (!strncmp(p,"seqsi",5))     /* state indicator matrix */
        err = seqsi();                 
    else if (!strncmp(p,"seqpe",5))     /* create sequence from episode data */
        err = seqpe();                 
    else if (!strncmp(p,"seqpd",5))     /* print sequence data */
        err = seqpd();                 
    else if (!strncmp(p,"seqevd",6))    /* data file with events */
        err = seqevd();                 
    else if (!strncmp(p,"seqev",5))     /* info about events */
        err = seqev();                 
    else if (!strncmp(p,"seqmd",5))     /* data for event models */
        err = seqmd();                 
    else if (!strncmp(p,"seqgm",5))     /* group membership */
        err = seqgm();                 
    else if (!strncmp(p,"seqm",4))      /* optimal matching */
        err = seqm();                 
    else if (!strncmp(p,"seqpm",5))     /* pattern matching */
        err = seqpm();                 
    else if (!strncmp(p,"seqsn",5))     /* states as new objects */
        err = seqsn();                 
    else if (!strncmp(p,"seqgsp",6))    /* generate plot data */
        err = seqgsp();                 

    else if (!strcmp(p,"seq"))          /* print info */
        err = prnsdi();                 
    else
        err = 1;

    return(err);
}

/* ------------------------------------------------------------------------ */
/*  check_sdj(j)    Check whether variable j is needed for sequence data.   */
/*                  Return 1 if needed, otherwise 0.                        */

int check_sdj(int j)
{
    register int i,k,n;

    if (SeqDN == 0)
        return(0);

    for (i = 0; i <= SeqDNH; ++i) {
        n = SeqDNV[i];
        for (k = 0; k < n; ++k) {
            if (j == SeqDV[i][k])
                return(1);
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  get_sdata(i,j,n)    get sequence data for variable i, case j, for       */
/*                      sequence n. Observe recoding of states.             */

int get_sdata(int i,int j,int n)
{
    register int s;

    s = (int)get_data(i,j);
    if (s >= -9 && SeqRCA[n] > 0)  
       s = SeqRC[n][s + 9];
    return(s);
}

/* ------------------------------------------------------------------------ */
/*  seq_sget(i,t,n)   Get state for case i at time t; -1 if not present.    */
/*                    n is number of sequence.                              */

int seq_sget(int i,int t,int n)
{
    register int j,s,nv;
    int ts,tf,ta,tb,sa,sb;

    if (n >= 0 && n <= SeqDNH) { 
        if (SeqDT[n] == 1) {
            if (t < 0 || t >= SeqDNV[n])
                return(-1);
            s = get_sdata(SeqDV[n][t],i,n);

            if (s < 0 && SeqGL[n] > 0) {
                ta = t - 1;
                sa = -1;
                while (ta >= 0) {
                    sa = get_sdata(SeqDV[n][ta],i,n);
                    if (sa >= 0)
                        break;
                    ta--;
                }   
                tb = t + 1;
                sb = -1;
                while (tb < SeqDNV[n]) {
                    sb = get_sdata(SeqDV[n][tb],i,n);
                    if (sb >= 0)
                        break;
                    tb++;
                }   
                if (sa >= 0 && sa == sb && tb - ta -1 <= SeqGL[n])
                    s = sa;
            }
            return(s);
        }
        else  {
            s = -1;
            ts = (int)get_data(SeqDV[n][1],i);
            if (t >= ts) {
                nv = SeqDNV[n];
                tf = (int)get_data(SeqDV[n][nv - 1],i);
                if (t == tf)
                    s = get_sdata(SeqDV[n][nv - 2],i,n);
                else if (t < tf) {
                    for (j = 3; j < nv; j += 2) {
                        if (t < (int)get_data(SeqDV[n][j],i)) {
                            s = get_sdata(SeqDV[n][j - 3],i,n);
                            break;
                        }
                    }
                }
            }
            return(s);
        }
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  seq_tsget(i,n)  Get first time point in sequence n for individual i.    */

int seq_tsget(int i,int n)
{
    register int t,typ;

    if (n >= 0 && n <= SeqDNH) { 
        typ = SeqDT[n];
        if (typ == 1) {
            for (t = 0; t < SeqDNV[n]; ++t) {
                if (get_sdata(SeqDV[n][t],i,n) >= 0)
                    return(t);
            }
        }
        else if (typ == 2) {
            t = (int)get_data(SeqDV[n][1],i);
            return(t);
        }
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  seq_tfget(i,n)  Get last time point in sequence n for individual i.     */

int seq_tfget(int i,int n)
{
    register int t,typ,nv;

    if (n >= 0 && n <= SeqDNH) { 
        typ = SeqDT[n];
        nv = SeqDNV[n];

        if (typ == 1) {
            for (t = nv - 1; t >= 0; --t) {
                if (get_sdata(SeqDV[n][t],i,n) >= 0)
                    return(t);
            }
        }
        else if (typ == 2) {
            t = (int)get_data(SeqDV[n][nv - 1],i);
            return(t);
        }
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  seq_getsn(sn,opt)                                                       */
/*                                                                          */
/*                  Get sequence number based on sn. If sn < 0 get number   */
/*                  number of first defined sequence. Otherwise check       */
/*                  whether sn is an existing sequence. Return sequence     */
/*                  number, or -1 if not available. If not available and    */
/*                  opt != 0, print error message.                          */
/*                                                                          */
/*  Note: sn is given as internal sequence number + 1.                      */

int seq_getsn(int sn,int opt)
{
    register int k;
    sn--;
    for (k = 0; k <= SeqDNH; ++k) {
        if (SeqDT[k]) {
            if (sn < 0 || sn == k) {
                if (opt) {
                    printf1("Using sequence data structure %d.\n",k + 1);
                }
                return(k);
            }
        }
    }
    if (opt && sn >= 0)  
        printf1("Error: sequence data structure %d not defined.\n",sn + 1);
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  seq_scheck(i,j)     Check whether sequences i and j have the same       */
/*                      state space. Return 0 if so, otherwise -1.          */

int seq_scheck(int i,int j)
{
    register int k,n;

    n = SeqSTN[i];
    if (n != SeqSTN[j])
        return(-1);
    for (k = 0; k < n; ++k) {
        if (SeqSTNI[i][k] != SeqSTNI[j][k])
            return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  seq_free(i)  free sequence data structure i.                            */
   
void seq_free(int i)                    
{
    if (SeqDNV[i] > 0) {
        free((char *)SeqDV[i]);
        memrq(-SeqDNV[i],sizeof(short));
        SeqDNV[i] = 0;
    }   
    if (SeqSTNH[i] > 0) {
        free((char *)SeqSTNII[i]);
        memrq(-SeqSTNH[i],sizeof(short));
        SeqSTNH[i] = 0;
    }   
    if (SeqSTN[i] > 0) {
        free((char *)SeqSTNI[i]);
        memrq(-SeqSTN[i],sizeof(short));
        SeqSTN[i] = 0;  
    }   
    if (SeqRCA[i] > 0) {
        free((char *)SeqRC[i]);
        memrq(-SeqRCA[i],sizeof(short));
        SeqRCA[i] = 0;
    }   
    SeqGL[i] = SeqDT[i] = 0;
    seqadj();           /* adjust number of sequences */
} 

/* ------------------------------------------------------------------------ */
/*  seq_afree(opt)  free all sequence data structures.                      */
/*                  if opt != 0 print message.                              */
   
void seq_afree(int opt)                 
{
    int i,n;

    n = 0;
    for (i = 0; i < SEQ_Max; ++i) {
        if (SeqDT[i]) {
            seq_free(i);
            n++;
        }
    }
    if (n > 0 && opt) 
        printf1("Deallocated %d sequence data structure(s).\n",n);
}

/* ------------------------------------------------------------------------ */
/*  seqdel()    delete sequence data.                                       */
/*              Return 0 if OK, -1 if error.                                */

int seqdel(void)                    
{
    int n;

    if (check_cmd(1))
        return(-1);

    if (SeqDN == 0) {
        printf1("No sequences defined.\n");
        return(0);
    }
    if (sscanf(CmdBuf,"seqdel=%d",&n) == 1) {
        if (n < 1 || n > SEQ_Max || SeqDT[n - 1] == 0)
            printf1("Error: there is no sequence data structure %d (command ignored).\n",n);
        else {
            seq_free(n - 1);
            printf1("Deallocated sequence data structure %d. Current memory: %d bytes.\n",n,MemReq);
        }
    }
    else if (!strcmp(CmdBuf,"seqdel")) {
        seq_afree(1);
        printf1("Current memory: %d bytes.\n",MemReq);
    }
    else {
        p_err(-1,1);
        return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  seqdef()    Create new sequence data. Command in CmdBuf.                */
/*                                                                          */
/*              seqdef(                                                     */
/*                  sn=...,         sequence number                         */
/*                  m=...,          type of definition                      */
/*                  rc=...,         recode info                             */
/*                  glen...,        max internal gap length, def. 0         */
/*              ) = varlist;                                                */
/*                                                                          */  
/*              m = 1 : varlist = Y0,Y1,... (discrete states)               */
/*              m = 2 : varlist = Y1,T1,... (discrete states)               */
/*              m = 3 : varlist = Y0,Y1,... (real-valued states)            */
/*              Note: m=3 not yet implemented.                              */
/*                                                                          */
/*  Return 0 if successful, -1 if error.                                    */

int seqdef(void)
{
    register int i,j;
    int err,s,smax,sn,nmiss,t,tl,tmin,tmax;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Creating a new sequence data structure. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 6,4,1))     /* get parameters */
        goto SEQNDFin;

    if (PMSN < 1)
        PMSN = 1;
    else if (PMSN > SEQ_Max) {
        printf1("Error: range of sequence data structures is: 1 - %d\n",SEQ_Max);
        goto SEQNDFin;
    }
    printf1("Sequence structure number: %d\n",PMSN);

    PMSN--;     /* internally we use sequence numbers 0,1,... */

    printf1("Sequence type: %d\n",PMM);
    if (PMM < 1 || PMM > 2) {
        printf1("Error: only m = 1 or 2 is possible.\n");
        goto SEQNDFin;
    }
    if (PMM == 2) {
        sn = PMNV / 2;
        if (2 * sn != PMNV) {
            printf1("Error: type 2 requires an even number of variables.\n");
            goto SEQNDFin;
        }
    }

    if (SeqDT[PMSN] != 0) {
        seq_free(PMSN);
        printf1("Previously defined sequence data structure %d will be deleted.\n",PMSN + 1);
    }

    /* PMNV = number of new variables on right-hand side */

    SeqDT[PMSN] = PMM;

    /* if PMGLEN > 0 set SeqGL */

    if (PMGLEN > 0) {
        if (PMM != 1) {
            printf1("Error: glen option only possible with m = 1.\n");
            goto SEQNDFin;
        }
        SeqGL[PMSN] = PMGLEN;
    }
    if (!(SeqDV[PMSN] = (short *)calloc(PMNV,sizeof(short)))) { 
        p_err(-2,1);
        goto SEQNDFin;
    }
    memrq(PMNV,sizeof(short));
    SeqDNV[PMSN] = PMNV;

    for (i = 0; i < PMNV; ++i)
        SeqDV[PMSN][i] = PMVIdx[i];

    tmax = 0;
    tmin = INTMAX;
    nmiss = 0;
    smax = -1;

    for (i = 0; i < NOC; ++i) {     /* get highest state number */
        tl = 0;                     /* and tmin,tmax for type 2 data */
        j = 0;                      /* also check ranges */
        while (j < PMNV) {
            s = (int)get_data(PMVIdx[j],i);
            if (s < 0)  
                nmiss++;
            else if (smax < s)
                smax = s;

            j++;
            if (PMM == 2) {
                t = (int)get_data(PMVIdx[j],i);
                if ((j == 1 && t < 0) || (j > 1 && t <= tl)) {
                    printf1("Error in case %d: need ascending event times.\n",i + 1);
                    goto SEQNDFin;
                }
                tl = t;
                if (tmin > t)
                    tmin = t;
                if (tmax < t)
                    tmax = t;
                j++;
            }
        }
    }
    if (smax < 0) {
        printf1("Error: cannot find any valid states.\n");
        goto SEQNDFin;
    }
    if (alloc_acn(smax + 1))
        goto SEQNDFin; 

    for (i = 0; i < NOC; ++i) {         /* get state numbers */
        for (j = 0; j < PMNV; ++j) {
            s = (int)get_data(PMVIdx[j],i);
            if (s >= 0)  
                AcN[s] += 1;
            if (PMM == 2)
                j++;
        }
    }

    if (PRCAlloc) {                  /* get recode information */
        if (seq_rc(PMSN,smax))
            goto SEQNDFin;
        j = 0;
        for (i = 0; i <= smax + 9; ++i) {
            if (j < SeqRC[PMSN][i])
                j = SeqRC[PMSN][i];
        }
        smax = j;
        if (smax < 0) {
            printf1("Error: no valid states after recoding.\n");
            goto SEQNDFin;
        }
        if (alloc_acn(smax + 1))
            goto SEQNDFin; 

        for (i = 0; i < NOC; ++i) {         /* get state numbers  */
            for (j = 0; j < PMNV; ++j) {
                s = (int)get_data(PMVIdx[j],i);
                if (s >= -9) {
                    s = SeqRC[PMSN][s + 9];
                        if (s >= 0) 
                              AcN[s] += 1;
                }
                if (PMM == 2)
                    j++;
            }
        }
    }
    if (!(SeqSTNII[PMSN] = (short *)calloc(smax + 1,sizeof(short)))) { 
        p_err(-2,1);
        goto SEQNDFin;
    }
    memrq(smax + 1,sizeof(short));
    SeqSTNH[PMSN] = smax + 1;

    sn = 0;                         /* number of different states */
    for (j = 0; j <= smax; ++j) {
        SeqSTNII[PMSN][j] = -1;
        if (AcN[j] > 0)  
            sn++;
    }
    if (!(SeqSTNI[PMSN] = (short *)calloc(sn,sizeof(short)))) { 
        p_err(-2,1);
        goto SEQNDFin;
    }
    memrq(sn,sizeof(short));
    SeqSTN[PMSN] = sn;

    i = 0;                          /* internal state number */
    for (j = 0; j <= smax; ++j) {
        if (AcN[j] > 0) {
            SeqSTNI[PMSN][i] = j;
            SeqSTNII[PMSN][j] = i;
            i++;
        }
    }
    /*************** 
    printf("SeqSTNI\n");
    for (i = 0; i < SeqSTN[PMSN]; ++i) {
        printf1("i=%d stni[][]=%d\n",i,SeqSTNI[PMSN][i]);
    }
    printf("SeqSTNII\n");
    for (i = 0; i < SeqSTNH[PMSN]; ++i) {
        printf1("i=%d stnii[][]=%d\n",i,SeqSTNII[PMSN][i]);
    }
    printf1("SeqRCA[PMSN] = %d\n",SeqRCA[PMSN]);
    printf("SeqRC\n");
    for (i = 0; i < SeqRCA[PMSN]; ++i)
        printf1("i=%d rc=%d\n",i,SeqRC[PMSN][i]);
    *************/
   
    if (PMM == 1) {                 /* set time axis */
        SeqTMin[PMSN] = 0;
        SeqTMax[PMSN] = PMNV - 1;
    }
    else if (PMM == 2) {
        SeqTMin[PMSN] = tmin;
        SeqTMax[PMSN] = tmax;
    }
    seqadj();           /* adjust number of sequences */
    prnsd();            /* print info */
    err = 0;

SEQNDFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  seq_rc(s,smax)      Get state recode information for sequence s.        */
/*                      smax is max number of states.                       */
/*                      Return 0 if OK, -1 if error.                        */

int seq_rc(int s,int smax)
{
    register int i;
    register char *p;
    int err,n,m;

    err = -2;
    if (!(SeqRC[s] = (short *)calloc(smax + 10,sizeof(short)))) { 
        p_err(-2,1);
        goto SEQRCFin;
    }
    memrq(smax + 10,sizeof(short));
    SeqRCA[s] = smax + 10;

    n = -9;
    for (i = 0; i <= smax + 9; ++i)  
        SeqRC[s][i] = n++;             
       
    p = PRC;
    err = -1;
               
    while (*p) {
        if (sscanf(p,"%d",&n) != 1 || n < -9) {
              goto SEQRCFin;
        }
        p = skip_int(p);
        if (!*p || *p++ != '[' || !*p)
            goto SEQRCFin;
        while (*p) {
            if (sscanf(p,"%d",&m) != 1 || m < -9 || m > smax)
                goto SEQRCFin;

            SeqRC[s][m + 9] = n; 
            p = skip_int(p);
            if (*p != ',')
                break;           
            p++;
        }
        if (*p++ != ']')
            goto SEQRCFin;
        if (!*p || *p != ',' || sscanf(p + 1,"%d",&n) != 1)
            break;
        p++;
    }
    if (!*p)
        err = 0;

SEQRCFin:
    if (err == -1)  
        printf1("Syntax error or unknown states in recode option.\n");
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  seqadj()    Adjust number of sequences and highest sequence number.     */

void seqadj(void)
{
    register int i;

    SeqTCMin = INTMAX;
    SeqTCMax = SeqDNH = SeqDN = 0;
    for (i = 0; i < SEQ_Max; ++i) {
        if (SeqDT[i]) {
            SeqDN++;
            if (SeqDNH < i)
                SeqDNH = i;
            if (SeqTCMin > SeqTMin[i])
                SeqTCMin = SeqTMin[i];
            if (SeqTCMax < SeqTMax[i])
                SeqTCMax = SeqTMax[i];
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  prnsd()     Print info about new sequence data.                         */

void prnsd(void)
{
    register int i,j,k;
    int n,s,ss,first;

    if (SeqDN == 0) {
        printf1("No sequences defined.\n");
        return;
    }
    printf1("Currently defined sequences:\n\n");
    printf1("Sequence          State       Time axis       Number\n");
    printf1("Structure Type  Variables  Minimum  Maximum  of States  States\n");
    prnchar('-',62,1);

    for (j = 0; j <= SeqDNH; ++j) {
        if (SeqDT[j]) {
            n = SeqDNV[j] / (int)SeqDT[j];
            printf1("%6d    %3d   %7d  %9d %8d  %8d   ",
                    j + 1,SeqDT[j],n,SeqTMin[j],SeqTMax[j],SeqSTN[j]);
                       
            for (i = 0; i < SeqSTN[j]; ++i) {
                s = SeqSTNI[j][i];
                printf1(" %d",s);
                if (SeqRCA[j] > 0) {
                    first = 1;
                    for (k = 0; k < SeqRCA[j]; ++k) {
                        ss = k - 9;              
                        if (SeqRC[j][k] == s && s != ss) {
                            if (first) {
                                printf1("[%d",ss);
                                first = 0;
                            }
                            else  
                                printf1(",%d",ss);
                        }
                    }
                    if (first == 0)
                        printf1("]");
                }
            }
            newline();
        }
    }
    printf1("\nRange of common time axis: %d to %d.\n",SeqTCMin,SeqTCMax);
}

/* ------------------------------------------------------------------------ */
/*  prnsdi()    Print info about currently defined sequences.               */
/*              Return 0 if OK, -1 if error.                                */

int prnsdi(void)                    
{
    if (check_cmd(1))
        return(-1);

    prnsd();
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  seqlg()         Sequence length and gaps.                               */
/*                  seqlg(sn=,sel=,v=,dtda=) = fname                        */
/*                  Return 0 if OK, -1 if error.                            */

int seqlg(void)
{
    register int i,j,l;
    int err,n,sn,a,b,ai,bi,t,len,gl,nrec,ng,gmin,gmax;

    err = -1;
    if (check_cmd(1))
        return(-1);
    printf1("Sequence length and gaps.\n");
         
    if (parm(CmdBuf + 5,1,1))    /* get parameters */
        goto SEQLFin;
        
    sn = seq_getsn(PMSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQLFin;

    a = SeqTMin[sn];
    b = SeqTMax[sn];

    fprintf(PMFd,"#       starting  ending  sequence     gap  min gap  max gap   number\n");
    fprintf(PMFd,"# Case      time    time    length  length   length   length  of gaps\n");

    nrec = 0;
    for (i = 0; i < NOC; ++i) {

        if (eval_sve(i) == 0)   
            continue;

        len = gl = 0;
        ai = bi = -1;

        ng = 0;         /* number of gaps */
        gmin = INTMAX;  /* min gap length */
        gmax = 0;       /* max gap length */

        for (t = a; t <= b; ++t) {
            n = seq_sget(i,t,sn);
            if (n >= 0) {
                ai = t;
                break;
            }
        }
        if (ai >= 0) {
            for (t = b; t >= ai; --t) {
                n = seq_sget(i,t,sn);
                if (n >= 0) {
                    bi = t;
                    break;
                }
            }
            len = bi - ai + 1;
            l = 0;
            for (t = ai + 1; t <= bi; ++t) {
                n = seq_sget(i,t,sn);
                if (n < 0) { 
                    gl++;
                    if (l == 0) {
                        ng++;
                        l = 1;
                    }
                    else
                        l++;
                }
                else if (l > 0) {
                    if (gmin > l)
                        gmin = l;
                    if (gmax < l)
                        gmax = l;
                    l = 0;
                }   
            }
        }
        if (gmin == INTMAX)
            gmin = 0;

        fprintf(PMFd,"%8d %7d %7d %9d %7d %8d %8d %8d ",
                                         i + 1,ai,bi,len,gl,gmin,gmax,ng);
        for (j = 0; j < PMNV; ++j) {
            n = (int)PMVIdx[j];
            fprintf(PMFd,VPFmtS[n],get_data(n,i));
        }
        fprintf(PMFd,"\n");
        nrec++;
    }
    prn_ssel(nrec);
    printf1("%d records written to: %s\n",nrec,PMFdName);

    if (PMTDAFDef && nrec > 0)                      /* write TDA description */
        seq_dtda(1,sn,PMFdName,nrec,PMNV,PMVIdx);

    err = 0;

SEQLFin:
    p_clean(); 
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqld()         Last occurrence of states.                              */
/*                                                                          */
/*                  seqld(                                                  */
/*                      s=...,      state number, def. 0                    */
/*                      sn= ...,    sequence number                         */
/*                      sel=,       selection                               */
/*                      v=,         add variables                           */
/*                      dtda=       TDA description of output file          */
/*                  ) = fname;                                              */
/*                                                                          */
/*                  Return 0 if OK, -1 if error.                            */

int seqld(void)
{
    register int i,j,l;
    int err,n,sn,a,b,ai,bi,t,len,gl,nrec;

    err = -1;
    if (check_cmd(1))
        return(-1);
    printf1("Last occurrence of states.\n");
         
    if (parm(CmdBuf + 5,1,1))    /* get parameters */
        goto SEQLDFin;
        
    sn = seq_getsn(PMSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQLDFin;

    a = SeqTMin[sn];
    b = SeqTMax[sn];

    fprintf(PMFd,"#       starting  ending  sequence    last        total\n");
    fprintf(PMFd,"# Case      time    time    length  occurrence  occurrence\n");

/**********

    nrec = 0;
    for (i = 0; i < NOC; ++i) {

        if (eval_sve(i) == 0)   
            continue;

        len = gl = 0;
        ai = bi = -1;

        for (t = a; t <= b; ++t) {
            n = seq_sget(i,t,sn);
            if (n >= 0) {
                ai = t;
                break;
            }
        }
        if (ai >= 0) {
            for (t = b; t >= ai; --t) {
                n = seq_sget(i,t,sn);
                if (n >= 0) {
                    bi = t;
                    break;
                }
            }
            len = bi - ai + 1;
            l = 0;
            for (t = ai + 1; t <= bi; ++t) {
                n = seq_sget(i,t,sn);
                if (n < 0) { 
                    gl++;
                    if (l == 0) {
                        ng++;
                        l = 1;
                    }
                    else
                        l++;
                }
                else if (l > 0) {
                    if (gmin > l)
                        gmin = l;
                    if (gmax < l)
                        gmax = l;
                    l = 0;
                }   
            }
        }
        if (gmin == INTMAX)
            gmin = 0;

        fprintf(PMFd,"%8d %7d %7d %9d %7d %8d %8d %8d ",
                                         i + 1,ai,bi,len,gl,gmin,gmax,ng);
        for (j = 0; j < PMNV; ++j) {
            n = (int)PMVIdx[j];
            fprintf(PMFd,VPFmtS[n],get_data(n,i));
        }
        fprintf(PMFd,"\n");
        nrec++;
    }
***********/


    prn_ssel(nrec);
    printf1("%d records written to: %s\n",nrec,PMFdName);

    if (PMTDAFDef && nrec > 0)                      /* write TDA description */
        seq_dtda(1,sn,PMFdName,nrec,PMNV,PMVIdx);

    err = 0;

SEQLDFin:
    p_clean(); 
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  prn_ssel(n)   If sel option print number of selected cases.             */

void prn_ssel(int n)
{
    if (SVEFlg)  
        printf1("Number of selected cases: %d\n",n);
}

/*--------------------------------------------------------------------------*/
/*  seq_csel()      return number of selected cases.                        */

int seq_csel(void)
{
    register int i;
    int n = 0;

    if (SVEFlg == 0)
        return(NOC);

    for (i = 0; i < NOC; ++i) {
        if (eval_sve(i) == 0)   
            continue;
        n++;
    }
    if (n < 1)  
        printf1("Error: no cases selected.\n");
    return(n);
}

/*--------------------------------------------------------------------------*/
/*  seq_dtda(typ,sn,fname,noc,nv,vidx)     Write TDA description file.      */
/*                                                                          */

void seq_dtda(int typ,int sn,char *fname,int noc,int nv,short *vidx)
{
    int i,j,n,ns,s,t;

    fprintf(PMTDAFd,"nvar(\n");
    fprintf(PMTDAFd,"  dfile = %s,\n",fname);
    fprintf(PMTDAFd,"  noc = %d,\n",noc);

    switch (typ) {

        case  1:                    /* seqlg */
            fprintf(PMTDAFd,"  CASE <5>[6.0] = c1, # case number\n");
            fprintf(PMTDAFd,"  TS   <5>[6.0] = c2, # starting time\n");
            fprintf(PMTDAFd,"  TF   <5>[6.0] = c3, # ending time\n");
            fprintf(PMTDAFd,"  SLEN <5>[6.0] = c4, # sequence length\n");
            fprintf(PMTDAFd,"  GLEN <5>[6.0] = c5, # total gap length\n");
            fprintf(PMTDAFd,"  GMIN <5>[6.0] = c6, # minimum gap length\n");
            fprintf(PMTDAFd,"  GMAX <5>[6.0] = c7, # maximum gap length\n");
            fprintf(PMTDAFd,"  NGAP <5>[6.0] = c8, # number of gaps\n");
            n = 9;
            break;

        case  2:                    /* seqgc */

            ns = SeqSTN[sn];        /* number of states */

            fprintf(PMTDAFd,"  CASE <5>[8.0] = c1 , # case number\n");
            fprintf(PMTDAFd,"  SLEN <2>[4.0] = c2 , # sequence length\n");
            fprintf(PMTDAFd,"  NDS  <2>[4.0] = c3 , # number of different states\n");
            fprintf(PMTDAFd,"  NEV  <2>[4.0] = c4 , # number state changes\n");
            n = 5;
            for (j = 0; j < ns; ++j) {
                s = SeqSTNI[sn][j];
                fprintf(PMTDAFd,"  DUR%-2d<2>[4.0] = c%-2d, # duration in state %d\n",s,n++,s);
            }
            fprintf(PMTDAFd,"  DURM <2>[4.0] = c%-2d, # duration in missing state\n",n++);

            for (j = 0; j < ns; ++j) {
                s = SeqSTNI[sn][j];
                fprintf(PMTDAFd,"  NEP%-2d<2>[4.0] = c%-2d, # number of episodes in state %d\n",s,n++,s);
            }
            fprintf(PMTDAFd,"  NEPM <2>[4.0] = c%-2d, # number of episodes in missing state\n",n++);
            break;

        case  3:                    /* seqsd */

            ns = SeqSTN[sn];        /* number of states */

            fprintf(PMTDAFd,"  TIME  <2>[6.0] = c1, # time\n");
            n = 2;
            for (j = 0; j < ns; ++j) {
                s = SeqSTNI[sn][j];
                fprintf(PMTDAFd,"  NST%-2d <2>[6.0] = c%d, # cases in state %d\n",s,n++,s);
            }
            fprintf(PMTDAFd,"  VALID <2>[6.0] = c%d, # cases in valid states\n",n++);
            fprintf(PMTDAFd,"  NMISS <2>[6.0] = c%d, # cases in missing state\n",n++);
            fprintf(PMTDAFd,"  TOTAL <2>[6.0] = c%d, # total number of cases\n",n++);
            break;

        case  4:                    /* seqen */
            fprintf(PMTDAFd,"  TIME <2>[ 6.0] = c1, # time\n");
            fprintf(PMTDAFd,"  N    <2>[ 6.0] = c2, # number of cases\n");
            fprintf(PMTDAFd,"  ENT  <4>[%d.%d] = c3, # entropy\n",PMFmt1,PMFmt2);
            nv = 0;
            break;

        case  5:                    /* seqtp */

            ns = SeqSTN[sn];        /* number of states */

            fprintf(PMTDAFd,"  TIME  <2>[ 6.0] = c1 , # time\n");
            fprintf(PMTDAFd,"  N     <2>[ 6.0] = c2 , # number of cases\n");
            n = 3;
            for (i = 0; i < ns; ++i) {
                s = SeqSTNI[sn][i];
                for (j = 0; j <= ns; ++j) {
                    if (j < ns) {
                        t = SeqSTNI[sn][j];
                        fprintf(PMTDAFd,"  TP%d_%d <4>[%d.%d] = c%-2d, # transition %d - %d\n",
                                                 s,t,PMFmt1,PMFmt2,n++,s,t);
                    }
                    else {
                        fprintf(PMTDAFd,"  TP%d_M <4>[%d.%d] = c%-2d, # transition %d - missing\n",
                                                 s,PMFmt1,PMFmt2,n++,s);
                    }
                }
            }
            nv = 0;
            break;

        case  6:                    /* seqsi */

            ns = SeqSTN[sn];        /* number of states */

            fprintf(PMTDAFd,"  Case <2>[6.0] = c1 , # case number\n");
            n = 2;
            for (i = 0; i < PMNTP; ++i) {
                for (j = 0; j < ns; ++j) {
                    s = SeqSTNI[sn][j];
                    t = (int)PMTP[i];
                    fprintf(PMTDAFd,"  S%d_%d <2>[2.0] = c%-2d, # state %d, time %d\n",s,t,n++,s,t);
                }
            }
            break;

        default:
            break;
    }
    for (i = 0; i < nv; ++i) {
        j = vidx[i];
        fprintf(PMTDAFd,"  %s <%d>[%d.%d]",VName[j],VSLen[j],VPFmt1[j],VPFmt2[j]);
        if (VLabel[j] != NULL)                                                      
            fprintf(PMTDAFd,"(%s)",VLabel[j]);                                           
        fprintf(PMTDAFd," = c%d,\n",n++);
    }
    fprintf(PMTDAFd,");\n");
    
    printf1("TDA description written to: %s\n",PMTDAFName);
}

/*--------------------------------------------------------------------------*/
/*  seqgc()         Characteristics of sequences                            */
/*                  seqgc(sn=,sel=,v=,dtda=) = fname                        */
/*                  Return 0 if OK, -1 if error.                            */

int seqgc(void)
{
    register int i,j,k;
    int err,nrec,m,n,n1,nc,ns,sn,a,b,ai,bi,t,dm;

    err = -1;         
    if (check_cmd(1))
        return(-1);
         
    printf1("Sequence characteristics.\n");
         
    if (parm(CmdBuf + 5,1,1))    /* get parameters */
        goto SEQFFin;
        
    sn = seq_getsn(PMSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQFFin;

    ns = SeqSTN[sn];        /* number of states */

    if (alloc_acn(ns + 1))
        goto SEQFFin;

    a = SeqTMin[sn];
    b = SeqTMax[sn];
    nrec = m = 0;
    for (i = 0; i < NOC; ++i) {

        if (eval_sve(i) == 0)   
            continue;

        ai = bi = -1;
        for (t = a; t <= b; ++t) {
            n = seq_sget(i,t,sn);
            if (n >= 0) {
                ai = t;
                break;
            }
        }
        if (ai < 0)
            continue;

        for (t = b; t >= ai; --t) {
            n = seq_sget(i,t,sn);
            if (n >= 0) {
                bi = t;
                break;
            }
        }
        fprintf(PMFd,"%6d ",i + 1);         /* case number */
        fprintf(PMFd,"%4d ",bi - ai + 1);   /* sequence length */
        m++;

        for (j = 0; j < ns; ++j)
            AcN[j] = 0;

        nc = 0;
        n1 = 0;
        dm = 0;
        for (t = ai; t <= bi; ++t) {
            n = seq_sget(i,t,sn);
            if (t > ai && n != n1)      /* state changes */
                nc++;
            n1 = n;
            if (n >= 0) {
                n = SeqSTNII[sn][n];
                AcN[n] += 1;
            }
            else
                dm++;                   /* duration in missing state */
        }
        n = 0; 
        for (j = 0; j < ns; ++j) {      /* # of different states */
            if (AcN[j] > 0)
                n++;
        }
        fprintf(PMFd,"%4d ",n);
        fprintf(PMFd,"%4d ",nc);

        for (j = 0; j < ns; ++j)        /* duration in state ... */
            fprintf(PMFd,"%4d ",AcN[j]);
        fprintf(PMFd,"%4d ",dm);

        /* number of episodes */

        for (j = 0; j <= ns; ++j)
            AcN[j] = 0;

        n1 = seq_sget(i,ai,sn);
        k = SeqSTNII[sn][n1];
        AcN[k] += 1;

        for (t = ai + 1; t <= bi; ++t) {
            n = seq_sget(i,t,sn);
            if (n != n1) {
                if (n >= 0) {
                    k = SeqSTNII[sn][n];
                    AcN[k] += 1;
                }
                else
                    AcN[ns] += 1;
            }
            n1 = n;
        }
        for (j = 0; j <= ns; ++j)        /* number of episodes */
            fprintf(PMFd,"%4d ",AcN[j]);

        /* additional variables */

        for (j = 0; j < PMNV; ++j) {
            n = (int)PMVIdx[j];
            fprintf(PMFd,VPFmtS[n],get_data(n,i));
        }
        fprintf(PMFd,"\n");
        nrec++;
    }
    prn_ssel(nrec);
    printf1("%d records written to: %s\n",nrec,PMFdName);

    if (PMTDAFDef && nrec > 0)                      /* write TDA description */
        seq_dtda(2,sn,PMFdName,nrec,PMNV,PMVIdx);
    err = 0;

SEQFFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqsd()         State distribution.                                     */
/*                  seqsd(sn=,sel=,dtda=) = fname                           */
/*                  Return 0 if OK, -1 if error.                            */

int seqsd(void)
{
    register int i,t;
    int err,n,a,b,sn,ns,nrec,nn;

    err = -1;         
    if (check_cmd(1))
        return(-1);
         
    printf1("State distributions.\n");
         
    if (parm(CmdBuf + 5,1,1))    /* get parameters */
        goto SEQSDFin;
        
    sn = seq_getsn(PMSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQSDFin;

    ns = SeqSTN[sn];        /* number of states */

    if (alloc_acn(ns + 1))
        goto SEQSDFin;

    a = SeqTMin[sn];
    b = SeqTMax[sn];

    nn = seq_csel();
    if (nn < 1)
        goto SEQSDFin;
         
    nrec = 0;
    for (t = a; t <= b; ++t) {

        for (i = 0; i <= ns; ++i)
            AcN[i] = 0;

        fprintf(PMFd,"%6d ",t);

        for (i = 0; i < NOC; ++i) {
            if (eval_sve(i) == 0)   
                continue;
            n = seq_sget(i,t,sn);
            if (n >= 0)
                n = SeqSTNII[sn][n];
            else
                n = ns;
            AcN[n] += 1;
        }
        n = 0;
        for (i = 0; i < ns; ++i) {
            fprintf(PMFd,"%6d ",AcN[i]);
            n += AcN[i];
        }
        fprintf(PMFd,"%6d %6d %6d\n",n,AcN[ns],nn);
        nrec++;
    }
    prn_ssel(nn);
    printf1("%d records written to: %s\n",nrec,PMFdName);

    if (PMTDAFDef && nrec > 0)                      /* write TDA description */
        seq_dtda(3,sn,PMFdName,nrec,PMNV,PMVIdx);
    err = 0;

SEQSDFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqen()         Entropy measures.                                       */
/*                  seqen(sn=,sel=,dtda=,fmt=) = fname                      */
/*                  Return 0 if OK, -1 if error.                            */

int seqen(void)
{
    register int i,t;
    int err,m,n,a,b,sn,ns,nn,nrec;
    double tmp,tmp1;

    err = -1;         
    if (check_cmd(1))
        return(-1);
         
    printf1("Entropy measures.\n");
         
    pmfmt(10,4);                /* default print format */
    if (parm(CmdBuf + 5,1,1))   /* get parameters */
        goto SEQENFin;
        
    if (PMFmt1 == 0)
        pmfmt(10,4);            /* default print format */

    sn = seq_getsn(PMSN,1);     /* get sequence number */
    if (sn < 0)  
        goto SEQENFin;

    ns = SeqSTN[sn];            /* number of states */

    if (alloc_acn(ns + 1))
        goto SEQENFin;

    a = SeqTMin[sn];
    b = SeqTMax[sn];

    nn = seq_csel();
    if (nn < 1)
        goto SEQENFin;
       
    fprintf(PMFd,"# Time      N    Entropy\n");

    nrec = 0;
    for (t = a; t <= b; ++t) {

        for (i = 0; i < ns; ++i)
            AcN[i] = 0;         

        m = 0;
        for (i = 0; i < NOC; ++i) {

            if (eval_sve(i) == 0)   
                continue;

            n = seq_sget(i,t,sn);
            if (n >= 0) {
                n = SeqSTNII[sn][n];
                AcN[n] += 1;
                m++;
            }
        }
        tmp = 0.0;
        if (m > 0) {
            for (i = 0; i < ns; ++i) {
                if (AcN[i] > 0) {
                    tmp1 = (double)AcN[i] / (double)m;
                    tmp -= tmp1 * rlog(tmp1);
                }
            }
        }
        fprintf(PMFd,"%6d %6d ",t,m);
        fprintf(PMFd,PMFmtS,tmp);
        fprintf(PMFd,"\n");
        nrec++;
    }
    prn_ssel(nn);
    printf1("%d records written to: %s\n",nrec,PMFdName);

    if (PMTDAFDef && nrec > 0)                      /* write TDA description */
        seq_dtda(4,sn,PMFdName,nrec,PMNV,PMVIdx);

    err = 0;

SEQENFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqtp()         Transition probabilities                                */
/*                  seqtp(sn=,sel=,fmt=,dtda=) = fname                      */
/*                  Return 0 if OK, -1 if error.                            */

int seqtp(void)
{
    register int i,j,t;
    int err,m,mm,mm0,n,n1,a,b,sn,ns,ns1,nss,nn,nrec;
    double tmp;

    err = -1;         
    if (check_cmd(1))
        return(-1);

    printf1("Transition probabilities.\n");
         
    if (parm(CmdBuf + 5,1,1))    /* get parameters */
        goto SEQTPFin;
        
    sn = seq_getsn(PMSN,1);     /* get sequence number */
    if (sn < 0)  
        goto SEQTPFin;

    if (PMFmt1 == 0)
        pmfmt(10,4);            /* default print format */

    ns = SeqSTN[sn];            /* number of states */
    ns1 = ns + 2;
    nss = ns * ns1;

    if (alloc_acn(nss + 1))
        goto SEQTPFin;

    a = SeqTMin[sn];
    b = SeqTMax[sn];

    nn = seq_csel();
    if (nn < 1)
        goto SEQTPFin;

    nrec = 0;
    for (t = a; t < b; ++t) {

        for (i = 0; i < nss; ++i)
            AcN[i] = 0;         

        mm0 = mm = m = 0;
        for (i = 0; i < NOC; ++i) {

            if (eval_sve(i) == 0)   
                continue;

            n = seq_sget(i,t,sn);
            if (n >= 0) {
                n = SeqSTNII[sn][n];
                m++;
                n1 = seq_sget(i,t + 1,sn);
                if (n1 >= 0) 
                    n1 = SeqSTNII[sn][n1];
                else  
                    n1 = ns;
                AcN[n * ns1 + n1] += 1;
                AcN[n * ns1 + ns + 1] += 1;
                if (n1 >= 0) {
                    mm0++;
                    if (n1 != n)
                        mm++;
                }
            }
        }
        fprintf(PMFd,"%6d %6d ",t,m);

        for (i = 0; i < ns; ++i) {
            n1 = AcN[i * ns1 + ns + 1];
            for (j = 0; j <= ns; ++j) {
                if (n1 <= 0)
                    tmp = 0.0;
                else 
                    tmp = (double)AcN[i * ns1 + j] / (double)n1;
                fprintf(PMFd,PMFmtS,tmp);
            }
        }
        if (mm0 > 0)
            tmp = (double)mm / (double)mm0;
        else
            tmp = 0.0;
        fprintf(PMFd,PMFmtS,tmp);
        fprintf(PMFd,"\n");
        nrec++;
    }
    prn_ssel(nn);
    printf1("%d records written to: %s\n",nrec,PMFdName);

    if (PMTDAFDef && nrec > 0)                      /* write TDA description */
        seq_dtda(5,sn,PMFdName,nrec,PMNV,PMVIdx);

    err = 0;

SEQTPFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqrd()     Create random sequences.                                    */
/*              seqrd(ns=,noc=,len=)=fname.                                 */
/*              Return 0 if OK, -1 if error.                                */

int seqrd(void)
{
    register int i,j;
    int err,m,n;

    err = -1;         
    if (check_cmd(1))
        return(-1);

    printf1("Creating random sequences.\n");
         
    if (parm(CmdBuf + 5,1,1))    /* get parameters */
        goto SEQRDFin;
        
    if (PMNS < 1)
        PMNS = 2;

    if (PMLEN < 1)
        PMLEN = 1;

    printf1("Sequence length: %d, number of states: %d\n",PMLEN,PMNS);

    m = 0;
    for (i = 1; i <= PMNOC; ++i) {
        fprintf(PMFd,"%6d ",i);
        for (j = 0; j < PMLEN; ++j) {
            n = (int)floor(random1() * (double)PMNS) + 1;
            fprintf(PMFd,"%2d ",n);
        }
        fprintf(PMFd,"\n");
        m++;
    }
    printf1("%d sequences written to output file: %s\n",m,PMFdName);
    err = 0;

SEQRDFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqsi()     Create state indicator matrix.                              */
/*              seqsi(m=,tp=,sel=,dtda=,v=)=fname.                          */
/*              Return 0 if OK, -1 if error.                                */

int seqsi(void)
{
    register int i,j;
    int err,m,n,nt,ns,nst,sn,t;

    err = -1;         
    if (check_cmd(1))
        return(-1);

    printf1("Printing state indicator matrix.\n");
         
    if (parm(CmdBuf + 5,1,1))    /* get parameters */
        goto SEQSIFin;
        
    if (PMM != 2)
        PMM = 1;

    if (PMNTP < 1) {
        p_err(-17,1);
        goto SEQSIFin;
    }
    sn = seq_getsn(PMSN,1);     /* get sequence number */
    if (sn < 0)  
        goto SEQSIFin;

    ns = SeqSTN[sn];    /* number of states */
    nst = ns * PMNTP;

    if (alloc_acn(nst + 1))
        goto SEQSIFin;

    m = 0;
    for (i = 0; i < NOC; ++i) {

        if (eval_sve(i) == 0)   
            continue;

        for (j = 0; j < nst; ++j)  
            AcN[j] = 0;

        nt = 0;
        for (j = 0; j < PMNTP; ++j) {
            t = (int)PMTP[j];
            n = seq_sget(i,t,sn);
            if (n >= 0) {
                n = SeqSTNII[sn][n];     /* internal state number */
                AcN[j * ns + n] = 1;
                nt++;
            }
            else if (PMM == 2)
                break;
        }
        if (PMM == 1 || nt == PMNTP) {
            fprintf(PMFd,"%6d ",i + 1);
            for (j = 0; j < nst; ++j)  
                fprintf(PMFd,"%1d ",AcN[j]);

            for (j = 0; j < PMNV; ++j) {
                n = (int)PMVIdx[j];
                fprintf(PMFd,VPFmtS[n],get_data(n,i));
            }
            fprintf(PMFd,"\n");
            m++;
        }
    }
    prn_ssel(m);
    printf1("%d records written to output file: %s\n",m,PMFdName);

    if (PMTDAFDef && m > 0)                      /* write TDA description */
        seq_dtda(6,sn,PMFdName,m,PMNV,PMVIdx);

    err = 0;

SEQSIFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqpe() ##  Create sequence data from episode data.                     */
/*                                                                          */
/*              seqpe(                                                      */
/*                  id=,                                                    */
/*                  org=,                                                   */
/*                  ts=,                                                    */
/*                  tf=,                                                    */
/*                  tp=,                                                    */
/*                  m=...,      missing value, def. -1                      */
/*                  m1=...,     missing values at the beginning, def. -1    */
/*                  m2=...,     missing values in between, def. -1          */
/*                  v=                                                      */
/*              ) = fname;                                                  */
/*                                                                          */  
/*  Return 0 if OK, -1 if error.                                            */

int seqpe(void)
{
    register int i,j;
    int err,m,a,b,i1,id,id1,s;
    double t,ts,tf;

    err = -1;         
    if (check_cmd(1))
        return(-1);

    printf1("Creating sequences from episodes.\n");
         
    if (parm(CmdBuf + 5,1,1))    /* get parameters */
        goto SEQPEFin;
        
    if (PMMFlg == 0) 
        PMM = -1;
    if (PMM1Flg == 0) 
        PMM1 = -1;
    if (PMM2Flg == 0) 
        PMM2 = -1;

    if (PMID < 0 || PMORG < 0 || PMTS < 0 || PMTF < 0 || PMNTP < 1) {
        p_err(-4,1);
        goto SEQPEFin;
    }
    if (alloc_acn(PMNTP + 1))
        goto SEQPEFin;

    for (j = 0; j < PMNTP; ++j)
        AcN[j] = PMM;

    id1 = (int)get_data(PMID,0);
    i1 = 0;

    a = PMTP[0];
    b = PMTP[PMNTP - 1];
    m = 0;
                   
    for (i = 0; i < NOC; ++i) {

        id = (int)get_data(PMID,i);
        s  = (int)get_data(PMORG,i);
        ts = get_data(PMTS,i);
        tf = get_data(PMTF,i);

        if (id != id1) {
            m++,
            seqpe1(i1,m,AcN);
            id1 = id;
            i1 = i;
            for (j = 0; j < PMNTP; ++j)
                AcN[j] = PMM;
        }
        if (ts <= tf && ts <= b && tf >= a) {       
            for (j = 0; j < PMNTP; ++j) {
                t = PMTP[j];
                if (t >= ts && t < tf)
                    AcN[j] = s;
            }
        }           
    }
    m++;
    seqpe1(i1,m,AcN);
    printf1("%d records written to output file: %s\n",m,PMFdName);
    err = 0;

SEQPEFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqpe1(i,id,st)     print sequence data for individual i with id and    */  
/*                      states st[].                                        */
/*                      called by seqpe().                                  */

void seqpe1(int i,int id,int *st)
{
    register int j,iv;

    /* insert new missing value codes */

    if (PMM1Flg) {
        for (j = 0; j < PMNTP; ++j) {
            if (st[j] == PMM)
                st[j] = PMM1;
            else
                break;
        }
    }
    if (PMM2Flg) {
        for (j = PMNTP - 1; j >= 0; --j) {
            if (st[j] == PMM)
                st[j] = PMM2;
            else
                break;
        }
    }

    fprintf(PMFd,"%6d ",id);
    for (j = 0; j < PMNTP; ++j)  
        fprintf(PMFd,"%2d ",st[j]);

    for (j = 0; j < PMNV; ++j) {
        iv = (int)PMVIdx[j];
        fprintf(PMFd,VPFmtS[iv],get_data(iv,i));
    }
    fprintf(PMFd,"\n");
}

/*--------------------------------------------------------------------------*/
/*  seqpd()         Print sequence data                                     */
/*                  seqpd(                                                  */
/*                    m=...,            option, def. 1.                     */
/*                    ns=...,           ns = 1                              */
/*                    s=...,            aggregation of time points          */
/*                    sel=...,          expression for case selection       */
/*                    v=...,            add variables                       */
/*                    dtda=...,         TDA description file                */
/*                  ) = fname;          output file                         */
/*                                                                          */
/*  m = 1                                                                   */
/*      2                                                                   */
/*      3                                                                   */
/*      4                                                                   */
/*      5                                                                   */
/*      6   print sequence of states, disregard durations.                  */
/*      7   same as 6                                                       */
/*      8   aggregate s time points.                                        */
/*      9   individual total durations in different states                  */
/*     10   first and last time point in each state                         */
/*                                                                          */
/*  If ns = 1, the command creates a new state space that distinguishes     */
/*  the number of occurrences in repeatable states. (Only with m = 1)       */
/*                                                                          */  
/*  Return 0 if OK, -1 if error.                                            */

int seqpd(void)
{
    register int i,j,k,l;
    int m,n,nm,err,org,ts,sn,iv,min,max,ta,tb,ns,ng,gmin,gmax,n1,len;
    int sa,sb,lcen,rcen;
   
    err = -1;         
    if (check_cmd(1))
        return(-1);

    printf1("Printing sequence data.\n");

    if (SeqDN == 0) {       /* check for sequence data */
        p_err(-13,1);
        goto SEQPDFin;
    }
    if (parm(CmdBuf + 5,1,1))    /* get parameters */
        goto SEQPDFin;

    min = SeqTCMin; 
    max = SeqTCMax; 
        
    m = 0;         
    switch (PMM) {

        case 1:    
            for (i = 0; i < NOC; ++i) {

                if (eval_sve(i) == 0)   
                    continue;

                fprintf(PMFd,"%6d ",i + 1);

                for (j = 0; j <= SeqDNH; ++j) {
                    if (SeqDT[j]) {

                        if (PMNS == 1) {
                            nm = SeqSTNH[j];
                            if (alloc_acn(nm + 1))
                                goto SEQPDFin;
                            for (k = 0; k <= nm; ++k)
                                AcN[k] = -10;
                            n1 = -1;
                            for (k = min; k <= max; ++k) {
                                n = seq_sget(i,k,j);
                                if (n >= 0 && n < nm && n != n1)
                                    AcN[n] += 10;
                                n1 = n;
                                if (n >= 0)
                                    n += AcN[n];
                                fprintf(PMFd,"%3d ",n);
                            }
                        }
                        else {
                            for (k = min; k <= max; ++k) {
                                n = seq_sget(i,k,j);
                                fprintf(PMFd,"%2d ",n);
                            }
                        }
                    }
                }
                for (k = 0; k < PMNV; ++k) {
                    iv = (int)PMVIdx[k];
                    fprintf(PMFd,VPFmtS[iv],get_data(iv,i));
                }
                fprintf(PMFd,"\n");
                m++;
            }
            break;

        case 2:
            for (i = 0; i < NOC; ++i) {
                if (eval_sve(i) == 0)   
                    continue;

                for (j = 0; j <= SeqDNH; ++j) {
                    if (SeqDT[j]) {
                        fprintf(PMFd,"%6d %2d ",i + 1,j + 1);
                        for (k = min; k <= max; ++k) {
                            n = seq_sget(i,k,j);
                            fprintf(PMFd,"%2d ",n);
                        }
                        for (k = 0; k < PMNV; ++k) {
                            iv = (int)PMVIdx[k];
                            fprintf(PMFd,VPFmtS[iv],get_data(iv,i));
                        }
                        fprintf(PMFd,"\n");
                        m++;
                    }
                }
            }
            break;

        case 3:
            for (i = 0; i < NOC; ++i) {

                if (eval_sve(i) == 0)   
                    continue;

                for (k = min; k <= max; ++k) {
                    fprintf(PMFd,"%6d %4d ",i + 1,k);
                    for (j = 0; j <= SeqDNH; ++j) {
                        if (SeqDT[j]) {
                            n = seq_sget(i,k,j);
                            fprintf(PMFd,"%2d ",n);
                        }
                    }
                    for (j = 0; j < PMNV; ++j) {
                        iv = (int)PMVIdx[j];
                        fprintf(PMFd,VPFmtS[iv],get_data(iv,i));
                    }
                    fprintf(PMFd,"\n");
                    m++;
                }
            }
            break;

        case 4:

            for (j = 0; j <= SeqDNH; ++j) {
                if (SeqDT[j])  
                    break;
            }
            for (i = 0; i < NOC; ++i) {

                if (eval_sve(i) == 0)   
                    continue;

                sn = 1;
                ts = k = min;
                org = seq_sget(i,k,j);
                while (++k <= max + 1) {
                    if (k <= max)
                        n = seq_sget(i,k,j);
                    else
                        n = org;

                    if (n != org || k == max + 1) {
                        fprintf(PMFd,"%6d %2d %2d %2d %4d %4d ",i + 1,sn++,org,n,ts,k);
                        for (l = 0; l < PMNV; ++l) {
                            iv = (int)PMVIdx[l];
                            fprintf(PMFd,VPFmtS[iv],get_data(iv,i));
                        }
                        fprintf(PMFd,"\n");
                        m++;
                        org = n;
                        ts = k;
                    }
                }
            }
            break;
         
        case 5:
            for (i = 0; i < NOC; ++i) {

                if (eval_sve(i) == 0)   
                    continue;

                fprintf(PMFd,"%6d ",i + 1);

                for (j = 0; j <= SeqDNH; ++j) {
                    if (SeqDT[j]) {

                        lcen = rcen = 0;
                        ta = tb = -1;
                        for (k = min; k <= max; ++k) {
                            n = seq_sget(i,k,j);
                            if (n >= 0) {
                                ta = k;
                                sa = n;
                                if (k == min)
                                    lcen = 1;
                                break;
                            }
                        }
                        for (k = max; k >= min; --k) {
                            n = seq_sget(i,k,j);
                            if (n >= 0) {
                                tb = k;
                                sb = n;
                                if (k == max)
                                    rcen = 1;
                                break;
                            }
                        }
                        len = gmax = n1 = ng = ns = 0;
                        gmin = tb - ta + 10;
                        for (k = ta; k <= tb; ++k) {
                            n = seq_sget(i,k,j);
                            if (n >= 0) {
                                ns++;   
                                if (len > 0) {
                                    gmin = imin(len,gmin);
                                    gmax = imax(len,gmax);
                                    len = 0;
                                }
                            }
                            else {
                                if (n1 >= 0) {
                                    ng++;
                                    len = 1;
                                }   
                                else
                                    len++;
                            }
                            n1 = n;
                        }
                        if (ng == 0)
                            gmin = gmax = 0;
    
                        fprintf(PMFd,"%5d %5d %2d %2d %1d %1d %5d %5d %5d %5d  ",
                            ta,tb,sa,sb,lcen,rcen,ns,ng,gmin,gmax);
                                     
                    }
                }
                for (k = 0; k < PMNV; ++k) {
                    iv = (int)PMVIdx[k];
                    fprintf(PMFd,VPFmtS[iv],get_data(iv,i));
                }
                fprintf(PMFd,"\n");
                m++;
            }
            break;
       
        case 6:
        case 7:

            for (j = 0; j <= SeqDNH; ++j) {
                if (SeqDT[j])  
                    break;
            }

            ns = 0;
            for (i = 0; i < NOC; ++i) {

                if (eval_sve(i) == 0)   
                    continue;

                ta = min;
                tb = max;
                for (k = min; k <= max; ++k) {
                    n = seq_sget(i,k,j);
                    if (n >= 0) {
                        ta = k;
                        break;
                    }
                }
                for (k = max; k >= min; --k) {
                    n = seq_sget(i,k,j);
                    if (n >= 0) {
                        tb = k;
                        break;
                    }
                }
                n1 = -1;
                ng = 0;
                for (k = ta; k <= tb; ++k) {
                    n = seq_sget(i,k,j);
                    if (PMM == 7 && n < 0)
                        continue;

                    if (n != n1)
                        ng++;
                    n1 = n;
                }
                ns = imax(ns,ng);
            }
            for (i = 0; i < NOC; ++i) {

                if (eval_sve(i) == 0)   
                    continue;

                fprintf(PMFd,"%6d ",i + 1);

                ta = min;
                tb = max;
                for (k = min; k <= max; ++k) {
                    n = seq_sget(i,k,j);
                    if (n >= 0) {
                        ta = k;
                        break;
                    }
                }
                for (k = max; k >= min; --k) {
                    n = seq_sget(i,k,j);
                    if (n >= 0) {
                        tb = k;
                        break;
                    }
                }
                n1 = -1;
                ng = 0;
                for (k = ta; k <= tb; ++k) {
                    n = seq_sget(i,k,j);
                    if (PMM == 7 && n < 0)
                        continue;

                    if (n != n1) {
                        ng++;
                        fprintf(PMFd,"%2d ",n);
                    }
                    n1 = n;
                }
                for (k = ng + 1; k <= ns; ++k)
                    fprintf(PMFd,"%2d ",-9);

                fprintf(PMFd,"%5d ",ng);
                  
                for (k = 0; k < PMNV; ++k) {
                    iv = (int)PMVIdx[k];
                    fprintf(PMFd,VPFmtS[iv],get_data(iv,i));
                }
                fprintf(PMFd,"\n");
                m++;
            }
            break;

        case 8:                        
            if (PMS < 1)
                PMS = 1;
            printf1("Aggregation with s = %d\n",PMS);

            for (j = 0; j <= SeqDNH; ++j) {         /* use first sequence */
                if (SeqDT[j])  
                    break;
            }
            ns = SeqSTN[j];        /* number of states */
            if (alloc_acn(ns + 1))
                goto SEQPDFin;

            min = SeqTMin[j];
            max = SeqTMax[j];

            for (i = 0; i < NOC; ++i) {

                if (eval_sve(i) == 0)   
                    continue;

                fprintf(PMFd,"%6d ",i + 1);

                k = min;
                while (k <= max) {

                    for (l = 0; l < ns; ++l)
                        AcN[l] = 0;
                    n1 = 0;
                    for (l = k; l < k + PMS; ++l) {
                        if (l > max)
                            break;
                        n = seq_sget(i,l,j);

                        if (n >= 0) {
                            iv = SeqSTNII[j][n];
                            AcN[iv] += 1;
                        }
                        else
                            n1++;                   /* duration in missing state */
                    }
                    n = -1;  
                    for (iv = 0; iv < ns; ++iv) {
                        if (AcN[iv] > n1) {
                            n1 = AcN[iv];
                            n  = SeqSTNI[j][iv];
                        }   
                    }
                    fprintf(PMFd,"%2d ",n);
                    k += PMS;
                }
                for (k = 0; k < PMNV; ++k) {
                    iv = (int)PMVIdx[k];
                    fprintf(PMFd,VPFmtS[iv],get_data(iv,i));
                }
                fprintf(PMFd,"\n");
                    m++;
            }
            break;

        case 9:                                               
            for (i = 0; i < NOC; ++i) {

                if (eval_sve(i) == 0)   
                    continue;

                fprintf(PMFd,"%6d ",i + 1);

                for (j = 0; j <= SeqDNH; ++j) {
                    if (SeqDT[j]) {
                        /*********************************
                        for (k = 0; k < SeqSTN[j]; ++k) {
                            n = SeqSTNI[j][k];
                            printf1(" %d",n);
                        }
                        printf1("\n");
                        ********************************/

                        n1 = SeqSTN[j];  
                        if (alloc_acn(n1 + 1))
                            goto SEQPDFin;

                        for (k = min; k <= max; ++k) {
                            n = seq_sget(i,k,j);
                            if (n >= 0)  
                                AcN[SeqSTNII[j][n]] += 1;
                            else
                                AcN[n1] += 1;
                        }   
                        for (k = 0; k < n1; ++k)
                            fprintf(PMFd,"%6d ",AcN[k]);
                        fprintf(PMFd,"%6d ",AcN[n1]);
                    }
                }
                for (k = 0; k < PMNV; ++k) {
                    iv = (int)PMVIdx[k];
                    fprintf(PMFd,VPFmtS[iv],get_data(iv,i));
                }
                fprintf(PMFd,"\n");
                m++;
            }
            break;

        case 10:                                     /* ### */
            for (i = 0; i < NOC; ++i) {

                if (eval_sve(i) == 0)   
                    continue;

                fprintf(PMFd,"%6d ",i + 1);

                for (j = 0; j <= SeqDNH; ++j) {
                    if (SeqDT[j]) {
                        /*********************************
                        for (k = 0; k < SeqSTN[j]; ++k) {
                            n = SeqSTNI[j][k];
                            printf1(" %d",n);
                        }
                        printf1("\n");
                        ********************************/

                        n1 = SeqSTN[j];  
                        if (alloc_acn(n1 + 1))
                            goto SEQPDFin;
                        if (alloc_acm(n1 + 1))
                            goto SEQPDFin;

                        for (k = 0; k <= n1; ++k)
                            AcN[k] = AcM[k] = -5;

                        for (k = min; k <= max; ++k) {
                            n = seq_sget(i,k,j);
                            if (n >= 0) {
                                n = SeqSTNII[j][n];
                                if (AcN[n] < 0)
                                    AcN[n] = k;
                                if (AcM[n] < k)
                                    AcM[n] = k;
                            }
                            else {
                                if (AcN[n1] < 0)
                                    AcN[n1] = k;
                                if (AcM[n1] < k)
                                    AcM[n1] = k;
                            }
                        }   
                        for (k = 0; k <= n1; ++k)
                            fprintf(PMFd,"%6d %6d  ",AcN[k],AcM[k]);
                    }
                }
                for (k = 0; k < PMNV; ++k) {
                    iv = (int)PMVIdx[k];
                    fprintf(PMFd,VPFmtS[iv],get_data(iv,i));
                }
                fprintf(PMFd,"\n");
                m++;
            }
            break;
        
        default:
            printf1("Error: option m=%d not available.\n",PMM);
            goto SEQPDFin;
    }
    prn_ssel(m);
    printf1("%d records written to output file: %s\n",m,PMFdName);

    if (PMTDAFDef && m > 0) {                    /* write TDA description */

        fprintf(PMTDAFd,"nvar(\n");
        fprintf(PMTDAFd,"  dfile = %s,\n",PMFdName);
        fprintf(PMTDAFd,"  noc = %d,\n",m);

        fprintf(PMTDAFd,"  CASE <5>[6.0] = c1 , # case number\n");
        n = 2;

        switch (PMM) {

            case 1:
                for (j = 0; j <= SeqDNH; ++j) {
                    if (SeqDT[j]) {
                        for (k = min; k <= max; ++k)  
                            fprintf(PMTDAFd,"  Y%d_%d <2>[2.0] = c%-2d, # state in seq %d, time %d\n",j,k,n++,j + 1,k);
                    }   
                }
                break;

            case  2:                               
                fprintf(PMTDAFd,"  SN   <2>[2.0] = c%-2d, # sequence number\n",n++);
                for (k = min; k <= max; ++k)  
                    fprintf(PMTDAFd,"  Y%-3d <2>[2.0] = c%-2d, # state at time %d\n",k,n++,k);
                break;

            case  3:                               
                fprintf(PMTDAFd,"  TIME <5>[4.0] = c%-2d, # time\n",n++);
                for (j = 0; j <= SeqDNH; ++j) {
                    if (SeqDT[j])  
                        fprintf(PMTDAFd,"  Y%-3d <2>[2.0] = c%-2d, # state in seq %d\n",j,n++,j + 1);
                }
                break;

            case  4:                               
                fprintf(PMTDAFd,"  SN   <2>[2.0] = c%-2d, # spell number\n",n++);
                fprintf(PMTDAFd,"  ORG  <5>[2.0] = c%-2d, # origin state\n",n++);
                fprintf(PMTDAFd,"  DES  <5>[2.0] = c%-2d, # destination state\n",n++);
                fprintf(PMTDAFd,"  TS   <5>[4.0] = c%-2d, # starting time\n",n++);
                fprintf(PMTDAFd,"  TF   <5>[4.0] = c%-2d, # ending time\n",n++);
                break;

            case 5:
                for (j = 0; j <= SeqDNH; ++j) {
                    if (SeqDT[j]) {
                        fprintf(PMTDAFd,"  TA%-2d   <5>[5.0] = c%-2d, # seq %d, first time point\n",j,n++,j + 1);
                        fprintf(PMTDAFd,"  TB%-2d   <5>[5.0] = c%-2d, # seq %d, last time point\n",j,n++,j + 1);
                        fprintf(PMTDAFd,"  SA%-2d   <2>[2.0] = c%-2d, # seq %d, state at TA\n",j,n++,j + 1);
                        fprintf(PMTDAFd,"  SB%-2d   <2>[2.0] = c%-2d, # seq %d, state at TB\n",j,n++,j + 1);
                        fprintf(PMTDAFd,"  LCen%-2d <2>[1.0] = c%-2d, # seq %d, 1 if left censored\n",j,n++,j + 1);
                        fprintf(PMTDAFd,"  RCen%-2d <2>[1.0] = c%-2d, # seq %d, 1 if right censored\n",j,n++,j + 1);
                        fprintf(PMTDAFd,"  NS%-2d   <5>[5.0] = c%-2d, # seq %d, number of valid states\n",j,n++,j + 1);
                        fprintf(PMTDAFd,"  NP%-2d   <5>[5.0] = c%-2d, # seq %d, number of parts\n",j,n++,j + 1);
                        fprintf(PMTDAFd,"  GLMin%-2d<5>[5.0] = c%-2d, # seq %d, minimal gap length\n",j,n++,j + 1);
                        fprintf(PMTDAFd,"  GLMax%-2d<5>[5.0] = c%-2d, # seq %d, maximal gap length\n",j,n++,j + 1);
                    }   
                }
                break;

            case  6:                               
            case  7:                               
                for (j = 1; j <= ns; ++j)  
                    fprintf(PMTDAFd,"  Y%-3d<2>[2.0] = c%-2d, # state\n",j,n++);
                fprintf(PMTDAFd,"  NS  <5>[2.0] = c%-2d, # number of states\n",n++);
                break;

            case 8:
                for (j = 0; j <= SeqDNH; ++j) {
                    if (SeqDT[j])  
                        break;
                }
                k = min;
                l = 0;
                while (k <= max) {
                    fprintf(PMTDAFd,"  Y%d <2>[2.0] = c%-2d, # state in time %d - %d\n",l++,n++,k,imin(max,k + PMS - 1));
                    k += PMS;
                }
                break;

            case 9:
                for (j = 0; j <= SeqDNH; ++j) {
                    if (SeqDT[j]) {
                        n1 = SeqSTN[j];  
                        for (k = 0; k < n1; ++k) {
                            ns = SeqSTNI[j][k];
                            fprintf(PMTDAFd,"  D%d_%d <5>[6.0] = c%-2d, # duration in seqence %d, state %d\n",j,ns,n++,j + 1,ns);
                        }
                        fprintf(PMTDAFd,"  D%d_m <5>[6.0] = c%-2d, # duration in seqence %d, state missing\n",j,n++,j + 1);
                    }   
                }
                break;

            case 10:
                for (j = 0; j <= SeqDNH; ++j) {
                    if (SeqDT[j]) {
                        n1 = SeqSTN[j];  
                        for (k = 0; k < n1; ++k) {
                            ns = SeqSTNI[j][k];
                            fprintf(PMTDAFd,"  F%d_%d <5>[6.0] = c%-2d, # first time point in seqence %d, state %d\n",j,ns,n++,j + 1,ns);
                            fprintf(PMTDAFd,"  L%d_%d <5>[6.0] = c%-2d, #  last time point in seqence %d, state %d\n",j,ns,n++,j + 1,ns);
                        }
                        fprintf(PMTDAFd,"  F%d_m <5>[6.0] = c%-2d, # first time point in seqence %d, state missing\n",j,n++,j + 1);
                        fprintf(PMTDAFd,"  L%d_m <5>[6.0] = c%-2d, #  last time point in seqence %d, state missing\n",j,n++,j + 1);
                    }   
                }
                break;

            default:
                break;
        }   

        for (i = 0; i < PMNV; ++i) {
            j = PMVIdx[i];
            fprintf(PMTDAFd,"  %s <%d>[%d.%d]",VName[j],VSLen[j],VPFmt1[j],VPFmt2[j]);
            if (VLabel[j] != NULL)                                                      
                fprintf(PMTDAFd,"(%s)",VLabel[j]);                                           
            fprintf(PMTDAFd," = c%-2d,\n",n++);
        }
        fprintf(PMTDAFd,");\n");
    
        printf1("TDA description written to: %s\n",PMTDAFName);
    }
    err = 0;

SEQPDFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqsn()         States as new objects.                                  */
/*                                                                          */
/*                  seqgc(                                                  */
/*                      m=...,      option, def. 1                          */
/*                      sn=...,     sequence number, def. 1                 */
/*                      dtda=...,   TDA description file                    */  
/*                      nfmt=...,   integer print format, def. 4            */
/*                      fmt=...,    floating point print format, def. 10.4  */
/*                  ) = fname;      output file                             */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int seqsn(void)
{
    register int i,j;
    int err,nrec,n,ni,nj,ns,ns1,sn,a,b,t,s,na,nb,nflag;
    double d;

    err = -1;         
    if (check_cmd(1))
        return(-1);
         
    printf1("States as new objects.\n");
         
    if (parm(CmdBuf + 5,1,1))    /* get parameters */
        goto SEQSNFin;
        
    if (PMFmt1 == 0)
        pmfmt(10,4);            /* default print format */

    sn = seq_getsn(PMSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQSNFin;

    ns = SeqSTN[sn];        /* number of states */

    if (alloc_acn(ns + 1))
        goto SEQSNFin;

    a = SeqTMin[sn];
    b = SeqTMax[sn];

    nrec = 0;

    switch (PMM) {

        case 1:    
        case 2:    
            for (t = a; t <= b; ++t) {

                na = nb = 0;
                for (j = 0; j < ns; ++j)
                    AcN[j] = 0;

                for (i = 0; i < NOC; ++i) {

                    n = seq_sget(i,t,sn);
                    if (n < 0)
                        continue;

                    ni = SeqSTNII[sn][n];
                    if (PMM == 1)
                        AcN[ni] += 1;
                    else if (PMM == 2)
                        AcN[ni] = i;      
    
                    if (t > a) {
                        if (seq_sget(i,t - 1,sn) < 0)
                            na++;
                    }
                    if (t < b) {
                        if (seq_sget(i,t + 1,sn) < 0)
                            nb++;
                    }
                }
                fprintf(PMFd,PMNFmtS,t);         /* time */
                n = 0;
                for (j = 0; j < ns; ++j) {
                    fprintf(PMFd,PMNFmtS,AcN[j]);  
                    if (PMM == 1)
                        n += AcN[j];
                    else if (AcN[j] > 0)
                        n++;
                }        
                fprintf(PMFd,PMNFmtS,n);  
                fprintf(PMFd,PMNFmtS,na);  
                fprintf(PMFd,PMNFmtS,nb);  
                fprintf(PMFd,"\n");
                nrec++;
            }
            break;

        case 3:                                 /* permutations */
            for (t = a + 1; t <= b; ++t) {

                for (j = 0; j < ns; ++j)
                    AcN[j] = -1;

                for (i = 0; i < NOC; ++i) {

                    na = seq_sget(i,t - 1,sn);
                    if (na < 0)
                        continue;

                    nb = seq_sget(i,t,sn);
                    if (nb < 0)
                        continue;

                    ni = SeqSTNII[sn][na];
                    nj = SeqSTNII[sn][nb];
                    AcN[ni] = nj;
                }
                fprintf(PMFd,PMNFmtS,t);         /* time */
                for (j = 0; j < ns; ++j)   
                    fprintf(PMFd,PMNFmtS,AcN[j] + 1);  
                fprintf(PMFd,"\n");
                nrec++;
            }
            break;

        case 4:    

            if (alloc_aci(ns + 1))
                goto SEQSNFin;
            if (alloc_acj(ns + 1))
                goto SEQSNFin;
            if (alloc_ack(ns + 1))
                goto SEQSNFin;
            if (alloc_acr(ns + 1))
                goto SEQSNFin;
            if (alloc_acs(ns + 1))
                goto SEQSNFin;

            for (i = 0; i < NOC; ++i) {

                for (j = 0; j < ns; ++j)
                    AcI[j] = 0;

                nflag = 1;
                for (t = a; t <= b; ++t) {
                    n = seq_sget(i,t,sn);
                    if (n < 0) {
                        nflag = 1;
                        continue;
                    }
                    ni = SeqSTNII[sn][n];
                    AcI[ni] = 1;

                    if (t > a && n != seq_sget(i,t - 1,sn))
                        AcR[ni] += 1;
                    if (t < b && n != seq_sget(i,t + 1,sn))
                        AcS[ni] += 1;

                    AcJ[ni] += 1;
                    if (nflag || (t > a && n != seq_sget(i,t - 1,sn)))
                        AcK[ni] += 1;
                    nflag = 0;
                }
                for (j = 0; j < ns; ++j)
                    AcN[j] += AcI[j];
            }

            for (j = 0; j < ns; ++j) {
                fprintf(PMFd,PMNFmtS,j + 1);  
                fprintf(PMFd,PMNFmtS,SeqSTNI[sn][j]);  
                fprintf(PMFd,PMNFmtS,AcN[j]);  
                fprintf(PMFd,PMNFmtS,AcR[j]);  
                fprintf(PMFd,PMNFmtS,AcS[j]);  
                fprintf(PMFd,PMNFmtS,AcJ[j]);  
                fprintf(PMFd,PMNFmtS,AcK[j]);  

                if (AcK[j] > 0)
                    d = (double)AcJ[j] / (double)AcK[j];
                else
                    d = 0.0;
                fprintf(PMFd,PMFmtS,d);  
                fprintf(PMFd,"\n");
                nrec++;
            }
            break;

        case 5:    

            ns1 = ns + 1;
            if (alloc_acn(ns1 * ns1 + 1))
                goto SEQSNFin;

            for (i = 0; i < NOC; ++i) {
                for (t = a + 1; t <= b; ++t) {
                    nb = seq_sget(i,t,sn);
                    na = seq_sget(i,t - 1,sn);
                    if (na < 0 && nb < 0)
                        continue;

                    if (na < 0) {
                        ni = SeqSTNII[sn][nb];
                        AcN[ns * ns1 + ni] += 1;
                    }
                    else if (nb < 0) {
                        ni = SeqSTNII[sn][na];
                        AcN[ni * ns1 + ns] += 1;
                    }
                    else {
                        ni = SeqSTNII[sn][na];
                        nj = SeqSTNII[sn][nb];
                        AcN[ni * ns1 + nj] += 1;
                    }
                }
            }
            for (i = 0; i <= ns; ++i) {

                fprintf(PMFd,PMNFmtS,i + 1);    

                if (i < ns)
                    n = SeqSTNI[sn][i];
                else
                    n = -1;

                fprintf(PMFd,PMNFmtS,n);    

                for (j = 0; j <= ns; ++j)  
                    fprintf(PMFd,PMNFmtS,AcN[i * ns1 + j]);  

                fprintf(PMFd,"\n");
                nrec++;
            }
            break;

        default:
            printf1("Error: option m=%d not available.\n",PMM);
            goto SEQSNFin;
    }
    printf1("%d records written to: %s\n",nrec,PMFdName);

    if (PMTDAFDef && nrec > 0) {                    /* write TDA description */

        fprintf(PMTDAFd,"nvar(\n");
        fprintf(PMTDAFd,"  dfile = %s,\n",PMFdName);
        fprintf(PMTDAFd,"  noc = %d,\n",nrec);

        n = 1;

        switch (PMM) {

            case 1:
            case 2:
                fprintf(PMTDAFd,"  T   <5>[%d.0] = c%-2d, # time\n",PMNFmt,n++);
                for (j = 0; j < ns; ++j) {
                    s = SeqSTNI[sn][j];
                    fprintf(PMTDAFd,"  N%-3d<5>[%d.0] = c%-2d, # number of sequences (or ID) in state %d\n",s,PMNFmt,n++,s);
                }
                fprintf(PMTDAFd,"  N   <5>[%d.0] = c%-2d, # total number of sequences\n",PMNFmt,n++);
                fprintf(PMTDAFd,"  NA  <5>[%d.0] = c%-2d, # number of new sequences\n",PMNFmt,n++);
                fprintf(PMTDAFd,"  NB  <5>[%d.0] = c%-2d, # number of dropped sequences\n",PMNFmt,n++);
                break;

            case 4:
                fprintf(PMTDAFd,"  RN  <5>[%d.0] = c%d, # record number\n",PMNFmt,n++);
                fprintf(PMTDAFd,"  S   <5>[%d.0] = c%d, # state\n",PMNFmt,n++);
                fprintf(PMTDAFd,"  NS  <5>[%d.0] = c%d, # number of sequences in state\n",PMNFmt,n++);
                fprintf(PMTDAFd,"  NA  <5>[%d.0] = c%d, # number of sequences entering state\n",PMNFmt,n++);
                fprintf(PMTDAFd,"  NB  <5>[%d.0] = c%d, # number of sequences leaving state\n",PMNFmt,n++);
                fprintf(PMTDAFd,"  TT  <5>[%d.0] = c%d, # total time in state\n",PMNFmt,n++);
                fprintf(PMTDAFd,"  NE  <5>[%d.0] = c%d, # number of episodes in state\n",PMNFmt,n++);
                fprintf(PMTDAFd,"  DUR <4>[%d.%d] = c%d, # mean duration (TT/NE)\n",PMFmt1,PMFmt2,n++);
                break;

            case 5:
                fprintf(PMTDAFd,"  RN  <5>[%d.0] = c%d, # record number\n",PMNFmt,n++);
                fprintf(PMTDAFd,"  S   <5>[%d.0] = c%d, # state\n",PMNFmt,n++);
                for (j = 1; j <= ns1; ++j)
                    fprintf(PMTDAFd,"  N%-3d<5>[%d.0] = c%d,\n",j,PMNFmt,n++);
                break;

            default:
                break;
        }   
        fprintf(PMTDAFd,");\n");
        printf1("TDA description written to: %s\n",PMTDAFName);
    }
    err = 0;

SEQSNFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqgm()         Indicators of group membership.                         */
/*                                                                          */
/*                  seqgm(                                                  */
/*                      sn=...,     sequence number, def. 1                 */
/*                      s=...,      state number, def. 0                    */
/*                      sk=...,     sequence of max 10 length values        */
/*                      dtda=...,   TDA description file                    */  
/*                      nfmt=...,   integer print format, def. 4            */
/*                  ) = fname;      output file                             */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int seqgm(void)
{
    register int i,j;
    int err,nrec,n,sn,a,b,t,na,nb,na1,na2,d;
    int nn[11];
    double dd;

    err = -1;         
    if (check_cmd(1))
        return(-1);
         
    printf1("Indicators of group membership.\n");
         
    if (parm(CmdBuf + 5,1,1))    /* get parameters */
        goto SEQMLFin;
        
    if (PMFmt1 == 0)
        pmfmt(10,4);            /* default print format */

    sn = seq_getsn(PMSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQMLFin;

    a = SeqTMin[sn];
    b = SeqTMax[sn];

    fprintf(PMFd,"# Time Number_in_state_%d",PMS);
    for (j = 1; j <= 10; ++j) {
        if (PMSK[j] > 0) 
            fprintf(PMFd," Len=%d ",PMSK[j]);  
    }
    fprintf(PMFd,"\n");

    nrec = 0;
    for (t = a; t <= b; ++t) {

        d = na = na1 = na2 = 0;
        for (j = 1; j <= 10; ++j)
            nn[j] = 0;

        for (i = 0; i < NOC; ++i) {

            n = seq_sget(i,t,sn);
            if (n != PMS)
                continue;
            na++;
            if (t == a)
                na1++;
            else if (seq_sget(i,t - 1,sn) != PMS)
                na1++;

            if (t == b)
                na2++;
            else if (seq_sget(i,t + 1,sn) != PMS)
                na2++;

            d++;
            for (j = t - 1; j >= a; --j) {
                if (seq_sget(i,j,sn) != PMS) 
                    break;
                d++;
            }
            nb = 1;   
            for (j = t - 1; j >= a; --j) {
                if (seq_sget(i,j,sn) != PMS)
                    break;
                nb++;
            }
            for (j = 1; j <= 10; ++j) {
                if (nb >= PMSK[j])
                    nn[j] += 1;
            }
        }
        fprintf(PMFd,PMNFmtS,t);         /* time */
        fprintf(PMFd,PMNFmtS,na);  
        fprintf(PMFd,PMNFmtS,na1);  
        fprintf(PMFd,PMNFmtS,na2);  

        dd = 0.0;
        if (na > 0)
            dd = (double)d / (double)na;

        fprintf(PMFd,PMFmtS,dd);  
        for (j = 1; j <= 10; ++j) {
            if (PMSK[j] > 0) 
                fprintf(PMFd,PMNFmtS,nn[j]);  
        }
        fprintf(PMFd,"\n");
        nrec++;
    }
    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

SEQMLFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqgsp      Generate data for sequence plots.                           */
/*                                                                          */
/*              seqggsp(                                                    */
/*                      sn=...,     sequence number, def. 1                 */
/*                      s=...,      state number, def. 0                    */
/*                      nfmt=...,   integer print format, def. 4            */
/*                  ) = fname;      output file                             */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int seqgsp(void)
{
    register int i,j;
    int err,nrec,n,sn,a,b,t,f,l;   

    err = -1;         
    if (check_cmd(1))
        return(-1);
         
    printf1("Generation of data for sequence plots.\n");
         
    if (parm(CmdBuf + 6,1,1))    /* get parameters */
        goto SEQGSPFin;
        
    sn = seq_getsn(PMSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQGSPFin;

    a = SeqTMin[sn];
    b = SeqTMax[sn];

    if (alloc_acn(NOC + 1))
        goto SEQGSPFin; 
    if (alloc_acm(NOC + 1))
        goto SEQGSPFin; 
    if (alloc_ack(NOC + 1))
        goto SEQGSPFin; 

    for (i = 0; i < NOC; ++i) {
        f = -1;
        l = 0;
        for (t = a; t <= b; ++t) {
            n = seq_sget(i,t,sn);
            if (n != PMS)
                continue;
            if (f == -1)
                f = t;
            l++;
        }
        AcN[i] = f;
        AcM[i] = l;
    }
    if (sortdpi2(NOC,AcN,AcM,AcK)) 
        goto SEQGSPFin;

    /******************* 
    printf("AcN: ");
    for (i = 0; i < NOC; ++i)
        printf(" %d ",AcN[i]);

    printf("\n AcM: ");
    for (i = 0; i < NOC; ++i)
        printf(" %d ",AcM[i]);
    printf("\n AcK: ");
    for (i = 0; i < NOC; ++i)
        printf(" %d ",AcK[i]);
    newline();
    *******************/

    nrec = 0;
    for (t = a; t <= b; ++t) {
        fprintf(PMFd,PMNFmtS,t);         /* time */
                 
        for (i = 0; i < NOC; ++i) {
            j = AcK[i];
            n = seq_sget(j,t,sn);
            if (n != PMS)
                n = 0;
            else
                n = NOC - i; 

            fprintf(PMFd,PMNFmtS,n);
        }
        fprintf(PMFd,"\n");
        nrec++;
    }
    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

SEQGSPFin:
    p_clean();
    return(err);
}

