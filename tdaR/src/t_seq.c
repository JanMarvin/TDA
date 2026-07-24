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
#include "tda_context.h"
#include "t_sort.h"

/*  functions in t_seq.c */

int t_seq(TDAContext *ctx);
int check_sdj(TDAContext *ctx, int j);
int get_sdata(TDAContext *ctx, int i,int j,int s);
int seq_sget(TDAContext *ctx, int i,int t,int n);
int seq_tsget(TDAContext *ctx, int i,int n);
int seq_tfget(TDAContext *ctx, int i,int n);
int seq_getsn(TDAContext *ctx, int sn,int opt);
int seq_scheck(TDAContext *ctx, int i,int j);
void seq_free(TDAContext *ctx, int i);                   
void seq_afree(TDAContext *ctx, int opt);                   
int seqdel(TDAContext *ctx);                   
int seqdef(TDAContext *ctx);
int seq_rc(TDAContext *ctx, int s,int smax);
void seqadj(TDAContext *ctx);
void prnsd(TDAContext *ctx);
int prnsdi(TDAContext *ctx);                   
int seqlg(TDAContext *ctx);
int seqld(TDAContext *ctx);
void prn_ssel(TDAContext *ctx, int n);
int seq_csel(TDAContext *ctx);
void seq_dtda(TDAContext *ctx, int typ,int sn,char *fname,int noc,int nv,short *vidx);
int seqgc(TDAContext *ctx);
int seqsd(TDAContext *ctx);
int seqen(TDAContext *ctx);
int seqtp(TDAContext *ctx);
int seqrd(TDAContext *ctx);
int seqsi(TDAContext *ctx);
int seqpe(TDAContext *ctx);
void seqpe1(TDAContext *ctx, int i,int id,int *st);
int seqpd(TDAContext *ctx);
int seqsn(TDAContext *ctx);
int seqgm(TDAContext *ctx);
int seqgsp(TDAContext *ctx);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */





/*--------------------------------------------------------------------------*/
/*  t_seq()     Entry point for sequence commands.                          */
/*              Return 0 if OK, -1 if error, 1 if command cannot be         */
/*              interpreted.                                                */

int t_seq(TDAContext *ctx)
{
    char *p;
    int err = 0;

    p = ctx->CmdBuf;

    if (!strncmp(p,"seqdef",6))         /* new sequence data */
        err = seqdef(ctx);                 
    else if (!strncmp(p,"seqdel",6))    /* delete sequence data */
        err = seqdel(ctx);                 
    else if (!strncmp(p,"seqlg",5))     /* sequence length */
        err = seqlg(ctx);                 
    else if (!strncmp(p,"seqld",5))     /* last occurrence */
        err = seqld(ctx);                 
    else if (!strncmp(p,"seqgc",5))     /* sequence characteristics */
        err = seqgc(ctx);                 
    else if (!strncmp(p,"seqsd",5))     /* state distributions */
        err = seqsd(ctx);                 
    else if (!strncmp(p,"seqen",5))     /* entropy */
        err = seqen(ctx);                 
    else if (!strncmp(p,"seqtp",5))     /* transition probabilities */
        err = seqtp(ctx);                 
    else if (!strncmp(p,"seqrd",5))     /* random sequences */
        err = seqrd(ctx);                 
    else if (!strncmp(p,"seqsi",5))     /* state indicator matrix */
        err = seqsi(ctx);                 
    else if (!strncmp(p,"seqpe",5))     /* create sequence from episode data */
        err = seqpe(ctx);                 
    else if (!strncmp(p,"seqpd",5))     /* print sequence data */
        err = seqpd(ctx);                 
    else if (!strncmp(p,"seqevd",6))    /* data file with events */
        err = seqevd(ctx);                 
    else if (!strncmp(p,"seqev",5))     /* info about events */
        err = seqev(ctx);                 
    else if (!strncmp(p,"seqmd",5))     /* data for event models */
        err = seqmd(ctx);                 
    else if (!strncmp(p,"seqgm",5))     /* group membership */
        err = seqgm(ctx);                 
    else if (!strncmp(p,"seqm",4))      /* optimal matching */
        err = seqm(ctx);                 
    else if (!strncmp(p,"seqpm",5))     /* pattern matching */
        err = seqpm(ctx);                 
    else if (!strncmp(p,"seqsn",5))     /* states as new objects */
        err = seqsn(ctx);                 
    else if (!strncmp(p,"seqgsp",6))    /* generate plot data */
        err = seqgsp(ctx);                 

    else if (!strcmp(p,"seq"))          /* print info */
        err = prnsdi(ctx);                 
    else
        err = 1;

    return(err);
}

/* ------------------------------------------------------------------------ */
/*  check_sdj(j)    Check whether variable j is needed for sequence data.   */
/*                  Return 1 if needed, otherwise 0.                        */

int check_sdj(TDAContext *ctx, int j)
{
    register int i,k,n;

    if (ctx->SeqDN == 0)
        return(0);

    for (i = 0; i <= ctx->SeqDNH; ++i) {
        n = ctx->SeqDNV[i];
        for (k = 0; k < n; ++k) {
            if (j == ctx->SeqDV[i][k])
                return(1);
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  get_sdata(i,j,n)    get sequence data for variable i, case j, for       */
/*                      sequence n. Observe recoding of states.             */

int get_sdata(TDAContext *ctx, int i,int j,int n)
{
    register int s;

    s = (int)get_data(ctx, i,j);
    if (s >= -9 && ctx->SeqRCA[n] > 0)  
       s = ctx->SeqRC[n][s + 9];
    return(s);
}

/* ------------------------------------------------------------------------ */
/*  seq_sget(i,t,n)   Get state for case i at time t; -1 if not present.    */
/*                    n is number of sequence.                              */

int seq_sget(TDAContext *ctx, int i,int t,int n)
{
    register int j,s,nv;
    int ts,tf,ta,tb,sa,sb;

    if (n >= 0 && n <= ctx->SeqDNH) { 
        if (ctx->SeqDT[n] == 1) {
            if (t < 0 || t >= ctx->SeqDNV[n])
                return(-1);
            s = get_sdata(ctx, ctx->SeqDV[n][t],i,n);

            if (s < 0 && ctx->SeqGL[n] > 0) {
                ta = t - 1;
                sa = -1;
                while (ta >= 0) {
                    sa = get_sdata(ctx, ctx->SeqDV[n][ta],i,n);
                    if (sa >= 0)
                        break;
                    ta--;
                }   
                tb = t + 1;
                sb = -1;
                while (tb < ctx->SeqDNV[n]) {
                    sb = get_sdata(ctx, ctx->SeqDV[n][tb],i,n);
                    if (sb >= 0)
                        break;
                    tb++;
                }   
                if (sa >= 0 && sa == sb && tb - ta -1 <= ctx->SeqGL[n])
                    s = sa;
            }
            return(s);
        }
        else  {
            s = -1;
            ts = (int)get_data(ctx, ctx->SeqDV[n][1],i);
            if (t >= ts) {
                nv = ctx->SeqDNV[n];
                tf = (int)get_data(ctx, ctx->SeqDV[n][nv - 1],i);
                if (t == tf)
                    s = get_sdata(ctx, ctx->SeqDV[n][nv - 2],i,n);
                else if (t < tf) {
                    for (j = 3; j < nv; j += 2) {
                        if (t < (int)get_data(ctx, ctx->SeqDV[n][j],i)) {
                            s = get_sdata(ctx, ctx->SeqDV[n][j - 3],i,n);
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

int seq_tsget(TDAContext *ctx, int i,int n)
{
    register int t,typ;

    if (n >= 0 && n <= ctx->SeqDNH) { 
        typ = ctx->SeqDT[n];
        if (typ == 1) {
            for (t = 0; t < ctx->SeqDNV[n]; ++t) {
                if (get_sdata(ctx, ctx->SeqDV[n][t],i,n) >= 0)
                    return(t);
            }
        }
        else if (typ == 2) {
            t = (int)get_data(ctx, ctx->SeqDV[n][1],i);
            return(t);
        }
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  seq_tfget(i,n)  Get last time point in sequence n for individual i.     */

int seq_tfget(TDAContext *ctx, int i,int n)
{
    register int t,typ,nv;

    if (n >= 0 && n <= ctx->SeqDNH) { 
        typ = ctx->SeqDT[n];
        nv = ctx->SeqDNV[n];

        if (typ == 1) {
            for (t = nv - 1; t >= 0; --t) {
                if (get_sdata(ctx, ctx->SeqDV[n][t],i,n) >= 0)
                    return(t);
            }
        }
        else if (typ == 2) {
            t = (int)get_data(ctx, ctx->SeqDV[n][nv - 1],i);
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

int seq_getsn(TDAContext *ctx, int sn,int opt)
{
    register int k;
    sn--;
    for (k = 0; k <= ctx->SeqDNH; ++k) {
        if (ctx->SeqDT[k]) {
            if (sn < 0 || sn == k) {
                if (opt) {
                    printf1(ctx, "Using sequence data structure %d.\n",k + 1);
                }
                return(k);
            }
        }
    }
    if (opt && sn >= 0)  
        printf1(ctx, "Error: sequence data structure %d not defined.\n",sn + 1);
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  seq_scheck(i,j)     Check whether sequences i and j have the same       */
/*                      state space. Return 0 if so, otherwise -1.          */

int seq_scheck(TDAContext *ctx, int i,int j)
{
    register int k,n;

    n = ctx->SeqSTN[i];
    if (n != ctx->SeqSTN[j])
        return(-1);
    for (k = 0; k < n; ++k) {
        if (ctx->SeqSTNI[i][k] != ctx->SeqSTNI[j][k])
            return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  seq_free(i)  free sequence data structure i.                            */
   
void seq_free(TDAContext *ctx, int i)                    
{
    if (ctx->SeqDNV[i] > 0) {
        free((char *)ctx->SeqDV[i]);
        memrq(ctx, -ctx->SeqDNV[i],sizeof(short));
        ctx->SeqDNV[i] = 0;
    }   
    if (ctx->SeqSTNH[i] > 0) {
        free((char *)ctx->SeqSTNII[i]);
        memrq(ctx, -ctx->SeqSTNH[i],sizeof(short));
        ctx->SeqSTNH[i] = 0;
    }   
    if (ctx->SeqSTN[i] > 0) {
        free((char *)ctx->SeqSTNI[i]);
        memrq(ctx, -ctx->SeqSTN[i],sizeof(short));
        ctx->SeqSTN[i] = 0;  
    }   
    if (ctx->SeqRCA[i] > 0) {
        free((char *)ctx->SeqRC[i]);
        memrq(ctx, -ctx->SeqRCA[i],sizeof(short));
        ctx->SeqRCA[i] = 0;
    }   
    ctx->SeqGL[i] = ctx->SeqDT[i] = 0;
    seqadj(ctx);           /* adjust number of sequences */
} 

/* ------------------------------------------------------------------------ */
/*  seq_afree(opt)  free all sequence data structures.                      */
/*                  if opt != 0 print message.                              */
   
void seq_afree(TDAContext *ctx, int opt)                 
{
    int i,n;

    n = 0;
    for (i = 0; i < SEQ_Max; ++i) {
        if (ctx->SeqDT[i]) {
            seq_free(ctx, i);
            n++;
        }
    }
    if (n > 0 && opt) 
        printf1(ctx, "Deallocated %d sequence data structure(s).\n",n);
}

/* ------------------------------------------------------------------------ */
/*  seqdel()    delete sequence data.                                       */
/*              Return 0 if OK, -1 if error.                                */

int seqdel(TDAContext *ctx)                    
{
    int n;

    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->SeqDN == 0) {
        printf1(ctx, "No sequences defined.\n");
        return(0);
    }
    if (sscanf(ctx->CmdBuf,"seqdel=%d",&n) == 1) {
        if (n < 1 || n > SEQ_Max || ctx->SeqDT[n - 1] == 0)
            printf1(ctx, "Error: there is no sequence data structure %d (command ignored).\n",n);
        else {
            seq_free(ctx, n - 1);
            printf1(ctx, "Deallocated sequence data structure %d. Current memory: %d bytes.\n",n,ctx->MemReq);
        }
    }
    else if (!strcmp(ctx->CmdBuf,"seqdel")) {
        seq_afree(ctx, 1);
        printf1(ctx, "Current memory: %d bytes.\n",ctx->MemReq);
    }
    else {
        p_err(ctx, -1,1);
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

int seqdef(TDAContext *ctx)
{
    register int i,j;
    int err,s,smax,sn,t,tl,tmin,tmax;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Creating a new sequence data structure. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 6,4,1))     /* get parameters */
        goto SEQNDFin;

    if (ctx->PMSN < 1)
        ctx->PMSN = 1;
    else if (ctx->PMSN > SEQ_Max) {
        printf1(ctx, "Error: range of sequence data structures is: 1 - %d\n",SEQ_Max);
        goto SEQNDFin;
    }
    printf1(ctx, "Sequence structure number: %d\n",ctx->PMSN);

    ctx->PMSN--;     /* internally we use sequence numbers 0,1,... */

    printf1(ctx, "Sequence type: %d\n",ctx->PMM);
    if (ctx->PMM < 1 || ctx->PMM > 2) {
        printf1(ctx, "Error: only m = 1 or 2 is possible.\n");
        goto SEQNDFin;
    }
    if (ctx->PMM == 2) {
        sn = ctx->PMNV / 2;
        if (2 * sn != ctx->PMNV) {
            printf1(ctx, "Error: type 2 requires an even number of variables.\n");
            goto SEQNDFin;
        }
    }

    if (ctx->SeqDT[ctx->PMSN] != 0) {
        seq_free(ctx, ctx->PMSN);
        printf1(ctx, "Previously defined sequence data structure %d will be deleted.\n",ctx->PMSN + 1);
    }

    /* PMNV = number of new variables on right-hand side */

    ctx->SeqDT[ctx->PMSN] = (char)(ctx->PMM);

    /* if PMGLEN > 0 set SeqGL */

    if (ctx->PMGLEN > 0) {
        if (ctx->PMM != 1) {
            printf1(ctx, "Error: glen option only possible with m = 1.\n");
            goto SEQNDFin;
        }
        ctx->SeqGL[ctx->PMSN] = ctx->PMGLEN;
    }
    if (!(ctx->SeqDV[ctx->PMSN] = (short *)calloc((size_t)(ctx->PMNV),sizeof(short)))) { 
        p_err(ctx, -2,1);
        goto SEQNDFin;
    }
    memrq(ctx, ctx->PMNV,sizeof(short));
    ctx->SeqDNV[ctx->PMSN] = (short)(ctx->PMNV);

    for (i = 0; i < ctx->PMNV; ++i)
        ctx->SeqDV[ctx->PMSN][i] = ctx->PMVIdx[i];

    tmax = 0;
    tmin = ctx->INTMAX;
    smax = -1;

    for (i = 0; i < ctx->NOC; ++i) {     /* get highest state number */
        tl = 0;                     /* and tmin,tmax for type 2 data */
        j = 0;                      /* also check ranges */
        while (j < ctx->PMNV) {
            s = (int)get_data(ctx, ctx->PMVIdx[j],i);
            if (s >= 0 && smax < s)
                smax = s;

            j++;
            if (ctx->PMM == 2) {
                t = (int)get_data(ctx, ctx->PMVIdx[j],i);
                if ((j == 1 && t < 0) || (j > 1 && t <= tl)) {
                    printf1(ctx, "Error in case %d: need ascending event times.\n",i + 1);
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
        printf1(ctx, "Error: cannot find any valid states.\n");
        goto SEQNDFin;
    }
    if (alloc_acn(ctx, smax + 1))
        goto SEQNDFin; 

    for (i = 0; i < ctx->NOC; ++i) {         /* get state numbers */
        for (j = 0; j < ctx->PMNV; ++j) {
            s = (int)get_data(ctx, ctx->PMVIdx[j],i);
            if (s >= 0)  
                ctx->AcN[s] += 1;
            if (ctx->PMM == 2)
                j++;
        }
    }

    if (ctx->PRCAlloc) {                  /* get recode information */
        if (seq_rc(ctx, ctx->PMSN,smax))
            goto SEQNDFin;
        j = 0;
        for (i = 0; i <= smax + 9; ++i) {
            if (j < ctx->SeqRC[ctx->PMSN][i])
                j = ctx->SeqRC[ctx->PMSN][i];
        }
        smax = j;
        if (smax < 0) {
            printf1(ctx, "Error: no valid states after recoding.\n");
            goto SEQNDFin;
        }
        if (alloc_acn(ctx, smax + 1))
            goto SEQNDFin; 

        for (i = 0; i < ctx->NOC; ++i) {         /* get state numbers  */
            for (j = 0; j < ctx->PMNV; ++j) {
                s = (int)get_data(ctx, ctx->PMVIdx[j],i);
                if (s >= -9) {
                    s = ctx->SeqRC[ctx->PMSN][s + 9];
                        if (s >= 0) 
                              ctx->AcN[s] += 1;
                }
                if (ctx->PMM == 2)
                    j++;
            }
        }
    }
    if (!(ctx->SeqSTNII[ctx->PMSN] = (short *)calloc((size_t)(smax + 1),sizeof(short)))) { 
        p_err(ctx, -2,1);
        goto SEQNDFin;
    }
    memrq(ctx, smax + 1,sizeof(short));
    ctx->SeqSTNH[ctx->PMSN] = (short)(smax + 1);

    sn = 0;                         /* number of different states */
    for (j = 0; j <= smax; ++j) {
        ctx->SeqSTNII[ctx->PMSN][j] = -1;
        if (ctx->AcN[j] > 0)  
            sn++;
    }
    if (!(ctx->SeqSTNI[ctx->PMSN] = (short *)calloc((size_t)(sn),sizeof(short)))) { 
        p_err(ctx, -2,1);
        goto SEQNDFin;
    }
    memrq(ctx, sn,sizeof(short));
    ctx->SeqSTN[ctx->PMSN] = (short)(sn);

    i = 0;                          /* internal state number */
    for (j = 0; j <= smax; ++j) {
        if (ctx->AcN[j] > 0) {
            ctx->SeqSTNI[ctx->PMSN][i] = (short)(j);
            ctx->SeqSTNII[ctx->PMSN][j] = (short)(i);
            i++;
        }
    }
    /*************** 
    tda_out("SeqSTNI\n");
    for (i = 0; i < SeqSTN[PMSN]; ++i) {
        printf1(ctx, "i=%d stni[][]=%d\n",i,SeqSTNI[PMSN][i]);
    }
    tda_out("SeqSTNII\n");
    for (i = 0; i < SeqSTNH[PMSN]; ++i) {
        printf1(ctx, "i=%d stnii[][]=%d\n",i,SeqSTNII[PMSN][i]);
    }
    printf1(ctx, "SeqRCA[PMSN] = %d\n",SeqRCA[PMSN]);
    tda_out("SeqRC\n");
    for (i = 0; i < SeqRCA[PMSN]; ++i)
        printf1(ctx, "i=%d rc=%d\n",i,SeqRC[PMSN][i]);
    *************/
   
    if (ctx->PMM == 1) {                 /* set time axis */
        ctx->SeqTMin[ctx->PMSN] = 0;
        ctx->SeqTMax[ctx->PMSN] = ctx->PMNV - 1;
    }
    else if (ctx->PMM == 2) {
        ctx->SeqTMin[ctx->PMSN] = tmin;
        ctx->SeqTMax[ctx->PMSN] = tmax;
    }
    seqadj(ctx);           /* adjust number of sequences */
    prnsd(ctx);            /* print info */
    err = 0;

SEQNDFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  seq_rc(s,smax)      Get state recode information for sequence s.        */
/*                      smax is max number of states.                       */
/*                      Return 0 if OK, -1 if error.                        */

int seq_rc(TDAContext *ctx, int s,int smax)
{
    register int i;
    register char *p;
    int err,n,m;

    err = -2;
    if (!(ctx->SeqRC[s] = (short *)calloc((size_t)(smax + 10),sizeof(short)))) { 
        p_err(ctx, -2,1);
        goto SEQRCFin;
    }
    memrq(ctx, smax + 10,sizeof(short));
    ctx->SeqRCA[s] = (short)(smax + 10);

    n = -9;
    for (i = 0; i <= smax + 9; ++i)  
        ctx->SeqRC[s][i] = (short)(n++);             
       
    p = ctx->PRC;
    err = -1;
               
    while (*p) {
        if (sscanf(p,"%d",&n) != 1 || n < -9) {
              goto SEQRCFin;
        }
        p = skip_int(ctx, p);
        if (!*p || *p++ != '[' || !*p)
            goto SEQRCFin;
        while (*p) {
            if (sscanf(p,"%d",&m) != 1 || m < -9 || m > smax)
                goto SEQRCFin;

            ctx->SeqRC[s][m + 9] = (short)(n); 
            p = skip_int(ctx, p);
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
        printf1(ctx, "Syntax error or unknown states in recode option.\n");
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  seqadj()    Adjust number of sequences and highest sequence number.     */

void seqadj(TDAContext *ctx)
{
    register int i;

    ctx->SeqTCMin = ctx->INTMAX;
    ctx->SeqTCMax = ctx->SeqDNH = ctx->SeqDN = 0;
    for (i = 0; i < SEQ_Max; ++i) {
        if (ctx->SeqDT[i]) {
            ctx->SeqDN++;
            if (ctx->SeqDNH < i)
                ctx->SeqDNH = i;
            if (ctx->SeqTCMin > ctx->SeqTMin[i])
                ctx->SeqTCMin = ctx->SeqTMin[i];
            if (ctx->SeqTCMax < ctx->SeqTMax[i])
                ctx->SeqTCMax = ctx->SeqTMax[i];
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  prnsd()     Print info about new sequence data.                         */

void prnsd(TDAContext *ctx)
{
    register int i,j,k;
    int n,s,ss,first;

    if (ctx->SeqDN == 0) {
        printf1(ctx, "No sequences defined.\n");
        return;
    }
    printf1(ctx, "Currently defined sequences:\n\n");
    printf1(ctx, "Sequence          State       Time axis       Number\n");
    printf1(ctx, "Structure Type  Variables  Minimum  Maximum  of States  States\n");
    prnchar(ctx, '-',62,1);

    for (j = 0; j <= ctx->SeqDNH; ++j) {
        if (ctx->SeqDT[j]) {
            n = ctx->SeqDNV[j] / (int)ctx->SeqDT[j];
            printf1(ctx, "%6d    %3d   %7d  %9d %8d  %8d   ",
                    j + 1,ctx->SeqDT[j],n,ctx->SeqTMin[j],ctx->SeqTMax[j],ctx->SeqSTN[j]);
#ifdef TDA_R_PACKAGE
            {
                double erow[6];
                erow[0] = (double)(j + 1);
                erow[1] = (double)ctx->SeqDT[j];
                erow[2] = (double)n;
                erow[3] = (double)ctx->SeqTMin[j];
                erow[4] = (double)ctx->SeqTMax[j];
                erow[5] = (double)ctx->SeqSTN[j];
                tda_export_row(ctx, "seq.defs", erow, 6);
            }
#endif
                       
            for (i = 0; i < ctx->SeqSTN[j]; ++i) {
                s = ctx->SeqSTNI[j][i];
                printf1(ctx, " %d",s);
#ifdef TDA_R_PACKAGE
                tda_export_cell(ctx, "seq.states", (double)s);
#endif
                if (ctx->SeqRCA[j] > 0) {
                    first = 1;
                    for (k = 0; k < ctx->SeqRCA[j]; ++k) {
                        ss = k - 9;              
                        if (ctx->SeqRC[j][k] == s && s != ss) {
                            if (first) {
                                printf1(ctx, "[%d",ss);
                                first = 0;
                            }
                            else  
                                printf1(ctx, ",%d",ss);
                        }
                    }
                    if (first == 0)
                        printf1(ctx, "]");
                }
            }
            newline(ctx);
#ifdef TDA_R_PACKAGE
            /* one structure's states are one row, and structures can
               have different numbers of them, so each is flushed on its
               own: seq.states, seq.states.2, ... in structure order */
            tda_export_endrow(ctx, "seq.states");
            tda_export_flush(ctx, "seq.states");
#endif
        }
    }
    printf1(ctx, "\nRange of common time axis: %d to %d.\n",ctx->SeqTCMin,ctx->SeqTCMax);
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "seq.defs");
    {
        double ax[2];
        ax[0] = (double)ctx->SeqTCMin;
        ax[1] = (double)ctx->SeqTCMax;
        tda_export_mat(ctx, "seq.timeaxis", ax, 1, 2);
    }
#endif
}

/* ------------------------------------------------------------------------ */
/*  prnsdi()    Print info about currently defined sequences.               */
/*              Return 0 if OK, -1 if error.                                */

int prnsdi(TDAContext *ctx)                    
{
    if (check_cmd(ctx, 1))
        return(-1);

    prnsd(ctx);
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  seqlg()         Sequence length and gaps.                               */
/*                  seqlg(sn=,sel=,v=,dtda=) = fname                        */
/*                  Return 0 if OK, -1 if error.                            */

int seqlg(TDAContext *ctx)
{
    register int i,j,l;
    int err,n,sn,a,b,ai,bi,t,len,gl,nrec,ng,gmin,gmax;
#ifdef TDA_R_PACKAGE
    double *erow = NULL;
    int ecol = 0;
#endif

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);
    printf1(ctx, "Sequence length and gaps.\n");
         
    if (parm(ctx, ctx->CmdBuf + 5,1,1))    /* get parameters */
        goto SEQLFin;
        
    sn = seq_getsn(ctx, ctx->PMSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQLFin;

    a = ctx->SeqTMin[sn];
    b = ctx->SeqTMax[sn];

    fprintf(ctx->PMFd,"#       starting  ending  sequence     gap  min gap  max gap   number\n");
    fprintf(ctx->PMFd,"# Case      time    time    length  length   length   length  of gaps\n");

    nrec = 0;
    for (i = 0; i < ctx->NOC; ++i) {

        if (eval_sve(ctx, i) == 0)   
            continue;

        len = gl = 0;
        ai = bi = -1;

        ng = 0;         /* number of gaps */
        gmin = ctx->INTMAX;  /* min gap length */
        gmax = 0;       /* max gap length */

        for (t = a; t <= b; ++t) {
            n = seq_sget(ctx, i,t,sn);
            if (n >= 0) {
                ai = t;
                break;
            }
        }
        if (ai >= 0) {
            for (t = b; t >= ai; --t) {
                n = seq_sget(ctx, i,t,sn);
                if (n >= 0) {
                    bi = t;
                    break;
                }
            }
            len = bi - ai + 1;
            l = 0;
            for (t = ai + 1; t <= bi; ++t) {
                n = seq_sget(ctx, i,t,sn);
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
        if (gmin == ctx->INTMAX)
            gmin = 0;

        fprintf(ctx->PMFd,"%8d %7d %7d %9d %7d %8d %8d %8d ",
                                         i + 1,ai,bi,len,gl,gmin,gmax,ng);
#ifdef TDA_R_PACKAGE
        if (erow == NULL) {
            ecol = 8 + ctx->PMNV;
            erow = (double *)malloc((size_t)ecol * sizeof(double));
        }
        if (erow != NULL) {
            erow[0] = (double)(i + 1); erow[1] = (double)ai;
            erow[2] = (double)bi;      erow[3] = (double)len;
            erow[4] = (double)gl;      erow[5] = (double)gmin;
            erow[6] = (double)gmax;    erow[7] = (double)ng;
        }
#endif
        for (j = 0; j < ctx->PMNV; ++j) {
            n = (int)ctx->PMVIdx[j];
            rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[n],get_data(ctx, n,i));
#ifdef TDA_R_PACKAGE
            if (erow != NULL)
                erow[8 + j] = get_data(ctx, n,i);
#endif
        }
        fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
        if (erow != NULL)
            tda_export_row(ctx, "seqlg.table", erow, ecol);
#endif
        nrec++;
    }
#ifdef TDA_R_PACKAGE
    free(erow);
    tda_export_flush(ctx, "seqlg.table");
#endif
    prn_ssel(ctx, nrec);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);

    if (ctx->PMTDAFDef && nrec > 0)                      /* write TDA description */
        seq_dtda(ctx, 1,sn,ctx->PMFdName,nrec,ctx->PMNV,ctx->PMVIdx);

    err = 0;

SEQLFin:
    p_clean(ctx); 
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

int seqld(TDAContext *ctx)
{
    int err,sn,nrec;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);
    printf1(ctx, "Last occurrence of states.\n");
         
    if (parm(ctx, ctx->CmdBuf + 5,1,1))    /* get parameters */
        goto SEQLDFin;
        
    sn = seq_getsn(ctx, ctx->PMSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQLDFin;


    fprintf(ctx->PMFd,"#       starting  ending  sequence    last        total\n");
    fprintf(ctx->PMFd,"# Case      time    time    length  occurrence  occurrence\n");

/**********

    nrec = 0;
    for (i = 0; i < NOC; ++i) {

        if (eval_sve(ctx, i) == 0)   
            continue;

        len = gl = 0;
        ai = bi = -1;

        for (t = a; t <= b; ++t) {
            n = seq_sget(ctx, i,t,sn);
            if (n >= 0) {
                ai = t;
                break;
            }
        }
        if (ai >= 0) {
            for (t = b; t >= ai; --t) {
                n = seq_sget(ctx, i,t,sn);
                if (n >= 0) {
                    bi = t;
                    break;
                }
            }
            len = bi - ai + 1;
            l = 0;
            for (t = ai + 1; t <= bi; ++t) {
                n = seq_sget(ctx, i,t,sn);
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
            fprintf(PMFd,VPFmtS[n],get_data(ctx, n,i));
        }
        fprintf(PMFd,"\n");
        nrec++;
    }
***********/


    nrec = 0;
    prn_ssel(ctx, nrec);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);

    if (ctx->PMTDAFDef && nrec > 0)                      /* write TDA description */
        seq_dtda(ctx, 1,sn,ctx->PMFdName,nrec,ctx->PMNV,ctx->PMVIdx);

    err = 0;

SEQLDFin:
    p_clean(ctx); 
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  prn_ssel(n)   If sel option print number of selected cases.             */

void prn_ssel(TDAContext *ctx, int n)
{
    if (ctx->SVEFlg)  
        printf1(ctx, "Number of selected cases: %d\n",n);
}

/*--------------------------------------------------------------------------*/
/*  seq_csel()      return number of selected cases.                        */

int seq_csel(TDAContext *ctx)
{
    register int i;
    int n = 0;

    if (ctx->SVEFlg == 0)
        return(ctx->NOC);

    for (i = 0; i < ctx->NOC; ++i) {
        if (eval_sve(ctx, i) == 0)   
            continue;
        n++;
    }
    if (n < 1)  
        printf1(ctx, "Error: no cases selected.\n");
    return(n);
}

/*--------------------------------------------------------------------------*/
/*  seq_dtda(typ,sn,fname,noc,nv,vidx)     Write TDA description file.      */
/*                                                                          */

void seq_dtda(TDAContext *ctx, int typ,int sn,char *fname,int noc,int nv,short *vidx)
{
    int i = 0,j = 0,n = 0,ns = 0,s = 0,t = 0;

    fprintf(ctx->PMTDAFd,"nvar(\n");
    fprintf(ctx->PMTDAFd,"  dfile = %s,\n",fname);
    fprintf(ctx->PMTDAFd,"  noc = %d,\n",noc);

    switch (typ) {

        case  1:                    /* seqlg */
            fprintf(ctx->PMTDAFd,"  CASE <5>[6.0] = c1, # case number\n");
            fprintf(ctx->PMTDAFd,"  TS   <5>[6.0] = c2, # starting time\n");
            fprintf(ctx->PMTDAFd,"  TF   <5>[6.0] = c3, # ending time\n");
            fprintf(ctx->PMTDAFd,"  SLEN <5>[6.0] = c4, # sequence length\n");
            fprintf(ctx->PMTDAFd,"  GLEN <5>[6.0] = c5, # total gap length\n");
            fprintf(ctx->PMTDAFd,"  GMIN <5>[6.0] = c6, # minimum gap length\n");
            fprintf(ctx->PMTDAFd,"  GMAX <5>[6.0] = c7, # maximum gap length\n");
            fprintf(ctx->PMTDAFd,"  NGAP <5>[6.0] = c8, # number of gaps\n");
            n = 9;
            break;

        case  2:                    /* seqgc */

            ns = ctx->SeqSTN[sn];        /* number of states */

            fprintf(ctx->PMTDAFd,"  CASE <5>[8.0] = c1 , # case number\n");
            fprintf(ctx->PMTDAFd,"  SLEN <2>[4.0] = c2 , # sequence length\n");
            fprintf(ctx->PMTDAFd,"  NDS  <2>[4.0] = c3 , # number of different states\n");
            fprintf(ctx->PMTDAFd,"  NEV  <2>[4.0] = c4 , # number state changes\n");
            n = 5;
            for (j = 0; j < ns; ++j) {
                s = ctx->SeqSTNI[sn][j];
                fprintf(ctx->PMTDAFd,"  DUR%-2d<2>[4.0] = c%-2d, # duration in state %d\n",s,n++,s);
            }
            fprintf(ctx->PMTDAFd,"  DURM <2>[4.0] = c%-2d, # duration in missing state\n",n++);

            for (j = 0; j < ns; ++j) {
                s = ctx->SeqSTNI[sn][j];
                fprintf(ctx->PMTDAFd,"  NEP%-2d<2>[4.0] = c%-2d, # number of episodes in state %d\n",s,n++,s);
            }
            fprintf(ctx->PMTDAFd,"  NEPM <2>[4.0] = c%-2d, # number of episodes in missing state\n",n++);
            break;

        case  3:                    /* seqsd */

            ns = ctx->SeqSTN[sn];        /* number of states */

            fprintf(ctx->PMTDAFd,"  TIME  <2>[6.0] = c1, # time\n");
            n = 2;
            for (j = 0; j < ns; ++j) {
                s = ctx->SeqSTNI[sn][j];
                fprintf(ctx->PMTDAFd,"  NST%-2d <2>[6.0] = c%d, # cases in state %d\n",s,n++,s);
            }
            fprintf(ctx->PMTDAFd,"  VALID <2>[6.0] = c%d, # cases in valid states\n",n++);
            fprintf(ctx->PMTDAFd,"  NMISS <2>[6.0] = c%d, # cases in missing state\n",n++);
            fprintf(ctx->PMTDAFd,"  TOTAL <2>[6.0] = c%d, # total number of cases\n",n++);
            break;

        case  4:                    /* seqen */
            fprintf(ctx->PMTDAFd,"  TIME <2>[ 6.0] = c1, # time\n");
            fprintf(ctx->PMTDAFd,"  N    <2>[ 6.0] = c2, # number of cases\n");
            fprintf(ctx->PMTDAFd,"  ENT  <4>[%d.%d] = c3, # entropy\n",ctx->PMFmt1,ctx->PMFmt2);
            nv = 0;
            break;

        case  5:                    /* seqtp */

            ns = ctx->SeqSTN[sn];        /* number of states */

            fprintf(ctx->PMTDAFd,"  TIME  <2>[ 6.0] = c1 , # time\n");
            fprintf(ctx->PMTDAFd,"  N     <2>[ 6.0] = c2 , # number of cases\n");
            n = 3;
            for (i = 0; i < ns; ++i) {
                s = ctx->SeqSTNI[sn][i];
                for (j = 0; j <= ns; ++j) {
                    if (j < ns) {
                        t = ctx->SeqSTNI[sn][j];
                        fprintf(ctx->PMTDAFd,"  TP%d_%d <4>[%d.%d] = c%-2d, # transition %d - %d\n",
                                                 s,t,ctx->PMFmt1,ctx->PMFmt2,n++,s,t);
                    }
                    else {
                        fprintf(ctx->PMTDAFd,"  TP%d_M <4>[%d.%d] = c%-2d, # transition %d - missing\n",
                                                 s,ctx->PMFmt1,ctx->PMFmt2,n++,s);
                    }
                }
            }
            nv = 0;
            break;

        case  6:                    /* seqsi */

            ns = ctx->SeqSTN[sn];        /* number of states */

            fprintf(ctx->PMTDAFd,"  Case <2>[6.0] = c1 , # case number\n");
            n = 2;
            for (i = 0; i < ctx->PMNTP; ++i) {
                for (j = 0; j < ns; ++j) {
                    s = ctx->SeqSTNI[sn][j];
                    t = (int)ctx->PMTP[i];
                    fprintf(ctx->PMTDAFd,"  S%d_%d <2>[2.0] = c%-2d, # state %d, time %d\n",s,t,n++,s,t);
                }
            }
            break;

        default:
            break;
    }
    for (i = 0; i < nv; ++i) {
        j = vidx[i];
        fprintf(ctx->PMTDAFd,"  %s <%d>[%d.%d]",ctx->VName[j],ctx->VSLen[j],ctx->VPFmt1[j],ctx->VPFmt2[j]);
        if (ctx->VLabel[j] != NULL)                                                      
            fprintf(ctx->PMTDAFd,"(%s)",ctx->VLabel[j]);                                           
        fprintf(ctx->PMTDAFd," = c%d,\n",n++);
    }
    fprintf(ctx->PMTDAFd,");\n");
    
    printf1(ctx, "TDA description written to: %s\n",ctx->PMTDAFName);
}

/*--------------------------------------------------------------------------*/
/*  seqgc()         Characteristics of sequences                            */
/*                  seqgc(sn=,sel=,v=,dtda=) = fname                        */
/*                  Return 0 if OK, -1 if error.                            */

int seqgc(TDAContext *ctx)
{
#ifdef TDA_R_PACKAGE
    double *grow = NULL;
    int gcol = 0;
#endif
    register int i,j,k;
    int err,nrec,m,n,n1,nc,ns,sn,a,b,ai,bi,t,dm;

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);
         
    printf1(ctx, "Sequence characteristics.\n");
         
    if (parm(ctx, ctx->CmdBuf + 5,1,1))    /* get parameters */
        goto SEQFFin;
        
    sn = seq_getsn(ctx, ctx->PMSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQFFin;

    ns = ctx->SeqSTN[sn];        /* number of states */

    if (alloc_acn(ctx, ns + 1))
        goto SEQFFin;

    a = ctx->SeqTMin[sn];
    b = ctx->SeqTMax[sn];
    nrec = m = 0;
    for (i = 0; i < ctx->NOC; ++i) {

        if (eval_sve(ctx, i) == 0)   
            continue;

        ai = bi = -1;
        for (t = a; t <= b; ++t) {
            n = seq_sget(ctx, i,t,sn);
            if (n >= 0) {
                ai = t;
                break;
            }
        }
        if (ai < 0)
            continue;

        for (t = b; t >= ai; --t) {
            n = seq_sget(ctx, i,t,sn);
            if (n >= 0) {
                bi = t;
                break;
            }
        }
        fprintf(ctx->PMFd,"%6d ",i + 1);         /* case number */
        fprintf(ctx->PMFd,"%4d ",bi - ai + 1);   /* sequence length */
#ifdef TDA_R_PACKAGE
        if (grow == NULL) {
            gcol = 4 + 2 * (ns + 1) + ctx->PMNV;
            grow = (double *)malloc((size_t)gcol * sizeof(double));
        }
        if (grow != NULL) {
            grow[0] = (double)(i + 1);
            grow[1] = (double)(bi - ai + 1);
        }
#endif
        m++;

        for (j = 0; j < ns; ++j)
            ctx->AcN[j] = 0;

        nc = 0;
        n1 = 0;
        dm = 0;
        for (t = ai; t <= bi; ++t) {
            n = seq_sget(ctx, i,t,sn);
            if (t > ai && n != n1)      /* state changes */
                nc++;
            n1 = n;
            if (n >= 0) {
                n = ctx->SeqSTNII[sn][n];
                ctx->AcN[n] += 1;
            }
            else
                dm++;                   /* duration in missing state */
        }
        n = 0; 
        for (j = 0; j < ns; ++j) {      /* # of different states */
            if (ctx->AcN[j] > 0)
                n++;
        }
        fprintf(ctx->PMFd,"%4d ",n);
        fprintf(ctx->PMFd,"%4d ",nc);
#ifdef TDA_R_PACKAGE
        if (grow != NULL) {
            grow[2] = (double)n;
            grow[3] = (double)nc;
        }
#endif

        for (j = 0; j < ns; ++j) {      /* duration in state ... */
            fprintf(ctx->PMFd,"%4d ",ctx->AcN[j]);
#ifdef TDA_R_PACKAGE
            if (grow != NULL)
                grow[4 + j] = (double)ctx->AcN[j];
#endif
        }
        fprintf(ctx->PMFd,"%4d ",dm);
#ifdef TDA_R_PACKAGE
        if (grow != NULL)
            grow[4 + ns] = (double)dm;
#endif

        /* number of episodes */

        for (j = 0; j <= ns; ++j)
            ctx->AcN[j] = 0;

        n1 = seq_sget(ctx, i,ai,sn);
        k = ctx->SeqSTNII[sn][n1];
        ctx->AcN[k] += 1;

        for (t = ai + 1; t <= bi; ++t) {
            n = seq_sget(ctx, i,t,sn);
            if (n != n1) {
                if (n >= 0) {
                    k = ctx->SeqSTNII[sn][n];
                    ctx->AcN[k] += 1;
                }
                else
                    ctx->AcN[ns] += 1;
            }
            n1 = n;
        }
        for (j = 0; j <= ns; ++j) {      /* number of episodes */
            fprintf(ctx->PMFd,"%4d ",ctx->AcN[j]);
#ifdef TDA_R_PACKAGE
            if (grow != NULL)
                grow[5 + ns + j] = (double)ctx->AcN[j];
#endif
        }

        /* additional variables */

        for (j = 0; j < ctx->PMNV; ++j) {
            n = (int)ctx->PMVIdx[j];
            rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[n],get_data(ctx, n,i));
#ifdef TDA_R_PACKAGE
            if (grow != NULL)
                grow[6 + 2 * ns + j] = get_data(ctx, n,i);
#endif
        }
        fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
        if (grow != NULL)
            tda_export_row(ctx, "seqgc.table", grow, gcol);
#endif
        nrec++;
    }
#ifdef TDA_R_PACKAGE
    free(grow);
    tda_export_flush(ctx, "seqgc.table");
#endif
    prn_ssel(ctx, nrec);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);

    if (ctx->PMTDAFDef && nrec > 0)                      /* write TDA description */
        seq_dtda(ctx, 2,sn,ctx->PMFdName,nrec,ctx->PMNV,ctx->PMVIdx);
    err = 0;

SEQFFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqsd()         State distribution.                                     */
/*                  seqsd(sn=,sel=,dtda=) = fname                           */
/*                  Return 0 if OK, -1 if error.                            */

int seqsd(TDAContext *ctx)
{
    register int i,t;
    int err,n,a,b,sn,ns,nrec,nn;
#ifdef TDA_R_PACKAGE
    double *srow = NULL;
#endif

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);
         
    printf1(ctx, "State distributions.\n");
         
    if (parm(ctx, ctx->CmdBuf + 5,1,1))    /* get parameters */
        goto SEQSDFin;
        
    sn = seq_getsn(ctx, ctx->PMSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQSDFin;

    ns = ctx->SeqSTN[sn];        /* number of states */

    if (alloc_acn(ctx, ns + 1))
        goto SEQSDFin;

    a = ctx->SeqTMin[sn];
    b = ctx->SeqTMax[sn];

    nn = seq_csel(ctx);
    if (nn < 1)
        goto SEQSDFin;
         
    nrec = 0;
    for (t = a; t <= b; ++t) {

        for (i = 0; i <= ns; ++i)
            ctx->AcN[i] = 0;

        fprintf(ctx->PMFd,"%6d ",t);
#ifdef TDA_R_PACKAGE
        if (srow == NULL)
            srow = (double *)malloc((size_t)(ns + 4) * sizeof(double));
        if (srow != NULL)
            srow[0] = (double)t;
#endif

        for (i = 0; i < ctx->NOC; ++i) {
            if (eval_sve(ctx, i) == 0)   
                continue;
            n = seq_sget(ctx, i,t,sn);
            if (n >= 0)
                n = ctx->SeqSTNII[sn][n];
            else
                n = ns;
            ctx->AcN[n] += 1;
        }
        n = 0;
        for (i = 0; i < ns; ++i) {
            fprintf(ctx->PMFd,"%6d ",ctx->AcN[i]);
#ifdef TDA_R_PACKAGE
            if (srow != NULL)
                srow[1 + i] = (double)ctx->AcN[i];
#endif
            n += ctx->AcN[i];
        }
        fprintf(ctx->PMFd,"%6d %6d %6d\n",n,ctx->AcN[ns],nn);
#ifdef TDA_R_PACKAGE
        if (srow != NULL) {
            srow[1 + ns] = (double)n;
            srow[2 + ns] = (double)ctx->AcN[ns];
            srow[3 + ns] = (double)nn;
            tda_export_row(ctx, "seqsd.table", srow, ns + 4);
        }
#endif
        nrec++;
    }
#ifdef TDA_R_PACKAGE
    free(srow);
    tda_export_flush(ctx, "seqsd.table");
#endif
    prn_ssel(ctx, nn);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);

    if (ctx->PMTDAFDef && nrec > 0)                      /* write TDA description */
        seq_dtda(ctx, 3,sn,ctx->PMFdName,nrec,ctx->PMNV,ctx->PMVIdx);
    err = 0;

SEQSDFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqen()         Entropy measures.                                       */
/*                  seqen(sn=,sel=,dtda=,fmt=) = fname                      */
/*                  Return 0 if OK, -1 if error.                            */

int seqen(TDAContext *ctx)
{
    register int i,t;
    int err,m,n,a,b,sn,ns,nn,nrec;
    double tmp,tmp1;

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);
         
    printf1(ctx, "Entropy measures.\n");
         
    pmfmt(ctx, 10,4);                /* default print format */
    if (parm(ctx, ctx->CmdBuf + 5,1,1))   /* get parameters */
        goto SEQENFin;
        
    if (ctx->PMFmt1 == 0)
        pmfmt(ctx, 10,4);            /* default print format */

    sn = seq_getsn(ctx, ctx->PMSN,1);     /* get sequence number */
    if (sn < 0)  
        goto SEQENFin;

    ns = ctx->SeqSTN[sn];            /* number of states */

    if (alloc_acn(ctx, ns + 1))
        goto SEQENFin;

    a = ctx->SeqTMin[sn];
    b = ctx->SeqTMax[sn];

    nn = seq_csel(ctx);
    if (nn < 1)
        goto SEQENFin;
       
    fprintf(ctx->PMFd,"# Time      N    Entropy\n");

    nrec = 0;
    for (t = a; t <= b; ++t) {

        for (i = 0; i < ns; ++i)
            ctx->AcN[i] = 0;         

        m = 0;
        for (i = 0; i < ctx->NOC; ++i) {

            if (eval_sve(ctx, i) == 0)   
                continue;

            n = seq_sget(ctx, i,t,sn);
            if (n >= 0) {
                n = ctx->SeqSTNII[sn][n];
                ctx->AcN[n] += 1;
                m++;
            }
        }
        tmp = 0.0;
        if (m > 0) {
            for (i = 0; i < ns; ++i) {
                if (ctx->AcN[i] > 0) {
                    tmp1 = (double)ctx->AcN[i] / (double)m;
                    tmp -= tmp1 * rlog(ctx, tmp1);
                }
            }
        }
        fprintf(ctx->PMFd,"%6d %6d ",t,m);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
        fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
        {
            double erow[3];
            erow[0] = (double)t;
            erow[1] = (double)m;
            erow[2] = tmp;
            tda_export_row(ctx, "seqen.table", erow, 3);
        }
#endif
        nrec++;
    }
    prn_ssel(ctx, nn);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);

    if (ctx->PMTDAFDef && nrec > 0)                      /* write TDA description */
        seq_dtda(ctx, 4,sn,ctx->PMFdName,nrec,ctx->PMNV,ctx->PMVIdx);

    err = 0;

SEQENFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqtp()         Transition probabilities                                */
/*                  seqtp(sn=,sel=,fmt=,dtda=) = fname                      */
/*                  Return 0 if OK, -1 if error.                            */

int seqtp(TDAContext *ctx)
{
    register int i,j,t;
    int err,m,mm,mm0,n,n1,a,b,sn,ns,ns1,nss,nn,nrec;
    double tmp;

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Transition probabilities.\n");
         
    if (parm(ctx, ctx->CmdBuf + 5,1,1))    /* get parameters */
        goto SEQTPFin;
        
    sn = seq_getsn(ctx, ctx->PMSN,1);     /* get sequence number */
    if (sn < 0)  
        goto SEQTPFin;

    if (ctx->PMFmt1 == 0)
        pmfmt(ctx, 10,4);            /* default print format */

    ns = ctx->SeqSTN[sn];            /* number of states */
    ns1 = ns + 2;
    nss = ns * ns1;

    if (alloc_acn(ctx, nss + 1))
        goto SEQTPFin;

    a = ctx->SeqTMin[sn];
    b = ctx->SeqTMax[sn];

    nn = seq_csel(ctx);
    if (nn < 1)
        goto SEQTPFin;

    nrec = 0;
    for (t = a; t < b; ++t) {

        for (i = 0; i < nss; ++i)
            ctx->AcN[i] = 0;         

        mm0 = mm = m = 0;
        for (i = 0; i < ctx->NOC; ++i) {

            if (eval_sve(ctx, i) == 0)   
                continue;

            n = seq_sget(ctx, i,t,sn);
            if (n >= 0) {
                n = ctx->SeqSTNII[sn][n];
                m++;
                n1 = seq_sget(ctx, i,t + 1,sn);
                if (n1 >= 0) 
                    n1 = ctx->SeqSTNII[sn][n1];
                else  
                    n1 = ns;
                ctx->AcN[n * ns1 + n1] += 1;
                ctx->AcN[n * ns1 + ns + 1] += 1;
                if (n1 >= 0) {
                    mm0++;
                    if (n1 != n)
                        mm++;
                }
            }
        }
        fprintf(ctx->PMFd,"%6d %6d ",t,m);
#ifdef TDA_R_PACKAGE
        tda_export_cell(ctx, "seqtp.table", (double)t);
        tda_export_cell(ctx, "seqtp.table", (double)m);
#endif

        for (i = 0; i < ns; ++i) {
            n1 = ctx->AcN[i * ns1 + ns + 1];
            for (j = 0; j <= ns; ++j) {
                if (n1 <= 0)
                    tmp = 0.0;
                else 
                    tmp = (double)ctx->AcN[i * ns1 + j] / (double)n1;
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
#ifdef TDA_R_PACKAGE
                tda_export_cell(ctx, "seqtp.table", tmp);
#endif
            }
        }
        if (mm0 > 0)
            tmp = (double)mm / (double)mm0;
        else
            tmp = 0.0;
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tmp);
        fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
        tda_export_cell(ctx, "seqtp.table", tmp);
        tda_export_endrow(ctx, "seqtp.table");
#endif
        nrec++;
    }
    prn_ssel(ctx, nn);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);

    if (ctx->PMTDAFDef && nrec > 0)                      /* write TDA description */
        seq_dtda(ctx, 5,sn,ctx->PMFdName,nrec,ctx->PMNV,ctx->PMVIdx);

    err = 0;

SEQTPFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqrd()     Create random sequences.                                    */
/*              seqrd(ns=,noc=,len=)=fname.                                 */
/*              Return 0 if OK, -1 if error.                                */

int seqrd(TDAContext *ctx)
{
    register int i,j;
    int err,m,n;

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Creating random sequences.\n");
         
    if (parm(ctx, ctx->CmdBuf + 5,1,1))    /* get parameters */
        goto SEQRDFin;
        
    if (ctx->PMNS < 1)
        ctx->PMNS = 2;

    if (ctx->PMLEN < 1)
        ctx->PMLEN = 1;

    printf1(ctx, "Sequence length: %d, number of states: %d\n",ctx->PMLEN,ctx->PMNS);

    m = 0;
    for (i = 1; i <= ctx->PMNOC; ++i) {
        fprintf(ctx->PMFd,"%6d ",i);
#ifdef TDA_R_PACKAGE
        tda_export_cell(ctx, "seqrd.table", (double)i);
#endif
        for (j = 0; j < ctx->PMLEN; ++j) {
            n = (int)floor(random1(ctx) * (double)ctx->PMNS) + 1;
            fprintf(ctx->PMFd,"%2d ",n);
#ifdef TDA_R_PACKAGE
            tda_export_cell(ctx, "seqrd.table", (double)n);
#endif
        }
        fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
        tda_export_endrow(ctx, "seqrd.table");
#endif
        m++;
    }
    printf1(ctx, "%d sequences written to output file: %s\n",m,ctx->PMFdName);
    err = 0;

SEQRDFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqsi()     Create state indicator matrix.                              */
/*              seqsi(m=,tp=,sel=,dtda=,v=)=fname.                          */
/*              Return 0 if OK, -1 if error.                                */

int seqsi(TDAContext *ctx)
{
    register int i,j;
    int err,m,n,nt,ns,nst,sn,t;

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Printing state indicator matrix.\n");
         
    if (parm(ctx, ctx->CmdBuf + 5,1,1))    /* get parameters */
        goto SEQSIFin;
        
    if (ctx->PMM != 2)
        ctx->PMM = 1;

    if (ctx->PMNTP < 1) {
        p_err(ctx, -17,1);
        goto SEQSIFin;
    }
    sn = seq_getsn(ctx, ctx->PMSN,1);     /* get sequence number */
    if (sn < 0)  
        goto SEQSIFin;

    ns = ctx->SeqSTN[sn];    /* number of states */
    nst = ns * ctx->PMNTP;

    if (alloc_acn(ctx, nst + 1))
        goto SEQSIFin;

    m = 0;
    for (i = 0; i < ctx->NOC; ++i) {

        if (eval_sve(ctx, i) == 0)   
            continue;

        for (j = 0; j < nst; ++j)  
            ctx->AcN[j] = 0;

        nt = 0;
        for (j = 0; j < ctx->PMNTP; ++j) {
            t = (int)ctx->PMTP[j];
            n = seq_sget(ctx, i,t,sn);
            if (n >= 0) {
                n = ctx->SeqSTNII[sn][n];     /* internal state number */
                ctx->AcN[j * ns + n] = 1;
                nt++;
            }
            else if (ctx->PMM == 2)
                break;
        }
        if (ctx->PMM == 1 || nt == ctx->PMNTP) {
            fprintf(ctx->PMFd,"%6d ",i + 1);
#ifdef TDA_R_PACKAGE
            tda_export_cell(ctx, "seqsi.table", (double)(i + 1));
#endif
            for (j = 0; j < nst; ++j) {  
                fprintf(ctx->PMFd,"%1d ",ctx->AcN[j]);
#ifdef TDA_R_PACKAGE
                tda_export_cell(ctx, "seqsi.table", (double)ctx->AcN[j]);
#endif
            }

            for (j = 0; j < ctx->PMNV; ++j) {
                n = (int)ctx->PMVIdx[j];
                rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[n],get_data(ctx, n,i));
#ifdef TDA_R_PACKAGE
                tda_export_cell(ctx, "seqsi.table", get_data(ctx, n,i));
#endif
            }
            fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
            tda_export_endrow(ctx, "seqsi.table");
#endif
            m++;
        }
    }
    prn_ssel(ctx, m);
    printf1(ctx, "%d records written to output file: %s\n",m,ctx->PMFdName);

    if (ctx->PMTDAFDef && m > 0)                      /* write TDA description */
        seq_dtda(ctx, 6,sn,ctx->PMFdName,m,ctx->PMNV,ctx->PMVIdx);

    err = 0;

SEQSIFin:
    p_clean(ctx);
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

int seqpe(TDAContext *ctx)
{
    register int i,j;
    int err,m,a,b,i1,id,id1,s;
    double t,ts,tf;

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Creating sequences from episodes.\n");
         
    if (parm(ctx, ctx->CmdBuf + 5,1,1))    /* get parameters */
        goto SEQPEFin;
        
    if (ctx->PMMFlg == 0) 
        ctx->PMM = -1;
    if (ctx->PMM1Flg == 0) 
        ctx->PMM1 = -1;
    if (ctx->PMM2Flg == 0) 
        ctx->PMM2 = -1;

    if (ctx->PMID < 0 || ctx->PMORG < 0 || ctx->PMTS < 0 || ctx->PMTF < 0 || ctx->PMNTP < 1) {
        p_err(ctx, -4,1);
        goto SEQPEFin;
    }
    if (alloc_acn(ctx, ctx->PMNTP + 1))
        goto SEQPEFin;

    for (j = 0; j < ctx->PMNTP; ++j)
        ctx->AcN[j] = ctx->PMM;

    id1 = (int)get_data(ctx, ctx->PMID,0);
    i1 = 0;

    a = (int)(ctx->PMTP[0]);
    b = (int)(ctx->PMTP[ctx->PMNTP - 1]);
    m = 0;
                   
    for (i = 0; i < ctx->NOC; ++i) {

        id = (int)get_data(ctx, ctx->PMID,i);
        s  = (int)get_data(ctx, ctx->PMORG,i);
        ts = get_data(ctx, ctx->PMTS,i);
        tf = get_data(ctx, ctx->PMTF,i);

        if (id != id1) {
            m++,
            seqpe1(ctx, i1,m,ctx->AcN);
            id1 = id;
            i1 = i;
            for (j = 0; j < ctx->PMNTP; ++j)
                ctx->AcN[j] = ctx->PMM;
        }
        if (ts <= tf && ts <= b && tf >= a) {       
            for (j = 0; j < ctx->PMNTP; ++j) {
                t = ctx->PMTP[j];
                if (t >= ts && t < tf)
                    ctx->AcN[j] = s;
            }
        }           
    }
    m++;
    seqpe1(ctx, i1,m,ctx->AcN);
    printf1(ctx, "%d records written to output file: %s\n",m,ctx->PMFdName);
    err = 0;

SEQPEFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  seqpe1(i,id,st)     print sequence data for individual i with id and    */  
/*                      states st[].                                        */
/*                      called by seqpe().                                  */

void seqpe1(TDAContext *ctx, int i,int id,int *st)
{
    register int j,iv;

    /* insert new missing value codes */

    if (ctx->PMM1Flg) {
        for (j = 0; j < ctx->PMNTP; ++j) {
            if (st[j] == ctx->PMM)
                st[j] = ctx->PMM1;
            else
                break;
        }
    }
    if (ctx->PMM2Flg) {
        for (j = ctx->PMNTP - 1; j >= 0; --j) {
            if (st[j] == ctx->PMM)
                st[j] = ctx->PMM2;
            else
                break;
        }
    }

    fprintf(ctx->PMFd,"%6d ",id);
#ifdef TDA_R_PACKAGE
    tda_export_cell(ctx, "seqpe.table", (double)id);
#endif
    for (j = 0; j < ctx->PMNTP; ++j) {  
        fprintf(ctx->PMFd,"%2d ",st[j]);
#ifdef TDA_R_PACKAGE
        tda_export_cell(ctx, "seqpe.table", (double)st[j]);
#endif
    }

    for (j = 0; j < ctx->PMNV; ++j) {
        iv = (int)ctx->PMVIdx[j];
        rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[iv],get_data(ctx, iv,i));
#ifdef TDA_R_PACKAGE
        tda_export_cell(ctx, "seqpe.table", get_data(ctx, iv,i));
#endif
    }
    fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
    tda_export_endrow(ctx, "seqpe.table");
#endif
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

int seqpd(TDAContext *ctx)
{
    register int i,j,k,l;
    int m,n,nm,err,org,ts,sn,iv,min,max,ta,tb,ns = 0,ng,gmin,gmax,n1,len;
    int sa = 0,sb = 0,lcen,rcen;
   
    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Printing sequence data.\n");

    if (ctx->SeqDN == 0) {       /* check for sequence data */
        p_err(ctx, -13,1);
        goto SEQPDFin;
    }
    if (parm(ctx, ctx->CmdBuf + 5,1,1))    /* get parameters */
        goto SEQPDFin;

    min = ctx->SeqTCMin; 
    max = ctx->SeqTCMax; 
        
    m = 0;         
    switch (ctx->PMM) {

        case 1:    
            for (i = 0; i < ctx->NOC; ++i) {

                if (eval_sve(ctx, i) == 0)   
                    continue;

                fprintf(ctx->PMFd,"%6d ",i + 1);

                for (j = 0; j <= ctx->SeqDNH; ++j) {
                    if (ctx->SeqDT[j]) {

                        if (ctx->PMNS == 1) {
                            nm = ctx->SeqSTNH[j];
                            if (alloc_acn(ctx, nm + 1))
                                goto SEQPDFin;
                            for (k = 0; k <= nm; ++k)
                                ctx->AcN[k] = -10;
                            n1 = -1;
                            for (k = min; k <= max; ++k) {
                                n = seq_sget(ctx, i,k,j);
                                if (n >= 0 && n < nm && n != n1)
                                    ctx->AcN[n] += 10;
                                n1 = n;
                                if (n >= 0)
                                    n += ctx->AcN[n];
                                fprintf(ctx->PMFd,"%3d ",n);
                            }
                        }
                        else {
                            for (k = min; k <= max; ++k) {
                                n = seq_sget(ctx, i,k,j);
                                fprintf(ctx->PMFd,"%2d ",n);
                            }
                        }
                    }
                }
                for (k = 0; k < ctx->PMNV; ++k) {
                    iv = (int)ctx->PMVIdx[k];
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[iv],get_data(ctx, iv,i));
                }
                fprintf(ctx->PMFd,"\n");
                m++;
            }
            break;

        case 2:
            for (i = 0; i < ctx->NOC; ++i) {
                if (eval_sve(ctx, i) == 0)   
                    continue;

                for (j = 0; j <= ctx->SeqDNH; ++j) {
                    if (ctx->SeqDT[j]) {
                        fprintf(ctx->PMFd,"%6d %2d ",i + 1,j + 1);
                        for (k = min; k <= max; ++k) {
                            n = seq_sget(ctx, i,k,j);
                            fprintf(ctx->PMFd,"%2d ",n);
                        }
                        for (k = 0; k < ctx->PMNV; ++k) {
                            iv = (int)ctx->PMVIdx[k];
                            rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[iv],get_data(ctx, iv,i));
                        }
                        fprintf(ctx->PMFd,"\n");
                        m++;
                    }
                }
            }
            break;

        case 3:
            for (i = 0; i < ctx->NOC; ++i) {

                if (eval_sve(ctx, i) == 0)   
                    continue;

                for (k = min; k <= max; ++k) {
                    fprintf(ctx->PMFd,"%6d %4d ",i + 1,k);
                    for (j = 0; j <= ctx->SeqDNH; ++j) {
                        if (ctx->SeqDT[j]) {
                            n = seq_sget(ctx, i,k,j);
                            fprintf(ctx->PMFd,"%2d ",n);
                        }
                    }
                    for (j = 0; j < ctx->PMNV; ++j) {
                        iv = (int)ctx->PMVIdx[j];
                        rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[iv],get_data(ctx, iv,i));
                    }
                    fprintf(ctx->PMFd,"\n");
                    m++;
                }
            }
            break;

        case 4:

            for (j = 0; j <= ctx->SeqDNH; ++j) {
                if (ctx->SeqDT[j])  
                    break;
            }
            for (i = 0; i < ctx->NOC; ++i) {

                if (eval_sve(ctx, i) == 0)   
                    continue;

                sn = 1;
                ts = k = min;
                org = seq_sget(ctx, i,k,j);
                while (++k <= max + 1) {
                    if (k <= max)
                        n = seq_sget(ctx, i,k,j);
                    else
                        n = org;

                    if (n != org || k == max + 1) {
                        fprintf(ctx->PMFd,"%6d %2d %2d %2d %4d %4d ",i + 1,sn++,org,n,ts,k);
                        for (l = 0; l < ctx->PMNV; ++l) {
                            iv = (int)ctx->PMVIdx[l];
                            rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[iv],get_data(ctx, iv,i));
                        }
                        fprintf(ctx->PMFd,"\n");
                        m++;
                        org = n;
                        ts = k;
                    }
                }
            }
            break;
         
        case 5:
            for (i = 0; i < ctx->NOC; ++i) {

                if (eval_sve(ctx, i) == 0)   
                    continue;

                fprintf(ctx->PMFd,"%6d ",i + 1);

                for (j = 0; j <= ctx->SeqDNH; ++j) {
                    if (ctx->SeqDT[j]) {

                        lcen = rcen = 0;
                        ta = tb = -1;
                        for (k = min; k <= max; ++k) {
                            n = seq_sget(ctx, i,k,j);
                            if (n >= 0) {
                                ta = k;
                                sa = n;
                                if (k == min)
                                    lcen = 1;
                                break;
                            }
                        }
                        for (k = max; k >= min; --k) {
                            n = seq_sget(ctx, i,k,j);
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
                            n = seq_sget(ctx, i,k,j);
                            if (n >= 0) {
                                ns++;   
                                if (len > 0) {
                                    gmin = imin(ctx, len,gmin);
                                    gmax = imax(ctx, len,gmax);
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
    
                        fprintf(ctx->PMFd,"%5d %5d %2d %2d %1d %1d %5d %5d %5d %5d  ",
                            ta,tb,sa,sb,lcen,rcen,ns,ng,gmin,gmax);
                                     
                    }
                }
                for (k = 0; k < ctx->PMNV; ++k) {
                    iv = (int)ctx->PMVIdx[k];
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[iv],get_data(ctx, iv,i));
                }
                fprintf(ctx->PMFd,"\n");
                m++;
            }
            break;
       
        case 6:
        case 7:

            for (j = 0; j <= ctx->SeqDNH; ++j) {
                if (ctx->SeqDT[j])  
                    break;
            }

            ns = 0;
            for (i = 0; i < ctx->NOC; ++i) {

                if (eval_sve(ctx, i) == 0)   
                    continue;

                ta = min;
                tb = max;
                for (k = min; k <= max; ++k) {
                    n = seq_sget(ctx, i,k,j);
                    if (n >= 0) {
                        ta = k;
                        break;
                    }
                }
                for (k = max; k >= min; --k) {
                    n = seq_sget(ctx, i,k,j);
                    if (n >= 0) {
                        tb = k;
                        break;
                    }
                }
                n1 = -1;
                ng = 0;
                for (k = ta; k <= tb; ++k) {
                    n = seq_sget(ctx, i,k,j);
                    if (ctx->PMM == 7 && n < 0)
                        continue;

                    if (n != n1)
                        ng++;
                    n1 = n;
                }
                ns = imax(ctx, ns,ng);
            }
            for (i = 0; i < ctx->NOC; ++i) {

                if (eval_sve(ctx, i) == 0)   
                    continue;

                fprintf(ctx->PMFd,"%6d ",i + 1);

                ta = min;
                tb = max;
                for (k = min; k <= max; ++k) {
                    n = seq_sget(ctx, i,k,j);
                    if (n >= 0) {
                        ta = k;
                        break;
                    }
                }
                for (k = max; k >= min; --k) {
                    n = seq_sget(ctx, i,k,j);
                    if (n >= 0) {
                        tb = k;
                        break;
                    }
                }
                n1 = -1;
                ng = 0;
                for (k = ta; k <= tb; ++k) {
                    n = seq_sget(ctx, i,k,j);
                    if (ctx->PMM == 7 && n < 0)
                        continue;

                    if (n != n1) {
                        ng++;
                        fprintf(ctx->PMFd,"%2d ",n);
                    }
                    n1 = n;
                }
                for (k = ng + 1; k <= ns; ++k)
                    fprintf(ctx->PMFd,"%2d ",-9);

                fprintf(ctx->PMFd,"%5d ",ng);
                  
                for (k = 0; k < ctx->PMNV; ++k) {
                    iv = (int)ctx->PMVIdx[k];
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[iv],get_data(ctx, iv,i));
                }
                fprintf(ctx->PMFd,"\n");
                m++;
            }
            break;

        case 8:                        
            if (ctx->PMS < 1)
                ctx->PMS = 1;
            printf1(ctx, "Aggregation with s = %d\n",ctx->PMS);

            for (j = 0; j <= ctx->SeqDNH; ++j) {         /* use first sequence */
                if (ctx->SeqDT[j])  
                    break;
            }
            ns = ctx->SeqSTN[j];        /* number of states */
            if (alloc_acn(ctx, ns + 1))
                goto SEQPDFin;

            min = ctx->SeqTMin[j];
            max = ctx->SeqTMax[j];

            for (i = 0; i < ctx->NOC; ++i) {

                if (eval_sve(ctx, i) == 0)   
                    continue;

                fprintf(ctx->PMFd,"%6d ",i + 1);

                k = min;
                while (k <= max) {

                    for (l = 0; l < ns; ++l)
                        ctx->AcN[l] = 0;
                    n1 = 0;
                    for (l = k; l < k + ctx->PMS; ++l) {
                        if (l > max)
                            break;
                        n = seq_sget(ctx, i,l,j);

                        if (n >= 0) {
                            iv = ctx->SeqSTNII[j][n];
                            ctx->AcN[iv] += 1;
                        }
                        else
                            n1++;                   /* duration in missing state */
                    }
                    n = -1;  
                    for (iv = 0; iv < ns; ++iv) {
                        if (ctx->AcN[iv] > n1) {
                            n1 = ctx->AcN[iv];
                            n  = ctx->SeqSTNI[j][iv];
                        }   
                    }
                    fprintf(ctx->PMFd,"%2d ",n);
                    k += ctx->PMS;
                }
                for (k = 0; k < ctx->PMNV; ++k) {
                    iv = (int)ctx->PMVIdx[k];
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[iv],get_data(ctx, iv,i));
                }
                fprintf(ctx->PMFd,"\n");
                    m++;
            }
            break;

        case 9:                                               
            for (i = 0; i < ctx->NOC; ++i) {

                if (eval_sve(ctx, i) == 0)   
                    continue;

                fprintf(ctx->PMFd,"%6d ",i + 1);

                for (j = 0; j <= ctx->SeqDNH; ++j) {
                    if (ctx->SeqDT[j]) {
                        /*********************************
                        for (k = 0; k < SeqSTN[j]; ++k) {
                            n = SeqSTNI[j][k];
                            printf1(ctx, " %d",n);
                        }
                        printf1(ctx, "\n");
                        ********************************/

                        n1 = ctx->SeqSTN[j];  
                        if (alloc_acn(ctx, n1 + 1))
                            goto SEQPDFin;

                        for (k = min; k <= max; ++k) {
                            n = seq_sget(ctx, i,k,j);
                            if (n >= 0)  
                                ctx->AcN[ctx->SeqSTNII[j][n]] += 1;
                            else
                                ctx->AcN[n1] += 1;
                        }   
                        for (k = 0; k < n1; ++k)
                            fprintf(ctx->PMFd,"%6d ",ctx->AcN[k]);
                        fprintf(ctx->PMFd,"%6d ",ctx->AcN[n1]);
                    }
                }
                for (k = 0; k < ctx->PMNV; ++k) {
                    iv = (int)ctx->PMVIdx[k];
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[iv],get_data(ctx, iv,i));
                }
                fprintf(ctx->PMFd,"\n");
                m++;
            }
            break;

        case 10:                                     /* ### */
            for (i = 0; i < ctx->NOC; ++i) {

                if (eval_sve(ctx, i) == 0)   
                    continue;

                fprintf(ctx->PMFd,"%6d ",i + 1);

                for (j = 0; j <= ctx->SeqDNH; ++j) {
                    if (ctx->SeqDT[j]) {
                        /*********************************
                        for (k = 0; k < SeqSTN[j]; ++k) {
                            n = SeqSTNI[j][k];
                            printf1(ctx, " %d",n);
                        }
                        printf1(ctx, "\n");
                        ********************************/

                        n1 = ctx->SeqSTN[j];  
                        if (alloc_acn(ctx, n1 + 1))
                            goto SEQPDFin;
                        if (alloc_acm(ctx, n1 + 1))
                            goto SEQPDFin;

                        for (k = 0; k <= n1; ++k)
                            ctx->AcN[k] = ctx->AcM[k] = -5;

                        for (k = min; k <= max; ++k) {
                            n = seq_sget(ctx, i,k,j);
                            if (n >= 0) {
                                n = ctx->SeqSTNII[j][n];
                                if (ctx->AcN[n] < 0)
                                    ctx->AcN[n] = k;
                                if (ctx->AcM[n] < k)
                                    ctx->AcM[n] = k;
                            }
                            else {
                                if (ctx->AcN[n1] < 0)
                                    ctx->AcN[n1] = k;
                                if (ctx->AcM[n1] < k)
                                    ctx->AcM[n1] = k;
                            }
                        }   
                        for (k = 0; k <= n1; ++k)
                            fprintf(ctx->PMFd,"%6d %6d  ",ctx->AcN[k],ctx->AcM[k]);
                    }
                }
                for (k = 0; k < ctx->PMNV; ++k) {
                    iv = (int)ctx->PMVIdx[k];
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[iv],get_data(ctx, iv,i));
                }
                fprintf(ctx->PMFd,"\n");
                m++;
            }
            break;
        
        default:
            printf1(ctx, "Error: option m=%d not available.\n",ctx->PMM);
            goto SEQPDFin;
    }
    prn_ssel(ctx, m);
    printf1(ctx, "%d records written to output file: %s\n",m,ctx->PMFdName);

    if (ctx->PMTDAFDef && m > 0) {                    /* write TDA description */

        fprintf(ctx->PMTDAFd,"nvar(\n");
        fprintf(ctx->PMTDAFd,"  dfile = %s,\n",ctx->PMFdName);
        fprintf(ctx->PMTDAFd,"  noc = %d,\n",m);

        fprintf(ctx->PMTDAFd,"  CASE <5>[6.0] = c1 , # case number\n");
        n = 2;

        switch (ctx->PMM) {

            case 1:
                for (j = 0; j <= ctx->SeqDNH; ++j) {
                    if (ctx->SeqDT[j]) {
                        for (k = min; k <= max; ++k)  
                            fprintf(ctx->PMTDAFd,"  Y%d_%d <2>[2.0] = c%-2d, # state in seq %d, time %d\n",j,k,n++,j + 1,k);
                    }   
                }
                break;

            case  2:                               
                fprintf(ctx->PMTDAFd,"  SN   <2>[2.0] = c%-2d, # sequence number\n",n++);
                for (k = min; k <= max; ++k)  
                    fprintf(ctx->PMTDAFd,"  Y%-3d <2>[2.0] = c%-2d, # state at time %d\n",k,n++,k);
                break;

            case  3:                               
                fprintf(ctx->PMTDAFd,"  TIME <5>[4.0] = c%-2d, # time\n",n++);
                for (j = 0; j <= ctx->SeqDNH; ++j) {
                    if (ctx->SeqDT[j])  
                        fprintf(ctx->PMTDAFd,"  Y%-3d <2>[2.0] = c%-2d, # state in seq %d\n",j,n++,j + 1);
                }
                break;

            case  4:                               
                fprintf(ctx->PMTDAFd,"  SN   <2>[2.0] = c%-2d, # spell number\n",n++);
                fprintf(ctx->PMTDAFd,"  ORG  <5>[2.0] = c%-2d, # origin state\n",n++);
                fprintf(ctx->PMTDAFd,"  DES  <5>[2.0] = c%-2d, # destination state\n",n++);
                fprintf(ctx->PMTDAFd,"  TS   <5>[4.0] = c%-2d, # starting time\n",n++);
                fprintf(ctx->PMTDAFd,"  TF   <5>[4.0] = c%-2d, # ending time\n",n++);
                break;

            case 5:
                for (j = 0; j <= ctx->SeqDNH; ++j) {
                    if (ctx->SeqDT[j]) {
                        fprintf(ctx->PMTDAFd,"  TA%-2d   <5>[5.0] = c%-2d, # seq %d, first time point\n",j,n++,j + 1);
                        fprintf(ctx->PMTDAFd,"  TB%-2d   <5>[5.0] = c%-2d, # seq %d, last time point\n",j,n++,j + 1);
                        fprintf(ctx->PMTDAFd,"  SA%-2d   <2>[2.0] = c%-2d, # seq %d, state at TA\n",j,n++,j + 1);
                        fprintf(ctx->PMTDAFd,"  SB%-2d   <2>[2.0] = c%-2d, # seq %d, state at TB\n",j,n++,j + 1);
                        fprintf(ctx->PMTDAFd,"  LCen%-2d <2>[1.0] = c%-2d, # seq %d, 1 if left censored\n",j,n++,j + 1);
                        fprintf(ctx->PMTDAFd,"  RCen%-2d <2>[1.0] = c%-2d, # seq %d, 1 if right censored\n",j,n++,j + 1);
                        fprintf(ctx->PMTDAFd,"  NS%-2d   <5>[5.0] = c%-2d, # seq %d, number of valid states\n",j,n++,j + 1);
                        fprintf(ctx->PMTDAFd,"  NP%-2d   <5>[5.0] = c%-2d, # seq %d, number of parts\n",j,n++,j + 1);
                        fprintf(ctx->PMTDAFd,"  GLMin%-2d<5>[5.0] = c%-2d, # seq %d, minimal gap length\n",j,n++,j + 1);
                        fprintf(ctx->PMTDAFd,"  GLMax%-2d<5>[5.0] = c%-2d, # seq %d, maximal gap length\n",j,n++,j + 1);
                    }   
                }
                break;

            case  6:                               
            case  7:                               
                for (j = 1; j <= ns; ++j)  
                    fprintf(ctx->PMTDAFd,"  Y%-3d<2>[2.0] = c%-2d, # state\n",j,n++);
                fprintf(ctx->PMTDAFd,"  NS  <5>[2.0] = c%-2d, # number of states\n",n++);
                break;

            case 8:
                for (j = 0; j <= ctx->SeqDNH; ++j) {
                    if (ctx->SeqDT[j])  
                        break;
                }
                k = min;
                l = 0;
                while (k <= max) {
                    fprintf(ctx->PMTDAFd,"  Y%d <2>[2.0] = c%-2d, # state in time %d - %d\n",l++,n++,k,imin(ctx, max,k + ctx->PMS - 1));
                    k += ctx->PMS;
                }
                break;

            case 9:
                for (j = 0; j <= ctx->SeqDNH; ++j) {
                    if (ctx->SeqDT[j]) {
                        n1 = ctx->SeqSTN[j];  
                        for (k = 0; k < n1; ++k) {
                            ns = ctx->SeqSTNI[j][k];
                            fprintf(ctx->PMTDAFd,"  D%d_%d <5>[6.0] = c%-2d, # duration in seqence %d, state %d\n",j,ns,n++,j + 1,ns);
                        }
                        fprintf(ctx->PMTDAFd,"  D%d_m <5>[6.0] = c%-2d, # duration in seqence %d, state missing\n",j,n++,j + 1);
                    }   
                }
                break;

            case 10:
                for (j = 0; j <= ctx->SeqDNH; ++j) {
                    if (ctx->SeqDT[j]) {
                        n1 = ctx->SeqSTN[j];  
                        for (k = 0; k < n1; ++k) {
                            ns = ctx->SeqSTNI[j][k];
                            fprintf(ctx->PMTDAFd,"  F%d_%d <5>[6.0] = c%-2d, # first time point in seqence %d, state %d\n",j,ns,n++,j + 1,ns);
                            fprintf(ctx->PMTDAFd,"  L%d_%d <5>[6.0] = c%-2d, #  last time point in seqence %d, state %d\n",j,ns,n++,j + 1,ns);
                        }
                        fprintf(ctx->PMTDAFd,"  F%d_m <5>[6.0] = c%-2d, # first time point in seqence %d, state missing\n",j,n++,j + 1);
                        fprintf(ctx->PMTDAFd,"  L%d_m <5>[6.0] = c%-2d, #  last time point in seqence %d, state missing\n",j,n++,j + 1);
                    }   
                }
                break;

            default:
                break;
        }   

        for (i = 0; i < ctx->PMNV; ++i) {
            j = ctx->PMVIdx[i];
            fprintf(ctx->PMTDAFd,"  %s <%d>[%d.%d]",ctx->VName[j],ctx->VSLen[j],ctx->VPFmt1[j],ctx->VPFmt2[j]);
            if (ctx->VLabel[j] != NULL)                                                      
                fprintf(ctx->PMTDAFd,"(%s)",ctx->VLabel[j]);                                           
            fprintf(ctx->PMTDAFd," = c%-2d,\n",n++);
        }
        fprintf(ctx->PMTDAFd,");\n");
    
        printf1(ctx, "TDA description written to: %s\n",ctx->PMTDAFName);
    }
    err = 0;

SEQPDFin:
    p_clean(ctx);
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

int seqsn(TDAContext *ctx)
{
    register int i,j;
    int err,nrec,n,ni,nj,ns,ns1 = 0,sn,a,b,t,s,na,nb,nflag;
    double d;

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);
         
    printf1(ctx, "States as new objects.\n");
         
    if (parm(ctx, ctx->CmdBuf + 5,1,1))    /* get parameters */
        goto SEQSNFin;
        
    if (ctx->PMFmt1 == 0)
        pmfmt(ctx, 10,4);            /* default print format */

    sn = seq_getsn(ctx, ctx->PMSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQSNFin;

    ns = ctx->SeqSTN[sn];        /* number of states */

    if (alloc_acn(ctx, ns + 1))
        goto SEQSNFin;

    a = ctx->SeqTMin[sn];
    b = ctx->SeqTMax[sn];

    nrec = 0;

    switch (ctx->PMM) {

        case 1:    
        case 2:    
            for (t = a; t <= b; ++t) {

                na = nb = 0;
                for (j = 0; j < ns; ++j)
                    ctx->AcN[j] = 0;

                for (i = 0; i < ctx->NOC; ++i) {

                    n = seq_sget(ctx, i,t,sn);
                    if (n < 0)
                        continue;

                    ni = ctx->SeqSTNII[sn][n];
                    if (ctx->PMM == 1)
                        ctx->AcN[ni] += 1;
                    else if (ctx->PMM == 2)
                        ctx->AcN[ni] = i;      
    
                    if (t > a) {
                        if (seq_sget(ctx, i,t - 1,sn) < 0)
                            na++;
                    }
                    if (t < b) {
                        if (seq_sget(ctx, i,t + 1,sn) < 0)
                            nb++;
                    }
                }
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,t);         /* time */
                n = 0;
                for (j = 0; j < ns; ++j) {
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcN[j]);  
                    if (ctx->PMM == 1)
                        n += ctx->AcN[j];
                    else if (ctx->AcN[j] > 0)
                        n++;
                }        
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,n);  
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,na);  
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nb);  
                fprintf(ctx->PMFd,"\n");
                nrec++;
            }
            break;

        case 3:                                 /* permutations */
            for (t = a + 1; t <= b; ++t) {

                for (j = 0; j < ns; ++j)
                    ctx->AcN[j] = -1;

                for (i = 0; i < ctx->NOC; ++i) {

                    na = seq_sget(ctx, i,t - 1,sn);
                    if (na < 0)
                        continue;

                    nb = seq_sget(ctx, i,t,sn);
                    if (nb < 0)
                        continue;

                    ni = ctx->SeqSTNII[sn][na];
                    nj = ctx->SeqSTNII[sn][nb];
                    ctx->AcN[ni] = nj;
                }
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,t);         /* time */
                for (j = 0; j < ns; ++j)   
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcN[j] + 1);  
                fprintf(ctx->PMFd,"\n");
                nrec++;
            }
            break;

        case 4:    

            if (alloc_aci(ctx, ns + 1))
                goto SEQSNFin;
            if (alloc_acj(ctx, ns + 1))
                goto SEQSNFin;
            if (alloc_ack(ctx, ns + 1))
                goto SEQSNFin;
            if (alloc_acr(ctx, ns + 1))
                goto SEQSNFin;
            if (alloc_acs(ctx, ns + 1))
                goto SEQSNFin;

            for (i = 0; i < ctx->NOC; ++i) {

                for (j = 0; j < ns; ++j)
                    ctx->AcI[j] = 0;

                nflag = 1;
                for (t = a; t <= b; ++t) {
                    n = seq_sget(ctx, i,t,sn);
                    if (n < 0) {
                        nflag = 1;
                        continue;
                    }
                    ni = ctx->SeqSTNII[sn][n];
                    ctx->AcI[ni] = 1;

                    if (t > a && n != seq_sget(ctx, i,t - 1,sn))
                        ctx->AcR[ni] += 1;
                    if (t < b && n != seq_sget(ctx, i,t + 1,sn))
                        ctx->AcS[ni] += 1;

                    ctx->AcJ[ni] += 1;
                    if (nflag || (t > a && n != seq_sget(ctx, i,t - 1,sn)))
                        ctx->AcK[ni] += 1;
                    nflag = 0;
                }
                for (j = 0; j < ns; ++j)
                    ctx->AcN[j] += ctx->AcI[j];
            }

            for (j = 0; j < ns; ++j) {
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j + 1);  
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->SeqSTNI[sn][j]);  
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcN[j]);  
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcR[j]);  
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcS[j]);  
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcJ[j]);  
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcK[j]);  

                if (ctx->AcK[j] > 0)
                    d = (double)ctx->AcJ[j] / (double)ctx->AcK[j];
                else
                    d = 0.0;
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,d);  
                fprintf(ctx->PMFd,"\n");
                nrec++;
            }
            break;

        case 5:    

            ns1 = ns + 1;
            if (alloc_acn(ctx, ns1 * ns1 + 1))
                goto SEQSNFin;

            for (i = 0; i < ctx->NOC; ++i) {
                for (t = a + 1; t <= b; ++t) {
                    nb = seq_sget(ctx, i,t,sn);
                    na = seq_sget(ctx, i,t - 1,sn);
                    if (na < 0 && nb < 0)
                        continue;

                    if (na < 0) {
                        ni = ctx->SeqSTNII[sn][nb];
                        ctx->AcN[ns * ns1 + ni] += 1;
                    }
                    else if (nb < 0) {
                        ni = ctx->SeqSTNII[sn][na];
                        ctx->AcN[ni * ns1 + ns] += 1;
                    }
                    else {
                        ni = ctx->SeqSTNII[sn][na];
                        nj = ctx->SeqSTNII[sn][nb];
                        ctx->AcN[ni * ns1 + nj] += 1;
                    }
                }
            }
            for (i = 0; i <= ns; ++i) {

                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);    

                if (i < ns)
                    n = ctx->SeqSTNI[sn][i];
                else
                    n = -1;

                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,n);    

                for (j = 0; j <= ns; ++j)  
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcN[i * ns1 + j]);  

                fprintf(ctx->PMFd,"\n");
                nrec++;
            }
            break;

        default:
            printf1(ctx, "Error: option m=%d not available.\n",ctx->PMM);
            goto SEQSNFin;
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);

    if (ctx->PMTDAFDef && nrec > 0) {                    /* write TDA description */

        fprintf(ctx->PMTDAFd,"nvar(\n");
        fprintf(ctx->PMTDAFd,"  dfile = %s,\n",ctx->PMFdName);
        fprintf(ctx->PMTDAFd,"  noc = %d,\n",nrec);

        n = 1;

        switch (ctx->PMM) {

            case 1:
            case 2:
                fprintf(ctx->PMTDAFd,"  T   <5>[%d.0] = c%-2d, # time\n",ctx->PMNFmt,n++);
                for (j = 0; j < ns; ++j) {
                    s = ctx->SeqSTNI[sn][j];
                    fprintf(ctx->PMTDAFd,"  N%-3d<5>[%d.0] = c%-2d, # number of sequences (or ID) in state %d\n",s,ctx->PMNFmt,n++,s);
                }
                fprintf(ctx->PMTDAFd,"  N   <5>[%d.0] = c%-2d, # total number of sequences\n",ctx->PMNFmt,n++);
                fprintf(ctx->PMTDAFd,"  NA  <5>[%d.0] = c%-2d, # number of new sequences\n",ctx->PMNFmt,n++);
                fprintf(ctx->PMTDAFd,"  NB  <5>[%d.0] = c%-2d, # number of dropped sequences\n",ctx->PMNFmt,n++);
                break;

            case 4:
                fprintf(ctx->PMTDAFd,"  RN  <5>[%d.0] = c%d, # record number\n",ctx->PMNFmt,n++);
                fprintf(ctx->PMTDAFd,"  S   <5>[%d.0] = c%d, # state\n",ctx->PMNFmt,n++);
                fprintf(ctx->PMTDAFd,"  NS  <5>[%d.0] = c%d, # number of sequences in state\n",ctx->PMNFmt,n++);
                fprintf(ctx->PMTDAFd,"  NA  <5>[%d.0] = c%d, # number of sequences entering state\n",ctx->PMNFmt,n++);
                fprintf(ctx->PMTDAFd,"  NB  <5>[%d.0] = c%d, # number of sequences leaving state\n",ctx->PMNFmt,n++);
                fprintf(ctx->PMTDAFd,"  TT  <5>[%d.0] = c%d, # total time in state\n",ctx->PMNFmt,n++);
                fprintf(ctx->PMTDAFd,"  NE  <5>[%d.0] = c%d, # number of episodes in state\n",ctx->PMNFmt,n++);
                fprintf(ctx->PMTDAFd,"  DUR <4>[%d.%d] = c%d, # mean duration (TT/NE)\n",ctx->PMFmt1,ctx->PMFmt2,n++);
                break;

            case 5:
                fprintf(ctx->PMTDAFd,"  RN  <5>[%d.0] = c%d, # record number\n",ctx->PMNFmt,n++);
                fprintf(ctx->PMTDAFd,"  S   <5>[%d.0] = c%d, # state\n",ctx->PMNFmt,n++);
                for (j = 1; j <= ns1; ++j)
                    fprintf(ctx->PMTDAFd,"  N%-3d<5>[%d.0] = c%d,\n",j,ctx->PMNFmt,n++);
                break;

            default:
                break;
        }   
        fprintf(ctx->PMTDAFd,");\n");
        printf1(ctx, "TDA description written to: %s\n",ctx->PMTDAFName);
    }
    err = 0;

SEQSNFin:
    p_clean(ctx);
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

int seqgm(TDAContext *ctx)
{
    register int i,j;
    int err,nrec,n,sn,a,b,t,na,nb,na1,na2,d;
    int nn[11];
    double dd;

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);
         
    printf1(ctx, "Indicators of group membership.\n");
         
    if (parm(ctx, ctx->CmdBuf + 5,1,1))    /* get parameters */
        goto SEQMLFin;
        
    if (ctx->PMFmt1 == 0)
        pmfmt(ctx, 10,4);            /* default print format */

    sn = seq_getsn(ctx, ctx->PMSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQMLFin;

    a = ctx->SeqTMin[sn];
    b = ctx->SeqTMax[sn];

    fprintf(ctx->PMFd,"# Time Number_in_state_%d",ctx->PMS);
    for (j = 1; j <= 10; ++j) {
        if (ctx->PMSK[j] > 0) 
            fprintf(ctx->PMFd," Len=%d ",ctx->PMSK[j]);  
    }
    fprintf(ctx->PMFd,"\n");

    nrec = 0;
    for (t = a; t <= b; ++t) {

        d = na = na1 = na2 = 0;
        for (j = 1; j <= 10; ++j)
            nn[j] = 0;

        for (i = 0; i < ctx->NOC; ++i) {

            n = seq_sget(ctx, i,t,sn);
            if (n != ctx->PMS)
                continue;
            na++;
            if (t == a)
                na1++;
            else if (seq_sget(ctx, i,t - 1,sn) != ctx->PMS)
                na1++;

            if (t == b)
                na2++;
            else if (seq_sget(ctx, i,t + 1,sn) != ctx->PMS)
                na2++;

            d++;
            for (j = t - 1; j >= a; --j) {
                if (seq_sget(ctx, i,j,sn) != ctx->PMS) 
                    break;
                d++;
            }
            nb = 1;   
            for (j = t - 1; j >= a; --j) {
                if (seq_sget(ctx, i,j,sn) != ctx->PMS)
                    break;
                nb++;
            }
            for (j = 1; j <= 10; ++j) {
                if (nb >= ctx->PMSK[j])
                    nn[j] += 1;
            }
        }
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,t);         /* time */
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,na);  
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,na1);  
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,na2);  

        dd = 0.0;
        if (na > 0)
            dd = (double)d / (double)na;

        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,dd);  
        for (j = 1; j <= 10; ++j) {
            if (ctx->PMSK[j] > 0) 
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nn[j]);  
        }
        fprintf(ctx->PMFd,"\n");
        nrec++;
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

SEQMLFin:
    p_clean(ctx);
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

int seqgsp(TDAContext *ctx)
{
    register int i,j;
    int err,nrec,n,sn,a,b,t,f,l;   

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);
         
    printf1(ctx, "Generation of data for sequence plots.\n");
         
    if (parm(ctx, ctx->CmdBuf + 6,1,1))    /* get parameters */
        goto SEQGSPFin;
        
    sn = seq_getsn(ctx, ctx->PMSN,1);   /* get sequence number */
    if (sn < 0)  
        goto SEQGSPFin;

    a = ctx->SeqTMin[sn];
    b = ctx->SeqTMax[sn];

    if (alloc_acn(ctx, ctx->NOC + 1))
        goto SEQGSPFin; 
    if (alloc_acm(ctx, ctx->NOC + 1))
        goto SEQGSPFin; 
    if (alloc_ack(ctx, ctx->NOC + 1))
        goto SEQGSPFin; 

    for (i = 0; i < ctx->NOC; ++i) {
        f = -1;
        l = 0;
        for (t = a; t <= b; ++t) {
            n = seq_sget(ctx, i,t,sn);
            if (n != ctx->PMS)
                continue;
            if (f == -1)
                f = t;
            l++;
        }
        ctx->AcN[i] = f;
        ctx->AcM[i] = l;
    }
    if (sortdpi2(ctx, ctx->NOC,ctx->AcN,ctx->AcM,ctx->AcK)) 
        goto SEQGSPFin;

    /******************* 
    tda_out("AcN: ");
    for (i = 0; i < NOC; ++i)
        tda_out(" %d ",AcN[i]);

    tda_out("\n AcM: ");
    for (i = 0; i < NOC; ++i)
        tda_out(" %d ",AcM[i]);
    tda_out("\n AcK: ");
    for (i = 0; i < NOC; ++i)
        tda_out(" %d ",AcK[i]);
    newline(ctx);
    *******************/

    nrec = 0;
    for (t = a; t <= b; ++t) {
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,t);         /* time */
                 
        for (i = 0; i < ctx->NOC; ++i) {
            j = ctx->AcK[i];
            n = seq_sget(ctx, j,t,sn);
            if (n != ctx->PMS)
                n = 0;
            else
                n = ctx->NOC - i; 

            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,n);
        }
        fprintf(ctx->PMFd,"\n");
        nrec++;
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

SEQGSPFin:
    p_clean(ctx);
    return(err);
}

