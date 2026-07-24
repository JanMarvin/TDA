/****************************************************************************/
/*  t_dmu                                                                   */
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
#include "t_freq.h"
#include "t_gdat.h"
#include "t_alloc.h"
#include "t_gf.h"
#include "t_sort.h"
#include "tda_context.h"
#include "tda_compat.h"

/* ------------------------------------------------------------------------ */
/*  functions in tda_dm.                                                    */

int dm_ccnt(TDAContext *ctx);
int dm_lcnt(TDAContext *ctx);
int dm_dump(TDAContext *ctx);
void prnbuf(TDAContext *ctx, unsigned char *buf, int n, int off);
int dm_dsplit(TDAContext *ctx);
int esort(TDAContext *ctx);
int escomp(const void *, const void *, void *);
long get_long(TDAContext *ctx, char *p);
void put_long(TDAContext *ctx, char *p,long n);
int eskip(TDAContext *ctx);
int eselect(TDAContext *ctx);
int eselcomp(const void *, const void *, void *);
int esel_c(TDAContext *ctx, int klen);
int emerge(TDAContext *ctx);
int ejoin(TDAContext *ctx);
int ejoin_nv(TDAContext *ctx, FILE *fd,char *buf,int *nrec);
int ejoin_rd(TDAContext *ctx, FILE *fd,char *fname,char *buf,int *eof,int *nrec, int *bvar,double *xvar,int nx,int *mxlev);
int ejoin_proc(TDAContext *ctx, int id,int r1,int r2,int l1,int l2,int nx1,int nx2);
void ejoin_prn(TDAContext *ctx, int nfmt,int val);
void ejoin_prn1(TDAContext *ctx, int nfmt,double val);
void ejoin_prn2(TDAContext *ctx, int nfmt,double val);

/*--------------------------------------------------------------------------*/
/*  Global variables                                                        */



/* ------------------------------------------------------------------------ */
/*  dm_ccnt()   Command: ccnt=fname. Print frequency distribution of        */
/*              characters in fname to standard output.                     */
/*              Return 0 if OK, -1 if error.                                */

int dm_ccnt(TDAContext *ctx)
{
    register int i,j;
    register char *p;
    int err,nn,cnt,rcnt;

    err = -1;
    nn = 20000;         /* size of read buffer */

    if (check_cmd(ctx, 1))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 4,10,1)) {   /* get parameters */
        p_clean(ctx);
        return(-1);   
    }
    printf1(ctx, "Frequency distribution of characters in file: %s\n",ctx->PMFdName);

    if (alloc_acn(ctx, 257))
        goto CCNTFin;
   
    if (alloc_acc(ctx, nn + 1))
        goto CCNTFin;
   
    rcnt = 0;
    while (((cnt = (int)fread(ctx->AcC,sizeof(char),(size_t)(nn),ctx->PMFd))) > 0) {

        rcnt += cnt;
        p = ctx->AcC;
        for (i = 0; i < cnt; ++i) {
            j = (int) *p++;
            if (j < 0)
                j += 256;
            ctx->AcN[j] += 1;
        }
    }
    printf1(ctx, "\nCharacter   Hexadecimal      Count\n");
    prnchar(ctx, '-',34,1);

    for (i = 0; i < 256; ++i) {
        if (ctx->AcN[i]) {
            if (isprint(i))  
                printf1(ctx, "    %c",(char)i);
            else  
                printf1(ctx, "     ");
            printf1(ctx, "            %02X     %10d\n",i,ctx->AcN[i]);
        }
    }
    prnchar(ctx, '-',34,1);
    printf1(ctx, "Sum                     %10d\n",rcnt);
    err = 0;

CCNTFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  dm_lcnt()   Command: lcnt(noc=)=fname. Print frequency distribution of  */
/*              length of records to standard output.                       */
/*              Return 0 if OK, -1 if error.                                */
/*                                                                          */
/*  Note: record length is printed without EOL!                             */  

int dm_lcnt(TDAContext *ctx)
{
    register int i,j,l;
    register char *p;
    int err,nn,n,rec,wflg,wflg1;
    double d1,d2;

    err = -1;
    wflg = wflg1 = 0;
    nn = 20000;         /* size of read buffer */

    if (check_cmd(ctx, 1))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 4,2,1)) {   /* get parameters */
        p_clean(ctx);
        return(-1);   
    }
    printf1(ctx, "Frequency distribution of record lengths (without EOL) in file: %s\n",ctx->PMFdName);
    printf1(ctx, "Max number of records: %d\n",ctx->PMNOC);

    if (alloc_acc(ctx, nn + 2))
        goto LCNTFin;

    if (alloc_acn(ctx, ctx->PMNOC + 1))
        goto LCNTFin;

    rec = 0;
    while (fgets(ctx->AcC,nn + 1,ctx->PMFd)) {
        if (rec >= ctx->PMNOC) {
            wflg = 1;
            break;
        }
        l = (int)(strlen(ctx->AcC));
        p = ctx->AcC + l;
        while (p > ctx->AcC && (*--p == LF || *p == CR))
            l--;

        ctx->AcN[++rec] = l;
        if (l >= nn)
            wflg1 = 1;
    }
    printf1(ctx, "Read %d records\n",rec);
    if (wflg)  
        printf1(ctx, "Warning: file has more than %d records\n",rec);
    if (wflg1)  
        printf1(ctx, "Warning: record length may have exceeded the maximum of %d bytes.\n",nn);

    if (rec <= 0)  
        goto LCNTFin;

    if (alloc_aci(ctx, rec + 1))
        goto LCNTFin;

    if (alloc_acj(ctx, rec + 1))
        goto LCNTFin;

    if (alloc_acm(ctx, rec + 1))
        goto LCNTFin;

    n = cfreq(ctx, 1,rec,ctx->AcN,rec,ctx->AcM,ctx->AcI,ctx->AcJ,0,&d1,&d2);
    if (n <= 0)  
        printf1(ctx, "Error: can't sort length of records.\n");
    else {
        printf1(ctx, "Index  Length  Frequency\n");
        for (i = 1; i <= n; ++i) {
            j = ctx->AcJ[i];
            printf1(ctx, "%5d %7d %10d\n",i,ctx->AcM[j],ctx->AcI[j]);
        }
    }
    err = 0;

LCNTFin:
    p_clean(ctx);    
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  dm_dump     Command: dump(nc=,s=)=fname.                                */
/*              hex dump of fname, nc is number of characters, s is offset  */
/*              Return 0 if OK, -1 if error.                                */

int dm_dump(TDAContext *ctx)
{
    int err,nn,n,off,rcnt,pcnt;
        
    err = -1;
    nn = 32000;         /* size of read buffer, should be multiple of 32 */

    if (check_cmd(ctx, 1))
        return(-1);
   
    if (parm(ctx, ctx->CmdBuf + 4,10,1)) {   /* get parameters */
        p_clean(ctx);
        return(-1);   
    }
    if (ctx->PMS < 0)
        ctx->PMS = 0;

    printf1(ctx, "Hex dump of file: %s [offset=%d]\n",ctx->PMFdName,ctx->PMS);
    if (ctx->PMS > 0) {
        if (fseek(ctx->PMFd,(long)ctx->PMS,0)) {
            printf1(ctx, "Can't seek to offset %d.\n",ctx->PMS);
            goto DUMPFin;
        }
        off = ctx->PMS;
    }
    else
        off = 0;

    if (alloc_acc(ctx, nn + 2))
        goto DUMPFin;

    rcnt = nn;         
    pcnt = 0;

    while (1) {

        if (ctx->PMNC > 0 && rcnt > ctx->PMNC)
            rcnt = ctx->PMNC;

        if (((n = (int)fread(ctx->AcC,sizeof(char),(size_t)(rcnt),ctx->PMFd))) < 0) {
            printf1(ctx, "Error in reading file: %s\n",ctx->PMFdName);
            goto DUMPFin;
        }
        if (!n)       
            break;

        prnbuf(ctx, (unsigned char *)ctx->AcC,n,off);
        off += n; 

        pcnt += n;
        if (ctx->PMNC > 0 && pcnt >= ctx->PMNC)
            break;
    }
    printf1(ctx, "\n");
    err = 0;

DUMPFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prnbuf(buf,n,off)                                                       */
/*      print n char of buf in hex and ASCII to stdout. off is used to      */
/*      count the header column.                                            */

void prnbuf(TDAContext *ctx, unsigned char *buf, int n, int off)
{
    register int ii,i;
    register unsigned char c;
    int llen = 16;                  /* characters per line */

    i = 0;
    while (1) {
        printf1(ctx, "\n%08x : ",off + i);
        for (ii = 0; ii < llen; ++ii) {
            if ((i + ii) >= n)
                break;
            printf1(ctx, "%02x ",buf[i + ii]);
        }
        for ( ; ii < llen; ++ii)  
            printf1(ctx, "   ");
        printf1(ctx, "  ");

        for (ii = 0; ii < llen; ++ii) {
            if ((i + ii) >= n)
                break;
            c = buf[i + ii];
            if (c >= 0x20 && c < 0x7f)  
                printf1(ctx, "%c",c);
            else  
                printf1(ctx, ".");
        }
        if ((i += llen) >= n)
            break;
    }
}

/* ------------------------------------------------------------------------ */
/*  dm_dsplit.  Command: dsplit(len=)=fname. Split file into parts.         */
/*              Return 0 if OK, -1 if error.                                */

int dm_dsplit(TDAContext *ctx)
{
    FILE *fd;
    int err,cnt,rcnt;
    char *p,fname[1000];

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 6,10,1)) {   /* get parameters */
        p_clean(ctx);
        return(-1);   
    }
    if (ctx->PMLEN < 1)
        ctx->PMLEN = 1000;

    printf1(ctx, "Split file %s into parts of length %d.\n",ctx->PMFdName,ctx->PMLEN);

    if (alloc_acc(ctx, ctx->PMLEN + 2))
        goto DSPLFin;

    strcpy(fname,ctx->PMFdName);
    strcat(fname,".a");
    p = fname + strlen(fname) - 1;

    rcnt = 0;
    while (((cnt = (int)fread(ctx->AcC,sizeof(char),(size_t)(ctx->PMLEN),ctx->PMFd))) > 0) {
        rcnt += cnt;
        printf1(ctx, "%s -- ",fname);

        if (!(fd = fopen(fname,OPEN_WB))) {
            printf1(ctx, "can't create this output file\n");
            break;         
        }
        if (fwrite(ctx->AcC,sizeof(char),(size_t)(cnt),fd) != (size_t)(cnt)) {
            printf1(ctx, "write error\n");
            break;
        }
        printf1(ctx, "written %d bytes\n",cnt);
        fclose(fd);
        if (*p == 'z') {
            p++;
            *p++ = 'a';
            *p-- = '\0';
        }
        else
            *p += 1;
    }
    printf1(ctx, "\nRead from input file %s: %d bytes.\n",ctx->PMFdName,rcnt);
    err = 0;
 
DSPLFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  esort()     Command: external sort.                                     */
/*                                                                          */
/*              esort(                                                      */
/*                  df = ...,               output file, required           */
/*                  sk = c1,...,c10,        up to 5 sort keys               */
/*                  noc=...,                max number of records, def.1000 */
/*                                          per temporary file              */
/*                  len=...,                max record length, def. 1000    */
/*                                                                          */
/*              ) = fname;  input file                                      */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

#define ESMAXNF 100         /* max number of temp. files */

int esort(TDAContext *ctx)
{
    FILE *ofd[ESMAXNF];
    register int i,j,k;
    int err,low[5],up[5],klen,kmax,keya,rcnt,tcnt,wcnt,nf,nn;
    long fptr,wptr;
    register char *p,*q,*s;
    char *kp,fname[500];

    keya = 0;
    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "External sort. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,2,1))     /* get parameters */
        goto ESFin;

    printf1(ctx, "Inputfile: %s\nOutputfile: ",ctx->PMFdName);
    if (ctx->PMF1Def)
        printf1(ctx, "%s\n",ctx->PMF1dName);
    else {
        printf1(ctx, "need output file.\n");
        goto ESFin;
    }

    printf1(ctx, "Sort keys:");
    klen = kmax = ctx->ESNK = 0;
    for (i = 1; i < 10; i += 2) {
        if (ctx->PMSK[i] > 0 && ctx->PMSK[i + 1] > 0) {
            if (ctx->PMSK[i + 1] < ctx->PMSK[i]) {
                printf1(ctx, "error in key definition.\n");
                goto ESFin;
            }
            low[ctx->ESNK] = ctx->PMSK[i] - 1;
            up [ctx->ESNK] = ctx->PMSK[i + 1];
            printf1(ctx, " [%d,%d]",low[ctx->ESNK] + 1,up[ctx->ESNK]);
            ctx->ESLEN[ctx->ESNK] = up[ctx->ESNK] - low[ctx->ESNK];
            klen += ctx->ESLEN[ctx->ESNK] + 1;
            if (kmax < up[ctx->ESNK])
                kmax = up[ctx->ESNK];

            ctx->ESNK++;
        }
    }
    if (ctx->ESNK == 0) {
        printf1(ctx, " need a valid definition with sk parameter.\n");
        goto ESFin;
    }
    printf1(ctx, "\nMax number of records: %d\n",ctx->PMNOC);

    if (ctx->PMLEN < 1)
        ctx->PMLEN = 1000;

    printf1(ctx, "Max record length: %d\n",ctx->PMLEN);
    if (ctx->PMLEN < kmax) {
        printf1(ctx, "Inconsistent with key definition.\n");
        goto ESFin;
    }
    if (ctx->PMNOC < ESMAXNF)
        ctx->PMNOC = ESMAXNF;

    klen += (int)sizeof(long);
    if (!(ctx->ESKEY  = (char *)calloc((size_t)(ctx->PMNOC) * (size_t)(klen),sizeof(char)))) {
        p_err(ctx, -2,1);             
        goto ESFin;    
    }         
    keya = 1;
    memrq(ctx, ctx->PMNOC * klen,sizeof(char));

    if (alloc_acc(ctx, ctx->PMLEN + 2))
        goto ESFin;   

    nf = tcnt = 0;
    fptr = ftell(ctx->PMFd);

    while (1) {

        if (fseek(ctx->PMFd,fptr,0)) {
            printf1(ctx, "Error: can't seek to offset %ld.\n",fptr);
            goto ESFin;
        }

        rcnt = 0;
        while (fgets(ctx->AcC,ctx->PMLEN + 1,ctx->PMFd)) {
           
            /* save key */

            k = (int)(strlen(ctx->AcC));
            p = ctx->AcC + k;
            for (i = k; i < kmax; ++i)
                *p++ = ' ';

        	   p = ctx->ESKEY + rcnt * klen;                 
            put_long(ctx, p,fptr);
            p += sizeof(long);     
	    
            for (k = 0; k < ctx->ESNK; ++k) {
                q = ctx->AcC + low[k];
                for (i = 0; i < ctx->ESLEN[k]; ++i)
                    *p++ = *q++;
                *p++ = '\0';
            }
       	    fptr = ftell(ctx->PMFd);
            tcnt++;
            prn_message(ctx, tcnt,0,0);
	  
            if (++rcnt >= ctx->PMNOC)
                break;
        }
        prn_message(ctx, tcnt,1,0);
         
        if (rcnt == 0)
            break;

        if (nf >= ESMAXNF) {
            printf1(ctx, "\nError: exceeded max number of temporary files.\n");
            goto ESFin;
        }
        strcpy(fname,ctx->PMF1dName);
        snprintf(fname + strlen(fname),sizeof(fname) - strlen(fname),".%d",nf);

        if (!(ofd[nf] = fopen(fname,OPEN_WR))) {
            printf1(ctx, "\nError: can't create a new temporary output file (%s)\n",fname);
            goto ESFin;    
        }
        if (ctx->SILENTFlg < 2)
            printfe(ctx, "Sorting: %s\n",fname);

        tda_qsort_r((char *)ctx->ESKEY,(size_t)(rcnt),(size_t)(klen), escomp, ctx);

        /* write to output file */

        if (ctx->SILENTFlg < 2)
            printfe(ctx, "Writing: %s\n",fname);

        for (i = 0; i < rcnt; ++i) {

            wptr = get_long(ctx, ctx->ESKEY + i * klen);

            if (fseek(ctx->PMFd,wptr,0)) {
                printf1(ctx, "Error: can't seek to offset %ld.\n",wptr);
                goto ESFin;
            }
            if (!fgets(ctx->AcC,ctx->PMLEN + 1,ctx->PMFd)) {
                printf1(ctx, "Error: can't read at offset %ld.\n",wptr);
                goto ESFin;
            }
            fprintf(ofd[nf],"%s",ctx->AcC);
        }
        fclose(ofd[nf]);
        nf++;

    }
    printf1(ctx, "\nRead %d records.\n",tcnt);
    printf1(ctx, "Number of temporary files: %d\n",nf);
          
    /* MERGING */

    if (ctx->SILENTFlg < 2)
        printfe(ctx, "\nMerging ...\n");

    if (alloc_acc(ctx, nf * (ctx->PMLEN + 2)))
        goto ESFin;   

    s = ctx->AcC;
    for (i = 0; i < nf; ++i) {
        strcpy(fname,ctx->PMF1dName);
        snprintf(fname + strlen(fname),sizeof(fname) - strlen(fname),".%d",i);
        if (!(ofd[i] = fopen(fname,OPEN_RD))) {
            printf1(ctx, "\nError: can't open temporary file (%s).\n",fname);
            goto ESFin;    
        }
        if (!fgets(s,ctx->PMLEN + 1,ofd[i])) {
            printf1(ctx, "Error: can't read temp file %s.\n",fname);
            goto ESFin;
        }
        p = ctx->ESKEY + i * klen;                 
        put_long(ctx, p,(long)i);
        p += sizeof(long);

        for (k = 0; k < ctx->ESNK; ++k) {
            q = s + low[k];
            for (j = 0; j < ctx->ESLEN[k]; ++j)
                *p++ = *q++;
            *p++ = '\0';
        }
        s += ctx->PMLEN + 2;
    }
    kp = ctx->ESKEY;
    nn = nf;
    wcnt = 0;
    while (nn > 0) {

        if (nn > 1)
            tda_qsort_r((char *)kp,(size_t)(nn),(size_t)(klen), escomp, ctx);
    
        wptr = get_long(ctx, kp);

        j = (int)wptr;
        s = ctx->AcC + (ctx->PMLEN + 2) * j;
        fprintf(ctx->PMF1d,"%s",s);
        wcnt++;
        prn_message(ctx, wcnt,0,1);

        if (!fgets(s,ctx->PMLEN + 1,ofd[j])) {
            *s = '\0';
            nn--;
            kp += klen;      
        }   
        else {
            p = kp + sizeof(long);

            for (k = 0; k < ctx->ESNK; ++k) {
                q = s + low[k];
                for (j = 0; j < ctx->ESLEN[k]; ++j)
                    *p++ = *q++;
                *p++ = '\0';
            }
        }
    }
    for (i = 0; i < nf; ++i)
        fclose(ofd[i]);

    prn_message(ctx, wcnt,1,1);
    printf1(ctx, "%d records written to: %s\n",wcnt,ctx->PMF1dName);
       
    err = 0;
 
ESFin:
    if (keya) {
        free((char *)ctx->ESKEY);
        memrq(ctx, -ctx->PMNOC * klen,sizeof(char));
    }
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  escomp()    compare function                                            */

int escomp(const void *arg1, const void *arg2, void *_ctx)
{
    TDAContext *ctx = (TDAContext *)_ctx;     
    register int k,n = 0;
    register char *p,*q;

    p = (char *)arg1 + sizeof(long);
    q = (char *)arg2 + sizeof(long);

    for (k = 0; k < ctx->ESNK; ++k) {
        n = strcmp(p,q);
        if (n)
            return(n);
        p += ctx->ESLEN[k] + 1;
        q += ctx->ESLEN[k] + 1;
    }
    return(n);
}

/* ------------------------------------------------------------------------ */
/*  get_long()                                                              */

long get_long(TDAContext *ctx, char *p)
{     
    (void)ctx;        /* unused: the signature is shared */
    register int i;
    char lc[sizeof(long)];
    register char *q;

    q = lc;
    for (i = 0; (size_t)i < sizeof(long); ++i)
        *q++ = *p++;
    {
        long lv;
        memcpy(&lv,lc,sizeof(lv));   /* a char array carries no alignment: copy the bytes */
        return(lv);
    }
}

/* ------------------------------------------------------------------------ */
/*  put_long()                                                              */

void put_long(TDAContext *ctx, char *p,long n)
{     
    (void)ctx;        /* unused: the signature is shared */
    register int i;
    char lc[sizeof(long)];
    register char *q;

    memcpy(lc,&n,sizeof(n));   /* a char array carries no alignment: copy the bytes */

    q = lc;
    for (i = 0; (size_t)i < sizeof(long); ++i)
        *p++ = *q++;
}

/* ------------------------------------------------------------------------ */
/*  eskip()     Command: delete selected columns                            */
/*                                                                          */
/*              eskip(                                                      */
/*                  df = ...,               output file, required           */
/*                  sk = ...,               up to 5 column ranges           */
/*                  len=...,                max record length, def. 1000    */
/*                                                                          */
/*              ) = fname;  input file                                      */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int eskip(TDAContext *ctx)
{
    register int i,k;
    int err,low[5],up[5],kmax,rcnt;
    register char *p;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Skip columns. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,2,1))     /* get parameters */
        goto ESKFin;

    printf1(ctx, "Inputfile: %s\nOutputfile: ",ctx->PMFdName);
    if (ctx->PMF1Def)
        printf1(ctx, "%s\n",ctx->PMF1dName);
    else {
        printf1(ctx, "need output file.\n");
        goto ESKFin;
    }
    printf1(ctx, "Column ranges:");
    kmax = ctx->ESNK = 0;
    for (i = 1; i < 10; i += 2) {
        if (ctx->PMSK[i] > 0 && ctx->PMSK[i + 1] > 0) {
            if (ctx->PMSK[i + 1] < ctx->PMSK[i]) {
                printf1(ctx, "error in sk definition.\n");
                goto ESKFin;
            }
            low[ctx->ESNK] = ctx->PMSK[i];
            up [ctx->ESNK] = ctx->PMSK[i + 1];
            printf1(ctx, " [%d,%d]",low[ctx->ESNK],up[ctx->ESNK]);
            if (kmax < up[ctx->ESNK])
                kmax = up[ctx->ESNK];

            ctx->ESNK++;
        }
    }
    if (ctx->ESNK == 0) {
        printf1(ctx, " need a valid definition with sk parameter.\n");
        goto ESKFin;
    }
    newline(ctx);
    if (ctx->PMLEN < 1)
        ctx->PMLEN = 1000;

    printf1(ctx, "Max record length: %d\n",ctx->PMLEN);
    if (ctx->PMLEN < kmax) {
        printf1(ctx, "Inconsistent with sk definition.\n");
        goto ESKFin;
    }
    if (alloc_acc(ctx, ctx->PMLEN + 2))
        goto ESKFin;   
    if (alloc_acn(ctx, ctx->PMLEN + 2))
        goto ESKFin;   

    for (k = 0; k < ctx->ESNK; ++k) {
        for (i = low[k]; i <= up[k]; ++i)
            ctx->AcN[i] = 1;
    }
    rcnt = 0;
    while (fgets(ctx->AcC,ctx->PMLEN + 1,ctx->PMFd)) {
           
        prn_message(ctx, ++rcnt,0,0);

        p = ctx->AcC;
        for (i = 1; i <= ctx->PMLEN + 1; ++i) {
            if (ctx->AcN[i] == 0 || *p == '\n')
                fprintf(ctx->PMF1d,"%c",*p);
            if (!*p || *p == '\n') 
                break;
            p++;
        }
    }
    prn_message(ctx, rcnt,1,0);
    printf1(ctx, "%d records written to: %s\n",rcnt,ctx->PMF1dName);
       
    err = 0;
 
ESKFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  eselect()   Command: select records.                                    */
/*                                                                          */
/*              eselect(                                                    */
/*                  df = ...,               output file, required           */
/*                  if = ...,               input file with list            */
/*                  sk = c1,...,c10,        2 keys                          */
/*                  noc=...,                max number of keys, def. 1000   */
/*                  len=...,                max record length, def. 1000    */
/*                  opt=...,                1 use all records               */
/*                                          2 max one record per key        */  
/*                                                                          */
/*              ) = fname;  input file                                      */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int eselect(TDAContext *ctx)
{
    register int i;
    int err,low[5],up[5],kmax,klen,rcnt,wcnt,keya;
    register char *p,*q,*s;

    keya = 0;
    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Selection of records. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 7,2,1))     /* get parameters */
        goto ESELFin;

    printf1(ctx, "Inputfile: %s\nOutputfile: ",ctx->PMFdName);
    if (ctx->PMF1Def)
        printf1(ctx, "%s\n",ctx->PMF1dName);
    else {
        printf1(ctx, "need output file.\n");
        goto ESELFin;
    }
    printf1(ctx, "Selection keys: ");
    if (ctx->PMF2Def)
        printf1(ctx, "%s\n",ctx->PMF2dName);
    else {
        printf1(ctx, "need a file.\n");
        goto ESELFin;
    }
    printf1(ctx, "Key definition:");
    kmax = ctx->ESNK = 0;
    for (i = 1; i < 10; i += 2) {
        if (ctx->PMSK[i] > 0 && ctx->PMSK[i + 1] > 0) {
            if (ctx->PMSK[i + 1] < ctx->PMSK[i]) {
                printf1(ctx, "error in sk definition.\n");
                goto ESELFin;
            }
            low[ctx->ESNK] = ctx->PMSK[i] - 1;
            up [ctx->ESNK] = ctx->PMSK[i + 1];
            printf1(ctx, " [%d,%d]",low[ctx->ESNK] + 1,up[ctx->ESNK]);
            if (kmax < up[ctx->ESNK])
                kmax = up[ctx->ESNK];

            ctx->ESNK++;
        }
    }
    if (ctx->ESNK != 2) {
        printf1(ctx, "\nError: need exactly two keys.\n");
        goto ESELFin;
    }
    newline(ctx);
    klen = up[0] - low[0];
    if (klen != up[1] - low[1]) {
        printf1(ctx, "Error: keys should have identical length.\n");
        goto ESELFin;
    }
    printf1(ctx, "Max number of keys: %d\n",ctx->PMNOC);
    if (ctx->PMLEN < 1)
        ctx->PMLEN = 1000;

    printf1(ctx, "Max record length: %d\n",ctx->PMLEN);
    if (ctx->PMLEN < kmax) {
        printf1(ctx, "Inconsistent with sk definition.\n");
        goto ESELFin;
    }
    if (alloc_acc(ctx, ctx->PMLEN + 2))
        goto ESELFin;   
    if (alloc_acd(ctx, klen + 1))
        goto ESELFin;   

    if (!(ctx->ESKEY  = (char *)calloc((size_t)(ctx->PMNOC) * (size_t)((klen + 1)),sizeof(char)))) {
        p_err(ctx, -2,1);             
        goto ESELFin;    
    }         
    keya = 1;
    memrq(ctx, ctx->PMNOC * (klen + 1),sizeof(char));

    /* read keys from if file */

    s = ctx->AcC + low[0];
    p = ctx->ESKEY;
    ctx->ESNK = 0;
    while (fgets(ctx->AcC,ctx->PMLEN + 1,ctx->PMF2d)) {
        if (ctx->ESNK >= ctx->PMNOC) {
            printf1(ctx, "Error: exceeded max number of keys.\n");
            goto ESELFin;
        }
        q = s;             
        for (i = 0; i < klen; ++i)
            *p++ = *q++;
        *p++ = '\0';
        prn_message(ctx, ++ctx->ESNK,0,0);
    }
    prn_message(ctx, ctx->ESNK,1,0);
    printf1(ctx, "Read %d keys from %s\n",ctx->ESNK,ctx->PMF2dName);

    if (ctx->SILENTFlg < 2)
        printfe(ctx, "Sorting keys ...\n");

    tda_qsort_r((char *)ctx->ESKEY,(size_t)(ctx->ESNK),(size_t)(klen + 1), eselcomp, ctx);

    s = ctx->AcC + low[1];
    wcnt = rcnt = 0;
    while (fgets(ctx->AcC,ctx->PMLEN + 1,ctx->PMFd)) {
           
        prn_message(ctx, ++rcnt,0,0);

        /* get key */

        p = ctx->AcD;
        q = s;
        for (i = 0; i < klen; ++i)
            *p++ = *q++;
        *p = '\0';

        if (esel_c(ctx, klen + 1)) {
            fprintf(ctx->PMF1d,"%s",ctx->AcC);
            wcnt++;
        }
    }
    prn_message(ctx, rcnt,1,0);
    printf1(ctx, "%d records written to: %s\n",wcnt,ctx->PMF1dName);
       
    err = 0;
 
ESELFin:
    if (keya) {
        free((char *)ctx->ESKEY);
        memrq(ctx, -ctx->PMNOC * (klen + 1),sizeof(char));
    }
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  eselcomp()    compare function                                          */

int eselcomp(const void *arg1, const void *arg2, void *_ctx)
{
    (void)_ctx;        /* unused: the signature is shared */
    register char *p,*q;

    p = (char *)arg1;
    q = (char *)arg2;
    return(strcmp(p,q));
}

/* ------------------------------------------------------------------------ */
/*  esel_c()    return 1 if key occurs in ESKEY[]                           */

int esel_c(TDAContext *ctx, int klen)
{     
    int n,low,mid,high;

    low = 1; high = ctx->ESNK;

    while (low <= high) {
        mid = (low + high) / 2;
        n = strcmp(ctx->AcD,ctx->ESKEY + (mid - 1) * klen);
        if (n < 0)
            high = mid - 1;
        else if (n > 0)
            low = mid + 1;
        else {
            if (ctx->PMOPT == 2)
                *(ctx->ESKEY + (mid - 1) * klen) = '\0';
            return(1);
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  emerge()    Command: external merge.                                    */
/*                                                                          */
/*              emerge(                                                     */
/*                  df = ...,               output file, required           */
/*                  mf = ...,               list with merge files           */
/*                  len=...,                max record length, def. 1000    */
/*                  sepc=...,               separation character, def ' '   */
/*                  m=...,                  0,...,9, default blanks         */  
/*              ) = fname;  input file                                      */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int emerge(TDAContext *ctx)
{
    register int i,k;
    int err,rcnt,wcnt,n,nl,kmax,nw,klen;
    register char sc,c,*p,*q,*s,*kp;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Merging files. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 6,2,1))     /* get parameters */
        goto EMFin;

    printf1(ctx, "Inputfile: %s\nOutputfile: ",ctx->PMFdName);
    if (ctx->PMF1Def)
        printf1(ctx, "%s\n",ctx->PMF1dName);
    else {
        printf1(ctx, "need output file.\n");
        goto EMFin;
    }
    if (ctx->MFN == 0) {
        printf1(ctx, "Error: need at least one file for merging.\n");
        goto EMFin;
    }  
    sc = ' ';
    if (ctx->PMMFlg)
        sc = (char)('0' + (ctx->PMM % 10));

    nl = 0;
    for (i = 0; i < ctx->MFN; ++i) {
        k = (int)(strlen(ctx->MFFNAME[i]));
        if (nl < k)
            nl = k;
    }
    klen = ctx->MFJ2[0] - ctx->MFJ1[0] + 1;
    kmax = 0;
    for (i = 0; i < ctx->MFN; ++i) {
        printf1(ctx, "Merging: %s ",ctx->MFFNAME[i]);
        prnchar(ctx, ' ',nl - (int)strlen(ctx->MFFNAME[i]),0);
        printf1(ctx, "[%2d-%2d] with %s [%2d-%2d]\n",
            ctx->MFI1[i],ctx->MFI2[i],ctx->PMFdName,ctx->MFJ1[i],ctx->MFJ2[i]);
        if (ctx->MFI1[i] < 1 || ctx->MFJ1[i] < 1 || klen != (ctx->MFI2[i] - ctx->MFI1[i] + 1)    
                                       || klen != (ctx->MFJ2[i] - ctx->MFJ1[i] + 1)) {
            printf1(ctx, "Error in key definitions.\n");
            goto EMFin;
        }
        if (kmax < ctx->MFI2[i])
            kmax = ctx->MFI2[i];
        if (kmax < ctx->MFJ2[i])
            kmax = ctx->MFJ2[i];
    }
    if (ctx->PMLEN < 1)
        ctx->PMLEN = 1000;

    printf1(ctx, "Max record length: %d\n",ctx->PMLEN);
    if (ctx->PMLEN < kmax) {
        printf1(ctx, "Inconsistent with key definitions.\n");
        goto EMFin;
    }
    if (alloc_acc(ctx, (ctx->PMLEN + 2) * (ctx->MFN + 1)))
        goto EMFin;   

    if (alloc_acd(ctx, (klen + 1) * (ctx->MFN + 1)))
        goto EMFin;   

    if (alloc_acn(ctx, ctx->MFN))
        goto EMFin;   

    if (alloc_acm(ctx, ctx->MFN))
        goto EMFin;   

    if (alloc_aci(ctx, ctx->MFN))
        goto EMFin;   

    if (alloc_acptr(ctx, ctx->MFN))
        goto EMFin;   

    for (i = 0; i < ctx->MFN; ++i) {
        p = ctx->AcPtr[i] = ctx->AcC + (i + 1) * (ctx->PMLEN + 2);
        q = p + ctx->MFI1[i] - 1;
        if (!fgets(p,ctx->PMLEN + 1,ctx->MFFD[i])) {
            printf1(ctx, "Error: can't read from %s.\n",ctx->MFFNAME[i]);
            goto EMFin;
        }
        n = 0;
        while (*p && *p != LF && *p != CR) {
            p++;
            n++;
        }
        *p = '\0';
        ctx->AcN[i] = n;               
        ctx->AcM[i] += 1;

        kp = ctx->AcD + (klen + 1) * (i + 1);
        for (k = 0; k < klen; ++k)
            *kp++ = *q++;
    }

    nw = wcnt = rcnt = 0;
    while (fgets(ctx->AcC,ctx->PMLEN + 1,ctx->PMFd)) {
           
        prn_message(ctx, ++rcnt,0,0);

        p = ctx->AcC;
        n = 0;
        while (*p && *p != LF && *p != CR) {
            p++;
            n++;
        }
        *p = '\0';
        if (n < ctx->MFJ2[0]) {
            nw++;
            continue;
        }
        fprintf(ctx->PMF1d,"%s",ctx->AcC);
        if (ctx->XSEPC)
            fprintf(ctx->PMF1d,"%c",ctx->XSEPC);

        /* get key from input file */

        *(ctx->AcC + ctx->MFJ2[0]) = '\0';
        kp = ctx->AcD;
        q = ctx->AcC + ctx->MFJ1[0] - 1;

        if (rcnt > 1) {             /* check ordering */
            if (strcmp(kp,q) > 0) {
                printf1(ctx, "\nError: %s not sorted in ascending order.\n",ctx->PMFdName);
                printf1(ctx, "Check records %d and %d.\n",rcnt-1,rcnt);
                goto EMFin;
            }
        }
        for (k = 0; k < klen; ++k)
            *kp++ = *q++;
   
        for (i = 0; i < ctx->MFN; ++i) {
            if (ctx->AcN[i] > 0) {
                s = ctx->AcPtr[i];
                p = s + ctx->MFI1[i] - 1;
                q = p + klen;
                c = *q;
                *q = '\0';

                while (strcmp(ctx->AcD,p) > 0) {

                    s = ctx->AcPtr[i];
                    if (!fgets(s,ctx->PMLEN + 1,ctx->MFFD[i])) {
                        ctx->AcN[i] = -ctx->AcN[i];
                        break;
                    }
                    n = 0;
                    while (*s && *s != LF && *s != CR) {
                        s++;
                        n++;
                    }
                    *s = '\0';

                    ctx->AcM[i] += 1;
                    c = *q;
                    *q = '\0';

                    kp = ctx->AcD + (klen + 1) * (i + 1);
                    if (strcmp(kp,p) > 0) {
                        printf1(ctx, "\nError: %s not sorted in ascending order.\n",ctx->MFFNAME[i]);
                        printf1(ctx, "Check records %d and %d.\n",ctx->AcM[i]-1,ctx->AcM[i]);
                        goto EMFin;
                    }
                    for (k = 0; k < klen; ++k)
                        *kp++ = *p++;
                    p = ctx->AcPtr[i] + ctx->MFI1[i] - 1;
                }
                n = strcmp(ctx->AcD,p);          
                *q = c;
                if (n == 0) {
                    fprintf(ctx->PMF1d,"%s",ctx->AcPtr[i]);
                    ctx->AcI[i] += 1;
                }
                else
                    fprnchar(ctx, ctx->PMF1d,(unsigned char)(sc),iabs(ctx, ctx->AcN[i]),0);
            }
            else  
                fprnchar(ctx, ctx->PMF1d,(unsigned char)(sc),-ctx->AcN[i],0);

            if (ctx->XSEPC)
                fprintf(ctx->PMF1d,"%c",ctx->XSEPC);
        }
        fprintf(ctx->PMF1d,"\n");
        wcnt++;
    }
    prn_message(ctx, rcnt,1,0);

    printf1(ctx, "\nRead %8d records from: %s\n",rcnt,ctx->PMFdName);
    for (i = 0; i < ctx->MFN; ++i) {
        printf1(ctx, "Read %8d records from: %s ",ctx->AcM[i],ctx->MFFNAME[i]);
        prnchar(ctx, ' ',nl - (int)strlen(ctx->MFFNAME[i]),0);
        printf1(ctx, " matched: %d\n",ctx->AcI[i]);
    }

    printf1(ctx, "\n%d records written to: %s\n",wcnt,ctx->PMF1dName);

    if (nw > 0) {
        printf1(ctx, "\nWarning: skipped %d records from %s\n",nw,ctx->PMFdName);
        printf1(ctx, "         that have record length less than %d.\n",ctx->MFJ2[0]);
    }
    err = 0;
 
EMFin:
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  ejoin()     Command: join two event data files.                         */
/*                                                                          */
/*              ejoin(                                                      */
/*                  if1= ...,               first file                      */
/*                  if2= ...,               second file                     */
/*                  len=...,                max record length, def. 1000    */
/*                  max=...,                max block size, def. 1000       */
/*                  nw=...,                 max number of levels, def. 1    */
/*                  noc-...,                read max noc records, def. all  */
/*                  fmt0=...,               print formats                   */
/*                  fmt1=...,               print formats                   */
/*                  fmt2=...,               print formats                   */
/*              ) = fname;  output file                                     */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int ejoin(TDAContext *ctx)
{
    int err,r,nv1,nv2,nx1,nx2,id1,id2,r1,r2;
    int wrec,idn1,idn2,njoin,l1,l2,rcnt;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Join two event data files. Current memory: %d bytes.\n\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,1,1))     /* get parameters */
        goto EJFin;

    printf1(ctx, "Outputfile : %s\n",ctx->PMFdName);

    if (ctx->PMMax < 1)
        ctx->PMMax = 1000;

    if (ctx->PMLEN < 1)
        ctx->PMLEN = 1000;

    printf1(ctx, "Max block size: %d\n",ctx->PMMax);
    printf1(ctx, "Max number of levels: %d\n",ctx->PMNW);
    printf1(ctx, "Max record length: %d\n",ctx->PMLEN);

    if (alloc_acs(ctx, 2 * ctx->PMMax + 1))
        goto EJFin;   
    if (alloc_aci(ctx, 2 * ctx->PMMax * ctx->PMNW + 1))
        goto EJFin;   
    if (alloc_acj(ctx, 2 * ctx->PMMax * ctx->PMNW + 1))
        goto EJFin;   

    njoin = idn1 = idn2 = nv1 = nv2 = nx1 = nx2 = 0;
    wrec = ctx->EJRECN1 = ctx->EJRECN2 = 0;
    ctx->EJEOF1 = ctx->EJEOF2 = 1;

    if (ctx->PMIF1Def) {
        if (alloc_acc(ctx, ctx->PMLEN + 2))
            goto EJFin;   
        nv1 = ejoin_nv(ctx, ctx->PMIF1d,ctx->AcC,&ctx->EJRECN1);
        printf1(ctx, "Inputfile 1: %s  (%d variables)\n",ctx->PMIF1Name,nv1);
        if (nv1 < 6) {
            printf1(ctx, "Error: need at least 6 variables.\n");
            goto EJFin;
        }
        nx1 = nv1 - 6;
        ctx->EJEOF1 = 0;
        if (alloc_acn(ctx, 6 * ctx->PMMax))
            goto EJFin;   
        if (alloc_acx(ctx, nx1 * ctx->PMMax + 1))
            goto EJFin;   
    }
    if (ctx->PMIF2Def) {
        if (alloc_acd(ctx, ctx->PMLEN + 2))
            goto EJFin;   
        nv2 = ejoin_nv(ctx, ctx->PMIF2d,ctx->AcD,&ctx->EJRECN2);
        printf1(ctx, "Inputfile 2: %s  (%d variables)\n",ctx->PMIF2Name,nv2);
        if (nv2 < 6) {
            printf1(ctx, "Error: need at least 6 variables.\n");
            goto EJFin;
        }
        nx2 = nv2 - 6;
        ctx->EJEOF2 = 0;
        if (alloc_acm(ctx, 6 * ctx->PMMax))
            goto EJFin;   
        if (alloc_acy(ctx, nx2 * ctx->PMMax + 1))
            goto EJFin;   
    }

    while (ctx->EJEOF1 == 0 || ctx->EJEOF2 == 0) {

        if (ctx->EJEOF1 == 0) {
            r1 = ejoin_rd(ctx, ctx->PMIF1d,ctx->PMIF1Name,ctx->AcC,&ctx->EJEOF1,&ctx->EJRECN1,ctx->AcN,ctx->AcX,nx1,&l1);
            if (r1 < 0)
                goto EJFin;
            id1 = ctx->AcN[0];
            idn1++;
        }
        else {
            id1 = ctx->INTMAX;
            r1 = 0;
            l1 = 1;
        }
        if (ctx->EJEOF2 == 0) {
            r2 = ejoin_rd(ctx, ctx->PMIF2d,ctx->PMIF2Name,ctx->AcD,&ctx->EJEOF2,&ctx->EJRECN2,ctx->AcM,ctx->AcY,nx2,&l2);
            if (r2 < 0)
                goto EJFin;
            id2 = ctx->AcM[0];
            idn2++;
        }
        else {
            id2 = ctx->INTMAX;
            r2 = 0;
            l2 = 1;
        }

EJCONT:
        while (id1 < id2) {
            if ((r = ejoin_proc(ctx, id1,r1,0,l1,1,nx1,nx2)) < 0)
                goto EJFin;
            wrec += r;

            rcnt = imax(ctx, ctx->EJRECN1,ctx->EJRECN2);
            prn_message(ctx, rcnt,0,0);
            if (ctx->PMNOCFlg && rcnt >= ctx->PMNOC)
                goto EJFin1;

            if (ctx->EJEOF1 == 0) {
                r1 = ejoin_rd(ctx, ctx->PMIF1d,ctx->PMIF1Name,ctx->AcC,&ctx->EJEOF1,&ctx->EJRECN1,ctx->AcN,ctx->AcX,nx1,&l1);
                if (r1 < 0)
                    goto EJFin;
                id1 = ctx->AcN[0];
                idn1++;
            }
            else {
                id1 = ctx->INTMAX;
                r1 = 0;
                l1 = 1;
            }
        }

        while (id2 < id1) {
            if ((r = ejoin_proc(ctx, id2,0,r2,1,l2,nx1,nx2)) < 0)
                goto EJFin;
            wrec += r;

            rcnt = imax(ctx, ctx->EJRECN1,ctx->EJRECN2);
            prn_message(ctx, rcnt,0,0);
            if (ctx->PMNOCFlg && rcnt >= ctx->PMNOC)
                goto EJFin1;

            if (ctx->EJEOF2 == 0) {
                r2 = ejoin_rd(ctx, ctx->PMIF2d,ctx->PMIF2Name,ctx->AcD,&ctx->EJEOF2,&ctx->EJRECN2,ctx->AcM,ctx->AcY,nx2,&l2);
                if (r2 < 0)
                    goto EJFin;
                id2 = ctx->AcM[0];
                idn2++;
            }
            else {
                id2 = ctx->INTMAX;
                r2 = 0;
                l2 = 1;
            }
        }
        if (id1 < id2)
            goto EJCONT;

        if (id1 == id2 && id1 < ctx->INTMAX) {

            if ((r = ejoin_proc(ctx, id1,r1,r2,l1,l2,nx1,nx2)) < 0)
                goto EJFin;
            wrec += r;
            njoin++;
        }
        rcnt = imax(ctx, ctx->EJRECN1,ctx->EJRECN2);
        prn_message(ctx, rcnt,0,0);
        if (ctx->PMNOCFlg && rcnt >= ctx->PMNOC)
            break;
    }

EJFin1:
    newline(ctx);
    if (ctx->PMIF1Def)
        printf1(ctx, "Read %d records (%d blocks) from: %s\n",ctx->EJRECN1,idn1,ctx->PMIF1Name);
    if (ctx->PMIF2Def)
        printf1(ctx, "Read %d records (%d blocks) from: %s\n",ctx->EJRECN2,idn2,ctx->PMIF2Name);
    if (ctx->PMIF1Def && ctx->PMIF2Def)
        printf1(ctx, "Number of common blocks: %d\n",njoin);
    printf1(ctx, "%d records written to %s\n",wrec,ctx->PMFdName);

    err = 0;
 
EJFin:
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  ejoin_nv()      return number of variables found in first record        */
/*                  of file fn.                                             */

int ejoin_nv(TDAContext *ctx, FILE *fd,char *buf,int *nrec)
{
    int nv;
    double x;
    register char *p;

    nv = 0;
    while (fgets(buf,ctx->PMLEN + 1,fd)) {
        *nrec += 1;
        if (check_drec(ctx, buf)) {      /* check for data records */
            p = buf;
            while (*p) {
                p = skip_b(ctx, p);
                if (sscanf(p,"%lf",&x) != 1)
                    break;
                nv++;
                p = skip_dbl(ctx, p);
            }
            break;
        }
    }
    return(nv);
}

/* -##--------------------------------------------------------------------- */
/*  ejoin_rd(fd,fname,buf,eof,nrec,bvar,xvar,nx,mxlev)                      */
/*                                                                          */
/*  Read one block of fd. Return number of records, or 0 if EOF.            */
/*  return -1 if number of records in block exceeds PMMax.                  */
/*  return -2 if less than nx variables.                                    */
/*  return -3 if max number of levels exceeded.                             */
/*  return max number of levels in mxlev                                    */

int ejoin_rd(TDAContext *ctx, FILE *fd,char *fname,char *buf,int *eof,int *nrec, int *bvar,double *xvar,int nx,int *mxlev)
{
    register int i,j,l;
    int n,nv,id;
    double x;
    register char *p;

    *mxlev = 0;

    if (*eof)
        return(0);

    n = nv = 0;

    while (1) {

        if (n >= ctx->PMMax) {
            printf1(ctx, "Error: exceeded max block size in file %s\n",fname);
            return(-1);
        }
        p = buf;
        nv = 0;

        while (*p) {
            p = skip_b(ctx, p);
            if (sscanf(p,"%lf",&x) != 1)
                break;
            if (nv < 6)  
                bvar[n * 6 + nv] = (int)x;
            else if (nv - 6 < nx)
                xvar[n * nx + nv - 6] = x;

            if (nv == 0 && n > 0) {
                if ((int)x != id)  
                    goto EJRDFin;
            }
            nv++;
            p = skip_dbl(ctx, p);
        }
        if (nv < nx + 6) {
            printf1(ctx, "Error: found less than %d variables in file %s.\n",nx + 6,fname);
            return(-2);
        }
        if (n == 0)
            id = bvar[0];
        n++;

        *eof = 1;

        while (fgets(buf,ctx->PMLEN + 1,fd)) {
            *nrec += 1;
            if (!check_drec(ctx, buf))        /* check for data records */
                continue;
            *eof = 0;
            break;
        }
        if (*eof)
            break;
    }

EJRDFin:
    j = 0;
    bvar[1] = 0;
    while (j + 1 < n) {
        l = 0;
        for (i = j + 1; i < n; ++i) {
            if (bvar[i * 6 + 3] < bvar[j * 6 + 4])  
                l++;
            else  
                l = 0;
            bvar[i * 6 + 1] = l;
            *mxlev = imax(ctx, *mxlev,l);
            if (l == 0)  
                break;
        }
        j++;
    }
    *mxlev += 1;

    if (*mxlev > ctx->PMNW) {
        printf1(ctx, "Error: exceeded max number of levels in file %s.\n",fname);
        return(-3);
    }
    return(n);
}

/* -##--------------------------------------------------------------------- */
/*  ejoin_proc()   Process one or two blocks. Return number of records      */
/*                 written to output file, or -1 if error.                  */

int ejoin_proc(TDAContext *ctx, int id,int r1,int r2,int l1,int l2,int nx1,int nx2)
{
    register int i,j,k,l;
    int n,nn,nnn,nfmt,t,ii,iii,ll,jj,wrec,ts,tf;

    n = wrec = 0;

    for (i = 0; i < r1; ++i) {
        ctx->AcS[n++] = ctx->AcN[i * 6 + 3];
        ctx->AcS[n++] = ctx->AcN[i * 6 + 4];
    }
    for (i = 0; i < r2; ++i) {
        ctx->AcS[n++] = ctx->AcM[i * 6 + 3];
        ctx->AcS[n++] = ctx->AcM[i * 6 + 4];
    }
    if (sorti(ctx, n,ctx->AcS,0))
        return(-1);  

    nn = 0;
    for (i = 1; i < n; ++i) {
        if (ctx->AcS[nn] < ctx->AcS[i]) 
            ctx->AcS[++nn] = ctx->AcS[i];
    }
    nn++;

    i = nn * imax(ctx, l1,l2);  
    for (k = 0; k < i; ++k)
        ctx->AcI[k] = ctx->AcJ[k] = ctx->EJMVAL;

    k = 0;
    t = ctx->AcS[0];
    for (i = 0; i < r1; ++i) {

        ts = ctx->AcN[i * 6 + 3];
        tf = ctx->AcN[i * 6 + 4];

        while (t < ts)  
            t = ctx->AcS[++k];

        while (t < tf) {
            l = 0;
            ctx->AcI[k] = i;

            ii = i;
            while (++ii < r1) {
                if (ctx->AcN[ii * 6 + 3] <= t) {
                    if (t < ctx->AcN[ii * 6 + 4]) {
                        l++;
                        ctx->AcI[l * nn + k] = ii;
                    }
                }   
                else
                    break;
            }
            t = ctx->AcS[++k];
        }
    }
    k = 0;
    t = ctx->AcS[0];
    for (i = 0; i < r2; ++i) {

        ts = ctx->AcM[i * 6 + 3];
        tf = ctx->AcM[i * 6 + 4];

        while (t < ts)  
            t = ctx->AcS[++k];

        while (t < tf) {
            l = 0;
            ctx->AcJ[k] = i;
            ii = i;
            while (++ii < r2) {
                if (ctx->AcM[ii * 6 + 3] <= t) {
                    if (t < ctx->AcM[ii * 6 + 4]) {
                        l++;
                        ctx->AcJ[l * nn + k] = ii;
                    }
                }   
                else
                    break;
            }
            t = ctx->AcS[++k];
        }
    }
    nnn = 0;
    for (k = 0; k < nn - 1; ++k) {
        ii = 0;
        for (l = 0; l < l1; ++l) {
            if (ctx->AcI[l * nn + k] < 0)
                break;
            ii = l;
        }
        iii = 0;
        for (l = 0; l < l2; ++l) {
            if (ctx->AcJ[l * nn + k] < 0)
                break;
            iii = l;
        }
        nnn += imax(ctx, ++ii,++iii);
    }
    for (k = 0; k < nn - 1; ++k) {
        ii = 0;
        for (l = 0; l < l1; ++l) {
            if (ctx->AcI[l * nn + k] < 0)
                break;
            ii = l;
        }
        iii = 0;
        for (l = 0; l < l2; ++l) {
            if (ctx->AcJ[l * nn + k] < 0)
                break;
            iii = l;
        }
        i = imax(ctx, ++ii,++iii);
  
        for (j = 0; j < i; ++j) {
            ejoin_prn(ctx, 0,id);
            ejoin_prn(ctx, 1,nnn);
            ejoin_prn(ctx, 2,k + 1);
            ejoin_prn(ctx, 3,j);
            ejoin_prn(ctx, 4,ctx->AcS[k]);
            ejoin_prn(ctx, 5,ctx->AcS[k + 1]);
            nfmt = 6;
            if (ctx->PMIF1Def) {
                l = ctx->AcI[j * nn + k];
                if (l >= 0)
                    ejoin_prn(ctx, nfmt++,ctx->AcN[l * 6 + 5]);
                else
                    ejoin_prn(ctx, nfmt++,ctx->EJMVAL);
            }
            if (ctx->PMIF2Def) {
                l = ctx->AcJ[j * nn + k];
                if (l >= 0)
                    ejoin_prn(ctx, nfmt++,ctx->AcM[l * 6 + 5]);
                else
                    ejoin_prn(ctx, nfmt++,ctx->EJMVAL);
            }
            if (ctx->PMIF1Def && nx1 > 0) {
                l = ctx->AcI[j * nn + k];
                if (l >= 0) {
                    ll = l * nx1;
                    for (jj = 0; jj < nx1; ++jj)  
                        ejoin_prn1(ctx, jj,ctx->AcX[ll++]);
                }
                else {
                    for (jj = 0; jj < nx1; ++jj)  
                        ejoin_prn1(ctx, jj,(double)ctx->EJMVAL);
                }
            }
            if (ctx->PMIF2Def && nx2 > 0) {
                l = ctx->AcJ[j * nn + k];
                if (l >= 0) {
                    ll = l * nx2;
                    for (jj = 0; jj < nx2; ++jj)  
                        ejoin_prn2(ctx, jj,ctx->AcY[ll++]);
                }
                else {
                    for (jj = 0; jj < nx2; ++jj)  
                        ejoin_prn2(ctx, jj,(double)ctx->EJMVAL);
                }
            }
            fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
            tda_export_endrow(ctx, "ejoin.table");
#endif
            wrec++;
        }
    }
    return(wrec);
}

/* ------------------------------------------------------------------------ */
/*  ejoin_prn()   print a number                                            */

void ejoin_prn(TDAContext *ctx, int nfmt,int val)
{
    int n;

    n = ctx->PMXFmtN[0];    
    if (n > 0) {            
        n = imin(ctx, nfmt,n - 1);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMXFmtS[0] + n * ctx->PMXFmtLen,(double)val);
    }
    else
        fprintf(ctx->PMFd,"%d ",val);
#ifdef TDA_R_PACKAGE
    /* every value ejoin writes goes through here; the row is closed by
       the newline the caller writes (see ejoin_prn1 below) */
    tda_export_cell(ctx, "ejoin.table", (double)val);
#endif
}

/* ------------------------------------------------------------------------ */
/*  ejoin_prn1()   print a number                                           */

void ejoin_prn1(TDAContext *ctx, int nfmt,double val)
{
    int n;

    n = ctx->PMXFmtN[1];    
    if (n > 0) {            
        n = imin(ctx, nfmt,n - 1);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMXFmtS[1] + n * ctx->PMXFmtLen,val);
    }
    else
        fprintf(ctx->PMFd,"%g ",val);
#ifdef TDA_R_PACKAGE
    tda_export_cell(ctx, "ejoin.table", val);
#endif
}

/* ------------------------------------------------------------------------ */
/*  ejoin_prn2()   print a number                                           */

void ejoin_prn2(TDAContext *ctx, int nfmt,double val)
{
    int n;

    n = ctx->PMXFmtN[2];    
    if (n > 0) {            
        n = imin(ctx, nfmt,n - 1);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMXFmtS[2] + n * ctx->PMXFmtLen,val);
    }
    else
        fprintf(ctx->PMFd,"%g ",val);
#ifdef TDA_R_PACKAGE
    tda_export_cell(ctx, "ejoin.table", val);
#endif
}



