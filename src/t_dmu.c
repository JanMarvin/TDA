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

/* ------------------------------------------------------------------------ */
/*  functions in tda_dm.                                                    */

int dm_ccnt(void);
int dm_lcnt(void);
int dm_dump(void);
void prnbuf(unsigned char *buf, int n, int off);
int dm_dsplit(void);
int esort(void);
int escomp(const void *arg1,const void *arg2);
long get_long(char *p);
void put_long(char *p,long n);
int eskip(void);
int eselect(void);
int eselcomp(const void *arg1,const void *arg2);
int esel_c(int klen);
int emerge(void);
int ejoin(void);
int ejoin_nv(FILE *fd,char *buf,int *nrec);
int ejoin_rd(FILE *fd,char *fname,char *buf,int *eof,int *nrec,
    int *bvar,double *xvar,int nx,int *mxlev);
int ejoin_proc(int id,int r1,int r2,int l1,int l2,int nx1,int nx2);
void ejoin_prn(int nfmt,int val);
void ejoin_prn1(int nfmt,double val);
void ejoin_prn2(int nfmt,double val);

/*--------------------------------------------------------------------------*/
/*  Global variables                                                        */

int DMNOC = 0;              /* max number of records                        */
int DMFFlg = 0;             /* Set for dmfile command                       */
char DMFNam[FNMaxLen];      /* File defined by dmfile command               */
FILE *DMFd;                 /* File descriptor for DMFNam                   */
char *ESKEY;                /* array with keys, used in esort()             */
int ESNK = 0;               /* number of keys                               */
int ESLEN[5];               /* length of keys                               */

int EJEOF1 = 0;             /* set if ejoin file 1 at eof                   */
int EJEOF2 = 0;             /* set if ejoin file 2 at eof                   */
int EJRECN1 = 0;
int EJRECN2 = 0;
int EJMVAL = -3;

/* ------------------------------------------------------------------------ */
/*  dm_ccnt()   Command: ccnt=fname. Print frequency distribution of        */
/*              characters in fname to standard output.                     */
/*              Return 0 if OK, -1 if error.                                */

int dm_ccnt(void)
{
    register int i,j;
    register char *p;
    int err,nn,cnt,rcnt;

    err = -1;
    nn = 20000;         /* size of read buffer */

    if (check_cmd(1))
        return(-1);

    if (parm(CmdBuf + 4,10,1)) {   /* get parameters */
        p_clean();
        return(-1);   
    }
    printf1("Frequency distribution of characters in file: %s\n",PMFdName);

    if (alloc_acn(257))
        goto CCNTFin;
   
    if (alloc_acc(nn + 1))
        goto CCNTFin;
   
    rcnt = 0;
    while ((cnt = fread(AcC,sizeof(char),nn,PMFd)) > 0) {

        rcnt += cnt;
        p = AcC;
        for (i = 0; i < cnt; ++i) {
            j = (int) *p++;
            if (j < 0)
                j += 256;
            AcN[j] += 1;
        }
    }
    printf1("\nCharacter   Hexadecimal      Count\n");
    prnchar('-',34,1);

    for (i = 0; i < 256; ++i) {
        if (AcN[i]) {
            if (isprint(i))  
                printf1("    %c",(char)i);
            else  
                printf1("     ");
            printf1("            %02X     %10d\n",i,AcN[i]);
        }
    }
    prnchar('-',34,1);
    printf1("Sum                     %10d\n",rcnt);
    err = 0;

CCNTFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  dm_lcnt()   Command: lcnt(noc=)=fname. Print frequency distribution of  */
/*              length of records to standard output.                       */
/*              Return 0 if OK, -1 if error.                                */
/*                                                                          */
/*  Note: record length is printed without EOL!                             */  

int dm_lcnt(void)
{
    register int i,j,l;
    register char *p;
    int err,nn,n,rec,wflg,wflg1;
    double d1,d2;

    err = -1;
    wflg = wflg1 = 0;
    nn = 20000;         /* size of read buffer */

    if (check_cmd(1))
        return(-1);

    if (parm(CmdBuf + 4,2,1)) {   /* get parameters */
        p_clean();
        return(-1);   
    }
    printf1("Frequency distribution of record lengths (without EOL) in file: %s\n",PMFdName);
    printf1("Max number of records: %d\n",PMNOC);

    if (alloc_acc(nn + 2))
        goto LCNTFin;

    if (alloc_acn(PMNOC + 1))
        goto LCNTFin;

    rec = 0;
    while (fgets(AcC,nn + 1,PMFd)) {
        if (rec >= PMNOC) {
            wflg = 1;
            break;
        }
        l = strlen(AcC);
        p = AcC + l;
        while (p > AcC && (*--p == LF || *p == CR))
            l--;

        AcN[++rec] = l;
        if (l >= nn)
            wflg1 = 1;
    }
    printf1("Read %d records\n",rec);
    if (wflg)  
        printf1("Warning: file has more than %d records\n",rec);
    if (wflg1)  
        printf1("Warning: record length may have exceeded the maximum of %d bytes.\n",nn);

    if (rec <= 0)  
        goto LCNTFin;

    if (alloc_aci(rec + 1))
        goto LCNTFin;

    if (alloc_acj(rec + 1))
        goto LCNTFin;

    if (alloc_acm(rec + 1))
        goto LCNTFin;

    n = cfreq(1,rec,AcN,rec,AcM,AcI,AcJ,0,&d1,&d2);
    if (n <= 0)  
        printf1("Error: can't sort length of records.\n");
    else {
        printf1("Index  Length  Frequency\n");
        for (i = 1; i <= n; ++i) {
            j = AcJ[i];
            printf1("%5d %7d %10d\n",i,AcM[j],AcI[j]);
        }
    }
    err = 0;

LCNTFin:
    p_clean();    
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  dm_dump     Command: dump(nc=,s=)=fname.                                */
/*              hex dump of fname, nc is number of characters, s is offset  */
/*              Return 0 if OK, -1 if error.                                */

int dm_dump(void)
{
    int err,nn,n,off,rcnt,pcnt;
        
    err = -1;
    nn = 32000;         /* size of read buffer, should be multiple of 32 */

    if (check_cmd(1))
        return(-1);
   
    if (parm(CmdBuf + 4,10,1)) {   /* get parameters */
        p_clean();
        return(-1);   
    }
    if (PMS < 0)
        PMS = 0;

    printf1("Hex dump of file: %s [offset=%d]\n",PMFdName,PMS);
    if (PMS > 0) {
        if (fseek(PMFd,(long)PMS,0)) {
            printf1("Can't seek to offset %d.\n",PMS);
            goto DUMPFin;
        }
        off = PMS;
    }
    else
        off = 0;

    if (alloc_acc(nn + 2))
        goto DUMPFin;

    rcnt = nn;         
    pcnt = 0;

    while (1) {

        if (PMNC > 0 && rcnt > PMNC)
            rcnt = PMNC;

        if ((n = fread(AcC,sizeof(char),rcnt,PMFd)) < 0) {
            printf1("Error in reading file: %s\n",PMFdName);
            goto DUMPFin;
        }
        if (!n)       
            break;

        prnbuf((unsigned char *)AcC,n,off);
        off += n; 

        pcnt += n;
        if (PMNC > 0 && pcnt >= PMNC)
            break;
    }
    printf1("\n");
    err = 0;

DUMPFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prnbuf(buf,n,off)                                                       */
/*      print n char of buf in hex and ASCII to stdout. off is used to      */
/*      count the header column.                                            */

void prnbuf(unsigned char *buf, int n, int off)
{
    register int ii,i;
    register unsigned char c;
    int llen = 16;                  /* characters per line */

    i = 0;
    while (1) {
        printf1("\n%08x : ",off + i);
        for (ii = 0; ii < llen; ++ii) {
            if ((i + ii) >= n)
                break;
            printf1("%02x ",buf[i + ii]);
        }
        for ( ; ii < llen; ++ii)  
            printf1("   ");
        printf1("  ");

        for (ii = 0; ii < llen; ++ii) {
            if ((i + ii) >= n)
                break;
            c = buf[i + ii];
            if (c >= 0x20 && c < 0x7f)  
                printf1("%c",c);
            else  
                printf1(".");
        }
        if ((i += llen) >= n)
            break;
    }
}

/* ------------------------------------------------------------------------ */
/*  dm_dsplit.  Command: dsplit(len=)=fname. Split file into parts.         */
/*              Return 0 if OK, -1 if error.                                */

int dm_dsplit(void)
{
    FILE *fd;
    int err,cnt,rcnt;
    char *p,fname[1000];

    err = -1;
    if (check_cmd(1))
        return(-1);

    if (parm(CmdBuf + 6,10,1)) {   /* get parameters */
        p_clean();
        return(-1);   
    }
    if (PMLEN < 1)
        PMLEN = 1000;

    printf1("Split file %s into parts of length %d.\n",PMFdName,PMLEN);

    if (alloc_acc(PMLEN + 2))
        goto DSPLFin;

    strcpy(fname,PMFdName);
    strcat(fname,".a");
    p = fname + strlen(fname) - 1;

    rcnt = 0;
    while ((cnt = fread(AcC,sizeof(char),PMLEN,PMFd)) > 0) {
        rcnt += cnt;
        printf1("%s -- ",fname);

        if (!(fd = fopen(fname,OPEN_WB))) {
            printf1("can't create this output file\n");
            break;         
        }
        if (fwrite(AcC,sizeof(char),cnt,fd) != cnt) {
            printf1("write error\n");
            break;
        }
        printf1("written %d bytes\n",cnt);
        fclose(fd);
        if (*p == 'z') {
            p++;
            *p++ = 'a';
            *p-- = '\0';
        }
        else
            *p += 1;
    }
    printf1("\nRead from input file %s: %d bytes.\n",PMFdName,rcnt);
    err = 0;
 
DSPLFin:
    p_clean();
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

int esort(void)
{
    FILE *ofd[ESMAXNF];
    register int i,j,k;
    int err,low[5],up[5],klen,kmax,keya,rcnt,tcnt,wcnt,nf,nn;
    long fptr,wptr;
    register char *p,*q,*s;
    char *kp,fname[500];

    keya = 0;
    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("External sort. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,2,1))     /* get parameters */
        goto ESFin;

    printf1("Inputfile: %s\nOutputfile: ",PMFdName);
    if (PMF1Def)
        printf1("%s\n",PMF1dName);
    else {
        printf1("need output file.\n");
        goto ESFin;
    }

    printf1("Sort keys:");
    klen = kmax = ESNK = 0;
    for (i = 1; i < 10; i += 2) {
        if (PMSK[i] > 0 && PMSK[i + 1] > 0) {
            if (PMSK[i + 1] < PMSK[i]) {
                printf1("error in key definition.\n");
                goto ESFin;
            }
            low[ESNK] = PMSK[i] - 1;
            up [ESNK] = PMSK[i + 1];
            printf1(" [%d,%d]",low[ESNK] + 1,up[ESNK]);
            ESLEN[ESNK] = up[ESNK] - low[ESNK];
            klen += ESLEN[ESNK] + 1;
            if (kmax < up[ESNK])
                kmax = up[ESNK];

            ESNK++;
        }
    }
    if (ESNK == 0) {
        printf1(" need a valid definition with sk parameter.\n");
        goto ESFin;
    }
    printf1("\nMax number of records: %d\n",PMNOC);

    if (PMLEN < 1)
        PMLEN = 1000;

    printf1("Max record length: %d\n",PMLEN);
    if (PMLEN < kmax) {
        printf1("Inconsistent with key definition.\n");
        goto ESFin;
    }
    if (PMNOC < ESMAXNF)
        PMNOC = ESMAXNF;

    klen += sizeof(long);
    if (!(ESKEY  = (char *)calloc(PMNOC * klen,sizeof(char)))) {
        p_err(-2,1);             
        goto ESFin;    
    }         
    keya = 1;
    memrq(PMNOC * klen,sizeof(char));

    if (alloc_acc(PMLEN + 2))
        goto ESFin;   

    nf = tcnt = 0;
    fptr = ftell(PMFd);

    while (1) {

        if (fseek(PMFd,fptr,0)) {
            printf1("Error: can't seek to offset %d.\n",fptr);
            goto ESFin;
        }

        rcnt = 0;
        while (fgets(AcC,PMLEN + 1,PMFd)) {
           
            /* save key */

            k = strlen(AcC);
            p = AcC + k;
            for (i = k; i < kmax; ++i)
                *p++ = ' ';

        	   p = ESKEY + rcnt * klen;                 
            put_long(p,fptr);
            p += sizeof(long);     
	    
            for (k = 0; k < ESNK; ++k) {
                q = AcC + low[k];
                for (i = 0; i < ESLEN[k]; ++i)
                    *p++ = *q++;
                *p++ = '\0';
            }
       	    fptr = ftell(PMFd);
            tcnt++;
            prn_message(tcnt,0,0);
	  
            if (++rcnt >= PMNOC)
                break;
        }
        prn_message(tcnt,1,0);
         
        if (rcnt == 0)
            break;

        if (nf >= ESMAXNF) {
            printf1("\nError: exceeded max number of temporary files.\n");
            goto ESFin;
        }
        strcpy(fname,PMF1dName);
        sprintf(fname + strlen(fname),".%d",nf);

        if (!(ofd[nf] = fopen(fname,OPEN_WR))) {
            printf1("\nError: can't create a new temporary output file (%s)\n",fname);
            goto ESFin;    
        }
        if (SILENTFlg < 2)
            printfe("Sorting: %s\n",fname);

        qsort((char *)ESKEY,rcnt,klen,escomp);

        /* write to output file */

        if (SILENTFlg < 2)
            printfe("Writing: %s\n",fname);

        for (i = 0; i < rcnt; ++i) {

            wptr = get_long(ESKEY + i * klen);

            if (fseek(PMFd,wptr,0)) {
                printf1("Error: can't seek to offset %d.\n",wptr);
                goto ESFin;
            }
            if (!fgets(AcC,PMLEN + 1,PMFd)) {
                printf1("Error: can't read at offset %d.\n",wptr);
                goto ESFin;
            }
            fprintf(ofd[nf],"%s",AcC);
        }
        fclose(ofd[nf]);
        nf++;

    }
    printf1("\nRead %d records.\n",tcnt);
    printf1("Number of temporary files: %d\n",nf);
          
    /* MERGING */

    if (SILENTFlg < 2)
        printfe("\nMerging ...\n");

    if (alloc_acc(nf * (PMLEN + 2)))
        goto ESFin;   

    s = AcC;
    for (i = 0; i < nf; ++i) {
        strcpy(fname,PMF1dName);
        sprintf(fname + strlen(fname),".%d",i);
        if (!(ofd[i] = fopen(fname,OPEN_RD))) {
            printf1("\nError: can't open temporary file (%s).\n",fname);
            goto ESFin;    
        }
        if (!fgets(s,PMLEN + 1,ofd[i])) {
            printf1("Error: can't read temp file %s.\n",fname);
            goto ESFin;
        }
        p = ESKEY + i * klen;                 
        put_long(p,(long)i);
        p += sizeof(long);

        for (k = 0; k < ESNK; ++k) {
            q = s + low[k];
            for (j = 0; j < ESLEN[k]; ++j)
                *p++ = *q++;
            *p++ = '\0';
        }
        s += PMLEN + 2;
    }
    kp = ESKEY;
    nn = nf;
    wcnt = 0;
    while (nn > 0) {

        if (nn > 1)
            qsort((char *)kp,nn,klen,escomp);
    
        wptr = get_long(kp);

        j = (int)wptr;
        s = AcC + (PMLEN + 2) * j;
        fprintf(PMF1d,"%s",s);
        wcnt++;
        prn_message(wcnt,0,1);

        if (!fgets(s,PMLEN + 1,ofd[j])) {
            *s = '\0';
            nn--;
            kp += klen;      
        }   
        else {
            p = kp + sizeof(long);

            for (k = 0; k < ESNK; ++k) {
                q = s + low[k];
                for (j = 0; j < ESLEN[k]; ++j)
                    *p++ = *q++;
                *p++ = '\0';
            }
        }
    }
    for (i = 0; i < nf; ++i)
        fclose(ofd[i]);

    prn_message(wcnt,1,1);
    printf1("%d records written to: %s\n",wcnt,PMF1dName);
       
    err = 0;
 
ESFin:
    if (keya) {
        free((char *)ESKEY);
        memrq(-PMNOC * klen,sizeof(char));
    }
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  escomp()    compare function                                            */

int escomp(const void *arg1,const void *arg2)
{     
    register int k,n;
    register char *p,*q;

    p = (char *)arg1 + sizeof(long);
    q = (char *)arg2 + sizeof(long);

    for (k = 0; k < ESNK; ++k) {
        n = strcmp(p,q);
        if (n)
            return(n);
        p += ESLEN[k] + 1;
        q += ESLEN[k] + 1;
    }
    return(n);
}

/* ------------------------------------------------------------------------ */
/*  get_long()                                                              */

long get_long(char *p)
{     
    register int i;
    char lc[sizeof(long)];
    register char *q;

    q = lc;
    for (i = 0; i < sizeof(long); ++i)
        *q++ = *p++;
    return(*(long *)lc);
}

/* ------------------------------------------------------------------------ */
/*  put_long()                                                              */

void put_long(char *p,long n)
{     
    register int i;
    char lc[sizeof(long)];
    register char *q;

    *(long *)lc = n;

    q = lc;
    for (i = 0; i < sizeof(long); ++i)
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

int eskip(void)
{
    register int i,k;
    int err,low[5],up[5],kmax,rcnt;
    register char *p;

    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("Skip columns. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,2,1))     /* get parameters */
        goto ESKFin;

    printf1("Inputfile: %s\nOutputfile: ",PMFdName);
    if (PMF1Def)
        printf1("%s\n",PMF1dName);
    else {
        printf1("need output file.\n");
        goto ESKFin;
    }
    printf1("Column ranges:");
    kmax = ESNK = 0;
    for (i = 1; i < 10; i += 2) {
        if (PMSK[i] > 0 && PMSK[i + 1] > 0) {
            if (PMSK[i + 1] < PMSK[i]) {
                printf1("error in sk definition.\n");
                goto ESKFin;
            }
            low[ESNK] = PMSK[i];
            up [ESNK] = PMSK[i + 1];
            printf1(" [%d,%d]",low[ESNK],up[ESNK]);
            if (kmax < up[ESNK])
                kmax = up[ESNK];

            ESNK++;
        }
    }
    if (ESNK == 0) {
        printf1(" need a valid definition with sk parameter.\n");
        goto ESKFin;
    }
    newline();
    if (PMLEN < 1)
        PMLEN = 1000;

    printf1("Max record length: %d\n",PMLEN);
    if (PMLEN < kmax) {
        printf1("Inconsistent with sk definition.\n");
        goto ESKFin;
    }
    if (alloc_acc(PMLEN + 2))
        goto ESKFin;   
    if (alloc_acn(PMLEN + 2))
        goto ESKFin;   

    for (k = 0; k < ESNK; ++k) {
        for (i = low[k]; i <= up[k]; ++i)
            AcN[i] = 1;
    }
    rcnt = 0;
    while (fgets(AcC,PMLEN + 1,PMFd)) {
           
        prn_message(++rcnt,0,0);

        p = AcC;
        for (i = 1; i <= PMLEN + 1; ++i) {
            if (AcN[i] == 0 || *p == '\n')
                fprintf(PMF1d,"%c",*p);
            if (!*p || *p == '\n') 
                break;
            p++;
        }
    }
    prn_message(rcnt,1,0);
    printf1("%d records written to: %s\n",rcnt,PMF1dName);
       
    err = 0;
 
ESKFin:
    p_clean();
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

int eselect(void)
{
    register int i;
    int err,low[5],up[5],kmax,klen,rcnt,wcnt,keya;
    register char *p,*q,*s;

    keya = 0;
    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("Selection of records. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 7,2,1))     /* get parameters */
        goto ESELFin;

    printf1("Inputfile: %s\nOutputfile: ",PMFdName);
    if (PMF1Def)
        printf1("%s\n",PMF1dName);
    else {
        printf1("need output file.\n");
        goto ESELFin;
    }
    printf1("Selection keys: ");
    if (PMF2Def)
        printf1("%s\n",PMF2dName);
    else {
        printf1("need a file.\n");
        goto ESELFin;
    }
    printf1("Key definition:");
    kmax = ESNK = 0;
    for (i = 1; i < 10; i += 2) {
        if (PMSK[i] > 0 && PMSK[i + 1] > 0) {
            if (PMSK[i + 1] < PMSK[i]) {
                printf1("error in sk definition.\n");
                goto ESELFin;
            }
            low[ESNK] = PMSK[i] - 1;
            up [ESNK] = PMSK[i + 1];
            printf1(" [%d,%d]",low[ESNK] + 1,up[ESNK]);
            if (kmax < up[ESNK])
                kmax = up[ESNK];

            ESNK++;
        }
    }
    if (ESNK != 2) {
        printf1("\nError: need exactly two keys.\n");
        goto ESELFin;
    }
    newline();
    klen = up[0] - low[0];
    if (klen != up[1] - low[1]) {
        printf1("Error: keys should have identical length.\n");
        goto ESELFin;
    }
    printf1("Max number of keys: %d\n",PMNOC);
    if (PMLEN < 1)
        PMLEN = 1000;

    printf1("Max record length: %d\n",PMLEN);
    if (PMLEN < kmax) {
        printf1("Inconsistent with sk definition.\n");
        goto ESELFin;
    }
    if (alloc_acc(PMLEN + 2))
        goto ESELFin;   
    if (alloc_acd(klen + 1))
        goto ESELFin;   

    if (!(ESKEY  = (char *)calloc(PMNOC * (klen + 1),sizeof(char)))) {
        p_err(-2,1);             
        goto ESELFin;    
    }         
    keya = 1;
    memrq(PMNOC * (klen + 1),sizeof(char));

    /* read keys from if file */

    s = AcC + low[0];
    p = ESKEY;
    ESNK = 0;
    while (fgets(AcC,PMLEN + 1,PMF2d)) {
        if (ESNK >= PMNOC) {
            printf1("Error: exceeded max number of keys.\n");
            goto ESELFin;
        }
        q = s;             
        for (i = 0; i < klen; ++i)
            *p++ = *q++;
        *p++ = '\0';
        prn_message(++ESNK,0,0);
    }
    prn_message(ESNK,1,0);
    printf1("Read %d keys from %s\n",ESNK,PMF2dName);

    if (SILENTFlg < 2)
        printfe("Sorting keys ...\n");

    qsort((char *)ESKEY,ESNK,klen + 1,eselcomp);

    s = AcC + low[1];
    wcnt = rcnt = 0;
    while (fgets(AcC,PMLEN + 1,PMFd)) {
           
        prn_message(++rcnt,0,0);

        /* get key */

        p = AcD;
        q = s;
        for (i = 0; i < klen; ++i)
            *p++ = *q++;
        *p = '\0';

        if (esel_c(klen + 1)) {
            fprintf(PMF1d,"%s",AcC);
            wcnt++;
        }
    }
    prn_message(rcnt,1,0);
    printf1("%d records written to: %s\n",wcnt,PMF1dName);
       
    err = 0;
 
ESELFin:
    if (keya) {
        free((char *)ESKEY);
        memrq(-PMNOC * (klen + 1),sizeof(char));
    }
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  eselcomp()    compare function                                          */

int eselcomp(const void *arg1,const void *arg2)
{     
    register char *p,*q;

    p = (char *)arg1;
    q = (char *)arg2;
    return(strcmp(p,q));
}

/* ------------------------------------------------------------------------ */
/*  esel_c()    return 1 if key occurs in ESKEY[]                           */

int esel_c(int klen)
{     
    int n,low,mid,high;

    low = 1; high = ESNK;

    while (low <= high) {
        mid = (low + high) / 2;
        n = strcmp(AcD,ESKEY + (mid - 1) * klen);
        if (n < 0)
            high = mid - 1;
        else if (n > 0)
            low = mid + 1;
        else {
            if (PMOPT == 2)
                *(ESKEY + (mid - 1) * klen) = '\0';
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

int emerge(void)
{
    register int i,k;
    int err,rcnt,wcnt,n,nl,kmax,nw,klen;
    register char sc,c,*p,*q,*s,*kp;

    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("Merging files. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 6,2,1))     /* get parameters */
        goto EMFin;

    printf1("Inputfile: %s\nOutputfile: ",PMFdName);
    if (PMF1Def)
        printf1("%s\n",PMF1dName);
    else {
        printf1("need output file.\n");
        goto EMFin;
    }
    if (MFN == 0) {
        printf1("Error: need at least one file for merging.\n");
        goto EMFin;
    }  
    sc = ' ';
    if (PMMFlg)
        sc = '0' + (PMM % 10);

    nl = 0;
    for (i = 0; i < MFN; ++i) {
        k = strlen(MFFNAME[i]);
        if (nl < k)
            nl = k;
    }
    klen = MFJ2[0] - MFJ1[0] + 1;
    kmax = 0;
    for (i = 0; i < MFN; ++i) {
        printf1("Merging: %s ",MFFNAME[i]);
        prnchar(' ',nl - strlen(MFFNAME[i]),0);
        printf1("[%2d-%2d] with %s [%2d-%2d]\n",
            MFI1[i],MFI2[i],PMFdName,MFJ1[i],MFJ2[i]);
        if (MFI1[i] < 1 || MFJ1[i] < 1 || klen != (MFI2[i] - MFI1[i] + 1)    
                                       || klen != (MFJ2[i] - MFJ1[i] + 1)) {
            printf1("Error in key definitions.\n");
            goto EMFin;
        }
        if (kmax < MFI2[i])
            kmax = MFI2[i];
        if (kmax < MFJ2[i])
            kmax = MFJ2[i];
    }
    if (PMLEN < 1)
        PMLEN = 1000;

    printf1("Max record length: %d\n",PMLEN);
    if (PMLEN < kmax) {
        printf1("Inconsistent with key definitions.\n");
        goto EMFin;
    }
    if (alloc_acc((PMLEN + 2) * (MFN + 1)))
        goto EMFin;   

    if (alloc_acd((klen + 1) * (MFN + 1)))
        goto EMFin;   

    if (alloc_acn(MFN))
        goto EMFin;   

    if (alloc_acm(MFN))
        goto EMFin;   

    if (alloc_aci(MFN))
        goto EMFin;   

    if (alloc_acptr(MFN))
        goto EMFin;   

    for (i = 0; i < MFN; ++i) {
        p = AcPtr[i] = AcC + (i + 1) * (PMLEN + 2);
        q = p + MFI1[i] - 1;
        if (!fgets(p,PMLEN + 1,MFFD[i])) {
            printf1("Error: can't read from %s.\n",MFFNAME[i]);
            goto EMFin;
        }
        n = 0;
        while (*p && *p != LF && *p != CR) {
            p++;
            n++;
        }
        *p = '\0';
        AcN[i] = n;               
        AcM[i] += 1;

        kp = AcD + (klen + 1) * (i + 1);
        for (k = 0; k < klen; ++k)
            *kp++ = *q++;
    }

    nw = wcnt = rcnt = 0;
    while (fgets(AcC,PMLEN + 1,PMFd)) {
           
        prn_message(++rcnt,0,0);

        p = AcC;
        n = 0;
        while (*p && *p != LF && *p != CR) {
            p++;
            n++;
        }
        *p = '\0';
        if (n < MFJ2[0]) {
            nw++;
            continue;
        }
        fprintf(PMF1d,"%s",AcC);
        if (XSEPC)
            fprintf(PMF1d,"%c",XSEPC);

        /* get key from input file */

        *(AcC + MFJ2[0]) = '\0';
        kp = AcD;
        q = AcC + MFJ1[0] - 1;

        if (rcnt > 1) {             /* check ordering */
            if (strcmp(kp,q) > 0) {
                printf1("\nError: %s not sorted in ascending order.\n",PMFdName);
                printf1("Check records %d and %d.\n",rcnt-1,rcnt);
                goto EMFin;
            }
        }
        for (k = 0; k < klen; ++k)
            *kp++ = *q++;
   
        for (i = 0; i < MFN; ++i) {
            if (AcN[i] > 0) {
                s = AcPtr[i];
                p = s + MFI1[i] - 1;
                q = p + klen;
                c = *q;
                *q = '\0';

                while (strcmp(AcD,p) > 0) {

                    s = AcPtr[i];
                    if (!fgets(s,PMLEN + 1,MFFD[i])) {
                        AcN[i] = -AcN[i];
                        break;
                    }
                    n = 0;
                    while (*s && *s != LF && *s != CR) {
                        s++;
                        n++;
                    }
                    *s = '\0';

                    AcM[i] += 1;
                    c = *q;
                    *q = '\0';

                    kp = AcD + (klen + 1) * (i + 1);
                    if (strcmp(kp,p) > 0) {
                        printf1("\nError: %s not sorted in ascending order.\n",MFFNAME[i]);
                        printf1("Check records %d and %d.\n",AcM[i]-1,AcM[i]);
                        goto EMFin;
                    }
                    for (k = 0; k < klen; ++k)
                        *kp++ = *p++;
                    p = AcPtr[i] + MFI1[i] - 1;
                }
                n = strcmp(AcD,p);          
                *q = c;
                if (n == 0) {
                    fprintf(PMF1d,"%s",AcPtr[i]);
                    AcI[i] += 1;
                }
                else
                    fprnchar(PMF1d,sc,iabs(AcN[i]),0);
            }
            else  
                fprnchar(PMF1d,sc,-AcN[i],0);

            if (XSEPC)
                fprintf(PMF1d,"%c",XSEPC);
        }
        fprintf(PMF1d,"\n");
        wcnt++;
    }
    prn_message(rcnt,1,0);

    printf1("\nRead %8d records from: %s\n",rcnt,PMFdName);
    for (i = 0; i < MFN; ++i) {
        printf1("Read %8d records from: %s ",AcM[i],MFFNAME[i]);
        prnchar(' ',nl - strlen(MFFNAME[i]),0);
        printf1(" matched: %d\n",AcI[i]);
    }

    printf1("\n%d records written to: %s\n",wcnt,PMF1dName);

    if (nw > 0) {
        printf1("\nWarning: skipped %d records from %s\n",nw,PMFdName);
        printf1("         that have record length less than %d.\n",MFJ2[0]);
    }
    err = 0;
 
EMFin:
    p_clean();
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

int ejoin(void)
{
    int err,r,nv1,nv2,nx1,nx2,id1,id2,r1,r2;
    int wrec,idn1,idn2,njoin,l1,l2,rcnt;

    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("Join two event data files. Current memory: %d bytes.\n\n",MemReq);

    if (parm(CmdBuf + 5,1,1))     /* get parameters */
        goto EJFin;

    printf1("Outputfile : %s\n",PMFdName);

    if (PMMax < 1)
        PMMax = 1000;

    if (PMLEN < 1)
        PMLEN = 1000;

    printf1("Max block size: %d\n",PMMax);
    printf1("Max number of levels: %d\n",PMNW);
    printf1("Max record length: %d\n",PMLEN);

    if (alloc_acs(2 * PMMax + 1))
        goto EJFin;   
    if (alloc_aci(2 * PMMax * PMNW + 1))
        goto EJFin;   
    if (alloc_acj(2 * PMMax * PMNW + 1))
        goto EJFin;   

    njoin = idn1 = idn2 = nv1 = nv2 = nx1 = nx2 = 0;
    wrec = EJRECN1 = EJRECN2 = 0;
    EJEOF1 = EJEOF2 = 1;

    if (PMIF1Def) {
        if (alloc_acc(PMLEN + 2))
            goto EJFin;   
        nv1 = ejoin_nv(PMIF1d,AcC,&EJRECN1);
        printf1("Inputfile 1: %s  (%d variables)\n",PMIF1Name,nv1);
        if (nv1 < 6) {
            printf1("Error: need at least 6 variables.\n");
            goto EJFin;
        }
        nx1 = nv1 - 6;
        EJEOF1 = 0;
        if (alloc_acn(6 * PMMax))
            goto EJFin;   
        if (alloc_acx(nx1 * PMMax + 1))
            goto EJFin;   
    }
    if (PMIF2Def) {
        if (alloc_acd(PMLEN + 2))
            goto EJFin;   
        nv2 = ejoin_nv(PMIF2d,AcD,&EJRECN2);
        printf1("Inputfile 2: %s  (%d variables)\n",PMIF2Name,nv2);
        if (nv2 < 6) {
            printf1("Error: need at least 6 variables.\n");
            goto EJFin;
        }
        nx2 = nv2 - 6;
        EJEOF2 = 0;
        if (alloc_acm(6 * PMMax))
            goto EJFin;   
        if (alloc_acy(nx2 * PMMax + 1))
            goto EJFin;   
    }

    while (EJEOF1 == 0 || EJEOF2 == 0) {

        if (EJEOF1 == 0) {
            r1 = ejoin_rd(PMIF1d,PMIF1Name,AcC,&EJEOF1,&EJRECN1,AcN,AcX,nx1,&l1);
            if (r1 < 0)
                goto EJFin;
            id1 = AcN[0];
            idn1++;
        }
        else {
            id1 = INTMAX;
            r1 = 0;
            l1 = 1;
        }
        if (EJEOF2 == 0) {
            r2 = ejoin_rd(PMIF2d,PMIF2Name,AcD,&EJEOF2,&EJRECN2,AcM,AcY,nx2,&l2);
            if (r2 < 0)
                goto EJFin;
            id2 = AcM[0];
            idn2++;
        }
        else {
            id2 = INTMAX;
            r2 = 0;
            l2 = 1;
        }

EJCONT:
        while (id1 < id2) {
            if ((r = ejoin_proc(id1,r1,0,l1,1,nx1,nx2)) < 0)
                goto EJFin;
            wrec += r;

            rcnt = imax(EJRECN1,EJRECN2);
            prn_message(rcnt,0,0);
            if (PMNOCFlg && rcnt >= PMNOC)
                goto EJFin1;

            if (EJEOF1 == 0) {
                r1 = ejoin_rd(PMIF1d,PMIF1Name,AcC,&EJEOF1,&EJRECN1,AcN,AcX,nx1,&l1);
                if (r1 < 0)
                    goto EJFin;
                id1 = AcN[0];
                idn1++;
            }
            else {
                id1 = INTMAX;
                r1 = 0;
                l1 = 1;
            }
        }

        while (id2 < id1) {
            if ((r = ejoin_proc(id2,0,r2,1,l2,nx1,nx2)) < 0)
                goto EJFin;
            wrec += r;

            rcnt = imax(EJRECN1,EJRECN2);
            prn_message(rcnt,0,0);
            if (PMNOCFlg && rcnt >= PMNOC)
                goto EJFin1;

            if (EJEOF2 == 0) {
                r2 = ejoin_rd(PMIF2d,PMIF2Name,AcD,&EJEOF2,&EJRECN2,AcM,AcY,nx2,&l2);
                if (r2 < 0)
                    goto EJFin;
                id2 = AcM[0];
                idn2++;
            }
            else {
                id2 = INTMAX;
                r2 = 0;
                l2 = 1;
            }
        }
        if (id1 < id2)
            goto EJCONT;

        if (id1 == id2 && id1 < INTMAX) {

            if ((r = ejoin_proc(id1,r1,r2,l1,l2,nx1,nx2)) < 0)
                goto EJFin;
            wrec += r;
            njoin++;
        }
        rcnt = imax(EJRECN1,EJRECN2);
        prn_message(rcnt,0,0);
        if (PMNOCFlg && rcnt >= PMNOC)
            break;
    }

EJFin1:
    newline();
    if (PMIF1Def)
        printf1("Read %d records (%d blocks) from: %s\n",EJRECN1,idn1,PMIF1Name);
    if (PMIF2Def)
        printf1("Read %d records (%d blocks) from: %s\n",EJRECN2,idn2,PMIF2Name);
    if (PMIF1Def && PMIF2Def)
        printf1("Number of common blocks: %d\n",njoin);
    printf1("%d records written to %s\n",wrec,PMFdName);

    err = 0;
 
EJFin:
    p_clean();
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  ejoin_nv()      return number of variables found in first record        */
/*                  of file fn.                                             */

int ejoin_nv(FILE *fd,char *buf,int *nrec)
{
    int nv;
    double x;
    register char *p;

    nv = 0;
    while (fgets(buf,PMLEN + 1,fd)) {
        *nrec += 1;
        if (check_drec(buf)) {      /* check for data records */
            p = buf;
            while (*p) {
                p = skip_b(p);
                if (sscanf(p,"%lf",&x) != 1)
                    break;
                nv++;
                p = skip_dbl(p);
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

int ejoin_rd(FILE *fd,char *fname,char *buf,int *eof,int *nrec,
    int *bvar,double *xvar,int nx,int *mxlev)
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

        if (n >= PMMax) {
            printf1("Error: exceeded max block size in file %s\n",fname);
            return(-1);
        }
        p = buf;
        nv = 0;

        while (*p) {
            p = skip_b(p);
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
            p = skip_dbl(p);
        }
        if (nv < nx + 6) {
            printf1("Error: found less than %d variables in file %s.\n",nx + 6,fname);
            return(-2);
        }
        if (n == 0)
            id = bvar[0];
        n++;

        *eof = 1;

        while (fgets(buf,PMLEN + 1,fd)) {
            *nrec += 1;
            if (!check_drec(buf))        /* check for data records */
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
            *mxlev = imax(*mxlev,l);
            if (l == 0)  
                break;
        }
        j++;
    }
    *mxlev += 1;

    if (*mxlev > PMNW) {
        printf1("Error: exceeded max number of levels in file %s.\n",fname);
        return(-3);
    }
    return(n);
}

/* -##--------------------------------------------------------------------- */
/*  ejoin_proc()   Process one or two blocks. Return number of records      */
/*                 written to output file, or -1 if error.                  */

int ejoin_proc(int id,int r1,int r2,int l1,int l2,int nx1,int nx2)
{
    register int i,j,k,l;
    int n,nn,nnn,nfmt,t,ii,iii,ll,jj,wrec,ts,tf;

    n = wrec = 0;

    for (i = 0; i < r1; ++i) {
        AcS[n++] = AcN[i * 6 + 3];
        AcS[n++] = AcN[i * 6 + 4];
    }
    for (i = 0; i < r2; ++i) {
        AcS[n++] = AcM[i * 6 + 3];
        AcS[n++] = AcM[i * 6 + 4];
    }
    if (sorti(n,AcS,0))
        return(-1);  

    nn = 0;
    for (i = 1; i < n; ++i) {
        if (AcS[nn] < AcS[i]) 
            AcS[++nn] = AcS[i];
    }
    nn++;

    i = nn * imax(l1,l2);  
    for (k = 0; k < i; ++k)
        AcI[k] = AcJ[k] = EJMVAL;

    k = 0;
    t = AcS[0];
    for (i = 0; i < r1; ++i) {

        ts = AcN[i * 6 + 3];
        tf = AcN[i * 6 + 4];

        while (t < ts)  
            t = AcS[++k];

        while (t < tf) {
            l = 0;
            AcI[k] = i;

            ii = i;
            while (++ii < r1) {
                if (AcN[ii * 6 + 3] <= t) {
                    if (t < AcN[ii * 6 + 4]) {
                        l++;
                        AcI[l * nn + k] = ii;
                    }
                }   
                else
                    break;
            }
            t = AcS[++k];
        }
    }
    k = 0;
    t = AcS[0];
    for (i = 0; i < r2; ++i) {

        ts = AcM[i * 6 + 3];
        tf = AcM[i * 6 + 4];

        while (t < ts)  
            t = AcS[++k];

        while (t < tf) {
            l = 0;
            AcJ[k] = i;
            ii = i;
            while (++ii < r2) {
                if (AcM[ii * 6 + 3] <= t) {
                    if (t < AcM[ii * 6 + 4]) {
                        l++;
                        AcJ[l * nn + k] = ii;
                    }
                }   
                else
                    break;
            }
            t = AcS[++k];
        }
    }
    nnn = 0;
    for (k = 0; k < nn - 1; ++k) {
        ii = 0;
        for (l = 0; l < l1; ++l) {
            if (AcI[l * nn + k] < 0)
                break;
            ii = l;
        }
        iii = 0;
        for (l = 0; l < l2; ++l) {
            if (AcJ[l * nn + k] < 0)
                break;
            iii = l;
        }
        nnn += imax(++ii,++iii);
    }
    for (k = 0; k < nn - 1; ++k) {
        ii = 0;
        for (l = 0; l < l1; ++l) {
            if (AcI[l * nn + k] < 0)
                break;
            ii = l;
        }
        iii = 0;
        for (l = 0; l < l2; ++l) {
            if (AcJ[l * nn + k] < 0)
                break;
            iii = l;
        }
        i = imax(++ii,++iii);
  
        for (j = 0; j < i; ++j) {
            ejoin_prn(0,id);
            ejoin_prn(1,nnn);
            ejoin_prn(2,k + 1);
            ejoin_prn(3,j);
            ejoin_prn(4,AcS[k]);
            ejoin_prn(5,AcS[k + 1]);
            nfmt = 6;
            if (PMIF1Def) {
                l = AcI[j * nn + k];
                if (l >= 0)
                    ejoin_prn(nfmt++,AcN[l * 6 + 5]);
                else
                    ejoin_prn(nfmt++,EJMVAL);
            }
            if (PMIF2Def) {
                l = AcJ[j * nn + k];
                if (l >= 0)
                    ejoin_prn(nfmt++,AcM[l * 6 + 5]);
                else
                    ejoin_prn(nfmt++,EJMVAL);
            }
            if (PMIF1Def && nx1 > 0) {
                l = AcI[j * nn + k];
                if (l >= 0) {
                    ll = l * nx1;
                    for (jj = 0; jj < nx1; ++jj)  
                        ejoin_prn1(jj,AcX[ll++]);
                }
                else {
                    for (jj = 0; jj < nx1; ++jj)  
                        ejoin_prn1(jj,(double)EJMVAL);
                }
            }
            if (PMIF2Def && nx2 > 0) {
                l = AcJ[j * nn + k];
                if (l >= 0) {
                    ll = l * nx2;
                    for (jj = 0; jj < nx2; ++jj)  
                        ejoin_prn2(jj,AcY[ll++]);
                }
                else {
                    for (jj = 0; jj < nx2; ++jj)  
                        ejoin_prn2(jj,(double)EJMVAL);
                }
            }
            fprintf(PMFd,"\n");
            wrec++;
        }
    }
    return(wrec);
}

/* ------------------------------------------------------------------------ */
/*  ejoin_prn()   print a number                                            */

void ejoin_prn(int nfmt,int val)
{
    int n;

    n = PMXFmtN[0];    
    if (n > 0) {            
        n = imin(nfmt,n - 1);
        fprintf(PMFd,PMXFmtS[0] + n * PMXFmtLen,(double)val);
    }
    else
        fprintf(PMFd,"%d ",val);
}

/* ------------------------------------------------------------------------ */
/*  ejoin_prn1()   print a number                                           */

void ejoin_prn1(int nfmt,double val)
{
    int n;

    n = PMXFmtN[1];    
    if (n > 0) {            
        n = imin(nfmt,n - 1);
        fprintf(PMFd,PMXFmtS[1] + n * PMXFmtLen,val);
    }
    else
        fprintf(PMFd,"%g ",val);
}

/* ------------------------------------------------------------------------ */
/*  ejoin_prn2()   print a number                                           */

void ejoin_prn2(int nfmt,double val)
{
    int n;

    n = PMXFmtN[2];    
    if (n > 0) {            
        n = imin(nfmt,n - 1);
        fprintf(PMFd,PMXFmtS[2] + n * PMXFmtLen,val);
    }
    else
        fprintf(PMFd,"%g ",val);
}



