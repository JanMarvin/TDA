/****************************************************************************/
/*  t_rzoo                                                                  */
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
#include "t_zoo.h"
#include "t_var.h"
#include "t_gdat.h"
#include "t_alloc.h"
#include "t_gf.h"

/*  functions in t_rzoo.c */

void arcd_off(void);
int arcd(void);
int arcd_alloc(int opt);
int get_len(char *p);
char *cpy_nam(char *s, char *p, int n);
void prn_afiles(void);
int get_drec(int fn, char *buf, int nmax);
int get_record(int fn, char *buf, int nmax);
void get_rdat(int fn);
int alloc_avar(int idx,int n,int opt);
int check_avar(int n,int adic);
int get_avar(void);
void get_astr(char *buf,int i);
int arcc(void);
int arcv(void);
int arcvc(void);
int get_afmt(int *w1,int *w2);
int iv_comp(const void *arg1,const void *arg2);

/* ------------------------------------------------------------------------ */
/*  Global variables                                                        */

int ARCDef = 0;             /* set if archive successfully opened           */
char ZADNam[FNMaxLen];      /* name of archive description file             */
char ZOONam[FNMaxLen];      /* name of zoo data archive                     */
FILE *ZOOFd;                /* file handle for zoo archive                  */
int ZOOFdo = 0;             /* set if opened                                */

int ZANF = 0;               /* number of files in description file          */

char **ZAFNam;              /* names of files                               */
int ZAFNamA = 0;        
int ZAFNamAI = 0;       
int *ZAFNRec;               /* number of records                            */
int ZAFNRecA = 0; 
short *ZAFNum;              /* logical file number                          */
int ZAFNumA = 0; 
short *ZAFTyp;              /* type of file                                 */
int ZAFTypA = 0; 
short *ZAFRLen;             /* record length                                */
int ZAFRLenA = 0; 
short *ZAFNVar;             /* number of variables                          */
int ZAFNVarA = 0; 
int ZABLen = 0;             /* max record length of files                   */
char *ZABuf;                /* Record buffer  for files                     */
int FDefLen = 0;            /* max length of file names                     */
char *ZAFReq;               /* flags for required files                     */
int ZAFReqA = 0; 
int VFN = -1;               /* Number of variable description file          */
/* ------------------------------------------------------------------------ */
int *ZAZOfs;                /* Direntry offset in archive                   */
int ZAZOfsA = 0; 
int *ZAZSiz;                /* Original size of file in bytes               */
int ZAZSizA = 0; 
char *ZAZTyp;               /* Packing method (1 - 2)                       */
int ZAZTypA = 0; 
char ZAInit;                /* Set if decompression initialized             */
int ZAAlloc;                /* Set if memory allocated                      */
char ZAEOF;                 /* Set if EOF reached                           */
long FPtr;                  /* File pointer for the files                   */

/* ------------------------------------------------------------------------ */
/*  For each file we allocate a read buffer where the decompression is      */
/*  done, and a line buffer.                                                */

char *Out_Buf_Adr;          /* Output buffer used for decompression         */
char *In_Buf_Adr;           /* Input buffer used for decompression          */
char *Out_Buf_Ptr;          /* Pointer to output buffer                     */
int Out_Buf_Cnt;            /* Character count of output buffer             */

/* ------------------------------------------------------------------------ */
/*  arcd_off    Turn off arrchive.                                          */

void arcd_off(void)
{
    if (ARCDef) {
        dbf_init(0);
        arcd_alloc(0);
        ARCDef = ZANF = 0;
        if (ZOOFdo) {
            fclose(ZOOFd);
            ZOOFdo = 0;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  arcd()      Command in CmdBuf:                                          */
/*              a) arcd = name_of_archive_description_file.                 */
/*                 then check this file and open access to archive; if      */
/*                 successful, set ARCDef = 1.                              */
/*              b) arcd = off                                               */
/*                 then turn off access to archive and free all previously  */
/*                 allocated memory and set ARCDef = 0.                     */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int arcd(void)
{
    FILE *fd;
    register int i,j;
    register char *p,*q;
    int err,m,n,first,rec,n1,n2,n3,n4;

    err = -1;
    if (check_cmd(1))
        return(-1);

    if (!strcmp(CmdBuf,"arcd")) {
        if (ARCDef == 0)
            printf1("Data archive not defined.\n");
        else {
            printf1("Archive: %s\n",ZOONam);
            printf1("Description file: %s\n",ZADNam);
            prn_afiles();
        }
        return(0);
    }

    if (ARCDef) {
        dbf_init(0);
        arcd_alloc(0);
        ARCDef = ZANF = 0;
        if (!strcmp(CmdBuf,"arcd=off")) {
            printf1("Access to data archive %s turned off.\nCurrent memory: %d bytes.\n",
                                                             ZOONam,MemReq);
            return(0);
        }
    }
    else if (!strcmp(CmdBuf,"arcd=off")) {
        printf1("No data archive defined (command ignored).\n");
        return(0);
    }
    strncpy(ZADNam,CmdBuf + 5,FNMaxLen);

    if (!(fd = fopen(ZADNam,OPEN_RD))) {
        printf1("Can't open archive description file: %s\n",ZADNam);
        return(-1);
    }
    printf1("Reading archive description file: %s\n",ZADNam);
        
    if (alloc_acc(RLMaxDef))
        goto AFin;

    /*  count number of files in description file in ZANF */

    ZANF = 0;
    while (fgets(AcC,RLMaxDef,fd)) {
        if ((p = check_comment(AcC)) != NULL)  
            ZANF++;
    }
    if (--ZANF < 1) {
        printf1("Error: no definition of files in: %s\n",ZADNam);
        goto AFin;
    }
    if (arcd_alloc(1)) {
        printf1("Insufficient memory for archive data structures.\n");
        goto AFin;
    }
    fseek(fd,0L,0);

    ZABLen = RLMaxDef;  /* max buffer size for data files, updated below */
    rec = m = 0;
    first = 1;
    err = -2;

    while (m < ZANF && fgets(AcC,RLMaxDef,fd)) {

        rec++;
        if ((p = check_comment(AcC)) != NULL) {

            q = p + strlen(p);
            while (--q > p && (*q == '\n' || *q == LF || *q == CR)) ;
            *++q = '\0';
   
            if (first) {              
                n = get_len(p);
                if (n >= FNMaxLen)  
                    goto AFin;
                cpy_nam(ZOONam,p,n);
                first = 0;
            }
            else {       
                if (!*p || sscanf(p,"%d",&n1) != 1 || n1 < 0)  
                    goto AFin;
                ZAFNum[m] = (short)n1;
                p = skip_cb(p);
                n = get_len(p);
                if (n >= FNMaxLen)  
                    goto AFin;
                if (FDefLen < n)
                    FDefLen = n;

                p = cpy_nam(ZAFNam[m],p,n);
                p = skip_b(p);
                n1 = n2 = n3 = n4 = -1;
                if (!*p || sscanf(p,"%d %d %d %d",&n1,&n2,&n3,&n4) != 4 ||
                    (n1 != 1 && n1 != 2) || n2 < 0 || n3 < 1 || n4 < 0) {
                    goto AFin;
                }
                ZAFTyp[m]  = (short)n1;
                ZAFRLen[m] = (short)n2;
                ZAFNRec[m] = n3;
                ZAFNVar[m] = (short)n4;    
                     
                if (ZABLen < ZAFRLen[m])
                    ZABLen = ZAFRLen[m];
                m++;
            }
        }
    }
    if (m != ZANF) {
        err = -2;
        goto AFin;
    }
    rec = 0;
    if (first) {
        printf1("Error: can't find name of ZOO archive.\n");
        goto AFin;
    }
    printf1("ZOO data archive: %s\n",ZOONam);

    /*  Open the ZOO archive and get header and directory information.      */

    if (!(ZOOFd = fopen(ZOONam,OPEN_RB))) {
        printf1("Can't open this file.\n");
        goto AFin;
    }
    ZOOFdo = 1;

    if (get_zoo())  
        goto AFin;

    printf1("Checking definition of files in archive.\n");
    err = m = n = 0;

    for (i = 0; i < ZANF; ++i) {
            
        if (!ZAZOfs[i]) {
            printf1("Can't find file %s in the archive.\n",ZAFNam[i]);
            err = -1;
        }
        else if (ZAZTyp[i] < 1 || ZAZTyp[i] > 2) {
            printf1("File %s has unknown packing method.\n",ZAFNam[i]);
            err = -1;
        }
        for (j = 0; j < i; ++j) {
            if (ZAFNum[i] == ZAFNum[j])  
                n++;
        }
        if (ZAFTyp[i] == 2) {
            VFN = i;
            m++;
        }
    }
    if (n > 0) {
        printf1("Logical file ID numbers are not unique.\n");
        err = -1;
    }
    if (m == 0) {
        printf1("Can't find a variable description file (type 2).\n");
        err = -1;
    }
    else if (m > 1) {
        printf1("Found more than one variable description file (type 2).\n");
        err = -1;
    }
    prn_afiles();           /* print list of files */

    if (err)
        goto AFin;

    ARCDef = 1;     /* set global flag for successfully opened archive */
    err = 0;      

AFin:
    fclose(fd);
    alloc_acc(0);

    if (err) {
        if (err == -2 && rec > 0)  
            printf1("Error in record %d of archive description file.\n",rec);
        arcd_alloc(0);
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  arcd_alloc(opt)     If opt != 0 allocate data structures for ZANF       */
/*                      data files, otherwise free the memory.              */
/*                      Return 0 if OK, -1 if error.                        */

int arcd_alloc(int opt)
{
    register int i;
    int err;

    if (opt == 0) {
        err = 0;
        goto AAFin;
    }   
    err = -1;

    if (!(ZAFNam = (char **) calloc(ZANF,sizeof(char *))))  
        goto AAFin; 
    ZAFNamA = ZANF;
    memrq(ZAFNamA,sizeof(char *));
   
    if (!(ZAFNum = (short *) calloc(ZANF,sizeof(short))))  
        goto AAFin; 
    ZAFNumA = ZANF;
    memrq(ZAFNumA,sizeof(short));

    if (!(ZAFTyp = (short *) calloc(ZANF,sizeof(short))))  
        goto AAFin; 
    ZAFTypA = ZANF;
    memrq(ZAFTypA,sizeof(short));

    if (!(ZAFRLen = (short *) calloc(ZANF,sizeof(short))))  
        goto AAFin; 
    ZAFRLenA = ZANF;
    memrq(ZAFRLenA,sizeof(short));

    if (!(ZAFNVar = (short *) calloc(ZANF,sizeof(short))))  
        goto AAFin; 
    ZAFNVarA = ZANF;
    memrq(ZAFNVarA,sizeof(short));

    if (!(ZAFNRec = (int *) calloc(ZANF,sizeof(int))))  
        goto AAFin; 
    ZAFNRecA = ZANF;
    memrq(ZAFNRecA,sizeof(int));

    if (!(ZAZSiz = (int *) calloc(ZANF,sizeof(int))))  
        goto AAFin; 
    ZAZSizA = ZANF;
    memrq(ZAZSizA,sizeof(int));

    if (!(ZAZTyp = (char *) calloc(ZANF,sizeof(char))))  
        goto AAFin; 
    ZAZTypA = ZANF;
    memrq(ZAZTypA,sizeof(char));

    if (!(ZAFReq = (char *) calloc(ZANF,sizeof(char))))  
        goto AAFin; 
    ZAFReqA = ZANF;
    memrq(ZAFReqA,sizeof(char));

    if (!(ZAZOfs = (int *) calloc(ZANF,sizeof(int))))  
        goto AAFin; 
    ZAZOfsA = ZANF;
    memrq(ZAZOfsA,sizeof(int));

    for (i = 0; i < ZANF; ++i) {
        if (!(ZAFNam[i] = (char *) calloc(FNMaxLen,sizeof(char))))  
            goto AAFin;
        memrq(FNMaxLen,1);
        ZAFNamAI++;
    }
    return(0);

AAFin:
    if (ZAFNamAI > 0) {
        for (i = 0; i < ZAFNamAI; ++i) {
            free(ZAFNam[i]);
            memrq(-FNMaxLen,1);
        }
        ZAFNamAI = 0;
    }
    if (ZAZOfsA > 0) {
        free((char *)ZAZOfs);
        memrq(-ZAZOfsA,sizeof(int));
        ZAZOfsA = 0;
    }      
    if (ZAFReqA > 0) {
        free((char *)ZAFReq);
        memrq(-ZAFReqA,sizeof(char));
        ZAFReqA = 0;
    }   
    if (ZAZTypA > 0) {
        free((char *)ZAZTyp);
        memrq(-ZAZTypA,sizeof(char));
        ZAZTypA = 0;
    }   
    if (ZAZSizA > 0) {
        free((char *)ZAZSiz);
        memrq(-ZAZSizA,sizeof(int));
        ZAZSizA = 0;
    }   
    if (ZAFNRecA > 0) {
        free((char *)ZAFNRec);
        memrq(-ZAFNRecA,sizeof(int));
        ZAFNRecA = 0;
    }
    if (ZAFNVarA > 0) {
        free((char *)ZAFNVar);
        memrq(-ZAFNVarA,sizeof(short));
        ZAFNVarA = 0;
    }
    if (ZAFRLenA > 0) {
        free((char *)ZAFRLen);
        memrq(-ZAFRLenA,sizeof(short));
        ZAFRLenA = 0;
    }
    if (ZAFTypA > 0) {
        free((char *)ZAFTyp);
        memrq(-ZAFTypA,sizeof(short));
        ZAFTypA = 0;
    }
    if (ZAFNumA > 0) {
        free((char *)ZAFNum);
        memrq(-ZAFNumA,sizeof(short));
        ZAFNumA = 0;
    }
    if (ZAFNamA > 0) {
        free((char *)ZAFNam);
        memrq(-ZAFNamA,sizeof(char *));
        ZAFNamA = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_len(p)    Get length of string at p until next blank or tab char.   */

int get_len(char *p)
{
    register int n = 0;

    while (*p && *p != ' ' && *p != '\t' && *p != LF && *p != CR && *p != '\n') {
        p++;
        n++;
    }
    return(n);
}

/* ------------------------------------------------------------------------ */
/*  cpy_nam(s,p,n)   Copy n characters from p to s. Append \0               */
/*                   Return pointer to next character.                      */

char *cpy_nam(char *s, char *p, int n)
{
    register int i;
    for (i = 0; i < n; ++i)
        *s++ = *p++;
    *s = '\0';
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  prn_afiles()  Print list of files described in archive description file */

void prn_afiles(void)
{
    register int i;

    printf1("\nFN  Type  RLen  Records  NVar    Size    M  Name\n");
    prnchar('-',48,1); 
                        
    for (i = 0; i < ZANF; ++i) {

        printf1("%3d %3d %6d %8d %5d ",
                    ZAFNum[i],ZAFTyp[i],ZAFRLen[i],ZAFNRec[i],ZAFNVar[i]);
        printf1("%9d %2d  %s\n",ZAZSiz[i],(int)ZAZTyp[i],ZAFNam[i]);
    }
    newline();
}

/* ------------------------------------------------------------------------ */
/*  get_drec(fn,buf,nmax)                                                   */
/*      Reads a record from file number fn into the buffer buf. Return the  */
/*      number of bytes in the buffer, or zero if EOF, or -1 if error.      */
/*      nmax is the maximum buffer length.                                  */
/*                                                                          */
/*      EOL characters are NOT counted and not put into the buffer; empty   */
/*      lines are skipped.                                                  */

int get_drec(int fn, char *buf, int nmax)
{
    register int n;
    register char *p;

    while ((n = get_record(fn,buf,nmax)) > 0) {
        p = buf + n - 1;
        while (n > 0 && (*p == LF || *p == CR)) {
            n--;
            *p-- = '\0';
        }
        if (n > 0)
            break;
    }
    return(n);
}

/* ------------------------------------------------------------------------ */
/*  get_record(fn,buf,nmax)                                                 */
/*      Reads a record from file number fn into the buffer buf. Return the  */
/*      number of bytes in the buffer, or zero if EOF, or -1 if error.      */
/*      nmax is the maximum buffer length.                                  */
/*                                                                          */
/*      The length of the record is determined according to ZAFRLen[fn].    */
/*      If ZAFRLen[fn] > 0 this length is taken. Otherwise the record is    */
/*      up to the next occurence of EOL. EOL is substituted by a single     */
/*      '\n' character. Always an \0 is appended!                           */

int get_record(int fn, char *buf, int nmax)
{
    register int n,i,j = 0;
    register char c,*q;

    q = buf;
    if (ZAFRLen[fn] <= 0) {
        while (j < nmax) {
            while (Out_Buf_Cnt > 0) {
                c = *Out_Buf_Ptr++;
                Out_Buf_Cnt--;
                if (c == CR || c == LF) {
                    if (Out_Buf_Cnt <= 0) 
                        get_rdat(fn);
                    if (Out_Buf_Cnt > 0) {
                        if (c == CR && *Out_Buf_Ptr == LF) {
                            Out_Buf_Ptr++;
                            Out_Buf_Cnt--;
                        }   
                    }
                    *q++ = '\n';
                    *q = '\0';
                    return(++j);
                }
                else {
                    *q++ = c;
                    if (++j >= nmax) {
                        printf1("\nError in retrieving file: %s\n",ZAFNam[fn]); 
                        printf1("Can't find an end-of-record character.\n");
                        printf1("Exceeded max record length: %d characters.\n",nmax);
                        return(-1);
                    }
                }
            }
            get_rdat(fn);
            if (Out_Buf_Cnt <= 0) {
                if (j > 0) {
                    printf1("\nError in retrieving file: %s\n",ZAFNam[fn]); 
                    printf1("Can't find an end-of-record character.\n");
                    printf1("Remaining characters: %d\nFound:",j);
                    q = buf;
                    for (i = 0; i < j; ++i)  
                        printf1(" %02x",*q++);
                    printf1("\n\n");
                    if (j > 1 || *buf != 0x1a)
                        return(-1);     
                }
                *q = '\0';
                return(0);
            }
        }
    }
    else {
        if ((n = ZAFRLen[fn]) > nmax)
            n = nmax;
        while (j < n) {
            if (!Out_Buf_Cnt) {
                get_rdat(fn);
                if (!Out_Buf_Cnt)  
                    break;
            } 
            if ((i = Out_Buf_Cnt) > n - j)
                i = n - j; 
            memcpy(q,Out_Buf_Ptr,i);
            q += i;
            Out_Buf_Ptr += i;
            Out_Buf_Cnt -= i;
            j += i;
        }
        *(buf + j) = '\0';
        return(j);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  get_rdat(fn)    Get new piece of uncompressed data for file fn          */
/*                  in the buffer Out_Buf_Adr. The number of bytes is       */
/*                  recorded in Out_Buf_Cnt. This function is simply        */
/*                  a switch in the functions lzd, and lzh_decode,          */
/*                  depending on the compression mode of the file.          */

void get_rdat(int fn)
{
    if (!ZAEOF) {
        switch (ZAZTyp[fn]) {
            case 1:   lzd(fn);
                      break; 
            case 2:   lzh_decode(fn);
                      break; 
            default:  printfe("\nUnknown compression type %d\n",(int)ZAZTyp[fn]);
                      gerr_exit(32);

        }
        Out_Buf_Ptr = Out_Buf_Adr;
    }
}

/* ------------------------------------------------------------------------ */
/*  alloc_avar(idx,n,opt)                                                   */
/*                                                                          */
/*                      If opt != 0, allocate the following arrays for n    */
/*                      archive variables, beginning with index idx.        */
/*                                                                          */
/*                      AVIdx       internal variable number                */
/*                      AVOff       offset                                  */
/*                      AVLen       length of variable                      */
/*                      AVFmt1      format                                  */
/*                      AVFmt2                                              */
/*                      AVMBlnk     missing values: blank                   */
/*                      AVMStar     missing values: star                    */
/*                      AVMPnt      missing values: point                   */
/*                      AVMGen      missing values: general                 */
/*                                                                          */
/*                      Set AVIdx to internal variable number and AVFmt1 to */
/*                      -length of archive variable name.                   */
/*                                                                          */
/*                      If opt == 0 free memory.                            */
/*                      Return 0 if OK, -1 if error.                        */

int alloc_avar(int idx,int n,int opt)
{
    register int i,j,l;
    int err = -1;

    if (opt == 0) {
        err = 0;
        goto AAVFin;
    }
    if (n <= 0)
        return(0);

    if (!(AVIdx = (short *)calloc(n,sizeof(short))))  
        goto AAVFin; 
    AVIdxA = n;
    memrq(n,sizeof(short));

    if (!(AVOff = (short *)calloc(n,sizeof(short))))  
        goto AAVFin; 
    AVOffA = n;
    memrq(n,sizeof(short));

    if (!(AVLen = (short *)calloc(n,sizeof(short))))  
        goto AAVFin; 
    AVLenA = n;
    memrq(n,sizeof(short));

    if (!(AVFmt1 = (short *)calloc(n,sizeof(short))))  
        goto AAVFin; 
    AVFmt1A = n;
    memrq(n,sizeof(short));

    if (!(AVFmt2 = (short *)calloc(n,sizeof(short))))  
        goto AAVFin; 
    AVFmt2A = n;
    memrq(n,sizeof(short));

    if (!(AVMBlnk = (int *)calloc(n,sizeof(int))))  
        goto AAVFin; 
    AVMBlnkA = n;
    memrq(n,sizeof(int));

    if (!(AVMStar = (int *)calloc(n,sizeof(int))))  
        goto AAVFin; 
    AVMStarA = n;
    memrq(n,sizeof(int));

    if (!(AVMPnt = (int *)calloc(n,sizeof(int))))  
        goto AAVFin; 
    AVMPntA = n;
    memrq(n,sizeof(int));

    if (!(AVMGen = (int *)calloc(n,sizeof(int))))  
        goto AAVFin; 
    AVMGenA = n;
    memrq(n,sizeof(int));

    i = 0;
    j = idx;
    while (j >= 0) {
        if (VTypA[j] == 2) {
            AVIdx[i] = j;
            l = strlen(VDef[j]) - 2;
            if (l > 0)
                AVFmt1[i] = -l;
            else  
                gerr_exit(40);
            
            if (++i >= n)
                break;
        }
        j = VNxt[j];
    }
    if (i != n)
        gerr_exit(46);
    return(0);

AAVFin:
    if (AVMGenA > 0) {
        free((char *)AVMGen);
        memrq(-AVMGenA,sizeof(int));
        AVMGenA = 0;
    }
    if (AVMPntA > 0) {
        free((char *)AVMPnt);
        memrq(-AVMPntA,sizeof(int));
        AVMPntA = 0;
    }
    if (AVMStarA > 0) {
        free((char *)AVMStar);
        memrq(-AVMStarA,sizeof(int));
        AVMStarA = 0;
    }
    if (AVMBlnkA > 0) {
        free((char *)AVMBlnk);
        memrq(-AVMBlnkA,sizeof(int));
        AVMBlnkA = 0;
    }
    if (AVFmt1A > 0) {
        free((char *)AVFmt1);
        memrq(-AVFmt1A,sizeof(short));
        AVFmt1A = 0;
    }
    if (AVFmt2A > 0) {
        free((char *)AVFmt2);
        memrq(-AVFmt2A,sizeof(short));
        AVFmt2A = 0;
    }
    if (AVLenA > 0) {
        free((char *)AVLen);
        memrq(-AVLenA,sizeof(short));
        AVLenA = 0;
    }
    if (AVOffA > 0) {
        free((char *)AVOff);
        memrq(-AVOffA,sizeof(short));
        AVOffA = 0;
    }
    if (AVIdxA > 0) {
        free((char *)AVIdx);
        memrq(-AVIdxA,sizeof(short));
        AVIdxA = 0;
    }
    if (err)  
        p_err(-2,1);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  check_avar(n,adic)  Search variable description file for the variables  */
/*                      in the list AVIdx, n variables. AVFmt1 contains     */
/*                      -length of variable name.                           */
/*                      Set: AVOff, AVLen, AVFmt1, AVFmt2.                  */
/*                      If adic != 0 recognize labels.                      */
/*                      Return index of data file if OK, -1 if error.       */

int check_avar(int n,int adic)
{
    register int i,j,l,ii;
    register char *p,*q;
    int w1,w2,m,fnni,fnn,fn,off,ll,len,err1,err2,err3;

    if (n <= 0)
        return(0);

    ll = 8;
    for (i = 0; i < n; ++i) {
        l = strlen(VDef[AVIdx[i]]) - 2;
        if (ll < l)
            ll = l;
    }
    printf1("\nVariable  ");
    prnchar(' ',ll - 8,0);
    printf1("File  ");
    prnchar(' ',FDefLen - 4,0);
    printf1("FN  Off  Len  Format  Label\n");
    prnchar('-',32 + ll + FDefLen,1);

    dbf_init(1);  /* init buffer for reading the var description file */

    /*  read the variable description file */

    err1 = err2 = err3 = m = 0;
    fnni = fnn = -1;

    while (get_record(VFN,ZABuf,ZABLen) > 0) {

        if ((p = check_comment(ZABuf)) != NULL) { /* skip comment lines */
NXT1:
            if (*ZABuf == ' ')
                goto NXT;

            for (i = 0; i < n; ++i) {

                if (AVFmt1[i] < 0.0) {

                    l = (int)(-AVFmt1[i]);
                    j = AVIdx[i];
                    q = VDef[j] + 2;

                    if (!strncmp(p,q,l) && *(p + l) == ' ') {

                        printf1("%s  ",q);
                        prnchar(' ',ll - l,0);

                        q = p + strlen(p);
                        while (q > p && (*--q == '\n' || *q == LF || *q == CR || *q == ' ')) ;
                        *++q = '\0';

                        w1 = w2 = 0;
                        p = skip_cb(p);
                        if (sscanf(p,"%d %d %d.%d",&fn,&off,&w1,&w2) == 4) ; 
                        else if (sscanf(p,"%d %d %d ",&fn,&off,&w1) == 3)   
                            w2 = 0;
                        else {
                            printf1(" cannot get information about variable\n");
                            fn = -1;
                        }

                        if (fnn >= 0 && fn == fnn)
                            j = fnni;
                        else {
                            j = -1;
                            for (ii = 0; ii < ZANF; ++ii) {
                                if (fn == (int)ZAFNum[ii]) {
                                    j = ii;
                                    break;
                                }
                            }
                            if (fnn < 0 && j >= 0) {
                                fnn = fn;
                                fnni = j;
                            }
                        }
                        if (j < 0) {
                            printf1("???? ");
                            prnchar(' ',FDefLen - 4,0);
                            err2 = 1;
                        }
                        else {
                            printf1("%s ",ZAFNam[j]);
                            prnchar(' ',FDefLen - strlen(ZAFNam[j]),0);
                            if (fn != fnn)
                                err3 = 1;
                        }
                        AVLen[i] = len = iabs(w1);

                        if (fn < 0 || off < 0 || len < 1 || w2 < 0 || w2 >= len)
                            err1 = 1;

                        if (w1 < 0) {       /* string variable */
                            VTyp[AVIdx[i]] = 1;
                            VSLen[AVIdx[i]] = w1;
                            w1 = w2 = 0;
                            if (len > SVBufLen) {
                                if (svb_alloc(len))
                                    return(-1);
                            }
                        }   
                        else {
                            VTyp[AVIdx[i]] = 3;
                            VSLen[AVIdx[i]] = get_afmt(&w1,&w2);
                        }
                        printf1("%3d %5d %4d %4d.%-3d ",fn,off,len,w1,w2);
                        AVOff[i] = off;
                        AVFmt1[i] = w1; 
                        AVFmt2[i] = w2; 

                        p = skip_cb(p);
                        p = skip_cb(p);
                        p = skip_cb(p);
                        if (*p)
                            printf1("%s\n",p);
                        else
                            printf1("\n");
                        m++;

                        if (adic) {
                            while (get_record(VFN,ZABuf,ZABLen) > 0) {
                                if ((p = check_comment(ZABuf)) != NULL) {
                                    if (*ZABuf == ' ') {

                                        q = p + strlen(p);
                                        while (q > p && (*--q == '\n' || *q == LF || *q == CR || *q == ' ')) ;
                                        *++q = '\0';
                                        printf1(" %s\n",p);

                                    }
                                    else if (m == n)
                                        break;
                                    else
                                        goto NXT1;
                                }
                            }
                        }
                        break;     
                    }
                }
                if (m == n)
                    break;
            }
        }
NXT:
        if (m == n)
            break;
    }
    if (m < n) {
        for (i = 0; i < n; ++i) {
            if (AVFmt1[i] < 0.0) {
                j = AVIdx[i];
                q = VDef[j] + 2;

                printf1("%s  ",q);
                prnchar(' ',ll + AVFmt1[i],0);
                printf1("can't find this variable.\n");
            }
        }
        return(-1);
    }
    if (err1) {
        printf1("\nIncorrect entries for at least one variable.\n");
        return(-1);
    }
    else if (err2 || fnni < 0) {
        printf1("\nError: there are unkown data files.\n");
        return(-1);
    }
    else if (err3) {
        printf1("\nError: variables do not belong to the same data file.\n");
        return(-1);
    }
    printf1("\nUsing archive data file: %s\n",ZAFNam[fnni]);
    printf1("Number of records: %d. Record length: %d",
                                         ZAFNRec[fnni],(int)ZAFRLen[fnni]);
    if ((int)ZAFRLen[fnni] == 0)  
        printf1(" (variable)");
    printf1(".\n");
    return(fnni);
}

/* ------------------------------------------------------------------------ */
/*  get_avar()      Read next record from archive file AVDFN and put values */
/*                  of requested NVArc variables into AVVAL[].              */
/*                  Also count missings in AVMBlnk, and so on.              */
/*                                                                          */
/*                  Return -1 if error, 0 if EOF, 1 if OK.                  */

int get_avar(void)
{
    register int i,len;
    int mval,rlen;
    double tmp;
    register char *p,*q,*s;

    rlen = get_drec(AVDFN,ZABuf,ZABLen);
    if (rlen <= 0)
        return(rlen);

    for (i = 0; i < NVArc; ++i) {

        if (VTyp[AVIdx[i]] == 1)    /* string variable */
            tmp = 0.0;
        else {
            q = ZABuf;     /* pointer to buffer with current record */

            len = AVLen[i];
            p = q + AVOff[i];
            s = p + len;

            if (s > q + rlen) {
                tmp = MBlnkVal;
                AVMBlnk[i] += 1;
            }
            else {
                while (p < s && *p == ' ') {    /* skip leading blanks */
                    p++;
                    len--;
                }
                if (len <= 0) {
                    tmp = MBlnkVal;
                    AVMBlnk[i] += 1;
                }
                else {
                    tmp = dscan(p,len,&mval);
                    if (mval) {
                        if (mval == 1)
                            AVMBlnk[i] += 1;
                        else if (mval == 2)
                            AVMStar[i] += 1;
                        else if (mval == 3)
                            AVMPnt[i] += 1;
                        else 
                            AVMGen[i] += 1;
                    }
                }
            }
        }
        AVVAL[AVIdx[i]] = tmp;
    }
    return(1);
}

/*--------------------------------------------------------------------------*/
/*  get_astr(buf,i)                                                         */
/*                                                                          */
/*  Get string from archive variable i into buf.                            */
          
void get_astr(char *buf,int i)
{
    register int j;
    register char *p,*q,*s;
    int len;

    p = ZABuf + AVOff[i];     /* pointer to buffer with current record */
    q = buf;
    s = ZABuf + ZABLen;
    len = AVLen[i];

    for (j = 0; j < len; ++j) {
        if (p < s)
            *q++ = *p++;
        else
            *q++ = ' ';
    }
    *q = '\0';
}

/* ------------------------------------------------------------------------ */
/*  arcc()      Test the archive and its description (arcc command).        */
/*              It is tried to read all files in the archive. The           */
/*              number of records read is checked against the number        */
/*              of records defined.                                         */
/*              Return 0 if OK, or -1 if errors.                            */

int arcc(void)
{
    register int i,j,n;
    int err,err1,terr,nrec,vdfi,vdfn,iv;
    register char *p,*q; 

    err = vdfn = vdfi = -1;

    if (check_cmd(1))
        return(-1);

    if (ARCDef == 0) {
        p_err(-9,1);
        return(-1);
    }
    printf1("Archive check. Reading files defined in: %s\n\n",ZADNam);

    terr = 0;
    for (i = 0; i < ZANF; ++i) {      /* read all files in the archive */

        printf1("File:%4d  %s",ZAFNum[i],ZAFNam[i]);
        prnchar(' ',FDefLen - strlen(ZAFNam[i]),0);
        /* fflush(stdout); */
  
        if (ZAFTyp[i] == 2 && vdfi < 0) {
            vdfi = i;
            vdfn = ZAFNRec[i];
            iv = 0;
            if (vdfn > 0) {
                if (alloc_acc((VNLMax + 1) * vdfn))
                    goto ACCFin;
            }
        }
        dbf_init(1);
        /*****
        ZAFReq[i] = 1;
        ****/
        err1 = nrec = 0;
        while ((n = get_record(i,ZABuf,ZABLen)) > 0) {
            if (ZAFRLen[i]) {
                if (n != ZAFRLen[i]) {
                    printf1(" error in record length or number of records\n");
                    err1 = 1;
                    break;
                }
            }
            if (i == vdfi) {
                if ((p = check_comment(ZABuf)) != NULL) {
                    if (*ZABuf != ' ') {
                        q = AcC + iv * (VNLMax + 1);
                        for (j = 0; j < VNLMax; ++j) {
                            if (*p && *p != ' ')
                                *q++ = *p++;
                            else
                                *q++ = ' ';
                        }
                        *q = '\0';
                        iv++;
                    }
                }
            }
            nrec++;
        }
        if (!err1) {
            printf1(" records: %7d\n",nrec);
            if (nrec != ZAFNRec[i]) {
                printf1("Inconsistent with defined number of records (%d)\n",
                                                                    ZAFNRec[i]);
                terr++;
            }
        }
        else
            terr++;
    }
    if (terr) {
        printf1("Found %d error(s) in archive %s\n",terr,ZOONam);
        goto ACCFin;
    }

    /* check variable description file */

    if (vdfi < 0) {
        printf1("Error: no variable description file.\n");
        goto ACCFin;
    }
    printf1("Found %d variables.\n",iv);

    qsort(AcC,iv,VNLMax + 1,iv_comp);

    for (i = 1; i < iv; ++i) {
        if (!strcmp(AcC + i * (VNLMax + 1),AcC + (i - 1) * (VNLMax + 1))) {
            printf1("Error: at least one variable name is not unique.\n");
            goto ACCFin;
        }
    }
    printf1("No errors found.\n");
    err = 0;
ACCFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  iv_comp()                                                               */

int iv_comp(const void *arg1,const void *arg2)
{     
    return(strcmp((char *)arg1,(char *)arg2));
}

/* ------------------------------------------------------------------------ */
/*  arcv()      arcv(fn=...)=fname. Print archive variables.                */
/*              Return 0 or -1 if syntax error.                             */

int arcv(void)
{
    register int i,j,k;
    int err,fnd,n,nn,fn,fna,fnn,off,cflag,w1,w2,rec,vtyp,vslen;
    register char *p,*q;

    err = -1;
    if (check_cmd(1))
        return(-1);

    if (ARCDef == 0) {
        p_err(-9,1);
        return(-1);
    }
    if (parm(CmdBuf + 4,1,0))     /* get parameters */
        goto PAFin;   

    printf1("Using variable description file: %s\n",ZAFNam[VFN]);

    err = n = 0;
    if (PMFN > 0) {

        if (alloc_acn(PMFN))
           goto PAFin;

        for (k = 0; k < PMFN; ++k) {
            j = -1;
            for (i = 0; i < ZANF; ++i) {
                if (!strcmp(PMFNam[k],ZAFNam[i]) && ZAFTyp[i] == 1) {
                    j = ZAFNum[i];
                    AcN[k] = j;
                    n++;
                    break;
                }
            }
            if (j < 0)  
                printf1("Data file %s not in archive (ignored).\n",PMFNam[k]);
        }
        if (n == 0) {
            printf1("No data files.\n");
            goto PAFin;
        }
    }

    dbf_init(1);    /* init reading of var description file */

    if (PMFDef)
        fprintf(PMFd,"# arcd = %s;\n# nvar( \n",ZADNam);
    else
        printf1("# arcd = %s;\n# nvar(\n",ZADNam);

    nn = 1;
    rec = fnd = 0;
    fnn = fna = -1;

    while (get_record(VFN,ZABuf,ZABLen) > 0) {

        rec++;
        if (*ZABuf == ' ')
            cflag = 1;
        else
            cflag = 0;

        if ((p = check_comment(ZABuf)) != NULL) {

            if (cflag && (PMArcDic == 0 || fnd == 0))
                continue; 

            q = p + strlen(p);
            while (q > p && (*--q == '\n' || *q == LF || *q == CR || *q == ' ')) ;
            *++q = '\0';

            if (fnd) {
                if (cflag && PMArcDic) {
                    if (PMFDef)
                        fprintf(PMFd,"#  %s\n",p);
                    else
                        printf1("#  %s\n",p);
                    nn++;
                    continue; 
                }
                fnd = 0;
            }
            q = skip_cb(p);
                                 
            if (sscanf(q,"%d %d %d.%d",&fn,&off,&w1,&w2) == 4) ;
            else if (sscanf(q,"%d %d %d.%d",&fn,&off,&w1,&w2) == 3)  
                w2 = 0;
            else  
                w1 = 0;  

            if (fn < 0 || off < 0 || iabs(w1) < 1 || w2 < 0 || w2 >= iabs(w1)) {
                printf1("Error in record %d of variable description file.\n",rec);
                err = -1;
                goto PAFin;
            }
            if (w1 < 0) {
                vtyp = 1;
                vslen = -w1;
                w1 = w2 = 0;
            }
            else {
                vtyp = 3;
            }
            j = 1;
            if (n > 0) {
                j = 0;
                for (i = 0; i < PMFN; ++i) {
                    if (AcN[i] == fn) {
                        j = 1;
                        break;
                    }
                }
            }
            if (j) {
                *(q - 1) = '\0';
                q = skip_cb(q);
                q = skip_cb(q);
                q = skip_cb(q);

                if (vtyp == 1)
                    j = -vslen;
                else
                    j = get_afmt(&w1,&w2);

                if (fnn != fn) {
                    for (k = 0; k < ZANF; ++k) {
                        if (ZAFNum[k] == fn) {
                            fnn = fn;
                            fna = k;
                            break;
                        }
                    }
                }
                if (fna < 0)
                    fna = 0;

                if (PMFDef)
                    fprintf(PMFd,"# %s<%d>[%d.%d] = A:%s, # [%s] %s\n",p,j,w1,w2,p,ZAFNam[fna],q);
                else
                    printf1("# %s<%d>[%d.%d] = A:%s, # [%s] %s\n",p,j,w1,w2,p,ZAFNam[fna],q);
                nn++;
                fnd = 1;
            }
        }
    }
    if (PMFDef)
        fprintf(PMFd,"# );\n");
    else
        printf1("# );\n");

    if (PMFDef && nn > 1)  
        printf1("%d records written to: %s\n",nn,PMFdName);

PAFin:
    p_clean();
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  arcvc   Check var description file.                                     */
/*                                                                          */
/*          arcvc(                                                          */
/*              df = ...,                                                   */
/*              opt = ...,  option, def. 1                                  */
/*                          1 = append file number                          */
/*                          2 = append file name (lower case)               */
/*                          3 = append file name (upper case)               */
/*          ) = file_name;                                                  */  
/*                                                                          */
/*          If df=... create new file with unique var names.                */
/*          Return 0 if OK, -1 if error.                                    */

int arcvc(void)
{
    register int i,j,l;
    int n,nn,m,fnd,nv,err,rec,vptra,cnta,fnna,fn,off,w1,w2,flen;
    short *count,*fnn;
    char fname[17],*p,*q,**vptr;

    err = -1;
    nn = rec = fnna = cnta = vptra = 0;

    if (check_cmd(1))
        return(-1);

    if (parm(CmdBuf + 5,2,1))     /* get parameters */
        goto PACFin;   

    printf1("Checking variable description file: %s\n",PMFdName);

    if (PMOPT > 3)
        PMOPT = 3;

    if (alloc_acc(RLMaxDef))
        goto PACFin;

    nv = 0;     /* count number of variables */

    while (fgets(AcC,RLMaxDef,PMFd)) {
        if (*AcC != ' ' && (p = check_comment(AcC)) != NULL)  
            nv++;
    }

    if (!(vptr = (char **) calloc(nv + 1,sizeof(char *)))) {
        p_err(-2,1);
        goto PACFin;
    }
    memrq(nv + 1,sizeof(char *));
    vptra = 1;

    if (!(count = (short *) calloc(nv + 1,sizeof(short)))) {
        p_err(-2,1);
        goto PACFin;
    }
    memrq(nv + 1,sizeof(short));
    cnta = 1;

    if (!(fnn = (short *) calloc(nv + 1,sizeof(short)))) {
        p_err(-2,1);
        goto PACFin;
    }
    memrq(nv + 1,sizeof(short));
    fnna = 1;

    fseek(PMFd,0,0);
    n = nn = rec = err = 0;
       
    while (fgets(AcC,RLMaxDef,PMFd)) {

        rec++;
        prn_message(rec,0,0);

        if (*AcC != ' ' && (p = check_comment(AcC)) != NULL) {

            p = skip_c(AcC);
            q = skip_b(p + 1);
            if (sscanf(q,"%d %d %d.%d",&fn,&off,&w1,&w2) == 4) ;
            else if (sscanf(q,"%d %d %d",&fn,&off,&w1) == 3)
                w2 = 0;
            else
                err = 1;

            if (err == 0) {
                if (fn < 1)
                    err = 2;
                else if (off < 0)
                    err = 3;
                else if (iabs(w1) < 1)
                    err = 4;
                else if (w2 < 0 || w2 >= iabs(w1))
                    err = 5;
            }
            if (err)  
                goto PACFin;

            *p = '\0';
            fnd = 0;
            for (j = 0; j < nn; ++j) {
                if (!strcmp(vptr[j],AcC)) {
                    count[j] += 1;
                    fnd = 1;
                    break;
                }
            }
            if (!fnd) {
                l = strlen(AcC) + 1;
                if (!(vptr[nn] = (char *) calloc(l,sizeof(char)))) {
                    p_err(-2,1);  
                    err = -1;
                    goto PACFin;
                }
                memrq(l,1);

                strcpy(vptr[nn],AcC);
                count[j] = 1;
                nn++;
            }
            fnn[n++] = fn;
        }
    }
    prn_message(rec,1,0);
    printf1("Found %d variables.\n",nv);
    j = 0;
    for (i = 0; i < nn; ++i) {
        if (count[i] > 1) {
            printf1("Found %2d times: %s\n",count[i],vptr[i]);
            j++;
        }
    }
    printf1("Number of variable names used more than once: %d\n",j);

    if (j > 0 && PMF1Def) {         /* create new file */

        printf1("Creating a new file with unique variable names.\n");

        fseek(PMFd,0,0);
        flen = n = rec = 0;
       
        while (fgets(AcC,RLMaxDef,PMFd)) {

            rec++;
            prn_message(rec,0,0);

            if (PMOPT > 1) {
                if (!strncmp(AcC,"# variable",10)) {
                    p = AcC + 38;
                    q = fname;
                    for (j = 0; j < 16; ++j) {
                        *q++ = *p++;
                        if (!*p || *p == '.' || *p == LF || *p == CR)
                            break;
                    }
                    *q = '\0';

                    if (PMOPT == 3) {
                        p = fname;
                        for (j = 0; j < 16; ++j) {
                            if (!*p)
                                break;
                            if (isalpha((int)*p))  
                                *p = (char)toupper((int)*p);
                            p++;
                        }
                    }
                    flen = strlen(fname);
                }
                if (flen == 0) {
                    printf1("Error: cannot create new file names with opt = %d.\n",PMOPT);
                    err = -1;
                    goto PACFin;
                }
            }

            m = 0;
            if (*AcC != ' ' && (p = check_comment(AcC)) != NULL) {

                p = skip_c(AcC);
                *p = '\0';

                m = 0;
                for (j = 0; j < nn; ++j) {
                    if (!strcmp(vptr[j],AcC)) {
                        m = count[j] - 1;
                        break;
                    }
                }
                if (m > 0) {
                    if (PMOPT == 1)
                        fprintf(PMF1d,"%s_%d %s",AcC,fnn[n],p + 1);
                    else
                        fprintf(PMF1d,"%s_%s %s",AcC,fname,p + 1);
                }
                *p = ' ';
                n++;
            }
            if (m == 0)
                fprintf(PMF1d,"%s",AcC);
        }
        prn_message(rec,1,0);
        printf1("%d records (%d variables) written to: %s\n",rec,n,PMF1dName);
    }
    err = 0;

PACFin:
    if (err > 0) {
        printf1("Error in record %d: %s\n",rec,AcC);
        if (err == 1)
            printf1("Need file number, offset and a correct format.\n");
        else if (err == 2)
            printf1("Logical file number must be positive.\n");
        else if (err == 3)
            printf1("Offset must not be negative.\n");
        else if (err == 4)
            printf1("Width must be at least one column.\n");
        else if (err == 5)
            printf1("Error in format specification.\n");
        err = -1;
    }
    for (i = 0; i < nn; ++i) {
        l = strlen(vptr[i]) + 1;
        free(vptr[i]);
        memrq(-l,1);
    }
    if (vptra) {  
        free((char *)vptr);
        memrq(-nv - 1,sizeof(char *));
    }
    if (cnta) {
        free((char *)count);
        memrq(-nv - 1,sizeof(short));
    }
    if (fnna) {
        free((char *)fnn);
        memrq(-nv - 1,sizeof(short));
    }
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_afmt(w1,w2)     create default print format and return storge size. */
/*                      not for string variables.                           */

int get_afmt(int *w1,int *w2)
{
    int s;

    if (*w1 < 2)
        *w1 = 2;

    if (*w2 == 0) {      /* integer */
        if (*w1 <= 2)
            s = 1;
        else if (*w1 <= 4)
            s = 2;
        else if (*w1 <= 9)
            s = 5;
        else            /* should be double */
            s = 8;
    }
    else {              /* float */
        if (*w1 <= 7)
            s = 4;
        else
            s = 8;
    }
    if (*w2 > 0 && *w1 < *w2 + 3)  /* adjust print format */
        *w1 = *w2 + 3;

    return(s);
}

