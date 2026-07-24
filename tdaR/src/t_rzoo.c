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
#include "tda_context.h"
#include "tda_compat.h"

/*  functions in t_rzoo.c */

void arcd_off(TDAContext *ctx);
int arcd(TDAContext *ctx);
int arcd_alloc(TDAContext *ctx, int opt);
int get_len(TDAContext *ctx, char *p);
char *cpy_nam(TDAContext *ctx, char *s, char *p, int n);
void prn_afiles(TDAContext *ctx);
int get_drec(TDAContext *ctx, int fn, char *buf, int nmax);
int get_record(TDAContext *ctx, int fn, char *buf, int nmax);
void get_rdat(TDAContext *ctx, int fn);
int alloc_avar(TDAContext *ctx, int idx,int n,int opt);
int check_avar(TDAContext *ctx, int n,int adic);
int get_avar(TDAContext *ctx);
void get_astr(TDAContext *ctx, char *buf,int i);
int arcc(TDAContext *ctx);
int arcv(TDAContext *ctx);
int arcvc(TDAContext *ctx);
int get_afmt(TDAContext *ctx, int *w1,int *w2);
int iv_comp(const void *, const void *, void *);

/* ------------------------------------------------------------------------ */
/*  Global variables                                                        */



/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------ */
/*  For each file we allocate a read buffer where the decompression is      */
/*  done, and a line buffer.                                                */


/* ------------------------------------------------------------------------ */
/*  arcd_off    Turn off arrchive.                                          */

void arcd_off(TDAContext *ctx)
{
    if (ctx->ARCDef) {
        dbf_init(ctx, 0);
        arcd_alloc(ctx, 0);
        ctx->ARCDef = ctx->ZANF = 0;
        if (ctx->ZOOFdo) {
            fclose(ctx->ZOOFd);
            ctx->ZOOFdo = 0;
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

int arcd(TDAContext *ctx)
{
    FILE *fd;
    register int i,j;
    register char *p,*q;
    int err,m,n,first,rec = 0,n1,n2,n3,n4;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    if (!strcmp(ctx->CmdBuf,"arcd")) {
        if (ctx->ARCDef == 0)
            printf1(ctx, "Data archive not defined.\n");
        else {
            printf1(ctx, "Archive: %s\n",ctx->ZOONam);
            printf1(ctx, "Description file: %s\n",ctx->ZADNam);
            prn_afiles(ctx);
        }
        return(0);
    }

    if (ctx->ARCDef) {
        dbf_init(ctx, 0);
        arcd_alloc(ctx, 0);
        ctx->ARCDef = ctx->ZANF = 0;
        if (!strcmp(ctx->CmdBuf,"arcd=off")) {
            printf1(ctx, "Access to data archive %s turned off.\nCurrent memory: %d bytes.\n",
                                                             ctx->ZOONam,ctx->MemReq);
            return(0);
        }
    }
    else if (!strcmp(ctx->CmdBuf,"arcd=off")) {
        printf1(ctx, "No data archive defined (command ignored).\n");
        return(0);
    }
    strncpy(ctx->ZADNam,ctx->CmdBuf + 5,FNMaxLen - 1);
    ctx->ZADNam[FNMaxLen - 1] = '\0';

    if (!(fd = fopen(ctx->ZADNam,OPEN_RD))) {
        printf1(ctx, "Can't open archive description file: %s\n",ctx->ZADNam);
        return(-1);
    }
    printf1(ctx, "Reading archive description file: %s\n",ctx->ZADNam);
        
    if (alloc_acc(ctx, RLMaxDef))
        goto AFin;

    /*  count number of files in description file in ZANF */

    ctx->ZANF = 0;
    while (fgets(ctx->AcC,RLMaxDef,fd)) {
        if ((p = check_comment(ctx, ctx->AcC)) != NULL)  
            ctx->ZANF++;
    }
    if (--ctx->ZANF < 1) {
        printf1(ctx, "Error: no definition of files in: %s\n",ctx->ZADNam);
        goto AFin;
    }
    if (arcd_alloc(ctx, 1)) {
        printf1(ctx, "Insufficient memory for archive data structures.\n");
        goto AFin;
    }
    fseek(fd,0L,0);

    ctx->ZABLen = RLMaxDef;  /* max buffer size for data files, updated below */
    rec = m = 0;
    first = 1;
    err = -2;

    while (m < ctx->ZANF && fgets(ctx->AcC,RLMaxDef,fd)) {

        rec++;
        if ((p = check_comment(ctx, ctx->AcC)) != NULL) {

            q = p + strlen(p);
            while (--q > p && (*q == '\n' || *q == LF || *q == CR)) ;
            *++q = '\0';
   
            if (first) {              
                n = get_len(ctx, p);
                if (n >= FNMaxLen)  
                    goto AFin;
                cpy_nam(ctx, ctx->ZOONam,p,n);
                first = 0;
            }
            else {       
                if (!*p || sscanf(p,"%d",&n1) != 1 || n1 < 0)  
                    goto AFin;
                ctx->ZAFNum[m] = (short)n1;
                p = skip_cb(ctx, p);
                n = get_len(ctx, p);
                if (n >= FNMaxLen)  
                    goto AFin;
                if (ctx->FDefLen < n)
                    ctx->FDefLen = n;

                p = cpy_nam(ctx, ctx->ZAFNam[m],p,n);
                p = skip_b(ctx, p);
                n1 = n2 = n3 = n4 = -1;
                if (!*p || sscanf(p,"%d %d %d %d",&n1,&n2,&n3,&n4) != 4 ||
                    (n1 != 1 && n1 != 2) || n2 < 0 || n3 < 1 || n4 < 0) {
                    goto AFin;
                }
                ctx->ZAFTyp[m]  = (short)n1;
                ctx->ZAFRLen[m] = (short)n2;
                ctx->ZAFNRec[m] = n3;
                ctx->ZAFNVar[m] = (short)n4;    
                     
                if (ctx->ZABLen < ctx->ZAFRLen[m])
                    ctx->ZABLen = ctx->ZAFRLen[m];
                m++;
            }
        }
    }
    if (m != ctx->ZANF) {
        err = -2;
        goto AFin;
    }
    rec = 0;
    if (first) {
        printf1(ctx, "Error: can't find name of ZOO archive.\n");
        goto AFin;
    }
    printf1(ctx, "ZOO data archive: %s\n",ctx->ZOONam);

    /*  Open the ZOO archive and get header and directory information.      */

    if (!(ctx->ZOOFd = fopen(ctx->ZOONam,OPEN_RB))) {
        printf1(ctx, "Can't open this file.\n");
        goto AFin;
    }
    ctx->ZOOFdo = 1;

    if (get_zoo(ctx))  
        goto AFin;

    printf1(ctx, "Checking definition of files in archive.\n");
    err = m = n = 0;

    for (i = 0; i < ctx->ZANF; ++i) {
            
        if (!ctx->ZAZOfs[i]) {
            printf1(ctx, "Can't find file %s in the archive.\n",ctx->ZAFNam[i]);
            err = -1;
        }
        /*  0 (stored) accepted as well as 1 (LZD) and 2 (LZH).  TDA only
            ever consumed archives made by Dhesi's zoo, which compresses,
            so it had never met a stored member and refused it as an
            "unknown packing method" -- while every other zoo reader takes
            one.  lzs() in t_zoo.c handles it.  */
        else if (ctx->ZAZTyp[i] < 0 || ctx->ZAZTyp[i] > 2) {
            printf1(ctx, "File %s has unknown packing method.\n",ctx->ZAFNam[i]);
            err = -1;
        }
        for (j = 0; j < i; ++j) {
            if (ctx->ZAFNum[i] == ctx->ZAFNum[j])  
                n++;
        }
        if (ctx->ZAFTyp[i] == 2) {
            ctx->VFN = i;
            m++;
        }
    }
    if (n > 0) {
        printf1(ctx, "Logical file ID numbers are not unique.\n");
        err = -1;
    }
    if (m == 0) {
        printf1(ctx, "Can't find a variable description file (type 2).\n");
        err = -1;
    }
    else if (m > 1) {
        printf1(ctx, "Found more than one variable description file (type 2).\n");
        err = -1;
    }
    prn_afiles(ctx);           /* print list of files */

    if (err)
        goto AFin;

    ctx->ARCDef = 1;     /* set global flag for successfully opened archive */
    err = 0;      

AFin:
    fclose(fd);
    alloc_acc(ctx, 0);

    if (err) {
        if (err == -2 && rec > 0)  
            printf1(ctx, "Error in record %d of archive description file.\n",rec);
        arcd_alloc(ctx, 0);
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  arcd_alloc(opt)     If opt != 0 allocate data structures for ZANF       */
/*                      data files, otherwise free the memory.              */
/*                      Return 0 if OK, -1 if error.                        */

int arcd_alloc(TDAContext *ctx, int opt)
{
    register int i;
    int err;

    if (opt == 0) {
        err = 0;
        goto AAFin;
    }   
    err = -1;

    if (!(ctx->ZAFNam = (char **) calloc((size_t)(ctx->ZANF),sizeof(char *))))  
        goto AAFin; 
    ctx->ZAFNamA = ctx->ZANF;
    memrq(ctx, ctx->ZAFNamA,sizeof(char *));
   
    if (!(ctx->ZAFNum = (short *) calloc((size_t)(ctx->ZANF),sizeof(short))))  
        goto AAFin; 
    ctx->ZAFNumA = ctx->ZANF;
    memrq(ctx, ctx->ZAFNumA,sizeof(short));

    if (!(ctx->ZAFTyp = (short *) calloc((size_t)(ctx->ZANF),sizeof(short))))  
        goto AAFin; 
    ctx->ZAFTypA = ctx->ZANF;
    memrq(ctx, ctx->ZAFTypA,sizeof(short));

    if (!(ctx->ZAFRLen = (short *) calloc((size_t)(ctx->ZANF),sizeof(short))))  
        goto AAFin; 
    ctx->ZAFRLenA = ctx->ZANF;
    memrq(ctx, ctx->ZAFRLenA,sizeof(short));

    if (!(ctx->ZAFNVar = (short *) calloc((size_t)(ctx->ZANF),sizeof(short))))  
        goto AAFin; 
    ctx->ZAFNVarA = ctx->ZANF;
    memrq(ctx, ctx->ZAFNVarA,sizeof(short));

    if (!(ctx->ZAFNRec = (int *) calloc((size_t)(ctx->ZANF),sizeof(int))))  
        goto AAFin; 
    ctx->ZAFNRecA = ctx->ZANF;
    memrq(ctx, ctx->ZAFNRecA,sizeof(int));

    if (!(ctx->ZAZSiz = (int *) calloc((size_t)(ctx->ZANF),sizeof(int))))  
        goto AAFin; 
    ctx->ZAZSizA = ctx->ZANF;
    memrq(ctx, ctx->ZAZSizA,sizeof(int));

    if (!(ctx->ZAZTyp = (char *) calloc((size_t)(ctx->ZANF),sizeof(char))))  
        goto AAFin; 
    ctx->ZAZTypA = ctx->ZANF;
    memrq(ctx, ctx->ZAZTypA,sizeof(char));

    if (!(ctx->ZAFReq = (char *) calloc((size_t)(ctx->ZANF),sizeof(char))))  
        goto AAFin; 
    ctx->ZAFReqA = ctx->ZANF;
    memrq(ctx, ctx->ZAFReqA,sizeof(char));

    if (!(ctx->ZAZOfs = (int *) calloc((size_t)(ctx->ZANF),sizeof(int))))  
        goto AAFin; 
    ctx->ZAZOfsA = ctx->ZANF;
    memrq(ctx, ctx->ZAZOfsA,sizeof(int));

    for (i = 0; i < ctx->ZANF; ++i) {
        if (!(ctx->ZAFNam[i] = (char *) calloc(FNMaxLen,sizeof(char))))  
            goto AAFin;
        memrq(ctx, FNMaxLen,1);
        ctx->ZAFNamAI++;
    }
    return(0);

AAFin:
    if (ctx->ZAFNamAI > 0) {
        for (i = 0; i < ctx->ZAFNamAI; ++i) {
            free(ctx->ZAFNam[i]);
            memrq(ctx, -FNMaxLen,1);
        }
        ctx->ZAFNamAI = 0;
    }
    if (ctx->ZAZOfsA > 0) {
        free((char *)ctx->ZAZOfs);
        memrq(ctx, -ctx->ZAZOfsA,sizeof(int));
        ctx->ZAZOfsA = 0;
    }      
    if (ctx->ZAFReqA > 0) {
        free((char *)ctx->ZAFReq);
        memrq(ctx, -ctx->ZAFReqA,sizeof(char));
        ctx->ZAFReqA = 0;
    }   
    if (ctx->ZAZTypA > 0) {
        free((char *)ctx->ZAZTyp);
        memrq(ctx, -ctx->ZAZTypA,sizeof(char));
        ctx->ZAZTypA = 0;
    }   
    if (ctx->ZAZSizA > 0) {
        free((char *)ctx->ZAZSiz);
        memrq(ctx, -ctx->ZAZSizA,sizeof(int));
        ctx->ZAZSizA = 0;
    }   
    if (ctx->ZAFNRecA > 0) {
        free((char *)ctx->ZAFNRec);
        memrq(ctx, -ctx->ZAFNRecA,sizeof(int));
        ctx->ZAFNRecA = 0;
    }
    if (ctx->ZAFNVarA > 0) {
        free((char *)ctx->ZAFNVar);
        memrq(ctx, -ctx->ZAFNVarA,sizeof(short));
        ctx->ZAFNVarA = 0;
    }
    if (ctx->ZAFRLenA > 0) {
        free((char *)ctx->ZAFRLen);
        memrq(ctx, -ctx->ZAFRLenA,sizeof(short));
        ctx->ZAFRLenA = 0;
    }
    if (ctx->ZAFTypA > 0) {
        free((char *)ctx->ZAFTyp);
        memrq(ctx, -ctx->ZAFTypA,sizeof(short));
        ctx->ZAFTypA = 0;
    }
    if (ctx->ZAFNumA > 0) {
        free((char *)ctx->ZAFNum);
        memrq(ctx, -ctx->ZAFNumA,sizeof(short));
        ctx->ZAFNumA = 0;
    }
    if (ctx->ZAFNamA > 0) {
        free((char *)ctx->ZAFNam);
        memrq(ctx, -ctx->ZAFNamA,sizeof(char *));
        ctx->ZAFNamA = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_len(p)    Get length of string at p until next blank or tab char.   */

int get_len(TDAContext *ctx, char *p)
{
    (void)ctx;        /* unused: the signature is shared */
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

char *cpy_nam(TDAContext *ctx, char *s, char *p, int n)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;
    for (i = 0; i < n; ++i)
        *s++ = *p++;
    *s = '\0';
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  prn_afiles()  Print list of files described in archive description file */

void prn_afiles(TDAContext *ctx)
{
    register int i;

    printf1(ctx, "\nFN  Type  RLen  Records  NVar    Size    M  Name\n");
    prnchar(ctx, '-',48,1); 
                        
    for (i = 0; i < ctx->ZANF; ++i) {

        printf1(ctx, "%3d %3d %6d %8d %5d ",
                    ctx->ZAFNum[i],ctx->ZAFTyp[i],ctx->ZAFRLen[i],ctx->ZAFNRec[i],ctx->ZAFNVar[i]);
        printf1(ctx, "%9d %2d  %s\n",ctx->ZAZSiz[i],(int)ctx->ZAZTyp[i],ctx->ZAFNam[i]);
    }
    newline(ctx);
}

/* ------------------------------------------------------------------------ */
/*  get_drec(fn,buf,nmax)                                                   */
/*      Reads a record from file number fn into the buffer buf. Return the  */
/*      number of bytes in the buffer, or zero if EOF, or -1 if error.      */
/*      nmax is the maximum buffer length.                                  */
/*                                                                          */
/*      EOL characters are NOT counted and not put into the buffer; empty   */
/*      lines are skipped.                                                  */

int get_drec(TDAContext *ctx, int fn, char *buf, int nmax)
{
    register int n;
    register char *p;

    while ((n = get_record(ctx, fn,buf,nmax)) > 0) {
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

int get_record(TDAContext *ctx, int fn, char *buf, int nmax)
{
    register int n,i,j = 0;
    register char c,*q;

    q = buf;
    if (ctx->ZAFRLen[fn] <= 0) {
        while (j < nmax) {
            while (ctx->Out_Buf_Cnt > 0) {
                c = *ctx->Out_Buf_Ptr++;
                ctx->Out_Buf_Cnt--;
                if (c == CR || c == LF) {
                    if (ctx->Out_Buf_Cnt <= 0) 
                        get_rdat(ctx, fn);
                    if (ctx->Out_Buf_Cnt > 0) {
                        if (c == CR && *ctx->Out_Buf_Ptr == LF) {
                            ctx->Out_Buf_Ptr++;
                            ctx->Out_Buf_Cnt--;
                        }   
                    }
                    *q++ = '\n';
                    *q = '\0';
                    return(++j);
                }
                else {
                    *q++ = c;
                    if (++j >= nmax) {
                        printf1(ctx, "\nError in retrieving file: %s\n",ctx->ZAFNam[fn]); 
                        printf1(ctx, "Can't find an end-of-record character.\n");
                        printf1(ctx, "Exceeded max record length: %d characters.\n",nmax);
                        return(-1);
                    }
                }
            }
            get_rdat(ctx, fn);
            if (ctx->Out_Buf_Cnt <= 0) {
                if (j > 0) {
                    printf1(ctx, "\nError in retrieving file: %s\n",ctx->ZAFNam[fn]); 
                    printf1(ctx, "Can't find an end-of-record character.\n");
                    printf1(ctx, "Remaining characters: %d\nFound:",j);
                    q = buf;
                    for (i = 0; i < j; ++i)  
                        printf1(ctx, " %02x",*q++);
                    printf1(ctx, "\n\n");
                    if (j > 1 || *buf != 0x1a)
                        return(-1);     
                }
                *q = '\0';
                return(0);
            }
        }
    }
    else {
        if ((n = ctx->ZAFRLen[fn]) > nmax)
            n = nmax;
        while (j < n) {
            if (!ctx->Out_Buf_Cnt) {
                get_rdat(ctx, fn);
                if (!ctx->Out_Buf_Cnt)  
                    break;
            } 
            if ((i = ctx->Out_Buf_Cnt) > n - j)
                i = n - j; 
            memcpy(q,ctx->Out_Buf_Ptr,(size_t)(i));
            q += i;
            ctx->Out_Buf_Ptr += i;
            ctx->Out_Buf_Cnt -= i;
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

void get_rdat(TDAContext *ctx, int fn)
{
    if (!ctx->ZAEOF) {
        switch (ctx->ZAZTyp[fn]) {
            case 0:   lzs(ctx, fn);
                      break;
            case 1:   lzd(ctx, fn);
                      break; 
            case 2:   lzh_decode(ctx, fn);
                      break; 
            default:  printfe(ctx, "\nUnknown compression type %d\n",(int)ctx->ZAZTyp[fn]);
                      gerr_exit(ctx, 32);

        }
        ctx->Out_Buf_Ptr = ctx->Out_Buf_Adr;
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

int alloc_avar(TDAContext *ctx, int idx,int n,int opt)
{
    register int i,j,l;
    int err = -1;

    if (opt == 0) {
        err = 0;
        goto AAVFin;
    }
    if (n <= 0)
        return(0);

    if (!(ctx->AVIdx = (short *)calloc((size_t)(n),sizeof(short))))  
        goto AAVFin; 
    ctx->AVIdxA = n;
    memrq(ctx, n,sizeof(short));

    if (!(ctx->AVOff = (short *)calloc((size_t)(n),sizeof(short))))  
        goto AAVFin; 
    ctx->AVOffA = n;
    memrq(ctx, n,sizeof(short));

    if (!(ctx->AVLen = (short *)calloc((size_t)(n),sizeof(short))))  
        goto AAVFin; 
    ctx->AVLenA = n;
    memrq(ctx, n,sizeof(short));

    if (!(ctx->AVFmt1 = (short *)calloc((size_t)(n),sizeof(short))))  
        goto AAVFin; 
    ctx->AVFmt1A = n;
    memrq(ctx, n,sizeof(short));

    if (!(ctx->AVFmt2 = (short *)calloc((size_t)(n),sizeof(short))))  
        goto AAVFin; 
    ctx->AVFmt2A = n;
    memrq(ctx, n,sizeof(short));

    if (!(ctx->AVMBlnk = (int *)calloc((size_t)(n),sizeof(int))))  
        goto AAVFin; 
    ctx->AVMBlnkA = n;
    memrq(ctx, n,sizeof(int));

    if (!(ctx->AVMStar = (int *)calloc((size_t)(n),sizeof(int))))  
        goto AAVFin; 
    ctx->AVMStarA = n;
    memrq(ctx, n,sizeof(int));

    if (!(ctx->AVMPnt = (int *)calloc((size_t)(n),sizeof(int))))  
        goto AAVFin; 
    ctx->AVMPntA = n;
    memrq(ctx, n,sizeof(int));

    if (!(ctx->AVMGen = (int *)calloc((size_t)(n),sizeof(int))))  
        goto AAVFin; 
    ctx->AVMGenA = n;
    memrq(ctx, n,sizeof(int));

    i = 0;
    j = idx;
    while (j >= 0) {
        if (ctx->VTypA[j] == 2) {
            ctx->AVIdx[i] = (short)(j);
            l = (int)(strlen(ctx->VDef[j]) - 2);
            if (l > 0)
                ctx->AVFmt1[i] = (short)(-l);
            else  
                gerr_exit(ctx, 40);
            
            if (++i >= n)
                break;
        }
        j = ctx->VNxt[j];
    }
    if (i != n)
        gerr_exit(ctx, 46);
    return(0);

AAVFin:
    if (ctx->AVMGenA > 0) {
        free((char *)ctx->AVMGen);
        memrq(ctx, -ctx->AVMGenA,sizeof(int));
        ctx->AVMGenA = 0;
    }
    if (ctx->AVMPntA > 0) {
        free((char *)ctx->AVMPnt);
        memrq(ctx, -ctx->AVMPntA,sizeof(int));
        ctx->AVMPntA = 0;
    }
    if (ctx->AVMStarA > 0) {
        free((char *)ctx->AVMStar);
        memrq(ctx, -ctx->AVMStarA,sizeof(int));
        ctx->AVMStarA = 0;
    }
    if (ctx->AVMBlnkA > 0) {
        free((char *)ctx->AVMBlnk);
        memrq(ctx, -ctx->AVMBlnkA,sizeof(int));
        ctx->AVMBlnkA = 0;
    }
    if (ctx->AVFmt1A > 0) {
        free((char *)ctx->AVFmt1);
        memrq(ctx, -ctx->AVFmt1A,sizeof(short));
        ctx->AVFmt1A = 0;
    }
    if (ctx->AVFmt2A > 0) {
        free((char *)ctx->AVFmt2);
        memrq(ctx, -ctx->AVFmt2A,sizeof(short));
        ctx->AVFmt2A = 0;
    }
    if (ctx->AVLenA > 0) {
        free((char *)ctx->AVLen);
        memrq(ctx, -ctx->AVLenA,sizeof(short));
        ctx->AVLenA = 0;
    }
    if (ctx->AVOffA > 0) {
        free((char *)ctx->AVOff);
        memrq(ctx, -ctx->AVOffA,sizeof(short));
        ctx->AVOffA = 0;
    }
    if (ctx->AVIdxA > 0) {
        free((char *)ctx->AVIdx);
        memrq(ctx, -ctx->AVIdxA,sizeof(short));
        ctx->AVIdxA = 0;
    }
    if (err)  
        p_err(ctx, -2,1);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  check_avar(n,adic)  Search variable description file for the variables  */
/*                      in the list AVIdx, n variables. AVFmt1 contains     */
/*                      -length of variable name.                           */
/*                      Set: AVOff, AVLen, AVFmt1, AVFmt2.                  */
/*                      If adic != 0 recognize labels.                      */
/*                      Return index of data file if OK, -1 if error.       */

int check_avar(TDAContext *ctx, int n,int adic)
{
    register int i,j,l,ii;
    register char *p,*q;
    int w1,w2,m,fnni,fnn,fn,off,ll,len,err1,err2,err3;

    if (n <= 0)
        return(0);

    ll = 8;
    for (i = 0; i < n; ++i) {
        l = (int)(strlen(ctx->VDef[ctx->AVIdx[i]]) - 2);
        if (ll < l)
            ll = l;
    }
    printf1(ctx, "\nVariable  ");
    prnchar(ctx, ' ',ll - 8,0);
    printf1(ctx, "File  ");
    prnchar(ctx, ' ',ctx->FDefLen - 4,0);
    printf1(ctx, "FN  Off  Len  Format  Label\n");
    prnchar(ctx, '-',32 + ll + ctx->FDefLen,1);

    dbf_init(ctx, 1);  /* init buffer for reading the var description file */

    /*  read the variable description file */

    err1 = err2 = err3 = m = 0;
    fnni = fnn = -1;

    while (get_record(ctx, ctx->VFN,ctx->ZABuf,ctx->ZABLen) > 0) {

        if ((p = check_comment(ctx, ctx->ZABuf)) != NULL) { /* skip comment lines */
NXT1:
            if (*ctx->ZABuf == ' ')
                goto NXT;

            for (i = 0; i < n; ++i) {

                if (ctx->AVFmt1[i] < 0.0) {

                    l = (int)(-ctx->AVFmt1[i]);
                    j = ctx->AVIdx[i];
                    q = ctx->VDef[j] + 2;

                    if (!strncmp(p,q,(size_t)(l)) && *(p + l) == ' ') {

                        printf1(ctx, "%s  ",q);
                        prnchar(ctx, ' ',ll - l,0);

                        q = p + strlen(p);
                        while (q > p && (*--q == '\n' || *q == LF || *q == CR || *q == ' ')) ;
                        *++q = '\0';

                        w1 = w2 = 0;
                        p = skip_cb(ctx, p);
                        if (sscanf(p,"%d %d %d.%d",&fn,&off,&w1,&w2) == 4) ; 
                        else if (sscanf(p,"%d %d %d ",&fn,&off,&w1) == 3)   
                            w2 = 0;
                        else {
                            printf1(ctx, " cannot get information about variable\n");
                            fn = -1;
                        }

                        if (fnn >= 0 && fn == fnn)
                            j = fnni;
                        else {
                            j = -1;
                            for (ii = 0; ii < ctx->ZANF; ++ii) {
                                if (fn == (int)ctx->ZAFNum[ii]) {
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
                            printf1(ctx, "???? ");
                            prnchar(ctx, ' ',ctx->FDefLen - 4,0);
                            err2 = 1;
                        }
                        else {
                            printf1(ctx, "%s ",ctx->ZAFNam[j]);
                            prnchar(ctx, ' ',ctx->FDefLen - (int)strlen(ctx->ZAFNam[j]),0);
                            if (fn != fnn)
                                err3 = 1;
                        }
                        ctx->AVLen[i] = (short)(len = iabs(ctx, w1));

                        if (fn < 0 || off < 0 || len < 1 || w2 < 0 || w2 >= len)
                            err1 = 1;

                        if (w1 < 0) {       /* string variable */
                            ctx->VTyp[ctx->AVIdx[i]] = 1;
                            ctx->VSLen[ctx->AVIdx[i]] = (short)(w1);
                            w1 = w2 = 0;
                            if (len > ctx->SVBufLen) {
                                if (svb_alloc(ctx, len))
                                    return(-1);
                            }
                        }   
                        else {
                            ctx->VTyp[ctx->AVIdx[i]] = 3;
                            ctx->VSLen[ctx->AVIdx[i]] = (short)(get_afmt(ctx, &w1,&w2));
                        }
                        printf1(ctx, "%3d %5d %4d %4d.%-3d ",fn,off,len,w1,w2);
                        ctx->AVOff[i] = (short)(off);
                        ctx->AVFmt1[i] = (short)(w1); 
                        ctx->AVFmt2[i] = (short)(w2); 

                        p = skip_cb(ctx, p);
                        p = skip_cb(ctx, p);
                        p = skip_cb(ctx, p);
                        if (*p)
                            printf1(ctx, "%s\n",p);
                        else
                            printf1(ctx, "\n");
                        m++;

                        if (adic) {
                            while (get_record(ctx, ctx->VFN,ctx->ZABuf,ctx->ZABLen) > 0) {
                                if ((p = check_comment(ctx, ctx->ZABuf)) != NULL) {
                                    if (*ctx->ZABuf == ' ') {

                                        q = p + strlen(p);
                                        while (q > p && (*--q == '\n' || *q == LF || *q == CR || *q == ' ')) ;
                                        *++q = '\0';
                                        printf1(ctx, " %s\n",p);

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
            if (ctx->AVFmt1[i] < 0.0) {
                j = ctx->AVIdx[i];
                q = ctx->VDef[j] + 2;

                printf1(ctx, "%s  ",q);
                prnchar(ctx, ' ',ll + ctx->AVFmt1[i],0);
                printf1(ctx, "can't find this variable.\n");
            }
        }
        return(-1);
    }
    if (err1) {
        printf1(ctx, "\nIncorrect entries for at least one variable.\n");
        return(-1);
    }
    else if (err2 || fnni < 0) {
        printf1(ctx, "\nError: there are unkown data files.\n");
        return(-1);
    }
    else if (err3) {
        printf1(ctx, "\nError: variables do not belong to the same data file.\n");
        return(-1);
    }
    printf1(ctx, "\nUsing archive data file: %s\n",ctx->ZAFNam[fnni]);
    printf1(ctx, "Number of records: %d. Record length: %d",
                                         ctx->ZAFNRec[fnni],(int)ctx->ZAFRLen[fnni]);
    if ((int)ctx->ZAFRLen[fnni] == 0)  
        printf1(ctx, " (variable)");
    printf1(ctx, ".\n");
    return(fnni);
}

/* ------------------------------------------------------------------------ */
/*  get_avar()      Read next record from archive file AVDFN and put values */
/*                  of requested NVArc variables into AVVAL[].              */
/*                  Also count missings in AVMBlnk, and so on.              */
/*                                                                          */
/*                  Return -1 if error, 0 if EOF, 1 if OK.                  */

int get_avar(TDAContext *ctx)
{
    register int i,len;
    int mval,rlen;
    double tmp;
    register char *p,*q,*s;

    rlen = get_drec(ctx, ctx->AVDFN,ctx->ZABuf,ctx->ZABLen);
    if (rlen <= 0)
        return(rlen);

    for (i = 0; i < ctx->NVArc; ++i) {

        if (ctx->VTyp[ctx->AVIdx[i]] == 1)    /* string variable */
            tmp = 0.0;
        else {
            q = ctx->ZABuf;     /* pointer to buffer with current record */

            len = ctx->AVLen[i];
            p = q + ctx->AVOff[i];
            s = p + len;

            if (s > q + rlen) {
                tmp = ctx->MBlnkVal;
                ctx->AVMBlnk[i] += 1;
            }
            else {
                while (p < s && *p == ' ') {    /* skip leading blanks */
                    p++;
                    len--;
                }
                if (len <= 0) {
                    tmp = ctx->MBlnkVal;
                    ctx->AVMBlnk[i] += 1;
                }
                else {
                    tmp = dscan(ctx, p,len,&mval);
                    if (mval) {
                        if (mval == 1)
                            ctx->AVMBlnk[i] += 1;
                        else if (mval == 2)
                            ctx->AVMStar[i] += 1;
                        else if (mval == 3)
                            ctx->AVMPnt[i] += 1;
                        else 
                            ctx->AVMGen[i] += 1;
                    }
                }
            }
        }
        ctx->AVVAL[ctx->AVIdx[i]] = tmp;
    }
    return(1);
}

/*--------------------------------------------------------------------------*/
/*  get_astr(buf,i)                                                         */
/*                                                                          */
/*  Get string from archive variable i into buf.                            */
          
void get_astr(TDAContext *ctx, char *buf,int i)
{
    register int j;
    register char *p,*q,*s;
    int len;

    p = ctx->ZABuf + ctx->AVOff[i];     /* pointer to buffer with current record */
    q = buf;
    s = ctx->ZABuf + ctx->ZABLen;
    len = ctx->AVLen[i];

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

int arcc(TDAContext *ctx)
{
    register int i,j,n;
    int err,err1,terr,nrec,vdfi,vdfn,iv = 0;
    register char *p,*q; 

    err = vdfn = vdfi = -1;

    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->ARCDef == 0) {
        p_err(ctx, -9,1);
        return(-1);
    }
    printf1(ctx, "Archive check. Reading files defined in: %s\n\n",ctx->ZADNam);

    terr = 0;
    for (i = 0; i < ctx->ZANF; ++i) {      /* read all files in the archive */

        printf1(ctx, "File:%4d  %s",ctx->ZAFNum[i],ctx->ZAFNam[i]);
        prnchar(ctx, ' ',ctx->FDefLen - (int)strlen(ctx->ZAFNam[i]),0);
        /* tda_out_flush(); */
  
        if (ctx->ZAFTyp[i] == 2 && vdfi < 0) {
            vdfi = i;
            vdfn = ctx->ZAFNRec[i];
            iv = 0;
            if (vdfn > 0) {
                if (alloc_acc(ctx, (VNLMax + 1) * vdfn))
                    goto ACCFin;
            }
        }
        dbf_init(ctx, 1);
        /*****
        ZAFReq[i] = 1;
        ****/
        err1 = nrec = 0;
        while ((n = get_record(ctx, i,ctx->ZABuf,ctx->ZABLen)) > 0) {
            if (ctx->ZAFRLen[i]) {
                if (n != ctx->ZAFRLen[i]) {
                    printf1(ctx, " error in record length or number of records\n");
                    err1 = 1;
                    break;
                }
            }
            if (i == vdfi) {
                if ((p = check_comment(ctx, ctx->ZABuf)) != NULL) {
                    if (*ctx->ZABuf != ' ') {
                        q = ctx->AcC + iv * (VNLMax + 1);
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
            printf1(ctx, " records: %7d\n",nrec);
            if (nrec != ctx->ZAFNRec[i]) {
                printf1(ctx, "Inconsistent with defined number of records (%d)\n",
                                                                    ctx->ZAFNRec[i]);
                terr++;
            }
        }
        else
            terr++;
    }
    if (terr) {
        printf1(ctx, "Found %d error(s) in archive %s\n",terr,ctx->ZOONam);
        goto ACCFin;
    }

    /* check variable description file */

    if (vdfi < 0) {
        printf1(ctx, "Error: no variable description file.\n");
        goto ACCFin;
    }
    printf1(ctx, "Found %d variables.\n",iv);

    tda_qsort_r(ctx->AcC,(size_t)(iv),VNLMax + 1, iv_comp, ctx);

    for (i = 1; i < iv; ++i) {
        if (!strcmp(ctx->AcC + i * (VNLMax + 1),ctx->AcC + (i - 1) * (VNLMax + 1))) {
            printf1(ctx, "Error: at least one variable name is not unique.\n");
            goto ACCFin;
        }
    }
    printf1(ctx, "No errors found.\n");
    err = 0;
ACCFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  iv_comp()                                                               */

int iv_comp(const void *arg1, const void *arg2, void *_ctx)
{
    (void)_ctx;        /* unused: the signature is shared */
    return(strcmp((char *)arg1,(char *)arg2));
}

/* ------------------------------------------------------------------------ */
/*  arcv()      arcv(fn=...)=fname. Print archive variables.                */
/*              Return 0 or -1 if syntax error.                             */

int arcv(TDAContext *ctx)
{
    register int i,j,k;
    int err,fnd,n,nn,fn,fna,fnn,off,cflag,w1,w2,rec,vtyp,vslen;
    register char *p,*q;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->ARCDef == 0) {
        p_err(ctx, -9,1);
        return(-1);
    }
    if (parm(ctx, ctx->CmdBuf + 4,1,0))     /* get parameters */
        goto PAFin;   

    printf1(ctx, "Using variable description file: %s\n",ctx->ZAFNam[ctx->VFN]);

    err = n = 0;
    if (ctx->PMFN > 0) {

        if (alloc_acn(ctx, ctx->PMFN))
           goto PAFin;

        for (k = 0; k < ctx->PMFN; ++k) {
            j = -1;
            for (i = 0; i < ctx->ZANF; ++i) {
                if (!strcmp(ctx->PMFNam[k],ctx->ZAFNam[i]) && ctx->ZAFTyp[i] == 1) {
                    j = ctx->ZAFNum[i];
                    ctx->AcN[k] = j;
                    n++;
                    break;
                }
            }
            if (j < 0)  
                printf1(ctx, "Data file %s not in archive (ignored).\n",ctx->PMFNam[k]);
        }
        if (n == 0) {
            printf1(ctx, "No data files.\n");
            goto PAFin;
        }
    }

    dbf_init(ctx, 1);    /* init reading of var description file */

    if (ctx->PMFDef)
        fprintf(ctx->PMFd,"# arcd = %s;\n# nvar( \n",ctx->ZADNam);
    else
        printf1(ctx, "# arcd = %s;\n# nvar(\n",ctx->ZADNam);

    nn = 1;
    rec = fnd = 0;
    fnn = fna = -1;

    while (get_record(ctx, ctx->VFN,ctx->ZABuf,ctx->ZABLen) > 0) {

        rec++;
        if (*ctx->ZABuf == ' ')
            cflag = 1;
        else
            cflag = 0;

        if ((p = check_comment(ctx, ctx->ZABuf)) != NULL) {

            if (cflag && (ctx->PMArcDic == 0 || fnd == 0))
                continue; 

            q = p + strlen(p);
            while (q > p && (*--q == '\n' || *q == LF || *q == CR || *q == ' ')) ;
            *++q = '\0';

            if (fnd) {
                if (cflag && ctx->PMArcDic) {
                    if (ctx->PMFDef)
                        fprintf(ctx->PMFd,"#  %s\n",p);
                    else
                        printf1(ctx, "#  %s\n",p);
                    nn++;
                    continue; 
                }
                fnd = 0;
            }
            q = skip_cb(ctx, p);
                                 
            if (sscanf(q,"%d %d %d.%d",&fn,&off,&w1,&w2) == 4) ;
            else if (sscanf(q,"%d %d %d.%d",&fn,&off,&w1,&w2) == 3)  
                w2 = 0;
            else  
                w1 = 0;  

            if (fn < 0 || off < 0 || iabs(ctx, w1) < 1 || w2 < 0 || w2 >= iabs(ctx, w1)) {
                printf1(ctx, "Error in record %d of variable description file.\n",rec);
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
                for (i = 0; i < ctx->PMFN; ++i) {
                    if (ctx->AcN[i] == fn) {
                        j = 1;
                        break;
                    }
                }
            }
            if (j) {
                *(q - 1) = '\0';
                q = skip_cb(ctx, q);
                q = skip_cb(ctx, q);
                q = skip_cb(ctx, q);

                if (vtyp == 1)
                    j = -vslen;
                else
                    j = get_afmt(ctx, &w1,&w2);

                if (fnn != fn) {
                    for (k = 0; k < ctx->ZANF; ++k) {
                        if (ctx->ZAFNum[k] == fn) {
                            fnn = fn;
                            fna = k;
                            break;
                        }
                    }
                }
                if (fna < 0)
                    fna = 0;

                if (ctx->PMFDef)
                    fprintf(ctx->PMFd,"# %s<%d>[%d.%d] = A:%s, # [%s] %s\n",p,j,w1,w2,p,ctx->ZAFNam[fna],q);
                else
                    printf1(ctx, "# %s<%d>[%d.%d] = A:%s, # [%s] %s\n",p,j,w1,w2,p,ctx->ZAFNam[fna],q);
                nn++;
                fnd = 1;
            }
        }
    }
    if (ctx->PMFDef)
        fprintf(ctx->PMFd,"# );\n");
    else
        printf1(ctx, "# );\n");

    if (ctx->PMFDef && nn > 1)  
        printf1(ctx, "%d records written to: %s\n",nn,ctx->PMFdName);

PAFin:
    p_clean(ctx);
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

int arcvc(TDAContext *ctx)
{
    register int i,j,l;
    int n,nn,m,fnd,nv,err,rec,vptra,cnta,fnna,fn,off,w1,w2,flen;
    short *count,*fnn;
    char fname[17],*p,*q,**vptr;

    err = -1;
    nn = rec = fnna = cnta = vptra = 0;

    if (check_cmd(ctx, 1))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 5,2,1))     /* get parameters */
        goto PACFin;   

    printf1(ctx, "Checking variable description file: %s\n",ctx->PMFdName);

    if (ctx->PMOPT > 3)
        ctx->PMOPT = 3;

    if (alloc_acc(ctx, RLMaxDef))
        goto PACFin;

    nv = 0;     /* count number of variables */

    while (fgets(ctx->AcC,RLMaxDef,ctx->PMFd)) {
        if (*ctx->AcC != ' ' && (p = check_comment(ctx, ctx->AcC)) != NULL)  
            nv++;
    }

    if (!(vptr = (char **) calloc((size_t)(nv + 1),sizeof(char *)))) {
        p_err(ctx, -2,1);
        goto PACFin;
    }
    memrq(ctx, nv + 1,sizeof(char *));
    vptra = 1;

    if (!(count = (short *) calloc((size_t)(nv + 1),sizeof(short)))) {
        p_err(ctx, -2,1);
        goto PACFin;
    }
    memrq(ctx, nv + 1,sizeof(short));
    cnta = 1;

    if (!(fnn = (short *) calloc((size_t)(nv + 1),sizeof(short)))) {
        p_err(ctx, -2,1);
        goto PACFin;
    }
    memrq(ctx, nv + 1,sizeof(short));
    fnna = 1;

    fseek(ctx->PMFd,0,0);
    n = nn = rec = err = 0;
       
    while (fgets(ctx->AcC,RLMaxDef,ctx->PMFd)) {

        rec++;
        prn_message(ctx, rec,0,0);

        if (*ctx->AcC != ' ' && (p = check_comment(ctx, ctx->AcC)) != NULL) {

            p = skip_c(ctx, ctx->AcC);
            q = skip_b(ctx, p + 1);
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
                else if (iabs(ctx, w1) < 1)
                    err = 4;
                else if (w2 < 0 || w2 >= iabs(ctx, w1))
                    err = 5;
            }
            if (err)  
                goto PACFin;

            *p = '\0';
            fnd = 0;
            for (j = 0; j < nn; ++j) {
                if (!strcmp(vptr[j],ctx->AcC)) {
                    count[j] += 1;
                    fnd = 1;
                    break;
                }
            }
            if (!fnd) {
                l = (int)(strlen(ctx->AcC) + 1);
                if (!(vptr[nn] = (char *) calloc((size_t)(l),sizeof(char)))) {
                    p_err(ctx, -2,1);  
                    err = -1;
                    goto PACFin;
                }
                memrq(ctx, l,1);

                strcpy(vptr[nn],ctx->AcC);
                count[j] = 1;
                nn++;
            }
            fnn[n++] = (short)(fn);
        }
    }
    prn_message(ctx, rec,1,0);
    printf1(ctx, "Found %d variables.\n",nv);
    j = 0;
    for (i = 0; i < nn; ++i) {
        if (count[i] > 1) {
            printf1(ctx, "Found %2d times: %s\n",count[i],vptr[i]);
            j++;
        }
    }
    printf1(ctx, "Number of variable names used more than once: %d\n",j);

    if (j > 0 && ctx->PMF1Def) {         /* create new file */

        printf1(ctx, "Creating a new file with unique variable names.\n");

        fseek(ctx->PMFd,0,0);
        flen = n = rec = 0;
       
        while (fgets(ctx->AcC,RLMaxDef,ctx->PMFd)) {

            rec++;
            prn_message(ctx, rec,0,0);

            if (ctx->PMOPT > 1) {
                if (!strncmp(ctx->AcC,"# variable",10)) {
                    p = ctx->AcC + 38;
                    q = fname;
                    for (j = 0; j < 16; ++j) {
                        *q++ = *p++;
                        if (!*p || *p == '.' || *p == LF || *p == CR)
                            break;
                    }
                    *q = '\0';

                    if (ctx->PMOPT == 3) {
                        p = fname;
                        for (j = 0; j < 16; ++j) {
                            if (!*p)
                                break;
                            if (isalpha((int)*p))  
                                *p = (char)toupper((int)*p);
                            p++;
                        }
                    }
                    flen = (int)(strlen(fname));
                }
                if (flen == 0) {
                    printf1(ctx, "Error: cannot create new file names with opt = %d.\n",ctx->PMOPT);
                    err = -1;
                    goto PACFin;
                }
            }

            m = 0;
            if (*ctx->AcC != ' ' && (p = check_comment(ctx, ctx->AcC)) != NULL) {

                p = skip_c(ctx, ctx->AcC);
                *p = '\0';

                m = 0;
                for (j = 0; j < nn; ++j) {
                    if (!strcmp(vptr[j],ctx->AcC)) {
                        m = count[j] - 1;
                        break;
                    }
                }
                if (m > 0) {
                    if (ctx->PMOPT == 1)
                        fprintf(ctx->PMF1d,"%s_%d %s",ctx->AcC,fnn[n],p + 1);
                    else
                        fprintf(ctx->PMF1d,"%s_%s %s",ctx->AcC,fname,p + 1);
                }
                *p = ' ';
                n++;
            }
            if (m == 0)
                fprintf(ctx->PMF1d,"%s",ctx->AcC);
        }
        prn_message(ctx, rec,1,0);
        printf1(ctx, "%d records (%d variables) written to: %s\n",rec,n,ctx->PMF1dName);
    }
    err = 0;

PACFin:
    if (err > 0) {
        printf1(ctx, "Error in record %d: %s\n",rec,ctx->AcC);
        if (err == 1)
            printf1(ctx, "Need file number, offset and a correct format.\n");
        else if (err == 2)
            printf1(ctx, "Logical file number must be positive.\n");
        else if (err == 3)
            printf1(ctx, "Offset must not be negative.\n");
        else if (err == 4)
            printf1(ctx, "Width must be at least one column.\n");
        else if (err == 5)
            printf1(ctx, "Error in format specification.\n");
        err = -1;
    }
    for (i = 0; i < nn; ++i) {
        l = (int)(strlen(vptr[i]) + 1);
        free(vptr[i]);
        memrq(ctx, -l,1);
    }
    if (vptra) {  
        free((char *)vptr);
        memrq(ctx, -nv - 1,sizeof(char *));
    }
    if (cnta) {
        free((char *)count);
        memrq(ctx, -nv - 1,sizeof(short));
    }
    if (fnna) {
        free((char *)fnn);
        memrq(ctx, -nv - 1,sizeof(short));
    }
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_afmt(w1,w2)     create default print format and return storge size. */
/*                      not for string variables.                           */

int get_afmt(TDAContext *ctx, int *w1,int *w2)
{
    (void)ctx;        /* unused: the signature is shared */
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

