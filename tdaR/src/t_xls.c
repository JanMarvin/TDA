/****************************************************************************/
/*  t_xls                                                                   */
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
#include "t_var.h"
#include "t_gdat.h"
#include "t_parm.h"
#include "t_alloc.h"
#include "t_sort.h"
#include "t_eval.h"
#include "t_eval1.h"
#include "t_gf.h"
#include "t_spss.h"
#include "tda_context.h"
#include "t_xls.h"

/*  functions in t_xls.c */

int rxls(TDAContext *ctx); 
int xls_ole_extract(TDAContext *ctx, char *tmp);
int xls_read(TDAContext *ctx);
int xls_short(TDAContext *ctx, short *n);
int xls_int(TDAContext *ctx, int *n);
int xls_dbl(TDAContext *ctx, double *x);
int xls_rk(TDAContext *ctx, double *x);
int xls_unknown(TDAContext *ctx);
int xls_boundsheet(TDAContext *ctx);
int xls_sst(TDAContext *ctx);
int xls_labelsst(TDAContext *ctx);
int xls_blank(TDAContext *ctx);
int xls_mulblank(TDAContext *ctx);
int xls_mulrk(TDAContext *ctx);
int xls_formula(TDAContext *ctx);
int xls_number(TDAContext *ctx);
void xls_rowcol(TDAContext *ctx, int row,int col);
void xls_putx(TDAContext *ctx, int row,int col,double x);
int xls_check(TDAContext *ctx, int j);
int rucinet(TDAContext *ctx); 
int xls_cellok(TDAContext *ctx, int row, int col);


/*--------------------------------------------------------------------------*/
/*  Global variables                                                        */

#define XLSHLEN 64          /* used to skip over OLE header                 */
#define XLSNFMax 50         /* max number of files                          */



/* ------------------------------------------------------------------------ */
/*  rxls            Read a xls (excel) file.                                */
/*                                                                          */
/*                  rxls(                                                   */
/*                      df=...,     write data to output file               */
/*                      fmt=...,    print format, def. 10.4                 */
/*                      prot=...,   protocol file                           */
/*                      nc=...,     1 don't write comment lines, def. 0     */
/*                      ns=...,     1 only write numerical columns, def. 0  */
/*                      ni=...,     1 ignore unknown record types, def. 0   */
/*                  ) = fname;      name of xls file                        */
/*                                                                          */
/*                                                                          */
/*  NOTE: the command only supports the BIFF8 format. Interpretation is     */
/*  based on the document "Microsoft Excel File Format" (Febr. 2002)        */
/*  written by Daniel Rentz and published as an internet document by        */  
/*  OpenOffice.org (http://sc.openoffice.org/excelfileformat.pdf).          */
/*                                                                          */
/*  NOTE: this command is not an "official" part of TDA and will not be     */
/*  maintained.                                                             */
/*                                                                          */
/*  Operation:                                                              */
/*                                                                          */  
/*  The command first reads the xls file and records the size of the        */
/*  table, the contents of the cells and the workbook globals file.         */
/*                                                                          */
/*  Then the workbook globals file is read again and a list of all strings  */
/*  (labels) is saved in internal memory.                                   */
/*                                                                          */
/*  Finally, each worksheet file is read again and the data are written     */
/*  to the output file specified with the df parameter. The data for        */
/*  different worksheets are written one after another into the output      */
/*  file. If there is more than one worksheet, there will be an additional  */
/*  first column that records the worksheet number. This column can be      */
/*  used to select data for specific worksheets.                            */
/*                                                                          */
/*  a)  If a row of the table contains at least one numerical entry (cell)  */
/*      it is written as a data record into the output file. Otherwise      */
/*      the record will begin with a comment sign (#) followed by all       */
/*      that are contained in cells of the current row.                     */
/*                                                                          */  
/*      If nc=1, the comment lines will not be written to the output file.  */
/*                                                                          */  
/*  b)  Columns that contain at least one string cell are taken to be       */
/*      string variables. If, in the same column, strings and numerical     */
/*      entries occur both (unfortunately, this is often the case with      */
/*      xls files), the numerical entries are written as strings with       */
/*      the free format "%-g".                                              */
/*                                                                          */
/*  c)  The standard output will contain information about the columns      */
/*      of the table, there type and additional labels, if present.         */
/*                                                                          */
/*  d)  If the prot parameter is used, a list of all strings from the       */
/*      workbook globals section is written into the protocol file.         */
/*                                                                          */  
/*  Return 0 if OK, -1 if error.                                            */

/* ------------------------------------------------------------------------ */
/*  xls_ole_extract     An .xls file is an OLE2 compound document, and the  */
/*  BIFF8 stream ("Workbook") inside it need not occupy consecutive         */
/*  sectors: a writer may interleave the summary streams or the FAT with    */
/*  it, and the sequential read that follows the BOF then runs into a       */
/*  foreign sector and reports a bad record length.  This follows the       */
/*  document's FAT (and the mini FAT for a stream below the cutoff), writes */
/*  the stream in order into the file named in tmp, and reopens it as       */
/*  ctx->PMFd; the parser then reads a plain BIFF stream from offset 0.     */
/*                                                                          */
/*  Layout, from the MS-CFB specification: 512-byte header with the sector  */
/*  shift at 30, the mini sector shift at 32, the directory start at 48,    */
/*  the mini stream cutoff at 56, the mini FAT start and count at 60/64,    */
/*  the DIFAT start and count at 68/72 and its first 109 entries at 76;     */
/*  sector n starts at (n + 1) * sector size; a directory entry is 128      */
/*  bytes with the UTF-16 name at 0, its byte length at 64, the type at 66  */
/*  (2 stream, 5 root), the start sector at 116 and the size at 120.        */
/*                                                                          */
/*  Return 0 if the stream was extracted, 1 if the file is not an OLE2      */
/*  document (left untouched, read as before), -1 on error.                 */

static unsigned int ole_u32(const unsigned char *b)
{
    return (unsigned int)b[0] | ((unsigned int)b[1] << 8) |
           ((unsigned int)b[2] << 16) | ((unsigned int)b[3] << 24);
}

static int ole_read_sector(FILE *f, unsigned int sec, unsigned int ssz,
                           unsigned char *buf)
{
    if (fseek(f,(long)((sec + 1U) * ssz),SEEK_SET))
        return(-1);
    return fread(buf,1,ssz,f) == ssz ? 0 : -1;
}

int xls_ole_extract(TDAContext *ctx, char *tmp)
{
    static const unsigned char magic[8] = {0xD0,0xCF,0x11,0xE0,0xA1,0xB1,0x1A,0xE1};
    unsigned char hdr[512];
    unsigned char *sec = NULL,*fat = NULL,*mfat = NULL,*mini = NULL;
    unsigned int ssz,msz,cutoff,dirstart,nfat,mfatstart,nmfat,difstart,ndif;
    unsigned int nfatent,nmfatent,i,k,s,wbstart,wbsize,rootstart,rootsize;
    unsigned int minisz,written,chunk;
    int found,err = -1;
    long nsec;
    FILE *out = NULL;

    rewind(ctx->PMFd);
    if (fread(hdr,1,512,ctx->PMFd) != 512 || memcmp(hdr,magic,8) != 0) {
        rewind(ctx->PMFd);
        return(1);
    }
    ssz = 1U << (hdr[30] | (hdr[31] << 8));
    msz = 1U << (hdr[32] | (hdr[33] << 8));
    nfat = ole_u32(hdr + 44);
    dirstart = ole_u32(hdr + 48);
    cutoff = ole_u32(hdr + 56);
    mfatstart = ole_u32(hdr + 60);
    nmfat = ole_u32(hdr + 64);
    difstart = ole_u32(hdr + 68);
    ndif = ole_u32(hdr + 72);
    if ((ssz != 512 && ssz != 4096) || msz > ssz || nfat == 0 || nfat > 100000U) {
        printf1(ctx, "Error: cannot interpret the OLE2 header of the xls file.\n");
        return(-1);
    }
    fseek(ctx->PMFd,0L,SEEK_END);
    nsec = ftell(ctx->PMFd) / (long)ssz;
    if (nsec < 2)
        return(-1);

    sec = (unsigned char *)malloc(ssz);
    fat = (unsigned char *)malloc((size_t)nfat * ssz);
    if (!sec || !fat) {
        p_err(ctx, -2,1);
        goto OLEFin;
    }

    /* the FAT: its sector numbers come from the header's DIFAT and, past
       109 of them, from the DIFAT chain */
    k = 0;
    for (i = 0; i < 109 && k < nfat; ++i) {
        s = ole_u32(hdr + 76 + 4 * i);
        if (s >= 0xFFFFFFFAU)
            break;
        if (ole_read_sector(ctx->PMFd,s,ssz,fat + (size_t)k * ssz))
            goto OLEBad;
        k++;
    }
    s = difstart;
    for (i = 0; i < ndif && k < nfat && s < 0xFFFFFFFAU; ++i) {
        if (ole_read_sector(ctx->PMFd,s,ssz,sec))
            goto OLEBad;
        for (unsigned int j = 0; j + 1 < ssz / 4 && k < nfat; ++j) {
            unsigned int t = ole_u32(sec + 4 * j);
            if (t >= 0xFFFFFFFAU)
                break;
            if (ole_read_sector(ctx->PMFd,t,ssz,fat + (size_t)k * ssz))
                goto OLEBad;
            k++;
        }
        s = ole_u32(sec + ssz - 4);
    }
    nfatent = k * (ssz / 4);

    /* the directory: the root entry, and the Workbook (or Book) stream */
    found = 0;
    wbstart = wbsize = rootstart = rootsize = 0;
    for (s = dirstart, i = 0; s < 0xFFFFFFFAU && i < nfatent; ++i, s = ole_u32(fat + 4 * s)) {
        if (s >= nfatent || ole_read_sector(ctx->PMFd,s,ssz,sec))
            goto OLEBad;
        for (k = 0; k + 128 <= ssz; k += 128) {
            unsigned char *e = sec + k;
            unsigned int nlen = (unsigned int)(e[64] | (e[65] << 8));
            if (i == 0 && k == 0) {
                rootstart = ole_u32(e + 116);
                rootsize = ole_u32(e + 120);
            }
            if (e[66] != 2 || found)
                continue;
            if ((nlen == 18 && memcmp(e,"W\0o\0r\0k\0b\0o\0o\0k\0",16) == 0) ||
                (nlen == 10 && memcmp(e,"B\0o\0o\0k\0",8) == 0)) {
                wbstart = ole_u32(e + 116);
                wbsize = ole_u32(e + 120);
                found = 1;
            }
        }
    }
    if (!found) {
        printf1(ctx, "Error: no Workbook stream in the xls file.\n");
        goto OLEFin;
    }

    if (!(out = fopen(tmp,OPEN_WB))) {
        printf1(ctx, "Error: cannot write %s\n",tmp);
        goto OLEFin;
    }
    written = 0;
    if (wbsize >= cutoff) {
        for (s = wbstart; s < 0xFFFFFFFAU && written < wbsize; s = ole_u32(fat + 4 * s)) {
            if (s >= nfatent || ole_read_sector(ctx->PMFd,s,ssz,sec))
                goto OLEBad;
            chunk = wbsize - written < ssz ? wbsize - written : ssz;
            if (fwrite(sec,1,chunk,out) != chunk)
                goto OLEBad;
            written += chunk;
        }
    }
    else {
        /* a small stream lives in the mini stream, itself a chain of
           ordinary sectors held by the root entry, and is addressed
           through the mini FAT in mini sectors */
        mfat = (unsigned char *)malloc((size_t)(nmfat ? nmfat : 1) * ssz);
        mini = (unsigned char *)malloc((size_t)(rootsize ? rootsize : 1) + ssz);
        if (!mfat || !mini) {
            p_err(ctx, -2,1);
            goto OLEFin;
        }
        for (s = mfatstart, k = 0; s < 0xFFFFFFFAU && k < nmfat; ++k, s = ole_u32(fat + 4 * s)) {
            if (s >= nfatent || ole_read_sector(ctx->PMFd,s,ssz,mfat + (size_t)k * ssz))
                goto OLEBad;
        }
        nmfatent = k * (ssz / 4);
        for (s = rootstart, minisz = 0; s < 0xFFFFFFFAU && minisz < rootsize; s = ole_u32(fat + 4 * s)) {
            if (s >= nfatent || ole_read_sector(ctx->PMFd,s,ssz,mini + minisz))
                goto OLEBad;
            minisz += ssz;
        }
        for (s = wbstart; s < 0xFFFFFFFAU && written < wbsize; s = ole_u32(mfat + 4 * s)) {
            if (s >= nmfatent || (s + 1U) * msz > minisz)
                goto OLEBad;
            chunk = wbsize - written < msz ? wbsize - written : msz;
            if (fwrite(mini + (size_t)s * msz,1,chunk,out) != chunk)
                goto OLEBad;
            written += chunk;
        }
    }
    if (written != wbsize)
        goto OLEBad;
    fclose(out);
    out = NULL;
    fclose(ctx->PMFd);
    if (!(ctx->PMFd = fopen(tmp,OPEN_RB))) {
        printf1(ctx, "Error: cannot reopen %s\n",tmp);
        ctx->PMFDef = 0;
        goto OLEFin;
    }
    err = 0;
    goto OLEFin;

OLEBad:
    printf1(ctx, "Error: the xls file's sector allocation table is damaged.\n");
OLEFin:
    if (out)
        fclose(out);
    free(sec);
    free(fat);
    free(mfat);
    free(mini);
    return(err);
}

int rxls(TDAContext *ctx)  
{
    register int i,j,k,jj;
    char buf[XLSHLEN + 1];
    short sn,len,ver,typ;
    int err,ii,l,n,nc,ne,nrow,ncol,eof,nrec;
    long fptr;
    char tmp[64];
    int extracted;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 4,10,1)) {   /* get parameters */
        p_clean(ctx);
        return(-1);   
    }   
    printf1(ctx, "Reading xls file: %s\n\n",ctx->PMFdName);

    /* pull the BIFF stream out of the OLE2 container first, so a
       stream stored in non-consecutive sectors reads in order */
    snprintf(tmp,sizeof(tmp),"%s","tda_rxls.tmp");
    extracted = xls_ole_extract(ctx, tmp);
    if (extracted < 0) {
        p_clean(ctx);
        remove(tmp);
        return(-1);
    }
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    ctx->ProtNRec = 0;
    ctx->XLSLNN = ctx->XLSLNA = ctx->XLSLN = ctx->XLSNF = ctx->XLSRF = 0;
    ctx->XLSLMin = ctx->XLSRowMin = ctx->XLSColMin = ctx->INTMAX;
    ctx->XLSLMax = ctx->XLSRowMax = ctx->XLSColMax = -1;
           
    /* skip the OLE header until one finds begin of BIFF8 */

    ii = 0;
    fptr = 0L;
    while (fread(buf,sizeof(char),XLSHLEN,ctx->PMFd) == XLSHLEN) {                   
        if (*buf == 0x09 && *(buf + 1) == 0x08) {
            ii = 1;
            fseek(ctx->PMFd,fptr,0);
            break;
        }
        fptr = ftell(ctx->PMFd);
    }
    if (ii == 0) {
        printf1(ctx, "Error: cannot find BIFF8 begin of file (0x0809).\n");
        goto RXLSFin;
    }
    eof = 0;
    ctx->XLSV = 0;

    while (1) {

        /* the extracted stream ends exactly after the last worksheet's
           EOF record: that is the end, not a read error */
        if (eof) {
            int c = fgetc(ctx->PMFd);
            if (c == EOF)
                break;
            ungetc(c,ctx->PMFd);
        }
        if (xls_short(ctx, &sn)) {           /* should be BOF */
            if (eof)
                break;
            goto RXLSFin;
        }
        if (sn != 0x0809) {
            if (eof == 0)
                printf1(ctx, "Error: cannot find BIFF8 begin of file (0x0809). Found: %04x\n",sn);
            break;             
        }
        if (xls_short(ctx, &len))     /* length */
            goto RXLSFin;

        if (xls_short(ctx, &ver))     /* version */
            goto RXLSFin;

        if (ctx->XLSV == 0) {
            printf1(ctx, "Version: %04x ",ver);
            if (ver == 0x0600 && len == 16) {
                ctx->XLSV = 8;
                printf1(ctx, "BIFF8\n");
            }
            else if (len == 8) {
                ctx->XLSV = 5;
                newline(ctx);
            }
            else {
                printf1(ctx, "[cannot use this version]\n");
                goto RXLSFin;
            }
        }
        if (xls_short(ctx, &typ))     /* type */
            goto RXLSFin;

        len -= 4;
        if (len < 4 || len > 12) {
            printf1(ctx, "Error: cannot interpret BOF.\n");
            goto RXLSFin;
        }
        if (fread(buf,sizeof(char),(size_t)(len),ctx->PMFd) != (size_t)(len)) {    /* skip rest of BOF */
            p_err(ctx, -7,1);
            goto RXLSFin;
        }
        printf1(ctx, "Begin of new file: ");
        fptr = ftell(ctx->PMFd);
        eof = 0;
        switch (typ) {
            case 0x0005:    if (ctx->XLSNF != 0) {
                                printf1(ctx, "Error: workbook globals should be the first file.\n");
                                goto RXLSFin;
                            }
                            ctx->XLSFPtr[ctx->XLSNF] = fptr;
                            ctx->XLSNF++;

                            printf1(ctx, "workbook globals.\n");
                            if (xls_read(ctx))
                                goto RXLSFin;
                            eof = 1;
                            break;

            case 0x0006:    printf1(ctx, "visual basic module.\n");
                            break;

            case 0x0010:    if (ctx->XLSNF == 0) {
                                printf1(ctx, "Error: worksheet should follow a globals section.\n");
                                goto RXLSFin;
                            }
                            if (ctx->XLSNF >= XLSNFMax) {
                                printf1(ctx, "Error: exceeded max number of files (%d).\n",XLSNFMax);
                                goto RXLSFin;
                            }
                            ctx->XLSFPtr[ctx->XLSNF] = fptr;
                            ctx->XLSNF++;

                            printf1(ctx, "worksheet.\n");
                            if (xls_read(ctx))
                                goto RXLSFin;
                            eof = 1;
                            break;

            case 0x0020:    printf1(ctx, "chart.\n");
                            break;

            case 0x0040:    printf1(ctx, "BIFF4 macro sheet.\n");
                            break;

            case 0x0100:    printf1(ctx, "BIFF4 workbooks globals.\n");
                            break;
            default:        printf1(ctx, "unknown.\n");
        }                   
        if (eof == 0) {
            printf1(ctx, "Not supported.\n");
            goto RXLSFin;
        }
    }
    printf1(ctx, "Number of worksheets: %d\n\n",ctx->XLSNF - 1);
    if (ctx->XLSNF < 2)
        goto RXLSFin;

    /* save strings from SST */

    if (!(ctx->XLSLabel = (char **)calloc((size_t)(ctx->XLSLN),sizeof(char *)))) {
        p_err(ctx, -2,1);
        ctx->XLSLN = 0;             
        goto RXLSFin;
    }
    ctx->XLSLNN = ctx->XLSLN;
    memrq(ctx, ctx->XLSLNN,sizeof(char *));

    ctx->XLSRF = 1;
    fseek(ctx->PMFd,ctx->XLSFPtr[0],0);
    if (xls_read(ctx))
        goto RXLSFin;

    /* get data for all worksheets */

    ctx->XLSNSkip = 0;
    ctx->XLSRow = ctx->XLSRowMax + 1;
    ctx->XLSCol = ctx->XLSColMax + 1;
    nrow = ctx->XLSRowMax - ctx->XLSRowMin + 1;
    ncol = ctx->XLSColMax - ctx->XLSColMin + 1;

    printf1(ctx, "Maximal size of table: %d (%d - %d) rows, %d (%d - %d) columns.\n",
        nrow,ctx->XLSRowMin,ctx->XLSRowMax,ncol,ctx->XLSColMin,ctx->XLSColMax);

    if (ctx->XLSRowMin > ctx->XLSRowMax || ctx->XLSColMin > ctx->XLSColMax)
        goto RXLSFin;

    if (alloc_acx(ctx, ctx->XLSRow * ctx->XLSCol + 1))
        goto RXLSFin;
    if (alloc_acn(ctx, ctx->XLSRow * ctx->XLSCol + 1))     /* >0 = number of label */
        goto RXLSFin;                       /* -1 not used */
                                            /* -2 numerical */
                                            /* -3 blank */
    if (alloc_aci(ctx, ctx->XLSRow + 1))
        goto RXLSFin;
    if (alloc_acj(ctx, ctx->XLSCol + 1))
        goto RXLSFin;
    if (alloc_ack(ctx, ctx->XLSCol + 1))
        goto RXLSFin;

    ctx->XLSRF = 2;
    nrec = 0;
    for (ii = 1; ii < ctx->XLSNF; ++ii) {

        printf1(ctx, "\nReading again: worksheet %d\n",ii);

        for (i = 0; i <= ctx->XLSRow * ctx->XLSCol; ++i) {
            ctx->AcN[i] = -1;
            ctx->AcX[i] = 0.0;
        }
        for (i = 0; i <= ctx->XLSRow; ++i)  
            ctx->AcI[i] =  0;
        for (j = 0; j <= ctx->XLSCol; ++j)  
            ctx->AcK[j] = ctx->AcJ[j] =  0;

        fseek(ctx->PMFd,ctx->XLSFPtr[ii],0);
        if (xls_read(ctx))             /* read worksheet ii */
            goto RXLSFin;

        /***
        for (i = XLSRowMin; i <= XLSRowMax; ++i) {
            tda_out("%5d : ",AcI[i]);
            for (j = XLSColMin; j <= XLSColMax; ++j)  
                tda_out("%3d ",AcN[i * XLSCol + j]);
            newline(ctx);
        }
        ***/

        ne = 0;
        for (i = ctx->XLSRowMin; i <= ctx->XLSRowMax; ++i) {
            if (ctx->AcI[i] > 0)
                ne++;
        }
        printf1(ctx, "%d rows contain at least one numerical entry and will be written\n",ne);
        printf1(ctx, "as data records to the output file. The following columns will be used.\n\n");
         
        printf1(ctx, "Column   New  Type         Row  Label\n");
        prnchar(ctx, '-',50,1);                       

        nc = 0;
        for (j = ctx->XLSColMin; j <= ctx->XLSColMax; ++j) {

            n = xls_check(ctx, j);
            if (n == 1) {
                printf1(ctx, "%6d %5d  numerical",j,nc);
                ctx->AcJ[nc++] = j;    
            }
            else if (n == 2) {
                if (ctx->PMNS != 1) {
                    printf1(ctx, "%6d %5d  string   ",j,nc);
                    ctx->AcJ[nc++] = j;
                    ctx->AcK[j] = 1;                   
                }
                else
                    printf1(ctx, "%6d        string   ",j);
            }
            else        
                printf1(ctx, "%6d        not used ",j);

            n = 0;
            for (i = ctx->XLSRowMin; i <= ctx->XLSRowMax; ++i) {
                jj = i * ctx->XLSCol + j;
                k = ctx->AcN[jj];        
                if (k >= 0 && k < ctx->XLSLNA) {
                    if (ctx->AcI[i] <= 0) {
                        if (n)
                            prnchar(ctx, ' ',23,0);
                        printf1(ctx, "  %5d  %s\n",i,ctx->XLSLabel[k]);
                        n++;
                    }
                }   
                if (ctx->AcI[i] > 0 && ctx->AcK[j]) {           /* string */ 
                    if (k == -2) {
                        snprintf(buf,sizeof(buf),"%-g",ctx->AcX[jj]);
                        ctx->AcK[j] = (int)(imax(ctx, ctx->AcK[j],(int)(strlen(buf))));
                    }
                    else if (k >= 0 && k < ctx->XLSLNA)
                        ctx->AcK[j] = (int)(imax(ctx, ctx->AcK[j],(int)(strlen(ctx->XLSLabel[k]))));
                }
            }
            if (n == 0)
                newline(ctx);
        }
        if (ne == 0)
            continue;

        if (ctx->PMF1Def) {  

            for (i = ctx->XLSRowMin; i <= ctx->XLSRowMax; ++i) {
                if (ctx->AcI[i] > 0) {
                    if (ctx->XLSNF > 2)
                        fprintf(ctx->PMF1d,"%4d ",ii);

                    for (jj = 0; jj < nc; ++jj) {
                        j = ctx->AcJ[jj];
                        k = i * ctx->XLSCol + j;
                        if (ctx->AcK[j] == 0) {
                            double xv = (ctx->AcN[k] == -2) ? ctx->AcX[k] : -1.0;
                            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,xv);
#ifdef TDA_R_PACKAGE
                            tda_export_str_row(ctx, "rxls.strings", "");
#endif
#ifdef TDA_R_PACKAGE
                            /*  Numeric cells as doubles.  Without this
                                the only route back was the df= file, so
                                tda_read_xls PARSED TEXT -- the one thing
                                this whole channel exists to remove.  A
                                text cell stays out: it is written as the
                                string itself.  */
                            tda_export_cell(ctx, "rxls.values", xv);
#endif
                        }
                        else {                        /* string */
#ifdef TDA_R_PACKAGE
                            /*  Text cells, one per cell in row-major
                                order so the block reshapes against
                                rxls.values; a numeric cell contributes
                                an empty string.  Without this a text
                                column came back as all-NA once the
                                frame was built from the exports --
                                clippy.xls lost "paperclip" entirely.  */
                            tda_export_str_row(ctx, "rxls.strings",
                                (ctx->AcN[k] >= 0 && ctx->AcN[k] < ctx->XLSLNA)
                                ? ctx->XLSLabel[ctx->AcN[k]] : "");
                            /*  and a NaN in the numeric block, so the
                                two blocks are the SAME SHAPE, cell for
                                cell.  Exporting numerics only for
                                numeric columns left the two with
                                different widths and nothing to line
                                them up by.  */
                            tda_export_cell(ctx, "rxls.values", (double)NAN);
#endif
                            if (ctx->AcN[k] >= 0 && ctx->AcN[k] < ctx->XLSLNA) {
                                fprintf(ctx->PMF1d,"%s ",ctx->XLSLabel[ctx->AcN[k]]);
                                l = (int)(strlen(ctx->XLSLabel[ctx->AcN[k]]));
                                fprnchar(ctx, ctx->PMF1d,' ',imax(ctx, 0,ctx->AcK[j] - l),0);
                            }
                            else if (ctx->AcN[k] == -2) {
                                snprintf(buf,sizeof(buf),"%-g",ctx->AcX[k]);
                                fprintf(ctx->PMF1d,"%s ",buf);
                                l = (int)(strlen(buf));
                                fprnchar(ctx, ctx->PMF1d,' ',imax(ctx, 0,ctx->AcK[j] - l),0);
                            }
                            else
                                fprnchar(ctx, ctx->PMF1d,' ',ctx->AcK[j] + 1,0);
                        }
                    }
                    fprintf(ctx->PMF1d,"\n");
#ifdef TDA_R_PACKAGE
                    tda_export_endrow(ctx, "rxls.values");
#endif
                    nrec++;
                }
                else if (ctx->PMNC != 1) {
                    n = 1;
                    for (j = ctx->XLSColMin; j <= ctx->XLSColMax; ++j) {
                        jj = i * ctx->XLSCol + j;
                        k = ctx->AcN[jj];
                        if (k >= 0 && k < ctx->XLSLNA) {
                            if (n) {
                                fprintf(ctx->PMF1d,"# ");
                                n = 0;
                            }
                            fprintf(ctx->PMF1d,"[%d] %s ",j,ctx->XLSLabel[k]);
                        }
                    }
                    if (n == 0) {
                        fprintf(ctx->PMF1d,"\n");
                        nrec++;
                    }
                }
            }
        }
#ifdef TDA_R_PACKAGE
        /*  One block per SHEET.  Sheets need not have the same
            width -- readxl's datasets.xls holds 11, 1 and 5 column
            ones -- so a single accumulator cannot hold them and only
            the first survived.  */
        tda_export_flush(ctx, "rxls.values");
        tda_export_str_flush(ctx, "rxls.strings");
#endif
    }
    /*  Braced deliberately.  The skipped-cell count was inserted between
        `if (ctx->PMF1Def)` and its body, which silently made it
        conditional on an output file having been requested AND made the
        records-written line unconditional -- it printed even when no
        file was written.  gcc's -Wmisleading-indentation catches this;
        the makefile's CFLAGS do not carry -Wall, so the build here did
        not.  */
    if (ctx->XLSNSkip) {
        printf1(ctx, "Cells outside the measured table (skipped): %d\n",
                ctx->XLSNSkip);
    }
    if (ctx->PMF1Def) {
        printf1(ctx, "\n%d records written to: %s\n",nrec,ctx->PMF1dName);
    }
    if (ctx->PMProtFd)
        printf1(ctx, "%d records written to: %s\n",ctx->ProtNRec,ctx->PMProtFName);
    err = 0;

RXLSFin:
    if (ctx->XLSLNN > 0) {
        for (i = 0; i < ctx->XLSLNA; ++i) {
            n = (int)(strlen(ctx->XLSLabel[i]) + 1);
            free(ctx->XLSLabel[i]);
            memrq(ctx, -n,sizeof(char));
        }
        free((char *)ctx->XLSLabel);
        memrq(ctx, -ctx->XLSLNN,sizeof(char *));
    }
    p_clean(ctx);
    if (extracted == 0)
        remove(tmp);
    return(err);
}


/* ------------------------------------------------------------------------ */
/*  xls_read            Read one file.                                      */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_read(TDAContext *ctx)
{
    short sn;
    long fptr;

    while (1) {

        if (xls_short(ctx, &sn))
            return(-1);     

        switch (sn) {
            case 0x0006:    if (xls_formula(ctx))          /* formula, value */
                                return(-1);    
                            break;

            case 0x0085:    if (xls_boundsheet(ctx))
                                return(-1);    
                            break;

            case 0x00fc:    if (xls_sst(ctx))
                                return(-1);     
                            break;

            case 0x000a:    if (xls_short(ctx, &sn))         /* EOF */
                                return(-1);     
                            return(0);

            case 0x00bd:    if (xls_mulrk(ctx))            /* multiple RK */
                                return(-1);     
                            break;       

            case 0x00be:    if (xls_mulblank(ctx))         /* multiple blanks */
                                return(-1);     
                            break;       

            case 0x00fd:    if (xls_labelsst(ctx))         /* label SST */
                                return(-1);     
                            break;       

            case 0x0201:    if (xls_blank(ctx))            /* blank cell */
                                return(-1);     
                            break;       

            case 0x0203:    if (xls_number(ctx))           /* number */
                                return(-1);     
                            break;       

            case 0x027e:    if (xls_rkvalue(ctx))          /* rk-value */
                                return(-1);     
                            break;       


            case 0x000c:    /* calccount */
            case 0x000d:    /* calcmode */
            case 0x000e:
            case 0x000f:
            case 0x0010:
            case 0x0011:
            case 0x0012:
            case 0x0013:    /* password */
            case 0x0014:                     
            case 0x0015:                     
            case 0x0017:    /* externsheet */
            case 0x0018:    /* name */
            case 0x0019:
            case 0x001a:
            case 0x001b:    /* ?-??? */
            case 0x001d:    /* ??? */
            case 0x0022:
            case 0x0026:
            case 0x0027:
            case 0x0028:
            case 0x0029:
            case 0x002a:
            case 0x002b:
            case 0x0031:    /* font */
            case 0x003d:
            case 0x0040:    /* backup */ 
            case 0x0041:    /* ??? */
            case 0x0042:
            case 0x004d:
            case 0x0055:    /* defcolwidth */
            case 0x005c:
            case 0x005d:    /* ??? */
            case 0x005f:
            case 0x007d:    /* colinfo */
            case 0x0080:                        
            case 0x0081:                        
            case 0x0082:                        
            case 0x0083:                        
            case 0x0084:                        
            case 0x008c:                        
            case 0x008d:
            case 0x0090:    /* ??? */
            case 0x0092:    /* palette */
            case 0x0094:    /* ??? */
            case 0x0099:    /* ??? */
            case 0x009c:
            case 0x00a0:    /* ??? */
            case 0x00a1:
            case 0x00ab:    /* ?-??? */
            case 0x00bf:    /* ?-??? */
            case 0x00c0:    /* ?-??? */
            case 0x00c1:
            case 0x00d7:
            case 0x00da:    /* bookbool */
            case 0x00e0:    /* XF */
            case 0x00e1:                
            case 0x00e2:                
            case 0x00e5:    /* merged cells */
            case 0x00eb:    /* ??? */   
            case 0x00ec:    /* ??? */   
            case 0x00ef:                 
            case 0x00ff:    /* extended SST */
            case 0x013d:
            case 0x0160:
            case 0x0161:
            case 0x01ae:    /* supblock */
            case 0x01af:
            case 0x01b7:
            case 0x01bc:
            case 0x01c0:
            case 0x01c1:
            case 0x0200:    /* dimensions */
            case 0x0204:    /* label version 5 */
            case 0x0208:
            case 0x020b:
            case 0x0225:
            case 0x023e:    /* ??? */
            case 0x0293:
            case 0x041e:
            case 0x04bc:
            case 0x105c:    /* ??? */

                            if (xls_unknown(ctx))
                                return(-1);    
                            break;
            default:        if (ctx->PMNI == 1) {
                                if (xls_unknown(ctx))
                                    return(-1);    
                                break;
                            }
                            printf1(ctx, "Error: unknown record type: %04x\n",sn);
                            fptr = ftell(ctx->PMFd);
                            printf1(ctx, "File position: %08lx\n",(unsigned long)(fptr - 2));
                            return(-1);   
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_short   Read next 2 bytes and return as short integer.              */
/*              Return 0 if OK, -1 if error.                                */  

int xls_short(TDAContext *ctx, short *sn)
{
    char buf[3];

    if (fread(buf,sizeof(char),2,ctx->PMFd) != 2) {
        p_err(ctx, -7,1);
        return(-1);       
    }
    memcpy(sn,buf,sizeof(*sn));   /* buf is a char array: no alignment */
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_int     Read next 4 bytes and return as integer.                    */
/*              Return 0 if OK, -1 if error.                                */  

int xls_int(TDAContext *ctx, int *n)
{
    char buf[5];

    if (fread(buf,sizeof(char),4,ctx->PMFd) != 4) {
        p_err(ctx, -7,1);
        return(-1);       
    }
    memcpy(n,buf,sizeof(*n));   /* a char array carries no alignment: copy the bytes */
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_dbl     Read next 8 bytes and return as double.                     */
/*              Return 0 if OK, -1 if error.                                */  

int xls_dbl(TDAContext *ctx, double *x)
{
    char buf[9];

    if (fread(buf,sizeof(char),8,ctx->PMFd) != 8) {
        p_err(ctx, -7,1);
        return(-1);       
    }
    memcpy(x,buf,sizeof(*x));   /* a char array carries no alignment: copy the bytes */
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  xls_rk      Read next 4 bytes and return as double (RK value).          */
/*              Return 0 if OK, -1 if error.                                */  

int xls_rk(TDAContext *ctx, double *x)
{
    int m,i;
    unsigned char buf[5],fbuf[11];
    unsigned int n;
                  
    if (fread(buf,sizeof(char),4,ctx->PMFd) != 4) {
        p_err(ctx, -7,1);
        return(-1);       
    }
    m = 0x01 & buf[0];
    i = 0x02 & buf[0];
    /**   
    tda_out("m=%d i=%d\n",m,i);        
    tda_out("buf: %02x %02x %02x %02x\n",buf[0],buf[1],buf[2],buf[3]);          
    **/   

    if (i) {                           /* integer */
        memcpy(&n,buf,sizeof(n));   /* a char array carries no alignment: copy the bytes */
        n = n >> 2;
        *x = (double)n;
    }
    else {                              /* float */
        fbuf[0] = 0x00;
        fbuf[1] = 0x00;
        fbuf[2] = 0x00;
        fbuf[3] = 0x00;
        fbuf[4] = buf[0] & 0xfc;
        fbuf[5] = buf[1];
        fbuf[6] = buf[2];
        fbuf[7] = buf[3];
        memcpy(x,fbuf,sizeof(*x));   /* a char array carries no alignment: copy the bytes */
    }
    if (m)   
        *x /= 100.0;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_unknown                                                             */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_unknown(TDAContext *ctx)
{
    char buf[10001];
    short sn;
    long fptr;

    if (xls_short(ctx, &sn))
        return(-1);       
           
    if (sn == 0)
        return(0);

    if (sn < 1 || sn > 10000) {
        printf1(ctx, "Error in length: %d\n",sn);
        fptr = ftell(ctx->PMFd);
        printf1(ctx, "File position: %08lx\n",(unsigned long)fptr);
        return(-1);
    }
    if (fread(buf,sizeof(char),(size_t)((int)sn),ctx->PMFd) != (size_t)((int)sn)) {
        p_err(ctx, -7,1);
        return(-1);    
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_boundsheet                                                          */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_boundsheet(TDAContext *ctx)
{
    unsigned char buf[257];
    short len,sn,rlen;

    if (xls_short(ctx, &len))
        return(-1);       
    if (xls_short(ctx, &sn))
        return(-1);       
    if (xls_short(ctx, &sn))
        return(-1);       
    if (xls_short(ctx, &sn))     /* options */
        return(-1);       

    if (ctx->XLSV == 8) {        /* get length of following string */
        if (xls_short(ctx, &sn))                    
            return(-1);       
        rlen = sn + 8l;
    }    
    else {
        if (fread(buf,sizeof(char),1,ctx->PMFd) != 1) {      
            p_err(ctx, -7,1);                                
            return(-1);       
        }
        sn = *buf;
        rlen = sn + 7; 
    }
    if (rlen != len) {
        printf1(ctx, "Boundsheet: error (%d, %d).\n",sn,len);
        return(-1);
    }
    if (sn < 1 || sn > 256) {
        printf1(ctx, "Boundsheet: error in length: %d\n",sn);
        return(-1);
    }
    if (fread(buf,sizeof(char),(size_t)((int)sn),ctx->PMFd) != (size_t)((int)sn)) {
        p_err(ctx, -7,1);
        return(-1);    
    }
    buf[(int)sn] = '\0';
    if (ctx->XLSRF == 0)
        printf1(ctx, "Boundsheet: %s\n",buf);
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  xls_sst                                                                 */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_sst(TDAContext *ctx)
{
    register char *p;
    char buf[10001];
    short len,sn,rtn;      
    int i,l,n,slen,rlen,rtlen;
    unsigned char c;

    if (xls_short(ctx, &len))
        return(-1);       

    if (xls_int(ctx, &n))
        return(-1);       

    if (ctx->PMProtFDef && ctx->XLSRF == 0) {
        fprintf(ctx->PMProtFd,"SST: total entries: %d\n",n);
        ctx->ProtNRec++;
    }
    if (xls_int(ctx, &n))
        return(-1);       

    if (ctx->PMProtFDef && ctx->XLSRF == 0) {
        fprintf(ctx->PMProtFd,"SST: different entries: %d\n\n",n);
        ctx->ProtNRec += 2;
    }
    rlen = 8;

    for (i = 0; i < n; ++i) {
        if (xls_short(ctx, &sn))
            return(-1);       
        slen = (int)sn;

        if (rlen >= len && sn == 0x003c) {             /* continue */
            if (xls_short(ctx, &sn))
                return(-1);       
            len += (int)sn;
            /** 
            tda_out("continue sn=%04x %d new len=%d rlen=%d i=%d n=%d\n",sn,sn,len,rlen,i,n);     
            **/ 
            if (xls_short(ctx, &sn))
                return(-1);       
            slen = (int)sn;
        }

        if (fread(buf,sizeof(char),1,ctx->PMFd) != 1) {
            p_err(ctx, -7,1);
            return(-1);    
        }
        c = (unsigned char)(buf[0]);
        if (c & 0x01)       /* 16 bit characters */
            slen *= 2;

        if (c & 0x04) {
            printf1(ctx, "SST error: cannot interpret far-east information (%04x).\n",c);
/*          return(-1);     */
        }
        rtlen = 0;
        if (c & 0x08) {
            if (xls_short(ctx, &rtn))
                return(-1);       
            rtlen = 4 * rtn;
            rlen += 2;
        }
        if (fread(buf,sizeof(char),(size_t)(slen),ctx->PMFd) != (size_t)(slen)) {
            p_err(ctx, -7,1);
            return(-1);    
        }
        buf[slen] = '\0';

        if (ctx->PMProtFDef && ctx->XLSRF == 0) {
            fprintf(ctx->PMProtFd,"%s\n",buf);
            ctx->ProtNRec++;
        }
        if (ctx->XLSRF == 0)
            ctx->XLSLN++;
    
        else if (ctx->XLSRF == 1) {
            if (ctx->XLSLNA >= ctx->XLSLN) {
                printf1(ctx, "SST error: exceeded max number of strings (%d).\n",ctx->XLSLN);
                exit(0);
            }
            l = (int)(strlen(buf));
            if (!(ctx->XLSLabel[ctx->XLSLNA] = (char *)calloc((size_t)(l + 1),sizeof(char)))) {
                p_err(ctx, -2,1);
                return(-1);            
            }
            memrq(ctx, l + 1,sizeof(char));

            p = buf;
            while (*p) {
                if (*p == '\n')
                    *p = ' ';
                p++;
            }
            strcpy(ctx->XLSLabel[ctx->XLSLNA],buf);
            ctx->XLSLNA++;
        }
        rlen += slen + 3;
        if (rtlen > 0) {
            if (fread(buf,sizeof(char),(size_t)(rtlen),ctx->PMFd) != (size_t)(rtlen)) {
                p_err(ctx, -7,1);
                return(-1);    
            }
            rlen += rtlen;
        }
    }
    if (rlen != len) {
        printf1(ctx, "SST error: %d, %d\n",rlen,len);
        return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_labelsst                                                            */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_labelsst(TDAContext *ctx)
{
    short sn,row,col,xf;
    int sst;

    if (xls_short(ctx, &sn))
        return(-1);       

    if (xls_short(ctx, &row))
        return(-1);       
    if (xls_short(ctx, &col))
        return(-1);       
    if (xls_short(ctx, &xf))
        return(-1);       
    if (xls_int(ctx, &sst))
        return(-1);       
           
    ctx->XLSLMin = imin(ctx, ctx->XLSLMin,sst);
    ctx->XLSLMax = imax(ctx, ctx->XLSLMax,sst);
       
    if (ctx->XLSRF == 2 && xls_cellok(ctx, row,col))  
        ctx->AcN[row * ctx->XLSCol + col] = sst;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_blank                                                               */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_blank(TDAContext *ctx)
{
    short sn,row,col;

    if (xls_short(ctx, &sn))
        return(-1);       
    if (xls_short(ctx, &row))
        return(-1);       
    if (xls_short(ctx, &col))
        return(-1);       
    if (xls_short(ctx, &sn))
        return(-1);       

    if (ctx->XLSRF == 0)
        xls_rowcol(ctx, row,col);
    else if (ctx->XLSRF == 2)  
        if (xls_cellok(ctx, row,col))
            ctx->AcN[row * ctx->XLSCol + col] = -3; 
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  xls_mulblank                                                            */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_mulblank(TDAContext *ctx)
{
    short sn,row,col,col1,xf;
    int i,n,len;

    if (xls_short(ctx, &sn))
        return(-1);       
    if (xls_short(ctx, &row))
        return(-1);       
    if (xls_short(ctx, &col))
        return(-1);       

    len = sn - 6;
    n = len / 2;
    if (n * 2 != len) {
        printf1(ctx, "MULBLANK error.\n");
        return(-1);
    }
    for (i = 0; i < n; ++i) {
        if (xls_short(ctx, &xf))
            return(-1);       

        if (ctx->XLSRF == 0)
            xls_rowcol(ctx, row,col + i);
        else if (ctx->XLSRF == 2)  
            if (xls_cellok(ctx, row,col))
            ctx->AcN[row * ctx->XLSCol + col] = -3; 
    }        
    if (xls_short(ctx, &col1))
        return(-1);       

    if (col + n - 1 != col1) {
        printf1(ctx, "MULBLANK1 error: %d, %d, %d.\n",col,col1,n);
        return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_mulrk                                                               */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_mulrk(TDAContext *ctx)
{
    short sn,xf,row,col,col1;
    int i,n,len;
    double x;

    if (xls_short(ctx, &sn))
        return(-1);       
    if (xls_short(ctx, &row))
        return(-1);       
    if (xls_short(ctx, &col))
        return(-1);       

    len = sn - 6;
    n = len / 6;
    if (n * 6 != len) {
        printf1(ctx, "MULRK error.\n");
        return(-1);
    }
    for (i = 0; i < n; ++i) {
        if (xls_short(ctx, &xf))
            return(-1);       
        if (xls_rk(ctx, &x))
            return(-1);       

        if (ctx->XLSRF == 0) 
            xls_rowcol(ctx, row,col + i);
        else if (ctx->XLSRF == 2)
            xls_putx(ctx, row,col + i,x); 
    }        
    if (xls_short(ctx, &col1))
        return(-1);       

    if (col + n - 1 != col1) {
        printf1(ctx, "MULRK1 error: %d, %d, %d.\n",col,col1,n);
        return(-1);
    }
    /**
    printf1(ctx, "MulRK:     %4d %4d %4d\n",row,col,col1);          
    **/
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_formula                                                             */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_formula(TDAContext *ctx)
{
    short sn,row,col,xf,opt;
    int len;
    char buf[100001];
    double x;

    if (xls_short(ctx, &sn))
        return(-1);       
    len = sn;

    if (xls_short(ctx, &row))
        return(-1);       
    if (xls_short(ctx, &col))
        return(-1);       
    if (xls_short(ctx, &xf))
        return(-1);       

    if (xls_dbl(ctx, &x))
        return(-1);       

    if (xls_short(ctx, &opt))
        return(-1);       

    len -= 16;
    if (len > 0) {
        if (fread(buf,sizeof(char),(size_t)(len),ctx->PMFd) != (size_t)(len)) {
            p_err(ctx, -7,1);
            return(-1);    
        }
    }
    if (ctx->XLSRF == 0)
        xls_rowcol(ctx, row,col);
    else if (ctx->XLSRF == 2)
        xls_putx(ctx, row,col,x);
    /**
    printf1(ctx, "Formula:   %4d %4d xf=%d opt=%04x value: %lg\n",row,col,xf,opt,x);
    **/
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_number                                                              */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_number(TDAContext *ctx)
{
    short sn,row,col,xf;
    double x;

    if (xls_short(ctx, &sn))
        return(-1);       

    if (xls_short(ctx, &row))
        return(-1);       
    if (xls_short(ctx, &col))
        return(-1);       
    if (xls_short(ctx, &xf))
        return(-1);       
    if (xls_dbl(ctx, &x))
        return(-1);       

    if (ctx->XLSRF == 0)
        xls_rowcol(ctx, row,col);
    else if (ctx->XLSRF == 2)
        xls_putx(ctx, row,col,x);

    /**
    printf1(ctx, "Number:    %4d %4d xf=%d value: %lg\n",row,col,xf,x);     
    **/
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_rkvalue                                                             */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_rkvalue(TDAContext *ctx)
{
    short sn,row,col,xf;
    double x;

    if (xls_short(ctx, &sn))
        return(-1);       

    if (xls_short(ctx, &row))
        return(-1);       
    if (xls_short(ctx, &col))
        return(-1);       
    if (xls_short(ctx, &xf))
        return(-1);       
    if (xls_rk(ctx, &x))
        return(-1);       

    if (ctx->XLSRF == 0)
        xls_rowcol(ctx, row,col);
    else if (ctx->XLSRF == 2)
        xls_putx(ctx, row,col,x);

    /**
    printf1(ctx, "RK-Value:  %4d %4d xf=%d value: %lg\n",row,col,xf,x);  
    **/
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_rowcol                                                              */
/*  Return 0 if OK, -1 if error.                                            */  

void xls_rowcol(TDAContext *ctx, int row,int col)
{
    if (ctx->XLSRF)
        return;

    ctx->XLSRowMin = imin(ctx, ctx->XLSRowMin,row);
    ctx->XLSRowMax = imax(ctx, ctx->XLSRowMax,row);
    ctx->XLSColMin = imin(ctx, ctx->XLSColMin,col);
    ctx->XLSColMax = imax(ctx, ctx->XLSColMax,col);
}

/* ------------------------------------------------------------------------ */
/*  xls_putx(x)     Put x into AcX.                                         */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */  

/* ------------------------------------------------------------------------ */
/*  xls_cellok(row,col)  Is this cell inside the table pass 1 measured?     */
/*                       Return 1 if it is, 0 if not.                      */

int xls_cellok(TDAContext *ctx, int row, int col)
{
    /*  AcN/AcX are allocated for XLSRow * XLSCol cells, where XLSRow and
        XLSCol come from the row/column maxima found in the FIRST pass.
        The second pass indexes them by absolute row and column, so a
        record naming a cell beyond those maxima writes outside the
        allocation.

        That is not hypothetical: readxl's type-me.xls does it, and the
        write corrupted the heap -- "munmap_chunk(): invalid pointer",
        abort.  Inside this R package TDA runs IN-PROCESS, so an abort
        takes the user's whole session down, with no error to catch.  A
        spreadsheet that cannot be read is a nuisance; one that kills
        the session while being read is not acceptable, so every write
        is checked and an out-of-range cell is skipped and counted.  */
    if (row < 0 || col < 0 || row >= ctx->XLSRow || col >= ctx->XLSCol) {
        ctx->XLSNSkip++;
        return(0);
    }
    return(1);
}

void xls_putx(TDAContext *ctx, int row,int col,double x)
{
    int i;

    if (ctx->XLSRF == 2) {
        if (!xls_cellok(ctx, row,col))
            return;
        i = row * ctx->XLSCol + col;
        ctx->AcX[i] = x;
        ctx->AcN[i] = -2;
        ctx->AcI[row] += 1;
    }
}

/* ------------------------------------------------------------------------ */
/*  xls_check(j)    Return 1 if column j contains only numerical entries,   */
/*                         2 if column j contains at least one string,      */
/*                         0 otherwise                                      */
/*                  Take into account only rows with at least one numerical */
/*                  entry.                                                  */

int xls_check(TDAContext *ctx, int j)
{
    int i,ij,nn,ns;

    nn = ns = 0;
    for (i = 0; i <= ctx->XLSRowMax; ++i) {
        if (ctx->AcI[i] > 0) {
            ij = i * ctx->XLSCol + j;
            if (ctx->AcN[ij] == -2)
                nn++;
            else if (ctx->AcN[ij] >= 0)
                ns++;
        }
    }
    if (ns == 0 && nn == 0)
        return(0);
    else if (ns > 0)
        return(2);
    else
        return(1);
}

/* -##--------------------------------------------------------------------- */
/*  rucinet         Read a Ucinet file.                                     */
/*                                                                          */
/*                  rxls(                                                   */
/*                      df=...,     write data to output file               */
/*                      fmt=...,    print format, def. 10.4                 */
/*                      n=...,..,   number of columns, def. 0               */
/*                      prn=...,    type of output, def. 0                  */
/*                                  0 = print as square matrix              */
/*                                  1 = print as edge list, only values     */
/*                                      greater than zero.                  */  
/*                  ) = fname;      name of Ucinet file.                    */
/*                                                                          */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int rucinet(TDAContext *ctx)  
{
    char buf[50];
    int err,nc,nrec,i,j;
    double x;    

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 7,10,1)) {   /* get parameters */
        p_clean(ctx);
        return(-1);   
    }   
    printf1(ctx, "Reading Ucinet file: %s\n\n",ctx->PMFdName);
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    /**********************
    fseek(PMFd,0x19da,0);
    i = 1;
    while (fread(buf,sizeof(char),1,PMFd) == 1) {                        
        if (buf[0])
            tda_out("%c",buf[0]);
        else {
            if (++i >= 43) {
                newline(ctx);
                i = 0;
            } 
            else
                tda_out(" ");
        }
    }
    *******************************/

    if (ctx->PMF1Def == 0) {
        printf1(ctx, "Error: need output file (df parameter).\n");
        goto RUCIFin;
    }            

    /*** 
    fread(buf,sizeof(char),27,PMFd);                               
    while (fread(buf,sizeof(char),20,PMFd) == 20)                          
        tda_out("%s,\n",buf);
    goto RUCIFin;
    ***/

    nc = nrec = 0;
    i = j = 1;

    while (fread(buf,sizeof(char),4,ctx->PMFd) == 4) {                        
        x = (double)(*(float *)buf);
        if (ctx->PMPRNO == 1) {
            if (x > 0.0) {
                fprintf(ctx->PMF1d,"%6d %6d ",i,j);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,(double)x);
                fprintf(ctx->PMF1d,"\n");
                nrec++;
            }
            if (++j > ctx->PMN) {
                i++;
                j = 1;
            }
        }
        else {
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,(double)x);
            if (++nc >= ctx->PMN) {
                fprintf(ctx->PMF1d,"\n");
                nrec++;
                nc = 0;
            }
        }
    }
    printf1(ctx, "\n%d records written to: %s\n",nrec,ctx->PMF1dName);
    err = 0;

RUCIFin:
    p_clean(ctx);
    return(err);
}





