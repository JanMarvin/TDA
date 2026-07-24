/****************************************************************************/
/*  t_spss                                                                  */
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
#include "t_gdat.h"
#include "t_sort.h"
#include "t_alloc.h"
#include "t_eval.h"
#include "t_eval1.h"
#include "t_gf.h"
#include "tda_context.h"

const char DIG30_data[] = {'0','1','2','3','4','5','6','7','8','9',
                'A','B','C','D','E','F','G','H','I','J',
                'K','L','M','N','O','P','Q','R','S','T'};

const char S5_data[] =
"uvwxyz .<(+0&[]!$*);^-/|,%_>?`:#@'=\"000000~000000000000000000000{}\\0000000000000";


const char S7_data[] =
"0001F/TDA Version 6.x";

const char S6_data[] =
"00000000000000000000000000000000000000000000000000000000SPSSPORTA8/190000006/000";

const char S4_data[] =
"0000000000000000000000000123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrst";

const char S3_data[] =
"0200002'220'&)3000#0000000000000000000000000000000000000000000000000000000000000";

const char S2_data[] =
"00000-0000-0000-0000--------------------!3#))0303300/240&),%00000000000000000000";

const char S1_data[] =
"00000@0000@0000@0000@@@@@@@@@@@@@@@@@@@@ASCII SPSS PORT FILE                    ";

/*  functions in t_spss.c */

int rd_spss(TDAContext *ctx); 
int buf_update(TDAContext *ctx, int init);
void prstr(TDAContext *ctx, char *p, int n);
char *getstr(TDAContext *ctx, int f);
int getdigit(TDAContext *ctx, char *p, int *err);
char *inum(TDAContext *ctx, char *p, int *n);
char *dnum(TDAContext *ctx, char *p, double *x, int *mv);
void digit_err(TDAContext *ctx, char *p);
int wr_spss(TDAContext *ctx);
void pnum1(TDAContext *ctx, FILE *fd,int n);
void pnum(TDAContext *ctx, FILE *fd,int n);
void pfnum(TDAContext *ctx, FILE *fd,double x);
void pstring(TDAContext *ctx, FILE *fd,char *s,int m);
void check_nl(TDAContext *ctx, FILE *fd);
int rd_spss1(TDAContext *ctx); 
int rd_dat(TDAContext *ctx, int n,char *buf);
int sav_sval(TDAContext *ctx, int *sval);
int sav_geti(TDAContext *ctx, int arch,int *err);
double sav_getd(TDAContext *ctx, int arch,int *err);
int wr_spss1(TDAContext *ctx);
void sav_puti(TDAContext *ctx, int n);
void sav_putd(TDAContext *ctx, double x);
int get_dvarp(TDAContext *ctx, int *np,char *vname);
void make_arcd(TDAContext *ctx, int fn,char *fname,int len,int nrec,int vn);

/*--------------------------------------------------------------------------*/
/*  Global variables                                                        */

#define RLEN 80             /* Record length                                */
#define HLEN 200            /* Header length                                */
#define TLEN 256            /* Length of translation table                  */
#define DTLEN 10            /* Length of creation date string               */
#define TTLEN  8            /* Length of creation time string               */




/* ------------------------------------------------------------------------ */
/*  rd_spss     Read SPSS Export file and create an internal data matrix.   */
/*                                                                          */
/*              rspss(                                                      */
/*                  len= ...,           record length of input file, def.80 */
/*                  noc=...,            # of records to read, def. all      */
/*                  msys=...,           new sys miss value code, def. -5    */  
/*                  fmt=...,            new print format                    */
/*                  df=...,             write data directly to output file  */
/*                  dvar(fn=...)=...,   create var description file         */
/*                  arcd=...,           archive description file            */
/*                  ns=...,             partition of string variables       */
/*              ) = file_name;                                              */
/*                                                                          */
/*              arcd can only be used in connection with df and dvar.       */
/*                                                                          */
/*              df=... can be used as dfa to append data.                   */
/*              dvara(...)=... can be used to append.                       */
/*              arcda=... can be used to append.                            */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int rd_spss(TDAContext *ctx)  
{
    register int i,j,k;
    register char *p;
    int err,idxn,nvar,nvar1,nvar2,m,n,vtyp,len,prec,mv,nmiss,nrec;
    int aptr,aoff,avint,avmax,avmin,avttyp,vlflag,pfmt1,pfmt2,nmiss1;
    int np,nrec1,nrec2,mxnoc,aflag,vrec,vn = 0,len1;
    short *vint,*vttyp;
    char vname[VNLMax + 1],vlabel[VLLMax + 1],vdef[400];
    double tol,a,x,y,*vmax,*vmin;

    tol = 10000.0 * ctx->EPSI;
    aptr = aoff = avint = avmin = avmax = avttyp = 0;
    idxn = get_nidx(ctx);      /* index to first new variable */
    ctx->VLabelLen = 0;
    err = -1;

    if (check_cmd(ctx, 0))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 5,10,1)) {   /* get parameters */
        p_clean(ctx);
        return(-1);   
    }   

    aflag = 0;
    if (ctx->PMARCFDef && ctx->PMDVARFDef && ctx->PMF1Def)
        aflag = 1;

    printf1(ctx, "Reading SPSS export file: %s\n",ctx->PMFdName);
    if (ctx->DMDef) {
        printf1(ctx, "Error: a data matrix already exists.\n");
        p_clean(ctx);
        return(-1);
    }
    if (ctx->PMLEN >= 80 && ctx->PMLEN <= 82) {
        ctx->RSRILen = ctx->PMLEN;
        printf1(ctx, "Assuming fixed record length: %d bytes.\n",ctx->RSRILen);
    }
    else  
        ctx->RSRILen = 0;

    if (ctx->PMNOCFlg)
        printf1(ctx, "Maximum number of cases: %d\n",ctx->PMNOC);

    if (ctx->PMRHSTRA && ctx->PMF1Def == 0)
        printf1(ctx, "Record selection (vsel): %s\n",ctx->PMRHSTR + 5);

    /* allocate read buffer */

    if (!(ctx->RSPBuf = (char *) calloc(2 * RLEN + 6,sizeof(char)))) {
        p_err(ctx, -2,1);
        goto RSPFin;
    }
    memrq(ctx, 2 * RLEN + 6,1);
    ctx->RSPBufA = 2 * RLEN + 6;
             
    ctx->RSPPtr = ctx->RSPBuf;
    ctx->RSPEOF = ctx->RSPCnt = 0;
    if (buf_update(ctx, 1)) 
        goto RSPFin;
           
    /*  Get the 5 * 40 bytes header */

    printf1(ctx, "\nHeader: ");          
    for (i = 1; i <= 5; ++i) {
        prstr(ctx, ctx->RSPPtr,40);
        ctx->RSPPtr += 40;
        ctx->RSPCnt -= 40;
        printf1(ctx, "\n        ");
        if (buf_update(ctx, 0))
            goto RSPFin;
    }

    /*  Skip translation table, length is HLEN = 256 chars */
  
    for (i = 1; i <= 3; ++i) {
        ctx->RSPPtr += RLEN;
        ctx->RSPCnt -= RLEN;
        if (buf_update(ctx, 0))
            goto RSPFin;
    }
    ctx->RSPPtr += 16;
    ctx->RSPCnt -= 16;
    if (buf_update(ctx, 0))
        goto RSPFin;
           
    /*  Read the SPSSPORT string, length is 8 characters */
   
    if (strncmp(ctx->RSPPtr,"SPSSPORT",8)) {
        printf1(ctx, "\nWarning: can't locate SPSSPORT string.\n");
        /*******************
        prstr(ctx, RSPPtr,40);                                       
        printf1(ctx, "\n");
        goto RSPFin;
        *******************/
    }
    ctx->RSPPtr += 8;
    ctx->RSPCnt -= 8;

    /*  Read the version code, length is 1 char */
  
    if (buf_update(ctx, 0))
        goto RSPFin;
    printf1(ctx, "\nFile-format version code: %c <%02x>",*ctx->RSPPtr,*ctx->RSPPtr);
    ctx->RSPPtr++;
    ctx->RSPCnt--;
  
    /*  Read creation date (10 char) and time */

    printf1(ctx, "\nCreation date: ");
    prstr(ctx, ctx->RSPPtr + 2,DTLEN - 2);
    ctx->RSPPtr += DTLEN;
    ctx->RSPCnt -= DTLEN;

    printf1(ctx, "\nCreation time: ");
    prstr(ctx, ctx->RSPPtr + 2,TTLEN - 2);
    ctx->RSPPtr += TTLEN;
    ctx->RSPCnt -= TTLEN;

    if (buf_update(ctx, 0))
        goto RSPFin;

    if (*ctx->RSPPtr == '1') {        /* get originating software */
        printf1(ctx, "\nSoftware: ");
        ctx->RSPCnt--;
        ctx->RSPPtr++;
        ctx->RSPPtr = getstr(ctx, 0);
        if (ctx->RSPPtr == NULL)
            goto RSPFin;
    }
    if (*ctx->RSPPtr == '2') {        /* get originating installation */
        if (buf_update(ctx, 0))
            goto RSPFin;
        printf1(ctx, "\nInstallation: ");
        ctx->RSPCnt--;
        ctx->RSPPtr++;
        ctx->RSPPtr = getstr(ctx, 0);
        if (ctx->RSPPtr == NULL)
            goto RSPFin;
    }
    if (*ctx->RSPPtr == '3') {        /* get file label */
        if (buf_update(ctx, 0))
            goto RSPFin;
        printf1(ctx, "\nLabel: ");
        ctx->RSPCnt--;
        ctx->RSPPtr++;
        ctx->RSPPtr = getstr(ctx, 0);
        if (ctx->RSPPtr == NULL)
            goto RSPFin;
    }
    printf1(ctx, "\n");           

    /*  get the number of variables */
   
    if (*ctx->RSPPtr != '4') {
        printf1(ctx, "\nError: can't find number of variables.\n");
        goto RSPFin;
    }
    if (buf_update(ctx, 0))
        goto RSPFin;

    ctx->RSPCnt--;
    ctx->RSPPtr = inum(ctx, ++ctx->RSPPtr,&nvar);
    if (ctx->RSPPtr == NULL)
        goto RSPFin;

    printf1(ctx, "\nNumber of variables: %d",nvar);
    if (nvar <= 0) {
        printf1(ctx, "\nError in number of variables.\n");
        goto RSPFin;
    }
    else if (nvar > ctx->MaxNV) {
        printf1(ctx, "\nError: exceeded max number of variables.\n");
        goto RSPFin;
    }
    if (*ctx->RSPPtr == '5') {        /* get number of base-30 digits */
        ctx->RSPCnt--;
        ctx->RSPPtr = inum(ctx, ++ctx->RSPPtr,&n);
        if (ctx->RSPPtr == NULL)
            goto RSPFin;
        printf1(ctx, "\nPrecison (base-30 digits): %d",n);
        if (n > 0)
            ctx->RSPMaxD = n;
    }

    printf1(ctx, "\nCase-weight variable: ");          

    if (*ctx->RSPPtr == '6') {        /* get case-weight variable */
        ctx->RSPCnt--;
        ctx->RSPPtr++;
        ctx->RSPPtr = getstr(ctx, 0);
        if (ctx->RSPPtr == NULL)
            goto RSPFin;
    }
    else  
        printf1(ctx, "not defined.");
    printf1(ctx, "\n");           

    if (!(vttyp = (short *) calloc((size_t)(nvar),sizeof(short)))) {
        p_err(ctx, -2,1);
        goto RSPFin;
    }
    avttyp = nvar;
    memrq(ctx, nvar,sizeof(short));

    /*  Read the variable descriptions */

    nvar1 = 0;      /* number of string variables */
    nvar2 = 0;      /* number of variables with unknown format */

    for (i = 0; i < nvar; ++i) {

        if (buf_update(ctx, 0))
            goto RSPFin;

        if (*ctx->RSPPtr != '7') {
            printf1(ctx, "\nError: can't find description of variable %d.\nFound: ",i);
            prstr(ctx, ctx->RSPPtr - 10,40);
            printf1(ctx, "\n");
            goto RSPFin;
        }

        /* get type of variable: 0 numerical, else string */

        ctx->RSPCnt--;
        ctx->RSPPtr = inum(ctx, ++ctx->RSPPtr,&vtyp);
        if (ctx->RSPPtr == NULL)
            goto RSPFin;

        vttyp[i] = (short)(vtyp); 

        if (vtyp > 0)       /* alphanumerical variable */
            nvar1++;

        /* get name of variable */
  
        ctx->RSPPtr = inum(ctx, ctx->RSPPtr,&n);
        if (ctx->RSPPtr == NULL)
            goto RSPFin;

        p = vname; 
        for (j = 0; j < n; ++j) {
            if (j < VNLMax)  
                *p++ = get_vnchar(ctx, *ctx->RSPPtr);
            ctx->RSPPtr++;
        }
        *p = '\0';
        ctx->RSPCnt -= n;

        /*  get print formats (not used) */

        for (j = 1; j <= 3; ++j) {
            ctx->RSPPtr = inum(ctx, ctx->RSPPtr,&n);
            if (ctx->RSPPtr == NULL)
                goto RSPFin;
        }

        /* get write formats; ignored */
   
        ctx->RSPPtr = inum(ctx, ctx->RSPPtr,&n);
        if (ctx->RSPPtr == NULL)
            goto RSPFin;

        if (n != 1 && n != 5)
            nvar2++;
           
        ctx->RSPPtr = inum(ctx, ctx->RSPPtr,&len);
        if (ctx->RSPPtr == NULL)
            goto RSPFin;

        ctx->RSPPtr = inum(ctx, ctx->RSPPtr,&prec);
        if (ctx->RSPPtr == NULL)
            goto RSPFin;

        /*  get missing value codes */
  
        n = 0;
        for (j = 0; j < 3; ++j) {

            switch (*ctx->RSPPtr) {

                case '8':   ctx->RSPCnt--;
                            ctx->RSPPtr = dnum(ctx, ++ctx->RSPPtr,&x,&mv);
                            n++;
                            break;
                case '9':   ctx->RSPCnt--;
                            ctx->RSPPtr = dnum(ctx, ++ctx->RSPPtr,&x,&mv);
                            n++;
                            break;
                case 'A':   ctx->RSPCnt--;
                            ctx->RSPPtr = dnum(ctx, ++ctx->RSPPtr,&x,&mv);
                            n++;
                            break;
                case 'B':   ctx->RSPCnt--;
                            ctx->RSPPtr = dnum(ctx, ++ctx->RSPPtr,&x,&mv);
                            if (ctx->RSPPtr == NULL)
                                break;
                            ctx->RSPPtr = dnum(ctx, ctx->RSPPtr,&y,&mv);
                            n++;
                            break;
                default:    goto CONT;
            }
            if (ctx->RSPPtr == NULL)
                goto RSPFin;
        }
CONT:
        snprintf(vdef,sizeof(vdef),"%s<4>[0.0]",vname);

        vlflag = 0;
        if (*ctx->RSPPtr == 'C') {   /* Var label is optional */

            ctx->RSPCnt--;
            ctx->RSPPtr = inum(ctx, ++ctx->RSPPtr,&n);
            if (ctx->RSPPtr == NULL)
                goto RSPFin;

            if (buf_update(ctx, 0))
                goto RSPFin;

            p = vlabel;
            for (j = 0; j < n; ++j) {
    
                if (ctx->RSPCnt < 10) {
                    if (buf_update(ctx, 0))
                        goto RSPFin;
                }
                if (j < VLLMax)
                    *p++ = *ctx->RSPPtr;
                ctx->RSPPtr++;
                ctx->RSPCnt--;
            }
            *p = '\0';
            p = vlabel;
            while (*p) {
                if (*p == '(')
                    *p = '[';
                else if (*p == ')')
                    *p = ']';
                p++;
            }
            p = vdef + strlen(vdef);
            snprintf(p,sizeof(vdef) - (size_t)(p - vdef),"(%s)",vlabel);
            vlflag = 1;
        }
        p = vdef + strlen(vdef);
        snprintf(p,sizeof(vdef) - (size_t)(p - vdef),"=spss(%d)",vtyp);
        if (save_var(ctx, vdef,0)) {       /* save variable definition */
            printf1(ctx, "\nError: can't save variable definitions.\n");
            goto RSPFin;
        }
    }
    printf1(ctx, "\nNumber of string variables: %d",nvar1);
    printf1(ctx, "\nNumber of variables with unknown format: %d\n",nvar2);
          
    /*  Read value labels (optional). */

    vlflag = 0;

    if (*ctx->RSPPtr == 'D') {

        printf1(ctx, "\nReading value labels.\n");
               
        /* allocate memory for file pointers */

        if (!(ctx->RSPVFPtr = (long *) calloc((size_t)(nvar),sizeof(long)))) {
            p_err(ctx, -2,1);
            goto RSPFin;
        }
        aptr = nvar;
        memrq(ctx, nvar,sizeof(long));

        if (!(ctx->RSPVFOff = (int *) calloc((size_t)(nvar),sizeof(int)))) {
            p_err(ctx, -2,1);
            goto RSPFin;
        }
        aoff = nvar;
        memrq(ctx, nvar,sizeof(int));

        vlflag = 1;

        while (*ctx->RSPPtr == 'D') {

            if (buf_update(ctx, 0))
                goto RSPFin;
            ctx->RSPPtr++;
            ctx->RSPCnt--;
                          
            ctx->RSPPtr = inum(ctx, ctx->RSPPtr,&n);
            if (ctx->RSPPtr == NULL)
                goto RSPFin;

            len = 0;
            for (j = 1; j <= n; ++j) {
                if (buf_update(ctx, 0))
                    goto RSPFin;

                ctx->RSPPtr = inum(ctx, ctx->RSPPtr,&m);
                if (ctx->RSPPtr == NULL)
                    goto RSPFin;

                p = vname;
                for (i = 0; i < m; ++i)  
                    *p++ = *ctx->RSPPtr++;
                *p = '\0';
                ctx->RSPCnt -= m;

                k = -1;
                for (i = 0; i < nvar; ++i) {
                    if (!strcmp(vname,ctx->VName[i])) {
                        k = i;
                        break;
                    }
                }
                if (k >= 0) {
                    ctx->RSPVFOff[k] = -1;
                    len = ctx->VTyp[k];
                }
            }
            if (buf_update(ctx, 0))
                goto RSPFin;

            for (i = 0; i < nvar; ++i) {
                if (ctx->RSPVFOff[i] < 0) {
                    ctx->RSPVFOff[i] = ctx->RSPFOff;
                    ctx->RSPVFPtr[i] = ctx->RSPFPtr;
                }
            }
            ctx->RSPPtr = inum(ctx, ctx->RSPPtr,&n);
            if (ctx->RSPPtr == NULL)
                goto RSPFin;

            for (j = 1; j <= n; ++j) {
                if (buf_update(ctx, 0))
                    goto RSPFin;

                if (len != 1) {       /* numerical variable */
                    ctx->RSPPtr = dnum(ctx, ctx->RSPPtr,&x,&mv);
                    if (ctx->RSPPtr == NULL)  
                        goto RSPFin;
                }
                else {              /* string variable */
                    ctx->RSPPtr = inum(ctx, ctx->RSPPtr,&m);
                    if (ctx->RSPPtr == NULL)
                        goto RSPFin;

                    ctx->RSPPtr += m;
                    ctx->RSPCnt -= m;
                }
                ctx->RSPPtr = inum(ctx, ctx->RSPPtr,&m);
                if (ctx->RSPPtr == NULL)
                    goto RSPFin;

                ctx->RSPPtr += m;
                ctx->RSPCnt -= m;
            }
        }
    }

    /* Read textual information documenting the file (optional) */
  
    if (*ctx->RSPPtr == 'E') {
        printf1(ctx, "\nFile documentation\n");
        ctx->RSPCnt--;
        ctx->RSPPtr = inum(ctx, ++ctx->RSPPtr,&n);
        if (ctx->RSPPtr == NULL)
            goto RSPFin;

        for (j = 1; j <= n; ++j) {
            if (buf_update(ctx, 0))
                goto RSPFin;
            if ((ctx->RSPPtr = getstr(ctx, 0)) == NULL)  
                goto RSPFin;
            printf1(ctx, "\n");          
        }
        printf1(ctx, "\n"); 
    }

    if (buf_update(ctx, 0))
        goto RSPFin;

    ctx->RSPFPtr1 = ctx->RSPFPtr;     /* remember file position for data */
    ctx->RSPFOff1 = ctx->RSPFOff;

    if (*ctx->RSPPtr != 'F') {               /* check for data stream */

        printf1(ctx, "Error: can't find data stream.\n");
        goto RSPFin;
    }
    ctx->RSPCnt--;
    ctx->RSPPtr++;
   
    /* allocate memory to check for maximum and integer values */
   
    if (!(vmax = (double *) calloc((size_t)(nvar),sizeof(double)))) {
        p_err(ctx, -2,1);
        goto RSPFin;
    }
    avmax = nvar;
    memrq(ctx, nvar,sizeof(double));
   
    if (!(vmin = (double *) calloc((size_t)(nvar),sizeof(double)))) {
        p_err(ctx, -2,1);
        goto RSPFin;
    }
    avmin = nvar;
    memrq(ctx, nvar,sizeof(double));

    if (!(vint = (short *) calloc((size_t)(nvar),sizeof(short)))) {
        p_err(ctx, -2,1);
        goto RSPFin;
    }
    avint = nvar;
    memrq(ctx, nvar,sizeof(short));
   
    printf1(ctx, "Reading data to check variables.\n");

    if (ctx->SILENTFlg < 2)
        printfe(ctx, "Reading: %s\n",ctx->PMFdName);

    nrec = 0;       /* number of records */
    nmiss = 0;      /* number of system missing values */
                     
    while (*ctx->RSPPtr && *ctx->RSPPtr != 'Z') {    /* this is the end of file marker */

        for (i = 0; i < nvar; ++i) {
            if (buf_update(ctx, 0))
                goto RSPFin;

            if (vttyp[i] > 0) {        /* string variable */
                ctx->RSPPtr = inum(ctx, ctx->RSPPtr,&m);
                if (ctx->RSPPtr == NULL)
                    goto RSPFin;
                if (m > vttyp[i])  
                    printf1(ctx, "Warning: variable %d in record %d (type=%d fnd=%d).\n",
                                            i+1,nrec+1,vttyp[i],m);

                while (m >= ctx->RSPCnt) {
                    m -= ctx->RSPCnt;
                    ctx->RSPPtr += ctx->RSPCnt;
                    ctx->RSPCnt = 0;
                    if (buf_update(ctx, 0))
                        goto RSPFin;
                }
                ctx->RSPPtr += m;          
                ctx->RSPCnt -= m;
            }
            else {
                ctx->RSPPtr = dnum(ctx, ctx->RSPPtr,&x,&mv);
                if (ctx->RSPPtr == NULL)
                    goto RSPFin;

                if (mv)  
                    nmiss++;    
                else {
                    if (vmax[i] < x)
                        vmax[i] = x;
                    if (vmin[i] > x)
                        vmin[i] = x;
           
                    y = (int)x;
                    if (fabs(y - x) > tol)      
                        vint[i] = 1;
                }
            }
        }
        nrec++;
        prn_message(ctx, nrec,0,0);
   
        if (ctx->PMNOCFlg && nrec >= ctx->PMNOC)
            break;
    
        if (ctx->RSPEOF > RLEN)
            break;
    }
    prn_message(ctx, nrec,1,0);

    printf1(ctx, "Read %d records.\n",nrec);
    if (nrec <= 0) {
        printf1(ctx, "Error.\n");
        goto RSPFin;
    }
    if (*ctx->RSPPtr != 'Z')  
        printf1(ctx, "Warning: can't find end of file character (Z).\n");
    printf1(ctx, "Number of system missing values: %d\n\n",nmiss);

    /* set appropriate storage size and write formats */

    for (i = 0; i < nvar; ++i) {

        if (vttyp[i] > 0) {
            pfmt1 = pfmt2 = 0; 
        }
        else {
            x = fabs(vmin[i]);
            y = fabs(vmax[i]);
            if (x < y)
                a = y;
            else
                a = x;

            if (vint[i])
                ctx->VSLen[i] = 8;
            else if ((int)a < 128)
                ctx->VSLen[i] = 1;
            else if ((int)a < 32000)
                ctx->VSLen[i] = 2;
            else 
                ctx->VSLen[i] = 5;

            x = vmin[i];
            y = vmax[i];

            if (vint[i]) {
                x = ceil(x);
                y = floor(y);
            }
            pfmt2 = 0;
            a = 10.0;
            for (j = 2; j <= 12; ++j) {
                if (-a < x && y < a * 10.0) {
                    pfmt1 = j;
                    break;
                }
                a *= 10.0;
            }
            if (vint[i]) {
                if (ctx->PMFmtF) {
                    pfmt1 = ctx->PMFmt1;
                    pfmt2 = ctx->PMFmt2;
                }
                else {
                    pfmt1 += 7;
                    pfmt2 = 6;
                }
            }
        }
        makefmt(ctx, &pfmt1,&pfmt2,ctx->VPFmtS[i],VPFmtSLen,0,ctx->SEPC,0);
        ctx->VPFmt1[i] = (short)pfmt1;
        ctx->VPFmt2[i] = (short)pfmt2;
    }
    prn_var(ctx, idxn);      /* print list of new variables */
    newline(ctx);     

    /* if requested create the variable description file */

    len = 0;
    if (ctx->PMDVARFDef) {
  
        printf1(ctx, "Creating variable description file: %s\n",ctx->PMDVARFName);
            
        fprintf(ctx->PMDVARFd,"# variable description file based on: %s\n",ctx->PMFdName);

        vrec = 2,
        len1 = vn = 0;

        for (i = 0; i < nvar; ++i) {

            fprintf(ctx->PMDVARFd,"%s ",ctx->VName[i]);
            fprnchar(ctx, ctx->PMDVARFd,' ',ctx->VNameLen - (int)strlen(ctx->VName[i]),0);
            fprintf(ctx->PMDVARFd," %3d %4d ",ctx->PMDVARFN,len);
            if (ctx->VTyp[i] != 1) {
                fprintf(ctx->PMDVARFd,"%4d",(int)ctx->VPFmt1[i]);
                if ((int)ctx->VPFmt2[i] > 0)
                    fprintf(ctx->PMDVARFd,".%-2d ",(int)ctx->VPFmt2[i]);
                else
                    fprintf(ctx->PMDVARFd,"    ");
            }
            else  
                fprintf(ctx->PMDVARFd,"%4d    ",-vttyp[i]);

            if (ctx->VLabel[i] != NULL)
                fprintf(ctx->PMDVARFd,"%s",ctx->VLabel[i]);
            fprintf(ctx->PMDVARFd,"\n");
            vrec++;
            vn++;
            if (ctx->VTyp[i] != 1)
                len += (int)ctx->VPFmt1[i] + 1;
            else {
                len1 = len;
                len += vttyp[i] + 1;
            }

            /* add value labels */

            if (vlflag && ctx->RSPVFPtr[i] > 0) {
                   
                ctx->RSPEOF = 0;
                fseek(ctx->PMFd,ctx->RSPVFPtr[i],0);
                ctx->RSPPtr = ctx->RSPBuf;
                ctx->RSPCnt = 0;
                if (buf_update(ctx, 0))
                    goto RSPFin;
                ctx->RSPPtr += ctx->RSPVFOff[i];
                ctx->RSPCnt -= ctx->RSPVFOff[i];
                if (buf_update(ctx, 0))
                    goto RSPFin;
                ctx->RSPPtr = inum(ctx, ctx->RSPPtr,&m);
                if (ctx->RSPPtr == NULL)
                    goto RSPFin;

                for (j = 1; j <= m; ++j) {
                    if (buf_update(ctx, 0))
                        goto RSPFin;

                    ctx->RSPPtr = dnum(ctx, ctx->RSPPtr,&x,&mv);
                    if (ctx->RSPPtr == NULL)  
                        goto RSPFin;

                    fprintf(ctx->PMDVARFd," %4g ",x);

                    ctx->RSPPtr = inum(ctx, ctx->RSPPtr,&n);
                    if (ctx->RSPPtr == NULL)
                        goto RSPFin;

                    for (k = 0; k < n; ++k)  
                        fprintf(ctx->PMDVARFd,"%c",*ctx->RSPPtr++);
                    fprintf(ctx->PMDVARFd,"\n");
                    ctx->RSPCnt -= n;
                    vrec++;
                }
            }
            if (ctx->VTyp[i] == 1) {     /* check for partition */
                                  
                n = get_dvarp(ctx, &np,ctx->VName[i]); 

                if (n > 0) {
                    if (n >= vttyp[i]) {
                        if (np == 2)
                            n = vttyp[i];
                        else
                            n = 0;
                    }
                }
                if (n > 0) {
                    m = j = 0;
                    while (m < vttyp[i]) {
                        fprintf(ctx->PMDVARFd,"%s%02d ",ctx->VName[i],j);
                        fprnchar(ctx, ctx->PMDVARFd,' ',ctx->VNameLen - (int)strlen(ctx->VName[i]) - 2,0);
                        fprintf(ctx->PMDVARFd," %3d %4d ",ctx->PMDVARFN,len1);
                        if (np == 1)
                            fprintf(ctx->PMDVARFd,"%4d    ",-n);
                        else
                            fprintf(ctx->PMDVARFd,"%4d    ",n);
                        if (ctx->VLabel[i] != NULL)
                            fprintf(ctx->PMDVARFd,"%s, part %02d",ctx->VLabel[i],j);
                        fprintf(ctx->PMDVARFd,"\n");
                        vn++;
                        len1 += n;
                        m += n;
                        j++;
                    }
                }
            }
        }
        printf1(ctx, "%d records (%d variables) written to: %s\n\n",vrec,vn,ctx->PMDVARFName);
    }
       
    /* if vsel */

    if (ctx->PMRHSTRA) {

        if (ctx->PMF1Def)  
            printf1(ctx, "Warning: vsel option will be ignored.\n");
        else {
            if ((n = v_parse(ctx, ctx->PMRHSTR + 5,0)) < 0 || ctx->ESCnt <= 0) {
                printf1(ctx, "Syntax error (%d) in vsel expression.\n",n);
                if (n < 0)
                    prn_emsg1(ctx, n);
                goto RSPFin;
            }
        }
    }
   
    /* reading data again to create internal data matrix, or directly
       writing to output file. */
    
    if (ctx->PMF1Def) {
        printf1(ctx, "Data will be directly written to output file: %s\n",ctx->PMF1dName);
        if (ctx->PMNOCFlg)
            mxnoc = ctx->PMNOC;
        else 
            mxnoc = ctx->INTMAX;
    }
    else {
        if (ctx->NOCMaxA <= 0)
            ctx->NOCMaxA = nrec;
        mxnoc = ctx->NOCMaxA;
       
        printf1(ctx, "Reading data again to create internal data matrix.\n");
        printf1(ctx, "Maximum number of cases: %d\n",ctx->NOCMaxA);
    }
    if (nmiss > 0)  
        printf1(ctx, "System missing values will be substituted by: %g\n",ctx->PMMSYS);

    ctx->RSPEOF = 0;                 /* seek to begin of data */
    fseek(ctx->PMFd,ctx->RSPFPtr1,0);
    ctx->RSPPtr = ctx->RSPBuf;
    ctx->RSPCnt = 0;
    if (buf_update(ctx, 0))
        goto RSPFin;
    ctx->RSPPtr += ctx->RSPFOff1;
    ctx->RSPCnt -= ctx->RSPFOff1;
    if (buf_update(ctx, 0))
        goto RSPFin;

    if (*ctx->RSPPtr != 'F') {               /* check for data stream */

        printf1(ctx, "Error: can't find data stream.\n");
        goto RSPFin;
    }
    ctx->RSPCnt--;
    ctx->RSPPtr++;

    /* allocate memory for data */
        
    if (ctx->PMF1Def == 0) {
        n = ctx->MemReq;
        if (alloc_vdat(ctx, idxn,1)) {
            printf1(ctx, "Error: insufficient memory for data matrix (%d cases).\n",ctx->NOCMaxA);
            goto RSPFin;
        }
        printf1(ctx, "Allocated %d bytes for data matrix.\n\n",ctx->MemReq - n);
    }
    nrec2 = nrec1 = 0;      /* number of records */
    nmiss = 0;              /* number of system missing values */

    while (*ctx->RSPPtr && *ctx->RSPPtr != 'Z') {    /* this is the end of file marker */

        if (nrec2 >= mxnoc)  
            break;
            
        nmiss1 = 0;
        for (i = 0; i < nvar; ++i) {
            if (buf_update(ctx, 0))
                goto RSPFin;

            if (vttyp[i] > 0) {        /* string variable */
                ctx->RSPPtr = inum(ctx, ctx->RSPPtr,&m);
                if (ctx->RSPPtr == NULL)
                    goto RSPFin;
                if (m > vttyp[i])  
                    printf1(ctx, "Warning: variable %d in record %d (type=%d fnd=%d).\n",
                                            i+1,nrec1+1,vttyp[i],m);
          
                p = ctx->VDPtr[i] - nrec1 * ctx->VSLen[i];                  
                for (k = 0; k < m; ++k) {

                    if (ctx->RSPCnt < 10) {
                        if (buf_update(ctx, 0))
                            goto RSPFin;
                    }
                    if (ctx->PMF1Def)
                        fprintf(ctx->PMF1d,"%c",*ctx->RSPPtr);
                    else
                        *p++ = *ctx->RSPPtr;

                    ctx->RSPPtr++;
                    ctx->RSPCnt--;
                }
                for (; k < vttyp[i]; ++k) {
                    if (ctx->PMF1Def)
                        fprintf(ctx->PMF1d," ");
                    else
                        *p++ = ' ';
                }
                if (ctx->PMF1Def)
                    fprintf(ctx->PMF1d," ");
            }
            else {
                ctx->RSPPtr = dnum(ctx, ctx->RSPPtr,&x,&mv);
                if (ctx->RSPPtr == NULL)
                    goto RSPFin;
    
                if (mv) {
                    nmiss1++;    
                    x = ctx->PMMSYS;
                }
                else if (vint[i] == 0) {
                    y = (int)x;
                    x = (int)y;
                }
                if (ctx->PMF1Def)
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->VPFmtS[i],x);
                else
                    put_data(ctx, x,i,nrec1);
     
            }
        }
        if (ctx->PMF1Def)
            fprintf(ctx->PMF1d,"\n");

        nrec2++;
        prn_message(ctx, nrec2,0,0);
        if (ctx->RSPEOF > RLEN)
            break;

        if (ctx->PMRHSTRA && ctx->PMF1Def == 0) {         /* vsel */

            n = v_eval1(ctx, nrec1,ctx->ESCnt,ctx->ESTyp,ctx->ESVal,ctx->ESIdx,&x,0,0,0,0,0);
            if (n) {
                printf1(ctx, "Can't evaluate vsel expression in record %d.\n",nrec2 + 1);
                prn_emsg2(ctx, n);
                goto RSPFin;
            }
            if (fabs(x) < ctx->EPSI2) {
                continue;     
            }
        }
        nrec1++;
        nmiss += nmiss1;
    }
    prn_message(ctx, nrec2,1,0);

    printf1(ctx, "Read %d records. ",nrec2);
    if (ctx->PMRHSTRA && ctx->PMF1Def == 0)
        printf1(ctx, "Selected: %d",nrec1);
    printf1(ctx, "\n");

    if (ctx->PMF1Def) {
        printf1(ctx, "%d records written to: %s\n",nrec1,ctx->PMF1dName);
        if (nrec1 <= 0)  
            goto RSPFin;
        if (len > 0)
            printf1(ctx, "Record length without EOL characters: %d bytes.\n",len);
    }
    else {
        if (nrec1 <= 0) {
            printf1(ctx, "Error: no data matrix created.\n");
            goto RSPFin;
        }
        ctx->NOCDM = ctx->NOC = nrec1;
        ctx->DMDef = 1;
       
        printf1(ctx, "Created a data matrix with %d variables and %d cases.\n",ctx->NVAR,ctx->NOC);
        if (nrec1 < nrec &&  ctx->PMRHSTRA == 0)
            printf1(ctx, "Warning: file contains more than %d records.\n",nrec1);
    }
    printf1(ctx, "Number of system missing values: %d\n",nmiss);

    if (aflag)          /* archive description file */
        make_arcd(ctx, ctx->PMDVARFN,ctx->PMF1dName,len,nrec1,vn);

    err = 0;

RSPFin:
    if (avttyp) {
        free((char *)vttyp);
        memrq(ctx, -avttyp,sizeof(short));
    }
    if (aptr) {
        free((char *)ctx->RSPVFPtr);
        memrq(ctx, -aptr,sizeof(long));
    }
    if (aoff) {
        free((char *)ctx->RSPVFOff);
        memrq(ctx, -aoff,sizeof(int));
    }
    if (avmax) {
        free((char *)vmax);
        memrq(ctx, -avmax,sizeof(double));
    }
    if (avmin) {
        free((char *)vmin);
        memrq(ctx, -avmin,sizeof(double));
    }
    if (avint) {
        free((char *)vint);
        memrq(ctx, -avint,sizeof(short));
    }
    if (ctx->RSPBufA > 0) {
        free(ctx->RSPBuf);
        memrq(ctx, -ctx->RSPBufA,1);
    }
    if (err || ctx->PMF1Def) {
        if (idxn < get_nidx(ctx))
            clear_avar(ctx, idxn);
        ctx->NVAR = ctx->NOCMaxA = ctx->NOCDM = ctx->NOC = 0;
    }
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  buf_update()    Read a new record into RSPBuf                           */
/*                  Return  0 if OK, -1 if error.                           */

int buf_update(TDAContext *ctx, int init)
{
    register int i;
    register char *p,*q;

    if (init)
        ctx->s_buf_update_fptr = 0L;

    if (ctx->RSPEOF || ctx->RSPCnt >= RLEN) {
        if (ctx->RSPEOF) {
            ctx->s_buf_update_fptr = 0L;
            ctx->RSPEOF++;
        }
        if (ctx->RSPCnt == RLEN) {
            ctx->RSPFPtr = ctx->s_buf_update_fptr;
            ctx->RSPFOff = 0;
        }
        else
            ctx->RSPFOff = 2 * RLEN - ctx->RSPCnt;
        return(0);
    }   
    ctx->RSPFOff = RLEN - ctx->RSPCnt;
    ctx->RSPFPtr = ctx->s_buf_update_fptr;

    ctx->s_buf_update_fptr = ftell(ctx->PMFd);

    p = ctx->RSPBuf;
    q = ctx->RSPPtr;
    for (i = 0; i < ctx->RSPCnt; ++i)
        *p++ = *q++;
    ctx->RSPPtr = ctx->RSPBuf;

    if (ctx->RSRILen == 0) {
        if (fgets(p,RLEN + 3,ctx->PMFd)) {
            q = p;
            for (i = 0; i < RLEN; ++i) {
                if (!*q || *q == CR || *q == LF)
                    break;
                q++;
            }
            for (; i < RLEN; ++i)  
                *q++ = ' ';
            *q = '\0';
            ctx->RSPCnt += RLEN;
        }
        else  
            ctx->RSPEOF = 1;
    }
    else {
        if (((i = (int)fread(p,sizeof(*p),(size_t)(ctx->RSRILen),ctx->PMFd))) < 0) {
            printf1(ctx, "\nError in reading the input file.\n");
            return(-1);
        }
        if (i < RLEN) {
            ctx->RSPCnt += i;
            ctx->RSPEOF = 1;
        }
        else  
            ctx->RSPCnt += RLEN;
    }
    if (ctx->RSPEOF)    
        *(ctx->RSPPtr + ctx->RSPCnt) = '\0';
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  prstr(p,n)      print a string at p with length n                       */

void prstr(TDAContext *ctx, char *p, int n)
{
    register int i;
    for (i = 0; i < n; ++i) {
        if (*p >= ' ' && *p <= 'z')
            printf1(ctx, "%c",*p);
        else
            printf1(ctx, "?");
        p++;
    }
}

/*--------------------------------------------------------------------------*/
/*  getstr(f)   Get a string, starting at RSPPtr. Return pointer to next    */
/*              character after the string.                                 */
/*              if f >= 0, the string is printed to stdout and              */
/*                             filled with up to f blanks.                  */
/*                                                                          */
/*                  If error, return NULL.                                  */

char *getstr(TDAContext *ctx, int f)
{           
    register int i;
    int n;

    ctx->RSPPtr = inum(ctx, ctx->RSPPtr,&n);
    if (ctx->RSPPtr == NULL)  
        return(NULL);

    if (n > 0) {

        if (ctx->RSPCnt <= n) {
            if (buf_update(ctx, 0))
                return(NULL);     
        }
        if (ctx->RSPCnt < n) {
            printf1(ctx, "Error in reading SPSS file.\n");
            return(NULL);
        }
        if (f >= 0) {
            for (i = 0; i < n; ++i)  
                printf1(ctx, "%c",*ctx->RSPPtr++);
            for (; i < f; ++i)  
                printf1(ctx, " ");
        }
        else
            ctx->RSPPtr += n;

        ctx->RSPCnt -= n;
    }
    return(ctx->RSPPtr);
}

/*--------------------------------------------------------------------------*/
/*  getdigit(p,err)      get and return a base30 digit at p.                */
/*                       in case of an error, return err = 1.               */

int getdigit(TDAContext *ctx, char *p, int *err)
{
    (void)ctx;        /* unused: the signature is shared */
    *err = 0;
    if (*p >= '0' && *p <= '9')
        return((int)(*p - '0'));
    else if (*p >= 'A' && *p <= 'T')
        return((int)(10 + *p - 'A'));
    else  
        *err = 1;
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  inum(p,n)       read an integer value at p and return it in n, return   */
/*                  pointer to next character. If error, return NULL        */

char *inum(TDAContext *ctx, char *p, int *n)
{
    int err;     
    register int ex = 0;
    register int nn = 0;
    register int neg = 0;
    register int len = 0;
    register char *q;
    double tmp;

    q = p;
    if (*p == '-') {
        neg = 1;
        ctx->RSPCnt--;
        p++;
    }
    *n = 0;
    while (*p && *p != '/') {
        if (*p == '+')
            ex = 1;
        else if (*p == '-')
            ex = -1;
        else {
            if (!ex) {
                *n *= 30;
                *n += getdigit(ctx, p,&err);
                if (err) {
                    digit_err(ctx, q);
                    return(NULL);
                }
                len++;
                if (len > 8)  
                    printf1(ctx, "Warning: found integer entry with %2d (base-30) digits.\n",len);
            }
            else {
                nn *= 30;
                nn += getdigit(ctx, p,&err);
                if (err) {
                    digit_err(ctx, q);
                    return(NULL);
                }
            }
        }
        p++;
        ctx->RSPCnt--;
    }
    if (neg)
        *n = -(*n);

    if (ex) {
        tmp = (double)*n;
        if (ex == 1)  
            tmp *= pow(30.0,(double)nn);
        else if (ex == -1)
            tmp /= pow(30.0,(double)nn);
        *n = (int) tmp;
    }
    ctx->RSPCnt--;
    return(++p);       
}

/*--------------------------------------------------------------------------*/
/*  dnum(p,x,mv)    read a double value at p and return it in x, return     */
/*                  pointer to next character. mv is set to 1 if a system   */
/*                  missing value code is found. If error, return NULL      */

char *dnum(TDAContext *ctx, char *p, double *x, int *mv)
{
    int err;
    register int neg = 0;
    register int pnt = 0;
    register int ex  = 0;
    register int k = 0;
    register int n = 0;
    double man = 0.0;
    double mex = 0.0;
    register char *q;

    q = p;
    while (*p && *p == ' ') {
        ctx->RSPCnt--;
        p++;
    }
    *mv = 0;
    if (*p == '*') {                /* check for internal missing value */
        *mv = 1;
        ctx->RSPCnt -= 2;
        p += 2;
        return(p);
    }
    if (*p == '-') {
        neg = 1;
        ctx->RSPCnt--;
        p++;
    }
    while (*p && *p != '/') {

        if (*p == '.') {
            pnt = 1;
            n = k;
        }
        else if (*p == '+')
            ex = 1;
        else if (*p == '-')  
            ex = -1;
        else {
            if (!ex) {
                man *= 30.0;
                man += (double) getdigit(ctx, p,&err);
                if (err) {
                    digit_err(ctx, q);
                    return(NULL);
                }
                k++;
                if (k > 13)  
                    printf1(ctx, "Warning: found entry with %2d (base-30) digits.\n",k);
            }
            else {
                mex *= 30.0;
                mex += (double) getdigit(ctx, p,&err);
                if (err) {
                    digit_err(ctx, q);
                    return(NULL);
                }
            }
        }
        ctx->RSPCnt--;
        p++;
    }
    if (neg)
        man = -man;

    if (pnt) {
        k -= n;
        while (k--)
            man /= 30.0;
    }
    if (ex == 1)  
        man *= pow(30.0,mex);
    else if (ex == -1)
        man /= pow(30.0,mex);
    *x = man;
    ctx->RSPCnt--;
    return(++p);
}

/*--------------------------------------------------------------------------*/
/*  digit_err(p)        print error message                                 */

void digit_err(TDAContext *ctx, char *p)
{
    register int i;
    printf1(ctx, "\nError: can't read a numerical value.");
    printf1(ctx, "\nFound: ");
    for (i = 0; i < 40; ++i)  
        printf1(ctx, "%c",*p++);
    printf1(ctx, "\n");
}

/*--------------------------------------------------------------------------*/
/*  wr_spss         Write SPSS export file.                                 */
/*                                                                          */
/*      wspss(                                                              */
/*              keep = varlist,                                             */
/*              drop = varlist,                                             */
/*              sort = varlist,                                             */
/*      ) = output-file-name;                                               */
/*                                                                          */
/*      Return 0 if OK, -1 if error.                                        */

char S1[] =
"00000@0000@0000@0000@@@@@@@@@@@@@@@@@@@@ASCII SPSS PORT FILE                    ";
char S2[] =
"00000-0000-0000-0000--------------------!3#))0303300/240&),%00000000000000000000";
char S3[] =
"0200002'220'&)3000#0000000000000000000000000000000000000000000000000000000000000";
char S4[] =
"0000000000000000000000000123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrst";
char S5[] =
"uvwxyz .<(+0&[]!$*);^-/|,%_>?`:#@'=\"000000~000000000000000000000{}\\0000000000000";
char S6[] =
"00000000000000000000000000000000000000000000000000000000SPSSPORTA8/190000006/000";
char S7[] =
"0001F/TDA Version 6.x";


int wr_spss(TDAContext *ctx)
{
    register int i,j,k,ii;
    int err,u,v,nv,sflag,l;
    double x;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 5,1,1)) {   /* get parameters */
        p_clean(ctx);
        return(-1);   
    }       
    printf1(ctx, "Writing SPSS export file: %s\n",ctx->PMFdName);

    if (ctx->PMKeep && ctx->PMDrop) {
        p_err(ctx, -16,1);
        goto WSPSSFin;
    }

    /* get list of variables in AcI[] */

    if (alloc_aci(ctx, imax(ctx, ctx->NVAR,ctx->PMNV)))
        goto WSPSSFin;

    nv = 0;
    if (ctx->PMKeep) {
        for (k = 0; k < ctx->PMNV; ++k)  
            ctx->AcI[nv++] = ctx->PMVIdx[k];
    }
    else {
        j = ctx->VIFirst;
        while (j >= 0) {
            i = 1;
            if (ctx->PMDrop) {
                for (k = 0; k < ctx->PMNV; ++k) {
                    if (ctx->PMVIdx[k] == j) {
                        i = 0;
                        break;
                    }   
                }
            }
            if (i) 
                ctx->AcI[nv++] = j;
            j = ctx->VNxt[j];
        }
    }
    if (nv == 0) {
        printf1(ctx, "No variables selected.\n"); 
        goto WSPSSFin;
    }
    sflag = 0;
    if (ctx->PM1NV > 0) {        /* sort */
        err = vsort(ctx, ctx->PM1NV,ctx->PM1VIdx,1,0,1);
        if (err)
            goto WSPSSFin;
        sflag = 1;
    }

    /* print header etc. */

    fprintf(ctx->PMFd,"%s\n",ctx->S1);
    fprintf(ctx->PMFd,"%s\n",ctx->S2);
    fprintf(ctx->PMFd,"%s\n",ctx->S3);
    fprintf(ctx->PMFd,"%s\n",ctx->S4);
    fprintf(ctx->PMFd,"%s\n",ctx->S5);
    fprintf(ctx->PMFd,"%s\n",ctx->S6);
    fprintf(ctx->PMFd,"%s",ctx->S7);

    /* print number of variables */

    fprintf(ctx->PMFd,"4");
    check_nl(ctx, ctx->PMFd);
    pnum(ctx, ctx->PMFd,nv);

    /* digits (optional) */

    fprintf(ctx->PMFd,"5");
    check_nl(ctx, ctx->PMFd);
    pnum(ctx, ctx->PMFd,10);

    /* case weight variable name (optional) not written */

    /* dictionary entries (required) */

    for (k = 0; k < nv; ++k) {
        j = ctx->AcI[k];
        fprintf(ctx->PMFd,"7");
        check_nl(ctx, ctx->PMFd);
        if (ctx->VTyp[j] == 1)
            pnum(ctx, ctx->PMFd,-ctx->VSLen[j]);
        else
            pnum(ctx, ctx->PMFd,0);
        pstring(ctx, ctx->PMFd,ctx->VName[j],8);
        u = (int)ctx->VPFmt1[j];
        v = (int)ctx->VPFmt2[j];
        if (u == 0) {
            u = 10;
            v = 4;
        }
        else if (u < 0) {
            u = 12;
            v = 4;
        }
        pnum(ctx, ctx->PMFd,5);
        pnum(ctx, ctx->PMFd,u);
        pnum(ctx, ctx->PMFd,v);
        pnum(ctx, ctx->PMFd,5);
        pnum(ctx, ctx->PMFd,u);
        pnum(ctx, ctx->PMFd,v);

        if (ctx->VLabelLen > 0) {
            if (ctx->VLabel[j] != NULL) {
                fprintf(ctx->PMFd,"C");
                check_nl(ctx, ctx->PMFd);
                pstring(ctx, ctx->PMFd,ctx->VLabel[j],40);
            }
        }
    }

    /* print data stream */

    if (ctx->SILENTFlg < 2)
        printfe(ctx, "Writing: %s\n",ctx->PMFdName);
   
    fprintf(ctx->PMFd,"F");
    check_nl(ctx, ctx->PMFd);

    for (i = 0; i < ctx->NOC; ++i) {
        ii = i;
        if (sflag)
            ii = ctx->VSORTPtr[i];

        for (k = 0; k < nv; ++k) {
            j = ctx->AcI[k];

            if (ctx->VTyp[j] == 1) {
                l = -ctx->VSLen[j];
                pstring(ctx, ctx->PMFd,ctx->VDPtr[j] + ii * l,l);
            }
            else {
                x = get_data(ctx, j,ii);           
                if (x < 0.0) {
                    fprintf(ctx->PMFd,"-");
                    check_nl(ctx, ctx->PMFd);
                    x = -x;
                }
                u = (int)x;
                if (x == (double)u)
                    pnum(ctx, ctx->PMFd,u);
                else  
                    pfnum(ctx, ctx->PMFd,x);
            }
        }
        prn_message(ctx, i + 1,0,1);
    }
    prn_message(ctx, ctx->NOC,1,1);
    fprintf(ctx->PMFd,"Z");            /* end of file */
    check_nl(ctx, ctx->PMFd);
    while (ctx->SPSSPtr++ < RLEN)  
        fprintf(ctx->PMFd,"Z");        
    fprintf(ctx->PMFd,"\n");

    printf1(ctx, "%d records with %d variables written to: %s\n",ctx->NOC,nv,ctx->PMFdName);
    err = 0;

WSPSSFin:
    if (ctx->PM1NV > 0)                     
        vsort(ctx, 0,ctx->PM1VIdx,0,0,1);
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pnum1(fd,n)     print integer n to fd (base 30). Update SPSSPtr.        */


void pnum1(TDAContext *ctx, FILE *fd,int n)
{
    int m,r;
    register char *p;
    char buf[100];

    if (n < 0) {
        fprintf(fd,"-");
        check_nl(ctx, fd);
        n = -n;
    }
    p = buf;
    while (n >= 30) {
        m = (int) (n / 30);
        r = n - 30 * m;
        *p++ = ctx->DIG30[r];
        n = m;
    }
    *p++ = ctx->DIG30[n];
    while (p > buf) {
        fprintf(fd,"%c",*--p);
        check_nl(ctx, fd);
    }
}

/*--------------------------------------------------------------------------*/
/*  pnum(fd,n)      print integer n to fd (base 30). Update SPSSPtr.        */

void pnum(TDAContext *ctx, FILE *fd,int n)
{
    pnum1(ctx, fd,n);
    fprintf(fd,"/");
    check_nl(ctx, fd);
}

/*--------------------------------------------------------------------------*/
/*  pfnum(fd,x)     print floating point number x to fd (base 30).          */
/*                  Update SPSSPtr. Note: x >= 0.0                          */

void pfnum(TDAContext *ctx, FILE *fd,double x)
{
    register int i;
    double a,b,c,d,e;

    e = floor(log(x) / log(30.0));
    b = x / pow(30.0,e);
    c = floor(b);
    if (c < 0.0 || c >= 30.0)  
        gerr_exit(ctx, 74);
         
    fprintf(fd,"%c",ctx->DIG30[(int)c]);
    check_nl(ctx, fd);
    b -= c;
    if (b > ctx->EPSI) {
        fprintf(fd,".");
        check_nl(ctx, fd);
        c = 30.0;
        for (i = 0; i < 10 ; ++i) {
            a = b * c;
            d = floor(a);
            fprintf(fd,"%c",ctx->DIG30[(int)d]);
            check_nl(ctx, fd);
            b -= d / c;
            if (b <= ctx->EPSI)  
                break;
            c *= 30;
        }
    }
    i = (int)e;
    if (i) {
        if (i < 0) {
            fprintf(fd,"-");
            check_nl(ctx, fd);
            i = -i;
        }
        else {
            fprintf(fd,"+");
            check_nl(ctx, fd);
        }
        pnum1(ctx, fd,i);
    }
    fprintf(fd,"/");
    check_nl(ctx, fd);
}

/*--------------------------------------------------------------------------*/
/*  pstring(fd,s,m)     print string s to fd. If m > 0, this is max length  */

void pstring(TDAContext *ctx, FILE *fd,char *s,int m)
{
    register int l;
    register char *p;

    p = s;
    l = (int)(strlen(p));
    if (m > 0 && l > m)  
        l = m;
    pnum(ctx, fd,l);
    while (l-- > 0) {
        fprintf(fd,"%c",*p++);
        check_nl(ctx, fd);
    }
}

/*--------------------------------------------------------------------------*/
/*  check_nl(fd)    check SPSSPtr. If >= RLEN print newline and set         */
/*                  SPSSPtr = 0;                                            */

void check_nl(TDAContext *ctx, FILE *fd)
{
    if (++ctx->SPSSPtr < 80)
        return;
    fprintf(fd,"\n");
    ctx->SPSSPtr = 0;
}

/* -##--------------------------------------------------------------------- */
/*  rd_spss1    Read SPSS sav file and create an internal data matrix.      */
/*                                                                          */
/*              rspss1(                                                     */
/*                  noc=...,            # of records to read, def. all      */
/*                  dvar=...,           file containing value label info    */
/*                  msys=...,           new sys miss value code, def. -5    */  
/*                  df=...,             write data directly to output file  */
/*              ) = file_name;                                              */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

#define SPBUFLEN 260        /* max length of read buffer */

int rd_spss1(TDAContext *ctx)  
{
    register int i,j;
    register char *p;
    int err,idxn,l,ln,n,n1,r,nr,nobs,cflag,widx,noc,nvar,nvar1,vlflag,nmiss,vtyp,rtyp;
    int slen,slen1,jlen,ii,jo,mxnoc,dvn,eof,nrec,nsys1,nsys2,pfmt1,pfmt2,sval[9];
    int arch;
    char buf[SPBUFLEN + 1],vname[VNLMax + 1],vlabel[VLLMax + 1],vdef[400];
    double a,x,y,bias,tol,smiss;
    long int fptr;
         
    arch = 1;
    slen = dvn = 0;
    tol = 10000.0 * ctx->EPSI;
    smiss = -(ctx->DBLMAX - 100.0);

    idxn = get_nidx(ctx);      /* index to first new variable */
    ctx->VLabelLen = 0;
    err = -1;

    if (check_cmd(ctx, 0))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 6,10,1)) {   /* get parameters */
        p_clean(ctx);
        return(-1);   
    }   
    printf1(ctx, "Reading SPSS sav file: %s\n",ctx->PMFdName);
    if (ctx->DMDef) {
        printf1(ctx, "Error: a data matrix already exists.\n");
        p_clean(ctx);
        return(-1);
    }
    if (ctx->PMNOCFlg)
        printf1(ctx, "Maximum number of cases: %d\n",ctx->PMNOC);
    newline(ctx);

    if ((r = rd_dat(ctx, 64,buf)) != 64)
        goto RSP1Fin; 

    printf1(ctx, "Identification: %s\n",buf);
    if (strncmp(buf,"$FL2",4)) {
        printf1(ctx, "Probably not an SPSS sav file.\n");
        err = -1;
        goto RSP1Fin;
    }
    n = sav_geti(ctx, arch,&r);
    if (r)
        goto RSP1Fin;
        
    if (n != 2)         /* file layout code should be 2 */           
        arch = 2;                  

    nobs = sav_geti(ctx, arch,&r);
    if (r)
        goto RSP1Fin;
    printf1(ctx, "Number of OBS elements per observation: %d\n",nobs);

    cflag = sav_geti(ctx, arch,&r);
    if (r)
        goto RSP1Fin;
    printf1(ctx, "Compression switch: %d\n",cflag);

    widx = sav_geti(ctx, arch,&r);
    if (r)
        goto RSP1Fin;
    printf1(ctx, "Index of case-weight variable: %d\n",widx);

    noc = sav_geti(ctx, arch,&r);
    if (r)
        goto RSP1Fin;
    printf1(ctx, "Number of cases: %d\n",noc);

    bias = sav_getd(ctx, arch,&r);
    if (r)
        goto RSP1Fin;
    printf1(ctx, "Compression bias: %g\n\n",bias);

    if ((r = rd_dat(ctx, 9,buf)) != 9)
        goto RSP1Fin; 
    printf1(ctx, "Creation date: %s\n",buf);

    if ((r = rd_dat(ctx, 8,buf)) != 8)
        goto RSP1Fin; 
    printf1(ctx, "Creation time: %s\n",buf);

    if ((r = rd_dat(ctx, 67,buf)) != 67)
        goto RSP1Fin; 
    printf1(ctx, "File label: %s\n\n",buf);

    if (nobs < 1) {
        printf1(ctx, "Error: need at least one variable description.\n");
        goto RSP1Fin;
    }
    if (alloc_acn(ctx, nobs + 1))        /* type of variable */
        goto RSP1Fin; 

    if (alloc_acm(ctx, nobs + 1))        /* spss-internal variable number */
        goto RSP1Fin; 

    if (alloc_aci(ctx, nobs + 1))        /* integer flags */
        goto RSP1Fin; 

    nvar1 = nvar = 0;
    for (i = 0; i < nobs; ++i) {

        n = sav_geti(ctx, arch,&r);
        if (r)
            goto RSP1Fin;
   
        if (n != 2) {
            printf1(ctx, "Error: can't read next dictionary entry.\n");
            goto RSP1Fin;
        }
        vtyp = sav_geti(ctx, arch,&r);
        if (r)
            goto RSP1Fin;

        ctx->AcN[i] = vtyp;

        if (vtyp >= 0)
            ctx->AcM[i] = nvar;
        else
            ctx->AcM[i] = -1;

        vlflag = sav_geti(ctx, arch,&r);
        if (r)
            goto RSP1Fin;

        nmiss = sav_geti(ctx, arch,&r);
        if (r)
            goto RSP1Fin;

        n = sav_geti(ctx, arch,&r);
        if (r)
            goto RSP1Fin;

        n = sav_geti(ctx, arch,&r);
        if (r)
            goto RSP1Fin;

        if ((r = rd_dat(ctx, 8,vname)) != 8) {
            printf1(ctx, "Error: can't read variable name.\n");
            goto RSP1Fin; 
        }
        p = vname;
        for (j = 0; j < VNLMax; ++j) {
            if (*p == ' ')
                break;
            p++;
        }
        *p = '\0';
        if (vlflag) {
            n1 = n = sav_geti(ctx, arch,&r);
            if (r)
                goto RSP1Fin;

            if (n > 0) {
                n = (n - 1) / 4 + 1;
                n *= 4;
                if (n > SPBUFLEN) {
                    printf1(ctx, "Error: exceeded max var label length.\n");
                    goto RSP1Fin;
                }
  
                if ((r = rd_dat(ctx, n,buf)) != n) {
                    printf1(ctx, "Error: can't read variable label.\n");
                    goto RSP1Fin; 
                }
                strncpy(vlabel,buf,VLLMax);
                vlabel[VLLMax] = '\0';
                vlabel[imin(ctx, VLLMax,n1)] = '\0';
   
                p = vlabel;
                while (*p) {
                    if (*p == '(')
                    *p = '[';
                    else if (*p == ')')
                        *p = ']';
                    p++;
                }
                snprintf(vdef,sizeof(vdef),"%s<4>[0.0](%s)=spss(%d)",vname,vlabel,vtyp);
            }
            else
                snprintf(vdef,sizeof(vdef),"%s<4>[0.0]=spss(%d)",vname,vtyp);
        }
        else
            snprintf(vdef,sizeof(vdef),"%s<4>[0.0]=spss(%d)",vname,vtyp);
    
        for (j = 0; j < iabs(ctx, nmiss); ++j) {
            x = sav_getd(ctx, arch,&r);
            if (r)
                goto RSP1Fin;
        }
        if (vtyp >= 0) {
            if (save_var(ctx, vdef,0)) {       /* save variable definition */

                printf1(ctx, "Error: can't save variable definitions.\n");
                goto RSP1Fin;
            }
            nvar++;
            if (vtyp > 0) {
                slen = imax(ctx, slen,vtyp);
                nvar1++;
            }
        }
    }
    printf1(ctx, "Number of variables: %d\n",nvar);
    printf1(ctx, "Number of string variables: %d\n\n",nvar1);

RSP1Cont:
    rtyp = sav_geti(ctx, arch,&r);
    if (r)
        goto RSP1Fin;

    if (rtyp == 3) {            /* value labels */
         
        if (dvn == 0 && ctx->PMDVARFDef) {
            fprintf(ctx->PMDVARFd,"Value labels found in: %s\n\n",ctx->PMFdName);
            dvn += 2;
        }
        n = sav_geti(ctx, arch,&r);       /* number of labels */
        if (r)
            goto RSP1Fin;

        for (i = 0; i < n; ++i) {

            x = sav_getd(ctx, arch,&r);
            if (r)
                goto RSP1Fin;

            if (ctx->PMDVARFDef) {
                l = (int)x;
                if ((double)l == x)
                    fprintf(ctx->PMDVARFd,"%12d  ",l);
                else
                    fprintf(ctx->PMDVARFd,"%12g  ",x);
            }
            if ((r = rd_dat(ctx, 1,buf)) != 1) {
                printf1(ctx, "Error in reading value labels.\n");
                goto RSP1Fin; 
            }
            ln = l = (int)buf[0];
            l = l / 8 + 1;
            l *= 8;

            if ((r = rd_dat(ctx, l - 1,buf)) != l - 1) {
                printf1(ctx, "Error in reading value labels.\n");
                goto RSP1Fin; 
            }
            if (ctx->PMDVARFDef) {
                *(buf + ln) = '\0';
                fprintf(ctx->PMDVARFd,"%s (%d)\n",buf,ln);
                dvn++;
            }
        }

        rtyp = sav_geti(ctx, arch,&r);
        if (r)
            goto RSP1Fin;

        if (rtyp != 4) {
            printf1(ctx, "Error: in reading value label entries.\n");
            goto RSP1Fin;
        }
        n = sav_geti(ctx, arch,&r);       /* number of variables */
        if (r)
            goto RSP1Fin;

        for (i = 0; i < n; ++i) {

            l = sav_geti(ctx, arch,&r);       /* number of variable */
            if (r)
                goto RSP1Fin;

            l--;
            if (l >= 0 && l < nobs) {
                ii = ctx->AcM[l];
                if (ctx->PMDVARFDef && ii >= 0 && ii < nvar)  
                    fprintf(ctx->PMDVARFd,"%s ",ctx->VName[ii]);
            }
        }
        if (ctx->PMDVARFDef) {
            fprintf(ctx->PMDVARFd,"\n");
            dvn++;
            fprnchar(ctx, ctx->PMDVARFd,'-',70,1);
            dvn++;
        }
        goto RSP1Cont;
    }
    if (ctx->PMDVARFDef) {
        printf1(ctx, "Value labels: ");
        if (dvn == 0)
            printf1(ctx, "not present.\n");
        else                
            printf1(ctx, "%d records written to: %s\n",dvn,ctx->PMDVARFName);
    }
    if (rtyp == 6) {        /* document record */

        printf1(ctx, "Reading document records (will be ignored).\n");

        n = sav_geti(ctx, arch,&r);       /* number of records */
        if (r)
            goto RSP1Fin;

        for (i = 0; i < n; ++i) {
            if ((r = rd_dat(ctx, 80,buf)) != 80)
                goto RSP1Fin; 
        }
        rtyp = sav_geti(ctx, arch,&r);
        if (r)
            goto RSP1Fin;
    }
    while (rtyp == 7) {        /* type 7 records */

        n = sav_geti(ctx, arch,&r);       /* subtype */
        if (r)
            goto RSP1Fin;

        l = sav_geti(ctx, arch,&r);       /* length */
        if (r)
            goto RSP1Fin;

        n = sav_geti(ctx, arch,&r);       /* number of elements */
        if (r)
            goto RSP1Fin;

        for (i = 0; i < n; ++i) {
            if ((r = rd_dat(ctx, l,buf)) != l)
                goto RSP1Fin; 
        }
        rtyp = sav_geti(ctx, arch,&r);
        if (r)
            goto RSP1Fin;

    }
    if (rtyp != 999) {
        printf1(ctx, "Error: can't find data stream.\n");
        goto RSP1Fin;
    }
    fptr = ftell(ctx->PMFd) - 4L;        /* save position to data stream */

    n = sav_geti(ctx, arch,&r);
    if (r)
        goto RSP1Fin;

    printf1(ctx, "Reading data to check variables.\n");

    if (ctx->SILENTFlg < 2)
        printfe(ctx, "Reading: %s\n",ctx->PMFdName);

    if (alloc_acx(ctx, nobs + 1)) 
        goto RSP1Fin; 
    if (alloc_acy(ctx, nobs + 1)) 
        goto RSP1Fin; 

    for (j = 0; j < nobs; ++j) {
        ctx->AcX[j] = ctx->DBLMAX;
        ctx->AcY[j] = ctx->DBLMIN;
    }
    eof = nrec = nsys1 = nsys2 = 0;

    if (cflag) {                /* if compressed data */

        j = 0;
        while (eof == 0) {
    
            nr = sav_sval(ctx, sval);         
            if (nr == 0) {
                eof = 1;
                break;
            }
            else if (nr < 0) {
                printf1(ctx, "Error in reading data stream.\n");
                goto RSP1Fin;
            }
            for (i = 0; i < nr; ++i) {
                if (sval[i] == 252) {        /* end of file */
                    eof = 1;
                    break;
                }
                else if (sval[i] == 253) {      /* uncompressed */
                    if (ctx->AcN[j] == 0) {
                        x = sav_getd(ctx, arch,&r);
                        if (r)
                            goto RSP1Fin;
                    }
                    else {
                        if ((r = rd_dat(ctx, 8,buf)) != 8)
                            goto RSP1Fin; 
                        x = 0.0;
                    }
                    ctx->AcX[j] = dmin(ctx, x,ctx->AcX[j]);
                    ctx->AcY[j] = dmax(ctx, x,ctx->AcY[j]);
                }
                else if (sval[i] == 254) {
                    x = ctx->MBlnkVal;
                    nsys1++;
                }
                else if (sval[i] == 255) {
                    x = ctx->PMMSYS;
                    nsys2++;
                }
                else  
                    x = (double)sval[i];

                ctx->AcX[j] = dmin(ctx, x,ctx->AcX[j]);
                ctx->AcY[j] = dmax(ctx, x,ctx->AcY[j]);
                y = floor(x + 0.5);
                if (fabs(y - x) > tol)      
                    ctx->AcI[j] = 1;

                if (++j >= nobs) {
                    j = 0;
                    nrec++;
                    prn_message(ctx, nrec,0,0);
   
                    if (ctx->PMNOCFlg && nrec >= ctx->PMNOC)
                        goto RSP1Nxt;
                }
            }
        }
    }
    else {                  /* uncompressed data */
       
        while (1) {

            for (j = 0; j < nobs; ++j) {

                x = sav_getd(ctx, arch,&r);
                if (r) {
                    if (j == 0)
                        goto RSP1Nxt;
                    goto RSP1Fin;
                }
                if (x <= smiss) {
                    nsys2++;
                    x = ctx->PMMSYS;
                }
                if (ctx->AcN[j] != 0)
                    continue;

                ctx->AcX[j] = dmin(ctx, x,ctx->AcX[j]);
                ctx->AcY[j] = dmax(ctx, x,ctx->AcY[j]);
                y = floor(x + 0.5);
                if (fabs(y - x) > tol)      
                    ctx->AcI[j] = 1;
            }
            j = 0;
            nrec++;
            prn_message(ctx, nrec,0,0);
   
            if (ctx->PMNOCFlg && nrec >= ctx->PMNOC)
                goto RSP1Nxt;
        }
    }

RSP1Nxt:
    prn_message(ctx, nrec,1,0);

    printf1(ctx, "Read %d records.\n",nrec);

    if (nrec <= 0) {
        printf1(ctx, "Error.\n");
        goto RSP1Fin;
    }
    if (j) {
        printf1(ctx, "Error: data stream ends with incomplete record.\n");
        goto RSP1Fin;
    }
    printf1(ctx, "Number of blank-type missing values: %d\n",nsys1);
    printf1(ctx, "Number of system-type missing values: %d\n\n",nsys2);

/*** 
for (i = 0; i < nobs; ++i)
    printf1(ctx, "i=%4d  %f %f aci=%4d acn=%4d   \n",i,AcX[i],AcY[i],AcI[i],AcN[i]);
***/            

    /* set appropriate storage size and write formats */

    ii = 0;
    for (i = 0; i < nobs; ++i) {

        if (ctx->AcN[i] < 0)
            continue;

        if (ctx->AcN[i] == 0) {
            x = fabs(ctx->AcX[i]);
            y = fabs(ctx->AcY[i]);
            if (x < y)
                a = y;
            else
                a = x;
    
            if (ctx->AcI[i])
                ctx->VSLen[ii] = 8;
            else if ((int)a < 128)
                ctx->VSLen[ii] = 1;
            else if ((int)a < 32000)
                ctx->VSLen[ii] = 2;
            else 
                ctx->VSLen[ii] = 5;

            x = ctx->AcX[i];
            y = ctx->AcY[i];

            if (ctx->AcI[i]) {
                x = ceil(x);
                y = floor(y);
            }
            pfmt2 = 0;
            a = 10.0;
            for (j = 2; j <= 12; ++j) {
                if (-a < x && y < a * 10.0) {
                    pfmt1 = j;
                    break;
                }
                a *= 10.0;
            }
            if (ctx->AcI[i]) {
                if (ctx->PMFmtF) {
                    pfmt1 = ctx->PMFmt1;
                    pfmt2 = ctx->PMFmt2;
                }
                else {
                    pfmt1 += 7;
                    pfmt2 = 6;
                }
            }
            makefmt(ctx, &pfmt1,&pfmt2,ctx->VPFmtS[ii],VPFmtSLen,0,ctx->SEPC,0);
            ctx->VPFmt1[ii] = (short)pfmt1;
            ctx->VPFmt2[ii] = (short)pfmt2;
        }
        ii++;
    }
    prn_var(ctx, idxn);    /* print list of new variables */
    newline(ctx);     
   
    if (ctx->PMF1Def) {
        printf1(ctx, "Data will be directly written to output file: %s\n",ctx->PMF1dName);
        if (ctx->PMNOCFlg)
            mxnoc = ctx->PMNOC;
        else 
            mxnoc = ctx->INTMAX;
    }
    else {
        if (ctx->NOCMaxA <= 0)
            ctx->NOCMaxA = nrec;
        mxnoc = imin(ctx, nrec,ctx->NOCMaxA);
       
        printf1(ctx, "Reading data again to create internal data matrix.\n");
        printf1(ctx, "Maximum number of cases: %d\n",ctx->NOCMaxA);
    }   
    fseek(ctx->PMFd,fptr,0);

    rtyp = sav_geti(ctx, arch,&r);
    if (r)
        goto RSP1Fin;

    if (rtyp != 999) {
        printf1(ctx, "Error: can't seek to data stream.\n");
        goto RSP1Fin;
    }
    n = sav_geti(ctx, arch,&r);
    if (r)
        goto RSP1Fin;
         
    /* allocate memory for data */
        
    if (ctx->PMF1Def == 0) {
        n = ctx->MemReq;
        if (alloc_vdat(ctx, idxn,1)) {
            printf1(ctx, "Error: insufficient memory for data matrix (%d cases).\n",ctx->NOCMaxA);
            goto RSP1Fin;
        }
        printf1(ctx, "Allocated %d bytes for data matrix.\n\n",ctx->MemReq - n);
    }
    slen1 = (slen - 1) / 8 + 1;
    slen1 *= 8;
    if (alloc_acc(ctx, slen1 + 10))        /* buffer for string variables */
        goto RSP1Fin; 

    eof = nrec = 0;
/* ### */
    if (cflag) {                /* if compressed data */

        ii = jlen = jo = j = 0;

        while (eof == 0) {  

            if (nrec >= mxnoc)  
                break;
            
            nr = sav_sval(ctx, sval);         

            if (nr == 0) {
                eof = 1;
                break;
            }   
            else if (nr < 0 || nr > 8) {
                printf1(ctx, "Error in reading data stream.\n");
                goto RSP1Fin;
            }

            for (i = 0; i < nr; ++i) {
                if (sval[i] == 252) {       /* end of file */
                    eof = 1;
                    break;
                }
                else if (sval[i] == 253) {  /* uncompressed */
                    if (ctx->AcN[j] == 0) {
                        x = sav_getd(ctx, arch,&r);
                        if (r)
                            goto RSP1Fin;
                        if (ctx->PMF1Def)
                            rt_fprintf_d(ctx, ctx->PMF1d,ctx->VPFmtS[ii],x);
                        else
                            put_data(ctx, x,ii,nrec);
                        ii++;
                    }
                    else {                
                        if (ctx->AcN[j] > 0)  
                            jlen = ctx->AcN[j];
                        if ((r = rd_dat(ctx, 8,ctx->AcC + jo)) != 8)
                            goto RSP1Fin; 
                        jo += 8;

                        if (jo >= jlen) {
                            if (ctx->PMF1Def) {
                                fwrite(ctx->AcC,(size_t)(jlen),1,ctx->PMF1d);
                                fprintf(ctx->PMF1d," ");
                            }
                            else
                                put_str(ctx, ctx->AcC,jlen,ii,nrec,0);
                            ii++;
                            jlen = jo = 0;
                        }
                    }
                }
                else if (sval[i] == 254) {
                    if (jo > 0 || ctx->AcN[j] > 0) {
                        if (ctx->AcN[j] > 0) {
                            jlen = ctx->AcN[j];
                            jo = 0;
                        }
                        snprintf(ctx->AcC + jo,9,"        ");
                        jo += 8;
                        if (jo >= jlen) {
                            if (ctx->PMF1Def) {
                                fwrite(ctx->AcC,(size_t)(jlen),1,ctx->PMF1d);
                                fprintf(ctx->PMF1d," ");
                            }
                            else
                                put_str(ctx, ctx->AcC,jlen,ii,nrec,0);
                            ii++;
                            jlen = jo = 0;
                        }
                    }
                    else {
                        x = ctx->MBlnkVal;
                        if (ctx->PMF1Def)
                            rt_fprintf_d(ctx, ctx->PMF1d,ctx->VPFmtS[ii],x);
                        else
                            put_data(ctx, x,ii,nrec);
                        ii++;
                    }
                }
                else if (sval[i] == 255) {
                    x = ctx->PMMSYS;
                    if (ctx->PMF1Def)
                        rt_fprintf_d(ctx, ctx->PMF1d,ctx->VPFmtS[ii],x);
                    else
                        put_data(ctx, x,ii,nrec);
                    ii++;
                }
                else { 
                    x = (double)sval[i];
                    if (ctx->PMF1Def)
                        rt_fprintf_d(ctx, ctx->PMF1d,ctx->VPFmtS[ii],x);
                    else
                        put_data(ctx, x,ii,nrec);
                    ii++;
                }

                if (++j >= nobs) {
                    ii = j = 0;
                    if (ctx->PMF1Def)
                        fprintf(ctx->PMF1d,"\n");
                    nrec++;
                    prn_message(ctx, nrec,0,0);
                }
            }
      
        }
    }
    else {                  /* uncompressed data */
       
        while (1) {
            ii = 0;
            for (j = 0; j < nobs; ++j) {

                if (ctx->AcN[j] == 0) {
                    x = sav_getd(ctx, arch,&r);
                    if (r) {
                        if (j == 0)
                            goto RSP1Nxt1;
                        goto RSP1Fin;
                    }
                    if (x <= smiss) {
                        nsys2++;
                        x = ctx->PMMSYS;
                    }
                    if (ctx->PMF1Def)
                        rt_fprintf_d(ctx, ctx->PMF1d,ctx->VPFmtS[ii],x);
                    else
                        put_data(ctx, x,ii,nrec);
                    ii++;
                }
                else {
                    if (ctx->AcN[j] > 0) {
                        jlen = ctx->AcN[j];
                        jo = (jlen - 1) / 8 + 1;
                        jo *= 8;
                        if ((r = rd_dat(ctx, jo,ctx->AcC)) != jo)
                            goto RSP1Fin; 
                        if (ctx->PMF1Def) {
                            fwrite(ctx->AcC,(size_t)(jlen),1,ctx->PMF1d);
                            fprintf(ctx->PMF1d," ");
                        }      
                        else
                            put_str(ctx, ctx->AcC,jlen,ii,nrec,0);
                        ii++;
                    }
                }
            }
            j = 0;
            nrec++;
            if (ctx->PMF1Def)
                fprintf(ctx->PMF1d,"\n");
            prn_message(ctx, nrec,0,0);
   
            if (nrec >= mxnoc)
                break;          
        }
    }

RSP1Nxt1:
    prn_message(ctx, nrec,1,0);
    printf1(ctx, "Read %d records.\n",nrec);

    if (ctx->PMF1Def) {
        printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);
        if (nrec <= 0)  
            goto RSP1Fin;
    }
    else {
        if (nrec <= 0) {
            printf1(ctx, "Error: no data matrix created.\n");
            goto RSP1Fin;
        }
        ctx->NOCDM = ctx->NOC = nrec;
        ctx->DMDef = 1;
        printf1(ctx, "Created a data matrix with %d variables and %d cases.\n",ctx->NVAR,ctx->NOC);
    }   
    err = 0;

RSP1Fin:
    if (err || ctx->PMF1Def) {
        if (idxn < get_nidx(ctx))
            clear_avar(ctx, idxn);
        ctx->NVAR = ctx->NOCMaxA = ctx->NOCDM = ctx->NOC = 0;
    }
    p_clean(ctx);
    return(err);
}


/* ------------------------------------------------------------------------ */
/*  rd_dat(n,buf)   read n bytes into buffer.                               */
/*                                                                          */

int rd_dat(TDAContext *ctx, int n,char *buf)
{
    int r;
    if (((r = (int)fread(buf,sizeof(char),(size_t)(n),ctx->PMFd))) != n) {     
        if (r == 0) {
            *buf = '\0';
            return(0);
        }
        printf1(ctx, "Error in reading the input file.\n");
        return(-1);
    }
    *(buf + r) = '\0';
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  sav_sval(sval)      get 8 bytes.                                        */

int sav_sval(TDAContext *ctx, int *sval)
{
    register int i,n;
    unsigned char buf[9];

    if ((n = rd_dat(ctx, 8,(char *)buf)) != 8) {
        if (n)
            return(-1);
        return(0);
    }
    n = 0;
    for (i = 0; i < 8; ++i) {
        if (buf[i] == 0)
            continue;
        else if (buf[i] == 0xff)
            sval[n++] = 255;   
        else if (buf[i] == 0xfe)
            sval[n++] = 254;    
        else if (buf[i] == 0xfd)
            sval[n++] = 253;
        else if (buf[i] == 0xfc)
            sval[n++] = 252;
        else  
            sval[n++] = (int)buf[i] - 100;
    }
    return(n);
}

/* -##--------------------------------------------------------------------- */
/*  sav_geti(arch,err)       return integer.                                */

int sav_geti(TDAContext *ctx, int arch,int *err)
{
    int r,n;
    char rbuf[6],buf[4];
    register char *p;

    *err = 0;
    if (((r = (int)fread(rbuf,sizeof(char),4,ctx->PMFd))) != 4) {           
        *err = -1;
        printf1(ctx, "Error: can't read next integer.\n");
        return(0);
    }
    p = rbuf;

    if (arch == 1) {
        if (ARCHTyp == 1) {
            buf[3] = *p++;
            buf[2] = *p++;
            buf[1] = *p++;
            buf[0] = *p;
        }
        else {
            buf[0] = *p++;
            buf[1] = *p++;
            buf[2] = *p++;
            buf[3] = *p;
        }
    }
    else {
        if (ARCHTyp == 2) {
            buf[3] = *p++;
            buf[2] = *p++;
            buf[1] = *p++;
            buf[0] = *p;
        }
        else {
            buf[0] = *p++;
            buf[1] = *p++;
            buf[2] = *p++;
            buf[3] = *p;
        }
    }
    /* buf is a char array and carries no alignment, so the bytes are
       copied into the int rather than read through a cast pointer */
    memcpy(&n,buf,sizeof(n));
    return(n);
}

/* ------------------------------------------------------------------------ */
/*  sav_getd(arch,err)      return double.                                  */
     
double sav_getd(TDAContext *ctx, int arch,int *err)
{
    double x;
    int r;
    char rbuf[10],buf[8];
    register char *p;

    *err = 0;
    if (((r = (int)fread(rbuf,sizeof(char),8,ctx->PMFd))) != 8) {           
        *err = -1;
        printf1(ctx, "Error (or eof): can't read next double.\n");
        return(0);
    }
    p = rbuf;

    if (arch == 1) {
        if (ARCHTyp == 1) {
            buf[7] = *p++;
            buf[6] = *p++;
            buf[5] = *p++;
            buf[4] = *p++;
            buf[3] = *p++;
            buf[2] = *p++;
            buf[1] = *p++;
            buf[0] = *p;
        }
        else {
            buf[0] = *p++;
            buf[1] = *p++;
            buf[2] = *p++;
            buf[3] = *p++;
            buf[4] = *p++;
            buf[5] = *p++;
            buf[6] = *p++;
            buf[7] = *p;
        }
    }
    else {
        if (ARCHTyp == 2) {
            buf[7] = *p++;
            buf[6] = *p++;
            buf[5] = *p++;
            buf[4] = *p++;
            buf[3] = *p++;
            buf[2] = *p++;
            buf[1] = *p++;
            buf[0] = *p;
        }
        else {
            buf[0] = *p++;
            buf[1] = *p++;
            buf[2] = *p++;
            buf[3] = *p++;
            buf[4] = *p++;
            buf[5] = *p++;
            buf[6] = *p++;
            buf[7] = *p;
        }
    }
    memcpy(&x,buf,sizeof(x));   /* see sav_geti(): buf is not aligned */
    return(x);
}

/*--------------------------------------------------------------------------*/
/*  wr_spss1        Write SPSS sav file.                                    */
/*  ####                                                                    */
/*      wspss1(                                                             */
/*              keep = varlist,                                             */
/*              drop = varlist,                                             */
/*              sort = varlist,                                             */
/*      ) = output-file-name;                                               */
/*                                                                          */
/*      Return 0 if OK, -1 if error.                                        */

int wr_spss1(TDAContext *ctx)
{
    register int i,j,k,l;
    int err,nv,sflag,nobs,n,ii,jj,olen = 0,lbuf[9];
    register char *p;
    unsigned char buf[9];
    double x,dbuf[9];
    char *pbuf[9];




    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 6,13,1)) {   /* get parameters */
        p_clean(ctx);
        return(-1);   
    }       
    printf1(ctx, "Writing SPSS sav file: %s\n",ctx->PMFdName);

    if (ctx->PMKeep && ctx->PMDrop) {
        p_err(ctx, -16,1);
        goto WSP1Fin;
    }

    /* get list of variables in AcI[] */

    if (alloc_aci(ctx, imax(ctx, ctx->NVAR,ctx->PMNV)))
        goto WSP1Fin;

    nv = 0;
    if (ctx->PMKeep) {
        for (k = 0; k < ctx->PMNV; ++k)  
            ctx->AcI[nv++] = ctx->PMVIdx[k];
    }
    else {
        j = ctx->VIFirst;
        while (j >= 0) {
            i = 1;
            if (ctx->PMDrop) {
                for (k = 0; k < ctx->PMNV; ++k) {
                    if (ctx->PMVIdx[k] == j) {
                        i = 0;
                        break;
                    }   
                }
            }
            if (i) 
                ctx->AcI[nv++] = j;
            j = ctx->VNxt[j];
        }
    }
    if (nv == 0) {
        printf1(ctx, "No variables selected.\n"); 
        goto WSP1Fin;
    }
    sflag = 0;
    if (ctx->PM1NV > 0) {        /* sort */
        err = vsort(ctx, ctx->PM1NV,ctx->PM1VIdx,1,0,1);
        if (err)
            goto WSP1Fin;
        sflag = 1;
    }

    nobs = 0;
    for (i = 0; i < nv; ++i) {
        j = ctx->AcI[i];
        if (ctx->VTyp[j] != 1) {
            nobs++;
        }
        else {
            l = -ctx->VSLen[j];
            n = (l - 1) / 8 + 1;
            nobs += n;
        }
    }
    if (alloc_ack(ctx, nobs + 1))    /* tda var index */
        goto WSP1Fin;
    if (alloc_acj(ctx, nobs + 1))    /* string length */
        goto WSP1Fin;
    if (alloc_acs(ctx, nobs + 1))    /* var label */
        goto WSP1Fin;
    if (alloc_acn(ctx, nobs + 1))    /* type */
        goto WSP1Fin;
    if (alloc_acm(ctx, nobs + 1))    /* format */
        goto WSP1Fin;

    nobs = 0;
    for (i = 0; i < nv; ++i) {
        j = ctx->AcI[i];
        if (ctx->VTyp[j] != 1) {
            ctx->AcK[nobs] = j;
            ctx->AcN[nobs] = 0;
            if (ctx->VPFmt1[i] > 0 && ctx->VPFmt2[i] < ctx->VPFmt1[i] - 1)
                ctx->AcM[nobs] = (int)ctx->VPFmt2[i] + (int)ctx->VPFmt1[i] * 256 + 5 * 256 * 256;
            else
                ctx->AcM[nobs] = 4 + 10 * 256 + 5 * 256 * 256;
            nobs++;
        }
        else {
            l = -ctx->VSLen[j];
            n = (l - 1) / 8 + 1;
            ctx->AcK[nobs] = j;
            ctx->AcN[nobs] = l;
            ctx->AcJ[nobs] = l;
            ctx->AcM[nobs] = l * 256 + 256 * 256;
            nobs++;
            for (k = 1; k < n; ++k) {
                ctx->AcK[nobs] = j;
                ctx->AcN[nobs] = -1;
                ctx->AcJ[nobs] = l;
                ctx->AcM[nobs] = l * 256 + 256 * 256;
                nobs++;
            }
        }
    }

    /* print header */

    fwrite("$FL2@(#)SPSS DATA FILE. SAV FILE CREATED BY TDA.",48,1,ctx->PMFd);
    fprnchar(ctx, ctx->PMFd,' ',16,0);

    sav_puti(ctx, 2);            /* file layout code */
    sav_puti(ctx, nobs);         /* number of OBS */
    sav_puti(ctx, 1);            /* compression code, always compressed */
    sav_puti(ctx, 0);            /* case weight variable index */
    sav_puti(ctx, ctx->NOC);          /* number of cases */
    sav_putd(ctx, 100.0);        /* compression bias */

    fprnchar(ctx, ctx->PMFd,' ',84,0);        /* no date, etc. */ 

    for (i = 0; i < nobs; ++i) {    /* write dictionary */

        j = ctx->AcK[i];
        sav_puti(ctx, 2);            /* record type */
        sav_puti(ctx, ctx->AcN[i]);       /* var type */

        ii = 0;
        if (ctx->VLabel[j] != NULL) {
            ii = (int)(strlen(ctx->VLabel[j]));
            n = (ii - 1) / 4;
            n++;
            n *= 4;
            sav_puti(ctx, 1);            /* label */
        }
        else
            sav_puti(ctx, 0);            /* no label */
        sav_puti(ctx, 0);            /* no missing value codes */
        sav_puti(ctx, ctx->AcM[i]);       /* print format */
        sav_puti(ctx, ctx->AcM[i]);       /* write format */

        p = ctx->VName[j];
        for (k = 0; k < 8; ++k) {
            if (*p)
                fprintf(ctx->PMFd,"%c",*p++);
            else
                fprintf(ctx->PMFd," ");
        }
        if (ii > 0) {
            sav_puti(ctx, ii);   
            fprintf(ctx->PMFd,"%s",ctx->VLabel[j]);
            for (k = ii; k < n; ++k)
                fprintf(ctx->PMFd," ");
        }
    }

    /* write data, always compressed */

    sav_puti(ctx, 999);          /* record type */
    sav_puti(ctx, 0);            /* filler */
         
    k = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        ii = i;
        if (sflag)
            ii = ctx->VSORTPtr[i];

        for (j = 0; j < nobs; ++j) {

            jj = ctx->AcK[j];
            if (ctx->AcN[j] == 0) {
                x = get_data(ctx, jj,ii);
                n = (int)x;
                if (x == (double)n && n >= -99 && n <= 151) {
                    n += 100;
                    buf[k] = (unsigned char)n;
                }   
                else {
                    buf[k] = 0xfd;
                    dbuf[k] = x;
                    pbuf[k] = NULL;
                }
            }
            else {
                buf[k] = 0xfd;
                if (ctx->AcN[j] > 0) {
                    pbuf[k] = ctx->VDPtr[jj] + ii * ctx->AcJ[j];                
                    olen = 8;
                }   
                else {
                    pbuf[k] = ctx->VDPtr[jj] + ii * ctx->AcJ[j] + olen;           
                    olen += 8;
                }   
                
                if (olen <= ctx->AcJ[j])  
                    lbuf[k] = 8;
                else
                    lbuf[k] = 8 - (olen - ctx->AcJ[j]);
            }
            if (++k >= 8) {
                fwrite(buf,8,1,ctx->PMFd);
                for (k = 0; k < 8; ++k) {
                    if (buf[k] == 0xfd) {
                        if (pbuf[k] == NULL) {
                            sav_putd(ctx, dbuf[k]);
                        }
                        else {
                            fwrite(pbuf[k],(size_t)(lbuf[k]),1,ctx->PMFd);
                            for (l = lbuf[k]; l < 8; ++l)
                                fprintf(ctx->PMFd," ");
                        }
                    }
                }
                k = 0;
            }
        }
        prn_message(ctx, i + 1,0,1);
    }
    if (k > 0) {
        while (k < 8)
            buf[k++] = 0xfc;          /* end of file */

        fwrite(buf,8,1,ctx->PMFd);
        for (k = 0; k < 8; ++k) {
            if (buf[k] == 0xfd) {
                if (pbuf[k] == NULL)
                    sav_putd(ctx, dbuf[k]);
                else {
                    fwrite(pbuf[k],(size_t)(lbuf[k]),1,ctx->PMFd);
                    for (l = lbuf[k]; l < 8; ++l)
                        fprintf(ctx->PMFd," ");
                }
            }
        }
    }
    prn_message(ctx, ctx->NOC,1,1);

    printf1(ctx, "%d records with %d variables written to: %s\n",ctx->NOC,nv,ctx->PMFdName);
    err = 0;
          
WSP1Fin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sav_puti(n)     write integer.                                          */

void sav_puti(TDAContext *ctx, int n)
{
    char buf[5];
    register char *p;

    p = (char *)&n;
    if (ARCHTyp == 1) {
        buf[3] = *p++;
        buf[2] = *p++;
        buf[1] = *p++;
        buf[0] = *p;
    }
    else {
        buf[0] = *p++;
        buf[1] = *p++;
        buf[2] = *p++;
        buf[3] = *p;
    }
    fwrite(buf,4,1,ctx->PMFd);
}

/* ------------------------------------------------------------------------ */
/*  sav_putd(x)     write double.                                           */
     
void sav_putd(TDAContext *ctx, double x)
{
    char buf[9];
    register char *p;

    p = (char *)&x;
    if (ARCHTyp == 1) {
        buf[7] = *p++;
        buf[6] = *p++;
        buf[5] = *p++;
        buf[4] = *p++;
        buf[3] = *p++;
        buf[2] = *p++;
        buf[1] = *p++;
        buf[0] = *p;
    }
    else {
        buf[0] = *p++;
        buf[1] = *p++;
        buf[2] = *p++;
        buf[3] = *p++;
        buf[4] = *p++;
        buf[5] = *p++;
        buf[6] = *p++;
        buf[7] = *p;
    }
    fwrite(buf,8,1,ctx->PMFd);
}

/* ------------------------------------------------------------------------ */
/*  get_dvarp(np)   get info about dvar partition.                          */
     
int get_dvarp(TDAContext *ctx, int *np,char *vname)
{
    register int j;   

    if (ctx->PMDVARP > 0) {
        *np = 1;
        return(ctx->PMDVARP);
    }
    else if (ctx->PMDVARPN > 0) {
        *np = 2;
        return(ctx->PMDVARPN);
    }
    else if (ctx->PMDVARVN > 0) {
        for (j = 0; j < ctx->PMDVARVN; ++j) {
            if (!strcmp(ctx->PMDVARVName[j],vname)) {
                *np = ctx->PMDVARVNP[j];
                return(ctx->PMDVARVNL[j]);
            }
        }
    }
    *np = 0;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  make_arcd()    update archive description file.                         */
     
void make_arcd(TDAContext *ctx, int fn,char *fname,int len,int nrec,int vn)
{
    int addc,n;

    addc = 1;
#if S_DOS               /* add to record length */
    addc++;
#endif
    n = 0;
    if (ctx->PMDVARFDef) {
        fclose(ctx->PMDVARFd);
        if ((ctx->PMDVARFd = fopen(ctx->PMDVARFName,OPEN_RD))) {   
            if (alloc_acc(ctx, 1001) == 0) {  
                while (fgets(ctx->AcC,1000,ctx->PMDVARFd))  
                    n++;
                alloc_acc(ctx, 0);  
            }
        }
    }
    if (ctx->PMARCFZOO)
        fprintf(ctx->PMARCFd,"%s\n",ctx->PMARCFZOOF);
    fprintf(ctx->PMARCFd,"%3d %s ",fn,fname);
    fprnchar(ctx, ctx->PMARCFd,' ',(int)(20 - strlen(fname)),0);
    fprintf(ctx->PMARCFd,"1 %6d %8d %6d\n",len + addc,nrec,vn);
    if (ctx->PMARCFVDF) {
        fprintf(ctx->PMARCFd,"999 %s ",ctx->PMARCFVDFF);
        fprnchar(ctx, ctx->PMARCFd,' ',(int)(20 - strlen(ctx->PMARCFVDFF)),0);
        fprintf(ctx->PMARCFd,"2 %6d %8d %6d\n",0,n,0);
    }
    printf1(ctx, "Added info to archive description file: %s\n",ctx->PMARCFName);
}


