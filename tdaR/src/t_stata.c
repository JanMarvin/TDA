/****************************************************************************/
/*  t_stata                                                                 */
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

/*  functions in t_stata.c */

int rd_stata(TDAContext *ctx); 
int st_gets(TDAContext *ctx, char *p);
int st_geti(TDAContext *ctx, char *p);
double st_getsf(TDAContext *ctx, char *p,int l); 
double st_getf(TDAContext *ctx, char *p);
double st_getd(TDAContext *ctx, char *p);
int wr_stata(TDAContext *ctx);

/*--------------------------------------------------------------------------*/
/*  Global variables                                                        */


/* ------------------------------------------------------------------------ */
/*  rd_stata        Read a STATA data file and create an internal data      */
/*                  matrix. Note that this function is only called if an    */  
/*                  an internal data matrix does not already exist.         */
/*                                                                          */
/*                  rstata(                                                 */
/*                      noc=...,     maximum number of cases                */
/*                      vsel=...,    case selection                         */
/*                      msys=...,    system missing value code, def. -5     */
/*                      dvar()=...,  create variable description file       */
/*                      df=...,      write data to output file              */
/*                      n = ...,     1 default var name convention          */         
/*                                   2 translate var names to upper case    */
/*                  ) = fname;                                              */
/*                                                                          */
/*  Version 3 :     without value labels                                    */
/*          4 :     without value labels                                    */
/*          6 :                                                             */
/*          7 :                                                             */
/*         10 :     format-114 datasets                                     */
/*                                                                          */
/*                  Return 0 if OK, -1 if error.                            */

int rd_stata(TDAContext *ctx)  
{
    register int i,j,k;
    register char *p,*q,*vl;
    char buf[20000],*typlist,*varlist,*srtlist,*fmtlist,*lbllist;
    char *varlab,*rbuf;
    unsigned char u;
    int err,idxn,nvar,nrec,nrec1,nrec2,n,rlen,nmiss,nmiss1,w,d,hlen,fmtlen,typ;
    int len,typlista,varlista,srtlista,fmtlista,lbllista,varlaba,rbufa;
    int vn = 0,m,nn,np,off = 0,off1,mxnoc,aflag,vlen,vlablen,fptr;
    double x;

    if (check_cmd(ctx, 0))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 6,10,1)) {   /* get parameters */
        p_clean(ctx);
        return(-1);   
    }   
    printf1(ctx, "Reading Stata file: %s\n",ctx->PMFdName);
    if (ctx->DMDef) {
        printf1(ctx, "Error: a data matrix already exists.\n");
        p_clean(ctx);
        return(-1);
    }
    aflag = 0;
    if (ctx->PMARCFDef && ctx->PMDVARFDef && ctx->PMF1Def)
        aflag = 1;

    if (ctx->PMNOCFlg)                     
        printf1(ctx, "Maximum number of cases: %d\n",ctx->PMNOC);

    if (ctx->PMRHSTRA && ctx->PMF1Def == 0)
        printf1(ctx, "Record selection (vsel): %s\n",ctx->PMRHSTR + 5);

    idxn = get_nidx(ctx);      /* index to first new variable */

    rbufa = varlaba = lbllista = fmtlista = srtlista = varlista = typlista = 0;
    err = -1;
           
    if (fread(buf,sizeof(char),10,ctx->PMFd) != 10) {
        p_err(ctx, -7,1);
        goto RSTATFin;
    }
    p = buf;
    printf1(ctx, "Release: 0x%02x ",*p);
    if (*p == 0x69) {
        typ = 4;
        hlen = 50;
        fmtlen = 12;
        vlen = 9;
        vlablen = 32;
    }
    else if (*p == 0x68) {
        typ = 3;
        hlen = 32;
        fmtlen = 7;
        vlen = 9;
        vlablen = 32;
    }
    else if (*p == 0x66) {
        /* format 102's data label is 30 bytes, not 32 (authoritative:
           readstata13's writer, save_pre13_dta.cpp) -- sharing 0x68's
           branch shifted every later section by two bytes and made
           genuine 102 files unreadable. */
        typ = 3;
        hlen = 30;
        fmtlen = 7;
        vlen = 9;
        vlablen = 32;
    }
    else if (*p == 0x6c) {
        typ = 6;
        hlen = 99;
        fmtlen = 12;
        vlen = 9;
        vlablen = 81;
    }
    else if (*p == 0x6e) {
        typ = 7;
        hlen = 99;
        fmtlen = 12;
        vlen = 33;
        vlablen = 81;
    }
    else if (*p == 0x72) {
        typ = 10;
        hlen = 99;
        fmtlen = 49;
        vlen = 33;
        vlablen = 81;
    }
    else {
        printf1(ctx, "Can't identify version.\n");
        goto RSTATFin;
    }
    printf1(ctx, "[version %d]\n",typ);

    if (*++p == 0x01)
        ctx->HILO = 1;
    else if (*p == 0x02)
        ctx->HILO = 2;
    else {
        printf1(ctx, "Can't identify HiLo indicator. Found: 0x%02x\n",(unsigned char)*p);
        goto RSTATFin;
    }
    if (*++p != 0x01) {
        printf1(ctx, "Probably not a stata data file (found: 0x%02x)\n",(unsigned char)*p);
        goto RSTATFin;
    }
    nvar = st_gets(ctx, buf + 4);
    nrec  = st_geti(ctx, buf + 6);

    printf1(ctx, "Number of variables: %d\n",nvar);
    printf1(ctx, "Number of cases: %d\n",nrec);

    if (nvar < 1 || nrec < 1) {
        printf1(ctx, "Probably an error.\n");
        goto RSTATFin;
    }
    if (nvar > ctx->MaxNV) {
        printf1(ctx, "Error: exceeded max number of variables.\n");
        goto RSTATFin;
    }
    if (fread(buf,sizeof(char),(size_t)(hlen),ctx->PMFd) != (size_t)(hlen)) {
        p_err(ctx, -7,1);
        goto RSTATFin;
    }
    printf1(ctx, "Data label: %s\n",buf);
    if (typ == 10) {
        printf1(ctx, "Time stamp: %s\n",buf + 81);
    }
    newline(ctx);

    if (!(typlist = (char *) calloc((size_t)(nvar + 1),sizeof(char)))) {
        p_err(ctx, -2,1);
        goto RSTATFin;
    }
    typlista = nvar + 1;
    memrq(ctx, typlista,1);
         
    if (!(varlist = (char *) calloc((size_t)(vlen) * (size_t)(nvar) + 1,sizeof(char)))) {
        p_err(ctx, -2,1);
        goto RSTATFin;
    }
    varlista = vlen * nvar + 1;
    memrq(ctx, varlista,1);

    if (!(srtlist = (char *) calloc((size_t)(2) * (size_t)((nvar + 1)) + 1,sizeof(char)))) {
        p_err(ctx, -2,1);
        goto RSTATFin;
    }
    srtlista = 2 * (nvar + 1) + 1;
    memrq(ctx, srtlista,1);

    if (!(fmtlist = (char *) calloc((size_t)(fmtlen) * (size_t)(nvar) + 1,sizeof(char)))) {
        p_err(ctx, -2,1);
        goto RSTATFin;
    }
    fmtlista = fmtlen * nvar + 1;
    memrq(ctx, fmtlista,1);

    if (!(lbllist = (char *) calloc((size_t)(vlen) * (size_t)(nvar) + 1,sizeof(char)))) {
        p_err(ctx, -2,1);
        goto RSTATFin;
    }
    lbllista = vlen * nvar + 1;
    memrq(ctx, lbllista,1);

    if (!(varlab = (char *) calloc((size_t)(vlablen) * (size_t)(nvar) + 1,sizeof(char)))) {
        p_err(ctx, -2,1);   
        goto RSTATFin;
    }
    varlaba = vlablen * nvar + 1;
    memrq(ctx, varlaba,1);
      
    if (fread(typlist,sizeof(char),(size_t)(nvar),ctx->PMFd) != (size_t)(nvar)) {
        p_err(ctx, -7,1);
        goto RSTATFin;
    }
    if (fread(varlist,sizeof(char),(size_t)(vlen * nvar),ctx->PMFd) != (size_t)(vlen * nvar)) {
        p_err(ctx, -7,1);
        goto RSTATFin;
    }
    if (fread(srtlist,sizeof(char),(size_t)(2 * (nvar + 1)),ctx->PMFd) != (size_t)(2 * (nvar + 1))) {
        p_err(ctx, -7,1);
        goto RSTATFin;
    }
    if (fread(fmtlist,sizeof(char),(size_t)(fmtlen * nvar),ctx->PMFd) != (size_t)(fmtlen * nvar)) {
        p_err(ctx, -7,1);
        goto RSTATFin;
    }
    if (fread(lbllist,sizeof(char),(size_t)(vlen * nvar),ctx->PMFd) != (size_t)(vlen * nvar)) {
        p_err(ctx, -7,1);
        goto RSTATFin;
    }
    if (fread(varlab,sizeof(char),(size_t)(vlablen * nvar),ctx->PMFd) != (size_t)(vlablen * nvar)) {
        p_err(ctx, -7,1);
        goto RSTATFin;
    }
    rlen = 0;   /* record length */

    if (alloc_acj(ctx, nvar))        /* used for storage size */
        goto RSTATFin;
    if (alloc_acr(ctx, nvar))        /* used for format w */
        goto RSTATFin;
    if (alloc_acs(ctx, nvar))        /* used for format d */
        goto RSTATFin;
         
    for (i = 0; i < nvar; ++i) {

        /********************************  
         printf1(ctx, "%d %c :%s:%s:%s:%s:\n",
         i,typlist[i],varlist + i * vlen,fmtlist + fmtlen * i,lbllist + i*vlen, 
         varlab + i * vlablen);
        **********************/
        p = varlist + i * vlen;

        if (typ >= 7)                   /* set max length of var name */
            *(p + VNLMax) = '\0';

        if (isalpha((int)*p)) {
            *p = (char)toupper((int)*p);
            p = buf;
        }
        else {
            buf[0] = '_';
            p = buf + 1;
        }

        /* get storage size */

        len = 1;

        if (typ == 10) {
            if ((unsigned char)typlist[i] == 0xfb) {
                n = 1;
                rlen += 1;
            }
            else if ((unsigned char)typlist[i] == 0xfc) {
                n = 2;
                rlen += 2;
            }
            else if ((unsigned char)typlist[i] == 0xfd) {
                n = 5;
                rlen += 4;
            }
            else if ((unsigned char)typlist[i] == 0xfe) {
                n = 4;
                rlen += 4;
            }
            else if ((unsigned char)typlist[i] == 0xff) {
                n = 8;
                rlen += 8;
            }
            else {
                u = (unsigned char)typlist[i];
                len = (int)u;
                rlen += len;        
                n = -len;
            }
        } 
        else {
            if (typlist[i] == 'b') {
                n = 1;
                rlen += 1;
            }
            else if (typlist[i] == 'i') {
                n = 2;
                rlen += 2;
            }
            else if (typlist[i] == 'l') {
                n = 5;
                rlen += 4;
            }
            else if (typlist[i] == 'f') {
                n = 4;
                rlen += 4;
            }
            else if (typlist[i] == 'd') {
                n = 8;
                rlen += 8;
            }
            else {
                u = ((unsigned char)typlist[i]) - 0x7f;
                len = (int)u;
                rlen += len;        
                n = -len;
            }
        }

        /* get print format */

        w = d = 0;
        if (n > 0) { /* numerical */
            q = fmtlist + fmtlen * i;
            if (sscanf(++q,"%d.%d",&w,&d) == 2 || sscanf(++q,"%d,%d",&w,&d) == 2) {
                q = skip_int(ctx, q);
                q = skip_int(ctx, ++q);
                if (*q == 'e')
                    w = -w;
                else if (*q == 'g')     /* free format */
                    w = d = 0;
            }
            len = iabs(ctx, w);
        }
        snprintf(p,sizeof(buf) - (size_t)(p - buf),"%s<%d>[%d.%d]",varlist + i * vlen,n,w,d);
    
        if (ctx->PMN == 2) {
            q = buf;
            while (*q != '<') {
                *q = (char)toupper((int)*q);
                q++;
            }
        }

        if (*(varlab + i * vlablen)) {
            p = varlab + i * vlablen;
            while (*p) {
                if (*p == '(')
                    *p = '[';
                else if (*p == ')')
                    *p = ']';
                p++;
            }
            p = buf + strlen(buf);
            snprintf(p,sizeof(buf) - (size_t)(p - buf),"(%s)",varlab + i * vlablen);
        }
        p = buf + strlen(buf);

        if (n < 0)
            snprintf(p,sizeof(buf) - (size_t)(p - buf),"=stata(%d)",len);     /* string variable */
        else if (typ == 10)
            snprintf(p,sizeof(buf) - (size_t)(p - buf),"=stata(0x%02x)",(unsigned char)typlist[i]);
        else                 
            snprintf(p,sizeof(buf) - (size_t)(p - buf),"=stata(%c)",typlist[i]);

        ctx->AcJ[i] = n;
        ctx->AcR[i] = w;
        ctx->AcS[i] = d;
                            
        if (save_var(ctx, buf,0)) {       /* save variable definition */

            printf1(ctx, "\nError: can't save variable definitions.\n");
            goto RSTATFin;              
        }
    }
    prn_var(ctx, idxn);      /* print list of new variables */
    printf1(ctx, "\n");

    /* skip expansion fields */
      
    while (typ >= 4) {      /* only for version 4 and higher */
         
        if (typ <= 6) {
            if (fread(buf,sizeof(char),3,ctx->PMFd) != 3) {
                p_err(ctx, -7,1);
                goto RSTATFin;
            }
            n = st_gets(ctx, buf + 1);
        }
        else {
            if (fread(buf,sizeof(char),5,ctx->PMFd) != 5) {
                p_err(ctx, -7,1);
                goto RSTATFin;
            }
            n = st_geti(ctx, buf + 1);
        }
        if (n == 0 && !*buf)
            break;

        if (fread(buf,sizeof(char),(size_t)(n),ctx->PMFd) != (size_t)(n)) {
            p_err(ctx, -7,1);
            goto RSTATFin;
        }
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
                goto RSTATFin;
            }
        }
    }

    /* read data to create internal data matrix, or directly writing
       to an output file. */

    if (rlen < 1) {
        printf1(ctx, "Error in record length.\n");
        goto RSTATFin;
    }
    if (!(rbuf = (char *) calloc((size_t)(rlen + 1),sizeof(char)))) {
        p_err(ctx, -2,1);
        goto RSTATFin;
    }
    rbufa = rlen + 1;
    memrq(ctx, rbufa,1);

    if (ctx->PMNOCFlg)
        mxnoc = imin(ctx, ctx->PMNOC,nrec);
    else 
        mxnoc = nrec;   

    if (ctx->PMF1Def) {
        printf1(ctx, "Data will be directly written to output file: %s\n",ctx->PMF1dName);
    }
    else {
        ctx->NOCMaxA = mxnoc;
        printf1(ctx, "Reading data to create internal data matrix.\n");
        printf1(ctx, "Maximum number of cases: %d\n",mxnoc);
        printf1(ctx, "Record length: %d\n",rlen);

        /* allocate memory for data */
        
        n = ctx->MemReq;
        if (alloc_vdat(ctx, idxn,1)) {
            printf1(ctx, "Error: insufficient memory for data matrix (%d cases).\n",ctx->NOCMaxA);
            goto RSTATFin;
        }
        printf1(ctx, "Allocated %d bytes for data matrix.\n\n",ctx->MemReq - n);
    }
    nrec2 = nrec1 = 0;      /* number of records */
    nmiss = 0;              /* number of system missing values */
           
    if (ctx->SILENTFlg < 2)
        printfe(ctx, "Reading: %s\n",ctx->PMFdName);

    for (nrec2 = 0; nrec2 < mxnoc; ++nrec2) {

        n = (int)(fread(rbuf,sizeof(char),(size_t)(rlen),ctx->PMFd));             
        if (n == 0)
            break;

        if (n != rlen) {
            printf1(ctx, "Warning: break with reading %d instead of %d bytes.\n",n,rlen);   
            break;          
        }

        q = p = rbuf;
        nmiss1 = 0;

        for (j = 0; j < nvar; ++j) {
            w = 0;
          
            if (typ == 10) {
                if  ((unsigned char)typlist[j] == 0xfb) {
                    if (*p >= 0x65) {
                        x = ctx->PMMSYS;
                        nmiss1++;
                    }
                    else
                        x = (double)*p;
                    p += 1;
                }
                else if  ((unsigned char)typlist[j] == 0xfc) {
                    n = st_gets(ctx, p);
                    if (n >= 0x7fe5) {
                        x = ctx->PMMSYS;
                        nmiss1++;
                    }
                    else
                        x = (double)n;
                    p += 2;
                }
                else if  ((unsigned char)typlist[j] == 0xfd) {
                    n = st_geti(ctx, p);
                    if (n == 0x7fffffe5) {
                        x = ctx->PMMSYS;
                        nmiss1++;
                    }
                    else
                        x = (double)n;
                    p += 4;
                }
                else if  ((unsigned char)typlist[j] == 0xfe) {
                    x = st_getf(ctx, p);
                    if (x == ctx->PMMSYS)
                        nmiss1++;
                    p += 4;
                }   
                else if  ((unsigned char)typlist[j] == 0xff) {
                    x = st_getd(ctx, p);
                    if (x == ctx->PMMSYS)
                        nmiss1++;
                    p += 8;
                }
                else {
                    u = (unsigned char)typlist[j];
                    w = n = (int)u;
                    q = p;
                    p += n;          
                }
            }
            else {
                switch (typlist[j]) {
                    case 'b':   if (*p == 0x7f) {
                                    x = ctx->PMMSYS;
                                    nmiss1++;
                                }
                                else
                                    x = (double)*p;
                                p += 1;
                                break;
                    case 'i':   n = st_gets(ctx, p);
                                if (n == 0x7fff) {
                                    x = ctx->PMMSYS;
                                    nmiss1++;
                                }
                                else
                                    x = (double)n;
                                p += 2;
                                break;
                    case 'l':   n = st_geti(ctx, p);
                                if (n == 0x7fffffff) {
                                    x = ctx->PMMSYS;
                                    nmiss1++;
                                }
                                else
                                    x = (double)n;
                                p += 4;
                                break;
                    case 'f':   x = st_getf(ctx, p);
                                if (x == ctx->PMMSYS)
                                    nmiss1++;
                                p += 4;
                                break;
                    case 'd':   x = st_getd(ctx, p);
                                if (x == ctx->PMMSYS)
                                    nmiss1++;
                                p += 8;
                                break;
                    default:    u = ((unsigned char)typlist[j]) - 0x7f;
                                w = n = (int)u;
                                q = p;
                                p += n;          
                                break;
                }
            }
            if (w > 0) {    /* string */     

                vl = ctx->SVBuf;
                for (i = 0; i < w; ++i) {
                    if (*q) 
                        *vl++ = *q++;
                    else
                        *vl++ = ' ';
                }
                if (ctx->PMF1Def) {
                    fwrite(ctx->SVBuf,(size_t)(w),1,ctx->PMF1d);
                    fprintf(ctx->PMF1d," ");
                }
                else  
                    put_str(ctx, ctx->SVBuf,w,j,nrec1,0);        
            }
            else {
                if (ctx->PMF1Def)
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->VPFmtS[j],x);
                else
                    put_data(ctx, x,j,nrec1);
            }
        }
        if (ctx->PMF1Def)
            fprintf(ctx->PMF1d,"\n");

        prn_message(ctx, nrec2,0,0);

        if (ctx->PMRHSTRA && ctx->PMF1Def == 0) {         /* vsel */

            n = v_eval1(ctx, nrec1,ctx->ESCnt,ctx->ESTyp,ctx->ESVal,ctx->ESIdx,&x,0,0,0,0,0);
            if (n) {
                printf1(ctx, "Can't evaluate vsel expression in record %d.\n",nrec2 + 1);
                prn_emsg2(ctx, n);
                goto RSTATFin;
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
        printf1(ctx, "Selected: %d.",nrec1);
    printf1(ctx, "\n");

    if (ctx->PMF1Def) {
        printf1(ctx, "%d records written to: %s\n",nrec1,ctx->PMF1dName);
        if (nrec1 <= 0)  
            goto RSTATFin;
    }
    else {
        if (nrec1 <= 0) {
            printf1(ctx, "Error: no data matrix created.\n");
            goto RSTATFin;
        }
        ctx->NOCDM = ctx->NOC = nrec1;
        ctx->DMDef = 1;
    
        printf1(ctx, "Created a data matrix with %d variables and %d cases.\n",ctx->NVAR,ctx->NOC);
    }
    if (nrec1 < nrec)  
        printf1(ctx, "Warning: file contains more than %d records.\n",nrec1);
    printf1(ctx, "Number of system missing values (substituted by %lg): %d\n",
                                                  ctx->PMMSYS,nmiss);
     
    /* try to read value labels */

    if (ctx->PMDVARFDef) {

        fprintf(ctx->PMDVARFd,"# Variable description file based on: %s\n",ctx->PMFdName);

        if (typ >= 6 && ctx->PMNOCFlg == 0) {    /* only if no limit to number of cases */

            if (alloc_aci(ctx, nvar))            /* used for file pointer */
                goto RSTATFin;

            nn = 10000;
            if (alloc_acn(ctx, nn + 1))          /* used for off[] */
                goto RSTATFin;
            if (alloc_acm(ctx, nn + 1))          /* used for val[] */
                goto RSTATFin;

            np = 10000;
            if (alloc_acc(ctx, np + 2))          /* used for read buffer */
                goto RSTATFin;

            while (1) {
                fptr = (int)(ftell(ctx->PMFd));
                if (fread(buf,sizeof(char),(size_t)(vlen + 7),ctx->PMFd) != (size_t)(vlen + 7))  
                    break;
      
                len = st_geti(ctx, buf);
                    if (len < 1) {
                    p_err(ctx, -44,1);
                    break;
                }
                p = buf + 4;
                        
                for (i = 0; i < nvar; ++i) {
                    if (!strcmp(p,lbllist + i * vlen))  
                        ctx->AcI[i] = fptr;
                }
                if (np < len) {
                    np = len;
                    if (alloc_acc(ctx, np + 2))      /* used for read buffer */
                        goto RSTATFin;
                }
                if (fread(ctx->AcC,sizeof(char),(size_t)(len),ctx->PMFd) != (size_t)(len)) {
                    p_err(ctx, -44,1);
                    goto RSTATFin;
                }
                n = st_geti(ctx, ctx->AcC);
                m = st_geti(ctx, ctx->AcC + 4);
                if (n < 1 || m < 1) {
                    p_err(ctx, -44,1);
                    break;
                }
                if (nn < n) {
                    nn = n;
                    if (alloc_acn(ctx, nn + 1))      /* used for off[] */
                        goto RSTATFin;
                    if (alloc_acm(ctx, nn + 1))      /* used for val[] */
                        goto RSTATFin;
                }
            }
        }
        vn = off = 0;
        for (i = 0; i < nvar; ++i) {

            if (ctx->PMN == 2) {
                q = varlist + i * vlen;
                while (*q) {
                    fprintf(ctx->PMDVARFd,"%c",(char)toupper((int)*q));
                    q++;
                }
                fprintf(ctx->PMDVARFd," ");
            }
            else
                fprintf(ctx->PMDVARFd,"%s ",varlist + i * vlen);

            fprnchar(ctx, ctx->PMDVARFd,' ',(int)(12 - strlen(varlist + i * vlen)),0);
            fprintf(ctx->PMDVARFd,"%3d %5d ",ctx->PMDVARFN,off);

            n = ctx->AcJ[i];
            w = ctx->AcR[i];
            d = ctx->AcS[i];
            if (n < 0)
                len = -n;
            else
                len = iabs(ctx, w);

            if (n < 0)
                fprintf(ctx->PMDVARFd,"%4d     ",n);
            else
                fprintf(ctx->PMDVARFd,"%4d.%-2d  ",w,d);

            vl = varlab + i * vlablen;       
            for (j = 0; j < vlablen; ++j) {
                if (!*vl)  
                    break;
                fprintf(ctx->PMDVARFd,"%c",*vl++);
            }   
            fprintf(ctx->PMDVARFd,"\n");
            vn++;
            off1 = off;
            off += len + 1;

            if (typ >= 6 && ctx->PMNOCFlg == 0 && ctx->AcI[i] != 0) {      /* read value labels */

                fseek(ctx->PMFd,ctx->AcI[i],0);

                if (fread(buf,sizeof(char),(size_t)(vlen + 7),ctx->PMFd) == (size_t)(vlen + 7)) {

                    len = st_geti(ctx, buf);
                    if (len > 0) {
                        if (fread(ctx->AcC,sizeof(char),(size_t)(len),ctx->PMFd) == (size_t)(len)) {
                            nn = st_geti(ctx, ctx->AcC);
                            m = st_geti(ctx, ctx->AcC + 4);
                            if (nn > 0 && m > 0) {
                                p = ctx->AcC + 8;
                                for (k = 0; k < nn; ++k) {
                                    ctx->AcN[k] = st_geti(ctx, p);
                                    p += 4;
                                }
                                for (k = 0; k < nn; ++k) {
                                    ctx->AcM[k] = st_geti(ctx, p);
                                    p += 4;
                                }
                                for (k = 0; k < nn; ++k)  
                                    fprintf(ctx->PMDVARFd," %6d  %s\n",ctx->AcM[k],p + ctx->AcN[k]);
                            }
                        }
                    }
                }
            }

            if (n < 0) {     /* check for partition */
                                  
                nn = get_dvarp(ctx, &np,varlist + i * vlen); 

                if (nn > 0) {
                    if (nn >= len) {
                        if (np == 2)
                            nn = len;
                        else
                            nn = 0;
                    }
                }
                if (nn > 0) {
                    m = j = 0;
                    while (m < len) {
                        fprintf(ctx->PMDVARFd,"%s%02d ",varlist + i * vlen,j);
                        fprnchar(ctx, ctx->PMDVARFd,' ',(int)(10 - strlen(varlist + i * vlen)),0);
                        fprintf(ctx->PMDVARFd,"%3d %5d ",ctx->PMDVARFN,off1);

                        if (np == 1)
                            fprintf(ctx->PMDVARFd,"%4d     ",-nn);
                        else
                            fprintf(ctx->PMDVARFd,"%4d     ",nn);

                        vl = varlab + i * vlablen;       
                        if (*vl) {
                            for (k = 0; k < vlablen; ++k) {
                                if (!*vl)  
                                    break;
                                fprintf(ctx->PMDVARFd,"%c",*vl++);
                            }   
                            fprintf(ctx->PMDVARFd,", part %02d",j);
                        }
                        fprintf(ctx->PMDVARFd,"\n");
                        vn++;
                        off1 += nn;
                        m += nn;
                        j++;
                    }
                }
            }
        }
        printf1(ctx, "Variable description file written to: %s\n",ctx->PMDVARFName);
    }

    if (aflag)          /* archive description file */
        make_arcd(ctx, ctx->PMDVARFN,ctx->PMF1dName,off,nrec1,vn);

    err = 0;

RSTATFin:
    if (typlista) {
        free(typlist);
        memrq(ctx, -typlista,1);
    }
    if (varlista) {
        free(varlist);
        memrq(ctx, -varlista,1);
    }
    if (srtlista) {
        free(srtlist);
        memrq(ctx, -srtlista,1);
    }
    if (fmtlista) {
        free(fmtlist);
        memrq(ctx, -fmtlista,1);
    }
    if (lbllista) {
        free(lbllist);
        memrq(ctx, -lbllista,1);
    }
    if (varlaba) {
        free(varlab);
        memrq(ctx, -varlaba,1);
    }
    if (rbufa) {
        free(rbuf);
        memrq(ctx, -rbufa,1);
    }
    if (err || ctx->PMF1Def) {
        clear_avar(ctx, idxn);
        ctx->NVAR = ctx->NOCMaxA = ctx->NOCDM = ctx->NOC = 0;
    }
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  st_getsf(p,l)       get integer from string variable beginning at p     */
/*                      and field width l.                                  */
/*                      If error, return -1.                                */

double st_getsf(TDAContext *ctx, char *p,int l)
{
    (void)ctx;        /* unused: the signature is shared */
    register int n = 0;
    register char *q;
    double tmp = 0.0;
    double mul = 1.0;

    q = p + l;
    while (l-- > 0) {
        if (*--q == ' ')
            break;
        if (isdigit((int)*q)) {
            n = (int)(*q - 0x30);
            tmp += mul * (double)n; 
            n++;
        }
        else
            return(-1.0);
        mul *= 10.0;
    }
    if (n)
        return(tmp); 
    else
        return(-1.0);
}

/* ------------------------------------------------------------------------ */
/*  st_gets(p)            return short integer beginning at p.              */

int st_gets(TDAContext *ctx, char *p)
{
    short int ns;
    char buf[2];

    if ((ctx->HILO == 1 && ARCHTyp == 2) || (ctx->HILO == 2 && ARCHTyp == 1)) {
        buf[1] = *p++;
        buf[0] = *p;
    }
    else {
        buf[0] = *p++;
        buf[1] = *p;
    }
    memcpy(&ns,buf,sizeof(ns));   /* a char array carries no alignment: copy the bytes */
    return((int)ns);
}

/* ------------------------------------------------------------------------ */
/*  st_geti(p)          return integer beginning at p.                      */

int st_geti(TDAContext *ctx, char *p)
{
    int n;
    char buf[4];
        
    if ((ctx->HILO == 1 && ARCHTyp == 2) || (ctx->HILO == 2 && ARCHTyp == 1)) {
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
    /* buf is a char array and carries no alignment, so the bytes are
       copied into the int rather than read through a cast pointer */
    memcpy(&n,buf,sizeof(n));
    return(n);
}

/* ------------------------------------------------------------------------ */
/*  st_getf(p)          return float beginning at p.                        */

double st_getf(TDAContext *ctx, char *p)
{
    float x;
    char buf[4];
    unsigned char *q;
         
    q = (unsigned char *)p;
    if (ctx->HILO == 1) {
        if (*q++ == 0x7f)
            return(ctx->PMMSYS);
        /**************
        if (*q++ == 0x7f && *q++ == 0x00 && *q++ == 0x00 && *q++ == 0x00)
            return(PMMSYS);
        **********************/
    }
    else if (ctx->HILO == 2) {
        if (*(q + 3) == 0x7f)
            return(ctx->PMMSYS);
        /************
        if (*q++ == 0x00 && *q++ == 0x00 && *q++ == 0x00 && *q++ == 0x7f)
            return(PMMSYS);
        ***********************/
    }
    if ((ctx->HILO == 1 && ARCHTyp == 2) || (ctx->HILO == 2 && ARCHTyp == 1)) {
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
    memcpy(&x,buf,sizeof(x));   /* a char array carries no alignment: copy the bytes */
    return((double)x);
}

/* ------------------------------------------------------------------------ */
/*  st_getd(p)         return double beginning at p.                        */

double st_getd(TDAContext *ctx, char *p)
{
    double x;
    char buf[8];
    unsigned char *q;
         
    q = (unsigned char *)p;
    if (ctx->HILO == 1) {
        if (*q++ == 0x7f && *q++ == 0xe0)
            return(ctx->PMMSYS);
        /********* 
        if (*q++ == 0x7f && *q++ == 0xe0 && *q++ == 0x00 && *q++ == 0x00 &&
            *q++ == 0x00 && *q++ == 0x00 && *q++ == 0x00 && *q++ == 0x00)
            return(PMMSYS);
        ********************/
    }
    else if (ctx->HILO == 2) {
        if (*(q + 7) == 0x7f && *(q + 6) == 0xe0)
            return(ctx->PMMSYS);
        /************
        if (*q++ == 0x00 && *q++ == 0x00 && *q++ == 0x00 && *q++ == 0x00 &&
            *q++ == 0x00 && *q++ == 0x00 && *q++ == 0xe0 && *q++ == 0x7f)
            return(PMMSYS);
        **************/
    }
    if ((ctx->HILO == 1 && ARCHTyp == 2) || (ctx->HILO == 2 && ARCHTyp == 1)) {
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
    memcpy(&x,buf,sizeof(x));   /* see st_geti(): buf is not aligned */
    return((double)x);
}

/* ------------------------------------------------------------------------ */
/*  wr_stata(idx)                                                           */
/*  ##                                                                      */
/*          wstata(                                                         */
/*              keep = varlist,                                             */
/*              drop = varlist,                                             */
/*              sort = varlist,                                             */
/*              ptyp = ...,      type of Stata release                      */
/*                               4, 6, 7 or 10 (default is 10)                   */
/*              n = ...,     1 default                                      */
/*                           2 translate var names to lower case            */
/*                           3 translate var names to upper case            */
/*          ) = fname;                                                      */
/*                                                                          */
/*  Write current data matrix to a stata output file.                       */
/*  Return 0 if OK, -1 if error.                                            */

int wr_stata(TDAContext *ctx)  
{
    register int i,j,k,ii;
    register char *p,*q;
    unsigned char uc;
    char c,buf[120];
    int w,d,nv,sflag,err,l,hlen,vlen,vlablen,elen,fmtlen;
    short int ns;
    float f;                    /* the Stata float field: written as 4
                                   raw bytes below, so its width is the
                                   file format, not a precision choice */
    double x;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 6,13,1)) {   /* get parameters */
        p_clean(ctx);
        return(-1);   
    }       
    printf1(ctx, "Writing Stata file: %s\n",ctx->PMFdName);
    if (ctx->PMPTyp == 4) {
        hlen = 50;
        vlen = 9;
        vlablen = 32;
        fmtlen = 12;
        elen = 3;
    }
    else if (ctx->PMPTyp == 6) {
        ctx->PMPTyp = 6;
        hlen = 99;
        vlen = 9;
        vlablen = 81;
        fmtlen = 12;
        elen = 3;
    }
    else if (ctx->PMPTyp == 7) {
        hlen = 99;
        vlen = 33;
        vlablen = 81;
        fmtlen = 12;
        elen = 5;
    }
    else {
        ctx->PMPTyp = 10;
        hlen = 99;
        vlen = 33;
        vlablen = 81;
        fmtlen = 49;
        elen = 5;
    }

    printf1(ctx, "For Stata release: %d\n",ctx->PMPTyp);
   
    if (ctx->PMKeep && ctx->PMDrop) {
        p_err(ctx, -16,1);
        goto WSTATFin;
    }

    /* get list of variables in AcI[] */

    if (alloc_aci(ctx, imax(ctx, ctx->NVAR,ctx->PMNV)))
        goto WSTATFin;

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
        goto WSTATFin;
    }
    sflag = 0;
    if (ctx->PM1NV > 0) {        /* sort */
        err = vsort(ctx, ctx->PM1NV,ctx->PM1VIdx,1,0,1);
        if (err)
            goto WSTATFin;
        sflag = 1;
    }

    if (ctx->SILENTFlg < 2)
        printfe(ctx, "Writing: %s\n",ctx->PMFdName);

    if (ctx->PMPTyp == 4)
        buf[0] = 0x69;      
    else if (ctx->PMPTyp == 6)
        buf[0] = 0x6c;      
    else if (ctx->PMPTyp == 7)
        buf[0] = 0x6e;      
    else
        buf[0] = 0x72;      
       
    if (ARCHTyp == 1)
        buf[1] = 0x01;
    else
        buf[1] = 0x02;
    buf[2] = 0x01;
    buf[3] = 0x01; 
    fwrite(buf,4,1,ctx->PMFd);
    ns = (short)nv;
    fwrite((char *)&ns,2,1,ctx->PMFd);   /* number of variables */                        
    fwrite((char *)&ctx->NOC,4,1,ctx->PMFd);  /* number of cases */
  
    for (i = 0; i < hlen; ++i)
        fprintf(ctx->PMFd,"%c",0x00);

    /* write storage size */

    for (k = 0; k < nv; ++k) {
        j = ctx->AcI[k];

        switch (ctx->VSLen[j]) {         /* storage size */
            case  0:
            case  1:    if (ctx->PMPTyp == 10)
                            uc = 0xfb;
                        else
                            uc = 'b';
                        break;
            case  2:    if (ctx->PMPTyp == 10)
                            uc = 0xfc;
                        else
                            uc = 'i';
                        break;
            case  5:    if (ctx->PMPTyp == 10)
                            uc = 0xfd;
                        else
                            uc = 'l';
                        break;
            case  4:    if (ctx->PMPTyp == 10)
                            uc = 0xfe;
                        else
                            uc = 'f';
                        break;
            case  8:    if (ctx->PMPTyp == 10)
                            uc = 0xff;
                        else
                            uc = 'd';
                        break;
            default:    l = -ctx->VSLen[j];
                        if (l < 1 || l > 128) {
                            printf1(ctx, "Error: exceeded max string length (128 bytes).\n");
                            goto WSTATFin;
                        }
                        if (ctx->PMPTyp == 10)
                            uc = (unsigned char)l;
                        else
                            uc = (unsigned char)l + 0x7f;
        }
        fprintf(ctx->PMFd,"%c",uc);
    }

    /* write variable names */

    buf[vlen - 1] = '\0';
    for (k = 0; k < nv; ++k) {
        j = ctx->AcI[k];
        /* the field is written whole, so the tail past the name has to
           be NUL and not whatever the stack held (valgrind) */
        memset(buf,0,(size_t)vlen);
        if (ctx->PMN == 2 || ctx->PMN == 3) {
            p = ctx->VName[j];
            q = buf;
            ii = (int)(imin(ctx, vlen - 1,(int)(strlen(ctx->VName[j]))));
            for (i = 0; i < ii; ++i) {           
                if (ctx->PMN == 2)
                    *q++ = (char)tolower((int)*p++);
                else
                    *q++ = (char)toupper((int)*p++);
            }
            if (ii < vlen)
                *q = '\0';
        }
        else
            strncpy(buf,ctx->VName[j],(size_t)(vlen - 1));    
        fwrite(buf,(size_t)(vlen),1,ctx->PMFd);   
    }

    /* srtlist will be empty */

    j = 2 * (nv + 1);
    for (i = 0; i < j; ++i)
        fprintf(ctx->PMFd,"%c",0x00);

    /* format list */

    buf[11] = '\0';

    for (k = 0; k < nv; ++k) {
        j = ctx->AcI[k];
        memset(buf,0,(size_t)fmtlen);           /* as above */
/* ## */
        if (ctx->VTyp[j] == 1) {
            snprintf(buf,sizeof(buf),"%%%ds",-ctx->VSLen[j]);
            fwrite(buf,(size_t)(fmtlen),1,ctx->PMFd);
        }
        else {
            w = (int)ctx->VPFmt1[j];
            d = (int)ctx->VPFmt2[j];
            if (w == 0)
                snprintf(buf,sizeof(buf),"%%8.0g");
            else if (w > 0)
                snprintf(buf,sizeof(buf),"%%%d.%df",w,d);
            else             
                snprintf(buf,sizeof(buf),"%%%d.%de",-w,d);
            fwrite(buf,(size_t)(fmtlen),1,ctx->PMFd);   
        }
    }

    /* lbllist will be empty */

    j = vlen * nv;
    for (i = 0; i < j; ++i)
        fprintf(ctx->PMFd,"%c",0x00);

    /* write variable labels */

    for (k = 0; k < vlablen; ++k)
        buf[k] = '\0';

    for (k = 0; k < nv; ++k) {
        j = ctx->AcI[k];
        memset(buf,0,(size_t)vlablen);          /* as above */

        if (ctx->VLabel[j] != NULL)  
            strncpy(buf,ctx->VLabel[j],(size_t)(vlablen - 1));
        else
            buf[0] = '\0';

        fwrite(buf,(size_t)(vlablen),1,ctx->PMFd);   
    }

    /* expansion field will be empty */

    for (i = 0; i < elen; ++i)
        fprintf(ctx->PMFd,"%c",0x00);

    /* write data */

    for (i = 0; i < ctx->NOC; ++i) {
        ii = i;
        if (sflag)
            ii = ctx->VSORTPtr[i];

        for (k = 0; k < nv; ++k) {
            j = ctx->AcI[k];

            x = get_data(ctx, j,ii);           
                     
            switch (ctx->VSLen[j]) {         /* storage size */
                case  0:
                case  1:    c = (char)x;
                            fprintf(ctx->PMFd,"%c",c);
                            break;
                case  2:    ns = (short)x;
                            fwrite((char *)&ns,2,1,ctx->PMFd);
                            break;
                case  5:    w = (int)x;
                            fwrite((char *)&w,4,1,ctx->PMFd);
                            break;
                case  4:    f = (float)x;
                            fwrite((char *)&f,4,1,ctx->PMFd);
                            break;
                case  8:    fwrite((char *)&x,8,1,ctx->PMFd);
                            break;
                default:    l = -ctx->VSLen[j];
                            fwrite(ctx->VDPtr[j] + ii * l,(size_t)(l),1,ctx->PMFd);

            }
        }
        prn_message(ctx, i+1,0,1);
    }
    prn_message(ctx, ctx->NOC,1,1);
    printf1(ctx, "%d records with %d variables written to: %s\n",ctx->NOC,nv,ctx->PMFdName);
    err = 0;

WSTATFin:
    if (ctx->PM1NV > 0)                     
        vsort(ctx, 0,ctx->PM1VIdx,0,0,1);
    p_clean(ctx);
    return(err);
}



