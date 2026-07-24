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

/*  functions in t_stata.c */

int rd_stata(void); 
int st_gets(char *p);
int st_geti(char *p);
double st_getsf(char *p,int l); 
double st_getf(char *p);
double st_getd(char *p);
int wr_stata(void);

/*--------------------------------------------------------------------------*/
/*  Global variables                                                        */

int HILO = 1;       /* 1 if hilo byte order, 2 if lohi byte order           */

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

int rd_stata(void)  
{
    register int i,j,k;
    register char *p,*q,*vl;
    char buf[20000],*typlist,*varlist,*srtlist,*fmtlist,*lbllist;
    char *varlab,*rbuf;
    unsigned char u;
    int err,idxn,nvar,nrec,nrec1,nrec2,n,rlen,nmiss,nmiss1,w,d,hlen,fmtlen,typ;
    int len,typlista,varlista,srtlista,fmtlista,lbllista,varlaba,rbufa,nskip;
    int vn,m,nn,np,off,off1,mxnoc,aflag,vlen,vlablen,fptr;
    double x;

    if (check_cmd(0))
        return(-1);

    if (parm(CmdBuf + 6,10,1)) {   /* get parameters */
        p_clean();
        return(-1);   
    }   
    printf1("Reading Stata file: %s\n",PMFdName);
    if (DMDef) {
        printf1("Error: a data matrix already exists.\n");
        p_clean();
        return(-1);
    }
    aflag = 0;
    if (PMARCFDef && PMDVARFDef && PMF1Def)
        aflag = 1;

    if (PMNOCFlg)                     
        printf1("Maximum number of cases: %d\n",PMNOC);

    if (PMRHSTRA && PMF1Def == 0)
        printf1("Record selection (vsel): %s\n",PMRHSTR + 5);

    idxn = get_nidx();      /* index to first new variable */

    rbufa = varlaba = lbllista = fmtlista = srtlista = varlista = typlista = 0;
    err = -1;
           
    if (fread(buf,sizeof(char),10,PMFd) != 10) {
        p_err(-7,1);
        goto RSTATFin;
    }
    p = buf;
    printf1("Release: 0x%02x ",*p);
    if (*p == 0x69) {
        typ = 4;
        hlen = 50;
        fmtlen = 12;
        vlen = 9;
        vlablen = 32;
    }
    else if (*p == 0x68 || *p == 0x66) {
        typ = 3;
        hlen = 32;
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
        printf1("Can't identify version.\n");
        goto RSTATFin;
    }
    printf1("[version %d]\n",typ);

    if (*++p == 0x01)
        HILO = 1;
    else if (*p == 0x02)
        HILO = 2;
    else {
        printf1("Can't identify HiLo indicator. Found: 0x%02x\n",(unsigned char)*p);
        goto RSTATFin;
    }
    if (*++p != 0x01) {
        printf1("Probably not a stata data file (found: 0x%02x)\n",(unsigned char)*p);
        goto RSTATFin;
    }
    nvar = st_gets(buf + 4);
    nrec  = st_geti(buf + 6);

    printf1("Number of variables: %d\n",nvar);
    printf1("Number of cases: %d\n",nrec);

    if (nvar < 1 || nrec < 1) {
        printf1("Probably an error.\n");
        goto RSTATFin;
    }
    if (nvar > MaxNV) {
        printf1("Error: exceeded max number of variables.\n");
        goto RSTATFin;
    }
    if (fread(buf,sizeof(char),hlen,PMFd) != hlen) {
        p_err(-7,1);
        goto RSTATFin;
    }
    printf1("Data label: %s\n",buf);
    if (typ == 10) {
        printf1("Time stamp: %s\n",buf + 81);
    }
    newline();

    if (!(typlist = (char *) calloc(nvar + 1,sizeof(char)))) {
        p_err(-2,1);
        goto RSTATFin;
    }
    typlista = nvar + 1;
    memrq(typlista,1);
         
    if (!(varlist = (char *) calloc(vlen * nvar + 1,sizeof(char)))) {
        p_err(-2,1);
        goto RSTATFin;
    }
    varlista = vlen * nvar + 1;
    memrq(varlista,1);

    if (!(srtlist = (char *) calloc(2 * (nvar + 1) + 1,sizeof(char)))) {
        p_err(-2,1);
        goto RSTATFin;
    }
    srtlista = 2 * (nvar + 1) + 1;
    memrq(srtlista,1);

    if (!(fmtlist = (char *) calloc(fmtlen * nvar + 1,sizeof(char)))) {
        p_err(-2,1);
        goto RSTATFin;
    }
    fmtlista = fmtlen * nvar + 1;
    memrq(fmtlista,1);

    if (!(lbllist = (char *) calloc(vlen * nvar + 1,sizeof(char)))) {
        p_err(-2,1);
        goto RSTATFin;
    }
    lbllista = vlen * nvar + 1;
    memrq(lbllista,1);

    if (!(varlab = (char *) calloc(vlablen * nvar + 1,sizeof(char)))) {
        p_err(-2,1);   
        goto RSTATFin;
    }
    varlaba = vlablen * nvar + 1;
    memrq(varlaba,1);
      
    if (fread(typlist,sizeof(char),nvar,PMFd) != nvar) {
        p_err(-7,1);
        goto RSTATFin;
    }
    if (fread(varlist,sizeof(char),vlen * nvar,PMFd) != vlen * nvar) {
        p_err(-7,1);
        goto RSTATFin;
    }
    if (fread(srtlist,sizeof(char),2 * (nvar + 1),PMFd) != 2 * (nvar + 1)) {
        p_err(-7,1);
        goto RSTATFin;
    }
    if (fread(fmtlist,sizeof(char),fmtlen * nvar,PMFd) != fmtlen * nvar) {
        p_err(-7,1);
        goto RSTATFin;
    }
    if (fread(lbllist,sizeof(char),vlen * nvar,PMFd) != vlen * nvar) {
        p_err(-7,1);
        goto RSTATFin;
    }
    if (fread(varlab,sizeof(char),vlablen * nvar,PMFd) != vlablen * nvar) {
        p_err(-7,1);
        goto RSTATFin;
    }
    rlen = 0;   /* record length */

    if (alloc_acj(nvar))        /* used for storage size */
        goto RSTATFin;
    if (alloc_acr(nvar))        /* used for format w */
        goto RSTATFin;
    if (alloc_acs(nvar))        /* used for format d */
        goto RSTATFin;
         
    for (i = 0; i < nvar; ++i) {

        /********************************  
         printf1("%d %c :%s:%s:%s:%s:\n",
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
                q = skip_int(q);
                q = skip_int(++q);
                if (*q == 'e')
                    w = -w;
                else if (*q == 'g')     /* free format */
                    w = d = 0;
            }
            len = iabs(w);
        }
        sprintf(p,"%s<%d>[%d.%d]",varlist + i * vlen,n,w,d);
    
        if (PMN == 2) {
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
            sprintf(p,"(%s)",varlab + i * vlablen);
        }
        p = buf + strlen(buf);

        if (n < 0)
            sprintf(p,"=stata(%d)",len);     /* string variable */
        else if (typ == 10)
            sprintf(p,"=stata(0x%02x)",(unsigned char)typlist[i]);
        else                 
            sprintf(p,"=stata(%c)",typlist[i]);

        AcJ[i] = n;
        AcR[i] = w;
        AcS[i] = d;
                            
        if (save_var(buf,0)) {       /* save variable definition */

            printf1("\nError: can't save variable definitions.\n");
            goto RSTATFin;              
        }
    }
    prn_var(idxn);      /* print list of new variables */
    printf1("\n");

    /* skip expansion fields */
      
    while (typ >= 4) {      /* only for version 4 and higher */
         
        if (typ <= 6) {
            if (fread(buf,sizeof(char),3,PMFd) != 3) {
                p_err(-7,1);
                goto RSTATFin;
            }
            n = st_gets(buf + 1);
        }
        else {
            if (fread(buf,sizeof(char),5,PMFd) != 5) {
                p_err(-7,1);
                goto RSTATFin;
            }
            n = st_geti(buf + 1);
        }
        if (n == 0 && !*buf)
            break;

        if (fread(buf,sizeof(char),n,PMFd) != n) {
            p_err(-7,1);
            goto RSTATFin;
        }
    }
              
    /* if vsel */

    if (PMRHSTRA) {
        if (PMF1Def)  
            printf1("Warning: vsel option will be ignored.\n");
        else {
            if ((n = v_parse(PMRHSTR + 5,0)) < 0 || ESCnt <= 0) {
                printf1("Syntax error (%d) in vsel expression.\n",n);
                if (n < 0)
                    prn_emsg1(n);
                goto RSTATFin;
            }
        }
    }

    /* read data to create internal data matrix, or directly writing
       to an output file. */

    if (rlen < 1) {
        printf1("Error in record length.\n");
        goto RSTATFin;
    }
    if (!(rbuf = (char *) calloc(rlen + 1,sizeof(char)))) {
        p_err(-2,1);
        goto RSTATFin;
    }
    rbufa = rlen + 1;
    memrq(rbufa,1);

    if (PMNOCFlg)
        mxnoc = imin(PMNOC,nrec);
    else 
        mxnoc = nrec;   

    if (PMF1Def) {
        printf1("Data will be directly written to output file: %s\n",PMF1dName);
    }
    else {
        NOCMaxA = mxnoc;
        printf1("Reading data to create internal data matrix.\n");
        printf1("Maximum number of cases: %d\n",mxnoc);
        printf1("Record length: %d\n",rlen);

        /* allocate memory for data */
        
        n = MemReq;
        if (alloc_vdat(idxn,1)) {
            printf1("Error: insufficient memory for data matrix (%d cases).\n",NOCMaxA);
            goto RSTATFin;
        }
        printf1("Allocated %d bytes for data matrix.\n\n",MemReq - n);
    }
    nrec2 = nrec1 = 0;      /* number of records */
    nmiss = 0;              /* number of system missing values */
    nskip = 0;              /* skipped by vsel */
           
    if (SILENTFlg < 2)
        printfe("Reading: %s\n",PMFdName);

    for (nrec2 = 0; nrec2 < mxnoc; ++nrec2) {

        n = fread(rbuf,sizeof(char),rlen,PMFd);             
        if (n == 0)
            break;

        if (n != rlen) {
            printf1("Warning: break with reading %d instead of %d bytes.\n",n,rlen);   
            break;          
        }

        q = p = rbuf;
        nmiss1 = 0;

        for (j = 0; j < nvar; ++j) {
            w = 0;
          
            if (typ == 10) {
                if  ((unsigned char)typlist[j] == 0xfb) {
                    if (*p >= 0x65) {
                        x = PMMSYS;
                        nmiss1++;
                    }
                    else
                        x = (double)*p;
                    p += 1;
                }
                else if  ((unsigned char)typlist[j] == 0xfc) {
                    n = st_gets(p);
                    if (n >= 0x7fe5) {
                        x = PMMSYS;
                        nmiss1++;
                    }
                    else
                        x = (double)n;
                        p += 2;
                }
                else if  ((unsigned char)typlist[j] == 0xfd) {
                    n = st_geti(p);
                    if (n == 0x7fffffe5) {
                        x = PMMSYS;
                        nmiss1++;
                    }
                    else
                        x = (double)n;
                    p += 4;
                }
                else if  ((unsigned char)typlist[j] == 0xfe) {
                    x = st_getf(p);
                    if (x == PMMSYS)
                        nmiss1++;
                    p += 4;
                }   
                else if  ((unsigned char)typlist[j] == 0xff) {
                    x = st_getd(p);
                    if (x == PMMSYS)
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
                                    x = PMMSYS;
                                    nmiss1++;
                                }
                                else
                                    x = (double)*p;
                                p += 1;
                                break;
                    case 'i':   n = st_gets(p);
                                if (n == 0x7fff) {
                                    x = PMMSYS;
                                    nmiss1++;
                                }
                                else
                                    x = (double)n;
                                p += 2;
                                break;
                    case 'l':   n = st_geti(p);
                                if (n == 0x7fffffff) {
                                    x = PMMSYS;
                                    nmiss1++;
                                }
                                else
                                    x = (double)n;
                                p += 4;
                                break;
                    case 'f':   x = st_getf(p);
                                if (x == PMMSYS)
                                    nmiss1++;
                                p += 4;
                                break;
                    case 'd':   x = st_getd(p);
                                if (x == PMMSYS)
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

                vl = SVBuf;
                for (i = 0; i < w; ++i) {
                    if (*q) 
                        *vl++ = *q++;
                    else
                        *vl++ = ' ';
                }
                if (PMF1Def) {
                    fwrite(SVBuf,w,1,PMF1d);
                    fprintf(PMF1d," ");
                }
                else  
                    put_str(SVBuf,w,j,nrec1,0);        
            }
            else {
                if (PMF1Def)
                    fprintf(PMF1d,VPFmtS[j],x);
                else
                    put_data(x,j,nrec1);
            }
        }
        if (PMF1Def)
            fprintf(PMF1d,"\n");

        prn_message(nrec2,0,0);

        if (PMRHSTRA && PMF1Def == 0) {         /* vsel */

            n = v_eval1(nrec1,ESCnt,ESTyp,ESVal,ESIdx,&x,0,0,0,0,0);
            if (n) {
                printf1("Can't evaluate vsel expression in record %d.\n",nrec2 + 1);
                prn_emsg2(n);
                goto RSTATFin;
            }
            if (fabs(x) < EPSI2) {
                nskip++;
                continue;     
            }
        }
        nrec1++;
        nmiss += nmiss1;
    }
    prn_message(nrec2,1,0);
    printf1("Read %d records. ",nrec2);
    if (PMRHSTRA && PMF1Def == 0)         
        printf1("Selected: %d.",nrec1);
    printf1("\n");

    if (PMF1Def) {
        printf1("%d records written to: %s\n",nrec1,PMF1dName);
        if (nrec1 <= 0)  
            goto RSTATFin;
    }
    else {
        if (nrec1 <= 0) {
            printf1("Error: no data matrix created.\n");
            goto RSTATFin;
        }
        NOCDM = NOC = nrec1;
        DMDef = 1;
    
        printf1("Created a data matrix with %d variables and %d cases.\n",NVAR,NOC);
    }
    if (nrec1 < nrec)  
        printf1("Warning: file contains more than %d records.\n",nrec1);
    printf1("Number of system missing values (substituted by %lg): %d\n",
                                                  PMMSYS,nmiss);
     
    /* try to read value labels */

    if (PMDVARFDef) {

        fprintf(PMDVARFd,"# Variable description file based on: %s\n",PMFdName);

        if (typ >= 6 && PMNOCFlg == 0) {    /* only if no limit to number of cases */

            if (alloc_aci(nvar))            /* used for file pointer */
                goto RSTATFin;

            nn = 10000;
            if (alloc_acn(nn + 1))          /* used for off[] */
                goto RSTATFin;
            if (alloc_acm(nn + 1))          /* used for val[] */
                goto RSTATFin;

            np = 10000;
            if (alloc_acc(np + 2))          /* used for read buffer */
                goto RSTATFin;

            while (1) {
                fptr = ftell(PMFd);
                if (fread(buf,sizeof(char),vlen + 7,PMFd) != vlen + 7)  
                    break;
      
                len = st_geti(buf);
                    if (len < 1) {
                    p_err(-44,1);
                    break;
                }
                p = buf + 4;
                        
                for (i = 0; i < nvar; ++i) {
                    if (!strcmp(p,lbllist + i * vlen))  
                        AcI[i] = fptr;
                }
                if (np < len) {
                    np = len;
                    if (alloc_acc(np + 2))      /* used for read buffer */
                        goto RSTATFin;
                }
                if (fread(AcC,sizeof(char),len,PMFd) != len) {
                    p_err(-44,1);
                    goto RSTATFin;
                }
                n = st_geti(AcC);
                m = st_geti(AcC + 4);
                if (n < 1 || m < 1) {
                    p_err(-44,1);
                    break;
                }
                if (nn < n) {
                    nn = n;
                    if (alloc_acn(nn + 1))      /* used for off[] */
                        goto RSTATFin;
                    if (alloc_acm(nn + 1))      /* used for val[] */
                        goto RSTATFin;
                }
            }
        }
        vn = off = 0;
        for (i = 0; i < nvar; ++i) {

            if (PMN == 2) {
                q = varlist + i * vlen;
                while (*q) {
                    fprintf(PMDVARFd,"%c",(char)toupper((int)*q));
                    q++;
                }
                fprintf(PMDVARFd," ");
            }
            else
                fprintf(PMDVARFd,"%s ",varlist + i * vlen);

            fprnchar(PMDVARFd,' ',12 - strlen(varlist + i * vlen),0);
            fprintf(PMDVARFd,"%3d %5d ",PMDVARFN,off);

            n = AcJ[i];
            w = AcR[i];
            d = AcS[i];
            if (n < 0)
                len = -n;
            else
                len = iabs(w);

            if (n < 0)
                fprintf(PMDVARFd,"%4d     ",n);
            else
                fprintf(PMDVARFd,"%4d.%-2d  ",w,d);

            vl = varlab + i * vlablen;       
            for (j = 0; j < vlablen; ++j) {
                if (!*vl)  
                    break;
                fprintf(PMDVARFd,"%c",*vl++);
            }   
            fprintf(PMDVARFd,"\n");
            vn++;
            off1 = off;
            off += len + 1;

            if (typ >= 6 && PMNOCFlg == 0 && AcI[i] != 0) {      /* read value labels */

                fseek(PMFd,AcI[i],0);

                if (fread(buf,sizeof(char),vlen + 7,PMFd) == vlen + 7) {

                    len = st_geti(buf);
                    if (len > 0) {
                        if (fread(AcC,sizeof(char),len,PMFd) == len) {
                            nn = st_geti(AcC);
                            m = st_geti(AcC + 4);
                            if (nn > 0 && m > 0) {
                                p = AcC + 8;
                                for (k = 0; k < nn; ++k) {
                                    AcN[k] = st_geti(p);
                                    p += 4;
                                }
                                for (k = 0; k < nn; ++k) {
                                    AcM[k] = st_geti(p);
                                    p += 4;
                                }
                                for (k = 0; k < nn; ++k)  
                                    fprintf(PMDVARFd," %6d  %s\n",AcM[k],p + AcN[k]);
                            }
                        }
                    }
                }
            }

            if (n < 0) {     /* check for partition */
                                  
                nn = get_dvarp(&np,varlist + i * vlen); 

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
                        fprintf(PMDVARFd,"%s%02d ",varlist + i * vlen,j);
                        fprnchar(PMDVARFd,' ',10 - strlen(varlist + i * vlen),0);
                        fprintf(PMDVARFd,"%3d %5d ",PMDVARFN,off1);

                        if (np == 1)
                            fprintf(PMDVARFd,"%4d     ",-nn);
                        else
                            fprintf(PMDVARFd,"%4d     ",nn);

                        vl = varlab + i * vlablen;       
                        if (*vl) {
                            for (k = 0; k < vlablen; ++k) {
                                if (!*vl)  
                                    break;
                                fprintf(PMDVARFd,"%c",*vl++);
                            }   
                            fprintf(PMDVARFd,", part %02d",j);
                        }
                        fprintf(PMDVARFd,"\n");
                        vn++;
                        off1 += nn;
                        m += nn;
                        j++;
                    }
                }
            }
        }
        printf1("Variable description file written to: %s\n",PMDVARFName);
    }

    if (aflag)          /* archive description file */
        make_arcd(PMDVARFN,PMF1dName,off,nrec1,vn);

    err = 0;

RSTATFin:
    if (typlista) {
        free(typlist);
        memrq(-typlista,1);
    }
    if (varlista) {
        free(varlist);
        memrq(-varlista,1);
    }
    if (srtlista) {
        free(srtlist);
        memrq(-srtlista,1);
    }
    if (fmtlista) {
        free(fmtlist);
        memrq(-fmtlista,1);
    }
    if (lbllista) {
        free(lbllist);
        memrq(-lbllista,1);
    }
    if (varlaba) {
        free(varlab);
        memrq(-varlaba,1);
    }
    if (rbufa) {
        free(rbuf);
        memrq(-rbufa,1);
    }
    if (err || PMF1Def) {
        clear_avar(idxn);
        NVAR = NOCMaxA = NOCDM = NOC = 0;
    }
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  st_getsf(p,l)       get integer from string variable beginning at p     */
/*                      and field width l.                                  */
/*                      If error, return -1.                                */

double st_getsf(char *p,int l)
{
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

int st_gets(char *p)
{
    short int *ns;
    char buf[2];

    if ((HILO == 1 && ARCHTyp == 2) || (HILO == 2 && ARCHTyp == 1)) {
        buf[1] = *p++;
        buf[0] = *p;
    }
    else {
        buf[0] = *p++;
        buf[1] = *p;
    }
    ns = (short int*)buf;
    return((int)*ns);
}

/* ------------------------------------------------------------------------ */
/*  st_geti(p)          return integer beginning at p.                      */

int st_geti(char *p)
{
    int *n;
    char buf[4];
        
    if ((HILO == 1 && ARCHTyp == 2) || (HILO == 2 && ARCHTyp == 1)) {
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
    n = (int *)buf;
    return(*n);
}

/* ------------------------------------------------------------------------ */
/*  st_getf(p)          return float beginning at p.                        */

double st_getf(char *p)
{
    float *x;
    char buf[4];
    unsigned char *q;
         
    q = (unsigned char *)p;
    if (HILO == 1) {
        if (*q++ == 0x7f)
            return(PMMSYS);
        /**************
        if (*q++ == 0x7f && *q++ == 0x00 && *q++ == 0x00 && *q++ == 0x00)
            return(PMMSYS);
        **********************/
    }
    else if (HILO == 2) {
        if (*(q + 3) == 0x7f)
            return(PMMSYS);
        /************
        if (*q++ == 0x00 && *q++ == 0x00 && *q++ == 0x00 && *q++ == 0x7f)
            return(PMMSYS);
        ***********************/
    }
    if ((HILO == 1 && ARCHTyp == 2) || (HILO == 2 && ARCHTyp == 1)) {
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
    x = (float *)buf;
    return(*x);
}

/* ------------------------------------------------------------------------ */
/*  st_getd(p)         return double beginning at p.                        */

double st_getd(char *p)
{
    double *x;
    char buf[8];
    unsigned char *q;
         
    q = (unsigned char *)p;
    if (HILO == 1) {
        if (*q++ == 0x7f && *q++ == 0xe0)
            return(PMMSYS);
        /********* 
        if (*q++ == 0x7f && *q++ == 0xe0 && *q++ == 0x00 && *q++ == 0x00 &&
            *q++ == 0x00 && *q++ == 0x00 && *q++ == 0x00 && *q++ == 0x00)
            return(PMMSYS);
        ********************/
    }
    else if (HILO == 2) {
        if (*(q + 7) == 0x7f && *(q + 6) == 0xe0)
            return(PMMSYS);
        /************
        if (*q++ == 0x00 && *q++ == 0x00 && *q++ == 0x00 && *q++ == 0x00 &&
            *q++ == 0x00 && *q++ == 0x00 && *q++ == 0xe0 && *q++ == 0x7f)
            return(PMMSYS);
        **************/
    }
    if ((HILO == 1 && ARCHTyp == 2) || (HILO == 2 && ARCHTyp == 1)) {
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
    x = (double *)buf;
    return(*x);
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

int wr_stata(void)  
{
    register int i,j,k,ii;
    register char *p,*q;
    unsigned char uc;
    char c,buf[120];
    int w,d,nv,sflag,err,l,hlen,vlen,vlablen,elen,fmtlen;
    short int ns;
    float f;
    double x;

    err = -1;
    if (check_cmd(0))
        return(-1);

    if (parm(CmdBuf + 6,13,1)) {   /* get parameters */
        p_clean();
        return(-1);   
    }       
    printf1("Writing Stata file: %s\n",PMFdName);
    if (PMPTyp == 4) {
        hlen = 50;
        vlen = 9;
        vlablen = 32;
        fmtlen = 12;
        elen = 3;
    }
    else if (PMPTyp == 6) {
        PMPTyp = 6;
        hlen = 99;
        vlen = 9;
        vlablen = 81;
        fmtlen = 12;
        elen = 3;
    }
    else if (PMPTyp == 7) {
        hlen = 99;
        vlen = 33;
        vlablen = 81;
        fmtlen = 12;
        elen = 5;
    }
    else {
        PMPTyp = 10;
        hlen = 99;
        vlen = 33;
        vlablen = 81;
        fmtlen = 49;
        elen = 5;
    }

    printf1("For Stata release: %d\n",PMPTyp);
   
    if (PMKeep && PMDrop) {
        p_err(-16,1);
        goto WSTATFin;
    }

    /* get list of variables in AcI[] */

    if (alloc_aci(imax(NVAR,PMNV)))
        goto WSTATFin;

    nv = 0;
    if (PMKeep) {
        for (k = 0; k < PMNV; ++k)  
            AcI[nv++] = PMVIdx[k];
    }
    else {
        j = VIFirst;
        while (j >= 0) {
            i = 1;
            if (PMDrop) {
                for (k = 0; k < PMNV; ++k) {
                    if (PMVIdx[k] == j) {
                        i = 0;
                        break;
                    }   
                }
            }
            if (i) 
                AcI[nv++] = j;
            j = VNxt[j];
        }
    }
    if (nv == 0) {
        printf1("No variables selected.\n"); 
        goto WSTATFin;
    }
    sflag = 0;
    if (PM1NV > 0) {        /* sort */
        err = vsort(PM1NV,PM1VIdx,1,0,1);
        if (err)
            goto WSTATFin;
        sflag = 1;
    }

    if (SILENTFlg < 2)
        printfe("Writing: %s\n",PMFdName);

    if (PMPTyp == 4)
        buf[0] = 0x69;      
    else if (PMPTyp == 6)
        buf[0] = 0x6c;      
    else if (PMPTyp == 7)
        buf[0] = 0x6e;      
    else
        buf[0] = 0x72;      
       
    if (ARCHTyp == 1)
        buf[1] = 0x01;
    else
        buf[1] = 0x02;
    buf[2] = 0x01;
    buf[3] = 0x01; 
    fwrite(buf,4,1,PMFd);
    ns = (short)nv;
    fwrite((char *)&ns,2,1,PMFd);   /* number of variables */                        
    fwrite((char *)&NOC,4,1,PMFd);  /* number of cases */
  
    for (i = 0; i < hlen; ++i)
        fprintf(PMFd,"%c",0x00);

    /* write storage size */

    for (k = 0; k < nv; ++k) {
        j = AcI[k];

        switch (VSLen[j]) {         /* storage size */
            case  0:
            case  1:    if (PMPTyp == 10)
                            uc = 0xfb;
                        else
                            uc = 'b';
                        break;
            case  2:    if (PMPTyp == 10)
                            uc = 0xfc;
                        else
                            uc = 'i';
                        break;
            case  5:    if (PMPTyp == 10)
                            uc = 0xfd;
                        else
                            uc = 'l';
                        break;
            case  4:    if (PMPTyp == 10)
                            uc = 0xfe;
                        else
                            uc = 'f';
                        break;
            case  8:    if (PMPTyp == 10)
                            uc = 0xff;
                        else
                            uc = 'd';
                        break;
            default:    l = -VSLen[j];
                        if (l < 1 || l > 128) {
                            printf1("Error: exceeded max string length (128 bytes).\n");
                            goto WSTATFin;
                        }
                        if (PMPTyp == 10)
                            uc = (unsigned char)l;
                        else
                            uc = (unsigned char)l + 0x7f;
        }
        fprintf(PMFd,"%c",uc);
    }

    /* write variable names */

    buf[vlen - 1] = '\0';
    for (k = 0; k < nv; ++k) {
        j = AcI[k];
        if (PMN == 2 || PMN == 3) {
            p = VName[j];
            q = buf;
            ii = imin(vlen - 1,strlen(VName[j]));
            for (i = 0; i < ii; ++i) {           
                if (PMN == 2)
                    *q++ = (char)tolower((int)*p++);
                else
                    *q++ = (char)toupper((int)*p++);
            }
            if (ii < vlen)
                *q = '\0';
        }
        else
            strncpy(buf,VName[j],vlen - 1);    
        fwrite(buf,vlen,1,PMFd);   
    }

    /* srtlist will be empty */

    j = 2 * (nv + 1);
    for (i = 0; i < j; ++i)
        fprintf(PMFd,"%c",0x00);

    /* format list */

    buf[11] = '\0';

    for (k = 0; k < nv; ++k) {
        j = AcI[k];
/* ## */
        if (VTyp[j] == 1) {
            sprintf(buf,"%%%ds",-VSLen[j]);
            fwrite(buf,fmtlen,1,PMFd);
        }
        else {
            w = (int)VPFmt1[j];
            d = (int)VPFmt2[j];
            if (w == 0)
                sprintf(buf,"%%8.0g");
            else if (w > 0)
                sprintf(buf,"%%%d.%df",w,d);
            else             
                sprintf(buf,"%%%d.%de",-w,d);
            fwrite(buf,fmtlen,1,PMFd);   
        }
    }

    /* lbllist will be empty */

    j = vlen * nv;
    for (i = 0; i < j; ++i)
        fprintf(PMFd,"%c",0x00);

    /* write variable labels */

    for (k = 0; k < vlablen; ++k)
        buf[k] = '\0';

    for (k = 0; k < nv; ++k) {
        j = AcI[k];

        if (VLabel[j] != NULL)  
            strncpy(buf,VLabel[j],vlablen - 1);
        else
            buf[0] = '\0';

        fwrite(buf,vlablen,1,PMFd);   
    }

    /* expansion field will be empty */

    for (i = 0; i < elen; ++i)
        fprintf(PMFd,"%c",0x00);

    /* write data */

    for (i = 0; i < NOC; ++i) {
        ii = i;
        if (sflag)
            ii = VSORTPtr[i];

        for (k = 0; k < nv; ++k) {
            j = AcI[k];

            x = get_data(j,ii);           
                     
            switch (VSLen[j]) {         /* storage size */
                case  0:
                case  1:    c = (char)x;
                            fprintf(PMFd,"%c",c);
                            break;
                case  2:    ns = (short)x;
                            fwrite((char *)&ns,2,1,PMFd);
                            break;
                case  5:    w = (int)x;
                            fwrite((char *)&w,4,1,PMFd);
                            break;
                case  4:    f = (float)x;
                            fwrite((char *)&f,4,1,PMFd);
                            break;
                case  8:    fwrite((char *)&x,8,1,PMFd);
                            break;
                default:    l = -VSLen[j];
                            fwrite(VDPtr[j] + ii * l,l,1,PMFd);

            }
        }
        prn_message(i+1,0,1);
    }
    prn_message(NOC,1,1);
    printf1("%d records with %d variables written to: %s\n",NOC,nv,PMFdName);
    err = 0;

WSTATFin:
    if (PM1NV > 0)                     
        vsort(0,PM1VIdx,0,0,1);
    p_clean();
    return(err);
}



