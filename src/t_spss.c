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

/*  functions in t_spss.c */

int rd_spss(void); 
int buf_update(int init);
void prstr(char *p, int n);
char *getstr(int f);
int getdigit(char *p, int *err);
char *inum(char *p, int *n);
char *dnum(char *p, double *x, int *mv);
void digit_err(char *p);
int wr_spss(void);
void pnum1(FILE *fd,int n);
void pnum(FILE *fd,int n);
void pfnum(FILE *fd,double x);
void pstring(FILE *fd,char *s,int m);
void check_nl(FILE *fd);
int rd_spss1(void); 
int rd_dat(int n,char *buf);
int sav_sval(int *sval);
int sav_geti(int arch,int *err);
double sav_getd(int arch,int *err);
int wr_spss1(void);
void sav_puti(int n);
void sav_putd(double x);
int get_dvarp(int *np,char *vname);
void make_arcd(int fn,char *fname,int len,int nrec,int vn);

/*--------------------------------------------------------------------------*/
/*  Global variables                                                        */

#define RLEN 80             /* Record length                                */
#define HLEN 200            /* Header length                                */
#define TLEN 256            /* Length of translation table                  */
#define DTLEN 10            /* Length of creation date string               */
#define TTLEN  8            /* Length of creation time string               */

int RSRILen = 0;            /* Input file record length                     */
int RSPEOF = 0;             /* set if eof in input file                     */
char *RSPBuf;               /* Read Buffer                                  */
int RSPBufA = 0;
char *RSPPtr;               /* Global pointer to read buffer                */
int RSPCnt = 0;             /* Number of char's in read buffer              */
int RSPMaxD = 10;           /* max digits                                   */

long RSPFPtr = 0;
int RSPFOff = 0;
long RSPFPtr1 = 0;
int RSPFOff1 = 0;
long RSPFPtr2 = 0;
int RSPFOff2 = 0;

long *RSPVFPtr; 
int *RSPVFOff; 

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

int rd_spss(void)  
{
    register int i,j,k;
    register char *p;
    int err,idxn,nvar,nvar1,nvar2,m,n,vtyp,len,prec,mv,nmiss,nrec;
    int aptr,aoff,avint,avmax,avmin,avttyp,vlflag,pfmt1,pfmt2,nskip,nmiss1;
    int np,nrec1,nrec2,mxnoc,aflag,vrec,vn,len1;
    short *vint,*vttyp;
    char vname[VNLMax + 1],vlabel[VLLMax + 1],vdef[400];
    double tol,a,x,y,*vmax,*vmin;

    tol = 10000.0 * EPSI;
    aptr = aoff = avint = avmin = avmax = avttyp = 0;
    idxn = get_nidx();      /* index to first new variable */
    VLabelLen = 0;
    err = -1;

    if (check_cmd(0))
        return(-1);

    if (parm(CmdBuf + 5,10,1)) {   /* get parameters */
        p_clean();
        return(-1);   
    }   

    aflag = 0;
    if (PMARCFDef && PMDVARFDef && PMF1Def)
        aflag = 1;

    printf1("Reading SPSS export file: %s\n",PMFdName);
    if (DMDef) {
        printf1("Error: a data matrix already exists.\n");
        p_clean();
        return(-1);
    }
    if (PMLEN >= 80 && PMLEN <= 82) {
        RSRILen = PMLEN;
        printf1("Assuming fixed record length: %d bytes.\n",RSRILen);
    }
    else  
        RSRILen = 0;

    if (PMNOCFlg)
        printf1("Maximum number of cases: %d\n",PMNOC);

    if (PMRHSTRA && PMF1Def == 0)
        printf1("Record selection (vsel): %s\n",PMRHSTR + 5);

    /* allocate read buffer */

    if (!(RSPBuf = (char *) calloc(2 * RLEN + 6,sizeof(char)))) {
        p_err(-2,1);
        goto RSPFin;
    }
    memrq(2 * RLEN + 6,1);
    RSPBufA = 2 * RLEN + 6;
             
    RSPPtr = RSPBuf;
    RSPEOF = RSPCnt = 0;
    if (buf_update(1)) 
        goto RSPFin;
           
    /*  Get the 5 * 40 bytes header */

    printf1("\nHeader: ");          
    for (i = 1; i <= 5; ++i) {
        prstr(RSPPtr,40);
        RSPPtr += 40;
        RSPCnt -= 40;
        printf1("\n        ");
        if (buf_update(0))
            goto RSPFin;
    }

    /*  Skip translation table, length is HLEN = 256 chars */
  
    for (i = 1; i <= 3; ++i) {
        RSPPtr += RLEN;
        RSPCnt -= RLEN;
        if (buf_update(0))
            goto RSPFin;
    }
    RSPPtr += 16;
    RSPCnt -= 16;
    if (buf_update(0))
        goto RSPFin;
           
    /*  Read the SPSSPORT string, length is 8 characters */
   
    if (strncmp(RSPPtr,"SPSSPORT",8)) {
        printf1("\nWarning: can't locate SPSSPORT string.\n");
        /*******************
        prstr(RSPPtr,40);                                       
        printf1("\n");
        goto RSPFin;
        *******************/
    }
    RSPPtr += 8;
    RSPCnt -= 8;

    /*  Read the version code, length is 1 char */
  
    if (buf_update(0))
        goto RSPFin;
    printf1("\nFile-format version code: %c <%02x>",*RSPPtr,*RSPPtr);
    RSPPtr++;
    RSPCnt--;
  
    /*  Read creation date (10 char) and time */

    printf1("\nCreation date: ");
    prstr(RSPPtr + 2,DTLEN - 2);
    RSPPtr += DTLEN;
    RSPCnt -= DTLEN;

    printf1("\nCreation time: ");
    prstr(RSPPtr + 2,TTLEN - 2);
    RSPPtr += TTLEN;
    RSPCnt -= TTLEN;

    if (buf_update(0))
        goto RSPFin;

    if (*RSPPtr == '1') {        /* get originating software */
        printf1("\nSoftware: ");
        RSPCnt--;
        RSPPtr++;
        RSPPtr = getstr(0);
        if (RSPPtr == NULL)
            goto RSPFin;
    }
    if (*RSPPtr == '2') {        /* get originating installation */
        if (buf_update(0))
            goto RSPFin;
        printf1("\nInstallation: ");
        RSPCnt--;
        RSPPtr++;
        RSPPtr = getstr(0);
        if (RSPPtr == NULL)
            goto RSPFin;
    }
    if (*RSPPtr == '3') {        /* get file label */
        if (buf_update(0))
            goto RSPFin;
        printf1("\nLabel: ");
        RSPCnt--;
        RSPPtr++;
        RSPPtr = getstr(0);
        if (RSPPtr == NULL)
            goto RSPFin;
    }
    printf1("\n");           

    /*  get the number of variables */
   
    if (*RSPPtr != '4') {
        printf1("\nError: can't find number of variables.\n");
        goto RSPFin;
    }
    if (buf_update(0))
        goto RSPFin;

    RSPCnt--;
    RSPPtr = inum(++RSPPtr,&nvar);
    if (RSPPtr == NULL)
        goto RSPFin;

    printf1("\nNumber of variables: %d",nvar);
    if (nvar <= 0) {
        printf1("\nError in number of variables.\n");
        goto RSPFin;
    }
    else if (nvar > MaxNV) {
        printf1("\nError: exceeded max number of variables.\n");
        goto RSPFin;
    }
    if (*RSPPtr == '5') {        /* get number of base-30 digits */
        RSPCnt--;
        RSPPtr = inum(++RSPPtr,&n);
        if (RSPPtr == NULL)
            goto RSPFin;
        printf1("\nPrecison (base-30 digits): %d",n);
        if (n > 0)
            RSPMaxD = n;
    }

    printf1("\nCase-weight variable: ");          

    if (*RSPPtr == '6') {        /* get case-weight variable */
        RSPCnt--;
        RSPPtr++;
        RSPPtr = getstr(0);
        if (RSPPtr == NULL)
            goto RSPFin;
    }
    else  
        printf1("not defined.");
    printf1("\n");           

    if (!(vttyp = (short *) calloc(nvar,sizeof(short)))) {
        p_err(-2,1);
        goto RSPFin;
    }
    avttyp = nvar;
    memrq(nvar,sizeof(short));

    /*  Read the variable descriptions */

    nvar1 = 0;      /* number of string variables */
    nvar2 = 0;      /* number of variables with unknown format */

    for (i = 0; i < nvar; ++i) {

        if (buf_update(0))
            goto RSPFin;

        if (*RSPPtr != '7') {
            printf1("\nError: can't find description of variable %d.\nFound: ",i);
            prstr(RSPPtr - 10,40);
            printf1("\n");
            goto RSPFin;
        }

        /* get type of variable: 0 numerical, else string */

        RSPCnt--;
        RSPPtr = inum(++RSPPtr,&vtyp);
        if (RSPPtr == NULL)
            goto RSPFin;

        vttyp[i] = vtyp; 

        if (vtyp > 0)       /* alphanumerical variable */
            nvar1++;

        /* get name of variable */
  
        RSPPtr = inum(RSPPtr,&n);
        if (RSPPtr == NULL)
            goto RSPFin;

        p = vname; 
        for (j = 0; j < n; ++j) {
            if (j < VNLMax)  
                *p++ = get_vnchar(*RSPPtr);
            RSPPtr++;
        }
        *p = '\0';
        RSPCnt -= n;

        /*  get print formats (not used) */

        for (j = 1; j <= 3; ++j) {
            RSPPtr = inum(RSPPtr,&n);
            if (RSPPtr == NULL)
                goto RSPFin;
        }

        /* get write formats; ignored */
   
        RSPPtr = inum(RSPPtr,&n);
        if (RSPPtr == NULL)
            goto RSPFin;

        if (n != 1 && n != 5)
            nvar2++;
           
        RSPPtr = inum(RSPPtr,&len);
        if (RSPPtr == NULL)
            goto RSPFin;

        RSPPtr = inum(RSPPtr,&prec);
        if (RSPPtr == NULL)
            goto RSPFin;

        /*  get missing value codes */
  
        n = 0;
        for (j = 0; j < 3; ++j) {

            switch (*RSPPtr) {

                case '8':   RSPCnt--;
                            RSPPtr = dnum(++RSPPtr,&x,&mv);
                            n++;
                            break;
                case '9':   RSPCnt--;
                            RSPPtr = dnum(++RSPPtr,&x,&mv);
                            n++;
                            break;
                case 'A':   RSPCnt--;
                            RSPPtr = dnum(++RSPPtr,&x,&mv);
                            n++;
                            break;
                case 'B':   RSPCnt--;
                            RSPPtr = dnum(++RSPPtr,&x,&mv);
                            if (RSPPtr == NULL)
                                break;
                            RSPPtr = dnum(RSPPtr,&y,&mv);
                            n++;
                            break;
                default:    goto CONT;
            }
            if (RSPPtr == NULL)
                goto RSPFin;
        }
CONT:
        sprintf(vdef,"%s<4>[0.0]",vname);

        vlflag = 0;
        if (*RSPPtr == 'C') {   /* Var label is optional */

            RSPCnt--;
            RSPPtr = inum(++RSPPtr,&n);
            if (RSPPtr == NULL)
                goto RSPFin;

            if (buf_update(0))
                goto RSPFin;

            p = vlabel;
            for (j = 0; j < n; ++j) {
    
                if (RSPCnt < 10) {
                    if (buf_update(0))
                        goto RSPFin;
                }
                if (j < VLLMax)
                    *p++ = *RSPPtr;
                RSPPtr++;
                RSPCnt--;
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
            sprintf(p,"(%s)",vlabel);
            vlflag = 1;
        }
        p = vdef + strlen(vdef);
        sprintf(p,"=spss(%d)",vtyp);
        if (save_var(vdef,0)) {       /* save variable definition */
            printf1("\nError: can't save variable definitions.\n");
            goto RSPFin;
        }
    }
    printf1("\nNumber of string variables: %d",nvar1);
    printf1("\nNumber of variables with unknown format: %d\n",nvar2);
          
    /*  Read value labels (optional). */

    vlflag = 0;

    if (*RSPPtr == 'D') {

        printf1("\nReading value labels.\n");
               
        /* allocate memory for file pointers */

        if (!(RSPVFPtr = (long *) calloc(nvar,sizeof(long)))) {
            p_err(-2,1);
            goto RSPFin;
        }
        aptr = nvar;
        memrq(nvar,sizeof(long));

        if (!(RSPVFOff = (int *) calloc(nvar,sizeof(int)))) {
            p_err(-2,1);
            goto RSPFin;
        }
        aoff = nvar;
        memrq(nvar,sizeof(int));

        vlflag = 1;

        while (*RSPPtr == 'D') {

            if (buf_update(0))
                goto RSPFin;
            RSPPtr++;
            RSPCnt--;
                          
            RSPPtr = inum(RSPPtr,&n);
            if (RSPPtr == NULL)
                goto RSPFin;

            len = 0;
            for (j = 1; j <= n; ++j) {
                if (buf_update(0))
                    goto RSPFin;

                RSPPtr = inum(RSPPtr,&m);
                if (RSPPtr == NULL)
                    goto RSPFin;

                p = vname;
                for (i = 0; i < m; ++i)  
                    *p++ = *RSPPtr++;
                *p = '\0';
                RSPCnt -= m;

                k = -1;
                for (i = 0; i < nvar; ++i) {
                    if (!strcmp(vname,VName[i])) {
                        k = i;
                        break;
                    }
                }
                if (k >= 0) {
                    RSPVFOff[k] = -1;
                    len = VTyp[k];
                }
            }
            if (buf_update(0))
                goto RSPFin;

            for (i = 0; i < nvar; ++i) {
                if (RSPVFOff[i] < 0) {
                    RSPVFOff[i] = RSPFOff;
                    RSPVFPtr[i] = RSPFPtr;
                }
            }
            RSPPtr = inum(RSPPtr,&n);
            if (RSPPtr == NULL)
                goto RSPFin;

            for (j = 1; j <= n; ++j) {
                if (buf_update(0))
                    goto RSPFin;

                if (len != 1) {       /* numerical variable */
                    RSPPtr = dnum(RSPPtr,&x,&mv);
                    if (RSPPtr == NULL)  
                        goto RSPFin;
                }
                else {              /* string variable */
                    RSPPtr = inum(RSPPtr,&m);
                    if (RSPPtr == NULL)
                        goto RSPFin;

                    RSPPtr += m;
                    RSPCnt -= m;
                }
                RSPPtr = inum(RSPPtr,&m);
                if (RSPPtr == NULL)
                    goto RSPFin;

                RSPPtr += m;
                RSPCnt -= m;
            }
        }
    }

    /* Read textual information documenting the file (optional) */
  
    if (*RSPPtr == 'E') {
        printf1("\nFile documentation\n");
        RSPCnt--;
        RSPPtr = inum(++RSPPtr,&n);
        if (RSPPtr == NULL)
            goto RSPFin;

        for (j = 1; j <= n; ++j) {
            if (buf_update(0))
                goto RSPFin;
            if ((RSPPtr = getstr(0)) == NULL)  
                goto RSPFin;
            printf1("\n");          
        }
        printf1("\n"); 
    }

    if (buf_update(0))
        goto RSPFin;

    RSPFPtr1 = RSPFPtr;     /* remember file position for data */
    RSPFOff1 = RSPFOff;

    if (*RSPPtr != 'F') {               /* check for data stream */

        printf1("Error: can't find data stream.\n");
        goto RSPFin;
    }
    RSPCnt--;
    RSPPtr++;
   
    /* allocate memory to check for maximum and integer values */
   
    if (!(vmax = (double *) calloc(nvar,sizeof(double)))) {
        p_err(-2,1);
        goto RSPFin;
    }
    avmax = nvar;
    memrq(nvar,sizeof(double));
   
    if (!(vmin = (double *) calloc(nvar,sizeof(double)))) {
        p_err(-2,1);
        goto RSPFin;
    }
    avmin = nvar;
    memrq(nvar,sizeof(double));

    if (!(vint = (short *) calloc(nvar,sizeof(short)))) {
        p_err(-2,1);
        goto RSPFin;
    }
    avint = nvar;
    memrq(nvar,sizeof(short));
   
    printf1("Reading data to check variables.\n");

    if (SILENTFlg < 2)
        printfe("Reading: %s\n",PMFdName);

    nrec = 0;       /* number of records */
    nmiss = 0;      /* number of system missing values */
                     
    while (*RSPPtr && *RSPPtr != 'Z') {    /* this is the end of file marker */

        for (i = 0; i < nvar; ++i) {
            if (buf_update(0))
                goto RSPFin;

            if (vttyp[i] > 0) {        /* string variable */
                RSPPtr = inum(RSPPtr,&m);
                if (RSPPtr == NULL)
                    goto RSPFin;
                if (m > vttyp[i])  
                    printf1("Warning: variable %d in record %d (type=%d fnd=%d).\n",
                                            i+1,nrec+1,vttyp[i],m);

                while (m >= RSPCnt) {
                    m -= RSPCnt;
                    RSPPtr += RSPCnt;
                    RSPCnt = 0;
                    if (buf_update(0))
                        goto RSPFin;
                }
                RSPPtr += m;          
                RSPCnt -= m;
            }
            else {
                RSPPtr = dnum(RSPPtr,&x,&mv);
                if (RSPPtr == NULL)
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
        prn_message(nrec,0,0);
   
        if (PMNOCFlg && nrec >= PMNOC)
            break;
    
        if (RSPEOF > RLEN)
            break;
    }
    prn_message(nrec,1,0);

    printf1("Read %d records.\n",nrec);
    if (nrec <= 0) {
        printf1("Error.\n");
        goto RSPFin;
    }
    if (*RSPPtr != 'Z')  
        printf1("Warning: can't find end of file character (Z).\n");
    printf1("Number of system missing values: %d\n\n",nmiss);

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
                VSLen[i] = 8;
            else if ((int)a < 128)
                VSLen[i] = 1;
            else if ((int)a < 32000)
                VSLen[i] = 2;
            else 
                VSLen[i] = 5;

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
                if (PMFmtF) {
                    pfmt1 = PMFmt1;
                    pfmt2 = PMFmt2;
                }
                else {
                    pfmt1 += 7;
                    pfmt2 = 6;
                }
            }
        }
        makefmt(&pfmt1,&pfmt2,VPFmtS[i],0,SEPC,0);
        VPFmt1[i] = (short)pfmt1;
        VPFmt2[i] = (short)pfmt2;
    }
    prn_var(idxn);      /* print list of new variables */
    newline();     

    /* if requested create the variable description file */

    len = 0;
    if (PMDVARFDef) {
  
        printf1("Creating variable description file: %s\n",PMDVARFName);
            
        fprintf(PMDVARFd,"# variable description file based on: %s\n",PMFdName);

        vrec = 2,
        len1 = vn = 0;

        for (i = 0; i < nvar; ++i) {

            fprintf(PMDVARFd,"%s ",VName[i]);
            fprnchar(PMDVARFd,' ',VNameLen - strlen(VName[i]),0);
            fprintf(PMDVARFd," %3d %4d ",PMDVARFN,len);
            if (VTyp[i] != 1) {
                fprintf(PMDVARFd,"%4d",(int)VPFmt1[i]);
                if ((int)VPFmt2[i] > 0)
                    fprintf(PMDVARFd,".%-2d ",(int)VPFmt2[i]);
                else
                    fprintf(PMDVARFd,"    ");
            }
            else  
                fprintf(PMDVARFd,"%4d    ",-vttyp[i]);

            if (VLabel[i] != NULL)
                fprintf(PMDVARFd,"%s",VLabel[i]);
            fprintf(PMDVARFd,"\n");
            vrec++;
            vn++;
            if (VTyp[i] != 1)
                len += (int)VPFmt1[i] + 1;
            else {
                len1 = len;
                len += vttyp[i] + 1;
            }

            /* add value labels */

            if (vlflag && RSPVFPtr[i] > 0) {
                   
                RSPEOF = 0;
                fseek(PMFd,RSPVFPtr[i],0);
                RSPPtr = RSPBuf;
                RSPCnt = 0;
                if (buf_update(0))
                    goto RSPFin;
                RSPPtr += RSPVFOff[i];
                RSPCnt -= RSPVFOff[i];
                if (buf_update(0))
                    goto RSPFin;
                RSPPtr = inum(RSPPtr,&m);
                if (RSPPtr == NULL)
                    goto RSPFin;

                for (j = 1; j <= m; ++j) {
                    if (buf_update(0))
                        goto RSPFin;

                    RSPPtr = dnum(RSPPtr,&x,&mv);
                    if (RSPPtr == NULL)  
                        goto RSPFin;

                    fprintf(PMDVARFd," %4g ",x);

                    RSPPtr = inum(RSPPtr,&n);
                    if (RSPPtr == NULL)
                        goto RSPFin;

                    for (k = 0; k < n; ++k)  
                        fprintf(PMDVARFd,"%c",*RSPPtr++);
                    fprintf(PMDVARFd,"\n");
                    RSPCnt -= n;
                    vrec++;
                }
            }
            if (VTyp[i] == 1) {     /* check for partition */
                                  
                n = get_dvarp(&np,VName[i]); 

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
                        fprintf(PMDVARFd,"%s%02d ",VName[i],j);
                        fprnchar(PMDVARFd,' ',VNameLen - strlen(VName[i]) - 2,0);
                        fprintf(PMDVARFd," %3d %4d ",PMDVARFN,len1);
                        if (np == 1)
                            fprintf(PMDVARFd,"%4d    ",-n);
                        else
                            fprintf(PMDVARFd,"%4d    ",n);
                        if (VLabel[i] != NULL)
                            fprintf(PMDVARFd,"%s, part %02d",VLabel[i],j);
                        fprintf(PMDVARFd,"\n");
                        vn++;
                        len1 += n;
                        m += n;
                        j++;
                    }
                }
            }
        }
        printf1("%d records (%d variables) written to: %s\n\n",vrec,vn,PMDVARFName);
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
                goto RSPFin;
            }
        }
    }
   
    /* reading data again to create internal data matrix, or directly
       writing to output file. */
    
    if (PMF1Def) {
        printf1("Data will be directly written to output file: %s\n",PMF1dName);
        if (PMNOCFlg)
            mxnoc = PMNOC;
        else 
            mxnoc = INTMAX;
    }
    else {
        if (NOCMaxA <= 0)
            NOCMaxA = nrec;
        mxnoc = NOCMaxA;
       
        printf1("Reading data again to create internal data matrix.\n");
        printf1("Maximum number of cases: %d\n",NOCMaxA);
    }
    if (nmiss > 0)  
        printf1("System missing values will be substituted by: %g\n",PMMSYS);

    RSPEOF = 0;                 /* seek to begin of data */
    fseek(PMFd,RSPFPtr1,0);
    RSPPtr = RSPBuf;
    RSPCnt = 0;
    if (buf_update(0))
        goto RSPFin;
    RSPPtr += RSPFOff1;
    RSPCnt -= RSPFOff1;
    if (buf_update(0))
        goto RSPFin;

    if (*RSPPtr != 'F') {               /* check for data stream */

        printf1("Error: can't find data stream.\n");
        goto RSPFin;
    }
    RSPCnt--;
    RSPPtr++;

    /* allocate memory for data */
        
    if (PMF1Def == 0) {
        n = MemReq;
        if (alloc_vdat(idxn,1)) {
            printf1("Error: insufficient memory for data matrix (%d cases).\n",NOCMaxA);
            goto RSPFin;
        }
        printf1("Allocated %d bytes for data matrix.\n\n",MemReq - n);
    }
    nrec2 = nrec1 = 0;      /* number of records */
    nmiss = 0;              /* number of system missing values */
    nskip = 0;              /* skipped by vsel */

    while (*RSPPtr && *RSPPtr != 'Z') {    /* this is the end of file marker */

        if (nrec2 >= mxnoc)  
            break;
            
        nmiss1 = 0;
        for (i = 0; i < nvar; ++i) {
            if (buf_update(0))
                goto RSPFin;

            if (vttyp[i] > 0) {        /* string variable */
                RSPPtr = inum(RSPPtr,&m);
                if (RSPPtr == NULL)
                    goto RSPFin;
                if (m > vttyp[i])  
                    printf1("Warning: variable %d in record %d (type=%d fnd=%d).\n",
                                            i+1,nrec1+1,vttyp[i],m);
          
                p = VDPtr[i] - nrec1 * VSLen[i];                  
                for (k = 0; k < m; ++k) {

                    if (RSPCnt < 10) {
                        if (buf_update(0))
                            goto RSPFin;
                    }
                    if (PMF1Def)
                        fprintf(PMF1d,"%c",*RSPPtr);
                    else
                        *p++ = *RSPPtr;

                    RSPPtr++;
                    RSPCnt--;
                }
                for (; k < vttyp[i]; ++k) {
                    if (PMF1Def)
                        fprintf(PMF1d," ");
                    else
                        *p++ = ' ';
                }
                if (PMF1Def)
                    fprintf(PMF1d," ");
            }
            else {
                RSPPtr = dnum(RSPPtr,&x,&mv);
                if (RSPPtr == NULL)
                    goto RSPFin;
    
                if (mv) {
                    nmiss1++;    
                    x = PMMSYS;
                }
                else if (vint[i] == 0) {
                    y = (int)x;
                    x = (int)y;
                }
                if (PMF1Def)
                    fprintf(PMF1d,VPFmtS[i],x);
                else
                    put_data(x,i,nrec1);
     
            }
        }
        if (PMF1Def)
            fprintf(PMF1d,"\n");

        nrec2++;
        prn_message(nrec2,0,0);
        if (RSPEOF > RLEN)
            break;

        if (PMRHSTRA && PMF1Def == 0) {         /* vsel */

            n = v_eval1(nrec1,ESCnt,ESTyp,ESVal,ESIdx,&x,0,0,0,0,0);
            if (n) {
                printf1("Can't evaluate vsel expression in record %d.\n",nrec2 + 1);
                prn_emsg2(n);
                goto RSPFin;
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
        printf1("Selected: %d",nrec1);
    printf1("\n");

    if (PMF1Def) {
        printf1("%d records written to: %s\n",nrec1,PMF1dName);
        if (nrec1 <= 0)  
            goto RSPFin;
        if (len > 0)
            printf1("Record length without EOL characters: %d bytes.\n",len);
    }
    else {
        if (nrec1 <= 0) {
            printf1("Error: no data matrix created.\n");
            goto RSPFin;
        }
        NOCDM = NOC = nrec1;
        DMDef = 1;
       
        printf1("Created a data matrix with %d variables and %d cases.\n",NVAR,NOC);
        if (nrec1 < nrec &&  PMRHSTRA == 0)
            printf1("Warning: file contains more than %d records.\n",nrec1);
    }
    printf1("Number of system missing values: %d\n",nmiss);

    if (aflag)          /* archive description file */
        make_arcd(PMDVARFN,PMF1dName,len,nrec1,vn);

    err = 0;

RSPFin:
    if (avttyp) {
        free((char *)vttyp);
        memrq(-avttyp,sizeof(short));
    }
    if (aptr) {
        free((char *)RSPVFPtr);
        memrq(-aptr,sizeof(long));
    }
    if (aoff) {
        free((char *)RSPVFOff);
        memrq(-aoff,sizeof(int));
    }
    if (avmax) {
        free((char *)vmax);
        memrq(-avmax,sizeof(double));
    }
    if (avmin) {
        free((char *)vmin);
        memrq(-avmin,sizeof(double));
    }
    if (avint) {
        free((char *)vint);
        memrq(-avint,sizeof(short));
    }
    if (RSPBufA > 0) {
        free(RSPBuf);
        memrq(-RSPBufA,1);
    }
    if (err || PMF1Def) {
        if (idxn < get_nidx())
            clear_avar(idxn);
        NVAR = NOCMaxA = NOCDM = NOC = 0;
    }
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  buf_update()    Read a new record into RSPBuf                           */
/*                  Return  0 if OK, -1 if error.                           */

int buf_update(int init)
{
    register int i;
    register char *p,*q;
    static long fptr = 0L;

    if (init)
        fptr = 0L;

    if (RSPEOF || RSPCnt >= RLEN) {
        if (RSPEOF) {
            fptr = 0L;
            RSPEOF++;
        }
        if (RSPCnt == RLEN) {
            RSPFPtr = fptr;
            RSPFOff = 0;
        }
        else
            RSPFOff = 2 * RLEN - RSPCnt;
        return(0);
    }   
    RSPFOff = RLEN - RSPCnt;
    RSPFPtr = fptr;

    fptr = ftell(PMFd);

    p = RSPBuf;
    q = RSPPtr;
    for (i = 0; i < RSPCnt; ++i)
        *p++ = *q++;
    RSPPtr = RSPBuf;

    if (RSRILen == 0) {
        if (fgets(p,RLEN + 3,PMFd)) {
            q = p;
            for (i = 0; i < RLEN; ++i) {
                if (!*q || *q == CR || *q == LF)
                    break;
                q++;
            }
            for (; i < RLEN; ++i)  
                *q++ = ' ';
            *q = '\0';
            RSPCnt += RLEN;
        }
        else  
            RSPEOF = 1;
    }
    else {
        if ((i = fread(p,sizeof(*p),RSRILen,PMFd)) < 0) {
            printf1("\nError in reading the input file.\n");
            return(-1);
        }
        if (i < RLEN) {
            RSPCnt += i;
            RSPEOF = 1;
        }
        else  
            RSPCnt += RLEN;
    }
    if (RSPEOF)    
        *(RSPPtr + RSPCnt) = '\0';
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  prstr(p,n)      print a string at p with length n                       */

void prstr(char *p, int n)
{
    register int i;
    for (i = 0; i < n; ++i) {
        if (*p >= ' ' && *p <= 'z')
            printf1("%c",*p);
        else
            printf1("?");
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

char *getstr(int f)
{           
    register int i;
    int n;

    RSPPtr = inum(RSPPtr,&n);
    if (RSPPtr == NULL)  
        return(NULL);

    if (n > 0) {

        if (RSPCnt <= n) {
            if (buf_update(0))
                return(NULL);     
        }
        if (RSPCnt < n) {
            printf1("Error in reading SPSS file.\n");
            return(NULL);
        }
        if (f >= 0) {
            for (i = 0; i < n; ++i)  
                printf1("%c",*RSPPtr++);
            for (; i < f; ++i)  
                printf1(" ");
        }
        else
            RSPPtr += n;

        RSPCnt -= n;
    }
    return(RSPPtr);
}

/*--------------------------------------------------------------------------*/
/*  getdigit(p,err)      get and return a base30 digit at p.                */
/*                       in case of an error, return err = 1.               */

int getdigit(char *p, int *err)
{
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

char *inum(char *p, int *n)
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
        RSPCnt--;
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
                *n += getdigit(p,&err);
                if (err) {
                    digit_err(q);
                    return(NULL);
                }
                len++;
                if (len > 8)  
                    printf1("Warning: found integer entry with %2d (base-30) digits.\n",len);
            }
            else {
                nn *= 30;
                nn += getdigit(p,&err);
                if (err) {
                    digit_err(q);
                    return(NULL);
                }
            }
        }
        p++;
        RSPCnt--;
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
    RSPCnt--;
    return(++p);       
}

/*--------------------------------------------------------------------------*/
/*  dnum(p,x,mv)    read a double value at p and return it in x, return     */
/*                  pointer to next character. mv is set to 1 if a system   */
/*                  missing value code is found. If error, return NULL      */

char *dnum(char *p, double *x, int *mv)
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
        RSPCnt--;
        p++;
    }
    *mv = 0;
    if (*p == '*') {                /* check for internal missing value */
        *mv = 1;
        RSPCnt -= 2;
        p += 2;
        return(p);
    }
    if (*p == '-') {
        neg = 1;
        RSPCnt--;
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
                man += (double) getdigit(p,&err);
                if (err) {
                    digit_err(q);
                    return(NULL);
                }
                k++;
                if (k > 13)  
                    printf1("Warning: found entry with %2d (base-30) digits.\n",k);
            }
            else {
                mex *= 30.0;
                mex += (double) getdigit(p,&err);
                if (err) {
                    digit_err(q);
                    return(NULL);
                }
            }
        }
        RSPCnt--;
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
    RSPCnt--;
    return(++p);
}

/*--------------------------------------------------------------------------*/
/*  digit_err(p)        print error message                                 */

void digit_err(char *p)
{
    register int i;
    printf1("\nError: can't read a numerical value.");
    printf1("\nFound: ");
    for (i = 0; i < 40; ++i)  
        printf1("%c",*p++);
    printf1("\n");
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

int SPSSPtr = 21;       /* pointer for writing records with RLEN rec length */

int wr_spss(void)
{
    register int i,j,k,ii;
    int err,u,v,nv,sflag,l;
    double x;

    err = -1;
    if (check_cmd(0))
        return(-1);

    if (parm(CmdBuf + 5,1,1)) {   /* get parameters */
        p_clean();
        return(-1);   
    }       
    printf1("Writing SPSS export file: %s\n",PMFdName);

    if (PMKeep && PMDrop) {
        p_err(-16,1);
        goto WSPSSFin;
    }

    /* get list of variables in AcI[] */

    if (alloc_aci(imax(NVAR,PMNV)))
        goto WSPSSFin;

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
        goto WSPSSFin;
    }
    sflag = 0;
    if (PM1NV > 0) {        /* sort */
        err = vsort(PM1NV,PM1VIdx,1,0,1);
        if (err)
            goto WSPSSFin;
        sflag = 1;
    }

    /* print header etc. */

    fprintf(PMFd,"%s\n",S1);
    fprintf(PMFd,"%s\n",S2);
    fprintf(PMFd,"%s\n",S3);
    fprintf(PMFd,"%s\n",S4);
    fprintf(PMFd,"%s\n",S5);
    fprintf(PMFd,"%s\n",S6);
    fprintf(PMFd,"%s",S7);

    /* print number of variables */

    fprintf(PMFd,"4");
    check_nl(PMFd);
    pnum(PMFd,nv);

    /* digits (optional) */

    fprintf(PMFd,"5");
    check_nl(PMFd);
    pnum(PMFd,10);

    /* case weight variable name (optional) not written */

    /* dictionary entries (required) */

    for (k = 0; k < nv; ++k) {
        j = AcI[k];
        fprintf(PMFd,"7");
        check_nl(PMFd);
        if (VTyp[j] == 1)
            pnum(PMFd,-VSLen[j]);
        else
            pnum(PMFd,0);
        pstring(PMFd,VName[j],8);
        u = (int)VPFmt1[j];
        v = (int)VPFmt2[j];
        if (u == 0) {
            u = 10;
            v = 4;
        }
        else if (u < 0) {
            u = 12;
            v = 4;
        }
        pnum(PMFd,5);
        pnum(PMFd,u);
        pnum(PMFd,v);
        pnum(PMFd,5);
        pnum(PMFd,u);
        pnum(PMFd,v);

        if (VLabelLen > 0) {
            if (VLabel[j] != NULL) {
                fprintf(PMFd,"C");
                check_nl(PMFd);
                pstring(PMFd,VLabel[j],40);
            }
        }
    }

    /* print data stream */

    if (SILENTFlg < 2)
        printfe("Writing: %s\n",PMFdName);
   
    fprintf(PMFd,"F");
    check_nl(PMFd);

    for (i = 0; i < NOC; ++i) {
        ii = i;
        if (sflag)
            ii = VSORTPtr[i];

        for (k = 0; k < nv; ++k) {
            j = AcI[k];

            if (VTyp[j] == 1) {
                l = -VSLen[j];
                pstring(PMFd,VDPtr[j] + ii * l,l);
            }
            else {
                x = get_data(j,ii);           
                if (x < 0.0) {
                    fprintf(PMFd,"-");
                    check_nl(PMFd);
                    x = -x;
                }
                u = (int)x;
                if (x == (double)u)
                    pnum(PMFd,u);
                else  
                    pfnum(PMFd,x);
            }
        }
        prn_message(i + 1,0,1);
    }
    prn_message(NOC,1,1);
    fprintf(PMFd,"Z");            /* end of file */
    check_nl(PMFd);
    while (SPSSPtr++ < RLEN)  
        fprintf(PMFd,"Z");        
    fprintf(PMFd,"\n");

    printf1("%d records with %d variables written to: %s\n",NOC,nv,PMFdName);
    err = 0;

WSPSSFin:
    if (PM1NV > 0)                     
        vsort(0,PM1VIdx,0,0,1);
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pnum1(fd,n)     print integer n to fd (base 30). Update SPSSPtr.        */

char DIG30[] = {'0','1','2','3','4','5','6','7','8','9',
                'A','B','C','D','E','F','G','H','I','J',
                'K','L','M','N','O','P','Q','R','S','T'};

void pnum1(FILE *fd,int n)
{
    int m,r;
    register char *p;
    char buf[100];

    if (n < 0) {
        fprintf(fd,"-");
        check_nl(fd);
        n = -n;
    }
    p = buf;
    while (n >= 30) {
        m = (int) (n / 30);
        r = n - 30 * m;
        sprintf(p++,"%c",DIG30[r]);
        n = m;
    }
    sprintf(p++,"%c",DIG30[n]);
    while (p > buf) {
        fprintf(fd,"%c",*--p);
        check_nl(fd);
    }
}

/*--------------------------------------------------------------------------*/
/*  pnum(fd,n)      print integer n to fd (base 30). Update SPSSPtr.        */

void pnum(FILE *fd,int n)
{
    pnum1(fd,n);
    fprintf(fd,"/");
    check_nl(fd);
}

/*--------------------------------------------------------------------------*/
/*  pfnum(fd,x)     print floating point number x to fd (base 30).          */
/*                  Update SPSSPtr. Note: x >= 0.0                          */

void pfnum(FILE *fd,double x)
{
    register int i;
    double a,b,c,d,e;

    e = floor(log(x) / log(30.0));
    b = x / pow(30.0,e);
    c = floor(b);
    if (c < 0.0 || c >= 30.0)  
        gerr_exit(74);
         
    fprintf(fd,"%c",DIG30[(int)c]);
    check_nl(fd);
    b -= c;
    if (b > EPSI) {
        fprintf(fd,".");
        check_nl(fd);
        c = 30.0;
        for (i = 0; i < 10 ; ++i) {
            a = b * c;
            d = floor(a);
            fprintf(fd,"%c",DIG30[(int)d]);
            check_nl(fd);
            b -= d / c;
            if (b <= EPSI)  
                break;
            c *= 30;
        }
    }
    i = (int)e;
    if (i) {
        if (i < 0) {
            fprintf(fd,"-");
            check_nl(fd);
            i = -i;
        }
        else {
            fprintf(fd,"+");
            check_nl(fd);
        }
        pnum1(fd,i);
    }
    fprintf(fd,"/");
    check_nl(fd);
}

/*--------------------------------------------------------------------------*/
/*  pstring(fd,s,m)     print string s to fd. If m > 0, this is max length  */

void pstring(FILE *fd,char *s,int m)
{
    register int l;
    register char *p;

    p = s;
    l = strlen(p);
    if (m > 0 && l > m)  
        l = m;
    pnum(fd,l);
    while (l-- > 0) {
        fprintf(fd,"%c",*p++);
        check_nl(fd);
    }
}

/*--------------------------------------------------------------------------*/
/*  check_nl(fd)    check SPSSPtr. If >= RLEN print newline and set         */
/*                  SPSSPtr = 0;                                            */

void check_nl(FILE *fd)
{
    if (++SPSSPtr < 80)
        return;
    fprintf(fd,"\n");
    SPSSPtr = 0;
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

int rd_spss1(void)  
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
    tol = 10000.0 * EPSI;
    smiss = -(DBLMAX - 100.0);

    idxn = get_nidx();      /* index to first new variable */
    VLabelLen = 0;
    err = -1;

    if (check_cmd(0))
        return(-1);

    if (parm(CmdBuf + 6,10,1)) {   /* get parameters */
        p_clean();
        return(-1);   
    }   
    printf1("Reading SPSS sav file: %s\n",PMFdName);
    if (DMDef) {
        printf1("Error: a data matrix already exists.\n");
        p_clean();
        return(-1);
    }
    if (PMNOCFlg)
        printf1("Maximum number of cases: %d\n",PMNOC);
    newline();

    if ((r = rd_dat(64,buf)) != 64)
        goto RSP1Fin; 

    printf1("Identification: %s\n",buf);
    if (strncmp(buf,"$FL2",4)) {
        printf1("Probably not an SPSS sav file.\n");
        err = -1;
        goto RSP1Fin;
    }
    n = sav_geti(arch,&r);
    if (r)
        goto RSP1Fin;
        
    if (n != 2)         /* file layout code should be 2 */           
        arch = 2;                  

    nobs = sav_geti(arch,&r);
    if (r)
        goto RSP1Fin;
    printf1("Number of OBS elements per observation: %d\n",nobs);

    cflag = sav_geti(arch,&r);
    if (r)
        goto RSP1Fin;
    printf1("Compression switch: %d\n",cflag);

    widx = sav_geti(arch,&r);
    if (r)
        goto RSP1Fin;
    printf1("Index of case-weight variable: %d\n",widx);

    noc = sav_geti(arch,&r);
    if (r)
        goto RSP1Fin;
    printf1("Number of cases: %d\n",noc);

    bias = sav_getd(arch,&r);
    if (r)
        goto RSP1Fin;
    printf1("Compression bias: %g\n\n",bias);

    if ((r = rd_dat(9,buf)) != 9)
        goto RSP1Fin; 
    printf1("Creation date: %s\n",buf);

    if ((r = rd_dat(8,buf)) != 8)
        goto RSP1Fin; 
    printf1("Creation time: %s\n",buf);

    if ((r = rd_dat(67,buf)) != 67)
        goto RSP1Fin; 
    printf1("File label: %s\n\n",buf);

    if (nobs < 1) {
        printf1("Error: need at least one variable description.\n");
        goto RSP1Fin;
    }
    if (alloc_acn(nobs + 1))        /* type of variable */
        goto RSP1Fin; 

    if (alloc_acm(nobs + 1))        /* spss-internal variable number */
        goto RSP1Fin; 

    if (alloc_aci(nobs + 1))        /* integer flags */
        goto RSP1Fin; 

    nvar1 = nvar = 0;
    for (i = 0; i < nobs; ++i) {

        n = sav_geti(arch,&r);
        if (r)
            goto RSP1Fin;
   
        if (n != 2) {
            printf1("Error: can't read next dictionary entry.\n");
            goto RSP1Fin;
        }
        vtyp = sav_geti(arch,&r);
        if (r)
            goto RSP1Fin;

        AcN[i] = vtyp;

        if (vtyp >= 0)
            AcM[i] = nvar;
        else
            AcM[i] = -1;

        vlflag = sav_geti(arch,&r);
        if (r)
            goto RSP1Fin;

        nmiss = sav_geti(arch,&r);
        if (r)
            goto RSP1Fin;

        n = sav_geti(arch,&r);
        if (r)
            goto RSP1Fin;

        n = sav_geti(arch,&r);
        if (r)
            goto RSP1Fin;

        if ((r = rd_dat(8,vname)) != 8) {
            printf1("Error: can't read variable name.\n");
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
            n1 = n = sav_geti(arch,&r);
            if (r)
                goto RSP1Fin;

            if (n > 0) {
                n = (n - 1) / 4 + 1;
                n *= 4;
                if (n > SPBUFLEN) {
                    printf1("Error: exceeded max var label length.\n");
                    goto RSP1Fin;
                }
  
                if ((r = rd_dat(n,buf)) != n) {
                    printf1("Error: can't read variable label.\n");
                    goto RSP1Fin; 
                }
                strncpy(vlabel,buf,VLLMax);
                vlabel[VLLMax] = '\0';
                vlabel[imin(VLLMax,n1)] = '\0';
   
                p = vlabel;
                while (*p) {
                    if (*p == '(')
                    *p = '[';
                    else if (*p == ')')
                        *p = ']';
                    p++;
                }
                sprintf(vdef,"%s<4>[0.0](%s)=spss(%d)",vname,vlabel,vtyp);
            }
            else
                sprintf(vdef,"%s<4>[0.0]=spss(%d)",vname,vtyp);
        }
        else
            sprintf(vdef,"%s<4>[0.0]=spss(%d)",vname,vtyp);
    
        for (j = 0; j < iabs(nmiss); ++j) {
            x = sav_getd(arch,&r);
            if (r)
                goto RSP1Fin;
        }
        if (vtyp >= 0) {
            if (save_var(vdef,0)) {       /* save variable definition */

                printf1("Error: can't save variable definitions.\n");
                goto RSP1Fin;
            }
            nvar++;
            if (vtyp > 0) {
                slen = imax(slen,vtyp);
                nvar1++;
            }
        }
    }
    printf1("Number of variables: %d\n",nvar);
    printf1("Number of string variables: %d\n\n",nvar1);

RSP1Cont:
    rtyp = sav_geti(arch,&r);
    if (r)
        goto RSP1Fin;

    if (rtyp == 3) {            /* value labels */
         
        if (dvn == 0 && PMDVARFDef) {
            fprintf(PMDVARFd,"Value labels found in: %s\n\n",PMFdName);
            dvn += 2;
        }
        n = sav_geti(arch,&r);       /* number of labels */
        if (r)
            goto RSP1Fin;

        for (i = 0; i < n; ++i) {

            x = sav_getd(arch,&r);
            if (r)
                goto RSP1Fin;

            if (PMDVARFDef) {
                l = (int)x;
                if ((double)l == x)
                    fprintf(PMDVARFd,"%12d  ",l);
                else
                    fprintf(PMDVARFd,"%12g  ",x);
            }
            if ((r = rd_dat(1,buf)) != 1) {
                printf1("Error in reading value labels.\n");
                goto RSP1Fin; 
            }
            ln = l = (int)buf[0];
            l = l / 8 + 1;
            l *= 8;

            if ((r = rd_dat(l - 1,buf)) != l - 1) {
                printf1("Error in reading value labels.\n");
                goto RSP1Fin; 
            }
            if (PMDVARFDef) {
                *(buf + ln) = '\0';
                fprintf(PMDVARFd,"%s (%d)\n",buf,ln);
                dvn++;
            }
        }

        rtyp = sav_geti(arch,&r);
        if (r)
            goto RSP1Fin;

        if (rtyp != 4) {
            printf1("Error: in reading value label entries.\n");
            goto RSP1Fin;
        }
        n = sav_geti(arch,&r);       /* number of variables */
        if (r)
            goto RSP1Fin;

        for (i = 0; i < n; ++i) {

            l = sav_geti(arch,&r);       /* number of variable */
            if (r)
                goto RSP1Fin;

            l--;
            if (l >= 0 && l < nobs) {
                ii = AcM[l];
                if (PMDVARFDef && ii >= 0 && ii < nvar)  
                    fprintf(PMDVARFd,"%s ",VName[ii]);
            }
        }
        if (PMDVARFDef) {
            fprintf(PMDVARFd,"\n");
            dvn++;
            fprnchar(PMDVARFd,'-',70,1);
            dvn++;
        }
        goto RSP1Cont;
    }
    if (PMDVARFDef) {
        printf1("Value labels: ");
        if (dvn == 0)
            printf1("not present.\n");
        else                
            printf1("%d records written to: %s\n",dvn,PMDVARFName);
    }
    if (rtyp == 6) {        /* document record */

        printf1("Reading document records (will be ignored).\n");

        n = sav_geti(arch,&r);       /* number of records */
        if (r)
            goto RSP1Fin;

        for (i = 0; i < n; ++i) {
            if ((r = rd_dat(80,buf)) != 80)
                goto RSP1Fin; 
        }
        rtyp = sav_geti(arch,&r);
        if (r)
            goto RSP1Fin;
    }
    while (rtyp == 7) {        /* type 7 records */

        n = sav_geti(arch,&r);       /* subtype */
        if (r)
            goto RSP1Fin;

        l = sav_geti(arch,&r);       /* length */
        if (r)
            goto RSP1Fin;

        n = sav_geti(arch,&r);       /* number of elements */
        if (r)
            goto RSP1Fin;

        for (i = 0; i < n; ++i) {
            if ((r = rd_dat(l,buf)) != l)
                goto RSP1Fin; 
        }
        rtyp = sav_geti(arch,&r);
        if (r)
            goto RSP1Fin;

    }
    if (rtyp != 999) {
        printf1("Error: can't find data stream.\n");
        goto RSP1Fin;
    }
    fptr = ftell(PMFd) - 4L;        /* save position to data stream */

    n = sav_geti(arch,&r);
    if (r)
        goto RSP1Fin;

    printf1("Reading data to check variables.\n");

    if (SILENTFlg < 2)
        printfe("Reading: %s\n",PMFdName);

    if (alloc_acx(nobs + 1)) 
        goto RSP1Fin; 
    if (alloc_acy(nobs + 1)) 
        goto RSP1Fin; 

    for (j = 0; j < nobs; ++j) {
        AcX[j] = DBLMAX;
        AcY[j] = DBLMIN;
    }
    eof = nrec = nsys1 = nsys2 = 0;

    if (cflag) {                /* if compressed data */

        j = 0;
        while (eof == 0) {
    
            nr = sav_sval(sval);         
            if (nr == 0) {
                eof = 1;
                break;
            }
            else if (nr < 0) {
                printf1("Error in reading data stream.\n");
                goto RSP1Fin;
            }
            for (i = 0; i < nr; ++i) {
                if (sval[i] == 252) {        /* end of file */
                    eof = 1;
                    break;
                }
                else if (sval[i] == 253) {      /* uncompressed */
                    if (AcN[j] == 0) {
                        x = sav_getd(arch,&r);
                        if (r)
                            goto RSP1Fin;
                    }
                    else {
                        if ((r = rd_dat(8,buf)) != 8)
                            goto RSP1Fin; 
                        x = 0.0;
                    }
                    AcX[j] = dmin(x,AcX[j]);
                    AcY[j] = dmax(x,AcY[j]);
                }
                else if (sval[i] == 254) {
                    x = MBlnkVal;
                    nsys1++;
                }
                else if (sval[i] == 255) {
                    x = PMMSYS;
                    nsys2++;
                }
                else  
                    x = (double)sval[i];

                AcX[j] = dmin(x,AcX[j]);
                AcY[j] = dmax(x,AcY[j]);
                y = floor(x + 0.5);
                if (fabs(y - x) > tol)      
                    AcI[j] = 1;

                if (++j >= nobs) {
                    j = 0;
                    nrec++;
                    prn_message(nrec,0,0);
   
                    if (PMNOCFlg && nrec >= PMNOC)
                        goto RSP1Nxt;
                }
            }
        }
    }
    else {                  /* uncompressed data */
       
        while (1) {

            for (j = 0; j < nobs; ++j) {

                x = sav_getd(arch,&r);
                if (r) {
                    if (j == 0)
                        goto RSP1Nxt;
                    goto RSP1Fin;
                }
                if (x <= smiss) {
                    nsys2++;
                    x = PMMSYS;
                }
                if (AcN[j] != 0)
                    continue;

                AcX[j] = dmin(x,AcX[j]);
                AcY[j] = dmax(x,AcY[j]);
                y = floor(x + 0.5);
                if (fabs(y - x) > tol)      
                    AcI[j] = 1;
            }
            j = 0;
            nrec++;
            prn_message(nrec,0,0);
   
            if (PMNOCFlg && nrec >= PMNOC)
                goto RSP1Nxt;
        }
    }

RSP1Nxt:
    prn_message(nrec,1,0);

    printf1("Read %d records.\n",nrec);

    if (nrec <= 0) {
        printf1("Error.\n");
        goto RSP1Fin;
    }
    if (j) {
        printf1("Error: data stream ends with incomplete record.\n");
        goto RSP1Fin;
    }
    printf1("Number of blank-type missing values: %d\n",nsys1);
    printf1("Number of system-type missing values: %d\n\n",nsys2);

/*** 
for (i = 0; i < nobs; ++i)
    printf1("i=%4d  %f %f aci=%4d acn=%4d   \n",i,AcX[i],AcY[i],AcI[i],AcN[i]);
***/            

    /* set appropriate storage size and write formats */

    ii = 0;
    for (i = 0; i < nobs; ++i) {

        if (AcN[i] < 0)
            continue;

        if (AcN[i] == 0) {
            x = fabs(AcX[i]);
            y = fabs(AcY[i]);
            if (x < y)
                a = y;
            else
                a = x;
    
            if (AcI[i])
                VSLen[ii] = 8;
            else if ((int)a < 128)
                VSLen[ii] = 1;
            else if ((int)a < 32000)
                VSLen[ii] = 2;
            else 
                VSLen[ii] = 5;

            x = AcX[i];
            y = AcY[i];

            if (AcI[i]) {
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
            if (AcI[i]) {
                if (PMFmtF) {
                    pfmt1 = PMFmt1;
                    pfmt2 = PMFmt2;
                }
                else {
                    pfmt1 += 7;
                    pfmt2 = 6;
                }
            }
            makefmt(&pfmt1,&pfmt2,VPFmtS[ii],0,SEPC,0);
            VPFmt1[ii] = (short)pfmt1;
            VPFmt2[ii] = (short)pfmt2;
        }
        ii++;
    }
    prn_var(idxn);    /* print list of new variables */
    newline();     
   
    if (PMF1Def) {
        printf1("Data will be directly written to output file: %s\n",PMF1dName);
        if (PMNOCFlg)
            mxnoc = PMNOC;
        else 
            mxnoc = INTMAX;
    }
    else {
        if (NOCMaxA <= 0)
            NOCMaxA = nrec;
        mxnoc = imin(nrec,NOCMaxA);
       
        printf1("Reading data again to create internal data matrix.\n");
        printf1("Maximum number of cases: %d\n",NOCMaxA);
    }   
    fseek(PMFd,fptr,0);

    rtyp = sav_geti(arch,&r);
    if (r)
        goto RSP1Fin;

    if (rtyp != 999) {
        printf1("Error: can't seek to data stream.\n");
        goto RSP1Fin;
    }
    n = sav_geti(arch,&r);
    if (r)
        goto RSP1Fin;
         
    /* allocate memory for data */
        
    if (PMF1Def == 0) {
        n = MemReq;
        if (alloc_vdat(idxn,1)) {
            printf1("Error: insufficient memory for data matrix (%d cases).\n",NOCMaxA);
            goto RSP1Fin;
        }
        printf1("Allocated %d bytes for data matrix.\n\n",MemReq - n);
    }
    slen1 = (slen - 1) / 8 + 1;
    slen1 *= 8;
    if (alloc_acc(slen1 + 10))        /* buffer for string variables */
        goto RSP1Fin; 

    eof = nrec = 0;
/* ### */
    if (cflag) {                /* if compressed data */

        ii = jlen = jo = j = 0;

        while (eof == 0) {  

            if (nrec >= mxnoc)  
                break;
            
            nr = sav_sval(sval);         

            if (nr == 0) {
                eof = 1;
                break;
            }   
            else if (nr < 0 || nr > 8) {
                printf1("Error in reading data stream.\n");
                goto RSP1Fin;
            }

            for (i = 0; i < nr; ++i) {
                if (sval[i] == 252) {       /* end of file */
                    eof = 1;
                    break;
                }
                else if (sval[i] == 253) {  /* uncompressed */
                    if (AcN[j] == 0) {
                        x = sav_getd(arch,&r);
                        if (r)
                            goto RSP1Fin;
                        if (PMF1Def)
                            fprintf(PMF1d,VPFmtS[ii],x);
                        else
                            put_data(x,ii,nrec);
                        ii++;
                    }
                    else {                
                        if (AcN[j] > 0)  
                            jlen = AcN[j];
                        if ((r = rd_dat(8,AcC + jo)) != 8)
                            goto RSP1Fin; 
                        jo += 8;

                        if (jo >= jlen) {
                            if (PMF1Def) {
                                fwrite(AcC,jlen,1,PMF1d);
                                fprintf(PMF1d," ");
                            }
                            else
                                put_str(AcC,jlen,ii,nrec,0);
                            ii++;
                            jlen = jo = 0;
                        }
                    }
                }
                else if (sval[i] == 254) {
                    if (jo > 0 || AcN[j] > 0) {
                        if (AcN[j] > 0) {
                            jlen = AcN[j];
                            jo = 0;
                        }
                        sprintf(AcC + jo,"        ");
                        jo += 8;
                        if (jo >= jlen) {
                            if (PMF1Def) {
                                fwrite(AcC,jlen,1,PMF1d);
                                fprintf(PMF1d," ");
                            }
                            else
                                put_str(AcC,jlen,ii,nrec,0);
                            ii++;
                            jlen = jo = 0;
                        }
                    }
                    else {
                        x = MBlnkVal;
                        if (PMF1Def)
                            fprintf(PMF1d,VPFmtS[ii],x);
                        else
                            put_data(x,ii,nrec);
                        ii++;
                    }
                }
                else if (sval[i] == 255) {
                    x = PMMSYS;
                    if (PMF1Def)
                        fprintf(PMF1d,VPFmtS[ii],x);
                    else
                        put_data(x,ii,nrec);
                    ii++;
                }
                else { 
                    x = (double)sval[i];
                    if (PMF1Def)
                        fprintf(PMF1d,VPFmtS[ii],x);
                    else
                        put_data(x,ii,nrec);
                    ii++;
                }

                if (++j >= nobs) {
                    ii = j = 0;
                    if (PMF1Def)
                        fprintf(PMF1d,"\n");
                    nrec++;
                    prn_message(nrec,0,0);
                }
            }
      
        }
    }
    else {                  /* uncompressed data */
       
        while (1) {
            ii = 0;
            for (j = 0; j < nobs; ++j) {

                if (AcN[j] == 0) {
                    x = sav_getd(arch,&r);
                    if (r) {
                        if (j == 0)
                            goto RSP1Nxt1;
                        goto RSP1Fin;
                    }
                    if (x <= smiss) {
                        nsys2++;
                        x = PMMSYS;
                    }
                    if (PMF1Def)
                        fprintf(PMF1d,VPFmtS[ii],x);
                    else
                        put_data(x,ii,nrec);
                    ii++;
                }
                else {
                    if (AcN[j] > 0) {
                        jlen = AcN[j];
                        jo = (jlen - 1) / 8 + 1;
                        jo *= 8;
                        if ((r = rd_dat(jo,AcC)) != jo)
                            goto RSP1Fin; 
                        if (PMF1Def) {
                            fwrite(AcC,jlen,1,PMF1d);
                            fprintf(PMF1d," ");
                        }      
                        else
                            put_str(AcC,jlen,ii,nrec,0);
                        ii++;
                    }
                }
            }
            j = 0;
            nrec++;
            if (PMF1Def)
                fprintf(PMF1d,"\n");
            prn_message(nrec,0,0);
   
            if (nrec >= mxnoc)
                break;          
        }
    }

RSP1Nxt1:
    prn_message(nrec,1,0);
    printf1("Read %d records.\n",nrec);

    if (PMF1Def) {
        printf1("%d records written to: %s\n",nrec,PMF1dName);
        if (nrec <= 0)  
            goto RSP1Fin;
    }
    else {
        if (nrec <= 0) {
            printf1("Error: no data matrix created.\n");
            goto RSP1Fin;
        }
        NOCDM = NOC = nrec;
        DMDef = 1;
        printf1("Created a data matrix with %d variables and %d cases.\n",NVAR,NOC);
    }   
    err = 0;

RSP1Fin:
    if (err || PMF1Def) {
        if (idxn < get_nidx())
            clear_avar(idxn);
        NVAR = NOCMaxA = NOCDM = NOC = 0;
    }
    p_clean();
    return(err);
}


/* ------------------------------------------------------------------------ */
/*  rd_dat(n,buf)   read n bytes into buffer.                               */
/*                                                                          */

int rd_dat(int n,char *buf)
{
    int r;
    if ((r = fread(buf,sizeof(char),n,PMFd)) != n) {     
        if (r == 0) {
            *buf = '\0';
            return(0);
        }
        printf1("Error in reading the input file.\n");
        return(-1);
    }
    *(buf + r) = '\0';
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  sav_sval(sval)      get 8 bytes.                                        */

int sav_sval(int *sval)
{
    register int i,n;
    unsigned char buf[9];

    if ((n = rd_dat(8,(char *)buf)) != 8) {
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

int sav_geti(int arch,int *err)
{
    int r,*n;
    char rbuf[6],buf[4];
    register char *p;

    *err = 0;
    if ((r = fread(rbuf,sizeof(char),4,PMFd)) != 4) {           
        *err = -1;
        printf1("Error: can't read next integer.\n");
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
    n = (int *)buf;
    return(*n);
}

/* ------------------------------------------------------------------------ */
/*  sav_getd(arch,err)      return double.                                  */
     
double sav_getd(int arch,int *err)
{
    double *x;
    int r;
    char rbuf[10],buf[8];
    register char *p;

    *err = 0;
    if ((r = fread(rbuf,sizeof(char),8,PMFd)) != 8) {           
        *err = -1;
        printf1("Error (or eof): can't read next double.\n");
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
    x = (double *)buf;
    return(*x);
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

int wr_spss1(void)
{
    register int i,j,k,l;
    int err,nv,sflag,nobs,n,ii,jj,olen,lbuf[9];
    register char *p;
    unsigned char buf[9];
    double x,dbuf[9];
    char *pbuf[9];




    err = -1;
    if (check_cmd(0))
        return(-1);

    if (parm(CmdBuf + 6,13,1)) {   /* get parameters */
        p_clean();
        return(-1);   
    }       
    printf1("Writing SPSS sav file: %s\n",PMFdName);

    if (PMKeep && PMDrop) {
        p_err(-16,1);
        goto WSP1Fin;
    }

    /* get list of variables in AcI[] */

    if (alloc_aci(imax(NVAR,PMNV)))
        goto WSP1Fin;

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
        goto WSP1Fin;
    }
    sflag = 0;
    if (PM1NV > 0) {        /* sort */
        err = vsort(PM1NV,PM1VIdx,1,0,1);
        if (err)
            goto WSP1Fin;
        sflag = 1;
    }

    nobs = 0;
    for (i = 0; i < nv; ++i) {
        j = AcI[i];
        if (VTyp[j] != 1) {
            nobs++;
        }
        else {
            l = -VSLen[j];
            n = (l - 1) / 8 + 1;
            nobs += n;
        }
    }
    if (alloc_ack(nobs + 1))    /* tda var index */
        goto WSP1Fin;
    if (alloc_acj(nobs + 1))    /* string length */
        goto WSP1Fin;
    if (alloc_acs(nobs + 1))    /* var label */
        goto WSP1Fin;
    if (alloc_acn(nobs + 1))    /* type */
        goto WSP1Fin;
    if (alloc_acm(nobs + 1))    /* format */
        goto WSP1Fin;

    nobs = 0;
    for (i = 0; i < nv; ++i) {
        j = AcI[i];
        if (VTyp[j] != 1) {
            AcK[nobs] = j;
            AcN[nobs] = 0;
            if (VPFmt1[i] > 0 && VPFmt2[i] < VPFmt1[i] - 1)
                AcM[nobs] = (int)VPFmt2[i] + (int)VPFmt1[i] * 256 + 5 * 256 * 256;
            else
                AcM[nobs] = 4 + 10 * 256 + 5 * 256 * 256;
            nobs++;
        }
        else {
            l = -VSLen[j];
            n = (l - 1) / 8 + 1;
            AcK[nobs] = j;
            AcN[nobs] = l;
            AcJ[nobs] = l;
            AcM[nobs] = l * 256 + 256 * 256;
            nobs++;
            for (k = 1; k < n; ++k) {
                AcK[nobs] = j;
                AcN[nobs] = -1;
                AcJ[nobs] = l;
                AcM[nobs] = l * 256 + 256 * 256;
                nobs++;
            }
        }
    }

    /* print header */

    fwrite("$FL2@(#)SPSS DATA FILE. SAV FILE CREATED BY TDA.",48,1,PMFd);
    fprnchar(PMFd,' ',16,0);

    sav_puti(2);            /* file layout code */
    sav_puti(nobs);         /* number of OBS */
    sav_puti(1);            /* compression code, always compressed */
    sav_puti(0);            /* case weight variable index */
    sav_puti(NOC);          /* number of cases */
    sav_putd(100.0);        /* compression bias */

    fprnchar(PMFd,' ',84,0);        /* no date, etc. */ 

    for (i = 0; i < nobs; ++i) {    /* write dictionary */

        j = AcK[i];
        sav_puti(2);            /* record type */
        sav_puti(AcN[i]);       /* var type */

        ii = 0;
        if (VLabel[j] != NULL) {
            ii = strlen(VLabel[j]);
            n = (ii - 1) / 4;
            n++;
            n *= 4;
            sav_puti(1);            /* label */
        }
        else
            sav_puti(0);            /* no label */
        sav_puti(0);            /* no missing value codes */
        sav_puti(AcM[i]);       /* print format */
        sav_puti(AcM[i]);       /* write format */

        p = VName[j];
        for (k = 0; k < 8; ++k) {
            if (*p)
                fprintf(PMFd,"%c",*p++);
            else
                fprintf(PMFd," ");
        }
        if (ii > 0) {
            sav_puti(ii);   
            fprintf(PMFd,"%s",VLabel[j]);
            for (k = ii; k < n; ++k)
                fprintf(PMFd," ");
        }
    }

    /* write data, always compressed */

    sav_puti(999);          /* record type */
    sav_puti(0);            /* filler */
         
    k = 0;
    for (i = 0; i < NOC; ++i) {
        ii = i;
        if (sflag)
            ii = VSORTPtr[i];

        for (j = 0; j < nobs; ++j) {

            jj = AcK[j];
            if (AcN[j] == 0) {
                x = get_data(jj,ii);
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
                if (AcN[j] > 0) {
                    pbuf[k] = VDPtr[jj] + ii * AcJ[j];                
                    olen = 8;
                }   
                else {
                    pbuf[k] = VDPtr[jj] + ii * AcJ[j] + olen;           
                    olen += 8;
                }   
                
                if (olen <= AcJ[j])  
                    lbuf[k] = 8;
                else
                    lbuf[k] = 8 - (olen - AcJ[j]);
            }
            if (++k >= 8) {
                fwrite(buf,8,1,PMFd);
                for (k = 0; k < 8; ++k) {
                    if (buf[k] == 0xfd) {
                        if (pbuf[k] == NULL) {
                            sav_putd(dbuf[k]);
                        }
                        else {
                            fwrite(pbuf[k],lbuf[k],1,PMFd);
                            for (l = lbuf[k]; l < 8; ++l)
                                fprintf(PMFd," ");
                        }
                    }
                }
                k = 0;
            }
        }
        prn_message(i + 1,0,1);
    }
    if (k > 0) {
        while (k < 8)
            buf[k++] = 0xfc;          /* end of file */

        fwrite(buf,8,1,PMFd);
        for (k = 0; k < 8; ++k) {
            if (buf[k] == 0xfd) {
                if (pbuf[k] == NULL)
                    sav_putd(dbuf[k]);
                else {
                    fwrite(pbuf[k],lbuf[k],1,PMFd);
                    for (l = lbuf[k]; l < 8; ++l)
                        fprintf(PMFd," ");
                }
            }
        }
    }
    prn_message(NOC,1,1);

    printf1("%d records with %d variables written to: %s\n",NOC,nv,PMFdName);
    err = 0;
          
WSP1Fin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sav_puti(n)     write integer.                                          */

void sav_puti(int n)
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
    fwrite(buf,4,1,PMFd);
}

/* ------------------------------------------------------------------------ */
/*  sav_putd(x)     write double.                                           */
     
void sav_putd(double x)
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
    fwrite(buf,8,1,PMFd);
}

/* ------------------------------------------------------------------------ */
/*  get_dvarp(np)   get info about dvar partition.                          */
     
int get_dvarp(int *np,char *vname)
{
    register int j;   

    if (PMDVARP > 0) {
        *np = 1;
        return(PMDVARP);
    }
    else if (PMDVARPN > 0) {
        *np = 2;
        return(PMDVARPN);
    }
    else if (PMDVARVN > 0) {
        for (j = 0; j < PMDVARVN; ++j) {
            if (!strcmp(PMDVARVName[j],vname)) {
                *np = PMDVARVNP[j];
                return(PMDVARVNL[j]);
            }
        }
    }
    *np = 0;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  make_arcd()    update archive description file.                         */
     
void make_arcd(int fn,char *fname,int len,int nrec,int vn)
{
    int addc,n;

    addc = 1;
#if S_DOS               /* add to record length */
    addc++;
#endif
    n = 0;
    if (PMDVARFDef) {
        fclose(PMDVARFd);
        if ((PMDVARFd = fopen(PMDVARFName,OPEN_RD))) {   
            if (alloc_acc(1001) == 0) {  
                while (fgets(AcC,1000,PMDVARFd))  
                    n++;
                alloc_acc(0);  
            }
        }
    }
    if (PMARCFZOO)
        fprintf(PMARCFd,"%s\n",PMARCFZOOF);
    fprintf(PMARCFd,"%3d %s ",fn,fname);
    fprnchar(PMARCFd,' ',20 - strlen(fname),0);
    fprintf(PMARCFd,"1 %6d %8d %6d\n",len + addc,nrec,vn);
    if (PMARCFVDF) {
        fprintf(PMARCFd,"999 %s ",PMARCFVDFF);
        fprnchar(PMARCFd,' ',20 - strlen(PMARCFVDFF),0);
        fprintf(PMARCFd,"2 %6d %8d %6d\n",0,n,0);
    }
    printf1("Added info to archive description file: %s\n",PMARCFName);
}


