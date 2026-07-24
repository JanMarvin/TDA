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

/*  functions in t_xls.c */

int rxls(void); 
int xls_read(void);
int xls_short(short *n);
int xls_int(int *n);
int xls_dbl(double *x);
int xls_rk(double *x);
int xls_unknown(void);
int xls_boundsheet(void);
int xls_sst(void);
int xls_labelsst(void);
int xls_blank(void);
int xls_mulblank(void);
int xls_mulrk(void);
int xls_formula(void);
int xls_number(void);
void xls_rowcol(int row,int col);
void xls_putx(int row,int col,double x);
int xls_check(int j);
int rucinet(void); 


/*--------------------------------------------------------------------------*/
/*  Global variables                                                        */

#define XLSHLEN 64          /* used to skip over OLE header                 */
#define XLSNFMax 50         /* max number of files                          */
int XLSV = 8;               /* Version: 5 or 8                              */
int XLSNF = 0;              /* actual number of files                       */
long XLSFPtr[XLSNFMax];     /* file pointer                                 */

int XLSRF = 0;              /* set to 1 for second read                     */
int ProtNRec = 0;           /* number of records writen to protocol file    */  
int XLSRowMin = 0;
int XLSRowMax = 0;
int XLSColMin = 0;
int XLSColMax = 0;
int XLSRow = 0;             /* number of rows in AcX                        */
int XLSCol = 0;             /* number of columns in AcX                     */
int XLSLMin = 0;            /* min number of label                          */
int XLSLMax = 0;            /* max number of label                          */

int XLSLN = 0;              /* number of labels                             */  
int XLSLNA = 0;             /* number of saved labels                       */
char **XLSLabel;            /* labels                                       */
int XLSLNN = 0;             /* allocated                                    */

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

int rxls(void)  
{
    register int i,j,k,jj;
    char buf[XLSHLEN + 1];
    short sn,len,ver,typ;
    int err,ii,l,n,nc,ne,nrow,ncol,eof,nrec;
    long fptr;

    err = -1;
    if (check_cmd(0))
        return(-1);

    if (parm(CmdBuf + 4,10,1)) {   /* get parameters */
        p_clean();
        return(-1);   
    }   
    printf1("Reading xls file: %s\n\n",PMFdName);
    if (PMFmtF == 0)
        pmfmt(10,4);

    ProtNRec = 0;
    XLSLNN = XLSLNA = XLSLN = XLSNF = XLSRF = 0;
    XLSLMin = XLSRowMin = XLSColMin = INTMAX;
    XLSLMax = XLSRowMax = XLSColMax = -1;
           
    /* skip the OLE header until one finds begin of BIFF8 */

    ii = 0;
    fptr = 0L;
    while (fread(buf,sizeof(char),XLSHLEN,PMFd) == XLSHLEN) {                   
        if (*buf == 0x09 && *(buf + 1) == 0x08) {
            ii = 1;
            fseek(PMFd,fptr,0);
            break;
        }
        fptr = ftell(PMFd);
    }
    if (ii == 0) {
        printf1("Error: cannot find BIFF8 begin of file (0x0809).\n");
        goto RXLSFin;
    }
    eof = 0;
    XLSV = 0;

    while (1) {

        if (xls_short(&sn)) {           /* should be BOF */
            if (eof)
                break;
            goto RXLSFin;
        }
        if (sn != 0x0809) {
            if (eof == 0)
                printf1("Error: cannot find BIFF8 begin of file (0x0809). Found: %04x\n",sn);
            break;             
        }
        if (xls_short(&len))     /* length */
            goto RXLSFin;

        if (xls_short(&ver))     /* version */
            goto RXLSFin;

        if (XLSV == 0) {
            printf1("Version: %04x ",ver);
            if (ver == 0x0600 && len == 16) {
                XLSV = 8;
                printf1("BIFF8\n");
            }
            else if (len == 8) {
                XLSV = 5;
                newline();
            }
            else {
                printf1("[cannot use this version]\n");
                goto RXLSFin;
            }
        }
        if (xls_short(&typ))     /* type */
            goto RXLSFin;

        len -= 4;
        if (len < 4 || len > 12) {
            printf1("Error: cannot interpret BOF.\n");
            goto RXLSFin;
        }
        if (fread(buf,sizeof(char),len,PMFd) != len) {    /* skip rest of BOF */
            p_err(-7,1);
            goto RXLSFin;
        }
        printf1("Begin of new file: ");
        fptr = ftell(PMFd);
        eof = 0;
        switch (typ) {
            case 0x0005:    if (XLSNF != 0) {
                                printf1("Error: workbook globals should be the first file.\n");
                                goto RXLSFin;
                            }
                            XLSFPtr[XLSNF] = fptr;
                            XLSNF++;

                            printf1("workbook globals.\n");
                            if (xls_read())
                                goto RXLSFin;
                            eof = 1;
                            break;

            case 0x0006:    printf1("visual basic module.\n");
                            break;

            case 0x0010:    if (XLSNF == 0) {
                                printf1("Error: worksheet should follow a globals section.\n");
                                goto RXLSFin;
                            }
                            if (XLSNF >= XLSNFMax) {
                                printf1("Error: exceeded max number of files (%d).\n",XLSNFMax);
                                goto RXLSFin;
                            }
                            XLSFPtr[XLSNF] = fptr;
                            XLSNF++;

                            printf1("worksheet.\n");
                            if (xls_read())
                                goto RXLSFin;
                            eof = 1;
                            break;

            case 0x0020:    printf1("chart.\n");
                            break;

            case 0x0040:    printf1("BIFF4 macro sheet.\n");
                            break;

            case 0x0100:    printf1("BIFF4 workbooks globals.\n");
                            break;
            default:        printf1("unknown.\n");
        }                   
        if (eof == 0) {
            printf1("Not supported.\n");
            goto RXLSFin;
        }
    }
    printf1("Number of worksheets: %d\n\n",XLSNF - 1);
    if (XLSNF < 2)
        goto RXLSFin;

    /* save strings from SST */

    if (!(XLSLabel = (char **)calloc(XLSLN,sizeof(char *)))) {
        p_err(-2,1);
        XLSLN = 0;             
        goto RXLSFin;
    }
    XLSLNN = XLSLN;
    memrq(XLSLNN,sizeof(char *));

    XLSRF = 1;
    fseek(PMFd,XLSFPtr[0],0);
    if (xls_read())
        goto RXLSFin;

    /* get data for all worksheets */

    XLSRow = XLSRowMax + 1;
    XLSCol = XLSColMax + 1;
    nrow = XLSRowMax - XLSRowMin + 1;
    ncol = XLSColMax - XLSColMin + 1;

    printf1("Maximal size of table: %d (%d - %d) rows, %d (%d - %d) columns.\n",
        nrow,XLSRowMin,XLSRowMax,ncol,XLSColMin,XLSColMax);

    if (XLSRowMin > XLSRowMax || XLSColMin > XLSColMax)
        goto RXLSFin;

    if (alloc_acx(XLSRow * XLSCol + 1))
        goto RXLSFin;
    if (alloc_acn(XLSRow * XLSCol + 1))     /* >0 = number of label */
        goto RXLSFin;                       /* -1 not used */
                                            /* -2 numerical */
                                            /* -3 blank */
    if (alloc_aci(XLSRow + 1))
        goto RXLSFin;
    if (alloc_acj(XLSCol + 1))
        goto RXLSFin;
    if (alloc_ack(XLSCol + 1))
        goto RXLSFin;

    XLSRF = 2;
    nrec = 0;
    for (ii = 1; ii < XLSNF; ++ii) {

        printf1("\nReading again: worksheet %d\n",ii);

        for (i = 0; i <= XLSRow * XLSCol; ++i) {
            AcN[i] = -1;
            AcX[i] = 0.0;
        }
        for (i = 0; i <= XLSRow; ++i)  
            AcI[i] =  0;
        for (j = 0; j <= XLSCol; ++j)  
            AcK[j] = AcJ[j] =  0;

        fseek(PMFd,XLSFPtr[ii],0);
        if (xls_read())             /* read worksheet ii */
            goto RXLSFin;

        /***
        for (i = XLSRowMin; i <= XLSRowMax; ++i) {
            printf("%5d : ",AcI[i]);
            for (j = XLSColMin; j <= XLSColMax; ++j)  
                printf("%3d ",AcN[i * XLSCol + j]);
            newline();
        }
        ***/

        ne = 0;
        for (i = XLSRowMin; i <= XLSRowMax; ++i) {
            if (AcI[i] > 0)
                ne++;
        }
        printf1("%d rows contain at least one numerical entry and will be written\n",ne);
        printf1("as data records to the output file. The following columns will be used.\n\n");
         
        printf1("Column   New  Type         Row  Label\n");
        prnchar('-',50,1);                       

        nc = 0;
        for (j = XLSColMin; j <= XLSColMax; ++j) {

            n = xls_check(j);
            if (n == 1) {
                printf1("%6d %5d  numerical",j,nc);
                AcJ[nc++] = j;    
            }
            else if (n == 2) {
                if (PMNS != 1) {
                    printf1("%6d %5d  string   ",j,nc);
                    AcJ[nc++] = j;
                    AcK[j] = 1;                   
                }
                else
                    printf1("%6d        string   ",j);
            }
            else        
                printf1("%6d        not used ",j);

            n = 0;
            for (i = XLSRowMin; i <= XLSRowMax; ++i) {
                jj = i * XLSCol + j;
                k = AcN[jj];        
                if (k >= 0 && k < XLSLNA) {
                    if (AcI[i] <= 0) {
                        if (n)
                            prnchar(' ',23,0);
                        printf1("  %5d  %s\n",i,XLSLabel[k]);
                        n++;
                    }
                }   
                if (AcI[i] > 0 && AcK[j]) {           /* string */ 
                    if (k == -2) {
                        sprintf(buf,"%-g",AcX[jj]);
                        AcK[j] = imax(AcK[j],strlen(buf));
                    }
                    else if (k >= 0 && k < XLSLNA)
                        AcK[j] = imax(AcK[j],strlen(XLSLabel[k]));
                }
            }
            if (n == 0)
                newline();
        }
        if (ne == 0)
            continue;

        if (PMF1Def) {  

            for (i = XLSRowMin; i <= XLSRowMax; ++i) {
                if (AcI[i] > 0) {
                    if (XLSNF > 2)
                        fprintf(PMF1d,"%4d ",ii);

                    for (jj = 0; jj < nc; ++jj) {
                        j = AcJ[jj];
                        k = i * XLSCol + j;
                        if (AcK[j] == 0) {
                            if (AcN[k] == -2)
                                fprintf(PMF1d,PMFmtS,AcX[k]);
                            else
                                fprintf(PMF1d,PMFmtS,-1.0);
                        }
                        else {                        /* string */
                            if (AcN[k] >= 0 && AcN[k] < XLSLNA) {
                                fprintf(PMF1d,"%s ",XLSLabel[AcN[k]]);
                                l = strlen(XLSLabel[AcN[k]]);
                                fprnchar(PMF1d,' ',imax(0,AcK[j] - l),0);
                            }
                            else if (AcN[k] == -2) {
                                sprintf(buf,"%-g",AcX[k]);
                                fprintf(PMF1d,"%s ",buf);
                                l = strlen(buf);
                                fprnchar(PMF1d,' ',imax(0,AcK[j] - l),0);
                            }
                            else
                                fprnchar(PMF1d,' ',AcK[j] + 1,0);
                        }
                    }
                    fprintf(PMF1d,"\n");
                    nrec++;
                }
                else if (PMNC != 1) {
                    n = 1;
                    for (j = XLSColMin; j <= XLSColMax; ++j) {
                        jj = i * XLSCol + j;
                        k = AcN[jj];
                        if (k >= 0 && k < XLSLNA) {
                            if (n) {
                                fprintf(PMF1d,"# ");
                                n = 0;
                            }
                            fprintf(PMF1d,"[%d] %s ",j,XLSLabel[k]);
                        }
                    }
                    if (n == 0) {
                        fprintf(PMF1d,"\n");
                        nrec++;
                    }
                }
            }
        }
    }
    if (PMF1Def)    
        printf1("\n%d records written to: %s\n",nrec,PMF1dName);
    if (PMProtFd)
        printf1("%d records written to: %s\n",ProtNRec,PMProtFName);
    err = 0;

RXLSFin:
    if (XLSLNN > 0) {
        for (i = 0; i < XLSLNA; ++i) {
            n = strlen(XLSLabel[i]) + 1;
            free(XLSLabel[i]);
            memrq(-n,sizeof(char));
        }
        free((char *)XLSLabel);
        memrq(-XLSLNN,sizeof(char *));
    }
    p_clean();
    return(err);
}


/* ------------------------------------------------------------------------ */
/*  xls_read            Read one file.                                      */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_read(void)
{
    short sn;
    long fptr;

    while (1) {

        if (xls_short(&sn))
            return(-1);     

        switch (sn) {
            case 0x0006:    if (xls_formula())          /* formula, value */
                                return(-1);    
                            break;

            case 0x0085:    if (xls_boundsheet())
                                return(-1);    
                            break;

            case 0x00fc:    if (xls_sst())
                                return(-1);     
                            break;

            case 0x000a:    if (xls_short(&sn))         /* EOF */
                                return(-1);     
                            return(0);

            case 0x00bd:    if (xls_mulrk())            /* multiple RK */
                                return(-1);     
                            break;       

            case 0x00be:    if (xls_mulblank())         /* multiple blanks */
                                return(-1);     
                            break;       

            case 0x00fd:    if (xls_labelsst())         /* label SST */
                                return(-1);     
                            break;       

            case 0x0201:    if (xls_blank())            /* blank cell */
                                return(-1);     
                            break;       

            case 0x0203:    if (xls_number())           /* number */
                                return(-1);     
                            break;       

            case 0x027e:    if (xls_rkvalue())          /* rk-value */
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

                            if (xls_unknown())
                                return(-1);    
                            break;
            default:        if (PMNI == 1) {
                                if (xls_unknown())
                                    return(-1);    
                                break;
                            }
                            printf1("Error: unknown record type: %04x\n",sn);
                            fptr = ftell(PMFd);
                            printf1("File position: %08x\n",fptr - 2);
                            return(-1);   
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_short   Read next 2 bytes and return as short integer.              */
/*              Return 0 if OK, -1 if error.                                */  

int xls_short(short *sn)
{
    char buf[3];

    if (fread(buf,sizeof(char),2,PMFd) != 2) {
        p_err(-7,1);
        return(-1);       
    }
    *sn = *(short *)buf;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_int     Read next 4 bytes and return as integer.                    */
/*              Return 0 if OK, -1 if error.                                */  

int xls_int(int *n)
{
    char buf[5];

    if (fread(buf,sizeof(char),4,PMFd) != 4) {
        p_err(-7,1);
        return(-1);       
    }
    *n = *(int *)buf;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_dbl     Read next 8 bytes and return as double.                     */
/*              Return 0 if OK, -1 if error.                                */  

int xls_dbl(double *x)
{
    char buf[9];

    if (fread(buf,sizeof(char),8,PMFd) != 8) {
        p_err(-7,1);
        return(-1);       
    }
    *x = *(double *)buf;
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  xls_rk      Read next 4 bytes and return as double (RK value).          */
/*              Return 0 if OK, -1 if error.                                */  

int xls_rk(double *x)
{
    int m,i;
    unsigned char buf[5],fbuf[11],*p;
    unsigned int n;
    float f;
                  
    if (fread(buf,sizeof(char),4,PMFd) != 4) {
        p_err(-7,1);
        return(-1);       
    }
    m = 0x01 & buf[0];
    i = 0x02 & buf[0];
    /**   
    printf("m=%d i=%d\n",m,i);        
    printf("buf: %02x %02x %02x %02x\n",buf[0],buf[1],buf[2],buf[3]);          
    **/   

    if (i) {                           /* integer */
        n = *(unsigned int *)buf;                                 
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
        *x = *(double *)fbuf;
    }
    if (m)   
        *x /= 100.0;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_unknown                                                             */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_unknown(void)
{
    char buf[10001];
    short sn;
    long fptr;

    if (xls_short(&sn))
        return(-1);       
           
    if (sn == 0)
        return(0);

    if (sn < 1 || sn > 10000) {
        printf1("Error in length: %d\n",sn);
        fptr = ftell(PMFd);
        printf1("File position: %08x\n",fptr);
        return(-1);
    }
    if (fread(buf,sizeof(char),(int)sn,PMFd) != (int)sn) {
        p_err(-7,1);
        return(-1);    
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_boundsheet                                                          */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_boundsheet(void)
{
    unsigned char buf[257];
    short len,sn,rlen;

    if (xls_short(&len))
        return(-1);       
    if (xls_short(&sn))
        return(-1);       
    if (xls_short(&sn))
        return(-1);       
    if (xls_short(&sn))     /* options */
        return(-1);       

    if (XLSV == 8) {        /* get length of following string */
        if (xls_short(&sn))                    
            return(-1);       
        rlen = sn + 8l;
    }    
    else {
        if (fread(buf,sizeof(char),1,PMFd) != 1) {      
            p_err(-7,1);                                
            return(-1);       
        }
        sn = *buf;
        rlen = sn + 7; 
    }
    if (rlen != len) {
        printf1("Boundsheet: error (%d, %d).\n",sn,len);
        return(-1);
    }
    if (sn < 1 || sn > 256) {
        printf1("Boundsheet: error in length: %d\n",sn);
        return(-1);
    }
    if (fread(buf,sizeof(char),(int)sn,PMFd) != (int)sn) {
        p_err(-7,1);
        return(-1);    
    }
    buf[(int)sn] = '\0';
    if (XLSRF == 0)
        printf1("Boundsheet: %s\n",buf);
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  xls_sst                                                                 */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_sst(void)
{
    register char *p;
    char buf[10001];
    short len,sn,rtn;      
    int i,l,n,slen,rlen,rtlen;
    unsigned char c;

    if (xls_short(&len))
        return(-1);       

    if (xls_int(&n))
        return(-1);       

    if (PMProtFDef && XLSRF == 0) {
        fprintf(PMProtFd,"SST: total entries: %d\n",n);
        ProtNRec++;
    }
    if (xls_int(&n))
        return(-1);       

    if (PMProtFDef && XLSRF == 0) {
        fprintf(PMProtFd,"SST: different entries: %d\n\n",n);
        ProtNRec += 2;
    }
    rlen = 8;

    for (i = 0; i < n; ++i) {
        if (xls_short(&sn))
            return(-1);       
        slen = (int)sn;

        if (rlen >= len && sn == 0x003c) {             /* continue */
            if (xls_short(&sn))
                return(-1);       
            len += (int)sn;
            /** 
            printf("continue sn=%04x %d new len=%d rlen=%d i=%d n=%d\n",sn,sn,len,rlen,i,n);     
            **/ 
            if (xls_short(&sn))
                return(-1);       
            slen = (int)sn;
        }

        if (fread(buf,sizeof(char),1,PMFd) != 1) {
            p_err(-7,1);
            return(-1);    
        }
        c = buf[0];
        if (c & 0x01)       /* 16 bit characters */
            slen *= 2;

        if (c & 0x04) {
            printf1("SST error: cannot interpret far-east information (%04x).\n",c);
/*          return(-1);     */
        }
        rtlen = 0;
        if (c & 0x08) {
            if (xls_short(&rtn))
                return(-1);       
            rtlen = 4 * rtn;
            rlen += 2;
        }
        if (fread(buf,sizeof(char),slen,PMFd) != slen) {
            p_err(-7,1);
            return(-1);    
        }
        buf[slen] = '\0';

        if (PMProtFDef && XLSRF == 0) {
            fprintf(PMProtFd,"%s\n",buf);
            ProtNRec++;
        }
        if (XLSRF == 0)
            XLSLN++;
    
        else if (XLSRF == 1) {
            if (XLSLNA >= XLSLN) {
                printf1("SST error: exceeded max number of strings (%d).\n",XLSLN);
                exit(0);
            }
            l = strlen(buf);
            if (!(XLSLabel[XLSLNA] = (char *)calloc(l + 1,sizeof(char)))) {
                p_err(-2,1);
                return(-1);            
            }
            memrq(l + 1,sizeof(char));

            p = buf;
            while (*p) {
                if (*p == '\n')
                    *p = ' ';
                p++;
            }
            strcpy(XLSLabel[XLSLNA],buf);
            XLSLNA++;
        }
        rlen += slen + 3;
        if (rtlen > 0) {
            if (fread(buf,sizeof(char),rtlen,PMFd) != rtlen) {
                p_err(-7,1);
                return(-1);    
            }
            rlen += rtlen;
        }
    }
    if (rlen != len) {
        printf1("SST error: %d, %d\n",rlen,len);
        return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_labelsst                                                            */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_labelsst(void)
{
    short sn,row,col,xf;
    int sst;

    if (xls_short(&sn))
        return(-1);       

    if (xls_short(&row))
        return(-1);       
    if (xls_short(&col))
        return(-1);       
    if (xls_short(&xf))
        return(-1);       
    if (xls_int(&sst))
        return(-1);       
           
    XLSLMin = imin(XLSLMin,sst);
    XLSLMax = imax(XLSLMax,sst);
       
    if (XLSRF == 2)  
        AcN[row * XLSCol + col] = sst;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_blank                                                               */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_blank(void)
{
    short sn,row,col;

    if (xls_short(&sn))
        return(-1);       
    if (xls_short(&row))
        return(-1);       
    if (xls_short(&col))
        return(-1);       
    if (xls_short(&sn))
        return(-1);       

    if (XLSRF == 0)
        xls_rowcol(row,col);
    else if (XLSRF == 2)  
        AcN[row * XLSCol + col] = -3; 
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  xls_mulblank                                                            */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_mulblank(void)
{
    short sn,row,col,col1,xf;
    int i,n,len;

    if (xls_short(&sn))
        return(-1);       
    if (xls_short(&row))
        return(-1);       
    if (xls_short(&col))
        return(-1);       

    len = sn - 6;
    n = len / 2;
    if (n * 2 != len) {
        printf1("MULBLANK error.\n");
        return(-1);
    }
    for (i = 0; i < n; ++i) {
        if (xls_short(&xf))
            return(-1);       

        if (XLSRF == 0)
            xls_rowcol(row,col + i);
        else if (XLSRF == 2)  
            AcN[row * XLSCol + col] = -3; 
    }        
    if (xls_short(&col1))
        return(-1);       

    if (col + n - 1 != col1) {
        printf1("MULBLANK1 error: %d, %d, %d.\n",col,col1,n);
        return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_mulrk                                                               */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_mulrk(void)
{
    short sn,xf,row,col,col1;
    int i,n,len;
    double x;

    if (xls_short(&sn))
        return(-1);       
    if (xls_short(&row))
        return(-1);       
    if (xls_short(&col))
        return(-1);       

    len = sn - 6;
    n = len / 6;
    if (n * 6 != len) {
        printf1("MULRK error.\n");
        return(-1);
    }
    for (i = 0; i < n; ++i) {
        if (xls_short(&xf))
            return(-1);       
        if (xls_rk(&x))
            return(-1);       

        if (XLSRF == 0) 
            xls_rowcol(row,col + i);
        else if (XLSRF == 2)
            xls_putx(row,col + i,x); 
    }        
    if (xls_short(&col1))
        return(-1);       

    if (col + n - 1 != col1) {
        printf1("MULRK1 error: %d, %d, %d.\n",col,col1,n);
        return(-1);
    }
    /**
    printf1("MulRK:     %4d %4d %4d\n",row,col,col1);          
    **/
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_formula                                                             */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_formula(void)
{
    short sn,row,col,xf,opt;
    int len;
    char buf[100001];
    double x;

    if (xls_short(&sn))
        return(-1);       
    len = sn;

    if (xls_short(&row))
        return(-1);       
    if (xls_short(&col))
        return(-1);       
    if (xls_short(&xf))
        return(-1);       

    if (xls_dbl(&x))
        return(-1);       

    if (xls_short(&opt))
        return(-1);       

    len -= 16;
    if (len > 0) {
        if (fread(buf,sizeof(char),len,PMFd) != len) {
            p_err(-7,1);
            return(-1);    
        }
    }
    if (XLSRF == 0)
        xls_rowcol(row,col);
    else if (XLSRF == 2)
        xls_putx(row,col,x);
    /**
    printf1("Formula:   %4d %4d xf=%d opt=%04x value: %lg\n",row,col,xf,opt,x);
    **/
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_number                                                              */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_number(void)
{
    short sn,row,col,xf;
    double x;

    if (xls_short(&sn))
        return(-1);       

    if (xls_short(&row))
        return(-1);       
    if (xls_short(&col))
        return(-1);       
    if (xls_short(&xf))
        return(-1);       
    if (xls_dbl(&x))
        return(-1);       

    if (XLSRF == 0)
        xls_rowcol(row,col);
    else if (XLSRF == 2)
        xls_putx(row,col,x);

    /**
    printf1("Number:    %4d %4d xf=%d value: %lg\n",row,col,xf,x);     
    **/
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_rkvalue                                                             */
/*  Return 0 if OK, -1 if error.                                            */  

int xls_rkvalue(void)
{
    short sn,row,col,xf;
    double x;

    if (xls_short(&sn))
        return(-1);       

    if (xls_short(&row))
        return(-1);       
    if (xls_short(&col))
        return(-1);       
    if (xls_short(&xf))
        return(-1);       
    if (xls_rk(&x))
        return(-1);       

    if (XLSRF == 0)
        xls_rowcol(row,col);
    else if (XLSRF == 2)
        xls_putx(row,col,x);

    /**
    printf1("RK-Value:  %4d %4d xf=%d value: %lg\n",row,col,xf,x);  
    **/
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xls_rowcol                                                              */
/*  Return 0 if OK, -1 if error.                                            */  

void xls_rowcol(int row,int col)
{
    if (XLSRF)
        return;

    XLSRowMin = imin(XLSRowMin,row);
    XLSRowMax = imax(XLSRowMax,row);
    XLSColMin = imin(XLSColMin,col);
    XLSColMax = imax(XLSColMax,col);
}

/* ------------------------------------------------------------------------ */
/*  xls_putx(x)     Put x into AcX.                                         */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */  

void xls_putx(int row,int col,double x)
{
    int i;

    if (XLSRF == 2) {
        i = row * XLSCol + col;
        AcX[i] = x;
        AcN[i] = -2;
        AcI[row] += 1;
    }
}

/* ------------------------------------------------------------------------ */
/*  xls_check(j)    Return 1 if column j contains only numerical entries,   */
/*                         2 if column j contains at least one string,      */
/*                         0 otherwise                                      */
/*                  Take into account only rows with at least one numerical */
/*                  entry.                                                  */

int xls_check(int j)
{
    int i,ij,nn,ns;

    nn = ns = 0;
    for (i = 0; i <= XLSRowMax; ++i) {
        if (AcI[i] > 0) {
            ij = i * XLSCol + j;
            if (AcN[ij] == -2)
                nn++;
            else if (AcN[ij] >= 0)
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

int rucinet(void)  
{
    char buf[50];
    int err,nc,nrec,i,j;
    float x;    

    err = -1;
    if (check_cmd(0))
        return(-1);

    if (parm(CmdBuf + 7,10,1)) {   /* get parameters */
        p_clean();
        return(-1);   
    }   
    printf1("Reading Ucinet file: %s\n\n",PMFdName);
    if (PMFmtF == 0)
        pmfmt(10,4);

    /**********************
    fseek(PMFd,0x19da,0);
    i = 1;
    while (fread(buf,sizeof(char),1,PMFd) == 1) {                        
        if (buf[0])
            printf("%c",buf[0]);
        else {
            if (++i >= 43) {
                newline();
                i = 0;
            } 
            else
                printf(" ");
        }
    }
    *******************************/

    if (PMF1Def == 0) {
        printf1("Error: need output file (df parameter).\n");
        goto RUCIFin;
    }            

    /*** 
    fread(buf,sizeof(char),27,PMFd);                               
    while (fread(buf,sizeof(char),20,PMFd) == 20)                          
        printf("%s,\n",buf);
    goto RUCIFin;
    ***/

    nc = nrec = 0;
    i = j = 1;

    while (fread(buf,sizeof(char),4,PMFd) == 4) {                        
        x = *(float *)buf;
        if (PMPRNO == 1) {
            if (x > 0.0) {
                fprintf(PMF1d,"%6d %6d ",i,j);
                fprintf(PMF1d,PMFmtS,(double)x);
                fprintf(PMF1d,"\n");
                nrec++;
            }
            if (++j > PMN) {
                i++;
                j = 1;
            }
        }
        else {
            fprintf(PMF1d,PMFmtS,(double)x);
            if (++nc >= PMN) {
                fprintf(PMF1d,"\n");
                nrec++;
                nc = 0;
            }
        }
    }
    printf1("\n%d records written to: %s\n",nrec,PMF1dName);
    err = 0;

RUCIFin:
    p_clean();
    return(err);
}





