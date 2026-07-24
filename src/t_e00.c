/****************************************************************************/
/*  t_e00                                                                   */
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
#include "t_gf.h"
#include "t_sd.h"

/*  functions in t_e00.c */

int sde00(void); 
void e00_err(int n);
char *e00_int(char *p,int *n);
char *e00_dbl(char *p,double *x,double *y,int n);
char *skip_dbl1(char *p);
void e00_bbarc(double x,double y,int first);
void e00_bblab(double x,double y,int first);
int e00_arc(int opt,int dflag,int *ne,int *np,double *xf,double *yf,int *cn,
    int *pn,int *ptr);
int e00_cnt(int dflag,int *ne);
int e00_lab(int opt,int dflag,int *ne,int *cn,double *xf,double *yf);
int e00_pal(int opt,int dflag,int *ne,int *na,int *arcs);
int e00_tol(int dflag,int *ne);
int e00_sin(int dflag,int *ne);
int e00_log(int dflag,int *ne);
int e00_prj(int dflag,int *ne);
int e00_ifo(int dflag,int *ne);
int e00_tx6(int dflag,int *ne);
int e00_tx7(int dflag,int *ne);
int e00_bnd(void);
int e00_prj_prn(void);
int e00_points(void);
int e00_lines(void);
int e00_polygons(void);
int e00_find(int u,int n,int *nlist);
int e00_poly_prn(int id,int id1,int id2,int nn,int n,int *idx,int *ptr,
    int *np,int nv,int nvd);
int e00_var_read(int opt,int ne);
void e00_var_prn(int nva);
int e00_var_buf(int nva,int irec);
void e00_dtda(int opt,int noc,int nv,int nvd);

/*--------------------------------------------------------------------------*/
/*  Global variables                                                        */

FILE *E00Fd;
int E00Fd_Open = 0;         /* set if file opened                           */
int E00DFlag = 0;           /* 2 if single, 3 if double precision           */

#define E00BufL 2000
char E00Buf[E00BufL + 2];      /* global read buffer                           */

long ARCFPtr = -1L;         /* file pointer for ARC section                 */
long CNTFPtr = -1L;         /* file pointer for CNT section                 */
long LABFPtr = -1L;         /* file pointer for LAB section                 */
long PALFPtr = -1L;         /* file pointer for PAL section                 */
long TOLFPtr = -1L;         /* file pointer for TOL section                 */
long SINFPtr = -1L;         /* file pointer for SIN section                 */
long LOGFPtr = -1L;         /* file pointer for LOG section                 */
long PRJFPtr = -1L;         /* file pointer for PRJ section                 */
long IFOFPtr = -1L;         /* file pointer for IFO section                 */
long TX6FPtr = -1L;         /* file pointer for TX6 section                 */
long TX7FPtr = -1L;         /* file pointer for TX7 section                 */

long BNDFPtr = -1L;         /* file pointer for BND section                 */
long AATFPtr = -1L;         /* file pointer for AAT section                 */
long PATFPtr = -1L;         /* file pointer for PAT section                 */

int NARC = 0;               /* number of entries in ARC                     */
int NARCP = 0;              /* number of points in ARC                      */
int NPAL = 0;               /* number of entries in PAL                     */
int NPALA = 0;              /* number of arcs in PAL                        */
int NCNT = 0;               /* number of entries in CNT                     */
int NLAB = 0;               /* number of entries in LAB                     */

double BNDXMin = 0.0;       /* bounding box from BND section                */
double BNDXMax = 0.0;
double BNDYMin = 0.0;
double BNDYMax = 0.0;

double ARCXMin = 0.0;       /* bounding box from ARC section                */
double ARCXMax = 0.0;
double ARCYMin = 0.0;
double ARCYMax = 0.0;

double LABXMin = 0.0;       /* bounding box from LAB section                */
double LABXMax = 0.0;
double LABYMin = 0.0;
double LABYMax = 0.0;

struct E00ATT {             /* structure for E00 attributes                 */
    char Name[17];
    int  Valid;
    int  Len;
    int  Typ;
    int  Fmt1;
    int  Fmt2; 
} *ATT;
int ATTN = 0;
int ATTLen = 0;

/* ------------------------------------------------------------------------ */
/*  sde00   Convert e00 file to spatial data file.                          */
/*                                                                          */
/*          sde00(                                                          */
/*              df=...,         output file (optional)                      */
/*              nc=...,         option for polygon coverages, def. 0        */
/*                              0 = drop first (universal) polygon          */
/*                              1 = keep this polygon                       */
/*              attr=...,       option for attributes, def. 0               */
/*                              0 = drop all attribute variables            */  
/*                              1 = drop only default variables             */
/*                              2 = include all variables                   */
/*              fmt=...,        floating point print format, def. 12.6      */
/*              dtda=...,       name of TDA description file                */
/*          ) = name_of_e00_file;                                           */
/*                                                                          */
/*  The command reads an uncompressed E00 ASCII file. The supported         */
/*  coverages are:                                                          */
/*                                                                          */
/*  1   point coverage      if no ARC section, but a LAB section            */
/*  2   line coverage       if ARC section, but no PAL section              */
/*  3   polygon coverage    if ARC section and PAL section                  */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int sde00(void)  
{
    int err,cflag,dflag,r,ne,np,n,patnv,patnva,patnrec,nrec;
    long fptr;

    E00DFlag = E00Fd_Open = 0;
    nrec = patnv = 0;
    err = -1;

    if (check_cmd(0))
        return(-1);
                    
    printf1("Reading an e00 file. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,8,1)) {     /* get parameters */
        goto SDE00Fin;
    }   
    if (PMNC < 0 || PMNC > 1)
        PMNC = 0;

    if (PMAttr > 2)
        PMAttr = 0;

    if (PMFmtF == 0)
        pmfmt(12,6);

    printf1("File: %s\n",PMRHSTR);

    if (!(E00Fd = fopen(PMRHSTR,OPEN_RD))) {   
        printf1("Error: can't open the file.\n");         
        goto SDE00Fin;
    }
    E00Fd_Open = 1;

    if (!fgets(E00Buf,E00BufL,E00Fd)) {
        printf1("Error: can't read the file.\n");
        goto SDE00Fin;
    }
    if (sscanf(E00Buf,"EXP %d",&cflag) != 1) {
        printf1("Error: probabily not an e00 ASCII file.\n");
        goto SDE00Fin;
    }
    if (cflag != 0) {
        printf1("Might be compressed. Can't continue.\n");
        goto SDE00Fin;
    }
    printf1("\nSection  records    entries\n");

    NARC = NPAL = NLAB = NCNT = 0;
    fptr = ftell(E00Fd);

    while (fgets(E00Buf,E00BufL,E00Fd)) {

        if (sscanf(E00Buf,"ARC %d ",&dflag) == 1) {
            printf1("ARC%2d ",dflag);
            if ((r = e00_arc(0,dflag,&ne,&np,AcX,AcY,AcI,AcN,AcK)) < 0) {
                e00_err(r);
                goto SDE00Fin;
            }
            printf1("%10d %10d\n",r,ne);
            NARC = ne;          
            NARCP = np;
            ARCFPtr = fptr;
            if (E00DFlag == 0)
                E00DFlag = dflag;
        }
        else if (sscanf(E00Buf,"CNT %d",&dflag) == 1) {
            printf1("CNT%2d ",dflag);
            if ((r = e00_cnt(dflag,&ne)) < 0) {
                e00_err(r);
                goto SDE00Fin;
            }
            printf1("%10d %10d\n",r,ne);
            NCNT = ne;
            CNTFPtr = fptr;
        }
        else if (sscanf(E00Buf,"LAB %d",&dflag) == 1) {
            printf1("LAB%2d ",dflag);
            if ((r = e00_lab(0,dflag,&ne,AcI,AcX,AcY)) < 0) {
                e00_err(r);
                goto SDE00Fin;
            }
            printf1("%10d %10d\n",r,ne);
            NLAB = ne;
            LABFPtr = fptr;
            if (E00DFlag == 0)
                E00DFlag = dflag;
        }
        else if (sscanf(E00Buf,"PAL %d",&dflag) == 1) {
            printf1("PAL%2d ",dflag);
            if ((r = e00_pal(0,dflag,&ne,&np,AcJ)) < 0) {
                e00_err(r);
                goto SDE00Fin;
            }
            printf1("%10d %10d\n",r,ne);
            NPAL = ne;
            NPALA = np;
            PALFPtr = fptr;
            if (E00DFlag == 0)
                E00DFlag = dflag;
        }
        else if (sscanf(E00Buf,"TOL %d",&dflag) == 1) {
            printf1("TOL%2d ",dflag);
            if ((r = e00_tol(dflag,&ne)) < 0) {
                e00_err(r);
                goto SDE00Fin;
            }
            printf1("%10d %10d\n",r,ne);
            TOLFPtr = fptr;
        }
        else if (sscanf(E00Buf,"SIN %d",&dflag) == 1) {
            printf1("SIN%2d ",dflag);
            if ((r = e00_sin(dflag,&ne)) < 0) {
                e00_err(r);
                goto SDE00Fin;
            }
            printf1("%10d %10d\n",r,ne);
            SINFPtr = fptr;
        }
        else if (sscanf(E00Buf,"LOG %d",&dflag) == 1) {
            printf1("LOG%2d ",dflag);
            if ((r = e00_log(dflag,&ne)) < 0) {
                e00_err(r);
                goto SDE00Fin;
            }
            printf1("%10d %10d\n",r,ne);
            LOGFPtr = fptr;
        }
        else if (sscanf(E00Buf,"PRJ %d",&dflag) == 1) {
            printf1("PRJ%2d ",dflag);
            if ((r = e00_prj(dflag,&ne)) < 0) {
                e00_err(r);
                goto SDE00Fin;
            }
            printf1("%10d %10d\n",r,ne);
            PRJFPtr = fptr;
        }
        else if (sscanf(E00Buf,"IFO %d",&dflag) == 1) {
            printf1("IFO%2d ",dflag);
            if ((r = e00_ifo(dflag,&ne)) < 0) {
                e00_err(r);
                goto SDE00Fin;
            }
            printf1("%10d %10d\n",r,ne);
        }
        else if (sscanf(E00Buf,"TX6 %d",&dflag) == 1) {
            printf1("TX6%2d ",dflag);
            if ((r = e00_tx6(dflag,&ne)) < 0) {
                e00_err(r);
                goto SDE00Fin;
            }
            printf1("%10d %10d\n",r,ne);
        }
        else if (sscanf(E00Buf,"TX7 %d",&dflag) == 1) {
            printf1("TX7%2d ",dflag);
            if ((r = e00_tx7(dflag,&ne)) < 0) {
                e00_err(r);
                goto SDE00Fin;
            }
            printf1("%10d %10d\n",r,ne);
        }
        else if (!(strncmp(E00Buf,"EOS",3)))
            break;
        else {
            printf1("\nUnknown record type: %s\n",E00Buf);
            goto SDE00Fin;
        }
        fptr = ftell(E00Fd);
    }
    if (PRJFPtr >= 0)                   /* print projection */
        e00_prj_prn();

    if (BNDFPtr >= 0 && e00_bnd() == 0) {
        printf1("\nBounding Box from BND section:\n");
        printf1("X: %20.8lf  %20.8lf\n",BNDXMin,BNDXMax);
        printf1("Y: %20.8lf  %20.8lf\n",BNDYMin,BNDYMax);
    }
    if (NARC > 0) {
        printf1("\nBounding Box from ARC section:\n");
        printf1("X: %20.8lf  %20.8lf\n",ARCXMin,ARCXMax);
        printf1("Y: %20.8lf  %20.8lf\n",ARCYMin,ARCYMax);
    }
    else if (NLAB > 0) {
        printf1("\nBounding Box from LAB section:\n");
        printf1("X: %20.8lf  %20.8lf\n",LABXMin,LABXMax);
        printf1("Y: %20.8lf  %20.8lf\n",LABYMin,LABYMax);
    }
    newline();

    if (PMF1Def == 0) {        
        err = 0;
        goto SDE00Fin;
    }
    prnchar('-',LLEN,1);

    /* writing the output file */

    if (NARC > 0) {
        if (NPAL > 0) {
            printf1("Interpretation: polygon coverage.\n");
            if (e00_polygons())
                goto SDE00Fin;
        }
        else {
            printf1("Interpretation: line coverage.\n");
            if (e00_lines())
                goto SDE00Fin;
        }
    }
    else if (NLAB > 0) {
        printf1("Interpretation: point coverage.\n");
        if (e00_points())
            goto SDE00Fin;
    }
    else {
        printf1("Cannot interpret this file.\n");
        goto SDE00Fin; 
    }
    err = 0;

SDE00Fin:
    E00DFlag = 0;
    if (E00Fd_Open) {   
        fclose(E00Fd);
        E00Fd_Open = 0;
    }
    if (ATTN > 0) {
        free((char *)ATT);
        memrq(-ATTN,sizeof(struct E00ATT));
        ATTN = 0;
    }
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  e00_err(n)      Print error message.                                    */

void e00_err(int n)
{
    printf1("\nError %d while reading the input file.\n",n);
    return;
}

/* ------------------------------------------------------------------------ */
/*  e00_int(p,n)    Read integer entry and return pointer to next location. */
/*                  If error, return NULL.                                  */

char *e00_int(char *p,int *n)
{
    p = skip_b(p);
    if (sscanf(p,"%d",n) != 1)   
        return(NULL);
    p = skip_int(p);
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  e00_dbl(p,x,y,n)  Read n (one or two) double entries and return         */
/*                    pointer to next location. If error return NULL.       */

char *e00_dbl(char *p,double *x,double *y,int n)
{
    p = skip_b(p);
    if (sscanf(p,"%lf",x) != 1)   
        return(NULL);

    if (n == 2) {
        p = skip_dbl1(p);
        p = skip_b(p);
        if (sscanf(p,"%lf",y) != 1)   
            return(NULL);
    }
    p = skip_dbl1(p);
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  skip_dbl1(p)                                                            */
/*      It is assumed that p is a pointer to a double or floating point     */
/*      value. The function returns a pointer to the next character after   */
/*      this value.                                                         */

char *skip_dbl1(char *p)
{
    register int pflag = 0;
    register int eflag = 0;

    while (*p && (isdigit((int)*p) || *p == '.' || *p == '-' || *p == '+' ||
                                                   *p == 'e' || *p == 'E')) {
        if (*p == 'e' || *p == 'E') {
            eflag++;
            pflag = 0;
        }
        if (*p == '+' || *p == '-') {
            pflag++;
            if (eflag && pflag >= 2)
                break;
        }
        p++;
    }
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  e00_bbarc(x,y,first)    bounding box: ARC                               */

void e00_bbarc(double x,double y,int first)
{
    if (first) {
        ARCXMin = ARCXMax = x;
        ARCYMin = ARCYMax = y;
        return;
    }
    ARCXMin = dmin(ARCXMin,x);
    ARCXMax = dmax(ARCXMax,x);
    ARCYMin = dmin(ARCYMin,y);
    ARCYMax = dmax(ARCYMax,y);
}

/* ------------------------------------------------------------------------ */
/*  e00_bblab(x,y,first)    bounding box: LAB                               */

void e00_bblab(double x,double y,int first)
{
    if (first) {
        LABXMin = LABXMax = x;
        LABYMin = LABYMax = y;
        return;
    }
    LABXMin = dmin(LABXMin,x);
    LABXMax = dmax(LABXMax,x);
    LABYMin = dmin(LABYMin,y);
    LABYMax = dmax(LABYMax,y);
}

/* ------------------------------------------------------------------------ */
/*  e00_arc(opt,dflag,ne,np,xf,yf,cn,pn,ptr)                                */
/*                                                                          */
/*  Read ARC records. Also calculate bounding box. Section might be empty!  */
/*  dflag = 2 (single precision) or 3 (double precision)                    */
/*                                                                          */
/*  If opt = 0 only read. If opt != 0 save coordinates in x and y, coverage */
/*  number in cn, number of points in pn, and pointer in ptr.               */
/*                                                                          */
/*  Return < 0 if error, or number of records. Return number of entries     */
/*  in ne and number of points in np.                                       */  

int e00_arc(int opt,int dflag,int *ne,int *np,double *xf,double *yf,int *cn,
    int *pn,int *ptr)
{
    register int i,j,k;
    register char *p;
    int m,n,r,first,m1,m2,m3,m4,m5;
    double x,y; 

    first = 1;
    *ne = 0;
    *np = 0;

    if (opt) {
        if (fseek(E00Fd,ARCFPtr,0) || !fgets(E00Buf,E00BufL,E00Fd) ||
            sscanf(E00Buf,"ARC %d",&n) != 1)  
        return(-1);
    }
    k = j = r = 0;
    while (fgets(E00Buf,E00BufL,E00Fd)) {
        r++;
        p = E00Buf;
        if ((p = e00_int(p,&n)) == NULL)
            return(-1);
        if ((p = e00_int(p,&m1)) == NULL)
            return(-1);
        if ((p = e00_int(p,&m2)) == NULL)
            return(-1);
        if ((p = e00_int(p,&m3)) == NULL)
            return(-1);
        if ((p = e00_int(p,&m4)) == NULL)
            return(-1);
        if ((p = e00_int(p,&m5)) == NULL)
            return(-1);
        if ((p = e00_int(p,&m)) == NULL)
            return(-1);

        if (n == -1 && m == 0)  
            return(r);
           
        if (m < 1)
            return(-2);

        *np += m;
        if (opt) {           
            cn[j] = n;      /* save coverage number */
            pn[j] = m;      /* save number of points */
            ptr[j] = k;     /* save pointers */
            j++;
        }
        while (m > 0) {
            if (!fgets(E00Buf,E00BufL,E00Fd))  
                return(-3);
            r++;
            p = E00Buf;
            if ((p = e00_dbl(p,&x,&y,2)) == NULL)  
                return(-3);
            m--;         
            if (opt) {
                xf[k] = x;
                yf[k] = y;
                k++;
            }       
            e00_bbarc(x,y,first);
            first = 0;

            if (m == 0)
                break;

            if (dflag == 2) {           /* single precision */
                if ((p = e00_dbl(p,&x,&y,2)) == NULL)
                    return(-3);
                m--;        
                if (opt) {
                    xf[k] = x;
                    yf[k] = y;
                    k++;
                }       
                e00_bbarc(x,y,0);
            }
        }
        *ne += 1;
    }
    return(-4);
}

/* ------------------------------------------------------------------------ */
/*  e00_cnt()       Read CNT records.                                       */
/*                                                                          */
/*  Return < 0 if error, or number of records. Return number of entries     */
/*  in ne.                                                                  */  

int e00_cnt(int dflag,int *ne)
{
    register char *p;
    int m,n,r,k;
    double x,y; 

    *ne = 0;
    r = 0;
    while (fgets(E00Buf,E00BufL,E00Fd)) {
        r++;
        p = E00Buf;
        if ((p = e00_int(p,&n)) == NULL)
            return(-1);

        if (n == -1)
            return(r);

        if ((p = e00_dbl(p,&x,&y,2)) == NULL)  
            return(-2);

        while (n > 0) {
            if (!(fgets(E00Buf,E00BufL,E00Fd))) 
                return(-3);
            r++;
            p = E00Buf;
            k = 8;
            while (n > 0 && k > 0) {
                if ((p = e00_int(p,&m)) == NULL)
                    return(-4);
                n--;
                k--;
            }
        }
        *ne += 1;
    }
    return(-5);
}

/* ---=-------------------------------------------------------------------- */
/*  e00_lab(opt,dflag,ne,cn,xf,yf)                                          */
/*                                                                          */
/*  Read the LAB section. Save number of entries in ne. Return number of    */
/*  records or a negative number if an error occurred.                      */
/*                                                                          */
/*  If opt != 0 save coverage number in cn and coordinates in xf and yf.    */

int e00_lab(int opt,int dflag,int *ne,int *cn,double *xf,double *yf)
{
    register char *p;
    int k,m,n,r,first;
    double x,y;
          
    first = 1;
    *ne = 0;

    if (opt) {
        if (fseek(E00Fd,LABFPtr,0) || !fgets(E00Buf,E00BufL,E00Fd) ||
            sscanf(E00Buf,"LAB %d",&n) != 1)  
        return(-2);
    }
    k = r = 0;
    while (fgets(E00Buf,E00BufL,E00Fd)) {
        r++;
        p = E00Buf;
        if ((p = e00_int(p,&n)) == NULL)  
            return(-3);
        if ((p = e00_int(p,&m)) == NULL)   
            return(-3);

        if (n == -1 && m == 0)  
            return(r);
 
        if ((p = e00_dbl(p,&x,&y,2)) == NULL)  
            return(-4);

        if (opt) {
            cn[k] = n;
            xf[k] = x;
            yf[k] = y;
            k++;
        }

        e00_bblab(x,y,first);
        first = 0;

        if (!(fgets(E00Buf,E00BufL,E00Fd))) 
            return(-5);
        r++;
        p = E00Buf;
        if ((p = e00_dbl(p,&x,&y,2)) == NULL)  
            return(-6);
          
        e00_bbarc(x,y,0);
                       
        if (dflag == 3) {           /* double precision */
            if (!(fgets(E00Buf,E00BufL,E00Fd))) 
                return(-7);
            r++;
            p = E00Buf;
        }
        if ((p = e00_dbl(p,&x,&y,2)) == NULL)        
            return(-8);
 
        e00_bbarc(x,y,0);
        *ne += 1;
    }
    return(-9);
}

/* ------------------------------------------------------------------------ */
/*  e00_pal(opt,dflag,ne,na,arcs)                                           */
/*                                                                          */
/*  If opt = 0 read the file and return number of records, also return      */
/*  number of polygons in ne and number of arcs in na.                      */
/*                                                                          */
/*  If opt != 0 save arc information in arc[].                              */
/*                                                                          */
/*  Return negative number if an error occurred.                            */

int e00_pal(int opt,int dflag,int *ne,int *na,int *arcs)
{
    register char *p;
    int k,n,r,m,m1,m2,m3,first;
    double x,y;

    first = 1;
    *ne = 0;
    *na = 0;

    if (opt) {
        if (fseek(E00Fd,PALFPtr,0) || !fgets(E00Buf,E00BufL,E00Fd) ||
            sscanf(E00Buf,"PAL %d",&n) != 1)  
        return(-1);
    }
    k = r = 0;
    while (fgets(E00Buf,E00BufL,E00Fd)) {
        r++;
        p = E00Buf;
        if ((p = e00_int(p,&n)) == NULL)
            return(-1);

        if (n == -1) {
            if (dflag == 3) {
                if (!(fgets(E00Buf,E00BufL,E00Fd)))       /* skip extra line */
                    return(-2);
            }
            return(r);
        }
        if (n < 1)
            return(-5);

        if ((p = e00_dbl(p,&x,&y,2)) == NULL)
            return(-3);

        if (dflag == 3) {                   /* double precision */
            if (!(fgets(E00Buf,E00BufL,E00Fd))) 
                return(-2);
            r++;
            p = E00Buf;
        }
        if ((p = e00_dbl(p,&x,&y,2)) == NULL)        
            return(-3);

        *na += n;           /* number of arcs */
        if (opt)
            arcs[k++] = n;

        while (n > 0) {
            if (!(fgets(E00Buf,E00BufL,E00Fd))) 
                return(-3);
            r++;    
            p = E00Buf;
            if ((p = e00_int(p,&m1)) == NULL)
                return(-4);
            if ((p = e00_int(p,&m2)) == NULL)
                return(-4);
            if ((p = e00_int(p,&m3)) == NULL)
                return(-4);
            n--;
            if (opt)
                arcs[k++] = m1;
    
            if (n == 0)
                break;

            if ((p = e00_int(p,&m1)) == NULL)
                return(-4);
            if ((p = e00_int(p,&m2)) == NULL)
                return(-4);
            if ((p = e00_int(p,&m3)) == NULL)
                return(-4);
            n--;
            if (opt)
                arcs[k++] = m1;
        }
        *ne += 1;
    }
    return(-5);
}

/* ------------------------------------------------------------------------ */
/*  e00_tol()       Read TOL records.                                       */
/*                                                                          */
/*  Return < 0 if error, or number of records. Return number of entries     */
/*  in ne.                                                                  */  

int e00_tol(int dflag,int *ne)
{
    register char *p;
    int n,m,r;
    double x;

    *ne = 0;
    r = 0;
    while (fgets(E00Buf,E00BufL,E00Fd)) {
        r++;
        p = E00Buf;
        if ((p = e00_int(p,&n)) == NULL)
            return(-1);
        if ((p = e00_int(p,&m)) == NULL)
            return(-1);

        if (n == -1 && m == 0)
            return(r);

        if ((p = e00_dbl(p,&x,&x,1)) == NULL)
            return(-2);

        *ne += 1;
    }
    return(-3);
}

/* ------------------------------------------------------------------------ */
/*  e00_sin()       Read SIN records.                                       */
/*                                                                          */
/*  Return < 0 if error, or number of records. Return number of entries     */
/*  in ne.                                                                  */  

int e00_sin(int dflag,int *ne)
{
    int r;

    *ne = 0;
    r = 0;
    while (fgets(E00Buf,E00BufL,E00Fd)) {
        r++;
        if (!strncmp(E00Buf,"EOX",3))
            return(r);
        *ne += 1;
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  e00_log()       Read LOG records.                                       */
/*                                                                          */
/*  Return < 0 if error, or number of records. Return number of entries     */
/*  in ne.                                                                  */  

int e00_log(int dflag,int *ne)
{
    int r;

    *ne = 0;
    r = 0;
    while (fgets(E00Buf,E00BufL,E00Fd)) {
        r++;
        if (!strncmp(E00Buf,"EOL",3))
            return(r);
        if (*E00Buf != '~')
            *ne += 1;
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  e00_prj()       Read PRJ records.                                       */
/*                                                                          */
/*  Return < 0 if error, or number of records. Return number of entries     */
/*  in ne.                                                                  */  

int e00_prj(int dflag,int *ne)
{
    int r;

    *ne = 0;
    r = 0;
    while (fgets(E00Buf,E00BufL,E00Fd)) {
        r++;
        if (!strncmp(E00Buf,"EOP",3))
            return(r);
        if (*E00Buf != '~')
            *ne += 1;
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  e00_ifo()       Read IFO records and set the following file pointers    */
/*                  if present:                                             */
/*                                                                          */
/*                  BNDFPtr         .BND                                    */
/*                  AATFPtr         .AAT                                    */
/*                  PATFPtr         .PAT                                    */
/*                                                                          */
/*  Return < 0 if error, or number of records. Return number of entries     */
/*  in ne.                                                                  */  

int e00_ifo(int dflag,int *ne)
{
    register char *p;
    int r;
    long fptr;

    BNDFPtr = AATFPtr = PATFPtr = -1L;

    *ne = 0;
    r = 0;

    fptr = ftell(E00Fd);
    while (fgets(E00Buf,E00BufL,E00Fd)) {
        r++;
        if (!strncmp(E00Buf,"EOI",3))
            return(r);
        p = E00Buf;
        while (*p) {
            if (*p == '.')
                break;
            p++;
        }
        if (!strncmp(p,".BND",4))
            BNDFPtr = fptr;
        else if (!strncmp(p,".AAT",4))
            AATFPtr = fptr;
        else if (!strncmp(p,".PAT",4))
            PATFPtr = fptr;

        *ne += 1;
        fptr = ftell(E00Fd);
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  e00_tx6()       Read TX6 records.                                       */
/*                                                                          */
/*  Return < 0 if error, or number of records. Return number of entries     */
/*  in ne.                                                                  */  

int e00_tx6(int dflag,int *ne)
{
    int r;

    *ne = 0;
    r = 0;
    while (fgets(E00Buf,E00BufL,E00Fd)) {
        r++;
        if (!strncmp(E00Buf,"JABBERWOCKY",11))
            return(r);
        *ne += 1;
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  e00_tx7()       Read TX7 records.                                       */
/*                                                                          */
/*  Return < 0 if error, or number of records. Return number of entries     */
/*  in ne.                                                                  */  

int e00_tx7(int dflag,int *ne)
{
    int r;

    *ne = 0;
    r = 0;
    while (fgets(E00Buf,E00BufL,E00Fd)) {
        r++;
        if (!strncmp(E00Buf,"JABBERWOCKY",11))
            return(r);
        *ne += 1;
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  e00_bnd()   Get bounding box from BND section.                          */
/*                                                                          */

int e00_bnd(void)
{
    register int i;
    register char *p;
    int nv,nva,nl;

    if (fseek(E00Fd,BNDFPtr,0))  
        return(-1);     

    if (!fgets(E00Buf,E00BufL,E00Fd))  
        return(-1);          

    p = E00Buf + 34;
    if (sscanf(p,"%4d%4d%4d",&nv,&nva,&nl) != 3)
        return(-1);

    for (i = 0; i < 5; ++i) {
        if (!fgets(E00Buf,E00BufL,E00Fd))  
            return(-1);        
    }
    if (nl == 32) {       /* double, one more line */
        if (!fgets(E00Buf + 80,E00BufL,E00Fd))  
            return(-1);          
    }
    p = E00Buf;
    if ((p = e00_dbl(p,&BNDXMin,&BNDYMin,2)) == NULL)        
        return(-1);       
    if ((p = e00_dbl(p,&BNDXMax,&BNDYMax,2)) == NULL)        
        return(-1);         
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  e00_prj_prn()   Print contents of projection section.                   */
/*                                                                          */

int e00_prj_prn(void)
{
    newline();
    if (fseek(E00Fd,PRJFPtr,0) || !fgets(E00Buf,E00BufL,E00Fd) || strncmp(E00Buf,"PRJ",3)) {
        printf1("Error: can't seek to PRJ section.\n");
        return(-1);
    }
    while (fgets(E00Buf,E00BufL,E00Fd)) {
        if (!strncmp(E00Buf,"EOP",3))
            return(0);

        if (*E00Buf != '~')
            printf1("%s",E00Buf);
    }
    return(-1);
}

/* -##--------------------------------------------------------------------- */
/*  e00_points()                                                            */
/*                                                                          */
/*  Create the output file based on a point coverage interpretation.        */
/*  NLAB  = number of entries in the LAB section.                           */
/*                                                                          */
/*  Return 0 if successful, or -1 if an error occurred.                     */

int e00_points(void)
{
    register int i,j,k;
    int r,nrec,id,n,nv,nvd;
    char *p;

    nrec = 0;
    printf1("Number of points: %d\n",NLAB);

    nv = -1;  
    nvd = 0;
    if (PMAttr > 0) {                     /* read variables from PAT section */
        if (PATFPtr >= 0) {     
            printf1("\nReading the PAT section.\n");
            nv = e00_var_read(1,NLAB);                 
            if (nv > 0)
                e00_var_prn(nv);
            else {
                PATFPtr = -1;
                nv = -1;
            }
        }
        else {
            printf1("Cannot find a PAT section.\n");
            nv = -1;
        }
    }
    if (PMAttr == 1) {
        nvd = 4;   
        if (nv <= nvd)
            nv = -1;
    }
    printf1("\nCreating the output file.\n");

    /* save point coordinates in AcX[i] and AcY[i], i = 0,...,NLAB - 1.
       coverage number in AcI[j]. */

    if (alloc_acx(NLAB + 1))
        return(-1);    
    if (alloc_acy(NLAB + 1))
        return(-1);    
    if (alloc_aci(NLAB + 1))
        return(-1);    
      
    if ((r = e00_lab(1,E00DFlag,&n,AcI,AcX,AcY)) < 0) {
        printf1("Error (%d) in reading the LAB section.\n",r);
        return(-1);       
    }
    id = 0;
    k = 0;
    for (j = 0; j < NLAB; ++j) {
        id++;
        fprintf(PMF1d,"%8d 1 1 ",id);
        if (nv > 0) {
            if (e00_var_buf(nv,k))
                return(-1);
            p = E00Buf;
            for (i = 0; i < nvd; ++i)  
                p += ATT[i].Len;
            fprintf(PMF1d,"%s",p);
        }
        fprintf(PMF1d,"\n");
        nrec++;

        fprintf(PMF1d,PMFmtS,AcX[k]);
        fprintf(PMF1d,PMFmtS,AcY[k]);
        fprintf(PMF1d,"\n");
        nrec++;
        k++;
    }
    printf1("%d records written to: %s\n",nrec,PMF1dName);
    if (PMTDAFDef)                  /* write TDA description file */
        e00_dtda(1,id,nv,nvd); 
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  e00_lines()                                                             */
/*                                                                          */
/*  Create the output file based on a line coverage interpretation.         */
/*  NARC  = number of entries in the ARC section.                           */
/*  NARCP = number of points in the ARC section.                            */
/*                                                                          */
/*  Return 0 if successful, or -1 if an error occurred.                     */

int e00_lines(void)
{
    register int i,j,k;
    int r,nrec,n,c,id,nv,nvd;
    char *p;

    nrec = 0;
    printf1("Number of lines: %d\n",NARC);
    printf1("Number of coordinates: %d\n",NARCP);

    nv = -1;  
    nvd = 0;
    if (PMAttr > 0) {                     /* read variables from AAT section */
        if (AATFPtr >= 0) {     
            printf1("\nReading the AAT section.\n");
            nv = e00_var_read(0,NARC);                 
            if (nv > 0)
                e00_var_prn(nv);
            else {
                AATFPtr = -1;
                nv = -1;
            }
        }
        else {
            printf1("Cannot find a AAT section.\n");
            nv = -1;
        }
    }
    if (PMAttr == 1) {
        nvd = 7;   
        if (nv <= nvd)
            nv = -1;
    }
    printf1("\nCreating the output file.\n");

    /* save line coordinates in AcX[i] and AcY[i], i = 0,...,NARC - 1.
       coverage number in AcI[j], number of points in line in AcN[j],
       j = 0,...,NARCP - 1, and pointer in AcK. */

    if (alloc_acx(NARCP + 1))
        return(-1);    
    if (alloc_acy(NARCP + 1))
        return(-1);    
    if (alloc_aci(NARC + 1))
        return(-1);    
    if (alloc_acn(NARC + 1))
        return(-1);    
    if (alloc_ack(NARC + 1))
        return(-1);    
      
    if ((r = e00_arc(1,E00DFlag,&n,&c,AcX,AcY,AcI,AcN,AcK)) < 0) {
        printf1("Error in reading the ARC section.\n");
        return(-1);       
    }
    id = 0;
    k = 0;
    for (j = 0; j < NARC; ++j) {
        id++;
        n = AcN[j];
        fprintf(PMF1d,"%8d 2 %7d ",id,n);
        if (nv > 0) {
            if (e00_var_buf(nv,id - 1))
                return(-1);
            p = E00Buf;
            for (i = 0; i < nvd; ++i)  
                p += ATT[i].Len;
            fprintf(PMF1d,"%s",p);
        }
        fprintf(PMF1d,"\n");
        nrec++;

        for (i = 0; i < n; ++i) {
            fprintf(PMF1d,PMFmtS,AcX[k]);
            fprintf(PMF1d,PMFmtS,AcY[k]);
            fprintf(PMF1d,"\n");
            nrec++;
            k++;
        }
    }
    printf1("%d records written to: %s\n",nrec,PMF1dName);
    if (PMTDAFDef)                  /* write TDA description file */
        e00_dtda(2,id,nv,nvd); 
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  e00_polygons()                                                          */
/*                                                                          */
/*  Create the output file based on a polygon coverage interpretation.      */
/*  NARC  = number of entries in the ARC section.                           */
/*  NARCP = number of points in the ARC section.                            */
/*                                                                          */
/*  If PMNC = 0 drop the first (universal) polygon. Print attribute         */
/*  variables according to PMAttr.                                          */
/*                                                                          */
/*  Return 0 if successful, or -1 if an error occurred.                     */

int e00_polygons(void)
{
    register int i,k,l;
    int r,nrec,n,c,na,nn,m,u,us,id,id1,id2,keep,nw,nw1,nv,nvd;

    nrec = 0;
    printf1("Number of polygons: %d\n",NPAL);
    printf1("Number of arcs: %d\n",NPALA);

    nv = -1;  
    nvd = 0;
    if (PMAttr > 0) {                     /* read variables from PAT section */
        if (PATFPtr >= 0) {     
            printf1("\nReading the PAT section.\n");
            nv = e00_var_read(1,NPAL);                 
            if (nv > 0)
                e00_var_prn(nv);
            else {
                PATFPtr = -1;
                nv = -1;
            }
        }
        else {
            printf1("Cannot find a PAT section.\n");
            nv = -1;
        }
    }
    if (PMAttr == 1) {
        nvd = 4;   
        if (nv <= nvd)
            nv = -1;
    }
    printf1("\nCreating the output file.\n");

    /* save polygon coordinates in AcX[i] and AcY[i], i = 0,...,NARC - 1.
       coverage number in AcI[j], number of points in line in AcN[j],
       j = 0,...,NARCP - 1, and pointer in AcK. */

    if (alloc_acx(NARCP + 1))
        return(-1);    
    if (alloc_acy(NARCP + 1))
        return(-1);    
    if (alloc_aci(NARC + 1))
        return(-1);    
    if (alloc_acn(NARC + 1))
        return(-1);    
    if (alloc_ack(NARC + 1))
        return(-1);    
      
    if ((r = e00_arc(1,E00DFlag,&n,&c,AcX,AcY,AcI,AcN,AcK)) < 0) {
        printf1("Error in reading the ARC section.\n");
        return(-1);       
    }

    /* save arc numbers from PAL section in AcJ */

    na = NPAL + NPALA;
    if (alloc_acj(na + 1))
        return(-1);    

    if ((r = e00_pal(1,E00DFlag,&n,&c,AcJ)) < 0) {
        e00_err(r);
        return(-1);    
    }

    /* check whether the entries in AcI are in ascending order */

    for (i = 1; i < NARC; ++i) {
        if (AcI[i] <= AcI[i - 1]) {
            printf1("Error: coverage numbers not in ascending order.\n");
            return(-1);
        }
    }

    /* write polygons, note that e00-polygons may consist of several parts */
     
    if (alloc_acm(NARC + 1))
        return(-1);    

    nw1 = nw = id = id1 = 0;
    if (PMNC == 0)
        keep = 0;
    else
        keep = 1;
    i = 0;
    while (i < na) {                
        id1++;
        id2 = 0;
        nn = 0;
        n = AcJ[i++];
        m = 0;
        for (k = 0; k < n; ++k) {
            u = AcJ[i + k];           
            us = 0;
            if (u < 0) {
                u = -u;
                us = 1;
            }
            if (u == 0) {
                if (m > 0 && keep) {
                    if (nn < 3)
                        nw1++;
                    else {
                        id++;
                        id2++;
                        r = e00_poly_prn(id,id1,id2,nn,m,AcM,AcK,AcN,nv,nvd);
                        if (r < 0)
                            return(-1);
                        nrec += r;
                    }
                    nn = m = 0;
                }
            }
            else {
                l = e00_find(u,NARC,AcI);
                if (l < 0)  
                    nw++;                          
                if (l >= 0) {
                    nn += AcN[l];
                    l++;
                    if (us)
                        l = -l;
                    AcM[m++] = l;
                }
            }
        }
        if (m > 0 && keep) {
            if (nn < 3)
                nw1++;
            else {
                id++;
                id2++;
                r = e00_poly_prn(id,id1,id2,nn,m,AcM,AcK,AcN,nv,nvd);
                if (r < 0)
                    return(-1);
                nrec += r;
            }
            nn = m = 0;
        }
        i += n;
        keep = 1;
    }
    printf1("Number of polygons in output file: %d\n",id);
    printf1("%d records written to: %s\n",nrec,PMF1dName);
    if (nw > 0)
        printf1("Warning: could not find %d arcs from PAL section.\n",nw);
    if (nw1 > 0)
        printf1("Warning: skipped %d polygon(s) with less than 3 points.\n",nw1);

    if (PMTDAFDef)                  /* write TDA description file */
        e00_dtda(3,id,nv,nvd); 
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  e00_find(u,n,nlist)                                                     */
/*                                                                          */
/*  Return index of u in list nlist[i], i = 0,...,n-1.                      */
/*  or -1 if not found.                                                     */
 
int e00_find(int u,int n,int *nlist)
{
    register int i,l,r,k;

    l = 0;
    r = n - 1;
    while (l <= r) {
        k = (l + r) / 2;
        if (u == nlist[k])
            return(k);
        if (u < nlist[k])
            r = k - 1;
        else
            l = k + 1;
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  e00_poly_prn(id,id1,id2,nn,n,idx,ptr,np,nv,nvd)                         */
/*                                                                          */
/*  Print the current polygon. Return number of records written to the      */
/*  output file, or -1 if error.                                            */
/*                                                                          */
/*  If nv > 0 add variables from PAT section. This depends on PMAttr.       */
/*  If nvd > 0 drop first nvd variables.                                    */
 
int e00_poly_prn(int id,int id1,int id2,int nn,int n,int *idx,int *ptr,
    int *np,int nv,int nvd)
{
    register int i,j,k,l;
    int m,s,nrec;
    char *p;
       
    fprintf(PMF1d,"%8d 3 %7d %8d %6d ",id,nn,id1,id2);
    if (nv > 0) {
        if (e00_var_buf(nv,id1 - 1))
            return(-1);
        p = E00Buf;
        for (i = 0; i < nvd; ++i)  
            p += ATT[i].Len;
        fprintf(PMF1d,"%s",p);
    }
    fprintf(PMF1d,"\n");
    nrec = 1;
        
    for (i = 0; i < n; ++i) {
        j = idx[i];
        s = 0;
        if (j < 0) {
            j = -j;
            s = 1;
        }
        j--;
        k = ptr[j];
        m = np[j];
        if (s == 0) {
            for (l = 0; l < m; ++l) {
                fprintf(PMF1d,PMFmtS,(double)AcX[k + l]);
                fprintf(PMF1d,PMFmtS,(double)AcY[k + l]);
                fprintf(PMF1d,"\n");
            }
        }
        else {
            for (l = m - 1; l >= 0; --l) {
                fprintf(PMF1d,PMFmtS,(double)AcX[k + l]);
                fprintf(PMF1d,PMFmtS,(double)AcY[k + l]);
                fprintf(PMF1d,"\n");
            }
        }
        nrec += m;
    }
    return(nrec);
}

/* ------------------------------------------------------------------------ */
/*  e00_var_read(opt,ne)                                                    */
/*                                                                          */
/*  If opt=0 read AAT section, else read PAT section. ne is the expected    */
/*  number of entries.                                                      */
/*  Save information about variables in E00ATT.                             */
/*  Return number of valid entries, or -1 if error.                         */
/*                                                                          */  
/*  Also create the array AcS with file pointers.                           */

int e00_var_read(int opt,int ne)
{
    register int i,j;
    int l,l1,n,nv,nva,nl,nrec,as,sp,nf,np,nt,nvl,vflag;
    long fptr;
    char *p,*q;

    if (opt == 0) {
        if (fseek(E00Fd,AATFPtr,0) || !fgets(E00Buf,E00BufL,E00Fd)) {
            printf1("Cannot read AAT section. Will be ignored.\n");
            return(-1);
        }
    }
    else {
        if (fseek(E00Fd,PATFPtr,0) || !fgets(E00Buf,E00BufL,E00Fd)) {
            printf1("Cannot read PAT section. Will be ignored.\n");
            return(-1);
        }
    }        

    /* nv  = number of valid attributes */
    /* nva = total number of attributes */
    /* nl  = length */
    /* nrec = number of records */

    p = E00Buf + 34;
    if (sscanf(p,"%4d%4d%4d%10d",&nv,&nva,&nl,&nrec) != 4 || nv < 1 ||  
        nva < 1 || nl < 1 || nrec < 1) {
        printf1("Cannot read PAT section. Will be ignored.\n");
        return(-1);
    }
    if (nrec != ne) {
        printf1("Inconsistent information in the AAT or PAT section.\n");
        printf1("PAT section will be ignored.\n");
        return(-1);
    }
    printf1("Number of valid entries: %d\n",nv);

    if (!(ATT = (struct E00ATT *)calloc(nva,sizeof(struct E00ATT)))) {
        p_err(-2,1);
        return(-1);
    }
    ATTN = nva;
    memrq(ATTN,sizeof(struct E00ATT));

    /* save variable descriptions in E00ATT */

    for (i = 0; i < nva; ++i) {

        ATT[i].Valid = 0;

        if (!fgets(E00Buf,E00BufL,E00Fd))  
            return(-1);
        p = E00Buf;
        q = ATT[i].Name;
        for (j = 0; j < 16; ++j) {
            if (*p == ' ')
                break;
            *q = *p++;
            if (*q == '#')
                *q = '_';
            q++;
        }
        *q = '\0';
        if (sscanf(E00Buf + 16,"%3d",&as) != 1)
            return(-1);
        if (sscanf(E00Buf + 21,"%4d",&sp) != 1)
            return(-1);
        if (sscanf(E00Buf + 28,"%4d",&nf) != 1)
            return(-1);
        if (nf < 0)
            nf = 0;
        ATT[i].Fmt1 = nf;

        if (sscanf(E00Buf + 32,"%2d",&np) != 1)
            return(-1);
        if (np < 0)
            np = 0;
        ATT[i].Fmt2 = np;
                        
        if (sscanf(E00Buf + 34,"%3d",&nt) != 1)
            return(-1);
        ATT[i].Typ = nt;

        if (sscanf(E00Buf + 66,"%3d",&nvl) != 1)
            return(-1);
 
        if (nvl > 0) {               /* only valid entries */
            vflag = 1;
            if (nt == 60 && as == 4)
                ATT[i].Len = 14;
            else if (nt == 60 && as == 8)
                ATT[i].Len = 24;
            else if (nt == 50 && as == 2)
                ATT[i].Len = 6;
            else if (nt == 50 && as == 4)
                ATT[i].Len = 11;
            else if (nt == 40)
                ATT[i].Len = 14;
            else if (nt == 30 && as > 0)
                ATT[i].Len = as;
            else if (nt == 20 && as > 0)
                ATT[i].Len = as;
            else if (nt == 10 && as > 0)
                ATT[i].Len = as;
            else  
                vflag = 0;
            if (vflag)
                ATT[i].Valid = 1;
        }
    }

    /* save file pointers in AcS */   

    if (nv < 3 || ATT[2].Valid == 0 || ATT[2].Typ != 50) {
        printf1("Cannot read identification numbers.\n");
        return(-1);
    }
    if (alloc_acs(ne + 1))
        return(-1);    

    ATTLen = 0;
    for (i = 0; i < nva; ++i) {
        if (ATT[i].Valid == 1)
            ATTLen += ATT[i].Len;
    }
    n = ATTLen / 80;
    n = (n + 1) * 80 + 2;
    if (n > E00BufL) {
        printf1("Exceeded maximal buffer length.\n");
        return(-1);
    }
    for (i = 0; i < nrec; ++i) {
        AcS[i] = (int)ftell(E00Fd);
        p = E00Buf;
        l = ATTLen;
        while (l > 0) {
            if (!fgets(p,l + 3,E00Fd))  
                return(-1);

            n = imin(l,80);
            for (j = 0; j < n; ++j) {
                if (*p == LF || *p == CR)
                    break;
                p++;
            }
            for (; j < n; ++j)
                *p++ = ' ';
            *p = '\0';
            l -= 80;
        }
        p = E00Buf;
        /****************
        printf("%s\n",p);
        ****************/
        p = E00Buf + ATT[0].Len + ATT[1].Len;
        *(p + ATT[2].Len) = '\0';
        q = skip_b(p);
        if (sscanf(q,"%d",&n) != 1) { 
            printf1("Cannot read identification numbers.\n");
            return(-1);
        }                            
    }
    return(nv);
}

/* ------------------------------------------------------------------------ */
/*  e00_var_prn(nva)    Print variable definitions from E00ATT to           */
/*                      standard output.                                    */

void e00_var_prn(int nva)
{
    register int i;
    int len,typ,n,m;   

    len = 8;
    for (i = 0; i < nva; ++i) {
        if (ATT[i].Valid == 0)
            continue;
        len = imax(len,strlen(ATT[i].Name));
    }
    printf1("\nIdx  Variable   ");
    prnchar(' ',len - 8,0);
    printf1("Type       Format\n");
    prnchar('-',len + 25,1);

    for (i = 0; i < nva; ++i) {
        if (ATT[i].Valid == 0)
            continue;

        printf1("%3d  %s   ",i + 1,ATT[i].Name);
        prnchar(' ',len - strlen(ATT[i].Name),0);
        typ = ATT[i].Typ;
        if (typ == 10 || typ == 20) 
            printf1("string      %d\n",ATT[i].Len);
        else {
            if (typ == 60 || typ == 40) {       /* floating point */
                m = ATT[i].Fmt2; 
                n = ATT[i].Fmt1 + m + 1;
            }
            else {                              /* integer */
                m = 0;
                n = ATT[i].Fmt1;                      
            }
            printf1("numerical   %d.%d\n",n,m);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  e00_var_buf(nva,irec)                                                   */
/*                                                                          */
/*  Create in E00Buf the variables for                                      */
/*  Return 0 if OK, -e if error.                                            */

int e00_var_buf(int nva,int irec)
{
    register int j;
    int n,l;
    char *p;
 
    if (fseek(E00Fd,(long)AcS[irec],0))  
        goto E00VWErr;
       
    n = ATTLen / 80;
    n = (n + 1) * 80 + 2;
         
    p = E00Buf;
    l = ATTLen;
    while (l > 0) {
        if (!fgets(p,l + 3,E00Fd))  
            goto E00VWErr;
                                     
        n = imin(l,80);
        for (j = 0; j < n; ++j) {
            if (*p == LF || *p == CR)
                break;
            p++;
        }
        for (; j < n; ++j)
            *p++ = ' ';
        *p = '\0';
        l -= 80;
    }
    return(0);

E00VWErr:
    printf1("Cannot read variables in AAT or PAT section.\n");
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  e00_dtda(opt,noc,nv,nvd)    Create TDA description file.                */
/*                                                                          */
/*  opt = 1 : point coverage                                                */
/*  opt = 2 : line coverage                                                 */
/*  opt = 3 : polygon coverage                                              */

void e00_dtda(int opt,int noc,int nv,int nvd)
{
    register int i,j,k,l;
    int n,m,typ,len;
    char *name;

    if (PMTDAFDef == 0) 
        return;

    fprintf(PMTDAFd,"sdnvar(\n");
    fprintf(PMTDAFd,"  dfile = %s,\n",PMF1dName);
    fprintf(PMTDAFd,"  noc = %d,\n",noc);

    if (opt == 1) {
        fprintf(PMTDAFd,"  ffmt=c1(1-8),c2(10),c3(12),");
        n = 3; j = 4; l = 14;
    }
    else if (opt == 2) {
        fprintf(PMTDAFd,"  ffmt=c1(1-8),c2(10),c3(12-18),");
        n = 3; j = 4; l = 20;
    }
    else if (opt == 3) {
        fprintf(PMTDAFd,"  ffmt=c1(1-8),c2(10),c3(12-18),c4(20-27),c5(29-34),");
        n = 5; j = 6; l = 36;
    }
    for (i = 0; i < nv; ++i) {
        if (i < nvd)
            continue;

        k = l + ATT[i].Len;
        if (ATT[i].Typ != 10 && ATT[i].Typ != 20) {
            if (++n > 5) {
                fprintf(PMTDAFd,"\n       ");
                n = 1;
            }
            fprintf(PMTDAFd,"c%d(%d-%d),",j,l,k - 1);
            j++;
        }
        l = k;
    }
    fprintf(PMTDAFd,"\n  SDID  <5>[ 8.0] = c1,\n");
    fprintf(PMTDAFd,"  SDTyp <1>[ 1.0] = c2,\n");
    if (opt == 1)
        fprintf(PMTDAFd,"  SDN   <1>[ 1.0] = c3,\n");
    else
        fprintf(PMTDAFd,"  SDN   <5>[ 7.0] = c3,\n");
    fprintf(PMTDAFd,"  SDPtr <5>[10.0] = rd,\n");

    if (opt == 1) {
        j = 4;
        l = 14;
    }
    else if (opt == 2) {
        j = 4;
        l = 20;
    }
    else if (opt == 3) {
        fprintf(PMTDAFd,"  SDID1 <5>[ 8.0] = c4,\n");
        fprintf(PMTDAFd,"  SDID2 <5>[ 6.0] = c5,\n");
        j = 6;
        l = 36;
    }
    for (i = 0; i < nv; ++i) {
        if (i < nvd)
            continue;

        name = ATT[i].Name;
        typ = ATT[i].Typ;
        len = ATT[i].Len;
        k = l + len;              

        if (typ == 10 || typ == 20) 
            fprintf(PMTDAFd,"  %s = str(%d,%d),\n",name,l,k - 1);
        else {
            if (typ == 60 || typ == 40) {       /* floating point */
                m = ATT[i].Fmt2; 
                n = ATT[i].Fmt1 + m + 1;
            }
            else {                              /* integer */
                m = 0;
                n = ATT[i].Fmt1;                      
            }
            fprintf(PMTDAFd,"  %s [%d.%d] = c%d,\n",name,n,m,j);
            j++;
        }
        l = k;
    }
    fprintf(PMTDAFd,");\n");
    printf1("TDA description written to: %s\n",PMTDAFName);
}



