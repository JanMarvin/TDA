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
#include "tda_context.h"

/*  functions in t_e00.c */

int sde00(TDAContext *ctx);
void e00_err(TDAContext *ctx, int n);
char *e00_int(TDAContext *ctx, char *p,int *n);
char *e00_dbl(TDAContext *ctx, char *p,double *x,double *y,int n);
char *skip_dbl1(TDAContext *ctx, char *p);
void e00_bbarc(TDAContext *ctx, double x,double y,int first);
void e00_bblab(TDAContext *ctx, double x,double y,int first);
int e00_arc(TDAContext *ctx, int opt,int dflag,int *ne,int *np,double *xf,double *yf,int *cn, int *pn,int *ptr);
int e00_cnt(TDAContext *ctx, int dflag,int *ne);
int e00_lab(TDAContext *ctx, int opt,int dflag,int *ne,int *cn,double *xf,double *yf);
int e00_pal(TDAContext *ctx, int opt,int dflag,int *ne,int *na,int *arcs);
int e00_tol(TDAContext *ctx, int dflag,int *ne);
int e00_sin(TDAContext *ctx, int dflag,int *ne);
int e00_log(TDAContext *ctx, int dflag,int *ne);
int e00_prj(TDAContext *ctx, int dflag,int *ne);
int e00_ifo(TDAContext *ctx, int dflag,int *ne);
int e00_tx6(TDAContext *ctx, int dflag,int *ne);
int e00_tx7(TDAContext *ctx, int dflag,int *ne);
int e00_bnd(TDAContext *ctx);
int e00_prj_prn(TDAContext *ctx);
int e00_points(TDAContext *ctx);
int e00_lines(TDAContext *ctx);
int e00_polygons(TDAContext *ctx);
int e00_find(TDAContext *ctx, int u,int n,int *nlist);
int e00_poly_prn(TDAContext *ctx, int id,int id1,int id2,int nn,int n,int *idx,int *ptr, int *np,int nv,int nvd);
int e00_var_read(TDAContext *ctx, int opt,int ne);
void e00_var_prn(TDAContext *ctx, int nva);
int e00_var_buf(TDAContext *ctx, int nva,int irec);
void e00_dtda(TDAContext *ctx, int opt,int noc,int nv,int nvd);

/*--------------------------------------------------------------------------*/
/*  Global variables                                                        */


#define E00BufL 2000







struct E00ATT {             /* structure for E00 attributes                 */
    char Name[17];
    int  Valid;
    int  Len;
    int  Typ;
    int  Fmt1;
    int  Fmt2;
} *ATT;

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

int sde00(TDAContext *ctx)
{
    int err,cflag,dflag,r,ne,np;
    long fptr;

    ctx->E00DFlag = ctx->E00Fd_Open = 0;
    err = -1;

    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Reading an e00 file. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,8,1)) {     /* get parameters */
        goto SDE00Fin;
    }
    if (ctx->PMNC < 0 || ctx->PMNC > 1)
        ctx->PMNC = 0;

    if (ctx->PMAttr > 2)
        ctx->PMAttr = 0;

    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 12,6);

    printf1(ctx, "File: %s\n",ctx->PMRHSTR);

    if (!(ctx->E00Fd = fopen(ctx->PMRHSTR,OPEN_RD))) {
        printf1(ctx, "Error: can't open the file.\n");
        goto SDE00Fin;
    }
    ctx->E00Fd_Open = 1;

    if (!fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)) {
        printf1(ctx, "Error: can't read the file.\n");
        goto SDE00Fin;
    }
    if (sscanf(ctx->E00Buf,"EXP %d",&cflag) != 1) {
        printf1(ctx, "Error: probabily not an e00 ASCII file.\n");
        goto SDE00Fin;
    }
    if (cflag != 0) {
        printf1(ctx, "Might be compressed. Can't continue.\n");
        goto SDE00Fin;
    }
    printf1(ctx, "\nSection  records    entries\n");

    ctx->NARC = ctx->NPAL = ctx->NLAB = ctx->NCNT = 0;
    fptr = ftell(ctx->E00Fd);

    while (fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)) {

        if (sscanf(ctx->E00Buf,"ARC %d ",&dflag) == 1) {
            printf1(ctx, "ARC%2d ",dflag);
            if ((r = e00_arc(ctx, 0,dflag,&ne,&np,ctx->AcX,ctx->AcY,ctx->AcI,ctx->AcN,ctx->AcK)) < 0) {
                e00_err(ctx, r);
                goto SDE00Fin;
            }
            printf1(ctx, "%10d %10d\n",r,ne);
            ctx->NARC = ne;
            ctx->NARCP = np;
            ctx->ARCFPtr = fptr;
            if (ctx->E00DFlag == 0)
                ctx->E00DFlag = dflag;
        }
        else if (sscanf(ctx->E00Buf,"CNT %d",&dflag) == 1) {
            printf1(ctx, "CNT%2d ",dflag);
            if ((r = e00_cnt(ctx, dflag,&ne)) < 0) {
                e00_err(ctx, r);
                goto SDE00Fin;
            }
            printf1(ctx, "%10d %10d\n",r,ne);
            ctx->NCNT = ne;
            ctx->CNTFPtr = fptr;
        }
        else if (sscanf(ctx->E00Buf,"LAB %d",&dflag) == 1) {
            printf1(ctx, "LAB%2d ",dflag);
            if ((r = e00_lab(ctx, 0,dflag,&ne,ctx->AcI,ctx->AcX,ctx->AcY)) < 0) {
                e00_err(ctx, r);
                goto SDE00Fin;
            }
            printf1(ctx, "%10d %10d\n",r,ne);
            ctx->NLAB = ne;
            ctx->LABFPtr = fptr;
            if (ctx->E00DFlag == 0)
                ctx->E00DFlag = dflag;
        }
        else if (sscanf(ctx->E00Buf,"PAL %d",&dflag) == 1) {
            printf1(ctx, "PAL%2d ",dflag);
            if ((r = e00_pal(ctx, 0,dflag,&ne,&np,ctx->AcJ)) < 0) {
                e00_err(ctx, r);
                goto SDE00Fin;
            }
            printf1(ctx, "%10d %10d\n",r,ne);
            ctx->NPAL = ne;
            ctx->NPALA = np;
            ctx->PALFPtr = fptr;
            if (ctx->E00DFlag == 0)
                ctx->E00DFlag = dflag;
        }
        else if (sscanf(ctx->E00Buf,"TOL %d",&dflag) == 1) {
            printf1(ctx, "TOL%2d ",dflag);
            if ((r = e00_tol(ctx, dflag,&ne)) < 0) {
                e00_err(ctx, r);
                goto SDE00Fin;
            }
            printf1(ctx, "%10d %10d\n",r,ne);
            ctx->TOLFPtr = fptr;
        }
        else if (sscanf(ctx->E00Buf,"SIN %d",&dflag) == 1) {
            printf1(ctx, "SIN%2d ",dflag);
            if ((r = e00_sin(ctx, dflag,&ne)) < 0) {
                e00_err(ctx, r);
                goto SDE00Fin;
            }
            printf1(ctx, "%10d %10d\n",r,ne);
            ctx->SINFPtr = fptr;
        }
        else if (sscanf(ctx->E00Buf,"LOG %d",&dflag) == 1) {
            printf1(ctx, "LOG%2d ",dflag);
            if ((r = e00_log(ctx, dflag,&ne)) < 0) {
                e00_err(ctx, r);
                goto SDE00Fin;
            }
            printf1(ctx, "%10d %10d\n",r,ne);
            ctx->LOGFPtr = fptr;
        }
        else if (sscanf(ctx->E00Buf,"PRJ %d",&dflag) == 1) {
            printf1(ctx, "PRJ%2d ",dflag);
            if ((r = e00_prj(ctx, dflag,&ne)) < 0) {
                e00_err(ctx, r);
                goto SDE00Fin;
            }
            printf1(ctx, "%10d %10d\n",r,ne);
            ctx->PRJFPtr = fptr;
        }
        else if (sscanf(ctx->E00Buf,"IFO %d",&dflag) == 1) {
            printf1(ctx, "IFO%2d ",dflag);
            if ((r = e00_ifo(ctx, dflag,&ne)) < 0) {
                e00_err(ctx, r);
                goto SDE00Fin;
            }
            printf1(ctx, "%10d %10d\n",r,ne);
        }
        else if (sscanf(ctx->E00Buf,"TX6 %d",&dflag) == 1) {
            printf1(ctx, "TX6%2d ",dflag);
            if ((r = e00_tx6(ctx, dflag,&ne)) < 0) {
                e00_err(ctx, r);
                goto SDE00Fin;
            }
            printf1(ctx, "%10d %10d\n",r,ne);
        }
        else if (sscanf(ctx->E00Buf,"TX7 %d",&dflag) == 1) {
            printf1(ctx, "TX7%2d ",dflag);
            if ((r = e00_tx7(ctx, dflag,&ne)) < 0) {
                e00_err(ctx, r);
                goto SDE00Fin;
            }
            printf1(ctx, "%10d %10d\n",r,ne);
        }
        else if (!(strncmp(ctx->E00Buf,"EOS",3)))
            break;
        else {
            printf1(ctx, "\nUnknown record type: %s\n",ctx->E00Buf);
            goto SDE00Fin;
        }
        fptr = ftell(ctx->E00Fd);
    }
    if (ctx->PRJFPtr >= 0)                   /* print projection */
        e00_prj_prn(ctx);

    if (ctx->BNDFPtr >= 0 && e00_bnd(ctx) == 0) {
        printf1(ctx, "\nBounding Box from BND section:\n");
        printf1(ctx, "X: %20.8lf  %20.8lf\n",ctx->BNDXMin,ctx->BNDXMax);
        printf1(ctx, "Y: %20.8lf  %20.8lf\n",ctx->BNDYMin,ctx->BNDYMax);
    }
    if (ctx->NARC > 0) {
        printf1(ctx, "\nBounding Box from ARC section:\n");
        printf1(ctx, "X: %20.8lf  %20.8lf\n",ctx->ARCXMin,ctx->ARCXMax);
        printf1(ctx, "Y: %20.8lf  %20.8lf\n",ctx->ARCYMin,ctx->ARCYMax);
    }
    else if (ctx->NLAB > 0) {
        printf1(ctx, "\nBounding Box from LAB section:\n");
        printf1(ctx, "X: %20.8lf  %20.8lf\n",ctx->LABXMin,ctx->LABXMax);
        printf1(ctx, "Y: %20.8lf  %20.8lf\n",ctx->LABYMin,ctx->LABYMax);
    }
    newline(ctx);

    if (ctx->PMF1Def == 0) {
        err = 0;
        goto SDE00Fin;
    }
    prnchar(ctx, '-',LLEN,1);

    /* writing the output file */

    if (ctx->NARC > 0) {
        if (ctx->NPAL > 0) {
            printf1(ctx, "Interpretation: polygon coverage.\n");
            if (e00_polygons(ctx))
                goto SDE00Fin;
        }
        else {
            printf1(ctx, "Interpretation: line coverage.\n");
            if (e00_lines(ctx))
                goto SDE00Fin;
        }
    }
    else if (ctx->NLAB > 0) {
        printf1(ctx, "Interpretation: point coverage.\n");
        if (e00_points(ctx))
            goto SDE00Fin;
    }
    else {
        printf1(ctx, "Cannot interpret this file.\n");
        goto SDE00Fin;
    }
    err = 0;

SDE00Fin:
    ctx->E00DFlag = 0;
    if (ctx->E00Fd_Open) {
        fclose(ctx->E00Fd);
        ctx->E00Fd_Open = 0;
    }
    if (ctx->ATTN > 0) {
        free((char *)ATT);
        memrq(ctx, -ctx->ATTN,sizeof(struct E00ATT));
        ctx->ATTN = 0;
    }
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  e00_err(n)      Print error message.                                    */

void e00_err(TDAContext *ctx, int n)
{
    printf1(ctx, "\nError %d while reading the input file.\n",n);
    return;
}

/* ------------------------------------------------------------------------ */
/*  e00_int(p,n)    Read integer entry and return pointer to next location. */
/*                  If error, return NULL.                                  */

char *e00_int(TDAContext *ctx, char *p,int *n)
{
    p = skip_b(ctx, p);
    if (sscanf(p,"%d",n) != 1)
        return(NULL);
    p = skip_int(ctx, p);
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  e00_dbl(p,x,y,n)  Read n (one or two) double entries and return         */
/*                    pointer to next location. If error return NULL.       */

char *e00_dbl(TDAContext *ctx, char *p,double *x,double *y,int n)
{
    p = skip_b(ctx, p);
    if (sscanf(p,"%lf",x) != 1)
        return(NULL);

    if (n == 2) {
        p = skip_dbl1(ctx, p);
        p = skip_b(ctx, p);
        if (sscanf(p,"%lf",y) != 1)
            return(NULL);
    }
    p = skip_dbl1(ctx, p);
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  skip_dbl1(p)                                                            */
/*      It is assumed that p is a pointer to a double or floating point     */
/*      value. The function returns a pointer to the next character after   */
/*      this value.                                                         */

char *skip_dbl1(TDAContext *ctx, char *p)
{
    (void)ctx;        /* unused: the signature is shared */
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

void e00_bbarc(TDAContext *ctx, double x,double y,int first)
{
    if (first) {
        ctx->ARCXMin = ctx->ARCXMax = x;
        ctx->ARCYMin = ctx->ARCYMax = y;
        return;
    }
    ctx->ARCXMin = dmin(ctx, ctx->ARCXMin,x);
    ctx->ARCXMax = dmax(ctx, ctx->ARCXMax,x);
    ctx->ARCYMin = dmin(ctx, ctx->ARCYMin,y);
    ctx->ARCYMax = dmax(ctx, ctx->ARCYMax,y);
}

/* ------------------------------------------------------------------------ */
/*  e00_bblab(x,y,first)    bounding box: LAB                               */

void e00_bblab(TDAContext *ctx, double x,double y,int first)
{
    if (first) {
        ctx->LABXMin = ctx->LABXMax = x;
        ctx->LABYMin = ctx->LABYMax = y;
        return;
    }
    ctx->LABXMin = dmin(ctx, ctx->LABXMin,x);
    ctx->LABXMax = dmax(ctx, ctx->LABXMax,x);
    ctx->LABYMin = dmin(ctx, ctx->LABYMin,y);
    ctx->LABYMax = dmax(ctx, ctx->LABYMax,y);
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

int e00_arc(TDAContext *ctx, int opt,int dflag,int *ne,int *np,double *xf,double *yf,int *cn, int *pn,int *ptr)
{
    register int j,k;
    register char *p;
    int m,n,r,first,m1,m2,m3,m4,m5;
    double x,y;

    first = 1;
    *ne = 0;
    *np = 0;

    if (opt) {
        if (fseek(ctx->E00Fd,ctx->ARCFPtr,0) || !fgets(ctx->E00Buf,E00BufL,ctx->E00Fd) ||
            sscanf(ctx->E00Buf,"ARC %d",&n) != 1)
        return(-1);
    }
    k = j = r = 0;
    while (fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)) {
        r++;
        p = ctx->E00Buf;
        if ((p = e00_int(ctx, p,&n)) == NULL)
            return(-1);
        if ((p = e00_int(ctx, p,&m1)) == NULL)
            return(-1);
        if ((p = e00_int(ctx, p,&m2)) == NULL)
            return(-1);
        if ((p = e00_int(ctx, p,&m3)) == NULL)
            return(-1);
        if ((p = e00_int(ctx, p,&m4)) == NULL)
            return(-1);
        if ((p = e00_int(ctx, p,&m5)) == NULL)
            return(-1);
        if ((p = e00_int(ctx, p,&m)) == NULL)
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
            if (!fgets(ctx->E00Buf,E00BufL,ctx->E00Fd))
                return(-3);
            r++;
            p = ctx->E00Buf;
            if ((p = e00_dbl(ctx, p,&x,&y,2)) == NULL)
                return(-3);
            m--;
            if (opt) {
                xf[k] = x;
                yf[k] = y;
                k++;
            }
            e00_bbarc(ctx, x,y,first);
            first = 0;

            if (m == 0)
                break;

            if (dflag == 2) {           /* single precision */
                if ((p = e00_dbl(ctx, p,&x,&y,2)) == NULL)
                    return(-3);
                m--;
                if (opt) {
                    xf[k] = x;
                    yf[k] = y;
                    k++;
                }
                e00_bbarc(ctx, x,y,0);
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

int e00_cnt(TDAContext *ctx, int dflag,int *ne)
{
    (void)dflag;        /* unused: the signature is shared */
    register char *p;
    int m,n,r,k;
    double x,y;

    *ne = 0;
    r = 0;
    while (fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)) {
        r++;
        p = ctx->E00Buf;
        if ((p = e00_int(ctx, p,&n)) == NULL)
            return(-1);

        if (n == -1)
            return(r);

        if ((p = e00_dbl(ctx, p,&x,&y,2)) == NULL)
            return(-2);

        while (n > 0) {
            if (!(fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)))
                return(-3);
            r++;
            p = ctx->E00Buf;
            k = 8;
            while (n > 0 && k > 0) {
                if ((p = e00_int(ctx, p,&m)) == NULL)
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

int e00_lab(TDAContext *ctx, int opt,int dflag,int *ne,int *cn,double *xf,double *yf)
{
    register char *p;
    int k,m,n,r,first;
    double x,y;

    first = 1;
    *ne = 0;

    if (opt) {
        if (fseek(ctx->E00Fd,ctx->LABFPtr,0) || !fgets(ctx->E00Buf,E00BufL,ctx->E00Fd) ||
            sscanf(ctx->E00Buf,"LAB %d",&n) != 1)
        return(-2);
    }
    k = r = 0;
    while (fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)) {
        r++;
        p = ctx->E00Buf;
        if ((p = e00_int(ctx, p,&n)) == NULL)
            return(-3);
        if ((p = e00_int(ctx, p,&m)) == NULL)
            return(-3);

        if (n == -1 && m == 0)
            return(r);

        if ((p = e00_dbl(ctx, p,&x,&y,2)) == NULL)
            return(-4);

        if (opt) {
            cn[k] = n;
            xf[k] = x;
            yf[k] = y;
            k++;
        }

        e00_bblab(ctx, x,y,first);
        first = 0;

        if (!(fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)))
            return(-5);
        r++;
        p = ctx->E00Buf;
        if ((p = e00_dbl(ctx, p,&x,&y,2)) == NULL)
            return(-6);

        e00_bbarc(ctx, x,y,0);

        if (dflag == 3) {           /* double precision */
            if (!(fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)))
                return(-7);
            r++;
            p = ctx->E00Buf;
        }
        if ((p = e00_dbl(ctx, p,&x,&y,2)) == NULL)
            return(-8);

        e00_bbarc(ctx, x,y,0);
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

int e00_pal(TDAContext *ctx, int opt,int dflag,int *ne,int *na,int *arcs)
{
    register char *p;
    int k,n,r,m1,m2,m3;
    double x,y;

    *ne = 0;
    *na = 0;

    if (opt) {
        if (fseek(ctx->E00Fd,ctx->PALFPtr,0) || !fgets(ctx->E00Buf,E00BufL,ctx->E00Fd) ||
            sscanf(ctx->E00Buf,"PAL %d",&n) != 1)
        return(-1);
    }
    k = r = 0;
    while (fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)) {
        r++;
        p = ctx->E00Buf;
        if ((p = e00_int(ctx, p,&n)) == NULL)
            return(-1);

        if (n == -1) {
            if (dflag == 3) {
                if (!(fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)))       /* skip extra line */
                    return(-2);
            }
            return(r);
        }
        if (n < 1)
            return(-5);

        if ((p = e00_dbl(ctx, p,&x,&y,2)) == NULL)
            return(-3);

        if (dflag == 3) {                   /* double precision */
            if (!(fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)))
                return(-2);
            r++;
            p = ctx->E00Buf;
        }
        if ((p = e00_dbl(ctx, p,&x,&y,2)) == NULL)
            return(-3);

        *na += n;           /* number of arcs */
        if (opt)
            arcs[k++] = n;

        while (n > 0) {
            if (!(fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)))
                return(-3);
            r++;
            p = ctx->E00Buf;
            if ((p = e00_int(ctx, p,&m1)) == NULL)
                return(-4);
            if ((p = e00_int(ctx, p,&m2)) == NULL)
                return(-4);
            if ((p = e00_int(ctx, p,&m3)) == NULL)
                return(-4);
            n--;
            if (opt)
                arcs[k++] = m1;

            if (n == 0)
                break;

            if ((p = e00_int(ctx, p,&m1)) == NULL)
                return(-4);
            if ((p = e00_int(ctx, p,&m2)) == NULL)
                return(-4);
            if ((p = e00_int(ctx, p,&m3)) == NULL)
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

int e00_tol(TDAContext *ctx, int dflag,int *ne)
{
    (void)dflag;        /* unused: the signature is shared */
    register char *p;
    int n,m,r;
    double x;

    *ne = 0;
    r = 0;
    while (fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)) {
        r++;
        p = ctx->E00Buf;
        if ((p = e00_int(ctx, p,&n)) == NULL)
            return(-1);
        if ((p = e00_int(ctx, p,&m)) == NULL)
            return(-1);

        if (n == -1 && m == 0)
            return(r);

        if ((p = e00_dbl(ctx, p,&x,&x,1)) == NULL)
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

int e00_sin(TDAContext *ctx, int dflag,int *ne)
{
    (void)dflag;        /* unused: the signature is shared */
    int r;

    *ne = 0;
    r = 0;
    while (fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)) {
        r++;
        if (!strncmp(ctx->E00Buf,"EOX",3))
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

int e00_log(TDAContext *ctx, int dflag,int *ne)
{
    (void)dflag;        /* unused: the signature is shared */
    int r;

    *ne = 0;
    r = 0;
    while (fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)) {
        r++;
        if (!strncmp(ctx->E00Buf,"EOL",3))
            return(r);
        if (*ctx->E00Buf != '~')
            *ne += 1;
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  e00_prj()       Read PRJ records.                                       */
/*                                                                          */
/*  Return < 0 if error, or number of records. Return number of entries     */
/*  in ne.                                                                  */

int e00_prj(TDAContext *ctx, int dflag,int *ne)
{
    (void)dflag;        /* unused: the signature is shared */
    int r;

    *ne = 0;
    r = 0;
    while (fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)) {
        r++;
        if (!strncmp(ctx->E00Buf,"EOP",3))
            return(r);
        if (*ctx->E00Buf != '~')
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

int e00_ifo(TDAContext *ctx, int dflag,int *ne)
{
    (void)dflag;        /* unused: the signature is shared */
    register char *p;
    int r;
    long fptr;

    ctx->BNDFPtr = ctx->AATFPtr = ctx->PATFPtr = -1L;

    *ne = 0;
    r = 0;

    fptr = ftell(ctx->E00Fd);
    while (fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)) {
        r++;
        if (!strncmp(ctx->E00Buf,"EOI",3))
            return(r);
        p = ctx->E00Buf;
        while (*p) {
            if (*p == '.')
                break;
            p++;
        }
        if (!strncmp(p,".BND",4))
            ctx->BNDFPtr = fptr;
        else if (!strncmp(p,".AAT",4))
            ctx->AATFPtr = fptr;
        else if (!strncmp(p,".PAT",4))
            ctx->PATFPtr = fptr;

        *ne += 1;
        fptr = ftell(ctx->E00Fd);
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  e00_tx6()       Read TX6 records.                                       */
/*                                                                          */
/*  Return < 0 if error, or number of records. Return number of entries     */
/*  in ne.                                                                  */

int e00_tx6(TDAContext *ctx, int dflag,int *ne)
{
    (void)dflag;        /* unused: the signature is shared */
    int r;

    *ne = 0;
    r = 0;
    while (fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)) {
        r++;
        if (!strncmp(ctx->E00Buf,"JABBERWOCKY",11))
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

int e00_tx7(TDAContext *ctx, int dflag,int *ne)
{
    (void)dflag;        /* unused: the signature is shared */
    int r;

    *ne = 0;
    r = 0;
    while (fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)) {
        r++;
        if (!strncmp(ctx->E00Buf,"JABBERWOCKY",11))
            return(r);
        *ne += 1;
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  e00_bnd()   Get bounding box from BND section.                          */
/*                                                                          */

int e00_bnd(TDAContext *ctx)
{
    register int i;
    register char *p;
    int nv,nva,nl;

    if (fseek(ctx->E00Fd,ctx->BNDFPtr,0))
        return(-1);

    if (!fgets(ctx->E00Buf,E00BufL,ctx->E00Fd))
        return(-1);

    p = ctx->E00Buf + 34;
    if (sscanf(p,"%4d%4d%4d",&nv,&nva,&nl) != 3)
        return(-1);

    for (i = 0; i < 5; ++i) {
        if (!fgets(ctx->E00Buf,E00BufL,ctx->E00Fd))
            return(-1);
    }
    if (nl == 32) {       /* double, one more line */
        if (!fgets(ctx->E00Buf + 80,E00BufL - 80,ctx->E00Fd))
            return(-1);
    }
    p = ctx->E00Buf;
    if ((p = e00_dbl(ctx, p,&ctx->BNDXMin,&ctx->BNDYMin,2)) == NULL)
        return(-1);
    if ((p = e00_dbl(ctx, p,&ctx->BNDXMax,&ctx->BNDYMax,2)) == NULL)
        return(-1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  e00_prj_prn()   Print contents of projection section.                   */
/*                                                                          */

int e00_prj_prn(TDAContext *ctx)
{
    newline(ctx);
    if (fseek(ctx->E00Fd,ctx->PRJFPtr,0) || !fgets(ctx->E00Buf,E00BufL,ctx->E00Fd) || strncmp(ctx->E00Buf,"PRJ",3)) {
        printf1(ctx, "Error: can't seek to PRJ section.\n");
        return(-1);
    }
    while (fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)) {
        if (!strncmp(ctx->E00Buf,"EOP",3))
            return(0);

        if (*ctx->E00Buf != '~')
            printf1(ctx, "%s",ctx->E00Buf);
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

int e00_points(TDAContext *ctx)
{
    register int i,j,k;
    int r,nrec,id,n,nv,nvd;
    char *p;

    nrec = 0;
    printf1(ctx, "Number of points: %d\n",ctx->NLAB);

    nv = -1;
    nvd = 0;
    if (ctx->PMAttr > 0) {                     /* read variables from PAT section */
        if (ctx->PATFPtr >= 0) {
            printf1(ctx, "\nReading the PAT section.\n");
            nv = e00_var_read(ctx, 1,ctx->NLAB);
            if (nv > 0)
                e00_var_prn(ctx, nv);
            else {
                ctx->PATFPtr = -1;
                nv = -1;
            }
        }
        else {
            printf1(ctx, "Cannot find a PAT section.\n");
            nv = -1;
        }
    }
    if (ctx->PMAttr == 1) {
        nvd = 4;
        if (nv <= nvd)
            nv = -1;
    }
    printf1(ctx, "\nCreating the output file.\n");

    /* save point coordinates in AcX[i] and AcY[i], i = 0,...,NLAB - 1.
       coverage number in AcI[j]. */

    if (alloc_acx(ctx, ctx->NLAB + 1))
        return(-1);
    if (alloc_acy(ctx, ctx->NLAB + 1))
        return(-1);
    if (alloc_aci(ctx, ctx->NLAB + 1))
        return(-1);

    if ((r = e00_lab(ctx, 1,ctx->E00DFlag,&n,ctx->AcI,ctx->AcX,ctx->AcY)) < 0) {
        printf1(ctx, "Error (%d) in reading the LAB section.\n",r);
        return(-1);
    }
    id = 0;
    k = 0;
    for (j = 0; j < ctx->NLAB; ++j) {
        id++;
        fprintf(ctx->PMF1d,"%8d 1 1 ",id);
        if (nv > 0) {
            if (e00_var_buf(ctx, nv,k))
                return(-1);
            p = ctx->E00Buf;
            for (i = 0; i < nvd; ++i)
                p += ATT[i].Len;
            fprintf(ctx->PMF1d,"%s",p);
        }
        fprintf(ctx->PMF1d,"\n");
        nrec++;

        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[k]);
        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[k]);
        fprintf(ctx->PMF1d,"\n");
        nrec++;
        k++;
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);
    if (ctx->PMTDAFDef)                  /* write TDA description file */
        e00_dtda(ctx, 1,id,nv,nvd);
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

int e00_lines(TDAContext *ctx)
{
    register int i,j,k;
    int r,nrec,n,c,id,nv,nvd;
    char *p;

    nrec = 0;
    printf1(ctx, "Number of lines: %d\n",ctx->NARC);
    printf1(ctx, "Number of coordinates: %d\n",ctx->NARCP);

    nv = -1;
    nvd = 0;
    if (ctx->PMAttr > 0) {                     /* read variables from AAT section */
        if (ctx->AATFPtr >= 0) {
            printf1(ctx, "\nReading the AAT section.\n");
            nv = e00_var_read(ctx, 0,ctx->NARC);
            if (nv > 0)
                e00_var_prn(ctx, nv);
            else {
                ctx->AATFPtr = -1;
                nv = -1;
            }
        }
        else {
            printf1(ctx, "Cannot find a AAT section.\n");
            nv = -1;
        }
    }
    if (ctx->PMAttr == 1) {
        nvd = 7;
        if (nv <= nvd)
            nv = -1;
    }
    printf1(ctx, "\nCreating the output file.\n");

    /* save line coordinates in AcX[i] and AcY[i], i = 0,...,NARC - 1.
       coverage number in AcI[j], number of points in line in AcN[j],
       j = 0,...,NARCP - 1, and pointer in AcK. */

    if (alloc_acx(ctx, ctx->NARCP + 1))
        return(-1);
    if (alloc_acy(ctx, ctx->NARCP + 1))
        return(-1);
    if (alloc_aci(ctx, ctx->NARC + 1))
        return(-1);
    if (alloc_acn(ctx, ctx->NARC + 1))
        return(-1);
    if (alloc_ack(ctx, ctx->NARC + 1))
        return(-1);

    if ((r = e00_arc(ctx, 1,ctx->E00DFlag,&n,&c,ctx->AcX,ctx->AcY,ctx->AcI,ctx->AcN,ctx->AcK)) < 0) {
        printf1(ctx, "Error in reading the ARC section.\n");
        return(-1);
    }
    id = 0;
    k = 0;
    for (j = 0; j < ctx->NARC; ++j) {
        id++;
        n = ctx->AcN[j];
        fprintf(ctx->PMF1d,"%8d 2 %7d ",id,n);
        if (nv > 0) {
            if (e00_var_buf(ctx, nv,id - 1))
                return(-1);
            p = ctx->E00Buf;
            for (i = 0; i < nvd; ++i)
                p += ATT[i].Len;
            fprintf(ctx->PMF1d,"%s",p);
        }
        fprintf(ctx->PMF1d,"\n");
        nrec++;

        for (i = 0; i < n; ++i) {
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[k]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[k]);
            fprintf(ctx->PMF1d,"\n");
            nrec++;
            k++;
        }
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);
    if (ctx->PMTDAFDef)                  /* write TDA description file */
        e00_dtda(ctx, 2,id,nv,nvd);
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

int e00_polygons(TDAContext *ctx)
{
    register int i,k,l;
    int r,nrec,n,c,na,nn,m,u,us,id,id1,id2,keep,nw,nw1,nv,nvd;

    nrec = 0;
    printf1(ctx, "Number of polygons: %d\n",ctx->NPAL);
    printf1(ctx, "Number of arcs: %d\n",ctx->NPALA);

    nv = -1;
    nvd = 0;
    if (ctx->PMAttr > 0) {                     /* read variables from PAT section */
        if (ctx->PATFPtr >= 0) {
            printf1(ctx, "\nReading the PAT section.\n");
            nv = e00_var_read(ctx, 1,ctx->NPAL);
            if (nv > 0)
                e00_var_prn(ctx, nv);
            else {
                ctx->PATFPtr = -1;
                nv = -1;
            }
        }
        else {
            printf1(ctx, "Cannot find a PAT section.\n");
            nv = -1;
        }
    }
    if (ctx->PMAttr == 1) {
        nvd = 4;
        if (nv <= nvd)
            nv = -1;
    }
    printf1(ctx, "\nCreating the output file.\n");

    /* save polygon coordinates in AcX[i] and AcY[i], i = 0,...,NARC - 1.
       coverage number in AcI[j], number of points in line in AcN[j],
       j = 0,...,NARCP - 1, and pointer in AcK. */

    if (alloc_acx(ctx, ctx->NARCP + 1))
        return(-1);
    if (alloc_acy(ctx, ctx->NARCP + 1))
        return(-1);
    if (alloc_aci(ctx, ctx->NARC + 1))
        return(-1);
    if (alloc_acn(ctx, ctx->NARC + 1))
        return(-1);
    if (alloc_ack(ctx, ctx->NARC + 1))
        return(-1);

    if ((r = e00_arc(ctx, 1,ctx->E00DFlag,&n,&c,ctx->AcX,ctx->AcY,ctx->AcI,ctx->AcN,ctx->AcK)) < 0) {
        printf1(ctx, "Error in reading the ARC section.\n");
        return(-1);
    }

    /* save arc numbers from PAL section in AcJ */

    na = ctx->NPAL + ctx->NPALA;
    if (alloc_acj(ctx, na + 1))
        return(-1);

    if ((r = e00_pal(ctx, 1,ctx->E00DFlag,&n,&c,ctx->AcJ)) < 0) {
        e00_err(ctx, r);
        return(-1);
    }

    /* check whether the entries in AcI are in ascending order */

    for (i = 1; i < ctx->NARC; ++i) {
        if (ctx->AcI[i] <= ctx->AcI[i - 1]) {
            printf1(ctx, "Error: coverage numbers not in ascending order.\n");
            return(-1);
        }
    }

    /* write polygons, note that e00-polygons may consist of several parts */

    if (alloc_acm(ctx, ctx->NARC + 1))
        return(-1);

    nw1 = nw = id = id1 = 0;
    if (ctx->PMNC == 0)
        keep = 0;
    else
        keep = 1;
    i = 0;
    while (i < na) {
        id1++;
        id2 = 0;
        nn = 0;
        n = ctx->AcJ[i++];
        m = 0;
        for (k = 0; k < n; ++k) {
            u = ctx->AcJ[i + k];
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
                        r = e00_poly_prn(ctx, id,id1,id2,nn,m,ctx->AcM,ctx->AcK,ctx->AcN,nv,nvd);
                        if (r < 0)
                            return(-1);
                        nrec += r;
                    }
                    nn = m = 0;
                }
            }
            else {
                l = e00_find(ctx, u,ctx->NARC,ctx->AcI);
                if (l < 0)
                    nw++;
                if (l >= 0) {
                    nn += ctx->AcN[l];
                    l++;
                    if (us)
                        l = -l;
                    ctx->AcM[m++] = l;
                }
            }
        }
        if (m > 0 && keep) {
            if (nn < 3)
                nw1++;
            else {
                id++;
                id2++;
                r = e00_poly_prn(ctx, id,id1,id2,nn,m,ctx->AcM,ctx->AcK,ctx->AcN,nv,nvd);
                if (r < 0)
                    return(-1);
                nrec += r;
            }
            nn = m = 0;
        }
        i += n;
        keep = 1;
    }
    printf1(ctx, "Number of polygons in output file: %d\n",id);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);
    if (nw > 0)
        printf1(ctx, "Warning: could not find %d arcs from PAL section.\n",nw);
    if (nw1 > 0)
        printf1(ctx, "Warning: skipped %d polygon(s) with less than 3 points.\n",nw1);

    if (ctx->PMTDAFDef)                  /* write TDA description file */
        e00_dtda(ctx, 3,id,nv,nvd);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  e00_find(u,n,nlist)                                                     */
/*                                                                          */
/*  Return index of u in list nlist[i], i = 0,...,n-1.                      */
/*  or -1 if not found.                                                     */

int e00_find(TDAContext *ctx, int u,int n,int *nlist)
{
    (void)ctx;        /* unused: the signature is shared */
    register int l,r,k;

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

int e00_poly_prn(TDAContext *ctx, int id,int id1,int id2,int nn,int n,int *idx,int *ptr, int *np,int nv,int nvd)
{
    register int i,j,k,l;
    int m,s,nrec;
    char *p;

    fprintf(ctx->PMF1d,"%8d 3 %7d %8d %6d ",id,nn,id1,id2);
    if (nv > 0) {
        if (e00_var_buf(ctx, nv,id1 - 1))
            return(-1);
        p = ctx->E00Buf;
        for (i = 0; i < nvd; ++i)
            p += ATT[i].Len;
        fprintf(ctx->PMF1d,"%s",p);
    }
    fprintf(ctx->PMF1d,"\n");
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
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,(double)ctx->AcX[k + l]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,(double)ctx->AcY[k + l]);
                fprintf(ctx->PMF1d,"\n");
            }
        }
        else {
            for (l = m - 1; l >= 0; --l) {
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,(double)ctx->AcX[k + l]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,(double)ctx->AcY[k + l]);
                fprintf(ctx->PMF1d,"\n");
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

int e00_var_read(TDAContext *ctx, int opt,int ne)
{
    register int i,j;
    int l,n,nv,nva,nl,nrec,as,sp,nf,np,nt,nvl,vflag;
    char *p,*q;

    if (opt == 0) {
        if (fseek(ctx->E00Fd,ctx->AATFPtr,0) || !fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)) {
            printf1(ctx, "Cannot read AAT section. Will be ignored.\n");
            return(-1);
        }
    }
    else {
        if (fseek(ctx->E00Fd,ctx->PATFPtr,0) || !fgets(ctx->E00Buf,E00BufL,ctx->E00Fd)) {
            printf1(ctx, "Cannot read PAT section. Will be ignored.\n");
            return(-1);
        }
    }

    /* nv  = number of valid attributes */
    /* nva = total number of attributes */
    /* nl  = length */
    /* nrec = number of records */

    p = ctx->E00Buf + 34;
    if (sscanf(p,"%4d%4d%4d%10d",&nv,&nva,&nl,&nrec) != 4 || nv < 1 ||
        nva < 1 || nl < 1 || nrec < 1) {
        printf1(ctx, "Cannot read PAT section. Will be ignored.\n");
        return(-1);
    }
    if (nrec != ne) {
        printf1(ctx, "Inconsistent information in the AAT or PAT section.\n");
        printf1(ctx, "PAT section will be ignored.\n");
        return(-1);
    }
    printf1(ctx, "Number of valid entries: %d\n",nv);

    if (!(ATT = (struct E00ATT *)calloc((size_t)(nva),sizeof(struct E00ATT)))) {
        p_err(ctx, -2,1);
        return(-1);
    }
    ctx->ATTN = nva;
    memrq(ctx, ctx->ATTN,sizeof(struct E00ATT));

    /* save variable descriptions in E00ATT */

    for (i = 0; i < nva; ++i) {

        ATT[i].Valid = 0;

        if (!fgets(ctx->E00Buf,E00BufL,ctx->E00Fd))
            return(-1);
        p = ctx->E00Buf;
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
        if (sscanf(ctx->E00Buf + 16,"%3d",&as) != 1)
            return(-1);
        if (sscanf(ctx->E00Buf + 21,"%4d",&sp) != 1)
            return(-1);
        if (sscanf(ctx->E00Buf + 28,"%4d",&nf) != 1)
            return(-1);
        if (nf < 0)
            nf = 0;
        ATT[i].Fmt1 = nf;

        if (sscanf(ctx->E00Buf + 32,"%2d",&np) != 1)
            return(-1);
        if (np < 0)
            np = 0;
        ATT[i].Fmt2 = np;

        if (sscanf(ctx->E00Buf + 34,"%3d",&nt) != 1)
            return(-1);
        ATT[i].Typ = nt;

        if (sscanf(ctx->E00Buf + 66,"%3d",&nvl) != 1)
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
        printf1(ctx, "Cannot read identification numbers.\n");
        return(-1);
    }
    if (alloc_acs(ctx, ne + 1))
        return(-1);

    ctx->ATTLen = 0;
    for (i = 0; i < nva; ++i) {
        if (ATT[i].Valid == 1)
            ctx->ATTLen += ATT[i].Len;
    }
    n = ctx->ATTLen / 80;
    n = (n + 1) * 80 + 2;
    if (n > E00BufL) {
        printf1(ctx, "Exceeded maximal buffer length.\n");
        return(-1);
    }
    for (i = 0; i < nrec; ++i) {
        ctx->AcS[i] = (int)ftell(ctx->E00Fd);
        p = ctx->E00Buf;
        l = ctx->ATTLen;
        while (l > 0) {
            if (!fgets(p,l + 3,ctx->E00Fd))
                return(-1);

            n = imin(ctx, l,80);
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
        p = ctx->E00Buf;
        /****************
        tda_out("%s\n",p);
        ****************/
        p = ctx->E00Buf + ATT[0].Len + ATT[1].Len;
        *(p + ATT[2].Len) = '\0';
        q = skip_b(ctx, p);
        if (sscanf(q,"%d",&n) != 1) {
            printf1(ctx, "Cannot read identification numbers.\n");
            return(-1);
        }
    }
    return(nv);
}

/* ------------------------------------------------------------------------ */
/*  e00_var_prn(nva)    Print variable definitions from E00ATT to           */
/*                      standard output.                                    */

void e00_var_prn(TDAContext *ctx, int nva)
{
    register int i;
    int len,typ,n,m;

    len = 8;
    for (i = 0; i < nva; ++i) {
        if (ATT[i].Valid == 0)
            continue;
        len = (int)(imax(ctx, len,(int)(strlen(ATT[i].Name))));
    }
    printf1(ctx, "\nIdx  Variable   ");
    prnchar(ctx, ' ',len - 8,0);
    printf1(ctx, "Type       Format\n");
    prnchar(ctx, '-',len + 25,1);

    for (i = 0; i < nva; ++i) {
        if (ATT[i].Valid == 0)
            continue;

        printf1(ctx, "%3d  %s   ",i + 1,ATT[i].Name);
        prnchar(ctx, ' ',len - (int)strlen(ATT[i].Name),0);
        typ = ATT[i].Typ;
        if (typ == 10 || typ == 20)
            printf1(ctx, "string      %d\n",ATT[i].Len);
        else {
            if (typ == 60 || typ == 40) {       /* floating point */
                m = ATT[i].Fmt2;
                n = ATT[i].Fmt1 + m + 1;
            }
            else {                              /* integer */
                m = 0;
                n = ATT[i].Fmt1;
            }
            printf1(ctx, "numerical   %d.%d\n",n,m);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  e00_var_buf(nva,irec)                                                   */
/*                                                                          */
/*  Create in E00Buf the variables for                                      */
/*  Return 0 if OK, -e if error.                                            */

int e00_var_buf(TDAContext *ctx, int nva,int irec)
{
    (void)nva;        /* unused: the signature is shared */
    register int j;
    int n,l;
    char *p;

    if (fseek(ctx->E00Fd,(long)ctx->AcS[irec],0))
        goto E00VWErr;

    n = ctx->ATTLen / 80;
    n = (n + 1) * 80 + 2;

    p = ctx->E00Buf;
    l = ctx->ATTLen;
    while (l > 0) {
        if (!fgets(p,l + 3,ctx->E00Fd))
            goto E00VWErr;

        n = imin(ctx, l,80);
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
    printf1(ctx, "Cannot read variables in AAT or PAT section.\n");
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  e00_dtda(opt,noc,nv,nvd)    Create TDA description file.                */
/*                                                                          */
/*  opt = 1 : point coverage                                                */
/*  opt = 2 : line coverage                                                 */
/*  opt = 3 : polygon coverage                                              */

void e00_dtda(TDAContext *ctx, int opt,int noc,int nv,int nvd)
{
    register int i,j = 0,k,l = 0;
    int n = 0,m,typ,len;
    char *name;

    if (ctx->PMTDAFDef == 0)
        return;

    fprintf(ctx->PMTDAFd,"sdnvar(\n");
    fprintf(ctx->PMTDAFd,"  dfile = %s,\n",ctx->PMF1dName);
    fprintf(ctx->PMTDAFd,"  noc = %d,\n",noc);

    if (opt == 1) {
        fprintf(ctx->PMTDAFd,"  ffmt=c1(1-8),c2(10),c3(12),");
        n = 3; j = 4; l = 14;
    }
    else if (opt == 2) {
        fprintf(ctx->PMTDAFd,"  ffmt=c1(1-8),c2(10),c3(12-18),");
        n = 3; j = 4; l = 20;
    }
    else if (opt == 3) {
        fprintf(ctx->PMTDAFd,"  ffmt=c1(1-8),c2(10),c3(12-18),c4(20-27),c5(29-34),");
        n = 5; j = 6; l = 36;
    }
    for (i = 0; i < nv; ++i) {
        if (i < nvd)
            continue;

        k = l + ATT[i].Len;
        if (ATT[i].Typ != 10 && ATT[i].Typ != 20) {
            if (++n > 5) {
                fprintf(ctx->PMTDAFd,"\n       ");
                n = 1;
            }
            fprintf(ctx->PMTDAFd,"c%d(%d-%d),",j,l,k - 1);
            j++;
        }
        l = k;
    }
    fprintf(ctx->PMTDAFd,"\n  SDID  <5>[ 8.0] = c1,\n");
    fprintf(ctx->PMTDAFd,"  SDTyp <1>[ 1.0] = c2,\n");
    if (opt == 1)
        fprintf(ctx->PMTDAFd,"  SDN   <1>[ 1.0] = c3,\n");
    else
        fprintf(ctx->PMTDAFd,"  SDN   <5>[ 7.0] = c3,\n");
    fprintf(ctx->PMTDAFd,"  SDPtr <5>[10.0] = rd,\n");

    if (opt == 1) {
        j = 4;
        l = 14;
    }
    else if (opt == 2) {
        j = 4;
        l = 20;
    }
    else if (opt == 3) {
        fprintf(ctx->PMTDAFd,"  SDID1 <5>[ 8.0] = c4,\n");
        fprintf(ctx->PMTDAFd,"  SDID2 <5>[ 6.0] = c5,\n");
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
            fprintf(ctx->PMTDAFd,"  %s = str(%d,%d),\n",name,l,k - 1);
        else {
            if (typ == 60 || typ == 40) {       /* floating point */
                m = ATT[i].Fmt2;
                n = ATT[i].Fmt1 + m + 1;
            }
            else {                              /* integer */
                m = 0;
                n = ATT[i].Fmt1;
            }
            fprintf(ctx->PMTDAFd,"  %s [%d.%d] = c%d,\n",name,n,m,j);
            j++;
        }
        l = k;
    }
    fprintf(ctx->PMTDAFd,");\n");
    printf1(ctx, "TDA description written to: %s\n",ctx->PMTDAFName);
}
