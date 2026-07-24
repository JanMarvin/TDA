#include "tda_rhooks.h"
/****************************************************************************/
/*  t_gen                                                                   */
/*                                                                          */
/*  TDA. Program for Transition Data Analysis, written by Goetz Rohwer.     */
/*  Copyright (C) 1989,1991-94 Goetz Rohwer. All rights reserved.           */
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
#include "t_var.h"
#include "t_gf.h"
#include <stdarg.h>
#include "tda_context.h"

/*  functions in t_gen.c */

void gerr_exit(TDAContext *ctx, int n);
void memrq(TDAContext *ctx, int n, int m);
void prn_mem(TDAContext *ctx);      
void newline(TDAContext *ctx);           
char *check_comment(TDAContext *ctx, char *p);
int check_drec(TDAContext *ctx, char *buf);
void prn_message(TDAContext *ctx, int rec,int opt,int wflag);
void get_sline(TDAContext *ctx, char *p); 
void prnchar(TDAContext *ctx, unsigned char c, int n, int mode);
void fprnchar(TDAContext *ctx, FILE *fd, unsigned char c, int n, int mode);
int mpr(TDAContext *ctx, int m,int n,double *a,int ivflg,double *b,char *fmt,char *fn,char *s,int aflag);
int mprf(TDAContext *ctx, int m,int n,float *a,char *fmt,FILE *fd);
int mprd(TDAContext *ctx, int m,int n,double *a,char *fmt,FILE *fd);
char *skip_b(TDAContext *ctx, char *p);
char *skip_c(TDAContext *ctx, char *p);
char *skip_cb(TDAContext *ctx, char *p);
char *skip_int(TDAContext *ctx, char *p);
char *skip_dbl(TDAContext *ctx, char *p);
char *skip_expr(TDAContext *ctx, char *p);
char *skip_blev(TDAContext *ctx, char *p);
char *skip_nc(TDAContext *ctx, char *p);
char *skip_com(TDAContext *ctx, char *p);
char *skip_xa(TDAContext *ctx, char *p);
int get_fname(TDAContext *ctx, char *p,char *q);
void printf1(TDAContext *ctx, const char *fmt, ...) TDA_PRINTF(2, 3);
void printf2(TDAContext *ctx, const char *fmt, ...) TDA_PRINTF(2, 3);
void printfe(TDAContext *ctx, const char *fmt, ...) TDA_PRINTF(2, 3);
void fflushe(TDAContext *ctx);
#if TIME_ON
void prn_time(TDAContext *ctx, FILE *fd);          
void prn1_time(TDAContext *ctx, char *buf); 
#endif

/*--------------------------------------------------------------------------*/
/*  Global variables                                                        */

                            /*  1 drop standard output                      */
                            /*  2 drop standard error output                */
                            /*  3 drop both                                 */
                            /* -1 include matrix and rep messages           */  





/****************************************************************************/
/*  gerr_exit(n)    Error message and exit.                                 */

void gerr_exit(TDAContext *ctx, int n)
{
    printfe(ctx, "ERROR %d. CAN'T CONTINUE.\n",n);
    printfe(ctx, "PLEASE REPORT THIS ERROR TO THE PROGRAM'S AUTHOR.\n");
    exit(0);
}

/*--------------------------------------------------------------------------*/
/*  memrq    Bookkeeping of memory usage.                                   */

void memrq(TDAContext *ctx, int n, int m)
{
    ctx->MemReq += n * m;
    if (ctx->MxMReq < ctx->MemReq)
        ctx->MxMReq = ctx->MemReq;
}

/* ------------------------------------------------------------------------ */
/*  prn_mem    print current memory.                                        */

void prn_mem(TDAContext *ctx)       
{
    printf1(ctx, "Current memory: %d bytes.\n",ctx->MemReq);
}

/* ------------------------------------------------------------------------ */
/*  prn_time(fd)  print current date and time to a file.                    */

#if TIME_ON

void prn_time(TDAContext *ctx, FILE *fd)       
{
    time_t t;
    if (time(&t) != -1) {
        if (fd == TDA_CONSOLE)
            printf1(ctx, "%s", ctime(&t));
        else
            fprintf(fd, "%s", ctime(&t));
    }
}
#endif

/* ------------------------------------------------------------------------ */
/*  prn1_time(buf)  print current date and time to buf.                     */

#if TIME_ON

void prn1_time(TDAContext *ctx, char *buf)  
{
    (void)ctx;        /* unused: the signature is shared */
    time_t t;

    *buf = '\0';
    if (time(&t) != -1)
        strcpy(buf,ctime(&t));
}
#endif

/* ------------------------------------------------------------------------ */
/*  newline()   print newline                                               */

void newline(TDAContext *ctx)            
{
    printf1(ctx, "\n");
}

/* ------------------------------------------------------------------------ */
/*  check_comment(p)    if string p is a comment line return NULL,          */
/*                      otherwise return pointer to first non-blank char.   */

char *check_comment(TDAContext *ctx, char *p)
{
    (void)ctx;        /* unused: the signature is shared */
    while (*p && *p == ' ')
        p++;

    if (*p && *p != '#' && *p != '*' && *p != ';' && *p != '\n'
                                                  && *p != LF && *p != CR)
        return(p);
    return(NULL);
}

/*--------------------------------------------------------------------------*/
/*  check_drec(buf)     return 1 if buf contains a data file record,        */
/*                      otherwise (if comment) return 0.                    */

int check_drec(TDAContext *ctx, char *buf)
{
    (void)ctx;        /* unused: the signature is shared */
    register char *p;

    p = buf;
    while (*p == ' ')
        p++;

    if (!*p || *p == '#' || *p == LF || *p == CR)
        return(0);
    return(1);
}

/*--------------------------------------------------------------------------*/
/*  prn_message(rec,opt,wflag)  print read message for current record rec.  */
/*                              if opt == 0 print according to PrnDfCnt,    */
/*                              otherwise print unconditionally.            */
/*                              If wflag != 0 print written instead of read */
          
void prn_message(TDAContext *ctx, int rec,int opt,int wflag)
{
    if (ctx->PrnDfCnt > 0 && ctx->SILENTFlg < 2) {
        if (opt || (rec / ctx->PrnDfCnt) * ctx->PrnDfCnt == rec) {
            if (wflag)
                printfe(ctx, "Records written: %7d     %c",rec,CR);
            else
                printfe(ctx, "Read records: %7d     %c",rec,CR);
            if (opt)
                printfe(ctx, "\n");
            fflushe(ctx);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  get_sline(p)   remove blanks from string p until end,                   */
/*                 except blanks in " " or ' ' or after help.               */

void get_sline(TDAContext *ctx, char *p)
{
    (void)ctx;        /* unused: the signature is shared */
    register char cflag,*q,*s;
    int hflag,cnt;

    s = q = p;   
    cflag = '\0';
    cnt = hflag = 0;
    while (*p && *p != '\n') {
/**     if (*p == '"' || *p == '`' || (*q == '\'' && *(q-1) != '>')) {  **/
        if (*p == '"' || *p == '`') {
            if (cflag == *p)
                cflag = '\0';
            else if (!cflag)
                cflag = *p;
            else {
                *q++ = *p;
                cnt++;
            }
        }
        else if ((*p != ' ' && *p != '\t') || cflag || hflag) {    
            *q++ = *p;
            cnt++;
        }
        if (cnt == 4 && !strncmp(s,"help",4))
            hflag = 1;
        p++;
    }
    *q = '\0';
}          

/* ------------------------------------------------------------------------ */
/*  prnchar(c,n,mode)                                                       */
/*      Print the character c n times to stdout. If mode != 0 append a      */
/*      newline.                                                            */

void prnchar(TDAContext *ctx, unsigned char c, int n, int mode)
{
    register int i;

    for (i = 0; i < n; ++i)  
        printf1(ctx, "%c",c);
    if (mode)  
        printf1(ctx, "\n");
}

/* ------------------------------------------------------------------------ */
/*  fprnchar(fd,c,n,mode)                                                   */
/*      Prints the character c n times into the file fd. If mode != 0       */
/*      append a newline.                                                   */

void fprnchar(TDAContext *ctx, FILE *fd, unsigned char c, int n, int mode)
{
    register int i;

    if (fd == TDA_CONSOLE)
        prnchar(ctx, c,n,mode);
    else {
        for (i = 0; i < n; ++i)
            fprintf(fd,"%c",c);
        if (mode)
            fprintf(fd,"\n");
    }
}

/*--------------------------------------------------------------------------*/
/*  mpr(m,n,a,ivflg,b,fmt,fn,s,aflag)                                       */
/*                                                                          */
/*  Print m,n matrix a with fmt to file fn. If ivflg > 0 also print b.      */
/*  If ivflg = 2 use square brackets.                                       */
/*  if s != NULL print string s first.                                      */
/*  if aflag != 0 append to file.                                           */
/*  return 0 if successful, -1 if error.                                    */

int mpr(TDAContext *ctx, int m,int n,double *a,int ivflg,double *b,char *fmt,char *fn,char *s,int aflag)
{
    FILE *fd = NULL;
    register int i,j,ii;

    if (fn != NULL) {
        if (aflag == 0) {
            if (!(fd = fopen(fn,OPEN_WR))) {
                printf1(ctx, "Error: can't open output file: %s\n",fn);
                return(-1);       
            }
        }
        else {
            if (!(fd = fopen(fn,OPEN_AP))) {
                printf1(ctx, "Error: can't append to output file: %s\n",fn);
                return(-1);       
            }
        }
    }
    if (s != NULL) {
        if (fn == NULL)
            printf1(ctx, "%s\n",s);
        else
            fprintf(fd,"%s\n",s);
    }
    for (i = 1; i <= m; ++i) {
        ii = (i - 1) * n;
        for (j = 1; j <= n; ++j) {
            if (fn == NULL) {
                if (ivflg == 2)
                    printf1(ctx, "[");
                rt_printf1_d(ctx, fmt,a[ii + j]);
                if (ivflg) {
                    rt_printf1_d(ctx, fmt,b[ii + j]);
                    if (ivflg == 2)
                        printf1(ctx, "]");
                }
            }
            else {
                if (ivflg == 2)
                    fprintf(fd,"[");
                rt_fprintf_d(ctx, fd,fmt,a[ii + j]);
                if (ivflg) {
                    rt_fprintf_d(ctx, fd,fmt,b[ii + j]);
                    if (ivflg == 2)
                        fprintf(fd,"]");
                }
            }
#ifdef TDA_R_PACKAGE
            /* both the console and the file branch: mpr() is the
               matrix printer behind mpr()= in a command file, which is
               how brr gets its S matrix */
            tda_export_cell(ctx, "mpr.matrix", a[ii + j]);
            if (ivflg)
                tda_export_cell(ctx, "mpr.matrix", b[ii + j]);
#endif
        }
        if (fn == NULL)  
            printf1(ctx, "\n");
        else
            fprintf(fd,"\n");
#ifdef TDA_R_PACKAGE
        tda_export_endrow(ctx, "mpr.matrix");
#endif

    }
#ifdef TDA_R_PACKAGE
    /* one flush per call, so several matrices become mpr.matrix,
       mpr.matrix.2, ... in call order */
    tda_export_flush(ctx, "mpr.matrix");
#endif
    if (fn != NULL)
        fclose(fd);

    return(0);
}

/*--------------------------------------------------------------------------*/
/*  mprf(m,n,a,fmt,fd)      print m,n matrix a with fmt to file fd          */
/*                          return 0 if successful, -1 if error.            */

int mprf(TDAContext *ctx, int m,int n,float *a,char *fmt,FILE *fd)
{
    register int i,j,ii;

    for (i = 1; i <= m; ++i) {
        ii = (i - 1) * n;
        for (j = 1; j <= n; ++j)  
            rt_fprintf_d(ctx, fd,fmt,(double)(a[ii + j]));
        fprintf(fd,"\n");
    }
    return(0);
}
   
/*--------------------------------------------------------------------------*/
/*  mprd(m,n,a,fmt,fd)      print m,n matrix a with fmt to file fd          */
/*                          return 0 if successful, -1 if error.            */

int mprd(TDAContext *ctx, int m,int n,double *a,char *fmt,FILE *fd)
{
    register int i,j,ii;

    for (i = 1; i <= m; ++i) {
        ii = (i - 1) * n;
        for (j = 1; j <= n; ++j)  
            rt_fprintf_d(ctx, fd,fmt,a[ii + j]);
        fprintf(fd,"\n");
    }
    return(0);
}
   
/* ------------------------------------------------------------------------ */
/*  skip_b(p) Skip blanks and tabs at pointer p and return pointer to next  */
/*            non-blank character.                                          */

char *skip_b(TDAContext *ctx, char *p)
{
    (void)ctx;        /* unused: the signature is shared */
    while (*p && (*p == ' ' || *p == '\t'))
        p++;
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  skip_c(p) Skip characters at pointer p until next blank or tab          */
/*            character.                                                    */

char *skip_c(TDAContext *ctx, char *p)
{
    (void)ctx;        /* unused: the signature is shared */
    while (*p && *p != ' ' && *p != '\t')
        p++;
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  skip_cb(p)  Skip characters at pointer p, first until next blank        */
/*              or tab character, then until next non-blank or non-tab.     */

char *skip_cb(TDAContext *ctx, char *p)
{
    (void)ctx;        /* unused: the signature is shared */
    while (*p && *p != ' ' && *p != '\t')
        p++;
    while (*p && (*p == ' ' || *p == '\t'))
        p++;
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  skip_int(p)                                                             */
/*      It is assumed that p is a pointer to an integer. The function       */
/*      returns a pointer to the next character after this integer.         */

char *skip_int(TDAContext *ctx, char *p)
{
    (void)ctx;        /* unused: the signature is shared */
    if (*p == '-' || *p == '+')
        p++;
   
    while (*p && isdigit((int)*p))
        p++;
   
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  skip_dbl(p)                                                             */
/*      It is assumed that p is a pointer to an double or floating point    */
/*      value. The function returns a pointer to the next character after   */
/*      this value.                                                         */

char *skip_dbl(TDAContext *ctx, char *p)
{
    (void)ctx;        /* unused: the signature is shared */
   
    while (*p && (isdigit((int)*p) || *p == '.' || *p == '-' || *p == '+' ||
                                                   *p == 'e' || *p == 'E'))
        p++;
   
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  skip_expr(p)                                                            */
/*      It is assumed that p is a pointer to an expression. The function    */
/*      returns a pointer to the first , or ) after the expression.         */
/*      end of the expression.                                              */

char *skip_expr(TDAContext *ctx, char *p)
{
    (void)ctx;        /* unused: the signature is shared */
    int n,m,m1;

    m1 = m = n = 0;

    while (*p) {
        if (*p == ',' && m1 == 0 && m == 0 && n == 0)
            break;
        if (*p == '(')
            n++;
        else if (*p == '[')
            m++;
        else if (*p == '<')
            m1++;
        else if (*p == ')') {
            if (m1 == 0 && m == 0 && n == 0)
                break;
            n--;
        }
        else if (*p == ']') {
            if (m1 == 0 && m == 0 && n == 0)
                break;
            m--;
        }
        else if (*p == '>') {
            if (m1 == 0 && m == 0 && n == 0)
                break;
            m1--;
        }
        p++;
    }
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  skip_blev(p)                                                            */
/*      It is assumed that p is a pointer to a '('  '[' or '{' or `<'       */
/*      The function returns a pointer to the first character after a       */  
/*      matching ')' or ']' or '}' or `>'.                                  */              

char *skip_blev(TDAContext *ctx, char *p)
{
    (void)ctx;        /* unused: the signature is shared */
    int m,n,mm,m1;
    char c;

    if (*p == '(') 
        c = ')';
    else if (*p == '[')
        c = ']';
    else if (*p == '{')
        c = '}';
    else if (*p == '<')
        c = '>';
    else
        return(p);

    p++;
    m1 = mm = m = n = 0;
    while (*p) {
        if (*p == c && m == 0 && n == 0 && mm == 0 && m1 == 0) {
            p++;
            break;
        }
        if (*p == '(')
            n++;
        else if (*p == '[')
            m++;
        else if (*p == '{')
            mm++;
        else if (*p == '<')
            m1++;
        else if (*p == ')') {
            if (*p == c && m == 0 && n == 0 && mm == 0 && m1 == 0) {
                p++;
                break;
            }
            n--;
        }
        else if (*p == ']') {
            if (*p == c && m == 0 && n == 0 && mm == 0 && m1 == 0) {
                p++;
                break;
            }
            m--;
        }
        else if (*p == '}') {
            if (*p == c && m == 0 && n == 0 && mm == 0 && m1 == 0) {
                p++;
                break;
            }
            mm--;
        }
        else if (*p == '>') {
            if (*p == c && m == 0 && n == 0 && mm == 0 && m1 == 0) {
                p++;
                break;
            }
            m1--;
        }
        p++;
    }
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  skip_nc(p)  skip until next comma or ')',                               */
/*              recognize (..) [..] {} <> pairs                             */

char *skip_nc(TDAContext *ctx, char *p)
{
    while (*p) {
        if (*p == ',' || *p == ')')
            return(p);
        else if (*p == '(' || *p == '[' || *p == '{' || *p == '<')
            p = skip_blev(ctx, p);
        else
            p++;
    }
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  skip_com(p)     skip until comma or ')'                                 */

char *skip_com(TDAContext *ctx, char *p)
{
    (void)ctx;        /* unused: the signature is shared */
    while (*p && *p != ',' && *p != ')')
        p++;
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  skip_xa(p)  skip xa(...)=varlist string                                 */

char *skip_xa(TDAContext *ctx, char *p)
{
    register int l;

    p += 2;
    p = skip_blev(ctx, p);
    if (*p != '=')
        return(p);

    if (*(p + 1) == '1' && *(p + 2) == ',')
        return(p + 2);
           
    while (*(p + 1)) {
        l = get_vnlen(ctx, p + 1);
        if (l <= 0)
            break;
        p += l + 1;
    }
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  get_fname(p,q)          get file name beginning at p and save in q.     */
/*                          if successful return 1, else 0.                 */

int get_fname(TDAContext *ctx, char *p,char *q)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;

    if (!*p)  
        return(0);
       
    for (i = 0; i < FNMaxLen; ++i) {
        if (!*p)
            break;
        *q++ = *p++;
    }
    *q = '\0';
    return(1);
}

#ifdef TDA_R_PACKAGE
/* ------------------------------------------------------------------------ */
/*  Generic print tap (package-only).                                       */
/*                                                                          */
/*  Rather than a producer beside every one of TDA's ~4300 numeric print    */
/*  sites, printf1() -- through which the console output all passes -- is   */
/*  tapped once here.  The format string is walked, every numeric           */
/*  conversion's own argument is taken from a va_copy of the same           */
/*  argument list the text path consumes, and the value is staged as        */
/*  (line, value).  Line numbers count newlines actually emitted, so a      */
/*  reader can map any output line back to the doubles behind it at full    */
/*  precision.                                                              */
/*                                                                          */
/*  The text path is untouched: the tap only reads a copy.  Arguments must  */
/*  be consumed in exactly printf's own order, including "*" width and      */
/*  precision, or the walk desynchronises -- which is why every conversion  */
/*  is handled, not only the numeric ones.                                  */

void tda_print_tap(TDAContext *ctx, const char *fmt, va_list ap);

void tda_print_tap(TDAContext *ctx, const char *fmt, va_list ap)
{
    const char *p;
    const char *lit_end;
    double v;
    int lmod, has;

    if (ctx == NULL || fmt == NULL)
        return;

    /* tda_out_line_start() is the count of lines already written by the
       output layer, so the line being written now is the next one:
       1-based, matching the index of res$output on the R side.  The tap
       cannot count fmt's own newlines instead -- a "%s" argument can
       carry one, which is how the command-file echo desynchronised it. */
    ctx->RExportLine = (int)tda_out_line_start() + 1;

    for (p = fmt; *p; ++p) {
        if (*p == '\n') {
            ctx->RExportLine++;
            continue;
        }
        if (*p != '%')
            continue;
        lit_end = p;              /* literal text ends where this spec starts */
        ++p;
        if (*p == '%')
            continue;
        while (*p && strchr("-+ #0'", *p))    /* flags */
            ++p;
        if (*p == '*') {                      /* width from an argument */
            (void)va_arg(ap, int);
            ++p;
        }
        else
            while (*p >= '0' && *p <= '9')
                ++p;
        if (*p == '.') {                      /* precision */
            ++p;
            if (*p == '*') {
                (void)va_arg(ap, int);
                ++p;
            }
            else
                while (*p >= '0' && *p <= '9')
                    ++p;
        }
        lmod = 0;                             /* length modifier */
        while (*p && strchr("hlLqjzt", *p)) {
            if (*p == 'l')
                lmod++;
            else if (*p == 'L')
                lmod = 3;
            else if (*p == 'h')
                lmod = -1;
            ++p;
        }
        if (!*p)
            break;
        has = 1;
        v = 0.0;
        switch (*p) {
            case 'd': case 'i':
                if (lmod >= 2)      v = (double)va_arg(ap, long long);
                else if (lmod == 1) v = (double)va_arg(ap, long);
                else                v = (double)va_arg(ap, int);
                break;
            case 'u': case 'o': case 'x': case 'X':
                if (lmod >= 2)      v = (double)va_arg(ap, unsigned long long);
                else if (lmod == 1) v = (double)va_arg(ap, unsigned long);
                else                v = (double)va_arg(ap, unsigned int);
                break;
            case 'e': case 'E': case 'f': case 'F':
            case 'g': case 'G': case 'a': case 'A':
                if (lmod == 3)      v = (double)va_arg(ap, long double);
                else                v = va_arg(ap, double);
                break;
            case 'c':
                (void)va_arg(ap, int);
                has = 0;
                break;
            case 's':
                (void)va_arg(ap, char *);
                has = 0;
                break;
            case 'p':
                (void)va_arg(ap, void *);
                has = 0;
                break;
            case 'n':
                (void)va_arg(ap, int *);
                has = 0;
                break;
            default:
                has = 0;
                break;
        }
        if (has) {
            double row[2];
            char lab[64];
            const char *b = fmt, *q2;
            size_t n2 = 0;

            row[0] = (double)ctx->RExportLine;
            row[1] = v;
            tda_export_row(ctx, "print.values", row, 2);
            /* The literal text of the format string in front of this
               conversion, so a value can be found by the LABEL TDA
               printed instead of by grepping the output for its line.
               That is what lets a reader take a scalar without the
               text at all.  A format whose label arrives through "%s"
               has no literal to take and stages an empty string. */
            for (q2 = fmt; q2 < lit_end; ++q2)
                if (*q2 == '\n')
                    b = q2 + 1;
            while (b < lit_end && (*b == ' ' || *b == '\t'))
                ++b;
            while (b < lit_end && n2 < sizeof(lab) - 1)
                lab[n2++] = *b++;
            while (n2 > 0 && (lab[n2 - 1] == ' ' || lab[n2 - 1] == '\t'))
                n2--;
            lab[n2] = '\0';
            tda_export_str_row(ctx, "print.labels", lab);
        }
    }
}
#endif

/* ------------------------------------------------------------------------ */
/*  printf1()       General printf function that recognizes SILENTFlg.      */

/*  Does this format string introduce an error message?  TDA has no error
    channel: an error is a printed line and the exit code stays 0 either
    way.  The wording is regular enough to recognise -- some 900 format
    strings begin with "Error" or "Syntax error" -- so a run can be
    counted in one place instead of at every site.

    What must NOT match is the coefficient table's own column header,
    which is printed as a bare "Error" or "Error " after its padding
    (t_pgen.c, t_gmin.c, t_nlreg.c, t_qrmod.c).  Hence the requirement
    that something follow the word.  Headers written as "  C/Error" keep
    their leading blanks and never matched in the first place. */

static int is_err_fmt(const char *f)
{
    while (*f == '\n' || *f == ' ' || *f == '\t')   /* "\nSyntax error." */
        ++f;
    if (!strncmp(f,"Syntax error",12) || !strncmp(f,"error",5) ||
        !strncmp(f,"Fatal",5) || !strncmp(f,"FATAL",5))
        return(1);
    if (strncmp(f,"Error",5) && strncmp(f,"ERROR",5))
        return(0);
    f += 5;
    if (*f == '\0' || (*f == ' ' && f[1] == '\0'))
        return(0);              /* the table header, "Error" or "Error " */
    return(1);
}

void printf1(TDAContext *ctx, const char *fmt, ...) 
{
    va_list ap;

    /* before the SILENTFlg return, so a silenced run still reports its
       errors -- the count is the only thing the caller can test */
    if (is_err_fmt(fmt))
        ++ctx->ErrCnt;

    if (ctx->SILENTFlg == 1 || ctx->SILENTFlg == 3)
        return;

    va_start(ap,fmt);
    tda_vout(fmt,ap);
    va_end(ap);
#ifdef TDA_R_PACKAGE
    va_start(ap,fmt);
    tda_print_tap(ctx, fmt,ap);
    va_end(ap);
#endif
}

/* printf1() with the text already formatted; the arguments after it are
   fmt's own and only feed the print tap. */
void printf1_rt(TDAContext *ctx, const char *fmt, const char *text, ...)
{
#ifdef TDA_R_PACKAGE
    va_list ap;
#endif

    if (is_err_fmt(fmt))
        ++ctx->ErrCnt;

    if (ctx->SILENTFlg == 1 || ctx->SILENTFlg == 3)
        return;

    tda_out("%s", text);
#ifdef TDA_R_PACKAGE
    va_start(ap,text);
    tda_print_tap(ctx, fmt,ap);
    va_end(ap);
#endif
}

/* ------------------------------------------------------------------------ */
/*  printf2()       General printf function that recognizes SILENTFlg.      */

void printf2(TDAContext *ctx, const char *fmt, ...) 
{
    va_list ap;

    if (ctx->SILENTFlg >= 0)
        return;

    va_start(ap,fmt);
    tda_vout(fmt,ap);
    va_end(ap);
}

/* ------------------------------------------------------------------------ */
/*  printfe()       fprintf to stderr.                                      */

#define FE_BufLen 132
char FE_Buf[FE_BufLen + 1];
int FE_BufPtr = 0;

void printfe(TDAContext *ctx, const char *fmt, ...) 
{
    va_list ap;

    /* the stderr twin of printf1(): gerr_exit() and the FATAL/DERIV
       reports go out this way, so they are counted here too */
    if (is_err_fmt(fmt))
        ++ctx->ErrCnt;

    va_start(ap,fmt);
    tda_vout_err(fmt,ap);
    va_end(ap);
}

/* ------------------------------------------------------------------------ */
/*  fflushe()     fflush stderr.                                            */

void fflushe(TDAContext *ctx)
{
    (void)ctx;        /* unused: the signature is shared */
    tda_err_flush();
}






void tda_reset_t_gen(void)
{
    memset(FE_Buf, 0, sizeof(FE_Buf));
    FE_BufPtr = 0;
}
