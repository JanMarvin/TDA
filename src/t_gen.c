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

/*  functions in t_gen.c */

void gerr_exit(int n);
void memrq(int n, int m);
void prn_mem(void);      
void newline(void);           
char *check_comment(char *p);
int check_drec(char *buf);
void prn_message(int rec,int opt,int wflag);
void get_sline(char *p); 
void prnchar(unsigned char c, int n, int mode);
void fprnchar(FILE *fd, unsigned char c, int n, int mode);
int mpr(int m,int n,double *a,int ivflg,double *b,char *fmt,char *fn,char *s,int aflag);
int mprf(int m,int n,float *a,char *fmt,FILE *fd);
int mprd(int m,int n,double *a,char *fmt,FILE *fd);
char *skip_b(char *p);
char *skip_c(char *p);
char *skip_cb(char *p);
char *skip_int(char *p);
char *skip_dbl(char *p);
char *skip_expr(char *p);
char *skip_blev(char *p);
char *skip_nc(char *p);
char *skip_com(char *p);
char *skip_xa(char *p);
int get_fname(char *p,char *q);
void printf1(const char *fmt, ...);
void printf2(const char *fmt, ...);
void printfe(const char *fmt, ...); 
void fflushe(void);
#if TIME_ON
void prn_time(FILE *fd);          
void prn1_time(char *buf); 
#endif

/*--------------------------------------------------------------------------*/
/*  Global variables                                                        */

int SILENTFlg = 0;          /*  0 default                                   */
                            /*  1 drop standard output                      */
                            /*  2 drop standard error output                */
                            /*  3 drop both                                 */
                            /* -1 include matrix and rep messages           */  

int MemReq = 0;             /* number of bytes of requested memory          */
int MxMReq = 0;             /* max memory request (bytes)                   */

int CmdBufL = 0;            /* actual command buffer length                 */
char *CmdBuf;               /* command buffer                               */

int INTMAX = 1;             /* largest integer                              */
double EPSI = 1.0;          /* machine's epsilon                            */
double EPSI1 = 1.0;         /* sqrt(EPSI)                                   */
double EPSI2 = 1.0;         /* 1000 * EPSI                                  */
double DBLMAX = 1.0;        /* largest number in machine                    */
double DBLMIN = 1.0;        /* smallest number in machine                   */
float MINFLOAT = 1.1754e-38;

int PrnDfCnt = 1000;        /* increment for read messages                  */

/****************************************************************************/
/*  gerr_exit(n)    Error message and exit.                                 */

void gerr_exit(int n)
{
    printfe("ERROR %d. CAN'T CONTINUE.\n",n);
    printfe("PLEASE REPORT THIS ERROR TO THE PROGRAM'S AUTHOR.\n");
    exit(0);
}

/*--------------------------------------------------------------------------*/
/*  memrq    Bookkeeping of memory usage.                                   */

void memrq(int n, int m)
{
    MemReq += n * m;
    if (MxMReq < MemReq)
        MxMReq = MemReq;
}

/* ------------------------------------------------------------------------ */
/*  prn_mem    print current memory.                                        */

void prn_mem(void)       
{
    printf1("Current memory: %d bytes.\n",MemReq);
}

/* ------------------------------------------------------------------------ */
/*  prn_time(fd)  print current date and time to a file.                    */

#if TIME_ON

void prn_time(FILE *fd)       
{
    time_t t;
    if (time(&t) != -1) {
        if (fd == stdout)
            printf1(ctime(&t));
        else
            fprintf(fd,ctime(&t));
    }
}
#endif

/* ------------------------------------------------------------------------ */
/*  prn1_time(buf)  print current date and time to buf.                     */

#if TIME_ON

void prn1_time(char *buf)  
{
    time_t t;

    *buf = '\0';
    if (time(&t) != -1)
        strcpy(buf,ctime(&t));
}
#endif

/* ------------------------------------------------------------------------ */
/*  newline()   print newline                                               */

void newline(void)            
{
    printf1("\n");
}

/* ------------------------------------------------------------------------ */
/*  check_comment(p)    if string p is a comment line return NULL,          */
/*                      otherwise return pointer to first non-blank char.   */

char *check_comment(char *p)
{
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

int check_drec(char *buf)
{
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
          
void prn_message(int rec,int opt,int wflag)
{
    if (PrnDfCnt > 0 && SILENTFlg < 2) {
        if (opt || (rec / PrnDfCnt) * PrnDfCnt == rec) {
            if (wflag)
                printfe("Records written: %7d     %c",rec,CR);
            else
                printfe("Read records: %7d     %c",rec,CR);
            if (opt)
                printfe("\n");
            fflushe();
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  get_sline(p)   remove blanks from string p until end,                   */
/*                 except blanks in " " or ' ' or after help.               */

void get_sline(char *p)
{
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

void prnchar(unsigned char c, int n, int mode)
{
    register int i;

    for (i = 0; i < n; ++i)  
        printf1("%c",c);
    if (mode)  
        printf1("\n");
}

/* ------------------------------------------------------------------------ */
/*  fprnchar(fd,c,n,mode)                                                   */
/*      Prints the character c n times into the file fd. If mode != 0       */
/*      append a newline.                                                   */

void fprnchar(FILE *fd, unsigned char c, int n, int mode)
{
    register int i;

    if (fd == stdout)
        prnchar(c,n,mode);
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

int mpr(int m,int n,double *a,int ivflg,double *b,char *fmt,char *fn,char *s,int aflag)
{
    FILE *fd;
    register int i,j,ii;

    if (fn != NULL) {
        if (aflag == 0) {
            if (!(fd = fopen(fn,OPEN_WR))) {
                printf1("Error: can't open output file: %s\n",fn);
                return(-1);       
            }
        }
        else {
            if (!(fd = fopen(fn,OPEN_AP))) {
                printf1("Error: can't append to output file: %s\n",fn);
                return(-1);       
            }
        }
    }
    if (s != NULL) {
        if (fn == NULL)
            printf1("%s\n",s);
        else
            fprintf(fd,"%s\n",s);
    }
    for (i = 1; i <= m; ++i) {
        ii = (i - 1) * n;
        for (j = 1; j <= n; ++j) {
            if (fn == NULL) {
                if (ivflg == 2)
                    printf1("[");
                printf1(fmt,a[ii + j]);
                if (ivflg) {
                    printf1(fmt,b[ii + j]);
                    if (ivflg == 2)
                        printf1("]");
                }
            }
            else {
                if (ivflg == 2)
                    fprintf(fd,"[");
                fprintf(fd,fmt,a[ii + j]);
                if (ivflg) {
                    fprintf(fd,fmt,b[ii + j]);
                    if (ivflg == 2)
                        fprintf(fd,"]");
                }
            }
        }
        if (fn == NULL)  
            printf1("\n");
        else
            fprintf(fd,"\n");

    }
    if (fn != NULL)
        fclose(fd);

    return(0);
}

/*--------------------------------------------------------------------------*/
/*  mprf(m,n,a,fmt,fd)      print m,n matrix a with fmt to file fd          */
/*                          return 0 if successful, -1 if error.            */

int mprf(int m,int n,float *a,char *fmt,FILE *fd)
{
    register int i,j,ii;

    for (i = 1; i <= m; ++i) {
        ii = (i - 1) * n;
        for (j = 1; j <= n; ++j)  
            fprintf(fd,fmt,a[ii + j]);
        fprintf(fd,"\n");
    }
    return(0);
}
   
/*--------------------------------------------------------------------------*/
/*  mprd(m,n,a,fmt,fd)      print m,n matrix a with fmt to file fd          */
/*                          return 0 if successful, -1 if error.            */

int mprd(int m,int n,double *a,char *fmt,FILE *fd)
{
    register int i,j,ii;

    for (i = 1; i <= m; ++i) {
        ii = (i - 1) * n;
        for (j = 1; j <= n; ++j)  
            fprintf(fd,fmt,a[ii + j]);
        fprintf(fd,"\n");
    }
    return(0);
}
   
/* ------------------------------------------------------------------------ */
/*  skip_b(p) Skip blanks and tabs at pointer p and return pointer to next  */
/*            non-blank character.                                          */

char *skip_b(char *p)
{
    while (*p && (*p == ' ' || *p == '\t'))
        p++;
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  skip_c(p) Skip characters at pointer p until next blank or tab          */
/*            character.                                                    */

char *skip_c(char *p)
{
    while (*p && *p != ' ' && *p != '\t')
        p++;
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  skip_cb(p)  Skip characters at pointer p, first until next blank        */
/*              or tab character, then until next non-blank or non-tab.     */

char *skip_cb(char *p)
{
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

char *skip_int(char *p)
{
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

char *skip_dbl(char *p)
{
   
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

char *skip_expr(char *p)
{
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

char *skip_blev(char *p)
{
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

char *skip_nc(char *p)
{
    while (*p) {
        if (*p == ',' || *p == ')')
            return(p);
        else if (*p == '(' || *p == '[' || *p == '{' || *p == '<')
            p = skip_blev(p);
        else
            p++;
    }
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  skip_com(p)     skip until comma or ')'                                 */

char *skip_com(char *p)
{
    while (*p && *p != ',' && *p != ')')
        p++;
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  skip_xa(p)  skip xa(...)=varlist string                                 */

char *skip_xa(char *p)
{
    register int l;

    p += 2;
    p = skip_blev(p);
    if (*p != '=')
        return(p);

    if (*(p + 1) == '1' && *(p + 2) == ',')
        return(p + 2);
           
    while (*(p + 1)) {
        l = get_vnlen(p + 1);
        if (l <= 0)
            break;
        p += l + 1;
    }
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  get_fname(p,q)          get file name beginning at p and save in q.     */
/*                          if successful return 1, else 0.                 */

int get_fname(char *p,char *q)
{
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

/* ------------------------------------------------------------------------ */
/*  printf1()       General printf function that recognizes SILENTFlg.      */

void printf1(const char *fmt, ...) 
{
    va_list ap;

    if (SILENTFlg == 1 || SILENTFlg == 3)
        return;

    va_start(ap,fmt);
    vprintf(fmt,ap);
    va_end(ap);
}

/* ------------------------------------------------------------------------ */
/*  printf2()       General printf function that recognizes SILENTFlg.      */

void printf2(const char *fmt, ...) 
{
    va_list ap;

    if (SILENTFlg >= 0)
        return;

    va_start(ap,fmt);
    vprintf(fmt,ap);
    va_end(ap);
}

/* ------------------------------------------------------------------------ */
/*  printfe()       fprintf to stderr.                                      */

#define FE_BufLen 132
char FE_Buf[FE_BufLen + 1];
int FE_BufPtr = 0;

void printfe(const char *fmt, ...) 
{
    va_list ap;

    va_start(ap,fmt);
    vfprintf(stderr,fmt,ap);
    va_end(ap);
}

/* ------------------------------------------------------------------------ */
/*  fflushe()     fflush stderr.                                            */

void fflushe(void)
{
    fflush(stderr);
}




