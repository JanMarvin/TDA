/****************************************************************************/
/*  t_help                                                                  */
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
#include "t_gf.h"
#include "t_pgen.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_help.c                                                   */

int help(char *key);
int prn_keys(char *key);
int init_hlp(int opt);
int alloc_hlp(int opt);
int s_match(char *p, char *s);
void hlp_norm(char *s);
int hlp_inp(void);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

#define MaxHlpKey 2000  /* max number of help keywords                      */
#define HlpPLen 24      /* number of lines per page                         */

#define HlpLLen 78      /* number of characters per line                    */
#define FNHLP "tda.hlp" /* name of tda help file                            */

FILE *FdHlp;            /* file handle for doc file                         */
int FdHlpFnd = 0;       /* set to 1 if successfully opened                  */
int HlpInit = 0;        /* set to 1 if doc keys initialized.                */
char **HlpKey;          /* array of key words                               */
int HlpKeyA = 0;        /* allocated                                        */
long *HlpPtr;           /* array of file pointers                           */
int HlpPtrA = 0;        /* allocated                                        */
int NHlpKey = 0;        /* number of key words.                             */
int MaxKeyLen = 0;      /* max key word length                              */
char HlpLKey[81];       /* last used key word                               */

/* ------------------------------------------------------------------------ */
/*  help(key)   look for key word in help file; if key contains wildcards   */
/*              print matching key words, otherwise the entry in doc file.  */
/*              Return 0 if OK, -1 if error.                                */
/*                                                                          */

int help(char *key)
{
    register int i;
    int n,nl,r,nlflg,ns,nn,fin,len;
    long fptr;
    char buf[201];
    register char *p;
       
    if (HlpInit == 0) {     /* initialization */
        if (init_hlp(1))
            return(-1);  
        printf1("Initialization of help, found %d key words.\n",NHlpKey);
        HlpLKey[0] = '\0';
    }
    hlp_norm(key);
    if (*key == '.')
        key = HlpLKey;
    else {
        if (!*key)
            strcpy(key,"help");
        strncpy(HlpLKey,key,80);
        *(HlpLKey + 80) = '\0';
    }
    n = prn_keys(key);
    if (n == 0) {
        strcat(key,"*");
        n = prn_keys(key);  
    }
    if (n == 0)
        printf1("No matching entry for: %s\n",key);
    if (n != 1)
        return(0);

    fptr = 0;
    len = fin = ns = nl = 0;
    nlflg = n = 0;
    for (i = 0; i < NHlpKey; ++i) {
        if (s_match(key,HlpKey[i])) {
            fflush(stdout);
            fseek(FdHlp,HlpPtr[i],0);

            while (fgets(buf,200,FdHlp)) {
                if (!strncmp(buf,"##",2)) { 
                    n = 1; 
                    fptr = ftell(FdHlp);
                }
                else if (*buf != '#') { 
                    n++;
                    break;
                }
            }
            if (n != 2)
                return(0);

            while (1) {
                if (n == 2) {
                    prnchar('-',LLEN,1);    
                    nlflg = 0;
                    nl = 1;
                    n = 3;
                    ns = 1;
                }
                p = buf + strlen(buf);
                while (--p >= buf && (*p == CR || *p == LF || *p == ' ')) ;
                if (p < buf) 
                    *++p = ' ';
                *++p = '\0';
                printf1("%s",buf);
                len = strlen(buf);
                nlflg = 1;
                nl++;
                ns++;
                fin = 1;
                while (fgets(buf,200,FdHlp)) {
                    if (*buf != '#') { 
                        fin = 0;
                        break;
                    }
                    else if (!strncmp(buf,"##",2))   
                        break;
                }
                if (fin || nl >= HlpPLen) {
                    if (fin && ns < HlpPLen)  
                        break;     
                          
                    prnchar(' ',HlpLLen - 4 - len,0);
                    r = hlp_inp();      /* wait for user input */
                    if (r == -1) {      /* previous page */
                        fseek(FdHlp,fptr,0);
                        nn = imax(0,ns - 2 * HlpPLen);
                        ns = 0;
                        while (fgets(buf,200,FdHlp)) {
                            if (*buf != '#') {
                                ns++;    
                                if (--nn < 0)
                                    break;
                            }   
                        }
                        prnchar('-',LLEN,1);    
                        nlflg = 0;
                        nl = 1;
                    }
                    else {
                        if (fin || r == 0)      /* quit */
                            return(0);
                        else if (r == 1)        /* next line */
                            nl = HlpPLen - 1;
                        else                    /* next page */
                            nl = 0;
                    }
                }
                else if (nlflg) {
                    newline();
                    nlflg = 0;
                }
            }
            if (nlflg)
                newline();
            break;
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prn_keys(key)   Print available key words that match key.               */
/*                  If just one key no print. Return number of keys.        */  

int prn_keys(char *key)
{
    register int i;
    int n,rcnt,len,cnt,nl,r,nlflg;
      
    cnt = len = 0;
    for (i = 0; i < NHlpKey; ++i) {
        if (s_match(key,HlpKey[i])) {
            cnt++;
            n = strlen(HlpKey[i]);
            if (len < n)
                len = n;
        }
    }
    if (cnt <= 1)  
        return(cnt);
      
    len++;       
    n = HlpLLen / (len + 1);
    nlflg = nl = 0;
    rcnt = 0;
    for (i = 0; i < NHlpKey; ++i) {
        if (s_match(key,HlpKey[i])) {
            printf1("%s ",HlpKey[i]);
            nlflg = 1;

            if (++rcnt >= n) {
                rcnt = 0;             
                if (++nl >= HlpPLen) {
                    r = hlp_inp();
                    if (r == 0)             /* quit */
                        return(cnt);
                    else if (r == 1)        /* next line */
                        nl = HlpPLen - 1;
                    else                    /* next page */
                        nl = 0;
                }
                else {
                    newline();
                    nlflg = 0;
                }
            }
            else
                prnchar(' ',len - strlen(HlpKey[i]),0);
        }
    }
    if (nlflg)
        newline();
    return(cnt);
}

/* ------------------------------------------------------------------------ */
/*  init_hlp(opt)   If opt != 0 initialization of help data structure.      */
/*                  Otherwise free previously allocated memory and close    */  
/*                  the help file.                                          */
/*                  Return 0 if OK, -1 if error.                            */

int init_hlp(int opt)
{
    register char *p,*q;
    int err,l;
    long fptr;
    char buf[501];

    if (opt == 0) {
        if (HlpInit) {
            alloc_hlp(0);
            if (FdHlpFnd)  
                fclose(FdHlp);
            FdHlpFnd = HlpInit = 0; 
        }
        return(0);
    }
    err = -1;

    if ((FdHlp = fopen(FNHLP,OPEN_RD)))   
        FdHlpFnd = 1;

    /* try to find the help file in any directory with executables. */
    
    if (FdHlpFnd == 0) {
     
        p = (char *) getenv("PATH");
        /* printf1("p=:%s:\n",p); */
        while (*p) {
            if (*p == PathSepC)
                p++;
            q = buf;
            while (*p && *p != PathSepC)      
                *q++ = *p++;
            *q++ = PSepC;
            *q = '\0';
            strcat(q,FNHLP);
            /** printf1("buf=:%s:\n",buf); **/
            if ((FdHlp = fopen(buf,OPEN_RD))) {
                FdHlpFnd = 1;
                break;    
            }
        }
        if (FdHlpFnd == 0) {
            printf1("Error: can't find help file: %s.\n",FNHLP);
            goto HLPIFin;
        }
    }
    if (alloc_hlp(1))       /* allocate memory */
        goto HLPIFin;

    /* read key words from doc file */

    fptr = ftell(FdHlp);
    NHlpKey = 0;

    while (fgets(buf,200,FdHlp)) {

        p = buf;
        if (*p++ == '#' && *p++ == '#') {
            hlp_norm(p);
            l = strlen(p);
            if (l > 0) {
                if (NHlpKey >= MaxHlpKey) {
                    printf1("Warning: exceeded max number of doc key words (%d).\n",
                                                                 MaxHlpKey);
                    break;
                }
                if (!(HlpKey[NHlpKey] = (char *)calloc(l + 1,sizeof(char)))) {
                    p_err(-2,1);
                    goto HLPIFin;
                }
                memrq(l + 1,sizeof(char));
                strcpy(HlpKey[NHlpKey],p);
                HlpPtr[NHlpKey] = fptr;
                NHlpKey++;
                if (MaxKeyLen < l)
                    MaxKeyLen = l;
            }
        }
        fptr = ftell(FdHlp);
    }
    HlpInit = 1;
    err = 0;

HLPIFin:
    if (err) {
        alloc_hlp(0);
        if (FdHlpFnd)  
            fclose(FdHlp);
        FdHlpFnd = HlpInit = 0; 
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  alloc_hlp(opt)  If opt != 0 allocate memory for MaxHlpKey keywords.     */
/*                  otherwise free previously allocated memory.             */
/*                  Return 0 if OK, -1 if error.                            */

int alloc_hlp(int opt)
{
    int i,l,err = 0;

    if (opt == 0)
        goto AHLPFin;

    err = -1;
    if (!(HlpKey = (char **)calloc(MaxHlpKey,sizeof(char *)))) {
        p_err(-2,1);
        goto AHLPFin;
    }
    HlpKeyA = MaxHlpKey;
    memrq(MaxHlpKey,sizeof(char *));

    if (!(HlpPtr = (long *)calloc(MaxHlpKey,sizeof(long)))) {
        p_err(-2,1);     
        goto AHLPFin;   
    }
    HlpPtrA = MaxHlpKey;
    memrq(MaxHlpKey,sizeof(long));
    return(0);

AHLPFin:
    if (HlpKeyA) {
        for (i = 0; i < NHlpKey; ++i) {
            l = strlen(HlpKey[i]) + 1;
            free((char *)HlpKey[i]);
            memrq(-l,sizeof(char));
        }
        free((char *)HlpKey);
        memrq(-HlpKeyA,sizeof(char *));
        HlpKeyA = 0;
    }
    if (HlpPtrA) {
        free((char *)HlpPtr);
        memrq(-HlpPtrA,sizeof(long));
        HlpPtrA = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  s_match(p,s)    checks whether the string p matches the string s.       */
/*                  p may contain wildcards (* and ?).                      */
/*                  return 1 if match, 0 otherwise.                         */
/*                                                                          */
/*  Code adapted from: Mike Cornelison, Two Wildcard Matching Utilities,    */
/*  in: C/C++ Users Journal, Vol. 13, No. 4 (April 1995), pp. 55--60.       */

int s_match(char *p, char *s)
{
    register int i,star;

SM_New:
    star = 0;
    while (*p == '*') {
        star = 1;
        p++;
    }
SM_Test:
    for (i = 0; p[i] && (p[i] != '*'); i++) {
        if (p[i] != s[i]) {
            if (!s[i])
                return(0);
            if (p[i] == '?')
                continue;
            if (!star)
                return(0);
            s++;
            goto SM_Test;
        }
    }
    if (p[i] == '*') {
        s += i;
        p += i;
        goto SM_New;
    }
    if (!s[i])
        return(1);
    if (i && p[i - 1] == '*')
        return(1);
    if (!star)
        return(0);
    s++;
    goto SM_Test;
}

/* ------------------------------------------------------------------------ */
/*  hlp_norm(s)     translate to lower case and remove multiple blanks.     */
/*                  skip blanks at beginning and end.                       */

void hlp_norm(char *s)
{
    register char *p,*q;
    register int n = 1;
          
    q = p = s;
    while (*q) {
        if (*q == LF || *q == CR) {
            *p = '\0';
            break;
        }
        else if (*q == ' ') {
            n++;
            if (n == 1)
                *p++ = ' ';
        }
        else {
            n = 0;
            *p++ = (char)tolower((int)*q);
        }
        q++;
    }
    if (p > s && *(p - 1) == ' ')
        p--;
    *p = '\0';
}

/* ------------------------------------------------------------------------ */
/*  hlp_inp()   Get user input.                                             */
/*  return:  -1  -, return                      previous page               */
/*            0  q, return                      quit                        */
/*            1  return                         next line                   */
/*            2  any character, return          next page                   */

int hlp_inp(void)
{
    int c;

    printf(" ");
    fflush(stdout);
    c = getc(stdin);
    fflush(stdin);
    if (c == '\n')
        return(1);
    else if (c == 'q')
        return(0);
    else if (c == '-')
        return(-1);
    return(2);
}


