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
#include "tda_context.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_help.c                                                   */

int help(TDAContext *ctx, char *key);
int prn_keys(TDAContext *ctx, char *key);
int init_hlp(TDAContext *ctx, int opt);
int alloc_hlp(TDAContext *ctx, int opt);
int s_match(TDAContext *ctx, char *p, char *s);
void hlp_norm(TDAContext *ctx, char *s);
int hlp_inp(TDAContext *ctx);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

#define MaxHlpKey 2000  /* max number of help keywords                      */
#define HlpPLen 24      /* number of lines per page                         */

#define HlpLLen 78      /* number of characters per line                    */
#define FNHLP "tda.hlp" /* name of tda help file                            */


/* ------------------------------------------------------------------------ */
/*  help(key)   look for key word in help file; if key contains wildcards   */
/*              print matching key words, otherwise the entry in doc file.  */
/*              Return 0 if OK, -1 if error.                                */
/*                                                                          */

int help(TDAContext *ctx, char *key)
{
    register int i;
    int n,nl,r,nlflg,ns,nn,fin,len;
    long fptr;
    char buf[201];
    register char *p;
       
    if (ctx->HlpInit == 0) {     /* initialization */
        if (init_hlp(ctx, 1))
            return(-1);  
        printf1(ctx, "Initialization of help, found %d key words.\n",ctx->NHlpKey);
        ctx->HlpLKey[0] = '\0';
    }
    hlp_norm(ctx, key);
    if (*key == '.')
        key = ctx->HlpLKey;
    else {
        if (!*key)
            strcpy(key,"help");
        strncpy(ctx->HlpLKey,key,80);
        *(ctx->HlpLKey + 80) = '\0';
    }
    n = prn_keys(ctx, key);
    if (n == 0) {
        strcat(key,"*");
        n = prn_keys(ctx, key);  
    }
    if (n == 0)
        printf1(ctx, "No matching entry for: %s\n",key);
    if (n != 1)
        return(0);

    fptr = 0;
    len = fin = ns = nl = 0;
    nlflg = n = 0;
    for (i = 0; i < ctx->NHlpKey; ++i) {
        if (s_match(ctx, key,ctx->HlpKey[i])) {
            tda_out_flush();
            fseek(ctx->FdHlp,ctx->HlpPtr[i],0);

            while (fgets(buf,200,ctx->FdHlp)) {
                if (!strncmp(buf,"##",2)) { 
                    n = 1; 
                    fptr = ftell(ctx->FdHlp);
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
                    prnchar(ctx, '-',LLEN,1);    
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
                printf1(ctx, "%s",buf);
                len = (int)(strlen(buf));
                nlflg = 1;
                nl++;
                ns++;
                fin = 1;
                while (fgets(buf,200,ctx->FdHlp)) {
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
                          
                    prnchar(ctx, ' ',HlpLLen - 4 - len,0);
                    r = hlp_inp(ctx);      /* wait for user input */
                    if (r == -1) {      /* previous page */
                        fseek(ctx->FdHlp,fptr,0);
                        nn = imax(ctx, 0,ns - 2 * HlpPLen);
                        ns = 0;
                        while (fgets(buf,200,ctx->FdHlp)) {
                            if (*buf != '#') {
                                ns++;    
                                if (--nn < 0)
                                    break;
                            }   
                        }
                        prnchar(ctx, '-',LLEN,1);    
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
                    newline(ctx);
                    nlflg = 0;
                }
            }
            if (nlflg)
                newline(ctx);
            break;
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prn_keys(key)   Print available key words that match key.               */
/*                  If just one key no print. Return number of keys.        */  

int prn_keys(TDAContext *ctx, char *key)
{
    register int i;
    int n,rcnt,len,cnt,nl,r,nlflg;
      
    cnt = len = 0;
    for (i = 0; i < ctx->NHlpKey; ++i) {
        if (s_match(ctx, key,ctx->HlpKey[i])) {
            cnt++;
            n = (int)(strlen(ctx->HlpKey[i]));
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
    for (i = 0; i < ctx->NHlpKey; ++i) {
        if (s_match(ctx, key,ctx->HlpKey[i])) {
            printf1(ctx, "%s ",ctx->HlpKey[i]);
            nlflg = 1;

            if (++rcnt >= n) {
                rcnt = 0;             
                if (++nl >= HlpPLen) {
                    r = hlp_inp(ctx);
                    if (r == 0)             /* quit */
                        return(cnt);
                    else if (r == 1)        /* next line */
                        nl = HlpPLen - 1;
                    else                    /* next page */
                        nl = 0;
                }
                else {
                    newline(ctx);
                    nlflg = 0;
                }
            }
            else
                prnchar(ctx, ' ',len - (int)strlen(ctx->HlpKey[i]),0);
        }
    }
    if (nlflg)
        newline(ctx);
    return(cnt);
}

/* ------------------------------------------------------------------------ */
/*  init_hlp(opt)   If opt != 0 initialization of help data structure.      */
/*                  Otherwise free previously allocated memory and close    */  
/*                  the help file.                                          */
/*                  Return 0 if OK, -1 if error.                            */

int init_hlp(TDAContext *ctx, int opt)
{
    register char *p,*q;
    int err,l;
    long fptr;
    char buf[501];

    if (opt == 0) {
        if (ctx->HlpInit) {
            alloc_hlp(ctx, 0);
            if (ctx->FdHlpFnd)  
                fclose(ctx->FdHlp);
            ctx->FdHlpFnd = ctx->HlpInit = 0; 
        }
        return(0);
    }
    err = -1;

    if ((ctx->FdHlp = fopen(FNHLP,OPEN_RD)))   
        ctx->FdHlpFnd = 1;

    /* TDA_HLP names the file outright, for anyone who keeps it elsewhere. */

    if (ctx->FdHlpFnd == 0) {
        p = (char *) getenv("TDA_HLP");
        if (p && *p && strlen(p) < sizeof(buf) - 1) {
            strcpy(buf,p);
            if ((ctx->FdHlp = fopen(buf,OPEN_RD)))
                ctx->FdHlpFnd = 1;
        }
    }

    /* Next to the executable, and in ../doc beside it.  tda.hlp historically
       sits with the binary, but only the current directory used to be
       searched, so starting TDA from anywhere else lost it. */

    if (ctx->FdHlpFnd == 0 && ctx->ExeDir[0]) {
        strcpy(buf,ctx->ExeDir);
        strcat(buf,FNHLP);
        if ((ctx->FdHlp = fopen(buf,OPEN_RD)))
            ctx->FdHlpFnd = 1;
    }
    if (ctx->FdHlpFnd == 0 && ctx->ExeDir[0]) {
        strcpy(buf,ctx->ExeDir);
        strcat(buf,"..");
        buf[strlen(buf) + 1] = '\0';
        buf[strlen(buf)] = PSepC;
        strcat(buf,"doc");
        buf[strlen(buf) + 1] = '\0';
        buf[strlen(buf)] = PSepC;
        strcat(buf,FNHLP);
        if ((ctx->FdHlp = fopen(buf,OPEN_RD)))
            ctx->FdHlpFnd = 1;
    }

    /* try to find the help file in any directory with executables. */
    
    if (ctx->FdHlpFnd == 0) {
     
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
            if ((ctx->FdHlp = fopen(buf,OPEN_RD))) {
                ctx->FdHlpFnd = 1;
                break;    
            }
        }
        if (ctx->FdHlpFnd == 0) {
            printf1(ctx, "Error: can't find help file: %s.\n",FNHLP);
            goto HLPIFin;
        }
    }
    if (alloc_hlp(ctx, 1))       /* allocate memory */
        goto HLPIFin;

    /* read key words from doc file */

    fptr = ftell(ctx->FdHlp);
    ctx->NHlpKey = 0;

    while (fgets(buf,200,ctx->FdHlp)) {

        p = buf;
        if (*p++ == '#' && *p++ == '#') {
            hlp_norm(ctx, p);
            l = (int)(strlen(p));
            if (l > 0) {
                if (ctx->NHlpKey >= MaxHlpKey) {
                    printf1(ctx, "Warning: exceeded max number of doc key words (%d).\n",
                                                                 MaxHlpKey);
                    break;
                }
                if (!(ctx->HlpKey[ctx->NHlpKey] = (char *)calloc((size_t)(l + 1),sizeof(char)))) {
                    p_err(ctx, -2,1);
                    goto HLPIFin;
                }
                memrq(ctx, l + 1,sizeof(char));
                strcpy(ctx->HlpKey[ctx->NHlpKey],p);
                ctx->HlpPtr[ctx->NHlpKey] = fptr;
                ctx->NHlpKey++;
                if (ctx->MaxKeyLen < l)
                    ctx->MaxKeyLen = l;
            }
        }
        fptr = ftell(ctx->FdHlp);
    }
    ctx->HlpInit = 1;
    err = 0;

HLPIFin:
    if (err) {
        alloc_hlp(ctx, 0);
        if (ctx->FdHlpFnd)  
            fclose(ctx->FdHlp);
        ctx->FdHlpFnd = ctx->HlpInit = 0; 
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  alloc_hlp(opt)  If opt != 0 allocate memory for MaxHlpKey keywords.     */
/*                  otherwise free previously allocated memory.             */
/*                  Return 0 if OK, -1 if error.                            */

int alloc_hlp(TDAContext *ctx, int opt)
{
    int i,l,err = 0;

    if (opt == 0)
        goto AHLPFin;

    err = -1;
    if (!(ctx->HlpKey = (char **)calloc(MaxHlpKey,sizeof(char *)))) {
        p_err(ctx, -2,1);
        goto AHLPFin;
    }
    ctx->HlpKeyA = MaxHlpKey;
    memrq(ctx, MaxHlpKey,sizeof(char *));

    if (!(ctx->HlpPtr = (long *)calloc(MaxHlpKey,sizeof(long)))) {
        p_err(ctx, -2,1);     
        goto AHLPFin;   
    }
    ctx->HlpPtrA = MaxHlpKey;
    memrq(ctx, MaxHlpKey,sizeof(long));
    return(0);

AHLPFin:
    if (ctx->HlpKeyA) {
        for (i = 0; i < ctx->NHlpKey; ++i) {
            l = (int)(strlen(ctx->HlpKey[i]) + 1);
            free((char *)ctx->HlpKey[i]);
            memrq(ctx, -l,sizeof(char));
        }
        free((char *)ctx->HlpKey);
        memrq(ctx, -ctx->HlpKeyA,sizeof(char *));
        ctx->HlpKeyA = 0;
    }
    if (ctx->HlpPtrA) {
        free((char *)ctx->HlpPtr);
        memrq(ctx, -ctx->HlpPtrA,sizeof(long));
        ctx->HlpPtrA = 0;
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

int s_match(TDAContext *ctx, char *p, char *s)
{
    (void)ctx;        /* unused: the signature is shared */
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

void hlp_norm(TDAContext *ctx, char *s)
{
    (void)ctx;        /* unused: the signature is shared */
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

int hlp_inp(TDAContext *ctx)
{
    (void)ctx;        /* unused: the signature is shared */
#ifdef TDA_R_PACKAGE
    /*  Paging waits on stdin, and under R there is nobody at a keyboard to
        press return: the run would hang.  Behave as if the reader had asked
        for the next page, so the whole entry comes out at once and the R
        side can hand it back as text. */
    return(1);
#else
    int c;

    tda_out(" ");
    tda_out_flush();
    c = getc(stdin);
    fflush(stdin);
    if (c == '\n')
        return(1);
    else if (c == 'q')
        return(0);
    else if (c == '-')
        return(-1);
    return(2);
#endif
}


