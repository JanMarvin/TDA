/****************************************************************************/
/*  t_var                                                                   */
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
#include "t_eval.h"
#include "t_eval2.h"
#include "t_eval3.h"
#include "t_gdat.h"
#include "t_parm.h"
#include "t_gf.h"
#include "t_alloc.h"
#include "t_freq.h"
#include "t_mdat.h"
#include "t_mat.h"
#include "t_sd.h"
#include "tda_context.h"

/*  functions in t_var.c */

char get_vnchar(TDAContext *ctx, char c);
int check_vname(TDAContext *ctx, char *p);
int get_vnlen(TDAContext *ctx, char *p);
int get_vnl(TDAContext *ctx, int n,short *vinum);
int alloc_vmax(TDAContext *ctx);
int get_nidx(TDAContext *ctx);
void clear_vidx(TDAContext *ctx, int i);
int get_vidx(TDAContext *ctx, char *p,char *name);
int get_vidx1(TDAContext *ctx, char *p,char *name);
int get_nvidx(TDAContext *ctx, char *p,short *vinum);
int save_var(TDAContext *ctx, char *vd,int opt);
int sub_ifs(TDAContext *ctx, char *vd,char *buf); 
void var_err(TDAContext *ctx, int n,char *vd);
int alloc_vdat(TDAContext *ctx, int idx,int opt);
int alloc_vsdat(TDAContext *ctx, int n,short *ivar);
int get_slen(TDAContext *ctx, int j,int noc);
void prn_var(TDAContext *ctx, int idx);
void prn_vname(TDAContext *ctx, int i);
void prn_vlabel(TDAContext *ctx, int i);
void prn_hvar(TDAContext *ctx);
void prn_hlabel(TDAContext *ctx);
void clear_dm(TDAContext *ctx);
void clear_var(TDAContext *ctx, int idx);
void clear_avar(TDAContext *ctx, int idx);
void free_var(TDAContext *ctx, int i);
int nlist(TDAContext *ctx);
int nl_check(TDAContext *ctx, char *nl);
char *get_nvi(TDAContext *ctx, char *s,int *n,int nv,short *vidx,int *nb);
char *get_nvia(TDAContext *ctx, char *s,int *n,int opt,int *nb);
void nl_free(TDAContext *ctx, int j);
int alloc_vl(TDAContext *ctx, int n);
int recode(TDAContext *ctx);
int ndvar(TDAContext *ctx);
int ndivar(TDAContext *ctx, char *vname,int n1,short *idx1,int n2,short *idx2);
int c_nlist(TDAContext *ctx, char *vname,int n,short *vidx);
int get_mxvlen(TDAContext *ctx, int n,short *vidx);
void prn_vlist(TDAContext *ctx, int n,short *vidx);
int svb_alloc(TDAContext *ctx, int n);

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */





                        /* available. 2 if archive variable and data not    */  
                        /* not already created.                             */



                        /* flag if at least one var label present.          */






/*--------------------------------------------------------------------------*/
/*  get_vnchar(c)   If c is valid character for variable name return c,     */
/*                  otherwise return underline.                             */

char get_vnchar(TDAContext *ctx, char c)
{
    (void)ctx;        /* unused: the signature is shared */
    if (c >= 'a' && c <= 'z')        
        return(c);
    if (c >= 'A' && c <= 'Z')        
        return(c);
    if (isdigit((int)c))
        return(c);
    if (c == '_' || c == '@' || c == '$')
        return(c);                 
    return('_');
}

/*--------------------------------------------------------------------------*/
/*  check_vname(p)  check whether p points to the beginning of a variable   */
/*                  name. Return 1 if this is the case, otherwise 0.        */
/*                                                                          */
/*                  A,...,Z,_,@,$                                           */

int check_vname(TDAContext *ctx, char *p)
{
    (void)ctx;        /* unused: the signature is shared */
    if ((*p >= 'A' && *p <= 'Z') || *p == '_' || *p == '@' || *p == '$')
        return(1);
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  get_vnlen(p)    get length of variable name pointed to by p.            */
/*                  return this length.                                     */

int get_vnlen(TDAContext *ctx, char *p)
{
    register int l;
           
    if (check_vname(ctx, p) == 0)
        return(0);
    p++;
    l = 1;
    while ((*p >= 'A' && *p <= 'Z') || (*p >= 'a' && *p <= 'z') ||  
                   isdigit((int)*p) || *p == '_' || *p == '@' || *p == '$') {
        l++;
        p++;
    }
    return(l);
}

/* ------------------------------------------------------------------------ */
/*  get_vnl(n,vinum)  return max length of variable names in vinum[]        */

int get_vnl(TDAContext *ctx, int n,short *vinum)
{
    register int i,j,l;

    l = 0;
    for (i = 0; i < n; ++i) {
        j = (int)(strlen(ctx->VName[vinum[i]]));
        if (l < j)
            l = j;
    }
    return(l);
}

/*--------------------------------------------------------------------------*/
/*  alloc_vmax()    Allocate memory for MaxNV variables. This memory is     */
/*                  only allocated once at the beginning of the program     */
/*                  and is never free'd. Note, we also allocate a global    */
/*                  stack for the parser.                                   */
/*                  return 0 if OK, -1 if error.                            */

int alloc_vmax(TDAContext *ctx)
{
    int i,n1,n2;

    n1 = (int)((size_t)(4 * ctx->MaxNV) * sizeof(char *) +
               (size_t)ctx->MaxNV * sizeof(short));  

    if (!(ctx->VDef   = (char **)calloc((size_t)(ctx->MaxNV),sizeof(char *))) ||
        !(ctx->VName  = (char **)calloc((size_t)(ctx->MaxNV),sizeof(char *))) ||
        !(ctx->VLabel = (char **)calloc((size_t)(ctx->MaxNV),sizeof(char *))) || 
        !(ctx->VDPtr  = (char **)calloc((size_t)(ctx->MaxNV),sizeof(char *))) || 
        !(ctx->VStrN  = (short *)calloc((size_t)(ctx->MaxNV),sizeof(short )))) {  

        p_err(ctx, -2,1);    
        return(-1);
    }         
    memrq(ctx, n1,1);

    n2 = (int)((size_t)ctx->MaxNV *
               (3 * sizeof(char) + 5 * sizeof(short) + sizeof(char **) +
                sizeof(int *) + sizeof(double) + sizeof(double *)));

    if (!(ctx->VAlloc = (char *)calloc((size_t)(ctx->MaxNV),sizeof(char))) ||
        !(ctx->VTyp   = (char *)calloc((size_t)(ctx->MaxNV),sizeof(char))) ||
        !(ctx->VTypA  = (char *)calloc((size_t)(ctx->MaxNV),sizeof(char))) ||
        !(ctx->VSLen  = (short *)calloc((size_t)(ctx->MaxNV),sizeof(short))) ||
        !(ctx->VPFmt1 = (short *)calloc((size_t)(ctx->MaxNV),sizeof(short))) ||
        !(ctx->VPFmt2 = (short *)calloc((size_t)(ctx->MaxNV),sizeof(short))) ||
        !(ctx->VPFmtS = (char **)calloc((size_t)(ctx->MaxNV),sizeof(char *))) ||
        !(ctx->VNxt   = (short *)calloc((size_t)(ctx->MaxNV),sizeof(short))) ||  
        !(ctx->VESCnt = (short *)calloc((size_t)(ctx->MaxNV),sizeof(short))) ||  
        !(ctx->VESTyp = (int **)calloc((size_t)(ctx->MaxNV),sizeof(int *))) ||  
        !(ctx->AVVAL  = (double *)calloc((size_t)(ctx->MaxNV),sizeof(double))) ||
        !(ctx->VESVal = (double **)calloc((size_t)(ctx->MaxNV),sizeof(double *)))) {  
        p_err(ctx, -2,1);   
        return(-1);
    }         
    memrq(ctx, n2,1);

    if (alloc_est(ctx, P_ALLOC,1))
        return(-1);

    for (i = 0; i < ctx->MaxNV; ++i) {

        if (!(ctx->VPFmtS[i] = (char *)calloc(VPFmtSLen,sizeof(char)))) {
            p_err(ctx, -2,1);             
            return(-1);
        }         
        memrq(ctx, VPFmtSLen,1);

        ctx->VLabel[i] = NULL;
        ctx->VNxt[i] = -1;
        ctx->VTyp[i] = ctx->VTypA[i] = 0;
        ctx->VAlloc[i] = (char)(ctx->VESCnt[i] = 0);
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  get_nidx()  get new variable index.                                     */
/*              return index, or -1 if maximum reached.                     */

int get_nidx(TDAContext *ctx)
{
    register int i,j;

    j = -1;
    for (i = 0; i < ctx->MaxNV; ++i) {
        if (ctx->VAlloc[i] == 0) {
            j = i;
            break;
        }
    }
    return(j);
}

/*--------------------------------------------------------------------------*/
/*  clear_vidx(i)   remove variable with index i (0 <= i < MaxNV) from      */
/*                  variable list.                                          */

void clear_vidx(TDAContext *ctx, int i)
{
    register int j,k;

    if (i < 0 || ctx->VIFirst < 0 || ctx->NVAR <= 0)  
        gerr_exit(ctx, 15);

    if (ctx->VAlloc[i] != 1)
        gerr_exit(ctx, 51);
      
    j = k = ctx->VIFirst;
    while (i != j) {
        k = j;
        if (k < 0)
            gerr_exit(ctx, 15);
        j = ctx->VNxt[k];
    }
    ctx->VNxt[k] = ctx->VNxt[i];
    if (ctx->VIFirst == i)
        ctx->VIFirst = ctx->VNxt[i];
    if (ctx->VILast == i) {
        if (ctx->VIFirst == -1)
            ctx->VILast = -1;
        else
            ctx->VILast = (short)(k);
    }
    ctx->VNxt[i] = -1;
    ctx->VAlloc[i] = 0;
    ctx->NVAR--;
}

/* ------------------------------------------------------------------------ */
/*  get_vidx(p,name)    p points to the begin of a variable name in a       */
/*                      string, possibly ended by: , {  < [ ( =             */
/*                      Return variable name in name[].                     */
/*                      Check whether this variable name is part of the     */
/*                      current data matrix. If this is the case return     */
/*                      internal variable number, otherwise -1.             */

int get_vidx(TDAContext *ctx, char *p,char *name)
{
    register int i;
    register char *q;

    q = name;
    for (i = 0; i < VNLMax; ++i) {
        if (!*p || *p == ',' || *p == '<' || *p == '[' || *p == '(' || *p == '=' || *p == '{')
            break;
        *q++ = *p++;
    }
    *q = '\0';

    i = ctx->VIFirst;
    while (i >= 0) {
        if (!strcmp(name,ctx->VName[i]))
            return(i);
        i = ctx->VNxt[i];
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  get_vidx1(p,name)   p points to the begin of a variable name in a       */
/*                      string, possibly ended by: , < [ ( = ) ] :          */
/*                      Return variable name in name[].                     */
/*                      Check whether this variable name is part of the     */
/*                      current data matrix. If this is the case return     */
/*                      internal variable number, otherwise -1.             */

int get_vidx1(TDAContext *ctx, char *p,char *name)
{
    register int i;
    register char *q;

    q = name;
    for (i = 0; i < VNLMax; ++i) {
        if (!*p || *p == ',' || *p == '<' || *p == '[' || *p == '(' || *p == '=' || *p == ')' || *p == ']' || *p == ':')
            break;
        *q++ = *p++;
    }
    *q = '\0';

    i = ctx->VIFirst;
    while (i >= 0) {
        if (!strcmp(name,ctx->VName[i]))
            return(i);
        i = ctx->VNxt[i];
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  get_nvidx(p,vinum)                                                      */
/*                                                                          */
/*  p points to a list of variables names. The function tries to find the   */
/*  corresponding internal variable numbers and puts these numbers into the */
/*  array vinum[]. It returns the number of variables, or -1 if an error,   */
/*  i.e. if a variable is undefined.                                        */

int get_nvidx(TDAContext *ctx, char *p,short *vinum)
{
    register int i,j,k,n;
    register char c,*q;
    int cflg,ia,ib;
    char name[VNLMax + 1];

    cflg = 0;
    ia = ib = -1;

    if (!*p || *p == ',')
        return(-1);

    k = n = 0;
    while (*p) {    
        q = name;
        for (i = 0; i < VNLMax; ++i) {
            if (!*p || *p == ',')
                break;
            *q++ = *p++;
        }
        c = *q;
        *q = '\0';

        j = 1;
        if (*name) {
            i = ctx->VIFirst;
            while (i >= 0) {
                if (!strcmp(name,ctx->VName[i])) {
                    j = 0;
                    break;    
                }
                i = ctx->VNxt[i];
            }
        }
        *q = c;
        if (j)          /* not found */
            return(-1);

        ib = i;
        if (cflg) {
            if (ia >= 0 && ia != ib) {

                i = ctx->VNxt[ia];
                while (i >= 0) {
                    if (!strcmp(name,ctx->VName[i]))  
                        break;    
                    vinum[k++] = (short)(i);
                    n++;
                    i = ctx->VNxt[i];
                }
            }
        }
        n++;
        vinum[k++] = (short)(ib);
        ia = ib;
           
        cflg = 0;
        if (*p == ',')
            p++;
        if (*p == ',') {
            cflg = 1;
            p++;
        }
    }
    return(n);
}

/*--##----------------------------------------------------------------------*/
/*  save_var(vd,opt)                                                        */
/*                                                                          */
/*                  Save information about variable defined by vd.          */
/*                  also check if variable definition could be correctly    */
/*                  parsed. If opt == 1 all  new variables are type 5.      */
/*                  Return 0 if OK,                                         */
/*                  -1 if syntax error,                                     */
/*                  -2 if type 5 variable and opt = 0.                      */
/*                  -3 if insufficient memory,                              */
/*                   1 if max number of variables reached.                  */
/*                                                                          */
/*  This function allocates/sets the following data structures.             */
/*  VName[i]        name of variable (allocated)                            */
/*  VLabel[i]       label, if defined, otherwise NULL (allocated, or NULL)  */
/*  VDef[i]         definition of variable (allocated)                      */
/*  VTyp[i]         type of variable                                        */
/*  VTypA[i]        1 if data already created, 2 if archive variable and    */
/*                  and data not already created.                           */
/*  VSLen[i]        storage size of variable                                */
/*  VPFmt1[i]       print format                                            */
/*  VPFmt2[i]       print format                                            */
/*  VPFmtS[i]       print format string                                     */
/*  VESTyp[i][]     allocated                                               */
/*  VESVal[i][]     allocated                                               */
/*  VESCnt[i]       set to 0 to flag that VESTyp and VESVal are allocated.  */
/*                  note: if VESCnt[i] >  0 at the beginning, this is an    */
/*                  allocation error.                                       */
/*                                                                          */
/*  Finally, if successful, one more variable is counted in NVAR.           */
/*  Note: this function also processes commands that define variables by    */
/*  reference to a data archive, the syntax is                              */
/*      VName ... = A:AVName                                                */
/*  The number of these variables is counted in NVArc1.                     */
/*                                                                          */
/*  We have the following types of variables:                               */
/*                                                                          */
/*  VTyp =  1   string variable                                             */
/*          2   numerical constant                                          */
/*          3   numerical variables which do not (directly or indirectly)   */
/*              invole type 2 operators.                                    */
/*          4   numerical variables that directly or indirectly involve     */
/*              type 2 operators                                            */
/*          5   Special variables created temporarily by edef().            */
/*                                                                          */
/*  Storage size:                                                           */
/*  VSLen = 0   Bit                                                         */
/*          1   Byte                                                        */
/*          2   short integer                                               */
/*          5   long integer                                                */
/*          4   single precision float                                      */  
/*          8   double precision float                                      */  
/*         <0   width of string variable                                    */
/*                                                                          */

int save_var(TDAContext *ctx, char *vd,int opt)
{
    register int j,l;
    register char *p,*q,*s,*pp;
    int err,err1,n,nc,nv,n2,n3,m,mm1,mm1s,mm2 = 0,mm3,blevel,ci,cii,w1,w2;
    int ni,ni1,ni2,ni3,ni3len,nibuflen;
    char *buf = NULL,*nibuf = NULL,name[VNLMax+100],vnum[60];

    err = err1 = 0;
    nibuflen = ni = ni1 = ni2 = ni3 = -1;

VDNXT:

    j = get_nidx(ctx);     /* get new variable index */

    if (j < 0) {
        printf1(ctx, "Error: exceeded max number of variables (%d).\n",ctx->MaxNV);
        return(1); 
    }
    if (ctx->VESCnt[j] > 0) {
        gerr_exit(ctx, 16);
    }
    if (ni == -1) {
        if (get_vidx(ctx, vd,name) >= 0) {
            printf1(ctx, "Already used: %s\n",vd);
            err = -1;
            goto SVFin;
        }
        mm1 = (int)strlen(name);
        mm1s = mm1;
        pp = p = vd + mm1;
         
        /* ### check for special {j,n,k} */

        if (*p == '{') {
            if (sscanf(p,"{%d,%d,%d}",&ni1,&ni2,&ni3) == 3 &&
                                 ni1 >= 0 && ni2 >= 1 && ni3 >= 1) {
                p = skip_int(ctx, p + 1);
                p = skip_int(ctx, p + 1);
                p = skip_int(ctx, p + 1);
            }
            else
                goto SVSERR; 

            if (*p++ != '}')
                goto SVSERR;
            pp = p;
            ni = ni1;
            ci = ni3;

            snprintf(vnum,sizeof(vnum),"c%d",ni3);
            ni3len = (int)(strlen(vnum));

            nibuflen = (int)(strlen(p) + 1000);
            if (!(nibuf = (char *)calloc((size_t)(nibuflen),sizeof(char)))) {
                err = -3;
                goto SVFin1;
            }
            memrq(ctx, nibuflen,1);
        }
    }
    else
        p = pp;

    if (ni >= 0) {
        snprintf(name + mm1s,sizeof(name) - (size_t)mm1s,"%d",ni);
        mm1 = (int)(strlen(name));
        if (mm1 > VNLMax) {
            err1 = 1;
            goto SVSERR;
        }
        if (get_vidx(ctx, name,name) >= 0) {
            printf1(ctx, "Already used: %s\n",name);
            err = -1;
            goto SVFin;
        }
        p = pp;
    }     
    if (ctx->VNameLen < mm1)
        ctx->VNameLen = mm1;
    mm1++;
    if (!(ctx->VName[j] = (char *)calloc((size_t)(mm1),sizeof(char)))) {
        err = -3;
        goto SVFin1;       
    }
    memrq(ctx, mm1,1);
    strcpy(ctx->VName[j],name);

    if (mat_getidx(ctx, ctx->VName[j],0) >= 0) {
        printf1(ctx, "Variable name %s already used for a matrix.\n",ctx->VName[j]);
        err = -1;
        goto SVFin2;
    }

    if (*p == '<') {
        if (sscanf(p,"<%d>",&n) == 1 &&            
                (n == 0 || n == 1 || n == 2 || n == 4 || n == 5 || n == 8 || n < 0))
            ctx->VSLen[j] = (short)(n);
        else {
            var_err(ctx, 1,vd);
            err = -1;
            goto SVFin2;
        }
        p = skip_int(ctx, ++p);
        p++;
    }
    else
        ctx->VSLen[j] = 4;   /* default storage size */

    ctx->VPFmt1[j] =  0;      
    ctx->VPFmt2[j] = -1;

    if (*p == '[') {
        if (sscanf(p,"[%d.%d]",&n,&m) == 2) {
            ctx->VPFmt1[j] = (short)(n);
            ctx->VPFmt2[j] = (short)(m);
            p = skip_int(ctx, ++p);
            p = skip_int(ctx, ++p);
        }
        else if (sscanf(p,"[%d]",&n) == 1 && n >= 0) {
            ctx->VPFmt1[j] = (short)(n);
            ctx->VPFmt2[j] = 0;
            p = skip_int(ctx, ++p);
        }
        else {
            var_err(ctx, 2,vd);
            err = -1;
            goto SVFin2;
        }
        p++;

        w1 = (int)ctx->VPFmt1[j];
        w2 = (int)ctx->VPFmt2[j];
    }
    else {              /* default print format */
        w2 = 0;
        switch (ctx->VSLen[j]) {
            case 0:   w1 =  2; break;
            case 1:   w1 =  4; break;
            case 2:   w1 =  6; break;
            case 5:   w1 = 11; break;
            default:  w1 = 0;  break;
        }
    }
    makefmt(ctx, &w1,&w2,ctx->VPFmtS[j],VPFmtSLen,0,' ',0);
    ctx->VPFmt1[j] = (short)w1;
    ctx->VPFmt2[j] = (short)w2;

    if (*p == '(') {                    /* check for label */
        p++;
        q = p;
        l = 0;
        blevel = 1;
        while (*q) {
            if (*q == '(')
                blevel++;
            else if (*q == ')') {
                if (--blevel == 0)
                    break;
            }
            l++;
            q++;
        }
        /************************
        while (*q && *q != ')') {
            l++;
            q++;
        }
        ***********************/
        if (*q != ')' || l > VLLMax) {
            var_err(ctx, 3,vd);
            err = -1;
            goto SVFin2;
        }
        mm2 = l + 1;
        if (!(ctx->VLabel[j] = (char *)calloc((size_t)(mm2),sizeof(char)))) {
            err = -3;           
            ctx->VLabel[j] = NULL;
            goto SVFin2;
        }
        memrq(ctx, mm2,1);

        /****************
        q = VLabel[j];
        while (*p && *p != ')')  
            *q++ = *p++;
        *q = '\0';
        p++;
        *********/
        s = ctx->VLabel[j];
        while (p < q)
            *s++ = *p++;
        *s = '\0';
        p++;
         
        if (ctx->VLabelLen < l)
            ctx->VLabelLen = l;
    }
    if (*p++ != '=') {
        var_err(ctx, 4,vd);
        err = -1;
        goto SVFin3;
    }
    l = (int)(strlen(p));      /* check variable's definition */
    mm3 = l + 1;
    if (l == 0) {
        var_err(ctx, 5,vd);
        err = -1;
        goto SVFin3;
    }

    l = (int)(strlen(p));      /* check variable's definition */
    mm3 = l + 1;
    if (l == 0) {
        var_err(ctx, 5,vd);
        err = -1;
        goto SVFin3;
    }
    if (ni >= 0) {      /* put var definition into nibuf */
        q = p;
        s = nibuf;

        while (*q) {
            if (sscanf(q,"c%d",&cii) == 1 && cii == ni3) {
                snprintf(s,(size_t)nibuflen - (size_t)(s - nibuf),"c%d",ci);
                s += strlen(s);
                q += ni3len;
            }
            else
                *s++ = *q++;
        }      
        *s = '\0';
        p = nibuf;
        l = (int)(strlen(p));  
        mm3 = l + 1;
    }

    q = p;              /* check for if-then-else expressions */
    n = 0;       
    while (*q) {
        if (!strncmp(q,"then",4)) {
            n = 1;
            break;
        }
        q++;
    }
    if (n) {
        if (!(buf = (char *)calloc((size_t)(2) * (size_t)(l),sizeof(char)))) {
            err = -3;
            goto SVFin3;
        }
        m = sub_ifs(ctx, p,buf);     /* change to if() operators */
        mm3 = m + 1;
        p = buf;
    }
                
    if (!(ctx->VDef[j] = (char *)calloc((size_t)(mm3),sizeof(char)))) {
        err = -3;
        if (n)
            free(buf);
        goto SVFin3;
    }
    memrq(ctx, mm3,1);
    strcpy(ctx->VDef[j],p);
    if (n)
        free(buf);
    /*  p pointed into buf, which has just been freed, and the fourteen
        strncmp/sscanf reads of p below then ran on freed memory --
        a use-after-free on ordinary nvar examples (ed1.cf, fml1.cf, frml1.cf).
        VDef[j] is the copy just taken, with the same contents and a
        lifetime that outlives this function, so the rest of the
        function reads that instead.  */
    p = ctx->VDef[j];

    if (!strncmp(p,"spss",4)) {   /* SPSS export/sav file variable */
        ctx->VTyp[j] = 3;
        if (sscanf(p,"spss(%d)",&n) == 1 && n > 0) {   /* string variable */
            ctx->VTyp[j] = 1;
            ctx->VSLen[j] = (short)(-n);           
            ctx->VStrN[j] = 0;
            if (n > ctx->SVBufLen) {
                if (svb_alloc(ctx, n))
                    goto SVFin3;
            }
        }
        ctx->ESCnt = 1;
        ctx->ESTyp[0] = -1;
        ctx->ESVal[0] = 0.0;
    }
    else if (!strncmp(p,"stata",5)) {   /* STATA file variable */
        ctx->VTyp[j] = 3;
        if (sscanf(p,"stata(%d)",&n) == 1 && n > 0) {   /* string variable */
            ctx->VTyp[j] = 1;
            ctx->VSLen[j] = (short)(-n);           
            ctx->VStrN[j] = 0;
            if (n > ctx->SVBufLen) {
                if (svb_alloc(ctx, n))
                    goto SVFin3;
            }
        }
        ctx->ESCnt = 1;
        ctx->ESTyp[0] = -1;
        ctx->ESVal[0] = 0.0;
    }
    else if (!strncmp(p,"str(",4)) {

        if (sscanf(p,"str(%d,%d)",&n,&m) != 2 || n < 1 || m < n) {
            var_err(ctx, 6,vd);
            err = -1;
            goto SVFin3;
        }
        p = skip_int(ctx, p + 4);
        p = skip_int(ctx, p + 1);
        if (*++p) {
            var_err(ctx, 6,vd);
            err = -1;
            goto SVFin3;
        }
        ctx->VTyp[j] = 1;
        ctx->VSLen[j] = (short)(n - m - 1);
        ctx->VStrN[j] = (short)(n - 1);

        if (-ctx->VSLen[j] > ctx->SVBufLen) {
            if (svb_alloc(ctx, -ctx->VSLen[j]))
                goto SVFin3;
        }
        ctx->ESCnt = 1;
        ctx->ESTyp[0] = ctx->VOFFS + j;
        ctx->ESVal[0] = 0.0;
    }
    else if (*p == 'A' && *(p + 1) == ':') {  /* check for archive variables */

        p += 2;
        n = (int)(strlen(p));
        /* n = get_vnlen(p); */

        if (n < 1 || n + 2 != l) {
            printf1(ctx, "Syntax error in archive variable: %s\n",vd);
            err = -1;
            goto SVFin4; 
        }
        if (ctx->VSLen[j] < 0) {
            ctx->VTyp[j] = 1;
            if (-ctx->VSLen[j] > ctx->SVBufLen) {
                if (svb_alloc(ctx, -ctx->VSLen[j]))
                    goto SVFin3;
            }
        }
        else
            ctx->VTyp[j] = 3;
        ctx->VTypA[j] = 2;
        ctx->ESCnt = 1;
        ctx->ESTyp[0] = ctx->VOFFS + j;
        ctx->ESVal[0] = 0.0;
        ctx->NVArc1++;
    }
    else if (*p == 'M' && *(p + 1) == ':') {  /* check for matrix variables */

        ctx->VTyp[j] = 3;
        ctx->ESCnt = 1;
        ctx->ESTyp[0] = ctx->VOFFS;
        ctx->ESVal[0] = 0.0;
    }
    else {

        ctx->VTyp[j] = 3;
        
        /* check whether the variable could be correctly parsed   
           and save results in a variable-specific parser stack. */          

        ctx->PREVName = ctx->VName[j];    /* special for pre operator */
        ctx->PREVNum = j;

        n = v_parse(ctx, ctx->VDef[j],0);

        ctx->PREVName = NULL;
        ctx->PREVNum = -1;
          
        if (n < 0 || ctx->ESCnt <= 0) {
            printf1(ctx, "Syntax error (%d) in variable: %s\n",-n,vd);
            if (n < 0)
                prn_emsg1(ctx, n);

            err = -1;
            goto SVFin4; 
        }

        /* check for a numerical constant */

        if (ctx->ESCnt == 1 && (ctx->ESTyp[0] == 0 || ctx->ESTyp[0] == 1902)) {
            ctx->VTyp[j] = 2;
            if (opt == 0)
                ctx->VSLen[j] = 8;
        }
        else {      /* check for type 2 and 3 operators */
            
            check_expr(ctx, ctx->ESCnt,ctx->ESTyp,&nv,&nc,&n2,&n3,0);
            if (n3 > 0) {
                if (opt == 0) {
                    p_err(ctx, -40,1);
                    err = -2;
                    goto SVFin4;
                }
                else if (n2 > 0) {
                    printf1(ctx, "Error: type 4 and type 5 variables not compatible.\n");
                    err = -2;
                    goto SVFin4;
                }
                else
                    ctx->VTyp[j] = 5;
            }
            else if (n2 > 0)
                ctx->VTyp[j] = 4;
        }
    }
    if (opt)
        ctx->VTyp[j] = 5;

    /* re-make print format for string variables */

    if (ctx->VTyp[j] == 1) {
        w1 = w2 = 0;
        makefmt(ctx, &w1,&w2,ctx->VPFmtS[j],VPFmtSLen,0,' ',0);
        ctx->VPFmt1[j] = (short)w1;
        ctx->VPFmt2[j] = (short)w2;
    }

    /* copy parser stack for variable */

    if (!(ctx->VESTyp[j] = (int *)calloc((size_t)(ctx->ESCnt),sizeof(int)))) {
        err = -3;         
        goto SVFin4;
    }
    memrq(ctx, ctx->ESCnt,sizeof(int));

    if (!(ctx->VESVal[j] = (double *)calloc((size_t)(ctx->ESCnt),sizeof(double)))) {
        err = -3;         
        goto SVFin5;
    }
    memrq(ctx, ctx->ESCnt,sizeof(double));

    ctx->VESCnt[j] = (short)ctx->ESCnt;
    for (l = 0; l < ctx->ESCnt; ++l) {
        ctx->VESTyp[j][l] = ctx->ESTyp[l];
        ctx->VESVal[j][l] = ctx->ESVal[l];
    }
    ctx->VAlloc[j] = 1;
    if (ctx->VIFirst < 0)
        ctx->VILast = ctx->VIFirst = (short)j;
    else {
        ctx->VNxt[ctx->VILast] = (short)(j);
        ctx->VILast = (short)(j);
    }
    ctx->NVAR++;
    if (ni >= 0 && ni < ni2) {
        ni++;
        ci++;                                
        goto VDNXT;
    }
    if (nibuflen > 0) {
        free(nibuf);
        memrq(ctx, -nibuflen,1);
    }
    return(0);

SVFin5:
    free((char *)ctx->VESTyp[j]);
    memrq(ctx, -ctx->ESCnt,sizeof(int));
SVFin4:
    free(ctx->VDef[j]);
    memrq(ctx, -mm3,1);
SVFin3:
    if (ctx->VLabel[j] != NULL) {
        free(ctx->VLabel[j]);
        memrq(ctx, -mm2,1);            
        ctx->VLabel[j] = NULL;
    }
SVFin2:
    free(ctx->VName[j]);
    memrq(ctx, -mm1,1);
SVFin1:
    if (err == -3)  
        printf1(ctx, "Insufficient memory for data structures of variable: %s\n",vd);
       
SVFin:
    if (nibuflen > 0) {
        free(nibuf);
        memrq(ctx, -nibuflen,1);
    }
    return(err);

SVSERR:
    if (nibuflen > 0) {
        free(nibuf);
        memrq(ctx, -nibuflen,1);
    }
    printf1(ctx, "Syntax error in variable: %s\n",vd);
    if (err1)   
        printf1(ctx, "Exceeded max var name length.\n");
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  sub_ifs(vd,buf)     Substitute if-then-else strings in variable         */
/*                      definition vd to if() notation. The modified        */
/*                      variable definition is put into buf. Return         */
/*                      length of buf.                                      */
            
int sub_ifs(TDAContext *ctx, char *vd,char *buf)
{
    register int level;
    register char *p,*q;
    int n,m,cnt[500];

    level = 0;
    cnt[level] = 0;
    p = vd;
    q = buf;

    while (*p) {

        if (!strncmp(p,"if",2) && (islower((int)*(p + 2)) == 0 || v_search(ctx, p + 2,&n,&m) >= 0)) {
            *q++ = 'i';
            *q++ = 'f';
            *q++ = '(';
            p += 2;
            cnt[level] += 1;
            level++;
            cnt[level] = 0;
        }
        else if (!strncmp(p,"then",4)) {
            *q++ = ',';
            p += 4;
        }
        else if (!strncmp(p,"else",4)) {
            while (cnt[level] > 0) {
                *q++ = ')';
                cnt[level] -= 1;
            }
            if (level > 0)
                level--;

            *q++ = ',';
            p += 4;
        }
        else   
            *q++ = *p++;
    }
    while (cnt[level] > 0) {
        *q++ = ')';
        cnt[level] -= 1;
    }
    *q = '\0';
    return ((int)(strlen(buf)));
}

/* ------------------------------------------------------------------------ */
/*  var_err(n,vd)     print syntax error for variable vd.                   */

void var_err(TDAContext *ctx, int n,char *vd)
{
    printf1(ctx, "Syntax error (%d): %s\n",n,vd);
}

/* ------------------------------------------------------------------------ */
/*  alloc_vdat(idx,opt)     If opt != 0  allocate memory for VDPtr and      */
/*                          NOCMaxA cases, beginning with variables at idx. */
/*                          Otherwise free the memory.                      */
/*                          Return 0 if OK, -1 if error.                    */
/*                                                                          */
/*                          If memory is successfully allocated, set        */ 
/*                          set VAlloc[j] = 2.                              */
/*                                                                          */
/*  VDPtr[i]        pointer to data field, allocated according to NOCMaxA   */
/*                  and VSLen.                                              */

int alloc_vdat(TDAContext *ctx, int idx,int opt)
{
    register int j,jj,l;
    int err = 0;

    jj = -1;
    if (opt) {
        j = idx;
        while (j >= 0) {
            if (ctx->VAlloc[j] != 1)
                gerr_exit(ctx, 48);

            l = get_slen(ctx, j,ctx->NOCMaxA);
            if (!(ctx->VDPtr[j] = (char *)calloc((size_t)(l),sizeof(char)))) {
                jj = j;         
                err = -1;
                goto AVErr; 
            }
            memrq(ctx, l,1);                    
            ctx->VAlloc[j] = 2;
            j = ctx->VNxt[j];
        }
        return(0);
    }

AVErr:
    j = idx;

    while (j >= 0) {

        if (err && j == jj)
            break;

        if (ctx->VAlloc[j] != 2)
            gerr_exit(ctx, 49);

        l = get_slen(ctx, j,ctx->NOCMaxA);
        free(ctx->VDPtr[j]);
        memrq(ctx, -l,1);
        ctx->VAlloc[j] = 1;
        j = ctx->VNxt[j];
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  alloc_vsdat(n,ivar)     Allocate memory for n variables with variable   */
/*                          number ivar[0],...,ivar[n-1].                   */
/*                          Return 0 if OK, -1 if error.                    */
/*                                                                          */
/*                          If memory is successfully allocated, set        */ 
/*                          set VAlloc[.] = 2.                              */
/*                                                                          */
/*  VDPtr[i]        pointer to data field, allocated according to NOCMaxA   */
/*                  and VSLen.                                              */

int alloc_vsdat(TDAContext *ctx, int n,short *ivar)
{
    register int i,j,ii,l;

    for (i = 0; i < n; ++i) {

        j = ivar[i];

        if (ctx->VAlloc[j] != 1)
            gerr_exit(ctx, 48);

        l = get_slen(ctx, j,ctx->NOCMaxA);

        if (!(ctx->VDPtr[j] = (char *)calloc((size_t)(l),sizeof(char)))) {

            for (ii = 0; ii < i; ++ii) {
                j = ivar[ii];

                l = get_slen(ctx, j,ctx->NOCMaxA);
                free(ctx->VDPtr[j]);
                memrq(ctx, -l,1);
                ctx->VAlloc[j] = 1;
            }
            return(-1);
        }
        memrq(ctx, l,1);                    
        ctx->VAlloc[j] = 2;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  get_slen(j,noc)     return required bytes for variable j and            */
/*                      number of cases noc.                                */

int get_slen(TDAContext *ctx, int j,int noc)
{
    int l;

    if (ctx->VTyp[j] == 2)
        l = 8;
    else if (ctx->VSLen[j] == 5)
        l = noc * 4;
    else if (ctx->VSLen[j] > 0)
        l = noc * ctx->VSLen[j];
    else if (ctx->VSLen[j] == 0)
        l = noc / 8 + 1;
    else
        l = -ctx->VSLen[j] * noc;
    return(l);
}

/* ------------------------------------------------------------------------ */
/*  prn_var(idx)    print variable definitions, begin with index idx.       */

void prn_var(TDAContext *ctx, int idx)
{
    register int i,j,cnt,k = 0;
    register char *p;
  
    if (ctx->VLabelLen > 0)
        k = 1;
    printf1(ctx, "Idx ");
    prn_hvar(ctx);
    prn_hlabel(ctx);
    printf1(ctx, " T   S  PFmt  Definition\n");
    prnchar(ctx, '-',29 + ctx->VNameLen + ctx->VLabelLen + k,1);
        
    if (idx < 0)
        i = ctx->VIFirst;
    else
        i = idx;

    j = 0;
    while (i >= 0) {              

        printf1(ctx, "%3d ",++j);
        prn_vname(ctx, i);
        prn_vlabel(ctx, i);
        printf1(ctx, "%2d %3d %3d.%-2d ",
                     ctx->VTyp[i],iabs(ctx, ctx->VSLen[i]),(int)ctx->VPFmt1[i],(int)ctx->VPFmt2[i]);
        p = ctx->VDef[i];
        cnt = 0;
        while (*p) {
            if (!strncmp(p,"if(",3)) {
                if (cnt) {
                    printf1(ctx, "\n ");
                    prnchar(ctx, ' ',16 + ctx->VNameLen + ctx->VLabelLen + k,0);
                }
                cnt++;
            }
            printf1(ctx, "%c",*p++);
        }
        printf1(ctx, "\n");
        i = ctx->VNxt[i];
    }
}

/* ------------------------------------------------------------------------ */
/*  prn_vname(i)    print name of variable i.                               */

void prn_vname(TDAContext *ctx, int i)
{
    printf1(ctx, "%s ",ctx->VName[i]);
    prnchar(ctx, ' ',ctx->VNameLen - (int)strlen(ctx->VName[i]),0);
}

/* ------------------------------------------------------------------------ */
/*  prn_vlabel(i)   print label of variable i.                              */

void prn_vlabel(TDAContext *ctx, int i)
{
    if (ctx->VLabelLen > 0) {
        if (ctx->VLabel[i] != NULL) {
            printf1(ctx, "%s ",ctx->VLabel[i]);
            prnchar(ctx, ' ',ctx->VLabelLen - (int)strlen(ctx->VLabel[i]),0);
        }
        else
            prnchar(ctx, ' ',ctx->VLabelLen + 1,0);
    }
}

/* ------------------------------------------------------------------------ */
/*  prn_hvar()   print variable string.                                     */

void prn_hvar(TDAContext *ctx)
{
    printf1(ctx, "Variable ");
    prnchar(ctx, ' ',ctx->VNameLen - 8,0);
}

/* ------------------------------------------------------------------------ */
/*  prn_hlabel()   print label string.                                      */

void prn_hlabel(TDAContext *ctx)
{
    if (ctx->VLabelLen > 0) {
        printf1(ctx, "Label  ");
        prnchar(ctx, ' ',ctx->VLabelLen - 6,0);
    }
}

/* ------------------------------------------------------------------------ */
/*  clear_dm()      clear all variables and data matrix.                    */

void clear_dm(TDAContext *ctx)
{
    register int i;

    nl_free(ctx, MaxNL);     /* clear all namelists */

    i = ctx->VIFirst;
    while (i >= 0) {
        clear_var(ctx, i);
        i = ctx->VIFirst;
    }
    svb_alloc(ctx, 0);
    sdnvar_close(ctx);
}

/* ------------------------------------------------------------------------ */
/*  clear_var(idx)      clear variable with index idx.                      */

void clear_var(TDAContext *ctx, int idx)
{
    register int l;

    /* first free data space */

    if (ctx->VAlloc[idx] == 2) {

        l = get_slen(ctx, idx,ctx->NOCMaxA);
        free(ctx->VDPtr[idx]);
        memrq(ctx, -l,1);
        ctx->VAlloc[idx] = 1;
    }
    free_var(ctx, idx);
    clear_vidx(ctx, idx);    /* sets VAlloc[idx] = 0 */
}

/* ------------------------------------------------------------------------ */
/*  clear_avar(idx)     clear all variables beginning with index idx.       */

void clear_avar(TDAContext *ctx, int idx)
{
    register int i,fnd;

    if (idx < 0)
        return;

    fnd = 0;
    i = ctx->VIFirst;
    while (i >= 0) {
        if (i == idx) {
            fnd = 1;
            break;
        }
        i = ctx->VNxt[i];
    }
    if (fnd == 0)
        return;

    i = ctx->VILast;
    while (i >= 0) {
        clear_var(ctx, i);
        if (i == idx)
            break;
        i = ctx->VILast;
    }
}

/* ------------------------------------------------------------------------ */
/*  free_var(i)     free data structures for variable i.                    */

void free_var(TDAContext *ctx, int i)
{
    if (i < 0)  
        gerr_exit(ctx, 17);
    if (ctx->VAlloc[i] != 1)
        gerr_exit(ctx, 50);
       
    memrq(ctx,(int)(-strlen(ctx->VName[i]) - 1),1);
    free(ctx->VName[i]);
    memrq(ctx,(int)(-strlen(ctx->VDef[i]) - 1),1);
    free(ctx->VDef[i]);
    if (ctx->VLabel[i] != NULL) {
        memrq(ctx,(int)(-strlen(ctx->VLabel[i]) - 1),1);
        free(ctx->VLabel[i]);
        ctx->VLabel[i] = NULL;
    }
    if (ctx->VESCnt[i] > 0) {
        free((char *)ctx->VESTyp[i]);
        memrq(ctx, -ctx->VESCnt[i],sizeof(int));
        free((char *)ctx->VESVal[i]);
        memrq(ctx, -ctx->VESCnt[i],sizeof(double));
        ctx->VESCnt[i] = 0;
    }
    ctx->VTypA[i] = ctx->VTyp[i] = 0;
}

/* ------------------------------------------------------------------------ */
/*  nlist()     set up a new namelist. nlist( Name = V1,... )               */
/*              Return 0 if OK, -1 if error.                                */

int nlist(TDAContext *ctx)
{
    register int i,j,jj;
    int err,l,n,vn,nb,r;          
    char *p,*q,*s,vname[VNLMax + 1];

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    s = ctx->CmdBuf + 6;
    l = get_vnlen(ctx, s);
    if (l < 1 || l > VNLMax) {
        printf1(ctx, "Error: no valid name for namelist.\n");
        goto NLSFin;
    }
    p = s + l;
    if (*p != '=') {
        printf1(ctx, "Error: %s\n",s);
        p_err(ctx, -1,1);
        goto NLSFin;
    }
    *p = '\0';
    vn = get_vidx1(ctx, s,vname);            

    if (vn >= 0 || nl_check(ctx, ctx->CmdBuf) >= 0) {
        printf1(ctx, "Error: name already used.\n");
        goto NLSFin;
    }
    j = -1;
    for (i = 0; i < MaxNL; ++i) {
        if (ctx->NLNV[i] == 0) {
            j = i;
            break;
        }
    }
    if (j < 0) {
        printf1(ctx, "Error: exceeded max number of namelists.\n");
        goto NLSFin;
    }
    strncpy(ctx->NLName[j],s,(size_t)(l));
    *p++ = '=';
    q = p;
    while (*q && *q != ')')
        q++;
    if (*q != ')') {
        printf1(ctx, "Error: %s\n",s);
        p_err(ctx, -1,1);
        goto NLSFin;
    }
    if (*(q - 1) == ',')
        q--;
    *q = '\0';

    q = get_nvi(ctx, p,&n,0,ctx->NLVIdx[j],&nb);
    if (nb) {
        p_err(ctx, -42,1);
        goto NLSFin;
    }
                             
    if (n == 0 || *q) {
        printf1(ctx, "Error: %s\n",s);
        p_err(ctx, -4,1);
        goto NLSFin;
    }
    if (n < 0) {
        printf1(ctx, "Error: reference to undefined variable(s).\n");
        goto NLSFin;
    }
    if (!(ctx->NLVIdx[j] = (short *)calloc((size_t)(n),sizeof(short)))) {  
        p_err(ctx, -2,1);
        goto NLSFin;
    }
    memrq(ctx, n,sizeof(short));
    ctx->NLNV[j] = (short)(n);
    q = get_nvi(ctx, p,&vn,n,ctx->NLVIdx[j],&nb);
    if (nb) {
        p_err(ctx, -42,1);
        goto NLSFin;
    }

    printf1(ctx, "New namelist: %s = %s",ctx->NLName[j],ctx->VName[ctx->NLVIdx[j][0]]);
    r = 0;
    for (i = 1; i < n; ++i) {
        jj = ctx->NLVIdx[j][i];
        printf1(ctx, ",%s",ctx->VName[jj]);
        if (ctx->VTyp[jj] == 1)
            r = 1;
    }
    printf1(ctx, "\n");
    if (r) {
        printf1(ctx, "Error: cannot use string variables in name lists.\n");
        goto NLSFin;
    }
    ctx->NNL++;
    err = 0;

NLSFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  nl_check(nl)    check whether namelist nl is defined. Return number     */
/*                  of namelist if defined, otherwise -1.                   */

int nl_check(TDAContext *ctx, char *nl)
{
    register int i;

    for (i = 0; i < MaxNL; ++i) {
        if (ctx->NLNV[i] > 0) {
            if (!strcmp(nl,ctx->NLName[i]))
                return(i);
        }
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  nl_free(j)      free memory for namelist j. If j == MaxNL free for      */
/*                  all name lists.                                         */

void nl_free(TDAContext *ctx, int j)
{
    register int i;

    if (j >= 0 && j < MaxNL) {
        if (ctx->NLNV[j] > 0) {
            free((char *)ctx->NLVIdx[j]);
            memrq(ctx, -ctx->NLNV[j],sizeof(short));
            ctx->NLNV[j] = 0;
            ctx->NLName[j][0] = '\0';
            ctx->NNL--;
        }
    }
    else if (j == MaxNL) {
        for (i = 0; i < MaxNL; ++i) {
            if (ctx->NLNV[i] > 0) {
                free((char *)ctx->NLVIdx[i]);
                memrq(ctx, -ctx->NLNV[i],sizeof(short));
                ctx->NLNV[i] = 0;
                ctx->NLName[i][0] = '\0';
            }
        }
        ctx->NNL = 0;
    }
}

/* ------------------------------------------------------------------------ */
/*  get_nvi(s,n,nv,vidx,nb)                                                 */
/*                                                                          */
/*                          s points to begin of a list of variable names.  */
/*                          Return pointer to next character after last     */
/*                          variable name. Return number of variables in    */  
/*                          n. If nv > 0 save up to nv var numbers in vidx  */
/*                          n = -1 if undefined var names. return number    */
/*                          of variable names in brackets in nb.            */
/*                                                                          */
/*  variables types are counted in global array VLTyp[].                    */

char *get_nvi(TDAContext *ctx, char *s,int *n,int nv,short *vidx,int *nb)
{
    register int i,j,k;
    register char c,*p,*q;
    int vn,vn1,l,nl,jj,bflag;
    char vname[VNLMax + 1],vname1[VNLMax + 1];

    for (i = 1; i <= 10; ++i)
        ctx->VLTyp[i] = 0;

    *nb = 0;
    bflag = i = 0;
    q = p = s;
    while (*p) {
        if (*p == '(') {
            if (bflag == 0) {
                bflag = 1;
                p++;
            }
            else {
                i = -1;
                break;
            }
        }
        l = get_vnlen(ctx, p);
        if (l <= 0)
            break;
        if (l > VNLMax) {
            i = -1;
            break;
        }
        q = p + l;
        c = *q;
        *q = '\0';
        vn = get_vidx(ctx, p,vname);
                 
        if (vn < 0) {
            j = nl_check(ctx, p);
            *q = c;
            if (j < 0) {
                i = -1;
                break;
            }
            nl = ctx->NLNV[j];   
            for (k = 0; k < nl; ++k) {
                vn = ctx->NLVIdx[j][k];
                if (i < nv) {
                    if (bflag == 0)
                        vidx[i] = (short)(vn);               
                    else
                        vidx[i] = (short)(-vn - 1);           
                }
                ctx->VLTyp[(int)ctx->VTyp[vn]] += 1;
                i++;
                if (bflag)
                    *nb += 1;
            }
        }
        else {
            *q = c;
            if (i < nv) {
                if (bflag == 0)
                    vidx[i] = (short)(vn);
                else
                    vidx[i] = (short)(-vn - 1);
            }
            ctx->VLTyp[(int)ctx->VTyp[vn]] += 1;
            i++;
            if (bflag)
                *nb += 1;
        }
        if (*q == ')') {
            if (bflag) {
                bflag = 0;
                q++;
            }
        }
        if (*q != ',') {
            p = q;
            break;
        }
        q++;
        p = q;
        if (*q == ',') {
            p++;
            l = get_vnlen(ctx, p);
            if (l <= 0 || l > VNLMax) {
                i = -1;
                break;
            }
            q = p + l;
            c = *q;
            *q = '\0';
            vn1 = get_vidx(ctx, p,vname1);
                     
            if (vn1 < 0) {
                jj = nl_check(ctx, p);
                *q = c;
                if (jj < 0) {
                    i = -1;
                    break;
                }
                nl = ctx->NLNV[jj];   
                vn1 = ctx->NLVIdx[jj][0];
                strcpy(vname1,ctx->VName[vn1]);
            }
            else {
                *q = c;
                nl = 1;
            }
            j = ctx->VNxt[vn];
            while (j >= 0 && j <= vn1) {
                if (!strcmp(vname1,ctx->VName[j]))  
                    break;    
                if (i < nv) {
                    if (bflag == 0)
                        vidx[i] = (short)(j);
                    else
                        vidx[i] = (short)(-j - 1);
                }
                ctx->VLTyp[(int)ctx->VTyp[j]] += 1;
                i++;
                if (bflag)
                    *nb += 1;
                j = ctx->VNxt[j];
            }
            if (i < nv) {
                if (bflag == 0)
                    vidx[i] = (short)(vn1);
                else
                    vidx[i] = (short)(-vn1 - 1);
            }
            ctx->VLTyp[(int)ctx->VTyp[vn1]] += 1;
            i++;
            if (bflag)
                *nb += 1;
            for (k = 1; k < nl; ++k) {
                vn1 = ctx->NLVIdx[jj][k];
                if (i < nv) {
                    if (bflag == 0)
                        vidx[i] = (short)(vn1);             
                    else
                        vidx[i] = (short)(-vn1 - 1);         
                }
                ctx->VLTyp[(int)ctx->VTyp[vn1]] += 1;
                i++;
                if (bflag)
                    *nb += 1;
            }
            if (bflag && *q == ')') {
                bflag = 0;
                q++;
            }
            if (*q != ',') {
                p = q;
                break;
            }
            q++;
            if (*q == ',') {
                i = -1;
                break;
            }           
            p = q;
        }
    }
    *n = i;
    if (*(p - 1) == '(')
        p--;
    if (*(p - 1) == ',')
        p--;
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  get_nvia(s,n,opt,nb)    s points to begin of a list of variable names.  */
/*                          Return pointer to next character after last     */
/*                          variable name. Return number of variables in    */  
/*                          n. If n > 0 create list of variable numbers in  */
/*                          VLNV, VLVIdx[]. If opt print error message.     */
/*                          Return number of variable names in brackets     */
/*                          in nb.                                          */
/*                                                                          */
/*  variables types are counted in global array VLTyp[].                    */

char *get_nvia(TDAContext *ctx, char *s,int *n,int opt,int *nb)
{
    char *p;

    p = get_nvi(ctx, s,n,0,ctx->VLVIdx,nb);
    if (*n < 1) {
        if (opt)
            printf1(ctx, "Syntax error or undefined variables.\n");
        return(p);
    }
    if (alloc_vl(ctx, *n)) {
        *n = -2;
        return(p);
    }
    p = get_nvi(ctx, s,n,ctx->VLNV,ctx->VLVIdx,nb);
    if (*n < 1) {
        if (opt)
            printf1(ctx, "Syntax error or undefined variables.\n");
        alloc_vl(ctx, 0);
    }
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  alloc_vl(n)     If n > 0 allocate VLVIdx[], otherwise free.             */
/*                  Return 0 if OK, -1 if error.                            */

int alloc_vl(TDAContext *ctx, int n)
{
    if (ctx->VLNV > 0) {
        free((char *)ctx->VLVIdx);
        memrq(ctx, -ctx->VLNV,sizeof(short));
        ctx->VLNV = 0;
    }
    if (n > 0) {
        if (!(ctx->VLVIdx = (short *)calloc((size_t)(n),sizeof(short)))) { 
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(short));
        ctx->VLNV = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  recode()    recode variables.                                           */
/*                                                                          */
/*              recode(                                                     */
/*                  dblock=...,         optional block mode                 */
/*                  VName[fmt]=expression,                                  */
/*              );                                                          */
/*                                                                          */
/*              Note: if tsel is active it is not deactivated. Recoding     */
/*              is only done for the currently selected cases, and a        */
/*              warning message is given.                                   */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int recode(TDAContext *ctx)
{
    register int i,j,ii;
    int err,l,m,n,vb,vn,w1,w2,wflg,iflg,nv,nc,n2,n3;
    register char c,*p,*q,*q1,*s;
    char vname[VNLMax + 1];
    double id,lastid;

    m = 0;
    vb = err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    if (alloc_acx(ctx, ctx->NOC))
        goto RECFin;

    p = ctx->CmdBuf + 6;
    if (*p++ != '(') {
        p_err(ctx, -1,1);
        goto RECFin;
    }
    while (*p) {
        if (!strncmp(p,"dblock=",7)) {
            q = p + 7;
            l = get_vnlen(ctx, q);
            if (l < 1 || l > VNLMax) {
                printf1(ctx, "Error: %s.\n",p);
                goto RECFin;
            }
            vb = get_vidx1(ctx, q,vname);            
            if (vb < 0) {
                printf1(ctx, "Undefined variable: %s.\n",p);
                goto RECFin;
            }   
            if (ctx->VTyp[vb] == 1) {
                printf1(ctx, "Error: cannot use string variables for dblock.\n");
                goto RECFin;
            }
            printf1(ctx, "Block mode with variable: %s\n",vname);
            p = q + l;
            if (*p == ',')
                p++;
        }
        else {
            l = get_vnlen(ctx, p);
            if (l < 1 || l > VNLMax) {
                printf1(ctx, "Error (in variable name?): %s.\n",p);
                goto RECFin;
            }
            q = p + l;
            wflg = 0;
            if (*q == '[') {
                if (sscanf(q,"[%d.%d]",&w1,&w2) != 2) {
                    p_err(ctx, -1,1);
                    goto RECFin;
                }
                wflg = 1;
                q = skip_dbl(ctx, q + 1);
                if (*q++ != ']') {
                    p_err(ctx, -1,1);
                    goto RECFin;
                }
            }
            if (*q++ != '=') {
                p_err(ctx, -1,1);
                goto RECFin;
            }
            vn = get_vidx1(ctx, p,vname);            

            if (vn < 0) {
                printf1(ctx, "Undefined variable: %s.\n",p);
                goto RECFin;
            }   
            q1 = skip_expr(ctx, q);
            c = *q1;
            *q1 = '\0';

            if (ctx->VTyp[vn] == 1 || ctx->VTyp[vn] == 2) {
                printf1(ctx, "Recoding: %s\n",p);
                printf1(ctx, "Error: cannot recode numerical constants and strings.\n");
                goto RECFin;
            }
            l = (int)(strlen(q));      /* check variable's definition */
            if (l == 0) {
                printf1(ctx, "Syntax error: %s\n",p);
                goto RECFin;
            }
            s = q;              /* check for if-then-else expressions */
            iflg = 0;
            while (*s) {
                if (!strncmp(s,"then",4)) {
                    iflg = 1;
                    break;
                }
                s++;
            }
            if (iflg) {
                if (alloc_acc(ctx, 2 * l + 100))
                    goto RECFin;
                    
                l = sub_ifs(ctx, q,ctx->AcC);     /* change to if() operators */
                q = ctx->AcC;
            }
            /* account for the old definition's length before freeing
               it -- the original order read the string after free */
            memrq(ctx,(int)(-strlen(ctx->VDef[vn]) - 1),1);
            free(ctx->VDef[vn]);

            if (!(ctx->VDef[vn] = (char *)calloc((size_t)(l + 1),sizeof(char)))) {
                p_err(ctx, -2,1);
                goto RECFin;
            }
            memrq(ctx, l + 1,sizeof(char));
            strcpy(ctx->VDef[vn],q);
                                  
            if ((n = v_parse(ctx, q,0)) < 0 || ctx->ESCnt <= 0) {
                printf1(ctx, "Syntax error (%d) in: %s\n",n,p);
                if (n < 0)
                    prn_emsg1(ctx, n);
                goto RECFin;
            }
            printf1(ctx, "Recoding: ");
            prn_vname(ctx, vn);
            printf1(ctx, " = %s\n",ctx->VDef[vn]);

            /* check type */

            check_expr(ctx, ctx->ESCnt,ctx->ESTyp,&nv,&nc,&n2,&n3,0);
            if (nc > 0) {
                p_err(ctx, -41,1); 
                goto RECFin;
            }
            if (n3 > 0) {
                p_err(ctx, -40,1);  
                goto RECFin;
            }
            if (n2 > 0)
                ctx->VTyp[vn] = 4;
            else ctx->VTyp[vn] = 3;

            *q1 = c,
             
            ii = 0;
            if (vb < 0)  
                i = ctx->NOC;
            else {
                i = 0;
                lastid = get_data(ctx, vb,i);
            }
            while (i < ctx->NOC) {
                id = get_data(ctx, vb,i);
                if (fabs(id - lastid) > ctx->EPSI1) {

                    n = v_eval2(ctx, -1,ctx->ESCnt,ctx->ESTyp,ctx->ESVal,ctx->ESIdx,ctx->AcX,ii,i,0,0);
                    if (n) {
                        printf1(ctx, "Error: can't evaluate expression.\n");
                        prn_emsg2(ctx, n);
                        goto RECFin;
                    }
                    for (j = ii; j < i; ++j)  
                        put_data(ctx, ctx->AcX[j - ii],vn,j);
 
                    ii = i;
                    lastid = id;
                }
                i++;
            }
            n = v_eval2(ctx, -1,ctx->ESCnt,ctx->ESTyp,ctx->ESVal,ctx->ESIdx,ctx->AcX,ii,i,0,0);
            if (n) {
                printf1(ctx, "Error: can't evaluate expression.\n");
                prn_emsg2(ctx, n);
                goto RECFin;
            }
            for (j = ii; j < i; ++j)     
                put_data(ctx, ctx->AcX[j - ii],vn,j);
  
            if (wflg) {
                makefmt(ctx, &w1,&w2,ctx->VPFmtS[vn],VPFmtSLen,0,' ',0);
                ctx->VPFmt1[vn] = (short)w1;
                ctx->VPFmt2[vn] = (short)w2;
            }
            m++;
            if (*q1 == ')')
                break;
            if (*q1++ != ',') {
                p_err(ctx, -1,1);
                goto RECFin;
            }
            p = q1;
        }
        if (*p == ')')
            break;
    }
    if (*p != ')') {
        p_err(ctx, -1,1);
        goto RECFin;
    }
    if (m == 0)
        printf1(ctx, "Nothing done.\n");
    else if (ctx->TSelFlg && ctx->NOC < ctx->NOCDM)
        printf1(ctx, "Warning: used only a subset of data matrix cases.\n");
    err = 0;

RECFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ndvar()     create dummy variables.                                     */
/*              Return 0 if OK, -1 if error.                                */

int ndvar(TDAContext *ctx)
{
    register int i,j;
    int err,l,vn,vx,r,n,idx,idxn,nd,nn,nb;
    short vv[2];
    register char *p,*q,*s;  
    char vname[VNLMax + 1],vname1[VNLMax + 1],vdef[500];
    double tmp;
               
    idxn = idx = get_nidx(ctx);      /* index to first new variable */
    if (idx < 0) {
        p_err(ctx, -5,1);
        return(-1); 
    }
    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Creating dummy variables. ");
    prn_mem(ctx);
    tsel_off(ctx, 1);        /* turn off temporary case selection */

    p = ctx->CmdBuf + 5;
    if (*p++ != '(' || !*p) {
        p_err(ctx, -1,1);
        goto NDVFin;
    }

    /* allocate memory */
   
    if (alloc_ack(ctx, ctx->NOC + 1))     /* value */
        goto NDVFin;

    if (alloc_acn(ctx, ctx->NOC + 1))     /* freq */
        goto NDVFin;

    if (alloc_aci(ctx, ctx->NOC + 1))     /* bf pointer */
        goto NDVFin;

    if (alloc_acr(ctx, ctx->NOC + 1))                         
        goto NDVFin;

    if (alloc_acs(ctx, ctx->NOC + 1))                         
        goto NDVFin;

    while (*p) {

        idxn = get_nidx(ctx);      /* index to first new variable */
        if (idxn < 0) {
            p_err(ctx, -5,1);
            goto NDVFin;
        }
        s = p;
        l = get_vnlen(ctx, p);
        if (l < 1 || l > VNLMax) {
            printf1(ctx, "Error (in variable name?): %s.\n",p);
            goto NDVFin;
        }
        vn = get_vidx1(ctx, p,vname);            
        if (vn >= 0) {
            printf1(ctx, "Error: %s already used.\n",ctx->VName[vn]);
            goto NDVFin;
        }
        q = p + l;
        if (!strncmp(q,"=dmul(",6)) {
            printf1(ctx, "\nCreating interaction effects: %s = dmul(...)\n",vname);
           
            q += 6;
            p = q;
            while (*p && *p != ',')
                p++;
            if (*p != ',') {
                err = -2;
                goto NDVFin;
            }
            *p = '\0';
            q = get_nvia(ctx, q,&n,1,&nb);
            if (nb) {
                p_err(ctx, -42,1);
                goto NDVFin;
            }
            *p++ = ',';
            if (n < 1 || ctx->VLNV < 1)
                goto NDVFin;

            nn = ctx->VLNV;
            if (alloc_acns(ctx, nn))                         
                goto NDVFin;

            for (j = 0; j < ctx->VLNV; ++j)  
                ctx->AcNS[j] = ctx->VLVIdx[j];
             
            p = get_nvia(ctx, p,&n,1,&nb);
            if (n < 1 || ctx->VLNV < 1)
                goto NDVFin;

            if (nb) {
                p_err(ctx, -42,1);
                goto NDVFin;
            }

            if (ndivar(ctx, vname,nn,ctx->AcNS,ctx->VLNV,ctx->VLVIdx))
                goto NDVFin;
           
            if (*p++ != ')') {
                err = -2;
                goto NDVFin;
            }
            goto NDVCONT;
        }
        else if (strncmp(q,"=dum(",5)) {
            err = -2;
            goto NDVFin;
        }
        q += 5;
        vx = get_vidx1(ctx, q,vname1);            
        if (vx < 0) {
            printf1(ctx, "Error: %s\n",p);
            printf1(ctx, "Undefined variable on right-hand side.\n");
            goto NDVFin;
        }
        q += strlen(ctx->VName[vx]);
        p = q;
        if (*p++ != ')') {
            err = -2;                         
            goto NDVFin;
        }

        /* get frequency distribution */

        vv[0] = (short)(vx);
        r = cfreq1(ctx, 1,vv,ctx->NOC,ctx->AcK,ctx->AcN,ctx->AcI,-1,&tmp);
        if (r <= 0) {
            printf1(ctx, "Error: cannot create frequency distribution for %s\n",ctx->VName[vx]);
            goto NDVFin;
        }
        nd = 0;
        i = 1;
        while (i <= r) {
            j = ctx->AcI[i];
            n = ctx->AcK[j];
            q = ctx->VName[vx];
            if (n >= 0) {
                snprintf(vdef,sizeof(vdef),"%s%d",vname,n);    
                ctx->AcR[nd] += ctx->AcN[j];
                ctx->AcS[nd++] = n;
            }
            else {                    
                snprintf(vdef,sizeof(vdef),"%sM",vname);    
                ctx->AcR[nd] += ctx->AcN[j];
                ctx->AcS[nd++] = -1;
            }    
            l = (int)(strlen(vdef));           
            if (l > VNLMax || get_vidx1(ctx, vdef,vname1) >= 0) {
                printf1(ctx, "Error: cannot create variable: %s\n",vdef);
                if (l > VNLMax)
                    printf1(ctx, "Exceeded max length of variable names.\n");
                else
                    printf1(ctx, "Variable name already used.\n");
                goto NDVFin;
            }
            snprintf(vdef,sizeof(vdef),"%s<0>[2.0]=%s[%d",vname1,q,n);    
            while (i < r && (n = ctx->AcK[ctx->AcI[i + 1]]) < 0) {
                ctx->AcR[nd - 1] += ctx->AcN[ctx->AcI[i + 1]];
                l = (int)(strlen(vdef));
                if (l > 450) {
                    printf1(ctx, "Error: exceeded max length of expression.\n");
                    goto NDVFin;
                }
                q = vdef + l;               
                snprintf(q,sizeof(vdef) - (size_t)(q - vdef),",%d",n);
                i++;
            }
            q = vdef + strlen(vdef);
            snprintf(q,sizeof(vdef) - (size_t)(q - vdef),"]");

            if (save_var(ctx, vdef,0)) {       /* save variable definition */
                printf1(ctx, "Error: cannot create: %s\n",vdef);
                goto NDVFin;
            }
            i++;
        }
        if (alloc_vdat(ctx, idxn,1)) {
            printf1(ctx, "Error: insufficient memory for new variables.\n");
            goto NDVFin;
        }
        for (i = 0; i < ctx->NOC; ++i) {
            n = (int)get_data(ctx, vx,i);
            l = idxn;
            if (n >= 0) {
                for (j = 0; j < nd; ++j) {
                    if (n == ctx->AcS[j]) {
                        l += j;
                        break;
                    }
                }
            }
            put_data(ctx, 1.0,l,i);
        }
        printf1(ctx, "\nVariable: ");
        prn_vname(ctx, vx);
        printf1(ctx, "\nValue        Cases   Pcnt  Dummy variable\n");
        prnchar(ctx, '-',41,1);

        if (alloc_acns(ctx, nd))                         
            goto NDVFin;

        nn = 0;
        for (j = 0; j < nd; ++j) {
            if (ctx->AcS[j] < 0)
                printf1(ctx, "negative  ");
            else {
                printf1(ctx, "%8d  ",ctx->AcS[j]);
                ctx->AcNS[nn++] = (short)(idxn + j);
            }

            printf1(ctx, "%8d %6.2lf  ",ctx->AcR[j],100.0 * (double)ctx->AcR[j] / (double)ctx->NOC);
            prn_vname(ctx, idxn + j);
            printf1(ctx, "= %s\n",ctx->VDef[idxn + j]);
        }
        if (nn > 0) {               /* create namelist */

            if (c_nlist(ctx, vname,nn,ctx->AcNS))
                goto NDVFin;
        }
        for (j = 0; j < nd; ++j)
            ctx->AcR[j] = 0;

NDVCONT:
        if (*p != ',')
            break;
        if (*++p == ')')
            break;
    }
    err = 0;

NDVFin:
    if (err) { 
        if (err == -2)
            printf1(ctx, "Syntax error: %s\n",s);
        clear_avar(ctx, idxn);
    }
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ndivar(vname,n1,idx1,n2,idx2)   Create interaction effects.             */
/*                                  Return 0 if OK, -1 if error.            */

int ndivar(TDAContext *ctx, char *vname,int n1,short *idx1,int n2,short *idx2)
{
    register int i,j,ii;
    int err,l,idx,idxn,ni,nj,nn,nv;
    char *p,vdef[500],vname1[VNLMax + 1];

    idx = get_nidx(ctx);      /* index to first new variable */
    if (idx < 0) {
        p_err(ctx, -5,1);
        return(-1); 
    }
    if (alloc_acms(ctx, n1 * n2 + 1))                         
        return(-1);  

    err = -1;

    printf1(ctx, "\nCases    Pcnt   Definition\n");
    prnchar(ctx, '-',43,1);
    nv = 0;
    for (i = 0; i < n1; ++i) {

        for (j = 0; j < n2; ++j) {

            idxn = get_nidx(ctx);      /* index to first new variable */
            if (idxn < 0) {
                p_err(ctx, -5,1);
                goto NDVIFin;
            }
            snprintf(vdef,sizeof(vdef),"%s_%s",ctx->VName[idx1[i]],ctx->VName[idx2[j]]);

            l = (int)(strlen(vdef));           
            if (l > VNLMax || get_vidx1(ctx, vdef,vname1) >= 0) {
                printf1(ctx, "Error: cannot create variable: %s\n",vdef);
                if (l > VNLMax)
                    printf1(ctx, "Exceeded max length of variable names.\n");
                else
                    printf1(ctx, "Variable name already used.\n");
                goto NDVIFin;
            }
            p = vdef + l;
            snprintf(p,sizeof(vdef) - (size_t)(p - vdef),"<0>[2.0]=%s&%s",ctx->VDef[idx1[i]],ctx->VDef[idx2[j]]);

            if (save_var(ctx, vdef,0)) {       /* save variable definition */
                printf1(ctx, "Error: cannot create: %s\n",vdef);
                goto NDVIFin;
            }
            if (alloc_vdat(ctx, idxn,1)) {
                printf1(ctx, "Error: insufficient memory for new variables.\n");
                goto NDVIFin;
            }
            nn = 0;
            for (ii = 0; ii < ctx->NOC; ++ii) {
                ni = (int)get_data(ctx, idx1[i],ii);
                if (ni) {
                    nj = (int)get_data(ctx, idx2[j],ii);
                    if (nj) {
                        put_data(ctx, 1.0,idxn,ii);
                        nn++;
                    }
                }
            }
            printf1(ctx, "%6d %6.2lf   ",nn,100.0 * (double)nn / (double)ctx->NOC);
            printf1(ctx, "%s",ctx->VName[idxn]);
            prnchar(ctx, ' ',(int)(16 - strlen(ctx->VName[idxn])),0);
            printf1(ctx, "= %s\n",ctx->VDef[idxn]);

            ctx->AcMS[nv++] = (short)(idxn);
        }
    }
    if (c_nlist(ctx, vname,nv,ctx->AcMS))
        goto NDVIFin;
    
    err = 0;

NDVIFin:
    if (err)   
        clear_avar(ctx, idx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  c_nlist(vname,n,vidx)   Create namelist                                 */
/*                          Return 0 if OK, -1 if error.                    */

int c_nlist(TDAContext *ctx, char *vname,int n,short *vidx)
{
    register int j,k;

    k = -1;
    for (j = 0; j < MaxNL; ++j) {
        if (ctx->NLNV[j] == 0) {
            k = j;
            break;
        }
    }
    if (k < 0) {
        printf1(ctx, "\nExceeded max number of namelists.\n");
        return(-1);
    }
    else {
        strcpy(ctx->NLName[k],vname);
        if (!(ctx->NLVIdx[k] = (short *)calloc((size_t)(n),sizeof(short)))) {  
            p_err(ctx, -2,1);
            return(-1);  
        }
        memrq(ctx, n,sizeof(short));
        ctx->NLNV[k] = (short)(n);
        for (j = 0; j < n; ++j)  
            ctx->NLVIdx[k][j] = vidx[j];
        ctx->NNL++;

        printf1(ctx, "\nNew namelist: %s = %s",vname,ctx->VName[vidx[0]]);
        if (n > 1)
            printf1(ctx, ",...,%s",ctx->VName[vidx[n - 1]]);
        newline(ctx);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  get_mxvlen(n,vidx)  Return max length of variable names in vidx[]       */

int get_mxvlen(TDAContext *ctx, int n,short *vidx)
{
    register int i,l,ml;

    ml = 0;
    for (i = 0; i < n; ++i) {
        l = (int)(strlen(ctx->VName[vidx[i]]));
        if (ml < l)
            ml = l;        
    }
    return(ml);
}

/* ------------------------------------------------------------------------ */
/*  prn_vlist(n,vidx)   Print names of variables in vidx[]                  */

void prn_vlist(TDAContext *ctx, int n,short *vidx)
{
    register int i;

    for (i = 0; i < n; ++i) {
        if (i)
            printf1(ctx, ",");
        printf1(ctx, "%s",ctx->VName[vidx[i]]);
    }
}

/* ------------------------------------------------------------------------ */
/*  svb_alloc(n)    If n > 0 allocate SVBuf, otherwise free.                */
/*                  Return 0 if OK, -1 if error.                            */

int svb_alloc(TDAContext *ctx, int n)
{
    if (ctx->SVBufA > 0) {
        free((char *)ctx->SVBuf);
        memrq(ctx, -ctx->SVBufA,sizeof(char));
        ctx->SVBufA = ctx->SVBufLen = 0;
    }
    if (n > 0) {
        if (!(ctx->SVBuf = (char *)calloc((size_t)(n + 2),sizeof(char)))) { 
            p_err(ctx, -2,1);    
            return(-1);
        }
        memrq(ctx, n + 2,sizeof(char));
        ctx->SVBufA = n + 2;
        ctx->SVBufLen = n;
    }         
    return(0);
}
