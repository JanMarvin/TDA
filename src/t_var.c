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

/*  functions in t_var.c */

char get_vnchar(char c);
int check_vname(char *p);
int get_vnlen(char *p);
int get_vnl(int n,short *vinum);
int alloc_vmax(void);
int get_nidx(void);
void clear_vidx(int i);
int get_vidx(char *p,char *name);
int get_vidx1(char *p,char *name);
int get_nvidx(char *p,short *vinum);
int save_var(char *vd,int opt);
int sub_ifs(char *vd,char *buf); 
void var_err(int n,char *vd);
int alloc_vdat(int idx,int opt);
int alloc_vsdat(int n,short *ivar);
int get_slen(int j,int noc);
void prn_var(int idx);
void prn_vname(int i);
void prn_vlabel(int i);
void prn_hvar(void);
void prn_hlabel(void);
void clear_dm(void);
void clear_var(int idx);
void clear_avar(int idx);
void free_var(int i);
int nlist(void);
int nl_check(char *nl);
char *get_nvi(char *s,int *n,int nv,short *vidx,int *nb);
char *get_nvia(char *s,int *n,int opt,int *nb);
void nl_free(int j);
int alloc_vl(int n);
int recode(void);
int ndvar(void);
int ndivar(char *vname,int n1,short *idx1,int n2,short *idx2);
int c_nlist(char *vname,int n,short *vidx);
int get_mxvlen(int n,short *vidx);
void prn_vlist(int n,short *vidx);
int svb_alloc(int n);

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */

#define VPFmtSLen 10    /* max length of print format strings               */
int MaxNV = 0;          /* max number of variables                          */
int NVAR = 0;           /* number of variables                              */

short VIFirst = -1;     /* index of first variable                          */
short VILast  = -1;     /* index of last variable                           */

char *VAlloc;           /* flags for allocated data structures              */
short *VNxt;            /* pointer to next variables                        */

char **VName;           /* names of variables                               */
char **VDef;            /* definitions of variables                         */
char **VLabel;          /* labels of variables                              */

short *VSLen;           /* storage length                                   */
char *VTyp;             /* type of variable                                 */
char *VTypA;            /* 1 if data for a variable in data matrix          */
                        /* available. 2 if archive variable and data not    */  
                        /* not already created.                             */
short *VPFmt1;          /* print format of variables                        */
short *VPFmt2;
char **VPFmtS;

char **VDPtr;           /* pointer to data fields                           */

int **VESTyp;           /* parser stack for variables                       */
double **VESVal;        /* parser value for variables                       */
short *VESCnt;          /* size of parser arrays                            */

int VNameLen = 8;       /* max length of variable names                     */
int VLabelLen = 0;      /* max length of variables labels, also used to     */
                        /* flag if at least one var label present.          */

int NVArc = 0;          /* number of archive variables                      */
int NVArc1 = 0;         /* number of new archive variables, set in save_var */
short *AVIdx;           /* internal variable number                         */
int AVIdxA = 0;
short *AVOff;           /* offset of variable                               */
int AVOffA = 0;
short *AVLen;           /* length of variable                               */
int AVLenA = 0;
short *AVFmt1;          /* format of variable                               */
int AVFmt1A = 0;
short *AVFmt2;          
int AVFmt2A = 0;
int *AVMBlnk;           /* missing values: blank                            */
int AVMBlnkA = 0;
int *AVMStar;           /* missing values: star                             */
int AVMStarA = 0;
int *AVMPnt ;           /* missing values: point                            */
int AVMPntA = 0;
int *AVMGen ;           /* missing values: general                          */
int AVMGenA = 0;
int AVDFN = -1;         /* index to data file for archive variables         */

short *VStrN;           /* first column of string variable (-1)             */

int NNL = 0;                    /* number of namelists                      */
char NLName[MaxNL][VNLMax+1];   /* names                                    */
short NLNV[MaxNL];              /* number of variables                      */
short *NLVIdx[MaxNL];           /* indices of variables                     */
short VLTyp[11];                /* counting variable types                  */

int VLNV = 0;           /* number of variables in VLVIdx[]                  */
short *VLVIdx;

char *SVBuf;            /* buffer for longest string                        */
int SVBufA = 0;         /* allocated                                        */
int SVBufLen = 0;       /* current maximal string length                    */

/*--------------------------------------------------------------------------*/
/*  get_vnchar(c)   If c is valid character for variable name return c,     */
/*                  otherwise return underline.                             */

char get_vnchar(char c)
{
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

int check_vname(char *p)
{
    if ((*p >= 'A' && *p <= 'Z') || *p == '_' || *p == '@' || *p == '$')
        return(1);
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  get_vnlen(p)    get length of variable name pointed to by p.            */
/*                  return this length.                                     */

int get_vnlen(char *p)
{
    register int l;
           
    if (check_vname(p) == 0)
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

int get_vnl(int n,short *vinum)
{
    register int i,j,l;

    l = 0;
    for (i = 0; i < n; ++i) {
        j = strlen(VName[vinum[i]]);
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

int alloc_vmax(void)
{
    int i,n1,n2;

    n1 = 4 * MaxNV * sizeof(char *) + MaxNV * sizeof(short);  

    if (!(VDef   = (char **)calloc(MaxNV,sizeof(char *))) ||
        !(VName  = (char **)calloc(MaxNV,sizeof(char *))) ||
        !(VLabel = (char **)calloc(MaxNV,sizeof(char *))) || 
        !(VDPtr  = (char **)calloc(MaxNV,sizeof(char *))) || 
        !(VStrN  = (short *)calloc(MaxNV,sizeof(short )))) {  

        p_err(-2,1);    
        return(-1);
    }         
    memrq(n1,1);

    n2 = MaxNV * (3 * sizeof(char) + 5 * sizeof(short) + sizeof(char **) +
                      sizeof(int *) + sizeof(double) + sizeof(double *));

    if (!(VAlloc = (char *)calloc(MaxNV,sizeof(char))) ||
        !(VTyp   = (char *)calloc(MaxNV,sizeof(char))) ||
        !(VTypA  = (char *)calloc(MaxNV,sizeof(char))) ||
        !(VSLen  = (short *)calloc(MaxNV,sizeof(short))) ||
        !(VPFmt1 = (short *)calloc(MaxNV,sizeof(short))) ||
        !(VPFmt2 = (short *)calloc(MaxNV,sizeof(short))) ||
        !(VPFmtS = (char **)calloc(MaxNV,sizeof(char *))) ||
        !(VNxt   = (short *)calloc(MaxNV,sizeof(short))) ||  
        !(VESCnt = (short *)calloc(MaxNV,sizeof(short))) ||  
        !(VESTyp = (int **)calloc(MaxNV,sizeof(int *))) ||  
        !(AVVAL  = (double *)calloc(MaxNV,sizeof(double))) ||
        !(VESVal = (double **)calloc(MaxNV,sizeof(double *)))) {  
        p_err(-2,1);   
        return(-1);
    }         
    memrq(n2,1);

    if (alloc_est(P_ALLOC,1))
        return(-1);

    for (i = 0; i < MaxNV; ++i) {

        if (!(VPFmtS[i] = (char *)calloc(VPFmtSLen,sizeof(char)))) {
            p_err(-2,1);             
            return(-1);
        }         
        memrq(VPFmtSLen,1);

        VLabel[i] = NULL;
        VNxt[i] = -1;
        VTyp[i] = VTypA[i] = 0;
        VAlloc[i] = VESCnt[i] = 0;
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  get_nidx()  get new variable index.                                     */
/*              return index, or -1 if maximum reached.                     */

int get_nidx(void)
{
    register int i,j;

    j = -1;
    for (i = 0; i < MaxNV; ++i) {
        if (VAlloc[i] == 0) {
            j = i;
            break;
        }
    }
    return(j);
}

/*--------------------------------------------------------------------------*/
/*  clear_vidx(i)   remove variable with index i (0 <= i < MaxNV) from      */
/*                  variable list.                                          */

void clear_vidx(int i)
{
    register int j,k;

    if (i < 0 || VIFirst < 0 || NVAR <= 0)  
        gerr_exit(15);

    if (VAlloc[i] != 1)
        gerr_exit(51);
      
    j = k = VIFirst;
    while (i != j) {
        k = j;
        if (k < 0)
            gerr_exit(15);
        j = VNxt[k];
    }
    VNxt[k] = VNxt[i];
    if (VIFirst == i)
        VIFirst = VNxt[i];
    if (VILast == i) {
        if (VIFirst == -1)
            VILast = -1;
        else
            VILast = k;
    }
    VNxt[i] = -1;
    VAlloc[i] = 0;
    NVAR--;
}

/* ------------------------------------------------------------------------ */
/*  get_vidx(p,name)    p points to the begin of a variable name in a       */
/*                      string, possibly ended by: , {  < [ ( =             */
/*                      Return variable name in name[].                     */
/*                      Check whether this variable name is part of the     */
/*                      current data matrix. If this is the case return     */
/*                      internal variable number, otherwise -1.             */

int get_vidx(char *p,char *name)
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

    i = VIFirst;
    while (i >= 0) {
        if (!strcmp(name,VName[i]))
            return(i);
        i = VNxt[i];
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

int get_vidx1(char *p,char *name)
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

    i = VIFirst;
    while (i >= 0) {
        if (!strcmp(name,VName[i]))
            return(i);
        i = VNxt[i];
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

int get_nvidx(char *p,short *vinum)
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
            i = VIFirst;
            while (i >= 0) {
                if (!strcmp(name,VName[i])) {
                    j = 0;
                    break;    
                }
                i = VNxt[i];
            }
        }
        *q = c;
        if (j)          /* not found */
            return(-1);

        ib = i;
        if (cflg) {
            if (ia >= 0 && ia != ib) {

                i = VNxt[ia];
                while (i >= 0) {
                    if (!strcmp(name,VName[i]))  
                        break;    
                    vinum[k++] = i;
                    n++;
                    i = VNxt[i];
                }
            }
        }
        n++;
        vinum[k++] = ib;
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

int save_var(char *vd,int opt)
{
    register int j,l;
    register char *p,*q,*s,*pp;
    int err,err1,n,nc,nv,n2,n3,m,mm1,mm1s,mm2,mm3,blevel,ci,cii,w1,w2;
    int ni,ni1,ni2,ni3,ni3len,nibuflen;
    char *buf,*nibuf,name[VNLMax+100],vnum[60];

    err = err1 = 0;
    nibuflen = ni = ni1 = ni2 = ni3 = -1;

VDNXT:

    j = get_nidx();     /* get new variable index */

    if (j < 0) {
        printf1("Error: exceeded max number of variables (%d).\n",MaxNV);
        return(1); 
    }
    if (VESCnt[j] > 0) {
        gerr_exit(16);
    }
    if (ni == -1) {
        if (get_vidx(vd,name) >= 0) {
            printf1("Already used: %s\n",vd);
            err = -1;
            goto SVFin;
        }
        mm1s = mm1 = strlen(name);
        pp = p = vd + mm1;
         
        /* ### check for special {j,n,k} */

        if (*p == '{') {
            if (sscanf(p,"{%d,%d,%d}",&ni1,&ni2,&ni3) == 3 &&
                                 ni1 >= 0 && ni2 >= 1 && ni3 >= 1) {
                p = skip_int(p + 1);
                p = skip_int(p + 1);
                p = skip_int(p + 1);
            }
            else
                goto SVSERR; 

            if (*p++ != '}')
                goto SVSERR;
            pp = p;
            ni = ni1;
            ci = ni3;

            sprintf(vnum,"c%d",ni3);
            ni3len = strlen(vnum);

            nibuflen = strlen(p) + 1000;
            if (!(nibuf = (char *)calloc(nibuflen,sizeof(char)))) {
                err = -3;
                goto SVFin1;
            }
            memrq(nibuflen,1);
        }
    }
    else
        p = pp;

    if (ni >= 0) {
        sprintf(name + mm1s,"%d",ni);
        mm1 = strlen(name);
        if (mm1 > VNLMax) {
            err1 = 1;
            goto SVSERR;
        }
        if (get_vidx(name,name) >= 0) {
            printf1("Already used: %s\n",name);
            err = -1;
            goto SVFin;
        }
        p = pp;
    }     
    if (VNameLen < mm1)
        VNameLen = mm1;
    mm1++;
    if (!(VName[j] = (char *)calloc(mm1,sizeof(char)))) {
        err = -3;
        goto SVFin1;       
    }
    memrq(mm1,1);
    strcpy(VName[j],name);

    if (mat_getidx(VName[j],0) >= 0) {
        printf1("Variable name %s already used for a matrix.\n",VName[j]);
        err = -1;
        goto SVFin2;
    }

    if (*p == '<') {
        if (sscanf(p,"<%d>",&n) == 1 &&            
                (n == 0 || n == 1 || n == 2 || n == 4 || n == 5 || n == 8 || n < 0))
            VSLen[j] = n;
        else {
            var_err(1,vd);
            err = -1;
            goto SVFin2;
        }
        p = skip_int(++p);
        p++;
    }
    else
        VSLen[j] = 4;   /* default storage size */

    VPFmt1[j] =  0;      
    VPFmt2[j] = -1;

    if (*p == '[') {
        if (sscanf(p,"[%d.%d]",&n,&m) == 2) {
            VPFmt1[j] = n;
            VPFmt2[j] = m;
            p = skip_int(++p);
            p = skip_int(++p);
        }
        else if (sscanf(p,"[%d]",&n) == 1 && n >= 0) {
            VPFmt1[j] = n;
            VPFmt2[j] = 0;
            p = skip_int(++p);
        }
        else {
            var_err(2,vd);
            err = -1;
            goto SVFin2;
        }
        p++;

        w1 = (int)VPFmt1[j];
        w2 = (int)VPFmt2[j];
    }
    else {              /* default print format */
        w2 = 0;
        switch (VSLen[j]) {
            case 0:   w1 =  2; break;
            case 1:   w1 =  4; break;
            case 2:   w1 =  6; break;
            case 5:   w1 = 11; break;
            default:  w1 = 0;  break;
        }
    }
    makefmt(&w1,&w2,VPFmtS[j],0,' ',0);
    VPFmt1[j] = (short)w1;
    VPFmt2[j] = (short)w2;

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
            var_err(3,vd);
            err = -1;
            goto SVFin2;
        }
        mm2 = l + 1;
        if (!(VLabel[j] = (char *)calloc(mm2,sizeof(char)))) {
            err = -3;           
            VLabel[j] = NULL;
            goto SVFin2;
        }
        memrq(mm2,1);

        /****************
        q = VLabel[j];
        while (*p && *p != ')')  
            *q++ = *p++;
        *q = '\0';
        p++;
        *********/
        s = VLabel[j];
        while (p < q)
            *s++ = *p++;
        *s = '\0';
        p++;
         
        if (VLabelLen < l)
            VLabelLen = l;
    }
    if (*p++ != '=') {
        var_err(4,vd);
        err = -1;
        goto SVFin3;
    }
    l = strlen(p);      /* check variable's definition */
    mm3 = l + 1;
    if (l == 0) {
        var_err(5,vd);
        err = -1;
        goto SVFin3;
    }

    l = strlen(p);      /* check variable's definition */
    mm3 = l + 1;
    if (l == 0) {
        var_err(5,vd);
        err = -1;
        goto SVFin3;
    }
    if (ni >= 0) {      /* put var definition into nibuf */
        q = p;
        s = nibuf;

        while (*q) {
            if (sscanf(q,"c%d",&cii) == 1 && cii == ni3) {
                sprintf(s,"c%d",ci);
                s += strlen(s);
                q += ni3len;
            }
            else
                *s++ = *q++;
        }      
        *s = '\0';
        p = nibuf;
        l = strlen(p);  
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
        if (!(buf = (char *)calloc(2 * l,sizeof(char)))) {
            err = -3;
            goto SVFin3;
        }
        m = sub_ifs(p,buf);     /* change to if() operators */
        mm3 = m + 1;
        p = buf;
    }
                
    if (!(VDef[j] = (char *)calloc(mm3,sizeof(char)))) {
        err = -3;
        if (n)
            free(buf);
        goto SVFin3;
    }
    memrq(mm3,1);
    strcpy(VDef[j],p);
    if (n)
        free(buf);

    if (!strncmp(p,"spss",4)) {   /* SPSS export/sav file variable */
        VTyp[j] = 3;
        if (sscanf(p,"spss(%d)",&n) == 1 && n > 0) {   /* string variable */
            VTyp[j] = 1;
            VSLen[j] = -n;           
            VStrN[j] = 0;
            if (n > SVBufLen) {
                if (svb_alloc(n))
                    goto SVFin3;
            }
        }
        ESCnt = 1;
        ESTyp[0] = -1;
        ESVal[0] = 0.0;
    }
    else if (!strncmp(p,"stata",5)) {   /* STATA file variable */
        VTyp[j] = 3;
        if (sscanf(p,"stata(%d)",&n) == 1 && n > 0) {   /* string variable */
            VTyp[j] = 1;
            VSLen[j] = -n;           
            VStrN[j] = 0;
            if (n > SVBufLen) {
                if (svb_alloc(n))
                    goto SVFin3;
            }
        }
        ESCnt = 1;
        ESTyp[0] = -1;
        ESVal[0] = 0.0;
    }
    else if (!strncmp(p,"str(",4)) {

        if (sscanf(p,"str(%d,%d)",&n,&m) != 2 || n < 1 || m < n) {
            var_err(6,vd);
            err = -1;
            goto SVFin3;
        }
        p = skip_int(p + 4);
        p = skip_int(p + 1);
        if (*++p) {
            var_err(6,vd);
            err = -1;
            goto SVFin3;
        }
        VTyp[j] = 1;
        VSLen[j] = n - m - 1;
        VStrN[j] = n - 1;

        if (-VSLen[j] > SVBufLen) {
            if (svb_alloc(-VSLen[j]))
                goto SVFin3;
        }
        ESCnt = 1;
        ESTyp[0] = VOFFS + j;
        ESVal[0] = 0.0;
    }
    else if (*p == 'A' && *(p + 1) == ':') {  /* check for archive variables */

        p += 2;
        n = strlen(p);
        /* n = get_vnlen(p); */

        if (n < 1 || n + 2 != l) {
            printf1("Syntax error in archive variable: %s\n",vd);
            err = -1;
            goto SVFin4; 
        }
        if (VSLen[j] < 0) {
            VTyp[j] = 1;
            if (-VSLen[j] > SVBufLen) {
                if (svb_alloc(-VSLen[j]))
                    goto SVFin3;
            }
        }
        else
            VTyp[j] = 3;
        VTypA[j] = 2;
        ESCnt = 1;
        ESTyp[0] = VOFFS + j;
        ESVal[0] = 0.0;
        NVArc1++;
    }
    else if (*p == 'M' && *(p + 1) == ':') {  /* check for matrix variables */

        VTyp[j] = 3;
        ESCnt = 1;
        ESTyp[0] = VOFFS;
        ESVal[0] = 0.0;
    }
    else {

        VTyp[j] = 3;
        
        /* check whether the variable could be correctly parsed   
           and save results in a variable-specific parser stack. */          

        PREVName = VName[j];    /* special for pre operator */
        PREVNum = j;

        n = v_parse(VDef[j],0);

        PREVName = NULL;
        PREVNum = -1;
          
        if (n < 0 || ESCnt <= 0) {
            printf1("Syntax error (%d) in variable: %s\n",-n,vd);
            if (n < 0)
                prn_emsg1(n);

            err = -1;
            goto SVFin4; 
        }

        /* check for a numerical constant */

        if (ESCnt == 1 && (ESTyp[0] == 0 || ESTyp[0] == 1902)) {
            VTyp[j] = 2;
            if (opt == 0)
                VSLen[j] = 8;
        }
        else {      /* check for type 2 and 3 operators */
            
            check_expr(ESCnt,ESTyp,&nv,&nc,&n2,&n3,0);
            if (n3 > 0) {
                if (opt == 0) {
                    p_err(-40,1);
                    err = -2;
                    goto SVFin4;
                }
                else if (n2 > 0) {
                    printf1("Error: type 4 and type 5 variables not compatible.\n");
                    err = -2;
                    goto SVFin4;
                }
                else
                    VTyp[j] = 5;
            }
            else if (n2 > 0)
                VTyp[j] = 4;
        }
    }
    if (opt)
        VTyp[j] = 5;

    /* re-make print format for string variables */

    if (VTyp[j] == 1) {
        w1 = w2 = 0;
        makefmt(&w1,&w2,VPFmtS[j],0,' ',0);
        VPFmt1[j] = (short)w1;
        VPFmt2[j] = (short)w2;
    }

    /* copy parser stack for variable */

    if (!(VESTyp[j] = (int *)calloc(ESCnt,sizeof(int)))) {
        err = -3;         
        goto SVFin4;
    }
    memrq(ESCnt,sizeof(int));

    if (!(VESVal[j] = (double *)calloc(ESCnt,sizeof(double)))) {
        err = -3;         
        goto SVFin5;
    }
    memrq(ESCnt,sizeof(double));

    VESCnt[j] = (short)ESCnt;
    for (l = 0; l < ESCnt; ++l) {
        VESTyp[j][l] = ESTyp[l];
        VESVal[j][l] = ESVal[l];
    }
    VAlloc[j] = 1;
    if (VIFirst < 0)
        VILast = VIFirst = j;
    else {
        VNxt[VILast] = j;
        VILast = j;
    }
    NVAR++;
    if (ni >= 0 && ni < ni2) {
        ni++;
        ci++;                                
        goto VDNXT;
    }
    if (nibuflen > 0) {
        free(nibuf);
        memrq(-nibuflen,1);
    }
    return(0);

SVFin5:
    free((char *)VESTyp[j]);
    memrq(-ESCnt,sizeof(int));
SVFin4:
    free(VDef[j]);
    memrq(-mm3,1);
SVFin3:
    if (VLabel[j] != NULL) {
        free(VLabel[j]);
        memrq(-mm2,1);            
        VLabel[j] = NULL;
    }
SVFin2:
    free(VName[j]);
    memrq(-mm1,1);
SVFin1:
    if (err == -3)  
        printf1("Insufficient memory for data structures of variable: %s\n",vd);
       
SVFin:
    if (nibuflen > 0) {
        free(nibuf);
        memrq(-nibuflen,1);
    }
    return(err);

SVSERR:
    if (nibuflen > 0) {
        free(nibuf);
        memrq(-nibuflen,1);
    }
    printf1("Syntax error in variable: %s\n",vd);
    if (err1)   
        printf1("Exceeded max var name length.\n");
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  sub_ifs(vd,buf)     Substitute if-then-else strings in variable         */
/*                      definition vd to if() notation. The modified        */
/*                      variable definition is put into buf. Return         */
/*                      length of buf.                                      */
            
int sub_ifs(char *vd,char *buf)
{
    register int level;
    register char *p,*q;
    int n,m,cnt[500];

    level = 0;
    cnt[level] = 0;
    p = vd;
    q = buf;

    while (*p) {

        if (!strncmp(p,"if",2) && (islower((int)*(p + 2)) == 0 || v_search(p + 2,&n,&m) >= 0)) {
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
    return(strlen(buf));
}

/* ------------------------------------------------------------------------ */
/*  var_err(n,vd)     print syntax error for variable vd.                   */

void var_err(int n,char *vd)
{
    printf1("Syntax error (%d): %s\n",n,vd);
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

int alloc_vdat(int idx,int opt)
{
    register int j,jj,l;
    int err = 0;

    jj = -1;
    if (opt) {
        j = idx;
        while (j >= 0) {
            if (VAlloc[j] != 1)
                gerr_exit(48);

            l = get_slen(j,NOCMaxA);
            if (!(VDPtr[j] = (char *)calloc(l,sizeof(char)))) {
                jj = j;         
                err = -1;
                goto AVErr; 
            }
            memrq(l,1);                    
            VAlloc[j] = 2;
            j = VNxt[j];
        }
        return(0);
    }

AVErr:
    j = idx;

    while (j >= 0) {

        if (err && j == jj)
            break;

        if (VAlloc[j] != 2)
            gerr_exit(49);

        l = get_slen(j,NOCMaxA);
        free(VDPtr[j]);
        memrq(-l,1);
        VAlloc[j] = 1;
        j = VNxt[j];
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

int alloc_vsdat(int n,short *ivar)
{
    register int i,j,ii,l;

    for (i = 0; i < n; ++i) {

        j = ivar[i];

        if (VAlloc[j] != 1)
            gerr_exit(48);

        l = get_slen(j,NOCMaxA);

        if (!(VDPtr[j] = (char *)calloc(l,sizeof(char)))) {

            for (ii = 0; ii < i; ++ii) {
                j = ivar[ii];

                l = get_slen(j,NOCMaxA);
                free(VDPtr[j]);
                memrq(-l,1);
                VAlloc[j] = 1;
            }
            return(-1);
        }
        memrq(l,1);                    
        VAlloc[j] = 2;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  get_slen(j,noc)     return required bytes for variable j and            */
/*                      number of cases noc.                                */

int get_slen(int j,int noc)
{
    int l;

    if (VTyp[j] == 2)
        l = 8;
    else if (VSLen[j] == 5)
        l = noc * 4;
    else if (VSLen[j] > 0)
        l = noc * VSLen[j];
    else if (VSLen[j] == 0)
        l = noc / 8 + 1;
    else
        l = -VSLen[j] * noc;
    return(l);
}

/* ------------------------------------------------------------------------ */
/*  prn_var(idx)    print variable definitions, begin with index idx.       */

void prn_var(int idx)
{
    register int i,j,cnt,k = 0;
    register char *p;
  
    if (VLabelLen > 0)
        k = 1;
    printf1("Idx ");
    prn_hvar();
    prn_hlabel();
    printf1(" T   S  PFmt  Definition\n");
    prnchar('-',29 + VNameLen + VLabelLen + k,1);
        
    if (idx < 0)
        i = VIFirst;
    else
        i = idx;

    j = 0;
    while (i >= 0) {              

        printf1("%3d ",++j);
        prn_vname(i);
        prn_vlabel(i);
        printf1("%2d %3d %3d.%-2d ",
                     VTyp[i],iabs(VSLen[i]),(int)VPFmt1[i],(int)VPFmt2[i]);
        p = VDef[i];
        cnt = 0;
        while (*p) {
            if (!strncmp(p,"if(",3)) {
                if (cnt) {
                    printf1("\n ");
                    prnchar(' ',16 + VNameLen + VLabelLen + k,0);
                }
                cnt++;
            }
            printf1("%c",*p++);
        }
        printf1("\n");
        i = VNxt[i];
    }
}

/* ------------------------------------------------------------------------ */
/*  prn_vname(i)    print name of variable i.                               */

void prn_vname(int i)
{
    printf1("%s ",VName[i]);
    prnchar(' ',VNameLen - strlen(VName[i]),0);
}

/* ------------------------------------------------------------------------ */
/*  prn_vlabel(i)   print label of variable i.                              */

void prn_vlabel(int i)
{
    if (VLabelLen > 0) {
        if (VLabel[i] != NULL) {
            printf1("%s ",VLabel[i]);
            prnchar(' ',VLabelLen - strlen(VLabel[i]),0);
        }
        else
            prnchar(' ',VLabelLen + 1,0);
    }
}

/* ------------------------------------------------------------------------ */
/*  prn_hvar()   print variable string.                                     */

void prn_hvar(void)
{
    printf1("Variable ");
    prnchar(' ',VNameLen - 8,0);
}

/* ------------------------------------------------------------------------ */
/*  prn_hlabel()   print label string.                                      */

void prn_hlabel(void)
{
    if (VLabelLen > 0) {
        printf1("Label  ");
        prnchar(' ',VLabelLen - 6,0);
    }
}

/* ------------------------------------------------------------------------ */
/*  clear_dm()      clear all variables and data matrix.                    */

void clear_dm(void)
{
    register int i;

    nl_free(MaxNL);     /* clear all namelists */

    i = VIFirst;
    while (i >= 0) {
        clear_var(i);
        i = VIFirst;
    }
    svb_alloc(0);
    sdnvar_close();
}

/* ------------------------------------------------------------------------ */
/*  clear_var(idx)      clear variable with index idx.                      */

void clear_var(int idx)
{
    register int l;

    /* first free data space */

    if (VAlloc[idx] == 2) {

        l = get_slen(idx,NOCMaxA);
        free(VDPtr[idx]);
        memrq(-l,1);
        VAlloc[idx] = 1;
    }
    free_var(idx);
    clear_vidx(idx);    /* sets VAlloc[idx] = 0 */
}

/* ------------------------------------------------------------------------ */
/*  clear_avar(idx)     clear all variables beginning with index idx.       */

void clear_avar(int idx)
{
    register int i,fnd;

    if (idx < 0)
        return;

    fnd = 0;
    i = VIFirst;
    while (i >= 0) {
        if (i == idx) {
            fnd = 1;
            break;
        }
        i = VNxt[i];
    }
    if (fnd == 0)
        return;

    i = VILast;
    while (i >= 0) {
        clear_var(i);
        if (i == idx)
            break;
        i = VILast;
    }
}

/* ------------------------------------------------------------------------ */
/*  free_var(i)     free data structures for variable i.                    */

void free_var(int i)
{
    if (i < 0)  
        gerr_exit(17);
    if (VAlloc[i] != 1)
        gerr_exit(50);
       
    memrq(-strlen(VName[i]) - 1,1);
    free(VName[i]);
    memrq(-strlen(VDef[i]) - 1,1);
    free(VDef[i]);
    if (VLabel[i] != NULL) {
        memrq(-strlen(VLabel[i]) - 1,1);
        free(VLabel[i]);
        VLabel[i] = NULL;
    }
    if (VESCnt[i] > 0) {
        free((char *)VESTyp[i]);
        memrq(-VESCnt[i],sizeof(int));
        free((char *)VESVal[i]);
        memrq(-VESCnt[i],sizeof(double));
        VESCnt[i] = 0;
    }
    VTypA[i] = VTyp[i] = 0;
}

/* ------------------------------------------------------------------------ */
/*  nlist()     set up a new namelist. nlist( Name = V1,... )               */
/*              Return 0 if OK, -1 if error.                                */

int nlist(void)
{
    register int i,j,jj;
    int err,l,n,vn,nb,r;          
    char *p,*q,*s,vname[VNLMax + 1];

    err = -1;
    if (check_cmd(0))
        return(-1);

    s = CmdBuf + 6;
    l = get_vnlen(s);
    if (l < 1 || l > VNLMax) {
        printf1("Error: no valid name for namelist.\n");
        goto NLSFin;
    }
    p = s + l;
    if (*p != '=') {
        printf1("Error: %s\n",s);
        p_err(-1,1);
        goto NLSFin;
    }
    *p = '\0';
    vn = get_vidx1(s,vname);            

    if (vn >= 0 || nl_check(CmdBuf) >= 0) {
        printf1("Error: name already used.\n");
        goto NLSFin;
    }
    j = -1;
    for (i = 0; i < MaxNL; ++i) {
        if (NLNV[i] == 0) {
            j = i;
            break;
        }
    }
    if (j < 0) {
        printf1("Error: exceeded max number of namelists.\n");
        goto NLSFin;
    }
    strncpy(NLName[j],s,l);
    *p++ = '=';
    q = p;
    while (*q && *q != ')')
        q++;
    if (*q != ')') {
        printf1("Error: %s\n",s);
        p_err(-1,1);
        goto NLSFin;
    }
    if (*(q - 1) == ',')
        q--;
    *q = '\0';

    q = get_nvi(p,&n,0,NLVIdx[j],&nb);
    if (nb) {
        p_err(-42,1);
        goto NLSFin;
    }
                             
    if (n == 0 || *q) {
        printf1("Error: %s\n",s);
        p_err(-4,1);
        goto NLSFin;
    }
    if (n < 0) {
        printf1("Error: reference to undefined variable(s).\n");
        goto NLSFin;
    }
    if (!(NLVIdx[j] = (short *)calloc(n,sizeof(short)))) {  
        p_err(-2,1);
        goto NLSFin;
    }
    memrq(n,sizeof(short));
    NLNV[j] = n;
    q = get_nvi(p,&vn,n,NLVIdx[j],&nb);
    if (nb) {
        p_err(-42,1);
        goto NLSFin;
    }

    printf1("New namelist: %s = %s",NLName[j],VName[NLVIdx[j][0]]);
    r = 0;
    for (i = 1; i < n; ++i) {
        jj = NLVIdx[j][i];
        printf1(",%s",VName[jj]);
        if (VTyp[jj] == 1)
            r = 1;
    }
    printf1("\n");
    if (r) {
        printf1("Error: cannot use string variables in name lists.\n");
        goto NLSFin;
    }
    NNL++;
    err = 0;

NLSFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  nl_check(nl)    check whether namelist nl is defined. Return number     */
/*                  of namelist if defined, otherwise -1.                   */

int nl_check(char *nl)
{
    register int i;

    for (i = 0; i < MaxNL; ++i) {
        if (NLNV[i] > 0) {
            if (!strcmp(nl,NLName[i]))
                return(i);
        }
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  nl_free(j)      free memory for namelist j. If j == MaxNL free for      */
/*                  all name lists.                                         */

void nl_free(int j)
{
    register int i;

    if (j >= 0 && j < MaxNL) {
        if (NLNV[j] > 0) {
            free((char *)NLVIdx[j]);
            memrq(-NLNV[j],sizeof(short));
            NLNV[j] = 0;
            NLName[j][0] = '\0';
            NNL--;
        }
    }
    else if (j == MaxNL) {
        for (i = 0; i < MaxNL; ++i) {
            if (NLNV[i] > 0) {
                free((char *)NLVIdx[i]);
                memrq(-NLNV[i],sizeof(short));
                NLNV[i] = 0;
                NLName[i][0] = '\0';
            }
        }
        NNL = 0;
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

char *get_nvi(char *s,int *n,int nv,short *vidx,int *nb)
{
    register int i,j,k;
    register char c,*p,*q;
    int vn,vn1,l,nl,jj,bflag;
    char vname[VNLMax + 1],vname1[VNLMax + 1];

    for (i = 1; i <= 10; ++i)
        VLTyp[i] = 0;

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
        l = get_vnlen(p);
        if (l <= 0)
            break;
        if (l > VNLMax) {
            i = -1;
            break;
        }
        q = p + l;
        c = *q;
        *q = '\0';
        vn = get_vidx(p,vname);
                 
        if (vn < 0) {
            j = nl_check(p);
            *q = c;
            if (j < 0) {
                i = -1;
                break;
            }
            nl = NLNV[j];   
            for (k = 0; k < nl; ++k) {
                vn = NLVIdx[j][k];
                if (i < nv) {
                    if (bflag == 0)
                        vidx[i] = vn;               
                    else
                        vidx[i] = -vn - 1;           
                }
                VLTyp[(int)VTyp[vn]] += 1;
                i++;
                if (bflag)
                    *nb += 1;
            }
        }
        else {
            *q = c;
            if (i < nv) {
                if (bflag == 0)
                    vidx[i] = vn;
                else
                    vidx[i] = -vn - 1;
            }
            VLTyp[(int)VTyp[vn]] += 1;
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
            l = get_vnlen(p);
            if (l <= 0 || l > VNLMax) {
                i = -1;
                break;
            }
            q = p + l;
            c = *q;
            *q = '\0';
            vn1 = get_vidx(p,vname1);
                     
            if (vn1 < 0) {
                jj = nl_check(p);
                *q = c;
                if (jj < 0) {
                    i = -1;
                    break;
                }
                nl = NLNV[jj];   
                vn1 = NLVIdx[jj][0];
                strcpy(vname1,VName[vn1]);
            }
            else {
                *q = c;
                nl = 1;
            }
            j = VNxt[vn];
            while (j >= 0 && j <= vn1) {
                if (!strcmp(vname1,VName[j]))  
                    break;    
                if (i < nv) {
                    if (bflag == 0)
                        vidx[i] = j;
                    else
                        vidx[i] = -j - 1;
                }
                VLTyp[(int)VTyp[j]] += 1;
                i++;
                if (bflag)
                    *nb += 1;
                j = VNxt[j];
            }
            if (i < nv) {
                if (bflag == 0)
                    vidx[i] = vn1;
                else
                    vidx[i] = -vn1 - 1;
            }
            VLTyp[(int)VTyp[vn1]] += 1;
            i++;
            if (bflag)
                *nb += 1;
            for (k = 1; k < nl; ++k) {
                vn1 = NLVIdx[jj][k];
                if (i < nv) {
                    if (bflag == 0)
                        vidx[i] = vn1;             
                    else
                        vidx[i] = -vn1 - 1;         
                }
                VLTyp[(int)VTyp[vn1]] += 1;
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

char *get_nvia(char *s,int *n,int opt,int *nb)
{
    char *p;

    p = get_nvi(s,n,0,VLVIdx,nb);
    if (*n < 1) {
        if (opt)
            printf1("Syntax error or undefined variables.\n");
        return(p);
    }
    if (alloc_vl(*n)) {
        *n = -2;
        return(p);
    }
    p = get_nvi(s,n,VLNV,VLVIdx,nb);
    if (*n < 1) {
        if (opt)
            printf1("Syntax error or undefined variables.\n");
        alloc_vl(0);
    }
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  alloc_vl(n)     If n > 0 allocate VLVIdx[], otherwise free.             */
/*                  Return 0 if OK, -1 if error.                            */

int alloc_vl(int n)
{
    if (VLNV > 0) {
        free((char *)VLVIdx);
        memrq(-VLNV,sizeof(short));
        VLNV = 0;
    }
    if (n > 0) {
        if (!(VLVIdx = (short *)calloc(n,sizeof(short)))) { 
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(short));
        VLNV = n;
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

int recode(void)
{
    register int i,j,ii;
    int err,l,m,n,vb,vn,w1,w2,wflg,iflg,nv,nc,n2,n3;
    register char c,*p,*q,*q1,*s;
    char vname[VNLMax + 1];
    double id,lastid;

    m = 0;
    vb = err = -1;
    if (check_cmd(0))
        return(-1);

    if (alloc_acx(NOC))
        goto RECFin;

    p = CmdBuf + 6;
    if (*p++ != '(') {
        p_err(-1,1);
        goto RECFin;
    }
    while (*p) {
        if (!strncmp(p,"dblock=",7)) {
            q = p + 7;
            l = get_vnlen(q);
            if (l < 1 || l > VNLMax) {
                printf1("Error: %s.\n",p);
                goto RECFin;
            }
            vb = get_vidx1(q,vname);            
            if (vb < 0) {
                printf1("Undefined variable: %s.\n",p);
                goto RECFin;
            }   
            if (VTyp[vb] == 1) {
                printf1("Error: cannot use string variables for dblock.\n");
                goto RECFin;
            }
            printf1("Block mode with variable: %s\n",vname);
            p = q + l;
            if (*p == ',')
                p++;
        }
        else {
            l = get_vnlen(p);
            if (l < 1 || l > VNLMax) {
                printf1("Error (in variable name?): %s.\n",p);
                goto RECFin;
            }
            q = p + l;
            wflg = 0;
            if (*q == '[') {
                if (sscanf(q,"[%d.%d]",&w1,&w2) != 2) {
                    p_err(-1,1);
                    goto RECFin;
                }
                wflg = 1;
                q = skip_dbl(q + 1);
                if (*q++ != ']') {
                    p_err(-1,1);
                    goto RECFin;
                }
            }
            if (*q++ != '=') {
                p_err(-1,1);
                goto RECFin;
            }
            vn = get_vidx1(p,vname);            

            if (vn < 0) {
                printf1("Undefined variable: %s.\n",p);
                goto RECFin;
            }   
            q1 = skip_expr(q);
            c = *q1;
            *q1 = '\0';

            if (VTyp[vn] == 1 || VTyp[vn] == 2) {
                printf1("Recoding: %s\n",p);
                printf1("Error: cannot recode numerical constants and strings.\n");
                goto RECFin;
            }
            l = strlen(q);      /* check variable's definition */
            if (l == 0) {
                printf1("Syntax error: %s\n",p);
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
                if (alloc_acc(2 * l + 100))
                    goto RECFin;
                    
                l = sub_ifs(q,AcC);     /* change to if() operators */
                q = AcC;
            }
            free(VDef[vn]);
            memrq(-strlen(VDef[vn]) - 1,1);

            if (!(VDef[vn] = (char *)calloc(l + 1,sizeof(char)))) {
                p_err(-2,1);
                goto RECFin;
            }
            memrq(l + 1,sizeof(char));
            strcpy(VDef[vn],q);
                                  
            if ((n = v_parse(q,0)) < 0 || ESCnt <= 0) {
                printf1("Syntax error (%d) in: %s\n",n,p);
                if (n < 0)
                    prn_emsg1(n);
                goto RECFin;
            }
            printf1("Recoding: ");
            prn_vname(vn);
            printf1(" = %s\n",VDef[vn]);

            /* check type */

            check_expr(ESCnt,ESTyp,&nv,&nc,&n2,&n3,0);
            if (nc > 0) {
                p_err(-41,1); 
                goto RECFin;
            }
            if (n3 > 0) {
                p_err(-40,1);  
                goto RECFin;
            }
            if (n2 > 0)
                VTyp[vn] = 4;
            else VTyp[vn] = 3;

            *q1 = c,
             
            ii = 0;
            if (vb < 0)  
                i = NOC;
            else {
                i = 0;
                lastid = get_data(vb,i);
            }
            while (i < NOC) {
                id = get_data(vb,i);
                if (fabs(id - lastid) > EPSI1) {

                    n = v_eval2(-1,ESCnt,ESTyp,ESVal,ESIdx,AcX,ii,i,0,0);
                    if (n) {
                        printf1("Error: can't evaluate expression.\n");
                        prn_emsg2(n);
                        goto RECFin;
                    }
                    for (j = ii; j < i; ++j)  
                        put_data(AcX[j - ii],vn,j);
 
                    ii = i;
                    lastid = id;
                }
                i++;
            }
            n = v_eval2(-1,ESCnt,ESTyp,ESVal,ESIdx,AcX,ii,i,0,0);
            if (n) {
                printf1("Error: can't evaluate expression.\n");
                prn_emsg2(n);
                goto RECFin;
            }
            for (j = ii; j < i; ++j)     
                put_data(AcX[j - ii],vn,j);
  
            if (wflg) {
                makefmt(&w1,&w2,VPFmtS[vn],0,' ',0);
                VPFmt1[vn] = (short)w1;
                VPFmt2[vn] = (short)w2;
            }
            m++;
            if (*q1 == ')')
                break;
            if (*q1++ != ',') {
                p_err(-1,1);
                goto RECFin;
            }
            p = q1;
        }
        if (*p == ')')
            break;
    }
    if (*p != ')') {
        p_err(-1,1);
        goto RECFin;
    }
    if (m == 0)
        printf1("Nothing done.\n");
    else if (TSelFlg && NOC < NOCDM)
        printf1("Warning: used only a subset of data matrix cases.\n");
    err = 0;

RECFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ndvar()     create dummy variables.                                     */
/*              Return 0 if OK, -1 if error.                                */

int ndvar(void)
{
    register int i,j;
    int err,l,vn,vx,r,n,idx,idxn,nd,nn,nb;
    short vv[2];
    register char *p,*q,*s;  
    char vname[VNLMax + 1],vname1[VNLMax + 1],vdef[500];
    double tmp;
               
    idxn = idx = get_nidx();      /* index to first new variable */
    if (idx < 0) {
        p_err(-5,1);
        return(-1); 
    }
    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("Creating dummy variables. ");
    prn_mem();
    tsel_off(1);        /* turn off temporary case selection */

    p = CmdBuf + 5;
    if (*p++ != '(' || !*p) {
        p_err(-1,1);
        goto NDVFin;
    }

    /* allocate memory */
   
    if (alloc_ack(NOC + 1))     /* value */
        goto NDVFin;

    if (alloc_acn(NOC + 1))     /* freq */
        goto NDVFin;

    if (alloc_aci(NOC + 1))     /* bf pointer */
        goto NDVFin;

    if (alloc_acr(NOC + 1))                         
        goto NDVFin;

    if (alloc_acs(NOC + 1))                         
        goto NDVFin;

    while (*p) {

        idxn = get_nidx();      /* index to first new variable */
        if (idxn < 0) {
            p_err(-5,1);
            goto NDVFin;
        }
        s = p;
        l = get_vnlen(p);
        if (l < 1 || l > VNLMax) {
            printf1("Error (in variable name?): %s.\n",p);
            goto NDVFin;
        }
        vn = get_vidx1(p,vname);            
        if (vn >= 0) {
            printf1("Error: %s already used.\n",VName[vn]);
            goto NDVFin;
        }
        q = p + l;
        if (!strncmp(q,"=dmul(",6)) {
            printf1("\nCreating interaction effects: %s = dmul(...)\n",vname);
           
            q += 6;
            p = q;
            while (*p && *p != ',')
                p++;
            if (*p != ',') {
                err = -2;
                goto NDVFin;
            }
            *p = '\0';
            q = get_nvia(q,&n,1,&nb);
            if (nb) {
                p_err(-42,1);
                goto NDVFin;
            }
            *p++ = ',';
            if (n < 1 || VLNV < 1)
                goto NDVFin;

            nn = VLNV;
            if (alloc_acns(nn))                         
                goto NDVFin;

            for (j = 0; j < VLNV; ++j)  
                AcNS[j] = VLVIdx[j];
             
            p = get_nvia(p,&n,1,&nb);
            if (n < 1 || VLNV < 1)
                goto NDVFin;

            if (nb) {
                p_err(-42,1);
                goto NDVFin;
            }

            if (ndivar(vname,nn,AcNS,VLNV,VLVIdx))
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
        vx = get_vidx1(q,vname1);            
        if (vx < 0) {
            printf1("Error: %s\n",p);
            printf1("Undefined variable on right-hand side.\n");
            goto NDVFin;
        }
        q += strlen(VName[vx]);
        p = q;
        if (*p++ != ')') {
            err = -2;                         
            goto NDVFin;
        }

        /* get frequency distribution */

        vv[0] = vx;
        r = cfreq1(1,vv,NOC,AcK,AcN,AcI,-1,&tmp);
        if (r <= 0) {
            printf1("Error: cannot create frequency distribution for %s\n",VName[vx]);
            goto NDVFin;
        }
        nd = 0;
        i = 1;
        while (i <= r) {
            j = AcI[i];
            n = AcK[j];
            q = VName[vx];
            if (n >= 0) {
                sprintf(vdef,"%s%d",vname,n);    
                AcR[nd] += AcN[j];
                AcS[nd++] = n;
            }
            else {                    
                sprintf(vdef,"%sM",vname);    
                AcR[nd] += AcN[j];
                AcS[nd++] = -1;
            }    
            l = strlen(vdef);           
            if (l > VNLMax || get_vidx1(vdef,vname1) >= 0) {
                printf1("Error: cannot create variable: %s\n",vdef);
                if (l > VNLMax)
                    printf1("Exceeded max length of variable names.\n");
                else
                    printf1("Variable name already used.\n");
                goto NDVFin;
            }
            sprintf(vdef,"%s<0>[2.0]=%s[%d",vname1,q,n);    
            while (i < r && (n = AcK[AcI[i + 1]]) < 0) {
                AcR[nd - 1] += AcN[AcI[i + 1]];
                l = strlen(vdef);
                if (l > 450) {
                    printf1("Error: exceeded max length of expression.\n");
                    goto NDVFin;
                }
                q = vdef + l;               
                sprintf(q,",%d",n);
                i++;
            }
            q = vdef + strlen(vdef);
            sprintf(q,"]");

            if (save_var(vdef,0)) {       /* save variable definition */
                printf1("Error: cannot create: %s\n",vdef);
                goto NDVFin;
            }
            i++;
        }
        if (alloc_vdat(idxn,1)) {
            printf1("Error: insufficient memory for new variables.\n");
            goto NDVFin;
        }
        for (i = 0; i < NOC; ++i) {
            n = (int)get_data(vx,i);
            l = idxn;
            if (n >= 0) {
                for (j = 0; j < nd; ++j) {
                    if (n == AcS[j]) {
                        l += j;
                        break;
                    }
                }
            }
            put_data(1.0,l,i);
        }
        printf1("\nVariable: ");
        prn_vname(vx);
        printf1("\nValue        Cases   Pcnt  Dummy variable\n");
        prnchar('-',41,1);

        if (alloc_acns(nd))                         
            goto NDVFin;

        nn = 0;
        for (j = 0; j < nd; ++j) {
            if (AcS[j] < 0)
                printf1("negative  ");
            else {
                printf1("%8d  ",AcS[j]);
                AcNS[nn++] = idxn + j;
            }

            printf1("%8d %6.2lf  ",AcR[j],100.0 * (double)AcR[j] / (double)NOC);
            prn_vname(idxn + j);
            printf1("= %s\n",VDef[idxn + j]);
        }
        if (nn > 0) {               /* create namelist */

            if (c_nlist(vname,nn,AcNS))
                goto NDVFin;
        }
        for (j = 0; j < nd; ++j)
            AcR[j] = 0;

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
            printf1("Syntax error: %s\n",s);
        clear_avar(idxn);
    }
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ndivar(vname,n1,idx1,n2,idx2)   Create interaction effects.             */
/*                                  Return 0 if OK, -1 if error.            */

int ndivar(char *vname,int n1,short *idx1,int n2,short *idx2)
{
    register int i,j,ii;
    int err,l,idx,idxn,ni,nj,nn,nv;
    char *p,vdef[500],vname1[VNLMax + 1];

    idx = get_nidx();      /* index to first new variable */
    if (idx < 0) {
        p_err(-5,1);
        return(-1); 
    }
    if (alloc_acms(n1 * n2 + 1))                         
        return(-1);  

    err = -1;

    printf1("\nCases    Pcnt   Definition\n");
    prnchar('-',43,1);
    nv = 0;
    for (i = 0; i < n1; ++i) {

        for (j = 0; j < n2; ++j) {

            idxn = get_nidx();      /* index to first new variable */
            if (idxn < 0) {
                p_err(-5,1);
                goto NDVIFin;
            }
            sprintf(vdef,"%s_%s",VName[idx1[i]],VName[idx2[j]]);

            l = strlen(vdef);           
            if (l > VNLMax || get_vidx1(vdef,vname1) >= 0) {
                printf1("Error: cannot create variable: %s\n",vdef);
                if (l > VNLMax)
                    printf1("Exceeded max length of variable names.\n");
                else
                    printf1("Variable name already used.\n");
                goto NDVIFin;
            }
            p = vdef + l;
            sprintf(p,"<0>[2.0]=%s&%s",VDef[idx1[i]],VDef[idx2[j]]);

            if (save_var(vdef,0)) {       /* save variable definition */
                printf1("Error: cannot create: %s\n",vdef);
                goto NDVIFin;
            }
            if (alloc_vdat(idxn,1)) {
                printf1("Error: insufficient memory for new variables.\n");
                goto NDVIFin;
            }
            nn = 0;
            for (ii = 0; ii < NOC; ++ii) {
                ni = (int)get_data(idx1[i],ii);
                if (ni) {
                    nj = (int)get_data(idx2[j],ii);
                    if (nj) {
                        put_data(1.0,idxn,ii);
                        nn++;
                    }
                }
            }
            printf1("%6d %6.2lf   ",nn,100.0 * (double)nn / (double)NOC);
            printf1("%s",VName[idxn]);
            prnchar(' ',16 - strlen(VName[idxn]),0);
            printf1("= %s\n",VDef[idxn]);

            AcMS[nv++] = idxn;
        }
    }
    if (c_nlist(vname,nv,AcMS))
        goto NDVIFin;
    
    err = 0;

NDVIFin:
    if (err)   
        clear_avar(idx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  c_nlist(vname,n,vidx)   Create namelist                                 */
/*                          Return 0 if OK, -1 if error.                    */

int c_nlist(char *vname,int n,short *vidx)
{
    register int j,k;

    k = -1;
    for (j = 0; j < MaxNL; ++j) {
        if (NLNV[j] == 0) {
            k = j;
            break;
        }
    }
    if (k < 0) {
        printf1("\nExceeded max number of namelists.\n");
        return(-1);
    }
    else {
        strcpy(NLName[k],vname);
        if (!(NLVIdx[k] = (short *)calloc(n,sizeof(short)))) {  
            p_err(-2,1);
            return(-1);  
        }
        memrq(n,sizeof(short));
        NLNV[k] = n;
        for (j = 0; j < n; ++j)  
            NLVIdx[k][j] = vidx[j];
        NNL++;

        printf1("\nNew namelist: %s = %s",vname,VName[vidx[0]]);
        if (n > 1)
            printf1(",...,%s",VName[vidx[n - 1]]);
        newline();
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  get_mxvlen(n,vidx)  Return max length of variable names in vidx[]       */

int get_mxvlen(int n,short *vidx)
{
    register int i,l,ml;

    ml = 0;
    for (i = 0; i < n; ++i) {
        l = strlen(VName[vidx[i]]);
        if (ml < l)
            ml = l;        
    }
    return(ml);
}

/* ------------------------------------------------------------------------ */
/*  prn_vlist(n,vidx)   Print names of variables in vidx[]                  */

void prn_vlist(int n,short *vidx)
{
    register int i;

    for (i = 0; i < n; ++i) {
        if (i)
            printf1(",");
        printf1("%s",VName[vidx[i]]);
    }
}

/* ------------------------------------------------------------------------ */
/*  svb_alloc(n)    If n > 0 allocate SVBuf, otherwise free.                */
/*                  Return 0 if OK, -1 if error.                            */

int svb_alloc(int n)
{
    if (SVBufA > 0) {
        free((char *)SVBuf);
        memrq(-SVBufA,sizeof(char));
        SVBufA = SVBufLen = 0;
    }
    if (n > 0) {
        if (!(SVBuf = (char *)calloc(n + 2,sizeof(char)))) { 
            p_err(-2,1);    
            return(-1);
        }
        memrq(n + 2,sizeof(char));
        SVBufA = n + 2;
        SVBufLen = n;
    }         
    return(0);
}
