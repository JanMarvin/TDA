/****************************************************************************/
/*  t_mat                                                                   */
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
#include "t_gf.h"
#include "t_alloc.h"
#include "t_matc.h"
#include "t_matf.h"
#include "t_imat.h"
#include "t_gcmd.h"
#include "t_gdd.h"
#include "t_mdat.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_mat                                                      */

int t_mat(void);
void mdefcpy(char *d,char *s);
void m_cmdmsg(void);
int alloc_local(char *cmd);
void free_local(void);
int check_local(char *name);
int alloc_mat(void);
void mat_free(void);
int mat_alloc(int idx,int opt);
int mat_getidx(char *mname,int opt);   
int mat_newidx(char *mname,int m,int n);
int mat_newmat(char *mname,int m,int n);
int mat_ncheck(char *p,int opt);
void mat_err(int opt);   
int mat_info(void);           
int mfree(char *def); 
int mdef(char *def); 
int m_mdefb(char *def); 
int mdeff(char *def); 
int mdefg(char *def); 
char *get_mname(char *p,char *mname,int opt);
int mat_ncopy(int row,int col,double *x,char *name);

/* ------------------------------------------------------------------------ */
/*  global variables.                                                       */

char **MatLocDef;       /* saves the local definitions                      */
int MaxMat = 0;         /* maximum number of matrices                       */
int *MatAlloc;          /* memory allocation for matrices                   */
double **MatVal;        /* matrix data                                      */
char **MatName;         /* names of matrices                                */
char **MatDef;          /* definition of matrices                           */
int *MatRow;            /* number of rows                                   */
int *MatCol;            /* number of columns                                */
int *MatLoc;            /* number of MacroExLevel if defined as local       */

int PMATFmt1 = 12;      /* print format: mfmt=                              */
int PMATFmt2 =  4;
int PMATFmtF = 0;
char PMATFmtS[20];        

char *MatCmdBuf;        /* saves command for error messages                 */
int MatCmdBufLen = 0;

/*--------------------------------------------------------------------------*/
/*  t_mat()     Entry point for matrix commands.                            */
/*              Return 0 if OK, -1 if error, 1 if command cannot be         */
/*              interpreted.                                                */

int t_mat(void)
{
    register char *p;
    int err,n,m;           

    if (PMATFmtF == 0) {
        PMMFmt1 = 10;                     /* make default print format */
        PMMFmt2 = 4; 
        makefmt(&PMATFmt1,&PMATFmt2,PMATFmtS,0,' ',0);
        PMATFmtF = 1;
    }
    err = 0; 
    p = CmdBuf;

    /* save command for error messages */

    n = strlen(p);
    if (!(MatCmdBuf = (char *)calloc(n+1,sizeof(char)))) {
        p_err(-2,1);
        return(-1);
    }         
    memrq(n + 1,sizeof(char));
    MatCmdBufLen = n + 1;
    strcpy(MatCmdBuf,p);

    if (sscanf(p,"mfmt=%d.%d",&n,&m) == 2) {
        PMATFmt1 = n;       
        PMATFmt2 = m; 
        makefmt(&PMATFmt1,&PMATFmt2,PMATFmtS,0,SEPC,0);
        printf2("%s\n",p);
    }

    else if (!strncmp(p,"mcent(",6))        /* mcent(A,R) */
        err = m_cent(p,0);

    else if (!strncmp(p,"mdcent(",7))       /* mdcent(A,R) */
        err = m_dcent(p);

    else if (!strncmp(p,"mstand(",7))       /* mstand(A,R) */
        err = m_cent(p,1);

    else if (!strncmp(p,"mcross(",7))       /* mcross(X,R) */
        err = m_cross(p);

    else if (!strncmp(p,"mchol(",6))        /* mchol(X,R) */
        err = m_chol(p);

    else if (!strncmp(p,"mag(",4))          /* mag(A,C,R,B) */
        err = m_agg(p);

    else if (!strncmp(p,"mcvec(",6))        /* mcvec(A,C) */
        err = m_vec(p,0);

    else if (!strncmp(p,"mrvec(",6))        /* mrvec(A,V) */
        err = m_vec(p,1);

    else if (!strncmp(p,"mivec(",6))        /* mivec(V,n,A) */
        err = m_ivec(p);

    else if (!strncmp(p,"mcsum(",6))        /* mcsum(A,U) */
        err = m_sum(p,1);

    else if (!strncmp(p,"mrsum(",6))        /* mrsum(A,V) */
        err = m_sum(p,0);

    else if (!strncmp(p,"mdrow(",6))        /* mdrow(X,R) */
        err = m_drow(p,0);

    else if (!strncmp(p,"mdcol(",6))        /* mdcol(X,R) */
        err = m_drow(p,1);

    else if (!strncmp(p,"mdiag(",6))        /* mdiag(A,R) */
        err = m_diag(p,0);

    else if (!strncmp(p,"mdiagd(",7))       /* mdiagd(A,R) */
        err = m_diag(p,1);

    else if (!strncmp(p,"mtransp(",8))      /* mtransp(X,R) */
         err = m_transp(p);

    else if (!strncmp(p,"msqrtd(",7))       /* msqrtd(X,R) */
        err = m_sqrt(p,0);

    else if (!strncmp(p,"msqrti(",7))       /* msqrti(X,R) */
        err = m_sqrt(p,1);

    else if (!strncmp(p,"mnrow(",6))        /* mnrow(A,R) */
        err = m_nrow(p,0);

    else if (!strncmp(p,"mncol(",6))        /* mncol(A,R) */
        err = m_nrow(p,1);

    else if (!strncmp(p,"mnc(",4))          /* mnc(A,x,opt,B) */
        err = m_mnc(p);

    else if (!strncmp(p,"mnorm(",6))        /* mnorm(A,R) */
        err = m_nrow(p,2);

    else if (!strncmp(p,"mnorm1(",7))       /* mnorm1(A,R) */
        err = m_nrow(p,3);

    else if (!strncmp(p,"mnorm2(",7))       /* mnorm2(A,R) */
        err = m_nrow(p,4);

    else if (!strncmp(p,"mtrace(",7))       /* mtrace(A,R) */
        err = m_nrow(p,5);

    else if (!strncmp(p,"mnum(",5))         /* mnum(x,d,n,A) */
        err = m_num(p);

    else if (!strncmp(p,"mevs(",5))         /* mevs(A,E,EV) */
        err = m_mevs(p);

    else if (!strncmp(p,"mev(",4))          /* mev(A,ER,EI,EV) */
        err = m_mev(p);

    else if (!strncmp(p,"mginv(",6))        /* mginv(A,R) */
        err = m_ginv(p);

    else if (!strncmp(p,"msvd1(",6))        /* msvd(A,Q,U,V*/
        err = m_svd(p,1);

    else if (!strncmp(p,"msvd(",5))         /* msvd(A,Q) */
        err = m_svd(p,0);

    else if (!strncmp(p,"mwvec1(",7))       /* mwvec1(A,W,T,R) */
        err = m_wvec1(p);

    else if (!strncmp(p,"mwvec(",6))        /* mwvec(A,W,R) */
        err = m_wvec(p);

    else if (!strncmp(p,"mscal1(",7))       /* mscal1(A,B) */
        err = m_scal1(p);

    else if (!strncmp(p,"mdefb(",6))        /* mdefb(A,bn) */
        err = m_mdefb(p);

    else if (!strncmp(p,"mdefc(",6))        /* mdefc(m,n,d,A) */
        err = m_mdefc(p);

    else if (!strncmp(p,"mdefi(",6))        /* mdefi(m,n,A) */
        err = m_mdefi(p);

    else if (!strncmp(p,"minvs(",6))        /* minvs(A,R) */
        err = m_invs(p);

    else if (!strncmp(p,"minvd(",6))        /* minvd(A,R) */
        err = m_invd(p);

    else if (!strncmp(p,"mple(",5))         /* mple(T,C,F,D) */
        err = m_ple(p);

    else if (!strncmp(p,"mnvar(",6))        /* mnvar(X) */
        err = m_nvar(p);

    else if (!strncmp(p,"mpsym(",6))        /* mpsym(A,P,B) */
        err = m_mperm(p,0);

    else if (!strncmp(p,"mprow(",6))        /* mprow(A,P,B) */
        err = m_mperm(p,1);

    else if (!strncmp(p,"mpcol(",6))        /* mpcol(A,P,B) */
        err = m_mperm(p,2);

    else if (!strncmp(p,"mpr",3))           /* print matrix */
        err = m_print(p);  

    else if (!strncmp(p,"mkp(",4))          /* mkp(A,B,R) */
        err = m_mkp(p);

    else if (!strncmp(p,"mlsei1(",7))       /* mlsei1(S,me,mi,R) */
        err = m_mlsei1(p);

    else if (!strncmp(p,"mnls(",5))         /* mnls(S,l,R) */
        err = m_mls(p,4);

    else if (!strncmp(p,"mlse(",5))         /* mlse(S,R) */
        err = m_mls(p,1);

    else if (!strncmp(p,"mlsi(",5))         /* mlsi(S,R) */
        err = m_mls(p,2);

    else if (!strncmp(p,"mlsei(",6))        /* mlsei(S,me,mi,R) */
        err = m_mls(p,3);

    else if (!strncmp(p,"mls(",4))          /* mls(S,R) */
        err = m_mls(p,0);

    else if (!strncmp(p,"mqpc(",5))         /* mqpc(C,D,A,B,me,X) */
        err = m_mqp(p,2);

    else if (!strncmp(p,"mqpb(",5))         /* mqpb(C,D,XL,XU,X) */
        err = m_mqp(p,1);

    else if (!strncmp(p,"mqp(",4))          /* mqp(C,D,X) */
        err = m_mqp(p,0);

    else if (!strncmp(p,"mlpi(",5))         /* mlpi(A,B,X) */
        err = m_mlpi(p);

    else if (!strncmp(p,"mlp1(",5))         /* mlp1(T,p,X,Y) */
        err = m_mlp(p,1);

    else if (!strncmp(p,"mlp(",4))          /* mlp(T,X,Y) */
        err = m_mlp(p,0);

    else if (!strncmp(p,"mmul(",5))         /* mmul(A,...,R,) */
        err = m_mul(p);
  
    else if (!strncmp(p,"mcath(",6))        /* mcath(A,...,R) */
        err = m_cat(p,0);

    else if (!strncmp(p,"mcatv(",6))        /* mcatv(A,...,R) */
        err = m_cat(p,1);

    else if (!strncmp(p,"mcathv(",7))       /* mcathv(A,...,R) */
        err = m_cat(p,2);

    else if (!strncmp(p,"msrow(",6))        /* msrow(A,D,R) */
        err = m_srow(p,0);

    else if (!strncmp(p,"mscol(",6))        /* mscol(A,D,R) */
        err = m_srow(p,1);

    else if (!strncmp(p,"msort(",6))        /* msort(A,D,R) */
        err = m_sort(p,0);

    else if (!strncmp(p,"msort1(",7))       /* msort1(A,D,R) */
        err = m_sort(p,2);

    else if (!strncmp(p,"mrank(",6))        /* mrank(A,D,R) */
        err = m_sort(p,1);

    else if (!strncmp(p,"msetv(",6))        /* msetv command */
        err = m_setv(p);  

    else if (!strncmp(p,"mdeff(",6))        /* matrix definition */
        err = mdeff(p);  

    else if (!strncmp(p,"mdefg",5))         /* matrix definition */
        err = mdefg(p);  

    else if (!strncmp(p,"mdef(",5))         /* matrix definition */
        err = mdef(p);  

    else if (!strcmp(p,"mdef"))             /* info about matrices */
        err = mat_info();

    else if (!strncmp(p,"mfree",5))         /* mfree(X) */
        err = mfree(p);  

    else if (!strncmp(p,"mexpr1(",7))       /* mexpr1(B,expression,A) */
        err = m_expr1(p);

    else if (!strncmp(p,"mexpr(",6))        /* mexpr(expression,A) */
        err = m_expr(p);

    else if (!strncmp(p,"mexp(",5))         /* mexp(expression,A) */
        err = m_exp(p);

    else if (!strncmp(p,"mbrr(",5))         /* mbrr(ns,nu,A) */
        err = m_brr(p);

    else if (!strncmp(p,"midf(",5))         /* midf(XL,XU,DL,DU,DM) */
        err = m_midf(p);

    else if (!strncmp(p,"midf1(",6))        /* midf1(XL,XU,F) */
        err = m_midf1(p);

    else if (!strncmp(p,"midf2(",6))        /* midf2(XL,XU,F) */
        err = m_midf2(p);

    else if (!strncmp(p,"midf3(",6))        /* midf3(XL,XU,XL1,XU1) */
        err = m_midf3(p);

    else if (!strncmp(p,"mtrim(",6))        /* mtrim(A,ca,ra,cb,rb,R) */
        err = m_mtrim(p);

    else if (!strncmp(p,"mmp1(",5))         /* mmp1(X,T,P) */
        err = m_mmp(p,1);

    else if (!strncmp(p,"mmp2(",5))         /* mmp2(X,T,P) */
        err = m_mmp(p,2);

    else if (!strncmp(p,"mmp(",4))          /* mmp(X,P) */
        err = m_mmp(p,0);

    else if (!strncmp(p,"mpz(",4))          /* mpz(A,B,P) */
        err = m_mpz(p);

    else if (!strncmp(p,"mpbl(",5))         /* mpbl(A,B,P,N,U) */
        err = m_mpb(p,0);

    else if (!strncmp(p,"mpbu(",5))         /* mpbu(A,B,P,N,U) */
        err = m_mpb(p,1);

    else if (!strncmp(p,"mqap(",5))         /* mqap(F,D,C,P) */
        err = m_mqap(p);

    else if (!strncmp(p,"mpinv(",6))        /* mpinv(P,Q) */
        err = m_mpinv(p);

    else if (!strncmp(p,"mcel(",5))         /* mcel(A,x,L) */
        err = m_mcel(p);

    else if (!strncmp(p,"mpfit(",6))        /* mpfit(A,R,C,B) */
        err = m_mpfit(p);

    else if (!strncmp(p,"mch(",4))          /* mch(A,B) */
        err = m_mch(p);

    else if (!strncmp(p,"mpit1(",6))        /* mpit1(F,N,Z,n,R) */
        err = m_mpit(p,1);

    else if (!strncmp(p,"mpit(",5))         /* mpit(F,N,n,R) */
        err = m_mpit(p,0);

    else if (!strncmp(p,"mldes(",6))        /* mldes(Z,G,D) */
        err = m_mldes(p);

    else if (!strncmp(p,"mkmet(",6))        /* mkmet(A,D) */
        err = m_kmet(p);

    else
        err = 1;

    if (MatCmdBufLen > 0) {    
        free(MatCmdBuf);
        memrq(-MatCmdBufLen,sizeof(char));
        MatCmdBufLen = 0;            
    }
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  mdefcpy(d,s)   copy matrix definition from s to d, maximal MatDefLen    */

void mdefcpy(char *d,char *s)
{
    strncpy(d,s,MatDefLen);
    *(d + MatDefLen) = '\0';
}

/*--------------------------------------------------------------------------*/
/*  m_cmdmsg    Print Error in command: ...                                 */

void m_cmdmsg(void)
{
    if (SILENTFlg >= 0 && MatCmdBufLen > 0)
        printf1("Error in command: %s\n",MatCmdBuf);
}

/*--------------------------------------------------------------------------*/
/*  alloc_local     If MacroExLevel > 0 and < MaxMat                        */
/*                  save string in MatLocDef[MacroExLevel].                 */
/*                  Syntax of local command is:                             */
/*                                                                          */
/*                  local(A,B,...);                                         */
/*                                                                          */
/*                  where A,B,... are matrix names.                         */
/*                                                                          */
/*                  Return 0 if OK, -1 if insufficient memory.              */

int alloc_local(char *cmd)
{
    int len;
    register char *p;

    printf2("local(%s [macro level %d]\n",cmd,MacroExLevel);

    if (MacroExLevel < 1 || MacroExLevel >= MaxMat)
        return(0);

    len = strlen(cmd);
    p = cmd + len;
    while (*--p == ')')
        len--;
    *++p = '\0';
    if (len < 1)   
        return(0);

    if (!(MatLocDef[MacroExLevel] = (char *)calloc(len+1,sizeof(char)))) {
        p_err(-2,1);
        return(-1);
    }         
    memrq(len + 1,sizeof(char));
    strcpy(MatLocDef[MacroExLevel],cmd);
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  free_local      If MacroExLevel > 0 and < MaxMat                        */
/*                  free local matrices.                                    */

void free_local(void)
{
    register int i;
          
    if (MacroExLevel > 0 && MacroExLevel < MaxMat) {
        for (i = 0; i < MaxMat; ++i) {
            if (MatAlloc[i] && MatLoc[i] == MacroExLevel)    
                mat_alloc(i,0);
        }
        if (MatLocDef[MacroExLevel] != NULL) {  
            i = strlen(MatLocDef[MacroExLevel]) + 1;
            free(MatLocDef[MacroExLevel]);
            memrq(-i,sizeof(char));
        }
    }
}

/*--------------------------------------------------------------------------*/
/*  check_local(A)  Return 1 if name A occurs in local list for current     */
/*                  MacroExLevel, otherwise return 0.                       */

int check_local(char *name)
{
    register char c,*p,*q;

    if (MacroExLevel < 1 || MacroExLevel >= MaxMat)
        return(0);

    p = MatLocDef[MacroExLevel];            
    if (p != NULL) {
        while (*p) {
            q = p;
            while (*q && *q != ',')
                q++;
            c = *q;
            *q = '\0';
            if (!strcmp(name,p)) {
                *q = c;
                return(1);
            }
            *q = c;
            p = q;
            if (*p++ != ',')
                break;
        }
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  alloc_mat()     Allocate memory for MaxNat matrices. This memory is     */
/*                  only allocated once at the beginning of the program     */
/*                  and is never free'd.                                    */
/*                  return 0 if OK, -1 if error.                            */

int alloc_mat(void)
{
    int i;

    if (!(MatAlloc = (int *)calloc(MaxMat,sizeof(int))) ||
        !(MatLoc   = (int *)calloc(MaxMat,sizeof(int))) ||
        !(MatRow   = (int *)calloc(MaxMat,sizeof(int))) ||
        !(MatCol   = (int *)calloc(MaxMat,sizeof(int))) ||
        !(MatLocDef = (char **)calloc(MaxMat,sizeof(char *))) ||
        !(MatName  = (char **)calloc(MaxMat,sizeof(char *))) ||
        !(MatDef   = (char **)calloc(MaxMat,sizeof(char *))) ||
        !(MatVal   = (double **)calloc(MaxMat,sizeof(double *))))   
        return(-1);
              
    memrq(MaxMat,4 * sizeof(int) + 3 * sizeof(char *) + sizeof(double *));

    for (i = 0; i < MaxMat; ++i) {
        MatLocDef[i] = NULL;
        if (!(MatName[i] = (char *)calloc(VNLMax + 1,sizeof(char))) ||
            !(MatDef[i] = (char *)calloc(MatDefLen + 1,sizeof(char))))  
            return(-1);
        *(MatDef[i] + MatDefLen) = '\0';
    }
    memrq(MaxMat * (VNLMax + MatDefLen + 2),sizeof(char));
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  mat_free()      free memory for all allocated matrices.                 */
   
void mat_free(void)
{
    int i;  
       
    for (i = 0; i < MaxMat; ++i) {
        if (MatAlloc[i])   
            mat_alloc(i,0);
    }
}

/*--------------------------------------------------------------------------*/
/*  mat_alloc(idx,opt)      allocate (opt = 1) or free (opt = 0) memory     */
/*                          for matrix with index idx.                      */
/*                          return 0 if OK, -1 if error.                    */

int mat_alloc(int idx,int opt)
{
    int i,mn;
        
    if (idx == MPLogIdx)  
        MPLogIdx = -1;
    else if (idx == MPParIdx)  
        MPParIdx = -1;
    else if (idx == MPCovIdx)  
        MPCovIdx = -1;
    else if (idx == MPGradIdx) {
        MPGradIdx = -1;
        MPGradRow = MPGradCol = 0;
    }
    else if (idx == MPResIdx) {
        MPResIdx = -1;
        MPResRow = MPGradCol = 0;
    }
    if (GD_TYP == 8) {                  /* check for graph definition */
        for (i = 0; i < GD_NG; ++i) {
            if (GD_EV[i] == idx) {
                gdd_free(1);
                break;
            }
        }
    }

    mn = MatRow[idx] * MatCol[idx] + 1;
    if (opt) {

        if (MatAlloc[idx] || mn <= 1) {
            printfe("MAT_ALLOC ERROR: mn = %d\n",mn);
            gerr_exit(71);
        }
        if (!(MatVal[idx] = (double *)calloc(mn,sizeof(double)))) {
            m_cmdmsg();
            printf1("Error: insufficient memory for (%d,%d) matrix.\n",        
                                                   MatRow[idx],MatCol[idx]);
            return(-1);
        }
        MatAlloc[idx] = 1;

        if (MacroExLevel > 0 && MacroExLevel < MaxMat) {
            if (check_local(MatName[idx]))
                MatLoc[idx] = MacroExLevel;
        }
    }
    else {
        if (MatAlloc[idx] == 0)
            gerr_exit(72);

        free((char *)MatVal[idx]);
        mn = -mn;
        MatAlloc[idx] = 0;
    }
    memrq(mn,sizeof(double));
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  mat_getidx(mname,opt)   return pointer to mname if defined,             */
/*                          otherwise -1.                                   */
/*                          if opt == 1 print Error: ... not defined        */

int mat_getidx(char *mname,int opt)
{
    register int i;

    for (i = 0; i < MaxMat; ++i) {
        if (MatAlloc[i] && !strcmp(mname,MatName[i])) {
            if (MacroExLevel > 0 && MacroExLevel < MaxMat && check_local(mname)) {
                if (MatLoc[i] == MacroExLevel)  
                    return(i);           
            }
            else  
                return(i);
        }
    }
    if (opt == 1) {
        m_cmdmsg();
        printf1("%s not defined.\n",mname);
    }
    return(-1);
}

/*--------------------------------------------------------------------------*/
/*  mat_newidx(mname,m,n)   Look for free pointer for a new matrix. If      */
/*                          found, set mname and m,n and return pointer,    */
/*                          otherwise -1.                                   */
/*                                                                          */
/*  Check also for name conflicts with variable and namelist names.         */
/*  In case of conflicts, make error message and return -1.                 */

int mat_newidx(char *mname,int m,int n)
{
    register int i;

    if (check_local(mname) == 0) {
        i = VIFirst;
        while (i >= 0) {
            if (!strcmp(mname,VName[i])) {
                m_cmdmsg();
                printf1("Error: name conflict with variable: %s\n",mname);
                return(-1);
            }
            i = VNxt[i];
        }
        if (nl_check(mname) >= 0) {
            m_cmdmsg();
            printf1("Error: name conflict with namelist: %s\n",mname);
            return(-1);
        }
    }
    for (i = 0; i < MaxMat; ++i) {
        if (MatAlloc[i] == 0) {
            strcpy(MatName[i],mname);
            MatRow[i] = m;
            MatCol[i] = n;
            return(i);           
        }
    }
    m_cmdmsg();
    printf1("Error: already allocated max number (%d) of matrices.\n",MaxMat);
    return(-1);
}

/*--------------------------------------------------------------------------*/
/*  mat_newmat(mname,m,n)   get a new matrix and allocate memory.           */
/*                          return index if OK, -1 if error.                */

int mat_newmat(char *mname,int m,int n)
{
    int i,idx;

    idx = mat_getidx(mname,0);  /* check if matrix already defined */
    if (idx >= 0) {
        if (GD_TYP == 8) {                   /* check if used for graph */
            for (i = 0; i < GD_NG; ++i) {
                if (GD_EV[i] == idx) {
                    printf1("Cannot automatically overwrite matrix %s\n",MatName[idx]);
                    printf1("This matrix is used for a graph definition.\n");
                    return(-1);
                }
            }
        }
        mat_alloc(idx,0);       /* free memory */
    }
    idx = mat_newidx(mname,m,n);      /* get new index */
    if (idx < 0)
        return(-1);      

    if (mat_alloc(idx,1)) {   /* allocate memory */
        m_cmdmsg();
        printf1("Error: insufficient memory for new matrix.\n");
        return(-1);
    }
    return(idx);
}

/*--------------------------------------------------------------------------*/
/*  mat_ncheck(p,opt)   check whether p points to a matrix name.            */
/*                      return 0 if OK, otherwise -1.                       */

int mat_ncheck(char *p,int opt)
{
    int l;
       
    if (check_vname(p) == 0)                            
        goto NCErr;
          
    l = get_vnlen(p);
    if (l < 1 || l > VNLMax)              
        goto NCErr;
    return(0);
   
NCErr:
    if (opt) {
        m_cmdmsg();
        printf1("Error in matrix name.\n");
    }
    return(-1);
}

/*--------------------------------------------------------------------------*/
/*  mat_err(opt)    print error message                                     */
/*                  opt 0   syntax error                                    */
/*                  opt 1   square matrix required                          */
/*                  opt 2   error in matrix dimensions                      */
/*                  opt 3   insufficient memory                             */
/*                  opt 4   matrix names must be different                  */
/*                  opt 5   equal number of rows required                   */
/*                  opt 6   error in index definition                       */
/*                  opt 7   column vector required                          */
/*                  opt 8   scalar expression(s) required                   */
/*                  opt 9   error in selection indices                      */
/*                  opt 10  error in matrix indices                         */
/*                  opt 11  need positive interval width                    */
/*                  opt 12  command requires block mode                     */
/*                  opt 13  block number out of range                       */
/*                  opt 14  error in number of constraints                  */
/*                  opt 15  error in resulting matrix dimensions            */
/*                  opt 16  error in index vector                           */
/*                  opt 17  error in permutation vector                     */
/*                  opt 18  inconsistent parametes                          */

void mat_err(int opt) 
{
    m_cmdmsg();

    if (opt == 0)
        printf1("Syntax error (or unknown strings in matrix expression).\n");
    else if (opt == 1)
        printf1("Error: square matrix required.\n");
    else if (opt == 2)
        printf1("Error in matrix dimensions.\n");
    else if (opt == 3)
        printf1("Error: insufficient memory.\n");
    else if (opt == 4)
        printf1("Error: matrix names must be different.\n");
    else if (opt == 5)
        printf1("Error: required equal number of rows.\n");
    else if (opt == 6)
        printf1("Error in definition of indices.\n");
    else if (opt == 7)
        printf1("Error: column vector required.\n");
    else if (opt == 8)
        printf1("Error: scalar expression(s) required.\n");
    else if (opt == 9)
        printf1("Error: in selection indices.\n");
    else if (opt == 10)
        printf1("Error: in matrix indices.\n");
    else if (opt == 11)
        printf1("Error: need positive interval width.\n");
    else if (opt == 12)
        printf1("Error: command requires block mode.\n");
    else if (opt == 13)
        printf1("Error: block number out of range (1 -- %d).\n",BNOC);
    else if (opt == 14)
        printf1("Error in number of constraints.\n");
    else if (opt == 15)
        printf1("Error in resulting matrix dimensions.\n");
    else if (opt == 16)
        printf1("Error in index vector (column or row).\n");
    else if (opt == 17)
        printf1("Error in permutation vector.\n");
    else if (opt == 18)
        printf1("Error: inconsistent parameters.\n");

}

/*--------------------------------------------------------------------------*/
/*  mat_info()  print info about currently defined matrices.                */

int mat_info(void)            
{
    register int i,n;
    int maxlen = 16;

    n = 0;
    for (i = 0; i < MaxMat; ++i) {
        if (MatAlloc[i])
            n++;
    }
    if (n == 0) {
        printf1("mdef: no matrices defined.\n");
        return(0);
    }
    printf1("mdef\nmatrix            loc  rows  columns  definition\n");
    prnchar('-',48,1);
                                 
    for (i = 0; i < MaxMat; ++i) {
        if (MatAlloc[i]) {
            printf1(MatName[i]);          
            prnchar(' ',maxlen - strlen(MatName[i]),0);
            printf1("  %3d %5d %6d    %s\n",MatLoc[i],
                                 MatRow[i],MatCol[i],MatDef[i]);
        }
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  mfree(def)  def is pointer to string mfree(X). Free memory for          */
/*              matrix X. Or mfree without arguments free's all matrices.   */
/*              return 0 if OK, -1 if error.                                */
        
int mfree(char *def)  
{
    register char *p,*q;
    int idx; 

    printf2("%s\n",def);
    if (!strcmp(def,"mfree")) {
        mat_free();
        return(0);
    }
    p = def + 5;
    if (*p++ != '(') {
        mat_err(0);
        return(-1);   
    }                                                   
    q = p;
    if (!*q) {
        mat_err(0);
        return(-1);   
    }
    while (*q && *q != ',' && *q != ')')
        q++;

    if (*q != ')') {
        mat_err(0);
        return(-1);   
    }
    *q = '\0';

    idx = mat_getidx(p,1);
    if (idx < 0)
        return(-1);

    mat_alloc(idx,0);
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  mdef(def)   Definition of a matrix. def is pointer to a string:         */
/*  ##          0   mdef(X,m,n);            use data matrix                 */
/*              1   mdef(X,m,n) = X1,...;   use variables                   */
/*              2   mdef(X,m,n) = x11,...,x12,      use data                */
/*                                ...        ,                              */
/*                                xm1,...,xmn;                              */
/*              3   mdef(X,m,n) = fname;    use fname                       */
/*              4   mdef(X);                use full data matrix            */
/*              5   mdef(X) = varlist;      use NOC x varlist               */
/*                                                                          */
/*              return 0 if OK, -1 if error.                                */

int mdef(char *def)  
{
    FILE *fd;
    register int i,j,k,ii;
    register char *p,*q;
    int r,m,n,idx,mtyp,dflag;
    char mtmp[MatDefLen + 1];
    double tmp;

    mdefcpy(mtmp,def);
    printf2("%s\n",mtmp);

    mtyp = 0;
    p = def + 4;
    dflag = 1;
    if (NVAR < 1 || NOC < 1 || !DMDef)  
        dflag = 0;

   
    if (*p++ != '(' || mat_ncheck(p,0)) {
        mat_err(0);                   
        return(-1);
    }
    q = p + get_vnlen(p);             
    if (*q == ',') {
        *q++ = '\0';
        if (sscanf(q,"%d,%d",&m,&n) != 2 || m < 1 || n < 1) {
            mat_err(0);    
            return(-1);
        }
        q = skip_int(q);
        q = skip_int(++q);
        if (*q++ != ')') {
            mat_err(0);       
            return(-1);
        }
    }
    else if (*q == ')') {
        *q++ = '\0';
        if (dflag == 0) {
            m_cmdmsg();
            printf1("Error: no data matrix.\n");
            return(-1);
        }
        m = NOC;
        n = 0;
    }
    else {
        mat_err(0);
        return(-1);
    }
    if (*q == '=') {
                     
        q++;
        if (dflag && check_vname(q)) {
            if (n > 0)
                mtyp = 1;
            else
                mtyp = 5;
        }
        else if (sscanf(q,"%lg",&tmp) == 1)  
            mtyp = 2;
        else if (*q)
            mtyp = 3;
        else {
            mat_err(0);    
            return(-1);
        }
    }
    else {
        if (*q) {
            mat_err(0);       
            return(-1);
        }
        if (dflag == 0) {
            printf1("Error: no data matrix.\n");
            return(-1);
        }
        if (n > 0)
            mtyp = 0;
        else
            mtyp = 4;
    }

    switch (mtyp) {
        case 4:                             /* full data matrix */
            n = NVAR;
        case 0:                             /* part of data matrix */
            idx = mat_newmat(p,m,n);        /* get new matrix */
            if (idx < 0)
                return(-1);        
            ii = m; 
            if (ii > NOC) 
                ii = NOC;

            for (i = 0; i < ii; ++i) {
                k = 1;
                j = VIFirst;
                while (j >= 0) {
                    MatVal[idx][i * n + k] = get_data(j,i);
                    if (++k > n)
                        break;
                    j = VNxt[j];
                }
            }
            break;

        case 2:                     /* defined by numerical entries */

            idx = mat_newmat(p,m,n);     /* get new matrix */
            if (idx < 0)
                return(-1);        

            for (i = 1; i <= m; ++i) {
                ii = (i - 1) * n;
                for (j = 1; j <= n; ++j) {

                    if (sscanf(q,"%lg",&tmp) != 1) {
                        m_cmdmsg();
                        printf1("Error: can't read entry (%d,%d).\n",i,j);
                        return(-1);
                    }
                    MatVal[idx][ii + j] = tmp;
                    q = skip_dval(q);
                    if (i == m && j == n)
                        break;
                    if (*q++ != ',') {
                        if (++j > n) {
                            ++i;
                            j = 1;
                        }
                        m_cmdmsg();  
                        printf1("Syntax error in entry (%d,%d).\n",i,j);
                        return(-1);
                    }
                }
            }
            if (*q)  
                printf1("Warning: data list contains more than %d entries.\n",m * n);

            break;

        case 3:                     /* defined by external file */

            idx = mat_newmat(p,m,n);         /* get new matrix */
            if (idx < 0)
                return(-1);        

            if (!(fd = fopen(q,"r"))) {
                m_cmdmsg();
                printf1("Error: can't open data file: %s\n",q);
                return(-1);     
            }
            if (alloc_acc(RLMaxDef + 1)) {
                fclose(fd);
                return(-1);
            }
            i = 0;
            while (i < m) {

                if (fgets(AcC,RLMaxDef,fd) == NULL || *AcC == 0x1a)
                    break;

                if (check_drec(AcC)) {      /* check for data records */

                    if (++i > m)
                        break;
                    ii = (i - 1) * n;

                    p = AcC;

                    for (j = 1; j <= n; ++j) {    /* read required variables */

                        p = skip_sep(p);    /* skip separator */

                        if (sscanf(p,"%lg",&tmp) != 1) {
                            m_cmdmsg();
                            printf1("Error: can't read entry (%d,%d).\n",i,j);
                            fclose(fd);
                            return(-1);
                        }
                        MatVal[idx][ii + j] = tmp;
                        p = skip_dval(p);
                    }
                }
            }
            fclose(fd);
            alloc_acc(0);
            if (i < m)  
                printf1("Warning: %s contains only %d data records.\n",q,i);
            break;

        case 1:                     /* based on varlist */
        case 5:     
            get_var(q,&r);
            if (r || PMNV < 1)
                return(-1);
            if (mtyp == 5)
                n = PMNV;

            idx = mat_newmat(p,m,n);     /* get new matrix */
            if (idx < 0)
                return(-1);        

            for (i = 0; i < m; ++i) {
                if (i >= NOC)
                    break;
                for (j = 0; j < PMNV; ++j) {
                    if (j >= n)
                        break;
                    MatVal[idx][i * n + j + 1] = get_data((int)PMVIdx[j],i);
                }           
            }
            break;
    }
    strcpy(MatDef[idx],mtmp);
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  m_mdefb(def,bn) Definition of a matrix. def is pointer to a string:     */
/*  ##              The command creates a matrix corresponding to block     */
/*                  bn (may be number or matrix). Command requiers          */
/*                  block mode.                                             */
/*                                                                          */
/*              return 0 if OK, -1 if error.                                */

int m_mdefb(char *def)  
{
    register int i,j,k,l;
    register char *p,*q;
    int err,n,bn,idx,ii,row,col,iv;
    char mtmp[MatDefLen + 1];

    err = -1;

    mdefcpy(mtmp,def);
    printf2("%s\n",mtmp);

    if (NVAR < 1 || NOC < 1 || !DMDef || BNOC <= 0) {
        mat_err(12);
        goto MDEFBFin;
    }
    p = def + 5;
   
    if (*p++ != '(' || mat_ncheck(p,0)) {
        mat_err(0);                   
        goto MDEFBFin;
    }
    q = p + get_vnlen(p);             

    if (*q != ',') {
        mat_err(0);
        goto MDEFBFin;
    }
    *q++ = '\0';

    if ((q = n_mexpr(q,1,0,&row,&col,&iv,NULL)) == NULL)  
        goto MDEFBFin;  

    if (row != 1 || col != 1) {
        mat_err(8);
        goto MDEFBFin;
    }
    bn = (int)get_mxval(0,1,1);

    if (bn < 1 || bn > BNOC) {
        mat_err(13);    
        goto MDEFBFin;
    }
    bn--;
    n = 0;
    ii = -1;
    for (i = 0; i < NOC; ++i) {
        if (DBlckPtr[i] == bn) {
            n++;
            if (ii < 0)
                ii = i;
        }
    }
    idx = mat_newmat(p,n,NVAR);    
    if (idx < 0)
        goto MDEFBFin;      
    
    l = 0;
    for (i = ii; i < ii + n; ++i) {
        k = 1;
        j = VIFirst;
        while (j >= 0) {
            MatVal[idx][l * NVAR + k++] = get_data(j,i);
            j = VNxt[j];
        }
        l++;
    }
    strcpy(MatDef[idx],mtmp);
    err = 0;   

MDEFBFin:
    mx_free();
    return(err);
}
  
/*--------------------------------------------------------------------------*/
/*  mdeff(def)  Definition of a matrix                                      */
/*              Syntax: mdeff(X) = fname.                                   */
/*              Create matrix X corresponding to data in fname.             */
/*                                                                          */
/*              return 0 if OK, -1 if error.                                */

int mdeff(char *def)  
{
    FILE *fd;
    int idx,m,n,i,j,ii;
    double tmp,*dp;
    char *p,*q,mname[VNLMax + 1];

    printf2("%s\n",def);
    p = def + 6;

    if ((p = get_mname(p,mname,1)) == NULL)
        return(-1);                       
    if (*p++ != ')' || *p++ != '=' || !*p) {
        mat_err(0);
        return(-1);
    }
    if (!(fd = fopen(p,"r"))) {
        m_cmdmsg();
        printf1("Error: can't open data file: %s\n",p);
        return(-1);     
    }
    if (alloc_acc(RLMaxDef + 1)) {
        fclose(fd);
        return(-1);
    }
    n = m = 0;
    while (fgets(AcC,RLMaxDef,fd) && *AcC != 0x1a) {
        if (check_drec(AcC)) {      /* check for data records */

            if (m == 0) {
                q = AcC;
                while (1) {
                    q = skip_sep(q);    /* skip separator */

                    if (sscanf(q,"%lg",&tmp) != 1)  
                        break;

                    n++;
                    q = skip_dval(q);
                }
            }
            m++;
        }
    }
    if (m == 0 || n == 0) {
        m_cmdmsg();
        printf1("Error: file %s does not contain readable data.\n",p);       
        return(-1);
    }
    if ((idx = mat_newmat(mname,m,n)) < 0)
        return(-1);       

    dp = MatVal[idx];
    strcpy(MatDef[idx],def);

    fseek(fd,0,0);
    i = 0;
    while (i < m) {

        if (fgets(AcC,RLMaxDef,fd) == NULL || *AcC == 0x1a)
            break;

        if (check_drec(AcC)) {      /* check for data records */

            if (++i > m)
                break;
            ii = (i - 1) * n;

            q = AcC;

            for (j = 1; j <= n; ++j) {    /* read required variables */

                q = skip_sep(q);    /* skip separator */

                if (sscanf(q,"%lg",&tmp) != 1) {
                    m_cmdmsg();
                    printf1("Error: can't read entry (%d,%d).\n",i,j);
                    fclose(fd);
                    alloc_acc(0);
                    return(-1);
                }
                MatVal[idx][ii + j] = tmp;
                q = skip_dval(q);
            }
        }
    }
    fclose(fd);
    alloc_acc(0);
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  mdefg(gn=...,sc=...) = Matrix name.                                     */
/*                                                                          */
/*  Create adjacency matrix for graph data, gn is graph number.             */
/*  Substitute missing values by sc.                                        */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int mdefg(char *def)  
{
    register int i,j;
    int idx,gn;           
    double sc,tmp;
    char *p,mname[VNLMax + 1];

    printf2("%s\n",def);
    p = def + 5;

    gn = 1;
    sc = -1.0;
    if (*p == '(') {
        p++;
        while (*p) {
            if (sscanf(p,"gn=%d",&gn) == 1 && gn >= 0)
                p = skip_int(p + 3);
            else if (sscanf(p,"sc=%lg",&tmp) == 1) {
                sc = tmp;
                p = skip_dbl(p + 3);
            }
            else if (*p == ')')
                break;
            else {
                mat_err(0);
                return(0);
            }
            if (*p == ',')  
                p++;
            else
                break;
        }
        if (*p++ != ')') {
            mat_err(0);
            return(-1);
        }
    }
    if (*p++ != '=') {
        mat_err(0);
        return(-1);
    }
    if ((p = get_mname(p,mname,1)) == NULL)
        return(-1);                       

    if (gdd_check(gn,1))
        return(-1);   

    if ((idx = mat_newmat(mname,GD_NP,GD_NP)) < 0)
        return(-1);       

    strcpy(MatDef[idx],def);

    for (i = 0; i < GD_NP; ++i) {
        for (j = 0; j < GD_NP; ++j) {
            tmp = gdd_adj(i,j,gn);
            if (tmp < 0.0)
                tmp = sc;
            MatVal[idx][i * GD_NP + j + 1] = tmp;
        }
    }
    return(0);
}
    
/*--------------------------------------------------------------------------*/
/*  get_mname(p,mname,opt)                                                  */
/*                  get mname from pointer p, return pointer to first       */
/*       character after mname. Print message and return NULL if err.      */

char *get_mname(char *p,char *mname,int opt)
{
    register int j;
    register char *q;

    q = mname;
    for (j = 0; j < VNLMax; ++j) {
        if (!*p || *p == ',' || *p == ')' || *p == '[' || *p == '(' || *p == ']')
            break;
        *q++ = *p++;
    }
    *q = '\0';
    if (mat_ncheck(mname,opt))   
        return(NULL);
    return(p);
}

/*--------------------------------------------------------------------------*/
/*  mat_ncopy(row,col,x,name)                                               */
/*                                                                          */  
/*  Create a new (row,col) matrix name and copy x[] into this matrix.       */
/*  Return index of new matrix, or -1 if error.                             */

int mat_ncopy(int row,int col,double *x,char *name)
{
    register int i;
    int idx,n;
    double *y;

    if ((idx = mat_newmat(name,row,col)) < 0)
        return(-1);       

    y = MatVal[idx];
    n = row * col;
    for (i = 1; i <= n; ++i)  
        y[i] = x[i];
    return(idx);
}


