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
#include "tda_context.h"
#include "t_mat.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_mat                                                      */

int t_mat(TDAContext *ctx);
void mdefcpy(TDAContext *ctx, char *d,char *s);
void m_cmdmsg(TDAContext *ctx);
int alloc_local(TDAContext *ctx, char *cmd);
void free_local(TDAContext *ctx);
int check_local(TDAContext *ctx, char *name);
int alloc_mat(TDAContext *ctx);
void mat_free(TDAContext *ctx);
int mat_alloc(TDAContext *ctx, int idx,int opt);
int mat_getidx(TDAContext *ctx, char *mname,int opt);   
int mat_newidx(TDAContext *ctx, char *mname,int m,int n);
int mat_newmat(TDAContext *ctx, char *mname,int m,int n);
int mat_ncheck(TDAContext *ctx, char *p,int opt);
void mat_err(TDAContext *ctx, int opt);   
int mat_info(TDAContext *ctx);           
int mfree(TDAContext *ctx, char *def); 
int mdef(TDAContext *ctx, char *def); 
int m_mdefb(TDAContext *ctx, char *def); 
int mdeff(TDAContext *ctx, char *def); 
int mdefg(TDAContext *ctx, char *def); 
char *get_mname(TDAContext *ctx, char *p,char *mname,int opt);
int mat_ncopy(TDAContext *ctx, int row,int col,double *x,char *name);

/* ------------------------------------------------------------------------ */
/*  global variables.                                                       */




/*--------------------------------------------------------------------------*/
/*  t_mat()     Entry point for matrix commands.                            */
/*              Return 0 if OK, -1 if error, 1 if command cannot be         */
/*              interpreted.                                                */

int t_mat(TDAContext *ctx)
{
    register char *p;
    int err,n,m;           

    if (ctx->PMATFmtF == 0) {
        ctx->PMMFmt1 = 10;                     /* make default print format */
        ctx->PMMFmt2 = 4; 
        makefmt(ctx, &ctx->PMATFmt1,&ctx->PMATFmt2,ctx->PMATFmtS,sizeof(ctx->PMATFmtS),0,' ',0);
        ctx->PMATFmtF = 1;
    }
    err = 0; 
    p = ctx->CmdBuf;

    /* save command for error messages */

    n = (int)(strlen(p));
    if (!(ctx->MatCmdBuf = (char *)calloc((size_t)(n+1),sizeof(char)))) {
        p_err(ctx, -2,1);
        return(-1);
    }         
    memrq(ctx, n + 1,sizeof(char));
    ctx->MatCmdBufLen = n + 1;
    strcpy(ctx->MatCmdBuf,p);

    if (sscanf(p,"mfmt=%d.%d",&n,&m) == 2) {
        ctx->PMATFmt1 = n;       
        ctx->PMATFmt2 = m; 
        makefmt(ctx, &ctx->PMATFmt1,&ctx->PMATFmt2,ctx->PMATFmtS,sizeof(ctx->PMATFmtS),0,ctx->SEPC,0);
        printf2(ctx, "%s\n",p);
    }

    else if (!strncmp(p,"mcent(",6))        /* mcent(A,R) */
        err = m_cent(ctx, p,0);

    else if (!strncmp(p,"mdcent(",7))       /* mdcent(A,R) */
        err = m_dcent(ctx, p);

    else if (!strncmp(p,"mstand(",7))       /* mstand(A,R) */
        err = m_cent(ctx, p,1);

    else if (!strncmp(p,"mcross(",7))       /* mcross(ctx, X,R) */
        err = m_cross(ctx, p);

    else if (!strncmp(p,"mchol(",6))        /* mchol(X,R) */
        err = m_chol(ctx, p);

    else if (!strncmp(p,"mag(",4))          /* mag(A,C,R,B) */
        err = m_agg(ctx, p);

    else if (!strncmp(p,"mcvec(",6))        /* mcvec(A,C) */
        err = m_vec(ctx, p,0);

    else if (!strncmp(p,"mrvec(",6))        /* mrvec(A,V) */
        err = m_vec(ctx, p,1);

    else if (!strncmp(p,"mivec(",6))        /* mivec(V,n,A) */
        err = m_ivec(ctx, p);

    else if (!strncmp(p,"mcsum(",6))        /* mcsum(A,U) */
        err = m_sum(ctx, p,1);

    else if (!strncmp(p,"mrsum(",6))        /* mrsum(A,V) */
        err = m_sum(ctx, p,0);

    else if (!strncmp(p,"mdrow(",6))        /* mdrow(X,R) */
        err = m_drow(ctx, p,0);

    else if (!strncmp(p,"mdcol(",6))        /* mdcol(X,R) */
        err = m_drow(ctx, p,1);

    else if (!strncmp(p,"mdiag(",6))        /* mdiag(A,R) */
        err = m_diag(ctx, p,0);

    else if (!strncmp(p,"mdiagd(",7))       /* mdiagd(A,R) */
        err = m_diag(ctx, p,1);

    else if (!strncmp(p,"mtransp(",8))      /* mtransp(X,R) */
         err = m_transp(ctx, p);

    else if (!strncmp(p,"msqrtd(",7))       /* msqrtd(X,R) */
        err = m_sqrt(ctx, p,0);

    else if (!strncmp(p,"msqrti(",7))       /* msqrti(X,R) */
        err = m_sqrt(ctx, p,1);

    else if (!strncmp(p,"mnrow(",6))        /* mnrow(A,R) */
        err = m_nrow(ctx, p,0);

    else if (!strncmp(p,"mncol(",6))        /* mncol(A,R) */
        err = m_nrow(ctx, p,1);

    else if (!strncmp(p,"mnc(",4))          /* mnc(A,x,opt,B) */
        err = m_mnc(ctx, p);

    else if (!strncmp(p,"mnorm(",6))        /* mnorm(A,R) */
        err = m_nrow(ctx, p,2);

    else if (!strncmp(p,"mnorm1(",7))       /* mnorm1(A,R) */
        err = m_nrow(ctx, p,3);

    else if (!strncmp(p,"mnorm2(",7))       /* mnorm2(A,R) */
        err = m_nrow(ctx, p,4);

    else if (!strncmp(p,"mtrace(",7))       /* mtrace(A,R) */
        err = m_nrow(ctx, p,5);

    else if (!strncmp(p,"mnum(",5))         /* mnum(x,d,n,A) */
        err = m_num(ctx, p);

    else if (!strncmp(p,"mevs(",5))         /* mevs(A,E,EV) */
        err = m_mevs(ctx, p);

    else if (!strncmp(p,"mev(",4))          /* mev(A,ER,EI,EV) */
        err = m_mev(ctx, p);

    else if (!strncmp(p,"mginv(",6))        /* mginv(A,R) */
        err = m_ginv(ctx, p);

    else if (!strncmp(p,"msvd1(",6))        /* msvd(A,Q,U,V*/
        err = m_svd(ctx, p,1);

    else if (!strncmp(p,"msvd(",5))         /* msvd(A,Q) */
        err = m_svd(ctx, p,0);

    else if (!strncmp(p,"mwvec1(",7))       /* mwvec1(A,W,T,R) */
        err = m_wvec1(ctx, p);

    else if (!strncmp(p,"mwvec(",6))        /* mwvec(A,W,R) */
        err = m_wvec(ctx, p);

    else if (!strncmp(p,"mscal1(",7))       /* mscal1(A,B) */
        err = m_scal1(ctx, p);

    else if (!strncmp(p,"mdefb(",6))        /* mdefb(A,bn) */
        err = m_mdefb(ctx, p);

    else if (!strncmp(p,"mdefc(",6))        /* mdefc(m,n,d,A) */
        err = m_mdefc(ctx, p);

    else if (!strncmp(p,"mdefi(",6))        /* mdefi(m,n,A) */
        err = m_mdefi(ctx, p);

    else if (!strncmp(p,"minvs(",6))        /* minvs(A,R) */
        err = m_invs(ctx, p);

    else if (!strncmp(p,"minvd(",6))        /* minvd(A,R) */
        err = m_invd(ctx, p);

    else if (!strncmp(p,"mple(",5))         /* mple(T,C,F,D) */
        err = m_ple(ctx, p);

    else if (!strncmp(p,"mnvar(",6))        /* mnvar(X) */
        err = m_nvar(ctx, p);

    else if (!strncmp(p,"mpsym(",6))        /* mpsym(A,P,B) */
        err = m_mperm(ctx, p,0);

    else if (!strncmp(p,"mprow(",6))        /* mprow(A,P,B) */
        err = m_mperm(ctx, p,1);

    else if (!strncmp(p,"mpcol(",6))        /* mpcol(A,P,B) */
        err = m_mperm(ctx, p,2);

    else if (!strncmp(p,"mpr",3))           /* print matrix */
        err = m_print(ctx, p);  

    else if (!strncmp(p,"mkp(",4))          /* mkp(A,B,R) */
        err = m_mkp(ctx, p);

    else if (!strncmp(p,"mlsei1(",7))       /* mlsei1(S,me,mi,R) */
        err = m_mlsei1(ctx, p);

    else if (!strncmp(p,"mnls(",5))         /* mnls(S,l,R) */
        err = m_mls(ctx, p,4);

    else if (!strncmp(p,"mlse(",5))         /* mlse(S,R) */
        err = m_mls(ctx, p,1);

    else if (!strncmp(p,"mlsi(",5))         /* mlsi(S,R) */
        err = m_mls(ctx, p,2);

    else if (!strncmp(p,"mlsei(",6))        /* mlsei(S,me,mi,R) */
        err = m_mls(ctx, p,3);

    else if (!strncmp(p,"mls(",4))          /* mls(S,R) */
        err = m_mls(ctx, p,0);

    else if (!strncmp(p,"mqpc(",5))         /* mqpc(C,D,A,B,me,X) */
        err = m_mqp(ctx, p,2);

    else if (!strncmp(p,"mqpb(",5))         /* mqpb(C,D,XL,XU,X) */
        err = m_mqp(ctx, p,1);

    else if (!strncmp(p,"mqp(",4))          /* mqp(C,D,X) */
        err = m_mqp(ctx, p,0);

    else if (!strncmp(p,"mlpi(",5))         /* mlpi(A,B,X) */
        err = m_mlpi(ctx, p);

    else if (!strncmp(p,"mlp1(",5))         /* mlp1(T,p,X,Y) */
        err = m_mlp(ctx, p,1);

    else if (!strncmp(p,"mlp(",4))          /* mlp(T,X,Y) */
        err = m_mlp(ctx, p,0);

    else if (!strncmp(p,"mmul(",5))         /* mmul(A,...,R,) */
        err = m_mul(ctx, p);
  
    else if (!strncmp(p,"mcath(",6))        /* mcath(A,...,R) */
        err = m_cat(ctx, p,0);

    else if (!strncmp(p,"mcatv(",6))        /* mcatv(A,...,R) */
        err = m_cat(ctx, p,1);

    else if (!strncmp(p,"mcathv(",7))       /* mcathv(A,...,R) */
        err = m_cat(ctx, p,2);

    else if (!strncmp(p,"msrow(",6))        /* msrow(A,D,R) */
        err = m_srow(ctx, p,0);

    else if (!strncmp(p,"mscol(",6))        /* mscol(A,D,R) */
        err = m_srow(ctx, p,1);

    else if (!strncmp(p,"msort(",6))        /* msort(A,D,R) */
        err = m_sort(ctx, p,0);

    else if (!strncmp(p,"msort1(",7))       /* msort1(A,D,R) */
        err = m_sort(ctx, p,2);

    else if (!strncmp(p,"mrank(",6))        /* mrank(A,D,R) */
        err = m_sort(ctx, p,1);

    else if (!strncmp(p,"msetv(",6))        /* msetv command */
        err = m_setv(ctx, p);  

    else if (!strncmp(p,"mdeff(",6))        /* matrix definition */
        err = mdeff(ctx, p);  

    else if (!strncmp(p,"mdefg",5))         /* matrix definition */
        err = mdefg(ctx, p);  

    else if (!strncmp(p,"mdef(",5))         /* matrix definition */
        err = mdef(ctx, p);  

    else if (!strcmp(p,"mdef"))             /* info about matrices */
        err = mat_info(ctx);

    else if (!strncmp(p,"mfree",5))         /* mfree(ctx, X) */
        err = mfree(ctx, p);  

    else if (!strncmp(p,"mexpr1(",7))       /* mexpr1(B,expression,A) */
        err = m_expr1(ctx, p);

    else if (!strncmp(p,"mexpr(",6))        /* mexpr(expression,A) */
        err = m_expr(ctx, p);

    else if (!strncmp(p,"mexp(",5))         /* mexp(expression,A) */
        err = m_exp(ctx, p);

    else if (!strncmp(p,"mbrr(",5))         /* mbrr(ns,nu,A) */
        err = m_brr(ctx, p);

    else if (!strncmp(p,"midf(",5))         /* midf(XL,XU,DL,DU,DM) */
        err = m_midf(ctx, p);

    else if (!strncmp(p,"midf1(",6))        /* midf1(XL,XU,F) */
        err = m_midf1(ctx, p);

    else if (!strncmp(p,"midf2(",6))        /* midf2(XL,XU,F) */
        err = m_midf2(ctx, p);

    else if (!strncmp(p,"midf3(",6))        /* midf3(XL,XU,XL1,XU1) */
        err = m_midf3(ctx, p);

    else if (!strncmp(p,"mtrim(",6))        /* mtrim(A,ca,ra,cb,rb,R) */
        err = m_mtrim(ctx, p);

    else if (!strncmp(p,"mmp1(",5))         /* mmp1(X,T,P) */
        err = m_mmp(ctx, p,1);

    else if (!strncmp(p,"mmp2(",5))         /* mmp2(X,T,P) */
        err = m_mmp(ctx, p,2);

    else if (!strncmp(p,"mmp(",4))          /* mmp(X,P) */
        err = m_mmp(ctx, p,0);

    else if (!strncmp(p,"mpz(",4))          /* mpz(A,B,P) */
        err = m_mpz(ctx, p);

    else if (!strncmp(p,"mpbl(",5))         /* mpbl(A,B,P,N,U) */
        err = m_mpb(ctx, p,0);

    else if (!strncmp(p,"mpbu(",5))         /* mpbu(A,B,P,N,U) */
        err = m_mpb(ctx, p,1);

    else if (!strncmp(p,"mqap(",5))         /* mqap(F,D,C,P) */
        err = m_mqap(ctx, p);

    else if (!strncmp(p,"mpinv(",6))        /* mpinv(P,Q) */
        err = m_mpinv(ctx, p);

    else if (!strncmp(p,"mcel(",5))         /* mcel(A,x,L) */
        err = m_mcel(ctx, p);

    else if (!strncmp(p,"mpfit(",6))        /* mpfit(ctx, A,R,C,B) */
        err = m_mpfit(ctx, p);

    else if (!strncmp(p,"mch(",4))          /* mch(A,B) */
        err = m_mch(ctx, p);

    else if (!strncmp(p,"mpit1(",6))        /* mpit1(F,N,Z,n,R) */
        err = m_mpit(ctx, p,1);

    else if (!strncmp(p,"mpit(",5))         /* mpit(F,N,n,R) */
        err = m_mpit(ctx, p,0);

    else if (!strncmp(p,"mldes(",6))        /* mldes(Z,G,D) */
        err = m_mldes(ctx, p);

    else if (!strncmp(p,"mkmet(",6))        /* mkmet(A,D) */
        err = m_kmet(ctx, p);

    else
        err = 1;

    if (ctx->MatCmdBufLen > 0) {    
        free(ctx->MatCmdBuf);
        memrq(ctx, -ctx->MatCmdBufLen,sizeof(char));
        ctx->MatCmdBufLen = 0;            
    }
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  mdefcpy(d,s)   copy matrix definition from s to d, maximal MatDefLen    */

void mdefcpy(TDAContext *ctx, char *d,char *s)
{
    (void)ctx;        /* unused: the signature is shared */
    strncpy(d,s,MatDefLen);
    *(d + MatDefLen) = '\0';
}

/*--------------------------------------------------------------------------*/
/*  m_cmdmsg    Print Error in command: ...                                 */

void m_cmdmsg(TDAContext *ctx)
{
    if (ctx->SILENTFlg >= 0 && ctx->MatCmdBufLen > 0)
        printf1(ctx, "Error in command: %s\n",ctx->MatCmdBuf);
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

int alloc_local(TDAContext *ctx, char *cmd)
{
    int len;
    register char *p;

    printf2(ctx, "local(%s [macro level %d]\n",cmd,ctx->MacroExLevel);

    if (ctx->MacroExLevel < 1 || ctx->MacroExLevel >= ctx->MaxMat)
        return(0);

    len = (int)(strlen(cmd));
    p = cmd + len;
    while (*--p == ')')
        len--;
    *++p = '\0';
    if (len < 1)   
        return(0);

    if (!(ctx->MatLocDef[ctx->MacroExLevel] = (char *)calloc((size_t)(len+1),sizeof(char)))) {
        p_err(ctx, -2,1);
        return(-1);
    }         
    memrq(ctx, len + 1,sizeof(char));
    strcpy(ctx->MatLocDef[ctx->MacroExLevel],cmd);
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  free_local      If MacroExLevel > 0 and < MaxMat                        */
/*                  free local matrices.                                    */

void free_local(TDAContext *ctx)
{
    register int i;
          
    if (ctx->MacroExLevel > 0 && ctx->MacroExLevel < ctx->MaxMat) {
        for (i = 0; i < ctx->MaxMat; ++i) {
            if (ctx->MatAlloc[i] && ctx->MatLoc[i] == ctx->MacroExLevel)    
                mat_alloc(ctx, i,0);
        }
        if (ctx->MatLocDef[ctx->MacroExLevel] != NULL) {  
            i = (int)(strlen(ctx->MatLocDef[ctx->MacroExLevel]) + 1);
            free(ctx->MatLocDef[ctx->MacroExLevel]);
            memrq(ctx, -i,sizeof(char));
        }
    }
}

/*--------------------------------------------------------------------------*/
/*  check_local(A)  Return 1 if name A occurs in local list for current     */
/*                  MacroExLevel, otherwise return 0.                       */

int check_local(TDAContext *ctx, char *name)
{
    register char c,*p,*q;

    if (ctx->MacroExLevel < 1 || ctx->MacroExLevel >= ctx->MaxMat)
        return(0);

    p = ctx->MatLocDef[ctx->MacroExLevel];            
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

int alloc_mat(TDAContext *ctx)
{
    int i;

    if (!(ctx->MatAlloc = (int *)calloc((size_t)(ctx->MaxMat),sizeof(int))) ||
        !(ctx->MatLoc   = (int *)calloc((size_t)(ctx->MaxMat),sizeof(int))) ||
        !(ctx->MatRow   = (int *)calloc((size_t)(ctx->MaxMat),sizeof(int))) ||
        !(ctx->MatCol   = (int *)calloc((size_t)(ctx->MaxMat),sizeof(int))) ||
        !(ctx->MatLocDef = (char **)calloc((size_t)(ctx->MaxMat),sizeof(char *))) ||
        !(ctx->MatName  = (char **)calloc((size_t)(ctx->MaxMat),sizeof(char *))) ||
        !(ctx->MatDef   = (char **)calloc((size_t)(ctx->MaxMat),sizeof(char *))) ||
        !(ctx->MatVal   = (double **)calloc((size_t)(ctx->MaxMat),sizeof(double *))))   
        return(-1);
              
    memrq(ctx, ctx->MaxMat,4 * sizeof(int) + 3 * sizeof(char *) + sizeof(double *));

    for (i = 0; i < ctx->MaxMat; ++i) {
        ctx->MatLocDef[i] = NULL;
        if (!(ctx->MatName[i] = (char *)calloc(VNLMax + 1,sizeof(char))) ||
            !(ctx->MatDef[i] = (char *)calloc(MatDefLen + 1,sizeof(char))))  
            return(-1);
        *(ctx->MatDef[i] + MatDefLen) = '\0';
    }
    memrq(ctx, ctx->MaxMat * (VNLMax + MatDefLen + 2),sizeof(char));
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  mat_free()      free memory for all allocated matrices.                 */
   
void mat_free(TDAContext *ctx)
{
    int i;  
       
    for (i = 0; i < ctx->MaxMat; ++i) {
        if (ctx->MatAlloc[i])   
            mat_alloc(ctx, i,0);
    }
}

/*--------------------------------------------------------------------------*/
/*  mat_alloc(idx,opt)      allocate (opt = 1) or free (opt = 0) memory     */
/*                          for matrix with index idx.                      */
/*                          return 0 if OK, -1 if error.                    */

int mat_alloc(TDAContext *ctx, int idx,int opt)
{
    int i,mn;
        
    if (idx == ctx->MPLogIdx)  
        ctx->MPLogIdx = -1;
    else if (idx == ctx->MPParIdx)  
        ctx->MPParIdx = -1;
    else if (idx == ctx->MPCovIdx)  
        ctx->MPCovIdx = -1;
    else if (idx == ctx->MPGradIdx) {
        ctx->MPGradIdx = -1;
        ctx->MPGradRow = ctx->MPGradCol = 0;
    }
    else if (idx == ctx->MPResIdx) {
        ctx->MPResIdx = -1;
        ctx->MPResRow = ctx->MPGradCol = 0;
    }
    if (ctx->GD_TYP == 8) {                  /* check for graph definition */
        for (i = 0; i < ctx->GD_NG; ++i) {
            if (ctx->GD_EV[i] == idx) {
                gdd_free(ctx, 1);
                break;
            }
        }
    }

    mn = ctx->MatRow[idx] * ctx->MatCol[idx] + 1;
    if (opt) {

        if (ctx->MatAlloc[idx] || mn <= 1) {
            printfe(ctx, "MAT_ALLOC ERROR: mn = %d\n",mn);
            gerr_exit(ctx, 71);
        }
        if (!(ctx->MatVal[idx] = (double *)calloc((size_t)(mn),sizeof(double)))) {
            m_cmdmsg(ctx);
            printf1(ctx, "Error: insufficient memory for (%d,%d) matrix.\n",        
                                                   ctx->MatRow[idx],ctx->MatCol[idx]);
            return(-1);
        }
        ctx->MatAlloc[idx] = 1;

        if (ctx->MacroExLevel > 0 && ctx->MacroExLevel < ctx->MaxMat) {
            if (check_local(ctx, ctx->MatName[idx]))
                ctx->MatLoc[idx] = ctx->MacroExLevel;
        }
    }
    else {
        if (ctx->MatAlloc[idx] == 0)
            gerr_exit(ctx, 72);

        free((char *)ctx->MatVal[idx]);
        mn = -mn;
        ctx->MatAlloc[idx] = 0;
    }
    memrq(ctx, mn,sizeof(double));
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  mat_getidx(mname,opt)   return pointer to mname if defined,             */
/*                          otherwise -1.                                   */
/*                          if opt == 1 print Error: ... not defined        */

int mat_getidx(TDAContext *ctx, char *mname,int opt)
{
    register int i;

    for (i = 0; i < ctx->MaxMat; ++i) {
        if (ctx->MatAlloc[i] && !strcmp(mname,ctx->MatName[i])) {
            if (ctx->MacroExLevel > 0 && ctx->MacroExLevel < ctx->MaxMat && check_local(ctx, mname)) {
                if (ctx->MatLoc[i] == ctx->MacroExLevel)  
                    return(i);           
            }
            else  
                return(i);
        }
    }
    if (opt == 1) {
        m_cmdmsg(ctx);
        printf1(ctx, "%s not defined.\n",mname);
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

int mat_newidx(TDAContext *ctx, char *mname,int m,int n)
{
    register int i;

    if (check_local(ctx, mname) == 0) {
        i = ctx->VIFirst;
        while (i >= 0) {
            if (!strcmp(mname,ctx->VName[i])) {
                m_cmdmsg(ctx);
                printf1(ctx, "Error: name conflict with variable: %s\n",mname);
                return(-1);
            }
            i = ctx->VNxt[i];
        }
        if (nl_check(ctx, mname) >= 0) {
            m_cmdmsg(ctx);
            printf1(ctx, "Error: name conflict with namelist: %s\n",mname);
            return(-1);
        }
    }
    for (i = 0; i < ctx->MaxMat; ++i) {
        if (ctx->MatAlloc[i] == 0) {
            strcpy(ctx->MatName[i],mname);
            ctx->MatRow[i] = m;
            ctx->MatCol[i] = n;
            return(i);           
        }
    }
    m_cmdmsg(ctx);
    printf1(ctx, "Error: already allocated max number (%d) of matrices.\n",ctx->MaxMat);
    return(-1);
}

/*--------------------------------------------------------------------------*/
/*  mat_newmat(mname,m,n)   get a new matrix and allocate memory.           */
/*                          return index if OK, -1 if error.                */

int mat_newmat(TDAContext *ctx, char *mname,int m,int n)
{
    int i,idx;

    idx = mat_getidx(ctx, mname,0);  /* check if matrix already defined */
    if (idx >= 0) {
        if (ctx->GD_TYP == 8) {                   /* check if used for graph */
            for (i = 0; i < ctx->GD_NG; ++i) {
                if (ctx->GD_EV[i] == idx) {
                    printf1(ctx, "Cannot automatically overwrite matrix %s\n",ctx->MatName[idx]);
                    printf1(ctx, "This matrix is used for a graph definition.\n");
                    return(-1);
                }
            }
        }
        mat_alloc(ctx, idx,0);       /* free memory */
    }
    idx = mat_newidx(ctx, mname,m,n);      /* get new index */
    if (idx < 0)
        return(-1);      

    if (mat_alloc(ctx, idx,1)) {   /* allocate memory */
        m_cmdmsg(ctx);
        printf1(ctx, "Error: insufficient memory for new matrix.\n");
        return(-1);
    }
    return(idx);
}

/*--------------------------------------------------------------------------*/
/*  mat_ncheck(p,opt)   check whether p points to a matrix name.            */
/*                      return 0 if OK, otherwise -1.                       */

int mat_ncheck(TDAContext *ctx, char *p,int opt)
{
    int l;
       
    if (check_vname(ctx, p) == 0)                            
        goto NCErr;
          
    l = get_vnlen(ctx, p);
    if (l < 1 || l > VNLMax)              
        goto NCErr;
    return(0);
   
NCErr:
    if (opt) {
        m_cmdmsg(ctx);
        printf1(ctx, "Error in matrix name.\n");
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

void mat_err(TDAContext *ctx, int opt) 
{
    m_cmdmsg(ctx);

    if (opt == 0)
        printf1(ctx, "Syntax error (or unknown strings in matrix expression).\n");
    else if (opt == 1)
        printf1(ctx, "Error: square matrix required.\n");
    else if (opt == 2)
        printf1(ctx, "Error in matrix dimensions.\n");
    else if (opt == 3)
        printf1(ctx, "Error: insufficient memory.\n");
    else if (opt == 4)
        printf1(ctx, "Error: matrix names must be different.\n");
    else if (opt == 5)
        printf1(ctx, "Error: required equal number of rows.\n");
    else if (opt == 6)
        printf1(ctx, "Error in definition of indices.\n");
    else if (opt == 7)
        printf1(ctx, "Error: column vector required.\n");
    else if (opt == 8)
        printf1(ctx, "Error: scalar expression(s) required.\n");
    else if (opt == 9)
        printf1(ctx, "Error: in selection indices.\n");
    else if (opt == 10)
        printf1(ctx, "Error: in matrix indices.\n");
    else if (opt == 11)
        printf1(ctx, "Error: need positive interval width.\n");
    else if (opt == 12)
        printf1(ctx, "Error: command requires block mode.\n");
    else if (opt == 13)
        printf1(ctx, "Error: block number out of range (1 -- %d).\n",ctx->BNOC);
    else if (opt == 14)
        printf1(ctx, "Error in number of constraints.\n");
    else if (opt == 15)
        printf1(ctx, "Error in resulting matrix dimensions.\n");
    else if (opt == 16)
        printf1(ctx, "Error in index vector (column or row).\n");
    else if (opt == 17)
        printf1(ctx, "Error in permutation vector.\n");
    else if (opt == 18)
        printf1(ctx, "Error: inconsistent parameters.\n");

}

/*--------------------------------------------------------------------------*/
/*  mat_info()  print info about currently defined matrices.                */

int mat_info(TDAContext *ctx)            
{
    register int i,n;
    int maxlen = 16;

    n = 0;
    for (i = 0; i < ctx->MaxMat; ++i) {
        if (ctx->MatAlloc[i])
            n++;
    }
    if (n == 0) {
        printf1(ctx, "mdef: no matrices defined.\n");
        return(0);
    }
    printf1(ctx, "mdef\nmatrix            loc  rows  columns  definition\n");
    prnchar(ctx, '-',48,1);
                                 
    for (i = 0; i < ctx->MaxMat; ++i) {
        if (ctx->MatAlloc[i]) {
            printf1(ctx, "%s", ctx->MatName[i]);          
            prnchar(ctx, ' ',maxlen - (int)strlen(ctx->MatName[i]),0);
            printf1(ctx, "  %3d %5d %6d    %s\n",ctx->MatLoc[i],
                                 ctx->MatRow[i],ctx->MatCol[i],ctx->MatDef[i]);
        }
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  mfree(def)  def is pointer to string mfree(X). Free memory for          */
/*              matrix X. Or mfree without arguments free's all matrices.   */
/*              return 0 if OK, -1 if error.                                */
        
int mfree(TDAContext *ctx, char *def)  
{
    register char *p,*q;
    int idx; 

    printf2(ctx, "%s\n",def);
    if (!strcmp(def,"mfree")) {
        mat_free(ctx);
        return(0);
    }
    p = def + 5;
    if (*p++ != '(') {
        mat_err(ctx, 0);
        return(-1);   
    }                                                   
    q = p;
    if (!*q) {
        mat_err(ctx, 0);
        return(-1);   
    }
    while (*q && *q != ',' && *q != ')')
        q++;

    if (*q != ')') {
        mat_err(ctx, 0);
        return(-1);   
    }
    *q = '\0';

    idx = mat_getidx(ctx, p,1);
    if (idx < 0)
        return(-1);

    mat_alloc(ctx, idx,0);
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

int mdef(TDAContext *ctx, char *def)  
{
    FILE *fd;
    register int i,j,k,ii;
    register char *p,*q;
    int r,m,n,idx,mtyp,dflag;
    char mtmp[MatDefLen + 1];
    double tmp;

    mdefcpy(ctx, mtmp,def);
    printf2(ctx, "%s\n",mtmp);

    mtyp = 0;
    p = def + 4;
    dflag = 1;
    if (ctx->NVAR < 1 || ctx->NOC < 1 || !ctx->DMDef)  
        dflag = 0;

   
    if (*p++ != '(' || mat_ncheck(ctx, p,0)) {
        mat_err(ctx, 0);                   
        return(-1);
    }
    q = p + get_vnlen(ctx, p);             
    if (*q == ',') {
        *q++ = '\0';
        if (sscanf(q,"%d,%d",&m,&n) != 2 || m < 1 || n < 1) {
            mat_err(ctx, 0);    
            return(-1);
        }
        q = skip_int(ctx, q);
        q = skip_int(ctx, ++q);
        if (*q++ != ')') {
            mat_err(ctx, 0);       
            return(-1);
        }
    }
    else if (*q == ')') {
        *q++ = '\0';
        if (dflag == 0) {
            m_cmdmsg(ctx);
            printf1(ctx, "Error: no data matrix.\n");
            return(-1);
        }
        m = ctx->NOC;
        n = 0;
    }
    else {
        mat_err(ctx, 0);
        return(-1);
    }
    if (*q == '=') {
                     
        q++;
        if (dflag && check_vname(ctx, q)) {
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
            mat_err(ctx, 0);    
            return(-1);
        }
    }
    else {
        if (*q) {
            mat_err(ctx, 0);       
            return(-1);
        }
        if (dflag == 0) {
            printf1(ctx, "Error: no data matrix.\n");
            return(-1);
        }
        if (n > 0)
            mtyp = 0;
        else
            mtyp = 4;
    }

    switch (mtyp) {
        case 4:                             /* full data matrix */
            n = ctx->NVAR;
            TDA_FALLTHROUGH;
        case 0:                             /* part of data matrix */
            idx = mat_newmat(ctx, p,m,n);        /* get new matrix */
            if (idx < 0)
                return(-1);        
            ii = m; 
            if (ii > ctx->NOC) 
                ii = ctx->NOC;

            for (i = 0; i < ii; ++i) {
                k = 1;
                j = ctx->VIFirst;
                while (j >= 0) {
                    ctx->MatVal[idx][i * n + k] = get_data(ctx, j,i);
                    if (++k > n)
                        break;
                    j = ctx->VNxt[j];
                }
            }
            break;

        case 2:                     /* defined by numerical entries */

            idx = mat_newmat(ctx, p,m,n);     /* get new matrix */
            if (idx < 0)
                return(-1);        

            for (i = 1; i <= m; ++i) {
                ii = (i - 1) * n;
                for (j = 1; j <= n; ++j) {

                    if (sscanf(q,"%lg",&tmp) != 1) {
                        m_cmdmsg(ctx);
                        printf1(ctx, "Error: can't read entry (%d,%d).\n",i,j);
                        return(-1);
                    }
                    ctx->MatVal[idx][ii + j] = tmp;
                    q = skip_dval(ctx, q);
                    if (i == m && j == n)
                        break;
                    if (*q++ != ',') {
                        if (++j > n) {
                            ++i;
                            j = 1;
                        }
                        m_cmdmsg(ctx);  
                        printf1(ctx, "Syntax error in entry (%d,%d).\n",i,j);
                        return(-1);
                    }
                }
            }
            if (*q)  
                printf1(ctx, "Warning: data list contains more than %d entries.\n",m * n);

            break;

        case 3:                     /* defined by external file */

            idx = mat_newmat(ctx, p,m,n);         /* get new matrix */
            if (idx < 0)
                return(-1);        

            if (!(fd = fopen(q,"r"))) {
                m_cmdmsg(ctx);
                printf1(ctx, "Error: can't open data file: %s\n",q);
                return(-1);     
            }
            if (alloc_acc(ctx, RLMaxDef + 1)) {
                fclose(fd);
                return(-1);
            }
            i = 0;
            while (i < m) {

                if (fgets(ctx->AcC,RLMaxDef,fd) == NULL || *ctx->AcC == 0x1a)
                    break;

                if (check_drec(ctx, ctx->AcC)) {      /* check for data records */

                    if (++i > m)
                        break;
                    ii = (i - 1) * n;

                    p = ctx->AcC;

                    for (j = 1; j <= n; ++j) {    /* read required variables */

                        p = skip_sep(ctx, p);    /* skip separator */

                        if (sscanf(p,"%lg",&tmp) != 1) {
                            m_cmdmsg(ctx);
                            printf1(ctx, "Error: can't read entry (%d,%d).\n",i,j);
                            fclose(fd);
                            return(-1);
                        }
                        ctx->MatVal[idx][ii + j] = tmp;
                        p = skip_dval(ctx, p);
                    }
                }
            }
            fclose(fd);
            alloc_acc(ctx, 0);
            if (i < m)  
                printf1(ctx, "Warning: %s contains only %d data records.\n",q,i);
            break;

        case 1:                     /* based on varlist */
        case 5:     
            get_var(ctx, q,&r);
            if (r || ctx->PMNV < 1)
                return(-1);
            if (mtyp == 5)
                n = ctx->PMNV;

            idx = mat_newmat(ctx, p,m,n);     /* get new matrix */
            if (idx < 0)
                return(-1);        

            for (i = 0; i < m; ++i) {
                if (i >= ctx->NOC)
                    break;
                for (j = 0; j < ctx->PMNV; ++j) {
                    if (j >= n)
                        break;
                    ctx->MatVal[idx][i * n + j + 1] = get_data(ctx, (int)ctx->PMVIdx[j],i);
                }           
            }
            break;
    }
    strcpy(ctx->MatDef[idx],mtmp);
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  m_mdefb(def,bn) Definition of a matrix. def is pointer to a string:     */
/*  ##              The command creates a matrix corresponding to block     */
/*                  bn (may be number or matrix). Command requiers          */
/*                  block mode.                                             */
/*                                                                          */
/*              return 0 if OK, -1 if error.                                */

int m_mdefb(TDAContext *ctx, char *def)  
{
    register int i,j,k,l;
    register char *p,*q;
    int err,n,bn,idx,ii,row,col,iv;
    char mtmp[MatDefLen + 1];

    err = -1;

    mdefcpy(ctx, mtmp,def);
    printf2(ctx, "%s\n",mtmp);

    if (ctx->NVAR < 1 || ctx->NOC < 1 || !ctx->DMDef || ctx->BNOC <= 0) {
        mat_err(ctx, 12);
        goto MDEFBFin;
    }
    p = def + 5;
   
    if (*p++ != '(' || mat_ncheck(ctx, p,0)) {
        mat_err(ctx, 0);                   
        goto MDEFBFin;
    }
    q = p + get_vnlen(ctx, p);             

    if (*q != ',') {
        mat_err(ctx, 0);
        goto MDEFBFin;
    }
    *q++ = '\0';

    if ((q = n_mexpr(ctx, q,1,0,&row,&col,&iv,NULL)) == NULL)  
        goto MDEFBFin;  

    if (row != 1 || col != 1) {
        mat_err(ctx, 8);
        goto MDEFBFin;
    }
    bn = (int)get_mxval(ctx, 0,1,1);

    if (bn < 1 || bn > ctx->BNOC) {
        mat_err(ctx, 13);    
        goto MDEFBFin;
    }
    bn--;
    n = 0;
    ii = -1;
    for (i = 0; i < ctx->NOC; ++i) {
        if (ctx->DBlckPtr[i] == bn) {
            n++;
            if (ii < 0)
                ii = i;
        }
    }
    idx = mat_newmat(ctx, p,n,ctx->NVAR);    
    if (idx < 0)
        goto MDEFBFin;      
    
    l = 0;
    for (i = ii; i < ii + n; ++i) {
        k = 1;
        j = ctx->VIFirst;
        while (j >= 0) {
            ctx->MatVal[idx][l * ctx->NVAR + k++] = get_data(ctx, j,i);
            j = ctx->VNxt[j];
        }
        l++;
    }
    strcpy(ctx->MatDef[idx],mtmp);
    err = 0;   

MDEFBFin:
    mx_free(ctx);
    return(err);
}
  
/*--------------------------------------------------------------------------*/
/*  mdeff(def)  Definition of a matrix                                      */
/*              Syntax: mdeff(X) = fname.                                   */
/*              Create matrix X corresponding to data in fname.             */
/*                                                                          */
/*              return 0 if OK, -1 if error.                                */

int mdeff(TDAContext *ctx, char *def)  
{
    FILE *fd;
    int idx,m,n,i,j,ii;
    double tmp;
    char *p,*q,mname[VNLMax + 1];

    printf2(ctx, "%s\n",def);
    p = def + 6;

    if ((p = get_mname(ctx, p,mname,1)) == NULL)
        return(-1);                       
    if (*p++ != ')' || *p++ != '=' || !*p) {
        mat_err(ctx, 0);
        return(-1);
    }
    if (!(fd = fopen(p,"r"))) {
        m_cmdmsg(ctx);
        printf1(ctx, "Error: can't open data file: %s\n",p);
        return(-1);     
    }
    if (alloc_acc(ctx, RLMaxDef + 1)) {
        fclose(fd);
        return(-1);
    }
    n = m = 0;
    while (fgets(ctx->AcC,RLMaxDef,fd) && *ctx->AcC != 0x1a) {
        if (check_drec(ctx, ctx->AcC)) {      /* check for data records */

            if (m == 0) {
                q = ctx->AcC;
                while (1) {
                    q = skip_sep(ctx, q);    /* skip separator */

                    if (sscanf(q,"%lg",&tmp) != 1)  
                        break;

                    n++;
                    q = skip_dval(ctx, q);
                }
            }
            m++;
        }
    }
    if (m == 0 || n == 0) {
        m_cmdmsg(ctx);
        printf1(ctx, "Error: file %s does not contain readable data.\n",p);       
        return(-1);
    }
    if ((idx = mat_newmat(ctx, mname,m,n)) < 0)
        return(-1);       

    strcpy(ctx->MatDef[idx],def);

    fseek(fd,0,0);
    i = 0;
    while (i < m) {

        if (fgets(ctx->AcC,RLMaxDef,fd) == NULL || *ctx->AcC == 0x1a)
            break;

        if (check_drec(ctx, ctx->AcC)) {      /* check for data records */

            if (++i > m)
                break;
            ii = (i - 1) * n;

            q = ctx->AcC;

            for (j = 1; j <= n; ++j) {    /* read required variables */

                q = skip_sep(ctx, q);    /* skip separator */

                if (sscanf(q,"%lg",&tmp) != 1) {
                    m_cmdmsg(ctx);
                    printf1(ctx, "Error: can't read entry (%d,%d).\n",i,j);
                    fclose(fd);
                    alloc_acc(ctx, 0);
                    return(-1);
                }
                ctx->MatVal[idx][ii + j] = tmp;
                q = skip_dval(ctx, q);
            }
        }
    }
    fclose(fd);
    alloc_acc(ctx, 0);
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  mdefg(gn=...,sc=...) = Matrix name.                                     */
/*                                                                          */
/*  Create adjacency matrix for graph data, gn is graph number.             */
/*  Substitute missing values by sc.                                        */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int mdefg(TDAContext *ctx, char *def)  
{
    register int i,j;
    int idx,gn;           
    double sc,tmp;
    char *p,mname[VNLMax + 1];

    printf2(ctx, "%s\n",def);
    p = def + 5;

    gn = 1;
    sc = -1.0;
    if (*p == '(') {
        p++;
        while (*p) {
            if (sscanf(p,"gn=%d",&gn) == 1 && gn >= 0)
                p = skip_int(ctx, p + 3);
            else if (sscanf(p,"sc=%lg",&tmp) == 1) {
                sc = tmp;
                p = skip_dbl(ctx, p + 3);
            }
            else if (*p == ')')
                break;
            else {
                mat_err(ctx, 0);
                return(0);
            }
            if (*p == ',')  
                p++;
            else
                break;
        }
        if (*p++ != ')') {
            mat_err(ctx, 0);
            return(-1);
        }
    }
    if (*p++ != '=') {
        mat_err(ctx, 0);
        return(-1);
    }
    if ((p = get_mname(ctx, p,mname,1)) == NULL)
        return(-1);                       

    if (gdd_check(ctx, gn,1))
        return(-1);   

    if ((idx = mat_newmat(ctx, mname,ctx->GD_NP,ctx->GD_NP)) < 0)
        return(-1);       

    strcpy(ctx->MatDef[idx],def);

    for (i = 0; i < ctx->GD_NP; ++i) {
        for (j = 0; j < ctx->GD_NP; ++j) {
            tmp = gdd_adj(ctx, i,j,gn);
            if (tmp < 0.0)
                tmp = sc;
            ctx->MatVal[idx][i * ctx->GD_NP + j + 1] = tmp;
        }
    }
    return(0);
}
    
/*--------------------------------------------------------------------------*/
/*  get_mname(p,mname,opt)                                                  */
/*                  get mname from pointer p, return pointer to first       */
/*       character after mname. Print message and return NULL if err.      */

char *get_mname(TDAContext *ctx, char *p,char *mname,int opt)
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
    if (mat_ncheck(ctx, mname,opt))   
        return(NULL);
    return(p);
}

/*--------------------------------------------------------------------------*/
/*  mat_ncopy(row,col,x,name)                                               */
/*                                                                          */  
/*  Create a new (row,col) matrix name and copy x[] into this matrix.       */
/*  Return index of new matrix, or -1 if error.                             */

int mat_ncopy(TDAContext *ctx, int row,int col,double *x,char *name)
{
    register int i;
    int idx,n;
    double *y;

    if ((idx = mat_newmat(ctx, name,row,col)) < 0)
        return(-1);       

    y = ctx->MatVal[idx];
    n = row * col;
    for (i = 1; i <= n; ++i)  
        y[i] = x[i];
    return(idx);
}


