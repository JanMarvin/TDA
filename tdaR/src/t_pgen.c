/****************************************************************************/
/*  t_pgen                                                                  */
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
#include "t_parm.h"
#include "t_gdat.h"
#include "t_eval.h"
#include "t_eval1.h"
#include "t_var.h"
#include "t_cdf.h"
#include "t_qrmod.h"
#include "t_gf.h"
#include "tda_context.h"

/*  functions in t_pgen.c */

void p_err(TDAContext *ctx, int typ,int n);
void ps_err(TDAContext *ctx, int typ,char *txt,int n);
void p_warn(TDAContext *ctx, int typ,int n);
void p_serr(TDAContext *ctx, char *s,int n);
void p_wmsg(TDAContext *ctx, int typ,char *fname,int n);
void p_vterr(TDAContext *ctx, int typ);
int check_cmd(TDAContext *ctx, int opt);  
void p_cmd(TDAContext *ctx, char *pcmd,int nw);           
void prn_sve(TDAContext *ctx);
int eval_sve(TDAContext *ctx, int i);
int getd1(TDAContext *ctx, double *x,int ix,int ig,int opt);
int getd2(TDAContext *ctx, double *x,double *y,int ix,int iy,int opt);
int getdxy(TDAContext *ctx, double *x,double *y,int nx,int nif,int opt);
int check_nvar(TDAContext *ctx, int opt);
void prn_nwvar(TDAContext *ctx, int opt);
void prn_data(TDAContext *ctx, int m,int nc,int n,double *x,FILE *fd,char *fmt);
void prn1_data(TDAContext *ctx, double *y,double *x,int m,int ny,int nx,FILE *fd,char *fmt);
void prn1_coeff(TDAContext *ctx, int n,double *x,double *se,int ni,int df,short *vidx,int seflg);
void makefmt(TDAContext *ctx, int *fmt1,int *fmt2,char *fmts,size_t fmtsz,int mlen,char sepc,int opt);
void makenfmt(TDAContext *ctx, int *fmt,char *fmts,size_t fmtsz,int mlen,char sepc);
void pmfmt(TDAContext *ctx, int n,int m);
void pmtfmt(TDAContext *ctx, int n,int m);
void prn_cwt(TDAContext *ctx);        
void prn_sfmt(TDAContext *ctx, char *s,int n,char *fmt,double x);
void prn_f1mat(TDAContext *ctx, FILE *fd,int prn,char *fmt,int n,double *x);

/* ------------------------------------------------------------------------ */
/*  p_err(typ,n)    Print error message, add n newlines                     */

void p_err(TDAContext *ctx, int typ,int n)
{
    int c0 = ctx->ErrCnt;

    switch (typ) {
        case -1:    printf1(ctx, "Syntax error.");
                    break;
        case -2:    printf1(ctx, "Insufficient memory.");
                    break;
        case -3:    printf1(ctx, "Error: no relational/graph data defined.");
                    break;
        case -4:    printf1(ctx, "Error: syntax error or undefined variables.");
                    break;
        case -5:    printf1(ctx, "Error: exceeded max number of variables.");
                    break;
        case -6:    printf1(ctx, "Error: no variables.");
                    break;
        case -7:    printf1(ctx, "Error: can't open or read input file.");
                    break;
        case -8:    printf1(ctx, "Error: command needs at least two variables.");
                    break;
        case -9:    printf1(ctx, "Error: no data archive defined.");
                    break;
        case -10:   printf1(ctx, "Error: no PostScript output file.");
                    break;
        case -11:   printf1(ctx, "Error: exceeded max number of sequence data structures.");
                    break;
        case -12:   printf1(ctx, "Error: need sequence 0 data structure.");
                    break;
        case -13:   printf1(ctx, "Error: no sequence data defined.");
                    break;
        case -14:   printf1(ctx, "Error: exceeded max number of model parameters.");
                    break;
        case -15:   printf1(ctx, "Error: no episode data defined.");
                    break;
        case -16:   printf1(ctx, "Error: keep and drop options not compatible.");
                    break;
        case -17:   printf1(ctx, "Error: need time points/periods.");
                    break;
        case -18:   printf1(ctx, "Error: no data matrix (command ignored).");
                    break;
        case -19:   printf1(ctx, "Error: no patterns defined.");
                    break;
        case -20:   printf1(ctx, "Error in pattern definition.");
                    break;
        case -21:   printf1(ctx, "Pattern definition inconsistent with selected state space.");
                    break;
        case -22:   printf1(ctx, "Error: exceeded max number of string variables.");
                    break;
        case -23:   printf1(ctx, "Error: insufficient memory for string variables.");
                    break;
        case -24:   printf1(ctx, "Error: need an rc parameter.");
                    break;
        case -25:   printf1(ctx, "Error: need even number of points.");
                    break;
        case -26:   printf1(ctx, "Error: need at least two points.");
                    break;
        case -27:   printf1(ctx, "Error: need an output file.");
                    break;
        case -28:   printf1(ctx, "Error: selection results in zero cases.");
                    break;
        case -29:   printf1(ctx, "Error: must be only a single wave.");
                    break;
        case -30:   printf1(ctx, "Error: number of cases is less than number of variables.");
                    break;
        case -31:   printf1(ctx, "Error: need at least two groups.");
                    break;
        case -32:   printf1(ctx, "Syntax error in function definition.");
                    break;
        case -33:   printf1(ctx, "Error: time periods must begin with zero.");
                    break;
        case -34:   printf1(ctx, "Syntax error in model definition.");
                    break;
        case -35:   printf1(ctx, "Error: exceeded maximum number of y terms.");
                    break;
        case -36:   printf1(ctx, "Error: y terms must be recursive.");
                    break;
        case -37:   printf1(ctx, "Error: reached max level of command files.");
                    break;
        case -38:   printf1(ctx, "Error: need multi-episode data.");
                    break;
        case -39:   printf1(ctx, "Error in creating episode data.");
                    break;
        case -40:   printf1(ctx, "Error: cannot use type 5 variables (type 3 operators).");
                    break;
        case -41:   printf1(ctx, "Error: reference to cj terms not possible.");
                    break;
        case -42:   printf1(ctx, "Error: varlist must not contain brackets.");
                    break;
        case -43:   printf1(ctx, "Error: exceeded max length of command buffer.");
                    break;
        case -44:   printf1(ctx, "Error: in reading value labels.");
                    break;
        case -45:   printf1(ctx, "Error: in number of coordinates.");
                    break;
        default:    printf1(ctx, "Undefined error.");
                    break;
    }
    /* printf1() counts the messages that begin with "Error" or "Syntax
       error"; the rest -- "Insufficient memory.", "Undefined error." and
       the two pattern messages -- are counted here. */
    if (ctx->ErrCnt == c0)
        ++ctx->ErrCnt;
    while (n-- > 0)
        newline(ctx);
}

/* ------------------------------------------------------------------------ */
/*  ps_err(typ,txt,n)   Print error message, add n newlines                 */

void ps_err(TDAContext *ctx, int typ,char *txt,int n)
{
    switch (typ) {
        case -1:    printf1(ctx, "Error: can't open: %s",txt);
                    break;

        default:    printf1(ctx, "Undefined error.");
                    break;
    }
    while (n-- > 0)
        newline(ctx);
}

/* ------------------------------------------------------------------------ */
/*  p_serr(s,n)     Print error message, add n newlines                     */

void p_serr(TDAContext *ctx, char *s,int n)
{
    printf1(ctx, "Error: %s",s);
    while (n-- > 0)
        newline(ctx);
}

/* ------------------------------------------------------------------------ */
/*  p_warn(typ,n)   Print warning message, add n newlines                   */

void p_warn(TDAContext *ctx, int typ,int n)
{
    printf1(ctx, "Warning: ");
    switch (typ) {
        case -1:    printf1(ctx, "sel parameter ignored.");
                    break;
        case -2:    printf1(ctx, "procedure is incompatible with episode splitting.");
                    break;
        case -3:    printf1(ctx, "varlist contains type 5 variables.");
                    break;
        default:    printf1(ctx, "Warning ...");
                    break;
    }
    while (n-- > 0)
        newline(ctx);
}

/* ------------------------------------------------------------------------ */
/*  p_wmsg(typ,fname,n)   print write message, and add n newlines.          */

void p_wmsg(TDAContext *ctx, int typ,char *fname,int n)
{
    switch (typ) {
        case  1:    printf1(ctx, "Data written to");
                    break;
        case  2:    printf1(ctx, "Covariance matrix written to");
                    break;
        case  3:    printf1(ctx, "Estimated values written to");
                    break;
        default:    printf1(ctx, "... ");
                    break;
    }
    printf1(ctx, ": %s",fname);
    while (n-- > 0)
        newline(ctx);
}

/* ------------------------------------------------------------------------ */
/*  p_vterr(typ)    Print vtyp error.                                       */

void p_vterr(TDAContext *ctx, int typ)
{
    printf1(ctx, "Error: cannot use type %d variables.\n",typ);
}

/* ------------------------------------------------------------------------ */
/*  check_cmd(opt)      Print command in CmdBuf and check syntax.           */
/*  ##                  If opt = 1 print full command, if 0 abbreviated.    */  
/*                      If opt = 2 print cmd()...                           */
/*                      Return 0 if OK, -1 if error.                        */

int check_cmd(TDAContext *ctx, int opt)   
{
    int err = -1;
    register char *p = ctx->CmdBuf;

    if (opt == 1)
        printf1(ctx, "%s",p);

    while (*p && *p != '(' && *p != '=') {
        if (opt == 0 || opt == 2)
            printf1(ctx, "%c",*p);
        p++;
    }        
    if (*p == '(') {
        if (opt != 1)
            printf1(ctx, "(...");
        p = skip_blev(ctx, p);
        if (*(p - 1) != ')')  
            goto PRNCMDFin;
        if (opt != 1)
            printf1(ctx, ")");
    }
    if (*p) {
        if (opt == 0)
            printf1(ctx, "%s",p);
        else if (opt == 2)
            printf1(ctx, "=...");

        if (*p++ != '=' || !*p)
            goto PRNCMDFin;
    }
    printf1(ctx, "\n");
    err = 0;

PRNCMDFin:
    if (err) {
        printf1(ctx, "\n");
        printf1(ctx, "Syntax error.\n");
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  p_cmd(pcmd,nw)      print pcmd. if nw > 0 assume right-hand side is a   */
/*                      list of variables for nw panel waves.               */  
/*  ##                                                                      */

void p_cmd(TDAContext *ctx, char *pcmd,int nw)            
{
    register int i;
    int l,nf,n,r,nn;
    register char c,*p,*q,*pp,*qq;

    r = nf = 0;
    p = pcmd;
    printf1(ctx, "Command: ");           
    l = (int)(strlen(pcmd));
    if (l < 60)  
        goto PCMDFin;
           
    q = p;
    while (*q) {
        if (*q == '=' || *q == '(')
            break;
        q++;
    }
    if (!*q)
        goto PCMDFin;

    c = *++q;
    *q = '\0';
    printf1(ctx, "%s\n",p);        
    *q = c;

    while (1) {
        p = q;
        if (!strncmp(p,"xa(",3) || !strncmp(p,"xb(",3) || 
            !strncmp(p,"xc(",3) || !strncmp(p,"xd(",3))   
            q = skip_xa(ctx, p);
           
        else if (!strncmp(p,"fn=",3)) 
            q = skip_expr(ctx, p + 3);   
        else if (sscanf(p,"y%d",&n) == 1) {
            q = skip_int(ctx, p + 1);
            q = skip_expr(ctx, q + 1);
        }
        else if (!strncmp(p,"grp=",4)) {
            q = p + 3;
            while (*(q + 1) && (l = get_vnlen(ctx, q + 1)) > 0)  
                q += l + 1;
        }
        else {
            pp = p;
            while (*pp) {
                if (nw > 0) {
                    if (sscanf(pp,"nw=%d",&nn) == 1) {    
                        if (nw < nn)
                            nw = nn;
                    }   
                }
                q = skip_nc(ctx, pp);

                if (r > 0 && nw > 0) {
                    if (nw == 1) {
                        prnchar(ctx, ' ',11,0);
                        goto PCMDFin;
                    }
                    for (i = 1; i < nw; ++i)
                        q = skip_nc(ctx, q + 1);
                    break;
                }
                qq = q + 1;
                if (isdigit((int)*qq) || *qq == '-' || *qq == '+')
                    pp = q + 1;
   
                else if (*q == ')' && *qq == '=') {
                    q = qq;
                    break;
                }
                else if (*qq == ')' && *++qq == '=') {
                    q = qq;
                    break;
                }
                else
                    break;
            }
        }
        if (*q == ',') {
            *q = '\0';
            printf1(ctx, "           %s,\n",p);
            *q++ = ',';
            nf++;
        }
        else if (*q == '=' && strlen(q) > 50) {
            *q = '\0';
            printf1(ctx, "           %s=\n",p);
            *q++ = '=';
            nf++;
            r++;
        }
        else {   
            prnchar(ctx, ' ',11,0);
            goto PCMDFin;
        }
    }

PCMDFin:
    printf1(ctx, "%s\n",p);
    /* if (nf)  */
    newline(ctx);
}

/* ------------------------------------------------------------------------ */
/*  prn_sve()       print info about case selection.                        */

void prn_sve(TDAContext *ctx)
{
    if (ctx->SVEFlg)  
        printf1(ctx, "Case selection: %s\n",ctx->SVESTR);
}

/* ------------------------------------------------------------------------ */
/*  eval_sve(i)     evaluate sel expression for case i. return 1 if         */
/*                  case to be selected, otherwise zero.                    */

int eval_sve(TDAContext *ctx, int i)
{
    int r;
    double x;

    if (ctx->SVEFlg == 0)
        return(1);

    r = v_eval1(ctx, i,ctx->SVECnt,ctx->SVETyp,ctx->SVEVal,ctx->ESIdx,&x,0,0,0,0,0);             

    if (r) {
        printf1(ctx, "Error (%d): can't evaluate sel expression in case %d (ignored).\n",r,i+1);
        /** prn_emsg2(r); **/
        return(0);
    }
    if (fabs(x) > ctx->EPSI1)
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  getd1(x,ix,ig,opt)      get data for var ix into x[].                   */  
/*                          recognize select cases option.                  */
/*                          if ig >= 0 take this as a grouping variable     */
/*                          and use only cases where value is != 0.         */
/*                                                                          */
/*                          if opt print cases and error message.           */
/*                                                                          */
/*      Return number of cases, or -1 if error.                             */

int getd1(TDAContext *ctx, double *x,int ix,int ig,int opt)
{
    register int i,j;
    double tmp;

    j = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        if (eval_sve(ctx, i) == 0)   
            continue;

        if (ig >= 0) {
            tmp = get_data(ctx, ig,i);
            if (fabs(tmp) < ctx->EPSI1)
                continue;
        }
        x[j] = get_data(ctx, ix,i);       
        j++;  
    }
    if (opt) {
        if (ctx->SVEFlg)  
            printf1(ctx, "Number of selected cases: %d\n",j);
        if (j == 0)  
            p_err(ctx, -28,1);  
    }
    return(j);
}

/* ------------------------------------------------------------------------ */
/*  getd2(x,y,ix,iy,opt)    get data for var ix and iy into x[] and y[].    */  
/*                          if opt print cases and error message.           */

int getd2(TDAContext *ctx, double *x,double *y,int ix,int iy,int opt)
{
    register int i,j;

    j = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        if (eval_sve(ctx, i) == 0)   
            continue;
       
        x[j] = get_data(ctx, ix,i);       
        y[j] = get_data(ctx, iy,i);       
        j++;  
    }
    if (opt) {
        if (ctx->SVEFlg)  
            printf1(ctx, "Number of selected cases: %d\n",j);
        if (j == 0)  
            p_err(ctx, -28,1);  
    }
    return(j);
}

/* ------------------------------------------------------------------------ */
/*  getdxy(x,y,nx,nif,opt)                                                  */
/*                                                                          */  
/*  get data into y[] and x[] based on variables in PMVIdx in standard      */
/*  order. If nif == 0 add intercept. Number of columns in x[] are nx.      */
/*  Recognize sel option.                                                   */
/*  Return number of cases.                                                 */
/*  If opt print cases and error message.                                   */
/*                                                                          */
/*  Note: y[1,...], x(i,j) = x[(i - 1) * nx + j]                            */

int getdxy(TDAContext *ctx, double *x,double *y,int nx,int nif,int opt)
{
    register int i,j,k,jj;

    k = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        if (eval_sve(ctx, i) == 0)   
            continue;
       
        y[k + 1] = get_data(ctx, ctx->PMVIdx[0],i);       

          
        if (nif == 0) {
            x[k * nx + 1] = 1.0;
            jj = 2;
        }
        else
            jj = 1;

        for (j = 1; j < ctx->PMNV; ++j) {
            x[k * nx + jj] = get_data(ctx, ctx->PMVIdx[j],i);       
            jj++;
        }
        k++;  
    }
    if (opt) {
        if (ctx->SVEFlg)  
            printf1(ctx, "Number of selected cases: %d\n",k);
        if (k == 0)  
            p_err(ctx, -28,1);  
    }
    return(k);
}

/* ------------------------------------------------------------------------ */
/*  check_nvar(opt) check number of variables, PMNV, and waves, PMNW.       */
/*                  set NPX = number of x variables.                        */
/*                  if opt != 1 also check z variables and set NPZ.         */
/*                  return number of variables, or 0 if error.              */

int check_nvar(TDAContext *ctx, int opt)
{
    int n;

    ctx->NPX = ctx->NPZ = 0;
    n = ctx->PMNV / ctx->PMNW;
    if (ctx->PMNW * n != ctx->PMNV) {
        printf1(ctx, "Error: number of variables/waves is inconsistent.\n");
        return(0);
    }
    ctx->NPX = n - 1;

    if (opt && ctx->PMNZ > 0 && ctx->PMNQ > 0) {
        ctx->NPZ = ctx->PMNZ / (ctx->PMNQ * ctx->PMNW);
        if (ctx->NPZ * ctx->PMNQ * ctx->PMNW != ctx->PMNZ) {
            printf1(ctx, "Error: number of z-variables is inconsistent with nq=%d.\n",ctx->PMNQ);
            ctx->NPZ = 0;
            return(0);
        }
    }
    return(n);
}

/* ------------------------------------------------------------------------ */
/*  prn_nwvar(opt)  print variables for regression etc model.               */
/*  ##              if opt = 1  print z variables.                          */
/*                  if opt = 2  called by mvreg (YL,YH)                     */

void prn_nwvar(TDAContext *ctx, int opt)
{
    register int i,j,k,l;
    int n,iv,len,iz,nw;

    nw = ctx->PMNW;
    if (opt == 2)
        nw = 1;
    len = 0;
    if (opt && ctx->NPZ > 0)  
        len = 4;

    l = ctx->PMNVLEN + 1;
    n = ctx->PMNV / nw;
    printf1(ctx, "Variables");
    k = 9;
    if (nw > 1) {
        printf1(ctx, " (%d waves)\n",nw);
        k += 11;
    }
    else if (opt <= 1) {
        printf1(ctx, " (cross-section)\n");
        k += 16;
    }
    if (opt <= 1)
        prnchar(ctx, '-',k,1);
    else
        newline(ctx);

    if (opt == 2)
        printf1(ctx, "YL  ");
    else           
        printf1(ctx, "Y   ");

    prnchar(ctx, ' ',len,0);
    printf1(ctx, ": ");
    k = 0;
    for (i = 0; i < nw; ++i) {
        iv = ctx->PMVIdx[k++];
        printf1(ctx, "%s ",ctx->VName[iv]);
        prnchar(ctx, ' ',l - (int)strlen(ctx->VName[iv]),0);
    }
    newline(ctx);

    for (j = 1; j < n; ++j) {
        if (opt == 2) {
            if (j == 1)
                printf1(ctx, "YH  ");
            else
                printf1(ctx, "X%-3d",j - 1);
        }
        else
            printf1(ctx, "X%-3d",j);

        prnchar(ctx, ' ',len,0);
        printf1(ctx, ": ");
        for (i = 0; i < nw; ++i) {
            iv = ctx->PMVIdx[k++];
            printf1(ctx, "%s ",ctx->VName[iv]);
            prnchar(ctx, ' ',l - (int)strlen(ctx->VName[iv]),0);
        }
        newline(ctx);
    }
    if (opt == 1 && ctx->NPZ > 0) {
                 
        for (iz = 0; iz < ctx->NPZ; ++iz) {
            k = iz * nw * ctx->PMNQ;
            for (j = 0; j < ctx->PMNQ; ++j) {
                if (j == 0)
                    printf1(ctx, "Z%-3d %2d : ",iz + 1,j + 1);
                else
                    printf1(ctx, "   %4d : ",j + 1);

                for (i = 0; i < nw; ++i) {
                    iv = ctx->PMZIdx[k + i * ctx->PMNQ + j];
                    printf1(ctx, "%s ",ctx->VName[iv]);
                    prnchar(ctx, ' ',l - (int)strlen(ctx->VName[iv]),0);
                }
                newline(ctx);
            }
        }
    }
    newline(ctx);
}

/* ------------------------------------------------------------------------ */
/*  prn_data(m,nc,n,x,fd,fmt)                                               */
/*                                                                          */
/*  print n columns from m,nc matrix x[] to fd with print format fmt.       */
/*  note: x(i,j) = x[(i - 1) * nc + j].                                     */

#ifdef TDA_R_PACKAGE
/*  export_prn(name,m,nc,n,x)   The matrix prn_data() would write, laid out
    the same way (row stride nc, columns 1..n), handed to R under `name`
    whether or not a file was asked for.  The lsreg family's covariance
    matrix reaches R this way; the ml families have ml.vcov. */
void export_prn(TDAContext *ctx, const char *name, int m,int nc,int n,double *x)
{
    int i,j;

    for (i = 0; i < m; ++i) {
        for (j = 1; j <= n; ++j)
            tda_export_cell(ctx, name, x[i * nc + j]);
        tda_export_endrow(ctx, name);
    }
    tda_export_flush(ctx, name);
}
#endif

void prn_data(TDAContext *ctx, int m,int nc,int n,double *x,FILE *fd,char *fmt)
{
    register int i,j;

    for (i = 0; i < m; ++i) {
        for (j = 1; j <= n; ++j) {
            rt_fprintf_d(ctx, fd,fmt,x[i * nc + j]);
#ifdef TDA_R_PACKAGE
            tda_export_cell(ctx, "prn.data", x[i * nc + j]);
#endif
        }
        fprintf(fd,"\n");
#ifdef TDA_R_PACKAGE
        tda_export_endrow(ctx, "prn.data");
#endif
    }
#ifdef TDA_R_PACKAGE
    /* one flush per call, so a run writing several matrices through
       this printer gets prn.data, prn.data.2, ... beside them.  The
       lsreg family's covariance matrix comes through here (t_lsreg.c
       pcov=), where the ml families have ml.vcov instead. */
    tda_export_flush(ctx, "prn.data");
#endif
}

/* ------------------------------------------------------------------------ */
/*  prn1_data(y,x,m,ny,nx,fd,fmt)                                           */
/*                                                                          */
/*  printf first ny columns from y[], then nx columns from x[] to fd        */
/*  with format fmt. m = number of cases.                                   */
/*  note: y(i,j) = y[(i - 1) * ny + j], and same for x[].                   */

void prn1_data(TDAContext *ctx, double *y,double *x,int m,int ny,int nx,FILE *fd,char *fmt)
{
    register int i,j;

    for (i = 0; i < m; ++i) {
        for (j = 1; j <= ny; ++j)
            rt_fprintf_d(ctx, fd,fmt,y[i * ny + j]);
        for (j = 1; j <= nx; ++j)
            rt_fprintf_d(ctx, fd,fmt,x[i * nx + j]);
        fprintf(fd,"\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  prn1_coeff(n,x,se,ni,df,vidx,seflg)                                     */
/*                                                                          */
/*  print estimated regression coefficients. x[1,...,n] = coefficients,     */
/*  se[] contains standard errors. if ni = 1 without intercept. df =        */  
/*  degrees of freedom for t distribution. If seflg = 1, se[] contains      */
/*  squared standard errors, otherwise stand. errors.                       */
/*                                                                          */
/*  it is assumed that variables are in standard format in vidx[].          */  
/*                                                                          */
/*  if PMPPFDef print estimated parameters to PMPPFd.                       */

void prn1_coeff(TDAContext *ctx, int n,double *x,double *se,int ni,int df,short *vidx,int seflg)
{
#ifdef TDA_R_PACKAGE
    /* the same coefficient/standard-error pairs the table below prints,
       handed to R directly (CONTRIBUTING.md); the printed table is
       untouched and stays the parsers' source of truth for now */
    if (n > 0 && x != NULL && se != NULL) {
        double *ev = (double *)malloc((size_t)n * 4u * sizeof(double));
        if (ev != NULL) {
            int ei;
            for (ei = 1; ei <= n; ++ei) {
                double ese, erat;
                ev[ei - 1] = x[ei];
                /* mirror the printed columns exactly: Error is sqrt of
                   the variance under seflg, and Error / Coeff-over-
                   Error / Signif are NaN wherever the table prints ---
                   (fixed-by-constraint parameter, no df, or a
                   non-positive variance); Signif carries the full
                   double the text truncates to four decimals */
                if (df > 0 && !(ctx->ParFixed && ctx->ParFixed[ei]) &&
                    se[ei] > 0.0) {
                    ese = seflg ? sqrt(se[ei]) : se[ei];
                    erat = x[ei] / ese;
                    ev[n + ei - 1] = ese;
                    ev[2 * n + ei - 1] = erat;
                    ev[3 * n + ei - 1] =
                        2.0 * cdtf(ctx, fabs(erat), df) - 1.0;
                }
                else {
                    ev[n + ei - 1] = (double)NAN;
                    ev[2 * n + ei - 1] = (double)NAN;
                    ev[3 * n + ei - 1] = (double)NAN;
                }
            }
            tda_export_mat(ctx, "coeff", ev, n, 4);
            free(ev);
        }
    }
#endif
    register int i,j,k,l,ii;
    double tmp,e;

    l = ctx->PMNVLEN + 1;
    if (l < 9)
        l = 9;

    printf1(ctx, "\nIdx  Wave  Variable ");
    prnchar(ctx, ' ',l - 8 + ctx->PMTFmt1 - 5,0); printf1(ctx, "Coeff ");          
    prnchar(ctx, ' ',ctx->PMTFmt1 - 5,0); printf1(ctx, "Error ");          
    prnchar(ctx, ' ',ctx->PMTFmt1 - 7,0); printf1(ctx, "Coeff/E  Signif\n");           
    prnchar(ctx, '-',19 + l + 3 * (ctx->PMTFmt1 + 1),1);

    k = 1;
    i = ctx->PMNW;

    for (j = 1; j <= n; ++j) {
        if (j == 1 && ni == 0) {
            printf1(ctx, "%3d     -  Intercept ",k);            
            prnchar(ctx, ' ',l - 9,0);
#ifdef TDA_R_PACKAGE
            tda_export_str_row(ctx, "coeff.names", "-|Intercept");
#endif
        }
        else {
            printf1(ctx, "%3d     1  %s ",k,ctx->VName[vidx[i]]);           
            prnchar(ctx, ' ',l - (int)strlen(ctx->VName[vidx[i]]),0);
#ifdef TDA_R_PACKAGE
            {
                char b[96];
                snprintf(b, sizeof(b), "1|%s", ctx->VName[vidx[i]]);
                tda_export_str_row(ctx, "coeff.names", b);
            }
#endif
        }
        rt_printf1_d(ctx, ctx->PMTFmtS,x[j]);            
                              
        /* PATCH (tdaR, 2026-08, v3): se[j] > ctx->EPSI1 rejected any
           coefficient's standard error whenever its raw *variance*
           fell below a fixed absolute floor, printing "---" instead --
           correct for a genuinely singular fit (caught separately,
           upstream, via ranka != nx1-ranke, unaffected by this
           change), wrong whenever a predictor's own large natural
           scale makes its variance legitimately tiny
           (examples/exam/lsreg3.cf's own Y, variance 8.9e-10) or, more
           to the point, whenever it is simply small
           (lsreg3.cf's own Y2, variance 6.3e-19 -- real, matching R's
           lm() exactly, not a magnitude any fixed floor can safely let
           through while also rejecting a constrained parameter's own
           noise).
           No magnitude threshold can tell those two situations apart
           by size alone (se[j] > 0.0 and se[j] > EPSI each fail one of the two cases). The only other real reason se[j]
           could land near zero is a parameter solely, directly fixed
           by one lsecon=/lsicon= (con=b3=1, not a multi-parameter
           con=b1+b2=1, where neither parameter alone has zero
           variance) -- tracked directly instead, in ctx->ParFixed[]
           (set by p_con() in t_con.c at the point the constraint
           itself is parsed, before it reaches the solver), so the
           check here no longer needs to guess from se[j]'s own
           magnitude at all. */
        if (df > 0 && !(ctx->ParFixed && ctx->ParFixed[j]) && se[j] > 0.0) {
            if (seflg)
                e = sqrt(se[j]);
            else
                e = se[j];

            rt_printf1_d(ctx, ctx->PMTFmtS,e);              
            tmp = x[j] / e;
            rt_printf1_d(ctx, ctx->PMTFmtS,tmp);          
            tmp = 2.0 * cdtf(ctx, fabs(tmp),df) - 1.0;
            printf1(ctx, "%7.4lf",tmp);          
        }
        else {
            prnchar(ctx, ' ',ctx->PMTFmt1 - 3,0); printf1(ctx, "--- ");          
            prnchar(ctx, ' ',ctx->PMTFmt1 - 3,0); printf1(ctx, "---     ---");          
        }
        newline(ctx);
        if (ni == 1 || j > 1) {
            for (ii = 2; ii <= ctx->PMNW; ++ii) {
                printf1(ctx, "      %3d  %s ",ii,ctx->VName[vidx[++i]]);          
#ifdef TDA_R_PACKAGE
                /* one entry per printed continuation line -- the wave-2+
                   variable names of the same coefficient; captured
                   before the pad's own ++i below, mirroring the printed
                   name, and never re-deriving TDA's own index walk */
                {
                    char b[96];
                    snprintf(b, sizeof(b), "%d|%s", ii,
                             ctx->VName[vidx[i]]);
                    tda_export_str_row(ctx, "coeff.names", b);
                }
#endif
                prnchar(ctx, ' ',l - (int)strlen(ctx->VName[vidx[++i]]),0);
                newline(ctx);
            }
            i++;
        }
        k++;
    }
#ifdef TDA_R_PACKAGE
    /* one label per printed name line; with a single wave (the usual
       case) this pairs 1:1 with the coeff matrix's rows */
    tda_export_str_flush(ctx, "coeff.names");
#endif
    if (ctx->PMPPFDef) {         /* print estimated parameters */

        for (j = 1; j <= n; ++j) {
            rt_fprintf_d(ctx, ctx->PMPPFd,ctx->PMTFmtS,x[j]);            
            /* Same reasoning as the console print above: a solely
               constrained parameter (ctx->ParFixed[j]) writes nothing
               here either, rather than whatever noise its own
               covariance computation happened to leave. */
            if (df > 0 && !(ctx->ParFixed && ctx->ParFixed[j]) && se[j] > 0.0) {
                if (seflg)
                    e = sqrt(se[j]);
                else
                    e = se[j];

                rt_fprintf_d(ctx, ctx->PMPPFd,ctx->PMTFmtS,e);              
            }
            fprintf(ctx->PMPPFd,"\n");
        }
        printf1(ctx, "\nParameter estimates written to: %s\n",ctx->PMPPFName);
    }
}

/* ------------------------------------------------------------------------ */
/*  makefmt(fmt1,fmt2,fmts,mlen,sepc,opt)                                   */
/*      Make format string with fieldlength fmt1, and precision fmt2.       */
/*      Put string in fmts. mlen is the minimum field length, sepc is       */
/*      the separation character at the end of the format string.           */
/*      If opt == 1 special format for plot labels.                         */
/*      if opt == 2 string variables                                        */
/*      if opt == -1 then format with leading zeros                         */

void makefmt(TDAContext *ctx, int *fmt1,int *fmt2,char *fmts,size_t fmtsz,int mlen,char sepc,int opt)
{
    if (opt == 2) {
        snprintf(fmts,fmtsz,"%%%ds%c",imax(ctx, iabs(ctx, *fmt1),mlen),sepc);
        return;
    }
    if (mlen && *fmt1 > -mlen) {          
        if (*fmt1 < 0)          
            *fmt1 = -mlen;        
        else if (*fmt1 < mlen)
            *fmt1 = mlen;
    }
    if (*fmt2 < 0)
        *fmt2 = 0;
    else if (*fmt2 >= 100)
        *fmt2 = 99;

    if (!*fmt1) {
        if (opt == 0)
            snprintf(fmts,fmtsz,"%%lg%c",sepc);
        else if (opt < 0)
            snprintf(fmts,fmtsz,"%%0lg%c",sepc);
        else
            snprintf(fmts,fmtsz,"(%%g)%c",sepc);
    }
    else if (*fmt1 > 0) {
        if (*fmt1 >= 100)
            *fmt1 = 99;
        if (*fmt2 > *fmt1)
            *fmt2 = *fmt1;

        if (opt == 0)
            snprintf(fmts,fmtsz,"%%%d.%df%c",*fmt1,*fmt2,sepc);
        else if (opt < 0)
            snprintf(fmts,fmtsz,"%%0%d.%df%c",*fmt1,*fmt2,sepc);
        else
            snprintf(fmts,fmtsz,"(%%%d.%df)%c",*fmt1,*fmt2,sepc);
    }
    else {
        *fmt1 = -*fmt1;
        if (*fmt1 >= 100)
            *fmt1 = 99;
        if (*fmt1 - *fmt2 < 8)
            *fmt1 = *fmt2 + 8;
        if (opt == 0)
            snprintf(fmts,fmtsz,"%%%d.%de%c",*fmt1,*fmt2,sepc);
        else if (opt < 0)
            snprintf(fmts,fmtsz,"%%0%d.%de%c",*fmt1,*fmt2,sepc);
        else
            snprintf(fmts,fmtsz,"(%%%d.%de)%c",*fmt1,*fmt2,sepc);
    }
}

/* ------------------------------------------------------------------------ */
/*  makenfmt(fmt,fmts,mlen,sepc)                                            */
/*      Make integer format string with fieldlength fmt.                    */
/*      Put string in fmts. mlen is the minimum field length, sepc is       */
/*      the separation character at the end of the format string.           */

void makenfmt(TDAContext *ctx, int *fmt,char *fmts,size_t fmtsz,int mlen,char sepc)
{
    (void)ctx;        /* unused: the signature is shared */
    if (*fmt < mlen)
        *fmt = mlen;
    else if (*fmt >= 100)
        *fmt = 99;

    snprintf(fmts,fmtsz,"%%%dd%c",*fmt,sepc);
}
      
/* ------------------------------------------------------------------------ */
/*  pmfmt(n,m)      make a new PMFmtS.                                      */

void pmfmt(TDAContext *ctx, int n,int m)
{
    ctx->PMFmt1 = n;
    ctx->PMFmt2 = m;  
    makefmt(ctx, &ctx->PMFmt1,&ctx->PMFmt2,ctx->PMFmtS,sizeof(ctx->PMFmtS),0,ctx->SEPC,0);
}

/* ------------------------------------------------------------------------ */
/*  pmtfmt(n,m)      make a new PMTFmtS.                                    */

void pmtfmt(TDAContext *ctx, int n,int m)
{
    ctx->PMTFmt1 = n;
    ctx->PMTFmt2 = m;  
    makefmt(ctx, &ctx->PMTFmt1,&ctx->PMTFmt2,ctx->PMTFmtS,sizeof(ctx->PMTFmtS),0,ctx->SEPC,0);
}

/* ------------------------------------------------------------------------ */
/*  prn_cwt()   Print case weight information.                              */

void prn_cwt(TDAContext *ctx)         
{
    if (ctx->WIVar >= 0) 
        printf1(ctx, "Using case weights defined by: %s\n",ctx->VName[ctx->WIVar]);
}

/* ------------------------------------------------------------------------ */
/*  prn_sfmt(s,n,fmt,x)     print string s, then value x with fmt.          */
/*                          if n > 0 adjust length.                         */

void prn_sfmt(TDAContext *ctx, char *s,int n,char *fmt,double x)
{
    int l;

    printf1(ctx, "%s ",s);
    if (n > 0) {
        l = (int)(strlen(s));
        if (n > l)
            prnchar(ctx, ' ',n - l,0);
    }
    rt_printf1_d(ctx, fmt,x);
    printf1(ctx, "\n");
#ifdef TDA_R_PACKAGE
    tda_export_row(ctx, "sfmt.values", &x, 1);
    tda_export_str_row(ctx, "sfmt.labels", s);
#endif
}

/* ------------------------------------------------------------------------ */
/*  prn_f1mat(fd,prn,fmt,n,x)   write lower triangle of x matrix to fd.     */
/*                                                                          */
/*                  prn= 0 : lower triangle, including main diagonal        */
/*                  prn= 1 : full square matrix                             */
/*                  prn= 2 : column vector: lower triangle                  */
/*                  prn= 3 : column vector: lower triangle incl. diag.      */
/*                  prn= 4 : column vector: full matrix                     */
/*                                                                          */

void prn_f1mat(TDAContext *ctx, FILE *fd,int prn,char *fmt,int n,double *x)
{
    register int i,j,k;

    if (prn >= 1 && prn <= 4) {

        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= n; ++j) {
                if (prn == 2 && j >= i)
                    break;
                else if (prn == 3 && j > i)
                     break;

                if (j <= i)
                    k = i * (i - 1) / 2 + j - 1;
                else
                    k = j * (j - 1) / 2 + i - 1;

                rt_fprintf_d(ctx, fd,fmt,x[k]);
                if (prn != 1)
                    fprintf(fd,"\n");
            }
            if (prn == 1)
                fprintf(fd,"\n");
        }
    }
    else {
        k = 0;
        for (i = 0; i < n; ++i) {
            for (j = 0; j <= i; ++j)  
                rt_fprintf_d(ctx, fd,fmt,x[k++]); 
            fprintf(fd,"\n");
        }
    }
}



