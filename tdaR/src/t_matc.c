/****************************************************************************/
/*  t_matc                                                                  */
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
#include "t_svd.h"
#include "t_lin.h"
#include "t_lsei.h"
#include "t_lp.h"
#include "t_gf.h"
#include "t_alloc.h"
#include "t_mat.h"
#include "t_mata.h"
#include "t_matf.h"
#include "t_sort.h"
#include "t_eval.h"
#include "t_brr.h"
#include "t_mes.h"
#include "t_eval4.h"
#include "t_ass.h"
#include "t_qp.h"
#include "tda_context.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_matc                                                     */

char *n_mexpr(TDAContext *ctx, char *p,int n,int idx,int *rmax,int *cmax,int *ivflg,char *sel);
int sel_expr(TDAContext *ctx, char *sel,int row,int col,int *rsel,int *csel,int *idx,double *dval);
void mx_free(TDAContext *ctx);
int mx_alloc(TDAContext *ctx, int idx,int row,int col);
int mx1_alloc(TDAContext *ctx, int idx,int row,int col);
double get_mxval(TDAContext *ctx, int idx,int row,int col);
double get_matval(TDAContext *ctx, int idx,int row,int col);
void m_copy(TDAContext *ctx, double *x,double *y,int n);
char *m_check1(TDAContext *ctx, char *p);                  
char *m_check2(TDAContext *ctx, char *p,int opt);                  
char *m_getmat(TDAContext *ctx, char *p,int row,int col,int *idx,char *cmd,int opt);
int m_cent(TDAContext *ctx, char *cmd,int opt);
int m_dcent(TDAContext *ctx, char *cmd);
int m_cross(TDAContext *ctx, char *cmd);
int m_chol(TDAContext *ctx, char *cmd);
int m_agg(TDAContext *ctx, char *cmd); 
int m_vec(TDAContext *ctx, char *cmd,int opt);
int m_ivec(TDAContext *ctx, char *cmd);
int m_mldes(TDAContext *ctx, char *cmd);
int m_sum(TDAContext *ctx, char *cmd,int opt);
int m_drow(TDAContext *ctx, char *cmd,int opt);
int m_diag(TDAContext *ctx, char *cmd,int opt);
int m_transp(TDAContext *ctx, char *cmd);
int m_sqrt(TDAContext *ctx, char *cmd,int opt);
int m_nrow(TDAContext *ctx, char *cmd,int opt);
int m_num(TDAContext *ctx, char *cmd);
int m_mch(TDAContext *ctx, char *cmd);
int m_mpfit(TDAContext *ctx, char *cmd);
int m_mpinv(TDAContext *ctx, char *cmd);
int m_mperm(TDAContext *ctx, char *cmd,int opt);
int m_mqap(TDAContext *ctx, char *cmd);
int m_mcel(TDAContext *ctx, char *cmd);
int m_mnc(TDAContext *ctx, char *cmd);
int m_mpz(TDAContext *ctx, char *cmd);
int m_mpb(TDAContext *ctx, char *cmd,int opt);
int m_mevs(TDAContext *ctx, char *cmd);
int m_mev(TDAContext *ctx, char *cmd);
int m_ginv(TDAContext *ctx, char *cmd);
int m_svd(TDAContext *ctx, char *cmd,int opt);
int m_wvec(TDAContext *ctx, char *cmd);
int m_wvec1(TDAContext *ctx, char *cmd);
int m_scal1(TDAContext *ctx, char *cmd);
int m_mdefc(TDAContext *ctx, char *cmd);
int m_mdefi(TDAContext *ctx, char *cmd);
int m_invs(TDAContext *ctx, char *cmd);
int m_invd(TDAContext *ctx, char *cmd);
int m_ple(TDAContext *ctx, char *cmd);
int m_nvar(TDAContext *ctx, char *cmd);
int m_print(TDAContext *ctx, char *cmd); 
int m_mkp(TDAContext *ctx, char *cmd);
int m_mls(TDAContext *ctx, char *cmd,int opt);
int m_mlp(TDAContext *ctx, char *cmd,int opt);
int m_mlpi(TDAContext *ctx, char *cmd);
int m_mul(TDAContext *ctx, char *cmd);
int m_cat(TDAContext *ctx, char *cmd,int opt);
int m_srow(TDAContext *ctx, char *cmd,int opt);
char *m_selidx(TDAContext *ctx, char *s,int opt,int *n,int *idx);
int m_sort(TDAContext *ctx, char *cmd,int opt);
char *get_mselect(TDAContext *ctx, char *p);
int m_exp(TDAContext *ctx, char *cmd);
int m_expr(TDAContext *ctx, char *cmd);
int m_expr1(TDAContext *ctx, char *cmd);
int m_setv(TDAContext *ctx, char *cmd); 
int m_brr(TDAContext *ctx, char *cmd); 
int m_mtrim(TDAContext *ctx, char *cmd); 
int m_mmp(TDAContext *ctx, char *cmd,int opt);
int m_mpit(TDAContext *ctx, char *cmd,int opt);
int m_mqp(TDAContext *ctx, char *cmd,int opt);
int m_mlsei1(TDAContext *ctx, char *cmd);
int m_kmet(TDAContext *ctx, char *cmd);

/* ------------------------------------------------------------------------ */
/*  gobal variables                                                         */

#define MXN 7

#define MXMaxLen 60

/*--------------------------------------------------------------------------*/
/*  n_mexpr(p,n,idx,rmax,cmax,ivflg,sel)                                    */
/*                                                                          */
/*  Get n matrix-expressions beginning at pointer p. Save values in         */
/*  MX[idx] (idx = 0,...,n-1). Return pointer to next character after       */
/*  last expression. Return NULL if error.                                  */
/*  Return max of rows in rmax, max of cols in cmax.                        */
/*                                                                          */
/*  If, at input, ivflg = 1 allow for intervals. Set                        */
/*  ivflg=1 if at least one expression contains an interval, otherwise set  */
/*  ivflg=0. Also set MXIV[]. Lower bounds in MX, upper bounds in MX1.      */
/*                                                                          */
/*  If sel != NULL interpret as string <rsel,csel,val> for selection of     */
/*  rows and columns in evaluating the matrix expression.                   */
/*                                                                          */
/*  Save name of expression in MXName                                       */

char *n_mexpr(TDAContext *ctx, char *p,int n,int idx,int *rmax,int *cmax,int *ivflg,char *sel)
{
    register int ii,i;
    int err,row,col,m,iv1flg,iv2flg,rsela,csela,sidx = 0;
    int *rsel,*csel;
    register char c,*q;
    double dval = 0.0; 

    rsela = csela = 0; 
    err = -1;
    *rmax = *cmax = 1;
    iv2flg = 0;

    for (ii = 0; ii < n; ++ii) {

        iv1flg = *ivflg;

        if (!*p) {
            mat_err(ctx, 0);
            break;
        }
        q = skip_expr(ctx, p);
        c = *q;
        *q = '\0';

        if (get_mexpr(ctx, p,&row,&col,&iv1flg) || row < 1 || col < 1)
             goto N_MEXPRFin;

        ctx->MXIV[idx + ii] = 0;

        strncpy(ctx->MXName,p,MXMaxLen);
        *(ctx->MXName + MXMaxLen) = '\0';

        *rmax = imax(ctx, *rmax,row);
        *cmax = imax(ctx, *cmax,col);

        rsel = csel = NULL;
 
        if (sel != NULL) {              /* get selection strings */
            if (!(rsel = (int *)calloc((size_t)(row  + 1),sizeof(int)))) {
                p_err(ctx, -2,1);
                goto N_MEXPRFin;
            }
            rsela = row + 1;
            memrq(ctx, rsela,sizeof(int));

            if (!(csel = (int *)calloc((size_t)(col  + 1),sizeof(int)))) {
                p_err(ctx, -2,1);
                goto N_MEXPRFin;
            }
            csela = col + 1;
            memrq(ctx, csela,sizeof(int));

            if (sel_expr(ctx, sel,row,col,rsel,csel,&sidx,&dval))
                goto N_MEXPRFin;
        }

        m = iabs(ctx, ctx->ESTyp[0]);
        if (ctx->ESCnt == 1 && m >= ctx->MOFFS && m < ctx->VOFFS) {
            m -= ctx->MOFFS;
        }
        else {
            m = -1;
            if (alloc_actmp(ctx, row * col + 1))
                goto N_MEXPRFin;

            if (iv1flg) {
                if (alloc_actmp1(ctx, row * col + 1))
                    goto N_MEXPRFin;
            }
            if (eval_mexpr(ctx, 0,row,col,ctx->AcTmp,iv1flg,ctx->AcTmp1,rsel,csel,sidx,dval))
                goto N_MEXPRFin;
        }
        if (mx_alloc(ctx, idx + ii,row,col))
            goto N_MEXPRFin;

        if (m >= 0) {
            m_copy(ctx, ctx->MX[idx + ii],ctx->MatVal[m],row * col);
            if (ctx->ESTyp[0] < 0) {
                for (i = 1; i <= row * col; ++i)
                    ctx->MX[idx + ii][i] *= -1;
            }
        }
        else {
            m_copy(ctx, ctx->MX[idx + ii],ctx->AcTmp,row * col);

            if (iv1flg) {
                if (mx1_alloc(ctx, idx + ii,row,col))
                    goto N_MEXPRFin;

                m_copy(ctx, ctx->MX1[idx + ii],ctx->AcTmp1,row * col);
                ctx->MXIV[idx + ii] = iv2flg = 1;
            }
        }

        *q = c;
        p = q;  
        if (ii == n - 1)    
            break;
        if (*p++ != ',' || !*p) {
            mat_err(ctx, 0);
            break;
        }
    }
    *ivflg = iv2flg;
    err = 0;

N_MEXPRFin:
    if (rsela > 0) {
        free((char *)rsel);
        memrq(ctx, -rsela,sizeof(int));
    }
    if (csela > 0) {
        free((char *)csel);
        memrq(ctx, -csela,sizeof(int));
    }
    alloc_actmp(ctx, 0);
    alloc_actmp1(ctx, 0);
    if (err) {
        mx_free(ctx);
        p = NULL;
    }
    return(p);
}

/*--------------------------------------------------------------------------*/
/*  sel_expr(sel,row,col,rsel,csel,idx,dval)                                */
/*                                                                          */
/*  check selection expression: <rowsel,colsel,x>                           */
/*  return selected row and column indices in rsel and csel.                */
/*  if x refers to a matrix return its index in idx, otherwise x must be    */
/*  a scalar constant and its value will be returned in dval.               */
/*                                                                          */
/*  return 0 if OK, -1 if error.                                            */

int sel_expr(TDAContext *ctx, char *sel,int row,int col,int *rsel,int *csel,int *idx,double *dval)
{
    int i,ii,cflag,mflag,n,a,b,d,*nsel;
    register char *p;

    *dval = 0.0;
    *idx = -1;

    p = sel;
    for (ii = 0; ii < 2; ++ii) {
        cflag = mflag = 0;
        if (ii == 0) {
            n = row;
            nsel = rsel;
        }
        else {
            n = col;
            nsel = csel;
        }
        for (i = 0; i < n; ++i)
            nsel[i] = 0;

        if (*p == '*') {
            for (i = 0; i < n; ++i)
                nsel[i] = 1;
            p++;
        }
        else {
            if (*p == '+')
                p++;
            else if (*p == '-') {
                mflag = 1;
                p++;
            }
            if (*p == '(') {
                cflag = 1;
                p++;
            }
            if (sscanf(p,"%d(%d)%d",&a,&d,&b) == 3) {
                while (a <= b) {
                    if (a >= 1 && a <= n)
                        nsel[a - 1] = 1;
                    a += d;
                }
                p = skip_int(ctx, p) + 1;
                p = skip_int(ctx, p) + 1;
                p = skip_int(ctx, p);
            }
            else {
                while (*p) {
                    if (sscanf(p,"%d",&a) != 1)  
                        goto SELEFin;

                    if (a >= 1 && a <= n)  
                        nsel[a - 1] = 1;
                    p = skip_int(ctx, p);
                    if (cflag == 0 || *p == ')')
                        break;
                    if (*p++ != ',')
                        goto SELEFin;

                }
            }
            if (cflag != 0) {
                if (*p != ')')
                    goto SELEFin;
                p++;
            }
            if (mflag) {
                for (i = 0; i < n; ++i) {
                    if (nsel[i] == 0)
                        nsel[i] = 1;
                    else
                        nsel[i] = 0;
                }   
            }
        }
        if (*p++ != ',')
            goto SELEFin;
    }
    if (sscanf(p,"%lf",dval) != 1) {        /* check for constant */
        if ((*idx = mat_getidx(ctx, p,1)) < 0)    /* check for matrix */
            goto SELEFin;
    }
    return(0);

SELEFin:
    printf1(ctx, "Error in selection string.\n");
    return(-1);
}

/*--------------------------------------------------------------------------*/
/*  mx_free()       Free MX[].                                              */

void mx_free(TDAContext *ctx)
{
    register int i;

    for (i = 0; i < MXN; ++i) {
        if (ctx->MXAlloc[i] > 0) {
            free((char *)ctx->MX[i]);
            memrq(ctx, -ctx->MXAlloc[i],sizeof(double));
            ctx->MXAlloc[i] = 0;
        }
        if (ctx->MX1Alloc[i] > 0) {
            free((char *)ctx->MX1[i]);
            memrq(ctx, -ctx->MX1Alloc[i],sizeof(double));
            ctx->MX1Alloc[i] = 0;
        }
    }
}

/*--------------------------------------------------------------------------*/
/*  mx_alloc(idx,row,col)   Allocate MX[idx]                                */
/*                          Return 0 if OK, -1 if error.                    */

int mx_alloc(TDAContext *ctx, int idx,int row,int col)
{

    if (idx < 0 || idx >= MXN || row < 1 || col < 1) {
        tda_err("ERROR in mx_alloc: idx=%d row=%d col=%d\n",idx,row,col);
        exit(0);
    }
    if (ctx->MXAlloc[idx] > 0) {
        free((char *)ctx->MX[idx]);
        memrq(ctx, -ctx->MXAlloc[idx],sizeof(double));
        ctx->MXAlloc[idx] = 0;
    }
    if (!(ctx->MX[idx] = (double *)calloc((size_t)(row) * (size_t)(col) + 1,sizeof(double)))) {
        m_cmdmsg(ctx);
        p_err(ctx, -2,1);
        return(-1);
    }         
    memrq(ctx, row * col + 1,sizeof(double));
    ctx->MXAlloc[idx] = row * col + 1;
    ctx->MXRow[idx] = row;
    ctx->MXCol[idx] = col;
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  mx1_alloc(idx,row,col)   Allocate MX1[idx]                              */
/*                          Return 0 if OK, -1 if error.                    */

int mx1_alloc(TDAContext *ctx, int idx,int row,int col)
{

    if (idx < 0 || idx >= MXN || row < 1 || col < 1) {
        tda_err("ERROR in mx1_alloc: idx=%d row=%d col=%d\n",idx,row,col);
        exit(0);
    }
    if (ctx->MX1Alloc[idx] > 0) {
        free((char *)ctx->MX1[idx]);
        memrq(ctx, -ctx->MX1Alloc[idx],sizeof(double));
        ctx->MX1Alloc[idx] = 0;
    }
    if (!(ctx->MX1[idx] = (double *)calloc((size_t)(row) * (size_t)(col) + 1,sizeof(double)))) {
        m_cmdmsg(ctx);
        p_err(ctx, -2,1);
        return(-1);
    }         
    memrq(ctx, row * col + 1,sizeof(double));
    ctx->MX1Alloc[idx] = row * col + 1;
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  get_mxval(idx,row,col)  Return MX[idx][row,col]                         */
/*                          row, col runs from 1,2,3,...                    */

double get_mxval(TDAContext *ctx, int idx,int row,int col)
{
    register int i,j;

    i = (--row) % ctx->MXRow[idx];
    j = (--col) % ctx->MXCol[idx] + 1;
    return(ctx->MX[idx][i * ctx->MXCol[idx] + j]);
}

/*--------------------------------------------------------------------------*/
/*  get_matval(idx,row,col) Return MatVal[idx][row,col]                     */
/*                          row, col runs from 1,2,3,...                    */

double get_matval(TDAContext *ctx, int idx,int row,int col)
{
    register int i,j;

    i = (--row) % ctx->MatRow[idx];
    j = (--col) % ctx->MatCol[idx] + 1;
    return(ctx->MatVal[idx][i * ctx->MatCol[idx] + j]);
}

/*--------------------------------------------------------------------------*/
/*  m_copy(x,y,n)       copy n elements from y to x.                        */

void m_copy(TDAContext *ctx, double *x,double *y,int n)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;
    for (i = 1; i <= n; ++i)
        x[i] = y[i];
}

/*--------------------------------------------------------------------------*/
/*  m_check1(p)     check whether *p is a comma. If it is return p++,       */
/*                  otherwise return NULL                                   */

char *m_check1(TDAContext *ctx, char *p)                   
{
    if (*p++ != ',' || !*p) {
        mat_err(ctx, 0);
        return(NULL);
    }
    return(p);
}

/*--------------------------------------------------------------------------*/
/*  m_check2(p,opt) check whether *p is a ")". If it is return p++,         */
/*                  otherwise return NULL                                   */
/*                  If opt != 0, also check that "0" is last character      */
/*                  in string.                                              */

char *m_check2(TDAContext *ctx, char *p,int opt)           
{
    if (*p++ != ')' || (opt && *p)) {
        mat_err(ctx, 0);
        return(NULL);
    }
    return(p);
}

/*--------------------------------------------------------------------------*/
/*  m_getmat(p,row,col,idx,cmd,opt)                                         */
/*                                                                          */
/*                              Get idx for new matrix pointed to by p.     */
/*                              Return pointer to next character after      */
/*                              name, or NULL if error.                     */
/*                              Copy cmd into matrix definition.            */
/*                                                                          */  
/*  If opt != 0 check that the string ends with ")".                        */  

char *m_getmat(TDAContext *ctx, char *p,int row,int col,int *idx,char *cmd,int opt)
{
    char mname[VNLMax + 1];  

    if ((p = get_mname(ctx, p,mname,1)) == NULL)
        return(p);
    if ((*idx = mat_newmat(ctx, mname,row,col)) < 0)
        return(NULL);
    mdefcpy(ctx, ctx->MatDef[*idx],cmd);
    if (opt)
        p = m_check2(ctx, p,1);
    return(p);
}

/*--------------------------------------------------------------------------*/
/*  m_cent(cmd,opt)  opt=0:  mcent(A,R)                                     */
/*                   opt=1:  mstand(A,R)                                    */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_cent(TDAContext *ctx, char *cmd,int opt)
{
    register int i,j;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp,tmp1;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = cmd + 6;
    if (opt)
        p++;

    p = n_mexpr(ctx, p,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MCENTFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MCENTFin;
    if ((p = m_getmat(ctx, p,row,col,&idx,cmd,1)) == NULL)
        goto MCENTFin;

    for (j = 1; j <= col; ++j) {
        tmp = 0.0;
        for (i = 1; i <= row; ++i)  
            tmp += get_mxval(ctx, 0,i,j);

        tmp /= (double)row;
        for (i = 0; i < row; ++i)  
            ctx->MatVal[idx][i * col + j] = get_mxval(ctx, 0,i + 1,j) - tmp;

        if (opt) {                      /* mstand */
            tmp = 0.0;
            for (i = 0; i < row; ++i) {
                tmp1 = ctx->MatVal[idx][i * col + j];
                tmp += tmp1 * tmp1;
            }
            tmp /= (double)row;

            if (tmp > 0.0) {
                tmp = sqrt(tmp);
                for (i = 0; i < row; ++i)    
                    ctx->MatVal[idx][i * col + j] /= tmp;
            }
        }
    }
    err = 0;

MCENTFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_dcent(cmd)       mdcent(R,A)                                          */
/*                     R should be a (symmetric) n x n matrix.              */
/*                     The command creates a double-centered matrix A       */
/*                     as used with classical metric scaling.               */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_dcent(TDAContext *ctx, char *cmd)
{
    register int i,j,k;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp,aij,r00,ri0,r0j;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = cmd + 7;

    p = n_mexpr(ctx, p,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MDCENTFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MDCENTFin;
    if (row != col) {
        mat_err(ctx, 1);
        goto MDCENTFin;
    }
    if ((p = m_getmat(ctx, p,row,col,&idx,cmd,1)) == NULL)
        goto MDCENTFin;

    r00 = 0.0;
    for (i = 1; i <= row; ++i) {
        for (j = 1; j <= col; ++j) {
            tmp = get_mxval(ctx, 0,i,j);
            r00 += tmp * tmp;
        }
    }
    r00 /= (double)(row * row);

    for (i = 1; i <= row; ++i) {
        ri0 = 0.0;
        for (j = 1; j <= col; ++j) {
            tmp = get_mxval(ctx, 0,i,j);
            ri0 += tmp * tmp;
        }
        ri0 /= (double)col;

        for (j = 1; j <= col; ++j) {
            tmp = get_mxval(ctx, 0,i,j);
            aij = r00 + tmp * tmp - ri0;

            r0j = 0.0;
            for (k = 1; k <= row; ++k) {
                tmp = get_mxval(ctx, 0,k,j);
                r0j += tmp * tmp;
            }
            r0j /= (double)row;
            aij -= r0j;
            ctx->MatVal[idx][(i - 1) * col + j] = -aij / 2.0;
        }
    }
    err = 0;

MDCENTFin:
    mx_free(ctx);
    return(err);
}
     
/*--------------------------------------------------------------------------*/
/*  m_cross(cmd)    mcross(A,R)          R = A'A                            */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_cross(TDAContext *ctx, char *cmd)
{
    register int i,j,l;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 7,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MCROSSFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MCROSSFin;
    if ((p = m_getmat(ctx, p,col,col,&idx,cmd,1)) == NULL)
        goto MCROSSFin;

    for (i = 0; i < col; ++i) {
        for (j = 1; j <= col; ++j) {
            tmp = 0.0;
            for (l = 1; l <= row; ++l)  
                tmp += get_mxval(ctx, 0,l,i + 1) * get_mxval(ctx, 0,l,j);
            ctx->MatVal[idx][i * col + j] = tmp;
        }
    }
    err = 0;

MCROSSFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_chol(cmd)      mchol(A,R)     Cholesky decomposition                  */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_chol(TDAContext *ctx, char *cmd)
{
    register int i,j,k;
    int err,row,col,idx,n,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 6,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MCHOLFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MCHOLFin;
    if (row != col) {
        mat_err(ctx, 1);
        goto MCHOLFin;
    }
    if ((p = m_getmat(ctx, p,row,col,&idx,cmd,1)) == NULL)
        goto MCHOLFin;

    n = (row * (row + 1)) / 2;
    if (alloc_actmp(ctx, n + 1))
        goto MCHOLFin;

    k = 1;
    for (i = 1; i <= row; ++i) {
        for (j = 1; j <= i; ++j)
            ctx->AcTmp[k++] = get_mxval(ctx, 0,i,j);  
    }
    if (decomp(ctx, row,ctx->AcTmp)) {
        m_cmdmsg(ctx);
        printf1(ctx, "Error: %s is not positive definite.\n",ctx->MXName);
        goto MCHOLFin;     
    }
    k = 1;
    for (i = 0; i < row; ++i) {
        for (j = 0; j <= i; ++j)
            ctx->MatVal[idx][i * row + j + 1] = ctx->AcTmp[k++];
    }
    err = 0;

MCHOLFin:
    alloc_actmp(ctx, 0);
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_agg(cmd)      Aggregation of a matrix.                                */
/*                                                                          */
/*  mag(A,C,R,B)                                                            */
/*                                                                          */
/*  Given: a (n,m)-matrix A, a (n,1) vector C and a (1,m)-vector R.         */
/*  The command creates a new (u,v)-matrix B that is an aggregate of A.     */
/*                                                                          */
/*  The vector C must contain indices from 1 to u, the vector R must        */
/*  contain indices from 1 to v. Then all rows of A having the same index   */
/*  in C and all columns of A having the same index in R are aggrated into  */
/*  a single element of B.                                                  */  
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_agg(TDAContext *ctx, char *cmd)
{
    register int i,j,ii,jj,k;
    int err,row,row1,row2,col,col1,col2,idx,ivflg,nr,nc,rmax,cmax,ir,ic;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    if ((p = n_mexpr(ctx, cmd + 4,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MAGFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MAGFin;
    if ((p = n_mexpr(ctx, p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
        goto MAGFin;
    if (row1 != row || col1 != 1) {
        mat_err(ctx, 2);
        goto MAGFin;
    }
    if ((p = m_check1(ctx, p)) == NULL)
        goto MAGFin;
    if ((p = n_mexpr(ctx, p,1,2,&row2,&col2,&ivflg,NULL)) == NULL)
        goto MAGFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MAGFin;
    if (row2 != 1 || col2 != col) {
        mat_err(ctx, 2);
        goto MAGFin;
    }
    rmax = 0;
    for (i = 1; i <= row; ++i)  
        rmax = imax(ctx, rmax,(int)ctx->MX[1][i]);
    cmax = 0;
    for (j = 1; j <= col; ++j)  
        cmax = imax(ctx, cmax,(int)ctx->MX[2][j]);

    if (rmax < 1 || cmax < 1) {
        mat_err(ctx, 15);
        goto MAGFin;
    }
    if (alloc_acr(ctx, rmax + 1))
        goto MAGFin;
    if (alloc_acs(ctx, cmax + 1))
        goto MAGFin;
      
    for (i = 1; i <= row; ++i) {
        if ((k = (int)ctx->MX[1][i]) > 0)  
            ctx->AcR[k] = 1;
    }
    for (j = 1; j <= col; ++j) {
        if ((k = (int)ctx->MX[2][j]) > 0)  
            ctx->AcS[k] = 1;
    }
    nr = nc = 0;
    for (i = 1; i <= rmax; ++i) {
        if (ctx->AcR[i] > 0)
            nr++;
    }
    for (j = 1; j <= cmax; ++j) {
        if (ctx->AcS[j] > 0)
            nc++;
    }

    if ((p = m_getmat(ctx, p,nr,nc,&idx,cmd,1)) == NULL)
        goto MAGFin;
      
    ir = 0;
    for (ii = 1; ii <= rmax; ++ii) {
        if (ctx->AcR[ii] == 0)
            continue;
        ir++;
        k = (ir - 1) * nc;
        ic = 0;

        for (jj = 1; jj <= cmax; ++jj) {
            if (ctx->AcS[jj] == 0)
                continue;
            ic++;

            tmp = 0.0;
            for (i = 1; i <= row; ++i) {
                if (ctx->MX[1][i] == ii) {
                    for (j = 1; j <= col; ++j) {
                        if (ctx->MX[2][j] == jj)  
                            tmp += get_mxval(ctx, 0,i,j);
                    }
                }
            }
            ctx->MatVal[idx][k + ic] = tmp;
        }
    }
    err = 0;

MAGFin:
    alloc_acr(ctx, 0);
    alloc_acs(ctx, 0);
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_vec(cmd,opt)      opt 0 :  mcvec(A,V)     stack columns               */
/*                      opt 1 :  mrvec(A,V)     stack rows                  */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_vec(TDAContext *ctx, char *cmd,int opt)
{
    register int i,j,k;
    int err,row,col,idx,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 6,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MVECFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MVECFin;

    if ((p = m_getmat(ctx, p,row * col,1,&idx,cmd,1)) == NULL)
        goto MVECFin;

    k = 1;
    if (opt == 0) {
        for (j = 1; j <= col; ++j) {
            for (i = 1; i <= row; ++i)     
                ctx->MatVal[idx][k++] = get_mxval(ctx, 0,i,j);
        }
    }
    else {
        for (i = 1; i <= row; ++i) {
            for (j = 1; j <= col; ++j)     
                ctx->MatVal[idx][k++] = get_mxval(ctx, 0,i,j);
        }
    }
    err = 0;

MVECFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_ivec(cmd)     mivec(V,n,A)   inverse vec operator                     */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_ivec(TDAContext *ctx, char *cmd)
{
    register int i,j,k;
    int err,row,col,n,m,q,idx,ivflg;
    register char *p;

    ivflg = 0;
    err = -1;
    printf2(ctx, "%s\n",cmd);

    if ((p = n_mexpr(ctx, cmd + 6,1,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto MIVECFin;

    if (col != 1) {
        mat_err(ctx, 7);
        goto MIVECFin;
    }
    q = row;

    if ((p = m_check1(ctx, p)) == NULL)
        goto MIVECFin;
    if ((p = n_mexpr(ctx, p,1,1,&row,&col,&ivflg,NULL)) == NULL) 
        goto MIVECFin;

    if (row != 1 || col != 1 || (n = (int)ctx->MX[1][1]) < 1) {
        mat_err(ctx, 18);
        goto MIVECFin;
    }
    m = q / n;
    if (n * m != q) {
        mat_err(ctx, 18);
        goto MIVECFin;
    }
    if ((p = m_check1(ctx, p)) == NULL)
        goto MIVECFin;
    if ((p = m_getmat(ctx, p,n,m,&idx,cmd,0)) == NULL)
        goto MIVECFin;

    k = 1;
    for (j = 1; j <= m; ++j) {
        for (i = 0; i < n; ++i)  
            ctx->MatVal[idx][i * m + j] = ctx->MX[0][k++];
    }   
    err = 0;

MIVECFin:
    mx_free(ctx);
    return(err);
}

/*--##----------------------------------------------------------------------*/
/*  m_mldes(cmd)    mldes(Z,G,D)   Design matrix for MLRC models            */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mldes(TDAContext *ctx, char *cmd)
{
    register int i1,i2,i,j,k,l,ii,kk,ll;
    int err,row,col,n,m,q,nn,idx,ivflg;
    register char *p;
    double *d,*zi,*zii;

    ivflg = 0;
    err = -1;
    printf2(ctx, "%s\n",cmd);

    if ((p = n_mexpr(ctx, cmd + 6,1,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto MLDESFin;

    n = row;
    q = col;
    m = q * q + 1;
    nn = n * n;

    if ((p = m_check1(ctx, p)) == NULL)
        goto MLDESFin;
    if ((p = n_mexpr(ctx, p,1,1,&row,&col,&ivflg,NULL)) == NULL) 
        goto MLDESFin;

    if (row != n || col != 1) {
        mat_err(ctx, 2);
        goto MLDESFin;
    }
    if ((p = m_check1(ctx, p)) == NULL)
        goto MLDESFin;
    if ((p = m_getmat(ctx, p,nn,m,&idx,cmd,0)) == NULL)
        goto MLDESFin;
       
    if (alloc_aci(ctx, n + 1))
        goto MLDESFin;
    if (alloc_acj(ctx, n + 1))
        goto MLDESFin;

    ctx->AcI[1] = ctx->AcJ[1] = 1;

    i = j = 1;
    for (k = 2; k <= n; ++k) {
        i++;
        if ((int)ctx->MX[1][k] != (int)ctx->MX[1][k - 1]) {
            i = 1;
            j++;
        }
        ctx->AcI[k] = i;
        ctx->AcJ[k] = j; 
    }
    l = 0;
    for (kk = 1; kk <= n; ++kk) {
        ii = ctx->AcI[kk];
        j = ctx->AcJ[kk];
        for (k = 1; k <= n; ++k) {
            i = ctx->AcI[k];

            /** tda_out("\ni=%2d j=%2d ii=%2d jj=%2d \n",i,AcJ[k],ii,j); **/

            if (ctx->AcJ[k] == j) {
                zi  = ctx->MX[0] + (k - 1) * q;
                zii = ctx->MX[0] + (kk - 1) * q;

                /** tda_out("zi=%lf %lf  zii=%lf %lf\n",zi[1],zi[2],zii[1],zii[2]); **/

                ll = 1;
                d = ctx->MatVal[idx] + l * m;
                for (i2 = 1; i2 <= q; ++i2) {
                    for (i1 = 1; i1 <= q; ++i1)  
                        d[ll++] = zi[i1] * zii[i2]; 
                }
                if (i == ii)
                    d[ll] = 1.0;                           
            }
            l++;
        }
    }
    err = 0;
      
MLDESFin:
    mx_free(ctx);
    alloc_aci(ctx, 0);
    alloc_acj(ctx, 0);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_sum(cmd,opt)      opt 0 :  mrsum(A,V)     row sums                    */
/*                      opt 1 :  mcsum(A,U)     column sums                 */
/*                                                                          */
/*  V(i) = A(i,1) + ... + A(i,m)                                            */
/*  U(j) = A(1,j) + ... + A(n,j)                                            */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_sum(TDAContext *ctx, char *cmd,int opt)
{
    register int i,j;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 6,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MSUMFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MSUMFin;

    if (opt == 0) {
        if ((p = m_getmat(ctx, p,row,1,&idx,cmd,1)) == NULL)
            goto MSUMFin;

        for (i = 1; i <= row; ++i) {
            tmp = 0.0;
            for (j = 1; j <= col; ++j) 
                tmp += get_mxval(ctx, 0,i,j);
            ctx->MatVal[idx][i] = tmp;
        }
    }
    else {
        if ((p = m_getmat(ctx, p,1,col,&idx,cmd,1)) == NULL)
            goto MSUMFin;

        for (j = 1; j <= col; ++j) {
            tmp = 0.0;
            for (i = 1; i <= row; ++i) 
                tmp += get_mxval(ctx, 0,i,j);
            ctx->MatVal[idx][j] = tmp;
        }
    }
    err = 0;

MSUMFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_drow(cmd,opt)     opt 0 :  mdrow(A,R)                                 */
/*                      opt 1 :  mdcol(A,R)                                 */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_drow(TDAContext *ctx, char *cmd,int opt)
{
    register int i,j;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 6,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MDROWFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MDROWFin;

    if (opt == 0) {
        if ((p = m_getmat(ctx, p,row,row,&idx,cmd,1)) == NULL)
            goto MDROWFin;

        for (i = 1; i <= row; ++i) {
            tmp = 0.0;
            for (j = 1; j <= col; ++j) 
                tmp += get_mxval(ctx, 0,i,j);
            ctx->MatVal[idx][(i - 1) * row + i] = tmp;
        }
    }
    else {
        if ((p = m_getmat(ctx, p,col,col,&idx,cmd,1)) == NULL)
            goto MDROWFin;

        for (j = 1; j <= col; ++j) {
            tmp = 0.0;
            for (i = 1; i <= row; ++i) 
                tmp += get_mxval(ctx, 0,i,j);
            ctx->MatVal[idx][(j - 1) * col + j] = tmp;
        }
    }
    err = 0;

MDROWFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_diag(cmd,opt)     opt 0 :  mdiag(A,R)                                 */
/*                      opt 1 :  mdiagd(A,R)                                */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_diag(TDAContext *ctx, char *cmd,int opt)
{
    register int i,j;
    int err,row,col,idx,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = cmd + 6;
    if (opt)
        p++;
    p = n_mexpr(ctx, p,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MDIAGFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MDIAGFin;

    if (opt == 0) {
        if (col != 1) {
            mat_err(ctx, 7);
            goto MDIAGFin;
        }
        if ((p = m_getmat(ctx, p,row,row,&idx,cmd,1)) == NULL)
            goto MDIAGFin;

        for (i = 1; i <= row; ++i)  
            ctx->MatVal[idx][(i - 1) * row + i] = get_mxval(ctx, 0,i,1);
    }
    else {
        j = imin(ctx, row,col);
        if ((p = m_getmat(ctx, p,j,1,&idx,cmd,1)) == NULL)
            goto MDIAGFin;

        for (i = 1; i <= j; ++i)  
            ctx->MatVal[idx][i] = get_mxval(ctx, 0,i,i);
    }
    err = 0;

MDIAGFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_transp(cmd)    mtransp(A,R)    Transposition                          */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_transp(TDAContext *ctx, char *cmd)
{
    register int i,j;
    int err,row,col,idx,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 8,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MTRANFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MTRANFin;

    if ((p = m_getmat(ctx, p,col,row,&idx,cmd,1)) == NULL)
        goto MTRANFin;

    for (i = 1; i <= row; ++i) {
        for (j = 1; j <= col; ++j)  
            ctx->MatVal[idx][(j - 1) * row + i] = get_mxval(ctx, 0,i,j);
    }
    err = 0;

MTRANFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_sqrt(cmd,opt)     opt 0 :  msqrtd(A,R)                                */
/*                      opt 1 :  msqrti(A,R)                                */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_sqrt(TDAContext *ctx, char *cmd,int opt)
{
    register int i,j;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 7,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MSQRTFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MSQRTFin;

    if (row != col) {
        mat_err(ctx, 1);
        goto MSQRTFin;
    }
    if ((p = m_getmat(ctx, p,row,col,&idx,cmd,1)) == NULL)
        goto MSQRTFin;

    j = 1;
    for (i = 1; i <= row; ++i) {
        tmp = get_mxval(ctx, 0,i,i); 
        if (tmp < 0.0) {
            m_cmdmsg(ctx);
            printf1(ctx, "Error: negative diagonal element in %s.\n",ctx->MXName);
            goto MSQRTFin;
        }
        if (tmp > 0.0)
            tmp = sqrt(tmp);
        if (opt == 0) 
            ctx->MatVal[idx][j] = tmp;
        else {
            if (tmp <= ctx->EPSI) {
                m_cmdmsg(ctx);
                printf1(ctx, "Error: (almost) zero diagonal element in %s.\n",ctx->MXName);
                goto MSQRTFin;
            }
            ctx->MatVal[idx][j] = 1.0 / tmp;
        }
        j += row + 1;
    }
    err = 0;

MSQRTFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_nrow(cmd,opt)     opt 0 :  mnrow(A,R)                                 */
/*                      opt 1 :  mncol(A,R)                                 */
/*                      opt 2 :  mnorm(A,R)                                 */
/*                      opt 3 :  mnorm1(A,R)                                */
/*                      opt 4 :  mnorm2(A,R)                                */
/*                      opt 5 :  mtrace(A,R)                                */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_nrow(TDAContext *ctx, char *cmd,int opt)
{
    register int i,n;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = cmd + 6;
    if (opt >= 3)
        p++;
    p = n_mexpr(ctx, p,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MNROWFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MNROWFin;

    if ((p = m_getmat(ctx, p,1,1,&idx,cmd,1)) == NULL)
        goto MNROWFin;

    switch (opt) {
        case  1:    tmp = (double)col;
                    break;
            
        case  2:    n = row * col;
                    tmp = 0.0;
                    for (i = 1; i <= n; ++i)
                        tmp = dmax(ctx, tmp,fabs(ctx->MX[0][i]));
                    break;

        case  3:    n = row * col;
                    tmp = 0.0;
                    for (i = 1; i <= n; ++i)
                        tmp += fabs(ctx->MX[0][i]);
                    break;

        case  4:    n = row * col;
                    tmp = 0.0;
                    for (i = 1; i <= n; ++i)
                        tmp += ctx->MX[0][i] * ctx->MX[0][i];
                    if (tmp > 0.0)
                        tmp = sqrt(tmp);
                    break;

        case  5:    n = imin(ctx, row,col);
                    tmp = 0.0;
                    for (i = 0; i < n; ++i)
                        tmp += ctx->MX[0][i * col + i + 1];
                    break;

        default:    tmp = (double)row;
                    break;
    }
    ctx->MatVal[idx][1] = tmp;
    err = 0;

MNROWFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_num(cmd)      mnum(x,d,n,R)                                           */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_num(TDAContext *ctx, char *cmd)
{
    register int i;
    int err,row,col,idx,n,ivflg;
    register char *p;
    double x,d;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 5,3,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MNUMFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MNUMFin;
    if (row != 1 || col != 1) {
        mat_err(ctx, 8);
        goto MNUMFin;
    }
    n = (int)get_mxval(ctx, 2,1,1);
    if (n < 1) {
        mat_err(ctx, 2);
        goto MNUMFin;
    }
    if ((p = m_getmat(ctx, p,n,1,&idx,cmd,1)) == NULL)
        goto MNUMFin;

    x = get_mxval(ctx, 0,1,1);
    d = get_mxval(ctx, 1,1,1);

    for (i = 1; i <= n; ++i)
        ctx->MatVal[idx][i] = x + d * (double)(i - 1);
    err = 0;

MNUMFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mch(cmd)  Check row and column sum conditions.                        */
/*                                                                          */
/*  mch(A,B)                                                                */
/*                                                                          */
/*  A (n,m), U (1,m), V (n,1)                                               */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mch(TDAContext *ctx, char *cmd)
{
    register int i,j;
    int err,row,col,idx,ivflg,n,nf;         
    register char *p;
    double z,s;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = cmd + 4;
    if ((p = n_mexpr(ctx, p,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MCHFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MCHFin;

    if (ctx->MXCol[0] != ctx->MXRow[0]) {       
        mat_err(ctx, 1);
        goto MCHFin;
    }
    n = ctx->MXRow[0];

    if ((p = m_getmat(ctx, p,n,1,&idx,cmd,1)) == NULL)
        goto MCHFin;

    n = ctx->MXRow[0];
    for (i = 0; i < n; ++i) {
        z = s = 0.0;
        nf = 0;
        for (j = i + 1; j < n; ++j) {
            z += ctx->MX[0][i * n + j + 1];
            s += ctx->MX[0][j * n + i + 1];
            /** tda_out("i=%d z=%g s=%g\n",i,z,s); **/
            if (z < s - ctx->EPSI1) {
                /* printf2("mch: i=%d z=%g s=%g\n",i + 1,z,s); */
                nf = 1;
                break;
            }
        }
        if (nf)
            ctx->MatVal[idx][i + 1] = 1.0;
    }
    for (i = n - 1; i > 0; --i) {
        z = s = 0.0;
        nf = 0;
        for (j = i - 1; j >= 0; --j) {
            z += ctx->MX[0][j * n + i + 1];
            s += ctx->MX[0][i * n + j + 1];
            /** tda_out("i=%d z=%g s=%g\n",i,z,s); **/
            if (z < s - ctx->EPSI1) {
                /* printf2("mch: i=%d z=%g s=%g\n",i + 1,z,s); */
                nf = 1;
                break;
            }
        }
        if (nf)
            ctx->MatVal[idx][i + 1] = 1.0;
    }
    err = 0;

MCHFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mpfit(cmd)    Iterative proportional fitting.                         */
/*                                                                          */
/*  mpfit(A,U,V,ITER,EPS,B)                                                 */
/*                                                                          */
/*  A (n,m), U (1,m), V (n,1)                                               */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mpfit(TDAContext *ctx, char *cmd)
{
    int err,row,col,idx,ivflg,iter;
    register char *p;
    double eps;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = cmd + 6;
    if ((p = n_mexpr(ctx, p,5,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MRASFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MRASFin;

    if (ctx->MXCol[0] != ctx->MXCol[1] || ctx->MXRow[1] != 1 ||
        ctx->MXRow[0] != ctx->MXRow[2] || ctx->MXCol[2] != 1 ||    
        ctx->MXRow[3] != 1        || ctx->MXCol[3] != 1 ||    
        ctx->MXRow[4] != 1        || ctx->MXCol[4] != 1) {    
        mat_err(ctx, 2);
        goto MRASFin;
    }
    if ((p = m_getmat(ctx, p,row,col,&idx,cmd,1)) == NULL)
        goto MRASFin;

    iter = (int)ctx->MX[3][1];
    eps  = ctx->MX[4][1];
    if (iter < 1) { 
        printf1(ctx, "mpfit command: error in number of iterations.\n");
        goto MRASFin;
    }
    if (eps < 0.0) { 
        printf1(ctx, "mpfit command: error in epsilon.\n");
        goto MRASFin;
    }
    if (mpfit(ctx, row,col,ctx->MX[0],ctx->MX[1],ctx->MX[2],ctx->MatVal[idx],&iter,&eps))
        goto MRASFin;
    printf2(ctx, "Number of iterations: %d. Final accuracy: %g\n",iter,eps);
    err = 0;

MRASFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mpinv(cmd)   Inverse permutation.                                     */
/*                                                                          */
/*  mpinv(P,Q)      Q(P(i)) = i.                                            */
/*                                                                          */
/*  P must be (n,1) or (1,n) vector                                         */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mpinv(TDAContext *ctx, char *cmd)
{
    register int i,k;
    int err,row,col,idx,ivflg,n;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = cmd + 6;
    if ((p = n_mexpr(ctx, p,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MPINVFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MPINVFin;

    if (row != 1 && col != 1) {
        mat_err(ctx, 2);
        goto MPINVFin;
    }
    if ((p = m_getmat(ctx, p,row,col,&idx,cmd,1)) == NULL)
        goto MPINVFin;

    n = imax(ctx, row,col);
    for (i = 1; i <= n; ++i) {
        k = (int)ctx->MX[0][i];
        if (k < 1 || k > n) {
            mat_err(ctx, 17);
            goto MPINVFin;
        }
        ctx->MatVal[idx][k] = (double)i;
    }
    err = 0;

MPINVFin:
    mx_free(ctx);
    return(err);
}


/*--------------------------------------------------------------------------*/
/*  m_mperm(cmd,opt)    Permutation of a matrix.                            */
/*                                                                          */
/*  opt 0:  mpsym(A,P,B)    b(i,j) = a(p(i),p(j))                           */
/*  opt 1:  mprow(A,P,B)    b(i,j) = a(p(i),j)                              */
/*  opt 2:  mpcol(A,P,B)    b(i,j) = a(i,p(j))                              */
/*                                                                          */
/*  opt 0:  A must be (n,n), p = (n,1) or (1,n)                             */
/*  opt 1:  A (n,m), p (n,1) or (1,n)                                       */ 
/*  opt 2:  A (n,m), p (1,m) or (m,1)                                       */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mperm(TDAContext *ctx, char *cmd,int opt)
{
    register int i,j,k,k1;
    int err,row,col,row1,col1,idx,ivflg,n;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = cmd + 6;
    if ((p = n_mexpr(ctx, p,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MPERMFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MPERMFin;

    if ((p = n_mexpr(ctx, p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
        goto MPERMFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MPERMFin;

    if (row1 != 1 && col1 != 1) {
        mat_err(ctx, 2);
        goto MPERMFin;
    }
    if ((p = m_getmat(ctx, p,row,col,&idx,cmd,1)) == NULL)
        goto MPERMFin;

    n = imax(ctx, row1,col1);
    if (opt == 0) {
        for (i = 1; i <= n; ++i) {
            k = (int)ctx->MX[1][i];
            if (k < 1 || k > row) {
                mat_err(ctx, 17);
                goto MPERMFin;
            }
        }
    }
    if (opt == 1) {
        if (n != row) {  
            mat_err(ctx, 2);
            goto MPERMFin;
        }
        for (i = 0; i < row; ++i) {
            k = (int)ctx->MX[1][i + 1] - 1;
            if (k < 0 || k >= row) {
                mat_err(ctx, 17);
                goto MPERMFin;
            }
            for (j = 1; j <= col; ++j)  
                ctx->MatVal[idx][i * col + j] = ctx->MX[0][k * col + j];
        }
    }
    else if (opt == 2) {
        if (n != col) {  
            mat_err(ctx, 2);
            goto MPERMFin;
        }
        for (j = 1; j <= col; ++j) {
            k = (int)ctx->MX[1][j];
            if (k < 1 || k > col) {
                mat_err(ctx, 17);
                goto MPERMFin;
            }
            for (i = 0; i < row; ++i)  
                ctx->MatVal[idx][i * col + j] = ctx->MX[0][i * col + k];
        }
    }
    else {
        if (row != col || n != row) {  
            mat_err(ctx, 2);
            goto MPERMFin;
        }
        for (i = 0; i < row; ++i) {
            k = (int)ctx->MX[1][i + 1] - 1;
            for (j = 1; j <= col; ++j) {
                k1 = (int)ctx->MX[1][j];
                ctx->MatVal[idx][i * col + j] = ctx->MX[0][k * col + k1];
            }
        }
    }
    err = 0;

MPERMFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mqap(cmd)     Quadratic assignment                                    */
/*                                                                          */
/*  mqap(F,D,C,P)   For a description see qap_w() in t_ass.c                */
/*                                                                          */
/*  Note: The main diagonal of F and D will be set to zero.                 */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mqap(TDAContext *ctx, char *cmd)
{
    register int i,j,n;
    int err,row,col,idx,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = cmd + 5;
    if ((p = n_mexpr(ctx, p,3,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MQAPFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MQAPFin;

    if (ctx->MXRow[0] != row || ctx->MXCol[0] != col ||
        ctx->MXRow[1] != row || ctx->MXCol[1] != col || 
        ctx->MXRow[2] != row || ctx->MXCol[2] != col) { 
        mat_err(ctx, 2);
        goto MQAPFin;
    }
    if ((p = m_getmat(ctx, p,row,1,&idx,cmd,1)) == NULL)
        goto MQAPFin;

    if (alloc_acn(ctx, 3 * row + 1))
        goto MQAPFin;

    if (alloc_actmp(ctx, 4 * row + 1))
        goto MQAPFin;

    n = 0;
    j = 1;
    for (i = 1; i <= row; ++i) {
        if (ctx->MX[0][j] != 0.0) {
            ctx->MX[0][j] = 0.0;
            n++;
        }
        if (ctx->MX[1][j] != 0.0) {
            ctx->MX[1][j] = 0.0;
            n++;
        }
        j += row + 1;
    }
    if (n > 0)  
        printf1(ctx, "Warning: changed %d main diagonal elements to zero.\n",n);

    if (qap_w(ctx, row,ctx->MX[2],ctx->MX[0],ctx->MX[1],ctx->AcN,ctx->AcTmp)) 
        goto MQAPFin;

    printf1(ctx, "Best value: ");
    rt_printf1_d(ctx, ctx->PMATFmtS,ctx->AcTmp[1]);
    printf1(ctx, "\n");
#ifdef TDA_R_PACKAGE
    /* the best objective value found, as a 1x1 */
    tda_export_cell(ctx, "mqap.value", ctx->AcTmp[1]);
    tda_export_endrow(ctx, "mqap.value");
#endif
    for (i = 1; i <= row; ++i)
        ctx->MatVal[idx][i] = (double)ctx->AcN[i];

    err = 0;

MQAPFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mcel(cmd)     Create edge list.                                       */
/*                                                                          */
/*  mnc(A,x,L)                                                              */
/*                                                                          */
/*  Given an (n,n)-matrix A, this command creates a (m,3)-matrix L,         */
/*  organized as an edge list. Valid edges are assumed if a(i,j) >= x.      */
/*  x must be scalar.                                                       */
/*                                                                          */  
/*  Return 0 if OK, -1 if error.                                            */

int m_mcel(TDAContext *ctx, char *cmd)
{
    register int i,j,k;
    int err,row,col,idx,ivflg,m,n;       
    register char *p;
    double x,aij;         

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = cmd + 5;
    if ((p = n_mexpr(ctx, p,2,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MCELFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MCELFin;

    if (ctx->MXRow[0] != ctx->MXCol[0] || ctx->MXRow[1] != 1 || ctx->MXCol[1] != 1) { 
        mat_err(ctx, 2);
        goto MCELFin;
    }
    row = ctx->MXRow[0];

    x = ctx->MX[1][1];
    m = 0;
    n = row * row;
    for (i = 1; i <= n; ++i) {
        if (ctx->MX[0][i] >= x)
            m++;
    }
    if (m == 0) {
        printf1(ctx, "mcel: number of edges is zero.\n");
        goto MCELFin;
    }
    if ((p = m_getmat(ctx, p,m,3,&idx,cmd,1)) == NULL)
        goto MCELFin;

    k = 0;
    for (i = 0; i < row; ++i) {
        for (j = 1; j <= row; ++j) {
            aij = ctx->MX[0][i * row + j];
            if (aij >= x) {         
                ctx->MatVal[idx][k * 3 + 1] = (double)(i + 1);
                ctx->MatVal[idx][k * 3 + 2] = (double)j;
                ctx->MatVal[idx][k * 3 + 3] = aij;       
                if (++k >= m)
                    break;
            }
        }
    }
    err = 0;

MCELFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mnc(cmd)                                                              */
/*                                                                          */
/*  mnc(A,x,opt,B)                                                          */
/*                                                                          */
/*  opt = 1:   bij = 0 if aij <= x ,  otherwise bij = aij                   */  
/*  opt = 2:   bij = 0 if aij <  x                                          */  
/*  opt = 3:   bij = 0 if aij >= x                                          */  
/*  opt = 4:   bij = 0 if aij >  x                                          */  
/*  opt = 5:   bij = 0 if abs(aij) <= x                                     */  
/*  opt = 6:   bij = 0 if abs(aij) <  x                                     */  
/*  opt = 7:   bij = 0 if abs(aij) >= x                                     */  
/*  opt = 8:   bij = 0 if abs(aij) >  x                                     */  
/*                                                                          */  
/*  Return 0 if OK, -1 if error.                                            */

int m_mnc(TDAContext *ctx, char *cmd)
{
    register int i;
    int err,row,col,idx,ivflg,n,opt;
    register char *p;
    double x;             

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = cmd + 4;
    if ((p = n_mexpr(ctx, p,3,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MNCFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MNCFin;

    if (ctx->MXRow[1] != 1 || ctx->MXCol[1] != 1 || ctx->MXRow[2] != 1 || ctx->MXCol[2] != 1) { 
        mat_err(ctx, 2);
        goto MNCFin;
    }
    row = ctx->MXRow[0];
    col = ctx->MXCol[0];

    if ((p = m_getmat(ctx, p,row,col,&idx,cmd,1)) == NULL)
        goto MNCFin;

    x   = ctx->MX[1][1];
    opt = (int)(ctx->MX[2][1]);
    n = row * col;

    switch (opt) {
        case 1:     for (i = 1; i <= n; ++i) {
                        if (ctx->MX[0][i] <= x)
                            ctx->MatVal[idx][i] = 0.0;
                        else
                            ctx->MatVal[idx][i] = ctx->MX[0][i];
                    }
                    break;

        case 2:     for (i = 1; i <= n; ++i) {
                        if (ctx->MX[0][i] < x)
                            ctx->MatVal[idx][i] = 0.0;
                        else
                            ctx->MatVal[idx][i] = ctx->MX[0][i];
                    }
                    break;

        case 3:     for (i = 1; i <= n; ++i) {
                        if (ctx->MX[0][i] >= x)
                            ctx->MatVal[idx][i] = 0.0;
                        else
                            ctx->MatVal[idx][i] = ctx->MX[0][i];
                    }
                    break;

        case 4:     for (i = 1; i <= n; ++i) {
                        if (ctx->MX[0][i] > x)
                            ctx->MatVal[idx][i] = 0.0;
                        else
                            ctx->MatVal[idx][i] = ctx->MX[0][i];
                    }
                    break;

        case 5:     for (i = 1; i <= n; ++i) {
                        if (fabs(ctx->MX[0][i]) <= x)
                            ctx->MatVal[idx][i] = 0.0;
                        else
                            ctx->MatVal[idx][i] = ctx->MX[0][i];
                    }
                    break;

        case 6:     for (i = 1; i <= n; ++i) {
                        if (fabs(ctx->MX[0][i]) < x)
                            ctx->MatVal[idx][i] = 0.0;
                        else
                            ctx->MatVal[idx][i] = ctx->MX[0][i];
                    }
                    break;

        case 7:     for (i = 1; i <= n; ++i) {
                        if (fabs(ctx->MX[0][i]) >= x)
                            ctx->MatVal[idx][i] = 0.0;
                        else
                            ctx->MatVal[idx][i] = ctx->MX[0][i];
                    }
                    break;

        case 8:     for (i = 1; i <= n; ++i) {
                        if (fabs(ctx->MX[0][i]) > x)
                            ctx->MatVal[idx][i] = 0.0;
                        else
                            ctx->MatVal[idx][i] = ctx->MX[0][i];
                    }
                    break;

        default:    mat_err(ctx, 0);
                    goto MNCFin;
    }
    err = 0;

MNCFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mpz(cmd)      mpz(A,B,P)                                              */
/*                                                                          */
/*  Row permutation for a zero-free diagonal. A must be (n,n). The command  */
/*  creates a (n,n)-matrix B and a (n,1)-vector P. P contains the           */
/*  permutation, B the permuted matrix.                                     */
/*                                                                          */  
/*  The algorithm is adapted from ACM algorithm 575 (I.S. Duff).            */
/*                                                                          */  
/*  Return 0 if OK, -1 if error.                                            */

int m_mpz(TDAContext *ctx, char *cmd)
{
    register int i,j,k;
    int err,row,col,idx,idx1,ivflg,n,nn,r;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 4,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MPZFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MPZFin;
    if (row != col) {
        mat_err(ctx, 1);
        goto MPZFin;
    }
    if ((p = m_getmat(ctx, p,row,col,&idx,cmd,0)) == NULL)
        goto MPZFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MPZFin;
    if ((p = m_getmat(ctx, p,row,1,&idx1,cmd,1)) == NULL)
        goto MPZFin;

    nn = 0;
    n = row * col;
    r = row + 1;

    for (i = 1; i <= n; ++i) {
        if (ctx->MX[0][i] != 0.0)
            nn++;
    }
    if (alloc_acn(ctx, nn + 1))
        goto MPZFin;
    if (alloc_aci(ctx, row + 1))
        goto MPZFin;
    if (alloc_acj(ctx, row + 1))
        goto MPZFin;
    if (alloc_acm(ctx, row + 1))
        goto MPZFin;
    if (alloc_acr(ctx, 4 * r + 1))
        goto MPZFin;

    k = 0;
    for (i = 0; i < row; ++i) {
        ctx->AcI[i + 1] = k + 1;
        n = 0;
        for (j = 1; j <= col; ++j) {
            if (ctx->MX[0][i * col + j] != 0.0) {
                ctx->AcN[++k] = j;
                n++;
            }
        }
        ctx->AcJ[i + 1] = n;
    }
    n = m_rperm(ctx, row,ctx->AcN,nn,ctx->AcI,ctx->AcJ,ctx->AcM,ctx->AcR,ctx->AcR + r,ctx->AcR + 2 * r,ctx->AcR + 3 * r);

    printf2(ctx, "Number of non-zeros: %d\n",n);

    for (i = 1; i <= row; ++i) {
        k = ctx->AcM[i];
        ctx->MatVal[idx1][i] = (double)k;          
        for (j = 1; j <= col; ++j)
            ctx->MatVal[idx][(i - 1) * col + j] = ctx->MX[0][(k - 1) * col + j];
    }
    err = 0;

MPZFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mpb(cmd,opt)  Permutation to block triangular form. The algorithm     */
/*                  is adapted from CACM 529 (I.S. Duff and J.K. Reid).     */
/*                                                                          */
/*  opt = 0 :  mpbl(A,B,P,N,U)    lower block diagonal                      */
/*  opt = 1 :  mpbu(A,B,P,N,U)    upper block diagonal                      */
/*                                                                          */
/*  A must be (n,n). The command creates a (n,n)-matrix B and a (n,1)-      */
/*  vector P. B contains the permuted matrix, P the permutation.            */
/*  N is (1,1) and contains the number of blocks. U is (k,1) with k =       */
/*  number of blocks. U(i,1) contains the row number where the i.th block   */
/*  begins.                                                                 */
/*                                                                          */  
/*  Return 0 if OK, -1 if error.                                            */

int m_mpb(TDAContext *ctx, char *cmd,int opt)
{
    register int i,j,k,l;
    int err,row,col,idx,idx1,idx2,idx3,ivflg,n,nn,r;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 5,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MPBFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MPBFin;
    if (row != col) {
        mat_err(ctx, 1);
        goto MPBFin;
    }
    if ((p = m_getmat(ctx, p,row,col,&idx,cmd,0)) == NULL)
        goto MPBFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MPBFin;

    if ((p = m_getmat(ctx, p,row,1,&idx1,cmd,0)) == NULL)
        goto MPBFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MPBFin;

    if ((p = m_getmat(ctx, p,1,1,&idx2,cmd,0)) == NULL)
        goto MPBFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MPBFin;

    nn = 0;
    n = row * col;
    r = row + 1;

    for (i = 1; i <= n; ++i) {
        if (ctx->MX[0][i] != 0.0)
            nn++;
    }
    if (alloc_acn(ctx, nn + 1))
        goto MPBFin;
    if (alloc_aci(ctx, row + 1))
        goto MPBFin;
    if (alloc_acj(ctx, row + 1))
        goto MPBFin;
    if (alloc_acm(ctx, row + 1))
        goto MPBFin;
    if (alloc_acs(ctx, row + 1))
        goto MPBFin;
    if (alloc_acr(ctx, 3 * r + 1))
        goto MPBFin;

    k = 0;
    for (i = 0; i < row; ++i) {
        ctx->AcI[i + 1] = k + 1;
        n = 0;
        for (j = 1; j <= col; ++j) {
            if (opt == 0) {
                if (ctx->MX[0][i * col + j] != 0.0) {
                    ctx->AcN[++k] = j;
                    n++;
                }
            }
            else {
                if (ctx->MX[0][(j - 1) * col + i + 1] != 0.0) {
                    ctx->AcN[++k] = j;
                    n++;
                }
            }
        }
        ctx->AcJ[i + 1] = n;
    }

    n = m_perm(ctx, row,ctx->AcN,nn,ctx->AcI,ctx->AcJ,ctx->AcM,ctx->AcS,ctx->AcR,ctx->AcR + r,ctx->AcR + 2 * r);

    printf2(ctx, "Number of blocks: %d\n",n);
    
    if ((p = m_getmat(ctx, p,n,1,&idx3,cmd,1)) == NULL)
        goto MPBFin;

    ctx->MatVal[idx2][1] = (double)n;

    for (i = 1; i <= n; ++i)  
        ctx->MatVal[idx3][i] = (double)ctx->AcS[i];

    for (i = 1; i <= row; ++i) {
        k = ctx->AcM[i];
        ctx->MatVal[idx1][i] = (double)k;          

        for (j = 1; j <= col; ++j) {
            l = ctx->AcM[j];
            ctx->MatVal[idx][(i - 1) * col + j] = ctx->MX[0][(k - 1) * col + l];
        }
    }
    err = 0;

MPBFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mevs(cmd)     mevs(A,E,V)                                             */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mevs(TDAContext *ctx, char *cmd)
{
    int err,row,col,idx,idx1,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 5,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MEVSFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MEVSFin;
    if (row != col) {
        mat_err(ctx, 1);
        goto MEVSFin;
    }
    if ((p = m_getmat(ctx, p,row,1,&idx,cmd,0)) == NULL)
        goto MEVSFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MEVSFin;
    if ((p = m_getmat(ctx, p,row,col,&idx1,cmd,1)) == NULL)
        goto MEVSFin;

    m_copy(ctx, ctx->MatVal[idx1],ctx->MX[0],row * col);
    err = evecf(ctx, row,ctx->MatVal[idx],ctx->MatVal[idx1]);
    if (err == -2) {
        mat_err(ctx, 3);
        goto MEVSFin;
    }
    else if (err) {  
        m_cmdmsg(ctx);
        printf1(ctx, "Error: no success in calculating eigenvalues.\n");
        goto MEVSFin;
    }
    err = 0;

MEVSFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mev(cmd)     mev(A,ER,EI,EVR,EVI)    (Grad/Brebner)                   */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mev(TDAContext *ctx, char *cmd)
{
    int err,row,col,idx,idx1,idx2,idx3,ivflg,indica;
    int *indic;
    register char *p;

    err = -1;
    indica = ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 4,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MEVFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MEVFin;
    if (row != col) {
        mat_err(ctx, 1);
        goto MEVFin;
    }
    if ((p = m_getmat(ctx, p,row,1,&idx,cmd,0)) == NULL)
        goto MEVFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MEVFin;
    if ((p = m_getmat(ctx, p,row,1,&idx1,cmd,0)) == NULL)
        goto MEVFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MEVFin;
    if ((p = m_getmat(ctx, p,row,col,&idx2,cmd,0)) == NULL)
        goto MEVFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MEVFin;
    if ((p = m_getmat(ctx, p,row,col,&idx3,cmd,1)) == NULL)
        goto MEVFin;

    if (!(indic = (int *)calloc((size_t)(row + 1),sizeof(int)))) {
        mat_err(ctx, 3);
        goto MEVFin;
    }                
    memrq(ctx, row + 1,sizeof(int));
    indica = row + 1;

    err = eigen1(ctx, row,ctx->MX[0],ctx->MatVal[idx],ctx->MatVal[idx1],ctx->MatVal[idx2],ctx->MatVal[idx3],indic);
     
    if (err > 0) {
        m_cmdmsg(ctx);
        printf1(ctx, "Error: no success in eigenvalue/vector calculations.\n");
        goto MEVFin;
    }
    else if (err < 0) {
        mat_err(ctx, 3);
        goto MEVFin;
    }
    err = 0;

MEVFin:
    if (indica > 0)
        free((char *)indic);
    memrq(ctx, -indica,sizeof(int));
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_ginv(cmd)     mginv(A,R)                                              */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_ginv(TDAContext *ctx, char *cmd)
{
    int err,row,col,idx,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 6,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MGINVFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MGINVFin;
    if (row < col) {
        mat_err(ctx, 2);
        goto MGINVFin;
    }
    if ((p = m_getmat(ctx, p,col,row,&idx,cmd,1)) == NULL)
        goto MGINVFin;

    m_copy(ctx, ctx->MatVal[idx],ctx->MX[0],row * col);
    err = ginv(ctx, row,col,ctx->MatVal[idx]);
    if (err < 0) {
        if (err == -1)
            mat_err(ctx, 3);
        else           
            mat_err(ctx, 2);
        goto MGINVFin;
    }
    printf1(ctx, "Pseudorank of %s: %d\n",ctx->MXName,err);
    err = 0;

MGINVFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_svd(cmd,opt)      opt 0 :  msvd(A,Q)                                  */
/*                      opt 1 :  msvd1(A,Q,U,V)                             */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_svd(TDAContext *ctx, char *cmd,int opt)
{
    int err,row,col,idx,idx1,idx2,no,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = cmd + 5;
    no = 1;
    if (opt) {
        p++;
        no = 0;
    }   
    p = n_mexpr(ctx, p,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MSVDFin;
    if (row < col) {
        mat_err(ctx, 2);
        goto MSVDFin;
    }
    if ((p = m_check1(ctx, p)) == NULL)
        goto MSVDFin;
    if ((p = m_getmat(ctx, p,col,1,&idx,cmd,no)) == NULL)
        goto MSVDFin;

    if (opt) {
        if ((p = m_check1(ctx, p)) == NULL)
            goto MSVDFin;
        if ((p = m_getmat(ctx, p,row,col,&idx1,cmd,0)) == NULL)
            goto MSVDFin;
        if ((p = m_check1(ctx, p)) == NULL)
            goto MSVDFin;
        if ((p = m_getmat(ctx, p,col,col,&idx2,cmd,1)) == NULL)
            goto MSVDFin;

        err = svdecomp(ctx, row,col,ctx->MX[0],ctx->MatVal[idx],ctx->MatVal[idx1],ctx->MatVal[idx2],3);
    }
    else {
        if (alloc_acu(ctx, row * col + 1))
            goto MSVDFin;
        if (alloc_acv(ctx, col * col + 1))
            goto MSVDFin;
        err = svdecomp(ctx, row,col,ctx->MX[0],ctx->MatVal[idx],ctx->AcU,ctx->AcV,0);    
        alloc_acu(ctx, 0);
        alloc_acv(ctx, 0);
    }
    if (err < 0) {
        if (err == -1)
            mat_err(ctx, 3);
        else           
            mat_err(ctx, 2);
        goto MSVDFin;
    }
    err = 0;

MSVDFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_wvec(cmd)              mwvec(A,W,R)                                   */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_wvec(TDAContext *ctx, char *cmd)
{
    register int i,j;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp,tmp1,x0,x1;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 6,2,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MWVECFin;
    if (col != 1) {
        mat_err(ctx, 7);
        goto MWVECFin;
    }
    if ((p = m_check1(ctx, p)) == NULL)
        goto MWVECFin;
    if ((p = m_getmat(ctx, p,row,1,&idx,cmd,1)) == NULL)
        goto MWVECFin;

    for (i = 1; i < row; ++i) {
        tmp = tmp1 = 0.0;
        for (j = i + 1; j <= row; ++j) {
            x0 = get_mxval(ctx, 0,j,1);
            x1 = get_mxval(ctx, 1,j,1);
            tmp += x0 * x1;                         
            tmp1 += x1;
        }
        if (tmp1 != 0.0)  
            tmp /= tmp1;
        else
            tmp = get_mxval(ctx, 0,i,1);
        ctx->MatVal[idx][i] = tmp;
    }
    ctx->MatVal[idx][row] = get_mxval(ctx, 0,row,1);
    err = 0;

MWVECFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_wvec1(cmd)            mwvec1(A,W,T,R)                                 */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_wvec1(TDAContext *ctx, char *cmd)
{
    register int i,j;
    int err,row,col,idx,ivflg;
    register char *p;
    double a,tmp,tmp1,tmp2,x0,x1;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 7,3,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MWVEC1Fin;
    if (col != 1) {
        mat_err(ctx, 7);
        goto MWVEC1Fin;
    }
    if ((p = m_check1(ctx, p)) == NULL)
        goto MWVEC1Fin;
    if ((p = m_getmat(ctx, p,row,1,&idx,cmd,1)) == NULL)
        goto MWVEC1Fin;

    for (i = 1; i <= row; ++i) {
       
        tmp = tmp1 = 0.0;
        a = get_mxval(ctx, 2,i,1);
        for (j = 1; j <= row; ++j) {
            tmp2 = get_mxval(ctx, 2,j,1);
            if (tmp2 <= a)
                continue;
            x0 = get_mxval(ctx, 0,j,1);
            x1 = get_mxval(ctx, 1,j,1);
            tmp += x0 * x1;
            tmp1 += x1;
        }
        if (tmp1 != 0.0)  
            tmp /= tmp1;
        else
            tmp = get_mxval(ctx, 0,i,1);
        ctx->MatVal[idx][i] = tmp;
    }
    err = 0;

MWVEC1Fin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_scal1(cmd)            mscal1(A,R)                                     */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_scal1(TDAContext *ctx, char *cmd)
{
    register int i,j;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 7,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MSCAL1Fin;
        
    if ((p = m_check1(ctx, p)) == NULL)
        goto MSCAL1Fin;
    if ((p = m_getmat(ctx, p,row,col,&idx,cmd,1)) == NULL)
        goto MSCAL1Fin;

    tmp = 0.0;
    for (i = 1; i <= row; ++i) {
        for (j = 1; j <= col; ++j)  
            tmp += get_mxval(ctx, 0,i,j);
    }
    if (tmp != 0.0) {
        for (i = 1; i <= row; ++i) {
            for (j = 1; j <= col; ++j)  
                ctx->MatVal[idx][(i - 1) * col + j] = get_mxval(ctx, 0,i,j) / tmp;
            }
    }
    err = 0;

MSCAL1Fin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mdefc(cmd)            mdefc(m,n,x,A)                                  */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mdefc(TDAContext *ctx, char *cmd)
{
    register int i;
    int err,row,col,idx,m,n,ivflg;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 6,3,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MDEFCFin;
    if (row != 1 || col != 1) {
        mat_err(ctx, 8);
        goto MDEFCFin;
    }
    m = (int)get_mxval(ctx, 0,1,1);
    n = (int)get_mxval(ctx, 1,1,1);
    if (m < 1 || n < 1) {
        mat_err(ctx, 2);
        goto MDEFCFin;
    }
    if ((p = m_check1(ctx, p)) == NULL)
        goto MDEFCFin;
    if ((p = m_getmat(ctx, p,m,n,&idx,cmd,1)) == NULL)
        goto MDEFCFin;

    m *= n;
    tmp = get_mxval(ctx, 2,1,1);
    for (i = 1; i <= m; ++i)
        ctx->MatVal[idx][i] = tmp;
    err = 0;

MDEFCFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mdefi(cmd)            mdefi(m,n,A)                                    */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mdefi(TDAContext *ctx, char *cmd)
{
    register int i,j;
    int err,row,col,idx,m,n,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 6,2,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MDEFIFin;
    if (row != 1 || col != 1) {
        mat_err(ctx, 8);
        goto MDEFIFin;
    }
    m = (int)get_mxval(ctx, 0,1,1);
    n = (int)get_mxval(ctx, 1,1,1);
    if (m < 1 || n < 1) {
        mat_err(ctx, 2);
        goto MDEFIFin;
    }
    if ((p = m_check1(ctx, p)) == NULL)
        goto MDEFIFin;
    if ((p = m_getmat(ctx, p,m,n,&idx,cmd,1)) == NULL)
        goto MDEFIFin;

    j = imin(ctx, m,n);
    for (i = 0; i < j; ++i)  
        ctx->MatVal[idx][i * n + i + 1] = 1.0;
    err = 0;

MDEFIFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_invs(cmd)             minvs(A,R)                                      */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_invs(TDAContext *ctx, char *cmd)
{
    int err,row,col,idx,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 6,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MINVSFin;
    if (row != col) {
        mat_err(ctx, 1);
        goto MINVSFin;
    }
    if ((p = m_check1(ctx, p)) == NULL)
        goto MINVSFin;
    if ((p = m_getmat(ctx, p,row,col,&idx,cmd,1)) == NULL)
        goto MINVSFin;

    m_copy(ctx, ctx->MatVal[idx],ctx->MX[0],row * col);
    col = syminv(ctx, row,ctx->MatVal[idx]);
    if (col < 0) {
        mat_err(ctx, 3);
        goto MINVSFin;
    }
    if (col < row) {
        m_cmdmsg(ctx);
        printf1(ctx, "%s not positive definite.\n",ctx->MXName);
        goto MINVSFin;
    }
    err = 0;

MINVSFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_invd(cmd)             minvd(A,R)                                      */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_invd(TDAContext *ctx, char *cmd)
{
    register int i,k;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 6,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MINVDFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MINVDFin;
    if ((p = m_getmat(ctx, p,row,col,&idx,cmd,1)) == NULL)
        goto MINVDFin;

    k = imin(ctx, row,col);
    for (i = 1; i <= k; ++i) {
        tmp = get_mxval(ctx, 0,i,i);
        if ((tmp) != 0.0)
            ctx->MatVal[idx][(i - 1) * col + i] = 1.0 / tmp;
    }
    err = 0;

MINVDFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_ple(cmd)             mple(T,C,F,D)                                    */
/*                                                                          */
/*  T is n-vector with values (times), C is n-vector with censoring         */
/*  information. T[i] = 1 if censored, T[i] = 0 if not censored!            */
/*  Return n-vector F containing the distribution function, and n-vector    */
/*  D with jumps.                                                           */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_ple(TDAContext *ctx, char *cmd)
{
    register int i,j,k;
    int err,m,row,col,row1,col1,idx,idx1,ivflg;
    register char *p;
    double tmp,tmp1;
    
    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    if ((p = n_mexpr(ctx, cmd + 5,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MPLEFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MPLEFin;
    if ((p = n_mexpr(ctx, p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
        goto MPLEFin;
    if (col != 1 || col1 != 1) {
        mat_err(ctx, 7);
        goto MPLEFin;
    }
    if (row != row1) {
        mat_err(ctx, 2);
        goto MPLEFin;
    }
    if ((p = m_check1(ctx, p)) == NULL)
        goto MPLEFin;
    if ((p = m_getmat(ctx, p,row,1,&idx,cmd,0)) == NULL)
        goto MPLEFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MPLEFin;
    if ((p = m_getmat(ctx, p,row,1,&idx1,cmd,1)) == NULL)
        goto MPLEFin;

    if (alloc_acn(ctx, row + 1))
        goto MPLEFin;
    if (alloc_acns(ctx, row + 1))
        goto MPLEFin;
    if (alloc_actmp(ctx, row + 1))
        goto MPLEFin;
    
    for (i = 1; i <= row; ++i) {
        if (ctx->MX[1][i] == 0.0)
            ctx->AcNS[i - 1] = 1;                /* not censored */
    }

    if (sortdp2a(ctx, row,ctx->MX[0] + 1,ctx->AcNS,ctx->AcN,0))      /* sort */
        return(-1);    

    tmp = 1.0 / (double)row;
    m = row - 1;
    for (i = 0; i < row; ++i) {
        j = ctx->AcN[i];
        ctx->AcTmp[j] += tmp;
        if (ctx->AcNS[j] == 0 && i < row - 1) {
            tmp1 = ctx->AcTmp[j] / (double)m;
            for (k = i + 1; k < row; ++k)  
                ctx->AcTmp[ctx->AcN[k]] += tmp1;
            ctx->AcTmp[j] = 0.0;
        }
        m--;
    }
    tmp = 0.0;
    for (i = 0; i < row; ++i) {
        k = ctx->AcN[i];
        tmp += ctx->AcTmp[k];
        ctx->MatVal[idx][k + 1] = tmp;
        ctx->MatVal[idx1][k + 1] = ctx->AcTmp[k];
    }
    err = 0;

MPLEFin:
    alloc_actmp(ctx, 0);
    alloc_acns(ctx, 0);
    alloc_acn(ctx, 0);
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_nvar(cmd)             mnvar(X)    Replace data matrix by X.           */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_nvar(TDAContext *ctx, char *cmd)
{
    register int i,j;
    int err,row,col,ivflg;
    register char *p;
    char vname1[4 * VNLMax+1];

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    if (ctx->DMDef) {
        m_cmdmsg(ctx);
        printf1(ctx, "Error: a data matrix already exists.\n");
        err = 0;
        goto MNVARFin;
    }
    p = n_mexpr(ctx, cmd + 6,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MNVARFin;
    if ((p = m_check2(ctx, p,0)) == NULL)
        goto MNVARFin;

    if (col > ctx->MaxNV) {
        m_cmdmsg(ctx);
        printf1(ctx, "Error: %d variables exceeds current maximum.\n",col);
        goto MNVARFin;
    }
    for (i = 0; i < col; ++i) {
        if (snprintf(vname1,sizeof(vname1),"%s%d<8>[%d.%d]=M:%s%d",ctx->MXName,i+1,ctx->PMMFmt1,ctx->PMMFmt2,ctx->MXName,i+1) >= (int)sizeof(vname1)) {
            m_cmdmsg(ctx);
            printf1(ctx, "Error: variable definition is too long.\n");
            goto MNVARFin;
        }
        if (save_var(ctx, vname1,0)) {       /* save variable definition */
            m_cmdmsg(ctx);
            printf1(ctx, "Error: can't save variable definitions.\n");
            goto MNVARFin;              
        }
    }
    /** prn_var(0); **/  /* print list of new variables */

    ctx->VLabelLen = 0;
    ctx->NOCMaxA = ctx->NOCDM = ctx->NOC = row;

    /* allocate memory for data */

    if (alloc_vdat(ctx, 0,1)) {
        m_cmdmsg(ctx);
        p_err(ctx, -2,1);
        goto MNVARFin;
    }
    for (i = 0; i < row; ++i) {
        for (j = 0; j < ctx->NVAR; ++j)  
            put_data(ctx, ctx->MX[0][i * col + j + 1],j,i);
    }
    ctx->DMDef = 1;
    printf2(ctx, "Created data matrix with %d variables and %d cases.\n",ctx->NVAR,ctx->NOCDM);
    err = 0;

MNVARFin:
    mx_free(ctx);
    if (err) {
        clear_dm(ctx);            /* clear all variables */
        ctx->NVAR = ctx->NOCMaxA = ctx->NOCDM = ctx->NOC = 0;
    }
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_print(cmd)        mpr[1][a](X [,string] ) [=fname]                    */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_print(TDAContext *ctx, char *cmd)
{
    int err,row,col,aflg,bflg,ivflg;
    register char *p,*q,*s;

    err = -1;
    ivflg = 1;

    printf2(ctx, "%s\n",cmd);
    bflg = aflg = 0;
    p = cmd + 3;
    if (*p == '1') {
        p++;
        bflg = 1;
    }
    if (*p == 'a') {
        p++;
        aflg = 1;
    }
    q = skip_blev(ctx, p);
    if (*--q != ')') {
        mat_err(ctx, 0);
        goto MPRFin;
    }
    p = n_mexpr(ctx, ++p,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MPRFin;

    ivflg = ctx->MXIV[0];
    if (ivflg && bflg)
        ivflg = 2;

    s = NULL;
    if (*p == ',') {
        p++;
        s = p;
        *q = '\0';
    }
    p = q + 1;
    if (*p++ == '=' && *p)
        err = mpr(ctx, row,col,ctx->MX[0],ivflg,ctx->MX1[0],ctx->PMATFmtS,p,s,aflg);
    else
        err = mpr(ctx, row,col,ctx->MX[0],ivflg,ctx->MX1[0],ctx->PMATFmtS,NULL,s,aflg);

MPRFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mkp(cmd)      mkp(A,B,R)      Kronecker product                       */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mkp(TDAContext *ctx, char *cmd)
{
    register int i,j,k,l;
    int err,row,col,row1,col1,idx,n,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    if ((p = n_mexpr(ctx, cmd + 4,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MKPFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MKPFin;
    if ((p = n_mexpr(ctx, p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
        goto MKPFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MKPFin;
    if ((p = m_getmat(ctx, p,row * row1,col * col1,&idx,cmd,1)) == NULL)
        goto MKPFin;
    n = col * col1;
    for (i = 0; i < row; ++i) {
        for (l = 0; l < row1; ++l) {
            for (j = 0; j < col; ++j) {
                for (k = 0; k < col1; ++k) {
                    ctx->MatVal[idx][(i * row1 +l) * n + j * col1 + k + 1] =
                        get_mxval(ctx, 0,i + 1,j + 1) * get_mxval(ctx, 1,l + 1,k + 1);
                }
            }
        }
    }
    err = 0;

MKPFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mls(cmd,opt)   opt=0:  mls (S,R)         least squares                */
/*                   opt=1:  mlse(S,R)                                      */
/*                   opt=2:  mlsi(S,R)                                      */
/*                   opt=3:  mlsei(S,me,mi,R)                               */
/*                   opt=4   mnls(S,me,l,R)                                 */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mls(TDAContext *ctx, char *cmd,int opt)
{
    int err,row,col,idx,ra,re,ivflg,me = 0,mi = 0,ma = 0,row1,col1;
    register char *p;
    double rne,rnl;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = cmd + 4;
    if (opt)
        p++;
    if (opt == 3)
        p++;
    if ((p = n_mexpr(ctx, p,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MLSFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MLSFin;
    if (col < 2) {
        mat_err(ctx, 2);
        goto MLSFin;
    }
    if (opt == 3 || opt == 4) {
        if ((p = n_mexpr(ctx, p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
            goto MLSFin;
        if ((p = m_check1(ctx, p)) == NULL)
            goto MLSFin;
        if (row1 != 1 || col1 != 1) {
            mat_err(ctx, 8);
            goto MLSFin;
        }
        me = (int)ctx->MX[1][1];

        if ((p = n_mexpr(ctx, p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
            goto MLSFin;
        if ((p = m_check1(ctx, p)) == NULL)
            goto MLSFin;
        if (row1 != 1 || col1 != 1) {
            mat_err(ctx, 8);
            goto MLSFin;
        }
        mi = (int)ctx->MX[1][1];

        if (me < 0 || mi < 0) {
            mat_err(ctx, 14);
            goto MLSFin;
        }
        if (opt == 3) {
            ma = row - me - mi;
            if (ma < 0) {
                mat_err(ctx, 14);
                goto MLSFin;
            }
        }
        else if (opt == 4) {
            ma = row - me;
            if (ma < 0) {
                mat_err(ctx, 14);
                goto MLSFin;
            }
            if (mi > col - 1)
                mi = col - 1;
        }
    }

    if ((p = m_getmat(ctx, p,col - 1,1,&idx,cmd,1)) == NULL)
        goto MLSFin;

    if (opt == 1)         /* mlse */
        err = lsei(ctx, ctx->MX[0],0,row,0,col - 1,ctx->MatVal[idx],0,&rne,&rnl,&ra,&re);
                  
    else if (opt == 2)    /* mlsi */
        err = lsei(ctx, ctx->MX[0],0,0,row,col - 1,ctx->MatVal[idx],0,&rne,&rnl,&ra,&re);
                  
    else if (opt == 3)    /* mlsei */
        err = lsei(ctx, ctx->MX[0],me,ma,mi,col - 1,ctx->MatVal[idx],0,&rne,&rnl,&ra,&re);
                  
    else if (opt == 4)    /* mnls */
        err = wnnls(ctx, ctx->MX[0],col,me,ma,col - 1,mi,ctx->MatVal[idx],&rnl);
                  
    else                  /* mls */
        err = lsei(ctx, ctx->MX[0],0,row,0,col - 1,ctx->MatVal[idx],0,&rne,&rnl,&ra,&re);
                  
                  
    switch (err) {
        case  0: break;
        case -1: m_cmdmsg(ctx);
                 printf1(ctx, "Equality constraints contradictory.\n");
                 goto MLSFin;
        case -2: m_cmdmsg(ctx);
                 printf1(ctx, "Cannot satisfy inequality constraints.\n");
                 goto MLSFin;
        case -3: m_cmdmsg(ctx);
                 printf1(ctx, "Cannot satisfy equality and inequality constraints.\n");
                 goto MLSFin;
        case -4: mat_err(ctx, 3);
                 goto MLSFin;
        default: m_cmdmsg(ctx);
                 printf1(ctx, "Undefined error (%d).\n",err);
                 goto MLSFin;
    }
    if (opt != 4)
        printf1(ctx, "Rank of left-hand side: %d\n",ra);
    printf2(ctx, "Norm of residuals: ");          
    rt_printf2_d(ctx, ctx->PMATFmtS,rnl);          
    printf2(ctx, "\n");          
#ifdef TDA_R_PACKAGE
    /* rank and residual norm, each as a 1x1 -- unconditional on
       SILENTFlg, unlike the printf2() line above. No rank for mnls
       (opt==4): lsei() never runs for it, so ra is never set. */
    if (opt != 4) {
        tda_export_cell(ctx, "mls.rank", (double)ra);
        tda_export_endrow(ctx, "mls.rank");
    }
    tda_export_cell(ctx, "mls.rnorm", rnl);
    tda_export_endrow(ctx, "mls.rnorm");
#endif
    err = 0;

MLSFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mlp(cmd,opt)  opt=0:  mlp (T,X,Y)     linear programming              */
/*                  opt=1:  mlp1(T,p,X,Y)                                   */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mlp(TDAContext *ctx, char *cmd,int opt)
{
    int err,row,col,row1,col1,idx,idx1,np,r,ivflg;
    register char *p;
    double tmp;

    ivflg = np = 0;
    err = -1;
    printf2(ctx, "%s\n",cmd);
    p = cmd + 4;
    if (opt)
        p++;
    if ((p = n_mexpr(ctx, p,1,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto MLPFin;
    if (row < 2 || col < 2) {
        mat_err(ctx, 2);
        goto MLPFin;
    }
    if (opt) {
        if ((p = m_check1(ctx, p)) == NULL)
            goto MLPFin;
        if ((p = n_mexpr(ctx, p,1,1,&row1,&col1,&ivflg,NULL)) == NULL) 
            goto MLPFin;
        if (row1 != 1 || col1 != 1 || (np = (int)ctx->MX[1][1]) < 0 || np >= row) {
            m_cmdmsg(ctx);
            printf1(ctx, "Error in second argument.\n");
            goto MLPFin;
        }
    }
    if ((p = m_check1(ctx, p)) == NULL)
        goto MLPFin;
    if ((p = m_getmat(ctx, p,col - 1,1,&idx,cmd,0)) == NULL)
        goto MLPFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MLPFin;
    if ((p = m_getmat(ctx, p,row - 1,1,&idx1,cmd,1)) == NULL)
        goto MLPFin;

    r = lpf1(ctx, row-1,col-1,np,ctx->MX[0],ctx->MatVal[idx],ctx->MatVal[idx1],&tmp,1.e-8);
    if (r) {
        m_cmdmsg(ctx);
        switch (r) {
            case -2:    printf1(ctx, "No solution.\n");
                        if (np >= 2) 
                            printf1(ctx, "(Possibly linear dependent equality constraints.)\n");
                        break;
            case -3:    printf1(ctx, "No solution. Dual objective function is unbounded.\n");
                        break;
            case -4:    printf1(ctx, "No solution. Primal objective function is unbounded.\n");
                        break;
            case -5:    printf1(ctx, "Error: insufficient memory.\n");
                        break;
            default:    printf1(ctx, "Error (%d).\n",r);
                        break;
        }   
        goto MLPFin;
    }
    printf2(ctx, "Value is: ");
    rt_printf2_d(ctx, ctx->PMATFmtS,tmp);
    printf2(ctx, "\n");
#ifdef TDA_R_PACKAGE
    /* the objective value at the returned solution, as a 1x1 --
       unconditional on SILENTFlg, unlike the printf2() lines above */
    tda_export_cell(ctx, "mlp.value", tmp);
    tda_export_endrow(ctx, "mlp.value");
#endif
    err = 0;

MLPFin:
    mx_free(ctx);
    return(err);
}
/*--------------------------------------------------------------------------*/
/*  m_mlpi(cmd)     mlpi(A,B,lmax)   0-1 linear programming (CACM 449).     */
/*                                                                          */
/*  The objective was stated nowhere in this source; established against    */
/*  exhaustive search (examples/coverage/mlpiops.cf):                       */
/*                                                                          */
/*      maximize A[1,]x + B[1]  over x in {0,1}^n                           */
/*      subject to A[i,]x >= B[i], i = 2..m.                                */
/*                                                                          */
/*  Integer coefficients, negative ones handled; B[1] is a constant added   */
/*  to the reported Value; up to lmax optimal solutions are printed.        */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mlpi(TDAContext *ctx, char *cmd)
{
    register int i,j;
    int err,row,col,r,ivflg,vmax,m,n,nest;
    register char *p;

    ivflg = 0;
    err = -1;
    printf2(ctx, "%s\n",cmd);
    p = cmd + 5;
    if ((p = n_mexpr(ctx, p,1,0,&m,&n,&ivflg,NULL)) == NULL) 
        goto MLPIFin;
    if (m < 2 || n < 2) {
        mat_err(ctx, 2);
        goto MLPIFin;
    }
    if ((p = m_check1(ctx, p)) == NULL)
        goto MLPIFin;
    if ((p = n_mexpr(ctx, p,1,1,&row,&col,&ivflg,NULL)) == NULL) 
        goto MLPIFin;

    if (row != m || col != 1) {
        mat_err(ctx, 2);  
        goto MLPIFin;
    }
    if ((p = m_check1(ctx, p)) == NULL)
        goto MLPIFin;

    if ((p = n_mexpr(ctx, p,1,2,&row,&col,&ivflg,NULL)) == NULL) 
        goto MLPIFin;

    if (row != 1 || col != 1) {
        mat_err(ctx, 8);
        goto MLPIFin;
    }
    nest = (int)ctx->MX[2][1];
          
    if (alloc_acn(ctx, nest * n + 1))
        goto MLPIFin;

    if (alloc_acm(ctx, 4 * m + 6 * n + m * n + 2))
        goto MLPIFin;

    if (alloc_acr(ctx, m * n + 1))
        goto MLPIFin;

    if (alloc_acs(ctx, m + 1))
        goto MLPIFin;

    for (i = 1; i <= m; ++i) {
        ctx->AcS[i] = (int)(ctx->MX[1][i]);
        for (j = 1; j <= n; ++j)
            ctx->AcR[(i - 1) * n + j] = (int)(ctx->MX[0][(i - 1) * n + j]);
    }
    r = lpi(ctx, n,m,nest,ctx->AcR,ctx->AcS,ctx->AcN,ctx->AcM,&vmax);

    if (r < 1) {
        m_cmdmsg(ctx);
        switch (r) {
            case -1:    printf1(ctx, "inconsistent constraints.\n");
                        break;
            case -2:    printf1(ctx, "exceeded maximal number of solutions.\n");
                        break;
            default:    printf1(ctx, "Error (%d).\n",r);
                        break;
        }   
        goto MLPIFin;
    }
    printf1(ctx, "Value: %d\nNumber of solutions: %d\n",vmax,r);
#ifdef TDA_R_PACKAGE
    /* the optimal value as a 1x1, and one row per optimal 0-1 vector */
    tda_export_cell(ctx, "mlpi.value", (double)vmax);
    tda_export_endrow(ctx, "mlpi.value");
#endif
    for (i = 0; i < r; ++i) {
        for (j = 1; j <= n; ++j) {
            tda_out("%2d ",ctx->AcN[i * n + j]);
#ifdef TDA_R_PACKAGE
            tda_export_cell(ctx, "mlpi.solutions", (double)ctx->AcN[i * n + j]);
#endif
        }
        newline(ctx);
#ifdef TDA_R_PACKAGE
        tda_export_endrow(ctx, "mlpi.solutions");
#endif
    }
    err = 0;

MLPIFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mul(cmd)      mmul(A1,A2,...,R)     multiplication                    */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mul(TDAContext *ctx, char *cmd)
{
    register int i,j,k,l;
    int err,row,col,row1,col1,idx,ivflg;
    register char *p,*q;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    if ((p = n_mexpr(ctx, cmd + 5,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MMULFin;

    while (1) {
        if ((p = m_check1(ctx, p)) == NULL)
            goto MMULFin;

        q = skip_expr(ctx, p);
        if (!*q || *q == ')')
            break;

        /* get next matrix into MX[1] */

        if ((p = n_mexpr(ctx, p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
            goto MMULFin;

        /* MX[0] * MX[1] -> V */

        if (alloc_acu(ctx, row * col1 + 1))
            goto MMULFin;

        k = imax(ctx, col,row1);
        for (i = 1; i <= row; ++i) {
            for (j = 1; j <= col1; ++j) {
                tmp = 0.0;
                for (l = 1; l <= k; ++l)
                    tmp += get_mxval(ctx, 0,i,l) * get_mxval(ctx, 1,l,j);
                ctx->AcU[(i - 1) * col1 + j] = tmp;
            }
        }

        /* copy U into MX[0], set dimension */

        col = col1;
        if (mx_alloc(ctx, 0,row,col))
            goto MMULFin;
        m_copy(ctx, ctx->MX[0],ctx->AcU,row * col1);
    }
    if ((p = m_getmat(ctx, p,row,col,&idx,cmd,1)) == NULL)
        goto MMULFin;
    m_copy(ctx, ctx->MatVal[idx],ctx->MX[0],row * col);
    err = 0;

MMULFin:
    alloc_acu(ctx, 0);
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_cat(cmd,opt)      opt 0 : mcath(A1,A2,...,R)      horizontal concat   */
/*                      opt 1 : mcatv(A1,A2,...,R)      vertical concat     */
/*                      opt 2 : mdcathv(A1,A2,...,R)    direct sum          */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_cat(TDAContext *ctx, char *cmd,int opt)
{
    register int i,j;
    int err,row,col,row1,col1,r,c,idx,ivflg;
    register char *p,*q;

    err = -1;
    ivflg = 0;

    r = 6;
    if (opt == 2)
        r++;       

    printf2(ctx, "%s\n",cmd);
    if ((p = n_mexpr(ctx, cmd + r,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MCATFin;

    while (1) {
        if ((p = m_check1(ctx, p)) == NULL)
            goto MCATFin;

        q = skip_expr(ctx, p);
        if (!*q || *q == ')')
            break;

        /* get next matrix into MX[1] */

        if ((p = n_mexpr(ctx, p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
            goto MCATFin;

        /* MX[0] + MX[1] -> V */

        if (opt == 0) {
            r = imax(ctx, row,row1);
            c = col + col1;
        }
        else if (opt == 1) {
            r = row + row1;
            c = imax(ctx, col,col1);
        }
        else  {
            r = row + row1;
            c = col + col1;
        }
        if (alloc_acu(ctx, r * c + 1))
            goto MCATFin;

        if (opt == 0) {
            for (i = 1; i <= r; ++i) {
                for (j = 1; j <= col; ++j) 
                    ctx->AcU[(i - 1) * c + j] = get_mxval(ctx, 0,i,j);

                for (j = col + 1; j <= c; ++j) 
                    ctx->AcU[(i - 1) * c + j] = get_mxval(ctx, 1,i,j - col);
            }
        }
        else if (opt == 1) {
            for (j = 1; j <= c; ++j) {
                for (i = 1; i <= row; ++i) 
                    ctx->AcU[(i - 1) * c + j] = get_mxval(ctx, 0,i,j);

                for (i = row + 1; i <= r; ++i) 
                    ctx->AcU[(i - 1) * c + j] = get_mxval(ctx, 1,i - row,j);
            }
        }
        else  {
            for (i = 1; i <= row; ++i) {
                for (j = 1; j <= col; ++j) 
                    ctx->AcU[(i - 1) * c + j] = get_mxval(ctx, 0,i,j);
            }
            for (i = row + 1; i <= r; ++i) {
                for (j = col + 1; j <= c; ++j) 
                    ctx->AcU[(i - 1) * c + j] = get_mxval(ctx, 1,i - row,j - col);
            }
        }

        /* copy U into MX[0], set dimension */

        row = r;   
        col = c;
        if (mx_alloc(ctx, 0,row,col))
            goto MCATFin;
        m_copy(ctx, ctx->MX[0],ctx->AcU,row * col);
    }
    if ((p = m_getmat(ctx, p,row,col,&idx,cmd,1)) == NULL)
        goto MCATFin;
    m_copy(ctx, ctx->MatVal[idx],ctx->MX[0],row * col);
    err = 0;

MCATFin:
    alloc_acu(ctx, 0);
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_srow(cmd,opt)  opt=0: msrow(X,D,R)    select rows                     */
/*                   opt=1: mscol(X,D,R)    select columns                  */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_srow(TDAContext *ctx, char *cmd,int opt)
{
    register int i,j,k,l;
    int err,row,col,row1,col1,m,idx,ivflg,sflag;             
    register char *p,*q;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    if ((p = n_mexpr(ctx, cmd + 6,1,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto MSROWFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MSROWFin;

    if (*p == '<') {
        sflag = 1;
        q = p;
        if (m_selidx(ctx, q,0,&m,ctx->AcI) == NULL || m < 1) {
            mat_err(ctx, 9);
            goto MSROWFin;
        }
        if (alloc_aci(ctx, m + 1))
            goto MSROWFin;

        if ((p = m_selidx(ctx, p,1,&m,ctx->AcI)) == NULL || m < 1) {
            mat_err(ctx, 9);
            goto MSROWFin;
        }
    }
    else {
        sflag = 0;

        if ((p = n_mexpr(ctx, p,1,1,&row1,&col1,&ivflg,NULL)) == NULL) 
            goto MSROWFin;

        if (row1 != 1 && col1 != 1) {
            mat_err(ctx, 9);
            goto MSROWFin;
        }
        m = imax(ctx, row1,col1);
    }
    if ((p = m_check1(ctx, p)) == NULL)
        goto MSROWFin;

    if (opt == 0) {
        if ((p = m_getmat(ctx, p,m,col,&idx,cmd,1)) == NULL)
            goto MSROWFin;

        l = 0;
        for (i = 1; i <= m; ++i) {
            if (sflag)
                k = ctx->AcI[i - 1];
            else
                k = (int)ctx->MX[1][i];

            if (k < 1 || (opt == 0 && k > row) || (opt == 1 && k > col)) {
                mat_err(ctx, 9);
                goto MSROWFin;
            }
            for (j = 1; j <= col; ++j)
                ctx->MatVal[idx][l * col + j] = get_mxval(ctx, 0,k,j);
            l++; 
        }
    }
    else {
        if ((p = m_getmat(ctx, p,row,m,&idx,cmd,1)) == NULL)
            goto MSROWFin;

        for (j = 1; j <= m; ++j) {
            if (sflag)
                k = ctx->AcI[j - 1];
            else
                k = (int)ctx->MX[1][j];

            if (k < 1 || (opt == 0 && k > row) || (opt == 1 && k > col)) {
                mat_err(ctx, 9);
                goto MSROWFin;
            }
            for (i = 0; i < row; ++i)
                ctx->MatVal[idx][i * m + j] = get_mxval(ctx, 0,i + 1,k);
        }
    }
    err = 0;

MSROWFin:
    alloc_aci(ctx, 0); 
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_selidx(s,opt,n,idx)   s points to a string of the form                */
/*                          <i1,i2,...> or <i1,i2,,i3,...>                  */
/*                                                                          */
/*  Return the number of indices in n. Return pointer to next character     */
/*  after >, or NULL if error.                                              */

char *m_selidx(TDAContext *ctx, char *s,int opt,int *n,int *idx)
{
    register int i,j;
    int nc,nca,ni;

    i = ni = 0;
    while (*s) {
        s++;
        if (sscanf(s,"%d,,%d",&nca,&nc) == 2 && nca >= 1 && nc >= 1) {
            if (nca <= nc) {
                for (j = nca; j <= nc; ++j) {
                    ni++;
                    if (opt)
                        idx[i++] = j;
                }
            }
            else {
                for (j = nca; j >= nc; --j) {
                    ni++;
                    if (opt)
                        idx[i++] = j;
                }
            }
            s = skip_int(ctx, s);
            s = skip_int(ctx, s + 2);
        }
        else if (sscanf(s,"%d",&nc) == 1 && nc >= 1) {
            ni++;
            s = skip_int(ctx, s);
            if (opt) 
                idx[i++] = nc;
        }
        else
            break;         

        if (*s == ',')
            continue;
        else if (*s == '>') {
            *n = ni;
            return(s + 1);
        }
        else
            break;
    }
    return(NULL);
}

/*--------------------------------------------------------------------------*/
/*  m_sort(cmd,opt)  opt=0: msort(X,D,R)                                    */
/*                   opt=1: mrank(X,D,R)                                    */
/*                   opt=2: msort1(X,D,R)                                   */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_sort(TDAContext *ctx, char *cmd,int opt)
{
    register int i,j,k,l,ll;
    int err,row,col,row1,col1,m,idx,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = cmd + 6;
    if (opt == 2)
        p++;
    if ((p = n_mexpr(ctx, p,1,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto MSORTFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MSORTFin;
    if ((p = n_mexpr(ctx, p,1,1,&row1,&col1,&ivflg,NULL)) == NULL) 
        goto MSORTFin;
    if (row1 != 1 && col1 != 1) {
        mat_err(ctx, 9);
        goto MSORTFin;
    }
    m = imax(ctx, row1,col1);
    for (i = 1; i <= m; ++i) {
        k = (int)ctx->MX[1][i];
        if (k < 1 || k > col) {
            mat_err(ctx, 9);
            goto MSORTFin;
        }
    }
    if ((p = m_check1(ctx, p)) == NULL)
        goto MSORTFin;

    if (alloc_acn(ctx, m + 1))        /* save indices for sort */
        goto MSORTFin;

    for (i = 1; i <= m; ++i)
        ctx->AcN[i] = (int)ctx->MX[1][i];

    if (alloc_aci(ctx, row + 1))         /* pointer for sorting */
        goto MSORTFin;

    if (sortdpn(ctx, row,col,ctx->MX[0],m,ctx->AcN,ctx->AcI))
        goto MSORTFin;      

    if (opt == 1) {                                         /* rank */
        if ((p = m_getmat(ctx, p,row,1,&idx,cmd,1)) == NULL)
            goto MSORTFin;

        for (i = 1; i <= row; ++i)
            ctx->MatVal[idx][i] = (double)ctx->AcI[i];
    }
    else {
        if (alloc_actmp(ctx, row * col + 1))     
            goto MSORTFin;

        k = 1;
        for (i = 1; i <= row; ++i) {
            l = (ctx->AcI[i] - 1) * col;
            ll = (k - 1) * col;
            for (j = 1; j <= col; ++j) 
                ctx->AcTmp[++ll] = ctx->MX[0][++l];
                       
            if (opt == 2 && k > 1) {
                k--;
                for (j = 1; j <= col; ++j) {
                    if (ctx->AcTmp[k * col + j] != ctx->AcTmp[(k-1) * col + j]) {
                        k++;
                        break;
                    }
                }   
            }
            k++;
        }
        k--;
        if ((p = m_getmat(ctx, p,k,col,&idx,cmd,1)) == NULL)
            goto MSORTFin;
        m_copy(ctx, ctx->MatVal[idx],ctx->AcTmp,k * col);
    }
    err = 0;

MSORTFin:
    alloc_actmp(ctx, 0);
    alloc_acn(ctx, 0);
    alloc_aci(ctx, 0);
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  get_mselect(p)  get <...> expression into AcC.                          */
/*                                                                          */
/*  Return pointer to next character, or NULL if error                      */

char *get_mselect(TDAContext *ctx, char *p)
{
    int l;
    register char *q;

    l = (int)(strlen(p));
    if (alloc_acc(ctx, l))   
        return(NULL);

    q = ctx->AcC;
    while (*p) { 
        if (*p == '>')
            break;
        *q++ = *p++;
    }
    *q = '\0';
    if (*p++ != '>' || *p++ != ',')
        return(NULL);
    return(p);
}

/*--------------------------------------------------------------------------*/
/*  m_exp(cmd)      mexp(expression,A)                                      */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_exp(TDAContext *ctx, char *cmd)
{
    register int i;
    int err,r,row,col,idx;
    char *p,*q;
    char mtmp[MatDefLen + 1];

    err = -1;
    printf2(ctx, "%s\n",cmd);
    mdefcpy(ctx, mtmp,cmd);

    p = cmd + 5;

    r = (int)(strlen(cmd));
    q = cmd + r;               
    if (*--q != ')') {
        mat_err(ctx, 0);
        goto MEXPFin;
    }
    *q = '\0';

    for (i = r; i >= 5; --i) {
        if (*--q == ',') 
            break;
    }
    if (*q != ',' || !*(q + 1)) {
        mat_err(ctx, 0);
        goto MEXPFin;
    }
    *q++ = '\0';

    if (mex_eval(ctx, p))
        goto MEXPFin;

    row = ctx->MEXOPRow[0];
    col = ctx->MEXOPCol[0];

    if ((idx = mat_newmat(ctx, q,row,col)) < 0)
        goto MEXPFin;

    m_copy(ctx, ctx->MatVal[idx],ctx->MEXOP[0],row * col);
    strcpy(ctx->MatDef[idx],mtmp);
    err = 0;

MEXPFin:       
    alloc_mex_eval(ctx, 0,0,0);                  
    alloc_mex(ctx, 0,0);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_expr(cmd)     mexpr(expression,A)                                     */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_expr(TDAContext *ctx, char *cmd)
{
    int err,row,col,idx,ivflg;
    register char *p;
    char *sel;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);

    sel = NULL;

    p = cmd + 6;
    if (*p == '<') {
        p = get_mselect(ctx, p + 1);
        if (p == NULL) {
            return(-1);
        }
        sel = ctx->AcC;
    }

    if ((p = n_mexpr(ctx, p,1,0,&row,&col,&ivflg,sel)) == NULL) 
        goto MEXPRFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MEXPRFin;
    if ((p = m_getmat(ctx, p,row,col,&idx,cmd,1)) == NULL)
        goto MEXPRFin;

    m_copy(ctx, ctx->MatVal[idx],ctx->MX[0],row * col);
    err = 0;

MEXPRFin:
    alloc_acc(ctx, 0);   
    mx_free(ctx);
    return(err);
}
/*--##----------------------------------------------------------------------*/
/*  m_expr1(cmd)     mexpr1(B,expression,A)                                 */
/*                                                                          */
/*  Evaluate expression in block mode, blocks defined by B.                 */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_expr1(TDAContext *ctx, char *cmd)
{
    register int j,k,k1,l;
    int err,ivflg,row,col,row1,col1,idx;
    register char c,*p,*q;
    double *tmp = NULL;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    if ((p = n_mexpr(ctx, cmd + 7,1,0,&row1,&col1,&ivflg,NULL)) == NULL) 
        goto MEXPR1Fin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MEXPR1Fin;

    q = skip_expr(ctx, p);
    c = *q;
    *q = '\0';
   
    if (get_mexpr(ctx, p,&row,&col,&ivflg) || row < 1 || col < 1) 
        goto MEXPR1Fin;

    row = imax(ctx, row,row1);
    if (row1 != row) {
        mat_err(ctx, 2);
        goto MEXPR1Fin;
    }
    *q = c;
    if ((q = m_check1(ctx, q)) == NULL)
        goto MEXPR1Fin;
    if ((q = m_getmat(ctx, q,row,col,&idx,cmd,1)) == NULL)
        goto MEXPR1Fin;

    if (alloc_actmp(ctx, col1 + 1))
        goto MEXPR1Fin;

    for (j = 1; j <= col1; ++j)
        ctx->AcTmp[j] = ctx->MX[0][j];         

    k1 = 0;
    k = 1;
    while (k <= row) {
        if (k == row)  
            l = 1;
        else {
            l = 0;
            for (j = 1; j <= col1; ++j) {
                if (ctx->AcTmp[j] != ctx->MX[0][k * col1 + j]) {
                    l = 1;
                    break;
                }
            }
        }
        if (l == 0)
            k++;
        else {
            if (eval_mexpr(ctx, k1,k,col,ctx->MatVal[idx],0,tmp,NULL,NULL,-1,0.0))
                goto MEXPR1Fin;

            for (j = 1; j <= col1; ++j)  
                ctx->AcTmp[j] = ctx->MX[0][k * col1 + j];   
            k1 = k;
            k++;
        }
    }
    err = 0;

MEXPR1Fin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_setv(cmd)  Syntax: msetv(expression,X(row,col))                       */
/*               Sets X(row,col) to the value of expression (must be a      */
/*               scalar expression). X must be a matrix name.               */
/*                                                                          */
/*               return 0 if OK, -1 if error.                               */

int m_setv(TDAContext *ctx, char *cmd)  
{
    int err,idx,row,col,ivflg;
    char *p,mname[VNLMax + 1];  
   
    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    if ((p = n_mexpr(ctx, cmd + 6,1,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto MSETVFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MSETVFin;
    if (row != 1 || col != 1) {
        mat_err(ctx, 8);
        goto MSETVFin;       
    }
    if ((p = get_mname(ctx, p,mname,1)) == NULL)
        return(-1);
    if ((idx = mat_getidx(ctx, mname,1)) < 0)
        return(-1);
    if (*p++ != '(') {
        mat_err(ctx, 0);
        goto MSETVFin;
    }
    if (sscanf(p,"%d,%d",&row,&col) == 2) {
        p = skip_int(ctx, p);
        p = skip_int(ctx, ++p);
    }
    else {
        if ((p = n_mexpr(ctx, p,2,1,&row,&col,&ivflg,NULL)) == NULL) 
            goto MSETVFin;
        if (row != 1 || col != 1) {
            mat_err(ctx, 8);
            goto MSETVFin;
        }
        row = (int)ctx->MX[1][1];
        col = (int)ctx->MX[2][1];
    }
    if ((p = m_check2(ctx, p,0)) == NULL)
        goto MSETVFin;
    if ((p = m_check2(ctx, p,1)) == NULL)
        goto MSETVFin;
    if (row < 1 || row > ctx->MatRow[idx] || col < 1 || col > ctx->MatCol[idx]) {
        mat_err(ctx, 10);
        goto MSETVFin;
    }
    ctx->MatVal[idx][(row - 1) * ctx->MatCol[idx] + col] = ctx->MX[0][1];
    err = 0;

MSETVFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_brr(cmd)    Create matrix with BRR indicators.                        */
/*                mbrr(ns,nu,A)                                             */
/*                where ns is number of strata, nu is number of secu.       */
/*                A has dimension ns x nr, nr is number of replications.    */
/*                Return 0 if OK, -1 if error.                              */

int m_brr(TDAContext *ctx, char *cmd)  
{
    register int i,j;
    int err,idx,row,col,ns,nu,nr,ivflg;
    char *p;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    if ((p = n_mexpr(ctx, cmd + 5,2,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto MBRRFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MBRRFin;
    if (row != 1 || col != 1) {
        mat_err(ctx, 8);
        goto MBRRFin;       
    }
    ns = (int)ctx->MX[0][1];
    nu = (int)ctx->MX[1][1];
    if (ctx->SILENTFlg < 0)  
        printf1(ctx, "Number of strata: %d, secu: %d\n",ns,nu);

    err = brr_gen(ctx, ns,nu,&nr,1);
    if (err || nr < 1)
        goto MBRRFin;

    if ((p = m_getmat(ctx, p,ns,nr,&idx,cmd,1)) == NULL)
        goto MBRRFin;

    for (i = 0; i < ns; ++i) {
        for (j = 0; j < nr; ++j)
            ctx->MatVal[idx][i * nr + j + 1] = (double)ctx->AcC[j * ns + i + 1];
    }
    err = 0;

MBRRFin:
    mx_free(ctx);
    alloc_acc(ctx, 0);           /* used by brr_gen(ctx) */
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mtrim(A,ca,ra,cb,rb,R)   Create new matrix.                           */
/*                                                                          */
/*      ca = leading columns,cb trailing columns                            */
/*      ra = leading rows, rb trailing rows.                                */
/*      if positive delete, if negative add zero rows/columns               */

int m_mtrim(TDAContext *ctx, char *cmd)  
{
    register int i,j;
    int err,idx,row,col,row1,col1,ivflg,ca,cb,ra,rb;
    int ia,ib,ii,ja,jb,jj;
    char *p;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    if ((p = n_mexpr(ctx, cmd + 6,1,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto MTRIMFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MTRIMFin;

    if ((p = n_mexpr(ctx, p,4,1,&row1,&col1,&ivflg,NULL)) == NULL) 
        goto MTRIMFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MTRIMFin;

    if (row1 != 1 || col1 != 1) {
        mat_err(ctx, 8);
        goto MTRIMFin;       
    }
    ca = (int)ctx->MX[1][1];
    ra = (int)ctx->MX[2][1];
    cb = (int)ctx->MX[3][1];
    rb = (int)ctx->MX[4][1];

    row1 = row - ra - rb;
    col1 = col - ca - cb;

    if (row1 <= 0 || col1 <= 0) {
        mat_err(ctx, 15);
        goto MTRIMFin;       
    }
    if ((p = m_getmat(ctx, p,row1,col1,&idx,cmd,1)) == NULL)
        goto MTRIMFin;
       

    ia = imax(ctx, 1 + ra,1);
    ib = imin(ctx, row - rb,row);
    ii = imax(ctx, 1 - ra,1);
    for (i = ia; i <= ib; ++i) {
        ja = imax(ctx, 1 + ca,1);
        jb = imin(ctx, col - cb,col);
        jj = imax(ctx, 1 - ca,1);
        for (j = ja; j <= jb; ++j)
            ctx->MatVal[idx][(ii - 1) * col1 + jj++] = ctx->MX[0][(i - 1) * col + j];
        ii++; 
    }
    err = 0;

MTRIMFin:
    mx_free(ctx);
    return(err);
}
/*--------------------------------------------------------------------------*/
/*  m_mmp(cmd,opt)   opt=0: mmp(X,P)      monotone projection (plain PAVA)  */
/*                   opt=1: mmp1(X,T,P)   Kruskal's PRIMARY tie approach:   */
/*                          within each tie group of T, sort X ascending,   */
/*                          PAVA the sequence, scatter back by rank slot    */
/*                          (result monotone up to reordering within ties)  */
/*                   opt=2: mmp2(X,T,P)   SECONDARY approach: each tie      */
/*                          group replaced by its mean, PAVA the means,     */
/*                          broadcast back                                  */
/*                                                                          */
/*  T's column defines tie groups by consecutive equal values (EPSI1).      */
/*  The comment shipped here misnumbered the opts (1/2/3) and called the    */
/*  third form mmp1; corrected against the dispatch in t_mat.c and the      */
/*  hand-verified control case examples/coverage/mmpops.cf.                 */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mmp(TDAContext *ctx, char *cmd,int opt) 
{
    register int i,j,k;
    int err,row,col,row1,col1,idx,ivflg;
    register char *p;
    double tmp,tmp1;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    i = 4;
    if (opt != 0)
        i++;

    p = n_mexpr(ctx, cmd + i,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MMPFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MMPFin;

    if (opt != 0) {
        p = n_mexpr(ctx, p,1,1,&row1,&col1,&ivflg,NULL);
        if (p == NULL)
            goto MMPFin;
        if ((p = m_check1(ctx, p)) == NULL)
            goto MMPFin;

        if (row1 != row || col1 != col) {
            mat_err(ctx, 0);
            goto MMPFin;
        }
        if (alloc_ack(ctx, row + 1))
            goto MMPFin;
        if (alloc_actmp1(ctx, row + 1))
            goto MMPFin;

        if (opt == 1) {
            if (alloc_acj(ctx, row + 1))
                goto MMPFin;
            if (alloc_acn(ctx, row + 1))
                goto MMPFin;
        }   
    }

    if ((p = m_getmat(ctx, p,row,col,&idx,cmd,1)) == NULL)
        goto MMPFin;

    if (alloc_actmp(ctx, row + 1))
        goto MMPFin;
         
    for (j = 1; j <= col; ++j) {

        for (i = 1; i <= row; ++i)  
            ctx->AcTmp[i - 1] = get_mxval(ctx, 0,i,j);
      
        if (opt == 0)
            mreg_proj(ctx, row,ctx->AcTmp);

        else {      
            k = 0;
            ctx->AcK[k] = 1;
            tmp = get_mxval(ctx, 1,1,j);

            for (i = 1; i < row; ++i) {
                tmp1 = get_mxval(ctx, 1,i + 1,j);
                if (fabs(tmp - tmp1) <= ctx->EPSI1)
                    ctx->AcK[k] += 1;            
                else {
                    k++;
                    ctx->AcK[k] = 1;
                }
                tmp = tmp1;
            }
            if (opt == 1) {
                if (mreg_proj1(ctx, row,k + 1,ctx->AcK,ctx->AcJ,ctx->AcN,ctx->AcTmp,ctx->AcTmp1)) {
                    mat_err(ctx, 3);
                    goto MMPFin;
                }
            }
            else if (opt == 2) {
                if (mreg_proj2(ctx, row,k + 1,ctx->AcK,ctx->AcTmp,ctx->AcTmp1)) {
                    mat_err(ctx, 3);
                    goto MMPFin;
                }
            }

        }
        for (i = 0; i < row; ++i)  
            ctx->MatVal[idx][i * col + j] = ctx->AcTmp[i];
        
    }
    err = 0;

MMPFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mpit(cmd,opt)      Iteration of Leslie matrix.                        */
/*                                                                          */
/*  opt = 0: mpit(F,N,n,R)                                                  */
/*  opt = 1: mpit1(F,N,Z,n,R)                                               */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mpit(TDAContext *ctx, char *cmd,int opt)
{
    register int i,j;
    int err,n,m,row,col,idx,ivflg;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);

    if (opt == 0) {
        p = cmd + 5;
        if ((p = n_mexpr(ctx, p,3,0,&row,&col,&ivflg,NULL)) == NULL) 
            goto MPITFin;
        if ((p = m_check1(ctx, p)) == NULL)
            goto MPITFin;

        n = ctx->MXRow[0];
        if (ctx->MXCol[0] != 2 || ctx->MXRow[1] != n || ctx->MXRow[2] != 1 || ctx->MXCol[2] != 1) {
            mat_err(ctx, 0);
            goto MPITFin;
        }
        m = (int)ctx->MX[2][1] + 1;          /* iterations */
    }
    else {
        p = cmd + 6;
        if ((p = n_mexpr(ctx, p,4,0,&row,&col,&ivflg,NULL)) == NULL) 
            goto MPITFin;
        if ((p = m_check1(ctx, p)) == NULL)
            goto MPITFin;

        n = ctx->MXRow[0];
        if (ctx->MXCol[0] != 2 || ctx->MXRow[1] != n || ctx->MXRow[2] != n || ctx->MXRow[3] != 1 || ctx->MXCol[3] != 1) {
            mat_err(ctx, 0);
            goto MPITFin;
        }
        m = (int)ctx->MX[3][1] + 1;          /* iterations */
    }
    if (m < 1)
        m = 1;

    if ((p = m_getmat(ctx, p,m,n,&idx,cmd,1)) == NULL)
        goto MPITFin;

    for (i = 1; i <= n; ++i)
        ctx->MatVal[idx][i] = ctx->MX[1][i];

    for (j = 1; j < m; ++j) {
        tmp = 0.0;
        for (i = 1; i <= n; ++i)  
            tmp += ctx->MX[0][(i - 1) * 2 + 1] * ctx->MatVal[idx][(j - 1) * n + i];
        ctx->MatVal[idx][j * n + 1] = tmp;
        for (i = 2; i <= n; ++i) 
            ctx->MatVal[idx][j * n + i] = ctx->MX[0][(i - 1) * 2] * ctx->MatVal[idx][(j - 1) * n + i - 1];

        if (opt == 1) {
            for (i = 1; i <= n; ++i)
                ctx->MatVal[idx][j * n + i] += ctx->MX[2][i];
        }
    }
    err = 0;

MPITFin:
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mqp(cmd,opt)   opt=0:  mqp (C,D,X)              quadratic programming */
/*                   opt=1:  mqpb(C,D,XL,XU,X)                              */
/*                   opt=2:  mqpc(C,D,A,B,me,X)                             */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mqp(TDAContext *ctx, char *cmd,int opt)
{
    register char *p;
    int i,m,me,err,r,row,col,idx,ivflg,row1,col1,row2,col2,row3,col3;
    int maxit,xla,xua,iaa,nact,*iact;
    double vsmall,diag,*xl = NULL,*xu = NULL,*a = NULL,*b = NULL;

    err = -1;
    me = m = iaa = xla = xua = ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = cmd + 4;
    if (opt)
        p++;

    if ((p = n_mexpr(ctx, p,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MPQFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MPQFin;

    if (row != col) {
        mat_err(ctx, 1);
        goto MPQFin;
    }
    if ((p = n_mexpr(ctx, p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
        goto MPQFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MPQFin;

    if (row1 != row || col1 != 1) {
        mat_err(ctx, 2);
        goto MPQFin;
    }
    if (opt == 1) {
        if ((p = n_mexpr(ctx, p,1,2,&row1,&col1,&ivflg,NULL)) == NULL)
            goto MPQFin;
        if ((p = m_check1(ctx, p)) == NULL)
            goto MPQFin;
        if (row1 != row || col1 != 1) {
            mat_err(ctx, 2);
            goto MPQFin;
        }
        if ((p = n_mexpr(ctx, p,1,3,&row1,&col1,&ivflg,NULL)) == NULL)
            goto MPQFin;
        if ((p = m_check1(ctx, p)) == NULL)
            goto MPQFin;
        if (row1 != row || col1 != 1) {
            mat_err(ctx, 2);
            goto MPQFin;
        }
    }
    else if (opt == 2) {
        if ((p = n_mexpr(ctx, p,1,2,&m,&col1,&ivflg,NULL)) == NULL)
            goto MPQFin;
        if ((p = m_check1(ctx, p)) == NULL)
            goto MPQFin;
        if (col1 != row) {
            mat_err(ctx, 2);
            goto MPQFin;
        }
        if ((p = n_mexpr(ctx, p,1,3,&row2,&col2,&ivflg,NULL)) == NULL)
            goto MPQFin;
        if ((p = m_check1(ctx, p)) == NULL)
            goto MPQFin;
        if (row2 != m || col2 != 1) {
            mat_err(ctx, 2);
            goto MPQFin;
        }
        if ((p = n_mexpr(ctx, p,1,4,&row3,&col3,&ivflg,NULL)) == NULL)
            goto MPQFin;
        if ((p = m_check1(ctx, p)) == NULL)
            goto MPQFin;
        if (row3 != 1 || col3 != 1) {
            mat_err(ctx, 2);
            goto MPQFin;
        }
        me = (int)ctx->MX[4][1];
        if (me < 0 || me > m) {
            printf1(ctx, "Error in number of equality constraints.\n");
            goto MPQFin;
        }
    }     
    if ((p = m_getmat(ctx, p,row,1,&idx,cmd,1)) == NULL)
        goto MPQFin;

    maxit = 50 * (row + m);
    vsmall = 2.0 * ctx->EPSI;

    if (opt == 0 || opt == 2) {
                     
        /* allocate and initialize xl and xu for lower and upper bounds */

        if (!(xl = (double *)calloc((size_t)(row + 1),sizeof(double)))) {
            mat_err(ctx, 3);
            goto MPQFin;
        }                
        xla = row + 1;
        memrq(ctx, xla,sizeof(double));

        if (!(xu = (double *)calloc((size_t)(row + 1),sizeof(double)))) {
            mat_err(ctx, 3);
            goto MPQFin;
        }                
        xua = row + 1;
        memrq(ctx, xua,sizeof(double));

        for (i = 1; i <= row; ++i) {
            xu[i] = ctx->DBLMAX / 1000.0;
            xl[i] = -xu[i];             
        }
    }

    /* allocate iact */

    if (!(iact = (int *)calloc((size_t)(row + 1),sizeof(int)))) {
        mat_err(ctx, 3);
        goto MPQFin;
    }                
    iaa = row + 1;
    memrq(ctx, iaa,sizeof(int));

    switch (opt) {
        case 0: r = qld(ctx, row,0,0,1,maxit,a,b,ctx->MX[1],ctx->MX[0],ctx->MatVal[idx],
                    xl,xu,iact,&nact,vsmall,&diag);
                break;
        case 1:
                r = qld(ctx, row,0,0,1,maxit,a,b,ctx->MX[1],ctx->MX[0],ctx->MatVal[idx],
                    ctx->MX[2],ctx->MX[3],iact,&nact,vsmall,&diag);
                break;

        case 2:
                r = qld(ctx, row,m,me,1,maxit,ctx->MX[2],ctx->MX[3],ctx->MX[1],ctx->MX[0],ctx->MatVal[idx],
                    xl,xu,iact,&nact,vsmall,&diag);
                break;

        default: goto MPQFin;
    }
    /*  tda_out("r=%d nact=%d diag=%g\n",r,nact,diag); */

    switch (r) {
        case  0: break;
        case  1: m_cmdmsg(ctx);
                 printf1(ctx, "Exceeded maximal number of iterations (%d).\n",maxit);
                 goto MPQFin;
        case  2: m_cmdmsg(ctx);
                 printf1(ctx, "Insufficient accuracy for convergence.\n");
                 goto MPQFin;
        case  3: m_cmdmsg(ctx);
                 printf1(ctx, "Insufficient memory.\n");
                 goto MPQFin;
        default: m_cmdmsg(ctx);
                 printf1(ctx, "Inconsistent constraints.\n");
                 goto MPQFin;
    }
    if (diag > 0.0)   
        printf2(ctx, "Coeff matrix was enlarged by %g times unitmatrix.\n",diag);
    err = 0;

MPQFin:
    if (xla > 0) {
        free((char *)xl);
        memrq(ctx, -xla,sizeof(double));
    }
    if (xua > 0) {
        free((char *)xu);
        memrq(ctx, -xua,sizeof(double));
    }
    if (iaa > 0) {
        free((char *)iact);
        memrq(ctx, -iaa,sizeof(int));
    }
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mlsei1(cmd)  mlsei1(S,me,mi,R)                                        */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mlsei1(TDAContext *ctx, char *cmd)
{
    register int i = 0,j = 0,k = 0;
    register char *p = NULL;
    int r = 0,err = 0,row = 0,col = 0,idx = 0,ivflg = 0,n = 0,me = 0,mi = 0,ma = 0,m = 0,row1 = 0,col1 = 0,xxa = 0,xya = 0,aa = 0,ba = 0;
    int maxit = 0,xla = 0,xua = 0,iaa = 0,nact = 0,*iact = NULL;
    double vsmall = 0.0,diag = 0.0,tmp = 0.0,*xl = NULL,*xu = NULL,*xx = NULL,*xy = NULL,*a = NULL,*b = NULL;

    err = -1;
    xla = xua = iaa = aa = ba = xxa = xya = ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = cmd + 7;

    if ((p = n_mexpr(ctx, p,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MLS1Fin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MLS1Fin;
    if (col < 2) {
        mat_err(ctx, 2);
        goto MLS1Fin;
    }
    if ((p = n_mexpr(ctx, p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
        goto MLS1Fin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MLS1Fin;
    if (row1 != 1 || col1 != 1) {
        mat_err(ctx, 8);
        goto MLS1Fin;
    }
    me = (int)ctx->MX[1][1];

    if ((p = n_mexpr(ctx, p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
        goto MLS1Fin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MLS1Fin;
    if (row1 != 1 || col1 != 1) {
        mat_err(ctx, 8);
        goto MLS1Fin;
    }
    mi = (int)ctx->MX[1][1];

    ma = row - me - mi;
    if (me < 0 || mi < 0 || ma < 1) {
        mat_err(ctx, 14);
        goto MLS1Fin;
    }
    n = col - 1;

    if (!(xx = (double *)calloc((size_t)(n) * (size_t)(n) + 1,sizeof(double)))) {
        mat_err(ctx, 3);
        goto MLS1Fin;
    }                
    xxa = n * n + 1;
    memrq(ctx, xxa,sizeof(double));

    if (!(xy = (double *)calloc((size_t)(n + 1),sizeof(double)))) {
        mat_err(ctx, 3);
        goto MLS1Fin;
    }                
    xya = n + 1;
    memrq(ctx, xya,sizeof(double));

    m = me + mi;

    if (m > 0) {
        if (!(a = (double *)calloc((size_t)(m) * (size_t)(n) + 1,sizeof(double)))) {
            mat_err(ctx, 3);
            goto MLS1Fin;
        }                
        aa = m * n + 1;
        memrq(ctx, aa,sizeof(double));

        if (!(b = (double *)calloc((size_t)(m + 1),sizeof(double)))) {
            mat_err(ctx, 3);
            goto MLS1Fin;
        }                
        ba = m + 1;
        memrq(ctx, ba,sizeof(double));
    }

    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) {
            tmp = 0.0;
            for (k = 1; k <= ma; ++k) {
                tmp += ctx->MX[0][(me + k - 1) * col + i] *
                       ctx->MX[0][(me + k - 1) * col + j];
            }
            xx[(i - 1) * n + j] = tmp;
        }   
        tmp = 0.0;
        for (k = 1; k <= ma; ++k) 
            tmp -= ctx->MX[0][(me + k - 1) * col + i] *
                       ctx->MX[0][(me + k - 1) * col + col];
        xy[i] = tmp;
    }
    if (me > 0) {
        for (i = 1; i <= me; ++i) {
            for (j = 1; j <= n; ++j)
                a[(i - 1) * n + j] = ctx->MX[0][(i - 1) * col + j];

            b[i] = ctx->MX[0][(i - 1) * col + col];
        }
    }
    if (mi > 0) {
        for (i = 1; i <= mi; ++i) {
            for (j = 1; j <= n; ++j)
                a[(me + i - 1) * n + j] = ctx->MX[0][(me + ma + i - 1) * col + j];

            b[me + i] = ctx->MX[0][(me + ma + i - 1) * col + col];
        }
    }
          
    if ((p = m_getmat(ctx, p,n,1,&idx,cmd,1)) == NULL)
        goto MLS1Fin;

    maxit = 50 * (n + m);
    vsmall = 2.0 * ctx->EPSI;
        
    /* allocate and initialize xl and xu for lower and upper bounds */

    if (!(xl = (double *)calloc((size_t)(row + 1),sizeof(double)))) {
        mat_err(ctx, 3);
        goto MLS1Fin;
    }                
    xla = row + 1;
    memrq(ctx, xla,sizeof(double));

    if (!(xu = (double *)calloc((size_t)(row + 1),sizeof(double)))) {
        mat_err(ctx, 3);
        goto MLS1Fin;
    }                
    xua = row + 1;
    memrq(ctx, xua,sizeof(double));

    for (i = 1; i <= row; ++i) {
        xu[i] = ctx->DBLMAX / 1000.0;
        xl[i] = -xu[i];             
    }

    /* allocate iact */

    if (!(iact = (int *)calloc((size_t)(row + 1),sizeof(int)))) {
        mat_err(ctx, 3);
        goto MLS1Fin;
    }                
    iaa = row + 1;
    memrq(ctx, iaa,sizeof(int));

    r = qld(ctx, n,m,me,1,maxit,a,b,xy,xx,ctx->MatVal[idx],xl,xu,iact,&nact,vsmall,&diag);

    /* tda_out("r=%d nact=%d diag=%g\n",r,nact,diag); */

    switch (r) {
        case  0: break;
        case  1: m_cmdmsg(ctx);
                 printf1(ctx, "Exceeded maximal number of iterations (%d).\n",maxit);
                 goto MLS1Fin;
        case  2: m_cmdmsg(ctx);
                 printf1(ctx, "Insufficient accuracy for convergence.\n");
                 goto MLS1Fin;
        case  3: m_cmdmsg(ctx);
                 printf1(ctx, "Insufficient memory.\n");
                 goto MLS1Fin;
        default: m_cmdmsg(ctx);
                 printf1(ctx, "Inconsistent constraints.\n");
                 goto MLS1Fin;
    }
    if (diag > 0.0)   
        printf2(ctx, "Coeff matrix was enlarged by %g times unitmatrix.\n",diag);
    err = 0;

MLS1Fin:
    if (xla > 0) {
        free((char *)xl);
        memrq(ctx, -xla,sizeof(double));
    }
    if (xua > 0) {
        free((char *)xu);
        memrq(ctx, -xua,sizeof(double));
    }
    if (iaa > 0) {
        free((char *)iact);
        memrq(ctx, -iaa,sizeof(int));
    }
    if (xxa > 0) {
        free((char *)xx);
        memrq(ctx, -xxa,sizeof(double));
    }
    if (xya > 0) {
        free((char *)xy);
        memrq(ctx, -xya,sizeof(double));
    }
    if (aa > 0) {
        free((char *)a);
        memrq(ctx, -aa,sizeof(double));
    }
    if (ba > 0) {
        free((char *)b);
        memrq(ctx, -ba,sizeof(double));
    }
    mx_free(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_kmet(cmd)      mkmet(A,D)     Kemeny metric                           */
/*                                                                          */
/*  Input:  (n,m) matrix A, each row contains a rank ordering.              */
/*  Output: (n,n) matrix D, with d(i,j) the Kemeny distance between row i   */
/*          and row j of A. (D symmetric with zero diagonal elements).      */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_kmet(TDAContext *ctx, char *cmd)
{
    register int i,j;
    int err,row,col,idx,ivflg,a,b;
    register char *p;
    double ai,aj,bi,bj,d;

    err = -1;
    ivflg = 0;

    printf2(ctx, "%s\n",cmd);
    p = n_mexpr(ctx, cmd + 6,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MKMETFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MKMETFin;
     
    if ((p = m_getmat(ctx, p,row,row,&idx,cmd,1)) == NULL)
        goto MKMETFin;

    for (a = 2; a <= row; ++a) {
        for (b = 1; b < a; ++b) {

            /* calculate Kemeny distance between row a and row b */

            d = 0.0;
            for (i = 2; i <= col; ++i) {
                ai = get_mxval(ctx, 0,a,i);  
                bi = get_mxval(ctx, 0,b,i);  
                for (j = 1; j < i; ++j) {
                    aj = get_mxval(ctx, 0,a,j);  
                    bj = get_mxval(ctx, 0,b,j);  
                    if (fabs(ai - aj) <= ctx->EPSI2) {
                        if (fabs(bi - bj) > ctx->EPSI2)  
                            d += 1.0;
                    }
                    else if (ai < aj) {
                        if (fabs(bi - bj) <= ctx->EPSI2)  
                            d += 1.0;
                        else if (bi > bj)
                            d += 2.0;
                    }
                    else {
                        if (fabs(bi - bj) <= ctx->EPSI2)  
                            d += 1.0;
                        else if (bi < bj)
                            d += 2.0;
                    }
                }
            }
            ctx->MatVal[idx][(a - 1) * row + b] = d;            
            ctx->MatVal[idx][(b - 1) * row + a] = d;            
        }
    }
    err = 0;

MKMETFin:
    mx_free(ctx);
    return(err);
}



