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

/* ------------------------------------------------------------------------ */
/*  functions in t_matc                                                     */

char *n_mexpr(char *p,int n,int idx,int *rmax,int *cmax,int *ivflg,char *sel);
int sel_expr(char *sel,int row,int col,int *rsel,int *csel,int *idx,double *dval);
void mx_free(void);
int mx_alloc(int idx,int row,int col);
int mx1_alloc(int idx,int row,int col);
double get_mxval(int idx,int row,int col);
double get_matval(int idx,int row,int col);
void m_copy(double *x,double *y,int n);
char *m_check1(char *p);                  
char *m_check2(char *p,int opt);                  
char *m_getmat(char *p,int row,int col,int *idx,char *cmd,int opt);
int m_cent(char *cmd,int opt);
int m_dcent(char *cmd);
int m_cross(char *cmd);
int m_chol(char *cmd);
int m_agg(char *cmd); 
int m_vec(char *cmd,int opt);
int m_ivec(char *cmd);
int m_mldes(char *cmd);
int m_sum(char *cmd,int opt);
int m_drow(char *cmd,int opt);
int m_diag(char *cmd,int opt);
int m_transp(char *cmd);
int m_sqrt(char *cmd,int opt);
int m_nrow(char *cmd,int opt);
int m_num(char *cmd);
int m_mch(char *cmd);
int m_mpfit(char *cmd);
int m_mpinv(char *cmd);
int m_mperm(char *cmd,int opt);
int m_mqap(char *cmd);
int m_mcel(char *cmd);
int m_mnc(char *cmd);
int m_mpz(char *cmd);
int m_mpb(char *cmd,int opt);
int m_mevs(char *cmd);
int m_mev(char *cmd);
int m_ginv(char *cmd);
int m_svd(char *cmd,int opt);
int m_wvec(char *cmd);
int m_wvec1(char *cmd);
int m_scal1(char *cmd);
int m_mdefc(char *cmd);
int m_mdefi(char *cmd);
int m_invs(char *cmd);
int m_invd(char *cmd);
int m_ple(char *cmd);
int m_nvar(char *cmd);
int m_print(char *cmd); 
int m_mkp(char *cmd);
int m_mls(char *cmd,int opt);
int m_mlp(char *cmd,int opt);
int m_mlpi(char *cmd);
int m_mul(char *cmd);
int m_cat(char *cmd,int opt);
int m_srow(char *cmd,int opt);
char *m_selidx(char *s,int opt,int *n,int *idx);
int m_sort(char *cmd,int opt);
char *get_mselect(char *p);
int m_exp(char *cmd);
int m_expr(char *cmd);
int m_expr1(char *cmd);
int m_setv(char *cmd); 
int m_brr(char *cmd); 
int m_mtrim(char *cmd); 
int m_mmp(char *cmd,int opt);
int m_mpit(char *cmd,int opt);
int m_mqp(char *cmd,int opt);
int m_mlsei1(char *cmd);
int m_kmet(char *cmd);

/* ------------------------------------------------------------------------ */
/*  gobal variables                                                         */

#define MXN 7

double *MX[MXN];
double *MX1[MXN];
int MXRow[MXN];
int MXCol[MXN];
int MXAlloc[MXN];
int MX1Alloc[MXN];
int MXIV[MXN];
#define MXMaxLen 60
char MXName[MXMaxLen + 1];      /* name of m-expression */

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

char *n_mexpr(char *p,int n,int idx,int *rmax,int *cmax,int *ivflg,char *sel)
{
    register int ii,i;
    int err,row,col,m,iv1flg,iv2flg,rsela,csela,sidx;
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
            mat_err(0);
            break;
        }
        q = skip_expr(p);
        c = *q;
        *q = '\0';

        if (get_mexpr(p,&row,&col,&iv1flg) || row < 1 || col < 1)
             goto N_MEXPRFin;

        MXIV[idx + ii] = 0;

        strncpy(MXName,p,MXMaxLen);
        *(MXName + MXMaxLen) = '\0';

        *rmax = imax(*rmax,row);
        *cmax = imax(*cmax,col);

        rsel = csel = NULL;
 
        if (sel != NULL) {              /* get selection strings */
            if (!(rsel = (int *)calloc(row  + 1,sizeof(int)))) {
                p_err(-2,1);
                goto N_MEXPRFin;
            }
            rsela = row + 1;
            memrq(rsela,sizeof(int));

            if (!(csel = (int *)calloc(col  + 1,sizeof(int)))) {
                p_err(-2,1);
                goto N_MEXPRFin;
            }
            csela = col + 1;
            memrq(csela,sizeof(int));

            if (sel_expr(sel,row,col,rsel,csel,&sidx,&dval))
                goto N_MEXPRFin;
        }

        m = iabs(ESTyp[0]);
        if (ESCnt == 1 && m >= MOFFS && m < VOFFS) {
            m -= MOFFS;
        }
        else {
            m = -1;
            if (alloc_actmp(row * col + 1))
                goto N_MEXPRFin;

            if (iv1flg) {
                if (alloc_actmp1(row * col + 1))
                    goto N_MEXPRFin;
            }
            if (eval_mexpr(0,row,col,AcTmp,iv1flg,AcTmp1,rsel,csel,sidx,dval))
                goto N_MEXPRFin;
        }
        if (mx_alloc(idx + ii,row,col))
            goto N_MEXPRFin;

        if (m >= 0) {
            m_copy(MX[idx + ii],MatVal[m],row * col);
            if (ESTyp[0] < 0) {
                for (i = 1; i <= row * col; ++i)
                    MX[idx + ii][i] *= -1;
            }
        }
        else {
            m_copy(MX[idx + ii],AcTmp,row * col);

            if (iv1flg) {
                if (mx1_alloc(idx + ii,row,col))
                    goto N_MEXPRFin;

                m_copy(MX1[idx + ii],AcTmp1,row * col);
                MXIV[idx + ii] = iv2flg = 1;
            }
        }

        *q = c;
        p = q;  
        if (ii == n - 1)    
            break;
        if (*p++ != ',' || !*p) {
            mat_err(0);
            break;
        }
    }
    *ivflg = iv2flg;
    err = 0;

N_MEXPRFin:
    if (rsela > 0) {
        free((char *)rsel);
        memrq(-rsela,sizeof(int));
    }
    if (csela > 0) {
        free((char *)csel);
        memrq(-csela,sizeof(int));
    }
    alloc_actmp(0);
    alloc_actmp1(0);
    if (err) {
        mx_free();
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

int sel_expr(char *sel,int row,int col,int *rsel,int *csel,int *idx,double *dval)
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
                p = skip_int(p) + 1;
                p = skip_int(p) + 1;
                p = skip_int(p);
            }
            else {
                while (*p) {
                    if (sscanf(p,"%d",&a) != 1)  
                        goto SELEFin;

                    if (a >= 1 && a <= n)  
                        nsel[a - 1] = 1;
                    p = skip_int(p);
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
        if ((*idx = mat_getidx(p,1)) < 0)    /* check for matrix */
            goto SELEFin;
    }
    return(0);

SELEFin:
    printf1("Error in selection string.\n");
    return(-1);
}

/*--------------------------------------------------------------------------*/
/*  mx_free()       Free MX[].                                              */

void mx_free(void)
{
    register int i;

    for (i = 0; i < MXN; ++i) {
        if (MXAlloc[i] > 0) {
            free((char *)MX[i]);
            memrq(-MXAlloc[i],sizeof(double));
            MXAlloc[i] = 0;
        }
        if (MX1Alloc[i] > 0) {
            free((char *)MX1[i]);
            memrq(-MX1Alloc[i],sizeof(double));
            MX1Alloc[i] = 0;
        }
    }
}

/*--------------------------------------------------------------------------*/
/*  mx_alloc(idx,row,col)   Allocate MX[idx]                                */
/*                          Return 0 if OK, -1 if error.                    */

int mx_alloc(int idx,int row,int col)
{

    if (idx < 0 || idx >= MXN || row < 1 || col < 1) {
        fprintf(stderr,"ERROR in mx_alloc: idx=%d row=%d col=%d\n",idx,row,col);
        exit(0);
    }
    if (MXAlloc[idx] > 0) {
        free((char *)MX[idx]);
        memrq(-MXAlloc[idx],sizeof(double));
        MXAlloc[idx] = 0;
    }
    if (!(MX[idx] = (double *)calloc(row * col + 1,sizeof(double)))) {
        m_cmdmsg();
        p_err(-2,1);
        return(-1);
    }         
    memrq(row * col + 1,sizeof(double));
    MXAlloc[idx] = row * col + 1;
    MXRow[idx] = row;
    MXCol[idx] = col;
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  mx1_alloc(idx,row,col)   Allocate MX1[idx]                              */
/*                          Return 0 if OK, -1 if error.                    */

int mx1_alloc(int idx,int row,int col)
{

    if (idx < 0 || idx >= MXN || row < 1 || col < 1) {
        fprintf(stderr,"ERROR in mx1_alloc: idx=%d row=%d col=%d\n",idx,row,col);
        exit(0);
    }
    if (MX1Alloc[idx] > 0) {
        free((char *)MX1[idx]);
        memrq(-MX1Alloc[idx],sizeof(double));
        MX1Alloc[idx] = 0;
    }
    if (!(MX1[idx] = (double *)calloc(row * col + 1,sizeof(double)))) {
        m_cmdmsg();
        p_err(-2,1);
        return(-1);
    }         
    memrq(row * col + 1,sizeof(double));
    MX1Alloc[idx] = row * col + 1;
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  get_mxval(idx,row,col)  Return MX[idx][row,col]                         */
/*                          row, col runs from 1,2,3,...                    */

double get_mxval(int idx,int row,int col)
{
    register int i,j;

    i = (--row) % MXRow[idx];
    j = (--col) % MXCol[idx] + 1;
    return(MX[idx][i * MXCol[idx] + j]);
}

/*--------------------------------------------------------------------------*/
/*  get_matval(idx,row,col) Return MatVal[idx][row,col]                     */
/*                          row, col runs from 1,2,3,...                    */

double get_matval(int idx,int row,int col)
{
    register int i,j;

    i = (--row) % MatRow[idx];
    j = (--col) % MatCol[idx] + 1;
    return(MatVal[idx][i * MatCol[idx] + j]);
}

/*--------------------------------------------------------------------------*/
/*  m_copy(x,y,n)       copy n elements from y to x.                        */

void m_copy(double *x,double *y,int n)
{
    register int i;
    for (i = 1; i <= n; ++i)
        x[i] = y[i];
}

/*--------------------------------------------------------------------------*/
/*  m_check1(p)     check whether *p is a comma. If it is return p++,       */
/*                  otherwise return NULL                                   */

char *m_check1(char *p)                   
{
    if (*p++ != ',' || !*p) {
        mat_err(0);
        return(NULL);
    }
    return(p);
}

/*--------------------------------------------------------------------------*/
/*  m_check2(p,opt) check whether *p is a ")". If it is return p++,         */
/*                  otherwise return NULL                                   */
/*                  If opt != 0, also check that "0" is last character      */
/*                  in string.                                              */

char *m_check2(char *p,int opt)           
{
    if (*p++ != ')' || (opt && *p)) {
        mat_err(0);
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

char *m_getmat(char *p,int row,int col,int *idx,char *cmd,int opt)
{
    char mname[VNLMax + 1];  

    if ((p = get_mname(p,mname,1)) == NULL)
        return(p);
    if ((*idx = mat_newmat(mname,row,col)) < 0)
        return(NULL);
    mdefcpy(MatDef[*idx],cmd);
    if (opt)
        p = m_check2(p,1);
    return(p);
}

/*--------------------------------------------------------------------------*/
/*  m_cent(cmd,opt)  opt=0:  mcent(A,R)                                     */
/*                   opt=1:  mstand(A,R)                                    */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_cent(char *cmd,int opt)
{
    register int i,j;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp,tmp1;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = cmd + 6;
    if (opt)
        p++;

    p = n_mexpr(p,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MCENTFin;
    if ((p = m_check1(p)) == NULL)
        goto MCENTFin;
    if ((p = m_getmat(p,row,col,&idx,cmd,1)) == NULL)
        goto MCENTFin;

    for (j = 1; j <= col; ++j) {
        tmp = 0.0;
        for (i = 1; i <= row; ++i)  
            tmp += get_mxval(0,i,j);

        tmp /= (double)row;
        for (i = 0; i < row; ++i)  
            MatVal[idx][i * col + j] = get_mxval(0,i + 1,j) - tmp;

        if (opt) {                      /* mstand */
            tmp = 0.0;
            for (i = 0; i < row; ++i) {
                tmp1 = MatVal[idx][i * col + j];
                tmp += tmp1 * tmp1;
            }
            tmp /= (double)row;

            if (tmp > 0.0) {
                tmp = sqrt(tmp);
                for (i = 0; i < row; ++i)    
                    MatVal[idx][i * col + j] /= tmp;
            }
        }
    }
    err = 0;

MCENTFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_dcent(cmd)       mdcent(R,A)                                          */
/*                     R should be a (symmetric) n x n matrix.              */
/*                     The command creates a double-centered matrix A       */
/*                     as used with classical metric scaling.               */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_dcent(char *cmd)
{
    register int i,j,k;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp,aij,r00,ri0,r0j;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = cmd + 7;

    p = n_mexpr(p,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MDCENTFin;
    if ((p = m_check1(p)) == NULL)
        goto MDCENTFin;
    if (row != col) {
        mat_err(1);
        goto MDCENTFin;
    }
    if ((p = m_getmat(p,row,col,&idx,cmd,1)) == NULL)
        goto MDCENTFin;

    r00 = 0.0;
    for (i = 1; i <= row; ++i) {
        for (j = 1; j <= col; ++j) {
            tmp = get_mxval(0,i,j);
            r00 += tmp * tmp;
        }
    }
    r00 /= (double)(row * row);

    for (i = 1; i <= row; ++i) {
        ri0 = 0.0;
        for (j = 1; j <= col; ++j) {
            tmp = get_mxval(0,i,j);
            ri0 += tmp * tmp;
        }
        ri0 /= (double)col;

        for (j = 1; j <= col; ++j) {
            tmp = get_mxval(0,i,j);
            aij = r00 + tmp * tmp - ri0;

            r0j = 0.0;
            for (k = 1; k <= row; ++k) {
                tmp = get_mxval(0,k,j);
                r0j += tmp * tmp;
            }
            r0j /= (double)row;
            aij -= r0j;
            MatVal[idx][(i - 1) * col + j] = -aij / 2.0;
        }
    }
    err = 0;

MDCENTFin:
    mx_free();
    return(err);
}
     
/*--------------------------------------------------------------------------*/
/*  m_cross(cmd)    mcross(A,R)          R = A'A                            */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_cross(char *cmd)
{
    register int i,j,l;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 7,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MCROSSFin;
    if ((p = m_check1(p)) == NULL)
        goto MCROSSFin;
    if ((p = m_getmat(p,col,col,&idx,cmd,1)) == NULL)
        goto MCROSSFin;

    for (i = 0; i < col; ++i) {
        for (j = 1; j <= col; ++j) {
            tmp = 0.0;
            for (l = 1; l <= row; ++l)  
                tmp += get_mxval(0,l,i + 1) * get_mxval(0,l,j);
            MatVal[idx][i * col + j] = tmp;
        }
    }
    err = 0;

MCROSSFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_chol(cmd)      mchol(A,R)     Cholesky decomposition                  */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_chol(char *cmd)
{
    register int i,j,k;
    int err,row,col,idx,n,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 6,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MCHOLFin;
    if ((p = m_check1(p)) == NULL)
        goto MCHOLFin;
    if (row != col) {
        mat_err(1);
        goto MCHOLFin;
    }
    if ((p = m_getmat(p,row,col,&idx,cmd,1)) == NULL)
        goto MCHOLFin;

    n = (row * (row + 1)) / 2;
    if (alloc_actmp(n + 1))
        goto MCHOLFin;

    k = 1;
    for (i = 1; i <= row; ++i) {
        for (j = 1; j <= i; ++j)
            AcTmp[k++] = get_mxval(0,i,j);  
    }
    if (decomp(row,AcTmp)) {
        m_cmdmsg();
        printf1("Error: %s is not positive definite.\n",MXName);
        goto MCHOLFin;     
    }
    k = 1;
    for (i = 0; i < row; ++i) {
        for (j = 0; j <= i; ++j)
            MatVal[idx][i * row + j + 1] = AcTmp[k++];
    }
    err = 0;

MCHOLFin:
    alloc_actmp(0);
    mx_free();
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

int m_agg(char *cmd)
{
    register int i,j,ii,jj,k;
    int err,row,row1,row2,col,col1,col2,idx,ivflg,nr,nc,rmax,cmax,ir,ic;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    if ((p = n_mexpr(cmd + 4,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MAGFin;
    if ((p = m_check1(p)) == NULL)
        goto MAGFin;
    if ((p = n_mexpr(p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
        goto MAGFin;
    if (row1 != row || col1 != 1) {
        mat_err(2);
        goto MAGFin;
    }
    if ((p = m_check1(p)) == NULL)
        goto MAGFin;
    if ((p = n_mexpr(p,1,2,&row2,&col2,&ivflg,NULL)) == NULL)
        goto MAGFin;
    if ((p = m_check1(p)) == NULL)
        goto MAGFin;
    if (row2 != 1 || col2 != col) {
        mat_err(2);
        goto MAGFin;
    }
    rmax = 0;
    for (i = 1; i <= row; ++i)  
        rmax = imax(rmax,(int)MatVal[1][i]);
    cmax = 0;
    for (j = 1; j <= col; ++j)  
        cmax = imax(cmax,(int)MatVal[2][j]);

    if (rmax < 1 || cmax < 1) {
        mat_err(15);
        goto MAGFin;
    }
    if (alloc_acr(rmax + 1))
        goto MAGFin;
    if (alloc_acs(cmax + 1))
        goto MAGFin;
      
    for (i = 1; i <= row; ++i) {
        if ((k = (int)MatVal[1][i]) > 0)  
            AcR[k] = 1;
    }
    for (j = 1; j <= col; ++j) {
        if ((k = (int)MatVal[2][j]) > 0)  
            AcS[k] = 1;
    }
    nr = nc = 0;
    for (i = 1; i <= rmax; ++i) {
        if (AcR[i] > 0)
            nr++;
    }
    for (j = 1; j <= cmax; ++j) {
        if (AcS[j] > 0)
            nc++;
    }

    if ((p = m_getmat(p,nr,nc,&idx,cmd,1)) == NULL)
        goto MAGFin;
      
    ir = 0;
    for (ii = 1; ii <= rmax; ++ii) {
        if (AcR[ii] == 0)
            continue;
        ir++;
        k = (ir - 1) * nc;
        ic = 0;

        for (jj = 1; jj <= cmax; ++jj) {
            if (AcS[jj] == 0)
                continue;
            ic++;

            tmp = 0.0;
            for (i = 1; i <= row; ++i) {
                if (MatVal[1][i] == ii) {
                    for (j = 1; j <= col; ++j) {
                        if (MatVal[2][j] == jj)  
                            tmp += get_mxval(0,i,j);
                    }
                }
            }
            MatVal[idx][k + ic] = tmp;
        }
    }
    err = 0;

MAGFin:
    alloc_acr(0);
    alloc_acs(0);
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_vec(cmd,opt)      opt 0 :  mcvec(A,V)     stack columns               */
/*                      opt 1 :  mrvec(A,V)     stack rows                  */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_vec(char *cmd,int opt)
{
    register int i,j,k;
    int err,row,col,idx,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 6,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MVECFin;
    if ((p = m_check1(p)) == NULL)
        goto MVECFin;

    if ((p = m_getmat(p,row * col,1,&idx,cmd,1)) == NULL)
        goto MVECFin;

    k = 1;
    if (opt == 0) {
        for (j = 1; j <= col; ++j) {
            for (i = 1; i <= row; ++i)     
                MatVal[idx][k++] = get_mxval(0,i,j);
        }
    }
    else {
        for (i = 1; i <= row; ++i) {
            for (j = 1; j <= col; ++j)     
                MatVal[idx][k++] = get_mxval(0,i,j);
        }
    }
    err = 0;

MVECFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_ivec(cmd)     mivec(V,n,A)   inverse vec operator                     */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_ivec(char *cmd)
{
    register int i,j,k;
    int err,row,col,n,m,q,idx,ivflg;
    register char *p;

    ivflg = 0;
    err = -1;
    printf2("%s\n",cmd);

    if ((p = n_mexpr(cmd + 6,1,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto MIVECFin;

    if (col != 1) {
        mat_err(7);
        goto MIVECFin;
    }
    q = row;

    if ((p = m_check1(p)) == NULL)
        goto MIVECFin;
    if ((p = n_mexpr(p,1,1,&row,&col,&ivflg,NULL)) == NULL) 
        goto MIVECFin;

    if (row != 1 || col != 1 || (n = (int)MX[1][1]) < 1) {
        mat_err(18);
        goto MIVECFin;
    }
    m = q / n;
    if (n * m != q) {
        mat_err(18);
        goto MIVECFin;
    }
    if ((p = m_check1(p)) == NULL)
        goto MIVECFin;
    if ((p = m_getmat(p,n,m,&idx,cmd,0)) == NULL)
        goto MIVECFin;

    k = 1;
    for (j = 1; j <= m; ++j) {
        for (i = 0; i < n; ++i)  
            MatVal[idx][i * m + j] = MX[0][k++];
    }   
    err = 0;

MIVECFin:
    mx_free();
    return(err);
}

/*--##----------------------------------------------------------------------*/
/*  m_mldes(cmd)    mldes(Z,G,D)   Design matrix for MLRC models            */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mldes(char *cmd)
{
    register int i1,i2,i,j,k,l,ii,kk,ll;
    int err,row,col,n,m,q,nn,idx,ivflg;
    register char *p;
    double *d,*zi,*zii;

    ivflg = 0;
    err = -1;
    printf2("%s\n",cmd);

    if ((p = n_mexpr(cmd + 6,1,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto MLDESFin;

    n = row;
    q = col;
    m = q * q + 1;
    nn = n * n;

    if ((p = m_check1(p)) == NULL)
        goto MLDESFin;
    if ((p = n_mexpr(p,1,1,&row,&col,&ivflg,NULL)) == NULL) 
        goto MLDESFin;

    if (row != n || col != 1) {
        mat_err(2);
        goto MLDESFin;
    }
    if ((p = m_check1(p)) == NULL)
        goto MLDESFin;
    if ((p = m_getmat(p,nn,m,&idx,cmd,0)) == NULL)
        goto MLDESFin;
       
    if (alloc_aci(n + 1))
        goto MLDESFin;
    if (alloc_acj(n + 1))
        goto MLDESFin;

    AcI[1] = AcJ[1] = 1;

    i = j = 1;
    for (k = 2; k <= n; ++k) {
        i++;
        if ((int)MX[1][k] != (int)MX[1][k - 1]) {
            i = 1;
            j++;
        }
        AcI[k] = i;
        AcJ[k] = j; 
    }
    l = 0;
    for (kk = 1; kk <= n; ++kk) {
        ii = AcI[kk];
        j = AcJ[kk];
        for (k = 1; k <= n; ++k) {
            i = AcI[k];

            /** printf("\ni=%2d j=%2d ii=%2d jj=%2d \n",i,AcJ[k],ii,j); **/

            if (AcJ[k] == j) {
                zi  = MX[0] + (k - 1) * q;
                zii = MX[0] + (kk - 1) * q;

                /** printf("zi=%lf %lf  zii=%lf %lf\n",zi[1],zi[2],zii[1],zii[2]); **/

                ll = 1;
                d = MatVal[idx] + l * m;
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
    mx_free();
    alloc_aci(0);
    alloc_acj(0);
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

int m_sum(char *cmd,int opt)
{
    register int i,j;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 6,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MSUMFin;
    if ((p = m_check1(p)) == NULL)
        goto MSUMFin;

    if (opt == 0) {
        if ((p = m_getmat(p,row,1,&idx,cmd,1)) == NULL)
            goto MSUMFin;

        for (i = 1; i <= row; ++i) {
            tmp = 0.0;
            for (j = 1; j <= col; ++j) 
                tmp += get_mxval(0,i,j);
            MatVal[idx][i] = tmp;
        }
    }
    else {
        if ((p = m_getmat(p,1,col,&idx,cmd,1)) == NULL)
            goto MSUMFin;

        for (j = 1; j <= col; ++j) {
            tmp = 0.0;
            for (i = 1; i <= row; ++i) 
                tmp += get_mxval(0,i,j);
            MatVal[idx][j] = tmp;
        }
    }
    err = 0;

MSUMFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_drow(cmd,opt)     opt 0 :  mdrow(A,R)                                 */
/*                      opt 1 :  mdcol(A,R)                                 */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_drow(char *cmd,int opt)
{
    register int i,j;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 6,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MDROWFin;
    if ((p = m_check1(p)) == NULL)
        goto MDROWFin;

    if (opt == 0) {
        if ((p = m_getmat(p,row,row,&idx,cmd,1)) == NULL)
            goto MDROWFin;

        for (i = 1; i <= row; ++i) {
            tmp = 0.0;
            for (j = 1; j <= col; ++j) 
                tmp += get_mxval(0,i,j);
            MatVal[idx][(i - 1) * row + i] = tmp;
        }
    }
    else {
        if ((p = m_getmat(p,col,col,&idx,cmd,1)) == NULL)
            goto MDROWFin;

        for (j = 1; j <= col; ++j) {
            tmp = 0.0;
            for (i = 1; i <= row; ++i) 
                tmp += get_mxval(0,i,j);
            MatVal[idx][(j - 1) * col + j] = tmp;
        }
    }
    err = 0;

MDROWFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_diag(cmd,opt)     opt 0 :  mdiag(A,R)                                 */
/*                      opt 1 :  mdiagd(A,R)                                */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_diag(char *cmd,int opt)
{
    register int i,j;
    int err,row,col,idx,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = cmd + 6;
    if (opt)
        p++;
    p = n_mexpr(p,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MDIAGFin;
    if ((p = m_check1(p)) == NULL)
        goto MDIAGFin;

    if (opt == 0) {
        if (col != 1) {
            mat_err(7);
            goto MDIAGFin;
        }
        if ((p = m_getmat(p,row,row,&idx,cmd,1)) == NULL)
            goto MDIAGFin;

        for (i = 1; i <= row; ++i)  
            MatVal[idx][(i - 1) * row + i] = get_mxval(0,i,1);
    }
    else {
        j = imin(row,col);
        if ((p = m_getmat(p,j,1,&idx,cmd,1)) == NULL)
            goto MDIAGFin;

        for (i = 1; i <= j; ++i)  
            MatVal[idx][i] = get_mxval(0,i,i);
    }
    err = 0;

MDIAGFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_transp(cmd)    mtransp(A,R)    Transposition                          */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_transp(char *cmd)
{
    register int i,j;
    int err,row,col,idx,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 8,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MTRANFin;
    if ((p = m_check1(p)) == NULL)
        goto MTRANFin;

    if ((p = m_getmat(p,col,row,&idx,cmd,1)) == NULL)
        goto MTRANFin;

    for (i = 1; i <= row; ++i) {
        for (j = 1; j <= col; ++j)  
            MatVal[idx][(j - 1) * row + i] = get_mxval(0,i,j);
    }
    err = 0;

MTRANFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_sqrt(cmd,opt)     opt 0 :  msqrtd(A,R)                                */
/*                      opt 1 :  msqrti(A,R)                                */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_sqrt(char *cmd,int opt)
{
    register int i,j;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 7,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MSQRTFin;
    if ((p = m_check1(p)) == NULL)
        goto MSQRTFin;

    if (row != col) {
        mat_err(1);
        goto MSQRTFin;
    }
    if ((p = m_getmat(p,row,col,&idx,cmd,1)) == NULL)
        goto MSQRTFin;

    j = 1;
    for (i = 1; i <= row; ++i) {
        tmp = get_mxval(0,i,i); 
        if (tmp < 0.0) {
            m_cmdmsg();
            printf1("Error: negative diagonal element in %s.\n",MXName);
            goto MSQRTFin;
        }
        if (tmp > 0.0)
            tmp = sqrt(tmp);
        if (opt == 0) 
            MatVal[idx][j] = tmp;
        else {
            if (tmp <= EPSI) {
                m_cmdmsg();
                printf1("Error: (almost) zero diagonal element in %s.\n",MXName);
                goto MSQRTFin;
            }
            MatVal[idx][j] = 1.0 / tmp;
        }
        j += row + 1;
    }
    err = 0;

MSQRTFin:
    mx_free();
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

int m_nrow(char *cmd,int opt)
{
    register int i,n;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = cmd + 6;
    if (opt >= 3)
        p++;
    p = n_mexpr(p,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MNROWFin;
    if ((p = m_check1(p)) == NULL)
        goto MNROWFin;

    if ((p = m_getmat(p,1,1,&idx,cmd,1)) == NULL)
        goto MNROWFin;

    switch (opt) {
        case  1:    tmp = (double)col;
                    break;
            
        case  2:    n = row * col;
                    tmp = 0.0;
                    for (i = 1; i <= n; ++i)
                        tmp = dmax(tmp,MX[0][i]);
                    break;

        case  3:    n = row * col;
                    tmp = 0.0;
                    for (i = 1; i <= n; ++i)
                        tmp += fabs(MX[0][i]);
                    break;

        case  4:    n = row * col;
                    tmp = 0.0;
                    for (i = 1; i <= n; ++i)
                        tmp += MX[0][i] * MX[0][i];
                    if (tmp > 0.0)
                        tmp = sqrt(tmp);
                    break;

        case  5:    n = imin(row,col);
                    tmp = 0.0;
                    for (i = 0; i < n; ++i)
                        tmp += MX[0][i * col + i + 1];
                    break;

        default:    tmp = (double)row;
                    break;
    }
    MatVal[idx][1] = tmp;
    err = 0;

MNROWFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_num(cmd)      mnum(x,d,n,R)                                           */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_num(char *cmd)
{
    register int i;
    int err,row,col,idx,n,ivflg;
    register char *p;
    double x,d;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 5,3,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MNUMFin;
    if ((p = m_check1(p)) == NULL)
        goto MNUMFin;
    if (row != 1 || col != 1) {
        mat_err(8);
        goto MNUMFin;
    }
    n = (int)get_mxval(2,1,1);
    if (n < 1) {
        mat_err(2);
        goto MNUMFin;
    }
    if ((p = m_getmat(p,n,1,&idx,cmd,1)) == NULL)
        goto MNUMFin;

    x = get_mxval(0,1,1);
    d = get_mxval(1,1,1);

    for (i = 1; i <= n; ++i)
        MatVal[idx][i] = x + d * (double)(i - 1);
    err = 0;

MNUMFin:
    mx_free();
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

int m_mch(char *cmd)
{
    register int i,j;
    int err,row,col,idx,ivflg,n,nf;         
    register char *p;
    double z,s;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = cmd + 4;
    if ((p = n_mexpr(p,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MCHFin;
    if ((p = m_check1(p)) == NULL)
        goto MCHFin;

    if (MXCol[0] != MXRow[0]) {       
        mat_err(1);
        goto MCHFin;
    }
    n = MXRow[0];

    if ((p = m_getmat(p,n,1,&idx,cmd,1)) == NULL)
        goto MCHFin;

    n = MXRow[0];
    for (i = 0; i < n; ++i) {
        z = s = 0.0;
        nf = 0;
        for (j = i + 1; j < n; ++j) {
            z += MX[0][i * n + j + 1];
            s += MX[0][j * n + i + 1];
            /** printf("i=%d z=%g s=%g\n",i,z,s); **/
            if (z < s - EPSI1) {
                /* printf2("mch: i=%d z=%g s=%g\n",i + 1,z,s); */
                nf = 1;
                break;
            }
        }
        if (nf)
            MatVal[idx][i + 1] = 1.0;
    }
    for (i = n - 1; i > 0; --i) {
        z = s = 0.0;
        nf = 0;
        for (j = i - 1; j >= 0; --j) {
            z += MX[0][j * n + i + 1];
            s += MX[0][i * n + j + 1];
            /** printf("i=%d z=%g s=%g\n",i,z,s); **/
            if (z < s - EPSI1) {
                /* printf2("mch: i=%d z=%g s=%g\n",i + 1,z,s); */
                nf = 1;
                break;
            }
        }
        if (nf)
            MatVal[idx][i + 1] = 1.0;
    }
    err = 0;

MCHFin:
    mx_free();
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

int m_mpfit(char *cmd)
{
    int err,row,col,idx,ivflg,iter;
    register char *p;
    double eps;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = cmd + 6;
    if ((p = n_mexpr(p,5,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MRASFin;
    if ((p = m_check1(p)) == NULL)
        goto MRASFin;

    if (MXCol[0] != MXCol[1] || MXRow[1] != 1 ||
        MXRow[0] != MXRow[2] || MXCol[2] != 1 ||    
        MXRow[3] != 1        || MXCol[3] != 1 ||    
        MXRow[4] != 1        || MXCol[4] != 1) {    
        mat_err(2);
        goto MRASFin;
    }
    if ((p = m_getmat(p,row,col,&idx,cmd,1)) == NULL)
        goto MRASFin;

    iter = (int)MX[3][1];
    eps  = MX[4][1];
    if (iter < 1) { 
        printf1("mpfit command: error in number of iterations.\n");
        goto MRASFin;
    }
    if (eps < 0.0) { 
        printf1("mpfit command: error in epsilon.\n");
        goto MRASFin;
    }
    if (mpfit(row,col,MX[0],MX[1],MX[2],MatVal[idx],&iter,&eps))
        goto MRASFin;
    printf2("Number of iterations: %d. Final accuracy: %g\n",iter,eps);
    err = 0;

MRASFin:
    mx_free();
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

int m_mpinv(char *cmd)
{
    register int i,k;
    int err,row,col,idx,ivflg,n;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = cmd + 6;
    if ((p = n_mexpr(p,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MPINVFin;
    if ((p = m_check1(p)) == NULL)
        goto MPINVFin;

    if (row != 1 && col != 1) {
        mat_err(2);
        goto MPINVFin;
    }
    if ((p = m_getmat(p,row,col,&idx,cmd,1)) == NULL)
        goto MPINVFin;

    n = imax(row,col);
    for (i = 1; i <= n; ++i) {
        k = (int)MX[0][i];
        if (k < 1 || k > n) {
            mat_err(17);
            goto MPINVFin;
        }
        MatVal[idx][k] = (double)i;
    }
    err = 0;

MPINVFin:
    mx_free();
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

int m_mperm(char *cmd,int opt)
{
    register int i,j,k,k1;
    int err,row,col,row1,col1,idx,ivflg,n;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = cmd + 6;
    if ((p = n_mexpr(p,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MPERMFin;
    if ((p = m_check1(p)) == NULL)
        goto MPERMFin;

    if ((p = n_mexpr(p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
        goto MPERMFin;
    if ((p = m_check1(p)) == NULL)
        goto MPERMFin;

    if (row1 != 1 && col1 != 1) {
        mat_err(2);
        goto MPERMFin;
    }
    if ((p = m_getmat(p,row,col,&idx,cmd,1)) == NULL)
        goto MPERMFin;

    n = imax(row1,col1);
    if (opt == 0) {
        for (i = 1; i <= n; ++i) {
            k = (int)MX[1][i];
            if (k < 1 || k > row) {
                mat_err(17);
                goto MPERMFin;
            }
        }
    }
    if (opt == 1) {
        if (n != row) {  
            mat_err(2);
            goto MPERMFin;
        }
        for (i = 0; i < row; ++i) {
            k = (int)MX[1][i + 1] - 1;
            if (k < 0 || k >= row) {
                mat_err(17);
                goto MPERMFin;
            }
            for (j = 1; j <= col; ++j)  
                MatVal[idx][i * col + j] = MX[0][k * col + j];
        }
    }
    else if (opt == 2) {
        if (n != col) {  
            mat_err(2);
            goto MPERMFin;
        }
        for (j = 1; j <= col; ++j) {
            k = (int)MX[1][j];
            if (k < 1 || k > col) {
                mat_err(17);
                goto MPERMFin;
            }
            for (i = 0; i < row; ++i)  
                MatVal[idx][i * col + j] = MX[0][i * col + k];
        }
    }
    else {
        if (row != col || n != row) {  
            mat_err(2);
            goto MPERMFin;
        }
        for (i = 0; i < row; ++i) {
            k = (int)MX[1][i + 1] - 1;
            for (j = 1; j <= col; ++j) {
                k1 = (int)MX[1][j];
                MatVal[idx][i * col + j] = MX[0][k * col + k1];
            }
        }
    }
    err = 0;

MPERMFin:
    mx_free();
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

int m_mqap(char *cmd)
{
    register int i,j,n;
    int err,row,col,idx,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = cmd + 5;
    if ((p = n_mexpr(p,3,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MQAPFin;
    if ((p = m_check1(p)) == NULL)
        goto MQAPFin;

    if (MXRow[0] != row || MXCol[0] != col ||
        MXRow[1] != row || MXCol[1] != col || 
        MXRow[2] != row || MXCol[2] != col) { 
        mat_err(2);
        goto MQAPFin;
    }
    if ((p = m_getmat(p,row,1,&idx,cmd,1)) == NULL)
        goto MQAPFin;

    if (alloc_acn(3 * row + 1))
        goto MQAPFin;

    if (alloc_actmp(4 * row + 1))
        goto MQAPFin;

    n = 0;
    j = 1;
    for (i = 1; i <= row; ++i) {
        if (MX[0][j] != 0.0) {
            MX[0][j] = 0.0;
            n++;
        }
        if (MX[1][j] != 0.0) {
            MX[1][j] = 0.0;
            n++;
        }
        j += row + 1;
    }
    if (n > 0)  
        printf1("Warning: changed %d main diagonal elements to zero.\n",n);

    if (qap_w(row,MX[2],MX[0],MX[1],AcN,AcTmp)) 
        goto MQAPFin;

    printf1("Best value: ",AcTmp[1]);
    printf1(PMATFmtS,AcTmp[1]);
    printf1("\n");
    for (i = 1; i <= row; ++i)
        MatVal[idx][i] = (double)AcN[i];

    err = 0;

MQAPFin:
    mx_free();
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

int m_mcel(char *cmd)
{
    register int i,j,k;
    int err,row,col,idx,ivflg,m,n;       
    register char *p;
    double x,aij;         

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = cmd + 5;
    if ((p = n_mexpr(p,2,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MCELFin;
    if ((p = m_check1(p)) == NULL)
        goto MCELFin;

    if (MXRow[0] != MXCol[0] || MXRow[1] != 1 || MXCol[1] != 1) { 
        mat_err(2);
        goto MCELFin;
    }
    row = MXRow[0];

    x = MX[1][1];
    m = 0;
    n = row * row;
    for (i = 1; i <= n; ++i) {
        if (MX[0][i] >= x)
            m++;
    }
    if (m == 0) {
        printf1("mcel: number of edges is zero.\n");
        goto MCELFin;
    }
    if ((p = m_getmat(p,m,3,&idx,cmd,1)) == NULL)
        goto MCELFin;

    k = 0;
    for (i = 0; i < row; ++i) {
        for (j = 1; j <= row; ++j) {
            aij = MX[0][i * row + j];
            if (aij >= x) {         
                MatVal[idx][k * 3 + 1] = (double)(i + 1);
                MatVal[idx][k * 3 + 2] = (double)j;
                MatVal[idx][k * 3 + 3] = aij;       
                if (++k >= m)
                    break;
            }
        }
    }
    err = 0;

MCELFin:
    mx_free();
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

int m_mnc(char *cmd)
{
    register int i;
    int err,row,col,idx,ivflg,n,opt;
    register char *p;
    double x;             

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = cmd + 4;
    if ((p = n_mexpr(p,3,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MNCFin;
    if ((p = m_check1(p)) == NULL)
        goto MNCFin;

    if (MXRow[1] != 1 || MXCol[1] != 1 || MXRow[2] != 1 || MXCol[2] != 1) { 
        mat_err(2);
        goto MNCFin;
    }
    row = MXRow[0];
    col = MXCol[0];

    if ((p = m_getmat(p,row,col,&idx,cmd,1)) == NULL)
        goto MNCFin;

    x   = MX[1][1];
    opt = MX[2][1];
    n = row * col;

    switch (opt) {
        case 1:     for (i = 1; i <= n; ++i) {
                        if (MX[0][i] <= x)
                            MatVal[idx][i] = 0.0;
                        else
                            MatVal[idx][i] = MX[0][i];
                    }
                    break;

        case 2:     for (i = 1; i <= n; ++i) {
                        if (MX[0][i] < x)
                            MatVal[idx][i] = 0.0;
                        else
                            MatVal[idx][i] = MX[0][i];
                    }
                    break;

        case 3:     for (i = 1; i <= n; ++i) {
                        if (MX[0][i] >= x)
                            MatVal[idx][i] = 0.0;
                        else
                            MatVal[idx][i] = MX[0][i];
                    }
                    break;

        case 4:     for (i = 1; i <= n; ++i) {
                        if (MX[0][i] > x)
                            MatVal[idx][i] = 0.0;
                        else
                            MatVal[idx][i] = MX[0][i];
                    }
                    break;

        case 5:     for (i = 1; i <= n; ++i) {
                        if (fabs(MX[0][i]) <= x)
                            MatVal[idx][i] = 0.0;
                        else
                            MatVal[idx][i] = MX[0][i];
                    }
                    break;

        case 6:     for (i = 1; i <= n; ++i) {
                        if (fabs(MX[0][i]) < x)
                            MatVal[idx][i] = 0.0;
                        else
                            MatVal[idx][i] = MX[0][i];
                    }
                    break;

        case 7:     for (i = 1; i <= n; ++i) {
                        if (fabs(MX[0][i]) >= x)
                            MatVal[idx][i] = 0.0;
                        else
                            MatVal[idx][i] = MX[0][i];
                    }
                    break;

        case 8:     for (i = 1; i <= n; ++i) {
                        if (fabs(MX[0][i]) > x)
                            MatVal[idx][i] = 0.0;
                        else
                            MatVal[idx][i] = MX[0][i];
                    }
                    break;

        default:    mat_err(0);
                    goto MNCFin;
    }
    err = 0;

MNCFin:
    mx_free();
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

int m_mpz(char *cmd)
{
    register int i,j,k;
    int err,row,col,idx,idx1,ivflg,n,nn,r;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 4,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MPZFin;
    if ((p = m_check1(p)) == NULL)
        goto MPZFin;
    if (row != col) {
        mat_err(1);
        goto MPZFin;
    }
    if ((p = m_getmat(p,row,col,&idx,cmd,0)) == NULL)
        goto MPZFin;
    if ((p = m_check1(p)) == NULL)
        goto MPZFin;
    if ((p = m_getmat(p,row,1,&idx1,cmd,1)) == NULL)
        goto MPZFin;

    nn = 0;
    n = row * col;
    r = row + 1;

    for (i = 1; i <= n; ++i) {
        if (MX[0][i] != 0.0)
            nn++;
    }
    if (alloc_acn(nn + 1))
        goto MPZFin;
    if (alloc_aci(row + 1))
        goto MPZFin;
    if (alloc_acj(row + 1))
        goto MPZFin;
    if (alloc_acm(row + 1))
        goto MPZFin;
    if (alloc_acr(4 * r + 1))
        goto MPZFin;

    k = 0;
    for (i = 0; i < row; ++i) {
        AcI[i + 1] = k + 1;
        n = 0;
        for (j = 1; j <= col; ++j) {
            if (MX[0][i * col + j] != 0.0) {
                AcN[++k] = j;
                n++;
            }
        }
        AcJ[i + 1] = n;
    }
    n = m_rperm(row,AcN,nn,AcI,AcJ,AcM,AcR,AcR + r,AcR + 2 * r,AcR + 3 * r);

    printf2("Number of non-zeros: %d\n",n);

    for (i = 1; i <= row; ++i) {
        k = AcM[i];
        MatVal[idx1][i] = (double)k;          
        for (j = 1; j <= col; ++j)
            MatVal[idx][(k - 1) * col + j] = MX[0][(i - 1) * col + j];
    }
    err = 0;

MPZFin:
    mx_free();
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

int m_mpb(char *cmd,int opt)
{
    register int i,j,k,l;
    int err,row,col,idx,idx1,idx2,idx3,ivflg,n,nn,r;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 5,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MPBFin;
    if ((p = m_check1(p)) == NULL)
        goto MPBFin;
    if (row != col) {
        mat_err(1);
        goto MPBFin;
    }
    if ((p = m_getmat(p,row,col,&idx,cmd,0)) == NULL)
        goto MPBFin;
    if ((p = m_check1(p)) == NULL)
        goto MPBFin;

    if ((p = m_getmat(p,row,1,&idx1,cmd,0)) == NULL)
        goto MPBFin;
    if ((p = m_check1(p)) == NULL)
        goto MPBFin;

    if ((p = m_getmat(p,1,1,&idx2,cmd,0)) == NULL)
        goto MPBFin;
    if ((p = m_check1(p)) == NULL)
        goto MPBFin;

    nn = 0;
    n = row * col;
    r = row + 1;

    for (i = 1; i <= n; ++i) {
        if (MX[0][i] != 0.0)
            nn++;
    }
    if (alloc_acn(nn + 1))
        goto MPBFin;
    if (alloc_aci(row + 1))
        goto MPBFin;
    if (alloc_acj(row + 1))
        goto MPBFin;
    if (alloc_acm(row + 1))
        goto MPBFin;
    if (alloc_acs(row + 1))
        goto MPBFin;
    if (alloc_acr(3 * r + 1))
        goto MPBFin;

    k = 0;
    for (i = 0; i < row; ++i) {
        AcI[i + 1] = k + 1;
        n = 0;
        for (j = 1; j <= col; ++j) {
            if (opt == 0) {
                if (MX[0][i * col + j] != 0.0) {
                    AcN[++k] = j;
                    n++;
                }
            }
            else {
                if (MX[0][(j - 1) * col + i + 1] != 0.0) {
                    AcN[++k] = j;
                    n++;
                }
            }
        }
        AcJ[i + 1] = n;
    }

    n = m_perm(row,AcN,nn,AcI,AcJ,AcM,AcS,AcR,AcR + r,AcR + 2 * r);

    printf2("Number of blocks: %d\n",n);
    
    if ((p = m_getmat(p,n,1,&idx3,cmd,1)) == NULL)
        goto MPBFin;

    MatVal[idx2][1] = (double)n;

    for (i = 1; i <= n; ++i)  
        MatVal[idx3][i] = (double)AcS[i];

    for (i = 1; i <= row; ++i) {
        k = AcM[i];
        MatVal[idx1][i] = (double)k;          

        for (j = 1; j <= col; ++j) {
            l = AcM[j];
            MatVal[idx][(i - 1) * col + j] = MX[0][(k - 1) * col + l];
        }
    }
    err = 0;

MPBFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mevs(cmd)     mevs(A,E,V)                                             */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mevs(char *cmd)
{
    int err,row,col,idx,idx1,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 5,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MEVSFin;
    if ((p = m_check1(p)) == NULL)
        goto MEVSFin;
    if (row != col) {
        mat_err(1);
        goto MEVSFin;
    }
    if ((p = m_getmat(p,row,1,&idx,cmd,0)) == NULL)
        goto MEVSFin;
    if ((p = m_check1(p)) == NULL)
        goto MEVSFin;
    if ((p = m_getmat(p,row,col,&idx1,cmd,1)) == NULL)
        goto MEVSFin;

    m_copy(MatVal[idx1],MX[0],row * col);
    err = evecf(row,MatVal[idx],MatVal[idx1]);
    if (err == -2) {
        mat_err(3);
        goto MEVSFin;
    }
    else if (err) {  
        m_cmdmsg();
        printf1("Error: no success in calculating eigenvalues.\n");
        goto MEVSFin;
    }
    err = 0;

MEVSFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mev(cmd)     mev(A,ER,EI,EVR,EVI)    (Grad/Brebner)                   */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mev(char *cmd)
{
    int err,row,col,idx,idx1,idx2,idx3,ivflg,indica;
    int *indic;
    register char *p;

    err = -1;
    indica = ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 4,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MEVFin;
    if ((p = m_check1(p)) == NULL)
        goto MEVFin;
    if (row != col) {
        mat_err(1);
        goto MEVFin;
    }
    if ((p = m_getmat(p,row,1,&idx,cmd,0)) == NULL)
        goto MEVFin;
    if ((p = m_check1(p)) == NULL)
        goto MEVFin;
    if ((p = m_getmat(p,row,1,&idx1,cmd,0)) == NULL)
        goto MEVFin;
    if ((p = m_check1(p)) == NULL)
        goto MEVFin;
    if ((p = m_getmat(p,row,col,&idx2,cmd,0)) == NULL)
        goto MEVFin;
    if ((p = m_check1(p)) == NULL)
        goto MEVFin;
    if ((p = m_getmat(p,row,col,&idx3,cmd,1)) == NULL)
        goto MEVFin;

    if (!(indic = (int *)calloc(row + 1,sizeof(int)))) {
        mat_err(3);
        goto MEVFin;
    }                
    memrq(row + 1,sizeof(int));
    indica = row + 1;

    err = eigen1(row,MX[0],MatVal[idx],MatVal[idx1],MatVal[idx2],MatVal[idx3],indic);
     
    if (err > 0) {
        m_cmdmsg();
        printf1("Error: no success in eigenvalue/vector calculations.\n");
        goto MEVFin;
    }
    else if (err < 0) {
        mat_err(3);
        goto MEVFin;
    }
    err = 0;

MEVFin:
    if (indica > 0)
        free((char *)indic);
    memrq(-indica,sizeof(int));
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_ginv(cmd)     mginv(A,R)                                              */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_ginv(char *cmd)
{
    int err,row,col,idx,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 6,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MGINVFin;
    if ((p = m_check1(p)) == NULL)
        goto MGINVFin;
    if (row < col) {
        mat_err(2);
        goto MGINVFin;
    }
    if ((p = m_getmat(p,col,row,&idx,cmd,1)) == NULL)
        goto MGINVFin;

    m_copy(MatVal[idx],MX[0],row * col);
    err = ginv(row,col,MatVal[idx]);
    if (err < 0) {
        if (err == -1)
            mat_err(3);
        else           
            mat_err(2);
        goto MGINVFin;
    }
    printf1("Pseudorank of %s: %d\n",MXName,err);
    err = 0;

MGINVFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_svd(cmd,opt)      opt 0 :  msvd(A,Q)                                  */
/*                      opt 1 :  msvd1(A,Q,U,V)                             */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_svd(char *cmd,int opt)
{
    int err,row,col,idx,idx1,idx2,no,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = cmd + 5;
    no = 1;
    if (opt) {
        p++;
        no = 0;
    }   
    p = n_mexpr(p,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MSVDFin;
    if (row < col) {
        mat_err(2);
        goto MSVDFin;
    }
    if ((p = m_check1(p)) == NULL)
        goto MSVDFin;
    if ((p = m_getmat(p,col,1,&idx,cmd,no)) == NULL)
        goto MSVDFin;

    if (opt) {
        if ((p = m_check1(p)) == NULL)
            goto MSVDFin;
        if ((p = m_getmat(p,row,col,&idx1,cmd,0)) == NULL)
            goto MSVDFin;
        if ((p = m_check1(p)) == NULL)
            goto MSVDFin;
        if ((p = m_getmat(p,col,col,&idx2,cmd,1)) == NULL)
            goto MSVDFin;

        err = svdecomp(row,col,MX[0],MatVal[idx],MatVal[idx1],MatVal[idx2],3);
    }
    else {
        if (alloc_acu(row * col + 1))
            goto MSVDFin;
        if (alloc_acv(col * col + 1))
            goto MSVDFin;
        err = svdecomp(row,col,MX[0],MatVal[idx],AcU,AcV,0);    
        alloc_acu(0);
        alloc_acv(0);
    }
    if (err < 0) {
        if (err == -1)
            mat_err(3);
        else           
            mat_err(2);
        goto MSVDFin;
    }
    err = 0;

MSVDFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_wvec(cmd)              mwvec(A,W,R)                                   */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_wvec(char *cmd)
{
    register int i,j;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp,tmp1,x0,x1;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 6,2,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MWVECFin;
    if (col != 1) {
        mat_err(7);
        goto MWVECFin;
    }
    if ((p = m_check1(p)) == NULL)
        goto MWVECFin;
    if ((p = m_getmat(p,row,1,&idx,cmd,1)) == NULL)
        goto MWVECFin;

    for (i = 1; i < row; ++i) {
        tmp = tmp1 = 0.0;
        for (j = i + 1; j <= row; ++j) {
            x0 = get_mxval(0,j,1);
            x1 = get_mxval(1,j,1);
            tmp += x0 * x1;                         
            tmp1 += x1;
        }
        if (tmp1 != 0.0)  
            tmp /= tmp1;
        else
            tmp = get_mxval(0,i,1);
        MatVal[idx][i] = tmp;
    }
    MatVal[idx][row] = get_mxval(0,row,1);
    err = 0;

MWVECFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_wvec1(cmd)            mwvec1(A,W,T,R)                                 */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_wvec1(char *cmd)
{
    register int i,j;
    int err,row,col,idx,ivflg;
    register char *p;
    double a,tmp,tmp1,tmp2,x0,x1;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 7,3,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MWVEC1Fin;
    if (col != 1) {
        mat_err(7);
        goto MWVEC1Fin;
    }
    if ((p = m_check1(p)) == NULL)
        goto MWVEC1Fin;
    if ((p = m_getmat(p,row,1,&idx,cmd,1)) == NULL)
        goto MWVEC1Fin;

    for (i = 1; i <= row; ++i) {
       
        tmp = tmp1 = 0.0;
        a = get_mxval(2,i,1);
        for (j = 1; j <= row; ++j) {
            tmp2 = get_mxval(2,j,1);
            if (tmp2 <= a)
                continue;
            x0 = get_mxval(0,j,1);
            x1 = get_mxval(1,j,1);
            tmp += x0 * x1;
            tmp1 += x1;
        }
        if (tmp1 != 0.0)  
            tmp /= tmp1;
        else
            tmp = get_mxval(0,i,1);
        MatVal[idx][i] = tmp;
    }
    err = 0;

MWVEC1Fin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_scal1(cmd)            mscal1(A,R)                                     */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_scal1(char *cmd)
{
    register int i,j;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 7,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MSCAL1Fin;
        
    if ((p = m_check1(p)) == NULL)
        goto MSCAL1Fin;
    if ((p = m_getmat(p,row,col,&idx,cmd,1)) == NULL)
        goto MSCAL1Fin;

    tmp = 0.0;
    for (i = 1; i <= row; ++i) {
        for (j = 1; j <= col; ++j)  
            tmp += get_mxval(0,i,j);
    }
    if (tmp != 0.0) {
        for (i = 1; i <= row; ++i) {
            for (j = 1; j <= col; ++j)  
                MatVal[idx][(i - 1) * col + j] = get_mxval(0,i,j) / tmp;
            }
    }
    err = 0;

MSCAL1Fin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mdefc(cmd)            mdefc(m,n,x,A)                                  */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mdefc(char *cmd)
{
    register int i;
    int err,row,col,idx,m,n,ivflg;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 6,3,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MDEFCFin;
    if (row != 1 || col != 1) {
        mat_err(8);
        goto MDEFCFin;
    }
    m = (int)get_mxval(0,1,1);
    n = (int)get_mxval(1,1,1);
    if (m < 1 || n < 1) {
        mat_err(2);
        goto MDEFCFin;
    }
    if ((p = m_check1(p)) == NULL)
        goto MDEFCFin;
    if ((p = m_getmat(p,m,n,&idx,cmd,1)) == NULL)
        goto MDEFCFin;

    m *= n;
    tmp = get_mxval(2,1,1);
    for (i = 1; i <= m; ++i)
        MatVal[idx][i] = tmp;
    err = 0;

MDEFCFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mdefi(cmd)            mdefi(m,n,A)                                    */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mdefi(char *cmd)
{
    register int i,j;
    int err,row,col,idx,m,n,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 6,2,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MDEFIFin;
    if (row != 1 || col != 1) {
        mat_err(8);
        goto MDEFIFin;
    }
    m = (int)get_mxval(0,1,1);
    n = (int)get_mxval(1,1,1);
    if (m < 1 || n < 1) {
        mat_err(2);
        goto MDEFIFin;
    }
    if ((p = m_check1(p)) == NULL)
        goto MDEFIFin;
    if ((p = m_getmat(p,m,n,&idx,cmd,1)) == NULL)
        goto MDEFIFin;

    j = imin(m,n);
    for (i = 0; i < j; ++i)  
        MatVal[idx][i * n + i + 1] = 1.0;
    err = 0;

MDEFIFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_invs(cmd)             minvs(A,R)                                      */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_invs(char *cmd)
{
    int err,row,col,idx,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 6,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MINVSFin;
    if (row != col) {
        mat_err(1);
        goto MINVSFin;
    }
    if ((p = m_check1(p)) == NULL)
        goto MINVSFin;
    if ((p = m_getmat(p,row,col,&idx,cmd,1)) == NULL)
        goto MINVSFin;

    m_copy(MatVal[idx],MX[0],row * col);
    col = syminv(row,MatVal[idx]);
    if (col < 0) {
        mat_err(3);
        goto MINVSFin;
    }
    if (col < row) {
        m_cmdmsg();
        printf1("%s not positive definite.\n",MXName);
        goto MINVSFin;
    }
    err = 0;

MINVSFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_invd(cmd)             minvd(A,R)                                      */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_invd(char *cmd)
{
    register int i,k;
    int err,row,col,idx,ivflg;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 6,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MINVDFin;
    if ((p = m_check1(p)) == NULL)
        goto MINVDFin;
    if ((p = m_getmat(p,row,col,&idx,cmd,1)) == NULL)
        goto MINVDFin;

    k = imin(row,col);
    for (i = 1; i <= k; ++i) {
        tmp = get_mxval(0,i,i);
        if (tmp)
            MatVal[idx][(i - 1) * col + i] = 1.0 / tmp;
    }
    err = 0;

MINVDFin:
    mx_free();
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

int m_ple(char *cmd)
{
    register int i,j,k;
    int err,m,row,col,row1,col1,idx,idx1,ivflg;
    register char *p;
    double tmp,tmp1;
    
    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    if ((p = n_mexpr(cmd + 5,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MPLEFin;
    if ((p = m_check1(p)) == NULL)
        goto MPLEFin;
    if ((p = n_mexpr(p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
        goto MPLEFin;
    if (col != 1 || col1 != 1) {
        mat_err(7);
        goto MPLEFin;
    }
    if (row != row1) {
        mat_err(2);
        goto MPLEFin;
    }
    if ((p = m_check1(p)) == NULL)
        goto MPLEFin;
    if ((p = m_getmat(p,row,1,&idx,cmd,0)) == NULL)
        goto MPLEFin;
    if ((p = m_check1(p)) == NULL)
        goto MPLEFin;
    if ((p = m_getmat(p,row,1,&idx1,cmd,1)) == NULL)
        goto MPLEFin;

    if (alloc_acn(row + 1))
        goto MPLEFin;
    if (alloc_acns(row + 1))
        goto MPLEFin;
    if (alloc_actmp(row + 1))
        goto MPLEFin;
    
    for (i = 1; i <= row; ++i) {
        if (MX[1][i] == 0.0)
            AcNS[i - 1] = 1;                /* not censored */
    }

    if (sortdp2a(row,MX[0] + 1,AcNS,AcN,0))      /* sort */
        return(-1);    

    tmp = 1.0 / (double)row;
    m = row - 1;
    for (i = 0; i < row; ++i) {
        j = AcN[i];
        AcTmp[j] += tmp;
        if (AcNS[j] == 0 && i < row - 1) {
            tmp1 = AcTmp[j] / (double)m;
            for (k = i + 1; k < row; ++k)  
                AcTmp[AcN[k]] += tmp1;
            AcTmp[j] = 0.0;
        }
        m--;
    }
    tmp = 0.0;
    for (i = 0; i < row; ++i) {
        k = AcN[i];
        tmp += AcTmp[k];
        MatVal[idx][k + 1] = tmp;
        MatVal[idx1][k + 1] = AcTmp[k];
    }
    err = 0;

MPLEFin:
    alloc_actmp(0);
    alloc_acns(0);
    alloc_acn(0);
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_nvar(cmd)             mnvar(X)    Replace data matrix by X.           */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_nvar(char *cmd)
{
    register int i,j;
    int err,row,col,ivflg;
    register char *p;
    char vname1[4 * VNLMax+1];

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    if (DMDef) {
        m_cmdmsg();
        printf1("Error: a data matrix already exists.\n");
        err = 0;
        goto MNVARFin;
    }
    p = n_mexpr(cmd + 6,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MNVARFin;
    if ((p = m_check2(p,0)) == NULL)
        goto MNVARFin;

    if (col > MaxNV) {
        m_cmdmsg();
        printf1("Error: %d variables exceeds current maximum.\n",col);
        goto MNVARFin;
    }
    for (i = 0; i < col; ++i) {
        sprintf(vname1,"%s%d<8>[%d.%d]=M:%s%d",MXName,i+1,PMMFmt1,PMMFmt2,MXName,i+1);
        if (save_var(vname1,0)) {       /* save variable definition */
            m_cmdmsg();
            printf1("Error: can't save variable definitions.\n");
            goto MNVARFin;              
        }
    }
    /** prn_var(0); **/  /* print list of new variables */

    VLabelLen = 0;
    NOCMaxA = NOCDM = NOC = row;

    /* allocate memory for data */

    if (alloc_vdat(0,1)) {
        m_cmdmsg();
        p_err(-2,1);
        goto MNVARFin;
    }
    for (i = 0; i < row; ++i) {
        for (j = 0; j < NVAR; ++j)  
            put_data(MX[0][i * col + j + 1],j,i);
    }
    DMDef = 1;
    printf2("Created data matrix with %d variables and %d cases.\n",NVAR,NOCDM);
    err = 0;

MNVARFin:
    mx_free();
    if (err) {
        clear_dm();            /* clear all variables */
        NVAR = NOCMaxA = NOCDM = NOC = 0;
    }
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_print(cmd)        mpr[1][a](X [,string] ) [=fname]                    */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_print(char *cmd)
{
    int err,row,col,aflg,bflg,ivflg;
    register char *p,*q,*s;

    err = -1;
    ivflg = 1;

    printf2("%s\n",cmd);
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
    q = skip_blev(p);
    if (*--q != ')') {
        mat_err(0);
        goto MPRFin;
    }
    p = n_mexpr(++p,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MPRFin;

    ivflg = MXIV[0];
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
        err = mpr(row,col,MX[0],ivflg,MX1[0],PMATFmtS,p,s,aflg);
    else
        err = mpr(row,col,MX[0],ivflg,MX1[0],PMATFmtS,NULL,s,aflg);

MPRFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mkp(cmd)      mkp(A,B,R)      Kronecker product                       */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mkp(char *cmd)
{
    register int i,j,k,l;
    int err,row,col,row1,col1,idx,n,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    if ((p = n_mexpr(cmd + 4,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MKPFin;
    if ((p = m_check1(p)) == NULL)
        goto MKPFin;
    if ((p = n_mexpr(p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
        goto MKPFin;
    if ((p = m_check1(p)) == NULL)
        goto MKPFin;
    if ((p = m_getmat(p,row * row1,col * col1,&idx,cmd,1)) == NULL)
        goto MKPFin;
    n = col * col1;
    for (i = 0; i < row; ++i) {
        for (l = 0; l < row1; ++l) {
            for (j = 0; j < col; ++j) {
                for (k = 0; k < col1; ++k) {
                    MatVal[idx][(i * row1 +l) * n + j * col1 + k + 1] =
                        get_mxval(0,i + 1,j + 1) * get_mxval(1,l + 1,k + 1);
                }
            }
        }
    }
    err = 0;

MKPFin:
    mx_free();
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

int m_mls(char *cmd,int opt)
{
    int err,row,col,idx,ra,re,ivflg,me,mi,ma,row1,col1;
    register char *p;
    double rne,rnl;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = cmd + 4;
    if (opt)
        p++;
    if (opt == 3)
        p++;
    if ((p = n_mexpr(p,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MLSFin;
    if ((p = m_check1(p)) == NULL)
        goto MLSFin;
    if (col < 2) {
        mat_err(2);
        goto MLSFin;
    }
    if (opt == 3 || opt == 4) {
        if ((p = n_mexpr(p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
            goto MLSFin;
        if ((p = m_check1(p)) == NULL)
            goto MLSFin;
        if (row1 != 1 || col1 != 1) {
            mat_err(8);
            goto MLSFin;
        }
        me = (int)MX[1][1];

        if ((p = n_mexpr(p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
            goto MLSFin;
        if ((p = m_check1(p)) == NULL)
            goto MLSFin;
        if (row1 != 1 || col1 != 1) {
            mat_err(8);
            goto MLSFin;
        }
        mi = (int)MX[1][1];

        if (me < 0 || mi < 0) {
            mat_err(14);
            goto MLSFin;
        }
        if (opt == 3) {
            ma = row - me - mi;
            if (ma < 0) {
                mat_err(14);
                goto MLSFin;
            }
        }
        else if (opt == 4) {
            ma = row - me;
            if (ma < 0) {
                mat_err(14);
                goto MLSFin;
            }
            if (mi > col - 1)
                mi = col - 1;
        }
    }

    if ((p = m_getmat(p,col - 1,1,&idx,cmd,1)) == NULL)
        goto MLSFin;

    if (opt == 1)         /* mlse */
        err = lsei(MX[0],0,row,0,col - 1,MatVal[idx],0,&rne,&rnl,&ra,&re);
                  
    else if (opt == 2)    /* mlsi */
        err = lsei(MX[0],0,0,row,col - 1,MatVal[idx],0,&rne,&rnl,&ra,&re);
                  
    else if (opt == 3)    /* mlsei */
        err = lsei(MX[0],me,ma,mi,col - 1,MatVal[idx],0,&rne,&rnl,&ra,&re);
                  
    else if (opt == 4)    /* mnls */
        err = wnnls(MX[0],col,me,ma,col - 1,mi,MatVal[idx],&rnl);
                  
    else                  /* mls */
        err = lsei(MX[0],0,row,0,col - 1,MatVal[idx],0,&rne,&rnl,&ra,&re);
                  
                  
    switch (err) {
        case  0: break;
        case -1: m_cmdmsg();
                 printf1("Equality constraints contradictory.\n");
                 goto MLSFin;
        case -2: m_cmdmsg();
                 printf1("Cannot satisfy inequality constraints.\n");
                 goto MLSFin;
        case -3: m_cmdmsg();
                 printf1("Cannot satisfy equality and inequality constraints.\n");
                 goto MLSFin;
        case -4: mat_err(3);
                 goto MLSFin;
        default: m_cmdmsg();
                 printf1("Undefined error (%d).\n",err);
                 goto MLSFin;
    }
    if (opt != 4)
        printf1("Rank of left-hand side: %d\n",ra);
    printf2("Norm of residuals: ");          
    printf2(PMATFmtS,rnl);          
    printf2("\n");          
    err = 0;

MLSFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mlp(cmd,opt)  opt=0:  mlp (T,X,Y)     linear programming              */
/*                  opt=1:  mlp1(T,p,X,Y)                                   */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mlp(char *cmd,int opt)
{
    int err,row,col,row1,col1,idx,idx1,np,r,ivflg;
    register char *p;
    double tmp;

    ivflg = np = 0;
    err = -1;
    printf2("%s\n",cmd);
    p = cmd + 4;
    if (opt)
        p++;
    if ((p = n_mexpr(p,1,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto MLPFin;
    if (row < 2 || col < 2) {
        mat_err(2);
        goto MLPFin;
    }
    if (opt) {
        if ((p = m_check1(p)) == NULL)
            goto MLPFin;
        if ((p = n_mexpr(p,1,1,&row1,&col1,&ivflg,NULL)) == NULL) 
            goto MLPFin;
        if (row1 != 1 || col1 != 1 || (np = (int)MX[1][1]) < 0 || np >= row) {
            m_cmdmsg();
            printf1("Error in second argument.\n");
            goto MLPFin;
        }
    }
    if ((p = m_check1(p)) == NULL)
        goto MLPFin;
    if ((p = m_getmat(p,col - 1,1,&idx,cmd,0)) == NULL)
        goto MLPFin;
    if ((p = m_check1(p)) == NULL)
        goto MLPFin;
    if ((p = m_getmat(p,row - 1,1,&idx1,cmd,1)) == NULL)
        goto MLPFin;

    r = lpf1(row-1,col-1,np,MX[0],MatVal[idx],MatVal[idx1],&tmp,1.e-8);
    if (r) {
        m_cmdmsg();
        switch (r) {
            case -2:    printf1("No solution.\n");
                        if (np >= 2) 
                            printf1("(Possibly linear dependent equality constraints.)\n");
                        break;
            case -3:    printf1("No solution. Dual objective function is unbounded.\n");
                        break;
            case -4:    printf1("No solution. Primal objective function is unbounded.\n");
                        break;
            case -5:    printf1("Error: insufficient memory.\n");
                        break;
            default:    printf1("Error (%d).\n",r);
                        break;
        }   
        goto MLPFin;
    }
    printf2("Value is: ");
    printf2(PMATFmtS,tmp);
    printf2("\n");
    err = 0;

MLPFin:
    mx_free();
    return(err);
}

/*--###---------------------------------------------------------------------*/
/*  m_mlpi(cmd)   linear programming with integer variables.                */
/*                                                                          */
/*  mlpi(A,B,lmax)                                                          */
/*                                                                          */
/*  lmax is maximal number of solutions. Solutions are always written       */ 
/*  into the standard output.                                               */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mlpi(char *cmd)
{
    register int i,j;
    int err,row,col,r,ivflg,vmax,m,n,nest;
    register char *p;

    ivflg = 0;
    err = -1;
    printf2("%s\n",cmd);
    p = cmd + 5;
    if ((p = n_mexpr(p,1,0,&m,&n,&ivflg,NULL)) == NULL) 
        goto MLPIFin;
    if (m < 2 || n < 2) {
        mat_err(2);
        goto MLPIFin;
    }
    if ((p = m_check1(p)) == NULL)
        goto MLPIFin;
    if ((p = n_mexpr(p,1,1,&row,&col,&ivflg,NULL)) == NULL) 
        goto MLPIFin;

    if (row != m || col != 1) {
        mat_err(2);  
        goto MLPIFin;
    }
    if ((p = m_check1(p)) == NULL)
        goto MLPIFin;

    if ((p = n_mexpr(p,1,2,&row,&col,&ivflg,NULL)) == NULL) 
        goto MLPIFin;

    if (row != 1 || col != 1) {
        mat_err(8);
        goto MLPIFin;
    }
    nest = (int)MX[2][1];
          
    if (alloc_acn(nest * n + 1))
        goto MLPIFin;

    if (alloc_acm(4 * m + 6 * n + m * n + 2))
        goto MLPIFin;

    if (alloc_acr(m * n + 1))
        goto MLPIFin;

    if (alloc_acs(m + 1))
        goto MLPIFin;

    for (i = 1; i <= m; ++i) {
        AcS[i] = MX[1][i];
        for (j = 1; j <= n; ++j)
            AcR[(i - 1) * n + j] = MX[0][(i - 1) * n + j];
    }
    r = lpi(n,m,nest,AcR,AcS,AcN,AcM,&vmax);

    if (r < 1) {
        m_cmdmsg();
        switch (r) {
            case -1:    printf1("inconsistent constraints.\n");
                        break;
            case -2:    printf1("exceeded maximal number of solutions.\n");
                        break;
            default:    printf1("Error (%d).\n",r);
                        break;
        }   
        goto MLPIFin;
    }
    printf1("Value: %d\nNumber of solutions: %d\n",vmax,r);
    for (i = 0; i < r; ++i) {
        for (j = 1; j <= n; ++j)
            printf("%2d ",AcN[i * n + j]);
        newline();
    }
    err = 0;

MLPIFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mul(cmd)      mmul(A1,A2,...,R)     multiplication                    */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mul(char *cmd)
{
    register int i,j,k,l;
    int err,row,col,row1,col1,idx,ivflg;
    register char *p,*q;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    if ((p = n_mexpr(cmd + 5,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MMULFin;

    while (1) {
        if ((p = m_check1(p)) == NULL)
            goto MMULFin;

        q = skip_expr(p);
        if (!*q || *q == ')')
            break;

        /* get next matrix into MX[1] */

        if ((p = n_mexpr(p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
            goto MMULFin;

        /* MX[0] * MX[1] -> V */

        if (alloc_acu(row * col1 + 1))
            goto MMULFin;

        k = imax(col,row1);
        for (i = 1; i <= row; ++i) {
            for (j = 1; j <= col1; ++j) {
                tmp = 0.0;
                for (l = 1; l <= k; ++l)
                    tmp += get_mxval(0,i,l) * get_mxval(1,l,j);
                AcU[(i - 1) * col1 + j] = tmp;
            }
        }

        /* copy U into MX[0], set dimension */

        col = col1;
        if (mx_alloc(0,row,col))
            goto MMULFin;
        m_copy(MX[0],AcU,row * col1);
    }
    if ((p = m_getmat(p,row,col,&idx,cmd,1)) == NULL)
        goto MMULFin;
    m_copy(MatVal[idx],MX[0],row * col);
    err = 0;

MMULFin:
    alloc_acu(0);
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_cat(cmd,opt)      opt 0 : mcath(A1,A2,...,R)      horizontal concat   */
/*                      opt 1 : mcatv(A1,A2,...,R)      vertical concat     */
/*                      opt 2 : mdcathv(A1,A2,...,R)    direct sum          */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_cat(char *cmd,int opt)
{
    register int i,j;
    int err,row,col,row1,col1,r,c,idx,ivflg;
    register char *p,*q;

    err = -1;
    ivflg = 0;

    r = 6;
    if (opt == 2)
        r++;       

    printf2("%s\n",cmd);
    if ((p = n_mexpr(cmd + r,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MCATFin;

    while (1) {
        if ((p = m_check1(p)) == NULL)
            goto MCATFin;

        q = skip_expr(p);
        if (!*q || *q == ')')
            break;

        /* get next matrix into MX[1] */

        if ((p = n_mexpr(p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
            goto MCATFin;

        /* MX[0] + MX[1] -> V */

        if (opt == 0) {
            r = imax(row,row1);
            c = col + col1;
        }
        else if (opt == 1) {
            r = row + row1;
            c = imax(col,col1);
        }
        else  {
            r = row + row1;
            c = col + col1;
        }
        if (alloc_acu(r * c + 1))
            goto MCATFin;

        if (opt == 0) {
            for (i = 1; i <= r; ++i) {
                for (j = 1; j <= col; ++j) 
                    AcU[(i - 1) * c + j] = get_mxval(0,i,j);

                for (j = col + 1; j <= c; ++j) 
                    AcU[(i - 1) * c + j] = get_mxval(1,i,j - col);
            }
        }
        else if (opt == 1) {
            for (j = 1; j <= c; ++j) {
                for (i = 1; i <= row; ++i) 
                    AcU[(i - 1) * c + j] = get_mxval(0,i,j);

                for (i = row + 1; i <= r; ++i) 
                    AcU[(i - 1) * c + j] = get_mxval(1,i - row,j);
            }
        }
        else  {
            for (i = 1; i <= row; ++i) {
                for (j = 1; j <= col; ++j) 
                    AcU[(i - 1) * c + j] = get_mxval(0,i,j);
            }
            for (i = row + 1; i <= r; ++i) {
                for (j = col + 1; j <= c; ++j) 
                    AcU[(i - 1) * c + j] = get_mxval(1,i - row,j - col);
            }
        }

        /* copy U into MX[0], set dimension */

        row = r;   
        col = c;
        if (mx_alloc(0,row,col))
            goto MCATFin;
        m_copy(MX[0],AcU,row * col);
    }
    if ((p = m_getmat(p,row,col,&idx,cmd,1)) == NULL)
        goto MCATFin;
    m_copy(MatVal[idx],MX[0],row * col);
    err = 0;

MCATFin:
    alloc_acu(0);
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_srow(cmd,opt)  opt=0: msrow(X,D,R)    select rows                     */
/*                   opt=1: mscol(X,D,R)    select columns                  */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_srow(char *cmd,int opt)
{
    register int i,j,k,l;
    int err,row,col,row1,col1,m,idx,ivflg,sflag;             
    register char *p,*q;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    if ((p = n_mexpr(cmd + 6,1,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto MSROWFin;
    if ((p = m_check1(p)) == NULL)
        goto MSROWFin;

    if (*p == '<') {
        sflag = 1;
        q = p;
        if (m_selidx(q,0,&m,AcI) == NULL || m < 1) {
            mat_err(9);
            goto MSROWFin;
        }
        if (alloc_aci(m + 1))
            goto MSROWFin;

        if ((p = m_selidx(p,1,&m,AcI)) == NULL || m < 1) {
            mat_err(9);
            goto MSROWFin;
        }
    }
    else {
        sflag = 0;

        if ((p = n_mexpr(p,1,1,&row1,&col1,&ivflg,NULL)) == NULL) 
            goto MSROWFin;

        if (row1 != 1 && col1 != 1) {
            mat_err(9);
            goto MSROWFin;
        }
        m = imax(row1,col1);
    }
    if ((p = m_check1(p)) == NULL)
        goto MSROWFin;

    if (opt == 0) {
        if ((p = m_getmat(p,m,col,&idx,cmd,1)) == NULL)
            goto MSROWFin;

        l = 0;
        for (i = 1; i <= m; ++i) {
            if (sflag)
                k = AcI[i - 1];
            else
                k = (int)MX[1][i];

            if (k < 1 || (opt == 0 && k > row) || (opt == 1 && k > col)) {
                mat_err(9);
                goto MSROWFin;
            }
            for (j = 1; j <= col; ++j)
                MatVal[idx][l * col + j] = get_mxval(0,k,j);
            l++; 
        }
    }
    else {
        if ((p = m_getmat(p,row,m,&idx,cmd,1)) == NULL)
            goto MSROWFin;

        for (j = 1; j <= m; ++j) {
            if (sflag)
                k = AcI[j - 1];
            else
                k = (int)MX[1][j];

            if (k < 1 || (opt == 0 && k > row) || (opt == 1 && k > col)) {
                mat_err(9);
                goto MSROWFin;
            }
            for (i = 0; i < row; ++i)
                MatVal[idx][i * m + j] = get_mxval(0,i + 1,k);
        }
    }
    err = 0;

MSROWFin:
    alloc_aci(0); 
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_selidx(s,opt,n,idx)   s points to a string of the form                */
/*                          <i1,i2,...> or <i1,i2,,i3,...>                  */
/*                                                                          */
/*  Return the number of indices in n. Return pointer to next character     */
/*  after >, or NULL if error.                                              */

char *m_selidx(char *s,int opt,int *n,int *idx)
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
            s = skip_int(s);
            s = skip_int(s + 2);
        }
        else if (sscanf(s,"%d",&nc) == 1 && nc >= 1) {
            ni++;
            s = skip_int(s);
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

int m_sort(char *cmd,int opt)
{
    register int i,j,k,l,ll;
    int err,row,col,row1,col1,m,idx,ivflg;
    register char *p;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = cmd + 6;
    if (opt == 2)
        p++;
    if ((p = n_mexpr(p,1,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto MSORTFin;
    if ((p = m_check1(p)) == NULL)
        goto MSORTFin;
    if ((p = n_mexpr(p,1,1,&row1,&col1,&ivflg,NULL)) == NULL) 
        goto MSORTFin;
    if (row1 != 1 && col1 != 1) {
        mat_err(9);
        goto MSORTFin;
    }
    m = imax(row1,col1);
    for (i = 1; i <= m; ++i) {
        k = (int)MX[1][i];
        if (k < 1 || k > col) {
            mat_err(9);
            goto MSORTFin;
        }
    }
    if ((p = m_check1(p)) == NULL)
        goto MSORTFin;

    if (alloc_acn(m + 1))        /* save indices for sort */
        goto MSORTFin;

    for (i = 1; i <= m; ++i)
        AcN[i] = (int)MX[1][i];

    if (alloc_aci(row + 1))         /* pointer for sorting */
        goto MSORTFin;

    if (sortdpn(row,col,MX[0],m,AcN,AcI))
        goto MSORTFin;      

    if (opt == 1) {                                         /* rank */
        if ((p = m_getmat(p,row,1,&idx,cmd,1)) == NULL)
            goto MSORTFin;

        for (i = 1; i <= row; ++i)
            MatVal[idx][i] = (double)AcI[i];
    }
    else {
        if (alloc_actmp(row * col + 1))     
            goto MSORTFin;

        k = 1;
        for (i = 1; i <= row; ++i) {
            l = (AcI[i] - 1) * col;
            ll = (k - 1) * col;
            for (j = 1; j <= col; ++j) 
                AcTmp[++ll] = MX[0][++l];
                       
            if (opt == 2 && k > 1) {
                k--;
                for (j = 1; j <= col; ++j) {
                    if (AcTmp[k * col + j] != AcTmp[(k-1) * col + j]) {
                        k++;
                        break;
                    }
                }   
            }
            k++;
        }
        k--;
        if ((p = m_getmat(p,k,col,&idx,cmd,1)) == NULL)
            goto MSORTFin;
        m_copy(MatVal[idx],AcTmp,k * col);
    }
    err = 0;

MSORTFin:
    alloc_actmp(0);
    alloc_acn(0);
    alloc_aci(0);
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  get_mselect(p)  get <...> expression into AcC.                          */
/*                                                                          */
/*  Return pointer to next character, or NULL if error                      */

char *get_mselect(char *p)
{
    int l;
    register char *q;

    l = strlen(p);
    if (alloc_acc(l))   
        return(NULL);

    q = AcC;
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

int m_exp(char *cmd)
{
    register int i;
    int err,r,row,col,idx;
    char *p,*q;
    char mtmp[MatDefLen + 1];

    err = -1;
    printf2("%s\n",cmd);
    mdefcpy(mtmp,cmd);

    p = cmd + 5;

    r = strlen(cmd);
    q = cmd + r;               
    if (*--q != ')') {
        mat_err(0);
        goto MEXPFin;
    }
    *q = '\0';

    for (i = r; i >= 5; --i) {
        if (*--q == ',') 
            break;
    }
    if (*q != ',' || !*(q + 1)) {
        mat_err(0);
        goto MEXPFin;
    }
    *q++ = '\0';

    if (mex_eval(p))
        goto MEXPFin;

    row = MEXOPRow[0];
    col = MEXOPCol[0];

    if ((idx = mat_newmat(q,row,col)) < 0)
        goto MEXPFin;

    m_copy(MatVal[idx],MEXOP[0],row * col);
    strcpy(MatDef[idx],mtmp);
    err = 0;

MEXPFin:       
    alloc_mex_eval(0,0,0);                  
    alloc_mex(0,0);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_expr(cmd)     mexpr(expression,A)                                     */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_expr(char *cmd)
{
    int err,row,col,idx,ivflg,sflag;
    register char *p;
    char *sel;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);

    sflag = 0;
    sel = NULL;

    p = cmd + 6;
    if (*p == '<') {
        p = get_mselect(p + 1);
        if (p == NULL) {
            return(-1);
        }
        sflag = 1;
        sel = AcC;
    }

    if ((p = n_mexpr(p,1,0,&row,&col,&ivflg,sel)) == NULL) 
        goto MEXPRFin;
    if ((p = m_check1(p)) == NULL)
        goto MEXPRFin;
    if ((p = m_getmat(p,row,col,&idx,cmd,1)) == NULL)
        goto MEXPRFin;

    m_copy(MatVal[idx],MX[0],row * col);
    err = 0;

MEXPRFin:
    alloc_acc(0);   
    mx_free();
    return(err);
}
/*--##----------------------------------------------------------------------*/
/*  m_expr1(cmd)     mexpr1(B,expression,A)                                 */
/*                                                                          */
/*  Evaluate expression in block mode, blocks defined by B.                 */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_expr1(char *cmd)
{
    register int j,k,k1,l;
    int err,ivflg,row,col,row1,col1,idx;
    register char c,*p,*q;
    double *tmp;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    if ((p = n_mexpr(cmd + 7,1,0,&row1,&col1,&ivflg,NULL)) == NULL) 
        goto MEXPR1Fin;
    if ((p = m_check1(p)) == NULL)
        goto MEXPR1Fin;

    q = skip_expr(p);
    c = *q;
    *q = '\0';
   
    if (get_mexpr(p,&row,&col,&ivflg) || row < 1 || col < 1) 
        goto MEXPR1Fin;

    row = imax(row,row1);
    if (row1 != row) {
        mat_err(2);
        goto MEXPR1Fin;
    }
    *q = c;
    if ((q = m_check1(q)) == NULL)
        goto MEXPR1Fin;
    if ((q = m_getmat(q,row,col,&idx,cmd,1)) == NULL)
        goto MEXPR1Fin;

    if (alloc_actmp(col1 + 1))
        goto MEXPR1Fin;

    for (j = 1; j <= col1; ++j)
        AcTmp[j] = MX[0][j];         

    k1 = 0;
    k = 1;
    while (k <= row) {
        if (k == row)  
            l = 1;
        else {
            l = 0;
            for (j = 1; j <= col1; ++j) {
                if (AcTmp[j] != MX[0][k * col1 + j]) {
                    l = 1;
                    break;
                }
            }
        }
        if (l == 0)
            k++;
        else {
            if (eval_mexpr(k1,k,col,MatVal[idx],0,tmp,NULL,NULL,-1,0.0))
                goto MEXPR1Fin;

            for (j = 1; j <= col1; ++j)  
                AcTmp[j] = MX[0][k * col1 + j];   
            k1 = k;
            k++;
        }
    }
    err = 0;

MEXPR1Fin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_setv(cmd)  Syntax: msetv(expression,X(row,col))                       */
/*               Sets X(row,col) to the value of expression (must be a      */
/*               scalar expression). X must be a matrix name.               */
/*                                                                          */
/*               return 0 if OK, -1 if error.                               */

int m_setv(char *cmd)  
{
    int err,idx,row,col,ivflg;
    char *p,mname[VNLMax + 1];  
   
    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    if ((p = n_mexpr(cmd + 6,1,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto MSETVFin;
    if ((p = m_check1(p)) == NULL)
        goto MSETVFin;
    if (row != 1 || col != 1) {
        mat_err(8);
        goto MSETVFin;       
    }
    if ((p = get_mname(p,mname,1)) == NULL)
        return(-1);
    if ((idx = mat_getidx(mname,1)) < 0)
        return(-1);
    if (*p++ != '(') {
        mat_err(0);
        goto MSETVFin;
    }
    if (sscanf(p,"%d,%d",&row,&col) == 2) {
        p = skip_int(p);
        p = skip_int(++p);
    }
    else {
        if ((p = n_mexpr(p,2,1,&row,&col,&ivflg,NULL)) == NULL) 
            goto MSETVFin;
        if (row != 1 || col != 1) {
            mat_err(8);
            goto MSETVFin;
        }
        row = (int)MX[1][1];
        col = (int)MX[2][1];
    }
    if ((p = m_check2(p,0)) == NULL)
        goto MSETVFin;
    if ((p = m_check2(p,1)) == NULL)
        goto MSETVFin;
    if (row < 1 || row > MatRow[idx] || col < 1 || col > MatCol[idx]) {
        mat_err(10);
        goto MSETVFin;
    }
    MatVal[idx][(row - 1) * MatCol[idx] + col] = MX[0][1];
    err = 0;

MSETVFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_brr(cmd)    Create matrix with BRR indicators.                        */
/*                mbrr(ns,nu,A)                                             */
/*                where ns is number of strata, nu is number of secu.       */
/*                A has dimension ns x nr, nr is number of replications.    */
/*                Return 0 if OK, -1 if error.                              */

int m_brr(char *cmd)  
{
    register int i,j;
    int err,idx,row,col,ns,nu,nr,ivflg;
    char *p;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    if ((p = n_mexpr(cmd + 5,2,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto MBRRFin;
    if ((p = m_check1(p)) == NULL)
        goto MBRRFin;
    if (row != 1 || col != 1) {
        mat_err(8);
        goto MBRRFin;       
    }
    ns = (int)MX[0][1];
    nu = (int)MX[1][1];
    if (SILENTFlg < 0)  
        printf1("Number of strata: %d, secu: %d\n",ns,nu);

    err = brr_gen(ns,nu,&nr,1);
    if (err || nr < 1)
        goto MBRRFin;

    if ((p = m_getmat(p,ns,nr,&idx,cmd,1)) == NULL)
        goto MBRRFin;

    for (i = 0; i < ns; ++i) {
        for (j = 0; j < nr; ++j)
            MatVal[idx][i * nr + j + 1] = (double)AcC[j * ns + i + 1];
    }
    err = 0;

MBRRFin:
    mx_free();
    alloc_acc(0);           /* used by brr_gen() */
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mtrim(A,ca,ra,cb,rb,R)   Create new matrix.                           */
/*                                                                          */
/*      ca = leading columns,cb trailing columns                            */
/*      ra = leading rows, rb trailing rows.                                */
/*      if positive delete, if negative add zero rows/columns               */

int m_mtrim(char *cmd)  
{
    register int i,j;
    int err,idx,row,col,row1,col1,ivflg,ca,cb,ra,rb;
    int ia,ib,ii,ja,jb,jj;
    char *p;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    if ((p = n_mexpr(cmd + 6,1,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto MTRIMFin;
    if ((p = m_check1(p)) == NULL)
        goto MTRIMFin;

    if ((p = n_mexpr(p,4,1,&row1,&col1,&ivflg,NULL)) == NULL) 
        goto MTRIMFin;
    if ((p = m_check1(p)) == NULL)
        goto MTRIMFin;

    if (row1 != 1 || col1 != 1) {
        mat_err(8);
        goto MTRIMFin;       
    }
    ca = (int)MX[1][1];
    ra = (int)MX[2][1];
    cb = (int)MX[3][1];
    rb = (int)MX[4][1];

    row1 = row - ra - rb;
    col1 = col - ca - cb;

    if (row1 <= 0 || col1 <= 0) {
        mat_err(15);
        goto MTRIMFin;       
    }
    if ((p = m_getmat(p,row1,col1,&idx,cmd,1)) == NULL)
        goto MTRIMFin;
       

    ia = imax(1 + ra,1);
    ib = imin(row - rb,row);
    ii = imax(1 - ra,1);
    for (i = ia; i <= ib; ++i) {
        ja = imax(1 + ca,1);
        jb = imin(col - cb,col);
        jj = imax(1 - ca,1);
        for (j = ja; j <= jb; ++j)
            MatVal[idx][(ii - 1) * col1 + jj++] = MX[0][(i - 1) * col + j];
        ii++; 
    }
    err = 0;

MTRIMFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mmp(cmd,opt)  Return monotone projection of X in Y.                   */
/*                                                                          */
/*  opt = 1:    mmp (X,Y)     ignore ties                                   */
/*  opt = 2:    mmp1(X,T,Y)   Kruskal's primary approach                    */
/*  opt = 3:    mmp1(X,T,Y)   Kruskal's secondary approach                  */
/*                                                                          */
/*  If opt = 2 or 3, then T contains information about ties.                */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mmp(char *cmd,int opt) 
{
    register int i,j,k;
    int err,row,col,row1,col1,idx,ivflg;
    register char *p;
    double tmp,tmp1;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    i = 4;
    if (opt != 0)
        i++;

    p = n_mexpr(cmd + i,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MMPFin;
    if ((p = m_check1(p)) == NULL)
        goto MMPFin;

    if (opt != 0) {
        p = n_mexpr(p,1,1,&row1,&col1,&ivflg,NULL);
        if (p == NULL)
            goto MMPFin;
        if ((p = m_check1(p)) == NULL)
            goto MMPFin;

        if (row1 != row || col1 != col) {
            mat_err(0);
            goto MMPFin;
        }
        if (alloc_ack(row + 1))
            goto MMPFin;
        if (alloc_actmp1(row + 1))
            goto MMPFin;

        if (opt == 1) {
            if (alloc_acj(row + 1))
                goto MMPFin;
            if (alloc_acn(row + 1))
                goto MMPFin;
        }   
    }

    if ((p = m_getmat(p,row,col,&idx,cmd,1)) == NULL)
        goto MMPFin;

    if (alloc_actmp(row + 1))
        goto MMPFin;
         
    for (j = 1; j <= col; ++j) {

        for (i = 1; i <= row; ++i)  
            AcTmp[i - 1] = get_mxval(0,i,j);
      
        if (opt == 0)
            mreg_proj(row,AcTmp);

        else {      
            k = 0;
            AcK[k] = 1;
            tmp = get_mxval(1,1,j);

            for (i = 1; i < row; ++i) {
                tmp1 = get_mxval(1,i + 1,j);
                if (fabs(tmp - tmp1) <= EPSI1)
                    AcK[k] += 1;            
                else {
                    k++;
                    AcK[k] = 1;
                }
                tmp = tmp1;
            }
            if (opt == 1) {
                if (mreg_proj1(row,k + 1,AcK,AcJ,AcN,AcTmp,AcTmp1)) {
                    mat_err(3);
                    goto MMPFin;
                }
            }
            else if (opt == 2) {
                if (mreg_proj2(row,k + 1,AcK,AcTmp,AcTmp1)) {
                    mat_err(3);
                    goto MMPFin;
                }
            }

        }
        for (i = 0; i < row; ++i)  
            MatVal[idx][i * col + j] = AcTmp[i];
        
    }
    err = 0;

MMPFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mpit(cmd,opt)      Iteration of Leslie matrix.                        */
/*                                                                          */
/*  opt = 0: mpit(F,N,n,R)                                                  */
/*  opt = 1: mpit1(F,N,Z,n,R)                                               */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mpit(char *cmd,int opt)
{
    register int i,j;
    int err,n,m,row,col,idx,ivflg;
    register char *p;
    double tmp;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);

    if (opt == 0) {
        p = cmd + 5;
        if ((p = n_mexpr(p,3,0,&row,&col,&ivflg,NULL)) == NULL) 
            goto MPITFin;
        if ((p = m_check1(p)) == NULL)
            goto MPITFin;

        n = MXRow[0];
        if (MXCol[0] != 2 || MXRow[1] != n || MXRow[2] != 1 || MXCol[2] != 1) {
            mat_err(0);
            goto MPITFin;
        }
        m = (int)MX[2][1] + 1;          /* iterations */
    }
    else {
        p = cmd + 6;
        if ((p = n_mexpr(p,4,0,&row,&col,&ivflg,NULL)) == NULL) 
            goto MPITFin;
        if ((p = m_check1(p)) == NULL)
            goto MPITFin;

        n = MXRow[0];
        if (MXCol[0] != 2 || MXRow[1] != n || MXRow[2] != n || MXRow[3] != 1 || MXCol[3] != 1) {
            mat_err(0);
            goto MPITFin;
        }
        m = (int)MX[3][1] + 1;          /* iterations */
    }
    if (m < 1)
        m = 1;

    if ((p = m_getmat(p,m,n,&idx,cmd,1)) == NULL)
        goto MPITFin;

    for (i = 1; i <= n; ++i)
        MatVal[idx][i] = MX[1][i];

    for (j = 1; j < m; ++j) {
        tmp = 0.0;
        for (i = 1; i <= n; ++i)  
            tmp += MX[0][(i - 1) * 2 + 1] * MatVal[idx][(j - 1) * n + i];
        MatVal[idx][j * n + 1] = tmp;
        for (i = 2; i <= n; ++i) 
            MatVal[idx][j * n + i] = MX[0][(i - 1) * 2] * MatVal[idx][(j - 1) * n + i - 1];

        if (opt == 1) {
            for (i = 1; i <= n; ++i)
                MatVal[idx][j * n + i] += MX[2][i];
        }
    }
    err = 0;

MPITFin:
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mqp(cmd,opt)   opt=0:  mqp (C,D,X)              quadratic programming */
/*                   opt=1:  mqpb(C,D,XL,XU,X)                              */
/*                   opt=2:  mqpc(C,D,A,B,me,X)                             */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mqp(char *cmd,int opt)
{
    register char *p;
    int i,j,m,me,err,r,row,col,idx,ivflg,row1,col1,row2,col2,row3,col3;
    int maxit,xla,xua,iaa,nact,*iact;
    double vsmall,diag,*xl,*xu,*a,*b;

    err = -1;
    me = m = iaa = xla = xua = ivflg = 0;

    printf2("%s\n",cmd);
    p = cmd + 4;
    if (opt)
        p++;

    if ((p = n_mexpr(p,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MPQFin;
    if ((p = m_check1(p)) == NULL)
        goto MPQFin;

    if (row != col) {
        mat_err(1);
        goto MPQFin;
    }
    if ((p = n_mexpr(p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
        goto MPQFin;
    if ((p = m_check1(p)) == NULL)
        goto MPQFin;

    if (row1 != row || col1 != 1) {
        mat_err(2);
        goto MPQFin;
    }
    if (opt == 1) {
        if ((p = n_mexpr(p,1,2,&row1,&col1,&ivflg,NULL)) == NULL)
            goto MPQFin;
        if ((p = m_check1(p)) == NULL)
            goto MPQFin;
        if (row1 != row || col1 != 1) {
            mat_err(2);
            goto MPQFin;
        }
        if ((p = n_mexpr(p,1,3,&row1,&col1,&ivflg,NULL)) == NULL)
            goto MPQFin;
        if ((p = m_check1(p)) == NULL)
            goto MPQFin;
        if (row1 != row || col1 != 1) {
            mat_err(2);
            goto MPQFin;
        }
    }
    else if (opt == 2) {
        if ((p = n_mexpr(p,1,2,&m,&col1,&ivflg,NULL)) == NULL)
            goto MPQFin;
        if ((p = m_check1(p)) == NULL)
            goto MPQFin;
        if (col1 != row) {
            mat_err(2);
            goto MPQFin;
        }
        if ((p = n_mexpr(p,1,3,&row2,&col2,&ivflg,NULL)) == NULL)
            goto MPQFin;
        if ((p = m_check1(p)) == NULL)
            goto MPQFin;
        if (row2 != m || col2 != 1) {
            mat_err(2);
            goto MPQFin;
        }
        if ((p = n_mexpr(p,1,4,&row3,&col3,&ivflg,NULL)) == NULL)
            goto MPQFin;
        if ((p = m_check1(p)) == NULL)
            goto MPQFin;
        if (row3 != 1 || col3 != 1) {
            mat_err(2);
            goto MPQFin;
        }
        me = (int)MX[4][1];
        if (me < 0 || me > m) {
            printf1("Error in number of equality constraints.\n");
            goto MPQFin;
        }
    }     
    if ((p = m_getmat(p,row,1,&idx,cmd,1)) == NULL)
        goto MPQFin;

    maxit = 50 * (row + m);
    vsmall = 2.0 * EPSI;

    if (opt == 0 || opt == 2) {
                     
        /* allocate and initialize xl and xu for lower and upper bounds */

        if (!(xl = (double *)calloc(row + 1,sizeof(double)))) {
            mat_err(3);
            goto MPQFin;
        }                
        xla = row + 1;
        memrq(xla,sizeof(double));

        if (!(xu = (double *)calloc(row + 1,sizeof(double)))) {
            mat_err(3);
            goto MPQFin;
        }                
        xua = row + 1;
        memrq(xua,sizeof(double));

        for (i = 1; i <= row; ++i) {
            xu[i] = DBLMAX / 1000.0;
            xl[i] = -xu[i];             
        }
    }

    /* allocate iact */

    if (!(iact = (int *)calloc(row + 1,sizeof(int)))) {
        mat_err(3);
        goto MPQFin;
    }                
    iaa = row + 1;
    memrq(iaa,sizeof(int));

    switch (opt) {
        case 0: r = qld(row,0,0,1,maxit,a,b,MX[1],MX[0],MatVal[idx],
                    xl,xu,iact,&nact,vsmall,&diag);
                break;
        case 1:
                r = qld(row,0,0,1,maxit,a,b,MX[1],MX[0],MatVal[idx],
                    MX[2],MX[3],iact,&nact,vsmall,&diag);
                break;

        case 2:
                r = qld(row,m,me,1,maxit,MX[2],MX[3],MX[1],MX[0],MatVal[idx],
                    xl,xu,iact,&nact,vsmall,&diag);
                break;

        default: goto MPQFin;
    }
    /*  printf("r=%d nact=%d diag=%g\n",r,nact,diag); */

    switch (r) {
        case  0: break;
        case  1: m_cmdmsg();
                 printf1("Exceeded maximal number of iterations (%d).\n",maxit);
                 goto MPQFin;
        case  2: m_cmdmsg();
                 printf1("Insufficient accuracy for convergence.\n");
                 goto MPQFin;
        case  3: m_cmdmsg();
                 printf1("Insufficient memory.\n");
                 goto MPQFin;
        default: m_cmdmsg();
                 printf1("Inconsistent constraints.\n");
                 goto MPQFin;
    }
    if (diag > 0.0)   
        printf2("Coeff matrix was enlarged by %g times unitmatrix.\n",diag);
    err = 0;

MPQFin:
    if (xla > 0) {
        free((char *)xl);
        memrq(-xla,sizeof(double));
    }
    if (xua > 0) {
        free((char *)xu);
        memrq(-xua,sizeof(double));
    }
    if (iaa > 0) {
        free((char *)iact);
        memrq(-iaa,sizeof(int));
    }
    mx_free();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_mlsei1(cmd)  mlsei1(S,me,mi,R)                                        */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int m_mlsei1(char *cmd)
{
    register int i,j,k;
    register char *p;
    int r,err,row,col,idx,ivflg,n,me,mi,ma,m,row1,col1,xxa,xya,aa,ba;
    int maxit,xla,xua,iaa,nact,*iact;
    double vsmall,diag,tmp,*xl,*xu,*xx,*xy,*a,*b;

    err = -1;
    xla = xua = iaa = aa = ba = xxa = xya = ivflg = 0;

    printf2("%s\n",cmd);
    p = cmd + 7;

    if ((p = n_mexpr(p,1,0,&row,&col,&ivflg,NULL)) == NULL)
        goto MLS1Fin;
    if ((p = m_check1(p)) == NULL)
        goto MLS1Fin;
    if (col < 2) {
        mat_err(2);
        goto MLS1Fin;
    }
    if ((p = n_mexpr(p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
        goto MLS1Fin;
    if ((p = m_check1(p)) == NULL)
        goto MLS1Fin;
    if (row1 != 1 || col1 != 1) {
        mat_err(8);
        goto MLS1Fin;
    }
    me = (int)MX[1][1];

    if ((p = n_mexpr(p,1,1,&row1,&col1,&ivflg,NULL)) == NULL)
        goto MLS1Fin;
    if ((p = m_check1(p)) == NULL)
        goto MLS1Fin;
    if (row1 != 1 || col1 != 1) {
        mat_err(8);
        goto MLS1Fin;
    }
    mi = (int)MX[1][1];

    ma = row - me - mi;
    if (me < 0 || mi < 0 || ma < 1) {
        mat_err(14);
        goto MLS1Fin;
    }
    n = col - 1;

    if (!(xx = (double *)calloc(n * n + 1,sizeof(double)))) {
        mat_err(3);
        goto MLS1Fin;
    }                
    xxa = n * n + 1;
    memrq(xxa,sizeof(double));

    if (!(xy = (double *)calloc(n + 1,sizeof(double)))) {
        mat_err(3);
        goto MLS1Fin;
    }                
    xya = n + 1;
    memrq(xya,sizeof(double));

    m = me + mi;

    if (m > 0) {
        if (!(a = (double *)calloc(m * n + 1,sizeof(double)))) {
            mat_err(3);
            goto MLS1Fin;
        }                
        aa = m * n + 1;
        memrq(aa,sizeof(double));

        if (!(b = (double *)calloc(m + 1,sizeof(double)))) {
            mat_err(3);
            goto MLS1Fin;
        }                
        ba = m + 1;
        memrq(ba,sizeof(double));
    }

    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) {
            tmp = 0.0;
            for (k = 1; k <= ma; ++k) {
                tmp += MX[0][(me + k - 1) * col + i] *
                       MX[0][(me + k - 1) * col + j];
            }
            xx[(i - 1) * n + j] = tmp;
        }   
        tmp = 0.0;
        for (k = 1; k <= ma; ++k) 
            tmp -= MX[0][(me + k - 1) * col + i] *
                       MX[0][(me + k - 1) * col + col];
        xy[i] = tmp;
    }
    if (me > 0) {
        for (i = 1; i <= me; ++i) {
            for (j = 1; j <= n; ++j)
                a[(i - 1) * n + j] = MX[0][(i - 1) * col + j];

            b[i] = MX[0][(i - 1) * col + col];
        }
    }
    if (mi > 0) {
        for (i = 1; i <= mi; ++i) {
            for (j = 1; j <= n; ++j)
                a[(me + i - 1) * n + j] = MX[0][(me + ma + i - 1) * col + j];

            b[me + i] = MX[0][(me + ma + i - 1) * col + col];
        }
    }
          
    if ((p = m_getmat(p,n,1,&idx,cmd,1)) == NULL)
        goto MLS1Fin;

    maxit = 50 * (n + m);
    vsmall = 2.0 * EPSI;
        
    /* allocate and initialize xl and xu for lower and upper bounds */

    if (!(xl = (double *)calloc(row + 1,sizeof(double)))) {
        mat_err(3);
        goto MLS1Fin;
    }                
    xla = row + 1;
    memrq(xla,sizeof(double));

    if (!(xu = (double *)calloc(row + 1,sizeof(double)))) {
        mat_err(3);
        goto MLS1Fin;
    }                
    xua = row + 1;
    memrq(xua,sizeof(double));

    for (i = 1; i <= row; ++i) {
        xu[i] = DBLMAX / 1000.0;
        xl[i] = -xu[i];             
    }

    /* allocate iact */

    if (!(iact = (int *)calloc(row + 1,sizeof(int)))) {
        mat_err(3);
        goto MLS1Fin;
    }                
    iaa = row + 1;
    memrq(iaa,sizeof(int));

    r = qld(n,m,me,1,maxit,a,b,xy,xx,MatVal[idx],xl,xu,iact,&nact,vsmall,&diag);

    /* printf("r=%d nact=%d diag=%g\n",r,nact,diag); */

    switch (r) {
        case  0: break;
        case  1: m_cmdmsg();
                 printf1("Exceeded maximal number of iterations (%d).\n",maxit);
                 goto MLS1Fin;
        case  2: m_cmdmsg();
                 printf1("Insufficient accuracy for convergence.\n");
                 goto MLS1Fin;
        case  3: m_cmdmsg();
                 printf1("Insufficient memory.\n");
                 goto MLS1Fin;
        default: m_cmdmsg();
                 printf1("Inconsistent constraints.\n");
                 goto MLS1Fin;
    }
    if (diag > 0.0)   
        printf2("Coeff matrix was enlarged by %g times unitmatrix.\n",diag);
    err = 0;

MLS1Fin:
    if (xla > 0) {
        free((char *)xl);
        memrq(-xla,sizeof(double));
    }
    if (xua > 0) {
        free((char *)xu);
        memrq(-xua,sizeof(double));
    }
    if (iaa > 0) {
        free((char *)iact);
        memrq(-iaa,sizeof(int));
    }
    if (xxa > 0) {
        free((char *)xx);
        memrq(-xxa,sizeof(double));
    }
    if (xya > 0) {
        free((char *)xy);
        memrq(-xya,sizeof(double));
    }
    if (aa > 0) {
        free((char *)a);
        memrq(-aa,sizeof(double));
    }
    if (ba > 0) {
        free((char *)b);
        memrq(-ba,sizeof(double));
    }
    mx_free();
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

int m_kmet(char *cmd)
{
    register int i,j;
    int err,row,col,idx,n,ivflg,a,b;
    register char *p;
    double ai,aj,bi,bj,d;

    err = -1;
    ivflg = 0;

    printf2("%s\n",cmd);
    p = n_mexpr(cmd + 6,1,0,&row,&col,&ivflg,NULL);
    if (p == NULL)
        goto MKMETFin;
    if ((p = m_check1(p)) == NULL)
        goto MKMETFin;
     
    if ((p = m_getmat(p,row,row,&idx,cmd,1)) == NULL)
        goto MKMETFin;

    for (a = 2; a <= row; ++a) {
        for (b = 1; b < a; ++b) {

            /* calculate Kemeny distance between row a and row b */

            d = 0.0;
            for (i = 2; i <= col; ++i) {
                ai = get_mxval(0,a,i);  
                bi = get_mxval(0,b,i);  
                for (j = 1; j < i; ++j) {
                    aj = get_mxval(0,a,j);  
                    bj = get_mxval(0,b,j);  
                    if (fabs(ai - aj) <= EPSI2) {
                        if (fabs(bi - bj) > EPSI2)  
                            d += 1.0;
                    }
                    else if (ai < aj) {
                        if (fabs(bi - bj) <= EPSI2)  
                            d += 1.0;
                        else if (bi > bj)
                            d += 2.0;
                    }
                    else {
                        if (fabs(bi - bj) <= EPSI2)  
                            d += 1.0;
                        else if (bi < bj)
                            d += 2.0;
                    }
                }
            }
            MatVal[idx][(a - 1) * row + b] = d;            
            MatVal[idx][(b - 1) * row + a] = d;            
        }
    }
    err = 0;

MKMETFin:
    mx_free();
    return(err);
}



