/****************************************************************************/
/*  t_matf                                                                  */
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
#include "t_eval.h"
#include "t_eval1.h"
#include "t_eval2.h"
#include "t_mat.h"
#include "t_matc.h"
#include "tda_context.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_matf                                                     */

void mp_info(TDAContext *ctx);                        
int mp_info1(TDAContext *ctx, int n,int t);                 
int mp_alloc(TDAContext *ctx, int t,int m,int n);                 
void mp_putvar(TDAContext *ctx, int i,int n,double *mat,int nv,short *vidx,int icase); 
int mp_putlog(TDAContext *ctx, double x);
int mp_putpar(TDAContext *ctx, int n,double *x);
int mp_putmpar(TDAContext *ctx, int m,int n,double *x);
int mp_putcov(TDAContext *ctx, int n,double *x);
int get_mexpr(TDAContext *ctx, char *exp,int *row,int *col,int *ivflg);
int eval_mexpr(TDAContext *ctx, int off,int row,int col,double *mat,int ivflg, double *mat1,int *rsel,int *csel,int idx,double dval);

/* ------------------------------------------------------------------------ */
/*  global variables.                                                       */


/* ------------------------------------------------------------------------ */
/*  mp_info()   Print info about matrices with estimation results.          */

void mp_info(TDAContext *ctx)                         
{
    int n = 0;
  
    n += mp_info1(ctx, ctx->PMMPLogDef,1);
    n += mp_info1(ctx, ctx->PMMPParDef,2);
    n += mp_info1(ctx, ctx->PMMPCovDef,3);
    n += mp_info1(ctx, ctx->PMMPGradDef,4);
    n += mp_info1(ctx, ctx->PMMPResDef,5);
    if (n) 
        newline(ctx);
}

/* ------------------------------------------------------------------------ */
/*  mp_info1(n,t)                                                           */

int mp_info1(TDAContext *ctx, int n,int t)                  
{
    int idx = 0;

    if (n == 0)
        return(0);
    if (n > 0)
        printf1(ctx, "Created");             
    else if (n < 0)
        printf1(ctx, "Could not create");
    printf1(ctx, " matrix (mp");
    switch (t) {
        case 1: printf1(ctx, "log ): "); idx = ctx->MPLogIdx; break;
        case 2: printf1(ctx, "par ): "); idx = ctx->MPParIdx; break;
        case 3: printf1(ctx, "cov ): "); idx = ctx->MPCovIdx; break;
        case 4: printf1(ctx, "grad): "); idx = ctx->MPGradIdx; break;
        case 5: printf1(ctx, "res ): "); idx = ctx->MPResIdx; break;
    }
    printf1(ctx, "%s", ctx->MatName[idx]);          
    /* prnchar(' ',VNLMax - strlen(MatName[idx]),0); */
    printf1(ctx, " [%d,%d]",ctx->MatRow[idx],ctx->MatCol[idx]);          
    if ((t == 4 && ctx->PM2NV > 0) || (t == 5 && ctx->PM3NV > 0))
        printf1(ctx, " %s",ctx->MatDef[idx]);          
    newline(ctx);
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  mp_alloc(t,m,n)         allocate special matrices                       */
/*                          t = 1 : MPLog                                   */
/*                          t = 2 : MPPar                                   */
/*                          t = 3 : MPCov                                   */
/*                          t = 4 : MPGrad                                  */
/*                          t = 5 : MPRes                                   */
/*  return 0 if OK, -1 if error.                                            */

int mp_alloc(TDAContext *ctx, int t,int m,int n)                  
{
    register int i,l;

    switch (t) {
        case 1: 
            if ((ctx->MPLogIdx = mat_newmat(ctx, ctx->PMMPLogName,m,n)) < 0) {
                ctx->MPLogIdx = -1;
                return(-1);
            }
            strcpy(ctx->MatDef[ctx->MPLogIdx],"mplog");
            break;

        case 2:
            if ((ctx->MPParIdx = mat_newmat(ctx, ctx->PMMPParName,m,n)) < 0) {
                ctx->MPParIdx = -1;
                return(-1);
            }
            strcpy(ctx->MatDef[ctx->MPParIdx],"mppar");
            break;

        case 3:
            if ((ctx->MPCovIdx = mat_newmat(ctx, ctx->PMMPCovName,m,n)) < 0) {
                ctx->MPCovIdx = -1;
                return(-1);
            }
            strcpy(ctx->MatDef[ctx->MPCovIdx],"mpcov");
            break;

        case 4:
            if ((ctx->MPGradIdx = mat_newmat(ctx, ctx->PMMPGradName,m,n)) < 0) {
                ctx->MPGradIdx = -1;
                ctx->MPGradRow = ctx->MPGradCol = 0;
                return(-1);
            }
            ctx->MPGradRow = m;
            ctx->MPGradCol = n;
            strcpy(ctx->MatDef[ctx->MPGradIdx],"mpgrad");
            l = 8;
            if (ctx->PM2NV > 0) {
                for (i = 0; i < ctx->PM2NV; ++i) {
                    if (l + (int)strlen(ctx->VName[ctx->PM2VIdx[i]]) + 1 < MatDefLen) {
                        if (i == 0)
                            strcat(ctx->MatDef[ctx->MPGradIdx],"(");
                        else       
                            strcat(ctx->MatDef[ctx->MPGradIdx],",");
                        strcat(ctx->MatDef[ctx->MPGradIdx],ctx->VName[ctx->PM2VIdx[i]]);
                    }
                }
                strcat(ctx->MatDef[ctx->MPGradIdx],")");
            }
            break;

        case 5:
            if ((ctx->MPResIdx = mat_newmat(ctx, ctx->PMMPResName,m,n)) < 0) {
                ctx->MPResIdx = -1;
                ctx->MPResRow = ctx->MPResCol = 0;
                return(-1);
            }
            ctx->MPResRow = m;
            ctx->MPResCol = n;
            strcpy(ctx->MatDef[ctx->MPResIdx],"mpres");
            l = 7;
            if (ctx->PM2NV > 0) {
                for (i = 0; i < ctx->PM3NV; ++i) {
                    if (l + (int)strlen(ctx->VName[ctx->PM3VIdx[i]]) + 1 < MatDefLen) {
                        if (i == 0)
                            strcat(ctx->MatDef[ctx->MPGradIdx],"(");
                        else       
                            strcat(ctx->MatDef[ctx->MPGradIdx],",");
                        strcat(ctx->MatDef[ctx->MPGradIdx],ctx->VName[ctx->PM3VIdx[i]]);
                    }
                }
                strcat(ctx->MatDef[ctx->MPGradIdx],")");
            }
            break;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mp_putvar(i,n,mat,nv,vidx,icase)                                        */
/*                                                                          */
/*  add values of variables to special matrices MPGrad, MPRes               */

void mp_putvar(TDAContext *ctx, int i,int n,double *mat,int nv,short *vidx,int icase)
{
    register int j,k;

    j = n - nv + 1;
    if (j < 1) {
        printfe(ctx, "ERROR in mp_putvar (n=%d, nv=%d)\n",n,nv);
        gerr_exit(ctx, 200);
    }
    for (k = 0; k < nv; ++k) {
        mat[i * n + j] = get_data(ctx, (int)vidx[k],icase);
        j++;
    }
}

/* ------------------------------------------------------------------------ */
/*  mp_putlog(x)        Create PMMPLog matrix and put x                     */
/*                      into the matrix.                                    */
/*  Return 0 if OK, -1 if error.                                            */

int mp_putlog(TDAContext *ctx, double x)
{
    if (ctx->PMMPLogDef == 1) {      
        if (mp_alloc(ctx, 1,1,1)) {          
            ctx->PMMPLogDef = -1;
            return(-1);
        }
        else  
            ctx->MatVal[ctx->MPLogIdx][1] = x;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mp_putpar(n,x)      Create PMMPPar matrix and put x[i] (i=1,n)          */
/*                      into the matrix.                                    */
/*  Return 0 if OK, -1 if error.                                            */

int mp_putpar(TDAContext *ctx, int n,double *x)
{
    register int i;

    if (n < 1)
        return(0);

    if (ctx->PMMPParDef == 1) {      
        if (mp_alloc(ctx, 2,n,1)) {          
            ctx->PMMPParDef = -1;
            return(-1);
        }
        else {
            for (i = 1; i <= n; ++i)   
                ctx->MatVal[ctx->MPParIdx][i] = x[i];
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mp_putmpar(m,n,x)   Create PMMPPar matrix with dimension (m,n) and      */
/*                      copy x[] into the matrix.                           */
/*  Return 0 if OK, -1 if error.                                            */

int mp_putmpar(TDAContext *ctx, int m,int n,double *x)
{
    register int i,mn;

    mn = m * n;
    if (mn < 1)
        return(0);

    if (ctx->PMMPParDef == 1) {      
        if (mp_alloc(ctx, 2,m,n)) {          
            ctx->PMMPParDef = -1;
            return(-1);
        }
        else {
            for (i = 1; i <= mn; ++i)   
                ctx->MatVal[ctx->MPParIdx][i] = x[i];
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mp_putcov(n,x)      Create PMMPCov matrix and put cov matrix into       */
/*                      the matrix.                                         */
/*  Return 0 if OK, -1 if error.                                            */

int mp_putcov(TDAContext *ctx, int n,double *x)
{
    register int i,j;

    if (n < 1)
        return(0);

    if (ctx->PMMPCovDef == 1) {      
        if (mp_alloc(ctx, 3,n,n)) {          
            ctx->PMMPCovDef = -1;
            return(-1);
        }
        else {
            for (i = 0; i < n; ++i) { 
                for (j = 1; j <= n; ++j)   
                    ctx->MatVal[ctx->MPCovIdx][i * n + j] = x[i * n + j];
            }
        }
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  get_mexpr(expr,row,col,ivflg)                                           */
/*                                                                          */
/*  Undocumented in the manual: an interval                                      */
/*  literal [lo,hi] may appear anywhere a scalar constant may.  It turns   */
/*  the whole expression interval-valued (*ivflg = 1): every element is    */
/*  then a pair, lower bounds in MX[], upper bounds in the second value    */
/*  plane MX1[] (allocated by mx1_alloc), evaluated with the i_* interval  */
/*  arithmetic.  mpr prints each element as its lo,hi pair.  Control       */
/*  case: examples/coverage/mivops.cf.                                     */
/*                                                                          */
/*  Get matrix expression.                                                  */
/*  Get parser stack and return also row and col dimension.                 */
/*  If ivflg != 0 allow for intervals. Return: ivflg = 1 if expression      */
/*  contains intervals, otherwise ivflg = 0.                                */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int get_mexpr(TDAContext *ctx, char *exp,int *row,int *col,int *ivflg)
{
    int err,idx;
       
    err = -1;

    ctx->MatExprFlg = 1;
    *row = *col = 0;

    if ((idx = v_parse(ctx, exp,*ivflg)) < 0 || ctx->ESCnt <= 0) {
        m_cmdmsg(ctx);
        printf1(ctx, "Syntax error (%d) in expression.\n",idx);
        if (idx < 0)
            prn_emsg1(ctx, idx);
        goto GMEXPRFin;
    }
    *ivflg = ctx->IVEXPRFlg;
/**
    printf1(ctx, "row=%d col=%d\n",MatExprRow,MatExprCol);
    for (i = 0; i < ESCnt; ++i) {
        printf1(ctx, "i=%d ESTyp=%d  ESVal=%lg ESIdx=%d \n",i,ESTyp[i],ESVal[i],ESIdx[i]);
    }
**/

    *row = ctx->MatExprRow;
    *col = ctx->MatExprCol;
    err = 0;

GMEXPRFin:
    ctx->MatExprFlg = 0;
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  eval_mexpr(off,row,col,mat,ivflg,mat1,rsel,csel,sidx,dval)              */
/*  ##                                                                      */  
/*  Evaluate matrix expression and save values in mat[].                    */
/*  for rows: i = off,...,row-1                                             */
/*                                                                          */
/*  If ivflg = 1 allow for intervals and save lower bounds in mat[],        */
/*  upper bounds in mat1[].                                                 */
/*                                                                          */
/*  If rsel != NULL and/or csel != NULL take this as selecting rows and     */
/*  columns for evaluating the matrix expression. If sidx >= 0 take this    */
/*  as index to a matrix to provide values for not selected elements. If    */
/*  sidx < 0 use dval.                                                      */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int eval_mexpr(TDAContext *ctx, int off,int row,int col,double *mat,int ivflg, double *mat1,int *rsel,int *csel,int sidx,double dval)
{
    register int i,j,k;
    int idx,r,err,kt,opt,nrow,isel,jsel;
    double tmp = 0.0,res[2];

    err = -1;
    nrow = row - off;

    if (alloc_aci(ctx, ctx->ESCnt))
        return(-1);

    if (alloc_acn(ctx, ctx->ESCnt))
        return(-1);
             
    if (alloc_acm(ctx, ctx->ESCnt))
        return(-1);
             
    if (alloc_acns(ctx, ctx->ESCnt))
        return(-1);      

    opt = 0;
    for (k = 0; k < ctx->ESCnt; ++k) {
        ctx->AcM[k] = ctx->ESTyp[k];                          /* need be saved */

        ctx->AcI[k] = 0;
        kt = ctx->ESTyp[k];
        if (kt < 0) {
            ctx->AcNS[k] = 1;
            kt = -kt;
        }
        if (ctx->OPT2A <= kt && kt < ctx->OPT2B)              /* contains typ 2 ops */
            opt = 1;

        if (kt >= ctx->MOFFS && kt < ctx->VOFFS) {
            ctx->AcI[k] = 1;
            ctx->AcN[k] = kt - ctx->MOFFS;
        }
        else if (kt >= ctx->NLOFFS && kt < ctx->MOFFS) {
            ctx->AcI[k] = 2;
            ctx->AcN[k] = kt - ctx->NLOFFS;
        }
        else if (kt >= ctx->VOFFS && kt < ctx->COFFS) {
            ctx->AcI[k] = 3;
            ctx->AcN[k] = kt - ctx->VOFFS;
        }
    }
    if (opt == 0) {                 /* without type 2 operators */

        for (k = 0; k < ctx->ESCnt; ++k) {
            if (ctx->AcI[k])  
                ctx->ESTyp[k] = 0;
        }
        for (i = off; i < row; ++i) {
            for (j = 0; j < col; ++j) {
                for (k = 0; k < ctx->ESCnt; ++k) {
                    if (ctx->AcI[k] == 0)
                        continue;
                    idx = ctx->AcN[k];
                    if (ctx->AcI[k] == 1)  
                        tmp = ctx->MatVal[idx][(i % ctx->MatRow[idx]) * ctx->MatCol[idx] + j % ctx->MatCol[idx] + 1];
                    else if (ctx->AcI[k] == 2)
                        tmp = get_data(ctx, ctx->NLVIdx[idx][j % ctx->NLNV[idx]],i % ctx->NOC);
                    else if (ctx->AcI[k] == 3)
                        tmp = get_data(ctx, idx,i % ctx->NOC);

                    if (ctx->AcNS[k])
                        tmp = -tmp;
                    ctx->ESVal[k] = tmp;
                }

                /*** 
                printf1(ctx, "iii=%d\n",i);   
                for (k = 0; k < ESCnt; ++k)  
                    printf1(ctx, "k=%d ESTyp=%d  ESVal=%lg ESIdx=%d\n",k,ESTyp[k],ESVal[k],ESIdx[k]);
                ***/

                isel = jsel = 1;
                if (rsel != NULL && rsel[i] == 0)
                    isel = 0;
                if (csel != NULL && csel[j] == 0)
                    jsel = 0;

                if (isel != 0 && jsel != 0) {
                    r = v_eval1(ctx, 0,ctx->ESCnt,ctx->ESTyp,ctx->ESVal,ctx->ESIdx,res,0,0,0,0,ivflg);
                    if (r) {
                        m_cmdmsg(ctx);
                        printf1(ctx, "Can't evaluate expression for matrix element (%d,%d).\n",i+1,j+1);
                        prn_emsg2(ctx, r);
                        goto EMEXPRFin;
                    }
                    mat[i * col + j + 1] = res[0];
                    if (ivflg)
                        mat1[i * col + j + 1] = res[1];
                }
                else {
                    if (sidx >= 0)  
                        dval = get_matval(ctx, sidx,i + 1,j + 1);
 
                    mat[i * col + j + 1] = dval;
                    if (ivflg)
                        mat1[i * col + j + 1] = dval;
                }
            }
        }
    }
    else {              /* with type 2 operators */

        if (alloc_acx(ctx, nrow + 1))
            goto EMEXPRFin;

        for (j = 0; j < col; ++j) {

            for (k = 0; k < ctx->ESCnt; ++k) {
                if (ctx->AcI[k] == 2) {     /* select variable from namelist */
                    ctx->ESTyp[k] = ctx->NLVIdx[ctx->AcN[k]][j % ctx->NLNV[ctx->AcN[k]]] + ctx->VOFFS;
                    if (ctx->AcNS[k])
                        ctx->ESTyp[k] = -ctx->ESTyp[k];
                }
            }
     
            /*** 
            printf1(ctx, "col=%d\n",j + 1);   
            for (k = 0; k < ESCnt; ++k)  
                printf1(ctx, "II: k=%d ESTyp=%d  ESVal=%lg ESIdx=%d AcN=%d\n",
                           k,ESTyp[k],ESVal[k],ESIdx[k],AcN[k]);
            ***/

            r = v_eval2(ctx, -1,ctx->ESCnt,ctx->ESTyp,ctx->ESVal,ctx->ESIdx,ctx->AcX,off,row,j,0);
            if (r) {
                m_cmdmsg(ctx);
                printf1(ctx, "Can't evaluate expression for matrix column %d.\n",j+1);
                prn_emsg2(ctx, r);
                goto EMEXPRFin;
            }
            for (k = off; k < row; ++k)
                mat[k * col + j + 1] = ctx->AcX[k - off];
        }
    }
    for (k = 0; k < ctx->ESCnt; ++k)  
        ctx->ESTyp[k] = ctx->AcM[k];

    err = 0;

EMEXPRFin:
    alloc_acns(ctx, 0);
    alloc_acm(ctx, 0);
    alloc_acn(ctx, 0);
    alloc_aci(ctx, 0);
    alloc_acx(ctx, 0);
    return(err);
}
