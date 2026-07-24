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

/* ------------------------------------------------------------------------ */
/*  functions in t_matf                                                     */

void mp_info(void);                        
int mp_info1(int n,int t);                 
int mp_alloc(int t,int m,int n);                 
void mp_putvar(int i,int n,double *mat,int nv,short *vidx,int icase); 
int mp_putlog(double x);
int mp_putpar(int n,double *x);
int mp_putmpar(int m,int n,double *x);
int mp_putcov(int n,double *x);
int get_mexpr(char *exp,int *row,int *col,int *ivflg);
int eval_mexpr(int off,int row,int col,double *mat,int ivflg,
    double *mat1,int *rsel,int *csel,int idx,double dval);

/* ------------------------------------------------------------------------ */
/*  global variables.                                                       */

int MPLogIdx = -1;          /* index of PMMPLog  matrix                     */
int MPParIdx = -1;          /* index of PMMPPar  matrix                     */
int MPCovIdx = -1;          /* index of PMMPCov  matrix                     */
int MPGradIdx = -1;         /* index of PMMPGrad matrix                     */
int MPResIdx = -1;          /* index of PMMPRes  matrix                     */
int MPGradRow = 0;          /* allocated rows in PMMPGrad                   */
int MPGradCol = 0;          /* allocated columns in PMMPGrad                */
int MPResRow = 0;           /* allocated rows in PMMPRes                    */
int MPResCol = 0;           /* allocated columns in PMMPRes                 */

/* ------------------------------------------------------------------------ */
/*  mp_info()   Print info about matrices with estimation results.          */

void mp_info(void)                         
{
    int n = 0;
  
    n += mp_info1(PMMPLogDef,1);
    n += mp_info1(PMMPParDef,2);
    n += mp_info1(PMMPCovDef,3);
    n += mp_info1(PMMPGradDef,4);
    n += mp_info1(PMMPResDef,5);
    if (n) 
        newline();
}

/* ------------------------------------------------------------------------ */
/*  mp_info1(n,t)                                                           */

int mp_info1(int n,int t)                  
{
    int idx;

    if (n == 0)
        return(0);
    if (n > 0)
        printf1("Created");             
    else if (n < 0)
        printf1("Could not create");
    printf1(" matrix (mp");
    switch (t) {
        case 1: printf1("log ): "); idx = MPLogIdx; break;
        case 2: printf1("par ): "); idx = MPParIdx; break;
        case 3: printf1("cov ): "); idx = MPCovIdx; break;
        case 4: printf1("grad): "); idx = MPGradIdx; break;
        case 5: printf1("res ): "); idx = MPResIdx; break;
    }
    printf1(MatName[idx]);          
    /* prnchar(' ',VNLMax - strlen(MatName[idx]),0); */
    printf1(" [%d,%d]",MatRow[idx],MatCol[idx]);          
    if ((t == 4 && PM2NV > 0) || (t == 5 && PM3NV > 0))
        printf1(" %s",MatDef[idx]);          
    newline();
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

int mp_alloc(int t,int m,int n)                  
{
    register int i,l;

    switch (t) {
        case 1: 
            if ((MPLogIdx = mat_newmat(PMMPLogName,m,n)) < 0) {
                MPLogIdx = -1;
                return(-1);
            }
            strcpy(MatDef[MPLogIdx],"mplog");
            break;

        case 2:
            if ((MPParIdx = mat_newmat(PMMPParName,m,n)) < 0) {
                MPParIdx = -1;
                return(-1);
            }
            strcpy(MatDef[MPParIdx],"mppar");
            break;

        case 3:
            if ((MPCovIdx = mat_newmat(PMMPCovName,m,n)) < 0) {
                MPCovIdx = -1;
                return(-1);
            }
            strcpy(MatDef[MPCovIdx],"mpcov");
            break;

        case 4:
            if ((MPGradIdx = mat_newmat(PMMPGradName,m,n)) < 0) {
                MPGradIdx = -1;
                MPGradRow = MPGradCol = 0;
                return(-1);
            }
            MPGradRow = m;
            MPGradCol = n;
            strcpy(MatDef[MPGradIdx],"mpgrad");
            l = 8;
            if (PM2NV > 0) {
                for (i = 0; i < PM2NV; ++i) {
                    if (l + strlen(VName[PM2VIdx[i]]) + 1 < MatDefLen) {
                        if (i == 0)
                            strcat(MatDef[MPGradIdx],"(");
                        else       
                            strcat(MatDef[MPGradIdx],",");
                        strcat(MatDef[MPGradIdx],VName[PM2VIdx[i]]);
                    }
                }
                strcat(MatDef[MPGradIdx],")");
            }
            break;

        case 5:
            if ((MPResIdx = mat_newmat(PMMPResName,m,n)) < 0) {
                MPResIdx = -1;
                MPResRow = MPResCol = 0;
                return(-1);
            }
            MPResRow = m;
            MPResCol = n;
            strcpy(MatDef[MPResIdx],"mpres");
            l = 7;
            if (PM2NV > 0) {
                for (i = 0; i < PM3NV; ++i) {
                    if (l + strlen(VName[PM3VIdx[i]]) + 1 < MatDefLen) {
                        if (i == 0)
                            strcat(MatDef[MPGradIdx],"(");
                        else       
                            strcat(MatDef[MPGradIdx],",");
                        strcat(MatDef[MPGradIdx],VName[PM3VIdx[i]]);
                    }
                }
                strcat(MatDef[MPGradIdx],")");
            }
            break;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mp_putvar(i,n,mat,nv,vidx,icase)                                        */
/*                                                                          */
/*  add values of variables to special matrices MPGrad, MPRes               */

void mp_putvar(int i,int n,double *mat,int nv,short *vidx,int icase)
{
    register int j,k;

    j = n - nv + 1;
    if (j < 1) {
        printfe("ERROR in mp_putvar (n=%d, nv=%d)\n",n,nv);
        gerr_exit(200);
    }
    for (k = 0; k < nv; ++k) {
        mat[i * n + j] = get_data((int)vidx[k],icase);
        j++;
    }
}

/* ------------------------------------------------------------------------ */
/*  mp_putlog(x)        Create PMMPLog matrix and put x                     */
/*                      into the matrix.                                    */
/*  Return 0 if OK, -1 if error.                                            */

int mp_putlog(double x)
{
    if (PMMPLogDef == 1) {      
        if (mp_alloc(1,1,1)) {          
            PMMPLogDef = -1;
            return(-1);
        }
        else  
            MatVal[MPLogIdx][1] = x;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mp_putpar(n,x)      Create PMMPPar matrix and put x[i] (i=1,n)          */
/*                      into the matrix.                                    */
/*  Return 0 if OK, -1 if error.                                            */

int mp_putpar(int n,double *x)
{
    register int i;

    if (n < 1)
        return(0);

    if (PMMPParDef == 1) {      
        if (mp_alloc(2,n,1)) {          
            PMMPParDef = -1;
            return(-1);
        }
        else {
            for (i = 1; i <= n; ++i)   
                MatVal[MPParIdx][i] = x[i];
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mp_putmpar(m,n,x)   Create PMMPPar matrix with dimension (m,n) and      */
/*                      copy x[] into the matrix.                           */
/*  Return 0 if OK, -1 if error.                                            */

int mp_putmpar(int m,int n,double *x)
{
    register int i,mn;

    mn = m * n;
    if (mn < 1)
        return(0);

    if (PMMPParDef == 1) {      
        if (mp_alloc(2,m,n)) {          
            PMMPParDef = -1;
            return(-1);
        }
        else {
            for (i = 1; i <= mn; ++i)   
                MatVal[MPParIdx][i] = x[i];
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mp_putcov(n,x)      Create PMMPCov matrix and put cov matrix into       */
/*                      the matrix.                                         */
/*  Return 0 if OK, -1 if error.                                            */

int mp_putcov(int n,double *x)
{
    register int i,j;

    if (n < 1)
        return(0);

    if (PMMPCovDef == 1) {      
        if (mp_alloc(3,n,n)) {          
            PMMPCovDef = -1;
            return(-1);
        }
        else {
            for (i = 0; i < n; ++i) { 
                for (j = 1; j <= n; ++j)   
                    MatVal[MPCovIdx][i * n + j] = x[i * n + j];
            }
        }
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  get_mexpr(expr,row,col,ivflg)                                           */
/*                                                                          */
/*  Get matrix expression.                                                  */
/*  Get parser stack and return also row and col dimension.                 */
/*  If ivflg != 0 allow for intervals. Return: ivflg = 1 if expression      */
/*  contains intervals, otherwise ivflg = 0.                                */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int get_mexpr(char *exp,int *row,int *col,int *ivflg)
{
    int err,idx;
       
    err = -1;

    MatExprFlg = 1;
    *row = *col = 0;

    if ((idx = v_parse(exp,*ivflg)) < 0 || ESCnt <= 0) {
        m_cmdmsg();
        printf1("Syntax error (%d) in expression.\n",idx);
        if (idx < 0)
            prn_emsg1(idx);
        goto GMEXPRFin;
    }
    *ivflg = IVEXPRFlg;
/**
    printf1("row=%d col=%d\n",MatExprRow,MatExprCol);
    for (i = 0; i < ESCnt; ++i) {
        printf1("i=%d ESTyp=%d  ESVal=%lg ESIdx=%d \n",i,ESTyp[i],ESVal[i],ESIdx[i]);
    }
**/

    *row = MatExprRow;
    *col = MatExprCol;
    err = 0;

GMEXPRFin:
    MatExprFlg = 0;
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

int eval_mexpr(int off,int row,int col,double *mat,int ivflg,
    double *mat1,int *rsel,int *csel,int sidx,double dval)
{
    register int i,j,k;
    int n,idx,r,err,kt,opt,nrow,isel,jsel;
    double tmp,res[2];

    err = -1;
    nrow = row - off;
    n = row * col;

    if (alloc_aci(ESCnt))
        return(-1);

    if (alloc_acn(ESCnt))
        return(-1);
             
    if (alloc_acm(ESCnt))
        return(-1);
             
    if (alloc_acns(ESCnt))
        return(-1);      

    opt = 0;
    for (k = 0; k < ESCnt; ++k) {
        AcM[k] = ESTyp[k];                          /* need be saved */

        AcI[k] = 0;
        kt = ESTyp[k];
        if (kt < 0) {
            AcNS[k] = 1;
            kt = -kt;
        }
        if (OPT2A <= kt && kt < OPT2B)              /* contains typ 2 ops */
            opt = 1;

        if (kt >= MOFFS && kt < VOFFS) {
            AcI[k] = 1;
            AcN[k] = kt - MOFFS;
        }
        else if (kt >= NLOFFS && kt < MOFFS) {
            AcI[k] = 2;
            AcN[k] = kt - NLOFFS;
        }
        else if (kt >= VOFFS && kt < COFFS) {
            AcI[k] = 3;
            AcN[k] = kt - VOFFS;
        }
    }
    if (opt == 0) {                 /* without type 2 operators */

        for (k = 0; k < ESCnt; ++k) {
            if (AcI[k])  
                ESTyp[k] = 0;
        }
        for (i = off; i < row; ++i) {
            for (j = 0; j < col; ++j) {
                for (k = 0; k < ESCnt; ++k) {
                    if (AcI[k] == 0)
                        continue;
                    idx = AcN[k];
                    if (AcI[k] == 1)  
                        tmp = MatVal[idx][(i % MatRow[idx]) * MatCol[idx] + j % MatCol[idx] + 1];
                    else if (AcI[k] == 2)
                        tmp = get_data(NLVIdx[idx][j % NLNV[idx]],i % NOC);
                    else if (AcI[k] == 3)
                        tmp = get_data(idx,i % NOC);

                    if (AcNS[k])
                        tmp = -tmp;
                    ESVal[k] = tmp;
                }

                /*** 
                printf1("iii=%d\n",i);   
                for (k = 0; k < ESCnt; ++k)  
                    printf1("k=%d ESTyp=%d  ESVal=%lg ESIdx=%d\n",k,ESTyp[k],ESVal[k],ESIdx[k]);
                ***/

                isel = jsel = 1;
                if (rsel != NULL && rsel[i] == 0)
                    isel = 0;
                if (csel != NULL && csel[j] == 0)
                    jsel = 0;

                if (isel != 0 && jsel != 0) {
                    r = v_eval1(0,ESCnt,ESTyp,ESVal,ESIdx,res,0,0,0,0,ivflg);
                    if (r) {
                        m_cmdmsg();
                        printf1("Can't evaluate expression for matrix element (%d,%d).\n",i+1,j+1);
                        prn_emsg2(r);
                        goto EMEXPRFin;
                    }
                    mat[i * col + j + 1] = res[0];
                    if (ivflg)
                        mat1[i * col + j + 1] = res[1];
                }
                else {
                    if (sidx >= 0)  
                        dval = get_matval(sidx,i + 1,j + 1);
 
                    mat[i * col + j + 1] = dval;
                    if (ivflg)
                        mat1[i * col + j + 1] = dval;
                }
            }
        }
    }
    else {              /* with type 2 operators */

        if (alloc_acx(nrow + 1))
            goto EMEXPRFin;

        for (j = 0; j < col; ++j) {

            for (k = 0; k < ESCnt; ++k) {
                if (AcI[k] == 2) {     /* select variable from namelist */
                    ESTyp[k] = NLVIdx[AcN[k]][j % NLNV[AcN[k]]] + VOFFS;
                    if (AcNS[k])
                        ESTyp[k] = -ESTyp[k];
                }
            }
     
            /*** 
            printf1("col=%d\n",j + 1);   
            for (k = 0; k < ESCnt; ++k)  
                printf1("II: k=%d ESTyp=%d  ESVal=%lg ESIdx=%d AcN=%d\n",
                           k,ESTyp[k],ESVal[k],ESIdx[k],AcN[k]);
            ***/

            r = v_eval2(-1,ESCnt,ESTyp,ESVal,ESIdx,AcX,off,row,j,0);
            if (r) {
                m_cmdmsg();
                printf1("Can't evaluate expression for matrix column %d.\n",j+1);
                prn_emsg2(r);
                goto EMEXPRFin;
            }
            for (k = off; k < row; ++k)
                mat[k * col + j + 1] = AcX[k - off];
        }
    }
    for (k = 0; k < ESCnt; ++k)  
        ESTyp[k] = AcM[k];

    err = 0;

EMEXPRFin:
    alloc_acns(0);
    alloc_acm(0);
    alloc_acn(0);
    alloc_aci(0);
    alloc_acx(0);
    return(err);
}
