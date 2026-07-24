/****************************************************************************/
/*  t_eval4                                                                 */
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
#include "t_var.h"
#include "t_gdat.h"
#include "t_rand.h"
#include "t_gf.h"
#include "t_parm.h"
#include "t_alloc.h"
#include "t_gmin.h"
#include "t_int.h"
#include "t_mat.h"
#include "t_matf.h"
#include "t_gcmd.h"
#include "t_eval.h"
#include "t_svd.h"
#include "tda_context.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_eval4.c                                                   */

int alloc_mex(TDAContext *ctx, int n,int opt);
int mex_parse(TDAContext *ctx, char *s);
void mex_expr(TDAContext *ctx);
void mex_expr1(TDAContext *ctx);
void mex_term(TDAContext *ctx);
void mex_factor(TDAContext *ctx);
void mex_match(TDAContext *ctx, int t);
int mex_lexan(TDAContext *ctx);             
int mparse(TDAContext *ctx);          
void m_parse(TDAContext *ctx, short cnt,int *typ,double *val);
int m_eval(TDAContext *ctx, int opt,int cnt,int *typ,double *val);                          
int check_dim_add(TDAContext *ctx, int i,int j);
int check_dim_mul(TDAContext *ctx, int i,int j);
int check_scalar(TDAContext *ctx, int i);
int check_square(TDAContext *ctx, int i);
int alloc_mex_eval(TDAContext *ctx, int opt,int n,int m);
void mex_error(TDAContext *ctx, int n);
int mex_eval(TDAContext *ctx, char *s);


/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

#define NUM 257         /* return code for floating point numbers           */
#define IDV 258         /* return code for V, C or X variables              */
#define IDM 259         /* return code for matrices                         */
#define IDG 260         /* return code for operators                        */
#define IDC 261         /* return code for constants, operators without     */


#define MEXOPMax 100

/*  The stack is organized as follows:                                      */
/*  MEXTyp =  0  floating point number (value in MEXVal)                    */
/*            1  +                                                          */
/*            2  -                                                          */
/*            3  *                                                          */
/*            5  /                                                          */
/*            6  ^                                                          */

/*  array with keywords for operators                                       */

struct MOPKEY {
    char *op;           /* keyword for operator                             */
    int  typ;           /* type number of operator                          */
    int  len;           /* length of keyword                                */
    int  narg;          /* number of arguments (0 if variable)              */
} MOPS[] = {
    { "cross",  101,  5,  1 },
    { "eexp",   105,  4,  1 },
    { "ginv",   100,  4,  1 },
    { "sqrt",   104,  4,  1 },
    { "trace",  103,  5,  1 },
    { "t",      102,  1,  1 },
    { "",         0,  0,  0 },
};


/* ------------------------------------------------------------------------ */
/*  alloc_mex(n,opt)                                                        */
/*                                                                          */
/*  Allocate (opt = 1) or free (opt = 0) MEXTyp with n elements.            */
/*  Return 0 if OK, -1 if error.                                            */

int alloc_mex(TDAContext *ctx, int n,int opt)
{
    if (ctx->MEXMaxLen > 0) {
        free((char *)ctx->MEXTyp);
        free((char *)ctx->MEXVal);
        memrq(ctx, -ctx->MEXMaxLen - 2,sizeof(int) + sizeof(double));
        ctx->MEXMaxLen = 0;
    }
    if (opt == 1) {
        if (!(ctx->MEXTyp = (int *)calloc((size_t)(n + 2),sizeof(int)))) { 
            p_err(ctx, -2,1);  
            return(-1);
        }
        if (!(ctx->MEXVal = (double *)calloc((size_t)(n + 2),sizeof(double)))) { 
            free((char *)ctx->MEXTyp);
            p_err(ctx, -2,1);  
            return(-1);
        }
        ctx->MEXMaxLen = n;
        memrq(ctx, n + 2,sizeof(int) + sizeof(double));
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  r = mex_parse(string)                                                   */
/*                                                                          */
/*  Transforms string to postfix notation. Results are put in the arrays    */
/*  MEXTyp[i]. Max stack size is assumed to be MEXMaxLen. After             */
/*  parsing, the number of entries in the stack is put into MEXCnt.         */
/*                                                                          */
/*  Return: 0 if ok, otherwise MEXErr.                                      */
/*                                                                          */
/*  MEXErr = 0  no error                                                    */
/*          -1  stack overflow                                              */
/*          -2  unknown operator or expression                              */
/*          -3  syntax error                                                */
/*          -6  scan error for floating point numbers                       */
/*          -7  misplaced minus sign                                        */
/*         -10  undefined name of matrix                                    */  
/*         -25  misplaced brackets                                          */
/*         -99  fatal error in stack organization                           */

int mex_parse(TDAContext *ctx, char *s)
{
    ctx->MEXLMin = ctx->MEXErr = ctx->MEXCnt = 0;
    ctx->MEXESP = s;
    ctx->MEXLLast = 1;
    ctx->MEX0Last = 0;

    ctx->MEXNxt = mex_lexan(ctx);
    if (ctx->MEXErr)
        return(ctx->MEXErr);
    while (ctx->MEXNxt) {
        mex_expr(ctx);
        if (ctx->MEXErr)
            return(ctx->MEXErr);
    }
    return(ctx->MEXErr);
}

/* ------------------------------------------------------------------------ */
/*  mex_expr                                                                */

void mex_expr(TDAContext *ctx)
{
    if (ctx->MEXCnt >= ctx->MEXMaxLen) {
        ctx->MEXErr = -1;
        return;
    }
    mex_expr1(ctx);
    if (ctx->MEXErr)
        return;

    while (1) {
        if (ctx->MEXCnt >= ctx->MEXMaxLen) {
            ctx->MEXErr = -1;
            break;  
        }
        switch (ctx->MEXNxt) {
            case '+':   mex_match(ctx, ctx->MEXNxt);
                        if (ctx->MEXErr)
                            return;
                        mex_expr1(ctx);
                        if (ctx->MEXErr)
                            return;
                        ctx->MEXTyp[ctx->MEXCnt] = 1;
                        ctx->MEXVal[ctx->MEXCnt] = 2.0;              
                        ctx->MEXCnt++;
                        break;
            case '-':   mex_match(ctx, ctx->MEXNxt);
                        if (ctx->MEXErr)
                            return;
                        mex_expr1(ctx);
                        if (ctx->MEXErr)
                            return;
                        ctx->MEXTyp[ctx->MEXCnt] = 2;
                        ctx->MEXVal[ctx->MEXCnt] = 2.0;              
                        ctx->MEXCnt++;
                        break;

            default:    return;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  mex_expr1                                                               */

void mex_expr1(TDAContext *ctx)
{
    if (ctx->MEXCnt >= ctx->MEXMaxLen) {
        ctx->MEXErr = -1;
        return;
    }
    mex_term(ctx);
    if (ctx->MEXErr)
        return;

    while (1) {
        if (ctx->MEXCnt >= ctx->MEXMaxLen) {
            ctx->MEXErr = -1;
            break;  
        }
        switch (ctx->MEXNxt) {
            case '*':   mex_match(ctx, ctx->MEXNxt);
                        if (ctx->MEXErr)
                            return;
                        mex_term(ctx);
                        if (ctx->MEXErr)
                            return;
                        ctx->MEXTyp[ctx->MEXCnt] = 3;
                        ctx->MEXVal[ctx->MEXCnt] = 2.0;              
                        ctx->MEXCnt++;
                        break;
            case '/':   mex_match(ctx, ctx->MEXNxt);
                        if (ctx->MEXErr)
                            return;
                        mex_term(ctx);
                        if (ctx->MEXErr)
                            return;
                        ctx->MEXTyp[ctx->MEXCnt] = 5;
                        ctx->MEXVal[ctx->MEXCnt] = 2.0;              
                        ctx->MEXCnt++;
                        break;
            default:    return;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  mex_term                                                                */

void mex_term(TDAContext *ctx)
{
    if (ctx->MEXCnt >= ctx->MEXMaxLen) {
        ctx->MEXErr = -1;
        return;
    }
    mex_factor(ctx);
    if (ctx->MEXErr)
        return;

    while (1) {
        if (ctx->MEXCnt >= ctx->MEXMaxLen) {
            ctx->MEXErr = -1;
            break;  
        }
        switch (ctx->MEXNxt) {

            case '^':   mex_match(ctx, ctx->MEXNxt);
                        if (ctx->MEXErr)
                            return;
                        mex_factor(ctx);
                        if (ctx->MEXErr)
                            return;
                        ctx->MEXTyp[ctx->MEXCnt] = 6;
                        ctx->MEXVal[ctx->MEXCnt] = 2.0;              
                        ctx->MEXCnt++;
                        break;
            default:    return;
        }
    }

}

/* ------------------------------------------------------------------------ */
/*  mex_factor                                                              */
/*                                                                          */
/*  MEXErr =  -2  if unknown operator or expression.                        */
/*           -25  if misplaced brackets.                                    */
    
void mex_factor(TDAContext *ctx)
{
    int i,tv,na;

    if (ctx->MEXCnt >= ctx->MEXMaxLen) {
        ctx->MEXErr = -1;
        return;
    }
    switch (ctx->MEXNxt) {
        case '(': 
                    mex_match(ctx, '(');
                    if (ctx->MEXErr)
                        return;
                    mex_expr(ctx);
                    if (ctx->MEXErr)
                        return;
                    mex_match(ctx, ')');
                    if (ctx->MEXNxt == '(')  
                        ctx->MEXErr = -25;
                    break;

        case NUM:   ctx->MEXTyp[ctx->MEXCnt] = 0;
                    ctx->MEXVal[ctx->MEXCnt] = ctx->MEXTVal;
                    ctx->MEXCnt++;
                    mex_match(ctx, NUM);
                    if (ctx->MEXNxt == '(')  
                        ctx->MEXErr = -25;
                    break;

        case IDV:                                               /* variable */
                    ctx->MEXTyp[ctx->MEXCnt] = (int)ctx->MEXTVal;  
                    ctx->MEXVal[ctx->MEXCnt++] = (double)(iabs(ctx, (int)ctx->MEXTVal) - ctx->VOFFS);
                    mex_match(ctx, IDV);
                    break;

        case IDM:                                               /* matrix */
                    ctx->MEXTyp[ctx->MEXCnt] = (int)ctx->MEXTVal;        
                    ctx->MEXVal[ctx->MEXCnt++] = (double)(iabs(ctx, (int)ctx->MEXTVal) - ctx->MOFFS);
                    mex_match(ctx, IDM);
                    break; 

        case IDG:   tv = (int)ctx->MEXTVal;
                    iabs(ctx, tv);
                    na = ctx->MEXNArg;

                    mex_match(ctx, IDG);
                    if (ctx->MEXErr)  
                        return;

                    mex_match(ctx, '(');
                    if (ctx->MEXErr)
                        return;

                    mex_expr(ctx);
                    if (ctx->MEXErr)
                        return;

                    for (i = 1; i < na; ++i) {
                        if (ctx->MEXCnt >= ctx->MEXMaxLen) {
                            ctx->MEXErr = -1;
                            return; 
                        }
                        mex_match(ctx, ',');
                        if (ctx->MEXErr)
                            return;

                        mex_expr(ctx);
                        if (ctx->MEXErr)
                            return;
                    }
                    mex_match(ctx, ')');
                    if (ctx->MEXCnt >= ctx->MEXMaxLen) {
                        ctx->MEXErr = -1;
                        return; 
                    }
                    ctx->MEXTyp[ctx->MEXCnt] = tv;
                    ctx->MEXVal[ctx->MEXCnt] = (double)na;

                    ctx->MEXCnt++;
                    if (ctx->MEXNxt == '(')  
                        ctx->MEXErr = -25;
                    break;

        default:    ctx->MEXErr = -2;
                    return;
    }
}

/* ------------------------------------------------------------------------ */
/*  mex_match                                                               */

void mex_match(TDAContext *ctx, int t)
{
    if (ctx->MEXNxt == t)
        ctx->MEXNxt = mex_lexan(ctx);
    else {
        ctx->MEXErr = -3;
        return; 
    }
}

/* ------------------------------------------------------------------------ */
/*  mex_lexan: scan the input stream.                                       */
/*                                                                          */
/*  Return: IDV   for variables                                             */
/*          IDM   for matrices                                              */
/*          IDG   for operators                                             */
/*                                                                          */
/*  Sets the following error codes in MEXErr:                               */
/*  MEXErr = 0  no error                                                    */
/*          -2  undefined operator                                          */
/*          -3  syntax error                                                */
/*         -10  undefined name of matrix                                    */  

int mex_lexan(TDAContext *ctx)              
{
    int i,fnd,l,len;  
    char c;
    double tmp;

    if (*ctx->MEXESP == '-') {
        if (ctx->MEXLLast == 1) {  /* if "-" should be part of following expression */
            ctx->MEXLMin = 1;
            ctx->MEXLLast = 0;
            ctx->MEXESP++;
        }
        else {
            ctx->MEXLLast = 0;             
            return((int)*ctx->MEXESP++);
        }
    }
    ctx->MEXLLast = 0;

    if (check_vname(ctx, ctx->MEXESP)) {     /* check for variable/matrix name */

        l = get_vnlen(ctx, ctx->MEXESP);
        c = *(ctx->MEXESP + l);
        *(ctx->MEXESP + l) = '\0';

        i = mat_getidx(ctx, ctx->MEXESP,0);
        *(ctx->MEXESP + l) = c;

        if (i >= 0) {               /* found matrix name */
            ctx->MEXESP += l;
            ctx->MEXTVal = (double)(ctx->MOFFS + i);
            if (ctx->MEXLMin) { 
                ctx->MEXTVal = -ctx->MEXTVal;
                ctx->MEXLMin = 0;
            }
            return(IDM);
        }
        if ((i = v_search1(ctx, ctx->MEXESP,&len)) >= 0) {   /* search for var name */
            ctx->MEXESP += len;
            ctx->MEXTVal = (double)(ctx->VOFFS + i);
            if (ctx->MEXLMin) { 
                ctx->MEXTVal = -ctx->MEXTVal;
                ctx->MEXLMin = 0;
            }
            return(IDV);
        }
        ctx->MEXESP += l;
        ctx->MEXErr = -10;
    }
    else if (*ctx->MEXESP >= 'a' && *ctx->MEXESP <= 'z') {    /* search for operator */

        ctx->MEXTVal = 0.0;
        fnd = i = 0;
        while (1) {
            if (MOPS[i].typ == 0)
                break;

            if (!strncmp(ctx->MEXESP,MOPS[i].op,(size_t)(MOPS[i].len))) {
                ctx->MEXNArg = MOPS[i].narg;
                ctx->MEXTVal = MOPS[i].typ;
                ctx->MEXESP += MOPS[i].len;
                fnd = 1;
                break;
            }
            i++;
        }
        if (fnd) {
            if (ctx->MEXLMin) { 
                ctx->MEXTVal = -ctx->MEXTVal;
                ctx->MEXLMin = 0;
            }
        }
        else
            ctx->MEXErr = -2;

        return(IDG);
    }
    else if (isdigit((int)*ctx->MEXESP) || *ctx->MEXESP == '.') {  /* number */
        if (ctx->MEXLMin) {
            ctx->MEXLMin = 0;
            ctx->MEXESP--;
        }
        if (sscanf(ctx->MEXESP,"%lf",&tmp) == 1) {
            ctx->MEXTVal = tmp; 
            ctx->MEXESP = skip_const(ctx, ctx->MEXESP);  
            if (ctx->MEXESP == NULL)
                ctx->MEXErr = -6;
        }
        else  
            ctx->MEXErr = -6;
 
        return(NUM);
    }
    else if (ctx->MEXLMin) {

        if (*ctx->MEXESP == '.') {
            ctx->MEXLMin = 0;
            ctx->MEXESP--;
            if (sscanf(ctx->MEXESP,"%lf",&tmp) == 1) {
                ctx->MEXTVal = tmp; 
                ctx->MEXESP = skip_const(ctx, ctx->MEXESP);  
                if (ctx->MEXESP == NULL)
                    ctx->MEXErr = -6;
            }
            else
                ctx->MEXErr = -6;
            return(NUM);
        }
        else if (*ctx->MEXESP == '(') {
            ctx->MEXLMin = 0;
            ctx->MEX0Last = 1;
            ctx->MEXTVal = 0.0;
            return(NUM);
        }
        else {
            ctx->MEXErr = -7;
            return((int)*ctx->MEXESP++);
        }
    }
    else if (ctx->MEX0Last) {
        ctx->MEX0Last = 0;
        return((int) '-');
    }
    else {
        switch (*ctx->MEXESP) {
            case '+':
            case '*':
            case '^':
            case '(':   
            case ',':   ctx->MEXLLast = 1;
                        break;
            default:    break;
        }
        return((int)*ctx->MEXESP++);
    }
    return(0);
}

/*--##----------------------------------------------------------------------*/
/*  mparse()    Execute parse command for matrix expressions.               */
/*                                                                          */  
/*              mparse(                                                     */
/*                  mfmt=...,               print format, def. 12.4         */  
/*              ) = matrix_expression.                                      */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int mparse(TDAContext *ctx)           
{
    register int i,j;
    int r,err = -1;

    if (check_cmd(ctx, 1))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 6,8,1))       /* get parameters */
        goto MPFin;

    if (alloc_mex(ctx, 1000,1))
        goto MPFin;

    if ((r = mex_parse(ctx, ctx->PMRHSTR)) || ctx->MEXCnt == 0) {
        mex_error(ctx, r); 
        goto MPFin;
    }
    m_parse(ctx,(short)(ctx->MEXCnt),ctx->MEXTyp,ctx->MEXVal);      /* show info */
           
    if ((r = m_eval(ctx, 0,ctx->MEXCnt,ctx->MEXTyp,ctx->MEXVal))) {                       
        mex_error(ctx, r); 
        goto MPFin;
    }
    printf1(ctx, "\nDimension: %d x %d. ",ctx->MEXOPRow[0],ctx->MEXOPCol[0]);
    printf1(ctx, "Max stack length: %d. ",ctx->MEXOPMaxL);
    printf1(ctx, "Max matrix size: %d\n\n",ctx->MEXOPMaxM);
          
    if (alloc_mex_eval(ctx, 1,ctx->MEXOPMaxL,ctx->MEXOPMaxM))
        goto MPFin;
             
    if ((r = m_eval(ctx, 1,ctx->MEXCnt,ctx->MEXTyp,ctx->MEXVal))) {                       
        mex_error(ctx, r); 
        goto MPFin;
    }
    for (i = 0; i < ctx->MEXOPRow[0]; ++i) {
        for (j = 1; j <= ctx->MEXOPCol[0]; ++j)  
            rt_printf1_d(ctx, ctx->PMMFmtS,ctx->MEXOP[0][i * ctx->MEXOPCol[0] + j]);
        newline(ctx);
    }
    newline(ctx);
    err = 0;

MPFin:       
    alloc_mex_eval(ctx, 0,0,0);                  
    alloc_mex(ctx, 0,0);
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_parse(cnt,typ)    Show parser stack for matrix expression on stdout.  */

void m_parse(TDAContext *ctx, short cnt,int *typ,double *val)
{
    register int i;
    int t,ta;

    printf1(ctx, "\n Cnt        Typ        Val   Dimension\n");
    prnchar(ctx, '-',38,1);

    for (i = 0; i < cnt; ++i) {
        t = typ[i];
        ta = iabs(ctx, t);
        printf1(ctx, "%4d %10d %10.4lf  ",i,t,val[i]);
        if (ta >= ctx->MOFFS)  
            printf1(ctx, "%5d x %d\n",ctx->MatRow[ta-ctx->MOFFS],ctx->MatCol[ta-ctx->MOFFS]);
        else
            printf1(ctx, "\n");
    }
}

/* -##--------------------------------------------------------------------- */
/*  m_eval(opt,cnt,typ,val)                                                 */
/*                                                                          */
/*  Evaluate the matrix expression defined in typ and val. Number of        */
/*  stack entries is cnt.                                                   */
/*                                                                          */
/*  If opt = 0, check for syntactical correctness and calculate:            */
/*  MEXOPMaxL       max stack length                                        */
/*  MEXOPMaxM       max matrix size                                         */
/*                                                                          */
/*  Otherwise: calculate result in MEXOP[0].                                */
/*                                                                          */
/*  Always calculate:                                                       */
/*  MEXOPRow[0]     number of rows of resulting matrix                      */
/*  MEXOPCol[0]     number of columns of resulting matrix                   */
/*                                                                          */
/*  Return 0 if ok, otherwise one of following error indicators:            */
/*                                                                          */
/* -99  =   fatal error in stack organization                               */
/*   1  =   out of stack space                                              */
/*   2  =   dimension error                                                 */
/*   3  =   cannot evaluate sqrt                                            */
/*   4  =   trace requires square matrix                                    */
/*   5  =   cannot evaluate power function                                  */
/*   6  =   division requires scalar                                        */
/*   7  =   division by zero                                                */
/*   8  =   exp out of range                                                */
/*   9  =   cannot calculate generalized inverse                            */

int m_eval(TDAContext *ctx, int opt,int cnt,int *typ,double *val)                           
{
    register int i,j,k;
    int is,jp,jp1,jp2,n,m,r,narg,utyp,atyp;
    double tmp,tmp1;
          
    /**
    for (i = 0; i < cnt; ++i) {
        printf1(ctx, "i=%d typ=%d val=%lf \n",i,typ[i],val[i]);
    }
    **/

    ctx->MEXOPMaxM = 1;    
          
    jp = 0;
    for (is = 0; is < cnt; ++is) {

        if (jp >= MEXOPMax)  
            return(1);

        utyp = typ[is];   
        atyp = iabs(ctx, utyp);
        narg = (int)val[is];

        if (atyp == 0) {                /* number */
            if (opt && jp >= MEXOPMax)  
                return(1);

            ctx->MEXOPRow[jp] = ctx->MEXOPCol[jp] = 1;            
            if (opt)
                ctx->MEXOP[jp][1] = val[is];
            jp++;
        }
        else if (atyp >= ctx->MOFFS) {            /* matrix */

            k = atyp - ctx->MOFFS;
            if (opt && jp >= MEXOPMax)  
                return(1);

            ctx->MEXOPRow[jp] = n = ctx->MatRow[k];
            ctx->MEXOPCol[jp] = m = ctx->MatCol[k];
            ctx->MEXOPMaxM = imax(ctx, ctx->MEXOPMaxM,n * m);
                
            if (opt) {
                r = n * m;                  
                for (i = 1; i <= r; ++i) {
                    if (utyp >= 0)
                        ctx->MEXOP[jp][i] = ctx->MatVal[k][i];
                    else
                        ctx->MEXOP[jp][i] = -ctx->MatVal[k][i];
                }
            }
            jp++;
        }
        else if (narg == 1) {                   /* one argument */

            if (jp < 1)  
                return(-99);
 
            jp1 = jp - 1;

            if (atyp == 104) {                      /* sqrt */
                if (opt) {
                    r = ctx->MEXOPRow[jp1] * ctx->MEXOPCol[jp1];
                    for (i = 1; i <= r; ++i) {
                        if (ctx->MEXOP[jp1][i] >= 0.0)
                            ctx->MEXOP[jp1][i] = sqrt(ctx->MEXOP[jp1][i]);
                        else
                            return(3);
                    }
                }
            }
            else if (atyp == 100) {                     /* ginv */
                n = ctx->MEXOPRow[jp1];
                m = ctx->MEXOPCol[jp1];

                if (n < m)
                    ctx->MEXOPMaxL = imax(ctx, ctx->MEXOPMaxL,jp + 1);

                if (opt) {
                    if (n >= m) {
                        r = ginv(ctx, n,m,ctx->MEXOP[jp1]);
                        if (r < 0)
                            return(9);
                    }
                    else {                      /* requires transposition */
                        r = n * m;
                        for (i = 0; i < m; ++i) {
                            for (j = 1; j <= n; ++j)
                                ctx->MEXOP[jp][i * n + j] = ctx->MEXOP[jp1][(j - 1) * m + i + 1];
                        }
                        r = ginv(ctx, m,n,ctx->MEXOP[jp]);
                        if (r < 0)
                            return(9);

                        for (i = 0; i < m; ++i) {
                            for (j = 1; j <= n; ++j)
                                ctx->MEXOP[jp1][i * n + j] = ctx->MEXOP[jp][(j - 1) * m + i + 1];
                        }
                    }
                }
                ctx->MEXOPRow[jp1] = m;
                ctx->MEXOPCol[jp1] = n; 
            }
            else if (atyp == 105) {                     /* eexp */
                if (opt) {
                    r = ctx->MEXOPRow[jp1] * ctx->MEXOPCol[jp1];
                    for (i = 1; i <= r; ++i) {
                        tmp = ctx->MEXOP[jp1][i];
                        if (tmp < ExpMin || tmp > ExpMax)
                            return(8);
                        ctx->MEXOP[jp1][i] = exp(tmp) / (1.0 + exp(tmp));
                    }           
                }
            }
            else if (atyp == 103) {                 /* trace */
                if (check_square(ctx, jp1)) {
                    if (opt) {
                        tmp = 0.0;
                        m = ctx->MEXOPCol[jp1];
                        for (i = 0; i < m; ++i)
                            tmp += ctx->MEXOP[jp1][i * m + i + 1];
                        ctx->MEXOP[jp1][1] = tmp;
                    }   
                    ctx->MEXOPRow[jp1] = ctx->MEXOPCol[jp1] = 1;
                }
                else
                    return(4);
            }
            else if (atyp == 102) {             /* transpose */

                ctx->MEXOPMaxL = imax(ctx, ctx->MEXOPMaxL,jp + 1);
                n = ctx->MEXOPRow[jp1];  
                m = ctx->MEXOPCol[jp1];  

                if (opt) {
                    r = n * m;
                    for (i = 1; i <= r; ++i)
                        ctx->MEXOP[jp][i] = ctx->MEXOP[jp1][i];

                    for (i = 0; i < m; ++i) {
                        for (j = 1; j <= n; ++j)
                            ctx->MEXOP[jp1][i * n + j] = ctx->MEXOP[jp][(j - 1) * m + i + 1];
                    }
                }
                ctx->MEXOPRow[jp1] = m;
                ctx->MEXOPCol[jp1] = n;
            }
            else if (atyp == 101) {                         /* cross product */

                ctx->MEXOPMaxL = imax(ctx, ctx->MEXOPMaxL,jp + 1);
                n = ctx->MEXOPRow[jp1];  
                m = ctx->MEXOPCol[jp1];  

                if (opt) {
                    for (i = 1; i <= m; ++i) {
                        for (j = 1; j <= m; ++j) {
                            tmp = 0.0;
                            for (k = 0; k < n; ++k) 
                                tmp += ctx->MEXOP[jp1][k * m + i] * ctx->MEXOP[jp1][k * m + j];
                            ctx->MEXOP[jp][(i - 1) * m + j] = tmp;
                        }
                    }
                    r = m * m;
                    for (i = 1; i <= r; ++i)
                        ctx->MEXOP[jp1][i] = ctx->MEXOP[jp][i];
                }
                ctx->MEXOPRow[jp1] = m;
                ctx->MEXOPMaxM = imax(ctx, ctx->MEXOPMaxM,m * m);
            }
            else                    /* undefined */
                gerr_exit(ctx, 4);

            if (utyp < 0) {                         /* change sign */
                r = ctx->MEXOPRow[jp1] * ctx->MEXOPCol[jp1];
                for (i = 1; i <= r; ++i)
                    ctx->MEXOP[jp1][i] = -ctx->MEXOP[jp1][i];
            }
        }   
        else if (narg == 2) {                   /* two arguments */

            if (jp < 2) 
                return(-99);
 
            jp1 = jp - 2;
            jp2 = jp - 1;

            if (atyp == 1) {                                /* addition */

                if (check_scalar(ctx, jp1)) {

                    n = ctx->MEXOPRow[jp2];  
                    m = ctx->MEXOPCol[jp2];  

                    if (opt) {
                        tmp = ctx->MEXOP[jp1][1];
                        r = n * m;
                        for (i = 1; i <= r; ++i)
                            ctx->MEXOP[jp1][i] = ctx->MEXOP[jp2][i] + tmp;
                    }
                    ctx->MEXOPRow[jp1] = n;              
                    ctx->MEXOPCol[jp1] = m;              
                    jp--;
                }
                else if (check_scalar(ctx, jp2)) {

                    n = ctx->MEXOPRow[jp1];  
                    m = ctx->MEXOPCol[jp1];  

                    if (opt) {
                        tmp = ctx->MEXOP[jp2][1];
                        r = n * m;
                        for (i = 1; i <= r; ++i)
                            ctx->MEXOP[jp1][i] += tmp;
                    }
                    jp--;
                }
                else {
                    if (check_dim_add(ctx, jp1,jp2))  
                        return(2);
     
                    if (opt) {
                        n = ctx->MEXOPRow[jp1] * ctx->MEXOPCol[jp1];  
                        for (i = 1; i <= n; ++i)
                            ctx->MEXOP[jp1][i] += ctx->MEXOP[jp2][i];
                    }
                    jp--;
                }
            }
            else if (atyp == 2) {                           /* subtraction */

                if (check_scalar(ctx, jp1)) {

                    n = ctx->MEXOPRow[jp2];  
                    m = ctx->MEXOPCol[jp2];  

                    if (opt) {
                        tmp = ctx->MEXOP[jp1][1];
                        r = n * m;
                        for (i = 1; i <= r; ++i)
                            ctx->MEXOP[jp1][i] = tmp - ctx->MEXOP[jp2][i];
                    }
                    ctx->MEXOPRow[jp1] = n;              
                    ctx->MEXOPCol[jp1] = m;              
                    jp--;
                }
                else if (check_scalar(ctx, jp2)) {

                    n = ctx->MEXOPRow[jp1];  
                    m = ctx->MEXOPCol[jp1];  

                    if (opt) {
                        tmp = ctx->MEXOP[jp2][1];
                        r = n * m;
                        for (i = 1; i <= r; ++i)
                            ctx->MEXOP[jp1][i] -= tmp;
                    }
                    jp--;
                }
                else {
                    if (check_dim_add(ctx, jp1,jp2))  
                        return(2);
     
                    if (opt) {
                        n = ctx->MEXOPRow[jp1] * ctx->MEXOPCol[jp1];  
                        for (i = 1; i <= n; ++i)
                            ctx->MEXOP[jp1][i] -= ctx->MEXOP[jp2][i];
                    }
                    jp--;
                }
            }
            else if (atyp == 3) {                           /* multiplication */

                if (check_scalar(ctx, jp2)) {
                    if (opt) {
                        n = ctx->MEXOPRow[jp1] * ctx->MEXOPCol[jp1];
                        tmp = ctx->MEXOP[jp2][1];
                        for (i = 1; i <= n; ++i)
                            ctx->MEXOP[jp1][i] *= tmp;           
                    }
                    jp--;
                }
                else if (check_scalar(ctx, jp1)) {
                    if (opt) {
                        n = ctx->MEXOPRow[jp2] * ctx->MEXOPCol[jp2];
                        tmp = ctx->MEXOP[jp1][1];
                        for (i = 1; i <= n; ++i)
                            ctx->MEXOP[jp1][i] = tmp * ctx->MEXOP[jp2][i];
                    }
                    ctx->MEXOPRow[jp1] = ctx->MEXOPRow[jp2];
                    ctx->MEXOPCol[jp1] = ctx->MEXOPCol[jp2];
                    jp--;
                }
                else {
                    if (check_dim_mul(ctx, jp1,jp2))  
                        return(2);
 
                    ctx->MEXOPMaxL = imax(ctx, ctx->MEXOPMaxL,jp + 1);
                    n = ctx->MEXOPRow[jp1];  
                    m = ctx->MEXOPCol[jp2];  

                    if (opt) {
                        if (jp >= MEXOPMax)  
                            return(1);

                        r = ctx->MEXOPCol[jp1];  
                        for (i = 0; i < n; ++i) {
                            for (j = 1; j <= m; ++j) {
                                tmp = 0.0;
                                for (k = 0; k < r; ++k) 
                                    tmp += ctx->MEXOP[jp1][i * r + k + 1] *
                                           ctx->MEXOP[jp2][k * m + j];
                                
                                ctx->MEXOP[jp][i * m + j] = tmp; 
                            }
                        }
                        r = n * m;
                        for (i = 1; i <= r; ++i)
                            ctx->MEXOP[jp1][i] = ctx->MEXOP[jp][i];
                    }
                    ctx->MEXOPRow[jp1] = n;              
                    ctx->MEXOPCol[jp1] = m;              
                    ctx->MEXOPMaxM = imax(ctx, ctx->MEXOPMaxM,n * m);
                    jp--;
                }
            }
            else if (atyp == 5) {                           /* division */

                if (check_scalar(ctx, jp2) == 0)  
                    return(6);

                if (opt) {
                    tmp = ctx->MEXOP[jp2][1];
                    if (tmp == 0.0)
                        return(7);
                    r = ctx->MEXOPRow[jp1] * ctx->MEXOPCol[jp1];
                    for (i = 1; i <= r; ++i)
                        ctx->MEXOP[jp1][i] /= tmp;           
                }
                jp--;
            }
            else if (atyp == 6) {                        /* ^ */
                if (check_scalar(ctx, jp2) == 0)  
                    return(5);

                if (opt) {
                    tmp1 = ctx->MEXOP[jp2][1];
                    r = ctx->MEXOPRow[jp1] * ctx->MEXOPCol[jp1];
                    for (i = 1; i <= r; ++i) {
                        tmp = ctx->MEXOP[jp1][i];
                        if (fabs(tmp) <= ctx->EPSI2)
                            ctx->MEXOP[jp1][i] = 0.0;
                        else if (tmp > 0.0)
                            ctx->MEXOP[jp1][i] = pow(tmp,tmp1);
                        else {
                            if (floor(tmp1) != tmp1)
                                return(5);
                            ctx->MEXOP[jp1][i] = pow(tmp,tmp1);
                        }
                    }   
                }
                jp--;
            }

            else                    /* undefined */
                gerr_exit(ctx, 4);

            if (utyp < 0) {                         /* change sign */
                r = ctx->MEXOPRow[jp1] * ctx->MEXOPCol[jp1];
                for (i = 1; i <= r; ++i)
                    ctx->MEXOP[jp1][i] = -ctx->MEXOP[jp1][i];
            }
        }
        ctx->MEXOPMaxL = imax(ctx, ctx->MEXOPMaxL,jp);
    } 
    if (jp != 1)           
        return(-99);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  check_dim_add(i,j)      Return 0 if OK, -1 if dimension error           */

int check_dim_add(TDAContext *ctx, int i,int j)
{
    if (ctx->MEXOPRow[i] != ctx->MEXOPRow[j] || ctx->MEXOPCol[i] != ctx->MEXOPCol[j])       
        return(-1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  check_dim_mul(i,j)      Return 0 if OK, -1 if dimension error           */

int check_dim_mul(TDAContext *ctx, int i,int j)
{
    if (ctx->MEXOPCol[i] != ctx->MEXOPRow[j])       
        return(-1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  check_scalar(i)     Return 1 if scalar, otherwise 0.                    */

int check_scalar(TDAContext *ctx, int i)
{
    if (ctx->MEXOPCol[i] == 1 && ctx->MEXOPRow[i] == 1)       
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  check_square(i)     Return 1 if square, otherwise 0.                    */

int check_square(TDAContext *ctx, int i)
{
    if (ctx->MEXOPCol[i] == ctx->MEXOPRow[i])       
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_mex_eval(opt,n,m)                                                 */
/*                                                                          */
/*  If opt != 0 allocate stack for the evaluation of a matrix expression.   */
/*  n is stack depth, m is max size of matrices.                            */
/*                                                                          */
/*  Otherwise free previously allocated memory.                             */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int alloc_mex_eval(TDAContext *ctx, int opt,int n,int m)
{
    int i;

    if (opt) {
        if (n >= MEXOPMax) {
            printf1(ctx, "Error: exceeding max stack size (%d).\n",MEXOPMax);
            return(-1);
        }
        ctx->MEXOPMaxB = m;

        for (i = 0; i < n; ++i) {
            if (!(ctx->MEXOP[i] = (double *)calloc((size_t)(m + 1),sizeof(double)))) { 
                p_err(ctx, -2,1);  
                goto AMFin;              
            }
            memrq(ctx, m + 1,sizeof(double));
            ctx->MEXOPMaxA++;
        }
        return(0);
    }

AMFin:
    for (i = 0; i < ctx->MEXOPMaxA; ++i) {
        free((char *)ctx->MEXOP[i]);
        memrq(ctx, -ctx->MEXOPMaxB - 1,sizeof(double));
    }
    ctx->MEXOPMaxA = ctx->MEXOPMaxB = 0;
    if (opt)
        return(-1);
    else
        return(0);
}

/* ------------------------------------------------------------------------ */
/*  mex_error(n)    Report error n.                                         */

void mex_error(TDAContext *ctx, int n)
{
    printf1(ctx, "\nError (%d) in evaluation of matrix expression.\n",n);
    printf1(ctx, "Error type: ");
    switch (n) {
        case -1:    printf1(ctx, "stack overflow.\n");
                    break;
        case -2:    printf1(ctx, "unknown operator or expression, or misplaced brackets.\n");
                    break;
        case -3:    printf1(ctx, "syntax error.\n");
                    break;
        case -6:    printf1(ctx, "scan error floating point number.\n");
                    break;
        case -7:    printf1(ctx, "misplaced minus sign.\n");
                    break;
        case -10:   printf1(ctx, "undefined matrix or variable name.\n");
                    break;
        case -25:   printf1(ctx, "misplaced brackets.\n");
                    break;

        case  1:    printf1(ctx, "out of stack space.\n");
                    break;
        case  2:    printf1(ctx, "error in matrix dimensions.\n");
                    break;
        case  3:    printf1(ctx, "cannot evaluate sqrt.\n");
                    break;
        case  4:    printf1(ctx, "trace requires square matrix.\n");
                    break;
        case  5:    printf1(ctx, "cannot evaluate power function.\n");
                    break;
        case  6:    printf1(ctx, "division requires scalar.\n");
                    break;
        case  7:    printf1(ctx, "division by zero.\n");
                    break;
        case  8:    printf1(ctx, "exp out of range.\n");
                    break;
        case  9:    printf1(ctx, "cannot calculate generalized inverse.\n");
                    break;
        default:    printf1(ctx, "unknown.\n");
                    break;
    }
}

/*--##----------------------------------------------------------------------*/
/*  mex_eval(s,row,col)     Evaluate matrix expression s. Put Result        */
/*                          into MEXOP[0], NEXOPRow[0], MEXOPCol[0].        */  
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int mex_eval(TDAContext *ctx, char *s)
{
    int err,r;

    err = -1;
    if (alloc_mex(ctx, 1000,1))
        goto MEVFin;

    if ((r = mex_parse(ctx, s)) || ctx->MEXCnt == 0) {
        mex_error(ctx, r); 
        goto MEVFin;
    }
    if ((r = m_eval(ctx, 0,ctx->MEXCnt,ctx->MEXTyp,ctx->MEXVal))) {                       
        mex_error(ctx, r); 
        goto MEVFin;
    }
    if (alloc_mex_eval(ctx, 1,ctx->MEXOPMaxL,ctx->MEXOPMaxM))
        goto MEVFin;
             
    if ((r = m_eval(ctx, 1,ctx->MEXCnt,ctx->MEXTyp,ctx->MEXVal))) {                       
        mex_error(ctx, r); 
        goto MEVFin;
    }
    err = 0;

MEVFin:       
    return(err);
}

