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

/* ------------------------------------------------------------------------ */
/*  functions in t_eval4.c                                                   */

int alloc_mex(int n,int opt);
int mex_parse(char *s);
void mex_expr(void);
void mex_expr1(void);
void mex_term(void);
void mex_factor(void);
void mex_match(int t);
int mex_lexan(void);             
int mparse(void);          
void m_parse(short cnt,int *typ,double *val);
int m_eval(int opt,int cnt,int *typ,double *val);                          
int check_dim_add(int i,int j);
int check_dim_mul(int i,int j);
int check_scalar(int i);
int check_square(int i);
int alloc_mex_eval(int opt,int n,int m);
void mex_error(int n);
int mex_eval(char *s);


/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

#define NUM 257         /* return code for floating point numbers           */
#define IDV 258         /* return code for V, C or X variables              */
#define IDM 259         /* return code for matrices                         */
#define IDG 260         /* return code for operators                        */
#define IDC 261         /* return code for constants, operators without     */

int MEXMaxLen = 0;      /* length of stack                                  */
int *MEXTyp;            /* type of stack entry                              */
double *MEXVal;         /* value of operand                                 */
int MEXCnt = 0;         /* number of entries in stack                       */
int MEXErr = 0;         /* error from parser                                */
int MEXNxt  = 0;        /* look ahead value in scanning input stream        */
double MEXTVal = 0.0;   /* value of operand produced by lexan()             */
int MEXNArg = 0;        /* number of arguments produced by lexan()          */
int MEX0Last = 0;       /* flags used by lexan()                            */
int MEXLLast = 1;       /* flags used by lexan()                            */
int MEXLMin  = 0;       /* flags used by lexan()                            */
char *MEXESP;           /* global pointer to string for evaluation          */

#define MEXOPMax 100
double *MEXOP[MEXOPMax + 2];
int MEXOPRow[MEXOPMax + 2];             
int MEXOPCol[MEXOPMax + 2];  
int MEXOPMaxA = 0;              /* allocated entries for MEXOP              */
int MEXOPMaxB = 0;              /* allocated entries for MEXOP[]            */
int MEXOPMaxL = 0;              /* max stack length                         */
int MEXOPMaxM = 0;              /* max matrix dimension                     */

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
    "cross",  101,  5,  1,
    "eexp",   105,  4,  1,
    "ginv",   100,  4,  1,
    "sqrt",   104,  4,  1,
    "trace",  103,  5,  1,
    "t",      102,  1,  1,
    "",         0,  0,  0,
};


/* ------------------------------------------------------------------------ */
/*  alloc_mex(n,opt)                                                        */
/*                                                                          */
/*  Allocate (opt = 1) or free (opt = 0) MEXTyp with n elements.            */
/*  Return 0 if OK, -1 if error.                                            */

int alloc_mex(int n,int opt)
{
    if (MEXMaxLen > 0) {
        free((char *)MEXTyp);
        free((char *)MEXVal);
        memrq(-MEXMaxLen - 2,sizeof(int) + sizeof(double));
        MEXMaxLen = 0;
    }
    if (opt == 1) {
        if (!(MEXTyp = (int *)calloc(n + 2,sizeof(int)))) { 
            p_err(-2,1);  
            return(-1);
        }
        if (!(MEXVal = (double *)calloc(n + 2,sizeof(double)))) { 
            free((char *)MEXTyp);
            p_err(-2,1);  
            return(-1);
        }
        MEXMaxLen = n;
        memrq(n + 2,sizeof(int) + sizeof(double));
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

int mex_parse(char *s)
{
    MEXLMin = MEXErr = MEXCnt = 0;
    MEXESP = s;
    MEXLLast = 1;
    MEX0Last = 0;

    MEXNxt = mex_lexan();
    if (MEXErr)
        return(MEXErr);
    while (MEXNxt) {
        mex_expr();
        if (MEXErr)
            return(MEXErr);
    }
    return(MEXErr);
}

/* ------------------------------------------------------------------------ */
/*  mex_expr                                                                */

void mex_expr(void)
{
    if (MEXCnt >= MEXMaxLen) {
        MEXErr = -1;
        return;
    }
    mex_expr1();
    if (MEXErr)
        return;

    while (1) {
        if (MEXCnt >= MEXMaxLen) {
            MEXErr = -1;
            break;  
        }
        switch (MEXNxt) {
            case '+':   mex_match(MEXNxt);
                        if (MEXErr)
                            return;
                        mex_expr1();
                        if (MEXErr)
                            return;
                        MEXTyp[MEXCnt] = 1;
                        MEXVal[MEXCnt] = 2.0;              
                        MEXCnt++;
                        break;
            case '-':   mex_match(MEXNxt);
                        if (MEXErr)
                            return;
                        mex_expr1();
                        if (MEXErr)
                            return;
                        MEXTyp[MEXCnt] = 2;
                        MEXVal[MEXCnt] = 2.0;              
                        MEXCnt++;
                        break;

            default:    return;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  mex_expr1                                                               */

void mex_expr1(void)
{
    if (MEXCnt >= MEXMaxLen) {
        MEXErr = -1;
        return;
    }
    mex_term();
    if (MEXErr)
        return;

    while (1) {
        if (MEXCnt >= MEXMaxLen) {
            MEXErr = -1;
            break;  
        }
        switch (MEXNxt) {
            case '*':   mex_match(MEXNxt);
                        if (MEXErr)
                            return;
                        mex_term();
                        if (MEXErr)
                            return;
                        MEXTyp[MEXCnt] = 3;
                        MEXVal[MEXCnt] = 2.0;              
                        MEXCnt++;
                        break;
            case '/':   mex_match(MEXNxt);
                        if (MEXErr)
                            return;
                        mex_term();
                        if (MEXErr)
                            return;
                        MEXTyp[MEXCnt] = 5;
                        MEXVal[MEXCnt] = 2.0;              
                        MEXCnt++;
                        break;
            default:    return;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  mex_term                                                                */

void mex_term(void)
{
    if (MEXCnt >= MEXMaxLen) {
        MEXErr = -1;
        return;
    }
    mex_factor();
    if (MEXErr)
        return;

    while (1) {
        if (MEXCnt >= MEXMaxLen) {
            MEXErr = -1;
            break;  
        }
        switch (MEXNxt) {

            case '^':   mex_match(MEXNxt);
                        if (MEXErr)
                            return;
                        mex_factor();
                        if (MEXErr)
                            return;
                        MEXTyp[MEXCnt] = 6;
                        MEXVal[MEXCnt] = 2.0;              
                        MEXCnt++;
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
    
void mex_factor(void)
{
    int i,tv,tva,na;

    tva = -1;
    if (MEXCnt >= MEXMaxLen) {
        MEXErr = -1;
        return;
    }
    switch (MEXNxt) {
        case '(': 
                    mex_match('(');
                    if (MEXErr)
                        return;
                    mex_expr();
                    if (MEXErr)
                        return;
                    mex_match(')');
                    if (MEXNxt == '(')  
                        MEXErr = -25;
                    break;

        case NUM:   MEXTyp[MEXCnt] = 0;
                    MEXVal[MEXCnt] = MEXTVal;
                    MEXCnt++;
                    mex_match(NUM);
                    if (MEXNxt == '(')  
                        MEXErr = -25;
                    break;

        case IDV:                                               /* variable */
                    MEXTyp[MEXCnt] = (int)MEXTVal;  
                    MEXVal[MEXCnt++] = (double)(iabs((int)MEXTVal) - VOFFS);
                    mex_match(IDV);
                    break;

        case IDM:                                               /* matrix */
                    MEXTyp[MEXCnt] = (int)MEXTVal;        
                    MEXVal[MEXCnt++] = (double)(iabs((int)MEXTVal) - MOFFS);
                    mex_match(IDM);
                    break; 

        case IDG:   tv = (int)MEXTVal;
                    tva = iabs(tv);
                    na = MEXNArg;

                    mex_match(IDG);
                    if (MEXErr)  
                        return;

                    mex_match('(');
                    if (MEXErr)
                        return;

                    mex_expr();
                    if (MEXErr)
                        return;

                    for (i = 1; i < na; ++i) {
                        if (MEXCnt >= MEXMaxLen) {
                            MEXErr = -1;
                            return; 
                        }
                        mex_match(',');
                        if (MEXErr)
                            return;

                        mex_expr();
                        if (MEXErr)
                            return;
                    }
                    mex_match(')');
                    if (MEXCnt >= MEXMaxLen) {
                        MEXErr = -1;
                        return; 
                    }
                    MEXTyp[MEXCnt] = tv;
                    MEXVal[MEXCnt] = (double)na;

                    MEXCnt++;
                    if (MEXNxt == '(')  
                        MEXErr = -25;
                    break;

        default:    MEXErr = -2;
                    return;
    }
}

/* ------------------------------------------------------------------------ */
/*  mex_match                                                               */

void mex_match(int t)
{
    if (MEXNxt == t)
        MEXNxt = mex_lexan();
    else {
        MEXErr = -3;
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

int mex_lexan(void)              
{
    int i,fnd,l,len;  
    char c;
    double tmp;

    if (*MEXESP == '-') {
        if (MEXLLast == 1) {  /* if "-" should be part of following expression */
            MEXLMin = 1;
            MEXLLast = 0;
            MEXESP++;
        }
        else {
            MEXLLast = 0;             
            return((int)*MEXESP++);
        }
    }
    MEXLLast = 0;

    if (check_vname(MEXESP)) {     /* check for variable/matrix name */

        l = get_vnlen(MEXESP);
        c = *(MEXESP + l);
        *(MEXESP + l) = '\0';

        i = mat_getidx(MEXESP,0);
        *(MEXESP + l) = c;

        if (i >= 0) {               /* found matrix name */
            MEXESP += l;
            MEXTVal = (double)(MOFFS + i);
            if (MEXLMin) { 
                MEXTVal = -MEXTVal;
                MEXLMin = 0;
            }
            return(IDM);
        }
        if ((i = v_search1(MEXESP,&len)) >= 0) {   /* search for var name */
            MEXESP += len;
            MEXTVal = (double)(VOFFS + i);
            if (MEXLMin) { 
                MEXTVal = -MEXTVal;
                MEXLMin = 0;
            }
            return(IDV);
        }
        MEXESP += l;
        MEXErr = -10;
    }
    else if (*MEXESP >= 'a' && *MEXESP <= 'z') {    /* search for operator */

        MEXTVal = 0.0;
        fnd = i = 0;
        while (1) {
            if (MOPS[i].typ == 0)
                break;

            if (!strncmp(MEXESP,MOPS[i].op,MOPS[i].len)) {
                MEXNArg = MOPS[i].narg;
                MEXTVal = MOPS[i].typ;
                MEXESP += MOPS[i].len;
                fnd = 1;
                break;
            }
            i++;
        }
        if (fnd) {
            if (MEXLMin) { 
                MEXTVal = -MEXTVal;
                MEXLMin = 0;
            }
        }
        else
            MEXErr = -2;

        return(IDG);
    }
    else if (isdigit((int)*MEXESP) || *MEXESP == '.') {  /* number */
        if (MEXLMin) {
            MEXLMin = 0;
            MEXESP--;
        }
        if (sscanf(MEXESP,"%lf",&tmp) == 1) {
            MEXTVal = tmp; 
            MEXESP = skip_const(MEXESP);  
            if (MEXESP == NULL)
                MEXErr = -6;
        }
        else  
            MEXErr = -6;
 
        return(NUM);
    }
    else if (MEXLMin) {

        if (*MEXESP == '.') {
            MEXLMin = 0;
            MEXESP--;
            if (sscanf(MEXESP,"%lf",&tmp) == 1) {
                MEXTVal = tmp; 
                MEXESP = skip_const(MEXESP);  
                if (MEXESP == NULL)
                    MEXErr = -6;
            }
            else
                MEXErr = -6;
            return(NUM);
        }
        else if (*MEXESP == '(') {
            MEXLMin = 0;
            MEX0Last = 1;
            MEXTVal = 0.0;
            return(NUM);
        }
        else {
            MEXErr = -7;
            return((int)*MEXESP++);
        }
    }
    else if (MEX0Last) {
        MEX0Last = 0;
        return((int) '-');
    }
    else {
        switch (*MEXESP) {
            case '+':
            case '*':
            case '^':
            case '(':   
            case ',':   MEXLLast = 1;
                        break;
            default:    break;
        }
        return((int)*MEXESP++);
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

int mparse(void)           
{
    register int i,j;
    int r,err = -1;

    if (check_cmd(1))
        return(-1);

    if (parm(CmdBuf + 6,8,1))       /* get parameters */
        goto MPFin;

    if (alloc_mex(1000,1))
        goto MPFin;

    if ((r = mex_parse(PMRHSTR)) || MEXCnt == 0) {
        mex_error(r); 
        goto MPFin;
    }
    m_parse(MEXCnt,MEXTyp,MEXVal);      /* show info */
           
    if ((r = m_eval(0,MEXCnt,MEXTyp,MEXVal))) {                       
        mex_error(r); 
        goto MPFin;
    }
    printf1("\nDimension: %d x %d. ",MEXOPRow[0],MEXOPCol[0]);
    printf1("Max stack length: %d. ",MEXOPMaxL);
    printf1("Max matrix size: %d\n\n",MEXOPMaxM);
          
    if (alloc_mex_eval(1,MEXOPMaxL,MEXOPMaxM))
        goto MPFin;
             
    if ((r = m_eval(1,MEXCnt,MEXTyp,MEXVal))) {                       
        mex_error(r); 
        goto MPFin;
    }
    for (i = 0; i < MEXOPRow[0]; ++i) {
        for (j = 1; j <= MEXOPCol[0]; ++j)  
            printf1(PMMFmtS,MEXOP[0][i * MEXOPCol[0] + j]);
        newline();
    }
    newline();
    err = 0;

MPFin:       
    alloc_mex_eval(0,0,0);                  
    alloc_mex(0,0);
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  m_parse(cnt,typ)    Show parser stack for matrix expression on stdout.  */

void m_parse(short cnt,int *typ,double *val)
{
    register int i;
    int t,ta;

    printf1("\n Cnt        Typ        Val   Dimension\n");
    prnchar('-',38,1);

    for (i = 0; i < cnt; ++i) {
        t = typ[i];
        ta = iabs(t);
        printf1("%4d %10d %10.4lf  ",i,t,val[i]);
        if (ta >= MOFFS)  
            printf1("%5d x %d\n",MatRow[ta-MOFFS],MatCol[ta-MOFFS]);
        else
            printf1("\n");
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

int m_eval(int opt,int cnt,int *typ,double *val)                           
{
    register int i,j,k;
    int is,jp,jp1,jp2,n,m,r,narg,utyp,atyp;
    double tmp,tmp1;
          
    /**
    for (i = 0; i < cnt; ++i) {
        printf1("i=%d typ=%d val=%lf \n",i,typ[i],val[i]);
    }
    **/

    MEXOPMaxM = 1;    
          
    jp = 0;
    for (is = 0; is < cnt; ++is) {

        if (jp >= MEXOPMax)  
            return(1);

        utyp = typ[is];   
        atyp = iabs(utyp);
        narg = (int)val[is];

        if (atyp == 0) {                /* number */
            if (opt && jp >= MEXOPMax)  
                return(1);

            MEXOPRow[jp] = MEXOPCol[jp] = 1;            
            if (opt)
                MEXOP[jp][1] = val[is];
            jp++;
        }
        else if (atyp >= MOFFS) {            /* matrix */

            k = atyp - MOFFS;
            if (opt && jp >= MEXOPMax)  
                return(1);

            MEXOPRow[jp] = n = MatRow[k];
            MEXOPCol[jp] = m = MatCol[k];
            MEXOPMaxM = imax(MEXOPMaxM,n * m);
                
            if (opt) {
                r = n * m;                  
                for (i = 1; i <= r; ++i) {
                    if (utyp >= 0)
                        MEXOP[jp][i] = MatVal[k][i];
                    else
                        MEXOP[jp][i] = -MatVal[k][i];
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
                    r = MEXOPRow[jp1] * MEXOPCol[jp1];
                    for (i = 1; i <= r; ++i) {
                        if (MEXOP[jp1][i] >= 0.0)
                            MEXOP[jp1][i] = sqrt(MEXOP[jp1][i]);
                        else
                            return(3);
                    }
                }
            }
            else if (atyp == 100) {                     /* ginv */
                n = MEXOPRow[jp1];
                m = MEXOPCol[jp1];

                if (n < m)
                    MEXOPMaxL = imax(MEXOPMaxL,jp + 1);

                if (opt) {
                    if (n >= m) {
                        r = ginv(n,m,MEXOP[jp1]);
                        if (r < 0)
                            return(9);
                    }
                    else {                      /* requires transposition */
                        r = n * m;
                        for (i = 0; i < m; ++i) {
                            for (j = 1; j <= n; ++j)
                                MEXOP[jp][i * n + j] = MEXOP[jp1][(j - 1) * m + i + 1];
                        }
                        r = ginv(m,n,MEXOP[jp]);
                        if (r < 0)
                            return(9);

                        for (i = 0; i < m; ++i) {
                            for (j = 1; j <= n; ++j)
                                MEXOP[jp1][i * n + j] = MEXOP[jp][(j - 1) * m + i + 1];
                        }
                    }
                }
                MEXOPRow[jp1] = m;
                MEXOPCol[jp1] = n; 
            }
            else if (atyp == 105) {                     /* eexp */
                if (opt) {
                    r = MEXOPRow[jp1] * MEXOPCol[jp1];
                    for (i = 1; i <= r; ++i) {
                        tmp = MEXOP[jp1][i];
                        if (tmp < ExpMin || tmp > ExpMax)
                            return(8);
                        MEXOP[jp1][i] = exp(tmp) / (1.0 + exp(tmp));
                    }           
                }
            }
            else if (atyp == 103) {                 /* trace */
                if (check_square(jp1)) {
                    if (opt) {
                        tmp = 0.0;
                        m = MEXOPCol[jp1];
                        for (i = 0; i < m; ++i)
                            tmp += MEXOP[jp1][i * m + i + 1];
                        MEXOP[jp1][1] = tmp;
                    }   
                    MEXOPRow[jp1] = MEXOPCol[jp1] = 1;
                }
                else
                    return(4);
            }
            else if (atyp == 102) {             /* transpose */

                MEXOPMaxL = imax(MEXOPMaxL,jp + 1);
                n = MEXOPRow[jp1];  
                m = MEXOPCol[jp1];  

                if (opt) {
                    r = n * m;
                    for (i = 1; i <= r; ++i)
                        MEXOP[jp][i] = MEXOP[jp1][i];

                    for (i = 0; i < m; ++i) {
                        for (j = 1; j <= n; ++j)
                            MEXOP[jp1][i * n + j] = MEXOP[jp][(j - 1) * m + i + 1];
                    }
                }
                MEXOPRow[jp1] = m;
                MEXOPCol[jp1] = n;
            }
            else if (atyp == 101) {                         /* cross product */

                MEXOPMaxL = imax(MEXOPMaxL,jp + 1);
                n = MEXOPRow[jp1];  
                m = MEXOPCol[jp1];  

                if (opt) {
                    for (i = 1; i <= m; ++i) {
                        for (j = 1; j <= m; ++j) {
                            tmp = 0.0;
                            for (k = 0; k < n; ++k) 
                                tmp += MEXOP[jp1][k * m + i] * MEXOP[jp1][k * m + j];
                            MEXOP[jp][(i - 1) * m + j] = tmp;
                        }
                    }
                    r = m * m;
                    for (i = 1; i <= r; ++i)
                        MEXOP[jp1][i] = MEXOP[jp][i];
                }
                MEXOPRow[jp1] = m;
                MEXOPMaxM = imax(MEXOPMaxM,m * m);
            }
            else                    /* undefined */
                gerr_exit(4);

            if (utyp < 0) {                         /* change sign */
                r = MEXOPRow[jp1] * MEXOPCol[jp1];
                for (i = 1; i <= r; ++i)
                    MEXOP[jp1][i] = -MEXOP[jp1][i];
            }
        }   
        else if (narg == 2) {                   /* two arguments */

            if (jp < 2) 
                return(-99);
 
            jp1 = jp - 2;
            jp2 = jp - 1;

            if (atyp == 1) {                                /* addition */

                if (check_scalar(jp1)) {

                    n = MEXOPRow[jp2];  
                    m = MEXOPCol[jp2];  

                    if (opt) {
                        tmp = MEXOP[jp1][1];
                        r = n * m;
                        for (i = 1; i <= r; ++i)
                            MEXOP[jp1][i] = MEXOP[jp2][i] + tmp;
                    }
                    MEXOPRow[jp1] = n;              
                    MEXOPCol[jp1] = m;              
                    jp--;
                }
                else if (check_scalar(jp2)) {

                    n = MEXOPRow[jp1];  
                    m = MEXOPCol[jp1];  

                    if (opt) {
                        tmp = MEXOP[jp2][1];
                        r = n * m;
                        for (i = 1; i <= r; ++i)
                            MEXOP[jp1][i] += tmp;
                    }
                    jp--;
                }
                else {
                    if (check_dim_add(jp1,jp2))  
                        return(2);
     
                    if (opt) {
                        n = MEXOPRow[jp1] * MEXOPCol[jp1];  
                        for (i = 1; i <= n; ++i)
                            MEXOP[jp1][i] += MEXOP[jp2][i];
                    }
                    jp--;
                }
            }
            else if (atyp == 2) {                           /* subtraction */

                if (check_scalar(jp1)) {

                    n = MEXOPRow[jp2];  
                    m = MEXOPCol[jp2];  

                    if (opt) {
                        tmp = MEXOP[jp1][1];
                        r = n * m;
                        for (i = 1; i <= r; ++i)
                            MEXOP[jp1][i] = tmp - MEXOP[jp2][i];
                    }
                    MEXOPRow[jp1] = n;              
                    MEXOPCol[jp1] = m;              
                    jp--;
                }
                else if (check_scalar(jp2)) {

                    n = MEXOPRow[jp1];  
                    m = MEXOPCol[jp1];  

                    if (opt) {
                        tmp = MEXOP[jp2][1];
                        r = n * m;
                        for (i = 1; i <= r; ++i)
                            MEXOP[jp1][i] -= tmp;
                    }
                    jp--;
                }
                else {
                    if (check_dim_add(jp1,jp2))  
                        return(2);
     
                    if (opt) {
                        n = MEXOPRow[jp1] * MEXOPCol[jp1];  
                        for (i = 1; i <= n; ++i)
                            MEXOP[jp1][i] -= MEXOP[jp2][i];
                    }
                    jp--;
                }
            }
            else if (atyp == 3) {                           /* multiplication */

                if (check_scalar(jp2)) {
                    if (opt) {
                        n = MEXOPRow[jp1] * MEXOPCol[jp1];
                        tmp = MEXOP[jp2][1];
                        for (i = 1; i <= n; ++i)
                            MEXOP[jp1][i] *= tmp;           
                    }
                    jp--;
                }
                else if (check_scalar(jp1)) {
                    if (opt) {
                        n = MEXOPRow[jp2] * MEXOPCol[jp2];
                        tmp = MEXOP[jp1][1];
                        for (i = 1; i <= n; ++i)
                            MEXOP[jp1][i] = tmp * MEXOP[jp2][i];
                    }
                    MEXOPRow[jp1] = MEXOPRow[jp2];
                    MEXOPCol[jp1] = MEXOPCol[jp2];
                    jp--;
                }
                else {
                    if (check_dim_mul(jp1,jp2))  
                        return(2);
 
                    MEXOPMaxL = imax(MEXOPMaxL,jp + 1);
                    n = MEXOPRow[jp1];  
                    m = MEXOPCol[jp2];  

                    if (opt) {
                        if (jp >= MEXOPMax)  
                            return(1);

                        r = MEXOPCol[jp1];  
                        for (i = 0; i < n; ++i) {
                            for (j = 1; j <= m; ++j) {
                                tmp = 0.0;
                                for (k = 0; k < r; ++k) 
                                    tmp += MEXOP[jp1][i * r + k + 1] *
                                           MEXOP[jp2][k * m + j];
                                
                                MEXOP[jp][i * m + j] = tmp; 
                            }
                        }
                        r = n * m;
                        for (i = 1; i <= r; ++i)
                            MEXOP[jp1][i] = MEXOP[jp][i];
                    }
                    MEXOPRow[jp1] = n;              
                    MEXOPCol[jp1] = m;              
                    MEXOPMaxM = imax(MEXOPMaxM,n * m);
                    jp--;
                }
            }
            else if (atyp == 5) {                           /* division */

                if (check_scalar(jp2) == 0)  
                    return(6);

                if (opt) {
                    tmp = MEXOP[jp2][1];
                    if (tmp == 0.0)
                        return(7);
                    r = MEXOPRow[jp1] * MEXOPCol[jp1];
                    for (i = 1; i <= r; ++i)
                        MEXOP[jp1][i] /= tmp;           
                }
                jp--;
            }
            else if (atyp == 6) {                        /* ^ */
                if (check_scalar(jp2) == 0)  
                    return(5);

                if (opt) {
                    tmp1 = MEXOP[jp2][1];
                    r = MEXOPRow[jp1] * MEXOPCol[jp1];
                    for (i = 1; i <= r; ++i) {
                        tmp = MEXOP[jp1][i];
                        if (fabs(tmp) <= EPSI2)
                            MEXOP[jp1][i] = 0.0;
                        else if (tmp > 0.0)
                            MEXOP[jp1][i] = pow(tmp,tmp1);
                        else {
                            if (floor(tmp1) != tmp1)
                                return(5);
                            MEXOP[jp1][i] = pow(tmp,tmp1);
                        }
                    }   
                }
                jp--;
            }

            else                    /* undefined */
                gerr_exit(4);

            if (utyp < 0) {                         /* change sign */
                r = MEXOPRow[jp1] * MEXOPCol[jp1];
                for (i = 1; i <= r; ++i)
                    MEXOP[jp1][i] = -MEXOP[jp1][i];
            }
        }
        MEXOPMaxL = imax(MEXOPMaxL,jp);
    } 
    if (jp != 1)           
        return(-99);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  check_dim_add(i,j)      Return 0 if OK, -1 if dimension error           */

int check_dim_add(int i,int j)
{
    if (MEXOPRow[i] != MEXOPRow[j] || MEXOPCol[i] != MEXOPCol[j])       
        return(-1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  check_dim_mul(i,j)      Return 0 if OK, -1 if dimension error           */

int check_dim_mul(int i,int j)
{
    if (MEXOPCol[i] != MEXOPRow[j])       
        return(-1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  check_scalar(i)     Return 1 if scalar, otherwise 0.                    */

int check_scalar(int i)
{
    if (MEXOPCol[i] == 1 && MEXOPRow[i] == 1)       
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  check_square(i)     Return 1 if square, otherwise 0.                    */

int check_square(int i)
{
    if (MEXOPCol[i] == MEXOPRow[i])       
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

int alloc_mex_eval(int opt,int n,int m)
{
    int i;

    if (opt) {
        if (n >= MEXOPMax) {
            printf1("Error: exceeding max stack size (%d).\n",MEXOPMax);
            return(-1);
        }
        MEXOPMaxB = m;

        for (i = 0; i < n; ++i) {
            if (!(MEXOP[i] = (double *)calloc(m + 1,sizeof(double)))) { 
                p_err(-2,1);  
                goto AMFin;              
            }
            memrq(m + 1,sizeof(double));
            MEXOPMaxA++;
        }
        return(0);
    }

AMFin:
    for (i = 0; i < MEXOPMaxA; ++i) {
        free((char *)MEXOP[i]);
        memrq(-MEXOPMaxB - 1,sizeof(double));
    }
    MEXOPMaxA = MEXOPMaxB = 0;
    if (opt)
        return(-1);
    else
        return(0);
}

/* ------------------------------------------------------------------------ */
/*  mex_error(n)    Report error n.                                         */

void mex_error(int n)
{
    printf1("\nError (%d) in evaluation of matrix expression.\n",n);
    printf1("Error type: ");
    switch (n) {
        case -1:    printf1("stack overflow.\n");
                    break;
        case -2:    printf1("unknown operator or expression, or misplaced brackets.\n");
                    break;
        case -3:    printf1("syntax error.\n");
                    break;
        case -6:    printf1("scan error floating point number.\n");
                    break;
        case -7:    printf1("misplaced minus sign.\n");
                    break;
        case -10:   printf1("undefined matrix or variable name.\n");
                    break;
        case -25:   printf1("misplaced brackets.\n");
                    break;

        case  1:    printf1("out of stack space.\n");
                    break;
        case  2:    printf1("error in matrix dimensions.\n");
                    break;
        case  3:    printf1("cannot evaluate sqrt.\n");
                    break;
        case  4:    printf1("trace requires square matrix.\n");
                    break;
        case  5:    printf1("cannot evaluate power function.\n");
                    break;
        case  6:    printf1("division requires scalar.\n");
                    break;
        case  7:    printf1("division by zero.\n");
                    break;
        case  8:    printf1("exp out of range.\n");
                    break;
        case  9:    printf1("cannot calculate generalized inverse.\n");
                    break;
        default:    printf1("unknown.\n");
                    break;
    }
}

/*--##----------------------------------------------------------------------*/
/*  mex_eval(s,row,col)     Evaluate matrix expression s. Put Result        */
/*                          into MEXOP[0], NEXOPRow[0], MEXOPCol[0].        */  
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int mex_eval(char *s)
{
    int err,r;

    err = -1;
    if (alloc_mex(1000,1))
        goto MEVFin;

    if ((r = mex_parse(s)) || MEXCnt == 0) {
        mex_error(r); 
        goto MEVFin;
    }
    if ((r = m_eval(0,MEXCnt,MEXTyp,MEXVal))) {                       
        mex_error(r); 
        goto MEVFin;
    }
    if (alloc_mex_eval(1,MEXOPMaxL,MEXOPMaxM))
        goto MEVFin;
             
    if ((r = m_eval(1,MEXCnt,MEXTyp,MEXVal))) {                       
        mex_error(r); 
        goto MEVFin;
    }
    err = 0;

MEVFin:       
    return(err);
}

