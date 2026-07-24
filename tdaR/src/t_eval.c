/****************************************************************************/
/*  t_eval                                                                  */
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
#include "t_cdf.h"
#include "t_rzoo.h"
#include "t_edat.h"
#include "t_parm.h"
#include "t_smo.h"
#include "t_alloc.h"
#include "t_gmin.h"
#include "t_int.h"
#include "t_mat.h"
#include "t_matf.h"
#include "t_gcmd.h"
#include "tda_context.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_eval.c                                                   */

int alloc_est(TDAContext *ctx, int n,int opt);
void parse(TDAContext *ctx, short cnt,int *typ,double *val,int *idx);
int v_parse(TDAContext *ctx, char *s,int ivflg);
void v_expr(TDAContext *ctx);
void v_expr1(TDAContext *ctx);
void v_term(TDAContext *ctx);
int v_getindex(TDAContext *ctx, int idx);
void v_factor(TDAContext *ctx);
void v_match(TDAContext *ctx, int t);
int v_lexan(TDAContext *ctx);             
char *skip_const(TDAContext *ctx, char *p);                   
int v_search(TDAContext *ctx, char *p,int *len,int *narg); 
int v_search1(TDAContext *ctx, char *p,int *len);
void prn_emsg1(TDAContext *ctx, int n);
void prn_emsg2(TDAContext *ctx, int n);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

#define NUM 257         /* return code for floating point numbers           */
#define IDV 258         /* return code for V, C or X variables              */
#define IDM 259         /* return code for matrices                         */
#define IDG 260         /* return code for operators                        */
#define IDC 261         /* return code for constants, operators without     */
                        /* arguments.                                       */
#define IDNL 262        /* return code for namelist                         */

#define DDTYP     400   /* type number for dummy operators                  */
#define CCTYP      99   /* type number for ,, in dummy operators            */
#define COLTYP     98   /* type number for :  in dummy operators            */
#define RVECTYP  1698   /* type number for row vectors <...>'               */
#define CVECTYP  1699   /* type number for column vectors <...>             */


                        /* must be the highest number                       */  




/*  ESErr =  error flag, see v_parse.                                       */



                            /* corresponding to time operator               */



/*  the stack is organized as follows:                                      */
/*  ESTyp =  0  floating point number (value in ESVal)                      */
/*           1  +                                                           */
/*           2  -                                                           */
/*           3  *                                                           */
/*           4  &                                                           */
/*           5  /                                                           */
/*           6  ^                                                           */
/*           7  %                                                           */
/*           8  |                                                           */
/*          20  ... operators, see OPKEY                                    */
/*              in this cases, ESVal contains number of arguments           */
/*                                                                          */
/*   NLOFFS+ i  i.th namelist                                               */
/*   MOFFS + i  i.th matrix                                                 */
/*   VOFFS + i  i.th V variable                                             */
/*   COFFS + i  i.th C variable                                             */
/*                                                                          */


/*  array with keywords for operators                                       */

struct OPKEY {
    char *op;           /* keyword for operator                             */
    short typ;          /* type number of operator                          */
    short len;          /* length of keyword                                */
    short narg;         /* number of arguments (0 if variable)              */
} OPS[] = {
    { "exp",    100,  3,  1 },
    { "log",    101,  3,  1 },
    { "rnd",    102,  3,  1 },
    { "abs",    103,  3,  1 },
    { "sign",   104,  4,  1 },
    { "floor",  105,  5,  1 },
    { "ceil",   106,  4,  1 },
    { "tr",     107,  2,  3 },
    { "min",    108,  3,  0 },
    { "max",    109,  3,  0 },
    { "sqrt",   110,  4,  1 },
    { "sin",    111,  3,  1 },
    { "cos",    112,  3,  1 },
    { "lgam",   113,  4,  1 },
    { "icg1",   114,  4,  2 },
    { "icb",    115,  3,  3 },
    { "eexp",   116,  4,  1 },
    { "digam",  117,  5,  1 },
    { "trigam", 118,  6,  1 },
    { "iv",     119,  2,  2 },
    { "icg",    120,  3,  2 },
    { "bc",     121,  2,  2 },
    { "bivn",   122,  4,  3 },
    { "mvn",    123,  3,  0 },
    { "poisson",124,  7,  2 },
    { "negbin", 125,  6,  3 },
    { "gclen",  126,  5,  4 },

    { "intn",   186,  4,  2 },      /* 184 - 186 reserved */
    { "inth",   189,  4,  2 },      /* 187 - 189 reserved */
    { "int",    193,  3,  3 },      /* 190 - 192 reserved */

                                /* 196 reserved for X(.,.) */
                                /* 197 reserved for X(.,j) */
                                /* 198 reserved for X(i,.) */
                                /* 199 reserved for X(i,j) */
               
    { "not",    200,  3,  1 },
    { "and",    201,  3,  2 },
    { "or",     202,  2,  2 },
    { "eq",     203,  2,  2 },
    { "ne",     204,  2,  2 },
    { "lt",     205,  2,  2 },
    { "le",     206,  2,  2 },
    { "gt",     207,  2,  2 },
    { "ge",     208,  2,  2 },
    { "if",     209,  2,  3 },

    { "ndf",    300,  3,  1 },
    { "ndi",    302,  3,  1 },
    { "nd",     301,  2,  1 },
    { "mr",     303,  2,  1 },
    { "td",     304,  2,  2 },
    { "chd",    305,  3,  2 },
    { "fd",     306,  2,  3 },
                                /* 444 reserved for [...] */
    { "jul",    402,  3,  3 },
    { "jul1",   403,  4,  3 },
    { "rd",     404,  2,  2 },
    { "rdmn",   405,  4,  2 },
    { "grd",    406,  3,  1 },
    { "july",   409,  4,  1 },
    { "julm",   410,  4,  1 },
    { "juld",   411,  4,  1 },
    { "num",    412,  3,  2 },
    { "recode", 413,  6,  2 },
    { "recode1",414,  7,  3 },
    { "row",    415,  3,  1 },
    { "col",    416,  3,  1 },      /* 417 reserved for begin of row/col expr */
    { "exists", 418,  6,  0 },
    { "rdp1",   419,  4,  1 },

    { "str",    450,  3,  2 },      /* string variable operator */
    { "strlen", 451,  6,  1 },      /* strlen(S) */
    { "strvp",  452,  5,  3 },      /* strvp(S,n,m) */
    { "strv",   453,  4,  1 },      /* strv(S) */

                                /* operators for sequences */
    { "slen",   501,  4, -1 },
    { "glen",   502,  4, -1 },

                                /* operators for spatial data */
    { "sdp",    601,  3,  2 },

    { "lag",   1500,  3,  2 },      /* begin of type 2 operators */
    { "sum",   1501,  3,  1 },
    { "mean",  1502,  4,  1 },
    { "std",   1503,  3,  1 },
    { "cum",   1504,  3,  1 },
    { "cd",    1505,  2,  1 },
    { "cdv",   1506,  3,  1 },
    { "vmin",  1507,  4,  1 },
    { "vmax",  1508,  4,  1 },
    { "suc",   1509,  3,  1 },
    { "pre",   1510,  3,  1 },
    { "mav",   1511,  3,  2 },
    { "vdif",  1512,  4,  2 },
    { "cmean", 1513,  5,  2 },
    { "cmean1",1514,  6,  2 },

    { "sort",  1520,  4,  1 },
    { "snum",  1521,  4,  1 },
    { "rank",  1522,  4,  1 },
    { "ndv",   1523,  3,  1 },
    { "ndv1",  1524,  4,  2 },
    { "ndv2",  1525,  4,  3 },
    { "strsp", 1526,  5,  1 },
    { "quant", 1527,  5,  2 },
    { "quant1",1528,  6,  4 },

    { "cnteq", 1530,  5,  2 },
    { "cntne", 1531,  5,  2 },
    { "cntlt", 1532,  5,  2 },
    { "cntgt", 1533,  5,  2 },
    { "cntle", 1534,  5,  2 },
    { "cntge", 1535,  5,  2 },

    { "gcnt",  1540,  4,  1 },
    { "grec",  1541,  4,  1 },
    { "gsn",   1542,  3,  1 },
    { "gfirst",1543,  6,  1 },
    { "glast", 1544,  5,  1 },

    { "gsum",  1550,  4,  2 },
    { "gmean", 1551,  5,  2 },
    { "gstd",  1552,  4,  2 },
    { "gmin",  1553,  4,  2 },
    { "gmax",  1554,  4,  2 },
    { "gsort", 1555,  5,  2 },
    { "gndv",  1556,  4,  2 },
    { "gndv1", 1557,  5,  3 },

    { "change",1600,  6,  1 },
    { "cntch", 1601,  5,  1 },
    { "ccntch",1602,  6,  1 },
    { "lagch", 1603,  5,  1 },
                            /* 1697 reserved for X(i) */
                            /* 1698, 1699 reserved for VECTYP */

    { "sma",   1700, -3,  1 },  /* special type 2 operators */
    { "smd",   1701, -3,  1 },

    { "bmin",  1880,  4,  1 },  /* block mode operators */
    { "bmax",  1881,  4,  1 },
    { "bfa",   1882,  3,  1 },
    { "bfr",   1883,  3,  1 },

    { "bnrec", 1890,  5,  0 },  /* here without arguments */
    { "brec",  1891,  4,  0 },
    { "bfirst",1892,  6,  0 },
    { "blast", 1893,  5,  0 },
    { "bnum",  1894,  4,  0 },
    { "brd",   1899,  3,  0 },
                            /* end of type 2 operators */

    { "case",  1900,  4,  0 },
    { "nvar",  1901,  4,  0 },
    { "pi",    1902,  2,  0 },
    { "rdn1",  1905,  4,  0 },
    { "rdn",   1904,  3,  0 },
    { "rd",    1903,  2,  0 },
    { "nocdm", 1906,  5,  0 },
    { "noc",   1907,  3,  0 },
    { "bnoc",  1908,  4,  0 },

    { "sdxmin",1910,  6,  0 },
    { "sdxmax",1911,  6,  0 },
    { "sdymin",1912,  6,  0 },
    { "sdymax",1913,  6,  0 },
    { "sdarea",1914,  6,  0 },
    { "sdlen", 1915,  5,  0 },

    { "t",     1969,  1,  0 },  /* special for integration */

                            /* begin of type 3 operators */
    { "time",  1970,  4,  0 },  /* special for Cox models */
    { "sn",    1980,  2,  0 },  /* special for episode data */
    { "org",   1981,  3,  0 },
    { "des",   1982,  3,  0 },
    { "ts",    1983,  2,  0 },
    { "tf",    1984,  2,  0 },
                            /* end of type 3 operators */
    { "",         0,  0,  0 },
};

/*  NOTE: type IAOFFS + i is used for intermediate function arguments, and
    FAOFFS + i is used for function arguments. */ 


/*  function definition */

                        /* FVFlg = 2 if using episode data                  */









/* this flag is zero. */

/* ------------------------------------------------------------------------ */
/*  alloc_est(n,opt)                                                        */
/*                                                                          */
/*  Allocate (opt = 1) or free (opt = 0) ESTyp and ESVal with n elements.   */
/*  Return 0 if OK, -1 if error.                                            */

int alloc_est(TDAContext *ctx, int n,int opt)
{

    if (opt == 1) {
        if (ctx->s_alloc_est_memcnt) {
            gerr_exit(ctx, 3);
        }
        if (!(ctx->ESTyp = (int *)calloc((size_t)(n + 2),sizeof(int)))) { 
            p_err(ctx, -2,1);  
            return(-1);
        }
        if (!(ctx->ESIdx = (int *)calloc((size_t)(n + 2),sizeof(int)))) { 
            p_err(ctx, -2,1);  
            return(-1);
        }
        if (!(ctx->ESVal = (double *)calloc((size_t)(n + 2),sizeof(double)))) { 
            p_err(ctx, -2,1);    
            return(-1);
        }
        ctx->ESMaxLen = n;
        ctx->s_alloc_est_memcnt = (int)((size_t)(n + 2) *
                                    (2 * sizeof(int) + sizeof(double)));
        memrq(ctx, ctx->s_alloc_est_memcnt,1);
        return(0);
    }
    else {
        free((char *)ctx->ESTyp);
        free((char *)ctx->ESIdx);
        free((char *)ctx->ESVal);
        memrq(ctx, -ctx->s_alloc_est_memcnt,1);
        ctx->ESMaxLen = ctx->s_alloc_est_memcnt = 0;
        return(0);
    }
}

/*--------------------------------------------------------------------------*/
/*  parse(cnt,typ,val,idx)  show parser stack on standard output.           */

void parse(TDAContext *ctx, short cnt,int *typ,double *val,int *idx)
{
    register int i;
    int t,ta;

    printf1(ctx, "\n Cnt         Typ         Val         Idx\n");
    prnchar(ctx, '-',40,1);

    for (i = 0; i < cnt; ++i) {
        t = typ[i];
        ta = iabs(ctx, t);
        printf1(ctx, "%4d  %10d  %10.1f  %10d  ",i,t,val[i],idx[i]);
        if (ta >= ctx->IFTYP)
            printf1(ctx, "special\n");
        else if (ta >= ctx->COFFS)
            printf1(ctx, "c%d\n",ta - ctx->COFFS);
        else if (ta >= ctx->VOFFS)
            printf1(ctx, "%s\n",ctx->VName[ta - ctx->VOFFS]);
        else if (ta >= ctx->MOFFS)
            printf1(ctx, "%s\n",ctx->MatName[ta - ctx->MOFFS]);
        else
            printf1(ctx, "\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  r = v_parse(string,ivflg)                                               */
/*                                                                          */
/*  Transforms string to postfix notation. Results are put in the arrays    */
/*  ESTyp[i] and ESVal[i]. Max stack size is assumed to be ESMaxLen. After  */
/*  parsing, the number of entries in the stack is put into ESCnt.          */
/*                                                                          */
/*  if ivflg != 0 interval operators are allowd, otherwise not.             */
/*  In any case, the following global variable is set:                      */
/*                                                                          */
/*  IVEXPRFlg = 1  if expression contains intervals                         */
/*              0  otherwise                                                */
/*                                                                          */
/*  Return: 0 if ok, otherwise ESErr.                                       */
/*  ESErr =  0  no error                                                    */
/*           1  no error, using data matrix variables                       */
/*          -1  stack overflow                                              */
/*          -2  if part of input stream cannot be interpreted               */
/*          -3  match error, syntax error                                   */
/*          -4  undefined operator                                          */
/*          -6  scan error for floating point numbers                       */
/*          -7  misplaced minus sign                                        */
/*          -8  undefined V variable                                        */
/*          -9  Xi variable with i <= 0                                     */
/*         -10  undefined name of matrix/var/namelist                       */  
/*         -11  if in pre(X), X is not a variable name                      */
/*         -12  error in ci reference                                       */
/*         -13  cannot use matrix reference                                 */
/*         -14  need reference to a matrix                                  */
/*         -15  matrix should have at least two columns                     */
/*         -16  integral must be first term in expression                   */
/*         -17  t operator may only occur inside integral                   */
/*         -18  if error in sequence definition                             */
/*         -19  if error in matrix dimension for matrix expressions         */
/*         -20  cannot use namelist                                         */
/*         -21  syntax error in matrix indices                              */
/*         -22  cannot use vector expression                                */
/*         -23  reference to undefined local matrix                         */
/*         -24  unvalid reference to macro argument                         */
/*         -25  misplaced bracket                                           */
/*         -26  expression mmust not contain intervals                      */
/*         -27  cannot combine str with other operators                     */
/*                                                                          */
/*  Note: basic ideas for this function are taken from                      */
/*  A.V. Aho, Compilers. Principles, Techniques and Tools, 1987.            */

int v_parse(TDAContext *ctx, char *s,int ivflg)
{
    register int i,iflg,t;

    ctx->IVEXPRFlg = 0;  
    ctx->MatExprCnt = 0;
    ctx->MatExprRow = ctx->MatExprCol = 1;
    ctx->ESP = s;
    ctx->EVFlag = ctx->ESTLMin = ctx->ESErr = ctx->ESCnt = 0;
    ctx->ESTLLast = 1;
    ctx->EST0Last = 0;
    ctx->IFLevel = 0;

    ctx->ESNxt = v_lexan(ctx);
    if (ctx->ESErr)
        return(ctx->ESErr);
    while (ctx->ESNxt) {
        v_expr(ctx);
        if (ctx->ESErr)
            return(ctx->ESErr);
    }
    if (ivflg == 0 && ctx->IVEXPRFlg)      
        ctx->ESErr = -26;

    if (ctx->ESErr == 0) {

        /* check for misplaced t operator (integration) */

        iflg = 0;
        for (i = 0; i < ctx->ESCnt; ++i) {
            t = ctx->ESTyp[i];
            if (t == 184 || t == 187 || t == 190)
                iflg = 1;
            else if (t == 186 || t == 189 || t == 193)
                iflg = 0;
            else if (t == 1969) {
                if (iflg == 0)
                    return(-17);
            }
        }   
        ctx->ESErr = ctx->EVFlag;
    }
    if (ctx->ESErr == 0) {       /* check for str-operator */
        iflg = 0;
        for (i = 0; i < ctx->ESCnt; ++i) {
            if (ctx->ESTyp[i] == 450) {
                iflg = 1;
                break;
            }
        }
        if (iflg && ctx->ESCnt != 3)
            ctx->ESErr = -27;
    }
    return(ctx->ESErr);
}

/* ------------------------------------------------------------------------ */
/*  expr                                                                    */

void v_expr(TDAContext *ctx)
{
    if (ctx->ESCnt >= ctx->ESMaxLen) {
        ctx->ESErr = -1;
        return;
    }
    v_expr1(ctx);
    if (ctx->ESErr)
        return;

    while (1) {
        if (ctx->ESCnt >= ctx->ESMaxLen) {
            ctx->ESErr = -1;
            break;  
        }
        switch (ctx->ESNxt) {
            case '+':   v_match(ctx, ctx->ESNxt);
                        if (ctx->ESErr)
                            return;
                        v_expr1(ctx);
                        if (ctx->ESErr)
                            return;
                        ctx->ESIdx[ctx->ESCnt] = 0;
                        ctx->ESVal[ctx->ESCnt] = 2.0;              
                        ctx->ESTyp[ctx->ESCnt++] = 1;
                        break;
            case '-':   v_match(ctx, ctx->ESNxt);
                        if (ctx->ESErr)
                            return;
                        v_expr1(ctx);
                        if (ctx->ESErr)
                            return;
                        ctx->ESIdx[ctx->ESCnt] = 0;
                        ctx->ESVal[ctx->ESCnt] = 2.0;              
                        ctx->ESTyp[ctx->ESCnt++] = 2;
                        break;

            case '|':   v_match(ctx, ctx->ESNxt);
                        if (ctx->ESErr)
                            return;
                        v_expr1(ctx);
                        if (ctx->ESErr)
                            return;
                        ctx->ESIdx[ctx->ESCnt] = 0;
                        ctx->ESVal[ctx->ESCnt] = 2.0;              
                        ctx->ESTyp[ctx->ESCnt++] = 8;
                        break;

            default:    return;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  expr1                                                                   */

void v_expr1(TDAContext *ctx)
{
    if (ctx->ESCnt >= ctx->ESMaxLen) {
        ctx->ESErr = -1;
        return;
    }
    v_term(ctx);
    if (ctx->ESErr)
        return;

    while (1) {
        if (ctx->ESCnt >= ctx->ESMaxLen) {
            ctx->ESErr = -1;
            break;  
        }
        switch (ctx->ESNxt) {
            case '*':   v_match(ctx, ctx->ESNxt);
                        if (ctx->ESErr)
                            return;
                        v_term(ctx);
                        if (ctx->ESErr)
                            return;
                        ctx->ESIdx[ctx->ESCnt] = 0;
                        ctx->ESVal[ctx->ESCnt] = 2.0;              
                        ctx->ESTyp[ctx->ESCnt++] = 3;
                        break;
            case '&':   v_match(ctx, ctx->ESNxt);
                        if (ctx->ESErr)
                            return;
                        v_term(ctx);
                        if (ctx->ESErr)
                            return;
                        ctx->ESIdx[ctx->ESCnt] = 0;
                        ctx->ESVal[ctx->ESCnt] = 2.0;              
                        ctx->ESTyp[ctx->ESCnt++] = 4;
                        break;
            case '/':   v_match(ctx, ctx->ESNxt);
                        if (ctx->ESErr)
                            return;
                        v_term(ctx);
                        if (ctx->ESErr)
                            return;
                        ctx->ESIdx[ctx->ESCnt] = 0;
                        ctx->ESVal[ctx->ESCnt] = 2.0;              
                        ctx->ESTyp[ctx->ESCnt++] = 5;
                        break;
            default:    return;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  term                                                                    */

void v_term(TDAContext *ctx)
{
    int nn;

    if (ctx->ESCnt >= ctx->ESMaxLen) {
        ctx->ESErr = -1;
        return;
    }
    v_factor(ctx);
    if (ctx->ESErr)
        return;

    while (1) {
        if (ctx->ESCnt >= ctx->ESMaxLen) {
            ctx->ESErr = -1;
            break;  
        }
        switch (ctx->ESNxt) {

            case '^':   v_match(ctx, ctx->ESNxt);
                        if (ctx->ESErr)
                            return;
                        v_factor(ctx);
                        if (ctx->ESErr)
                            return;
                        ctx->ESIdx[ctx->ESCnt] = 0;
                        ctx->ESVal[ctx->ESCnt] = 2.0;              
                        ctx->ESTyp[ctx->ESCnt++] = 6;
                        break;
            case '%':
                        v_match(ctx, ctx->ESNxt);
                        if (ctx->ESErr)
                            return;
                        v_factor(ctx);
                        if (ctx->ESErr)
                            return;
                        ctx->ESIdx[ctx->ESCnt] = 0;
                        ctx->ESVal[ctx->ESCnt] = 2.0;              
                        ctx->ESTyp[ctx->ESCnt++] = 7;
                        break;
            case '[': 
                        v_match(ctx, '[');
                        if (ctx->ESErr)
                            return;
                        v_expr(ctx);
                        if (ctx->ESErr)
                            return;
                        nn = 1;
                        while (1) {
                            if (ctx->ESNxt == ',') {
                                nn++;
                                v_match(ctx, ',');
                                if (ctx->ESErr)   
                                    return;
                                if (ctx->ESNxt == ',') {
                                    ctx->ESIdx[ctx->ESCnt] = 0;
                                    ctx->ESTyp[ctx->ESCnt] = CCTYP;
                                    ctx->ESVal[ctx->ESCnt++] = 0.0;
                                    nn++;
                                    v_match(ctx, ',');
                                    if (ctx->ESErr)
                                        return;
                                }  
                                v_expr(ctx);
                            }
                            else if (ctx->ESNxt == ':') {
                                nn++;
                                v_match(ctx, ':');
                                if (ctx->ESErr)   
                                    return;
                                ctx->ESIdx[ctx->ESCnt] = 0;
                                ctx->ESTyp[ctx->ESCnt] = COLTYP;
                                ctx->ESVal[ctx->ESCnt++] = 0.0;
                                nn++;
                                v_expr(ctx);
                            }
                            else
                                break;
                            if (ctx->ESErr)
                                return;
                        }
                        v_match(ctx, ']');
                        ctx->ESIdx[ctx->ESCnt] = 0;
                        ctx->ESTyp[ctx->ESCnt] = DDTYP;
                        ctx->ESVal[ctx->ESCnt++] = (double)nn;
                        break;

            default:    return;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  v_getindex(idx)         Get matrix, var or namelist indices.            */
/*                          199   X(i,j)                                    */
/*                          1697  X(i)                                      */
/*                          198   X(i,.)                                    */
/*                          197   X(.,j)                                    */
/*                          196   X(.,.)                                    */
/*                                                                          */
/*  Return 1 if full row dimension, otherwise 0.                            */

int v_getindex(TDAContext *ctx, int idx)
{
    int t = 199;
    int na = 3;
    int r = 0;

    if (*ctx->ESP == '.') {
        t = 197;
        na = 2;
        ctx->ESP++;
        v_match(ctx, ctx->ESNxt);
        if (ctx->ESNxt != ',') {
            ctx->ESErr = -21;
            return(r);
        }
        goto GICONT0;
    }
    v_match(ctx, ctx->ESNxt);
    if (ctx->ESErr)    
        return(r);
    v_expr(ctx);
    if (ctx->ESErr)   
        return(r);
GICONT0:
    if (ctx->ESNxt != ',') {
        if (ctx->ESNxt == ')') {
            t = 1697;
            na = 2;
            r = 1;
            goto GICONT1;
        }
        ctx->ESErr = -21;
        return(r);
    }
    if (*ctx->ESP == '.') {
        if (t == 197) {
            t = 196;
            na = 1;
        }
        else {
            t = 198;
            na = 2;
        }
        ctx->ESP++;
        v_match(ctx, ctx->ESNxt);
        goto GICONT1;
    }
    v_match(ctx, ctx->ESNxt);
    if (ctx->ESErr)   
        return(r);
    v_expr(ctx);
    if (ctx->ESErr)   
        return(r);
GICONT1:
    v_match(ctx, ')');
    if (ctx->ESErr) { 
        ctx->ESErr = -21;
        return(r);
    }
    if (ctx->ESCnt >= ctx->ESMaxLen) {
        ctx->ESErr = -1;
        return(r);
    }
    ctx->ESVal[ctx->ESCnt] = (double)na;       
    ctx->ESTyp[ctx->ESCnt] = t;
    ctx->ESIdx[ctx->ESCnt] = idx;             
    ctx->ESCnt++;
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  factor                                                                  */
/*                                                                          */
/*  ESErr = -11  if in pre(X), X is not a variable name                     */
/*                                                                          */
/*          -13  cannot use matrix reference                                */
/*          -14  need reference to a matrix                                 */
/*          -15  matrix should have at least two columns                    */
/*          -16  integral must be first term in expression                  */
/*          -18  if error in sequence definition                            */
/*          -19  if error in matrix dimension for matrix expressions        */
/*          -20  cannot use namelist                                        */
/*          -21  syntax error in matrix indices                             */
/*          -22  cannot use vector expression                               */
/*          -25  misplaced bracket                                          */
    
void v_factor(TDAContext *ctx)
{
    int i,len,tv,tva,tvaa,na,nn,v1,v2;

    tva = -1;

    if (ctx->ESCnt >= ctx->ESMaxLen) {
        ctx->ESErr = -1;
        return;
    }
    switch (ctx->ESNxt) {
        case '(': 
                    v_match(ctx, '(');
                    if (ctx->ESErr)
                        return;
                    v_expr(ctx);
                    if (ctx->ESErr)
                        return;
                    v_match(ctx, ')');
                    if (ctx->ESNxt == '(')  
                        ctx->ESErr = -25;
                    break;

        case '[':                       /* interval */
                    ctx->IVEXPRFlg = 1;  
                    v_match(ctx, '[');
                    if (ctx->ESErr)
                        return;
                    v_expr(ctx);
                    if (ctx->ESErr)
                        return;
                    v_match(ctx, ',');
                    v_expr(ctx);
                    if (ctx->ESErr)
                        return;
                    v_match(ctx, ']');
                    if (ctx->ESNxt == '(')  
                        ctx->ESErr = -25;

                    ctx->ESIdx[ctx->ESCnt] = 0;
                    ctx->ESTyp[ctx->ESCnt] = 119;
                    ctx->ESVal[ctx->ESCnt] = (double)2;
                    ctx->ESCnt++;
                    break;

        case NUM:   ctx->ESTyp[ctx->ESCnt] = 0;
                    ctx->ESIdx[ctx->ESCnt] = 0;
                    ctx->ESVal[ctx->ESCnt++] = ctx->ESTVal;
                    v_match(ctx, NUM);
                    if (ctx->ESNxt == '(')  
                        ctx->ESErr = -25;
                    break;
        case IDV: 
                    ctx->ESTyp[ctx->ESCnt] = (int)ctx->ESTVal;
                    ctx->ESIdx[ctx->ESCnt] = nn = iabs(ctx, ctx->ESTyp[ctx->ESCnt]);
                    ctx->ESVal[ctx->ESCnt] = (double)(iabs(ctx, (int)ctx->ESTVal) - ctx->VOFFS);
                    ctx->ESCnt++;
                    v_match(ctx, IDV);

                    if (ctx->MatExprFlg) {
                        if (ctx->ESNxt == '(') {
                            if (v_getindex(ctx, ctx->ESTyp[ctx->ESCnt - 1])) {
                                if (ctx->MatExprCnt == 0)  
                                    ctx->MatExprRow = imax(ctx, ctx->MatExprRow,ctx->NOC);
                            }
                        }
                        else if (ctx->MatExprCnt == 0) {
                            ctx->MatExprRow = imax(ctx, ctx->MatExprRow,ctx->NOC);
                            ctx->MatExprCol = imax(ctx, ctx->MatExprCol,1);
                        }
                    }
                    break;
        case IDM:                   /* ## matrix */

                    if (ctx->MatExprFlg == 0) {
                        ctx->ESErr = -13;
                        return;
                    }
                    nn = iabs(ctx, (int)ctx->ESTVal) - ctx->MOFFS;

                    ctx->ESTyp[ctx->ESCnt] = (int)ctx->ESTVal;
                    ctx->ESIdx[ctx->ESCnt] = iabs(ctx, ctx->ESTyp[ctx->ESCnt]);
                    ctx->ESVal[ctx->ESCnt] = (double)(iabs(ctx, (int)ctx->ESTVal) - ctx->MOFFS);
                    ctx->ESCnt++;
                    v_match(ctx, IDM);

                    if (ctx->ESNxt == '(') {
                        if (v_getindex(ctx, ctx->ESTyp[ctx->ESCnt - 1])) {
                            if (ctx->MatExprCnt == 0)  
                                ctx->MatExprRow = imax(ctx, ctx->MatExprRow,ctx->MatRow[nn]);
                        }
                    }
                    else if (ctx->MatExprCnt == 0) {
                        ctx->MatExprRow = imax(ctx, ctx->MatExprRow,ctx->MatRow[nn]);
                        ctx->MatExprCol = imax(ctx, ctx->MatExprCol,ctx->MatCol[nn]);
                    }
                    return;

        case IDNL:                  /* ## namelist */

                    if (ctx->MatExprFlg == 0 || ctx->NOC == 0) {
                        ctx->ESErr = -20;
                        return;
                    }
                    nn = iabs(ctx, (int)ctx->ESTVal) - ctx->NLOFFS;

                    ctx->ESTyp[ctx->ESCnt] = (int)ctx->ESTVal;
                    ctx->ESIdx[ctx->ESCnt] = iabs(ctx, ctx->ESTyp[ctx->ESCnt]);
                    ctx->ESVal[ctx->ESCnt] = (double)(iabs(ctx, (int)ctx->ESTVal) - ctx->NLOFFS);
                    ctx->ESCnt++;
                    v_match(ctx, IDNL);

                    if (ctx->ESNxt == '(') {
                        if (v_getindex(ctx, ctx->ESTyp[ctx->ESCnt - 1])) {
                            if (ctx->MatExprCnt == 0)     
                                ctx->MatExprRow = imax(ctx, ctx->MatExprRow,ctx->NOC);
                        }
                    }
                    else if (ctx->MatExprCnt == 0) {
                        ctx->MatExprRow = imax(ctx, ctx->MatExprRow,ctx->NOC);
                        ctx->MatExprCol = imax(ctx, ctx->MatExprCol,ctx->NLNV[nn]);
                    }
                    return;

        case IDC:   ctx->ESIdx[ctx->ESCnt] = 0;
                    ctx->ESTyp[ctx->ESCnt] = (int)ctx->ESTVal;
                    ctx->ESVal[ctx->ESCnt++] = 0.0;
                    v_match(ctx, IDC);
                    if (ctx->ESNxt == '(')  
                        ctx->ESErr = -25;
                    break;

        case IDG:   tv = (int)ctx->ESTVal;
   
                    if (tv == 119)
                        ctx->IVEXPRFlg = 1;  

                    tva = iabs(ctx, tv);
                    na = ctx->ESNArg;
                    v_match(ctx, IDG);
                    if (ctx->ESErr)  
                        return;

                    if (tva == 415 || tva == 416) {     /* row and col */
                        ctx->ESIdx[ctx->ESCnt] = 0;
                        ctx->ESTyp[ctx->ESCnt] = 417;             
                        ctx->ESVal[ctx->ESCnt++] = 0.0;
                        if (ctx->ESCnt >= ctx->ESMaxLen) {
                            ctx->ESErr = -1;
                            return; 
                        }
                        ctx->MatExprCnt = 1;   
                    }

                    /* special for int (integration) */

                    if (tva == 186 || tva == 189 || tva == 193) {
                        if (ctx->ESCnt != 0) {
                            ctx->ESErr = -16;
                            return;
                        }
                        if (tva == 186)  
                            tvaa = 184;
                        else if (tva == 189)  
                            tvaa = 187;
                        else                     
                            tvaa = 190;

                        ctx->ESIdx[ctx->ESCnt] = 0;
                        ctx->ESTyp[ctx->ESCnt] = 184;
                        ctx->ESVal[ctx->ESCnt] = (double)tvaa;
                        ctx->ESCnt++;
                    }   

                    /* special for pre and suc operators */

                    if (tva == 1510 || tva == 1509) {
              
                        if (ctx->ESNxt != '(') {
                            ctx->ESErr = -3;
                            return;
                        }
                        if ((i = v_search1(ctx, ctx->ESP,&len)) >= 0) {
                            ctx->ESP += len;
                            if (*ctx->ESP != ')') {
                                ctx->ESErr = -11;
                                return;
                            }
                            ctx->ESIdx[ctx->ESCnt] = 0;
                            ctx->ESTyp[ctx->ESCnt] = tv;
                            ctx->ESVal[ctx->ESCnt++] = (double)i;
                            v_match(ctx, ctx->ESNxt);
                            v_match(ctx, ')');
                            return;
                        }
                        else if (tva == 1509) {
                            ctx->ESErr = -11;
                        }
                        else if (ctx->PREVName != NULL && ctx->PREVNum >= 0) {

                            /* check for self-referencing pre() */

                            len = (int)(strlen(ctx->PREVName));
                            if (!strncmp(ctx->ESP,ctx->PREVName,(size_t)(len))) {
                                ctx->ESP += len;
                                if (*ctx->ESP != ')') {
                                    ctx->ESErr = -11;
                                    return;
                                }
                                ctx->ESIdx[ctx->ESCnt] = 0;
                                ctx->ESTyp[ctx->ESCnt] = tv;
                                ctx->ESVal[ctx->ESCnt++] = (double)ctx->PREVNum;
                                v_match(ctx, ctx->ESNxt);
                                v_match(ctx, ')');
                                return;
                            }
                            else  
                                ctx->ESErr = -11;
                        }
                        else  
                            ctx->ESErr = -11;
                        return;
                    }
                    v_match(ctx, '(');

                    if (tva == 418) {               /* exists(X) */
                        ctx->ESIdx[ctx->ESCnt] = 0;
                        ctx->ESTyp[ctx->ESCnt] = 418;
                        if (ctx->ESErr || ctx->ESNxt == ')')
                            ctx->ESVal[ctx->ESCnt++] = -1.0;
                        else  
                            ctx->ESVal[ctx->ESCnt++] = 0.0;
                        if (ctx->ESNxt != ')') 
                            v_match(ctx, ctx->ESNxt);
                        ctx->ESErr = 0;
                        v_match(ctx, ')');
                        return;
                    }
                    if (ctx->ESErr)
                        return;

                    if (tva == 123) {           /* mvn(A,k,eps,x1,...,xn) */

                        if (ctx->ESNxt != IDM)
                            ctx->ESErr = -14;
                        else {
                            nn = iabs(ctx, (int)ctx->ESTVal) - ctx->MOFFS;
                            if (ctx->MatCol[nn] < 2)  
                                ctx->ESErr = -15;
                            else {
                                ctx->ESIdx[ctx->ESCnt] = 0;
                                ctx->ESTyp[ctx->ESCnt] = (int)ctx->ESTVal;
                                ctx->ESVal[ctx->ESCnt++] = (double)nn;
                                v_match(ctx, IDM);
                            }
                        }
                    }
                    else
                        v_expr(ctx);

                    if (ctx->ESErr)
                        return;

                    if (na == -1) {     /* special for sequences */
                        v_match(ctx, ',');
                        if (ctx->ESErr)
                            return;
                        v_match(ctx, ',');
                        if (ctx->ESErr)
                            return;
                        v_factor(ctx);         /* must be a variable */
                        if (ctx->ESErr)
                            return;

                        v1 = (int)ctx->ESVal[ctx->ESCnt - 2];
                        v2 = (int)ctx->ESVal[ctx->ESCnt - 1];
                        if (v1 != ctx->ESTyp[ctx->ESCnt - 2] - ctx->VOFFS ||
                            v2 != ctx->ESTyp[ctx->ESCnt - 1] - ctx->VOFFS || v1 > v2) {
                            ctx->ESErr = -18;
                            return;
                        }
                        nn = 2;
                        ctx->ESCnt--;
                        i = ctx->VNxt[v1];
                        while (i < v2) {
                            if (ctx->ESCnt + 2 >= ctx->ESMaxLen) {
                                ctx->ESErr = -1;
                                return;
                            }
                            ctx->ESIdx[ctx->ESCnt] = 0;
                            ctx->ESVal[ctx->ESCnt] = (double)i;
                            ctx->ESTyp[ctx->ESCnt++] = ctx->VOFFS + i;
                            nn++;
                            i = ctx->VNxt[i];
                        }
                        ctx->ESIdx[ctx->ESCnt] = 0;
                        ctx->ESVal[ctx->ESCnt] = (double)v2;
                        ctx->ESTyp[ctx->ESCnt++] = ctx->VOFFS + v2;
                        na = nn;
                    }
                    else if (na == 0) {  /* variable number of arguments */
                        nn = 1;
                        while (ctx->ESNxt == ',') {

                            if (ctx->ESCnt >= ctx->ESMaxLen) {
                                ctx->ESErr = -1;
                                return; 
                            }
                            nn++;
                            v_match(ctx, ',');
                            if (ctx->ESErr)
                                return;

                            v_expr(ctx);
                            if (ctx->ESErr)
                                return;
                        }
                        na = nn;
                    }
                    else {
                        if (tva == 209)
                            ctx->IFLevel++;

                        for (i = 1; i < na; ++i) {

                            if (ctx->ESCnt >= ctx->ESMaxLen) {
                                ctx->ESErr = -1;
                                return; 
                            }
                            v_match(ctx, ',');
                            if (ctx->ESErr)
                                return;

                            if (tva == 209) {       /* special for if */
                                ctx->ESIdx[ctx->ESCnt] = 0;
                                ctx->ESTyp[ctx->ESCnt] = ctx->IFTYP + ctx->IFLevel;
                                ctx->ESVal[ctx->ESCnt++] = (double)i;
                            }
                            else if (tva == 186) {  /* special for integration */
                                ctx->ESIdx[ctx->ESCnt] = 0;
                                ctx->ESTyp[ctx->ESCnt] = 184 + i;
                                ctx->ESVal[ctx->ESCnt] = (double)ctx->ESCnt;
                                ctx->ESCnt++;
                            }
                            else if (tva == 189) {  /* special for integration */
                                ctx->ESIdx[ctx->ESCnt] = 0;
                                ctx->ESTyp[ctx->ESCnt] = 187 + i;
                                ctx->ESVal[ctx->ESCnt] = (double)ctx->ESCnt;
                                ctx->ESCnt++;
                            }
                            else if (tva == 193) {  /* special for integration */
                                ctx->ESIdx[ctx->ESCnt] = 0;
                                ctx->ESTyp[ctx->ESCnt] = 190 + i;
                                ctx->ESVal[ctx->ESCnt] = (double)ctx->ESCnt;
                                ctx->ESCnt++;
                            }
                            if (tva == 413 || (i == 1 && tva == 414) || /* recode */
                                tva == 405) {                           /* rdmn */

                                if (ctx->ESNxt != IDM)
                                    ctx->ESErr = -14;
                                else {
                                    nn = iabs(ctx, (int)ctx->ESTVal) - ctx->MOFFS;
                                    if (ctx->MatCol[nn] < 2)
                                        ctx->ESErr = -15;
                                    else {
                                        ctx->ESIdx[ctx->ESCnt] = 0;
                                        ctx->ESTyp[ctx->ESCnt] = (int)ctx->ESTVal;
                                        ctx->ESVal[ctx->ESCnt++] = (double)nn;
                                        v_match(ctx, IDM);
                                    }
                                }
                            }
                            else
                                v_expr(ctx);
                            if (ctx->ESErr)
                                return;
                        }
                        if (tva == 209) {        /* special for if */
                            if (ctx->ESCnt >= ctx->ESMaxLen) {
                                ctx->ESErr = -1;
                                return; 
                            }
                            ctx->ESIdx[ctx->ESCnt] = 0;
                            ctx->ESTyp[ctx->ESCnt] = ctx->IFTYP + ctx->IFLevel;
                            ctx->ESVal[ctx->ESCnt++] = (double)i;
                            ctx->IFLevel--;
                        }
                    }
                    v_match(ctx, ')');
                    if (ctx->ESCnt >= ctx->ESMaxLen) {
                        ctx->ESErr = -1;
                        return; 
                    }
                    ctx->ESIdx[ctx->ESCnt] = 0;
                    ctx->ESTyp[ctx->ESCnt] = tv;
                    ctx->ESVal[ctx->ESCnt] = (double)na;

                    if ((tv >= 451 && tv <= 453) || tv == 1526) {
                        for (i = ctx->ESCnt - 1; i >= 0; --i) {
                            nn = ctx->ESTyp[i];    
                            if (nn >= ctx->VOFFS && nn < ctx->COFFS) {
                                ctx->ESVal[ctx->ESCnt] = (double)nn;
                                break;
                            }
                        }
                    }
                    ctx->ESCnt++;

                    if (tva == 415 || tva == 416) 
                        ctx->MatExprCnt = 0;    
                    if (ctx->ESNxt == '(')  
                        ctx->ESErr = -25;
                    break;

        case '<':                           /* vector operators */
                    if (ctx->MatExprFlg == 0) {
                        ctx->ESErr = -22;
                        return;
                    }
                    v_match(ctx, '<');
                    if (ctx->ESErr)
                        return;
                    v_expr(ctx);
                    if (ctx->ESErr)
                        return;
                    nn = 1;
                    while (1) {
                        if (ctx->ESNxt == ',') {
                            nn++;
                            v_match(ctx, ',');
                            if (ctx->ESErr)   
                                return;
                            v_expr(ctx);
                        }
                        else
                            break;
                        if (ctx->ESErr)
                            return;
                    }
                    v_match(ctx, '>');
                    ctx->ESIdx[ctx->ESCnt] = 0;
                    ctx->ESVal[ctx->ESCnt] = (double)nn;
                    if (ctx->ESNxt == '\'') {
                        v_match(ctx, '\'');
                        ctx->ESTyp[ctx->ESCnt] = RVECTYP;
                        if (ctx->MatExprCnt == 0) {
                            ctx->MatExprRow = imax(ctx, ctx->MatExprRow,1);
                            ctx->MatExprCol = imax(ctx, ctx->MatExprCol,nn);
                        }
                    }
                    else {
                        ctx->ESTyp[ctx->ESCnt] = CVECTYP;
                        if (ctx->MatExprCnt == 0) {
                            ctx->MatExprRow = imax(ctx, ctx->MatExprRow,nn);
                            ctx->MatExprCol = imax(ctx, ctx->MatExprCol,1);
                        }
                    }
                    ctx->ESCnt++;
                    break;

        default:    ctx->ESErr = -2;
                    return;
    }
}

/* ------------------------------------------------------------------------ */
/*  match                                                                   */

void v_match(TDAContext *ctx, int t)
{
    if (ctx->ESNxt == t)
        ctx->ESNxt = v_lexan(ctx);
    else {
        ctx->ESErr = -3;
        return; 
    }
}

/* ------------------------------------------------------------------------ */
/*  lexan: scan the input stream.                                           */
/*                                                                          */
/*  Return: IDV   for V, C, or X variables                                  */
/*          IDM   for matrices                                              */
/*          IDG   for operators                                             */
/*          NUM   floating point numbers                                    */
/*          IDC   constants and operators without arguments                 */
/*                                                                          */
/*  Sets the following error codes in ESErr:                                */
/*  ESErr =  0  no error                                                    */
/*          -3  syntax error                                                */
/*          -4  undefined operator                                          */
/*          -6  scan error for floating point numbers                       */
/*          -7  misplaced minus sign                                        */
/*          -8  undefined V variable                                        */
/*          -9  Xi variable with i <= 0                                     */
/*         -10  undefined name of matrix/variable/namelist.                 */  
/*         -12  if error in ci reference                                    */  

int v_lexan(TDAContext *ctx)              
{
    register int i;
    register char c,*p;
    int n,l,len,narg,lflag;
    double tmp;
         
    if (*ctx->ESP == '-') {
        if (ctx->ESTLLast == 1) {  /* if "-" should be part of following expression */
            ctx->ESTLMin = 1;
            ctx->ESTLLast = 0;
            ctx->ESP++;
        }
        else {
            ctx->ESTLLast = 0;             
            return((int)*ctx->ESP++);
        }
    }
    ctx->ESTLLast = 0;

    /* ### */
    
    if (*ctx->ESP == '$') {
        printf1(ctx, "currmac=%d\n",ctx->CurMacroNum);
        /**
        if (MacroExLevel < 1 || MacroExLevel >= MaxMat || sscanf(ESP,"$%d",&n) != 1) {
            ESErr = -24;
            return(0);
        }
        if (n < 1 || n > MacroNP[.....])  
            -1
        else if (MacroArg[...][n] == NULL)
            0
        else
            1
        **/

    }
     
    if (check_vname(ctx, ctx->ESP)) {     /* check for variable/matrix name */

        l = get_vnlen(ctx, ctx->ESP);
        c = *(ctx->ESP + l);
        *(ctx->ESP + l) = '\0';

        lflag = check_local(ctx, ctx->ESP);         

        i = mat_getidx(ctx, ctx->ESP,0);
        *(ctx->ESP + l) = c;
        if (i < 0 && lflag) {       /* local matrix name not found */
            ctx->ESErr = -23;
            return(IDM);
        }
        if (i >= 0) {               /* found matrix name */
            ctx->ESP += l;
            ctx->ESTVal = (double)(ctx->MOFFS + i);
            if (ctx->ESTLMin) { 
                ctx->ESTVal = -ctx->ESTVal;
                ctx->ESTLMin = 0;
            }
            return(IDM);
        }
        if ((i = v_search1(ctx, ctx->ESP,&len)) >= 0) {   /* search for var name */
            ctx->ESP += len;
            ctx->ESTVal = (double)(ctx->VOFFS + i);
            if (ctx->ESTLMin) { 
                ctx->ESTVal = -ctx->ESTVal;
                ctx->ESTLMin = 0;
            }
            ctx->EVFlag = 1;
            return(IDV);
        }
        *(ctx->ESP + l) = '\0';
        i = nl_check(ctx, ctx->ESP);          /* check for namelist */
        *(ctx->ESP + l) = c;
        if (i >= 0) {
            ctx->ESP += l;
            ctx->ESTVal = (double)(ctx->NLOFFS + i);
            if (ctx->ESTLMin) { 
                ctx->ESTVal = -ctx->ESTVal;
                ctx->ESTLMin = 0;
            }
            return(IDNL);
        }
        ctx->ESP += l;
        ctx->ESErr = -10;
    }
    else if (*ctx->ESP >= 'a' && *ctx->ESP <= 'z') {

        /* search for C or X variable */

        ctx->ESTVal = 0.0;
        if (*ctx->ESP == 'c' && sscanf(ctx->ESP,"c%d",&n) == 1) {
            if (n < 1 || n > ctx->CIdxMax)
                ctx->ESErr = -12;
            else  
                ctx->ESTVal = (double)(ctx->COFFS + n);
        }
        else
            n = 0;

        if (ctx->ESErr)
            return(0);

        if (n > 0) {
            if (ctx->ESTLMin) { 
                ctx->ESTVal = -ctx->ESTVal;
                ctx->ESTLMin = 0;
            }
            ctx->ESP++;
            while (*ctx->ESP && isdigit((int)*ctx->ESP)) ctx->ESP++;
            return(IDV);
        }
        else {      /* search for an operator */

            i = v_search(ctx, ctx->ESP,&len,&narg);

            if (i < 0) {
                ctx->ESErr = -4;
            }
            else {
                if (len >= 0)
                    ctx->ESP += len;
                else {              /* special type 2 operators */
                    ctx->ESP -= len;
                    if (*ctx->ESP != '[') {
                        ctx->ESErr = -3;
                        return((int)*ctx->ESP);
                    }
                    ctx->ESP = skip_blev(ctx, ctx->ESP);
                }
                ctx->ESTVal = (double)i;
                ctx->ESNArg  = narg;
            }
            if (ctx->ESTLMin) { 
                ctx->ESTVal = -ctx->ESTVal;
                ctx->ESTLMin = 0;
            }
            if (i >= 1890)
                return(IDC);
            else
                return(IDG);
        }
    }

    /* check for numerical entry (integer or floating point) */

    else if (isdigit((int)*ctx->ESP) || *ctx->ESP == '.') {       /* changed 8/5/97 */
        if (ctx->ESTLMin) {
            ctx->ESTLMin = 0;
            ctx->ESP--;
        }
        if (sscanf(ctx->ESP,"%lf",&tmp) == 1) {
            ctx->ESTVal = tmp; 
            p = skip_const(ctx, ctx->ESP);  
            if (p == NULL)
                ctx->ESErr = -6;
            else
                ctx->ESP = p;
        }
        else
            ctx->ESErr = -6;
        return(NUM);
    }
    else if (*ctx->ESP == '#') {                   /* hexadecimal number */
        unsigned int hx;

        if (sscanf(++ctx->ESP,"%x",&hx) == 1) {
            /* the literal keeps its old value: glibc wrote the bit
               pattern into the int, so #ffffffff was -1, and it stays
               -1 -- the unsigned scan only removes the UB above
               0x7fffffff, it does not change what TDA computed */
            ctx->ESTVal = (double)(int)hx;
            if (ctx->ESTLMin) {
                ctx->ESTLMin = 0;
                ctx->ESTVal = -ctx->ESTVal;
            }
            while (*ctx->ESP && (isdigit((int)*ctx->ESP) || (*ctx->ESP >= 'a' && *ctx->ESP <= 'f')))
                ctx->ESP++;
        }
        else
            ctx->ESErr = -6;
        return(NUM);
    }
    else if (ctx->ESTLMin) {

        if (*ctx->ESP == '.') {
            ctx->ESTLMin = 0;
            ctx->ESP--;
            if (sscanf(ctx->ESP,"%lf",&tmp) == 1) {
                ctx->ESTVal = tmp; 
                p = skip_const(ctx, ctx->ESP);  
                if (p == NULL)
                    ctx->ESErr = -6;
                else
                    ctx->ESP = p;
            }
            else
                ctx->ESErr = -6;
            return(NUM);
        }
        else if (*ctx->ESP == '(') {
            ctx->ESTLMin = 0;
            ctx->EST0Last = 1;
            ctx->ESTVal = 0.0;
            return(NUM);
        }
        else {
            ctx->ESErr = -7;
            return((int)*ctx->ESP++);
        }
    }
    else if (ctx->EST0Last) {
        ctx->EST0Last = 0;
        return((int) '-');
    }
    else {
        switch (*ctx->ESP) {
            case '+':
            case '*':
            case '&':
            case '/':
            case '^':
            case '%':   
            case '(':   
            case '[': 
            case '<': 
            case ',':   ctx->ESTLLast = 1;
                        break;
            default:    break;
        }
        return((int)*ctx->ESP++);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  p = skip_const(p)                                                       */
/*      It is assumed that p points to a string with integer or floating    */
/*      point number. p is then set to the next character after the end     */
/*      of the numerical expression. If an error occurs, the function       */
/*      returns NULL.                                                       */

char *skip_const(TDAContext *ctx, char *p)                    
{
    (void)ctx;        /* unused: the signature is shared */
    register int pcnt = 0, pcnt1 = 0, pcnt2 = 0, pcnt3 = 0;

    if (*p == '+' || *p == '-')  
        p++;

    if (isdigit((int)*p))
        pcnt = 1;
    else
        p--;

    while (*++p) {
        if (isdigit((int)*p))
            pcnt++;
        else if (*p == '.') {
            if (pcnt1++)  
                break;
        }
        else if (*p == 'E' || *p == 'e') {
            pcnt = 0;
            if (pcnt2++)
                return(NULL);
        }
        else if (*p == '+' || *p == '-') {
            if (!pcnt2 || pcnt3++)
                break;
        }
        else
            break;
    }
    if (!pcnt)  
        return(NULL);
    else  
        return(p);
}

/* ------------------------------------------------------------------------ */
/*  v_search  Check for defined operators. If FNPN > 0 also check the       */
/*            entries in FNPDef,FNPLen. Then return IAOFFS + i. And entries */
/*            in FNArgDef,FNArgLen. Then return FAOFFS + i.                 */

int v_search(TDAContext *ctx, char *p,int *len,int *narg)
{
    register int l,ll,i;
    register char *q;

    ll = i = 0;
    q = p;
    while (*q) {
        if (islower((int)*q) || isdigit((int)*q))
            ll++;
        else
            break;
        q++;
    }
    if (ll == 0)
        return(-1);

    while (1) {
        if (OPS[i].typ == 0)
            break;
        l = OPS[i].len;
        if (ll == l) {
            if (!strncmp(p,OPS[i].op,(size_t)(l))) {
                if (!OPS[i].narg || *(p + l) == '(') {
                    *len = OPS[i].len;
                    *narg = OPS[i].narg;
                    return((int)OPS[i].typ);
                }
            }
        }
        else if (ll == -l) {
            if (!strncmp(p,OPS[i].op,(size_t)(ll))) {
                if (!OPS[i].narg || *(p + ll) == '[') {
                    *len = OPS[i].len;    
                    *narg = OPS[i].narg;
                    return((int)OPS[i].typ);
                }
            }
        }
        i++;
    }
    if (ctx->FNPN > 0) {
        for (i = 0; i < ctx->FNPN; ++i) {
            l = ctx->FNPLen[i];
            if (l == ll && !strncmp(p,ctx->FNPDef[i],(size_t)(l))) {
                *len = l;
                *narg = 0;
                return(ctx->IAOFFS + i);
            }
        }
    }
    if (ctx->FNArgN > 0) {
        for (i = 0; i < ctx->FNArgN; ++i) {
            l = ctx->FNArgLen[i];
            if (l == ll && !strncmp(p,ctx->FNArgDef[i],(size_t)(l))) {
                *len = l;
                *narg = 0;
                return(ctx->FAOFFS + i);
            }
        }
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  v_search1  check for data matrix variables. Return internal variable    */
/*             number if found, otherwise -1.                               */

int v_search1(TDAContext *ctx, char *p,int *len)
{
    register int i,l;        

    l = get_vnlen(ctx, p);       /* get length of var name at p */

    if (l > 0) {
        i = ctx->VIFirst;
        while (i >= 0) {
            if (!strncmp(ctx->VName[i],p,(size_t)(l))) {
                if ((size_t)l == strlen(ctx->VName[i])) {
                    *len = l;
                    return(i);
                }
            }
            i = ctx->VNxt[i];
        }
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  prn_emsg1(n)    print error message for v_parse errors.                 */
/*                                                                          */
/*          -1  stack overflow                                              */
/*          -2  if part of input stream cannot be interpreted               */
/*          -3  match error                                                 */
/*          -4  undefined operator                                          */
/*          -6  scan error for floating point numbers                       */
/*          -7  misplaced minus sign                                        */
/*          -8  undefined V variable                                        */
/*          -9  Xi variable with i <= 0                                     */
/*         -10  undefined name of matrix/var/namelist                       */  
/*         -11  if in pre(X) or suc(X), X is not a variable name            */
/*         -12  error in ci reference                                       */
/*         -13  cannot use matrix reference                                 */
/*         -14  need reference to a matrix                                  */
/*         -15  matrix should have at least two columns                     */
/*         -16  integral must be first term in expression                   */
/*         -17  t operator may only occur inside integral                   */
/*         -18  error in sequence definition                                */
/*         -19  inconsistent matrix dimensions                              */
/*         -20  cannot use namelist                                         */
/*         -21  error in matrix indices                                     */
/*         -22  cannot use vector expression                                */
/*         -23  reference to undefined local matrix                         */
/*         -24  unvalid reference to macro argument                         */
/*         -25  misplaced bracket                                           */
/*         -26  expressoin must not contain intervals                       */
/*         -27  cannot combine str with other operators                     */
    
void prn_emsg1(TDAContext *ctx, int n)
{
    printf1(ctx, "Error type: ");
    switch (n) {
        case -1:    printf1(ctx, "parser stack overflow.\n");
                    break;
        case -2:    printf1(ctx, "expression contains unknown operators.\n");
                    break;
        case -3:    printf1(ctx, "syntax error.\n");
                    break;
        case -4:    printf1(ctx, "expression contains undefined operator or argument.\n");
                    break;
        case -6:    printf1(ctx, "scan error for floating point number.\n");
                    break;
        case -7:    printf1(ctx, "misplaced minus sign.\n");
                    break;
        case -8:    printf1(ctx, "reference to an undefined variable.\n");
                    break;
        case -9:    printf1(ctx, "wrong reference to an Xi term.\n");
                    break;
        case -10:   printf1(ctx, "probably wrong reference to a variable or matrix.\n");
                    break;
        case -11:   printf1(ctx, "error in pre or suc operator, syntax is: pre(VName), suc(VName).\n");
                    printf1(ctx, "suc operator must refer to a previously defined variable.\n");
                    break;
        case -12:   printf1(ctx, "ci term needs positive integer i.\n");
                    break;
        case -13:   printf1(ctx, "cannot use matrix reference.\n");
                    break;
        case -14:   printf1(ctx, "need reference to a matrix.\n");
                    break;
        case -15:   printf1(ctx, "matrix should have at least two columns.\n");
                    break;
        case -16:   printf1(ctx, "integral must be first term in expression.\n");
                    break;
        case -17:   printf1(ctx, "t operator may only occur inside integral.\n");
                    break;
        case -18:   printf1(ctx, "error in sequence definition.\n");
                    break;
        case -19:   printf1(ctx, "inconsistent matrix dimensions.\n");
                    break;
        case -20:   printf1(ctx, "cannot use reference to namelist.\n");
                    break;
        case -21:   printf1(ctx, "syntax error in matrix indices.\n");
                    break;
        case -22:   printf1(ctx, "cannot use vector expression.\n");
                    break;
        case -23:   printf1(ctx, "reference to undefined local matrix.\n");
                    break;
        case -24:   printf1(ctx, "unvalid reference to macro argument.\n");
                    break;
        case -25:   printf1(ctx, "misplaced bracket.\n");
                    break;
        case -26:   printf1(ctx, "expression may not contain intervals.\n");
                    break;
        case -27:   printf1(ctx, "cannot combine str with other operators.\n");
                    break;
        default:    printf1(ctx, "unspecified.\n");
                    break;
    }
}

/* ------------------------------------------------------------------------ */
/*  prn_emsg2(n)    print error message for v_eval errors.                  */
/*                                                                          */
/*  -1  =   undefined X or C or V variable                                  */
/*   1  =   undefined operator                                              */
/*   2  =   division by zero                                                */
/*   3  =   % operator finds non-integers                                   */
/*   4  =   exp out of range                                                */
/*   5  =   argument of log is too small (< LogMin)                         */
/*   6  =   negative argument in sqrt                                       */
/*   7  =   negative or zero argument for lgam                              */
/*   8  =   arguments for icg out of range                                  */
/*   9  =   arguments for icb out of range                                  */
/*  10  =   incorrect range in tr(x,a,b)  (a > b)                           */
/*  11  =   argument for ndi() out ouf range [0,1]                          */
/*  12  =   incorrect arguments for TD                                      */
/*  13  =   incorrect arguments for CHD                                     */
/*  14  =   incorrect arguments for FD                                      */
/*  15  =   invalid quantile request                                        */
/*  16  =   rdmn(i,A), i out of range                                       */
/*  17  =   undefined intermediate function parameter                       */
/*  18  =   undefined function argument                                     */
/*  24  =   error: operator cannot be used for derivatives                  */
/*  27  =   error: integration interval not positive                        */
/*  28  =   error: no success in numerical integration                      */
/*  29  =   error: cannot use operator for interval expression              */
/*  30  =   error: in number of Hermite integration points                  */
/*  31  =   error: misplaced integration variable (t)                       */
/*  32  =   error: cannot calculate incomplete gamma integral               */
/*  33  =   error: cannot calculate binomial coefficient                    */
/*  34  =   error: cannot calculate biv normal distribution                 */
/*  35  =   error: cannot init rdmn generator                               */
/*  36  =   error: wrong index in rdmn operator                             */
/*  37  =   error: need a square matrix                                     */
/*  38  =   error: exceeded max matrix dimensions                           */
/*  39  =   error: in number of arguments                                   */
/*  40  =   error: in mnv evaluation                                        */
/*  41  =   error: in correlation matrix                                    */
/*  42  =   error: in poisson operator                                      */
/*  43  =   error: in negbin operator                                       */
/*  44  =   error: can only calculate first derivative                      */
/*  45  =   error: in matrix indices                                        */
/*  46  =   error: argument of rdp1 out of range.                           */
/*  47  =   error: need reference to a string variable                      */
/*  48  =   error: in range of arguments                                    */
/*  49  =   error: operator requires spatial data                           */
/*  50  =   error: cannot read spatial data                                 */
/*                                                                          */
/*  90  =   syntax error (missing operator?)                                */
/*  92  =   insufficient stack space for derivatives                        */
/*  95  =   insufficient memory                                             */
/*  96  =   no cases in data matrix                                         */
/*  97  =   out of stack space                                              */
/*  98  =   if expression contains references to c1,c2,...                  */
/*  99  =   operator that cannot be used in connection with type 2 ops.     */

void prn_emsg2(TDAContext *ctx, int n)
{
    printf1(ctx, "Error type: ");
    switch (n) {
        case -1:    printf1(ctx, "undefined reference (variable, x, y, or c term).\n");
                    break;
        case  1:    printf1(ctx, "expression contains undefined operator or argument.\n");
                    break;
        case  2:    printf1(ctx, "division by zero.\n");
                    break;
        case  3:    printf1(ctx, "non-integers in modulus operator.\n");
                    break;
        case  4:    printf1(ctx, "argument of exp exceeds limits.\n");
                    break;
        case  5:    printf1(ctx, "argument of log is less than or equal to zero.\n");
                    break;
        case  6:    printf1(ctx, "negative argument in sqrt.\n");
                    break;
        case  7:    printf1(ctx, "negative or zero argument in lgam, digam or trigam.\n");
                    break;
        case  8:    printf1(ctx, "icg arguments out of range.\n");
                    break;
        case  9:    printf1(ctx, "icb arguments out of range.\n");
                    break;
        case 10:    printf1(ctx, "incorrect arguments for tr(x,a,b).\n");
                    break;
        case 11:    printf1(ctx, "argument for ndi out of range.\n");
                    break;
        case 12:    printf1(ctx, "argument for td out of range.\n");
                    break;
        case 13:    printf1(ctx, "argument for chd out of range.\n");
                    break;
        case 14:    printf1(ctx, "argument for fd out of range.\n");
                    break;
        case 15:    printf1(ctx, "invalid quantile request.\n");
                    break;
        case 16:    printf1(ctx, "error in argument to rdmn.\n");
                    break;
        case 17:    printf1(ctx, "error: undefined intermediate function parameter.\n");                   
                    break;
        case 18:    printf1(ctx, "error: undefined function argument.\n");                                 
                    break;
        case 24:    printf1(ctx, "expression contains an operator that cannot\n"); 
                    printf1(ctx, "be used for derivatives.\n");
                    break;
        case 27:    printf1(ctx, "integration interval not positive.\n");
                    break;
        case 28:    printf1(ctx, "no success in numerical integration.\n");
                    break;
        case 29:    printf1(ctx, "at least one operator is not allowed in interval expression.\n");
                    break;
        case 30:    printf1(ctx, "error in number of Hermite integration points.\n");
                    break;
        case 31:    printf1(ctx, "misplaced integration variable t.\n");
                    break;
        case 32:    printf1(ctx, "cannot calculate incomplete gamma integral.\n");
                    break;
        case 33:    printf1(ctx, "cannot evaluate bc().\n");
                    break;
        case 34:    printf1(ctx, "cannot evaluate bivn().\n");
                    break;
        case 35:    printf1(ctx, "cannot init rdmn generator.\n");
                    break;
        case 36:    printf1(ctx, "wrong index in rdmn operator.\n");
                    break;
        case 37:    printf1(ctx, "need a square matrix.\n");
                    break;
        case 38:    printf1(ctx, "exceeded maximal matrix dimension.\n");
                    break;
        case 39:    printf1(ctx, "wrong number of arguments.\n");
                    break;
        case 40:    printf1(ctx, "cannot evaluate mvn (for required accuracy).\n");
                    break;
        case 41:    printf1(ctx, "improper correlation matrix.\n");
                    break;
        case 42:    printf1(ctx, "error in poisson operator.\n");
                    break;
        case 43:    printf1(ctx, "error in negbin operator.\n");
                    break;
        case 44:    printf1(ctx, "can only calculate first derivative.\n");
                    break;
        case 45:    printf1(ctx, "matrix indices not suitable for current object.\n");
                    break;
        case 46:    printf1(ctx, "argument for rdp1 out of range.\n");
                    break;
        case 47:    printf1(ctx, "need reference to a string variable.\n");
                    break;
        case 48:    printf1(ctx, "error in range of arguments.\n");
                    break;
        case 49:    printf1(ctx, "operator requires spatial data.\n");
                    break;
        case 50:    printf1(ctx, "cannot read spatial data.\n");
                    break;

        case 90:    printf1(ctx, "syntax error (missing operator?).\n");
                    break;
        case 92:    printf1(ctx, "insufficient stack space for derivatives.\n");
                    break;
        case 93:    printf1(ctx, "if-then-else incompatible with type 2 operators.\n");
                    break;
        case 94:    printf1(ctx, "can't evaluate this if-then-else construction\n");
                    printf1(ctx, "with type 2 operators and/or type 4 variables.\n");
                    break;
        case 95:    printf1(ctx, "insufficient memory.\n");
                    break;
        case 96:    printf1(ctx, "no cases in data matrix.\n");
                    break;
        case 97:    printf1(ctx, "insufficient stack space.\n");
                    break;
        case 98:    printf1(ctx, "type 2 operator contains c terms.\n");
                    break;
        case 99:    printf1(ctx, "wrong combination of operator types.\n");
                    break;
        default:    printf1(ctx, "unspecified (%d).\n",n);
                    break;
    }
}


