/****************************************************************************/
/*  t_eval3                                                                 */
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
#include "t_eval.h"
#include "t_eval1.h"
#include "tda_context.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_eval3.c                                                  */

int eval_specf(TDAContext *ctx, int atyp,int vn,int noc,int sl,double **buf);
void check_expr(TDAContext *ctx, short cnt,int *ptyp,int *nv,int *nc,int *n2,int *n3,int opt);
int evalf(TDAContext *ctx, int iflag);
int fnd_alloc(TDAContext *ctx, int opt,int deriv,int fnn,int pn,int iflag);
int deriv0(TDAContext *ctx, int opt,int deriv,int j,int iflag);
int deriv1(TDAContext *ctx, int deriv,int ix,int j,int iflag);
int deriv2(TDAContext *ctx, int op,int deriv,int j,double *val,int iflag);
int deriv3(TDAContext *ctx, int op,int deriv,int j,double *val,int iflag);
int derivf1(TDAContext *ctx, int op,int deriv,int j,double *val,int iflag);
int derivf2(TDAContext *ctx, int op,int deriv,int j,double *val,double *dval);
int derivf3(TDAContext *ctx, int op,int deriv,int j,double *val,int iflag);
double get_integrand(TDAContext *ctx, double x,int *err);
int ia_op(TDAContext *ctx, int typ,int j,double *buf);
void ia_mul(TDAContext *ctx, double a,double b,double c,double d,double *ra,double *rb);
void ia_iv(TDAContext *ctx, double a,double b,double c,double d,double *ra,double *rb);
void ia_sin(TDAContext *ctx, double a,double b,double *ra,double *rb);
void ia_cos(TDAContext *ctx, double a,double b,double *ra,double *rb);

/* ------------------------------------------------------------------------ */
/*  eval_specf(atyp,vn,noc,sl,buf)                                          */
/*  ##                      Evaluation of special type 2 operators          */
/*                          for variable iv. sl is stacklevel for buf[sl]   */
/*                          return 0 if OK, else error code.                */
/*                                                                          */
/*  case 1700:              sma[](X)                                        */
/*  case 1701:              smd[](X)                                        */
/*                                                                          */

int eval_specf(TDAContext *ctx, int atyp,int vn,int noc,int sl,double **buf)
{     
    register int j;
    register char c,*p,*q;
    int err = 0;

    if (vn < 0)
        return(26);

    p = ctx->VDef[vn];
    while (*p && *p != '[')
        p++;
    q = skip_blev(ctx, p);
    if (*p != '[' || *--q != ']')  
        return(26);
         
    *p = '(';
    *q = ')';
    c = *++q;
    *q = '\0';

    if (parm(ctx, p,0,0)) {  /* get parameters */
        err = 25;
        goto ESPECFin;
    }
    if (atyp == 1700) {
        if (ctx->PMR < 1)
            ctx->PMR = 1;
        if (ctx->PMNTP < 2 || ctx->PMOPT < 1 || ctx->PMOPT > 2) {
            err = 25;
            goto ESPECFin;
        }
        j = sl - 1;
        if (s_sma(ctx, noc,buf[j],ctx->PMNTP,ctx->PMTP,ctx->PMR,ctx->PMOPT)) {
            err = 25;
            goto ESPECFin;
        }
    }
    else if (atyp == 1701) {
        if (ctx->PMRHSTRA < 1 || ctx->PMOPT < 1 || ctx->PMOPT > 3) {
            err = 25;
            goto ESPECFin;
        }
        j = sl - 1;
        if (s_smd(ctx, ctx->PMRHSTR,ctx->PMOPT,noc,buf[j])) { 
            err = 25;
            goto ESPECFin;
        }
    }
ESPECFin:
    *q-- = c;
    *q = ']';
    *p = '[';
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  check_expr(cnt,ptyp,nv,nc,n2,n3,opt)                                    */
/*                                                                          */
/*  check expression, return                                                */
/*  nv = number of variable references                                      */
/*  nc = number of c term referendes                                        */
/*  n2 = number of type 2 operators or type 4 variables (if opt = 0)        */
/*  n3 = number of type 3 operators or type 5 variables (if opt = 0)        */
/*  If opt = 1 check only operators.                                        */ 
         
void check_expr(TDAContext *ctx, short cnt,int *ptyp,int *nv,int *nc,int *n2,int *n3,int opt)
{
    register int i,t;
                 
    *nc = 0;         
    *nv = 0;
    *n2 = 0;
    *n3 = 0;

    for (i = 0; i < cnt; ++i) {
        t = iabs(ctx, ptyp[i]);
        if (t >= ctx->VOFFS && t < ctx->COFFS) {
            *nv += 1;
            if (opt == 0) {
                if (ctx->VTyp[t - ctx->VOFFS] == 4)
                    *n2 += 1;
                else if (ctx->VTyp[t - ctx->VOFFS] == 5)
                    *n3 += 1;
            }
        }
        else if (t >= ctx->COFFS && t < ctx->COFFMAX)
            *nc += 1;
        else if (t >= 1500 && t < 1900)
            *n2 += 1;
        else if (t >= 1970 && t < ctx->IAOFFS)
            *n3 += 1;
    }
}

/* ------------------------------------------------------------------------ */
/*  evalf(iflag)    Command for function evaluation.                        */
/*  ##                                                                      */
/*          evalf(                                                          */
/*              fmt=...,            print format, def. 12.4                 */
/*              arg1=...,                                                   */
/*              arg2=...,                                                   */
/*          ) = function;                                                   */
/*                                                                          */
/*          arguments must be given as                                      */
/*                                                                          */
/*          arg = value,                for standard functions, or          */
/*          arg = lvalue,uvalue,        for inclusion functions.            */
/*                                                                          */
/*          evalf()   only function value                                   */
/*          evalf1()  with first derivaties                                 */
/*          evalf2()  with first and second derivatives                     */  
/*          evalfi()  inclusion function (if iflag = 1)                     */
/*          evalfi1() inclusion function first derivative (if iflag = 1)    */
/*                                                                          */
/*          Return 0 if OK, otherwise -1.                                   */

int evalf(TDAContext *ctx, int iflag)
{
    register int i,j,k;
    int err,aflg,r,n,m,deriv,ii,jj,n3,len;
    register char *p,*q,*s;
    double x,xu,tmp;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    len = 5;
    printf1(ctx, "Function evaluation. Current memory: %d bytes.\n",ctx->MemReq);
    if (iflag) {
        printf1(ctx, "Calculating an inclusion function.\n");
        len++;
    }
    deriv = 0;
    p = ctx->CmdBuf + len;
    if (*p == '1') {
        deriv = 1;
        p++;
    }
    else if (*p == '2') {
        deriv = 2;
        p++;
    }
    if (iflag && deriv > 1) {
        printf1(ctx, "Error: cannot calculate inclusion functions for second derivatives.\n");
        return(-1);
    }
    q = p;

    aflg = 0;
    if (*q == '(') {
        p = skip_blev(ctx, q);
        aflg = 1;
    }
    if (*p++ != '=' || !*p) {
        p_err(ctx, -1,1);
        goto EVALFFin;
    }
    r = get_func(ctx, p,1,&n3,iflag,NULL);
    if (r || ctx->FNFlg == 0)
        goto EVALFFin;

    if (n3 > 0) {
        p_err(ctx, -40,1);
        goto EVALFFin;
    }

    for (i = 0; i < ctx->FNArgN; ++i) {
        ctx->FNArgVal[i] = 0.0;              /* default arguments */
        if (iflag)
            ctx->FNArgVal1[i] = 0.0;
    }
    pmfmt(ctx, 12,4);                        /* default print format */

    if (aflg) {                         /* check arguments */
        while (*++q && *q != ')') {

            if (sscanf(q,"fmt=%d.%d",&n,&m) == 2) {
                pmfmt(ctx, n,m);
                q = skip_dbl(ctx, q + 4);
                continue;
            }
            k = -1;
            s = q;
            for (i = 0; i < ctx->FNArgN; ++i) {
                j = ctx->FNArgSP[i];
                if (!strncmp(q,ctx->FNArgDef[j],(size_t)(ctx->FNArgLen[j]))) {
                    q += ctx->FNArgLen[j];
                    if (sscanf(q,"=%lf",&x) != 1) {
                        err = -9;
                        goto EVALFFin;
                    }   
                    q = skip_dbl(ctx, q + 1);
                    xu = x;
                    if (iflag && sscanf(q,",%lf",&tmp) == 1) {
                        xu = tmp; 
                        q = skip_dbl(ctx, q + 1);
                    }
                    k = i;
                    break;
                }
            }
            if (k < 0) {
                err = -9;
                goto EVALFFin;
            }
            ctx->FNArgVal[k] = x;
            if (iflag)
                ctx->FNArgVal1[k] = xu;

            if (*q == ')')
                break;
        }
    }
                 
    if (ctx->FNArgN == 0) {
        printf1(ctx, "\nFunction without arguments.");
        deriv = 0;
    }
    else {
        printf1(ctx, "\nArguments:");
        r = 0;
        for (i = 0; i < ctx->FNArgN; ++i) {
            j = ctx->FNArgSP[i];
            if (iflag == 0)
                printf1(ctx, " %s=%lg",ctx->FNArgDef[j],ctx->FNArgVal[i]);
            else {
                printf1(ctx, " %s=[%lg,%lg]",ctx->FNArgDef[j],ctx->FNArgVal[i],ctx->FNArgVal1[i]);
                if (ctx->FNArgVal1[i] < ctx->FNArgVal[i])  
                    r++;
            }
        }
        newline(ctx);
        if (r) {
            printf1(ctx, "Error in argument intervals.\n");
            goto EVALFFin;
        }
        if (deriv > 0) {
            if (fnd_alloc(ctx, 1,deriv,ctx->FNArgN,ctx->FNPN,iflag))
                goto EVALFFin;
        }
    }
    newline(ctx);
    if (ctx->FVFlg)
        prn_feval(ctx);

    if (deriv) {        /* need memory for gradient and hessian */
          
        if (alloc_acu(ctx, ctx->FNArgN + 1))
            goto EVALFFin;
  
        if (alloc_acv(ctx, ctx->FNArgN * ctx->FNArgN + 1))
            goto EVALFFin;
    }
    ctx->NINTMUsed = -1;
    if (iflag == 0)
        r = get_flval(ctx, &x,ctx->FNArgN,ctx->FNArgVal,deriv,1,ctx->AcU,ctx->AcV,ctx->AcV);
    else
        r = get_ifval(ctx, &x,&xu,ctx->FNArgN,ctx->FNArgVal,ctx->FNArgVal1,deriv,ctx->AcU,ctx->AcV);

    if (r) {      /* r from v_eval1(ctx) */

        printf1(ctx, "Can't evaluate function or derivatives.\n");
        prn_emsg2(ctx, r);
        goto EVALFFin;
    }
    if (ctx->NINTMUsed >= 0)  
        ni_info(ctx);

    newline(ctx);
    if (iflag)
        printf1(ctx, "[");
    rt_printf1_d(ctx, ctx->PMFmtS,x);
    if (iflag) {
        printf1(ctx, ",");
        rt_printf1_d(ctx, ctx->PMFmtS,xu);
        printf1(ctx, "]");
    }
    printf1(ctx, "   function value\n");
#ifdef TDA_R_PACKAGE
    /* the interval the evaluation produced; iflag == 0 means a point
       value, which is exported as a degenerate [x,x] so the reader has
       one shape to handle */
    {
        double erow[2];
        erow[0] = x;
        erow[1] = iflag ? xu : x;
        tda_export_row(ctx, "evalfi.value", erow, 2);
    }
#endif

    if (deriv > 0) {
        for (i = 1; i <= ctx->FNArgN; ++i) {
            ii = ctx->FNArgSP[i - 1];

            if (iflag == 0) 
                rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcU[i]);
            else {
                printf1(ctx, "[");
                rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcU[i - 1]);
                printf1(ctx, ",");
                rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcV[i - 1]);
                printf1(ctx, "]");
            }
            printf1(ctx, "   gradient %s\n",ctx->FNArgDef[ii]);
#ifdef TDA_R_PACKAGE
            {
                double grow[2];
                grow[0] = iflag ? ctx->AcU[i - 1] : ctx->AcU[i];
                grow[1] = iflag ? ctx->AcV[i - 1] : ctx->AcU[i];
                tda_export_row(ctx, "evalfi.gradient", grow, 2);
            }
#endif
        }
        if (deriv > 1) {
            k = 0;
            for (i = 0; i < ctx->FNArgN; ++i) {
                ii = ctx->FNArgSP[i];
                for (j = 0; j <= i; ++j) {
                    jj = ctx->FNArgSP[j];
                    rt_printf1_d(ctx, ctx->PMFmtS,ctx->AcV[i * ctx->FNArgN + j + 1]);
                    printf1(ctx, "   hessian  %s %s\n",ctx->FNArgDef[ii],ctx->FNArgDef[jj]);
                }
            }
        }
    }
    err = 0;         

EVALFFin:
    if (err == -9) {
        p = skip_nc(ctx, s);
        *p = '\0';
        printf1(ctx, "Error in scanning function arguments: %s\n",s);
        err = -1;
    }
    fnd_alloc(ctx, 0,0,0,0,0);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  fnd_alloc(opt,deriv,fnn,pn,iflag)                                       */
/*                                                                          */
/*  If opt != 0 allocate memory for gradients and hessians, depending on    */  
/*  deriv, fnn,pn. If iflag != 0 allocate also FNPGrad1.                    */
/*                                                                          */
/*  If opt = 0 free previously allocated memory.                            */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int fnd_alloc(TDAContext *ctx, int opt,int deriv,int fnn,int pn,int iflag)
{
    int err,i;

    err = -1;
    if (opt == 0) {
        err = 0;
        goto FNDAFin;
    }
    if (deriv < 1)
        return(0);

    ctx->FNGradL = fnn;
    ctx->FNHessL = fnn * (fnn + 1) / 2;

    if (fnn <= 0)
        return(0);

    if (!(ctx->FNATyp = (short *)calloc((size_t)(ctx->MaxSTD),sizeof(short)))) { 
        p_err(ctx, -2,1);  
        goto FNDAFin;
    }
    memrq(ctx, ctx->MaxSTD,sizeof(short));
    ctx->FNATypA = ctx->MaxSTD;

    if (!(ctx->FNGrad = (double **)calloc((size_t)(ctx->MaxSTD),sizeof(double *)))) { 
        p_err(ctx, -2,1);  
        goto FNDAFin;
    }
    memrq(ctx, ctx->MaxSTD,sizeof(double *));
    ctx->FNGradAA = ctx->MaxSTD;

    if (deriv > 1) {
        if (!(ctx->FNHess = (double **)calloc((size_t)(ctx->MaxSTD),sizeof(double *)))) { 
            free((char *)ctx->FNGrad);
            p_err(ctx, -2,1);  
            goto FNDAFin;
        }
        memrq(ctx, ctx->MaxSTD,sizeof(double *));
        ctx->FNHessAA = ctx->MaxSTD;
    }
    for (i = 0; i < ctx->MaxSTD; ++i) {
        if (!(ctx->FNGrad[i] = (double *)calloc((size_t)(ctx->FNGradL),sizeof(double)))) { 
            p_err(ctx, -2,1);  
            goto FNDAFin;
        }
        memrq(ctx, ctx->FNGradL,sizeof(double));
        ctx->FNGradA++;
    }
    for (i = 0; i < pn; ++i) {
        if (!(ctx->FNPGrad[i] = (double *)calloc((size_t)(ctx->FNGradL),sizeof(double)))) { 
            p_err(ctx, -2,1);  
            goto FNDAFin;
        }
        memrq(ctx, ctx->FNGradL,sizeof(double));
        ctx->FNPGradA++;

        if (iflag) {
            if (!(ctx->FNPGrad1[i] = (double *)calloc((size_t)(ctx->FNGradL),sizeof(double)))) { 
                p_err(ctx, -2,1);  
                goto FNDAFin;
            }
            memrq(ctx, ctx->FNGradL,sizeof(double));
            ctx->FNPGrad1A++;
        }
    }
    if (deriv > 1) {
        for (i = 0; i < ctx->MaxSTD; ++i) {
            if (!(ctx->FNHess[i] = (double *)calloc((size_t)(ctx->FNHessL),sizeof(double)))) { 
                p_err(ctx, -2,1);  
                goto FNDAFin;
            }
            memrq(ctx, ctx->FNHessL,sizeof(double));
            ctx->FNHessA++;
        }
        for (i = 0; i < pn; ++i) {
            if (!(ctx->FNPHess[i] = (double *)calloc((size_t)(ctx->FNHessL),sizeof(double)))) { 
                p_err(ctx, -2,1);  
                goto FNDAFin;
            }
            memrq(ctx, ctx->FNHessL,sizeof(double));
            ctx->FNPHessA++;
        }
    }
    return(0);

FNDAFin:
    if (ctx->FNATypA > 0) {
        free((char *)ctx->FNATyp);
        memrq(ctx, -ctx->FNATypA,sizeof(short));
        ctx->FNATypA = 0;         
    }
    if (ctx->FNGradL > 0) {
        for (i = 0; i < ctx->FNGradA; ++i) {
            free((char *)ctx->FNGrad[i]);
            memrq(ctx, -ctx->FNGradL,sizeof(double));
        }
        if (ctx->FNPGradA > 0) {
            for (i = 0; i < ctx->FNPGradA; ++i) {
                free((char *)ctx->FNPGrad[i]);
                memrq(ctx, -ctx->FNGradL,sizeof(double));
            }
        }
        if (ctx->FNPGrad1A > 0) {
            for (i = 0; i < ctx->FNPGrad1A; ++i) {
                free((char *)ctx->FNPGrad1[i]);
                memrq(ctx, -ctx->FNGradL,sizeof(double));
            }
        }
        ctx->FNPGrad1A = ctx->FNPGradA = ctx->FNGradN = ctx->FNGradL = ctx->FNGradA = 0;
    }
    if (ctx->FNHessL > 0) {
        for (i = 0; i < ctx->FNHessA; ++i) {
            free((char *)ctx->FNHess[i]);
            memrq(ctx, -ctx->FNHessL,sizeof(double));
        }
        if (ctx->FNPHessA > 0) {
            for (i = 0; i < ctx->FNPHessA; ++i) {
                free((char *)ctx->FNPHess[i]);
                memrq(ctx, -ctx->FNHessL,sizeof(double));
            }
        }
        ctx->FNPHessA = ctx->FNHessN = ctx->FNHessL = ctx->FNHessA = 0;
    }
    if (ctx->FNGradAA > 0) {
        free((char *)ctx->FNGrad); 
        memrq(ctx, -ctx->MaxSTD,sizeof(double *));
        ctx->FNGradAA = 0;         
    }
    if (ctx->FNHessAA > 0) {
        free((char *)ctx->FNHess); 
        memrq(ctx, -ctx->MaxSTD,sizeof(double *));
        ctx->FNHessAA = 0;         
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  deriv0(opt,deriv,j,iflag)                                               */
/*                                                                          */
/*  opt = 0   set to 0                                                      */
/*  opt = 1   change sign                                                   */
/*                                                                          */
/*  Return 0 if OK,   92 if out of stack space.                             */

int deriv0(TDAContext *ctx, int opt,int deriv,int j,int iflag)
{
    register int i;
    double tmp;

    if (j != ctx->FNGradN || (deriv > 1 && ctx->FNGradN != ctx->FNHessN))
        gerr_exit(ctx, 84);
   
    if (opt == 0) {

        if (j + iflag >= ctx->MaxSTD)
            return(92);

        ctx->FNATyp[j] = 1;
        for (i = 0; i < ctx->FNGradL; ++i)
            ctx->FNGrad[j][i] = 0.0;
        ctx->FNGradN++;         
        if (iflag) {
            j++;
            ctx->FNATyp[j] = 1;
            for (i = 0; i < ctx->FNGradL; ++i)
                ctx->FNGrad[j][i] = 0.0;
            ctx->FNGradN++;         
        }    
        if (deriv > 1) {
            for (i = 0; i < ctx->FNHessL; ++i)
                ctx->FNHess[j][i] = 0.0;
            ctx->FNHessN++;         
        }
    }
    else if (opt == 1) {

        j--;
        if (j < 0)
            gerr_exit(ctx, 78);

        if (iflag == 0) {
            for (i = 0; i < ctx->FNGradL; ++i)
                ctx->FNGrad[j][i] *= -1;
        }
        else {
            for (i = 0; i < ctx->FNGradL; ++i) {
                tmp = ctx->FNGrad[j][i];
                ctx->FNGrad[j][i] = -ctx->FNGrad[j - 1][i];
                ctx->FNGrad[j - 1][i] = -tmp;
            }
        }
        if (deriv > 1) {
            for (i = 0; i < ctx->FNHessL; ++i)
                ctx->FNHess[j][i] *= -1;
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  deriv1(deriv,ix,j,iflag)                                                */

int deriv1(TDAContext *ctx, int deriv,int ix,int j,int iflag)
{
    register int i;

    if (j != ctx->FNGradN || (deriv > 1 && ctx->FNGradN != ctx->FNHessN)) {
        printfe(ctx, "j=%d fngradn=%d deriv=%d\n",j,ctx->FNGradN,deriv);
        gerr_exit(ctx, 85);
    }
    if (j + iflag >= ctx->MaxSTD)
        return(92);
  
    ctx->FNATyp[j] = 0;
    for (i = 0; i < ctx->FNGradL; ++i)
        ctx->FNGrad[j][i] = 0.0;
    ctx->FNGrad[j][ix] = 1.0;
    ctx->FNGradN++;         

    if (iflag) {
        j++;
        ctx->FNATyp[j] = 0;
        for (i = 0; i < ctx->FNGradL; ++i)
            ctx->FNGrad[j][i] = 0.0;
        ctx->FNGrad[j][ix] = 1.0;
        ctx->FNGradN++;         
    }

    if (deriv > 1) {
        for (i = 0; i < ctx->FNHessL; ++i)
            ctx->FNHess[j][i] = 0.0;
        ctx->FNHessN++;         
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  deriv2(op,deriv,j,val,iflag)                                            */
/*  ##                                                                      */
/*  op =  1 : addition                                                      */
/*        2 : subtraction                                                   */
/*        3 : multiplication                                                */
/*        4 : division                                                      */
/*        5 : power                                                         */
/*        6 : poisson                                                       */
/*                                                                          */
/*  return:  0  if OK                                                       */
/*           2  division by zero                                            */
/*           5  argument of log is too small (< LogMin)                     */
/*          42  if error in poisson operator                                */

int deriv2(TDAContext *ctx, int op,int deriv,int j,double *val,int iflag)
{
    register int i,ii,k;
    int j1,j2,j1u = 0,j2u = 0;
    double lg,tmp,tmp1,tmp2,tmpa,tmpb;

    if (j != ctx->FNGradN || (deriv > 1 && ctx->FNGradN != ctx->FNHessN))
        gerr_exit(ctx, 86);

    if (iflag == 0) {
        j1 = j - 1;
        j2 = j - 2;
    }
    else {
        j1u = j - 1;
        j1  = j - 2;
        j2u = j - 3;
        j2  = j - 4;
    }

    if (j2 < 0)           
        gerr_exit(ctx, 81);
 
    if (op == 0 || (ctx->FNATyp[j1] && ctx->FNATyp[j2])) {

        for (i = 0; i < ctx->FNGradL; ++i)
            ctx->FNGrad[j2][i] = 0.0;             
        ctx->FNGradN--;

        if (iflag) {
            for (i = 0; i < ctx->FNGradL; ++i)
                ctx->FNGrad[j2u][i] = 0.0;             
            ctx->FNGradN--;
        }
                 
        if (deriv > 1) {
            for (i = 0; i < ctx->FNHessL; ++i)
                ctx->FNHess[j2][i] = 0.0;               
            ctx->FNHessN--;
        }
        return(0);
    }
    switch (op) {

        case 1:                                     /* addition */

            for (i = 0; i < ctx->FNGradL; ++i)
                ctx->FNGrad[j2][i] += ctx->FNGrad[j1][i];
            ctx->FNGradN--;
               
            if (iflag) {
                for (i = 0; i < ctx->FNGradL; ++i)
                    ctx->FNGrad[j2u][i] += ctx->FNGrad[j1u][i];
                ctx->FNGradN--;
            }  
            if (deriv > 1) {
                for (i = 0; i < ctx->FNHessL; ++i)
                    ctx->FNHess[j2][i] += ctx->FNHess[j1][i];
                ctx->FNHessN--;
            }
            break;

        case 2:                                     /* subtraction */

            if (iflag == 0) {
                for (i = 0; i < ctx->FNGradL; ++i)
                    ctx->FNGrad[j2][i] -= ctx->FNGrad[j1][i];
                ctx->FNGradN--;
            }
            else {
                for (i = 0; i < ctx->FNGradL; ++i) {
                    ctx->FNGrad[j2][i] -= ctx->FNGrad[j1u][i];
                    ctx->FNGrad[j2u][i] -= ctx->FNGrad[j1][i];
                }   
                ctx->FNGradN -= 2;
            }
                 
            if (deriv > 1) {
                for (i = 0; i < ctx->FNHessL; ++i)
                    ctx->FNHess[j2][i] -= ctx->FNHess[j1][i];
                ctx->FNHessN--;
            }
            break;

        case 3:                                     /* multiplication */

            if (deriv > 1) {
                k = 0;
                for (i = 0; i < ctx->FNGradL; ++i) {
                    for (ii = 0; ii <= i; ++ii) {
                        tmp = 0.0;
                        if (ctx->FNATyp[j2] == 0)  
                            tmp += ctx->FNHess[j2][k] * val[j1];
                        if (ctx->FNATyp[j1] == 0)  
                            tmp += ctx->FNHess[j1][k] * val[j2];
                        if (ctx->FNATyp[j1] == 0 && ctx->FNATyp[j2] == 0)  
                            tmp += ctx->FNGrad[j2][i] * ctx->FNGrad[j1][ii] +
                                   ctx->FNGrad[j2][ii] * ctx->FNGrad[j1][i];
                        ctx->FNHess[j2][k] = tmp;
                        k++;
                    }
                }
                ctx->FNHessN--;
            }
            if (iflag == 0) {
                for (i = 0; i < ctx->FNGradL; ++i) {
                    tmp = 0.0;
                    if (ctx->FNATyp[j2] == 0)
                        tmp += ctx->FNGrad[j2][i] * val[j1];
                    if (ctx->FNATyp[j1] == 0)
                        tmp += ctx->FNGrad[j1][i] * val[j2];
                    ctx->FNGrad[j2][i] = tmp;
                }
                ctx->FNGradN--;
            }
            else {
                for (i = 0; i < ctx->FNGradL; ++i) {

                    ia_mul(ctx, ctx->FNGrad[j2][i],ctx->FNGrad[j2u][i],val[j1],val[j1u],&tmpa,&tmpb);
                    ia_mul(ctx, val[j2],val[j2u],ctx->FNGrad[j1][i],ctx->FNGrad[j1u][i],&tmp1,&tmp2);
                    ctx->FNGrad[j2][i]  = tmpa + tmp1;
                    ctx->FNGrad[j2u][i] = tmpb + tmp2;
                }
                ctx->FNGradN -= 2;
            }
            break;

        case 4:                                     /* division */

            if (deriv > 1) {
    
                k = 0;
                for (i = 0; i < ctx->FNGradL; ++i) {
                    for (ii = 0; ii <= i; ++ii) {
                        tmp = 0.0;
                        if (ctx->FNATyp[j2] == 0) {
                            if ((ctx->FNHess[j2][k]) != 0.0)
                                tmp += ctx->FNHess[j2][k] * val[j1];
                        }
                        if (ctx->FNATyp[j1] == 0) {
                            if ((ctx->FNHess[j1][k]) != 0.0 && (val[j2]) != 0.0)
                                tmp -= ctx->FNHess[j1][k] * val[j2];
                            if (((tmp2 = ctx->FNGrad[j1][i] * ctx->FNGrad[j1][ii] * val[j2] * 2.0)) != 0.0) {
                                if (val[j1] == 0.0)  
                                    return(2);
                                tmp += tmp2 / val[j1];
                            }
                        }
                        if (ctx->FNATyp[j1] == 0 && ctx->FNATyp[j2] == 0) {
                            tmp -= ctx->FNGrad[j2][ii] * ctx->FNGrad[j1][i];
                            tmp -= ctx->FNGrad[j2][i] * ctx->FNGrad[j1][ii];
                        }           
                        if ((tmp) != 0.0) {
                            tmp1 = val[j1] * val[j1];
                            if (tmp1 == 0.0)  
                                return(2);
                            tmp /= tmp1;
                        }
                        ctx->FNHess[j2][k] = tmp;
                        k++;
                    }
                }
                ctx->FNHessN--;
            }
            if (iflag == 0) {
                for (i = 0; i < ctx->FNGradL; ++i) {
                    tmp = 0.0;
                    if (ctx->FNATyp[j2] == 0)
                        tmp += ctx->FNGrad[j2][i] * val[j1];
                    if (ctx->FNATyp[j1] == 0)
                        tmp -= ctx->FNGrad[j1][i] * val[j2];
                    if ((tmp) != 0.0) {
                        tmp1 = val[j1] * val[j1];
                        if (tmp1 == 0.0)  
                            return(2);
                        tmp /= tmp1;
                    }
                    ctx->FNGrad[j2][i] = tmp;
                }
                ctx->FNGradN--;
            }
            else {
                for (i = 0; i < ctx->FNGradL; ++i) {

                    ia_mul(ctx, ctx->FNGrad[j2][i],ctx->FNGrad[j2u][i],val[j1],val[j1u],&tmpa,&tmpb);
                    ia_mul(ctx, val[j2],val[j2u],ctx->FNGrad[j1][i],ctx->FNGrad[j1u][i],&tmp1,&tmp2);
                    ctx->FNGrad[j2][i]  = tmpa - tmp2;
                    ctx->FNGrad[j2u][i] = tmpb - tmp1;
                    if (val[j1] <= 0.0 && val[j1u] >= 0.0)
                        return(2);
                    tmp1 = 1.0 / val[j1];
                    tmp2 = 1.0 / val[j1u];
                    ia_mul(ctx, ctx->FNGrad[j2][i],ctx->FNGrad[j2u][i],tmp2,tmp1,&tmpa,&tmpb);
                    ctx->FNGrad[j2][i]  = tmpa;
                    ctx->FNGrad[j2u][i] = tmpb;
                    ia_mul(ctx, ctx->FNGrad[j2][i],ctx->FNGrad[j2u][i],tmp2,tmp1,&tmpa,&tmpb);
                    ctx->FNGrad[j2][i]  = tmpa;
                    ctx->FNGrad[j2u][i] = tmpb;
                }
                ctx->FNGradN -= 2;
            }
            break;

        case 5:                                     /* ## ^ */
      
            if (iflag)
                return(29);

            tmp = val[j2];
            if (fabs(tmp) <= ctx->EPSI2) {
                for (i = 0; i < ctx->FNGradL; ++i)  
                    ctx->FNGrad[j2][i] = 0.0;
                ctx->FNGradN--;
                if (deriv > 1) {
                    for (i = 0; i < ctx->FNHessL; ++i)  
                        ctx->FNHess[j2][i] = 0.0;
                    ctx->FNHessN--;
                }
                break;         
            }
            if (tmp < 0.0)
                tmp1 = pow(tmp,floor(val[j1]));
            else
                tmp1 = pow(tmp,val[j1]);

            if (deriv > 1) {
                k = 0;
                for (i = 0; i < ctx->FNGradL; ++i) {
                    for (ii = 0; ii <= i; ++ii) {
                        tmp = 0.0;
                        if (ctx->FNATyp[j1] == 0) {
                            if (val[j2] <= 0.0)      
                                return(5);

                            lg = rlog(ctx, val[j2]);
                            tmp += lg * (ctx->FNHess[j1][k] + lg * ctx->FNGrad[j1][i] * ctx->FNGrad[j1][ii]);
                        }   
                        if (ctx->FNATyp[j2] == 0) {
                            if (val[j2] == 0.0)      
                                return(2);
 
                            lg = val[j1] / val[j2];
                            tmp += lg * (ctx->FNHess[j2][k] + ctx->FNGrad[j2][i] * ctx->FNGrad[j2][ii] * (val[j1] - 1.0) / val[j2]);
                        }
                        if (ctx->FNATyp[j1] == 0 && ctx->FNATyp[j2] == 0) {
                            if (val[j2] <= 0.0)      
                                return(5);

                            tmp2 = (val[j1] * rlog(ctx, val[j2]) + 1.0) / val[j2];
                            tmp += tmp2 * (ctx->FNGrad[j2][i] * ctx->FNGrad[j1][ii] +
                                           ctx->FNGrad[j1][i] * ctx->FNGrad[j2][ii]);  
                        }
                        if ((tmp) != 0.0)
                            tmp *= tmp1;
                        ctx->FNHess[j2][k] = tmp;
                        k++;
                    }
                }
                ctx->FNHessN--;
            }
            for (i = 0; i < ctx->FNGradL; ++i) {
                tmp = 0.0;
                if (ctx->FNATyp[j1] == 0 && (ctx->FNGrad[j1][i]) != 0.0) {
                    if (val[j2] <= 0.0)      
                        return(5);
 
                    tmp += ctx->FNGrad[j1][i] * rlog(ctx, val[j2]);
                }
                if (ctx->FNATyp[j2] == 0 && ((tmp2 = ctx->FNGrad[j2][i] * val[j1])) != 0.0) {
                    if (val[j2] == 0.0)      
                        return(2);
 
                    tmp += tmp2 / val[j2];
                }
                if ((tmp) != 0.0)
                    tmp *= tmp1;
                ctx->FNGrad[j2][i] = tmp;
            }
            ctx->FNGradN--;
            break;

        case 6:                         /* poisson(theta,k) */
    
            if (l_poisson(ctx, val[j2],val[j1],&tmp,&tmp1,&tmp2))
                return(42);

            val[j2] = tmp;

            if (deriv > 1) {
                k = 0;
                for (i = 0; i < ctx->FNGradL; ++i) {
                    for (ii = 0; ii <= i; ++ii) {
                        tmp = 0.0;
                        if (ctx->FNATyp[j2] == 0) {
                            tmp += tmp1 * ctx->FNHess[j2][k] +
                                   tmp2 * ctx->FNGrad[j2][i] * ctx->FNGrad[j2][ii]; 
                        }           
                        if (ctx->FNATyp[j1] == 0)  
                            return(42);

                        ctx->FNHess[j2][k] = tmp;
                        k++;
                    }
                }
                ctx->FNHessN--;
            }
            if (deriv > 0) {
                for (i = 0; i < ctx->FNGradL; ++i) {
                    tmp = 0.0;
                    if (ctx->FNATyp[j2] == 0)  
                        tmp += ctx->FNGrad[j2][i] * tmp1;    
                       
                    if (ctx->FNATyp[j1] == 0)  
                        return(42);
                       
                    ctx->FNGrad[j2][i] = tmp;
                }
                ctx->FNGradN--;
            }
            break;

        default: printfe(ctx, "ERROR IN DERIV2. CASE %d\n",op);
                 gerr_exit(ctx, 206);
    }
    ctx->FNATyp[j2] *= ctx->FNATyp[j1];
    if (iflag)
        ctx->FNATyp[j2u] *= ctx->FNATyp[j1u];
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  deriv3(op,deriv,j,val,iflag)                                            */
/*  ##                                                                      */
/*  op =  1 : negbin                                                        */
/*                                                                          */
/*  return:  0  if OK                                                       */
/*          43  if error in negbin operator                                 */

int deriv3(TDAContext *ctx, int op,int deriv,int j,double *val,int iflag)
{
    register int i,ii,k;
    int j1,j2,j3,j1u,j2u,j3u;
    double tmp;
    double dval[6];

    if (j != ctx->FNGradN || (deriv > 1 && ctx->FNGradN != ctx->FNHessN))
        gerr_exit(ctx, 86);

    if (iflag == 0) {
        j1 = j - 1;
        j2 = j - 2;
        j3 = j - 3;
    }
    else {
        j1u = j - 1;
        j1  = j - 2;
        j2u = j - 3;
        j2  = j - 4;
        j3u = j - 5;
        j3  = j - 6;
    }

    if (j3 < 0)           
        gerr_exit(ctx, 81);
 
    switch (op) {

        case 1:                         /* negbin(alpha,gamma,k) */
    
            if (ctx->FNATyp[j1] == 0)        /* derivatives only for alpha,gamma */
                return(43);

            if (l_negbin(ctx, val[j3],val[j2],val[j1],deriv,dval))
                return(43);

            val[j3] = dval[0];

            if (deriv > 1) {
                k = 0;
                for (i = 0; i < ctx->FNGradL; ++i) {
                    for (ii = 0; ii <= i; ++ii) {
                        tmp = 0.0;
                        if (ctx->FNATyp[j3] == 0) {
                            tmp += dval[1] * ctx->FNHess[j3][k] +
                                   dval[2] * ctx->FNGrad[j3][i] * ctx->FNGrad[j3][ii]; 
                        }           
                        if (ctx->FNATyp[j2] == 0) {
                            tmp += dval[3] * ctx->FNHess[j2][k] +
                                   dval[4] * ctx->FNGrad[j2][i] * ctx->FNGrad[j2][ii]; 
                                    
                            if (ctx->FNATyp[j3] == 0)  
                                tmp += dval[5] * (ctx->FNGrad[j3][i] * ctx->FNGrad[j2][ii] +
                                                  ctx->FNGrad[j3][ii] * ctx->FNGrad[j2][i]);
                        }           
                        ctx->FNHess[j3][k] = tmp;
                        k++;
                    }
                }
                ctx->FNHessN -= 2;
            }
            if (deriv > 0) {
                for (i = 0; i < ctx->FNGradL; ++i) {
                    tmp = 0.0;
                    if (ctx->FNATyp[j3] == 0)  
                        tmp += ctx->FNGrad[j3][i] * dval[1]; 
                       
                    if (ctx->FNATyp[j2] == 0)  
                        tmp += ctx->FNGrad[j2][i] * dval[3]; 
                       
                    ctx->FNGrad[j3][i] = tmp;
                }
                ctx->FNGradN -= 2;
            }
            break;

        default: printfe(ctx, "ERROR IN DERIV3. CASE %d\n",op);
                 gerr_exit(ctx, 207);
    }
    ctx->FNATyp[j2] *= ctx->FNATyp[j1];
    ctx->FNATyp[j3] *= ctx->FNATyp[j2];
    if (iflag) {
        ctx->FNATyp[j2u] *= ctx->FNATyp[j1u];
        ctx->FNATyp[j3u] *= ctx->FNATyp[j2u];
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  derivf1(op,deriv,j,val,iflag)                                           */
/*                                                                          */
/*  op =  1 :   abs                                                         */
/*        3 :   sqrt                                                        */
/*        4 :   exp                                                         */
/*        5 :   eexp                                                        */
/*        6 :   log                                                         */
/*        7 :   sin                                                         */
/*        8 :   cos                                                         */
/*        9 :   lgam                                                        */
/*       10 :   ndf                                                         */
/*       11 :   nd                                                          */
/*       12 :   ndi                                                         */
/*                                                                          */
/*  return:  0  if OK                                                       */
/*           2  division by zero                                            */
/*           6  negative argument in sqrt()                                 */
/*           7  negative or zero argument for lgam, digam                   */

int derivf1(TDAContext *ctx, int op,int deriv,int j,double *val,int iflag)
{
    register int i,ii,k;
    int err,j1,j1u = 0;       
    double tmp = 0.0,tmp1,tmp2,tmp3,tmpa,tmpb;

    if (j != ctx->FNGradN || (deriv > 1 && ctx->FNGradN != ctx->FNHessN))
        gerr_exit(ctx, 87);
   
    if (iflag == 0)
        j1 = j - 1;
    else {
        j1 = j - 2;
        j1u = j - 1;
    }
    if (j1 < 0)
        gerr_exit(ctx, 82);

    if (op == 0 || ctx->FNATyp[j1] == 1) {
                            
        if (deriv > 0) {
            for (i = 0; i < ctx->FNGradL; ++i)
                ctx->FNGrad[j1][i] = 0.0;             
        }
        if (deriv > 1) {
            for (i = 0; i < ctx->FNHessL; ++i)
                ctx->FNHess[j1][i] = 0.0;             
        }
        return(0);
    }                 

    switch (op) {

        case 1:                                     /* abs */

            if (iflag)
                return(29);

            if (val[j1] < 0.0) {
                if (deriv > 0) {
                    for (i = 0; i < ctx->FNGradL; ++i)
                        ctx->FNGrad[j1][i] *= -1.0;           
                }
                if (deriv > 1) {
                    for (i = 0; i < ctx->FNHessL; ++i)
                        ctx->FNHess[j1][i] *= -1.0;             
                }
            }
            break;

        case 3:                                     /* sqrt */

            if (iflag)
                return(29);

            if (deriv > 1) {
                k = 0;
                for (i = 0; i < ctx->FNGradL; ++i) {
                    for (ii = 0; ii <= i; ++ii) {
                        tmp = 0.0;
                        if (ctx->FNATyp[j1] == 0) {
                            tmp = ctx->FNHess[j1][k];
                            if (((tmp1 = ctx->FNGrad[j1][i] * ctx->FNGrad[j1][ii])) != 0.0) {
                                if (val[j1] == 0.0)  
                                    return(2);
                                tmp -= tmp1 / (2.0 * val[j1]);
                            }
                            if ((tmp) != 0.0) {
                                if (val[j1] <= 0.0)         
                                    return(2);
                                tmp /= (2.0 * sqrt(val[j1]));
                            }   
                        }
                        ctx->FNHess[j1][k] = tmp;               
                        k++;
                    }
                }
            }
            for (i = 0; i < ctx->FNGradL; ++i) {
                tmp = 0.0;
                if (ctx->FNATyp[j1] == 0 && (ctx->FNGrad[j1][i]) != 0.0) {
                    if (val[j1] <= 0.0)      
                        return(2);
                    tmp = 0.5 * ctx->FNGrad[j1][i] / sqrt(val[j1]);
                }
                ctx->FNGrad[j1][i] = tmp;            
            }
            break;

        case 4:                                     /* exp */

            if (iflag == 0)
                tmp = rexp(ctx, val[j1]);

            if (deriv > 1) {
                k = 0;
                for (i = 0; i < ctx->FNGradL; ++i) {
                    tmp1 = ctx->FNGrad[j1][i];
                    for (ii = 0; ii <= i; ++ii) {
                        tmp2 = ctx->FNGrad[j1][ii];
                        ctx->FNHess[j1][k] = tmp * (tmp1 * tmp2 + ctx->FNHess[j1][k]);
                        k++;
                    }   
                }
            }
            if (deriv > 0) {
                if (iflag == 0) {
                    for (i = 0; i < ctx->FNGradL; ++i)  
                        ctx->FNGrad[j1][i] *= tmp;       
                }
                else {
                    for (i = 0; i < ctx->FNGradL; ++i) {
                        tmpa = rexp(ctx, val[j1]);
                        tmpb = rexp(ctx, val[j1u]);
                        ia_mul(ctx, ctx->FNGrad[j1][i],ctx->FNGrad[j1u][i],tmpa,tmpb,&tmp1,&tmp2);
                        ctx->FNGrad[j1][i] = tmp1;
                        ctx->FNGrad[j1u][i] = tmp2;
                    }
                }
            }
            break;

        case 5:                                     /* ##  eexp */

            if (iflag)
                return(29);

            tmp = rexp(ctx, val[j1]);
 
            tmp1 = tmp / (1.0 + tmp);
            tmp2 = tmp1 * (1.0 - tmp1);
            tmp3 = tmp2 * (1.0 - 2.0 * tmp1);

            if (deriv > 1) {
                k = 0;
                for (i = 0; i < ctx->FNGradL; ++i) {
                    for (ii = 0; ii <= i; ++ii) {
                        ctx->FNHess[j1][k] = tmp2 * ctx->FNHess[j1][k] +
                            tmp3 * ctx->FNGrad[j1][i] * ctx->FNGrad[j1][ii];
                        k++; 
                    }
                }
            }
            if (deriv > 0) {
                for (i = 0; i < ctx->FNGradL; ++i)  
                    ctx->FNGrad[j1][i] *= tmp2;
            }
            break;

        case 6:                                     /* log */

            if (deriv > 1) {
                k = 0;
                for (i = 0; i < ctx->FNGradL; ++i) {
                    for (ii = 0; ii <= i; ++ii) {
                        tmp = 0.0;
                        if ((ctx->FNHess[j1][k]) != 0.0) {
                            if (val[j1] == 0)  
                                return(2);
 
                            tmp += ctx->FNHess[j1][k] / val[j1];
                        }
                        if (((tmp2 = ctx->FNGrad[j1][i] * ctx->FNGrad[j1][ii])) != 0.0) {
                            tmp1 = val[j1] * val[j1];
                            if (tmp1 == 0)  
                                return(2);
 
                            tmp -= tmp2 / tmp1;        
                        }
                        ctx->FNHess[j1][k] = tmp;               
                        k++; 
                    }
                }
            }
            if (deriv > 0) {
                if (iflag == 0) {
                    for (i = 0; i < ctx->FNGradL; ++i) {
                        tmp = 0.0;

                        if ((ctx->FNGrad[j1][i]) != 0.0) {
                            if (val[j1] == 0.0)      
                                return(2);
 
                            tmp = ctx->FNGrad[j1][i] / val[j1];
                        }
                        ctx->FNGrad[j1][i] = tmp;       
                    }
                }
                else {
                    for (i = 0; i < ctx->FNGradL; ++i) {
                        if (val[j1] <= 0.0 && val[j1u] >= 0.0)
                            return(2);
   
                        ia_mul(ctx, 1.0,1.0,1.0 / val[j1u],1.0 / val[j1],&tmpa,&tmpb);
                        ia_mul(ctx, ctx->FNGrad[j1][i],ctx->FNGrad[j1u][i],tmpa,tmpb,&tmp1,&tmp2);
                        ctx->FNGrad[j1][i] = tmp1;
                        ctx->FNGrad[j1u][i] = tmp2;
                    }
                }
            }
            break;

        case 7:                                     /* sin */

            if (iflag == 0)
                tmp = cos(val[j1]);

            if (deriv > 1) {
                k = 0;
                tmp1 = sin(val[j1]);
                for (i = 0; i < ctx->FNGradL; ++i) {
                    for (ii = 0; ii <= i; ++ii) {
                        ctx->FNHess[j1][k] = tmp * ctx->FNHess[j1][k] -
                                  tmp1 * ctx->FNGrad[j1][i] * ctx->FNGrad[j1][ii];
                        k++; 
                    }
                }
            }
            if (deriv > 0) {

                if (iflag == 0) {
                    for (i = 0; i < ctx->FNGradL; ++i)  
                        ctx->FNGrad[j1][i] *= tmp;       
                }
                else {
                    for (i = 0; i < ctx->FNGradL; ++i) {
                        ia_cos(ctx, val[j1],val[j1u],&tmpa,&tmpb);
                        ia_mul(ctx, ctx->FNGrad[j1][i],ctx->FNGrad[j1u][i],tmpa,tmpb,&tmp1,&tmp2);
                        ctx->FNGrad[j1][i] = tmp1;
                        ctx->FNGrad[j1u][i] = tmp2;
                    }
                }
            }
            break;

        case 8:                                     /* cos */

            if (iflag)
                return(29);

            tmp = -sin(val[j1]);

            if (deriv > 1) {
                k = 0;
                tmp1 = cos(val[j1]);
                for (i = 0; i < ctx->FNGradL; ++i) {
                    for (ii = 0; ii <= i; ++ii) {
                        ctx->FNHess[j1][k] = tmp * ctx->FNHess[j1][k] -
                                  tmp1 * ctx->FNGrad[j1][i] * ctx->FNGrad[j1][ii];
                        k++; 
                    }
                }
            }
            if (deriv > 0) {
                for (i = 0; i < ctx->FNGradL; ++i)  
                    ctx->FNGrad[j1][i] *= tmp;       
            }
            break;

        case 9:                                     /* lgam */

            if (iflag)
                return(29);

            tmp = val[j1];
            if (tmp <= 0.0)
                return(7);
            tmp1 = digam(ctx, tmp,&err);
            if (err)
                return(7);

            if (deriv > 1) {
                tmp2 = trigam(ctx, tmp,&err);
                if (err)
                    return(7);

                k = 0;
                for (i = 0; i < ctx->FNGradL; ++i) {
                    for (ii = 0; ii <= i; ++ii) {
                        ctx->FNHess[j1][k] = tmp1 * ctx->FNHess[j1][k] +
                                  tmp2 * ctx->FNGrad[j1][i] * ctx->FNGrad[j1][ii];
                        k++; 
                    }
                }
            }
            if (deriv > 0) {
                for (i = 0; i < ctx->FNGradL; ++i)  
                    ctx->FNGrad[j1][i] *= tmp1;       
            }
            break;

        case 10:                                    /* ndf */

            if (iflag)
                return(29);

            tmp1 = dnf(ctx, val[j1]);

            if (deriv > 1) {
                k = 0;
                for (i = 0; i < ctx->FNGradL; ++i) {
                    for (ii = 0; ii <= i; ++ii) {
                        tmp = 0.0;
                        if ((ctx->FNHess[j1][k]) != 0.0)
                            tmp -= val[j1] * ctx->FNHess[j1][k];
                        if (((tmp2 = ctx->FNGrad[j1][i] * ctx->FNGrad[j1][ii])) != 0.0)
                            tmp += tmp2 * (val[j1] * val[j1] - 1.0);
                        if ((tmp) != 0.0)
                            tmp *= tmp1;

                        ctx->FNHess[j1][k] = tmp;
                        k++; 
                    }
                }
            }
            for (i = 0; i < ctx->FNGradL; ++i) {
                if ((ctx->FNGrad[j1][i]) != 0.0)
                    ctx->FNGrad[j1][i] *= -val[j1] * tmp1;
            }
            break;

        case 11:                                    /* nd */
        case 12:                                    /* ndi */

            if (iflag)
                return(29);

            tmp1 = dnf(ctx, val[j1]);
            if (op == 12) {
                if (deriv > 1)
                    return(44);
                if (tmp1 == 0.0)
                    return(2);
                tmp1 = 1.0 / tmp1;
            }
       
            if (deriv > 1) {
                k = 0;
                for (i = 0; i < ctx->FNGradL; ++i) {
                    for (ii = 0; ii <= i; ++ii) {
                        tmp = ctx->FNHess[j1][k];
                        if (((tmp2 = ctx->FNGrad[j1][i] * ctx->FNGrad[j1][ii])) != 0.0)
                            tmp -= tmp2 * val[j1];
                        if ((tmp) != 0.0)
                            tmp *= tmp1;

                        ctx->FNHess[j1][k] = tmp;
                        k++; 
                    }
                }
            }
            for (i = 0; i < ctx->FNGradL; ++i) {
                if ((ctx->FNGrad[j1][i]) != 0.0)
                    ctx->FNGrad[j1][i] *= tmp1;
            }
            break;

        default: printfe(ctx, "ERROR IN DERIVF1. CASE %d\n",op);
                 gerr_exit(ctx, 208);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  derivf2(op,deriv,j,val,dval)                                            */
/*  ##                                                                      */
/*  op  =  0    : set zero                                                  */
/*         1    : incomplete gamma integral                                 */
/*                                                                          */
/*  return:  0  if OK                                                       */
/*          32  if error in calculation of incomplete gamma integral.       */  
/*          42  if error in evaluation of log of poisson distribution.      */  

int derivf2(TDAContext *ctx, int op,int deriv,int j,double *val,double *dval)
{
    (void)val;        /* unused: the signature is shared */
    register int i,ii,k;
    int j1,j2;
    double tmp;

    if (j != ctx->FNGradN || (deriv > 1 && ctx->FNGradN != ctx->FNHessN))
        gerr_exit(ctx, 86);

    j1 = j - 1;
    j2 = j - 2;

    if (j2 < 0)           
        gerr_exit(ctx, 81);
 
    if (op == 0 || (ctx->FNATyp[j1] && ctx->FNATyp[j2])) {

        for (i = 0; i < ctx->FNGradL; ++i)
            ctx->FNGrad[j2][i] = 0.0;             
        ctx->FNGradN--;

        if (deriv > 1) {
            for (i = 0; i < ctx->FNHessL; ++i)
                ctx->FNHess[j2][i] = 0.0;               
            ctx->FNHessN--;
        }
        return(0);
    }
    switch (op) {

        case 1:                         /* incomplete gamma integral */
   
            if (deriv > 1) {
                k = 0;
                for (i = 0; i < ctx->FNGradL; ++i) {
                    for (ii = 0; ii <= i; ++ii) {
                        tmp = 0.0;
                        if (ctx->FNATyp[j2] == 0) {
                            tmp += dval[1] * ctx->FNHess[j2][k] +
                                   dval[2] * ctx->FNGrad[j2][i] * ctx->FNGrad[j2][ii]; 
                        }           
                        if (ctx->FNATyp[j1] == 0) {
                            tmp += dval[3] * ctx->FNHess[j1][k] +
                                   dval[4] * ctx->FNGrad[j1][i] * ctx->FNGrad[j1][ii]; 
                                    
                            if (ctx->FNATyp[j2] == 0) {

                                tmp += dval[5] * (ctx->FNGrad[j1][i] * ctx->FNGrad[j2][ii] +
                                                   ctx->FNGrad[j2][i] * ctx->FNGrad[j1][ii]);
                            }
                        }
                        ctx->FNHess[j2][k] = tmp;
                        k++;
                    }
                }
                ctx->FNHessN--;
            }
            if (deriv > 0) {
                for (i = 0; i < ctx->FNGradL; ++i) {
                    tmp = 0.0;
                    if (ctx->FNATyp[j2] == 0)  
                        tmp += ctx->FNGrad[j2][i] * dval[1];
                       
                    if (ctx->FNATyp[j1] == 0)  
                        tmp += ctx->FNGrad[j1][i] * dval[3];
                       
                    ctx->FNGrad[j2][i] = tmp;
                }
                ctx->FNGradN--;
            }
            break;

        default: printfe(ctx, "ERROR IN DERIVF2. CASE %d\n",op);
                 gerr_exit(ctx, 209);
    }
    ctx->FNATyp[j2] *= ctx->FNATyp[j1];
    return(0);
}
      
/* ------------------------------------------------------------------------ */
/*  derivf3(op,deriv,j,val,iflag)                                           */
/*                                                                          */
/*  op  1 : tr()                                                            */
/*      2 : bivn()                                                          */
/*                                                                          */  
/*  return:  0  if OK                                                       */
/*           2  division by zero                                            */
/*           5  argument of log is too small (< LogMin)                     */
/*          34  rho for bivn() is out of range.                             */

int derivf3(TDAContext *ctx, int op,int deriv,int j,double *val,int iflag)
{
    register int i,k,ii;      
    int j1,j2,j3;
    double x,y,a,aa,b,tmp,r,xy,yx,phix,phiy,phixy,phiyx,pphixy,pphiyx,phi2;
    double tmp11,tmp12,tmp13,tmp21,tmp22,tmp23,tmp31,tmp32,tmp33;

    if (j != ctx->FNGradN || (deriv > 1 && ctx->FNGradN != ctx->FNHessN))
        gerr_exit(ctx, 87);

    j1 = j - 1;
    j2 = j - 2;
    j3 = j - 3;

    if (j3 < 0)
        gerr_exit(ctx, 83);

    switch (op) {

        case 1:                                     /* tr(x,a,b) */

            if (iflag)
                return(29);

            x = val[j3];
            a = val[j2];
            b = val[j1];

            if (deriv > 0) {
                if (x < a || x > b) {
                    for (i = 0; i < ctx->FNGradL; ++i)
                        ctx->FNGrad[j3][i] = 0.0;           
                }
                ctx->FNGradN -= 2;
            }
            if (deriv > 1) {
                if (x < a || x > b) {
                    for (i = 0; i < ctx->FNHessL; ++i)
                        ctx->FNHess[j3][i] = 0.0;           
                }
                ctx->FNHessN -= 2;
            }
            break;

        case 2:                                     /* bivn(ctx, x1,x2,r) */
            if (iflag)
                return(29);

            x = val[j3];
            y = val[j2];
            r = val[j1];



            if (r <= -1.0 || r >= 1.0)
                return(34);

            aa = 1.0 - r * r;
            a = sqrt(aa);

            xy = (x - r * y) / a;
            yx = (y - r * x) / a;

            phix = dnf(ctx, x);
            phiy = dnf(ctx, y);
            phixy = dnf(ctx, xy);
            phiyx = dnf(ctx, yx);
            pphixy = cdnf(ctx, xy);
            pphiyx = cdnf(ctx, yx);
            phi2 = phix * phiyx / a;

            if (deriv > 1) {

   
                tmp11 = -phix * (r * phiyx / a + x * pphiyx);
                tmp12 =  phi2; 
                tmp13 =  -phi2 * xy / a;
                tmp21 =  tmp12;    

                tmp22 = -phiy * (r * phixy / a + y * pphixy);
                tmp23 =  -phi2 * yx / a;

                tmp31 =  tmp13;
                tmp32 =  tmp23;
                tmp33 =  phix * phiyx * ((xy * yx + r) / aa) / a;

                k = 0;
                for (i = 0; i < ctx->FNGradL; ++i) {
                    for (ii = 0; ii <= i; ++ii) {
                        tmp = 0.0;
                        if (ctx->FNATyp[j3] == 0) {        
                            tmp += phix * phiyx * ctx->FNHess[j3][k];  
                            tmp += ctx->FNGrad[j3][i] * (tmp11 * ctx->FNGrad[j3][i] +
                                                    tmp12 * ctx->FNGrad[j2][i] +
                                                    tmp13 * ctx->FNGrad[j1][i]); 
                        }     
                        if (ctx->FNATyp[j2] == 0) {        
                            tmp += phiy * phixy * ctx->FNHess[j2][k];  
                            tmp += ctx->FNGrad[j2][i] * (tmp21 * ctx->FNGrad[j3][i] +
                                                    tmp22 * ctx->FNGrad[j2][i] +
                                                    tmp23 * ctx->FNGrad[j1][i]); 
                        }      
                        if (ctx->FNATyp[j1] == 0) {
                            tmp += phi2 * ctx->FNHess[j1][k];
                            tmp += ctx->FNGrad[j1][i] * (tmp31 * ctx->FNGrad[j3][i] +
                                                    tmp32 * ctx->FNGrad[j2][i] +
                                                    tmp33 * ctx->FNGrad[j1][i]); 
                        }
                        ctx->FNHess[j3][k] = tmp;
                        k++;
                    }
                }
                ctx->FNHessN -= 2;
            }
            if (deriv > 0) {
                for (i = 0; i < ctx->FNGradL; ++i) {
                    tmp = 0.0;
                    if (ctx->FNATyp[j3] == 0)      
                        tmp += ctx->FNGrad[j3][i] * phix * pphiyx;

                    if (ctx->FNATyp[j2] == 0)      
                        tmp += ctx->FNGrad[j2][i] * phiy * pphixy;

                    if (ctx->FNATyp[j1] == 0)       
                        tmp += ctx->FNGrad[j1][i] * phi2;             
   
                    ctx->FNGrad[j3][i] = tmp;
                }
                ctx->FNGradN -= 2;
            }
            break;

        default: printfe(ctx, "ERROR IN DERIVF3. CASE %d\n",op);
                 gerr_exit(ctx, 210);
    }
    ctx->FNATyp[j2] *= ctx->FNATyp[j1];
    ctx->FNATyp[j3] *= ctx->FNATyp[j2];
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  get_integrand(x,err)                                                    */
/*  ##                                                                      */
/*  get integrand for user defined function and derivatives.                */
/*                                                                          */
/*  uses parameters set in v_eval1()                                        */
/*  ESCnt,ESTyp,ESVal = parser stack                                        */
/*  ESCase = current case                                                   */
/*  ESIVP = index to val[]                                                  */
/*  ESDeriv = 0, 1, or 2 for derivatives                                    */
/*  ESDerivI = index of argument for derivatives                            */

double get_integrand(TDAContext *ctx, double x,int *err)
{
    double f;

    ctx->NINTTVAR = x;

    *err = v_eval1(ctx, ctx->ESCase,ctx->ESCnt,ctx->ESTyp,ctx->ESVal,ctx->ESIdx,&f,0,ctx->ESDeriv,0,ctx->ESIVP,0);
    if (*err)
        return(0.0);

    if (ctx->ESDeriv == 1)
        return(ctx->FNGrad[ctx->ESIVP][ctx->ESDerivI]);
    else if (ctx->ESDeriv == 2)
        return(ctx->FNHess[ctx->ESIVP][ctx->ESDerivI]);
    return(f);
}
   
/* ------------------------------------------------------------------------ */
/*  ia_op(typ,j,buf)                                                        */
/*                                                                          */
/*  Basic interval operations. Return j >= 0 if OK, -1 if error.            */

int ia_op(TDAContext *ctx, int typ,int j,double *buf)
{
    register int ja,jb,jc,jd;
    double a,b,c,d;

    ja = j - 4;
    jb = j - 3;
    jc = j - 2;
    jd = j - 1;

    switch (typ) {
        case 1:                         /* addition */
            buf[ja] += buf[jc];
            buf[jb] += buf[jd];
            j -= 2;
            break;

        case 2:                         /* subtraction */
            buf[ja] -= buf[jd];
            buf[jb] -= buf[jc];
            j -= 2;
            break;

        case 3:                         /* division */
            a = buf[ja];
            b = buf[jb];
            c = buf[jc];
            d = buf[jd];
            if (c <= 0.0 && d >= 0.0)   
                return(-1);
               
            buf[jc] = 1.0 / d;
            buf[jd] = 1.0 / c;          /* continue with multiplication */
            TDA_FALLTHROUGH;

        case 4:                         /* multiplication */
            a = buf[ja];
            b = buf[jb];
            c = buf[jc];
            d = buf[jd];

            ia_mul(ctx, a,b,c,d,&buf[ja],&buf[jb]);
            j -= 2;
            break;

        case 5:                         /* exp */
            buf[jc] = rexp(ctx, buf[jc]);
            buf[jd] = rexp(ctx, buf[jd]);
            break;

        case 6:                         /* log */
            if (buf[jc] <= 0.0 || buf[jd] <= 0.0)
                return(-1);
            buf[jc] = rlog(ctx, buf[jc]);
            buf[jd] = rlog(ctx, buf[jd]);
            break;

        case 10:                        /* sine */
            a = buf[jc];
            b = buf[jd];
            ia_sin(ctx, a,b,&buf[jc],&buf[jd]);
            break;

        case 11:                        /* cosine */
            a = buf[jc];
            b = buf[jd];
            ia_cos(ctx, a,b,&buf[jc],&buf[jd]);
            break;

        case 53:                         /* le */
            if (buf[jd] < buf[ja]) {
                buf[ja] = 0.0;
                buf[jb] = 0.0;
            }
            else if (buf[jb] <= buf[jc]) {
                buf[ja] = 1.0;
                buf[jb] = 1.0;
            }
            else {
                buf[ja] = 0.0;
                buf[jb] = 1.0;
            }
            j -= 2;
            break;

        case 54:                         /* lt */
            if (buf[jd] <= buf[ja]) {
                buf[ja] = 0.0;
                buf[jb] = 0.0;
            }
            else if (buf[jb] < buf[jc]) {
                buf[ja] = 1.0;
                buf[jb] = 1.0;
            }
            else {
                buf[ja] = 0.0;
                buf[jb] = 1.0;
            }
            j -= 2;
            break;

        case 55:                         /* ge */
            if (buf[jb] < buf[jc]) {
                buf[ja] = 0.0;
                buf[jb] = 0.0;
            }
            else if (buf[jd] <= buf[ja]) {
                buf[ja] = 1.0;
                buf[jb] = 1.0;
            }
            else {
                buf[ja] = 0.0;
                buf[jb] = 1.0;
            }
            j -= 2;
            break;

        case 56:                         /* gt */
            if (buf[jb] <= buf[jc]) {
                buf[ja] = 0.0;
                buf[jb] = 0.0;
            }
            else if (buf[jd] < buf[ja]) {
                buf[ja] = 1.0;
                buf[jb] = 1.0;
            }
            else {
                buf[ja] = 0.0;
                buf[jb] = 1.0;
            }
            j -= 2;
            break;

        case 57:                         /* eq */

            if (buf[jb] < buf[jc] || buf[jd] < buf[ja]) {
                buf[ja] = 0.0;
                buf[jb] = 0.0;
            }
            else if (buf[ja] == buf[jb] && buf[jb] == buf[jc] && buf[jc] == buf[jd]) {  
                buf[ja] = 1.0;
                buf[jb] = 1.0;
            }
            else {
                buf[ja] = 0.0;
                buf[jb] = 1.0;
            }
            j -= 2;
            break;

        case 58:                         /* ne */

            if (buf[jb] < buf[jc] || buf[jd] < buf[ja]) {
                buf[ja] = 1.0;
                buf[jb] = 1.0;
            }
            else if (buf[ja] == buf[jb] && buf[jb] == buf[jc] && buf[jc] == buf[jd]) {  
                buf[ja] = 0.0;
                buf[jb] = 0.0;
            }
            else {
                buf[ja] = 0.0;
                buf[jb] = 1.0;
            }
            j -= 2;
            break;

        case 60:                         /* and */

            if ((buf[ja] == buf[jb] && buf[jb] == 0.0) || 
                                      (buf[jc] == buf[jd] && buf[jd] == 0.0)) {  
                buf[ja] = 0.0;
                buf[jb] = 0.0;
            }
            else if (buf[ja] == buf[jb] && buf[jb] == buf[jc] && buf[jc] == buf[jd] && buf[jd] == 1.0) {
                buf[ja] = 1.0;
                buf[jb] = 1.0;
            }
            else {
                buf[ja] = 0.0;
                buf[jb] = 1.0;
            }
            j -= 2;
            break;

        case 61:                         /* or */

            if ((buf[ja] == buf[jb] && buf[jb] == 1.0) || 
                                      (buf[jc] == buf[jd] && buf[jd] == 1.0)) {  
                buf[ja] = 1.0;
                buf[jb] = 1.0;
            }
            else if (buf[ja] == buf[jb] && buf[jb] == buf[jc] && buf[jc] == buf[jd] && buf[jd] == 0.0) {
                buf[ja] = 0.0;
                buf[jb] = 0.0;
            }
            else {
                buf[ja] = 0.0;
                buf[jb] = 1.0;
            }
            j -= 2;
            break;


        default:
            gerr_exit(ctx, 90);
    }
    return(j);
}

/* ------------------------------------------------------------------------ */
/*  ia_mul(a,b,c,d,ra,rb)   Multiplication [a,b] * [c,d] -> [ra,rb]         */

void ia_mul(TDAContext *ctx, double a,double b,double c,double d,double *ra,double *rb)
{
    (void)ctx;        /* unused: the signature is shared */
    double ad,bc,bd;

    *ra = *rb = a * c;
    ad = a * d;
    bc = b * c;
    bd = b * d;
    if (*ra > ad)
        *ra = ad;
    if (*ra > bc)
        *ra = bc;
    if (*ra > bd)
        *ra = bd;

    if (*rb < ad)
        *rb = ad;
    if (*rb < bc)
        *rb = bc;
    if (*rb < bd)
        *rb = bd;
}

/* ------------------------------------------------------------------------ */
/*  ia_iv(a,b,c,d,ra,rb)    Interval  ra = min(a,c), rb = max(b,d)          */

void ia_iv(TDAContext *ctx, double a,double b,double c,double d,double *ra,double *rb)
{
    (void)ctx;        /* unused: the signature is shared */
    if (a < c)
        *ra = a;
    else
        *ra = c;

    if (b > d)
        *rb = b;
    else
        *rb = d;
}                     

/* ------------------------------------------------------------------------ */
/*  ia_sin(a,b,ra,rb)   sin([a,b]) -> [ra,rb]                               */

void ia_sin(TDAContext *ctx, double a,double b,double *ra,double *rb)
{
    double tmp;

    if (a == b)  
        *ra = *rb = sin(a);

    else {
        *ra = -1.0;
        *rb =  1.0;
                 
        if (b - a < 2.0 * Pi) {

            tmp = floor(a / (2.0 * Pi));
            tmp *= 2.0 * Pi;
            a -= tmp;
            b -= tmp;

            if (a <= Pi / 2.0) {
                if (b <= Pi / 2.0) {
                    *ra = sin(a);
                    *rb = sin(b);
                }
                else if (b < 1.5 * Pi) {
                    *ra = dmin(ctx, sin(a),sin(b));
                    *rb = 1.0;
                }
            }
            else if (a <= 1.5 * Pi) {
                if (b <= 1.5 * Pi) {
                    *ra = sin(b);
                    *rb = sin(a);
                }
                else {
                    *ra = -1.0;
                    *rb = dmax(ctx, sin(a),sin(b));
                }
            }
            else {
                if (b <= 2.5 * Pi) {
                    *ra = sin(a);
                    *rb = sin(b);
                }
                else if (b <= 3.5 * Pi) {
                    *ra = dmin(ctx, sin(a),sin(b));
                    *rb = 1.0;
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  ia_cos(a,b,ra,rb)   cos([a,b]) -> [ra,rb]                               */

void ia_cos(TDAContext *ctx, double a,double b,double *ra,double *rb)
{
    double tmp;

    if (a == b)  
        *ra = *rb = cos(a);

    else {
        *ra = -1.0;
        *rb =  1.0;
                 
        if (b - a < 2.0 * Pi) {

            tmp = floor(a / (2.0 * Pi));
            tmp *= 2.0 * Pi;
            a -= tmp;
            b -= tmp;

            if (a <= Pi) {
                if (b <= Pi) {
                    *ra = cos(b);
                    *rb = cos(a);
                }
                else if (b < 2.0 * Pi) {
                    *ra = -1.0;
                    *rb = dmax(ctx, cos(a),cos(b));
                }
            }
            else if (a <= 2.0 * Pi) {
                if (b <= 2.0 * Pi) {
                    *ra = cos(a);
                    *rb = cos(b);
                }
                else if (b <= 3.0 * Pi) {
                    *ra = dmin(ctx, cos(a),cos(b));
                    *rb = 1.0;
                }
            }
        }
    }
}


