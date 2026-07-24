/****************************************************************************/
/*  t_int                                                                   */
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
#include "t_alloc.h"
#include "t_gmin.h"
#include "t_fnrtb.h"
#include "t_gf.h"
#include "t_eval.h"
#include "t_eval3.h"
#include "tda_context.h"

/*  functions in t_int.c */

int niset(TDAContext *ctx);
void ni_info(TDAContext *ctx);
int t_int(TDAContext *ctx);
int prn_int(TDAContext *ctx, int m,double rerr);
double ni_gen(TDAContext *ctx, double a,double b,int ftyp,int nhp,int htyp,int *err);
double ni_func(TDAContext *ctx, double x,int *err,int ftyp);
int qsniff(TDAContext *ctx, double a,double b,double rerr,double *val,int *nf,int *nfu,int ftyp);
double qsniff1(TDAContext *ctx, double u,double v,int ftyp,int *err);
int qng(TDAContext *ctx, double a,double b,double epsabs,double epsrel,double *res, double *abserr,int *nf,int ftyp);
int qsub(TDAContext *ctx, double a,double b,double rerr,double *val,double *relerr,int *nfs, int *rtyp,int ftyp);
int qsubf(TDAContext *ctx, double a,double b,double rerr,double *val,double *val1,int *nf,int ftyp);
int qtrap(TDAContext *ctx, double a,double b,double rerr,double *val,int *nf,int ftyp);
int qsimp(TDAContext *ctx, double a,double b,double rerr,double *val,int *nf,int ftyp);
double trapzd(TDAContext *ctx, double a,double b,int n,int *nf,int *err,int ftyp);


/* parameter for Hermite integration */

           
/* ------------------------------------------------------------------------ */
/*  niset()         Set options for numerical integration                   */
/*                  Return 0 if OK, otherwise -1.                           */

int niset(TDAContext *ctx)
{
    int err;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);
         
    if (parm(ctx, ctx->CmdBuf + 5,3,1))   /* get parameters */
        goto NISETFin;

    if (ctx->PMRHSI < 1 || ctx->PMRHSI > 5) {
        p_err(ctx, -1,1);
        goto NISETFin;
    }
    ctx->NINTMETH = ctx->PMRHSI;
    ctx->NINTRERR = ctx->PMRERR;
    ctx->NINTAERR = ctx->PMAERR;
    printf1(ctx, "Numerical integration with method %d. Relative error: %lg (%lg).\n",
                                                    ctx->NINTMETH,ctx->NINTRERR,ctx->NINTAERR);
    err = 0;

NISETFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ni_info()       Print method about numerical integration                */

void ni_info(TDAContext *ctx)
{
    if (ctx->NINTMUsed >= 1 && ctx->NINTMUsed <= 5)  
        printf1(ctx, "Numerical integration with method %d. Relative error: %lg (%lg).\n",
                                                    ctx->NINTMETH,ctx->NINTRERR,ctx->NINTAERR);
}

/* ------------------------------------------------------------------------ */
/*  t_int()         Numerical integreation.                                 */
/*                  Return 0 if OK, otherwise -1.                           */

int t_int(TDAContext *ctx)
{
    int err = 0,nf = 0,nfu = 0,r = 0,rtyp = 0;
    double res = 0.0,aerr = 0.0;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "Numerical integration. Current memory: %d bytes.\n",ctx->MemReq);
         
    if (parm(ctx, ctx->CmdBuf + 3,9,1))   /* get parameters */
        goto TIFin;

    if (ctx->PMFTYP5) {
        p_err(ctx, -40,1);
        goto TIFin;
    }
    if (ctx->FNArgN != 1) {
        printf1(ctx, "Error: function should contain exactly one argument.\n");
        goto TIFin;
    }
    printf1(ctx, "Function argument: %s\n",ctx->FNArgDef[0]);

    if (ctx->PMXYFlg == 0 || ctx->PMX >= ctx->PMY) {
        printf1(ctx, "Error: need a valid integration interval.\n");
        goto TIFin;
    }    
    printf1(ctx, "Integration interval: %lg -- %lg\n",ctx->PMX,ctx->PMY);

    if (prn_int(ctx, ctx->NINTMETH,ctx->NINTRERR))       /* print method */
        goto TIFin;

    if (ctx->NINTMETH == 4)  
        printf1(ctx, "Absolute error: %12.5e\n",ctx->NINTAERR);           

    switch (ctx->NINTMETH) {
        case 1:     r = qng(ctx, ctx->PMX,ctx->PMY,ctx->NINTAERR,ctx->NINTRERR,&res,&aerr,&nf,0);
                    break;
        case 2:     r = qsub(ctx, ctx->PMX,ctx->PMY,ctx->NINTRERR,&res,&aerr,&nf,&rtyp,0);
                    break;
        case 3:     r = qsniff(ctx, ctx->PMX,ctx->PMY,ctx->NINTRERR,&res,&nf,&nfu,0);
                    break;
        case 4:     r = qtrap(ctx, ctx->PMX,ctx->PMY,ctx->NINTRERR,&res,&nf,0);
                    break;
        case 5:     r = qsimp(ctx, ctx->PMX,ctx->PMY,ctx->NINTRERR,&res,&nf,0);
                    break;
    }
    if (r < 0) {
        printf1(ctx, "Error in function evaluation.\n");
        goto TIFin;
    }
    printf1(ctx, "Approximation: "); 
    rt_printf1_d(ctx, ctx->PMFmtS,res); 

    if (r != 0)
        printf1(ctx, " (not successful)");

    printf1(ctx, "\nNumber of function calls: %d\n",nf);
    err = 0;

TIFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_int(m,rerr)     print method etc for numerical integration.         */
/*                      m = method, rerr = relative error.                  */
/*                      return 0 if OK, -1 if error.                        */

int prn_int(TDAContext *ctx, int m,double rerr)
{

    printf1(ctx, "Method "); 
    switch (m) {
        case 1:     printf1(ctx, "1 (QNG).\n"); 
                    break;
        case 2:     printf1(ctx, "2 (QSUB).\n");
                    break;
        case 3:     printf1(ctx, "3 (SNIFF).\n");
                    break;
        case 4:     printf1(ctx, "4 (QTRAP).\n");
                    break;
        case 5:     printf1(ctx, "5 (QSIMP).\n");
                    break;
        default:    printf1(ctx, "%d (not available).\n",m); 
                    return(-1);
    }
    printf1(ctx, "Relative error: %12.5e\n",rerr);          
    return(0);
}

/****************************************************************************/
/*  ni_gen(a,b,ftyp,nhp,htyp,err)                                           */
/*                                                                          */
/*      Numerical integration with method NINTMETH. Integration interval is */
/*      [a,b]. Function depends on ftyp.                                    */
/*                                                                          */
/*      ftyp = 1 :  user-defined function (get_integrand)                   */  
/*             2 :  MPOL1                                                   */
/*             3 :  MPOL1                                                   */
/*             4 :  MPOL1                                                   */
/*                                                                          */
/*      If nhp > 0 Hermite integration, ignore a, b, and NINTMETH.          */
/*      If htyp = 0 standard Hermite integration, if htyp = 1 integration   */
/*      with normal density.                                                */
/*                                                                          */
/*      err = 0 if successful                                               */
/*            1 if not successful                                           */
/*           -1 if error in function evaluation.                            */
/*           -2 if error in Hermite integration                             */
/*                                                                          */
/*      Return: value of integral                                           */

double ni_gen(TDAContext *ctx, double a,double b,int ftyp,int nhp,int htyp,int *err)
{
    int r = 0,nf = 0,nfu = 0,rtyp = 0,n = 0,i = 0,ii = 0;
    double res = 0.0,aerr = 0.0,x = 0.0,tmp = 0.0;
    double tp = 0.5641895835477563; /* 1 / sqrt(pi) */
    double s2 = 1.4142135623730950; /* sqrt(2) */

    if (nhp > 0) {      /* ## Hermite integration */

        if (nhp < 2 || nhp > 12 || nhp == 8 || nhp == 10 || nhp == 11) {
            *err = -2;
            return(0.0);
        }
        res = 0.0;
        if (ctx->HIW[nhp][0] > 0.0) {
            tmp = ni_func(ctx, 0.0,err,ftyp);
/*    printf1("x=%lg nhp=%d tmp=%lg err=%d \n",0.0,nhp,tmp,*err);      
*/
            if (*err)
                return(0.0);
            if (htyp == 1)
                res = ctx->HIW[nhp][0] * tmp;
            else
                res = ctx->HIW1[nhp][0] * tmp;                     
        }  
        n = nhp / 2;

        for (ii = 1; ii <= 2; ++ii) {
            for (i = 1; i <= n; ++i) {

                if (ii == 1)
                    x = ctx->HIZ[nhp][i];
                else
                    x = -ctx->HIZ[nhp][i];
/**
                tmp = ni_func(ctx, x,err,ftyp);
printf1(ctx, "x=%lg nhp=%d tmp=%lg err=%d ftyp=%d ii=%d i=%d \n",x,nhp,tmp,*err,ftyp,ii,i);
                if (*err)
                    return(0.0);
**/


                if (htyp == 1)  {
                    res += ctx->HIW[nhp][i] * ni_func(ctx, x * s2,err,ftyp);
/*
                    tmp = ni_func(ctx, x * s2,err,ftyp);
printf1(ctx, "x=%lg nhp=%d tmp=%lg err=%d ftyp=%d ii=%d i=%d \n",x,nhp,tmp,*err,ftyp,ii,i);
*/
}

                else    {
                    res += ctx->HIW1[nhp][i] * ni_func(ctx, x,err,ftyp);
/*
                    tmp  =  ni_func(ctx, x,err,ftyp);
printf1(ctx, "x=%lg nhp=%d tmp=%lg err=%d ftyp=%d ii=%d i=%d \n",x,nhp,tmp,*err,ftyp,ii,i);
*/
}
                if (*err)
                    return(0.0);
            }
        }
        if (htyp == 1)
            res *= tp;
        return(res);
    }
         
    switch (ctx->NINTMETH) {
        case 1:     r = qng(ctx, a,b,ctx->NINTAERR,ctx->NINTRERR,&res,&aerr,&nf,ftyp);
                    break;
        case 2:     r = qsub(ctx, a,b,ctx->NINTRERR,&res,&aerr,&nf,&rtyp,ftyp);
                    break;
        case 3:     r = qsniff(ctx, a,b,ctx->NINTRERR,&res,&nf,&nfu,ftyp);
                    break;
        case 4:     r = qtrap(ctx, a,b,ctx->NINTRERR,&res,&nf,ftyp);
                    break;
        case 5:     r = qsimp(ctx, a,b,ctx->NINTRERR,&res,&nf,ftyp);
                    break;
    }
    if (r)
        res = 0.0;
    *err = r;
    return(res);
}

/* ------------------------------------------------------------------------ */
/*  ni_func(x,err,ftyp)  function used for numerical integration.           */

double ni_func(TDAContext *ctx, double x,int *err,int ftyp)
{
    double f = 0.0;
        
    switch (ftyp) {
        case  0:
            *err = get_flval(ctx, &f,1,&x,0,0,&f,&f,&f);
            break;
             
        case  1:    /* special for user-defined functions */

            f = get_integrand(ctx, x,err);
            break;      

        case  2:
            f = get_mpol1(ctx, x,err);
            break;
        case  3:
            f = get_mpol2(ctx, x,err);
            break;
        case  4:
            f = get_mpol3(ctx, x,err);
            break;
    }
    return(f);
}

/****************************************************************************/
/*  qsniff                                                                  */
/*                                                                          */
/*      Numerical integration, adapted from:rischen Integration             */
/*                                                                          */
/*      Ref. S. Garribba, L. Quartapelle, G. Reina, SNIFF: Efficient        */
/*      Self-Tuning Algorithm for Numerical Integration, in: Computing      */
/*      20 (1978), 363 - 375.                                               */
/*                                                                          */
/*      int qsniff(a,b,rerr,val,nf,nfu,ftyp)                                */
/*                                                                          */
/*      double a,b      input:  integration interval.                       */
/*      double rerr     input:  relative error.                             */
/*      double *val     output: value of integral.                          */
/*      int *nf         output: number of function calls                    */
/*      int *nfu        output: number of finally used function values      */
/*                                                                          */
/*      if ftyp = 0  use standard function,                                 */
/*                1  use user-defined rate models.                          */
/*               >1  other functions                                        */
/*                                                                          */
/*      Return:         0  if successful (or mostly successful)             */
/*                      1  not successful                                   */
/*                     -1  error in function evaluation.                    */

int qsniff(TDAContext *ctx, double a,double b,double rerr,double *val,int *nf,int *nfu,int ftyp)
{
    int i,j,n2,n3,mar,err;
    double ac,ap,ba,bc,c,ca,cp,ds,e,ea,h,r,s,sa,sc,sp,ss,sss,t,tmp,tmp1;

    int no = 6;      /* number of points of the elementary integrator */

    *nf = *nfu = 0;

    t  = 1.2;
    sp = 4.0;

    tmp = log10(rerr);
    tmp1 = -sqrt(tmp * tmp + sp * sp);
    e = pow(10.0,tmp1) / (2.0 * pow(t,12.0));

    n2 = 2 * no;
    n3 = 3 * no;

    c = b;
    sss = 0.0;
    ds  = 1.0;
    mar = i = j = 0;
  
SNIter:
    while (1) {
        if (a == b) {
            *val = sss;     /* estimated value of integral */
            *nf = i;        /* number of function calls */
            *nfu = j;       /* finally used function values */
            return(mar);
        }
        if (a != c)
            break;

        /**  mar = 2; **/
        mar = 1;
        a = b;
    }
    ac = (a + c) / 2.0;
    sa = qsniff1(ctx, a,ac,ftyp,&err);
    if (err)
        return(-1);

    sc = qsniff1(ctx, ac,c,ftyp,&err);
    if (err)
        return(-1);

    s  = qsniff1(ctx, a,c,ftyp,&err);
    if (err)
        return(-1);

    ss = sa + sc;
    i += n3;

    if (c == b && fabs(sa) < fabs(sc)) {
        ds = -ds;
        s = -s;
        ss = -ss;
        c = a;
        a = b;
        b = c;
    }
    ba = b - a;
    bc = b - c;
    ca = c - a;

    if (ss == s) {
        sss += ds * ss;
        j += n2;
        e *= fabs(bc / ca);
        a = c;
        c = b;
        goto SNIter;
    }
    h = 1.0;
    if (ss != 0.0) {
        if ((h = fabs(sss / ss)) < 1.0)
            h = 1.0;
    }
    ea = h * e;
    if (ea < ctx->EPSI) {
        ea = ctx->EPSI;
        /** mar = 1; **/
        mar = 0;
    }

    tmp = ((ss - s) / (4096.0 * ss - s)) / ea;
    r = pow(fabs(tmp),1.0 / 12.0);

    if (r > t) {
        e /= r;
        c = a + ca / r;
        goto SNIter;
    }
    sss += ds * ss;
    j += n2;
    ap = c;
    cp = c + ca / r;

    if ((a < c && cp >= b) || (a >= c && cp < b)) {
        e = (e / r) * fabs(ba / (cp - a));
        a = ap;
        c = b;
    }
    else {
        e /= r;
        a = ap;
        c = cp;
    }
    goto SNIter;
}

/* ------------------------------------------------------------------------ */
/*  qsniff1(u,v,ftyp,err)                                                   */
/*                                                                          */
/*                  Elementary integrator: 6-point Gauss Legendre Rule.     */
/*                  Used by qsniff().                                       */

double qsniff1(TDAContext *ctx, double u,double v,int ftyp,int *err)
{
    int k = 0;
    double du = 0.0,uo = 0.0,p = 0.0,uk = 0.0,f1 = 0.0,f2 = 0.0;

    static double w[3] = { 0.171324492379170,
                           0.360761573048139,
                           0.467913934572691};
    static double x[3] = { 0.932469514203152,
                           0.661209386466265,
                           0.238619186083196};
    du = (v - u) / 2.0;
    uo = (u + v) / 2.0;
    p  = 0.0;

    for (k = 0; k < 3; ++k) {
        uk = x[k] * du;
        f1 = ni_func(ctx, uo + uk,err,ftyp);
        if (*err)
            return(0.0);

        f2 = ni_func(ctx, uo - uk,err,ftyp);
        if (*err)
            return(0.0);

        p += w[k] * (f1 + f2);           
    }
    p *= du;
    return(p);
}

/****************************************************************************/
/*  qng                                                                     */
/*      Numerical integration:                                              */
/*                                                                          */
/*      Calculates: I = Integral of F over (a,b)                            */
/*      hopefully with                                                      */
/*                                                                          */
/*          ABS(I - Result) <= MAX(epsabs,epsrel * ABS(I))                  */
/*                                                                          */
/*      Algorithm adopted from:                                             */
/*      R. Piessens, E. de Doncker-Kapenga, C.W. Ueberhuber, D.K.           */
/*      Kahaner, Quadpack. A Subroutine Package for Automatic Integration,  */
/*      Springer 1983, p. 130 - 136.                                        */
/*                                                                          */
/*      int qng(a,b,epsabs,epsrel,res,abserr,nf,ftyp)                       */
/*                                                                          */
/*      double a,b      input:  integration interval                        */
/*      double epsabs   input:  absolute error                              */
/*      double epsrel   input:  relative error                              */
/*      double *res     output: approximation to integral                   */
/*      double *abserr  output: estimate of resulting abs error             */
/*      int *nf         output: number of function calls                    */
/*                                                                          */
/*      if ftyp = 0  use standard function,                                 */
/*                1  use user-defined rate models.                          */
/*               >1  other functions                                        */
/*                                                                          */
/*      Return:         0  successful                                       */
/*                      1  not successful                                   */
/*                     -1  if error in function evaluation                  */
/*                     -2  error in input parameter                         */

int qng(TDAContext *ctx, double a,double b,double epsabs,double epsrel,double *res, double *abserr,int *nf,int ftyp)
{
    int l,k,ipx,err;
    double tmp,res10,res21,res43 = 0.0,res87,hlgth,dhlgth,centr,fcentr;
    double absc,fval1,fval2,fval,reskh,resasc = 0.0,resabs = 0.0,f1,f2;
    double fv1[5],fv2[5],fv3[5],fv4[5],savfun[21];

    static double x1[5] = {    9.739065285171717e-1, 8.650633666889845e-1,
                               6.794095682990244e-1, 4.333953941292472e-1,
                               1.488743389816312e-1 };

    static double x2[5] = {    9.956571630258081e-1, 9.301574913557082e-1,
                               7.808177265864169e-1, 5.627571346686047e-1,
                               2.943928627014602e-1 };

    static double x3[11] = {   9.993333609019321e-1, 9.874334029080889e-1,
                               9.548079348142663e-1, 9.001486957483283e-1,
                               8.251983149831142e-1, 7.321483889893050e-1,
                               6.228479705377252e-1, 4.994795740710565e-1,
                               3.649016613465808e-1, 2.222549197766013e-1,
                               7.465061746138332e-2 };

    static double x4[22] = {   9.999029772627292e-1, 9.979898959866787e-1, 
                               9.921754978606872e-1, 9.813581635727128e-1, 
                               9.650576238583846e-1, 9.431676131336706e-1,
                               9.158064146855072e-1, 8.832216577713165e-1,
                               8.457107484624157e-1, 8.035576580352310e-1,
                               7.570057306854956e-1, 7.062732097873218e-1,
                               6.515894665011779e-1, 5.932233740579611e-1,
                               5.314936059708319e-1, 4.667636230420228e-1,
                               3.994248478592188e-1, 3.298748771061883e-1,
                               2.585035592021616e-1, 1.856953965683467e-1,
                               1.118422131799075e-1, 3.735212339461987e-2 };

    static double w10[5] = {   6.667134430868814e-2, 1.494513491505806e-1,
                               2.190863625159820e-1, 2.692667193099964e-1,
                               2.955242247147529e-1 };

    static double w21a[5] = {  3.255816230796473e-2, 7.503967481091995e-2,
                               1.093871588022976e-1, 1.347092173114733e-1,
                               1.477391049013385e-1 };

    static double w21b[6] = {  1.169463886737187e-2, 5.475589657435200e-2,
                               9.312545458369761e-2, 1.234919762620659e-1,
                               1.427759385770601e-1, 1.494455540029169e-1 };

    static double w43a[10] = { 1.629673428966656e-2, 3.752287612086950e-2,
                               5.469490205825544e-2, 6.735541460947809e-2,
                               7.387019963239395e-2, 5.768556059769796e-3, 
                               2.737189059324884e-2, 4.656082691042883e-2,
                               6.174499520144256e-2, 7.138726726869340e-2 };

    static double w43b[12] = { 1.844477640212414e-3, 1.079868958589165e-2,
                               2.189536386779543e-2, 3.259746397534569e-2,
                               4.216313793519181e-2, 5.074193960018458e-2,
                               5.837939554261925e-2, 6.474640495144589e-2,
                               6.956619791235648e-2, 7.282444147183321e-2,
                               7.450775101417512e-2, 7.472214751740301e-2 };

    static double w87a[21] = { 8.148377384149173e-3, 1.876143820156282e-2,
                               2.734745105005229e-2, 3.367770731163793e-2,
                               3.693509982042791e-2, 2.884872430211531e-3,
                               1.368594602271270e-2, 2.328041350288831e-2,
                               3.087249761171336e-2, 3.569363363941877e-2,
                               9.152833452022414e-4, 5.399280219300471e-3,
                               1.094767960111893e-2, 1.629873169678734e-2,
                               2.108156888920384e-2, 2.537096976925383e-2,
                               2.918969775647575e-2, 3.237320246720279e-2,
                               3.478309895036514e-2, 3.641222073135179e-2,
                               3.725387550304771e-2 };

    static double w87b[23] = { 2.741455637620724e-4, 1.807124155057943e-3, 
                               4.096869282759165e-3, 6.758290051847379e-3, 
                               9.549957672201647e-3, 1.232944765224485e-2,
                               1.501044734638895e-2, 1.754896798624319e-2, 
                               1.993803778644089e-2, 2.219493596101229e-2,
                               2.433914712600081e-2, 2.637450541483921e-2, 
                               2.828691078877120e-2, 3.005258112809270e-2,
                               3.164675137143993e-2, 3.305041341997850e-2,
                               3.425509970422606e-2, 3.526241266015668e-2,
                               3.607698962288870e-2, 3.669860449845609e-2,
                               3.712054926983258e-2, 3.733422875193504e-2,
                               3.736107376267902e-2 };

    *res = *abserr = 0.0;
    *nf = 0;

    if (epsabs < 0.0 && epsrel < 0.0)       /* check parameter */
        return(-2);

    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);
    centr = 0.5 * (b + a);
    fcentr = ni_func(ctx, centr,&err,ftyp);
    if (err)
        return(-1);
    *nf = 21;

    for (l = 1; l <= 3; ++l) {

        if (l == 1) {
            res10 = 0.0;
            res21 = w21b[5] * fcentr;
            resabs = w21b[5] * fabs(fcentr);

            for (k = 0; k < 5; ++k) {
                absc = hlgth * x1[k];
                fval1 = ni_func(ctx, centr + absc,&err,ftyp);
                if (err)
                    return(-1);
                fval2 = ni_func(ctx, centr - absc,&err,ftyp);
                if (err)
                    return(-1);

                fval = fval1 + fval2;
                res10 += w10[k] * fval;
                res21 += w21a[k] * fval;
                resabs += w21a[k] * (fabs(fval1) + fabs(fval2));
                savfun[k] = fval;
                fv1[k] = fval1;
                fv2[k] = fval2;
            }
            ipx = 5;
            for (k = 0; k < 5; ++k) {
                absc = hlgth * x2[k];
                fval1 = ni_func(ctx, centr + absc,&err,ftyp);
                if (err)
                    return(-1);
                fval2 = ni_func(ctx, centr - absc,&err,ftyp);
                if (err)
                    return(-1);

                fval  = fval1 + fval2;
                res21 += w21b[k] * fval;
                resabs += w21b[k] * (fabs(fval1) + fabs(fval2));
                savfun[ipx++] = fval;
                fv3[k] = fval1;
                fv4[k] = fval2;
            }
            *res = res21 * hlgth;
            resabs *= dhlgth;
            reskh = 0.5 * res21;
            resasc = w21b[5] * fabs(fcentr - reskh);
            for (k = 0; k < 5; ++k) {
              resasc += w21a[k] * (fabs(fv1[k] - reskh) + fabs(fv2[k] - reskh))
                      + w21b[k] * (fabs(fv3[k] - reskh) + fabs(fv4[k] - reskh));
            }
            *abserr = fabs((res21 - res10) * hlgth);
            resasc *= dhlgth;
        }
        else if (l == 2) {

            res43 = w43b[11] * fcentr;
            *nf = 43;
            for (k = 0; k < 10; ++k)  
                res43 += savfun[k] * w43a[k];

            for (k = 0; k < 11; ++k) {
                absc = hlgth * x3[k];

                f1 = ni_func(ctx, absc + centr,&err,ftyp);
                if (err)
                    return(-1);
                f2 = ni_func(ctx, centr - absc,&err,ftyp);
                if (err)
                    return(-1);

                fval = f1 + f2;
                res43 += fval * w43b[k];
                savfun[ipx++] = fval;
            }
            *res = res43 * hlgth;
            *abserr = fabs((res43 - res21) * hlgth);
        }
        else {
            res87 = w87b[22] * fcentr;
            *nf = 87;
            for (k = 0; k < 21; ++k)
                res87 += savfun[k] * w87a[k];

            for (k = 0; k < 22; ++k) {
                absc = hlgth * x4[k];

                f1 = ni_func(ctx, absc + centr,&err,ftyp);
                if (err)
                    return(-1);
                f2 = ni_func(ctx, centr - absc,&err,ftyp);
                if (err)
                    return(-1);

                res87 += w87b[k] * (f1 + f2);
            }
            *res = res87 * hlgth;
            *abserr = fabs((res87 - res43) * hlgth);
        }
        if (resasc != 0.0 && *abserr != 0.0) {
            tmp = pow(200.0  * *abserr / resasc,1.5);
            if (tmp > 1.0)
                tmp = 1.0;

            *abserr = resasc * tmp;
        }
        tmp = 50.0 * ctx->EPSI;
        if (resabs > ctx->DBLMIN / tmp) {
            *abserr = dmax(ctx, *abserr,tmp * resabs);
        }
        tmp = dmax(ctx, epsrel * fabs(*res),epsabs);
        if (*abserr <= tmp)
            return(0);
    }
    return(1);  /* no convergence */
}

/****************************************************************************/
/*  qsub                                                                    */
/*                                                                          */
/*      Numerical integreation, adapted from:                               */
/*      Ref. T.N.L. Patterson, Algorithm for Automatic Numerical            */
/*      Integration Over a Finite Intervall (Algorithm 468), Comm. of the   */
/*      ACM 16 (1973), 694 - 699.                                           */
/*                                                                          */
/*      int qsub(a,b,rerr,val,relerr,nfs,rtyp,ftyp)                         */
/*                                                                          */
/*      double a,b      input:  integration interval                        */
/*      double rerr     input:  relative error                              */
/*      double *val     output: approximation to integral                   */
/*      double *relerr  output: estimated rel error.                        */
/*      int *nfs        output: number of function calls                    */
/*      int *rtyp       output: 0  successful (without subdivision)         */
/*                              1  successful (with subdivision)            */
/*                              2  successful (but relaxed criterion)       */
/*                            < 0  probably not reliable                    */
/*                                                                          */
/*      if ftyp = 0  use standard function,                                 */
/*                1  use user-defined rate models.                          */
/*               >1  other functions                                        */
/*                                                                          */
/*      Return:         0  successful                                       */
/*                      1  not successful                                   */
/*                     -1  error in function evaluation.                    */

int qsub(TDAContext *ctx, double a,double b,double rerr,double *val,double *relerr,int *nfs, int *rtyp,int ftyp)
{
    int ic,j,jj,m1,m2,n,nf,bad,out,r,rhs;
    double alpha,beta,ival,ival1,estim,comp,comp1,h;
    int nmax = 4096;

    *nfs = 0;
    *rtyp = 0;

    r = qsubf(ctx, a,b,rerr,&ival,&ival1,&nf,ftyp);   

    *nfs += nf;

    if (r < 0)
        return(-1);

    if (r == 0) {
        if (ival != 0.0)
            *relerr = fabs((ival - ival1) / ival);
        else
            *relerr = 0.0;

        *val = ival;
        return(0);
    }
    estim = fabs(ival * rerr);
    ic = 1;
    rhs = 0;
    n = 1;
    h = b - a;
    bad = 1;

Lab10:
    *val = 0.0;
    *relerr = 0.0;
    h *= 0.5;
    n += n;
   
    m1 = bad;
    m2 = bad + 1;
    out = 1;

Lab50:
    if (m1 > m2)
        goto Lab90;

    comp1 = 0.0;

    for (jj = m1; jj <= m2; ++jj) {

        if (rhs)
            j = m2 + m1 - jj;
        else
            j = jj;

        alpha = a + h * ((double)j - 1.0);
        beta = alpha + h;

        r = qsubf(ctx, alpha,beta,rerr,&ival,&ival1,&nf,ftyp);
        *nfs += nf;

        if (r < 0)
            return(-1);

        comp = fabs(ival - ival1);
        comp1 += comp;

        if (r) {
            if (comp <= estim) {
                if (ic >= 0)
                    ic = 2;
                else
                    ic = -2;
            }
            else if (n < nmax) {
                bad = 2 * j - 1;
                if (((j / 2) * 2) == j)
                    rhs = 1;
                else
                    rhs = 0;

                goto Lab10;
            }
            else if (ic > 0)   
                ic = -ic;
        }
        *val += ival;
    }
    *relerr += comp1;

Lab90:
    if (out == 1) {
        m1 = 1;
        m2 = bad - 1;
        rhs = 0;
        out = 2;
    }
    else if (out == 2) {
        m1 = bad + 2;
        m2 = n;
        out = 3;
    }
    else {
        *relerr /= fabs(*val);
        *rtyp = ic;
        if (ic >= 0)
            return(0);
        else
            return(1);
    }
    goto Lab50;
}

/****************************************************************************/
/*  qsubf                                                                   */
/*                                                                          */
/*      Called by qsub(). Adapted from                                      */
/*      Ref. T.N.L. Patterson, Algorithm for Automatic Numerical            */
/*      Integration Over a Finite Intervall (Algorithm 468), Comm. of the   */
/*      ACM 16 (1973), 694 - 699.                                           */
/*                                                                          */
/*      int qsubf(a,b,rerr,val,val1,nf,ftyp)                                */
/*                                                                          */
/*      double a,b          input:  integration interval                    */
/*      double rerr         input:  relative error                          */
/*      double *val,*val1   output: value of integral in last two           */
/*                                  approximations                          */
/*      int *nf             output: number of function calls                */
/*                                                                          */
/*      if ftyp = 0  use standard function,                                 */
/*                1  use user-defined rate models.                          */
/*               >1  other functions                                        */
/*                                                                          */
/*      Return:         0  if successful                                    */
/*                      1  if not successful                                */
/*                     -1  if error in function evaluation                  */

int qsubf(TDAContext *ctx, double a,double b,double rerr,double *val,double *val1,int *nf,int ftyp)
{
    int i,j,k,iold,inew,err;
    double f1,f2,sum,diff,fzero,acum,x,funct[130];

    static double p[382] = {                   0.0000000000000000E00,
                        0.7745966692414833E00, 0.5555555555555555E00,
                        0.8888888888888888E00, 0.2684880898683334E00,
                        0.9604912687080202E00, 0.1046562260264672E00,
                        0.4342437493468025E00, 0.4013974147759622E00,
                        0.4509165386584741E00, 0.1344152552437842E00,
                        0.5160328299707973E-1, 0.2006285293769890E00,
                        0.9938319632127550E00, 0.1700171962994026E-1,
                        0.8884592328722569E00, 0.9292719531512453E-1,
                        0.6211029467372264E00, 0.1715119091363913E00,
                        0.2233866864289668E00, 0.2191568584015874E00,
                        0.2255104997982066E00, 0.6720775429599070E-1,
                        0.2580759809617665E-1, 0.1003142786117955E00,
                        0.8434565739321106E-2, 0.4646289326175798E-1,
                        0.8575592004999035E-1, 0.1095784210559246E00,
                          
                        0.9990981249676675E00, 0.2544780791561874E-2,
                        0.9815311495537401E00, 0.1644604985438781E-1,
                        0.9296548574297400E00, 0.3595710330712932E-1,
                        0.8367259381688687E00, 0.5697950949412335E-1,
                        0.7024962064915270E00, 0.7687962049900353E-1,
                        0.5313197436443756E00, 0.9362710998126447E-1,
                        0.3311353932579768E00, 0.1056698935802348E00,
                        0.1124889431331866E00, 0.1119568730209534E00,
                        0.1127552567207686E00, 0.3360387714820773E-1,
                        0.1290380010035126E-1, 0.5015713930589953E-1,
                        0.4217630441558854E-2, 0.2323144663991026E-1,
                        0.4287796002500773E-1, 0.5478921052796286E-1,
                        0.1265156556230068E-2, 0.8223007957235929E-2,
                        0.1797855156812827E-1, 0.2848975474583354E-1,

                        0.3843981024945553E-1, 0.4681355499062801E-1,
                        0.5283494679011651E-1, 0.5597843651047631E-1,
                        0.9998728881203576E00, 0.3632214818455306E-3,
                        0.9972062593722219E00, 0.2579049794685688E-2,
                        0.9886847575474294E00, 0.6115506822117246E-2,
                        0.9721828747485817E00, 0.1049824690962132E-1,
                        0.9463428583734029E00, 0.1540675046655949E-1,
                        0.9103711569570042E00, 0.2059423391591271E-1,
                        0.8639079381936904E00, 0.2586967932721474E-1,
                        0.8069405319502176E00, 0.3107355111168796E-1,
                        0.7397560443526947E00, 0.3606443278078257E-1,
                        0.6629096600247805E00, 0.4071551011694431E-1,
                        0.5771957100520458E00, 0.4491453165363219E-1,
                        0.4836180269458410E00, 0.4856433040667319E-1,

                        0.3833593241987303E00, 0.5158325395204845E-1,
                        0.2777498220218243E00, 0.5390549933526606E-1,
                        0.1632352515522074E00, 0.5548140435655936E-1,
                        0.5634431304659278E-1, 0.5627769983125430E-1,
                        0.5637762836038471E-1, 0.1680193857410386E-1,
                        0.6451900050175736E-2, 0.2507856965294976E-1,
                        0.2108815245726632E-2, 0.1161572331995513E-1,
                        0.2143898001250386E-1, 0.2739460526398143E-1,
                        0.6326073193626335E-3, 0.4111503978654693E-2,
                        0.8989275784064135E-2, 0.1424487737291677E-1,
                        0.1921990512472776E-1, 0.2340677749531400E-1,
                        0.2641747339505825E-1, 0.2798921825523815E-1,
                        0.1807395644453883E-3, 0.1289524082610417E-2,
                        0.3057753410175531E-2, 0.5249123454808859E-2,

                        0.7703375233279741E-2, 0.1029711695795635E-1,
                        0.1293483966360737E-1, 0.1553677555584398E-1,
                        0.1803221639039128E-1, 0.2035775505847215E-1,
                        0.2245726582681609E-1, 0.2428216520333659E-1,
                        0.2579162697602422E-1, 0.2695274966763303E-1,
                        0.2774070217827968E-1, 0.2813884991562715E-1,
                        0.9999824303548915E00, 0.5053609520786251E-4,
                        0.9995987996719106E00, 0.3777466463269846E-3,
                        0.9983166353184073E00, 0.9383698485423815E-3,
                        0.9957241046984071E00, 0.1681142865421469E-2,
                        0.9914957211781061E00, 0.2568764943794020E-2,
                        0.9853714995985203E00, 0.3572892783517299E-2,
                        0.9771415146397057E00, 0.4671050372114321E-2,
                        0.9666378515584165E00, 0.5843449875835639E-2,

                        0.9537300064257611E00, 0.7072489995433555E-2,
                        0.9383203977795928E00, 0.8342838753968157E-2,
                        0.9203400254700124E00, 0.9641177729702536E-2,
                        0.8997448997769400E00, 0.1095573338783790E-1,
                        0.8765134144847052E00, 0.1227583056008277E-1,
                        0.8506444947683502E00, 0.1359157100976554E-1,
                        0.8221562543649804E00, 0.1489364166481518E-1,
                        0.7910849337998483E00, 0.1617321872957771E-1,
                        0.7574839663805136E00, 0.1742193015946417E-1,
                        0.7214230853700989E00, 0.1863184825413879E-1,
                        0.6829874310910792E00, 0.1979549504809749E-1,
                        0.6422766425097595E00, 0.2090585144581202E-1,
                        0.5994039302422428E00, 0.2195636630531782E-1,
                        0.5544951326319325E00, 0.2294096422938774E-1,

                        0.5076877575337166E00, 0.2385405210603854E-1,
                        0.4591300119898323E00, 0.2469052474448767E-1,
                        0.4089798212298886E00, 0.2544576996546476E-1,
                        0.3574038378315321E00, 0.2611567337670609E-1,
                        0.3045764415567140E00, 0.2669662292745035E-1,
                        0.2506787303034831E00, 0.2718551322962479E-1,
                        0.1958975027111001E00, 0.2757974956648187E-1,
                        0.1404242331525601E00, 0.2787725147661370E-1,
                        0.8445404008371088E-1, 0.2807645579381724E-1,
                        0.2818464894974569E-1, 0.2817631903301660E-1,
                        0.2818881418019235E-1, 0.8400969287051932E-2,
                        0.3225950025087868E-2, 0.1253928482647488E-1,
                        0.1054407622863316E-2, 0.5807861659977567E-2,
                        0.1071949000625193E-1, 0.1369730263199071E-1,

                        0.3163036608222644E-3, 0.2055751989327346E-2,
                        0.4494637892032067E-2, 0.7122438686458387E-2,
                        0.9609952562363883E-2, 0.1170338874765700E-1,
                        0.1320873669752912E-1, 0.1399460912761907E-1,
                        0.9037273465875114E-4, 0.6447620413057247E-3,
                        0.1528876705087765E-2, 0.2624561727404429E-2,
                        0.3851687616639870E-2, 0.5148558478978177E-2,
                        0.6467419831803686E-2, 0.7768387777921991E-2,
                        0.9016108195195643E-2, 0.1017887752923607E-1,
                        0.1122863291340804E-1, 0.1214108260166829E-1,
                        0.1289581348801211E-1, 0.1347637483381651E-1,
                        0.1387035108913984E-1, 0.1406942495781357E-1,
                        0.2515787038428066E-4, 0.1888732645065049E-3,
                        0.4691849242478504E-3, 0.8405714327107224E-3,

                        0.1284382471897010E-2, 0.1786446391758649E-2,
                        0.2335525186057160E-2, 0.2921724937917819E-2,
                        0.3536244997716777E-2, 0.4171419376984078E-2,
                        0.4820588864851268E-2, 0.5477866693918950E-2,
                        0.6137915280041385E-2, 0.6795785504882773E-2,
                        0.7446820832407591E-2, 0.8086609364788859E-2,
                        0.8710965079732086E-2, 0.9315924128069395E-2,
                        0.9897747524048749E-2, 0.1045292575290601E-1,
                        0.1097818315265891E-1, 0.1147048211469387E-1,
                        0.1192702605301927E-1, 0.1234526237224383E-1,
                        0.1272288498273238E-1, 0.1305783668835304E-1,
                        0.1334831146372517E-1, 0.1359275661481239E-1,
                        0.1378987478324093E-1, 0.1393862573830685E-1,
                        0.1403822789690862E-1, 0.1408815951650830E-1,

                        0.9999975963797484E00, 0.6937936432410826E-5,
                        0.9999439962070543E00, 0.5327529366978061E-4,
                        0.9997604909244320E00, 0.1357549109492287E-3,
                        0.9993803380250235E00, 0.2492124004829972E-3,
                        0.9987456144680951E00, 0.3897452844732822E-3,
                        0.9978053544959572E00, 0.5542953149303747E-3,
                        0.9965141459148902E00, 0.7402828042445033E-3,
                        0.9948315028006210E00, 0.9453615168585253E-3,
                        0.9927213442827886E00, 0.1167484117429959E-2,
                        0.9901513704007701E00, 0.1404907995655144E-2,
                        0.9870925279540340E00, 0.1656112728154452E-2,
                        0.9835186975786327E00, 0.1919712971013872E-2,
                        0.9794062816708626E00, 0.2194406925363838E-2,
                        0.9747344597524026E00, 0.2478958226657567E-2,

                        0.9694846595024592E00, 0.2772195764593450E-2,
                        0.9636406215698121E00, 0.3073018434702578E-2,
                        0.9571882161098609E00, 0.3380397991086920E-2,
                        0.9501152975212948E00, 0.3693377917025650E-2,
                        0.9424115651910830E00, 0.4011068724075023E-2,
                        0.9340684361577257E00, 0.4332640968092982E-2,
                        0.9250789329070756E00, 0.4657317299756854E-2,
                        0.9154375871557650E00, 0.4984364564765538E-2,
                        0.9051403588132615E00, 0.5313086605187056E-2,
                        0.8941845683355590E00, 0.5642818101384444E-2,
                        0.8825688402473419E00, 0.5972919565508165E-2,
                        0.8702930555481139E00, 0.6302773449085758E-2,
                        0.8573583108862321E00, 0.6631781242901887E-2,
                        0.8437668826727086E00, 0.6959361409390422E-2,

                        0.8295221946374014E00, 0.7284947980553807E-2,
                        0.8146287876551374E00, 0.7607989665719056E-2,
                        0.7990922909608414E00, 0.7927949334294849E-2,
                        0.7829193941182830E00, 0.8244303763032868E-2,
                        0.7661178193037600E00, 0.8556543561307689E-2,
                        0.7486962936169366E00, 0.8864173209482494E-2,
                        0.7306645212421812E00, 0.9166711163560788E-2,
                        0.7120331553622520E00, 0.9463689993830065E-2,
                        0.6962137697791147E00, 0.9754656536317411E-2,
                        0.6730188302304184E00, 0.1003917204405684E-1,
                        0.6526616654100174E00, 0.1031681233094762E-1,
                        0.6317564377111942E00, 0.1058716790488519E-1,
                        0.6103181137151864E00, 0.1084984408933731E-1,
                        0.5883624344476625E00, 0.1110446113400692E-1,

                        0.5659058854236544E00, 0.1135065431598059E-1,
                        0.5429656664983114E00, 0.1158807403304395E-1,
                        0.5195596615374570E00, 0.1181638589083023E-1,
                        0.4957064079187614E00, 0.1203527078527956E-1,
                        0.4714250658716588E00, 0.1224442498161198E-1,
                        0.4467353876620284E00, 0.1244356019071403E-1,
                        0.4216576866261633E00, 0.1263240364354207E-1,
                        0.3962128060576159E00, 0.1281069816387736E-1,
                        0.3704220879500782E00, 0.1297820223953739E-1,
                        0.3443073415994380E00, 0.1313469009196015E-1,
                        0.3178908120684766E00, 0.1327995174393053E-1,
                        0.2911951485182466E00, 0.1341379308511009E-1,
                        0.2642433724109267E00, 0.1353603593495621E-1,
                        0.2370588455898297E00, 0.1364651810257129E-1,

                        0.2096652382431811E00, 0.1374509344300189E-1,
                        0.1820864967592521E00, 0.1383163190950642E-1,
                        0.1543468114813781E00, 0.1390601960132546E-1,
                        0.1264705843723019E00, 0.1396815880651693E-1,
                        0.9848239659811920E-1, 0.1401796803945660E-1,
                        0.7040697604285517E-1, 0.1405538207264996E-1,
                        0.4226916476536360E-1, 0.1408035196255366E-1,
                        0.1409388641078246E-1, 0.1409284506916040E-1,
                        0.1409440709009617E-1 };

    *nf = 0;
    *val = *val1 = 0.0;
  
    if (rerr <= 0.0 || a >= b)  
        return(0);

    sum   = (b + a) / 2.0;
    diff  = (b - a) / 2.0;
    fzero = ni_func(ctx, sum,&err,ftyp);
    if (err)
        return(-1);

    *val  = 2.0 * fzero * diff;

    i = iold = 0; inew = 1; acum = 0.0; 

    for (k = 2; k <= 8; ++k) {

        *val1 = *val;

        if (k > 2) {
            acum = 0.0;
            for (j = 1; j <= iold; ++j)  
                acum += p[++i] * funct[j];
        }
        iold += inew;
        for (j = inew; j <= iold; ++j) {
            x = p[++i] * diff;
            f1 = ni_func(ctx, sum + x,&err,ftyp);
            if (err)
                return(-1);
            f2 = ni_func(ctx, sum - x,&err,ftyp);
            if (err)
                return(-1);

            funct[j] = f1 + f2;                       
            acum += p[++i] * funct[j];
        }
        inew = iold + 1;

        *val = (acum + p[++i] * fzero) * diff;

        if (fabs(*val - *val1) <= rerr * fabs(*val)) {
            *nf = inew + iold;
            return(0);           /* convergence reached to requ. accuracy */
        }
    }
    *nf = inew + iold;
    return(1);                   /* no convergence to requested accuracy  */
}

/****************************************************************************/
/*  qtrap                                                                   */
/*                                                                          */
/*      Numerical integration, adapted from:                                */
/*      W.H. Press, B.P. Flannery, S.A. Teukolsky, W.T. Vetterling,         */
/*      Numerical Recipes in C. Cambridge Univ. Press 1988, p. 121          */
/*                                                                          */
/*      int qtrap(a,b,rerr,val,nf,ftyp)                                     */
/*                                                                          */
/*      double a,b          input:  integration interval                    */
/*      double rerr         input:  relative error                          */
/*      double *val         output: value of integral                       */
/*      int *nf             output: number of function calls                */
/*                                                                          */
/*      if ftyp = 0  use standard function,                                 */
/*                1  use user-defined rate models.                          */
/*               >1  other functions                                        */
/*                                                                          */
/*      Return:         0  if successful                                    */
/*                      1  if not successful                                */
/*                     -1  if error in function evaluation                  */

#define QTJMAX 50

int qtrap(TDAContext *ctx, double a,double b,double rerr,double *val,int *nf,int ftyp)
{
    int j,nnf,err;
    double olds,s;

    *nf = 0;
    olds = ctx->DBLMIN;    
    for (j = 1; j <= QTJMAX; ++j) {
        s = trapzd(ctx, a,b,j,&nnf,&err,ftyp);
        *nf += nnf;
        if (err)
            return(-1);
        if (fabs(s - olds) < rerr * fabs(olds)) {
            *val = s;
            return(0);
        }
        olds = s;
    }
    return(1);
}

/****************************************************************************/
/*  qsimp                                                                   */
/*                                                                          */
/*      Numerical integration, adapted from:                                */
/*      W.H. Press, B.P. Flannery, S.A. Teukolsky, W.T. Vetterling,         */
/*      Numerical Recipes in C. Cambridge Univ. Press 1988, p. 123          */
/*                                                                          */
/*      int qsimp(a,b,rerr,val,nf,ftyp)                                     */
/*                                                                          */
/*      double a,b          input:  integration interval                    */
/*      double rerr         input:  relative error                          */
/*      double *val         output: value of integral                       */
/*      int *nf             output: number of function calls                */
/*                                                                          */
/*      if ftyp = 0  use standard function,                                 */
/*                1  use user-defined rate models.                          */
/*                                                                          */
/*      Return:         0  if successful                                    */
/*                      1  if not successful                                */
/*                     -1  if error in function evaluation                  */

int qsimp(TDAContext *ctx, double a,double b,double rerr,double *val,int *nf,int ftyp)
{
    int j,err,nnf;
    double s,st,ost,os;

    *nf = 0;
    ost = os = ctx->DBLMIN;    
    for (j = 1; j <= QTJMAX; ++j) {
        st = trapzd(ctx, a,b,j,&nnf,&err,ftyp);
        *nf += nnf;
        if (err)
            return(-1);
        s = (4.0 * st - ost) / 3.0;
        if (fabs(s - os) < rerr * fabs(os)) {
            *val = s;
            return(0);
        }
        os = s;
        ost = st;
    }
    return(1);
}

/****************************************************************************/
/*  trapzd                                                                  */
/*                                                                          */
/*      Numerical integration, adapted from:                                */
/*      W.H. Press, B.P. Flannery, S.A. Teukolsky, W.T. Vetterling,         */
/*      Numerical Recipes in C. Cambridge Univ. Press 1988, p. 120          */
/*                                                                          */
/*      double trapzd(a,b,n,nf,err,ftyp)                                    */
/*                                                                          */
/*      double a,b          input:  integration interval                    */
/*      int n               input:  number of call                          */  
/*      int *nf             output: number of function calls                */
/*      int *err            output: 0 if OK, -1 if error in func evaluation */
/*                                                                          */
/*      if ftyp = 0  use standard function,                                 */
/*                1  use user-defined rate models.                          */
/*               >1  other functions                                        */
/*                                                                          */
/*      Return: value of integral                                           */

double trapzd(TDAContext *ctx, double a,double b,int n,int *nf,int *err,int ftyp)
{
    int j = 0;
    double x = 0.0,tnm = 0.0,sum = 0.0,del = 0.0,fa = 0.0,fb = 0.0;

    *nf = 0;

    if (n == 1) {
        ctx->s_trapzd_it = 1;
        fa = ni_func(ctx, a,err,ftyp);
        *nf += 1;
        if (*err)
            return(-1);
        fb = ni_func(ctx, b,err,ftyp);
        *nf += 1;
        if (*err)
            return(-1);

        ctx->s_trapzd_s = 0.5 * (b - a) * (fa + fb);
        return(ctx->s_trapzd_s);
    }
    else {
        tnm = (double)ctx->s_trapzd_it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        for (sum = 0.0, j = 1; j <= ctx->s_trapzd_it; j++, x += del) {
            fa = ni_func(ctx, x,err,ftyp);
            *nf += 1;
            if (*err)
                return(-1);
            sum += fa;
        }
        ctx->s_trapzd_it *= 2;
        ctx->s_trapzd_s = 0.5 * (ctx->s_trapzd_s + (b - a) * sum / tnm);
        return(ctx->s_trapzd_s);
    }
}

