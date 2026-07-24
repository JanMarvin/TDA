/****************************************************************************/
/*  t_gf                                                                    */
/*                                                                          */
/*  TDA. Program for Transition Data Analysis, written by Goetz Rohwer.     */
/*  Copyright (C) 1989,1991-94 Goetz Rohwer. All rights reserved.           */
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
#include "tda_context.h"

/*  functions in t_gf.c */

void clear_npflg(TDAContext *ctx);
void prn_npflg(TDAContext *ctx);
int imax(TDAContext *ctx, int i,int j);
int imin(TDAContext *ctx, int i,int j);
int iabs(TDAContext *ctx, int x);
double dmax(TDAContext *ctx, double x,double y);
double dmin(TDAContext *ctx, double x,double y);
double dsign(TDAContext *ctx, double x);
int isign(TDAContext *ctx, int x);
double rexp(TDAContext *ctx, double x);
double rlog(TDAContext *ctx, double x);
double loggam(TDAContext *ctx, double x,int *err);
double icgam(TDAContext *ctx, double x, double alpha,int *err);
double beta(TDAContext *ctx, double a, double b,int *err);
double incbeta(TDAContext *ctx, double x, double aa, double bb,int *err);
double digam(TDAContext *ctx, double x,int *err);
double trigam(TDAContext *ctx, double x,int *err);
int icdgam(TDAContext *ctx, double x,double p,double *d);
double bincoeff(TDAContext *ctx, int n,int  m);
int l_poisson(TDAContext *ctx, double theta,double k,double *val,double *val1,double *val2);
int l_negbin(TDAContext *ctx, double alpha,double gamma,double k,int opt,double *d);

/* ------------------------------------------------------------------------ */
/*  global variables defined in t_gf.c                                      */

                            /*  1 = overflow, underflow with exp, log.      */
                            /*  2 = overflow, underflow with fn             */
                            /*  3 = log likelihood positive                 */
                            /*  4 = Hessian not pos def.                    */
                            /*  5 = problem when evaluating time in Cox     */
                            /*      models.                                 */
                            /*  6 = error in range of probabilities         */
                            /*  7 = inaccuracy in incomplete gamma integral */     
                            /*  8 = inaccuracy in incomplete beta integral  */     
                            /*  9 = error in covariance matrix              */
                            /* 10 = integration did not always achieve      */
                            /*      required accuracy                       */

/* ------------------------------------------------------------------------ */
/*  clear_npflg()   Clear numerical problem flags.                          */

void clear_npflg(TDAContext *ctx)
{
    register int i;
    for (i = 1; i <= 20; ++i)
        ctx->NPFlgs[i] = 0;
}

/* ------------------------------------------------------------------------ */
/*  prn_npflgs()        Print information about numerical problems.         */

void prn_npflg(TDAContext *ctx)
{
    register int i,j,first = 1;

    j = 0;
    for (i = 1; i <= 20; ++i) {
        if (ctx->NPFlgs[i]) {
            if (first) {
                first = 0;
                printf1(ctx, "Numerical problems.\n");
            }
            switch (i) {
                case  1:            
                    printf1(ctx, "  1 : overflow/underflow in exp() or log().\n"); 
                    break;
                case  2:            
                    printf1(ctx, "  2 : overflow/underflow in log-likelihood.\n"); 
                    break;
                /**********************************************
                case  3:            
                    printf1(ctx, "  3 : log-likelihood positive.\n"); 
                    break;
                ************************************************/
                case  4:            
                    printf1(ctx, "  4 : Hessian not positive definite.\n"); 
                    break;
                case  5:            
                    printf1(ctx, "  5 : error in evaluating type 4 variable.\n"); 
                    break;
                case  6:            
                    printf1(ctx, "  6 : error in range of probabilities.\n"); 
                    break;
                case  7:            
                    printf1(ctx, "  7 : inaccuracy in incomplete gamma integral.\n"); 
                    break;
                case  8:            
                    printf1(ctx, "  8 : inaccuracy in incomplete beta integral.\n"); 
                    break;
                case  9:            
                    printf1(ctx, "  9 : error in covariance matrix.\n"); 
                    break;
                case 10:            
                    printf1(ctx, " 10 : integration not always achieved required accuracy.\n");
                    break;
                default:
                    break;
            }
            j++;
        }
    }
    if (j)
        newline(ctx);
}

/* ------------------------------------------------------------------------ */
/*  imax(i,j)    return maximum of i and j                                  */
                     
int imax(TDAContext *ctx, int i,int j)
{
    (void)ctx;        /* unused: the signature is shared */
    if (i >= j)
        return(i);
    else
        return(j);
}

/* ------------------------------------------------------------------------ */
/*  imin(i,j)   return minimum of i and j                                   */
                     
int imin(TDAContext *ctx, int i,int j)
{
    (void)ctx;        /* unused: the signature is shared */
    if (i < j)
        return(i);
    else
        return(j);
}

/* ------------------------------------------------------------------------ */
/*  dmax(x,y)   return maximum of x and y                                   */
                     
double dmax(TDAContext *ctx, double x,double y)
{
    (void)ctx;        /* unused: the signature is shared */
    if (x >= y)
        return(x);
    else
        return(y);
}

/* ------------------------------------------------------------------------ */
/*  dmin(x,y)   return minimum of x and y                                   */

double dmin(TDAContext *ctx, double x,double y)
{
    (void)ctx;        /* unused: the signature is shared */
    if (x < y)
        return(x);
    else
        return(y);
}

/* ------------------------------------------------------------------------ */
/*  dsign(x)    return sign(x).                                             */

double dsign(TDAContext *ctx, double x)
{
    (void)ctx;        /* unused: the signature is shared */
    if (x > 0.0)
        return(1.0);
    else if (x < 0.0)
        return(-1.0);
    else
        return(0.0);
}

/* ------------------------------------------------------------------------ */
/*  isign(x)    return sign(x).                                             */

int isign(TDAContext *ctx, int x)
{
    (void)ctx;        /* unused: the signature is shared */
    if (x > 0)
        return(1);
    else if (x < 0)
        return(-1);
    else
        return(0);
}

/* ------------------------------------------------------------------------ */
/*  iabs(x)                                                                 */

int iabs(TDAContext *ctx, int x)
{
    (void)ctx;        /* unused: the signature is shared */
    if (x >= 0)
        return(x);
    else
        return(-x);
}

/* ------------------------------------------------------------------------ */
/*  rlog(), rexp()                                                          */
/*      Local versions with range checking and corrections.                 */

double rexp(TDAContext *ctx, double x)
{
    if (x < ExpMin) {
        x = ExpMin;
        ctx->NPFlgs[1] += 1;
    }
    else if (x > ExpMax) {
        x = ExpMax;
        ctx->NPFlgs[1] += 1;
    }
    return(exp(x));
}

double rlog(TDAContext *ctx, double x)
{
    if (x < LogMin) {
        x = LogMin;
        ctx->NPFlgs[1] += 1;
    }
    return(log(x));
}

/* ------------------------------------------------------------------------ */
/*  loggam                                                                  */
/*      Log of the gamma function.                                          */
/*      Ref.: M.C. Pike, I.D. Hill, Logarithm of Gamma Function,            */
/*      Algorithm 291, Comm. ACM 9 (1966), p. 684                           */
/*                                                                          */
/*      Call: double loggam(x,err)                                          */
/*      Input: double x > 0.0                                               */
/*      Return: value of log gamma at x                                     */
/*                                                                          */
/*      Note: if x out of range err return a value != 0.                    */

double loggam(TDAContext *ctx, double x,int *err)
{
    (void)ctx;        /* unused: the signature is shared */
    double loggamma,f,z;

    if (x <= 0.0) {
        *err = -1;
        return(0.0);
    }
    *err = 0;
    if (x < 7.0) {
        f = 1.0;
        for (z = x; z < 7.0; z += 1.0) {
            x = z;
            f *= z;
        }
        x += 1.0;
        f = -log(f);
    }
    else {
        f = 0.0;
    }
    z = 1.0 / x * 1.0 / x;
    loggamma = f + (x - 0.5) * log(x) - x +
               0.918938533204673  +
           (((-0.000595238095238  * z +
               0.000793650793651) * z -
               0.002777777777778) * z +
               0.083333333333333) / x;

    return(loggamma);
}

/* ------------------------------------------------------------------------ */
/*  icgam                                                                   */
/*      Incomplete gamma integral.                                          */
/*      Ref.: Chi-Leung Lau, A Simple Series for the Incomplete Gamma       */
/*      Integral, AS 147, Applied Statistics 29 (1980), 113 - 114           */
/*                                                                          */
/*      Call: double icgam(x,alpha,&err)                                    */
/*      Input: double x     (argument, x > 0.0)                             */
/*      Input: double alpha (index of gamma function, alpha > 0.0)          */
/*      Return: value of the integral.                                      */
/*                                                                          */
/*      Note: if x out of range the function returns with err != 0.         */

double icgam(TDAContext *ctx, double x, double alpha,int *err)
{
    double lgam,c,f;   
    double eps = 1.e-6;     /* controls accuracy */

    *err = 0;
    if (x <= 0.0 || alpha <=  0.0) {
        *err = -1;
        return(0.0);
    }
    f = rexp(ctx, alpha * rlog(ctx, x) - loggam(ctx, alpha + 1.0,err) - x);
    if (*err)
        return(0.0);

    if (f < ctx->EPSI) {
        ctx->NPFlgs[7] += 1;
        return(1.0 - ctx->EPSI);
    }
    c = lgam = 1.0;

    while (1) {
        alpha += 1.0;
        c *= x / alpha;
        lgam += c;
        if ((c / lgam) <= eps)
            break;
    }
    return(lgam * f);
}

/****************************************************************************/
/*  beta    Beta function                                                   */

double beta(TDAContext *ctx, double a, double b,int *err)
{
    double x = 0.0;

    *err = 0;
    x = rexp(ctx, loggam(ctx, a,err) + loggam(ctx, b,err) - loggam(ctx, a + b,err));
    return(x);
}

/* ------------------------------------------------------------------------ */
/*  incbeta  Incomplete beta integral.                                      */
/*           Ref: L. Baker, C Mathematical Function Handbook,               */
/*           McGraw-Hill 1992, p. 580.                                      */
/*                                                                          */
/*  Note: accuracy is controlled with tol, currently set to 1.e-8.          */
/*  If the max number of iterations, itmax, is not sufficient to reach      */
/*  this accuracy, the numerical problems flag NPFlgs[8] is set.            */

double incbeta(TDAContext *ctx, double x, double aa, double bb,int *err)
{
    register int m;
    int itmax = 25;
    double dc,a,b,z,bmult,offset,ab,am1,ap1,d,h,c,fm,fm1,delta;
    double tol = 1.e-8;
    double small = 1.e-10;

    if (x < 0.0 || x > 1.0) {
        *err = -1;
        return(0.0);    
    }
    *err = 0;
    if (x == 0.0 || x == 1.0)
        return(x);

    bmult = rexp(ctx, loggam(ctx, aa + bb,err) - loggam(ctx, aa,err) - loggam(ctx, bb,err)
                                 + aa * rlog(ctx, x) + bb * rlog(ctx, 1.0 - x));
    /* was `if (err)`, which tests the pointer (never NULL here), so
       the function returned 0 for interior x.  With `*err` the result agrees with R's pbeta at every tested point.  Practical impact
       is likely small -- the operator apparently went unused. */
    if (*err)
        return(0);

    if (x < ((aa + 1.0) / (aa + bb + 2.0))) {
        a = aa;
        b = bb;
        z = x;
        offset = 0.0;
        bmult /= aa;
    }
    else {
        a = bb;
        b = aa;
        z = 1.0 - x;
        bmult = -bmult / bb;
        offset = 1.0;
    }
    ab = a + b;
    am1 = a - 1.0;
    ap1 = a + 1.0;
    d = 0.0;
    h = small;
    c = h;
    for (m = 0; m <= itmax; m++) {
        fm = (double)m;
        fm1 = 2.0 * fm;
        if (m)
            dc = fm * (b - fm) * z / ((am1 + fm1) * (a + fm1));
        else
            dc = 1.0;
        d = 1.0 + d * dc;
        c = 1.0 + dc / c;
        if (d == 0.0)
            d = small;
        if (c == 0.0)
            c = small;
        d = 1.0 / d;
        delta = d * c;
        h *= delta;
        dc = -(a + fm) * (ab + fm) * z / ((a + fm1) * (ap1 + fm1));
        d = 1.0 + d * dc;
        c = 1.0 + dc / c;
        if (d == 0.0)
            d = small;
        if (c == 0.0)
            c = small;
        d = 1.0 / d;
        delta = d * c;
        h *= delta;
        if (fabs(delta - 1.0) < tol) {
           return(h * bmult + offset);
        }
    }
    ctx->NPFlgs[8] += 1;
    return(h * bmult + offset);
}

/* ------------------------------------------------------------------------ */
/*  digam                                                                   */
/*      First derivative of log of gamma function.                          */
/*      Ref. J. M. Bernardo, Psi (Digamma) Function. Algorithm AS 103.      */
/*      Applied Statistics 25 (1976), pp. 315 - 317.                        */
/*                                                                          */
/*      Call: double digam(x,&err)                                          */
/*      Input: double x > 0.0                                               */
/*      Return: value of derivative of loggam(x).                           */
/*                                                                          */
/*      Note: if x out of range the function returns with err != 0.         */

double digam(TDAContext *ctx, double x,int *err)
{
    double r,y,dg;
    static double s = 1.e-5;
    static double c = 8.5;
    static double s3 = 8.333333333e-2;
    static double s4 = 8.333333333e-3;
    static double s5 = 3.968253968e-3;
    static double d1 = -0.5772156649;

    if (x <= 0.0) {
        *err = -1; 
        return(0.0);   
    }
    *err = 0;

    dg = 0.0;
    y = x;
    if (y <= s) {
        dg = d1 - 1.0 / y;
        return(dg);
    }

    while (y < c) {
        dg -= 1.0 / y;
        y += 1.0;
    }
    r = 1.0 / y;
    dg += rlog(ctx, y) - 0.5 * r;
    r *= r;
    dg -= r * (s3 - r * (s4 - r * s5));
    return(dg);
}

/* ------------------------------------------------------------------------ */
/*  trigam                                                                  */
/*      Second derivative of log of gamma function.                         */
/*      Ref. B. E. Schneider, Trigamma Function. Algorithm AS 121.          */
/*      Applied Statistics 27 (1978), pp. 97 - 99.                          */
/*                                                                          */
/*      Call: double trigam(x,&err)                                         */
/*      Input: double x > 0.0                                               */
/*      Return: value of scnd derivative of loggam(x).                      */
/*                                                                          */
/*      Note: if x out of range the function returns with err != 0.         */
/*      the program is exited.                                              */

double trigam(TDAContext *ctx, double x,int *err)
{
    (void)ctx;        /* unused: the signature is shared */
    double tg,y,z;
    static double a = 1.e-4;
    static double b = 5.0;
    static double b2 =  0.1666666667;
    static double b4 = -0.03333333333;
    static double b6 =  0.02380952381;
    static double b8 = -0.03333333333;

    if (x <= 0.0) {
        *err = -1;
        return(0.0);    
    }
    *err = 0;  

    tg = 0.0;
    z = x;
    if (z <= a) {
        tg = 1.0 / (z * z);
        return(tg);
    }
    while (z < b) {
        tg += 1.0 / (z * z);
        z += 1.0;
    }
    y = 1.0 / (z * z);
    tg += 0.5 * y + (1.0 + y * (b2 + y * (b4 + y * (b6 + y * b8)))) / z;
    return(tg);
}

/* ------------------------------------------------------------------------ */
/*  icdgam                                                                  */
/*      Derivatives of incomplete gamma integral.                           */
/*      Ref. R. J. Moore, Derivatives of the Incomplete Gamma Integral.     */
/*      Applied Statistics 31 (1982), pp. 330 - 335.                        */
/*                                                                          */
/*      Call:    int icdgam(x,p,d)                                          */
/*      Input: double x > 0.0                                               */
/*      Input: double p > 0.0                                               */
/*      Output: d[1] = deriv  I(x,p) / d  x                                 */
/*              d[2] = deriv2 I(x,p) / d2 x                                 */
/*              d[3] = deriv  I(x,p) / d  p                                 */
/*              d[4] = deriv2 I(x,p) / d2 p                                 */
/*              d[5] = deriv2 I(x,p) / dx dp                                */
/*              d[6] = I(x,p)                                               */
/*                                                                          */
/*      Return: 0 if OK, -1 if no convergence                               */
/*                                                                          */
/*      Note: if x out of range the function returns with -1                */
        
#define ICDTMAX 200

int icdgam(TDAContext *ctx, double x,double p,double *d)
{
    register int i,i2;
    int err;
    double pm1,xlog,gplog,psip,f,gp1log,dfp,dfpp,psip1,tmaxp,c,s;               
    double cp,cpp,dsp,dspp,a,b,cpc,term,so,an,psidp,psidp1;
    double pn[7],dp[7],dpp[7];

    static double e = 1.e-6;
    static double oflo = 1.e30;
    static double zero = 1.e-30;

    if (x <= 0.0 || p <= 0.0) {
        return(-1);   
    }
    err = 0;

    gplog  = loggam(ctx, p,&err);
    if (err)
        return(-1);

    gp1log = gplog + rlog(ctx, p);
    psip   = digam(ctx, p,&err);
    if (err)
        return(-1);

    psip1  = psip + 1.0 / p;
    psidp  = trigam(ctx, p,&err);
    if (err)
        return(-1);

    psidp1 = psidp - 1.0 / (p * p);

    pm1 = p - 1.0;
    xlog = rlog(ctx, x);
    d[1] = rexp(ctx, -gplog + pm1 * xlog - x);
    d[2] = d[1] * (pm1 / x - 1.0);
    d[5] = d[1] * (xlog - psip);

    if (x <= 1.0 || x < p) {

        f = rexp(ctx, p * xlog - gp1log - x);
        dfp = f * (xlog - psip1);
        dfpp = dfp * dfp / f - f * psidp1;

        tmaxp = (double)ICDTMAX + p;
        c = s = 1.0;
        cp = cpp = dsp = dspp = 0.0;
        a = p;

        while (1) {
            a += 1.0;
            cpc = cp / c;
            cp = cpc - 1.0 / a;
            cpp = cpp / c - cpc * cpc + 1.0 / (a * a);
            c = c * x / a;
            cp *= c;
            cpp = cpp * c + cp * cp / c;
            s += c;
            dsp += cp;
            dspp += cpp;
            if (a > tmaxp)
                return(-2);

            if (c <= e * s)
                break;
        }
        d[6] = s * f;
        d[3] = s * dfp + f * dsp;
        d[4] = s * dfpp + 2.0 * dfp * dsp + f * dspp;
        return(0);
    }

    f = rexp(ctx, p * xlog - gplog - x);
    dfp = f * (xlog - psip);
    dfpp = dfp * dfp / f - f * psidp;
    a = pm1;
    b = x + 1.0 - a;
    term = 0.0;
    pn[1] = 1.0;
    pn[2] = x;
    pn[3] = x + 1.0;
    pn[4] = x * b;
    so = pn[3] / pn[4];

    for (i = 1; i <= 4; ++i)  
        dp[i] = dpp[i] = 0.0;

    dp[4] = -x;

    while (1) {
        a -= 1.0;
        b += 2.0;
        term += 1.0;
        an = a * term;
        pn[5] = b * pn[3] + an * pn[1];
        pn[6] = b * pn[4] + an * pn[2];
        dp[5] = b * dp[3] - pn[3] + an * dp[1] + pn[1] * term;
        dp[6] = b * dp[4] - pn[4] + an * dp[2] + pn[2] * term;
        dpp[5] = b * dpp[3] + an * dpp[1] + 2.0 * (term * dp[1] - dp[3]);
        dpp[6] = b * dpp[4] + an * dpp[2] + 2.0 * (term * dp[2] - dp[4]);

        if (fabs(pn[6]) >= zero) {
            s = pn[5] / pn[6];
            c = fabs(s - so);
            if (c * p <= e) {
                if (c <= e * s)
                    break;
            }
            so = s;
        }
        for (i = 1; i <= 4; ++i) {
            i2 = i + 2;
            dp[i] = dp[i2];
            dpp[i] = dpp[i2];
            pn[i] = pn[i2];
        }
        if (term > (double)ICDTMAX)
            return(-1);

        if (fabs(pn[5]) >= oflo) {
            for (i = 1; i <= 4; ++i) {
                dp[i] /= oflo;
                dpp[i] /= oflo;
                pn[i] /= oflo;
            }
        }
    }

    d[6] = 1.0 - f * s;
    dsp = (dp[5] - s * dp[6]) / pn[6];
    dspp = (dpp[5] - s * dpp[6] - 2.0 * dsp * dp[6]) / pn[6];
    d[3] = -f * dsp - s * dfp;
    d[4] = -f * dspp - 2.0 * dsp * dfp - s * dfpp;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  bincoeff(n,m)       Binomial coefficient                                */
/*                      see: alg. ACM 19, 1960.                             */
 
double bincoeff(TDAContext *ctx, int n,int  m)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,b;
    double a;

    if (n < 0 || m < 0 || m > n)
        return(-1.0);
    if (m == 0 || m == n)
        return(1.0);

    a = 1.0;
    if (2 * m > n)
        b = n - m;
    else
        b = m;
    for (i = 0; i < b; ++i) {
        a *= (double)(n - i);
        if (a <= 0.0)
            return(-2);
        a /= (double)(i + 1);
    }
    return(a);
}

/* ------------------------------------------------------------------------ */
/*  l_poisson(theta,k,val,val1,val2)                                        */
/*                                                                          */
/*  Calculate logarithm of Poisson distribution.                            */
/*  Return: val  = value                                                    */
/*          val1 = first derivative wrt theta.                              */
/*          val2 = second derivative wrt theta.                             */
/*                                                                          */
/*  Return 0 if OK, -1 if theta <= 0 or k < 0.                              */
 
int l_poisson(TDAContext *ctx, double theta,double k,double *val,double *val1,double *val2)
{
    int err;

    if (theta <= 0.0 || k < 0.0)
        return(-1);

    *val = k * rlog(ctx, theta) - loggam(ctx, k + 1.0,&err) - theta;
    if (err)
        return(-1);

    *val1 = k / theta - 1.0;
    *val2 = -k / (theta * theta);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  l_negbin(alpha,gamma,k,opt,d)                                           */
/*                                                                          */
/*  Calculate logarithm of negbin distribution and, optionally, first       */
/*  and second derivatives wrt to alpha and gamma.                          */
/*                                                                          */
/*  d[0] = value                                                            */
/*  d[1] = derivative wrt alpha                                             */
/*  d[2] = derivative wrt alpha,alpha                                       */
/*  d[3] = derivative wrt gamma                                             */
/*  d[4] = derivative wrt gamma,gamma                                       */
/*  d[5] = derivative wrt alpha,gamma                                       */
/*                                                                          */
/*  First derivatives, if opt >= 1, first and second derivatives,           */
/*  if opt >= 2.                                                            */
/*                                                                          */
/*  Return 0 if OK, -1 if alpha <= 0 or gamma <= 0 or k < 0                 */
 
int l_negbin(TDAContext *ctx, double alpha,double gamma,double k,int opt,double *d)
{
    int err;
    double ak,ag,la,lg;

    if (alpha <= 0.0 || gamma <= 0.0 || k < 0.0)
        return(-1);

    ak = alpha + k;
    ag = alpha + gamma;
    la = rlog(ctx, alpha / ag);
    lg = rlog(ctx, gamma / ag);

    d[0] = loggam(ctx, ak,&err) - loggam(ctx, k + 1.0,&err) - loggam(ctx, alpha,&err);
    if (err)
        return(-1);
    d[0] += k * lg + alpha * la;

    if (opt >= 1) {
        d[1] = digam(ctx, ak,&err) - digam(ctx, alpha,&err);
        if (err)
            return(-1);
        d[1] += la + (gamma - k) / ag;

        d[3] = alpha * (k / gamma - 1.0) / ag;

        if (opt >= 2) {
            d[2] = trigam(ctx, ak,&err) - trigam(ctx, alpha,&err);
            if (err)
                return(-1);
            d[2] += (gamma / alpha - (gamma - k) / ag) / ag;

            d[4] = alpha * ((1.0 - k / gamma) / ag - k / (gamma * gamma)) / ag;
            d[5] = (k - gamma) / (ag * ag);
        }
    }
    return(0);
}




