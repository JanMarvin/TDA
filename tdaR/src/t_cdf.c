/****************************************************************************/
/*  t_cdf                                                                   */
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
#include "t_gf.h"
#include "t_lin.h"
#include "tda_context.h"

/*  functions in t_cdf */

double dnf(TDAContext *ctx, double x);
double cdnf(TDAContext *ctx, double z);
double cdnf1(TDAContext *ctx, double z);
double cdnif(TDAContext *ctx, double p);
double cdnif1(TDAContext *ctx, double p);
double cdchif(TDAContext *ctx, double x, int df);
double cdff(TDAContext *ctx, double x, int m, int n);
double cdtf(TDAContext *ctx, double t, int df);
double cdnm(TDAContext *ctx, double x);
int mecdf(TDAContext *ctx, int ndim,double *d,double *rho,double *prob);
double xphi(TDAContext *ctx, double x,double y);
int dmv(TDAContext *ctx, int m,int k,double *h,double *r,double *prob,double eps,double *err);
int dmv1(TDAContext *ctx, int m,int k,double *h,double *r,double *prob);
int dmvsolve(TDAContext *ctx, int mr,int m,int k,double *h,double *r,double *prob);
double bivn(TDAContext *ctx, double ah,double ak,double r,int *err);

/* ------------------------------------------------------------------------ */
/*  Global variables                                                        */


/* ------------------------------------------------------------------------ */
/*  dnf   Density of the standard normal distribution                       */
/*                                                                          */
/*        Call: double rdnf(x)                                              */
/*        Input: double x  argument                                         */
/*        Return: density at x.                                             */

static double SQ2Pi = 2.506628274631;    /* Sqrt of 2 * pi */ 

double dnf(TDAContext *ctx, double x)
{
    x *= x;
    x /= -2.0;
    x = rexp(ctx, x);
    x /= SQ2Pi;     
    return(x);
}

/* ------------------------------------------------------------------------ */
/*  cdnf    Normal distribution function.                                   */
/*          Ref. I.D. Hill, The Normal Integral (Algorithm AS 66), Applied  */
/*          Statistics 22 (1973), 424 - 427.                                */

double cdnf(TDAContext *ctx, double z)
{
    int up;
    double y,alnorm;

    static double con    =  1.28;
    static double ltone  =  7.70;
    static double utzero = 37.00;

    /*
    if (NDTyp) 
        return(cdnf1(ctx, z));
    */
    up = 0;     
    if (z < 0.0) {
        up = 1;
        z = -z;
    }
    if (z  <= ltone || (up && z <= utzero)) {
        y = 0.5 * z * z;
        if (z > con) {
            alnorm =  0.398942280385 * exp(-y) /
                     (z - 3.8052e-8      +  1.00000615302 /
                     (z + 3.98064794e-4  +  1.98615381364 /
                     (z - 0.151679116635 +  5.29330324926 /
                     (z + 4.8385912808   - 15.1508972451 /
                     (z + 0.742380924027 + 30.789933034 /
                     (z + 3.99019417011))))));
        }
        else {
            alnorm =  0.5 - z * (0.398942280444 - 0.399903438504 * y /
                     (y + 5.75885480458  - 29.8213557808 /
                     (y + 2.62433121679  + 48.6959930692 /
                     (y + 5.92885724438))));
        }
    }
    else
        alnorm = ctx->EPSI1;

    if (!up)
        alnorm = 1.0 - alnorm;
    if (alnorm < ctx->EPSI1)
        alnorm = ctx->EPSI1;   
  
    if (alnorm >= 1.0)
        alnorm = 1.0 - ctx->EPSI1;
    return(alnorm);
}

/* ------------------------------------------------------------------------ */
/*  cdnf1   Normal distribution function.                                   */
/*          Ref. ACM Algorithm 304.                                         */

double cdnf1(TDAContext *ctx, double z)
{
    (void)ctx;        /* unused: the signature is shared */
    register int upper;
    double y,z2,q1,q2,p1,p2,m,n,s,t,norm;

    if (z > 0.0)
        upper = 0;
    else
        upper = 1;

    z = fabs(z);
    z2 = z * z;
    y = 0.3989422804014327 * exp(-0.5 * z2);
    n = y / z;
    if (!upper && 1 - n == 1.0)
        norm = 1.0;
    else if (upper && n == 0.0)
        norm = 0.0;
    else {
        if ((upper && z > 2.32) || (!upper && z > 3.5)) {
            q1 = z;
            p2 = y * z;
            n = 1.0;
            p1 = y;
            q2 = z2 + 1.0;
            if (upper) {
                s = m = p1 / q1;
                t = p2 / q2;
            }
            else {
                s = m = 1.0 - p1 / q1;
                t = 1.0 - p2 / q2;
            }
            while (m != t && s != t) {
                n += 1.0;
                s = z * p2 + n * p1;
                p1 = p2;
                p2 = s;
                s = z * q2 + n * q1;
                q1 = q2;
                q2 = s;
                s = m;
                m = t;
                if (upper)
                    t = p2 / q2;
                else
                    t = 1.0 - p2 / q2;
            }
            norm = t;
        }
        else {
            z *= y;
            s = z;
            n = 1.0;
            t = 0.0;
            while (s != t) {
                n += 2.0;
                t = s;
                z = z * z2 / n;
                s += z;
            }
            if (upper)
                norm = 0.5 - s;
            else
                norm = 0.5 + s;
        }
    }
    return(norm);
}

/* ------------------------------------------------------------------------ */
/*  cdnif  Inverse of the normal distribution function.                     */
/*         Ref. J D. Beasley, S.G. Springer, The Percentage Points of the   */
/*         Normal Distribution (Algorithm AS 111), Applied Statistics 26,   */
/*         1977, 118 - 121.                                                 */
/*                                                                          */
/*  Note: we now use cdnif1(), see below, instead of cdnif().               */

double cdnif(TDAContext *ctx, double p)
{
    double q,r,pfnd;
    static double split =   0.42;
    static double a0    =   2.50662823884;
    static double a1    = -18.61500062529;
    static double a2    =  41.39119773534;
    static double a3    = -25.44106049637;
    static double b1    =  -8.47351093090;
    static double b2    =  23.08336743743;
    static double b3    = -21.06224101826;
    static double b4    =   3.13082909833;
    static double c0    =  -2.78718931138;
    static double c1    =  -2.29796479134;
    static double c2    =   4.85014127135;
    static double c3    =   2.32121276858;
    static double d1    =   3.54388924762;
    static double d2    =   1.63706781897;

    if (p <= 0.0 || p >= 1.0) { 
        printfe(ctx, "Error in cdnif with argument %g\n",p);
        gerr_exit(ctx, 27);
    }  
    q = p - 0.5;
    if (fabs(q) <= split) {
        r = q * q;
        pfnd = q * (((a3 * r + a2) * r + a1) * r + a0);
        pfnd /= ((((b4 * r + b3) * r + b2) * r + b1) * r + 1.0);
    }
    else {
        r = p;
        if (q > 0.0) r = 1.0 - p;
        r = sqrt(-log(r));
        pfnd = (((c3 * r + c2) * r + c1) * r + c0);
        pfnd /= ((d2 * r + d1) * r + 1.0);
        if (q < 0.0)
            pfnd = -pfnd;
    }
    return(pfnd);
}

/* ------------------------------------------------------------------------ */
/*  cdnif1   Inverse normal distribution function.                          */
/*           Ref. ACM Algorithm 442.                                        */

double cdnif1(TDAContext *ctx, double p)
{
    double x,y,s,ndev;

    if (p > 0.5)
        x = 1.0 - p;
    else
        x = p;
    if (x < 0.0) {
        printfe(ctx, "Error in cdnif1 with argument %g\n",p);
        gerr_exit(ctx, 28);
    }  
    s = sqrt(-2.0 * log(x));
    x = ((-7.49101 * s - 448.047) * s - 1266.846) /
            (((s + 109.8371) * s + 748.189) * s + 498.003) + s;
    if (p < 0.5)
        x = -x;
    y = p - cdnf(ctx, x);
    s = dnf(ctx, x);
    y /= s;
    s = x * x;
    /*****************************************************************
    ndev = (((((((((720.0 * s + 2556.0) * s + 1740.0) * s + 127.0) *
           (y / 7.0) +
                ((120.0 * s + 326.0) * s + 127.0) * x) * (y / 6.0) +
                (24.0 * s + 46.0) * s + 7.0) * (y / 40.0) +
                (0.75 * s + 0.875) * x) * y + s + 0.5) * (y / 3.0) +
                + x * 0.5) * y + 1.0) * y + x + 0.832E-24 * x;
    *****************************************************************/
    ndev = ((((((s * 0.6 + 1.15) * s + 0.175) * y + (s * 0.75 + 0.875) * x) *
            y + s + 0.5) * (y / 3.0) + x * 0.5) * y + 1.0) * y + x +
            0.42E-19 * x;

    return(ndev);
}

/* ------------------------------------------------------------------------ */
/*  cdchif  Chi square distribution function                                */
/*                                                                          */
/*      Ref. I.D. Hill, M.C. Pike, Chi-Squared Integral, Algorithm 299,     */
/*      Communications of the ACM 10 (1967), No. 4, 243 - 244.              */
/*                                                                          */
/*      double cdchif(x,df)                                                 */
/*      double x    argument                                                */
/*      int df      degrees of freedom                                      */
/*      Return      value of distribution function                          */
/*                                                                          */

#define const1 0.572364942925 
#define const2 0.564189583548 
#define max1   1000.0 
#define max2    200.0

#define ETA1 1.e-12

double cdchif(TDAContext *ctx, double x, int df)
{
    int even,bigx;
    double a,c,e,s,y = 0.0,z,chiprob;

    if (x < 0.0 || df < 1)
        return(-1.0);

    a = 0.5 * x;
    if (a < ETA1)
        return(0.0);

    if (a > max1)
        return(1.0);

    if (a > max2)
        bigx = 1;
    else
        bigx = 0;

    if (((df / 2) * 2) == df)
        even = 1;
    else
        even = 0;

    if (even || (df > 2 && !bigx))
        y = exp(-a);

    if (even)
        s = y;
    else
        s = 2.0 * cdnf(ctx, -sqrt(x));

    if (df > 2) {
        x = 0.5 * ((double)df - 1.0);
        if (even)
            z = 1.0;
        else
            z = 0.5;
        if (bigx) {
            if (even)
                e = 0.0;
            else
                e = const1;
            c = log(a);
            while (z <= x) {
                e += log(z);
                s += exp(c * z - a - e);
                z += 1.0;
            }
            chiprob = 1.0 - s;
        }
        else {
            if (even)
                e = 1.0;
            else
                e = const2 / sqrt(a);
            c = 0.0;
            while (z <= x) {
                e *= a / z;
                c += e;
                z += 1.0;
            }
            chiprob = 1.0 - (c * y + s);
        }
    }
    else
        chiprob = 1.0 - s;
    return(chiprob);
}

/* ------------------------------------------------------------------------ */
/*  cdff(x,m,n)   Distribution function of F distribution at value x        */
/*                with m,n degrees of freedom.                              */
/*                                                                          */
/*      Ref. E. Dorrer, F-Distribution, Algorithm 322, Communications       */
/*      of the ACM 11 (1968), 115. Certification: Comm. ACM (1969).         */

double cdff(TDAContext *ctx, double x, int m, int n)
{
    (void)ctx;        /* unused: the signature is shared */
    int i,j,a,b;
    double d,p,w,y,z;

    a = 2 * (m / 2) - m + 2;
    b = 2 * (n / 2) - n + 2;
    w = x * (double)m / (double)n;
    z = 1.0 / (1.0 + w);
    if (a == 1) {
        if (b == 1) {
            p = sqrt(w);
            y = 0.3183098862;   /* 1/pi */
            d = y * z / p;
            p = 2.0 * y * atan(p);
        }
        else {
            p = sqrt(w * z);
            d = 0.5 * p * z / w;
        }
    }
    else {
        if (b == 1) {
            p = sqrt(z);
            d = 0.5 * z * p;
            p = 1.0 - p;
        }
        else {
            d = z * z;
            p = w * z;
        }
    }
    y = 2.0 * w / z;
    for (j = b + 2; j <= n; j += 2) {
        d *= z * (1.0 + (double)a / ((double)j - 2.0));
        if (a == 1)
            p += d * y / ((double)j - 1.0);
        else
            p = (p + w) * z;
    }
    y = w * z;
    z = 2.0 / z;
    b = n - 2;
    for (i = a + 2; i <= m; i += 2) {
        j = i + b;
        d *= y * ((double)j / ((double)i - 2.0));
        p -= z * (d / (double)j);
    }
    if (p > 1.0)
        return(1.0);
    if (p < 0.0)
        return(0.0);
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  cdtf(t,df)  Distribution function of t distribution at value t with df  */
/*              degrees of freedom.                                         */
/*                                                                          */
/*      Ref. D.A. Levine, Student's t-Distribution, Algorithm 344, in:      */
/*      Communications of the ACM 12 (1969), 37 - 38. Vgl. auch die         */
/*      Remarks in: Comm. ACM 13 (1970), 124.                               */
/*                                                                          */
/*      Accurracy is about 6 significant digits.                            */

double cdtf(TDAContext *ctx, double t, int df)
{
    int i,n;
    double ans,d2,f1,f2,t1,t2;

    static double d1 = 0.63661977236758134;  /* 2/pi */

    if (df < 1) {
        printfe(ctx, "Error in cdtf: %d degrees of freedom not possible.\n",df);
        gerr_exit(ctx, 29);
    }
    t = fabs(t);
    t1 = t / sqrt((double)df);
    t2 = 1.0 / (1.0 + t1 * t1);

    if ((df / 2) * 2 == df) {
        d2 = t1 * sqrt(t2);
        ans = 1.0 - d2;
        if (df == 2)
            goto CDTFin;
        f1 = 1.0;
    }
    else {
        ans = 1.0 - d1 * atan(t1);
        if (df == 1)
            goto CDTFin;
        d2 = d1 * t1 * t2;
        ans -= d2;
        if (df == 3)
            goto CDTFin;

        f1 = 0.0;
    }
    n = (df - 2) / 2;
    for (i = 1; i <= n; ++i) {
        f2 = 2.0 * (double)i - f1;
        d2 *= t2 * f2 / (f2 + 1.0);
        ans -= d2;
    }
CDTFin:
    if (ans <= 0.0)
        return(1.0);
    else
        return(1.0 - ans / 2.0);
}

/* ------------------------------------------------------------------------ */
/*  cdnm                                                                    */
/*      Reciprocal of Mill's ratio: f(x) / (1 - F(x)), where f is the       */
/*      density and F the distribution function of the standard normal      */
/*      distribution.                                                       */
/*                                                                          */
/*      Ref. A. V. Swan, The Reciprocal of Mill's Ratio, Algorithm AS 17,   */
/*      Applied Statistics 18 (1969), 115 - 116.                            */
/*                                                                          */
/*      Call: double cdnm(x)                                                */
/*      Input: double x                                                     */
/*      Return: value of function at x.                                     */
/*                                                                          */
/*      If x < -22.9, this function is of order 10.E-114 and the return     */
/*      value is set to zero.                                               */
/*                                                                          */

static double FPi = 1.2533141373155;  /* sqrt(pi/2) */
static double Least = -22.9;

double cdnm(TDAContext *ctx, double x)
{
    double a,a0,a1,a2,b,b0,b1,b2,r,s,t,d,rmill;

    if (x == 0.0) 
        rmill = 1.0 / FPi;
    else if (x < Least)
        rmill = 0.0;
    else {
        d = 1.0;
        if (x < 0.0)
            d = -1.0;
        x = fabs(x);
        if (x <= 2.0) {
            s = 0.0;
            a = 1.0;
            r = t = x;
            b = x * x;
            while (fabs(s-t) > ctx->EPSI2) {
                a += 2.0;
                s = t;
                r *= b / a;
                t += r;
            }
            rmill = 1.0 / (FPi * exp(0.5 * b) - d * t);
        }
        else {
            a = 2.0;
            r = s = b1 = x;
            a1 = x * x + 1.0;
            a2 = x * (a1 + 2.0);
            b2 = a1 + 1.0;
            t  = a2 / b2;
            while (fabs(r-t) > ctx->EPSI2 && fabs(s-t) > ctx->EPSI2) {
                a += 1.0;
                a0 = a1;
                a1 = a2;
                a2 = x * a1 + a * a0;
                b0 = b1;
                b1 = b2;
                b2 = x * b1 + a * b0;
                r  = s;
                s  = t;
                t  = a2 / b2;
            }
            rmill = t;
            if (d < 0.0)
                rmill /= (2.0 * FPi * exp(0.5 * x * x) * t - 1.0);
        }
    }
    return(rmill);
}

/* ------------------------------------------------------------------------ */
/*  mecdf(ndim,d,rho,prob,ier)                                              */
/*                                                                          */
/*  Calculates the cum distribution function for a multivariate normal      */
/*  distribution with ndim dimensions, that is:                             */
/*                                                                          */
/*  prob = Prob(X1 > d1,...,Xn > dn), n = ndim.                             */
/*                                                                          */
/*  rho[] is the lower triangle of the correlation matrix.                  */
/*                                                                          */
/*  The function uses the Mendell-Elston procedure as described in          */
/*  Kamakura (1989). We follow the code given in ACM algorithm 717.         */
/*  There is the following note:                                            */
/*  NOTE:  Equation (15) in Kamakura has an error.                          */
/*                                                                          */

#define MVNNDMax 21     /* max number of dimensions is MVNNDMax - 1 */

int mecdf(TDAContext *ctx, int ndim,double *d,double *rho,double *prob)
{
    register int i,j,k,ir;
    double probi,tmp;
    double u[MVNNDMax],uumz[MVNNDMax],sig[MVNNDMax][MVNNDMax];
    double z[MVNNDMax][MVNNDMax],r[MVNNDMax][MVNNDMax][MVNNDMax];

    if (ndim >= MVNNDMax)      
        return(-1);

    ir = 0;
    for (i = 1; i <= ndim; ++i) {
        z[i][0] = d[i];
        for (j = 1; j < i; ++j)  
            r[j][i][0] = rho[++ir];
    }
    *prob = 1.0 - cdnf(ctx, z[1][0]);           
    if (*prob <= 0.0)         
        return(1);
                   
    u[1] = xphi(ctx, z[1][0],0.0) / *prob;
    uumz[1] = u[1] * (u[1] - z[1][0]);

    for (i = 2; i <= ndim; ++i) {
        for (j = 1; j < i; ++j) {
            for (k = 1; k < j; ++k) {
                tmp = r[j][i][k - 1] - r[k][j][k - 1] * r[k][i][k - 1] * uumz[k];
                r[j][i][k] = tmp / sig[j][k] / sig[i][k]; 
            }
            sig[i][j] = sqrt(1.0 - uumz[j] * r[j][i][j - 1] * r[j][i][j - 1]);
            z[i][j] = (z[i][j - 1] - u[j] * r[j][i][j - 1]) / sig[i][j];
        }
        probi = 1.0 - cdnf(ctx, z[i][i - 1]);
        if (probi <= 0.0)           
            return(i);
                  
        *prob *= probi;         
        if (i < ndim) {       
            u[i] = xphi(ctx, z[i][i - 1],0.0) / probi;
            uumz[i] = u[i] * (u[i] - z[i][i - 1]);
        }       
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  xphi(x,y)                                                               */
       

double xphi(TDAContext *ctx, double x,double y)
{
    double arg,phi = 0.0;

    arg = -0.5 * x * x - ctx->SQ2P - y;
    phi = rexp(ctx, arg);
    return(phi);
}

/* ------------------------------------------------------------------------ */
/*  dmv                                                                     */
/*  ##                                                                      */
/*  Multivariate normal integral, adapted from:                             */
/*  Z. Drezner, Computation of the Multivariate Normal Integral.            */
/*  ACM Transactions on Mathematical Software 18 (1992), 470 - 480.         */
/*  ACM Algorithm 725.                                                      */
/*                                                                          */
/*  err = dmv (m,k,h,r,prob,eps,err)                                        */
/*                                                                          */
/*  m      input   =  dimension (2 <= m <= DMVMax)                          */        
/*  k      input   =  specifies type of algorithm.                          */
/*                    if 2 <= k <= 10, one Gaussian quadrature with k       */
/*                    points for each variable is performed.                */
/*                    eps is disregarded in this case.                      */
/*                                                                          */
/*                    if 12 <= k <= 20, progressive Gaussian quadratures    */  
/*                    are performed for ki = 2,...,kmax = k - 10.           */
/*                    until the difference between two successive           */
/*                    calculations does not exceede eps. The last value     */
/*                    of k is then returned. If no convergence, then        */
/*                    the last difference is returned in err.               */
/*                                                                          */
/*  h      input   =  vector of upper limits                                */
/*  r      input   =  correlation matrix                                    */
/*  prob   output  =  probability                                           */
/*  eps    input   =  required accuracy.                                    */
/*  err    output  =  difference in prob between two successive             */  
/*                    calculations, only used if 12 <= k <= 20.             */
/*                                                                          */
/*  return ier =      0 if successful return (2 <= k <= 10, or if           */
/*                      12 <= k <= 20, then if eps is satisfied).           */
/*                    1 if  12 <= k <= 20 and eps not satisfied.            */
/*                   -1 if error in dimension m                             */
/*                   -2 if r is singular or has negative determinant.       */  

#define DMVCUT -15.0    /* if x < CUT) then exp(x) = 0 */

int dmv(TDAContext *ctx, int m,int k,double *h,double *r,double *prob,double eps,double *err)
{
    int ier,kmax,ker;
    double pold;

    *err = *prob = 0.0;
    if (k >= 2 && k <= 10) {
        ier = dmv1(ctx, m,k,h,r,prob);
        return(ier);
    }  
    kmax = 10;
    if (k > 11 && k < 21)
        kmax = k - 10;

    pold = -1.0;
    for (ker = 2; ker <= kmax; ++ker) {
        ier = dmv1(ctx, m,ker,h,r,prob);
        if (ier)         
            return(ier);

        *err = *prob - pold;
        if (fabs(*err) <= eps)
            return(0);              
        pold = *prob;
    }
    return(1);
}           

/* ------------------------------------------------------------------------ */
/*  dmv1(m,k,h,r,prob)                                                      */        
/*                                                                          */
/*  Calculates the probability for a given k.                               */
/*  k is the number of points in the Gauss quadrature.                      */
/*  if k < 2 or k > 10 then k is set to 10.                                 */
/*                                                                          */
/*  Return:      0 if successful                                            */
/*              -1 if error in dimension.                                   */
/*              -2 if r singular or negative determinant.                   */

int dmv1(TDAContext *ctx, int m,int k,double *h,double *r,double *prob)
{
    register int i,j,l,ij;
    int ier,ipos,l1,l2;
    int list[DMVMax + 1],list1[DMVMax + 1];
    double s,h1[DMVMax + 1],r1[DMVMax * DMVMax + 1];

    if (m < 0 || m > 20)           
        return(-1);

    *prob = 0.0;
    for (i = 1; i <= m; ++i) {
        list[i] = 0;
        if (h[i] > 0.0)
            list[i] = 1;
    }

L20:
    ipos = 0;
    l = 0;
    for (i = 1; i <= m; ++i) {
        if (list[i] != 2) {           
            l++;
            list1[l] = i;
        }
    }
    for (i = 1; i <= l; ++i) { 
        l1 = list1[i];
        h1[i] = h[l1];
        if (list[l1] == 1) {
            h1[i] = -h1[i];
            ipos = 1 - ipos;
        }
        for (j = 1; j <= l; ++j) {
            l2 = list1[j];
            ij = (i - 1) * m + j;
            r1[ij] = r[(l1 - 1) * m + l2];
            if (list[l1] == 1)
                r1[ij] = -r1[ij];
            if (list[l2] == 1)
                r1[ij] = -r1[ij];
        }
    }
    ier = dmvsolve(ctx, m,l,k,h1,r1,&s);
    if (ier)
        return(ier);

    if (ipos == 1)
        s = -s;
    *prob += s;     
    for (i = 1; i <= m; ++i) {
        if (list[i] != 0) {           
            list[i] += 1;         
            if (list[i] == 2)
                goto L20;
            list[i] = 1;
        }
    }
    return(0);
}          

static double MVNCOEF[10][10] = {
   {0.0,              0.0,              0.0,
    0.0,              0.0,              0.0,
    0.0,              0.0,              0.0,
    0.0},
   {0.6405291796843790,0.245697745768379,0.0,           
    0.0,              0.0,               0.0,
    0.0,              0.0,               0.0,
    0.0},
   {0.446029770466658,0.396468266998335,4.37288879877644e-2,
    0.0,              0.0,              0.0,
    0.0,              0.0,              0.0,
    0.0},                                         
   {0.325302999756919,0.421107101852062,0.133442500357520,
    0.637432348625728e-2,0.0,           0.0,
    0.0,                 0.0,           0.0,
    0.0},
   {0.248406152028443,0.392331066652399,0.211418193076057,
    0.332466603513439e-1,0.824853344515628e-3,0.0,         
    0.0,                 0.0,                 0.0,
    0.0},
   {0.196849675488598,0.349154201525395,0.257259520584421,
    0.0760131375840057,0.685191862513596e-2,9.84716452019267e-5,
    0.0,               0.0,                 0.0,
    0.0},                                                       
   {0.160609965149261,0.306319808158099,0.275527141784905,
    0.120630193130784,0.0218922863438067,0.00123644672831056,
    1.10841575911059e-5,0.0,             0.0,
    0.0},
   {0.13410918845336,0.26833075447264,0.275953397988422,   
    0.15744828261879,0.0448141099174625,0.00536793575602526,  
    2.02063649132407e-4,1.19259692659532e-6,0.0,
    0.0},
   {0.114088970242118,0.235940791223685,0.266425473630253,
    0.183251679101663,0.0713440493066916,0.0139814184155604,
    0.00116385272078519,0.305670214897831e-4,1.23790511337496e-7,
    0.0},
   {0.0985520975191087,0.208678066608185,0.252051688403761,
    0.198684340038387,0.097198422760062,0.0270244164355446,
    0.00380464962249537,2.28886243044656e-4,4.34534479844469e-6,
    1.24773714817825e-8} 
    };

static double MVNX[10][10] = {
   {0.0,              0.0,              0.0,
    0.0,              0.0,              0.0,
    0.0,              0.0,              0.0,
    0.0},
   {0.300193931060839,1.25242104533372,0.0,       
    0.0,              0.0,             0.0,
    0.0,              0.0,             0.0,
    0.0},
   {0.190554149798192,0.848251867544577,1.79977657841573,         
    0.0,              0.0,               0.0,
    0.0,              0.0,               0.0,
    0.0},
   {0.133776446996068,0.624324690187190,1.34253782564499,
    2.26266447701036,0.0,               0.0,
    0.0,             0.0,               0.0,
    0.0},
   {0.100242151968216,0.482813966046201,1.06094982152572, 
    1.77972941852026,2.66976035608766,0.0,
    0.0,             0.0,             0.0,
    0.0},
   {0.0786006594130979,0.386739410270631,0.866429471682044,
    1.46569804966352,2.172707796939,3.03682016932287,          
    0.0,             0.0,           0.0,
    0.0},
   {0.0637164846067008,0.318192018888619,0.724198989258373,
    1.23803559921509,1.83852822027095,2.53148815132768,
    3.37345643012458,0.0,             0.0,
    0.0},
   {0.0529786439318514,0.267398372167767,0.616302884182402,  
    1.06424631211623,1.58885586227006,2.18392115309586,   
    2.86313388370808,3.6860071627244,0.0,
    0.0},
   {0.0449390308011934,0.228605305560535,0.532195844331646,
    0.927280745338081,1.39292385519588,1.91884309919743,
    2.50624783400574,3.17269213348124,3.97889886978978,         
    0.0},
   {0.0387385243257289,0.198233304013083,0.465201111814767,  
    0.816861885592273,1.23454132402818,1.70679814968913,
    2.22994008892494,2.80910374689875,3.46387241949586,
    4.25536180636608}   
};


/* ------------------------------------------------------------------------ */
/*  ier = dmvsolve(mr,m,k,h,r,prob)                                         */
/*  ##                                                                      */
/*  calculates the multivariate probability of a given vector h             */
/*  and for given k; mr is dimension of r[]                                 */
/*                                                                          */
/*  Return:      0 if successful                                            */
/*              -1 if error in dimension.                                   */
/*              -2 if r singular or negative determinant.                   */
 
int dmvsolve(TDAContext *ctx, int mr,int m,int k,double *h,double *r,double *prob)
{
    register int i,j,ij,ii;
    int ier,l;
    int list[DMVMax + 1];
    double pi,solve,det,prod,sum,tmp;
    double h1[DMVMax + 1],h2[DMVMax + 1],y[DMVMax + 1];
    double d[3],r1[DMVMax * DMVMax + 1];

    *prob = 0.0;
    pi = 4.0 * atan(1.0);

    ier = 0;
    if (m < 0 || m > DMVMax)
        return(-1);

    if (k < 2 || k > 10)
        k = 10;
    if (m == 0) {
        *prob = 1.0;
        return(0);
    }
       
    for (i = 1; i <= m; ++i) {
        for (j = 1; j <= m; ++j) {
            ij = (i - 1) * mr + j;
            r1[ij] = r[ij];
        }
    }
    ier = dpofa(ctx, r1,mr,m);
    if (ier)         
        return(-2);

    ier = dpodi(ctx, r1,mr,m,d,11);
    if (ier || d[1] <= 0.0 || d[2] <= -30)           
        return(-2);

    det = d[1] * pow(10.0,d[2]);
    for (i = 1; i <= m; ++i) {
        for (j = i; j <= m; ++j)  
            r1[(j - 1) * mr + i] = r1[(i - 1) * mr + j];
    }
    prod = 1.0;
    for (i = 1; i <= m; ++i) {
        ii = (i - 1) * mr + i;
        if (r1[ii] == 0.0)
            return(-2);
        h1[i] = sqrt(2.0 / r1[ii]);
        if (h1[i] == 0.0)
            return(-2);
        h2[i] = h[i] / h1[i];
        prod *= pi * r1[ii];  
        list[i] = 1;
    }
    prod = 1.0 / sqrt(det * prod); 

    for (i = 1; i <= m; ++i) {
        for (j = 1; j <= m; ++j)  
            r1[(i - 1) * mr + j] *= h1[i] * h1[j] * 0.5;   
    }
    solve = 0.0;

L60:
    sum = 0.0;
    for (i = 1; i <= m; ++i) {
        l = list[i];
        tmp = MVNX[k - 1][l - 1];     
        sum += tmp * tmp;                        
        y[i] = tmp - h2[i];
    }
    for (i = 1; i <= m; ++i) {
        for (j = 1; j <= m; ++j)  
            sum -= y[i] * y[j] * r1[(i - 1) * mr + j];
    }
    if (sum >= DMVCUT) {           
        sum = exp(sum);
        for (i = 1; i <= m; ++i) {
            sum *= MVNCOEF[k - 1][list[i] - 1];       
        }
        solve += sum;       
    }
    for (i = 1; i <= m; ++i) {
        list[i] += 1;         
        if (list[i] <= k)              
            goto L60;
        list[i] = 1;
    }
    solve *= prod;
    *prob = solve;
    return(0);
}           

/* ------------------------------------------------------------------------ */
/*  bivn(ah,ak,r,err)   Calculates the value of the bivariate normal        */  
/*                      distribution with correlation r, at values ah,ak.   */
/*                                                                          */
/*  Source: T.G. Donnelly, Algorithm 462. Bivariate Normal Distribution.    */
/*                                                                          */
/*  Return: lower tail value; *err = 0 if OK, -1 if error.                  */
  
double bivn(TDAContext *ctx, double ah,double ak,double r,int *err)
{
    int is;
    double gh,gk,rr,b,wh,wk,gw,sgn,t,g2,h2,a2,h4,ex,w2,ap,sp,s2,s1,sn;
    double con,sqr,conex,cn;
    double pi2 = 6.283185307179587;

    /* int idig = 15; */

    if (r < -1.0 || r > 1.0) {
        *err = -1;
        return(0.0);
    }
    *err = 0;
    b = 0.0;
    gh = cdnf1(ctx, ah) / 2.0;
    gk = cdnf1(ctx, ak) / 2.0;
    ah = -ah;
    ak = -ak;

    if (r == 0.0) {
        b = 4.0 * gh * gk;
        goto Fin;
    }
    rr = 1.0 - r * r;
    if (rr == 0.0) {
        if (r < 0.0) {
            if (ah + ak < 0.0)  
                b = 2.0 * (gh + gk) - 1.0;
        }
        else {
            if (ah < ak)  
                b = 2.0 * gk;
            else
                b = 2.0 * gh;
        }
        goto Fin;
    }
    sqr = sqrt(rr);
    con = pi2 * 1.0e-15 / 2.0;

    /*************************
    if (idig == 15)  
        con = pi2 * 1.0e-15 / 2.0;
    else {
        con = pi2 / 2.0;
        for (i = 1; i <= idig; ++i)
            con /= 10.0;
    }
    ******************************/

    if (ah == 0.0) {
        if (ak == 0.0) {
            b = atan(r / sqr) / pi2 + 0.25;
            goto Fin;
        }
    }
    else {
        b = gh;
        if (ah * ak == 0.0)
            goto L200;
   
        else if (ah * ak < 0.0)
            b -= 0.5;
    }
    b += gk;
    if (ah == 0.0) {
        wh = -ak;
        wk = (ah / ak - r) / sqr;
        gw = 2.0 * gk;
        is = 1;
        goto L210;
    }
L200:
    wh = -ah;
    wk = (ak / ah - r) / sqr;
    gw = 2.0 * gh;
    is = -1;
L210:
    sgn = -1.0;
    t = 0.0;
    if (wk == 0.0)
        goto L320;

    if (fabs(wk) == 1.0) {
        t = wk * gw * (1.0 - gw) / 2.0;
        goto L310;
    }
    else if (fabs(wk) > 1.0) {
        sgn = -sgn;
        wh *= wk;
        g2 = cdnf1(ctx, wh);
        wk = 1.0 / wk;

        if (wk < 0.0)
            b += 0.5;
     
        b = b - (gw + g2) / 2.0 + gw * g2;
    }
    h2 = wh * wh;
    a2 = wk * wk;
    h4 = h2 / 2.0;
    ex = exp(-h4);
    w2 = h4 * ex;
    sp = ap = 1.0;
    s2 = ap - ex;
    sn = s1 = 0.0;
    conex = fabs(con / wk);
    goto L290;

L280:
    sn = sp;
    sp += 1.0;
    s2 -= w2;
    w2 *= h4 / sp;
    ap *= -a2;
L290:
    cn = ap * s2 / (sn + sp);
    s1 += cn;
    if (fabs(cn) > conex)
        goto L280;

    t = (atan(wk) - wk * s1) / pi2;
L310:
    b += sgn * t;
L320:
    if (is < 0) {
        if (ak != 0.0) {
            wh = -ak;
            wk = (ah / ak - r) / sqr;
            gw = 2.0 * gk;
            is = 1;
            goto L210;
        }
    }
Fin:
    if (b <= 0.0)
        b = ctx->EPSI;
    if (b >= 1.0)
        b = 1.0 - ctx->EPSI;
    return(b);
}


