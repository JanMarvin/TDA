#include "tda_rhooks.h"
/****************************************************************************/
/*  t_spl                                                                   */
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
#include "tda_context.h"

/*  functions in t_spl.c */

int spltena(TDAContext *ctx, int n,double *x,double *y,double *xp,double *yp, int mode,double slp1,double slpn,double sigma);
void spltenad(TDAContext *ctx, int n,double *x,double *y,double *xp,double *yp, double sigma,double t,double *tx,double *ty,int first);
int spltenc(TDAContext *ctx, int n,double *x,double *y,double *xp,double *yp,double sigma);
void spltencd(TDAContext *ctx, int n,double *x,double *y,double *xp,double *yp, double sigma,double t,double *tx,double *ty,int first);
int splakima(TDAContext *ctx, int mode,int l,double *x,double *y,int m,int n,double *u,double *v);
int splf(TDAContext *ctx, int m,double *x,double *y,double *w,double xa,double xb, int k,double s,int *nn,double * t,double *c,int maxit);
int bandet(TDAContext *ctx, double *g, int ls, int nk1);
void bansol(TDAContext *ctx, double *g, int ls, int nk1, double *z, double *c);
double spld(TDAContext *ctx, int n, double *t,int k,double *c,int nu,double arg,int l,double *tmp);
int interv(TDAContext *ctx, int n,double *x,double arg, int init);

/* ------------------------------------------------------------------------ */
/*  spltena                                                                 */
/*                                                                          */
/*      Calculates a spline under tension for an arbitrary series of        */
/*      points (x,y).                                                       */
/*                                                                          */
/*      Adopted from:                                                       */
/*      A.K. Cline, Six Subprograms for Curve Fitting Using Splines         */
/*      Under Tension (Algorithm 476), Communications of the ACM, Vol. 17,  */
/*      1974, No. 4, 220 - 223.                                             */
/*                                                                          */
/*      Note: The curve is mapped onto the interval [0,1]. To calculate     */
/*      smoothed points used the function splenad().                        */
/*                                                                          */
/*      Call:  int spltena(n,x,y,xp,yp,mode,slp1,slp2,sigma)                */
/*                                                                          */
/*      int n        input:  number of data points, n >= 2                  */
/*      double *x    input:  array with x coordinates, x(i), i = 1,n        */
/*                   output: unchanged                                      */
/*      double *y    input:  array with y coordinates, y(i), i = 1,n        */
/*                   output: unchanged                                      */
/*      double *xp   input:  array xp(i), 0 = 1,n                           */
/*                   output: description of smooth curve                    */
/*                           note: xp[0] contains the polygonal arclength   */
/*                           of the spline, to be used by spltenad().       */
/*      double *yp   input:  array yp(i), 0 = 1,n                           */
/*                   output: description of smooth curve                    */
/*      int mode     input:  mode                                           */
/*                           0 internal calculation of slp1 and slpn.       */
/*                           1 slp1 und slpn are given as input parameters. */
/*      double spl1  input:  derivative at x(1), if mode = 1                */
/*      double spln  input:  derivative at x(n), if mode = 1                */
/*      double sigma input:  tension factor, sigma > 0.0                    */
/*                                                                          */
/*      Return:      0  if successful                                       */
/*                  -1  error in input parameters                           */
/*                  -2  insufficient memory                                 */

/* ------------------------------------------------------------------------ */
/*  spltenad                                                                */
/*      Calculation of smooths points; to be used in connection with        */
/*      spltena().                                                          */
/*                                                                          */
/*      Adopted from:                                                       */
/*      A.K. Cline, Six Subprograms for Curve Fitting Using Splines         */
/*      Under Tension (Algorithm 476), Communications of the ACM, Vol. 17,  */
/*      1974, No. 4, 220 - 223.                                             */
/*                                                                          */
/*      Call: spltenad(n,x,y,xp,yp,sigma,t,tx,ty,first)                     */
/*                                                                          */
/*      int n        input:  number of data points, n >= 2                  */
/*      double *x    input:  array with x coordinates, x(i), i = 1,n        */
/*                   output: unchanged                                      */
/*      double *y    input:  array with y coordinates, y(i), i = 1,n        */
/*                   output: unchanged                                      */
/*      double *xp   input:  description of smooth curve                    */
/*                           xp[0] must contain the polygonal arclength     */
/*                           of the spline.                                 */
/*                   output: unchanged                                      */
/*      double *yp   input:  description of smooth curve                    */
/*                   output: unchanged                                      */
/*      double sigma input:  tension factor, sigma > 0.0                    */
/*      double t     input:  abcissa in [0,1] to be used for calculation    */
/*                           of corresponding (xt,yt)                       */
/*      double *tx   output: x coordinate of smoothed curve at t            */
/*      double *ty   output: y coordinate of smoothed curve at t            */
/*      int first    input:  first = 1 when the function is called first;   */
/*                           afterwards: first = 0.                         */

/* ------------------------------------------------------------------------ */
/*  spltenc                                                                 */
/*      Calculates a closed spline under tension for an arbitrary series    */
/*      of points (x,y).                                                    */
/*                                                                          */
/*      Adopted from:                                                       */
/*      A.K. Cline, Six Subprograms for Curve Fitting Using Splines         */
/*      Under Tension (Algorithm 476), Communications of the ACM, Vol. 17,  */
/*      1974, No. 4, 220 - 223.                                             */
/*                                                                          */
/*      Note: The curve is mapped onto the interval [0,1]. To calculate     */
/*      smoothed points used the function splenac().                        */
/*                                                                          */
/*      Call: int spltenc(n,x,y,xp,yp,sigma)                                */
/*                                                                          */
/*      int n        input:  number of data points, n >= 2                  */
/*      double *x    input:  array with x coordinates, x(i), i = 1,n        */
/*                   output: unchanged                                      */
/*      double *y    input:  array with y coordinates, y(i), i = 1,n        */
/*                   output: unchanged                                      */
/*      double *xp   input:  array xp(i), 0 = 1,n                           */
/*                   output: description of smooth curve                    */
/*                           note: xp[0] contains the polygonal arclength   */
/*                           of the spline, to be used by spltenad().       */
/*      double *yp   input:  array yp(i), 0 = 1,n                           */
/*                   output: description of smooth curve                    */
/*      double sigma input:  tension factor, sigma > 0.0                    */
/*                                                                          */
/*      Return:      0  if successful                                       */
/*                  -1  error in input parameters                           */
/*                  -2  insufficient memory                                 */

/* ------------------------------------------------------------------------ */
/*  spltencd                                                                */
/*      Calculation of smooths points; to be used in connection with        */
/*      spltenc().                                                          */
/*                                                                          */
/*      Adopted from:                                                       */
/*      A.K. Cline, Six Subprograms for Curve Fitting Using Splines         */
/*      Under Tension (Algorithm 476), Communications of the ACM, Vol. 17,  */
/*      1974, No. 4, 220 - 223.                                             */
/*                                                                          */
/*      Call: spltencd(n,x,y,xp,yp,sigma,t,tx,ty,first)                     */
/*                                                                          */
/*      int n        input:  number of data points, n >= 2                  */
/*      double *x    input:  array with x coordinates, x(i), i = 1,n        */
/*                   output: unchanged                                      */
/*      double *y    input:  array with y coordinates, y(i), i = 1,n        */
/*                   output: unchanged                                      */
/*      double *xp   input:  description of smooth curve                    */
/*                           xp[0] must contain the polygonal arclength     */
/*                           of the spline.                                 */
/*                   output: unchanged                                      */
/*      double *yp   input:  description of smooth curve                    */
/*                   output: unchanged                                      */
/*      double sigma input:  tension factor, sigma > 0.0                    */
/*      double t     input:  abcissa in [0,1] to be used for calculation    */
/*                           of corresponding (xt,yt)                       */
/*      double *tx   output: x coordinate of smoothed curve at t            */
/*      double *ty   output: y coordinate of smoothed curve at t            */
/*      int first    input:  first = 1 when the function is called first;   */
/*                           afterwards: first = 0.                         */

/* ------------------------------------------------------------------------ */
/*  splakima                                                                */
/*      Calculation of a smooth curve, addopted from:                       */
/*      H. Akima, Interpolation and smooth curve fitting based on local     */
/*      procedures (algorithm 433), Communications of the ACM 15 (1972),    */
/*      no. 10, 914 -918                                                    */
/*                                                                          */
/*      Call: int splakima (mode,l,x,y,m,n,u,v)                             */
/*                                                                          */
/*      int mode    input:   1  x coordinates must be strictly ascending    */
/*                           2  arbitrary points allowed                    */
/*      int l       input:  number of data points, l >= 2                   */
/*      double *x   input:  array x(i), i = 1,l  with x coordinates         */
/*                  output: unchanged                                       */
/*      double *y   input:  array y(i), i = 1,l  with y coordinates         */
/*                  output: unchanged                                       */
/*      int m       input:  number of subintervals, m >= 2                  */
/*      int n       input:  length of the arrays u and v with               */
/*                          n = (l - 1) * m + 1                             */
/*      double *u   input:  array u(i) of length n                          */
/*                  output: contains the x coordinates of the smooth curve  */
/*      double *v   input:  array v[i] of length n                          */
/*                  output: contains the y coordinates of the smooth curve  */
/*                                                                          */
/*      Return:     0  successful                                           */
/*                 -1  error in input parameters                            */
/*                 -2  if x coordinates not increasing                      */
/*                 -3  if identical points                                  */

int splakima(TDAContext *ctx, int mode,int l,double *x,double *y,int m,int n,double *u,double *v)
{
    (void)ctx;        /* unused: the signature is shared */
    int i,j,md0,mdm1,l0,m0,mm1,n0,lm1,k,k5;
    double a1 = 0.0,a2 = 0.0,a3,a4,m1 = 0.0,m2 = 0.0,m3,m4,sw,w2,w3;
    double t2 = 0.0,t3,x2,x3,x4,x5,y2,y3,y4,y5,sin2 = 0.0,sin3,cos2 = 0.0,cos3,p1,rm;

    x5 = y5 = t3 = sin3 = cos3 = 0.0;
    md0 = mode;
    mdm1 = md0 - 1;
    l0 = l;
    lm1 = l0 - 1;
    m0 = m;
    mm1 = m0 - 1;
    n0 = n;

    if (md0 <= 0 || md0 >= 3 || lm1 <= 0 || mm1 <= 0 || n0 != (lm1 * m0 + 1))  
        return(-1);
 
    if (md0 == 1) {
        i = 2;

        if (x[1] < x[2]) {  
            for (i = 3; i <= l0; ++i) {
                if (x[i] <= x[i - 1])  
                    return(-2);
            }
        }
        else if (x[1] > x[2]) {
            for (i = 3; i <= l0; ++i) {
                if (x[i] >= x[i - 1])  
                    return(-2);
            }
        }
        else  
            return(-2);
    }
    else {
        for (i = 2; i <= l0; ++i) {
            if (x[i - 1] == x[i] && y[i - 1] == y[i])  
                return(-3);
        }
    }
    k = n0 + m0;
    i = l0 + 1;

    for (j=1; j<=l0; ++j) {
        k -= m0;
        i--;
        u[k] = x[i];
        v[k] = y[i];  
    }
    rm = m0;
    rm = 1.0 / rm; 
    k5 = m0 + 1;

    for (i = 1; i <= l0; ++i) {

        if (i > 1)
            goto Lab40;

        x3 = u[1];
        y3 = v[1];
        x4 = u[m0 + 1];
        y4 = v[m0 + 1];
        a3 = x4 - x3;
        m3 = y4 - y3;

        if (mdm1 == 0)
            m3 /= a3;

        if ( l0 != 2) 
            goto Lab41;

        a4 = a3;
        m4 = m3;

Lab31:
        if (md0 == 2) {
            a2 = a3 + a3 - a4;
            a1 = a2 + a2 - a3;
        }
        m2 = m3 + m3 - m4;
        m1 = m2 + m2 - m3;

        if (md0 == 1)   
            goto Lab51;
        else
            goto Lab56;

Lab40:
        x2 = x3;
        y2 = y3; 
        x3 = x4;
        y3 = y4; 
        x4 = x5; 
        y4 = y5; 

        a1 = a2;
        m1 = m2; 
        a2 = a3; 
        m2 = m3; 
        a3 = a4;
        m3 = m4;

        if (i >= lm1)
            goto Lab42;

Lab41:
        k5 = k5 + m0;
        x5 = u[k5];
        y5 = v[k5]; 
        a4 = x5 - x4;
        m4 = y5 - y4;

        if (mdm1 == 0)
            m4 /= a4;
        goto Lab43;

Lab42:
        if (mdm1 != 0)
            a4 = a3 + a3 - a2;
        m4 = m3 + m3 - m2;

Lab43:
        if (i == 1) 
            goto Lab31;
        if (md0 == 1) 
            goto Lab50;
        else
            goto Lab55;

Lab50:
        t2 = t3;

Lab51:
        w2 = fabs(m4 - m3);
        w3 = fabs(m2 - m1); 
        sw = w2 + w3;

        if (sw == 0.0) {
            w2 = w3 = 0.5;
            sw = 1.0;
        }
        t3 = (w2 * m2 + w3 * m3) / sw;

        if (i > 1)
            goto Lab60;
        else
            goto Lab80;

Lab55:
        cos2 = cos3;
        sin2 = sin3;

Lab56:
        w2 = fabs(a3 * m4 - a4 * m3);
        w3 = fabs(a1 * m2 - a2 * m1);

        if ((w2 + w3) != 0.0)
            goto Lab57;

        w2 = sqrt(a3 * a3 + m3 * m3);
        w3 = sqrt(a2 * a2 + m2 * m2);

Lab57:
        cos3 = w2 * a2 + w3 * a3;
        sin3 = w2 * m2 + w3 * m3;
        sw = cos3 * cos3 + sin3 * sin3;
        if (sw == 0.0) 
            goto Lab58;

        sw = sqrt(sw);
        cos3 /= sw;
        sin3 /= sw;

Lab58:
        if (i > 1)  
            goto Lab65;
        else
            goto Lab80;

Lab60:
        w2 = (2.0 * (m2 - t2) + m2 - t3) / a2;
        w3 = (-m2 - m2 + t2 + t3) / (a2 * a2);
        goto Lab70;

Lab65:
        sw = sqrt(a2 * a2 + m2 * m2);
        p1 = sw * cos2;
        a1 = 3.0 * a2 - sw * (cos2 + cos2 + cos3);
        m1 = a2 - p1 - a1;

        t2 = sw * sin2;
        w2 = 3.0 * m2 - sw * (sin2 + sin2 + sin3);
        w3 = m2 - t2 - w2;

        goto Lab75;

Lab70:
        a2 = a2 * rm;
        sw = 0.0;

        for (j = 1; j <= mm1; ++j) {
            k++;
            sw += a2; 
            u[k] = x2 + sw;
            v[k] = y2 + sw * (t2 + sw * (w2 + sw * w3));
        }
        goto Lab79;

Lab75:
        sw = 0.0;

        for (j = 1; j <= mm1; ++j) {
            k++;
            sw += rm;
            u[k] = x2 + sw * (p1 + sw * (a1 + sw * m1));
            v[k] = y2 + sw * (t2 + sw * (w2 + sw * w3)); 
        }

Lab79:
        k++;
Lab80: ;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  splf    Calculation of a smoothing spline based on:                     */
/*  ##                                                                      */
/*      P. Dierckx, An Algorithm for Smoothing, Differentiation and         */
/*      Integration of Experimental Data Using Spline Functions, in:        */
/*      Journal of Computational and Applied Mathematics, Vol. 1 (1975),    */
/*      No. 3, 165 - 184.                                                   */
/*                                                                          */
/*      Gegeben seien m Datenpunkte mit Abszissen x(i), Ordinaten y(i)      */
/*      und Gewichten w(i), i = 1,m. Beachte, dass m > 2 gelten muss und    */
/*      die Abszissenwerte eine streng aufsteigende Folge bilden muessen.   */
/*      Ausserdem sei ein Glaettungsparameter S >= 0.0 gegeben.             */
/*                                                                          */
/*      Es wird dann eine Splinefunktion s(x) des Grades k berechnet; mit   */
/*      Knoten t(j) (j = 1,n). Die berechnete Splinefunktion erfuellt die   */
/*      Bedingung, dass es sich um die "glatteste" Funktion handelt, die    */
/*      die Bedingung erfuellt:                                             */
/*                                                                          */
/*          SUM (i = 1,m) wi * (yi - s(xi)) ** 2 <= S                       */
/*                                                                          */
/*      Beachte: wenn S == 0.0, wird eine interpolierende Splinefunktion    */
/*      berechnet (Return-Code 1). Wenn S >> 0.0, wird naeherungsweise      */
/*      eine OLS-Approximation mit einem Polynom k.ten Grades berechnet     */
/*      (Return-Code 2).                                                    */
/*                                                                          */
/*      Der Output dieser Funktion besteht in den Knoten der Splinefunktion */
/*      sowie in den Koeffizienten c(j) der normalisierten B-Spline-        */
/*      Repraesentation der Splinefunktion.                                 */
/*                                                                          */
/*      Zur Weiterverarbeitung dieses Outputs stehen die Funktionen:        */
/*                                                                          */
/*          spldf   :   Berechnung von Funktionswerten und Ableitungen      */
/*          splif   :   Berechnung von Integralen                           */
/*                                                                          */
/*      zur Verfuegung.                                                     */
/*                                                                          */
/*      Beachte: damit die Splinefunktion berechnet werden kann, muss ein   */
/*      Schaetzintervall [xa,xb] angegeben werden, wobei gelten muss:       */
/*                                                                          */
/*          xa <= x[1] < x[m] <= xb                                         */
/*                                                                          */
/*      Beachte: Da die Anzahl der Knoten nicht vorab bestimmt werden kann, */
/*      muss die Funktion mit einem Schaetzwert fuer die Anzahl der Knoten  */
/*      aufgerufen werden. Fuer die Bestimmung dieses Schaetzwerts (nn)     */
/*      kann man sich an folgende Regeln halten:                            */
/*                                                                          */
/*      (1) 3 * k <= nn <= m + k + 1                                        */
/*      (2) je kleiner der Wert von S, desto groesser muss nn sein          */
/*      (3) normalerweise ist nn = m / 2 ein ausreichender Schaetzwert      */
/*                                                                          */
/*                                                                          */
/*      Call:  int splf(m,x,y,w,xa,xb,k,s,nn,t,c,maxit)                     */
/*                                                                          */
/*      int m         input:  number of data points, m >= 2 * k             */
/*      double *x     input:  array x(i), i = 1,m with x coordinates        */
/*                    output: not changed                                   */
/*      double *y     input:  array y(i), i = 1,m with y coordinates        */
/*                    output: not changed                                   */
/*      double *w     input:  array w(i), i = 1,m with weights, > 0.0       */
/*                    output: unchanged                                     */
/*      double xa,xb  input:  end points of the interval for approximation  */
/*                            with: xa <= x[1] < x[m] <= xb                 */
/*      int k         input:  degree of spline: k >= 2                      */
/*      double s      input:  smoothing factor s >= 0.0                     */
/*      int *nn       input:  upper bound for number of knots.              */
/*                    output: actual number of knots.                       */
/*      double *t     input:  array t(i) (i = 1,nn)                         */
/*                    output: knots of the spline function                  */
/*      double *c     input:  array c(i) (i = 1,nn - k - 1)                 */
/*                    output: coefficients of the B-spline representation.  */
/*      int maxit     input:  max number of iterations.                     */
/*                                                                          */
/*      Return:     0  successful calculation of smoothing spline           */
/*                  1  successful calculation of interpolating spline       */
/*                  2  successful OLS approximation                         */
/*                 -1  error in input parameters                            */
/*                 -2  if nn is not large enough                            */
/*                 -3  if exceeded max number of iterations                 */
/*                 -4  unexpected error (ier = 3)                           */
/*                 -5  unexpected error (ier = 4)                           */
/*                 -6  abscissae not strictly increasing                    */
/*                 -7  if weights are not strictly positive                 */
/*                 -8  if insufficient memory.                              */

static int ad,bd,gd;       /* column dimension of a, b, g */

int splf(TDAContext *ctx, int m,double *x,double *y,double *w,double xa,double xb, int k,double s,int *nn,double * t,double *c,int maxit)
{
    register int i,j,jj;
    int i1,ik,i2,i3,ir,ir1,ir2 = 0,in,it,iter,ifib1,ifib2,ifib3,j1,ni,ne,rflg;
    int k1,k5,k6,l,li,lj,lk,ll,ls,l1,l2,l3,l4,m1,mi,me,m2,n,nnn,nmin,nmax,nk1;
    double wmax,s2,f,f1,f2,fp,fp2,arg,p,pinv,dfp;
    double *a,*b,*g,*v,*z,*h,*h1,*h2;
    int a_a,b_a,g_a,v_a,z_a,h_a,h1_a,h2_a;

    double tol = 0.01;          /* tolerance for newton iteration   */

    a_a = b_a = g_a = v_a = z_a = h_a = h1_a = h2_a = 0;

    rflg = -1;
    p = pinv = 0.0;

    if (k < 2 || m < 2 * k || s < 0.0)  
        return(-1);

    if (xa > x[1] || x[m] > xb)
        return(-1);

    if ((wmax = w[1]) <= 0.0)
        return(-7);

    for (i = 2; i <= m; ++i) {
        if (x[i - 1] >= x[i])
            return(-6);
        if (w[i] <= 0.0)
            return(-7);
        if (w[i] > wmax)
            wmax = w[i];
    }
    s2 = sqrt(s);
    ad = k1 = k + 1;
    bd = gd = k + 2;
    nnn = *nn;

    rflg = -8;
    if (!(a  = (double *) calloc ((size_t)((nnn - k1) * ad + 1),sizeof(double))))   
        goto Splf_Fin;
    a_a = (nnn - k1) * ad + 1;
    memrq(ctx, a_a,sizeof(double));

    if (!(b  = (double *) calloc ((size_t)((nnn - k1) * bd + 1),sizeof(double))))   
        goto Splf_Fin;
    b_a = (nnn - k1) * bd + 1;
    memrq(ctx, b_a,sizeof(double));

    if (!(g  = (double *) calloc ((size_t)((nnn - k1) * gd + 1),sizeof(double))))   
        goto Splf_Fin;
    g_a = (nnn - k1) * gd + 1;
    memrq(ctx, g_a,sizeof(double));

    if (!(v  = (double *) calloc ((size_t)(nnn - k),sizeof(double))))   
        goto Splf_Fin;
    v_a = nnn - k;
    memrq(ctx, v_a,sizeof(double));

    if (!(z  = (double *) calloc ((size_t)(nnn - k),sizeof(double))))   
        goto Splf_Fin;
    z_a = nnn - k;
    memrq(ctx, z_a,sizeof(double));

    if (!(h  = (double *) calloc ((size_t)(2 * k + 3),sizeof(double))))   
        goto Splf_Fin;
    h_a = 2 * k + 3;
    memrq(ctx, h_a,sizeof(double));

    if (!(h1 = (double *) calloc ((size_t)(k + 2),sizeof(double))))   
        goto Splf_Fin;
    h1_a = k + 2;
    memrq(ctx, h1_a,sizeof(double));

    if (!(h2 = (double *) calloc ((size_t)(k + 2),sizeof(double))))   
        goto Splf_Fin;
    h2_a = k + 2;
    memrq(ctx, h2_a,sizeof(double));
    rflg = 0;

    m1 = m - 1;
    ifib1 = 1;
    ifib2 = 2;
    nmin = 3 * k + 1;
    nmax = m + k1;
    n = nmin;
    if (s == 0.0)
        n = nmax;
    else
        n = nmin;

Lab100:
    *nn = n;

    if (n > nnn) {      /* storage space for number of knots insufficient */
        rflg = -2;
        goto Splf_Fin;
    }
    nk1 = n - k1;
    m2 = n - 2 * k1 + 1;
    l = m1 / m2 + 1;
    ir = m1 - m2 * (l - 1);
    mi = l + 1;
    me = m - l;
    ni = k1 + 1;
    ne = nk1;

    for (jj = 1; jj <= 2; ++jj) {
        ir1 = ir / 2;

        for (j = 1; j <= ir1; ++j) {
             t[ni] = x[mi];
             t[ne] = x[me];
             ni++;
             mi += l;
             ne--;
             me -= l;
        }
        if (jj == 2) 
            break;

        ir2 = 2 * ir1 - ir + 1;
        l--;
        ir = m2 - ir - ir2;
        mi--;
        me++;
    }
    if (ir1 * 2 != ir) {
        t[ni] = x[mi];
        if (ir2 != 1)
            t[ni] = (x[mi] + x[mi + 1]) * 0.5;
    }
    t[k1] = xa;
    t[nk1 + 1] = xb;
    f1 = t[k1 + 1] - xa;
    f2 = xb - t[nk1];

    for (j = 1; j <= k; ++j) {
        i = k1 - j;
        in = n - i;
        t[i] = t[i + 1] - f1;
        t[in + 1] = t[in] + f2;
    }
    for (ir = 1; ir <= nk1; ++ir) {
        z[ir] = 0.0;
        for (ik = 1; ik <= k1; ++ik)
            a[(ir - 1) * ad + ik] = 0.0;
    }
    l = k1;

    for (it = 1; it <= m; ++it) {
        arg = x[it];
        if (arg >= t[l + 1] && l != nk1)
            l++;

        h1[1] = 1.0;
        for (j = 1; j <= k; ++j) {
            h2[1] = 0.0;
            j1 = j + 1;
            for (i = 1; i <= j; ++i) {
                li = l + i;
                lj = li - j;
                f1 = h1[i] / (t[li] - t[lj]);
                h2[i] += f1 * (t[li] - arg);
                h2[i + 1] = f1 * (arg - t[lj]);
            }
            for (i = 1; i <= j1; ++i)  
                h1[i] = h2[i];
        }
        lk = l - k;

        for (l1 = lk; l1 <= l; ++l1) {
            k5 = l1 - l + k1;
            z[l1] += h1[k5] * y[it] * w[it];

            for (l2 = l1; l2 <= l; ++l2) {
                k6 = l2 - l + k1;
                ir = l2;
                ik = k1 - l2 + l1;
                a[(ir - 1) * ad + ik] += h1[k5] * h1[k6] * w[it];
            }
        }
    }
    iter = 0;

    for (ir = 1; ir <= nk1; ++ir) {
        for (ik = 1; ik <= k1; ++ik)  
            g[(ir - 1) * gd + ik] = a[(ir - 1) * ad + ik];
    }
    ls = k1;

    while (1) {     /* Label 400 */

        if (bandet(ctx, g,ls,nk1)) { 
            if (iter == 0) {
                rflg = -5;
                goto Splf_Fin;
            }
            p *= 10.0;
            goto Lab700;
        }
        bansol(ctx, g,ls,nk1,z,c);
     
        if (s == 0.0) {
            rflg = 1;
            goto Splf_Fin;           /* interpolating spline */
        }
        fp = 0.0;
        l = k1;
     
        for (it = 1; it <= m; ++it) {
            arg = x[it];
            if (arg >= t[l + 1] && l != nk1)
                l++;
     
            for (i = 1; i <= k1; ++i)  
                h[i] = c[l + i - k1];
     
            for (j = 2; j <= k1; ++j) {
                for (jj = j; jj <= k1; ++jj) {
                    i = j + k1 - jj;
                    li = l + i - k1;
                    lj = l + i - j + 1;
     
                    h[i] = ((arg - t[li]) * h[i] + (t[lj] - arg) *
                                                h[i - 1]) / (t[lj] - t[li]);
                }
            }
            fp += w[it] * (y[it] - h[k1]) * (y[it] - h[k1]);
        }
        if (fabs((fp - s) / s) < tol) {     /* check convergence */
            rflg = 0;
            goto Splf_Fin;
        }
        if (iter != 0) {
            if (s > fp) {
                if (iter < 2)       /* OLS approximation */
                    rflg = 2;
                else
                    rflg = -4;
                goto Splf_Fin;
            }
            if (iter >= maxit) {      /* max number of iterations */
                rflg = -3;
                goto Splf_Fin;
            }
            for (i = 1; i <= nk1; ++i) {
                if ((l1 = i - k1) < 1)
                    l1 = 1;
     
                f = 0.0;
                for (i2 = l1; i2 <= i; ++i2) {
                    i3 = i2 - i + ls;
                    f += c[i2] * b[(i - 1) * bd + i3];
                }
                if (i != nk1) {
                    if ((l1 = i + k1) > nk1)
                        l1 = nk1;
                    i1 = i + 1;
                    for (i2 = i1; i2 <= l1; ++i2) {
                        i3 = i + ls - i2;
                        f += c[i2] * b[(i2 - 1) * bd + i3];
                    }
                }
                v[i] = f;
            }
            bansol(ctx, g,ls,nk1,v,c);
            dfp = 0.0;
     
            for (i = 1; i <= nk1; ++i)  
                dfp += c[i] * v[i];

            dfp = -2.0 * dfp * pinv * pinv * pinv;
            fp2 = sqrt(fp);
            p  += 2.0 * fp2 / s2 * (s2 * fp2 - fp) / dfp;
            iter++;
        }
        else {
            if (fp > s) 
                break;
            ls++;        
     
            for (ir = 1; ir <= nk1; ++ir) {
                for (ik = 1; ik <= ls; ++ik)
                    b[(ir - 1) * bd + ik] = 0.0;
            }
            for (l = ls; l <= nk1; ++l) {
                for (j = 1; j <= k1; ++j) {
                    l1 = l + j;
                    l2 = l1 - ls;
                    k5 = k1 + j;
                    h[j] = t[l] - t[l2];
                    h[k5] = t[l] - t[l1];
                }
                f1 = -h[ls] * h[k1];
     
                for (j = 1; j <= ls; ++j) {
                    for (i = j; i <= ls; ++i) {
                        f = 1.0;
                        for (ll = 1; ll <= k1; ++ll) {
                            l1 = j + ll - 1;
                            l2 = i + ll - 1;
                            f2 = h[l1] * h[l2] / f1;
                            f  = f * f2;
                        }
                        ir = i - 1 + l - k1;
                        ik = ls - i + j;
                        l1 = l + j - 1;
                        l2 = l + i - 1;
                        l3 = l1 - k1;
                        l4 = l2 - k1;
                        b[(ir - 1) * bd + ik] += 
                                (t[l1] - t[l3]) * (t[l2] - t[l4]) / (f * f1);
                    }
                }
            }
            iter = 1;
            p = 0.0001 / wmax;
        }
Lab700:
        pinv = 1.0 / p;
     
        for (ir = 1; ir <= nk1; ++ir) {
            g[(ir - 1) * gd + 1] = pinv * b[(ir - 1) * bd + 1];
            for (ik = 2; ik <= ls; ++ik) {
                g[(ir - 1) * gd + ik] = a[(ir - 1) * ad + ik - 1] + pinv 
                                * b[(ir - 1) * bd + ik];
            }
        }
    }       /* end of while: Label 400 */

    if (n == nmax) {
        rflg = 1;           /* interpolating spline */
        goto Splf_Fin;
    }
    ifib3 = ifib1 + ifib2;
    ifib1 = ifib2;
    ifib2 = ifib3;
    n = nmin + ifib3;
    if (n > nmax) 
        n = nmax;

    goto Lab100;

Splf_Fin:
    if (a_a > 0) {
        free((char *)a);
        memrq(ctx, -a_a,sizeof(double));
    }
    if (b_a > 0) {
        free((char *)b);
        memrq(ctx, -b_a,sizeof(double));
    }
    if (g_a > 0) {
        free((char *)g);
        memrq(ctx, -g_a,sizeof(double));
    }
    if (v_a > 0) {
        free((char *)v);
        memrq(ctx, -v_a,sizeof(double));
    }
    if (z_a > 0) {
        free((char *)z);
        memrq(ctx, -z_a,sizeof(double));
    }
    if (h_a > 0) {
        free((char *)h);
        memrq(ctx, -h_a,sizeof(double));
    }
    if (h1_a > 0) {
        free((char *)h1);
        memrq(ctx, -h1_a,sizeof(double));
    }
    if (h2_a > 0) {
        free((char *)h2);
        memrq(ctx, -h2_a,sizeof(double));
    }
    return(rflg);
}

/* ------------------------------------------------------------------------ */

int bandet(TDAContext *ctx, double *g, int ls, int nk1)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j,k;
    int p,q,r,s;
    double y = 0.0;

    for (i = 1; i <= nk1; ++i) {
        p = 1;
        if (i <= ls)
            p = ls - i + 1;

        r = i - ls + p;

        for (j = p; j <= ls; ++j) {
            s = j - 1;
            q = ls - j + p;
            y = g[(i - 1) * gd + j];

            for (k = p; k <= s; ++k) {
                y -= g[(i - 1) * gd + k] * g[(r - 1) * gd + q];
                q++;
            }
            if (j == ls)
                break;
            g[(i - 1) * gd + j] = y * g[(r - 1) * gd + ls];
            r++;
        }
        if (y <= 0.0)
            return(-1);     /* g not positive definite */

        g[(i - 1) * gd + j] = 1.0 / sqrt(y);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */

void bansol(TDAContext *ctx, double *g, int ls, int nk1, double *z, double *c)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j;
    int p,q,r,l,k;
    double y;

    l = ls - 1;

    for (i = 1; i <= nk1; ++i) {
        y = z[i];
        if (i != 1) {
            p = 1;
            if (i <= ls)
                p = ls - i + 1;

            q = i;

            for (j = p; j <= l; ++j) {
                k = p + l - j;
                q--;
                y -= g[(i - 1) * gd + k] * c[q];
            }
        }
        c[i] = y * g[(i - 1) * gd + ls];
    }
    for (i = 1; i <= nk1; ++i) {
        r = nk1 + 1 - i;
        y = c[r];
        if (r != nk1) {
            p = 1;
            if (nk1 - r < ls)
                p = ls - nk1 + r;

            q = r;
            for (j = p; j <= l; ++j) {
                k = p + l - j;
                q++;
                y -= g[(q - 1) * gd + k] * c[q];
            }
        }
        c[r] = y * g[(r - 1) * gd + ls];
    }
}

/* ------------------------------------------------------------------------ */
/*  spld                                                                    */
/*                                                                          */
/*      Berechnung der Werte oder der Ableitungen einer Splinefunktion,     */
/*      die zuvor mit der Funktion splf() berechnet worden ist. Zur         */
/*      naeheren Beschreibung vgl. splf.                                    */
/*                                                                          */
/*      Es sei s(x) eine Splinefunktion vom Grad k, die mit splf()          */
/*      berechnet worden ist. Als Output von splf() stehen zur Verfuegung   */
/*      die Folge der Knoten t[j] (j = 1,n) und die Koeffizienten der       */
/*      B-Spline-Darstellung der Splinefunktion c[j] (j = 1, n - k - 1).    */
/*      Diese Funktion berechnet dann den Wert der nu.ten Ableitung der     */
/*      Splinefunktion an der Stelle arg. Dabei muss gelten: nu >= 0.       */
/*      Wenn nu = 0, wird der Wert der Splinefunktion an der Stelle arg     */
/*      berechnet.                                                          */
/*                                                                          */
/*      Als Input muss weiterhin zur Verfuegung gestellt werden ein Index l */
/*      so dass:                                                            */
/*                  t[l] <= x < t[l + 1]                                    */
/*                                                                          */
/*      oder l = n - k - 1, wenn arg = t[n - k]                             */
/*                                                                          */
/*      Call: double spld(n,t,k,c,nu,arg,l,tmp)                             */
/*                                                                          */
/*      int n        input:  Anzahl der Knoten                              */
/*      double *t    input:  Feld mit den Knotenwerten t(i), i = 1,n        */
/*      int k        input:  Grad der Splinefunktion                        */
/*      double *c    input:  Feld der Laenge: 1,n - k - 1; enthaelt die     */
/*                           Koeffizienten der B-Spline-Darstellung         */
/*      int nu       input:  Ordnung der Ableitung, nu >= 0                 */
/*      double arg   input:  das Argument, bei dem der Wert der Spline-     */
/*                           funktion berechnet werden soll                 */
/*      int l        input:  Parameter, der die Position des Arguments in   */
/*                           der Menge der Knoten bestimmt; es muss gelten: */
/*                           t[l] <= arg < t[l + 1] oder l = n - k - 1,     */
/*                           wenn arg = t[n - k]                            */
/*      double *tmp  input:  work array of length >= k + 2                  */
/*                                                                          */
/*      Return:      Wert der Splinefunktion an der Stelle arg              */
/*                                                                          */

double spld(TDAContext *ctx, int n, double *t,int k,double *c,int nu,double arg,int l,double *tmp) 
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j,jj;
    int li,lj,ik,k1,nu1,nu2;
    double deriv;

    i = n;
    deriv = 0.0;
    k1 = k + 1;

    if (nu < 0 || nu >= k1)
        return(deriv);

    for (i = 1; i <= k1; ++i) {
        ik = l + i - k1;
        tmp[i] = c[ik];
    }
    if (nu > 0) {
        nu1 = nu + 1;
        for (j = 2; j <= nu1; ++j) {
            for (jj = j; jj <= k1; ++jj) {
                i = j + k1 - jj;
                li = l + i - k1;
                lj = l + i - j + 1;
                tmp[i] = (tmp[i] - tmp[i - 1]) / (t[lj] - t[li]);
            }
        }
    }
    if (nu == 0 || nu < k1 - 1) {
        nu2 = nu + 2;

        for (j = nu2; j <= k1; ++j) {
            for (jj = j; jj <= k1; ++jj) {
                i = j + k1 - jj;
                li = l + i - k1;
                lj = l + i - j + 1;
                tmp[i] = ((arg-t[li]) * tmp[i] + 
                                   (t[lj]-arg) * tmp[i-1])/(t[lj]-t[li]);
            }
        }
    }
    deriv = tmp[k1];
    if (nu) {
        for (i = 1; i <= nu; ++i) 
            deriv *= (double)(k1 - i);
    }
    return (deriv);
}

/* ------------------------------------------------------------------------ */
/*  interv  Used to calculate the place of an entry into an ascending       */
/*          series.                                                         */
/*                                                                          */
/*      Ref. C. de Boor, A Practical Guide to Splines, NY-Berlin 1978, 92   */
/*                                                                          */
/*      int interv(n,xt,x,init)                                             */
/*                                                                          */
/*      int n        input:  number of values in the series x[].            */
/*      double *x    input:  array x[i] (i = 1,n), must be ascending.       */
/*      double arg   input:  an entry in the series.                        */
/*      int init     input:  1 for initialization                           */
/*                                                                          */
/*      Return:      0  if           arg < x[1]                             */
/*                   i  if   x[i] <= arg < x[i + 1] (i = 1,n - 1)           */
/*                   n  if   x[n] <= arg                                    */
/*                                                                          */

int interv(TDAContext *ctx, int n,double *x,double arg, int init)
{
    int ihi,middle,istep;

    if (init) {
        ctx->s_interv_ilo = 1;
        return(0);
    }
    ihi = ctx->s_interv_ilo + 1;

    if (ihi >= n) {
        if (arg >= x[n])
            return(n);
        if (n <= 1)
            return(0);
        ctx->s_interv_ilo = n - 1;
        ihi = n;
    }
    if (arg >= x[ihi]) {
        istep = 1;
        while(1) {
            ctx->s_interv_ilo = ihi; ihi = ctx->s_interv_ilo + istep;
            if (ihi >= n)
                break;
            if (arg < x[ihi])
                goto IV_Fin;
            istep *= 2;
        }
        if (arg >= x[n])
            return(n);
        ihi = n;
    }
    else {
        if (arg >= x[ctx->s_interv_ilo]) 
            return(ctx->s_interv_ilo);
        istep = 1;
        while(1) {
            ihi = ctx->s_interv_ilo;
            ctx->s_interv_ilo = ihi - istep;
            if (ctx->s_interv_ilo <= 1)
                break;
            if (arg >= x[ctx->s_interv_ilo])
                goto IV_Fin;
            istep *= 2;
        }
        ctx->s_interv_ilo = 1;
        if (arg < x[1])
            return(0);
    }
IV_Fin:
    while(1) {
        if ((middle = (ctx->s_interv_ilo + ihi) / 2) == ctx->s_interv_ilo)
            return(ctx->s_interv_ilo);
        if (arg < x[middle])
            ihi = middle;
        else 
            ctx->s_interv_ilo = middle;
    }
}



void tda_reset_t_spl(void)
{
    ad = bd = gd = 0;
}
