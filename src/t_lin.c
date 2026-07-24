/****************************************************************************/
/*  t_lin                                                                   */
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
#include "t_gf.h"

/*  functions in t_lin.c */

int decomp(int n, double *x);
void dclear(int n, double *da);
void dcopy(int n,double *dx,int incx, double *dy,int incy);
double ddot(int n,double *dx,int incx,double *dy,int incy);
void dscal(int n,double da,double *dx,int incx);
double fsign(double a,double b);
void daxpy(int n,double da,double *dx,int incx,double *dy,int incy);
double twonrm(int n,double *v); 
int dpodi(double *a,int ma,int n,double *det,int job);
int dpofa(double *a,int ma,int n);
void drot(int n,double *dx,int incx,double *dy,int incy,double dc,double ds);
void drotg(double *da,double *db,double *dc,double *ds);
void dxpy(int n,int m,double *x,double *y,double *xpy);
void dzero(int n,int m,double *a,int cda);
double dnrm2(int n,double *dx,int incx);
void dqrdc(double *x,int cdx,int n,int p,double *qraux,int *jpvt,      
    double *work,int job);
void dqrsl(double *x,int cdx,int n,int k,double *qraux,double *y,double *qy,
    double *qty,double *b,double *rsd,double *xb,int job,int *info);
void dtrco(double *t,int cdt,int n,double *rcond,double *z,int job);
int idamax(int n,double *dx,int incx);
void dswap(int n,double *dx,int incx,double *dy,int incy);
double dasum(int n,double *dx,int incx);
void dtrsl(double *t,int cdt,int n,double *b,int job,int *info);                                                

/* ------------------------------------------------------------------------ */
/*  decomp(n,x)                                                             */
/*      Cholesky decomposition of a (n,n) matrix x. It is assumed that x    */
/*      is given as an lower triangular matrix with n (n + 1) / 2 elements. */
/*      x is substituted by its lower Cholesky factor.                      */
/*      Return 0 if positive definite, else 1.                              */

int decomp(int n, double *x)
{
    register int i,j,k,l,m;
    double tmp,tmp1;

    for (i = 1; i <= n; ++i) {
        m = i * (i - 1) / 2;    
        tmp1 = x[m + i];
        for (j = 1; j < i; ++j) {
            l = j * (j - 1) / 2;
            tmp = x[m + j];
            for (k = 1; k < j; ++k)  
                tmp -= x[m + k] * x[++l];
            tmp /= x[j * (j + 1) / 2];
            x[m + j] = tmp;
            tmp1 -= tmp * tmp;
        }
        if (tmp1 < EPSI1)   
            return(i);
        x[m + i] = sqrt(tmp1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  dclear.   Sets a vector to zero.                                        */

void dclear(int n, double *da)
{
    register int i;
    for (i = 1; i <= n; ++i)  
        da[i] = 0.0;
}

/* ------------------------------------------------------------------------ */
/*  dcopy(n,dx,incx,dy,incy)                                                */
/*  qq                                                                      */
/*  COPIES A VECTOR, X, TO A VECTOR, Y.                                     */
/*  USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.                        */
/*  JACK DONGARRA, LINPACK, 3/11/78.                                        */
/*                                                                          */
 
void dcopy(int n,double *dx,int incx, double *dy,int incy)
{
    register int i,ix,iy,m;

    if (n <= 0)         
        return;

    if (incx == 1 && incy == 1)             
        goto L20;
 
    /* CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
    /* NOT EQUAL TO 1 */
 
    ix = iy = 1;
    if (incx < 0)
        ix = (-n + 1) * incx + 1;
    if (incy < 0)
        iy = (-n + 1) * incy + 1;

    for (i = 1; i <= n; ++i) {
        dy[iy] = dx[ix];
        ix = ix + incx;
        iy = iy + incy;
    }
    return;   
 
    /* CODE FOR BOTH INCREMENTS EQUAL TO 1 */
    /* CLEAN-UP LOOP */
 
L20:
    m = n % 7;   
    if (m == 0)              
        goto L40;
    for (i = 1; i <= m; ++i) 
        dy[i] = dx[i];
    if (n <  7)             
        return;
L40:
    for (i = m + 1; i <= n; i += 7) {
        dy[i] = dx[i];
        dy[i + 1] = dx[i + 1];
        dy[i + 2] = dx[i + 2];
        dy[i + 3] = dx[i + 3];
        dy[i + 4] = dx[i + 4];
        dy[i + 5] = dx[i + 5];
        dy[i + 6] = dx[i + 6];
    }
}
     
/* ------------------------------------------------------------------------ */
/*  ddot(n,dx,incx,dy,incy)                                                 */
/*                                                                          */
/*     forms the dot product of two vectors.                                */
/*     uses unrolled loops for increments equal to one.                     */
/*     jack dongarra, linpack, 3/11/78.                                     */
/*                                                                          */ 

double ddot(int n,double *dx,int incx,double *dy,int incy)
{
    register int i,ix,iy;
    int m;
    double ddot;

    ddot = 0.0;
    if (n <= 0)
        return(ddot);

    if (incx == 1 && incy == 1)           
        goto L20;

    /* code for unequal increments or equal increments not equal to 1 */
    
    ix = iy = 1;
    if (incx < 0)
        ix = (-n + 1) * incx + 1;
    if (incy < 0)
        iy = (-n + 1) * incy + 1;

    for (i = 1; i <= n; ++i) {
        ddot += dx[ix] * dy[iy];
        ix += incx;
        iy += incy;
    }
    return(ddot);
   
    /* code for both increments equal to 1 */ 
     
L20:
    m = n % 5;    
    if (m != 0) {
        for (i = 1; i <= m; ++i)
            ddot += dx[i] * dy[i];

        if (n < 5)
            return(ddot);
    }
    for (i = m + 1; i <= n; i += 5) {
        ddot += dx[i] * dy[i] + dx[i + 1] * dy[i + 1] +
                dx[i + 2] * dy[i + 2] + dx[i + 3] * dy[i + 3] +
                dx[i + 4] * dy[i + 4];
    }
    return(ddot);
}             
   
/* ------------------------------------------------------------------------ */
/*  dscal(n,da,dx,incx)                                                     */
/*                                                                          */
/*    scales a vector by a constant.                                        */
/*    uses unrolled loops for increment equal to one.                       */
/*    jack dongarra, linpack, 3/11/78.                                      */
/*                                                                          */

void dscal(int n,double da,double *dx,int incx)
{
    register int i;
    int nincx,m;
 
    if (n <= 0)
        return;
    if (incx == 1)
        goto L20;
 
    /* code for increment not equal to 1 */
 
    nincx = n * incx;
    for (i = 1; i <= nincx; i += incx)
        dx[i] *= da;         
    return;
 
    /* code for increment equal to 1 */
 
L20:
    m = n % 5;
    if (m != 0) {
        for (i = 1; i <= m; ++i) 
            dx[i] *= da;         
        if (n < 5)
            return;
    } 
    for (i = m + 1; i <= n; i += 5) {
        dx[i] *= da;         
        dx[i + 1] *= da;         
        dx[i + 2] *= da;           
        dx[i + 3] *= da;         
        dx[i + 4] *= da;         
    }
}           

/* ------------------------------------------------------------------------ */
/*  fsign(a,b)                                                              */

double fsign(double a,double b)
{
    double tmp;

    if (b >= 0.0)
        tmp = 1.0;
    else
        tmp = -1.0;

    tmp *= fabs(a);
    return(tmp);
}

/* ------------------------------------------------------------------------ */
/*  daxpy(n,da,dx,incx,dy,incy)                                             */
/*                                                                          */
/*    constant times a vector plus a vector.                                */
/*    uses unrolled loops for increments equal to one.                      */
/*    jack dongarra, linpack, 3/11/78.                                      */
/*                                                                          */
 
void daxpy(int n,double da,double *dx,int incx,double *dy,int incy)
{
    register int i,ix,iy;
    int m;

    if (n <= 0 || da == 0.0)
        return;

    if (incx == 1 && incy == 1)           
        goto L20;
 
    /* code for unequal increments or equal increments not equal to 1 */
                           
    ix = iy = 1;
    if (incx < 0)
        ix = (-n + 1) * incx + 1;
    if (incy < 0)
        iy = (-n + 1) * incy + 1;

    for (i = 1; i <= n; ++i) {
        dy[iy] += da * dx[ix];
        ix += incx;
        iy += incy;
    }
    return;
 
    /* code for both increments equal to 1 */
 
L20:
    m = n % 4;

    if (m != 0) {           
        for (i = 1; i <= m; ++i)  
            dy[i] += da * dx[i];
        if(n < 4)
            return;
    } 
    for (i = m + 1; i <= n; i += 4) {
        dy[i] += da * dx[i];
        dy[i + 1] += da * dx[i + 1];
        dy[i + 2] += da * dx[i + 2];
        dy[i + 3] += da * dx[i + 3];
    }
}         

/* ------------------------------------------------------------------------ */
/*  twonrm()                                                                */
/*  compute L2 norm.                                                        */

double twonrm(int n,double *v)
{
    double temp;

    temp = ddot(n,v,1,v,1);
    return(sqrt(temp));
}             

/* ------------------------------------------------------------------------ */
/*  dpodi(a,ma,n,det,job)                                                   */
/*                                                                          */
/*    dpodi computes the determinant and inverse of a certain               */
/*    double precision symmetric positive definite matrix (see below)       */
/*    using the factors computed by dpoco, dpofa or dqrdc.                  */
/*                                                                          */
/*     on entry                                                             */
/*                                                                          */
/*        a       double precision(n,ma)                                    */
/*                the output  a  from dpoco or dpofa                        */
/*                or the output  x  from dqrdc.                             */
/*                                                                          */
/*        ma      number of columns in a                                    */
/*                                                                          */
/*        n       number of rows in a                                       */
/*                                                                          */
/*        job     integer                                                   */
/*                = 11   both determinant and inverse.                      */
/*                = 01   inverse only.                                      */
/*                = 10   determinant only.                                  */
/*                                                                          */
/*     on return                                                            */
/*                                                                          */
/*        a       if dpoco or dpofa was used to factor  a  then             */
/*                dpodi produces the upper half of inverse(a) .             */
/*                if dqrdc was used to decompose  x  then                   */
/*                dpodi produces the upper half of inverse(trans(x)*x)      */
/*                where trans(x) is the transpose.                          */
/*                elements of  a  below the diagonal are unchanged.         */
/*                if the units digit of job is zero,  a  is unchanged.      */
/*                                                                          */
/*        det     double precision(2)                                       */
/*                determinant of  a  or of  trans(x)*x  if requested.       */
/*                otherwise not referenced.                                 */
/*                determinant = det(1) * 10.0**det(2)                       */
/*                with  1.0 .le. det(1) .lt. 10.0                           */
/*                or  det(1) .eq. 0.0 .                                     */
/*                                                                          */
/*     error condition                                                      */
/*                                                                          */
/*        a division by zero will occur if the input factor contains        */
/*        a zero on the diagonal and the inverse is requested.              */
/*        it will not occur if the subroutines are called correctly         */
/*        and if dpoco or dpofa has set info .eq. 0 .                       */
/*                                                                          */
/*     Return 0 if OK, -1 if error.                                         */
/*                                                                          */
/*     linpack.  this version dated 08/14/78 .                              */
/*     cleve moler, university of new mexico, argonne national lab.         */
/*                                                                          */
  
int dpodi(double *a,int ma,int n,double *det,int job)
{
    register int i,j,k,kk;
    double s,t;

    if (job / 10 == 0)             
        goto L70;

    det[1] = 1.0;  
    det[2] = 0.0; 
    s = 10.0;  

    for (i = 1; i <= n; ++i) {
        det[1] *= pow(a[(i - 1) * ma + i],2.0);
        if (det[1] == 0.0)              
            break;     
       
        while (1) {
            if (det[1] >= 1.0)            
                break;       
            det[1] *= s;             
            det[2] -= 1.0;   
        }
        while (1) {
            if (det[1] < s)             
                break;
            if (s == 0.0) 
                return(-1);
            det[1] /= s;
            det[2] += 1.0;   
        }
    }
                 
L70:    /* compute inverse(r) */
 
    if (job % 10 == 0)              
        return(0); 

    for (k = 1; k <= n; ++k) {
        kk = (k - 1) * ma + k;
        if (a[kk] == 0.0)
            return(-1);
        a[kk] = 1.0 / a[kk];
        t = -a[kk];

        dscal(k - 1,t,a + k - 1,ma);

        if (k + 1 <= n) {        
            for (j = k + 1; j <= n; ++j) {
                t = a[(k - 1) * ma + j];
                a[(k - 1) * ma + j] = 0.0;   
                daxpy(k,t,a + k - 1,ma,a + j - 1,ma);
            }
        }
    }
 
    /* form  inverse(r) * trans(inverse(r)) */
  
    for (j = 1; j <= n; ++j) {
        if (j > 1) {             
            for (k = 1; k < j; ++k) {
                t = a[(k - 1) * ma + j];
                daxpy(k,t,a + j - 1,ma,a + k - 1,ma);
            }
        }
        t = a[(j - 1) * ma + j];
        dscal(j,t,a + j - 1,ma);
    }
    return(0);
}           
  
/* ------------------------------------------------------------------------ */
/*  info = dpofa(a,ma,n)                                                    */
/*                                                                          */
/*     dpofa factors a double precision symmetric positive definite         */
/*     matrix.                                                              */
/*                                                                          */
/*     dpofa is usually called by dpoco, but it can be called               */
/*     directly with a saving in time if  rcond  is not needed.             */
/*     (time for dpoco) = (1 + 18/n)*(time for dpofa) .                     */
/*                                                                          */
/*     on entry                                                             */
/*                                                                          */
/*        a       double precision(n,ma)                                    */
/*                the symmetric matrix to be factored.  only the            */
/*                diagonal and upper triangle are used.                     */
/*                                                                          */
/*        ma      the column dimension of the array  a .                    */
/*                                                                          */
/*        n       number of rows in a.                                      */
/*                                                                          */
/*     on return                                                            */
/*                                                                          */
/*        a       an upper triangular matrix  r  so that  a = trans(r)*r    */
/*                where  trans(r)  is the transpose.                        */
/*                the strict lower triangle is unaltered.                   */
/*                if  info .ne. 0 , the factorization is not complete.      */
/*                                                                          */
/*        info    integer                                                   */
/*                = 0  for normal return.                                   */
/*                = k  signals an error condition.  the leading minor       */
/*                     of order  k  is not positive definite.               */
/*                                                                          */
/*     linpack.  this version dated 08/14/78 .                              */
/*     cleve moler, university of new mexico, argonne national lab.         */
 
int dpofa(double *a,int ma,int n)
{
    register int j,k,jj,kk;
    double s,t;

    for (j = 1; j <= n; ++j) {
        s = 0.0;   
        if (j > 1) {
            for (k = 1; k < j; ++k) {
                kk = (k - 1) * ma + k;
                t = a[(k - 1) * ma + j] - ddot(k - 1,a + k - 1,ma,a + j - 1,ma);
                if (a[kk] == 0.0)
                    return(-1);
                t /= a[kk];
                a[(k - 1) * ma + j] = t;
                s += t * t;
            }
        }
        jj = (j - 1) * ma + j;
        s = a[jj] - s;
        if (s <= 0.0)             
            return(j);   
        a[jj] = sqrt(s);
    }
    return(0);
}           

/*  ----------------------------------------------------------------------- */    
/*  drot()                                                                  */
/*                                                                          */
/*  drot(n,dx,incx,dy,incy,dc,ds)                                           */
/*                                                                          */
/*  PURPOSE  APPLY D.P. GIVENS ROTATION                                     */
/*  BLAS routine.                                                           */
/*                                                                          */

void drot(int n,double *dx,int incx,double *dy,int incy,double dc,double ds)
{
    int i,kx,ky,nsteps;
    double w,z;
 
    if(n <= 0 || (ds == 0.0 && dc == 1.0))
        return;  

    if (incx == incy && incx > 0) {             
 
        nsteps = incx * n;
        for (i = 1; i <= nsteps; i += incx) {
            w = dx[i];
            z = dy[i];
            dx[i] = dc * w + ds * z;
            dy[i] = -ds * w + dc * z;
        }
    }
    else {   
        kx = ky = 1;
        if (incx < 0)
            kx = 1 - (n-1) * incx;
        if (incy < 0)
            ky = 1 - (n-1) * incy;
 
        for (i = 1; i <= n; ++i) {
            w = dx[kx];
            z = dy[ky];
            dx[kx] = dc * w + ds * z;
            dy[ky] = -ds * w + dc * z;
            kx += incx;
            ky += incy;
        }
    }
}

/*  ----------------------------------------------------------------------- */    
/*  drotg()                                                                 */
/*                                                                          */
/*  drotg(da,db,dc,ds)                                                      */
/*                                                                          */
/*  BLAS routine.                                                           */
/*                                                                          */

void drotg(double *da,double *db,double *dc,double *ds)
{
    double r,u,v;
  
    if (fabs(*da) <= fabs(*db))
        goto L10;                         
 
    u = *da + *da;
    v = *db / u;
  
    r = sqrt(0.25 + v * v) * u;
 
    *dc = *da / r;
    *ds = v * (*dc + *dc);
    *db = *ds;
    *da = r;
    return;
 
L10:
    if (*db == 0.0)
        goto L20;           

    u = *db + *db;
    v = *da / u;
 
    *da = sqrt(0.25 + v * v) * u;
 
    *ds = *db / *da;
    *dc = v * (*ds + *ds);
    if (*dc == 0.0)
        goto L15; 

    *db = 1.0 / *dc;
    return;

L15:
    *db = 1.0;
    return;  
 
L20:
    *dc = 1.0;
    *ds = 0.0;
    return;  
} 

/*  ----------------------------------------------------------------------- */    
/*  dxpy()                                                                  */
/*                                                                          */
/*  dxpy(n,m,x,y,xpy)                                         */
/*                                                                          */
/*  PURPOSE  COMPUTE XPY = X + Y                                            */
/*                                                                          */

void dxpy(int n,int m,double *x,double *y,double *xpy)
{
    register int i,nm;   
 
    nm = n * m;
    for (i = 1; i <= nm; ++i)  
        xpy[i] = x[i] + y[i];  
}

/*  ----------------------------------------------------------------------- */    
/*  dzero()                                                                 */
/*                                                                          */
/*  dzero(n,m,a,cda)                                                        */
/*                                                                          */
/*  PURPOSE  set a[(i - 1) * cda + j] = 0.                                  */
/*                                                                          */

void dzero(int n,int m,double *a,int cda)
{
    register int i,j;

    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= m; ++j)
            a[(i - 1) * cda + j] = 0.0;
    }
}

/*  ----------------------------------------------------------------------- */    
/*  dnrm2()                                                                 */
/*                                                                          */
/*  double precision function dnrm2(n,dx,incx)                              */
/*  Euclidean norm of dx[].                                                 */
/*  BLAS.                                                                   */
/*                                                                          */

double dnrm2(int n,double *dx,int incx)
{
    register int i,j;
    int next,nn;
    double hitest,sum,xmax,nrm;  

    double cutlo = 8.232-11;
    double cuthi = 1.30419;
 
    if(n <= 0)           
        return(0.0);

    next = 30;
    xmax = sum = 0.0;
    nn = n * incx;
    i = 1;

L20:
    if (next == 50)
        goto L50;
    else if (next == 70)
        goto L70;
    else if (next == 110)
        goto L110;

    if(fabs(dx[i]) > cutlo)
        goto L85;

    next = 50;
    xmax = 0.0;  
 
L50:
    if (dx[i] == 0.0)
        goto L200;

    if (fabs(dx[i]) > cutlo)
        goto L85;
 
    next = 70;
    goto L105;
 
L100:
    i = j;
    next = 110;
    sum = (sum / dx[i]) / dx[i];

L105:
    xmax = fabs(dx[i]);
    goto L115;
 
L70:
    if (fabs(dx[i]) > cutlo)
        goto L75;
 
L110:
    if (fabs(dx[i]) <= xmax)
        goto L115;

    sum = 1.0 + sum * (xmax / dx[i]) * (xmax / dx[i]);
    xmax = fabs(dx[i]);
    goto L200;
 
L115:
    sum += (dx[i] / xmax) * (dx[i] / xmax);
    goto L200;
 
L75:
    sum = (sum * xmax) * xmax;
 
L85:
    hitest = cuthi / (double)n;      
 
    for (j = i; j <= nn; j += incx) {
        if (fabs(dx[j]) >= hitest)
            goto L100;
        sum += dx[j] * dx[j];
    }
    nrm = sqrt(sum);
    return(nrm);
 
L200:
    i += incx;
    if (i <= nn)
        goto L20;   
 
    nrm = xmax * sqrt(sum);
    return(nrm);
}             

/*  ----------------------------------------------------------------------- */  
/*  dqrdc()                                                                 */
/*                                                                          */
/*  dqrdc(x,cdx,n,p,qraux,jpvt,work,job)                                    */
/*                                                                          */
/*  PURPOSE  USES HOUSEHOLDER TRANSFORMATIONS TO COMPUTE THE QR FACTORI-    */
/*           ZATION OF N BY P MATRIX X.  COLUMN PIVOTING IS OPTIONAL.       */
/*                                                                          */
/*                                                                          */
/*  ON ENTRY                                                                */
/*        X       DOUBLE PRECISION(LDX,P), WHERE LDX .GE. N.                */
/*                X CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO BE        */   
/*                COMPUTED.                                                 */
/*        CDX     INTEGER.                                                  */
/*                CDX IS THE COLUMN DIMENSION OF THE ARRAY X.               */
/*                This has been changed to accomodate for C conventions!    */
/*                                                                          */
/*        N       INTEGER.                                                  */
/*                N IS THE NUMBER OF ROWS OF THE MATRIX X.                  */
/*        P       INTEGER.                                                  */
/*                P IS THE NUMBER OF COLUMNS OF THE MATRIX X.               */
/*        JPVT    INTEGER(P).                                               */
/*                JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION         */
/*                OF THE PIVOT COLUMNS.  THE K-TH COLUMN X(K) OF X          */
/*                IS PLACED IN ONE OF THREE CLASSES ACCORDING TO THE        */
/*                VALUE OF JPVT(K).                                         */
/*                   IF JPVT(K) .GT. 0, THEN X(K) IS AN INITIAL             */
/*                                      COLUMN.                             */
/*                   IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE COLUMN.         */
/*                   IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL COLUMN.        */
/*                BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL COLUMNS     */
/*                ARE MOVED TO THE BEGINNING OF THE ARRAY X AND FINAL       */
/*                COLUMNS TO THE END.  BOTH INITIAL AND FINAL COLUMNS       */
/*                ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY       */
/*                FREE COLUMNS ARE MOVED.  AT THE K-TH STAGE OF THE         */
/*                REDUCTION, IF X(K) IS OCCUPIED BY A FREE COLUMN           */
/*                IT IS INTERCHANGED WITH THE FREE COLUMN OF LARGEST        */
/*                REDUCED NORM.  JPVT IS NOT REFERENCED IF                  */
/*                JOB .EQ. 0.                                               */
/*        WORK    DOUBLE PRECISION(P).                                      */
/*                WORK IS A WORK ARRAY.  WORK IS NOT REFERENCED IF          */
/*                JOB .EQ. 0.                                               */
/*        JOB     INTEGER.                                                  */
/*                JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING.         */
/*                IF JOB .EQ. 0, NO PIVOTING IS DONE.                       */
/*                IF JOB .NE. 0, PIVOTING IS DONE.                          */
/*     ON RETURN                                                            */
/*        X       X CONTAINS IN ITS UPPER TRIANGLE THE UPPER                */
/*                TRIANGULAR MATRIX R OF THE QR FACTORIZATION.              */
/*                BELOW ITS DIAGONAL X CONTAINS INFORMATION FROM            */
/*                WHICH THE ORTHOGONAL PART OF THE DECOMPOSITION            */
/*                CAN BE RECOVERED.  NOTE THAT IF PIVOTING HAS              */
/*                BEEN REQUESTED, THE DECOMPOSITION IS NOT THAT             */
/*                OF THE ORIGINAL MATRIX X BUT THAT OF X                    */
/*                WITH ITS COLUMNS PERMUTED AS DESCRIBED BY JPVT.           */
/*        QRAUX   DOUBLE PRECISION(P).                                      */
/*                QRAUX CONTAINS FURTHER INFORMATION REQUIRED TO RECOVER    */
/*                THE ORTHOGONAL PART OF THE DECOMPOSITION.                 */
/*        JPVT    JPVT(K) CONTAINS THE INDEX OF THE COLUMN OF THE           */
/*                ORIGINAL MATRIX THAT HAS BEEN INTERCHANGED INTO           */
/*                THE K-TH COLUMN, IF PIVOTING WAS REQUESTED.               */
/*     LINPACK.  THIS VERSION DATED 08/14/78 .                              */
/*     G. W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.         */
/*                                                                          */ 
/*  REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,        */
/*                 *LINPACK USERS  GUIDE*, SIAM, 1979.                      */

void dqrdc(double *x,int cdx,int n,int p,double *qraux,int *jpvt,      
    double *work,int job)
{
    register int j,jj; 
    int jp,l,lp1,lup,maxj,pl,pu,negj,swapj;
    double maxnrm,nrmxl,t,tt;
 
    pl = 1;
    pu = 0;
    if (job == 0)            
        goto L60;
 
    /*  PIVOTING HAS BEEN REQUESTED.  REARRANGE THE COLUMNS
        ACCORDING TO JPVT. */

    for (j = 1; j <= p; ++j) {
        swapj = 0;
        if (jpvt[j] > 0)
            swapj = 1;

        negj = 0;
        if (jpvt[j] < 0)
            negj = 1;

        jpvt[j] = j;
        if (negj)
            jpvt[j] = -j;

        if (swapj) {           
            if (j != pl)
                dswap(n,x+pl-1,cdx,x+j-1,cdx);
            jpvt[j] = jpvt[pl];
            jpvt[pl] = j;
            pl++;       
        }             
    }
    pu = p;
    for (jj = 1; jj <= p; ++jj) {
        j = p - jj + 1;
        if (jpvt[j] < 0) {         
            jpvt[j] = -jpvt[j];
            if (j != pu) {              
                dswap(n,x+pu-1,cdx,x+j-1,cdx);
                jp = jpvt[pu];
                jpvt[pu] = jpvt[j];
                jpvt[j] = jp;
            }
            pu--;        
        }
    }

L60:   /* COMPUTE THE NORMS OF THE FREE COLUMNS. */
 
    for (j = pl; j <= pu; ++j) {
        qraux[j] = dnrm2(n,x + j - 1,cdx);                  
        work[j] = qraux[j];
    }
 
    /* PERFORM THE HOUSEHOLDER REDUCTION OF X. */
 
    lup = imin(n,p);

    for (l = 1; l <= lup; ++l) {
        if (l < pl || l >= pu)
            goto L120;

        /* LOCATE THE COLUMN OF LARGEST NORM AND BRING IT
           INTO THE PIVOT POSITION. */
 
        maxnrm = 0.0;  
        maxj = l;
        for (j = l; j <= pu; ++j) {
            if (qraux[j] > maxnrm) {           
                maxnrm = qraux[j];
                maxj = j;
            }
        }
        if (maxj != l) {           
            dswap(n,x+l-1,cdx,x+maxj-1,cdx);
            qraux[maxj] = qraux[l];
            work[maxj] = work[l];
            jp = jpvt[maxj];
            jpvt[maxj] = jpvt[l];
            jpvt[l] = jp;
        }
L120:
        qraux[l] = 0.0;  
        if (l == n)
            continue; 
 
        /* COMPUTE THE HOUSEHOLDER TRANSFORMATION FOR COLUMN L. */

        nrmxl = dnrm2(n-l+1,x+(l-1)*cdx+l-1,cdx);
        if (nrmxl == 0.0)
            continue; 

        if (x[(l-1)*cdx+l] != 0.0)
            nrmxl = fsign(nrmxl,x[(l-1)*cdx + l]);

        dscal(n-l+1,1.0/nrmxl,x+(l-1)*cdx+l-1,cdx);
        x[(l-1)*cdx+l] += 1.0;

        /* APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS,
           UPDATING THE NORMS. */
 
        lp1 = l + 1;
        if (p >= lp1) {           
            for (j = lp1; j <= p; ++j) {
                t = -ddot(n-l+1,x+(l-1)*cdx+l-1,cdx,x+(l-1)*cdx+j-1,cdx)/
                                                                x[(l-1)*cdx+l];                

                daxpy(n-l+1,t,x+(l-1)*cdx+l-1,cdx,x+(l-1)*cdx+j-1,cdx);
                if (j < pl || j > pu)
                    goto L150;
                if (qraux[j] != 0.0) {          
                    tt = 1.0 - pow(fabs(x[(l-1)*cdx+j]) / qraux[j],2.0);
                    tt = dmax(tt,0.0);
                    t = tt;
                    tt = 1.0 + 0.05 * tt * pow(qraux[j] / work[j],2.0);
                    if (tt != 1.0)
                        qraux[j] *= sqrt(t);
                    else {
                        qraux[j] = dnrm2(n-l,x+l*cdx+j-1,cdx);
                        work[j] = qraux[j];
                    }
                }
L150:           ;          
            }
        }

        /* SAVE THE TRANSFORMATION. */
 
        qraux[l] = x[(l-1)*cdx+l];
        x[(l-1)*cdx+l] = -nrmxl;
    }
}

/*  ----------------------------------------------------------------------- */  
/*  dqrsl()                                                                 */
/*                                                                          */
/*  dqrsl(x,cdx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)                       */
/*                                                                          */
/*  PURPOSE  APPLIES THE OUTPUT OF DQRDC TO COMPUTE COORDINATE              */
/*            TRANSFORMATIONS, PROJECTIONS, AND LEAST SQUARES SOLUTIONS.    */
/*                                                                          */
/*     DQRSL APPLIES THE OUTPUT OF DQRDC TO COMPUTE COORDINATE              */
/*     TRANSFORMATIONS, PROJECTIONS, AND LEAST SQUARES SOLUTIONS.           */
/*     FOR K .LE. MIN(N,P), LET XK BE THE MATRIX                            */
/*            XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K)))                 */
/*     FORMED FROM COLUMNNS JPVT(1), ... ,JPVT(K) OF THE ORIGINAL           */
/*     N X P MATRIX X THAT WAS INPUT TO DQRDC (IF NO PIVOTING WAS           */
/*     DONE, XK CONSISTS OF THE FIRST K COLUMNS OF X IN THEIR               */
/*     ORIGINAL ORDER).  DQRDC PRODUCES A FACTORED ORTHOGONAL MATRIX Q      */
/*     AND AN UPPER TRIANGULAR MATRIX R SUCH THAT                           */
/*              XK = Q * (R)                                                */
/*                       (0)                                                */
/*     THIS INFORMATION IS CONTAINED IN CODED FORM IN THE ARRAYS            */
/*     X AND QRAUX.                                                         */
/*     ON ENTRY                                                             */
/*        X      DOUBLE PRECISION(LDX,P).                                   */
/*               X CONTAINS THE OUTPUT OF DQRDC.                            */
/*        CDX    INTEGER.                                                   */
/*               CDX IS THE COLUMN DIMENSION OF THE ARRAY X.                */
/*        N      INTEGER.                                                   */
/*               N IS THE NUMBER OF ROWS OF THE MATRIX XK.  IT MUST         */
/*               HAVE THE SAME VALUE AS N IN DQRDC.                         */
/*        K      INTEGER.                                                   */
/*               K IS THE NUMBER OF COLUMNS OF THE MATRIX XK.  K            */
/*               MUST NOT BE GREATER THAN MIN(N,P), WHERE P IS THE          */
/*               SAME AS IN THE CALLING SEQUENCE TO DQRDC.                  */
/*        QRAUX  DOUBLE PRECISION(P).                                       */
/*               QRAUX CONTAINS THE AUXILIARY OUTPUT FROM DQRDC.            */
/*        Y      DOUBLE PRECISION(N)                                        */
/*               Y CONTAINS AN N-VECTOR THAT IS TO BE MANIPULATED           */
/*               BY DQRSL.                                                  */
/*        JOB    INTEGER.                                                   */
/*               JOB SPECIFIES WHAT IS TO BE COMPUTED.  JOB HAS             */
/*               THE DECIMAL EXPANSION ABCDE, WITH THE FOLLOWING            */
/*               MEANING.                                                   */
/*                    IF A .NE. 0, COMPUTE QY.                              */
/*                    IF B,C,D, OR E .NE. 0, COMPUTE QTY.                   */
/*                    IF C .NE. 0, COMPUTE B.                               */
/*                    IF D .NE. 0, COMPUTE RSD.                             */
/*                    IF E .NE. 0, COMPUTE XB.                              */
/*               NOTE THAT A REQUEST TO COMPUTE B, RSD, OR XB               */
/*               AUTOMATICALLY TRIGGERS THE COMPUTATION OF QTY, FOR         */
/*               WHICH AN ARRAY MUST BE PROVIDED IN THE CALLING             */
/*               SEQUENCE.                                                  */
/*     ON RETURN                                                            */
/*        QY     DOUBLE PRECISION(N).                                       */
/*               QY CONTAINS Q*Y, IF ITS COMPUTATION HAS BEEN               */
/*               REQUESTED.                                                 */
/*        QTY    DOUBLE PRECISION(N).                                       */
/*               QTY CONTAINS TRANS(Q)*Y, IF ITS COMPUTATION HAS            */
/*               BEEN REQUESTED.  HERE TRANS(Q) IS THE                      */
/*               TRANSPOSE OF THE MATRIX Q.                                 */
/*        B      DOUBLE PRECISION(K)                                        */
/*               B CONTAINS THE SOLUTION OF THE LEAST SQUARES PROBLEM       */
/*                    MINIMIZE NORM2(Y - XK*B),                             */
/*               IF ITS COMPUTATION HAS BEEN REQUESTED.  (NOTE THAT         */
/*               IF PIVOTING WAS REQUESTED IN DQRDC, THE J-TH               */
/*               COMPONENT OF B WILL BE ASSOCIATED WITH COLUMN JPVT(J)      */
/*               OF THE ORIGINAL MATRIX X THAT WAS INPUT INTO DQRDC.)       */
/*        RSD    DOUBLE PRECISION(N).                                       */
/*               RSD CONTAINS THE LEAST SQUARES RESIDUAL Y - XK*B,          */
/*               IF ITS COMPUTATION HAS BEEN REQUESTED.  RSD IS             */
/*               ALSO THE ORTHOGONAL PROJECTION OF Y ONTO THE               */
/*               ORTHOGONAL COMPLEMENT OF THE COLUMN SPACE OF XK.           */
/*        XB     DOUBLE PRECISION(N).                                       */
/*               XB CONTAINS THE LEAST SQUARES APPROXIMATION XK*B,          */
/*               IF ITS COMPUTATION HAS BEEN REQUESTED.  XB IS ALSO         */
/*               THE ORTHOGONAL PROJECTION OF Y ONTO THE COLUMN SPACE       */
/*               OF X.                                                      */
/*        INFO   INTEGER.                                                   */
/*               INFO IS ZERO UNLESS THE COMPUTATION OF B HAS               */
/*               BEEN REQUESTED AND R IS EXACTLY SINGULAR.  IN              */
/*               THIS CASE, INFO IS THE INDEX OF THE FIRST ZERO             */
/*               DIAGONAL ELEMENT OF R AND B IS LEFT UNALTERED.             */
/*     THE PARAMETERS QY, QTY, B, RSD, AND XB ARE NOT REFERENCED            */
/*     IF THEIR COMPUTATION IS NOT REQUESTED AND IN THIS CASE               */
/*     CAN BE REPLACED BY DUMMY VARIABLES IN THE CALLING PROGRAM.           */
/*     TO SAVE STORAGE, THE USER MAY IN SOME CASES USE THE SAME             */
/*     ARRAY FOR DIFFERENT PARAMETERS IN THE CALLING SEQUENCE.  A           */
/*     FREQUENTLY OCCURING EXAMPLE IS WHEN ONE WISHES TO COMPUTE            */
/*     ANY OF B, RSD, OR XB AND DOES NOT NEED Y OR QTY.  IN THIS            */
/*     CASE ONE MAY IDENTIFY Y, QTY, AND ONE OF B, RSD, OR XB, WHILE        */
/*     PROVIDING SEPARATE ARRAYS FOR ANYTHING ELSE THAT IS TO BE            */
/*     COMPUTED.  THUS THE CALLING SEQUENCE                                 */
/*          CALL DQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO)            */
/*     WILL RESULT IN THE COMPUTATION OF B AND RSD, WITH RSD                */
/*     OVERWRITING Y.  MORE GENERALLY, EACH ITEM IN THE FOLLOWING           */
/*     LIST CONTAINS GROUPS OF PERMISSIBLE IDENTIFICATIONS FOR              */
/*     A SINGLE CALLING SEQUENCE.                                           */
/*          1. (Y,QTY,B) (RSD) (XB) (QY)                                    */
/*          2. (Y,QTY,RSD) (B) (XB) (QY)                                    */
/*          3. (Y,QTY,XB) (B) (RSD) (QY)                                    */
/*          4. (Y,QY) (QTY,B) (RSD) (XB)                                    */
/*          5. (Y,QY) (QTY,RSD) (B) (XB)                                    */
/*          6. (Y,QY) (QTY,XB) (B) (RSD)                                    */
/*     IN ANY GROUP THE VALUE RETURNED IN THE ARRAY ALLOCATED TO            */
/*     THE GROUP CORRESPONDS TO THE LAST MEMBER OF THE GROUP.               */
/*     LINPACK.  THIS VERSION DATED 08/14/78 .                              */
/*     G. W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.         */
/*     REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,     */
/*                 *LINPACK USERS  GUIDE*, SIAM, 1979.                      */
/*     ROUTINES CALLED  DAXPY,DCOPY,DDOT                                    */

void dqrsl(double *x,int cdx,int n,int k,double *qraux,double *y,double *qy,
    double *qty,double *b,double *rsd,double *xb,int job,int *info)
{
    register int i,j,jj;
    int ju,kp1,cb,cqty,cqy,cr,cxb;
    double t,temp;
 
    /* DETERMINE WHAT IS TO BE COMPUTED. */
 
    *info = cqy = cqty = cb = cr = cxb = 0;

    if (job/10000 != 0)  
        cqy = 1;
    if ((job % 10000) != 0)
        cqty = 1;
    if ((job % 1000)/100 != 0)
        cb = 1;
    if ((job % 100)/10 != 0)
        cr = 1;
    if ((job % 10) != 0)
        cxb = 1;
    ju = imin(k,n-1);
 
    /* SPECIAL ACTION WHEN N=1. */
 
    if (ju == 0) {
        if (cqy)
            qy[1] = y[1];
        if (cqty)
            qty[1] = y[1];
        if (cxb)
            xb[1] = y[1];
        if (cb) {
            if (x[1] != 0.0)
                b[1] = y[1] / x[1];
            else
                *info = 1;
        }
        if (cr)
            rsd[1] = 0.0;  
        return;
    }
 
    /* SET UP TO COMPUTE QY OR QTY. */
 
    if (cqy)
        dcopy(n,y,1,qy,1);
    if (cqty)
        dcopy(n,y,1,qty,1);
    if (!cqy)
        goto L70;
 
    /* COMPUTE QY. */
  
    for (jj = 1; jj <= ju; ++jj) {
        j = ju - jj + 1;
        if (qraux[j] == 0.0)              
            continue;
        temp = x[(j-1)*cdx+j];
        x[(j-1)*cdx+j] = qraux[j];
        t = -ddot(n-j+1,x+(j-1)*cdx+j-1,cdx,qy+j-1,1) / x[(j-1)*cdx+j];
        daxpy(n-j+1,t,x+(j-1)*cdx+j-1,cdx,qy+j-1,1);
        x[(j-1)*cdx+j] = temp;
    }
     
L70:
    if (!cqty)
        goto L100;
 
    /* COMPUTE TRANS(Q)*Y. */
 
    for (j = 1; j <= ju; ++j) {
        if (qraux[j] == 0.0)            
            continue;
        temp = x[(j-1)*cdx+j];
        x[(j-1)*cdx+j] = qraux[j];
        t = -ddot(n-j+1,x+(j-1)*cdx+j-1,cdx,qty+j-1,1) / x[(j-1)*cdx+j];
        daxpy(n-j+1,t,x+(j-1)*cdx+j-1,cdx,qty+j-1,1);
        x[(j-1)*cdx+j] = temp;
    }

L100:  /* SET UP TO COMPUTE B, RSD, OR XB. */
 
    if (cb)
        dcopy(k,qty,1,b,1);
    kp1 = k + 1;
    if (cxb)
        dcopy(k,qty,1,xb,1);
    if (cr && k < n)
        dcopy(n-k,qty+kp1-1,1,rsd+kp1-1,1);

    if (cxb && kp1 <= n) {
        for (i = kp1; i <= n; ++i) {
            xb[i] = 0.0;  
        }
    }
    if (cr) {            
        for (i = 1; i <= k; ++i)  
            rsd[i] = 0.0;  
    }
    if (cb) {       /* compute B */
 
        for (jj = 1; jj <= k; ++jj) {
            j = k - jj + 1;
            if (x[(j-1)*cdx+j] == 0.0) {            
                *info = j;
                break;
            }
            b[j] = b[j]/x[(j-1)*cdx+j];
            if (j != 1) {           
                t = -b[j];
                daxpy(j-1,t,x+j-1,cdx,b,1);
            }
        }
    }
    if (cr || cxb) {    /* COMPUTE RSD OR XB AS REQUIRED. */
 
        for (jj = 1; jj <= ju; ++jj) {
            j = ju - jj + 1;
            if (qraux[j] != 0.0) {           
                temp = x[(j-1)*cdx+j];
                x[(j-1)*cdx+j] = qraux[j];
                if (cr) {              
                    t = -ddot(n-j+1,x+(j-1)*cdx+j-1,cdx,rsd+j-1,1) / x[(j-1)*cdx+j];
                    daxpy(n-j+1,t,x+(j-1)*cdx+j-1,1,rsd+j-1,1);
                }
                if (cxb) {             
                    t = -ddot(n-j+1,x+(j-1)*cdx+j-1,cdx,xb+j-1,1) / x[(j-1)*cdx+j];
                    daxpy(n-j+1,t,x+(j-1)*cdx+j-1,cdx,xb+j-1,1);
                }
                x[(j-1)*cdx+j] = temp;
            }
        }
    }
}

/*  ----------------------------------------------------------------------- */  
/*  dtrco()                                                                 */
/*                                                                          */
/*  dtrco(t,cdt,n,rcond,z,job)                                              */
/*                                                                          */
/*   PURPOSE  ESTIMATES THE CONDITION OF A DOUBLE PRECISION TRIANGULAR      */
/*            MATRIX.                                                       */
/*   DESCRIPTION                                                            */
/*     DTRCO ESTIMATES THE CONDITION OF A DOUBLE PRECISION TRIANGULAR       */
/*     MATRIX.                                                              */
/*     ON ENTRY                                                             */
/*        T       DOUBLE PRECISION(LDT,N)                                   */
/*                T CONTAINS THE TRIANGULAR MATRIX.  THE ZERO               */
/*                ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND            */
/*                THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE            */
/*                USED TO STORE OTHER INFORMATION.                          */
/*        CDT     INTEGER                                                   */
/*                CDT IS THE COLUMN DIMENSION OF THE ARRAY T.               */
/*        N       INTEGER                                                   */
/*                N IS THE ORDER OF THE SYSTEM.                             */
/*        JOB     INTEGER                                                   */
/*                = 0         T  IS LOWER TRIANGULAR.                       */
/*                = NONZERO   T  IS UPPER TRIANGULAR.                       */
/*     ON RETURN                                                            */
/*        RCOND   DOUBLE PRECISION                                          */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  T .           */
/*                FOR THE SYSTEM  T*X = B , RELATIVE PERTURBATIONS          */
/*                IN  T  AND  B  OF SIZE  EPSILON  MAY CAUSE                */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .    */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION        */
/*                           1.0 + RCOND .EQ. 1.0                           */
/*                IS TRUE, THEN  T  MAY BE SINGULAR TO WORKING              */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF            */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE             */
/*                UNDERFLOWS.                                               */
/*        Z       DOUBLE PRECISION(N)                                       */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.     */
/*                IF  T  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS         */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT              */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .                       */
/*     LINPACK.  THIS VERSION DATED 08/14/78 .                              */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.         */
/*     REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,     */
/*                 *LINPACK USERS  GUIDE*, SIAM, 1979.                      */
/*     ROUTINES CALLED  DASUM,DAXPY,DSCAL                                   */

void dtrco(double *t,int cdt,int n,double *rcond,double *z,int job)
{
    int i1,j,j1,j2,k,kk,l,lower;
    double ek,s,sm,tnorm,w,wk,wkm,ynorm;
  
    lower = 0;
    if (job == 0)          
        lower = 1;
 
    /* COMPUTE 1-NORM OF T */
 
    tnorm = 0.0;  
    for (j = 1; j <= n; ++j) {
        l = j;
        if (lower)
            l = n + 1 - j;
        i1 = 1;
        if (lower)
            i1 = j;
        tnorm = dmax(tnorm,dasum(l,t+(i1-1)*cdt+j-1,cdt));
    }
 
    /* RCOND = 1/(NORM(T)*(ESTIMATE OF NORM(INVERSE(T)))) .
       ESTIMATE = NORM(Z)/NORM(Y) WHERE  T*Z = Y  AND  TRANS(T)*Y = E .
       TRANS(T)  IS THE TRANSPOSE OF T .
       THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
       GROWTH IN THE ELEMENTS OF Y .
       THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */
 
    /* SOLVE TRANS(T)*Y = E */
 
    ek = 1.0;  
    for (j = 1; j <= n; ++j) 
        z[j] = 0.0;  

    for (kk = 1; kk <= n; ++kk) {
        k = kk;
        if (lower)
            k = n + 1 - kk;
        if (z[k] != 0.0)
            ek = fsign(ek,-z[k]);
  
        if (fabs(ek-z[k]) > fabs(t[(k-1)*cdt+k])) {          
            s = fabs(t[(k-1)*cdt+k]) / fabs(ek-z[k]);
            dscal(n,s,z,1);
            ek *= s;   
        }
   
        wk = ek - z[k];
        wkm = -ek - z[k];
        s = fabs(wk);
        sm = fabs(wkm);
  
        if (t[(k-1)*cdt+k] != 0.0) {           
            wk = wk / t[(k-1)*cdt+k];
            wkm = wkm / t[(k-1)*cdt+k];
        }
        else {
            wk = 1.0;  
            wkm = 1.0;  
        }
        if (kk != n) {        
            j1 = k + 1;
            if (lower)
                j1 = 1;
            j2 = n;
            if (lower)
                j2 = k - 1;
   
            for (j = j1; j <= j2; ++j) {
                sm += fabs(z[j]+wkm * t[(k-1)*cdt+j]);            
                z[j] += wk * t[(k-1)*cdt+j];        
                s += fabs(z[j]);
            }
  
            if (s < sm) {        
                w = wkm - wk;
                wk = wkm;
                for (j = j1; j <= j2; ++j) {
                       z[j] += w * t[(k-1)*cdt+j];                   
                }
            }
        }
        z[k] = wk;
    }
    s = 1.0 / dasum(n,z,1);
    dscal(n,s,z,1);
    ynorm = 1.0;  
 
    /* SOLVE T*Z = Y */
 
    for (kk = 1; kk <= n; ++kk) {
        k = n + 1 - kk;
        if (lower)
            k = kk;
   
        if (fabs(z[k]) > fabs(t[(k-1)*cdt+k])) {
            s = fabs(t[(k-1)*cdt+k]) / fabs(z[k]);
            dscal(n,s,z,1);
            ynorm *= s;      
        }
        if (t[(k-1)*cdt+k] != 0.0)
            z[k] /=  t[(k-1)*cdt+k];
        if (t[(k-1)*cdt+k] == 0.0)
            z[k] = 1.0;  
        i1 = 1;
        if (lower)
            i1 = k + 1;
        if (kk < n) {          
            w = -z[k];
            daxpy(n-kk,w,t+(i1-1)*cdt+k-1,cdt,z+i1-1,1);
        }
    }

    /* MAKE ZNORM = 1.0 */

    s = 1.0 / dasum(n,z,1);
    dscal(n,s,z,1);
    ynorm *= s;          
    if (tnorm != 0.0)
        *rcond = ynorm/tnorm;
    else                  
        *rcond = 0.0;  
}

/*  ----------------------------------------------------------------------- */  
/*  idamax()                                                                */
/*                                                                          */
/*  int idamax(n,dx,incx)                                                   */
/*                                                                          */
/*  PURPOSE  FIND LARGEST COMPONENT OF D.P. VECTOR                          */
/*  BLAS.                                                                   */

int idamax(int n,double *dx,int incx)
{
    int i,ii,ns,imax;
    double dmax,xmag;
 
    if (n <= 0)
        return(0);
    if (n <= 1)
        return(1);

    imax = 1;
    if (incx != 1) {       
        dmax = fabs(dx[1]);
        ns = n * incx;
        ii = 1;
        for (i = 1; i <= ns; i += incx) {
            xmag = fabs(dx[i]);
            if (xmag > dmax) {        
                imax = ii;
                dmax = xmag;
            }
            ii++;        
        }
    }
    else {
        dmax = fabs(dx[1]);
        for (i = 2; i <= n; ++i) {
            xmag = fabs(dx[i]);
            if (xmag <= dmax)
                continue;             
            imax = i;
            dmax = xmag; 
        }
    }
    return(imax);
}

/*  ----------------------------------------------------------------------- */  
/*  dswap()                                                                 */
/*                                                                          */
/*  dswap(n,dx,incx,dy,incy)                                                */
/*                                                                          */  
/*  PURPOSE  INTERCHANGE D.P. VECTORS                                       */
/*  BLAS.                                                                   */

void dswap(int n,double *dx,int incx,double *dy,int incy)
{
    int i,ix,iy,m,ns;
    double dtemp1,dtemp2,dtemp3;
 
    if (n <= 0)
        return;     
    if (incx == incy) {
        if (incx == 1)
            goto L20;
        else if (incx > 1)
            goto L60;
    }
    ix = iy = 1;
    if (incx < 0)
        ix = (-n+1)*incx + 1;
    if (incy < 0)
        iy = (-n+1)*incy + 1;

    for (i = 1; i <= n; ++i) {
        dtemp1 = dx[ix];
        dx[ix] = dy[iy];
        dy[iy] = dtemp1;
        ix = ix + incx;
        iy = iy + incy;
    }
    return;

L20:
    m = n % 3;
    if (m == 0)
        goto L40;             
    for (i = 1; i <= m; ++i) {
        dtemp1 = dx[i];
        dx[i] = dy[i];
        dy[i] = dtemp1;
    }
    if (n < 3)
        return;                

L40:
    for (i = m + 1; i <= n; i += 3) {
        dtemp1 = dx[i];
        dtemp2 = dx[i+1];
        dtemp3 = dx[i+2];
        dx[i] = dy[i];
        dx[i+1] = dy[i+1];
        dx[i+2] = dy[i+2];
        dy[i] = dtemp1;
        dy[i+1] = dtemp2;
        dy[i+2] = dtemp3;
    }
    return;  

L60:
    ns = n * incx;
    for (i = 1; i <= ns; i += incx) {
        dtemp1 = dx[i];
        dx[i] = dy[i];
        dy[i] = dtemp1;
    }
    return;
}            

/*  ----------------------------------------------------------------------- */  
/*  dasum()                                                                 */
/*                                                                          */
/*  double dasum(n,dx,incx)                                                 */
/*                                                                          */
/*  PURPOSE  SUM OF MAGNITUDES OF D.P. VECTOR COMPONENTS                    */
/*  BLAS.                                                                   */

double dasum(int n,double *dx,int incx)
{
    register int i,m,ns;
    double ds;

    ds = 0.0;
    if (n <= 0)
        return(ds);    

    if (incx != 1) {       
        ns = n * incx;
        for (i = 1; i <= ns; i += incx)  
            ds += fabs(dx[i]);
        return(ds);
    }
    m = n % 6;
    if (m != 0) {
        for (i = 1; i <= m; ++i) 
            ds += fabs(dx[i]);
        if (n < 6) 
            return(ds);
    }
    for (i = m + 1; i <= n; i += 6) {
        ds += fabs(dx[i]) + fabs(dx[i+1]) + fabs(dx[i+2]) +
              fabs(dx[i+3]) + fabs(dx[i+4]) + fabs(dx[i+5]);
    }
    return(ds);
}

/*  ----------------------------------------------------------------------- */  
/*  dtrsl()                                                                 */
/*                                                                          */
/*  dtrsl(t,cdt,n,b,job,info)                                               */
/*                                                                          */
/*  PURPOSE  SOLVES SYSTEMS OF THE FORM  T*X=B OR  TRANS(T)*X=B  WHERE T    */
/*            IS A TRIANGULAR MATRIX OF ORDER N.                            */
/*                                                                          */
/*     DTRSL SOLVES SYSTEMS OF THE FORM                                     */
/*                   T * X = B                                              */
/*     OR                                                                   */
/*                   TRANS(T) * X = B                                       */
/*     WHERE T IS A TRIANGULAR MATRIX OF ORDER N.  HERE TRANS(T)            */
/*     DENOTES THE TRANSPOSE OF THE MATRIX T.                               */
/*     ON ENTRY                                                             */
/*         T         DOUBLE PRECISION(LDT,N)                                */
/*                   T CONTAINS THE MATRIX OF THE SYSTEM.  THE ZERO         */
/*                   ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND         */
/*                   THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE         */
/*                   USED TO STORE OTHER INFORMATION.                       */
/*         CDT       INTEGER                                                */
/*                   CDT IS THE COLUMN  DIMENSION OF THE ARRAY T.           */
/*         N         INTEGER                                                */
/*                   N IS THE ORDER OF THE SYSTEM.                          */
/*         B         DOUBLE PRECISION(N).                                   */
/*                   B CONTAINS THE RIGHT HAND SIDE OF THE SYSTEM.          */
/*         JOB       INTEGER                                                */
/*                   JOB SPECIFIES WHAT KIND OF SYSTEM IS TO BE SOLVED.     */
/*                   IF JOB IS                                              */
/*                        00   SOLVE T*X=B, T LOWER TRIANGULAR,             */
/*                        01   SOLVE T*X=B, T UPPER TRIANGULAR,             */
/*                        10   SOLVE TRANS(T)*X=B, T LOWER TRIANGULAR,      */
/*                        11   SOLVE TRANS(T)*X=B, T UPPER TRIANGULAR.      */
/*     ON RETURN                                                            */
/*         B         B CONTAINS THE SOLUTION, IF INFO .EQ. 0.               */
/*                   OTHERWISE B IS UNALTERED.                              */
/*         INFO      INTEGER                                                */
/*                   INFO CONTAINS ZERO IF THE SYSTEM IS NONSINGULAR.       */
/*                   OTHERWISE INFO CONTAINS THE INDEX OF                   */
/*                   THE FIRST ZERO DIAGONAL ELEMENT OF T.                  */
/*     LINPACK.  THIS VERSION DATED 08/14/78 .                              */
/*     G. W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.         */
/*     REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,     */
/*                 *LINPACK USERS  GUIDE*, SIAM, 1979.                      */
/*                                                                          */

void dtrsl(double *t,int cdt,int n,double *b,int job,int *info)                                                 
{
    int cas,j,jj;
    double temp;       

    for (*info = 1; *info <= n; *info += 1) {
        if (t[(*info - 1) * cdt + *info] == 0.0)
            return;                              
    }
    *info = 0;
 
    cas = 1;
    if ((job % 10) != 0)
        cas = 2;
    if ((job % 100)/10 != 0)
        cas += 2;

    if (cas == 2)
        goto L50;
    else if (cas == 3)
        goto L80;
    else if (cas == 4)
        goto L110;
    /*  GO TO (20,50,80,110), cas   */
    
    b[1] = b[1] / t[1];   
    for (j = 2; j <= n; ++j) {
        temp = -b[j-1];
        daxpy(n-j+1,temp,t+(j-1)*cdt+j-2,cdt,b+j-1,1);
        b[j] = b[j]/t[(j-1)*cdt+j];        
    }
    return;
 
L50:
    b[n] = b[n] / t[(n-1)*cdt+n];          
    for (jj = 2; jj <= n; ++jj) {
        j = n - jj + 1;
        temp = -b[j+1];
        daxpy(j,temp,t+j,cdt,b,1);
        b[j] = b[j] / t[(j-1)*cdt+j];
    }
    return;
 
L80:
    b[n] = b[n]/ t[(n-1)*cdt+n];
    for (jj = 2; jj <= n; ++jj) {
        j = n - jj + 1;
        b[j] -= ddot(jj-1,t+j*cdt+j-1,cdt,b+j,1);
        b[j] = b[j]/t[(j-1)*cdt+j];
    }
    return;   
 
L110:
    b[1] = b[1]/t[1];  
    for (j = 2; j <= n; ++j) {
        b[j] -= ddot(j-1,t+j-1,cdt,b,1);
        b[j] = b[j]/t[(j-1)*cdt+j];
    }
}           




