/****************************************************************************/
/*  t_qp                                                                    */
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
#include "t_alloc.h"

/*  functions in t_qp.c */

int qld(int n,int m,int meq,int lql,int maxit,double *a,double *b,double *grad,
    double *g,double *x,double *xl,double *xu,int *iact,int *nact,
    double vsmall,double *diag);

/* ------------------------------------------------------------------------ */
/*                                                                          */
/*  int qld(int n,int m,int meq,int lql,int maxit,double *a,double *b,      */
/*      double *grad,double *g,double *x,double *xl,double *xu,             */
/*      int *iact,int *nact,double vsmall,double *diag)                     */
/*                                                                          */
/*  THIS SUBROUTINE SOLVES THE QUADRATIC PROGRAMMING PROBLEM                */
/*                                                                          */
/*      MINIMIZE      GRAD'*X  +  0.5 * X*G*X                               */
/*      SUBJECT TO    A(K)*X  =  B(K)   K=1,2,...,MEQ,                      */
/*                    A(K)*X >=  B(K)   K=MEQ+1,...,M,                      */
/*                    XL  <=  X  <=  XU                                     */
/*                                                                          */
/*  THE QUADRATIC PROGRAMMING METHOD PROCEEDS FROM AN INITIAL CHOLESKY-     */
/*  DECOMPOSITION OF THE OBJECTIVE FUNCTION MATRIX, TO CALCULATE THE        */
/*  UNIQUELY DETERMINED MINIMIZER OF THE UNCONSTRAINED PROBLEM.             */
/*  SUCCESSIVELY ALL VIOLATED CONSTRAINTS ARE ADDED TO A WORKING SET        */
/*  AND A MINIMIZER OF THE OBJECTIVE FUNCTION SUBJECT TO ALL CONSTRAINTS    */
/*  IN THIS WORKING SET IS COMPUTED. IT IS POSSIBLE THAT CONSTRAINTS        */
/*  HAVE TO LEAVE THE WORKING SET.                                          */
/*                                                                          */
/*                                                                          */
/*  DESCRIPTION OF PARAMETERS:                                              */
/*                                                                          */
/*    N        : IS THE NUMBER OF VARIABLES.                                */
/*    M        : TOTAL NUMBER OF CONSTRAINTS.                               */
/*    MEQ      : NUMBER OF EQUALITY CONTRAINTS.                             */
/*    MMAX     : ROW DIMENSION OF A, DIMENSION OF B. MMAX MUST BE AT        */
/*               LEAST ONE AND GREATER OR EQUAL TO M.                       */
/*    MN       : MUST BE EQUAL M + N.                                       */
/*    MNN      : MUST BE EQUAL M + N + N.                                   */
/*    NMAX     : ROW DIMENSION OF G. MUST BE AT LEAST N.                    */
/*    LQL      : DETERMINES INITIAL DECOMPOSITION.                          */
/*       lql = 0 : THE UPPER TRIANGULAR PART OF THE MATRIX G                */
/*                 CONTAINS INITIALLY THE CHOLESKY-FACTOR OF A SUITABLE     */
/*                 DECOMPOSITION.                                           */
/*       lql = 1 : THE INITIAL CHOLESKY-FACTORISATION OF G IS TO BE         */
/*                 PERFORMED BY THE ALGORITHM.                              */
/*    A(MMAX,NMAX) : A IS A MATRIX WHOSE COLUMNS ARE THE CONSTRAINTS NORMALS*/
/*    B(MMAX)  : CONTAINS THE RIGHT HAND SIDES OF THE CONSTRAINTS.          */
/*    grad(N)  : CONTAINS THE OBJECTIVE FUNCTION VECTOR grad.               */
/*    G(NMAX,N): CONTAINS THE SYMMETRIC OBJECTIVE FUNCTION MATRIX.          */
/*    XL(N), XU(N): CONTAIN THE LOWER AND UPPER BOUNDS FOR X.               */
/*    X(N)     : VECTOR OF VARIABLES.                                       */
/*    NACT     : FINAL NUMBER OF ACTIVE CONSTRAINTS.                        */
/*    IACT[K] (K=1,2,...,nact): INDICES OF THE FINAL ACTIVE CONSTRAINTS.    */
/*    INFO     : REASON FOR THE RETURN FROM THE SUBROUTINE.                 */
/*        info = 0 : CALCULATION WAS TERMINATED SUCCESSFULLY.               */
/*        info = 1 : MAXIMUM NUMBER OF ITERATIONS ATTAINED.                 */
/*        info = 2 : ACCURACY IS INSUFFICIENT TO MAINTAIN INCREASING        */
/*                   FUNCTION VALUES.                                       */
/*        info < 0 : THE CONSTRAINT WITH INDEX ABS(info) AND THE CON-       */
/*                   STRAINTS WHOSE INDICES ARE iact[k], K=1,2,...,nact,    */
/*                   ARE INCONSISTENT.                                      */
/*    MAXIT    : MAXIMUM NUMBER OF ITERATIONS.                              */
/*    VSMALL   : REQUIRED ACCURACY TO BE ACHIEVED (E.G. IN THE ORDER OF THE */
/*               MACHINE PRECISION FOR SMALL AND WELL-CONDITIONED PROBLEMS).*/
/*    DIAG     : ON RETURN DIAG IS EQUAL TO THE MULTIPLE OF THE UNIT MATRIX */
/*               THAT WAS ADDED TO G TO ACHIEVE POSITIVE DEFINITENESS.      */
/*    W(LW)    : THE ELEMENTS OF W(.) ARE USED FOR WORKING SPACE. THE LENGTH*/
/*               OF W MUST NOT BE LESS THAN (1.5*NMAX*NMAX + 10*NMAX + M).  */
/*               WHEN info = 0 ON RETURN, THE LAGRANGE MULTIPLIERS OF THE   */
/*               FINAL ACTIVE CONSTRAINTS ARE HELD IN W(K), K=1,2,...,nact. */
/*                                                                          */
/*  THE VALUES OF N, M, meq, MMAX, MN, MNN AND NMAX AND THE ELEMENTS OF     */
/*  A, B, grad AND G ARE NOT ALTERED.                                       */
/*                                                                          */
/*  THE FOLLOWING INTEGERS ARE USED TO PARTITION W:                         */
/*    THE FIRST N ELEMENTS OF W HOLD LAGRANGE MULTIPLIER ESTIMATES.         */
/*    W(IWZ+I+(N-1)*J) HOLDS THE MATRIX ELEMENT Z(I,J).                     */
/*    W(IWR+I+0.5*J*(J-1)) HOLDS THE UPPER TRIANGULAR MATRIX                */
/*      ELEMENT R(I,J). THE SUBSEQUENT N COMPONENTS OF W MAY BE             */
/*      TREATED AS AN EXTRA COLUMN OF R(.,.).                               */
/*    W(IWW-N+I) (I=1,2,...,N) ARE USED FOR TEMPORARY STORAGE.              */
/*    W(IWW+I) (I=1,2,...,N) ARE USED FOR TEMPORARY STORAGE.                */
/*    W(IWD+I) (I=1,2,...,N) HOLDS G(I,I) DURING THE CALCULATION.           */
/*    W(IWX+I) (I=1,2,...,N) HOLDS VARIABLES THAT WILL BE USED TO           */
/*      TEST THAT THE ITERATIONS INCREASE THE OBJECTIVE FUNCTION.           */
/*    W(IWA+K) (K=1,2,...,M) USUALLY HOLDS THE RECIPROCAL OF THE            */
/*      LENGTH OF THE K-TH CONSTRAINT, BUT ITS SIGN INDICATES               */
/*      WHETHER THE CONSTRAINT IS ACTIVE.                                   */
/*                                                                          */
/*                                                                          */
/*  AUTHOR:    K. SCHITTKOWSKI,                                             */
/*             MATHEMATISCHES INSTITUT,                                     */
/*             UNIVERSITAET BAYREUTH,                                       */
/*             8580 BAYREUTH,                                               */
/*             GERMANY, F.R.                                                */
/*                                                                          */
/*  AUTHOR OF ORIGINAL VERSION:                                             */
/*             M.J.D. POWELL, DAMTP,                                        */
/*             UNIVERSITY OF CAMBRIDGE, SILVER STREET                       */
/*             CAMBRIDGE,                                                   */
/*             ENGLAND                                                      */
/*                                                                          */
/*  REFERENCE: M.J.D. POWELL: ZQPCVX, A FORTRAN SUBROUTINE FOR CONVEX       */
/*             PROGRAMMING, REPORT DAMTP/1983/NA17, UNIVERSITY OF           */
/*             CAMBRIDGE, ENGLAND, 1983.                                    */
/*                                                                          */
/*  VERSION :  2.0 (MARCH, 1987)                                            */
/*                                                                          */
/* ------------------------------------------------------------------------ */
/*                                                                          */
/*  Adapted into C, Dec. 2003 (G.R.)                                        */
/*                                                                          */
/*  The function returns:                                                   */
/*                                                                          */ 
/*      info = 0 : CALCULATION WAS TERMINATED SUCCESSFULLY.                 */
/*      info = 1 : MAXIMUM NUMBER OF ITERATIONS ATTAINED.                   */
/*      info = 2 : ACCURACY IS INSUFFICIENT TO MAINTAIN INCREASING          */
/*                 FUNCTION VALUES.                                         */
/*      info = 3 : insufficient memory.                                     */
/*      info < 0 : THE CONSTRAINT WITH INDEX ABS(info) AND THE CON-         */
/*                 STRAINTS WHOSE INDICES ARE iact[k], K=1,2,...,nact,      */
/*                 ARE INCONSISTENT.                                        */
/*                                                                          */
/*  Note:                                                                   */
/*      1) W is locally allocated.                                          */
/*      2) nmax always equals n.                                            */
/*      3) mmax always equals m.                                            */
/*                                                                          */
/* ------------------------------------------------------------------------ */
 
int qld(int n,int m,int meq,int lql,int maxit,double *a,double *b,double *grad,
    double *g,double *x,double *xl,double *xu,int *iact,int *nact,
    double vsmall,double *diag)
{
int iii,jjj;

    int i,j,k,k1,kk,id,ir,ip,il,iu,ix,iy,iw,is,ju,nu;
    int iwz,iwr,iww,iwd,iwa,ifinc,kfinc,jfinc,iterc,itref,knext,kflag;
    int iflag,lflag,jflag,kdrop,mflag,iwx,iwy,iws,iza,ipp,ira,irb;
    int ia,ii,jl,iwwn,iz,nflag,nm,lower,lw,info,mn; 
             
    double sum,suma,sumb,sumc,sumx,sumy,cvmax,diagr,fdiffa,fdiff;
    double temp,tempa,step,ratio,parinc,parnew,xmagr,xmag,ga,gb,res;
    double *w; 
    double vfact = 1.0;

    info = 3;
    lw = (3 * n * n) / 2 + 10 * n + 2 * m + 10;
    if (!(w = (double *) calloc(lw,sizeof(double))))  
        return(3);
    memrq(lw,sizeof(double));

    mn = m + n;
    iwz = n;
    iwr = iwz + n * n;
    iww = iwr + (n * (n + 3)) / 2;
    iwd = iww + n;
    iwx = iwd + n;
    iwa = iwx + n;
  
/*  SET SOME PARAMETERS.                                        
    NUMBER LESS THAN VSMALL ARE ASSUMED TO BE NEGLIGIBLE.       
    THE MULTIPLE OF I THAT IS ADDED TO G IS AT MOST DIAGR TIMES 
    THE LEAST MULTIPLE OF I THAT GIVES POSITIVE DEFINITENESS.   
    X IS RE-INITIALISED IF ITS MAGNITUDE IS REDUCED BY THE      
    FACTOR XMAGR.                                               
    A CHECK IS MADE FOR AN INCREASE IN F EVERY IFINC ITERATIONS,  
    AFTER KFINC ITERATIONS ARE COMPLETED. */
 
    diagr = 2.0; 
    xmagr = 1.0e-2;
    ifinc = 3;
    kfinc = imax(10,n);
 
/*  FIND THE RECIPROCALS OF THE LENGTHS OF THE CONSTRAINT NORMALS. 
    RETURN IF A CONSTRAINT IS INFEASIBLE DUE TO A ZERO NORMAL. */
 
    *nact = 0;
    if (m <= 0)
        goto L45;

    for (k = 1; k <= m; ++k) {
        sum = 0.0;
        for (i = 1; i <= n; ++i) {
            sum = sum + pow(a[(k - 1) * n + i],2.0);             
        }
        if (sum > 0.0)
            goto L20;
        if (b[k] == 0.0)
            goto L30;
        info = -k;
        if (k <= meq)
            goto L730;

        if (b[k] <= 0.0)
            goto L30;
        else
            goto L730;
L20:
        sum = 1.0 / sqrt(sum);
L30:
        ia = iwa + k;
L40:
        w[ia] = sum;
    }
L45:
    for (k = 1; k <= n; ++k) {
        ia = iwa + m + k;
        w[ia] = 1.0;
    }
 
/* IF NECESSARY INCREASE THE DIAGONAL ELEMENTS OF G. */
 
    if (lql == 0)
        goto L165;
    
    *diag = 0.0;
    for (i = 1; i <= n; ++i) {
        id = iwd + i;
        w[id] = g[(i - 1) * n + i];           
        *diag = dmax(*diag,vsmall - w[id]);
        if (i == n)
            goto L60;
        ii = i + 1;
        for (j = ii; j <= n; ++j) {
            ga = -dmin(w[id],g[(j-1) * n + j]);
            gb = fabs(w[id]-g[(j - 1) * n + j]) + fabs(g[(i - 1) * n + j]);
            if (gb > 0.0)
                ga = ga + pow(g[(i - 1) * n + j],2.0) / gb;
            *diag = dmax(*diag,ga);
        }
L60:
        continue;
    }
    if (*diag <= 0.0)
        goto L90;

L70:
    *diag *= diagr;
    for (i = 1; i <= n; ++i) {
        id = iwd + i;
        g[(i - 1) * n + i] = *diag + w[id];
    }
 
/*  FORM THE CHOLESKY FACTORISATION OF G. THE TRANSPOSE   */
/*  OF THE FACTOR WILL BE PLACED IN THE R-PARTITION OF W. */
 
L90:
    ir = iwr;
    for (j = 1; j <= n; ++j) {
        ira = iwr;
        irb = ir + 1;
        for (i = 1; i <= j; ++i) {
                        
            temp = g[(i - 1) * n + j];
            if (i == 1)
                goto L110;

            for (k = irb; k <= ir; ++k) {
                ira++;        
                temp = temp - w[k] * w[ira];
            }
L110:
            ir++;    
            ira++;        
            if (i < j)
                w[ir] = temp / w[ira];
        }
        if (temp < vsmall)
            goto L140;
        w[ir] = sqrt(temp);
    }
    goto L170;
 
/*  INCREASE FURTHER THE DIAGONAL ELEMENT OF G. */
  
L140:
    w[j] = 1.0;
    sumx = 1.0;
    k = j;
L150:
    sum = 0.0;
    ira = ir - 1;

    for (i = k; i <= j; ++i) {
        sum = sum - w[ira] * w[i];
        ira += i; 
    }
    ir = ir - k;
    k--;  
    w[k] = sum / w[ir];
    sumx = sumx + w[k] * w[k];
    if (k >= 2)
        goto L150;

    *diag = *diag + vsmall - temp / sumx;
    goto L70;
 
/* STORE THE CHOLESKY FACTORISATION IN THE R-PARTITION OF W */
  
L165:
    ir = iwr;
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= i; ++j) {
            ir++;     
            w[ir] = g[(j - 1) * n + i];
        }
    }
 
/* SET Z THE INVERSE OF THE MATRIX IN R. */
  
L170:
    nm = n - 1;
    for (i = 1; i <= n; ++i) {
        iz = iwz + i;
        if (i == 1)
            goto L190;
        for (j = 2; j <= i; ++j) {
            w[iz] = 0.0;
            iz += n;
        }
L190:
        ir = iwr + (i + i*i) / 2;
        w[iz] = 1.0 / w[ir];
        if (i == n)
            continue;  
        iza = iz;
        for (j = i; j <= nm; ++j) {
            ir += i;
            sum = 0.0;
            for (k = iza; k <= iz; k += n) {
                sum = sum + w[k] * w[ir];
                ir++;    
            }
            iz += n;
            w[iz] = -sum / w[ir];
        }
    }
 
/* SET THE INITIAL VALUES OF SOME VARIABLES.                     
   ITERC COUNTS THE NUMBER OF ITERATIONS.                        
   ITREF IS SET TO ONE WHEN ITERATIVE REFINEMENT IS REQUIRED.    
   JFINC INDICATES WHEN TO TEST FOR AN INCREASE IN F. */
                                                                  
    iterc = 1;
    itref = 0;
    jfinc = -kfinc;
 
/* SET X TO ZERO AND SET THE CORRESPONDING RESIDUALS OF THE    
   KUHN-TUCKER CONDITIONS. */
 
L230:
    iflag = 1;
    iws = iww - n;

    for (i = 1; i <= n; ++i) {
        x[i] = 0.0;
        iw = iww + i;
        w[iw] = grad[i];
        if (i > *nact)
            goto L240;
        w[i] = 0.0;
        is = iws + i;
        k = iact[i];
        if (k <= m)
            goto L235;
        if (k > mn)
            goto L234;
        k1 = k - m;
        w[is] = xl[k1];
        goto L240;
L234:
        k1 = k - mn;
        w[is] = -xu[k1];
        goto L240;
L235:
        w[is] = b[k];
L240:   continue;
    }
    xmag = 0.0;
    vfact = 1.0; 
    if (*nact <= 0) 
        goto L340;
    else 
        goto L280;

    /* if (*nact) 340,340,280 */
  
/* SET THE RESIDUALS OF THE KUHN-TUCKER CONDITIONS FOR GENERAL X. */
  
L250:
    iflag = 2;
    iws = iww - n; 

    for (i = 1; i <= n; ++i) {
        iw = iww + i;
        w[iw] = grad[i];
        if (lql)
            goto L259;
        id = iwd + i;
        w[id] = 0.0;
        for (j = i; j <= n; ++j) {
            w[id] = w[id] + g[(i - 1) * n + j] * x[j];
        }
        for (j = 1; j <= i; ++j) {
            id = iwd + j;
            w[iw] = w[iw] + g[(j - 1) * n + i] * w[id];
        }
        goto L260;
L259:
        for (j = 1; j <= n; ++j) {
            w[iw] = w[iw] + g[(i - 1) * n + j] * x[j];
        }
L260:   continue;
    }
    if (*nact == 0)
        goto L340;

    for (k = 1; k <= *nact; ++k) {
        kk = iact[k];
        is = iws + k;
        if (kk > m)
            goto L265;
        w[is] = b[kk];

        for (i = 1; i <= n; ++i) {
            iw = iww + i; 
            w[iw] = w[iw] - w[k] * a[(kk - 1) * n + i];
            w[is] = w[is] - x[i] * a[(kk - 1) * n + i];
        }
        goto L270;
L265:
        if (kk > mn)
            goto L266;
        k1 = kk - m;
        iw = iww + k1;
        w[iw] = w[iw]-w[k];
        w[is] = xl[k1] - x[k1];
        goto L270;
L266:
        k1 = kk - mn;
        iw = iww + k1;
        w[iw] = w[iw] + w[k];
        w[is] = -xu[k1] + x[k1];
L270:   continue;   
    }

/*  PRE-MULTIPLY THE VECTOR IN THE S-PARTITION OF W BY THE  */
/*  INVERS OF R TRANSPOSE. */
  
L280:
    ir = iwr;
    ip = iww + 1;
    ipp = iww + n;
    il = iws + 1;
    iu = iws + *nact;

    for (i = il; i <= iu; ++i) {
        sum = 0.0;
        if (i == il)
            goto L300;
        ju = i - 1;
        for (j = il; j <= ju; ++j) {
            ir++;             
            sum = sum + w[ir] * w[j];
        }
L300:
        ir++;           
        w[i] = (w[i] - sum) / w[ir];
    }

/* SHIFT X TO SATISFY THE ACTIVE CONSTRAINTS AND MAKE THE */
/* CORRESPONDING CHANGE TO THE GRADIENT RESIDUALS. */
 
    for (i = 1; i <= n; ++i) {
        iz = iwz + i;
        sum = 0.0;
        for (j = il; j <= iu; ++j) {
            sum = sum + w[j] * w[iz];
            iz = iz + n;
        }
        x[i] = x[i] + sum;
        if (lql)
            goto L329;
        id = iwd + i;
        w[id] = 0.0;
        for (j = i; j <= n; ++j) {
            w[id] = w[id] + g[(i - 1) * n + j] * sum;
        }
        iw = iww + i;
        for (j = 1; j <= i; ++j) {
            id = iwd + j;
            w[iw] = w[iw] + g[(j - 1) * n + i] * w[id];
        }
        goto L330;
L329:
        for (j = 1; j <= n; ++j) {
            iw = iww + j;
            w[iw] = w[iw] + sum * g[(i - 1) * n + j];
        }
L330:
        continue; 
    }
   
/* FORM THE SCALAR PRODUCT OF THE CURRENT GRADIENT RESIDUALS */
/* WITH EACH COLUMN OF Z. */
  
L340:
    kflag = 1;
    goto L930;
L350:
    if (*nact == n)
        goto L380;
 
/* SHIFT X SO THAT IT SATISFIES THE REMAINING KUHN-TUCKER CONDITIONS. */
  
    il = iws + *nact + 1;
    iza = iwz + *nact * n;

    for (i = 1; i <= n; ++i) {
        sum = 0.0;
        iz = iza + i;
        for (j = il; j <= iww; ++j) {
            sum = sum + w[iz] * w[j];
            iz = iz + n;
        }
        x[i] = x[i] - sum;
    }
    info = 0;
    if (*nact == 0)
        goto L410;
 
/* UPDATE THE LAGRANGE MULTIPLIERS. */
  
L380:
    lflag = 3;
    goto L740;

L390:
    for (k = 1; k <= *nact; ++k) {
        iw = iww + k;
        w[k] = w[k] + w[iw]; 
    }
 
/* REVISE THE VALUES OF XMAG. */
/* BRANCH IF ITERATIVE REFINEMENT IS REQUIRED. */
 
L410:
    jflag = 1;
    goto L910;

L420:
    if (iflag == itref) 
        goto L250;
 
/* DELETE A CONSTRAINT IF A LAGRANGE MULTIPLIER OF AN */
/* INEQUALITY CONSTRAINT IS NEGATIVE. */
 
    kdrop = 0;
    goto L440;

L430:
    kdrop = kdrop + 1;
    if (w[kdrop] >= 0.0)
        goto L440;
    if (iact[kdrop] <= meq) 
        goto L440;
    nu = *nact; 
    mflag = 1;
    goto L800;

L440:
    if (kdrop < *nact)
        goto L430;
 
/* SEEK THE GREATEAST NORMALISED CONSTRAINT VIOLATION, DISREGARDING */
/* ANY THAT MAY BE DUE TO COMPUTER ROUNDING ERRORS. */
   
L450:
    cvmax = 0.0;
    if (m <= 0)
        goto L481;

    for (k = 1; k <= m; ++k) {
        ia = iwa + k;
        if (w[ia] <= 0.0)
            goto L480;
        sum = -b[k];

        for (i = 1; i <= n; ++i) {
            sum = sum + x[i] * a[(k - 1) * n + i];
        }
        sumx = -sum * w[ia];
        if (k <= meq)
            sumx = fabs(sumx);
        if (sumx <= cvmax)
            goto L480;
        temp = fabs(b[k]);

        for (i = 1; i <= n; ++i) {
            temp = temp + fabs(x[i] * a[(k - 1) * n + i]);
        }
        tempa = temp + fabs(sum);
        if (tempa <= temp)
            goto L480;
        temp = temp + 1.5 * fabs(sum);
        if (temp <= tempa)
            goto L480;
        cvmax = sumx;
        res = sum;
        knext = k;
L480:   continue;
    }

L481:
    for (k = 1; k <= n; ++k) {
        lower = 1;     
        ia = iwa + m + k;
        if (w[ia] <= 0.0)
            goto L485;
        sum = xl[k] - x[k];

        if (sum < 0.0)
            goto L482;
        else if (sum > 0.0)
            goto L483;
        else
            goto L485;

L482:
        sum = x[k] - xu[k];
        lower = 0;       
L483:
        if (sum <= cvmax)
            goto L485;
        cvmax = sum;
        res = -sum;
        knext = k + m;
        if (lower)
            goto L485;
        knext = k + mn;
L485:
        continue;
    }

/* TEST FOR CONVERGENCE */
 
    info = 0;
    if (cvmax <= vsmall)
        goto L700;
 
/* RETURN IF, DUE TO ROUNDING ERRORS, THE ACTUAL CHANGE IN */
/* X MAY NOT INCREASE THE OBJECTIVE FUNCTION */
 
    jfinc = jfinc + 1;
    if (jfinc == 0)
        goto L510;
    if (jfinc != ifinc)
        goto L530;
    fdiff = 0.0;
    fdiffa = 0.0;

    for (i = 1; i <= n; ++i) {
        sum = 2.0 * grad[i];
        sumx = fabs(sum);
        if (lql)
            goto L489;
        id = iwd + i;
        w[id] = 0.0;

        for (j = i; j <= n; ++j) {
            ix = iwx + j;
            w[id] = w[id] + g[(i - 1) * n + j] * (w[ix] + x[j]);
        }
        for (j = 1; j <= i; ++j) {
            id = iwd + j;
            temp = g[(j - 1) * n + i] * w[id];
            sum = sum + temp;
            sumx = sumx + fabs(temp);
        }
        goto L495;
L489:
        for (j = 1; j <= n; ++j) {
            ix = iwx + j; 
            temp = g[(i - 1) * n + j] * (w[ix] + x[j]);
            sum = sum + temp;
            sumx = sumx + fabs(temp);
        }
L495:
        ix = iwx + i;
        fdiff = fdiff + sum * (x[i] - w[ix]);
L500:
        fdiffa = fdiffa + sumx * fabs(x[i] - w[ix]);
    }

    info = 2;
    sum = fdiffa + fdiff;
    if (sum <= fdiffa)
        goto L700;
    temp = fdiffa + 1.5 * fdiff;
    if (temp <= sum)
        goto L700;
    jfinc = 0;
    info = 0;
   
L510:
    for (i = 1; i <= n; ++i) {
        ix = iwx + i;
        w[ix] = x[i];
    }
 
/* FORM THE SCALAR PRODUCT OF THE NEW CONSTRAINT NORMAL WITH EACH */
/* COLUMN OF Z. PARNEW WILL BECOME THE LAGRANGE MULTIPLIER OF */
/* THE NEW CONSTRAINT. */
  
L530:
    iterc = iterc + 1;
    if (iterc <= maxit)
        goto L531;

    info = 1;
    goto L710;
L531:
    iws = iwr + (*nact + *nact * *nact) / 2;
    if (knext > m)
        goto L541;

    for (i = 1; i <= n; ++i) {
        iw = iww + i;
        w[iw] = a[(knext - 1) * n + i];
    }
    goto L549;

L541:
    for (i = 1; i <= n; ++i) {
        iw = iww + i;
        w[iw] = 0.0;
    }
    k1 = knext - m;
    if (k1 > n)
        goto L545;

    iw = iww + k1;
    w[iw] = 1.0;
    iz = iwz + k1;

    for (i = 1; i <= n; ++i) {
        is = iws + i;
        w[is] = w[iz];
        iz = iz + n;
    }
    goto L550;
   
L545:
    k1 = knext - mn;
    iw = iww + k1;
    w[iw] = -1.0;
    iz = iwz + k1;

    for (i = 1; i <= n; ++i) {
        is = iws + i;
        w[is] = -w[iz];
        iz = iz + n;
    }
    goto L550;

L549:
    kflag = 2;
    goto L930;

L550:
    parnew = 0.0;
 
/* APPLY GIVENS ROTATIONS TO MAKE THE LAST (N-NACT-2) SCALAR */
/* PRODUCTS EQUAL TO ZERO. */
 
    if (*nact == n)
        goto L570;
    nu = n;
    nflag = 1;
    goto L860;
 
/* BRANCH IF THERE IS NO NEED TO DELETE A CONSTRAINT. */
   
L560:
    is = iws + *nact;
    if (*nact == 0)
        goto L640;

    suma = 0.0;
    sumb = 0.0;
    sumc = 0.0;
    iz = iwz + *nact * n;

    for (i = 1; i <= n; ++i) {
        iz = iz + 1;
        iw = iww + i;
        suma = suma + w[iw] * w[iz];
        sumb = sumb + fabs(w[iw] * w[iz]);
        sumc = sumc + w[iz] * w[iz];
    }
    temp = sumb + 0.1 * fabs(suma);
    tempa = sumb + 0.2 * fabs(suma);
    if (temp <= sumb)
        goto L570;
    if (tempa <= temp)
        goto L570;
    if (sumb > vsmall)
        goto L5;
    goto L570;

L5:
    sumc = sqrt(sumc);
    ia = iwa + knext;
    if (knext <= m)
        sumc = sumc / w[ia];
    temp = sumc + 0.1 * fabs(suma);
    tempa = sumc + 0.2 * fabs(suma);
    if (temp <= sumc)
        goto L567;
    if (tempa <= temp)
        goto L567;
    goto L640;
 
/* CALCULATE THE MULTIPLIERS FOR THE NEW CONSTRAINT NORMAL
   EXPRESSED IN TERMS OF THE ACTIVE CONSTRAINT NORMALS.
   THEN WORK OUT WHICH CONTRAINT TO DROP. */
 
L567:
    lflag = 4;
    goto L740;
L570:
    lflag = 1;
    goto L740;
 
/* COMPLETE THE TEST FOR LINEARLY DEPENDENT CONSTRAINTS. */
 
L571:
    if (knext > m)
        goto L574;
               
    for (i = 1; i <= n; ++i) {
        suma = a[(knext - 1) * n + i];
        sumb = fabs(suma);
        if (*nact == 0)
            goto L581;

        for (k = 1; k <= *nact; ++k) {
            kk = iact[k];
            if (kk <= m)
                goto L568;
            kk = kk - m;
            temp = 0.0;
            if (kk == i)
                temp = w[iww + kk];
            kk = kk - n;
            if (kk == i)
                temp = -w[iww + kk];
            goto L569;
L568:       continue;
            iw = iww + k;
            temp = w[iw] * a[(kk - 1) * n + i];
L569:       continue;
            suma = suma - temp;
            sumb = sumb + fabs(temp);
        }
L581:
        if (suma <= vsmall)
            goto L573;

        temp = sumb + 0.1 * fabs(suma);
        tempa = sumb + 0.2 * fabs(suma);
        if (temp <= sumb)
            goto L573;
        if (tempa <= temp)
            goto L573;
        goto L630;
L573:   continue;
    }
    lflag = 1;
    goto L775;

L574:
    k1 = knext - m;
    if (k1  > n)
        k1 = k1 - n;

    for (i = 1; i <= n; ++i) {
        suma = 0.0;
        if (i != k1)
            goto L575;
        suma = 1.0;
        if (knext > mn)
            suma = -1.0;
L575:   sumb = fabs(suma);
        if (*nact == 0)
            goto L582;

        for (k = 1; k <= *nact; ++k) {
            kk = iact[k];
            if (kk <= m)
                goto L579;
            kk = kk - m;
            temp = 0.0;
            if (kk == i)
                temp = w[iww + kk];
            kk = kk - n;
            if (kk == i)
                temp = -w[iww + kk];
            goto L576;
L579:       iw = iww + k;
            temp = w[iw] * a[(kk - 1) * n + i];
L576:       suma = suma - temp;
L577:       sumb = sumb + fabs(temp);
        }
L582:
        temp = sumb + 0.1 * fabs(suma);
        tempa = sumb + 0.2 * fabs(suma);
        if (temp <= sumb)
            goto L578;
        if (tempa <= temp)
            goto L578;
        goto L630;
L578:   continue;
    }
    lflag = 1;
    goto L775;
   
/* BRANCH IF THE CONTRAINTS ARE INCONSISTENT. */
  
L580:
    info = -knext;
    if (kdrop == 0)
        goto L700;
    parinc = ratio;
    parnew = parinc;
 
/* REVISE THE LAGRANGE MULTIPLIERS OF THE ACTIVE CONSTRAINTS. */
 
L590:
    if (*nact == 0)
        goto L601;

    for (k = 1; k <= *nact; ++k) {
        iw = iww + k;
        w[k] = w[k] - parinc * w[iw];
        if (iact[k] > meq)
            w[k] = dmax(0.0,w[k]);
    }
L601:
    if (kdrop  == 0)
        goto L680;
 
/* DELETE THE CONSTRAINT TO BE DROPPED.
   SHIFT THE VECTOR OF SCALAR PRODUCTS.
   THEN, IF APPROPRIATE, MAKE ONE MORE SCALAR PRODUCT ZERO. */
 
    nu = *nact + 1;
    mflag = 2;
    goto L800;
L610:
    iws = iws - *nact - 1;
    nu = imin(n,nu);

    for (i = 1; i <= nu; ++i) {
        is = iws + i;
        j = is + *nact;
        w[is] = w[j + 1];
    }
    nflag = 2;
    goto L860;
 
/* CALCULATE THE STEP TO THE VIOLATED CONSTRAINT. */
 
L630:
    is = iws + *nact;
L640:
    sumy = w[is + 1];
    step = -res / sumy;
    parinc = step / sumy;
    if (*nact == 0)
        goto L660;
   
/* CALCULATE THE CHANGES TO THE LAGRANGE MULTIPLIERS, AND REDUCE
   THE STEP ALONG THE NEW SEARCH DIRECTION IF NECESSARY. */
 
    lflag = 2;
    goto L740;

L650:
    if (kdrop  == 0)
        goto L660;
    temp = 1.0 - ratio / parinc;
    if (temp <= 0.0)
        kdrop = 0;
    if (kdrop == 0)
        goto L660;
    step = ratio * sumy;
    parinc = ratio;
    res = temp * res;
 
/* UPDATE X AND THE LAGRANGE MULTIPIERS.
   DROP A CONSTRAINT IF THE FULL STEP IS NOT TAKEN. */
 
L660:
    iwy = iwz + *nact * n;

    for (i = 1; i <= n; ++i) {
        iy = iwy + i;
        x[i] = x[i] + step * w[iy];
    }
    parnew = parnew + parinc;
    if (*nact >= 1)
        goto L590;
 
/* ADD THE NEW CONSTRAINT TO THE ACTIVE SET. */
 
L680:
    *nact = *nact + 1;
    w[*nact] = parnew;
    iact[*nact] = knext;
    ia = iwa + knext;
    if (knext > mn)
        ia = ia - n;
    w[ia] = -w[ia];
 
/* ESTIMATE THE MAGNITUDE OF X. THEN BEGIN A NEW ITERATION,
   RE-INITILISING X IF THIS MAGNITUDE IS SMALL. */
 
    jflag = 2;
    goto L910;

L690:
    if (sum < (xmagr * xmag))
        goto L230;
    /* if (itref) 450,450,250 */

    if (itref <= 0)
        goto L450;
    else
        goto L250;
 
/* INITIATE ITERATIVE REFINEMENT IF IT HAS NOT YET BEEN USED,
   OR RETURN AFTER RESTORING THE DIAGONAL ELEMENTS OF G. */
 
L700:
    if (iterc == 0) 
        goto L710;
    itref = itref + 1;
    jfinc = -1;
    if (itref == 1)
        goto L250;
L710:
    if (lql == 0)
        goto LFIN;   

    for (i = 1; i <= n; ++i) {
        id = iwd + i;
        g[(i - 1) * n + i] = w[id];
    }
L730:
    goto LFIN;   
 
/*  THE REMAINIG INSTRUCTIONS ARE USED AS SUBROUTINES. */
 
/*  CALCULATE THE LAGRANGE MULTIPLIERS BY PRE-MULTIPLYING THE  */
/*  VECTOR IN THE S-PARTITION OF W BY THE INVERSE OF R.        */
   
L740:
    ir = iwr + (*nact + *nact * *nact) / 2;
    i = *nact;
    sum = 0.0;
    goto L770;

L750:
    ira = ir - 1;
    sum = 0.0;
    if (*nact == 0)
        goto L761;

    for (j = i; j <= *nact; ++j) {
        iw = iww + j;
        sum = sum + w[ira] * w[iw];
        ira += j;
    }
L761:
    ir -= i;
    i--;    
L770:
    iw = iww + i;
    is = iws + i;
    w[iw] = (w[is] - sum) / w[ir];
    if (i > 1)
        goto L750;
    if (lflag == 3)
        goto L390;
    if (lflag == 4)
        goto L571;
 
/* CALCULATE THE NEXT CONSTRAINT TO DROP. */
 
L775:
    ip = iww + 1;
    ipp = iww + *nact;
    kdrop = 0;
    if (*nact == 0)
        goto L791;

    for (k = 1; k <= *nact; ++k) {
        if (iact[k] <= meq)
            goto L790;
        iw = iww + k;
        if ((res * w[iw]) >= 0.0)
            goto L790;
        temp = w[k] / w[iw];
        if (kdrop == 0)
            goto L780; 
        if (fabs(temp) >= fabs(ratio))
            goto L790;
L780:
        kdrop = k;
        ratio = temp;
L790:
        continue;
    }
L791:
    if (lflag == 1)
        goto L580;
    else
        goto L650;
 
/*  DROP THE CONSTRAINT IN POSITION KDROP IN THE ACTIVE SET. */
 
L800:
    ia = iwa + iact[kdrop];
    if (iact[kdrop] > mn)
        ia -= n;   
    w[ia] = -w[ia];
    if (kdrop == *nact)
        goto L850;
 
/* SET SOME INDICES AND CALCULATE THE ELEMENTS OF THE NEXT GIVENS ROTATION.*/
 
    iz = iwz + kdrop * n;
    ir = iwr + (kdrop + kdrop * kdrop) / 2;
L810:
    ira = ir;
    ir = ir + kdrop + 1;
    temp = dmax(fabs(w[ir - 1]),fabs(w[ir]));
    sum = temp * sqrt(pow((w[ir - 1] / temp),2.0) + pow((w[ir] / temp),2.0)); 
    ga = w[ir - 1] / sum;
    gb = w[ir] / sum;
 
/* EXCHANGE THE COLUMNS OF R. */
 
    for (i = 1; i <= kdrop; ++i) {
        ira++;            
        j = ira - kdrop;
        temp = w[ira];
        w[ira] = w[j];
        w[j] = temp;
    }
    w[ir] = 0.0;
 
/* APPLY THE ROTATION TO THE ROWS OF R. */
 
    w[j] = sum;
    kdrop++;             
    for (i = kdrop; i <= nu; ++i) {
        temp = ga * w[ira] + gb * w[ira + 1];
        w[ira + 1] = ga * w[ira + 1] - gb * w[ira];
        w[ira] = temp;
        ira += i;           
    }

/* APPLY THE ROTATION TO THE COLUMNS OF Z. */
 
    for (i = 1; i <= n; ++i) {
        iz++;           
        j = iz - n;
        temp = ga * w[j] + gb * w[iz];
        w[iz] = ga * w[iz] - gb * w[j];
        w[j] = temp;
    }
 
/*  REVISE iact AND THE LAGRANGE MULTIPLIERS. */
 
    iact[kdrop - 1] = iact[kdrop];
    w[kdrop - 1] = w[kdrop];
    if (kdrop < *nact)
        goto L810;
L850:
    *nact -= 1;        
    if (mflag == 1)
        goto L250;
    else
        goto L610;
 
/*  APPLY GIVENS ROTATION TO REDUCE SOME OF THE SCALAR  */
/*  PRODUCTS IN THE S-PARTITION OF W TO ZERO.           */
 
L860:
    iz = iwz + nu * n;
L870:
    iz = iz - n;
L880:
    is = iws + nu;
    nu--;         
    if (nu == *nact)
        goto L900;
    if (w[is] == 0.0)
        goto L870;

    temp = dmax(fabs(w[is - 1]),fabs(w[is]));
    sum = temp * sqrt(pow((w[is - 1] / temp),2.0) + pow((w[is] / temp),2.0));
    ga = w[is - 1] / sum;
    gb = w[is] / sum;
    w[is - 1] = sum;

    for (i = 1; i <= n; ++i) {
        k = iz + n;
        temp = ga * w[iz] + gb * w[k];
        w[k] = ga * w[k] - gb * w[iz];
        w[iz] = temp;
        iz--;        
    }
    goto L880;
L900:
    if (nflag == 1)
        goto L560;
    else
        goto L630;
 
/* CALCULATE THE MAGNITUDE OF X AN REVISE XMAG. */
   
L910:
    sum = 0.0;
    for (i  =  1; i <= n; ++i) {
        sum = sum + fabs(x[i]) * vfact * (fabs(grad[i]) + fabs(g[(i-1)*n+i] * x[i]));
        if (lql)
            goto L920;
        if (sum < 1.e-30)
            goto L920;
        vfact = 1.e-10 * vfact;
        sum = 1.e-10 * sum;
        xmag = 1.e-10 * xmag;
L920:   continue;
    }
L925:
    xmag = dmax(xmag,sum);
    if (jflag == 1)
        goto L420;
    else
        goto L690;
 
/* PRE-MULTIPLY THE VECTOR IN THE W-PARTITION OF W BY Z TRANSPOSE. */
   
L930:
    jl = iww + 1;
    iz = iwz; 
    for (i = 1; i <= n; ++i) {
        is = iws + i;
        w[is] = 0.0;
        iwwn = iww + n;
        for (j = jl; j <= iwwn; ++j) {
            iz++;            
            w[is] = w[is] + w[iz] * w[j];
        }
    }
    if (kflag == 1)
        goto L350;
    else
        goto L550;

LFIN:
    free((char *)w);
    memrq(-lw,sizeof(double));
    return(info);
}            

