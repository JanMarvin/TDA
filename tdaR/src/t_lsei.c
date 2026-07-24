/****************************************************************************/
/*  t_lsei                                                                  */
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
#include "tda_context.h"

/*  functions in t_lsei.c */

int lsei(TDAContext *ctx, double *w, int me, int ma, int mg, int n, double *x, int cov, double *rnorme, double *rnorml, int *ranka, int *ranke);
int lsi(TDAContext *ctx, double *w, int nw, int ma, int mg, int n, double *x, double *rnorm, int *ranka, int cov);
int lpdp(TDAContext *ctx, double *a, int na, int m, int n1, int n2, double *x, double *wnorm);
int wnnls(TDAContext *ctx, double *w, int nw, int me, int ma, int n, int l, double *x, double *rnorm);
int wnlsm(TDAContext *ctx, double *w, int nw, int mme, int ma, int n, int l, double *x, double *rnorm, int *ipivot, int *itype, double *wd, double *h, double *scale, double *z, double *temp, double *d);
int wnlit(TDAContext *ctx, double *w, int nw, int m, int n, int l, int *ipivot, int *itype, double *h, double *scale, double *rnorm, int *idope, double *dope, int done);
void srotmg(TDAContext *ctx, double *sd1, double *sd2, double *sx1, double *sy1, double *sparam);
void srotm(TDAContext *ctx, int n, double *sx, double *sy, double *sparam);
int lhhfti(TDAContext *ctx, int m, int n, int na, int nb, double *a, double *b, double *rnorm, int mode, double *d, double *h, int *ip, double tau);
void lhh12(TDAContext *ctx, int mode, int lpivot, int l1, int m, double *u, int iue, double *up, double *c, int ice, int icv, int ncv);
int syminv(TDAContext *ctx, int n, double *x);

/* ------------------------------------------------------------------------ */
/*  lsei                                                                    */
/*                                                                          */
/*      Solves a least squares problem with linear constraints.             */
/*      The code is adapted from: R.J. Hanson, K.H. Haskell, Algorithm 587  */
/*      from Collected Algorithms of ACM.                                   */
/*                                                                          */
/*      The problem solved by this function is as follows:                  */
/*                                                                          */
/*          Ex =  F     E = (me,n)-matrix, F = (me)-vector                  */
/*          Ax =  B     A = (ma,n)-matrix, B = (ma)-vector                  */
/*          Gx >= H     G = (mg,n)-matrix, H = (mg)-vector                  */
/*                                                                          */
/*      Where Ex = F are equality constraints, Ax = B is solved in the      */
/*      least squares sense, and Gx >= H are inequality constraints.        */
/*                                                                          */
/*      Possible dimensions: me >= 0, ma >= 0, mg >= 0, n >= 0              */
/*                                                                          */
/*      Call: int lsei(w,me,ma,mg,n,x,cov,&rnorme,&rnorml,&ranka,&ranke)    */
/*                                                                          */
/*      Input data are expected in a two-dimensional array w(i,j), stored   */
/*      in a one-dimensional array by w(i,j) = w[(i - 1) * nw + j] with     */
/*      nw = n + 1, as follows:                                             */
/*                                                                          */
/*              (E,F)   me rows )                                           */
/*          w = (A,B)   ma rows )   m = me + ma + mg rows                   */
/*              (G,H)   mg rows )                                           */
/*                                                                          */
/*      The solution, if found, is returned in the (n)-vector x. On         */
/*      successful return, rnorme contains the Euclidean norm of the        */
/*      residuals F - Ex, and rnorml contains the norm of the residuals     */
/*      B - Ax; also ranka and ranke contain the rank of the A-matrix and   */
/*      of the E-matrix, respectively.                                      */
/*                                                                          */
/*      cov is a logical flag. If set, the function tries to calculate      */
/*      the covariance matrix. If successful, the covariance matrix is      */
/*      returned in w(i,j); i,j = 1,n.                                      */
/*                                                                          */
/*      Note: for covariance matrix calculation, the array w[] must have    */
/*      at least n rows!                                                    */
/*                                                                          */
/*      Return codes:     = 0   Successful operation.                       */
/*                         -1   Equality constraints are contradictory      */
/*                         -2   Inequality constraints are incompatible,    */
/*                         -3   Equality and inequality constraints are     */
/*                              contradictory.                              */
/*                         -4   Insufficient memory.                        */
/*                                                                          */
/*      Note: For calculation of tolerances, the square root of the         */
/*      machine precision (EPSI) is used.                                   */
/*                                                                          */
/*      Note: the input data are always scaled such that the norm of        */
/*      the columns in w[] becomes 1.                                       */
/*                                                                          */

int lsei(TDAContext *ctx, double *w,int me,int ma,int mg,int n,double *x,int cov, double *rnorme,double *rnorml,int *ranka,int *ranke)
{
    register int i,j,k,l,ll;
    int ierr,m,kranke,n1,nw,imaxx,jp1,mapke1,mend,wsn;
    double tmp,enorm,fnorm,xnorm,xnrme,snmax,rnmax,sn,rn,up,tau,tau2;
    double gam,uj,vj,rb,size,*ws;

    tau = sqrt(ctx->EPSI);
    tau2 = tau * tau;

    *rnorme = 0.0;
    *rnorml = 0.0;
    *ranka = 0;
    *ranke = 0;
    nw = n + 1;
    m = ma + me + mg;

    if (cov && m < n)  
        gerr_exit(ctx, 43);

    ierr = 0;
    if (n <= 0 || m <= 0)
        return(ierr);

    kranke = n;
    if (me < n)
        kranke = me;
    n1 = 2 * kranke + 1;

    /*  allocate ws for temporary storage */

    wsn = n + n1 + 1;
    if (!(ws = (double *) calloc((size_t)(wsn),sizeof(double)))) {
        ierr = -4;
        goto LSEIFin;
    }
    memrq(ctx, wsn,sizeof(double));

    for (j = 1; j <= n; ++j) {      /* scaling */
        sn = 0.0;
        for (i = 1; i <= m; ++i) {
            tmp = w[(i - 1) * nw + j];         
            sn += tmp * tmp;
        }
        if (sn > 0) {
            sn = 1.0 / sqrt(sn);
            for (i = 1; i <= m; ++i)
                w[(i - 1) * nw + j] *= sn;   
        }
        else
            sn = 1.0; 

        ws[n1 + j - 1] = sn;                
    }

    /*  compute norm of equality constraints and right side */

    l = kranke;
    enorm = 0.0;
    for (j = 1; j <= n; ++j) {
        tmp = 0.0;
        for (ll = 1; ll <= me; ++ll)
            tmp += fabs(w[(ll - 1) * nw + j]);
        if (enorm < tmp)
            enorm = tmp;
    }
    fnorm = 0;
    for (ll = 1; ll <= me; ++ll)
        fnorm += fabs(w[(ll - 1) * nw + nw]);

    /*  compute maximum ratio of vector lengths, partition at col i */

    if (l > 0) {
        snmax = rnmax = 0.0;
        imaxx = 1;
        for (i = 1; i <= l; ++i) {
            for (k = i; k <= me; ++k) {
                sn = 0.0;
                for (ll = i; ll <= n; ++ll)
                    sn += w[(k - 1) * nw + ll] * w[(k - 1) * nw + ll];
                rn = 0.0;
                for (ll = 1; ll < i; ++ll)
                    rn += w[(k - 1) * nw + ll] * w[(k - 1) * nw + ll];

                if (rn == 0.0 && sn > snmax) {
                    snmax = sn;
                    imaxx = k;
                }
                else if (k == i || sn * rnmax > rn * snmax) {
                    snmax = sn;
                    rnmax = rn;
                    imaxx = k;
                }
            }
            if (i != imaxx) {
                for (ll = 1; ll <= nw; ++ll) {
                    tmp = w[(i - 1) * nw + ll];
                    w[(i - 1) * nw + ll] = w[(imaxx - 1) * nw + ll];
                    w[(imaxx - 1) * nw + ll] = tmp;
                }
            }
            if (snmax > tau2 * rnmax) {
                lhh12(ctx, 1,i,i+1,n,w+(i-1)*nw,1,ws+i,w+i*nw,1,nw,m-i);
            }
            else {
                kranke = i - 1;
                break;
            }
        }
    }
    *ranke = kranke;        /* save rank of E matrix */

    /*  save diagonal elements */

    for (ll = 1; ll <= kranke; ++ll)
        ws[kranke + ll] = w[(ll - 1) * nw + ll];

    /*  apply householder transformations */

    if (kranke > 0 && kranke < me) {
        for (k = kranke; k >= 1; --k) {
            lhh12(ctx, 1,k,kranke+1,me,w+k-1,nw,&up,w,nw,1,k-1);
            lhh12(ctx, 2,k,kranke+1,me,w+k-1,nw,&up,w+nw-1,nw,1,1);
        }
    }
    if (kranke > 0) {
        for (ll = 1; ll <= kranke; ++ll)    {
            x[ll] = w[(ll - 1) * nw + nw];
        }
        for ( i = 1; i <= kranke; ++i) {
            tmp = 0.0;
            for (ll = 1; ll < i; ++ll)
                tmp += w[(i - 1) * nw + ll] * x[ll];
            x[i] = (x[i] - tmp) / w[(i - 1) * nw + i];
        }
    }

    /*  compute residuals for reduced problem */

    if (me < m) {
        for (i = me + 1; i <= m; ++i) {
            tmp = 0.0;
            for (ll = 1; ll <= kranke; ++ll)
                tmp += w[(i - 1) * nw + ll] * x[ll];
            w[(i - 1) * nw + nw] -= tmp;

            sn = 0.0;
            for (ll = 1; ll <= kranke; ++ll)
                sn += w[(i - 1) * nw + ll] * w[(i - 1) * nw + ll];

            rn = 0.0;
            for (ll = kranke + 1; ll <= n; ++ll)
                rn += w[(i - 1) * nw + ll] * w[(i - 1) * nw + ll];

            if (rn <= tau2 * sn && kranke < n) {
                for (ll = kranke + 1; ll <= n; ++ll)
                    w[(i - 1) * nw + ll] = 0.0;                                 
            }
        }
    }

    /*  norm of equality constraints */

    tmp = 0.0;
    for (ll = kranke + 1; ll <= me; ++ll)
        tmp += w[(ll - 1) * nw + nw] * w[(ll - 1) * nw + nw];
    *rnorme = sqrt(tmp);

    /*  move reduced problem data upwards if kranke less than me */

    if (kranke < me) {
        for (j = 1; j <= nw; ++j) {
            for (ll = 0; ll < m - me; ++ll)
                w[(kranke + ll) * nw + j] = w[(me + ll) * nw + j];
        }
    }
    
    /*  call lsi to solve the reduced problem */

    ierr = lsi(ctx, w + kranke * nw + kranke,nw,ma,mg,
                                      n - kranke,x + kranke,rnorml,ranka,cov);
    if (ierr)
        goto LSEIFin1;
   
    /*  test for consistency of equality constraints */

    if (me > 0) {
        xnrme = 0.0;
        for (ll = 1; ll <= kranke; ++ll) 
            xnrme += fabs(w[(ll - 1) * nw + nw]);

        if (*rnorme > tau * (enorm * xnrme + fnorm)) {
            ierr = -1;
        }
        if (kranke == n && mg > 0) {
            xnorm = 0.0;
            for (ll = 1; ll <= n; ++ll)
                xnorm += fabs(x[ll]);

            mapke1 = ma + kranke + 1;
            mend = ma + kranke + mg;

            for (i = mapke1; i <= mend; ++i) {
                tmp = 0.0;
                for (ll = 1; ll <= n; ++ll)
                    tmp += fabs(w[(i - 1) * nw + ll]);
                size = tmp * xnorm + fabs(w[(i - 1) * nw + nw]);

                if ((w[(i - 1) * nw + nw] > tau * size)) {
                    if (ierr == -1)
                        ierr = -3;
                    else
                        ierr = -2;
                    break;
                }
            }
        }
    }
    if (kranke > 0) {

        /*  replace diagonal elements */

        for (ll = 1; ll <= kranke; ++ll)
            w[(ll - 1) * nw + ll] = ws[kranke + ll];

        /*  reapply transformations to put solution in orig coordinates */

        for (i = kranke; i >= 1; --i) {
            lhh12(ctx, 2,i,i+1,n,w+(i-1)*nw,1,ws+i,x,1,1,1);
        }

        /*  adjust covariance matrix for equality constraints */

        if (cov) {
            for (j = kranke; j >= 1; --j) {
                if (j < n) {
                    rb = ws[j] * w[(j - 1) * nw + j];
                    if (rb != 0.0)
                        rb = 1.0 / rb;
                    jp1 = j + 1;

                    for (i = jp1; i <= n; ++i) {
                        tmp = 0.0;
                        for (ll = 0; ll < n - j; ++ll)
                            tmp += w[(i - 1) * nw + jp1 + ll] *
                                                     w[(j - 1) * nw + jp1 + ll];
                        w[(i - 1) * nw + j] = tmp * rb;
                    }
                    tmp = 0.0;
                    for (ll = 0; ll < n - j; ++ll)
                        tmp += w[(jp1 + ll - 1) * nw + j] *  
                                                     w[(j - 1) * nw + jp1 + ll];
                    gam = tmp * rb / 2.0;

                    for (ll = 0; ll < n - j; ++ll)
                        w[(jp1 + ll - 1) * nw + j] +=
                            gam * w[(j - 1) * nw + jp1 + ll];

                    for (i = jp1; i <= n; ++i) {
                        for (k = i; k <= n; ++k) {
                            w[(i - 1) * nw + k] +=
                                 w[(j - 1) * nw + i] * w[(k - 1) * nw + j] +
                                 w[(i - 1) * nw + j] * w[(j - 1) * nw + k];
                            w[(k - 1) * nw + i] = w[(i - 1) * nw + k];   
                        }
                    }
                    uj = ws[j];
                    vj = gam * uj;
                    w[(j - 1) * nw + j] = uj * vj + uj * vj;
                
                    for (i = jp1; i <= n; ++i)  {
                        w[(j - 1) * nw + i] = uj * w[(i - 1) * nw + j] +
                                              vj * w[(j - 1) * nw + i];
                    }
                    for (ll = 0; ll < n - j; ++ll)
                        w[(jp1 + ll - 1) * nw + j] = w[(j - 1) * nw + jp1 + ll];
                }
            }
        }
    }

    if (cov) {         /* rescale covariance matrix */ 
        for (i = 1; i <= n; ++i) {
            l = n1 + i;
            for (ll = 1; ll <= n; ++ll)
                w[(i - 1) * nw + ll] *= ws[l - 1];
            for (ll = 1; ll <= n; ++ll)
                w[(ll - 1) * nw + i] *= ws[l - 1];
        }
    }
    if (!ierr || ierr == -1) {                  /*  rescale the solution vector */
        for (j = 1; j <= n; ++j)  
            x[j] *= ws[n1 + j - 1];
    }

LSEIFin1:
    free((char *)ws);
    memrq(ctx, -wsn,sizeof(double));
LSEIFin:
    return(ierr);
}

/* ------------------------------------------------------------------------ */
/*  lsi                                                                     */
/*                                                                          */
/*      Solves a least squares problem with linear inequality constraints.  */
/*      The code is adapted from: R.J. Hanson, K.H. Haskell, Algorithm 587  */
/*      from Collected Algorithms of ACM.                                   */
/*                                                                          */
/*      The problem solved by this function is as follows:                  */
/*                                                                          */
/*          Ax =  B     A = (ma,n)-matrix, B = (ma)-vector                  */
/*          Gx >= H     G = (mg,n)-matrix, H = (mg)-vector                  */
/*                                                                          */
/*      Where Ex = F are equality constraints, Ax = B is solved in the      */
/*      least squares sense, and Gx >= H are inequality constraints.        */
/*                                                                          */
/*      Possible dimensions: ma >= 0, mg >= 0, n >= 0                       */
/*                                                                          */
/*      Call: int lsi(w,nw,ma,mg,n,x,&rnorm,&ranka,cov)                     */
/*                                                                          */
/*      Input data are expected in a two-dimensional array w(i,j), stored   */
/*      in a one-dimensional array by w(i,j) = w[(i - 1) * nw + j] with     */
/*      nw = number of columns in w (not neccessarily equal to n + 1) as    */
/*      follows:                                                            */
/*                                                                          */
/*          w = (A,B)   ma rows )   m = ma + mg rows                        */
/*              (G,H)   mg rows )                                           */
/*                                                                          */
/*      The vectors B and H are to be found in the n+1.th column of w.      */
/*                                                                          */
/*      The solution, if found, is returned in the (n)-vector x. On         */
/*      successful return, rnorm contains the Euclidean norm of the         */
/*      residuals B - Ax. The rank of A is returned in ranka.               */
/*                                                                          */
/*      cov is a logical flag. If set, the function tries to calculate      */
/*      the covariance matrix. If successful, the covariance matrix is      */
/*      returned in w(i,j); i,j = 1,n.                                      */
/*                                                                          */
/*      Note: for covariance matrix calculation, the array w[] must have    */
/*      at least n rows!                                                    */
/*                                                                          */
/*      Return codes:       0   Successful operation.                       */
/*                         -2   No success. Inequality constraints are      */
/*                              incompatible.                               */
/*                         -4   Insufficient memory.                        */
/*                                                                          */
/*      Note: For calculation of tolerances, an externally defined          */
/*      EPSI is used.                                                       */
    
int lsi(TDAContext *ctx, double *w,int nw,int ma,int mg,int n,double *x,double *rnorm, int *ranka,int cov)
{
    register int i,j,k,ll,ii,l;
    int lw,il,np1,ierr,m,krank,minman,n1,n2,n3,krm1,krp1,map1;
    int *ip;
    double rb,tau,tol,anorm,xnorm,fac,gam,tmp,tmp1,*ws;

    m = ma + mg;
    ierr = 0;
    *rnorm = 0.0;
    m = ma + mg;
    np1 = n + 1; 
    krank = 0;
    *ranka = 0;
    if (n <= 0 || m <= 0)
        return(ierr);

    /*  allocate ip to store permutation indices */

    if (!(ip = (int *) calloc((size_t)(n + 1),sizeof(int)))) {
        ierr = -4;
        goto LSIFin;
    }
    memrq(ctx, n + 1,sizeof(int));

    /*  compute matrix norm of least squares equations */

    anorm = 0.0;
    for (j = 1; j <= n; ++j) {
        tmp = 0.0;
        for (ll = 1; ll <= ma; ++ll)
            tmp += fabs(w[(ll - 1) * nw + j]);
        if (anorm < tmp)
            anorm = tmp;
    }
    tol = sqrt(ctx->EPSI);
    tau = tol * anorm;

    k = m;
    if (n > m)
        k = n;

    minman = ma;
    if (n < ma)
        minman = n;

    n1 = k + 1;
    n2 = n1 + n;

    /*  allocate ws for temporary storage */

    lw = n2 + 2 * k + 1;

    if (!(ws = (double *) calloc((size_t)(lw + 1),sizeof(double)))) {
        ierr = -4;
        goto LSIFin1;
    }
    memrq(ctx, lw + 1,sizeof(double));

    /*  compute Householder orthogonal decomposition of matrix */

    for (ll = 1; ll <= n; ++ll)
        ws[ll] = 0.0;

    for (ll = 1; ll <= ma; ++ll)
        ws[ll] = w[(ll - 1) * nw + np1];
        
    krank = lhhfti(ctx, ma,n,nw,1,w,ws,rnorm,1, ws + n2 - 1,ws + n1 - 1,ip,tau);

    *ranka = krank;
    fac = 1.0;
    gam = (double)(ma - krank);
    if (krank < ma)
        fac = *rnorm * *rnorm / gam;

    /*  reduce problem to be solved with lpdp */

    map1 = ma + 1;
    if (ma < m) {
        if (minman > 0) {
            for (i = map1; i <= m; ++i) {
                tmp = 0.0;
                for (ll = 1; ll <= n; ++ll)
                    tmp += w[(i - 1) * nw + ll] * ws[ll];
                w[(i - 1) * nw + np1] -= tmp;
            }
            for (i = 1; i <= minman; ++i) {
                j = ip[i];
                if (i != j) {
                    for (ll = 0; ll < mg; ++ll) {
                        tmp = w[(map1 + ll - 1) * nw + i];
                        w[(map1 + ll - 1) * nw + i] = w[(map1 + ll - 1) * nw + j];
                        w[(map1 + ll - 1) * nw + j] = tmp;
                    }
                }
            }
            if (0 < krank && krank < n) {
                for (ii = 1; ii <= krank; ++ii) {
                    i = krank + 1 - ii;
                    l = n1 + i;
                    lhh12(ctx, 2,i,krank+1,n,w+(i-1)*nw,1,ws+l-1,w+(map1-1)*nw,1,nw,mg);
                }
            }
            for (i = map1; i <= m; ++i) {
                if (0 < krank) {
                    for (j = 1; j <= krank; ++j) {
                        tmp = 0.0;
                        for (ll = 1; ll < j; ++ll)
                            tmp += w[(ll - 1) * nw + j] * w[(i - 1) * nw + ll];
                        w[(i - 1) * nw + j] =
                             (w[(i - 1) * nw + j] - tmp) / w[(j - 1) * nw + j];
                    }
                }
            }
        }
        
        /*  call lpdp to solve the reduced constrained problem.
            return codes are  0  if successful 
                             -2  if a solution was not found
                             -4  if insufficient memory
        */
   
        ierr = lpdp(ctx, w + (map1 -1) * nw,nw,mg,krank,n-krank,x,&xnorm);

        if (!ierr) {  /*  if successful */

            if (krank > 0) {
                for (ii = 1; ii <= krank; ++ii) {
                    i = krank + 1 - ii;
                    tmp = 0.0;
                    for (ll = 1; ll < ii; ++ll)
                        tmp += w[(i - 1) * nw + i + ll] * x[i + ll];
                    x[i] = (x[i] - tmp) / w[(i - 1) * nw + i];
                }
            }
            if (minman > 0) {

                /*  repermute variables to input order */

                for (i = minman; i >= 1; --i) {
                    j = ip[i];
                    if (i != j) {
                        tmp = x[i];
                        x[i] = x[j];
                        x[j] = tmp;
                    }
                }

                /*  add solution of unconstrained problem */

                for (i = 1; i <= n; ++i)  
                    x[i] += ws[i];

                /*  compute norm of residuals */

                *rnorm = sqrt(*rnorm * *rnorm + xnorm * xnorm);
            }
        }
    }
    else {      /* ma == m, no inequality constraints */

        for (ll = 1; ll <= n; ++ll)
            x[ll] = ws[ll];
    }
       
    if (ierr)
        goto LSIFin2;

    /*  calculation of the covariance matrix */

    if (cov && krank > 0) {

        krm1 = krank - 1;
        krp1 = krank + 1;

        for (ll = 0; ll < krank; ++ll)
            ws[n2 + ll] = w[ll * nw + ll + 1];

        for (j = 1; j <= krank; ++j) 
            w[(j - 1) * nw + j] = 1.0 / w[(j - 1) * nw + j];

        if (krank > 1) {
            for (i = 1; i <= krm1; ++i) {
                for (j = i + 1; j <= krank; ++j) {
                    tmp = 0.0;
                    for (ll = 0; ll < j - i; ++ll)
                        tmp += w[(i - 1) * nw + i + ll] *
                                        w[(i + ll - 1) * nw + j];
                    w[(i - 1) * nw + j] = -tmp * w[(j - 1) * nw + j];
                }
            }
        }
        for (i = 1; i <= krank; ++i) {
            for (j = i; j <= krank; ++j) {
                tmp = 0.0;
                for (ll = 0; ll < krank + 1 - j; ++ll)
                  tmp += w[(i - 1) * nw + j + ll] * w[(j - 1) * nw + j + ll];
                w[(i - 1) * nw + j] = tmp;
            }
        }
        if (krank < n) {
            for (j = 1; j <= krank; ++j) {
                for (ll = 1; ll <= j; ++ll)
                    w[(j - 1) * nw + ll] = w[(ll - 1) * nw + j];
            }
            for (i = krp1; i <= n; ++i) {
                for (ll = 1; ll <= i; ++ll)
                    w[(i - 1) * nw + ll] = 0.0;
            }
            n3 = n2 + krp1;
            for (i = 1; i <= krank; ++i) {
                l = n1 + i;
                k = n2 + i;
                rb = ws[l-1] * ws[k-1];

                if (rb < 0.0) {
 
                    rb = 1.0 / rb;
                    for (ll = 0; ll < n; ++ll)
                        ws[n3 + ll] = 0.0;

                    l = n1 + i;
                    k = n3 + i;
                    ws[k-1] = ws[l-1];

                    for (j = krp1; j <= n; ++j) {
                        k = n3 + j;
                        ws[k-1] = w[(i - 1) * nw + j];
                    }
                    for (j = 1; j <= n; ++j) {
                        l = n3 + i;
                        k = n3 + j;
                        tmp = 0.0;
                        for (ll = 0; ll < j - i; ++ll)
                            tmp += w[(j - 1) * nw + i + ll] * ws[l-1+ll];
                        tmp1 = 0.0;
                        for (ll = 0; ll < n - j + 1; ++ll)
                            tmp1 += w[(j + ll - 1) * nw + j] * ws[k-1+ll];
                        ws[j] = tmp + tmp1;
                        ws[j] *= rb;
                    }
                    l = n3 + i;
                    gam = 0.0;
                    for (ll = 0; ll < n - i + 1; ++ll)
                        gam += ws[l - 1 + ll] * ws[i + ll];
                    gam *= rb / 2.0;
                    for (ll = 0; ll < n - i + 1; ++ll)
                        ws[i + ll] += gam * ws[l-1 + ll];

                    for (j = i; j <= n; ++j) {
                        if (i > 1) {
                            k = n3 + j;
                            for (l = 1; l < i; ++l)  
                                w[(j - 1) * nw + l]  += ws[k-1] * ws[l];   
                        }
                        k = n3 + j;
                        for (l = i; l <= j; ++l) {
                            il = n3 + l;
                            w[(j - 1) * nw + l] += ws[j] * ws[il - 1] +
                                                   ws[l] * ws[k-1];
                        }
                    }
                }
            }
            for (i = 1; i <= n; ++i) {
                for (ll = 1; ll <= i; ++ll)
                    w[(ll - 1) * nw + i] = w[(i - 1) * nw + ll];
            }
        }

        /*  repermute rows and columns */

        for (i = minman; i >= 1; --i) {
            k = ip[i];
            if (i != k) {
                tmp = w[(i - 1) * nw + i];
                w[(i - 1) * nw + i] = w[(k - 1) * nw + k];
                w[(k - 1) * nw + k] = tmp;
                for (ll = 1; ll < i; ++ll) {
                    tmp = w[(ll - 1) * nw + i];
                    w[(ll - 1) * nw + i] = w[(ll - 1) * nw + k];
                    w[(ll - 1) * nw + k] = tmp;
                }
                for (ll = 1; ll < k - i; ++ll) {
                    tmp = w[(i - 1) * nw + i + ll];
                    w[(i - 1) * nw + i + ll] = w[(i + ll - 1) * nw + k];
                    w[(i + ll - 1) * nw + k] = tmp;
                }
                for (ll = k + 1; ll <= n; ++ll) {
                    tmp = w[(i - 1) * nw + ll];
                    w[(i - 1) * nw + ll] = w[(k - 1) * nw + ll];
                    w[(k - 1) * nw + ll] = tmp;
                }
            }
        }
        for (j = 1; j <= n; ++j) {
            for (ll = 1; ll <= j; ++ll)
                w[(ll - 1) * nw + j] *= fac;
            for (ll = 1; ll <= j; ++ll)
                w[(j - 1) * nw + ll] = w[(ll - 1) * nw + j];
        }
    }

LSIFin2:
    free((char *)ws);
    memrq(ctx, -lw - 1,sizeof(double));
LSIFin1:
    free((char *)ip);
    memrq(ctx, -n - 1,sizeof(int));
LSIFin:
    return(ierr);
}

/* ------------------------------------------------------------------------ */
/*  lpdp                                                                    */
/*                                                                          */
/*      This function is part of the lsei/lsi package to solve constrained  */
/*      least square problems, adapted from: R.J. Hanson, K.H. Haskell,     */
/*      Algorithm 587 from Collected Algorithms of ACM.                     */
/*                                                                          */
/*      The problem solved by this function is: determine a n1 vector w     */
/*      and a n2 vector z which minimize the Eucledian length of w s.t.     */
/*      Gw + Hz >= Y. G is a (m,n1)-matrix, H is a (m,n2)-matrix, Y is a    */
/*      (m)-vector.                                                         */
/*                                                                          */
/*      Call: int lpdp(a,na,m,n1,n2,x,wnorm)                                */
/*                                                                          */
/*      Input data are expected in a two-dimensional array a(i,j), stored   */
/*      in a one-dimensional array by a(i,j) = a[(i - 1) * na + j], as      */
/*      follows:                                                            */
/*                                                                          */
/*              A = (G,H,Y)   m rows, n1 + n2 + 1 columns.                  */
/*                                                                          */
/*      The solution, if found, is returned in the (n1 + n2)-vector x       */
/*      with x = (w,z). If successful, the norm of w is returned in wnorm.  */
/*                                                                          */
/*      Return codes:     = 0   Successful operation.                       */
/*                         -2   Inequality constraints are incompatible,    */
/*                         -4   Insufficient memory.                        */
        
int lpdp(TDAContext *ctx, double *a,int na,int m,int n1,int n2,double *x, double *wnorm)
{
    register int i,j,ll,l,np1;
    int nws,n,ierr;
    double sc,fac,rnorm,ynorm,tmp,*w,*wz;

    fac = 0.1;
    n = n1 + n2; 
    if (m <= 0) {
        for (i = 1; i <= n; ++i)
            x[i] = 0.0;
        *wnorm = 0.0;
        return(0);
    }
    ierr = 0;
    np1 = n + 1;
    nws = m + 1;

    /*  allocate  temporary storage */

    if (!(w = (double *) calloc((size_t)(nws) * (size_t)(np1) + 1,sizeof(double)))) {
        ierr = -4;
        goto LPDPFin;
    }
    memrq(ctx, nws * np1 + 1,sizeof(double));

    if (!(wz = (double *) calloc((size_t)(m + 1),sizeof(double)))) {
        ierr = -4;
        goto LPDPFin1;
    }
    memrq(ctx, m + 1,sizeof(double));
   
    for (i = 1; i <= m; ++i) {
        tmp = 0.0;
        for (ll = 1; ll <= n; ++ll)
            tmp += a[(i - 1) * na + ll] * a[(i - 1) * na + ll];
        sc = sqrt(tmp);

        if (sc != 0.0) {
            sc = 1.0 / sc;
            for (ll = 1; ll <= np1; ++ll)
                a[(i - 1) * na + ll] *= sc;
        }
    }
    tmp = 0.0;
    for (ll = 1; ll <= m; ++ll)
        tmp += a[(ll - 1) * na + np1] * a[(ll - 1) * na + np1];
    ynorm = sqrt(tmp);
    if (ynorm != 0.0) {
        sc = 1.0 / ynorm;
        for (ll = 1; ll <= m; ++ll)
            a[(ll - 1) * na + np1] *= sc;
    }
    j = n1 + 1;

    while (j <= n) {
        tmp = 0.0;
        for (ll = 1; ll <= m; ++ll)
            tmp += a[(ll - 1) * na + j] * a[(ll - 1) * na + j];
        sc = sqrt(tmp);

        if (sc != 0.0)
            sc = 1.0 / sc;
        for (ll = 1; ll <= m; ++ll)
            a[(ll - 1) * na + j] *= sc;
        x[j] = sc;
        j++;
    }
    if (n1 > 0) {
        for (i = 1; i <= m; ++i) {
            for (ll = 1; ll <= n2; ++ll)
                w[(ll - 1) * nws + i] = a[(i - 1) * na + n1 + ll];  

            for (ll = 1; ll <= n1; ++ll)
                w[(n2 + ll - 1) * nws + i] = a[(i - 1) * na + ll];  

            w[(n1 + n2) * nws + i] = a[(i - 1) * na + np1];  
        }
        for (ll = 1; ll <= n; ++ll)  
            w[(ll - 1) * nws + nws] = 0.0;
        w[(n1 + n2) * nws + nws] = 1.0;

        ierr = wnnls(ctx, w,m+1,n2,np1 - n2,m,0,wz,&rnorm);
        if (ierr)
            goto LPDPFin2;

        tmp = 0.0;
        for (ll = 1; ll <= m; ++ll)
            tmp += a[(ll - 1) * na + np1] * wz[ll]; 
        sc = 1.0 - tmp;

        if (1.0 + fac * fabs(sc) != 1.0 && rnorm > 0.0) {
            sc = 1.0 / sc;
            for (j = 1; j <= n1; ++j) {
                tmp = 0.0;
                for (ll = 1; ll <= m; ++ll)
                    tmp += a[(ll - 1) * na + j] * wz[ll];
                x[j] = sc * tmp;
            }
            for (i = 1; i <= m; ++i) {
                tmp = 0.0;
                for (ll = 1; ll <= n1; ++ll)
                    tmp += a[(i - 1) * na + ll] * x[ll];
                a[(i - 1) * na + np1] -= tmp;
            }
        }
        else {
            ierr = -2;
            goto LPDPFin2;
        }
    }
    if (n2 > 0) {
        for (i = 1; i <= m; ++i) {
            for (ll = 1; ll <= n2; ++ll)  
                w[(ll - 1) * nws + i] = a[(i - 1) * na + n1 + ll];  
            w[n2 * nws + i] = a[(i - 1) * na + np1];  
        }
        for (ll = 1; ll <= n2; ++ll)  
            w[(ll - 1) * nws + nws] = 0.0;
        w[n2 * nws + nws] = 1.0;

        ierr = wnnls(ctx, w,m+1,0,n2 +1,m,0,wz,&rnorm);
        if (ierr)
            goto LPDPFin2;

        tmp = 0.0;
        for (ll = 1; ll <= m; ++ll)
            tmp += a[(ll - 1) * na + np1] * wz[ll];
        sc = 1.0 - tmp;

        if ((1.0 + fac * fabs(sc) != 1.0 && rnorm > 0.0)) {
            sc = 1.0 / sc;
            for (j = 1; j <= n2; ++j) {
                l = n1 + j;
                tmp = 0.0;
                for (ll = 1; ll <= m; ++ll)  
                    tmp += a[(ll - 1) * na + l] * wz[ll];
                x[l] *= sc * tmp;
            }
        }
        else {
            ierr = -2;
            goto LPDPFin2;
        }
    }
    for (ll = 1; ll <= n; ++ll) 
        x[ll] *= ynorm;
    tmp = 0.0;
    for (ll = 1; ll <= n1; ++ll) 
        tmp += x[ll] * x[ll];
    *wnorm = sqrt(tmp);

LPDPFin2:
    free((char *)wz);
    memrq(ctx, -m - 1,sizeof(double));
LPDPFin1:
    free((char *)w);
    memrq(ctx, -nws * np1 - 1,sizeof(double));
LPDPFin:
    return(ierr);
}   

/* ------------------------------------------------------------------------ */
/*  wnnls                                                                   */
/*      This function is part of the lsei/lsi package to solve constrained  */
/*      least square problems, adapted from: R.J. Hanson, K.H. Haskell,     */
/*      Algorithm 587 from Collected Algorithms of ACM.                     */
/*                                                                          */
/*      Return:     0  if successful                                        */
/*                 -4  if insufficient memory                               */
/*                 -5  if exceeded max number of iterations in wnlsm        */
        
int wnnls(TDAContext *ctx, double *w,int nw,int me,int ma,int n,int l,double *x,double *rnorm)
{
    int ierr,k,*ipivot,*itype;
    double *wd,*h,*scale,*z,*temp,*d;

    if (ma + me <= 0 || n <= 0)  
        return(0);

    /*  allocate  temporary storage */

    k = me + ma;
    if (k < n)
        k = n;

    ierr = -4;
    if (!(ipivot = (int *) calloc((size_t)(n + 1),sizeof(int))))   
        goto WNNLSFin;
    memrq(ctx, n + 1,sizeof(int));
   
    if (!(itype = (int *) calloc((size_t)(k + 1),sizeof(int))))  
        goto WNNLSFin1;
    memrq(ctx, k + 1,sizeof(int));
   
    if (!(wd = (double *) calloc((size_t)(n + 1),sizeof(double))))  
        goto WNNLSFin2;
    memrq(ctx, n + 1,sizeof(double));
   
    if (!(h = (double *) calloc((size_t)(n + 1),sizeof(double))))  
        goto WNNLSFin3;
    memrq(ctx, n + 1,sizeof(double));
   
    if (!(scale = (double *) calloc((size_t)(k + 1),sizeof(double))))  
        goto WNNLSFin4;
    memrq(ctx, k + 1,sizeof(double));

    if (!(z = (double *) calloc((size_t)(n + 1),sizeof(double))))  
        goto WNNLSFin5;
    memrq(ctx, n + 1,sizeof(double));
  
    if (!(temp = (double *) calloc((size_t)(n + 1),sizeof(double))))  
        goto WNNLSFin6;
    memrq(ctx, n + 1,sizeof(double));

    if (!(d = (double *) calloc((size_t)(n + 1),sizeof(double))))  
        goto WNNLSFin7;
    memrq(ctx, n + 1,sizeof(double));

    ierr = wnlsm(ctx, w,nw,me,ma,n,l,x,rnorm,ipivot,itype,wd,h,scale,z,temp,d);

    free((char *)d);
    memrq(ctx, -n - 1,sizeof(double));
WNNLSFin7:
    free((char *)temp);
    memrq(ctx, -n - 1,sizeof(double));
WNNLSFin6:
    free((char *)z);
    memrq(ctx, -n - 1,sizeof(double));
WNNLSFin5:
    free((char *)scale);
    memrq(ctx, -k - 1,sizeof(double));
WNNLSFin4:
    free((char *)h);
    memrq(ctx, -n - 1,sizeof(double));
WNNLSFin3:
    free((char *)wd);
    memrq(ctx, -n - 1,sizeof(double));
WNNLSFin2:
    free((char *)itype);
    memrq(ctx, -k - 1,sizeof(int));
WNNLSFin1:
    free((char *)ipivot);
    memrq(ctx, -n - 1,sizeof(int));
WNNLSFin:
    return(ierr);
}

/* ------------------------------------------------------------------------ */
/*  wnlsm                                                                   */
/*      This function is part of the lsei/lsi package to solve constrained  */
/*      least square problems, adapted from: R.J. Hanson, K.H. Haskell,     */
/*      Algorithm 587 from Collected Algorithms of ACM.                     */
/*                                                                          */
/*      Return:     0  if successful                                        */
/*                 -5  if exceeded max number of iterations                 */
           
int wnlsm(TDAContext *ctx, double *w,int nw,int mme,int ma,int n,int l, double *x,double *rnorm,int *ipivot,int *itype, double *wd,double *h,double *scale,double *z,double *temp,double *d)
{
    register int i,j,l1,ll,jj;
    int m,me,done,hitcon,iter,itmax,np1,lp1,itemp,mep1,nsoln,krank,krp1;
    int imaxx,nm1,niv,niv1,feasbl,jcon,nsp1,iwmax,jm1,jp,isol,pos;
    int ierr,idope[9];
    double tau,fac,blowup,eanorm,bnorm,alamda,alsq,t,zz,alpha,sm;
    double amax,wmax,z2,tmp,tmp1,dope[5],sparam[6];

    pos = hitcon = done = ierr = 0;
    m = ma + mme;
    me = mme;
    mep1 = me + 1;
    fac = 1.E-4;

    tau = sqrt(ctx->EPSI);
    blowup = tau;
    for (ll = 1; ll <= n; ++ll)
        d[ll] = 1.0;

    for (j = 1; j <= n; ++j) {
        for (ll = 1; ll <= m; ++ll)
            w[(ll - 1) * nw + j] *= d[j];
    }
    iter = 0;
    itmax = 3 * (n - l);
    lp1 = l + 1;
    nsoln = l;
    nsp1 = nsoln + 1;
    np1 = n + 1;
    nm1 = n - 1;
    l1 = m;
    if (m > l)
        l1 = l;

    for (j = 1; j <= n; ++j) {
        tmp = 0.0;
        for (ll = 1; ll <= m; ++ll)
            tmp += fabs(w[(ll - 1) * nw + j]);
        wd[j] = tmp;
    }
    imaxx = 1;
    tmp = fabs(wd[1]);
    for (ll = 2; ll <= n; ++ll) {
        tmp1 = fabs(wd[ll]);
        if (tmp1 > tmp) {
            imaxx = ll;
            tmp = tmp1;
        }
    }
    eanorm = wd[imaxx];
    bnorm = 0.0;
    for (ll = 1; ll <= m; ++ll)
        bnorm += fabs(w[(ll - 1) * nw + np1]);

    alamda = eanorm / (ctx->EPSI * fac);
    alsq = alamda * alamda;

    for (i = 1; i <= m; ++i) {
        if (i <= me) {
            t = alsq;
            itemp = 0;
        }
        else {
            t = 1.0;
            itemp = 1;
        }
        scale[i] = t;
        itype[i] = itemp;
    }
    for (ll = 1; ll <= n; ++ll) {
        x[ll] = 0.0;
        ipivot[ll] = ll;
    }
    for (ll = 1; ll <= l; ++ll)
        wd[ll] = 0.0;

    idope[1] = me;
    idope[2] = mep1;
    idope[3] = 0;
    idope[4] = 1;
    idope[5] = nsoln;
    idope[6] = 0;
    idope[7] = 1;
    idope[8] = l1;
 
    dope[1] = alsq;
    dope[2] = eanorm;
    dope[3] = fac;
    dope[4] = tau;

    done = wnlit(ctx, w,nw,m,n,l,ipivot,itype,h,scale,rnorm,idope,dope,done);

    me = idope[1];
    mep1 = idope[2];
    krank = idope[3];
    krp1 = idope[4];
    nsoln = idope[5];
    niv = idope[6];
    niv1 = idope[7];
    l1 = idope[8];

L20:
    if (done)   
        isol = 1;
    else
        isol = lp1;

    if (nsoln >= isol) {
        for (ll = 1; ll <= niv; ++ll)
            temp[ll] = w[(ll - 1) * nw + np1];

        for (jj = isol; jj <= nsoln; ++jj) {
            j = nsoln - jj + isol;
            if (j > krank) {
                i = niv - jj + isol;
            }
            else
                i = j;
            if (j > krank && j <= l) {
                z[j] = 0.0;
            }
            else {
                z[j] = temp[i] / w[(i - 1) * nw + j];
                for (ll = 1; ll <= i - 1; ++ll)
                    temp[ll] -= z[j] * w[(ll - 1) * nw + j];
            }
        }
    }
    if (done) {
        for (ll = 1; ll <= nsoln; ++ll)  
            x[ll] = z[ll];

        if (0 < krank && krank < l) {
            for (i = 1; i <= krank; ++i) {
                lhh12(ctx, 2,i,krp1,l,w+(i-1)*nw,1,h+i,x,1,1,1);
            }
        }
        if (nsoln < n) {
            for (ll = 0; ll < n - nsoln; ++ll)   
                x[nsp1 + ll] = 0.0;
        }
        for (i = 1; i <= n; ++i) {
            j = i;
            while (ipivot[j] != i) {
                j++;
                if (j > n) {
                    printfe(ctx, "Error in wnnls (pivot search)\n");
                    gerr_exit(ctx, 990);
                }
            }
            ipivot[j] = ipivot[i];
            ipivot[i] = j;
            tmp = x[j];
            x[j] = x[i];
            x[i] = tmp;
        }
        for (j = 1; j <= n; ++j)  
            x[j] *= d[j];

        if (niv < m) {
            for (i = niv1; i <= m; ++i) {
                t = w[(i - 1) * nw + np1];
                if (i <= me)
                    t /= alamda;
                t = (scale[i] * t) * t;
                *rnorm += t;
            }
        }
        *rnorm = sqrt(*rnorm);
        return(ierr);
    }
    if (++iter > itmax) {
        ierr = -5;
        done = 1;
    }
    alpha = 2.0;
    hitcon = 0;
    if (l < nsoln) {
        for (j = lp1; j <= nsoln; ++j) {
            zz = z[j];
            if (zz <= 0.0) {
                t = x[j] / (x[j] - zz);
                if (t < alpha) {
                    alpha = t;
                    jcon = j;
                }
                hitcon = 1;
            }
        }
    }
    if (hitcon) {
        if (lp1 <= nsoln) {
            for (j = lp1; j <= nsoln; ++j) {
                x[j] = x[j] + alpha * (z[j] - x[j]);
            }
        }
        feasbl = 0;
        while (!feasbl) {
            for (i = 1; i <= m; ++i) {
                t = w[(i - 1) * nw + jcon];
                for (ll = 0; ll < n - jcon; ++ll)
                    w[(i - 1) * nw + jcon + ll] = w[(i - 1) * nw + jcon + 1+ll];
                w[(i - 1) * nw + n] = t;
            }
            itemp = ipivot[jcon];
            if (jcon < n) {
                for (i = jcon; i <= nm1; ++i) 
                    ipivot[i] = ipivot[i + 1];
            }
            ipivot[n] = itemp;
            for (ll = 0; ll < n - jcon; ++ll)
                x[jcon + ll] = x[jcon + 1 + ll];
    
            x[n] = 0.0;
            nsp1 = nsoln;
            nsoln--;
            niv1 = niv;
            niv--;
     
            j = jcon;
            i = krank + jcon - l;
    
            while (j <= nsoln) {
                if (itype[i] == 0 && itype[i+1] == 0) {
                    if (w[i * nw + j] != 0.0) {
                        srotmg(ctx, scale+i,scale+i+1,w+(i-1)*nw+j,w+i*nw+j,sparam);
                        w[i * nw + j] = 0.0;
                        srotm(ctx, np1-j,w+(i-1)*nw+j,w+i*nw+j,sparam);
                    }
                }
                else if (itype[i] == 1 && itype[i+1] == 1) {
                    if (w[i * nw + j] != 0.0) {
                        srotmg(ctx, scale+i,scale+i+1,w+(i-1)*nw+j,w+i*nw+j,sparam);
                        w[i * nw + j] = 0.0;
                        srotm(ctx, np1-j,w+(i-1)*nw+j,w+i*nw+j,sparam);
                    }
                }
                else if (itype[i] == 1 && itype[i+1] == 0) {
                    for (ll = 1; ll <= np1; ++ll) {
                        tmp = w[(i - 1) * nw + ll];
                        w[(i - 1) * nw + ll] = w[(i + 1 - 1) * nw + ll];
                        w[(i + 1 - 1) * nw + ll] = tmp;
                    }
                    tmp = scale[i];
                    scale[i] = scale[i + 1];
                    scale[i + 1] = tmp;
    
                    itemp = itype[i+1];
                    itype[i+1] = itype[i];
                    itype[i] = itemp;
      
                    if (w[i * nw + j] != 0.0) {
                        srotmg(ctx, scale+i,scale+i+1,w+(i-1)*nw+j,w+i*nw+j,sparam);
                        w[i * nw + j] = 0.0;
                        srotm(ctx, np1-j,w+(i-1)*nw+j,w+i*nw+j,sparam);
                    }
                }
                else if (itype[i] == 0 && itype[i+1] == 1) {
                    t = scale[i] * w[(i - 1) * nw + j] * w[(i - 1) * nw + j] / alsq;
                    if (t > tau * tau * eanorm * eanorm) {
                        if (w[i * nw + j] != 0.0) {
                            srotmg(ctx, scale+i,scale+i+1,w+(i-1)*nw+j,w+i*nw+j,sparam);
                            w[i * nw + j] = 0.0;
                            srotm(ctx, np1-j,w+(i-1)*nw+j,w+i*nw+j,sparam);
                        }
                    }
                    else {
                        for (ll = 1; ll <= np1; ++ll) {
                            tmp = w[(i - 1) * nw + ll];
                            w[(i - 1) * nw + ll] = w[(i + 1 - 1) * nw + ll];
                            w[(i + 1 - 1) * nw + ll] = tmp;
                        }
                        tmp = scale[i];
                        scale[i] = scale[i + 1];
                        scale[i + 1] = tmp;
                        itemp = itype[i+1];
                        itype[i+1] = itype[i];
                        itype[i] = itemp;
                        w[i * nw + j] = 0.0;
                    }
                }
                i++;
                j++;
            }
            if (lp1 <= nsoln) {
                for (jcon = lp1; jcon <= nsoln; ++jcon) {
                    if (x[jcon] <= 0.0)
                        goto L600;
                }
            }
            feasbl = 1;
L600:       ;
        }
    }
    else {
        for (ll = 1; ll <= nsoln; ++ll)
            x[ll] = z[ll];
    
        if (nsoln < n) {
            for (ll = 0; ll < n - nsoln; ++ll)
                x[nsp1 + ll] = 0.0;
        }
        i = niv1;
    
        while (i <= me) {
     
            if (itype[i] == 0) {
                i++;
            }
            else {
                for (ll = 1; ll <= np1; ++ll) {
                    tmp = w[(i - 1) * nw + ll];
                    w[(i - 1) * nw + ll] = w[(me - 1) * nw + ll];
                    w[(me - 1) * nw + ll] = tmp;
                }
                tmp = scale[i];
                scale[i] = scale[me];
                scale[me] = tmp;
    
                itemp = itype[i];
                itype[i] = itype[me];
                itype[me] = itemp;
                mep1 = me;
                me--;
            }
        }
        if (nsp1 <= n) {
            for (j = nsp1; j <= n; ++j) {
                sm = 0.0;
                if (nsoln < m) {
                    for (i = nsp1; i <= m; ++i)  
                        sm += scale[i] * w[(i - 1) * nw + j] * w[(i - 1) * nw + np1];
                }
                wd[j] = sm;
            }
        }
        goto L750;
L740:
        if (!pos && !done) {
L750:
            wmax = 0.0;
            iwmax = nsp1;
            if (nsp1 <= n) {
                for (j = nsp1; j <= n; ++j) {
                    if (wd[j] > wmax) {
                        wmax = wd[j];
                        iwmax = j;
                    }
                }
            }
            if (wmax <= 0.0) {
                done = 1;
            }
            else {
                wd[iwmax] = 0.0;
                nsoln = nsp1;
                nsp1 = nsoln + 1;
                niv = niv1;
                niv1 = niv + 1;
                if (nsoln != iwmax) {
                    for (ll = 1; ll <= m; ++ll) {
                        tmp = w[(ll - 1) * nw + nsoln];
                        w[(ll - 1) * nw + nsoln] = w[(ll - 1) * nw + iwmax];
                        w[(ll - 1) * nw + iwmax] = tmp;
                    }
                    wd[iwmax] = wd[nsoln];
                    wd[nsoln] = 0.0;
                    itemp = ipivot[nsoln];
                    ipivot[nsoln] = ipivot[iwmax];
                    ipivot[iwmax] = itemp;
                }
                j = m;
                while (j > niv) {
                    jm1 = j - 1;
                    jp = jm1;
         
                    if (j == mep1) {
                        imaxx = me;
                        amax = scale[me] * w[(me - 1) * nw + nsoln] *
                                                          w[(me - 1) * nw + nsoln];
                        while (jp >= niv) {
                            t = scale[jp] * w[(jp - 1) * nw + nsoln] *
                                                          w[(jp - 1) * nw + nsoln];
                            if (t > amax) {
                                imaxx = jp;
                                amax = t;
                            }
                            jp--;
                        }
                        jp = imaxx;
                    }
                    if (w[(j - 1) * nw + nsoln] != 0.0) {
                        srotmg(ctx, scale+jp,scale+j,w+(jp-1)*nw+nsoln,
                                                         w+(j-1)*nw+nsoln,sparam);
                        w[(j - 1) * nw + nsoln] = 0.0;
                        srotm(ctx, np1-nsoln,w+(jp-1)*nw+nsp1-1,w+(j-1)*nw+nsp1-1,sparam);
                    }
                    j = jm1;
                }
                if (w[(niv - 1) * nw + nsoln] != 0.0) {
                    isol = niv;
                    z2 = w[(isol - 1) * nw + np1] / w[(isol - 1) * nw + nsoln];
                    z[nsoln] = z2;
                    if (z2 > 0.0)
                        pos = 1;
                    else
                        pos = 0;
    
                    if (z2 * eanorm >= bnorm && pos) {
                        if (blowup * z2 * eanorm >= bnorm)
                            pos = 0;
                        else
                            pos = 1;
                    }
                }
                else if (niv <= me && w[(mep1 - 1) * nw + nsoln] != 0.0) {
                    isol = mep1;
                    z2 = w[(isol - 1) * nw + np1] / w[(isol - 1) * nw + nsoln];
                    z[nsoln] = z2;
                    if (z2 > 0.0)
                        pos = 1;
                    else
                        pos = 0;
                    if (z2 * eanorm >= bnorm && pos) {
                        if (blowup * z2 * eanorm >= bnorm)
                            pos = 0;
                        else
                            pos = 1;
                    }
                    if (pos) {
                        for (ll = 1; ll <= np1; ++ll) {
                            tmp = w[(mep1 - 1) * nw + ll];
                            w[(mep1 - 1) * nw + ll] = w[(niv - 1) * nw + ll];
                            w[(niv - 1) * nw + ll] = tmp;
                        }
                        tmp = scale[mep1];
                        scale[mep1] = scale[niv];
                        scale[niv] = tmp;
        
                        itemp = itype[mep1];
                        itype[mep1] = itype[niv];
                        itype[niv] = itemp;
                        me = mep1;
                        mep1 = me + 1;
                    }
                }
                else
                    pos = 0;

                if (!pos) {
                    nsp1 = nsoln;
                    nsoln--;
                    niv1 = niv;
                    niv--;
                }
            }
            goto L740;
        }
    }
    goto L20;
}

/* ------------------------------------------------------------------------ */
/*  wnlit                                                                   */
/*      This function is part of the lsei/lsi package to solve constrained  */
/*      least square problems, adapted from: R.J. Hanson, K.H. Haskell,     */
/*      Algorithm 587 from Collected Algorithms of ACM.                     */

int wnlit(TDAContext *ctx, double *w,int nw,int m,int n,int l,int *ipivot,int *itype, double *h,double *scale,double *rnorm,int *idope,double *dope,int done)
{
    register int i,j,ll,ip1,ir,ic,jj,k;
    int lb,me,mep1,krank,krp1,nsoln,niv,niv1,l1,recalc,jp,jm1; 
    int max,lend,mend,indep,kk,np1,lp1,irp1,i1,j1,itemp;
    double alsq,eanorm,tau,factor,t,hbar,rn,sn,tmp,tmp1,sparam[6];
    double tenm3 = 1.e-3;

    hbar = 0.0;
    me = idope[1];
    mep1 = idope[2];
    krank = idope[3];
    krp1 = idope[4];
    nsoln = idope[5];
    niv = idope[6];
    niv1 = idope[7];
    l1 = idope[8];
    alsq = dope[1];
    eanorm = dope[2];

    tau = dope[4];
    np1 = n + 1;
    lb = m - 1;
    if (l < lb)
        lb = l;

    recalc = 1;
    *rnorm = 0.0; 
    krank = 0;
    factor = 1.0;
    i = 1;
    ip1 = 2;
    lend = l;

    while (i <= lb) {
        if (i <= me) {
            ir = i;
            mend = m;
            if (ir != 1 && !recalc) {
 
                for (j = i; j <= lend; ++j)  
                    h[j] -= scale[ir-1] * w[(ir - 2) * nw + j] *
                                                         w[(ir - 2) * nw + j];
                max = 1;
                tmp = fabs(h[i]);
                for (ll = 2; ll <= lend - i + 1; ++ll ) {
                    tmp1 = fabs(h[i - 1 + ll]);
                    if (tmp < tmp1) {
                        max = ll;
                        tmp = tmp1;
                    }
                }
                max = max + i - 1;
                tmp = hbar + tenm3 * h[max];
                if (tmp == hbar)
                    recalc = 1;
                else
                    recalc = 0;
            }
            if (recalc) {
                for (j = i; j <= lend; ++j) {
                    h[j] = 0.0;
                    for (k = ir; k <= mend; ++k)  
                        h[j] += scale[k] * w[(k - 1) * nw + j] *
                                                       w[(k - 1) * nw + j];
                }
                max = 1;
                tmp = fabs(h[i]);
                for (ll = 2; ll <= lend - i + 1; ++ll ) {
                    tmp1 = fabs(h[i - 1 + ll]);
                    if (tmp < tmp1) {
                        max = ll;
                        tmp = tmp1;
                    }
                }
                max = max + i - 1;
                hbar = h[max];
            }
            if (max != i) {
                itemp = ipivot[i];
                ipivot[i] = ipivot[max];
                ipivot[max] = itemp;
                for (ll = 1; ll <= m; ++ll) {
                    tmp = w[(ll - 1) * nw + max];
                    w[(ll - 1) * nw + max] = w[(ll - 1) * nw + i];
                    w[(ll - 1) * nw + i] = tmp;
                }
                t = h[max];
                h[max] = h[i];
                h[i] = t;
            }
L30:
            ic = i;
            sn = rn = 0.0;
            for (j = 1; j <= mend; ++j) {
                t = scale[j];
                if (j <= me)
                    t /= factor;
                t *= w[(j - 1) * nw + ic] * w[(j - 1) * nw + ic];
                if (j < ir)  
                    sn += t;
                else 
                    rn += t;
            }
            tmp = tau * tau * sn;
            if (rn > tmp)
                indep = 1;
            else
                indep = 0;

            if (indep) {
                j = m;
                for (jj = ip1; jj <= m; ++jj) {
                    jm1 = j - 1;
                    jp = jm1;
                    if (jj == m) {
                        if (i <  mep1) {
                            j = mep1;
                            jp = i;
                            t = scale[jp] * w[(jp-1) * nw + i] *
                                            w[(jp-1) * nw + i] * tau * tau;
                            tmp = scale[j] * w[(j - 1) * nw + i] *
                                             w[(j - 1) * nw + i];
                            if (t <= tmp)
                                goto L130;
                        }
                    }
                    else if (j == mep1) {
                        j = jm1;
                        jm1 = j - 1;
                        jp = jm1;
                    }
                    if (w[(j - 1) * nw + i] != 0.0) {
                        srotmg(ctx, scale+jp,scale+j,w+(jp-1)*nw+i,
                                                          w+(j-1)*nw+i,sparam);
                        w[(j - 1) * nw + i] = 0.0;
                        srotm(ctx, np1-i,w+(jp-1)*nw+ip1-1,w+(j-1)*nw+ip1-1,sparam);
                    }
                    j = jm1;
                }
                goto L140;
            } 
            if (lend > i) {
                max = lend;
                if (max != i) {
                    itemp = ipivot[i];
                    ipivot[i] = ipivot[max];
                    ipivot[max] = itemp;
                    for (ll = 1; ll <= m; ++ll) {
                        tmp = w[(ll - 1) * nw + max];
                        w[(ll - 1) * nw + max] = w[(ll - 1) * nw + i];
                        w[(ll - 1) * nw + i] = tmp;
                    }
                    t = h[max];
                    h[max] = h[i];
                    h[i] = t;
                }
                lend--;
                max = 1;
                tmp = fabs(h[i]);
                for (ll = 2; ll <= lend - i + 1; ++ll  ) {
                    tmp1 = fabs(h[i - 1 + ll]);
                    if (tmp < tmp1) {
                        max = ll;
                        tmp = tmp1;
                    }
                }
                max = max + i - 1;
                hbar = h[max];
                goto L30;
            }
        }
L130:
        krank = i - 1;
        goto L160;
L140:
        i = ip1;
        ip1++;
    }
    krank = l1;
L160:
    krp1 = krank + 1;
    if (krank < me) {
        factor = alsq;
        for (i = krp1; i <= me; ++i) {
            if (l > 0)  
                tmp = 0.0;
            else
                tmp = w[(i - 1) * nw + 1];
            for (ll = 2; ll <= l; ++ll)
                w[(i - 1) * nw + ll] = tmp;
        }
        ir = krp1;
        if (l < n) {
            lp1 = l + 1;
            recalc = 1;
            lb = l + me - krank;
            if (lb > n)
                lb = n;
            i = lp1;
            ip1 = i + 1;
L180:
            if (i > lb)
                goto L280;

            ir = krank + i - l;
            lend = n;
            mend = me;
            if (ir != 1 && !recalc) {
 
                for (j = i; j <= lend; ++j)  
                    h[j] -= scale[ir-1] * w[(ir - 2) * nw + j] *
                                                        w[(ir - 2) * nw + j];
                max = 1;
                tmp = fabs(h[i]);
                for (ll = 2; ll <= lend - i + 1; ++ll ) {
                    tmp1 = fabs(h[i - 1 + ll]);
                    if (tmp < tmp1) {
                        max = ll;
                        tmp = tmp1;
                    }
                }
                max = max + i - 1;
                tmp = hbar + tenm3 * h[max];
                if (tmp == hbar)
                    recalc = 1;
                else
                    recalc = 0;
            }
            if (recalc) {
                for (j = i; j <= lend; ++j) {
                    h[j] = 0.0;
                    for (k = ir; k <= mend; ++k)  
                        h[j] += scale[k] * w[(k - 1) * nw + j] *
                                                          w[(k - 1) * nw + j];
                }
                max = 1;
                tmp = fabs(h[i]);
                for (ll = 2; ll <= lend - i + 1; ++ll ) {
                    tmp1 = fabs(h[i - 1 + ll]);
                    if (tmp < tmp1) {
                        max = ll;
                        tmp = tmp1;
                    }
                }
                max = max + i - 1;
                hbar = h[max];
            }
            if (max != i) {
                itemp = ipivot[i];
                ipivot[i] = ipivot[max];
                ipivot[max] = itemp;
                for (ll = 1; ll <= m; ++ll) {
                    tmp = w[(ll - 1) * nw + max];
                    w[(ll - 1) * nw + max] = w[(ll - 1) * nw + i];
                    w[(ll - 1) * nw + i] = tmp;
                }
                t = h[max];
                h[max] = h[i];
                h[i] = t;
            }
            j = me;
            while (j > ir) {
                jm1 = j - 1;
                if (w[(j - 1) * nw + i] != 0.0) {
                    srotmg(ctx, scale+jm1,scale+j,w+(jm1-1)*nw+i,w+(j-1)*nw+i,sparam);
                    w[(j - 1) * nw + i] = 0.0;
                    srotm(ctx, np1-i,w+(jm1-1)*nw+ip1-1,w+(j-1)*nw+ip1-1,sparam);
                }
                j = jm1;
            }
            ic = i;
            sn = rn = 0.0;
            for (j = 1; j <= mend; ++j) {
                t = scale[j];
                if (j <= me)
                    t /= factor;
                t *= w[(j - 1) * nw + ic] * w[(j - 1) * nw + ic];
                if (j < ir)  
                    sn += t;
                else 
                    rn += t;
            }
            tmp = tau * tau * sn;
            if (rn > tmp)
                indep = 1;
            else
                indep = 0;

            if (indep)
                goto L270;
        }
        jj = ir;
        while (ir <= me) {
            for (ll = 1; ll <= n; ++ll)  
                w[(ir - 1) * nw + ll] = 0.0;
            *rnorm += (scale[ir] * w[(ir - 1) * nw + np1] / alsq) *
                                                     w[(ir - 1) * nw + np1];
            w[(ir - 1) * nw + np1] = 0.0;
            scale[ir] = 1.0;
            itype[ir] = 1;
            ir++;
        }
        me = jj - 1;
        mep1 = me + 1;
        goto L300;
L270:
        i = ip1;
        ip1++;
        goto L180;
L280:       ;
    }
L300:
    if (krank < l1) {
 
        recalc = 1;
        factor = alsq;
        kk = krp1;
        i = kk;
        ip1 = i + 1;

        while (i <= l1) {
            ir = mep1;
            lend = l;
            mend = m;
            if (ir != 1 && !recalc) {
 
                for (j = i; j <= lend; ++j)  
                    h[j] -= scale[ir-1] * w[(ir - 2) * nw + j] *
                                                         w[(ir - 2) * nw + j];
                max = 1;
                tmp = fabs(h[i]);
                for (ll = 2; ll <= lend - i + 1; ++ll ) {
                    tmp1 = fabs(h[i - 1 + ll]);
                    if (tmp < tmp1) {
                        max = ll;
                        tmp = tmp1;
                    }
                }
                max = max + i - 1;
                tmp = hbar + tenm3 * h[max];
                if (tmp == hbar)
                    recalc = 1;
                else
                    recalc = 0;
            }
            if (recalc) {
                for (j = i; j <= lend; ++j) {
                    h[j] = 0.0;
                    for (k = ir; k <= mend; ++k)  
                        h[j] += scale[k] * w[(k - 1) * nw + j] *
                                                         w[(k - 1) * nw + j];
                }
                max = 1;
                tmp = fabs(h[i]);
                for (ll = 2; ll <= lend - i + 1; ++ll ) {
                    tmp1 = fabs(h[i - 1 + ll]);
                    if (tmp < tmp1) {
                        max = ll;
                        tmp = tmp1;
                    }
                }
                max = max + i - 1;
                hbar = h[max];
            }
            if (max != i) {
                itemp = ipivot[i];
                ipivot[i] = ipivot[max];
                ipivot[max] = itemp;
                for (ll = 1; ll <= m; ++ll) {
                    tmp = w[(ll - 1) * nw + max];
                    w[(ll - 1) * nw + max] = w[(ll - 1) * nw + i];
                    w[(ll - 1) * nw + i] = tmp;
                }
                t = h[max];
                h[max] = h[i];
                h[i] = t;
            }
            irp1 = ir + 1;
            if (irp1 <= m) {
                j = m;
                for (jj = irp1; jj <= m; ++jj) {
                    jm1 = j - 1;
                    if (w[(j - 1) * nw + i] != 0.0) {
                        srotmg(ctx, scale+jm1,scale+j,w+(jm1-1)*nw+i,w+(j-1)*nw+i,sparam);
                        w[(j - 1) * nw + i] = 0.0;
                        srotm(ctx, np1-i,w+(jm1-1)*nw+ip1-1,w+(j-1)*nw+ip1-1,sparam);
                    }
                    j = jm1;
                }
            }
            t = scale[ir] * w[(ir - 1) * nw + i] *
                            w[(ir - 1) * nw + i];
            tmp = tau * tau * eanorm * eanorm;
            if (t > tmp)
                indep = 1;
            else
                indep = 0;

            if (indep) {
                rn = 0.0;
                for (i1 = ir; i1 <= m; ++i1) {
                    for (j1 = ip1; j1 <= n; ++j1) {
                        tmp = scale[i1] * w[(i1 - 1) * nw + j1] *
                                                        w[(i1 - 1) * nw + j1];
                        if (rn < tmp)
                            rn = tmp;

                    }
                }
                tmp = tau * tau * rn;
                if (t > tmp)
                    indep = 1;
                else
                    indep = 0;
            }
            if (indep) {
                for (ll = 1; ll <= np1; ++ll) {
                    tmp = w[(krp1 - 1) * nw + ll];
                    w[(krp1 - 1) * nw + ll] = w[(ir - 1) * nw + ll];
                    w[(ir - 1) * nw + ll] = tmp;
                }
                tmp = scale[krp1];
                scale[krp1] = scale[ir];
                scale[ir] = tmp;
                itype[ir] = 0;
                t = sqrt(scale[krp1]);
                for (ll = 1; ll <= np1; ++ll)  
                    w[(krp1 - 1) * nw + ll] *= t;

                scale[krp1] = alsq;
                me = mep1;
                mep1 = me + 1;
                krank = krp1;
                krp1 = krank + 1;
            }
            else
                break;

            i = ip1;
            ip1++;
        }
    }
    if (krank < l) {
        for (i = 1; i <= krank; ++i) {
            j = krp1 - i;
            lhh12(ctx, 1,j,krp1,l,w+(j-1)*nw,1,h+j,w,1,nw,j-1);
        }
    }
    niv = krank + nsoln - l;
    niv1 = niv + 1;
    if (l == n)  
        done = 1;

    idope[1] = me;
    idope[2] = mep1;
    idope[3] = krank;
    idope[4] = krp1;
    idope[5] = nsoln;
    idope[6] = niv;
    idope[7] = niv1;
    idope[8] = l1;
    return(done);
}
 
/* ------------------------------------------------------------------------ */
/*  srotmg                                                                  */
/*      This function is part of the lsei/lsi package to solve constrained  */
/*      least square problems, adapted from: R.J. Hanson, K.H. Haskell,     */
/*      Algorithm 587 from Collected Algorithms of ACM.                     */

void srotmg(TDAContext *ctx, double *sd1,double *sd2,double *sx1,double *sy1,double *sparam)
{
    (void)ctx;        /* unused: the signature is shared */
    int iflag = 1;
    double sflag,sp1,sp2,su,sq1,sq2,stemp,sh11 = 0.0,sh12 = 0.0,sh21 = 0.0,sh22 = 0.0;
    double gam = 4096.E0;
    double gamsq = 1.678E7;
    double rgam = 2.441E-4;
    double rgamsq = 5.960E-8;

    if (*sd1 >= 0.0) {
        sp2 = *sd2 * *sy1;
        if (sp2 == 0.0) {
            sflag = -2.0;
            sparam[1] = sflag;
            return;
        }
        sp1 = *sd1 * *sx1;
        sq2 = sp2 * *sy1;
        sq1 = sp1 * *sx1;

        if (fabs(sq1) > fabs(sq2)) {
            sh21 = -(*sy1) / *sx1;
            sh12 =  sp2 / sp1;
            su = 1.0 - sh12 * sh21;
            if (su > 0.0) {
                sflag = 0.0;
                *sd1 = *sd1 / su;
                *sd2 = *sd2 / su;
                *sx1 = *sx1 * su;
                goto L100;
            }
        }
        else {
            if (sq2 >= 0.0) {
                sflag = 1.0;
                sh11 = sp1 / sp2;
                sh22 = *sx1 / *sy1;
                su = 1.0 + sh11 * sh22;
                stemp = *sd2 / su;
                *sd2 = *sd1 / su;
                *sd1 = stemp;
                *sx1 = *sy1 * su;
                goto L100;
            }
        }
    }
    sflag = -1.0;
    sh11 = sh12 = sh21 = sh22 = 0.0;
    *sd1  = *sd2  = *sx1  = 1.0;
    goto L220;
L100:
    if (iflag == 1) {
        rgam = 1.0 / gam;
        gamsq = gam * gam;
        rgamsq = rgam * rgam;
        iflag = 2;
    }
    while (*sd1 <= rgamsq) {
        if (*sd1 == 0.0)
            goto L160;

        if (sflag >= 0.0) {
            if (sflag == 0.0) {
                sh11  =  1.0;
                sh22  =  1.0;
                sflag = -1.0;
            }
            else {
                sh21  = -1.0;
                sh12  =  1.0;
                sflag = -1.0;
            }
        }
        *sd1 = *sd1 * gamsq;
        *sx1 = *sx1 * rgam;
        sh11 = sh11 * rgam;
        sh12 = sh12 * rgam;
    }
    while (*sd1 >= gamsq) {
        if (sflag >= 0.0) {
            if (sflag == 0.0) {
                sh11  =  1.0;
                sh22  =  1.0;
                sflag = -1.0;
            }
            else {
                sh21  = -1.0;
                sh12  =  1.0;
                sflag = -1.0;
            }
        }
        *sd1 = *sd1 * rgamsq;
        *sx1 = *sx1 * gam;
        sh11 = sh11 * gam;
        sh12 = sh12 * gam;
    }
L160:
    while (fabs(*sd2) <= rgamsq) {
        if (*sd2 == 0.0)
            goto L220;

        if (sflag >= 0.0) {
            if (sflag == 0.0) {
                sh11  =  1.0;
                sh22  =  1.0;
                sflag = -1.0;
            }
            else {
                sh21  = -1.0;
                sh12  =  1.0;
                sflag = -1.0;
            }
        }
        *sd2  = *sd2 * gamsq;
        sh21 = sh21 * rgam;
        sh22 = sh22 * rgam;
    }
    while (fabs(*sd2) >= gamsq) {
        if (sflag >= 0.0) {
            if (sflag == 0.0) {
                sh11  =  1.0;
                sh22  =  1.0;
                sflag = -1.0;
            }
            else {
                sh21  = -1.0;
                sh12  =  1.0;
                sflag = -1.0;
            }
        }
        *sd2  = *sd2 * rgamsq;
        sh21 = sh21 * gam;
        sh22 = sh22 * gam;
    }
L220:
    if (sflag < 0.0) {
        sparam[2] = sh11;
        sparam[3] = sh21;
        sparam[4] = sh12;
        sparam[5] = sh22;
    }
    else if (sflag > 0.0) {
        sparam[2] = sh11;
        sparam[5] = sh22;
    }
    else {
        sparam[3] = sh21;
        sparam[4] = sh12;
    }
    sparam[1] = sflag;
}

/* ------------------------------------------------------------------------ */
/*  srotm                                                                   */
/*      This function is part of the lsei/lsi package to solve constrained  */
/*      least square problems, adapted from: R.J. Hanson, K.H. Haskell,     */
/*      Algorithm 587 from Collected Algorithms of ACM.                     */

void srotm(TDAContext *ctx, int n,double *sx,double *sy,double *sparam)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;
    double sflag,sh11,sh12,sh21,sh22,w,z;

    sflag = sparam[1];
    if (n > 0 && sflag != -2.0) {
        if (sflag < 0.0) {
            sh11 = sparam[2];
            sh12 = sparam[4];
            sh21 = sparam[3];
            sh22 = sparam[5];
            for (i = 1; i <= n;++i) {
                w = sx[i];
                z = sy[i];
                sx[i] = w * sh11 + z * sh12;
                sy[i] = w * sh21 + z * sh22;
            }
        }
        else if (sflag > 0.0) {
            sh11 = sparam[2];
            sh22 = sparam[5];
            for (i = 1; i <= n;++i) {
                w = sx[i];
                z = sy[i];
                sx[i] =  w * sh11 + z;
                sy[i] = -w + sh22 * z;
            }
        }
        else {
            sh12 = sparam[4];
            sh21 = sparam[3];
            for (i = 1; i <= n;++i) {
                w = sx[i];
                z = sy[i];
                sx[i] = w + z * sh12;
                sy[i] = w * sh21 + z;
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  lhhfti                                                                  */
/*                                                                          */
/*      Solves a general linear system: a x = b in the least squares sense. */
/*                                                                          */
/*      The code is adapted from:                                           */
/*      C.L. Lawson, R.J. Hanson, Solving Least Squares Problems,           */
/*      Englewood Cliffs (Prentice-Hall) 1974                               */
/*                                                                          */
/*      Call: int lhhfti(m,n,na,nb,a,b,rnorm,mode,d,h,ip,tau)               */
/*                                                                          */
/*      a is a general (m,n) matrix, stored in a one-dimensional array,     */
/*      with a(i,j) = a[(i - 1) * na + j], na >= n.                         */
/*      b is a (m,nb) matrix, stored in a one-dimensional array             */
/*      with b(i,j) = b[(i - 1) * nb + j].                                  */
/*                                                                          */
/*      d, h, and ip are one-dimensional working arrays of length n.        */
/*                                                                          */
/*      Is is necessary that m >= nb.                                       */
/*                                                                          */
/*      If mode == 0, the function determines the pseudo rank of a, using   */
/*      the tolerance tau, and then returns this value. In this case, the   */
/*      b matrix is not used.                                               */
/*                                                                          */
/*      If mode == 1, the function calculates a solution of                 */
/*      a * x = b in the least square sense. The solution is returned in    */
/*      the field of b. The function returns again the pseudo rank of a.    */
/*                                                                          */
/*      If mode == 2, additional to mode == 1 the covariance matrix, i.e.   */
/*      the inverse of cross product matrix a'a is calculated and returned  */
/*      in the field of a. Note: this is done only if the pseudo rank of    */
/*      a is equal to n.                                                    */
/*                                                                          */
/*      If mode == 1 or 2, rnorm contains on return the norm of the         */
/*      least squares residuals.                                            */
/*                                                                          */

int lhhfti(TDAContext *ctx, int m,int n,int na,int nb,double *a,double *b,double *rnorm, int mode,double *d,double *h,int *ip,double tau)
{
    register int i,j,jb,l,ii;
    int jj,k,kp1,ldiag,lmax = 0,rk;
    double hmax,factor,sm,sm1,temp,tmp; 
  
    hmax = 0;
    factor = 0.001;
    rk = k = 0;  
    ldiag = n;
    if (m < n)
        ldiag = m;
    if (ldiag <= 0)
        return(0);

    for (j = 1; j <= ldiag; ++j) {

        /* Update squared column lengths and find lmax. */

        if (j > 1) {
            lmax = j;
            for (l = j; l <= n; ++l) {
                temp = a[(j - 2) * na + l];
                d[l] -= temp * temp;
                if (d[l] > d[lmax])
                    lmax = l;
            }
        }

        /* Compute squared column lengths and find lmax. */

        if (j == 1 || hmax + factor * d[lmax] <= hmax) {
            lmax = j;
            for (l = j; l <= n; ++l) {
                d[l] = 0.0;
                for (i = j; i <= m; ++i) {
                    temp = a[(i - 1) * na + l];
                    d[l] += temp * temp;
                }
               if (d[l] > d[lmax]) lmax = l;
            }
            hmax = d[lmax];
        }

        /* lmax has been determined. Do column interchanges if needed. */

        ip[j] = lmax;
        if (ip[j] != j) {
            for (i = 1; i <= m; ++i) {
                tmp = a[(i - 1) * na + j];
                a[(i - 1) * na + j] = a[(i - 1) * na + lmax];
                a[(i - 1) * na + lmax] = tmp;
            }
            d[lmax] = d[j];
        }

        /* Compute the j.th transformation and apply it to a[] and b[]. */

        lhh12(ctx, 1,j,j + 1,m,a + j - 1,na,d + j,a + j,na,1,n - j);
        lhh12(ctx, 2,j,j + 1,m,a + j - 1,na,d + j,b,nb,1,nb);
    }

    /*  Determine the pseudo rank rk, using the tolerance tau  */

    k = ldiag;
    for (j = 1; j <= ldiag; ++j) {
        if (fabs(a[(j - 1) * na + j]) <= tau) {
            k = j - 1;
            break;
        }
    }
    rk = k;
   
    if (mode > 0) {

        /* Compute the norms of the residual vectors. */

        kp1 = k + 1;
        for (jb = 1; jb <= nb; ++jb) {
            tmp = 0.0;
            for (i = kp1; i <= m; ++i) {
                temp = b[(i - 1) * nb + jb];
                tmp += temp * temp;
            }
            rnorm[jb - 1] = sqrt(tmp);
        }

        /* If the pseudo rank equal (or less than) 0, set b[] = 0.0 */

        if (k <= 0) {
            for (jb = 1; jb <= nb; ++jb)  
                for (i = 1; i <= n; ++i) 
                    b[(i - 1) * nb + jb] = 0.0;
        }
        else {

            /* If the pseudo rank is less than n, compute householder
               decomposition of first k rows. */

            if (k != n) {
                for (ii = 1; ii <= k; ++ii) {
                    i = kp1 - ii;
                    lhh12(ctx, 1,i,kp1,n,a + (i - 1) * na,1,h + i,a,1,na,i - 1);
                }
            }

            /* Now do for all columns of the matrix b[] */

            for (jb = 1; jb <= nb; ++jb) {

                /* Solve the k by k triangular system. */

                for (l = 1; l <= k; ++l) {
                    sm = 0.0; i = kp1 - l;
                    if (i != k) {
                        for (j = i + 1; j <= k; ++j)
                            sm += a[(i - 1) * na + j] * b[(j - 1) * nb + jb];
                    }
                    sm1 = sm;
                    temp = b[(i - 1) * nb + jb];
                    b[(i - 1) * nb + jb] = (temp - sm1) / a[(i - 1) * na + i];
                }

                /* Complete computation of solution vector. */

                if (k != n) {
                    for (j = kp1; j <= n; ++j)
                        b[(j - 1) * nb + jb] = 0.0;
                    for (i = 1; i <= k; ++i) 
                    lhh12(ctx, 2,i,kp1,n,a + (i - 1) * na,1,h + i,b + jb - 1,nb,1,1);
                }
   
                /* Reorder the solution vector to compensate for the
                   column interchange. */

                for (jj = 1; jj <= ldiag; ++jj) {
                    j = ldiag + 1 - jj;
                    if (ip[j] != j) {
                        l = ip[j];
                        tmp = b[(l - 1) * nb + jb];
                        b[(l - 1) * nb + jb] = b[(j - 1) * nb + jb];
                        b[(j - 1) * nb + jb] = tmp;
                    }
                }
            }
        }
    }
    if (mode == 2 && rk == n) {

        /* Compute the (unscaled) covariance matrix (i.e. the inverse of
           the cross product of the input matrix a[]. */

        for (j = 1; j <= n; ++j)
            a[(j - 1) * na + j] = 1.0 / a[(j - 1) * na + j];

        for (i = 1; i <= n - 1; ++i) {
            for (j = i + 1; j <= n; ++j) {
                sm = 0.0;
                for (l = i; l <= j - 1; ++l)  
                    sm += a[(i - 1) * na + l] * a[(l - 1) * na + j];
                a[(i - 1) * na + j] = -sm * a[(j - 1) * na + j];
            }
        }
        for (i = 1; i <= n; ++i) {
            for (j = i; j <= n; ++j) {
                sm = 0.0;
                for (l = j; l <= n; ++l)
                    sm += a[(i - 1) * na + l] * a[(j - 1) * na + l];
                a[(i - 1) * na + j] = sm;
            }
        }
        for (ii = 2; ii <= n; ++ii) {
            i = n + 1 - ii;
            if (ip[i] != i) {
                k = ip[i];
                tmp = a[(i - 1) * na + i];
                a[(i - 1) * na + i] = a[(k - 1) * na + k];
                a[(k - 1) * na + k] = tmp;
                for (l = 2; l <= i; ++l) {
                    tmp = a[(l - 2) * na + i];
                    a[(l - 2) * na + i] = a[(l - 2) * na + k];
                    a[(l - 2) * na + k] = tmp;
                }
                for (l = i + 1; l <= k - 1; ++l) {
                    tmp = a[(i - 1) * na + l];
                    a[(i - 1) * na + l] = a[(l - 1) * na + k];
                    a[(l - 1) * na + k] = tmp;
                }
                for (l = k + 1; l <= n; ++l) {
                    tmp = a[(i - 1) * na + l];
                    a[(i - 1) * na + l] = a[(k - 1) * na + l];
                    a[(k - 1) * na + l] = tmp;
                }
            }
        }
        for (i = 1; i <= n; ++i)
            for (j = i + 1; j <= n; ++j)
                a[(j - 1) * na + i] = a[(i - 1) * na + j];
    }
    return(rk);
}

/* ------------------------------------------------------------------------ */
/*  lhh12   Performs Householder transformations used by lhhfti, lsif,      */
/*          lsei, and lsi.                                                  */
/*                                                                          */
/*          The code is adapted from:                                       */
/*          C.L. Lawson, R.J. Hanson, Solving Least Squares Problems,       */
/*          Englewood Cliffs (Prentice-Hall) 1974                           */

void lhh12(TDAContext *ctx, int mode,int lpivot,int l1,int m,double *u,int iue,double *up, double *c,int ice,int icv,int ncv)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j,i2,i3,i4;
    int incr; double b,cl,clinv,sm,sm1,temp,tpiv;

    if (lpivot <= 0 || l1 <= lpivot || l1 > m)
        return;

    tpiv = u[1 + iue * (lpivot - 1)];
    cl = fabs(tpiv);

    if (mode == 2) {
        if (cl <= 0.0)
            return;
    }
    else {
        for (j = l1; j <= m; ++j) {
            temp = fabs(u[1 + iue * (j - 1)]);
            if (temp > cl) cl = temp;
        }
        if (cl <= 0.0) 
            return;

        clinv = 1.0 / cl;
        sm = tpiv * clinv * tpiv * clinv;
        for (j = l1; j <= m; ++j) {
            temp = u[1 + iue * (j - 1)];
            sm += temp * clinv * temp * clinv;
        }
        sm1 = sm; cl *= sqrt(sm1);
        if (tpiv > 0.0) cl = -cl;
        *up = tpiv - cl;
        tpiv = u[1 + iue * (lpivot - 1)] = cl;
    }
    if (ncv <= 0) 
        return;

    b = *up * tpiv;
    if (b >= 0.0)
        return;

    b = 1.0 / b;
    i2 = 1 - icv + ice * (lpivot - 1);
    incr = ice * (l1 - lpivot);
    for (j = 1; j <= ncv; ++j) {
        i2 += icv; i3 = i2 + incr; i4 = i3;
        sm = c[i2] * *up;
        for (i = l1; i <= m; ++i) {
            sm += c[i3] * u[1 + iue * (i - 1)];
            i3 += ice;
        }
        if (sm != 0.0) {
            sm *= b;
            c[i2] += sm * *up;
            for (i = l1; i <= m; ++i) {
                c[i4] += sm * u[1 + iue * (i - 1)];
                i4 += ice;
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  syminv(n,x)                                                             */
/*      Inversion of a symmetrical (positive definite) matrix x with        */
/*      dimension n. The inversion is performed directly in the lower       */
/*      triangular, then it is copied to the upper one. The function        */
/*      returns the rank of the matrix. For the determination of the rank   */
/*      the tolerance sqrt(EPSI) is used. Return -1 if error (alloc).       */

int syminv(TDAContext *ctx, int n, double *x)
{
    register int i,j,k,kk;
    int ii,jj,aflag,r;
    double tau,crit,tmp1,tmp2,tmp3,*wrk;

    tau = sqrt(ctx->EPSI);

    if (!(wrk = (double *)calloc((size_t)(n + 1),sizeof(double))))
        return(-1);
    memrq(ctx, n + 1,sizeof(double));
                
    r = n;      
    k = 1;      
    for (i = 1; i <= n; ++i) {          
        wrk[i] = x[k++];
        k += n;
    }
    for (k = 1; k <= n; ++k) {

        aflag = 0;                          /* flag for aliased parameters  */

        /* check for current sweep-column */

        kk = (k - 1) * n;
        if ((tmp1 = x[kk + k]) <= tau * wrk[k])   
            aflag = 1;

        else {
 
            j = 0;
            for (i = 1; i < k; ++i) {
                tmp2 = x[++j];

                if (tmp2 != 0.0) {
                    tmp3 = x[kk + i];
                    crit = 1.0 / (tmp2 + (tmp3 * tmp3) / tmp1);
                    if (crit < tau * wrk[i]) {
                        aflag = 1;
                        break;
                    }
                }
                j += n;
            }
        }

        if (!aflag) {

            /* collinearity check passed, so updating triangle */
 
            for (j = 1; j < k; ++j) {
                jj = (j - 1) * n;
                for (i = 1; i <= j; ++i)  
                    x[++jj] += x[kk + i] * x[kk + j] / tmp1;
            }
            for (i = k + 1; i <= n; ++i) {
                ii = (i - 1) * n;
                for (j = k + 1; j <= i; ++j)  
                    x[ii + j] -= x[ii + k] * x[(j - 1) * n + k] / tmp1;
            }
            for (i = k + 1; i <= n; ++i) {
                ii = (i - 1) * n;
                for (j = 1; j < k; ++j)  
                    x[ii + j] -= x[ii + k] * x[kk + j] / tmp1;
            }
            for (j = 1; j < k; ++j)
                x[kk + j] /= tmp1;
     
            for (i = k + 1; i <= n; ++i)  
                x[(i - 1) * n + k] /= -tmp1;

            x[kk + k] = 1.0 / tmp1;
        }
        else {
 
            /*  collinearity check not passed, so updating the k column     */
            /*  and k row with zeroes                                       */
 
            for (i = 1; i <= k; ++i)
                x[kk + i] = 0.0;

            for (i = k + 1; i <= n; ++i)  
                x[(i - 1) * n + k] = 0.0;
            r--;
        }
    }
    k = 1;                                 /* copy lower to upper triangle  */
    for (i = 1; i <= n; ++i) {          
        ii = (i - 1) * n;
        for (j = 1; j < i; ++j)  
            x[(j - 1) * n + i] = x[ii + j];
    }
    memrq(ctx, -n - 1,sizeof(double));
    free((char *)wrk);
    return(r);
}


