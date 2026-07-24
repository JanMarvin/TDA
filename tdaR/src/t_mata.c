/****************************************************************************/
/*  t_mata                                                                  */
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
#include "t_gf.h"
#include "t_alloc.h"
#include "t_mat.h"
#include "t_matc.h"
#include "tda_context.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_mata                                                     */

int m_rperm(TDAContext *ctx, int n,int *icn,int licn,int *ip,int *lenr,int *iperm, int *pr,int *arp,int *cv,int *out);
int m_perm(TDAContext *ctx, int n,int *icn,int licn,int *ip,int *lenr,int *arp,    int *ib,int *lowl,int *numb,int *prev);
int mpfit(TDAContext *ctx, int n,int m,double *a,double *u,double *v,double *b,int *iter,double *eps);
int mpfit_check(TDAContext *ctx, int n,int m,double *a,double *u,double *v);

/* ------------------------------------------------------------------------ */
/*  m_rperm() Row permutation for a zero-free diagonal. Adapted from        */
/*            CACM 575 (I.S. Duff).                                         */
/*                                                                          */
/*  n         order of matrix                                               */
/*  icn       array containing the column indices of the non-zeroes. Those  */
/*            belonging to a single row must be contiguous but the ordering */
/*            of column indices within each row is unimportant and wasted   */
/*            space between rows is permitted.                              */
/*  licn      length of array icn.                                          */
/*  ip        ip[i], i=1,2,...n, is the position in array icn of the first  */
/*            column index of a non-zero in row i.                          */
/*  lenr      lenr[i] is the number of non-zeros in row i, i=1,2,..n.       */
/*  iperm     contains permutation to make diagonal have the smallest       */
/*            number of zeros on it.  Elements (iperm[i],i) i=1, ... n are  */
/*            non-zero at the end of the algorithm unless matrix is         */
/*            structurally singular. In this case, (iperm[i],i) will be     */
/*            zero for n - numnz entries.                                   */
/*  numnz     number of non-zeros on diagonal of permuted matrix            */
/*                                                                          */
/*  pr, arp, cv and out are work arrays of length 1,...n.                   */
/*                                                                          */
/*  Return: numnz.                                                          */
  
int m_rperm(TDAContext *ctx, int n,int *icn,int licn,int *ip,int *lenr,int *iperm, int *pr,int *arp,int *cv,int *out)
{
    (void)ctx; (void)licn;        /* unused: the signature is shared */
    register int i,ii = 0,k,kk;
    int in1,in2 = 0,j,j1,jord,numnz;

    for (i = 1; i <= n; ++i) {
        arp[i] = lenr[i] - 1;
        cv[i] = 0;
        iperm[i] = 0;
    }
    numnz = 0;
 
    for (jord = 1; jord <= n; ++jord) {

        j = jord;
        pr[j] = -1;

        for (k = 1; k <= jord; ++k) {

            in1 = arp[j];
            if (in1 < 0)
                goto L60;             
            in2 = ip[j] + lenr[j] - 1;
            in1 = in2 - in1;

            for (ii = in1; ii <= in2; ++ii) {
                i = icn[ii];
                if (iperm[i] == 0)
                    goto L110;
            }
            arp[j] = -1;

L60:        out[j] = lenr[j] - 1;

            for (kk = 1; kk <= jord; ++kk) {

                in1 = out[j];
                if (in1 < 0)               
                    goto L80;

                in2 = ip[j] + lenr[j] - 1;
                in1 = in2 - in1;

                for (ii = in1; ii <= in2; ++ii) {
                    i = icn[ii];
                    if (cv[i] == jord)                
                        continue;

                    j1 = j;
                    j = iperm[i];
                    cv[i] = jord;
                    pr[j] = j1;
                    out[j1] = in2 - ii - 1;
                    goto L100; 
                }
L80:            j = pr[j];
                if (j == -1)                    
                    goto L130;
            }
L100:       ;
        }

L110:
        iperm[i] = j;
        arp[j] = in2 - ii - 1;
        numnz++;            

        for (k = 1; k <= jord; ++k) {
            j = pr[j];
            if (j == -1)
                goto L130;

            ii = ip[j] + lenr[j] - out[j] - 2;
            i = icn[ii];
            iperm[i] = j;
        }
L130:   ;
    }
    if (numnz == n)                
        goto L500;

    for (i = 1; i <= n; ++i)  
        arp[i] = 0;
    k = 0;
    for (i = 1; i <= n; ++i) {
        if (iperm[i] != 0) {
            j = iperm[i];
            arp[j] = i;
        }
        else {
            k++;
            out[k] = i;
        }
    }
    k = 0;
    for (i = 1; i <= n; ++i) {
        if (arp[i] != 0)                   
            continue;
        k++;      
        iperm[out[k]] = i;
    }
L500:
    return(numnz);
}

/* ------------------------------------------------------------------------ */
/*  m_perm()  Permutation to block triangular form. Adapted from            */
/*            CACM 529 (I.S. Duff and J.K. Reid)                            */
/*                                                                          */
/*  n         order of matrix                                               */
/*  icn       array containing the column indices of the non-zeroes. Those  */
/*            belonging to a single row must be contiguous but the ordering */
/*            of column indices within each row is unimportant and wasted   */
/*            space between rows is permitted.                              */
/*  licn      length of array icn.                                          */
/*  ip        ip[i], i=1,2,...n, is the position in array icn of the first  */
/*            column index of a non-zero in row i.                          */
/*  lenr      lenr[i] is the number of non-zeros in row i, i=1,2,..n.       */
/*  arp       ior[i] gives the position in the original ordering of the row */
/*            or column which is in position i in the permuted form, i=1,n  */
/*  ib        ib[i] is the row number in the permuted matrix of the         */ 
/*            beginning of block i, i=1,2,...num.                           */
/*  num       number of blocks found.                                       */
/*                                                                          */
/*  lowl, numb, prev are work arrays of length n.                           */
/*                                                                          */
/*  Return: num.                                                            */
  
int m_perm(TDAContext *ctx, int n,int *icn,int licn,int *ip,int *lenr,int *arp,    int *ib,int *lowl,int *numb,int *prev)
{
    (void)ctx; (void)licn;        /* unused: the signature is shared */
    register int i,j,k,ii;
    int num,icnt,lcnt,nnm1,iv,iw,isn,ist,ist1,i1,i2,dummy,stp;

    icnt = 0;
    num = 0;
    nnm1 = n + n - 1;

    for (j = 1; j <= n; ++j) {
        numb[j] = 0;
        arp[j] = lenr[j] - 1;
    }
    for (isn = 1; isn <= n; ++isn) {

        if (numb[isn] != 0)
            continue;                        
        iv = isn;
        ist = 1;
        lowl[iv] = 1;
        numb[iv] = 1;
        ib[n] = iv;
 
        for (dummy = 1; dummy <= nnm1; ++dummy) {

            i1 = arp[iv];
            if (i1 < 0)
                goto L30;

            i2 = ip[iv] + lenr[iv] - 1;
            i1 = i2 - i1;
 
            for (ii = i1; ii <= i2; ++ii) {
                iw = icn[ii];
                if (numb[iw] == 0)
                    goto L70;                       

                if (lowl[iw] < lowl[iv])
                    lowl[iv] = lowl[iw];
            }
            arp[iv] = -1;

L30:        if (lowl[iv] < numb[iv])
                goto L60;                          
 
            num++;             
            ist1 = n + 1 - ist;
            lcnt = icnt + 1;

            for (stp = ist1; stp <= n; ++stp) {
                iw = ib[stp];
                lowl[iw] = n + 1;
                icnt++;           
                numb[iw] = icnt;
                if (iw == iv)
                    break;                   
            }
            ist = n - stp;
            ib[num] = lcnt;

            if (ist != 0)
                goto L60;

            if (icnt < n)
                goto L90;               
            goto L100;
 
L60:        iw = iv;
            iv = prev[iv]; 

            if (lowl[iw] < lowl[iv])
                lowl[iv] = lowl[iw];
            goto L80;
 
L70:        arp[iv] = i2 - ii - 1;
            prev[iw] = iv;
            iv = iw;
            ist++;            
            lowl[iv] = ist;
            numb[iv] = ist;
            k = n + 1 - ist;
            ib[k] = iv;
L80:        ; 
        }
L90:    ;    
    }
L100:
    for (i = 1; i <= n; ++i) {
        ii = numb[i];
        arp[ii] = i;
    }
    return(num);
}

/* ------------------------------------------------------------------------ */
/*  mpfit(n,m,a,u,v,b,iter,eps)  Iterative proportional fitting.            */
/*                                                                          */
/*  a is the input table (n x m). u is the new column sums (1 x m), c is    */
/*  the new row sums (n x 1). b (n,m) is to be estimated.                   */  
/*                                                                          */
/*  The input table a must not contain negative entries. It is also         */
/*  required that all row and column sums are positive. The same is         */
/*  required for the new row and new column sums. Furthermore, the          */
/*  sum of prescribed row sums and the sum of the prescribed column sums    */
/*  should be equal. Otherwise the command will issue a warnings messsage.  */
/*                                                                          */  
/*  The algorithm performs a maximum of iter iterations. It is terminated   */  
/*  when the maximum deviation between the estimated and prescribed row     */
/*  and column sums is not greater than eps.                                */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */
  
int mpfit(TDAContext *ctx, int n,int m,double *a,double *u,double *v,double *b,int *iter,double *eps)
{
    register int i,j;
    int it;
    double dev = 0.0,s,tmp;

    if (mpfit_check(ctx, n,m,a,u,v))
        return(-1);

    for (i = 1; i <= n; ++i) {          /* set starting values */
        for (j = 1; j <= m; ++j)  
            b[(i - 1) * m + j] = a[(i - 1) * m + j];
    }
              
    for (it = 1; it <= *iter; ++it) {
        for (j = 1; j <= m; ++j) {
            s = 0.0;
            for (i = 1; i <= n; ++i)
                s += b[(i - 1) * m + j];

            s = u[j] / s;
            for (i = 1; i <= n; ++i)  
                b[(i - 1) * m + j] *= s;
        }
        for (i = 1; i <= n; ++i) {
            s = 0.0;
            for (j = 1; j <= m; ++j)
                s += b[(i - 1) * m + j];

            s = v[i] / s;
            for (j = 1; j <= m; ++j)   
                b[(i - 1) * m + j] *= s;
        }
        dev = 0.0;
        for (i = 1; i <= n; ++i) {
            tmp = 0.0;
            for (j = 1; j <= m; ++j)
                tmp += b[(i - 1) * m + j];
            dev = dmax(ctx, dev,fabs(tmp - v[i]));
        }
        for (j = 1; j <= m; ++j) {
            tmp = 0.0;
            for (i = 1; i <= n; ++i)
                tmp += b[(i - 1) * m + j];
            dev = dmax(ctx, dev,fabs(tmp - u[j]));
        }
        if (dev <= *eps)
            break;

    }
    *iter = it;
    *eps  = dev;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mpfit_check(n,m,a,u,v)  Check for correct input values.                 */
/*                          Return 0 if OK, -1 if error.                    */
 
int mpfit_check(TDAContext *ctx, int n,int m,double *a,double *u,double *v)
{
    register int i,j;
    double s,aij,usum,vsum;

    usum = vsum = 0.0;
    for (i = 1; i <= n; ++i) {

        if (v[i] <= 0.0) {
            printf1(ctx, "Error: requested row sum %d is negative or zero.\n",i);
            return(-1);
        }
        vsum += v[i];

        s = 0.0;
        for (j = 1; j <= m; ++j) {
            if ((aij = a[(i - 1) * m + j]) < -ctx->EPSI1) {
                printf1(ctx, "Error: negative coefficient in table (%d,%d).\n",i,j);
                return(-1);
            }
            s += aij;
        }
        if (s <= 0.0) {
            printf1(ctx, "Error: column sum %d is zero.\n",i);
            return(-1);
        }
    }
    for (j = 1; j <= m; ++j) {

        if (u[j] <= 0.0) {
            printf1(ctx, "Error: requested column sum %d is negative or zero.\n",j);
            return(-1);
        }
        usum += u[j];
        s = 0.0;
        for (i = 1; i <= n; ++i) {
            if ((aij = a[(i - 1) * m + j]) > 0.0)
                s += aij;
        }
        if (s <= 0.0) {
            printf1(ctx, "Error: row sum %d is zero.\n",j);
            return(-1);
        }
    }
    if (fabs(usum - vsum) > ctx->EPSI1)  
        printf1(ctx, "Warning: unequal sums of marginals (%g, %g).\n",usum,vsum); 
    return(0);
}



