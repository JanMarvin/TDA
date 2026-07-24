/****************************************************************************/
/*  t_ds                                                                    */
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
#include "t_var.h"
#include "t_gdat.h"
#include "t_sort.h"
#include "t_cdf.h"
#include "t_gf.h"

/*  functions in t_ds.c */

int mean(int n,short *vinum,double *m);
int std(int n,short *vinum,double *m,double *s,int opt);
int corr(int n,short *vinum,double *m,double *s,double *x);
void rcorr(int n,short *vinum,double *x);
int cov(int n,short *vinum,double *m,double *x);
double xmean(int n,double *x,double *wt);
double xstd(int n,double *x,double *wt,double *std);
double quantf(int n,double *x,double q);

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */

float *DE_X;            /* variable for density estimation                  */

/*--------------------------------------------------------------------------*/
/*  mean(n,vinum,m)     Calculate mean of variables vinum[i] (i=0,n-1)      */
/*                      Return values in m[i].                              */
/*                      If defined use case weights.                        */
/*  Return 0 if OK, -1 if error.                                            */

int mean(int n,short *vinum,double *m)
{
    register int i,j;
    double wt,wsum;

    wt = 1.0;
    if (WIVar >= 0)
        wsum = 0.0;
    else
        wsum = (double)NOC;

    for (j = 0; j < n; ++j)
        m[j] = 0.0;

    for (i = 0; i < NOC; ++i) {

        if (WIVar >= 0) {
            wt = get_data(WIVar,i) * WNorm;
            wsum += wt;
        }
        for (j = 0; j < n; ++j) {
            if (WIVar >= 0)
                m[j] += get_data(vinum[j],i) * wt;
            else
                m[j] += get_data(vinum[j],i);
        }
    }
    for (j = 0; j < n; ++j) {
        if (wsum > 0.0)
            m[j] /= wsum;
        else
            m[j] = 0.0;
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  std(n,vinum,m,s,opt)    Calculate mean for variables vinum[i] in m[]    */
/*                          and std deviation in s[].                       */
/*                          Use case weights if defined.                    */
/*                          if opt, do not divide by n - 1.                 */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int std(int n,short *vinum,double *m,double *s,int opt)
{
    register int i,j;
    double wt,tmp,wsum;

    wt = 1.0;
    if (WIVar >= 0)
        wsum = 0.0;
    else
        wsum = (double)NOC;

    if (mean(n,vinum,m))  
        return(-1);

    for (j = 0; j < n; ++j)
        s[j] = 0.0;

    for (i = 0; i < NOC; ++i) {

        if (WIVar >= 0) {
            wt = get_data(WIVar,i) * WNorm;
            wsum += wt;
        }
        for (j = 0; j < n; ++j) {
            tmp = get_data(vinum[j],i) - m[j];
            tmp *= tmp;
            if (WIVar >= 0)
                s[j] += tmp * wt;
            else
                s[j] += tmp;
        }
    }
    if (wsum <= 1.0) {
        for (j = 0; j < n; ++j)
            s[j] = 0.0;
    }
    else {
        if (opt == 0) {
            wsum -= 1.0;
            for (j = 0; j < n; ++j)
                s[j] /= tmp;
        } 
        for (j = 0; j < n; ++j) {
            tmp = s[j];
            if (tmp > 0.0)
                s[j] = sqrt(tmp);
        }
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  corr(n,vinum,m,s,x)                                                     */
/*                                                                          */
/*  calculate mean of variables vinum[i] (i = 0,n-1) in m[], standard       */
/*  deviation in s[], and lower triangle, including main diagonal, of       */
/*  corr matrix in x[].                                                     */
/*                                                                          */
/*  Return 0 if OK, -1 if error, -2 if insufficient memory.                 */

int corr(int n,short *vinum,double *m,double *s,double *x)
{
    register int i,j,k,l;
    double wt,xy,tmp;

    if (std(n,vinum,m,s,1))  
        return(-1);
      
    wt = 1.0;
    l = 0;
    for (j = 0; j < n; ++j) {
        for (k = 0; k <= j; ++k) {

            if (k == j) {
                x[l++] = 1.0;
                continue;
            }
            xy = 0.0;
            for (i = 0; i < NOC; ++i) {

                tmp = (get_data(vinum[j],i) - m[j]) *
                      (get_data(vinum[k],i) - m[k]);  

                if (WIVar >= 0) {
                    wt = get_data(WIVar,i) * WNorm;
                    xy += wt * tmp;
                }
                else
                    xy += tmp;
            }
            tmp = s[j] * s[k];
            if (tmp <= 0.0)
                x[l] = -2.0;
            else
                x[l] = xy / tmp;
            l++;
        }
    }
    return(0); 
}

/*--##----------------------------------------------------------------------*/
/*  rcorr(n,vinum,x)                                                        */
/*                                                                          */
/*  return rank correlation (tau) of variables with indices vinum[i],       */
/*  i = 0,...,n-1, in the lower triangle of x[], including main diagonal.   */

void rcorr(int n,short *vinum,double *x)
{
    register int ii,jj,l,j,k;
    double nn,xy,xi,yi,xj,yj;

    nn = (double)(NOC * (NOC - 1) / 2);

    l = 0;
    for (j = 0; j < n; ++j) {
        for (k = 0; k <= j; ++k) {

            if (k == j) {
                x[l++] = 1.0;
                continue;
            }
            xy = 0.0;
            for (ii = 1; ii < NOC; ++ii) {

                xi = get_data(vinum[j],ii);            
                yi = get_data(vinum[k],ii);            

                for (jj = 0; jj < ii; ++jj) {

                    xj = get_data(vinum[j],jj);            
                    yj = get_data(vinum[k],jj);            

                    xy += dsign(xi - xj) * dsign(yi - yj);
                }
            }
            x[l] = xy / nn;
            l++;
        }
    }
    return; 
}

/*--------------------------------------------------------------------------*/
/*  cov(n,vinum,m,x)                                                        */
/*                                                                          */
/*  calculate mean of variables vinum[i] (i = 0,n-1) in m[], and            */
/*  lower triangle, including diagonal, of cov matrix in x[].               */
/*                                                                          */
/*  Return 0 if OK, -1 if error                                             */

int cov(int n,short *vinum,double *m,double *x)
{
    register int i,j,k,l;
    double wt,xy,tmp,f;

    if (mean(n,vinum,m))  
        return(-1);
      
    if (NOC == 1)
        f = 1.0;
    else
        f = (double)(NOC - 1);

    l = 0;
    for (j = 0; j < n; ++j) {
        for (k = 0; k <= j; ++k) {

            xy = 0.0;
            for (i = 0; i < NOC; ++i) {

                tmp = (get_data(vinum[j],i) - m[j]) *
                      (get_data(vinum[k],i) - m[k]);  

                if (WIVar >= 0) {
                    wt = get_data(WIVar,i) * WNorm;
                    xy += wt * tmp;
                }
                else
                    xy += tmp;
            }
            x[l++] = xy / f;
        }
    }
    return(0); 
}

/*--------------------------------------------------------------------------*/
/*  xmean(n,x,wt)   Calculate mean of x[i] i = 0,...,n-1 with weights wt.   */    
/*                  Return: mean value.                                     */

double xmean(int n,double *x,double *wt)
{
    register int i;
    double wsum,mean;

    wsum = mean = 0.0;

    for (i = 0; i < n; ++i) {
        mean += x[i] * wt[i];
        wsum += wt[i];
    }
    if (wsum > 0.0)
        mean /= wsum;
    else
        mean = 0.0;
    return(mean);
}

/*--------------------------------------------------------------------------*/
/*  xstd(n,x,wt,std)    Calculate mean of x[i] (i = 0,...,n-1) and          */
/*                      standard deviation, using weights wt[].             */
/*                      REturn mean, std. deviation in std.                 */  

double xstd(int n,double *x,double *wt,double *std)
{
    register int i;
    double tmp,wsum,mean,s;

    wsum = mean = 0.0;

    for (i = 0; i < n; ++i) {
        mean += x[i] * wt[i];
        wsum += wt[i];
    }
    if (wsum <= 0.0) {
        mean = 0.0;
        *std = 0.0;
    }     
    else {
        mean /= wsum;
        s = 0.0;
        for (i = 0; i < n; ++i) {
            tmp = x[i] - mean;
            s += tmp * tmp * wt[i];
        }
        if (wsum > 1.0)
            s /= (wsum - 1.0);
        if (s > 0.0)
            s = sqrt(s);
        *std = s;
    }
    return(mean);
}
   
/* ------------------------------------------------------------------------ */
/*  quantf(n,x,q)   return q-quantile of x[0,...,n-1], assuming that x[]    */
/*                  is sorted in ascending order.                           */

double quantf(int n,double *x,double q)
{
    double qn,qnf,tmp;

    if (n < 1)
        return(0.0);

    if (n == 1)
        return(x[0]);

    qn = q * (double)(n + 1);
    if (qn <= 1.0)
        return(x[0]);
    else if (qn >= (double)n)
        return(x[n - 1]);

    qnf = floor(qn);
    tmp = (1.0 - (qn - qnf)) * x[(int)qnf - 1] + (qn - qnf) * x[(int)qnf];
    return(tmp);
}


