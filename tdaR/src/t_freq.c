/****************************************************************************/
/*  t_freq                                                                  */
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
#include "t_gdat.h"
#include "tda_context.h"

/*  functions in t_freq.c */

int cfreq(TDAContext *ctx, int n, int m, int *x, int mcel, int *value, int *freq, int *bf, int wflg, double *wt, double *wfreq);
int cfreq1(TDAContext *ctx, int n,short *idx,int mcel,int *value,int *freq,int *bf, int widx,double *wfreq);

/* ------------------------------------------------------------------------ */
/*  cfreq   Joint frequency distribution.                                   */
/*                                                                          */
/*      The code of this function is adapted from:                          */
/*      B.L. Leathers, Tabulating Sparse Joint Frequency Distributions,     */
/*      Algorithm AS 119, Applied Statistics 26 (1977), pp 364 - 368.       */
/*                                                                          */
/*      int cfreq(n,m,x,mcel,value,freq,bf,wflg,wt,wfreq)                   */
/*                                                                          */
/*      x is a field of input data with n columns and m rows. mcel is the   */
/*      max number of different categories. value is a field of n columns   */
/*      and mcel rows. It is used to build the combinations of categories   */
/*      in the contingency table. freq is a one-dimensional field with      */
/*      mcel entries, at return it contains the frequencies. bf is a        */
/*      one dimensional field of mcel entries, at return it contains        */
/*      pointers so that: value[bf[i],j] and freq[bf[i]] are sorted in      */
/*      ascending order.                                                    */
/*                                                                          */
/*      Return: if > 0 it is the number of combinations of categories       */
/*      in the contingency table, i.e. the effective number of rows in      */
/*      value, freq and bf. Return is -1 if mcel is too small, -2 if        */
/*      the (predefined) stack size is too small, or -3 if insufficient     */
/*      memory.                                                             */
/*                                                                          */
/*      Added following option:                                             */
/*                                                                          */
/*      int wflg        if 1 then on input                                  */
/*      double *wt      is a field of length m with weights, and on output  */
/*      double *wfreq   is, analogously to freq, a field containing         */
/*                      weighted frequencies.                               */
/*                                                                          */
/*      if wflg == 0, then wt and wfreq are not used.                       */
/*                                                                          */

int cfreq(TDAContext *ctx, int n, int m, int *x, int mcel, int *value, int *freq, int *bf, int wflg, double *wt, double *wfreq)
{
    register int i,j,k,nc;
    int ctype,lroot,nuroot,slevel,sub,subsub,mstk,lthis;
    int *lchild,*rchild,*stack,*test;
    int mm1,mm2,mm3,mm4;
    double w = 0.0;

    mm1 = mm2 = mm3 = mm4 = lroot = 0;

    mstk = (int)log((double)mcel);
    mstk = 4 * (mstk + 2);

    nc = -3;
    if (!(test = (int *) calloc((size_t)(n + 1),sizeof(int)))) { 
        p_err(ctx, -2,1);     
        goto FFin;
    }
    mm1 = n + 1;
    memrq(ctx, mm1,sizeof(int));

    if (!(lchild = (int *) calloc((size_t)(mcel + 1),sizeof(int)))) {
        p_err(ctx, -2,1);    
        goto FFin;
    }
    mm2 = mcel + 1;
    memrq(ctx, mm2,sizeof(int));

    if (!(rchild = (int *) calloc((size_t)(mcel + 1),sizeof(int)))) {
        p_err(ctx, -2,1);      
        goto FFin;
    }
    mm3 = mcel + 1;
    memrq(ctx, mm3,sizeof(int));

    if (!(stack = (int *) calloc((size_t)(mstk + 1),sizeof(int)))) {
        p_err(ctx, -2,1);     
        goto FFin;
    }
    mm4 = mstk + 1;
    memrq(ctx, mm4,sizeof(int));

    nc = 0;
    for (k = 1; k <= m; ++k) {

        for (i = 1; i <= n; ++i)  
            test[i] = x[(k - 1) * n + i];

        if (wflg)
            w = wt[k];

        lthis = lroot;
        slevel = 0;

FCont1:
        if (lthis != 0) {

            for (i = 1; i <= n; ++i) {
                ctype = test[i] - value[(lthis - 1) * n + i];
                if (ctype) {
                    if (ctype >= 0)
                        ctype = 1;
                    else
                        ctype = -1;
                    slevel += 2;
                    if (slevel > mstk) {
                        nc = -2;
                        goto FFin;
                    }
                    stack[slevel - 1] = ctype;
                    stack[slevel] = lthis;
         
                    if (ctype == -1)
                    lthis = lchild[lthis];
                    if (ctype == 1)
                        lthis = rchild[lthis];
                    goto FCont1;
                }
            }
            freq[lthis]++;
            if (wflg)
                wfreq[lthis] += w;

            goto FCont2;
        }
        if (nc >= mcel) {
            nc = -1;
            goto FFin;
        }
        nc++;
         
        for (i = 1; i <= n; ++i)  
            value[(nc - 1) * n + i] = test[i]; 
         
        freq[nc] = 1;
        if (wflg)
            wfreq[nc] = w;

        lchild[nc] = rchild[nc] = bf[nc] = 0;
         
        if (lroot) {
            lthis = stack[slevel];
            ctype = stack[slevel - 1];
            if (ctype == -1)
                lchild[lthis] = nc;
            if (ctype == 1)
                rchild[lthis] = nc;
        }
        else  
            lroot = nc;
         
        while (1) {
            if (slevel <= 0)
                goto FCont2;
         
            lthis = stack[slevel];
            ctype = stack[slevel - 1];
            slevel -= 2;
            if (bf[lthis])
                break;
         
            bf[lthis] = ctype;
        }
        if (ctype != bf[lthis]) {
            bf[lthis] = 0;
            goto FCont2;
        }
        if (ctype != -1) {  
            sub = rchild[lthis];
         
            if (bf[sub] != 1) {
                subsub = lchild[sub];
                nuroot = subsub;
                lchild[sub] = rchild[subsub];
                rchild[lthis] = lchild[subsub];
                lchild[subsub] = lthis;
                rchild[subsub] = sub;
                bf[lthis] = bf[sub] = 0;
         
                if (bf[subsub] == -1)
                    bf[sub] = 1;
                if (bf[subsub] == 1)
                    bf[lthis] = -1;
                bf[subsub] = 0;
            }
            else {
                nuroot = sub;
                rchild[lthis] = lchild[sub];
                lchild[sub] = lthis;
                bf[lthis] = bf[sub] = 0;
            }
        }
        else {
            sub = lchild[lthis];
         
            if (bf[sub] != -1) {
                subsub = rchild[sub];
                nuroot = subsub;
                rchild[sub] = lchild[subsub];
                lchild[lthis] = rchild[subsub];
                rchild[subsub] = lthis;
                lchild[subsub] = sub;
                bf[lthis] = bf[sub] = 0;
         
                if (bf[subsub] == -1)
                    bf[lthis] = 1;
                if (bf[subsub] == 1)
                    bf[sub] = -1;
                bf[subsub] = 0;
            }
            else {
                nuroot = sub;
                lchild[lthis] = rchild[sub];
                rchild[sub] = lthis;
                bf[lthis] = bf[sub] = 0;
            }
        }
        if (slevel) {
            lthis = stack[slevel];
            ctype = stack[slevel - 1];
            if (ctype == -1)
                lchild[lthis] = nuroot;
            if (ctype == 1)
                rchild[lthis] = nuroot;
        }
        else
            lroot = nuroot;
FCont2: ;
    }
    slevel = 0;
    j = 1;
    lthis = lroot;
    for (i = 1; i <= nc; ++i) {

        while (1) {
            if (!lthis)
                break;
            slevel++;
            stack[slevel] = lthis;
            lthis = lchild[lthis];
        }
        lthis = stack[slevel];
        slevel--;
        bf[j] = lthis;
        j++;
        lthis = rchild[lthis];
    }

FFin:
    if (mm4) {
        free((char *)stack);
        memrq(ctx, -mm4,sizeof(int));
    }  
    if (mm3) {
        free((char *)rchild);
        memrq(ctx, -mm3,sizeof(int));
    }   
    if (mm2) {
        free((char *)lchild);
        memrq(ctx, -mm2,sizeof(int));
    }
    if (mm1) {
        free((char *)test);
        memrq(ctx, -mm1,sizeof(int));
    }
    return(nc);
}

/* ------------------------------------------------------------------------ */
/*  cfreq1  Joint frequency distribution.                                   */
/*                                                                          */
/*      This function is very similar to cfreq. The only difference is      */
/*      that cfreq1 reads its input data directly from the data matrix.     */
/*                                                                          */
/*      int cfreq1(n,idx,mcel,value,freq,bf,widx,wfreq)                     */
/*                                                                          */
/*      int n           number of variables                                 */
/*      short *idx      internal numbers of variables idx[j], j=0,n-1.      */
/*      int mcel        max number of different combinations                */
/*      int *value      mcel x n array for different combinations           */
/*      int *freq       mcel x 1 array with frequencies                     */
/*      int *bf         mcel x 1 array with pointers for sorting            */
/*      int widx        if >=0 this is the internal number of a variable    */
/*                      used for weighting.                                 */
/*      double *wfreq   mcel x 1 array with weighted frequencies            */
/*                                                                          */
/*      On output, freq[bf[i]] (i=1,n), or wfreq[bf[i]], is the             */
/*      frequency of value[bf[i],j], with values in ascending order.        */
/*                                                                          */
/*      Return: if > 0 it is the number of combinations of categories       */
/*      in the contingency table, i.e. the effective number of rows in      */
/*      value, freq and bf. Return is -1 if mcel is too small, -2 if        */
/*      the (predefined) stack size is too small, or -3 if insufficient     */
/*      memory.                                                             */
/*  ###                                                                     */
/*      If widx < 0, the weighting variable, and wfreq, are not used.       */

int cfreq1(TDAContext *ctx, int n,short *idx,int mcel,int *value,int *freq,int *bf, int widx,double *wfreq)
{
    register int i,j,k,nc;
    int ctype,lroot,nuroot,slevel,sub,subsub,mstk,lthis;
    int *lchild,*rchild,*stack,*test;
    int mm1,mm2,mm3,mm4;
    double w = 0.0;

    lroot = mm1 = mm2 = mm3 = mm4 = 0;
   
    mstk = (int)log((double)mcel);
    mstk = 4 * (mstk + 2);
    nc = -3;

    if (!(test = (int *) calloc((size_t)(n + 1),sizeof(int)))) { 
        p_err(ctx, -2,1);    
        goto F1Fin;
    }
    mm1 = n + 1;
    memrq(ctx, mm1,sizeof(int));

    if (!(lchild = (int *) calloc((size_t)(mcel + 1),sizeof(int)))) {
        p_err(ctx, -2,1);    
        goto F1Fin;
    }
    mm2 = mcel + 1;
    memrq(ctx, mm2,sizeof(int));

    if (!(rchild = (int *) calloc((size_t)(mcel + 1),sizeof(int)))) {
        p_err(ctx, -2,1);    
        goto F1Fin;
    }
    mm3 = mcel + 1;
    memrq(ctx, mm3,sizeof(int));

    if (!(stack = (int *) calloc((size_t)(mstk + 1),sizeof(int)))) {
        p_err(ctx, -2,1);     
        goto F1Fin;
    }
    mm4 = mstk + 1;
    memrq(ctx, mm4,sizeof(int));

    nc = 0;
    for (k = 1; k <= ctx->NOC; ++k) {

        for (i = 1; i <= n; ++i)  
            test[i] = (int)get_data(ctx, (int)idx[i - 1],k - 1);
         
        if (widx >= 0)  
            w = get_data(ctx, widx,k - 1); 

        lthis = lroot;
        slevel = 0;

F1Cont1:
        if (lthis != 0) {

            for (i = 1; i <= n; ++i) {
                ctype = test[i] - value[(lthis - 1) * n + i];
                if (ctype) {
                    if (ctype >= 0)
                        ctype = 1;
                    else
                        ctype = -1;
                    slevel += 2;
                    if (slevel > mstk) {
                        nc = -2;
                        goto F1Fin;
                    }
                    stack[slevel - 1] = ctype;
                    stack[slevel] = lthis;
         
                    if (ctype == -1)
                    lthis = lchild[lthis];
                    if (ctype == 1)
                        lthis = rchild[lthis];
                    goto F1Cont1;
                }
            }
            freq[lthis]++;
            if (widx >= 0)
                wfreq[lthis] += w;

            goto FCont2;
        }
        if (nc >= mcel) {
            nc = -1;
            goto F1Fin;
        }
        nc++;
         
        for (i = 1; i <= n; ++i)  
            value[(nc - 1) * n + i] = test[i]; 
         
        freq[nc] = 1;
        if (widx >= 0)
            wfreq[nc] = w;

        lchild[nc] = rchild[nc] = bf[nc] = 0;
         
        if (lroot) {
            lthis = stack[slevel];
            ctype = stack[slevel - 1];
            if (ctype == -1)
                lchild[lthis] = nc;
            if (ctype == 1)
                rchild[lthis] = nc;
        }
        else  
            lroot = nc;
         
        while (1) {
            if (slevel <= 0)
                goto FCont2;
         
            lthis = stack[slevel];
            ctype = stack[slevel - 1];
            slevel -= 2;
            if (bf[lthis])
                break;
         
            bf[lthis] = ctype;
        }
        if (ctype != bf[lthis]) {
            bf[lthis] = 0;
            goto FCont2;
        }
        if (ctype != -1) {  
            sub = rchild[lthis];
         
            if (bf[sub] != 1) {
                subsub = lchild[sub];
                nuroot = subsub;
                lchild[sub] = rchild[subsub];
                rchild[lthis] = lchild[subsub];
                lchild[subsub] = lthis;
                rchild[subsub] = sub;
                bf[lthis] = bf[sub] = 0;
         
                if (bf[subsub] == -1)
                    bf[sub] = 1;
                if (bf[subsub] == 1)
                    bf[lthis] = -1;
                bf[subsub] = 0;
            }
            else {
                nuroot = sub;
                rchild[lthis] = lchild[sub];
                lchild[sub] = lthis;
                bf[lthis] = bf[sub] = 0;
            }
        }
        else {
            sub = lchild[lthis];
         
            if (bf[sub] != -1) {
                subsub = rchild[sub];
                nuroot = subsub;
                rchild[sub] = lchild[subsub];
                lchild[lthis] = rchild[subsub];
                rchild[subsub] = lthis;
                lchild[subsub] = sub;
                bf[lthis] = bf[sub] = 0;
         
                if (bf[subsub] == -1)
                    bf[lthis] = 1;
                if (bf[subsub] == 1)
                    bf[sub] = -1;
                bf[subsub] = 0;
            }
            else {
                nuroot = sub;
                lchild[lthis] = rchild[sub];
                rchild[sub] = lthis;
                bf[lthis] = bf[sub] = 0;
            }
        }
        if (slevel) {
            lthis = stack[slevel];
            ctype = stack[slevel - 1];
            if (ctype == -1)
                lchild[lthis] = nuroot;
            if (ctype == 1)
                rchild[lthis] = nuroot;
        }
        else
            lroot = nuroot;
FCont2: ;
    }
    slevel = 0;
    j = 1;
    lthis = lroot;
    for (i = 1; i <= nc; ++i) {

        while (1) {
            if (!lthis)
                break;
            slevel++;
            stack[slevel] = lthis;
            lthis = lchild[lthis];
        }
        lthis = stack[slevel];
        slevel--;
        bf[j] = lthis;
        j++;
        lthis = rchild[lthis];
    }

F1Fin:
    if (mm4) {
        free((char *)stack);
        memrq(ctx, -mm4,sizeof(int));
    }
    if (mm3) {
        free((char *)rchild);
        memrq(ctx, -mm3,sizeof(int));
    }
    if (mm2) {
        free((char *)lchild);
        memrq(ctx, -mm2,sizeof(int));
    }
    if (mm1) {
        free((char *)test);
        memrq(ctx, -mm1,sizeof(int));
    }  
    return(nc);
}

/*  end of t_freq.c */















    

