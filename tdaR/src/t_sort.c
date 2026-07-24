/****************************************************************************/
/*  t_sort                                                                  */
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
#include "t_gdat.h"
#include "t_var.h"
#include "t_gm.h"
#include "tda_context.h"
#include "tda_compat.h"
#include <stdlib.h>
#include <string.h>

static void tda_merge(char *a, char *tmp, size_t lo, size_t mid, size_t hi,
                      size_t sz, int (*cmp)(const void *, const void *, void *),
                      void *thunk)
{
    size_t i = lo, j = mid, k = lo;

    while (i < mid && j < hi) {
        if (cmp(a + j * sz, a + i * sz, thunk) < 0) {
            memcpy(tmp + k * sz, a + j * sz, sz); ++j;
        } else {
            memcpy(tmp + k * sz, a + i * sz, sz); ++i;
        }
        ++k;
    }
    if (i < mid) memcpy(tmp + k * sz, a + i * sz, (mid - i) * sz);
    if (j < hi)  memcpy(tmp + k * sz, a + j * sz, (hi - j) * sz);
    memcpy(a + lo * sz, tmp + lo * sz, (hi - lo) * sz);
}

static void tda_msort(char *a, char *tmp, size_t lo, size_t hi, size_t sz,
                      int (*cmp)(const void *, const void *, void *), void *thunk)
{
    size_t mid;

    if (hi - lo < 2) return;
    mid = lo + (hi - lo) / 2;
    tda_msort(a, tmp, lo, mid, sz, cmp, thunk);
    tda_msort(a, tmp, mid, hi, sz, cmp, thunk);
    if (cmp(a + mid * sz, a + (mid - 1) * sz, thunk) < 0)
        tda_merge(a, tmp, lo, mid, hi, sz, cmp, thunk);
}

void tda_qsort_r(void *base, size_t n, size_t sz,
                 int (*cmp)(const void *, const void *, void *), void *thunk)
{
    char *tmp;

    if (n < 2 || sz == 0) return;
    tmp = malloc(n * sz);
    if (tmp) {
        tda_msort((char *)base, tmp, 0, n, sz, cmp, thunk);
        free(tmp);
        return;
    }
    /* out of memory: stable insertion sort by adjacent swaps */
    {
        char *a = (char *)base;
        size_t i, j, b;
        for (i = 1; i < n; ++i)
            for (j = i; j > 0 && cmp(a + j * sz, a + (j - 1) * sz, thunk) < 0; --j)
                for (b = 0; b < sz; ++b) {
                    char c = a[j * sz + b];
                    a[j * sz + b] = a[(j - 1) * sz + b];
                    a[(j - 1) * sz + b] = c;
                }
    }
}

/*  functions in t_sort.c */

int sortd(TDAContext *ctx, int n,double *x,int opt);
int sortd2(TDAContext *ctx, int n,double *x,double *y);
int sortd3(TDAContext *ctx, int n,double *x,double *y,double *z);
int s1comp(const void *, const void *, void *);
int vsort(TDAContext *ctx, int n,short *vidx,int opt,int cflg,int pflg);   
int dvcomp(const void *, const void *, void *);
int sorti(TDAContext *ctx, int n,int *x,int opt);
int sorti2(TDAContext *ctx, int n,int *x,int *y);
int sicomp(const void *, const void *, void *);
int sortdpi(TDAContext *ctx, int n,int *x,int *ptr);
int sortdpi2(TDAContext *ctx, int n,int *x,int *y,int *ptr);
int si2comp(const void *, const void *, void *);
int sortdp(TDAContext *ctx, int n,double *x,int *ptr);
int sortdp2(TDAContext *ctx, int n,double *x,double *y,int *ptr);
int s2comp(const void *, const void *, void *);
int sortdp2f(TDAContext *ctx, int n,float *x,float *y,int *ptr);
int s2fcomp(const void *, const void *, void *);
int sortdp4f(TDAContext *ctx, int n,float *xa,float *ya,float *xb,float *yb,int *ptr);
int s4fcomp(const void *, const void *, void *);
int sortdp2a(TDAContext *ctx, int n,double *x,short *c,int *ptr,int opt);
int s2acomp(const void *, const void *, void *);
int s2bcomp(const void *, const void *, void *);
int sortdpn(TDAContext *ctx, int m,int n,double *x,int nc,int *col,int *ptr);
int sncomp(const void *, const void *, void *);
int sortdpn1(TDAContext *ctx, int m,int n,double *x,int nc,int *col,int *ptr,short *cen);
int sn1comp(const void *, const void *, void *);
int sortdp2c(TDAContext *ctx, int n,double *x,double *y,int *ptr,double xmin,double ymin, double xmax,double ymax);
int s2ccomp(const void *, const void *, void *);


/* ------------------------------------------------------------------------ */
/*  sortd(n,x,opt)  sort x[] containing n values in ascending order.        */
/*                  if opt=0 ascending order, else descending.              */
/*                  Return 0 if OK, otherwise -1 (insuff memory).           */

int sortd(TDAContext *ctx, int n,double *x,int opt)
{
    register int i,j;
    int err,ptr_a,srt_a;
    int *ptr;

    if (n < 2)
        return(0);
    ptr_a = srt_a = 0;
    err = -1;

    /* allocate a pointer for sorting */

    if (!(ptr = (int *)calloc((size_t)(n),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto SRTDFin;
    }
    memrq(ctx, n,sizeof(int));
    ptr_a = 1;

    if (!(ctx->SORTV = (double *)calloc((size_t)(n),sizeof(double)))) { 
        p_err(ctx, -2,1);
        goto SRTDFin;
    }
    memrq(ctx, n,sizeof(double));
    srt_a = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        ctx->SORTV[i] = x[i];
    }

    tda_qsort_r((char *)ptr,(size_t)(n),sizeof(int), s1comp, ctx);

    if (opt == 0) {
        for (i = 0; i < n; ++i)  
            x[i] = ctx->SORTV[ptr[i]];
    }
    else {
        j = n;  
        for (i = 0; i < n; ++i)  
            x[--j] = ctx->SORTV[ptr[i]];
    }
    err = 0;

SRTDFin:
    if (srt_a) {
        free((char *)ctx->SORTV);
        memrq(ctx, -n,sizeof(double));
    }
    if (ptr_a) {
        free((char *)ptr);
        memrq(ctx, -n,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sortd2(n,x,y)       sort x[] containing n values in ascending order.    */
/*                      change values in y[] in the same order.             */
/*                      Return 0 if OK, otherwise -1 (insuff memory).       */

int sortd2(TDAContext *ctx, int n,double *x,double *y)
{
    register int i;
    int err,ptr_a,srt_a;
    int *ptr;

    if (n < 2)
        return(0);
    ptr_a = srt_a = 0;
    err = -1;

    /* allocate a pointer for sorting */

    if (!(ptr = (int *)calloc((size_t)(n),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto SRTD2Fin;
    }
    memrq(ctx, n,sizeof(int));
    ptr_a = 1;

    if (!(ctx->SORTV = (double *)calloc((size_t)(n),sizeof(double)))) { 
        p_err(ctx, -2,1);
        goto SRTD2Fin;
    }
    memrq(ctx, n,sizeof(double));
    srt_a = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        ctx->SORTV[i] = x[i];
    }

    tda_qsort_r((char *)ptr,(size_t)(n),sizeof(int), s1comp, ctx);

    for (i = 0; i < n; ++i)  
        x[i] = ctx->SORTV[ptr[i]];
    for (i = 0; i < n; ++i)  
        ctx->SORTV[i] = y[ptr[i]];
    for (i = 0; i < n; ++i)  
        y[i] = ctx->SORTV[i];               

    err = 0;

SRTD2Fin:
    if (srt_a) {
        free((char *)ctx->SORTV);
        memrq(ctx, -n,sizeof(double));
    }
    if (ptr_a) {
        free((char *)ptr);
        memrq(ctx, -n,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sortd3(n,x,y,z)     sort x[] containing n values in ascending order.    */
/*                      change values in y[] and z[] in the same order.     */
/*                      Return 0 if OK, otherwise -1 (insuff memory).       */

int sortd3(TDAContext *ctx, int n,double *x,double *y,double *z)
{
    register int i;
    int err,ptr_a,srt_a;
    int *ptr;

    if (n < 2)
        return(0);
    ptr_a = srt_a = 0;
    err = -1;

    /* allocate a pointer for sorting */

    if (!(ptr = (int *)calloc((size_t)(n),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto SRTD3Fin;
    }
    memrq(ctx, n,sizeof(int));
    ptr_a = 1;

    if (!(ctx->SORTV = (double *)calloc((size_t)(n),sizeof(double)))) { 
        p_err(ctx, -2,1);
        goto SRTD3Fin;
    }
    memrq(ctx, n,sizeof(double));
    srt_a = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        ctx->SORTV[i] = x[i];
    }

    tda_qsort_r((char *)ptr,(size_t)(n),sizeof(int), s1comp, ctx);

    for (i = 0; i < n; ++i)  
        x[i] = ctx->SORTV[ptr[i]];
    for (i = 0; i < n; ++i)  
        ctx->SORTV[i] = y[ptr[i]];
    for (i = 0; i < n; ++i)  
        y[i] = ctx->SORTV[i];               
    for (i = 0; i < n; ++i)  
        ctx->SORTV[i] = z[ptr[i]];
    for (i = 0; i < n; ++i)  
        z[i] = ctx->SORTV[i];               

    err = 0;

SRTD3Fin:
    if (srt_a) {
        free((char *)ctx->SORTV);
        memrq(ctx, -n,sizeof(double));
    }
    if (ptr_a) {
        free((char *)ptr);
        memrq(ctx, -n,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  s1comp()    compare function                                            */

int s1comp(const void *arg1, const void *arg2, void *_ctx)
{
    TDAContext *ctx = (TDAContext *)_ctx;     
    double x; 

    x = ctx->SORTV[*(int *)arg1] - ctx->SORTV[*(int *)arg2];
    if (x > 0.0)
        return(1);
    else if (x < 0.0)
        return(-1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  vsort(n,vidx,opt,cflg,pflg)                                             */
/*                                                                          */
/*          If opt != 0 sort data matrix with respect to variables in       */
/*          vidx[i] (i=0,n-1) in ascending oder. Create pointer             */
/*          VSORTPtr[j] (j=0,NOC-1) for sorted cases. If opt = 0 free       */
/*          previously allocated memory. If cflg = 1 check whether data     */
/*          are already sorted. If pflg != 0 print info.                    */
/*                                                                          */
/*  Return:  0 if OK                                                        */
/*          -1 if insuff memory                                             */

int vsort(TDAContext *ctx, int n,short *vidx,int opt,int cflg,int pflg)
{     
    register int i,j,k;
    int err,r;
    double tmp;

    err = -1;
    if (opt == 0) {
        err = 0;
        goto VSORTFin;
    }
    if (pflg) {
        printf1(ctx, "Sort with variable(s): %s",ctx->VName[vidx[0]]);
        for (i = 1; i < n; ++i)  
            printf1(ctx, ",%s",ctx->VName[vidx[i]]);
        printf1(ctx, "\n");
        prn_mem(ctx);
    }
    r = 0;
    for (i = 0; i < n; ++i) {
        if (ctx->VTyp[vidx[i]] == 1) {
            r = 1;
            break;
        }
    }
    if (r) {
        printf1(ctx, "Error: can't sort string variables.\n");
        goto VSORTFin;
    }

    if (n < 1)
        goto VSORTFin;
    if (!(ctx->VSORTIdx = (int *)calloc((size_t)(n),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto VSORTFin;
    }
    memrq(ctx, n,sizeof(int));
    ctx->VSORTIdxA = ctx->VSORTIdxN = n;
    for (i = 0; i < n; ++i)
        ctx->VSORTIdx[i] = vidx[i];

    if (!(ctx->VSORTPtr = (int *)calloc((size_t)(ctx->NOC),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto VSORTFin;
    }
    memrq(ctx, ctx->NOC,sizeof(int));
    ctx->VSORTPtrA = ctx->NOC;

    for (i = 0; i < ctx->NOC; ++i)
        ctx->VSORTPtr[i] = i;

    if (cflg) {
        for (i = 1; i < ctx->NOC; ++i) {

            for (k = 0; k < ctx->VSORTIdxN; ++k) {
                j = ctx->VSORTIdx[k];
                tmp = (double)((float)(get_data(ctx, j,i) - get_data(ctx, j,i - 1)));
                if (tmp > 0.0)
                    break;
                else if (tmp < 0.0) {
                    cflg = 0;
                    break;
                }
            }         
            if (cflg == 0)
                break;
        }
        if (cflg) {
            if (pflg)
                printf1(ctx, "Data already sorted.\n");
            return(0);
        }
        if (pflg)
            printf1(ctx, "Data not already sorted, see records %d and %d.\n",i,i + 1);
    }
    tda_qsort_r((char *)ctx->VSORTPtr,(size_t)(ctx->NOC),sizeof(int), dvcomp, ctx);
    return(0);

VSORTFin:
    if (ctx->VSORTIdxA > 0) {
        free((char *)ctx->VSORTIdx);
        memrq(ctx, -ctx->VSORTIdxA,sizeof(int));
        ctx->VSORTIdxA = ctx->VSORTIdxN = 0;
    }
    if (ctx->VSORTPtrA > 0) {
        free((char *)ctx->VSORTPtr);
        memrq(ctx, -ctx->VSORTPtrA,sizeof(int));
        ctx->VSORTPtrA = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  dvcomp()   Compare function for vsort                                   */

int dvcomp(const void *arg1, const void *arg2, void *_ctx)
{
    TDAContext *ctx = (TDAContext *)_ctx;     
    register int i,j;
    double tmp;

    for (i = 0; i < ctx->VSORTIdxN; ++i) {
        j = ctx->VSORTIdx[i];
        tmp = (double)((float)(get_data(ctx, j,*(int *)arg1) - get_data(ctx, j,*(int *)arg2)));
        if (tmp > 0.0)
            return(1);
        else if (tmp < 0.0)
            return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  sorti(n,x,opt)  sort x[] containing n values in ascending order.        */
/*                  if opt=0 ascending order, else descending.              */
/*                  Return 0 if OK, otherwise -1 (insuff memory).           */

int sorti(TDAContext *ctx, int n,int *x,int opt)
{
    register int i,j;
    int err,ptr_a,srt_a;
    int *ptr;

    if (n < 2)
        return(0);
    ptr_a = srt_a = 0;
    err = -1;

    /* allocate a pointer for sorting */
     
    if (!(ptr = (int *)calloc((size_t)(n + 3),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto SRTIFin;
    }
    memrq(ctx, n,sizeof(int));
    ptr_a = 1;

    if (!(ctx->SORTVI = (int *)calloc((size_t)(n),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto SRTIFin;
    }
    memrq(ctx, n,sizeof(int));
    srt_a = 1;
    
    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        ctx->SORTVI[i] = x[i];
    }
            
    tda_qsort_r((char *)ptr,(size_t)(n),sizeof(int), sicomp, ctx);

    if (opt == 0) {
        for (i = 0; i < n; ++i)       
            x[i] = ctx->SORTVI[ptr[i]];
    }
    else {
        j = n;  
        for (i = 0; i < n; ++i)  
            x[--j] = ctx->SORTVI[ptr[i]];
    }
    err = 0;

SRTIFin:
    if (srt_a) {
        free((char *)ctx->SORTVI);
        memrq(ctx, -n,sizeof(int));
    }
    if (ptr_a) {
        free((char *)ptr);
        memrq(ctx, -n,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sorti2(n,x,y)       sort x[] containing n values in ascending order.    */
/*                      change values in y[] in the same order.             */
/*                      Return 0 if OK, otherwise -1 (insuff memory).       */

int sorti2(TDAContext *ctx, int n,int *x,int *y)
{
    register int i;
    int err,ptr_a,srt_a;
    int *ptr;

    if (n < 2)
        return(0);
    ptr_a = srt_a = 0;
    err = -1;

    /* allocate a pointer for sorting */

    if (!(ptr = (int *)calloc((size_t)(n),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto SRTI2Fin;
    }
    memrq(ctx, n,sizeof(int));
    ptr_a = 1;

    if (!(ctx->SORTVI = (int *)calloc((size_t)(n),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto SRTI2Fin;
    }
    memrq(ctx, n,sizeof(int));
    srt_a = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        ctx->SORTVI[i] = x[i];
    }

    tda_qsort_r((char *)ptr,(size_t)(n),sizeof(int), sicomp, ctx);

    for (i = 0; i < n; ++i)  
        x[i] = ctx->SORTVI[ptr[i]];
    for (i = 0; i < n; ++i)  
        ctx->SORTVI[i] = y[ptr[i]];
    for (i = 0; i < n; ++i)  
        y[i] = ctx->SORTVI[i];               

    err = 0;

SRTI2Fin:
    if (srt_a) {
        free((char *)ctx->SORTVI);
        memrq(ctx, -n,sizeof(int));
    }
    if (ptr_a) {
        free((char *)ptr);
        memrq(ctx, -n,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sicomp()    compare function                                            */

int sicomp(const void *arg1, const void *arg2, void *_ctx)
{
    TDAContext *ctx = (TDAContext *)_ctx;     
    int x; 

    x = ctx->SORTVI[*(int *)arg1] - ctx->SORTVI[*(int *)arg2];
    if (x > 0)
        return(1);
    else if (x < 0)
        return(-1);
    return(0);
}
   
/* ------------------------------------------------------------------------ */
/*  sortdpi(n,x,ptr) create pointer ptr for sorting x[] in ascending order. */
/*                   Return 0 if OK, otherwise -1 (insuff memory).          */

int sortdpi(TDAContext *ctx, int n,int *x,int *ptr)
{
    register int i;
    int err,srt_a;

    if (n < 2) {
        ptr[0] = 0;
        return(0);
    }
    srt_a = 0;
    err = -1;

    if (!(ctx->SORTVI = (int *)calloc((size_t)(n),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto SRTDPIFin;
    }
    memrq(ctx, n,sizeof(int));
    srt_a = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        ctx->SORTVI[i] = x[i];
    }

    tda_qsort_r((char *)ptr,(size_t)(n),sizeof(int), sicomp, ctx);
    err = 0;

SRTDPIFin:
    if (srt_a) {
        free((char *)ctx->SORTVI);
        memrq(ctx, -n,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sortdpi2(n,x,y,ptr)     create pointer ptr for sorting x[] and y[]      */
/*                          ascending order. First wrt x, then wrt y.       */
/*                          Return 0 if OK, otherwise -1 (insuff memory).   */

int sortdpi2(TDAContext *ctx, int n,int *x,int *y,int *ptr)
{
    register int i;

    if (n < 2) {
        ptr[0] = 0;
        return(0);
    }
    ctx->SORTVI = x;
    ctx->SORTVI1 = y;
    for (i = 0; i < n; ++i)  
        ptr[i] = i;

    tda_qsort_r((char *)ptr,(size_t)(n),sizeof(int), si2comp, ctx);

    return(0);
}

/* ------------------------------------------------------------------------ */
/*  si2comp()    compare function                                           */

int si2comp(const void *arg1, const void *arg2, void *_ctx)
{
    TDAContext *ctx = (TDAContext *)_ctx;     
    int x; 

    x = ctx->SORTVI[*(int *)arg1] - ctx->SORTVI[*(int *)arg2];
    if (x > 0)
        return(1);
    else if (x < 0)
        return(-1);

    x = ctx->SORTVI1[*(int *)arg1] - ctx->SORTVI1[*(int *)arg2];
    if (x > 0)
        return(1);
    else if (x < 0)
        return(-1);
    return(0);
}
   
/* ------------------------------------------------------------------------ */
/*  sortdp(n,x,ptr) create pointer ptr for sorting x[] in ascending order.  */
/*                  Return 0 if OK, otherwise -1 (insuff memory).           */

int sortdp(TDAContext *ctx, int n,double *x,int *ptr)
{
    register int i;

    if (n < 2) {
        ptr[0] = 0;
        return(0);
    }
    for (i = 0; i < n; ++i)  
        ptr[i] = i;
       
    ctx->SORTV = x;
    tda_qsort_r((char *)ptr,(size_t)(n),sizeof(int), s1comp, ctx);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  sortdp2(n,x,y,ptr)  create pointer ptr for sorting x[], y[] in          */
/*                      ascending order. First wrt x, then wrt y.           */
/*                      Return 0 if OK, otherwise -1 (insuff memory).       */

int sortdp2(TDAContext *ctx, int n,double *x,double *y,int *ptr)
{
    register int i;
    int err,srt_a,srt_a1;

    if (n < 2) {
        ptr[0] = 0;
        return(0);
    }
    srt_a1 = srt_a = 0;
    err = -1;

    if (!(ctx->SORTV = (double *)calloc((size_t)(n),sizeof(double)))) { 
        p_err(ctx, -2,1);
        goto SRTDFin;
    }
    memrq(ctx, n,sizeof(double));
    srt_a = 1;

    if (!(ctx->SORTV1 = (double *)calloc((size_t)(n),sizeof(double)))) { 
        p_err(ctx, -2,1);
        goto SRTDFin;
    }
    memrq(ctx, n,sizeof(double));
    srt_a1 = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        ctx->SORTV[i] = x[i];
        ctx->SORTV1[i] = y[i];
    }

    tda_qsort_r((char *)ptr,(size_t)(n),sizeof(int), s2comp, ctx);
    err = 0;

SRTDFin:
    if (srt_a) {
        free((char *)ctx->SORTV);
        memrq(ctx, -n,sizeof(double));
    }
    if (srt_a1) {
        free((char *)ctx->SORTV1);
        memrq(ctx, -n,sizeof(double));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  s2comp()    compare function                                            */

int s2comp(const void *arg1, const void *arg2, void *_ctx)
{
    TDAContext *ctx = (TDAContext *)_ctx;     
    double x; 

    x = ctx->SORTV[*(int *)arg1] - ctx->SORTV[*(int *)arg2];
    if (x > 0.0)
        return(1);
    else if (x < 0.0)
        return(-1);
    else {
        x = ctx->SORTV1[*(int *)arg1] - ctx->SORTV1[*(int *)arg2];
        if (x > 0.0)
            return(1);
        else if (x < 0.0)
            return(-1);
    }
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  sortdp2f(n,x,y,ptr) create pointer ptr for sorting x[], y[] in          */
/*                      ascending order. First wrt x, then wrt y.           */
/*                      Return 0 if OK, otherwise -1 (insuff memory).       */

int sortdp2f(TDAContext *ctx, int n,float *x,float *y,int *ptr)
{
    register int i;
    int err,srt_a,srt_a1;

    if (n < 2) {
        ptr[0] = 0;
        return(0);
    }
    srt_a1 = srt_a = 0;
    err = -1;

    if (!(ctx->SORTVF = (float *)calloc((size_t)(n),sizeof(float)))) { 
        p_err(ctx, -2,1);
        goto SRTDFFin;
    }
    memrq(ctx, n,sizeof(float));
    srt_a = 1;

    if (!(ctx->SORTVF1 = (float *)calloc((size_t)(n),sizeof(float)))) { 
        p_err(ctx, -2,1);
        goto SRTDFFin;
    }
    memrq(ctx, n,sizeof(float));
    srt_a1 = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        ctx->SORTVF[i] = x[i];
        ctx->SORTVF1[i] = y[i];
    }
    tda_qsort_r((char *)ptr,(size_t)(n),sizeof(int), s2fcomp, ctx);
    err = 0;

SRTDFFin:
    if (srt_a) {
        free((char *)ctx->SORTVF);
        memrq(ctx, -n,sizeof(float));
    }
    if (srt_a1) {
        free((char *)ctx->SORTVF1);
        memrq(ctx, -n,sizeof(float));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  s2fcomp()   compare function                                            */

int s2fcomp(const void *arg1, const void *arg2, void *_ctx)
{
    TDAContext *ctx = (TDAContext *)_ctx;     
    double x; 
    
    x = (double)(ctx->SORTVF[*(int *)arg1] - ctx->SORTVF[*(int *)arg2]);
    if (x > 0.0)
        return(1);
    else if (x < 0.0)
        return(-1);
    else {
        x = (double)(ctx->SORTVF1[*(int *)arg1] - ctx->SORTVF1[*(int *)arg2]);
        if (x > 0.0)
            return(1);
        else if (x < 0.0)
            return(-1);
    }
    return(0);
}

/* -###-------------------------------------------------------------------- */
/*  sortdp4f(n,xa,ya,xb,yb,ptr)                                             */
/*                                                                          */
/*  Create pointer ptr for sorting xa[], ya[], xb[] and yb[] in ascending   */
/*  order. Return 0 if OK, otherwise -1 (insuff memory).                    */

int sortdp4f(TDAContext *ctx, int n,float *xa,float *ya,float *xb,float *yb,int *ptr)
{
    register int i;
    int err,srt_a,srt_a1,srt_a2,srt_a3;

    if (n < 2) {
        ptr[0] = 0;
        return(0);
    }
    srt_a3 = srt_a2 = srt_a1 = srt_a = 0;
    err = -1;

    if (!(ctx->SORTVF = (float *)calloc((size_t)(n),sizeof(float)))) { 
        p_err(ctx, -2,1);
        goto SRTDFFin;
    }
    memrq(ctx, n,sizeof(float));
    srt_a = 1;

    if (!(ctx->SORTVF1 = (float *)calloc((size_t)(n),sizeof(float)))) { 
        p_err(ctx, -2,1);
        goto SRTDFFin;
    }
    memrq(ctx, n,sizeof(float));
    srt_a1 = 1;

    if (!(ctx->SORTVF2 = (float *)calloc((size_t)(n),sizeof(float)))) { 
        p_err(ctx, -2,1);
        goto SRTDFFin;
    }
    memrq(ctx, n,sizeof(float));
    srt_a2 = 1;

    if (!(ctx->SORTVF3 = (float *)calloc((size_t)(n),sizeof(float)))) { 
        p_err(ctx, -2,1);
        goto SRTDFFin;
    }
    memrq(ctx, n,sizeof(float));
    srt_a3 = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        ctx->SORTVF[i] = xa[i];
        ctx->SORTVF1[i] = ya[i];
        ctx->SORTVF2[i] = xb[i];
        ctx->SORTVF3[i] = yb[i];
    }
    tda_qsort_r((char *)ptr,(size_t)(n),sizeof(int), s4fcomp, ctx);
    err = 0;

SRTDFFin:
    if (srt_a) {
        free((char *)ctx->SORTVF);
        memrq(ctx, -n,sizeof(float));
    }
    if (srt_a1) {
        free((char *)ctx->SORTVF1);
        memrq(ctx, -n,sizeof(float));
    }
    if (srt_a2) {
        free((char *)ctx->SORTVF2);
        memrq(ctx, -n,sizeof(float));
    }
    if (srt_a3) {
        free((char *)ctx->SORTVF3);
        memrq(ctx, -n,sizeof(float));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  s4fcomp()   compare function                                            */

int s4fcomp(const void *arg1, const void *arg2, void *_ctx)
{
    TDAContext *ctx = (TDAContext *)_ctx;     
    double x; 
    
    x = (double)(ctx->SORTVF[*(int *)arg1] - ctx->SORTVF[*(int *)arg2]);
    if (x > 0.0)
        return(1);
    else if (x < 0.0)
        return(-1);
    else {
        x = (double)(ctx->SORTVF1[*(int *)arg1] - ctx->SORTVF1[*(int *)arg2]);
        if (x > 0.0)
            return(1);
        else if (x < 0.0)
            return(-1);
        else {
            x = (double)(ctx->SORTVF2[*(int *)arg1] - ctx->SORTVF2[*(int *)arg2]);
            if (x > 0.0)
                return(1);
            else if (x < 0.0)
                return(-1);
            else {
                x = (double)(ctx->SORTVF3[*(int *)arg1] - ctx->SORTVF3[*(int *)arg2]);
                if (x > 0.0)
                    return(1);
                else if (x < 0.0)
                    return(-1);
            }
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  sortdp2a(n,x,c,ptr,opt)                                                 */
/*                                                                          */
/*                       create pointer ptr for sorting x[], c[] in         */
/*                       ascending order. First wrt x in ascending order,   */
/*                       then wrt c in descending order (if opt == 0) or    */
/*                       in ascending order (if opt != 0).                  */
/*                       Return 0 if OK, otherwise -1 (insuff memory).      */

int sortdp2a(TDAContext *ctx, int n,double *x,short *c,int *ptr,int opt)
{
    register int i;
    int err,srt_a,srt_a1;

    if (n < 2) {
        ptr[0] = 0;
        return(0);
    }
    srt_a1 = srt_a = 0;
    err = -1;

    if (!(ctx->SORTV = (double *)calloc((size_t)(n),sizeof(double)))) { 
        p_err(ctx, -2,1);
        goto SRTD2AFin;
    }
    memrq(ctx, n,sizeof(double));
    srt_a = 1;

    if (!(ctx->SORTVS = (short *)calloc((size_t)(n),sizeof(short)))) { 
        p_err(ctx, -2,1);
        goto SRTD2AFin;
    }
    memrq(ctx, n,sizeof(short));
    srt_a1 = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        ctx->SORTV[i] = x[i];
        ctx->SORTVS[i] = c[i];
    }
    if (opt == 0)
        tda_qsort_r((char *)ptr,(size_t)(n),sizeof(int), s2acomp, ctx);
    else
        tda_qsort_r((char *)ptr,(size_t)(n),sizeof(int), s2bcomp, ctx);
    err = 0;

SRTD2AFin:
    if (srt_a) {
        free((char *)ctx->SORTV);
        memrq(ctx, -n,sizeof(double));
    }
    if (srt_a1) {
        free((char *)ctx->SORTVS);
        memrq(ctx, -n,sizeof(short));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  s2acomp()   compare function                                            */

int s2acomp(const void *arg1, const void *arg2, void *_ctx)
{
    TDAContext *ctx = (TDAContext *)_ctx;     
    double x; 
    short c;

    x = ctx->SORTV[*(int *)arg1] - ctx->SORTV[*(int *)arg2];
    if (x > ctx->EPSI1)
        return(1);
    else if (x < -ctx->EPSI1)
        return(-1);
    else {
        c = ctx->SORTVS[*(int *)arg1] - ctx->SORTVS[*(int *)arg2];
        if (c > 0)
            return(-1);
        else if (c < 0)
            return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  s2bcomp()   compare function                                            */

int s2bcomp(const void *arg1, const void *arg2, void *_ctx)
{
    TDAContext *ctx = (TDAContext *)_ctx;     
    double x; 
    short c;

    x = ctx->SORTV[*(int *)arg1] - ctx->SORTV[*(int *)arg2];
    if (x > ctx->EPSI1)
        return(1);
    else if (x < -ctx->EPSI1)
        return(-1);
    else {
        c = ctx->SORTVS[*(int *)arg1] - ctx->SORTVS[*(int *)arg2];
        if (c > 0)
            return(1);
        else if (c < 0)
            return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  sortdpn(m,n,x,nc,col,ptr)                                               */
/*                                                                          */
/*  x is m x n matrix. col is array with nc elements: col[i], i=1,nc        */  
/*  referring to column numbers of x. Create a pointer ptr that sort x      */
/*  in ascending order, first with respect to col[1], then col[2], etc.     */
/*  ptr[j] j = 1,...,m                                                      */
/*                                                                          */
/*  Return 0 if OK, otherwise -1                                            */

int sortdpn(TDAContext *ctx, int m,int n,double *x,int nc,int *col,int *ptr)
{
    register int i;

    if (nc < 1)
        return(0);

    if (m < 2) {
        ptr[1] = 1;
        return(0);
    }
    ctx->SORTV = x;
    ctx->SORTVN = n;
    ctx->SORTVI = col;
    ctx->SORTVIN = nc;

    for (i = 1; i <= m; ++i)  
        ptr[i] = i;

    tda_qsort_r((char *)(ptr + 1),(size_t)(m),sizeof(int), sncomp, ctx);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  sncomp()    compare function                                            */

int sncomp(const void *arg1, const void *arg2, void *_ctx)
{
    TDAContext *ctx = (TDAContext *)_ctx;     
    register int i,j;
    double x; 

    for (i = 1; i <= ctx->SORTVIN; ++i) {
        j = ctx->SORTVI[i];
        x = ctx->SORTV[(*(int *)arg1 - 1) * ctx->SORTVN + j] -
            ctx->SORTV[(*(int *)arg2 - 1) * ctx->SORTVN + j];

        if (x > 0.0)
            return(1);
        else if (x < 0.0)
            return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  sortdpn1(m,n,x,nc,col,ptr,cen)                                          */
/*                                                                          */
/*  x is m x n matrix. col is array with nc elements: col[i], i=1,nc        */  
/*  referring to column numbers of x. Create a pointer ptr that sort x      */
/*  in ascending order, first with respect to col[1], then col[2], etc.     */
/*  In case of ties, sort in ascending order according to cen.              */
/*  ptr[j] j = 1,...,m                                                      */
/*                                                                          */
/*  Return 0 if OK, otherwise -1                                            */

int sortdpn1(TDAContext *ctx, int m,int n,double *x,int nc,int *col,int *ptr,short *cen)
{
    register int i;

    if (nc < 1)
        return(0);

    if (m < 2) {
        ptr[1] = 1;
        return(0);
    }
    ctx->SORTV = x;
    ctx->SORTVN = n;
    ctx->SORTVI = col;
    ctx->SORTVIN = nc;
    ctx->SORTVS = cen;

    for (i = 1; i <= m; ++i)  
        ptr[i] = i;

    tda_qsort_r((char *)(ptr + 1),(size_t)(m),sizeof(int), sn1comp, ctx);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  sn1comp()    compare function                                           */

int sn1comp(const void *arg1, const void *arg2, void *_ctx)
{
    TDAContext *ctx = (TDAContext *)_ctx;     
    register int i,j;
    short c;
    double x; 

    for (i = 1; i <= ctx->SORTVIN; ++i) {
        j = ctx->SORTVI[i];
        x = ctx->SORTV[(*(int *)arg1 - 1) * ctx->SORTVN + j] -
            ctx->SORTV[(*(int *)arg2 - 1) * ctx->SORTVN + j];

        if (x > 0.0)
            return(1);
        else if (x < 0.0)
            return(-1);
    }
    c = ctx->SORTVS[*(int *)arg1] - ctx->SORTVS[*(int *)arg2];
    if (c > 0)
        return(1);
    else if (c < 0)
        return(-1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  sortdp2c(n,x,y,ptr,xmin,ymin,xmax,ymax)                                 */
/*                                                                          */
/*  create pointer ptr for sorting x and y in cyclical order around the     */
/*  borders of the rectangle defined by xmin,ymin,xmax,ymax.                */
/*                                                                          */
/*  Return 0 if OK, otherwise -1 (insuff memory).                           */

int sortdp2c(TDAContext *ctx, int n,double *x,double *y,int *ptr,double xmin,double ymin, double xmax,double ymax)
{
    register int i;
    int err,srt_a,srt_a1;

    if (n < 2) {
        ptr[0] = 0;
        return(0);
    }
    srt_a1 = srt_a = 0;
    err = -1;

    if (!(ctx->SORTV = (double *)calloc((size_t)(n),sizeof(double)))) { 
        p_err(ctx, -2,1);
        goto SRTDCFin;
    }
    memrq(ctx, n,sizeof(double));
    srt_a = 1;

    if (!(ctx->SORTV1 = (double *)calloc((size_t)(n),sizeof(double)))) { 
        p_err(ctx, -2,1);
        goto SRTDCFin;
    }
    memrq(ctx, n,sizeof(double));
    srt_a1 = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        ctx->SORTV[i] = x[i];
        ctx->SORTV1[i] = y[i];
    }
    ctx->SORTXMin = xmin;
    ctx->SORTXMax = xmax;
    ctx->SORTYMin = ymin;
    ctx->SORTYMax = ymax;

    tda_qsort_r((char *)ptr,(size_t)(n),sizeof(int), s2ccomp, ctx);
    err = 0;

SRTDCFin:
    if (srt_a) {
        free((char *)ctx->SORTV);
        memrq(ctx, -n,sizeof(double));
    }
    if (srt_a1) {
        free((char *)ctx->SORTV1);
        memrq(ctx, -n,sizeof(double));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  s2ccomp()   compare function                                            */

int s2ccomp(const void *arg1, const void *arg2, void *_ctx)
{
    TDAContext *ctx = (TDAContext *)_ctx;     
    int ra,rb;
    double xa,ya,xb,yb;

    xa = ctx->SORTV[*(int *)arg1];
    xb = ctx->SORTV[*(int *)arg2];
    ya = ctx->SORTV1[*(int *)arg1];
    yb = ctx->SORTV1[*(int *)arg2];

    ra = g_left(ctx, ctx->SORTXMin,ctx->SORTYMin,ctx->SORTXMax,ctx->SORTYMax,xa,ya);
    rb = g_left(ctx, ctx->SORTXMin,ctx->SORTYMin,ctx->SORTXMax,ctx->SORTYMax,xb,yb);
    /**  
    tda_out("xa=%lg,%lg xb=%lg %lg ra=%d rb=%d\n",xa,ya,xb,yb,ra,rb);
    **/    

    if (ra == 0) {
        if (rb ==  1)   
            return(-1); 
        else {
            if (xa < xb - ctx->EPSI1) 
                return(-1);       
            else if (xa > xb + ctx->EPSI1)
                return(1);           
            else if (ya < yb)        
                return(-1);          
            else if (ya > yb)        
                return(1);           
            return(0);
        }
    }
    else { 
        if (rb == 0)   
            return(1); 
        else {
            if (xa < xb - ctx->EPSI1) 
                return(1);       
            else if (xa > xb + ctx->EPSI1)
                return(-1);          
            else if (ya < yb)        
                return(1);           
            else if (ya > yb)        
                return(-1);          
            return(0);
        }
    }
}



