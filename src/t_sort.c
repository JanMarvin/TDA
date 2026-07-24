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

/*  functions in t_sort.c */

int sortd(int n,double *x,int opt);
int sortd2(int n,double *x,double *y);
int sortd3(int n,double *x,double *y,double *z);
int s1comp(const void *arg1,const void *arg2);
int vsort(int n,short *vidx,int opt,int cflg,int pflg);   
int dvcomp(const void *arg1,const void *arg2);
int sorti(int n,int *x,int opt);
int sorti2(int n,int *x,int *y);
int sicomp(const void *arg1,const void *arg2);
int sortdpi(int n,int *x,int *ptr);
int sortdpi2(int n,int *x,int *y,int *ptr);
int si2comp(const void *arg1,const void *arg2);
int sortdp(int n,double *x,int *ptr);
int sortdp2(int n,double *x,double *y,int *ptr);
int s2comp(const void *arg1,const void *arg2);
int sortdp2f(int n,float *x,float *y,int *ptr);
int s2fcomp(const void *arg1,const void *arg2);
int sortdp4f(int n,float *xa,float *ya,float *xb,float *yb,int *ptr);
int s4fcomp(const void *arg1,const void *arg2);
int sortdp2a(int n,double *x,short *c,int *ptr,int opt);
int s2acomp(const void *arg1,const void *arg2);
int s2bcomp(const void *arg1,const void *arg2);
int sortdpn(int m,int n,double *x,int nc,int *col,int *ptr);
int sncomp(const void *arg1,const void *arg2);
int sortdpn1(int m,int n,double *x,int nc,int *col,int *ptr,short *cen);
int sn1comp(const void *arg1,const void *arg2);
int sortdp2c(int n,double *x,double *y,int *ptr,double xmin,double ymin,
    double xmax,double ymax);
int s2ccomp(const void *arg1,const void *arg2);

double *SORTV;
double *SORTV1;
float *SORTVF;
float *SORTVF1;
float *SORTVF2;
float *SORTVF3;
short *SORTVS;
int SORTVN = 0;
int *SORTVI;
int *SORTVI1;
int SORTVIN = 0;
int *VSORTPtr;      /* pointer to sorted data matrix cases                  */
int VSORTPtrA = 0;  /* if allocated                                         */
int *VSORTIdx;      /* variables for sorting                                */
int VSORTIdxA = 0;  /* if allocated                                         */
int VSORTIdxN = 0;  /* number of variables in VSORTIdx[]                    */
double SORTXMin = 0.0;
double SORTXMax = 0.0;
double SORTYMin = 0.0;
double SORTYMax = 0.0;

/* ------------------------------------------------------------------------ */
/*  sortd(n,x,opt)  sort x[] containing n values in ascending order.        */
/*                  if opt=0 ascending order, else descending.              */
/*                  Return 0 if OK, otherwise -1 (insuff memory).           */

int sortd(int n,double *x,int opt)
{
    register int i,j;
    int err,ptr_a,srt_a;
    int *ptr;

    if (n < 2)
        return(0);
    ptr_a = srt_a = 0;
    err = -1;

    /* allocate a pointer for sorting */

    if (!(ptr = (int *)calloc(n,sizeof(int)))) { 
        p_err(-2,1);
        goto SRTDFin;
    }
    memrq(n,sizeof(int));
    ptr_a = 1;

    if (!(SORTV = (double *)calloc(n,sizeof(double)))) { 
        p_err(-2,1);
        goto SRTDFin;
    }
    memrq(n,sizeof(double));
    srt_a = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        SORTV[i] = x[i];
    }

    qsort((char *)ptr,n,sizeof(int),s1comp);

    if (opt == 0) {
        for (i = 0; i < n; ++i)  
            x[i] = SORTV[ptr[i]];
    }
    else {
        j = n;  
        for (i = 0; i < n; ++i)  
            x[--j] = SORTV[ptr[i]];
    }
    err = 0;

SRTDFin:
    if (srt_a) {
        free((char *)SORTV);
        memrq(-n,sizeof(double));
    }
    if (ptr_a) {
        free((char *)ptr);
        memrq(-n,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sortd2(n,x,y)       sort x[] containing n values in ascending order.    */
/*                      change values in y[] in the same order.             */
/*                      Return 0 if OK, otherwise -1 (insuff memory).       */

int sortd2(int n,double *x,double *y)
{
    register int i;
    int err,ptr_a,srt_a;
    int *ptr;

    if (n < 2)
        return(0);
    ptr_a = srt_a = 0;
    err = -1;

    /* allocate a pointer for sorting */

    if (!(ptr = (int *)calloc(n,sizeof(int)))) { 
        p_err(-2,1);
        goto SRTD2Fin;
    }
    memrq(n,sizeof(int));
    ptr_a = 1;

    if (!(SORTV = (double *)calloc(n,sizeof(double)))) { 
        p_err(-2,1);
        goto SRTD2Fin;
    }
    memrq(n,sizeof(double));
    srt_a = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        SORTV[i] = x[i];
    }

    qsort((char *)ptr,n,sizeof(int),s1comp);

    for (i = 0; i < n; ++i)  
        x[i] = SORTV[ptr[i]];
    for (i = 0; i < n; ++i)  
        SORTV[i] = y[ptr[i]];
    for (i = 0; i < n; ++i)  
        y[i] = SORTV[i];               

    err = 0;

SRTD2Fin:
    if (srt_a) {
        free((char *)SORTV);
        memrq(-n,sizeof(double));
    }
    if (ptr_a) {
        free((char *)ptr);
        memrq(-n,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sortd3(n,x,y,z)     sort x[] containing n values in ascending order.    */
/*                      change values in y[] and z[] in the same order.     */
/*                      Return 0 if OK, otherwise -1 (insuff memory).       */

int sortd3(int n,double *x,double *y,double *z)
{
    register int i;
    int err,ptr_a,srt_a;
    int *ptr;

    if (n < 2)
        return(0);
    ptr_a = srt_a = 0;
    err = -1;

    /* allocate a pointer for sorting */

    if (!(ptr = (int *)calloc(n,sizeof(int)))) { 
        p_err(-2,1);
        goto SRTD3Fin;
    }
    memrq(n,sizeof(int));
    ptr_a = 1;

    if (!(SORTV = (double *)calloc(n,sizeof(double)))) { 
        p_err(-2,1);
        goto SRTD3Fin;
    }
    memrq(n,sizeof(double));
    srt_a = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        SORTV[i] = x[i];
    }

    qsort((char *)ptr,n,sizeof(int),s1comp);

    for (i = 0; i < n; ++i)  
        x[i] = SORTV[ptr[i]];
    for (i = 0; i < n; ++i)  
        SORTV[i] = y[ptr[i]];
    for (i = 0; i < n; ++i)  
        y[i] = SORTV[i];               
    for (i = 0; i < n; ++i)  
        SORTV[i] = z[ptr[i]];
    for (i = 0; i < n; ++i)  
        z[i] = SORTV[i];               

    err = 0;

SRTD3Fin:
    if (srt_a) {
        free((char *)SORTV);
        memrq(-n,sizeof(double));
    }
    if (ptr_a) {
        free((char *)ptr);
        memrq(-n,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  s1comp()    compare function                                            */

int s1comp(const void *arg1,const void *arg2)
{     
    double x; 

    x = SORTV[*(int *)arg1] - SORTV[*(int *)arg2];
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

int vsort(int n,short *vidx,int opt,int cflg,int pflg)
{     
    register int i,j,k;
    int err,r;
    float tmp;

    err = -1;
    if (opt == 0) {
        err = 0;
        goto VSORTFin;
    }
    if (pflg) {
        printf1("Sort with variable(s): %s",VName[vidx[0]]);
        for (i = 1; i < n; ++i)  
            printf1(",%s",VName[vidx[i]]);
        printf1("\n");
        prn_mem();
    }
    r = 0;
    for (i = 0; i < n; ++i) {
        if (VTyp[vidx[i]] == 1) {
            r = 1;
            break;
        }
    }
    if (r) {
        printf1("Error: can't sort string variables.\n");
        goto VSORTFin;
    }

    if (!(VSORTIdx = (int *)calloc(n,sizeof(int)))) { 
        p_err(-2,1);
        goto VSORTFin;
    }
    memrq(n,sizeof(int));
    VSORTIdxA = VSORTIdxN = n;
    for (i = 0; i < n; ++i)
        VSORTIdx[i] = vidx[i];

    if (!(VSORTPtr = (int *)calloc(NOC,sizeof(int)))) { 
        p_err(-2,1);
        goto VSORTFin;
    }
    memrq(NOC,sizeof(int));
    VSORTPtrA = NOC;

    for (i = 0; i < NOC; ++i)
        VSORTPtr[i] = i;

    if (cflg) {
        for (i = 1; i < NOC; ++i) {

            for (k = 0; k < VSORTIdxN; ++k) {
                j = VSORTIdx[k];
                tmp = (float)(get_data(j,i) - get_data(j,i - 1));
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
                printf1("Data already sorted.\n");
            return(0);
        }
        if (pflg)
            printf1("Data not already sorted, see records %d and %d.\n",i,i + 1);
    }
    qsort((char *)VSORTPtr,NOC,sizeof(int),dvcomp);
    return(0);

VSORTFin:
    if (VSORTIdxA > 0) {
        free((char *)VSORTIdx);
        memrq(-VSORTIdxA,sizeof(int));
        VSORTIdxA = VSORTIdxN = 0;
    }
    if (VSORTPtrA > 0) {
        free((char *)VSORTPtr);
        memrq(-VSORTPtrA,sizeof(int));
        VSORTPtrA = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  dvcomp()   Compare function for vsort                                   */

int dvcomp(const void *arg1,const void *arg2)
{     
    register int i,j;
    float tmp;

    for (i = 0; i < VSORTIdxN; ++i) {
        j = VSORTIdx[i];
        tmp = (float)(get_data(j,*(int *)arg1) - get_data(j,*(int *)arg2));
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

int sorti(int n,int *x,int opt)
{
    register int i,j;
    int err,ptr_a,srt_a;
    int *ptr;

    if (n < 2)
        return(0);
    ptr_a = srt_a = 0;
    err = -1;

    /* allocate a pointer for sorting */
     
    if (!(ptr = (int *)calloc(n + 3 ,sizeof(int)))) { 
        p_err(-2,1);
        goto SRTIFin;
    }
    memrq(n,sizeof(int));
    ptr_a = 1;

    if (!(SORTVI = (int *)calloc(n,sizeof(int)))) { 
        p_err(-2,1);
        goto SRTIFin;
    }
    memrq(n,sizeof(int));
    srt_a = 1;
    
    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        SORTVI[i] = x[i];
    }
            
    qsort((char *)ptr,n,sizeof(int),sicomp);

    if (opt == 0) {
        for (i = 0; i < n; ++i)       
            x[i] = SORTVI[ptr[i]];
    }
    else {
        j = n;  
        for (i = 0; i < n; ++i)  
            x[--j] = SORTVI[ptr[i]];
    }
    err = 0;

SRTIFin:
    if (srt_a) {
        free((char *)SORTVI);
        memrq(-n,sizeof(int));
    }
    if (ptr_a) {
        free((char *)ptr);
        memrq(-n,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sorti2(n,x,y)       sort x[] containing n values in ascending order.    */
/*                      change values in y[] in the same order.             */
/*                      Return 0 if OK, otherwise -1 (insuff memory).       */

int sorti2(int n,int *x,int *y)
{
    register int i;
    int err,ptr_a,srt_a;
    int *ptr;

    if (n < 2)
        return(0);
    ptr_a = srt_a = 0;
    err = -1;

    /* allocate a pointer for sorting */

    if (!(ptr = (int *)calloc(n,sizeof(int)))) { 
        p_err(-2,1);
        goto SRTI2Fin;
    }
    memrq(n,sizeof(int));
    ptr_a = 1;

    if (!(SORTVI = (int *)calloc(n,sizeof(int)))) { 
        p_err(-2,1);
        goto SRTI2Fin;
    }
    memrq(n,sizeof(int));
    srt_a = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        SORTVI[i] = x[i];
    }

    qsort((char *)ptr,n,sizeof(int),sicomp);

    for (i = 0; i < n; ++i)  
        x[i] = SORTVI[ptr[i]];
    for (i = 0; i < n; ++i)  
        SORTVI[i] = y[ptr[i]];
    for (i = 0; i < n; ++i)  
        y[i] = SORTVI[i];               

    err = 0;

SRTI2Fin:
    if (srt_a) {
        free((char *)SORTVI);
        memrq(-n,sizeof(int));
    }
    if (ptr_a) {
        free((char *)ptr);
        memrq(-n,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sicomp()    compare function                                            */

int sicomp(const void *arg1,const void *arg2)
{     
    int x; 

    x = SORTVI[*(int *)arg1] - SORTVI[*(int *)arg2];
    if (x > 0)
        return(1);
    else if (x < 0)
        return(-1);
    return(0);
}
   
/* ------------------------------------------------------------------------ */
/*  sortdpi(n,x,ptr) create pointer ptr for sorting x[] in ascending order. */
/*                   Return 0 if OK, otherwise -1 (insuff memory).          */

int sortdpi(int n,int *x,int *ptr)
{
    register int i;
    int err,srt_a;

    if (n < 2) {
        ptr[0] = 0;
        return(0);
    }
    srt_a = 0;
    err = -1;

    if (!(SORTVI = (int *)calloc(n,sizeof(int)))) { 
        p_err(-2,1);
        goto SRTDPIFin;
    }
    memrq(n,sizeof(int));
    srt_a = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        SORTVI[i] = x[i];
    }

    qsort((char *)ptr,n,sizeof(int),sicomp);
    err = 0;

SRTDPIFin:
    if (srt_a) {
        free((char *)SORTVI);
        memrq(-n,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sortdpi2(n,x,y,ptr)     create pointer ptr for sorting x[] and y[]      */
/*                          ascending order. First wrt x, then wrt y.       */
/*                          Return 0 if OK, otherwise -1 (insuff memory).   */

int sortdpi2(int n,int *x,int *y,int *ptr)
{
    register int i;

    if (n < 2) {
        ptr[0] = 0;
        return(0);
    }
    SORTVI = x;
    SORTVI1 = y;
    for (i = 0; i < n; ++i)  
        ptr[i] = i;

    qsort((char *)ptr,n,sizeof(int),si2comp);

    return(0);
}

/* ------------------------------------------------------------------------ */
/*  si2comp()    compare function                                           */

int si2comp(const void *arg1,const void *arg2)
{     
    int x; 

    x = SORTVI[*(int *)arg1] - SORTVI[*(int *)arg2];
    if (x > 0)
        return(1);
    else if (x < 0)
        return(-1);

    x = SORTVI1[*(int *)arg1] - SORTVI1[*(int *)arg2];
    if (x > 0)
        return(1);
    else if (x < 0)
        return(-1);
    return(0);
}
   
/* ------------------------------------------------------------------------ */
/*  sortdp(n,x,ptr) create pointer ptr for sorting x[] in ascending order.  */
/*                  Return 0 if OK, otherwise -1 (insuff memory).           */

int sortdp(int n,double *x,int *ptr)
{
    register int i;

    if (n < 2) {
        ptr[0] = 0;
        return(0);
    }
    for (i = 0; i < n; ++i)  
        ptr[i] = i;
       
    SORTV = x;
    qsort((char *)ptr,n,sizeof(int),s1comp);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  sortdp2(n,x,y,ptr)  create pointer ptr for sorting x[], y[] in          */
/*                      ascending order. First wrt x, then wrt y.           */
/*                      Return 0 if OK, otherwise -1 (insuff memory).       */

int sortdp2(int n,double *x,double *y,int *ptr)
{
    register int i;
    int err,srt_a,srt_a1;

    if (n < 2) {
        ptr[0] = 0;
        return(0);
    }
    srt_a1 = srt_a = 0;
    err = -1;

    if (!(SORTV = (double *)calloc(n,sizeof(double)))) { 
        p_err(-2,1);
        goto SRTDFin;
    }
    memrq(n,sizeof(double));
    srt_a = 1;

    if (!(SORTV1 = (double *)calloc(n,sizeof(double)))) { 
        p_err(-2,1);
        goto SRTDFin;
    }
    memrq(n,sizeof(double));
    srt_a1 = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        SORTV[i] = x[i];
        SORTV1[i] = y[i];
    }

    qsort((char *)ptr,n,sizeof(int),s2comp);
    err = 0;

SRTDFin:
    if (srt_a) {
        free((char *)SORTV);
        memrq(-n,sizeof(double));
    }
    if (srt_a1) {
        free((char *)SORTV1);
        memrq(-n,sizeof(double));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  s2comp()    compare function                                            */

int s2comp(const void *arg1,const void *arg2)
{     
    double x; 

    x = SORTV[*(int *)arg1] - SORTV[*(int *)arg2];
    if (x > 0.0)
        return(1);
    else if (x < 0.0)
        return(-1);
    else {
        x = SORTV1[*(int *)arg1] - SORTV1[*(int *)arg2];
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

int sortdp2f(int n,float *x,float *y,int *ptr)
{
    register int i;
    int err,srt_a,srt_a1;

    if (n < 2) {
        ptr[0] = 0;
        return(0);
    }
    srt_a1 = srt_a = 0;
    err = -1;

    if (!(SORTVF = (float *)calloc(n,sizeof(float)))) { 
        p_err(-2,1);
        goto SRTDFFin;
    }
    memrq(n,sizeof(float));
    srt_a = 1;

    if (!(SORTVF1 = (float *)calloc(n,sizeof(float)))) { 
        p_err(-2,1);
        goto SRTDFFin;
    }
    memrq(n,sizeof(float));
    srt_a1 = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        SORTVF[i] = x[i];
        SORTVF1[i] = y[i];
    }
    qsort((char *)ptr,n,sizeof(int),s2fcomp);
    err = 0;

SRTDFFin:
    if (srt_a) {
        free((char *)SORTVF);
        memrq(-n,sizeof(float));
    }
    if (srt_a1) {
        free((char *)SORTVF1);
        memrq(-n,sizeof(float));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  s2fcomp()   compare function                                            */

int s2fcomp(const void *arg1,const void *arg2)
{     
    float x; 
    
    x = SORTVF[*(int *)arg1] - SORTVF[*(int *)arg2];
    if (x > 0.0)
        return(1);
    else if (x < 0.0)
        return(-1);
    else {
        x = SORTVF1[*(int *)arg1] - SORTVF1[*(int *)arg2];
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

int sortdp4f(int n,float *xa,float *ya,float *xb,float *yb,int *ptr)
{
    register int i;
    int err,srt_a,srt_a1,srt_a2,srt_a3;

    if (n < 2) {
        ptr[0] = 0;
        return(0);
    }
    srt_a3 = srt_a2 = srt_a1 = srt_a = 0;
    err = -1;

    if (!(SORTVF = (float *)calloc(n,sizeof(float)))) { 
        p_err(-2,1);
        goto SRTDFFin;
    }
    memrq(n,sizeof(float));
    srt_a = 1;

    if (!(SORTVF1 = (float *)calloc(n,sizeof(float)))) { 
        p_err(-2,1);
        goto SRTDFFin;
    }
    memrq(n,sizeof(float));
    srt_a1 = 1;

    if (!(SORTVF2 = (float *)calloc(n,sizeof(float)))) { 
        p_err(-2,1);
        goto SRTDFFin;
    }
    memrq(n,sizeof(float));
    srt_a2 = 1;

    if (!(SORTVF3 = (float *)calloc(n,sizeof(float)))) { 
        p_err(-2,1);
        goto SRTDFFin;
    }
    memrq(n,sizeof(float));
    srt_a3 = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        SORTVF[i] = xa[i];
        SORTVF1[i] = ya[i];
        SORTVF2[i] = xb[i];
        SORTVF3[i] = yb[i];
    }
    qsort((char *)ptr,n,sizeof(int),s4fcomp);
    err = 0;

SRTDFFin:
    if (srt_a) {
        free((char *)SORTVF);
        memrq(-n,sizeof(float));
    }
    if (srt_a1) {
        free((char *)SORTVF1);
        memrq(-n,sizeof(float));
    }
    if (srt_a2) {
        free((char *)SORTVF2);
        memrq(-n,sizeof(float));
    }
    if (srt_a3) {
        free((char *)SORTVF3);
        memrq(-n,sizeof(float));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  s4fcomp()   compare function                                            */

int s4fcomp(const void *arg1,const void *arg2)
{     
    float x; 
    
    x = SORTVF[*(int *)arg1] - SORTVF[*(int *)arg2];
    if (x > 0.0)
        return(1);
    else if (x < 0.0)
        return(-1);
    else {
        x = SORTVF1[*(int *)arg1] - SORTVF1[*(int *)arg2];
        if (x > 0.0)
            return(1);
        else if (x < 0.0)
            return(-1);
        else {
            x = SORTVF2[*(int *)arg1] - SORTVF2[*(int *)arg2];
            if (x > 0.0)
                return(1);
            else if (x < 0.0)
                return(-1);
            else {
                x = SORTVF3[*(int *)arg1] - SORTVF3[*(int *)arg2];
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

int sortdp2a(int n,double *x,short *c,int *ptr,int opt)
{
    register int i;
    int err,srt_a,srt_a1;

    if (n < 2) {
        ptr[0] = 0;
        return(0);
    }
    srt_a1 = srt_a = 0;
    err = -1;

    if (!(SORTV = (double *)calloc(n,sizeof(double)))) { 
        p_err(-2,1);
        goto SRTD2AFin;
    }
    memrq(n,sizeof(double));
    srt_a = 1;

    if (!(SORTVS = (short *)calloc(n,sizeof(short)))) { 
        p_err(-2,1);
        goto SRTD2AFin;
    }
    memrq(n,sizeof(short));
    srt_a1 = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        SORTV[i] = x[i];
        SORTVS[i] = c[i];
    }
    if (opt == 0)
        qsort((char *)ptr,n,sizeof(int),s2acomp);
    else
        qsort((char *)ptr,n,sizeof(int),s2bcomp);
    err = 0;

SRTD2AFin:
    if (srt_a) {
        free((char *)SORTV);
        memrq(-n,sizeof(double));
    }
    if (srt_a1) {
        free((char *)SORTVS);
        memrq(-n,sizeof(short));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  s2acomp()   compare function                                            */

int s2acomp(const void *arg1,const void *arg2)
{     
    double x; 
    short c;

    x = SORTV[*(int *)arg1] - SORTV[*(int *)arg2];
    if (x > EPSI1)
        return(1);
    else if (x < -EPSI1)
        return(-1);
    else {
        c = SORTVS[*(int *)arg1] - SORTVS[*(int *)arg2];
        if (c > 0)
            return(-1);
        else if (c < 0)
            return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  s2bcomp()   compare function                                            */

int s2bcomp(const void *arg1,const void *arg2)
{     
    double x; 
    short c;

    x = SORTV[*(int *)arg1] - SORTV[*(int *)arg2];
    if (x > EPSI1)
        return(1);
    else if (x < -EPSI1)
        return(-1);
    else {
        c = SORTVS[*(int *)arg1] - SORTVS[*(int *)arg2];
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

int sortdpn(int m,int n,double *x,int nc,int *col,int *ptr)
{
    register int i;

    if (nc < 1)
        return(0);

    if (m < 2) {
        ptr[1] = 1;
        return(0);
    }
    SORTV = x;
    SORTVN = n;
    SORTVI = col;
    SORTVIN = nc;

    for (i = 1; i <= m; ++i)  
        ptr[i] = i;

    qsort((char *)(ptr + 1),m,sizeof(int),sncomp);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  sncomp()    compare function                                            */

int sncomp(const void *arg1,const void *arg2)
{     
    register int i,j;
    double x; 

    for (i = 1; i <= SORTVIN; ++i) {
        j = SORTVI[i];
        x = SORTV[(*(int *)arg1 - 1) * SORTVN + j] -
            SORTV[(*(int *)arg2 - 1) * SORTVN + j];

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

int sortdpn1(int m,int n,double *x,int nc,int *col,int *ptr,short *cen)
{
    register int i;

    if (nc < 1)
        return(0);

    if (m < 2) {
        ptr[1] = 1;
        return(0);
    }
    SORTV = x;
    SORTVN = n;
    SORTVI = col;
    SORTVIN = nc;
    SORTVS = cen;

    for (i = 1; i <= m; ++i)  
        ptr[i] = i;

    qsort((char *)(ptr + 1),m,sizeof(int),sn1comp);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  sn1comp()    compare function                                           */

int sn1comp(const void *arg1,const void *arg2)
{     
    register int i,j;
    short c;
    double x; 

    for (i = 1; i <= SORTVIN; ++i) {
        j = SORTVI[i];
        x = SORTV[(*(int *)arg1 - 1) * SORTVN + j] -
            SORTV[(*(int *)arg2 - 1) * SORTVN + j];

        if (x > 0.0)
            return(1);
        else if (x < 0.0)
            return(-1);
    }
    c = SORTVS[*(int *)arg1] - SORTVS[*(int *)arg2];
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

int sortdp2c(int n,double *x,double *y,int *ptr,double xmin,double ymin,
    double xmax,double ymax)
{
    register int i;
    int err,srt_a,srt_a1;

    if (n < 2) {
        ptr[0] = 0;
        return(0);
    }
    srt_a1 = srt_a = 0;
    err = -1;

    if (!(SORTV = (double *)calloc(n,sizeof(double)))) { 
        p_err(-2,1);
        goto SRTDCFin;
    }
    memrq(n,sizeof(double));
    srt_a = 1;

    if (!(SORTV1 = (double *)calloc(n,sizeof(double)))) { 
        p_err(-2,1);
        goto SRTDCFin;
    }
    memrq(n,sizeof(double));
    srt_a1 = 1;

    for (i = 0; i < n; ++i) {
        ptr[i] = i;
        SORTV[i] = x[i];
        SORTV1[i] = y[i];
    }
    SORTXMin = xmin;
    SORTXMax = xmax;
    SORTYMin = ymin;
    SORTYMax = ymax;

    qsort((char *)ptr,n,sizeof(int),s2ccomp);
    err = 0;

SRTDCFin:
    if (srt_a) {
        free((char *)SORTV);
        memrq(-n,sizeof(double));
    }
    if (srt_a1) {
        free((char *)SORTV1);
        memrq(-n,sizeof(double));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  s2ccomp()   compare function                                            */

int s2ccomp(const void *arg1,const void *arg2)
{     
    int ra,rb;
    double xa,ya,xb,yb;

    xa = SORTV[*(int *)arg1];
    xb = SORTV[*(int *)arg2];
    ya = SORTV1[*(int *)arg1];
    yb = SORTV1[*(int *)arg2];

    ra = g_left(SORTXMin,SORTYMin,SORTXMax,SORTYMax,xa,ya);
    rb = g_left(SORTXMin,SORTYMin,SORTXMax,SORTYMax,xb,yb);
    /**  
    printf("xa=%lg,%lg xb=%lg %lg ra=%d rb=%d\n",xa,ya,xb,yb,ra,rb);
    **/    

    if (ra == 0) {
        if (rb ==  1)   
            return(-1); 
        else {
            if (xa < xb - EPSI1) 
                return(-1);       
            else if (xa > xb + EPSI1)
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
            if (xa < xb - EPSI1) 
                return(1);       
            else if (xa > xb + EPSI1)
                return(-1);          
            else if (ya < yb)        
                return(1);           
            else if (ya > yb)        
                return(-1);          
            return(0);
        }
    }
}



