/****************************************************************************/
/*  t_alloc                                                                 */
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

/*  functions in t_alloc.c */

void a_clean(void);
int alloc_actmp(int n);
int alloc_actmp1(int n);
int alloc_acx(int n);
int alloc_acy(int n);
int alloc_acy1(int n);
int alloc_acz(int n);
int alloc_act(int n);
int alloc_acu(int n);
int alloc_acv(int n);
int alloc_acw(int n);
int alloc_acxf(int n);
int alloc_acyf(int n);
int alloc_aczf(int n);
int alloc_acuf(int n);
int alloc_acvf(int n);
int alloc_acm(int n);
int alloc_acn(int n);
int alloc_aci(int n);
int alloc_acj(int n);
int alloc_acr(int n);
int alloc_acs(int n);
int alloc_ack(int n);
int alloc_acl(int n);
int alloc_acc(int n);
int alloc_acd(int n);
int alloc_ace(int n);
int alloc_acf(int n);
int alloc_acns(int n);
int alloc_acms(int n);
int alloc_acptr(int n);
int alloc_aciptr(int n);
int alloc_acjptr(int n);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

int AcTmpN = 0;
double *AcTmp;
int AcTmp1N = 0;
double *AcTmp1;
int AcXN = 0;
double *AcX;
int AcYN = 0;
double *AcY;
int AcY1N = 0;
double *AcY1;
int AcZN = 0;
double *AcZ;
int AcTN = 0;
double *AcT;
int AcUN = 0;
double *AcU;
int AcVN = 0;
double *AcV;
int AcWN = 0;
double *AcW;
int AcXFN = 0;
float *AcXF;
int AcYFN = 0;
float *AcYF;
int AcZFN = 0;
float *AcZF;
int AcUFN = 0;
float *AcUF;
int AcVFN = 0;
float *AcVF;
int AcNM = 0;
int *AcM;
int AcNN = 0;
int *AcN;
int AcIN = 0;
int *AcI;
int AcJN = 0;
int *AcJ;
int AcRN = 0;
int *AcR;
int AcSN = 0;
int *AcS;
int AcKN = 0;
int *AcK;
int AcLN = 0;
int *AcL;
int AcCN = 0;
char *AcC;
int AcDN = 0;
char *AcD;
int AcEN = 0;
char *AcE;
int AcFN = 0;
char *AcF;
int AcNNS = 0;
short *AcNS;
int AcNMS = 0;
short *AcMS;
int AcNPtr = 0;
char **AcPtr;
int AcNIPtr = 0;
int **AcIPtr;
int AcNJPtr = 0;
int **AcJPtr;

/* ------------------------------------------------------------------------ */
/*  a_clean()   Free all previously allocated memory.                       */

void a_clean(void)
{
    alloc_actmp(0);
    alloc_actmp1(0);
    alloc_acx(0);
    alloc_acy(0);
    alloc_acy1(0);
    alloc_acz(0);
    alloc_act(0);
    alloc_acu(0);
    alloc_acv(0);
    alloc_acw(0);
    alloc_acxf(0);
    alloc_acyf(0);
    alloc_aczf(0);
    alloc_acuf(0);
    alloc_acvf(0);
    alloc_acm(0);
    alloc_acn(0);
    alloc_aci(0);
    alloc_acj(0);
    alloc_acr(0);
    alloc_acs(0);
    alloc_ack(0);
    alloc_acl(0);
    alloc_acc(0);
    alloc_acd(0);
    alloc_ace(0);
    alloc_acf(0);
    alloc_acns(0);
    alloc_acms(0);
    alloc_acptr(0);     
    alloc_aciptr(0);     
    alloc_acjptr(0);     
}

/* ------------------------------------------------------------------------ */
/*  alloc_actmp(n)  If n > 0 allocate AcTmp otherwise free.                 */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_actmp(int n)
{
    if (AcTmpN > 0) {
        free((char *)AcTmp);
        memrq(-AcTmpN,sizeof(double));
        AcTmpN = 0;
    }
    if (n > 0) {
        if (!(AcTmp = (double *)calloc(n,sizeof(double)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(double));
        AcTmpN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_actmp1(n)  If n > 0 allocate AcTmp1 otherwise free.               */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_actmp1(int n)
{
    if (AcTmp1N > 0) {
        free((char *)AcTmp1);
        memrq(-AcTmp1N,sizeof(double));
        AcTmp1N = 0;
    }
    if (n > 0) {
        if (!(AcTmp1 = (double *)calloc(n,sizeof(double)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(double));
        AcTmp1N = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acx(n)    If n > 0 allocate AcX otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acx(int n)
{
    if (AcXN > 0) {
        free((char *)AcX);
        memrq(-AcXN,sizeof(double));
        AcXN = 0;
    }
    if (n > 0) {
        if (!(AcX = (double *)calloc(n,sizeof(double)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(double));
        AcXN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acy(n)    If n > 0 allocate AcY otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acy(int n)
{
    if (AcYN > 0) {
        free((char *)AcY);
        memrq(-AcYN,sizeof(double));
        AcYN = 0;
    }
    if (n > 0) {
        if (!(AcY = (double *)calloc(n,sizeof(double)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(double));
        AcYN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acy1(n)    If n > 0 allocate AcY1 otherwise free.                 */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acy1(int n)
{
    if (AcY1N > 0) {
        free((char *)AcY1);
        memrq(-AcY1N,sizeof(double));
        AcY1N = 0;
    }
    if (n > 0) {
        if (!(AcY1 = (double *)calloc(n,sizeof(double)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(double));
        AcY1N = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acz(n)    If n > 0 allocate AcZ otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acz(int n)
{
    if (AcZN > 0) {
        free((char *)AcZ);
        memrq(-AcZN,sizeof(double));
        AcZN = 0;
    }
    if (n > 0) {
        if (!(AcZ = (double *)calloc(n,sizeof(double)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(double));
        AcZN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_act(n)    If n > 0 allocate AcT otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_act(int n)
{
    if (AcTN > 0) {
        free((char *)AcT);
        memrq(-AcTN,sizeof(double));
        AcTN = 0;
    }
    if (n > 0) {
        if (!(AcT = (double *)calloc(n,sizeof(double)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(double));
        AcTN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acu(n)    If n > 0 allocate AcU otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acu(int n)
{
    if (AcUN > 0) {
        free((char *)AcU);
        memrq(-AcUN,sizeof(double));
        AcUN = 0;
    }
    if (n > 0) {
        if (!(AcU = (double *)calloc(n,sizeof(double)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(double));
        AcUN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acv(n)    If n > 0 allocate AcV otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acv(int n)
{
    if (AcVN > 0) {
        free((char *)AcV);
        memrq(-AcVN,sizeof(double));
        AcVN = 0;
    }
    if (n > 0) {
        if (!(AcV = (double *)calloc(n,sizeof(double)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(double));
        AcVN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acw(n)    If n > 0 allocate AcW otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acw(int n)
{
    if (AcWN > 0) {
        free((char *)AcW);
        memrq(-AcWN,sizeof(double));
        AcWN = 0;
    }
    if (n > 0) {
        if (!(AcW = (double *)calloc(n,sizeof(double)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(double));
        AcWN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acxf(n)   If n > 0 allocate AcXF otherwise free.                  */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acxf(int n)
{
    if (AcXFN > 0) {
        free((char *)AcXF);
        memrq(-AcXFN,sizeof(float));
        AcXFN = 0;
    }
    if (n > 0) {
        if (!(AcXF = (float *)calloc(n,sizeof(float)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(float));
        AcXFN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acyf(n)   If n > 0 allocate AcYF otherwise free.                  */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acyf(int n)
{
    if (AcYFN > 0) {
        free((char *)AcYF);
        memrq(-AcYFN,sizeof(float));
        AcYFN = 0;
    }
    if (n > 0) {
        if (!(AcYF = (float *)calloc(n,sizeof(float)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(float));
        AcYFN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_aczf(n)   If n > 0 allocate AcZF otherwise free.                  */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_aczf(int n)
{
    if (AcZFN > 0) {
        free((char *)AcZF);
        memrq(-AcZFN,sizeof(float));
        AcZFN = 0;
    }
    if (n > 0) {
        if (!(AcZF = (float *)calloc(n,sizeof(float)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(float));
        AcZFN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acuf(n)   If n > 0 allocate AcUF otherwise free.                  */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acuf(int n)
{
    if (AcUFN > 0) {
        free((char *)AcUF);
        memrq(-AcUFN,sizeof(float));
        AcUFN = 0;
    }
    if (n > 0) {
        if (!(AcUF = (float *)calloc(n,sizeof(float)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(float));
        AcUFN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acvf(n)   If n > 0 allocate AcVF otherwise free.                  */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acvf(int n)
{
    if (AcVFN > 0) {
        free((char *)AcVF);
        memrq(-AcVFN,sizeof(float));
        AcVFN = 0;
    }
    if (n > 0) {
        if (!(AcVF = (float *)calloc(n,sizeof(float)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(float));
        AcVFN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acm(n)    If n > 0 allocate AcM otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acm(int n)
{
    if (AcNM > 0) {
        free((char *)AcM);
        memrq(-AcNM,sizeof(int));
        AcNM = 0;
    }
    if (n > 0) {
        if (!(AcM = (int *)calloc(n,sizeof(int)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(int));
        AcNM = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acn(n)    If n > 0 allocate AcN otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acn(int n)
{
    if (AcNN > 0) {
        free((char *)AcN);
        memrq(-AcNN,sizeof(int));
        AcNN = 0;
    }
    if (n > 0) {
        if (!(AcN = (int *)calloc(n,sizeof(int)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(int));
        AcNN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_aci(n)    If n > 0 allocate AcI otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_aci(int n)
{
    if (AcIN > 0) {
        free((char *)AcI);
        memrq(-AcIN,sizeof(int));
        AcIN = 0;
    }
    if (n > 0) {
        if (!(AcI = (int *)calloc(n,sizeof(int)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(int));
        AcIN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acj(n)    If n > 0 allocate AcJ otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acj(int n)
{
    if (AcJN > 0) {
        free((char *)AcJ);
        memrq(-AcJN,sizeof(int));
        AcJN = 0;
    }
    if (n > 0) {
        if (!(AcJ = (int *)calloc(n,sizeof(int)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(int));
        AcJN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acr(n)    If n > 0 allocate AcR otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acr(int n)
{
    if (AcRN > 0) {
        free((char *)AcR);
        memrq(-AcRN,sizeof(int));
        AcRN = 0;
    }
    if (n > 0) {
        if (!(AcR = (int *)calloc(n,sizeof(int)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(int));
        AcRN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acs(n)    If n > 0 allocate AcS otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acs(int n)
{
    if (AcSN > 0) {
        free((char *)AcS);
        memrq(-AcSN,sizeof(int));
        AcSN = 0;
    }
    if (n > 0) {
        if (!(AcS = (int *)calloc(n,sizeof(int)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(int));
        AcSN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_ack(n)    If n > 0 allocate AcK otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_ack(int n)
{
    if (AcKN > 0) {
        free((char *)AcK);
        memrq(-AcKN,sizeof(int));
        AcKN = 0;
    }
    if (n > 0) {
        if (!(AcK = (int *)calloc(n,sizeof(int)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(int));
        AcKN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acl(n)    If n > 0 allocate AcL otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acl(int n)
{
    if (AcLN > 0) {
        free((char *)AcL);
        memrq(-AcLN,sizeof(int));
        AcLN = 0;
    }
    if (n > 0) {
        if (!(AcL = (int *)calloc(n,sizeof(int)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(int));
        AcLN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acc(n)    If n > 0 allocate AcC otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acc(int n)
{
    if (AcCN > 0) {
        free((char *)AcC);
        memrq(-AcCN,sizeof(char));
        AcCN = 0;
    }
    if (n > 0) {
        if (!(AcC = (char *)calloc(n,sizeof(char)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(char));
        AcCN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acd(n)    If n > 0 allocate AcD otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acd(int n)
{
    if (AcDN > 0) {
        free((char *)AcD);
        memrq(-AcDN,sizeof(char));
        AcDN = 0;
    }
    if (n > 0) {
        if (!(AcD = (char *)calloc(n,sizeof(char)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(char));
        AcDN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_ace(n)    If n > 0 allocate AcE otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_ace(int n)
{
    if (AcEN > 0) {
        free((char *)AcE);
        memrq(-AcEN,sizeof(char));
        AcEN = 0;
    }
    if (n > 0) {
        if (!(AcE = (char *)calloc(n,sizeof(char)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(char));
        AcEN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acf(n)    If n > 0 allocate AcF otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acf(int n)
{
    if (AcFN > 0) {
        free((char *)AcF);
        memrq(-AcFN,sizeof(char));
        AcFN = 0;
    }
    if (n > 0) {
        if (!(AcF = (char *)calloc(n,sizeof(char)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(char));
        AcFN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acns(n)   If n > 0 allocate AcNS otherwise free.                  */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acns(int n)
{
    if (AcNNS > 0) {
        free((char *)AcNS);
        memrq(-AcNNS,sizeof(short));
        AcNNS = 0;
    }
    if (n > 0) {
        if (!(AcNS = (short *)calloc(n,sizeof(short)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(short));
        AcNNS = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acms(n)   If n > 0 allocate AcMS otherwise free.                  */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acms(int n)
{
    if (AcNMS > 0) {
        free((char *)AcMS);
        memrq(-AcNMS,sizeof(short));
        AcNMS = 0;
    }
    if (n > 0) {
        if (!(AcMS = (short *)calloc(n,sizeof(short)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(short));
        AcNMS = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acptr(n)  If n > 0 allocate AcPtr otherwise free.                 */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acptr(int n)
{
    register int i;

    if (AcNPtr > 0) {
        free((char **)AcPtr);
        memrq(-AcNPtr,sizeof(char *));
        AcNPtr = 0;
    }
    if (n > 0) {
        if (!(AcPtr = (char **)calloc(n,sizeof(char *)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(char *));
        AcNPtr = n;
        for (i = 0; i < n; ++i)
            AcPtr[i] = NULL;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_aciptr(n)  If n > 0 allocate AcIPtr otherwise free.               */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_aciptr(int n)
{
    register int i;

    if (AcNIPtr > 0) {
        free((char **)AcIPtr);
        memrq(-AcNIPtr,sizeof(int *));
        AcNIPtr = 0;
    }
    if (n > 0) {
        if (!(AcIPtr = (int **)calloc(n,sizeof(int *)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(int *));
        AcNIPtr = n;
        for (i = 0; i < n; ++i)
            AcIPtr[i] = NULL;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acjptr(n)  If n > 0 allocate AcJPtr otherwise free.               */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acjptr(int n)
{
    register int i;

    if (AcNJPtr > 0) {
        free((char **)AcJPtr);
        memrq(-AcNJPtr,sizeof(int *));
        AcNJPtr = 0;
    }
    if (n > 0) {
        if (!(AcJPtr = (int **)calloc(n,sizeof(int *)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(int *));
        AcNJPtr = n;
        for (i = 0; i < n; ++i)
            AcJPtr[i] = NULL;
    }
    return(0);
}













