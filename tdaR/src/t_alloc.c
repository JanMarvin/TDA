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
#include "tda_context.h"

/*  functions in t_alloc.c */

void a_clean(TDAContext *ctx);
int alloc_actmp(TDAContext *ctx, int n);
int alloc_actmp1(TDAContext *ctx, int n);
int alloc_acx(TDAContext *ctx, int n);
int alloc_acy(TDAContext *ctx, int n);
int alloc_acy1(TDAContext *ctx, int n);
int alloc_acz(TDAContext *ctx, int n);
int alloc_act(TDAContext *ctx, int n);
int alloc_acu(TDAContext *ctx, int n);
int alloc_acv(TDAContext *ctx, int n);
int alloc_acw(TDAContext *ctx, int n);
int alloc_acxf(TDAContext *ctx, int n);
int alloc_acyf(TDAContext *ctx, int n);
int alloc_aczf(TDAContext *ctx, int n);
int alloc_acuf(TDAContext *ctx, int n);
int alloc_acvf(TDAContext *ctx, int n);
int alloc_acm(TDAContext *ctx, int n);
int alloc_acn(TDAContext *ctx, int n);
int alloc_aci(TDAContext *ctx, int n);
int alloc_acj(TDAContext *ctx, int n);
int alloc_acr(TDAContext *ctx, int n);
int alloc_acs(TDAContext *ctx, int n);
int alloc_ack(TDAContext *ctx, int n);
int alloc_acl(TDAContext *ctx, int n);
int alloc_acc(TDAContext *ctx, int n);
int alloc_acd(TDAContext *ctx, int n);
int alloc_ace(TDAContext *ctx, int n);
int alloc_acf(TDAContext *ctx, int n);
int alloc_acns(TDAContext *ctx, int n);
int alloc_acms(TDAContext *ctx, int n);
int alloc_acptr(TDAContext *ctx, int n);
int alloc_aciptr(TDAContext *ctx, int n);
int alloc_acjptr(TDAContext *ctx, int n);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */


/* ------------------------------------------------------------------------ */
/*  a_clean()   Free all previously allocated memory.                       */

void a_clean(TDAContext *ctx)
{
    alloc_actmp(ctx, 0);
    alloc_actmp1(ctx, 0);
    alloc_acx(ctx, 0);
    alloc_acy(ctx, 0);
    alloc_acy1(ctx, 0);
    alloc_acz(ctx, 0);
    alloc_act(ctx, 0);
    alloc_acu(ctx, 0);
    alloc_acv(ctx, 0);
    alloc_acw(ctx, 0);
    alloc_acxf(ctx, 0);
    alloc_acyf(ctx, 0);
    alloc_aczf(ctx, 0);
    alloc_acuf(ctx, 0);
    alloc_acvf(ctx, 0);
    alloc_acm(ctx, 0);
    alloc_acn(ctx, 0);
    alloc_aci(ctx, 0);
    alloc_acj(ctx, 0);
    alloc_acr(ctx, 0);
    alloc_acs(ctx, 0);
    alloc_ack(ctx, 0);
    alloc_acl(ctx, 0);
    alloc_acc(ctx, 0);
    alloc_acd(ctx, 0);
    alloc_ace(ctx, 0);
    alloc_acf(ctx, 0);
    alloc_acns(ctx, 0);
    alloc_acms(ctx, 0);
    alloc_acptr(ctx, 0);     
    alloc_aciptr(ctx, 0);     
    alloc_acjptr(ctx, 0);     
}

/* ------------------------------------------------------------------------ */
/*  alloc_actmp(n)  If n > 0 allocate AcTmp otherwise free.                 */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_actmp(TDAContext *ctx, int n)
{
    if (ctx->AcTmpN > 0) {
        free((char *)ctx->AcTmp);
        memrq(ctx, -ctx->AcTmpN,sizeof(double));
        ctx->AcTmpN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcTmp = (double *)calloc((size_t)(n),sizeof(double)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(double));
        ctx->AcTmpN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_actmp1(n)  If n > 0 allocate AcTmp1 otherwise free.               */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_actmp1(TDAContext *ctx, int n)
{
    if (ctx->AcTmp1N > 0) {
        free((char *)ctx->AcTmp1);
        memrq(ctx, -ctx->AcTmp1N,sizeof(double));
        ctx->AcTmp1N = 0;
    }
    if (n > 0) {
        if (!(ctx->AcTmp1 = (double *)calloc((size_t)(n),sizeof(double)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(double));
        ctx->AcTmp1N = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acx(n)    If n > 0 allocate AcX otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acx(TDAContext *ctx, int n)
{
    if (ctx->AcXN > 0) {
        free((char *)ctx->AcX);
        memrq(ctx, -ctx->AcXN,sizeof(double));
        ctx->AcXN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcX = (double *)calloc((size_t)(n),sizeof(double)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(double));
        ctx->AcXN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acy(n)    If n > 0 allocate AcY otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acy(TDAContext *ctx, int n)
{
    if (ctx->AcYN > 0) {
        free((char *)ctx->AcY);
        memrq(ctx, -ctx->AcYN,sizeof(double));
        ctx->AcYN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcY = (double *)calloc((size_t)(n),sizeof(double)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(double));
        ctx->AcYN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acy1(n)    If n > 0 allocate AcY1 otherwise free.                 */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acy1(TDAContext *ctx, int n)
{
    if (ctx->AcY1N > 0) {
        free((char *)ctx->AcY1);
        memrq(ctx, -ctx->AcY1N,sizeof(double));
        ctx->AcY1N = 0;
    }
    if (n > 0) {
        if (!(ctx->AcY1 = (double *)calloc((size_t)(n),sizeof(double)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(double));
        ctx->AcY1N = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acz(n)    If n > 0 allocate AcZ otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acz(TDAContext *ctx, int n)
{
    if (ctx->AcZN > 0) {
        free((char *)ctx->AcZ);
        memrq(ctx, -ctx->AcZN,sizeof(double));
        ctx->AcZN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcZ = (double *)calloc((size_t)(n),sizeof(double)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(double));
        ctx->AcZN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_act(n)    If n > 0 allocate AcT otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_act(TDAContext *ctx, int n)
{
    if (ctx->AcTN > 0) {
        free((char *)ctx->AcT);
        memrq(ctx, -ctx->AcTN,sizeof(double));
        ctx->AcTN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcT = (double *)calloc((size_t)(n),sizeof(double)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(double));
        ctx->AcTN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acu(n)    If n > 0 allocate AcU otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acu(TDAContext *ctx, int n)
{
    if (ctx->AcUN > 0) {
        free((char *)ctx->AcU);
        memrq(ctx, -ctx->AcUN,sizeof(double));
        ctx->AcUN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcU = (double *)calloc((size_t)(n),sizeof(double)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(double));
        ctx->AcUN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acv(n)    If n > 0 allocate AcV otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acv(TDAContext *ctx, int n)
{
    if (ctx->AcVN > 0) {
        free((char *)ctx->AcV);
        memrq(ctx, -ctx->AcVN,sizeof(double));
        ctx->AcVN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcV = (double *)calloc((size_t)(n),sizeof(double)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(double));
        ctx->AcVN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acw(n)    If n > 0 allocate AcW otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acw(TDAContext *ctx, int n)
{
    if (ctx->AcWN > 0) {
        free((char *)ctx->AcW);
        memrq(ctx, -ctx->AcWN,sizeof(double));
        ctx->AcWN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcW = (double *)calloc((size_t)(n),sizeof(double)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(double));
        ctx->AcWN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acxf(n)   If n > 0 allocate AcXF otherwise free.                  */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acxf(TDAContext *ctx, int n)
{
    if (ctx->AcXFN > 0) {
        free((char *)ctx->AcXF);
        memrq(ctx, -ctx->AcXFN,sizeof(float));
        ctx->AcXFN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcXF = (float *)calloc((size_t)(n),sizeof(float)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(float));
        ctx->AcXFN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acyf(n)   If n > 0 allocate AcYF otherwise free.                  */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acyf(TDAContext *ctx, int n)
{
    if (ctx->AcYFN > 0) {
        free((char *)ctx->AcYF);
        memrq(ctx, -ctx->AcYFN,sizeof(float));
        ctx->AcYFN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcYF = (float *)calloc((size_t)(n),sizeof(float)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(float));
        ctx->AcYFN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_aczf(n)   If n > 0 allocate AcZF otherwise free.                  */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_aczf(TDAContext *ctx, int n)
{
    if (ctx->AcZFN > 0) {
        free((char *)ctx->AcZF);
        memrq(ctx, -ctx->AcZFN,sizeof(float));
        ctx->AcZFN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcZF = (float *)calloc((size_t)(n),sizeof(float)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(float));
        ctx->AcZFN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acuf(n)   If n > 0 allocate AcUF otherwise free.                  */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acuf(TDAContext *ctx, int n)
{
    if (ctx->AcUFN > 0) {
        free((char *)ctx->AcUF);
        memrq(ctx, -ctx->AcUFN,sizeof(float));
        ctx->AcUFN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcUF = (float *)calloc((size_t)(n),sizeof(float)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(float));
        ctx->AcUFN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acvf(n)   If n > 0 allocate AcVF otherwise free.                  */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acvf(TDAContext *ctx, int n)
{
    if (ctx->AcVFN > 0) {
        free((char *)ctx->AcVF);
        memrq(ctx, -ctx->AcVFN,sizeof(float));
        ctx->AcVFN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcVF = (float *)calloc((size_t)(n),sizeof(float)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(float));
        ctx->AcVFN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acm(n)    If n > 0 allocate AcM otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acm(TDAContext *ctx, int n)
{
    if (ctx->AcNM > 0) {
        free((char *)ctx->AcM);
        memrq(ctx, -ctx->AcNM,sizeof(int));
        ctx->AcNM = 0;
    }
    if (n > 0) {
        if (!(ctx->AcM = (int *)calloc((size_t)(n),sizeof(int)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(int));
        ctx->AcNM = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acn(n)    If n > 0 allocate AcN otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acn(TDAContext *ctx, int n)
{
    if (ctx->AcNN > 0) {
        free((char *)ctx->AcN);
        memrq(ctx, -ctx->AcNN,sizeof(int));
        ctx->AcNN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcN = (int *)calloc((size_t)(n),sizeof(int)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(int));
        ctx->AcNN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_aci(n)    If n > 0 allocate AcI otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_aci(TDAContext *ctx, int n)
{
    if (ctx->AcIN > 0) {
        free((char *)ctx->AcI);
        memrq(ctx, -ctx->AcIN,sizeof(int));
        ctx->AcIN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcI = (int *)calloc((size_t)(n),sizeof(int)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(int));
        ctx->AcIN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acj(n)    If n > 0 allocate AcJ otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acj(TDAContext *ctx, int n)
{
    if (ctx->AcJN > 0) {
        free((char *)ctx->AcJ);
        memrq(ctx, -ctx->AcJN,sizeof(int));
        ctx->AcJN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcJ = (int *)calloc((size_t)(n),sizeof(int)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(int));
        ctx->AcJN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acr(n)    If n > 0 allocate AcR otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acr(TDAContext *ctx, int n)
{
    if (ctx->AcRN > 0) {
        free((char *)ctx->AcR);
        memrq(ctx, -ctx->AcRN,sizeof(int));
        ctx->AcRN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcR = (int *)calloc((size_t)(n),sizeof(int)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(int));
        ctx->AcRN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acs(n)    If n > 0 allocate AcS otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acs(TDAContext *ctx, int n)
{
    if (ctx->AcSN > 0) {
        free((char *)ctx->AcS);
        memrq(ctx, -ctx->AcSN,sizeof(int));
        ctx->AcSN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcS = (int *)calloc((size_t)(n),sizeof(int)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(int));
        ctx->AcSN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_ack(n)    If n > 0 allocate AcK otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_ack(TDAContext *ctx, int n)
{
    if (ctx->AcKN > 0) {
        free((char *)ctx->AcK);
        memrq(ctx, -ctx->AcKN,sizeof(int));
        ctx->AcKN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcK = (int *)calloc((size_t)(n),sizeof(int)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(int));
        ctx->AcKN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acl(n)    If n > 0 allocate AcL otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acl(TDAContext *ctx, int n)
{
    if (ctx->AcLN > 0) {
        free((char *)ctx->AcL);
        memrq(ctx, -ctx->AcLN,sizeof(int));
        ctx->AcLN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcL = (int *)calloc((size_t)(n),sizeof(int)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(int));
        ctx->AcLN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acc(n)    If n > 0 allocate AcC otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acc(TDAContext *ctx, int n)
{
    if (ctx->AcCN > 0) {
        free((char *)ctx->AcC);
        memrq(ctx, -ctx->AcCN,sizeof(char));
        ctx->AcCN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcC = (char *)calloc((size_t)(n),sizeof(char)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(char));
        ctx->AcCN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acd(n)    If n > 0 allocate AcD otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acd(TDAContext *ctx, int n)
{
    if (ctx->AcDN > 0) {
        free((char *)ctx->AcD);
        memrq(ctx, -ctx->AcDN,sizeof(char));
        ctx->AcDN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcD = (char *)calloc((size_t)(n),sizeof(char)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(char));
        ctx->AcDN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_ace(n)    If n > 0 allocate AcE otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_ace(TDAContext *ctx, int n)
{
    if (ctx->AcEN > 0) {
        free((char *)ctx->AcE);
        memrq(ctx, -ctx->AcEN,sizeof(char));
        ctx->AcEN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcE = (char *)calloc((size_t)(n),sizeof(char)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(char));
        ctx->AcEN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acf(n)    If n > 0 allocate AcF otherwise free.                   */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acf(TDAContext *ctx, int n)
{
    if (ctx->AcFN > 0) {
        free((char *)ctx->AcF);
        memrq(ctx, -ctx->AcFN,sizeof(char));
        ctx->AcFN = 0;
    }
    if (n > 0) {
        if (!(ctx->AcF = (char *)calloc((size_t)(n),sizeof(char)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(char));
        ctx->AcFN = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acns(n)   If n > 0 allocate AcNS otherwise free.                  */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acns(TDAContext *ctx, int n)
{
    if (ctx->AcNNS > 0) {
        free((char *)ctx->AcNS);
        memrq(ctx, -ctx->AcNNS,sizeof(short));
        ctx->AcNNS = 0;
    }
    if (n > 0) {
        if (!(ctx->AcNS = (short *)calloc((size_t)(n),sizeof(short)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(short));
        ctx->AcNNS = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acms(n)   If n > 0 allocate AcMS otherwise free.                  */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acms(TDAContext *ctx, int n)
{
    if (ctx->AcNMS > 0) {
        free((char *)ctx->AcMS);
        memrq(ctx, -ctx->AcNMS,sizeof(short));
        ctx->AcNMS = 0;
    }
    if (n > 0) {
        if (!(ctx->AcMS = (short *)calloc((size_t)(n),sizeof(short)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(short));
        ctx->AcNMS = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acptr(n)  If n > 0 allocate AcPtr otherwise free.                 */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acptr(TDAContext *ctx, int n)
{
    register int i;

    if (ctx->AcNPtr > 0) {
        free((char **)ctx->AcPtr);
        memrq(ctx, -ctx->AcNPtr,sizeof(char *));
        ctx->AcNPtr = 0;
    }
    if (n > 0) {
        if (!(ctx->AcPtr = (char **)calloc((size_t)(n),sizeof(char *)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(char *));
        ctx->AcNPtr = n;
        for (i = 0; i < n; ++i)
            ctx->AcPtr[i] = NULL;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_aciptr(n)  If n > 0 allocate AcIPtr otherwise free.               */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_aciptr(TDAContext *ctx, int n)
{
    register int i;

    if (ctx->AcNIPtr > 0) {
        free((char **)ctx->AcIPtr);
        memrq(ctx, -ctx->AcNIPtr,sizeof(int *));
        ctx->AcNIPtr = 0;
    }
    if (n > 0) {
        if (!(ctx->AcIPtr = (int **)calloc((size_t)(n),sizeof(int *)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(int *));
        ctx->AcNIPtr = n;
        for (i = 0; i < n; ++i)
            ctx->AcIPtr[i] = NULL;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_acjptr(n)  If n > 0 allocate AcJPtr otherwise free.               */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_acjptr(TDAContext *ctx, int n)
{
    register int i;

    if (ctx->AcNJPtr > 0) {
        free((char **)ctx->AcJPtr);
        memrq(ctx, -ctx->AcNJPtr,sizeof(int *));
        ctx->AcNJPtr = 0;
    }
    if (n > 0) {
        if (!(ctx->AcJPtr = (int **)calloc((size_t)(n),sizeof(int *)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(int *));
        ctx->AcNJPtr = n;
        for (i = 0; i < n; ++i)
            ctx->AcJPtr[i] = NULL;
    }
    return(0);
}













