/****************************************************************************/
/*  t_com                                                                   */
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
#include "t_var.h"
#include "t_gdat.h"
#include "tda_context.h"

/*  functions in t_com.c */

int com(TDAContext *ctx);
int comb_n(TDAContext *ctx, int n,int *com,int first); 
int comb_nm(TDAContext *ctx, int n,int m,int *com,int first);
int comb_nm1(TDAContext *ctx, int n,int m,int *com,int first);
int perm(TDAContext *ctx, int n,int *com,int *p,int *d,int first);
int partit(TDAContext *ctx, int *k,int *p,int last);
int partit1(TDAContext *ctx, int n,int k,int *c,int *d,int last);
int partit2(TDAContext *ctx, int n,int r,int *a,int *p,int *t,int first);
int pcyc(TDAContext *ctx);
int indep(TDAContext *ctx);


/* ------------------------------------------------------------------------ */
/*  com         Combinatorial patterns.                                     */
/*  ##                                                                      */
/*              com(                                                        */
/*                  opt = ...,  1 = all ordered n-tuples                    */
/*                              2 = all m-sets from {0,...,n-1}             */
/*                              3 = all m-tuples from ...                   */
/*                              4 = all permutations of 0,...,n-1           */
/*                              5 = all partitions of n                     */
/*                              6 = all partitions of n into m parts (m>1)  */
/*                              7 = all partitions of n into m 2 subsets    */
/*                  n=...,      dimension, def. 1                           */
/*                  m=...,      def. 1                                      */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int com(TDAContext *ctx)
{
    register int j;
    int err,first,last,n,m = 0,r,k;

    err = -1;
    n = 0;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Combinatorial patterns. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 3,1,1))       /* get parameters */
        goto COMFin;

    if (ctx->PMN < 1)
        ctx->PMN = 1;
    if (ctx->PMM < 1)
        ctx->PMM = 1;

    if (ctx->PMOPT > 7)
        ctx->PMOPT = 7;

    if (alloc_acn(ctx, imax(ctx, ctx->PMN,ctx->PMM) + 1))  
        goto COMFin;

    if (ctx->PMOPT == 1) {
        printf1(ctx, "Option 1: all n-tuples (n=%d).\n",ctx->PMN);
        m = ctx->PMN;
    }
    else if (ctx->PMOPT == 2) {
        if (ctx->PMM > ctx->PMN)
            ctx->PMM = ctx->PMN;
        m = ctx->PMM;
        printf1(ctx, "Option 2: all m-sets (n=%d, m=%d).\n",ctx->PMN,ctx->PMM);
    }
    else if (ctx->PMOPT == 3) {
        printf1(ctx, "Option 3: all m-tuples (n=%d, m=%d).\n",ctx->PMN,ctx->PMM);
        m = ctx->PMM;
    }
    else if (ctx->PMOPT == 4) {
        printf1(ctx, "Option 4: all permutations (n=%d).\n",ctx->PMN);
        if (alloc_aci(ctx, ctx->PMN + 1))  
            goto COMFin;
        if (alloc_acj(ctx, ctx->PMN + 1))  
            goto COMFin;
        m = ctx->PMN;
    }
    else if (ctx->PMOPT == 5) {
        printf1(ctx, "Option 5: all partitions of n=%d.\n",ctx->PMN);
    }
    else if (ctx->PMOPT == 6) {
        if (ctx->PMM > ctx->PMN)
            ctx->PMM = ctx->PMN;
        else if (ctx->PMM < 2)
            ctx->PMM = 2;
        printf1(ctx, "Option 6: all partitions of n=%d into m=%d parts.\n",ctx->PMN,ctx->PMM);
        if (alloc_aci(ctx, ctx->PMN + 1))  
            goto COMFin;
    }
    else if (ctx->PMOPT == 7) {
        if (ctx->PMM > ctx->PMN)
            ctx->PMM = ctx->PMN;
        else if (ctx->PMM < 1)
            ctx->PMM = 1;
        printf1(ctx, "Option 7: all partitions of n=%d into m=%d subsets.\n",ctx->PMN,ctx->PMM);
        if (alloc_aci(ctx, ctx->PMN + 1))  
            goto COMFin;
        if (alloc_acj(ctx, ctx->PMN + 1))  
            goto COMFin;
    }

    if (ctx->PMOPT < 5) {
        n = 0;
        first = 1;

        while (1) {

            if (ctx->PMOPT == 1)
                r = comb_n(ctx, ctx->PMN,ctx->AcN,first);
            else if (ctx->PMOPT == 2)
                r = comb_nm(ctx, ctx->PMN,ctx->PMM,ctx->AcN,first);
            else if (ctx->PMOPT == 3)
                r = comb_nm1(ctx, ctx->PMN,ctx->PMM,ctx->AcN,first);
            else if (ctx->PMOPT == 4)
                r = perm(ctx, ctx->PMN,ctx->AcN,ctx->AcI,ctx->AcJ,first);

            if (r == 0)
                break;

            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,++n);
            for (j = 0; j < m; ++j)
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcN[j]);
            fprintf(ctx->PMFd,"\n");

            first = 0;
        }
    }
    else if (ctx->PMOPT == 5) {
        ctx->AcN[1] = ctx->PMN;
        k = 1;
        last = 1;
        n = 0;
        while (1) {
            r = partit(ctx, &k,ctx->AcN,last);           
         
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,++n);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,k);
            for (j = 1; j <= ctx->PMN; ++j)
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcN[j]);
            fprintf(ctx->PMFd,"\n");

            if (r)
                break;
            last = 0;
        }
    }
    else if (ctx->PMOPT == 6) {
        last = 1;
        n = 0;
        while (1) {
            r = partit1(ctx, ctx->PMN,ctx->PMM,ctx->AcN,ctx->AcI,last);
            if (r)
                break;
         
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,++n);
            for (j = 1; j <= ctx->PMM; ++j)
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcN[j]);
            fprintf(ctx->PMFd,"\n");
            last = 0;
        }
    }
    else if (ctx->PMOPT == 7) {
        last = 1;
        n = 0;
        while (1) {
            r = partit2(ctx, ctx->PMN,ctx->PMM,ctx->AcN,ctx->AcI,ctx->AcJ,last);
            if (r)
                break;
         
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,++n);
            for (j = 1; j <= ctx->PMN; ++j)
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcN[j]);
            fprintf(ctx->PMFd,"\n");
            last = 0;
        }
    }

    printf1(ctx, "%d records written to: %s\n",n,ctx->PMFdName);
    err = 0;
           
COMFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  comb_n(n,com,first)                                                     */
/*                                                                          */  
/*  Enumerates all orderes n-tuples. At first call, first = 1, afterwards   */
/*  first = 0. Returns 1 if a new pattern found, otherwise 0.               */

int comb_n(TDAContext *ctx, int n,int *com,int first)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j,fin;

    if (first) {
        for (i = 0; i < n; ++i)
            com[i] = 0;
        return(1);
    }
    fin = 0;
    for (i = 0; i < n; ++i) {
        if (com[i] < n - 1) {
            com[i] += 1;
            for (j = 0; j < i; ++j)  
                com[j] = 0;
            fin = 1;
            break;
        }
    }
    return(fin);
}

/* ------------------------------------------------------------------------ */
/*  comb_nm(n,m,com,first)                                                  */
/*                                                                          */  
/*  Generate all subsets of {0,...,n-1} with m elements.                    */
/*  Return 1 if new subset found, otherwise return 0.                       */

int comb_nm(TDAContext *ctx, int n,int m,int *com,int first)
{
    (void)ctx;        /* unused: the signature is shared */
    register int k,j,l;

    if (first) {
        for (k = 0; k < m; ++k)
            com[k] = n;
    }
    for (k = 1; k <= m; ++k) {
        j = m - k;
        if (com[j] < n - k) {
            l = com[j];
            for (k = j; k < m; ++k)  
                com[k] = ++l;
            return(1);
        }
    }
    for (k = 0; k < m; ++k) 
        com[k] = k;
    if (first)
        return(1);
    else
        return(0);
}

/* ------------------------------------------------------------------------ */
/*  comb_nm1(n,m,com,first)                                                 */
/*                                                                          */  
/*  Enumerates all m-tuples from {0,...,n-1} x {0,...,n-1} x ...            */
/*                                                                          */
/*  At first call, first = 1, afterwards                                    */
/*  first = 0. Returns 1 if a new pattern found, otherwise 0.               */

int comb_nm1(TDAContext *ctx, int n,int m,int *com,int first)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j,fin;

    if (first) {
        for (i = 0; i < m; ++i)
            com[i] = 0;
        return(1);
    }
    fin = 0;
    for (i = 0; i < m; ++i) {
        if (com[i] < n - 1) {
            com[i] += 1;
            for (j = 0; j < i; ++j)  
                com[j] = 0;
            fin = 1;
            break;
        }
    }
    return(fin);
}

/* ------------------------------------------------------------------------ */
/*  perm(n,com,p,d,first)                                                   */
/*                                                                          */  
/*  Creates all permutations of {0,...,n-1}                                 */
/*                                                                          */
/*  At first call, first = 1, afterwards                                    */
/*  first = 0. Returns 1 if a new pattern found, otherwise 0.               */
/*  p[] and d[] are integer arrays of dimension n.                          */

int perm(TDAContext *ctx, int n,int *com,int *p,int *d,int first)
{
    (void)ctx;        /* unused: the signature is shared */
    register int k,m,q,t;

    if (first) {
        for (k = 0; k < n; ++k) {
            com[k] = k;
            p[k] = 0;
            d[k] = 1;
        }
        return(1);
    }
    m = n - 1;
    k = 0;

    while (1) {
        q = p[m] + d[m];
        p[m] = q;
        if (q == m + 1)  
            d[m] = -1;
        else {
            if (q != 0)  
                break;  
            d[m] = 1;
            k++;
        }
        if (m <= 1) {
            q = 1;
            return(0);
        }               
        m--;
    }
    q += k;
    t = com[q - 1];
    com[q - 1] = com[q];
    com[q] = t;
    return(1);
}

/* -##--------------------------------------------------------------------- */
/*  partit(k,p,last)        creates all partitions of n.                    */
/*                                                                          */  
/*  Algorithm adapted from: CACM 371 (Partitions in natural order).         */
/*                                                                          */
/*  In order to create all partitions, partit() must first be called        */
/*  with p[1] = n, k = 1, last = 1. Afterwards with last = 0.               */
/*  When the last partition is generated, the functions returns 1.          */
/*  p[i] is indexed by i = 1,...n.                                          */

int partit(TDAContext *ctx, int *k,int *p,int last)
{
    register int t;

    if (last) {
        for (ctx->s_partit_m = 1; ctx->s_partit_m <= *k; ++ctx->s_partit_m) {
            if (p[ctx->s_partit_m] == 1) 
                goto PC;
        }
        ctx->s_partit_m = *k;
        goto PC;
    }
    t = *k - ctx->s_partit_m;
    *k = ctx->s_partit_m;
    p[ctx->s_partit_m] -= 1;
PA:
    if (p[*k] > t)
        goto PB;

    t -= p[*k];
    *k += 1;
    p[*k] = p[*k - 1];
    goto PA;

PB:
    *k += 1;
    p[*k] = t + 1;
    if (p[ctx->s_partit_m] != 1)
        ctx->s_partit_m = *k;
PC:
    if (p[ctx->s_partit_m] == 1)  
        ctx->s_partit_m--;
    if (ctx->s_partit_m == 0)
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  partit1(n,k,c,d,last)   creates all partitions                          */
/*                          n = c[1] + c[2] + ... + c[k], k > 1.            */
/*                          c[] and d[] are arrays indexed 1,...,k          */
/*                                                                          */  
/*  Algorithm adapted from: CACM 72.                                        */
/*                                                                          */
/*  In order to create all partitions, partit1() must first be called       */
/*  with last = 1. Then in all further calls last must be 0. If no          */
/*  more partitions the function returns 1.                                 */

int partit1(TDAContext *ctx, int n,int k,int *c,int *d,int last)
{
    (void)ctx;        /* unused: the signature is shared */
    register int j;

    if (last == 1) {
        c[k] = n - k + 1;
        for (j = 1; j < k; ++j)  
            c[j] = 1;
        return(0);
    }
    for (j = 1; j <= k; ++j)  
        d[j] = c[j] - 1;
    j = k;
  
    while (d[j] <= 0) {
        j--;
        if (j == 1)  
            return(1);
    }
    d[j] = 0;
    d[j - 1] += 1;
    d[k] = c[j] - 2;
    for (j = 1; j <= k; ++j)
        c[j] = d[j] + 1;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  partit2(n,r,a,p,t,last)                                                 */
/*                                                                          */
/*  creates all partitions of (1,...,n) into r subsets. r > 3 required.     */
/*                                                                          */
/*  First call with first = 1, afterwards first = 0. The function returns   */
/*  0 if a new partition is found, otherwise 1.                             */
/*                                                                          */  
/*  Algorithm adapted from: CACM 477.                                       */
/*                                                                          */

int partit2(TDAContext *ctx, int n,int r,int *a,int *p,int *t,int first)
{
    register int j,l;
    int s,u,ns;

    if (first) {
        ctx->s_partit2_last = 0;
             
        if (r == 1) {
            for (j = 1; j <= n; ++j)
                a[j] = 1;
            ctx->s_partit2_last = 1;
            return(0);
        }
        else if (r == 2) {
            ctx->s_partit2_z = n / 2;
            if (2 * ctx->s_partit2_z == n)
                ctx->s_partit2_even = 1;                       
            else
                ctx->s_partit2_even = 0;                           

            ctx->s_partit2_k = 1;
            p[1] = 1;
            a[1] = 1;
            for (j = 2; j <= n; ++j)
                a[j] = 2;

            if (n == 2)
                ctx->s_partit2_last = 1;
            return(0);
        }
        else {
            ctx->s_partit2_k = n - r + 1;
            for (j = 1; j <= ctx->s_partit2_k; ++j)  
                a[j] = p[j] = 1;
            for (j = ctx->s_partit2_k + 1; j <= n; ++j)
                a[j] = 1 + j - ctx->s_partit2_k;
            ctx->s_partit2_i = ctx->s_partit2_k;
            t[ctx->s_partit2_k] = ctx->s_partit2_k - 1;
            t[ctx->s_partit2_k - 1] = 0;
            ctx->s_partit2_z = 1;
        }                  
        goto PFIN;
    }
    if (ctx->s_partit2_last)  
        return(1);
                            
    if (r == 2) {                     
        while (1) {
            s = n;
            for (j = ctx->s_partit2_k; j >= 1; --j) {
                if (ctx->s_partit2_even && ctx->s_partit2_k == ctx->s_partit2_z && j == 1)
                    break;

                p[j] += 1;
                if (p[j] <= s) {
                    for (l = j + 1; l <= ctx->s_partit2_k; ++l)
                        p[l] = p[l - 1] + 1;
                    goto P2NXT;
                }
                s--;
            }
            if (ctx->s_partit2_k == ctx->s_partit2_z) {
                ctx->s_partit2_last = 1;
                return(1);
            }
            ctx->s_partit2_k++;
            p[1] = 1;
            for (l = 2; l <= ctx->s_partit2_k; ++l)
                p[l] = p[l - 1] + 1;
            break;        
        }
P2NXT:
        for (l = 1; l <= n; ++l)
            a[l] = 2; 
        for (l = 1; l <= ctx->s_partit2_k; ++l)
            a[p[l]] = 1;
        return(0);
    }
          
    /*  algorithm for r > 2 */

    s = ctx->s_partit2_i;
    u = 0;

    switch (p[ctx->s_partit2_i]) {
      case 1:
        ns = a[ctx->s_partit2_i] = ctx->s_partit2_z = 2;
        p[ctx->s_partit2_i] = 7;
        if (ctx->s_partit2_i == ctx->s_partit2_k) {
            u = ctx->s_partit2_k = ctx->s_partit2_k + 1;
            a[u] = 1;
            p[ctx->s_partit2_k] = 6;
        }
        goto NOGA;

      case 2:
        ns = a[ctx->s_partit2_i] = ctx->s_partit2_z = ctx->s_partit2_z - 1;
        u = ctx->s_partit2_k;
        a[ctx->s_partit2_k] = ctx->s_partit2_z + 1;
        ctx->s_partit2_k--;
        if (ctx->s_partit2_z == 2) {
            p[ctx->s_partit2_i] = 7;
            goto NOGA;
        }
        p[ctx->s_partit2_i] = 3;
        goto OFRA;
      

      case 3:
        ns = a[ctx->s_partit2_i] = a[ctx->s_partit2_i] - 1;
        if (ns != 2)
            goto OFRA;
        p[ctx->s_partit2_i] = 7;
        goto NOGA;

      case 4:
        u = ctx->s_partit2_k;
        a[u] = ctx->s_partit2_z;
        ctx->s_partit2_z--;
        ctx->s_partit2_k--;
        TDA_FALLTHROUGH;

      case 5:
        ns = a[ctx->s_partit2_i] = 1;
        p[ctx->s_partit2_i] = 6;
        goto NOGA;

      case 6:
        if (ctx->s_partit2_z == r) {
            ns = a[ctx->s_partit2_i] = r;
            p[ctx->s_partit2_i] = 3;
        }
        else {
            ns = a[ctx->s_partit2_i] = ctx->s_partit2_z = ctx->s_partit2_z + 1;
            p[ctx->s_partit2_i] = 2;
            u = ctx->s_partit2_k = ctx->s_partit2_k + 1;
            a[ctx->s_partit2_k] = 1;
            p[ctx->s_partit2_k] = 6;
        }
        goto OFRA;

      case 7:
        ns = a[ctx->s_partit2_i] = a[ctx->s_partit2_i] + 1;
        if (ns >= ctx->s_partit2_z) {
            if (ctx->s_partit2_z == r)
                p[ctx->s_partit2_i] = 5;
            else if (a[ctx->s_partit2_i] == ctx->s_partit2_z + 1) {
                ctx->s_partit2_z++;
                p[ctx->s_partit2_i] = 4;
                u = ctx->s_partit2_k = ctx->s_partit2_k + 1;
                a[ctx->s_partit2_k] = 1;
                p[ctx->s_partit2_k] = 6;
            }
        }
        goto OFRA;
    }

NOGA:
    if (ctx->s_partit2_i == ctx->s_partit2_k) {
        ctx->s_partit2_i = t[ctx->s_partit2_i];
        goto PFIN;     
    }
    if (t[ctx->s_partit2_i] < 1) {
        if (-t[ctx->s_partit2_i] != ctx->s_partit2_i - 1)  
            t[ctx->s_partit2_i - 1] = t[ctx->s_partit2_i];
        t[ctx->s_partit2_i] = ctx->s_partit2_i - 1;
    }
    if (ctx->s_partit2_i != ctx->s_partit2_k - 1) {
        t[ctx->s_partit2_k] = ctx->s_partit2_k - 1;
        t[ctx->s_partit2_k - 1] = -ctx->s_partit2_i - 1;
    }   
    t[ctx->s_partit2_i + 1] = t[ctx->s_partit2_i];
    ctx->s_partit2_i = ctx->s_partit2_k; 
    goto PFIN;       

OFRA:
    if (ctx->s_partit2_i == ctx->s_partit2_k)
        goto PFIN;       
    t[ctx->s_partit2_k] = ctx->s_partit2_k - 1;
    if (ctx->s_partit2_i != ctx->s_partit2_k - 1)
        t[ctx->s_partit2_k - 1] = -ctx->s_partit2_i;
    ctx->s_partit2_i = ctx->s_partit2_k;

PFIN:
    if (ctx->s_partit2_i == 1)
        ctx->s_partit2_last = 1;
    return(0);
}


/* ------------------------------------------------------------------------ */
/*  pcyc    Calculating cycles of permutations.                             */
/*                                                                          */
/*          pcyc(                                                           */
/*              df=...,         output file (required)                      */
/*              nfmt=...,       integer print format, def. 4                */
/*          ) = varlist;        varlist (required)                          */  
/*                                                                          */
/*  Return 0 if successful, otherwise -1.                                   */

int pcyc(TDAContext *ctx)
{
    register int i,j,k,l,ll;
    int err,n,s,nc,ncc,nic,iflag,lcc,lic,nrec;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Cycles of permutations. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,4,1))     /* get parameters */
        goto PCYCFin;

    printf1(ctx, "Number of positions: %d\n",ctx->PMNV);
    if (ctx->PMF1Def == 0) { 
        printf1(ctx, "Error: need an output file.\n");
        goto PCYCFin;
    }

    if (alloc_acn(ctx, ctx->PMNV + 1))  
        goto PCYCFin;
    if (alloc_ack(ctx, ctx->PMNV + 1))  
        goto PCYCFin;

    if (alloc_acr(ctx, ctx->PMNV + 1))  
        goto PCYCFin;
    if (alloc_acs(ctx, ctx->PMNV * ctx->PMNV + 1))  
        goto PCYCFin;

    nrec = 0;
    for (i = 0; i < ctx->NOC; ++i) {

        iflag = n = 0;
        for (j = 1; j <= ctx->PMNV; ++j)  
            ctx->AcK[j] = 0;

        for (j = 1; j <= ctx->PMNV; ++j) {
            s = (int)get_data(ctx, ctx->PMVIdx[j - 1],i);
            if (s >= 1 && s <= ctx->PMNV) {
                n++;
                if (ctx->AcK[s]) {
                    iflag = 1;
                    break;
                }
                ctx->AcN[j] = s;
                ctx->AcK[s] = 1;
            }
            else
                ctx->AcN[j] = 0;
        }
        if (iflag) {
            printf1(ctx, "Inconsistent data in case %d.\n",i + 1);
            continue;
        }
        nc = ncc = nic = 0;

        for (j = 1; j <= ctx->PMNV; ++j)  
            ctx->AcK[j] = 0;

        iflag = 0;                      /* begin with incomplete cycles */
        while (iflag == 0) {
            iflag = 1;
            for (j = 1; j <= ctx->PMNV; ++j) {
                if (ctx->AcK[j])
                    continue;

                if (ctx->AcN[j] == 0) {
                    iflag = 0;
                    ll = 0;
                    ctx->AcK[j] = 1;
                    ctx->AcS[nc * ctx->PMNV + ll] = j;
                    ll++;
                    while (1) {
                        k = 0;
                        for (l = 1; l <= ctx->PMNV; ++l) {
                            if (ctx->AcK[l] == 0 && ctx->AcN[l] == j) {
                                k = l;
                                break;
                            }
                        }
                        if (k == 0)    
                            break;
                        
                        ctx->AcK[k] = 1;
                        ctx->AcS[nc * ctx->PMNV + ll] = k;
                        ll++;
                        j = k;
                    }
                    ctx->AcR[nc] = ll;
                    nic++;
                    nc++;
                    break;
                }
            }
        }
        iflag = 0;                      /* continue with complete cycles */
        while (iflag == 0) {           
            iflag = 1;
            for (j = 1; j <= ctx->PMNV; ++j) {
                if (ctx->AcK[j])
                    continue;
                iflag = 0;
                ll = 0;
                ctx->AcK[j] = 1;
                ctx->AcS[nc * ctx->PMNV + ll] = j;
                ll++;
                while (1) {
                    k = ctx->AcN[j];
                    if (ctx->AcK[k])    
                        break;
                    ctx->AcK[k] = 1;
                    ctx->AcS[nc * ctx->PMNV + ll] = k;
                    ll++;
                    j = k;
                }
                ctx->AcR[nc] = ll;
                ncc++;
                nc++;
                break;
            }
        }

        /* change order of incomplete cycles */

        for (j = 0; j < nic; ++j) {
            ll = ctx->AcR[j];
            for (l = 0; l < ll; ++l)
                ctx->AcK[l] = ctx->AcS[j * ctx->PMNV + l];
            for (l = 0; l < ll; ++l)
                ctx->AcS[j * ctx->PMNV + l] = ctx->AcK[ll - l - 1];
        }

        /* calculate length of cycles */

        lcc = lic = 0;
        for (j = 0; j < nic; ++j)   
            lic += (ctx->AcR[j] - 1);
        for (j = nic; j < nc; ++j)   
            lcc += (ctx->AcR[j] - 1);
    
        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,i + 1);
        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,n);
        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ncc);
        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,nic);
        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,lcc);
        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,lic);

        for (j = 1; j <= ctx->PMNV; ++j)
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->AcN[j]);

        fprintf(ctx->PMF1d,"  ");
        for (j = nic; j < nc; ++j) {
            fprintf(ctx->PMF1d,"(%d",ctx->AcS[j * ctx->PMNV]);
            ll = ctx->AcR[j];
            for (l = 1; l < ll; ++l)
                fprintf(ctx->PMF1d,",%d",ctx->AcS[j * ctx->PMNV + l]);
            fprintf(ctx->PMF1d,")");
        }
        for (j = 0; j < nic; ++j) {
            fprintf(ctx->PMF1d,"[%d",ctx->AcS[j * ctx->PMNV]);
            ll = ctx->AcR[j];
            for (l = 1; l < ll; ++l)
                fprintf(ctx->PMF1d,",%d",ctx->AcS[j * ctx->PMNV + l]);
            fprintf(ctx->PMF1d,"]");
        }
        fprintf(ctx->PMF1d,"\n");
        nrec++;
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);
    err = 0;

PCYCFin:
    p_clean(ctx);      
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  indep   Delta-independence of two variables.                            */
/*                                                                          */
/*          indep(                                                          */
/*              w=...,          variables specifying weights                */
/*              y=...,          partition of Y's property space             */

/*                                                                          */
/*          ) = X,Y;            two variables (integer-valued)              */  
/*                                                                          */
/*  Return 0 if successful, otherwise -1.                                   */

int indep(TDAContext *ctx)
{
    register int i,j;
    int err,ix,iy,iw,x,y,xmin,xmax,ymin,ymax,nx,ny,np,ip; 
    double w,wsum,delta,fx,fy,fysum;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Delta-independence of two variables. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,4,1))     /* get parameters */
        goto INDEPFin;

    if (ctx->PMNV != 2) {
        printf1(ctx, "Error: need exactly two variables on right-hand side.\n");
        goto INDEPFin;
    }
    ix = ctx->PMVIdx[0];
    iy = ctx->PMVIdx[1];
    iw = ctx->PMWVar;

    xmin = xmax = (int)get_data(ctx, ix,0);
    ymin = ymax = (int)get_data(ctx, iy,0);

    for (i = 1; i < ctx->NOC; ++i) {
        x = (int)get_data(ctx, ix,i);
        y = (int)get_data(ctx, iy,i);
        xmin = imin(ctx, xmin,x);
        xmax = imax(ctx, xmax,x);
        ymin = imin(ctx, ymin,y);
        ymax = imax(ctx, ymax,y);
    }
    printf1(ctx, "Range of X: %d to %d\n",xmin,xmax);
    printf1(ctx, "Range of Y: %d to %d\n",ymin,ymax);

    if (xmin < 0 || ymin < 0) {
        printf1(ctx, "Error: values of X and Y must not be negative.\n");
        goto INDEPFin;
    }
    nx = xmax - xmin + 1;
    ny = ymax - ymin + 1;

    /* create partition of Y's property space */ 

    if (alloc_ack(ctx, ny + 1))  
        goto INDEPFin;
                 
    if (ctx->PMNTP > 0) {
        np = 0;
        j = ymin;
        for (i = 0; i < ctx->PMNTP; ++i) {
            y = (int)ctx->PMTP[i];
            if (y < ymin)
                continue;
            np++;
            while (j <= y) {
                if (j > ymax)
                    break;
                ctx->AcK[j++ - ymin] = np;
            }
            if (j > ymax)
                break;
        }
    }
    else {
        np = 0;
        for (j = ymin; j <= ymax; ++j)
            ctx->AcK[j - ymin] = ++np;
    }
    printf1(ctx, "Partition of Y's property space into %d subsets\n",np);

    for (i = 1; i <= np; ++i) {
        printf1(ctx, "%5d : ",i);
        for (j = ymin; j <= ymax; ++j) {
            if (ctx->AcK[j - ymin] == i)
                printf1(ctx, "%d ",j);
        }
        newline(ctx);
    }



    if (alloc_acx(ctx, nx * ny + 1))  
        goto INDEPFin;
    if (alloc_acy(ctx, ny + 1))  
        goto INDEPFin;

    wsum = 0.0;
    w = 1.0;
    for (i = 0; i < ctx->NOC; ++i) {
        x = (int)get_data(ctx, ix,i);
        y = (int)get_data(ctx, iy,i);

        if (iw >= 0) {
            w = get_data(ctx, iw,i);
            if (w < 0.0) {
                printf1(ctx, "Error: weights must not be negative.\n");
                goto INDEPFin;
            }
        }
        ctx->AcX[(x - xmin) * ny + (y - ymin)] += w;
        ctx->AcY[y - ymin] += w;
        wsum += w;
    }
    printf1(ctx, "Sum of frequencies: %g\n",wsum);
    if (wsum <= ctx->EPSI1) {
        err = 0;
        goto INDEPFin;
    }



    /* do separately for all subsets of partition of Y's property space */ 

    for (ip = 1; ip <= np; ++ip) {

        printf1(ctx, "Y-subset %d : ",ip);
        for (j = ymin; j <= ymax; ++j) {
            if (ctx->AcK[j - ymin] == ip)
                printf1(ctx, "%d ",j);
        }
        newline(ctx);

        /* begin with one-element X-sets */
       
        delta = 0.0;

        for (i = xmin; i <= xmax; ++i) {
            fx = 0.0;
            for (j = ymin; j <= ymax; ++j) 
                fx += ctx->AcX[(i - xmin) * ny + (j - ymin)];
            fx /= wsum;

            fy = fysum = 0.0;

            for (j = ymin; j <= ymax; ++j) {
                if (ctx->AcK[j - ymin] != ip)
                    continue;

                fy += ctx->AcX[(i - xmin) * ny + (j - ymin)];
                fysum += ctx->AcY[j - ymin];
            }
            if (fysum == 0.0)
                continue;
            fy /= fysum;

            delta = dmax(ctx, delta,fabs(fy - fx));

        }
        printf1(ctx, "Maximal delta for 1-element sets: %g\n",delta);

#ifdef TDA_R_PACKAGE
        /* one row per Y-subset: subset number, smallest and largest y
           value in it, delta -- next to the print so text and export
           can never drift apart */
        {
            int jj, ylo1 = 0, yhi1 = 0, seen = 0;
            for (jj = ymin; jj <= ymax; ++jj)
                if (ctx->AcK[jj - ymin] == ip) {
                    if (!seen) { ylo1 = jj; seen = 1; }
                    yhi1 = jj;
                }
            tda_export_cell(ctx, "indep.deltas", (double)ip);
            tda_export_cell(ctx, "indep.deltas", (double)ylo1);
            tda_export_cell(ctx, "indep.deltas", (double)yhi1);
            tda_export_cell(ctx, "indep.deltas", delta);
            tda_export_endrow(ctx, "indep.deltas");
        }
#endif

    }



    err = 0;

INDEPFin:
    p_clean(ctx);      
    return(err);
}



