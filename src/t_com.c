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

/*  functions in t_com.c */

int com(void);
int comb_n(int n,int *com,int first); 
int comb_nm(int n,int m,int *com,int first);
int comb_nm1(int n,int m,int *com,int first);
int perm(int n,int *com,int *p,int *d,int first);
int partit(int *k,int *p,int last);
int partit1(int n,int k,int *c,int *d,int last);
int partit2(int n,int r,int *a,int *p,int *t,int first);
int pcyc(void);
int indep(void);


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

int com(void)
{
    register int j;
    int err,first,last,n,m,r,k;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Combinatorial patterns. Current memory: %d bytes.\n",MemReq);
        
    if (parm(CmdBuf + 3,1,1))       /* get parameters */
        goto COMFin;

    if (PMN < 1)
        PMN = 1;
    if (PMM < 1)
        PMM = 1;

    if (PMOPT > 7)
        PMOPT = 7;

    if (alloc_acn(imax(PMN,PMM) + 1))  
        goto COMFin;

    if (PMOPT == 1) {
        printf1("Option 1: all n-tuples (n=%d).\n",PMN);
        m = PMN;
    }
    else if (PMOPT == 2) {
        if (PMM > PMN)
            PMM = PMN;
        m = PMM;
        printf1("Option 2: all m-sets (n=%d, m=%d).\n",PMN,PMM);
    }
    else if (PMOPT == 3) {
        printf1("Option 3: all m-tuples (n=%d, m=%d).\n",PMN,PMM);
        m = PMM;
    }
    else if (PMOPT == 4) {
        printf1("Option 4: all permutations (n=%d).\n",PMN);
        if (alloc_aci(PMN + 1))  
            goto COMFin;
        if (alloc_acj(PMN + 1))  
            goto COMFin;
        m = PMN;
    }
    else if (PMOPT == 5) {
        printf1("Option 5: all partitions of n=%d.\n",PMN);
    }
    else if (PMOPT == 6) {
        if (PMM > PMN)
            PMM = PMN;
        else if (PMM < 2)
            PMM = 2;
        printf1("Option 6: all partitions of n=%d into m=%d parts.\n",PMN,PMM);
        if (alloc_aci(PMN + 1))  
            goto COMFin;
    }
    else if (PMOPT == 7) {
        if (PMM > PMN)
            PMM = PMN;
        else if (PMM < 1)
            PMM = 1;
        printf1("Option 7: all partitions of n=%d into m=%d subsets.\n",PMN,PMM);
        if (alloc_aci(PMN + 1))  
            goto COMFin;
        if (alloc_acj(PMN + 1))  
            goto COMFin;
    }

    if (PMOPT < 5) {
        n = 0;
        first = 1;

        while (1) {

            if (PMOPT == 1)
                r = comb_n(PMN,AcN,first);
            else if (PMOPT == 2)
                r = comb_nm(PMN,PMM,AcN,first);
            else if (PMOPT == 3)
                r = comb_nm1(PMN,PMM,AcN,first);
            else if (PMOPT == 4)
                r = perm(PMN,AcN,AcI,AcJ,first);

            if (r == 0)
                break;

            fprintf(PMFd,PMNFmtS,++n);
            for (j = 0; j < m; ++j)
                fprintf(PMFd,PMNFmtS,AcN[j]);
            fprintf(PMFd,"\n");

            first = 0;
        }
    }
    else if (PMOPT == 5) {
        AcN[1] = PMN;
        k = 1;
        last = 1;
        n = 0;
        while (1) {
            r = partit(&k,AcN,last);           
         
            fprintf(PMFd,PMNFmtS,++n);
            fprintf(PMFd,PMNFmtS,k);
            for (j = 1; j <= PMN; ++j)
                fprintf(PMFd,PMNFmtS,AcN[j]);
            fprintf(PMFd,"\n");

            if (r)
                break;
            last = 0;
        }
    }
    else if (PMOPT == 6) {
        last = 1;
        n = 0;
        while (1) {
            r = partit1(PMN,PMM,AcN,AcI,last);
            if (r)
                break;
         
            fprintf(PMFd,PMNFmtS,++n);
            for (j = 1; j <= PMM; ++j)
                fprintf(PMFd,PMNFmtS,AcN[j]);
            fprintf(PMFd,"\n");
            last = 0;
        }
    }
    else if (PMOPT == 7) {
        last = 1;
        n = 0;
        while (1) {
            r = partit2(PMN,PMM,AcN,AcI,AcJ,last);
            if (r)
                break;
         
            fprintf(PMFd,PMNFmtS,++n);
            for (j = 1; j <= PMN; ++j)
                fprintf(PMFd,PMNFmtS,AcN[j]);
            fprintf(PMFd,"\n");
            last = 0;
        }
    }

    printf1("%d records written to: %s\n",n,PMFdName);
    err = 0;
           
COMFin:       
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  comb_n(n,com,first)                                                     */
/*                                                                          */  
/*  Enumerates all orderes n-tuples. At first call, first = 1, afterwards   */
/*  first = 0. Returns 1 if a new pattern found, otherwise 0.               */

int comb_n(int n,int *com,int first)
{
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

int comb_nm(int n,int m,int *com,int first)
{
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

int comb_nm1(int n,int m,int *com,int first)
{
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

int perm(int n,int *com,int *p,int *d,int first)
{
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

int partit(int *k,int *p,int last)
{
    register int t;
    static int m;

    if (last) {
        for (m = 1; m <= *k; ++m) {
            if (p[m] == 1) 
                goto PC;
        }
        m = *k;
        goto PC;
    }
    t = *k - m;
    *k = m;
    p[m] -= 1;
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
    if (p[m] != 1)
        m = *k;
PC:
    if (p[m] == 1)  
        m--;
    if (m == 0)
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

int partit1(int n,int k,int *c,int *d,int last)
{
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

int partit2(int n,int r,int *a,int *p,int *t,int first)
{
    register int j,l;
    int s,u,os,ou,ns,nu;
    static int k,z,i,last,even;

    if (first) {
        last = 0;
             
        if (r == 1) {
            for (j = 1; j <= n; ++j)
                a[j] = 1;
            last = 1;
            return(0);
        }
        else if (r == 2) {
            z = n / 2;
            if (2 * z == n)
                even = 1;                       
            else
                even = 0;                           

            k = 1;
            p[1] = 1;
            a[1] = 1;
            for (j = 2; j <= n; ++j)
                a[j] = 2;

            if (n == 2)
                last = 1;
            return(0);
        }
        else {
            k = n - r + 1;
            for (j = 1; j <= k; ++j)  
                a[j] = p[j] = 1;
            for (j = k + 1; j <= n; ++j)
                a[j] = 1 + j - k;
            i = k;
            t[k] = k - 1;
            t[k - 1] = 0;
            z = 1;
        }                  
        goto PFIN;
    }
    if (last)  
        return(1);
                            
    if (r == 2) {                     
        while (1) {
            s = n;
            for (j = k; j >= 1; --j) {
                if (even && k == z && j == 1)
                    break;

                p[j] += 1;
                if (p[j] <= s) {
                    for (l = j + 1; l <= k; ++l)
                        p[l] = p[l - 1] + 1;
                    goto P2NXT;
                }
                s--;
            }
            if (k == z) {
                last = 1;
                return(1);
            }
            k++;
            p[1] = 1;
            for (l = 2; l <= k; ++l)
                p[l] = p[l - 1] + 1;
            break;        
        }
P2NXT:
        for (l = 1; l <= n; ++l)
            a[l] = 2; 
        for (l = 1; l <= k; ++l)
            a[p[l]] = 1;
        return(0);
    }
          
    /*  algorithm for r > 2 */

    s = i;
    os = a[s];
    u = 0;

    switch (p[i]) {
      case 1:
        ns = a[i] = z = 2;
        p[i] = 7;
        if (i == k) {
            u = k = k + 1;
            ou = a[u];
            nu = a[u] = 1;
            p[k] = 6;
        }
        goto NOGA;

      case 2:
        ns = a[i] = z = z - 1;
        u = k;
        ou = a[u];
        nu = a[k] = z + 1;
        k--;
        if (z == 2) {
            p[i] = 7;
            goto NOGA;
        }
        p[i] = 3;
        goto OFRA;
      

      case 3:
        ns = a[i] = a[i] - 1;
        if (ns != 2)
            goto OFRA;
        p[i] = 7;
        goto NOGA;

      case 4:
        u = k;
        ou = a[u];
        nu = a[u] = z;
        z--;
        k--;

      case 5:
        ns = a[i] = 1;
        p[i] = 6;
        goto NOGA;

      case 6:
        if (z == r) {
            ns = a[i] = r;
            p[i] = 3;
        }
        else {
            ns = a[i] = z = z + 1;
            p[i] = 2;
            u = k = k + 1;
            ou = a[u];
            nu = a[k] = 1;
            p[k] = 6;
        }
        goto OFRA;

      case 7:
        ns = a[i] = a[i] + 1;
        if (ns >= z) {
            if (z == r)
                p[i] = 5;
            else if (a[i] == z + 1) {
                z++;
                p[i] = 4;
                u = k = k + 1;
                ou = a[u];
                nu = a[k] = 1;
                p[k] = 6;
            }
        }
        goto OFRA;
    }

NOGA:
    if (i == k) {
        i = t[i];
        goto PFIN;     
    }
    if (t[i] < 1) {
        if (-t[i] != i - 1)  
            t[i - 1] = t[i];
        t[i] = i - 1;
    }
    if (i != k - 1) {
        t[k] = k - 1;
        t[k - 1] = -i - 1;
    }   
    t[i + 1] = t[i];
    i = k; 
    goto PFIN;       

OFRA:
    if (i == k)
        goto PFIN;       
    t[k] = k - 1;
    if (i != k - 1)
        t[k - 1] = -i;
    i = k;

PFIN:
    if (i == 1)
        last = 1;
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

int pcyc(void)
{
    register int i,j,k,l,ll;
    int err,n,s,nc,ncc,nic,iflag,lcc,lic,nrec;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Cycles of permutations. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,4,1))     /* get parameters */
        goto PCYCFin;

    printf1("Number of positions: %d\n",PMNV);
    if (PMF1Def == 0) { 
        printf1("Error: need an output file.\n");
        goto PCYCFin;
    }

    if (alloc_acn(PMNV + 1))  
        goto PCYCFin;
    if (alloc_ack(PMNV + 1))  
        goto PCYCFin;

    if (alloc_acr(PMNV + 1))  
        goto PCYCFin;
    if (alloc_acs(PMNV * PMNV + 1))  
        goto PCYCFin;

    nrec = 0;
    for (i = 0; i < NOC; ++i) {

        iflag = n = 0;
        for (j = 1; j <= PMNV; ++j)  
            AcK[j] = 0;

        for (j = 1; j <= PMNV; ++j) {
            s = (int)get_data(PMVIdx[j - 1],i);
            if (s >= 1 && s <= PMNV) {
                n++;
                if (AcK[s]) {
                    iflag = 1;
                    break;
                }
                AcN[j] = s;
                AcK[s] = 1;
            }
            else
                AcN[j] = 0;
        }
        if (iflag) {
            printf1("Inconsistent data in case %d.\n",i + 1);
            continue;
        }
        nc = ncc = nic = 0;

        for (j = 1; j <= PMNV; ++j)  
            AcK[j] = 0;

        iflag = 0;                      /* begin with incomplete cycles */
        while (iflag == 0) {
            iflag = 1;
            for (j = 1; j <= PMNV; ++j) {
                if (AcK[j])
                    continue;

                if (AcN[j] == 0) {
                    iflag = 0;
                    ll = 0;
                    AcK[j] = 1;
                    AcS[nc * PMNV + ll] = j;
                    ll++;
                    while (1) {
                        k = 0;
                        for (l = 1; l <= PMNV; ++l) {
                            if (AcK[l] == 0 && AcN[l] == j) {
                                k = l;
                                break;
                            }
                        }
                        if (k == 0)    
                            break;
                        
                        AcK[k] = 1;
                        AcS[nc * PMNV + ll] = k;
                        ll++;
                        j = k;
                    }
                    AcR[nc] = ll;
                    nic++;
                    nc++;
                    break;
                }
            }
        }
        iflag = 0;                      /* continue with complete cycles */
        while (iflag == 0) {           
            iflag = 1;
            for (j = 1; j <= PMNV; ++j) {
                if (AcK[j])
                    continue;
                iflag = 0;
                ll = 0;
                AcK[j] = 1;
                AcS[nc * PMNV + ll] = j;
                ll++;
                while (1) {
                    k = AcN[j];
                    if (AcK[k])    
                        break;
                    AcK[k] = 1;
                    AcS[nc * PMNV + ll] = k;
                    ll++;
                    j = k;
                }
                AcR[nc] = ll;
                ncc++;
                nc++;
                break;
            }
        }

        /* change order of incomplete cycles */

        for (j = 0; j < nic; ++j) {
            ll = AcR[j];
            for (l = 0; l < ll; ++l)
                AcK[l] = AcS[j * PMNV + l];
            for (l = 0; l < ll; ++l)
                AcS[j * PMNV + l] = AcK[ll - l - 1];
        }

        /* calculate length of cycles */

        lcc = lic = 0;
        for (j = 0; j < nic; ++j)   
            lic += (AcR[j] - 1);
        for (j = nic; j < nc; ++j)   
            lcc += (AcR[j] - 1);
    
        fprintf(PMF1d,PMNFmtS,i + 1);
        fprintf(PMF1d,PMNFmtS,n);
        fprintf(PMF1d,PMNFmtS,ncc);
        fprintf(PMF1d,PMNFmtS,nic);
        fprintf(PMF1d,PMNFmtS,lcc);
        fprintf(PMF1d,PMNFmtS,lic);

        for (j = 1; j <= PMNV; ++j)
            fprintf(PMF1d,PMNFmtS,AcN[j]);

        fprintf(PMF1d,"  ");
        for (j = nic; j < nc; ++j) {
            fprintf(PMF1d,"(%d",AcS[j * PMNV]);
            ll = AcR[j];
            for (l = 1; l < ll; ++l)
                fprintf(PMF1d,",%d",AcS[j * PMNV + l]);
            fprintf(PMF1d,")");
        }
        for (j = 0; j < nic; ++j) {
            fprintf(PMF1d,"[%d",AcS[j * PMNV]);
            ll = AcR[j];
            for (l = 1; l < ll; ++l)
                fprintf(PMF1d,",%d",AcS[j * PMNV + l]);
            fprintf(PMF1d,"]");
        }
        fprintf(PMF1d,"\n");
        nrec++;
    }
    printf1("%d records written to: %s\n",nrec,PMF1dName);
    err = 0;

PCYCFin:
    p_clean();      
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

int indep(void)
{
    register int i,j;
    int err,ix,iy,iw,x,y,xmin,xmax,ymin,ymax,nx,ny,np,ip; 
    double w,wsum,delta,fx,fy,fysum;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Delta-independence of two variables. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,4,1))     /* get parameters */
        goto INDEPFin;

    if (PMNV != 2) {
        printf1("Error: need exactly two variables on right-hand side.\n");
        goto INDEPFin;
    }
    ix = PMVIdx[0];
    iy = PMVIdx[1];
    iw = PMWVar;

    xmin = xmax = (int)get_data(ix,0);
    ymin = ymax = (int)get_data(iy,0);

    for (i = 1; i < NOC; ++i) {
        x = (int)get_data(ix,i);
        y = (int)get_data(iy,i);
        xmin = imin(xmin,x);
        xmax = imax(xmax,x);
        ymin = imin(ymin,y);
        ymax = imax(ymax,y);
    }
    printf1("Range of X: %d to %d\n",xmin,xmax);
    printf1("Range of Y: %d to %d\n",ymin,ymax);

    if (xmin < 0 || ymin < 0) {
        printf1("Error: values of X and Y must not be negative.\n");
        goto INDEPFin;
    }
    nx = xmax - xmin + 1;
    ny = ymax - ymin + 1;

    /* create partition of Y's property space */ 

    if (alloc_ack(ny + 1))  
        goto INDEPFin;
                 
    if (PMNTP > 0) {
        np = 0;
        j = ymin;
        for (i = 0; i < PMNTP; ++i) {
            y = (int)PMTP[i];
            if (y < ymin)
                continue;
            np++;
            while (j <= y) {
                if (j > ny)
                    break;
                AcK[j++] = np;
            }
            if (j > ny)
                break;
        }
    }
    else {
        np = 0;
        for (j = ymin; j <= ymax; ++j)
            AcK[j] = ++np;
    }
    printf1("Partition of Y's property space into %d subsets\n",np);
    for (i = 0; i <= ny; ++i)
        printf("%d ",AcK[i]);
    newline();

    for (i = 1; i <= np; ++i) {
        printf1("%5d : ",i);
        for (j = ymin; j <= ymax; ++j) {
            if (AcK[j] == i)
                printf("%d ",j);
        }
        newline();
    }



    if (alloc_acx(nx * ny + 1))  
        goto INDEPFin;
    if (alloc_acy(ny + 1))  
        goto INDEPFin;

    wsum = 0.0;
    w = 1.0;
    for (i = 0; i < NOC; ++i) {
        x = (int)get_data(ix,i);
        y = (int)get_data(iy,i);

        if (iw >= 0) {
            w = get_data(iw,i);
            if (w < 0.0) {
                printf1("Error: weights must not be negative.\n");
                goto INDEPFin;
            }
        }
        AcX[(x - 1) * nx + y] += w;
        AcY[y] += w;
        wsum += w;
    }
    printf1("Sum of frequencies: %g\n",wsum);
    if (wsum <= EPSI1) {
        err = 0;
        goto INDEPFin;
    }

    for (i = xmin; i <= xmax; ++i) {
        printf("i=%3d ",i);
        for (j = ymin; j <= ymax; ++j) {
            printf("%4d ",(int)AcX[(i - 1) * nx + j]);
        }
        printf("\n");
    }
    for (j = ymin; j <= ymax; ++j) {
        printf("%4d ",(int)AcY[j]);
    }
    printf("\n");


    /* do separately for all subsets of partition of Y's property space */ 

    for (ip = 1; ip <= np; ++ip) {

        printf1("Y-subset %d : ",ip);
        for (j = ymin; j <= ymax; ++j) {
            if (AcK[j] == ip)
                printf("%d ",j);
        }
        newline();

        /* begin with one-element X-sets */
       
        delta = 0.0;

        for (i = xmin; i <= xmax; ++i) {
            fx = 0.0;
            for (j = ymin; j <= ymax; ++j) 
                fx += AcX[(i - 1) * nx + j];
            fx /= wsum;

            fy = fysum = 0.0;

            for (j = ymin; j <= ymax; ++j) {
                if (AcK[j] != ip)
                    continue;

                fy += AcX[(i - 1) * nx + j];
                fysum += AcY[j];
            }
            if (fysum == 0.0)
                continue;
            fy /= fysum;

            delta = dmax(delta,fabs(fy - fx));

        }
        printf1("Maximal delta for 1-element sets: %g\n",delta);

    }



    err = 0;

INDEPFin:
    p_clean();      
    return(err);
}



