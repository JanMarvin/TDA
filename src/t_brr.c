/****************************************************************************/
/*  t_brr                                                                   */
/*                                                                          */
/*  TDA. Program for Transition Data Analysis, written by Goetz Rohwer.     */
/*  Copyright (C) 1989,1991-94 Goetz Rohwer. All rights reserved.           */
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
#include "t_gf.h"
#include "t_alloc.h"
#include "t_parm.h"

#define BRRTEST 0

/*  functions in t_brr.c */

#ifndef LINT
int brr(void);
int brr_gen(int ns, int nu,int *nrep,int opt);
int brr_mod(int a, int p);
int brr_pp(int t, int n, int *a);
#endif

/* ------------------------------------------------------------------------ */
/*  Pointer to starting element of polynomial used for arithmetic in        */
/*  GP(p**n), i.e. pointers to the fields BHHMPC, also to BHHGPC            */

short BHHP[13][14] = {
/*   0   1   2   3   4   5   6   7   8   9  10  11  12  13 */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1,  0, 88, -1,121, -1,139, -1, -1, -1,157, -1,169,
    -1, -1,  3, 91, -1,124, -1,142, -1, -1, -1,160, -1,172,
    -1, -1,  7, 95, -1,128, -1,146, -1, -1, -1, -1, -1, -1,
    -1, -1, 12,100, -1,133, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, 18,106, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, 25,113, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, 33, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, 42, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, 52, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, 63, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 
};

/*  Coefficients of polynomials used for arithmetic in GP(p**n) */

#define MaxPDEG 12

char BHHMPC[] = {
    /* p =  2  n =  2 */  1,1,1,                        /*  0       */
    /* p =  2  n =  3 */  1,1,0,1,                      /*  3       */
    /* p =  2  n =  4 */  1,1,0,0,1,                    /*  7       */
    /* p =  2  n =  5 */  1,0,1,0,0,1,                  /*  12      */
    /* p =  2  n =  6 */  1,1,0,0,0,0,1,                /*  18      */
    /* p =  2  n =  7 */  1,1,0,0,0,0,0,1,              /*  25      */
    /* p =  2  n =  8 */  1,1,0,1,1,0,0,0,1,            /*  33      */
    /* p =  2  n =  9 */  1,1,0,0,0,0,0,0,0,1,          /*  42      */
    /* p =  2  n = 10 */  1,0,0,1,0,0,0,0,0,0,1,        /*  52      */
    /* p =  2  n = 11 */  1,0,1,0,0,0,0,0,0,0,0,1,      /*  63      */
    /* p =  2  n = 12 */  1,0,0,1,0,0,0,0,0,0,0,0,1,    /*  75      */

    /* p =  3  n =  2 */  1,0,1,                        /*  88      */
    /* p =  3  n =  3 */  1,2,0,1,                      /*  91      */
    /* p =  3  n =  4 */  2,1,0,0,1,                    /*  95      */
    /* p =  3  n =  5 */  2,1,0,1,0,1,                  /*  100     */
    /* p =  3  n =  6 */  2,1,0,0,0,0,1,                /*  106     */
    /* p =  3  n =  7 */  1,0,2,0,0,0,0,1,              /*  113     */

    /* p =  5  n =  2 */  2,0,1,                        /*  121     */
    /* p =  5  n =  3 */  4,1,1,1,                      /*  124     */
    /* p =  5  n =  4 */  2,0,0,0,1,                    /*  128     */
    /* p =  5  n =  5 */  1,4,0,0,0,1,                  /*  133     */

    /* p =  7  n =  2 */  1,0,1,                        /*  139     */
    /* p =  7  n =  3 */  2,0,0,1,                      /*  142     */
    /* p =  7  n =  4 */  1,6,0,0,1,                    /*  146     */
    /* p =  7  n =  5 */  5,5,2,0,0,1,                  /*  151     */

    /* p = 11  n =  2 */  1,0,1,                        /*  157     */
    /* p = 11  n =  3 */  5,8,0,1,                      /*  160     */
    /* p = 11  n =  4 */  6,9,0,0,1,                    /*  164     */

    /* p = 13  n =  2 */  2,0,1,                        /*  169     */
    /* p = 13  n =  3 */  2,0,0,1                       /*  172     */
};

/*  Coefficients of generating polynomials for GP(p**n) */
 
char BHHGPC[] = {
    /* p =  2  n =  2 */  0,0,1,
    /* p =  2  n =  3 */  0,1,0,0,
    /* p =  2  n =  4 */  0,1,0,0,0,
    /* p =  2  n =  5 */  0,1,0,0,0,0,
    /* p =  2  n =  6 */  0,1,0,0,0,0,0,
    /* p =  2  n =  7 */  0,1,0,0,0,0,0,0,
    /* p =  2  n =  8 */  1,1,0,0,0,0,0,0,0,
    /* p =  2  n =  9 */  1,1,1,0,0,0,0,0,0,0,
    /* p =  2  n = 10 */  0,1,0,0,0,0,0,0,0,0,0,
    /* p =  2  n = 11 */  0,1,0,0,0,0,0,0,0,0,0,0,
    /* p =  2  n = 12 */  0,0,0,0,0,0,0,0,0,0,0,0,0,

    /* p =  3  n =  2 */  1,1,0,
    /* p =  3  n =  3 */  0,1,0,0,
    /* p =  3  n =  4 */  0,1,0,0,0,
    /* p =  3  n =  5 */  0,2,0,0,0,0,
    /* p =  3  n =  6 */  0,1,0,0,0,0,0,
    /* p =  3  n =  7 */  0,1,0,0,0,0,0,0,

    /* p =  5  n =  2 */  1,1,0,
    /* p =  5  n =  3 */  1,1,0,0,
    /* p =  5  n =  4 */  1,1,0,0,0,
    /* p =  5  n =  5 */  0,2,0,0,0,0,

    /* p =  7  n =  2 */  2,1,0,
    /* p =  7  n =  3 */  1,3,0,0,
    /* p =  7  n =  4 */  2,1,0,0,0,
    /* p =  7  n =  5 */  0,0,0,0,0,0,

    /* p = 11  n =  2 */  4,1,0,
    /* p = 11  n =  3 */  0,1,0,0,
    /* p = 11  n =  4 */  0,0,0,0,0,

    /* p = 13  n =  2 */  2,1,0,
    /* p = 13  n =  3 */  2,1,0,0
};

/* ------------------------------------------------------------------------ */
/*  brr         Calculate number of replications for BRR method.            */
/*                                                                          */
/*              brr = ns,nu;                                                */
/*                                                                          */
/*              ns = number of strate                                       */
/*              nu = number of secu                                         */
/*                                                                          */
/*              Return 0 if OK, otherwise -1.                               */

int brr(void)
{
    int err,ns,nu,nr,sflg,opt;

    opt = 0;
    err = -1;
    if (check_cmd(1))
        return(-1);

    if (parm(CmdBuf + 3,5,1))     /* get parameters */
        goto BRR1Fin;

    if (PMTabFDef)
        opt = 2;

    ns = (int)PMRHSA;
    nu = (int)PMRHSB;
    sflg = SILENTFlg;
    SILENTFlg = -1;
    err = brr_gen(ns,nu,&nr,opt);
    SILENTFlg = sflg;
    if (err == 0 && opt)
        printf1("Matrix written to: %s\n",PMTabFName);

BRR1Fin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  brr_gen(ns,nu,opt)  Generates orthogonal weights for ns strata and      */
/*                      nu Secu's.                                          */
/*                                                                          */
/*  The number of rows in AcC (= number of replications) is returned        */
/*  in nrep.                                                                */
/*                                                                          */
/*  If opt != 0, the function creates nr x ns matrix AcC containing         */  
/*  references for the secus. If opt == 2 the matrix is written into        */
/*  PMTabFd.                                                                */
/*                                                                          */
/*  Return 0 if OK, -1 if error (makes local error messages)                */

int brr_gen(int ns, int nu,int *nrep,int opt)
{
    register int i,j,k;
    int n,m,np,nr,na,nb,nc,nd,err;
    int a[MaxPDEG + 2],b[MaxPDEG + 2];
    int c[2 * MaxPDEG + 4],d[2 * MaxPDEG + 4];
    double x;

    err = -1;
    *nrep = 0;
    if (ns < 2) { 
        printf1("Error: number of strata less than 2.\n");
        return(-1);
    } 
    if (nu < 2 || nu > 13) {
        printf1("Error: secu=%d not possible.\n",nu);
        return(-1);
    } 
    for (i = 2; i < nu; ++i) {  /* check for prime */ 
        n = nu / i;
        if (n * i == nu) {
            printf1("Error: number of secu must be prime.\n");
            return(-1);
        }
    }

    /*  calculate np, modulus for polynomial ring */

    np = (int) ceil( log((double)(ns * (nu - 1) + 1)) / log((double)nu) -
                                                                      1.e-12);

    /*  calculate nr, the number of replications */

    nr = (int)pow((double)nu,(double)np);

    printf2("Creating orthogonal weights.\n");
    printf2("Power of Galois field: %d\n",np);
    printf2("Number of replications: %d\n",nr);
         
    if (opt) {
        if (alloc_acc(nr * ns + 1))
            return(-1);
    }      

    /*  get polynomial for arithmetic in GP(p**n) in b[], the generating
        polynomial in a[] */

    if (np < 2 || np > 12 || nu < 2 || nu > 13 || (j = BHHP[np][nu]) < 0) {
        printf1("Error: can't generate orthogonal weights for strata=%d, secu=%d\n",
                                                                    ns,nu);
        goto BRRFin;
    }
    for (k = 0; k <= np; ++k) {
        a[k] = BHHGPC[j];
        b[k] = BHHMPC[j];
        j++;
    }
    na = nb = np;
    while (na && a[na] == 0)
        na--;

    printf2("Generating polynomial:     ");
    for (k = 0; k <= na; ++k)
        printf2("%d ",a[k]);
    printf2("\nPolynomial for arithmetic: ");
    for (k = 0; k <= nb; ++k)
        printf2("%d ",b[k]);
    printf2("\n");
          
    if (na == 0) {
        printf1("Error in generating polynomial.\n");
        goto BRRFin;
    }

    /* calculate leading column in AcC: the polynomial c[], at beginning
       equal to a[], is successively multiplied by a[] */

    if (opt)
        AcC[1] = (char)a[0];
    for (j = 0; j <= na; ++j)
        c[j] = a[j];
    nc = na;

#if BRRTEST
        brr_pp(3,nc,c);
#endif

    for (i = 2; i < nr; ++i) {

        /* multiplication of c[] with a[], result in d[] */

#if BRRTEST
            printf1("\nreplication: %d\n",i);
            printf1("multiplicate:             ");
            brr_pp(3,nc,c);
            printf1("with:                     ");
            brr_pp(1,na,a);
#endif    
        nd = na + nc;
        for (j = 0; j <= nd; ++j)
            d[j] = 0;

        for (j = 0; j <= na; ++j) {
            for (k = 0; k <= nc; ++k)  
                d[j + k] += a[j] * c[k];
        }
#if BRRTEST
            printf1("result:                   ");
            brr_pp(4,nd,d);
#endif      

        /* change d[] to modulus nu */

        for (j = 0; j <= nd; ++j)
            d[j] = brr_mod(d[j],nu);
        while (nd && d[nd] == 0)
            nd--;

#if BRRTEST
            printf1("result after modulus p:   ");
            brr_pp(4,nd,d);
#endif     

        /* change d[] to modulus of b[] */

        while (nd >= nb) {
            n = nd - nb;
            x = (double)d[nd] / (double)b[nb];
            m = (int)x;
            if ((double)m != x) {
                printfe("FATAL ERROR IN BRR (m=%d x=%lf)\n",m,x);
                gerr_exit(224);
            }
            for (j = n; j <= nd; ++j)
                d[j] -= m * b[j - n];
            while (nd && d[nd] == 0)
                nd--;
        }
#if BRRTEST
            printf1("result after modulus b[]: ");
            brr_pp(4,nd,d);
#endif     

        /* again change d[] with modulus nu */

        for (j = 0; j <= nd; ++j)  
            d[j] = brr_mod(d[j],nu);
           
        while (nd && d[nd] == 0)
            nd--;

#if BRRTEST
            printf1("result after modulus p:   ");
            brr_pp(4,nd,d);
#endif     

        /* copy d[] to c[] */

        nc = nd;
        for (j = 0; j <= nc; ++j)
            c[j] = d[j];


        /* save zero order coefficient of c[] in AcC */

        if (opt)
            AcC[(i - 1) * ns + 1] = (char)c[0];

    }
    if (nc != 0 || c[0] != 1) {
        printf1("Error: no success in calculating orthogonal weights.\n");
        goto BRRFin;
    }
            
    /* create matrix of orthogonal weights */

    if (opt) {
        for (j = 2; j <= ns; ++j) {
            for (i = 2; i < nr; ++i)
                AcC[(i - 2) * ns + j] = AcC[(i - 1) * ns + j - 1];
            AcC[(nr - 2) * ns + j] = AcC[j - 1];
        }
        if (opt == 2) {         /* write into PMTabFd */

            for (i = 1; i <= ns; ++i) {
                for (j = 0; j < nr; ++j)
                    fprintf(PMTabFd,"%d ",(int)AcC[j * ns + i]);
                fprintf(PMTabFd,"\n");
            }
        }
    }
    err = 0;

BRRFin:
    if (err) {
        if (opt)
            alloc_acc(0);
    }
    else
        *nrep = nr;
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  brr_mod(a,p)    return a modulus p                                      */

int brr_mod(int a, int p) 
{
    int m;
    if (a > 0)
        m = a % p;
    else if (a < 0) {
        m = p - ((-a) % p);
        if (m >= p)
            m = m % p;
    }
    else
        m = 0;
    return(m);
}

/*--------------------------------------------------------------------------*/
/*  brr_pp(t,n,a)   Print a polynomial, only for testing.                   */

#if BRRTEST
int brr_pp(int t, int n, int *a) 
{
    register int j;
    if (t == 1)
        printf1("a: ");
    else if (t == 2)
        printf1("b: ");
    else if (t == 3)
        printf1("c: ");
    else if (t == 4)
        printf1("d: ");

    for (j = 0; j <= n; ++j)
        printf1("%2d ",a[j]);
    printf1("\n");
}
#endif

