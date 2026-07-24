/****************************************************************************/
/*  t_gcmd                                                                  */
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
#include "t_eval.h"  
#include "t_gdat.h"  
#include "t_gf.h"  
#include "t_var.h"  

/*  functions in t_gcmd.c */

int mem_use(void);
int prntime(void);
int prn_txt(void);
int machp(int opt);
int machp1(double start, int base);
double machf(int n);
double dsum(double a,double b);
int machcomp (double a, double b);
int t_parse(void);          
int data(int opt);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

int MacroN = 0;             /* number of macros                             */
char *MacroName[MaxMacro];  /* names of macros                              */
char *MacroDef[MaxMacro];   /* definition of macros                         */
char *MacroArgs[MaxMacro];  /* macro arguments                              */
short MacroNP[MaxMacro];    /* number of parameters                         */
int MacroNLen = 0;          /* max length of macro names                    */
int MacroExLevel = 0;       /* execution level of macro, used for local     */
int CurMacroNum  = -1;      /* current macor number, used for exists()      */

/* ------------------------------------------------------------------------ */
/*  mem_use()       print current use of memory.                            */
/*                  return 0 if OK, -1 if error.                            */
 
int mem_use(void)
{
    /****
    if (check_cmd(1))
        return(-1);
    ***/ 
    printf1("Currently requested memory: %d (%d) bytes.\n",MemReq,MxMReq);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prntime()       print current time.                                     */
/*                  return 0 if OK, -1 if error.                            */
 
int prntime(void)
{
    /***
    if (check_cmd(1))
        return(-1);
    ***/
#if TIME_ON
    printf1("Current time: ");
    prn_time(stdout);
#endif
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prn_txt()       print text: print(...).                                 */
/*                  return 0 if OK, -1 if error.                            */
 
int prn_txt(void)
{
    register char *p,*q;
    
    p = CmdBuf + 5;
    q = skip_blev(p);
    if (*--q != ')') {
        p_err(-1,1);
        return(-1);
    }
    *q = '\0';
    printf1("%s\n",p + 1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  machp(opt)  if opt == 0 set basic parameters                            */
/*                 opt == 1 print to standard out                           */
/*                                                                          */
/*      This function calculates the following machine parameters           */
/*                                                                          */
/*      beta  - The base of the machine.                                    */
/*      t     - The number of (beta) digits in the mantissa.                */
/*      rnd   - Whether proper rounding  (RND = 1)  or chopping (rnd = 0)   */
/*              occurs in addition. This may not be a reliable guide to     */
/*              the way in which the machine perfoms its arithmetic.        */
/*      eps   - The smallest positive number such that                      */
/*              fl(1.0 - eps) < 1.0, where fl denotes the computed value.   */
/*      emin  - The minimum exponent before (gradual) underflow occurs.     */
/*      rmin  - The smallest normalized number for the machine given by     */
/*              base ** (emin - 1), where base is the floating point        */
/*              value of beta.                                              */
/*      emax  - The maximum exponent before overflow occurs.                */
/*      rmax  - The largest positive number for the machine given by        */
/*              base ** emax * (1 - EPS), where base is the floating        */
/*              point value of beta.                                        */
/*      sfmin - a safe minimum, such that 1/sfmin does not overflow.        */
/*      imax  - largest integer                                             */
/*                                                                          */
/*      The code is adapted from the FORTRAN function DLAMCH (CMACH) which  */
/*      is part of the PAPACK package. The preliminary version dates from:  */
/*      March 26, 1990                                                      */
/*      Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,   */
/*      Courant Institute, NAG Ltd., and Rice University                    */
/*                                                                          */
/*      The  routine  is  based  on  the  routine  ENVRON  by  Malcolm and  */
/*      incorporates suggestions by Gentleman and Marovich. See             */
/*      Malcolm M. A. (1972) Algorithms to reveal properties of             */
/*          floating-point arithmetic. Comms. of the ACM, 15, 949-951.      */
/*      Gentleman W. M. and Marovich S. B. (1974) More on algorithms        */
/*          that reveal properties of floating point arithmetic units.      */
/*          Comms. of the ACM, 17, 276-277.                                 */
/*                                                                          */
/*      The  computation  of  EPS  is based  on  a  routine,  PARANOIA by   */
/*      W. Kahan of the University of California at Berkeley.               */
/*                                                                          */
/*      If opt == 0 the following global variables are set:                 */
/*      EPSI = eps                                                          */
/*      EPSI1 = sqrt(EPSI)                                                  */
/*      EPSI2 = 1000.0 * EPSI                                               */
/*      DBLMAX = rmax                                                       */
/*      DBLMIN = sfmin                                                      */
/*      INTMAX = imax                                                       */

int machp(int opt)
{
    int iwarn = 0;
    int i,beta,emax,emin,t,gnmin,gpmin,ngnmin,ngpmin,ieee,lieee1,rnd;
    int imax,tmax,lexp,uexp,ltry,exbits,nbits,expsum,itmp;
    double eps,rmax,rmin,recbas,y,z,oldy,sfmin;
    double a,b,c,savec,f,one,base,small,two,zero,tmp,t1,t2;
    /****************
    double half,sixth,third;
    ***********************/
    zero = 0.0;
    one = 1.0;
    two = 2.0;

    if (opt) {
        if (check_cmd(1))
            return(-1);
        printf1("Some parameters of the current machine:\n");
    }

    /*  Compute beta, the base of the machine. */
    /*  First, compute  a = 2.0**m  with the  smallest positive integer m such
        that fl( a + 1.0 ) = a. */

    a = one;
    c = one;

    while (machcomp(c,one)) {
        a *= 2.0;
	c = dsum(dsum(a,one),-a);
    }

    /*  Now compute  b = 2.0**m  with the smallest positive integer m
        such that fl( a + b ) .gt. a. */

    b = one;
    c = dsum(a,b);

    while (machcomp(c,a)) {
        b *= 2.0;
        c = dsum(a,b);
    }

    /*  Now compute the base.  a and c  are neighbouring floating point
        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
        their difference is beta. Adding 0.25 to c is to ensure that it
        is truncated to beta and not ( beta - 1 ). */

    tmp = one / 4.0;
    savec = c;
    c = dsum(c,-a);
    beta = (int) (c + tmp);

    if (opt)  
        printf1("Base of the machine: %d\n",beta);

    /*  Compute the largest integer, assume it is a power of beta */

    y = 1.0;
    tmax = imax = (int)y;
    while (1) {
        y *= 2.0;
        imax = (int)y;
        /*********************************************
        printf1("y=%lf imax=%d tmax=%d\n",y,imax,tmax);
        if (imax == tmax)
            break;
        tmax = imax;
        *****************/
        if (machf(imax) != y)
            break;
    }
    while (1) {
        y -= 1.0;
        imax = (int)y;
        if (machf(imax) == y)
            break;
    }
    INTMAX = imax;

    if (opt)  
        printf1("Largest integer: %d\n",imax);

    /* Now determine whether rounding or chopping occurs,  by adding a
       bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a. */
 
    b = (double)beta;
    f = dsum(b / 2.0, -b / 100.0);
    c = dsum(f,a);
    if (machcomp(c,a))
        rnd = 1;
    else
        rnd = 0;

    f = dsum(b / 2.0,b / 100.0);
    c = dsum(f,a);
    if (rnd && machcomp(c,a))
        rnd = 0;

    if (opt) {
        if (rnd)  
            printf1("Rounding seems proper\n");
        else
            printf1("Rounding seems to be by chopping\n");
    }

    /* Try and decide whether rounding is done in the  IEEE  'round to
       nearest' style. B/2 is half a unit in the last place of the two
       numbers A and savec. Furthermore, A is even, i.e. has last  bit
       zero, and savec is odd. Thus adding B/2 to A should not  change
       A, but adding B/2 to savec should change savec. */

    t1 = dsum(b / 2.0,a);
    t2 = dsum(b / 2.0,savec);
    if (machcomp(t1,a) && t2 > savec && rnd) {
        lieee1 = 1;
        if (opt)
            printf1("Seems to be: IEEE 'round to nearest' style\n");
    }
    else {
        lieee1 = 0;
        if (opt)
            printf1("Seems to be: NOT IEEE 'round to nearest' style\n");
    }

    /* Now find  the  mantissa, t.  It should  be the  integer part of
       log to the base beta of a,  however it is safer to determine  t
       by powering.  So we find t as the smallest positive integer for
       which fl( beta**t + 1.0 ) = 1.0. */

    t = 0;
    a = one;
    c = one;
    while (machcomp(one,c)) {
        t++;
        a *= (double)beta;
        c = dsum(dsum(one,a),-a);
    }

    if (opt)  
        printf1("Bits used for the mantissa: %d\n",t);

    /*  Start to find eps, i.e. the smallest positive number such that
        fl(1.0 - eps) less than 1.0. */

    b = (double)beta;
    a = pow(b,-(double)(t - 1));
    eps = a;
 
    /*  Try some tricks to see whether or not this is the correct eps. */
    /******************************************************************** 
    b = two / 3.0;
    half = one / 2.0;
    sixth = dsum(b,-half);
    third = dsum(sixth,sixth);
    b = dsum(third,-half);
    b = dsum(b,sixth);
    b = fabs(b);
    if (b < eps)
        b = eps;
 
    eps = 1.0;
    while (eps > b && b > zero) {
        eps = b;
        c = dsum(half * eps,pow(two,5.0) * pow(eps,2.0));
        c = dsum(half,-c);
        b = dsum(half,c);
        c = dsum(half,-b);
        b = dsum(half,c);
    }
    if (a < eps)
        eps = a;
    ******************************************************************/

    if (opt)  
        printf1("Machine epsilon: %20.15e\n",eps);
    else {
        EPSI = eps;
        EPSI1 = sqrt(EPSI);
        EPSI2 = 1000.0 * EPSI;
    }

    if (rnd)  
        tmp = pow((double)beta,1.0 - (double)(t - 1)) / 2.0;
    else  
        tmp = pow((double)beta,1.0 - (double)(t - 1));

    if (opt)  
        printf1("Machine epsilon (rounded): %20.15e\n",tmp);

    /*  Now find emin, i.e. the minimum exponent before (gradual)
        underflow occurs.
                        
        Let a = + or - 1, and + or - (1 + base ** (-3)), with base =
        (double)beta. Keep dividing a by beta until (gradual) underflow 
        occurs. This is detected when we cannot recover the previous a. */
 
    base = one / (double)beta;
    small = one;
    for (i = 1; i <= 3; ++i) 
        small = dsum(small * base,zero);

    a = dsum(one,small);
    ngpmin = machp1(one,beta);
    ngnmin = machp1(-one,beta);
    gpmin  = machp1(a,beta);
    gnmin  = machp1(-a,beta);
    ieee = 0;
 
    if ((ngpmin == ngnmin) && (gpmin == gnmin)) {

        if (ngpmin == gpmin) {
            if (opt)
                printf1("Not a twos-complement machine, no gradual underflow.\n");
            emin = ngpmin;
        }
        else if ((gpmin - ngpmin) == 3) {
            if (opt)
                printf1("Not a twos-complement machine, with gradual underflow.\n");
            emin = ngpmin - 1 + t;
            ieee = 1;
        }
        else {  
            if (opt)
                printf1("No known machine.\n");
            emin = ngpmin;
            if (emin > gpmin)
                emin = gpmin;
            iwarn = 1;
        }
    }
    else if ((ngpmin == gpmin) && (ngnmin == gnmin)) {

        if (fabs((double)(ngpmin - ngnmin)) == 1.0) {
            if (opt)
                printf1("A twos-complement machine, no gradual underflow.\n");
            emin = ngpmin;
            if (emin < ngnmin)
                emin = ngnmin;
        }
        else {  
            if (opt)
                printf1("No known machine.\n");
            emin = ngpmin;
            if (emin > ngnmin)
                emin = ngnmin;
            iwarn = 1;
        }
    }
    else if ((fabs((double)(ngpmin - ngnmin)) == 1.0) && (gpmin == gnmin)) {
        itmp = ngpmin;
        if (itmp > ngnmin)
            itmp = ngnmin;
        if ((gpmin - itmp) == 3) {
            if (opt)  
                printf1("A twos-complement machine, with gradual underflow.\n");
                     
            itmp = ngpmin;
            if (itmp < ngnmin)
                itmp = ngnmin;
            emin = itmp - 1 + t;
        }
        else {  
            if (opt)  
                printf1("No known machine.\n");
                     
            emin = ngpmin;
            if (emin > ngnmin)
                emin = ngnmin;
            iwarn = 1;
        }
    }
    else {  
        if (opt)   
            printf1("No known machine.\n");

        emin = ngpmin;
        if (emin > ngnmin)
            emin = ngnmin;
        if (emin > gpmin)
            emin = gpmin;
        if (emin > gnmin)
            emin = gnmin;
        iwarn = 1;
    }
 
    /* Assume IEEE arithmetic if we found denormalised numbers above,
       or if arithmetic seems to round in the IEEE style, determined
       above. A true IEEE machine should have both things true;
       however, faulty machines may have one or the other. */

    if (ieee || lieee1) {
        ieee = 1;
        if (opt)
            printf1("Seems to be the IEEE standard.\n");
    }
    else {
        ieee = 0;
        if (opt)
            printf1("Seems NOT to be the IEEE standard.\n");
    }
    if (opt) {
        printf1("Minimum exponent before underflow: %d\n",emin);
        if (iwarn)   
            printf1("This value may be incorrect!\n");
    }

    /*  Compute rmin, i.e. the smallest normalized number for the machine
        given by base ** (emin - 1), where base = (double)beta. Calculation
        is done by successive division by beta to avoid underflow during
        the computation. */

    rmin = 1.0;
    for (i = 1; i <= 1 - emin; ++i)
        rmin = dsum(rmin * base,zero);

    if (opt)  
        printf1("Smallest normalized number: %20.15e\n",rmin);

    /*  Compute emax, i.e. the maximum exponent before overflow occurs;
        and rmax, i.e. the largest positive number for the machine given
        by base ** emax  *  (1 - eps), where base = (double) beta.

        First compute lexp and uexp, two powers of 2 that bound
        abs(emin). We then assume that emax + abs(emin) will sum
        approximately to the bound that is closest to abs(emin).
        (emax is the exponent of the required number rmax). */

    lexp = 1;
    exbits = 1;
    while (1) {
        ltry = lexp * 2;
        if (ltry > -emin)  
            break;
        lexp = ltry;
        exbits++;
    }
    if (lexp == -emin) 
        uexp = lexp;
    else {
        uexp = ltry;
        exbits++;
    }

    /*  Now -lexp is less than or equal to emin, and -uexp is greater
        than or equal to emin. exbits is the number of bits needed to
        store the exponent. */

    if ((uexp + emin) > (-lexp - emin)) 
        expsum = 2 * lexp;
    else
        expsum = 2 * uexp;

    /*  expsum is the exponent range, approximately equal to
        emax - emin + 1. */

    emax = expsum + emin - 1;
    nbits = 1 + exbits + t;
 
    /* nbits is the total number of bits needed to store a
       floating-point number. */
 
    if ((nbits % 2) == 1 && (beta == 2)) {
 
        /* Either there are an odd number of bits used to store a
           floating-point number, which is unlikely, or some bits are
           not used in the representation of numbers, which is possible,
           (e.g. Cray machines) or the mantissa has an implicit bit,
           (e.g. IEEE machines, Dec Vax machines), which is perhaps the
           most likely. We have to assume the last alternative.
           If this is true, then we need to reduce emax by one because
           there must be some way of representing zero in an implicit-bit
           system. On machines like Cray, we are reducing emax by one
           unnecessarily. */
 
        emax--;
    }
    if (ieee) {
 
        /* Assume we are on an IEEE machine which reserves one exponent
           for infinity and NaN. */
  
        emax--;
    }
    if (opt)  
        printf1("Maximum exponent before overflow: %d\n",emax);

    /*  Now create rmax, the largest machine number, which should
        be equal to (1.0 - beta ** (-t)) * beta ** emax. */
 
    /*  First compute 1.0 - beta ** (-t), being careful that the
        result is less than 1.0. */
 
    recbas = one / (double)beta;
    z = one;
    y = zero;
    for (i = 1; i <= t; ++i) {
        z = z * recbas;
        if (y < one)
            oldy = y;
        y = dsum(y,z);
    }
    if (y >= one)
        y = oldy;
 
    /*  Now multiply by beta ** emax to get rmax. */
 
    for (i = 1; i <= emax; ++i)
        y = dsum(y * (double)beta,zero);

    rmax = y;
    if (opt)  
        printf1("Largest normalized number: %20.15e\n",rmax);
    else
        DBLMAX = rmax;

    /*  Calculate a safe minimum, sfmin, such that 1 / sfmin does not
        overflow. */

    sfmin = rmin;
    small = one / rmax;

    if (small >= sfmin) {

        /* Use small plus a bit, to avoid the possibility of rounding
           causing overflow when computing  1 / sfmin. */
 
        sfmin = small * (one + eps);
    }
    if (opt)  
        printf1("Safe minimum: %20.15e\n",sfmin);
    else
        DBLMIN = sfmin;

    return(0);
}

/* ------------------------------------------------------------------------ */
/*  machp1 (start,base)                                                     */
/*      Service routine for machp, returns emin.                            */

int machp1(double start, int base)
{
    int i,emin;
    double a,b1,b2,c1,c2,d1,d2,one,rbase,zero;

    a = start;
    one = 1.0;
    rbase = one / (double)base;
    zero = 0.0;
    emin = 1;
    b1 = dsum(a * rbase,zero);
    c1 = a;
    c2 = a;
    d1 = a;
    d2 = a;

    while (machcomp(c1,a) && machcomp(c2,a) && machcomp(d1,a) && machcomp(d2,a)) {
        emin--;
        a = b1;
        b1 = dsum(a / (double)base,zero);
        c1 = dsum(b1 * (double)base,zero);
        d1 = zero;
        for (i = 1; i <= base; ++i)
            d1 += b1;
        b2 = dsum(a * rbase,zero);
        c2 = dsum(b2 / rbase,zero);
        d2 = zero;
        for (i = 1; i <= base; ++i)
            d2 += b2;
    }
    return(emin);
}

double machf(int n)
{
    return((double)n);
}

/* ------------------------------------------------------------------------ */
/*  double dsum(a,b)                                                        */
/*      dsum is intended to force  a and b to be stored prior to doing      */
/*      the addition of a and b.  For use in situations where  optimizers   */
/*      might hold one of these in a register.                              */

double dsum (double a, double b) 
{
    return(a + b);
}

/* ------------------------------------------------------------------------ */
/*  int  machcomp(a,b)                                                      */
/*      machcomp is intended to force  a and b to be compared explicitly.   */

int machcomp (double a, double b) 
{
  if (a==b)
    return(1);
  return(0);
}

/*--------------------------------------------------------------------------*/
/*  t_parse()   execute parse command.                                      */
/*              return 0 if OK, -1 if error.                                */

int t_parse(void)           
{
    int r;
    char *p;

    if (check_cmd(1))
        return(-1);

    p = CmdBuf + 5;
    if (*p++ != '=' || !*p) {
        p_err(-1,1);
        return(-1);
    }             
    MatExprFlg = 1;

    r = v_parse(p,1);
    if (r < 0 || ESCnt == 0) {
        printf1("Syntax error (%d).\n",r);
        if (r < 0)
            prn_emsg1(r);
    }
    else {
        parse(ESCnt,ESTyp,ESVal,ESIdx);
        printf1("Dimension: %d x %d\n",MatExprRow,MatExprCol);
    }
    MatExprFlg = 0;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  data(opt)   print info about current data, if opt != 0 show also        */
/*              list of variables.                                          */
 
int data(int opt)
{
    register int i,j,k;
    int n,l;

    /***
    if (check_cmd(1))
        return(-1);
    ***/
    if (DMDef == 0) {
        printf1("No data matrix defined.\n");
        return(0);
    }
    printf1("Number of cases in data matrix: %d\n",NOCDM);
    if (TSelFlg)
        printf1("Selected with tsel command: %d\n",NOC);
    printf1("Number of variables: %d\n",NVAR);
    printf1("Number of namelists: %d\n",NNL);
    if (opt) {
        printf1("\n");
        prn_var(VIFirst);

        if (NNL > 0) {
            printf1("\nNamelist");
            if (NNL > 1)
                printf1("s");
            newline();
            n = 1;
            for (i = 0; i < MaxNL; ++i) {
                if (NLNV[i] > 0) {
                    l = strlen(NLName[i]);
                    if (n < l)
                        n = l;
                }
            }
            for (i = 0; i < MaxNL; ++i) {
                j = NLNV[i];
                if (j > 0) {
                    printf1("%s",NLName[i]);
                    prnchar(' ',n - strlen(NLName[i]),0);
                    printf1(" = %s",VName[NLVIdx[i][0]]);
                    for (k = 1; k < j; ++k)
                        printf1(",%s",VName[NLVIdx[i][k]]);
                    printf1("\n");
                }
            }
        }
    }
    return(0);
}






