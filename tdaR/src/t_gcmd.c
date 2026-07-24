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
#include "tda_context.h"
#include <float.h>
#include <limits.h>

/*  functions in t_gcmd.c */

int mem_use(TDAContext *ctx);
int prntime(TDAContext *ctx);
int prn_txt(TDAContext *ctx);
int machp(TDAContext *ctx, int opt);
int machp1(TDAContext *ctx, double start, int base);
double machf(TDAContext *ctx, int n);
double dsum(TDAContext *ctx, double a,double b);
void mach_ftz_note(TDAContext *ctx);
int machcomp(TDAContext *ctx, double a, double b);
int t_parse(TDAContext *ctx);
int data(TDAContext *ctx, int opt);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */


/* ------------------------------------------------------------------------ */
/*  mem_use()       print current use of memory.                            */
/*                  return 0 if OK, -1 if error.                            */

int mem_use(TDAContext *ctx)
{
    /****
    if (check_cmd(ctx, 1))
        return(-1);
    ***/
    printf1(ctx, "Currently requested memory: %d (%d) bytes.\n",ctx->MemReq,ctx->MxMReq);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prntime()       print current time.                                     */
/*                  return 0 if OK, -1 if error.                            */

int prntime(TDAContext *ctx)
{
    /***
    if (check_cmd(ctx, 1))
        return(-1);
    ***/
#if TIME_ON
    printf1(ctx, "Current time: ");
    prn_time(ctx, TDA_CONSOLE);
#endif
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prn_txt()       print text: print(...).                                 */
/*                  return 0 if OK, -1 if error.                            */

int prn_txt(TDAContext *ctx)
{
    register char *p,*q;

    p = ctx->CmdBuf + 5;
    q = skip_blev(ctx, p);
    if (*--q != ')') {
        p_err(ctx, -1,1);
        return(-1);
    }
    *q = '\0';
    printf1(ctx, "%s\n",p + 1);
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

int machp(TDAContext *ctx, int opt)
{
    int iwarn = 0;
    int i = 0,beta = 0,emax = 0,emin = 0,t = 0,gnmin = 0,gpmin = 0,ngnmin = 0,ngpmin = 0,ieee = 0,lieee1 = 0,rnd = 0;
    int imax = 0,lexp = 0,uexp = 0,ltry = 0,exbits = 0,nbits = 0,expsum = 0,itmp = 0;
    double eps = 0.0,rmax = 0.0,rmin = 0.0,recbas = 0.0,y = 0.0,z = 0.0,oldy = 0.0,sfmin = 0.0;
    double a = 0.0,b = 0.0,c = 0.0,savec = 0.0,f = 0.0,one = 0.0,base = 0.0,small = 0.0,zero = 0.0,tmp = 0.0,t1 = 0.0,t2 = 0.0;
    /****************
    double half = 0.0,sixth = 0.0,third = 0.0;
    ***********************/
    zero = 0.0;
    one = 1.0;

    if (opt) {
        if (check_cmd(ctx, 1))
            return(-1);
        printf1(ctx, "Some parameters of the current machine:\n");
    }

    /* On any IEEE 754 machine (all modern hardware) the values machp
       detects through floating-point experiments are exactly known from
       the C standard headers.  Use them directly to avoid undefined
       behaviour in the detection loops (float->int cast overflow,
       etc.) that trap on ARM64 / Apple Silicon with strict clang. */
    if (!opt) {
        ctx->EPSI   = DBL_EPSILON;
        ctx->EPSI1  = DBL_EPSILON < 1e-10 ? 1.4901161193847656e-08  /* sqrt(DBL_EPSILON) */
                                           : 1.0;
        ctx->EPSI2  = 1000.0 * DBL_EPSILON;
        ctx->DBLMAX = DBL_MAX;
        ctx->DBLMIN = DBL_MIN;
        ctx->INTMAX = INT_MAX;
        return(0);
    }

    /*  Compute beta, the base of the machine. */
    /*  First, compute  a = 2.0**m  with the  smallest positive integer m such
        that fl( a + 1.0 ) = a. */

    a = one;
    c = one;

    while (machcomp(ctx, c,one)) {
        a *= 2.0;
	c = dsum(ctx, dsum(ctx, a,one),-a);
    }

    /*  Now compute  b = 2.0**m  with the smallest positive integer m
        such that fl( a + b ) .gt. a. */

    b = one;
    c = dsum(ctx, a,b);

    while (machcomp(ctx, c,a)) {
        b *= 2.0;
        c = dsum(ctx, a,b);
    }

    /*  Now compute the base.  a and c  are neighbouring floating point
        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
        their difference is beta. Adding 0.25 to c is to ensure that it
        is truncated to beta and not ( beta - 1 ). */

    tmp = one / 4.0;
    savec = c;
    c = dsum(ctx, c,-a);
    beta = (int) (c + tmp);

    if (opt)
        printf1(ctx, "Base of the machine: %d\n",beta);

    /*  Compute the largest integer, assume it is a power of beta */

    y = 1.0;
    imax = (int)y;
    while (1) {
        y *= 2.0;
        /* guard against UB: (int)y is undefined if y > INT_MAX.
           On ARM64/Apple clang this traps. Stop before overflow. */
        if (y > 2147483647.0) {
            imax = 2147483647;
            break;
        }
        imax = (int)y;
        if (machf(ctx, imax) != y)
            break;
    }
    while (1) {
        y -= 1.0;
        if (y < 0.0) { imax = 0; break; }
        imax = (int)y;
        if (machf(ctx, imax) == y)
            break;
    }
    ctx->INTMAX = imax;

    if (opt)
        printf1(ctx, "Largest integer: %d\n",imax);

    /* Now determine whether rounding or chopping occurs,  by adding a
       bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a. */

    b = (double)beta;
    f = dsum(ctx, b / 2.0, -b / 100.0);
    c = dsum(ctx, f,a);
    if (machcomp(ctx, c,a))
        rnd = 1;
    else
        rnd = 0;

    f = dsum(ctx, b / 2.0,b / 100.0);
    c = dsum(ctx, f,a);
    if (rnd && machcomp(ctx, c,a))
        rnd = 0;

    if (opt) {
        if (rnd)
            printf1(ctx, "Rounding seems proper\n");
        else
            printf1(ctx, "Rounding seems to be by chopping\n");
    }

    /* Try and decide whether rounding is done in the  IEEE  'round to
       nearest' style. B/2 is half a unit in the last place of the two
       numbers A and savec. Furthermore, A is even, i.e. has last  bit
       zero, and savec is odd. Thus adding B/2 to A should not  change
       A, but adding B/2 to savec should change savec. */

    t1 = dsum(ctx, b / 2.0,a);
    t2 = dsum(ctx, b / 2.0,savec);
    if (machcomp(ctx, t1,a) && t2 > savec && rnd) {
        lieee1 = 1;
        if (opt)
            printf1(ctx, "Seems to be: IEEE 'round to nearest' style\n");
    }
    else {
        lieee1 = 0;
        if (opt)
            printf1(ctx, "Seems to be: NOT IEEE 'round to nearest' style\n");
    }

    /* Now find  the  mantissa, t.  It should  be the  integer part of
       log to the base beta of a,  however it is safer to determine  t
       by powering.  So we find t as the smallest positive integer for
       which fl( beta**t + 1.0 ) = 1.0. */

    t = 0;
    a = one;
    c = one;
    while (machcomp(ctx, one,c)) {
        t++;
        a *= (double)beta;
        c = dsum(ctx, dsum(ctx, one,a),-a);
    }

    if (opt)
        printf1(ctx, "Bits used for the mantissa: %d\n",t);

    /*  Start to find eps, i.e. the smallest positive number such that
        fl(1.0 - eps) less than 1.0. */

    b = (double)beta;
    a = pow(b,-(double)(t - 1));
    eps = a;

    /*  Try some tricks to see whether or not this is the correct eps. */
    /********************************************************************
    b = two / 3.0;
    half = one / 2.0;
    sixth = dsum(ctx, b,-half);
    third = dsum(ctx, sixth,sixth);
    b = dsum(ctx, third,-half);
    b = dsum(ctx, b,sixth);
    b = fabs(b);
    if (b < eps)
        b = eps;

    eps = 1.0;
    while (eps > b && b > zero) {
        eps = b;
        c = dsum(ctx, half * eps,pow(two,5.0) * pow(eps,2.0));
        c = dsum(ctx, half,-c);
        b = dsum(ctx, half,c);
        c = dsum(ctx, half,-b);
        b = dsum(ctx, half,c);
    }
    if (a < eps)
        eps = a;
    ******************************************************************/

    if (opt)
        printf1(ctx, "Machine epsilon: %20.15e\n",eps);
    else {
        ctx->EPSI = eps;
        ctx->EPSI1 = sqrt(ctx->EPSI);
        ctx->EPSI2 = 1000.0 * ctx->EPSI;
    }

    if (rnd)
        tmp = pow((double)beta,1.0 - (double)(t - 1)) / 2.0;
    else
        tmp = pow((double)beta,1.0 - (double)(t - 1));

    if (opt)
        printf1(ctx, "Machine epsilon (rounded): %20.15e\n",tmp);

    /*  Now find emin, i.e. the minimum exponent before (gradual)
        underflow occurs.

        Let a = + or - 1, and + or - (1 + base ** (-3)), with base =
        (double)beta. Keep dividing a by beta until (gradual) underflow
        occurs. This is detected when we cannot recover the previous a. */

    base = one / (double)beta;
    small = one;
    for (i = 1; i <= 3; ++i)
        small = dsum(ctx, small * base,zero);

    a = dsum(ctx, one,small);
    ngpmin = machp1(ctx, one,beta);
    ngnmin = machp1(ctx, -one,beta);
    gpmin  = machp1(ctx, a,beta);
    gnmin  = machp1(ctx, -a,beta);
    ieee = 0;

    if ((ngpmin == ngnmin) && (gpmin == gnmin)) {

        if (ngpmin == gpmin) {
            if (opt) {
                printf1(ctx, "Not a twos-complement machine, no gradual underflow.\n");
                mach_ftz_note(ctx);
            }
            emin = ngpmin;
        }
        else if ((gpmin - ngpmin) == 3) {
            if (opt)
                printf1(ctx, "Not a twos-complement machine, with gradual underflow.\n");
            emin = ngpmin - 1 + t;
            ieee = 1;
        }
        else {
            if (opt)
                printf1(ctx, "No known machine.\n");
            emin = ngpmin;
            if (emin > gpmin)
                emin = gpmin;
            iwarn = 1;
        }
    }
    else if ((ngpmin == gpmin) && (ngnmin == gnmin)) {

        if (fabs((double)(ngpmin - ngnmin)) == 1.0) {
            if (opt)
                printf1(ctx, "A twos-complement machine, no gradual underflow.\n");
            emin = ngpmin;
            if (emin < ngnmin)
                emin = ngnmin;
        }
        else {
            if (opt)
                printf1(ctx, "No known machine.\n");
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
                printf1(ctx, "A twos-complement machine, with gradual underflow.\n");

            itmp = ngpmin;
            if (itmp < ngnmin)
                itmp = ngnmin;
            emin = itmp - 1 + t;
        }
        else {
            if (opt)
                printf1(ctx, "No known machine.\n");

            emin = ngpmin;
            if (emin > ngnmin)
                emin = ngnmin;
            iwarn = 1;
        }
    }
    else {
        if (opt)
            printf1(ctx, "No known machine.\n");

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
            printf1(ctx, "Seems to be the IEEE standard.\n");
    }
    else {
        ieee = 0;
        if (opt)
            printf1(ctx, "Seems NOT to be the IEEE standard.\n");
    }
    if (opt) {
        printf1(ctx, "Minimum exponent before underflow: %d\n",emin);
        if (iwarn)
            printf1(ctx, "This value may be incorrect!\n");
    }

    /*  Compute rmin, i.e. the smallest normalized number for the machine
        given by base ** (emin - 1), where base = (double)beta. Calculation
        is done by successive division by beta to avoid underflow during
        the computation. */

    rmin = 1.0;
    for (i = 1; i <= 1 - emin; ++i)
        rmin = dsum(ctx, rmin * base,zero);

    if (opt)
        printf1(ctx, "Smallest normalized number: %20.15e\n",rmin);

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
        printf1(ctx, "Maximum exponent before overflow: %d\n",emax);

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
        y = dsum(ctx, y,z);
    }
    if (y >= one)
        y = oldy;

    /*  Now multiply by beta ** emax to get rmax. */

    for (i = 1; i <= emax; ++i)
        y = dsum(ctx, y * (double)beta,zero);

    rmax = y;
    if (opt)
        printf1(ctx, "Largest normalized number: %20.15e\n",rmax);
    else
        ctx->DBLMAX = rmax;

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
        printf1(ctx, "Safe minimum: %20.15e\n",sfmin);
    else
        ctx->DBLMIN = sfmin;

    return(0);
}

/* ------------------------------------------------------------------------ */
/*  machp1 (start,base)                                                     */
/*      Service routine for machp, returns emin.                            */

int machp1(TDAContext *ctx, double start, int base)
{
    int i,emin;
    double a,b1,b2,c1,c2,d1,d2,one,rbase,zero;

    a = start;
    one = 1.0;
    rbase = one / (double)base;
    zero = 0.0;
    emin = 1;
    b1 = dsum(ctx, a * rbase,zero);
    c1 = a;
    c2 = a;
    d1 = a;
    d2 = a;

    while (machcomp(ctx, c1,a) && machcomp(ctx, c2,a) && machcomp(ctx, d1,a) && machcomp(ctx, d2,a)) {
        emin--;
        a = b1;
        b1 = dsum(ctx, a / (double)base,zero);
        c1 = dsum(ctx, b1 * (double)base,zero);
        d1 = zero;
        for (i = 1; i <= base; ++i)
            d1 += b1;
        b2 = dsum(ctx, a * rbase,zero);
        c2 = dsum(ctx, b2 / rbase,zero);
        d2 = zero;
        for (i = 1; i <= base; ++i)
            d2 += b2;
    }
    return(emin);
}

double machf(TDAContext *ctx, int n)
{
    (void)ctx;        /* unused: the signature is shared */
    return((double)n);
}

/* ------------------------------------------------------------------------ */
/*  mach_ftz_note()                                                         */
/*                                                                          */
/*  Absence of gradual underflow on hardware that has it almost always      */
/*  means the build flushes denormals to zero, which -ffast-math does by    */
/*  linking crtfastmath.o -- it sets flush-to-zero in the FPU control       */
/*  register before main() runs.  The probe is reporting the truth about    */
/*  the runtime, but a user comparing results against another build needs   */
/*  to know why the two differ, so say it here rather than leave them to    */
/*  find it.                                                                */

void mach_ftz_note(TDAContext *ctx)
{
    volatile double tiny = 1.0e-308;
    volatile double d = tiny / 1.0e10;

    if (d == 0.0)
        printf1(ctx, "Note: this build flushes denormals to zero (as -ffast-math does).\n"
                     "      Results may differ slightly from a build that does not.\n");
}

/* ------------------------------------------------------------------------ */
/*  double dsum(a,b)                                                        */
/*      dsum is intended to force  a and b to be stored prior to doing      */
/*      the addition of a and b.  For use in situations where  optimizers   */
/*      might hold one of these in a register.                              */

double dsum(TDAContext *ctx, double a, double b)
{
    /*  The comment above says this forces a and b to be stored before
        the addition.  A plain "return a + b" does not: every modern
        compiler inlines it and folds the arithmetic, and the machine
        probe that depends on it then measures the compiler rather than
        the hardware.  volatile forces real stores and loads, which
        defeats constant folding, fused multiply-add contraction and
        reassociation alike.  This is the same device LAPACK uses in
        dlamc3, for the same reason.

        It does NOT rescue a -ffast-math build, and nothing here could:
        that flag links crtfastmath.o, which sets flush-to-zero in the
        FPU control register at startup, so denormals genuinely do not
        exist at runtime.  The probe then reports "no gradual
        underflow", and it is right to -- measured directly, a denormal
        computed under -ffast-math really does come back as 0.

        All thirty callers are in the machine-characteristics probe, so
        the cost is a handful of stores, once.  */
    volatile double x = a;
    volatile double y = b;
    volatile double s = x + y;

    (void)ctx;
    return(s);
}

/* ------------------------------------------------------------------------ */
/*  int  machcomp(a,b)                                                      */
/*      machcomp is intended to force  a and b to be compared explicitly.   */

int machcomp(TDAContext *ctx, double a, double b)
{
    (void)ctx;        /* unused: the signature is shared */
  if (a==b)
    return(1);
  return(0);
}

/*--------------------------------------------------------------------------*/
/*  t_parse()   execute parse command.                                      */
/*              return 0 if OK, -1 if error.                                */

int t_parse(TDAContext *ctx)
{
    int r;
    char *p;

    if (check_cmd(ctx, 1))
        return(-1);

    p = ctx->CmdBuf + 5;
    if (*p++ != '=' || !*p) {
        p_err(ctx, -1,1);
        return(-1);
    }
    ctx->MatExprFlg = 1;

    r = v_parse(ctx, p,1);
    if (r < 0 || ctx->ESCnt == 0) {
        printf1(ctx, "Syntax error (%d).\n",r);
        if (r < 0)
            prn_emsg1(ctx, r);
    }
    else {
        parse(ctx, ctx->ESCnt,ctx->ESTyp,ctx->ESVal,ctx->ESIdx);
        printf1(ctx, "Dimension: %d x %d\n",ctx->MatExprRow,ctx->MatExprCol);
    }
    ctx->MatExprFlg = 0;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  data(opt)   print info about current data, if opt != 0 show also        */
/*              list of variables.                                          */

int data(TDAContext *ctx, int opt)
{
    register int i,j,k;
    int n,l;

    /***
    if (check_cmd(ctx, 1))
        return(-1);
    ***/
    if (ctx->DMDef == 0) {
        printf1(ctx, "No data matrix defined.\n");
        return(0);
    }
    printf1(ctx, "Number of cases in data matrix: %d\n",ctx->NOCDM);
    if (ctx->TSelFlg)
        printf1(ctx, "Selected with tsel command: %d\n",ctx->NOC);
    printf1(ctx, "Number of variables: %d\n",ctx->NVAR);
    printf1(ctx, "Number of namelists: %d\n",ctx->NNL);
    if (opt) {
        printf1(ctx, "\n");
        prn_var(ctx, ctx->VIFirst);

        if (ctx->NNL > 0) {
            printf1(ctx, "\nNamelist");
            if (ctx->NNL > 1)
                printf1(ctx, "s");
            newline(ctx);
            n = 1;
            for (i = 0; i < MaxNL; ++i) {
                if (ctx->NLNV[i] > 0) {
                    l = (int)(strlen(ctx->NLName[i]));
                    if (n < l)
                        n = l;
                }
            }
            for (i = 0; i < MaxNL; ++i) {
                j = ctx->NLNV[i];
                if (j > 0) {
                    printf1(ctx, "%s",ctx->NLName[i]);
                    prnchar(ctx, ' ',n - (int)strlen(ctx->NLName[i]),0);
                    printf1(ctx, " = %s",ctx->VName[ctx->NLVIdx[i][0]]);
                    for (k = 1; k < j; ++k)
                        printf1(ctx, ",%s",ctx->VName[ctx->NLVIdx[i][k]]);
                    printf1(ctx, "\n");
                }
            }
        }
    }
    return(0);
}
