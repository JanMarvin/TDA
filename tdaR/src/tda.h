/****************************************************************************/
/* tda.h                                                                   */
/*           Definitions to include in all TDA modules.                     */
/*                                                                          */

#ifndef _TDA_H
#define _TDA_H

#define TDA_Version 6.4

/* Marks a switch case that deliberately runs into the next one, for
   -Wimplicit-fallthrough.  GCC also accepts a /-* FALLTHROUGH *-/
   comment, clang does not; the attribute works on both.  The [[ ]]
   spelling is only used where the standard has it, since -Wpedantic
   rejects it before C23. */
#if defined(__has_attribute)
#  if __has_attribute(fallthrough)
#    define TDA_FALLTHROUGH __attribute__((fallthrough))
#  endif
#elif defined(__GNUC__) && __GNUC__ >= 7
#  define TDA_FALLTHROUGH __attribute__((fallthrough))
#endif
#if !defined(TDA_FALLTHROUGH) && defined(__STDC_VERSION__) && \
    __STDC_VERSION__ >= 202311L
#  define TDA_FALLTHROUGH [[fallthrough]]
#endif
#ifndef TDA_FALLTHROUGH
#  define TDA_FALLTHROUGH ((void)0)
#endif

/*  set following flags:
    for DOS set: S_DOS to 1, S_UNIX to 0
    for UNIX set: S_UNIX to 1, S_DOS 0
*/

#ifndef S_DOS
#define S_DOS       0
#endif
#ifndef S_UNIX
#define S_UNIX      1
#endif



/* set ARCHTyp to 1 for 32 bit IEEE with standard byte order, e.g.,
   SPARC, RS/6000.

   set ARCHTyp to 2 for 32 bit IEEE with reversed byte order, e.g.,
   most INTEL processors. In particular, when compiling TDA for DOS-
   like platforms, ARCHTyp should normally be set to 2.
*/

#define  ARCHTyp    2

/* ALPHA64 selects a 32-bit TLONG.  It was written for the DEC Alpha, where
   long is 64 bits and int 32 -- which is exactly the LP64 model every current
   Linux and macOS build uses, so it has to be on there too.  TLONG carries the
   4-byte fields of the Zoo archive format: with a 64-bit TLONG, to_long()
   reassembles a negative 32-bit value as a large positive one and the
   zoo_start + zoo_minus == 0 check in get_zoo() fails for every archive. */
#if defined(__LP64__) || defined(_LP64) || defined(__alpha)
#define ALPHA64     1
#else
#define ALPHA64     0
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#define TIME_ON     1       /* set TIME_ON to 0 if problems with time.h */
#if TIME_ON
#include <time.h>
#include "tda_const.h"
#include "tda_context.h"
#endif

#if ALPHA64                 /* for 64 architecture  */
#define TLONG int
#else
#define TLONG long
#endif

#if S_UNIX
#define  OPEN_RD    "r"
#define  OPEN_WR    "w"
#define  OPEN_AP    "a"
/* "b" is a no-op on POSIX and required on Windows, which the portable build
   reaches through S_UNIX: without it the C runtime translates CRLF inside
   binary files and the Zoo archive, shapefile, SPSS and dBase readers get
   corrupted input.
   OPEN_WR is deliberately left as text mode.  That means TDA's own output
   files end their lines CRLF on Windows and LF elsewhere, so anything that
   reports their size or dumps their bytes differs by one byte per line
   between platforms -- which is why dsplit, dump and rplz are excluded from
   the generated tests.  Changing OPEN_WR to "wb" would make those portable
   and every output file byte-identical across platforms, at the cost of
   giving Windows users LF-only files. */
#define  OPEN_RB    "rb"
#define  OPEN_WB    "wb+"
#define  OPEN_RU    "r+"
#define  PSepC      '/'     /* path name separator */
#define  PSepS      "/"
#define  PathSepC   ':'
#endif

#if S_DOS
#define  OPEN_RD    "rt"
#define  OPEN_WR    "wt"
#define  OPEN_AP    "at"
#define  OPEN_RB    "rb"
#define  OPEN_WB    "wb+"
#define  OPEN_RU    "rt+"
#define  PSepC      '\\'    /* path name separator */
#define  PSepS      "\\"
#define  PathSepC   ';'
#endif

/* ------------------------------------------------------------------------ */

#define Pi  3.141592653589793

#define CR  0x0d            /* Carriage return                              */
#define LF  0x0a            /* Line feed                                    */
#define ExpMin  -350.0      /* Max and Min arguments for rexp-function      */
#define ExpMax   350.0      /* -350, 350                                    */
#define LogMin 10e-153      /* Min value for rlog function                  */

#define CmdBufLen   20000   /* length of command buffer                     */
#define FNMaxLen      120   /* max file name length                         */
#define LLMaxLen      460   /* max length of interactive lines              */
#define RLMaxDef    20000   /* default max record length for input files    */
#define OBufLen     21000   /* max length of output lines                   */
                            /* note: OBufLen >= RLMaxDef required           */
#define MaxIFLEV      100   /* max nesting of if/endif                      */
#define MaxREP        100   /* max nesting of repeat and while              */
#define MaxREPCMD    1000   /* max commands in each repeat and while level  */

#define P_ALLOC     10000   /* stack length for parser                      */
#define LLEN           76   /* length of separator lines                    */
#define VNLMax         32   /* max length of variable names                 */
#define VLLMax        120   /* max length of variable labels                */
#define MatDefLen      78   /* max length of matrix definition              */
#define SFMTLEN        16   /* max length of format strings                 */
#define MaxCTDim       26   /* max number of contingency table dimensions   */
#define MaxCTCat     1000   /* max number of categories per dimension       */
#define MaxLL          50   /* max number of loglinear models per section   */
#define MaxP         1000   /* max number of model parameters               */
#define MaxW          200   /* max number of panel waves                    */
#define MaxYC         100   /* Max number of categories in dep. variables   */
#define MaxMDC         10   /* Max number of mdcon commands                 */
#define PLMaxLen      120   /* max length of plot labels                    */
#define MaxMatchV       5   /* max number of matching variable (pairs)      */
#define MaxNVRec      100   /* max number of vrec commands                  */
#define MaxXARG      1000   /* max number of x arguments for atab classes   */
#define MaxFNP         50   /* max temp parameters in function definition   */
#define MaxDGRP        10   /* max number of dgrp groups                    */

#define MaxMatDef     200   /* def. max number of matrices                  */
#define MaxMatFunc    100   /* max number of matrix functions               */
#define MaxNVDef     2000   /* def. max number of variables                 */
#define NOCDef       1000   /* default (max) number of cases in data matrix */
#define MaxSTDDef      20   /* default max stack size for derivatives       */
#define MaxNL          50   /* max number of namelists                      */

/* ------------------------------------------------------------------------ */
#define MaxSPL       1000   /* Max number of split variables                */
#define MaxTranDef    100   /* Default max number of transitions            */
#define MaxNT         250   /* Max number of time periods (<= 256)          */
#define MaxVP         250   /* Max number of vector definitions             */
#define MaxVT         500   /* Max number of variables per vector           */
#define MaxG          100   /* Max number of groups                         */
#define MaxSt         255   /* Max number of O- and D- states               */
#define MaxQ          200   /* Max number of quantile definitions           */
#define DMVMax         20   /* Max number of dimensions for dmv()           */

#define MaxNodeN    32000   /* highest node number                          */
#define RDP_VMax     1000   /* max number of variables in RDP_VIdx          */
#define DS_VMax      1000   /* max number of variables in DS_VIdx           */
#define SEQ_Max       100   /* max number of seqence data structures        */
#define PMPSMax        20   /* max number of test patterns                  */
#define MFMAX          50   /* max number of merge files                    */
#define MaxMacro      100   /* max number of macros                         */

/* ------------------------------------------------------------------------ */
/*  Model numbers                                                           */

#define MCOX       1        /* Cox Model                                    */
#define MEXP       2        /* Exponential Model                            */
#define MEXP1      3        /* Exponential Model with Time-periods          */
#define MPOL       4        /* Model with Polynomial Rates (I)              */
#define MPOL1      5        /* Model with Polynomial Rates (II)             */
#define MGM        6        /* Gompertz-Makeham Model                       */
#define MWEI       7        /* Weibull Model                                */
#define MSIC       8        /* Sickle Model                                 */
#define MLL        9        /* Log-logistic Model (I)                       */
#define MLL2      10        /* Log-logistic Model (II)                      */
#define MLLA      11        /* Log-logistic Model (IIa)                     */
#define MLN       12        /* Log-normal Model (I)                         */
#define MGAM      13        /* Gamma Model                                  */
#define MIG       14        /* Inverse Gaussian Model                       */
#define MEXP2     16        /* Exponential Model with Time-periods, type II */

#define DLR       20        /* Discrete Time Logistic Regression            */
#define CLL       21        /* Discrete Time Complementary Log-Log          */

#define QRLOG1    41        /* Binary logit                                 */
#define QRLOG2    42        /* Ordered logit (kumulative logits)            */
#define QRLOG3    43        /* Multinomial logit                            */
#define QRLOG4    44        /* Conditional logit                            */

#define QRPROB1   51        /* Binary probit                                */
#define QRPROB2   52        /* Ordered probit                               */
#define QRPROB3   53        /* Multivariate probit                          */
#define QRPROB4   54        /* Sim. binary probits                          */

#define QRLOG1A   101       /* simple logit model (test of fnml())          */
#define QRCLOG    102       /* complementary log-log                        */

#define RMOD1     111       /* Rasch model 1                                */

#define CLMOD1    121       /* LS fit of ultrametric distances              */

#define MDSMOD1   201       /* multi-dimensional scaling                    */
#define MDSMOD2   202       /* multi-dimensional scaling                    */
#define MDSMOD3   203       /* multi-dimensional scaling                    */
#define MDSMOD4   204       /* multi-dimensional scaling                    */
#define MDSMOD5   205       /* multi-dimensional scaling                    */

#define IVRMOD1   221       /* interval regression                          */
#define IVRMOD2   222       /* interval regression                          */
#define IVRMOD3   223       /* interval regression                          */

#define UMOD      999       /* user-defined function                        */

/* ------------------------------------------------------------------------ */

/* Under the package build every fprintf() to a context-held stream is
   routed through tda_fprintf() (tda_export.c), which writes exactly the
   same bytes and additionally hands the numbers to R.  ctx is in scope
   at all 3247 call sites -- checked, not assumed -- so the macro needs
   no argument beyond what the calls already pass.  Files that must
   reach the real fprintf define TDA_NO_FPRINTF_TAP before including
   this header. */
#if defined(TDA_R_PACKAGE) && !defined(TDA_NO_FPRINTF_TAP)
#define fprintf(fd, ...) tda_fprintf(ctx, fd, __VA_ARGS__)
#endif

#endif /* _TDA_H */

/* end of tda.h */
