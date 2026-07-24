#ifndef TDA_CONST_H
#define TDA_CONST_H

/****************************************************************************/
/* tda.h                                                                   */
/*           Definitions to include in all TDA modules.                     */
/*                                                                          */

#define TDA_Version 6.4 

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

#define FNMaxLen      120   /* max file name length                         */
#define MaxIFLEV      100   /* max nesting of if/endif                      */
#define MaxREP        100   /* max nesting of repeat and while              */
#define VNLMax         32   /* max length of variable names                 */
#define VPFmtSLen      10   /* size of ctx->VPFmtS[] print format strings   */
#define NVPFmtSLen     10   /* size of ctx->NVPFmtS[] print format strings  */
#define MaxLL          50   /* max number of loglinear models per section   */
#define MaxP         1000   /* max number of model parameters               */
#define MaxMatchV       5   /* max number of matching variable (pairs)      */
#define MaxFNP         50   /* max temp parameters in function definition   */
#define MaxDGRP        10   /* max number of dgrp groups                    */
#define MaxNL          50   /* max number of namelists                      */
#define MaxSPL       1000   /* Max number of split variables                */
#define MaxG          100   /* Max number of groups                         */
#define DMVMax         20   /* Max number of dimensions for dmv()           */
#define SEQ_Max       100   /* max number of seqence data structures        */
#define PMPSMax        20   /* max number of test patterns                  */
#define MFMAX          50   /* max number of merge files                    */
#define MaxMacro      100   /* max number of macros                         */

/* Module-private constants hoisted for struct sizing */
#ifndef MEXOPMax
#  define MEXOPMax   100
#endif
#ifndef MXN
#  define MXN          7
#endif
#ifndef MXMaxLen
#  define MXMaxLen    60
#endif
#ifndef PMFNMax
#  define PMFNMax    100
#endif
#ifndef E00BufL
#  define E00BufL   2000
#endif
#ifndef DFILMax
#  define DFILMax    100
#endif
#ifndef FE_BufLen
#  define FE_BufLen  132
#endif
#ifndef XLSNFMax
#  define XLSNFMax    50
#endif
#ifndef ESMAXNF
#  define ESMAXNF    100
#endif
#ifndef SIZ_TEXT
#  define SIZ_TEXT    20
#  define PATHSIZE   256
#  define LFNAMESIZE 256
#endif

#endif /* TDA_CONST_H */
