#ifndef TDA_RHOOKS_H
#define TDA_RHOOKS_H

#include <stdio.h>
#include <stdarg.h>
#include "tda_attr.h"

/* TDA keeps some state in file-scope variables rather than in TDAContext.
   Running tda_main() more than once in a process -- which the R binding does,
   but the command line program never did -- would carry that state over, so
   each owning file exposes a reset that puts it back to its start-up value.
   Every one of these is zero-initialised at start-up, so zeroing reproduces a
   fresh process exactly. */

void tda_reset_t_gen(void);
void tda_reset_t_plot(void);
void tda_reset_t_tri(void);
void tda_reset_t_lp(void);
void tda_reset_t_spl(void);
void tda_reset_t_zoo(void);
void tda_reset_globals(void);

#ifdef TDA_R_PACKAGE

/* exit() has 141 call sites.  Under R the process must survive, so it becomes
   a longjmp back to the .Call entry point.  It never returns -- it either
   longjmps or raises an R error -- and saying so keeps the compiler from
   treating every `exit()` as a path that carries on into the code below it,
   which is what made gd_proj_gety() look as though it could return an
   uninitialised y. */
#if defined(__GNUC__) || defined(__clang__)
void tda_r_exit(int status) __attribute__((noreturn));
#else
void tda_r_exit(int status);
#endif

#ifdef exit
#undef exit
#endif
#define exit(x) tda_r_exit(x)

/* R will not have a package writing to the process's own stdout or stderr:
   the connection has to be R's, so that output lands where R is sending it
   and the check tools find no printf, puts, stdout or stderr in the shared
   object.

   Redefining the C names themselves does not work -- the macros reach the
   declarations in <stdio.h> and glibc's inline functions -- so TDA gets its
   own three primitives instead.  They behave exactly like printf, vprintf
   and fprintf(stderr, ...), and outside the package that is what they are. */

#endif /* TDA_R_PACKAGE */

/* Console output.  tda_out() is printf, tda_vout() is vprintf and tda_err()
   is fprintf(stderr, ...); under R all three go through R's connection.
   TDA_CONSOLE is what it passes around when a FILE * argument means "the
   screen", and what `fd == TDA_CONSOLE` tests against. */
void tda_out(const char *fmt, ...) TDA_PRINTF(1, 2);
void tda_vout(const char *fmt, va_list ap) TDA_PRINTF(1, 0);
void tda_err(const char *fmt, ...) TDA_PRINTF(1, 2);
void tda_vout_err(const char *fmt, va_list ap) TDA_PRINTF(1, 0);
void tda_err_flush(void);
void tda_out_flush(void);
void tda_console_reset(void);
void tda_cf_set(const char *text);
const char *tda_console_text(int which, size_t *len);
/* long tda_out_line(void);  -- unused, body #if 0'd in tda_out.c */
long tda_out_line_start(void);
void tda_out_line_reset(void);

#ifdef TDA_R_PACKAGE
FILE *tda_console(void);

#define TDA_CONSOLE tda_console()
#else
#define TDA_CONSOLE stdout
#endif

#endif /* TDA_RHOOKS_H */
