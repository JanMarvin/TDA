/* tda_interrupt.c -- user-interrupt support for the R package build.
 *
 * Long branch-and-bound runs (igmin, ivreg_min) poll
 * tda_check_interrupt() from their iteration loops.  The check uses
 * R_ToplevelExec around R_CheckUserInterrupt -- Urbanek's idiom -- so a
 * pending interrupt is DETECTED without the longjmp that
 * R_CheckUserInterrupt itself would perform.  The loops then unwind
 * through TDA's ordinary error paths, freeing everything, and the R
 * wrapper re-raises a proper interrupt condition.
 *
 * The actual check is throttled: R_ToplevelExec is not free, so it runs
 * every 1024th call; once an interrupt is seen the flag is sticky until
 * tda_interrupt_reset() at the start of the next run.
 *
 * Package-only: the standalone binary gets Ctrl-C from the OS.
 *
 * API status (checked against WRE for R 4.6): both entry points are
 * documented API -- R_CheckUserInterrupt in the "The R API" chapter's
 * interrupt section, R_ToplevelExec in "Condition handling and cleanup
 * code".  WRE steers new code toward R_UnwindProtect and friends for
 * most uses of R_ToplevelExec, but detecting a pending interrupt
 * WITHOUT unwinding is exactly the case those don't cover, and this
 * idiom is in long-standing use across CRAN packages.  If R core ever
 * demotes R_ToplevelExec, this file is the only place to touch.
 */

#include <Rinternals.h>
#include <R_ext/Utils.h>

static int tda_intr_seen = 0;
static unsigned tda_intr_tick = 0;

static void tda_intr_probe(void *dummy)
{
    (void)dummy;
    R_CheckUserInterrupt();
}

void tda_interrupt_reset(void)
{
    tda_intr_seen = 0;
    tda_intr_tick = 0;
}

int tda_interrupt_seen(void)
{
    return tda_intr_seen;
}

int tda_check_interrupt(void)
{
    if (tda_intr_seen)
        return 1;
    if ((++tda_intr_tick & 1023u) != 0u)
        return 0;
    if (R_ToplevelExec(tda_intr_probe, NULL) == FALSE)
        tda_intr_seen = 1;
    return tda_intr_seen;
}
