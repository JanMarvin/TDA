#include "tda.h"
#include "tda_context.h"
#include "t_gen.h"

/*  xshow opened a plot window through t_xwin.c, an X11 front end that
    is NOT part of this source tree and never was: the only definition
    of x_show has always been this stub.  The call site was wrapped in
    "#if S_XWIN", which tda.h defaulted to 1 and tda_const.h to 0 --
    hence the "S_XWIN redefined" warning on every file -- and the
    guarded line called x_show() with no argument while the stub takes
    a context, so it could not have compiled had anyone enabled it.

    The conditional is gone.  The command still exists, because TDA's
    own help documents it, and now says why it does nothing rather than
    returning silently.  */

int x_show(TDAContext *ctx)
{
    printf1(ctx, "Error: xshow needs the X11 front end, which is not part "
                 "of this build.\n");
    return(-1);
}
