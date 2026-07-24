#include "tda_rhooks.h"

void tda_reset_globals(void)
{
    tda_reset_t_gen();
    tda_reset_t_plot();
    tda_reset_t_tri();
    tda_reset_t_lp();
    tda_reset_t_spl();
    tda_reset_t_zoo();
}
