#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP C_tda_run(SEXP, SEXP, SEXP);
extern SEXP C_tda_unzoo(SEXP);
extern SEXP C_tda_zoo(SEXP, SEXP, SEXP, SEXP);

/* void (*)(void) is the generic function type for -Wcast-function-type */
#define TDA_DL(f) ((DL_FUNC)(void (*)(void))(f))

static const R_CallMethodDef CallEntries[] = {
    {"C_tda_run", TDA_DL(&C_tda_run), 3},
    {"C_tda_unzoo", TDA_DL(&C_tda_unzoo), 1},
    {"C_tda_zoo", TDA_DL(&C_tda_zoo), 4},
    {NULL, NULL, 0}
};
void R_init_tdaR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
