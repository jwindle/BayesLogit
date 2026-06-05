#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP C_rpg_hybrid(SEXP, SEXP, SEXP);
extern SEXP C_rpg_hybrid_fill(SEXP, SEXP);
extern SEXP C_rpg_gamma(SEXP, SEXP, SEXP);
extern SEXP C_rpg_devroye(SEXP, SEXP);
extern SEXP C_rpg_sp(SEXP, SEXP);

static const R_CallMethodDef callMethods[] = {
    {"C_rpg_hybrid",      (DL_FUNC) &C_rpg_hybrid,      3},
    {"C_rpg_hybrid_fill", (DL_FUNC) &C_rpg_hybrid_fill, 2},
    {"C_rpg_gamma",       (DL_FUNC) &C_rpg_gamma,       3},
    {"C_rpg_devroye",     (DL_FUNC) &C_rpg_devroye,     2},
    {"C_rpg_sp",          (DL_FUNC) &C_rpg_sp,          2},
    {NULL, NULL, 0}
};

void R_init_BayesLogitCallableTest(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
