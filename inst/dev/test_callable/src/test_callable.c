/* Exercises the BayesLogit C-callable API via R_GetCCallable, exactly as a
 * downstream package would use it via LinkingTo: BayesLogit and BayesLogit.h.
 */

#include <R.h>
#include <Rinternals.h>
#include <BayesLogit.h>   /* from inst/include/, resolved via LinkingTo */

/* rpg_hybrid — scalar draw */
SEXP C_rpg_hybrid(SEXP h_, SEXP z_, SEXP n_)
{
    BayesLogit_rpg_hybrid_t fn = BayesLogit_rpg_hybrid();
    if (!fn) error("R_GetCCallable returned NULL for rpg_hybrid");

    int n = asInteger(n_);
    SEXP out = PROTECT(allocVector(REALSXP, n));
    GetRNGstate();
    for (int i = 0; i < n; i++)
        REAL(out)[i] = fn(asReal(h_), asReal(z_));
    PutRNGstate();
    UNPROTECT(1);
    return out;
}

/* rpg_hybrid_fill — vectorised draw */
SEXP C_rpg_hybrid_fill(SEXP h_, SEXP z_)
{
    BayesLogit_rpg_hybrid_fill_t fn = BayesLogit_rpg_hybrid_fill();
    if (!fn) error("R_GetCCallable returned NULL for rpg_hybrid_fill");

    int n = length(h_);
    SEXP out = PROTECT(allocVector(REALSXP, n));
    GetRNGstate();
    fn(n, REAL(h_), REAL(z_), REAL(out));
    PutRNGstate();
    UNPROTECT(1);
    return out;
}

/* rpg_gamma — scalar draw */
SEXP C_rpg_gamma(SEXP h_, SEXP z_, SEXP trunc_)
{
    BayesLogit_rpg_gamma_t fn = BayesLogit_rpg_gamma();
    if (!fn) error("R_GetCCallable returned NULL for rpg_gamma");

    GetRNGstate();
    double val = fn(asReal(h_), asReal(z_), asInteger(trunc_));
    PutRNGstate();
    return ScalarReal(val);
}

/* rpg_devroye — scalar draw */
SEXP C_rpg_devroye(SEXP h_, SEXP z_)
{
    BayesLogit_rpg_devroye_t fn = BayesLogit_rpg_devroye();
    if (!fn) error("R_GetCCallable returned NULL for rpg_devroye");

    GetRNGstate();
    double val = fn(asInteger(h_), asReal(z_));
    PutRNGstate();
    return ScalarReal(val);
}

/* rpg_sp — scalar draw */
SEXP C_rpg_sp(SEXP h_, SEXP z_)
{
    BayesLogit_rpg_sp_t fn = BayesLogit_rpg_sp();
    if (!fn) error("R_GetCCallable returned NULL for rpg_sp");

    int iter = 0;
    GetRNGstate();
    double val = fn(asReal(h_), asReal(z_), &iter);
    PutRNGstate();
    return ScalarReal(val);
}
