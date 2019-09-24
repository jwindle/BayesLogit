// Jari Oksanen, 2019

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> /* for NULL */
#include <R_ext/Rdynload.h>
#include "polyagamma_wrapper.h"

/* Skeletons of declarations of public C functions */

/* // extern void EM(void *, void *, void *, void *, void *, void *, void *, void *); */
/* extern void rpg_alt(void *, void *, void *, void *); */
/* extern void rpg_devroye(void *, void *, void *, void *); */
/* extern void rpg_gamma(void *, void *, void *, void *, void *); */
/* extern void rpg_hybrid(void *, void *, void *, void *); */
/* extern void rpg_sp(void *, void *, void *, void *, void*); */

/* table of C calls */

/* NB. The Makevars of 0.6 release exclude AR1.o and DynExpFamMH.o
 * from the DLL, and the functions in these are therefore commented
 * out from the following CMethodDef list (they are undefined). */

static const R_CMethodDef cMethods[] = {
  //    {"EM",                       (DL_FUNC) &EM,                       8},
    {"rpg_alt",                  (DL_FUNC) &rpg_alt,                  4},
    {"rpg_devroye",              (DL_FUNC) &rpg_devroye,              4},
    {"rpg_gamma",                (DL_FUNC) &rpg_gamma,                5},
    {"rpg_hybrid",               (DL_FUNC) &rpg_hybrid,               4},
    {"rpg_sp",                   (DL_FUNC) &rpg_sp,                   5},
    {NULL, NULL, 0}
};

/* initialize registration */

void R_init_BayesLogit2(DllInfo *dll)
{
    R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
