// -*- mode: c++; c-basic-offset: 4 -*-
// (C) Nicholas Polson, James Scott, Jesse Windle, 2012-2019
// Contributor: Jari Oksanen

// This file is part of BayesLogit.

// BayesLogit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.

// BayesLogit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with
// BayesLogit.  If not, see <https://www.gnu.org/licenses/>.


#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> /* for NULL */
#include <R_ext/Rdynload.h>
// #include "polyagamma_wrapper.h"

/* Skeletons of declarations of public C functions */

// extern void EM(void *, void *, void *, void *, void *, void *, void *, void *);
extern void rpg_alt(void *, void *, void *, void *);
extern void rpg_devroye(void *, void *, void *, void *);
extern void rpg_gamma(void *, void *, void *, void *, void *);
extern void rpg_hybrid(void *, void *, void *, void *);
extern void rpg_sp(void *, void *, void *, void *, void*);

extern double BayesLogit_rpg_gamma(double, double, int);
extern void BayesLogit_rpg_gamma_fill(int, const double *, const double *,
                                      int, double *);
extern double BayesLogit_rpg_devroye(int, double);
extern void BayesLogit_rpg_devroye_fill(int, const int *, const double *,
                                        double *);
extern double BayesLogit_rpg_sp(double, double, int *);
extern void BayesLogit_rpg_sp_fill(int, const double *, const double *,
                                   double *, int *);
extern double BayesLogit_rpg_hybrid(double, double);
extern void BayesLogit_rpg_hybrid_fill(int, const double *, const double *,
                                       double *);

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

void R_init_BayesLogit(DllInfo *dll)
{
    R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);

    R_RegisterCCallable("BayesLogit", "rpg_gamma",
                        (DL_FUNC) &BayesLogit_rpg_gamma);
    R_RegisterCCallable("BayesLogit", "rpg_gamma_fill",
                        (DL_FUNC) &BayesLogit_rpg_gamma_fill);
    R_RegisterCCallable("BayesLogit", "rpg_devroye",
                        (DL_FUNC) &BayesLogit_rpg_devroye);
    R_RegisterCCallable("BayesLogit", "rpg_devroye_fill",
                        (DL_FUNC) &BayesLogit_rpg_devroye_fill);
    R_RegisterCCallable("BayesLogit", "rpg_sp",
                        (DL_FUNC) &BayesLogit_rpg_sp);
    R_RegisterCCallable("BayesLogit", "rpg_sp_fill",
                        (DL_FUNC) &BayesLogit_rpg_sp_fill);
    R_RegisterCCallable("BayesLogit", "rpg_hybrid",
                        (DL_FUNC) &BayesLogit_rpg_hybrid);
    R_RegisterCCallable("BayesLogit", "rpg_hybrid_fill",
                        (DL_FUNC) &BayesLogit_rpg_hybrid_fill);
}
