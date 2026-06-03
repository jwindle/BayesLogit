// -*- mode: c++; c-basic-offset: 4 -*-
// Cross-package C-callable API for BayesLogit Polya-Gamma RNGs.

#ifndef BAYESLOGIT_API_H
#define BAYESLOGIT_API_H

#include <R_ext/Rdynload.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * These functions use R's RNG but do not call GetRNGstate() or PutRNGstate().
 * Callers must own R's RNG state before using them, following the same pattern
 * as native R RNG calls such as unif_rand(), rnorm(), or rgamma().
 */

typedef double (*BayesLogit_rpg_gamma_t)(double h, double z, int trunc);
typedef void (*BayesLogit_rpg_gamma_fill_t)(int n, const double *h,
                                            const double *z, int trunc,
                                            double *out);

typedef double (*BayesLogit_rpg_devroye_t)(int h, double z);
typedef void (*BayesLogit_rpg_devroye_fill_t)(int n, const int *h,
                                              const double *z, double *out);

typedef double (*BayesLogit_rpg_sp_t)(double h, double z, int *iter);
typedef void (*BayesLogit_rpg_sp_fill_t)(int n, const double *h,
                                         const double *z, double *out,
                                         int *iter);

typedef double (*BayesLogit_rpg_hybrid_t)(double h, double z);
typedef void (*BayesLogit_rpg_hybrid_fill_t)(int n, const double *h,
                                             const double *z, double *out);

static inline BayesLogit_rpg_gamma_t BayesLogit_rpg_gamma(void)
{
    return (BayesLogit_rpg_gamma_t)
        R_GetCCallable("BayesLogit", "rpg_gamma");
}

static inline BayesLogit_rpg_gamma_fill_t BayesLogit_rpg_gamma_fill(void)
{
    return (BayesLogit_rpg_gamma_fill_t)
        R_GetCCallable("BayesLogit", "rpg_gamma_fill");
}

static inline BayesLogit_rpg_devroye_t BayesLogit_rpg_devroye(void)
{
    return (BayesLogit_rpg_devroye_t)
        R_GetCCallable("BayesLogit", "rpg_devroye");
}

static inline BayesLogit_rpg_devroye_fill_t BayesLogit_rpg_devroye_fill(void)
{
    return (BayesLogit_rpg_devroye_fill_t)
        R_GetCCallable("BayesLogit", "rpg_devroye_fill");
}

static inline BayesLogit_rpg_sp_t BayesLogit_rpg_sp(void)
{
    return (BayesLogit_rpg_sp_t)
        R_GetCCallable("BayesLogit", "rpg_sp");
}

static inline BayesLogit_rpg_sp_fill_t BayesLogit_rpg_sp_fill(void)
{
    return (BayesLogit_rpg_sp_fill_t)
        R_GetCCallable("BayesLogit", "rpg_sp_fill");
}

static inline BayesLogit_rpg_hybrid_t BayesLogit_rpg_hybrid(void)
{
    return (BayesLogit_rpg_hybrid_t)
        R_GetCCallable("BayesLogit", "rpg_hybrid");
}

static inline BayesLogit_rpg_hybrid_fill_t BayesLogit_rpg_hybrid_fill(void)
{
    return (BayesLogit_rpg_hybrid_fill_t)
        R_GetCCallable("BayesLogit", "rpg_hybrid_fill");
}

#ifdef __cplusplus
}
#endif

#endif
