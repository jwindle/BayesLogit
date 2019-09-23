// Jesse Windle, 2019
// Code modified from <https://github.com/jwindle/BayesLogit>

#include "R.h"
#include "Rmath.h"

// Define MACROS

#ifndef SIMPLE_RNG_WRAPPER
#define SIMPLE_RNG_WRAPPER

#define RCHECK 1000

const double SQRT2PI = 2.50662827;

inline void check_R_interupt(int count)
{
    #ifdef USE_R
    if (count % RCHECK == 0) R_CheckUserInterrupt();
    #endif
}

// CDF
#define p_norm(X, USE_LOG) pnorm((X), 0.0, 1.0, true, (USE_LOG))
#define p_gamma_scale(X, SHAPE, SCALE, USE_LOG) pgamma((X), (SHAPE), (SCALE), 1., (USE_LOG))
#define p_gamma_rate(X, SHAPE, RATE, USE_LOG) pgamma((X), (SHAPE), (1./RATE), 1., (USE_LOG))

// Random variates
#define expon_mean(MEAN) rexp((MEAN))
#define expon_rate(RATE) rexp((1./RATE))
#define unif() runif(0.0, 1.0)
#define norm(MEAN, SD) rnorm((MEAN), (SD))
#define gamma_scale(SHAPE, SCALE) rgamma((SHAPE), (SCALE))
#define gamma_rate(SHAPE, RATE) rgamma((SHAPE), (1./RATE))

// Scientific functions
#define lgamma(X) lgammafn((X))

#endif
