// Jesse Windle, 2019
// Code modified from <https://github.com/jwindle/BayesLogit>

#include <stdio.h>
#include <stdexcept>
#include <cmath>

#ifndef TRUNCATED_GAMMA
#define TRUNCATED_GAMMA

double right_tgamma_reject(double shape, double rate);
double omega_k(int k, double a, double b);
double right_tgamma_beta(double shape, double rate);
double rtgamma_rate(double shape, double rate, double right_t);
double ltgamma(double shape, double rate, double trunc);

#endif
