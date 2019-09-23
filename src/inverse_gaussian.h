// Jesse Windle, 2019
// Code modified from <https://github.com/jwindle/BayesLogit>

#include <stdio.h>
#include <stdexcept>
#include <cmath>

#ifndef INVERSE_GAUSSIAN
#define INVERSE_GAUSSIAN

double igauss(double mu, double lambda);
double p_igauss(double x, double mu, double lambda);

#endif
