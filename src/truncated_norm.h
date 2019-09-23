// Jesse Windle, 2019
// Method by Christian Robert
// Code modified from <https://github.com/jwindle/BayesLogit>

#include <stdio.h>
#include <stdexcept>
#include <cmath>

#ifndef TRUNCATED_NORM
#define TRUNCATED_NORM

// Throw runtime exception or return.
#ifndef TREOR
#ifndef NTHROW
#define TREOR(MESS, VAL) throw std::runtime_error(MESS);
#else
#define TREOR(MESS, VAL) {fprintf(stderr, MESS); return VAL;}
#endif
#endif

double flat(double left, double right);
double texpon_rate(double left, double rate);
double texpon_rate(double left, double right, double rate);
double alphastar(double left);
double lowerbound(double left);
double tnorm(double left);
double tnorm(double left, double mu, double sd);
double tnorm(double left, double right, double mu, double sd);
double rtinvchi2(double scale, double trunc);

#endif
