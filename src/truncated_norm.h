// -*- mode: c++; c-basic-offset: 4 -*-
// (C) Nicholas Polson, James Scott, Jesse Windle, 2012-2019

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



// Method from Christian Robert, *Simulation of truncated normal variables*

// Code modified from <https://github.com/jwindle/BayesLogit>

#include <stdio.h>
#include <stdexcept>
#include <cmath>
#include "R.h"

#ifndef TRUNCATED_NORM
#define TRUNCATED_NORM

// Throw runtime exception or print message and return.
#ifndef TREOR
#ifndef NTHROW
#define TREOR(MESS, VAL) throw std::runtime_error(MESS);
#else
#ifndef USE_R
#define TREOR(MESS, VAL) {fprintf(stderr, MESS); return VAL;}
#else
#define TREOR(MESS, VAL) {Rprintf(MESS); return VAL;}
#endif
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
