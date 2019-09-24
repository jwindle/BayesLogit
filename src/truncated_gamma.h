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



// See Philippe, *Simulation of right and left truncated gamma
// distributions by mixtures* and Dagpunar *Sampling of variates from
// a truncated gamma distribution*

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
