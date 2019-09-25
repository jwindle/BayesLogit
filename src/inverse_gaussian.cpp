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



#include "inverse_gaussian.h"
#include "simple_RNG_wrapper.h"

double igauss(double mu, double lambda)
{
    // See R code for specifics.
    double mu2 = mu * mu;
    double Y = norm(0.0, 1.0);
    Y *= Y;
    double W = mu + 0.5 * mu2 * Y / lambda;
    double X = W - sqrt(W*W - mu2);
    if (unif() > mu / (mu + X))
        X = mu2 / X;
    return X;
}


double p_igauss(double x, double mu, double lambda)
{
    // z = 1 / mean
    double z = 1 / mu;
    double b = sqrt(lambda / x) * (x * z - 1);
    double a = sqrt(lambda / x) * (x * z + 1) * -1.0;
    // double y = p_norm(b, false) + exp(2 * lambda * z) * p_norm(a, false);
    double y = p_norm(b, false) + exp(2 * lambda * z + p_norm(a, true));
    return y;
}
