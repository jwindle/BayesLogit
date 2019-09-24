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



#ifndef __POLYAGAMMAAPPROXSP__
#define __POLYAGAMMAAPPROXSP__


#include "simple_RNG_wrapper.h"
#include "truncated_norm.h"
#include "truncated_gamma.h"
#include "inverse_gaussian.h"
#include "PolyaGamma.h"

#include <cmath>
#include <vector>
#include <stdio.h>
#include <stdexcept>


// Function and derivative.
struct FD {
  double val;
  double der;
};

struct Line {
  double slope;
  double icept;
};


// PolyaGamma approximation by SP.
class PolyaGammaApproxSP
{

public:

  int draw(double& d, double h, double z, int max_iter=200);

protected:

  // Helper.
  
  double w_left (double trunc, double h, double z);
  double w_right(double trunc, double h, double z);

  void   delta_func(double x, double z  , FD& delta);
  double phi_func  (double x, double mid, FD& phi);
  
  double tangent_to_eta(double x, double z, double mid, Line& tl);

  double sp_approx(double x, double n, double z);

  double cos_rt(double v);

  // YV yv;

  double rtigauss(double mu, double lambda, double trunc);
  double y_func(double v); // y = tan(sqrt(v)) / sqrt(v);

};

//------------------------------------------------------------------------------

// double v_secant(double y, double vb, double va, double tol=1e-8, int maxiter=100);
// double v_func(double y);

#endif
