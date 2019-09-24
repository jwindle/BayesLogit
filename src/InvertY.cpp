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



#include "InvertY.h"
#include <stdio.h>

#ifdef USE_R
#include "R.h"
#endif

//------------------------------------------------------------------------------

double y_eval(double v)
{
  double y   = 0.0;
  double r   = sqrt(fabs(v));
  if (v > tol)
    y = tan(r) / r;
  else if (v < -1*tol)
    y = tanh(r) / r;
  else
    y = 1 + (1/3) * v + (2/15) * v * v + (17/315) * v * v * v;
  return y;
}

void ydy_eval(double v, double* yp, double* dyp)
{
  // double r   = sqrt(fabs(v));

  double y = y_eval(v);
  *yp = y;

  if (fabs(v) >= tol) 
    *dyp = 0.5 * (y*y + (1-y) / v);
  else
    *dyp = 0.5 * (y*y - 1/3 - (2/15) * v);

}

double f_eval(double v, void * params)
{
  double y = *((double*) params);
  return y_eval(v) - y;
}

void fdf_eval(double v, void* params, double* fp, double* dfp)
{
  double y = *((double*)params);
  ydy_eval(v, fp, dfp);
  *fp  -= y;
}

double df_eval(double v, void * params)
{
  double f, df;
  ydy_eval(v, &f, &df);
  return df;
}

double v_eval(double y, double tol, int max_iter)
{
  double ylower = ygrid[0];
  double yupper = ygrid[grid_size-1];

  if (y < ylower) {
    return -1. / (y*y);
  } else if (y > yupper) {
    double v = atan(0.5 * y * IYPI);
    return v*v;
  }
  else if (y==1) return 0.0;
    
  double id = (log(y) / log(2.0) + 4.0) / 0.1;
  // printf("y, id, y[id], v[id]: %g, %g, %g, %g\n", y, id, ygrid[(int)id], vgrid[(int)id]);

  // C++ default is truncate decimal portion.
  int idlow  = (int)id;
  int idhigh = (int)id + 1;
  double vl  = vgrid[idlow];  // lower bound
  double vh  = vgrid[idhigh]; // upper bound

  int    iter = 0;
  double diff = tol + 1.0;
  double vnew = vl;
  double vold = vl;
  double f0, f1;

  while (diff > tol && iter < max_iter) {
    iter++;
    vold = vnew;
    fdf_eval(vold, &y, &f0, &f1);
    vnew = vold - f0 / f1;
    vnew = vnew > vh ? vh : vnew;
    vnew = vnew < vl ? vl : vnew;
    diff = fabs(vnew - vold);
    // printf("iter: %i, v: %g, diff: %g\n", iter, vnew, diff);
  }

  if (iter >= max_iter) {
    #ifndef USE_R
    fprintf(stderr, "InvertY.cpp, v_eval: reached max_iter: %i\n", iter);
    #else
    Rprintf("InvertY.cpp, v_eval: reached max_iter: %i\n", iter);
    #endif
  }

  return vnew;
}
