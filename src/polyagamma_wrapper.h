// -*- mode: c++; c-basic-offset: 4 -*-
// (C) Nicholas Polson, James Scott, Jesse Windle, 2012-2019
// Contributor: Jari Oksanen

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



#ifndef POLYAGAMMA_HYBRID
#define POLYAGAMMA_HYBRID

#ifdef __cplusplus
extern "C" {
#endif

extern void rpg_gamma(double *x, double *n, double *z, int *num, int *trunc);
extern void rpg_devroye(double *x, int *n, double *z, int *num);
extern void rpg_alt(double *x, double *h, double *z, int* num);
extern void rpg_sp(double *x, double *h, double *z, int* num, int *iter);
extern void rpg_hybrid(double *x, double *h, double *z, int* num);

#ifdef __cplusplus
}
#endif

#endif
