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



// MAKE SURE YOU USE GetRNGState() and PutRNGState() !!!

#include "polyagamma_wrapper.h"
#include "PolyaGamma.h"
#include "PolyaGammaApproxAlt.h"
#include "PolyaGammaApproxSP.h"
#include "simple_RNG_wrapper.h"


void rpg_gamma(double *x, double *n, double *z, int *num, int *trunc)
{
    PolyaGamma pg(*trunc);

#ifdef USE_R
    GetRNGstate();
#endif

    for(int i=0; i < *num; ++i){
#ifdef USE_R
        if (i % 1000 == 0) R_CheckUserInterrupt();
#endif
        if (n[i]!=0.0)
            x[i] = pg.draw_sum_of_gammas(n[i], z[i]);
        else
            x[i] = 0.0;

    }

#ifdef USE_R
    PutRNGstate();
#endif
} // rpg



void rpg_devroye(double *x, int *n, double *z, int *num)
{
    PolyaGamma pg(1);

#ifdef USE_R
    GetRNGstate();
#endif

    for(int i=0; i < *num; ++i){
        if (n[i]!=0)
            x[i] = pg.draw(n[i], z[i]);
        else
            x[i] = 0.0;
    }

#ifdef USE_R
    PutRNGstate();
#endif
} 

void rpg_alt(double *x, double *h, double *z, int* num)
{
    PolyaGammaApproxAlt pg;

#ifdef USE_R
    GetRNGstate();
#endif

    for(int i=0; i < *num; ++i){
        if (h[i]!=0)
            x[i] = pg.draw(h[i], z[i]);
        else
            x[i] = 0.0;
    }

#ifdef USE_R
    PutRNGstate();
#endif
}

void rpg_sp(double *x, double *h, double *z, int* num, int *iter)
{
    PolyaGammaApproxSP pg;

#ifdef USE_R
    GetRNGstate();
#endif

    for(int i=0; i < *num; ++i){
        if (h[i]!=0)
            iter[i] = pg.draw(x[i], h[i], z[i]);
        else
            x[i] = 0.0;
    }

#ifdef USE_R
    PutRNGstate();
#endif
}

void rpg_hybrid(double *x, double *h, double *z, int* num)
{
    PolyaGamma dv(1000);
    // PolyaGammaApproxAlt alt;
    PolyaGammaApproxSP sp;

#ifdef USE_R
    GetRNGstate();
#endif

    for(int i=0; i < *num; ++i){
        double b = h[i];
        if (b > 170) {
            double m = dv.pg_m1(b,z[i]);
            double v = dv.pg_m2(b,z[i]) - m*m;
            x[i] = norm(m, sqrt(v));
        }
        else if (b > 13) {
            sp.draw(x[i], b, z[i]);
        }
        else if (b==1 || b==2) {
            x[i] = dv.draw((int)b, z[i]);
        }
	// Need to review "alt" sampler.
        // else if (b > 1) {
        //     x[i] = alt.draw(b, z[i]);
        // }
        else if (b > 0) {
            x[i] = dv.draw_sum_of_gammas(b, z[i]);
        }
        else {
            x[i] = 0.0;
        }
    }

#ifdef USE_R
    PutRNGstate();
#endif
}



// double rpg_hybrid(double b_, double z_)
// {
//     PolyaGamma          dv;
//     PolyaGammaApproxAlt al;
//     PolyaGammaApproxSP  sp;

//     double x;

//     double b = (double) b_;
//     double z = (double) z_;

//     if (b > 170) {
// 	double m = dv.pg_m1(b,z);
// 	double v = dv.pg_m2(b,z) - m*m;
// 	x = (double) norm(m, sqrt(v));
//     }
//     else if (b > 13) {
// 	sp.draw(x, b, z);
//     }
//     else if (b==1 || b==2) {
// 	x = dv.draw((int)b, z);
//     }
//     else if (b > 1) {
// 	x = al.draw(b, z);
//     }
//     else if (b > 0) {
// 	x = dv.draw_sum_of_gammas(b, z);
//     }
//     else {
// 	x = 0.0;
//     }

//     return (double) x;
// }
