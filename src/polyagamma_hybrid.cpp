// -*- c-basic-offset: 4 -*-

// Jesse Windle, 2019

#include "polyagamma_hybrid.h"
#include "PolyaGamma.h"
#include "PolyaGammaApproxAlt.h"
#include "PolyaGammaApproxSP.h"
#include "simple_RNG_wrapper.h"


double rpg_hybrid(double b_, double z_)
{
    PolyaGamma          dv;
    PolyaGammaApproxAlt al;
    PolyaGammaApproxSP  sp;

    double x;

    double b = (double) b_;
    double z = (double) z_;

    if (b > 170) {
	double m = dv.pg_m1(b,z);
	double v = dv.pg_m2(b,z) - m*m;
	x = (double) norm(m, sqrt(v));
    }
    else if (b > 13) {
	sp.draw(x, b, z);
    }
    else if (b==1 || b==2) {
	x = dv.draw((int)b, z);
    }
    else if (b > 1) {
	x = al.draw(b, z, r);
    }
    else if (b > 0) {
	x = dv.draw_sum_of_gammas(b, z);
    }
    else {
	x = 0.0;
    }

    return (double) x;
}
