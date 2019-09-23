// -*- c-basic-offset: 4; -*-

// Jesse Windle, 2019
// Method from Christian Robert
// Code modified from <https://github.com/jwindle/BayesLogit>

#include "truncated_norm.h"
#include "simple_RNG_wrapper.h"


double flat(double a, double b) {
    return a + unif() * (b-a);
}

double texpon_rate(double left, double rate){
    if (rate < 0) TREOR("texpon_rate: rate < 0, return 0\n", 0.0);
    // return left - log(unif()) / rate;
    return expon_rate(rate) + left;
}

double texpon_rate(double left, double right, double rate)
{
    if (left == right) return left;
    if (left > right) TREOR("texpon_rate: left > right, return 0.\n", 0.0);
    if (rate < 0) TREOR("texpon_rate: rate < 0, return 0\n", 0.0);

    double b = 1 - exp(rate * (left - right));
    double y = 1 - b * unif();
    return left - log(y) / rate;
}

double alphastar(double left)
{
    return 0.5 * (left + sqrt(left*left + 4));
}

double lowerbound(double left)
{
    double astar  = alphastar(left);
    double lbound = left + exp(0.5 * + 0.5 * left * (left - astar)) / astar;
    return lbound;
} 

double tnorm(double left)
{
    double rho, ppsl;
    int count = 1;

    if (left < 0) { // Accept/Reject Normal
        while (true) {
            ppsl = norm(0.0, 1.0);
            if (ppsl > left) return ppsl;
            check_R_interupt(count++);
            #ifndef NDEBUG
            if (count > RCHECK * 1000) fprintf(stderr, "left < 0; count: %i\n", count);
            #endif
        }
    }
    else { // Accept/Reject Exponential
        // return tnorm_tail(left); // Use Devroye.
        double astar = alphastar(left);
        while (true) {
            ppsl = texpon_rate(left, astar);
            rho  = exp( -0.5 * (ppsl - astar) * (ppsl - astar) );
            if (unif() < rho) return ppsl;
            check_R_interupt(count++);
            #ifndef NDEBUG
            if (count > RCHECK * 1000) fprintf(stderr, "left > 0; count: %i\n", count);
            #endif
        }
    }

}

double tnorm(double left, double right)
{
    // The most difficult part of this algorithm is figuring out all the
    // various cases.  An outline is summarized in the Appendix.

    // Check input
    #ifdef USE_R
    if (ISNAN(right) || ISNAN(left))
    #else
    if (std::isnan(right) || std::isnan(left))
    #endif
	{
	    fprintf(stderr, "Warning: nan sent to tnorm: left=%g, right=%g\n", left, right);
	    TREOR("tnorm: parameter problem.\n", 0.5 * (left + right));
	    // throw std::runtime_error("tnorm: parameter problem.\n");
	}
    
    if (right < left) {
        fprintf(stderr, "Warning: left: %g, right:%g.\n", left, right);
        TREOR("tnorm: parameter problem.\n", 0.5 * (left + right));
    }
    
    double rho, ppsl;
    int count = 1;
    
    if (left >= 0) {
        double lbound = lowerbound(left);
        if (right > lbound) { // Truncated Exponential.
            double astar = alphastar(left);
            while (true) {
		ppsl = texpon_rate(left, right, astar);
                rho  = exp(-0.5*(ppsl - astar)*(ppsl-astar));
                if (unif() < rho) return ppsl;
		if (count > RCHECK * 10) fprintf(stderr, "left >= 0, right > lbound; count: %i\n", count);
                // if (ppsl < right) return ppsl;
            }
        }
        else {
            while (true) {
                ppsl = flat(left, right);
                rho  = exp(0.5 * (left*left - ppsl*ppsl));
                if (unif() < rho) return ppsl;
                check_R_interupt(count++);
                #ifndef NDEBUG
                if (count > RCHECK * 10) fprintf(stderr, "left >= 0, right <= lbound; count: %i\n", count);
                #endif
            }
        }
    }
    else if (right >= 0) {
        if ( (right - left) < SQRT2PI ){
            while (true) {
                ppsl = flat(left, right);
                rho  = exp(-0.5 * ppsl * ppsl);
                if (unif() < rho) return ppsl;
                check_R_interupt(count++);
                #ifndef NDEBUG
                if (count > RCHECK * 10) fprintf(stderr, "First, left < 0, right >= 0, count: %i\n", count);
                #endif
            }
        }
        else{
            while (true) {
                ppsl = norm(0., 1.);
                if (left < ppsl && ppsl < right) return ppsl;
                check_R_interupt(count++);
                #ifndef NDEBUG
                if (count > RCHECK * 10) fprintf(stderr, "Second, left < 0, right > 0, count: %i\n", count);
                #endif
            }
        }
    }
    else {
        return -1. * tnorm(-1.0 * right, -1.0 * left);
    }
} 

double tnorm(double left, double mu, double sd)
{
    double newleft = (left - mu) / sd;
    return mu + tnorm(newleft) * sd;
}

double tnorm(double left, double right, double mu, double sd)
{
    if (left==right) return left;

    double newleft  = (left - mu) / sd;
    double newright = (right - mu) / sd;

    // I want to check this here as well so we can see what the input was.
    // It may be more elegant to try and catch tdraw.
    if (newright < newleft) {
        fprintf(stderr, "left, right, mu, sd: %g, %g, %g, %g \n", left, right, mu, sd);
        fprintf(stderr, "nleft, nright: %g, %g\n", newleft, newright);
        TREOR("tnorm: parameter problem.\n", 0.5 * (left + right));
    }

    double tdraw = tnorm(newleft, newright);
    double draw = mu + tdraw * sd;

    // It may be the case that there is some numerical error and that the draw
    // ends up out of bounds.
    if (draw < left || draw > right){
        fprintf(stderr, "Error in tnorm: draw not in bounds.\n");
        fprintf(stderr, "left, right, mu, sd: %g, %g, %g, %g\n", left, right, mu, sd);
        fprintf(stderr, "nleft, nright, tdraw, draw: %g, %g, %g, %g\n", newleft, newright, tdraw, draw);
        TREOR("Aborting and returning average of left and right.\n",  0.5 * (left + right));
    }

    return draw;
}

double rtinvchi2(double scale, double trunc)
{
    double R = trunc / scale;
    // double X = 0.0;
    // // I need to consider using a different truncated normal sampler.
    // double E1 = r.expon_rate(1.0); double E2 = r.expon_rate(1.0);
    // while ( (E1*E1) > (2 * E2 / R)) {
    //   // printf("E %g %g %g %g\n", E1, E2, E1*E1, 2*E2/R);
    //   E1 = r.expon_rate(1.0); E2 = r.expon_rate(1.0);
    // }
    // // printf("E %g %g \n", E1, E2);
    // X = 1 + E1 * R;
    // X = R / (X * X);
    // X = scale * X;
    double E = tnorm(1/sqrt(R));
    double X = scale / (E*E);
    return X;
}
