// Jesse Windle, 2019
// Code modified from <https://github.com/jwindle/BayesLogit>

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
    double y = p_norm(b, false) + exp(2 * lambda * z) * p_norm(a, false);
    return y;
}
