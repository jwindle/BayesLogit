% package: BayesLogit
\name{rpg}

\alias{rpg}
\alias{rpg.devroye}
\alias{rpg.gamma}
%\alias{rpg.alt}
\alias{rpg.sp}
\alias{rpg.gamma.R}
\alias{rpg.devroye.R}
\alias{rpg.sp.R}

\title{Polya-Gamma Random Variates}

\description{
  Generate random variates from the Polya-Gamma distribution.
}

\details{

  A random variable X with distribution PG(h,z) is distributed like

  \deqn{ X \sim \sum_{k=1}^\infty G(h,1) / ( 2 \pi^2 (k-1/2)^2 + z^2/2). }

  The density for X may be derived by exponentially tilting the
  PG(h,0) density:

  \deqn{ p(x|h,z) \propto \exp(-x z^2/2) p(x|h,0).}

  Different methods for generating this random variable are
  implemented, each of which is useful for certain parameters.  The
  parameters supplied by the user automatically determine which method
  is used.  One may manually call each routine using
  \code{rpg.METHOD}.  Functions ending in ".R" are pure R
  implementations.

  You may call rpg when n and z are vectors.

}

\usage{

rpg(num=1, h=1, z=0.0)

rpg.gamma(num=1, h=1, z=0.0, trunc=200)

rpg.devroye(num=1, h=1, z=0.0)

rpg.sp(num=1, h=1, z=0.0)

rpg.gamma.R(num=1, h=1, z=0.0, trunc=200)

rpg.devroye.R(num=1, h=1, z=0.0)

}

% rpg.sp.R(num=1, h=1, z=0.0)

% rpg.alt(num=1, h=1, z=0.0)


\arguments{

  \item{num}{The number of random variates to simulate.}
  \item{h}{Shape parameter.  \code{h} must be >= 1 if not using sum of gammas method.}
  \item{z}{Parameter associated with tilting.}
  \item{trunc}{The number of elements used in sum of gammas approximation.}
}

%  \item{track.iter}{The number of proposals made before accepting.}

\value{

  This function returns \code{num} Polya-Gamma samples.

}

\examples{

h = c(1, 2, 3);
z = c(4, 5, 6);

## Devroye-like method -- only use if h contains integers, preferably small integers.
X = rpg.devroye(100, h, z);

h = c(1.2, 2.3, 3.2);
z = c(4, 5, 6);

## Sum of gammas method -- this is slow.
X = rpg.gamma(100, h, z);

h = c(1, 4, 2.3);
z = c(4, 5, 6);

## Hybrid method -- automatically chooses best procedure.
X = rpg(100, h, z);

}

%\seealso{
%  \code{\link{logit.EM}}, \code{\link{logit}}
%}

\references{
  Nicholas G. Polson, James G. Scott, and Jesse Windle.  Bayesian inference for
  logistic models using Polya-Gamma latent variables.
  \url{http://arxiv.org/abs/1205.0310}
}

\keyword{rpg}
\keyword{Polya-Gamma}
\keyword{polyagamma}
