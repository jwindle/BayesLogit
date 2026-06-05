% package: BayesLogit
\name{pg.m1}

\alias{pg.m1}
\alias{pg.m2}
\alias{pg.var}

\title{Polya-Gamma Distribution Moments}

\description{
  Compute the first moment (mean), second moment, and variance of the
  Polya-Gamma distribution PG(b, z).
}

\usage{
pg.m1(b, z)
pg.m2(b, z)
pg.var(b, z)
}

\arguments{
  \item{b}{Shape parameter. Must be positive.}
  \item{z}{Tilt parameter. May be any real number.}
}

\details{

  For \eqn{X \sim \mathrm{PG}(b, z)}, the moments are derived from
  the moment generating function

  \deqn{ E[e^{tX}] = \left( \frac{\cosh(z/2)}{\cosh(\sqrt{z^2/4 - t/2})} \right)^b. }

  The first moment is

  \deqn{ E[X] = \frac{b \tanh(z/2)}{2z}, \quad z \neq 0, }

  with \eqn{E[X] = b/4} when \eqn{z = 0}.  Taylor series are used
  near \eqn{z = 0} for numerical stability.

  \code{pg.m2} returns \eqn{E[X^2]}, and \code{pg.var} returns
  \eqn{\mathrm{Var}(X) = E[X^2] - E[X]^2}.

}

\value{
  A scalar giving the requested moment.
}

\examples{
## Mean and variance of PG(1, 0)
pg.m1(1, 0)   # 0.25
pg.var(1, 0)  # 1/24

## Compare sample moments to theoretical values
set.seed(1)
x <- rpg(10000, h = 2, z = 1)
mean(x)
pg.m1(2, 1)

var(x)
pg.var(2, 1)
}

\seealso{
  \code{\link{rpg}}
}

\references{
  Nicholas G. Polson, James G. Scott, and Jesse Windle.  Bayesian inference for
  logistic models using Polya-Gamma latent variables.
  \url{https://arxiv.org/abs/1205.0310}
}

\keyword{Polya-Gamma}
\keyword{polyagamma}
