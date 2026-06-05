library(BayesLogit)

N        <- 10000
MEAN_TOL <- 0.05  # 5% relative tolerance
VAR_TOL  <- 0.10  # 10% relative tolerance
TRUNC    <- 0.64  # truncation point used by rtigauss

################################################################################
## rigauss: Inverse-Gaussian(mu, lambda)
## E[X] = mu,  Var[X] = mu^3 / lambda
################################################################################

rigauss_cases <- list(
  list(mu = 1.0, lambda = 1.0),
  list(mu = 0.5, lambda = 2.0),
  list(mu = 2.0, lambda = 1.0),
  list(mu = 1.0, lambda = 4.0)
)

for (tc in rigauss_cases) {
  local({
    mu     <- tc$mu
    lambda <- tc$lambda
    test_that(paste("rigauss mean matches IG(mu =", mu, ", lambda =", lambda, ")"), {
      set.seed(1)
      x <- replicate(N, BayesLogit:::rigauss(mu, lambda))
      expect_true(all(x > 0))
      expect_equal(mean(x), mu, tolerance = MEAN_TOL)
    })
    test_that(paste("rigauss variance matches IG(mu =", mu, ", lambda =", lambda, ")"), {
      set.seed(1)
      x <- replicate(N, BayesLogit:::rigauss(mu, lambda))
      expect_equal(var(x), mu^3 / lambda, tolerance = VAR_TOL)
    })
  })
}

################################################################################
## rtigauss: truncated Inverse-Gaussian(1/|Z|, 1) on (0, TRUNC)
##
## No closed-form moments for the truncated distribution, so we validate by:
##   1. Checking all samples fall in (0, TRUNC)
##   2. Comparing mean and variance against a rejection sampler from the
##      same distribution (draw from IG(1/Z, 1), discard if > TRUNC)
################################################################################

rejection_rtigauss <- function(Z, n) {
  replicate(n, {
    x <- TRUNC + 1
    while (x > TRUNC) x <- BayesLogit:::rigauss(1 / abs(Z), 1.0)
    x
  })
}

rtigauss_Z_values <- c(0.5, 1.0, 2.0, 4.0)

for (Z in rtigauss_Z_values) {
  local({
    z <- Z
    test_that(paste("rtigauss samples lie in (0, TRUNC) for Z =", z), {
      set.seed(1)
      x <- replicate(N, BayesLogit:::rtigauss(z))
      expect_true(all(x > 0))
      expect_true(all(x < TRUNC))
    })
    test_that(paste("rtigauss mean matches rejection sampler for Z =", z), {
      set.seed(1); x <- replicate(N, BayesLogit:::rtigauss(z))
      set.seed(1); y <- rejection_rtigauss(z, N)
      expect_equal(mean(x), mean(y), tolerance = MEAN_TOL)
    })
    test_that(paste("rtigauss variance matches rejection sampler for Z =", z), {
      set.seed(1); x <- replicate(N, BayesLogit:::rtigauss(z))
      set.seed(1); y <- rejection_rtigauss(z, N)
      expect_equal(var(x), var(y), tolerance = VAR_TOL)
    })
  })
}
