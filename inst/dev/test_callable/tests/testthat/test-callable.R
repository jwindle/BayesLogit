library(BayesLogit)
library(BayesLogitCallableTest)

N        <- 5000
MEAN_TOL <- 0.05
VAR_TOL  <- 0.10

## Each test calls through R_GetCCallable, verifying both pointer resolution
## and statistical correctness of the samples.

# --- rpg_hybrid (scalar loop) ------------------------------------------------

test_that("callable rpg_hybrid returns positive finite values", {
  set.seed(1)
  x <- callable_rpg_hybrid(50, h = 2, z = 1)
  expect_length(x, 50)
  expect_true(all(x > 0))
  expect_true(all(is.finite(x)))
})

test_that("callable rpg_hybrid mean matches PG(b=2, z=1)", {
  set.seed(1)
  x <- callable_rpg_hybrid(N, h = 2, z = 1)
  expect_equal(mean(x), pg.m1(2, 1), tolerance = MEAN_TOL)
})

test_that("callable rpg_hybrid variance matches PG(b=2, z=1)", {
  set.seed(1)
  x <- callable_rpg_hybrid(N, h = 2, z = 1)
  expect_equal(var(x), pg.var(2, 1), tolerance = VAR_TOL)
})

# --- rpg_hybrid_fill (vectorised) --------------------------------------------

test_that("callable rpg_hybrid_fill returns correct-length positive vector", {
  set.seed(1)
  x <- callable_rpg_hybrid_fill(rep(2, 50), rep(1, 50))
  expect_length(x, 50)
  expect_true(all(x > 0))
  expect_true(all(is.finite(x)))
})

test_that("callable rpg_hybrid_fill mean matches PG(b=2, z=1)", {
  set.seed(1)
  x <- callable_rpg_hybrid_fill(rep(2, N), rep(1, N))
  expect_equal(mean(x), pg.m1(2, 1), tolerance = MEAN_TOL)
})

# --- rpg_gamma ---------------------------------------------------------------

test_that("callable rpg_gamma mean matches PG(b=1, z=0)", {
  set.seed(1)
  x <- callable_rpg_gamma(N, h = 1, z = 0)
  expect_equal(mean(x), pg.m1(1, 0), tolerance = MEAN_TOL)
})

# --- rpg_devroye -------------------------------------------------------------

test_that("callable rpg_devroye mean matches PG(b=2, z=1)", {
  set.seed(1)
  x <- callable_rpg_devroye(N, h = 2L, z = 1)
  expect_equal(mean(x), pg.m1(2, 1), tolerance = MEAN_TOL)
})

# --- rpg_sp ------------------------------------------------------------------

test_that("callable rpg_sp mean matches PG(b=4, z=1)", {
  set.seed(1)
  x <- callable_rpg_sp(N, h = 4, z = 1)
  expect_equal(mean(x), pg.m1(4, 1), tolerance = MEAN_TOL)
})
