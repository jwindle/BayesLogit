library(BayesLogit)

N         <- 10000
MEAN_TOL  <- 0.05   # 5% relative tolerance on sample mean
VAR_TOL   <- 0.10   # 10% relative tolerance on sample variance

################################################################################
## Input validation
################################################################################

test_that("rpg errors on non-positive h", {
  expect_error(rpg(10, h =  0), regexp = "h must be > 0")
  expect_error(rpg(10, h = -1), regexp = "h must be > 0")
})

test_that("rpg.gamma errors on negative h or bad trunc", {
  expect_error(rpg.gamma(10, h = -1),          regexp = "h must be greater than zero")
  expect_error(rpg.gamma(10, h = 1, trunc = 0), regexp = "xtrunc must be > 0")
})

test_that("rpg.devroye errors on negative h", {
  expect_error(rpg.devroye(10, h = -1), regexp = "h must be greater than zero")
})

test_that("rpg.sp errors when h < 1", {
  expect_error(rpg.sp(10, h = 0.5), regexp = "h must be >= 1")
  expect_error(rpg.sp(10, h = 0),   regexp = "h must be >= 1")
})

test_that("rpg.sp.R errors when h < 1", {
  expect_error(rpg.sp.R(10, h = 0.5), regexp = "h must be >= 1")
  expect_error(rpg.sp.R(10, h = 0),   regexp = "h must be >= 1")
})

################################################################################
## Output sanity: length, positivity, finiteness
################################################################################

sanity_cases <- list(
  list(fn = function(n) rpg(n, h = 2, z = 1),           label = "rpg"),
  list(fn = function(n) rpg.gamma(n, h = 2, z = 1),      label = "rpg.gamma"),
  list(fn = function(n) rpg.devroye(n, h = 2L, z = 1),   label = "rpg.devroye"),
  list(fn = function(n) rpg.sp(n, h = 2, z = 1),         label = "rpg.sp"),
  list(fn = function(n) rpg.sp.R(n, h = 2, z = 1),       label = "rpg.sp.R"),
  list(fn = function(n) rpg.devroye.R(n, h = 2L, z = 1), label = "rpg.devroye.R"),
  list(fn = function(n) rpg.gamma.R(n, h = 2, z = 1),    label = "rpg.gamma.R")
)

for (sc in sanity_cases) {
  local({
    fn    <- sc$fn
    label <- sc$label
    test_that(paste(label, "returns correct-length, positive, finite vector"), {
      set.seed(42)
      x <- fn(50)
      expect_length(x, 50)
      expect_true(all(x > 0))
      expect_true(all(is.finite(x)))
    })
  })
}

################################################################################
## Moment accuracy: sample mean and variance vs. pg.m1 / pg.var
################################################################################

moment_cases <- list(
  list(b = 1, z = 0),
  list(b = 1, z = 2),
  list(b = 2, z = 0),
  list(b = 4, z = 1)
)

moment_samplers <- list(
  list(fn = function(b, z) rpg(N, h = b, z = z),            label = "rpg"),
  list(fn = function(b, z) rpg.gamma(N, h = b, z = z),       label = "rpg.gamma"),
  list(fn = function(b, z) rpg.devroye(N, h = b, z = z),     label = "rpg.devroye"),
  list(fn = function(b, z) rpg.sp(N, h = b, z = z),          label = "rpg.sp"),
  list(fn = function(b, z) rpg.sp.R(N, h = b, z = z),        label = "rpg.sp.R"),
  list(fn = function(b, z) rpg.devroye.R(N, h = b, z = z),   label = "rpg.devroye.R"),
  list(fn = function(b, z) rpg.gamma.R(N, h = b, z = z),     label = "rpg.gamma.R")
)

for (ms in moment_samplers) {
  for (tc in moment_cases) {
    local({
      fn    <- ms$fn
      label <- ms$label
      b     <- tc$b
      z     <- tc$z
      test_that(paste(label, "mean matches PG(b =", b, ", z =", z, ")"), {
        set.seed(1)
        x <- fn(b, z)
        expect_equal(mean(x), pg.m1(b, z), tolerance = MEAN_TOL)
      })
      test_that(paste(label, "variance matches PG(b =", b, ", z =", z, ")"), {
        set.seed(1)
        x <- fn(b, z)
        expect_equal(var(x), pg.var(b, z), tolerance = VAR_TOL)
      })
    })
  }
}

################################################################################
## Saddlepoint sampler: large-h accuracy
## The SP approximation tightens as h grows (CLT), so these are the most
## meaningful accuracy checks for rpg.sp specifically.
################################################################################

sp_large_h_cases <- list(
  list(b = 20, z = 0),
  list(b = 50, z = 2),
  list(b = 100, z = 4)
)

for (tc in sp_large_h_cases) {
  local({
    b <- tc$b
    z <- tc$z
    for (lbl in c("rpg.sp", "rpg.sp.R")) {
      local({
        sampler <- if (lbl == "rpg.sp") rpg.sp else rpg.sp.R
        nm      <- lbl
        test_that(paste(nm, "mean matches PG(b =", b, ", z =", z, ") [large h]"), {
          set.seed(1)
          x <- sampler(N, h = b, z = z)
          expect_equal(mean(x), pg.m1(b, z), tolerance = MEAN_TOL)
        })
        test_that(paste(nm, "variance matches PG(b =", b, ", z =", z, ") [large h]"), {
          set.seed(1)
          x <- sampler(N, h = b, z = z)
          expect_equal(var(x), pg.var(b, z), tolerance = VAR_TOL)
        })
      })
    }
  })
}
