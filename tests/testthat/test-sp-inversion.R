test_that("saddlepoint inversion converges near y = 1", {
  # The midpoint 1.1 * tanh(z / 2) / (z / 2) is approximately
  # 0.9999999982216754, whose inverse lies in the Taylor-series branch.
  # Integer division in that branch previously exhausted 1,000 iterations.
  z <- 1.106469275917348
  for (sampler in list(rpg.sp, rpg)) {
    set.seed(1)
    expect_output(x <- sampler(5, h = 50, z = z), NA)
    expect_true(all(is.finite(x) & x > 0))
  }
})
