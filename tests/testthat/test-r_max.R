# Tests for r_max()

test_that("r_max returns the values reported in Section 3.3", {
  expect_equal(r_max(kappa = 1, alpha = 0.025, beta = 0.3), 0.2863,
               tolerance = 1e-3)
  expect_equal(r_max(kappa = 1, alpha = 0.025, beta = 0.2), 0.2385,
               tolerance = 1e-3)
  expect_equal(r_max(kappa = 1, alpha = 0.025, beta = 0.1), 0.1857,
               tolerance = 1e-3)
})

test_that("r_max at kappa = 1 has the explicit form", {
  for (beta in c(0.3, 0.2, 0.15, 0.1)) {
    z_alpha <- qnorm(1 - 0.025)
    z_beta  <- qnorm(1 - beta)
    expected <- 1 - (z_alpha + z_beta) ^ 2 /
      (z_alpha + qnorm(sqrt(1 - beta))) ^ 2
    expect_equal(r_max(kappa = 1, beta = beta), expected, tolerance = 1e-12)
  }
})

test_that("r_max uses alpha = 0.025 and beta = 0.2 by default", {
  expect_equal(r_max(kappa = 1.3), r_max(kappa = 1.3, alpha = 0.025,
                                         beta = 0.2))
})

test_that("r_max is non-increasing in kappa and bounded by r_max(1)", {
  kappa_seq <- seq(1, 10, by = 0.01)
  for (beta in c(0.3, 0.2, 0.13, 0.1)) {
    vals <- sapply(kappa_seq, r_max, beta = beta)
    # Allow for rounding at the order of machine precision
    expect_true(all(diff(vals) <= 1e-12))
    expect_true(all(vals >= -1e-12))
    expect_true(all(vals <= r_max(kappa = 1, beta = beta)))
  }
})

test_that("r_max decreases as the target power increases", {
  for (kappa in c(1, 1.3, 1.6)) {
    vals <- sapply(c(0.3, 0.2, 0.1), function(b) r_max(kappa, beta = b))
    expect_true(all(diff(vals) < 0))
  }
})

test_that("r_max reproduces the Doody 2014 example (about 4.6%)", {
  kappa <- (3.1 / 12) / (1.8 / 9)
  expect_equal(r_max(kappa = kappa, alpha = 0.025, beta = 0.13), 0.046,
               tolerance = 0.01)
})

test_that("r_max rejects kappa < 1", {
  expect_error(r_max(kappa = 0.9), "kappa must be >= 1")
})
