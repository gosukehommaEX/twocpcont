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

test_that("r_max uses alpha = 0.025, beta = 0.2, exact = TRUE by default", {
  expect_equal(r_max(kappa = 1.3), r_max(kappa = 1.3, alpha = 0.025,
                                         beta = 0.2, exact = TRUE))
})

test_that("r_max is non-increasing in kappa and bounded by r_max(1)", {
  kappa_seq <- seq(1, 10, by = 0.01)
  for (beta in c(0.3, 0.2, 0.13, 0.1)) {
    for (ex in c(TRUE, FALSE)) {
      vals <- sapply(kappa_seq, r_max, beta = beta, exact = ex)
      # Allow for rounding at the order of machine precision
      expect_true(all(diff(vals) <= 1e-12))
      expect_true(all(vals >= -1e-12))
      expect_true(all(vals <= r_max(kappa = 1, beta = beta) + 1e-12))
    }
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

test_that("r_max with exact = TRUE solves the rho = 0 equation", {
  z_alpha <- qnorm(1 - 0.025)
  for (beta in c(0.3, 0.2, 0.1)) {
    for (kappa in c(1.05, 1.3, 1.6, 2.5)) {
      z_beta <- qnorm(1 - beta)
      lam <- (z_alpha + z_beta) / sqrt(1 - r_max(kappa, beta = beta))
      expect_equal(pnorm(lam - z_alpha) * pnorm(kappa * lam - z_alpha),
                   1 - beta, tolerance = 1e-10)
    }
  }
})

test_that("r_max: exact and closed-form versions agree closely", {
  kappa_seq <- seq(1, 3, by = 0.05)
  for (beta in c(0.3, 0.2, 0.15, 0.13, 0.1)) {
    v_exact  <- sapply(kappa_seq, r_max, beta = beta, exact = TRUE)
    v_approx <- sapply(kappa_seq, r_max, beta = beta, exact = FALSE)
    expect_lt(max(abs(v_exact - v_approx)), 5e-5)
  }
  expect_equal(r_max(1, beta = 0.2, exact = FALSE),
               r_max(1, beta = 0.2, exact = TRUE), tolerance = 1e-12)
})

test_that("r_max with exact = TRUE returns 0 for very large kappa", {
  expect_identical(r_max(kappa = 10, beta = 0.2), 0)
})

test_that("r_max rejects kappa < 1", {
  expect_error(r_max(kappa = 0.9), "kappa must be >= 1")
})

test_that("r_max rejects an invalid exact argument", {
  expect_error(r_max(kappa = 1.2, exact = NA), "exact must be")
  expect_error(r_max(kappa = 1.2, exact = "yes"), "exact must be")
  expect_error(r_max(kappa = 1.2, exact = c(TRUE, FALSE)), "exact must be")
})

test_that("r_max(1, beta) reproduces the sample size ratio of Hung and Wang (2009)", {
  # Hung and Wang (2009, Section 3) give, for K independent endpoints with
  # equal effect sizes, n / m1 = [{qnorm(1 - alpha) + qnorm((1 - beta)^(1/K))} /
  # {qnorm(1 - alpha) + qnorm(1 - beta)}]^2, where m1 is the single-endpoint
  # sample size.  For K = 2 this equals 1 / {1 - r_max(1, beta)}.
  z_a <- qnorm(1 - 0.025)
  for (beta in c(0.20, 0.15, 0.10)) {
    ratio_hw <- ((z_a + qnorm(sqrt(1 - beta))) / (z_a + qnorm(1 - beta))) ^ 2
    expect_equal(1 / (1 - r_max(1, alpha = 0.025, beta = beta)), ratio_hw,
                 tolerance = 1e-10)
  }
  # Values at K = 2 read from Figure 1 of Hung and Wang (2009):
  # about 1.31, 1.27, and 1.23 for power 0.80, 0.85, and 0.90
  ratio <- sapply(c(0.20, 0.15, 0.10), function(b) 1 / (1 - r_max(1, beta = b)))
  # Absolute difference: expect_equal() uses a relative tolerance
  expect_lt(max(abs(ratio - c(1.31, 1.27, 1.23))), 0.01)
})
