# Tests for kappa_star()

test_that("kappa_star returns the values reported in the manuscript", {
  expect_equal(kappa_star(epsilon = 0.010, beta = 0.2), 1.614,
               tolerance = 1e-3)
  expect_equal(kappa_star(epsilon = 0.050, beta = 0.2), 1.369,
               tolerance = 1e-3)
  expect_equal(kappa_star(epsilon = 0.005, beta = 0.2), 1.700,
               tolerance = 1e-3)
  expect_equal(kappa_star(epsilon = 0.005, beta = 0.1), 1.512,
               tolerance = 1e-3)
  expect_equal(kappa_star(epsilon = 0.050, beta = 0.1), 1.239,
               tolerance = 1e-3)
})

test_that("kappa_star uses alpha = 0.025 and beta = 0.2 by default", {
  expect_equal(kappa_star(epsilon = 0.01),
               kappa_star(epsilon = 0.01, alpha = 0.025, beta = 0.2))
})

test_that("kappa_star solves the exact rho = 0 relation", {
  # At kappa = kappa_star, lambda = (z_alpha + z_beta) / sqrt(1 - eps)
  # satisfies Phi(lambda - z_alpha) Phi(kappa lambda - z_alpha) = 1 - beta.
  for (beta in c(0.3, 0.2, 0.1)) {
    for (eps in c(0.005, 0.01, 0.05, 0.1)) {
      z_alpha <- qnorm(1 - 0.025)
      lam <- (z_alpha + qnorm(1 - beta)) / sqrt(1 - eps)
      ks  <- kappa_star(epsilon = eps, beta = beta)
      expect_equal(pnorm(lam - z_alpha) * pnorm(ks * lam - z_alpha),
                   1 - beta, tolerance = 1e-10)
    }
  }
})

test_that("kappa_star is the exact inverse of r_max", {
  for (beta in c(0.3, 0.2, 0.15, 0.1)) {
    for (eps in c(0.005, 0.01, 0.02, 0.05, 0.1, 0.15)) {
      ks <- kappa_star(epsilon = eps, beta = beta)
      expect_equal(r_max(kappa = ks, beta = beta), eps, tolerance = 1e-8)
      expect_equal(r_max(kappa = ks, beta = beta, exact = FALSE), eps,
                   tolerance = 1e-3)
    }
  }
})

test_that("kappa_star exceeds 1 and decreases in epsilon and in power", {
  eps_seq <- c(0.005, 0.01, 0.02, 0.05, 0.1, 0.15)
  ks_80 <- sapply(eps_seq, kappa_star, beta = 0.2)
  ks_90 <- sapply(eps_seq, kappa_star, beta = 0.1)
  expect_true(all(ks_80 > 1))
  expect_true(all(ks_90 > 1))
  expect_true(all(diff(ks_80) < 0))
  expect_true(all(ks_90 < ks_80))
})

test_that("kappa_star equals 1 when epsilon is at least r_max(1)", {
  for (beta in c(0.2, 0.1)) {
    r1 <- r_max(kappa = 1, beta = beta)
    expect_equal(kappa_star(epsilon = r1, beta = beta), 1, tolerance = 1e-8)
    expect_identical(kappa_star(epsilon = 0.5, beta = beta), 1)
    expect_gt(kappa_star(epsilon = r1 - 0.01, beta = beta), 1)
  }
})

test_that("kappa_star rejects epsilon outside (0, 1)", {
  expect_error(kappa_star(epsilon = 0), "epsilon must lie in")
  expect_error(kappa_star(epsilon = 1), "epsilon must lie in")
  expect_error(kappa_star(epsilon = -0.1), "epsilon must lie in")
})
