# Tests for twocpcont_ss()

# Scenario grid of the numerical study (Section 3 of the manuscript)
scenarios <- list(
  list(d1 = 0.3, d2 = 0.3, s1 = 1, s2 = 1.0, r = 1),
  list(d1 = 0.3, d2 = 0.3, s1 = 1, s2 = 1.0, r = 2),
  list(d1 = 0.3, d2 = 0.3, s1 = 1, s2 = 1.5, r = 1),
  list(d1 = 0.3, d2 = 0.3, s1 = 1, s2 = 1.5, r = 2),
  list(d1 = 0.7, d2 = 0.5, s1 = 1, s2 = 1.0, r = 1),
  list(d1 = 0.7, d2 = 0.5, s1 = 1, s2 = 1.0, r = 2),
  list(d1 = 0.7, d2 = 0.5, s1 = 1, s2 = 1.5, r = 1),
  list(d1 = 0.7, d2 = 0.5, s1 = 1, s2 = 1.5, r = 2)
)
rho_vals <- c(0.0, 0.3, 0.5, 0.8)

# Helper: sample size for one scenario
ss_for <- function(sc, rho, beta, method, gl_nodes = 5L) {
  twocpcont_ss(delta1 = sc$d1, delta2 = sc$d2, sd1 = sc$s1, sd2 = sc$s2,
               rho = rho, r = sc$r, alpha = 0.025, beta = beta,
               method = method, gl_nodes = gl_nodes)
}

test_that("twocpcont_ss returns an object of the correct class", {
  ss <- twocpcont_ss(delta1 = 1, delta2 = 1, sd1 = 2, sd2 = 2,
                     rho = 0.5, r = 1, alpha = 0.025, beta = 0.2,
                     method = "univariate")
  expect_s3_class(ss, "twocpcont_ss")
  expect_s3_class(ss, "data.frame")
  expect_named(ss, c("delta1", "delta2", "sd1", "sd2", "rho", "r", "alpha",
                     "beta", "method", "n1", "n2", "N"))
  expect_equal(ss$N, ss$n1 + ss$n2)
})

test_that("twocpcont_ss uses the bivariate method by default", {
  ss <- twocpcont_ss(delta1 = 1, delta2 = 1, sd1 = 2, sd2 = 2, rho = 0.5)
  expect_equal(ss$method, "bivariate")
  expect_equal(c(ss$r, ss$alpha, ss$beta), c(1, 0.025, 0.2))
})

test_that("twocpcont_ss rejects an unknown method", {
  expect_error(twocpcont_ss(delta1 = 1, delta2 = 1, sd1 = 2, sd2 = 2,
                            rho = 0.5, method = "exact"))
})

test_that("twocpcont_ss: methods agree on all 32 cells at 80% and 85% power", {
  for (beta in c(0.20, 0.15)) {
    for (sc in scenarios) {
      for (rho in rho_vals) {
        expect_equal(ss_for(sc, rho, beta, "univariate")$N,
                     ss_for(sc, rho, beta, "bivariate")$N)
      }
    }
  }
})

test_that("twocpcont_ss: methods agree at 90% power except one boundary cell", {
  # At (0.3, 0.3, 1, 1, r = 1, rho = 0.5) the exact continuous n2 is
  # 278.0009 and the closed-form value is 277.9970, so the ceiling
  # gives n2 = 278 instead of 279.  The power at n2 = 278 falls short
  # of 0.90 by about 1e-6.
  for (sc in scenarios) {
    for (rho in rho_vals) {
      N_uni <- ss_for(sc, rho, 0.10, "univariate")$N
      N_biv <- ss_for(sc, rho, 0.10, "bivariate")$N
      is_boundary <- sc$d1 == 0.3 && sc$d2 == 0.3 && sc$s2 == 1 &&
        sc$r == 1 && rho == 0.5
      if (is_boundary) {
        expect_equal(c(N_uni, N_biv), c(556, 558))
      } else {
        expect_equal(N_uni, N_biv)
      }
    }
  }
  pw <- twocpcont_power(n1 = 278, n2 = 278,
                        delta1 = 0.3, delta2 = 0.3, sd1 = 1, sd2 = 1,
                        rho = 0.5, method = "bivariate")
  expect_lt(pw$powerCoprimary, 0.90)
  expect_gt(pw$powerCoprimary, 0.90 - 1e-5)
})

test_that("twocpcont_ss: 7 quadrature nodes give the same N as 5 nodes", {
  for (sc in scenarios) {
    for (rho in rho_vals) {
      expect_equal(ss_for(sc, rho, 0.2, "univariate", gl_nodes = 7)$N,
                   ss_for(sc, rho, 0.2, "univariate", gl_nodes = 5)$N)
    }
  }
})

test_that("twocpcont_ss: bivariate search returns the smallest adequate n2", {
  for (sc in scenarios[c(1, 4, 7)]) {
    ss <- ss_for(sc, 0.5, 0.2, "bivariate")
    pw_at <- function(n2) {
      twocpcont_power(n1 = ceiling(sc$r * n2), n2 = n2,
                      delta1 = sc$d1, delta2 = sc$d2,
                      sd1 = sc$s1, sd2 = sc$s2, rho = 0.5,
                      method = "bivariate")$powerCoprimary
    }
    expect_gte(pw_at(ss$n2), 0.8)
    expect_lt(pw_at(ss$n2 - 1), 0.8)
  }
})

test_that("twocpcont_ss: n1 equals ceiling(r n2)", {
  for (m in c("univariate", "bivariate")) {
    ss <- twocpcont_ss(delta1 = 0.7, delta2 = 0.5, sd1 = 1, sd2 = 1.5,
                       rho = 0.5, r = 2, method = m)
    expect_equal(ss$n1, ceiling(2 * ss$n2))
    ss <- twocpcont_ss(delta1 = 0.7, delta2 = 0.5, sd1 = 1, sd2 = 1.5,
                       rho = 0.5, r = 1.5, method = m)
    expect_equal(ss$n1, ceiling(1.5 * ss$n2))
  }
})

test_that("twocpcont_ss: N is non-increasing in rho", {
  for (sc in scenarios) {
    N_vec <- sapply(rho_vals, function(rr) {
      ss_for(sc, rr, 0.2, "univariate")$N
    })
    expect_true(all(diff(N_vec) <= 0))
  }
})

test_that("twocpcont_ss: N does not depend on the order of the endpoints", {
  for (m in c("univariate", "bivariate")) {
    ss_a <- twocpcont_ss(delta1 = 0.7, delta2 = 0.5, sd1 = 1, sd2 = 1.5,
                         rho = 0.5, r = 2, method = m)
    ss_b <- twocpcont_ss(delta1 = 0.5, delta2 = 0.7, sd1 = 1.5, sd2 = 1,
                         rho = 0.5, r = 2, method = m)
    expect_equal(ss_a$N, ss_b$N)
  }
})

test_that("twocpcont_ss: N depends on delta and sd only through delta / sd", {
  ss_a <- twocpcont_ss(delta1 = 0.7, delta2 = 0.5, sd1 = 1, sd2 = 1.5,
                       rho = 0.5, r = 2, method = "univariate")
  ss_b <- twocpcont_ss(delta1 = 7, delta2 = 5, sd1 = 10, sd2 = 15,
                       rho = 0.5, r = 2, method = "univariate")
  expect_equal(ss_a$N, ss_b$N)
})

test_that("twocpcont_ss at kappa = 1 and rho = 0 has the explicit form", {
  # lambda_w* = z_alpha + qnorm(sqrt(1 - beta)) and c = 2 / delta^2
  lam <- qnorm(1 - 0.025) + qnorm(sqrt(0.8))
  n2  <- ceiling(2 / 0.3 ^ 2 * lam ^ 2)
  ss  <- twocpcont_ss(delta1 = 0.3, delta2 = 0.3, sd1 = 1, sd2 = 1,
                      rho = 0, r = 1, method = "univariate")
  expect_equal(ss$n2, n2)
  expect_equal(ss$N, 460)
})

test_that("Doody 2014 example reproduces N = 1000 at rho = 0", {
  ss <- twocpcont_ss(delta1 = 1.8, delta2 = 3.1, sd1 = 9, sd2 = 12,
                     rho = 0, r = 1, alpha = 0.025, beta = 0.13,
                     method = "univariate")
  expect_equal(ss$N, 1000)
})

test_that("Doody 2014 example: N over rho matches both methods", {
  N_uni <- sapply(c(0, 0.3, 0.5, 0.8), function(rr) {
    twocpcont_ss(delta1 = 1.8, delta2 = 3.1, sd1 = 9, sd2 = 12,
                 rho = rr, r = 1, alpha = 0.025, beta = 0.13,
                 method = "univariate")$N
  })
  N_biv <- sapply(c(0, 0.3, 0.5, 0.8), function(rr) {
    twocpcont_ss(delta1 = 1.8, delta2 = 3.1, sd1 = 9, sd2 = 12,
                 rho = rr, r = 1, alpha = 0.025, beta = 0.13,
                 method = "bivariate")$N
  })
  expect_equal(N_uni, c(1000, 990, 980, 960))
  expect_equal(N_biv, N_uni)
})

test_that("twocpcont_ss runs at the boundary values of rho", {
  for (rho in c(0, 0.99)) {
    ss <- twocpcont_ss(delta1 = 0.3, delta2 = 0.3, sd1 = 1, sd2 = 1,
                       rho = rho, r = 1, method = "univariate")
    expect_true(is.finite(ss$N))
    expect_gt(ss$N, 0)
  }
})
