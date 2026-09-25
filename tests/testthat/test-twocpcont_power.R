# Tests for twocpcont_power()

test_that("twocpcont_power returns an object of the correct class", {
  pw <- twocpcont_power(n1 = 100, n2 = 100,
                        delta1 = 1, delta2 = 1, sd1 = 2, sd2 = 2,
                        rho = 0.5, alpha = 0.025,
                        method = "univariate")
  expect_s3_class(pw, "twocpcont_power")
  expect_s3_class(pw, "data.frame")
  expect_named(pw, c("n1", "n2", "delta1", "delta2", "sd1", "sd2", "rho",
                     "alpha", "method", "power1", "power2",
                     "powerCoprimary"))
  expect_equal(nrow(pw), 1L)
})

test_that("twocpcont_power uses the bivariate method by default", {
  pw <- twocpcont_power(n1 = 100, n2 = 100,
                        delta1 = 1, delta2 = 1, sd1 = 2, sd2 = 2,
                        rho = 0.5)
  expect_equal(pw$method, "bivariate")
  expect_equal(pw$alpha, 0.025)
})

test_that("twocpcont_power rejects an unknown method", {
  expect_error(twocpcont_power(n1 = 100, n2 = 100,
                               delta1 = 1, delta2 = 1, sd1 = 2, sd2 = 2,
                               rho = 0.5, method = "exact"))
})

test_that("twocpcont_power: univariate and bivariate methods agree", {
  cases <- list(
    list(n1 = 230, n2 = 230, d1 = 0.3, d2 = 0.3, s1 = 1, s2 = 1.0,
         rho = 0.5, tol = 1e-6),
    list(n1 = 100, n2 = 50,  d1 = 0.7, d2 = 0.5, s1 = 1, s2 = 1.0,
         rho = 0.3, tol = 1e-6),
    list(n1 = 150, n2 = 150, d1 = 0.7, d2 = 0.5, s1 = 1, s2 = 1.5,
         rho = -0.3, tol = 1e-6),
    list(n1 = 230, n2 = 230, d1 = 0.3, d2 = 0.3, s1 = 1, s2 = 1.0,
         rho = 0.8, tol = 1e-5)
  )
  for (cs in cases) {
    args <- list(n1 = cs$n1, n2 = cs$n2,
                 delta1 = cs$d1, delta2 = cs$d2, sd1 = cs$s1, sd2 = cs$s2,
                 rho = cs$rho, alpha = 0.025)
    pw_uni <- do.call(twocpcont_power, c(args, method = "univariate"))
    pw_biv <- do.call(twocpcont_power, c(args, method = "bivariate"))
    expect_lt(abs(pw_uni$powerCoprimary - pw_biv$powerCoprimary), cs$tol)
    expect_equal(pw_uni$power1, pw_biv$power1)
    expect_equal(pw_uni$power2, pw_biv$power2)
  }
})

test_that("twocpcont_power: marginal powers follow the normal formula", {
  n1 <- 120; n2 <- 80
  pw <- twocpcont_power(n1 = n1, n2 = n2,
                        delta1 = 0.7, delta2 = 0.5, sd1 = 1, sd2 = 1.5,
                        rho = 0.3, alpha = 0.025, method = "bivariate")
  se_inv  <- 1 / sqrt(1 / n1 + 1 / n2)
  z_alpha <- qnorm(1 - 0.025)
  expect_equal(pw$power1, pnorm(0.7 / 1 * se_inv - z_alpha))
  expect_equal(pw$power2, pnorm(0.5 / 1.5 * se_inv - z_alpha))
})

test_that("twocpcont_power at rho = 0 equals the product of marginals", {
  args <- list(n1 = 100, n2 = 100,
               delta1 = 0.7, delta2 = 0.5, sd1 = 1, sd2 = 1,
               rho = 0, alpha = 0.025)
  pw_uni <- do.call(twocpcont_power, c(args, method = "univariate"))
  pw_biv <- do.call(twocpcont_power, c(args, method = "bivariate"))
  expect_equal(pw_uni$powerCoprimary, pw_uni$power1 * pw_uni$power2)
  expect_equal(pw_biv$powerCoprimary, pw_biv$power1 * pw_biv$power2,
               tolerance = 1e-8)
})

test_that("twocpcont_power keeps the input order of the endpoints", {
  for (m in c("univariate", "bivariate")) {
    pw_a <- twocpcont_power(n1 = 100, n2 = 100,
                            delta1 = 0.7, delta2 = 0.5, sd1 = 1, sd2 = 1.5,
                            rho = 0.5, method = m)
    pw_b <- twocpcont_power(n1 = 100, n2 = 100,
                            delta1 = 0.5, delta2 = 0.7, sd1 = 1.5, sd2 = 1,
                            rho = 0.5, method = m)
    expect_equal(pw_a$power1, pw_b$power2)
    expect_equal(pw_a$power2, pw_b$power1)
    expect_equal(pw_a$powerCoprimary, pw_b$powerCoprimary, tolerance = 1e-10)
  }
})

test_that("twocpcont_power: co-primary power does not exceed the marginals", {
  for (rho in c(-0.5, 0, 0.3, 0.8)) {
    pw <- twocpcont_power(n1 = 100, n2 = 100,
                          delta1 = 0.4, delta2 = 0.5, sd1 = 1, sd2 = 1,
                          rho = rho, method = "bivariate")
    expect_lte(pw$powerCoprimary, min(pw$power1, pw$power2))
  }
})

test_that("twocpcont_power increases in rho and in the sample size", {
  for (m in c("univariate", "bivariate")) {
    pw_rho <- sapply(c(-0.5, 0, 0.3, 0.5, 0.8), function(rr) {
      twocpcont_power(n1 = 100, n2 = 100,
                      delta1 = 0.4, delta2 = 0.4, sd1 = 1, sd2 = 1,
                      rho = rr, method = m)$powerCoprimary
    })
    expect_true(all(diff(pw_rho) > 0))
    pw_n <- sapply(c(50, 100, 150, 200), function(nn) {
      twocpcont_power(n1 = nn, n2 = nn,
                      delta1 = 0.4, delta2 = 0.4, sd1 = 1, sd2 = 1,
                      rho = 0.5, method = m)$powerCoprimary
    })
    expect_true(all(diff(pw_n) > 0))
  }
})
