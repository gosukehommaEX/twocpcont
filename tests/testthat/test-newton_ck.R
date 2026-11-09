# Tests for newton_ck()

test_that("newton_ck returns a one-row data frame with the documented columns", {
  ck <- newton_ck(kappa = 1.2, rho = 0.3)
  expect_s3_class(ck, "data.frame")
  expect_equal(nrow(ck), 1L)
  expect_named(ck, c("kappa", "rho", "alpha", "beta", "C", "lambda",
                     "iterations", "n_eval", "converged"))
  expect_equal(ck$lambda, ck$C + qnorm(1 - 0.025))
  expect_equal(c(ck$alpha, ck$beta), c(0.025, 0.2))
})

test_that("newton_ck reproduces the example of Sugimoto et al. (2012)", {
  # Supporting Information, Appendix A.2.1: gamma = 8/7, rho = 0.5,
  # alpha = 0.025 and 1 - beta = 0.8 give C = 0.9988124, and with
  # delta = (0.4, 0.35) and p = 1, n = 2 (C + z_0.025)^2 / 0.35^2 = 142.93
  ck <- newton_ck(kappa = 8 / 7, rho = 0.5, alpha = 0.025, beta = 0.2)
  expect_true(ck$converged)
  expect_equal(round(ck$C, 7), 0.9988124, tolerance = 1e-12)
  expect_equal(round(2 * ck$lambda ^ 2 / 0.35 ^ 2, 2), 142.93,
               tolerance = 1e-12)
})

test_that("newton_ck meets its stopping rule near the exact root", {
  z_a <- qnorm(1 - 0.025)
  for (k in c(1, 1.3, 2)) {
    for (rho in c(0, 0.5, 0.8)) {
      for (beta in c(0.2, 0.1)) {
        ck <- newton_ck(kappa = k, rho = rho, alpha = 0.025, beta = beta)
        pw <- pbivnorm::pbivnorm(k * ck$C + z_a * (k - 1), ck$C, rho = rho)
        C_exact <- uniroot(function(cc) {
          pbivnorm::pbivnorm(k * cc + z_a * (k - 1), cc, rho = rho) -
            (1 - beta)
        }, lower = -2, upper = 5, tol = 1e-13)$root
        expect_true(ck$converged)
        # |G| < 1e-5 and G <= 0: the power is at least the target
        expect_gte(pw, 1 - beta)
        expect_lt(pw - (1 - beta), 1e-5)
        expect_gte(ck$C, C_exact - 1e-10)
        expect_lt(ck$C - C_exact, 1e-4)
      }
    }
  }
})

test_that("newton_ck counts the evaluations of the bivariate normal CDF", {
  ck <- newton_ck(kappa = 1.3, rho = 0.5, alpha = 0.025, beta = 0.2)
  # Three evaluations per update, one for the final check and at least
  # two for the initial value (the two ends of the uniroot() interval)
  expect_gte(ck$n_eval - 3L * ck$iterations - 1L, 2L)
  expect_gte(ck$iterations, 1L)
})

test_that("newton_ck reports that the stopping rule was not met", {
  # The initial value is the solution for kappa = 1, so at kappa = 2 the
  # stopping rule is not met without an update
  ck <- newton_ck(kappa = 2, rho = 0.5, max_iter = 0)
  expect_false(ck$converged)
  expect_equal(ck$iterations, 0L)
})

test_that("newton_ck rejects invalid arguments", {
  expect_error(newton_ck(kappa = 0.9, rho = 0.5), "kappa")
  expect_error(newton_ck(kappa = 1.2, rho = 1), "rho")
  expect_error(newton_ck(kappa = 1.2, rho = -1), "rho")
})
