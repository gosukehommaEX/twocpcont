# Unit tests for the twocpcont package
# All tests are kept in a single file as requested.

# ============================================================
# (a) Mathematical equivalence of the two methods
# ============================================================
test_that("twocpcont_ss: univariate and bivariate methods agree", {
  # 8 representative patterns from the manuscript Section 3
  patterns <- list(
    list(d1 = 0.3, d2 = 0.3, s1 = 1, s2 = 1,   r = 1, rho = 0.0),
    list(d1 = 0.3, d2 = 0.3, s1 = 1, s2 = 1,   r = 1, rho = 0.5),
    list(d1 = 0.3, d2 = 0.3, s1 = 1, s2 = 1,   r = 1, rho = 0.8),
    list(d1 = 0.7, d2 = 0.5, s1 = 1, s2 = 1,   r = 1, rho = 0.5),
    list(d1 = 0.7, d2 = 0.5, s1 = 1, s2 = 1,   r = 2, rho = 0.5),
    list(d1 = 0.3, d2 = 0.3, s1 = 1, s2 = 1.5, r = 1, rho = 0.5),
    list(d1 = 0.7, d2 = 0.5, s1 = 1, s2 = 1.5, r = 1, rho = 0.5),
    list(d1 = 0.7, d2 = 0.5, s1 = 1, s2 = 1.5, r = 2, rho = 0.8)
  )
  for (p in patterns) {
    ss_uni <- twocpcont_ss(delta1 = p$d1, delta2 = p$d2,
                           sd1 = p$s1, sd2 = p$s2,
                           rho = p$rho, r = p$r,
                           alpha = 0.025, beta = 0.2,
                           method = "univariate")
    ss_biv <- twocpcont_ss(delta1 = p$d1, delta2 = p$d2,
                           sd1 = p$s1, sd2 = p$s2,
                           rho = p$rho, r = p$r,
                           alpha = 0.025, beta = 0.2,
                           method = "bivariate")
    expect_equal(ss_uni$N, ss_biv$N)
  }
})

test_that("twocpcont_power: univariate and bivariate methods agree", {
  # Power evaluation at fixed n1, n2 should match across methods
  pw_uni <- twocpcont_power(n1 = 230, n2 = 230,
                            delta1 = 0.3, delta2 = 0.3,
                            sd1 = 1, sd2 = 1,
                            rho = 0.5, alpha = 0.025,
                            method = "univariate")
  pw_biv <- twocpcont_power(n1 = 230, n2 = 230,
                            delta1 = 0.3, delta2 = 0.3,
                            sd1 = 1, sd2 = 1,
                            rho = 0.5, alpha = 0.025,
                            method = "bivariate")
  expect_equal(pw_uni$power, pw_biv$power, tolerance = 1e-4)
})

# ============================================================
# (b) Known values from the manuscript
# ============================================================
test_that("r_max returns the values reported in Section 3.3", {
  # R_max(kappa = 1, 1 - beta = 0.80) = 23.85%
  expect_equal(r_max(kappa = 1, alpha = 0.025, beta = 0.2),
               0.2385, tolerance = 1e-3)
  # R_max(kappa = 1, 1 - beta = 0.70) = 28.63%
  expect_equal(r_max(kappa = 1, alpha = 0.025, beta = 0.3),
               0.2863, tolerance = 1e-3)
  # R_max(kappa = 1, 1 - beta = 0.90) = 18.57%
  expect_equal(r_max(kappa = 1, alpha = 0.025, beta = 0.1),
               0.1857, tolerance = 1e-3)
})

test_that("kappa_star returns the value reported in Section 2", {
  # kappa_star(eps = 0.01, 1 - beta = 0.80) = 1.614
  expect_equal(kappa_star(epsilon = 0.01, alpha = 0.025, beta = 0.2),
               1.614, tolerance = 1e-3)
  # kappa_star(eps = 0.05, 1 - beta = 0.80) = 1.369
  expect_equal(kappa_star(epsilon = 0.05, alpha = 0.025, beta = 0.2),
               1.369, tolerance = 1e-3)
})

test_that("Doody 2014 example reproduces N = 1000 at rho = 0", {
  # EXPEDITION 1 design parameters
  ss <- twocpcont_ss(delta1 = 1.8, delta2 = 3.1,
                     sd1 = 9, sd2 = 12,
                     rho = 0, r = 1,
                     alpha = 0.025, beta = 0.13,
                     method = "univariate")
  expect_equal(ss$N, 1000L)
})

test_that("Doody 2014 example: N is monotonically non-increasing in rho", {
  N_vec <- sapply(c(0, 0.3, 0.5, 0.8), function(rr) {
    twocpcont_ss(delta1 = 1.8, delta2 = 3.1,
                 sd1 = 9, sd2 = 12,
                 rho = rr, r = 1,
                 alpha = 0.025, beta = 0.13,
                 method = "univariate")$N
  })
  expect_equal(N_vec, c(1000L, 990L, 980L, 960L))
  expect_true(all(diff(N_vec) <= 0))
})

# ============================================================
# (c) Edge cases (boundary inputs)
# ============================================================
test_that("twocpcont_ss runs at boundary rho = 0", {
  ss <- twocpcont_ss(delta1 = 0.3, delta2 = 0.3,
                     sd1 = 1, sd2 = 1,
                     rho = 0, r = 1,
                     alpha = 0.025, beta = 0.2,
                     method = "univariate")
  expect_true(is.numeric(ss$N))
  expect_true(ss$N > 0)
})

test_that("twocpcont_ss runs at boundary rho = 0.99", {
  ss <- twocpcont_ss(delta1 = 0.3, delta2 = 0.3,
                     sd1 = 1, sd2 = 1,
                     rho = 0.99, r = 1,
                     alpha = 0.025, beta = 0.2,
                     method = "univariate")
  expect_true(is.numeric(ss$N))
  expect_true(ss$N > 0)
})

test_that("twocpcont_ss handles allocation ratio r = 2", {
  ss <- twocpcont_ss(delta1 = 0.3, delta2 = 0.3,
                     sd1 = 1, sd2 = 1,
                     rho = 0.5, r = 2,
                     alpha = 0.025, beta = 0.2,
                     method = "univariate")
  expect_true(is.numeric(ss$N))
  expect_true(ss$N > 0)
})

# ============================================================
# (d) Class and print method tests
# ============================================================
test_that("twocpcont_ss returns an object of the correct class", {
  ss <- twocpcont_ss(delta1 = 1, delta2 = 1, sd1 = 2, sd2 = 2,
                     rho = 0.5, r = 1,
                     alpha = 0.025, beta = 0.2,
                     method = "univariate")
  expect_s3_class(ss, "twocpcont_ss")
})

test_that("twocpcont_power returns an object of the correct class", {
  pw <- twocpcont_power(n1 = 100, n2 = 100,
                        delta1 = 1, delta2 = 1, sd1 = 2, sd2 = 2,
                        rho = 0.5, alpha = 0.025,
                        method = "univariate")
  expect_s3_class(pw, "twocpcont_power")
})

test_that("print method for twocpcont_ss runs without error", {
  ss <- twocpcont_ss(delta1 = 1, delta2 = 1, sd1 = 2, sd2 = 2,
                     rho = 0.5, r = 1,
                     alpha = 0.025, beta = 0.2,
                     method = "univariate")
  expect_output(print(ss))
})

test_that("print method for twocpcont_power runs without error", {
  pw <- twocpcont_power(n1 = 100, n2 = 100,
                        delta1 = 1, delta2 = 1, sd1 = 2, sd2 = 2,
                        rho = 0.5, alpha = 0.025,
                        method = "univariate")
  expect_output(print(pw))
})

# ============================================================
# (e) GL_nodes_and_weights basic verification
# ============================================================
test_that("GL_nodes_and_weights returns the requested number of nodes", {
  gl <- GL_nodes_and_weights(gl_nodes = 5)
  expect_equal(length(gl$x), 5)
  expect_equal(length(gl$w), 5)
})

test_that("GL_nodes_and_weights: weights sum to 2 on [-1, 1]", {
  for (k in c(3, 5, 7, 10)) {
    gl <- GL_nodes_and_weights(gl_nodes = k)
    expect_equal(sum(gl$w), 2, tolerance = 1e-10)
  }
})

test_that("GL_nodes_and_weights: all nodes lie in [-1, 1]", {
  gl <- GL_nodes_and_weights(gl_nodes = 5)
  expect_true(all(gl$x >= -1))
  expect_true(all(gl$x <= 1))
})

test_that("GL_nodes_and_weights: nodes are in increasing order", {
  gl <- GL_nodes_and_weights(gl_nodes = 5)
  expect_true(all(diff(gl$x) > 0))
})

test_that("GL_nodes_and_weights: 5-point rule integrates x^8 correctly", {
  # Exact: integral of x^8 over [-1, 1] = 2 / 9
  gl <- GL_nodes_and_weights(gl_nodes = 5)
  expect_equal(sum(gl$w * gl$x ^ 8), 2 / 9, tolerance = 1e-10)
})
