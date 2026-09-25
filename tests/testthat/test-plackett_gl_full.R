# Tests for plackett_gl_full()

# Helper: arguments u_w and u_s as functions of lambda_w
u_args <- function(lam_w, kappa, alpha = 0.025) {
  z_alpha <- qnorm(1 - alpha)
  list(u_w = lam_w - z_alpha, u_s = kappa * lam_w - z_alpha)
}

# Helper: I_GL, I_GL_1 or I_GL_2 as a function of lambda_w
pk_at <- function(lam_w, kappa, rho, component, gl_nodes = 5L) {
  u <- u_args(lam_w, kappa)
  plackett_gl_full(u$u_w, u$u_s, rho, kappa, gl_nodes = gl_nodes)[[component]]
}

test_that("plackett_gl_full returns a list with three components", {
  pk <- plackett_gl_full(u_w = 0.84, u_s = 0.84, rho = 0.3, kappa = 1)
  expect_named(pk, c("I_GL", "I_GL_1", "I_GL_2"))
})

test_that("plackett_gl_full vanishes at rho = 0", {
  pk <- plackett_gl_full(u_w = 0.84, u_s = 1.2, rho = 0, kappa = 1.4)
  expect_equal(pk, list(I_GL = 0, I_GL_1 = 0, I_GL_2 = 0))
})

test_that("I_GL matches the exact bivariate normal correction", {
  # Exact value: Phi_2(u_w, u_s; rho) - Phi(u_w) Phi(u_s)
  cases <- list(
    list(lam = 2.8, kappa = 1.0, rho = 0.3,  tol = 1e-8),
    list(lam = 2.8, kappa = 1.0, rho = 0.5,  tol = 1e-8),
    list(lam = 2.5, kappa = 1.4, rho = 0.5,  tol = 1e-8),
    list(lam = 2.8, kappa = 1.0, rho = -0.5, tol = 1e-8),
    list(lam = 2.8, kappa = 1.0, rho = 0.8,  tol = 1e-5),
    list(lam = 2.5, kappa = 2.0, rho = 0.8,  tol = 1e-5)
  )
  for (cs in cases) {
    u <- u_args(cs$lam, cs$kappa)
    exact <- pbivnorm::pbivnorm(u$u_w, u$u_s, rho = cs$rho) -
      pnorm(u$u_w) * pnorm(u$u_s)
    pk <- plackett_gl_full(u$u_w, u$u_s, cs$rho, cs$kappa)
    expect_lt(abs(pk$I_GL - exact), cs$tol)
  }
})

test_that("More quadrature nodes reduce the error at rho = 0.8", {
  u <- u_args(2.8, 1)
  exact <- pbivnorm::pbivnorm(u$u_w, u$u_s, rho = 0.8) -
    pnorm(u$u_w) * pnorm(u$u_s)
  err5  <- abs(plackett_gl_full(u$u_w, u$u_s, 0.8, 1, gl_nodes = 5)$I_GL -
                 exact)
  err10 <- abs(plackett_gl_full(u$u_w, u$u_s, 0.8, 1, gl_nodes = 10)$I_GL -
                 exact)
  expect_lt(err10, err5)
})

test_that("I_GL_1 and I_GL_2 match finite differences in lambda_w", {
  h <- 1e-4
  cases <- list(
    list(lam = 2.8, kappa = 1.0, rho = 0.5),
    list(lam = 2.5, kappa = 1.4, rho = 0.3),
    list(lam = 2.5, kappa = 2.0, rho = 0.8),
    list(lam = 2.8, kappa = 1.2, rho = -0.5)
  )
  for (cs in cases) {
    d1_fd <- (pk_at(cs$lam + h, cs$kappa, cs$rho, "I_GL") -
                pk_at(cs$lam - h, cs$kappa, cs$rho, "I_GL")) / (2 * h)
    d2_fd <- (pk_at(cs$lam + h, cs$kappa, cs$rho, "I_GL_1") -
                pk_at(cs$lam - h, cs$kappa, cs$rho, "I_GL_1")) / (2 * h)
    expect_equal(pk_at(cs$lam, cs$kappa, cs$rho, "I_GL_1"), d1_fd,
                 tolerance = 1e-6)
    expect_equal(pk_at(cs$lam, cs$kappa, cs$rho, "I_GL_2"), d2_fd,
                 tolerance = 1e-6)
  }
})

test_that("plackett_gl_full returns the analytic limit at rho = 1", {
  u <- u_args(2.5, 1.3)
  pk <- plackett_gl_full(u$u_w, u$u_s, rho = 1, kappa = 1.3)
  expect_equal(pk$I_GL, pnorm(u$u_w) * (1 - pnorm(u$u_s)))
})

test_that("Derivatives at rho = 1 match finite differences", {
  h <- 1e-4
  for (kappa in c(1, 1.3, 2)) {
    d1_fd <- (pk_at(2.5 + h, kappa, 1, "I_GL") -
                pk_at(2.5 - h, kappa, 1, "I_GL")) / (2 * h)
    d2_fd <- (pk_at(2.5 + h, kappa, 1, "I_GL_1") -
                pk_at(2.5 - h, kappa, 1, "I_GL_1")) / (2 * h)
    expect_equal(pk_at(2.5, kappa, 1, "I_GL_1"), d1_fd, tolerance = 1e-6)
    expect_equal(pk_at(2.5, kappa, 1, "I_GL_2"), d2_fd, tolerance = 1e-6)
  }
})

test_that("plackett_gl_full returns the analytic limits at rho = -1", {
  # lam_w = 0.5 gives Phi(u_w) + Phi(u_s) <= 1; lam_w = 3 gives > 1
  h <- 1e-4
  for (lam in c(0.5, 3)) {
    u  <- u_args(lam, 1.2)
    pk <- plackett_gl_full(u$u_w, u$u_s, rho = -1, kappa = 1.2)
    Phi_w <- pnorm(u$u_w)
    Phi_s <- pnorm(u$u_s)
    expect_equal(pk$I_GL, max(0, Phi_w + Phi_s - 1) - Phi_w * Phi_s)
    d1_fd <- (pk_at(lam + h, 1.2, -1, "I_GL") -
                pk_at(lam - h, 1.2, -1, "I_GL")) / (2 * h)
    d2_fd <- (pk_at(lam + h, 1.2, -1, "I_GL_1") -
                pk_at(lam - h, 1.2, -1, "I_GL_1")) / (2 * h)
    expect_equal(pk$I_GL_1, d1_fd, tolerance = 1e-6)
    expect_equal(pk$I_GL_2, d2_fd, tolerance = 1e-6)
  }
})
