# Tests for the internal function .correction_term()

test_that(".correction_term equals the root with the smaller absolute value", {
  C0  <- c(-0.01, 0.02, -0.03)
  P_1 <- c(0.30, 0.25, 0.40)
  P_2 <- c(-0.20, 0.15, 0.30)
  for (i in seq_along(C0)) {
    d_quad <- (-P_1[i] + sqrt(P_1[i] ^ 2 - 2 * P_2[i] * C0[i])) / P_2[i]
    d <- .correction_term(C0[i], P_1[i], P_2[i])
    expect_equal(d, d_quad, tolerance = 1e-12)
    # The value solves the quadratic Taylor equation
    expect_lt(abs(C0[i] + P_1[i] * d + P_2[i] * d ^ 2 / 2), 1e-14)
  }
})

test_that(".correction_term remains defined when P_2 is zero", {
  expect_equal(.correction_term(0.01, 0.3, 0), -0.01 / 0.3,
               tolerance = 1e-14)
  expect_equal(.correction_term(0, 0.3, 0.2), 0)
})

test_that(".correction_term uses the first-order correction for a negative discriminant", {
  # P_1^2 - 2 P_2 C0 = 0.01 - 2 < 0
  expect_warning(d <- .correction_term(1, 0.1, 1), "discriminant")
  expect_equal(d, -1 / 0.1)
})
