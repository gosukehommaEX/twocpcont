# Tests for GL_nodes_and_weights()

test_that("GL_nodes_and_weights returns the requested number of nodes", {
  for (k in c(1, 2, 5, 7, 10)) {
    gl <- GL_nodes_and_weights(gl_nodes = k)
    expect_named(gl, c("x", "w"))
    expect_length(gl$x, k)
    expect_length(gl$w, k)
  }
})

test_that("GL_nodes_and_weights uses 5 nodes by default", {
  expect_equal(GL_nodes_and_weights(), GL_nodes_and_weights(gl_nodes = 5))
})

test_that("GL_nodes_and_weights handles the 1-point rule", {
  expect_equal(GL_nodes_and_weights(gl_nodes = 1), list(x = 0, w = 2))
})

test_that("GL_nodes_and_weights reproduces the known 2-point rule", {
  gl <- GL_nodes_and_weights(gl_nodes = 2)
  expect_equal(gl$x, c(-1, 1) / sqrt(3), tolerance = 1e-12)
  expect_equal(gl$w, c(1, 1), tolerance = 1e-12)
})

test_that("GL_nodes_and_weights: weights sum to 2 on [-1, 1]", {
  for (k in c(3, 5, 7, 10)) {
    gl <- GL_nodes_and_weights(gl_nodes = k)
    expect_equal(sum(gl$w), 2, tolerance = 1e-10)
  }
})

test_that("GL_nodes_and_weights: nodes lie in (-1, 1) in increasing order", {
  for (k in c(2, 5, 10)) {
    gl <- GL_nodes_and_weights(gl_nodes = k)
    expect_true(all(gl$x > -1 & gl$x < 1))
    expect_true(all(diff(gl$x) > 0))
    expect_true(all(gl$w > 0))
  }
})

test_that("GL_nodes_and_weights: nodes and weights are symmetric about 0", {
  for (k in c(4, 5, 10)) {
    gl <- GL_nodes_and_weights(gl_nodes = k)
    expect_equal(gl$x, -rev(gl$x), tolerance = 1e-12)
    expect_equal(gl$w, rev(gl$w), tolerance = 1e-12)
  }
})

test_that("GL_nodes_and_weights is exact up to degree 2 gl_nodes - 1", {
  # Integral of x^m over [-1, 1]: 2 / (m + 1) for even m, 0 for odd m
  for (k in 1:8) {
    gl <- GL_nodes_and_weights(gl_nodes = k)
    for (m in 0:(2 * k - 1)) {
      exact <- if (m %% 2 == 0) 2 / (m + 1) else 0
      expect_equal(sum(gl$w * gl$x ^ m), exact, tolerance = 1e-10)
    }
  }
})

test_that("GL_nodes_and_weights: 5-point rule integrates x^8 correctly", {
  gl <- GL_nodes_and_weights(gl_nodes = 5)
  expect_equal(sum(gl$w * gl$x ^ 8), 2 / 9, tolerance = 1e-10)
})

test_that("GL_nodes_and_weights rejects a non-positive number of nodes", {
  expect_error(GL_nodes_and_weights(gl_nodes = 0), "positive integer")
  expect_error(GL_nodes_and_weights(gl_nodes = -3), "positive integer")
})
