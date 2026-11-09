# Tests for print.twocpcont_power()

test_that("print.twocpcont_power shows the three sections and the method", {
  pw_uni <- twocpcont_power(n1 = 100, n2 = 100, delta1 = 1, delta2 = 1,
                            sd1 = 2, sd2 = 2, rho = 0.5,
                            method = "univariate")
  pw_biv <- twocpcont_power(n1 = 100, n2 = 100, delta1 = 1, delta2 = 1,
                            sd1 = 2, sd2 = 2, rho = 0.5,
                            method = "bivariate")
  out_uni <- capture.output(print(pw_uni))
  out_biv <- capture.output(print(pw_biv))
  expect_true(any(grepl("Power calculation", out_uni)))
  for (s in c("[Design]", "[Computation]", "[Results]")) {
    expect_true(any(grepl(s, out_uni, fixed = TRUE)))
  }
  expect_true(any(grepl("Plackett-GL approx. (proposed)", out_uni,
                        fixed = TRUE)))
  expect_true(any(grepl("Bivariate normal (conventional)", out_biv,
                        fixed = TRUE)))
})

test_that("print.twocpcont_power returns its argument invisibly", {
  pw <- twocpcont_power(n1 = 100, n2 = 100, delta1 = 1, delta2 = 1,
                        sd1 = 2, sd2 = 2, rho = 0.5,
                        method = "univariate")
  capture.output(res <- withVisible(print(pw)))
  expect_false(res$visible)
  expect_identical(res$value, pw)
})

test_that("print.twocpcont_power accepts the digits argument", {
  pw <- twocpcont_power(n1 = 100, n2 = 100, delta1 = 1, delta2 = 1,
                        sd1 = 2, sd2 = 2, rho = 0.5,
                        method = "univariate")
  expect_output(print(pw, digits = 2), "Co-primary power")
})

test_that("print.twocpcont_power labels delta1 and delta2 as mean differences", {
  pw <- twocpcont_power(n1 = 100, n2 = 100, delta1 = 1.8, delta2 = 3.1,
                        sd1 = 9, sd2 = 12, rho = 0.5)
  expect_output(print(pw), "Mean difference \\(d1, d2\\) : 1.800, 3.100")
})
