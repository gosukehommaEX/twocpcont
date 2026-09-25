# Tests for print.twocpcont_ss()

test_that("print.twocpcont_ss shows the three sections and the method", {
  ss_uni <- twocpcont_ss(delta1 = 1, delta2 = 1, sd1 = 2, sd2 = 2,
                         rho = 0.5, method = "univariate")
  ss_biv <- twocpcont_ss(delta1 = 1, delta2 = 1, sd1 = 2, sd2 = 2,
                         rho = 0.5, method = "bivariate")
  out_uni <- capture.output(print(ss_uni))
  out_biv <- capture.output(print(ss_biv))
  expect_true(any(grepl("Sample size calculation", out_uni)))
  for (s in c("[Design]", "[Computation]", "[Results]")) {
    expect_true(any(grepl(s, out_uni, fixed = TRUE)))
  }
  expect_true(any(grepl("Noniterative closed-form (proposed)", out_uni,
                        fixed = TRUE)))
  expect_true(any(grepl("Bivariate normal (conventional)", out_biv,
                        fixed = TRUE)))
  expect_true(any(grepl(paste0(": ", ss_uni$N, "$"), out_uni)))
})

test_that("print.twocpcont_ss returns its argument invisibly", {
  ss <- twocpcont_ss(delta1 = 1, delta2 = 1, sd1 = 2, sd2 = 2,
                     rho = 0.5, method = "univariate")
  capture.output(res <- withVisible(print(ss)))
  expect_false(res$visible)
  expect_identical(res$value, ss)
})

test_that("print.twocpcont_ss accepts the digits argument", {
  ss <- twocpcont_ss(delta1 = 0.7, delta2 = 0.5, sd1 = 1, sd2 = 1.5,
                     rho = 0.35, method = "univariate")
  expect_output(print(ss, digits = 2), "0.35")
})
