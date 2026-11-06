test_that("PBV standard normal CDF wrapper matches pnorm", {
  x <- c(-2, -0.5, 0, 0.75, 2)
  expect_equal(vapply(x, pbv_rcpp_pnorm0, numeric(1)), pnorm(x), tolerance = 1e-12)
})

test_that("PBV bivariate normal CDF matches the orthant identity", {
  rho <- c(-0.9, -0.5, 0, 0.5, 0.9)
  expected <- 0.25 + asin(rho) / (2 * pi)
  actual <- vapply(rho, function(r) pbv_rcpp_pbvnorm0(0, 0, r), numeric(1))

  expect_equal(actual, expected, tolerance = 1e-7)
})

test_that("Gaussian copula marginal CDF is a product at zero correlation", {
  x <- c(-1, 0.25, 1.5)
  y <- c(0.5, -0.75, 2)
  mar_param <- c(0, 1, 0, 1)

  actual <- vapply(seq_along(x), function(i) {
    mtd:::pBiCopMar_cpp(x[i], y[i], 1L, 0, 3L, mar_param)
  }, numeric(1))
  expect_equal(actual, pnorm(x) * pnorm(y), tolerance = 1e-7)
})

test_that("vectorized Gaussian copula marginal CDF agrees with scalar calls", {
  x <- c(-1, 0.25, 1.5)
  y <- c(0.5, -0.75, 2)
  rho <- c(-0.3, 0, 0.6)
  mar_param <- matrix(rep(c(0, 1, 0, 1), length(x)), nrow = length(x), byrow = TRUE)

  vectorized <- mtd:::pBiCopMar_cpp2(x, y, 1L, rho, 3L, mar_param)
  scalar <- vapply(seq_along(x), function(i) {
    mtd:::pBiCopMar_cpp(x[i], y[i], 1L, rho[i], 3L, mar_param[i, ])
  }, numeric(1))

  expect_equal(as.numeric(vectorized), scalar, tolerance = 1e-12)
})

test_that("Gamma-copula randomized residual uses the marginal CDF at independence", {
  mtdorder <- 2
  obs <- c(0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8)
  shape <- 2
  rate <- 1
  weight <- matrix(c(0.4, 0.6), nrow = mtdorder, ncol = 1)
  rho <- matrix(0, nrow = mtdorder, ncol = 1)

  residual <- mtd:::rqrCopGammaMTD(obs, weight, shape, rate, 1L, rho, 2L)
  expected <- qnorm(pgamma(obs[(mtdorder + 1):length(obs)], shape = shape, rate = rate))

  expect_equal(dim(residual), c(length(obs) - mtdorder, 1))
  expect_equal(as.numeric(residual), expected, tolerance = 1e-7)
})