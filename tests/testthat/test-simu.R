w <- c(0.5, 0.5)

test_that("rMTD validates the model order", {
  expect_error(rMTD(2.5, w, "poisson", list(), 10), "must be an integer")
  expect_error(rMTD(1, w, "poisson", list(), 10), "greater than two")
})

test_that("rMTD validates the weights", {
  expect_error(rMTD(2, 1, "poisson", list(), 10), "at least two")
  expect_error(rMTD(3, w, "poisson", list(), 10), "same as the mtdorder")
  expect_error(rMTD(2, c(1, 0), "poisson", list(), 10), "positive vector")
})

test_that("rMTD validates the family", {
  expect_error(rMTD(2, w, param = list(), size = 10), "family must be specified")
  expect_error(rMTD(2, w, "foo", list(), 10), "not a valid option")
})

test_that("rMTD validates truncation settings", {
  expect_error(rMTD(2, w, "poisson", list(), size = 2, trun = TRUE),
               "greater than the model order")
  expect_error(rMTD(2, w, "poisson", list(), size = 10, trun = TRUE),
               "initial values must be specified")
  expect_error(rMTD(2, w, "poisson", list(), size = 10, trun = TRUE, init_data = 1),
               "same as mtdorder")
})

test_that("rMTD requires family parameters", {
  expect_error(rMTD(2, w, "gaussian", list(), 10), "rho must be specified")
})

test_that("rMTD simulates Gaussian data of the requested size", {
  param <- list(rho = c(0.8, 0.5, 0.2), mu = 2, sigma2 = 1)
  set.seed(1)
  out <- rMTD(3, c(0.5, 0.3, 0.2), "gaussian", param, size = 100)

  expect_named(out, c("y", "label", "lags"))
  expect_length(out$y, 100)
  expect_false(anyNA(out$y))
  expect_true(all(out$label[4:100] %in% 1:3))
})

test_that("rMTD Gaussian simulation is reproducible with a seed", {
  param <- list(rho = c(0.8, 0.5, 0.2), mu = 2, sigma2 = 1)
  set.seed(1)
  a <- rMTD(3, c(0.5, 0.3, 0.2), "gaussian", param, size = 50)
  set.seed(1)
  b <- rMTD(3, c(0.5, 0.3, 0.2), "gaussian", param, size = 50)
  expect_identical(a, b)
})

test_that("rMTD Gaussian simulation keeps init_data when trun = TRUE", {
  param <- list(rho = c(0.8, 0.5, 0.2), mu = 2, sigma2 = 1)
  set.seed(1)
  out <- rMTD(3, c(0.5, 0.3, 0.2), "gaussian", param, size = 50,
              trun = TRUE, init_data = c(1, 2, 3))

  expect_equal(as.numeric(out$y[1:3]), c(1, 2, 3))
  expect_length(out$y, 50)
  expect_length(out$label, 47)
})
