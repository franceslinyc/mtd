test_that("lagmat with k = 1 drops the last observation", {
  expect_equal(lagmat(1:5, 1), 1:4)
})

test_that("lagmat returns an (n - k) x k matrix of lagged values", {
  out <- lagmat(1:6, 2)
  expect_equal(dim(out), c(4, 2))
  expect_equal(unname(out), matrix(c(2:5, 1:4), ncol = 2))

  out <- lagmat(1:6, 3)
  expect_equal(dim(out), c(3, 3))
  expect_equal(unname(out), matrix(c(3:5, 2:4, 1:3), ncol = 3))
})

test_that("specRound pads to k decimal places", {
  expect_equal(specRound(3.14159, 2), "3.14")
  expect_equal(specRound(2, 2), "2.00")
})

test_that("specCI formats quantile limits", {
  expect_equal(specCI(0:100), "(2.50, 97.50)")
})

test_that("summarize.numeric returns point estimate and interval", {
  x <- 0:100
  expect_equal(summarize(x, "mean"), "50.00 (2.50, 97.50)")
  expect_equal(summarize(x, "median"), "50.00 (2.50, 97.50)")
  expect_length(summarize(x, "both"), 2)
})

test_that("summarize.numeric is case-insensitive for point_est", {
  expect_equal(summarize(0:100, "MEAN"), summarize(0:100, "mean"))
})

test_that("summarize.numeric rejects invalid arguments", {
  x <- 0:100
  expect_error(summarize(x, "mode"), "not a valid option")
  expect_error(summarize(x, "mean", lower_prob = 0), "lower_prob must be between 0 and 1")
  expect_error(summarize(x, "mean", upper_prob = 1), "upper_prob must be between 0 and 1")
  expect_error(summarize(x, "mean", lower_prob = 0.9, upper_prob = 0.1),
               "lower_prob must be smaller than upper_prob")
})

test_that("summarize.matrix summarizes each row", {
  x <- matrix(rep(0:100, 2), nrow = 2, byrow = TRUE)
  expect_equal(summarize(x, "mean"), rep("50.00 (2.50, 97.50)", 2))
})

test_that("summarize.matrix rejects invalid arguments", {
  x <- matrix(rnorm(20), nrow = 2)
  expect_error(summarize(x, "mode"), "not a valid option")
  expect_error(summarize(x, "mean", lower_prob = 0.9, upper_prob = 0.1),
               "lower_prob must be smaller than upper_prob")
})
