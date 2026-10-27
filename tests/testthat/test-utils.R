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

test_that("psgamma is a shifted Gamma CDF", {
  x <- c(0.1, 0.5, 1, 2)
  expect_equal(psgamma(x, mu = 4, beta = 2, epsilon = 0.5),
               c(0, 0, pgamma(0.5, shape = 2, scale = 2),
                 pgamma(1.5, shape = 2, scale = 2)))
})

test_that("pzig combines uniform and shifted Gamma CDFs", {
  x <- c(0, 0.25, 0.5, 1.5)
  mu <- 4
  beta <- 2
  epsilon <- 0.5
  Pi <- 0.3
  expected <- Pi * punif(x, 0, epsilon) +
    (1 - Pi) * ifelse(x >= epsilon,
                      pgamma(x - epsilon, shape = mu / beta, scale = beta),
                      0)
  expect_equal(pzig(x, mu, beta, epsilon, Pi), expected)
})

test_that("qsgamma shifts Gamma quantiles by epsilon", {
  probs <- c(0.1, 0.5, 0.9)
  expect_equal(qsgamma(probs, mu = 4, beta = 2, epsilon = 0.5),
               qgamma(probs, shape = 2, scale = 2) + 0.5)
})

test_that("qzig inverts pzig across both mixture components", {
  mu <- 4
  beta <- 2
  epsilon <- 0.5
  Pi <- 0.3
  x <- c(0.05, 0.2, 0.45, 0.5, 0.8, 2)
  result <- suppressWarnings(qzig(pzig(x, mu, beta, epsilon, Pi),
                                  mu, beta, epsilon, Pi))
  expect_equal(result, x, tolerance = 1e-8)
})

test_that("rsgamma produces shifted positive Gamma draws reproducibly", {
  set.seed(24)
  x <- rsgamma(100, mu = 4, beta = 2, epsilon = 0.5)
  set.seed(24)
  y <- rsgamma(100, mu = 4, beta = 2, epsilon = 0.5)

  expect_identical(x, y)
  expect_length(x, 100)
  expect_true(all(x >= 0.5))
})

test_that("rzig produces draws from both mixture components reproducibly", {
  set.seed(24)
  x <- rzig(1000, mu = 4, beta = 2, epsilon = 0.5, Pi = 0.3)
  set.seed(24)
  y <- rzig(1000, mu = 4, beta = 2, epsilon = 0.5, Pi = 0.3)

  expect_identical(x, y)
  expect_length(x, 1000)
  expect_true(all(x >= 0))
  expect_gt(sum(x < 0.5), 200)
  expect_lt(sum(x < 0.5), 400)
  expect_true(any(x > 0.5))
})

test_that("R simulation helpers do not mask the Rcpp scalar pzig and qzig exports", {
  expect_length(pzig(0.2, 4, 2, 0.5, 0.3), 1)
  expect_length(qzig(0.5, 4, 2, 0.5, 0.3), 1)
  expect_true(is.finite(psgamma(1, 4, 2, 0.5)))
})
