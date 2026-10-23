gaussian_setup <- function() {
  set.seed(1)
  param <- list(rho = c(0.8, 0.5, 0.2), mu = 2, sigma2 = 1)
  obs <- rMTD(3, c(0.5, 0.3, 0.2), "gaussian", param, size = 200)$y

  list(
    obs = obs,
    prior = list(mu_0 = 0, sigma2_0 = 10, u_0 = 2, v_0 = 2, alpha = 1),
    tuning = list(step_size = 0.1),
    starting = list(rho = rep(0.5, 3), sigma2 = 1),
    mcmc_param = list(niter = 200, nburn = 100, nthin = 2)
  )
}

test_that("tsMTD fits a Gaussian MTD model with the expected output shape", {
  s <- gaussian_setup()
  set.seed(2)
  fit <- tsMTD(s$obs, 3, "gaussian", "sb", s$prior, s$tuning, s$starting, s$mcmc_param)

  expect_named(fit, c("weight", "rho", "mu", "sigma2"))
  expect_equal(dim(fit$weight), c(3, 50))
  expect_equal(dim(fit$rho), c(3, 50))
  expect_length(fit$mu, 50)
  expect_length(fit$sigma2, 50)
})

test_that("tsMTD Gaussian samples are valid", {
  s <- gaussian_setup()
  set.seed(2)
  fit <- tsMTD(s$obs, 3, "gaussian", "sb", s$prior, s$tuning, s$starting, s$mcmc_param)

  expect_false(anyNA(unlist(fit)))
  expect_true(all(fit$sigma2 > 0))
  expect_equal(colSums(fit$weight), rep(1, 50))
})

test_that("tsMTD Gaussian fit is reproducible with a seed", {
  s <- gaussian_setup()
  set.seed(2)
  a <- tsMTD(s$obs, 3, "gaussian", "sb", s$prior, s$tuning, s$starting, s$mcmc_param)
  set.seed(2)
  b <- tsMTD(s$obs, 3, "gaussian", "sb", s$prior, s$tuning, s$starting, s$mcmc_param)
  expect_identical(a, b)
})

test_that("tsMTD validates its basic arguments", {
  s <- gaussian_setup()
  expect_error(
    tsMTD(s$obs[1:3], 3, "gaussian", "sb", s$prior, s$tuning, s$starting, s$mcmc_param),
    "data size must be greater"
  )
  expect_error(
    tsMTD(s$obs, 3, "foo", "sb", s$prior, s$tuning, s$starting, s$mcmc_param),
    "not a valid option"
  )
  s$prior$mu_0 <- NULL
  expect_error(
    tsMTD(s$obs, 3, "gaussian", "sb", s$prior, s$tuning, s$starting, s$mcmc_param),
    "mu_0 must be specified"
  )
})
