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

gamma_copula_setup <- function() {
  set.seed(24)
  J <- 5
  expDecay <- function(x, la) exp(-la * x)
  true_weight <- expDecay(1:J, 1)
  true_weight <- true_weight / sum(true_weight)
  param <- list(rho = c(0.7, 0.5, 0.3, 0.1, 0.1), shape = 7, rate = 1)
  obs <- rMTD(J, true_weight, "gamma", param, size = 80, trun = FALSE)$y

  list(
    obs = obs,
    prior = list(alpha = 1, u_alp = 2, v_alp = 1, u_be = 2, v_be = 1),
    tuning = list(step_size = 0.15, upper_alp = 12, upper_be = 4),
    starting = list(rho = rep(0.3, J), alp = 7, be = 1),
    mcmc_param = list(niter = 12, nburn = 6, nthin = 2)
  )
}

zigamma_copula_setup <- function() {
  set.seed(24)
  J <- 5
  expDecay <- function(x, la) exp(-la * x)
  true_weight <- expDecay(1:J, 1)
  true_weight <- true_weight / sum(true_weight)
  param <- list(rho = c(0.7, 0.5, 0.3, 0.1, 0.1),
                mu = 7, beta = 1, epsilon = 0.3, Pi = 0.25)
  obs <- suppressWarnings(
    rMTD(J, true_weight, "zigamma", param, size = 80, trun = FALSE)$y
  )

  list(
    obs = obs,
    prior = list(alpha = 1, u_alp = 2, v_alp = 1, u_be = 2, v_be = 1),
    tuning = list(step_size = 0.15, upper_alp = 12, upper_be = 4),
    starting = list(rho = rep(0.3, J), alp = 7, be = 1, eps = param$epsilon, Pi = param$Pi),
    mcmc_param = list(niter = 12, nburn = 6, nthin = 2)
  )
}

test_that("tsMTD runs the Gamma copula sampler", {
  s <- gamma_copula_setup()
  set.seed(11)
  fit <- tsMTD(s$obs, 5, "gamma", "sb", s$prior, s$tuning, s$starting,
               s$mcmc_param, cop_family = 1)

  expect_named(fit, c("weight", "rho", "shape1", "shape2"))
  expect_equal(dim(fit$weight), c(5, 3))
  expect_equal(dim(fit$rho), c(5, 3))
  expect_true(all(is.finite(unlist(fit))))
  expect_equal(colSums(fit$weight), rep(1, 3))
  expect_true(all(fit$shape1 > 0))
  expect_true(all(fit$shape2 > 0))
})

test_that("tsMTD runs the ZI Gamma copula sampler with fixed eps", {
  s <- zigamma_copula_setup()
  fit <- tsMTD(s$obs, 5, "zigamma", "sb", s$prior, s$tuning, s$starting,
               s$mcmc_param, cop_family = 1)

  expect_named(fit, c("weight", "rho", "shape1", "shape2", "Pi", "eps"))
  expect_equal(dim(fit$weight), c(5, 3))
  expect_equal(dim(fit$rho), c(5, 3))
  expect_true(all(is.finite(unlist(fit))))
  expect_equal(colSums(fit$weight), rep(1, 3))
  expect_true(all(fit$shape1 > 0))
  expect_true(all(fit$shape2 > 0))
  expect_true(all(fit$Pi > 0 & fit$Pi < 1))
  expect_equal(as.numeric(fit$eps), rep(0.3, 3))
})

test_that("tsMTD requires supported copula settings for Gamma families", {
  s <- gamma_copula_setup()
  expect_error(
    tsMTD(s$obs, 5, "gamma", "sb", s$prior, s$tuning, s$starting, s$mcmc_param),
    "cop_family must be specified"
  )
  expect_error(
    tsMTD(s$obs, 5, "zigamma", "sb", s$prior, s$tuning, s$starting, s$mcmc_param,
          cop_family = 2),
    "only cop_family = 1"
  )
})
