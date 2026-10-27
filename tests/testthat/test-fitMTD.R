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
  set.seed(10)
  list(
    obs = rgamma(80, shape = 2, scale = 1.5),
    prior = list(alpha = 1, u_alp = 2, v_alp = 1, u_be = 2, v_be = 1),
    tuning = list(step_size = 0.15, upper_alp = 8, upper_be = 4),
    starting = list(rho = c(0.2, 0.35), alp = 2, be = 1.5),
    mcmc_param = list(niter = 12, nburn = 6, nthin = 2)
  )
}

test_that("tsMTD runs the Gamma copula sampler", {
  s <- gamma_copula_setup()
  set.seed(11)
  fit <- tsMTD(s$obs, 2, "gamma", "sb", s$prior, s$tuning, s$starting,
               s$mcmc_param, cop_family = 1)

  expect_named(fit, c("weight", "rho", "shape1", "shape2"))
  expect_equal(dim(fit$weight), c(2, 3))
  expect_equal(dim(fit$rho), c(2, 3))
  expect_true(all(is.finite(unlist(fit))))
  expect_equal(colSums(fit$weight), rep(1, 3))
  expect_true(all(fit$shape1 > 0))
  expect_true(all(fit$shape2 > 0))
})

test_that("tsMTD runs the ZI Gamma copula sampler with fixed eps", {
  s <- gamma_copula_setup()
  set.seed(12)
  s$obs <- ifelse(runif(80) < 0.25, runif(80, 0, 0.3), rgamma(80, 2 / 1.5, scale = 1.5) + 0.3)
  s$starting$eps <- 0.3
  s$starting$Pi <- 0.25
  set.seed(13)
  fit <- tsMTD(s$obs, 2, "zigamma", "sb", s$prior, s$tuning, s$starting,
               s$mcmc_param, cop_family = 1)

  expect_named(fit, c("weight", "rho", "shape1", "shape2", "Pi", "eps"))
  expect_equal(dim(fit$weight), c(2, 3))
  expect_equal(dim(fit$rho), c(2, 3))
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
    tsMTD(s$obs, 2, "gamma", "sb", s$prior, s$tuning, s$starting, s$mcmc_param),
    "cop_family must be specified"
  )
  expect_error(
    tsMTD(s$obs, 2, "zigamma", "sb", s$prior, s$tuning, s$starting, s$mcmc_param,
          cop_family = 2),
    "only cop_family = 1"
  )
})
