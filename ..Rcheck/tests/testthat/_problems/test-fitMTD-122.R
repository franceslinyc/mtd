# Extracted from test-fitMTD.R:122

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "mtd", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
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

# test -------------------------------------------------------------------------
s <- zigamma_copula_setup()
fit <- tsMTD(s$obs, 5, "zigamma", "sb", s$prior, s$tuning, s$starting,
               s$mcmc_param, cop_family = 1)
expect_named(fit, c("weight", "rho", "shape1", "shape2", "Pi", "eps"))
