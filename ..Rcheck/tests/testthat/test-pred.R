prediction_copula_fit <- function(family) {
  set.seed(24)
  mtdorder <- 5
  expDecay <- function(x, la) exp(-la * x)
  weight <- expDecay(1:mtdorder, 1)
  weight <- weight / sum(weight)
  rho <- c(0.7, 0.5, 0.3, 0.1, 0.1)

  if (family == "gamma") {
    param <- list(rho = rho, shape = 7, rate = 1)
    obs <- rMTD(mtdorder, weight, "gamma", param, size = 80)$y
    prior <- list(alpha = 1, u_alp = 2, v_alp = 1, u_be = 2, v_be = 1)
    tuning <- list(step_size = 0.15, upper_alp = 12, upper_be = 4)
    starting <- list(rho = rep(0.3, mtdorder), alp = 7, be = 1)
  } else {
    param <- list(rho = rho, mu = 7, beta = 1, epsilon = 0.3, Pi = 0.25)
    obs <- suppressWarnings(rMTD(mtdorder, weight, "zigamma", param, size = 80)$y)
    prior <- list(alpha = 1, u_alp = 2, v_alp = 1, u_be = 2, v_be = 1)
    tuning <- list(step_size = 0.15, upper_alp = 12, upper_be = 4)
    starting <- list(rho = rep(0.3, mtdorder), alp = 7, be = 1,
                     eps = param$epsilon, Pi = param$Pi)
  }

  mcmc_param <- list(niter = 12, nburn = 6, nthin = 2)
  fit <- tsMTD(obs, mtdorder, family, "sb", prior, tuning, starting,
               mcmc_param, cop_family = 1)

  list(fit = fit, lags = as.matrix(lagmat(obs, mtdorder)))
}

test_that("predMTD predicts Gamma copula quantiles", {
  s <- prediction_copula_fit("gamma")
  pred <- predMTD(s$fit, "gamma", s$lags, probs = c(0.1, 0.9))

  expect_equal(dim(pred), c(2, nrow(s$lags)))
  expect_true(all(is.finite(pred)))
})

test_that("predMTD predicts ZI Gamma copula quantiles", {
  s <- prediction_copula_fit("zigamma")
  pred <- predMTD(s$fit, "zigamma", s$lags, probs = c(0.1, 0.9))

  expect_equal(dim(pred), c(2, nrow(s$lags)))
  expect_true(all(is.finite(pred)))
})

test_that("predMTD_2 returns Gamma copula posterior predictive means", {
  s <- prediction_copula_fit("gamma")
  pred <- predMTD_2(s$fit, "gamma", s$lags, probs = c(0.1, 0.9))

  expect_length(pred, nrow(s$lags))
  expect_true(all(is.finite(pred)))
})

test_that("predMTD_2 returns ZI Gamma copula posterior predictive means", {
  s <- prediction_copula_fit("zigamma")
  pred <- predMTD_2(s$fit, "zigamma", s$lags, probs = c(0.1, 0.9))

  expect_length(pred, nrow(s$lags))
  expect_true(all(is.finite(pred)))
})
