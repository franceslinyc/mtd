##############################################################################
### Gibbs sampler for Gaussian MTD
##############################################################################
gibbsGMTD <- function(obs, mtdorder, weight_param, prior, tuning, starting, updateWeight, niter, nburn, nthin){
  
  # prior
  mu_0 <- prior$mu_0
  sigma2_0 <- prior$sigma2_0
  u_0 <- prior$u_0
  v_0 <- prior$v_0
  
  # tuning parameters
  step_size <- tuning$step_size
  
  #################
  # initialization
  #################
  ## initial values of the parameters
  rho <- starting$rho
  sigma2 <- starting$sigma2
  lag_mat <- as.matrix(lagmat(obs, mtdorder))
  dat <- obs[-(1:mtdorder)]
  trun_size <- length(obs) - mtdorder
  weight <- extraDistr::rdirichlet(1, rep(1 / mtdorder, mtdorder))
  data_label <- sample(1:mtdorder, trun_size, replace = TRUE,  weight)
  unique_label <- as.numeric(names(table(data_label)))
  n_j <- rep(0, mtdorder)
  n_j[unique_label] <- as.numeric(table(data_label))
  dep_lag <- as.numeric(array(NA, dim = trun_size))
  for (t in 1:trun_size){
    dep_lag[t] <- lag_mat[t, data_label[t]]
  }
  dat_rho <- rho[data_label]
  logprho <- function(x, dd, lags, mu, sigma2) {
    logp <- dnorm(dd, (1 - x) * mu + x * lags, sqrt(sigma2 * (1 - x^2)), log = TRUE)
    sum(logp)
  }
  
  ## empty stuffs to save samples
  weight_save <- array(NA, dim = c(mtdorder, niter - nburn))
  rho_save <- array(NA, dim = c(mtdorder, niter - nburn))
  mu_save <- array(NA, dim = niter - nburn)
  sigma2_save <- array(NA, dim = niter - nburn)
  
  #######
  # mcmc
  #######
  for (iter in 1:niter){
    
    ## update mu
    mu <- updateGMTDmu(dat, dep_lag, dat_rho, sigma2, mu_0, sigma2_0)

    ## update sigma2
    sigma2 <- updateGMTDsigma2(dat, dep_lag, dat_rho, mu, trun_size, u_0, v_0)
    
    ## update rho
    rho <- updateGMTDrho(rho, dat, dep_lag, data_label, mtdorder, step_size, mu, sigma2, logprho)
    
    ## update labels
    labels <- updateGMTDlabel(dat, lag_mat, mtdorder, weight, rep(mu, mtdorder), rep(sigma2, mtdorder), rho)
    data_label <- as.numeric(labels$data_label)
    dep_lag <- as.numeric(labels$dep_lag)
    unique_label <- as.numeric(names(table(data_label)))
    n_j <- rep(0, mtdorder)
    n_j[unique_label] <- as.numeric(table(data_label))
    dat_rho <- rho[data_label]
    
    ## sample weights
    weight <- updateWeight(n_j, weight_param)
    
    ## save samples
    if (iter > nburn){
      weight_save[, iter - nburn] <- weight
      rho_save[, iter - nburn] <- rho
      mu_save[iter - nburn] <- mu
      sigma2_save[iter - nburn] <- sigma2
    }
  }
  
  ######################
  # thinning and output
  ######################
  selc_index <- seq(1, niter - nburn, by = nthin)
  list(weight = weight_save[, selc_index], 
       rho = rho_save[, selc_index],
       mu = mu_save[selc_index],
       sigma2 = sigma2_save[selc_index])
  
}

##############################################################################
### Gibbs sampler for Copula Gaussian MTD
##############################################################################
gibbsCopGMTD <- function(obs, mtdorder, weight_param, prior, tuning, starting, updateWeight, niter, nburn, nthin, nreport, 
                         cop_family, marg_family, verbose){ 
  
  step_size <- tuning$step_size
  rho <- starting$rho         

  mu_0 <- prior$mu_0
  sigma2_0 <- prior$sigma2_0
  u_0 <- prior$u_0
  v_0 <- prior$v_0

  upper_mu <- tuning$upper_mu
  lower_mu <- tuning$lower_mu
  upper_sigma2 <- tuning$upper_sigma2
  lower_sigma2 <- tuning$lower_sigma2

  mu <- starting$mu
  sigma2 <- starting$sigma2

  if (verbose) {
    
    if (is.null(nreport)) nreport <- 1000
    
    cat("--------------------------------------------\n")
    cat("\t  Running MCMC\n");
    cat("--------------------------------------------\n")
    
  }

  lag_mat <- as.matrix(lagmat(obs, mtdorder))
  dat <- obs[-(1:mtdorder)]
  trun_size <- length(obs) - mtdorder
  weight <- extraDistr::rdirichlet(1, rep(1 / mtdorder, mtdorder))

  data_label <- sample(1:mtdorder, trun_size, replace = TRUE,  weight)
  
  unique_label <- as.numeric(names(table(data_label)))
  n_j <- rep(0, mtdorder)
  n_j[unique_label] <- as.numeric(table(data_label))

  dep_lag <- as.numeric(array(NA, dim = trun_size))
  for (t in 1:trun_size){
    dep_lag[t] <- lag_mat[t, data_label[t]]
  }
  dat_rho <- rho[data_label]

  weight_save <- array(NA, dim = c(mtdorder, niter - nburn))
  rho_save <- array(NA, dim = c(mtdorder, niter - nburn))
  mu_save <- array(NA, dim = niter - nburn)
  sigma2_save <- array(NA, dim = niter - nburn)
  data_label_save <- array(NA, dim = c(trun_size, niter - nburn))
  dep_lag_save <- array(NA, dim = c(trun_size, niter - nburn))

  block_runtime <- 0

  for (iter in 1:niter){
    
    start_time <- Sys.time()

    mu <- updateCopGMTDmu(mu, dat, dep_lag, mu_0, sigma2_0, dat_rho, sigma2, trun_size, step_size, 
                          cop_family, marg_family, lower_mu, upper_mu)
    
    sigma2 <- updateCopGMTDsigma2(sigma2, dat, dep_lag, u_0, v_0, dat_rho, mu, trun_size, step_size,
                                  cop_family, marg_family, lower_sigma2, upper_sigma2)

    rho <- updateCopGMTDrho(rho, dat, dep_lag, data_label, mtdorder, step_size, mu, sigma2, 
                            cop_family, marg_family)
    
    yy_param1 <- rep(mu, trun_size)
    yy_param2 <- rep(sqrt(sigma2), trun_size)                           
    yy_ne_param1 <- matrix(rep(yy_param1, mtdorder), ncol = mtdorder)
    yy_ne_param2 <- matrix(rep(yy_param2, mtdorder), ncol = mtdorder)    
    labels <- updateCopGMTDlabel(dat, lag_mat, mtdorder, weight, rho, 
                                 cop_family, marg_family,
                                 yy_param1, yy_param2, yy_ne_param1, yy_ne_param2)
    data_label <- as.numeric(labels$data_label)
    dep_lag <- as.numeric(labels$dep_lag)
    unique_label <- as.numeric(names(table(data_label)))
    n_j <- rep(0, mtdorder)
    n_j[unique_label] <- as.numeric(table(data_label))
    dat_rho <- rho[data_label]
    
    weight <- updateWeight(n_j, weight_param)
    
    end_time <- Sys.time()
    block_runtime <- block_runtime + as.numeric(end_time - start_time)
    
    if (verbose) {
      if (iter %% nreport == 0) {
        cat(paste0("Iterations: ", iter, "/", niter, 
                   "  Percentage: ", specRound(iter / niter * 100, 2), "%\n"))
        ert <- (block_runtime / nreport) * (niter - iter)
        cat(paste0("Estimated remaining time: ", 
                   specRound(ert / 60, 2), " minutes \n"))
        block_runtime <- 0
        cat("--------------------------------------------\n")
      }
    }  

    if (iter > nburn){
      weight_save[, iter - nburn] <- weight
      rho_save[, iter - nburn] <- rho
      mu_save[iter - nburn] <- mu
      sigma2_save[iter - nburn] <- sigma2
      data_label_save[, iter - nburn] <- data_label
      dep_lag_save[, iter - nburn] <- dep_lag
    }
  }
  
  selc_index <- seq(1, niter - nburn, by = nthin)
  list(weight = weight_save[, selc_index], 
       rho = rho_save[, selc_index],
       shape1 = mu_save[selc_index],
       shape2 = sigma2_save[selc_index], 
       data_label = data_label_save[, selc_index],
       dep_lag = dep_lag_save[, selc_index]
       )
  
}

##############################################################################
### Gibbs sampler for Copula Gamma MTD
##############################################################################
gibbsCopGammaMTD <- function(obs, mtdorder, weight_param, prior, tuning, starting, updateWeight, niter, nburn, nthin, nreport, 
                             cop_family, marg_family, verbose){

  # prior
  u_alp <- prior$u_alp
  v_alp <- prior$v_alp
  u_be <- prior$u_be
  v_be <- prior$v_be

  # tuning parameters
  step_size <- tuning$step_size
  upper_alp <- tuning$upper_alp
  upper_be <- tuning$upper_be

  # starting values
  alp <- starting$alp
  be <- starting$be
  rho <- starting$rho        

  #--------------------------
  # MCMC settings
  #--------------------------

  if (verbose) {
    
    if (is.null(nreport)) nreport <- 1000
    
    cat("--------------------------------------------\n")
    cat("\t  Running MCMC\n");
    cat("--------------------------------------------\n")
    
  }

  #################
  # initialization
  #################  

  lag_mat <- as.matrix(lagmat(obs, mtdorder))
  dat <- obs[-(1:mtdorder)]
  trun_size <- length(obs) - mtdorder
  weight <- extraDistr::rdirichlet(1, rep(1 / mtdorder, mtdorder))

  data_label <- sample(1:mtdorder, trun_size, replace = TRUE, weight)
  
  unique_label <- as.numeric(names(table(data_label)))
  n_j <- rep(0, mtdorder)
  n_j[unique_label] <- as.numeric(table(data_label))

  dep_lag <- as.numeric(array(NA, dim = trun_size))
  for (t in 1:trun_size){
    dep_lag[t] <- lag_mat[t, data_label[t]]
  }
  dat_rho <- rho[data_label]

  ## empty stuffs to save samples
  weight_save <- array(NA, dim = c(mtdorder, niter - nburn))
  rho_save <- array(NA, dim = c(mtdorder, niter - nburn))
  alp_save <- array(NA, dim = niter - nburn)
  be_save <- array(NA, dim = niter - nburn)

  #######
  # mcmc
  #######

  block_runtime <- 0

  for (iter in 1:niter){
    
    start_time <- Sys.time()
  
    ## update alp
    alp <- updateCopGammaMTDalp(alp, dat, dep_lag, u_alp, v_alp, dat_rho, be, trun_size, step_size, 
                                cop_family, marg_family, upper_alp)

    ## update be (rate parameter)
    be <- updateCopGammaMTDbe(be, dat, dep_lag, u_be, v_be, dat_rho, alp, trun_size, step_size, 
                              cop_family, marg_family, upper_be)

    ## update rho  
    rho <- updateCopGammaMTDrho(rho, dat, dep_lag, data_label, mtdorder, step_size, alp, be,
                                cop_family, marg_family)

    ## update labels    
    yy_param1 <- rep(alp, trun_size)
    yy_param2 <- rep(be, trun_size)                           
    yy_ne_param1 <- matrix(rep(yy_param1, mtdorder), ncol = mtdorder) 
    yy_ne_param2 <- matrix(rep(yy_param2, mtdorder), ncol = mtdorder)    
    labels <- updateCopGMTDlabel(dat, lag_mat, mtdorder, weight, rho,
                                 cop_family, marg_family,
                                 yy_param1, yy_param2, yy_ne_param1, yy_ne_param2)
    data_label <- as.numeric(labels$data_label)
    dep_lag <- as.numeric(labels$dep_lag)
    unique_label <- as.numeric(names(table(data_label)))
    n_j <- rep(0, mtdorder)
    n_j[unique_label] <- as.numeric(table(data_label))
    dat_rho <- rho[data_label]

    ## sample weights    
    weight <- updateWeight(n_j, weight_param)

    #--------------------------
    # Print MCMC progress
    #--------------------------

    end_time <- Sys.time()
    block_runtime <- block_runtime + as.numeric(end_time - start_time)
    
    if (verbose) {
      if (iter %% nreport == 0) {
        cat(paste0("Iterations: ", iter, "/", niter, 
                   "  Percentage: ", specRound(iter / niter * 100, 2), "%\n"))
        ert <- (block_runtime / nreport) * (niter - iter)
        cat(paste0("Estimated remaining time: ", 
                   specRound(ert / 60, 2), " minutes \n"))
        block_runtime <- 0
        cat("--------------------------------------------\n")
      }
    }  

    ## save samples
    if (iter > nburn){
      weight_save[, iter - nburn] <- weight
      rho_save[, iter - nburn] <- rho
      alp_save[iter - nburn] <- alp
      be_save[iter - nburn] <- be
    }
  }

  ######################
  # thinning and output
  ######################

  selc_index <- seq(1, niter - nburn, by = nthin)
  list(weight = weight_save[, selc_index], 
       rho = rho_save[, selc_index],
       shape1 = alp_save[selc_index],
       shape2 = be_save[selc_index]
       )
  
}  

##############################################################################
### Gibbs sampler for Copula ZI Gamma MTD
##############################################################################
gibbsCopZIGammaMTD <- function(obs, mtdorder, weight_param, prior, tuning, starting, updateWeight, niter, nburn, nthin, nreport, 
                               cop_family, marg_family, verbose){

  # prior
  u_alp <- prior$u_alp
  v_alp <- prior$v_alp
  u_be <- prior$u_be
  v_be <- prior$v_be

  # tuning parameters
  step_size <- tuning$step_size
  upper_alp <- tuning$upper_alp
  upper_be <- tuning$upper_be

  # starting values
  alp <- starting$alp
  be <- starting$be 
  rho <- starting$rho       
  eps <- starting$eps
  Pi <- starting$Pi

  #--------------------------
  # MCMC settings
  #--------------------------

  if (verbose) {
    
    if (is.null(nreport)) nreport <- 1000
    
    cat("--------------------------------------------\n")
    cat("\t  Running MCMC\n");
    cat("--------------------------------------------\n")
    
  }

  #################
  # initialization
  #################  

  lag_mat <- as.matrix(lagmat(obs, mtdorder))
  dat <- obs[-(1:mtdorder)]
  trun_size <- length(obs) - mtdorder
  weight <- extraDistr::rdirichlet(1, rep(1 / mtdorder, mtdorder))

  data_label <- sample(1:mtdorder, trun_size, replace = TRUE, weight)
  
  unique_label <- as.numeric(names(table(data_label)))
  n_j <- rep(0, mtdorder)
  n_j[unique_label] <- as.numeric(table(data_label))

  dep_lag <- as.numeric(array(NA, dim = trun_size))
  for (t in 1:trun_size){
    dep_lag[t] <- lag_mat[t, data_label[t]]
  }
  dat_rho <- rho[data_label]

  ## empty stuffs to save samples
  weight_save <- array(NA, dim = c(mtdorder, niter - nburn))
  rho_save <- array(NA, dim = c(mtdorder, niter - nburn))
  alp_save <- array(NA, dim = niter - nburn)
  be_save <- array(NA, dim = niter - nburn)
  Pi_save <- array(NA, dim = niter - nburn)
  #eps_save <- array(NA, dim = niter - nburn)

  #######
  # mcmc
  #######

  block_runtime <- 0

  for (iter in 1:niter){
    
    start_time <- Sys.time()

    ## update mu
    alp <- updateZICopGammaMTDalp(alp, dat, dep_lag, u_alp, v_alp, dat_rho, be, eps, Pi, trun_size, step_size, 
                                  cop_family, marg_family, upper_alp)

    ## update be (scale parameter)
    be <- updateZICopGammaMTDbe(be, dat, dep_lag, u_be, v_be, dat_rho, alp, eps, Pi, trun_size, step_size, 
                                cop_family, marg_family, upper_be)

    ## update eps
    # eps <- updateZICopGammaMTDeps(eps, dat, dep_lag, eps_0, dat_rho, alp, be, Pi, trun_size, step_size, 
    #                               cop_family, marg_family, upper_eps)

    ## update pi 
    Pi <- updateZICopGammaMTDPi(Pi, dat, dep_lag, dat_rho, alp, be, eps, trun_size, step_size, 
                                cop_family, marg_family)    

    ## update rho                                
    rho <- updateCopZIGammaMTDrho(rho, dat, dep_lag, data_label, mtdorder, step_size, alp, be, eps, Pi, 
                                  cop_family, marg_family)
    
    ## update labels    
    yy_param1 <- rep(alp, trun_size)
    yy_param2 <- rep(be, trun_size) 
    yy_param3 <- rep(eps, trun_size)
    yy_param4 <- rep(Pi, trun_size)
    yy_ne_param1 <- matrix(rep(yy_param1, mtdorder), ncol = mtdorder) 
    yy_ne_param2 <- matrix(rep(yy_param2, mtdorder), ncol = mtdorder) 
    yy_ne_param3 <- matrix(rep(yy_param3, mtdorder), ncol = mtdorder)
    yy_ne_param4 <- matrix(rep(yy_param4, mtdorder), ncol = mtdorder)
    labels <- updateZIGMTDlabel(dat, lag_mat, mtdorder, weight, rho,
                                cop_family, marg_family,
                                yy_param1, yy_param2, yy_ne_param1, yy_ne_param2, 
                                yy_param3, yy_param4, yy_ne_param3, yy_ne_param4)
    data_label <- as.numeric(labels$data_label)
    dep_lag <- as.numeric(labels$dep_lag)
    unique_label <- as.numeric(names(table(data_label)))
    n_j <- rep(0, mtdorder)
    n_j[unique_label] <- as.numeric(table(data_label))
    dat_rho <- rho[data_label]

    ## sample weights
    weight <- updateWeight(n_j, weight_param)

    #--------------------------
    # Print MCMC progress
    #--------------------------

    end_time <- Sys.time()
    block_runtime <- block_runtime + as.numeric(end_time - start_time)
    
    if (verbose) {
      if (iter %% nreport == 0) {
        cat(paste0("Iterations: ", iter, "/", niter, 
                   "  Percentage: ", specRound(iter / niter * 100, 2), "%\n"))
        ert <- (block_runtime / nreport) * (niter - iter)
        cat(paste0("Estimated remaining time: ", 
                   specRound(ert / 60, 2), " minutes \n"))
        block_runtime <- 0
        cat("--------------------------------------------\n")
      }
    }  

    ## save samples
    if (iter > nburn){
      weight_save[, iter - nburn] <- weight
      rho_save[, iter - nburn] <- rho
      alp_save[iter - nburn] <- alp
      be_save[iter - nburn] <- be
      Pi_save[iter - nburn] <- Pi
      #eps_save[iter - nburn] <- eps
    }
  }

  ######################
  # thinning and output
  ######################

  selc_index <- seq(1, niter - nburn, by = nthin)
  list(weight = weight_save[, selc_index], 
       rho = rho_save[, selc_index],
       shape1 = alp_save[selc_index],
       shape2 = be_save[selc_index],
       Pi = Pi_save[selc_index]#,
       #eps = eps_save[selc_index]
       )
  
} 

##############################################################################
### Gibbs sampler for Poisson MTD
##############################################################################
gibbsPMTD <- function(obs, mtdorder, weight_param, prior, tuning, starting, updateWeight, niter, nburn, nthin){
  
  # Prior
  u_la <- prior$u_la
  v_la<- prior$v_la
  u_th <- prior$u_th
  v_th <- prior$v_th
  
  #################
  # initialization
  #################
  ## initial values of the parameters
  lag_mat <- as.matrix(lagmat(obs, mtdorder))
  dat <- obs[-(1:mtdorder)]
  trun_size <- length(obs) - mtdorder
  weight <- extraDistr::rdirichlet(1, rep(1 / mtdorder, mtdorder))
  data_label <- sample(1:mtdorder, trun_size, replace = TRUE,  weight)
  unique_label <- as.numeric(names(table(data_label)))
  n_j <- as.numeric(table(data_label))
  dep_lag <- array(NA, dim = trun_size)
  for (t in 1:trun_size) dep_lag[t] <- lag_mat[t, data_label[t]]
  temp_ut <- rpois(trun_size, 1)
  temp_ut2 <- apply(cbind(temp_ut, dat - dep_lag), 1, max)
  ut <- apply(cbind(temp_ut2, dat), 1, min)
  
  ## empty stuffs to save samples
  weight_save <- array(NA, dim = c(mtdorder, niter - nburn))
  ut_save <- array(NA, dim = c(trun_size, niter - nburn))
  la_save <- array(NA, dim = niter - nburn)
  th_save <- array(NA, dim = niter - nburn)
  
  #######
  # mcmc
  #######
  for (iter in 1:niter){
    
    ## update th
    th <- updatePMTDth(dat, dep_lag, ut, u_th, v_th)
 
    ## update lambda
    la <- updatePMTDla(ut, trun_size, u_la, v_la)

    ## update ut
    ut <- as.numeric(updatePMTDut(ut, dat, dep_lag, la, la/(1-th)-la))
    
    ## update labels
    labels <- updatePMTDlabel(dat, lag_mat, mtdorder, weight, th, ut)
    data_label <- as.numeric(labels$data_label)
    dep_lag <- as.numeric(labels$dep_lag)
    unique_label <- as.numeric(names(table(data_label)))
    n_j <- rep(0, mtdorder)
    n_j[unique_label] <- as.numeric(table(data_label))
    
    ## update weights
    weight <- updateWeight(n_j, weight_param)
    
    ## save samples
    if (iter > nburn){
      weight_save[, iter - nburn] <- weight
      ut_save[, iter - nburn] <- ut
      la_save[iter - nburn] <- la
      th_save[iter - nburn] <- th
    }
  }
  
  ####################
  # thinning and output
  ######################  
  selc_index <- seq(1, niter - nburn, by = nthin)
  list(weight = weight_save[, selc_index],
       ut = ut_save[, selc_index],
       la = la_save[selc_index],
       th = th_save[selc_index])
}

##############################################################################
### Gibbs sampler for negative binomial MTD
#############################################################################
gibbsNBMTD <- function(obs, mtdorder, weight_param, prior, tuning, starting, updateWeight, niter, nburn, nthin) {
  
  # prior
  u_th <- prior$u_th
  v_th <- prior$v_th
  u_psi <- prior$u_psi
  v_psi <- prior$v_psi
  u_kap <- prior$u_kap
  v_kap <- prior$v_kap
  
  # tuning parameters
  se_kap <- tuning$se_kap
  
  #################
  # initialization
  #################
  ## initial values of the parameters
  kap <- starting$kap
  lag_mat <- as.matrix(lagmat(obs, mtdorder))
  dat <- obs[-(1:mtdorder)]
  trun_size <- length(obs) - mtdorder
  weight <- extraDistr::rdirichlet(1, rep(1 / mtdorder, mtdorder))
  data_label <- sample(1:mtdorder, trun_size, replace = TRUE,  weight)
  unique_label <- as.numeric(names(table(data_label)))
  n_j <- as.numeric(table(data_label))
  dep_lag <- array(NA, dim = trun_size)
  for (t in 1:trun_size) dep_lag[t] <- lag_mat[t, data_label[t]]
  temp_ut <- rpois(trun_size, 1)
  temp_ut2 <- apply(cbind(temp_ut, dat - dep_lag), 1, max)
  ut <- apply(cbind(temp_ut2, dat), 1, min)
  
  ## empty stuffs to save samples
  weight_save <- array(NA, dim = c(mtdorder, niter - nburn))
  ut_save <- array(NA, dim = c(trun_size, niter - nburn))
  psi_save <- array(NA, dim = niter - nburn)
  th_save <- array(NA, dim = niter - nburn)
  kap_save <- array(NA, dim = niter - nburn)
  kap_accept_count <- array(0, dim = niter)
  
  #######
  # mcmc
  #######
  for (iter in 1:niter){
    
    ## update th
    th <- updateNBMTDth(dat, dep_lag, ut, u_th, v_th)

    ## update psi
    psi <- updateNBMTDpsi(dep_lag, ut, trun_size, kap, u_psi, v_psi)

    ## update ut
    ut <- as.numeric(updateNBMTDut(ut, dat, dep_lag, th, psi, kap))
    
    ## update kap
    kap_res <- updateNBMTDkap(kap, se_kap, ut, dep_lag, psi, u_kap, v_kap)
    kap <- kap_res$kap
    if (kap_res$accept) kap_accept_count[iter] <- 1
    
    ## update labels
    labels <- updateNBMTDlabel(dat, lag_mat, mtdorder, weight, th, ut, psi, kap)
    data_label <- as.numeric(labels$data_label)
    dep_lag <- as.numeric(labels$dep_lag)
    unique_label <- as.numeric(names(table(data_label)))
    n_j <- rep(0, mtdorder)
    n_j[unique_label] <- as.numeric(table(data_label))
    
    ## update weights
    weight <- updateWeight(n_j, weight_param)
    
    ## save samples
    if (iter > nburn){
      weight_save[, iter - nburn] <- weight
      ut_save[, iter - nburn] <- ut
      psi_save[iter - nburn] <- psi
      th_save[iter - nburn] <- th
      kap_save[iter - nburn] <- kap
    }
  }
  
  ####################
  # thinning and output
  ######################  
  selc_index <- seq(1, niter - nburn, by = nthin)
  list(weight = weight_save[, selc_index],
       ut = ut_save[, selc_index],
       psi = psi_save[selc_index],
       th = th_save[selc_index],
       kap = kap_save[selc_index],
       accept_cout = kap_accept_count)
}

##############################################################################
### Gibbs sampler for Lomax MTD regression model 
#############################################################################
gibbsRegLMTD <- function(obs, XX, mtdorder, weight_param, prior, tuning, starting, updateWeight, niter, nburn, nthin){
  
  # prior
  u_alpha <- prior$u_alpha
  v_alpha <- prior$v_alpha
  u_phi <- prior$u_phi
  v_phi <- prior$v_phi
  
  # tuning parameter
  se_phi <- tuning$se_phi
  se_bb <- tuning$se_bb
  
  #################
  # initialization
  ################
  ## initial values of the parameters
  phi <- starting$phi
  bb <- starting$bb
  mu_obs <- exp(XX %*% bb)
  eps_lag_mat <- as.matrix(lagmat(obs / mu_obs, mtdorder))
  eps_dat <- (obs / mu_obs)[-(1:mtdorder)]
  trun_size <- length(obs) - mtdorder
  weight <- extraDistr::rdirichlet(1, rep(1 / mtdorder, mtdorder))
  data_label <- sample(1:mtdorder, trun_size, replace = TRUE,  weight)
  unique_label <- as.numeric(names(table(data_label)))
  n_j <- as.numeric(table(data_label))
  eps_dep_lag <- array(NA, dim = trun_size)
  for (t in 1:trun_size){
    eps_dep_lag[t] <- eps_lag_mat[t, data_label[t]]
  }
  
  ## empty stuffs to save samples
  weight_save <- array(NA, dim = c(mtdorder, niter - nburn))
  alpha_save <- array(NA, dim = niter - nburn)
  phi_save <- array(NA, dim = niter - nburn)
  bb_save <- array(NA, dim = c(ncol(XX), niter - nburn))
  phi_accept_count <- array(0, dim = niter)
  bb_accept_count <- array(0, dim = c(length(bb), niter))
  
  # mcmc
  
  for (iter in 1:niter){
    
    ## update alpha
    alpha <- updateRegLMTDalpha(eps_dat, eps_dep_lag, trun_size, phi, u_alpha, v_alpha)
    
    ## update phi
    phi_res <- updateRegLMTDphi(phi, se_phi, eps_dat, eps_dep_lag, trun_size, u_alpha, v_alpha, u_phi, v_phi)
    phi <- phi_res$phi
    if (phi_res$accept) phi_accept_count[iter] <- 1
    
    ## update beta
    bb_res <- updateRegLMTDbeta(bb, se_bb, bb_accept_count, eps_dat, eps_dep_lag, eps_lag_mat,
                                mtdorder, obs, XX, data_label, trun_size, phi, u_alpha, v_alpha, iter)
    
    bb <- bb_res$bb
    bb_accept_count <- bb_res$bb_accept_count
    eps_dat <- bb_res$eps_dat
    eps_dep_lag <- bb_res$eps_dep_lag
    eps_lag_mat <- bb_res$eps_lag_mat
    
    ## update labels
    labels <- updateLMTDlabel(eps_dat, eps_lag_mat, mtdorder, weight, rep(alpha, mtdorder), rep(phi, mtdorder))
    data_label <- as.numeric(labels$data_label)
    eps_dep_lag <- as.numeric(labels$dep_lag)
    unique_label <- as.numeric(names(table(data_label)))
    n_j <- rep(0, mtdorder)
    n_j[unique_label] <- as.numeric(table(data_label))
    
    ## update weights
    weight <- updateWeight(n_j, weight_param)
    
    ## save samples
    if (iter > nburn){
      weight_save[, iter - nburn] <- weight
      alpha_save[iter - nburn] <- alpha
      phi_save[iter - nburn] <- phi
      bb_save[, iter - nburn] <- bb
    }
  }
  
  ######################
  # thinning and output
  ######################
  selc_index <- seq(1, niter - nburn, by = nthin)
  list(weight = weight_save[, selc_index],
       alpha = alpha_save[selc_index],
       phi = phi_save[selc_index],
       phi_accept_count = phi_accept_count,
       bb = bb_save[, selc_index])
}