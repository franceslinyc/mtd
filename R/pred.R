#' One-step ahead prediction
#'
#' This function generates one-step-ahead posterior predictive distribution 
#' given a fitted MTD model.
#
#' @param res a list of MCMC results from the MTD model; an output from functions 
#'        either \code{\link{tsMTD}} or \code{\link{tsRegMTD}}.
#' @param family a quoted keyword that specifies the stationary marginal distribution.
#'        Supported keywords are: "poisson", "negative binomial", "gamma", "zigamma".
#' @param lags a numeric matrix of the lags of the observation. Each row 
#'        corresponds to an observation. The dimension of the matrix is 
#'        truncated data size x mtdorder.
#' @param probs a numeric vector of probabilities with values in \eqn{[0, 1]}.
#'
#' @references 
#' Zheng, X., Kottas, A., and Sansó, B. (2021),
#' "On construction and estimation of stationary mixture transition distribution models,"
#' \emph{arXiv:2010.12696}.
#' 
#' @export
#'
predMTD <- function(res, family, lags, probs){
  
  if (missing(family)) {
    
    stop("error: family must be specified")
    
  }
  
  family <- tolower(family)
  family_names <- c("poisson", "negative binomial", "gamma", "zigamma")
  
  if (!family %in% family_names) {
    
    stop("error: specified family '", family, "' is not a valid option;
         available families are ", 
         paste(family_names, collapse = ", ", sep = "") ,".")
    
  }
  
  if (any(probs <= 0, probs >= 1)) {
    
    stop("error: probs must be a positive vector and each element must be between 0 and 1.")
    
  }
  
  if (family == "poisson") {
    
    weight <- res$weight
    la <- res$la
    th <- res$th
    
    if (nrow(weight) != ncol(lags)) {
      
      stop("error: the model order used to fit the MTD is different from the number of lags of each observaton.")
      
    }
    
    if (is.null(la) || is.null(th) || is.null(weight)) {
      
      stop("error: lack of posterior samples.")
      
    }
    
    mtdorder <- nrow(weight)
    
    out <- predPMTD(mtdorder, weight, lags, la, th, probs) 
      
  }
  else if (family == "negative binomial") {
    
    weight <- res$weight
    th <- res$th
    psi <- res$psi
    kap <- res$kap
    
    if (nrow(weight) != ncol(lags)) {
      
      stop("error: the model order used to fit the MTD is different from the number of lags of each observaton.")
      
    }
    
    
    if (is.null(th) || is.null(psi) || is.null(kap) || is.null(weight)) {
      
      stop("error: lack of posterior samples.")
      
    }
    
    mtdorder <- nrow(weight)
    
    out <- predNBMTD(mtdorder, weight, lags, th, psi, kap, probs) 
    
  }
  else if (family == "gamma") {
    weight <- res$weight
    rho <- res$rho
    alp <- res$shape1
    be <- res$shape2

    if (nrow(weight) != ncol(lags) || nrow(rho) != ncol(lags)) {
      stop("error: the model order used to fit the MTD is different from the number of lags of each observation.")
    }

    if (is.null(rho) || is.null(alp) || is.null(be) || is.null(weight)) {
      stop("error: lack of posterior samples.")
    }

    mtdorder <- nrow(weight)
    cop_family <- 1
    marg_family <- 2
    out <- predCopGammaMTD(mtdorder, weight, lags, alp, be, cop_family, rho, marg_family, probs)
  }
  else if (family == "zigamma") {
    weight <- res$weight
    rho <- res$rho
    alp <- res$shape1
    be <- res$shape2
    Pi <- res$Pi
    eps <- res$eps

    if (nrow(weight) != ncol(lags) || nrow(rho) != ncol(lags)) {
      stop("error: the model order used to fit the MTD is different from the number of lags of each observation.")
    }

    if (is.null(rho) || is.null(alp) || is.null(be) || is.null(weight) || is.null(Pi) || is.null(eps)) {
      stop("error: lack of posterior samples.")
    }

    mtdorder <- nrow(weight)
    cop_family <- 1
    marg_family <- 4
    out <- predCopZIGammaMTD(mtdorder, weight, lags, alp, be, Pi, eps, cop_family, rho, marg_family, probs)
  }
  
  out
  
}

predMTD_2 <- function(res, family, lags, probs){
  if (family == "gamma") {
    weight <- res$weight
    rho <- res$rho
    alp <- res$shape1
    be <- res$shape2

    if (nrow(weight) != ncol(lags) || nrow(rho) != ncol(lags)) {
      stop("error: the model order used to fit the MTD is different from the number of lags of each observation.")
    }

    if (is.null(rho) || is.null(alp) || is.null(be) || is.null(weight)) {
      stop("error: lack of posterior samples.")
    }

    mtdorder <- nrow(weight)
    cop_family <- 1
    marg_family <- 2
    out <- predCopGammaMTD_2(mtdorder, weight, lags, alp, be, cop_family, rho, marg_family, probs)
  }
  else if (family == "zigamma") {
    weight <- res$weight
    rho <- res$rho
    alp <- res$shape1
    be <- res$shape2
    Pi <- res$Pi
    eps <- res$eps

    if (nrow(weight) != ncol(lags) || nrow(rho) != ncol(lags)) {
      stop("error: the model order used to fit the MTD is different from the number of lags of each observation.")
    }

    if (is.null(rho) || is.null(alp) || is.null(be) || is.null(weight) || is.null(Pi) || is.null(eps)) {
      stop("error: lack of posterior samples.")
    }

    mtdorder <- nrow(weight)
    cop_family <- 1
    marg_family <- 4
    out <- predCopZIGammaMTD_2(mtdorder, weight, lags, alp, be, Pi, eps, cop_family, rho, marg_family, probs)
  }

  out
}