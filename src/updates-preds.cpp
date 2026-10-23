#include <iostream>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
#include "utils.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

#define PMAX  1.0 - 1e-10
#define PMIN  1e-10

// Adapted from nnmp (X. Zheng), copulas.cpp
double dBiCop_cpp(const double& uu,
                  const double& vv,
                  const int& cop_family,
                  const double& cop_param,
                  const bool& logd = false) {
  
  double logdens;
  
  if (cop_family == 1) {
    double xx = R::qnorm(uu, 0.0, 1.0, true, false);
    double yy = R::qnorm(vv, 0.0, 1.0, true, false);
    double rho = cop_param;
    double rho2 = pow(rho, 2);
    double xx2 = pow(xx, 2);
    double yy2 = pow(yy, 2);
    double ee = rho2 * (xx2 + yy2) - 2.0 * rho * xx * yy;
    logdens = -.5 * ee / (1.0 - rho2) - .5 * log(1.0 - rho2);
  } 
  else if (cop_family == 2) {
    double eta = cop_param;
    if (eta < 1) {
      stop("error: gumbel copula parameter must be greater than 1.");
    }
    if (eta > 50.0) {
      eta = 50.0;
    }
    double xx = -log(uu);
    double yy = -log(vv);
    double tt = pow(xx, eta) + pow(yy, eta);
    logdens = -pow(tt, 1.0 / eta) + log(pow(tt, 1.0 / eta) + eta - 1.0) + 
      (1.0 / eta - 2.0) * log(tt) + (eta - 1.0) * (log(xx) + log(yy)) + xx + yy;
  }
  else if (cop_family == 3) {
    double delta = cop_param;
    if (delta > 98) {
      delta = 98;
    }
    double ee = pow(uu, -delta) + pow(vv, -delta) - 1.0;
    logdens = log(1.0 + delta) - (delta + 1.0) * (log(uu) + log(vv)) - (2.0 + 1.0 / delta) * log(ee);
  }
  else {
    stop("error: no families are found.");
  }
  
  if (logd) {
    return logdens;
  } else {
    return exp(logdens);
  }
  
}

// Adapted from nnmp (X. Zheng), copulas.cpp
arma::colvec dBiCop_cpp2(const arma::colvec& uu,
                         const arma::colvec& vv,
                         const int& cop_family,
                         const arma::colvec& cop_param,
                         const bool& logd = false) {
  
  int nn = uu.n_rows;
  arma::colvec dens(nn);
  
  for (int i = 0; i < nn; ++i) {
    dens(i) = dBiCop_cpp(uu(i), vv(i), cop_family, cop_param(i), logd);
  } 
  
  return dens;
  
}

// [[Rcpp::export]]
double psgamma(double x, double mu, double beta, double eps){
  return (x >= eps) * R::pgamma(x - eps, mu / beta, beta, true, false);
}

// [[Rcpp::export]]
double pzig(double x, double mu, double beta, double eps, double pi){
  return pi * R::punif(x, 0, eps, true, false) + (1 - pi) * psgamma(x, mu, beta, eps); 
}

// [[Rcpp::export]]
double qsgamma(double x, double mu, double beta, double eps) {
  return R::qgamma(x, mu / beta, beta, true, false) + eps;
}

// [[Rcpp::export]]
double qzig(double u, double mu, double beta, double pi, double eps) {
  if (u <= pi) {
    return R::qunif(u / pi, 0.0, eps, true, false);
  } else {
    double adj_u = (u - pi) / (1.0 - pi);
    return qsgamma(adj_u, mu, beta, eps);
  }
}

// Adapted from nnmp (X. Zheng), copulas.cpp; marg_family 3 and 4 added
// [[Rcpp::export]]
double dBiCopMar_cpp(const double& xx,
                     const double& yy,
                     const int& cop_family,
                     const double& cop_param,
                     const int& marg_family,
                     const arma::colvec& mar_param,
                     const bool& ce = false,
                     const bool& logd = false) {
  
  double uu, vv;
  
  if (marg_family == 1) {
    uu = R::pbeta(xx, mar_param(0), mar_param(1), true, false);
    vv = R::pbeta(yy, mar_param(2), mar_param(3), true, false);
  }
  else if (marg_family == 2) {
    uu = R::pgamma(xx, mar_param(0), 1.0 / mar_param(1), true, false);
    vv = R::pgamma(yy, mar_param(2), 1.0 / mar_param(3), true, false);
  } 
  else if (marg_family == 3) {
    uu = R::pnorm(xx, mar_param(0), mar_param(1), true, false);
    vv = R::pnorm(yy, mar_param(2), mar_param(3), true, false);
  }   
  else if (marg_family == 4) {
    uu = pzig(xx, mar_param(0), mar_param(1), mar_param(4), mar_param(5));
    vv = pzig(yy, mar_param(2), mar_param(3), mar_param(6), mar_param(7));
  } 
  else {
    stop("error: no such a family is available.");
  }
  
  if (uu < PMIN) {
    uu = PMIN;
  } else if (uu > PMAX) {
    uu = PMAX;
  }
  if (vv < PMIN) {
    vv = PMIN;
  } else if (vv > PMAX) {
    vv = PMAX;
  }
  
  return dBiCop_cpp(uu, vv, cop_family, cop_param, logd);
  
}

// Adapted from nnmp (X. Zheng), copulas.cpp; marg_family 3 and 4 added
// [[Rcpp::export]]
arma::colvec dBiCopMar_cpp2(const arma::colvec& xx,
                            const arma::colvec& yy,
                            const int& cop_family,
                            const arma::colvec cop_param,
                            const int& marg_family,
                            const arma::mat& mar_param,
                            const bool& ce = false,
                            const bool& logd = false) {
  
  int nn = xx.n_rows;    
  arma::colvec uu(nn), vv(nn);
  
  if (marg_family == 1) {
    for (int i = 0; i < nn; ++i) {
      uu(i) = R::pbeta(xx(i), mar_param(i,0), mar_param(i,1), true, false);
      vv(i) = R::pbeta(yy(i), mar_param(i,2), mar_param(i,3), true, false);
    }
  }
  else if (marg_family == 2) {
    for (int i = 0; i < nn; ++i) {
      uu(i) = R::pgamma(xx(i), mar_param(i,0), 1.0 / mar_param(i,1), true, false);
      vv(i) = R::pgamma(yy(i), mar_param(i,2), 1.0 / mar_param(i,3), true, false);
    }
  }                                 
  else if (marg_family == 3) {
    for (int i = 0; i < nn; ++i) {
      uu(i) = R::pnorm(xx(i), mar_param(i,0), mar_param(i,1), true, false);
      vv(i) = R::pnorm(yy(i), mar_param(i,2), mar_param(i,3), true, false);
    }
  }
  else if (marg_family == 4) {
    for (int i = 0; i < nn; ++i) {
      uu(i) = pzig(xx(i), mar_param(i,0), mar_param(i,1), mar_param(i, 4), mar_param(i, 5));
      vv(i) = pzig(yy(i), mar_param(i,2), mar_param(i,3), mar_param(i, 6), mar_param(i, 7)); 
    }
  } 
  else {
    stop("error: no such a family is available.");
  }
  
  for (int i = 0; i < nn; ++i) {
    if (uu(i) < PMIN) {
      uu(i) = PMIN;
    } else if (uu(i) > PMAX) {
      uu(i) = PMAX;
    }
    if (vv(i) < PMIN) {
      vv(i) = PMIN;
    } else if (vv(i) > PMAX) {
      vv(i) = PMAX;
    }
  }
  
  return dBiCop_cpp2(uu, vv, cop_family, cop_param, logd);    
  
}

// Adapted from nnmp (X. Zheng), copulas.cpp
double rConGum_mnr(const double& vv,
                   const double& eta,
                   const double& tol,
                   const int& maxiter) {

  double ww = -log(vv);
  double pp = arma::randu();
  double upb = 100.0;
  double lob = ww;
  double zz0 = (lob + upb) / 2.0;

  double rr, tt, tt1, tt2, zz1, uu;

  for (int i = 0; i < maxiter; ++i) {    

    tt = zz0 + (eta - 1.0) * log(zz0) - (ww + (eta - 1.0) * log(ww) - log(pp));
    tt1 = 1.0 + (eta - 1.0) / zz0;
    tt2 = (1.0 - eta) / pow(zz0, 2);
    zz1 = zz0 - tt * tt1 / (tt1 * tt1  - tt * tt2);

    if (zz1 < lob || zz1 > upb) {
      if (tt < 0) {
        lob = zz0;
      } else {
        upb = zz0;
      }
      zz1 = (lob + upb) / 2.0;
    } 

    rr = zz1 + (eta - 1.0) * log(zz1) - (ww + (eta - 1.0) * log(ww) - log(pp));

    if (abs(rr) <= tol) {
      break;
    }

    zz0 = zz1;

    if (isnan(rr)) {
      stop("Numerical issue occurred: the function outputs nan");
    }
    if (i == (maxiter - 1)){
      stop("Failed to converge!");
    }
  }

  uu = exp(-pow((pow(zz1, eta) - pow(ww, eta)), 1.0/eta));

  if (uu == 1){
    stop("Numerical issue occurred: uu = 1");
  }

  return uu;

}

// Adapted from nnmp (X. Zheng), copulas.cpp
double rConCop(const double& vv,
               const int& cop_family,
               const double& cop_param,
               const double& tol = 1e-5,
               const int& maxiter = 1000) {

  double zz = arma::randu();
  double uu;

  if (cop_family == 1) {
    double rho = cop_param;
    double qzz = R::qnorm(zz, 0.0, 1.0, true, false);
    double qvv = R::qnorm(vv, 0.0, 1.0, true, false);
    uu = R::pnorm(sqrt(1.0 - rho * rho) * qzz + rho * qvv, 0.0, 1.0, true, false);
  }
  else if (cop_family == 2) {
    double eta = cop_param;
    if (eta > 50.0) {
      eta = 50.0;
    }
    uu = rConGum_mnr(vv, eta, tol, maxiter);
  }
  else if (cop_family == 3) {
    double delta = cop_param;
    if (delta > 98.0) {
      delta = 98.0;
    }
    double tt = exp(-delta / (1.0 + delta) * log(zz)) - 1.0;
    uu = pow(exp(log(tt) - delta * log(vv)) + 1.0, -1.0 / delta);
  }
  else {
    stop("error: no families are found.");
  }

  return uu;

}

// update auxiliary variables of the Poisson MTD model
// [[Rcpp::export]]
arma::colvec updatePMTDut(const arma::colvec& ut,
                          const arma::colvec& dat,
                          const arma::colvec& dep_lag,
                          const double& la,
                          const double& ga) {
  
  int trun_size = dat.n_rows;
  double th = ga / (la + ga);
  double ut_lb, val_len, prop_ut, cur_ut, idat, ilag_dep, prop_loglik, cur_loglik, diff_loglik;
  arma::Col<int> prop_ut_vals;
  arma::colvec idistr;
  arma::colvec new_ut(trun_size);

  for (int i = 0; i < trun_size; ++i) {
    
    idat = dat(i);
    ilag_dep = dep_lag(i);
    cur_ut = ut(i);
    ut_lb = max(0.0, idat - ilag_dep);
    val_len = idat - ut_lb + 1.0;
    idistr.set_size(val_len);
    idistr.fill(1.0 / val_len);
    
    if (ut_lb == idat) {
      prop_ut = idat;
    } else {
      prop_ut_vals = arma::linspace<arma::Col<int>>(ut_lb, idat, val_len);
      prop_ut = as_scalar(Rcpp::RcppArmadillo::sample(prop_ut_vals, 1, 1, idistr));
    }
    
    prop_loglik = R::dbinom(idat - prop_ut, ilag_dep, th, true) + R::dpois(prop_ut, la, true);
    cur_loglik = R::dbinom(idat - cur_ut, ilag_dep, th, true) + R::dpois(cur_ut, la, true);
    diff_loglik = prop_loglik - cur_loglik;
    if (diff_loglik > log(arma::randu())) {
      new_ut(i) = prop_ut;
    } else {
      new_ut(i) = cur_ut;
    }
  }

  return new_ut;
  
}

// update auxiliary variables of the negative binomial MTD model
// [[Rcpp::export]]
arma::colvec updateNBMTDut(const arma::colvec& ut,
                           const arma::colvec& dat,
                           const arma::colvec& dep_lag,
                           const double& th,
                           const double& alp,
                           const double& kap) {
  
  int trun_size = dat.n_rows;
  double ut_lb, val_len, prop_ut, cur_ut, idat, ilag_dep, prop_loglik, cur_loglik, diff_loglik;
  arma::Col<int> prop_ut_vals;
  arma::colvec idistr;
  arma::colvec new_ut(trun_size);
  
  for (int i = 0; i < trun_size; ++i) {
    idat = dat(i);
    ilag_dep = dep_lag(i);
    cur_ut = ut(i);
    ut_lb = max(0.0, idat - ilag_dep);
    if (ut_lb == idat) {
      prop_ut = idat;
    } else {
      val_len = idat - ut_lb + 1.0;
      idistr.set_size(val_len);
      idistr.fill(1.0 / val_len);
      prop_ut_vals = arma::linspace<arma::Col<int>>(ut_lb, idat, val_len);
      prop_ut = as_scalar(Rcpp::RcppArmadillo::sample(prop_ut_vals, 1, 1, idistr));
    }
    prop_loglik = R::dbinom(idat - prop_ut, ilag_dep, th, true) + R::dnbinom(prop_ut, kap + ilag_dep, alp, true);
    cur_loglik = R::dbinom(idat - cur_ut, ilag_dep, th, true) + R::dnbinom(cur_ut, kap + ilag_dep, alp, true);
    diff_loglik = prop_loglik - cur_loglik;
    if (diff_loglik > log(arma::randu())) {
      new_ut(i) = prop_ut;
    } else {
      new_ut(i) = cur_ut;
    }
  }
  
  return new_ut;
  
}

// Adapted from nnmp (X. Zheng), funcs-nnmp-cop.cpp
// [[Rcpp::export]]
List updateCopGMTDlabel(const arma::colvec& data,
                        const arma::mat& lags,
                        const int& mtdorder, 
                        const arma::colvec& weight,
                        const arma::colvec& rho,
                        const int& cop_family,
                        const int& marg_family, 
                        const arma::colvec& yy_param1,
                        const arma::colvec& yy_param2,
                        const arma::mat& yy_ne_param1,
                        const arma::mat& yy_ne_param2) {

  int data_size = data.n_rows;
  int label;
  arma::Col<int> labels(data_size);
  arma::Col<int> label_vals = arma::linspace<arma::Col<int>>(1, mtdorder, mtdorder); 
  arma::colvec dep_lag(data_size), idata(mtdorder);
  arma::colvec ilag, iloglik, idistr;
  arma::colvec iyy_param1, iyy_param2;
  arma::mat mar_param;
  
  iyy_param1.set_size(mtdorder);
  iyy_param2.set_size(mtdorder);
  
  for (int i = 0; i < data_size; ++i){
    
    idata.fill(data(i));
    ilag = lags.row(i).t();
    
    mar_param.set_size(mtdorder, 4);
    iyy_param1.fill(yy_param1(i));
    iyy_param2.fill(yy_param2(i));
    mar_param.col(0) = iyy_param1;
    mar_param.col(1) = iyy_param2;
    mar_param.col(2) = yy_ne_param1.row(i).t();
    mar_param.col(3) = yy_ne_param2.row(i).t();

    iloglik = dBiCopMar_cpp2(idata, ilag, cop_family, rho, marg_family, mar_param, false, true);
    idistr = log(weight) + iloglik;                                  
    idistr = exp(idistr - idistr.max());
    idistr = idistr / sum(idistr);
    label = as_scalar(Rcpp::RcppArmadillo::sample(label_vals, 1, 1, idistr));
    labels(i) = label;
    dep_lag(i) = ilag(label - 1);
  }

  return List::create(Named("data_label") = labels, Named("dep_lag") = dep_lag);

}

// Adapted from nnmp (X. Zheng), funcs-nnmp-cop.cpp
// [[Rcpp::export]]
List updateZIGMTDlabel(const arma::colvec& data,
                       const arma::mat& lags,
                       const int& mtdorder, 
                       const arma::colvec& weight,
                       const arma::colvec& rho,
                       const int& cop_family,
                       const int& marg_family, 
                       const arma::colvec& yy_param1,
                       const arma::colvec& yy_param2,
                       const arma::mat& yy_ne_param1,
                       const arma::mat& yy_ne_param2,
                       const arma::colvec& yy_param3,
                       const arma::colvec& yy_param4,
                       const arma::mat& yy_ne_param3,
                       const arma::mat& yy_ne_param4) {

  int data_size = data.n_rows;
  int label;
  arma::Col<int> labels(data_size);
  arma::Col<int> label_vals = arma::linspace<arma::Col<int>>(1, mtdorder, mtdorder); 
  arma::colvec dep_lag(data_size), idata(mtdorder);
  arma::colvec ilag, iloglik, idistr;
  arma::colvec iyy_param1, iyy_param2, iyy_param3, iyy_param4;
  arma::mat mar_param;
  
  iyy_param1.set_size(mtdorder);
  iyy_param2.set_size(mtdorder);
  iyy_param3.set_size(mtdorder);
  iyy_param4.set_size(mtdorder);
  
  for (int i = 0; i < data_size; ++i){
    
    idata.fill(data(i));
    ilag = lags.row(i).t();
    
    mar_param.set_size(mtdorder, 8);
    iyy_param1.fill(yy_param1(i));
    iyy_param2.fill(yy_param2(i));
    iyy_param3.fill(yy_param3(i));
    iyy_param4.fill(yy_param4(i));
    mar_param.col(0) = iyy_param1;
    mar_param.col(1) = iyy_param2;
    mar_param.col(2) = yy_ne_param1.row(i).t();
    mar_param.col(3) = yy_ne_param2.row(i).t();
    mar_param.col(4) = iyy_param3;
    mar_param.col(5) = iyy_param4;
    mar_param.col(6) = yy_ne_param3.row(i).t();
    mar_param.col(7) = yy_ne_param4.row(i).t();

    iloglik = dBiCopMar_cpp2(idata, ilag, cop_family, rho, marg_family, mar_param, false, true);
    idistr = log(weight) + iloglik;                                  
    idistr = exp(idistr - idistr.max());
    idistr = idistr / sum(idistr);
    label = as_scalar(Rcpp::RcppArmadillo::sample(label_vals, 1, 1, idistr));
    labels(i) = label;
    dep_lag(i) = ilag(label - 1);
  }

  return List::create(Named("data_label") = labels, Named("dep_lag") = dep_lag);

}

// update labels of the Gaussian MTD model
// [[Rcpp::export]]
List updateGMTDlabel(const arma::colvec& data, 
                     const arma::mat& lags, 
                     const int& mtdorder, 
                     const arma::colvec& weight, 
                     const arma::colvec& mu, 
                     const arma::colvec& sigma2,
                     const arma::colvec& rho) {

  int data_size = data.n_rows;
  int label;
  arma::Col<int> label_vals = arma::linspace<arma::Col<int>>(1, mtdorder, mtdorder); 
  arma::Col<int> labels(data_size);
  arma::colvec dep_lag(data_size), idata(mtdorder);
  arma::colvec ilag, imu, isigma2, iloglik, idistr;
  
  for (int i = 0; i < data_size; ++i){
    idata.fill(data(i));
    ilag = lags.row(i).t();
    imu = mu + rho % (ilag - mu);
    isigma2 = sigma2 % (1.0 - arma::pow(rho, 2));
    iloglik = log(arma::normpdf(idata, imu, arma::sqrt(isigma2)));
    idistr = log(weight) + iloglik;
    idistr = exp(idistr - idistr.max());
    idistr = idistr / sum(idistr);
    label = as_scalar(Rcpp::RcppArmadillo::sample(label_vals, 1, 1, idistr));
    labels(i) = label;
    dep_lag(i) = ilag(label - 1);
  }

  return List::create(Named("data_label") = labels, Named("dep_lag") = dep_lag);

}

// update labels of the Poisson MTD model
// [[Rcpp::export]]
List updatePMTDlabel(const arma::colvec& data, 
                     const arma::mat& lags, 
                     const int& mtdorder, 
                     const arma::colvec& weight, 
                     const double& th, 
                     const arma::colvec& ut) {

  int data_size = data.n_rows;
  int label, idata;
  arma::Col<int> label_vals = arma::linspace<arma::Col<int>>(1, mtdorder, mtdorder); 
  arma::Col<int> labels(data_size);
  arma::colvec dep_lag(data_size), iloglik(mtdorder);
  arma::colvec ilag, idistr;
  
  for (int i = 0; i < data_size; ++i){
    idata = data(i) - ut(i);
    ilag = lags.row(i).t();
    iloglik.fill(NA_REAL);
    for (int j = 0; j < mtdorder; ++j){
      iloglik(j) = R::dbinom(idata, ilag(j), th, true);
    }
    idistr = log(weight) + iloglik;
    idistr = exp(idistr - idistr.max());
    idistr = idistr / sum(idistr);
    label = as_scalar(Rcpp::RcppArmadillo::sample(label_vals, 1, 1, idistr));
    labels(i) = label;
    dep_lag(i) = ilag(label - 1);
  }

  return List::create(Named("data_label") = labels, Named("dep_lag") = dep_lag);

}

// update labels of the negative binomial MTD model
// [[Rcpp::export]]
List updateNBMTDlabel(const arma::colvec& data, 
                      const arma::mat& lags, 
                      const int& mtdorder, 
                      const arma::colvec& weight, 
                      const double& th, 
                      const arma::colvec& ut,
                      const double& alp,
                      const double& kap) {

  int data_size = data.n_rows;
  int label, idata, iut;
  arma::Col<int> label_vals = arma::linspace<arma::Col<int>>(1, mtdorder, mtdorder); 
  arma::Col<int> labels(data_size);
  arma::colvec dep_lag(data_size), iloglik(mtdorder);
  arma::colvec ilag, idistr;
  
  for (int i = 0; i < data_size; ++i){
    iut = ut(i);
    idata = data(i) - iut;
    ilag = lags.row(i).t();
    iloglik.fill(NA_REAL);
    for (int j = 0; j < mtdorder; ++j){
      iloglik(j) = R::dbinom(idata, ilag(j), th, true) + R::dnbinom(iut, kap + ilag(j), alp, true);
    }
    idistr = log(weight) + iloglik;
    idistr = exp(idistr - idistr.max());
    idistr = idistr / sum(idistr);
    label = as_scalar(Rcpp::RcppArmadillo::sample(label_vals, 1, 1, idistr));
    labels(i) = label;
    dep_lag(i) = ilag(label - 1);
  }

  return List::create(Named("data_label") = labels, Named("dep_lag") = dep_lag);

}

// update labels of the Lomax MTD regression model
// [[Rcpp::export]]
List updateLMTDlabel(const arma::colvec& data, 
                     const arma::mat& lags, 
                     const int& mtdorder, 
                     const arma::colvec& weight, 
                     const arma::colvec& alpha, 
                     const arma::colvec& phi) {

  int data_size = data.n_rows;
  int label;
  arma::Col<int> label_vals = arma::linspace<arma::Col<int>>(1, mtdorder, mtdorder); 
  arma::Col<int> labels(data_size);
  arma::colvec dep_lag(data_size), idata(mtdorder);
  arma::colvec ilag, iloglik, idistr;

  for (int i = 0; i < data_size; ++i){
    idata.fill(data(i));
    ilag = lags.row(i).t();
    iloglik = log(alpha) + alpha % log(phi + ilag) - (alpha + 1.0) % log(phi + ilag + idata);
    idistr = log(weight) + iloglik;
    idistr = exp(idistr - idistr.max());
    idistr = idistr / sum(idistr);
    label = as_scalar(Rcpp::RcppArmadillo::sample(label_vals, 1, 1, idistr));
    labels(i) = label;
    dep_lag(i) = ilag(label - 1);
  }

  return List::create(Named("data_label") = labels, Named("dep_lag") = dep_lag);

}

// one-step ahead posterior predictive distribution of the Poisson MTD model
// [[Rcpp::export]]
arma::mat predPMTD(const int& mtdorder, 
                   const arma::mat& weight, 
                   const arma::mat& lags,
                   const arma::colvec& la, 
                   const arma::colvec& th,
                   const arma::colvec& probs) {

  int sample_size = weight.n_cols;
  int trun_size = lags.n_rows;
  int probs_size = probs.n_rows;
  int label;
  double iter_la, iter_th, iter_ut, dep_lag;
  arma::mat preds(probs_size, trun_size);
  arma::colvec ipred(sample_size);
  arma::Col<int> label_vals = arma::linspace<arma::Col<int>>(0, mtdorder - 1, mtdorder);
  arma::colvec ilags, iter_weight;

  for (int i = 0; i < trun_size; ++ i) {
    
    ilags = lags.row(i).t();
    
    for (int iter = 0; iter < sample_size; ++iter) {
      iter_la = la(iter);
      iter_th = th(iter);
      iter_weight = weight.col(iter);
      label = arma::as_scalar(Rcpp::RcppArmadillo::sample(label_vals, 1, 1, iter_weight));
      dep_lag = ilags(label);
      iter_ut = R::rpois(iter_la);
      ipred(iter) = iter_ut + R::rbinom(dep_lag, iter_th);
    }
    
    preds.col(i) = arma::quantile(ipred, probs);
    
  }

  return preds;

}

// one-step ahead posterior predictive distribution of the negative binomial MTD model
// [[Rcpp::export]]
arma::mat predNBMTD(const int& mtdorder, 
                    const arma::mat& weight, 
                    const arma::mat& lags,
                    const arma::colvec& th, 
                    const arma::colvec& psi,
                    const arma::colvec& kap,
                    const arma::colvec& probs) {
  
  int sample_size = weight.n_cols;
  int trun_size = lags.n_rows;
  int probs_size = probs.n_rows;
  int label;
  double iter_th, iter_psi, iter_kap, iter_ut, dep_lag;
  arma::mat preds(probs_size, trun_size);
  arma::colvec ipred(sample_size);
  arma::Col<int> label_vals = arma::linspace<arma::Col<int>>(0, mtdorder - 1, mtdorder);
  arma::colvec ilags, iter_weight;
  
  for (int i = 0; i < trun_size; ++ i) {
    
    ilags = lags.row(i).t();
    
    for (int iter = 0; iter < sample_size; ++iter) {
      
      iter_th = th(iter);
      iter_psi = psi(iter);
      iter_kap = kap(iter);
      iter_weight = weight.col(iter);
      label = arma::as_scalar(Rcpp::RcppArmadillo::sample(label_vals, 1, 1, iter_weight));
      dep_lag = ilags(label);
      iter_ut = R::rnbinom(iter_kap + dep_lag, iter_psi);
      ipred(iter) = iter_ut + R::rbinom(dep_lag, iter_th);
      
    }
    
    preds.col(i) = arma::quantile(ipred, probs);
    
  }
  
  return preds;
  
}

// [[Rcpp::export]]
arma::mat predCopGammaMTD(const int& mtdorder, 
                          const arma::mat& weight, 
                          const arma::mat& lags,
                          const arma::colvec& alp,
                          const arma::colvec& be,
                          const int& cop_family,
                          const arma::mat& rho,
                          const int& marg_family,                     
                          const arma::colvec& probs) {

  int sample_size = weight.n_cols;
  int trun_size = lags.n_rows;
  int probs_size = probs.n_rows;
  int label;
  double iter_alp, iter_be, dep_lag, iter_vv, iter_uu, iter_cop_param; 
  arma::mat preds(probs_size, trun_size);
  arma::colvec ipred(sample_size);
  arma::Col<int> label_vals = arma::linspace<arma::Col<int>>(0, mtdorder - 1, mtdorder);
  arma::colvec ilags, iter_weight, iter_rho;

  for (int i = 0; i < trun_size; ++ i) {

    ilags = lags.row(i).t();

    for (int iter = 0; iter < sample_size; ++iter) {

      iter_rho = rho.col(iter);
      iter_alp = alp(iter);
      iter_be = be(iter);
      iter_weight = weight.col(iter);
      
      label = arma::as_scalar(Rcpp::RcppArmadillo::sample(label_vals, 1, 1, iter_weight));
      dep_lag = ilags(label);
      iter_cop_param = iter_rho(label);
      
      iter_vv = pgamma(dep_lag, iter_alp, iter_be);
      iter_uu = rConCop(iter_vv, cop_family, iter_cop_param); 
      ipred(iter) = qgamma(iter_uu, iter_alp, iter_be);       

    }

    preds.col(i) = arma::quantile(ipred, probs);

  }

  return preds; 

}

// [[Rcpp::export]]
arma::mat predCopGammaMTD_2(const int& mtdorder, 
                            const arma::mat& weight, 
                            const arma::mat& lags,
                            const arma::colvec& alp,
                            const arma::colvec& be,
                            const int& cop_family,
                            const arma::mat& rho,
                            const int& marg_family,                     
                            const arma::colvec& probs) {

  int sample_size = weight.n_cols;
  int trun_size = lags.n_rows;
  int probs_size = probs.n_rows;
  int label;
  double iter_alp, iter_be, dep_lag, iter_vv, iter_uu, iter_cop_param; 
  arma::colvec pred_mean(trun_size);
  arma::mat preds(probs_size, trun_size);
  arma::colvec ipred(sample_size);
  arma::Col<int> label_vals = arma::linspace<arma::Col<int>>(0, mtdorder - 1, mtdorder);
  arma::colvec ilags, iter_weight, iter_rho;

  for (int i = 0; i < trun_size; ++ i) {

    ilags = lags.row(i).t();

    for (int iter = 0; iter < sample_size; ++iter) {

      iter_rho = rho.col(iter);
      iter_alp = alp(iter);
      iter_be = be(iter);
      iter_weight = weight.col(iter);

      label = arma::as_scalar(Rcpp::RcppArmadillo::sample(label_vals, 1, 1, iter_weight));
      dep_lag = ilags(label);
      iter_cop_param = iter_rho(label);

      iter_vv = pgamma(dep_lag, iter_alp, iter_be);
      iter_uu = rConCop(iter_vv, cop_family, iter_cop_param); 
      ipred(iter) = qgamma(iter_uu, iter_alp, iter_be);       

    }

    preds.col(i) = arma::quantile(ipred, probs);
    pred_mean(i) = mean(ipred);

  }

  return pred_mean; 

}

// [[Rcpp::export]]
arma::mat predCopZIGammaMTD(const int& mtdorder, 
                            const arma::mat& weight, 
                            const arma::mat& lags,
                            const arma::colvec& alp,
                            const arma::colvec& be,
                            const arma::colvec& Pi,
                            const arma::colvec& eps,
                            const int& cop_family,
                            const arma::mat& rho,
                            const int& marg_family,                     
                            const arma::colvec& probs) {

  int sample_size = weight.n_cols;
  int trun_size = lags.n_rows;
  int probs_size = probs.n_rows;
  int label;
  double iter_alp, iter_be, iter_Pi, iter_eps, dep_lag, iter_vv, iter_uu, iter_cop_param;  
  arma::mat preds(probs_size, trun_size);
  arma::colvec ipred(sample_size);
  arma::Col<int> label_vals = arma::linspace<arma::Col<int>>(0, mtdorder - 1, mtdorder);
  arma::colvec ilags, iter_weight, iter_rho;

  for (int i = 0; i < trun_size; ++ i) {

    ilags = lags.row(i).t();

    for (int iter = 0; iter < sample_size; ++iter) {

      iter_rho = rho.col(iter);
      iter_alp = alp(iter);
      iter_be = be(iter);
      iter_Pi = Pi(iter);
      iter_eps = eps(iter);
      iter_weight = weight.col(iter);

      label = arma::as_scalar(Rcpp::RcppArmadillo::sample(label_vals, 1, 1, iter_weight));
      dep_lag = ilags(label);
      iter_cop_param = iter_rho(label);

      iter_vv = pzig(dep_lag, iter_alp, iter_be, iter_eps, iter_Pi);
      iter_uu = rConCop(iter_vv, cop_family, iter_cop_param);
      ipred(iter) = qzig(iter_uu, iter_alp, iter_be, iter_Pi, iter_eps);
    
    }

    preds.col(i) = arma::quantile(ipred, probs);

  }

  return preds;

}

// [[Rcpp::export]]
arma::mat predCopZIGammaMTD_2(const int& mtdorder, 
                              const arma::mat& weight, 
                              const arma::mat& lags,
                              const arma::colvec& alp,
                              const arma::colvec& be,
                              const arma::colvec& Pi,
                              const arma::colvec& eps,
                              const int& cop_family,
                              const arma::mat& rho,
                              const int& marg_family,                     
                              const arma::colvec& probs) {

  int sample_size = weight.n_cols;
  int trun_size = lags.n_rows;
  int probs_size = probs.n_rows;
  int label;
  double iter_alp, iter_be, iter_Pi, iter_eps, dep_lag, iter_vv, iter_uu, iter_cop_param;  
  arma::colvec pred_mean(trun_size);
  arma::mat preds(probs_size, trun_size);
  arma::colvec ipred(sample_size);
  arma::Col<int> label_vals = arma::linspace<arma::Col<int>>(0, mtdorder - 1, mtdorder);
  arma::colvec ilags, iter_weight, iter_rho;

  for (int i = 0; i < trun_size; ++ i) {

    ilags = lags.row(i).t();

    for (int iter = 0; iter < sample_size; ++iter) {

      iter_rho = rho.col(iter);
      iter_alp = alp(iter);
      iter_be = be(iter);
      iter_Pi = Pi(iter);
      iter_eps = eps(iter);
      iter_weight = weight.col(iter);

      label = arma::as_scalar(Rcpp::RcppArmadillo::sample(label_vals, 1, 1, iter_weight));
      dep_lag = ilags(label);
      iter_cop_param = iter_rho(label);

      iter_vv = pzig(dep_lag, iter_alp, iter_be, iter_eps, iter_Pi);
      iter_uu = rConCop(iter_vv, cop_family, iter_cop_param);
      ipred(iter) = qzig(iter_uu, iter_alp, iter_be, iter_Pi, iter_eps);

    }

    preds.col(i) = arma::quantile(ipred, probs);
    pred_mean(i) = mean(ipred);
  
  }

  return pred_mean; 

}


