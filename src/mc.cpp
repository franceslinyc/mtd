#include <iostream>
#include <RcppArmadillo.h>
#include <math.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double dCondBiPois(const int& uu,
                   const int& vv,
                   const double& th1,
                   const double& th2,
                   const double& th12) {
  
  int zz = min(uu, vv);                    
  double logp1 = log(th12) - log(th2 + th12);
  double logp2 = log(th2) - log(th2 + th12);
  double cc1, cc2, dens;
  arma::Col<int> grid = arma::linspace<arma::Col<int>>(0, zz, zz + 1); 
  arma::colvec yy(zz + 1);
  
  for (int j = 0; j < zz + 1; ++j) {
    cc1 = lgamma(vv + 1) - lgamma(j + 1) - lgamma(vv - j + 1);
    cc2 = j * logp1 + (vv - j) * logp2 + (uu - j) * log(th1) - lgamma(uu - j + 1);
    yy(j) = exp(cc1 + cc2);
  }
  
  dens = exp(-th1) * sum(yy);
  
  return dens;
  
}

// [[Rcpp::export]]
double pCondBiPois(const int& uu,
                   const int& vv,
                   const double& th1,
                   const double& th2,
                   const double& th12) {
  
  arma::colvec yy(uu + 1);  
  
  for (int i = 0; i < uu + 1; ++i) {
    yy(i) = dCondBiPois(i, vv, th1, th2, th12);
  }
  
  return sum(yy);
  
}

// [[Rcpp::export]]
double pPMTD(const int& obs,
             const arma::Col<int>& lags,
             const arma::colvec& weight,
             const double& la,
             const double& ga) {
  
  int mtdorder = lags.n_rows;
  arma::colvec yy(mtdorder);
  
  for (int j = 0; j < mtdorder; ++j) {
    yy(j) = weight(j) * pCondBiPois(obs, lags(j), la, la, ga);
  }
  
  return sum(yy);
  
}

// [[Rcpp::export]]
arma::mat rqrPMTD(const arma::Col<int>& obs,
                  const arma::mat& weight,
                  const arma::colvec& la,
                  const arma::colvec& ga) {
  
  int mtdorder = weight.n_rows;
  int trun_size = obs.n_rows - mtdorder;
  int sample_size = la.n_rows; 
  int iobs;
  double iter_la, iter_ga, at, bt, pt;
  arma::colvec iter_weight;
  arma::Col<int> ilags_rev, ilags;
  arma::mat rt(trun_size, sample_size);
  
  for (int iter = 0; iter < sample_size; ++iter) {
    
    iter_la = la(iter);
    iter_ga = ga(iter);
    iter_weight = weight.col(iter);
    
    for (int i = 0; i < trun_size; ++i) {
      
      iobs = obs(i + mtdorder);
      ilags_rev = obs.rows(i, i + mtdorder - 1);
      ilags = arma::reverse(ilags_rev);
      bt = pPMTD(iobs, ilags, iter_weight, iter_la, iter_ga);
      at = pPMTD(iobs - 1, ilags, iter_weight, iter_la, iter_ga);
      pt = R::runif(at, bt);
      rt(i, iter) = R::qnorm(pt, 0.0, 1.0, true, false);
      
    }
    
  }
  
  return rt;
  
}

// [[Rcpp::export]]
double pLMTD(const double& obs,
             const arma::colvec& lags,
             const arma::colvec& weight,
             const double& alpha,
             const double& phi) {
  
  int mtdorder = lags.n_rows;
  arma::colvec yy(mtdorder);
  
  for (int j = 0; j < mtdorder; ++j) {
    yy(j) = weight(j) * (1.0 - pow((1.0 + obs / (lags(j) + phi)), -alpha));
  }
  
  return sum(yy);
  
}

// [[Rcpp::export]]
arma::mat rqrRegLMTD(const arma::colvec& obs,
                     const arma::mat& xx,
                     const arma::mat& weight,
                     const arma::mat& bb, 
                     const arma::colvec& alpha,
                     const arma::colvec& phi) {
  
  int mtdorder = weight.n_rows;
  int trun_size = obs.n_rows - mtdorder;
  int sample_size = alpha.n_rows;
  arma::rowvec ixx;
  arma::colvec eps, ilags_rev, ilags, ipp, iter_bb, iter_weight;
  double iEps, iter_alpha, iter_phi, pt;
  arma::mat rt(trun_size, sample_size);
  
  for (int iter = 0; iter < sample_size; ++iter) {
    
    iter_bb = bb.col(iter);
    iter_alpha = alpha(iter);
    iter_phi = phi(iter);
    iter_weight = weight.col(iter);
    eps = obs % arma::exp(-xx * iter_bb);
    
    for (int i = 0; i < trun_size; ++i) {
      iEps = eps(i + mtdorder);
      ilags_rev = eps.rows(i, i + mtdorder - 1);
      ilags = arma::reverse(ilags_rev);
      pt = pLMTD(iEps, ilags, iter_weight, iter_alpha, iter_phi);
      rt(i, iter) = R::qnorm(pt, 0.0, 1.0, true, false);
    }
    
  }
  
  return rt;
  
}

// Adapted from pbv (Alexander Robitzsch), GPL (>= 2): pbv_rcpp_bvnorm.cpp.
// [[Rcpp::export]]
double pbv_rcpp_pnorm0(double z) {
  return ::Rf_pnorm5(z, 0.0, 1.0, 1, 0);
}

// Adapted from pbv (Alexander Robitzsch), GPL (>= 2): pbv_rcpp_bvnorm.cpp.
// [[Rcpp::export]]
double pbv_rcpp_pbvnorm0(double h1, double hk, double r) {
  int nx = 5;
  Rcpp::NumericVector x(nx);
  Rcpp::NumericVector w(nx);
  x[0] = .04691008;
  x[1] = .23076534;
  x[2] = .5;
  x[3] = .76923466;
  x[4] = .95308992;
  w[0] = .018854042;
  w[1] = .038088059;
  w[2] = .0452707394;
  w[3] = .038088059;
  w[4] = .018854042;

  double bv = 0;
  double r1, r2, rr, rr2, r3, h3, h5, h6, h7, aa, ab, h11;
  double cor_max = 0.7;
  double bv_fac1 = 0.13298076;
  double bv_fac2 = 0.053051647;
  double h2 = hk;
  double h12 = (h1 * h1 + h2 * h2) / 2;
  double r_abs = std::abs(r);

  if (r_abs > cor_max) {
    r2 = 1.0 - r * r;
    r3 = std::sqrt(r2);
    if (r < 0) h2 = -h2;
    h3 = h1 * h2;
    h7 = std::exp(-h3 / 2.0);
    if (r_abs < 1) {
      h6 = std::abs(h1 - h2);
      h5 = h6 * h6 / 2.0;
      h6 = h6 / r3;
      aa = 0.5 - h3 / 8.0;
      ab = 3.0 - 2.0 * aa * h5;
      bv = bv_fac1 * h6 * ab * (1 - pbv_rcpp_pnorm0(h6)) -
        std::exp(-h5 / r2) * (ab + aa * r2) * bv_fac2;
      for (int ii = 0; ii < nx; ii++) {
        r1 = r3 * x[ii];
        rr = r1 * r1;
        r2 = std::sqrt(1.0 - rr);
        bv += -w[ii] * std::exp(-h5 / rr) *
          (std::exp(-h3 / (1.0 + r2)) / r2 / h7 - 1.0 - aa * rr);
      }
    }
    h11 = std::min(h1, h2);
    bv = bv * r3 * h7 + pbv_rcpp_pnorm0(h11);
    if (r < 0) bv = pbv_rcpp_pnorm0(h1) - bv;
  } else {
    h3 = h1 * h2;
    for (int ii = 0; ii < nx; ii++) {
      r1 = r * x[ii];
      rr2 = 1.0 - r1 * r1;
      bv += w[ii] * std::exp((r1 * h3 - h12) / rr2) / std::sqrt(rr2);
    }
    bv = pbv_rcpp_pnorm0(h1) * pbv_rcpp_pnorm0(h2) + r * bv;
  }
  return bv;
}

#define PMAX (1.0 - 1e-10)
#define PMIN 1e-10

double pBiCop_cpp(const double& uu, const double& vv, const int& cop_family,
                  const double& cop_param, const bool& logd = false) {
  double logdens;
  if (cop_family == 1) {
    double xx = R::qnorm(uu, 0.0, 1.0, true, false);
    double yy = R::qnorm(vv, 0.0, 1.0, true, false);
    logdens = pbv_rcpp_pbvnorm0(xx, yy, cop_param);
  } else if (cop_family == 2) {
    double eta = cop_param;
    if (eta < 1) stop("error: gumbel copula parameter must be greater than 1.");
    if (eta > 50.0) eta = 50.0;
    double xx = -log(uu);
    double yy = -log(vv);
    double tt = pow(xx, eta) + pow(yy, eta);
    logdens = -pow(tt, 1.0 / eta) + log(pow(tt, 1.0 / eta) + eta - 1.0) +
      (1.0 / eta - 2.0) * log(tt) + (eta - 1.0) * (log(xx) + log(yy)) + xx + yy;
  } else if (cop_family == 3) {
    double delta = cop_param;
    if (delta > 98) delta = 98;
    double ee = pow(uu, -delta) + pow(vv, -delta) - 1.0;
    logdens = log(1.0 + delta) - (delta + 1.0) * (log(uu) + log(vv)) -
      (2.0 + 1.0 / delta) * log(ee);
  } else {
    stop("error: no families are found.");
  }

  if (cop_family == 1 || logd) return logdens;
  return exp(logdens);
}

arma::colvec pBiCop_cpp2(const arma::colvec& uu, const arma::colvec& vv,
                         const int& cop_family, const arma::colvec& cop_param,
                         const bool& logd = false) {
  int nn = uu.n_rows;
  arma::colvec result(nn);
  for (int i = 0; i < nn; ++i) {
    result(i) = pBiCop_cpp(uu(i), vv(i), cop_family, cop_param(i), logd);
  }
  return result;
}

// [[Rcpp::export]]
double psgamma_v2(double x, double mu, double beta, double eps) {
  return (x >= eps) * R::pgamma(x - eps, mu / beta, beta, true, false);
}

// [[Rcpp::export]]
double pzig_v2(double x, double mu, double beta, double eps, double pi) {
  return pi * R::punif(x, 0, eps, true, false) +
    (1 - pi) * psgamma_v2(x, mu, beta, eps);
}

// [[Rcpp::export]]
double pBiCopMar_cpp(const double& xx, const double& yy, const int& cop_family,
                     const double& cop_param, const int& marg_family,
                     const arma::colvec& mar_param, const bool& ce = false,
                     const bool& logd = false) {
  double uu, vv;
  if (marg_family == 1) {
    uu = R::pbeta(xx, mar_param(0), mar_param(1), true, false);
    vv = R::pbeta(yy, mar_param(2), mar_param(3), true, false);
  } else if (marg_family == 2) {
    uu = R::pgamma(xx, mar_param(0), 1.0 / mar_param(1), true, false);
    vv = R::pgamma(yy, mar_param(2), 1.0 / mar_param(3), true, false);
  } else if (marg_family == 3) {
    uu = R::pnorm(xx, mar_param(0), mar_param(1), true, false);
    vv = R::pnorm(yy, mar_param(2), mar_param(3), true, false);
  } else if (marg_family == 4) {
    uu = pzig_v2(xx, mar_param(0), mar_param(1), mar_param(4), mar_param(5));
    vv = pzig_v2(yy, mar_param(2), mar_param(3), mar_param(6), mar_param(7));
  } else {
    stop("error: no such a family is available.");
  }

  uu = std::max(PMIN, std::min(PMAX, uu));
  vv = std::max(PMIN, std::min(PMAX, vv));
  return pBiCop_cpp(uu, vv, cop_family, cop_param, logd);
}

// [[Rcpp::export]]
arma::colvec pBiCopMar_cpp2(const arma::colvec& xx, const arma::colvec& yy,
                            const int& cop_family, const arma::colvec cop_param,
                            const int& marg_family, const arma::mat& mar_param,
                            const bool& ce = false, const bool& logd = false) {
  int nn = xx.n_rows;
  arma::colvec uu(nn), vv(nn);
  if (marg_family == 1) {
    for (int i = 0; i < nn; ++i) {
      uu(i) = R::pbeta(xx(i), mar_param(i, 0), mar_param(i, 1), true, false);
      vv(i) = R::pbeta(yy(i), mar_param(i, 2), mar_param(i, 3), true, false);
    }
  } else if (marg_family == 2) {
    for (int i = 0; i < nn; ++i) {
      uu(i) = R::pgamma(xx(i), mar_param(i, 0), 1.0 / mar_param(i, 1), true, false);
      vv(i) = R::pgamma(yy(i), mar_param(i, 2), 1.0 / mar_param(i, 3), true, false);
    }
  } else if (marg_family == 3) {
    for (int i = 0; i < nn; ++i) {
      uu(i) = R::pnorm(xx(i), mar_param(i, 0), mar_param(i, 1), true, false);
      vv(i) = R::pnorm(yy(i), mar_param(i, 2), mar_param(i, 3), true, false);
    }
  } else if (marg_family == 4) {
    for (int i = 0; i < nn; ++i) {
      uu(i) = pzig_v2(xx(i), mar_param(i, 0), mar_param(i, 1), mar_param(i, 4), mar_param(i, 5));
      vv(i) = pzig_v2(yy(i), mar_param(i, 2), mar_param(i, 3), mar_param(i, 6), mar_param(i, 7));
    }
  } else {
    stop("error: no such a family is available.");
  }

  for (int i = 0; i < nn; ++i) {
    uu(i) = std::max(PMIN, std::min(PMAX, uu(i)));
    vv(i) = std::max(PMIN, std::min(PMAX, vv(i)));
  }
  return pBiCop_cpp2(uu, vv, cop_family, cop_param, logd);
}

// [[Rcpp::export]]
double pCopGammaMTD(const double& obs, const arma::colvec& lags,
                    const arma::colvec& weight, const double& alp,
                    const double& be, const int& cop_family,
                    const arma::colvec& rho, const int& marg_family,
                    const arma::colvec& mar_param) {
  int mtdorder = lags.n_rows;
  arma::colvec result(mtdorder);
  for (int j = 0; j < mtdorder; ++j) {
    double copula_cdf = pBiCopMar_cpp(obs, lags(j), cop_family, rho(j),
                                     marg_family, mar_param, false, false);
    double lag_cdf = R::pgamma(lags(j), alp, 1.0 / be, true, false);
    result(j) = weight(j) * copula_cdf / lag_cdf;
  }
  return sum(result);
}

// [[Rcpp::export]]
arma::mat rqrCopGammaMTD(const arma::colvec& obs, const arma::mat& weight,
                         const arma::colvec& alp, const arma::colvec& be,
                         const int& cop_family, const arma::mat& rho,
                         const int& marg_family) {
  int mtdorder = weight.n_rows;
  int trun_size = obs.n_rows - mtdorder;
  int sample_size = alp.n_rows;
  arma::colvec ilags_rev, ilags, iter_weight, iter_rho, iter_mar;
  double iobs, iter_alp, iter_be, pt;
  arma::mat rt(trun_size, sample_size);

  for (int iter = 0; iter < sample_size; ++iter) {
    iter_rho = rho.col(iter);
    iter_alp = alp(iter);
    iter_be = be(iter);
    iter_mar.set_size(4);
    iter_mar(0) = iter_alp;
    iter_mar(1) = iter_be;
    iter_mar(2) = iter_alp;
    iter_mar(3) = iter_be;
    iter_weight = weight.col(iter);

    for (int i = 0; i < trun_size; ++i) {
      iobs = obs(i + mtdorder);
      ilags_rev = obs.rows(i, i + mtdorder - 1);
      ilags = arma::reverse(ilags_rev);
      pt = pCopGammaMTD(iobs, ilags, iter_weight, iter_alp, iter_be,
                        cop_family, iter_rho, marg_family, iter_mar);
      rt(i, iter) = R::qnorm(pt, 0.0, 1.0, true, false);
    }
  }
  return rt;
}
