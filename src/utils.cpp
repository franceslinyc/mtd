#include <iostream>
#include <RcppArmadillo.h>
#include <math.h>

using namespace std;
using namespace Rcpp;

// Adapted from nnmp (X. Zheng), utils.cpp
#define EPSILON DBL_EPSILON

double qbeta(double uu, double shape1, double shape2) {
  return R::qbeta(uu, shape1, shape2, true, false);
}

double pbeta(double uu, double shape1, double shape2) {
  return R::pbeta(uu, shape1, shape2, true, false);
}

double qgamma(double uu, double shape1, double shape2) {
  return R::qgamma(uu, shape1, 1.0 / shape2, true, false);
}

double pgamma(double uu, double shape1, double shape2) {
  return R::pgamma(uu, shape1, 1.0 / shape2, true, false);
}
