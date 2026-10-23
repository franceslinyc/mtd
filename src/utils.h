#ifndef UTILS_H
#define UTILS_H
#include <RcppArmadillo.h>

double qbeta(double uu, double shape1, double shape2);

double pbeta(double uu, double shape1, double shape2);

double qgamma(double uu, double shape1, double shape2);

double pgamma(double uu, double shape1, double shape2);
  
#endif 
