#include "BaseFunc.hpp"


// [[Rcpp::export]]
double test_getLogDist(NumericVector& x, NumericVector& y, double s)
{
  return getLogDist(x, y, s);
}

// [[Rcpp::export]]
NumericMatrix testInverse(NumericMatrix& x)
{
  MapMatd y(as<MapMatd> (x));
  return wrap(y.inverse());
}

// [[Rcpp::export]]
NumericVector test_getLogDistVector(NumericMatrix& x, NumericVector& y, double s)
{
  return getLogDistVector(x, y, s);
}

/***R
# x <- c(1, 2, 3)
# y <- c(4, 5, 6)
# test_getLogDist(x, y, 0)
# x <- matrix(c(7, 10, 15, 22), nrow = 2, ncol = 2)
# x %*% testInverse(x)
x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3, byrow = T)
y <- c(4, 5, 6)
test_getLogDistVector(x, y, 1)
*/
