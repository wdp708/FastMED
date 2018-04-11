#include "BaseFunc.hpp"

// [[Rcpp::export]]
int test_which_max(NumericVector& x)
{
  return which_max(x);
}

// [[Rcpp::export]]
NumericMatrix test_fastpdist(NumericMatrix& x, NumericMatrix& y)
{
  return fastpdist(x, y);
}

// [[Rcpp::export]]
NumericVector test_fastpdist2(NumericMatrix& x, NumericVector& y)
{
  return fastpdist2(x, y);
}

// [[Rcpp::export]]
NumericMatrix test_rbindNullMatrix(NumericMatrix& x)
{
  NumericMatrix res;
  return rbindM(res, x);
}

/***R
 test_which_max(c(1, 5, 6, 8, 1, 2))
 */

