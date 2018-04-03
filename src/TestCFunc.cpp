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

// [[Rcpp::export]]
int testAssign(NumericMatrix& x, NumericMatrix& y)
{
  cout<<x<<endl;
  cout<<y<<endl;
  x(_, 1) = y(_, 1);
  x(_, 1) = y(_, 2);
  cout<<x<<endl;
  cout<<y<<endl;
  MatrixXd m_test(x.nrow(), x.ncol());
  NumericVector tempRow;
  tempRow = x(_, 1);
  m_test.col(1) = as<MapVectd> (tempRow);
  tempRow = y(_, 1);
  m_test.col(1) = as<MapVectd> (tempRow);
  cout<<m_test<<endl;
  return 0;
}

// [[Rcpp::export]]
NumericMatrix test_transMatrix(MatrixXd& x, NumericMatrix& y)
{
  return transMatrix(x, y);
}


// [[Rcpp::export]]
NumericVector test_transVector(MatrixXd& x, NumericVector& y)
{
  return transVector(x, y);
}

// [[Rcpp::export]]
NumericVector test_rowMin(NumericMatrix& x)
{
  return rowMin(x);
}

// [[Rcpp::export]]
NumericVector test_compareMin(NumericVector& x, NumericVector& y)
{
  return compareMin(x, y);
}

/***R
# x <- c(1, 2, 3)
# y <- c(4, 5, 6)
# test_getLogDist(x, y, 0)
# x <- matrix(c(7, 10, 15, 22), nrow = 2, ncol = 2)
# x %*% testInverse(x)
# x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3, byrow = T)
# y <- c(4, 5, 6)
# test_getLogDistVector(x, y, 1)
# x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3, byrow = T)
# y <- matrix(c(1, 2, 3, 4, 3, 6, 0, 6, 9), nrow = 3, ncol = 3, byrow = T)
# testAssign(x, y)
# x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3, byrow = T)
# y <- matrix(c(1, 2, 3, 4, 3, 6, 0, 6, 9), nrow = 3, ncol = 3, byrow = T)
# test_transMatrix(x, y)
# test_transVector(x, y[1,])
x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3, byrow = T)
test_rowMin(x)
test_compareMin(x[1, ], x[2, ])
*/
