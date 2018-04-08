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

// [[Rcpp::export]]
NumericMatrix test_cbindV(NumericMatrix& x, NumericVector& y)
{
  return cbindV(x, y);
}

// [[Rcpp::export]]
NumericVector test_sum(NumericMatrix& x)
{
  return rowSums(x);
}

// [[Rcpp::export]]
NumericMatrix test_pdist(NumericMatrix& x, NumericMatrix& y)
{
  return fastpdist(x, y);
}


// [[Rcpp::export]]
NumericMatrix test_rowErase(NumericMatrix& x, int rowID)
{
  x = rowErase(x, rowID);
  return x;
}

// [[Rcpp::export]]
NumericMatrix test_colErase(NumericMatrix& x, int colID)
{
  return colErase(x, colID);
}

// [[Rcpp::export]]
int test_reportMaxPrime(int n)
{
  return reportMaxPrime(n);
}

// [[Rcpp::export]]
NumericMatrix test_varCPP(NumericMatrix& x)
{
  return varCPP(x);
}

// [[Rcpp::export]]
NumericMatrix test_sqrtVarMatrix(NumericMatrix& x)
{
  return sqrtVarMatrix(x);
}


// [[Rcpp::export]]
NumericVector test_quantileCPP(NumericVector& x, NumericVector& q)
{
  return quantileCPP(x, q);
}

// [[Rcpp::export]]
NumericVector test_sampleCPP(int dim)
{
  return sampleCPP(dim);
}

// [[Rcpp::export]]
NumericMatrix test_subMatrixCols(NumericMatrix& x, NumericVector& cols)
{
  return subMatrixCols(x, cols);
}

// [[Rcpp::export]]
IntegerVector test_orderCPP(NumericVector x)
{
  return orderCPP(x);
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
# x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3, byrow = T)
# test_rowMin(x)
# test_compareMin(x[1, ], x[2, ])

# x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3, byrow = T)
# y <- c(10, 15, 11)
# z <- test_cbindV(x, y)
# y <- c(1, -100, -50)
# print(z)
# print(y)
# x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3, byrow = T)
# y <- matrix(c(1, 2, 3, 4, 3, 6, 0, 6, 9), nrow = 3, ncol = 3, byrow = T)
# test_pdist(x, y)
# test_rowErase(x, 0)
# test_colErase(test_rowErase(x, 0), 0)
# test_reportMaxPrime(110)

# x <- MaxPro::MaxProLHD(20, 2)$Design
# varMatrix <- test_varCPP(x)
# test_sqrtVarMatrix(varMatrix)
# eig <- eigen(varMatrix)
# sqrtSinv <- eig$vec%*%diag(1/sqrt(eig$val))%*%t(eig$vectors)
# print(sqrtSinv)
*/
