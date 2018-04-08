#ifndef BASEFUNC__
#define BASEFUNC__


#include <Rcpp.h>
#include <math.h>
#include <unordered_set>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace RcppEigen;
using namespace Eigen;
using namespace std;
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<Eigen::VectorXd> MapVectd;
// [[Rcpp::plugins(cpp11)]]


#define minTwo(a, b) (a<b? a:b)
#define MaxTwo(a, b) (a>b? a:b)


/*
 * an alternative function is sugar::cumprod
 */
// [[Rcpp::export]]
NumericVector cumprodC(NumericVector& x)
{
  NumericVector res(x.size());
  double temp = 1.0;
  for(int i = 0; i < x.size(); i++)
  {
    temp *= x[i];
    res[i] = temp;
  }
  return res;
}

/*
 * an alternative functin is sugar::which_min
 */
// [[Rcpp::export]]
int argMin(NumericVector& x)
{
  return std::min_element(x.begin(), x.end()) - x.begin();
}

/*
 * an alternative function is sugar::which_max
 */
// [[Rcpp::export]]
int argMax(NumericVector& x)
{
  return std::max_element(x.begin(), x.end()) - x.begin();
}

/*
 * compute the distance between two vectors
 */
// [[Rcpp::export]]
double getLogDist(NumericVector& x, NumericVector& y, double s)
{
  double res(1.0);
  int n(x.size());
  try
  {
    if(n != y.size())
    {
      throw "vectors must have the same length.";
    }
    if(s > 1e-6)
    {
      NumericVector temp = pow(abs(x - y), s);
      res = accumulate(temp.begin(), temp.end(), 0.0);
      res = log(res / n ) / s;
    }
    else
    {
      NumericVector temp = log(abs(x - y));
      res = 1.0  / n * accumulate(temp.begin(), temp.end(), 0.0);
    }
    return res;
  }catch(const char* e)
  {
    cout<<e<<endl;
    return -1.0;
  }
}

// [[Rcpp::export]]
NumericVector getLogDistVector(NumericMatrix& x, NumericVector& y, double s)
{
  NumericVector res(x.nrow(), -1.0);
  int n(y.size());
  try
  {
    if(x.ncol() != n)
    {
      throw "The columns of x must be same as the size of vector y.";
    }
    NumericVector tempRow;
    for(int i = 0; i < x.nrow(); i++)
    {
      tempRow = x(i, _);
      res[i] = getLogDist(tempRow, y, s);
    }
    return res;
  }catch(const char* e)
  {
    cout<< e <<endl;
    return res;
  }
}

// [[Rcpp::export]]
NumericMatrix transMatrix(MatrixXd& SigmaMatrix, NumericMatrix& x)
{
  MapMatd tempX(as<MapMatd> (x));
  MatrixXd res(x.nrow(), x.ncol());
  res = tempX * SigmaMatrix;
  return wrap(res);
}

// [[Rcpp::export]]
NumericVector transVector(MatrixXd& SigmaMatrix, NumericVector& x)
{
  MapMatd tempX(as<MapMatd> (x));
  MatrixXd res = tempX.transpose() * SigmaMatrix;
  return wrap(res.row(0));
}

/*
 * return the minimum value of each row
 */
// [[Rcpp::export]]
NumericVector rowMin(NumericMatrix& x)
{
  NumericVector res(x.nrow());
  for(int i = 0; i < x.ncol(); i++)
  {
    res[i] = min(x(i, _));
  }
  return res;
}

/*
 * return the element-wise minimum values between two vectors
 */
// [[Rcpp::export]]
NumericVector compareMin(NumericVector& x, NumericVector& y)
{
  NumericVector res(x.size());
  for(int i = 0; i < x.size(); i++)
  {
    res[i] = minTwo(x[i], y[i]);
  }
  return res;
}


/*
 * bind a vector to matrix by column.
 * an alternative function is sugar::cbind
 */
NumericMatrix cbindV(NumericMatrix& x, NumericVector& y)
{
  NumericMatrix res(x.nrow(), x.ncol() + 1);
  for(int i = 0; i < x.ncol(); i++)
  {
    res(_, i) = x(_, i);
  }
  res(_, x.ncol()) = y;
  return res;
}

NumericMatrix rbindV(NumericMatrix& x, NumericVector& y)
{
  NumericMatrix res(x.nrow() + 1, x.ncol());
  for(int i = 0; i < x.nrow(); i++)
  {
    res(i, _) = x(i, _);
  }
  res(x.nrow(), _) = y;
  return res;
}

/*
 * compute the pair-wise distance between two observations.
 */
NumericMatrix fastpdist(NumericMatrix& X, NumericMatrix& Y)
{
  int n(X.nrow()), m(Y.nrow()), k(Y.ncol());
  NumericMatrix tempX = NumericMatrix(n, k, (X * X).begin());
  NumericMatrix tempY = NumericMatrix(m, k, (Y * Y).begin());
  NumericVector Xn = rowSums(tempX);
  NumericVector Yn = rowSums(tempY);
  MatrixXd C = as<MapMatd> (X) * as<MapMatd> (Y).transpose();
  NumericMatrix Cn = wrap(-2.0 * C);
  for(int i = 0; i < m; i++)
  {
    Cn(_, i) = Cn(_, i) + Xn;
  }
  for(int i = 0; i < n; i++)
  {
    Cn(i, _) = Cn(i, _) + Yn;
  }
  return NumericMatrix(n, m, sqrt(Cn).begin());
}


/*
 * erase a row from one matrix
 */
NumericMatrix rowErase(NumericMatrix& x, int rowID)
{
  NumericMatrix res(x.nrow() - 1, x.ncol());
  int iter(0);
  for(int i = 0; i < x.nrow(); i++)
  {
    if(i != rowID)
    {
      res.row(iter) = x.row(i);
      iter++;
    }
  }
  return res;
}


/*
 * erase a column from one matrix
 */
NumericMatrix colErase(NumericMatrix& x, int colID)
{
  NumericMatrix res(x.nrow(), x.ncol() - 1);
  int iter(0);
  for(int i = 0; i < x.ncol(); i++)
  {
    if(i != colID)
    {
      res.column(iter) = x.column(i);
      iter++;
    }
  }
  return res;
}


// [[Rcpp::export]]
bool isPrime(int n)
{
  if(n <= 1) return false;
  if(n <= 3) return true;
  if (n%2 == 0 || n%3 == 0) return false;
  for (int i=5; i*i<=n; i=i+6)
    if (n%i == 0 || n%(i+2) == 0)
      return false;
    return true;
}

int reportMaxPrime(int n)
{
  int MaxPrime(0);
  for(int i = 0; i <= n; i++)
  {
    if(isPrime(i))
    {
      MaxPrime = MaxTwo(MaxPrime, i);
    }
  }
  return MaxPrime;
}


NumericMatrix varCPP(NumericMatrix& x)
{
  NumericVector Means = colMeans(x);
  NumericMatrix tempXx(clone(x));
  for(int i = 0; i < x.nrow(); i++)
  {
    tempXx(i, _) = tempXx(i, _) - Means;
  }
  MapMatd tempX(as<MapMatd> (tempXx));
  return wrap((tempX.transpose() * tempX) / (x.nrow() - 1.0));
}


NumericMatrix sqrtVarMatrix(NumericMatrix& x)
{
  MapMatd tempX(as<MapMatd> (x));
  SelfAdjointEigenSolver<MatrixXd> res(tempX);
  // VectorXd eigenvalues = res.eigenvalues();
  // MatrixXd eigenvectors = res.eigenvectors();
  VectorXd eigenvalues = res.eigenvalues();
  for(int i = 0; i < eigenvalues.size(); i++)
  {
    eigenvalues[i] = 1.0 / sqrt(eigenvalues[i]);
  }
  return wrap(res.eigenvectors() * eigenvalues.asDiagonal()  * res.eigenvectors().transpose());
}


NumericVector quantileCPP(NumericVector& x, NumericVector& q)
{
  // NumericVector y(clone(x));
  // std::sort(y.begin(), y.end());
  // return y[floor(x.size() * (q))];
  Environment stats("package:stats");
  Function quantile = stats["quantile"];
  int npr = q.size();
  NumericVector ans(npr);
  ans = quantile(x, q);
  return ans;
}

NumericVector sampleCPP(int dim)
{
  Environment stats("package:base");
  Function sample = stats["sample"];
  return sample(dim);
}

NumericMatrix subMatrixCols(NumericMatrix& x, NumericVector& cols)
{
  NumericMatrix res(x.nrow(), cols.size());
  for(int i = 0; i < cols.size(); i++)
  {
    res(_, i) = x(_, cols[i]);
  }
  return res;
}

IntegerVector orderCPP(NumericVector x) {
  if (is_true(any(duplicated(x)))) {
    Rf_warning("There are duplicates in 'x'; order not guaranteed to match that of R's base::order");
  }
  NumericVector sorted = clone(x).sort();
  return match(sorted, x);
}


#endif
