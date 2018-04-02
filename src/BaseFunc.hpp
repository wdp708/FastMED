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

// [[Rcpp::export]]
int argMin(NumericVector& x)
{
  return std::min_element(x.begin(), x.end()) - x.begin();
}

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
      res = 1.0  / n * accumulate(temp.begin(), temp.end(), 1.0, multiplies<double>());
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
