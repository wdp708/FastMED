#include <Rcpp.h>
#include <math.h>
#include <unordered_set>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace RcppEigen;
using namespace Eigen;
using namespace std;
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
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
