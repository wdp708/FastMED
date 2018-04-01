/*
 * generate MED points.
 */
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

#include "BaseFunc.hpp"

/*
 * Choose $n$ MED points from candidates
 */
// [[Rcpp::export]]
NumericMatrix chooseMED(NumericMatrix Cand, // candidate points
                        NumericVector lfCand, // logarithms of candidate points
                        int n, // the number of MED points
                        NumericMatrix SigmaK, // Covariance matrix of candidate points
                        double gamma, //
                        double s)
{
  int N(Cand.nrow());
  int MaxLfIndex(0);
  MatrixXd m_MED(n, Cand.ncol());
  MaxLfIndex = argMax(lfCand);
  MapMatd m_SigmaK(as<MapMatd> (SigmaK));
  MatrixXd SigmaInv = m_SigmaK.inverse();

  return NumericMatrix(wrap(m_MED));
}
