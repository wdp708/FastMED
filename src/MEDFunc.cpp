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
  int N(Cand.nrow());// number of candidate list
  int MaxLfIndex(0); // the index of maximum
  MatrixXd m_MED(n, Cand.ncol()); // MED points
  MapMatd m_SigmaK(as<MapMatd> (SigmaK));
  MatrixXd SigmaInv = m_SigmaK.inverse();
  NumericVector current_point;

  // find out the first MED point which with maximum logarithm function value
  MaxLfIndex = argMax(lfCand);
  current_point = Cand(MaxlfIndex, _);
  m_MED.row(0) = as<MapVectd> (current_point);

  return NumericMatrix(wrap(m_MED));
}
