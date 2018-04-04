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
List chooseMED(NumericMatrix Cand, // candidate points
                        NumericVector lfCand, // logarithms of candidate points
                        int n, // the number of MED points
                        NumericMatrix SigmaK, // Covariance matrix of candidate points
                        double gamma, //
                        double s)
{
  int N(Cand.nrow());// number of candidate list
  int MaxLfIndex(0); // the index of maximum
  MatrixXd m_MED(n, Cand.ncol()); // MED points
  NumericVector lfVector(n);
  MapMatd m_SigmaK(as<MapMatd> (SigmaK));
  MatrixXd SigmaInv = m_SigmaK.inverse();
  NumericVector current_point;
  NumericVector currentCriterion(N, 0.0);
  NumericVector minCriterion(N, 1e16);
  double current_lf(0.0);

  // find out the first MED point which with maximum logarithm function value
  MaxLfIndex = argMax(lfCand);
  current_point = Cand(MaxLfIndex, _);
  current_lf = lfCand[MaxLfIndex];
  m_MED.row(0) = as<MapVectd> (current_point);
  lfVector[0] = current_lf;

  // search sequential MED points using criterion
  NumericMatrix transCand = transMatrix(SigmaInv, Cand);
  NumericVector transCurrentPoint = transVector(SigmaInv, current_point);
  NumericMatrix MEDCriterion(N, n); // define the criterion matrix
  std::fill(MEDCriterion.begin(), MEDCriterion.end(), 1e16); // fill with large enough number
  for(int i = 1; i < n; i++)
  {
    currentCriterion = 0.5 * gamma * (lfCand + current_lf) + Cand.ncol() *
      getLogDistVector(transCand, transCurrentPoint, s);
    minCriterion = compareMin(minCriterion, currentCriterion);
    MaxLfIndex = argMax(minCriterion);
    current_point = Cand(MaxLfIndex, _);
    current_lf = lfCand[MaxLfIndex];
    lfVector[i] = current_lf;
    m_MED.row(i) = as<MapVectd> (current_point);
    transCurrentPoint = transVector(SigmaInv, current_point);
  }
  return List::create(_["Design"] = wrap(m_MED), _["lfVector"] = lfVector);
}


/*
 * Augment the candidate points.
 */
// [[Rcpp::expoer]]
NumericMatrix augCandidate(NumericMatrix& AugMatrix, NumericMatrix& ExistMatrix)
{
  int N1(AugMatrix.nrow());
  int N2(ExistMatrix.nrow());
  NumericMatrix DistMat = fastpdist(AugMatrix, ExistMatrix);
  NumericMatrix DistAug = fastpdist(AugMatrix, AugMatrix);
  NumericMatrix res(round(N1 / 2.0), AugMatrix.ncol());
  int MaxIndex(0);
  for(int i = 0; i < round(N1 / 2.0); i++)
  {
    MaxIndex = which_max(rowMin(DistMat));
    res(i, _) = AugMatrix(MaxIndex, _);
    AugMatrix = rowErase(AugMatrix, MaxIndex);
    DistMat = rowErase(DistMat, MaxIndex);
    NumericVector DistV = rowErase(DistAug, MaxIndex)(_, MaxIndex);
    DistAug = rowErase(DistAug, MaxIndex);
    DistAug = colErase(DistAug, MaxIndex);
    DistMat = cbind(DistMat, DistV);
  }
  return res;
}


/***R
Cand <- matrix(rnorm(500), ncol = 1)
lfCand <- log(dnorm(Cand))
SigmaK <- matrix(c(1))
y <- chooseMED(Cand, lfCand, 50, SigmaK, 1.0, 2.0)
hist(y$Design)
*/
