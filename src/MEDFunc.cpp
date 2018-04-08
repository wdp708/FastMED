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
#include "FastLattice.hpp"

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
  // MapMatd m_SigmaK(as<MapMatd> (SigmaK));
  MatrixXd SigmaInv = as<MapMatd> (sqrtVarMatrix(SigmaK));
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
// [[Rcpp::export]]
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

// [[Rcpp::export]]
List generateMEDPoints(int dim, NumericMatrix ExploreDesign, Function logFunc)
{
  int DesignSize = reportMaxPrime(100 + 5 * dim); // the number of MED points
  NumericMatrix LatticDesign = generateLattice(DesignSize, dim);
  NumericMatrix InitialDesign(clone(LatticDesign));
  NumericVector LfVector(DesignSize, 0.0);
  NumericVector FuncRes;
  for(int i = 0; i < DesignSize; i++)
  {
    FuncRes = logFunc(InitialDesign(i, _));
    LfVector[i] = FuncRes[0];
  }
  int K = ceil(4 * sqrt(dim)); // the number is SA
  double gamma = 1.0 / K;
  NumericMatrix SigmaK = varCPP(InitialDesign);
  double s(0.0);
  List m_MED = chooseMED(InitialDesign, LfVector, DesignSize, SigmaK, gamma, s);

  /*
   *
   */
  NumericVector q = NumericVector::create(0.9, 0.1);
  NumericVector q_vector(q.size());
  NumericVector lfD(m_MED[1]);
  NumericMatrix MEDesign = m_MED[0];
  NumericVector perIndex(dim);
  NumericMatrix pairwiseDist(DesignSize, DesignSize);
  NumericVector DistRowMin(DesignSize);
  NumericVector RadiusVect(DesignSize);
  for(int loop_k = 2; loop_k <= K; loop_k++)
  {
    gamma = (loop_k - 1.0) / (K - 1.0);
    SigmaK = (loop_k - 1.0) / loop_k * SigmaK;
    q_vector = quantileCPP(lfD, q);
    s = round(2.0 * (1.0 - exp(-gamma * (q_vector[0] - q_vector[1]))));
    perIndex = sampleCPP(dim) - 1;
    ExploreDesign = subMatrixCols(ExploreDesign, perIndex);
    // compute the pairwise distance between MED points.
    pairwiseDist = fastpdist(MEDesign, MEDesign);
    pairwiseDist.fill_diag(10 * dim);
    DistRowMin = rowMin(pairwiseDist);
    for(int i = 0; i < DesignSize; i++)
    {
      RadiusVect[i] = pairwiseDist(i, orderCPP(pairwiseDist(i, _))[1] - 1);
    }
    // loop for adding MED
    for(int j = 0; j < DesignSize; j++)
    {
      //TODO:
    }
  }


  return List::create(_["Design"]=MEDesign, _["LfuncVec"]=lfD, _["PairWiseDist"]=pairwiseDist,
                      _["Radius"]=RadiusVect);
  //return List::create(_["Design"]=InitialDesign, _["LfuncVec"]=LfVector);
}

/***R
x <- matrix(c(1, 2, 3, 4), nrow = 2)
m_result = generateMEDPoints(2, x, lf)
*/
