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
 * Choose $n$ MED points from candidates.
 */
// [[Rcpp::export]]
List chooseMED(NumericMatrix Cand, // candidate points
                        NumericVector lfCand, // logarithms of candidate points
                        int n, // the number of MED points
                        NumericMatrix SigmaK, // Covariance matrix of candidate points
                        double gamma, //
                        double s)
{
  int N(Cand.nrow());// number of candidates
  int MaxLfIndex(0); // the index of maximum log density function values
  MatrixXd m_MED(n, Cand.ncol()); // MED points
  NumericVector lfVector(n); // log density function values at corresponding MED points
  MatrixXd SigmaInv = as<MapMatd> (sqrtVarMatrix(SigmaK)); // the sqrt inverse of variance-covariance matrix
  NumericVector current_point; // new added MED point
  NumericVector currentCriterion(N, 0.0); // MED criterion with current point and candidates
  NumericVector minCriterion(N, 1e16); // minimum MED criterion along row direction
  double current_lf(0.0);

  /*
   * find out the first MED point which with maximum logarithm function value
   */
  MaxLfIndex = argMax(lfCand); // index of point with maximum log density function value
  current_point = Cand(MaxLfIndex, _); // point with maximum log density function value
  current_lf = lfCand[MaxLfIndex]; // maximum log density function value
  m_MED.row(0) = as<MapVectd> (current_point); // set the point with maximum log density function value as first MED point
  lfVector[0] = current_lf; // store the corresponding log density function value

  /*
   * search the next MED point
   */
  NumericMatrix transCand = transMatrix(SigmaInv, Cand); // translate andidates by multipling sqrt inverse of variance matrix
  NumericVector transCurrentPoint = transVector(SigmaInv, current_point); // translate current point
  NumericMatrix MEDCriterion(N, n); // criterion matrix
  std::fill(MEDCriterion.begin(), MEDCriterion.end(), 1e16); // fill with large enough number
  for(int i = 1; i < n; i++)
  {
    currentCriterion = 0.5 * gamma * (lfCand + current_lf) + Cand.ncol() *
      getLogDistVector(transCand, transCurrentPoint, s); // compute MED criterion with current point and candidates
    minCriterion = compareMin(minCriterion, currentCriterion); // comapre the MED criterion
    MaxLfIndex = argMax(minCriterion); // index of best point within candidates
    current_point = Cand(MaxLfIndex, _); // choose the best point as the next MED point
    current_lf = lfCand[MaxLfIndex]; // record the corresponding log density function value
    lfVector[i] = current_lf;
    m_MED.row(i) = as<MapVectd> (current_point); // extend MED matrix
    transCurrentPoint = transVector(SigmaInv, current_point); // update current point
  }
  return List::create(_["Design"] = wrap(m_MED), _["lfVector"] = lfVector);
}


/*
 * Augment the candidate points.
 * In this step, for each candiate one new point is selected in its neighborhood which is defined by the nearest
 */
// [[Rcpp::export]]
NumericMatrix augCandidate(NumericMatrix& AugMatrix, NumericMatrix& ExistMatrix)
{
  int N1(AugMatrix.nrow());
  //int N2(ExistMatrix.nrow());
  NumericMatrix DistMat = fastpdist(AugMatrix, ExistMatrix);
  NumericMatrix DistAug = fastpdist(AugMatrix, AugMatrix);
  NumericMatrix res(round(N1 / 2.0), AugMatrix.ncol());
  int MaxIndex(0);
  for(int i = 0; i < round(N1 / 2.0); i++)
  {
    //cout<<rowMin(DistMat)<<endl;
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

class LKPredictor
{
public:
  LKPredictor(){};
  LKPredictor(NumericMatrix& DesignPoints, NumericVector& LfVector)
  {
    m_DesignPoints = clone(DesignPoints);
    m_LfVector = clone(LfVector);
    this->N = DesignPoints.nrow();
    this->Dim = DesignPoints.ncol();
  }
  void setData(NumericMatrix& DesignPoints, NumericVector& LfVector)
  {
    m_DesignPoints = clone(DesignPoints);
    m_LfVector = clone(LfVector);
    this->N = DesignPoints.nrow();
    this->Dim = DesignPoints.ncol();
  }
  void fit()
  {
    NumericMatrix res = fastpdist(this->m_DesignPoints, this->m_DesignPoints);
    double quant = MaxTwo(1.0 - 5.0 * this->Dim / this->N, .5);
    // fill diagnal
    res.fill_diag(10.0 * this->Dim);
    NumericVector tempVector = rowMin(res);
    NumericVector quantVect = NumericVector::create(quant);
    this->dbar = 2.0 * this->Dim * quant * quantileCPP(tempVector, quantVect)[0];
    res.fill_diag(0.0);
    int Xsize = res.ncol() * res.nrow();
    for(int i = 0; i < Xsize; i++)
    {
      res[i] = 1.0 / (1.0 + (res[i] / this->dbar) * (res[i] / this->dbar));
    }
    res.fill_diag(1.0 + 1e-6);
    cout<<dbar<<endl;
    this->corMatrixInv = as<MapMatd> (res).inverse();
    this->coef = this->corMatrixInv * as<MapVectd> (this->m_LfVector);
    this->dnom = this->corMatrixInv.rowwise().sum();
  }
  List predictor(NumericVector& x)
  {
    NumericVector distVect = fastpdist2(this->m_DesignPoints, x);
    for(int i = 0; i < distVect.size(); i++)
    {
      distVect[i] = 1.0 / (1.0 + (distVect[i] / this->dbar) * (distVect[i] / this->dbar));
    }
    MapVectd temp1(as<MapVectd> (distVect));
    double mind = min(distVect);
    double pred = (temp1.transpose() * this->coef)[0] / (temp1.transpose() * this->dnom)[0];
    return List::create(_["mind"]=mind, _["pred"]=pred);
  }
public:
  NumericMatrix m_DesignPoints;
  NumericVector m_LfVector;
private:
  MatrixXd corMatrixInv;
  VectorXd coef;
  VectorXd dnom;
  double dbar;
  int N;
  int Dim;
};

// [[Rcpp::export]]
List LKpred(NumericVector& x, NumericMatrix& DesignPoints, NumericVector& LfVector)
{
  LKPredictor m_pred(DesignPoints, LfVector);
  m_pred.fit();
  return m_pred.predictor(x);
}

// [[Rcpp::export]]
List generateMEDPoints(int dim, NumericMatrix ExploreDesign, Function logFunc)
{
  /*
   * In the initial step, generate samples from uniform distribution by latttice rules
   */
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
  int ncl = ExploreDesign.nrow();
  int nc = round(ncl / 2.0);
  NumericMatrix AugMatrix;
  NumericMatrix M1, M2, M0; // temp matrix
  LKPredictor m_pred;
  /*
   *
   */
  NumericVector q = NumericVector::create(0.9, 0.1);
  NumericVector q_vector(q.size());
  NumericVector lfD(m_MED[1]); // vector of log density function value at MED points which are chosen in last step
  NumericMatrix MEDesign = m_MED[0]; // MED points which are chosen in last step
  NumericVector perIndex(dim);
  NumericMatrix pairwiseDist(DesignSize, DesignSize);
  NumericVector DistRowMin(DesignSize);
  NumericVector RadiusVect(DesignSize);
  /*
   * Definitions of variables
   */
  NumericVector InitialDistVect(InitialDesign.nrow());
  NumericVector tempCol(InitialDesign.ncol());
  NumericVector tempOrder;
  NumericVector OrderCl(DesignSize);
  NumericVector OrderClD(nc+1);
  NumericVector OrderClD2(round(nc / 2.0));
  NumericVector dD;
  NumericVector ru = NumericVector::create(0.5, 0.5);
  /*
   * define the design in K^{th} step and the corresponding log density function value
   */
  NumericMatrix MEDk(DesignSize, MEDesign.ncol());
  std::fill(MEDk.begin(), MEDk.end(), 0.0);
  NumericVector lfDk(DesignSize, 0.0);
  int MaxIndex(0);
  NumericMatrix newPoints(DesignSize, MEDesign.ncol());
  std::fill(newPoints.begin(), newPoints.end(), 0.0);
  NumericVector newlf(DesignSize, 0.0);
  NumericVector tempPoint;

  for(int loop_k = 2; loop_k <= 2; loop_k++)
  {
    gamma = (loop_k - 1.0) / (K - 1.0);
    SigmaK = (loop_k - 1.0) / loop_k * SigmaK;
    // TODO: need to update SigmaK
    q_vector = quantileCPP(lfD, q);
    s = round(2.0 * (1.0 - exp(-gamma * (q_vector[0] - q_vector[1])))); // the s value in distance measure
    perIndex = sampleCPP(dim) - 1;
    ExploreDesign = subMatrixCols(ExploreDesign, perIndex); // permutation of matrix columns
    pairwiseDist = fastpdist(MEDesign, MEDesign); // compute the pair-wise log distances between MED points
    pairwiseDist.fill_diag(10 * dim); // fill the diagonal elments of distance matrix with 10*p
    DistRowMin = rowMin(pairwiseDist); // find out the nearest distance within MED points for each point
    /*
     * compute the search radius for each MED point
     */
    for(int i = 0; i < DesignSize; i++)
    {
      RadiusVect[i] = pairwiseDist(i, orderCPP(pairwiseDist(i, _))[1] - 1);
    }
    /*
     * For each MED point, search one new candidate within its neighborhood
     */
    for(int j = 0; j < 1; j++)
    {
      tempCol = MEDesign(j, _); // center of neighborhood
      InitialDistVect = fastpdist2(InitialDesign, tempCol); // compute the distances between center and candidates
      // find out the first $n$ indexes according to distance measures
      tempOrder = orderCPP(InitialDistVect);
      OrderCl = head(tempOrder, DesignSize);
      dD = pairwiseDist(j, _); // distances between center and MED points
      dD[j] = 0.0; // re-fill the diagonal elements

      /*
       * rescale the space-filling design
       */
      AugMatrix = RadiusVect[j] / sqrt(dim) * 2.0 * (ExploreDesign - 0.5);
      for(int l = 0; l < AugMatrix.nrow(); l++)
      {
        AugMatrix(l, _) = tempCol + AugMatrix(l, _);
      }

      /*
       * linear combinations of the adjacent points
       */
      tempOrder = orderCPP(dD);
      OrderClD = head(tempOrder, nc + 1);
      if(j == 1)
      {
        //TODO: re-set the random weights
        ru[1] = 0.45;//runif(1, 0.25, 0.75)[0];
      }
      OrderClD2 = tail(head(OrderClD, round(nc / 2.0) + 1), round(nc / 2.0));
      NumericVector tempRowID = tail(OrderClD, OrderClD.size() - 1) - 1;
      M1 = ru[0] * subMatrixRows(MEDesign,  tempRowID);
      for(int l = 0; l < M1.nrow(); l++)
      {
        M1(l, _) = M1(l, _) + (1 - ru[0]) * MEDesign(j, _);
      }
      tempRowID = OrderClD2 - 1;
      M2 = subMatrixRows(MEDesign, tempRowID);
      for(int l = 0; l < M2.nrow(); l++)
      {
        M2(l, _) = (1.0 + ru[1]) * MEDesign(j, _) - ru[1] * M2(l, _);
      }
      M0 = rbindM(M1, M2);


      // find out the first $n$ nearest points in candiates to current center
      tempRowID = OrderCl - 1;
      NumericMatrix tempInitialDes = subMatrixRows(InitialDesign, tempRowID);
      NumericVector tempInitialLfVect = subElements(LfVector, tempRowID);
      NumericVector lBound = colMin(tempInitialDes);
      NumericVector uBound = colMax(tempInitialDes);
      lBound = lBound - (gamma - 1.0 / (K - 1)) * (uBound - lBound) / 6;
      uBound = uBound + (gamma - 1.0 / (K - 1)) * (uBound - lBound) / 6;
      // trunctate rescaled space-filling design points by bounds which are defined using nearest neighborhood points
      for(int l = 0; l < dim; l++)
      {
        AugMatrix(_, l) = pmin(pmax(AugMatrix(_, l), lBound[l]), uBound[l]);
      }
      AugMatrix = rbindM(AugMatrix, M0);
      AugMatrix = augCandidate(AugMatrix, tempInitialDes); // filter new points by distance
      /*
       * fit surrogate model using limit kriging model
       */
      m_pred.setData(tempInitialDes, tempInitialLfVect);
      m_pred.fit();

      NumericVector predVector(AugMatrix.nrow()); // prediction of log density function values at new points
      for(int l = 0; l < AugMatrix.nrow(); l++)
      {
        NumericVector tempAugRow = AugMatrix(l, _);
        List predRes = m_pred.predictor(tempAugRow);
        double mind = predRes[0];
        double predVal = predRes[1];
        double penalty(0.0);
        if(mind < .1 / DesignSize)
        {
          penalty = log(mind);
        }
        else
        {
          penalty = 0.0;
        }
        if(j == 1)
        {
          predVector[l] = gamma * predVal + penalty;
        }
        else
        {
          NumericVector tempCriterion(j - 1, 1e16);
          NumericVector searchPoint = AugMatrix(l, _);
          NumericVector tempMEDPoint;
          for(int lll = 0; lll < j - 1; lll++)
          {
            tempMEDPoint = MEDk(lll, _);
            tempCriterion[lll] = .5 * gamma * predVal + .5 * gamma * lfDk[lll]
            + getLogDist(tempMEDPoint, searchPoint, s) + penalty;
          }
          predVector[l] = min(tempCriterion);
        }
      }
      MaxIndex = argMax(predVector);
      tempPoint = AugMatrix(MaxIndex, _);
      newPoints(j, _) = tempPoint;
      FuncRes = logFunc(tempPoint);
      newlf[j] = FuncRes[0];
      /*
       * update MED points in augmentation of k^{th} step
       */
      // if(j == 1)
      // {
      //   if
      // }

    }
    InitialDesign = rbindM(InitialDesign, newPoints);
    LfVector = extendV(LfVector, newlf);
    m_MED = chooseMED(InitialDesign, LfVector, DesignSize, SigmaK, gamma, s);
    lfD = m_MED[1];
    MEDesign = wrap(m_MED[0]);
    SigmaK = varCPP(MEDesign);
  }
  return List::create(_["Design"]=MEDesign, _["LfuncVec"]=lfD, _["PairWiseDist"]=pairwiseDist,
                      _["Radius"]=RadiusVect, _["InitialD"]=InitialDesign, _["debug"]=AugMatrix);
  //return List::create(_["Design"]=InitialDesign, _["LfuncVec"]=LfVector);
}

/***R
x <- matrix(c(1, 2, 3, 4), nrow = 2)
m_result = generateMEDPoints(2, D1, lf)
*/
