/*
 * Generate MED points by using fast algorithm which will be published.
 * Copyright: Dianpeng Wang <wdp@bit.edu.cn>
 * Date: 20180524
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
 * This is the main function which is used to generate MED from objective density.
 * @ dim: the dimension of distribution
 * @ ExploreDesign: one given space-filling design which will be used to search new candidate in a neighhorhood
 * @ logFunc: Logarithm of density function
 */
List generateMED(int dim, NumericMatrix ExploreDesign, Function logFunc)
{
  int DesignSize = report
}
