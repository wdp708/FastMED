#include "FastLattice.hpp"

// [[Rcpp::export]]
NumericMatrix test_generateLattice(int n, int p)
{
  return generateLattice(n, p);
}

/***R
test_generateLattice(11, 2)
*/
