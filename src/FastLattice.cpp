/*
 * Generate scalar rank-1 lattice rules using fast algorithm.
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

// [[Rcpp::export]]
ComplexVector FFT(NumericVector& x){
  Environment base("package:stats");
  Function fft = base["fft"];
  ComplexVector res = fft(_["z"] = x);
  return res;
}

// [[Rcpp::export]]
ComplexVector IFFT(ComplexVector& x){
  Environment base("package:stats");
  Function fft = base["fft"];
  ComplexVector res = fft(_["z"] = x, _["inverse"] = true);
  return res;
}

double calOmega(double x){
  return 2.0 * PI * PI * (x * x - x + 1.0 / 6.0);
}

NumericVector calOmega(NumericVector& x){
  NumericVector res(x.size());
  for(int i = 0; i < x.size(); i++)
  {
    res[i] = calOmega(x[i]);
  }
  return res;
}

// Calculates x ^ a (mod n) using the well known “Russian peasant method”
// [[Rcpp::export]]
int POWMOD(int x, int a, int n)
{
  int y(1);
  int u(x);
  while(a > 0){
    if(a % 2 == 1){
      y = (y * u) % n;
    }
    u = (u * u) % n;
    a = floor(a / 2.0);
  }
  return y;
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
int findPrimefactors(unordered_set<int> &s, int n)
{
  while(n%2 == 0)
  {
    s.insert(2);
    n = n / 2;
  }

  for(int i = 3; i < sqrt(n); i = i + 2)
  {
    while(n % i == 0)
    {
      s.insert(i);
      n = n / i;
    }
  }
  if(n > 2)
    s.insert(n);
  return 1;
}

// [[Rcpp::export]]
// store the prime factors of number n.
NumericVector findFactorize(int n)
{
  unordered_set<int> s;
  findPrimefactors(s, n);
  NumericVector res(s.size());
  int i(0);
  for(auto it = s.begin(); it != s.end(); it++)
  {
    res[i] = (*it);
    i++;
  }
  return res;
}

// [[Rcpp::export]]
int generateOrp(int n)
{
  try{
    if(!isPrime(n))
    {
      throw "n is not a prime";
    }
    NumericVector primef = findFactorize(n - 1);
    int g(2), i(1);
    while(i <= primef.size())
    {
      if(POWMOD(g, (n - 1) / primef[i - 1], n) == 1)
      {
        g++;
        i = 0;
      }
      i++;
    }
    return g;
  }catch(const char* & e){
    cout<<e<<endl;
    return -1;
  }
}


// [[Rcpp::export]]
NumericMatrix generateLattice(int n, int p){
  int s_max(p); // number of dimensions
  NumericVector gamma(s_max, 1.0 / s_max);
  NumericVector beta(s_max, 1.0);
  NumericMatrix res(n, p);
  std::fill(res.begin(), res.end(), 0.0);
  int m((n - 1) / 2);
  NumericVector E2(m, 0.0);
  NumericVector cumbeta = cumprodC(beta);
  int g(0);
  try
  {
    g = generateOrp(n);
    if(g == -1){
      throw "n must be prime!";
    }
  }catch(const char* & e)
  {
    cout<<e<<endl;
    return res;
  }
  NumericVector perm(m, 0.0);
  perm[0] = 1.0;
  for(int j = 1; j < perm.size(); j++)
  {
    perm[j] = ((int)perm[j - 1] * g) % n;
    perm[j] = minTwo(n - perm[j], perm[j]);
  }
  NumericVector temp = perm / n;
  NumericVector psi = calOmega(temp);
  double psi0 = calOmega(0.0);
  ComplexVector fft_psi = FFT(psi);
  NumericVector z(s_max, 0.0);
  NumericVector e2(s_max, 0.0);
  NumericVector q(m, 1.0);
  double q0 = 1;
  //
  Environment base("package:base");
  Function Re = base["Re"];
  int min_index;
  double min_E2;
  //double noise;
  NumericVector temp_psi(m, 1.0);
  int temp_index(0);
  for(int s = 1; s <= s_max; s++)
  {
    ComplexVector fft_q = FFT(q);
    fft_q = fft_q * fft_psi;
    fft_q = IFFT(fft_q);
    E2 = Re(_["z"] = fft_q);
    E2 = E2 / m;
    min_index = argMin(E2);
    min_E2 = E2[min_index];
    if(s==1)
    {
      min_index = 0;
      //noise = abs(E2[0] - min_E2);
    }
    z[s-1] = perm[min_index];
    e2[s-1] = -cumbeta[s-1] +(beta[s-1] * (q0 + 2 * sum(q)) + gamma[s-1]
                            * (psi0 * q0 + 2 * min_E2)) / n;
    temp_index = min_index;
    for(int i = 0; i <= min_index; i++)
    {
      temp_psi[i] = psi[min_index - i];
    }
    temp_index = m - 1;
    for(int i = min_index + 1; i < psi.size(); i++)
    {
      temp_psi[i] = psi[temp_index];
      temp_index--;
    }
    q = (beta[s] + gamma[s] * temp_psi) * q;
    q0 = (beta[s] + gamma[s] * psi0) * q0;
  }
  MatrixXd m_res(n, p);
  VectorXd l_bound(p);
  VectorXd u_bound(p);
  for(int i = 0; i < p; i++)
  {
    l_bound(i) = 1e60;
    u_bound(i) = -1e60;
    for(int j = 0; j < n; j++)
    {
      m_res(j, i) = ((int) z[i] * (j + 1)) % n;
      m_res(j, i) = m_res(j, i) * 1.0 / n;
      l_bound(i) = minTwo(l_bound(i), m_res(j, i));
      u_bound(i) = MaxTwo(u_bound(i), m_res(j, i));
    }
  }
  m_res = 0.5 / n + (1 - 1.0 / n) * (m_res - MatrixXd::Ones(n, 1) * l_bound.transpose()).array() /
    (MatrixXd::Ones(n, 1) * (u_bound - l_bound).transpose()).array();
  return NumericMatrix(wrap(m_res));
  // return res;
}



/***R

# x <- 1:4
# fft(x)
# Re(FFT(x))
# IFFT(FFT(x)) / length(x)
# POWMOD(2, 3, 10)
# isPrime(10)
# findFactorize(511)
# generateOrp(11)
# cumprodC(c(1.0, 2.0, 1.0))
# argMin(c(2.0, 1.5, 2.5, 3.5))
generateLattice(11, 2)
*/
