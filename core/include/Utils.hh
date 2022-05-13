#pragma once

#include <limits>
#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>

typedef std::vector<std::complex<double>> vec_cmplx;
typedef std::vector<vec_cmplx> mat_cmplx;
typedef std::vector<mat_cmplx> tens_cmplx;
typedef std::vector<tens_cmplx> qtens_cmplx;

typedef std::vector<double> vec_double;
typedef std::vector<vec_double> mat_double;
typedef std::vector<mat_double> tens_double;
typedef std::vector<tens_double> tens2_double;
typedef std::vector<tens2_double> jtens2_double;
typedef std::vector<jtens2_double> ijtens2_double;

  template<int n, typename RF> 
constexpr RF powr(const RF& x)
{
  if constexpr (n == 0)
    return 1.;
  else if constexpr (n < 0)
    return 1.0 / powr<-n>(x);
  else
    return x * powr<n-1>(x);
}

  template<typename T>
bool isEqual(T a, T b, T eps_ = std::numeric_limits<T>::epsilon())
{
  T diff = std::fabs(a - b);
  if (diff <= eps_)
    return true;
  if (diff < std::fmax(std::fabs(a), std::fabs(b)) * eps_)
    return true;
  return false;
}

  template<typename RF>
RF linearMapTo(const RF& val, const RF& A, const RF& B, const RF& a, const RF& b)
{
  return (val - A)*(b-a)/(B-A) + a;
}

  template<typename RF>
std::vector<RF> linearMapTo(const std::vector<RF>& val, const RF& A, const RF& B, const RF& a, const RF& b)
{
  std::vector<RF> ret(val);
  std::for_each(ret.begin(), ret.end(), [&](RF& x){ x = linearMapTo(x, A, B, a, b);});
  return ret;
}
