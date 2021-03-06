#pragma once

#include <limits>
#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>
#include <iostream>

  template<int n, typename RF> 
constexpr RF powr(const RF& x)
{
  if constexpr (n == 0)
    return 1.;
  else if constexpr (n < 0)
    return 1. / powr<-n>(x);
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

void debug_out(const std::string& msg, bool debug)
{
  if (debug)
    std::cout << msg << std::flush;
}

int locate(const std::vector<double> &xx, const double& x)
{
  /*
   * This function has been taken and slightly been modified from the
   * "Numerical Recipes in C" book.
   *
   * Given a vector xx of size n and given a value x, it returns an integer j,
   * such that x is between xx[j] and xx[j+1]. xx must be sorted. Returns -1
   * or n when is x is out of range.
   */

  const unsigned n = xx.size();

  if (isEqual(x, xx[0]))
    return 0;
  else if (isEqual(x, xx[n-1]))
    return n-2;

  int ju, jm;
  int jl;
  bool ascnd;

  jl = -1;
  ju = n;

  ascnd = (xx[n-1] >= xx[0]);

  while (ju - jl > 1) {
    jm = (ju + jl) / 2;
    if ((x >= xx[jm]) == ascnd) {
      jl = jm;
    } else {
      ju = jm;
    }
  }

  return jl;
}
