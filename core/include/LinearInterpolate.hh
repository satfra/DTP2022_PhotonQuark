#pragma once

#include <vector>

#include "Utils.hh"

template<typename _RF, typename _RF_f>
class lInterpolator
{
  public:
    using RF = _RF;
    using RF_f = _RF_f;
    using Range = std::vector<RF_f>;
    using Grid = std::vector<RF>;

    lInterpolator(const Range& _x, const Grid& _f)
      : x(_x), f(_f), a(x.front()), b(x.back()) {}

    RF_f operator()(const RF& y) const
    {
      if(y > b || y < a)
        throw std::runtime_error("Interpolating outside bounds");

      const auto idx = locate(x, y);

      const RF t = (y - x[idx]) / (x[idx+1] - x[idx]);
      return t*f[idx+1] + (1.-t)*f[idx];
    }

  private:
    const Range& x;
    const Grid& f;
    const RF& a,b;
};

template<typename _RF, typename _RF_f>
class lInterpolator2d
{
  public:
    using RF = _RF;
    using RF_f = _RF_f;
    using Range = std::vector<RF>;
    using Grid = std::vector<std::vector<RF_f>>;

    lInterpolator2d(const Range& _x1, const Range& _x2, const Grid& _f)
      : x1(_x1),x2(_x2), f(_f), a(x1.front()), b(x1.back()), c(x2.front()), d(x2.back()) {}

    RF_f operator()(const RF& y,const RF& z) const
    {
      if (y > b || y < a || z > d || z < c)
        throw std::runtime_error("Interpolating outside bounds");

      const auto idx1 = locate(x1, y);
      const auto idx2 = locate(x2, z);

      const RF t1 = (y - x1[idx1]) / (x1[idx1 + 1] - x1[idx1]);
      const RF t2 = (z - x2[idx2]) / (x2[idx2 + 1] - x2[idx2]);

      return t2 * t1 * f[idx1 + 1][idx2 + 1]
        + t2 * (1. - t1) * f[idx1][idx2 + 1]
        + (1. - t2) * t1 * f[idx1 + 1][idx2]
        + (1. - t2) * (1. - t1)*f[idx1][idx2];
    }

  private:
    const Range& x1, x2;
    const Grid& f;
    const RF& a, b, c, d;
};
