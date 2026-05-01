#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
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

    RF_f operator()(const RF& y_in) const
    {
      // Tolerate sub-ULP excursions past the grid endpoints (a few bits
      // of roundoff are unavoidable when the caller computes the query
      // point via floating-point arithmetic, especially under -ffast-math).
      const RF tol = 16. * std::numeric_limits<RF>::epsilon()
                   * std::max<RF>(std::abs(a), std::abs(b));
      if(y_in > b + tol || y_in < a - tol)
        throw std::runtime_error("Interpolating outside bounds");
      const RF y = std::clamp(y_in, a, b);

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

    RF_f operator()(const RF& y_in,const RF& z_in) const
    {
      // Tolerate sub-ULP excursions past the grid endpoints (-ffast-math
      // and identical-formula-but-different-translation-unit can produce
      // bit-different results for the same mathematical value).
      const RF tol_y = 16. * std::numeric_limits<RF>::epsilon()
                     * std::max<RF>(std::abs(a), std::abs(b));
      const RF tol_z = 16. * std::numeric_limits<RF>::epsilon()
                     * std::max<RF>(std::abs(c), std::abs(d));
      if (y_in > b + tol_y || y_in < a - tol_y
       || z_in > d + tol_z || z_in < c - tol_z)
        throw std::runtime_error("Interpolating outside bounds");
      const RF y = std::clamp(y_in, a, b);
      const RF z = std::clamp(z_in, c, d);

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
