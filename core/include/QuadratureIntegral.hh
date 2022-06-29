#pragma once

#include <vector>
#include <algorithm>
#include "Utils.hh"
#include "LegendrePolynomials.hh"
#include "parameters.hh"

template<typename POL, typename RF = typename POL::RF>
class yIntegral
{
  private:
    POL polynomials;

  public:
    auto operator()(const std::vector<double> &fun) const
    {
      const auto& w = polynomials.weights();

      auto res = 0.;
      for(unsigned i = 0; i < POL::order; ++i)
        res += w[i] * fun[i];
      return res;
    }
};

// performs a 2d gauss-legendre quadrature integral on a two-dimensional array, which is already discretized on the legendre nodes of an arbitrary range
template<typename POL1, typename POL2, typename RF = typename POL1::RF>
class gaulegIntegral2d
{
  private:
    POL1 polynomials_1;
    POL2 polynomials_2;

  public:
    auto operator()(const std::vector<std::vector<std::complex<double>>> &fun, const RF& a1, const RF& b1, const RF& a2, const RF& b2) const
    {
      const auto& w1 = polynomials_1.weights();
      const RF dx1 = (b1-a1)/2.;

      const auto& w2 = polynomials_2.weights();
      const RF dx2 = (b2-a2)/2.;

      std::complex<double> res = 0.;

      for(unsigned j = 0; j < POL1::order; ++j)
        for(unsigned i = 0; i < POL2::order; ++i)
          res += dx1 * w1[j] * ( dx2 * w2[i] * fun[j][i] );

      return res;
    }
};

template<typename POL, typename RF = typename POL::RF>
class qIntegral
{
  private:
    POL polynomials;

  public:
    template<typename FUN>
    auto operator()(FUN& fun, const RF& a, const RF& b) const
    {
      const auto& z = linearMapTo(polynomials.zeroes(), RF(-1.), RF(1.), a, b);
      const auto& w = polynomials.weights();
      const RF dx = (b-a)/2.;

      decltype(fun(a)) res(0.);
      for(unsigned i = 0; i < POL::order; ++i)
        res += dx * w[i] * fun(z[i]);
      return res;
    }
};

template<typename POL1, typename POL2, typename RF = typename POL1::RF>
class qIntegral2d
{
  private:
    POL1 polynomials_1;
    POL2 polynomials_2;

  public:
    template<typename FUN>
    auto operator()(FUN& fun, const RF& a1, const RF& b1, const RF& a2, const RF& b2) const
    {
      const auto& z1 = linearMapTo(polynomials_1.zeroes(), RF(-1.), RF(1.), a1, b1);
      const auto& w1 = polynomials_1.weights();
      const RF dx1 = (b1-a1)/2.;

      const auto& z2 = linearMapTo(polynomials_2.zeroes(), RF(-1.), RF(1.), a2, b2);
      const auto& w2 = polynomials_2.weights();
      const RF dx2 = (b2-a2)/2.;

      decltype(fun(a1,a2)) result = 0.;

      for(unsigned j = 0; j < POL1::order; ++j)
        for(unsigned i = 0; i < POL2::order; ++i)
          result += dx1 * w1[j] * ( dx2 * w2[i] * fun(z1[j], z2[i]));

      return result;
    }
};

template<typename POL1, typename POL2, typename POL3, typename RF = typename POL1::RF>
class qIntegral3d
{
  private:
    POL1 polynomials_1;
    POL2 polynomials_2;
    POL3 polynomials_3;

  public:
    template<typename FUN>
    auto operator()(FUN& fun, const RF& a1, const RF& b1, const RF& a2, const RF& b2, const RF& a3, const RF& b3) const
    {
      const auto& z1 = linearMapTo(polynomials_1.zeroes(), RF(-1.), RF(1.), a1, b1);
      const auto& w1 = polynomials_1.weights();
      const RF dx1 = (b1-a1)/2.;

      const auto& z2 = linearMapTo(polynomials_2.zeroes(), RF(-1.), RF(1.), a2, b2);
      const auto& w2 = polynomials_2.weights();
      const RF dx2 = (b2-a2)/2.;

      const auto& z3 = linearMapTo(polynomials_3.zeroes(), RF(-1.), RF(1.), a3, b3);
      const auto& w3 = polynomials_3.weights();
      const RF dx3 = (b3-a3)/2.;

      decltype(fun(a1,a2,a3)) result = 0.;

      for(unsigned j = 0; j < POL1::order; ++j)
        for(unsigned i = 0; i < POL2::order; ++i)
          for(unsigned k = 0; k < POL3::order; ++k)
            result += dx1 * w1[j] * (dx2 * w2[i] * ( dx3 * w3[k] * fun(z1[j], z2[i], z3[k])));

      return result;
    }
};
