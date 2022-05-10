#pragma once

#include <vector>
#include <algorithm>
#include <Utils.hh>

#include <tbb/parallel_for.h>

template<typename POL, typename RF = typename POL::RF>
class qIntegral
{
  private:
    POL polynomials;

  public:
    qIntegral() {}
    template<typename FUN>
    RF operator()(FUN& fun, const RF& a, const RF& b)
    {
      const auto& z = linearMapTo(polynomials.zeroes(), RF(-1.), RF(1.), a, b);
      const auto& w = polynomials.weights();
      const RF dx = (b-a)/2.;

      RF res(0.);
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
    qIntegral2d() {}

    template<typename FUN>
    RF operator()(FUN& fun, const RF& a1, const RF& b1, const RF& a2, const RF& b2)
    {
      const auto& z1 = linearMapTo(polynomials_1.zeroes(), RF(-1.), RF(1.), a1, b1);
      const auto& w1 = polynomials_1.weights();
      const RF dx1 = (b1-a1)/2.;

      const auto& z2 = linearMapTo(polynomials_2.zeroes(), RF(-1.), RF(1.), a2, b2);
      const auto& w2 = polynomials_2.weights();
      const RF dx2 = (b2-a2)/2.;

      std::vector<RF> res(POL1::order);
      auto inner_for_loop = [&](size_t j) {
        for(unsigned i = 0; i < POL2::order; ++i)
          res[j] += dx1 * w1[j] * ( dx2 * w2[i] * fun(z1[j], z2[i]));
      };
      tbb::parallel_for(unsigned(0), POL1::order, inner_for_loop);

      RF result = 0.;
      for(unsigned j = 0; j < POL1::order; ++j)
        result += res[j];

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
    qIntegral3d() {}

    template<typename FUN>
    RF operator()(FUN& fun, const RF& a1, const RF& b1, const RF& a2, const RF& b2, const RF& a3, const RF& b3)
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

      std::vector<RF> res(POL1::order);
      auto inner_for_loop = [&](size_t j) {
        for(unsigned i = 0; i < POL2::order; ++i)
          for(unsigned k = 0; k < POL3::order; ++k)
            res[j] += dx1 * w1[j] * (dx2 * w2[j] * ( dx3 * w3[i] * fun(z1[j], z2[i], z3[k])));
      };
      tbb::parallel_for(unsigned(0), POL1::order, inner_for_loop);

      RF result = 0.;
      for(unsigned j = 0; j < POL1::order; ++j)
        result += res[j];

      return result;
    }
};
