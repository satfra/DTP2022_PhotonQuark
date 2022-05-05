#pragma once

#include <vector>
#include <algorithm>
#include <Utils.hh>

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
