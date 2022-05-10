#pragma once

#include <cmath>
#include <Utils.hh>
#include <PolynomialBase.hh>

template<unsigned o, typename _RF = double>
class LegendrePolynomial : public PolynomialBase<LegendrePolynomial<o>,_RF>
{
  public:
    static constexpr unsigned order = o;
    using RF = _RF;
    LegendrePolynomial() {}

  private:
    template<unsigned n>
      RF _P(const RF& x) const
      {
        if constexpr (n == 0)
          return 1.;
        else if constexpr (n == 1)
          return x;
        else
          return (2.*n-1)/RF(n) * x * _P<n-1>(x) - (n-1.)/RF(n) * _P<n-2>(x);
      }

    template<unsigned n>
      RF _dP(const RF& x) const
      {
        if constexpr (n == 0)
          return 0.;
        else if constexpr (n == 1)
          return 1.;
        else
          return ( (n+1.) * x * _P<n>(x) - (n+1.) * _P<n+1>(x) ) / (1. - powr<2>(x));
      }
  protected:
    virtual RF P(const RF& x) const override
    {
      std::vector<RF> res(order+1);
      res[0] = 1.;
      res[1] = x;
      for(unsigned int n = 2; n <= order; ++n)
          res[n] = (2.*n-1)/RF(n) * x * res[n-1] - (n-1.)/RF(n) * res[n-2];
      return res[order];
    }

    virtual RF dP(const RF& x) const override
    {
      std::vector<RF> der(order+1);
      std::vector<RF> res(order+1);
      res[0] = 1.;
      der[0] = 0.;
      res[1] = x;
      der[1] = 1.;
      for(unsigned int n = 2; n <= order; ++n)
      {
          res[n] = (2.*n-1)/RF(n) * x * res[n-1] - (n-1.)/RF(n) * res[n-2];
          der[n] = n * res[n-1] + x * der[n-1];
      }
      return der[order];
    }
};
