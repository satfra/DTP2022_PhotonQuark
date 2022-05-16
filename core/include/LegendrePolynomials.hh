#pragma once

#include <cmath>
#include <Utils.hh>
#include <PolynomialBase.hh>

template<unsigned o, typename _RF = double>
class LegendrePolynomial : public PolynomialBase<LegendrePolynomial<o>,_RF>
{
  using PB = PolynomialBase<LegendrePolynomial<o>,_RF>;
  public:
    static constexpr unsigned order = o;
    using RF = _RF;
    LegendrePolynomial() 
    {
      PB::calcZeros();
      PB::calcWeights();
    }

  protected:
    virtual RF P(const RF& x) const override
    {
      std::vector<RF> res(order+1);
      res[0] = 1.;
      res[1] = x;
      for(unsigned n = 2; n <= order; ++n)
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
      for(unsigned n = 2; n <= order; ++n)
      {
          res[n] = (2.*n-1)/RF(n) * x * res[n-1] - (n-1.)/RF(n) * res[n-2];
          der[n] = n * res[n-1] + x * der[n-1];
      }
      return der[order];
    }

    virtual RF initialGuess(unsigned j) const override
    {
      return std::cos(M_PI * (RF(j+1) - 0.25) / (order + 0.5));
    }
};
