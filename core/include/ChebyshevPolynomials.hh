#pragma once

#include <cmath>
#include <Utils.hh>
#include <PolynomialBase.hh>

template<unsigned o, typename _RF = double>
class ChebyshevPolynomial : public PolynomialBase<ChebyshevPolynomial<o>,_RF>
{
  public:
    static constexpr unsigned order = o;
    using RF = _RF;
    ChebyshevPolynomial() {}

  protected:
    virtual RF P(const RF& x) const override
    {
      std::vector<RF> res(order+1);
      res[0] = 1.;
      res[1] = x;
      for(unsigned n = 2; n <= order; ++n)
          res[n] = 2.*x * res[n-1] - res[n-2];
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
          res[n] = 2.*x * res[n-1] - res[n-2];
          der[n] = 2.* res[n-1] + 2.*x * der[n-1] - der[n-2];
      }
      return der[order];
    }
};
