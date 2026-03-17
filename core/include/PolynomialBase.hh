#pragma once

#include <vector>
#include <cmath>

#include "Utils.hh"

template<typename POL, typename _RF>
class PolynomialBase
{
  public:
    using RF = _RF;

  protected:
    using LI = unsigned long long;

    static constexpr RF NEWTON_PRECISION = 1e-14;
    static constexpr LI NEWTON_MAX_STEPS = 1e+5;

    RF getZero(const RF& x_0) const
    {
      RF error(1.);
      unsigned long stepnumber(0);
      RF x_new = x_0;
      while (stepnumber < 3 || (error > NEWTON_PRECISION && stepnumber < NEWTON_MAX_STEPS))
      {
        const RF x_next = x_new - P(x_new) / dP(x_new);
        error = std::abs(x_next - x_new);
        x_new = x_next;
        stepnumber++;
      }
      return x_new;
    }

    std::vector<RF> _z;
    std::vector<RF> _w;

    virtual void calcZeros()
    {
      _z.resize(order);
      for(unsigned i = 0; i < order; ++i)
      {
        unsigned j = order - 1 - i;
        const RF iG = initialGuess(j);
        _z[i] = getZero(iG);
      }
    }

    virtual void calcWeights()
    {
      if(_z.size() == 0)
        calcZeros();
      _w.resize(order);
      for(unsigned i = 0; i < order; ++i)
      {
        const RF& x_i = _z[i];
        _w[i] = 2. / ((1. - powr<2>(x_i)) * powr<2>(dP(x_i)));
      }
    }

    virtual RF P(const RF& x) const = 0;
    virtual RF dP(const RF& x) const = 0;
    virtual RF initialGuess(unsigned j) const = 0;

  public:
    static constexpr unsigned order = POL::order;

    PolynomialBase()
    {
    }

    const std::vector<RF>& zeroes() const
    {
      return _z;
    }

    const std::vector<RF>& weights() const
    {
      return _w;
    }

    RF operator()(const RF& x) const
    {
      return P(x);
    }
};
