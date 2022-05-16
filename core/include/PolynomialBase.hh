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

    static constexpr RF __PRECISION = 1e-6;
    static constexpr LI __MAXSTEPS = 1e+5;

    RF getZero(const RF& x_0) const
    {
      RF error(1.);
      unsigned long stepnumber(0);
      RF x_old = x_0;
      RF x_new = x_0;
      while (stepnumber < 3 || (error > __PRECISION && stepnumber < __MAXSTEPS))
      {
        RF buf = x_new;
        x_new = x_old - P(x_old) / dP(x_old);
        x_old = buf;
        error = std::abs(x_old - x_new);
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
