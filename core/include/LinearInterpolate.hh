#pragma once

#include <vector>

template<typename _RF = double>
class lInterpolator
{
    using LI = unsigned long long;
  public:
    using RF = _RF;
    using Range = std::vector<RF>;

    lInterpolator(const Range& _x, const Range& _f)
      : x(_x), f(_f), a(x.front()), b(x.back())
    {
      
    }

    RF operator()(const RF& y)
    {
      if(y > b || y < a)
        throw std::runtime_error("Interpolating outside bounds");

      LI idx = x.size() + 1;
      for(LI i = 0; i < x.size()-1; ++i)
        if(y >= x[i] && y <= x[i+1])
        {
          idx = i;
          break;
        }
      
      const RF t = (y - x[idx]) / (x[idx+1] - x[idx]);
      return t*f[idx+1] + (1.-t)*f[idx];
    }

  private:
    Range x,f;
    RF a,b;
};
