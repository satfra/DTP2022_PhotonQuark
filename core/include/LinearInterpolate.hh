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

template<typename _RF = double>
class lInterpolator2d
{
	using LI = unsigned long long;
public:
	using RF = _RF;
	using Range = std::vector<RF>;
	using Grid = std::vector<std::vector<RF>>;

	lInterpolator2d(const Range& _x1, const Range& _x2, const Grid& _f)
		: x1(_x1),x2(_x2), f(_f), a(x1.front()), b(x1.back()), c(x2.front()), d(x2.back())
	{
	}

	RF operator()(const RF& y,const RF& z)
	{
		if (y > b || y < a || z > d || z < c)
			throw std::runtime_error("Interpolating outside bounds");
	
		LI idx1 = x1.size() + 1;
		for (LI i = 0; i < x1.size() - 1; ++i)
			if (y >= x1[i] && y <= x1[i + 1])
			{
				idx1 = i;
				break;
			}

		LI idx2 = x2.size() + 1;
		for (LI i = 0; i < x2.size() - 1; ++i)
			if (z >= x2[i] && z <= x2[i + 1])
			{
				idx2 = i;
				break;
			}

		const RF t1 = (y - x1[idx1]) / (x1[idx1 + 1] - x1[idx1]);
		const RF t2 = (y - x2[idx2]) / (x2[idx2 + 1] - x2[idx2]);
		return t2 * (t1 *f[idx1 + 1][idx2] + (1. - t1)*f[idx1][idx2]) + (1. - t2)*(t1 *f[idx1 + 1][idx2+1] + (1. - t1)*f[idx1][idx2+1]);
	}

private:
	Range x1, x2;
	RF a, b, c, d;
	Grid f;
};
