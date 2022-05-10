#pragma once

#include <vector>

// TODO: Put me at a nice home please :)
unsigned int locate(const std::vector<double> &xx, double x)
{
    /*
     * This function has been taken and slightly been modified from the
     * "Numerical Recipes in C" book.
     *
     * Given a vector xx of size n and given a value x, it returns an integer j,
     * such that x is between xx[j] and xx[j+1]. xx must be sorted. Returns -1
     * or n when is x is out of range.
     */

    const unsigned int n = xx.size();

    if (x == xx[0]) return 0;
    else if (x == xx[n-1]) return n-2;

    unsigned int ju, jm, jl;
    bool ascnd;

    jl = -1;
    ju = n;

    ascnd = (xx[n-1] >= xx[0]);

    while (ju - jl > 1) {
        jm = (ju + jl) / 2;
        if ((x >= xx[jm]) == ascnd) {
            // TODO: Test, if the innermost parentheses are correct. Might be a source of errors.
            jl = jm;
        } else {
            ju = jm;
        }
    }

    return jl;
}


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

      // Search for y in x
      idx = locate(x, y);

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
