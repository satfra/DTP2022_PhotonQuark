#pragma once

#include <cmath>
#include <optional>

#include "Utils.hh"
#include "types.hh"
#include "parameters.hh"
#include "LinearInterpolate.hh"

#include "quark_dse.hh"

class quark_model
{
  private:
    static constexpr double z_2 = 0.97;
    static constexpr double ScaleFactor_AM = 1./(0.7*0.7);

  public:
    // The model for A(x) given in the project description
    double A(const double& p_sq) const
    {
      const double x = ScaleFactor_AM * p_sq;
      return 0.95 + 0.3 / std::log(x + 2.0) + 0.1 / (1.0 + x) + 0.29 * std::exp(-0.1 * x)
        - 0.18 * std::exp(-3.0 * x);
    }

    // The model for B(x) given in the project description
    double M(const double& p_sq) const
    {
      const double x = ScaleFactor_AM * p_sq;
      return 0.06 / (1.0 + x) + 0.44 * std::exp(-0.66 * x) + 0.009 / pow(std::log(x + 2.0), 0.48);
    }

    double z2() const
    {
      return z_2;
    }
};

class quark_DSE
{
  public:
    double A(const double& p_sq) const
    {
      return (*ip_a)(std::log(p_sq));
    }

    double M(const double& p_sq) const
    {
      const double q = std::log(p_sq);
      return (*ip_b)(q) / (*ip_a)(q);
    }

    double z2() const
    {
      return quark_z2;
    }

    quark_DSE()
    {
      const mat_double quark_a_and_b = quark_iterate_dressing_functions(
          parameters::physical::quark_a0,
          parameters::physical::m_c,
          parameters::physical::m_c,
          parameters::physical::mu);
      quark_a    = quark_a_and_b[0];
      quark_b    = quark_a_and_b[1];
      quark_grid = quark_a_and_b[2]; // logarithmic grid in p²
      quark_z2   = quark_a_and_b[3][0];
      // Construct interpolators once, referencing the now-stable member vectors
      ip_a.emplace(quark_grid, quark_a);
      ip_b.emplace(quark_grid, quark_b);
    }

  private:
    vec_double quark_a, quark_b, quark_grid;
    double quark_z2 = 0.;
    // Interpolators stored as members so they are constructed only once.
    // std::optional is used because lInterpolator holds references to the
    // member vectors, which are not ready until the constructor body runs.
    std::optional<lInterpolator<double, double>> ip_a, ip_b;
};
