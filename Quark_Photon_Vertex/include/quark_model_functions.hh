#pragma once

#include <cmath>
#include "parameters.hh"
#include "quark_dse.hh"
#include "Utils.hh"
#include "LinearInterpolate.hh"

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
    // The model for A(x) given in the project description
    double A(const double& p_sq) const
    {
      const double q = std::log(p_sq);
      lInterpolator ip_a(quark_grid, quark_a);
      return ip_a(q);
    }

    // The model for B(x) given in the project description
    double M(const double& p_sq) const
    {
      const double q = std::log(p_sq);
      lInterpolator ip_a(quark_grid, quark_a);
      lInterpolator ip_b(quark_grid, quark_b);
      return ip_b(q)/ip_a(q);
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
      quark_a = quark_a_and_b[0];
      quark_b = quark_a_and_b[1];
      quark_grid = quark_a_and_b[2]; // Logarithmic grid in p^2
      quark_z2 = quark_a_and_b[3][0];
    }

  private:
    vec_double quark_a, quark_b, quark_grid;
    double quark_z2;
};
