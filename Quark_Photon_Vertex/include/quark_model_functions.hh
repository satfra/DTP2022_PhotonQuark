#pragma once

#include <cmath>

#include "Utils.hh"
#include "types.hh"
#include "parameters.hh"
#include "spline.h"

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
      return ip_a(std::log(p_sq));
    }

    double M(const double& p_sq) const
    {
      const double q = std::log(p_sq);
      return ip_b(q) / ip_a(q);
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
      // Cubic spline (vs piecewise-linear) is required so that the WTI
      // finite-difference Δ_A = (A(k+²)−A(k−²))/(k+²−k−²) does not pick up
      // step jumps at DSE-grid nodes — see phd-work-horak commit c4f017f.
      ip_a.set_points(quark_grid, quark_a);
      ip_b.set_points(quark_grid, quark_b);
    }

  private:
    vec_double quark_a, quark_b, quark_grid;
    double quark_z2 = 0.;
    tk::spline ip_a, ip_b;
};
