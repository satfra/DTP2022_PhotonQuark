#pragma once

#include <cmath>

#include "Utils.hh"
#include "types.hh"
#include "parameters.hh"
#include "flavor.hh"
#include "spline.h"

#include "quark_dse.hh"

class quark_model
{
  private:
    static constexpr double z_2 = 0.97;
    static constexpr double ScaleFactor_AM = 1./(0.7*0.7);
    // Maris-Tandy IR parameters carried from the selected flavour. The
    // analytic model is mass-independent, so m_c/mu/a0 are ignored here —
    // this is why non-light flavours require the quark DSE (-d).
    double eta_mt_, lambda_mt_;

  public:
    explicit quark_model(const FlavorParams& fp)
      : eta_mt_(fp.eta_mt), lambda_mt_(fp.lambda_mt) {}

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

    double z2() const { return z_2; }
    double eta_mt() const { return eta_mt_; }
    double lambda_mt() const { return lambda_mt_; }
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

    double z2() const { return quark_z2; }
    double eta_mt() const { return eta_mt_; }
    double lambda_mt() const { return lambda_mt_; }

    // Raw DSE solution, for output. quark_grid holds log(p²); quark_a/quark_b
    // are A(p²) and B(p²) sampled on that grid.
    const vec_double& dse_log_p_sq() const { return quark_grid; }
    const vec_double& dse_A() const { return quark_a; }
    const vec_double& dse_B() const { return quark_b; }

    explicit quark_DSE(const FlavorParams& fp)
      : eta_mt_(fp.eta_mt), lambda_mt_(fp.lambda_mt)
    {
      const mat_double quark_a_and_b = quark_iterate_dressing_functions(
          fp.quark_a0,
          fp.m_c,
          fp.m_c,
          fp.mu,
          fp.eta_mt,
          fp.lambda_mt);
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
    double eta_mt_, lambda_mt_;
    tk::spline ip_a, ip_b;
};
