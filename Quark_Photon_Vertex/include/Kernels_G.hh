#pragma once

#include <array>
#include <complex>
#include <cstdint>
#include <iostream>
#include <stdexcept>

#include "quark_model_functions.hh"

template<typename Quark>
class G
{
  private:
    static constexpr std::complex<double> II = {0.0, 1.0};
    static constexpr std::complex<double> zero = {0.0, 0.0};

    using submatrix = std::array<std::array<std::complex<double>, 4>, 4>;
    submatrix G_upper{};
    submatrix G_middle{};
    submatrix G_lower{};

    // Block-tensor selector that maps (i,j) ∈ [0,12)² to one of the ten
    // distinct G̃_{ij} components computed once in the constructor. The
    // dispatch matches kernel_G_snake's switch in the previous version.
    static std::complex<double> select_g_snake(
        unsigned i, unsigned j,
        double g11, std::complex<double> g12, std::complex<double> g13, double g14,
        double g22, double g23, std::complex<double> g24,
        double g33, std::complex<double> g34, double g44)
    {
      const unsigned super_idx = 100u * (i + 1u) + (j + 1u);
      switch (super_idx) {
        case 101: case 808: case 1212:
          return g11;
        case 102: case 201: case 708: case 807: case 1211: case 1112:
          return g12;
        case 103: case 301:
          return g13;
        case 608: case 806: case 1012: case 1210:
          return -g13;
        case 104: case 401:
          return g14;
        case 508: case 805: case 912: case 1209:
          return -g14;
        case 202: case 707: case 1111:
          return g22;
        case 203: case 302:
          return g23;
        case 607: case 706: case 1011: case 1110:
          return -g23;
        case 204: case 402:
          return g24;
        case 507: case 705: case 911: case 1109:
          return -g24;
        case 303: case 606: case 1010:
          return g33;
        case 304: case 403: case 506: case 605: case 910: case 1009:
          return g34;
        case 404: case 505: case 909:
          return g44;
        default:
          return 0.0;
      }
    }

  public:
    // Implements the 12×12 G_{ij} kernel (G̃ + propagator factor) from
    // Eqs. (56)-(57) of the project description. Every (k_sq, z, q²) point
    // produces ten distinct G̃ components — the previous version recomputed
    // them all 48× per construction (3 sub-blocks × 16 entries), and called
    // sigma_v / quark.A / quark.M four times each. Now: precompute once, then
    // a single switch dispatch fills the three 4×4 sub-blocks.
    G(const double& k_sq, const double& z, const double& q_sq, const Quark& quark)
    {
      const double sqrt_kq = std::sqrt(k_sq * q_sq);
      const double kp_sq = k_sq + 0.25 * q_sq + sqrt_kq * z;
      const double km_sq = k_sq + 0.25 * q_sq - sqrt_kq * z;

      const double m_kp = quark.M(kp_sq);
      const double m_km = quark.M(km_sq);
      const double sig_kp = 1.0 / (quark.A(kp_sq) * (kp_sq + m_kp * m_kp));
      const double sig_km = 1.0 / (quark.A(km_sq) * (km_sq + m_km * m_km));
      const double sig_factor = -sig_kp * sig_km;

      const double sigma_m = 0.5 * (m_kp + m_km);
      const double delta_m = (m_kp - m_km) / (kp_sq - km_sq);
      const double m_bar_sq = m_kp * m_km;

      const double q = std::sqrt(q_sq);
      const double k = std::sqrt(k_sq);
      const double z_sq = z * z;
      const double y = std::sqrt(1.0 - z_sq);

      const double                g11 = m_bar_sq + k_sq - 0.25 * q_sq;
      const std::complex<double>  g12 = II * q * (sigma_m - 2.0 * k_sq * z_sq * delta_m);
      const std::complex<double>  g13 = -2.0 * II * k_sq * q * z * y * delta_m;
      const double                g14 = -k * q * y;
      const double                g22 = m_bar_sq - (1.0 - 2.0 * z_sq) * k_sq - 0.25 * q_sq;
      const double                g23 = 2.0 * k_sq * z * y;
      const std::complex<double>  g24 = 2.0 * II * k * y * sigma_m;
      const double                g33 = m_bar_sq + (1.0 - 2.0 * z_sq) * k_sq + 0.25 * q_sq;
      const std::complex<double>  g34 = II * k * z * (q_sq * delta_m - 2.0 * sigma_m);
      const double                g44 = m_bar_sq - k_sq + 0.25 * q_sq;

      for (unsigned i = 0; i < 4; ++i) {
        for (unsigned j = 0; j < 4; ++j) {
          G_upper[i][j]  = sig_factor * select_g_snake(i,     j,     g11, g12, g13, g14, g22, g23, g24, g33, g34, g44);
          G_middle[i][j] = sig_factor * select_g_snake(i + 4, j + 4, g11, g12, g13, g14, g22, g23, g24, g33, g34, g44);
          G_lower[i][j]  = sig_factor * select_g_snake(i + 8, j + 8, g11, g12, g13, g14, g22, g23, g24, g33, g34, g44);
        }
      }
    }

    const std::complex<double>& get(const unsigned& i, const unsigned& j) const
    {
      if (i > 11 || j > 11)
        throw std::runtime_error("Function get(..) out of range in Kernels_G");
      else if (i < 4 && j < 4)
        return G_upper[i][j];
      else if (i >= 4 && i < 8 && j >= 4 && j < 8)
        return G_middle[i - 4][j - 4];
      else if (i >= 8 && i < 12 && j >= 8 && j < 12)
        return G_lower[i - 8][j - 8];
      return zero;
    }
};
