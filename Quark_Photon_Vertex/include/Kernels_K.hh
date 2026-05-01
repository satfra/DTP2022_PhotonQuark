#pragma once

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <momentumtransform.hh>

class K
{
  private:
    double y, l2, u, uprime, V, w, wprime, X;

    // 12×12 sparsity mask — true means K_{ij}(...) is identically zero in the
    // basis. Built once at compile time so the BSE hot loop does a flat table
    // lookup instead of walking the explicit if-chain.
    static constexpr std::array<bool, 144> kZeroMask = []() {
      std::array<bool, 144> m{};
      for (unsigned i = 0; i < 12; ++i) {
        for (unsigned j = 0; j < 12; ++j) {
          bool zero;
          if ((i == 0 && j == 1) || (i == 1 && j == 0) ||
              (i == 2 && j == 3) || (i == 3 && j == 2) ||
              (i == 3 && j == 4) || (i == 4 && j == 3) ||
              (i == 4 && j == 5) || (i == 5 && j == 4) ||
              (i == 6 && j == 7) || (i == 7 && j == 6) ||
              (i == 7 && j == 8) || (i == 8 && j == 7) ||
              (i == 8 && j == 9) || (i == 9 && j == 8) ||
              (i == 11 && j == 10) || (i == 10 && j == 11) ||
              (i == 0 && j == 2) || (i == 2 && j == 0) ||
              (i == 0 && j == 3) || (i == 3 && j == 0) ||
              (i == 0 && j == 4) || (i == 4 && j == 0) ||
              (i == 0 && j == 7) || (i == 7 && j == 0) ||
              (i == 1 && j == 3) || (i == 3 && j == 1) ||
              (i == 1 && j == 4) || (i == 4 && j == 1) ||
              (i == 1 && j == 5) || (i == 5 && j == 1) ||
              (i == 1 && j == 6) || (i == 6 && j == 1) ||
              (i == 2 && j == 4) || (i == 4 && j == 2) ||
              (i == 2 && j == 5) || (i == 5 && j == 2) ||
              (i == 2 && j == 6) || (i == 6 && j == 2) ||
              (i == 3 && j == 5) || (i == 5 && j == 3) ||
              (i == 3 && j == 6) || (i == 6 && j == 3) ||
              (i == 3 && j == 7) || (i == 7 && j == 3) ||
              (i == 4 && j == 6) || (i == 6 && j == 4) ||
              (i == 4 && j == 7) || (i == 7 && j == 4) ||
              (i == 5 && j == 7) || (i == 7 && j == 5) ||
              (i == 8 && j == 10) || (i == 10 && j == 8) ||
              (i == 8 && j == 11) || (i == 11 && j == 8) ||
              (i == 9 && j == 11) || (i == 11 && j == 9))
            zero = true;
          else if (i < 8 && j < 8)
            zero = false;
          else if (i >= 8 && j >= 8)
            zero = false;
          else
            zero = true;
          m[i * 12 + j] = zero;
        }
      }
      return m;
    }();

  public:
    // Hot-path constructor used by precalculate_K_kernel. Caller hoists the
    // y-independent quantities sqrt(k²·k'²), sqrt(1-z²), sqrt(1-z'²),
    // sqrt(k²), sqrt(k'²) so the y-quadrature lambda runs them once across
    // the entire grid instead of once per integrand evaluation.
    K(const double& k_sq, const double& k_sq_prime,
      const double& z, const double& z_prime,
      const double& y_,
      const double& k_kp_sqrt,
      const double& s_z, const double& s_z_prime,
      const double& k_v, const double& k_prime_v)
    {
      y = y_;
      l2 = k_sq + k_sq_prime
            - 2. * k_kp_sqrt * (z * z_prime + y * s_z * s_z_prime);
      u = k_v * s_z;
      uprime = k_prime_v * s_z_prime;
      const double inv_l2 = 1.0 / l2;
      V = (k_v * z - k_prime_v * z_prime) * inv_l2;
      w = u * u * inv_l2;
      wprime = uprime * uprime * inv_l2;
      X = u * uprime * inv_l2;
    }

    // Backwards-compatible constructor — delegates to the hot-path version
    // after computing the precomputed quantities itself. q_sq is unused but
    // retained for API stability.
    K(const double& k_sq, const double& k_sq_prime, const double& z,
      const double& z_prime, const double& y_, const double& /*q_sq*/)
        : K(k_sq, k_sq_prime, z, z_prime, y_,
            std::sqrt(k_sq * k_sq_prime),
            std::sqrt(1. - z * z),
            std::sqrt(1. - z_prime * z_prime),
            std::sqrt(k_sq),
            std::sqrt(k_sq_prime))
    {}

    static bool isZeroIndex(const unsigned& i, const unsigned& j)
    {
      if (i > 11 || j > 11)
        throw std::runtime_error("Function get(..) out of range in Kernels_K");
      return kZeroMask[i * 12 + j];
    }

    // Returns K_{ij}(y, l², u, u', V, w, w', X).
    // Implements the Bethe-Salpeter kernel in the 12-component basis.
    // See Eqs. (51)-(52) in the project description.
    double get(const unsigned& i, const unsigned& j) const
    {
      if (i > 11 || j > 11)
        throw std::runtime_error("Function get(..) out of range in Kernels_K");

      const unsigned i_prime = i+1;
      const unsigned j_prime = j+1;
      const unsigned super_idx = 100 * i_prime + j_prime;

      switch(super_idx)
      {
        case 101:
          return K11(y, l2, u, uprime, V, w, wprime, X);
        case 202:
          return K22(y, l2, u, uprime, V, w, wprime, X);
        case 303:
          return K33(y, l2, u, uprime, V, w, wprime, X);
        case 404:
          return K44(y, l2, u, uprime, V, w, wprime, X);
        case 505:
          return K55(y, l2, u, uprime, V, w, wprime, X);
        case 606:
          return K66(y, l2, u, uprime, V, w, wprime, X);
        case 707:
          return K77(y, l2, u, uprime, V, w, wprime, X);
        case 808:
          return K88(y, l2, u, uprime, V, w, wprime, X);
        case 106:
          return K16(y, l2, u, uprime, V, w, wprime, X);
        case 601:
          return K61(y, l2, u, uprime, V, w, wprime, X);
        case 107:
          return K17(y, l2, u, uprime, V, w, wprime, X);
        case 701:
          return K71(y, l2, u, uprime, V, w, wprime, X); 
        case 607:
          return K67(y, l2, u, uprime, V, w, wprime, X);
        case 706:
          return K76(y, l2, u, uprime, V, w, wprime, X);
        case 203:
          return K23(y, l2, u, uprime, V, w, wprime, X);
        case 302:
          return K32(y, l2, u, uprime, V, w, wprime, X);
        case 208:
          return K28(y, l2, u, uprime, V, w, wprime, X);
        case 802:
          return K82(y, l2, u, uprime, V, w, wprime, X);
        case 308:
          return K38(y, l2, u, uprime, V, w, wprime, X);
        case 803:
          return K83(y, l2, u, uprime, V, w, wprime, X);

        case 909:
          return K99(y, l2, u, uprime, V, w, wprime, X);
        case 1010:
          return K1010(y, l2, u, uprime, V, w, wprime, X);
        case 1011:
          return K1011(y, l2, u, uprime, V, w, wprime, X);
        case 1110:
          return K1110(y, l2, u, uprime, V, w, wprime, X);
        case 1111:
          return K1111(y, l2, u, uprime, V, w, wprime, X);
        case 1212:
          return K1212(y, l2, u, uprime, V, w, wprime, X);
      }
      return 0.;
    }

  private:
    double K11(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return -(1. + y * y) / 2. - y * (1. - y * y) * X;
    }

    double K22(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return -(1. + y * y) / 2. * (1. - 2. * l2 * V * V) + y * (1. - y * y) * X;
    }

    double K33(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return  y * (1. - 2. * l2 * V * V) - (1. - y * y) * X;
    }

    double K44(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return y + (1. - y * y) * X;
    }

    double K55(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return 3. * y;
    }

    double K66(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return -y * (1. + 2. * l2 * V * V);
    }

    double K77(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return -y * y * (3. - 2. * l2 * V*V) + 2. * y * (1. - y * y) * X;
    }

    double K88(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return y * y - 2. * y * (1. - y * y) * X;
    }

    double K16(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return std::sqrt(2.) * (1. - y * y) * uprime * V;
    }

    double K61(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return -std::sqrt(2.) * (1. - y * y) * u * V;
    }

    double K17(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return -(1. - y * y) / std::sqrt(2.) * (1. + 2. * wprime - 2. * y * X);
    }

    double K71(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return -(1. - y * y) / std::sqrt(2.) * (1. + 2 * w - 2. * y * X);
    }

    double K23(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return (2. * y * u - (1. + y * y) * uprime) * V;
    }

    double K32(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return -(2. * y * uprime - (1. + y * y) * u) * V;
    }

    double K67(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return 2. * y * (uprime - y * u) * V;
    }

    double K76(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return -2. * y * (u - y * uprime) * V;
    }

    double K28(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return K71(y, l2, u, uprime, V, w, wprime, X)
        + std::sqrt(2.) * (1. - y * y);
    }

    double K82(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return K17(y, l2, u, uprime, V, w, wprime, X)
        + std::sqrt(2.) * (1. - y * y);
    }

    double K38(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return -K61(y, l2, u, uprime, V, w, wprime, X);
    }

    double K83(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return -K16(y, l2, u, uprime, V, w, wprime, X);
    }

    double K99(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return K55(y, l2, u, uprime, V, w, wprime, X) / y;
    }

    double K1010(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return K66(y, l2, u, uprime, V, w, wprime, X) / y;
    }

    double K1011(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return K67(y, l2, u, uprime, V, w, wprime, X) / y;
    }

    double K1110(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return K76(y, l2, u, uprime, V, w, wprime, X) / y;
    }

    double K1111(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return K77(y, l2, u, uprime, V, w, wprime, X) / y;
    }

    double K1212(const double& y, const double& l2, const double& u, const double& uprime, const double& V, const double& w, const double& wprime, const double& X) const
    {
      return K88(y, l2, u, uprime, V, w, wprime, X) / y;
    }
};
