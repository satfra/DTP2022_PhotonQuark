#pragma once

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <momentumtransform.hh>
#include <acc.hh>

// 12×12 sparsity mask — true means K_{ij}(...) is identically zero in the basis.
// Built once at compile time so the BSE hot loop does a flat table lookup
// instead of walking the explicit if-chain. Kept at namespace scope (not as a
// static member of K): with nvc -acc, a `static constexpr std::array` data
// member inside a class whose methods are `acc routine seq` made the whole
// class's device method dispatch return zero. This mask is host-only anyway
// (the offloaded K' precalc prunes zero (i,j) via the sparse slot map).
inline constexpr std::array<bool, 144> kK_ZeroMask = []() {
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

// The Bethe-Salpeter kernel components K_{ij} — the single definition of the
// formulas (Eqs. (51)-(52)); K::get() delegates here. Comments name the
// component each case implements; the previously separate K11..K1212 methods
// have been folded in, with their cross-references (K28 = K71 + √2(1-y²),
// K99 = K55/y, ...) resolved.
//
// This MUST stay a free function with the formulas inlined in the switch: nvc++
// 26.5 miscompiles a device routine whose switch dispatches to *member method*
// calls once there are more than ~20 cases — it silently falls through and
// returns 0, while the identical arithmetic inlined is correct. See the
// standalone reproducer in tools/nvc_openacc_switch_bug.cpp. The eight
// quantities are the ones the K constructor builds; the offloaded caller
// (precalculate_K_kernel) computes them the same way.
ACC_ROUTINE_SEQ
inline double k_component(unsigned i, unsigned j,
    double y, double l2, double u, double uprime,
    double V, double w, double wprime, double X)
{
  const double sqrt2 = 1.4142135623730950488;  // std::sqrt(2.)
  const double omy2 = 1. - y * y;              // (1 - y²)
  const double opy2 = 1. + y * y;              // (1 + y²)
  const unsigned super_idx = 100u * (i + 1u) + (j + 1u);
  switch (super_idx)
  {
    case 101:  return -opy2 / 2. - y * omy2 * X;                          // K11
    case 202:  return -opy2 / 2. * (1. - 2. * l2 * V * V) + y * omy2 * X;  // K22
    case 303:  return  y * (1. - 2. * l2 * V * V) - omy2 * X;             // K33
    case 404:  return  y + omy2 * X;                                     // K44
    case 505:  return  3. * y;                                          // K55
    case 606:  return -y * (1. + 2. * l2 * V * V);                       // K66
    case 707:  return -y * y * (3. - 2. * l2 * V * V) + 2. * y * omy2 * X; // K77
    case 808:  return  y * y - 2. * y * omy2 * X;                        // K88
    case 106:  return  sqrt2 * omy2 * uprime * V;                        // K16
    case 601:  return -sqrt2 * omy2 * u * V;                             // K61
    case 107:  return -omy2 / sqrt2 * (1. + 2. * wprime - 2. * y * X);    // K17
    case 701:  return -omy2 / sqrt2 * (1. + 2. * w - 2. * y * X);         // K71
    case 607:  return  2. * y * (uprime - y * u) * V;                    // K67
    case 706:  return -2. * y * (u - y * uprime) * V;                    // K76
    case 203:  return  (2. * y * u - opy2 * uprime) * V;                 // K23
    case 302:  return -(2. * y * uprime - opy2 * u) * V;                 // K32
    case 208:  return -omy2 / sqrt2 * (1. + 2. * w - 2. * y * X) + sqrt2 * omy2;       // K28 = K71 + √2(1-y²)
    case 802:  return -omy2 / sqrt2 * (1. + 2. * wprime - 2. * y * X) + sqrt2 * omy2;  // K82 = K17 + √2(1-y²)
    case 308:  return  sqrt2 * omy2 * u * V;                             // K38 = -K61
    case 803:  return -sqrt2 * omy2 * uprime * V;                        // K83 = -K16
    case 909:  return  3.;                                             // K99  = K55/y
    case 1010: return -(1. + 2. * l2 * V * V);                           // K1010= K66/y
    case 1011: return  2. * (uprime - y * u) * V;                        // K1011= K67/y
    case 1110: return -2. * (u - y * uprime) * V;                        // K1110= K76/y
    case 1111: return -y * (3. - 2. * l2 * V * V) + 2. * omy2 * X;        // K1111= K77/y
    case 1212: return  y - 2. * omy2 * X;                               // K1212= K88/y
  }
  return 0.;
}

class K
{
  private:
    double y, l2, u, uprime, V, w, wprime, X;

  public:
    // Hot-path constructor used by precalculate_K_kernel. Caller hoists the
    // y-independent quantities sqrt(k²·k'²), sqrt(1-z²), sqrt(1-z'²),
    // sqrt(k²), sqrt(k'²) so the y-quadrature lambda runs them once across
    // the entire grid instead of once per integrand evaluation.
    ACC_ROUTINE_SEQ
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

    // Constexpr noexcept lookup — the hot path skips 80% of (i,j) pairs via
    // this check, so we want it inlined to a single array load. Bounds
    // guarded by an assert in debug builds; production builds trust the
    // caller (always n_structs = 12 in our codebase).
    // Host-only (reads the static constexpr kZeroMask): the offloaded K'
    // precalc tests emptiness via the sparse slot map instead, so this needs
    // no device version.
    static constexpr bool isZeroIndex(unsigned i, unsigned j) noexcept
    {
      return kK_ZeroMask[i * 12 + j];
    }

    // Returns K_{ij}(y, l², u, u', V, w, w', X).
    // Implements the Bethe-Salpeter kernel in the 12-component basis.
    // See Eqs. (51)-(52) in the project description.
    ACC_ROUTINE_SEQ
    double get(unsigned i, unsigned j) const
    {
      // Device code (OpenACC) cannot throw; the hot loop only ever passes
      // valid 0..11 indices, so the bounds guard is host-only.
#ifndef _OPENACC
      if (i > 11 || j > 11)
        throw std::runtime_error("Function get(..) out of range in Kernels_K");
#endif

      return k_component(i, j, y, l2, u, uprime, V, w, wprime, X);
    }

};
