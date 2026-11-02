#pragma once

#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Utils.hh"
#include "types.hh"
#include "hdf5IO.hh"
#include "iteration.hh"
#include "parameters.hh"
#include "QuadratureIntegral.hh"
#include "LegendrePolynomials.hh"
#include "ChebyshevPolynomial2.hh"
#include "complex_spline.hh"

namespace hvp
{
  // Computes the renormalised hadronic-vacuum-polarisation loop on the
  // existing q_grid. The b's already contain the quark propagators
  // (b ≡ S·Γ·S, see literature/3-Quark-Photon-Vertex.pdf), so the
  // integrand is simply the Dirac trace times the 4D Euclidean measure
  // — no σ_v / A / M factors are added here.
  //
  // From the notebook UnneccecaryNotebook.nb the trace contracted with
  // γ^μ collapses (in the b-basis) to
  //   tr = 4·√2·b₁ + 4·b₇ + 4·b₁₀
  // ↔ b[0], b[6], b[9] in the code's 0-indexed tensor.
  //
  // The integrand subtracts tr(k, p²_min) at p²_min = q_grid[0] so the
  // quadratically divergent parts cancel pointwise; the remaining constant
  // is absorbed into C by fit_low_p_sq.
  //
  // Quadrature: Legendre in log(k'²), Cheb2 in z' — the latter coincides
  // with z_grid by construction. b1/b7/b10 are cubic-splined in log(k'²)
  // at each z'-grid index (the same fix applied to a_iteration_step in
  // the BSE; linear interp was the dominant integrand error). Direct
  // (zp_idx, k_quad_idx) loops let us index splines by zp_idx instead of
  // a value-to-index lookup.
  inline vec_double compute_hvp_pi(
      const tens_cmplx& b1, const tens_cmplx& b7, const tens_cmplx& b10,
      const vec_double& q_grid, const vec_double& k_grid, const vec_double& /*z_grid*/)
  {
    using namespace parameters::numerical;

    vec_double pi_values(q_steps, 0.0);

    constexpr double pref_b1  = 4.0 * M_SQRT2;
    constexpr double pref_b7  = 4.0;
    constexpr double pref_b10 = 4.0;

    constexpr unsigned q_ren = 0;
    // b*[q_ren] returns a Tensor3 row view (Row2D); take it by value
    // since it's a small pointer+strides struct, not a heap-backed vector.
    const auto b1_ren  = b1[q_ren];
    const auto b7_ren  = b7[q_ren];
    const auto b10_ren = b10[q_ren];

    // Legendre / Cheb2 quadrature data (built once).
    static const LegendrePolynomial<k_steps> lp_k;
    static const auto& k_unit_zeros   = lp_k.zeroes();
    static const auto& k_unit_weights = lp_k.weights();
    static const ChebyshevPolynomial2<z_steps> cp_z;
    static const auto z_weights = cp_z.weights();

    const double k_a = k_grid.front();
    const double k_b = k_grid.back();
    const double dk = 0.5 * (k_b - k_a);
    const double k_mid = 0.5 * (k_a + k_b);

    vec_double k_quad_log(k_steps), k_quad_powr2(k_steps);
    for (unsigned q = 0; q < k_steps; ++q) {
      k_quad_log[q] = k_mid + dk * k_unit_zeros[q];
      const double k_sq = std::exp(k_quad_log[q]);
      k_quad_powr2[q] = k_sq * k_sq;
    }

    // Build z_steps cubic splines in log(k'²) from a 2D (k, z) view.
    // Works for any view exposing [k_idx][zp_idx] (Row2D from types.hh).
    auto build_splines = [&](const auto& view2d) {
      std::vector<ComplexSpline> splines(z_steps);
      for (unsigned zp = 0; zp < z_steps; ++zp) {
        vec_cmplx slice(k_steps);
        for (unsigned k = 0; k < k_steps; ++k)
          slice[k] = view2d[k][zp];
        splines[zp].set_points(k_grid, slice);
      }
      return splines;
    };

    // Renormalisation-point splines: built once, reused across all q_iter.
    const auto sp_b1_ren  = build_splines(b1_ren);
    const auto sp_b7_ren  = build_splines(b7_ren);
    const auto sp_b10_ren = build_splines(b10_ren);

    #pragma omp parallel for
    for (unsigned q_iter = 0; q_iter < q_steps; ++q_iter)
    {
      const auto sp_b1  = build_splines(b1[q_iter]);
      const auto sp_b7  = build_splines(b7[q_iter]);
      const auto sp_b10 = build_splines(b10[q_iter]);

      std::complex<double> integral{0., 0.};
      for (unsigned zp = 0; zp < z_steps; ++zp)
      {
        std::complex<double> partial{0., 0.};
        for (unsigned q = 0; q < k_steps; ++q)
        {
          const double log_k_sq = k_quad_log[q];

          const std::complex<double> trace_p =
              pref_b1  * sp_b1 [zp](log_k_sq) +
              pref_b7  * sp_b7 [zp](log_k_sq) +
              pref_b10 * sp_b10[zp](log_k_sq);

          const std::complex<double> trace_ren =
              pref_b1  * sp_b1_ren [zp](log_k_sq) +
              pref_b7  * sp_b7_ren [zp](log_k_sq) +
              pref_b10 * sp_b10_ren[zp](log_k_sq);

          // 4D Euclidean measure with log-k² Jacobian (k'⁴ from
          // (1/2)·k² dk² → (k²)² d(log k²)). √(1−z'²) is the Cheb2 weight,
          // baked into z_weights; do not include it here. The integrand is
          // independent of the remaining angles: ∫dφ = 2π and ∫dy = 2.
          const double measure = int_factors * 4. * M_PI * k_quad_powr2[q];
          partial += k_unit_weights[q] * (trace_p - trace_ren) * measure;
        }
        partial *= dk;
        integral += z_weights[zp] * partial;
      }
      pi_values[q_iter] = integral.real();
    }

    return pi_values;
  }

  // Low-p² fit  T(p²) = C + s·p² + O(p⁴)  of the traced loop T computed by
  // compute_hvp_pi. C is the (regulator-dependent) quadratic divergence left
  // in the trace, s = −3Π(0)/Z₂ the log-divergent Π(0). A cubic in p² is
  // least-squares fitted to all points with p² ≤ hvp_fit_max_p_sq; p² is
  // rescaled to [0, 1] for conditioning.
  struct LowPSqFit { double C, s; };

  inline LowPSqFit fit_low_p_sq(const vec_double& p_sq, const vec_double& T)
  {
    using parameters::numerical::hvp_fit_max_p_sq;
    constexpr unsigned n_coeff = 4;

    unsigned n = 0;
    while (n < p_sq.size() && p_sq[n] <= hvp_fit_max_p_sq) ++n;
    if (n < n_coeff + 2)
      throw std::runtime_error("hvp: need at least 6 q-grid points with p² <= "
                               "hvp_fit_max_p_sq for the p² -> 0 fit");
    const double scale = p_sq[n - 1];

    // normal equations  (VᵀV) c = Vᵀ T  for the Vandermonde matrix V
    double M[n_coeff][n_coeff + 1] = {};
    for (unsigned r = 0; r < n; ++r) {
      double pw[n_coeff];
      pw[0] = 1.;
      for (unsigned i = 1; i < n_coeff; ++i) pw[i] = pw[i - 1] * p_sq[r] / scale;
      for (unsigned i = 0; i < n_coeff; ++i) {
        for (unsigned j = 0; j < n_coeff; ++j) M[i][j] += pw[i] * pw[j];
        M[i][n_coeff] += pw[i] * T[r];
      }
    }
    // Gaussian elimination with partial pivoting
    for (unsigned c = 0; c < n_coeff; ++c) {
      unsigned piv = c;
      for (unsigned r = c + 1; r < n_coeff; ++r)
        if (std::abs(M[r][c]) > std::abs(M[piv][c])) piv = r;
      std::swap(M[c], M[piv]);
      for (unsigned r = 0; r < n_coeff; ++r) {
        if (r == c) continue;
        const double f = M[r][c] / M[c][c];
        for (unsigned j = c; j <= n_coeff; ++j) M[r][j] -= f * M[c][j];
      }
    }
    return {M[0][n_coeff] / M[0][0], M[1][n_coeff] / M[1][1] / scale};
  }

  // Top-level driver: take b₁,b₇,b₁₀ (in memory, from iterate_a_and_b),
  // compute the traced loop T(p²), and store the renormalised HVP
  //   Π̂(p²) = Π(p²) − Π(0) = (Z₂/3) · [ s − (T(p²) − C) / p² ]
  // (one colour, unit charge; see fit_low_p_sq for C and s) under /hvp.
  // Normalisation: Π_μν = −Z₂ ∫ tr[γ_μ S Γ_ν S] = (p²δ_μν − p_μp_ν) Π(p²), and
  // T = ∫ tr[γ_μ S Γ_μ S] = −3p²Π/Z₂ + C. At large p², Π̂ → ln(p²)/(12π²).
  inline void hvp_driver(const HvpBInput& b_in,
                         const vec_double& q_grid,
                         const vec_double& k_grid,
                         const vec_double& z_grid,
                         const std::string& h5_path)
  {
    std::cout << "\n_____________________________________\n\n"
              << "Calculating the hadronic vacuum polarisation...\n";

    std::cout << " Integrating the loop on the q grid..." << std::flush;
    const vec_double trace =
        compute_hvp_pi(b_in.b1, b_in.b7, b_in.b10, q_grid, k_grid, z_grid);
    std::cout << " done\n";

    const LowPSqFit fit = fit_low_p_sq(q_grid, trace);
    vec_double pi_hat(q_grid.size());
    for (std::size_t i = 0; i < q_grid.size(); ++i)
      pi_hat[i] = b_in.z2 / 3. * (fit.s - (trace[i] - fit.C) / q_grid[i]);
    const double pi_at_zero = -b_in.z2 / 3. * fit.s;

    qpv_hdf5::append_hvp(h5_path, q_grid, pi_hat, pi_at_zero, fit.C);

    std::cout << " HVP written to " << h5_path << " (/hvp)"
              << "  (Pi(p^2) - Pi(0); unrenormalised Pi(0) = " << pi_at_zero << ")\n";
  }
} // namespace hvp
