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
#include "fileIO.hh"
#include "parameters.hh"
#include "QuadratureIntegral.hh"
#include "LegendrePolynomials.hh"
#include "ChebyshevPolynomial2.hh"
#include "complex_spline.hh"

namespace hvp
{
  // Reads b_file_idx_<idx>.dat (written by saveToFile_withGrids in fileIO.hh)
  // and returns the b-tensor indexed [q_iter][k_idx][z_idx]. Mirrors the
  // file layout: header line starting with '#', then for each q_sq one
  // block of (k_steps × z_steps) lines "q_sq i k_sq z Re Im", separated
  // from the next block by a blank line.
  inline tens_cmplx load_b_from_file(unsigned idx)
  {
    using namespace parameters::numerical;
    const std::string fname = "b_file_idx_" + std::to_string(idx) + ".dat";
    std::ifstream in(fname);
    if (!in)
      throw std::runtime_error("hvp: cannot open " + fname);

    tens_cmplx out(q_steps, k_steps, z_steps);

    unsigned q_iter = 0;
    unsigned points_in_block = 0;
    std::string line;
    while (std::getline(in, line))
    {
      if (line.empty() || line[0] == '#')
      {
        if (points_in_block > 0)
        {
          if (points_in_block != k_steps * z_steps)
            throw std::runtime_error("hvp: malformed block in " + fname);
          ++q_iter;
          points_in_block = 0;
        }
        continue;
      }

      std::istringstream iss(line);
      double q_sq, k_sq, z, re, im;
      unsigned i_dummy;
      if (!(iss >> q_sq >> i_dummy >> k_sq >> z >> re >> im))
        throw std::runtime_error("hvp: parse error in " + fname);

      const unsigned k_idx = points_in_block / z_steps;
      const unsigned z_idx = points_in_block % z_steps;
      if (q_iter >= q_steps)
        throw std::runtime_error("hvp: more q-blocks in " + fname + " than q_steps");
      out[q_iter][k_idx][z_idx] = std::complex<double>(re, im);
      ++points_in_block;
    }
    if (points_in_block == k_steps * z_steps)
      ++q_iter;

    if (q_iter != q_steps)
      throw std::runtime_error("hvp: " + fname + " has " + std::to_string(q_iter)
          + " q-blocks, expected " + std::to_string(q_steps));
    return out;
  }

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
  // The subtraction Π̃(p²) = ∫ [tr(k,p²) − tr(k,p²_min)] · measure
  // is performed *inside* the integral so the divergent parts cancel
  // at the integrand level, leaving a well-defined finite integrand.
  // The renormalisation point is q_grid[0] (= parameters::numerical::min_q_sq).
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
          // baked into z_weights; do not include it here.
          const double measure = int_factors * 2. * M_PI * k_quad_powr2[q];
          partial += k_unit_weights[q] * (trace_p - trace_ren) * measure;
        }
        partial *= dk;
        integral += z_weights[zp] * partial;
      }
      pi_values[q_iter] = integral.real();
    }

    return pi_values;
  }

  // Top-level driver: load b₁,b₇,b₁₀ from disk, compute Π̃(p²) with
  // the subtraction baked into the integrand, and write hvp.dat.
  inline void hvp_driver(const vec_double& q_grid,
                         const vec_double& k_grid,
                         const vec_double& z_grid)
  {
    std::cout << "\n_____________________________________\n\n"
              << "Calculating the hadronic vacuum polarisation...\n";

    std::cout << " Loading b_1, b_7, b_10 from disk..." << std::flush;
    const tens_cmplx b1  = load_b_from_file(0);
    const tens_cmplx b7  = load_b_from_file(6);
    const tens_cmplx b10 = load_b_from_file(9);
    std::cout << " done\n";

    std::cout << " Integrating the loop on the q grid..." << std::flush;
    const vec_double pi_renorm = compute_hvp_pi(b1, b7, b10, q_grid, k_grid, z_grid);
    std::cout << " done\n";

    const std::string out_file = "hvp.dat";
    std::ofstream fs(out_file, std::ofstream::out | std::ofstream::trunc);
    fs << "#p_sq Pi_renormalised\n";
    for (unsigned i = 0; i < q_grid.size(); ++i)
      fs << q_grid[i] << " " << pi_renorm[i] << "\n";
    fs.close();

    std::cout << " HVP written to " << out_file
              << "  (renormalised at p² = " << q_grid.front() << ")\n";
  }
} // namespace hvp
