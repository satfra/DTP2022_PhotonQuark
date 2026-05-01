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
#include "LinearInterpolate.hh"

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

    tens_cmplx out(q_steps,
        mat_cmplx(k_steps, vec_cmplx(z_steps, 0.0)));

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

  // z dimension uses Gauss–Chebyshev type 2 so the √(1−z²) Jacobian of
  // the 4D loop measure is absorbed into the quadrature weights.
  using HvpIntegrator = qIntegral2d<
      LegendrePolynomial<parameters::numerical::k_steps>,
      ChebyshevPolynomial2<parameters::numerical::z_steps>>;

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
  inline vec_double compute_hvp_pi(
      const tens_cmplx& b1, const tens_cmplx& b7, const tens_cmplx& b10,
      const vec_double& q_grid, const vec_double& k_grid, const vec_double& z_grid)
  {
    using namespace parameters::numerical;

    HvpIntegrator qint2d;
    vec_double pi_values(q_steps, 0.0);

    constexpr double pref_b1  = 4.0 * M_SQRT2;
    constexpr double pref_b7  = 4.0;
    constexpr double pref_b10 = 4.0;

    constexpr unsigned q_ren = 0;
    const mat_cmplx& b1_ren  = b1[q_ren];
    const mat_cmplx& b7_ren  = b7[q_ren];
    const mat_cmplx& b10_ren = b10[q_ren];

    #pragma omp parallel for
    for (unsigned q_iter = 0; q_iter < q_steps; ++q_iter)
    {
      lInterpolator2d ip_b1     (k_grid, z_grid, b1[q_iter]);
      lInterpolator2d ip_b7     (k_grid, z_grid, b7[q_iter]);
      lInterpolator2d ip_b10    (k_grid, z_grid, b10[q_iter]);
      lInterpolator2d ip_b1_ren (k_grid, z_grid, b1_ren);
      lInterpolator2d ip_b7_ren (k_grid, z_grid, b7_ren);
      lInterpolator2d ip_b10_ren(k_grid, z_grid, b10_ren);

      auto integrand = [&](const double& log_k_sq, const double& z)
      {
        const std::complex<double> trace_p =
            pref_b1  * ip_b1 (log_k_sq, z) +
            pref_b7  * ip_b7 (log_k_sq, z) +
            pref_b10 * ip_b10(log_k_sq, z);

        const std::complex<double> trace_ren =
            pref_b1  * ip_b1_ren (log_k_sq, z) +
            pref_b7  * ip_b7_ren (log_k_sq, z) +
            pref_b10 * ip_b10_ren(log_k_sq, z);

        (void)z;  // z enters only through the trace; no explicit Jacobian.
        const double k_sq = std::exp(log_k_sq);
        // 4D Euclidean measure with the log-k² Jacobian:
        //   d⁴k/(2π)⁴ · f  =  int_factors · 2π · √(1−z²) · k² · d(k²) dz · f
        //                  =  int_factors · 2π · √(1−z²) · (k²)² · d(log k²) dz · f
        // The √(1−z²) factor is the Chebyshev type 2 weight (see
        // HvpIntegrator above), so it must NOT be multiplied in here.
        const double measure = int_factors * 2. * M_PI * powr<2>(k_sq);

        return measure * (trace_p - trace_ren);
      };

      const std::complex<double> integral = qint2d(integrand,
          k_grid.front(), k_grid.back(),
          -1.0, 1.0);

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
