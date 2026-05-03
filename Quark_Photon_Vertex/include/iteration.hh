#pragma once

#include <complex>
#include <chrono>
#include "omp.h"

#include "Utils.hh"
#include "types.hh"
#include "QuadratureIntegral.hh"
#include "ChebyshevPolynomial2.hh"
#include "LinearInterpolate.hh"
#include "spline.h"
#include "fileIO.hh"

#include "parameters.hh"
#include "Kernels_G.hh"
#include "Kernels_K.hh"
#include "momentumtransform.hh"
#include "basistransform.hh"
#include "maris_tandy.hh"
#include "WTI.hh"

using Integrator1d = qIntegral<LegendrePolynomial<parameters::numerical::y_steps>>;
// z dimension uses Gauss–Chebyshev type 2: the √(1−z²) Jacobian of the
// 4D loop measure is absorbed into the quadrature weights, so it must
// NOT appear in the integrand and the z-bounds must be (-1, 1).
using Integrator2d = qIntegral2d<LegendrePolynomial<parameters::numerical::k_steps>, ChebyshevPolynomial2<parameters::numerical::z_steps>>;

double update_accuracy(const unsigned z_0, const tens_cmplx &a, const mat_cmplx &a_old)
{
  double current_acc = 0.;

  using namespace parameters::numerical;
  for (unsigned i = 0; i < n_structs; ++i)
    for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
    {
      std::complex<double> current_a = 0.;
      for (unsigned z_idx = 0; z_idx < a[i][k_idx].size(); ++z_idx)
        current_a += std::abs(a[i][k_idx][z_idx]) / double(a[i][k_idx].size());
      const auto diff = current_a - a_old[i][k_idx];
      const auto sum = current_a + a_old[i][k_idx];
      if(!isEqual(std::abs(sum), 0.))
        current_acc = std::max(current_acc, abs(diff) / std::abs(sum));
    }
  return current_acc;
}

mat_cmplx average_array_full(const tens_cmplx &a)
{
  using namespace parameters::numerical;
  mat_cmplx a_z0(a.size(), vec_cmplx(a[0].size(), 0.0));

  for (unsigned i = 0; i < a.size(); ++i)
    for (unsigned k_idx = 0; k_idx < a[i].size(); ++k_idx)
      for (unsigned z_idx = 0; z_idx < a[i][k_idx].size(); ++z_idx)
        a_z0[i][k_idx] += std::abs(a[i][k_idx][z_idx]) / double(a[i][k_idx].size());
  return a_z0;
}

mat_cmplx average_array_z0(const tens_cmplx &a, const unsigned z_0)
{
  using namespace parameters::numerical;
  mat_cmplx a_z0(a.size(), vec_cmplx(a[0].size(), 0.0));

  for (unsigned i = 0; i < a.size(); ++i)
    for (unsigned k_idx = 0; k_idx < a[i].size(); ++k_idx)
      a_z0[i][k_idx] = 0.5 * (a[i][k_idx][z_0-1] + a[i][k_idx][z_0]);
  return a_z0;
}

double a0 (const unsigned& i)
{
  if(i == 0)
    return std::sqrt(2.);
  else if(i == 6 || i == 9)
    return 1.;
  return 0.;
}

// Initializes a_i with the bare vertex (inhomogeneous term a_i^0).
// See Eq. (48) in the project description.
template<typename Quark>
void a_initialize(tens_cmplx &a, const Quark& quark)
{
  using namespace parameters::numerical;
  #pragma omp parallel for collapse(2)
  for (unsigned i = 0; i < n_structs; ++i)
    for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
      for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx)
        a[i][k_idx][z_idx] = quark.z2() * a0(i);
}

// Cubic spline of a complex-valued grid function in one real coordinate.
// Used to interpolate b[j](log k'², z'_idx) along log(k'²) at each fixed
// Cheb2 z'-grid index. Replaces the linear-in-log(k'²) interpolation that
// was the dominant interp-error contributor at the IR k² corner where
// WTI2 max-error sits.
struct ComplexSpline {
  tk::spline re_, im_;

  void set_points(const std::vector<double>& x, const vec_cmplx& y) {
    std::vector<double> y_re(y.size()), y_im(y.size());
    for (std::size_t i = 0; i < y.size(); ++i) {
      y_re[i] = y[i].real();
      y_im[i] = y[i].imag();
    }
    re_.set_points(x, y_re);
    im_.set_points(x, y_im);
  }

  std::complex<double> operator()(double x) const {
    return std::complex<double>(re_(x), im_(x));
  }
};

// Build cubic splines of b[j][:, zp_idx] over log(k'²) for each (j, zp_idx).
// Done once per BSE iteration after b_iteration_step refreshes b. The 384
// splines (n_structs × z_steps at default sizes) are cheap: each is an
// O(k_steps) tridiagonal solve. Storage ≈ 4·k_steps doubles per spline,
// total ~1.6 MB at the default grid.
std::vector<std::vector<ComplexSpline>> build_b_splines(
    const tens_cmplx& b, const vec_double& k_grid)
{
  using namespace parameters::numerical;
  std::vector<std::vector<ComplexSpline>> spline_b(
      n_structs, std::vector<ComplexSpline>(z_steps));
  for (unsigned j = 0; j < n_structs; ++j) {
    for (unsigned zp_idx = 0; zp_idx < z_steps; ++zp_idx) {
      vec_cmplx b_slice(k_steps);
      for (unsigned k = 0; k < k_steps; ++k)
        b_slice[k] = b[j][k][zp_idx];
      spline_b[j][zp_idx].set_points(k_grid, b_slice);
    }
  }
  return spline_b;
}

template<typename Quark>
void a_iteration_step(const tens_cmplx & /*b*/,
    const ijtens2_double &K_prime, const vec_double &z_grid, const vec_double &k_grid, tens_cmplx &a,
    const Integrator2d& /*qint2d*/, const Quark& quark,
    const std::vector<std::vector<ComplexSpline>>& spline_b)
{
  using namespace parameters::numerical;

  // Cheb2 z'-quadrature: nodes coincide with z_grid by construction (the
  // Integrator2d uses ChebyshevPolynomial2<z_steps> on (−1, 1) and so does
  // Simulation.cpp when building z_grid), so no z-interpolation of b/K'
  // is needed — we just sample the stored grid value at zp_idx.
  static const ChebyshevPolynomial2<z_steps> cp_z;
  static const auto z_weights = cp_z.weights();

  // Legendre quadrature in log(k'²), mapped from [−1,1] to
  // [k_grid[0], k_grid[k_steps-1]] via dk = (b−a)/2, k_mid = (a+b)/2.
  static const LegendrePolynomial<k_steps> lp_k;
  static const auto& k_unit_zeros = lp_k.zeroes();
  static const auto& k_unit_weights = lp_k.weights();

  const double k_a = k_grid[0];
  const double k_b = k_grid[k_steps - 1];
  const double dk = 0.5 * (k_b - k_a);
  const double k_mid = 0.5 * (k_a + k_b);

  std::vector<double> k_quad_log(k_steps), k_quad_powr2(k_steps);
  for (unsigned q = 0; q < k_steps; ++q) {
    k_quad_log[q] = k_mid + dk * k_unit_zeros[q];
    const double k_prime_sq = std::exp(k_quad_log[q]);
    k_quad_powr2[q] = k_prime_sq * k_prime_sq;
  }

  #pragma omp parallel for collapse(3)
  for (unsigned i = 0; i < n_structs; ++i)
    for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
      for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx)
      {
        // Inhomogeneous term: a^0_i (Eq. 48), scaled by Z_2.
        a[i][k_idx][z_idx] = quark.z2() * a0(i);

        for (unsigned j = 0; j < n_structs; ++j)
        {
          if (K::isZeroIndex(i, j))
            continue;

          // K_prime stays linear-in-log(k'²) for now: per-(i,k_idx,z_idx,j)
          // spline rebuilds would cost ~3M tridiagonal solves per BSE step
          // and dominate the run. lInterpolator2d returns the exact stored
          // value when z' is at a grid point (which it is here), so the
          // 2D interp degenerates to 1D linear-in-log(k'²) cleanly.
          const lInterpolator2d interp_K(k_grid, z_grid, K_prime.slice2d(i, k_idx, z_idx, j));

          std::complex<double> integral{0., 0.};
          for (unsigned zp_idx = 0; zp_idx < z_steps; ++zp_idx)
          {
            const double z_prime = z_grid[zp_idx];
            const auto& spline_bj = spline_b[j][zp_idx];

            std::complex<double> partial{0., 0.};
            for (unsigned q = 0; q < k_steps; ++q)
            {
              const double K_at = interp_K.unchecked(k_quad_log[q], z_prime);
              const auto b_at = spline_bj(k_quad_log[q]);
              partial += k_unit_weights[q] * K_at * b_at * k_quad_powr2[q];
            }
            partial *= dk;
            // Cheb2 weight already absorbs √(1−z'²); no extra factor here.
            integral += z_weights[zp_idx] * partial;
          }

          a[i][k_idx][z_idx] += integral * 2.0 * M_PI * int_factors;
        }
      }
}

template<typename Quark>
void b_iteration_step(const tens_cmplx &a, const double &q_sq,
    const vec_double &z_grid, const vec_double &k_grid, tens_cmplx &b, const Quark& quark)
{
  using namespace parameters::numerical;
  #pragma omp parallel for collapse(2)
  for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
    for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx) 
    {
      const double k_sq = std::exp(k_grid[k_idx]);
      const double &z = z_grid[z_idx];

      // Evaluate Gij
      const G<Quark> g_kernel(k_sq, z, q_sq, quark);
      for (unsigned i = 0; i < n_structs; ++i)
      {
        // Initialize the b's to 0
        b[i][k_idx][z_idx] = 0.0;

        // Add stuff to the b's
        for (unsigned j = 0; j < n_structs; ++j)
          b[i][k_idx][z_idx] += g_kernel.get(i, j) * a[j][k_idx][z_idx];
      }
    }
}

template<typename Quark>
void precalculate_K_kernel(const vec_double &y_grid,
    const Integrator1d &qint1d, const double &q_sq,
    const vec_double &z_grid, const vec_double &k_grid,
    const vec_double &k_sq_table, const vec_double &k_table,
    const vec_double &s_z_table,
    ijtens2_double &K_prime, const Quark& quark, const bool use_PauliVillars)
{
  using namespace parameters::numerical;
  // collapse(5) widens the OpenMP iteration space by 144× over the previous
  // collapse(3); pruned (i,j) pairs (≈40% via K::isZeroIndex) become cheap
  // continues but threads stay saturated. z_prime_idx remains serial inside
  // each task so the inner write to K_prime walks a contiguous row.
  #pragma omp parallel for collapse(5)
  for (unsigned i = 0; i < n_structs; ++i)
    for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
      for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx)
        for (unsigned j = 0; j < n_structs; ++j)
          for (unsigned k_prime_idx = 0; k_prime_idx < k_steps; ++k_prime_idx)
          {
            if (K::isZeroIndex(i, j))
              continue;

            const double k_sq = k_sq_table[k_idx];
            const double k_prime_sq = k_sq_table[k_prime_idx];
            const double k_v = k_table[k_idx];
            const double k_prime_v = k_table[k_prime_idx];
            const double k_kp_sqrt = k_v * k_prime_v;
            const double& z = z_grid[z_idx];
            const double s_z = s_z_table[z_idx];

            for (unsigned z_prime_idx = 0; z_prime_idx < z_steps; ++z_prime_idx)
            {
              const double& z_prime = z_grid[z_prime_idx];
              const double s_z_prime = s_z_table[z_prime_idx];

              auto f = [&](const double &y)
              {
                // Inlined momentumtransform::l2 with the precomputed
                // sqrt(k²·k'²) and sin terms; matches the old formula bit-for-
                // bit modulo associativity (≤1 ULP under -ffast-math).
                const double l_sq = k_sq + k_prime_sq
                    - 2. * k_kp_sqrt * (z * z_prime + y * s_z * s_z_prime);

                const double gl = use_PauliVillars
                    ? pauli_villars_g(l_sq, quark)
                    : maris_tandy_g(l_sq, quark);

                K k_kernel(k_sq, k_prime_sq, z, z_prime, y,
                           k_kp_sqrt, s_z, s_z_prime, k_v, k_prime_v);
                return gl * k_kernel.get(i, j);
              };

              K_prime(i, k_idx, z_idx, j, k_prime_idx, z_prime_idx) =
                  qint1d(f, y_grid[0], y_grid[y_steps - 1]);
            }
          }
}

void transform_a_to_fg(tens_cmplx &a, const double& q_sq, const vec_double &k_grid, const vec_double &z_grid)
{
  using namespace parameters::numerical;
  using namespace basistransform;

  for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
    for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx)
    {
      vec_cmplx a_copy(n_structs);
      for(unsigned i = 0; i < n_structs; ++i)
        a_copy[i] = a[i][k_idx][z_idx];

      const double Q = std::sqrt(q_sq);
      const double k_sq = std::exp(k_grid[k_idx]);
      const double k = std::sqrt(k_sq);
      const double& z = z_grid[z_idx];
      const double s = std::sqrt(1. - powr<2>(z));

      a[0][k_idx][z_idx] = f1(Q, s, z, k, a_copy);
      a[1][k_idx][z_idx] = f2(Q, s, z, k, a_copy);
      a[2][k_idx][z_idx] = f3(Q, s, z, k, a_copy);
      a[3][k_idx][z_idx] = f4(Q, s, z, k, a_copy);
      a[4][k_idx][z_idx] = f5(Q, s, z, k, a_copy);
      a[5][k_idx][z_idx] = f6(Q, s, z, k, a_copy);
      a[6][k_idx][z_idx] = f7(Q, s, z, k, a_copy);
      a[7][k_idx][z_idx] = f8(Q, s, z, k, a_copy);

      a[8][k_idx][z_idx] = g1(Q, s, z, k, a_copy);
      a[9][k_idx][z_idx] = g2(Q, s, z, k, a_copy);
      a[10][k_idx][z_idx] = g3(Q, s, z, k, a_copy);
      a[11][k_idx][z_idx] = g4(Q, s, z, k, a_copy);
    }
}

template<typename Quark>
void iterate_a_and_b(const vec_double &q_grid, const vec_double &z_grid, const vec_double &k_grid, const vec_double &y_grid, const bool use_PauliVillars, const bool debug)
{
  using namespace std::chrono;
  const auto start_time = steady_clock::now();

  using namespace parameters::numerical;
  const unsigned z_0 = z_grid.size() / 2;

  // Do some Legendre Magic
  Integrator1d qint1d;
  Integrator2d qint2d;

  const Quark quark;

  // Hoist the k² = exp(k_grid), √k² and √(1-z²) tables out of the inner
  // loops — they only depend on the (fixed) k/z grids, so the K-kernel
  // precalculation goes from ~10⁹ exp/sqrt evaluations per q-point to a few
  // hundred table lookups.
  vec_double k_sq_table(k_steps);
  vec_double k_table(k_steps);
  for (unsigned i = 0; i < k_steps; ++i) {
    k_sq_table[i] = std::exp(k_grid[i]);
    k_table[i] = std::sqrt(k_sq_table[i]);
  }
  vec_double s_z_table(z_steps);
  for (unsigned i = 0; i < z_steps; ++i)
    s_z_table[i] = std::sqrt(1. - z_grid[i] * z_grid[i]);

  // prepare output files
  emptyIdxFile<12>("fg_file", "#q_sq i k_sq z Re(fg) Im(fg)");
  emptyIdxFile<12>("fg_z0_file", "#q_sq i k_sq Re(fg) Im(fg)");
  emptyIdxFile<12>("b_file", "#q_sq i k_sq z Re(b) Im(b)");
  emptyIdxFile<3>("w_file", "#q_sq i k_sq z Re(w) Im(w)");
  emptyIdxFile<3>("w_z0_file", "#q_sq i k_sq Re(w) Im(w)");

  // loop over q
  for (unsigned q_iter = 0; q_iter < q_steps; q_iter++)
  {
    const auto iter_start_time = steady_clock::now();

    const double &q_sq = q_grid[q_iter];
    std::cout << "\n_____________________________________\n\n"
      << "Calculation for q^2 = " << q_sq << "\n";

    // define new a, b
    tens_cmplx a(n_structs, k_steps, z_steps);
    tens_cmplx b(n_structs, k_steps, z_steps);

    // Precalculate the K kernel
    std::cout << " Calculating K'_ij..." << std::flush;
    // K_prime is sparse over (i,j): only 26 of 144 (i,j) entries are non-zero
    // (8+4 block decoupling + within-block sparsity from K::isZeroIndex).
    // SparseTensor6 only allocates the 26 slots, ~5.5× memory reduction.
    ijtens2_double K_prime(k_steps, z_steps, k_steps, z_steps,
                           [](unsigned i, unsigned j) { return K::isZeroIndex(i, j); });
    precalculate_K_kernel(y_grid, qint1d, q_sq, z_grid, k_grid,
                          k_sq_table, k_table, s_z_table,
                          K_prime, quark, use_PauliVillars);
    std::cout << " done\n";

    // Initialize a with bare vertex
    a_initialize(a, quark);

    // Main self-consistent iteration: a → b → a until convergence.
    // Implements the coupled system in Eq. (43) of the project description.
    std::cout << "  Starting iteration...\n";
    double current_acc = 1.0;
    unsigned current_step = 0;
    while (max_steps > current_step++ && current_acc > target_acc)
    {
      debug_out("\n    Started a step...\n", debug);
      
      // copy for checking the convergence
      const auto a_old = average_array_full(a);

      debug_out("    Calculating b_i...", debug);
      b_iteration_step(a, q_sq, z_grid, k_grid, b, quark);
      debug_out(" done\n", debug);

      // Refresh the cubic splines of b[j](log k'²) at each z'-grid point.
      // Used inside a_iteration_step to interpolate b at the Legendre
      // quadrature nodes in k' — replaces the linear log-k' interpolation
      // that limited WTI2 accuracy at the IR k corner.
      const auto spline_b = build_b_splines(b, k_grid);

      debug_out("    Calculating a_i...", debug);
      a_iteration_step(b, K_prime, z_grid, k_grid, a, qint2d, quark, spline_b);
      debug_out(" done\n", debug);

      // check the convergence
      current_acc = update_accuracy(z_0, a, a_old);
      debug_out("    current_step = " + std::to_string(current_step) + "\n    current_acc = " + std::to_string(current_acc) + "\n", debug);
    }
    if (current_acc < target_acc)
      std::cout << "  Converged after " << current_step << " steps.\n";
    else
      std::cout << "  ! Did not converge !\n";

    std::cout << " Saving results..." << std::flush;
    // resync b with the converged a, since the last loop step recomputes a
    // from the previous-iteration b. With this call b = G_kernel · a holds for
    // the converged a, which is what the HVP loop expects (b ≡ S·Γ·S).
    b_iteration_step(a, q_sq, z_grid, k_grid, b, quark);
    saveToFile_withGrids<n_structs>(b, "b_file", q_sq, k_grid, z_grid);
    // transform to g,f (almost in place!)
    transform_a_to_fg(a, q_sq, k_grid, z_grid);
    const auto& fg = a;
    // save to the prepared file
    saveToFile_withGrids<n_structs>(fg, "fg_file", q_sq, k_grid, z_grid);
    // save to the prepared file
    const auto fg_z0 = average_array_z0(fg, z_0);
    saveToFile_withGrids<n_structs>(fg_z0, "fg_z0_file", q_sq, k_grid);
    std::cout << "  done\n";

    const auto iter_end_time = steady_clock::now();
    std::cout << "Calculation finished after " << duration_cast<milliseconds>(iter_end_time - iter_start_time).count()/1000.<< "s\n";
  }

  std::cout << "\n_____________________________________\n\n"
    << "\nCalculating the WTIs...\n";

  // check WTI
  for (unsigned q_iter = 0; q_iter < q_steps; q_iter++)
  {
    const double &q_sq = q_grid[q_iter];
    tens_cmplx w(3, k_steps, z_steps);
    for (unsigned k_idx = 0; k_idx < parameters::numerical::k_steps; ++k_idx)
    {
      for (unsigned z_idx = 0; z_idx < parameters::numerical::z_steps; ++z_idx)
      {
        const double q_sq = q_grid[q_iter];
        const double Q = std::sqrt(q_sq);
        const double k_sq = std::exp(k_grid[k_idx]);
        const double k = std::sqrt(k_sq);
        const double& z = z_grid[z_idx];

        double kplus2 = k_sq+ q_sq/4. + k*Q*z;
        double kminus2 = k_sq + q_sq/4. - k*Q*z;

        w[0][k_idx][z_idx] = Sigma_A<Quark>(kminus2,kplus2,quark);
        w[1][k_idx][z_idx] = Delta_A<Quark>(kminus2,kplus2,quark);
        w[2][k_idx][z_idx] = Delta_B<Quark>(kminus2,kplus2,quark);
      }
    }
    saveToFile_withGrids<3>(w, "w_file", q_sq, k_grid, z_grid);
    const auto w_z0 = average_array_z0(w, z_0);
    saveToFile_withGrids<3>(w_z0, "w_z0_file", q_sq, k_grid);
  }

  auto end_time = steady_clock::now();
  std::cout << "\nProgram finished after " << duration_cast<milliseconds>(end_time - start_time).count()/1000.<< "s\n";
}
