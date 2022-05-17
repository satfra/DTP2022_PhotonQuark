#pragma once

#include <complex>
#include <chrono>
#include "omp.h"

#include "Utils.hh"
#include "types.hh"
#include "QuadratureIntegral.hh"
#include "fileIO.hh"

#include "parameters.hh"
#include "Kernels_G.hh"
#include "Kernels_K.hh"
#include "momentumtransform.hh"
#include "basistransform.hh"
#include "maris_tandy.hh"
#include "WTI.hh"

using Integrator1d = qIntegral<LegendrePolynomial<parameters::numerical::y_steps>>;
using Integrator2d = qIntegral2d<LegendrePolynomial<parameters::numerical::k_steps>, LegendrePolynomial<parameters::numerical::z_steps>>;

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

mat_cmplx average_array_full(const tens_cmplx &a, const unsigned z_0)
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

template<typename Quark>
void a_iteration_step(const tens_cmplx &b,
    const ijtens2_double &K_prime, const vec_double &z_grid, const vec_double &k_grid, tens_cmplx &a,
    const Integrator2d& qint2d, const Quark& quark)
{
  using namespace parameters::numerical;
  #pragma omp parallel for collapse(3)
  for (unsigned i = 0; i < n_structs; ++i)
    for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
      for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx) 
      {
        // Initialize the a's with the inhomogeneous term
        a[i][k_idx][z_idx] = quark.z2() * a0(i);

        for (unsigned j = 0; j < n_structs; ++j)
        {
          if (K::isZeroIndex(i, j))
            continue;

          // The function to integrate
          auto f = [&](const double &k_prime_sq_log, const double &z_prime) {
            const double k_prime_sq = std::exp(k_prime_sq_log);
            lInterpolator2d interpolate2d_b(k_grid, z_grid, b[j]);
            const auto b_j = interpolate2d_b(k_prime_sq_log, z_prime);

            lInterpolator2d interpolate2d_K(k_grid, z_grid, K_prime[i][k_idx][z_idx][j]);
            const auto K_prime_ij = interpolate2d_K(k_prime_sq_log, z_prime);

            return K_prime_ij * b_j * std::sqrt(1. - powr<2>(z_prime)) * powr<2>(k_prime_sq);
          };

          // Evaluate the integral
          const std::complex<double> integral = qint2d(f, 
              k_grid[0], k_grid[k_steps - 1], 
              z_grid[0], z_grid[z_steps - 1]);

          // Add this to the a's
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
    const vec_double &z_grid, const vec_double &k_grid, ijtens2_double &K_prime, const Quark& quark, const bool use_PauliVillars)
{
  using namespace parameters::numerical;
  #pragma omp parallel for collapse(3)
  for (unsigned i = 0; i < n_structs; ++i)
    for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
      for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx)
        for (unsigned j = 0; j < n_structs; ++j)
        {
          if (K::isZeroIndex(i, j))
            continue;
          for (unsigned k_prime_idx = 0; k_prime_idx < k_steps; ++k_prime_idx)
            for (unsigned z_prime_idx = 0; z_prime_idx < z_steps; ++z_prime_idx)
            {
              const double &z = z_grid[z_idx];
              const double k_sq = std::exp(k_grid[k_idx]);
              const double &z_prime = z_grid[z_prime_idx];
              const double k_prime_sq = std::exp(k_grid[k_prime_idx]);

              K_prime[i][k_idx][z_idx][j][k_prime_idx][z_prime_idx] = 0.;

              auto f = [&](const double &y)
              {
                const double l_sq = momentumtransform::l2(k_sq, k_prime_sq, z, z_prime, y);

                const double gl = use_PauliVillars ? pauli_villars_g(l_sq, quark): maris_tandy_g(l_sq, quark);

                K k_kernel(k_sq, k_prime_sq, z, z_prime, y, q_sq);
                return gl * k_kernel.get(i, j);
              };

              // Evaluate the integral
              const double integral = qint1d(f, y_grid[0], y_grid[y_steps - 1]);

              // Add this to the a's
              K_prime[i][k_idx][z_idx][j][k_prime_idx][z_prime_idx] = integral;
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

  const vec_cmplx temp0(z_steps, 0.0);
  const mat_cmplx temp1(k_steps, temp0);

  const vec_double temp0_d(z_steps, 0.0);
  const mat_double temp1_d(k_steps, temp0_d);
  const tens_double temp2_d(n_structs, temp1_d);
  const tens2_double temp3_d(z_steps, temp2_d);
  const jtens2_double temp4_d(k_steps, temp3_d);

  // Do some Legendre Magic
  Integrator1d qint1d;
  Integrator2d qint2d;

  const Quark quark;

  // prepare output files
  emptyIdxFile<12>("fg_file", "#q_sq i k_sq z Re(fg) Im(fg)");
  emptyIdxFile<12>("fg_z0_file", "#q_sq i k_sq Re(fg) Im(fg)");
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
    tens_cmplx a(n_structs, temp1);
    tens_cmplx b(n_structs, temp1);

    // Precalculate the K kernel
    std::cout << " Calculating K'_ij..." << std::flush;
    ijtens2_double K_prime(n_structs, temp4_d);
    precalculate_K_kernel(y_grid, qint1d, q_sq, z_grid, k_grid, K_prime, quark, use_PauliVillars);
    std::cout << " done\n";

    // Initialize a with bare vertex
    a_initialize(a, quark);

    std::cout << "  Starting iteration...\n";
    double current_acc = 1.0;
    unsigned current_step = 0;
    while (max_steps > current_step++ && current_acc > target_acc)
    {
      debug_out("\n    Started a step...\n", debug);
      
      // copy for checking the convergence
      const auto a_old = average_array_full(a, z_0);

      debug_out("    Calculating b_i...", debug);
      b_iteration_step(a, q_sq, z_grid, k_grid, b, quark);
      debug_out(" done\n", debug);

      debug_out("    Calculating a_i...", debug);
      a_iteration_step(b, K_prime, z_grid, k_grid, a, qint2d, quark);
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
    // transform to g,f (almost in place!)
    transform_a_to_fg(a, q_sq, k_grid, z_grid);
    const auto& fg = a;
    // save to the prepared file
    saveToFile_withGrids<n_structs>(fg, "fg_file.dat", q_sq, k_grid, z_grid);
    // save to the prepared file
    const auto fg_z0 = average_array_z0(fg, z_0);
    saveToFile_withGrids<n_structs>(fg_z0, "fg_z0_file.dat", q_sq, k_grid);
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
    tens_cmplx w(3, temp1);
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
    saveToFile_withGrids<3>(w, "w_file.dat", q_sq, k_grid, z_grid);
    const auto w_z0 = average_array_z0(w, z_0);
    saveToFile_withGrids<3>(w_z0, "w_z0_file.dat", q_sq, k_grid);
  }

  auto end_time = steady_clock::now();
  std::cout << "\nProgram finished after " << duration_cast<milliseconds>(end_time - start_time).count()/1000.<< "s\n";
}
