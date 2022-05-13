#pragma once

#include <vector>
#include <complex>
#include "Utils.hh"
#include "coeff.hh"
#include "Kernels_G.hh"
#include "Kernels_K.hh"
#include "quark_model_functions.hh"
#include "momentumtransform.hh"
#include "fileIO.hh"
#include "omp.h"

void iterate_a_and_b(const vec_double &q_grid, const vec_double &z_grid, const vec_double &k_grid,
    const vec_double &y_grid) {
  // Model value for Z_2. Must be updated once we use a real quark.
  constexpr double z_2 = 0.97;

  constexpr unsigned int n_structs = 12;
  const unsigned int k_steps = k_grid.size();
  const unsigned int z_steps = z_grid.size();
  const unsigned int q_steps = q_grid.size();
  const unsigned int y_steps = y_grid.size();
  const vec_cmplx temp0(z_steps, 0.0);
  const mat_cmplx temp1(k_steps, temp0);
  const tens_cmplx temp2(n_structs, temp1);
  qtens_cmplx a(q_steps, temp2);
  qtens_cmplx b(q_steps, temp2);

  const vec_double temp0_d(z_steps, 0.0);
  const mat_double temp1_d(k_steps, temp0_d);
  const tens_double temp2_d(n_structs, temp1_d);
  const tens2_double temp3_d(z_steps, temp2_d);
  const jtens2_double temp4_d(k_steps, temp3_d);

  constexpr double target_acc = 1e-3;
  constexpr unsigned int max_steps = 30;

  // Do some Legendre Magic
  constexpr unsigned order_z_prime = 4;
  constexpr unsigned order_k_prime = 10;
  constexpr unsigned order_y = 4;
  qIntegral2d<LegendrePolynomial<order_k_prime>, LegendrePolynomial<order_z_prime>> qint2d;
  qIntegral<LegendrePolynomial<order_y>> qint1d;

  // loop over q
  for (unsigned int q_iter = 0; q_iter < q_steps; q_iter++) {
    const double& q_sq = q_grid[q_iter];
    double current_acc = 1.0;
    unsigned current_step = 0;

    ijtens2_double K_prime(n_structs, temp4_d);
    #pragma omp parallel for collapse(2)
    for (unsigned i = 0; i < n_structs; ++i)
    {
      for (unsigned j = 0; j < n_structs; ++j)
      {
        if (K::isZeroIndex(i,j))
          continue;
        for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx) 
        {
          const double& k_sq = k_grid[k_idx];
          for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx)
          {
            const double& z = z_grid[z_idx];
            for (unsigned k_prime_idx = 0; k_prime_idx < k_steps; ++k_prime_idx)
            {
              const double& k_prime_sq = k_grid[k_prime_idx];
              for (unsigned z_prime_idx = 0; z_prime_idx < z_steps; ++z_prime_idx)
              {
                const double& z_prime = z_grid[z_prime_idx];

                // The function to integrate
                auto f = [&](const double &y) {
                  const double l_sq = momentumtransform::l2(k_sq, k_prime_sq, z, z_prime, y);
                  const double gl = maris_tandy_g(l_sq, 1.8, 0.72);
                  K k_kernel(k_sq, k_prime_sq, z, z_prime, y, q_sq);
                  return gl * k_kernel.get(i, j);
                };

                // Evaluate the integral
                const double integral = qint1d(f, y_grid[0], y_grid[y_steps - 1]);

                // Add this to the a's
                K_prime[i][k_idx][z_idx][j][k_prime_idx][z_prime_idx] += integral;
              }
            }
          }
        }
        std::cout << "K_" << i << j << " = " << K_prime[i][0][0][j][0][0] << "\n";
      }
    }
    std::cout << "Calculated K'_ij...\n";

    // Initialize the a's with the inhomogeneous term
    for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
      for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx)
        for (unsigned i = 0; i < n_structs; ++i)
          a[q_iter][i][k_idx][z_idx] = z_2 * a0(i);

    std::cout << "Initialized a_i...\n";

    while (max_steps > current_step && current_acc > target_acc) {
      std::cout << "Started a step...\n";

      const std::complex<double> a_old = a[q_iter][0][0][0];

      #pragma omp parallel for collapse(2)
      for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx) {
        for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx) {
          const double& k_sq = k_grid[k_idx];
          const double& z = z_grid[z_idx];

          // Evaluate Gij
          G g_kernel(k_sq, z, q_sq);
          for (unsigned i = 0; i < n_structs; ++i) {
            // Initialize the b's to 0
            b[q_iter][i][k_idx][z_idx] = 0.0;

            for (unsigned j = 0; j < n_structs; ++j) {
              // Add stuff to the b's
              b[q_iter][i][k_idx][z_idx] += g_kernel.get(i, j) * a[q_iter][j][k_idx][z_idx];
              //std::cout << "g_kernel[" << i << "," << j << "] = " << g_kernel.get(i, j) << "\n";
            }
            //std::cout << "b[" << i << "]≠" << b[q_iter][i][k_idx][z_idx] << "\n";
          }
        }
      }

      std::cout << "  Calculated b_i...\n";

      #pragma omp parallel for// collapse(2)
      for (unsigned i = 0; i < n_structs; ++i) {
        for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx) {
          for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx) {
            //const double& k_sq = k_grid[k_idx];
            //const double& z = z_grid[z_idx];
            // Initialize the a's with the inhomogeneous term
            a[q_iter][i][k_idx][z_idx] = z_2 * a0(i);

            for (unsigned j = 0; j < n_structs; ++j) {

              if (K::isZeroIndex(i,j))
                continue;
              // The function to integrate
              auto f = [=](const double &k_prime_sq, const double &z_prime) {
                lInterpolator2d interpolate2d_b(k_grid, z_grid, b[q_iter][j]);
                const auto b_j = interpolate2d_b(k_prime_sq, z_prime);

                lInterpolator2d interpolate2d_K(k_grid, z_grid, K_prime[i][k_idx][z_idx][j]);
                const auto K_prime_ij = interpolate2d_K(k_prime_sq, z_prime);

                return 2.0 * M_PI * K_prime_ij * b_j;
              };

              // Evaluate the integral
              const std::complex<double> integral = qint2d(f, 
                  k_grid[0], k_grid[k_steps - 1], 
                  z_grid[0], z_grid[z_steps - 1]);

              //std::cout << "integral = " << integral << " at k2=" << k_sq << " z=" << z << " i=" << i << " j=" << j << "\n";

              // Add this to the a's
              a[q_iter][i][k_idx][z_idx] += integral;
            }
          }
        }
      }

      std::cout << "  Calculated a_i...\n";
  
      // TODO better convergence test
      const std::complex<double> a_new = a[q_iter][0][0][0];
      current_acc = abs(a_new - a_old) / abs(a_new + a_old);
      ++current_step;
      if (current_step == max_steps) {
        std::cout << "Maximum iterations reached!" << std::endl;
      }

      std::cout << "current_step = " << current_step << "\n"
        << "a[q_iter][0][0][0] = " << a[q_iter][0][0][3] << "\n"
        << "b[q_iter][0][0][0] = " << b[q_iter][0][0][3] << "\n";
    }
  }
  saveToFile(a, "file_a", "q_sq i k_sq z");
  saveToFile(b, "file_b", "q_sq i k_sq z");
}
