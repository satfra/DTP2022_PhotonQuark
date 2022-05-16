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
#include "parameters.hh"
#include "basistransform.hh"
#include "omp.h"

using Integrator1d = qIntegral<LegendrePolynomial<parameters::numerical::y_steps>>;
using Integrator2d = qIntegral2d<LegendrePolynomial<parameters::numerical::k_steps>, LegendrePolynomial<parameters::numerical::z_steps>>;

double update_accuracy(const unsigned int z_0, const qtens_cmplx &a, unsigned int q_iter,
    double current_acc, const std::vector<mat_cmplx> &a_old) {
  using namespace parameters::numerical;
  for (unsigned i = 0; i < n_structs; ++i) {
    for (unsigned k_idx = 0; k_idx < parameters::numerical::k_steps; ++k_idx) {
      const auto diff = a[q_iter][i][k_idx][z_0] - a_old[i][k_idx][z_0];
      const auto sum = a[q_iter][i][k_idx][z_0] + a_old[i][k_idx][z_0];
      current_acc = std::max(current_acc, abs(diff) / abs(sum));
    }
  }
  return current_acc;
}

void a_initialize(qtens_cmplx &a, unsigned int q_iter) {
  using namespace parameters::numerical;
#pragma omp parallel for collapse(2)
  for (unsigned i = 0; i < n_structs; ++i) {
    for (unsigned k_idx = 0; k_idx < parameters::numerical::k_steps; ++k_idx) {
      for (unsigned z_idx = 0; z_idx < parameters::numerical::z_steps; ++z_idx) {
        // Initialize the b's to 0
        a[q_iter][i][k_idx][z_idx] = parameters::physical::z_2 * a0(i);
      }
    }
  }
}

void a_iteration_step(const qtens_cmplx &b, unsigned int q_iter,
    const ijtens2_double &K_prime, const vec_double &z_grid, const vec_double &k_grid, qtens_cmplx &a,
    const Integrator2d& qint2d) {
  using namespace parameters::numerical;
#pragma omp parallel for collapse(2)
  for (unsigned i = 0; i < n_structs; ++i) {
    for (unsigned k_idx = 0; k_idx < parameters::numerical::k_steps; ++k_idx) {
      for (unsigned z_idx = 0; z_idx < parameters::numerical::z_steps; ++z_idx) {
        // Initialize the a's with the inhomogeneous term
        a[q_iter][i][k_idx][z_idx] = parameters::physical::z_2 * a0(i);

        for (unsigned j = 0; j < n_structs; ++j) {
          if (K::isZeroIndex(i, j))
            continue;

          // The function to integrate
          auto f = [&](const double &k_prime_sq_log, const double &z_prime) {
            const double k_prime_sq = exp(k_prime_sq_log);
            lInterpolator2d interpolate2d_b(k_grid, z_grid, b[q_iter][j]);
            const auto b_j = interpolate2d_b(k_prime_sq_log, z_prime);

            lInterpolator2d interpolate2d_K(k_grid, z_grid, K_prime[i][k_idx][z_idx][j]);
            const auto K_prime_ij = interpolate2d_K(k_prime_sq_log, z_prime);

            return K_prime_ij * b_j * sqrt(1. - powr<2>(z_prime)) * k_prime_sq * k_prime_sq;
          };

          // Evaluate the integral
          const std::complex<double> integral = qint2d(f, k_grid[0],
              k_grid[parameters::numerical::k_steps - 1], z_grid[0],
              z_grid[parameters::numerical::z_steps - 1]);

          // Add this to the a's
          a[q_iter][i][k_idx][z_idx] += integral * 2.0 * M_PI * int_factors;
        }
      }
    }
  }
}

void b_iteration_step(const qtens_cmplx &a, unsigned int q_iter, const double &q_sq,
    const vec_double &z_grid, const vec_double &k_grid, qtens_cmplx &b) {
  using namespace parameters::numerical;
#pragma omp parallel for collapse(2)
  for (unsigned k_idx = 0; k_idx < parameters::numerical::k_steps; ++k_idx) {
    for (unsigned z_idx = 0; z_idx < parameters::numerical::z_steps; ++z_idx) {
      const double k_sq = std::exp(k_grid[k_idx]);
      const double &z = z_grid[z_idx];

      // Evaluate Gij
      G g_kernel(k_sq, z, q_sq);
      for (unsigned i = 0; i < n_structs; ++i) {
        // Initialize the b's to 0
        b[q_iter][i][k_idx][z_idx] = 0.0;

        // Add stuff to the b's
        for (unsigned j = 0; j < n_structs; ++j) {
          b[q_iter][i][k_idx][z_idx] += g_kernel.get(i, j) * a[q_iter][j][k_idx][z_idx];
        }
      }
    }
  }
}

void precalculate_K_kernel(const vec_double &y_grid,
    const Integrator1d &qint1d, const double &q_sq,
    const vec_double &z_grid, const vec_double &k_grid, ijtens2_double &K_prime) {
  using namespace parameters::numerical;
#pragma omp parallel for collapse(2)
  for (unsigned i = 0; i < n_structs; ++i) {
    for (unsigned j = 0; j < n_structs; ++j) {
      if (K::isZeroIndex(i, j))
        continue;

      for (unsigned k_idx = 0; k_idx < parameters::numerical::k_steps; ++k_idx) {
        const double k_sq = std::exp(k_grid[k_idx]);

        for (unsigned z_idx = 0; z_idx < parameters::numerical::z_steps; ++z_idx) {
          const double &z = z_grid[z_idx];

          for (unsigned k_prime_idx = 0; k_prime_idx < parameters::numerical::k_steps; ++k_prime_idx) {
            const double &k_prime_sq = exp(k_grid[k_prime_idx]);

            for (unsigned z_prime_idx = 0; z_prime_idx < parameters::numerical::z_steps; ++z_prime_idx) {
              const double &z_prime = z_grid[z_prime_idx];

              K_prime[i][k_idx][z_idx][j][k_prime_idx][z_prime_idx] = 0.;

              auto f = [&](const double &y) {
                const double l_sq = momentumtransform::l2(k_sq, k_prime_sq, z, z_prime, y);

                const double gl = parameters::physical::pauliVillars ? pauli_villars_g(l_sq): maris_tandy_g(l_sq);
                K k_kernel(k_sq, k_prime_sq, z, z_prime, y, q_sq);
                return gl * k_kernel.get(i, j);
              };

              // Evaluate the integral
              const double integral = qint1d(f, y_grid[0], y_grid[parameters::numerical::y_steps - 1]);

              // Add this to the a's
              K_prime[i][k_idx][z_idx][j][k_prime_idx][z_prime_idx] = integral;
            }
          }
        }
      }
    }
  }
}

void transform_a_to_fg(qtens_cmplx &a, 
    const vec_double &q_grid, const vec_double &k_grid, const vec_double &z_grid) {
  using namespace parameters::numerical;
  using namespace basistransform;

  for (unsigned int q_iter = 0; q_iter < q_steps; q_iter++) {
    for (unsigned k_idx = 0; k_idx < parameters::numerical::k_steps; ++k_idx) {
      for (unsigned z_idx = 0; z_idx < parameters::numerical::z_steps; ++z_idx) {
        vec_cmplx a_copy(n_structs);
        for(unsigned i = 0; i < n_structs; ++i)
          a_copy[i] = a[q_iter][i][k_idx][z_idx];

        const double Q = std::sqrt(q_grid[q_iter]);
        const double k_sq = std::exp(k_grid[k_idx]);
        const double k = std::sqrt(k_sq);
        const double& z = z_grid[z_idx];
        const double s = std::sqrt(1. - powr<2>(z));

        a[q_iter][0][k_idx][z_idx] = f1(Q, s, z, k, a_copy);
        a[q_iter][1][k_idx][z_idx] = f2(Q, s, z, k, a_copy);
        a[q_iter][2][k_idx][z_idx] = f3(Q, s, z, k, a_copy);
        a[q_iter][3][k_idx][z_idx] = f4(Q, s, z, k, a_copy);
        a[q_iter][4][k_idx][z_idx] = f5(Q, s, z, k, a_copy);
        a[q_iter][5][k_idx][z_idx] = f6(Q, s, z, k, a_copy);
        a[q_iter][6][k_idx][z_idx] = f7(Q, s, z, k, a_copy);
        a[q_iter][7][k_idx][z_idx] = f8(Q, s, z, k, a_copy);

        a[q_iter][8][k_idx][z_idx] = g1(Q, s, z, k, a_copy);
        a[q_iter][9][k_idx][z_idx] = g2(Q, s, z, k, a_copy);
        a[q_iter][10][k_idx][z_idx] = g3(Q, s, z, k, a_copy);
        a[q_iter][11][k_idx][z_idx] = g4(Q, s, z, k, a_copy);
      }
    }
  }
}

void iterate_a_and_b(const vec_double &q_grid, const vec_double &z_grid, const vec_double &k_grid, const vec_double &y_grid)
{
  using namespace parameters::numerical;
  const unsigned z_0 = z_grid.size() / 2;

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

  // Do some Legendre Magic
  Integrator1d qint1d;
  Integrator2d qint2d;

  // loop over q
  for (unsigned int q_iter = 0; q_iter < q_steps; q_iter++) {
    const double &q_sq = q_grid[q_iter];
    double current_acc = 1.0;
    unsigned current_step = 0;

    ijtens2_double K_prime(n_structs, temp4_d);

    // Precalculate the K kernel
    precalculate_K_kernel(y_grid, qint1d, q_sq, z_grid, k_grid, K_prime);
    std::cout << "Calculated K'_ij...\n";

    // Initialize a with bare vertex
    a_initialize(a, q_iter);

    std::cout << "Starting iteration...\n";
    while (max_steps > current_step && current_acc > target_acc) {
      std::cout << "\nStarted a step...\n";

      // consider NOT copying... RAM lol
      const auto a_old = a[q_iter];

      b_iteration_step(a, q_iter, q_sq, z_grid, k_grid, b);

      std::cout << "  Calculated b_i...\n";

      a_iteration_step(b, q_iter, K_prime, z_grid, k_grid, a, qint2d);

      std::cout << "  Calculated a_i...\n";

      current_acc = 0.;
      current_acc = update_accuracy(z_0, a, q_iter, current_acc, a_old);

      ++current_step;
      if (current_step == max_steps) {
        std::cout << "Maximum iterations reached!" << std::endl;
      }
      std::cout << "  current_step = " << current_step << "\n" << "  current_acc = " << current_acc << "\n";
      if (current_acc < target_acc) {
        std::cout << "----> Converged!" << std::endl;
      }
    }
  }

  // transform to g,f (almost in place!)
  transform_a_to_fg(a, q_grid, k_grid, z_grid);
  const auto& fg = a;
 
  // TODO check WTI

  saveToFile_withGrids(fg, "fg_file.dat", "#q_sq i k_sq z Re(fg) Im(fg)", 
      q_grid, k_grid, z_grid);
}
