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

void iterate_a_and_b(const vec_double &q_grid, const vec_double &z_grid, const vec_double &k_grid,
                     const vec_double &y_grid) {
    // Model value for Z_2. Must be updated once we use a real quark.
    constexpr double z_2 = 0.97;

    constexpr unsigned int n_structs = 12;
    const unsigned int k_steps = k_grid.size();
    const unsigned int z_steps = z_grid.size();
    const unsigned int q_steps = q_grid.size();
    const unsigned int y_steps = y_grid.size();
    const vec_cmplx temp1(z_steps, 1.0);
    const vec_cmplx temp0(z_steps, 1.0);
    const mat_cmplx temp2(k_steps, temp1);
    const mat_cmplx temp3(k_steps, temp0);
    const tens_cmplx temp4(n_structs, temp2);
    const tens_cmplx temp5(n_structs, temp3);
    qtens_cmplx a(q_steps, temp4);
    qtens_cmplx b(q_steps, temp5);

    constexpr double target_acc = 1e-3;
    constexpr unsigned int max_steps = 30;

    // Do some Legendre Magic
    constexpr unsigned order_z_prime = 8;
    constexpr unsigned order_k_prime = 20;
    constexpr unsigned order_y = 6;
    qIntegral3d<LegendrePolynomial<order_k_prime>, LegendrePolynomial<order_z_prime>, LegendrePolynomial<order_y>> qint;

    // loop over q
    for (unsigned int q_iter = 0; q_iter < q_steps; q_iter++) {
        const double q_sq = q_grid[q_iter];
        double current_acc = 1.0;
        unsigned int current_step = 0;
        while (max_steps > current_step && current_acc > target_acc) {

            const std::complex<double> a_old = a[0][0][0][0];

            // loop over i
            for (unsigned int i = 0; i < n_structs; ++i) {
                // loop over k
                for (unsigned int k_idx = 0; k_idx < k_steps; ++k_idx) {
                    const double k_sq = k_grid[k_idx];
                    // loop over z
                    for (unsigned int z_idx = 0; z_idx < z_steps; ++z_idx) {
                        const double z = z_grid[z_idx];

                        // Initialize the b's to 0
                        b[q_iter][i][k_idx][z_idx] = 0.0;

                        // Evaluate Gij
                        G g_kernel(k_sq, z, q_sq);
                        for (unsigned int j = 0; j < n_structs; ++j) {
                            // Add stuff to the b's
                            b[q_iter][i][k_idx][z_idx] += g_kernel.get(i, j) * a[q_iter][j][k_idx][z_idx];
                        }
                    }
                }
            }

            // loop over i
            for (unsigned int i = 0; i < n_structs; ++i) {
                // loop over k
                for (unsigned int k_idx = 0; k_idx < k_steps; ++k_idx) {
                    const double k_sq = k_grid[k_idx];
                    // loop over z
                    for (unsigned int z_idx = 0; z_idx < z_steps; ++z_idx) {
                        const double z = z_grid[z_idx];


                        // Initialize the a's with the inhomogeneous term
                        a[q_iter][i][k_idx][z_idx] = z_2 * a0(i);

                        // loop over j
                        for (unsigned int j = 0; j < n_structs; ++j) {
                            // The function to integrate
                            auto f = [&](const double &k_prime_sq, const double &z_prime, const double &y) {
                                const double l_sq = momentumtransform::l2(k_sq, k_prime_sq, z, z_prime, y);
                                const double gl = maris_tandy_g(l_sq, 1.8, 0.72);
                                K k_kernel(k_sq, k_prime_sq, z, z_prime, y, q_sq);
                                lInterpolator2d interpolate2d(k_grid, z_grid, b[q_iter][j]);
                                const auto b_j = interpolate2d(k_prime_sq, z_prime);

                                return 2.0 * M_PI * gl * k_kernel.get(i, j) * b_j;
                            };

                            // Evaluate the integral
                            const std::complex<double> integral = qint(f, k_grid[0], k_grid[k_steps - 1], z_grid[0],
                                                                       z_grid[z_steps - 1], y_grid[0],
                                                                       y_grid[y_steps - 1]);

                            // Add this to the a's
                            a[q_iter][i][k_idx][z_idx] += integral;
                        }
                    }
                }
            }

            const std::complex<double> a_new = a[0][0][0][0];
            current_acc = abs(a_new - a_old) / abs(a_new + a_old);
            ++current_step;
            if (current_step == max_steps) {
                std::cout << "Maximum iterations reached!" << current_step << " " << max_steps << std::endl;
            }
        }
    }
    saveToFile(a, "file_a");
    saveToFile(b, "file_b");
}
