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
    const vec_cmplx temp0(z_steps, 0.0);
    const mat_cmplx temp2(k_steps, temp1);
    const mat_cmplx temp3(k_steps, temp0);
    const tens_cmplx temp4(n_structs, temp2);
    const tens_cmplx temp5(n_structs, temp3);
    qtens_cmplx a(q_steps, temp4);
    qtens_cmplx b(q_steps, temp5);

    constexpr double target_acc = 1e-5;
    constexpr unsigned int max_steps = 25;
    double current_acc_a;
    double current_acc_b;
    bool a_converged = false;
    bool b_converged = false;

    // Do some Legendre Magic
    constexpr unsigned order_z_prime = 8;
    constexpr unsigned order_k_prime = 20;
    constexpr unsigned order_y = 6;
    qIntegral3d<LegendrePolynomial<order_k_prime>, LegendrePolynomial<order_z_prime>, LegendrePolynomial<order_y>> qint;

    // If only one of a and b is converged, this counter will start counting the amount of iterations this is the
    // case for. If this happens for too many iterations, it might be a good idea to check, if the other functions
    // will change again. Thus, they converged flags will be turned off again, once this reaches a certain value.
    unsigned int suspicion_counter = 0;
    // loop over q
    for (unsigned int q_iter = 0; q_iter < q_steps; q_iter++) {
        const double q_sq = q_grid[q_iter];
        constexpr unsigned int max_sus_counter = 10;
        // loop over i
        for (unsigned int i = 0; i < n_structs; ++i) {
            unsigned int current_step = 0;
            while ((!a_converged || !b_converged) && max_steps > current_step) {
                if (a_converged != b_converged) {
                    ++suspicion_counter;
                } else {
                    suspicion_counter = 0;
                }

                const std::complex<double> a_old = a[0][0][0][0];
                const std::complex<double> b_old = b[0][0][0][0];

                // loop over k
                for (unsigned int k_idx = 0; k_idx < k_steps; ++k_idx) {
                    const double k_sq = k_grid[k_idx];
                    // loop over z
                    for (unsigned int z_idx = 0; z_idx < z_steps; ++z_idx) {
                        const double z = z_grid[z_idx];

                        if (!a_converged || suspicion_counter >= max_sus_counter) {
                            suspicion_counter = 0;
                            /*
                 _                              _        _            _            _           _                 _          _          _            _
                / /\                           /\ \     /\ \         /\ \         /\ \        / /\              /\ \       /\ \       /\ \         /\ \     _
               / /  \                          \ \ \    \_\ \       /  \ \       /  \ \      / /  \             \_\ \      \ \ \     /  \ \       /  \ \   /\_\
              / / /\ \                         /\ \_\   /\__ \     / /\ \ \     / /\ \ \    / / /\ \            /\__ \     /\ \_\   / /\ \ \     / /\ \ \_/ / /
             / / /\ \ \        ____           / /\/_/  / /_ \ \   / / /\ \_\   / / /\ \_\  / / /\ \ \          / /_ \ \   / /\/_/  / / /\ \ \   / / /\ \___/ /
            / / /  \ \ \     /\____/\        / / /    / / /\ \ \ / /_/_ \/_/  / / /_/ / / / / /  \ \ \        / / /\ \ \ / / /    / / /  \ \_\ / / /  \/____/
           / / /___/ /\ \    \/____\/       / / /    / / /  \/_// /____/\    / / /__\/ / / / /___/ /\ \      / / /  \/_// / /    / / /   / / // / /    / / /
          / / /_____/ /\ \                 / / /    / / /      / /\____\/   / / /_____/ / / /_____/ /\ \    / / /      / / /    / / /   / / // / /    / / /
         / /_________/\ \ \            ___/ / /__  / / /      / / /______  / / /\ \ \  / /_________/\ \ \  / / /   ___/ / /__  / / /___/ / // / /    / / /
        / / /_       __\ \_\          /\__\/_/___\/_/ /      / / /_______\/ / /  \ \ \/ / /_       __\ \_\/_/ /   /\__\/_/___\/ / /____\/ // / /    / / /
        \_\___\     /____/_/          \/_________/\_\/       \/__________/\/_/    \_\/\_\___\     /____/_/\_\/    \/_________/\/_________/ \/_/     \/_/
                             */
                            // Initialize the a's with the inhomogeneous term
                            a[q_iter][i][k_idx][z_idx] = z_2 * a0(i);

                            // loop over j
                            for (unsigned int j = 0; j < n_structs; ++j) {
                                // The function to integrate
                                auto f = [&](const double &k_prime_sq, const double &z_prime, const double &y) {
                                    const double l_sq = momentumtransform::l2(k_sq, k_prime_sq, z, z_prime, y);
                                    const double gl = maris_tandy_g(l_sq, 1.8, 0.72);
                                    K k_kernel(k_sq, k_prime_sq, z, z_prime, y, q_sq);
                                    // TODO: Make this work with the interface provided by Jonas
                                    lInterpolator2d interpolate2d(k_grid, z_grid, b[q_iter][j]);
                                    const auto b_j = interpolate2d(k_prime_sq, z_prime);

                                    return gl * k_kernel.get(i, j) * b_j;
                                };

                                // Evaluate the integral
                                const std::complex<double> integral = qint(f, k_grid[0], k_grid[k_steps - 1], z_grid[0],
                                                                           z_grid[z_steps - 1], y_grid[0],
                                                                           y_grid[y_steps - 1]);

                                // Add this to the a's
                                a[q_iter][i][k_idx][z_idx] += integral;
                            }
                        } else if (suspicion_counter < max_sus_counter) ++suspicion_counter;

                        if (!b_converged || suspicion_counter >= max_sus_counter) {
                            suspicion_counter = 0;
                            /*
               _                           _        _            _            _           _                 _          _          _            _
              / /\                        /\ \     /\ \         /\ \         /\ \        / /\              /\ \       /\ \       /\ \         /\ \     _
             / /  \                       \ \ \    \_\ \       /  \ \       /  \ \      / /  \             \_\ \      \ \ \     /  \ \       /  \ \   /\_\
            / / /\ \                      /\ \_\   /\__ \     / /\ \ \     / /\ \ \    / / /\ \            /\__ \     /\ \_\   / /\ \ \     / /\ \ \_/ / /
           / / /\ \ \     ____           / /\/_/  / /_ \ \   / / /\ \_\   / / /\ \_\  / / /\ \ \          / /_ \ \   / /\/_/  / / /\ \ \   / / /\ \___/ /
          / / /\ \_\ \  /\____/\        / / /    / / /\ \ \ / /_/_ \/_/  / / /_/ / / / / /  \ \ \        / / /\ \ \ / / /    / / /  \ \_\ / / /  \/____/
         / / /\ \ \___\ \/____\/       / / /    / / /  \/_// /____/\    / / /__\/ / / / /___/ /\ \      / / /  \/_// / /    / / /   / / // / /    / / /
        / / /  \ \ \__/               / / /    / / /      / /\____\/   / / /_____/ / / /_____/ /\ \    / / /      / / /    / / /   / / // / /    / / /
       / / /____\_\ \             ___/ / /__  / / /      / / /______  / / /\ \ \  / /_________/\ \ \  / / /   ___/ / /__  / / /___/ / // / /    / / /
      / / /__________\           /\__\/_/___\/_/ /      / / /_______\/ / /  \ \ \/ / /_       __\ \_\/_/ /   /\__\/_/___\/ / /____\/ // / /    / / /
      \/_____________/           \/_________/\_\/       \/__________/\/_/    \_\/\_\___\     /____/_/\_\/    \/_________/\/_________/ \/_/     \/_/

                             */
                            // Initialize the b's to 0
                            b[q_iter][i][k_idx][z_idx] = 0.0;

                            // Evaluate Gij
                            G g_kernel(k_sq, z, q_sq);
                            for (unsigned int j = 0; j < n_structs; ++j) {
                                // Add stuff to the b's
                                b[q_iter][i][k_idx][z_idx] += g_kernel.get(i, j) * a[q_iter][j][k_idx][z_idx];
                            }
                        } else if (suspicion_counter < max_sus_counter) ++suspicion_counter;
                    }
                }
                // TODO: Do this estimation and switching on and off for all i and q separately
                // TODO: Do a better error estimation
                const std::complex<double> a_new = a[0][0][0][0];
                const std::complex<double> b_new = b[0][0][0][0];
                current_acc_a = abs(a_new - a_old) / abs(a_new + a_old);
                current_acc_b = abs(b_new - b_old) / abs(b_new + b_old);
                if (current_acc_a < target_acc) a_converged = true;
                if (current_acc_b < target_acc) b_converged = true;
                ++current_step;
                if (current_step == max_steps) std::cout << "Maximum iterations reached!" << std::endl;
            }
        }
    }
    saveToFile(a, "file_a");
    saveToFile(b, "file_b");
}
