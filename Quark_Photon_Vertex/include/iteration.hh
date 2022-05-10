#include <vector>
#include <complex>
#include "coeff.hh"
#include "Kernels_G.hh"
#include "Kernels_K.hh"
#include "quark_model_functions.hh"

typedef std::vector<std::complex<double>> vec_cmplx;
typedef std::vector<vec_cmplx> mat_cmplx;
typedef std::vector<mat_cmplx> tens_cmplx;
typedef std::vector<double> vec_double;

// TODO: Write this function and put it in an appropriate place.
double interpolate(const vec_double &grid, const double x) {
    return 1.0;
}

void iterate_a_and_b(const vec_double &q_grid, const vec_double &z_grid, const vec_double &k_grid, const vec_double & y_grid)
{
    constexpr double z_2 = 0.97; // Model value for Z_2. Must be updated once we use a real quark.

    constexpr unsigned int n_structs = 12;
    const unsigned int k_steps = k_grid.size();
    const unsigned int z_steps = z_grid.size();
    const unsigned int q_steps = q_grid.size();
    const unsigned int y_steps = y_grid.size();
    const vec_cmplx temp2(z_steps, 1.0);
    const mat_cmplx temp(k_steps, temp2);
    tens_cmplx a(n_structs, temp);
    tens_cmplx b(n_structs, temp);

    constexpr double target_acc = 1e-5;
    constexpr unsigned int max_steps = 100;
    double current_acc = 1.0;
    unsigned int current_step = 0;

    // Generate the K kernel, so we don't have to do a lot of function calls
    // Inside the iterations
    const unsigned int k_dimension = n_structs * n_structs * k_steps * k_steps * z_steps * z_steps * y_steps * q_steps;

    // Do some Legendre Magic
    constexpr unsigned order_z_prime = 10;
    constexpr unsigned order_k_prime = 30;
    qIntegral2d<LegendrePolynomial<order_k_prime>, LegendrePolynomial<order_z_prime>> qint;

    while (current_acc > target_acc && max_steps > current_step) {
        ++current_step;
        // For now this will be for a fixed value of q_sq. TODO: Generalize this
        const double q_sq = 0.0;

        // loop over i
        for (unsigned int i = 0; i < n_structs; ++i) {
            // loop over k
            for (unsigned int k_idx = 0; k_idx < k_steps; ++k_idx) {
                const double k_sq = k_grid[k_idx];
                // loop over z
                for (unsigned int z_idx = 0; z_idx < z_steps; ++z_idx) {
                    const double z = z_grid[z_idx];

                    // Initialize them with the inhomogeneous term
                    a[i][k_idx][z_idx] = z_2 * a0(i);

                    // loop over j
                    for (unsigned int j = 0; j < n_structs; ++j) {
                        // The function to integrate
                        auto f = [const &](const double& k_prime_sq, const double& z_prime){
                            const double l_sq = k_sq + k_prime_sq - 2.0 * std::sqrt(k_prime_sq * k_sq);
                            const double gl = maris_tandy_g(l_sq);
                            const K k_kernel = new K(k_sq, k_prime_sq, z, z_prime, y, q_sq);
                            const double b_j = interpolate2d(b[j], k_prime_sq, z_prime);

                            return gl * k_kernel.get(i, j) * b_j;
                        };

                        // Evaluate the integral
                        const std::complex<double> integral = qint(f, z_grid[0], z_grid[z_steps-1], k_grid[0],
                                                                                        k_grid[k_steps-1]);

                        // Add this to the a's
                        a[i][k_idx][z_idx] += integral;
                    }
                }
            }
        }
    }
}
