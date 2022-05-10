#include <vector>
#include <complex>
#include "coeff.hh"

typedef std::vector<std::complex<double>> vec_cmplx;
typedef std::vector<vec_cmplx> mat_cmplx;
typedef std::vector<mat_cmplx> tens_cmplx;
typedef std::vector<double> vec_double;

void iterate_a_and_b(const vec_double &q_grid, const vec_double &z_grid, const vec_double &k_grid, const vec_double & y_grid)
{
    constexpr double z_2 = 0.97; // Model value for Z_2. Must be updated once we use a real quark.

    constexpr unsigned int n_structs = 12;
    const unsigned int k_steps = k_grid.size();
    const unsigned int z_steps = z_grid.size();
    const unsigned int q_steps = q_grid.size();
    const unsigned int y_steps = y_grid.size();
    const vec_cmplx temp2(z_steps, 0.0);
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
    vec_double k_kernel(k_dimension, 0.0);

    // TODO: Parallelize me!
    // ... but I might not even be needed ...
    for (int super_idx = 0; super_idx < k_dimension; ++super_idx) {
        const unsigned int i = k_dimension % n_structs;
        const unsigned int j = (k_dimension / n_structs) % n_structs;
        const unsigned int k_idx = (k_dimension / (n_structs * n_structs)) % k_steps;
        const unsigned int k_prime_idx = (k_dimension / (n_structs * n_structs * k_steps)) % k_steps;
        const unsigned int z_idx = (k_dimension / (n_structs * n_structs * k_steps * k_steps)) % z_steps;
        const unsigned int z_prime_idx = (k_dimension / (n_structs * n_structs * k_steps * k_steps * z_steps)) % z_steps;
        const unsigned int y_idx = (k_dimension / (n_structs * n_structs * k_steps * k_steps * z_steps * z_steps)) % y_steps;
        const unsigned int q_idx = (k_dimension / (n_structs * n_structs * k_steps * k_steps * z_steps * z_steps * y_steps)) % q_steps;

        const double k = k_grid[k_idx];
        const double k_prime = k_grid[k_prime_idx];
        const double z = z_grid[z_idx];
        const double z_prime = z_grid[z_prime_idx];
        const double y = y_grid[y_idx];
        const double q = q_grid[q_idx];

        // TODO: Transform the variables to the ones used in the K kernel functions

        k_kernel[super_idx] = 0.0; // TODO: Call the K kernel functions here
    }

    while (current_acc > target_acc && max_steps > current_step) {
        ++current_step;
        // TODO: Loop over Q
        // Iterate the a_i's
        for (int i = 0; i < n_structs; ++i) {
            for (int z_idx = 0; z_idx < z_steps; ++z_idx) {
                const double z = z_grid[z_idx];
                for (int k_idx = 0; k_idx < k_steps; ++k_idx) {
                    const double k = k_grid[k_idx];
                    a[i][k_idx][z_idx] = z_2 * a0(i);
                    for (int j = 0; j < n_structs; ++j) {
                        a[i][k_idx][z_idx] += Integrate[
                                g((k-kprime)*(k-kprime)) * K(i, j, k*k, kprime*kprime, z, zprime, y, q*q) * b[j](kprime*kprime, zprime, q*q)
                                ,kprime, zprime
                                ] // TODO: Finish this up
                    }
                }
            }
        }
    }
}
