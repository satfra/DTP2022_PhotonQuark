#include <vector>
#include <complex>
#include "coeff.hh"

typedef std::vector<std::complex<double>> vec_cmplx;
typedef std::vector<double> vec_double;

void iterate_a_and_b(const vec_double &q_grid, const vec_double &z_grid, const vec_double &k_grid)
{
    constexpr double z_2 = 0.97; // Model value for Z_2. Must be updated once we use a real quark.

    constexpr unsigned int n_structs = 12;
    const vec_cmplx temp(k_grid.size(), 0.0);
    std::vector<vec_cmplx> a(n_structs, temp);
    std::vector<vec_cmplx> b(n_structs, temp);

    constexpr double target_acc = 1e-5;
    constexpr unsigned int max_steps = 100;
    double current_acc = 1.0;
    unsigned int current_step = 0;

    while (current_acc > target_acc && max_steps > current_step) {
        ++current_step;
        // TODO: Loop over Q
        // Iterate the a_i's
        for (int i = 0; i < n_structs; ++i) {
            for (int z_idx = 0; z_idx < ; ++z_idx) {
                const double z = z_grid[z_idx];
                for (int k_idx = 0; k < k_grid.size(); ++k) {
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
