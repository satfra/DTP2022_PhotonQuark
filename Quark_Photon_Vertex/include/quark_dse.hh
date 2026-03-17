#pragma once

#include <algorithm>
#include <cmath>

#include "Utils.hh"
#include "types.hh"
#include "maris_tandy.hh"
#include <LegendrePolynomials.hh>

// ---------------------------------------------------------------------------
// Angular integration matrix: pre-computes the combined momentum+angle
// integrands for A and B, reused across DSE iterations to save time.
// ---------------------------------------------------------------------------
inline tens_double init_brl_angular_matrix(double mu,
    const vec_double &dse_absci_q, const vec_double &dse_weights_q,
    const vec_double &dse_absci_ang, const vec_double &dse_weights_ang)
{
    tens_double temp(2);
    mat_double temp_mat(parameters::numerical::quark_dse_steps_q);
    vec_double temp_vec(parameters::numerical::quark_dse_steps_q + 1);
    for (unsigned int i = 0; i < parameters::numerical::quark_dse_steps_q; ++i)
        temp_mat[i] = temp_vec;
    temp[0] = temp_mat;
    temp[1] = temp_mat;

#pragma omp parallel for
    for (unsigned int qidx = 0; qidx < parameters::numerical::quark_dse_steps_q; ++qidx) {
        const double q = exp(dse_absci_q[qidx] / 2.0);
        for (unsigned int pidx = 0; pidx < parameters::numerical::quark_dse_steps_q + 1; ++pidx) {
            double p;
            if (pidx == parameters::numerical::quark_dse_steps_q) {
                p = mu;
            } else {
                p = exp(dse_absci_q[pidx] / 2.0);
            }
            double sa = 0.0;
            double sb = 0.0;

            for (unsigned int ang_i = 0; ang_i < parameters::numerical::quark_dse_steps_z; ++ang_i) {
                const double psi = dse_absci_ang[ang_i];
                const double z = cos(psi);
                const double sin2psi = 1.0 - z * z; // sin²(ψ) = angular_function(ψ)
                const double pq = p * q * z;
                const double kk = p * p + q * q - 2.0 * p * q * z;

                const double mt = maris_tandy_alpha(kk);

                sa += dse_weights_ang[ang_i] *
                        dse_weights_q[qidx] * sin2psi *
                        (2.0 * (M_1_PI * M_1_PI) / (3.0 * p * p)) *
                        q * q * q * q * (1.0 / kk) *
                        (mt * (
                                -2*p*p + (2*(p*p - pq)*(p*p - pq))/kk + 3*pq
                            )
                        );

                sb += dse_weights_ang[ang_i] *
                        dse_weights_q[qidx] * sin2psi *
                        (2.0 * M_1_PI * M_1_PI) * (1.0 / 3.0) *
                        q * q * q * q * (1.0 / kk) *
                        (mt * 3.0);
            }
            temp[0][qidx][pidx] = sa;
            temp[1][qidx][pidx] = sb;
        }
    }
    return temp;
}

inline double brl_integrate_coupled_a(int pidx, const vec_double &a_values,
                               const vec_double &b_values, const tens_double &angular_matrix,
                               const vec_double &dse_absci_q)
{
    /*
     * Integration needed to iterate A. The angular integration is already
     * baked into angular_matrix; only the momentum sum remains.
     */
    double s_a = 0.0;

#pragma omp parallel for reduction(+:s_a)
    for (unsigned int j = 0; j < parameters::numerical::quark_dse_steps_q; ++j) {
        const double q = exp(dse_absci_q[j] / 2.0);
        s_a += angular_matrix[0][j][pidx] * a_values[j] /
               ((q * a_values[j] * q * a_values[j]) +
                (b_values[j] * b_values[j]));
    }
    return s_a;
}

inline double brl_integrate_coupled_b(int pidx, const vec_double &a_values,
        const vec_double &b_values, const tens_double &angular_matrix,
        const vec_double &absci_x)
{
    /*
     * Integration needed to iterate B. The angular integration is already
     * baked into angular_matrix; only the momentum sum remains.
     */
    double s_b = 0.0;

#pragma omp parallel for reduction(+:s_b)
    for (unsigned int j = 0; j < parameters::numerical::quark_dse_steps_q; ++j) {
        const double q = exp(absci_x[j] / 2.0);
        s_b += angular_matrix[1][j][pidx] * b_values[j] /
                ((q * a_values[j] * q * a_values[j]) +
                (b_values[j] * b_values[j]));
    }
    return s_b;
}

inline mat_double quark_iterate_dressing_functions(double a0, double b0, double mc, double mu)
{
    /*
     * Iterates the quark dressing functions A and B until convergence.
     * Returns [a_values, b_values, dse_absci_q, renorm] as a 2d vector.
     *
     * Gauss-Legendre grids are built using the existing LegendrePolynomial
     * infrastructure (replaces the C-style gauleg() from Numerical Recipes).
     * NOTE: do NOT parallelize the grid construction — see commit fd56693.
     */
    vec_double a_values(parameters::numerical::quark_dse_steps_q, a0);
    vec_double b_values(parameters::numerical::quark_dse_steps_q, b0);

    // Build GL momentum grid on [2·log(λ_IR), 2·log(λ_UV)]
    LegendrePolynomial<parameters::numerical::quark_dse_steps_q> lp_q;
    const double q_lo = 2.0 * std::log(parameters::physical::lambda_IR);
    const double q_hi = 2.0 * std::log(parameters::physical::lambda_UV);
    const vec_double dse_absci_q = linearMapTo(lp_q.zeroes(), -1., 1., q_lo, q_hi);
    vec_double dse_weights_q = lp_q.weights();
    for (auto& w : dse_weights_q) w *= (q_hi - q_lo) / 2.;

    // Build GL angular grid on [0, π]
    LegendrePolynomial<parameters::numerical::quark_dse_steps_z> lp_z;
    const vec_double dse_absci_z = linearMapTo(lp_z.zeroes(), -1., 1., 0., M_PI);
    vec_double dse_weights_z = lp_z.weights();
    for (auto& w : dse_weights_z) w *= M_PI / 2.;

#pragma omp parallel for
    for (unsigned int i = 0; i < parameters::numerical::quark_dse_steps_q; ++i) {
        double p = exp(dse_absci_q[i] / 2.0);
        if (p > 1.0) {
            a_values[i] = 1;
        } else {
            b_values[i] = 0.8;
        }
    }

    /*
     * Pre-compute the angular matrix once; it's the same for every iteration.
     */
    tens_double angular_matrix = init_brl_angular_matrix(
        mu, dse_absci_q, dse_weights_q, dse_absci_z, dse_weights_z);

    double current_acc, a_start, b_start, a_end, b_end, sigma_a, sigma_b;
    vec_double renorm = {1.0, 1.0};
    vec_double new_a(parameters::numerical::quark_dse_steps_q);
    vec_double new_b(parameters::numerical::quark_dse_steps_q);
    unsigned int k = 0;

    do {
        a_start = a_values[parameters::numerical::quark_dse_steps_q - 1];
        b_start = b_values[parameters::numerical::quark_dse_steps_q - 1];

#pragma omp parallel for
        for (unsigned int pidx = 0; pidx < parameters::numerical::quark_dse_steps_q; ++pidx) {
            new_a[pidx] = renorm[0] * (1.0 + renorm[0] *
                    brl_integrate_coupled_a(pidx, a_values, b_values,
                            angular_matrix, dse_absci_q));
            new_b[pidx] = renorm[0] * (mc * renorm[1] + renorm[0] *
                 brl_integrate_coupled_b(pidx, a_values, b_values,
                            angular_matrix, dse_absci_q));
        }

        sigma_a = brl_integrate_coupled_a(parameters::numerical::quark_dse_steps_q,
                a_values, b_values, angular_matrix, dse_absci_q);
        renorm[0] = 1.0 / (1.0 + renorm[0] * sigma_a);
        sigma_b = brl_integrate_coupled_b(parameters::numerical::quark_dse_steps_q,
                a_values, b_values, angular_matrix, dse_absci_q);
        renorm[1] = 1.0 / renorm[0] - renorm[0] * sigma_b / mc;

        for (unsigned int l = 0; l < parameters::numerical::quark_dse_steps_q; ++l) {
            a_values[l] = new_a[l];
            b_values[l] = new_b[l];
        }

        a_end = a_values[parameters::numerical::quark_dse_steps_q - 1];
        b_end = b_values[parameters::numerical::quark_dse_steps_q - 1];

        current_acc = std::max(fabs((b_end - b_start) / (b_end + b_start)),
                fabs((a_end - a_start) / (a_end + a_start)));
        ++k;
    } while (k < parameters::numerical::quark_dse_max_steps &&
             current_acc > parameters::numerical::quark_dse_acc);

    std::cout << "\nQuark DSE iteration converged after " << k << " iterations.\n";

    mat_double a2d(4);
    a2d[0] = a_values;
    a2d[1] = b_values;
    a2d[2] = dse_absci_q;
    a2d[3] = renorm;

    return a2d;
}
