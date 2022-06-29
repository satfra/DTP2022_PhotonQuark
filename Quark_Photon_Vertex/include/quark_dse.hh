#pragma once

#include "Utils.hh"
#include "types.hh"
#include "maris_tandy.hh"

struct weight_struc {
    const vec_double &dse_absci_q;
    const vec_double &dse_weights_q;
    const vec_double &dse_absci_ang;
    const vec_double &dse_weights_ang;
};

mat_double gauleg(double x1, double x2, int n)
{
    /*
     * This is the Gauß-Legendre method taken from the
     * "Numerical Recipes in C" Book. It has only been modified to
     * return the generated abscissas and weights as a 2d vector starting at 0.
     */
    int m,j,i;
    double xm,xl;
    vec_double x(n+1);
    vec_double w(n+1);

    m=(n+1)/2;
    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    for (i=1;i<=m;i++) {
        double p1, p2, p3, pp, z, z1;
        z=cos(M_PI*(i-0.25)/(n+0.5));
        do {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
            }
            pp=n*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;
        } while (fabs(z-z1) > 1e-9);
        x[i]=xm-xl*z;
        x[n+1-i]=xm+xl*z;
        w[i]=2.0*xl/((1.0-z*z)*pp*pp);
        w[n+1-i]=w[i];
    }

    // Preparing vector for return...
    for (int k = 0; k < n; ++k) {
        x[k] = x[k+1];
        w[k] = w[k+1];
    }
    x.pop_back();
    w.pop_back();

    mat_double a2d(2);
    a2d[0] = x;
    a2d[1] = w;

    return a2d;
}

inline double angular_function(double psi)
{
    /*
     * The angular part of the function, which is being integrated
     * to get A and B.
     */
    const double z = sin(psi);
    return z * z;
}


// TODO: Do the same for the precalculation and BSE functions.
// TODO: Check the exact momentum argument in the equations.

tens_double init_brl_angular_matrix(double mu, const weight_struc &weights)
{
    tens_double temp(2);
    mat_double temp_mat(parameters::numerical::quark_dse_steps_q);
    vec_double temp_vec(parameters::numerical::quark_dse_steps_q + 1);
    //#pragma omp parallel for
    for (unsigned int i = 0; i < parameters::numerical::quark_dse_steps_q; ++i) {
        temp_mat[i] = temp_vec;
    }
    temp[0] = temp_mat;
    temp[1] = temp_mat;

#pragma omp parallel for
    for (unsigned int qidx = 0; qidx < parameters::numerical::quark_dse_steps_q; ++qidx) {
        const double q = exp(weights.dse_absci_q[qidx] / 2.0);
        for (unsigned int pidx = 0; pidx < parameters::numerical::quark_dse_steps_q + 1; ++pidx) {
            double p;
            if (pidx == parameters::numerical::quark_dse_steps_q){
                p = mu;
            } else {
                p = exp(weights.dse_absci_q[pidx] / 2.0);
            }
            double sa = 0.0;
            double sb = 0.0;

            for (unsigned int ang_i = 0; ang_i < parameters::numerical::quark_dse_steps_z; ++ang_i) {
                const double psi = weights.dse_absci_ang[ang_i];
                const double z = cos(psi);
                const double pq = p * q * z;
                const double kk = p * p + q * q - 2.0 * p * q * z;

                const double mt = maris_tandy_alpha(kk);

                sa += weights.dse_weights_ang[ang_i] *
                        weights.dse_weights_q[qidx] * angular_function(psi) *
                        (2.0 * (M_1_PI * M_1_PI) / (3.0 * p * p)) *
                        q * q * q * q * (1.0 / kk) *
                        (mt * (
                                -2*p*p + (2*(p*p - pq)*(p*p - pq))/kk + 3*pq
                            )
                        );

                sb += weights.dse_weights_ang[ang_i] *
                        weights.dse_weights_q[qidx] * angular_function(psi) *
                        (2.0 * M_1_PI * M_1_PI) * (1.0 / 3.0) *
                        q * q * q * q * (1.0 / kk) *
                        (mt * (
                                3.0
                            )
                        );

            }
            temp[0][qidx][pidx] = sa;
            temp[1][qidx][pidx] = sb;
        }
    }
    return temp;
}

double brl_integrate_coupled_a(int pidx, const vec_double &a_values,
                               const vec_double &b_values, const tens_double &angular_matrix,
                               const weight_struc &weights)
{
    /*
     * This function does the integration, which is needed to iterate A. It uses
     * a Gaußian Quadrature using the abscissas and weights calculated earlier.
     */
    double s_a = 0.0;

    #pragma omp parallel for reduction(+:s_a)
    for (unsigned int j = 0; j < parameters::numerical::quark_dse_steps_q; ++j) {
        const double q = exp(weights.dse_absci_q[j] / 2.0);

        s_a += angular_matrix[0][j][pidx] * a_values[j] /
               ((q * a_values[j] * q * a_values[j]) +
                (b_values[j] * b_values[j]));

    }
    return s_a;
}

double brl_integrate_coupled_b(int pidx, const vec_double &a_values,
        const vec_double &b_values, const tens_double &angular_matrix,
        const vec_double &absci_x)
{
    /*
     * This function does the integration, which is needed to iterate B. It uses
     * a Gaußian Quadrature using the abscissas and weights calculated earlier.
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

mat_double quark_iterate_dressing_functions(double a0, double b0, double mc, double mu)
{
    std::cout << "Solving equation...";
     /*
     * The iterator for the dressing functions A and B. It iterates for
     * MAX_ITERATIONS steps, or until the wanted accuracy is reached. The
     * resulting dressing functions are calculated at the abscissas given in
     * weights->dse_absci_q and are returned 2d vector. At first, the values for
     * A and B are initialized to be a0 and b0 for all values of p.
     */
    vec_double a_values(parameters::numerical::quark_dse_steps_q, a0);
    vec_double b_values(parameters::numerical::quark_dse_steps_q, b0);

    const mat_double temp_q = gauleg(2.0 * std::log(parameters::physical::lambda_IR_DSE),
                               2.0 * std::log(parameters::physical::lambda_UV_DSE),
                               parameters::numerical::quark_dse_steps_q);
    const vec_double dse_absci_q = temp_q[0];
    const vec_double dse_weights_q = temp_q[1];

    const mat_double temp_z = gauleg(0.0, M_PI, parameters::numerical::quark_dse_steps_z);
    const vec_double dse_absci_z = temp_z[0];
    const vec_double dse_weights_z = temp_z[1];

    const weight_struc weights = weight_struc {
        dse_absci_q,
        dse_weights_q,
        dse_absci_z,
        dse_weights_z
    };

#pragma omp parallel for
    for (unsigned int i = 0; i < parameters::numerical::quark_dse_steps_q; ++i) {
        double p = exp(weights.dse_absci_q[i] / 2.0);
        if (p > 1.0) {
            a_values[i] = 1;
        } else {
            b_values[i] = 0.8;
        }
    }

    /*
     * Since the angular integration is the same for each iteration, it is done
     * once and will then be used in each iteration to save a lot of time.
     */
    tens_double angular_matrix = init_brl_angular_matrix(mu, weights);
    double current_acc, a_start, b_start, a_end, b_end, sigma_a, sigma_b;
    vec_double renorm = {1.0, 1.0};
    vec_double new_a(parameters::numerical::quark_dse_steps_q);
    vec_double new_b(parameters::numerical::quark_dse_steps_q);
    unsigned int k = 0;

    do{
        a_start = a_values[parameters::numerical::quark_dse_steps_q - 1];
        b_start = b_values[parameters::numerical::quark_dse_steps_q - 1];
        /*
         * This do-while loop does the iterations. In each iteration, the new
         * values for A and B are generated by integrating the previous values
         * (or initial values for the 0th step) and storing them in new_a and
         * new_b respectively. Then, the values stored in a_values and b_values
         * get updated. If the update is small enough, the iterator stops.
         */
#pragma omp parallel for
        for (unsigned int pidx = 0; pidx < parameters::numerical::quark_dse_steps_q; ++pidx) {
            /*
             * This for loop loops over all momenta in absci_x and updates the
             * values for new_a and new_b at each momentum.
             */
            new_a[pidx] = renorm[0] * (1.0 + renorm[0] *
                    brl_integrate_coupled_a(pidx, a_values, b_values,
                            angular_matrix, weights));
            new_b[pidx] = renorm[0] * (mc * renorm[1] + renorm[0] *
                 brl_integrate_coupled_b(pidx, a_values, b_values,
                            angular_matrix, weights.dse_absci_q));
        }

        /*
         * Afterwards, the new values for A and B are used to calculate the
         * renormalization constants Z2 and Zm for the current iteration.
         */
        sigma_a = brl_integrate_coupled_a(parameters::numerical::quark_dse_steps_q, a_values, b_values,
                angular_matrix, weights);
        renorm[0] = 1.0 / (1.0 + renorm[0] * sigma_a);
        sigma_b = brl_integrate_coupled_b(parameters::numerical::quark_dse_steps_q, a_values, b_values,
                                          angular_matrix, weights.dse_absci_q);
        renorm[1] = 1.0 / renorm[0] - renorm[0] * sigma_b / mc;

        for (unsigned int l = 0; l < parameters::numerical::quark_dse_steps_q; ++l) {
            /*
             * Now we update the values stored in a_values and b_values. It is
             * important to update them by value and not by reference, so that
             * the correct values are used during integration.
             */
            a_values[l] = new_a[l];
            b_values[l] = new_b[l];
        }

        // Then we check, if the error is small enough
        a_end = a_values[parameters::numerical::quark_dse_steps_q - 1];
        b_end = b_values[parameters::numerical::quark_dse_steps_q - 1];

        current_acc = std::max(fabs((b_end - b_start) / (b_end + b_start)),
                fabs((a_end - a_start) / (a_end + a_start)));

        ++k;
    } while(k < parameters::numerical::quark_dse_max_steps && current_acc > parameters::numerical::quark_dse_acc);

    std::cout << " done (" << k << " iterations";

    // Preparing array for return...
    mat_double a2d(4);
    a2d[0] = a_values;
    a2d[1] = b_values;
    a2d[2] = dse_absci_q;
    a2d[3] = renorm;

    return a2d;
}
