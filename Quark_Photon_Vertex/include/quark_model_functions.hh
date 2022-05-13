#pragma once

#include <cmath>

// The model for A(x) given in the project description
double quark_a_function_model(const double& x)
{
    return ( 0.95 + 0.3 / std::log(x + 2.0) + 0.1 / (1.0 + x) + 0.29 * std::exp(-0.1 * x)
        - 0.18 * std::exp(-3.0 * x) );
}

// The model for B(x) given in the project description
double quark_m_function_model(const double& x)
{
    return 0.06 / (1.0 + x) + 0.44 * std::exp(-0.66 * x) + 0.009 / pow(std::log(x + 2.0), 0.48);
}

/*
 * This is the running coupling used in the Maris-Tandy model. This version
 * has been taken from 1606.09602v2.
 */
double maris_tandy_alpha(const double& p_squared, const double& eta, const double& lambda)
{
    constexpr double lambda_qcd = 0.234;
    constexpr double lambda_0 = 1.0;
    constexpr double gamma_m = 0.48;
    const double x = p_squared / (lambda * lambda);
    const double irterm = eta * eta *
            eta * eta * eta *
            eta * eta * M_PI * x * x *
            exp(-(eta * eta) * x);
    const double uvterm = (2.0 * M_PI * gamma_m * (1.0 - exp(-p_squared /
            (lambda_0 * lambda_0)))) / log(M_E * M_E - 1.0 + (1.0 + p_squared /
            (lambda_qcd * lambda_qcd)) * (1.0 + p_squared /
            (lambda_qcd * lambda_qcd)));
    return irterm + uvterm;
}

// This is just the function g(k^2) from eq. 19 in the project description,
// i.e. the function alpha(k^2) with some constant factors
double maris_tandy_g(const double& p_squared, const double& eta, const double& lambda)
{
    constexpr double z_2 = 0.97; // For now this is our go-to value. Has to be
                                 // adjusted when we use the real propagator.

    return z_2 * z_2 * 16.0 * M_PI * maris_tandy_alpha(p_squared, eta, lambda) / (3.0 * p_squared);
}
