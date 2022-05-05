#include "quark_model_functions.h"
#include "cmath"

// The model for A(x) given in the project description
double quark_a_function_model(double x)
{
    return 0.95 + 0.3 / std::log(x + 2.0) + 0.1 / (1.0 + x) + 0.29 * std::exp(-0.1 * x)
        - 0.18 * std::exp(-3.0 * x);
}

// The model for B(x) given in the project description
double quark_m_function_model(double x)
{
    return 0.06 / (1.0 + x) + 0.44 * std::exp(-0.66 * x) + 0.009 / pow(std::log(x + 2.0, 0.48));
}

/*
 * This is the running coupling used in the Maris-Tandy model. This version
 * has been taken from 1606.09602v2.
 */
double dressing_model_gamma_gamma(double p_squared, double eta, double lambda)
{
    const double lambda_qcd = 0.234;
    const double lambda_0 = 1.0;
    const double gamma_m = 0.48;
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
