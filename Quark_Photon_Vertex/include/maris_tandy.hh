#pragma once

#include "parameters.hh"
#include "Utils.hh"

/*
 * This is the running coupling used in the Maris-Tandy model. This version
 * has been taken from 1606.09602v2.
 */
double maris_tandy_alpha(const double& p_squared)
{
  using namespace parameters::physical;
  const double x = p_squared / (lambda_mt * lambda_mt);
  const double irterm = eta_mt * eta_mt *
    eta_mt * eta_mt * eta_mt *
    eta_mt * eta_mt * M_PI * x * x *
    exp(-(eta_mt * eta_mt) * x);
  const double uvterm = (2.0 * M_PI * gamma_m * (1.0 - exp(-p_squared /
          (lambda_0 * lambda_0)))) / log(M_E * M_E - 1.0 + (1.0 + p_squared /
          (lambda_qcd * lambda_qcd)) * (1.0 + p_squared /
            (lambda_qcd * lambda_qcd)));
  return irterm + uvterm;
}

double pauli_villars_alpha(const double& p_squared)
{
  double const lambda_sq=parameters::physical::lambda_pv * parameters::physical::lambda_pv;
   using namespace parameters::physical;
    const double x = p_squared / (lambda_mt * lambda_mt);
    const double irterm = eta_mt * eta_mt *
            eta_mt * eta_mt * eta_mt *
            eta_mt * eta_mt * M_PI * x * x *
            exp(-(eta_mt * eta_mt) * x);
    const double uvterm = (2.0 * M_PI * gamma_m * (1.0 - exp(-p_squared /
            (lambda_0 * lambda_0)))) / log(M_E * M_E - 1.0 + (1.0 + p_squared /
            (lambda_qcd * lambda_qcd)) * (1.0 + p_squared /
            (lambda_qcd * lambda_qcd)));
    return (irterm + uvterm)/(1. + p_squared / lambda_sq);
}

template<typename Quark>
double pauli_villars_g(const double& p_squared, const Quark& quark)
{
  return powr<2>(quark.z2()) * 16.0 * M_PI * pauli_villars_alpha(p_squared) /
    (3.0 * p_squared);
}

// This is just the function g(k^2) from eq. 19 in the project description,
// i.e. the function alpha(k^2) with some constant factors
template<typename Quark>
double maris_tandy_g(const double& p_squared, const Quark& quark)
{
  return powr<2>(quark.z2()) * 16.0 * M_PI * maris_tandy_alpha(p_squared) /
    (3.0 * p_squared);
}
