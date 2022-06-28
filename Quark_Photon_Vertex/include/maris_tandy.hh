#pragma once

#include "parameters.hh"
#include "Utils.hh"

// BUG CHECKED

/*
 * This is the running coupling used in the Maris-Tandy model. This version
 * has been taken from 1606.09602v2.
 */
double maris_tandy_alpha(const double& p_squared)
{
  using namespace parameters::physical;

  const double x = p_squared / lambda_mt / lambda_mt;

  const double irterm = powr<7>(eta_mt) * M_PI * x * x * std::exp( - eta_mt * eta_mt * x);
  const double uvterm = (2.0 * M_PI * gamma_m * (1.0 - exp(- x))) 
   / log( M_E * M_E - 1.0 + (1.0 + p_squared / lambda_qcd / lambda_qcd) * (1.0 + p_squared / lambda_qcd / lambda_qcd));  

  return uvterm + irterm;
}

/* 
 * This is just the function g(k^2) from eq. 19 in the project description,
 * i.e. the function alpha(k^2) with some constant factors
 */
//template<typename Quark>
double maris_tandy_g(const double& p_squared, const double& z2)
{
  return z2 * z2 * 16.0 * M_PI * maris_tandy_alpha(p_squared) / (3.0 * p_squared);
}

/* 
 * This is just the function g(k^2) from eq. 19 in the project description,
 * i.e. the function alpha(k^2) with some constant factors and modified with
 * Pauli-Villars regularization
 */
//template<typename Quark>
double pauli_villars_g(const double& p_squared, const double& z2)
{
  using namespace parameters::physical;

  return maris_tandy_g(p_squared, z2) / (1. + p_squared / lambda_pv / lambda_pv);
}
