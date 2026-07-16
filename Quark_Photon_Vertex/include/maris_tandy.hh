#pragma once

#include "parameters.hh"
#include "Utils.hh"
#include "acc.hh"

/*
 * This is the running coupling used in the Maris-Tandy model. This version
 * has been taken from 1606.09602v2.
 * Implements Eq. (19) in the project description.
 */
// eta_mt and lambda_mt are the per-flavour Maris-Tandy IR parameters, passed
// in at runtime (see flavor.hh). The UV parameters (gamma_m, lambda_0,
// lambda_qcd) are flavour-independent and stay in parameters::physical.
ACC_ROUTINE_SEQ
inline double maris_tandy_alpha(const double& p_squared, const double eta_mt, const double lambda_mt)
{
  using namespace parameters::physical;
  const double x = p_squared / powr<2>(lambda_mt);
  const double irterm = powr<7>(eta_mt) * M_PI * powr<2>(x) * std::exp( - powr<2>(eta_mt) * x);
  const double uvterm = (2.0 * M_PI * gamma_m * (1.0 - exp(-p_squared / powr<2>(lambda_0))) )
    / std::log(powr<2>(M_E) - 1.0 + (1.0 + p_squared / powr<2>(lambda_qcd)) * (1.0 + p_squared / powr<2>(lambda_qcd)));
  return irterm + uvterm;
}

/*
 * This is just the function g(k^2) from Eq. (19) in the project description,
 * i.e. the function alpha(k^2) with some constant factors.
 * Used as the gluon interaction kernel in the BSE, see Eq. (64).
 *
 * Scalar overload: takes the three quark quantities (Z_2, eta, lambda) as plain
 * doubles so it can be called from an OpenACC device region without mapping a
 * quark object (which, in the DSE case, holds host-side splines). The templated
 * overload below just extracts these three scalars and delegates, so host
 * callers are unchanged and the arithmetic is identical.
 */
ACC_ROUTINE_SEQ
inline double maris_tandy_g(const double& p_squared, const double z2,
    const double eta_mt, const double lambda_mt)
{
  return powr<2>(z2) * 16.0 * M_PI *
    maris_tandy_alpha(p_squared, eta_mt, lambda_mt)
    / (3.0 * p_squared);
}

template<typename Quark>
double maris_tandy_g(const double& p_squared, const Quark& quark)
{
  return maris_tandy_g(p_squared, quark.z2(), quark.eta_mt(), quark.lambda_mt());
}

/*
 * This is just the function g(k^2) from Eq. (19) in the project description,
 * i.e. the function alpha(k^2) with some constant factors and modified with
 * Pauli-Villars regularization. See Eq. (64). Scalar/templated overload pair as
 * for maris_tandy_g above.
 */
ACC_ROUTINE_SEQ
inline double pauli_villars_g(const double& p_squared, const double z2,
    const double eta_mt, const double lambda_mt)
{
  using namespace parameters::physical;
  constexpr double lambda_sq = powr<2>(lambda_pv);

  return powr<2>(z2) * 16.0 * M_PI *
    ( maris_tandy_alpha(p_squared, eta_mt, lambda_mt) / (1. + p_squared / lambda_sq) )
    / (3.0 * p_squared);
}

template<typename Quark>
double pauli_villars_g(const double& p_squared, const Quark& quark)
{
  return pauli_villars_g(p_squared, quark.z2(), quark.eta_mt(), quark.lambda_mt());
}
