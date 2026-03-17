#pragma once //include stuff only once

#include <cmath>
#include <Utils.hh>

namespace momentumtransform
{
  // Momentum variable transformations for the quark-photon vertex kinematics.
  // See Eq. (49)-(50) in the project description.

  // l² = |k - k'|²: relative squared momentum between k and k'
  // Eq. (49)
  inline double l2 (const double& k_sq, const double& k_sq_prime, const double& z, const double& z1, const double& y)
  {
    return k_sq + k_sq_prime - 2.*std::sqrt(k_sq*k_sq_prime) * (z*z1 + y*std::sqrt(1.-z*z)*std::sqrt(1.-z1*z1));
  }

  // u = k * sin(θ_k): transverse component of k w.r.t. the photon momentum
  // Eq. (50)
  inline double u (const double& k_sq, const double& z)
  {
    return std::sqrt(k_sq) * std::sqrt(1. - powr<2>(z));
  }

  // V = (k_z - k'_z) / l: longitudinal projection difference normalized by l
  // Eq. (50)
  inline double V (const double& k_sq, const double& k_sq_prime, const double& z, const double& z1, const double& l2)
  {
    return (std::sqrt(k_sq) * z - std::sqrt(k_sq_prime) * z1) / l2;
  }

  // w = u² / l²: dimensionless transverse variable for k
  // Eq. (50)
  inline double w (const double& u, const double& l2)
  {
    return powr<2>(u)/l2;
  }

  // X = u * u' / l²: mixed transverse variable coupling k and k'
  // Eq. (50)
  inline double X (const double& u, const double& u1, const double& l2)
  {
    return u*u1/l2;
  }
}
