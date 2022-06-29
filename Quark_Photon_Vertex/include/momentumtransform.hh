#pragma once //include stuff only once 

#include <cmath>
#include "Utils.hh"

// BUG CHECKED

namespace momentumtransform
{
  double l_sq (const double& k_sq, const double& k_prime_sq, const double& z, const double& z_prime, const double& y)
  {
    return k_sq + k_prime_sq - 2.*std::sqrt(k_sq * k_prime_sq) * (z*z_prime + y*std::sqrt(1.-z*z)*std::sqrt(1.-z_prime*z_prime));
  } 

  double u (const double& k_sq, const double& z)
  {
    return std::sqrt(k_sq) * std::sqrt(1. - z*z);
  }

  double V (const double& k_sq, const double& k_sq_prime, const double& z, const double& z1, const double& l_sq)
  {
    return (std::sqrt(k_sq) * z - std::sqrt(k_sq_prime) * z1) / l_sq;
  }

  double w (const double& u, const double& l_sq)
  {
    return u*u / l_sq;
  }

  double X (const double& u, const double& u_prime, const double& l_sq)
  {
    return u*u_prime / l_sq;
  }
}
