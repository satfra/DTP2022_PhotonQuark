#pragma once //include stuff only once 

#include <cmath>
#include <Utils.hh>

namespace momentumtransform
{
  double l2 (const double& k_sq, const double& k_sq_prime, const double& z, const double& z1, const double& y)
  {
    return k_sq + k_sq_prime - 2.*std::sqrt(k_sq*k_sq_prime) * (z*z1 + y*std::sqrt(1.-z*z)*std::sqrt(1.-z1*z1));
  } 

  double u (const double& k_sq, const double& z)
  {
    return std::sqrt(k_sq) * std::sqrt(1. - powr<2>(z));
  }

  double V (const double& k_sq, const double& k_sq_prime, const double& z, const double& z1, const double& l2)
  {
    return (std::sqrt(k_sq) * z - std::sqrt(k_sq_prime) * z1) / l2;
  }

  double w (const double& u, const double& l2)
  {
    return powr<2>(u)/l2;
  }

  double X (const double& u, const double& u1, const double& l2)
  {
    return u*u1/l2;
  }
}
