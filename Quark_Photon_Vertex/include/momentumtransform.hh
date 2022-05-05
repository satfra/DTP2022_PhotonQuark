#pragma once //include stuff only once 
 
#include <cmath>
#include <Utils.hh>

double l2 (const double& k, const double& k1, const double& z, const double& z1, const double& y)
{
  return k*k+k1*k1-2.*k*k1*(z*z1+y*std::sqrt(1.-z*z)*std::sqrt(1.-z1*z1));
} 

double u (const double& k, const double& z)
{
  return k * sqrt (1. - powr<2>(z));
}

double V (const double& k, const double& z, const double& k1, const double& z1, const double& l)
{
  return (k * z - k1 * z1) / (powr<2>(l))
}

double w (const double& u, const double& l)
{
  return powr<2>(u/l);
}

double X (const double& u, const double& u1, const double& l)
{
  return u*u1/powr<2>(l);
}
