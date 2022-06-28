#pragma once

#include <cmath>
#include <vector>
#include "momentumtransform.hh"

using namespace momentumtransform;

class K
{
  private:
    //double y, l_sq(k_sq, k_prime_sq, z, z_prime, y), u, u(k_prime_sq, z_prime), V, w, w_prime, X(u, u(k_prime_sq, z_prime), l_sq(k_sq, k_prime_sq, z, z_prime, y));
    double k_sq, k_prime_sq, z, z_prime, y;
  public:
    K(const double& k_sq_, const double& k_prime_sq_, const double& z_, const double& z_prime_, const double& y_)
    {
      k_sq = k_sq_;
      k_prime_sq = k_prime_sq_;
      z = z_;
      z_prime = z_prime_;
      y = y_;
    }

    static bool isZeroIndex(const unsigned& i, const unsigned& j)
    {
      if (i > 11 || j > 11)
        throw std::runtime_error("Function get(..) out of range in Kernels_K");
      else if ((i == 0 && j == 1) || (i == 1 && j == 0) ||
          (i == 2 && j == 3) || (i == 3 && j == 2) ||
          (i == 3 && j == 4) || (i == 4 && j == 3) ||
          (i == 4 && j == 5) || (i == 5 && j == 4) ||
          (i == 6 && j == 7) || (i == 7 && j == 6) ||
          (i == 7 && j == 8) || (i == 8 && j == 7) ||
          (i == 8 && j == 9) || (i == 9 && j == 8) ||
          (i == 11 && j == 10) || (i == 10 && j == 11) ||
          (i == 0 && j == 2) || (i == 2 && j == 0) ||
          (i == 0 && j == 3) || (i == 3 && j == 0) ||
          (i == 0 && j == 4) || (i == 4 && j == 0) ||
          (i == 0 && j == 7) || (i == 7 && j == 0) ||
          (i == 1 && j == 3) || (i == 3 && j == 1) ||
          (i == 1 && j == 4) || (i == 4 && j == 1) ||
          (i == 1 && j == 5) || (i == 5 && j == 1) ||
          (i == 1 && j == 6) || (i == 6 && j == 1) ||
          (i == 2 && j == 4) || (i == 4 && j == 2) ||
          (i == 2 && j == 5) || (i == 5 && j == 2) ||
          (i == 2 && j == 6) || (i == 6 && j == 2) ||
          (i == 3 && j == 5) || (i == 5 && j == 3) ||
          (i == 3 && j == 6) || (i == 6 && j == 3) ||
          (i == 3 && j == 7) || (i == 7 && j == 3) ||
          (i == 4 && j == 6) || (i == 6 && j == 4) ||
          (i == 4 && j == 7) || (i == 7 && j == 4) ||
          (i == 5 && j == 7) || (i == 7 && j == 5) ||
          (i == 8 && j == 10) || (i == 10 && j == 8) ||
          (i == 8 && j == 11) || (i == 11 && j == 8) ||
          (i == 9 && j == 11) || (i == 11 && j == 9))
          return true;
      else if (i < 8 && j < 8)
        return false;
      else if (i >= 8 && j >= 8)
        return false;
      return true;
    }

    double get(const unsigned& i, const unsigned& j) const
    {
      if (i > 11 || j > 11)
        throw std::runtime_error("Function get(..) out of range in Kernels_K");

      const unsigned i_prime = i+1;
      const unsigned j_prime = j+1;
      const unsigned super_idx = 100 * i_prime + j_prime;

      switch(super_idx)
      {
        // case 101:
        //   return K11;
        // case 202:
        //   return K22;
        // case 303:
        //   return K33;
         case 101:
          return K11(k_sq, k_prime_sq, z, z_prime, y);
        case 202:
          return K22(k_sq, k_prime_sq, z, z_prime, y);
        case 303:
          return K33(k_sq, k_prime_sq, z, z_prime, y);
        case 404:
          return K44(k_sq, k_prime_sq, z, z_prime, y);
        case 505:
          return K55(k_sq, k_prime_sq, z, z_prime, y);
        case 606:
          return K66(k_sq, k_prime_sq, z, z_prime, y);
        case 707:
          return K77(k_sq, k_prime_sq, z, z_prime, y);
        case 808:
          return K88(k_sq, k_prime_sq, z, z_prime, y);
        case 106:
          return K16(k_sq, k_prime_sq, z, z_prime, y);
        case 601:
          return K61(k_sq, k_prime_sq, z, z_prime, y);
        case 107:
          return K17(k_sq, k_prime_sq, z, z_prime, y);
        case 701:
          return K71(k_sq, k_prime_sq, z, z_prime, y); 
        case 607:
          return K67(k_sq, k_prime_sq, z, z_prime, y);
        case 706:
          return K76(k_sq, k_prime_sq, z, z_prime, y);
        case 203:
          return K23(k_sq, k_prime_sq, z, z_prime, y);
        case 302:
          return K32(k_sq, k_prime_sq, z, z_prime, y);
        case 208:
          return K71(k_sq, k_prime_sq, z, z_prime, y) + std::sqrt(2.) * (1. - y * y);
        case 802:
          return K17(k_sq, k_prime_sq, z, z_prime, y) + std::sqrt(2.) * (1. - y * y);
        case 308:
          return -K61(k_sq, k_prime_sq, z, z_prime, y);
        case 803:
          return -K16(k_sq, k_prime_sq, z, z_prime, y);

        case 909:
          return K55(k_sq, k_prime_sq, z, z_prime, y)/y;
        case 1010:
          return K66(k_sq, k_prime_sq, z, z_prime, y)/y;
        case 1011:
          return K67(k_sq, k_prime_sq, z, z_prime, y)/y;
        case 1110:
          return K76(k_sq, k_prime_sq, z, z_prime, y)/y;
        case 1111:
          return K77(k_sq, k_prime_sq, z, z_prime, y)/y;
        case 1212:
          return K88(k_sq, k_prime_sq, z, z_prime, y)/y;
      }
      return 0.;
    }

  public:
    double K11(const double& k_sq, const double& k_prime_sq, const double& z, const double& z_prime, const double& y) const
    {
      return -(1. + y * y) / 2. - y * (1. - y * y) * X(u(k_sq, z), u(k_prime_sq, z_prime), l_sq(k_sq, k_prime_sq, z, z_prime, y));
    }

    double K22(const double& k_sq, const double& k_prime_sq, const double& z, const double& z_prime, const double& y) const
    {
      return -(1. + y * y) / 2. * (1. - 2. * l_sq(k_sq, k_prime_sq, z, z_prime, y) * V(k_sq, k_prime_sq, z, z_prime, l_sq(k_sq, k_prime_sq, z, z_prime, y)) * V(k_sq, k_prime_sq, z, z_prime, l_sq(k_sq, k_prime_sq, z, z_prime, y))) + y * (1. - y * y) * X(u(k_sq, z), u(k_prime_sq, z_prime), l_sq(k_sq, k_prime_sq, z, z_prime, y));
    }

    double K33(const double& k_sq, const double& k_prime_sq, const double& z, const double& z_prime, const double& y) const
    {
      return  y * (1. - 2. * l_sq(k_sq, k_prime_sq, z, z_prime, y) * V(k_sq, k_prime_sq, z, z_prime, l_sq(k_sq, k_prime_sq, z, z_prime, y)) * V(k_sq, k_prime_sq, z, z_prime, l_sq(k_sq, k_prime_sq, z, z_prime, y))) - (1. - y * y) * X(u(k_sq, z), u(k_prime_sq, z_prime), l_sq(k_sq, k_prime_sq, z, z_prime, y));
    }
    double K44(const double& k_sq, const double& k_prime_sq, const double& z, const double& z_prime, const double& y) const
    {
      return y + (1. - y * y) * X(u(k_sq, z), u(k_prime_sq, z_prime), l_sq(k_sq, k_prime_sq, z, z_prime, y));
    }

    double K55(const double& k_sq, const double& k_prime_sq, const double& z, const double& z_prime, const double& y) const
    {
      return 3. * y;
    }

    double K66(const double& k_sq, const double& k_prime_sq, const double& z, const double& z_prime, const double& y) const
    {
      return -y * (1. + 2. * l_sq(k_sq, k_prime_sq, z, z_prime, y) * V(k_sq, k_prime_sq, z, z_prime, l_sq(k_sq, k_prime_sq, z, z_prime, y)) * V(k_sq, k_prime_sq, z, z_prime, l_sq(k_sq, k_prime_sq, z, z_prime, y)));
    }

    double K77(const double& k_sq, const double& k_prime_sq, const double& z, const double& z_prime, const double& y) const
    {
      return -y * y * (3. - 2. * l_sq(k_sq, k_prime_sq, z, z_prime, y) * V(k_sq, k_prime_sq, z, z_prime, l_sq(k_sq, k_prime_sq, z, z_prime, y)) * V(k_sq, k_prime_sq, z, z_prime, l_sq(k_sq, k_prime_sq, z, z_prime, y))) + 2. * y * (1. - y * y) * X(u(k_sq, z), u(k_prime_sq, z_prime), l_sq(k_sq, k_prime_sq, z, z_prime, y));
    }

    double K88(const double& k_sq, const double& k_prime_sq, const double& z, const double& z_prime, const double& y) const
    {
      return y * y - 2. * y * (1. - y * y) * X(u(k_sq, z), u(k_prime_sq, z_prime), l_sq(k_sq, k_prime_sq, z, z_prime, y));
    }

    double K16(const double& k_sq, const double& k_prime_sq, const double& z, const double& z_prime, const double& y) const
    {
      return std::sqrt(2.) * (1. - y * y) * u(k_prime_sq, z_prime) * V(k_sq, k_prime_sq, z, z_prime, l_sq(k_sq, k_prime_sq, z, z_prime, y));
    }

    double K61(const double& k_sq, const double& k_prime_sq, const double& z, const double& z_prime, const double& y) const
    {
      return -std::sqrt(2.) * (1. - y * y) * u(k_sq, z) * V(k_sq, k_prime_sq, z, z_prime, l_sq(k_sq, k_prime_sq, z, z_prime, y));
    }

    double K17(const double& k_sq, const double& k_prime_sq, const double& z, const double& z_prime, const double& y) const
    {
      return -(1. - y * y) / std::sqrt(2.) * (1. + 2. * w(u(k_prime_sq, z_prime), l_sq(k_sq, k_prime_sq, z, z_prime, y)) - 2. * y * X(u(k_sq, z), u(k_prime_sq, z_prime), l_sq(k_sq, k_prime_sq, z, z_prime, y)));
    }

    double K71(const double& k_sq, const double& k_prime_sq, const double& z, const double& z_prime, const double& y) const
    {
      return -(1. - y * y) / std::sqrt(2.) * (1. + 2 * w(u(k_sq, z), l_sq(k_sq, k_prime_sq, z, z_prime, y)) - 2. * y * X(u(k_sq, z), u(k_prime_sq, z_prime), l_sq(k_sq, k_prime_sq, z, z_prime, y)));
    }

    double K23(const double& k_sq, const double& k_prime_sq, const double& z, const double& z_prime, const double& y) const
    {
      return (2. * y * u(k_sq, z) - (1. + y * y) * u(k_prime_sq, z_prime)) * V(k_sq, k_prime_sq, z, z_prime, l_sq(k_sq, k_prime_sq, z, z_prime, y));
    }

    double K32(const double& k_sq, const double& k_prime_sq, const double& z, const double& z_prime, const double& y) const
    {
      return -(2. * y * u(k_prime_sq, z_prime) - (1. + y * y) * u(k_sq, z)) * V(k_sq, k_prime_sq, z, z_prime, l_sq(k_sq, k_prime_sq, z, z_prime, y));
    }

    double K67(const double& k_sq, const double& k_prime_sq, const double& z, const double& z_prime, const double& y) const
    {
      return 2. * y * (u(k_prime_sq, z_prime) - y * u(k_sq, z)) * V(k_sq, k_prime_sq, z, z_prime, l_sq(k_sq, k_prime_sq, z, z_prime, y));
    }

    double K76(const double& k_sq, const double& k_prime_sq, const double& z, const double& z_prime, const double& y) const
    {
      return -2. * y * (u(k_sq, z) - y * u(k_prime_sq, z_prime)) * V(k_sq, k_prime_sq, z, z_prime, l_sq(k_sq, k_prime_sq, z, z_prime, y));
    }

};
