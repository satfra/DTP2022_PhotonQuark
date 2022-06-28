#pragma once

#include <iostream>
#include <complex>
#include <vector>

#include "quark_model_functions.hh"

// BUG CHECKED

// template<typename Quark>
class G
{
  private:
    static constexpr std::complex<double> II = {0.0, 1.0}; // NOLINT(cert-err58-cpp)

  public:
    // const Quark& quark;
    // kinematic variables and variables built from the quark mass function
    const double k_plus_sq, k_minus_sq, y, Sigma_M, Delta_M, M_bar_sq, sigma_p, sigma_m;  

    // kernel matrix entries
    const std::complex<double> g11, g12, g13, g14, g22, g23, g24, g33, g34, g44;

    G(const double& k_sq, const double& z, const double& q_sq, const double& quark_M_p, const double& quark_M_m, const double& quark_A_p, const double& quark_A_m) : 
    k_plus_sq(k_sq + 0.25 * q_sq + std::sqrt(k_sq * q_sq) * z),
    k_minus_sq(k_sq + 0.25 * q_sq - std::sqrt(k_sq * q_sq) * z),
    y(std::sqrt(1.0 - z*z)),
    Sigma_M(0.5 * ( quark_M_p + quark_M_m )),
    Delta_M( (quark_M_p - quark_M_m) / (k_plus_sq - k_minus_sq)),
    M_bar_sq(quark_M_p * quark_M_m),
    sigma_p(1.0 / quark_A_p / (k_plus_sq + quark_M_p * quark_M_p)),
    sigma_m(1.0 / quark_A_m / (k_minus_sq + quark_M_m * quark_M_m)),

    g11(M_bar_sq + k_sq - 0.25 * q_sq),
    g12(II * std::sqrt(q_sq) * (Sigma_M - 2.0 * k_sq * z*z * Delta_M)),
    g13(-2.0 * II * k_sq * std::sqrt(q_sq) * z * y * Delta_M),
    g14(- std::sqrt(k_sq) * std::sqrt(q_sq) * y),
    g22(M_bar_sq - (1.0 - 2.0 * z*z) * k_sq - 0.25 * q_sq),
    g23(2.0 * k_sq * z * y),
    g24(2.0 * II * std::sqrt(k_sq) * y * Sigma_M),
    g33(M_bar_sq + (1.0 - 2.0 * z*z) * k_sq + 0.25 * q_sq),
    g34(II * std::sqrt(k_sq) * z * (q_sq * Delta_M - 2.0 * Sigma_M)),
    g44(M_bar_sq - k_sq + 0.25 * q_sq)
    {

    }

    // double sigma_v(const double& k_sq) const
    // {
    //   return 1.0 / quark.A(k_sq) / (k_sq + quark.M(k_sq) * quark.M(k_sq));
    // }

    const std::complex<double> get(const unsigned& i, const unsigned& j) const
    {
      if (i > 11 || j > 11)
        throw std::runtime_error("ERROR: Index out of scope in G kernel. i = " + std::to_string(i) + " j = " + std::to_string(j) + ".");
      
      std::complex<double> G_tilde_ij;

      const unsigned i_prime = i+1;
      const unsigned j_prime = j+1;
      const unsigned int super_idx = 100 * i_prime + j_prime;

      switch (super_idx) {
        case 101:
        case 808:
        case 1212:
          G_tilde_ij = g11;
          break;
        case 102:
        case 201:
        case 708:
        case 807:
        case 1211:
        case 1112:
          G_tilde_ij = g12;
          break;
        case 103:
        case 301:
          G_tilde_ij = g13;
          break;
        case 608:
        case 806:
        case 1012:
        case 1210:
          G_tilde_ij = -g13;
          break;
        case 104:
        case 401:
          G_tilde_ij = g14;
          break;
        case 508:
        case 805:
        case 912:
        case 1209:
          G_tilde_ij = -g14;
          break;
        case 202:
        case 707:
        case 1111:
          G_tilde_ij = g22;
          break;
        case 203:
        case 302:
          G_tilde_ij = g23;
          break;
        case 607:
        case 706:
        case 1011:
        case 1110:
          G_tilde_ij = -g23;
          break;
        case 204:
        case 402:
          G_tilde_ij = g24;
          break;
        case 507:
        case 705:
        case 911:
        case 1109:
          G_tilde_ij = -g24;
          break;
        case 303:
        case 606:
        case 1010:
          G_tilde_ij = g33;
          break;
        case 304:
        case 403:
        case 506:
        case 605:
        case 910:
        case 1009:
          G_tilde_ij = g34;
          break;
        case 404:
        case 505:
        case 909:
          G_tilde_ij = g44;
          break;
        default:
          G_tilde_ij = 0.0;
      };

    // return - sigma_v(k_plus_sq) * sigma_v(k_minus_sq) * G_tilde_ij;
    return - sigma_p * sigma_m * G_tilde_ij;
    //return G_tilde_ij;
    }
};


template<typename Quark>
class G_old
{
  private:
    static constexpr std::complex<double> II = {0.0, 1.0}; // NOLINT(cert-err58-cpp)

  public:
    const Quark& quark;
    // kinematic variables and variables built from the quark mass function
    const double k_plus_sq, k_minus_sq, y, Sigma_M, Delta_M, M_bar_sq;  

    // kernel matrix entries
    const std::complex<double> g11, g12, g13, g14, g22, g23, g24, g33, g34, g44;

    G_old(const double& k_sq, const double& z, const double& q_sq, const Quark& quark_) : 
    quark(quark_),
    k_plus_sq(k_sq + 0.25 * q_sq + std::sqrt(k_sq * q_sq) * z),
    k_minus_sq(k_sq + 0.25 * q_sq - std::sqrt(k_sq * q_sq) * z),
    y(std::sqrt(1.0 - z*z)),
    Sigma_M(0.5 * ( quark_.M(k_plus_sq) + quark_.M(k_minus_sq) )),
    Delta_M( (quark_.M(k_plus_sq) - quark_.M(k_minus_sq)) / (k_plus_sq - k_minus_sq)),
    M_bar_sq(quark_.M(k_plus_sq) * quark_.M(k_minus_sq)),
    g11(M_bar_sq + k_sq - 0.25 * q_sq),
    g12(II * std::sqrt(q_sq) * (Sigma_M - 2.0 * k_sq * z*z * Delta_M)),
    g13(-2.0 * II * k_sq * std::sqrt(q_sq) * z * y * Delta_M),
    g14(- std::sqrt(k_sq) * std::sqrt(q_sq) * y),
    g22(M_bar_sq - (1.0 - 2.0 * z*z) * k_sq - 0.25 * q_sq),
    g23(2.0 * k_sq * z * y),
    g24(2.0 * II * std::sqrt(k_sq) * y * Sigma_M),
    g33(M_bar_sq + (1.0 - 2.0 * z*z) * k_sq + 0.25 * q_sq),
    g34(II * std::sqrt(k_sq) * z * (q_sq * Delta_M - 2.0 * Sigma_M)),
    g44(M_bar_sq - k_sq + 0.25 * q_sq)
    {

    }

    double sigma_v(const double& k_sq) const
    {
      return 1.0 / quark.A(k_sq) / (k_sq + quark.M(k_sq) * quark.M(k_sq));
    }

    const std::complex<double> get(const unsigned& i, const unsigned& j) const
    {
      if (i > 11 || j > 11)
        throw std::runtime_error("ERROR: Index out of scope in G kernel. i = " + std::to_string(i) + " j = " + std::to_string(j) + ".");
      
      std::complex<double> G_tilde_ij;

      const unsigned i_prime = i+1;
      const unsigned j_prime = j+1;
      const unsigned int super_idx = 100 * i_prime + j_prime;

      switch (super_idx) {
        case 101:
        case 808:
        case 1212:
          G_tilde_ij = g11;
          break;
        case 102:
        case 201:
        case 708:
        case 807:
        case 1211:
        case 1112:
          G_tilde_ij = g12;
          break;
        case 103:
        case 301:
          G_tilde_ij = g13;
          break;
        case 608:
        case 806:
        case 1012:
        case 1210:
          G_tilde_ij = -g13;
          break;
        case 104:
        case 401:
          G_tilde_ij = g14;
          break;
        case 508:
        case 805:
        case 912:
        case 1209:
          G_tilde_ij = -g14;
          break;
        case 202:
        case 707:
        case 1111:
          G_tilde_ij = g22;
          break;
        case 203:
        case 302:
          G_tilde_ij = g23;
          break;
        case 607:
        case 706:
        case 1011:
        case 1110:
          G_tilde_ij = -g23;
          break;
        case 204:
        case 402:
          G_tilde_ij = g24;
          break;
        case 507:
        case 705:
        case 911:
        case 1109:
          G_tilde_ij = -g24;
          break;
        case 303:
        case 606:
        case 1010:
          G_tilde_ij = g33;
          break;
        case 304:
        case 403:
        case 506:
        case 605:
        case 910:
        case 1009:
          G_tilde_ij = g34;
          break;
        case 404:
        case 505:
        case 909:
          G_tilde_ij = g44;
          break;
        default:
          G_tilde_ij = 0.0;
      };

    return - sigma_v(k_plus_sq) * sigma_v(k_minus_sq) * G_tilde_ij;
    
    //return G_tilde_ij;
    }
};