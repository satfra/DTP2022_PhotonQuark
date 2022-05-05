#include "quark_model_functions.hh"
#include "iostream"
#include "complex"
#include "vector"

// The factor of 0.7 is put into place for now, but might need double-checking.
#define FACTOR_ZERO_POINT_SEVEN 0.7
namespace sthagel_utility {
    const std::complex<double> II = {0.0, 1.0}; // NOLINT(cert-err58-cpp)
}

double sigma_v(const double k_sq)
{
    const double z = 1.0 / quark_a_function_model(FACTOR_ZERO_POINT_SEVEN * k_sq);
    const double m = quark_m_function_model(FACTOR_ZERO_POINT_SEVEN * k_sq);
    const double denom = 1.0 / (k_sq + m * m);

    return z * denom;
}

std::complex<double> kernel_G_snake(const uint8_t i, const uint8_t j, const double k_sq, const double z, const double q_sq)
{
    if (i > 11 || j > 11) {
        std::cout << "ERROR: Index out of scope in G kernel. i = " << i << " j = " << j << "." << std::endl;

    }
    // Those four values can also be passed directly, which might be a bit faster.
    // But maybe the compiler does that for us anyway.
    const double kp_sq = k_sq + 0.25 * q_sq + std::sqrt(k_sq * q_sq) * z;
    const double km_sq = k_sq + 0.25 * q_sq - std::sqrt(k_sq * q_sq) * z;
    const double m_kp = quark_m_function_model(FACTOR_ZERO_POINT_SEVEN * kp_sq);
    const double m_km = quark_m_function_model(FACTOR_ZERO_POINT_SEVEN * km_sq);

    const double sigma_m = 0.5 * (m_kp + m_km);
    const double delta_m = (m_kp - m_km) / (kp_sq - km_sq);
    const double m_bar_sq = m_kp * m_km;

    const std::complex<double> II = sthagel_utility::II;

    const double q = std::sqrt(q_sq);
    const double k = std::sqrt(k_sq);
    const double z_sq = z*z;
    const double y = std::sqrt(1.0 - z_sq);

    const double g11 = m_bar_sq + k_sq - 0.25 * q_sq;
    const std::complex<double> g12 = II * q * (sigma_m - 2.0 * k_sq * z_sq * delta_m);
    const std::complex<double> g13 = -2.0 * II * k_sq * q * z * y * delta_m;
    const double g14 = -k * q * y;
    const double g22 = m_bar_sq - (1.0 - 2.0 * z_sq) * k_sq - 0.25 * q_sq;
    const double g23 = 2.0 * k_sq * z * y;
    const std::complex<double> g24 = 2.0 * II * k * y * sigma_m;
    const double g33 = m_bar_sq + (1.0 - 2.0 * z_sq) * k_sq + 0.25 * q_sq;
    const std::complex<double> g34 =  II * k * z * (q_sq * delta_m - 2.0 * sigma_m);
    const double g44 = m_bar_sq - k_sq + 0.25 * q_sq;

    const uint8_t i_prime = i+1;
    const uint8_t j_prime = j+1;
    const unsigned int super_idx = i_prime + 100 * j_prime;

    switch (super_idx) {
    case 101:
    case 808:
    case 1212:
        return g11;
    case 102:
    case 201:
    case 708:
    case 807:
    case 1211:
    case 1112:
        return g12;
    case 103:
    case 301:
        return g13;
    case 608:
    case 806:
    case 1012:
    case 1210:
        return -g13;
    case 104:
    case 401:
        return g14;
    case 508:
    case 805:
    case 912:
    case 1209:
        return -g14;
    case 202:
    case 707:
    case 1111:
        return g22;
    case 203:
    case 302:
        return g23;
    case 607:
    case 706:
    case 1011:
    case 1110:
        return -g23;
    case 204:
    case 402:
        return g24;
    case 507:
    case 705:
    case 911:
    case 1109:
        return -g24;
    case 303:
    case 606:
    case 1010:
        return g33;
    case 304:
    case 403:
    case 506:
    case 605:
    case 910:
    case 1009:
        return g34;
    case 404:
    case 505:
    case 909:
        return g44;

    default:
        return 0.0;
    }
    return 0.0;
}

std::complex<double> kernel_G(const u_int8_t i, const u_int8_t j, const double k_sq, const double z, const double q_sq)
{
    const double kp_sq = k_sq + 0.25 * q_sq + std::sqrt(k_sq * q_sq) * z;
    const double km_sq = k_sq + 0.25 * q_sq - std::sqrt(k_sq * q_sq) * z;

    const double sig_kp = sigma_v(kp_sq);
    const double sig_km = sigma_v(km_sq);

    return sig_kp * sig_km * kernel_G_snake(i, j, k_sq, z, q_sq);
}

#undef FACTOR_ZERO_POINT_SEVEN
