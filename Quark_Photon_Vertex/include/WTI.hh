#pragma once

#include <quark_model_functions.hh>
#include <cmath>

// Ward-Takahashi identity helper functions.
// Sigma_A, Delta_A, Delta_B are the symmetrized combinations of quark
// dressing functions appearing in the WTI constraint on the quark-photon vertex.
// See Eqs. (12)-(13) and (55) in the project description.

template<typename Quark>
double Sigma_A (const double &kminus2, const double &kplus2, const Quark &quark)
{
    // Eq. (12): Σ_A(k+,k-) = ½(A(k+²) + A(k-²))
    return 0.5 * (quark.A(kplus2) + quark.A(kminus2));
}

template<typename Quark>
double Delta_A (const double &kminus2, const double &kplus2, const Quark &quark)
{
    // Eq. (12): Δ_A(k+,k-) = (A(k+²) - A(k-²)) / (k+² - k-²)
   return (quark.A(kplus2) - quark.A(kminus2)) / (kplus2 - kminus2);
}

template<typename Quark>
double Delta_B (const double &kminus2, const double &kplus2, const Quark &quark)
{
    // Eq. (13): Δ_B(k+,k-) = (A(k+²)M(k+²) - A(k-²)M(k-²)) / (k+² - k-²)
   return (quark.A(kplus2)*quark.M(kplus2) - quark.A(kminus2)*quark.M(kminus2)) / (kplus2 - kminus2);
}
