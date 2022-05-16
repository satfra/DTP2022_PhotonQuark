#pragma once

#include <quark_model_functions.hh>
#include <cmath>

template<typename Quark>
double Sigma_A (double &kminus2, double &kplus2, const Quark &quark)
{
    return 0.5 * (quark.A(kplus2) + quark.A(kminus2));
}

template<typename Quark>
double Delta_A (double &kminus2, double &kplus2, const Quark &quark)
{
   return (quark.A(kplus2) - quark.A(kminus2)) / (kplus2 - kminus2);
}

template<typename Quark>
double Delta_B (double &kminus2, double &kplus2, const Quark &quark)
{
   return (quark.A(kplus2)*quark.M(kplus2) - quark.A(kminus2)*quark.M(kminus2)) / (kplus2 - kminus2);
}
