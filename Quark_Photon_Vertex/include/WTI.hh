#pragma once

#include "quark_model_functions.hh"
#include <cmath>

double Sigma_A (double &kminus2, double &kplus2, const double& quark_M_p, const double& quark_M_m, const double& quark_A_p, const double& quark_A_m)
{
    return 0.5 * (quark_A_p + quark_A_m);
}

double Delta_A (double &kminus2, double &kplus2, const double& quark_M_p, const double& quark_M_m, const double& quark_A_p, const double& quark_A_m)
{
   return (quark_A_p - quark_A_m) / (kplus2 - kminus2);
}

double Delta_B (double &kminus2, double &kplus2, const double& quark_M_p, const double& quark_M_m, const double& quark_A_p, const double& quark_A_m)
{
   return (quark_A_p * quark_M_p - quark_A_m * quark_M_m) / (kplus2 - kminus2);
}
