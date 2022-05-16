#pragma once

#include <iostream>
#include <quark_model_functions.hh>
#include <basistransform.hh>
#include <cmath>
#include <Utils.hh>

double Sigma_A (double kminus2, double kplus2)

{

    double xplus = kplus2/pow(0.7,2);
    double xminus=kminus2/pow(0.7,2);

    return 0.5 * (quark_a_function_model(xplus)+quark_a_function_model(xminus));
}

double Delta_A (double kminus2, double kplus2)

{

    double xplus = kplus2/pow(0.7,2);
    double xminus=kminus2/pow(0.7,2);

   return (quark_a_function_model(xplus)-quark_a_function_model(xminus))/(kplus2 - kminus2);
}

double Delta_B (double kminus2, double kplus2)

{

    double xplus = kplus2/pow(0.7,2);
    double xminus=kminus2/pow(0.7,2);

   return (quark_a_function_model(xplus)*quark_m_function_model(xplus) - quark_a_function_model(xminus)*quark_m_function_model(xminus))/(kplus2 - kminus2);
}

//g1(const double Q, const double s, const double z,const double k, std::vector< std::complex<double> > a)

void wti (const double Q, const double s, const double z, const double k,  std::complex<double> g1 ,std::complex<double> g2, std::complex<double> g3)
{
    double kplus2 =  powr<2>(k) + powr<2>(Q)/4 + k*Q*z;
    double kminus2 =  powr<2>(k) + powr<2>(Q)/4 - k*Q*z;

    std::complex<double> wti1= g1 - std::complex<double> (Sigma_A(kminus2,kplus2));
    std::complex<double> wti2= g2 - std::complex<double>(Delta_A(kminus2,kplus2));
    std::complex<double> wti3= g3 - std::complex<double> (Delta_B(kminus2,kplus2));
   
     std::cout << "deviation of WTI given by:" << std::endl;
     std::cout << "devg1 = " << wti1 << std::endl;
     std:: cout << "devg2= " << wti2  << std::endl;
     std:: cout << "devg3= " << wti3 << std::endl;
}

