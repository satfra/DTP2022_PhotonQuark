#pragma once

#include <iostream>
#include <quark_model_functions.hh>
#include <basistransform.hh>
#include <cmath>
#include <Utils.hh>

template<typename Quark>
double Sigma_A (double &kminus2, double &kplus2, const Quark &quark)

{

    
    return 0.5 * (quark.A(kplus2)+quark.A(kminus2));
}

template<typename Quark>
double Delta_A (double &kminus2, double &kplus2, const Quark &quark)

{

    

   return (quark.A(kplus2)-quark.A(kminus2))/(kplus2 - kminus2);
}


template<typename Quark>
double Delta_B (double &kminus2, double &kplus2, const Quark &quark)


{

   return (quark.A(kplus2)*quark.M(kplus2) - quark.A(kminus2)*quark.M(kminus2))/(kplus2 - kminus2);
}

//g1(const double Q, const double s, const double z,const double k, std::vector< std::complex<double> > a)
template<typename Quark> 
void wti (const double &Q, const double &s, const double &z, const double &k,  std::complex<double> &g1 ,std::complex<double> &g2, std::complex<double> &g3, const Quark &quark)
{
    double kplus2 =  powr<2>(k) + powr<2>(Q)/4 + k*Q*z;
    double kminus2 =  powr<2>(k) + powr<2>(Q)/4 - k*Q*z;

    std::complex<double> wti1= g1 - std::complex<double> (Sigma_A<Quark>(kminus2,kplus2,quark));
    std::complex<double> wti2= g2 - std::complex<double>(Delta_A<Quark>(kminus2,kplus2,quark));
    std::complex<double> wti3= g3 - std::complex<double> (Delta_B<Quark>(kminus2,kplus2,quark));

   
     std::cout << "for parameters (Q,s,z,k) = " << "( " << Q << "," << s << "," << z << ", " << k << ")" << std::endl;
     std::cout << "deviation of WTI given by:" << std::endl;
     std::cout << "devg1 = " << wti1 << std::endl;
     std:: cout << "devg2= " << wti2  << std::endl;
     std:: cout << "devg3= " << wti3 << std::endl;
}

