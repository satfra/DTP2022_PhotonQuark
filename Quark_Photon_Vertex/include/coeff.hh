#pragma once
#include <cmath>

double a0 (const unsigned& i)
{
  if (i==1)
    return std::sqrt(2.);
  else if (i==7 || i==10)
    return 1.;
  else 
    return 0.;
}
