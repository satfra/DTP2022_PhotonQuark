#pragma once
#include <cmath>

double a0 (const unsigned& i_)
{
  unsigned i = i_+1;
  if (i==1)
    return std::sqrt(2.);
  else if (i==7 || i==10)
    return 1.;
  return 0.;
}
