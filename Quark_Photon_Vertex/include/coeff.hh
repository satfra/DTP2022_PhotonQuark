#pragma once
#include <cmath>

// your code

double a0(unsigned int i)
{
    if (i==1) {
        return std::sqrt(2);
    }
    else if (i==7 || i==10) {
        return 1;
    }
    else {
        return 0;
    }
}
