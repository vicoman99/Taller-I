/** @file Frac.h
 *  @brief Implementation of Frac
 *  @author Victor Coman
 *  @date May 2021
*/
#include "Frac.h"
#include <math.h>

double Frac(double x)
{
    double res = x - floor(x);
    return res;
}
