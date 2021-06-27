/** @file sign.h
 *  @brief Specification of sign
 *
 *  @author Victor Coman
 *  @date May 2021
*/

#include "sign.h"
#include <math.h>

double sign(double x, double y)
{
    double result;
    if (y > 0.0)
    {
        result = fabs(x);
    }
    else
    {
        result = -fabs(x);
    }

    return result;
}