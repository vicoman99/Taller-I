/** @file unit.h
 *  @brief Implementation of unit
 *
 *  @author Victor Coman
 *  @date May 2021
 */

#include "unit.h"
#include "arrays.h"

double *unit(double *vec, int n)
{
    double *outvec, small, magv;

    outvec = vector(n);
    small = 0.000001;
    magv = norma(vec, 3);

    if (magv > small)
        for (int i = 0; i < 3; i++)
        {
            outvec[i] = vec[i] / magv;
        }
    else
        for (int i = 0; i < 3; i++)
        {
            outvec[i] = 0.0;
        }
    return outvec;
}
