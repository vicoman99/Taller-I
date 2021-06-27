/** @file angl.h
 *  @brief Implementation of angl
 *
 *  @author Victor Coman
 *  @date May 2021
*/

#include "angl.h"
#include "arrays.h"
#include <math.h>
#include "sign.h"

double angl(double *vec1, double *vec2)
{
    double magv1, magv2, temp, theta;
    float small = 0.00000001;
    float undefined = 999999.1;

    magv1 = norma(vec1, 3);
    magv2 = norma(vec2, 3);

    if (magv1 * magv2 > small * small)
    {
        temp = dot(vec1, 3, vec2, 3) / (magv1 * magv2);
        if (fabs(temp) > 1.0)
        {
            temp = sign(temp, 3) * 1.0;
        }

        theta = acos(temp);
    }
    else
    {
        theta = undefined;
    }

    return theta;
}
