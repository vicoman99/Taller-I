/** @file position.c
 *  @brief Implementation of position
 *  @author Victor Coman
 *  @date May 2021
*/
#include <math.h>
#include "sat_const.h"
#include "arrays.h"
#include "position.h"

double *position(double lon, double lat, double h)
{
    double R_equ, f, e2, CosLat, SinLat, N, *r;

    r = vector(3);

    R_equ = R_Earth;
    f = f_Earth;

    e2 = f * (2.0 - f); // Square of eccentricity
    CosLat = cos(lat);  // (Co)sine of geodetic latitude
    SinLat = sin(lat);

    // Position vector
    N = R_equ / sqrt(1.0 - e2 * SinLat * SinLat);

    r[0] = (N + h) * CosLat * cos(lon);
    r[1] = (N + h) * CosLat * sin(lon);
    r[2] = ((1.0 - e2) * N + h) * SinLat;

    return r;
}
