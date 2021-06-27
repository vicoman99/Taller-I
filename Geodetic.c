/** @file Geodetic.h
 *  @brief Implementation of Geodetic
 *  @author Victor Coman
 *  @date May 2021
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arrays.h"
#include "sat_const.h"
#include "Geodetic.h"

void Geodetic(double *r, double *lon, double *lat, double *h)
{
    double epsRequ, e2, X, Y, Z, rho2, dZ, ZdZ, Nh, SinPhi, N, dZ_new;
    epsRequ = (EPSILON) * (R_Earth);    // Convergence criterion
    e2 = (f_Earth) * (2.0 - (f_Earth)); // Square of eccentricity

    X = r[1]; // Cartesian coordinates
    Y = r[2];
    Z = r[3];
    rho2 = X * X + Y * Y; // Square of distance from z-axis

    // Check validity of input data
    if (norma(r, 3) == 0.0)
    {
        printf("invalid input in Geodetic constructor\n");
        *lon = 0.0;
        *lat = 0.0;
        *h = -(R_Earth);
    }

    // Iteration
    dZ = e2 * Z;

    while (1)
    {
        ZdZ = Z + dZ;
        Nh = sqrt(rho2 + ZdZ * ZdZ);
        SinPhi = ZdZ / Nh; // Sine of geodetic latitude
        N = (R_Earth) / sqrt(1.0 - e2 * SinPhi * SinPhi);
        dZ_new = N * e2 * SinPhi;
        if (abs(dZ - dZ_new) < epsRequ)
        {
            break;
        }
        dZ = dZ_new;
    }

    // Longitude, latitude, altitude
    *lon = atan2(Y, X);
    *lat = atan2(ZdZ, sqrt(rho2));
    *h = Nh - N;

}
