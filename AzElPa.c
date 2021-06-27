/** @file AzElPa.h
 *  @brief Implementation of AzElPa
 *  @author Victor Coman
 *  @date May 2021
*/
#include <stdio.h>
#include <math.h>
#include "arrays.h"
#include "sat_const.h"
#include "AzElPa.h"

void AzElPa(double *s, double Az, double El, double *dAds, double *dEds)
{
    double rho;

    rho = sqrt(s[0] * s[0] + s[1] * s[1]);

    // Angles
    Az = atan2(s[0], s[1]);

    if (Az < 0.0)
    {
        Az = Az + pi2;
    }

    El = atan(s[3] / rho);

    // Partials
    dAds[0] = s[1] / (rho * rho);
    dAds[1] = -s[0] / (rho * rho);
    dAds[2] = 0.0;

    dEds[0] = (-s[0]) * s[2] / (rho * dot(s, 3, s, 3));
    dEds[1] = (-s[1]) * s[2] / (rho * dot(s, 3, s, 3));
    dEds[2] = rho / dot(s, 3, s, 3);

}