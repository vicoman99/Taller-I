/** @file elements.h
 *  @brief Implementation of elements
 *  @author Victor Coman
 *  @date May 2021
*/
#include <stdlib.h>
#include <math.h>
#include "arrays.h"
#include "sat_const.h"
#include "elements.h"

void elements(double *y, double *p, double *a, double *e, double *i, double *Omega, double *omega , double *M)
{
    double *r, *v, *h;
    double magh, H, u, R, eCosE, eSinE, e2, E, nu;

    r = vector(3);
    v = vector(3);
    for (int j = 0; j < 3; j++)
    {
        r[j] = y[j];     // Position
        v[j] = y[j + 3]; // Velocity
    }
    h = crossProd(r, 3, v, 3); // Areal velocity

    magh = norma(h, 3);
    *p = magh * magh / GM_Earth;
    H = norma(h, 3);
    *Omega = atan2(h[0], -h[1]); // Long. ascend. node
    *Omega = fmod(*Omega, pi2);
    *i = atan2(sqrt(h[0] * h[0] + h[1] * h[1]), h[2]); // Inclination
    u = atan2(r[2] * H, -r[0] * h[1] + r[1] * h[0]);  // Arg. of latitude
    R = norma(r, 3);                                  // Distance
    *a = 1 / (2 / R - dot(v, 3, v, 3) / GM_Earth);     // Semi-major axis
    eCosE = 1 - R / *a;                                // e*cos(E)
    eSinE = dot(r, 3, v, 3) / sqrt(GM_Earth * *a);     // e*sin(E)
    e2 = eCosE * eCosE + eSinE * eSinE;
    *e = sqrt(e2);                                   // Eccentricity
    E = atan2(eSinE, eCosE);                        // Eccentric anomaly
    *M = fmod(E - eSinE, pi2);                       // Mean anomaly
    nu = atan2(sqrt(1.0 - e2) * eSinE, eCosE - e2); // True anomaly
    *omega = fmod(u - nu, pi2);                      // Arg. of perihelion
}
