/** @file AccelHarmonic.h
 *  @brief Implementation of AccelHarmonic
 *  @author Victor Coman
 *  @date May 2021
*/
#include <stdio.h>
#include <math.h>
#include "arrays.h"
#include "sat_const.h"
#include "AccelHarmonic.h"
#include "Legendre.h"
#include "globales.h"

double *AccelHarmonic(double *r, double **E, int nE, int n_max, int m_max)
{
    extern double **Cnm, **Snm;
    double r_ref, gm, d, latgc, lon, dUdr, dUdlatgc, dUdlon, b1, b2, b3, r2xy;
    double *r_bf = vector(3), *a_bf = vector(3), *a = vector(3), **pnm = array(3,3), **dpnm = array(3,3);
    int q1, q2, q3;

    r_ref = 6378.1363e3; // Earths radius [m]; GGM03S
    gm = 398600.4415e9;  // [m^3/s^2]; GGM03S

    // Body-fixed position
    r_bf = matXvec(E, 3, 3, r, 3);

    // Auxiliary quantities
    d = norma(r_bf, 3); // distance
    latgc = asin(r_bf[2] / d);
    lon = atan2(r_bf[1], r_bf[0]);

    Legendre(n_max, m_max, latgc, pnm, dpnm);

    dUdr = 0;
    dUdlatgc = 0;
    dUdlon = 0;
    q3 = 0;
    q2 = q3;
    q1 = q2;

    for (int n = 0; n < n_max; n++)
    {
        b1 = (-gm / pow(d, 2.0)) * pow((r_ref / d), n) * (n + 1);
        b2 = (gm / d) * pow((r_ref / d), n);
        b3 = (gm / d) * pow((r_ref / d), n);
        for (int m = 0; m < m_max; m++)
        {
            q1 = q1 + pnm[n + 1][m + 1] * (Cnm[n + 1][m + 1] * cos(m * lon) + Snm[n + 1][m + 1] * sin(m * lon));
            q2 = q2 + dpnm[n + 1][m + 1] * (Cnm[n + 1][m + 1] * cos(m * lon) + Snm[n + 1][m + 1] * sin(m * lon));
            q3 = q3 + m * pnm[n + 1][m + 1] * (Snm[n + 1][m + 1] * cos(m * lon) - Cnm[n + 1][m + 1] * sin(m * lon));
        }
        dUdr = dUdr + q1 * b1;
        dUdlatgc = dUdlatgc + q2 * b2;
        dUdlon = dUdlon + q3 * b3;
        q3 = 0;
        q2 = q3;
        q1 = q2;
    }

    // Body-fixed acceleration
    r2xy = pow(r_bf[0], 2) + pow(r_bf[1], 2);

    a_bf = vector(3);
    a_bf[0] = (1 / d * dUdr - r_bf[2] / (pow(d, 2) * sqrt(r2xy)) * dUdlatgc) * r_bf[0] - (1 / r2xy * dUdlon) * r_bf[1];
    a_bf[1] = (1 / d * dUdr - r_bf[2] / (pow(d, 2) * sqrt(r2xy)) * dUdlatgc) * r_bf[1] + (1 / r2xy * dUdlon) * r_bf[0];
    a_bf[2] = 1 / d * dUdr * r_bf[2] + sqrt(r2xy) / pow(d, 2) * dUdlatgc;

    // Inertial acceleration
    a = matXvec(trasp(E, nE), nE, nE, a_bf, 3);

    freeArray(pnm, n_max + 2, m_max + 2);
    freeArray(dpnm, n_max + 2, m_max + 2);
    freeVector(a_bf, 3);

    return a;
}
