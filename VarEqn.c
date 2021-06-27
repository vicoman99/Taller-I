/** @file PoleMatrix.h
 *  @brief Implementation of PoleMatrix
 *  @author Victor Coman
 *  @date May 2021
*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include "arrays.h"
#include "sat_const.h"
#include "globales.h"
#include "AccelHarmonic.h"
#include "G_AccelHarmonic.h"
#include "VarEqn.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "timediff.h"
#include "IERS.h"

void VarEqn(double x, double *yPhi, double *yPhip)
{
    extern Param AuxParam;

    int i, j;
    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC, Mjd_UT1, Mjd_TT;
    double **P, **N, **T, **E, **Phi, **G, **dfdy, **Phip;
    double *r, *v, *a;

    Phi = array(6, 6);
    r = vector(3);
    v = vector(3);

    IERS(AuxParam.Mjd_UTC, 'l', &x_pole, &y_pole, &UT1_UTC, &LOD, &dpsi, &deps, &dx_pole, &dy_pole, &TAI_UTC);
    timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);

    Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;

    // Transformation matrix
    P = PrecMatrix(MJD_J2000, AuxParam.Mjd_TT + x / 86400);
    N = NutMatrix(AuxParam.Mjd_TT + x / 86400);
    T = prod(N, 3, 3, P, 3, 3);
    E = prod(PoleMatrix(x_pole, y_pole), 3, 3, prod(GHAMatrix(Mjd_UT1), 3, 3, T, 3, 3), 3, 3);

    // State vector components
    for (i = 0; i < 3; i++)
    {
        r[i] = yPhi[i];
        r[i + 3] = yPhi[i + 3];
    }

    Phi = zeros(6, 6);

    // State transition matrix
    for (j = 0; j < 6; j++)
    {
        for (i = 0; i < 6; i++)
        {
            Phi[i][j] = yPhi[5 * j + 1 + i];
        }
    }

    // Acceleration and gradient
    a = AccelHarmonic(r, E, 3, AuxParam.n, AuxParam.m);
    G = G_AccelHarmonic(r, E, 3, AuxParam.n, AuxParam.m);

    // Time derivative of state transition matrix
    dfdy = array(6, 6);

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            dfdy[i][j] = 0.0;         // dv/dr(i,j)
            dfdy[i + 3][j] = G[i][j]; // da/dr(i,j)

            if (i == j)
            {
                dfdy[i][j + 3] = 1;
            }
            else
            {
                dfdy[i][j + 3] = 0;
                dfdy[i + 3][j + 3] = 0.0;
            }
        }
    }

    Phip = prod(dfdy, 3, 3, Phi, 3, 3);

    // Derivative of combined state vector and state transition matrix
    for (i = 0; i < 3; i++)
    {
        yPhip[i] = v[i];     // dr/dt(i)
        yPhip[i + 3] = a[i]; // dv/dt(i)
    }

    for (i = 0; i < 6; i++)
    {
        for (j = 0; j < 6; j++)
        {
            yPhip[6 * j + i] = Phip[i][j]; // dPhi/dt(i,j)
        }
    }

}