/** @file gast.h
 *  @brief Implementation of gast
 *  @author Victor Coman
 *  @date May 2021
*/
#include "Accel.h"
#include "sat_const.h"
#include "timediff.h"
#include "arrays.h"
#include "IERS.h"
#include "Mjday_TBD.h"
#include "AccelHarmonic.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "AccelPointMass.h"
#include "globales.h"
#include "NutMatrix.h"
#include "JPL_Eph_DE430.h"
#include "PrecMatrix.h"

void Accel(double x, double *Y, double *dY)
{
    extern Param AuxParam;
    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC, Mjd_UT1, Mjd_TT, MJD_TDB;
    double *r_Mercury = vector(3), *r_Venus = vector(3), *r_Earth = vector(3), *r_Mars = vector(3), *r_Jupiter = vector(3), *r_Saturn = vector(3), *r_Uranus = vector(3), *r_Neptune = vector(3), *r_Pluto = vector(3), *r_Moon = vector(3), *r_Sun = vector(3), *a = vector(3);
    double **P = array(3,3), **N = array(3,3), **T = array(3,3), **E = array(3,3);

    IERS(AuxParam.Mjd_UTC + x / 86400.0, 'l', &x_pole, &y_pole, &UT1_UTC, &LOD, &dpsi, &deps, &dx_pole, &dy_pole, &TAI_UTC);

    timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
    Mjd_UT1 = AuxParam.Mjd_UTC + x / 86400.0 + UT1_UTC / 86400.0;
    Mjd_TT = AuxParam.Mjd_UTC + x / 86400.0 + TT_UTC / 86400.0;

    P = PrecMatrix(MJD_J2000, Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = prod(N, 3, 3, P, 3, 3);

    JPL_Eph_DE430(MJD_TDB, r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun);

    E = prod(PoleMatrix(x_pole, y_pole), 3, 3, prod(GHAMatrix(Mjd_UT1), 3, 3, T, 3, 3), 3, 3);

    // Acceleration due to harmonic gravity field
    a = AccelHarmonic(Y, E, 3, AuxParam.n, AuxParam.m);

    // Luni-solar perturbations
    if (AuxParam.sun)
    {
        a = sumV(a, 3, AccelPointMass(Y, 3, r_Sun, 3, GM_Sun), 3);
    }

    if (AuxParam.moon)
    {
        a = sumV(a, 3, AccelPointMass(Y, 3, r_Moon, 3, GM_Moon), 3);
    }

    // Planetary perturbations
    if (AuxParam.planets)
    {
        a = sumV(a, 3, AccelPointMass(Y, 3, r_Mercury, 3, GM_Mercury), 3);
        a = sumV(a, 3, AccelPointMass(Y, 3, r_Venus, 3, GM_Venus), 3);
        a = sumV(a, 3, AccelPointMass(Y, 3, r_Mars, 3, GM_Mars), 3);
        a = sumV(a, 3, AccelPointMass(Y, 3, r_Jupiter, 3, GM_Jupiter), 3);
        a = sumV(a, 3, AccelPointMass(Y, 3, r_Saturn, 3, GM_Saturn), 3);
        a = sumV(a, 3, AccelPointMass(Y, 3, r_Uranus, 3, GM_Uranus), 3);
        a = sumV(a, 3, AccelPointMass(Y, 3, r_Neptune, 3, GM_Neptune), 3);
        a = sumV(a, 3, AccelPointMass(Y, 3, r_Pluto, 3, GM_Pluto), 3);
    }

    int i;
    for (i = 0; i < 3; i++)
    {
        dY[i] = Y[3 + i];
    }
    for (i = 0; i < 3; i++)
    {
        dY[i + 3] = a[i];
    }
}
