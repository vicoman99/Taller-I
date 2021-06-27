/** @file anglesdr.h
 *  @brief Implementation of anglesdr
 *  @author Victor Coman
 *  @date May 2021
*/
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include "Anglesg.h"
#include "rpoly.h"
#include "globales.h"
#include "Geodetic.h"
#include "hgibbs.h"
#include "sat_const.h"
#include "arrays.h"
#include "NutMatrix.h"
#include "GHAMatrix.h"
#include "LTC.h"
#include "angl.h"
#include "elements.h"
#include "PrecMatrix.h"
#include "PoleMatrix.h"
#include "gibbs.h"
#include "IERS.h"
#include "timediff.h"
#include "doubler.h"

void anglesdr (double *r2, double *v2, double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, double *rsite1, double *rsite2, double *rsite3)
{
    extern double **eopdata;

    double magr1in = 1.1 * R_Earth;
    double magr2in = 1.11* R_Earth;
    char direct  = 'y';

    double tol    = 1e-8 * R_Earth;
    double pctchg = 0.005;

    double t1 = (Mjd1 - Mjd2)*86400.0;
    double t3 = (Mjd3 - Mjd2)*86400.0;

    double *los1 = vector(3);
    double *los2 = vector(3);
    double *los3 = vector(3);

    los1[0] = cos(el1) * sin(az1);
    los1[1] = cos(el1) * cos(az1);
    los1[2] = sin(el1);

    los2[0] = cos(el2) * sin(az2);
    los2[1] = cos(el2) * cos(az2);
    los2[2] = sin(el2);

    los3[0] = cos(el3) * sin(az3);
    los3[1] = cos(el3) * cos(az3);
    los3[2] = sin(el3);

    double lon1, lat1, h1, lon2, lat2, h2, lon3, lat3, h3;

    Geodetic(rsite1, &lon1, &lat1, &h1);
    Geodetic(rsite2, &lon2, &lat2, &h2);
    Geodetic(rsite3, &lon3, &lat3, &h3);

    double **M1 = array(3,3);
    double **M2 = array(3,3);
    double **M3 = array(3,3);

    M1 = LTC(lon1, lat1);
    M2 = LTC(lon2, lat2);
    M3 = LTC(lon3, lat3);

    //body-fixed system

    los1 = matXvec(trasp(M1, 3), 3, 3, los1, 3);
    los2 = matXvec(trasp(M2, 3), 3, 3, los2, 3);
    los3 = matXvec(trasp(M3, 3), 3, 3, los3, 3);

    //mean of date system (J2000)
    double Mjd_UTC = Mjd1;
    double *x_pole, *y_pole, UT1_UTC, *LOD, *dpsi,  *deps, *dx_pole, *dy_pole, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;

    IERS(Mjd_UTC, 'l', x_pole, y_pole, &UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, &TAI_UTC);

    timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);

    double Mjd_TT = Mjd_UTC + TT_UTC / 86400.0;
    double Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;

    double **P = PrecMatrix(MJD_J2000, Mjd_TT);
    double **N = NutMatrix(Mjd_TT);
    double **T = prod(N, 3, 3, P, 3, 3);
    double **E = prod(prod(PoleMatrix(*x_pole, *y_pole), 3, 3, GHAMatrix(Mjd_UT1), 3, 3), 3, 3, T, 3, 3);

    los1 = matXvec(trasp(E, 3), 3, 3, los1, 3);
    rsite1 = matXvec(trasp(E, 3), 3, 3, rsite1, 3);

    Mjd_UTC = Mjd2;

    IERS(Mjd_UTC, 'l', x_pole, y_pole, &UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, &TAI_UTC);
    timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC / 86400.0;
    Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;

    P = PrecMatrix(MJD_J2000, Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = prod(N, 3, 3, P, 3, 3);
    E = prod(prod(PoleMatrix(*x_pole, *y_pole), 3, 3, GHAMatrix(Mjd_UT1), 3, 3), 3, 3, T, 3, 3);

    los2 = matXvec(trasp(E, 3), 3, 3, los2, 3);
    rsite2 = matXvec(trasp(E, 3), 3, 3, rsite2, 3);

    Mjd_UTC = Mjd3;

    IERS(Mjd_UTC, 'l', x_pole, y_pole, &UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, &TAI_UTC);
    timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC / 86400.0;
    Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;

    P = PrecMatrix(MJD_J2000, Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = prod(N, 3, 3, P, 3, 3);
    E = prod(prod(PoleMatrix(*x_pole, *y_pole), 3, 3, GHAMatrix(Mjd_UT1), 3, 3), 3, 3, T, 3, 3);

    los3 = matXvec(trasp(E, 3), 3, 3, los3, 3);
    rsite3 = matXvec(trasp(E, 3), 3, 3, rsite3, 3);

    double magr1old  = 99999999.9;
    double magr2old  = 99999999.9;
    double magrsite1 = norma(rsite1, 3);
    double magrsite2 = norma(rsite2, 3);
    double magrsite3 = norma(rsite3, 3);

    double cc1 = 2.0*dot(los1, 3, rsite1, 3);
    double cc2 = 2.0*dot(los2, 3, rsite2, 3);
    double ktr = 0;

    double *r3, f1, f2, q1, magr1, magr2, a, deltae32, f, g, magr1o, deltar1, pf1pr1, pf2pr1, f1delr1, f2delr1, magr2o, f1delr2, f2delr2, pf1pr2, pf2pr2, deltar2, delta, delta1, delta2;

    while (fabs(magr1in-magr1old) > tol || fabs(magr2in-magr2old) > tol)
    {
        ktr++;
        doubler(cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct,r2,r3,&f1,&f2,&q1,&magr1,&magr2,&a,&deltae32);

        f  = 1.0 - a/magr2*(1.0-cos(deltae32));
        g  = t3 - sqrt(pow(a,3)/GM_Earth)*(deltae32-sin(deltae32));
        v2 = vec_x_esc(sumV(r3, 3, vec_x_esc(r2, 3, -f), 3), 3, 1/g);

        magr1o = magr1in;
        magr1in = (1.0+pctchg)*magr1in;
        deltar1 = pctchg*magr1in;

        doubler(cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct,r2,r3,&f1delr1,&f2delr1,&q1,&magr1,&magr2,&a,&deltae32);

        pf1pr1 = (f1delr1-f1)/deltar1;
        pf2pr1 = (f2delr1-f2)/deltar1;

        magr1in = magr1o;
        deltar1 = pctchg*magr1in;
        magr2o = magr2in;
        magr2in = (1.0+pctchg)*magr2in;
        deltar2 = pctchg*magr2in;

        doubler(cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct,r2,r3,&f1delr1,&f2delr1,&q1,&magr1,&magr2,&a,&deltae32);

        pf1pr2 = (f1delr2-f1)/deltar2;
        pf2pr2 = (f2delr2-f2)/deltar2;

        magr2in = magr2o;
        deltar2 = pctchg*magr2in;

        delta  = pf1pr1*pf2pr2 - pf2pr1*pf1pr2;
        delta1 = pf2pr2*f1 - pf1pr2*f2;
        delta2 = pf1pr1*f2 - pf2pr1*f1;

        deltar1 = -delta1/delta;
        deltar2 = -delta2/delta;

        magr1old = magr1in;
        magr2old = magr2in;

        magr1in = magr1in + deltar1;
        magr2in = magr2in + deltar2;
    }

    doubler(cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct,r2,r3,&f1,&f2,&q1,&magr1,&magr2,&a,&deltae32);

    f  = 1.0 - a/magr2*(1.0-cos(deltae32));
    g  = t3 - sqrt(pow(a,3)/GM_Earth)*(deltae32-sin(deltae32));
    v2 = vec_x_esc(sumV(r3, 3, vec_x_esc(r2, 3, -f), 3), 3, 1/g);
}
