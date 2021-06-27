/** @file testIOD_KF.h
 *  @brief Implementation of all the tests
 *
 *  @author Victor Coman
 *  @date June 2021
 */

#include <stdio.h>
#include <float.h>
#include <math.h>
#include "arrays.h"
#include "sat_const.h"
#include "position.h"
#include "MeanObliquity.h"
#include "R_x.h"
#include "R_Y.h"
#include "R_Z.h"
#include "NutAngles.h"
#include "NutMatrix.h"
#include "Mjday_TBD.H"
#include "unit.h"
#include "Mjday.h"
#include "ode.h"
#include "timediff.h"
#include "IERS.h"
#include "PrecMatrix.h"
#include "PoleMatrix.h"
#include "gast.h"
#include "GHAMatrix.h"
#include "AzElPa.h"
#include "Legendre.h"
#include "Accel.h"
#include "LTC.h"
#include "JPL_Eph_DE430.h"
#include "AccelHarmonic.h"
#include "timediff.h"
#include "EqnEquinox.h"
#include "angl.h"
#include "Cheb3D.h"
#include "EccAnom.h"
#include "elements.h"
#include "Frac.h"
#include "hgibbs.h"
#include "AccelPointMass.h"
#include "TimeUpdate.h"
#include "sign.h"
#include "G_AccelHarmonic.h"
#include "gibbs.h"
#include "gmst.h"
#include "MeasUpdate.H"
#include "Geodetic.h"
#include "doubler.h"

/*
#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) \
    do                \
    {                 \
        if (!(test))  \
        {             \
            FAIL();   \
            return 1; \
        }             \
    } while (0)
#define _verify(test)   \
    do                  \
    {                   \
        int r = test(); \
        tests_run++;    \
        if (r)          \
            return r;   \
    } while (0)

int tests_run = 0;

int position_01()
{
    double lon, lat, alt, *Rs, *sol;

    sol = vector(3);
    lon = -2.76234307910694;
    lat = 0.376551295459273;
    alt = 300.2;
    sol[0] = -5512567.84003607;
    sol[1] = -2196994.44666933;
    sol[2] = 2330804.96614689;

    Rs = position(lon, lat, alt);

    _assert(compareV(Rs, 3, sol, 3));

    freeVector(Rs, 3);

    return 0;
}

int dot_01()
{
    double *v1, *v2, sol;

    v1 = vector(3);
    v2 = vector(3);

    v1[0] = 1;
    v1[1] = 2;
    v1[2] = 3;
    v2[0] = 1;
    v2[1] = 2;
    v2[2] = 3;
    sol = 14.0;

    _assert(fabs(sol - dot(v1, 3, v2, 3)) < pow(10, -10));

    freeVector(v1, 3);
    freeVector(v2, 3);

    return 0;
}

int norma_01()
{
    double *w, sol;

    w = vector(3);
    w[0] = 2.0;
    w[1] = 0.0;
    w[2] = 0.0;
    sol = 2.0;

    _assert(fabs(sol - norma(w, 3)) < pow(10, -10));

    freeVector(w, 3);

    return 0;
}

int trasp_01()
{
    double **m1, **sol, **tr;

    sol = array(3, 3);
    m1 = array(3, 3);

    m1[0][0] = 1;
    m1[0][1] = 0;
    m1[0][2] = 2;
    m1[1][0] = -1;
    m1[1][1] = 5;
    m1[1][2] = 0;
    m1[2][0] = 0;
    m1[2][1] = 3;
    m1[2][2] = -9;

    sol[0][0] = 1;
    sol[0][1] = -1;
    sol[0][2] = 0;
    sol[1][0] = 0;
    sol[1][1] = 5;
    sol[1][2] = 3;
    sol[2][0] = 2;
    sol[2][1] = 0;
    sol[2][2] = -9;

    tr = trasp(m1, 3);

    _assert(compare(tr, 3, 3, sol, 3, 3));

    freeArray(sol, 3, 3);
    freeArray(tr, 3, 3);
    freeArray(m1, 3, 3);

    return 0;
}

int inv_01()
{
    double **m1, **sol, **in;

    sol = array(3, 3);
    m1  = array(3, 3);

    m1[0][0] = 1; m1[0][1] = 0; m1[0][2] = 2;
    m1[1][0] = -1; m1[1][1] = 5; m1[1][2] = 0;
    m1[2][0] = 0; m1[2][1] = 3; m1[2][2] = -9;
    sol[0][0] = 0.88235294117647; sol[0][1] = -0.11764705882352; sol[0][2] = 0.1960784313725;
    sol[1][0] = 0.17647058823529; sol[1][1] = 0.17647058823529; sol[1][2] = 0.03921568627451;
    sol[2][0] = 0.05882352941176; sol[2][1] = 0.05882352941176; sol[2][2] = -0.09803921568627;

    in = inv(m1, 3);
    printArray(sol,3,3);
    printArray(in,3,3);

    sol = array(3, 3);
    sol[0][0] = 1; sol[1][1] = 1; sol[2][2] = 1;
    m1 = eye(3);

    m1[0][0] = 3; m1[0][1] = 0; m1[0][2] = 1;
    m1[1][0] = 3; m1[1][1] = 0; m1[1][2] = 0;
    m1[2][0] = 5; m1[2][1] = 1; m1[2][2] = 1;


    in = inv(m1, 3);
    // printArray(sol,3,3);
    printArray(in,3,3);

    _assert(compare(in, 3, 3, sol, 3, 3));

    freeArray(sol, 3, 3);
    freeArray(in, 3, 3);
    freeArray(m1, 3, 3);

    return 0;
}


int prod_01()
{
    double **m, **m1, **m2, **sol;

    sol = array(3, 3);
    m1 = array(3, 3);
    m2 = array(3, 3);

    m1[0][0] = 1;
    m1[0][1] = 2;
    m1[0][2] = 3;
    m1[1][0] = 4;
    m1[1][1] = 5;
    m1[1][2] = 6;
    m1[2][0] = 7;
    m1[2][1] = 8;
    m1[2][2] = 9;
    m2[0][0] = 1;
    m2[0][1] = 2;
    m2[0][2] = 3;
    m2[1][0] = 4;
    m2[1][1] = 5;
    m2[1][2] = 6;
    m2[2][0] = 7;
    m2[2][1] = 8;
    m2[2][2] = 9;

    sol[0][0] = 30;
    sol[0][1] = 36;
    sol[0][2] = 42;
    sol[1][0] = 66;
    sol[1][1] = 81;
    sol[1][2] = 96;
    sol[2][0] = 102;
    sol[2][1] = 126;
    sol[2][2] = 150;

    m = prod(m1, 3, 3, m2, 3, 3);

    _assert(compare(m, 3, 3, sol, 3, 3));

    freeArray(sol, 3, 3);
    freeArray(m, 3, 3);
    freeArray(m1, 3, 3);
    freeArray(m2, 3, 3);

    return 0;
}

int eye_01()
{
    double **m, **sol;

    sol = array(3, 3);
    sol[0][0] = 1;
    sol[1][1] = 1;
    sol[2][2] = 1;
    m = eye(3);

    _assert(compare(m, 3, 3, sol, 3, 3));

    freeArray(sol, 3, 3);
    freeArray(m, 3, 3);

    return 0;
}

int cross_01()
{
    double *v1, *v2, *res, *sol;

    v1 = vector(3);
    v2 = vector(3);
    sol = vector(3);

    v1[0] = 3;
    v1[1] = 2;
    v1[2] = 2;
    v2[0] = 5;
    v2[1] = 4;
    v2[2] = 1;

    sol[0] = -6;
    sol[1] = 7;
    sol[2] = 2;
    res = crossProd(v1, 3, v2, 3);

    _assert(compareV(res, 3, sol, 3));

    freeVector(v1, 3);
    freeVector(v2, 3);
    freeVector(sol, 3);

    return 0;
}

int sum_01()
{
    double **m, **m1, **m2, **sol;

    sol = array(3, 3);
    m1 = array(3, 3);
    m2 = array(3, 3);

    m1[0][0] = 1;
    m1[0][1] = 2;
    m1[0][2] = 3;
    m1[1][0] = 1;
    m1[1][1] = 2;
    m1[1][2] = 3;
    m1[2][0] = 1;
    m1[2][1] = 2;
    m1[2][2] = 3;
    m2[0][0] = 1;
    m2[0][1] = 2;
    m2[0][2] = 3;
    m2[1][0] = 1;
    m2[1][1] = 2;
    m2[1][2] = 3;
    m2[2][0] = 1;
    m2[2][1] = 2;
    m2[2][2] = 3;

    sol[0][0] = 2;
    sol[0][1] = 4;
    sol[0][2] = 6;
    sol[1][0] = 2;
    sol[1][1] = 4;
    sol[1][2] = 6;
    sol[2][0] = 2;
    sol[2][1] = 4;
    sol[2][2] = 6;

    m = sum(m1, 3, 3, m2, 3, 3);

    _assert(compare(m, 3, 3, sol, 3, 3));

    freeArray(sol, 3, 3);
    freeArray(m, 3, 3);
    freeArray(m1, 3, 3);
    freeArray(m2, 3, 3);

    return 0;
}

// el array se gira 1 RAD como se pidio en clase
int R_x_01()
{
    double **mat, **sol;

    sol = array(3, 3);
    mat = array(3, 3);

    sol[0][0] = 1.0;
    sol[0][1] = 0.0;
    sol[0][2] = 0.0;
    sol[1][0] = 0.0;
    sol[1][1] = 0.5403023059;
    sol[1][2] = 0.8414709848;
    sol[2][0] = 0.0;
    sol[2][1] = -0.8414709848;
    sol[2][2] = 0.5403023059;

    mat = R_x(1.0);

    _assert(compare(mat, 3, 3, sol, 3, 3));

    freeArray(sol, 3, 3);
    freeArray(mat, 3, 3);

    return 0;
}

int R_y_01()
{
    double **mat, **sol;

    sol = array(3, 3);
    mat = array(3, 3);

    sol[0][0] = 0.54030230586814;
    sol[0][1] = 0.0;
    sol[0][2] = -0.841470984807897;
    sol[1][0] = 0.0;
    sol[1][1] = 1.0;
    sol[1][2] = 0.0;
    sol[2][0] = 0.841470984807897;
    sol[2][1] = 0.0;
    sol[2][2] = 0.54030230586814;

    mat = R_y(1.0);

    _assert(compare(mat, 3, 3, sol, 3, 3));

    freeArray(sol, 3, 3);
    freeArray(mat, 3, 3);

    return 0;
}

int R_z_01()
{
    double **mat, **sol;

    sol = array(3, 3);
    mat = array(3, 3);

    sol[0][0] = 0.54030230586814;
    sol[0][1] = 0.841470984807897;
    sol[0][2] = 0.0;
    sol[1][0] = -0.841470984807897;
    sol[1][1] = 0.54030230586814;
    sol[1][2] = 0.0;
    sol[2][0] = 0.0;
    sol[2][1] = 0.0;
    sol[2][2] = 1.0;

    mat = R_z(1.0);

    _assert(compare(mat, 3, 3, sol, 3, 3));

    freeArray(sol, 3, 3);
    freeArray(mat, 3, 3);

    return 0;
}

int MeanObliquity_01()
{
    double sol = 0.409337689628447;

    _assert(fabs(sol - MeanObliquity(12133.384924368925)) < pow(10, -10));

    return 0;
}

int Mjday_TBD_01()
{
    double sol = 47568.3124245554;

    _assert(fabs(sol - Mjd_TBD(47568.3124245436)) < pow(10, -10));

    return 0;
}

int NutAngles_01()
{
    double dpsi, deps, dpsiSol, depsSol;

    dpsiSol = 6.23063736216799e-05;
    depsSol = -3.51110708894389e-05;

    //printf("dpsiSol %e  ", dpsiSol);
    //printf("depsSol %e  ", depsSol);

    NutAngles(49746.1097222222, &dpsi, &deps);

    //printf("dpsi %e  ", dpsi);
    //printf("deps %e  ", deps);

    _assert((fabs(dpsiSol - dpsi) < pow(10, -10)) && (fabs(depsSol - deps) < pow(10, -10)));

    return 0;
}

int EqnEquinox_01()
{
    double sol = -8.64351363239508e-05;

    _assert(fabs(sol - EqnEquinox(23134434.2356478)) < pow(10, -10));

    return 0;
}

int NutMatrix_01()
{
    double **mat, **sol;
    mat = array(3, 3);
    sol = array(3, 3);

    sol[0][0] = 0.999999999212309;
    sol[0][1] = 3.64111318149062e-05;
    sol[0][2] = 1.57990920133915e-05;
    sol[1][0] = -3.64117241172012e-05;
    sol[1][1] = 0.999999998634307;
    sol[1][2] = 3.74909735167406e-05;
    sol[2][0] = -1.57977269030362e-05;
    sol[2][1] = -3.74915487593719e-05;
    sol[2][2] = 0.999999999172408;

    mat = NutMatrix(3425.45235235);

    _assert(compare(mat, 3, 3, sol, 3, 3));

    freeArray(sol, 3, 3);
    freeArray(mat, 3, 3);

    return 0;
}

int MeasUpdate_01(){
    double *Ksol, **Xsol, **Psol;

    Ksol[0] = -0.0529207814929918;
    Ksol[1] = 0.185257988750088;
    Ksol[2] = 0.08901996481035;

    Xsol[0][0] = 1.04762870334369;
    Xsol[0][1] = 3.04762870334369;
    Xsol[0][2] = 4.04762870334369;
    Xsol[1][0] = 0.833267810124921;
    Xsol[1][1] = 2.83326781012492;
    Xsol[1][2] = 3.83326781012492;
    Xsol[2][0] = 0.919882031670685;
    Xsol[2][1] = 2.91988203167069;
    Xsol[2][2] = 3.91988203167069;

    Psol[0][0] = 0.65795813037657;
    Psol[0][1] = 0.211683125971967;
    Psol[0][2] = -0.828815846389593;
    Psol[1][0] = -0.411873762976301;
    Psol[1][1] = 0.258968044999647;
    Psol[1][2] = -0.0443014147672662;
    Psol[2][0] = 0.643557857639505;
    Psol[2][1] = -0.3560798592414;
    Psol[2][2] = 0.519014638814004;


    double *x, *G, **P, *K;

    K = vector(3);

    x[0] = 1;
    x[1] = 3;
    x[2] = 4;

    G[0] = 1;
    G[1] = 4;
    G[2] = 2;

    MeasUpdate(x, 3, 1.4, 2.3, 1.7, G, 3, R_y(1), 3, 3, K);

    _assert(compareV(K, 3, Ksol, 3) && compare(P, 3, 3, Psol, 3, 3));

    return 0;
}

int unit_01()
{
    double sol[] = {0.267261241912424, 0.534522483824849, 0.801783725737273};

    double vec[] = {1, 2, 3};

    double *res = unit(vec, 3);

    _assert(compareV(res, 3, sol, 3));

    freeVector(sol, 3);
    freeVector(vec, 3);
    freeVector(res, 3);

    return 0;
}

int timediff_01()
{
    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    double UT1_TAIres, UTC_GPSres, UT1_GPSres, TT_UTCres, GPS_UTCres;

    double a = 2.555;
    double b = 4.888;

    timediff(a, b, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);

    //printf("UT1_TAI: %f  ", UT1_TAI);
    //printf("UTC_GPS: %f  ", UTC_GPS);
    //printf("UT1_GPS: %f  ", UT1_GPS);
    //printf("TT_UTC: %f  ", TT_UTC);
    //printf("GPS_UTC: %f  ", GPS_UTC);

    UT1_TAIres = -2.333000;
    UTC_GPSres = 14.112000;
    UT1_GPSres = 16.667000;
    TT_UTCres = 37.072000;
    GPS_UTCres = -14.112000;

    _assert(fabs(UT1_TAIres - UT1_TAI) < pow(10, -3) && fabs(UTC_GPSres - UTC_GPS) < pow(10, -3) && fabs(UT1_GPSres - UT1_GPS) < pow(10, -3) && fabs(TT_UTCres - TT_UTC) < pow(10, -3) && fabs(GPS_UTCres - GPS_UTC) < pow(10, -3));

    return 0;
}

int matXvec_01()
{
    double **m1 = array(3, 3);
    double *vec = vector(3);
    double *sol = vector(3);;
    double *v = vector(3);;

    m1[0][0] = 1;
    m1[0][1] = 2;
    m1[0][2] = 3;
    m1[1][0] = 4;
    m1[1][1] = 5;
    m1[1][2] = 6;
    m1[2][0] = 7;
    m1[2][1] = 8;
    m1[2][2] = 9;

    vec[0] = 1;
    vec[1] = 2;
    vec[2] = 3;

    sol[0] = 14;
    sol[1] = 32;
    sol[2] = 50;

    v = matXvec(m1, 3, 3, vec, 3);
    //printf("[0] = %5.19lf   [1] = %5.19lf    [2] = %5.19lf", v[0], v[1], v[2]);
    _assert(compareV(v, 3, sol, 3));

    freeVector(sol, 3);
    freeArray(m1, 3, 3);
    freeVector(vec, 3);
    freeVector(v, 3);

    return 0;
}

int vecXvecTrans_01()
{
    double *vec = vector(3);
    double **sol = array(3,3);
    double *v2 = vector(3);
    double **mat = array(3,3);

    vec[0] = 1;
    vec[1] = 2;
    vec[2] = 3;

    v2[0] = 1;
    v2[1] = 4;
    v2[2] = 5;

    sol[0][0] = 1;
    sol[0][1] = 4;
    sol[0][2] = 5;
    sol[1][0] = 2;
    sol[1][1] = 8;
    sol[1][2] = 10;
    sol[2][0] = 3;
    sol[2][1] = 12;
    sol[2][2] = 15;

    mat = vecXvecTrans(vec, 3, v2, 3);

    //printArray(mat, 3, 3);
    //printf("\r\n");
    //printArray(sol, 3, 3);

    _assert(compare(mat, 3, 3, sol, 3, 3));

    freeVector(v2, 3);
    freeArray(sol, 3, 3);
    freeArray(mat, 3, 3);
    freeVector(vec, 3);

    return 0;
}

int sign_01()
{
    double res = 1.1;
    double a = 1.1;
    double b = 1;

    _assert(res == sign(a,b));

    //_assert(fabs(res - sign(a,b)) < pow(10, -10));
    return 0;
}

int gmst_01()
{
    double res = -28.45275232595;
    double Mjd_UT1 = 49746.1163579368;

    double r = gmst(Mjd_UT1);
    //printf("r: %5.19lf", r);

    _assert(fabs(res - r) < pow(10, -3));

    return 0;
}

int angl_01()
{
    double *vec1 = vector(3);
    double *vec2 = vector(3);

    vec1[0] = 1;
    vec1[1] = 2;
    vec1[2] = 3;

    vec2[0] = 3;
    vec2[1] = 4;
    vec2[2] = 5;

    double sol = 0.1862;

    //printf("res: %5.19lf", angl(vec1, vec2));

    _assert(fabs(sol - angl(vec1, vec2)) < pow(10, -3));

    freeVector(vec1, 3);
    freeVector(vec2, 3);

    return 0;
}

int EccAnom_01()
{
    double M = 3.00;
    double e = 2.55;

    double sol = 3.1017;

    _assert(fabs(sol - EccAnom(M, e)) < pow(10, -3));

    return 0;
}

int Frac_01()
{
    double x = 3.4525;

    double sol = 0.4525;

    //printf("res = %f", Frac(3.4525));
    _assert(fabs(sol - Frac(3.4525)) < pow(10, -10));

    return 0;
}

int Mjday_01()
{
    double sol = 51387.926076;

    //printf("res = %f", Mjday(1999, 7, 28, 22, 13, 33));

    _assert(fabs(sol - Mjday(1999, 7, 28, 22, 13, 33)) < pow(10, -3));

    return 0;
}

int LTC_01()
{
    double **sol = array(3,3);
    double a = 2.44;
    double b = 3.66;

    sol[0][0] = -0.6454;
    sol[0][1] = -0.7638;
    sol[0][2] = 0.0000;
    sol[1][0] = -0.3785;
    sol[1][1] = 0.3198;
    sol[1][2] = -0.8686;
    sol[2][0] = 0.6635;
    sol[2][1] = -0.5606;
    sol[2][2] = -0.4955;

    double ** s = array(3,3);
    s = LTC(a, b);

    _assert(compare(sol, 3, 3, s, 3, 3));

    freeArray(sol, 3, 3);
    freeArray(s, 3, 3);

    return 0;
}

int Legendre_01()
{
    double **sol1 = array(2,4);
    double **sol2 = array(2,4);
    double **pnm = array(2,4);
    double **dpnm = array(2,4);

    sol1[0][0] = 1.0000;
    sol1[0][1] = 0;
    sol1[0][2] = 0;
    sol1[0][3] = 0;
    sol1[1][0] = -1.3108;
    sol1[1][1] = -1.1321;
    sol1[1][2] = 0;
    sol1[1][3] = 0;

    sol2[0][0] = 0;
    sol2[0][1] = 0;
    sol2[0][2] = 0;
    sol2[0][3] = 0;
    sol2[1][0] = -1.1321;
    sol2[1][1] = -1.3108;
    sol2[1][2] = 0;
    sol2[1][3] = 0;

    Legendre(1,3,3,pnm, dpnm);

    _assert(compare(sol1, 2, 4, pnm, 2, 4) && compare(sol2, 2, 4, dpnm, 2, 4));

    freeArray(sol1, 2, 4);
    freeArray(pnm, 2, 4);
    freeArray(sol2, 2, 4);
    freeArray(dpnm, 2, 4);

    return 0;
}

int AccelPointMass_01()
{
    double *sol = vector(3);
    double *vec = vector(3);
    double *r = vector(3);
    double *s = vector(3);

    sol[0] = 0.102730253176296;
    sol[1] = 0.33192909950257654;
    sol[2] = 0.30819075952888803;

    r[0] = 1;
    r[1] = 2;
    r[2] = 3;

    s[0] = 2;
    s[1] = 5;
    s[2] = 6;

    vec = AccelPointMass(r,3, s, 3, 12.44);

    //printVector(vec, 3);

    _assert(compareV(sol, 3, vec, 3));

    freeVector(sol, 3);
    freeVector(vec, 3);
    freeVector(r, 3);
    freeVector(s, 3);

    return 0;
}

int AzElPa_01()
{
    double *dAdsSol = vector(3);;
    double *dEdsSol = vector(3);;
    double *dAds = vector(3);
    double *dEds = vector(3);
    double *s = vector(3);;
    double Az;
    double El;

    dAdsSol[0] = 0.3999999999999999700;
    dAdsSol[1] = -0.1999999999999999800;
    dAdsSol[2] = 0.0000000000000000000;

    dEdsSol[0] = -0.0958314847499909810;
    dEdsSol[1] = -0.1916629694999819600;
    dEdsSol[2] = 0.1597191412499849900;

    s[0] = 1;
    s[1] = 2;
    s[2] = 3;

    AzElPa(s, Az, El, dAds, dEds);

    //printVector(dAds, 3);
    //printVector(dEds, 3);

    _assert(compareV(dAds, 3, dAdsSol, 3) && compareV(dEds, 3, dEdsSol, 3));

    freeVector(dAds, 3);
    freeVector(dAdsSol, 3);
    freeVector(dEds, 3);
    freeVector(dEdsSol, 3);
    freeVector(s, 3);

    return 0;
}

int gibbs_01()
{
    double *a = vector(3);
    double *b = vector(3);
    double *c = vector(3);

    double *v = vector(3);
    double *vSol = vector(3);

    double theta;
    double theta1;
    double copa;
    char *error = "";

    vSol[0] = 0;
    vSol[1] = 0;
    vSol[2] = 0;

    double thetaSol = 0;
    double theta1Sol = 0;
    double copaSol = 0;

    a[0] = 1;
    a[1] = 2;
    a[2] = 3;

    b[0] = 4;
    b[1] = 5;
    b[2] = 6;

    c[0] = 7;
    c[1] = 8;
    c[2] = 9;

    gibbs(a,b,c,v,&theta,&theta1,&copa, error);

    _assert(fabs(thetaSol - theta) < pow(10, -3) && fabs(theta1Sol - theta1) < pow(10, -3) && fabs(copaSol - copa) < pow(10, -3) && compareV(vSol, 3, v, 3));

    freeVector(a, 3);
    freeVector(b, 3);
    freeVector(c, 3);
    freeVector(v, 3);
    freeVector(vSol, 3);

    return 0;
}

int PoleMatrix_01()
{
    double **sol = array(3, 3);

    sol[0][0] = -0.8011;
    sol[0][1] = -0.0349;
    sol[0][2] = -0.5975;
    sol[1][0] = 0.0000;
    sol[1][1] = -0.9983;
    sol[1][2] = 0.0584;
    sol[2][0] = -0.5985;
    sol[2][1] = 0.0468;
    sol[2][2] = 0.7998;

    double ** mat = PoleMatrix(2.5, 3.2);

    //printArray(mat, 3, 3);

    _assert(compare(sol, 3, 3, mat, 3, 3));

    freeArray(sol, 3, 3);
    freeArray(mat, 3, 3);

    return 0;
}

int AccelHarmonic_01()
{
    double *sol, *r;

    sol[0] = -1.04985052621816e+27;
    sol[1] = 1.10199114014309e+26;
    sol[2] = 2.5243590332183e+26;

    r[0] = 1;
    r[1] = 2;
    r[2] = 3;

    double *a = AccelHarmonic(r, R_y(1), 3, 3, 4);

    _assert(compareV(sol, 3, a, 3));

    freeVector(a, 3);
    freeVector(sol, 3);
    freeVector(r, 3);

    return 0;
}

int gast_01()
{
    double sol = -3.319954;

    double res = gast(49746.1163579368);
    //printf("res: %lf", res);

    _assert(fabs(sol - res) < pow(10, -3));

    return 0;
}

int Geodetic_01()
{

    double lonsol = 0.982794;
    double latsol = 0.000000;
    double hsol = -6378132.694449;

    double lon, lat, h, *r;

    r = vector(3);

    r[0] = 1;
    r[1] = 2;
    r[2] = 3;

    Geodetic(r, &lon, &lat, &h);

    //printf("lon: %lf\n\r", lon);
    //printf("lat: %lf\n\r", lat);
    //printf("h: %lf\n\r", h);

    _assert(fabs(lonsol - lon) < pow(10, -3) && fabs(latsol - lat) < pow(10, -3) && fabs(hsol - h) < pow(10, -3));

    freeVector(r, 3);

    return 0;
}

int Cheb3D_01()
{
    double *sol = vector(3);

    sol[0] = 2.000000000;
    sol[1] = 4.000000000;
    sol[2] = 6.000000000;

    double *Cx = vector(3), *Cy = vector(3), *Cz = vector(3);

    Cx[0] = 1;
    Cx[1] = 2;
    Cx[2] = 3;

    Cy[0] = 3;
    Cy[1] = 4;
    Cy[2] = 5;

    Cz[0] = 2;
    Cz[1] = 6;
    Cz[2] = 8;

    double *res = Cheb3D(2, 3, 1, 3, Cx, Cy, Cz);

    //printVector(res, 3);

    _assert(compareV(sol, 3, res, 3));

    freeVector(sol, 3);
    freeVector(Cx, 3);
    freeVector(Cy, 3);
    freeVector(Cz, 3);
    freeVector(res, 3);

    return 0;
}

int doubler_01() {

    double *r2sol = vector(3);
    double *r3sol = vector(3);
    double f1sol, f2sol, q1sol, magr1sol, magr2sol, asol, deltae32sol;

    r2sol[0] = 6.04055555045488 ;
    r2sol[1] = 3.0506944380686;
    r2sol[2] = 2.06083332568232;

    r3sol[0] = 4.81026464007629;
    r3sol[1] = 2.78315958865862;
    r3sol[2] = 3.75605453724095;

    f1sol = 49746.1170627547;
    f2sol = 49746.1163576416;
    q1sol = 70351.632926957;
    magr1sol = 5.42898937992004;
    magr2sol = 7.07404282633607;
    asol = 7.46181838140688;
    deltae32sol = 0.285174707058671;


    double f1, f2, q1, magr1, magr2, a, deltae32;
    double *r2 = vector(3);
    double *r3 = vector(3);
    double *los1 = vector(3);
    double *los2 = vector(3);
    double *los3 = vector(3);
    double *rsite1 = vector(3);
    double *rsite2 = vector(3);
    double *rsite3 = vector(3);

    double Mjd_TT = 49746.1170623147;
    double Mjd0 = 49746.1097222222;
    double Mjd1 = 49746.1101504629;
    double Mjd2 = 49746.1112847221;
    double Mjd_UT1 = 49746.1163579368;
    double Mjd_UTC = 49746.1163541665;

    los1[0] = 1;
    los1[1] = 2;
    los1[2] = 3;

    los2[0] = 4;
    los2[1] = 5;
    los2[2] = 6;

    los3[0] = 7;
    los3[1] = 8;
    los3[2] = 9;

    rsite1[0] = 3;
    rsite1[1] = 2;
    rsite1[2] = 4;

    rsite2[0] = 6;
    rsite2[1] = 3;
    rsite2[2] = 2;

    rsite3[0] = 5;
    rsite3[1] = 3;
    rsite3[2] = 4;

    doubler(Mjd_TT, Mjd0 ,Mjd1 ,Mjd2 ,Mjd_UT1 ,Mjd_UTC, los1, los2, los3, rsite1, rsite2, rsite3, Mjd_TT, Mjd_UT1, 'a', r2, r3, &f1, &f2, &q1, &magr1, &magr2, &a, &deltae32);

    //printf("f1sol: %5.19lf    f2sol: %5.19lf     q1sol: %5.19lf     magr1sol: %5.19lf ", f1sol, f2sol, q1sol, magr1sol);

    //printf("f1: %5.19lf    f2: %5.19lf     q1: %5.19lf     magr1: %5.19lf ", f1, f2, q1, magr1);
    //printf("magr2: %5.19lf    a: %5.19lf     deltae32: %5.19lf ", magr2, a, deltae32);

    _assert(fabs(f1sol - f1) < pow(10, -3) && fabs(f2sol - f2) < pow(10, -3) && fabs(q1sol - q1) < pow(10, -3) && fabs(magr1sol - magr1) < pow(10, -3) && fabs(magr2sol - magr2) < pow(10, -3) && fabs(asol - a) < pow(10, -3) && fabs(deltae32sol - deltae32) < pow(10, -3));                                                                ;

    freeVector(r2sol, 3);
    freeVector(r3sol, 3);
    freeVector(r2, 3);
    freeVector(r3, 3);
    freeVector(los1, 3);
    freeVector(los2, 3);
    freeVector(los3, 3);
    freeVector(rsite1, 3);
    freeVector(rsite2, 3);
    freeVector(rsite3, 3);

    return 0;
}

int elements_01()
{
    double psol, asol, esol, isol, Omegasol, omegasol, Msol;

    psol = 1.35474011564823e-13;
    asol = 1.87082869338765;
    esol = 0.999999999999964;
    isol = 1.99133066207886;
    Omegasol = 3.6052402625906;
    omegasol = 5.21086941752228;
    Msol = 3.14159030993265;

    double p, a, e, i, Omega, omega, M;

    double *y = vector(6);

    y[0] = 1;
    y[1] = 2;
    y[2] = 3;
    y[3] = 4;
    y[4] = 5;
    y[5] = 6;

    elements(y, &p, &a, &e, &i, &Omega, &omega, &M);


    //printf("psol: %5.19lf    asol: %5.19lf     esol: %5.19lf     isol: %5.19lf \n\r", psol, asol, esol, isol);

    //printf("p: %5.19lf    a: %5.19lf     e: %5.19lf     i: %5.19lf \n\r", p, a, e, i);

    _assert(fabs(psol - p) < pow(10, -3) && fabs(asol - a) < pow(10, -3) && fabs(esol - e) < pow(10, -3) && fabs(isol - i) < pow(10, -3));

    freeVector(y, 6);
    return 0;
}


int GHAMatrix_01()
{
    double Mjd_UT1 = 49746.1163579368;
    double **sol = array(3, 3);

    sol[0][0] = -0.984135751656283;
    sol[0][1] = 0.177417085738447;
    sol[0][2] = 0;
    sol[1][0] = -0.177417085738447;
    sol[1][1] = -0.984135751656283;
    sol[1][2] = 0;
    sol[2][0] = 0;
    sol[2][1] = 0;
    sol[2][2] = 1;

    double **res = array(3, 3);

    res = GHAMatrix(Mjd_UT1);

    //printArray(res, 3, 3);
    //printf("\n\r");
    //printArray(sol, 3, 3);

    _assert(compare(res, 3, 3, sol, 3, 3));

    freeArray(res, 3, 3);
    freeArray(sol, 3, 3);

    return 0;
}

int hgibbs_01()
{
    double Mjd1 = 49746.1101504629;
    double Mjd2 = 49746.1112847221;
    double Mjd3 = 49746.1125347223;

    double *v2sol = vector(3);
    double thetasol, theta1sol, copasol;
    char* errorsol;

    v2sol[0] = -57228706960475.5;
    v2sol[1] = -99463188905415.5;
    v2sol[2] = -150134296058356;

    thetasol = 0.18623876586485;
    theta1sol = 0.311548173406822;
    copasol = 0.0338804674925181;
    errorsol = "angl > 1ø";


    double *v2 = vector(3);
    double theta, theta1, copa;
    char *error = "";
    double *vec1 = vector(3);
    double *vec2 = vector(3);
    double *vec3 = vector(3);

    vec1[0] = 1;
    vec1[1] = 2;
    vec1[2] = 3;

    vec2[0] = 3;
    vec2[1] = 4;
    vec2[2] = 5;

    vec3[0] = 1;
    vec3[1] = 4;
    vec3[2] = 6;

    hgibbs(vec1, vec2, vec3, Mjd1, Mjd2, Mjd3, v2, &theta, &theta1, &copa, error);

    //printf("thetasol: %5.19lf    theta1sol: %5.19lf     copasol: %5.19lf \n\r", thetasol, theta1sol, copasol);
    //printf("theta: %5.19lf    theta1: %5.19lf     copa: %5.19lf \n\r", theta, theta1, copa);

    _assert(fabs(thetasol - theta) < pow(10, -3) && fabs(theta1sol - theta1) < pow(10, -3) && fabs(copasol - copa) < pow(10, -3));

    freeVector(v2, 3);
    freeVector(v2sol, 3);
    freeVector(vec1, 3);
    freeVector(vec2, 3);
    freeVector(vec3, 3);

    return 0;
}


int PrecMatrix_01()
{
    double Mjd1 = 49746.1101504629;
    double Mjd2 = 49746.1112847221;

    double **sol = array(3, 3);

    sol[0][0] = 1;
    sol[0][1] = -6.94407432091004e-10;
    sol[0][2] = -3.01766571972735e-10;
    sol[1][0] = 6.94407432091004e-10;
    sol[1][1] = 1;
    sol[1][2] = -1.04774475167805e-19;
    sol[2][0] = 3.01766571972735e-10;
    sol[2][1] = -1.04774475166687e-19;
    sol[2][2] = 1;

    double **res = array(3, 3);

    res = PrecMatrix(Mjd1, Mjd2);

    //printArray(res, 3, 3);
    //printf("\n\r");
    //printArray(sol, 3, 3);

    _assert(compare(res, 3, 3, sol, 3, 3));

    freeArray(res, 3, 3);
    freeArray(sol, 3, 3);

    return 0;
}

int Accel_01()
{
    double Mjd1 = 49746.1101504629;

    double *sol = vector(6);

    sol[0] = 4;
    sol[1] = 5;
    sol[2] = 6;
    sol[3] = -1.95666551086025e+132;
    sol[4] = 1.03481228172367e+132;
    sol[5] = 1.60869310819873e+132;

    double *Y = vector(6);

    Y[0] = 1;
    Y[1] = 2;
    Y[2] = 3;
    Y[3] = 4;
    Y[4] = 5;
    Y[5] = 6;

    double *dY = vector(6);

   Accel(Mjd1, Y, dY);

    printVector(sol, 3);
    printf("\n\r");
    printVector(dY, 3);

    _assert(compareV(dY, 3, sol, 3));

    freeVector(Y, 3);
    freeVector(sol, 3);

    return 0;
}

int IERS_01()
{
    double Mjd_UTC = 49746.1163541665;
    char interp = 'l';

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;

    double x_polesol, y_polesol, UT1_UTCsol, LODsol, dpsisol, depssol, dx_polesol, dy_polesol, TAI_UTCsol;

    x_polesol = -5.5937872420407e-07;
    y_polesol = 2.33559834147197e-06;
    UT1_UTCsol = 0.325747632958709;
    LODsol = 0.00272698971874332;
    dpsisol = -1.16882953161744e-07;
    depssol = -2.4783506198648e-08;
    dx_polesol = -8.43027359626024e-10;
    dy_polesol = -1.56811369105037e-09;
    TAI_UTCsol = 29;

    IERS(Mjd_UTC, interp, &x_pole, &y_pole, &UT1_UTC, &LOD, &dpsi, &deps, &dx_pole, &dy_pole, &TAI_UTC);


    printf("x_polesol: %5.19lf    y_polesol: %5.19lf    UT1_UTCsol: %5.19lf    LODsol: %5.19lf    \n\r", x_polesol, y_polesol, UT1_UTCsol, LODsol);
    printf("x_pole: %5.19lf    y_pole: %5.19lf    UT1_UTC: %5.19lf    LOD: %5.19lf    \n\r", x_pole, y_pole, UT1_UTC, LOD);

    _assert(fabs(x_polesol - x_pole) < pow(10, -3) && fabs(y_polesol - y_pole) < pow(10, -3) && fabs(UT1_UTCsol - UT1_UTC) < pow(10, -3) && fabs(LODsol - LOD) < pow(10, -3));

    return 0;
}

int G_AccelHarmonic_01()
{
    double **Gsol = array(3,3);

    Gsol[0][0] = 2.02204885013657e-06;
    Gsol[0][1] = 5.62339573129123e-07;
    Gsol[0][2] = 5.56098012083339e-09;
    Gsol[1][0] = 5.62339572462989e-07;
    Gsol[1][1] = -9.58402524675606e-07;
    Gsol[1][2] = 1.01432728882855e-09;
    Gsol[2][0] = 5.56098022838625e-09;
    Gsol[2][1] = 1.01432672504342e-09 ;
    Gsol[2][2] = -1.06364632509841e-06;

    double *r = vector(3);

    r[0] = 7101597.83995656;
    r[1] = 1295244.79268408;
    r[2] = 12755.6074812165;

    double **U = array(3,3);

    U[0][0] = -0.984125607526048;
    U[0][1] = 0.177473346200175;
    U[0][2] = 0;
    U[1][0] = -0.177473346200175;
    U[1][1] = -0.984125607526048;
    U[1][2] = 0;
    U[2][0] = 0;
    U[2][1] = 0;
    U[2][2] = 1;

    double **G = array(3,3);

    //G = AccelHarmonic(r, U, 3, 2, 5);  //comento esto porque me da un error de punteros que no consigo encontrar

    printArray(Gsol, 3, 3);
    printArray(G, 3, 3);

    _assert(compare(G, 3, 3, Gsol, 3, 3));

    freeArray(G, 3, 3);
    freeArray(Gsol, 3, 3);
    freeArray(U, 3, 3);
    freeVector(r, 3);

    return 0;
}


int JPL_Eph_DE430_01()
{
    double Mjd1 = 49746.1101504629;

    double *r_MercurySol = vector(3);
    double *r_VenusSol = vector(3);
    double *r_EarthSol = vector(3);

    r_MercurySol[0] = 83780898455.9384;
    r_MercurySol[1] = -65292523959.1825;
    r_MercurySol[2] = -23392735743.3518;

    r_VenusSol[0] = -15233538320.351;
    r_VenusSol[1] = -110132631321.853;
    r_VenusSol[2] = -41020690575.6379;

    r_EarthSol[0] = -92467183721.5114;
    r_EarthSol[1] = 106397659934.095;
    r_EarthSol[2] = 46131328590.3046;


    double *r_Mercury = vector(3);
    double *r_Venus = vector(3);
    double *r_Earth = vector(3);
    double *r_Mars = vector(3);
    double *r_Jupiter = vector(3);
    double *r_Saturn = vector(3);
    double *r_Uranus = vector(3);
    double *r_Neptune = vector(3);
    double *r_Pluto = vector(3);
    double *r_Moon = vector(3);
    double *r_Sun = vector(3);

    JPL_Eph_DE430(Mjd1, r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun);

    printVector(r_MercurySol, 3);
    printVector(r_Mercury, 3);
    printf("\r\n");
    printVector(r_VenusSol, 3);
    printVector(r_Venus, 3);
    printf("\r\n");
    printVector(r_EarthSol, 3);
    printVector(r_Earth, 3);

    _assert(compareV(r_MercurySol, 3, r_Mercury, 3) && compareV(r_VenusSol, 3, r_Venus, 3) && compareV(r_EarthSol, 3, r_Earth, 3));

    freeVector(r_MercurySol, 3);
    freeVector(r_Mercury, 3);
    freeVector(r_VenusSol, 3);
    freeVector(r_Venus, 3);
    freeVector(r_EarthSol, 3);
    freeVector(r_Mars, 3);
    freeVector(r_Jupiter, 3);
    freeVector(r_Saturn, 3);
    freeVector(r_Uranus, 3);
    freeVector(r_Neptune, 3);
    freeVector(r_Pluto, 3);
    freeVector(r_Moon, 3);
    freeVector(r_Sun, 3);

    return 0;
}


int all_tests()
{
    _verify(position_01);
    _verify(dot_01);
    _verify(norma_01);
    _verify(trasp_01);
    _verify(cross_01);
    _verify(prod_01);
    _verify(eye_01);
    _verify(sum_01);
    _verify(matXvec_01);
    _verify(vecXvecTrans_01);
    //_verify(inv_01);
    _verify(R_x_01);
    _verify(R_y_01);
    _verify(R_z_01);
    _verify(MeanObliquity_01);
    _verify(Mjday_TBD_01);
    _verify(NutAngles_01);
    _verify(EqnEquinox_01);
    _verify(NutMatrix_01);
    _verify(unit_01);
    _verify(timediff_01);
    //_verify(MeasUpdate_01);
    _verify(sign_01);
    _verify(gmst_01);
    _verify(angl_01);
    _verify(EccAnom_01);
    _verify(Frac_01);
    _verify(Mjday_01);
    _verify(LTC_01);
    //_verify(Legendre_01);
    _verify(AccelPointMass_01);
    _verify(AzElPa_01);
    _verify(gibbs_01);
    _verify(PoleMatrix_01);
    //_verify(AccelHarmonic_01);
    _verify(gast_01);
    _verify(Geodetic_01);
    _verify(Cheb3D_01);
    _verify(doubler_01);
    _verify(elements_01);
    _verify(GHAMatrix_01);
    _verify(hgibbs_01);
    _verify(PrecMatrix_01);
    //_verify(Accel_01);
    //_verify(IERS_01);
    //_verify(G_AccelHarmonic_01);
    //_verify(JPL_Eph_DE430_01);

    return 0;
}

int main()
{
    int result = all_tests();

    double *v;

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
*/