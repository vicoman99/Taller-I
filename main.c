/** @file main.c
 *  @brief Initial Orbit Determination using Gauss and Extended Kalman Filter methods
 *  @author Victor Coman
 *  @date May 2021
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globales.h"
#include "arrays.h"
#include "sat_const.h"
#include "Mjday.h"
#include "position.h"
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
#include "R_x.h"
#include "R_y.h"
#include "R_z.h"
#include "gmst.h"
#include "TimeUpdate.h"
#include "VarEqn.h"
#include "MeasUpdate.H"
#include "Anglesg.h"
#include "Anglesdr.h"

int main()
{
    extern double **PC, **Cnm, **Snm, **eopdata;
    extern int n_eqn;
    extern Param AuxParam;

    double **P, **obs, aux1, aux2, ss, az, el, Dist, sigma_range, sigma_az, sigma_el, lat, lon, alt, *Rs, *Y0_apr, Mjd0, Mjd_UTC;
    FILE *fp;
    int f, c, n, m, nobs, i, Y, M, D, hh, mm;
    char line[55], y[5], mo[3], d[3], h[3], mi[3], s[7], a[9], e[9], di[10];

    PC = array(2285, 1020);
    fp = fopen("../Data/DE430Coeff.txt", "r");

    if (fp == NULL)
    {
        printf("Fail open DE430Coeff.txt file \n");
        exit(EXIT_FAILURE);
    }

    for (n = 0; n <= 180; ++n)
    {
        for (m = 0; m <= n; ++m)
        {
            fscanf(fp, "%d%d%lf%lf%lf%lf", &f, &c, &Cnm[n + 1][m + 1], &Snm[n + 1][m + 1], &aux1, aux2);
        }
    }
    fclose(fp);

    eopdata = array(13, 21413);
    fp = fopen("../Data/eop19620101.txt", "r");

    if (fp == NULL)
    {
        printf("Fail open eop19620101.txt file \n");
        exit(EXIT_FAILURE);
    }

    //  ----------------------------------------------------------------------------------------------------
    // |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
    // |(0h UTC)           "         "          s          s          "        "          "         "     s
    //  ----------------------------------------------------------------------------------------------------

    for (f = 0; f <= 21413; ++f)
    {
        fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &eopdata[0][f], &eopdata[1][f], &eopdata[2][f], &eopdata[3][f], &eopdata[4][f], &eopdata[5][f], &eopdata[6][f], &eopdata[7][f], &eopdata[8][f], &eopdata[9][f], &eopdata[10][f], &eopdata[11][f], &eopdata[12][f]);
    }
    fclose(fp);

    nobs = 46;
    obs = array(nobs, 4);

    //read observations
    fp = fopen("../Data/GEOS3.txt", "r");

    if (fp == NULL)
    {
        printf("Fail open GEOS3.txt file \n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i <= nobs; ++i)
    {
        fgets(line, sizeof(line) + 2, fp);
        //print("%s",line);
        strncpy(y, &line[0], 4);
        y[4] = '\0';
        Y = atoi(y);
        //print("%d\n",atoi(y));
        strncpy(mo, &line[5], 2);
        mo[2] = '\0';
        M = atoi(mo);
        //print("%d\n",atoi(mo));
        strncpy(d, &line[8], 2);
        d[2] = '\0';
        D = atoi(d);
        //print("%d\n",atoi(d));
        strncpy(h, &line[12], 2);
        h[2] = '\0';
        hh = atoi(h);
        //print("%d\n",atoi(h));
        strncpy(mi, &line[15], 2);
        mi[2] = '\0';
        mm = atoi(mi);
        //print("%d\n",atoi(mi));
        strncpy(s, &line[18], 6);
        s[6] = '\0';
        ss = atof(s);
        //print("%lf\n",atof(s));
        strncpy(a, &line[25], 8);
        a[8] = '\0';
        az = atof(a);
        //print("%lf\n",atof(a));
        strncpy(e, &line[35], 8);
        e[8] = '\0';
        el = atof(e);
        //print("%lf\n",atof(e));
        strncpy(di, &line[44], 9);
        di[9] = '\0';
        Dist = atof(di);
        //print("%lf\n",atof(di));
        obs[i][0] = Mjday(Y, M, D, hh, mm, ss);
        obs[i][1] = Rad * az;
        obs[i][2] = Rad * el;
        obs[i][3] = 1e3 * Dist;
    }

    fclose(fp);

    sigma_range = 92.5;      //[m]
    sigma_az = 0.0224 * Rad; //[rad]
    sigma_el = 0.0139 * Rad; //[rad]

    //Kaena point station
    lat = Rad * 21.5748;     //[rad]
    lon = Rad * (-158.2706); //[rad]
    alt = 300.20;            //m

    Rs = vector(3);
    Rs = position(lon, lat, alt);

    double Mjd1 = obs[0][0];
    double Mjd2 = obs[8][0];
    double Mjd3 = obs[17][0];
    double *r2 = vector(3);
    double *v2 = vector(3);

    //JF[r2,v2] = anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);
    anglesg(r2, v2, obs[0][1], obs[8][1], obs[17][1], obs[0][2], obs[8][2], obs[17][2], Mjd1, Mjd2, Mjd3, Rs, Rs, Rs);

    // [r2,v2] = anglesdr(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);
    anglesdr(r2, v2, obs[0][1], obs[8][1], obs[17][1], obs[0][2], obs[8][2], obs[17][2], Mjd1, Mjd2, Mjd3, Rs, Rs, Rs);

    Y0_apr = vector(6);
    Y0_apr[0] = 6221397.62857869;
    Y0_apr[1] = 2867713.77965741;
    Y0_apr[2] = 3006155.9850995;
    Y0_apr[3] = 4645.0472516175;
    Y0_apr[4] = -2752.21591588182;
    Y0_apr[5] = -7507.99940986939;

    Mjd0 = Mjday(1995, 1, 29, 02, 38, 0);
    Mjd_UTC = obs[8][0];

    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;

    n_eqn = 6;
    double t = 0.0;
    double relerr = 1e-13;
    double abserr = 1e-6;
    int iflag = 1;
    double *work = vector(100 + 21 * n_eqn);
    double YY[6];
    int *iwork;

    for (i = 0; i < 6; i++)
    {
        YY[i] = Y0_apr[i];
    }

    ode(Accel, n_eqn, YY, &t, -(obs[8][0] - Mjd0) * 86400.0, relerr, abserr, &iflag, work, iwork);

    P = array(6, 6);

    for (int i = 0; i < 3; i++)
    {
        P[i][i] = 1e8;
        P[i + 3][i + 3] = 1e3;
    }

    double **LT = LTC(lon, lat);

    double *yPhi = vector(42);
    double **Phi = array(6, 6);

    //Measurement loop
    t = 0;

    double t_old, *Y_old;

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC, Mjd_TT, Mjd_UT1, theta, **U, *r;
    double Azim, Elev, *dAds, *dEds, *dAdY, *dEdY, *dDdY;

    dAds = vector(3);
    dEds = vector(3);
    dAdY = vector(6);
    dEdY = vector(6);
    dDdY = vector(6);

    for (i = 0; i < nobs; i++)
    {
        //Previous step
        t_old = t;
        Y_old = YY;

        //Time increment and propagation
        Mjd_UTC = obs[i][1];            // Modified Julian Date
        t = (Mjd_UTC - Mjd0) * 86400.0; // Time since epoch [s]

        IERS(Mjd0, 'l', &x_pole, &y_pole, &UT1_UTC, &LOD, &dpsi, &deps, &dx_pole, &dy_pole, &TAI_UTC);
        timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);

        Mjd_TT = Mjd_UTC + TT_UTC / 86400;
        Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;

        AuxParam.Mjd_UTC = Mjd_UTC;
        AuxParam.Mjd_TT = Mjd_TT;

        for (int ii = 0; i < 6; i++)
        {
            yPhi[ii] = Y_old[ii];
            for (int j = 0; j < 6; j++)
            {
                if (ii == j)
                {
                    yPhi[6 * j + ii] = 1;
                }
                else
                {
                    yPhi[6 * j + ii] = 0;
                }
            }
        }

        n_eqn = 42;

        //yPhi = DEInteg (@VarEqn,0,t-t_old,1e-13,1e-6,42,yPhi);
        ode(VarEqn, n_eqn, YY, &t, t-t_old, relerr, abserr, &iflag, work, iwork);

        // Extract state transition matrices
        for (int j = 0; j < 6; j++)
        {
            for (int i = 0; i < 6; i++)
            {
                Phi[i][j] = yPhi[6 * j + i];
            }
        }

        n_eqn = 6;

        //Y = DEInteg (@Accel,0,t-t_old,1e-13,1e-6,6,Y_old);
        ode(Accel, n_eqn, YY, &t, t-t_old, relerr, abserr, &iflag, work, iwork);

        // Topocentric coordinates
        theta = gmst(Mjd_UT1); // Earth rotation
        U = R_z(theta);
        r[0] = YY[0];
        r[1] = YY[1];
        r[2] = YY[2];
        double *ss = matXvec(LT, 3, 3, sumV(matXvec(U, 3, 3, r, 3), 3, vec_x_esc(Rs, 3, -1), 3), 3); // Topocentric position [m]

        // Time update
        P = TimeUpdate(P, 6, Phi, 6, 0);

        // Azimuth and partials
        AzElPa(ss, Azim, Elev, dAds, dEds); //Azimuth, Elevation

        double **parte1 = prod(LT, 3, 3, U, 3, 3);
        double *parte2 = vecXmat(dAds, 3, parte1, 3, 3);
        for (int i = 0; i < 3; i++)
        {
            dAdY[i] = parte2[i];
            dAdY[i + 3] = 0;
        }

        double *K;
        // Measurement update
        MeasUpdate(YY, 6, obs[i][2], Azim, sigma_az, dAdY, 6, P, 6, 6, K);

        // Elevation and partials
        r[0] = YY[0];
        r[1] = YY[1];
        r[2] = YY[2];

        ss = matXvec(LT, 3, 3, sumV(matXvec(U, 3, 3, r, 3), 3, vec_x_esc(Rs, 3, -1), 3), 3); // Topocentric position [m]
        AzElPa(ss, Azim, Elev, dAds, dEds);                                                  //Azimuth, Elevation

        parte1 = prod(LT, 3, 3, U, 3, 3);
        parte2 = vecXmat(dEds, 3, parte1, 3, 3);
        for (int i = 0; i < 3; i++)
        {
            dEdY[i] = parte2[i];
            dEdY[i + 3] = 0;
        }

        // Measurement update
        MeasUpdate(YY,6, obs[i][3], Elev, sigma_el, dEdY, 6, P, 6, 6, K);

        // Range and partials
        r[0] = YY[0];
        r[1] = YY[1];
        r[2] = YY[2];

        ss = matXvec(LT, 3, 3, sumV(matXvec(U, 3, 3, r, 3), 3, vec_x_esc(Rs, 3, -1), 3), 3); // Topocentric position [m]
        Dist = norma(ss, 3);
        double *dDds = vec_x_esc(ss, 3, 1 / Dist);

        parte1 = prod(LT, 3, 3, U, 3, 3);
        parte2 = vecXmat(dDds, 3, parte1, 3, 3);
        for (int i = 0; i < 3; i++)
        {
            dDdY[i] = parte2[i];
            dDdY[i + 3] = 0;
        }

        // Measurement update
        MeasUpdate(YY,6,  obs[i][4], Dist, sigma_range, dDdY, 6, P, 6, 6, K);
    }

    IERS(Mjd0, 'l', &x_pole, &y_pole, &UT1_UTC, &LOD, &dpsi, &deps, &dx_pole, &dy_pole, &TAI_UTC);
    timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);

    Mjd_TT = Mjd_UTC + TT_UTC / 86400;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;


    double *r_Mercury, *r_Venus, *r_Earth, *r_Mars, *r_Jupiter, *r_Saturn, *r_Uranus, *r_Neptune, *r_Pluto, *r_Moon, *r_Sun;
    
    JPL_Eph_DE430(Mjd0, r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun);

    double *vv;
    vv = vector(21413);

    for (int i = 0; i < 21413; i++)
    {
        vv[i] = eopdata[3][i];
    }

    find2(vv, 21413, Mjd0);

    double dYY[6];
    Accel(1.0, YY, dYY);

    double *v3, **m3;
    v3 = vector(3);
    v3[0] = 72000;
    v3[1] = 0;
    v3[2] = 0;
    m3 = R_y(1.0);

    r_Sun[0] = 92293598495.9694;
    r_Sun[1] = -105378575309.884;
    r_Sun[2] = -45687832459.835;


    //Y0 = DEInteg (Accel,0,-(obs[46][1]-obs[1][1])*86400.0,1e-13,1e-6,6,YY);
    ode(Accel, n_eqn, YY, &t, -(obs[45][0] - obs[0][0]) * 86400.0, relerr, abserr, &iflag, work, iwork);

    double *Y_true = vector(6);
    Y_true[0] = 5753.173e3;
    Y_true[1] = 2673.361e3;
    Y_true[2] = 3440.304e3;
    Y_true[3] = 4.324207e3;
    Y_true[4] = -1.924299e3;
    Y_true[5] = -5.728216e3;


    printf("\n Error of Position Estimation \n");
    printf("dX%10.1f [m] \n",YY[1]-Y_true[1]);
    printf("dY%10.1f [m] \n",YY[2]-Y_true[2]);
    printf("dZ%10.1f [m] \n",YY[3]-Y_true[3]);
    printf("\n Error of Velocity Estimation \n");
    printf("dVx%8.1f [m/s] \n",YY[4]-Y_true[4]);
    printf("'dVy%8.1f [m/s] \n",YY[5]-Y_true[5]);
    printf("dVz%8.1f [m/s] \n",YY[6]-Y_true[6]);


    freeArray(PC, 2285, 1020);
    freeArray(Cnm, 182, 182);
    freeArray(Snm, 182, 182);
    freeArray(eopdata, 13, 21413);
    freeArray(obs, nobs, 4);

    return 0;
}