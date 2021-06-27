/** @file PrecMatrix.h
 *  @brief Implementation of PrecMatrix
 *  @author Victor Coman
 *  @date May 2021
*/

#include <stdio.h>
#include <math.h>
#include "arrays.h"
#include "sat_const.h"
#include "R_y.h"
#include "R_z.h"

double **PrecMatrix(double Mjd_1, double Mjd_2)
{
    double T, dT, zeta, z, theta;
    double **PrecMat = array(3,3);

    T = (Mjd_1 - MJD_J2000) / 36525.0;
    dT = (Mjd_2 - Mjd_1) / 36525.0;

    // Precession angles
    zeta = ((2306.2181 + (1.39656 - 0.000139 * T) * T) + ((0.30188 - 0.000344 * T) + 0.017998 * dT) * dT) * dT / (Arcs);
    z = zeta + ((0.79280 + 0.000411 * T) + 0.000205 * dT) * dT * dT / (Arcs);
    theta = ((2004.3109 - (0.85330 + 0.000217 * T) * T) - ((0.42665 + 0.000217 * T) + 0.041833 * dT) * dT) * dT / (Arcs);

    //Precession matrix
    PrecMat = prod(R_z(-z), 3, 3, prod(R_y(theta), 3, 3, R_z(-zeta), 3, 3), 3, 3);

    return PrecMat;
}
