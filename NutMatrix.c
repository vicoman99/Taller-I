/** @file NutMatrix.h
 *  @brief Implementation of NutMatrix
 *
 *  @author Victor Coman
 *  @date May 2021
*/

#include "NutMatrix.h"
#include <math.h>
#include "NutAngles.h"
#include "R_x.h"
#include "R_z.h"
#include "arrays.h"
#include "MeanObliquity.h"

double **NutMatrix(double Mjd_TT)
{
    double eps, dpsi, deps;
    // Mean obliquity of the ecliptic
    eps = MeanObliquity(Mjd_TT);

    // Nutation in longitude and obliquity
    NutAngles(Mjd_TT, &dpsi, &deps);

    // Transformation from mean to true equator and equinox
    double **NutMat = prod( R_x(-eps - deps), 3, 3, prod(R_z(-dpsi), 3, 3, R_x(eps), 3, 3), 3, 3);
    return NutMat;
}
