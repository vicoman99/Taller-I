/** @file GHAMatrix.h
 *  @brief Implementation of GHAMatrix
 *  @author Victor Coman
 *  @date May 2021
*/

#include <stdio.h>
#include <math.h>
#include "arrays.h"
#include "sat_const.h"
#include "GHAMatrix.h"
#include "R_z.h"
#include "gast.h"

double **GHAMatrix(double Mjd_UT1)
{
    double **GHAmat = R_z(gast(Mjd_UT1));
    return GHAmat;
}
