/** @file PoleMatrix.c
 *  @brief Implementation of PoleMatrix
 *  @author Victor Coman
 *  @date May 2021
*/

#include <stdio.h>
#include <math.h>
#include "arrays.h"
#include "sat_const.h"
#include "GHAMatrix.h"
#include "R_y.h"
#include "R_x.h"
#include "PoleMatrix.h"

double **PoleMatrix(double xp, double yp)
{
    double **PoleMat = prod(R_y(-xp), 3, 3, R_x(-yp), 3, 3);
    return PoleMat;
}
