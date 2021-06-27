
/** @file R_x.c
 *  @brief Implementation of R_x
 * 
 *  @author Victor Coman
 *  @date May 2021
 */

#include <math.h>
#include "R_x.h"
#include "arrays.h"

double **R_x(double angle)
{
    double C = cos(angle);
    double S = sin(angle);
    double **rotmat = array(3, 3);

    rotmat[0][0] = 1.0;
    rotmat[1][0] = 0.0;
    rotmat[2][0] = 0.0;
    rotmat[0][1] = 0.0;
    rotmat[1][1] = C;
    rotmat[1][2] = S;
    rotmat[2][0] = 0.0;
    rotmat[2][1] = -1.0 * S;
    rotmat[2][2] = C;

    return rotmat;
}
