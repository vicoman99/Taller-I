/** @file R_z.c
 *  @brief Implementation of R_z
 * 
 *  @author Victor Coman
 *  @date May 2021
 */

#include "R_z.h"
#include <math.h>
#include "arrays.h"

double **R_z(double angle)
{
    double C = cos(angle);
    double S = sin(angle);
    double **rotmat = array(3, 3);

    rotmat[0][0] = C;
    rotmat[1][0] = -1.0 * S;
    rotmat[2][0] = 0.0;
    rotmat[0][1] = S;
    rotmat[1][1] = C;
    rotmat[1][2] = 0.0;
    rotmat[2][0] = 0.0;
    rotmat[2][1] = 0.0;
    rotmat[2][2] = 1.0;

    return rotmat;
}