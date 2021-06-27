/** @file LTC.h
 *  @brief Implementation of LTC
 *  @author Victor Coman
 *  @date May 2021
*/

#include <stdio.h>
#include <math.h>
#include "arrays.h"
#include "sat_const.h"
#include "LTC.h"
#include "R_y.h"
#include "R_z.h"

double **LTC(double lon, double lat)
{
    double **M =  array(3,3), aux;

    M = prod(R_y(-1.0 * lat), 3, 3, R_z(lon), 3, 3);

    for (int i = 0; i < 3; i++)
    {
        aux = M[0][i];
        M[0][i] = M[1][i];
        M[1][i] = M[2][i];
        M[2][i] = aux;
    }

    return M;
}
