/** @file G_AccelHarmonic.h
 *  @brief Implementation of G_AccelHarmonic
 *  @author Victor Coman
 *  @date May 2021
*/
#include <stdio.h>
#include <math.h>
#include "arrays.h"
#include "sat_const.h"
#include "G_AccelHarmonic.h"
#include "AccelHarmonic.h"

double **G_AccelHarmonic(double *r, double **U, int nU, int n_max, int m_max)
{
    double d;
    double **G, *dr, *da, *da1, *da2;

    d = 1.0; // Position increment [m]
    G = array(3, 3);
    dr = vector(3);
    da = vector(3);
    da1 = vector(3);
    da2 = vector(3);

    // Gradient
    for (int i = 0; i < 3; i++)
    {
        // Set offset in i-th component of the position vector
        for (int j = 0; j < 3; j++)
        {
            dr[j] = 0.0; //dr[:] = 0.0;
        }

        dr[i] = d / 2;

        // Acceleration difference
        da1 = AccelHarmonic(sumV(r, 3, dr, 3), U, nU, n_max, m_max);
        da1 = AccelHarmonic(sumV(r, 3, vec_x_esc(dr, 3, -1), 3), U, nU, n_max, m_max);
        da = sumV(da1, 3, vec_x_esc(da2, 3, -1), 3);

        // Derivative with respect to i-th axis
        for (int j = 0; j < 3; j++)
        {
            G[j][i] = da[j] / d; //G(:, i) = da / d;
        }
    }

    freeVector(da, 3);
    freeVector(da1, 3);
    freeVector(da2, 3);
    freeVector(dr, 3);

    return G;
}
