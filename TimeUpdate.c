/** @file TimeUpdate.h
 *  @brief Implementation of TimeUpdate
 *  @author Victor Coman
 *  @date May 2021
*/

#include <stdio.h>
#include <math.h>
#include "arrays.h"
#include "sat_const.h"

double **TimeUpdate(double **P, int np, double **Phi, int nphi, int Qdt)
{
    P = prod(Phi, nphi, nphi, prod(P, np, np, trasp(Phi, nphi), nphi, nphi), nphi, nphi);
    for (int i = 0; i < np; i++)
    {
        for (int j = 0; j < np; j++)
        {
            P[i][j] = P[i][j] + Qdt;
        }
    }
    return P;
}
