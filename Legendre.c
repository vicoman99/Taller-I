/** @file Legendre.h
 *  @brief Implementation of Legendre
 *  @author Victor Coman
 *  @date May 2021
*/
#include <stdio.h>
#include <math.h>
#include "arrays.h"
#include "sat_const.h"
#include "Legendre.h"

void Legendre(int n, int m, double fi, double **pnm, double **dpnm)
{
    pnm[0][0] = 1;
    dpnm[0][0] = 0;
    pnm[2][2] = sqrt(3) * cos(fi);
    dpnm[2][2] = -sqrt(3) * sin(fi);

    //diagonal coefficients
    for (int i = 0; i < 2; i++)
    {
        pnm[i + 1][i + 1] = sqrt((2 * i + 1) / (2 * i)) * cos(fi) * pnm[i][i];
    }

    for (int i = 0; i < 2; i++)
    {
        dpnm[i + 1][i + 1] = sqrt((2 * i + 1) / (2 * i)) * ((cos(fi) * dpnm[i][i]) - (sin(fi) * pnm[i][i]));
    }

    //horizontal first step coefficients
    for (int i = 0; i < 1; i++)
    {
        pnm[i + 1][i] = sqrt(2 * i + 1) * sin(fi) * pnm[i][i];
    }

    for (int i = 0; i < 1; i++)
    {
        dpnm[i + 1][i] = sqrt(2 * i + 1) * ((cos(fi) * pnm[i][i]) + (sin(fi) * dpnm[i][i]));
    }

    //horizontal second step coefficients
    int j = 0;
    int k = 2;

    while (1)
    {
        for (int i = 0; i < k; i++)
        {
            pnm[i + 1][j + 1] = sqrt((2 * i + 1) / ((i - j) * (i + j))) * ((sqrt(2 * i - 1) * sin(fi) * pnm[i][j + 1]) - (sqrt(((i + j - 1) * (i - j - 1)) / (2 * i - 3)) * pnm[i - 1][j + 1]));
        }

        j = j + 1;
        k = k + 1;

        if (j > m)
        {
            break;
        }
    }

    j = 0;
    k = 2;

    while (1)
    {
        for (int i = 0; i < k; i++)
        {
            dpnm[i + 1][j + 1] = sqrt((2 * i + 1) / ((i - j) * (i + j))) * ((sqrt(2 * i - 1) * sin(fi) * dpnm[i][j + 1]) + (sqrt(2 * i - 1) * cos(fi) * pnm[i][j + 1]) - (sqrt(((i + j - 1) * (i - j - 1)) / (2 * i - 3)) * dpnm[i - 1][j + 1]));
        }

        j = j + 1;
        k = k + 1;

        if (j > m)
        {
            break;
        }
    }

}
