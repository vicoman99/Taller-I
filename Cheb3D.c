/** @file Cheb3D.h
 *  @brief Implementation of Cheb3D
 *  @author Victor Coman
 *  @date May 2021
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arrays.h"
#include "Cheb3D.h"

double *Cheb3D(double t, double N, double Ta, double Tb, double *Cx, double *Cy, double *Cz)
{
    double tau;
    double *f1, *f2, *old_f1, *ChebApp;

    // Check validity
    if ((t < Ta) || (Tb < t))
    {
        printf("error");
        exit(EXIT_FAILURE);
    }

    // Clenshaw algorithm
    tau = (2 * t - Ta - Tb) / (Tb - Ta);

    f1 = zerosV(3);
    f2 = zerosV(3);
    old_f1 = zerosV(3);
    ChebApp = zerosV(3);

    for (int i = N; i >= 2; i--)
    {
        for (int j = 0; j < 3; j++)
        {
            old_f1[j] = f1[j];
            f1[j] = 2 * tau * f1[j] - f2[j];
            f2 = old_f1;
        }
        f1[0] += Cx[i];
        f1[1] += Cy[i];
        f1[2] += Cz[i];
    }

    for (int j = 0; j < 3; j++)
    {
        ChebApp[j] = 2 * tau * f1[j] - f2[j];
    }
    ChebApp[0] += Cx[1];
    ChebApp[1] += Cy[1];
    ChebApp[2] += Cz[1];

    return ChebApp;
}
