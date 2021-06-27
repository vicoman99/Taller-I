/** @file EccAnom.h
 *  @brief Specification of EccAnom
 *
 *  @author Victor Coman
 *  @date May 2021
*/

#include "sat_const.h"
#include "EccAnom.h"
#include <stdio.h>

double EccAnom(double M, double e)
{
    double E;
    int maxit = 15;
    int i = 1;

    // Valor de comienzo
    M = fmod(M, 2.0 * M_PI);

    if (e < 0.8)
    {
        E = M;
    }
    else
    {
        E = M_PI;
    }
    double f = E - e * sin(E) - M;
    E = E - f / (1.0 - e * cos(E));

    // iteraciones
    while (fabs(f) > 1e2 * EPSILON)
    {
        f = E - e * sin(E) - M;
        E = E - f / (1.0 - e * cos(E));
        i = i + 1;
        if (i == maxit)
        {
            printf("Problemas de convergencia en EccAnom");
        }
    }
    return E;
}