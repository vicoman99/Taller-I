/** @file globales.h
 *  @brief Specification of globales
 *
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef GLOBALES_H
#define GLOBALES_H

typedef struct
{
    double Mjd_UTC, Mjd_TT;
    int n, m, sun, moon, planets;
} Param;

double **PC, **Cnm, **Snm, **eopdata;
int n_eqn;
Param AuxParam;

#endif