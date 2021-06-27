/** @file EccAnom.h
 *  @brief Specification of EccAnom
 *
 *  @author Victor Coman
 *  @date May 2021
 */

#ifndef ECCANOM_H
#define ECCANOM_H
#include <math.h>

/** Computes the eccentric anomaly for elliptic orbits
 *  @param [in] M         Mean anomaly in [rad]
 *  @param [in] e         Eccentricity of the orbit [0,1]
 *  @return Eccentric anomaly in [rad]
 */
double EccAnom(double M, double e);


#endif