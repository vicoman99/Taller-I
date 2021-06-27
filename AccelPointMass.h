/** @file AccelPointMass.h
 *  @brief Specification of AccelPointMass
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef ACCELPOINTMASS_H
#define ACCELPOINTMASS_H

/** Computes the perturbational acceleration due to a point mass
 *  @param [in] r           Satellite position vector
 *  @param [in] s           Point mass position vector
 *  @param [in] GM          Gravitational coefficient of point mass
 *  @return a    		Acceleration(a = d ^ 2r / dt ^ 2)
 */
double *AccelPointMass(double *r, int nr, double *s, int ns, double GM);

#endif
