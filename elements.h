/** @file elements.h
 *  @brief Specification of elements
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef ELEMENTS_H
#define ELEMENTS_H

/** Computes the osculating Keplerian elements from the satellite state vector for elliptic orbits
 *  @brief Notes: The function cannot be used with state vectors describing a circular or non-inclined orbit.
 *  @param [in] y        State vector (x,y,z,vx,vy,vz)
 *  @param [out] p        semilatus rectum [m]
 *  @param [out] a        Semimajor axis
 *  @param [out] e        Eccentricity
 *  @param [out] i        Inclination [rad]
 *  @param [out] Omega    Longitude of the ascending node [rad]
 *  @param [out] omega    Argument of pericenter [rad]
 *  @param [out] M        Mean anomaly [rad]
 */
void elements(double *y, double *p, double *a, double *e, double *i, double *Omega, double *omega , double *M);

#endif