/** @file Accel.h
 *  @brief Specification of Accel
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef ACCEL_H
#define ACCEL_H

/** Computes the acceleration of an Earth orbiting satellite due to the Earth's harmonic gravity field, the gravitational perturbations of the Sun and Moon, the solar radiation pressure and the atmospheric drag
 *  @param [in] Mjd_TT      Terrestrial Time (Modified Julian Date)
 *  @param [in] Y           Satellite state vector in the ICRF/EME2000 system
 *  @param [out] dY		    Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
 */
void Accel(double x, double *Y, double *dY);

#endif
