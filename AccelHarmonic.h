/** @file AccelHarmonic.h
 *  @brief Specification of AccelHarmonic
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef ACCELHARMONIC_H
#define ACCELHARMONIC_H

/** Computes the acceleration due to the harmonic gravity field of the central body
 *  @param [in] r           Satellite position vector in the inertial system
 *  @param [in] E           Transformation matrix to body-fixed system
 *  @param [in] n_max       Maximum degree
 *  @param [in] m_max       Maximum order (m_max<=n_max; m_max=0 for zonals, only)
 *  @return a           Acceleration (a=d^2r/dt^2)
 */
double *AccelHarmonic(double *r, double **E, int nE, int n_max, int m_max);

#endif
