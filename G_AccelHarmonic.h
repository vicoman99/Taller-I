/** @file G_AccelHarmonic.h
 *  @brief Specification of G_AccelHarmonic
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef G_ACCELHARMONIC_H
#define G_ACCELHARMONIC_H

/** Computes the gradient of the Earth's harmonic gravity field
 *  @param [in] r           Satellite position vector in the true-of-date system
 *  @param [in] U           Transformation matrix to body-fixed system
 *  @param [in] n           Gravity model degree
 *  @param [in] m 			Gravity model order
 *  @return  G    		Gradient (G=da/dr) in the true-of-date system
 */
double **G_AccelHarmonic(double *r, double **U, int nU, int n_max, int m_max);

#endif
