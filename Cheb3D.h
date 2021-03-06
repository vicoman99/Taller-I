/** @file Cheb3D.h
 *  @brief Specification of Cheb3D
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef CHEB3D_H
#define CHEB3D_H

/** Chebyshev approximation of 3-dimensional vectors
 *  @param [in] N       Number of coefficients
 *  @param [in] Ta      Begin interval
 *  @param [in] Tb      End interval
 *  @param [in] Cx      Coefficients of Chebyshev polyomial (x-coordinate)
 *  @param [in] Cy      Coefficients of Chebyshev polyomial (y-coordinate)
 *  @param [in] Cz      Coefficients of Chebyshev polyomial (z-coordinate)
 *  @return vector
 */
double *Cheb3D(double t, double N, double Ta, double Tb, double *Cx, double *Cy, double *Cz);

#endif
