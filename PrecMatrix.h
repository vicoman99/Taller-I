/** @file PrecMatrix.h
 *  @brief Specification of PrecMatrix
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef PRECMATRIX_H
#define PRECMATRIX_H

/** Precession transformation of equatorial coordinates
 *  @param [in] Mjd_1     Epoch given (Modified Julian Date TT)
 *  @param [in] MjD_2     Epoch to precess to (Modified Julian Date TT)
 *  @param [out] PrecMat   Precession transformation matrix
 */
double **PrecMatrix(double Mjd_1, double Mjd_2);

#endif
