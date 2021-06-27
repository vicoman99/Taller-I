/** @file NutMatrix.h
 *  @brief Specification of NutMatrix
 *
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef NUTMATRIX_H__
#define NUTMATRIX_H__

/** Transformation from mean to true equator and equinox
 *  @param [in] Mjd_TT    Modified Julian Date (Terrestrial Time)
 *  @return NutMat    Nutation matrix
 */
double** NutMatrix(double Mjd_TT);

#endif