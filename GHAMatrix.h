/** @file GHAMatrix.h
 *  @brief Specification of GHAMatrix
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef GHAMATRIX_H
#define GHAMATRIX_H

/** Transformation from true equator and equinox to Earth equator and Greenwich meridian system
 *  @param [in] Mjd_UT1   Modified Julian Date UT1
 *  @param [out] GHAmat    Greenwich Hour Angle matrix
 */
double **GHAMatrix(double Mjd_UT1);

#endif
