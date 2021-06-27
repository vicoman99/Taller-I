/** @file PoleMatrix.h
 *  @brief Specification of PoleMatrix
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef POLEMATRIX_H
#define POLEMATRIX_H

/** Transformation from pseudo Earth-fixed to Earth-fixed coordinates for a given date
 *  @param [in] Pole coordinte(xp,yp)
 *  @param [out] PoleMat   Pole matrix
 */
double **PoleMatrix(double xp, double yp);

#endif/* POLEMATRIX_H */
