/** @file hgibbs.h
 *  @brief Specification of hgibbs
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef HGIBBS_H
#define HGIBBS_H

/** this function implements the herrick-gibbs approximation for orbit determination, and finds the middle velocity vector for the 3 given position vectors.
 *  @param [in] r1          - ijk position vector #1         m
 *  @param [in] r2          - ijk position vector #2         m
 *  @param [in] r3          - ijk position vector #3         m
 *  @param [in] Mjd1        - julian date of 1st sighting    days from 4713 bc
 *  @param [in] Mjd2        - julian date of 2nd sighting    days from 4713 bc
 *  @param [in] Mjd3        - julian date of 3rd sighting    days from 4713 bc
 *  @param [out] v2          - ijk velocity vector for r2     m/s
 *  @param [out] theta       - angl between vectors           rad
 */
void hgibbs(double *r1, double *r2, double *r3, double Mjd1, double Mjd2, double Mjd3, double *v2, double *theta, double *theta1, double *copa, char *error);

#endif
