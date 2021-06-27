/** @file gibbs.h
 *  @brief Specification of gibbs
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef GIBBS_H
#define GIBBS_H

/** this function performs the gibbs method of orbit determination. this method determines the velocity at the middle point of the 3 given position vectors.
 *  @param [in] r1          - ijk position vector #1         m
 *  @param [in] r2          - ijk position vector #2         m
 *  @param [in] r3          - ijk position vector #3         m
 *  @param [out] v2          - ijk velocity vector for r2     m/s
 *  @param [out] theta       - angl between vectors           rad
 */
void gibbs(double *r1, double *r2, double *r3, double *v2, double *theta, double *theta1, double *copa, char *error);

#endif