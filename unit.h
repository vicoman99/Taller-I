/** @file unit.h
 *  @brief Specification of unit
 *  @author Victor Coman
 *  @date May 2021
 */

#ifndef UNIT_H
#define UNIT_H

/** this function calculates a unit vector given the original vector.if a // zero vector is input, the vector is set to zero
 *  @param [in] vec - vector
 *  @return outvec - unit vector
 */
double *unit(double *vec, int n);

#endif