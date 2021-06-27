/** @file angl.h
 *  @brief Specification of angl
 *
 *  @author Victor Coman
 *  @date May 2021
 */

#ifndef ANGL_H
#define ANGL_H

/** Greenwich Mean Sidereal Time
 *  @param [in] vec1         - vector 1
 *  @param [in] vec2         - vector 2
 *  @return theta        - angle between the two vectors  -pi to pi
 */
double angl(double *vec1, double *vec2);

#endif