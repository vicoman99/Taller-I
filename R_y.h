/** @file R_y.h
 *  @brief Specification of R_y
 *
 *  @author Victor Coman
 *  @date May 2021
 */

#ifndef R_Y_H
#define R_Y_H

/** Matrix rotation on the y axis
 *  @param [in] angle angle of rotation [rad]
 *  @return rotmat - vector result
 */
double **R_y(double angle);

#endif /* R_Y_H */