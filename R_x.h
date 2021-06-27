/** @file R_x.h
 *  @brief Specification of R_x
 *
 *  @author Victor Coman
 *  @date May 2021
 */

#ifndef R_X_H
#define R_X_H

/** Matrix rotation on the x axis
 *  @param [in] angle angle of rotation [rad]
 *  @return rotmat - vector result
 */
double **R_x(double angle);

#endif /* R_X_H */