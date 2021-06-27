/** @file R_z.h
 *  @brief Specification of R_z
 *
 *  @author Victor Coman
 *  @date May 2021
 */

#ifndef R_Z_H
#define R_Z_H

/** Matrix rotation on the y axis
 *  @param [in] angle angle of rotation [rad]
 *  @return rotmat - vector result
 */
double **R_z(double angle);

#endif /* R_Z_H */