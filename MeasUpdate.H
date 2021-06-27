/** @file MeasUpdate.h
 *  @brief Specification of MeasUpdate
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef MEASUPDATE_H
#define MEASUPDATE_H

/** update of the measurement
 *  @param [in] x
 *  @param [in] z
 *  @param [in] g
 *  @param [in] s
 *  @param [in] G
 *  @param [in] P
 *  @param [in] N
 *  @param [out] K
 *  @param [out] x
 *  @param [in, out] P
 */
void MeasUpdate(double *x, int nx, double z, double g, double s, double *G, int ng, double **P, int np, int n, double *K);

#endif