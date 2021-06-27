/** @file TimeUpdate.h
 *  @brief Specification of TimeUpdate
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef TIMEUPDATE_H
#define TIMEUPDATE_H

/** Update of time
 *  @param [in] P
 *  @param [in] Phi
 *  @param [in] Qdt
 *  @return array
 */
double **TimeUpdate(double **P, int np, double **Phi, int nphi, int Qdt);

#endif
