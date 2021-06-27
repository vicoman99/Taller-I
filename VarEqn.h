/** @file VarEqn.h
 *  @brief Specification of VarEqn
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef VAREQN_H
#define VAREQN_H

/** Computes the variational equations, i.e. the derivative of the state vector and the state transition matrix
 *  @param [in] x           Time since epoch in [s]
 *  @param [in] yPhi        (6+36)-dim vector comprising the state vector (y) and the state transition matrix (Phi) in column wise storage order
 *  @param [out] yPhip       Derivative of yPhi
 */
void VarEqn(double x, double *yPhi, double *yPhip);

#endif
