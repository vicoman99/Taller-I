/** @file Legendre.h
 *  @brief Specification of Legendre
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef LEGENDRE_H
#define LEGENDRE_H

/**
 *  @param [in] n
 *  @param [in] m
 *  @param [in] fi   [rad]
 *  @param [out] pnm
 *  @param [out] dpnm
 */
void Legendre(int n, int m, double fi, double **pnm, double **dpnm);

#endif