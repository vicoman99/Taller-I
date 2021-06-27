/** @file gast.h
 *  @brief Implementation of gast
 *  @author Victor Coman
 *  @date May 2021
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arrays.h"
#include "sat_const.h"
#include "AccelPointMass.h"

double *AccelPointMass(double *r, int nr, double *s, int ns, double GM)
{
	if (nr != ns)
	{
		printf("Error");
		exit(EXIT_FAILURE);
	}
	// Relative position vector of satellite w.r.t.point mass
	double *d = sumV(r, nr, vec_x_esc(s, ns, -1.0), ns);

	// Acceleration
	double *parte1 = vec_x_esc(d, nr, 1 / pow(norma(d, 3), 3.0));
	double *parte2 = vec_x_esc(s, ns, 1 / pow(norma(s, 3), 3.0));
	double *a = vec_x_esc(sumV(parte1, nr, parte2, nr), nr, -GM);

	return a;
}
