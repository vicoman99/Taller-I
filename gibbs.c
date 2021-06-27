/** @file gibbs.h
 *  @brief Implementation of gibbs
 *  @author Victor Coman
 *  @date May 2021
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "arrays.h"
#include "gibbs.h"
#include "sat_const.h"
#include "unit.h"
#include "angl.h"

void gibbs(double *r1, double *r2, double *r3, double *v2, double *theta, double *theta1, double *copa, char *error)
{
    double small, magr1, magr2, magr3, magd, magn, r1mr2, r3mr1, r2mr3, l, tover2;
    double *p, *q, *w, *pn, *r1n, *d, *n, *nn, *dn, *s, *b, **result;
    small = 0.00000001;
    *theta = 0.0;
    *theta1 = 0.0;

    magr1 = norma(r1, 3);
    magr2 = norma(r2, 3);
    magr3 = norma(r3, 3);

    for (int i = 0; i < 3; i++)
    {
        v2[i] = 0.0;
    }

    p = crossProd(r2, 3, r3, 3);
    q = crossProd(r3, 3, r1, 3);
    w = crossProd(r1, 3, r2, 3);

    pn = unit(p, 3);
    r1n = unit(r1, 3);
    *copa = asin(dot(pn, 3, r1n, 3));

    if (fabs(dot(r1n, 3, pn, 3)) > 0.017452406)
    {
        error = "not coplanar";
    }

    d = sumV(p, 3, sumV(q, 3, w, 3), 3);
    magd = norma(d, 3);

    n = sumV(vec_x_esc(p, 3, magr1), 3, sumV(vec_x_esc(q, 3, magr2), 3, vec_x_esc(w, 3, magr3), 3), 3);
    magn = norma(n, 3);
    nn = unit(n, 3);
    dn = unit(d, 3);

    // -------------------------------------------------------------
    // determine if  the orbit is possible. both d and n must be in
    // the same direction, and non-zero.
    // -------------------------------------------------------------

    if ((fabs(magd) < small) || (fabs(magn) < small) || (dot(nn, 3, dn, 3) < small))
    {
        error = "impossible";
    }
    else
    {
        *theta = angl(r1, r2);
        *theta1 = angl(r2, r3);
    }

    // ----------- perform gibbs method to find v2 -----------

    r1mr2 = magr1 - magr2;
    r3mr1 = magr3 - magr1;
    r2mr3 = magr2 - magr3;

    s = sumV(vec_x_esc(r3, 3, r1mr2), 3, sumV(vec_x_esc(r2, 3, r3mr1), 3, vec_x_esc(r1, 3, r2mr3), 3), 3);
    b = crossProd(d, 3, r2, 3);
    l = sqrt(GM_Earth / (magd * magn));
    tover2 = l / magr2;
    v2 = sumV(vec_x_esc(b, 3, tover2), 3, vec_x_esc(s, 3, l), 3);
}
