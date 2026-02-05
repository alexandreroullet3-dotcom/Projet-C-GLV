#ifndef GLV_CURVES_H
#define GLV_CURVES_H

#include "glv_decompose.h"
#include "quadratic_solver.h"

/*
 * Structure pour une courbe GLV.
 */
typedef struct {
    ECCurve E;
    ECPointProj P;
    ECPointProj phiP;
    mpz_t n;
    mpz_t lambda;
    mpz_t beta;
    Z2 v1;
    Z2 v2;
} GLVCurve;

/*
 * Fonctions d'initialisation/clear.
 */
void init_secp256k1_curve(GLVCurve *curve);
void init_example2_curve(GLVCurve *curve);
void init_example3_curve(GLVCurve *curve);
void clear_curve(GLVCurve *curve);

#endif
