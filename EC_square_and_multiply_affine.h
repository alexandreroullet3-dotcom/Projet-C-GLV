#ifndef EC_SQUARE_AND_MULTIPLY_AFFINE_H
#define EC_SQUARE_AND_MULTIPLY_AFFINE_H


#include "EC_add_affine.h"

/*
 * =========================
 * Multiplication scalaire sur points affines
 * Algorithme square-and-multiply classique
 * =========================
 *
 * R = k * P, où k est un mpz_t
 */
void ec_scalar_mul_affine(ECPointAffine *R, 
                          const ECPointAffine *P, 
                          const mpz_t k, 
                          const ECCurve *E);

#endif /* EC_SQUARE_AND_MULTIPLY_AFFINE_H */
