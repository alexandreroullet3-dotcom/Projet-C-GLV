#ifndef EC_SQUARE_AND_MULTIPLY_AFFINE_H
#define EC_SQUARE_AND_MULTIPLY_AFFINE_H

#include "EC_struct.h"
#include "EC_add_affine.h"

// multiplication scalaire : R = k * P
// k est un mpz_t
void ec_scalar_mul_affine(ECPointAffine *R, const ECPointAffine *P, const mpz_t k, const ECCurve *E);

#endif