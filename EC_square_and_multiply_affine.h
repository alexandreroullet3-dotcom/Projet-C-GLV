#ifndef EC_SQUARE_AND_MULTIPLY_H
#define EC_SQUARE_AND_MULTIPLY_H

#include "EC_struct.h"
#include "EC_add.h"

// multiplication scalaire : R = k * P
// k est un mpz_t
void ec_scalar_mul(ECPoint *R, const ECPoint *P, const mpz_t k, const ECCurve *E);

#endif