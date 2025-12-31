#ifndef DOUBLE_SCALAR_MULTIPLICATION
#define DOUBLE_SCALAR_MULTIPLICATION
#include "EC_add_proj.h"
#include "precompute_table.h"

void ec_double_scalar_multiplication(ECPointProj *R, const ECPointProj *P, 
        const ECPointProj *Q, const mpz_t k, const mpz_t l, unsigned int w, const ECCurve *E);

void ec_double_scalar_multiplication_wrapper(ECPointProj *R, const ECPointProj *P, 
        const ECPointProj *Q, const mpz_t k, const mpz_t l, const ECCurve *E);
#endif