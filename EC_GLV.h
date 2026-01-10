#ifndef EC_GLV_H
#define EC_GLV_H
#include "EC_struct.h"

void ec_scal_mul_glv(ECPointProj *R, const ECPointProj *P, const mpz_t k, const ECCurve * E, const mpz_t x1, 
                    const mpz_t y1, const mpz_t x2, const mpz_t y2, const mpz_t beta);
#endif