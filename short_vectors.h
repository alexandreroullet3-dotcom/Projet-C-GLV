#ifndef SHORT_VECTORS
#define SHORT_VECTORS

#include "EC_struct.h"

typedef struct {
    mpz_t x;
    mpz_t y;
} Z2;

void z2_init(Z2 *v);

void z2_clear(Z2 *v);

void glv_basis(Z2 *v1, Z2 *v2, const mpz_t n, const mpz_t lambda);

void mpz_round_div(mpz_t r, const mpz_t num, const mpz_t den);

#endif