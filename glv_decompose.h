#ifndef GLV_DECOMPOSE_H
#define GLV_DECOMPOSE_H

#include "short_vectors.h"

void mpz_round_div(mpz_t r, const mpz_t num, const mpz_t den);

void glv_nearest_vector(Z2 *v,
                        const mpz_t k,
                        const Z2 *v1,
                        const Z2 *v2);

#endif /* GLV_DECOMPOSE_H */