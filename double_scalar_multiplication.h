#ifndef DOUBLE_SCALAR_MULTIPLICATION_H
#define DOUBLE_SCALAR_MULTIPLICATION_H

#include "EC_add_proj.h"
#include "precompute_table.h"

/*
 * =========================
 * Double multiplication scalaire
 * =========================
 *
 * R = k*P + l*Q en utilisant la technique de fenêtres simultanées
 * w = taille de la fenêtre
 */
void ec_double_scalar_multiplication(ECPointProj *R,
                                     const ECPointProj *P,
                                     const ECPointProj *Q,
                                     const mpz_t k,
                                     const mpz_t l,
                                     unsigned int w,
                                     const ECCurve *E);

#endif /* DOUBLE_SCALAR_MULTIPLICATION_H */