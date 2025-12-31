#ifndef DOUBLE_SCALAR_MULTIPLICATION
#define DOUBLE_SCALAR_MULTIPLICATION

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
void ec_double_scalar_multiplication(ECPointProj *R, const ECPointProj *P, const ECPointProj *Q, 
                                    const mpz_t k, const mpz_t l, unsigned int w, const ECCurve *E);

/*
 * Wrapper pour double multiplication scalaire
 * Choisit une taille de fenêtre par défaut (ex: w=2)
 */
void ec_double_scalar_multiplication_wrapper(ECPointProj *R, const ECPointProj *P, const ECPointProj *Q, 
                                            const mpz_t k, const mpz_t l, const ECCurve *E);

#endif /* DOUBLE_SCALAR_MULTIPLICATION */
