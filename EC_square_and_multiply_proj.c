#include "EC_square_and_multiply_proj.h"

/*
 * =========================
 * Multiplication scalaire en coordonnées projectives
 * Algorithme square-and-multiply
 * =========================
 *
 * R = k * P, k est un mpz_t
 */
void ec_scalar_mul_proj(ECPointProj *R,
                        const ECPointProj *P,
                        const mpz_t k,
                        const ECCurve *E)
{
    ECPointProj Q; // Résultat temporaire initialisé au point à l'infini.
    ec_point_proj_init(&Q);
    Q.infinity = 1;

    ECPointProj Pcopy; // Copie de P pour ne pas modifier l'original.
    ec_point_proj_init(&Pcopy);
    ec_point_proj_copy(&Pcopy, P);

    ECPointProj tmp; // Variable temporaire pour doublement/ajout.
    ec_point_proj_init(&tmp);

    size_t nbits = mpz_sizeinbase(k, 2);
    if (nbits == 0) {     // k = 0 -> R = O
        ec_point_proj_copy(R, &Q);
        goto cleanup;
    }

    // Boucle de multiplication MSB → LSB
    for (size_t i = 0; i < nbits; i++) {
        size_t bit_index = nbits - 1 - i;

        // Doublement du résultat courant
        ec_point_double_proj(&tmp, &Q, E);
        ec_point_proj_copy(&Q, &tmp);

        // Ajout de P si le bit courant est 1
        if (mpz_tstbit(k, bit_index)) {
            ec_point_add_proj(&tmp, &Q, &Pcopy, E);
            ec_point_proj_copy(&Q, &tmp);
        }
    }

    // Copie finale dans R
    ec_point_proj_copy(R, &Q);

cleanup:
    ec_point_proj_clear(&Q);
    ec_point_proj_clear(&Pcopy);
    ec_point_proj_clear(&tmp);
}