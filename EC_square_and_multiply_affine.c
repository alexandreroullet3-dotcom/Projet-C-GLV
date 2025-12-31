#include "EC_square_and_multiply_affine.h"

/*
 * R = kP en coordonnées affines, MSB → LSB
 */
void ec_scalar_mul_affine(ECPointAffine *R,
                          const ECPointAffine *P,
                          const mpz_t k,
                          const ECCurve *E)
{
    ECPointAffine Q, Pcopy, tmp;
    ec_point_affine_init(&Q);
    ec_point_affine_init(&Pcopy);
    ec_point_affine_init(&tmp);

    Q.infinity = 1;                 // Q = O
    ec_point_affine_copy(&Pcopy, P);

    size_t nbits = mpz_sizeinbase(k, 2);
    if (nbits == 0) { // k = 0
        ec_point_affine_copy(R, &Q);
        goto cleanup;
    }

    // Boucle MSB → LSB
    for (size_t i = 0; i < nbits; i++) {
        size_t bit_index = nbits - 1 - i;

        // Doublement
        ec_point_double_affine(&tmp, &Q, E);
        ec_point_affine_copy(&Q, &tmp);

        // Ajout si le bit courant est 1
        if (mpz_tstbit(k, bit_index)) {
            ec_point_add_affine(&tmp, &Q, &Pcopy, E);
            ec_point_affine_copy(&Q, &tmp);
        }
    }

    // Copie finale
    ec_point_affine_copy(R, &Q);

cleanup:
    ec_point_affine_clear(&Q);
    ec_point_affine_clear(&Pcopy);
    ec_point_affine_clear(&tmp);
}