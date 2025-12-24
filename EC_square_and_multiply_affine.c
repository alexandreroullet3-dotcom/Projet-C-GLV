#include "EC_square_and_multiply_affine.h"

void ec_scalar_mul_affine(ECPointAffine *R,
                          const ECPointAffine *P,
                          const mpz_t k,
                          const ECCurve *E)
{
    ECPointAffine Q, Pcopy;
    ec_point_affine_init(&Q);
    ec_point_affine_init(&Pcopy);

    Q.infinity = 1;               // Q = O
    ec_point_affine_copy(&Pcopy, P);

    mpz_t n;
    mpz_init_set(n, k);

    while (mpz_cmp_ui(n, 0) > 0) {
        if (mpz_odd_p(n)) {
            ec_point_add_affine(&Q, &Q, &Pcopy, E);
        }
        ec_point_double_affine(&Pcopy, &Pcopy, E);
        mpz_fdiv_q_2exp(n, n, 1);
    }

    ec_point_affine_copy(R, &Q);

    ec_point_affine_clear(&Q);
    ec_point_affine_clear(&Pcopy);
    mpz_clear(n);
}
