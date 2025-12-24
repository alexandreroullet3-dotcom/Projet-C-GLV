#include "EC_square_and_multiply_affine.h"

void ec_scalar_mul_affine(ECPointAffine *R, const ECPointAffine *P, const mpz_t k, const ECCurve *E)
{
    ECPointAffine Q;
    ec_point_affine_init(&Q);
    Q.infinity = 1; // Q = point à l'infini
    ECPointAffine tmp;
    ec_point_affine_init(&tmp);

    mpz_t n;
    mpz_init_set(n, k);

    ECPointAffine P_copy;
    ec_point_affine_init(&P_copy);
    mpz_set(P_copy.x, P->x);
    mpz_set(P_copy.y, P->y);
    P_copy.infinity = P->infinity;

    while (mpz_cmp_ui(n, 0) > 0) {
        if (mpz_odd_p(n)) {
            ec_point_add_affine(&tmp, &Q, &P_copy, E);
            ec_point_affine_init(&Q);
            mpz_set(Q.x, tmp.x);
            mpz_set(Q.y, tmp.y);
            Q.infinity = tmp.infinity;
        }
        ec_point_double_affine(&tmp, &P_copy, E);
        mpz_set(P_copy.x, tmp.x);
        mpz_set(P_copy.y, tmp.y);
        P_copy.infinity = tmp.infinity;

        mpz_fdiv_q_2exp(n, n, 1); // n = n / 2
    }

    mpz_set(R->x, Q.x);
    mpz_set(R->y, Q.y);
    R->infinity = Q.infinity;

    ec_point_affine_clear(&Q);
    ec_point_affine_clear(&tmp);
    ec_point_affine_clear(&P_copy);
    mpz_clear(n);
}