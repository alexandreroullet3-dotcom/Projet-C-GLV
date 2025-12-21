#include "EC_square_and_multiply.h"

void ec_scalar_mul(ECPoint *R, const ECPoint *P, const mpz_t k, const ECCurve *E)
{
    ECPoint Q;
    ec_point_init(&Q);
    Q.infinity = 1; // Q = point à l'infini
    ECPoint tmp;
    ec_point_init(&tmp);

    mpz_t n;
    mpz_init_set(n, k);

    ECPoint P_copy;
    ec_point_init(&P_copy);
    mpz_set(P_copy.x, P->x);
    mpz_set(P_copy.y, P->y);
    P_copy.infinity = P->infinity;

    while (mpz_cmp_ui(n, 0) > 0) {
        if (mpz_odd_p(n)) {
            ec_point_add(&tmp, &Q, &P_copy, E);
            ec_point_init(&Q);
            mpz_set(Q.x, tmp.x);
            mpz_set(Q.y, tmp.y);
            Q.infinity = tmp.infinity;
        }
        ec_point_double(&tmp, &P_copy, E);
        mpz_set(P_copy.x, tmp.x);
        mpz_set(P_copy.y, tmp.y);
        P_copy.infinity = tmp.infinity;

        mpz_fdiv_q_2exp(n, n, 1); // n = n / 2
    }

    mpz_set(R->x, Q.x);
    mpz_set(R->y, Q.y);
    R->infinity = Q.infinity;

    ec_point_clear(&Q);
    ec_point_clear(&tmp);
    ec_point_clear(&P_copy);
    mpz_clear(n);
}