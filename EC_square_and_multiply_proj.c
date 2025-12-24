#include "EC_square_and_multiply_proj.h"

void ec_scalar_mul_proj(ECPointProj *R,
                        const ECPointProj *P,
                        const mpz_t k,
                        const ECCurve *E)
{
    ECPointProj Q, Pcopy;
    ec_point_proj_init(&Q);
    ec_point_proj_init(&Pcopy);

    Q.infinity = 1;            // Q = O
    ec_point_proj_copy(&Pcopy, P);

    mpz_t n;
    mpz_init_set(n, k);

    while (mpz_cmp_ui(n, 0) > 0) {
        if (mpz_odd_p(n)) {
            ec_point_add_proj(&Q, &Q, &Pcopy, E);
        }
        ec_point_double_proj(&Pcopy, &Pcopy, E);
        mpz_fdiv_q_2exp(n, n, 1);
    }

    ec_point_proj_copy(R, &Q);

    ec_point_proj_clear(&Q);
    ec_point_proj_clear(&Pcopy);
    mpz_clear(n);
}
