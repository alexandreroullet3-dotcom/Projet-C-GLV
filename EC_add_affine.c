#include "EC_add.h"

void ec_point_add(ECPoint *R,
                  const ECPoint *P,
                  const ECPoint *Q,
                  const ECCurve *E)
{
    if (P->infinity) { *R = *Q; return; }
    if (Q->infinity) { *R = *P; return; }

    // P == -Q
    if (mpz_cmp(P->x, Q->x) == 0) {
        mpz_t tmp;
        mpz_init(tmp);
        mpz_add(tmp, P->y, Q->y);
        mpz_mod(tmp, tmp, E->p);
        if (mpz_cmp_ui(tmp, 0) == 0) {
            R->infinity = 1;
            mpz_clear(tmp);
            return;
        }
        mpz_clear(tmp);
    }

    mpz_t lambda, tmp;
    mpz_inits(lambda, tmp, NULL);

    // lambda = (y2 - y1)/(x2 - x1)
    mpz_sub(lambda, Q->y, P->y);
    mpz_sub(tmp, Q->x, P->x);
    mpz_mod(tmp, tmp, E->p);
    mpz_invert(tmp, tmp, E->p);
    mpz_mul(lambda, lambda, tmp);
    mpz_mod(lambda, lambda, E->p);

    // x3 = lambda^2 - x1 - x2
    mpz_mul(R->x, lambda, lambda);
    mpz_sub(R->x, R->x, P->x);
    mpz_sub(R->x, R->x, Q->x);
    mpz_mod(R->x, R->x, E->p);

    // y3 = lambda(x1 - x3) - y1
    mpz_sub(tmp, P->x, R->x);
    mpz_mul(tmp, lambda, tmp);
    mpz_sub(tmp, tmp, P->y);
    mpz_mod(R->y, tmp, E->p);

    R->infinity = 0;

    mpz_clears(lambda, tmp, NULL);
}


void ec_point_double(ECPoint *R,
                     const ECPoint *P,
                     const ECCurve *E)
{
    if (P->infinity) { R->infinity = 1; return; }

    mpz_t lambda, tmp;
    mpz_inits(lambda, tmp, NULL);

    // lambda = (3*x^2 + a)/(2*y)
    mpz_mul(tmp, P->x, P->x);
    mpz_mul_ui(tmp, tmp, 3);
    mpz_add(tmp, tmp, E->a);

    mpz_mul_ui(lambda, P->y, 2);
    mpz_invert(lambda, lambda, E->p);
    mpz_mul(lambda, tmp, lambda);
    mpz_mod(lambda, lambda, E->p);

    // x3
    mpz_mul(R->x, lambda, lambda);
    mpz_sub(R->x, R->x, P->x);
    mpz_sub(R->x, R->x, P->x);
    mpz_mod(R->x, R->x, E->p);

    // y3
    mpz_sub(tmp, P->x, R->x);
    mpz_mul(tmp, lambda, tmp);
    mpz_sub(tmp, tmp, P->y);
    mpz_mod(R->y, tmp, E->p);

    R->infinity = 0;

    mpz_clears(lambda, tmp, NULL);
}