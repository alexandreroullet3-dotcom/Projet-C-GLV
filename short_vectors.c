#include "short_vectors.h"

void z2_init(Z2 *v)
{
    mpz_init(v->x);
    mpz_init(v->y);
}

void z2_clear(Z2 *v)
{
    mpz_clear(v->x);
    mpz_clear(v->y);
}

void glv_basis(Z2 *v1, Z2 *v2, const mpz_t n, const mpz_t lambda)
{
    mpz_t r0, r1, r2;
    mpz_t t0, t1, t2;
    mpz_t q, tmp, sqrt_n;

    mpz_inits(r0, r1, r2, t0, t1, t2, q, tmp, sqrt_n, NULL);

    mpz_set(r0, n);
    mpz_set(r1, lambda);
    mpz_set_ui(t0, 0);
    mpz_set_ui(t1, 1);

    mpz_sqrt(sqrt_n, n);

    mpz_t rm, tm, rm1, tm1, rm2, tm2;
    mpz_inits(rm, tm, rm1, tm1, rm2, tm2, NULL);

    while (mpz_cmp(r1, sqrt_n) > 0) {
        mpz_fdiv_q(q, r0, r1);

        mpz_mul(tmp, q, r1);
        mpz_sub(r2, r0, tmp);

        mpz_mul(tmp, q, t1);
        mpz_sub(t2, t0, tmp);

        mpz_set(r0, r1);
        mpz_set(t0, t1);
        mpz_set(r1, r2);
        mpz_set(t1, t2);
    }

    /* À ce point :
       r0 = r_m
       r1 = r_{m+1}
       r2 = r_{m+2} (virtuel si besoin)
    */

    mpz_set(rm, r0);
    mpz_set(tm, t0);
    mpz_set(rm1, r1);
    mpz_set(tm1, t1);

    /* calcul r_{m+2} */
    mpz_fdiv_q(q, r0, r1);
    mpz_mul(tmp, q, r1);
    mpz_sub(rm2, r0, tmp);

    mpz_mul(tmp, q, t1);
    mpz_sub(tm2, t0, tmp);

    /* v1 = (r_{m+1}, -t_{m+1}) */
    z2_init(v1);
    mpz_set(v1->x, rm1);
    mpz_neg(v1->y, tm1);

    /* v2 = plus court entre (r_m, -t_m) et (r_{m+2}, -t_{m+2}) */
    mpz_t n1, n2;
    mpz_inits(n1, n2, NULL);

    mpz_mul(n1, rm, rm);
    mpz_mul(tmp, tm, tm);
    mpz_add(n1, n1, tmp);

    mpz_mul(n2, rm2, rm2);
    mpz_mul(tmp, tm2, tm2);
    mpz_add(n2, n2, tmp);

    z2_init(v2);
    if (mpz_cmp(n1, n2) <= 0) {
        mpz_set(v2->x, rm);
        mpz_neg(v2->y, tm);
    } else {
        mpz_set(v2->x, rm2);
        mpz_neg(v2->y, tm2);
    }

    mpz_clears(r0, r1, r2, t0, t1, t2, q, tmp, sqrt_n,
               rm, tm, rm1, tm1, rm2, tm2, n1, n2, NULL);
}