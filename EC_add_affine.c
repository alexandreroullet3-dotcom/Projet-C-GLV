#include "EC_add_affine.h"


void ec_point_double_affine(ECPointAffine *R,
                     const ECPointAffine *P,
                     const ECCurve *E)
{
    if (P->infinity) { R->infinity = 1; return; }
    // Cas y == 0 -> Infini
    if (mpz_cmp_ui(P->y, 0) == 0) { R->infinity = 1; return; }

    mpz_t lambda, num, den, x3, y3;
    mpz_inits(lambda, num, den, x3, y3, NULL);

    // lambda = (3x^2 + a) / 2y
    mpz_mul(num, P->x, P->x);
    mpz_mul_ui(num, num, 3);
    mpz_add(num, num, E->a); // 3x^2 + a

    mpz_mul_ui(den, P->y, 2);
    mpz_invert(den, den, E->p); // 1 / 2y

    mpz_mul(lambda, num, den);
    mpz_mod(lambda, lambda, E->p);

    // x3 = lambda^2 - 2x
    mpz_mul(x3, lambda, lambda);
    mpz_sub(x3, x3, P->x);
    mpz_sub(x3, x3, P->x);
    mpz_mod(x3, x3, E->p);

    // y3 = lambda(x - x3) - y
    mpz_sub(num, P->x, x3);
    mpz_mul(y3, lambda, num);
    mpz_sub(y3, y3, P->y);
    mpz_mod(y3, y3, E->p);

    // Copie finale
    mpz_set(R->x, x3);
    mpz_set(R->y, y3);
    R->infinity = 0;

    mpz_clears(lambda, num, den, x3, y3, NULL);
}


void ec_point_add_affine(ECPointAffine *R,
                  const ECPointAffine *P,
                  const ECPointAffine *Q,
                  const ECCurve *E)
{
    // 1. Gestion du point à l'infini
    if (P->infinity) { ec_point_affine_copy(R, Q); return; }
    if (Q->infinity) { ec_point_affine_copy(R, P); return; }

    // 2. Gestion des abscisses égales (x1 == x2)
    if (mpz_cmp(P->x, Q->x) == 0) {
        // Cas P == Q
        if (mpz_cmp(P->y, Q->y) == 0) {
            ec_point_double_affine(R, P, E);
            return;
        }
        
        // P == -Q -> C'est l'INFINI
        R->infinity = 1;
        return;
    }

    // 3. Variables temporaires pour le calcul
    mpz_t lambda, num, den, x3, y3;
    mpz_inits(lambda, num, den, x3, y3, NULL);

    // lambda = (y2 - y1) / (x2 - x1)
    mpz_sub(num, Q->y, P->y);   // y2 - y1
    mpz_sub(den, Q->x, P->x);   // x2 - x1
    mpz_mod(den, den, E->p);
    
    mpz_invert(den, den, E->p); // 1 / (x2 - x1)
    
    mpz_mul(lambda, num, den);
    mpz_mod(lambda, lambda, E->p); // lambda final

    // x3 = lambda^2 - x1 - x2
    mpz_mul(x3, lambda, lambda);
    mpz_sub(x3, x3, P->x);
    mpz_sub(x3, x3, Q->x);
    mpz_mod(x3, x3, E->p);

    // y3 = lambda(x1 - x3) - y1
    mpz_sub(num, P->x, x3);     // x1 - x3
    mpz_mul(y3, lambda, num);
    mpz_sub(y3, y3, P->y);
    mpz_mod(y3, y3, E->p);

    // 4. Copie finale dans R
    mpz_set(R->x, x3);
    mpz_set(R->y, y3);
    R->infinity = 0;

    // 5. Nettoyage
    mpz_clears(lambda, num, den, x3, y3, NULL);
}