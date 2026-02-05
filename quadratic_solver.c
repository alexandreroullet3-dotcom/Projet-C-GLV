#include "quadratic_solver.h"

#include <stdio.h>
#include <stdlib.h>

/*
 * Tonelli–Shanks
 * Résout r^2 ≡ n (mod p), p premier impair
 * Retourne 1 si solution, 0 sinon
 */
int mpz_tonelli_shanks(mpz_t r, const mpz_t n, const mpz_t p)
{
    /* Cas triviaux */
    if (mpz_cmp_ui(n, 0) == 0) {
        mpz_set_ui(r, 0);
        return 1;
    }

    /* Vérification résidu quadratique */
    if (mpz_legendre(n, p) != 1) {
        return 0;
    }
    
    mpz_t tmp;
    mpz_init(tmp);
    /* Cas rapide p ≡ 3 (mod 4) */
    if (mpz_mod_ui(tmp, p, 4) == 3) {
        mpz_t e;
        mpz_init(e);
        mpz_add_ui(e, p, 1);
        mpz_fdiv_q_2exp(e, e, 2);  // (p+1)/4
        mpz_powm(r, n, e, p);
        mpz_clears(e, tmp, NULL);
        return 1;
    }

    /* Décomposition p-1 = q * 2^s avec q impair */
    mpz_t q;
    mpz_init(q);
    mpz_sub_ui(q, p, 1);

    unsigned long s = 0;
    while (mpz_even_p(q)) {
        mpz_fdiv_q_2exp(q, q, 1);
        s++;
    }

    /* Trouver z non-résidu quadratique */
    mpz_t z;
    mpz_init_set_ui(z, 2);
    while (mpz_legendre(z, p) != -1) {
        mpz_add_ui(z, z, 1);
    }

    /* Initialisations */
    mpz_t c, t;
    mpz_inits(c, t, NULL);

    mpz_powm(c, z, q, p);    // c = z^q
    mpz_powm(t, n, q, p);    // t = n^q

    mpz_t e;
    mpz_init(e);
    mpz_add_ui(e, q, 1);
    mpz_fdiv_q_2exp(e, e, 1); // (q+1)/2
    mpz_powm(r, n, e, p);     // r = n^((q+1)/2)

    unsigned long m = s;

    /* Boucle principale */
    while (mpz_cmp_ui(t, 1) != 0) {
        
        mpz_set(tmp, t);

        unsigned long i;
        for (i = 1; i < m; i++) {
            mpz_mul(tmp, tmp, tmp);
            mpz_mod(tmp, tmp, p);
            if (mpz_cmp_ui(tmp, 1) == 0) {
                break;
            }
        }

        if (i == m) {
            mpz_clears(q, z, c, t, e, tmp, NULL);
            return 0;
        }

        mpz_t b;
        mpz_init(b);
        mpz_powm_ui(b, c, 1UL << (m - i - 1), p);

        mpz_mul(r, r, b);
        mpz_mod(r, r, p);

        mpz_mul(c, b, b);
        mpz_mod(c, c, p);

        mpz_mul(t, t, c);
        mpz_mod(t, t, p);

        m = i;

        mpz_clear(b);
    }

    mpz_clears(q, z, c, t, e, tmp, NULL);
    return 1;
}

// focntion qui Résout une equation x^2 +Ax +B =0 mod m et stock les solution dasn r1 et r2
void solve_quadratic_equation(mpz_t r1, mpz_t r2, mpz_t A, mpz_t B, const mpz_t m)
{
    mpz_t D, sqrt_D, inv_2, t;
    mpz_inits(D, sqrt_D, inv_2, t, NULL);

    // D = A^2 - 4B
    mpz_mul(D, A, A);
    mpz_set_ui(t, 4);
    mpz_mul(t, t, B);
    mpz_sub(D, D, t);
    mpz_mod(D, D, m);

    // Racine carrée de D (utilise la fonction mpz_tonelli_shanks)
    if (mpz_tonelli_shanks(sqrt_D, D, m) == 0) {
        printf("Erreur : Pas de solution quadratique ! Les entrées étaient\n");
        gmp_printf("A = %Zd, B = %Zd, m = %Zx\n", A, B, m);
        exit(1);
    }
    // Calcul de : (-A ± sqrt) / 2
    mpz_set_ui(inv_2, 2);
    mpz_invert(inv_2, inv_2, m);
    mpz_neg(t, A);
    
    // r1
    mpz_add(r1, t, sqrt_D);
    mpz_mul(r1, r1, inv_2);
    mpz_mod(r1, r1, m);

    // r2
    mpz_sub(r2, t, sqrt_D);
    mpz_mul(r2, r2, inv_2);
    mpz_mod(r2, r2, m);

    mpz_clears(D, sqrt_D, inv_2, t, NULL);
}

/*
 * Calcule beta en testant les racines de l'équation caractéristique.
 * La signature inclut P et lambda pour vérifier la correspondance.
 */
void trouver_constantes_glv(mpz_t beta,
                            const ECCurve *E,
                            const ECPointProj *P,
                            const mpz_t lambda,
                            int type)
{
    mpz_t A, B, b1, b2;
    mpz_inits(A, B, b1, b2, NULL);

    // 1. Définition de l'équation caractéristique selon le type.
    if (type == 1) {        // y^2 = x^3 + b : x^2 + x + 1 = 0
        mpz_set_ui(A, 1);
        mpz_set_ui(B, 1);
    } else if (type == 2) { // y^2 = x^3 + ax : x^2 + 1 = 0
        mpz_set_ui(A, 0);
        mpz_set_ui(B, 1);
    } else if (type == 3) { // y^2 = x^3 - 3/4x^2 - 2x - 1 : x^2 - x + 2 = 0
        mpz_sub_ui(A, E->p, 1); // A = -1 mod p
        mpz_set_ui(B, 2);
    }

    // 2. Calcul des deux racines possibles pour beta mod p
    solve_quadratic_equation(b1, b2, A, B, E->p);

    // 3. VÉRIFICATION : On teste laquelle des deux racines correspond à lambda
    ECPointProj R_test, phiP_test;
    ec_point_proj_init(&R_test);
    ec_point_proj_init(&phiP_test);

    // On calcule R = [lambda]P de manière classique [cite: 10]
    ec_scalar_mul_proj(&R_test, P, lambda, E);

    // On teste b1 avec l'endomorphisme
    ECPointAffine P_aff, phiP_aff;
    ec_point_affine_init(&P_aff);
    ec_point_affine_init(&phiP_aff);
    proj_to_affine(&P_aff, P, E);

    if (type == 1) {
        ec_endo_phi1_affine(&phiP_aff, &P_aff, E, b1);
    } else if (type == 2) {
        ec_endo_phi2_affine(&phiP_aff, &P_aff, E, b1);
    } else if (type == 3) {
        ec_endo_phi3_affine(&phiP_aff, &P_aff, E, b1);
    }

    affine_to_proj(&phiP_test, &phiP_aff);

    // Si phi_b1(P) == [lambda]P, alors beta = b1. Sinon, beta = b2.
    if (ec_cmp_proj(&phiP_test, &R_test, E) == 0) {
        mpz_set(beta, b1);
    } else {
        mpz_set(beta, b2);
    }

    // Nettoyage
    ec_point_proj_clear(&R_test);
    ec_point_proj_clear(&phiP_test);
    ec_point_affine_clear(&P_aff);
    ec_point_affine_clear(&phiP_aff);
    mpz_clears(A, B, b1, b2, NULL);
}