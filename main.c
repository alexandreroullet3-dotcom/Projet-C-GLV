#include <stdio.h>
#include <stdlib.h>
#include "EC_struct.h"
#include "double_scalar_multiplication.h"
#include "EC_square_and_multiply_proj.h" // Les additions sont deja dans square and multiply
#include "EC_square_and_multiply_affine.h"
#include "EC_conversions.h"
#include "EC_glv.h"


int main() {
    /*// ---------------------------
    // 1. Définir la courbe
    // ---------------------------
    ECCurve E;
    mpz_init_set_str(E.p, "9739", 10);
    mpz_init_set_ui(E.a, 497);
    mpz_init_set_ui(E.b, 1768);

    // ---------------------------
    // 2. Point affine de test
    // ---------------------------
    ECPointAffine Pa;
    ec_point_affine_init(&Pa);
    mpz_set_ui(Pa.x, 493);
    mpz_set_ui(Pa.y, 5564);
    Pa.infinity = 0;

    ECPointAffine Qaff;
    ec_point_affine_init(&Qaff);
    mpz_set_ui(Qaff.x, 1539);
    mpz_set_ui(Qaff.y, 4742);
    Qaff.infinity = 0;
    
    // Conversion en projectif pour les tests
    ECPointProj Pp;
    ec_point_proj_init(&Pp);
    affine_to_proj(&Pp, &Pa);
    
    ECPointProj Qp;
    ec_point_proj_init(&Qp);
    affine_to_proj(&Qp, &Qaff);

    gmp_printf("Point P départ : (%Zd, %Zd)\n", Pa.x, Pa.y);

    // =================================================================
    // TEST 1 : DOUBLEMENT (2P)
    // =================================================================

    // A. Calcul Affine
    ECPointAffine R_double_aff;
    ec_point_affine_init(&R_double_aff);
    ec_point_double_affine(&R_double_aff, &Pa, &E);
    gmp_printf("[Affine]   2P = (%Zd, %Zd)\n", R_double_aff.x, R_double_aff.y);

    // B. Calcul Projectif
    ECPointProj R_double_proj;
    ec_point_proj_init(&R_double_proj);
    ec_point_double_proj(&R_double_proj, &Pp, &E);
    
    // Conversion résultat projectif -> affine pour affichage
    ECPointAffine R_double_proj_conv;
    ec_point_affine_init(&R_double_proj_conv);
    proj_to_affine(&R_double_proj_conv, &R_double_proj, &E);
    gmp_printf("[Projectif] 2P = (%Zd, %Zd)\n", R_double_proj_conv.x, R_double_proj_conv.y);

    // =================================================================
    // TEST 2 : ADDITION (P + P)
    // =================================================================

    // A. Calcul Affine
    ECPointAffine R_add_aff;
    ec_point_affine_init(&R_add_aff);
    ec_point_add_affine(&R_add_aff, &Pa, &Pa, &E);
    gmp_printf("[Affine]   P+P = (%Zd, %Zd)\n", R_add_aff.x, R_add_aff.y);

    // B. Calcul Projectif
    ECPointProj R_add_proj;
    ec_point_proj_init(&R_add_proj);
    ec_point_add_proj(&R_add_proj, &Pp, &Pp, &E);

    ECPointAffine R_add_proj_conv;
    ec_point_affine_init(&R_add_proj_conv);
    proj_to_affine(&R_add_proj_conv, &R_add_proj, &E);
    gmp_printf("[Projectif] P+P = (%Zd, %Zd)\n", R_add_proj_conv.x, R_add_proj_conv.y);

    // =================================================================
    // TEST 3 : MULTIPLICATION SCALAIRE (20P)
    // =================================================================

    mpz_t k;
    mpz_init_set_ui(k, 20);

    // A. Calcul Affine
    ECPointAffine R_mul_aff;
    ec_point_affine_init(&R_mul_aff);
    ec_scalar_mul_affine(&R_mul_aff, &Pa, k, &E);
    gmp_printf("[Affine]   20P = (%Zd, %Zd)\n", R_mul_aff.x, R_mul_aff.y);

    // B. Calcul Projectif
    ECPointProj R_mul_proj;
    ec_point_proj_init(&R_mul_proj);
    ec_scalar_mul_proj(&R_mul_proj, &Pp, k, &E);

    ECPointAffine R_mul_proj_conv;
    ec_point_affine_init(&R_mul_proj_conv);
    proj_to_affine(&R_mul_proj_conv, &R_mul_proj, &E);
    gmp_printf("[Projectif] 20P = (%Zd, %Zd)\n", R_mul_proj_conv.x, R_mul_proj_conv.y);

    //------------------------------------------
    //Calcul de 20P+15Q avec l'algorithme 1
    //-----------------------------------------

    //calcul naif

    ECPointProj Rnaif, tmp1, tmp2;
    ec_point_proj_init(&Rnaif);
    ec_point_proj_init(&tmp1);
    ec_point_proj_init(&tmp2);

    ECPointAffine Rnaifaff;
    ec_point_affine_init(&Rnaifaff);

    mpz_t l;
    mpz_init_set_ui(l, 15);
    
    ec_scalar_mul_proj(&tmp1, &Pp, k, &E);
    ec_scalar_mul_proj(&tmp2, &Qp, l, &E);
    ec_point_add_proj(&Rnaif, &tmp1, &tmp2, &E);
    proj_to_affine(&Rnaifaff, &Rnaif, &E);
    gmp_printf("Manière naïve: 20P+ + 15Q = (%Zd, %Zd)\n", Rnaifaff.x, Rnaifaff.y);

    //Avec l'algo 1
    ECPointProj Rshamir;
    ec_point_proj_init(&Rshamir);

    ECPointAffine Rshamiraff;
    ec_point_affine_init(&Rshamiraff);

    ec_double_scalar_multiplication(&Rshamir, &Pp, &Qp, k, l, 1, &E);
    proj_to_affine(&Rshamiraff, &Rshamir, &E);
    gmp_printf("Methode shamir: 20P+ + 15Q = (%Zd, %Zd)\n", Rshamiraff.x, Rshamiraff.y);

    // --- Nettoyage ---
    ec_curve_clear(&E);
    mpz_clear(k);
    
    ec_point_affine_clear(&Pa);
    ec_point_affine_clear(&R_double_aff);
    ec_point_affine_clear(&R_double_proj_conv);
    ec_point_affine_clear(&R_add_aff);
    ec_point_affine_clear(&R_add_proj_conv);
    ec_point_affine_clear(&R_mul_aff);
    ec_point_affine_clear(&R_mul_proj_conv);
    ec_point_affine_clear(&Rnaifaff);
    ec_point_affine_clear(&Rshamiraff);

    ec_point_proj_clear(&Pp);
    ec_point_proj_clear(&R_double_proj);
    ec_point_proj_clear(&R_add_proj);
    ec_point_proj_clear(&R_mul_proj);
    ec_point_proj_clear(&Rnaif);
    ec_point_proj_clear(&Rshamir);
    ec_point_proj_clear(&tmp1);
    ec_point_proj_clear(&tmp2);*/
    
//////////////////  TEST de la fonction GLV //////////////////////////

// 1. Initialisation des variables GMP
    mpz_t p, n, lambda, beta, k, x1, y1, x2, y2;
    mpz_inits(p, n, lambda, beta, k, x1, y1, x2, y2, NULL);

    // 2. Initialisation des Points
    ECPointProj P, R_glv;
    ec_point_proj_init(&P);
    ec_point_proj_init(&R_glv);

    // 3. Initialisation de la Courbe
    ECCurve E;
    mpz_inits(E.a, E.b, E.p, NULL);

    // =================================================================
    // PARAMÈTRES DE LA COURBE secp256k1
    // =================================================================
    
    // Corps p = 2^256 - 2^32 - 977
    mpz_set_str(p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    mpz_set(E.p, p);
    
    // Equation y^2 = x^3 + 7
    mpz_set_ui(E.a, 0);
    mpz_set_ui(E.b, 7);
    
    // Ordre n
    mpz_set_str(n, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16);

    // Point Générateur P
    mpz_set_str(P.X, "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_set_str(P.Y, "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);
    mpz_set_ui(P.Z, 1);
    P.infinity = 0;

    // =================================================================
    // PARAMÈTRES GLV (Pré-calculés pour secp256k1)
    // =================================================================
    
    // Beta (Racine cubique de l'unité mod p)
    mpz_set_str(beta, "7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee", 16);
    
    // Lambda (Racine cubique de l'unité mod n)
    mpz_set_str(lambda, "5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72", 16);

    // Vecteurs de la base réduite (Attention aux signes !)
    mpz_set_str(x1, "3086d221a7d46bcde86c90e49284eb15", 16);
    
    mpz_set_str(y1, "e4437ed6010e88286f547fa90abfe4c3", 16);
    mpz_neg(y1, y1); // y1 est négatif !

    mpz_set_str(x2, "114ca50f7a8e2f3f657c1108d9d44cfd8", 16); 
    
    mpz_set_str(y2, "3086d221a7d46bcde86c90e49284eb15", 16);

    // =================================================================
    // TEST
    // =================================================================
    
    // Choix d'un scalaire k aléatoire
    // Ici un nombre arbitraire de 256 bits
    mpz_set_str(k, "123456789ABCDEF123456789ABCDEF123456789ABCDEF123456789ABCDEF1234", 16);

    printf("=== Lancement du Test GLV ===\n");
    gmp_printf("k = %Zx\n", k);

    // Je les laisse par défaut pour correspondre à ton prototype actuel.
    ec_scal_mul_glv(&R_glv, &P, k, &E, x1, y1, x2, y2, beta);

    printf("\nResultat Projectif (X, Y, Z) :\n");
    gmp_printf("X: %Zx\n", R_glv.X);
    gmp_printf("Y: %Zx\n", R_glv.Y);
    gmp_printf("Z: %Zx\n", R_glv.Z);

    // VÉRIFICATION si le point est sur la courbe
    
    // Conversion en Affine pour vérifier l'équation y^2 = x^3 + 7
    ECPointAffine Check;
    ec_point_affine_init(&Check);
    
    if (R_glv.infinity) {
        printf("\n[INFO] Le résultat est le point à l'infini (0).\n");
    } else {
        proj_to_affine(&Check, &R_glv, &E);

        mpz_t lhs, rhs;
        mpz_inits(lhs, rhs, NULL);

        // Calcule y^2
        mpz_mul(lhs, Check.y, Check.y);
        mpz_mod(lhs, lhs, p);

        // Calcule x^3 + 7
        mpz_mul(rhs, Check.x, Check.x); // x^2
        mpz_mod(rhs, rhs, p);
        mpz_mul(rhs, rhs, Check.x); // x^3
        mpz_mod(rhs, rhs, p);
        mpz_add(rhs, rhs, E.b); // + 7
        mpz_mod(rhs, rhs, p);

        // Comparaison
        if (mpz_cmp(lhs, rhs) == 0) {
            printf("\n[SUCCESS] Le point est sur la courbe.\n");
        } else {
            printf("\n[ERROR] Le point n'est pas sur la courbe.\n");
            gmp_printf("Attendu (x^3+7): %Zx\n", rhs);
            gmp_printf("Obtenu (y^2)   : %Zx\n", lhs);
        }
        mpz_clears(lhs, rhs, NULL);
    }

    // Nettoyage final
    mpz_clears(p, n, lambda, beta, k, x1, y1, x2, y2, E.a, E.b, E.p, NULL);
    ec_point_proj_clear(&P);
    ec_point_proj_clear(&R_glv);
    ec_point_affine_clear(&Check);

    return 0;

}
