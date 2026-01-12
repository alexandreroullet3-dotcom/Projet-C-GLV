#include <stdio.h>
#include <stdlib.h>
#include "EC_struct.h"
#include "double_scalar_multiplication.h"
#include "EC_square_and_multiply_proj.h" // Les additions sont deja dans square and multiply
#include "EC_square_and_multiply_affine.h"
#include "EC_conversions.h"
#include "EC_GLV.h"


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

    // 1. Initialisation des variables

    mpz_t p, n, lambda, beta, k_big, x1, y1, x2, y2;
    mpz_inits(p, n, lambda, beta, k_big, x1, y1, x2, y2, NULL);

    mpz_t k_test; // Variable pour les petits tests (2 et 5)
    mpz_init(k_test);

    // 2. Initialisation des Points et Courbe

    ECCurve E;
    mpz_inits(E.a, E.b, E.p, NULL);

    ECPointProj P, R_temp, R_glv;
    ec_point_proj_init(&P);
    ec_point_proj_init(&R_temp);
    ec_point_proj_init(&R_glv);

    ECPointAffine Affine_Disp; // Pour afficher (x,y)
    ec_point_affine_init(&Affine_Disp);

    // 3. PARAMÈTRES courbe SECP256K1 & GLV

    mpz_set_str(p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    mpz_set(E.p, p);
    mpz_set_ui(E.a, 0);
    mpz_set_ui(E.b, 7);
    mpz_set_str(n, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16);

    // Point P
    mpz_set_str(P.X, "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_set_str(P.Y, "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);
    mpz_set_ui(P.Z, 1);
    P.infinity = 0;

    // Paramètres GLV
    mpz_set_str(beta, "7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee", 16);
    mpz_set_str(lambda, "5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72", 16);
    mpz_set_str(x1, "3086d221a7d46bcde86c90e49284eb15", 16);
    mpz_set_str(y1, "e4437ed6010e88286f547fa90abfe4c3", 16);
    mpz_neg(y1, y1); 
    mpz_set_str(x2, "114ca50f7a8e2f3f657c1108d9d44cfd8", 16); 
    mpz_set_str(y2, "3086d221a7d46bcde86c90e49284eb15", 16);


    // TESTS DE BASE (P, 2P, 5P)

    printf("\n=== verification valeur simple ===\n");

    // AFFICHER P 
    proj_to_affine(&Affine_Disp, &P, &E);
    printf("1. Point de départ P :\n");
    gmp_printf("   X: %Zx\n   Y: %Zx\n", Affine_Disp.x, Affine_Disp.y);

    printf("\n----------------------------------------------------------\n");
    printf("Test 1 : Calcul de 2P (doit être pareil partout)\n");
    
    // 2P VIA DOUBLEMENT (ec_point_double_proj)
    
    ec_point_double_proj(&R_temp, &P, &E);
    proj_to_affine(&Affine_Disp, &R_temp, &E);
    printf("\n>> Méthode Doublement :\n");
    gmp_printf("   X: %Zx\n   Y: %Zx\n", Affine_Disp.x, Affine_Disp.y);

    // 2P VIA SCALAIRE CLASSIQUE (ec_scalar_mul_proj avec k=2)
    mpz_set_ui(k_test, 2);
    ec_scalar_mul_proj(&R_temp, &P, k_test, &E); 
    proj_to_affine(&Affine_Disp, &R_temp, &E);
    printf("\n>> Méthode Somme Classique (k=2) :\n");
    gmp_printf("   X: %Zx\n   Y: %Zx\n", Affine_Disp.x, Affine_Disp.y);

    // 2P VIA GLV 
    ec_scal_mul_glv(&R_glv, &P, k_test, &E, x1, y1, x2, y2, beta);
    proj_to_affine(&Affine_Disp, &R_glv, &E);
    printf("\n>> Méthode GLV (k=2) :\n");
    gmp_printf("   X: %Zx\n   Y: %Zx\n", Affine_Disp.x, Affine_Disp.y);


    printf("\n----------------------------------------------------------\n");
    printf("Test B : Calcul de 5P\n");

    // 5P VIA SCALAIRE CLASSIQUE 
    mpz_set_ui(k_test, 5);
    ec_scalar_mul_proj(&R_temp, &P, k_test, &E); 
    proj_to_affine(&Affine_Disp, &R_temp, &E);
    printf("\n>> Méthode Somme Classique (k=5) :\n");
    gmp_printf("   X: %Zx\n   Y: %Zx\n", Affine_Disp.x, Affine_Disp.y);

    // 5P VIA GLV 
    ec_scal_mul_glv(&R_glv, &P, k_test, &E, x1, y1, x2, y2, beta);
    
    if (R_glv.infinity) {
        printf("\n>> Méthode GLV (k=5) : POINT A L'INFINI\n");
    } else {
        proj_to_affine(&Affine_Disp, &R_glv, &E);
        printf("\n>> Méthode GLV (k=5) :\n");
        gmp_printf("   X: %Zx\n   Y: %Zx\n", Affine_Disp.x, Affine_Disp.y);
    }
    
    // Comparaison 
    ECPointAffine TempAff;
    ec_point_affine_init(&TempAff);
    proj_to_affine(&TempAff, &R_temp, &E);
    
    if (mpz_cmp(TempAff.x, Affine_Disp.x) == 0 && mpz_cmp(TempAff.y, Affine_Disp.y) == 0) {
         printf("\n SUCCÈS : GLV et Classique donnent le même résultat pour 5P.\n");
    } else {
         printf("\n ÉCHEC : Les résultats pour 5P sont différents !\n");
    }
    ec_point_affine_clear(&TempAff);

    // TEST élevé (Grand k)

    printf("\n\n=== TEST GLV GRAND ENTIER ===\n");
    
    mpz_set_str(k_big, "123456789ABCDEF123456789ABCDEF123456789ABCDEF123456789ABCDEF1234", 16);
    gmp_printf("k = %Zx\n", k_big);

    ec_scal_mul_glv(&R_glv, &P, k_big, &E, x1, y1, x2, y2, beta);

    printf("\nResultat Projectif (X, Y, Z) :\n");
    gmp_printf("X: %Zx\n", R_glv.X);
    gmp_printf("Y: %Zx\n", R_glv.Y);
    gmp_printf("Z: %Zx\n", R_glv.Z);

    // VÉRIFICATION COURBE
    if (R_glv.infinity) {
        printf("\n Le résultat est le point à l'infini (0).\n");
    } else {
        proj_to_affine(&Affine_Disp, &R_glv, &E);

        mpz_t lhs, rhs;
        mpz_inits(lhs, rhs, NULL);

        // y^2
        mpz_mul(lhs, Affine_Disp.y, Affine_Disp.y);
        mpz_mod(lhs, lhs, p);

        // x^3 + 7
        mpz_mul(rhs, Affine_Disp.x, Affine_Disp.x);
        mpz_mod(rhs, rhs, p);
        mpz_mul(rhs, rhs, Affine_Disp.x);
        mpz_mod(rhs, rhs, p);
        mpz_add(rhs, rhs, E.b);
        mpz_mod(rhs, rhs, p);

        if (mpz_cmp(lhs, rhs) == 0) {
            printf("\n[SUCCESS] Le point final est bien sur la courbe.\n");
        } else {
            printf("\n[ERROR] Le point n'est pas sur la courbe.\n");
        }
        mpz_clears(lhs, rhs, NULL);
    }

    // Nettoyage final

    mpz_clears(p, n, lambda, beta, k_big, k_test, x1, y1, x2, y2, E.a, E.b, E.p, NULL);
    ec_point_proj_clear(&P);
    ec_point_proj_clear(&R_glv);
    ec_point_proj_clear(&R_temp);
    ec_point_affine_clear(&Affine_Disp);
    return 0;

}
