#include <stdio.h>
#include "EC_struct.h"
#include "double_scalar_multiplication.h"
#include "EC_square_and_multiply_proj.h" // Les additions sont deja dans square and multiply
#include "EC_square_and_multiply_affine.h"


int main() {
    // ---------------------------
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
    ec_point_proj_clear(&tmp2);
    
    
    return 0;
}
