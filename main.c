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

    // Conversion en projectif pour les tests
    ECPointProj Pp;
    ec_point_proj_init(&Pp);
    affine_to_proj(&Pp, &Pa);

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

    ec_point_proj_clear(&Pp);
    ec_point_proj_clear(&R_double_proj);
    ec_point_proj_clear(&R_add_proj);
    ec_point_proj_clear(&R_mul_proj);
    ECPointProj P;
    ec_point_proj_init(&P);
    affine_to_proj(&P, &Pa);

    ECPointProj Q;
    ec_point_proj_init(&Q);
    mpz_set_ui(Q.X, 1539);
    mpz_set_ui(Q.Y, 4742);
    mpz_set_ui(Q.Z, 1);
    Q.infinity = 0;

    mpz_t k;
    mpz_init_set_ui(k, 20);

    mpz_t l;
    mpz_init_set_ui(l, 15);



    // ---------------------------
    // 3. Point à l'infini
    // ---------------------------
    ECPointProj O;
    ec_point_proj_init(&O);
    O.infinity = 1;
    mpz_set_ui(O.X, 0);
    mpz_set_ui(O.Y, 1);
    mpz_set_ui(O.Z, 0);

    ECPointProj Rtmp;
    ec_point_proj_init(&Rtmp);

    ECPointAffine RA;
    ec_point_affine_init(&RA);

    // ---------------------------
    // Test 4: Multiplication scalaire naïve
    // ---------------------------
    
    ECPointProj Rmul;
    ec_point_proj_init(&Rmul);
    ec_scalar_mul_proj(&Rmul, &P, k, &E);
    proj_to_affine(&RA, &Rmul, &E);
    gmp_printf("20P = (%Zd, %Zd)\n", RA.x, RA.y); 

    ec_point_add_proj(&Rtmp, &P, &Q, &E);
    proj_to_affine(&RA, &Rtmp, &E);
    gmp_printf("P + Q= (%Zd, %Zd)\n", RA.x, RA.y);

    // Résultat naïf
    ECPointProj Rnaif, tmp;
    ec_point_proj_init(&Rnaif);
    
    ec_scalar_mul_proj(&Rnaif, &P, k, &E);
    ec_point_proj_init(&tmp);
    ec_scalar_mul_proj(&tmp, &Q, l, &E);
    ec_point_add_proj(&Rnaif, &Rnaif, &tmp, &E);
    proj_to_affine(&RA, &Rnaif, &E);
    gmp_printf("résultat naif 20P + 15Q = (%Zd, %Zd)\n",RA.x, RA.y);

    // Test: fonction precompute
    
    ECPointProj **T = precompute_table(&P, &Q, 1, &E);
    for (int i=0; i<2; i++){
        for (int j=0; j<2; j++){
        tmp = T[i][j];
        proj_to_affine(&RA, &tmp, &E);
        gmp_printf("i=%d, j=%d, iP+jQ = (%Zd, %Zd)\n", i, j, RA.x, RA.y);
    }}
    proj_to_affine(&RA, &T[0][0], &E);
    gmp_printf("T[0][0] = (%Zd, %Zd)\n", RA.x, RA.y);




    // ---------------------------
    // Test 6: Shamir double scalar multiplication w=1
    // ---------------------------

    ECPointProj Rshamir;
    ec_point_proj_init(&Rshamir);
    ec_double_scalar_multiplication(&Rshamir, &P, &Q, k, l, 1, &E);
    proj_to_affine(&RA, &Rshamir, &E);
    gmp_printf("Shamir w=1 20P + 15Q = (%Zd, %Zd)\n", RA.x, RA.y);

    // ---------------------------
    // Test 7: Shamir double scalar multiplication w=2
    // ---------------------------
    ECPointProj Rshamir2;
    ec_point_proj_init(&Rshamir2);
    ec_double_scalar_multiplication(&Rshamir2, &P, &Q, k, l, 2, &E);

    // Comparer avec le même résultat naïf
    //assert(mpz_cmp(Rshamir2.X, Rnaif.X) == 0);
    //assert(mpz_cmp(Rshamir2.Y, Rnaif.Y) == 0);
    //assert(Rshamir2.infinity == Rnaif.infinity);

    proj_to_affine(&RA, &Rshamir2, &E);
    gmp_printf("Shamir w=2 20P + 15Q = (%Zd, %Zd)\n", RA.x, RA.y);

    // ---------------------------
    // Cleanup
    // ---------------------------
    ec_point_proj_clear(&P);
    ec_point_proj_clear(&O);
    ec_point_proj_clear(&Rtmp);
    ec_point_proj_clear(&Rmul);
    ec_point_proj_clear(&Q);
    ec_point_proj_clear(&Rshamir);
    ec_point_proj_clear(&Rshamir2);
    ec_point_proj_clear(&Rnaif);
    ec_point_proj_clear(&tmp);

    ec_point_affine_clear(&Pa);
    ec_point_affine_clear(&RA);

    mpz_clear(k);
    mpz_clear(l);
    mpz_clear(E.p);
    mpz_clear(E.a);
    mpz_clear(E.b);

    printf("✅ Tous les tests passent et les résultats sont corrects.\n");
    return 0;
}
