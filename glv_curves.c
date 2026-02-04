#include "glv_curves.h"
#include "EC_endo_phi_GLV.h" // Assurez-vous que vos fonctions phi y sont déclarées

void calculer_phiP_optimise(GLVCurve *curve, int type) {
    ECPointAffine P_aff, phiP_aff;
    ec_point_affine_init(&P_aff);
    ec_point_affine_init(&phiP_aff);

    // 1. Passage en affine (adapter selon votre signature réelle de proj_to_affine)
    // Si proj_to_affine prend aussi 2 arguments, retirez &curve->E ici aussi.
    proj_to_affine(&P_aff, &curve->P, &curve->E);

    // 2. Choix de l'endomorphisme
    if (type == 1) {
        ec_endo_phi1_affine(&phiP_aff, &P_aff, &curve->E, curve->beta);
    } else if (type == 2) {
        ec_endo_phi2_affine(&phiP_aff, &P_aff, &curve->E, curve->beta);
    } else if (type == 3) {
        ec_endo_phi3_affine(&phiP_aff, &P_aff, &curve->E, curve->beta);
    }

    // 3. RÉPARATION : On retire le 3ème argument (la courbe)
    affine_to_proj(&curve->phiP, &phiP_aff); 

    ec_point_affine_clear(&P_aff);
    ec_point_affine_clear(&phiP_aff);
}

// Valeurs tabulées et trouvées en ligne

void init_secp256k1_curve(GLVCurve *curve) {
    mpz_inits(curve->E.a, curve->E.b, curve->E.p, curve->n, curve->lambda, curve->beta, NULL);
    mpz_set_ui(curve->E.a, 0);
    mpz_set_ui(curve->E.b, 7);
    mpz_set_str(curve->E.p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    mpz_set_str(curve->n, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16);

    ec_point_proj_init(&curve->P);
    ec_point_proj_init(&curve->phiP);
    mpz_set_str(curve->P.X, "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_set_str(curve->P.Y, "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);
    mpz_set_ui(curve->P.Z, 1);
    curve->P.infinity = 0;

    mpz_set_str(curve->lambda, "5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72", 16);
    trouver_constantes_glv(curve->beta, &curve->E, &curve->P, curve->lambda, 1);

    // RÉPARATION : Utilisation de l'endomorphisme phi1 au lieu de la mul_scalaire
    calculer_phiP_optimise(curve, 1);

    z2_init(&curve->v1); z2_init(&curve->v2);
    glv_basis(&curve->v1, &curve->v2, curve->n, curve->lambda);
}

//Valeurs calculées avec gen_const_example2.sage

void init_example2_curve(GLVCurve *curve) {
    mpz_inits(curve->E.a, curve->E.b, curve->E.p, curve->n, curve->lambda, curve->beta, NULL);
    mpz_set_str(curve->E.a, "2e818c97303c2c8e6ee49e2b6cacc754eff27e8346e4a23da9b7b882685b7c72", 16);
    mpz_set_ui(curve->E.b, 0);
    mpz_set_str(curve->E.p, "b4fc23d83418e4d099141c1a435cbb663817e03477f8f84f3afd51e63e89ef31", 16);
    mpz_set_str(curve->n, "b4fc23d83418e4d099141c1a435cbb67505e70fbe7dfa084bca52c64f59eca7a", 16);

    ec_point_proj_init(&curve->P);
    ec_point_proj_init(&curve->phiP);
    mpz_set_str(curve->P.X, "a3cf5403796bcb5f78979ac1eac261d7d7c5785a712d01a52519542c2e6f0dc6", 16);
    mpz_set_str(curve->P.Y, "5a4073851f2f574df70c8c84906546b6ca76421fcc7db0901ffc77fa9a2ed540", 16);
    mpz_set_ui(curve->P.Z, 1);
    curve->P.infinity = 0;

    // --- CORRECTION TYPE 2 : Initialiser lambda AVANT trouver_constantes_glv ---
    mpz_set_str(curve->lambda, "81e16b4d3131f1322cf0ab2ba439286b3962df578bf08b5dc04a37e1c8e31333", 16);
    trouver_constantes_glv(curve->beta, &curve->E, &curve->P, curve->lambda, 2);

    calculer_phiP_optimise(curve, 2);

    z2_init(&curve->v1); z2_init(&curve->v2);
    glv_basis(&curve->v1, &curve->v2, curve->n, curve->lambda);
}

//Valeurs calculées avec gen_const_example3.sage

void init_example3_curve(GLVCurve *curve) {
    mpz_inits(curve->E.a, curve->E.b, curve->E.p, curve->n, curve->lambda, curve->beta, curve->E.a2, NULL);
    mpz_t t; mpz_init_set_ui(t, 4);

    mpz_set_si(curve->E.a, -2);
    mpz_set_si(curve->E.b, -1);
    
    // --- CORRECTION TYPE 3 : Suppression des "0x" pour p et n ---
    // ici on a n = h*r avec h=4, mais on ne fait GLV que sur le sous-groupe d'ordre r
    mpz_set_str(curve->E.p, "ba41f30a2231820b3fa2957bae0b9fed1801dff14e9cded971efdb63aeb369d5", 16);
    mpz_set_str(curve->n, "2e907cc2888c6082cfe8a55eeb82e7faf521130fcfa0ca47fa6ecebf5bc4759b", 16);

    mpz_invert(curve->E.a2, t, curve->E.p);
    mpz_mul_ui(curve->E.a2, curve->E.a2, 3);
    mpz_neg(curve->E.a2, curve->E.a2);
    mpz_mod(curve->E.a2, curve->E.a2, curve->E.p);

    ec_point_proj_init(&curve->P);
    ec_point_proj_init(&curve->phiP);
    
    // --- CORRECTION TYPE 3 : Suppression des "0x" pour les coordonnées ---
    mpz_set_str(curve->P.X, "b926c0c1e0bb366767cef2ecb7b6b0363611073bb6caeec225040dbd6986d21c", 16);
    mpz_set_str(curve->P.Y, "583d94cf7256ce2c92dbd84aac314a0c699a4e3b94b75baf83d67cd5dd54d11d", 16);
    mpz_set_ui(curve->P.Z, 1);
    curve->P.infinity = 0;
    mpz_clear(t);

    // --- CORRECTION TYPE 3 : Initialiser lambda AVANT trouver_constantes_glv (et sans 0x) ---
    mpz_set_str(curve->lambda, "1ea1f990009d0fe14b48af8f43d1a1bc32a0e228508b71c353c0be059898928", 16);
    mpz_set_str(curve->beta, "660bd4979c65fd51364f5530f1ef8e4ceaeaccdb85033f184764b3a1551c8c5c", 16);
    //trouver_constantes_glv(curve->beta, &curve->E, &curve->P, curve->lambda, 3);

    calculer_phiP_optimise(curve, 3);
    
// TEST DE VALIDATION (à retirer après)
    ECPointProj test_phi;
    ec_point_proj_init(&test_phi);
    ec_scalar_mul_proj(&test_phi, &curve->P, curve->lambda, &curve->E);
    if (ec_cmp_proj(&test_phi, &curve->phiP, &curve->E) != 0) {
        printf("ALERTE : phi(P) != lambda*P ! Le probleme vient de beta ou de l'arithmetique.\n");
    }
    ec_point_proj_clear(&test_phi);
///
    z2_init(&curve->v1); z2_init(&curve->v2);
    glv_basis(&curve->v1, &curve->v2, curve->n, curve->lambda);

}