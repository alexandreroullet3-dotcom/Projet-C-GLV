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

void init_example3_curve(GLVCurve *curve) {
    mpz_inits(curve->E.a, curve->E.b, curve->E.p, curve->n, curve->lambda, curve->beta, curve->E.a2, NULL);
    mpz_t t; mpz_init_set_ui(t, 4);

    mpz_set_si(curve->E.a, -2);
    mpz_set_si(curve->E.b, -1);
    
    // --- CORRECTION TYPE 3 : Suppression des "0x" pour p et n ---
    mpz_set_str(curve->E.p, "b67e57937fad77013fa77527b678e2d9056747a359a88a044c9c5a431c85372b", 16);
    mpz_set_str(curve->n, "b67e57937fad77013fa77527b678e2d97d0bdb060e17a8cd13ae56ddeb4e66e0", 16);

    mpz_invert(curve->E.a2, t, curve->E.p);
    mpz_mul_ui(curve->E.a2, curve->E.a2, 3);
    mpz_neg(curve->E.a2, curve->E.a2);
    mpz_mod(curve->E.a2, curve->E.a2, curve->E.p);

    ec_point_proj_init(&curve->P);
    ec_point_proj_init(&curve->phiP);
    
    // --- CORRECTION TYPE 3 : Suppression des "0x" pour les coordonnées ---
    mpz_set_str(curve->P.X, "7a4ec4f1a9c0b1e3a12f1044099f5ef280eddb93240f9c49f6963f158b5d9a1e", 16);
    mpz_set_str(curve->P.Y, "7d156c0124cc06b1ae06e598925c39f5cd9c4985aa46e93b6e59919d55531890", 16);
    mpz_set_ui(curve->P.Z, 1);
    curve->P.infinity = 0;

    // --- CORRECTION TYPE 3 : Initialiser lambda AVANT trouver_constantes_glv (et sans 0x) ---
    mpz_set_str(curve->lambda, "3cd4c7dbd539d255bfe27c6292284b9dd459490204b28d99b13a1cf4a3c4ccf5", 16);
    mpz_set_str(curve->beta, "687d2bd31b3b44c496ef57ab0219526715054aed8e7c22e6e9ed9310f5f90bb6", 16);
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

    mpz_clear(t);
}