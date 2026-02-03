#include "glv_curves.h"

void init_secp256k1_curve(GLVCurve *curve) {
    // Initialisation de la courbe
    mpz_inits(curve->E.a, curve->E.b, curve->E.p, curve->n, NULL);
    mpz_set_ui(curve->E.a, 0);
    mpz_set_ui(curve->E.b, 7);
    mpz_set_str(curve->E.p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    mpz_set_str(curve->n, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16);

    // Point générateur
    ec_point_proj_init(&curve->P);
    mpz_set_str(curve->P.X, "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_set_str(curve->P.Y, "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);
    mpz_set_ui(curve->P.Z, 1);
    curve->P.infinity = 0;

    // GLV paramètres
    mpz_inits(curve->lambda, curve->beta, NULL);
    trouver_constantes_glv(curve->beta, &curve->E, 1);
    mpz_set_str(curve->lambda, "5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72", 16);

    // Base du réseau GLV
    z2_init(&curve->v1);
    z2_init(&curve->v2);
    glv_basis(&curve->v1, &curve->v2, curve->n, curve->lambda);
}

/*Les valeurs ont été calculées avec le fichier gen_const_example2.sage*/
void init_example2_curve(GLVCurve *curve) {
    // Initialisation de la courbe
    mpz_inits(curve->E.a, curve->E.b, curve->E.p, curve->n, curve->lambda, NULL);
    mpz_set_str(curve->E.a, "2e818c97303c2c8e6ee49e2b6cacc754eff27e8346e4a23da9b7b882685b7c72", 16);
    mpz_set_ui(curve->E.b, 0);
    mpz_set_str(curve->E.p, "b4fc23d83418e4d099141c1a435cbb663817e03477f8f84f3afd51e63e89ef31", 16);
    mpz_set_str(curve->n, "b4fc23d83418e4d099141c1a435cbb67505e70fbe7dfa084bca52c64f59eca7a", 16);

    // Point générateur
    ec_point_proj_init(&curve->P);
    mpz_set_str(curve->P.X, "1e57bf13b5d9247da3cfdc0d46fe2b3924fd3698b201ab650ba511c9fdecf107", 16);
    mpz_set_str(curve->P.Y, "7aac691affa92ae823ad1af5cf523f0dce0c5f471527c1b4901b042fb3057409", 16);
    mpz_set_ui(curve->P.Z, 1);
    curve->P.infinity = 0;

    //trouver beta (une racine carré de -1 mod p)
    trouver_constantes_glv(curve->beta, &curve->E, 2);

    //trouver lambda (une racine carré de -1 mod n), calculé avec gen_const_example2.sage
    mpz_set_str(curve->lambda, "81e16b4d3131f1322cf0ab2ba439286b3962df578bf08b5dc04a37e1c8e31333", 16);

    // Base du réseau GLV
    z2_init(&curve->v1);
    z2_init(&curve->v2);
    glv_basis(&curve->v1, &curve->v2, curve->n, curve->lambda);
}


/*Les valeurs ont été calculées avec le fichier gen_const_example3.sage*/

void init_example3_curve(GLVCurve *curve){
    // Initialisation de la courbe
    mpz_inits(curve->E.a, curve->E.b, curve->E.p,
              curve->n, curve->lambda, curve->beta,
              curve->E.a2, NULL);
    mpz_t t;
    mpz_init_set_ui(t, 4);
    mpz_set_si(curve->E.a, -2);
    mpz_set_ui(curve->E.b, -1);
    mpz_set_str(curve->E.p, "e6adac14aa1890c61edfaeb1f66359f6468e064f93ce403cfb73b9e3cb7f44ad", 16);
    mpz_set_str(curve->n, "e6adac14aa1890c61edfaeb1f66359f59c5613554c28c9ab168f0b1f3c52fca4", 16);
    mpz_invert(curve->E.a2, t, curve->E.p);
    mpz_mul_ui(curve->E.a2, curve->E.a2, 3);
    mpz_neg(curve->E.a2, curve->E.a2);

    // Point générateur
    ec_point_proj_init(&curve->P);
    mpz_set_str(curve->P.X, "2958133073b80f31070226c37e132acccafc05892f2fd67faa556b8c5dacd8d", 16);
    mpz_set_str(curve->P.Y, "a92358ec0c20bdaf22e615ab89169b0158e4912b6e2c107b43c169cf773716bf", 16);
    mpz_set_ui(curve->P.Z, 1);
    curve->P.infinity = 0;

    // GLV paramètres
    trouver_constantes_glv(curve->beta, &curve->E, 3);

    // lambda 
    mpz_set_str(curve->lambda, "cec272884084085b2e7c660e7e5a27cc0d98b9741cba044bf0f30f059e313209", 16);

    // Base du réseau GLV
    z2_init(&curve->v1);
    z2_init(&curve->v2);
    glv_basis(&curve->v1, &curve->v2, curve->n, curve->lambda);

    mpz_clear(t);
}

void clear_curve(GLVCurve *curve){
    mpz_clears(curve->beta, curve->lambda, curve->n, NULL);
    z2_clear(&curve->v1);
    z2_clear(&curve->v2);
    ec_curve_clear(&curve->E);
    ec_point_proj_clear(&curve->P);
}