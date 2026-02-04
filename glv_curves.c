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
    ec_point_proj_init(&curve->phiP);
    mpz_set_str(curve->P.X, "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_set_str(curve->P.Y, "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);
    mpz_set_ui(curve->P.Z, 1);
    curve->P.infinity = 0;

    // GLV paramètres
    mpz_inits(curve->lambda, curve->beta, NULL);
    trouver_constantes_glv(curve->beta, &curve->E, 1);
    mpz_set_str(curve->lambda, "5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72", 16);
    ec_scalar_mul_proj(&curve->phiP, &curve->P, curve->lambda, &curve->E);

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
    ec_point_proj_init(&curve->phiP);
    mpz_set_str(curve->P.X, "a3cf5403796bcb5f78979ac1eac261d7d7c5785a712d01a52519542c2e6f0dc6", 16);
    mpz_set_str(curve->P.Y, "5a4073851f2f574df70c8c84906546b6ca76421fcc7db0901ffc77fa9a2ed540", 16);
    mpz_set_ui(curve->P.Z, 1);
    curve->P.infinity = 0;

    //trouver beta (une racine carré de -1 mod p)
    trouver_constantes_glv(curve->beta, &curve->E, 2);

    //trouver lambda (une racine carré de -1 mod n), calculé avec gen_const_example2.sage
    mpz_set_str(curve->lambda, "81e16b4d3131f1322cf0ab2ba439286b3962df578bf08b5dc04a37e1c8e31333", 16);
    ec_scalar_mul_proj(&curve->phiP, &curve->P, curve->lambda, &curve->E);

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
    mpz_set_str(curve->E.p, "c1da32bc97dae988c398d3d419ce5210634811ec6bd31474803febcf727f2129", 16);
    mpz_set_str(curve->n, "c1da32bc97dae988c398d3d419ce5211028b5b2bf6a01e44aa7fa1a81f9a90c4", 16);

    /*mpz_t t;
    mpz_init_set_si(t, -16);
    mpz_invert(t, t, curve->E.p);
    mpz_mul_ui(t, t, 35);
    mpz_mod(t, t, curve->E.p);
    mpz_set(curve->E.a, t);
    mpz_set_si(t, -64);
    mpz_invert(t, t, curve->E.p);
    mpz_mul_ui(t, t, 49);
    mpz_mod(t, t, curve->E.p);
    mpz_set(curve->E.b, t);*/
    


    
    mpz_t t;
    mpz_init_set_ui(t, 4);
    mpz_set_si(curve->E.a, -2);
    mpz_set_si(curve->E.b, -1);
    mpz_invert(curve->E.a2, t, curve->E.p);
    mpz_mul_ui(curve->E.a2, curve->E.a2, 3);
    mpz_neg(curve->E.a2, curve->E.a2);

    // Point générateur
    ec_point_proj_init(&curve->P);
    ec_point_proj_init(&curve->phiP);
    mpz_set_str(curve->P.X, "46d91a4b0d3f776d93bbac95a9754fa097cc72296faa177698d7f78d968c9399", 16);
    mpz_set_str(curve->P.Y, "7b817bd75ab53007dd2e82ce2a8f792204df0bd46995d5d31c355fb86d515431", 16);
    mpz_set_ui(curve->P.Z, 1);
    curve->P.infinity = 0;
    mpz_clear(t);

    // GLV paramètres
    trouver_constantes_glv(curve->beta, &curve->E, 3);

    // lambda 
    mpz_set_str(curve->lambda, "0x81c1ea47be450c44006169334b5c00996fca8415616e952dbca4e83b126568aa", 16);
    ec_scalar_mul_proj(&curve->phiP, &curve->P, curve->lambda, &curve->E);
    
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
    ec_point_proj_clear(&curve->phiP);
}