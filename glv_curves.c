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

    ECPointAffine P_aff;
    ec_point_affine_init(&P_aff);
    proj_to_affine(&P_aff, &curve->P,&curve->E);

    // GLV paramètres
    mpz_inits(curve->lambda, curve->beta, NULL);
    trouver_constantes_glv(curve->beta, curve->lambda, &curve->E, curve->n, &P_aff, 1);
    mpz_set_str(curve->lambda, "5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72", 16);

    // Base du réseau GLV
    z2_init(&curve->v1);
    z2_init(&curve->v2);
    glv_basis(&curve->v1, &curve->v2, curve->n, curve->lambda);

    ec_point_affine_clear(&P_aff);
}

/*Les valeurs ont été calculées avec le fichier gen_const_example2.sage*/
void init_example2_curve(GLVCurve *curve) {
    // Initialisation de la courbe
    mpz_inits(curve->E.a, curve->E.b, curve->E.p, curve->n, NULL);
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

    ECPointAffine P_aff;
    ec_point_affine_init(&P_aff);
    proj_to_affine(&P_aff, &curve->P,&curve->E);

    //On calcule ici lambda à la main avec le crt en utilisant la factorisation de n donné par sage
    // n = 2 * 184239961 * 3248770255115366281 * 868559398343199438216463774116001181
    //La fonction mpz_tonelli_shanks ne marche que pour n premier
    mpz_inits(curve->lambda, curve->beta, NULL);

    mpz_t t, a_1, a_2, a_3, a_4, b, c, d, e;
    mpz_inits(t, a_1, a_2, a_3, a_4, b, c, d, e, NULL);

    mpz_set_ui(b, 2);
    mpz_set_ui(c, 184239961);
    mpz_set_ui(d, 3248770255115366281);
    mpz_set_str(e, "868559398343199438216463774116001181", 10);

    mpz_set_ui(a_1, 1);

    mpz_sub_ui(t, c, 1);    
    mpz_tonelli_shanks(a_2, t, c);

    mpz_sub_ui(t, d, 1);
    mpz_tonelli_shanks(a_3, t, d);

    mpz_sub_ui(t, e, 1);
    mpz_tonelli_shanks(a_4, t, e);

    mpz_t n_1, n_2, n_3, n_4, e_1, e_2, e_3, e_4;
    mpz_inits(n_1, n_2, n_3, n_4, e_1, e_2, e_3, e_4, NULL);

    mpz_divexact(n_1, curve->n, b);
    mpz_divexact(n_2, curve->n, c);
    mpz_divexact(n_3, curve->n, d);
    mpz_divexact(n_4, curve->n, e);

    mpz_invert(e_1, b, n_1);
    mpz_invert(e_2, c, n_2);
    mpz_invert(e_3, d, n_3);
    mpz_invert(e_4, e, n_4);

    mpz_mul(a_1, a_1, n_1);
    mpz_mul(a_1, a_1, e_1);

    mpz_mul(a_2, a_2, n_2);
    mpz_mul(a_2, a_2, e_2);

    mpz_mul(a_3, a_3, n_3);
    mpz_mul(a_3, a_3, e_3);

    mpz_mul(a_4, a_4, n_4);
    mpz_mul(a_4, a_4, e_4);

    mpz_add(t, a_1, a_2);
    mpz_add(t, t, a_3);
    mpz_add(t, t, a_4);

    mpz_mod(curve->lambda, t, curve->n);


    //Constantes GLV
    //jai pas modifié la signature de cette fonction mais elle ne modifie pas lambda

    trouver_constantes_glv(curve->beta, curve->lambda, &curve->E, curve->n, &P_aff, 2);
    
    gmp_printf("lambda = %Zx\n", curve->lambda);
    gmp_printf("beta = %Zx\n", curve->beta);

    // Base du réseau GLV
    z2_init(&curve->v1);
    z2_init(&curve->v2);
    glv_basis(&curve->v1, &curve->v2, curve->n, curve->lambda);

    mpz_clears(t, a_1, a_2, a_3, a_4, n_1, n_2, n_3, n_4, e_1, e_2, e_3, e_4, b, c, d, e, NULL);
}



/*Les valeurs ont été calculées avec le fichier gen_const_example3.sage*/

void init_example3_curve(GLVCurve *curve){
    // Initialisation de la courbe
    mpz_inits(curve->E.a, curve->E.b, curve->E.p, curve->n, NULL);
    
    //probleme: nous on a que des courbes y^3 = x^3 +ax + b
    //mais dans cet exemple on a coefficient x^2
    //jsp comment faire, si faut changer ECCurve en ajoutant un champ a2 par exemple, mais jsp si
    //faudrait tout changer du coup?
    
    //mpz_set_ui(curve->E.a, 0);
    //mpz_set_ui(curve->E.b, 7);
    mpz_set_str(curve->E.p, "910a26fd7f781ab40c311a6d5ca610036508d3f5e9899d9d49f5a6cccd103a6d", 16);
    mpz_set_str(curve->n, "910a26fd7f781ab40c311a6d5ca610038b33954211fffaf5d82f7a5b9f698092", 16);

    // Point générateur
    ec_point_proj_init(&curve->P);
    mpz_set_str(curve->P.X, "3F2D711A892A5B5074A60ED58BF5CBD13D9334A7F4B2E5265DB5DBE0C851B3B5", 16);
    mpz_set_str(curve->P.Y, "503AD2618C647B7D52BDCA4656350F71EEC2271168C2C9C09885794C2ADBB035", 16);
    mpz_set_ui(curve->P.Z, 1);
    curve->P.infinity = 0;

    ECPointAffine P_aff;
    ec_point_affine_init(&P_aff);
    proj_to_affine(&P_aff, &curve->P,&curve->E);

    // GLV paramètres
    mpz_inits(curve->lambda, curve->beta, NULL);
    trouver_constantes_glv(curve->beta, curve->lambda, &curve->E, curve->n, &P_aff, 3);

    // Base du réseau GLV
    z2_init(&curve->v1);
    z2_init(&curve->v2);
    glv_basis(&curve->v1, &curve->v2, curve->n, curve->lambda);


    ec_point_affine_clear(&P_aff);
}