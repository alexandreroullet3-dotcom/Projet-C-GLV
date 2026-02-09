#include "glv_curves.h"

#include <stdio.h>

#include "EC_endo_phi_GLV.h"

static void print_curve_params(const GLVCurve *curve) {
    printf("Paramètres publics :\n");
    gmp_printf(
        "\nNombre premier utilisé:\np = %Zx\nOrdre de la courbe :\nn = %Zx\n"
        "Point générateur :\nP = (%Zx,%Zx)\n\n",
        curve->E.p,
        curve->n,
        curve->P.X,
        curve->P.Y);
}


/*
 * Calcule phi(P) en affine puis reconvertit en projectif.
 */
void calculer_phiP_optimise(GLVCurve *curve, int type)
{
    ECPointAffine P_aff, phiP_aff;
    ec_point_affine_init(&P_aff);
    ec_point_affine_init(&phiP_aff);

    // Passage en affine.
    proj_to_affine(&P_aff, &curve->P, &curve->E);

    // Choix de l'endomorphisme selon le type de courbe.
    switch (type) {
        case 1:
            ec_endo_phi1_affine(&phiP_aff, &P_aff, &curve->E, curve->beta);
            break;
        case 2:
            ec_endo_phi2_affine(&phiP_aff, &P_aff, &curve->E, curve->beta);
            break;
        case 3:
            ec_endo_phi3_affine(&phiP_aff, &P_aff, &curve->E, curve->beta);
            break;
        default:
            break;
    }

    affine_to_proj(&curve->phiP, &phiP_aff);

    ec_point_affine_clear(&P_aff);
    ec_point_affine_clear(&phiP_aff);
}

/*
 * Valeurs tabulées (secp256k1).
 */
void init_secp256k1_curve(GLVCurve *curve)
{
    mpz_inits(curve->E.a, curve->E.b, curve->E.p, curve->n, curve->lambda, curve->beta, curve->E.a2, NULL);
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

    calculer_phiP_optimise(curve, 1);

    z2_init(&curve->v1);
    z2_init(&curve->v2);
    glv_basis(&curve->v1, &curve->v2, curve->n, curve->lambda);
    print_curve_params(curve);
}

/*
 * Valeurs calculées avec gen_const_example2.sage.
 */
void init_example2_curve(GLVCurve *curve)
{
    mpz_inits(curve->E.a, curve->E.b, curve->E.p, curve->n, curve->lambda, curve->beta, curve->E.a2, NULL);
    mpz_set_str(curve->E.a, "5d00518733d00bbff1eb8461f994c688e1230b843372b49515829a7ccbcaf04e", 16);
    mpz_set_ui(curve->E.b, 0);

    // Ici on a n = h*r avec h=2, mais on ne fait GLV que sur le sous-groupe d'ordre r.
    mpz_set_str(curve->E.p, "819b395e84dc993cb56563295dd41ddf3f3679adc076bcd2b50be89f583b511d", 16);
    mpz_set_str(curve->n, "40cd9caf426e4c9e5ab2b194aeea0ef03b5024b0514a3c3c205cb23214f5f259", 16);

    ec_point_proj_init(&curve->P);
    ec_point_proj_init(&curve->phiP);
    mpz_set_str(curve->P.X, "924916153b530a9d5b459742ec5b797c47d4d7bcd53f8b2b8735011cac52d1b", 16);
    mpz_set_str(curve->P.Y, "4068d45bb0b8e75b5ba5a22e340737d25acd6ae0570e2bf545f36d1d31cd709a", 16);
    mpz_set_ui(curve->P.Z, 1);
    curve->P.infinity = 0;

    // Initialiser lambda avant trouver_constantes_glv.
    mpz_set_str(curve->lambda, "60a9f11a57b6d81a918bccc547e18eb730747c1c7262ee63025e410caf7b6357", 16);
    trouver_constantes_glv(curve->beta, &curve->E, &curve->P, curve->lambda, 2);

    calculer_phiP_optimise(curve, 2);

    z2_init(&curve->v1);
    z2_init(&curve->v2);
    glv_basis(&curve->v1, &curve->v2, curve->n, curve->lambda);
    print_curve_params(curve);
}

/*
 * Valeurs calculées avec gen_const_example3.sage.
 */
void init_example3_curve(GLVCurve *curve)
{
    mpz_inits(curve->E.a, curve->E.b, curve->E.p, curve->n, curve->lambda, curve->beta, curve->E.a2, NULL);
    mpz_t t;
    mpz_init_set_ui(t, 4);

    mpz_set_si(curve->E.a, -2);
    mpz_set_si(curve->E.b, -1);
    
    // Ici on a n = h*r avec h=4, mais on ne fait GLV que sur le sous-groupe d'ordre r.
    mpz_set_str(curve->E.p, "ba41f30a2231820b3fa2957bae0b9fed1801dff14e9cded971efdb63aeb369d5", 16);
    mpz_set_str(curve->n, "2e907cc2888c6082cfe8a55eeb82e7faf521130fcfa0ca47fa6ecebf5bc4759b", 16);

    mpz_invert(curve->E.a2, t, curve->E.p);
    mpz_mul_ui(curve->E.a2, curve->E.a2, 3);
    mpz_neg(curve->E.a2, curve->E.a2);
    mpz_mod(curve->E.a2, curve->E.a2, curve->E.p);

    ec_point_proj_init(&curve->P);
    ec_point_proj_init(&curve->phiP);
    
    mpz_set_str(curve->P.X, "b926c0c1e0bb366767cef2ecb7b6b0363611073bb6caeec225040dbd6986d21c", 16);
    mpz_set_str(curve->P.Y, "583d94cf7256ce2c92dbd84aac314a0c699a4e3b94b75baf83d67cd5dd54d11d", 16);
    mpz_set_ui(curve->P.Z, 1);
    curve->P.infinity = 0;
    mpz_clear(t);

    mpz_set_str(curve->lambda, "1ea1f990009d0fe14b48af8f43d1a1bc32a0e228508b71c353c0be059898928", 16);
    mpz_set_str(curve->beta, "660bd4979c65fd51364f5530f1ef8e4ceaeaccdb85033f184764b3a1551c8c5c", 16);
    //trouver_constantes_glv(curve->beta, &curve->E, &curve->P, curve->lambda, 3);

    calculer_phiP_optimise(curve, 3);
    
    // Test de validation interne.
    ECPointProj test_phi;
    ec_point_proj_init(&test_phi);
    ec_scalar_mul_proj(&test_phi, &curve->P, curve->lambda, &curve->E);
    if (ec_cmp_proj(&test_phi, &curve->phiP, &curve->E) != 0) {
        printf("ALERTE : phi(P) != lambda*P ! Le probleme vient de beta ou de l'arithmetique.\n");
    }
    ec_point_proj_clear(&test_phi);

    z2_init(&curve->v1);
    z2_init(&curve->v2);
    glv_basis(&curve->v1, &curve->v2, curve->n, curve->lambda);
    print_curve_params(curve);

}

void clear_curve(GLVCurve *curve)
{   mpz_clears(curve->beta, curve->lambda, curve->n, NULL); 
    ec_curve_clear(&curve->E);
    ec_point_proj_clear(&curve->P);
    ec_point_proj_clear(&curve->phiP);
    z2_clear(&curve->v1);
    z2_clear(&curve->v2);
}
