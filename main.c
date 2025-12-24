#include <stdio.h>
#include "EC_struct.h"
#include "EC_square_and_multiply_proj.h"
#include "EC_square_and_multiply_affine.h"

int main() {
    // --- 1. Définir une petite courbe pour test ---
    ECCurve E;
    mpz_init_set_str(E.p, "9739", 10);  // petit premier pour test
    mpz_init_set_ui(E.a, 497);
    mpz_init_set_ui(E.b, 1768);

    // --- 2. Définir un point affine ---
    ECPointAffine Pa;
    ec_point_affine_init(&Pa);
    mpz_set_ui(Pa.x, 493);
    mpz_set_ui(Pa.y, 5564);
    Pa.infinity = 0;

    // --- 3. Convertir en projectif ---
    ECPointProj Pp;
    ec_point_proj_init(&Pp);
    affine_to_proj(&Pp, &Pa);

    // --- 4. Tester doublement ---
    ECPointProj Rdouble;
    ec_point_proj_init(&Rdouble);
    ec_point_double_proj(&Rdouble, &Pp, &E);

    ECPointAffine Rdouble_aff;
    ec_point_affine_init(&Rdouble_aff);
    proj_to_affine(&Rdouble_aff, &Rdouble, &E);

    gmp_printf("2P = (%Zd, %Zd)\n", Rdouble_aff.x, Rdouble_aff.y);

    // --- 5. Tester addition ---
    ECPointProj Radd;
    ec_point_proj_init(&Radd);
    ec_point_add_proj(&Radd, &Pp, &Pp, &E);  // P + P

    ECPointAffine Radd_aff;
    ec_point_affine_init(&Radd_aff);
    proj_to_affine(&Radd_aff, &Radd, &E);

    gmp_printf("P + P = (%Zd, %Zd)\n", Radd_aff.x, Radd_aff.y);

    // --- 6. Tester multiplication scalaire ---
    ECPointProj Rmul;
    ec_point_proj_init(&Rmul);
    mpz_t k;
    mpz_init_set_ui(k, 20);

    ec_scalar_mul_proj(&Rmul, &Pp, k, &E);

    ECPointAffine Rmul_aff;
    ec_point_affine_init(&Rmul_aff);
    proj_to_affine(&Rmul_aff, &Rmul, &E);

    gmp_printf("20P = (%Zd, %Zd)\n", Rmul_aff.x, Rmul_aff.y);

    // --- 7. Nettoyage ---
    ec_point_affine_clear(&Pa);
    ec_point_proj_clear(&Pp);
    ec_point_proj_clear(&Rdouble);
    ec_point_affine_clear(&Rdouble_aff);
    ec_point_proj_clear(&Radd);
    ec_point_affine_clear(&Radd_aff);
    ec_point_proj_clear(&Rmul);
    ec_point_affine_clear(&Rmul_aff);
    mpz_clear(k);
    mpz_clear(E.p);
    mpz_clear(E.a);
    mpz_clear(E.b);

    return 0;
}
