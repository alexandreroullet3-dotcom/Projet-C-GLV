#include "EC_GLV_mpz_t.h"
void EC_GLV_mpz_t(const GLVCurve *curve){
    mpz_t k;
    mpz_init(k);
    ECPointAffine Raff;
    ec_point_affine_init(&Raff);
    ECPointProj R_glv, R_classic;
    ec_point_proj_init(&R_glv);
    ec_point_proj_init(&R_classic);
    
    printf("Entrez l'entier k, en hex, sans le 0x\n");
    char hex_str[256]; // plus grand pour éviter overflow
    if (scanf("%255s", hex_str) != 1) {
        printf("Erreur de lecture de k\n");
        return ;
    }

    if (mpz_set_str(k, hex_str, 16) != 0) {
        printf("Erreur : k n'est pas un hex valide !\n");
        return ;
        }

    ec_scal_mul_glv(&R_glv, &curve->P, &curve->phiP, k, &curve->E, &curve->v1, &curve->v2);

    ec_scalar_mul_proj(&R_classic, &curve->P, k, &curve->E);

    if (ec_cmp_proj(&R_classic, &R_glv, &curve->E) == 0){
        printf("Le calcul avec GLV coïncide avec l'exponentiation rapide et le résultat est \n");
        proj_to_affine(&Raff, &R_glv, &curve->E);
        gmp_printf("k*P = (%Zx, %Zx)\n", Raff.x, Raff.y);
    }
    else{
        printf("Les résultats GLV et exponentiation rapide ne coïncident pas\n");
    }
    mpz_clear(k);
    ec_point_affine_clear(&Raff);
    ec_point_proj_clear(&R_glv);
    ec_point_proj_clear(&R_classic);

}