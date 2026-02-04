#include "EC_DH.h"

void ec_dh(GLVCurve *curve){
    mpz_t k;
    mpz_init(k);
    ECPointAffine Raff;
    ECPointProj R_glv;
    ec_point_affine_init(&Raff);
    ec_point_proj_init(&R_glv);
    srand(time(NULL));
    gmp_randstate_t state;
    gmp_randinit_default(state); // Initialise le générateur par défaut

    // On peut le "seeder" avec le temps ou un autre nombre aléatoire
    unsigned long seed = time(NULL);
    gmp_randseed_ui(state, seed);
    // k aléatoire modulo n
    mpz_urandomm(k, state, curve->n);
    ec_scal_mul_glv(&R_glv, &curve->P, &curve->phiP, k, &curve->E, &curve->v1, &curve->v2);
    proj_to_affine(&Raff, &R_glv, &curve->E);
    
    gmp_printf("Votre clé publique est :\n(%Zx, %Zx)\nVotre clé privé est :\n%Zx\n", Raff.x, Raff.y, k);
    gmp_randclear(state);
    mpz_clear(k);
    ec_point_affine_clear(&Raff);
    ec_point_proj_clear(&R_glv);
    }