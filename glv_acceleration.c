#include "glv_acceleration.h"

void glv_acceleration(const GLVCurve *curve){
    // ---------- Boucle de tests ----------
    const int N_TESTS = 1000;
    mpz_t k;
    mpz_init(k);

    ECPointProj R_classic, R_glv;
    ec_point_proj_init(&R_classic);
    ec_point_proj_init(&R_glv);
    clock_t start_classic, end_classic;
    clock_t start_glv, end_glv;

    double time_classic = 0.0;
    double time_glv = 0.0;

    srand(time(NULL));
    gmp_randstate_t state;
    gmp_randinit_default(state); // Initialise le générateur par défaut

    // On peut le "seeder" avec le temps ou un autre nombre aléatoire
    unsigned long seed = time(NULL);
    gmp_randseed_ui(state, seed);

    for (int i = 0; i < N_TESTS; i++) {
        // k aléatoire modulo n
        mpz_urandomm(k, state, curve->n);

        // -- Multiplication classique --
        start_classic = clock();
        ec_scalar_mul_proj(&R_classic, &curve->P, k, &curve->E);
        end_classic = clock();
        time_classic += (double)(end_classic - start_classic) / CLOCKS_PER_SEC;

        // -- Multiplication GLV --
        start_glv = clock();
        ec_scal_mul_glv(&R_glv, &curve->P, &curve->phiP, k, &curve->E, &curve->v1, &curve->v2);
        end_glv = clock();
        time_glv += (double)(end_glv - start_glv) / CLOCKS_PER_SEC;

        // Vérification facultative
        if (ec_cmp_proj(&R_classic, &R_glv, &curve->E) == 1) {
            printf("Mismatch à l'itération %d !\n", i);
            break;
        }
    }

    // ---------- 3. Résultat ----------
    printf("\nTemps total pour %d multiplications scalaires sur des entiers aléatoires:\n", N_TESTS);
    printf(" - Méthode classique : %.6f s\n", time_classic);
    printf(" - Méthode GLV       : %.6f s\n", time_glv);
    printf(" - Gain GLV          : %.2f x\n", time_classic / time_glv);
    gmp_randclear(state);
    mpz_clear(k);
    ec_point_proj_clear(&R_classic);
    ec_point_proj_clear(&R_glv);
    }