#include <time.h>

#include "EC_square_and_multiply_proj.h" // Les additions sont deja dans square and multiply
#include "EC_GLV.h"
#include "quadratic_solver.h"
#include "glv_curves.h"


int main() {
    
//////////////////  TEST de la fonction GLV //////////////////////////
    
    GLVCurve C;
    init_secp256k1_curve(&C);
    ECPointProj P, R_classic, R_glv;
    ec_point_proj_init(&P);
    ec_point_proj_init(&R_classic);
    ec_point_proj_init(&R_glv);
    mpz_t beta, lambda, n, p;
    mpz_init_set(beta, C.beta);
    mpz_init_set(lambda, C.lambda);
    mpz_init_set(n, C.n);

    Z2 v1, v2;
    z2_init(&v1);
    z2_init(&v2);
    mpz_set(v1.x, C.v1.x);
    mpz_set(v1.y, C.v1.y);
    mpz_set(v2.x, C.v2.x);
    mpz_set(v2.y, C.v2.y);
    ECCurve E;
    ec_curve_init(&E);
    mpz_set(E.a, C.E.a);
    mpz_set(E.b, C.E.b);
    mpz_set(E.p, C.E.p);
    mpz_init_set(p, E.p);
    /*
    // 1. Initialisation des variables

    mpz_t p, n, lambda, beta;
    mpz_inits(p, n, lambda, beta, NULL);

    // 2. Initialisation des Points et Courbe

    ECCurve E;
    mpz_inits(E.a, E.b, E.p, NULL);

    ECPointProj P, R_classic, R_glv;
    ec_point_proj_init(&P);
    ec_point_proj_init(&R_classic);
    ec_point_proj_init(&R_glv);


    // 3. PARAMÈTRES courbe SECP256K1 & GLV

    mpz_set_str(p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    mpz_set(E.p, p);
    mpz_set_ui(E.a, 0);
    mpz_set_ui(E.b, 7);
    mpz_set_str(n, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16);

    // Point P
    mpz_set_str(P.X, "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_set_str(P.Y, "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);
    mpz_set_ui(P.Z, 1);
    P.infinity = 0;

    // Paramètres GLV
    mpz_set_str(beta, "7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee", 16);
    mpz_set_str(lambda, "5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72", 16);
    //mpz_set_str(x1, "3086d221a7d46bcde86c90e49284eb15", 16);
    //mpz_set_str(y1, "e4437ed6010e88286f547fa90abfe4c3", 16);
    //mpz_neg(y1, y1); 
    //mpz_set_str(x2, "114ca50f7a8e2f3f657c1108d9d44cfd8", 16); 
    //mpz_set_str(y2, "3086d221a7d46bcde86c90e49284eb15", 16);

    Z2 v1, v2;
    z2_init(&v1);
    z2_init(&v2);
    glv_basis(&v1, &v2, n, lambda);
    gmp_printf("v1: \n x:%Zx \n y:%Zx\n",v1.x, v1.y);
    gmp_printf("v2: \n x:%Zx \n y:%Zx\n",v2.x, v2.y);*/
    

    // ---------- 2. Boucle de tests ----------
    const int N_TESTS = 1000;
    mpz_t k;
    mpz_init(k);

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
        mpz_urandomm(k, state, n);

        // -- Multiplication classique --
        start_classic = clock();
        ec_scalar_mul_proj(&R_classic, &P, k, &E);
        end_classic = clock();
        time_classic += (double)(end_classic - start_classic) / CLOCKS_PER_SEC;

        // -- Multiplication GLV --
        start_glv = clock();
        ec_scal_mul_glv(&R_glv, &P, k, &E, &v1, &v2, beta, n);
        end_glv = clock();
        time_glv += (double)(end_glv - start_glv) / CLOCKS_PER_SEC;

        // Vérification facultative
        if (!ec_cmp_proj(&R_classic, &R_glv, &E)) {
            printf("Mismatch à l'itération %d !\n", i);
            break;
        }
    }

    // ---------- 3. Résultat ----------
    printf("\nTemps total pour %d multiplications scalaires :\n", N_TESTS);
    printf(" - Méthode classique : %.6f s\n", time_classic);
    printf(" - Méthode GLV      : %.6f s\n", time_glv);
    printf(" - Gain GLV         : %.2f x\n", time_classic / time_glv);

    // ---------- 4. Nettoyage ----------
    mpz_clears(p, n, lambda, beta, k, state, NULL);
    ec_point_proj_clear(&P);
    ec_point_proj_clear(&R_classic);
    ec_point_proj_clear(&R_glv);
    z2_clear(&v1);
    z2_clear(&v2);
    ec_curve_clear(&E);

    return 0;

}