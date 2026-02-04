#include "EC_square_and_multiply_proj.h"
#include "EC_GLV.h"
#include "quadratic_solver.h"
#include "glv_curves.h"
#include "glv_acceleration.h"
#include "EC_DH.h"
#include "EC_GLV_mpz_t.h"


int main() {
    int type;
    printf("Entrez un nombre :\n1 pour faire GLV sur secp256k1, équation de la forme y^2 = x^3 + b avec b = 7\n");
    printf("2 pour faire GLV sur une équation de la forme y^2 = x^3 + a*x avec a un entier aléatoire mod p\n");
    printf("3 pour faire GLV sur une équation de la forme y^2 = x^3 - 3/4*x^2 - 2*x - 1\n");
    scanf("%d", &type);
    
    //////////////////  TEST de la fonction GLV //////////////////////////
    GLVCurve C;
    if (type == 1){
        init_secp256k1_curve(&C);
        printf("Paramètres publique :\n");
        gmp_printf("\n Nombre premier utilisé:\np = %Zx\nOrdre de la courbe :\nn = %Zx\nPoint générateur :\nP = (%Zx,%Zx)\n\n",C.E.p, C.n, C.P.X, C.P.Y);
    }
    else if (type == 2){
        init_example2_curve(&C);
        printf("Paramètres publique :\n");
        gmp_printf("\nNombre premier utilisé:\np = %Zx\nOrdre de la courbe :\nn = %Zx\nPoint générateur :\nP = (%Zx,%Zx)\n\n",C.E.p, C.n, C.P.X, C.P.Y);
    }
    else if (type == 3){
        init_example3_curve(&C);
        printf("Paramètres publique :\n");
        gmp_printf("\nNombre premier utilisé:\np = %Zx\nOrdre de la courbe :\nn = %Zx\nPoint générateur :\nP = (%Zx,%Zx)\n\n",C.E.p, C.n, C.P.X, C.P.Y);
    }
    else{
        printf("Veuillez entrer un nombre entre 1 et 3\n");
        return 0;
    }
    ECPointProj P, phiP, R_classic, R_glv;
    ec_point_proj_init(&P);
    ec_point_proj_init(&phiP);
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
    mpz_set(E.a2, C.E.a2);
    mpz_init_set(p, E.p);
    ec_point_proj_copy(&P, &C.P);
    ec_point_proj_copy(&phiP, &C.phiP);

    // ==========================================
    //      BLOC DE DIAGNOSTIC DE SÉCURITÉ
    // ==========================================
    printf("\n--- Vérification de la cohérence GLV ---\n");
    
    ECPointProj test_phi; // <--- C'est cette ligne qui manquait !
    ec_point_proj_init(&test_phi);
    
    // Test 1 : Est-ce que phi(P) == [lambda]P ?
    ec_scalar_mul_proj(&test_phi, &P, lambda, &E);
    
    if (ec_cmp_proj(&test_phi, &phiP, &E) == 0) {
        printf("[OK] : phi(P) == [lambda]P. L'endomorphisme est parfait.\n");
    } else {
        // Test 2 : Est-ce que phi(P) == [1 - lambda]P ?
        mpz_t other_lambda;
        mpz_init_set_ui(other_lambda, 1);
        mpz_sub(other_lambda, other_lambda, lambda);
        mpz_mod(other_lambda, other_lambda, n);
        
        ec_scalar_mul_proj(&test_phi, &P, other_lambda, &E);
        
        if (ec_cmp_proj(&test_phi, &phiP, &E) == 0) {
            printf("[ALERTE] : phi(P) == [1 - lambda]P.\n");
            printf("-> Action : Remplace lambda par (1 - lambda) mod n dans glv_curves.c\n");
        } else {
            // Test 3 : Est-ce que phi(P) == [-lambda]P ?
            ec_point_proj_neg(&test_phi, &test_phi); // On ne fait que tester le signe
            ec_scalar_mul_proj(&test_phi, &P, lambda, &E);
            ec_point_proj_neg(&test_phi, &test_phi);
            
            if (ec_cmp_proj(&test_phi, &phiP, &E) == 0) {
                printf("[ALERTE] : phi(P) == [-lambda]P.\n");
                printf("-> Action : Verifie le signe de l'ordonnee dans ton endomorphisme.\n");
            } else {
                printf("[ERREUR CRITIQUE] : phi(P) ne correspond a rien.\n");
                printf("-> Action : Verifie a2 dans EC_add_proj.c et ta formule de phi3.\n");
            }
        }
        mpz_clear(other_lambda);
    }
    ec_point_proj_clear(&test_phi);
    printf("----------------------------------------\n\n");

    
    printf("Entrez un nombre :\n");
    printf("1 pour utiliser GLV sur 1000 entiers aléatoires et constater l'accélération\n");
    printf("2 pour utiliser GLV sur un exemple spécifique\n");
    printf("3 pour une génération de clé Diffie-Hellmann\n");
    int option;
    scanf("%d", &option);

    if (option == 1){
        glv_acceleration(&C);
    }

    if (option == 2){
        EC_GLV_mpz_t(&C);
    }

    if (option == 3){
        ec_dh(&C);
    }

    // ---------- 4. Nettoyage ----------
    mpz_clears(p, n, lambda, beta, NULL);
    ec_point_proj_clear(&P);
    ec_point_proj_clear(&phiP);
    ec_point_proj_clear(&R_classic);
    ec_point_proj_clear(&R_glv);
    z2_clear(&v1);
    z2_clear(&v2);
    ec_curve_clear(&E);
    clear_curve(&C);

    return 0;
}