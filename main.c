#include <stdio.h>

#include "EC_square_and_multiply_proj.h"
#include "glv_curves.h"
#include "glv_acceleration.h"
#include "EC_DH.h"
#include "EC_GLV_demo.h"

enum {
    kCurveMinChoice = 1,
    kCurveMaxChoice = 3,
    kActionMinChoice = 1,
    kActionMaxChoice = 3
};

static int read_choice(const char *prompt, int min_choice, int max_choice) {
    int choice = 0;

    for (;;) {
        printf("%s", prompt);
        if (scanf("%d", &choice) == 1 && choice >= min_choice && choice <= max_choice) {
            return choice;
        }
        printf("Veuillez entrer un nombre entre %d et %d.\n", min_choice, max_choice);
        while (getchar() != '\n') {
        }
    }
}


static int init_curve_from_choice(GLVCurve *curve, int choice) {
    switch (choice) {
        case 1:
            init_secp256k1_curve(curve);
            return 1;
        case 2:
            init_example2_curve(curve);
            return 1;
        case 3:
            init_example3_curve(curve);
            return 1;
        default:
            return 0;
    }
}

int main() {
    const char *curve_prompt =
        "Entrez un nombre :\n"
        "1 pour faire GLV sur secp256k1, équation de la forme y^2 = x^3 + b avec b = 7\n"
        "2 pour faire GLV sur une équation de la forme y^2 = x^3 + a*x avec a un entier "
        "aléatoire mod p\n"
        "3 pour faire GLV sur une équation de la forme y^2 = x^3 - 3/4*x^2 - 2*x - 1\n";
    int type = read_choice(curve_prompt, kCurveMinChoice, kCurveMaxChoice);

    // ======================== TEST de la fonction GLV ========================
    GLVCurve C;
    if (!init_curve_from_choice(&C, type)) {
        printf("Veuillez entrer un nombre entre 1 et 3\n");
        return 0;
    }

    const char *action_prompt =
        "Entrez un nombre :\n"
        "1 pour utiliser GLV sur 1000 entiers aléatoires et constater l'accélération\n"
        "2 pour utiliser GLV sur un exemple spécifique\n"
        "3 pour une génération de clé Diffie-Hellmann\n";
    int option = read_choice(action_prompt, kActionMinChoice, kActionMaxChoice);

    switch (option) {
        case 1:
            const int N_tests = 1000;
            glv_acceleration(&C, N_tests);
            break;
        case 2:
            EC_GLV_demo(&C);
            break;
        case 3:
            ec_dh(&C);
            break;
        default:
            break;
    }

    // ----------  Nettoyage ----------
    
    clear_curve(&C);

    return 0;
}