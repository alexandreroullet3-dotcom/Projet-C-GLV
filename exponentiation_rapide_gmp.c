#include <gmp.h>
void mpz_pow_fast(mpz_t result, const mpz_t base, unsigned long exp) {
    mpz_set_ui(result, 1);      // result = 1
    mpz_t b;
    mpz_init_set(b, base);      // b = base

    while (exp > 0) {
        if (exp % 2 == 1) {     // si le bit courant est 1
            mpz_mul(result, result, b);  // result *= b
        }
        mpz_mul(b, b, b);       // b = b^2
        exp /= 2;               // passer au bit suivant
    }

    mpz_clear(b);
}