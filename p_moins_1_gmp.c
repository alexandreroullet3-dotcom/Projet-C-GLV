#include <gmp.h>
#include <stdio.h>
void factoriser_rsa(mpz_t result, const mpz_t N){
    mpz_t a;
    mpz_init_set_ui(a, 2); // choix de la base arbitraire, 2 marche
    mpz_t temp;
    mpz_init(temp);
    for (int i = 2; i < 10000; i++){ // le 10000 est arbitraire, on considère la base de factorisation des premiers <10000
        mpz_sub_ui(temp, a, 1);
        mpz_gcd(result, temp, N); // r = pcgd(a^k!-1, N)
        if ((mpz_cmp_ui(result, 1) != 0) && (mpz_cmp(result, N) != 0)){//si on trouve un facteur non trivial
            mpz_clear(temp);
            return ;
        }
        mpz_powm_ui(a, a, i, N); // a = a^i mod N, ça donne donc a^(i!) mod N puisqu'on itère 
    }
    mpz_clear(temp);
    return;
}


int main(){
    mpz_t N, result, temp;
    mpz_init_set_str(N, "199214358783833785496649131630759414803916321139456200129431155042143170897974614023327", 10);
    mpz_init(result);
    mpz_init(temp);
    factoriser_rsa(result, N);
    mpz_out_str(stdout, 10, result);
    printf("\n");
    mpz_cdiv_q(temp, N, result);
    mpz_out_str(stdout, 10, temp);
    printf("\n");
    mpz_mul(temp, temp, result);
    printf("%d\n",mpz_cmp(N, temp));
}