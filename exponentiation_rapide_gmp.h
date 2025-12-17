/*
    Calcule result = base^exp
    base et result sont des mpz_t
    exp est un unsigned long (entier non négatif)
*/
void mpz_pow_fast(mpz_t result, const mpz_t base, unsigned long exp); 