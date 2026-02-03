p = random_prime(2^256, lbound=2^255)
while kronecker(-7, p) == -1:
    p = random_prime(2^256, lbound=2^255)

F = GF(p)
E = EllipticCurve(F, [0, 0, -3/4, -2, -1])

n = E.order()

print("p =", hex(p))
print("n =", hex(n))

fac = n.factor()
for p,e in fac:
    if p == 2:
        continue
    print(p, kronecker(-7, p))
print(fac)

P = E.random_point()
while P.order() != n:
    P = E.random_point()
print("P =", hex(P[0]), hex(P[1]), hex(P[2]))

L_moduli = []
L_roots = []

for prime, e in fac:
    modulus = prime**e
    R = Integers(modulus)

    if prime == 2:
        # pour 2^e, on teste toutes les valeurs possibles pour x^2-x+2 mod 2^e
        roots = [x for x in range(modulus) if (x**2 - x + 2) % modulus == 0]
        L_moduli.append(modulus)
        L_roots.append(roots[0])  # on prend une solution parmi celles possibles
        continue

    # pour les autres nombres premiers
    s = R(-7).sqrt(all=True, extend=False)[0]  # racine de -7 mod p^e
    inv2 = R(2)**-1                         # inverse de 2 modulo p^e
    x = (1 + s) * inv2                     # solution x = (1 + sqrt(-7))/2 mod p^e
    L_moduli.append(modulus)
    L_roots.append(Integer(x))

# Combiner toutes les solutions avec CRT
racine = crt(L_roots, L_moduli)
print("lambda = ", hex(racine))