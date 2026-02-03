p = random_prime(2^256, lbound=2^255)
while mod(p, 4) != 1:
    p = random_prime(2^256, lbound=2^255)

F = GF(p)

a = int(F.random_element())
E = EllipticCurve(F, [a, 0])

n = E.order()

P = E.random_point()
while P.order() != n:
    P = E.random_point()

print("a =", hex(a))
print("p =", hex(p))
print("n =", hex(n))
print("P =", hex(P[0]),hex(P[1]), hex(P[2]))
fac = n.factor()
print(fac)

L_moduli = []
L_roots = []

for prime, e in fac:
    modulus = prime**e
    R = Integers(modulus)
    
    if prime == 2:
        # pour 2^e, on teste toutes les valeurs possibles pour x^2+1 mod 2^e
        roots = [x for x in range(modulus) if (x**2 + 1) % modulus == 0]
        L_moduli.append(modulus)
        L_roots.append(roots[0])  # on prend une solution parmi celles possibles
        continue
    
    # pour les autres nombres premiers
    x = R(-1).sqrt(all=True, extend=False)[0]  # racine de -1 mod p^e
    L_moduli.append(modulus)
    L_roots.append(Integer(x))

# Combiner toutes les solutions avec CRT
racine = crt(L_roots, L_moduli)
print("lambda = ", hex(racine))