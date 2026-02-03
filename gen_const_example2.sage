p = random_prime(2^256, lbound=2^255)
while mod(p, 4) != 1:
    p = random_prime(2^256, lbound=2^255)
p = 0xb4fc23d83418e4d099141c1a435cbb663817e03477f8f84f3afd51e63e89ef31
F = GF(p)

a = int(F.random_element())
a = 0x2e818c97303c2c8e6ee49e2b6cacc754eff27e8346e4a23da9b7b882685b7c72
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