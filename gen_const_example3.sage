i = 1
while i:
    p = random_prime(2^256, lbound=2^255)
    if kronecker(-7, p) == -1:
        continue
    F = GF(p)
    E = EllipticCurve(F, [0, -3/4, 0, -2, -1])

    n = E.order()
    fac = n.factor()

    # on cherche un facteur premier ~256 bits
    for q,e in fac:
        if q.nbits() > 250 and q.is_prime():
            r = q
            h = n // r
            if h <= 16  :
                print("GOOD CURVE FOUND")
                print("p =", hex(p))
                print("r =", hex(r))
                print("h =", h)
                i = 0


print("n =", hex(n))
print(fac)

P = E.random_point()
print("n//P.order() ", hex(n//P.order()))
print("P =", hex(P[0]), hex(P[1]), hex(P[2]))

r = fac[1][0]
R = Integers(r)
s = R(-7).sqrt(all=True, extend=False)[0]  # racine de -7 mod p^e
inv2 = R(2)**-1                         # inverse de 2 modulo p^e
x = (1 + s) * inv2
print("lambda = ", hex(x))

R = Integers(p)
t = R(-7).sqrt(all=True, extend=False)[0]  # racine de -7 mod p^e
inv2 = R(2)**-1                         # inverse de 2 modulo p^e
y = (1 + t) * inv2
print("beta = ", hex(y))