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
print("P=", P)
print("P=", hex(P[0]),hex(P[1]), hex(P[2]))
print(n.factor())