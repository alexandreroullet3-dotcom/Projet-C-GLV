p = random_prime(2^256, lbound=2^255)
while not is_square(-7, p):
    p = random_prime(2^256, lbound=2^255)

F = GF(p)
E = EllipticCurve(F, [0, 0, -3/4, -2, -1])

n = E.order()

P = E.random_point()
while P.order() != n:
    P = E.random_point()

print("p =", hex(p))
print("n =", hex(n))
print("P=", hex(P[0]), hex(P[1]), hex(P[2]))
print(n.factor())