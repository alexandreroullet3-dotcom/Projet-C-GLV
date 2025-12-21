#include "EC_struct.h"

void ec_point_add(ECPoint *R, const ECPoint *P, const ECPoint *Q, const ECCurve *E);

void ec_point_double(ECPoint *R, const ECPoint *P, const ECCurve *E);