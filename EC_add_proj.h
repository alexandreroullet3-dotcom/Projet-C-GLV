#ifndef EC_ADD_PROJ_H
#define EC_ADD_PROJ_H 
#include "EC_struct.h"

void ec_point_double_proj(ECPointProj *R, const ECPointProj *P, const ECCurve *E);

void ec_point_add_proj(ECPointProj *R, const ECPointProj *P, const ECPointProj *Q, const ECCurve *E);

#endif