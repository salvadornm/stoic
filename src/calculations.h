#ifndef _calculations_h
#define _calculations_h

#include "global.h"

void calcPressure(particleset & vd);
double viscous(const Point<3,double> & dr, double rr2, Point<3,double> & dv, double rhoa, double rhob, double massb, double & visc);


#endif // _calculations_h