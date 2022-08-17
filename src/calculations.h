#ifndef _calculations_h
#define _calculations_h

#include "global.h"

void updateEqtnState(particleset & vd);
void updateThermalProperties2(particleset & vd, int a);
void updateThermalProperties1(particleset & vd, int a);
double viscous(const Point<3,double> & dr, double rr2, Point<3,double> & dv, double rhoa, double rhob, double massb, double visc); //double & visc)
double calculateMass(particleset & vd, engine eng);


#endif // _calculations_h