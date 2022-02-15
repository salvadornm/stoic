
#include "kernel.h"

const double H = 0.0147224318643;
//define cubic SPH Kernel
const double a2 = 1.0;// 1.0/M_PI/H/H/H;
double Wab(double r)
{
    r /= H;
    if (r < 1.0)
        return (1.0 - 3.0/2.0*r*r + 3.0/4.0*r*r*r)*a2;
    else if (r < 2.0)
        return (0.25*(2.0 - r)*(2.0 - r)*(2.0 - r))*a2;
    else
        return 0.0;
}

// define gradient of cubic kernel function
const double c1 = -3.0/M_PI/H/H/H/H;
const double d1 = 9.0/4.0/M_PI/H/H/H/H;
const double c2 = -3.0/4.0/M_PI/H/H/H/H;
const double a2_4 = 0.25*a2;


void DWab(Point<3,double> & dx, Point<3,double> & DW, double r, bool print)
{
    const double qq=r/H;
    double qq2 = qq * qq;
    double fac1 = (c1*qq + d1*qq2)/r;
    double b1 = (qq < 1.0)?1.0f:0.0f;
    double wqq = (2.0 - qq);
    double fac2 = c2 * wqq * wqq / r;
    double b2 = (qq >= 1.0 && qq < 2.0)?1.0f:0.0f;
    double factor = (b1*fac1 + b2*fac2);
    DW.get(0) = factor * dx.get(0);
    DW.get(1) = factor * dx.get(1);
    DW.get(2) = factor * dx.get(2);
}
