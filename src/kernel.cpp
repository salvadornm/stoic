
#include "kernel.h"

//define cubic SPH Kernel where r is distance btwn them
double Wab(double r)
{
    r /= H;
    //std::cout << "r/h = " << r << std::endl;
    if (r < 1.0)
        return (1.0 - (3.0/2.0)*r*r + (3.0/4.0)*r*r*r)*a2;
    else if (r < 2.0)
        return (0.25*(2.0 - r)*(2.0 - r)*(2.0 - r))*a2;
    else
        return 0.0;
}

// define gradient of cubic kernel function
double DWab(Point<3,double> & dx, Point<3,double> & DW, double r, bool print)
{
    r /= H;
    double factor = 0;
    //std::cout << "r/h = " << r << std::endl;
    if (r < 1.0)
        factor = (c1*r + d1*r*r);  //*(a2/H));
    else if (r < 2.0)
        factor = (c2*(2.0 - r)*(2.0 - r));  //*(a2/H);
    else
        factor = 0.0;

    DW.get(0) = factor * dx.get(0);
    DW.get(1) = factor * dx.get(1);
    DW.get(2) = factor * dx.get(2);
    return factor;
}

