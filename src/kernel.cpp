
#include "kernel.h"

//define cubic SPH Kernel where r is distance btwn them
double Wab(double r, double H)
{   double a2 = 1.0/(pi*pow(H,3));
    //for 1D: a2 = 2.0;     //2/3H;
    //for 2D: a2 = 10/7;    //10/7/M_PI/H/H;

    double W = 0;
    r /= H;
    //std::cout << "r/h = " << r << std::endl;
    if (r < 1.0){
        W = (1.0 - (3.0/2.0)*r*r + (3.0/4.0)*r*r*r);}
    else if (r < 2.0){
        W = (0.25*(2.0 - r)*(2.0 - r)*(2.0 - r));}
    else{
        W = 0.0;}
    
    //std::cout << "W: " << W << std::endl;
    //std::cout << "W*a2: " << W*a2 << std::endl;
    return W * a2;
}

// define gradient of cubic kernel function
double DWab(Point<3,double> & dx, Point<3,double> & DW, double r, double H)
{   
    double a2 = 1.0/(pi*pow(H,3));
    //for 1D: a2 = 2.0;     //2/3H;
    //for 2D: a2 = 10/7;    //10/7/M_PI/H/H;
    
    const double c1 = -3.0; //-3.0/M_PI/H/H/H/H;
    const double d1 = 9.0/4.0;  //9.0/4.0/M_PI/H/H/H/H;
    const double c2 = -3.0/4.0; //-3.0/4.0/M_PI/H/H/H/H;

    r /= H;
    double factor = 0;
    std::cout << "r/h = " << r << std::endl;
    if (r < 1.0)
        factor = (c1*r + d1*r*r);
    else if (r < 2.0)
        factor = (c2*(2.0 - r)*(2.0 - r));
    else
        factor = 0.0;

    factor = factor * (a2/H);
    
    //std::cout << "a/h = " << a2/H << std::endl;
    //std::cout << "factor = " << factor << std::endl;
    //std::cout << "dx_i = " << dx.get(0) << std::endl;
    DW.get(0) = factor * dx.get(0);
    DW.get(1) = factor * dx.get(1);
    DW.get(2) = factor * dx.get(2);
    return factor;
}

