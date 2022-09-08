
#include "kernel.h"

//define cubic SPH Kernel where r is distance btwn them
double Wab(double r, double H)
{   
    //W* = W * H^d
    //double a2 = 1.0/(pi*pow(H,3));
    double a2 = 1.0/pi;    

    //for 1D: a2 = 2.0/3H;   
    //for 2D: a2 = 10/(7*pi*pow(H,2));
    //for 3D: a2 = 1.0/(pi*pow(H,3));

    double W = 0;
    r /= H;
    //std::cout << "r/h = " << r << std::endl;
    if (r < 1.0){
        W = (1.0 - (3.0/2.0)*r*r + (3.0/4.0)*r*r*r);}
    else if (r < 2.0){
        W = (0.25*(2.0 - r)*(2.0 - r)*(2.0 - r));}
    else{
        W = 0.0;}

    return W * a2;
}

// define gradient of cubic kernel function
double DWab(Point<3,double> & dx, Point<3,double> & DW, double r, double H)
{   
    //W* = W * H^d
    //double a2 = 1.0/(pi*pow(H,3));
    double a2 = 1.0/pi;  

    //for 1D: a2 = 2.0/3H;   
    //for 2D: a2 = 10/(7*pi*pow(H,2));
    //for 3D: a2 = 1.0/(pi*pow(H,3));
    
    const double c1 = -3.0; 
    const double d1 = 9.0/4.0;  
    const double c2 = -3.0/4.0; 

    double r_h = r / H;
    double factor = 0;
    
    if (r_h < 1.0)
        factor = (c1*r_h + d1*r_h*r_h);
    else if (r_h < 2.0)
        factor = (c2*(2.0 - r_h)*(2.0 - r_h));
    else
        factor = 0.0;

    //factor = factor*a2;
    factor = factor * (a2/H);
    
    DW.get(0) = factor * dx.get(0);
    DW.get(1) = factor * dx.get(1);
    DW.get(2) = factor * dx.get(2);
    return factor;
}

