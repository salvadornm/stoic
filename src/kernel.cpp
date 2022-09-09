
#include "kernel.h"

double kern_a2 = 1.0/pi;    
//W* = W * H^NDIM -> a2 = a2 * H^NDIM
//for 1D: a2 = 2.0/3H;   
//for 2D: a2 = 10/(7*pi*pow(H,2));
//for 3D: a2 = 1.0/(pi*pow(H,3));

//define cubic SPH Kernel where r is distance btwn points
double Wab(double r, double H)
{   
    double W = 0;
    double r_h = r / H;
    
    //std::cout << "r/h = " << r_h << std::endl;
    if (r < 1.0){
        W = (1.0 - (3.0/2.0)*r_h*r_h + (3.0/4.0)*r_h*r_h*r_h);}
    else if (r < 2.0){
        W = (0.25*(2.0 - r_h)*(2.0 - r_h)*(2.0 - r_h));}
    else{
        W = 0.0;}

    return W * kern_a2;
}

// define gradient of cubic kernel function
double * DWab(Point<3,double> dr, double r, double H)
{       
    static double dW[NDIM];
    double factor = 0;
    double r_h = r / H;
    
    if (r_h < 1.0)
        factor = -3.0*r_h + (9.0/4.0)*r_h*r_h;
    else if (r_h < 2.0)
        factor = (-3.0/4.0)*(2.0 - r_h)*(2.0 - r_h);
    else
        factor = 0.0;

    factor = factor * (kern_a2/H);
    
    for (int i = 0; i < NDIM ; i++) 
    {
        //dW[i] =   -dr.get(i)*factor/abs(dr.get(i))/h;
        dW[i] =   1.0/(dr.get(i) + 1.E-08);
    }

    return dW;
}

