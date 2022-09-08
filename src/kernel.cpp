
#include "kernel.h"

double kern_a2 = 1.0/pi ;

//define cubic SPH Kernel where r is distance btwn them
double Wab(double r, double h)
{   
    double W = 0;
    r /= h;
    //std::cout << "r/h = " << r << std::endl;
    if (r < 1.0){
        W = 1.0 - 3.0/2.0*r*r + (3.0/4.0)*r*r*r;}
    else if (r < 2.0){
        W  =0.25*(2.0 - r)*(2.0 - r)*(2.0 - r);}
    else{
        W = 0.0;}

    return W *  kern_a2;
}

// define gradient of cubic kernel function
double * DWab(double r, double h,Point<3,double> dr)
{   
    static double dW[NDIM];
    double factor = 0.0;
    double q = r / h;
    
    if (q < 1.0)
        factor = -3.0*q + 9.0/4.0*q*q;
    else if (q < 2.0)
        factor = -3.0/4.0*(2.0 - q)*(2.0 - q);
    else
        factor = 0.0;

    factor = factor * kern_a2;

    for (int i = 0; i < NDIM ; i++) 
    {
        //dW[i] =   -dr.get(i)*factor/abs(dr.get(i))/h;
        dW[i] =   1.0/(dr.get(i) + 1.E-08);
    }        
  
    


    return dW;
}

