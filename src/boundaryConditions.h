#ifndef _boundaryConditions_h
#define _boundaryConditions_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"
#include <math.h>

//determine if point p is inside the cylinder
//return 0 if outside cylinder, 1 if inside cylinder
int inCylinder(double posx, double posy, double posz, int p, engine eng, double & psi)
{  // Get the distance between a and b
    double r_cyl = eng.bore/2;
    Point<2,double> xa {posx,posy};
    Point<2,double> center {r_cyl, r_cyl};

    Point<2,double> dr = xa - center;
    double r2 = norm2(dr);  //norm2 = (sum of the squares)
    double R = sqrt(r2);

    //is the particle within the z(height) boundaries...
    //assume piston @ BDC and not moving (if moving: vd.getpos(p)[2] < eng.stroke - y) , where y is instantaneous distance from TDC  
    if (posz < 0 || posz - eng.stroke > 0){
        return 0;   //out of bounds
    }    

    //is particle within xy plane boundaries....
    if (r_cyl - R < 0) {
        psi = r_cyl - R;
        return 0;   //out of bounds
    }
    else{   // if (r_cyl - R > 0)
        return 1;   //in bounds
    }
}

//determine if point p is inside the cylinder
//return 0 if outside cylinder, 1 if inside cylinder
void initialBoundary(double & posx, double & posy, engine eng, Cfd sim)
{  // Get the distance between a and b
    double r_cyl = eng.bore/2;
    double r = r_cyl * sqrt(((double)rand() / RAND_MAX));
    double theta = ((double)rand()/ RAND_MAX) * 2 * pi;

    posx = r_cyl + r * cos(theta);
    posy = r_cyl + r * sin(theta);
}
#endif // _boundaryConditions_h