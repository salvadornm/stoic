
#ifndef _engineKinematics_h
#define _engineKinematics_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"

void initialize_geometry(engine &eng){

    eng.Nrps = eng.Nrpm/60;   //[rps]
    eng.crankRad = eng.stroke/2;
    eng.smp = 2*eng.stroke*eng.Nrps;
    eng.ca = eng.ca_init;

    //volumes
    eng.volumeC = (pi/4)*eng.bore*eng.bore*eng.stroke;
    eng.Vdisp = (pi/4)*(eng.bore*eng.bore*eng.stroke);
    eng.VBDC = eng.Vdisp * (eng.Rcomp / (eng.Rcomp-1));
    eng.VTDC = eng.VBDC - eng.Vdisp;
}

void update_CA(double dt, engine &eng){
    double N = eng.Nrps;    //revolutions per second
    double rotation = N*dt;     //xx revolutions

    eng.ca += rotation*360;
}

//update instantaneous volume,stroke. return piston height
double movePiston(engine &eng){
    double a = eng.crankRad;    //crank radius
    double l = eng.conRod;
    double theta = eng.ca*(pi/180);     //convert to radians
    
    eng.s_inst = a - l + sqrt(l*l - a*a*sin(theta)*sin(theta)) + a*cos(theta);
    
    double yt = 2*a - eng.s_inst;                           //distance of piston from top
    eng.volumeC = eng.VTDC + (pi/4)*eng.bore*eng.bore*yt;   //updated volume of cylinder

    return yt;
}
 

#endif // _engineKinematics_h