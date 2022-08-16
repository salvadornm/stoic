
#ifndef _engineKinematics_h
#define _engineKinematics_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"

void initialize_geometry(Cfd &simulation, engine &eng){

    eng.Nrps = eng.Nrpm/60;   //[rps]
    eng.crankRad = eng.stroke/2;
    eng.smp = 2*eng.stroke*eng.Nrps;

    //initialize
    eng.ca = eng.ca_init;
    eng.s_inst = 0;

    //volumes
    eng.Vdisp = (pi/4) * eng.bore * eng.bore * eng.stroke;  //[m3]
    eng.VBDC = eng.Vdisp * (eng.Rcomp / (eng.Rcomp-1));
    eng.VTDC = eng.VBDC - eng.Vdisp;
    eng.height = eng.stroke + ((4*eng.VTDC) / (pi*eng.bore*eng.bore));
    eng.volumeC = (pi/4)*eng.bore*eng.bore*eng.height;
    cout << "vol_c: " << eng.volumeC << " VTDC: " << eng.VTDC << " VBDC: " << eng.VBDC << endl;
    cout << "eng height: " << eng.height << endl;


    //calculations
    simulation.ppv = simulation.nparticles/eng.volumeC;
    simulation.dp = std::cbrt(1/simulation.ppv);
    simulation.H = 3*simulation.dp;
    std::cout << "simulation.H: " << simulation.H << std::endl; //eventually move global.h H value to this... placeholder - not done yet
    simulation.lx = 2 * eng.bore;
    simulation.ly = 2 * eng.bore;
    simulation.lz = 2 * eng.height;
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
    eng.volumeC = eng.VTDC + (pi/4)*eng.bore*eng.bore*(yt);   //updated volume of cylinder

    //st = eng.stroke - yt;
    return eng.s_inst;
}
 

#endif // _engineKinematics_h