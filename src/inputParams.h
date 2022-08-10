
#ifndef _inputParams_h
#define _inputParams_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"


void inputUserParams(Cfd &simulation, engine &eng)
{
    //simulation parameters
    simulation.nparticles = 5000; //if you change this change the H value in global.h!
    simulation.nsteps = 100; //100
    simulation.dt = 0.003;   //0.01
    simulation.frame = 2;   //10
    simulation.rad = 2;
    simulation.dp = 1/sqrt(simulation.nparticles);

    //engine geometery parameters
    eng.bore = 12.065/100;
    eng.stroke = 14.0/100;
    eng.volumeC = (pi/4)*eng.bore*eng.bore*eng.stroke;

    //calculations
    simulation.ppv = simulation.nparticles/eng.volumeC;
    simulation.dp = std::cbrt(1/simulation.ppv);
    simulation.H = 3*simulation.dp;
    //SNM
    //simulation.H = 0.005;   //0.01  //changed in global.h, simulation.H does nothing right now...
    
    
    std::cout << "simulation.H: " << simulation.H << std::endl;
    simulation.lx = 2 * eng.bore;
    simulation.ly = 2 * eng.bore;
    simulation.lz = 2 * eng.stroke;

}
#endif // _inputParams_h