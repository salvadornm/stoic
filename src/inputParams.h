
#ifndef _inputParams_h
#define _inputParams_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"


void inputUserParams(Cfd &simulation, engine &eng )
{
    //simulation parameters
    simulation.nparticles = 10000; //if you change this change the H value in global.h!
    simulation.nsteps = 50; //100
    simulation.dt = 0.003;   //0.01
    simulation.frame = 2;   //10
    simulation.rad = 2;
    simulation.dp = 1/sqrt(simulation.nparticles);

    //engine geometery parameters (for the volvo td100)
    eng.bore = 12.065/100;  //[m]
    eng.stroke = 14.0/100;  //[m]
    eng.conRod = 26.0/100;  //[m]
    eng.Rcomp = 15;
    eng.Nrpm = 1500;   //[rpm]
    eng.ca_init = 180.0;  //start at BDC

    //calculations
    simulation.ppv = simulation.nparticles/eng.volumeC;
    simulation.dp = std::cbrt(1/simulation.ppv);
    simulation.H = 3*simulation.dp;
    std::cout << "simulation.H: " << simulation.H << std::endl; //eventually move global.h H value to this... placeholder - not done yet
    simulation.lx = 2 * eng.bore;
    simulation.ly = 2 * eng.bore;
    simulation.lz = 2 * eng.stroke;

}

void initialize_vel(particleset &vd, double key, double lx)
{
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::normal_distribution<double> distribution(0.0, lx);

    //set random velocity of the particles
    double numberx = distribution(generator);
    double numbery = distribution(generator);
    double numberz = distribution(generator);

    //set the property of the particles : eventually velocity will be initialized from turbulence files
    vd.template getProp<i_velocity>(key)[0] = numberx;
    vd.template getProp<i_velocity>(key)[1] = numbery;
    vd.template getProp<i_velocity>(key)[2] = numberz;
}

void initialize_temp(particleset &vd, Cfd &simulation, double key, engine eng)
{
    // Get the distance between a and b
    double r_cyl = eng.bore/2;
    Point<3,double> cyl_center {r_cyl, r_cyl};
    Point<3,double> pos {vd.getPos(key)[0], vd.getPos(key)[1]};
    //double radius = sqrt(norm2(pos - r_cyl));
    double radius = sqrt(pow((vd.getPos(key)[0]-r_cyl),2)+pow((vd.getPos(key)[1]-r_cyl),2));

    if(radius < 0.1*eng.bore){
        vd.template getProp<i_temperature>(key) = 2500.0;   //hot product temperature
    } else{
        vd.template getProp<i_temperature>(key) = 700.0;    //cool reactant temperature
    }
}


#endif // _inputParams_h