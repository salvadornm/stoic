
#ifndef _inputParams_h
#define _inputParams_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"


void inputUserParams(Cfd &simulation, engine &eng )
{
    //simulation parameters
    simulation.nparticles = 10000; //
    simulation.nsteps = 3; //100
    simulation.dt = 0.00005;   //0.01
    simulation.frame = 1;   //10
    simulation.rad = 2;
    simulation.dp = 1/sqrt(simulation.nparticles);

    //engine geometery parameters (for the volvo td100)
    eng.bore = 12.065/100;  //[m]
    eng.stroke = 14.0/100;  //[m]
    eng.conRod = 26.0/100;  //[m]
    eng.Rcomp = 15;
    eng.Nrpm = 1500;   //[rpm]
    eng.ca_init = 180.0;  //start at BDC

    eng.Twall = 300; //[K]
}

void initialize_vel(particleset &vd, Cfd &simulation, double key, engine eng)
{
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::normal_distribution<double> distribution(0.0, eng.smp);

    //set random velocity of the particles
    double numberx = distribution(generator);
    double numbery = distribution(generator);
    double numberz = distribution(generator);

    //set the property of the particles : eventually velocity will be initialized from turbulence files
    vd.template getProp<i_velocity>(key)[0] = 0;// numberx;
    vd.template getProp<i_velocity>(key)[1] = 0;// numbery;
    vd.template getProp<i_velocity>(key)[2] = 0;// numberz;
    
    //use to test different velocity inputs per location in cylinder  
    /*
    // Get the distance between a and b
    double r_cyl = eng.bore/2;
    Point<3,double> cyl_center {r_cyl, r_cyl};
    Point<3,double> pos {vd.getPos(key)[0], vd.getPos(key)[1]};
    //double radius = sqrt(norm2(pos - r_cyl));
    double radius = sqrt(pow((vd.getPos(key)[0]-r_cyl),2)+pow((vd.getPos(key)[1]-r_cyl),2));

    if(radius > 0.2*eng.bore){
        vd.template getProp<i_velocity>(key)[0] = 0;
        vd.template getProp<i_velocity>(key)[1] = 0;
        vd.template getProp<i_velocity>(key)[2] = 0;
    }
    */
}

void initialize_temp(particleset &vd, Cfd &simulation, double key, engine eng)
{
    // Get the distance between a and b
    double r_cyl = eng.bore/2;
    Point<3,double> cyl_center {r_cyl, r_cyl};
    Point<3,double> pos {vd.getPos(key)[0], vd.getPos(key)[1]};
    
    double radius = sqrt(pow((vd.getPos(key)[0]-r_cyl),2)+pow((vd.getPos(key)[1]-r_cyl),2));

    if(radius < 0.2*eng.bore){
        vd.template getProp<i_temperature>(key) = 2500.0;   //hot product temperature
    } else{
        vd.template getProp<i_temperature>(key) = 700.0;    //cool reactant temperature
    }
}

void initialize_pres(particleset &vd, Cfd &simulation, double key, engine eng)
{
    // Get the distance between a and b
    double r_cyl = eng.bore/2;
    Point<3,double> cyl_center {r_cyl, r_cyl};
    Point<3,double> pos {vd.getPos(key)[0], vd.getPos(key)[1]};
    
    double radius = sqrt(pow((vd.getPos(key)[0]-r_cyl),2)+pow((vd.getPos(key)[1]-r_cyl),2));

    if(radius < 0.2*eng.bore){
        vd.template getProp<i_pressure>(key) = 10*101300;   //hot product temperature
    } else{
        vd.template getProp<i_pressure>(key) = 101300;
    }

    //initialize pressure for kernel testing in 1 Dimension only!
    //vd.template getProp<i_pressure>(key) = 4*radius;
    //vd.template getProp<i_pressure>(key) = 101300*vd.getPos(key)[0];
}

void initialize_dvdmean(particleset &vd, double key)
{
    //initialize dvdmean particles
    for (size_t j = 0; j < 7 ; j++)
    {  
        vd.template getProp<i_vdmean>(key)[j] = 0.0; 
    }
    for (size_t j = 0; j < 3 ; j++)
    { 
        vd.template getProp<i_dvdmean>(key)[j][i_momentum]  = 0.0;
        vd.template getProp<i_dvdmean>(key)[j][i_rho]  = 0.0;
        vd.template getProp<i_dvdmean>(key)[j][i_energy]  = 0.0;
        vd.template getProp<i_dvdmean>(key)[j][i_pressure]  = 0.0;
        vd.template getProp<i_dvdmean>(key)[j][i_temperature]  = 0.0;
    }
}


#endif // _inputParams_h