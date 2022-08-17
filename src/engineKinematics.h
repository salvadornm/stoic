
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

void updateSimulation(Cfd &simulation, engine eng){

    simulation.ppv = simulation.nparticles/eng.volumeC;
    simulation.dp = std::cbrt(1/simulation.ppv);
    simulation.H = 3*simulation.dp;
    std::cout << "simulation.H: " << simulation.H << std::endl; 

}

void pistonInteraction(particleset &vd, Cfd &simulation, engine eng){

    auto it = vd.getDomainIterator();
    double Tmean;
    int count = 0;

    //calculate mean temperature from previous iteration
    while (it.isNext())
    {
        auto a = it.get();
        Tmean += vd.template getProp<i_temperature>(a);
        
        ++it; ++count;
    }
    Tmean = Tmean/count;
    
    //fluid mean density
    double rhomean = simulation.m_tot/eng.volumeC;
    double Pmean = rhomean * Tmean;
    cout << "rhomean: " << rhomean << endl;

    auto it2 = vd.getDomainIterator();
    while (it2.isNext())
    {   
        auto p = it2.get();    //contains (i,j,k) index of grid
        int p1 = p.getKey();
        vd.template getProp<i_rho>(p) = rhomean;  //[pa] atmospheric pressure <- EQTN TO UPDATE THIS?

        updateThermalProperties1(vd, p1);    //equation of state: update pressure and temperature
        //output_vd(vd,p1);

        //update particle position/velocity IF impacted by piston move
        // update temporary particle position
        Point <3,double> pos {vd.getPos(p)[0],vd.getPos(p)[1],vd.getPos(p)[2]};
        vector <double> vel {vd.template getProp<i_velocity>(p)[0],vd.template getProp<i_velocity>(p)[1],vd.template getProp<i_velocity>(p)[2]};

        //check particle inside cylinder
        int bc_flag = inCylinder(pos, eng);
        while (bc_flag == 0){
            pistonBound(pos, vel, eng);
            topBound(pos, vel, eng);

            bc_flag = inCylinder(pos, eng);
        }
    
        vd.getPos(p)[0] = pos[0];
        vd.getPos(p)[1] = pos[1];
        vd.getPos(p)[2] = pos[2];
        vd.template getProp<i_velocity>(p)[0] = vel[0];
        vd.template getProp<i_velocity>(p)[1] = vel[1];
        vd.template getProp<i_velocity>(p)[2] = vel[2];   

        ++it2;
    }
}
 

#endif // _engineKinematics_h