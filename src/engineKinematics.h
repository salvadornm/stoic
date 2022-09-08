
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
    eng.flag = 0;

    //volumes
    eng.Vdisp = (pi/4) * eng.bore * eng.bore * eng.stroke;  //[m3]
    eng.VBDC = eng.Vdisp * (eng.Rcomp / (eng.Rcomp-1));
    eng.VTDC = eng.VBDC - eng.Vdisp;
    eng.height = eng.stroke + ((4*eng.VTDC) / (pi*eng.bore*eng.bore));
    eng.volumeC = (pi/4)*eng.bore*eng.bore*eng.height;

    //calculations
    simulation.ppv = simulation.nparticles/eng.volumeC;
    simulation.dp = std::cbrt(1/simulation.ppv);
    simulation.H = (5/3)*simulation.dp;
  
    simulation.H = 0.02;

    std::cout << "simulation.H: " << simulation.H << std::endl; 
    simulation.lx = 2 * eng.bore;
    simulation.ly = 2 * eng.bore;
    simulation.lz = 2 * eng.height;
    std::cout << "Simulation.dt: " << simulation.dt << std::endl;
}

void update_CA(double dt, engine &eng){
    double N = eng.Nrps;    //revolutions per second
    double rotation = N*dt;     //xx revolutions

    eng.ca += rotation*360;
}

void updateSimulation(Cfd &simulation, engine eng){

    simulation.ppv = simulation.nparticles/eng.volumeC;
    simulation.dp = std::cbrt(1/simulation.ppv);
    simulation.H = (5/3)*simulation.dp; //*simulation.dp;
    simulation.r_cut = 2*simulation.H;
    simulation.Eta2 = 0.01 * simulation.H * simulation.H;
    //std::cout << "volumec: " << eng.volumeC << std::endl;
    //std::cout << "simulation.H: " << simulation.H << std::endl; 

}

//update instantaneous volume,stroke. return piston height
void movePiston(engine &eng){
    double a = eng.crankRad;    //crank radius
    double l = eng.conRod;
    double theta = eng.ca*(pi/180);     //convert to radians
    double s_temp = 0; double vol_temp = 0;
    
    s_temp = a - l + sqrt(l*l - a*a*sin(theta)*sin(theta)) + a*cos(theta);

    if (s_temp < eng.s_inst){eng.flag = 0;} //expansion stroke
    else{eng.flag = 1;} //compression stroke

    eng.dStime = s_temp - eng.s_inst;
    
    eng.s_inst = s_temp;
    double yt = 2*a - eng.s_inst;                           //distance of piston from top
    vol_temp = eng.VTDC + (pi/4)*pow(eng.bore,2)*(yt);   //updated volume of cylinder
    
    eng.dVol = vol_temp - eng.volumeC;
    eng.volumeC = vol_temp;
    //st = eng.stroke - yt;
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
    //cout << "pmean: " << Pmean << endl;

    auto it2 = vd.getDomainIterator();
    while (it2.isNext())
    {   
        auto p = it2.get();    //contains (i,j,k) index of grid
        int p1 = p.getKey();
        vd.template getProp<i_rho>(p) = rhomean;  //[pa] 

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