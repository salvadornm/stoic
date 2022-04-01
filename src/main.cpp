#include <stddef.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "Vector/vector_dist.hpp"
#include "timer.hpp"
#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"

using namespace std;

//include functions
#include "global.h"
#include "findNeighbors.h"
#include "calculations.h"
#include "updateProps.h"
#include "moveParticles.h"

#define BOUNDARY 0 // A constant to indicate boundary particles
#define FLUID 1 // A constant to indicate fluid particles

int main(int argc, char* argv[])
{
    cout << "Hello World \n" << endl;
    
    //** VARIABLES **//
    // to do: create separate input parameter function/file
    Cfd simulation;
    std::default_random_engine generator;

    //simulation parameters
    simulation.nparticles = 100; //1000
    simulation.nsteps = 50; //100
    simulation.dt = 0.01;
    simulation.frame = 10;
    simulation.rad = 2;

    //basic engine params
    engine eng;
    eng.bore = 12.065;
    eng.stroke = 14.0;

    //time step variables
    double Vol = 1; //[m3]
    double mass_mix, mass_p; 

    // implements a 1D std:vector like structure to create grid   
    openfpm::vector<double> x;
    openfpm::vector<openfpm::vector<double>> y;
    std::normal_distribution<double> distribution(0.0,1.0);

    //** initialize openfpm **//
    openfpm_init(&argc,&argv);

    Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});  //define 3D box
    size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};    // boundary conditions: Periodic means the 1 boundary is equal to the 0 boundary
    //Ghost<3,double> ghost(simulation.rad);   // extended to contain interaction radius

    // extended boundary around the domain, and the processor domain
    Ghost<3,double> g(2*H);

    particleset vd(0,domain,bc,g,DEC_GRAN(512)); 
    particleset vdmean(0,domain,bc,g,DEC_GRAN(512)); 
    particleset dvdmean(0,domain,bc,g,DEC_GRAN(512)); 
    particleset dvdmeany(0,domain,bc,g,DEC_GRAN(512)); 
    particleset dvdmeanz(0,domain,bc,g,DEC_GRAN(512)); //added today

    size_t cnt = 0; //used later    
    
    openfpm::vector<std::string> names({"velocity","rho","energy","Pressure","Temperature","scalars","species"});
    vd.setPropNames(names);
    vdmean.setPropNames(names);
    dvdmean.setPropNames(names);
    dvdmeany.setPropNames(names);
    dvdmeanz.setPropNames(names);
    

    //**ADD PARTICLES**//
    auto it = vd.getDomainIterator();
    for (size_t i = 0; i < simulation.nparticles ; i++)
    {
        vd.add();
        auto key = it.get();    //contains (i,j,k) index of grid

        // we define x, assign a random position between 0.0 and 1.0
        vd.getPos(key)[0] = (double)rand() / RAND_MAX;  //rand_max just normalizes it to between 0 and 1
        vd.getPos(key)[1] = (double)rand() / RAND_MAX;
        vd.getPos(key)[2] = (double)rand() / RAND_MAX;

        //set random velocity of the particles
        double numberx = distribution(generator);
        double numbery = distribution(generator);
        double numberz = distribution(generator);

        //set the property of the particles
        vd.template getProp<i_velocity>(key)[0] = numberx;
        vd.template getProp<i_velocity>(key)[1] = numbery;
        vd.template getProp<i_velocity>(key)[2] = numberz;

        //set remaining properties
        vd.template getProp<i_temperature>(key) = 300; //[K]
        vd.template getProp<i_pressure>(key) = 101300;  //[pa] atmospheric pressure
        vd.template getProp<i_rho>(key) = 850; //[kg/m3] temporary placeholder
        
        // next particle
        ++it;
        cnt++;
    }

    //** Map Particles to grid **/
    vd.map();       //distribute particle positions across the processors       
    vd.ghost_get<>();   //syncs the ghost with the newly mapped particles

    // initialise mean
    vdmean = vd;
    // initilaise grads (values no meaning)
    dvdmean = vd;
    dvdmeany = vd;
    dvdmeanz = vd;

    cout << " ---------  particles initialized  ------- " << cnt << endl;


    //**ASSIGNING VALUES**//

    timer tsim;
    tsim.start();
    double dt = simulation.dt;
    unsigned long int f = 0;

        
    auto NN = vd.getCellList(4*H);
    // Time loop
    for (size_t i = 0; i < simulation.nsteps ; i++)
    {
        auto it3 = vd.getDomainIterator();  //iterator that traverses the particles in the domain

        //find new volume of cylinder
        // crank angle;
        Vol = 1; //[m3]
        mass_mix = 0;
        //calculate mass of particle
        mass_p = mass_mix/simulation.nparticles;

        find_neighbors(vd, vdmean,dvdmean, NN, H); //contaions properties of neighbors
 

        // Particle loop...
        while (it3.isNext())
        {
            auto p = it3.get();
            int place = p.getKey();
            
            //updateEqtnState(vd);    //calc pressure based on local density

            updateParticleProperties(vd, vdmean, dvdmean, place, dt, H);
            moveParticles(vd, place, dt);
            
            vdmean.getPos(p)[0] = vd.getPos(p)[0];
            vdmean.getPos(p)[1] = vd.getPos(p)[1];
            vdmean.getPos(p)[2] = vd.getPos(p)[2];
            
            ++it3;
        }

        //MOVE BOUNDARY/UPDATE PISTON LOCATION

        // Map particles and re-sync the ghost
        vd.map();
        vd.template ghost_get<>();        
        NN = vd.getCellList(4*H);

        // collect some statistic about the configuration
        if (i % simulation.frame == 0)
        {
            // Write the particle position for visualization (Without ghost)
            vd.deleteGhost();
            vd.write_frame("particles_",i);
            vdmean.write_frame("partmean_",i);
            
            // resync the ghost
            vd.ghost_get<>();

            //**Reduce**//
            auto & v_cl = create_vcluster();
            v_cl.sum(cnt);
            v_cl.execute();
        }

        if (i % 10 ==0 )
        {
          cout << " --------- STEP   i="  << i  << std::endl;  
        }
    }

    tsim.stop();
    std::cout << "Time: " << tsim.getwct() << std::endl;

    cout << " ---------  END  ------- "  << endl;
    openfpm_finalize(); //De - initialize openfpm
}