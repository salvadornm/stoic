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
    simulation.nsteps = 100; //100
    simulation.dt = 0.01;
    simulation.frame = 10;
    simulation.rad = 2;

    //basic engine params
    engine eng;
    eng.bore = 12.065;
    eng.stroke = 14.0;

    // implements a 1D std:vector like structure to create grid   
    openfpm::vector<double> x;
    openfpm::vector<openfpm::vector<double>> y;

    // test gaussian distribtuion (0,1)
    std::normal_distribution<double> distribution(0.0,1.0);

    //** initialize openfpm **//
    openfpm_init(&argc,&argv);

    Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});  //define 3D box
    size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};    // boundary conditions: Periodic means the 1 boundary is equal to the 0 boundary
    //Ghost<3,double> ghost(simulation.rad);   // extended to contain interaction radius

    // extended boundary around the domain, and the processor domain
    Ghost<3,double> g(2*H);

    particleset vd(0,domain,bc,g,DEC_GRAN(512)); 
    particleset vdmean(0,domain,bc,g,DEC_GRAN(512)); //added today


    size_t cnt = 0; //used later    
    
    openfpm::vector<std::string> names({"velocity","rho","energy","Pressure","Temperature","scalars","species"});
    vd.setPropNames(names);
    vdmean.setPropNames(names);

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
        
        // next particle
        ++it;
        cnt++;
    }

    //** Map Particles to grid **/
    vd.map();       //distribute particle positions across the processors       
    vd.ghost_get<>();   //syncs the ghost with the newly mapped particles

    vdmean = vd;

    cout << " ---------  particles initialized  ------- " << cnt << endl;


    //**ASSIGNING VALUES**//

    timer tsim;
    tsim.start();
    double dt = simulation.dt;
    unsigned long int f = 0;

    //in the time loop:
    // in the particle loop:
        //find nearest neighbors

        //compute averages
        //compute gradient

    //advance particles
        //updateParticleProperties
        //moveParticles
    
    //check boundary

    //move boundary (movePiston)

   // const double H = simulation.rad;
    auto NN = vd.getCellList(2*H);
    // Time loop
    for (size_t i = 0; i < simulation.nsteps ; i++)
    {
        auto it3 = vd.getDomainIterator();  //iterator that traverses the particles in the domain

        // For each particle...
        while (it3.isNext())
        {
            auto p = it3.get();
            
            // (this may break due to int overflow with high numbers of particles, change to double)
            int place = p.getKey();
            
            updateEqtnState(vd);    //calc pressure based on local density

            find_neighbors(vd, vdmean, NN, H); //contaions properties of neighbors
            updateParticleProperties(vd, vdmean, place, dt);

            moveParticles(vd, vdmean, place, dt);
            
            ++it3;
        }

        // Map particles and re-sync the ghost
        vd.map();
        vd.template ghost_get<>();

        // collect some statistic about the configuration
        if (i % simulation.frame == 0)
        {
            // Write the particle position for visualization (Without ghost)
            vd.deleteGhost();
            vd.write_frame("particles_",i);

            // resync the ghost
            vd.ghost_get<>();

            //**Reduce**//
            auto & v_cl = create_vcluster();
            v_cl.sum(cnt);
            v_cl.execute();
        }
    }

    tsim.stop();
    std::cout << "Time: " << tsim.getwct() << std::endl;

    cout << " ---------  END  ------- "  << endl;

    //**VISUALIZATION**// <--- NEEDS UPDATED STILL

    //De - initialize openfpm
    openfpm_finalize();
}