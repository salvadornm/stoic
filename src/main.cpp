#include <stddef.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "Vector/vector_dist.hpp"
#include "Draw/DrawParticles.hpp"
#include "timer.hpp"
#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"

using namespace std;

//include functions
#include "global.h"
#include "vector_math.h"
#include "inputParams.h"
#include "findNeighbors.h"
#include "calculations.h"
#include "updateProps.h"
#include "boundaryConditions.h"
#include "moveParticles.h"
#include "test.h"


//selections for run setup ------------------------------------
#define BOUNDARY 0 // A constant to indicate boundary particles
#define FLUID 1 // A constant to indicate fluid particles

int main(int argc, char* argv[])
{
    cout << "Hello World \n" << endl;
    cout << "H = " << H << endl;
    
    //-- VARIABLES --// ------------------------------------

    Cfd simulation;
    engine eng;
    std::default_random_engine generator;

    //turbulences
    turbulence turb;
    thermal therm;

    //-- initialize user inputs --//
    inputUserParams(simulation,eng);

    //-- initialize openfpm --//
    openfpm_init(&argc,&argv);

    // implements a 1D std:vector like structure to create grid   
    openfpm::vector<double> x;
    openfpm::vector<openfpm::vector<double>> y;
    std::normal_distribution<double> distribution(0.0, simulation.lx);//(0.0,1.0);

    Box<3,double> domain({0.0,0.0,0.0},{simulation.lx,simulation.ly,simulation.lz});
    size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};    // boundary conditions: Periodic means the 1 boundary is equal to the 0 boundary

    // extended boundary around the domain, and the processor domain
    Ghost<3,double> g(0.1); //g(2*H);

    // initialize particlesets
    particleset vd(0,domain,bc,g,DEC_GRAN(512));  

    openfpm::vector<std::string> names({"velocity","rho","energy","Pressure","Temperature","scalars","species","vdmean","dvdmean"});
    openfpm::vector<std::string> grad_names({"momentum","density","energy","Pressure","Temperature"});
    vd.setPropNames(names);
    
    size_t cnt = 0; //used later  
    
    //--ADD PARTICLES--//  ------------------------------------
    auto it = vd.getDomainIterator();
    
    for (size_t i = 0; i < simulation.nparticles ; i++)
    {
        vd.add();
        auto key = it.get();    //contains (i,j,k) index of grid
        int key1 = key.getKey();

        // Assign a random position within engine cylinder
        vd.getPos(key)[2] = ((double)rand() / RAND_MAX) * eng.stroke;

        initialBoundary(vd.getPos(key)[0], vd.getPos(key)[1], eng, simulation);

        //set random velocity of the particles
        double numberx = distribution(generator);
        double numbery = distribution(generator);
        double numberz = distribution(generator);

        //set the property of the particles : eventually velocity will be initialized from turbulence files
        vd.template getProp<i_velocity>(key)[0] = numberx;
        vd.template getProp<i_velocity>(key)[1] = numbery;
        vd.template getProp<i_velocity>(key)[2] = numberz;


        //initialize remaining properties (placeholder values for now)
        vd.template getProp<i_temperature>(key) = 300; //[K]
        vd.template getProp<i_pressure>(key) = 101300;  //[pa] atmospheric pressure
        vd.template getProp<i_energy>(key) = 1; //temporary placeholder
        vd.template getProp<i_rho>(key) = .1; //temporary placeholder
        
        updateDensity(vd, key1);    //equation of state
        //cout << "initial energy: " << vd.template getProp<i_energy>(key) << " density: " << vd.template getProp<i_rho>(key) << endl;

        //vary initial pressure/density (fx in test.cpp)
        //vary_initialization(vd, simulation, key1);
        
        //initialize dvdmean particles
        for (size_t j = 0; j < 7 ; j++)
        { vd.template getProp<i_vdmean>(key)[j] = 0.0; }
        for (size_t j = 0; j < 3 ; j++)
        { vd.template getProp<i_dvdmean>(key)[j][i_momentum]  = 0.0;
          vd.template getProp<i_dvdmean>(key)[j][i_rho]  = 0.0;
          vd.template getProp<i_dvdmean>(key)[j][i_energy]  = 0.0;
          vd.template getProp<i_dvdmean>(key)[j][i_pressure]  = 0.0;
          vd.template getProp<i_dvdmean>(key)[j][i_temperature]  = 0.0;
        }
        
        //updateChemicalProperties(vd); //initialize temperature and composition
        //updateThermalProperties1(vd);   //update conductivity,diffusivity,specific haet capacity, based on T/Y 
        //updateDensity(vd);
        //checkContinuity
                
        // next particle
        ++it;
        cnt++;
    }

    //-- Map Particles to grid --/
    vd.map();       //distribute particle positions across the processors
    vd.write_frame("particles_",0);       
    vd.ghost_get<>();   //syncs the ghost with the newly mapped particles

    cout << " ---------  particles initialized  ------- " << cnt << endl;


    //--PARTICLE TIME LOOP--//  ------------------------------------

    timer tsim;
    tsim.start();
    double dt = simulation.dt;
    unsigned long int f = 0;
    int count = 0;
    double cfl = 0;
    
    auto NN = vd.getCellList(r_cut);  //define neighborhoods with radius (repeated at end of time loop)

    // Time loop
    for (size_t i = 0; i < simulation.nsteps ; i++)
    {
        auto it3 = vd.getDomainIterator();  //iterator that traverses the particles in the domain 
        std::cout << "--------step: " << i << " ------" << std::endl;
        find_neighbors(vd, NN); //contaions properties of neighbors

        //function to solve for new cylinder geometry
        count = 0;

        // Particle loop...
        while (it3.isNext())
        {
            count++;
            auto p = it3.get();
            int place = p.getKey();

           // std::cout << count << " particle " << std::endl;
           // output_vd(vd,place);    //output particle properties
            
            updateParticleProperties(vd, place, dt, H, turb, simulation);

            // smooth P

            double cfl_temp = (dt/H) * (vd.template getProp<i_velocity>(p)[0]+vd.template getProp<i_velocity>(p)[1]+vd.template getProp<i_velocity>(p)[2]);
            cfl = std::max(abs(cfl_temp),cfl);       //look for max not average. needs to be absolute.  should hold
       
            moveParticles(vd, place, dt, eng);
            
            //updateThermalProperties1(cd, place);   //update pressure/temperature equation of state

            ++it3;
        }

        
        std::cout << "cfl: " << cfl << endl;
        cfl = 0;

        //MOVE BOUNDARY/UPDATE PISTON LOCATION
        //moveCrank(eng);

        // Map particles and re-sync the ghost
        vd.map();
        vd.template ghost_get<>();        
        NN = vd.getCellList(r_cut);

        // collect statistics on the configuration
        if ((i+1) % simulation.frame == 0)
        {
            // Write the particle position for visualization (Without ghost)
            vd.deleteGhost();
            vd.write_frame("particles_",i+1);
            
            // resync the ghost
            vd.ghost_get<>();

            //--Reduce--//
            auto & v_cl = create_vcluster();
            v_cl.sum(cnt);
            v_cl.execute();
        }

        if (i % 10 ==0 )
        {
          cout << " Step "  << i  << std::endl;  
        }  
         
    } 

    tsim.stop();
    std::cout << "Total Time: " << tsim.getwct() << std::endl;

    cout << " ---------  END  ------- "  << endl;
    openfpm_finalize(); //De - initialize openfpm
}