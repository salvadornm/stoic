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
#include "findNeighbors.h"
#include "calculations.h"
#include "updateProps.h"
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
    // to do: create separate input parameter function/file
    Cfd simulation;
    std::default_random_engine generator;

    //simulation parameters
    simulation.nparticles = 1000; //1000
    simulation.nsteps = 10; //100
    simulation.dt = 0.01;
    simulation.frame = 1;   //10
    simulation.rad = 2;
    simulation.dp = 1/sqrt(simulation.nparticles);

    //basic engine params --> make simulation dx, dy, dz
    engine eng;
    eng.bore = 12.065/100;
    eng.stroke = 14.0;
    eng.volumeC = 1; //eng.bore*eng.bore*eng.bore;

    simulation.ppv = simulation.nparticles/eng.volumeC;
    simulation.dp = std::cbrt(1/simulation.ppv);
    simulation.H = 3*simulation.dp;
    cout << "simulation.H: " << simulation.H << endl;
    simulation.lx = 1;  //eng.bore;
    simulation.ly = 1;  //eng.bore;
    simulation.lz = 1;  //eng.bore;

    //turbulences
    turbulence turb;
    thermal therm;
    
    //time step variables
    double Vol = 1; //[m3]
    double mass_mix, mass_p; 

    //-- initialize openfpm --//
    openfpm_init(&argc,&argv);

    // implements a 1D std:vector like structure to create grid   
    openfpm::vector<double> x;
    openfpm::vector<openfpm::vector<double>> y;
    std::normal_distribution<double> distribution(0.0, 1.0);//(0.0,1.0);

    //Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
    Box<3,double> domain({0.0,0.0,0.0},{simulation.lx,simulation.ly,simulation.lz});
    size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};    // boundary conditions: Periodic means the 1 boundary is equal to the 0 boundary

    // extended boundary around the domain, and the processor domain
    Ghost<3,double> g(0.1); //g(2*H); //why using g(2*H)??

    // initialize particlesets
    particleset vd(0,domain,bc,g,DEC_GRAN(512));  

    openfpm::vector<std::string> names({"velocity","rho","energy","Pressure","Temperature","scalars","species","vdmean","dvdmean"});
    openfpm::vector<std::string> grad_names({"momentum","density","energy","Pressure","Temperature"});
    vd.setPropNames(names);
    
    size_t cnt = 0; //used later  
    //auto obstacle_box = DrawParticles::DrawSkin(vd,bc,domain,piston);

    
    //--ADD PARTICLES--//  ------------------------------------
    auto it = vd.getDomainIterator();
    
    for (size_t i = 0; i < simulation.nparticles ; i++)
    {
        vd.add();
        auto key = it.get();    //contains (i,j,k) index of grid

        // we define x, assign a random position between 0.0 and 1.0
        vd.getPos(key)[0] = ((double)rand() / RAND_MAX) * simulation.lx;  //rand_max just normalizes it to between 0 and 1
        vd.getPos(key)[1] = ((double)rand() / RAND_MAX) * simulation.ly;
        vd.getPos(key)[2] = ((double)rand() / RAND_MAX) * simulation.lz;

        //set random velocity of the particles
        double numberx = 1; //distribution(generator);
        double numbery = 1; //distribution(generator);
        double numberz = 1; //distribution(generator);

        //set the property of the particles : eventually velocity will be initialized from turbulence files
        vd.template getProp<i_velocity>(key)[0] = numberx;
        vd.template getProp<i_velocity>(key)[1] = numbery;
        vd.template getProp<i_velocity>(key)[2] = numberz;

        //initialize remaining properties (placeholder values for now)
        vd.template getProp<i_temperature>(key) = 300; //[K]
        vd.template getProp<i_pressure>(key) = 101300;  //[pa] atmospheric pressure
        vd.template getProp<i_energy>(key) = 1; //temporary placeholder
        vd.template getProp<i_rho>(key) = .1; //temporary placeholder
        
        //initialize dvdmean particles
        for (size_t j = 0; j < 7.0 ; j++)
        { vd.template getProp<i_vdmean>(key)[j] = 0.0; }
        for (size_t j = 0; j < 3.0 ; j++)
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

    
    auto NN = vd.getCellList(r_cut);  //define neighborhoods with radius (repeated at end of time loop)
    // x = NN.getNNIteratorRadius();

    // Time loop
    for (size_t i = 0; i < simulation.nsteps ; i++)
    {
        auto it3 = vd.getDomainIterator();  //iterator that traverses the particles in the domain 
        std::cout << "step: " << i << std::endl;
        find_neighbors(vd, NN); //contaions properties of neighbors

        //function to solve for new cylinder geometry
        count = 0;

        // Particle loop...
        while (it3.isNext())
        {
            count++;
            auto p = it3.get();
            int place = p.getKey();

            std::cout << count << " og vd particle " << std::endl;
            output_vd(vd,place);
            
            //updateEqtnState(vd);    //calc pressure based on local density
            double a5 = vd.getProp<i_velocity>(p)[0];
            double a6 = vd.getProp<i_velocity>(p)[1];
            double a7 = vd.getProp<i_velocity>(p)[2];
            double b5 = vd.getProp<i_vdmean>(p)[i_velx];
            double b6 = vd.getProp<i_vdmean>(p)[i_vely];
            double b7 = vd.getProp<i_vdmean>(p)[i_velz];

            //std::cout << count << " og vd particle " << std::endl;
            //std::cout << " vx = " << a5 << "    vy = " << a6 << "   vz = " << a7 << std::endl;
            //std::cout << " density = " << vd.getProp<i_rho>(p) << std::endl;

            updateParticleProperties(vd, place, dt, H, turb);
            //update thermal properties; pressure, total energy
            //updateThermalProperties2(vd, vdmean, place);
            //updateThermoDynamics
            //advanceNavierStokes
            //cout << "New pressure = " << vd.template getProp<i_pressure>(p) << " new temp = "<< vd.template getProp<i_temperature>(p)  << endl;
                        
            moveParticles(vd, place, dt, eng);
            
            std::cout << "updated vd particle --------" << std::endl;
            std::cout << " vx = " << vd.getProp<i_velocity>(p)[0] << "      vy = " << vd.getProp<i_velocity>(p)[1] << "     vz = " << vd.getProp<i_velocity>(p)[2] << std::endl;
            std::cout << " density = " << vd.getProp<i_rho>(p) << std::endl;
            
            /*
            std::cout << "(vdmean particle) --------" << std::endl;
            std::cout << " vx = " << b5 << "    vy = " << b6 << "   vz = " << b7 << std::endl;
            */
            ++it3;

        }
        //cout << "post move" << endl;
        //stateOfNeighbors(vd, NN);
        
        //moveParticles(vd, place, dt);
        //applyBC(vd,place,dt);
        //updateThermalProperties
        //continuityEquation

        //MOVE BOUNDARY/UPDATE PISTON LOCATION
        //moveCrank(eng);

        // Map particles and re-sync the ghost
        vd.map();
        vd.template ghost_get<>();        
        NN = vd.getCellList(r_cut);

        // collect statistics on the configuration
        //if (i+1 % simulation.frame == 0)
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
          cout << " --------- STEP   i="  << i  << std::endl;  
        }  
         
    } 

    tsim.stop();
    std::cout << "Total Time: " << tsim.getwct() << std::endl;

    cout << " ---------  END  ------- "  << endl;
    openfpm_finalize(); //De - initialize openfpm
}