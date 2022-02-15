#include <stddef.h>
#include "Vector/vector_dist.hpp"
#include "timer.hpp"
#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"
#include <math.h>

#include <iostream>
#include <stdlib.h>
using namespace std;

//include functions
#include "kernel.h"
#include "findNeighbors.h"

//Create global variable class
 class Cfd
  {
    public:
    int nsteps;
    int nparticles;
    int frame;
    double dt,dx,dy,dz;
    double rad;
  };
  class engine
  { 
      public:
      double bore,stroke,conRod,crankRad,Rcomp;
      int rpm,rps;
      double Vdisp,VBDC,VTDC;
  };


// Initialize global velocity/force
constexpr int velocity = 0;
constexpr int force = 0;
const double pi = 3.14159265358979323846;

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
    simulation.nparticles = 1000;
    simulation.nsteps = 100;
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
    Ghost<3,double> ghost(simulation.rad);   // extended to contain interaction radius

    //** VECTOR INSTANTIATION **//
    //contains two vectoral properties
    vector_dist<3,double, aggregate<double[3],double[3]> > vd(0,domain,bc,ghost);
    vector_dist<3,double,aggregate<size_t,double,  double,    double,     double,     double[3], double[3], double[3]>> particleset;

    size_t cnt = 0; //used later    
    
    openfpm::vector<std::string> names1({"velocity","force"});
    openfpm::vector<std::string> names2({"type","rho","rho_prev","Pressure","drho","force","velocity","velocity_prev"});
    vd.setPropNames(names1);
    particleset.setPropNames(names2);

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
        vd.template getProp<velocity>(key)[0] = numberx;
        vd.template getProp<velocity>(key)[1] = numbery;
        vd.template getProp<velocity>(key)[2] = numberz;
        
        // next particle
        ++it;
        cnt++;
    }

    //** Map Particles to grid **/
    vd.map();       //distribute particle positions across the processors       
    vd.ghost_get<>();   //syncs the ghost with the newly mapped particles

    cout << " ---------  particles initialized  ------- " << cnt << endl;
 
    // test kernel function
    double test12 = Wab(0.001);
    cout << " ---------  kernel output 0.15  ------- " << test12 << endl;

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

    const double H = simulation.rad;
    auto NN = vd.getCellList(2*H);
    double temp = 0;
    // Time loop
    for (size_t i = 0; i < simulation.nsteps ; i++)
    {
        auto it3 = vd.getDomainIterator();  //iterator that traverses the particles in the domain

        // Calculate velocities to update positions for each particle
        while (it3.isNext())
        {
            auto p = it3.get();

            //find_neighbors(vd, NN, temp, H);

            // v = v + .5dt calculate v(tn + 0.5) += 0.5*dt;
            // velocity is always dependent on the previous velocity (getProp)
            vd.template getProp<velocity>(p)[0] += (0.5-vd.template getProp<velocity>(p)[0])*dt/0.1;
            vd.template getProp<velocity>(p)[1] += (0.5-vd.template getProp<velocity>(p)[1])*dt/0.1;
            vd.template getProp<velocity>(p)[2] += (0.5-vd.template getProp<velocity>(p)[2])*dt/0.1;

            // update particle position based on velocity
            vd.getPos(p)[0] += vd.template getProp<velocity>(p)[0]*dt;
            vd.getPos(p)[1] += vd.template getProp<velocity>(p)[1]*dt;
            vd.getPos(p)[2] += vd.template getProp<velocity>(p)[2]*dt;

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