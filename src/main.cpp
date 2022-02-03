#include <stddef.h>
#include "Vector/vector_dist.hpp"
#include "timer.hpp"
#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"

#include <iostream>
#include <stdlib.h>
using namespace std;

//Create global variable class
 class Cfd
  {
    public:
    int nsteps;
    int nparticles;
    int frame;
    double dt,dx,dy,dz;
    double delta;
  };
  class engine
  { 
      public:
      double bore,stroke,conRod,crankRad,Rcomp;
      int rpm,rps;
      double Vdisp,VBDC,VTDC;
  };

  class grid_prop
  {
      public:
      int nx,ny,nz;
      int nr,nphi;
      int dim;
      double lx,ly,lz;
      double dx,dy,dz;
      double ang;
  };

// Initialize global velocity/force
constexpr int velocity = 0;
constexpr int force = 0;
const double pi = 3.14159265358979323846;


int main(int argc, char* argv[])
{
    cout << "Hello World \n" << endl;
    
    //** VARIABLES **//
    Cfd simulation;
    std::default_random_engine generator;

    //simulation parameters
    simulation.nparticles = 1000;
    simulation.nsteps = 100;
    simulation.dt = 0.01;
    simulation.frame = 10;
    simulation.delta = 0.01;

    //basic engine params
    engine eng;
    eng.bore = 12.065;
    eng.stroke = 14.0;
    
    double lx = eng.bore/2;
    double ly = eng.stroke + eng.VTDC/(pi*pow(eng.bore,2)/4);

    // implements a 1D std:vector like structure to create grid   
    openfpm::vector<double> x;
    openfpm::vector<openfpm::vector<double>> y;

    //gaussian distribtuion (0,1)
    std::normal_distribution<double> distribution(0.0,1.0);

    //** INITIALIZATION **//
    //initialize openfpm
    openfpm_init(&argc,&argv);

    Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});  //define 3D box
    size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};    // boundary conditions: Periodic means the 1 boundary is equal to the 0 boundary
    Ghost<3,double> ghost(simulation.delta);   // extended to contain interaction radius

    //** VECTOR INSTANTIATION **//
    //contains two vectoral properties
    vector_dist<3,double, aggregate<double[3],double[3]> > vd(0,domain,bc,ghost);

    size_t cnt = 0; //used later    
    
    openfpm::vector<std::string> names({"velocity","force"});
    vd.setPropNames(names);

    //**ASSIGN Random POSITION**//
    //should generate all values
    auto it = vd.getDomainIterator();
    for (size_t i = 0; i < simulation.nparticles ; i++)
    {
        vd.add();
        auto key = it.get();    //contains (i,j,k) index of grid

        // we define x, assign a random position between 0.0 and 1.0
        vd.getPos(key)[0] = (double)rand() / RAND_MAX;  //rand_max just normalizes it to between 0 and 1
        vd.getPos(key)[1] = (double)rand() / RAND_MAX;
        vd.getPos(key)[2] = (double)rand() / RAND_MAX;

        //set random property of the particles
        double numberx = distribution(generator);
        double numbery = distribution(generator);
        double numberz = distribution(generator);

        //set the property of the particles
        vd.template getProp<velocity>(key)[0] = 2.0;
        vd.template getProp<velocity>(key)[1] = 1.0;
        vd.template getProp<velocity>(key)[2] = 0.0;
        
        // next particle
        ++it;
        cnt++;
    }

    //** Map Particles to grid **/
    vd.map();       //distribute particle positions across the processors       
    vd.ghost_get<>();   //syncs the ghost with the newly mapped particles

    cout << " ---------  particles initialized  ------- " << cnt << endl;
 
    //**ASSIGNING VALUES**//

    timer tsim;
    tsim.start();
    double dt = simulation.dt;
    unsigned long int f = 0;

    auto NN = vd.getCellList<CELL_MEMBAL(3,double)>(simulation.delta); //cell list structure

    // MD time stepping
    for (size_t i = 0; i < simulation.nsteps ; i++)
    {
        auto it3 = vd.getDomainIterator();  //iterator that traverses the particles in the domain

        // Calculate velocities to update positions
        while (it3.isNext())
        {
            auto p = it3.get();

            // v = v + .5dt calculate v(tn + 0.5)
            // velocity is always dependent on the previous velocity (getProp)
            //vd.template getProp<velocity>(p)[0] += 0.5*dt;
            //vd.template getProp<velocity>(p)[1] += 0.5*dt;
            //vd.template getProp<velocity>(p)[2] += 0.5*dt;

            // calculate x(tn + 1)
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
    //A VTK file contains information about particle position and properties


    // save vtk format (vtk is always the default)
    //vd.write("particles_moving");
    // save in vtk format with time
    //vd.write("particles_moving_with_time","time=1.234");
    // save in vtk format with time
    //vd.write("particles_moving_with_time_bin","time=1.234",VTK_WRITER | FORMAT_BINARY);
    
    //De - initialize openfpm
    openfpm_finalize();
}