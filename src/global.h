#ifndef _global_h
#define _global_h

#include "Vector/vector_dist.hpp"

// // properties index
// const int i_rho         = 1;
// const int i_energy      = 2;
// const int i_pressure    = 3;
// const int i_temperature = 4;
// const int i_velocity    = 5;
// const int i_scalars     = 6;
// const int i_species     = 7;

// Initialize global velocity/force
constexpr int velocity = 0;
constexpr int force = 0;
const double pi = 3.14159265358979323846;

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

// particle structure
typedef vector_dist<3,double,aggregate<size_t,double,  double,    double,  double, double[3], double, double>> particleset;

#endif // _global_h