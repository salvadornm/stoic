#ifndef _global_h
#define _global_h

#include "Vector/vector_dist.hpp"

// // properties index
constexpr int i_velocity    = 0;
constexpr int i_rho         = 1;
constexpr int i_energy      = 2;
constexpr int i_pressure    = 3;
constexpr int i_temperature = 4;
constexpr int i_scalars     = 5;
constexpr int i_species     = 6;

// Initialize global velocity/force
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
typedef vector_dist<3,double,aggregate<double[3], double,  double,    double,  double, double, double>> particleset;

#endif // _global_h