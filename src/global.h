#ifndef _global_h
#define _global_h

#include "Vector/vector_dist.hpp"

// // properties index
constexpr int i_velocity    = 0;
constexpr int i_momentum    = 0;
constexpr int i_rho         = 1;
constexpr int i_energy      = 2;
constexpr int i_pressure    = 3;
constexpr int i_temperature = 4;
constexpr int i_scalars     = 5;
constexpr int i_species     = 6;
constexpr int i_vdmean     = 7;
constexpr int i_dvdmean     = 8;


// Initialize global vars
const double pi = 3.14159265358979323846;
const double H = 0.0147224318643; //for kernel
const double Eta2 = 0.01 * H*H;
const double R_air = 287; //[J/kg/K]
const double R_global = .831; //[J/mol/K]

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

  class thermal
  { 
      public:
      double k, cp, cv, R;
  };
  
  class turbulence
  { 
      public:
      double k_sgs, Eps_sgs;
  };

// particle structure
typedef vector_dist<3,double,aggregate<double[3], double,  double,    double,  double, double, double, double[7], double[5][4]>> particleset;
//typedef vector_dist<3,double,aggregate<double[3], double,  double,    double,  double, double, double>> particleset;
//typedef vector_dist<1,double,aggregate<double,  double,    double,  double, double>> gradientset;

#endif // _global_h