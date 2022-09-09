#ifndef _global_h
#define _global_h

#include "Vector/vector_dist.hpp"

// // properties index
constexpr int i_velocity    = 0;
constexpr int i_momentum    = 0;
constexpr int i_rho         = 1;
constexpr int i_energy      = 2;
constexpr int i_visc        = 2;
constexpr int i_pressure    = 3;
constexpr int i_temperature = 4;
constexpr int i_scalars     = 5;
constexpr int i_species     = 6;
constexpr int i_vdmean      = 7;
constexpr int i_dvdmean     = 8;
constexpr int i_velx        = 0;
constexpr int i_vely        = 5;
constexpr int i_velz        = 6;

constexpr int NDIM          = 3;
constexpr int NVARSOLVE     = 5;
constexpr int NVAR          = 7;

// Initialize global vars
const double pi = 3.14159265358979323846;
const double R_air = 287; //[J/kgK]
const double R_global = 8.31; //[J/mol/K]
const double cp_global = 1000; //[J/kgK]
const double cv_global = 718; //[J/kgK]

//Create global variable class
 class Cfd
  {
    public:
    int nsteps;
    int nparticles;
    int frame;
    double dt,dx,dy,dz, dp;
    double rad, ppv, H, r_cut;
    double lx, ly, lz;
    double m_tot, Eta2;
    double P0,T0, Pmean, Tmean, Rhomean;
  };
  class engine
  { 
      public:
      double bore,stroke,conRod,crankRad,Rcomp;
      double Nrpm,Nrps, smp;
      double Vdisp,VBDC,VTDC;
      double volumeC, height;
      double s_inst,V_inst, dStime, dVol; //instantaneous stroke and volume
      double ca,ca_init;  //crank angle
      int flag; //compression or expansion stroke
      double Twall;
  };

  class thermal
  { 
      public:
      double k, cp, cv, R;
  };
  
  class turbulence
  { 
      public:
      double k_sgs, Eps_sgs, C0, delta;
  };

// particle structure
typedef vector_dist<3,double,aggregate<double[3], double,  double,    double,  double, double, double, double[7], double[3][5]>> particleset;
//typedef vector_dist<3,double,aggregate<double[3], double,  double,    double,  double, double, double>> particleset;
//typedef vector_dist<1,double,aggregate<double,  double,    double,  double, double>> gradientset;

#endif // _global_h