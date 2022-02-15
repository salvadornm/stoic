#include "Vector/vector_dist.hpp"

// // properties index
// const int i_rho         = 1;
// const int i_energy      = 2;
// const int i_pressure    = 3;
// const int i_temperature = 4;
// const int i_velocity    = 5;
// const int i_scalars     = 6;
// const int i_species     = 7;

// particle structure
typedef vector_dist<3,double,aggregate<size_t,double,  double,    double,  double, double[3], double, double>> particleset;
