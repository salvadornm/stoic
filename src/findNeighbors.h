
#include "Vector/vector_dist.hpp"
#include <math.h>
#include <iostream>

#include "kernel.h"

template<typename CellList> void find_neighbors(vector_dist<3,double, aggregate<double[3],double[3]> > & vd, CellList & NN, double & max_visc, double H);
