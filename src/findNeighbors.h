#ifndef _neighbors_h
#define _neighbors_h

#include "Vector/vector_dist.hpp"
#include <math.h>
#include <iostream>
#include "kernel.h"
#include "global.h"

template<typename CellList> void find_neighbors(particleset  & vd, CellList & NN, double & max_visc, double H);

#endif // _neighbors_h