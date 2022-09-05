#ifndef _kernel_h
#define _kernel_h

#include "Vector/vector_dist.hpp"
#include <math.h>
#include "global.h"


//SPH Kernel --> kernel.cpp
double Wab(double r, double H);
double DWab(Point<3,double> & dx, Point<3,double> & DW, double r, double H);

#endif // _kernel_h