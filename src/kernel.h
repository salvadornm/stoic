#ifndef _kernel_h
#define _kernel_h

#include "Vector/vector_dist.hpp"
#include <math.h>
#include "global.h"


//SPH Kernel --> kernel.cpp
double Wab(double r, double h);
//double DWab(double r, double h,Point<3,double> dr, Point<3,double> & DW );
double * DWab(double r, double h,Point<3,double> dr);

#endif // _kernel_h