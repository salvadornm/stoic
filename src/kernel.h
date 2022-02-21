#ifndef _kernel_h
#define _kernel_h

#include "Vector/vector_dist.hpp"
#include <math.h>
#include <iostream>

//const double H = 0.0147224318643;

//SPH Kernel --> kernel.cpp
double Wab(double r);
void DWab(Point<3,double> & dx, Point<3,double> & DW, double r, bool print);

#endif // _kernel_h