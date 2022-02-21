#ifndef _kernel_h
#define _kernel_h

#include "Vector/vector_dist.hpp"

//SPH Kernel --> kernel.cpp
double Wab(double r);
void DWab(Point<3,double> & dx, Point<3,double> & DW, double r, bool print);

#endif // _kernel_h