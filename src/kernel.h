
#include "Vector/vector_dist.hpp"
#include <math.h>
#include <iostream>

//SPH Kernel --> kernel.cpp
double Wab(double r);
void DWab(Point<3,double> & dx, Point<3,double> & DW, double r, bool print);
