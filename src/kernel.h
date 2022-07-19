#ifndef _kernel_h
#define _kernel_h

#include "Vector/vector_dist.hpp"
#include <math.h>
#include "global.h"

const double a2 = 1.0; //1/M_PI/H/H/H;
//for 1D: a2 = 2.0;     //2/3/H;
//for 2D: a2 = 10/7;    //10/7/M_PI/H/H;

const double c1 = -3.0; //-3.0/M_PI/H/H/H/H;
const double d1 = 9.0/4.0;  //9.0/4.0/M_PI/H/H/H/H;
const double c2 = -3.0/4.0; //-3.0/4.0/M_PI/H/H/H/H;
const double a2_4 = 0.25*a2;

//SPH Kernel --> kernel.cpp
double Wab(double r);
double DWab(Point<3,double> & dx, Point<3,double> & DW, double r, bool print);

#endif // _kernel_h